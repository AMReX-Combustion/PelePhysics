#ifndef MECHANISM_CPP
#define MECHANISM_CPP

#include "chemistry_file.H"

/* PURE CPU stuff  */
#ifndef AMREX_USE_GPU
namespace thermo
{
    amrex::Real fwd_A[84], fwd_beta[84], fwd_Ea[84];
    amrex::Real low_A[84], low_beta[84], low_Ea[84];
    amrex::Real rev_A[84], rev_beta[84], rev_Ea[84];
    amrex::Real troe_a[84],troe_Ts[84], troe_Tss[84], troe_Tsss[84];
    amrex::Real sri_a[84], sri_b[84], sri_c[84], sri_d[84], sri_e[84];
    amrex::Real activation_units[84], prefactor_units[84], phase_units[84];
    int is_PD[84], troe_len[84], sri_len[84], nTB[84], *TBid[84];
    amrex::Real *TB[84];
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
    amrex::Real k_f_s[84*npt], Kc_s[84*npt], mixture[npt], g_RT[21*npt];
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
    vcomp_wdot_51_84(npt, wdot, mixture, sc, k_f_s, Kc_s, tc, invT, T);
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

        Kc_s[0*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[7*npt+i]) - (g_RT[9*npt+i]));
        Kc_s[1*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[9*npt+i]) - (g_RT[10*npt+i]));
        Kc_s[2*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[13*npt+i]) - (g_RT[14*npt+i]));
        Kc_s[3*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[14*npt+i]) - (g_RT[15*npt+i]));
        Kc_s[4*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[16*npt+i]) - (g_RT[17*npt+i]));
        Kc_s[5*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[17*npt+i]) - (g_RT[18*npt+i]));
        Kc_s[6*npt+i] = refCinv * exp((g_RT[0*npt+i] + g_RT[11*npt+i]) - (g_RT[14*npt+i]));
        Kc_s[7*npt+i] = refCinv * exp((2.000000 * g_RT[9*npt+i]) - (g_RT[18*npt+i]));
        Kc_s[8*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[2*npt+i]) - (g_RT[4*npt+i]));
        Kc_s[9*npt+i] = refCinv * exp((g_RT[2*npt+i] + g_RT[11*npt+i]) - (g_RT[12*npt+i]));
        Kc_s[10*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[3*npt+i]) - (g_RT[6*npt+i]));
        Kc_s[11*npt+i] = refCinv * exp((2.000000 * g_RT[1*npt+i]) - (g_RT[0*npt+i]));
        Kc_s[12*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[4*npt+i]) - (g_RT[5*npt+i]));
        Kc_s[13*npt+i] = refC * exp((g_RT[13*npt+i]) - (g_RT[1*npt+i] + g_RT[11*npt+i]));
        Kc_s[14*npt+i] = exp((g_RT[0*npt+i] + g_RT[2*npt+i]) - (g_RT[1*npt+i] + g_RT[4*npt+i]));
        Kc_s[15*npt+i] = exp((g_RT[2*npt+i] + g_RT[6*npt+i]) - (g_RT[3*npt+i] + g_RT[4*npt+i]));
        Kc_s[16*npt+i] = exp((g_RT[2*npt+i] + g_RT[7*npt+i]) - (g_RT[1*npt+i] + g_RT[13*npt+i]));
        Kc_s[17*npt+i] = exp((g_RT[2*npt+i] + g_RT[8*npt+i]) - (g_RT[1*npt+i] + g_RT[13*npt+i]));
        Kc_s[18*npt+i] = exp((g_RT[2*npt+i] + g_RT[9*npt+i]) - (g_RT[1*npt+i] + g_RT[14*npt+i]));
        Kc_s[19*npt+i] = exp((g_RT[2*npt+i] + g_RT[10*npt+i]) - (g_RT[4*npt+i] + g_RT[9*npt+i]));
        Kc_s[20*npt+i] = exp((g_RT[2*npt+i] + g_RT[13*npt+i]) - (g_RT[4*npt+i] + g_RT[11*npt+i]));
        Kc_s[21*npt+i] = exp((g_RT[2*npt+i] + g_RT[13*npt+i]) - (g_RT[1*npt+i] + g_RT[12*npt+i]));
        Kc_s[22*npt+i] = exp((g_RT[2*npt+i] + g_RT[14*npt+i]) - (g_RT[4*npt+i] + g_RT[13*npt+i]));
        Kc_s[23*npt+i] = exp((g_RT[2*npt+i] + g_RT[16*npt+i]) - (g_RT[9*npt+i] + g_RT[13*npt+i]));
        Kc_s[24*npt+i] = exp((g_RT[2*npt+i] + g_RT[17*npt+i]) - (g_RT[9*npt+i] + g_RT[14*npt+i]));
        Kc_s[25*npt+i] = exp((g_RT[2*npt+i] + g_RT[18*npt+i]) - (g_RT[4*npt+i] + g_RT[17*npt+i]));
        Kc_s[26*npt+i] = exp((g_RT[3*npt+i] + g_RT[11*npt+i]) - (g_RT[2*npt+i] + g_RT[12*npt+i]));
        Kc_s[27*npt+i] = exp((g_RT[3*npt+i] + g_RT[14*npt+i]) - (g_RT[6*npt+i] + g_RT[13*npt+i]));
        Kc_s[28*npt+i] = refCinv * exp((g_RT[1*npt+i] + 2.000000 * g_RT[3*npt+i]) - (g_RT[3*npt+i] + g_RT[6*npt+i]));
        Kc_s[29*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[3*npt+i] + g_RT[5*npt+i]) - (g_RT[5*npt+i] + g_RT[6*npt+i]));
        Kc_s[30*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[3*npt+i] + g_RT[19*npt+i]) - (g_RT[6*npt+i] + g_RT[19*npt+i]));
        Kc_s[31*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[3*npt+i] + g_RT[20*npt+i]) - (g_RT[6*npt+i] + g_RT[20*npt+i]));
        Kc_s[32*npt+i] = exp((g_RT[1*npt+i] + g_RT[3*npt+i]) - (g_RT[2*npt+i] + g_RT[4*npt+i]));
        Kc_s[33*npt+i] = refCinv * exp((g_RT[0*npt+i] + 2.000000 * g_RT[1*npt+i]) - (2.000000 * g_RT[0*npt+i]));
        Kc_s[34*npt+i] = refCinv * exp((2.000000 * g_RT[1*npt+i] + g_RT[5*npt+i]) - (g_RT[0*npt+i] + g_RT[5*npt+i]));
        Kc_s[35*npt+i] = refCinv * exp((2.000000 * g_RT[1*npt+i] + g_RT[12*npt+i]) - (g_RT[0*npt+i] + g_RT[12*npt+i]));
        Kc_s[36*npt+i] = exp((g_RT[1*npt+i] + g_RT[6*npt+i]) - (g_RT[0*npt+i] + g_RT[3*npt+i]));
        Kc_s[37*npt+i] = exp((g_RT[1*npt+i] + g_RT[6*npt+i]) - (2.000000 * g_RT[4*npt+i]));
        Kc_s[38*npt+i] = exp((g_RT[1*npt+i] + g_RT[10*npt+i]) - (g_RT[0*npt+i] + g_RT[9*npt+i]));
        Kc_s[39*npt+i] = exp((g_RT[1*npt+i] + g_RT[13*npt+i]) - (g_RT[0*npt+i] + g_RT[11*npt+i]));
        Kc_s[40*npt+i] = exp((g_RT[1*npt+i] + g_RT[14*npt+i]) - (g_RT[0*npt+i] + g_RT[13*npt+i]));
        Kc_s[41*npt+i] = exp((g_RT[1*npt+i] + g_RT[15*npt+i]) - (g_RT[4*npt+i] + g_RT[9*npt+i]));
        Kc_s[42*npt+i] = exp((g_RT[1*npt+i] + g_RT[18*npt+i]) - (g_RT[0*npt+i] + g_RT[17*npt+i]));
        Kc_s[43*npt+i] = exp((g_RT[0*npt+i] + g_RT[4*npt+i]) - (g_RT[1*npt+i] + g_RT[5*npt+i]));
        Kc_s[44*npt+i] = exp((2.000000 * g_RT[4*npt+i]) - (g_RT[2*npt+i] + g_RT[5*npt+i]));
        Kc_s[45*npt+i] = exp((g_RT[4*npt+i] + g_RT[6*npt+i]) - (g_RT[3*npt+i] + g_RT[5*npt+i]));
        Kc_s[46*npt+i] = exp((g_RT[4*npt+i] + g_RT[7*npt+i]) - (g_RT[1*npt+i] + g_RT[14*npt+i]));
        Kc_s[47*npt+i] = exp((g_RT[4*npt+i] + g_RT[8*npt+i]) - (g_RT[1*npt+i] + g_RT[14*npt+i]));
        Kc_s[48*npt+i] = exp((g_RT[4*npt+i] + g_RT[9*npt+i]) - (g_RT[5*npt+i] + g_RT[7*npt+i]));
        Kc_s[49*npt+i] = exp((g_RT[4*npt+i] + g_RT[9*npt+i]) - (g_RT[5*npt+i] + g_RT[8*npt+i]));
        Kc_s[50*npt+i] = exp((g_RT[4*npt+i] + g_RT[10*npt+i]) - (g_RT[5*npt+i] + g_RT[9*npt+i]));
        Kc_s[51*npt+i] = exp((g_RT[4*npt+i] + g_RT[11*npt+i]) - (g_RT[1*npt+i] + g_RT[12*npt+i]));
        Kc_s[52*npt+i] = exp((g_RT[4*npt+i] + g_RT[13*npt+i]) - (g_RT[5*npt+i] + g_RT[11*npt+i]));
        Kc_s[53*npt+i] = exp((g_RT[4*npt+i] + g_RT[14*npt+i]) - (g_RT[5*npt+i] + g_RT[13*npt+i]));
        Kc_s[54*npt+i] = exp((g_RT[4*npt+i] + g_RT[18*npt+i]) - (g_RT[5*npt+i] + g_RT[17*npt+i]));
        Kc_s[55*npt+i] = exp((g_RT[6*npt+i] + g_RT[7*npt+i]) - (g_RT[4*npt+i] + g_RT[14*npt+i]));
        Kc_s[56*npt+i] = exp((g_RT[6*npt+i] + g_RT[9*npt+i]) - (g_RT[3*npt+i] + g_RT[10*npt+i]));
        Kc_s[57*npt+i] = exp((g_RT[6*npt+i] + g_RT[9*npt+i]) - (g_RT[4*npt+i] + g_RT[15*npt+i]));
        Kc_s[58*npt+i] = exp((g_RT[6*npt+i] + g_RT[11*npt+i]) - (g_RT[4*npt+i] + g_RT[12*npt+i]));
        Kc_s[59*npt+i] = exp((g_RT[3*npt+i] + g_RT[7*npt+i]) - (g_RT[4*npt+i] + g_RT[13*npt+i]));
        Kc_s[60*npt+i] = exp((g_RT[0*npt+i] + g_RT[7*npt+i]) - (g_RT[1*npt+i] + g_RT[9*npt+i]));
        Kc_s[61*npt+i] = exp((g_RT[7*npt+i] + g_RT[9*npt+i]) - (g_RT[1*npt+i] + g_RT[16*npt+i]));
        Kc_s[62*npt+i] = exp((g_RT[7*npt+i] + g_RT[10*npt+i]) - (2.000000 * g_RT[9*npt+i]));
        Kc_s[63*npt+i] = exp((g_RT[8*npt+i] + g_RT[19*npt+i]) - (g_RT[7*npt+i] + g_RT[19*npt+i]));
        Kc_s[64*npt+i] = exp((g_RT[8*npt+i] + g_RT[20*npt+i]) - (g_RT[7*npt+i] + g_RT[20*npt+i]));
        Kc_s[65*npt+i] = refC * exp((g_RT[3*npt+i] + g_RT[8*npt+i]) - (g_RT[1*npt+i] + g_RT[4*npt+i] + g_RT[11*npt+i]));
        Kc_s[66*npt+i] = exp((g_RT[3*npt+i] + g_RT[8*npt+i]) - (g_RT[5*npt+i] + g_RT[11*npt+i]));
        Kc_s[67*npt+i] = exp((g_RT[0*npt+i] + g_RT[8*npt+i]) - (g_RT[1*npt+i] + g_RT[9*npt+i]));
        Kc_s[68*npt+i] = exp((g_RT[5*npt+i] + g_RT[8*npt+i]) - (g_RT[5*npt+i] + g_RT[7*npt+i]));
        Kc_s[69*npt+i] = exp((g_RT[8*npt+i] + g_RT[9*npt+i]) - (g_RT[1*npt+i] + g_RT[16*npt+i]));
        Kc_s[70*npt+i] = exp((g_RT[8*npt+i] + g_RT[10*npt+i]) - (2.000000 * g_RT[9*npt+i]));
        Kc_s[71*npt+i] = exp((g_RT[8*npt+i] + g_RT[11*npt+i]) - (g_RT[7*npt+i] + g_RT[11*npt+i]));
        Kc_s[72*npt+i] = exp((g_RT[8*npt+i] + g_RT[12*npt+i]) - (g_RT[7*npt+i] + g_RT[12*npt+i]));
        Kc_s[73*npt+i] = exp((g_RT[8*npt+i] + g_RT[12*npt+i]) - (g_RT[11*npt+i] + g_RT[14*npt+i]));
        Kc_s[74*npt+i] = exp((g_RT[3*npt+i] + g_RT[9*npt+i]) - (g_RT[2*npt+i] + g_RT[15*npt+i]));
        Kc_s[75*npt+i] = exp((g_RT[3*npt+i] + g_RT[9*npt+i]) - (g_RT[4*npt+i] + g_RT[14*npt+i]));
        Kc_s[76*npt+i] = exp((2.000000 * g_RT[9*npt+i]) - (g_RT[1*npt+i] + g_RT[17*npt+i]));
        Kc_s[77*npt+i] = exp((g_RT[9*npt+i] + g_RT[13*npt+i]) - (g_RT[10*npt+i] + g_RT[11*npt+i]));
        Kc_s[78*npt+i] = exp((g_RT[9*npt+i] + g_RT[14*npt+i]) - (g_RT[10*npt+i] + g_RT[13*npt+i]));
        Kc_s[79*npt+i] = exp((g_RT[9*npt+i] + g_RT[18*npt+i]) - (g_RT[10*npt+i] + g_RT[17*npt+i]));
        Kc_s[80*npt+i] = refC * exp((g_RT[5*npt+i] + g_RT[13*npt+i]) - (g_RT[1*npt+i] + g_RT[5*npt+i] + g_RT[11*npt+i]));
        Kc_s[81*npt+i] = exp((g_RT[3*npt+i] + g_RT[13*npt+i]) - (g_RT[6*npt+i] + g_RT[11*npt+i]));
        Kc_s[82*npt+i] = exp((g_RT[3*npt+i] + g_RT[15*npt+i]) - (g_RT[6*npt+i] + g_RT[14*npt+i]));
        Kc_s[83*npt+i] = exp((g_RT[3*npt+i] + g_RT[17*npt+i]) - (g_RT[6*npt+i] + g_RT[16*npt+i]));
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

        /*reaction 1: H + CH2 (+M) <=> CH3 (+M) */
        phi_f = sc[1*npt+i]*sc[7*npt+i];
        alpha = mixture[i] + (TB[0][0] - 1)*sc[0*npt+i] + (TB[0][1] - 1)*sc[5*npt+i] + (TB[0][2] - 1)*sc[10*npt+i] + (TB[0][3] - 1)*sc[11*npt+i] + (TB[0][4] - 1)*sc[12*npt+i] + (TB[0][5] - 1)*sc[18*npt+i] + (TB[0][6] - 1)*sc[20*npt+i];
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
        phi_r = sc[9*npt+i];
        Kc = Kc_s[0*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[7*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;

        /*reaction 2: H + CH3 (+M) <=> CH4 (+M) */
        phi_f = sc[1*npt+i]*sc[9*npt+i];
        alpha = mixture[i] + (TB[1][0] - 1)*sc[0*npt+i] + (TB[1][1] - 1)*sc[5*npt+i] + (TB[1][2] - 1)*sc[10*npt+i] + (TB[1][3] - 1)*sc[11*npt+i] + (TB[1][4] - 1)*sc[12*npt+i] + (TB[1][5] - 1)*sc[18*npt+i] + (TB[1][6] - 1)*sc[20*npt+i];
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
        phi_r = sc[10*npt+i];
        Kc = Kc_s[1*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[9*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;

        /*reaction 3: H + HCO (+M) <=> CH2O (+M) */
        phi_f = sc[1*npt+i]*sc[13*npt+i];
        alpha = mixture[i] + (TB[2][0] - 1)*sc[0*npt+i] + (TB[2][1] - 1)*sc[5*npt+i] + (TB[2][2] - 1)*sc[10*npt+i] + (TB[2][3] - 1)*sc[11*npt+i] + (TB[2][4] - 1)*sc[12*npt+i] + (TB[2][5] - 1)*sc[18*npt+i] + (TB[2][6] - 1)*sc[20*npt+i];
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
        phi_r = sc[14*npt+i];
        Kc = Kc_s[2*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[13*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;

        /*reaction 4: H + CH2O (+M) <=> CH3O (+M) */
        phi_f = sc[1*npt+i]*sc[14*npt+i];
        alpha = mixture[i] + (TB[3][0] - 1)*sc[0*npt+i] + (TB[3][1] - 1)*sc[5*npt+i] + (TB[3][2] - 1)*sc[10*npt+i] + (TB[3][3] - 1)*sc[11*npt+i] + (TB[3][4] - 1)*sc[12*npt+i] + (TB[3][5] - 1)*sc[18*npt+i];
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
        phi_r = sc[15*npt+i];
        Kc = Kc_s[3*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[14*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;

        /*reaction 5: H + C2H4 (+M) <=> C2H5 (+M) */
        phi_f = sc[1*npt+i]*sc[16*npt+i];
        alpha = mixture[i] + (TB[4][0] - 1)*sc[0*npt+i] + (TB[4][1] - 1)*sc[5*npt+i] + (TB[4][2] - 1)*sc[10*npt+i] + (TB[4][3] - 1)*sc[11*npt+i] + (TB[4][4] - 1)*sc[12*npt+i] + (TB[4][5] - 1)*sc[18*npt+i] + (TB[4][6] - 1)*sc[20*npt+i];
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
        phi_r = sc[17*npt+i];
        Kc = Kc_s[4*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[16*npt+i] -= qdot;
        wdot[17*npt+i] += qdot;

        /*reaction 6: H + C2H5 (+M) <=> C2H6 (+M) */
        phi_f = sc[1*npt+i]*sc[17*npt+i];
        alpha = mixture[i] + (TB[5][0] - 1)*sc[0*npt+i] + (TB[5][1] - 1)*sc[5*npt+i] + (TB[5][2] - 1)*sc[10*npt+i] + (TB[5][3] - 1)*sc[11*npt+i] + (TB[5][4] - 1)*sc[12*npt+i] + (TB[5][5] - 1)*sc[18*npt+i] + (TB[5][6] - 1)*sc[20*npt+i];
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
        phi_r = sc[18*npt+i];
        Kc = Kc_s[5*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[17*npt+i] -= qdot;
        wdot[18*npt+i] += qdot;

        /*reaction 7: H2 + CO (+M) <=> CH2O (+M) */
        phi_f = sc[0*npt+i]*sc[11*npt+i];
        alpha = mixture[i] + (TB[6][0] - 1)*sc[0*npt+i] + (TB[6][1] - 1)*sc[5*npt+i] + (TB[6][2] - 1)*sc[10*npt+i] + (TB[6][3] - 1)*sc[11*npt+i] + (TB[6][4] - 1)*sc[12*npt+i] + (TB[6][5] - 1)*sc[18*npt+i] + (TB[6][6] - 1)*sc[20*npt+i];
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
        wdot[0*npt+i] -= qdot;
        wdot[11*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;

        /*reaction 8: 2.000000 CH3 (+M) <=> C2H6 (+M) */
        phi_f = pow(sc[9*npt+i], 2.000000);
        alpha = mixture[i] + (TB[7][0] - 1)*sc[0*npt+i] + (TB[7][1] - 1)*sc[5*npt+i] + (TB[7][2] - 1)*sc[10*npt+i] + (TB[7][3] - 1)*sc[11*npt+i] + (TB[7][4] - 1)*sc[12*npt+i] + (TB[7][5] - 1)*sc[18*npt+i] + (TB[7][6] - 1)*sc[20*npt+i];
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
        phi_r = sc[18*npt+i];
        Kc = Kc_s[7*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[9*npt+i] -= 2.000000 * qdot;
        wdot[18*npt+i] += qdot;

        /*reaction 9: O + H + M <=> OH + M */
        phi_f = sc[1*npt+i]*sc[2*npt+i];
        alpha = mixture[i] + (TB[8][0] - 1)*sc[0*npt+i] + (TB[8][1] - 1)*sc[5*npt+i] + (TB[8][2] - 1)*sc[10*npt+i] + (TB[8][3] - 1)*sc[11*npt+i] + (TB[8][4] - 1)*sc[12*npt+i] + (TB[8][5] - 1)*sc[18*npt+i] + (TB[8][6] - 1)*sc[20*npt+i];
        k_f = alpha * k_f_s[8*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i];
        Kc = Kc_s[8*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;

        /*reaction 10: O + CO + M <=> CO2 + M */
        phi_f = sc[2*npt+i]*sc[11*npt+i];
        alpha = mixture[i] + (TB[9][0] - 1)*sc[0*npt+i] + (TB[9][1] - 1)*sc[3*npt+i] + (TB[9][2] - 1)*sc[5*npt+i] + (TB[9][3] - 1)*sc[10*npt+i] + (TB[9][4] - 1)*sc[11*npt+i] + (TB[9][5] - 1)*sc[12*npt+i] + (TB[9][6] - 1)*sc[18*npt+i] + (TB[9][7] - 1)*sc[20*npt+i];
        k_f = alpha * k_f_s[9*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[12*npt+i];
        Kc = Kc_s[9*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[11*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;

        /*reaction 11: H + O2 + M <=> HO2 + M */
        phi_f = sc[1*npt+i]*sc[3*npt+i];
        alpha = mixture[i] + (TB[10][0] - 1)*sc[3*npt+i] + (TB[10][1] - 1)*sc[5*npt+i] + (TB[10][2] - 1)*sc[11*npt+i] + (TB[10][3] - 1)*sc[12*npt+i] + (TB[10][4] - 1)*sc[18*npt+i] + (TB[10][5] - 1)*sc[19*npt+i] + (TB[10][6] - 1)*sc[20*npt+i];
        k_f = alpha * k_f_s[10*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i];
        Kc = Kc_s[10*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;

        /*reaction 12: 2.000000 H + M <=> H2 + M */
        phi_f = pow(sc[1*npt+i], 2.000000);
        alpha = mixture[i] + (TB[11][0] - 1)*sc[0*npt+i] + (TB[11][1] - 1)*sc[5*npt+i] + (TB[11][2] - 1)*sc[10*npt+i] + (TB[11][3] - 1)*sc[12*npt+i] + (TB[11][4] - 1)*sc[18*npt+i] + (TB[11][5] - 1)*sc[20*npt+i];
        k_f = alpha * k_f_s[11*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i];
        Kc = Kc_s[11*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= 2.000000 * qdot;

        /*reaction 13: H + OH + M <=> H2O + M */
        phi_f = sc[1*npt+i]*sc[4*npt+i];
        alpha = mixture[i] + (TB[12][0] - 1)*sc[0*npt+i] + (TB[12][1] - 1)*sc[5*npt+i] + (TB[12][2] - 1)*sc[10*npt+i] + (TB[12][3] - 1)*sc[18*npt+i] + (TB[12][4] - 1)*sc[20*npt+i];
        k_f = alpha * k_f_s[12*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i];
        Kc = Kc_s[12*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;

        /*reaction 14: HCO + M <=> H + CO + M */
        phi_f = sc[13*npt+i];
        alpha = mixture[i] + (TB[13][0] - 1)*sc[0*npt+i] + (TB[13][1] - 1)*sc[5*npt+i] + (TB[13][2] - 1)*sc[10*npt+i] + (TB[13][3] - 1)*sc[11*npt+i] + (TB[13][4] - 1)*sc[12*npt+i] + (TB[13][5] - 1)*sc[18*npt+i];
        k_f = alpha * k_f_s[13*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[11*npt+i];
        Kc = Kc_s[13*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[11*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;

        /*reaction 15: O + H2 <=> H + OH */
        phi_f = sc[0*npt+i]*sc[2*npt+i];
        k_f = k_f_s[14*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[4*npt+i];
        Kc = Kc_s[14*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;

        /*reaction 16: O + HO2 <=> OH + O2 */
        phi_f = sc[2*npt+i]*sc[6*npt+i];
        k_f = k_f_s[15*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[4*npt+i];
        Kc = Kc_s[15*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[4*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 17: O + CH2 <=> H + HCO */
        phi_f = sc[2*npt+i]*sc[7*npt+i];
        k_f = k_f_s[16*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[13*npt+i];
        Kc = Kc_s[16*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[7*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;

        /*reaction 18: O + CH2(S) <=> H + HCO */
        phi_f = sc[2*npt+i]*sc[8*npt+i];
        k_f = k_f_s[17*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[13*npt+i];
        Kc = Kc_s[17*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[8*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;

        /*reaction 19: O + CH3 <=> H + CH2O */
        phi_f = sc[2*npt+i]*sc[9*npt+i];
        k_f = k_f_s[18*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[14*npt+i];
        Kc = Kc_s[18*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[9*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;

        /*reaction 20: O + CH4 <=> OH + CH3 */
        phi_f = sc[2*npt+i]*sc[10*npt+i];
        k_f = k_f_s[19*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[9*npt+i];
        Kc = Kc_s[19*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[9*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;

        /*reaction 21: O + HCO <=> OH + CO */
        phi_f = sc[2*npt+i]*sc[13*npt+i];
        k_f = k_f_s[20*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[11*npt+i];
        Kc = Kc_s[20*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[11*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;

        /*reaction 22: O + HCO <=> H + CO2 */
        phi_f = sc[2*npt+i]*sc[13*npt+i];
        k_f = k_f_s[21*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[12*npt+i];
        Kc = Kc_s[21*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;

        /*reaction 23: O + CH2O <=> OH + HCO */
        phi_f = sc[2*npt+i]*sc[14*npt+i];
        k_f = k_f_s[22*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[13*npt+i];
        Kc = Kc_s[22*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[13*npt+i] += qdot;
        wdot[14*npt+i] -= qdot;

        /*reaction 24: O + C2H4 <=> CH3 + HCO */
        phi_f = sc[2*npt+i]*sc[16*npt+i];
        k_f = k_f_s[23*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[9*npt+i]*sc[13*npt+i];
        Kc = Kc_s[23*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;
        wdot[13*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;

        /*reaction 25: O + C2H5 <=> CH3 + CH2O */
        phi_f = sc[2*npt+i]*sc[17*npt+i];
        k_f = k_f_s[24*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[9*npt+i]*sc[14*npt+i];
        Kc = Kc_s[24*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;
        wdot[14*npt+i] += qdot;
        wdot[17*npt+i] -= qdot;

        /*reaction 26: O + C2H6 <=> OH + C2H5 */
        phi_f = sc[2*npt+i]*sc[18*npt+i];
        k_f = k_f_s[25*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[17*npt+i];
        Kc = Kc_s[25*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[17*npt+i] += qdot;
        wdot[18*npt+i] -= qdot;

        /*reaction 27: O2 + CO <=> O + CO2 */
        phi_f = sc[3*npt+i]*sc[11*npt+i];
        k_f = k_f_s[26*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[12*npt+i];
        Kc = Kc_s[26*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[11*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;

        /*reaction 28: O2 + CH2O <=> HO2 + HCO */
        phi_f = sc[3*npt+i]*sc[14*npt+i];
        k_f = k_f_s[27*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[13*npt+i];
        Kc = Kc_s[27*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[13*npt+i] += qdot;
        wdot[14*npt+i] -= qdot;

        /*reaction 29: H + 2.000000 O2 <=> HO2 + O2 */
        phi_f = sc[1*npt+i]*pow(sc[3*npt+i], 2.000000);
        k_f = k_f_s[28*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[6*npt+i];
        Kc = Kc_s[28*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] -= 2.000000 * qdot;
        wdot[3*npt+i] += qdot;
        wdot[6*npt+i] += qdot;

        /*reaction 30: H + O2 + H2O <=> HO2 + H2O */
        phi_f = sc[1*npt+i]*sc[3*npt+i]*sc[5*npt+i];
        k_f = k_f_s[29*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[6*npt+i];
        Kc = Kc_s[29*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] -= qdot;
        wdot[5*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] += qdot;

        /*reaction 31: H + O2 + N2 <=> HO2 + N2 */
        phi_f = sc[1*npt+i]*sc[3*npt+i]*sc[19*npt+i];
        k_f = k_f_s[30*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[19*npt+i];
        Kc = Kc_s[30*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;
        wdot[19*npt+i] += qdot;

        /*reaction 32: H + O2 + AR <=> HO2 + AR */
        phi_f = sc[1*npt+i]*sc[3*npt+i]*sc[20*npt+i];
        k_f = k_f_s[31*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[20*npt+i];
        Kc = Kc_s[31*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[20*npt+i] -= qdot;
        wdot[20*npt+i] += qdot;

        /*reaction 33: H + O2 <=> O + OH */
        phi_f = sc[1*npt+i]*sc[3*npt+i];
        k_f = k_f_s[32*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[4*npt+i];
        Kc = Kc_s[32*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;

        /*reaction 34: 2.000000 H + H2 <=> 2.000000 H2 */
        phi_f = sc[0*npt+i]*pow(sc[1*npt+i], 2.000000);
        k_f = k_f_s[33*npt+i];
        q_f = phi_f * k_f;
        phi_r = pow(sc[0*npt+i], 2.000000);
        Kc = Kc_s[33*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[0*npt+i] += 2.000000 * qdot;
        wdot[1*npt+i] -= 2.000000 * qdot;

        /*reaction 35: 2.000000 H + H2O <=> H2 + H2O */
        phi_f = pow(sc[1*npt+i], 2.000000)*sc[5*npt+i];
        k_f = k_f_s[34*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[5*npt+i];
        Kc = Kc_s[34*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= 2.000000 * qdot;
        wdot[5*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;

        /*reaction 36: 2.000000 H + CO2 <=> H2 + CO2 */
        phi_f = pow(sc[1*npt+i], 2.000000)*sc[12*npt+i];
        k_f = k_f_s[35*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[12*npt+i];
        Kc = Kc_s[35*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= 2.000000 * qdot;
        wdot[12*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;

        /*reaction 37: H + HO2 <=> O2 + H2 */
        phi_f = sc[1*npt+i]*sc[6*npt+i];
        k_f = k_f_s[36*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[3*npt+i];
        Kc = Kc_s[36*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 38: H + HO2 <=> 2.000000 OH */
        phi_f = sc[1*npt+i]*sc[6*npt+i];
        k_f = k_f_s[37*npt+i];
        q_f = phi_f * k_f;
        phi_r = pow(sc[4*npt+i], 2.000000);
        Kc = Kc_s[37*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[4*npt+i] += 2.000000 * qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 39: H + CH4 <=> CH3 + H2 */
        phi_f = sc[1*npt+i]*sc[10*npt+i];
        k_f = k_f_s[38*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[9*npt+i];
        Kc = Kc_s[38*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;

        /*reaction 40: H + HCO <=> H2 + CO */
        phi_f = sc[1*npt+i]*sc[13*npt+i];
        k_f = k_f_s[39*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[11*npt+i];
        Kc = Kc_s[39*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;

        /*reaction 41: H + CH2O <=> HCO + H2 */
        phi_f = sc[1*npt+i]*sc[14*npt+i];
        k_f = k_f_s[40*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[13*npt+i];
        Kc = Kc_s[40*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;
        wdot[14*npt+i] -= qdot;

        /*reaction 42: H + CH3O <=> OH + CH3 */
        phi_f = sc[1*npt+i]*sc[15*npt+i];
        k_f = k_f_s[41*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[9*npt+i];
        Kc = Kc_s[41*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[9*npt+i] += qdot;
        wdot[15*npt+i] -= qdot;

        /*reaction 43: H + C2H6 <=> C2H5 + H2 */
        phi_f = sc[1*npt+i]*sc[18*npt+i];
        k_f = k_f_s[42*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[17*npt+i];
        Kc = Kc_s[42*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[17*npt+i] += qdot;
        wdot[18*npt+i] -= qdot;

        /*reaction 44: OH + H2 <=> H + H2O */
        phi_f = sc[0*npt+i]*sc[4*npt+i];
        k_f = k_f_s[43*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[5*npt+i];
        Kc = Kc_s[43*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;

        /*reaction 45: 2.000000 OH <=> O + H2O */
        phi_f = pow(sc[4*npt+i], 2.000000);
        k_f = k_f_s[44*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[5*npt+i];
        Kc = Kc_s[44*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[4*npt+i] -= 2.000000 * qdot;
        wdot[5*npt+i] += qdot;

        /*reaction 46: OH + HO2 <=> O2 + H2O */
        phi_f = sc[4*npt+i]*sc[6*npt+i];
        k_f = k_f_s[45*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[5*npt+i];
        Kc = Kc_s[45*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 47: OH + CH2 <=> H + CH2O */
        phi_f = sc[4*npt+i]*sc[7*npt+i];
        k_f = k_f_s[46*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[14*npt+i];
        Kc = Kc_s[46*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[7*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;

        /*reaction 48: OH + CH2(S) <=> H + CH2O */
        phi_f = sc[4*npt+i]*sc[8*npt+i];
        k_f = k_f_s[47*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[14*npt+i];
        Kc = Kc_s[47*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[8*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;

        /*reaction 49: OH + CH3 <=> CH2 + H2O */
        phi_f = sc[4*npt+i]*sc[9*npt+i];
        k_f = k_f_s[48*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[7*npt+i];
        Kc = Kc_s[48*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[7*npt+i] += qdot;
        wdot[9*npt+i] -= qdot;

        /*reaction 50: OH + CH3 <=> CH2(S) + H2O */
        phi_f = sc[4*npt+i]*sc[9*npt+i];
        k_f = k_f_s[49*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[8*npt+i];
        Kc = Kc_s[49*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[8*npt+i] += qdot;
        wdot[9*npt+i] -= qdot;
    }
}

void vcomp_wdot_51_84(int npt, amrex::Real *  wdot, amrex::Real *  mixture, amrex::Real *  sc,
		amrex::Real *  k_f_s, amrex::Real *  Kc_s,
		amrex::Real *  tc, amrex::Real *  invT, amrex::Real *  T)
{
    for (int i=0; i<npt; i++) {
        amrex::Real qdot, q_f, q_r, phi_f, phi_r, k_f, k_r, Kc;
        amrex::Real alpha;

        /*reaction 51: OH + CH4 <=> CH3 + H2O */
        phi_f = sc[4*npt+i]*sc[10*npt+i];
        k_f = k_f_s[50*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[9*npt+i];
        Kc = Kc_s[50*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[9*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;

        /*reaction 52: OH + CO <=> H + CO2 */
        phi_f = sc[4*npt+i]*sc[11*npt+i];
        k_f = k_f_s[51*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[12*npt+i];
        Kc = Kc_s[51*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[11*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;

        /*reaction 53: OH + HCO <=> H2O + CO */
        phi_f = sc[4*npt+i]*sc[13*npt+i];
        k_f = k_f_s[52*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[11*npt+i];
        Kc = Kc_s[52*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[11*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;

        /*reaction 54: OH + CH2O <=> HCO + H2O */
        phi_f = sc[4*npt+i]*sc[14*npt+i];
        k_f = k_f_s[53*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[13*npt+i];
        Kc = Kc_s[53*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[13*npt+i] += qdot;
        wdot[14*npt+i] -= qdot;

        /*reaction 55: OH + C2H6 <=> C2H5 + H2O */
        phi_f = sc[4*npt+i]*sc[18*npt+i];
        k_f = k_f_s[54*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[17*npt+i];
        Kc = Kc_s[54*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[17*npt+i] += qdot;
        wdot[18*npt+i] -= qdot;

        /*reaction 56: HO2 + CH2 <=> OH + CH2O */
        phi_f = sc[6*npt+i]*sc[7*npt+i];
        k_f = k_f_s[55*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[14*npt+i];
        Kc = Kc_s[55*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[7*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;

        /*reaction 57: HO2 + CH3 <=> O2 + CH4 */
        phi_f = sc[6*npt+i]*sc[9*npt+i];
        k_f = k_f_s[56*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[10*npt+i];
        Kc = Kc_s[56*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[9*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;

        /*reaction 58: HO2 + CH3 <=> OH + CH3O */
        phi_f = sc[6*npt+i]*sc[9*npt+i];
        k_f = k_f_s[57*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[15*npt+i];
        Kc = Kc_s[57*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[9*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;

        /*reaction 59: HO2 + CO <=> OH + CO2 */
        phi_f = sc[6*npt+i]*sc[11*npt+i];
        k_f = k_f_s[58*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[12*npt+i];
        Kc = Kc_s[58*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[11*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;

        /*reaction 60: CH2 + O2 <=> OH + HCO */
        phi_f = sc[3*npt+i]*sc[7*npt+i];
        k_f = k_f_s[59*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[13*npt+i];
        Kc = Kc_s[59*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;

        /*reaction 61: CH2 + H2 <=> H + CH3 */
        phi_f = sc[0*npt+i]*sc[7*npt+i];
        k_f = k_f_s[60*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[9*npt+i];
        Kc = Kc_s[60*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;

        /*reaction 62: CH2 + CH3 <=> H + C2H4 */
        phi_f = sc[7*npt+i]*sc[9*npt+i];
        k_f = k_f_s[61*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[16*npt+i];
        Kc = Kc_s[61*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;
        wdot[9*npt+i] -= qdot;
        wdot[16*npt+i] += qdot;

        /*reaction 63: CH2 + CH4 <=> 2.000000 CH3 */
        phi_f = sc[7*npt+i]*sc[10*npt+i];
        k_f = k_f_s[62*npt+i];
        q_f = phi_f * k_f;
        phi_r = pow(sc[9*npt+i], 2.000000);
        Kc = Kc_s[62*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[7*npt+i] -= qdot;
        wdot[9*npt+i] += 2.000000 * qdot;
        wdot[10*npt+i] -= qdot;

        /*reaction 64: CH2(S) + N2 <=> CH2 + N2 */
        phi_f = sc[8*npt+i]*sc[19*npt+i];
        k_f = k_f_s[63*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[19*npt+i];
        Kc = Kc_s[63*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[7*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[19*npt+i] -= qdot;
        wdot[19*npt+i] += qdot;

        /*reaction 65: CH2(S) + AR <=> CH2 + AR */
        phi_f = sc[8*npt+i]*sc[20*npt+i];
        k_f = k_f_s[64*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[20*npt+i];
        Kc = Kc_s[64*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[7*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[20*npt+i] -= qdot;
        wdot[20*npt+i] += qdot;

        /*reaction 66: CH2(S) + O2 <=> H + OH + CO */
        phi_f = sc[3*npt+i]*sc[8*npt+i];
        k_f = k_f_s[65*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[4*npt+i]*sc[11*npt+i];
        Kc = Kc_s[65*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;

        /*reaction 67: CH2(S) + O2 <=> CO + H2O */
        phi_f = sc[3*npt+i]*sc[8*npt+i];
        k_f = k_f_s[66*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[11*npt+i];
        Kc = Kc_s[66*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;

        /*reaction 68: CH2(S) + H2 <=> CH3 + H */
        phi_f = sc[0*npt+i]*sc[8*npt+i];
        k_f = k_f_s[67*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[9*npt+i];
        Kc = Kc_s[67*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;

        /*reaction 69: CH2(S) + H2O <=> CH2 + H2O */
        phi_f = sc[5*npt+i]*sc[8*npt+i];
        k_f = k_f_s[68*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[7*npt+i];
        Kc = Kc_s[68*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[7*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;

        /*reaction 70: CH2(S) + CH3 <=> H + C2H4 */
        phi_f = sc[8*npt+i]*sc[9*npt+i];
        k_f = k_f_s[69*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[16*npt+i];
        Kc = Kc_s[69*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[9*npt+i] -= qdot;
        wdot[16*npt+i] += qdot;

        /*reaction 71: CH2(S) + CH4 <=> 2.000000 CH3 */
        phi_f = sc[8*npt+i]*sc[10*npt+i];
        k_f = k_f_s[70*npt+i];
        q_f = phi_f * k_f;
        phi_r = pow(sc[9*npt+i], 2.000000);
        Kc = Kc_s[70*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[8*npt+i] -= qdot;
        wdot[9*npt+i] += 2.000000 * qdot;
        wdot[10*npt+i] -= qdot;

        /*reaction 72: CH2(S) + CO <=> CH2 + CO */
        phi_f = sc[8*npt+i]*sc[11*npt+i];
        k_f = k_f_s[71*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[11*npt+i];
        Kc = Kc_s[71*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[7*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[11*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;

        /*reaction 73: CH2(S) + CO2 <=> CH2 + CO2 */
        phi_f = sc[8*npt+i]*sc[12*npt+i];
        k_f = k_f_s[72*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[12*npt+i];
        Kc = Kc_s[72*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[7*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[12*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;

        /*reaction 74: CH2(S) + CO2 <=> CO + CH2O */
        phi_f = sc[8*npt+i]*sc[12*npt+i];
        k_f = k_f_s[73*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[11*npt+i]*sc[14*npt+i];
        Kc = Kc_s[73*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[8*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;
        wdot[12*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;

        /*reaction 75: CH3 + O2 <=> O + CH3O */
        phi_f = sc[3*npt+i]*sc[9*npt+i];
        k_f = k_f_s[74*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[15*npt+i];
        Kc = Kc_s[74*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[9*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;

        /*reaction 76: CH3 + O2 <=> OH + CH2O */
        phi_f = sc[3*npt+i]*sc[9*npt+i];
        k_f = k_f_s[75*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[14*npt+i];
        Kc = Kc_s[75*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[9*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;

        /*reaction 77: 2.000000 CH3 <=> H + C2H5 */
        phi_f = pow(sc[9*npt+i], 2.000000);
        k_f = k_f_s[76*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[17*npt+i];
        Kc = Kc_s[76*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[9*npt+i] -= 2.000000 * qdot;
        wdot[17*npt+i] += qdot;

        /*reaction 78: CH3 + HCO <=> CH4 + CO */
        phi_f = sc[9*npt+i]*sc[13*npt+i];
        k_f = k_f_s[77*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[10*npt+i]*sc[11*npt+i];
        Kc = Kc_s[77*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[9*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;
        wdot[11*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;

        /*reaction 79: CH3 + CH2O <=> HCO + CH4 */
        phi_f = sc[9*npt+i]*sc[14*npt+i];
        k_f = k_f_s[78*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[10*npt+i]*sc[13*npt+i];
        Kc = Kc_s[78*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[9*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;
        wdot[13*npt+i] += qdot;
        wdot[14*npt+i] -= qdot;

        /*reaction 80: CH3 + C2H6 <=> C2H5 + CH4 */
        phi_f = sc[9*npt+i]*sc[18*npt+i];
        k_f = k_f_s[79*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[10*npt+i]*sc[17*npt+i];
        Kc = Kc_s[79*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[9*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;
        wdot[17*npt+i] += qdot;
        wdot[18*npt+i] -= qdot;

        /*reaction 81: HCO + H2O <=> H + CO + H2O */
        phi_f = sc[5*npt+i]*sc[13*npt+i];
        k_f = k_f_s[80*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[5*npt+i]*sc[11*npt+i];
        Kc = Kc_s[80*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[11*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;

        /*reaction 82: HCO + O2 <=> HO2 + CO */
        phi_f = sc[3*npt+i]*sc[13*npt+i];
        k_f = k_f_s[81*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[11*npt+i];
        Kc = Kc_s[81*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[11*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;

        /*reaction 83: CH3O + O2 <=> HO2 + CH2O */
        phi_f = sc[3*npt+i]*sc[15*npt+i];
        k_f = k_f_s[82*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[14*npt+i];
        Kc = Kc_s[82*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[14*npt+i] += qdot;
        wdot[15*npt+i] -= qdot;

        /*reaction 84: C2H5 + O2 <=> HO2 + C2H4 */
        phi_f = sc[3*npt+i]*sc[17*npt+i];
        k_f = k_f_s[83*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[16*npt+i];
        Kc = Kc_s[83*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[16*npt+i] += qdot;
        wdot[17*npt+i] -= qdot;
    }
}


static amrex::Real T_save = -1;
#ifdef _OPENMP
#pragma omp threadprivate(T_save)
#endif

static amrex::Real k_f_save[84];
#ifdef _OPENMP
#pragma omp threadprivate(k_f_save)
#endif

static amrex::Real Kc_save[84];
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

    amrex::Real qdot, q_f[84], q_r[84];
    amrex::Real sc_qss[1];
    comp_qfqr(q_f, q_r, sc, sc_qss, tc, invT);

    for (int i = 0; i < 21; ++i) {
        wdot[i] = 0.0;
    }

    qdot = q_f[0]-q_r[0];
    wdot[1] -= qdot;
    wdot[7] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[1]-q_r[1];
    wdot[1] -= qdot;
    wdot[9] -= qdot;
    wdot[10] += qdot;

    qdot = q_f[2]-q_r[2];
    wdot[1] -= qdot;
    wdot[13] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[3]-q_r[3];
    wdot[1] -= qdot;
    wdot[14] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[4]-q_r[4];
    wdot[1] -= qdot;
    wdot[16] -= qdot;
    wdot[17] += qdot;

    qdot = q_f[5]-q_r[5];
    wdot[1] -= qdot;
    wdot[17] -= qdot;
    wdot[18] += qdot;

    qdot = q_f[6]-q_r[6];
    wdot[0] -= qdot;
    wdot[11] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[7]-q_r[7];
    wdot[9] -= 2.000000 * qdot;
    wdot[18] += qdot;

    qdot = q_f[8]-q_r[8];
    wdot[1] -= qdot;
    wdot[2] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[9]-q_r[9];
    wdot[2] -= qdot;
    wdot[11] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[10]-q_r[10];
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[6] += qdot;

    qdot = q_f[11]-q_r[11];
    wdot[0] += qdot;
    wdot[1] -= 2.000000 * qdot;

    qdot = q_f[12]-q_r[12];
    wdot[1] -= qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[13]-q_r[13];
    wdot[1] += qdot;
    wdot[11] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[14]-q_r[14];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[15]-q_r[15];
    wdot[2] -= qdot;
    wdot[3] += qdot;
    wdot[4] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[16]-q_r[16];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[7] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[17]-q_r[17];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[8] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[18]-q_r[18];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[9] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[19]-q_r[19];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[9] += qdot;
    wdot[10] -= qdot;

    qdot = q_f[20]-q_r[20];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[11] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[21]-q_r[21];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[12] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[22]-q_r[22];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[13] += qdot;
    wdot[14] -= qdot;

    qdot = q_f[23]-q_r[23];
    wdot[2] -= qdot;
    wdot[9] += qdot;
    wdot[13] += qdot;
    wdot[16] -= qdot;

    qdot = q_f[24]-q_r[24];
    wdot[2] -= qdot;
    wdot[9] += qdot;
    wdot[14] += qdot;
    wdot[17] -= qdot;

    qdot = q_f[25]-q_r[25];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[17] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[26]-q_r[26];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[11] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[27]-q_r[27];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[13] += qdot;
    wdot[14] -= qdot;

    qdot = q_f[28]-q_r[28];
    wdot[1] -= qdot;
    wdot[3] -= 2.000000 * qdot;
    wdot[3] += qdot;
    wdot[6] += qdot;

    qdot = q_f[29]-q_r[29];
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[5] -= qdot;
    wdot[5] += qdot;
    wdot[6] += qdot;

    qdot = q_f[30]-q_r[30];
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[19] -= qdot;
    wdot[19] += qdot;

    qdot = q_f[31]-q_r[31];
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[20] -= qdot;
    wdot[20] += qdot;

    qdot = q_f[32]-q_r[32];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[33]-q_r[33];
    wdot[0] -= qdot;
    wdot[0] += 2.000000 * qdot;
    wdot[1] -= 2.000000 * qdot;

    qdot = q_f[34]-q_r[34];
    wdot[0] += qdot;
    wdot[1] -= 2.000000 * qdot;
    wdot[5] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[35]-q_r[35];
    wdot[0] += qdot;
    wdot[1] -= 2.000000 * qdot;
    wdot[12] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[36]-q_r[36];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[3] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[37]-q_r[37];
    wdot[1] -= qdot;
    wdot[4] += 2.000000 * qdot;
    wdot[6] -= qdot;

    qdot = q_f[38]-q_r[38];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[9] += qdot;
    wdot[10] -= qdot;

    qdot = q_f[39]-q_r[39];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[11] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[40]-q_r[40];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[13] += qdot;
    wdot[14] -= qdot;

    qdot = q_f[41]-q_r[41];
    wdot[1] -= qdot;
    wdot[4] += qdot;
    wdot[9] += qdot;
    wdot[15] -= qdot;

    qdot = q_f[42]-q_r[42];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[17] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[43]-q_r[43];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[44]-q_r[44];
    wdot[2] += qdot;
    wdot[4] -= 2.000000 * qdot;
    wdot[5] += qdot;

    qdot = q_f[45]-q_r[45];
    wdot[3] += qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[46]-q_r[46];
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[7] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[47]-q_r[47];
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[8] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[48]-q_r[48];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[7] += qdot;
    wdot[9] -= qdot;

    qdot = q_f[49]-q_r[49];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[8] += qdot;
    wdot[9] -= qdot;

    qdot = q_f[50]-q_r[50];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[9] += qdot;
    wdot[10] -= qdot;

    qdot = q_f[51]-q_r[51];
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[11] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[52]-q_r[52];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[11] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[53]-q_r[53];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[13] += qdot;
    wdot[14] -= qdot;

    qdot = q_f[54]-q_r[54];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[17] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[55]-q_r[55];
    wdot[4] += qdot;
    wdot[6] -= qdot;
    wdot[7] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[56]-q_r[56];
    wdot[3] += qdot;
    wdot[6] -= qdot;
    wdot[9] -= qdot;
    wdot[10] += qdot;

    qdot = q_f[57]-q_r[57];
    wdot[4] += qdot;
    wdot[6] -= qdot;
    wdot[9] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[58]-q_r[58];
    wdot[4] += qdot;
    wdot[6] -= qdot;
    wdot[11] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[59]-q_r[59];
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[7] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[60]-q_r[60];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[7] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[61]-q_r[61];
    wdot[1] += qdot;
    wdot[7] -= qdot;
    wdot[9] -= qdot;
    wdot[16] += qdot;

    qdot = q_f[62]-q_r[62];
    wdot[7] -= qdot;
    wdot[9] += 2.000000 * qdot;
    wdot[10] -= qdot;

    qdot = q_f[63]-q_r[63];
    wdot[7] += qdot;
    wdot[8] -= qdot;
    wdot[19] -= qdot;
    wdot[19] += qdot;

    qdot = q_f[64]-q_r[64];
    wdot[7] += qdot;
    wdot[8] -= qdot;
    wdot[20] -= qdot;
    wdot[20] += qdot;

    qdot = q_f[65]-q_r[65];
    wdot[1] += qdot;
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[8] -= qdot;
    wdot[11] += qdot;

    qdot = q_f[66]-q_r[66];
    wdot[3] -= qdot;
    wdot[5] += qdot;
    wdot[8] -= qdot;
    wdot[11] += qdot;

    qdot = q_f[67]-q_r[67];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[8] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[68]-q_r[68];
    wdot[5] -= qdot;
    wdot[5] += qdot;
    wdot[7] += qdot;
    wdot[8] -= qdot;

    qdot = q_f[69]-q_r[69];
    wdot[1] += qdot;
    wdot[8] -= qdot;
    wdot[9] -= qdot;
    wdot[16] += qdot;

    qdot = q_f[70]-q_r[70];
    wdot[8] -= qdot;
    wdot[9] += 2.000000 * qdot;
    wdot[10] -= qdot;

    qdot = q_f[71]-q_r[71];
    wdot[7] += qdot;
    wdot[8] -= qdot;
    wdot[11] -= qdot;
    wdot[11] += qdot;

    qdot = q_f[72]-q_r[72];
    wdot[7] += qdot;
    wdot[8] -= qdot;
    wdot[12] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[73]-q_r[73];
    wdot[8] -= qdot;
    wdot[11] += qdot;
    wdot[12] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[74]-q_r[74];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[9] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[75]-q_r[75];
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[9] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[76]-q_r[76];
    wdot[1] += qdot;
    wdot[9] -= 2.000000 * qdot;
    wdot[17] += qdot;

    qdot = q_f[77]-q_r[77];
    wdot[9] -= qdot;
    wdot[10] += qdot;
    wdot[11] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[78]-q_r[78];
    wdot[9] -= qdot;
    wdot[10] += qdot;
    wdot[13] += qdot;
    wdot[14] -= qdot;

    qdot = q_f[79]-q_r[79];
    wdot[9] -= qdot;
    wdot[10] += qdot;
    wdot[17] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[80]-q_r[80];
    wdot[1] += qdot;
    wdot[5] -= qdot;
    wdot[5] += qdot;
    wdot[11] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[81]-q_r[81];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[11] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[82]-q_r[82];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[14] += qdot;
    wdot[15] -= qdot;

    qdot = q_f[83]-q_r[83];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[16] += qdot;
    wdot[17] -= qdot;

    return;
}

void comp_k_f(amrex::Real *  tc, amrex::Real invT, amrex::Real *  k_f)
{
    for (int i=0; i<84; ++i) {
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

    Kc[0] = g_RT[1] + g_RT[7] - g_RT[9];
    Kc[1] = g_RT[1] + g_RT[9] - g_RT[10];
    Kc[2] = g_RT[1] + g_RT[13] - g_RT[14];
    Kc[3] = g_RT[1] + g_RT[14] - g_RT[15];
    Kc[4] = g_RT[1] + g_RT[16] - g_RT[17];
    Kc[5] = g_RT[1] + g_RT[17] - g_RT[18];
    Kc[6] = g_RT[0] + g_RT[11] - g_RT[14];
    Kc[7] = 2.000000*g_RT[9] - g_RT[18];
    Kc[8] = g_RT[1] + g_RT[2] - g_RT[4];
    Kc[9] = g_RT[2] + g_RT[11] - g_RT[12];
    Kc[10] = g_RT[1] + g_RT[3] - g_RT[6];
    Kc[11] = -g_RT[0] + 2.000000*g_RT[1];
    Kc[12] = g_RT[1] + g_RT[4] - g_RT[5];
    Kc[13] = -g_RT[1] - g_RT[11] + g_RT[13];
    Kc[14] = g_RT[0] - g_RT[1] + g_RT[2] - g_RT[4];
    Kc[15] = g_RT[2] - g_RT[3] - g_RT[4] + g_RT[6];
    Kc[16] = -g_RT[1] + g_RT[2] + g_RT[7] - g_RT[13];
    Kc[17] = -g_RT[1] + g_RT[2] + g_RT[8] - g_RT[13];
    Kc[18] = -g_RT[1] + g_RT[2] + g_RT[9] - g_RT[14];
    Kc[19] = g_RT[2] - g_RT[4] - g_RT[9] + g_RT[10];
    Kc[20] = g_RT[2] - g_RT[4] - g_RT[11] + g_RT[13];
    Kc[21] = -g_RT[1] + g_RT[2] - g_RT[12] + g_RT[13];
    Kc[22] = g_RT[2] - g_RT[4] - g_RT[13] + g_RT[14];
    Kc[23] = g_RT[2] - g_RT[9] - g_RT[13] + g_RT[16];
    Kc[24] = g_RT[2] - g_RT[9] - g_RT[14] + g_RT[17];
    Kc[25] = g_RT[2] - g_RT[4] - g_RT[17] + g_RT[18];
    Kc[26] = -g_RT[2] + g_RT[3] + g_RT[11] - g_RT[12];
    Kc[27] = g_RT[3] - g_RT[6] - g_RT[13] + g_RT[14];
    Kc[28] = g_RT[1] + 2.000000*g_RT[3] - g_RT[3] - g_RT[6];
    Kc[29] = g_RT[1] + g_RT[3] + g_RT[5] - g_RT[5] - g_RT[6];
    Kc[30] = g_RT[1] + g_RT[3] - g_RT[6] + g_RT[19] - g_RT[19];
    Kc[31] = g_RT[1] + g_RT[3] - g_RT[6] + g_RT[20] - g_RT[20];
    Kc[32] = g_RT[1] - g_RT[2] + g_RT[3] - g_RT[4];
    Kc[33] = g_RT[0] - 2.000000*g_RT[0] + 2.000000*g_RT[1];
    Kc[34] = -g_RT[0] + 2.000000*g_RT[1] + g_RT[5] - g_RT[5];
    Kc[35] = -g_RT[0] + 2.000000*g_RT[1] + g_RT[12] - g_RT[12];
    Kc[36] = -g_RT[0] + g_RT[1] - g_RT[3] + g_RT[6];
    Kc[37] = g_RT[1] - 2.000000*g_RT[4] + g_RT[6];
    Kc[38] = -g_RT[0] + g_RT[1] - g_RT[9] + g_RT[10];
    Kc[39] = -g_RT[0] + g_RT[1] - g_RT[11] + g_RT[13];
    Kc[40] = -g_RT[0] + g_RT[1] - g_RT[13] + g_RT[14];
    Kc[41] = g_RT[1] - g_RT[4] - g_RT[9] + g_RT[15];
    Kc[42] = -g_RT[0] + g_RT[1] - g_RT[17] + g_RT[18];
    Kc[43] = g_RT[0] - g_RT[1] + g_RT[4] - g_RT[5];
    Kc[44] = -g_RT[2] + 2.000000*g_RT[4] - g_RT[5];
    Kc[45] = -g_RT[3] + g_RT[4] - g_RT[5] + g_RT[6];
    Kc[46] = -g_RT[1] + g_RT[4] + g_RT[7] - g_RT[14];
    Kc[47] = -g_RT[1] + g_RT[4] + g_RT[8] - g_RT[14];
    Kc[48] = g_RT[4] - g_RT[5] - g_RT[7] + g_RT[9];
    Kc[49] = g_RT[4] - g_RT[5] - g_RT[8] + g_RT[9];
    Kc[50] = g_RT[4] - g_RT[5] - g_RT[9] + g_RT[10];
    Kc[51] = -g_RT[1] + g_RT[4] + g_RT[11] - g_RT[12];
    Kc[52] = g_RT[4] - g_RT[5] - g_RT[11] + g_RT[13];
    Kc[53] = g_RT[4] - g_RT[5] - g_RT[13] + g_RT[14];
    Kc[54] = g_RT[4] - g_RT[5] - g_RT[17] + g_RT[18];
    Kc[55] = -g_RT[4] + g_RT[6] + g_RT[7] - g_RT[14];
    Kc[56] = -g_RT[3] + g_RT[6] + g_RT[9] - g_RT[10];
    Kc[57] = -g_RT[4] + g_RT[6] + g_RT[9] - g_RT[15];
    Kc[58] = -g_RT[4] + g_RT[6] + g_RT[11] - g_RT[12];
    Kc[59] = g_RT[3] - g_RT[4] + g_RT[7] - g_RT[13];
    Kc[60] = g_RT[0] - g_RT[1] + g_RT[7] - g_RT[9];
    Kc[61] = -g_RT[1] + g_RT[7] + g_RT[9] - g_RT[16];
    Kc[62] = g_RT[7] - 2.000000*g_RT[9] + g_RT[10];
    Kc[63] = -g_RT[7] + g_RT[8] + g_RT[19] - g_RT[19];
    Kc[64] = -g_RT[7] + g_RT[8] + g_RT[20] - g_RT[20];
    Kc[65] = -g_RT[1] + g_RT[3] - g_RT[4] + g_RT[8] - g_RT[11];
    Kc[66] = g_RT[3] - g_RT[5] + g_RT[8] - g_RT[11];
    Kc[67] = g_RT[0] - g_RT[1] + g_RT[8] - g_RT[9];
    Kc[68] = g_RT[5] - g_RT[5] - g_RT[7] + g_RT[8];
    Kc[69] = -g_RT[1] + g_RT[8] + g_RT[9] - g_RT[16];
    Kc[70] = g_RT[8] - 2.000000*g_RT[9] + g_RT[10];
    Kc[71] = -g_RT[7] + g_RT[8] + g_RT[11] - g_RT[11];
    Kc[72] = -g_RT[7] + g_RT[8] + g_RT[12] - g_RT[12];
    Kc[73] = g_RT[8] - g_RT[11] + g_RT[12] - g_RT[14];
    Kc[74] = -g_RT[2] + g_RT[3] + g_RT[9] - g_RT[15];
    Kc[75] = g_RT[3] - g_RT[4] + g_RT[9] - g_RT[14];
    Kc[76] = -g_RT[1] + 2.000000*g_RT[9] - g_RT[17];
    Kc[77] = g_RT[9] - g_RT[10] - g_RT[11] + g_RT[13];
    Kc[78] = g_RT[9] - g_RT[10] - g_RT[13] + g_RT[14];
    Kc[79] = g_RT[9] - g_RT[10] - g_RT[17] + g_RT[18];
    Kc[80] = -g_RT[1] + g_RT[5] - g_RT[5] - g_RT[11] + g_RT[13];
    Kc[81] = g_RT[3] - g_RT[6] - g_RT[11] + g_RT[13];
    Kc[82] = g_RT[3] - g_RT[6] - g_RT[14] + g_RT[15];
    Kc[83] = g_RT[3] - g_RT[6] - g_RT[16] + g_RT[17];

    for (int i=0; i<84; ++i) {
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
    Kc[7] *= refCinv;
    Kc[8] *= refCinv;
    Kc[9] *= refCinv;
    Kc[10] *= refCinv;
    Kc[11] *= refCinv;
    Kc[12] *= refCinv;
    Kc[13] *= refC;
    Kc[28] *= refCinv;
    Kc[29] *= refCinv;
    Kc[30] *= refCinv;
    Kc[31] *= refCinv;
    Kc[33] *= refCinv;
    Kc[34] *= refCinv;
    Kc[35] *= refCinv;
    Kc[65] *= refC;
    Kc[80] *= refC;

    return;
}

void comp_qfqr(amrex::Real *  qf, amrex::Real *  qr, amrex::Real *  sc, amrex::Real * qss_sc, amrex::Real *  tc, amrex::Real invT)
{

    /*reaction 1: H + CH2 (+M) <=> CH3 (+M) */
    qf[0] = sc[1]*sc[7];
    qr[0] = sc[9];

    /*reaction 2: H + CH3 (+M) <=> CH4 (+M) */
    qf[1] = sc[1]*sc[9];
    qr[1] = sc[10];

    /*reaction 3: H + HCO (+M) <=> CH2O (+M) */
    qf[2] = sc[1]*sc[13];
    qr[2] = sc[14];

    /*reaction 4: H + CH2O (+M) <=> CH3O (+M) */
    qf[3] = sc[1]*sc[14];
    qr[3] = sc[15];

    /*reaction 5: H + C2H4 (+M) <=> C2H5 (+M) */
    qf[4] = sc[1]*sc[16];
    qr[4] = sc[17];

    /*reaction 6: H + C2H5 (+M) <=> C2H6 (+M) */
    qf[5] = sc[1]*sc[17];
    qr[5] = sc[18];

    /*reaction 7: H2 + CO (+M) <=> CH2O (+M) */
    qf[6] = sc[0]*sc[11];
    qr[6] = sc[14];

    /*reaction 8: 2.000000 CH3 (+M) <=> C2H6 (+M) */
    qf[7] = pow(sc[9], 2.000000);
    qr[7] = sc[18];

    /*reaction 9: O + H + M <=> OH + M */
    qf[8] = sc[1]*sc[2];
    qr[8] = sc[4];

    /*reaction 10: O + CO + M <=> CO2 + M */
    qf[9] = sc[2]*sc[11];
    qr[9] = sc[12];

    /*reaction 11: H + O2 + M <=> HO2 + M */
    qf[10] = sc[1]*sc[3];
    qr[10] = sc[6];

    /*reaction 12: 2.000000 H + M <=> H2 + M */
    qf[11] = pow(sc[1], 2.000000);
    qr[11] = sc[0];

    /*reaction 13: H + OH + M <=> H2O + M */
    qf[12] = sc[1]*sc[4];
    qr[12] = sc[5];

    /*reaction 14: HCO + M <=> H + CO + M */
    qf[13] = sc[13];
    qr[13] = sc[1]*sc[11];

    /*reaction 15: O + H2 <=> H + OH */
    qf[14] = sc[0]*sc[2];
    qr[14] = sc[1]*sc[4];

    /*reaction 16: O + HO2 <=> OH + O2 */
    qf[15] = sc[2]*sc[6];
    qr[15] = sc[3]*sc[4];

    /*reaction 17: O + CH2 <=> H + HCO */
    qf[16] = sc[2]*sc[7];
    qr[16] = sc[1]*sc[13];

    /*reaction 18: O + CH2(S) <=> H + HCO */
    qf[17] = sc[2]*sc[8];
    qr[17] = sc[1]*sc[13];

    /*reaction 19: O + CH3 <=> H + CH2O */
    qf[18] = sc[2]*sc[9];
    qr[18] = sc[1]*sc[14];

    /*reaction 20: O + CH4 <=> OH + CH3 */
    qf[19] = sc[2]*sc[10];
    qr[19] = sc[4]*sc[9];

    /*reaction 21: O + HCO <=> OH + CO */
    qf[20] = sc[2]*sc[13];
    qr[20] = sc[4]*sc[11];

    /*reaction 22: O + HCO <=> H + CO2 */
    qf[21] = sc[2]*sc[13];
    qr[21] = sc[1]*sc[12];

    /*reaction 23: O + CH2O <=> OH + HCO */
    qf[22] = sc[2]*sc[14];
    qr[22] = sc[4]*sc[13];

    /*reaction 24: O + C2H4 <=> CH3 + HCO */
    qf[23] = sc[2]*sc[16];
    qr[23] = sc[9]*sc[13];

    /*reaction 25: O + C2H5 <=> CH3 + CH2O */
    qf[24] = sc[2]*sc[17];
    qr[24] = sc[9]*sc[14];

    /*reaction 26: O + C2H6 <=> OH + C2H5 */
    qf[25] = sc[2]*sc[18];
    qr[25] = sc[4]*sc[17];

    /*reaction 27: O2 + CO <=> O + CO2 */
    qf[26] = sc[3]*sc[11];
    qr[26] = sc[2]*sc[12];

    /*reaction 28: O2 + CH2O <=> HO2 + HCO */
    qf[27] = sc[3]*sc[14];
    qr[27] = sc[6]*sc[13];

    /*reaction 29: H + 2.000000 O2 <=> HO2 + O2 */
    qf[28] = sc[1]*pow(sc[3], 2.000000);
    qr[28] = sc[3]*sc[6];

    /*reaction 30: H + O2 + H2O <=> HO2 + H2O */
    qf[29] = sc[1]*sc[3]*sc[5];
    qr[29] = sc[5]*sc[6];

    /*reaction 31: H + O2 + N2 <=> HO2 + N2 */
    qf[30] = sc[1]*sc[3]*sc[19];
    qr[30] = sc[6]*sc[19];

    /*reaction 32: H + O2 + AR <=> HO2 + AR */
    qf[31] = sc[1]*sc[3]*sc[20];
    qr[31] = sc[6]*sc[20];

    /*reaction 33: H + O2 <=> O + OH */
    qf[32] = sc[1]*sc[3];
    qr[32] = sc[2]*sc[4];

    /*reaction 34: 2.000000 H + H2 <=> 2.000000 H2 */
    qf[33] = sc[0]*pow(sc[1], 2.000000);
    qr[33] = pow(sc[0], 2.000000);

    /*reaction 35: 2.000000 H + H2O <=> H2 + H2O */
    qf[34] = pow(sc[1], 2.000000)*sc[5];
    qr[34] = sc[0]*sc[5];

    /*reaction 36: 2.000000 H + CO2 <=> H2 + CO2 */
    qf[35] = pow(sc[1], 2.000000)*sc[12];
    qr[35] = sc[0]*sc[12];

    /*reaction 37: H + HO2 <=> O2 + H2 */
    qf[36] = sc[1]*sc[6];
    qr[36] = sc[0]*sc[3];

    /*reaction 38: H + HO2 <=> 2.000000 OH */
    qf[37] = sc[1]*sc[6];
    qr[37] = pow(sc[4], 2.000000);

    /*reaction 39: H + CH4 <=> CH3 + H2 */
    qf[38] = sc[1]*sc[10];
    qr[38] = sc[0]*sc[9];

    /*reaction 40: H + HCO <=> H2 + CO */
    qf[39] = sc[1]*sc[13];
    qr[39] = sc[0]*sc[11];

    /*reaction 41: H + CH2O <=> HCO + H2 */
    qf[40] = sc[1]*sc[14];
    qr[40] = sc[0]*sc[13];

    /*reaction 42: H + CH3O <=> OH + CH3 */
    qf[41] = sc[1]*sc[15];
    qr[41] = sc[4]*sc[9];

    /*reaction 43: H + C2H6 <=> C2H5 + H2 */
    qf[42] = sc[1]*sc[18];
    qr[42] = sc[0]*sc[17];

    /*reaction 44: OH + H2 <=> H + H2O */
    qf[43] = sc[0]*sc[4];
    qr[43] = sc[1]*sc[5];

    /*reaction 45: 2.000000 OH <=> O + H2O */
    qf[44] = pow(sc[4], 2.000000);
    qr[44] = sc[2]*sc[5];

    /*reaction 46: OH + HO2 <=> O2 + H2O */
    qf[45] = sc[4]*sc[6];
    qr[45] = sc[3]*sc[5];

    /*reaction 47: OH + CH2 <=> H + CH2O */
    qf[46] = sc[4]*sc[7];
    qr[46] = sc[1]*sc[14];

    /*reaction 48: OH + CH2(S) <=> H + CH2O */
    qf[47] = sc[4]*sc[8];
    qr[47] = sc[1]*sc[14];

    /*reaction 49: OH + CH3 <=> CH2 + H2O */
    qf[48] = sc[4]*sc[9];
    qr[48] = sc[5]*sc[7];

    /*reaction 50: OH + CH3 <=> CH2(S) + H2O */
    qf[49] = sc[4]*sc[9];
    qr[49] = sc[5]*sc[8];

    /*reaction 51: OH + CH4 <=> CH3 + H2O */
    qf[50] = sc[4]*sc[10];
    qr[50] = sc[5]*sc[9];

    /*reaction 52: OH + CO <=> H + CO2 */
    qf[51] = sc[4]*sc[11];
    qr[51] = sc[1]*sc[12];

    /*reaction 53: OH + HCO <=> H2O + CO */
    qf[52] = sc[4]*sc[13];
    qr[52] = sc[5]*sc[11];

    /*reaction 54: OH + CH2O <=> HCO + H2O */
    qf[53] = sc[4]*sc[14];
    qr[53] = sc[5]*sc[13];

    /*reaction 55: OH + C2H6 <=> C2H5 + H2O */
    qf[54] = sc[4]*sc[18];
    qr[54] = sc[5]*sc[17];

    /*reaction 56: HO2 + CH2 <=> OH + CH2O */
    qf[55] = sc[6]*sc[7];
    qr[55] = sc[4]*sc[14];

    /*reaction 57: HO2 + CH3 <=> O2 + CH4 */
    qf[56] = sc[6]*sc[9];
    qr[56] = sc[3]*sc[10];

    /*reaction 58: HO2 + CH3 <=> OH + CH3O */
    qf[57] = sc[6]*sc[9];
    qr[57] = sc[4]*sc[15];

    /*reaction 59: HO2 + CO <=> OH + CO2 */
    qf[58] = sc[6]*sc[11];
    qr[58] = sc[4]*sc[12];

    /*reaction 60: CH2 + O2 <=> OH + HCO */
    qf[59] = sc[3]*sc[7];
    qr[59] = sc[4]*sc[13];

    /*reaction 61: CH2 + H2 <=> H + CH3 */
    qf[60] = sc[0]*sc[7];
    qr[60] = sc[1]*sc[9];

    /*reaction 62: CH2 + CH3 <=> H + C2H4 */
    qf[61] = sc[7]*sc[9];
    qr[61] = sc[1]*sc[16];

    /*reaction 63: CH2 + CH4 <=> 2.000000 CH3 */
    qf[62] = sc[7]*sc[10];
    qr[62] = pow(sc[9], 2.000000);

    /*reaction 64: CH2(S) + N2 <=> CH2 + N2 */
    qf[63] = sc[8]*sc[19];
    qr[63] = sc[7]*sc[19];

    /*reaction 65: CH2(S) + AR <=> CH2 + AR */
    qf[64] = sc[8]*sc[20];
    qr[64] = sc[7]*sc[20];

    /*reaction 66: CH2(S) + O2 <=> H + OH + CO */
    qf[65] = sc[3]*sc[8];
    qr[65] = sc[1]*sc[4]*sc[11];

    /*reaction 67: CH2(S) + O2 <=> CO + H2O */
    qf[66] = sc[3]*sc[8];
    qr[66] = sc[5]*sc[11];

    /*reaction 68: CH2(S) + H2 <=> CH3 + H */
    qf[67] = sc[0]*sc[8];
    qr[67] = sc[1]*sc[9];

    /*reaction 69: CH2(S) + H2O <=> CH2 + H2O */
    qf[68] = sc[5]*sc[8];
    qr[68] = sc[5]*sc[7];

    /*reaction 70: CH2(S) + CH3 <=> H + C2H4 */
    qf[69] = sc[8]*sc[9];
    qr[69] = sc[1]*sc[16];

    /*reaction 71: CH2(S) + CH4 <=> 2.000000 CH3 */
    qf[70] = sc[8]*sc[10];
    qr[70] = pow(sc[9], 2.000000);

    /*reaction 72: CH2(S) + CO <=> CH2 + CO */
    qf[71] = sc[8]*sc[11];
    qr[71] = sc[7]*sc[11];

    /*reaction 73: CH2(S) + CO2 <=> CH2 + CO2 */
    qf[72] = sc[8]*sc[12];
    qr[72] = sc[7]*sc[12];

    /*reaction 74: CH2(S) + CO2 <=> CO + CH2O */
    qf[73] = sc[8]*sc[12];
    qr[73] = sc[11]*sc[14];

    /*reaction 75: CH3 + O2 <=> O + CH3O */
    qf[74] = sc[3]*sc[9];
    qr[74] = sc[2]*sc[15];

    /*reaction 76: CH3 + O2 <=> OH + CH2O */
    qf[75] = sc[3]*sc[9];
    qr[75] = sc[4]*sc[14];

    /*reaction 77: 2.000000 CH3 <=> H + C2H5 */
    qf[76] = pow(sc[9], 2.000000);
    qr[76] = sc[1]*sc[17];

    /*reaction 78: CH3 + HCO <=> CH4 + CO */
    qf[77] = sc[9]*sc[13];
    qr[77] = sc[10]*sc[11];

    /*reaction 79: CH3 + CH2O <=> HCO + CH4 */
    qf[78] = sc[9]*sc[14];
    qr[78] = sc[10]*sc[13];

    /*reaction 80: CH3 + C2H6 <=> C2H5 + CH4 */
    qf[79] = sc[9]*sc[18];
    qr[79] = sc[10]*sc[17];

    /*reaction 81: HCO + H2O <=> H + CO + H2O */
    qf[80] = sc[5]*sc[13];
    qr[80] = sc[1]*sc[5]*sc[11];

    /*reaction 82: HCO + O2 <=> HO2 + CO */
    qf[81] = sc[3]*sc[13];
    qr[81] = sc[6]*sc[11];

    /*reaction 83: CH3O + O2 <=> HO2 + CH2O */
    qf[82] = sc[3]*sc[15];
    qr[82] = sc[6]*sc[14];

    /*reaction 84: C2H5 + O2 <=> HO2 + C2H4 */
    qf[83] = sc[3]*sc[17];
    qr[83] = sc[6]*sc[16];

    amrex::Real T = tc[1];

    /*compute the mixture concentration */
    amrex::Real mixture = 0.0;
    for (int i = 0; i < 21; ++i) {
        mixture += sc[i];
    }

    amrex::Real Corr[84];
    for (int i = 0; i < 84; ++i) {
        Corr[i] = 1.0;
    }

    /* troe */
    {
        amrex::Real alpha[8];
        alpha[0] = mixture + (TB[0][0] - 1)*sc[0] + (TB[0][1] - 1)*sc[5] + (TB[0][2] - 1)*sc[10] + (TB[0][3] - 1)*sc[11] + (TB[0][4] - 1)*sc[12] + (TB[0][5] - 1)*sc[18] + (TB[0][6] - 1)*sc[20];
        alpha[1] = mixture + (TB[1][0] - 1)*sc[0] + (TB[1][1] - 1)*sc[5] + (TB[1][2] - 1)*sc[10] + (TB[1][3] - 1)*sc[11] + (TB[1][4] - 1)*sc[12] + (TB[1][5] - 1)*sc[18] + (TB[1][6] - 1)*sc[20];
        alpha[2] = mixture + (TB[2][0] - 1)*sc[0] + (TB[2][1] - 1)*sc[5] + (TB[2][2] - 1)*sc[10] + (TB[2][3] - 1)*sc[11] + (TB[2][4] - 1)*sc[12] + (TB[2][5] - 1)*sc[18] + (TB[2][6] - 1)*sc[20];
        alpha[3] = mixture + (TB[3][0] - 1)*sc[0] + (TB[3][1] - 1)*sc[5] + (TB[3][2] - 1)*sc[10] + (TB[3][3] - 1)*sc[11] + (TB[3][4] - 1)*sc[12] + (TB[3][5] - 1)*sc[18];
        alpha[4] = mixture + (TB[4][0] - 1)*sc[0] + (TB[4][1] - 1)*sc[5] + (TB[4][2] - 1)*sc[10] + (TB[4][3] - 1)*sc[11] + (TB[4][4] - 1)*sc[12] + (TB[4][5] - 1)*sc[18] + (TB[4][6] - 1)*sc[20];
        alpha[5] = mixture + (TB[5][0] - 1)*sc[0] + (TB[5][1] - 1)*sc[5] + (TB[5][2] - 1)*sc[10] + (TB[5][3] - 1)*sc[11] + (TB[5][4] - 1)*sc[12] + (TB[5][5] - 1)*sc[18] + (TB[5][6] - 1)*sc[20];
        alpha[6] = mixture + (TB[6][0] - 1)*sc[0] + (TB[6][1] - 1)*sc[5] + (TB[6][2] - 1)*sc[10] + (TB[6][3] - 1)*sc[11] + (TB[6][4] - 1)*sc[12] + (TB[6][5] - 1)*sc[18] + (TB[6][6] - 1)*sc[20];
        alpha[7] = mixture + (TB[7][0] - 1)*sc[0] + (TB[7][1] - 1)*sc[5] + (TB[7][2] - 1)*sc[10] + (TB[7][3] - 1)*sc[11] + (TB[7][4] - 1)*sc[12] + (TB[7][5] - 1)*sc[18] + (TB[7][6] - 1)*sc[20];
        for (int i=0; i<8; i++)
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
        alpha = mixture + (TB[8][0] - 1)*sc[0] + (TB[8][1] - 1)*sc[5] + (TB[8][2] - 1)*sc[10] + (TB[8][3] - 1)*sc[11] + (TB[8][4] - 1)*sc[12] + (TB[8][5] - 1)*sc[18] + (TB[8][6] - 1)*sc[20];
        Corr[8] = alpha;
        alpha = mixture + (TB[9][0] - 1)*sc[0] + (TB[9][1] - 1)*sc[3] + (TB[9][2] - 1)*sc[5] + (TB[9][3] - 1)*sc[10] + (TB[9][4] - 1)*sc[11] + (TB[9][5] - 1)*sc[12] + (TB[9][6] - 1)*sc[18] + (TB[9][7] - 1)*sc[20];
        Corr[9] = alpha;
        alpha = mixture + (TB[10][0] - 1)*sc[3] + (TB[10][1] - 1)*sc[5] + (TB[10][2] - 1)*sc[11] + (TB[10][3] - 1)*sc[12] + (TB[10][4] - 1)*sc[18] + (TB[10][5] - 1)*sc[19] + (TB[10][6] - 1)*sc[20];
        Corr[10] = alpha;
        alpha = mixture + (TB[11][0] - 1)*sc[0] + (TB[11][1] - 1)*sc[5] + (TB[11][2] - 1)*sc[10] + (TB[11][3] - 1)*sc[12] + (TB[11][4] - 1)*sc[18] + (TB[11][5] - 1)*sc[20];
        Corr[11] = alpha;
        alpha = mixture + (TB[12][0] - 1)*sc[0] + (TB[12][1] - 1)*sc[5] + (TB[12][2] - 1)*sc[10] + (TB[12][3] - 1)*sc[18] + (TB[12][4] - 1)*sc[20];
        Corr[12] = alpha;
        alpha = mixture + (TB[13][0] - 1)*sc[0] + (TB[13][1] - 1)*sc[5] + (TB[13][2] - 1)*sc[10] + (TB[13][3] - 1)*sc[11] + (TB[13][4] - 1)*sc[12] + (TB[13][5] - 1)*sc[18];
        Corr[13] = alpha;
    }

    for (int i=0; i<84; i++)
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

    amrex::Real q_f[84], q_r[84];
    amrex::Real sc_qss[1];
    comp_qfqr(q_f, q_r, sc, sc_qss, tc, invT);

    for (int i = 0; i < 84; ++i) {
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
    for (id = 0; id < 84; ++id) {
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

    for (id = 0; id < 84; ++id) {
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
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*H */
    YOW += y[2]*imw[2]; /*O */
    YOW += y[3]*imw[3]; /*O2 */
    YOW += y[4]*imw[4]; /*OH */
    YOW += y[5]*imw[5]; /*H2O */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*CH2 */
    YOW += y[8]*imw[8]; /*CH2(S) */
    YOW += y[9]*imw[9]; /*CH3 */
    YOW += y[10]*imw[10]; /*CH4 */
    YOW += y[11]*imw[11]; /*CO */
    YOW += y[12]*imw[12]; /*CO2 */
    YOW += y[13]*imw[13]; /*HCO */
    YOW += y[14]*imw[14]; /*CH2O */
    YOW += y[15]*imw[15]; /*CH3O */
    YOW += y[16]*imw[16]; /*C2H4 */
    YOW += y[17]*imw[17]; /*C2H5 */
    YOW += y[18]*imw[18]; /*C2H6 */
    YOW += y[19]*imw[19]; /*N2 */
    YOW += y[20]*imw[20]; /*AR */
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
    for (id = 0; id < 84; ++id) {
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
    for (id = 0; id < 84; ++id) {
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
    for (id = 0; id < 84; ++id) {
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
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*1.007970; /*H */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*18.015340; /*H2O */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*14.027090; /*CH2 */
    XW += x[8]*14.027090; /*CH2(S) */
    XW += x[9]*15.035060; /*CH3 */
    XW += x[10]*16.043030; /*CH4 */
    XW += x[11]*28.010550; /*CO */
    XW += x[12]*44.009950; /*CO2 */
    XW += x[13]*29.018520; /*HCO */
    XW += x[14]*30.026490; /*CH2O */
    XW += x[15]*31.034460; /*CH3O */
    XW += x[16]*28.054180; /*C2H4 */
    XW += x[17]*29.062150; /*C2H5 */
    XW += x[18]*30.070120; /*C2H6 */
    XW += x[19]*28.013400; /*N2 */
    XW += x[20]*39.948000; /*AR */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 21; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 84; ++id) {
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
    /*reaction 1: H + CH2 (+M) <=> CH3 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[0][0] - 1)*sc[0] + (TB[0][1] - 1)*sc[5] + (TB[0][2] - 1)*sc[10] + (TB[0][3] - 1)*sc[11] + (TB[0][4] - 1)*sc[12] + (TB[0][5] - 1)*sc[18] + (TB[0][6] - 1)*sc[20];
    /* forward */
    phi_f = sc[1]*sc[7];
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
    phi_r = sc[9];
    Kc = refCinv * exp(g_RT[1] + g_RT[7] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[7]) + (h_RT[9]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[7] -= q; /* CH2 */
    wdot[9] += q; /* CH3 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[0][0] - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[7] -= dqdci;                /* dwdot[CH2]/d[H2] */
        J[9] += dqdci;                /* dwdot[CH3]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[7];
        J[23] -= dqdci;               /* dwdot[H]/d[H] */
        J[29] -= dqdci;               /* dwdot[CH2]/d[H] */
        J[31] += dqdci;               /* dwdot[CH3]/d[H] */
        /* d()/d[H2O] */
        dqdci = (TB[0][1] - 1)*dcdc_fac;
        J[111] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[117] -= dqdci;              /* dwdot[CH2]/d[H2O] */
        J[119] += dqdci;              /* dwdot[CH3]/d[H2O] */
        /* d()/d[CH2] */
        dqdci =  + k_f*sc[1];
        J[155] -= dqdci;              /* dwdot[H]/d[CH2] */
        J[161] -= dqdci;              /* dwdot[CH2]/d[CH2] */
        J[163] += dqdci;              /* dwdot[CH3]/d[CH2] */
        /* d()/d[CH3] */
        dqdci =  - k_r;
        J[199] -= dqdci;              /* dwdot[H]/d[CH3] */
        J[205] -= dqdci;              /* dwdot[CH2]/d[CH3] */
        J[207] += dqdci;              /* dwdot[CH3]/d[CH3] */
        /* d()/d[CH4] */
        dqdci = (TB[0][2] - 1)*dcdc_fac;
        J[221] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[227] -= dqdci;              /* dwdot[CH2]/d[CH4] */
        J[229] += dqdci;              /* dwdot[CH3]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[0][3] - 1)*dcdc_fac;
        J[243] -= dqdci;              /* dwdot[H]/d[CO] */
        J[249] -= dqdci;              /* dwdot[CH2]/d[CO] */
        J[251] += dqdci;              /* dwdot[CH3]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[0][4] - 1)*dcdc_fac;
        J[265] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[271] -= dqdci;              /* dwdot[CH2]/d[CO2] */
        J[273] += dqdci;              /* dwdot[CH3]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (TB[0][5] - 1)*dcdc_fac;
        J[397] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[403] -= dqdci;              /* dwdot[CH2]/d[C2H6] */
        J[405] += dqdci;              /* dwdot[CH3]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (TB[0][6] - 1)*dcdc_fac;
        J[441] -= dqdci;              /* dwdot[H]/d[AR] */
        J[447] -= dqdci;              /* dwdot[CH2]/d[AR] */
        J[449] += dqdci;              /* dwdot[CH3]/d[AR] */
    }
    else {
        dqdc[0] = TB[0][0]*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[7];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = TB[0][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac + k_f*sc[1];
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac - k_r;
        dqdc[10] = TB[0][2]*dcdc_fac;
        dqdc[11] = TB[0][3]*dcdc_fac;
        dqdc[12] = TB[0][4]*dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = TB[0][5]*dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = TB[0][6]*dcdc_fac;
        for (int k=0; k<21; k++) {
            J[22*k+1] -= dqdc[k];
            J[22*k+7] -= dqdc[k];
            J[22*k+9] += dqdc[k];
        }
    }
    J[463] -= dqdT; /* dwdot[H]/dT */
    J[469] -= dqdT; /* dwdot[CH2]/dT */
    J[471] += dqdT; /* dwdot[CH3]/dT */

    /*reaction 2: H + CH3 (+M) <=> CH4 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[1][0] - 1)*sc[0] + (TB[1][1] - 1)*sc[5] + (TB[1][2] - 1)*sc[10] + (TB[1][3] - 1)*sc[11] + (TB[1][4] - 1)*sc[12] + (TB[1][5] - 1)*sc[18] + (TB[1][6] - 1)*sc[20];
    /* forward */
    phi_f = sc[1]*sc[9];
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
    phi_r = sc[10];
    Kc = refCinv * exp(g_RT[1] + g_RT[9] - g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[9]) + (h_RT[10]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[9] -= q; /* CH3 */
    wdot[10] += q; /* CH4 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[1][0] - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[9] -= dqdci;                /* dwdot[CH3]/d[H2] */
        J[10] += dqdci;               /* dwdot[CH4]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[9];
        J[23] -= dqdci;               /* dwdot[H]/d[H] */
        J[31] -= dqdci;               /* dwdot[CH3]/d[H] */
        J[32] += dqdci;               /* dwdot[CH4]/d[H] */
        /* d()/d[H2O] */
        dqdci = (TB[1][1] - 1)*dcdc_fac;
        J[111] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[119] -= dqdci;              /* dwdot[CH3]/d[H2O] */
        J[120] += dqdci;              /* dwdot[CH4]/d[H2O] */
        /* d()/d[CH3] */
        dqdci =  + k_f*sc[1];
        J[199] -= dqdci;              /* dwdot[H]/d[CH3] */
        J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
        J[208] += dqdci;              /* dwdot[CH4]/d[CH3] */
        /* d()/d[CH4] */
        dqdci = (TB[1][2] - 1)*dcdc_fac - k_r;
        J[221] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[229] -= dqdci;              /* dwdot[CH3]/d[CH4] */
        J[230] += dqdci;              /* dwdot[CH4]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[1][3] - 1)*dcdc_fac;
        J[243] -= dqdci;              /* dwdot[H]/d[CO] */
        J[251] -= dqdci;              /* dwdot[CH3]/d[CO] */
        J[252] += dqdci;              /* dwdot[CH4]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[1][4] - 1)*dcdc_fac;
        J[265] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[273] -= dqdci;              /* dwdot[CH3]/d[CO2] */
        J[274] += dqdci;              /* dwdot[CH4]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (TB[1][5] - 1)*dcdc_fac;
        J[397] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[405] -= dqdci;              /* dwdot[CH3]/d[C2H6] */
        J[406] += dqdci;              /* dwdot[CH4]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (TB[1][6] - 1)*dcdc_fac;
        J[441] -= dqdci;              /* dwdot[H]/d[AR] */
        J[449] -= dqdci;              /* dwdot[CH3]/d[AR] */
        J[450] += dqdci;              /* dwdot[CH4]/d[AR] */
    }
    else {
        dqdc[0] = TB[1][0]*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[9];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = TB[1][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac + k_f*sc[1];
        dqdc[10] = TB[1][2]*dcdc_fac - k_r;
        dqdc[11] = TB[1][3]*dcdc_fac;
        dqdc[12] = TB[1][4]*dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = TB[1][5]*dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = TB[1][6]*dcdc_fac;
        for (int k=0; k<21; k++) {
            J[22*k+1] -= dqdc[k];
            J[22*k+9] -= dqdc[k];
            J[22*k+10] += dqdc[k];
        }
    }
    J[463] -= dqdT; /* dwdot[H]/dT */
    J[471] -= dqdT; /* dwdot[CH3]/dT */
    J[472] += dqdT; /* dwdot[CH4]/dT */

    /*reaction 3: H + HCO (+M) <=> CH2O (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[2][0] - 1)*sc[0] + (TB[2][1] - 1)*sc[5] + (TB[2][2] - 1)*sc[10] + (TB[2][3] - 1)*sc[11] + (TB[2][4] - 1)*sc[12] + (TB[2][5] - 1)*sc[18] + (TB[2][6] - 1)*sc[20];
    /* forward */
    phi_f = sc[1]*sc[13];
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
    phi_r = sc[14];
    Kc = refCinv * exp(g_RT[1] + g_RT[13] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[13]) + (h_RT[14]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[13] -= q; /* HCO */
    wdot[14] += q; /* CH2O */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[2][0] - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[13] -= dqdci;               /* dwdot[HCO]/d[H2] */
        J[14] += dqdci;               /* dwdot[CH2O]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[13];
        J[23] -= dqdci;               /* dwdot[H]/d[H] */
        J[35] -= dqdci;               /* dwdot[HCO]/d[H] */
        J[36] += dqdci;               /* dwdot[CH2O]/d[H] */
        /* d()/d[H2O] */
        dqdci = (TB[2][1] - 1)*dcdc_fac;
        J[111] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[123] -= dqdci;              /* dwdot[HCO]/d[H2O] */
        J[124] += dqdci;              /* dwdot[CH2O]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (TB[2][2] - 1)*dcdc_fac;
        J[221] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[233] -= dqdci;              /* dwdot[HCO]/d[CH4] */
        J[234] += dqdci;              /* dwdot[CH2O]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[2][3] - 1)*dcdc_fac;
        J[243] -= dqdci;              /* dwdot[H]/d[CO] */
        J[255] -= dqdci;              /* dwdot[HCO]/d[CO] */
        J[256] += dqdci;              /* dwdot[CH2O]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[2][4] - 1)*dcdc_fac;
        J[265] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[277] -= dqdci;              /* dwdot[HCO]/d[CO2] */
        J[278] += dqdci;              /* dwdot[CH2O]/d[CO2] */
        /* d()/d[HCO] */
        dqdci =  + k_f*sc[1];
        J[287] -= dqdci;              /* dwdot[H]/d[HCO] */
        J[299] -= dqdci;              /* dwdot[HCO]/d[HCO] */
        J[300] += dqdci;              /* dwdot[CH2O]/d[HCO] */
        /* d()/d[CH2O] */
        dqdci =  - k_r;
        J[309] -= dqdci;              /* dwdot[H]/d[CH2O] */
        J[321] -= dqdci;              /* dwdot[HCO]/d[CH2O] */
        J[322] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
        /* d()/d[C2H6] */
        dqdci = (TB[2][5] - 1)*dcdc_fac;
        J[397] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[409] -= dqdci;              /* dwdot[HCO]/d[C2H6] */
        J[410] += dqdci;              /* dwdot[CH2O]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (TB[2][6] - 1)*dcdc_fac;
        J[441] -= dqdci;              /* dwdot[H]/d[AR] */
        J[453] -= dqdci;              /* dwdot[HCO]/d[AR] */
        J[454] += dqdci;              /* dwdot[CH2O]/d[AR] */
    }
    else {
        dqdc[0] = TB[2][0]*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[13];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = TB[2][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = TB[2][2]*dcdc_fac;
        dqdc[11] = TB[2][3]*dcdc_fac;
        dqdc[12] = TB[2][4]*dcdc_fac;
        dqdc[13] = dcdc_fac + k_f*sc[1];
        dqdc[14] = dcdc_fac - k_r;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = TB[2][5]*dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = TB[2][6]*dcdc_fac;
        for (int k=0; k<21; k++) {
            J[22*k+1] -= dqdc[k];
            J[22*k+13] -= dqdc[k];
            J[22*k+14] += dqdc[k];
        }
    }
    J[463] -= dqdT; /* dwdot[H]/dT */
    J[475] -= dqdT; /* dwdot[HCO]/dT */
    J[476] += dqdT; /* dwdot[CH2O]/dT */

    /*reaction 4: H + CH2O (+M) <=> CH3O (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[3][0] - 1)*sc[0] + (TB[3][1] - 1)*sc[5] + (TB[3][2] - 1)*sc[10] + (TB[3][3] - 1)*sc[11] + (TB[3][4] - 1)*sc[12] + (TB[3][5] - 1)*sc[18];
    /* forward */
    phi_f = sc[1]*sc[14];
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
    phi_r = sc[15];
    Kc = refCinv * exp(g_RT[1] + g_RT[14] - g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[14]) + (h_RT[15]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[14] -= q; /* CH2O */
    wdot[15] += q; /* CH3O */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[3][0] - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[14] -= dqdci;               /* dwdot[CH2O]/d[H2] */
        J[15] += dqdci;               /* dwdot[CH3O]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[14];
        J[23] -= dqdci;               /* dwdot[H]/d[H] */
        J[36] -= dqdci;               /* dwdot[CH2O]/d[H] */
        J[37] += dqdci;               /* dwdot[CH3O]/d[H] */
        /* d()/d[H2O] */
        dqdci = (TB[3][1] - 1)*dcdc_fac;
        J[111] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[124] -= dqdci;              /* dwdot[CH2O]/d[H2O] */
        J[125] += dqdci;              /* dwdot[CH3O]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (TB[3][2] - 1)*dcdc_fac;
        J[221] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[234] -= dqdci;              /* dwdot[CH2O]/d[CH4] */
        J[235] += dqdci;              /* dwdot[CH3O]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[3][3] - 1)*dcdc_fac;
        J[243] -= dqdci;              /* dwdot[H]/d[CO] */
        J[256] -= dqdci;              /* dwdot[CH2O]/d[CO] */
        J[257] += dqdci;              /* dwdot[CH3O]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[3][4] - 1)*dcdc_fac;
        J[265] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[278] -= dqdci;              /* dwdot[CH2O]/d[CO2] */
        J[279] += dqdci;              /* dwdot[CH3O]/d[CO2] */
        /* d()/d[CH2O] */
        dqdci =  + k_f*sc[1];
        J[309] -= dqdci;              /* dwdot[H]/d[CH2O] */
        J[322] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
        J[323] += dqdci;              /* dwdot[CH3O]/d[CH2O] */
        /* d()/d[CH3O] */
        dqdci =  - k_r;
        J[331] -= dqdci;              /* dwdot[H]/d[CH3O] */
        J[344] -= dqdci;              /* dwdot[CH2O]/d[CH3O] */
        J[345] += dqdci;              /* dwdot[CH3O]/d[CH3O] */
        /* d()/d[C2H6] */
        dqdci = (TB[3][5] - 1)*dcdc_fac;
        J[397] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[410] -= dqdci;              /* dwdot[CH2O]/d[C2H6] */
        J[411] += dqdci;              /* dwdot[CH3O]/d[C2H6] */
    }
    else {
        dqdc[0] = TB[3][0]*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[14];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = TB[3][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = TB[3][2]*dcdc_fac;
        dqdc[11] = TB[3][3]*dcdc_fac;
        dqdc[12] = TB[3][4]*dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac + k_f*sc[1];
        dqdc[15] = dcdc_fac - k_r;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = TB[3][5]*dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        for (int k=0; k<21; k++) {
            J[22*k+1] -= dqdc[k];
            J[22*k+14] -= dqdc[k];
            J[22*k+15] += dqdc[k];
        }
    }
    J[463] -= dqdT; /* dwdot[H]/dT */
    J[476] -= dqdT; /* dwdot[CH2O]/dT */
    J[477] += dqdT; /* dwdot[CH3O]/dT */

    /*reaction 5: H + C2H4 (+M) <=> C2H5 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[4][0] - 1)*sc[0] + (TB[4][1] - 1)*sc[5] + (TB[4][2] - 1)*sc[10] + (TB[4][3] - 1)*sc[11] + (TB[4][4] - 1)*sc[12] + (TB[4][5] - 1)*sc[18] + (TB[4][6] - 1)*sc[20];
    /* forward */
    phi_f = sc[1]*sc[16];
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
    phi_r = sc[17];
    Kc = refCinv * exp(g_RT[1] + g_RT[16] - g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[16]) + (h_RT[17]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[16] -= q; /* C2H4 */
    wdot[17] += q; /* C2H5 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[4][0] - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[16] -= dqdci;               /* dwdot[C2H4]/d[H2] */
        J[17] += dqdci;               /* dwdot[C2H5]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[16];
        J[23] -= dqdci;               /* dwdot[H]/d[H] */
        J[38] -= dqdci;               /* dwdot[C2H4]/d[H] */
        J[39] += dqdci;               /* dwdot[C2H5]/d[H] */
        /* d()/d[H2O] */
        dqdci = (TB[4][1] - 1)*dcdc_fac;
        J[111] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[126] -= dqdci;              /* dwdot[C2H4]/d[H2O] */
        J[127] += dqdci;              /* dwdot[C2H5]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (TB[4][2] - 1)*dcdc_fac;
        J[221] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[236] -= dqdci;              /* dwdot[C2H4]/d[CH4] */
        J[237] += dqdci;              /* dwdot[C2H5]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[4][3] - 1)*dcdc_fac;
        J[243] -= dqdci;              /* dwdot[H]/d[CO] */
        J[258] -= dqdci;              /* dwdot[C2H4]/d[CO] */
        J[259] += dqdci;              /* dwdot[C2H5]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[4][4] - 1)*dcdc_fac;
        J[265] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[280] -= dqdci;              /* dwdot[C2H4]/d[CO2] */
        J[281] += dqdci;              /* dwdot[C2H5]/d[CO2] */
        /* d()/d[C2H4] */
        dqdci =  + k_f*sc[1];
        J[353] -= dqdci;              /* dwdot[H]/d[C2H4] */
        J[368] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
        J[369] += dqdci;              /* dwdot[C2H5]/d[C2H4] */
        /* d()/d[C2H5] */
        dqdci =  - k_r;
        J[375] -= dqdci;              /* dwdot[H]/d[C2H5] */
        J[390] -= dqdci;              /* dwdot[C2H4]/d[C2H5] */
        J[391] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
        /* d()/d[C2H6] */
        dqdci = (TB[4][5] - 1)*dcdc_fac;
        J[397] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[412] -= dqdci;              /* dwdot[C2H4]/d[C2H6] */
        J[413] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (TB[4][6] - 1)*dcdc_fac;
        J[441] -= dqdci;              /* dwdot[H]/d[AR] */
        J[456] -= dqdci;              /* dwdot[C2H4]/d[AR] */
        J[457] += dqdci;              /* dwdot[C2H5]/d[AR] */
    }
    else {
        dqdc[0] = TB[4][0]*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[16];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = TB[4][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = TB[4][2]*dcdc_fac;
        dqdc[11] = TB[4][3]*dcdc_fac;
        dqdc[12] = TB[4][4]*dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac + k_f*sc[1];
        dqdc[17] = dcdc_fac - k_r;
        dqdc[18] = TB[4][5]*dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = TB[4][6]*dcdc_fac;
        for (int k=0; k<21; k++) {
            J[22*k+1] -= dqdc[k];
            J[22*k+16] -= dqdc[k];
            J[22*k+17] += dqdc[k];
        }
    }
    J[463] -= dqdT; /* dwdot[H]/dT */
    J[478] -= dqdT; /* dwdot[C2H4]/dT */
    J[479] += dqdT; /* dwdot[C2H5]/dT */

    /*reaction 6: H + C2H5 (+M) <=> C2H6 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[5][0] - 1)*sc[0] + (TB[5][1] - 1)*sc[5] + (TB[5][2] - 1)*sc[10] + (TB[5][3] - 1)*sc[11] + (TB[5][4] - 1)*sc[12] + (TB[5][5] - 1)*sc[18] + (TB[5][6] - 1)*sc[20];
    /* forward */
    phi_f = sc[1]*sc[17];
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
    phi_r = sc[18];
    Kc = refCinv * exp(g_RT[1] + g_RT[17] - g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[17]) + (h_RT[18]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[17] -= q; /* C2H5 */
    wdot[18] += q; /* C2H6 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[5][0] - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[17] -= dqdci;               /* dwdot[C2H5]/d[H2] */
        J[18] += dqdci;               /* dwdot[C2H6]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[17];
        J[23] -= dqdci;               /* dwdot[H]/d[H] */
        J[39] -= dqdci;               /* dwdot[C2H5]/d[H] */
        J[40] += dqdci;               /* dwdot[C2H6]/d[H] */
        /* d()/d[H2O] */
        dqdci = (TB[5][1] - 1)*dcdc_fac;
        J[111] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[127] -= dqdci;              /* dwdot[C2H5]/d[H2O] */
        J[128] += dqdci;              /* dwdot[C2H6]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (TB[5][2] - 1)*dcdc_fac;
        J[221] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[237] -= dqdci;              /* dwdot[C2H5]/d[CH4] */
        J[238] += dqdci;              /* dwdot[C2H6]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[5][3] - 1)*dcdc_fac;
        J[243] -= dqdci;              /* dwdot[H]/d[CO] */
        J[259] -= dqdci;              /* dwdot[C2H5]/d[CO] */
        J[260] += dqdci;              /* dwdot[C2H6]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[5][4] - 1)*dcdc_fac;
        J[265] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[281] -= dqdci;              /* dwdot[C2H5]/d[CO2] */
        J[282] += dqdci;              /* dwdot[C2H6]/d[CO2] */
        /* d()/d[C2H5] */
        dqdci =  + k_f*sc[1];
        J[375] -= dqdci;              /* dwdot[H]/d[C2H5] */
        J[391] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
        J[392] += dqdci;              /* dwdot[C2H6]/d[C2H5] */
        /* d()/d[C2H6] */
        dqdci = (TB[5][5] - 1)*dcdc_fac - k_r;
        J[397] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[413] -= dqdci;              /* dwdot[C2H5]/d[C2H6] */
        J[414] += dqdci;              /* dwdot[C2H6]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (TB[5][6] - 1)*dcdc_fac;
        J[441] -= dqdci;              /* dwdot[H]/d[AR] */
        J[457] -= dqdci;              /* dwdot[C2H5]/d[AR] */
        J[458] += dqdci;              /* dwdot[C2H6]/d[AR] */
    }
    else {
        dqdc[0] = TB[5][0]*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[17];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = TB[5][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = TB[5][2]*dcdc_fac;
        dqdc[11] = TB[5][3]*dcdc_fac;
        dqdc[12] = TB[5][4]*dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac + k_f*sc[1];
        dqdc[18] = TB[5][5]*dcdc_fac - k_r;
        dqdc[19] = dcdc_fac;
        dqdc[20] = TB[5][6]*dcdc_fac;
        for (int k=0; k<21; k++) {
            J[22*k+1] -= dqdc[k];
            J[22*k+17] -= dqdc[k];
            J[22*k+18] += dqdc[k];
        }
    }
    J[463] -= dqdT; /* dwdot[H]/dT */
    J[479] -= dqdT; /* dwdot[C2H5]/dT */
    J[480] += dqdT; /* dwdot[C2H6]/dT */

    /*reaction 7: H2 + CO (+M) <=> CH2O (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[6][0] - 1)*sc[0] + (TB[6][1] - 1)*sc[5] + (TB[6][2] - 1)*sc[10] + (TB[6][3] - 1)*sc[11] + (TB[6][4] - 1)*sc[12] + (TB[6][5] - 1)*sc[18] + (TB[6][6] - 1)*sc[20];
    /* forward */
    phi_f = sc[0]*sc[11];
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
    Kc = refCinv * exp(g_RT[0] + g_RT[11] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[11]) + (h_RT[14]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[11] -= q; /* CO */
    wdot[14] += q; /* CH2O */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[6][0] - 1)*dcdc_fac + k_f*sc[11];
        J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
        J[11] -= dqdci;               /* dwdot[CO]/d[H2] */
        J[14] += dqdci;               /* dwdot[CH2O]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[6][1] - 1)*dcdc_fac;
        J[110] -= dqdci;              /* dwdot[H2]/d[H2O] */
        J[121] -= dqdci;              /* dwdot[CO]/d[H2O] */
        J[124] += dqdci;              /* dwdot[CH2O]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (TB[6][2] - 1)*dcdc_fac;
        J[220] -= dqdci;              /* dwdot[H2]/d[CH4] */
        J[231] -= dqdci;              /* dwdot[CO]/d[CH4] */
        J[234] += dqdci;              /* dwdot[CH2O]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[6][3] - 1)*dcdc_fac + k_f*sc[0];
        J[242] -= dqdci;              /* dwdot[H2]/d[CO] */
        J[253] -= dqdci;              /* dwdot[CO]/d[CO] */
        J[256] += dqdci;              /* dwdot[CH2O]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[6][4] - 1)*dcdc_fac;
        J[264] -= dqdci;              /* dwdot[H2]/d[CO2] */
        J[275] -= dqdci;              /* dwdot[CO]/d[CO2] */
        J[278] += dqdci;              /* dwdot[CH2O]/d[CO2] */
        /* d()/d[CH2O] */
        dqdci =  - k_r;
        J[308] -= dqdci;              /* dwdot[H2]/d[CH2O] */
        J[319] -= dqdci;              /* dwdot[CO]/d[CH2O] */
        J[322] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
        /* d()/d[C2H6] */
        dqdci = (TB[6][5] - 1)*dcdc_fac;
        J[396] -= dqdci;              /* dwdot[H2]/d[C2H6] */
        J[407] -= dqdci;              /* dwdot[CO]/d[C2H6] */
        J[410] += dqdci;              /* dwdot[CH2O]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (TB[6][6] - 1)*dcdc_fac;
        J[440] -= dqdci;              /* dwdot[H2]/d[AR] */
        J[451] -= dqdci;              /* dwdot[CO]/d[AR] */
        J[454] += dqdci;              /* dwdot[CH2O]/d[AR] */
    }
    else {
        dqdc[0] = TB[6][0]*dcdc_fac + k_f*sc[11];
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = TB[6][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = TB[6][2]*dcdc_fac;
        dqdc[11] = TB[6][3]*dcdc_fac + k_f*sc[0];
        dqdc[12] = TB[6][4]*dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac - k_r;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = TB[6][5]*dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = TB[6][6]*dcdc_fac;
        for (int k=0; k<21; k++) {
            J[22*k+0] -= dqdc[k];
            J[22*k+11] -= dqdc[k];
            J[22*k+14] += dqdc[k];
        }
    }
    J[462] -= dqdT; /* dwdot[H2]/dT */
    J[473] -= dqdT; /* dwdot[CO]/dT */
    J[476] += dqdT; /* dwdot[CH2O]/dT */

    /*reaction 8: 2.000000 CH3 (+M) <=> C2H6 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[7][0] - 1)*sc[0] + (TB[7][1] - 1)*sc[5] + (TB[7][2] - 1)*sc[10] + (TB[7][3] - 1)*sc[11] + (TB[7][4] - 1)*sc[12] + (TB[7][5] - 1)*sc[18] + (TB[7][6] - 1)*sc[20];
    /* forward */
    phi_f = pow(sc[9], 2.000000);
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
    phi_r = sc[18];
    Kc = refCinv * exp(2.000000*g_RT[9] - g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[9]) + (h_RT[18]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[9] -= 2 * q; /* CH3 */
    wdot[18] += q; /* C2H6 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[7][0] - 1)*dcdc_fac;
        J[9] += -2 * dqdci;           /* dwdot[CH3]/d[H2] */
        J[18] += dqdci;               /* dwdot[C2H6]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[7][1] - 1)*dcdc_fac;
        J[119] += -2 * dqdci;         /* dwdot[CH3]/d[H2O] */
        J[128] += dqdci;              /* dwdot[C2H6]/d[H2O] */
        /* d()/d[CH3] */
        dqdci =  + k_f*2.000000*sc[9];
        J[207] += -2 * dqdci;         /* dwdot[CH3]/d[CH3] */
        J[216] += dqdci;              /* dwdot[C2H6]/d[CH3] */
        /* d()/d[CH4] */
        dqdci = (TB[7][2] - 1)*dcdc_fac;
        J[229] += -2 * dqdci;         /* dwdot[CH3]/d[CH4] */
        J[238] += dqdci;              /* dwdot[C2H6]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[7][3] - 1)*dcdc_fac;
        J[251] += -2 * dqdci;         /* dwdot[CH3]/d[CO] */
        J[260] += dqdci;              /* dwdot[C2H6]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[7][4] - 1)*dcdc_fac;
        J[273] += -2 * dqdci;         /* dwdot[CH3]/d[CO2] */
        J[282] += dqdci;              /* dwdot[C2H6]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (TB[7][5] - 1)*dcdc_fac - k_r;
        J[405] += -2 * dqdci;         /* dwdot[CH3]/d[C2H6] */
        J[414] += dqdci;              /* dwdot[C2H6]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (TB[7][6] - 1)*dcdc_fac;
        J[449] += -2 * dqdci;         /* dwdot[CH3]/d[AR] */
        J[458] += dqdci;              /* dwdot[C2H6]/d[AR] */
    }
    else {
        dqdc[0] = TB[7][0]*dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = TB[7][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac + k_f*2.000000*sc[9];
        dqdc[10] = TB[7][2]*dcdc_fac;
        dqdc[11] = TB[7][3]*dcdc_fac;
        dqdc[12] = TB[7][4]*dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = TB[7][5]*dcdc_fac - k_r;
        dqdc[19] = dcdc_fac;
        dqdc[20] = TB[7][6]*dcdc_fac;
        for (int k=0; k<21; k++) {
            J[22*k+9] += -2 * dqdc[k];
            J[22*k+18] += dqdc[k];
        }
    }
    J[471] += -2 * dqdT; /* dwdot[CH3]/dT */
    J[480] += dqdT; /* dwdot[C2H6]/dT */

    /*reaction 9: O + H + M <=> OH + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[8][0] - 1)*sc[0] + (TB[8][1] - 1)*sc[5] + (TB[8][2] - 1)*sc[10] + (TB[8][3] - 1)*sc[11] + (TB[8][4] - 1)*sc[12] + (TB[8][5] - 1)*sc[18] + (TB[8][6] - 1)*sc[20];
    /* forward */
    phi_f = sc[1]*sc[2];
    k_f = prefactor_units[8] * fwd_A[8]
                * exp(fwd_beta[8] * tc[0] - activation_units[8] * fwd_Ea[8] * invT);
    dlnkfdT = fwd_beta[8] * invT + activation_units[8] * fwd_Ea[8] * invT2;
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
        dqdci = (TB[8][0] - 1)*q_nocor;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[2] -= dqdci;                /* dwdot[O]/d[H2] */
        J[4] += dqdci;                /* dwdot[OH]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[2];
        J[23] -= dqdci;               /* dwdot[H]/d[H] */
        J[24] -= dqdci;               /* dwdot[O]/d[H] */
        J[26] += dqdci;               /* dwdot[OH]/d[H] */
        /* d()/d[O] */
        dqdci =  + k_f*sc[1];
        J[45] -= dqdci;               /* dwdot[H]/d[O] */
        J[46] -= dqdci;               /* dwdot[O]/d[O] */
        J[48] += dqdci;               /* dwdot[OH]/d[O] */
        /* d()/d[OH] */
        dqdci =  - k_r;
        J[89] -= dqdci;               /* dwdot[H]/d[OH] */
        J[90] -= dqdci;               /* dwdot[O]/d[OH] */
        J[92] += dqdci;               /* dwdot[OH]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (TB[8][1] - 1)*q_nocor;
        J[111] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[112] -= dqdci;              /* dwdot[O]/d[H2O] */
        J[114] += dqdci;              /* dwdot[OH]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (TB[8][2] - 1)*q_nocor;
        J[221] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[222] -= dqdci;              /* dwdot[O]/d[CH4] */
        J[224] += dqdci;              /* dwdot[OH]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[8][3] - 1)*q_nocor;
        J[243] -= dqdci;              /* dwdot[H]/d[CO] */
        J[244] -= dqdci;              /* dwdot[O]/d[CO] */
        J[246] += dqdci;              /* dwdot[OH]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[8][4] - 1)*q_nocor;
        J[265] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[266] -= dqdci;              /* dwdot[O]/d[CO2] */
        J[268] += dqdci;              /* dwdot[OH]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (TB[8][5] - 1)*q_nocor;
        J[397] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[398] -= dqdci;              /* dwdot[O]/d[C2H6] */
        J[400] += dqdci;              /* dwdot[OH]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (TB[8][6] - 1)*q_nocor;
        J[441] -= dqdci;              /* dwdot[H]/d[AR] */
        J[442] -= dqdci;              /* dwdot[O]/d[AR] */
        J[444] += dqdci;              /* dwdot[OH]/d[AR] */
    }
    else {
        dqdc[0] = TB[8][0]*q_nocor;
        dqdc[1] = q_nocor + k_f*sc[2];
        dqdc[2] = q_nocor + k_f*sc[1];
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor - k_r;
        dqdc[5] = TB[8][1]*q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = TB[8][2]*q_nocor;
        dqdc[11] = TB[8][3]*q_nocor;
        dqdc[12] = TB[8][4]*q_nocor;
        dqdc[13] = q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = TB[8][5]*q_nocor;
        dqdc[19] = q_nocor;
        dqdc[20] = TB[8][6]*q_nocor;
        for (int k=0; k<21; k++) {
            J[22*k+1] -= dqdc[k];
            J[22*k+2] -= dqdc[k];
            J[22*k+4] += dqdc[k];
        }
    }
    J[463] -= dqdT; /* dwdot[H]/dT */
    J[464] -= dqdT; /* dwdot[O]/dT */
    J[466] += dqdT; /* dwdot[OH]/dT */

    /*reaction 10: O + CO + M <=> CO2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[9][0] - 1)*sc[0] + (TB[9][1] - 1)*sc[3] + (TB[9][2] - 1)*sc[5] + (TB[9][3] - 1)*sc[10] + (TB[9][4] - 1)*sc[11] + (TB[9][5] - 1)*sc[12] + (TB[9][6] - 1)*sc[18] + (TB[9][7] - 1)*sc[20];
    /* forward */
    phi_f = sc[2]*sc[11];
    k_f = prefactor_units[9] * fwd_A[9]
                * exp(fwd_beta[9] * tc[0] - activation_units[9] * fwd_Ea[9] * invT);
    dlnkfdT = fwd_beta[9] * invT + activation_units[9] * fwd_Ea[9] * invT2;
    /* reverse */
    phi_r = sc[12];
    Kc = refCinv * exp(g_RT[2] + g_RT[11] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[11]) + (h_RT[12]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[11] -= q; /* CO */
    wdot[12] += q; /* CO2 */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[9][0] - 1)*q_nocor;
        J[2] -= dqdci;                /* dwdot[O]/d[H2] */
        J[11] -= dqdci;               /* dwdot[CO]/d[H2] */
        J[12] += dqdci;               /* dwdot[CO2]/d[H2] */
        /* d()/d[O] */
        dqdci =  + k_f*sc[11];
        J[46] -= dqdci;               /* dwdot[O]/d[O] */
        J[55] -= dqdci;               /* dwdot[CO]/d[O] */
        J[56] += dqdci;               /* dwdot[CO2]/d[O] */
        /* d()/d[O2] */
        dqdci = (TB[9][1] - 1)*q_nocor;
        J[68] -= dqdci;               /* dwdot[O]/d[O2] */
        J[77] -= dqdci;               /* dwdot[CO]/d[O2] */
        J[78] += dqdci;               /* dwdot[CO2]/d[O2] */
        /* d()/d[H2O] */
        dqdci = (TB[9][2] - 1)*q_nocor;
        J[112] -= dqdci;              /* dwdot[O]/d[H2O] */
        J[121] -= dqdci;              /* dwdot[CO]/d[H2O] */
        J[122] += dqdci;              /* dwdot[CO2]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (TB[9][3] - 1)*q_nocor;
        J[222] -= dqdci;              /* dwdot[O]/d[CH4] */
        J[231] -= dqdci;              /* dwdot[CO]/d[CH4] */
        J[232] += dqdci;              /* dwdot[CO2]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[9][4] - 1)*q_nocor + k_f*sc[2];
        J[244] -= dqdci;              /* dwdot[O]/d[CO] */
        J[253] -= dqdci;              /* dwdot[CO]/d[CO] */
        J[254] += dqdci;              /* dwdot[CO2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[9][5] - 1)*q_nocor - k_r;
        J[266] -= dqdci;              /* dwdot[O]/d[CO2] */
        J[275] -= dqdci;              /* dwdot[CO]/d[CO2] */
        J[276] += dqdci;              /* dwdot[CO2]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (TB[9][6] - 1)*q_nocor;
        J[398] -= dqdci;              /* dwdot[O]/d[C2H6] */
        J[407] -= dqdci;              /* dwdot[CO]/d[C2H6] */
        J[408] += dqdci;              /* dwdot[CO2]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (TB[9][7] - 1)*q_nocor;
        J[442] -= dqdci;              /* dwdot[O]/d[AR] */
        J[451] -= dqdci;              /* dwdot[CO]/d[AR] */
        J[452] += dqdci;              /* dwdot[CO2]/d[AR] */
    }
    else {
        dqdc[0] = TB[9][0]*q_nocor;
        dqdc[1] = q_nocor;
        dqdc[2] = q_nocor + k_f*sc[11];
        dqdc[3] = TB[9][1]*q_nocor;
        dqdc[4] = q_nocor;
        dqdc[5] = TB[9][2]*q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = TB[9][3]*q_nocor;
        dqdc[11] = TB[9][4]*q_nocor + k_f*sc[2];
        dqdc[12] = TB[9][5]*q_nocor - k_r;
        dqdc[13] = q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = TB[9][6]*q_nocor;
        dqdc[19] = q_nocor;
        dqdc[20] = TB[9][7]*q_nocor;
        for (int k=0; k<21; k++) {
            J[22*k+2] -= dqdc[k];
            J[22*k+11] -= dqdc[k];
            J[22*k+12] += dqdc[k];
        }
    }
    J[464] -= dqdT; /* dwdot[O]/dT */
    J[473] -= dqdT; /* dwdot[CO]/dT */
    J[474] += dqdT; /* dwdot[CO2]/dT */

    /*reaction 11: H + O2 + M <=> HO2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[10][0] - 1)*sc[3] + (TB[10][1] - 1)*sc[5] + (TB[10][2] - 1)*sc[11] + (TB[10][3] - 1)*sc[12] + (TB[10][4] - 1)*sc[18] + (TB[10][5] - 1)*sc[19] + (TB[10][6] - 1)*sc[20];
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = prefactor_units[10] * fwd_A[10]
                * exp(fwd_beta[10] * tc[0] - activation_units[10] * fwd_Ea[10] * invT);
    dlnkfdT = fwd_beta[10] * invT + activation_units[10] * fwd_Ea[10] * invT2;
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
        J[23] -= dqdci;               /* dwdot[H]/d[H] */
        J[25] -= dqdci;               /* dwdot[O2]/d[H] */
        J[28] += dqdci;               /* dwdot[HO2]/d[H] */
        /* d()/d[O2] */
        dqdci = (TB[10][0] - 1)*q_nocor + k_f*sc[1];
        J[67] -= dqdci;               /* dwdot[H]/d[O2] */
        J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
        J[72] += dqdci;               /* dwdot[HO2]/d[O2] */
        /* d()/d[H2O] */
        dqdci = (TB[10][1] - 1)*q_nocor;
        J[111] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[113] -= dqdci;              /* dwdot[O2]/d[H2O] */
        J[116] += dqdci;              /* dwdot[HO2]/d[H2O] */
        /* d()/d[HO2] */
        dqdci =  - k_r;
        J[133] -= dqdci;              /* dwdot[H]/d[HO2] */
        J[135] -= dqdci;              /* dwdot[O2]/d[HO2] */
        J[138] += dqdci;              /* dwdot[HO2]/d[HO2] */
        /* d()/d[CO] */
        dqdci = (TB[10][2] - 1)*q_nocor;
        J[243] -= dqdci;              /* dwdot[H]/d[CO] */
        J[245] -= dqdci;              /* dwdot[O2]/d[CO] */
        J[248] += dqdci;              /* dwdot[HO2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[10][3] - 1)*q_nocor;
        J[265] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[267] -= dqdci;              /* dwdot[O2]/d[CO2] */
        J[270] += dqdci;              /* dwdot[HO2]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (TB[10][4] - 1)*q_nocor;
        J[397] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[399] -= dqdci;              /* dwdot[O2]/d[C2H6] */
        J[402] += dqdci;              /* dwdot[HO2]/d[C2H6] */
        /* d()/d[N2] */
        dqdci = (TB[10][5] - 1)*q_nocor;
        J[419] -= dqdci;              /* dwdot[H]/d[N2] */
        J[421] -= dqdci;              /* dwdot[O2]/d[N2] */
        J[424] += dqdci;              /* dwdot[HO2]/d[N2] */
        /* d()/d[AR] */
        dqdci = (TB[10][6] - 1)*q_nocor;
        J[441] -= dqdci;              /* dwdot[H]/d[AR] */
        J[443] -= dqdci;              /* dwdot[O2]/d[AR] */
        J[446] += dqdci;              /* dwdot[HO2]/d[AR] */
    }
    else {
        dqdc[0] = q_nocor;
        dqdc[1] = q_nocor + k_f*sc[3];
        dqdc[2] = q_nocor;
        dqdc[3] = TB[10][0]*q_nocor + k_f*sc[1];
        dqdc[4] = q_nocor;
        dqdc[5] = TB[10][1]*q_nocor;
        dqdc[6] = q_nocor - k_r;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = TB[10][2]*q_nocor;
        dqdc[12] = TB[10][3]*q_nocor;
        dqdc[13] = q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = TB[10][4]*q_nocor;
        dqdc[19] = TB[10][5]*q_nocor;
        dqdc[20] = TB[10][6]*q_nocor;
        for (int k=0; k<21; k++) {
            J[22*k+1] -= dqdc[k];
            J[22*k+3] -= dqdc[k];
            J[22*k+6] += dqdc[k];
        }
    }
    J[463] -= dqdT; /* dwdot[H]/dT */
    J[465] -= dqdT; /* dwdot[O2]/dT */
    J[468] += dqdT; /* dwdot[HO2]/dT */

    /*reaction 12: 2.000000 H + M <=> H2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[11][0] - 1)*sc[0] + (TB[11][1] - 1)*sc[5] + (TB[11][2] - 1)*sc[10] + (TB[11][3] - 1)*sc[12] + (TB[11][4] - 1)*sc[18] + (TB[11][5] - 1)*sc[20];
    /* forward */
    phi_f = pow(sc[1], 2.000000);
    k_f = prefactor_units[11] * fwd_A[11]
                * exp(fwd_beta[11] * tc[0] - activation_units[11] * fwd_Ea[11] * invT);
    dlnkfdT = fwd_beta[11] * invT + activation_units[11] * fwd_Ea[11] * invT2;
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
        dqdci = (TB[11][0] - 1)*q_nocor - k_r;
        J[0] += dqdci;                /* dwdot[H2]/d[H2] */
        J[1] += -2 * dqdci;           /* dwdot[H]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*2.000000*sc[1];
        J[22] += dqdci;               /* dwdot[H2]/d[H] */
        J[23] += -2 * dqdci;          /* dwdot[H]/d[H] */
        /* d()/d[H2O] */
        dqdci = (TB[11][1] - 1)*q_nocor;
        J[110] += dqdci;              /* dwdot[H2]/d[H2O] */
        J[111] += -2 * dqdci;         /* dwdot[H]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (TB[11][2] - 1)*q_nocor;
        J[220] += dqdci;              /* dwdot[H2]/d[CH4] */
        J[221] += -2 * dqdci;         /* dwdot[H]/d[CH4] */
        /* d()/d[CO2] */
        dqdci = (TB[11][3] - 1)*q_nocor;
        J[264] += dqdci;              /* dwdot[H2]/d[CO2] */
        J[265] += -2 * dqdci;         /* dwdot[H]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (TB[11][4] - 1)*q_nocor;
        J[396] += dqdci;              /* dwdot[H2]/d[C2H6] */
        J[397] += -2 * dqdci;         /* dwdot[H]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (TB[11][5] - 1)*q_nocor;
        J[440] += dqdci;              /* dwdot[H2]/d[AR] */
        J[441] += -2 * dqdci;         /* dwdot[H]/d[AR] */
    }
    else {
        dqdc[0] = TB[11][0]*q_nocor - k_r;
        dqdc[1] = q_nocor + k_f*2.000000*sc[1];
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor;
        dqdc[5] = TB[11][1]*q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = TB[11][2]*q_nocor;
        dqdc[11] = q_nocor;
        dqdc[12] = TB[11][3]*q_nocor;
        dqdc[13] = q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = TB[11][4]*q_nocor;
        dqdc[19] = q_nocor;
        dqdc[20] = TB[11][5]*q_nocor;
        for (int k=0; k<21; k++) {
            J[22*k+0] += dqdc[k];
            J[22*k+1] += -2 * dqdc[k];
        }
    }
    J[462] += dqdT; /* dwdot[H2]/dT */
    J[463] += -2 * dqdT; /* dwdot[H]/dT */

    /*reaction 13: H + OH + M <=> H2O + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[12][0] - 1)*sc[0] + (TB[12][1] - 1)*sc[5] + (TB[12][2] - 1)*sc[10] + (TB[12][3] - 1)*sc[18] + (TB[12][4] - 1)*sc[20];
    /* forward */
    phi_f = sc[1]*sc[4];
    k_f = prefactor_units[12] * fwd_A[12]
                * exp(fwd_beta[12] * tc[0] - activation_units[12] * fwd_Ea[12] * invT);
    dlnkfdT = fwd_beta[12] * invT + activation_units[12] * fwd_Ea[12] * invT2;
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
        dqdci = (TB[12][0] - 1)*q_nocor;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[4] -= dqdci;                /* dwdot[OH]/d[H2] */
        J[5] += dqdci;                /* dwdot[H2O]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[4];
        J[23] -= dqdci;               /* dwdot[H]/d[H] */
        J[26] -= dqdci;               /* dwdot[OH]/d[H] */
        J[27] += dqdci;               /* dwdot[H2O]/d[H] */
        /* d()/d[OH] */
        dqdci =  + k_f*sc[1];
        J[89] -= dqdci;               /* dwdot[H]/d[OH] */
        J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
        J[93] += dqdci;               /* dwdot[H2O]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (TB[12][1] - 1)*q_nocor - k_r;
        J[111] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[114] -= dqdci;              /* dwdot[OH]/d[H2O] */
        J[115] += dqdci;              /* dwdot[H2O]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (TB[12][2] - 1)*q_nocor;
        J[221] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[224] -= dqdci;              /* dwdot[OH]/d[CH4] */
        J[225] += dqdci;              /* dwdot[H2O]/d[CH4] */
        /* d()/d[C2H6] */
        dqdci = (TB[12][3] - 1)*q_nocor;
        J[397] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[400] -= dqdci;              /* dwdot[OH]/d[C2H6] */
        J[401] += dqdci;              /* dwdot[H2O]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (TB[12][4] - 1)*q_nocor;
        J[441] -= dqdci;              /* dwdot[H]/d[AR] */
        J[444] -= dqdci;              /* dwdot[OH]/d[AR] */
        J[445] += dqdci;              /* dwdot[H2O]/d[AR] */
    }
    else {
        dqdc[0] = TB[12][0]*q_nocor;
        dqdc[1] = q_nocor + k_f*sc[4];
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor + k_f*sc[1];
        dqdc[5] = TB[12][1]*q_nocor - k_r;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = TB[12][2]*q_nocor;
        dqdc[11] = q_nocor;
        dqdc[12] = q_nocor;
        dqdc[13] = q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = TB[12][3]*q_nocor;
        dqdc[19] = q_nocor;
        dqdc[20] = TB[12][4]*q_nocor;
        for (int k=0; k<21; k++) {
            J[22*k+1] -= dqdc[k];
            J[22*k+4] -= dqdc[k];
            J[22*k+5] += dqdc[k];
        }
    }
    J[463] -= dqdT; /* dwdot[H]/dT */
    J[466] -= dqdT; /* dwdot[OH]/dT */
    J[467] += dqdT; /* dwdot[H2O]/dT */

    /*reaction 14: HCO + M <=> H + CO + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[13][0] - 1)*sc[0] + (TB[13][1] - 1)*sc[5] + (TB[13][2] - 1)*sc[10] + (TB[13][3] - 1)*sc[11] + (TB[13][4] - 1)*sc[12] + (TB[13][5] - 1)*sc[18];
    /* forward */
    phi_f = sc[13];
    k_f = prefactor_units[13] * fwd_A[13]
                * exp(fwd_beta[13] * tc[0] - activation_units[13] * fwd_Ea[13] * invT);
    dlnkfdT = fwd_beta[13] * invT + activation_units[13] * fwd_Ea[13] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[11];
    Kc = refC * exp(-g_RT[1] - g_RT[11] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[13]) + (h_RT[1] + h_RT[11]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[11] += q; /* CO */
    wdot[13] -= q; /* HCO */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[13][0] - 1)*q_nocor;
        J[1] += dqdci;                /* dwdot[H]/d[H2] */
        J[11] += dqdci;               /* dwdot[CO]/d[H2] */
        J[13] -= dqdci;               /* dwdot[HCO]/d[H2] */
        /* d()/d[H] */
        dqdci =  - k_r*sc[11];
        J[23] += dqdci;               /* dwdot[H]/d[H] */
        J[33] += dqdci;               /* dwdot[CO]/d[H] */
        J[35] -= dqdci;               /* dwdot[HCO]/d[H] */
        /* d()/d[H2O] */
        dqdci = (TB[13][1] - 1)*q_nocor;
        J[111] += dqdci;              /* dwdot[H]/d[H2O] */
        J[121] += dqdci;              /* dwdot[CO]/d[H2O] */
        J[123] -= dqdci;              /* dwdot[HCO]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (TB[13][2] - 1)*q_nocor;
        J[221] += dqdci;              /* dwdot[H]/d[CH4] */
        J[231] += dqdci;              /* dwdot[CO]/d[CH4] */
        J[233] -= dqdci;              /* dwdot[HCO]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[13][3] - 1)*q_nocor - k_r*sc[1];
        J[243] += dqdci;              /* dwdot[H]/d[CO] */
        J[253] += dqdci;              /* dwdot[CO]/d[CO] */
        J[255] -= dqdci;              /* dwdot[HCO]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[13][4] - 1)*q_nocor;
        J[265] += dqdci;              /* dwdot[H]/d[CO2] */
        J[275] += dqdci;              /* dwdot[CO]/d[CO2] */
        J[277] -= dqdci;              /* dwdot[HCO]/d[CO2] */
        /* d()/d[HCO] */
        dqdci =  + k_f;
        J[287] += dqdci;              /* dwdot[H]/d[HCO] */
        J[297] += dqdci;              /* dwdot[CO]/d[HCO] */
        J[299] -= dqdci;              /* dwdot[HCO]/d[HCO] */
        /* d()/d[C2H6] */
        dqdci = (TB[13][5] - 1)*q_nocor;
        J[397] += dqdci;              /* dwdot[H]/d[C2H6] */
        J[407] += dqdci;              /* dwdot[CO]/d[C2H6] */
        J[409] -= dqdci;              /* dwdot[HCO]/d[C2H6] */
    }
    else {
        dqdc[0] = TB[13][0]*q_nocor;
        dqdc[1] = q_nocor - k_r*sc[11];
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor;
        dqdc[5] = TB[13][1]*q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = TB[13][2]*q_nocor;
        dqdc[11] = TB[13][3]*q_nocor - k_r*sc[1];
        dqdc[12] = TB[13][4]*q_nocor;
        dqdc[13] = q_nocor + k_f;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = TB[13][5]*q_nocor;
        dqdc[19] = q_nocor;
        dqdc[20] = q_nocor;
        for (int k=0; k<21; k++) {
            J[22*k+1] += dqdc[k];
            J[22*k+11] += dqdc[k];
            J[22*k+13] -= dqdc[k];
        }
    }
    J[463] += dqdT; /* dwdot[H]/dT */
    J[473] += dqdT; /* dwdot[CO]/dT */
    J[475] -= dqdT; /* dwdot[HCO]/dT */

    /*reaction 15: O + H2 <=> H + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[2];
    k_f = prefactor_units[14] * fwd_A[14]
                * exp(fwd_beta[14] * tc[0] - activation_units[14] * fwd_Ea[14] * invT);
    dlnkfdT = fwd_beta[14] * invT + activation_units[14] * fwd_Ea[14] * invT2;
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
    J[22] -= dqdci;               /* dwdot[H2]/d[H] */
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[24] -= dqdci;               /* dwdot[O]/d[H] */
    J[26] += dqdci;               /* dwdot[OH]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[0];
    J[44] -= dqdci;               /* dwdot[H2]/d[O] */
    J[45] += dqdci;               /* dwdot[H]/d[O] */
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[48] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[1];
    J[88] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[89] += dqdci;               /* dwdot[H]/d[OH] */
    J[90] -= dqdci;               /* dwdot[O]/d[OH] */
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[462] -= dqdT;               /* dwdot[H2]/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */

    /*reaction 16: O + HO2 <=> OH + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[6];
    k_f = prefactor_units[15] * fwd_A[15]
                * exp(fwd_beta[15] * tc[0] - activation_units[15] * fwd_Ea[15] * invT);
    dlnkfdT = fwd_beta[15] * invT + activation_units[15] * fwd_Ea[15] * invT2;
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
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[47] += dqdci;               /* dwdot[O2]/d[O] */
    J[48] += dqdci;               /* dwdot[OH]/d[O] */
    J[50] -= dqdci;               /* dwdot[HO2]/d[O] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[4];
    J[68] -= dqdci;               /* dwdot[O]/d[O2] */
    J[69] += dqdci;               /* dwdot[O2]/d[O2] */
    J[70] += dqdci;               /* dwdot[OH]/d[O2] */
    J[72] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[3];
    J[90] -= dqdci;               /* dwdot[O]/d[OH] */
    J[91] += dqdci;               /* dwdot[O2]/d[OH] */
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[94] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[2];
    J[134] -= dqdci;              /* dwdot[O]/d[HO2] */
    J[135] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[136] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[138] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[465] += dqdT;               /* dwdot[O2]/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[468] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 17: O + CH2 <=> H + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[7];
    k_f = prefactor_units[16] * fwd_A[16]
                * exp(fwd_beta[16] * tc[0] - activation_units[16] * fwd_Ea[16] * invT);
    dlnkfdT = fwd_beta[16] * invT + activation_units[16] * fwd_Ea[16] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[13];
    Kc = exp(-g_RT[1] + g_RT[2] + g_RT[7] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[7]) + (h_RT[1] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[2] -= q; /* O */
    wdot[7] -= q; /* CH2 */
    wdot[13] += q; /* HCO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[13];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[24] -= dqdci;               /* dwdot[O]/d[H] */
    J[29] -= dqdci;               /* dwdot[CH2]/d[H] */
    J[35] += dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[7];
    J[45] += dqdci;               /* dwdot[H]/d[O] */
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[51] -= dqdci;               /* dwdot[CH2]/d[O] */
    J[57] += dqdci;               /* dwdot[HCO]/d[O] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[2];
    J[155] += dqdci;              /* dwdot[H]/d[CH2] */
    J[156] -= dqdci;              /* dwdot[O]/d[CH2] */
    J[161] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[167] += dqdci;              /* dwdot[HCO]/d[CH2] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[1];
    J[287] += dqdci;              /* dwdot[H]/d[HCO] */
    J[288] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[293] -= dqdci;              /* dwdot[CH2]/d[HCO] */
    J[299] += dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[469] -= dqdT;               /* dwdot[CH2]/dT */
    J[475] += dqdT;               /* dwdot[HCO]/dT */

    /*reaction 18: O + CH2(S) <=> H + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[8];
    k_f = prefactor_units[17] * fwd_A[17]
                * exp(fwd_beta[17] * tc[0] - activation_units[17] * fwd_Ea[17] * invT);
    dlnkfdT = fwd_beta[17] * invT + activation_units[17] * fwd_Ea[17] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[13];
    Kc = exp(-g_RT[1] + g_RT[2] + g_RT[8] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[8]) + (h_RT[1] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[2] -= q; /* O */
    wdot[8] -= q; /* CH2(S) */
    wdot[13] += q; /* HCO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[13];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[24] -= dqdci;               /* dwdot[O]/d[H] */
    J[30] -= dqdci;               /* dwdot[CH2(S)]/d[H] */
    J[35] += dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[8];
    J[45] += dqdci;               /* dwdot[H]/d[O] */
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[52] -= dqdci;               /* dwdot[CH2(S)]/d[O] */
    J[57] += dqdci;               /* dwdot[HCO]/d[O] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[2];
    J[177] += dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[178] -= dqdci;              /* dwdot[O]/d[CH2(S)] */
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[189] += dqdci;              /* dwdot[HCO]/d[CH2(S)] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[1];
    J[287] += dqdci;              /* dwdot[H]/d[HCO] */
    J[288] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[294] -= dqdci;              /* dwdot[CH2(S)]/d[HCO] */
    J[299] += dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[475] += dqdT;               /* dwdot[HCO]/dT */

    /*reaction 19: O + CH3 <=> H + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[9];
    k_f = prefactor_units[18] * fwd_A[18]
                * exp(fwd_beta[18] * tc[0] - activation_units[18] * fwd_Ea[18] * invT);
    dlnkfdT = fwd_beta[18] * invT + activation_units[18] * fwd_Ea[18] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[14];
    Kc = exp(-g_RT[1] + g_RT[2] + g_RT[9] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[9]) + (h_RT[1] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[2] -= q; /* O */
    wdot[9] -= q; /* CH3 */
    wdot[14] += q; /* CH2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[14];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[24] -= dqdci;               /* dwdot[O]/d[H] */
    J[31] -= dqdci;               /* dwdot[CH3]/d[H] */
    J[36] += dqdci;               /* dwdot[CH2O]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[9];
    J[45] += dqdci;               /* dwdot[H]/d[O] */
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[53] -= dqdci;               /* dwdot[CH3]/d[O] */
    J[58] += dqdci;               /* dwdot[CH2O]/d[O] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[2];
    J[199] += dqdci;              /* dwdot[H]/d[CH3] */
    J[200] -= dqdci;              /* dwdot[O]/d[CH3] */
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[212] += dqdci;              /* dwdot[CH2O]/d[CH3] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[1];
    J[309] += dqdci;              /* dwdot[H]/d[CH2O] */
    J[310] -= dqdci;              /* dwdot[O]/d[CH2O] */
    J[317] -= dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[322] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */
    J[476] += dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 20: O + CH4 <=> OH + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[10];
    k_f = prefactor_units[19] * fwd_A[19]
                * exp(fwd_beta[19] * tc[0] - activation_units[19] * fwd_Ea[19] * invT);
    dlnkfdT = fwd_beta[19] * invT + activation_units[19] * fwd_Ea[19] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[9];
    Kc = exp(g_RT[2] - g_RT[4] - g_RT[9] + g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[10]) + (h_RT[4] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    wdot[9] += q; /* CH3 */
    wdot[10] -= q; /* CH4 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[10];
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[48] += dqdci;               /* dwdot[OH]/d[O] */
    J[53] += dqdci;               /* dwdot[CH3]/d[O] */
    J[54] -= dqdci;               /* dwdot[CH4]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[9];
    J[90] -= dqdci;               /* dwdot[O]/d[OH] */
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[97] += dqdci;               /* dwdot[CH3]/d[OH] */
    J[98] -= dqdci;               /* dwdot[CH4]/d[OH] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[4];
    J[200] -= dqdci;              /* dwdot[O]/d[CH3] */
    J[202] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[207] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[208] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[2];
    J[222] -= dqdci;              /* dwdot[O]/d[CH4] */
    J[224] += dqdci;              /* dwdot[OH]/d[CH4] */
    J[229] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[230] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[471] += dqdT;               /* dwdot[CH3]/dT */
    J[472] -= dqdT;               /* dwdot[CH4]/dT */

    /*reaction 21: O + HCO <=> OH + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[13];
    k_f = prefactor_units[20] * fwd_A[20]
                * exp(fwd_beta[20] * tc[0] - activation_units[20] * fwd_Ea[20] * invT);
    dlnkfdT = fwd_beta[20] * invT + activation_units[20] * fwd_Ea[20] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[11];
    Kc = exp(g_RT[2] - g_RT[4] - g_RT[11] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[13]) + (h_RT[4] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    wdot[11] += q; /* CO */
    wdot[13] -= q; /* HCO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[13];
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[48] += dqdci;               /* dwdot[OH]/d[O] */
    J[55] += dqdci;               /* dwdot[CO]/d[O] */
    J[57] -= dqdci;               /* dwdot[HCO]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[11];
    J[90] -= dqdci;               /* dwdot[O]/d[OH] */
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[99] += dqdci;               /* dwdot[CO]/d[OH] */
    J[101] -= dqdci;              /* dwdot[HCO]/d[OH] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[4];
    J[244] -= dqdci;              /* dwdot[O]/d[CO] */
    J[246] += dqdci;              /* dwdot[OH]/d[CO] */
    J[253] += dqdci;              /* dwdot[CO]/d[CO] */
    J[255] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[2];
    J[288] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[290] += dqdci;              /* dwdot[OH]/d[HCO] */
    J[297] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[299] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[473] += dqdT;               /* dwdot[CO]/dT */
    J[475] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 22: O + HCO <=> H + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[13];
    k_f = prefactor_units[21] * fwd_A[21]
                * exp(fwd_beta[21] * tc[0] - activation_units[21] * fwd_Ea[21] * invT);
    dlnkfdT = fwd_beta[21] * invT + activation_units[21] * fwd_Ea[21] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[12];
    Kc = exp(-g_RT[1] + g_RT[2] - g_RT[12] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[13]) + (h_RT[1] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[2] -= q; /* O */
    wdot[12] += q; /* CO2 */
    wdot[13] -= q; /* HCO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[12];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[24] -= dqdci;               /* dwdot[O]/d[H] */
    J[34] += dqdci;               /* dwdot[CO2]/d[H] */
    J[35] -= dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[13];
    J[45] += dqdci;               /* dwdot[H]/d[O] */
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[56] += dqdci;               /* dwdot[CO2]/d[O] */
    J[57] -= dqdci;               /* dwdot[HCO]/d[O] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[1];
    J[265] += dqdci;              /* dwdot[H]/d[CO2] */
    J[266] -= dqdci;              /* dwdot[O]/d[CO2] */
    J[276] += dqdci;              /* dwdot[CO2]/d[CO2] */
    J[277] -= dqdci;              /* dwdot[HCO]/d[CO2] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[2];
    J[287] += dqdci;              /* dwdot[H]/d[HCO] */
    J[288] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[298] += dqdci;              /* dwdot[CO2]/d[HCO] */
    J[299] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[474] += dqdT;               /* dwdot[CO2]/dT */
    J[475] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 23: O + CH2O <=> OH + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[14];
    k_f = prefactor_units[22] * fwd_A[22]
                * exp(fwd_beta[22] * tc[0] - activation_units[22] * fwd_Ea[22] * invT);
    dlnkfdT = fwd_beta[22] * invT + activation_units[22] * fwd_Ea[22] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[13];
    Kc = exp(g_RT[2] - g_RT[4] - g_RT[13] + g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[14]) + (h_RT[4] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    wdot[13] += q; /* HCO */
    wdot[14] -= q; /* CH2O */
    /* d()/d[O] */
    dqdci =  + k_f*sc[14];
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[48] += dqdci;               /* dwdot[OH]/d[O] */
    J[57] += dqdci;               /* dwdot[HCO]/d[O] */
    J[58] -= dqdci;               /* dwdot[CH2O]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[13];
    J[90] -= dqdci;               /* dwdot[O]/d[OH] */
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[101] += dqdci;              /* dwdot[HCO]/d[OH] */
    J[102] -= dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[4];
    J[288] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[290] += dqdci;              /* dwdot[OH]/d[HCO] */
    J[299] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[300] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[2];
    J[310] -= dqdci;              /* dwdot[O]/d[CH2O] */
    J[312] += dqdci;              /* dwdot[OH]/d[CH2O] */
    J[321] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[322] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[475] += dqdT;               /* dwdot[HCO]/dT */
    J[476] -= dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 24: O + C2H4 <=> CH3 + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[16];
    k_f = prefactor_units[23] * fwd_A[23]
                * exp(fwd_beta[23] * tc[0] - activation_units[23] * fwd_Ea[23] * invT);
    dlnkfdT = fwd_beta[23] * invT + activation_units[23] * fwd_Ea[23] * invT2;
    /* reverse */
    phi_r = sc[9]*sc[13];
    Kc = exp(g_RT[2] - g_RT[9] - g_RT[13] + g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[16]) + (h_RT[9] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[9] += q; /* CH3 */
    wdot[13] += q; /* HCO */
    wdot[16] -= q; /* C2H4 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[16];
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[53] += dqdci;               /* dwdot[CH3]/d[O] */
    J[57] += dqdci;               /* dwdot[HCO]/d[O] */
    J[60] -= dqdci;               /* dwdot[C2H4]/d[O] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[13];
    J[200] -= dqdci;              /* dwdot[O]/d[CH3] */
    J[207] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[211] += dqdci;              /* dwdot[HCO]/d[CH3] */
    J[214] -= dqdci;              /* dwdot[C2H4]/d[CH3] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[9];
    J[288] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[295] += dqdci;              /* dwdot[CH3]/d[HCO] */
    J[299] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[302] -= dqdci;              /* dwdot[C2H4]/d[HCO] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[2];
    J[354] -= dqdci;              /* dwdot[O]/d[C2H4] */
    J[361] += dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[365] += dqdci;              /* dwdot[HCO]/d[C2H4] */
    J[368] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    /* d()/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[471] += dqdT;               /* dwdot[CH3]/dT */
    J[475] += dqdT;               /* dwdot[HCO]/dT */
    J[478] -= dqdT;               /* dwdot[C2H4]/dT */

    /*reaction 25: O + C2H5 <=> CH3 + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[17];
    k_f = prefactor_units[24] * fwd_A[24]
                * exp(fwd_beta[24] * tc[0] - activation_units[24] * fwd_Ea[24] * invT);
    dlnkfdT = fwd_beta[24] * invT + activation_units[24] * fwd_Ea[24] * invT2;
    /* reverse */
    phi_r = sc[9]*sc[14];
    Kc = exp(g_RT[2] - g_RT[9] - g_RT[14] + g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[17]) + (h_RT[9] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[9] += q; /* CH3 */
    wdot[14] += q; /* CH2O */
    wdot[17] -= q; /* C2H5 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[17];
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[53] += dqdci;               /* dwdot[CH3]/d[O] */
    J[58] += dqdci;               /* dwdot[CH2O]/d[O] */
    J[61] -= dqdci;               /* dwdot[C2H5]/d[O] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[14];
    J[200] -= dqdci;              /* dwdot[O]/d[CH3] */
    J[207] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[212] += dqdci;              /* dwdot[CH2O]/d[CH3] */
    J[215] -= dqdci;              /* dwdot[C2H5]/d[CH3] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[9];
    J[310] -= dqdci;              /* dwdot[O]/d[CH2O] */
    J[317] += dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[322] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[325] -= dqdci;              /* dwdot[C2H5]/d[CH2O] */
    /* d()/d[C2H5] */
    dqdci =  + k_f*sc[2];
    J[376] -= dqdci;              /* dwdot[O]/d[C2H5] */
    J[383] += dqdci;              /* dwdot[CH3]/d[C2H5] */
    J[388] += dqdci;              /* dwdot[CH2O]/d[C2H5] */
    J[391] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
    /* d()/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[471] += dqdT;               /* dwdot[CH3]/dT */
    J[476] += dqdT;               /* dwdot[CH2O]/dT */
    J[479] -= dqdT;               /* dwdot[C2H5]/dT */

    /*reaction 26: O + C2H6 <=> OH + C2H5 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[18];
    k_f = prefactor_units[25] * fwd_A[25]
                * exp(fwd_beta[25] * tc[0] - activation_units[25] * fwd_Ea[25] * invT);
    dlnkfdT = fwd_beta[25] * invT + activation_units[25] * fwd_Ea[25] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[17];
    Kc = exp(g_RT[2] - g_RT[4] - g_RT[17] + g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[18]) + (h_RT[4] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    wdot[17] += q; /* C2H5 */
    wdot[18] -= q; /* C2H6 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[18];
    J[46] -= dqdci;               /* dwdot[O]/d[O] */
    J[48] += dqdci;               /* dwdot[OH]/d[O] */
    J[61] += dqdci;               /* dwdot[C2H5]/d[O] */
    J[62] -= dqdci;               /* dwdot[C2H6]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[17];
    J[90] -= dqdci;               /* dwdot[O]/d[OH] */
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[105] += dqdci;              /* dwdot[C2H5]/d[OH] */
    J[106] -= dqdci;              /* dwdot[C2H6]/d[OH] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[4];
    J[376] -= dqdci;              /* dwdot[O]/d[C2H5] */
    J[378] += dqdci;              /* dwdot[OH]/d[C2H5] */
    J[391] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[392] -= dqdci;              /* dwdot[C2H6]/d[C2H5] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[2];
    J[398] -= dqdci;              /* dwdot[O]/d[C2H6] */
    J[400] += dqdci;              /* dwdot[OH]/d[C2H6] */
    J[413] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[414] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[464] -= dqdT;               /* dwdot[O]/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[479] += dqdT;               /* dwdot[C2H5]/dT */
    J[480] -= dqdT;               /* dwdot[C2H6]/dT */

    /*reaction 27: O2 + CO <=> O + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[11];
    k_f = prefactor_units[26] * fwd_A[26]
                * exp(fwd_beta[26] * tc[0] - activation_units[26] * fwd_Ea[26] * invT);
    dlnkfdT = fwd_beta[26] * invT + activation_units[26] * fwd_Ea[26] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[12];
    Kc = exp(-g_RT[2] + g_RT[3] + g_RT[11] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[11]) + (h_RT[2] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* O */
    wdot[3] -= q; /* O2 */
    wdot[11] -= q; /* CO */
    wdot[12] += q; /* CO2 */
    /* d()/d[O] */
    dqdci =  - k_r*sc[12];
    J[46] += dqdci;               /* dwdot[O]/d[O] */
    J[47] -= dqdci;               /* dwdot[O2]/d[O] */
    J[55] -= dqdci;               /* dwdot[CO]/d[O] */
    J[56] += dqdci;               /* dwdot[CO2]/d[O] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[11];
    J[68] += dqdci;               /* dwdot[O]/d[O2] */
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[77] -= dqdci;               /* dwdot[CO]/d[O2] */
    J[78] += dqdci;               /* dwdot[CO2]/d[O2] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[3];
    J[244] += dqdci;              /* dwdot[O]/d[CO] */
    J[245] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[253] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[254] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[2];
    J[266] += dqdci;              /* dwdot[O]/d[CO2] */
    J[267] -= dqdci;              /* dwdot[O2]/d[CO2] */
    J[275] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[276] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[464] += dqdT;               /* dwdot[O]/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[473] -= dqdT;               /* dwdot[CO]/dT */
    J[474] += dqdT;               /* dwdot[CO2]/dT */

    /*reaction 28: O2 + CH2O <=> HO2 + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[14];
    k_f = prefactor_units[27] * fwd_A[27]
                * exp(fwd_beta[27] * tc[0] - activation_units[27] * fwd_Ea[27] * invT);
    dlnkfdT = fwd_beta[27] * invT + activation_units[27] * fwd_Ea[27] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[13];
    Kc = exp(g_RT[3] - g_RT[6] - g_RT[13] + g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[14]) + (h_RT[6] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    wdot[13] += q; /* HCO */
    wdot[14] -= q; /* CH2O */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[14];
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[72] += dqdci;               /* dwdot[HO2]/d[O2] */
    J[79] += dqdci;               /* dwdot[HCO]/d[O2] */
    J[80] -= dqdci;               /* dwdot[CH2O]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[13];
    J[135] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[138] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[145] += dqdci;              /* dwdot[HCO]/d[HO2] */
    J[146] -= dqdci;              /* dwdot[CH2O]/d[HO2] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[6];
    J[289] -= dqdci;              /* dwdot[O2]/d[HCO] */
    J[292] += dqdci;              /* dwdot[HO2]/d[HCO] */
    J[299] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[300] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[3];
    J[311] -= dqdci;              /* dwdot[O2]/d[CH2O] */
    J[314] += dqdci;              /* dwdot[HO2]/d[CH2O] */
    J[321] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[322] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[468] += dqdT;               /* dwdot[HO2]/dT */
    J[475] += dqdT;               /* dwdot[HCO]/dT */
    J[476] -= dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 29: H + 2.000000 O2 <=> HO2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*pow(sc[3], 2.000000);
    k_f = prefactor_units[28] * fwd_A[28]
                * exp(fwd_beta[28] * tc[0] - activation_units[28] * fwd_Ea[28] * invT);
    dlnkfdT = fwd_beta[28] * invT + activation_units[28] * fwd_Ea[28] * invT2;
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
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[25] -= dqdci;               /* dwdot[O2]/d[H] */
    J[28] += dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1]*2.000000*sc[3] - k_r*sc[6];
    J[67] -= dqdci;               /* dwdot[H]/d[O2] */
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[72] += dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[3];
    J[133] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[135] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[138] += dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[468] += dqdT;               /* dwdot[HO2]/dT */

    /*reaction 30: H + O2 + H2O <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3]*sc[5];
    k_f = prefactor_units[29] * fwd_A[29]
                * exp(fwd_beta[29] * tc[0] - activation_units[29] * fwd_Ea[29] * invT);
    dlnkfdT = fwd_beta[29] * invT + activation_units[29] * fwd_Ea[29] * invT2;
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
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[25] -= dqdci;               /* dwdot[O2]/d[H] */
    J[28] += dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1]*sc[5];
    J[67] -= dqdci;               /* dwdot[H]/d[O2] */
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[72] += dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[1]*sc[3] - k_r*sc[6];
    J[111] -= dqdci;              /* dwdot[H]/d[H2O] */
    J[113] -= dqdci;              /* dwdot[O2]/d[H2O] */
    J[116] += dqdci;              /* dwdot[HO2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[5];
    J[133] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[135] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[138] += dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[468] += dqdT;               /* dwdot[HO2]/dT */

    /*reaction 31: H + O2 + N2 <=> HO2 + N2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3]*sc[19];
    k_f = prefactor_units[30] * fwd_A[30]
                * exp(fwd_beta[30] * tc[0] - activation_units[30] * fwd_Ea[30] * invT);
    dlnkfdT = fwd_beta[30] * invT + activation_units[30] * fwd_Ea[30] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[19];
    Kc = refCinv * exp(g_RT[1] + g_RT[3] - g_RT[6] + g_RT[19] - g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3] + h_RT[19]) + (h_RT[6] + h_RT[19]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[3]*sc[19];
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[25] -= dqdci;               /* dwdot[O2]/d[H] */
    J[28] += dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1]*sc[19];
    J[67] -= dqdci;               /* dwdot[H]/d[O2] */
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[72] += dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[19];
    J[133] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[135] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[138] += dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/d[N2] */
    dqdci =  + k_f*sc[1]*sc[3] - k_r*sc[6];
    J[419] -= dqdci;              /* dwdot[H]/d[N2] */
    J[421] -= dqdci;              /* dwdot[O2]/d[N2] */
    J[424] += dqdci;              /* dwdot[HO2]/d[N2] */
    /* d()/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[468] += dqdT;               /* dwdot[HO2]/dT */

    /*reaction 32: H + O2 + AR <=> HO2 + AR */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3]*sc[20];
    k_f = prefactor_units[31] * fwd_A[31]
                * exp(fwd_beta[31] * tc[0] - activation_units[31] * fwd_Ea[31] * invT);
    dlnkfdT = fwd_beta[31] * invT + activation_units[31] * fwd_Ea[31] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[20];
    Kc = refCinv * exp(g_RT[1] + g_RT[3] - g_RT[6] + g_RT[20] - g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3] + h_RT[20]) + (h_RT[6] + h_RT[20]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[3]*sc[20];
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[25] -= dqdci;               /* dwdot[O2]/d[H] */
    J[28] += dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1]*sc[20];
    J[67] -= dqdci;               /* dwdot[H]/d[O2] */
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[72] += dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[20];
    J[133] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[135] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[138] += dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/d[AR] */
    dqdci =  + k_f*sc[1]*sc[3] - k_r*sc[6];
    J[441] -= dqdci;              /* dwdot[H]/d[AR] */
    J[443] -= dqdci;              /* dwdot[O2]/d[AR] */
    J[446] += dqdci;              /* dwdot[HO2]/d[AR] */
    /* d()/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[468] += dqdT;               /* dwdot[HO2]/dT */

    /*reaction 33: H + O2 <=> O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = prefactor_units[32] * fwd_A[32]
                * exp(fwd_beta[32] * tc[0] - activation_units[32] * fwd_Ea[32] * invT);
    dlnkfdT = fwd_beta[32] * invT + activation_units[32] * fwd_Ea[32] * invT2;
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
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[24] += dqdci;               /* dwdot[O]/d[H] */
    J[25] -= dqdci;               /* dwdot[O2]/d[H] */
    J[26] += dqdci;               /* dwdot[OH]/d[H] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[4];
    J[45] -= dqdci;               /* dwdot[H]/d[O] */
    J[46] += dqdci;               /* dwdot[O]/d[O] */
    J[47] -= dqdci;               /* dwdot[O2]/d[O] */
    J[48] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1];
    J[67] -= dqdci;               /* dwdot[H]/d[O2] */
    J[68] += dqdci;               /* dwdot[O]/d[O2] */
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[70] += dqdci;               /* dwdot[OH]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[2];
    J[89] -= dqdci;               /* dwdot[H]/d[OH] */
    J[90] += dqdci;               /* dwdot[O]/d[OH] */
    J[91] -= dqdci;               /* dwdot[O2]/d[OH] */
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[464] += dqdT;               /* dwdot[O]/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */

    /*reaction 34: 2.000000 H + H2 <=> 2.000000 H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*pow(sc[1], 2.000000);
    k_f = prefactor_units[33] * fwd_A[33]
                * exp(fwd_beta[33] * tc[0] - activation_units[33] * fwd_Ea[33] * invT);
    dlnkfdT = fwd_beta[33] * invT + activation_units[33] * fwd_Ea[33] * invT2;
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
    J[22] += dqdci;               /* dwdot[H2]/d[H] */
    J[23] += -2 * dqdci;          /* dwdot[H]/d[H] */
    /* d()/dT */
    J[462] += dqdT;               /* dwdot[H2]/dT */
    J[463] += -2 * dqdT;          /* dwdot[H]/dT */

    /*reaction 35: 2.000000 H + H2O <=> H2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[1], 2.000000)*sc[5];
    k_f = prefactor_units[34] * fwd_A[34]
                * exp(fwd_beta[34] * tc[0] - activation_units[34] * fwd_Ea[34] * invT);
    dlnkfdT = fwd_beta[34] * invT + activation_units[34] * fwd_Ea[34] * invT2;
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
    J[22] += dqdci;               /* dwdot[H2]/d[H] */
    J[23] += -2 * dqdci;          /* dwdot[H]/d[H] */
    /* d()/d[H2O] */
    dqdci =  + k_f*pow(sc[1], 2.000000) - k_r*sc[0];
    J[110] += dqdci;              /* dwdot[H2]/d[H2O] */
    J[111] += -2 * dqdci;         /* dwdot[H]/d[H2O] */
    /* d()/dT */
    J[462] += dqdT;               /* dwdot[H2]/dT */
    J[463] += -2 * dqdT;          /* dwdot[H]/dT */

    /*reaction 36: 2.000000 H + CO2 <=> H2 + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[1], 2.000000)*sc[12];
    k_f = prefactor_units[35] * fwd_A[35]
                * exp(fwd_beta[35] * tc[0] - activation_units[35] * fwd_Ea[35] * invT);
    dlnkfdT = fwd_beta[35] * invT + activation_units[35] * fwd_Ea[35] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[12];
    Kc = refCinv * exp(-g_RT[0] + 2.000000*g_RT[1] + g_RT[12] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[1] + h_RT[12]) + (h_RT[0] + h_RT[12]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= 2 * q; /* H */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[12];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += -2 * dqdci;           /* dwdot[H]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*2.000000*sc[1]*sc[12];
    J[22] += dqdci;               /* dwdot[H2]/d[H] */
    J[23] += -2 * dqdci;          /* dwdot[H]/d[H] */
    /* d()/d[CO2] */
    dqdci =  + k_f*pow(sc[1], 2.000000) - k_r*sc[0];
    J[264] += dqdci;              /* dwdot[H2]/d[CO2] */
    J[265] += -2 * dqdci;         /* dwdot[H]/d[CO2] */
    /* d()/dT */
    J[462] += dqdT;               /* dwdot[H2]/dT */
    J[463] += -2 * dqdT;          /* dwdot[H]/dT */

    /*reaction 37: H + HO2 <=> O2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[6];
    k_f = prefactor_units[36] * fwd_A[36]
                * exp(fwd_beta[36] * tc[0] - activation_units[36] * fwd_Ea[36] * invT);
    dlnkfdT = fwd_beta[36] * invT + activation_units[36] * fwd_Ea[36] * invT2;
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
    J[22] += dqdci;               /* dwdot[H2]/d[H] */
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[25] += dqdci;               /* dwdot[O2]/d[H] */
    J[28] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[0];
    J[66] += dqdci;               /* dwdot[H2]/d[O2] */
    J[67] -= dqdci;               /* dwdot[H]/d[O2] */
    J[69] += dqdci;               /* dwdot[O2]/d[O2] */
    J[72] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[1];
    J[132] += dqdci;              /* dwdot[H2]/d[HO2] */
    J[133] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[135] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[138] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[462] += dqdT;               /* dwdot[H2]/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[465] += dqdT;               /* dwdot[O2]/dT */
    J[468] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 38: H + HO2 <=> 2.000000 OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[6];
    k_f = prefactor_units[37] * fwd_A[37]
                * exp(fwd_beta[37] * tc[0] - activation_units[37] * fwd_Ea[37] * invT);
    dlnkfdT = fwd_beta[37] * invT + activation_units[37] * fwd_Ea[37] * invT2;
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
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[26] += 2 * dqdci;           /* dwdot[OH]/d[H] */
    J[28] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*2.000000*sc[4];
    J[89] -= dqdci;               /* dwdot[H]/d[OH] */
    J[92] += 2 * dqdci;           /* dwdot[OH]/d[OH] */
    J[94] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[1];
    J[133] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[136] += 2 * dqdci;          /* dwdot[OH]/d[HO2] */
    J[138] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[466] += 2 * dqdT;           /* dwdot[OH]/dT */
    J[468] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 39: H + CH4 <=> CH3 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[10];
    k_f = prefactor_units[38] * fwd_A[38]
                * exp(fwd_beta[38] * tc[0] - activation_units[38] * fwd_Ea[38] * invT);
    dlnkfdT = fwd_beta[38] * invT + activation_units[38] * fwd_Ea[38] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[9];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[9] + g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[10]) + (h_RT[0] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[9] += q; /* CH3 */
    wdot[10] -= q; /* CH4 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[9];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[9] += dqdci;                /* dwdot[CH3]/d[H2] */
    J[10] -= dqdci;               /* dwdot[CH4]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[10];
    J[22] += dqdci;               /* dwdot[H2]/d[H] */
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[31] += dqdci;               /* dwdot[CH3]/d[H] */
    J[32] -= dqdci;               /* dwdot[CH4]/d[H] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[0];
    J[198] += dqdci;              /* dwdot[H2]/d[CH3] */
    J[199] -= dqdci;              /* dwdot[H]/d[CH3] */
    J[207] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[208] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[1];
    J[220] += dqdci;              /* dwdot[H2]/d[CH4] */
    J[221] -= dqdci;              /* dwdot[H]/d[CH4] */
    J[229] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[230] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[462] += dqdT;               /* dwdot[H2]/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[471] += dqdT;               /* dwdot[CH3]/dT */
    J[472] -= dqdT;               /* dwdot[CH4]/dT */

    /*reaction 40: H + HCO <=> H2 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[13];
    k_f = prefactor_units[39] * fwd_A[39]
                * exp(fwd_beta[39] * tc[0] - activation_units[39] * fwd_Ea[39] * invT);
    dlnkfdT = fwd_beta[39] * invT + activation_units[39] * fwd_Ea[39] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[11];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[11] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[13]) + (h_RT[0] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[11] += q; /* CO */
    wdot[13] -= q; /* HCO */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[11];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[11] += dqdci;               /* dwdot[CO]/d[H2] */
    J[13] -= dqdci;               /* dwdot[HCO]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[13];
    J[22] += dqdci;               /* dwdot[H2]/d[H] */
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[33] += dqdci;               /* dwdot[CO]/d[H] */
    J[35] -= dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[0];
    J[242] += dqdci;              /* dwdot[H2]/d[CO] */
    J[243] -= dqdci;              /* dwdot[H]/d[CO] */
    J[253] += dqdci;              /* dwdot[CO]/d[CO] */
    J[255] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[1];
    J[286] += dqdci;              /* dwdot[H2]/d[HCO] */
    J[287] -= dqdci;              /* dwdot[H]/d[HCO] */
    J[297] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[299] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[462] += dqdT;               /* dwdot[H2]/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[473] += dqdT;               /* dwdot[CO]/dT */
    J[475] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 41: H + CH2O <=> HCO + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[14];
    k_f = prefactor_units[40] * fwd_A[40]
                * exp(fwd_beta[40] * tc[0] - activation_units[40] * fwd_Ea[40] * invT);
    dlnkfdT = fwd_beta[40] * invT + activation_units[40] * fwd_Ea[40] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[13];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[13] + g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[14]) + (h_RT[0] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[13] += q; /* HCO */
    wdot[14] -= q; /* CH2O */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[13];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[13] += dqdci;               /* dwdot[HCO]/d[H2] */
    J[14] -= dqdci;               /* dwdot[CH2O]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[14];
    J[22] += dqdci;               /* dwdot[H2]/d[H] */
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[35] += dqdci;               /* dwdot[HCO]/d[H] */
    J[36] -= dqdci;               /* dwdot[CH2O]/d[H] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[0];
    J[286] += dqdci;              /* dwdot[H2]/d[HCO] */
    J[287] -= dqdci;              /* dwdot[H]/d[HCO] */
    J[299] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[300] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[1];
    J[308] += dqdci;              /* dwdot[H2]/d[CH2O] */
    J[309] -= dqdci;              /* dwdot[H]/d[CH2O] */
    J[321] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[322] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[462] += dqdT;               /* dwdot[H2]/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[475] += dqdT;               /* dwdot[HCO]/dT */
    J[476] -= dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 42: H + CH3O <=> OH + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[15];
    k_f = prefactor_units[41] * fwd_A[41]
                * exp(fwd_beta[41] * tc[0] - activation_units[41] * fwd_Ea[41] * invT);
    dlnkfdT = fwd_beta[41] * invT + activation_units[41] * fwd_Ea[41] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[9];
    Kc = exp(g_RT[1] - g_RT[4] - g_RT[9] + g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[15]) + (h_RT[4] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[4] += q; /* OH */
    wdot[9] += q; /* CH3 */
    wdot[15] -= q; /* CH3O */
    /* d()/d[H] */
    dqdci =  + k_f*sc[15];
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[26] += dqdci;               /* dwdot[OH]/d[H] */
    J[31] += dqdci;               /* dwdot[CH3]/d[H] */
    J[37] -= dqdci;               /* dwdot[CH3O]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[9];
    J[89] -= dqdci;               /* dwdot[H]/d[OH] */
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[97] += dqdci;               /* dwdot[CH3]/d[OH] */
    J[103] -= dqdci;              /* dwdot[CH3O]/d[OH] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[4];
    J[199] -= dqdci;              /* dwdot[H]/d[CH3] */
    J[202] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[207] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[213] -= dqdci;              /* dwdot[CH3O]/d[CH3] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[1];
    J[331] -= dqdci;              /* dwdot[H]/d[CH3O] */
    J[334] += dqdci;              /* dwdot[OH]/d[CH3O] */
    J[339] += dqdci;              /* dwdot[CH3]/d[CH3O] */
    J[345] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[471] += dqdT;               /* dwdot[CH3]/dT */
    J[477] -= dqdT;               /* dwdot[CH3O]/dT */

    /*reaction 43: H + C2H6 <=> C2H5 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[18];
    k_f = prefactor_units[42] * fwd_A[42]
                * exp(fwd_beta[42] * tc[0] - activation_units[42] * fwd_Ea[42] * invT);
    dlnkfdT = fwd_beta[42] * invT + activation_units[42] * fwd_Ea[42] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[17];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[17] + g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[18]) + (h_RT[0] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[17] += q; /* C2H5 */
    wdot[18] -= q; /* C2H6 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[17];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[17] += dqdci;               /* dwdot[C2H5]/d[H2] */
    J[18] -= dqdci;               /* dwdot[C2H6]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[18];
    J[22] += dqdci;               /* dwdot[H2]/d[H] */
    J[23] -= dqdci;               /* dwdot[H]/d[H] */
    J[39] += dqdci;               /* dwdot[C2H5]/d[H] */
    J[40] -= dqdci;               /* dwdot[C2H6]/d[H] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[0];
    J[374] += dqdci;              /* dwdot[H2]/d[C2H5] */
    J[375] -= dqdci;              /* dwdot[H]/d[C2H5] */
    J[391] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[392] -= dqdci;              /* dwdot[C2H6]/d[C2H5] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[1];
    J[396] += dqdci;              /* dwdot[H2]/d[C2H6] */
    J[397] -= dqdci;              /* dwdot[H]/d[C2H6] */
    J[413] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[414] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[462] += dqdT;               /* dwdot[H2]/dT */
    J[463] -= dqdT;               /* dwdot[H]/dT */
    J[479] += dqdT;               /* dwdot[C2H5]/dT */
    J[480] -= dqdT;               /* dwdot[C2H6]/dT */

    /*reaction 44: OH + H2 <=> H + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[4];
    k_f = prefactor_units[43] * fwd_A[43]
                * exp(fwd_beta[43] * tc[0] - activation_units[43] * fwd_Ea[43] * invT);
    dlnkfdT = fwd_beta[43] * invT + activation_units[43] * fwd_Ea[43] * invT2;
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
    J[22] -= dqdci;               /* dwdot[H2]/d[H] */
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[26] -= dqdci;               /* dwdot[OH]/d[H] */
    J[27] += dqdci;               /* dwdot[H2O]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[0];
    J[88] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[89] += dqdci;               /* dwdot[H]/d[OH] */
    J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[93] += dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[1];
    J[110] -= dqdci;              /* dwdot[H2]/d[H2O] */
    J[111] += dqdci;              /* dwdot[H]/d[H2O] */
    J[114] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[115] += dqdci;              /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[462] -= dqdT;               /* dwdot[H2]/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[466] -= dqdT;               /* dwdot[OH]/dT */
    J[467] += dqdT;               /* dwdot[H2O]/dT */

    /*reaction 45: 2.000000 OH <=> O + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[4], 2.000000);
    k_f = prefactor_units[44] * fwd_A[44]
                * exp(fwd_beta[44] * tc[0] - activation_units[44] * fwd_Ea[44] * invT);
    dlnkfdT = fwd_beta[44] * invT + activation_units[44] * fwd_Ea[44] * invT2;
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
    J[46] += dqdci;               /* dwdot[O]/d[O] */
    J[48] += -2 * dqdci;          /* dwdot[OH]/d[O] */
    J[49] += dqdci;               /* dwdot[H2O]/d[O] */
    /* d()/d[OH] */
    dqdci =  + k_f*2.000000*sc[4];
    J[90] += dqdci;               /* dwdot[O]/d[OH] */
    J[92] += -2 * dqdci;          /* dwdot[OH]/d[OH] */
    J[93] += dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[2];
    J[112] += dqdci;              /* dwdot[O]/d[H2O] */
    J[114] += -2 * dqdci;         /* dwdot[OH]/d[H2O] */
    J[115] += dqdci;              /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[464] += dqdT;               /* dwdot[O]/dT */
    J[466] += -2 * dqdT;          /* dwdot[OH]/dT */
    J[467] += dqdT;               /* dwdot[H2O]/dT */

    /*reaction 46: OH + HO2 <=> O2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[6];
    k_f = prefactor_units[45] * fwd_A[45]
                * exp(fwd_beta[45] * tc[0] - activation_units[45] * fwd_Ea[45] * invT);
    dlnkfdT = fwd_beta[45] * invT + activation_units[45] * fwd_Ea[45] * invT2;
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
    J[69] += dqdci;               /* dwdot[O2]/d[O2] */
    J[70] -= dqdci;               /* dwdot[OH]/d[O2] */
    J[71] += dqdci;               /* dwdot[H2O]/d[O2] */
    J[72] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[6];
    J[91] += dqdci;               /* dwdot[O2]/d[OH] */
    J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[93] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[94] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[3];
    J[113] += dqdci;              /* dwdot[O2]/d[H2O] */
    J[114] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[115] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[116] -= dqdci;              /* dwdot[HO2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[4];
    J[135] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[136] -= dqdci;              /* dwdot[OH]/d[HO2] */
    J[137] += dqdci;              /* dwdot[H2O]/d[HO2] */
    J[138] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[465] += dqdT;               /* dwdot[O2]/dT */
    J[466] -= dqdT;               /* dwdot[OH]/dT */
    J[467] += dqdT;               /* dwdot[H2O]/dT */
    J[468] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 47: OH + CH2 <=> H + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[7];
    k_f = prefactor_units[46] * fwd_A[46]
                * exp(fwd_beta[46] * tc[0] - activation_units[46] * fwd_Ea[46] * invT);
    dlnkfdT = fwd_beta[46] * invT + activation_units[46] * fwd_Ea[46] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[14];
    Kc = exp(-g_RT[1] + g_RT[4] + g_RT[7] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[7]) + (h_RT[1] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[4] -= q; /* OH */
    wdot[7] -= q; /* CH2 */
    wdot[14] += q; /* CH2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[14];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[26] -= dqdci;               /* dwdot[OH]/d[H] */
    J[29] -= dqdci;               /* dwdot[CH2]/d[H] */
    J[36] += dqdci;               /* dwdot[CH2O]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[7];
    J[89] += dqdci;               /* dwdot[H]/d[OH] */
    J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[95] -= dqdci;               /* dwdot[CH2]/d[OH] */
    J[102] += dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[4];
    J[155] += dqdci;              /* dwdot[H]/d[CH2] */
    J[158] -= dqdci;              /* dwdot[OH]/d[CH2] */
    J[161] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[168] += dqdci;              /* dwdot[CH2O]/d[CH2] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[1];
    J[309] += dqdci;              /* dwdot[H]/d[CH2O] */
    J[312] -= dqdci;              /* dwdot[OH]/d[CH2O] */
    J[315] -= dqdci;              /* dwdot[CH2]/d[CH2O] */
    J[322] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[466] -= dqdT;               /* dwdot[OH]/dT */
    J[469] -= dqdT;               /* dwdot[CH2]/dT */
    J[476] += dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 48: OH + CH2(S) <=> H + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[8];
    k_f = prefactor_units[47] * fwd_A[47]
                * exp(fwd_beta[47] * tc[0] - activation_units[47] * fwd_Ea[47] * invT);
    dlnkfdT = fwd_beta[47] * invT + activation_units[47] * fwd_Ea[47] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[14];
    Kc = exp(-g_RT[1] + g_RT[4] + g_RT[8] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[8]) + (h_RT[1] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[4] -= q; /* OH */
    wdot[8] -= q; /* CH2(S) */
    wdot[14] += q; /* CH2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[14];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[26] -= dqdci;               /* dwdot[OH]/d[H] */
    J[30] -= dqdci;               /* dwdot[CH2(S)]/d[H] */
    J[36] += dqdci;               /* dwdot[CH2O]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[8];
    J[89] += dqdci;               /* dwdot[H]/d[OH] */
    J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[96] -= dqdci;               /* dwdot[CH2(S)]/d[OH] */
    J[102] += dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[4];
    J[177] += dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[180] -= dqdci;              /* dwdot[OH]/d[CH2(S)] */
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[190] += dqdci;              /* dwdot[CH2O]/d[CH2(S)] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[1];
    J[309] += dqdci;              /* dwdot[H]/d[CH2O] */
    J[312] -= dqdci;              /* dwdot[OH]/d[CH2O] */
    J[316] -= dqdci;              /* dwdot[CH2(S)]/d[CH2O] */
    J[322] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[466] -= dqdT;               /* dwdot[OH]/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[476] += dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 49: OH + CH3 <=> CH2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[9];
    k_f = prefactor_units[48] * fwd_A[48]
                * exp(fwd_beta[48] * tc[0] - activation_units[48] * fwd_Ea[48] * invT);
    dlnkfdT = fwd_beta[48] * invT + activation_units[48] * fwd_Ea[48] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[7];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[7] + g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[9]) + (h_RT[5] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[7] += q; /* CH2 */
    wdot[9] -= q; /* CH3 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[9];
    J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[93] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[95] += dqdci;               /* dwdot[CH2]/d[OH] */
    J[97] -= dqdci;               /* dwdot[CH3]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[7];
    J[114] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[115] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[117] += dqdci;              /* dwdot[CH2]/d[H2O] */
    J[119] -= dqdci;              /* dwdot[CH3]/d[H2O] */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[5];
    J[158] -= dqdci;              /* dwdot[OH]/d[CH2] */
    J[159] += dqdci;              /* dwdot[H2O]/d[CH2] */
    J[161] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[163] -= dqdci;              /* dwdot[CH3]/d[CH2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[4];
    J[202] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[203] += dqdci;              /* dwdot[H2O]/d[CH3] */
    J[205] += dqdci;              /* dwdot[CH2]/d[CH3] */
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[466] -= dqdT;               /* dwdot[OH]/dT */
    J[467] += dqdT;               /* dwdot[H2O]/dT */
    J[469] += dqdT;               /* dwdot[CH2]/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */

    /*reaction 50: OH + CH3 <=> CH2(S) + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[9];
    k_f = prefactor_units[49] * fwd_A[49]
                * exp(fwd_beta[49] * tc[0] - activation_units[49] * fwd_Ea[49] * invT);
    dlnkfdT = fwd_beta[49] * invT + activation_units[49] * fwd_Ea[49] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[8];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[8] + g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[9]) + (h_RT[5] + h_RT[8]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[8] += q; /* CH2(S) */
    wdot[9] -= q; /* CH3 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[9];
    J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[93] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[96] += dqdci;               /* dwdot[CH2(S)]/d[OH] */
    J[97] -= dqdci;               /* dwdot[CH3]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[8];
    J[114] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[115] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[118] += dqdci;              /* dwdot[CH2(S)]/d[H2O] */
    J[119] -= dqdci;              /* dwdot[CH3]/d[H2O] */
    /* d()/d[CH2(S)] */
    dqdci =  - k_r*sc[5];
    J[180] -= dqdci;              /* dwdot[OH]/d[CH2(S)] */
    J[181] += dqdci;              /* dwdot[H2O]/d[CH2(S)] */
    J[184] += dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[185] -= dqdci;              /* dwdot[CH3]/d[CH2(S)] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[4];
    J[202] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[203] += dqdci;              /* dwdot[H2O]/d[CH3] */
    J[206] += dqdci;              /* dwdot[CH2(S)]/d[CH3] */
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[466] -= dqdT;               /* dwdot[OH]/dT */
    J[467] += dqdT;               /* dwdot[H2O]/dT */
    J[470] += dqdT;               /* dwdot[CH2(S)]/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */

    /*reaction 51: OH + CH4 <=> CH3 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[10];
    k_f = prefactor_units[50] * fwd_A[50]
                * exp(fwd_beta[50] * tc[0] - activation_units[50] * fwd_Ea[50] * invT);
    dlnkfdT = fwd_beta[50] * invT + activation_units[50] * fwd_Ea[50] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[9];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[9] + g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[10]) + (h_RT[5] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[9] += q; /* CH3 */
    wdot[10] -= q; /* CH4 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[10];
    J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[93] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[97] += dqdci;               /* dwdot[CH3]/d[OH] */
    J[98] -= dqdci;               /* dwdot[CH4]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[9];
    J[114] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[115] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[119] += dqdci;              /* dwdot[CH3]/d[H2O] */
    J[120] -= dqdci;              /* dwdot[CH4]/d[H2O] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[5];
    J[202] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[203] += dqdci;              /* dwdot[H2O]/d[CH3] */
    J[207] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[208] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[4];
    J[224] -= dqdci;              /* dwdot[OH]/d[CH4] */
    J[225] += dqdci;              /* dwdot[H2O]/d[CH4] */
    J[229] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[230] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[466] -= dqdT;               /* dwdot[OH]/dT */
    J[467] += dqdT;               /* dwdot[H2O]/dT */
    J[471] += dqdT;               /* dwdot[CH3]/dT */
    J[472] -= dqdT;               /* dwdot[CH4]/dT */

    /*reaction 52: OH + CO <=> H + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[11];
    k_f = prefactor_units[51] * fwd_A[51]
                * exp(fwd_beta[51] * tc[0] - activation_units[51] * fwd_Ea[51] * invT);
    dlnkfdT = fwd_beta[51] * invT + activation_units[51] * fwd_Ea[51] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[12];
    Kc = exp(-g_RT[1] + g_RT[4] + g_RT[11] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[11]) + (h_RT[1] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[4] -= q; /* OH */
    wdot[11] -= q; /* CO */
    wdot[12] += q; /* CO2 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[12];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[26] -= dqdci;               /* dwdot[OH]/d[H] */
    J[33] -= dqdci;               /* dwdot[CO]/d[H] */
    J[34] += dqdci;               /* dwdot[CO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[11];
    J[89] += dqdci;               /* dwdot[H]/d[OH] */
    J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[99] -= dqdci;               /* dwdot[CO]/d[OH] */
    J[100] += dqdci;              /* dwdot[CO2]/d[OH] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[4];
    J[243] += dqdci;              /* dwdot[H]/d[CO] */
    J[246] -= dqdci;              /* dwdot[OH]/d[CO] */
    J[253] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[254] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[1];
    J[265] += dqdci;              /* dwdot[H]/d[CO2] */
    J[268] -= dqdci;              /* dwdot[OH]/d[CO2] */
    J[275] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[276] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[466] -= dqdT;               /* dwdot[OH]/dT */
    J[473] -= dqdT;               /* dwdot[CO]/dT */
    J[474] += dqdT;               /* dwdot[CO2]/dT */

    /*reaction 53: OH + HCO <=> H2O + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[13];
    k_f = prefactor_units[52] * fwd_A[52]
                * exp(fwd_beta[52] * tc[0] - activation_units[52] * fwd_Ea[52] * invT);
    dlnkfdT = fwd_beta[52] * invT + activation_units[52] * fwd_Ea[52] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[11];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[11] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[13]) + (h_RT[5] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[11] += q; /* CO */
    wdot[13] -= q; /* HCO */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[13];
    J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[93] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[99] += dqdci;               /* dwdot[CO]/d[OH] */
    J[101] -= dqdci;              /* dwdot[HCO]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[11];
    J[114] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[115] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[121] += dqdci;              /* dwdot[CO]/d[H2O] */
    J[123] -= dqdci;              /* dwdot[HCO]/d[H2O] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[5];
    J[246] -= dqdci;              /* dwdot[OH]/d[CO] */
    J[247] += dqdci;              /* dwdot[H2O]/d[CO] */
    J[253] += dqdci;              /* dwdot[CO]/d[CO] */
    J[255] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[4];
    J[290] -= dqdci;              /* dwdot[OH]/d[HCO] */
    J[291] += dqdci;              /* dwdot[H2O]/d[HCO] */
    J[297] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[299] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[466] -= dqdT;               /* dwdot[OH]/dT */
    J[467] += dqdT;               /* dwdot[H2O]/dT */
    J[473] += dqdT;               /* dwdot[CO]/dT */
    J[475] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 54: OH + CH2O <=> HCO + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[14];
    k_f = prefactor_units[53] * fwd_A[53]
                * exp(fwd_beta[53] * tc[0] - activation_units[53] * fwd_Ea[53] * invT);
    dlnkfdT = fwd_beta[53] * invT + activation_units[53] * fwd_Ea[53] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[13];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[13] + g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[14]) + (h_RT[5] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[13] += q; /* HCO */
    wdot[14] -= q; /* CH2O */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[14];
    J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[93] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[101] += dqdci;              /* dwdot[HCO]/d[OH] */
    J[102] -= dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[13];
    J[114] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[115] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[123] += dqdci;              /* dwdot[HCO]/d[H2O] */
    J[124] -= dqdci;              /* dwdot[CH2O]/d[H2O] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[5];
    J[290] -= dqdci;              /* dwdot[OH]/d[HCO] */
    J[291] += dqdci;              /* dwdot[H2O]/d[HCO] */
    J[299] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[300] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[4];
    J[312] -= dqdci;              /* dwdot[OH]/d[CH2O] */
    J[313] += dqdci;              /* dwdot[H2O]/d[CH2O] */
    J[321] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[322] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[466] -= dqdT;               /* dwdot[OH]/dT */
    J[467] += dqdT;               /* dwdot[H2O]/dT */
    J[475] += dqdT;               /* dwdot[HCO]/dT */
    J[476] -= dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 55: OH + C2H6 <=> C2H5 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[18];
    k_f = prefactor_units[54] * fwd_A[54]
                * exp(fwd_beta[54] * tc[0] - activation_units[54] * fwd_Ea[54] * invT);
    dlnkfdT = fwd_beta[54] * invT + activation_units[54] * fwd_Ea[54] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[17];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[17] + g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[18]) + (h_RT[5] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[17] += q; /* C2H5 */
    wdot[18] -= q; /* C2H6 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[18];
    J[92] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[93] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[105] += dqdci;              /* dwdot[C2H5]/d[OH] */
    J[106] -= dqdci;              /* dwdot[C2H6]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[17];
    J[114] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[115] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[127] += dqdci;              /* dwdot[C2H5]/d[H2O] */
    J[128] -= dqdci;              /* dwdot[C2H6]/d[H2O] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[5];
    J[378] -= dqdci;              /* dwdot[OH]/d[C2H5] */
    J[379] += dqdci;              /* dwdot[H2O]/d[C2H5] */
    J[391] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[392] -= dqdci;              /* dwdot[C2H6]/d[C2H5] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[4];
    J[400] -= dqdci;              /* dwdot[OH]/d[C2H6] */
    J[401] += dqdci;              /* dwdot[H2O]/d[C2H6] */
    J[413] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[414] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[466] -= dqdT;               /* dwdot[OH]/dT */
    J[467] += dqdT;               /* dwdot[H2O]/dT */
    J[479] += dqdT;               /* dwdot[C2H5]/dT */
    J[480] -= dqdT;               /* dwdot[C2H6]/dT */

    /*reaction 56: HO2 + CH2 <=> OH + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[7];
    k_f = prefactor_units[55] * fwd_A[55]
                * exp(fwd_beta[55] * tc[0] - activation_units[55] * fwd_Ea[55] * invT);
    dlnkfdT = fwd_beta[55] * invT + activation_units[55] * fwd_Ea[55] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[14];
    Kc = exp(-g_RT[4] + g_RT[6] + g_RT[7] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[7]) + (h_RT[4] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] += q; /* OH */
    wdot[6] -= q; /* HO2 */
    wdot[7] -= q; /* CH2 */
    wdot[14] += q; /* CH2O */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[14];
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[94] -= dqdci;               /* dwdot[HO2]/d[OH] */
    J[95] -= dqdci;               /* dwdot[CH2]/d[OH] */
    J[102] += dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[7];
    J[136] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[138] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[139] -= dqdci;              /* dwdot[CH2]/d[HO2] */
    J[146] += dqdci;              /* dwdot[CH2O]/d[HO2] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[6];
    J[158] += dqdci;              /* dwdot[OH]/d[CH2] */
    J[160] -= dqdci;              /* dwdot[HO2]/d[CH2] */
    J[161] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[168] += dqdci;              /* dwdot[CH2O]/d[CH2] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[4];
    J[312] += dqdci;              /* dwdot[OH]/d[CH2O] */
    J[314] -= dqdci;              /* dwdot[HO2]/d[CH2O] */
    J[315] -= dqdci;              /* dwdot[CH2]/d[CH2O] */
    J[322] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[468] -= dqdT;               /* dwdot[HO2]/dT */
    J[469] -= dqdT;               /* dwdot[CH2]/dT */
    J[476] += dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 57: HO2 + CH3 <=> O2 + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[9];
    k_f = prefactor_units[56] * fwd_A[56]
                * exp(fwd_beta[56] * tc[0] - activation_units[56] * fwd_Ea[56] * invT);
    dlnkfdT = fwd_beta[56] * invT + activation_units[56] * fwd_Ea[56] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[10];
    Kc = exp(-g_RT[3] + g_RT[6] + g_RT[9] - g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[9]) + (h_RT[3] + h_RT[10]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += q; /* O2 */
    wdot[6] -= q; /* HO2 */
    wdot[9] -= q; /* CH3 */
    wdot[10] += q; /* CH4 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[10];
    J[69] += dqdci;               /* dwdot[O2]/d[O2] */
    J[72] -= dqdci;               /* dwdot[HO2]/d[O2] */
    J[75] -= dqdci;               /* dwdot[CH3]/d[O2] */
    J[76] += dqdci;               /* dwdot[CH4]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[9];
    J[135] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[138] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[141] -= dqdci;              /* dwdot[CH3]/d[HO2] */
    J[142] += dqdci;              /* dwdot[CH4]/d[HO2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[6];
    J[201] += dqdci;              /* dwdot[O2]/d[CH3] */
    J[204] -= dqdci;              /* dwdot[HO2]/d[CH3] */
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[208] += dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[3];
    J[223] += dqdci;              /* dwdot[O2]/d[CH4] */
    J[226] -= dqdci;              /* dwdot[HO2]/d[CH4] */
    J[229] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[230] += dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[465] += dqdT;               /* dwdot[O2]/dT */
    J[468] -= dqdT;               /* dwdot[HO2]/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */
    J[472] += dqdT;               /* dwdot[CH4]/dT */

    /*reaction 58: HO2 + CH3 <=> OH + CH3O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[9];
    k_f = prefactor_units[57] * fwd_A[57]
                * exp(fwd_beta[57] * tc[0] - activation_units[57] * fwd_Ea[57] * invT);
    dlnkfdT = fwd_beta[57] * invT + activation_units[57] * fwd_Ea[57] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[15];
    Kc = exp(-g_RT[4] + g_RT[6] + g_RT[9] - g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[9]) + (h_RT[4] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] += q; /* OH */
    wdot[6] -= q; /* HO2 */
    wdot[9] -= q; /* CH3 */
    wdot[15] += q; /* CH3O */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[15];
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[94] -= dqdci;               /* dwdot[HO2]/d[OH] */
    J[97] -= dqdci;               /* dwdot[CH3]/d[OH] */
    J[103] += dqdci;              /* dwdot[CH3O]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[9];
    J[136] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[138] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[141] -= dqdci;              /* dwdot[CH3]/d[HO2] */
    J[147] += dqdci;              /* dwdot[CH3O]/d[HO2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[6];
    J[202] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[204] -= dqdci;              /* dwdot[HO2]/d[CH3] */
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[213] += dqdci;              /* dwdot[CH3O]/d[CH3] */
    /* d()/d[CH3O] */
    dqdci =  - k_r*sc[4];
    J[334] += dqdci;              /* dwdot[OH]/d[CH3O] */
    J[336] -= dqdci;              /* dwdot[HO2]/d[CH3O] */
    J[339] -= dqdci;              /* dwdot[CH3]/d[CH3O] */
    J[345] += dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[468] -= dqdT;               /* dwdot[HO2]/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */
    J[477] += dqdT;               /* dwdot[CH3O]/dT */

    /*reaction 59: HO2 + CO <=> OH + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[11];
    k_f = prefactor_units[58] * fwd_A[58]
                * exp(fwd_beta[58] * tc[0] - activation_units[58] * fwd_Ea[58] * invT);
    dlnkfdT = fwd_beta[58] * invT + activation_units[58] * fwd_Ea[58] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[12];
    Kc = exp(-g_RT[4] + g_RT[6] + g_RT[11] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[11]) + (h_RT[4] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] += q; /* OH */
    wdot[6] -= q; /* HO2 */
    wdot[11] -= q; /* CO */
    wdot[12] += q; /* CO2 */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[12];
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[94] -= dqdci;               /* dwdot[HO2]/d[OH] */
    J[99] -= dqdci;               /* dwdot[CO]/d[OH] */
    J[100] += dqdci;              /* dwdot[CO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[11];
    J[136] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[138] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[143] -= dqdci;              /* dwdot[CO]/d[HO2] */
    J[144] += dqdci;              /* dwdot[CO2]/d[HO2] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[6];
    J[246] += dqdci;              /* dwdot[OH]/d[CO] */
    J[248] -= dqdci;              /* dwdot[HO2]/d[CO] */
    J[253] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[254] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[4];
    J[268] += dqdci;              /* dwdot[OH]/d[CO2] */
    J[270] -= dqdci;              /* dwdot[HO2]/d[CO2] */
    J[275] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[276] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[468] -= dqdT;               /* dwdot[HO2]/dT */
    J[473] -= dqdT;               /* dwdot[CO]/dT */
    J[474] += dqdT;               /* dwdot[CO2]/dT */

    /*reaction 60: CH2 + O2 <=> OH + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[7];
    k_f = prefactor_units[59] * fwd_A[59]
                * exp(fwd_beta[59] * tc[0] - activation_units[59] * fwd_Ea[59] * invT);
    dlnkfdT = fwd_beta[59] * invT + activation_units[59] * fwd_Ea[59] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[13];
    Kc = exp(g_RT[3] - g_RT[4] + g_RT[7] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[7]) + (h_RT[4] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[4] += q; /* OH */
    wdot[7] -= q; /* CH2 */
    wdot[13] += q; /* HCO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[7];
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[70] += dqdci;               /* dwdot[OH]/d[O2] */
    J[73] -= dqdci;               /* dwdot[CH2]/d[O2] */
    J[79] += dqdci;               /* dwdot[HCO]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[13];
    J[91] -= dqdci;               /* dwdot[O2]/d[OH] */
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[95] -= dqdci;               /* dwdot[CH2]/d[OH] */
    J[101] += dqdci;              /* dwdot[HCO]/d[OH] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[3];
    J[157] -= dqdci;              /* dwdot[O2]/d[CH2] */
    J[158] += dqdci;              /* dwdot[OH]/d[CH2] */
    J[161] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[167] += dqdci;              /* dwdot[HCO]/d[CH2] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[4];
    J[289] -= dqdci;              /* dwdot[O2]/d[HCO] */
    J[290] += dqdci;              /* dwdot[OH]/d[HCO] */
    J[293] -= dqdci;              /* dwdot[CH2]/d[HCO] */
    J[299] += dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[469] -= dqdT;               /* dwdot[CH2]/dT */
    J[475] += dqdT;               /* dwdot[HCO]/dT */

    /*reaction 61: CH2 + H2 <=> H + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[7];
    k_f = prefactor_units[60] * fwd_A[60]
                * exp(fwd_beta[60] * tc[0] - activation_units[60] * fwd_Ea[60] * invT);
    dlnkfdT = fwd_beta[60] * invT + activation_units[60] * fwd_Ea[60] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[9];
    Kc = exp(g_RT[0] - g_RT[1] + g_RT[7] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[7]) + (h_RT[1] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[1] += q; /* H */
    wdot[7] -= q; /* CH2 */
    wdot[9] += q; /* CH3 */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[7];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += dqdci;                /* dwdot[H]/d[H2] */
    J[7] -= dqdci;                /* dwdot[CH2]/d[H2] */
    J[9] += dqdci;                /* dwdot[CH3]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[9];
    J[22] -= dqdci;               /* dwdot[H2]/d[H] */
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[29] -= dqdci;               /* dwdot[CH2]/d[H] */
    J[31] += dqdci;               /* dwdot[CH3]/d[H] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[0];
    J[154] -= dqdci;              /* dwdot[H2]/d[CH2] */
    J[155] += dqdci;              /* dwdot[H]/d[CH2] */
    J[161] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[163] += dqdci;              /* dwdot[CH3]/d[CH2] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[1];
    J[198] -= dqdci;              /* dwdot[H2]/d[CH3] */
    J[199] += dqdci;              /* dwdot[H]/d[CH3] */
    J[205] -= dqdci;              /* dwdot[CH2]/d[CH3] */
    J[207] += dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[462] -= dqdT;               /* dwdot[H2]/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[469] -= dqdT;               /* dwdot[CH2]/dT */
    J[471] += dqdT;               /* dwdot[CH3]/dT */

    /*reaction 62: CH2 + CH3 <=> H + C2H4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[7]*sc[9];
    k_f = prefactor_units[61] * fwd_A[61]
                * exp(fwd_beta[61] * tc[0] - activation_units[61] * fwd_Ea[61] * invT);
    dlnkfdT = fwd_beta[61] * invT + activation_units[61] * fwd_Ea[61] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[16];
    Kc = exp(-g_RT[1] + g_RT[7] + g_RT[9] - g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7] + h_RT[9]) + (h_RT[1] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[7] -= q; /* CH2 */
    wdot[9] -= q; /* CH3 */
    wdot[16] += q; /* C2H4 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[16];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[29] -= dqdci;               /* dwdot[CH2]/d[H] */
    J[31] -= dqdci;               /* dwdot[CH3]/d[H] */
    J[38] += dqdci;               /* dwdot[C2H4]/d[H] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[9];
    J[155] += dqdci;              /* dwdot[H]/d[CH2] */
    J[161] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[163] -= dqdci;              /* dwdot[CH3]/d[CH2] */
    J[170] += dqdci;              /* dwdot[C2H4]/d[CH2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[7];
    J[199] += dqdci;              /* dwdot[H]/d[CH3] */
    J[205] -= dqdci;              /* dwdot[CH2]/d[CH3] */
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[214] += dqdci;              /* dwdot[C2H4]/d[CH3] */
    /* d()/d[C2H4] */
    dqdci =  - k_r*sc[1];
    J[353] += dqdci;              /* dwdot[H]/d[C2H4] */
    J[359] -= dqdci;              /* dwdot[CH2]/d[C2H4] */
    J[361] -= dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[368] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[469] -= dqdT;               /* dwdot[CH2]/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */
    J[478] += dqdT;               /* dwdot[C2H4]/dT */

    /*reaction 63: CH2 + CH4 <=> 2.000000 CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[7]*sc[10];
    k_f = prefactor_units[62] * fwd_A[62]
                * exp(fwd_beta[62] * tc[0] - activation_units[62] * fwd_Ea[62] * invT);
    dlnkfdT = fwd_beta[62] * invT + activation_units[62] * fwd_Ea[62] * invT2;
    /* reverse */
    phi_r = pow(sc[9], 2.000000);
    Kc = exp(g_RT[7] - 2.000000*g_RT[9] + g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7] + h_RT[10]) + (2.000000*h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] -= q; /* CH2 */
    wdot[9] += 2 * q; /* CH3 */
    wdot[10] -= q; /* CH4 */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[10];
    J[161] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[163] += 2 * dqdci;          /* dwdot[CH3]/d[CH2] */
    J[164] -= dqdci;              /* dwdot[CH4]/d[CH2] */
    /* d()/d[CH3] */
    dqdci =  - k_r*2.000000*sc[9];
    J[205] -= dqdci;              /* dwdot[CH2]/d[CH3] */
    J[207] += 2 * dqdci;          /* dwdot[CH3]/d[CH3] */
    J[208] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[7];
    J[227] -= dqdci;              /* dwdot[CH2]/d[CH4] */
    J[229] += 2 * dqdci;          /* dwdot[CH3]/d[CH4] */
    J[230] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[469] -= dqdT;               /* dwdot[CH2]/dT */
    J[471] += 2 * dqdT;           /* dwdot[CH3]/dT */
    J[472] -= dqdT;               /* dwdot[CH4]/dT */

    /*reaction 64: CH2(S) + N2 <=> CH2 + N2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[19];
    k_f = prefactor_units[63] * fwd_A[63]
                * exp(fwd_beta[63] * tc[0] - activation_units[63] * fwd_Ea[63] * invT);
    dlnkfdT = fwd_beta[63] * invT + activation_units[63] * fwd_Ea[63] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[19];
    Kc = exp(-g_RT[7] + g_RT[8] + g_RT[19] - g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[19]) + (h_RT[7] + h_RT[19]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] += q; /* CH2 */
    wdot[8] -= q; /* CH2(S) */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[19];
    J[161] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[162] -= dqdci;              /* dwdot[CH2(S)]/d[CH2] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[19];
    J[183] += dqdci;              /* dwdot[CH2]/d[CH2(S)] */
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    /* d()/d[N2] */
    dqdci =  + k_f*sc[8] - k_r*sc[7];
    J[425] += dqdci;              /* dwdot[CH2]/d[N2] */
    J[426] -= dqdci;              /* dwdot[CH2(S)]/d[N2] */
    /* d()/dT */
    J[469] += dqdT;               /* dwdot[CH2]/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */

    /*reaction 65: CH2(S) + AR <=> CH2 + AR */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[20];
    k_f = prefactor_units[64] * fwd_A[64]
                * exp(fwd_beta[64] * tc[0] - activation_units[64] * fwd_Ea[64] * invT);
    dlnkfdT = fwd_beta[64] * invT + activation_units[64] * fwd_Ea[64] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[20];
    Kc = exp(-g_RT[7] + g_RT[8] + g_RT[20] - g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[20]) + (h_RT[7] + h_RT[20]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] += q; /* CH2 */
    wdot[8] -= q; /* CH2(S) */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[20];
    J[161] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[162] -= dqdci;              /* dwdot[CH2(S)]/d[CH2] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[20];
    J[183] += dqdci;              /* dwdot[CH2]/d[CH2(S)] */
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    /* d()/d[AR] */
    dqdci =  + k_f*sc[8] - k_r*sc[7];
    J[447] += dqdci;              /* dwdot[CH2]/d[AR] */
    J[448] -= dqdci;              /* dwdot[CH2(S)]/d[AR] */
    /* d()/dT */
    J[469] += dqdT;               /* dwdot[CH2]/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */

    /*reaction 66: CH2(S) + O2 <=> H + OH + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[8];
    k_f = prefactor_units[65] * fwd_A[65]
                * exp(fwd_beta[65] * tc[0] - activation_units[65] * fwd_Ea[65] * invT);
    dlnkfdT = fwd_beta[65] * invT + activation_units[65] * fwd_Ea[65] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[4]*sc[11];
    Kc = refC * exp(-g_RT[1] + g_RT[3] - g_RT[4] + g_RT[8] - g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[8]) + (h_RT[1] + h_RT[4] + h_RT[11]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[3] -= q; /* O2 */
    wdot[4] += q; /* OH */
    wdot[8] -= q; /* CH2(S) */
    wdot[11] += q; /* CO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[4]*sc[11];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[25] -= dqdci;               /* dwdot[O2]/d[H] */
    J[26] += dqdci;               /* dwdot[OH]/d[H] */
    J[30] -= dqdci;               /* dwdot[CH2(S)]/d[H] */
    J[33] += dqdci;               /* dwdot[CO]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[8];
    J[67] += dqdci;               /* dwdot[H]/d[O2] */
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[70] += dqdci;               /* dwdot[OH]/d[O2] */
    J[74] -= dqdci;               /* dwdot[CH2(S)]/d[O2] */
    J[77] += dqdci;               /* dwdot[CO]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[1]*sc[11];
    J[89] += dqdci;               /* dwdot[H]/d[OH] */
    J[91] -= dqdci;               /* dwdot[O2]/d[OH] */
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[96] -= dqdci;               /* dwdot[CH2(S)]/d[OH] */
    J[99] += dqdci;               /* dwdot[CO]/d[OH] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[3];
    J[177] += dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[179] -= dqdci;              /* dwdot[O2]/d[CH2(S)] */
    J[180] += dqdci;              /* dwdot[OH]/d[CH2(S)] */
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[187] += dqdci;              /* dwdot[CO]/d[CH2(S)] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[1]*sc[4];
    J[243] += dqdci;              /* dwdot[H]/d[CO] */
    J[245] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[246] += dqdci;              /* dwdot[OH]/d[CO] */
    J[250] -= dqdci;              /* dwdot[CH2(S)]/d[CO] */
    J[253] += dqdci;              /* dwdot[CO]/d[CO] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[473] += dqdT;               /* dwdot[CO]/dT */

    /*reaction 67: CH2(S) + O2 <=> CO + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[8];
    k_f = prefactor_units[66] * fwd_A[66]
                * exp(fwd_beta[66] * tc[0] - activation_units[66] * fwd_Ea[66] * invT);
    dlnkfdT = fwd_beta[66] * invT + activation_units[66] * fwd_Ea[66] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[11];
    Kc = exp(g_RT[3] - g_RT[5] + g_RT[8] - g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[8]) + (h_RT[5] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[5] += q; /* H2O */
    wdot[8] -= q; /* CH2(S) */
    wdot[11] += q; /* CO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[8];
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[71] += dqdci;               /* dwdot[H2O]/d[O2] */
    J[74] -= dqdci;               /* dwdot[CH2(S)]/d[O2] */
    J[77] += dqdci;               /* dwdot[CO]/d[O2] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[11];
    J[113] -= dqdci;              /* dwdot[O2]/d[H2O] */
    J[115] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[118] -= dqdci;              /* dwdot[CH2(S)]/d[H2O] */
    J[121] += dqdci;              /* dwdot[CO]/d[H2O] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[3];
    J[179] -= dqdci;              /* dwdot[O2]/d[CH2(S)] */
    J[181] += dqdci;              /* dwdot[H2O]/d[CH2(S)] */
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[187] += dqdci;              /* dwdot[CO]/d[CH2(S)] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[5];
    J[245] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[247] += dqdci;              /* dwdot[H2O]/d[CO] */
    J[250] -= dqdci;              /* dwdot[CH2(S)]/d[CO] */
    J[253] += dqdci;              /* dwdot[CO]/d[CO] */
    /* d()/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[467] += dqdT;               /* dwdot[H2O]/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[473] += dqdT;               /* dwdot[CO]/dT */

    /*reaction 68: CH2(S) + H2 <=> CH3 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[8];
    k_f = prefactor_units[67] * fwd_A[67]
                * exp(fwd_beta[67] * tc[0] - activation_units[67] * fwd_Ea[67] * invT);
    dlnkfdT = fwd_beta[67] * invT + activation_units[67] * fwd_Ea[67] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[9];
    Kc = exp(g_RT[0] - g_RT[1] + g_RT[8] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[8]) + (h_RT[1] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[1] += q; /* H */
    wdot[8] -= q; /* CH2(S) */
    wdot[9] += q; /* CH3 */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[8];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += dqdci;                /* dwdot[H]/d[H2] */
    J[8] -= dqdci;                /* dwdot[CH2(S)]/d[H2] */
    J[9] += dqdci;                /* dwdot[CH3]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[9];
    J[22] -= dqdci;               /* dwdot[H2]/d[H] */
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[30] -= dqdci;               /* dwdot[CH2(S)]/d[H] */
    J[31] += dqdci;               /* dwdot[CH3]/d[H] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[0];
    J[176] -= dqdci;              /* dwdot[H2]/d[CH2(S)] */
    J[177] += dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[185] += dqdci;              /* dwdot[CH3]/d[CH2(S)] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[1];
    J[198] -= dqdci;              /* dwdot[H2]/d[CH3] */
    J[199] += dqdci;              /* dwdot[H]/d[CH3] */
    J[206] -= dqdci;              /* dwdot[CH2(S)]/d[CH3] */
    J[207] += dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[462] -= dqdT;               /* dwdot[H2]/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[471] += dqdT;               /* dwdot[CH3]/dT */

    /*reaction 69: CH2(S) + H2O <=> CH2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[8];
    k_f = prefactor_units[68] * fwd_A[68]
                * exp(fwd_beta[68] * tc[0] - activation_units[68] * fwd_Ea[68] * invT);
    dlnkfdT = fwd_beta[68] * invT + activation_units[68] * fwd_Ea[68] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[7];
    Kc = exp(g_RT[5] - g_RT[5] - g_RT[7] + g_RT[8]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[8]) + (h_RT[5] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] += q; /* CH2 */
    wdot[8] -= q; /* CH2(S) */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[8] - k_r*sc[7];
    J[117] += dqdci;              /* dwdot[CH2]/d[H2O] */
    J[118] -= dqdci;              /* dwdot[CH2(S)]/d[H2O] */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[5];
    J[161] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[162] -= dqdci;              /* dwdot[CH2(S)]/d[CH2] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[5];
    J[183] += dqdci;              /* dwdot[CH2]/d[CH2(S)] */
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    /* d()/dT */
    J[469] += dqdT;               /* dwdot[CH2]/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */

    /*reaction 70: CH2(S) + CH3 <=> H + C2H4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[9];
    k_f = prefactor_units[69] * fwd_A[69]
                * exp(fwd_beta[69] * tc[0] - activation_units[69] * fwd_Ea[69] * invT);
    dlnkfdT = fwd_beta[69] * invT + activation_units[69] * fwd_Ea[69] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[16];
    Kc = exp(-g_RT[1] + g_RT[8] + g_RT[9] - g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[9]) + (h_RT[1] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[8] -= q; /* CH2(S) */
    wdot[9] -= q; /* CH3 */
    wdot[16] += q; /* C2H4 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[16];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[30] -= dqdci;               /* dwdot[CH2(S)]/d[H] */
    J[31] -= dqdci;               /* dwdot[CH3]/d[H] */
    J[38] += dqdci;               /* dwdot[C2H4]/d[H] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[9];
    J[177] += dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[185] -= dqdci;              /* dwdot[CH3]/d[CH2(S)] */
    J[192] += dqdci;              /* dwdot[C2H4]/d[CH2(S)] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[8];
    J[199] += dqdci;              /* dwdot[H]/d[CH3] */
    J[206] -= dqdci;              /* dwdot[CH2(S)]/d[CH3] */
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[214] += dqdci;              /* dwdot[C2H4]/d[CH3] */
    /* d()/d[C2H4] */
    dqdci =  - k_r*sc[1];
    J[353] += dqdci;              /* dwdot[H]/d[C2H4] */
    J[360] -= dqdci;              /* dwdot[CH2(S)]/d[C2H4] */
    J[361] -= dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[368] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */
    J[478] += dqdT;               /* dwdot[C2H4]/dT */

    /*reaction 71: CH2(S) + CH4 <=> 2.000000 CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[10];
    k_f = prefactor_units[70] * fwd_A[70]
                * exp(fwd_beta[70] * tc[0] - activation_units[70] * fwd_Ea[70] * invT);
    dlnkfdT = fwd_beta[70] * invT + activation_units[70] * fwd_Ea[70] * invT2;
    /* reverse */
    phi_r = pow(sc[9], 2.000000);
    Kc = exp(g_RT[8] - 2.000000*g_RT[9] + g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[10]) + (2.000000*h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[8] -= q; /* CH2(S) */
    wdot[9] += 2 * q; /* CH3 */
    wdot[10] -= q; /* CH4 */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[10];
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[185] += 2 * dqdci;          /* dwdot[CH3]/d[CH2(S)] */
    J[186] -= dqdci;              /* dwdot[CH4]/d[CH2(S)] */
    /* d()/d[CH3] */
    dqdci =  - k_r*2.000000*sc[9];
    J[206] -= dqdci;              /* dwdot[CH2(S)]/d[CH3] */
    J[207] += 2 * dqdci;          /* dwdot[CH3]/d[CH3] */
    J[208] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[8];
    J[228] -= dqdci;              /* dwdot[CH2(S)]/d[CH4] */
    J[229] += 2 * dqdci;          /* dwdot[CH3]/d[CH4] */
    J[230] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[471] += 2 * dqdT;           /* dwdot[CH3]/dT */
    J[472] -= dqdT;               /* dwdot[CH4]/dT */

    /*reaction 72: CH2(S) + CO <=> CH2 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[11];
    k_f = prefactor_units[71] * fwd_A[71]
                * exp(fwd_beta[71] * tc[0] - activation_units[71] * fwd_Ea[71] * invT);
    dlnkfdT = fwd_beta[71] * invT + activation_units[71] * fwd_Ea[71] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[11];
    Kc = exp(-g_RT[7] + g_RT[8] + g_RT[11] - g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[11]) + (h_RT[7] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] += q; /* CH2 */
    wdot[8] -= q; /* CH2(S) */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[11];
    J[161] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[162] -= dqdci;              /* dwdot[CH2(S)]/d[CH2] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[11];
    J[183] += dqdci;              /* dwdot[CH2]/d[CH2(S)] */
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[8] - k_r*sc[7];
    J[249] += dqdci;              /* dwdot[CH2]/d[CO] */
    J[250] -= dqdci;              /* dwdot[CH2(S)]/d[CO] */
    /* d()/dT */
    J[469] += dqdT;               /* dwdot[CH2]/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */

    /*reaction 73: CH2(S) + CO2 <=> CH2 + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[12];
    k_f = prefactor_units[72] * fwd_A[72]
                * exp(fwd_beta[72] * tc[0] - activation_units[72] * fwd_Ea[72] * invT);
    dlnkfdT = fwd_beta[72] * invT + activation_units[72] * fwd_Ea[72] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[12];
    Kc = exp(-g_RT[7] + g_RT[8] + g_RT[12] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[12]) + (h_RT[7] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] += q; /* CH2 */
    wdot[8] -= q; /* CH2(S) */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[12];
    J[161] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[162] -= dqdci;              /* dwdot[CH2(S)]/d[CH2] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[12];
    J[183] += dqdci;              /* dwdot[CH2]/d[CH2(S)] */
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    /* d()/d[CO2] */
    dqdci =  + k_f*sc[8] - k_r*sc[7];
    J[271] += dqdci;              /* dwdot[CH2]/d[CO2] */
    J[272] -= dqdci;              /* dwdot[CH2(S)]/d[CO2] */
    /* d()/dT */
    J[469] += dqdT;               /* dwdot[CH2]/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */

    /*reaction 74: CH2(S) + CO2 <=> CO + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[12];
    k_f = prefactor_units[73] * fwd_A[73]
                * exp(fwd_beta[73] * tc[0] - activation_units[73] * fwd_Ea[73] * invT);
    dlnkfdT = fwd_beta[73] * invT + activation_units[73] * fwd_Ea[73] * invT2;
    /* reverse */
    phi_r = sc[11]*sc[14];
    Kc = exp(g_RT[8] - g_RT[11] + g_RT[12] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[12]) + (h_RT[11] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[8] -= q; /* CH2(S) */
    wdot[11] += q; /* CO */
    wdot[12] -= q; /* CO2 */
    wdot[14] += q; /* CH2O */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[12];
    J[184] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[187] += dqdci;              /* dwdot[CO]/d[CH2(S)] */
    J[188] -= dqdci;              /* dwdot[CO2]/d[CH2(S)] */
    J[190] += dqdci;              /* dwdot[CH2O]/d[CH2(S)] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[14];
    J[250] -= dqdci;              /* dwdot[CH2(S)]/d[CO] */
    J[253] += dqdci;              /* dwdot[CO]/d[CO] */
    J[254] -= dqdci;              /* dwdot[CO2]/d[CO] */
    J[256] += dqdci;              /* dwdot[CH2O]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  + k_f*sc[8];
    J[272] -= dqdci;              /* dwdot[CH2(S)]/d[CO2] */
    J[275] += dqdci;              /* dwdot[CO]/d[CO2] */
    J[276] -= dqdci;              /* dwdot[CO2]/d[CO2] */
    J[278] += dqdci;              /* dwdot[CH2O]/d[CO2] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[11];
    J[316] -= dqdci;              /* dwdot[CH2(S)]/d[CH2O] */
    J[319] += dqdci;              /* dwdot[CO]/d[CH2O] */
    J[320] -= dqdci;              /* dwdot[CO2]/d[CH2O] */
    J[322] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[470] -= dqdT;               /* dwdot[CH2(S)]/dT */
    J[473] += dqdT;               /* dwdot[CO]/dT */
    J[474] -= dqdT;               /* dwdot[CO2]/dT */
    J[476] += dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 75: CH3 + O2 <=> O + CH3O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[9];
    k_f = prefactor_units[74] * fwd_A[74]
                * exp(fwd_beta[74] * tc[0] - activation_units[74] * fwd_Ea[74] * invT);
    dlnkfdT = fwd_beta[74] * invT + activation_units[74] * fwd_Ea[74] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[15];
    Kc = exp(-g_RT[2] + g_RT[3] + g_RT[9] - g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[9]) + (h_RT[2] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* O */
    wdot[3] -= q; /* O2 */
    wdot[9] -= q; /* CH3 */
    wdot[15] += q; /* CH3O */
    /* d()/d[O] */
    dqdci =  - k_r*sc[15];
    J[46] += dqdci;               /* dwdot[O]/d[O] */
    J[47] -= dqdci;               /* dwdot[O2]/d[O] */
    J[53] -= dqdci;               /* dwdot[CH3]/d[O] */
    J[59] += dqdci;               /* dwdot[CH3O]/d[O] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[9];
    J[68] += dqdci;               /* dwdot[O]/d[O2] */
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[75] -= dqdci;               /* dwdot[CH3]/d[O2] */
    J[81] += dqdci;               /* dwdot[CH3O]/d[O2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[3];
    J[200] += dqdci;              /* dwdot[O]/d[CH3] */
    J[201] -= dqdci;              /* dwdot[O2]/d[CH3] */
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[213] += dqdci;              /* dwdot[CH3O]/d[CH3] */
    /* d()/d[CH3O] */
    dqdci =  - k_r*sc[2];
    J[332] += dqdci;              /* dwdot[O]/d[CH3O] */
    J[333] -= dqdci;              /* dwdot[O2]/d[CH3O] */
    J[339] -= dqdci;              /* dwdot[CH3]/d[CH3O] */
    J[345] += dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[464] += dqdT;               /* dwdot[O]/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */
    J[477] += dqdT;               /* dwdot[CH3O]/dT */

    /*reaction 76: CH3 + O2 <=> OH + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[9];
    k_f = prefactor_units[75] * fwd_A[75]
                * exp(fwd_beta[75] * tc[0] - activation_units[75] * fwd_Ea[75] * invT);
    dlnkfdT = fwd_beta[75] * invT + activation_units[75] * fwd_Ea[75] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[14];
    Kc = exp(g_RT[3] - g_RT[4] + g_RT[9] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[9]) + (h_RT[4] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[4] += q; /* OH */
    wdot[9] -= q; /* CH3 */
    wdot[14] += q; /* CH2O */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[9];
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[70] += dqdci;               /* dwdot[OH]/d[O2] */
    J[75] -= dqdci;               /* dwdot[CH3]/d[O2] */
    J[80] += dqdci;               /* dwdot[CH2O]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[14];
    J[91] -= dqdci;               /* dwdot[O2]/d[OH] */
    J[92] += dqdci;               /* dwdot[OH]/d[OH] */
    J[97] -= dqdci;               /* dwdot[CH3]/d[OH] */
    J[102] += dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[3];
    J[201] -= dqdci;              /* dwdot[O2]/d[CH3] */
    J[202] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[212] += dqdci;              /* dwdot[CH2O]/d[CH3] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[4];
    J[311] -= dqdci;              /* dwdot[O2]/d[CH2O] */
    J[312] += dqdci;              /* dwdot[OH]/d[CH2O] */
    J[317] -= dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[322] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[466] += dqdT;               /* dwdot[OH]/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */
    J[476] += dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 77: 2.000000 CH3 <=> H + C2H5 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[9], 2.000000);
    k_f = prefactor_units[76] * fwd_A[76]
                * exp(fwd_beta[76] * tc[0] - activation_units[76] * fwd_Ea[76] * invT);
    dlnkfdT = fwd_beta[76] * invT + activation_units[76] * fwd_Ea[76] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[17];
    Kc = exp(-g_RT[1] + 2.000000*g_RT[9] - g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[9]) + (h_RT[1] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[9] -= 2 * q; /* CH3 */
    wdot[17] += q; /* C2H5 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[17];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[31] += -2 * dqdci;          /* dwdot[CH3]/d[H] */
    J[39] += dqdci;               /* dwdot[C2H5]/d[H] */
    /* d()/d[CH3] */
    dqdci =  + k_f*2.000000*sc[9];
    J[199] += dqdci;              /* dwdot[H]/d[CH3] */
    J[207] += -2 * dqdci;         /* dwdot[CH3]/d[CH3] */
    J[215] += dqdci;              /* dwdot[C2H5]/d[CH3] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[1];
    J[375] += dqdci;              /* dwdot[H]/d[C2H5] */
    J[383] += -2 * dqdci;         /* dwdot[CH3]/d[C2H5] */
    J[391] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[471] += -2 * dqdT;          /* dwdot[CH3]/dT */
    J[479] += dqdT;               /* dwdot[C2H5]/dT */

    /*reaction 78: CH3 + HCO <=> CH4 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[13];
    k_f = prefactor_units[77] * fwd_A[77]
                * exp(fwd_beta[77] * tc[0] - activation_units[77] * fwd_Ea[77] * invT);
    dlnkfdT = fwd_beta[77] * invT + activation_units[77] * fwd_Ea[77] * invT2;
    /* reverse */
    phi_r = sc[10]*sc[11];
    Kc = exp(g_RT[9] - g_RT[10] - g_RT[11] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[9] + h_RT[13]) + (h_RT[10] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[9] -= q; /* CH3 */
    wdot[10] += q; /* CH4 */
    wdot[11] += q; /* CO */
    wdot[13] -= q; /* HCO */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[13];
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[208] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[209] += dqdci;              /* dwdot[CO]/d[CH3] */
    J[211] -= dqdci;              /* dwdot[HCO]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[11];
    J[229] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[230] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[231] += dqdci;              /* dwdot[CO]/d[CH4] */
    J[233] -= dqdci;              /* dwdot[HCO]/d[CH4] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[10];
    J[251] -= dqdci;              /* dwdot[CH3]/d[CO] */
    J[252] += dqdci;              /* dwdot[CH4]/d[CO] */
    J[253] += dqdci;              /* dwdot[CO]/d[CO] */
    J[255] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[9];
    J[295] -= dqdci;              /* dwdot[CH3]/d[HCO] */
    J[296] += dqdci;              /* dwdot[CH4]/d[HCO] */
    J[297] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[299] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */
    J[472] += dqdT;               /* dwdot[CH4]/dT */
    J[473] += dqdT;               /* dwdot[CO]/dT */
    J[475] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 79: CH3 + CH2O <=> HCO + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[14];
    k_f = prefactor_units[78] * fwd_A[78]
                * exp(fwd_beta[78] * tc[0] - activation_units[78] * fwd_Ea[78] * invT);
    dlnkfdT = fwd_beta[78] * invT + activation_units[78] * fwd_Ea[78] * invT2;
    /* reverse */
    phi_r = sc[10]*sc[13];
    Kc = exp(g_RT[9] - g_RT[10] - g_RT[13] + g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[9] + h_RT[14]) + (h_RT[10] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[9] -= q; /* CH3 */
    wdot[10] += q; /* CH4 */
    wdot[13] += q; /* HCO */
    wdot[14] -= q; /* CH2O */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[14];
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[208] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[211] += dqdci;              /* dwdot[HCO]/d[CH3] */
    J[212] -= dqdci;              /* dwdot[CH2O]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[13];
    J[229] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[230] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[233] += dqdci;              /* dwdot[HCO]/d[CH4] */
    J[234] -= dqdci;              /* dwdot[CH2O]/d[CH4] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[10];
    J[295] -= dqdci;              /* dwdot[CH3]/d[HCO] */
    J[296] += dqdci;              /* dwdot[CH4]/d[HCO] */
    J[299] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[300] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[9];
    J[317] -= dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[318] += dqdci;              /* dwdot[CH4]/d[CH2O] */
    J[321] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[322] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */
    J[472] += dqdT;               /* dwdot[CH4]/dT */
    J[475] += dqdT;               /* dwdot[HCO]/dT */
    J[476] -= dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 80: CH3 + C2H6 <=> C2H5 + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[18];
    k_f = prefactor_units[79] * fwd_A[79]
                * exp(fwd_beta[79] * tc[0] - activation_units[79] * fwd_Ea[79] * invT);
    dlnkfdT = fwd_beta[79] * invT + activation_units[79] * fwd_Ea[79] * invT2;
    /* reverse */
    phi_r = sc[10]*sc[17];
    Kc = exp(g_RT[9] - g_RT[10] - g_RT[17] + g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[9] + h_RT[18]) + (h_RT[10] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[9] -= q; /* CH3 */
    wdot[10] += q; /* CH4 */
    wdot[17] += q; /* C2H5 */
    wdot[18] -= q; /* C2H6 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[18];
    J[207] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[208] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[215] += dqdci;              /* dwdot[C2H5]/d[CH3] */
    J[216] -= dqdci;              /* dwdot[C2H6]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[17];
    J[229] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[230] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[237] += dqdci;              /* dwdot[C2H5]/d[CH4] */
    J[238] -= dqdci;              /* dwdot[C2H6]/d[CH4] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[10];
    J[383] -= dqdci;              /* dwdot[CH3]/d[C2H5] */
    J[384] += dqdci;              /* dwdot[CH4]/d[C2H5] */
    J[391] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[392] -= dqdci;              /* dwdot[C2H6]/d[C2H5] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[9];
    J[405] -= dqdci;              /* dwdot[CH3]/d[C2H6] */
    J[406] += dqdci;              /* dwdot[CH4]/d[C2H6] */
    J[413] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[414] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[471] -= dqdT;               /* dwdot[CH3]/dT */
    J[472] += dqdT;               /* dwdot[CH4]/dT */
    J[479] += dqdT;               /* dwdot[C2H5]/dT */
    J[480] -= dqdT;               /* dwdot[C2H6]/dT */

    /*reaction 81: HCO + H2O <=> H + CO + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[13];
    k_f = prefactor_units[80] * fwd_A[80]
                * exp(fwd_beta[80] * tc[0] - activation_units[80] * fwd_Ea[80] * invT);
    dlnkfdT = fwd_beta[80] * invT + activation_units[80] * fwd_Ea[80] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[5]*sc[11];
    Kc = refC * exp(-g_RT[1] + g_RT[5] - g_RT[5] - g_RT[11] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[13]) + (h_RT[1] + h_RT[5] + h_RT[11]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[11] += q; /* CO */
    wdot[13] -= q; /* HCO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[5]*sc[11];
    J[23] += dqdci;               /* dwdot[H]/d[H] */
    J[33] += dqdci;               /* dwdot[CO]/d[H] */
    J[35] -= dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[13] - k_r*sc[1]*sc[11];
    J[111] += dqdci;              /* dwdot[H]/d[H2O] */
    J[121] += dqdci;              /* dwdot[CO]/d[H2O] */
    J[123] -= dqdci;              /* dwdot[HCO]/d[H2O] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[1]*sc[5];
    J[243] += dqdci;              /* dwdot[H]/d[CO] */
    J[253] += dqdci;              /* dwdot[CO]/d[CO] */
    J[255] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[5];
    J[287] += dqdci;              /* dwdot[H]/d[HCO] */
    J[297] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[299] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[463] += dqdT;               /* dwdot[H]/dT */
    J[473] += dqdT;               /* dwdot[CO]/dT */
    J[475] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 82: HCO + O2 <=> HO2 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[13];
    k_f = prefactor_units[81] * fwd_A[81]
                * exp(fwd_beta[81] * tc[0] - activation_units[81] * fwd_Ea[81] * invT);
    dlnkfdT = fwd_beta[81] * invT + activation_units[81] * fwd_Ea[81] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[11];
    Kc = exp(g_RT[3] - g_RT[6] - g_RT[11] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[13]) + (h_RT[6] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    wdot[11] += q; /* CO */
    wdot[13] -= q; /* HCO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[13];
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[72] += dqdci;               /* dwdot[HO2]/d[O2] */
    J[77] += dqdci;               /* dwdot[CO]/d[O2] */
    J[79] -= dqdci;               /* dwdot[HCO]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[11];
    J[135] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[138] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[143] += dqdci;              /* dwdot[CO]/d[HO2] */
    J[145] -= dqdci;              /* dwdot[HCO]/d[HO2] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[6];
    J[245] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[248] += dqdci;              /* dwdot[HO2]/d[CO] */
    J[253] += dqdci;              /* dwdot[CO]/d[CO] */
    J[255] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[3];
    J[289] -= dqdci;              /* dwdot[O2]/d[HCO] */
    J[292] += dqdci;              /* dwdot[HO2]/d[HCO] */
    J[297] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[299] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[468] += dqdT;               /* dwdot[HO2]/dT */
    J[473] += dqdT;               /* dwdot[CO]/dT */
    J[475] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 83: CH3O + O2 <=> HO2 + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[15];
    k_f = prefactor_units[82] * fwd_A[82]
                * exp(fwd_beta[82] * tc[0] - activation_units[82] * fwd_Ea[82] * invT);
    dlnkfdT = fwd_beta[82] * invT + activation_units[82] * fwd_Ea[82] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[14];
    Kc = exp(g_RT[3] - g_RT[6] - g_RT[14] + g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[15]) + (h_RT[6] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    wdot[14] += q; /* CH2O */
    wdot[15] -= q; /* CH3O */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[15];
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[72] += dqdci;               /* dwdot[HO2]/d[O2] */
    J[80] += dqdci;               /* dwdot[CH2O]/d[O2] */
    J[81] -= dqdci;               /* dwdot[CH3O]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[14];
    J[135] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[138] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[146] += dqdci;              /* dwdot[CH2O]/d[HO2] */
    J[147] -= dqdci;              /* dwdot[CH3O]/d[HO2] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[6];
    J[311] -= dqdci;              /* dwdot[O2]/d[CH2O] */
    J[314] += dqdci;              /* dwdot[HO2]/d[CH2O] */
    J[322] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[323] -= dqdci;              /* dwdot[CH3O]/d[CH2O] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[3];
    J[333] -= dqdci;              /* dwdot[O2]/d[CH3O] */
    J[336] += dqdci;              /* dwdot[HO2]/d[CH3O] */
    J[344] += dqdci;              /* dwdot[CH2O]/d[CH3O] */
    J[345] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[468] += dqdT;               /* dwdot[HO2]/dT */
    J[476] += dqdT;               /* dwdot[CH2O]/dT */
    J[477] -= dqdT;               /* dwdot[CH3O]/dT */

    /*reaction 84: C2H5 + O2 <=> HO2 + C2H4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[17];
    k_f = prefactor_units[83] * fwd_A[83]
                * exp(fwd_beta[83] * tc[0] - activation_units[83] * fwd_Ea[83] * invT);
    dlnkfdT = fwd_beta[83] * invT + activation_units[83] * fwd_Ea[83] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[16];
    Kc = exp(g_RT[3] - g_RT[6] - g_RT[16] + g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[17]) + (h_RT[6] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    wdot[16] += q; /* C2H4 */
    wdot[17] -= q; /* C2H5 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[17];
    J[69] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[72] += dqdci;               /* dwdot[HO2]/d[O2] */
    J[82] += dqdci;               /* dwdot[C2H4]/d[O2] */
    J[83] -= dqdci;               /* dwdot[C2H5]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[16];
    J[135] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[138] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[148] += dqdci;              /* dwdot[C2H4]/d[HO2] */
    J[149] -= dqdci;              /* dwdot[C2H5]/d[HO2] */
    /* d()/d[C2H4] */
    dqdci =  - k_r*sc[6];
    J[355] -= dqdci;              /* dwdot[O2]/d[C2H4] */
    J[358] += dqdci;              /* dwdot[HO2]/d[C2H4] */
    J[368] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[369] -= dqdci;              /* dwdot[C2H5]/d[C2H4] */
    /* d()/d[C2H5] */
    dqdci =  + k_f*sc[3];
    J[377] -= dqdci;              /* dwdot[O2]/d[C2H5] */
    J[380] += dqdci;              /* dwdot[HO2]/d[C2H5] */
    J[390] += dqdci;              /* dwdot[C2H4]/d[C2H5] */
    J[391] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
    /* d()/dT */
    J[465] -= dqdT;               /* dwdot[O2]/dT */
    J[468] += dqdT;               /* dwdot[HO2]/dT */
    J[478] += dqdT;               /* dwdot[C2H4]/dT */
    J[479] -= dqdT;               /* dwdot[C2H5]/dT */

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
    *LENIMC = 86;}


void egtransetLENRMC(int* LENRMC ) {
    *LENRMC = 9114;}


void egtransetNO(int* NO ) {
    *NO = 4;}


void egtransetKK(int* KK ) {
    *KK = 21;}


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
    WT[7] = 1.40270900E+01;
    WT[8] = 1.40270900E+01;
    WT[9] = 1.50350600E+01;
    WT[10] = 1.60430300E+01;
    WT[11] = 2.80105500E+01;
    WT[12] = 4.40099500E+01;
    WT[13] = 2.90185200E+01;
    WT[14] = 3.00264900E+01;
    WT[15] = 3.10344600E+01;
    WT[16] = 2.80541800E+01;
    WT[17] = 2.90621500E+01;
    WT[18] = 3.00701200E+01;
    WT[19] = 2.80134000E+01;
    WT[20] = 3.99480000E+01;
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
    EPS[7] = 1.44000000E+02;
    EPS[8] = 1.44000000E+02;
    EPS[9] = 1.44000000E+02;
    EPS[10] = 1.41400000E+02;
    EPS[11] = 9.81000000E+01;
    EPS[12] = 2.44000000E+02;
    EPS[13] = 4.98000000E+02;
    EPS[14] = 4.98000000E+02;
    EPS[15] = 4.17000000E+02;
    EPS[16] = 2.80800000E+02;
    EPS[17] = 2.52300000E+02;
    EPS[18] = 2.52300000E+02;
    EPS[19] = 9.75300000E+01;
    EPS[20] = 1.36500000E+02;
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
    SIG[7] = 3.80000000E+00;
    SIG[8] = 3.80000000E+00;
    SIG[9] = 3.80000000E+00;
    SIG[10] = 3.74600000E+00;
    SIG[11] = 3.65000000E+00;
    SIG[12] = 3.76300000E+00;
    SIG[13] = 3.59000000E+00;
    SIG[14] = 3.59000000E+00;
    SIG[15] = 3.69000000E+00;
    SIG[16] = 3.97100000E+00;
    SIG[17] = 4.30200000E+00;
    SIG[18] = 4.30200000E+00;
    SIG[19] = 3.62100000E+00;
    SIG[20] = 3.33000000E+00;
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
    DIP[9] = 0.00000000E+00;
    DIP[10] = 0.00000000E+00;
    DIP[11] = 0.00000000E+00;
    DIP[12] = 0.00000000E+00;
    DIP[13] = 0.00000000E+00;
    DIP[14] = 0.00000000E+00;
    DIP[15] = 1.70000000E+00;
    DIP[16] = 0.00000000E+00;
    DIP[17] = 0.00000000E+00;
    DIP[18] = 0.00000000E+00;
    DIP[19] = 0.00000000E+00;
    DIP[20] = 0.00000000E+00;
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
    POL[8] = 0.00000000E+00;
    POL[9] = 0.00000000E+00;
    POL[10] = 2.60000000E+00;
    POL[11] = 1.95000000E+00;
    POL[12] = 2.65000000E+00;
    POL[13] = 0.00000000E+00;
    POL[14] = 0.00000000E+00;
    POL[15] = 0.00000000E+00;
    POL[16] = 0.00000000E+00;
    POL[17] = 0.00000000E+00;
    POL[18] = 0.00000000E+00;
    POL[19] = 1.76000000E+00;
    POL[20] = 0.00000000E+00;
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
    ZROT[7] = 0.00000000E+00;
    ZROT[8] = 0.00000000E+00;
    ZROT[9] = 0.00000000E+00;
    ZROT[10] = 1.30000000E+01;
    ZROT[11] = 1.80000000E+00;
    ZROT[12] = 2.10000000E+00;
    ZROT[13] = 0.00000000E+00;
    ZROT[14] = 2.00000000E+00;
    ZROT[15] = 2.00000000E+00;
    ZROT[16] = 1.50000000E+00;
    ZROT[17] = 1.50000000E+00;
    ZROT[18] = 1.50000000E+00;
    ZROT[19] = 4.00000000E+00;
    ZROT[20] = 0.00000000E+00;
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
    NLIN[7] = 1;
    NLIN[8] = 1;
    NLIN[9] = 1;
    NLIN[10] = 2;
    NLIN[11] = 1;
    NLIN[12] = 1;
    NLIN[13] = 2;
    NLIN[14] = 2;
    NLIN[15] = 2;
    NLIN[16] = 2;
    NLIN[17] = 2;
    NLIN[18] = 2;
    NLIN[19] = 1;
    NLIN[20] = 0;
}


/*Poly fits for the viscosities, dim NO*KK */
void egtransetCOFETA(amrex::Real* COFETA) {
    COFETA[0] = -1.38347699E+01;
    COFETA[1] = 1.00106621E+00;
    COFETA[2] = -4.98105694E-02;
    COFETA[3] = 2.31450475E-03;
    COFETA[4] = -2.04078397E+01;
    COFETA[5] = 3.65436395E+00;
    COFETA[6] = -3.98339635E-01;
    COFETA[7] = 1.75883009E-02;
    COFETA[8] = -1.50926240E+01;
    COFETA[9] = 1.92606504E+00;
    COFETA[10] = -1.73487476E-01;
    COFETA[11] = 7.82572931E-03;
    COFETA[12] = -1.71618309E+01;
    COFETA[13] = 2.68036374E+00;
    COFETA[14] = -2.72570227E-01;
    COFETA[15] = 1.21650964E-02;
    COFETA[16] = -1.50620763E+01;
    COFETA[17] = 1.92606504E+00;
    COFETA[18] = -1.73487476E-01;
    COFETA[19] = 7.82572931E-03;
    COFETA[20] = -1.05420863E+01;
    COFETA[21] = -1.37777096E+00;
    COFETA[22] = 4.20502308E-01;
    COFETA[23] = -2.40627230E-02;
    COFETA[24] = -1.71463238E+01;
    COFETA[25] = 2.68036374E+00;
    COFETA[26] = -2.72570227E-01;
    COFETA[27] = 1.21650964E-02;
    COFETA[28] = -2.02663469E+01;
    COFETA[29] = 3.63241793E+00;
    COFETA[30] = -3.95581049E-01;
    COFETA[31] = 1.74725495E-02;
    COFETA[32] = -2.02663469E+01;
    COFETA[33] = 3.63241793E+00;
    COFETA[34] = -3.95581049E-01;
    COFETA[35] = 1.74725495E-02;
    COFETA[36] = -2.02316497E+01;
    COFETA[37] = 3.63241793E+00;
    COFETA[38] = -3.95581049E-01;
    COFETA[39] = 1.74725495E-02;
    COFETA[40] = -2.00094664E+01;
    COFETA[41] = 3.57220167E+00;
    COFETA[42] = -3.87936446E-01;
    COFETA[43] = 1.71483254E-02;
    COFETA[44] = -1.66188336E+01;
    COFETA[45] = 2.40307799E+00;
    COFETA[46] = -2.36167638E-01;
    COFETA[47] = 1.05714061E-02;
    COFETA[48] = -2.40014975E+01;
    COFETA[49] = 5.14359547E+00;
    COFETA[50] = -5.74269731E-01;
    COFETA[51] = 2.44937679E-02;
    COFETA[52] = -1.98501306E+01;
    COFETA[53] = 2.69480162E+00;
    COFETA[54] = -1.65880845E-01;
    COFETA[55] = 3.14504769E-03;
    COFETA[56] = -1.98330577E+01;
    COFETA[57] = 2.69480162E+00;
    COFETA[58] = -1.65880845E-01;
    COFETA[59] = 3.14504769E-03;
    COFETA[60] = -1.99945919E+01;
    COFETA[61] = 2.86923313E+00;
    COFETA[62] = -2.03325661E-01;
    COFETA[63] = 5.39056989E-03;
    COFETA[64] = -2.50655444E+01;
    COFETA[65] = 5.33982977E+00;
    COFETA[66] = -5.89982992E-01;
    COFETA[67] = 2.47780650E-02;
    COFETA[68] = -2.46581414E+01;
    COFETA[69] = 5.19497183E+00;
    COFETA[70] = -5.78827339E-01;
    COFETA[71] = 2.46050921E-02;
    COFETA[72] = -2.46410937E+01;
    COFETA[73] = 5.19497183E+00;
    COFETA[74] = -5.78827339E-01;
    COFETA[75] = 2.46050921E-02;
    COFETA[76] = -1.65695594E+01;
    COFETA[77] = 2.39056562E+00;
    COFETA[78] = -2.34558144E-01;
    COFETA[79] = 1.05024037E-02;
    COFETA[80] = -1.90422907E+01;
    COFETA[81] = 3.47025711E+00;
    COFETA[82] = -3.75102111E-01;
    COFETA[83] = 1.66086076E-02;
}


/*Poly fits for the conductivities, dim NO*KK */
void egtransetCOFLAM(amrex::Real* COFLAM) {
    COFLAM[0] = 9.24084392E+00;
    COFLAM[1] = -4.69568931E-01;
    COFLAM[2] = 1.15980279E-01;
    COFLAM[3] = -2.61033830E-03;
    COFLAM[4] = -8.57929284E-01;
    COFLAM[5] = 3.65436395E+00;
    COFLAM[6] = -3.98339635E-01;
    COFLAM[7] = 1.75883009E-02;
    COFLAM[8] = 1.69267361E+00;
    COFLAM[9] = 1.92606504E+00;
    COFLAM[10] = -1.73487476E-01;
    COFLAM[11] = 7.82572931E-03;
    COFLAM[12] = -1.93718739E+00;
    COFLAM[13] = 2.89110219E+00;
    COFLAM[14] = -2.71096923E-01;
    COFLAM[15] = 1.15344907E-02;
    COFLAM[16] = 1.41666185E+01;
    COFLAM[17] = -3.24630496E+00;
    COFLAM[18] = 5.33878183E-01;
    COFLAM[19] = -2.32905165E-02;
    COFLAM[20] = 2.33729817E+01;
    COFLAM[21] = -8.96536433E+00;
    COFLAM[22] = 1.52828665E+00;
    COFLAM[23] = -7.58551778E-02;
    COFLAM[24] = -1.12960913E+00;
    COFLAM[25] = 2.34014305E+00;
    COFLAM[26] = -1.63245030E-01;
    COFLAM[27] = 5.80319600E-03;
    COFLAM[28] = 1.29177902E+01;
    COFLAM[29] = -3.73745535E+00;
    COFLAM[30] = 7.15831021E-01;
    COFLAM[31] = -3.63846910E-02;
    COFLAM[32] = 1.89383266E+01;
    COFLAM[33] = -6.51018128E+00;
    COFLAM[34] = 1.13292060E+00;
    COFLAM[35] = -5.69603226E-02;
    COFLAM[36] = 1.39937901E+01;
    COFLAM[37] = -4.64256494E+00;
    COFLAM[38] = 9.07728674E-01;
    COFLAM[39] = -4.77274469E-02;
    COFLAM[40] = 1.33091602E+01;
    COFLAM[41] = -4.96140261E+00;
    COFLAM[42] = 1.03295505E+00;
    COFLAM[43] = -5.63420090E-02;
    COFLAM[44] = 1.18777264E+01;
    COFLAM[45] = -3.15463949E+00;
    COFLAM[46] = 6.01973268E-01;
    COFLAM[47] = -3.03211261E-02;
    COFLAM[48] = -1.13649314E+01;
    COFLAM[49] = 5.88177395E+00;
    COFLAM[50] = -5.68651819E-01;
    COFLAM[51] = 2.03561485E-02;
    COFLAM[52] = 6.30243508E+00;
    COFLAM[53] = -2.22810801E+00;
    COFLAM[54] = 6.37340514E-01;
    COFLAM[55] = -3.81056018E-02;
    COFLAM[56] = 5.39305086E+00;
    COFLAM[57] = -2.39312375E+00;
    COFLAM[58] = 7.39585221E-01;
    COFLAM[59] = -4.58435589E-02;
    COFLAM[60] = -6.14588187E+00;
    COFLAM[61] = 2.47428873E+00;
    COFLAM[62] = 6.43999571E-02;
    COFLAM[63] = -1.45368336E-02;
    COFLAM[64] = -1.46152839E+01;
    COFLAM[65] = 6.36251406E+00;
    COFLAM[66] = -5.03832130E-01;
    COFLAM[67] = 1.26121050E-02;
    COFLAM[68] = -8.95009705E+00;
    COFLAM[69] = 4.02515080E+00;
    COFLAM[70] = -1.84063946E-01;
    COFLAM[71] = -1.94054752E-03;
    COFLAM[72] = -1.09902209E+01;
    COFLAM[73] = 4.70647707E+00;
    COFLAM[74] = -2.52272495E-01;
    COFLAM[75] = 1.75193258E-04;
    COFLAM[76] = 1.29306158E+01;
    COFLAM[77] = -3.52817362E+00;
    COFLAM[78] = 6.45499013E-01;
    COFLAM[79] = -3.19375299E-02;
    COFLAM[80] = -3.17202048E+00;
    COFLAM[81] = 3.47025711E+00;
    COFLAM[82] = -3.75102111E-01;
    COFLAM[83] = 1.66086076E-02;
}


/*Poly fits for the diffusion coefficients, dim NO*KK*KK */
void egtransetCOFD(amrex::Real* COFD) {
    COFD[0] = -1.03270606E+01;
    COFD[1] = 2.19285409E+00;
    COFD[2] = -7.54492786E-02;
    COFD[3] = 3.51398213E-03;
    COFD[4] = -1.14366381E+01;
    COFD[5] = 2.78323501E+00;
    COFD[6] = -1.51214064E-01;
    COFD[7] = 6.75150012E-03;
    COFD[8] = -1.09595712E+01;
    COFD[9] = 2.30836460E+00;
    COFD[10] = -8.76339315E-02;
    COFD[11] = 3.90878445E-03;
    COFD[12] = -1.18988955E+01;
    COFD[13] = 2.57507000E+00;
    COFD[14] = -1.24033737E-01;
    COFD[15] = 5.56694959E-03;
    COFD[16] = -1.09628982E+01;
    COFD[17] = 2.30836460E+00;
    COFD[18] = -8.76339315E-02;
    COFD[19] = 3.90878445E-03;
    COFD[20] = -1.71982995E+01;
    COFD[21] = 4.63881404E+00;
    COFD[22] = -3.86139633E-01;
    COFD[23] = 1.66955081E-02;
    COFD[24] = -1.18998012E+01;
    COFD[25] = 2.57507000E+00;
    COFD[26] = -1.24033737E-01;
    COFD[27] = 5.56694959E-03;
    COFD[28] = -1.25098960E+01;
    COFD[29] = 2.77873601E+00;
    COFD[30] = -1.50637360E-01;
    COFD[31] = 6.72684281E-03;
    COFD[32] = -1.25098960E+01;
    COFD[33] = 2.77873601E+00;
    COFD[34] = -1.50637360E-01;
    COFD[35] = 6.72684281E-03;
    COFD[36] = -1.25141260E+01;
    COFD[37] = 2.77873601E+00;
    COFD[38] = -1.50637360E-01;
    COFD[39] = 6.72684281E-03;
    COFD[40] = -1.24693568E+01;
    COFD[41] = 2.76686648E+00;
    COFD[42] = -1.49120141E-01;
    COFD[43] = 6.66220432E-03;
    COFD[44] = -1.17159737E+01;
    COFD[45] = 2.48123210E+00;
    COFD[46] = -1.11322604E-01;
    COFD[47] = 4.99282389E-03;
    COFD[48] = -1.37794315E+01;
    COFD[49] = 3.23973858E+00;
    COFD[50] = -2.09989036E-01;
    COFD[51] = 9.27667906E-03;
    COFD[52] = -1.60517370E+01;
    COFD[53] = 4.11188603E+00;
    COFD[54] = -3.21540884E-01;
    COFD[55] = 1.40482564E-02;
    COFD[56] = -1.60528285E+01;
    COFD[57] = 4.11188603E+00;
    COFD[58] = -3.21540884E-01;
    COFD[59] = 1.40482564E-02;
    COFD[60] = -1.58456300E+01;
    COFD[61] = 4.02074783E+00;
    COFD[62] = -3.10018522E-01;
    COFD[63] = 1.35599552E-02;
    COFD[64] = -1.42229194E+01;
    COFD[65] = 3.38669384E+00;
    COFD[66] = -2.28784122E-01;
    COFD[67] = 1.00790953E-02;
    COFD[68] = -1.39913897E+01;
    COFD[69] = 3.26384506E+00;
    COFD[70] = -2.12947087E-01;
    COFD[71] = 9.39743888E-03;
    COFD[72] = -1.39924781E+01;
    COFD[73] = 3.26384506E+00;
    COFD[74] = -2.12947087E-01;
    COFD[75] = 9.39743888E-03;
    COFD[76] = -1.16906297E+01;
    COFD[77] = 2.47469981E+00;
    COFD[78] = -1.10436257E-01;
    COFD[79] = 4.95273813E-03;
    COFD[80] = -1.23130152E+01;
    COFD[81] = 2.74418790E+00;
    COFD[82] = -1.46230156E-01;
    COFD[83] = 6.53948886E-03;
    COFD[84] = -1.14366381E+01;
    COFD[85] = 2.78323501E+00;
    COFD[86] = -1.51214064E-01;
    COFD[87] = 6.75150012E-03;
    COFD[88] = -1.47968712E+01;
    COFD[89] = 4.23027636E+00;
    COFD[90] = -3.36139991E-01;
    COFD[91] = 1.46507621E-02;
    COFD[92] = -1.34230272E+01;
    COFD[93] = 3.48624238E+00;
    COFD[94] = -2.41554467E-01;
    COFD[95] = 1.06263545E-02;
    COFD[96] = -1.46550083E+01;
    COFD[97] = 3.83606243E+00;
    COFD[98] = -2.86076532E-01;
    COFD[99] = 1.25205829E-02;
    COFD[100] = -1.34247866E+01;
    COFD[101] = 3.48624238E+00;
    COFD[102] = -2.41554467E-01;
    COFD[103] = 1.06263545E-02;
    COFD[104] = -1.95739570E+01;
    COFD[105] = 5.61113230E+00;
    COFD[106] = -4.90190187E-01;
    COFD[107] = 2.03260675E-02;
    COFD[108] = -1.46554748E+01;
    COFD[109] = 3.83606243E+00;
    COFD[110] = -2.86076532E-01;
    COFD[111] = 1.25205829E-02;
    COFD[112] = -1.57972369E+01;
    COFD[113] = 4.22225052E+00;
    COFD[114] = -3.35156428E-01;
    COFD[115] = 1.46104855E-02;
    COFD[116] = -1.57972369E+01;
    COFD[117] = 4.22225052E+00;
    COFD[118] = -3.35156428E-01;
    COFD[119] = 1.46104855E-02;
    COFD[120] = -1.57994893E+01;
    COFD[121] = 4.22225052E+00;
    COFD[122] = -3.35156428E-01;
    COFD[123] = 1.46104855E-02;
    COFD[124] = -1.57199037E+01;
    COFD[125] = 4.19936335E+00;
    COFD[126] = -3.32311009E-01;
    COFD[127] = 1.44921003E-02;
    COFD[128] = -1.43151174E+01;
    COFD[129] = 3.68038508E+00;
    COFD[130] = -2.65779346E-01;
    COFD[131] = 1.16360771E-02;
    COFD[132] = -1.76147026E+01;
    COFD[133] = 4.86049500E+00;
    COFD[134] = -4.12200578E-01;
    COFD[135] = 1.77160971E-02;
    COFD[136] = -1.97544450E+01;
    COFD[137] = 5.56931926E+00;
    COFD[138] = -4.89105511E-01;
    COFD[139] = 2.04493129E-02;
    COFD[140] = -1.97550088E+01;
    COFD[141] = 5.56931926E+00;
    COFD[142] = -4.89105511E-01;
    COFD[143] = 2.04493129E-02;
    COFD[144] = -1.92718582E+01;
    COFD[145] = 5.41172124E+00;
    COFD[146] = -4.73213887E-01;
    COFD[147] = 1.99405473E-02;
    COFD[148] = -1.82251914E+01;
    COFD[149] = 5.05237312E+00;
    COFD[150] = -4.35182396E-01;
    COFD[151] = 1.86363074E-02;
    COFD[152] = -1.79339327E+01;
    COFD[153] = 4.91373893E+00;
    COFD[154] = -4.18747629E-01;
    COFD[155] = 1.79856610E-02;
    COFD[156] = -1.79344949E+01;
    COFD[157] = 4.91373893E+00;
    COFD[158] = -4.18747629E-01;
    COFD[159] = 1.79856610E-02;
    COFD[160] = -1.42894441E+01;
    COFD[161] = 3.67490723E+00;
    COFD[162] = -2.65114792E-01;
    COFD[163] = 1.16092671E-02;
    COFD[164] = -1.54738604E+01;
    COFD[165] = 4.15765300E+00;
    COFD[166] = -3.27126237E-01;
    COFD[167] = 1.42762611E-02;
    COFD[168] = -1.09595712E+01;
    COFD[169] = 2.30836460E+00;
    COFD[170] = -8.76339315E-02;
    COFD[171] = 3.90878445E-03;
    COFD[172] = -1.34230272E+01;
    COFD[173] = 3.48624238E+00;
    COFD[174] = -2.41554467E-01;
    COFD[175] = 1.06263545E-02;
    COFD[176] = -1.32093628E+01;
    COFD[177] = 2.90778936E+00;
    COFD[178] = -1.67388544E-01;
    COFD[179] = 7.45220609E-03;
    COFD[180] = -1.43139231E+01;
    COFD[181] = 3.17651319E+00;
    COFD[182] = -2.02028974E-01;
    COFD[183] = 8.94232502E-03;
    COFD[184] = -1.32244035E+01;
    COFD[185] = 2.90778936E+00;
    COFD[186] = -1.67388544E-01;
    COFD[187] = 7.45220609E-03;
    COFD[188] = -1.94093572E+01;
    COFD[189] = 5.16013126E+00;
    COFD[190] = -4.46824543E-01;
    COFD[191] = 1.90464887E-02;
    COFD[192] = -1.43190389E+01;
    COFD[193] = 3.17651319E+00;
    COFD[194] = -2.02028974E-01;
    COFD[195] = 8.94232502E-03;
    COFD[196] = -1.50584249E+01;
    COFD[197] = 3.47945612E+00;
    COFD[198] = -2.40703722E-01;
    COFD[199] = 1.05907441E-02;
    COFD[200] = -1.50584249E+01;
    COFD[201] = 3.47945612E+00;
    COFD[202] = -2.40703722E-01;
    COFD[203] = 1.05907441E-02;
    COFD[204] = -1.50766130E+01;
    COFD[205] = 3.47945612E+00;
    COFD[206] = -2.40703722E-01;
    COFD[207] = 1.05907441E-02;
    COFD[208] = -1.50270339E+01;
    COFD[209] = 3.46140064E+00;
    COFD[210] = -2.38440092E-01;
    COFD[211] = 1.04960087E-02;
    COFD[212] = -1.40999008E+01;
    COFD[213] = 3.08120012E+00;
    COFD[214] = -1.89629903E-01;
    COFD[215] = 8.40361952E-03;
    COFD[216] = -1.70534856E+01;
    COFD[217] = 4.14240922E+00;
    COFD[218] = -3.25239774E-01;
    COFD[219] = 1.41980687E-02;
    COFD[220] = -1.94313116E+01;
    COFD[221] = 5.02567894E+00;
    COFD[222] = -4.32045169E-01;
    COFD[223] = 1.85132214E-02;
    COFD[224] = -1.94373127E+01;
    COFD[225] = 5.02567894E+00;
    COFD[226] = -4.32045169E-01;
    COFD[227] = 1.85132214E-02;
    COFD[228] = -1.88179418E+01;
    COFD[229] = 4.79683898E+00;
    COFD[230] = -4.04829719E-01;
    COFD[231] = 1.74325475E-02;
    COFD[232] = -1.74792112E+01;
    COFD[233] = 4.29676909E+00;
    COFD[234] = -3.44085306E-01;
    COFD[235] = 1.49671135E-02;
    COFD[236] = -1.72496634E+01;
    COFD[237] = 4.17889917E+00;
    COFD[238] = -3.29752510E-01;
    COFD[239] = 1.43850275E-02;
    COFD[240] = -1.72556499E+01;
    COFD[241] = 4.17889917E+00;
    COFD[242] = -3.29752510E-01;
    COFD[243] = 1.43850275E-02;
    COFD[244] = -1.40756935E+01;
    COFD[245] = 3.07549274E+00;
    COFD[246] = -1.88889344E-01;
    COFD[247] = 8.37152866E-03;
    COFD[248] = -1.49610527E+01;
    COFD[249] = 3.41988961E+00;
    COFD[250] = -2.33128386E-01;
    COFD[251] = 1.02689994E-02;
    COFD[252] = -1.18988955E+01;
    COFD[253] = 2.57507000E+00;
    COFD[254] = -1.24033737E-01;
    COFD[255] = 5.56694959E-03;
    COFD[256] = -1.46550083E+01;
    COFD[257] = 3.83606243E+00;
    COFD[258] = -2.86076532E-01;
    COFD[259] = 1.25205829E-02;
    COFD[260] = -1.43139231E+01;
    COFD[261] = 3.17651319E+00;
    COFD[262] = -2.02028974E-01;
    COFD[263] = 8.94232502E-03;
    COFD[264] = -1.55511344E+01;
    COFD[265] = 3.48070094E+00;
    COFD[266] = -2.40859499E-01;
    COFD[267] = 1.05972514E-02;
    COFD[268] = -1.43340796E+01;
    COFD[269] = 3.17651319E+00;
    COFD[270] = -2.02028974E-01;
    COFD[271] = 8.94232502E-03;
    COFD[272] = -2.12652533E+01;
    COFD[273] = 5.59961818E+00;
    COFD[274] = -4.91624858E-01;
    COFD[275] = 2.05035550E-02;
    COFD[276] = -1.55588279E+01;
    COFD[277] = 3.48070094E+00;
    COFD[278] = -2.40859499E-01;
    COFD[279] = 1.05972514E-02;
    COFD[280] = -1.63254691E+01;
    COFD[281] = 3.82388595E+00;
    COFD[282] = -2.84480724E-01;
    COFD[283] = 1.24506311E-02;
    COFD[284] = -1.63254691E+01;
    COFD[285] = 3.82388595E+00;
    COFD[286] = -2.84480724E-01;
    COFD[287] = 1.24506311E-02;
    COFD[288] = -1.63493345E+01;
    COFD[289] = 3.82388595E+00;
    COFD[290] = -2.84480724E-01;
    COFD[291] = 1.24506311E-02;
    COFD[292] = -1.62724462E+01;
    COFD[293] = 3.79163564E+00;
    COFD[294] = -2.80257365E-01;
    COFD[295] = 1.22656902E-02;
    COFD[296] = -1.52721107E+01;
    COFD[297] = 3.36790500E+00;
    COFD[298] = -2.26321740E-01;
    COFD[299] = 9.97135055E-03;
    COFD[300] = -1.84688406E+01;
    COFD[301] = 4.49330851E+00;
    COFD[302] = -3.68208715E-01;
    COFD[303] = 1.59565402E-02;
    COFD[304] = -2.08204449E+01;
    COFD[305] = 5.35267674E+00;
    COFD[306] = -4.69010505E-01;
    COFD[307] = 1.98979152E-02;
    COFD[308] = -2.08293255E+01;
    COFD[309] = 5.35267674E+00;
    COFD[310] = -4.69010505E-01;
    COFD[311] = 1.98979152E-02;
    COFD[312] = -2.04928958E+01;
    COFD[313] = 5.22397933E+00;
    COFD[314] = -4.54138171E-01;
    COFD[315] = 1.93249285E-02;
    COFD[316] = -1.89544778E+01;
    COFD[317] = 4.68595732E+00;
    COFD[318] = -3.91842840E-01;
    COFD[319] = 1.69262542E-02;
    COFD[320] = -1.86335932E+01;
    COFD[321] = 4.53572533E+00;
    COFD[322] = -3.73386925E-01;
    COFD[323] = 1.61678881E-02;
    COFD[324] = -1.86424545E+01;
    COFD[325] = 4.53572533E+00;
    COFD[326] = -3.73386925E-01;
    COFD[327] = 1.61678881E-02;
    COFD[328] = -1.52414485E+01;
    COFD[329] = 3.35922578E+00;
    COFD[330] = -2.25181399E-01;
    COFD[331] = 9.92132878E-03;
    COFD[332] = -1.62380075E+01;
    COFD[333] = 3.72612300E+00;
    COFD[334] = -2.71663673E-01;
    COFD[335] = 1.18889643E-02;
    COFD[336] = -1.09628982E+01;
    COFD[337] = 2.30836460E+00;
    COFD[338] = -8.76339315E-02;
    COFD[339] = 3.90878445E-03;
    COFD[340] = -1.34247866E+01;
    COFD[341] = 3.48624238E+00;
    COFD[342] = -2.41554467E-01;
    COFD[343] = 1.06263545E-02;
    COFD[344] = -1.32244035E+01;
    COFD[345] = 2.90778936E+00;
    COFD[346] = -1.67388544E-01;
    COFD[347] = 7.45220609E-03;
    COFD[348] = -1.43340796E+01;
    COFD[349] = 3.17651319E+00;
    COFD[350] = -2.02028974E-01;
    COFD[351] = 8.94232502E-03;
    COFD[352] = -1.32399106E+01;
    COFD[353] = 2.90778936E+00;
    COFD[354] = -1.67388544E-01;
    COFD[355] = 7.45220609E-03;
    COFD[356] = -1.94253036E+01;
    COFD[357] = 5.16013126E+00;
    COFD[358] = -4.46824543E-01;
    COFD[359] = 1.90464887E-02;
    COFD[360] = -1.43394069E+01;
    COFD[361] = 3.17651319E+00;
    COFD[362] = -2.02028974E-01;
    COFD[363] = 8.94232502E-03;
    COFD[364] = -1.50724636E+01;
    COFD[365] = 3.47945612E+00;
    COFD[366] = -2.40703722E-01;
    COFD[367] = 1.05907441E-02;
    COFD[368] = -1.50724636E+01;
    COFD[369] = 3.47945612E+00;
    COFD[370] = -2.40703722E-01;
    COFD[371] = 1.05907441E-02;
    COFD[372] = -1.50911794E+01;
    COFD[373] = 3.47945612E+00;
    COFD[374] = -2.40703722E-01;
    COFD[375] = 1.05907441E-02;
    COFD[376] = -1.50420953E+01;
    COFD[377] = 3.46140064E+00;
    COFD[378] = -2.38440092E-01;
    COFD[379] = 1.04960087E-02;
    COFD[380] = -1.41191261E+01;
    COFD[381] = 3.08120012E+00;
    COFD[382] = -1.89629903E-01;
    COFD[383] = 8.40361952E-03;
    COFD[384] = -1.70757047E+01;
    COFD[385] = 4.14240922E+00;
    COFD[386] = -3.25239774E-01;
    COFD[387] = 1.41980687E-02;
    COFD[388] = -1.94507876E+01;
    COFD[389] = 5.02567894E+00;
    COFD[390] = -4.32045169E-01;
    COFD[391] = 1.85132214E-02;
    COFD[392] = -1.94570287E+01;
    COFD[393] = 5.02567894E+00;
    COFD[394] = -4.32045169E-01;
    COFD[395] = 1.85132214E-02;
    COFD[396] = -1.88378874E+01;
    COFD[397] = 4.79683898E+00;
    COFD[398] = -4.04829719E-01;
    COFD[399] = 1.74325475E-02;
    COFD[400] = -1.74984476E+01;
    COFD[401] = 4.29676909E+00;
    COFD[402] = -3.44085306E-01;
    COFD[403] = 1.49671135E-02;
    COFD[404] = -1.72691500E+01;
    COFD[405] = 4.17889917E+00;
    COFD[406] = -3.29752510E-01;
    COFD[407] = 1.43850275E-02;
    COFD[408] = -1.72753760E+01;
    COFD[409] = 4.17889917E+00;
    COFD[410] = -3.29752510E-01;
    COFD[411] = 1.43850275E-02;
    COFD[412] = -1.40949196E+01;
    COFD[413] = 3.07549274E+00;
    COFD[414] = -1.88889344E-01;
    COFD[415] = 8.37152866E-03;
    COFD[416] = -1.49826725E+01;
    COFD[417] = 3.41988961E+00;
    COFD[418] = -2.33128386E-01;
    COFD[419] = 1.02689994E-02;
    COFD[420] = -1.71982995E+01;
    COFD[421] = 4.63881404E+00;
    COFD[422] = -3.86139633E-01;
    COFD[423] = 1.66955081E-02;
    COFD[424] = -1.95739570E+01;
    COFD[425] = 5.61113230E+00;
    COFD[426] = -4.90190187E-01;
    COFD[427] = 2.03260675E-02;
    COFD[428] = -1.94093572E+01;
    COFD[429] = 5.16013126E+00;
    COFD[430] = -4.46824543E-01;
    COFD[431] = 1.90464887E-02;
    COFD[432] = -2.12652533E+01;
    COFD[433] = 5.59961818E+00;
    COFD[434] = -4.91624858E-01;
    COFD[435] = 2.05035550E-02;
    COFD[436] = -1.94253036E+01;
    COFD[437] = 5.16013126E+00;
    COFD[438] = -4.46824543E-01;
    COFD[439] = 1.90464887E-02;
    COFD[440] = -1.19157919E+01;
    COFD[441] = 9.28955130E-01;
    COFD[442] = 2.42107090E-01;
    COFD[443] = -1.59823963E-02;
    COFD[444] = -2.06463744E+01;
    COFD[445] = 5.41688482E+00;
    COFD[446] = -4.73387188E-01;
    COFD[447] = 1.99280175E-02;
    COFD[448] = -2.12639214E+01;
    COFD[449] = 5.61184117E+00;
    COFD[450] = -4.90532156E-01;
    COFD[451] = 2.03507922E-02;
    COFD[452] = -2.12639214E+01;
    COFD[453] = 5.61184117E+00;
    COFD[454] = -4.90532156E-01;
    COFD[455] = 2.03507922E-02;
    COFD[456] = -2.12831323E+01;
    COFD[457] = 5.61184117E+00;
    COFD[458] = -4.90532156E-01;
    COFD[459] = 2.03507922E-02;
    COFD[460] = -2.14087397E+01;
    COFD[461] = 5.57282008E+00;
    COFD[462] = -4.76690890E-01;
    COFD[463] = 1.94000719E-02;
    COFD[464] = -2.11388331E+01;
    COFD[465] = 5.55529675E+00;
    COFD[466] = -4.87942518E-01;
    COFD[467] = 2.04249054E-02;
    COFD[468] = -2.07653719E+01;
    COFD[469] = 5.01092022E+00;
    COFD[470] = -3.77985635E-01;
    COFD[471] = 1.40968645E-02;
    COFD[472] = -1.77498543E+01;
    COFD[473] = 3.57475686E+00;
    COFD[474] = -1.56396297E-01;
    COFD[475] = 3.12157721E-03;
    COFD[476] = -1.77563250E+01;
    COFD[477] = 3.57475686E+00;
    COFD[478] = -1.56396297E-01;
    COFD[479] = 3.12157721E-03;
    COFD[480] = -1.65295288E+01;
    COFD[481] = 2.97569206E+00;
    COFD[482] = -6.75652842E-02;
    COFD[483] = -1.08648422E-03;
    COFD[484] = -2.08812333E+01;
    COFD[485] = 5.08859217E+00;
    COFD[486] = -3.90525428E-01;
    COFD[487] = 1.47376395E-02;
    COFD[488] = -2.12597312E+01;
    COFD[489] = 5.24930667E+00;
    COFD[490] = -4.17435088E-01;
    COFD[491] = 1.61434424E-02;
    COFD[492] = -2.12661865E+01;
    COFD[493] = 5.24930667E+00;
    COFD[494] = -4.17435088E-01;
    COFD[495] = 1.61434424E-02;
    COFD[496] = -2.10643259E+01;
    COFD[497] = 5.53614847E+00;
    COFD[498] = -4.86046736E-01;
    COFD[499] = 2.03659188E-02;
    COFD[500] = -2.12755888E+01;
    COFD[501] = 5.60381989E+00;
    COFD[502] = -4.91225459E-01;
    COFD[503] = 2.04487844E-02;
    COFD[504] = -1.18998012E+01;
    COFD[505] = 2.57507000E+00;
    COFD[506] = -1.24033737E-01;
    COFD[507] = 5.56694959E-03;
    COFD[508] = -1.46554748E+01;
    COFD[509] = 3.83606243E+00;
    COFD[510] = -2.86076532E-01;
    COFD[511] = 1.25205829E-02;
    COFD[512] = -1.43190389E+01;
    COFD[513] = 3.17651319E+00;
    COFD[514] = -2.02028974E-01;
    COFD[515] = 8.94232502E-03;
    COFD[516] = -1.55588279E+01;
    COFD[517] = 3.48070094E+00;
    COFD[518] = -2.40859499E-01;
    COFD[519] = 1.05972514E-02;
    COFD[520] = -1.43394069E+01;
    COFD[521] = 3.17651319E+00;
    COFD[522] = -2.02028974E-01;
    COFD[523] = 8.94232502E-03;
    COFD[524] = -2.06463744E+01;
    COFD[525] = 5.41688482E+00;
    COFD[526] = -4.73387188E-01;
    COFD[527] = 1.99280175E-02;
    COFD[528] = -1.55666415E+01;
    COFD[529] = 3.48070094E+00;
    COFD[530] = -2.40859499E-01;
    COFD[531] = 1.05972514E-02;
    COFD[532] = -1.63301444E+01;
    COFD[533] = 3.82388595E+00;
    COFD[534] = -2.84480724E-01;
    COFD[535] = 1.24506311E-02;
    COFD[536] = -1.63301444E+01;
    COFD[537] = 3.82388595E+00;
    COFD[538] = -2.84480724E-01;
    COFD[539] = 1.24506311E-02;
    COFD[540] = -1.63542394E+01;
    COFD[541] = 3.82388595E+00;
    COFD[542] = -2.84480724E-01;
    COFD[543] = 1.24506311E-02;
    COFD[544] = -1.62775714E+01;
    COFD[545] = 3.79163564E+00;
    COFD[546] = -2.80257365E-01;
    COFD[547] = 1.22656902E-02;
    COFD[548] = -1.52792891E+01;
    COFD[549] = 3.36790500E+00;
    COFD[550] = -2.26321740E-01;
    COFD[551] = 9.97135055E-03;
    COFD[552] = -1.84777607E+01;
    COFD[553] = 4.49330851E+00;
    COFD[554] = -3.68208715E-01;
    COFD[555] = 1.59565402E-02;
    COFD[556] = -2.08277598E+01;
    COFD[557] = 5.35267674E+00;
    COFD[558] = -4.69010505E-01;
    COFD[559] = 1.98979152E-02;
    COFD[560] = -2.08367725E+01;
    COFD[561] = 5.35267674E+00;
    COFD[562] = -4.69010505E-01;
    COFD[563] = 1.98979152E-02;
    COFD[564] = -2.02637994E+01;
    COFD[565] = 5.14984081E+00;
    COFD[566] = -4.46093018E-01;
    COFD[567] = 1.90396647E-02;
    COFD[568] = -1.89616623E+01;
    COFD[569] = 4.68595732E+00;
    COFD[570] = -3.91842840E-01;
    COFD[571] = 1.69262542E-02;
    COFD[572] = -1.86409139E+01;
    COFD[573] = 4.53572533E+00;
    COFD[574] = -3.73386925E-01;
    COFD[575] = 1.61678881E-02;
    COFD[576] = -1.86499071E+01;
    COFD[577] = 4.53572533E+00;
    COFD[578] = -3.73386925E-01;
    COFD[579] = 1.61678881E-02;
    COFD[580] = -1.52486273E+01;
    COFD[581] = 3.35922578E+00;
    COFD[582] = -2.25181399E-01;
    COFD[583] = 9.92132878E-03;
    COFD[584] = -1.62465583E+01;
    COFD[585] = 3.72612300E+00;
    COFD[586] = -2.71663673E-01;
    COFD[587] = 1.18889643E-02;
    COFD[588] = -1.25098960E+01;
    COFD[589] = 2.77873601E+00;
    COFD[590] = -1.50637360E-01;
    COFD[591] = 6.72684281E-03;
    COFD[592] = -1.57972369E+01;
    COFD[593] = 4.22225052E+00;
    COFD[594] = -3.35156428E-01;
    COFD[595] = 1.46104855E-02;
    COFD[596] = -1.50584249E+01;
    COFD[597] = 3.47945612E+00;
    COFD[598] = -2.40703722E-01;
    COFD[599] = 1.05907441E-02;
    COFD[600] = -1.63254691E+01;
    COFD[601] = 3.82388595E+00;
    COFD[602] = -2.84480724E-01;
    COFD[603] = 1.24506311E-02;
    COFD[604] = -1.50724636E+01;
    COFD[605] = 3.47945612E+00;
    COFD[606] = -2.40703722E-01;
    COFD[607] = 1.05907441E-02;
    COFD[608] = -2.12639214E+01;
    COFD[609] = 5.61184117E+00;
    COFD[610] = -4.90532156E-01;
    COFD[611] = 2.03507922E-02;
    COFD[612] = -1.63301444E+01;
    COFD[613] = 3.82388595E+00;
    COFD[614] = -2.84480724E-01;
    COFD[615] = 1.24506311E-02;
    COFD[616] = -1.73027557E+01;
    COFD[617] = 4.21416723E+00;
    COFD[618] = -3.34163932E-01;
    COFD[619] = 1.45697432E-02;
    COFD[620] = -1.73027557E+01;
    COFD[621] = 4.21416723E+00;
    COFD[622] = -3.34163932E-01;
    COFD[623] = 1.45697432E-02;
    COFD[624] = -1.73198034E+01;
    COFD[625] = 4.21416723E+00;
    COFD[626] = -3.34163932E-01;
    COFD[627] = 1.45697432E-02;
    COFD[628] = -1.72556729E+01;
    COFD[629] = 4.19029808E+00;
    COFD[630] = -3.31177076E-01;
    COFD[631] = 1.44446234E-02;
    COFD[632] = -1.59634533E+01;
    COFD[633] = 3.67388294E+00;
    COFD[634] = -2.64990709E-01;
    COFD[635] = 1.16042706E-02;
    COFD[636] = -1.93015555E+01;
    COFD[637] = 4.85015581E+00;
    COFD[638] = -4.10945109E-01;
    COFD[639] = 1.76651398E-02;
    COFD[640] = -2.14160703E+01;
    COFD[641] = 5.56531152E+00;
    COFD[642] = -4.88789821E-01;
    COFD[643] = 2.04437116E-02;
    COFD[644] = -2.14215700E+01;
    COFD[645] = 5.56531152E+00;
    COFD[646] = -4.88789821E-01;
    COFD[647] = 2.04437116E-02;
    COFD[648] = -2.09376196E+01;
    COFD[649] = 5.40870099E+00;
    COFD[650] = -4.73017610E-01;
    COFD[651] = 1.99399066E-02;
    COFD[652] = -1.98418115E+01;
    COFD[653] = 5.04367502E+00;
    COFD[654] = -4.34153325E-01;
    COFD[655] = 1.85956055E-02;
    COFD[656] = -1.95263312E+01;
    COFD[657] = 4.90255048E+00;
    COFD[658] = -4.17368501E-01;
    COFD[659] = 1.79287358E-02;
    COFD[660] = -1.95318173E+01;
    COFD[661] = 4.90255048E+00;
    COFD[662] = -4.17368501E-01;
    COFD[663] = 1.79287358E-02;
    COFD[664] = -1.59404882E+01;
    COFD[665] = 3.66853818E+00;
    COFD[666] = -2.64346221E-01;
    COFD[667] = 1.15784613E-02;
    COFD[668] = -1.71942502E+01;
    COFD[669] = 4.14993355E+00;
    COFD[670] = -3.26168062E-01;
    COFD[671] = 1.42364115E-02;
    COFD[672] = -1.25098960E+01;
    COFD[673] = 2.77873601E+00;
    COFD[674] = -1.50637360E-01;
    COFD[675] = 6.72684281E-03;
    COFD[676] = -1.57972369E+01;
    COFD[677] = 4.22225052E+00;
    COFD[678] = -3.35156428E-01;
    COFD[679] = 1.46104855E-02;
    COFD[680] = -1.50584249E+01;
    COFD[681] = 3.47945612E+00;
    COFD[682] = -2.40703722E-01;
    COFD[683] = 1.05907441E-02;
    COFD[684] = -1.63254691E+01;
    COFD[685] = 3.82388595E+00;
    COFD[686] = -2.84480724E-01;
    COFD[687] = 1.24506311E-02;
    COFD[688] = -1.50724636E+01;
    COFD[689] = 3.47945612E+00;
    COFD[690] = -2.40703722E-01;
    COFD[691] = 1.05907441E-02;
    COFD[692] = -2.12639214E+01;
    COFD[693] = 5.61184117E+00;
    COFD[694] = -4.90532156E-01;
    COFD[695] = 2.03507922E-02;
    COFD[696] = -1.63301444E+01;
    COFD[697] = 3.82388595E+00;
    COFD[698] = -2.84480724E-01;
    COFD[699] = 1.24506311E-02;
    COFD[700] = -1.73027557E+01;
    COFD[701] = 4.21416723E+00;
    COFD[702] = -3.34163932E-01;
    COFD[703] = 1.45697432E-02;
    COFD[704] = -1.73027557E+01;
    COFD[705] = 4.21416723E+00;
    COFD[706] = -3.34163932E-01;
    COFD[707] = 1.45697432E-02;
    COFD[708] = -1.73198034E+01;
    COFD[709] = 4.21416723E+00;
    COFD[710] = -3.34163932E-01;
    COFD[711] = 1.45697432E-02;
    COFD[712] = -1.72556729E+01;
    COFD[713] = 4.19029808E+00;
    COFD[714] = -3.31177076E-01;
    COFD[715] = 1.44446234E-02;
    COFD[716] = -1.59634533E+01;
    COFD[717] = 3.67388294E+00;
    COFD[718] = -2.64990709E-01;
    COFD[719] = 1.16042706E-02;
    COFD[720] = -1.93015555E+01;
    COFD[721] = 4.85015581E+00;
    COFD[722] = -4.10945109E-01;
    COFD[723] = 1.76651398E-02;
    COFD[724] = -2.14160703E+01;
    COFD[725] = 5.56531152E+00;
    COFD[726] = -4.88789821E-01;
    COFD[727] = 2.04437116E-02;
    COFD[728] = -2.14215700E+01;
    COFD[729] = 5.56531152E+00;
    COFD[730] = -4.88789821E-01;
    COFD[731] = 2.04437116E-02;
    COFD[732] = -2.09376196E+01;
    COFD[733] = 5.40870099E+00;
    COFD[734] = -4.73017610E-01;
    COFD[735] = 1.99399066E-02;
    COFD[736] = -1.98418115E+01;
    COFD[737] = 5.04367502E+00;
    COFD[738] = -4.34153325E-01;
    COFD[739] = 1.85956055E-02;
    COFD[740] = -1.95263312E+01;
    COFD[741] = 4.90255048E+00;
    COFD[742] = -4.17368501E-01;
    COFD[743] = 1.79287358E-02;
    COFD[744] = -1.95318173E+01;
    COFD[745] = 4.90255048E+00;
    COFD[746] = -4.17368501E-01;
    COFD[747] = 1.79287358E-02;
    COFD[748] = -1.59404882E+01;
    COFD[749] = 3.66853818E+00;
    COFD[750] = -2.64346221E-01;
    COFD[751] = 1.15784613E-02;
    COFD[752] = -1.71942502E+01;
    COFD[753] = 4.14993355E+00;
    COFD[754] = -3.26168062E-01;
    COFD[755] = 1.42364115E-02;
    COFD[756] = -1.25141260E+01;
    COFD[757] = 2.77873601E+00;
    COFD[758] = -1.50637360E-01;
    COFD[759] = 6.72684281E-03;
    COFD[760] = -1.57994893E+01;
    COFD[761] = 4.22225052E+00;
    COFD[762] = -3.35156428E-01;
    COFD[763] = 1.46104855E-02;
    COFD[764] = -1.50766130E+01;
    COFD[765] = 3.47945612E+00;
    COFD[766] = -2.40703722E-01;
    COFD[767] = 1.05907441E-02;
    COFD[768] = -1.63493345E+01;
    COFD[769] = 3.82388595E+00;
    COFD[770] = -2.84480724E-01;
    COFD[771] = 1.24506311E-02;
    COFD[772] = -1.50911794E+01;
    COFD[773] = 3.47945612E+00;
    COFD[774] = -2.40703722E-01;
    COFD[775] = 1.05907441E-02;
    COFD[776] = -2.12831323E+01;
    COFD[777] = 5.61184117E+00;
    COFD[778] = -4.90532156E-01;
    COFD[779] = 2.03507922E-02;
    COFD[780] = -1.63542394E+01;
    COFD[781] = 3.82388595E+00;
    COFD[782] = -2.84480724E-01;
    COFD[783] = 1.24506311E-02;
    COFD[784] = -1.73198034E+01;
    COFD[785] = 4.21416723E+00;
    COFD[786] = -3.34163932E-01;
    COFD[787] = 1.45697432E-02;
    COFD[788] = -1.73198034E+01;
    COFD[789] = 4.21416723E+00;
    COFD[790] = -3.34163932E-01;
    COFD[791] = 1.45697432E-02;
    COFD[792] = -1.73374529E+01;
    COFD[793] = 4.21416723E+00;
    COFD[794] = -3.34163932E-01;
    COFD[795] = 1.45697432E-02;
    COFD[796] = -1.72738845E+01;
    COFD[797] = 4.19029808E+00;
    COFD[798] = -3.31177076E-01;
    COFD[799] = 1.44446234E-02;
    COFD[800] = -1.59863030E+01;
    COFD[801] = 3.67388294E+00;
    COFD[802] = -2.64990709E-01;
    COFD[803] = 1.16042706E-02;
    COFD[804] = -1.93276434E+01;
    COFD[805] = 4.85015581E+00;
    COFD[806] = -4.10945109E-01;
    COFD[807] = 1.76651398E-02;
    COFD[808] = -2.14391943E+01;
    COFD[809] = 5.56531152E+00;
    COFD[810] = -4.88789821E-01;
    COFD[811] = 2.04437116E-02;
    COFD[812] = -2.14449559E+01;
    COFD[813] = 5.56531152E+00;
    COFD[814] = -4.88789821E-01;
    COFD[815] = 2.04437116E-02;
    COFD[816] = -2.09612557E+01;
    COFD[817] = 5.40870099E+00;
    COFD[818] = -4.73017610E-01;
    COFD[819] = 1.99399066E-02;
    COFD[820] = -1.98646734E+01;
    COFD[821] = 5.04367502E+00;
    COFD[822] = -4.34153325E-01;
    COFD[823] = 1.85956055E-02;
    COFD[824] = -1.95494668E+01;
    COFD[825] = 4.90255048E+00;
    COFD[826] = -4.17368501E-01;
    COFD[827] = 1.79287358E-02;
    COFD[828] = -1.95552142E+01;
    COFD[829] = 4.90255048E+00;
    COFD[830] = -4.17368501E-01;
    COFD[831] = 1.79287358E-02;
    COFD[832] = -1.59633387E+01;
    COFD[833] = 3.66853818E+00;
    COFD[834] = -2.64346221E-01;
    COFD[835] = 1.15784613E-02;
    COFD[836] = -1.72196961E+01;
    COFD[837] = 4.14993355E+00;
    COFD[838] = -3.26168062E-01;
    COFD[839] = 1.42364115E-02;
    COFD[840] = -1.24693568E+01;
    COFD[841] = 2.76686648E+00;
    COFD[842] = -1.49120141E-01;
    COFD[843] = 6.66220432E-03;
    COFD[844] = -1.57199037E+01;
    COFD[845] = 4.19936335E+00;
    COFD[846] = -3.32311009E-01;
    COFD[847] = 1.44921003E-02;
    COFD[848] = -1.50270339E+01;
    COFD[849] = 3.46140064E+00;
    COFD[850] = -2.38440092E-01;
    COFD[851] = 1.04960087E-02;
    COFD[852] = -1.62724462E+01;
    COFD[853] = 3.79163564E+00;
    COFD[854] = -2.80257365E-01;
    COFD[855] = 1.22656902E-02;
    COFD[856] = -1.50420953E+01;
    COFD[857] = 3.46140064E+00;
    COFD[858] = -2.38440092E-01;
    COFD[859] = 1.04960087E-02;
    COFD[860] = -2.14087397E+01;
    COFD[861] = 5.57282008E+00;
    COFD[862] = -4.76690890E-01;
    COFD[863] = 1.94000719E-02;
    COFD[864] = -1.62775714E+01;
    COFD[865] = 3.79163564E+00;
    COFD[866] = -2.80257365E-01;
    COFD[867] = 1.22656902E-02;
    COFD[868] = -1.72556729E+01;
    COFD[869] = 4.19029808E+00;
    COFD[870] = -3.31177076E-01;
    COFD[871] = 1.44446234E-02;
    COFD[872] = -1.72556729E+01;
    COFD[873] = 4.19029808E+00;
    COFD[874] = -3.31177076E-01;
    COFD[875] = 1.44446234E-02;
    COFD[876] = -1.72738845E+01;
    COFD[877] = 4.19029808E+00;
    COFD[878] = -3.31177076E-01;
    COFD[879] = 1.44446234E-02;
    COFD[880] = -1.72167708E+01;
    COFD[881] = 4.16886779E+00;
    COFD[882] = -3.28518156E-01;
    COFD[883] = 1.43341626E-02;
    COFD[884] = -1.59525102E+01;
    COFD[885] = 3.66023858E+00;
    COFD[886] = -2.63401043E-01;
    COFD[887] = 1.15432000E-02;
    COFD[888] = -1.92867554E+01;
    COFD[889] = 4.83375900E+00;
    COFD[890] = -4.09146560E-01;
    COFD[891] = 1.76006599E-02;
    COFD[892] = -2.14022336E+01;
    COFD[893] = 5.55346617E+00;
    COFD[894] = -4.87783156E-01;
    COFD[895] = 2.04210886E-02;
    COFD[896] = -2.14082453E+01;
    COFD[897] = 5.55346617E+00;
    COFD[898] = -4.87783156E-01;
    COFD[899] = 2.04210886E-02;
    COFD[900] = -2.11381508E+01;
    COFD[901] = 5.45574440E+00;
    COFD[902] = -4.77436155E-01;
    COFD[903] = 2.00644596E-02;
    COFD[904] = -1.98075055E+01;
    COFD[905] = 5.02169524E+00;
    COFD[906] = -4.31582804E-01;
    COFD[907] = 1.84953568E-02;
    COFD[908] = -1.94763688E+01;
    COFD[909] = 4.87333294E+00;
    COFD[910] = -4.13769241E-01;
    COFD[911] = 1.77802244E-02;
    COFD[912] = -1.94823660E+01;
    COFD[913] = 4.87333294E+00;
    COFD[914] = -4.13769241E-01;
    COFD[915] = 1.77802244E-02;
    COFD[916] = -1.59327297E+01;
    COFD[917] = 3.65620899E+00;
    COFD[918] = -2.62933804E-01;
    COFD[919] = 1.15253223E-02;
    COFD[920] = -1.71754154E+01;
    COFD[921] = 4.13131681E+00;
    COFD[922] = -3.23897559E-01;
    COFD[923] = 1.41438222E-02;
    COFD[924] = -1.17159737E+01;
    COFD[925] = 2.48123210E+00;
    COFD[926] = -1.11322604E-01;
    COFD[927] = 4.99282389E-03;
    COFD[928] = -1.43151174E+01;
    COFD[929] = 3.68038508E+00;
    COFD[930] = -2.65779346E-01;
    COFD[931] = 1.16360771E-02;
    COFD[932] = -1.40999008E+01;
    COFD[933] = 3.08120012E+00;
    COFD[934] = -1.89629903E-01;
    COFD[935] = 8.40361952E-03;
    COFD[936] = -1.52721107E+01;
    COFD[937] = 3.36790500E+00;
    COFD[938] = -2.26321740E-01;
    COFD[939] = 9.97135055E-03;
    COFD[940] = -1.41191261E+01;
    COFD[941] = 3.08120012E+00;
    COFD[942] = -1.89629903E-01;
    COFD[943] = 8.40361952E-03;
    COFD[944] = -2.11388331E+01;
    COFD[945] = 5.55529675E+00;
    COFD[946] = -4.87942518E-01;
    COFD[947] = 2.04249054E-02;
    COFD[948] = -1.52792891E+01;
    COFD[949] = 3.36790500E+00;
    COFD[950] = -2.26321740E-01;
    COFD[951] = 9.97135055E-03;
    COFD[952] = -1.59634533E+01;
    COFD[953] = 3.67388294E+00;
    COFD[954] = -2.64990709E-01;
    COFD[955] = 1.16042706E-02;
    COFD[956] = -1.59634533E+01;
    COFD[957] = 3.67388294E+00;
    COFD[958] = -2.64990709E-01;
    COFD[959] = 1.16042706E-02;
    COFD[960] = -1.59863030E+01;
    COFD[961] = 3.67388294E+00;
    COFD[962] = -2.64990709E-01;
    COFD[963] = 1.16042706E-02;
    COFD[964] = -1.59525102E+01;
    COFD[965] = 3.66023858E+00;
    COFD[966] = -2.63401043E-01;
    COFD[967] = 1.15432000E-02;
    COFD[968] = -1.50233475E+01;
    COFD[969] = 3.26660767E+00;
    COFD[970] = -2.13287177E-01;
    COFD[971] = 9.41137857E-03;
    COFD[972] = -1.81735763E+01;
    COFD[973] = 4.38391495E+00;
    COFD[974] = -3.54941287E-01;
    COFD[975] = 1.54195107E-02;
    COFD[976] = -2.05045578E+01;
    COFD[977] = 5.23843909E+00;
    COFD[978] = -4.55815614E-01;
    COFD[979] = 1.93898040E-02;
    COFD[980] = -2.05128705E+01;
    COFD[981] = 5.23843909E+00;
    COFD[982] = -4.55815614E-01;
    COFD[983] = 1.93898040E-02;
    COFD[984] = -2.02642227E+01;
    COFD[985] = 5.14499740E+00;
    COFD[986] = -4.45694430E-01;
    COFD[987] = 1.90318646E-02;
    COFD[988] = -1.86157761E+01;
    COFD[989] = 4.55689508E+00;
    COFD[990] = -3.75937921E-01;
    COFD[991] = 1.62703488E-02;
    COFD[992] = -1.83455435E+01;
    COFD[993] = 4.42828044E+00;
    COFD[994] = -3.60417833E-01;
    COFD[995] = 1.56455103E-02;
    COFD[996] = -1.83538377E+01;
    COFD[997] = 4.42828044E+00;
    COFD[998] = -3.60417833E-01;
    COFD[999] = 1.56455103E-02;
    COFD[1000] = -1.50031687E+01;
    COFD[1001] = 3.26223357E+00;
    COFD[1002] = -2.12746642E-01;
    COFD[1003] = 9.38912883E-03;
    COFD[1004] = -1.60074211E+01;
    COFD[1005] = 3.63723937E+00;
    COFD[1006] = -2.60754222E-01;
    COFD[1007] = 1.14428814E-02;
    COFD[1008] = -1.37794315E+01;
    COFD[1009] = 3.23973858E+00;
    COFD[1010] = -2.09989036E-01;
    COFD[1011] = 9.27667906E-03;
    COFD[1012] = -1.76147026E+01;
    COFD[1013] = 4.86049500E+00;
    COFD[1014] = -4.12200578E-01;
    COFD[1015] = 1.77160971E-02;
    COFD[1016] = -1.70534856E+01;
    COFD[1017] = 4.14240922E+00;
    COFD[1018] = -3.25239774E-01;
    COFD[1019] = 1.41980687E-02;
    COFD[1020] = -1.84688406E+01;
    COFD[1021] = 4.49330851E+00;
    COFD[1022] = -3.68208715E-01;
    COFD[1023] = 1.59565402E-02;
    COFD[1024] = -1.70757047E+01;
    COFD[1025] = 4.14240922E+00;
    COFD[1026] = -3.25239774E-01;
    COFD[1027] = 1.41980687E-02;
    COFD[1028] = -2.07653719E+01;
    COFD[1029] = 5.01092022E+00;
    COFD[1030] = -3.77985635E-01;
    COFD[1031] = 1.40968645E-02;
    COFD[1032] = -1.84777607E+01;
    COFD[1033] = 4.49330851E+00;
    COFD[1034] = -3.68208715E-01;
    COFD[1035] = 1.59565402E-02;
    COFD[1036] = -1.93015555E+01;
    COFD[1037] = 4.85015581E+00;
    COFD[1038] = -4.10945109E-01;
    COFD[1039] = 1.76651398E-02;
    COFD[1040] = -1.93015555E+01;
    COFD[1041] = 4.85015581E+00;
    COFD[1042] = -4.10945109E-01;
    COFD[1043] = 1.76651398E-02;
    COFD[1044] = -1.93276434E+01;
    COFD[1045] = 4.85015581E+00;
    COFD[1046] = -4.10945109E-01;
    COFD[1047] = 1.76651398E-02;
    COFD[1048] = -1.92867554E+01;
    COFD[1049] = 4.83375900E+00;
    COFD[1050] = -4.09146560E-01;
    COFD[1051] = 1.76006599E-02;
    COFD[1052] = -1.81735763E+01;
    COFD[1053] = 4.38391495E+00;
    COFD[1054] = -3.54941287E-01;
    COFD[1055] = 1.54195107E-02;
    COFD[1056] = -2.13425698E+01;
    COFD[1057] = 5.40460130E+00;
    COFD[1058] = -4.72718910E-01;
    COFD[1059] = 1.99362717E-02;
    COFD[1060] = -2.19215555E+01;
    COFD[1061] = 5.45216133E+00;
    COFD[1062] = -4.52916925E-01;
    COFD[1063] = 1.80456400E-02;
    COFD[1064] = -2.19317743E+01;
    COFD[1065] = 5.45216133E+00;
    COFD[1066] = -4.52916925E-01;
    COFD[1067] = 1.80456400E-02;
    COFD[1068] = -2.20421041E+01;
    COFD[1069] = 5.52708332E+00;
    COFD[1070] = -4.68000808E-01;
    COFD[1071] = 1.89131908E-02;
    COFD[1072] = -2.16802612E+01;
    COFD[1073] = 5.52918296E+00;
    COFD[1074] = -4.85360709E-01;
    COFD[1075] = 2.03448006E-02;
    COFD[1076] = -2.14224484E+01;
    COFD[1077] = 5.41729961E+00;
    COFD[1078] = -4.73400281E-01;
    COFD[1079] = 1.99269567E-02;
    COFD[1080] = -2.14326461E+01;
    COFD[1081] = 5.41729961E+00;
    COFD[1082] = -4.73400281E-01;
    COFD[1083] = 1.99269567E-02;
    COFD[1084] = -1.81432461E+01;
    COFD[1085] = 4.37565431E+00;
    COFD[1086] = -3.53906025E-01;
    COFD[1087] = 1.53760786E-02;
    COFD[1088] = -1.93483692E+01;
    COFD[1089] = 4.79506290E+00;
    COFD[1090] = -4.04621659E-01;
    COFD[1091] = 1.74244230E-02;
    COFD[1092] = -1.60517370E+01;
    COFD[1093] = 4.11188603E+00;
    COFD[1094] = -3.21540884E-01;
    COFD[1095] = 1.40482564E-02;
    COFD[1096] = -1.97544450E+01;
    COFD[1097] = 5.56931926E+00;
    COFD[1098] = -4.89105511E-01;
    COFD[1099] = 2.04493129E-02;
    COFD[1100] = -1.94313116E+01;
    COFD[1101] = 5.02567894E+00;
    COFD[1102] = -4.32045169E-01;
    COFD[1103] = 1.85132214E-02;
    COFD[1104] = -2.08204449E+01;
    COFD[1105] = 5.35267674E+00;
    COFD[1106] = -4.69010505E-01;
    COFD[1107] = 1.98979152E-02;
    COFD[1108] = -1.94507876E+01;
    COFD[1109] = 5.02567894E+00;
    COFD[1110] = -4.32045169E-01;
    COFD[1111] = 1.85132214E-02;
    COFD[1112] = -1.77498543E+01;
    COFD[1113] = 3.57475686E+00;
    COFD[1114] = -1.56396297E-01;
    COFD[1115] = 3.12157721E-03;
    COFD[1116] = -2.08277598E+01;
    COFD[1117] = 5.35267674E+00;
    COFD[1118] = -4.69010505E-01;
    COFD[1119] = 1.98979152E-02;
    COFD[1120] = -2.14160703E+01;
    COFD[1121] = 5.56531152E+00;
    COFD[1122] = -4.88789821E-01;
    COFD[1123] = 2.04437116E-02;
    COFD[1124] = -2.14160703E+01;
    COFD[1125] = 5.56531152E+00;
    COFD[1126] = -4.88789821E-01;
    COFD[1127] = 2.04437116E-02;
    COFD[1128] = -2.14391943E+01;
    COFD[1129] = 5.56531152E+00;
    COFD[1130] = -4.88789821E-01;
    COFD[1131] = 2.04437116E-02;
    COFD[1132] = -2.14022336E+01;
    COFD[1133] = 5.55346617E+00;
    COFD[1134] = -4.87783156E-01;
    COFD[1135] = 2.04210886E-02;
    COFD[1136] = -2.05045578E+01;
    COFD[1137] = 5.23843909E+00;
    COFD[1138] = -4.55815614E-01;
    COFD[1139] = 1.93898040E-02;
    COFD[1140] = -2.19215555E+01;
    COFD[1141] = 5.45216133E+00;
    COFD[1142] = -4.52916925E-01;
    COFD[1143] = 1.80456400E-02;
    COFD[1144] = -1.90328712E+01;
    COFD[1145] = 3.99221757E+00;
    COFD[1146] = -2.19854880E-01;
    COFD[1147] = 6.22736279E-03;
    COFD[1148] = -1.90413348E+01;
    COFD[1149] = 3.99221757E+00;
    COFD[1150] = -2.19854880E-01;
    COFD[1151] = 6.22736279E-03;
    COFD[1152] = -2.01801667E+01;
    COFD[1153] = 4.53183330E+00;
    COFD[1154] = -3.02186760E-01;
    COFD[1155] = 1.02756490E-02;
    COFD[1156] = -2.16296373E+01;
    COFD[1157] = 5.29019717E+00;
    COFD[1158] = -4.24502606E-01;
    COFD[1159] = 1.65197343E-02;
    COFD[1160] = -2.19229190E+01;
    COFD[1161] = 5.41841631E+00;
    COFD[1162] = -4.46818971E-01;
    COFD[1163] = 1.77127652E-02;
    COFD[1164] = -2.19313638E+01;
    COFD[1165] = 5.41841631E+00;
    COFD[1166] = -4.46818971E-01;
    COFD[1167] = 1.77127652E-02;
    COFD[1168] = -2.04750581E+01;
    COFD[1169] = 5.23112374E+00;
    COFD[1170] = -4.54967682E-01;
    COFD[1171] = 1.93570423E-02;
    COFD[1172] = -2.14255087E+01;
    COFD[1173] = 5.52240865E+00;
    COFD[1174] = -4.84699537E-01;
    COFD[1175] = 2.03247833E-02;
    COFD[1176] = -1.60528285E+01;
    COFD[1177] = 4.11188603E+00;
    COFD[1178] = -3.21540884E-01;
    COFD[1179] = 1.40482564E-02;
    COFD[1180] = -1.97550088E+01;
    COFD[1181] = 5.56931926E+00;
    COFD[1182] = -4.89105511E-01;
    COFD[1183] = 2.04493129E-02;
    COFD[1184] = -1.94373127E+01;
    COFD[1185] = 5.02567894E+00;
    COFD[1186] = -4.32045169E-01;
    COFD[1187] = 1.85132214E-02;
    COFD[1188] = -2.08293255E+01;
    COFD[1189] = 5.35267674E+00;
    COFD[1190] = -4.69010505E-01;
    COFD[1191] = 1.98979152E-02;
    COFD[1192] = -1.94570287E+01;
    COFD[1193] = 5.02567894E+00;
    COFD[1194] = -4.32045169E-01;
    COFD[1195] = 1.85132214E-02;
    COFD[1196] = -1.77563250E+01;
    COFD[1197] = 3.57475686E+00;
    COFD[1198] = -1.56396297E-01;
    COFD[1199] = 3.12157721E-03;
    COFD[1200] = -2.08367725E+01;
    COFD[1201] = 5.35267674E+00;
    COFD[1202] = -4.69010505E-01;
    COFD[1203] = 1.98979152E-02;
    COFD[1204] = -2.14215700E+01;
    COFD[1205] = 5.56531152E+00;
    COFD[1206] = -4.88789821E-01;
    COFD[1207] = 2.04437116E-02;
    COFD[1208] = -2.14215700E+01;
    COFD[1209] = 5.56531152E+00;
    COFD[1210] = -4.88789821E-01;
    COFD[1211] = 2.04437116E-02;
    COFD[1212] = -2.14449559E+01;
    COFD[1213] = 5.56531152E+00;
    COFD[1214] = -4.88789821E-01;
    COFD[1215] = 2.04437116E-02;
    COFD[1216] = -2.14082453E+01;
    COFD[1217] = 5.55346617E+00;
    COFD[1218] = -4.87783156E-01;
    COFD[1219] = 2.04210886E-02;
    COFD[1220] = -2.05128705E+01;
    COFD[1221] = 5.23843909E+00;
    COFD[1222] = -4.55815614E-01;
    COFD[1223] = 1.93898040E-02;
    COFD[1224] = -2.19317743E+01;
    COFD[1225] = 5.45216133E+00;
    COFD[1226] = -4.52916925E-01;
    COFD[1227] = 1.80456400E-02;
    COFD[1228] = -1.90413348E+01;
    COFD[1229] = 3.99221757E+00;
    COFD[1230] = -2.19854880E-01;
    COFD[1231] = 6.22736279E-03;
    COFD[1232] = -1.90499441E+01;
    COFD[1233] = 3.99221757E+00;
    COFD[1234] = -2.19854880E-01;
    COFD[1235] = 6.22736279E-03;
    COFD[1236] = -2.01889168E+01;
    COFD[1237] = 4.53183330E+00;
    COFD[1238] = -3.02186760E-01;
    COFD[1239] = 1.02756490E-02;
    COFD[1240] = -2.16379567E+01;
    COFD[1241] = 5.29019717E+00;
    COFD[1242] = -4.24502606E-01;
    COFD[1243] = 1.65197343E-02;
    COFD[1244] = -2.19313890E+01;
    COFD[1245] = 5.41841631E+00;
    COFD[1246] = -4.46818971E-01;
    COFD[1247] = 1.77127652E-02;
    COFD[1248] = -2.19399793E+01;
    COFD[1249] = 5.41841631E+00;
    COFD[1250] = -4.46818971E-01;
    COFD[1251] = 1.77127652E-02;
    COFD[1252] = -2.04833713E+01;
    COFD[1253] = 5.23112374E+00;
    COFD[1254] = -4.54967682E-01;
    COFD[1255] = 1.93570423E-02;
    COFD[1256] = -2.14353267E+01;
    COFD[1257] = 5.52240865E+00;
    COFD[1258] = -4.84699537E-01;
    COFD[1259] = 2.03247833E-02;
    COFD[1260] = -1.58456300E+01;
    COFD[1261] = 4.02074783E+00;
    COFD[1262] = -3.10018522E-01;
    COFD[1263] = 1.35599552E-02;
    COFD[1264] = -1.92718582E+01;
    COFD[1265] = 5.41172124E+00;
    COFD[1266] = -4.73213887E-01;
    COFD[1267] = 1.99405473E-02;
    COFD[1268] = -1.88179418E+01;
    COFD[1269] = 4.79683898E+00;
    COFD[1270] = -4.04829719E-01;
    COFD[1271] = 1.74325475E-02;
    COFD[1272] = -2.04928958E+01;
    COFD[1273] = 5.22397933E+00;
    COFD[1274] = -4.54138171E-01;
    COFD[1275] = 1.93249285E-02;
    COFD[1276] = -1.88378874E+01;
    COFD[1277] = 4.79683898E+00;
    COFD[1278] = -4.04829719E-01;
    COFD[1279] = 1.74325475E-02;
    COFD[1280] = -1.65295288E+01;
    COFD[1281] = 2.97569206E+00;
    COFD[1282] = -6.75652842E-02;
    COFD[1283] = -1.08648422E-03;
    COFD[1284] = -2.02637994E+01;
    COFD[1285] = 5.14984081E+00;
    COFD[1286] = -4.46093018E-01;
    COFD[1287] = 1.90396647E-02;
    COFD[1288] = -2.09376196E+01;
    COFD[1289] = 5.40870099E+00;
    COFD[1290] = -4.73017610E-01;
    COFD[1291] = 1.99399066E-02;
    COFD[1292] = -2.09376196E+01;
    COFD[1293] = 5.40870099E+00;
    COFD[1294] = -4.73017610E-01;
    COFD[1295] = 1.99399066E-02;
    COFD[1296] = -2.09612557E+01;
    COFD[1297] = 5.40870099E+00;
    COFD[1298] = -4.73017610E-01;
    COFD[1299] = 1.99399066E-02;
    COFD[1300] = -2.11381508E+01;
    COFD[1301] = 5.45574440E+00;
    COFD[1302] = -4.77436155E-01;
    COFD[1303] = 2.00644596E-02;
    COFD[1304] = -2.02642227E+01;
    COFD[1305] = 5.14499740E+00;
    COFD[1306] = -4.45694430E-01;
    COFD[1307] = 1.90318646E-02;
    COFD[1308] = -2.20421041E+01;
    COFD[1309] = 5.52708332E+00;
    COFD[1310] = -4.68000808E-01;
    COFD[1311] = 1.89131908E-02;
    COFD[1312] = -2.01801667E+01;
    COFD[1313] = 4.53183330E+00;
    COFD[1314] = -3.02186760E-01;
    COFD[1315] = 1.02756490E-02;
    COFD[1316] = -2.01889168E+01;
    COFD[1317] = 4.53183330E+00;
    COFD[1318] = -3.02186760E-01;
    COFD[1319] = 1.02756490E-02;
    COFD[1320] = -1.95877017E+01;
    COFD[1321] = 4.27643051E+00;
    COFD[1322] = -2.68040901E-01;
    COFD[1323] = 8.77650113E-03;
    COFD[1324] = -2.19670848E+01;
    COFD[1325] = 5.48847873E+00;
    COFD[1326] = -4.59558930E-01;
    COFD[1327] = 1.84107961E-02;
    COFD[1328] = -2.21070030E+01;
    COFD[1329] = 5.55072945E+00;
    COFD[1330] = -4.72525345E-01;
    COFD[1331] = 1.91674202E-02;
    COFD[1332] = -2.21157340E+01;
    COFD[1333] = 5.55072945E+00;
    COFD[1334] = -4.72525345E-01;
    COFD[1335] = 1.91674202E-02;
    COFD[1336] = -2.02268902E+01;
    COFD[1337] = 5.13632093E+00;
    COFD[1338] = -4.44839124E-01;
    COFD[1339] = 1.90058354E-02;
    COFD[1340] = -2.10026861E+01;
    COFD[1341] = 5.38326647E+00;
    COFD[1342] = -4.71201048E-01;
    COFD[1343] = 1.99207516E-02;
    COFD[1344] = -1.42229194E+01;
    COFD[1345] = 3.38669384E+00;
    COFD[1346] = -2.28784122E-01;
    COFD[1347] = 1.00790953E-02;
    COFD[1348] = -1.82251914E+01;
    COFD[1349] = 5.05237312E+00;
    COFD[1350] = -4.35182396E-01;
    COFD[1351] = 1.86363074E-02;
    COFD[1352] = -1.74792112E+01;
    COFD[1353] = 4.29676909E+00;
    COFD[1354] = -3.44085306E-01;
    COFD[1355] = 1.49671135E-02;
    COFD[1356] = -1.89544778E+01;
    COFD[1357] = 4.68595732E+00;
    COFD[1358] = -3.91842840E-01;
    COFD[1359] = 1.69262542E-02;
    COFD[1360] = -1.74984476E+01;
    COFD[1361] = 4.29676909E+00;
    COFD[1362] = -3.44085306E-01;
    COFD[1363] = 1.49671135E-02;
    COFD[1364] = -2.08812333E+01;
    COFD[1365] = 5.08859217E+00;
    COFD[1366] = -3.90525428E-01;
    COFD[1367] = 1.47376395E-02;
    COFD[1368] = -1.89616623E+01;
    COFD[1369] = 4.68595732E+00;
    COFD[1370] = -3.91842840E-01;
    COFD[1371] = 1.69262542E-02;
    COFD[1372] = -1.98418115E+01;
    COFD[1373] = 5.04367502E+00;
    COFD[1374] = -4.34153325E-01;
    COFD[1375] = 1.85956055E-02;
    COFD[1376] = -1.98418115E+01;
    COFD[1377] = 5.04367502E+00;
    COFD[1378] = -4.34153325E-01;
    COFD[1379] = 1.85956055E-02;
    COFD[1380] = -1.98646734E+01;
    COFD[1381] = 5.04367502E+00;
    COFD[1382] = -4.34153325E-01;
    COFD[1383] = 1.85956055E-02;
    COFD[1384] = -1.98075055E+01;
    COFD[1385] = 5.02169524E+00;
    COFD[1386] = -4.31582804E-01;
    COFD[1387] = 1.84953568E-02;
    COFD[1388] = -1.86157761E+01;
    COFD[1389] = 4.55689508E+00;
    COFD[1390] = -3.75937921E-01;
    COFD[1391] = 1.62703488E-02;
    COFD[1392] = -2.16802612E+01;
    COFD[1393] = 5.52918296E+00;
    COFD[1394] = -4.85360709E-01;
    COFD[1395] = 2.03448006E-02;
    COFD[1396] = -2.16296373E+01;
    COFD[1397] = 5.29019717E+00;
    COFD[1398] = -4.24502606E-01;
    COFD[1399] = 1.65197343E-02;
    COFD[1400] = -2.16379567E+01;
    COFD[1401] = 5.29019717E+00;
    COFD[1402] = -4.24502606E-01;
    COFD[1403] = 1.65197343E-02;
    COFD[1404] = -2.19670848E+01;
    COFD[1405] = 5.48847873E+00;
    COFD[1406] = -4.59558930E-01;
    COFD[1407] = 1.84107961E-02;
    COFD[1408] = -2.19327397E+01;
    COFD[1409] = 5.60638188E+00;
    COFD[1410] = -4.91272522E-01;
    COFD[1411] = 2.04396264E-02;
    COFD[1412] = -2.18190539E+01;
    COFD[1413] = 5.55753905E+00;
    COFD[1414] = -4.88136714E-01;
    COFD[1415] = 2.04294957E-02;
    COFD[1416] = -2.18273547E+01;
    COFD[1417] = 5.55753905E+00;
    COFD[1418] = -4.88136714E-01;
    COFD[1419] = 2.04294957E-02;
    COFD[1420] = -1.85864144E+01;
    COFD[1421] = 4.54915847E+00;
    COFD[1422] = -3.75000738E-01;
    COFD[1423] = 1.62324821E-02;
    COFD[1424] = -1.98040322E+01;
    COFD[1425] = 4.97569695E+00;
    COFD[1426] = -4.26123307E-01;
    COFD[1427] = 1.82788664E-02;
    COFD[1428] = -1.39913897E+01;
    COFD[1429] = 3.26384506E+00;
    COFD[1430] = -2.12947087E-01;
    COFD[1431] = 9.39743888E-03;
    COFD[1432] = -1.79339327E+01;
    COFD[1433] = 4.91373893E+00;
    COFD[1434] = -4.18747629E-01;
    COFD[1435] = 1.79856610E-02;
    COFD[1436] = -1.72496634E+01;
    COFD[1437] = 4.17889917E+00;
    COFD[1438] = -3.29752510E-01;
    COFD[1439] = 1.43850275E-02;
    COFD[1440] = -1.86335932E+01;
    COFD[1441] = 4.53572533E+00;
    COFD[1442] = -3.73386925E-01;
    COFD[1443] = 1.61678881E-02;
    COFD[1444] = -1.72691500E+01;
    COFD[1445] = 4.17889917E+00;
    COFD[1446] = -3.29752510E-01;
    COFD[1447] = 1.43850275E-02;
    COFD[1448] = -2.12597312E+01;
    COFD[1449] = 5.24930667E+00;
    COFD[1450] = -4.17435088E-01;
    COFD[1451] = 1.61434424E-02;
    COFD[1452] = -1.86409139E+01;
    COFD[1453] = 4.53572533E+00;
    COFD[1454] = -3.73386925E-01;
    COFD[1455] = 1.61678881E-02;
    COFD[1456] = -1.95263312E+01;
    COFD[1457] = 4.90255048E+00;
    COFD[1458] = -4.17368501E-01;
    COFD[1459] = 1.79287358E-02;
    COFD[1460] = -1.95263312E+01;
    COFD[1461] = 4.90255048E+00;
    COFD[1462] = -4.17368501E-01;
    COFD[1463] = 1.79287358E-02;
    COFD[1464] = -1.95494668E+01;
    COFD[1465] = 4.90255048E+00;
    COFD[1466] = -4.17368501E-01;
    COFD[1467] = 1.79287358E-02;
    COFD[1468] = -1.94763688E+01;
    COFD[1469] = 4.87333294E+00;
    COFD[1470] = -4.13769241E-01;
    COFD[1471] = 1.77802244E-02;
    COFD[1472] = -1.83455435E+01;
    COFD[1473] = 4.42828044E+00;
    COFD[1474] = -3.60417833E-01;
    COFD[1475] = 1.56455103E-02;
    COFD[1476] = -2.14224484E+01;
    COFD[1477] = 5.41729961E+00;
    COFD[1478] = -4.73400281E-01;
    COFD[1479] = 1.99269567E-02;
    COFD[1480] = -2.19229190E+01;
    COFD[1481] = 5.41841631E+00;
    COFD[1482] = -4.46818971E-01;
    COFD[1483] = 1.77127652E-02;
    COFD[1484] = -2.19313890E+01;
    COFD[1485] = 5.41841631E+00;
    COFD[1486] = -4.46818971E-01;
    COFD[1487] = 1.77127652E-02;
    COFD[1488] = -2.21070030E+01;
    COFD[1489] = 5.55072945E+00;
    COFD[1490] = -4.72525345E-01;
    COFD[1491] = 1.91674202E-02;
    COFD[1492] = -2.18190539E+01;
    COFD[1493] = 5.55753905E+00;
    COFD[1494] = -4.88136714E-01;
    COFD[1495] = 2.04294957E-02;
    COFD[1496] = -2.15575659E+01;
    COFD[1497] = 5.44803850E+00;
    COFD[1498] = -4.76610560E-01;
    COFD[1499] = 2.00355294E-02;
    COFD[1500] = -2.15660171E+01;
    COFD[1501] = 5.44803850E+00;
    COFD[1502] = -4.76610560E-01;
    COFD[1503] = 2.00355294E-02;
    COFD[1504] = -1.83166353E+01;
    COFD[1505] = 4.42045763E+00;
    COFD[1506] = -3.59451578E-01;
    COFD[1507] = 1.56056164E-02;
    COFD[1508] = -1.94928815E+01;
    COFD[1509] = 4.83189721E+00;
    COFD[1510] = -4.08932249E-01;
    COFD[1511] = 1.75924650E-02;
    COFD[1512] = -1.39924781E+01;
    COFD[1513] = 3.26384506E+00;
    COFD[1514] = -2.12947087E-01;
    COFD[1515] = 9.39743888E-03;
    COFD[1516] = -1.79344949E+01;
    COFD[1517] = 4.91373893E+00;
    COFD[1518] = -4.18747629E-01;
    COFD[1519] = 1.79856610E-02;
    COFD[1520] = -1.72556499E+01;
    COFD[1521] = 4.17889917E+00;
    COFD[1522] = -3.29752510E-01;
    COFD[1523] = 1.43850275E-02;
    COFD[1524] = -1.86424545E+01;
    COFD[1525] = 4.53572533E+00;
    COFD[1526] = -3.73386925E-01;
    COFD[1527] = 1.61678881E-02;
    COFD[1528] = -1.72753760E+01;
    COFD[1529] = 4.17889917E+00;
    COFD[1530] = -3.29752510E-01;
    COFD[1531] = 1.43850275E-02;
    COFD[1532] = -2.12661865E+01;
    COFD[1533] = 5.24930667E+00;
    COFD[1534] = -4.17435088E-01;
    COFD[1535] = 1.61434424E-02;
    COFD[1536] = -1.86499071E+01;
    COFD[1537] = 4.53572533E+00;
    COFD[1538] = -3.73386925E-01;
    COFD[1539] = 1.61678881E-02;
    COFD[1540] = -1.95318173E+01;
    COFD[1541] = 4.90255048E+00;
    COFD[1542] = -4.17368501E-01;
    COFD[1543] = 1.79287358E-02;
    COFD[1544] = -1.95318173E+01;
    COFD[1545] = 4.90255048E+00;
    COFD[1546] = -4.17368501E-01;
    COFD[1547] = 1.79287358E-02;
    COFD[1548] = -1.95552142E+01;
    COFD[1549] = 4.90255048E+00;
    COFD[1550] = -4.17368501E-01;
    COFD[1551] = 1.79287358E-02;
    COFD[1552] = -1.94823660E+01;
    COFD[1553] = 4.87333294E+00;
    COFD[1554] = -4.13769241E-01;
    COFD[1555] = 1.77802244E-02;
    COFD[1556] = -1.83538377E+01;
    COFD[1557] = 4.42828044E+00;
    COFD[1558] = -3.60417833E-01;
    COFD[1559] = 1.56455103E-02;
    COFD[1560] = -2.14326461E+01;
    COFD[1561] = 5.41729961E+00;
    COFD[1562] = -4.73400281E-01;
    COFD[1563] = 1.99269567E-02;
    COFD[1564] = -2.19313638E+01;
    COFD[1565] = 5.41841631E+00;
    COFD[1566] = -4.46818971E-01;
    COFD[1567] = 1.77127652E-02;
    COFD[1568] = -2.19399793E+01;
    COFD[1569] = 5.41841631E+00;
    COFD[1570] = -4.46818971E-01;
    COFD[1571] = 1.77127652E-02;
    COFD[1572] = -2.21157340E+01;
    COFD[1573] = 5.55072945E+00;
    COFD[1574] = -4.72525345E-01;
    COFD[1575] = 1.91674202E-02;
    COFD[1576] = -2.18273547E+01;
    COFD[1577] = 5.55753905E+00;
    COFD[1578] = -4.88136714E-01;
    COFD[1579] = 2.04294957E-02;
    COFD[1580] = -2.15660171E+01;
    COFD[1581] = 5.44803850E+00;
    COFD[1582] = -4.76610560E-01;
    COFD[1583] = 2.00355294E-02;
    COFD[1584] = -2.15746136E+01;
    COFD[1585] = 5.44803850E+00;
    COFD[1586] = -4.76610560E-01;
    COFD[1587] = 2.00355294E-02;
    COFD[1588] = -1.83249299E+01;
    COFD[1589] = 4.42045763E+00;
    COFD[1590] = -3.59451578E-01;
    COFD[1591] = 1.56056164E-02;
    COFD[1592] = -1.95026789E+01;
    COFD[1593] = 4.83189721E+00;
    COFD[1594] = -4.08932249E-01;
    COFD[1595] = 1.75924650E-02;
    COFD[1596] = -1.16906297E+01;
    COFD[1597] = 2.47469981E+00;
    COFD[1598] = -1.10436257E-01;
    COFD[1599] = 4.95273813E-03;
    COFD[1600] = -1.42894441E+01;
    COFD[1601] = 3.67490723E+00;
    COFD[1602] = -2.65114792E-01;
    COFD[1603] = 1.16092671E-02;
    COFD[1604] = -1.40756935E+01;
    COFD[1605] = 3.07549274E+00;
    COFD[1606] = -1.88889344E-01;
    COFD[1607] = 8.37152866E-03;
    COFD[1608] = -1.52414485E+01;
    COFD[1609] = 3.35922578E+00;
    COFD[1610] = -2.25181399E-01;
    COFD[1611] = 9.92132878E-03;
    COFD[1612] = -1.40949196E+01;
    COFD[1613] = 3.07549274E+00;
    COFD[1614] = -1.88889344E-01;
    COFD[1615] = 8.37152866E-03;
    COFD[1616] = -2.10643259E+01;
    COFD[1617] = 5.53614847E+00;
    COFD[1618] = -4.86046736E-01;
    COFD[1619] = 2.03659188E-02;
    COFD[1620] = -1.52486273E+01;
    COFD[1621] = 3.35922578E+00;
    COFD[1622] = -2.25181399E-01;
    COFD[1623] = 9.92132878E-03;
    COFD[1624] = -1.59404882E+01;
    COFD[1625] = 3.66853818E+00;
    COFD[1626] = -2.64346221E-01;
    COFD[1627] = 1.15784613E-02;
    COFD[1628] = -1.59404882E+01;
    COFD[1629] = 3.66853818E+00;
    COFD[1630] = -2.64346221E-01;
    COFD[1631] = 1.15784613E-02;
    COFD[1632] = -1.59633387E+01;
    COFD[1633] = 3.66853818E+00;
    COFD[1634] = -2.64346221E-01;
    COFD[1635] = 1.15784613E-02;
    COFD[1636] = -1.59327297E+01;
    COFD[1637] = 3.65620899E+00;
    COFD[1638] = -2.62933804E-01;
    COFD[1639] = 1.15253223E-02;
    COFD[1640] = -1.50031687E+01;
    COFD[1641] = 3.26223357E+00;
    COFD[1642] = -2.12746642E-01;
    COFD[1643] = 9.38912883E-03;
    COFD[1644] = -1.81432461E+01;
    COFD[1645] = 4.37565431E+00;
    COFD[1646] = -3.53906025E-01;
    COFD[1647] = 1.53760786E-02;
    COFD[1648] = -2.04750581E+01;
    COFD[1649] = 5.23112374E+00;
    COFD[1650] = -4.54967682E-01;
    COFD[1651] = 1.93570423E-02;
    COFD[1652] = -2.04833713E+01;
    COFD[1653] = 5.23112374E+00;
    COFD[1654] = -4.54967682E-01;
    COFD[1655] = 1.93570423E-02;
    COFD[1656] = -2.02268902E+01;
    COFD[1657] = 5.13632093E+00;
    COFD[1658] = -4.44839124E-01;
    COFD[1659] = 1.90058354E-02;
    COFD[1660] = -1.85864144E+01;
    COFD[1661] = 4.54915847E+00;
    COFD[1662] = -3.75000738E-01;
    COFD[1663] = 1.62324821E-02;
    COFD[1664] = -1.83166353E+01;
    COFD[1665] = 4.42045763E+00;
    COFD[1666] = -3.59451578E-01;
    COFD[1667] = 1.56056164E-02;
    COFD[1668] = -1.83249299E+01;
    COFD[1669] = 4.42045763E+00;
    COFD[1670] = -3.59451578E-01;
    COFD[1671] = 1.56056164E-02;
    COFD[1672] = -1.49828430E+01;
    COFD[1673] = 3.25781069E+00;
    COFD[1674] = -2.12199367E-01;
    COFD[1675] = 9.36657283E-03;
    COFD[1676] = -1.59877782E+01;
    COFD[1677] = 3.63340763E+00;
    COFD[1678] = -2.60307961E-01;
    COFD[1679] = 1.14256954E-02;
    COFD[1680] = -1.23130152E+01;
    COFD[1681] = 2.74418790E+00;
    COFD[1682] = -1.46230156E-01;
    COFD[1683] = 6.53948886E-03;
    COFD[1684] = -1.54738604E+01;
    COFD[1685] = 4.15765300E+00;
    COFD[1686] = -3.27126237E-01;
    COFD[1687] = 1.42762611E-02;
    COFD[1688] = -1.49610527E+01;
    COFD[1689] = 3.41988961E+00;
    COFD[1690] = -2.33128386E-01;
    COFD[1691] = 1.02689994E-02;
    COFD[1692] = -1.62380075E+01;
    COFD[1693] = 3.72612300E+00;
    COFD[1694] = -2.71663673E-01;
    COFD[1695] = 1.18889643E-02;
    COFD[1696] = -1.49826725E+01;
    COFD[1697] = 3.41988961E+00;
    COFD[1698] = -2.33128386E-01;
    COFD[1699] = 1.02689994E-02;
    COFD[1700] = -2.12755888E+01;
    COFD[1701] = 5.60381989E+00;
    COFD[1702] = -4.91225459E-01;
    COFD[1703] = 2.04487844E-02;
    COFD[1704] = -1.62465583E+01;
    COFD[1705] = 3.72612300E+00;
    COFD[1706] = -2.71663673E-01;
    COFD[1707] = 1.18889643E-02;
    COFD[1708] = -1.71942502E+01;
    COFD[1709] = 4.14993355E+00;
    COFD[1710] = -3.26168062E-01;
    COFD[1711] = 1.42364115E-02;
    COFD[1712] = -1.71942502E+01;
    COFD[1713] = 4.14993355E+00;
    COFD[1714] = -3.26168062E-01;
    COFD[1715] = 1.42364115E-02;
    COFD[1716] = -1.72196961E+01;
    COFD[1717] = 4.14993355E+00;
    COFD[1718] = -3.26168062E-01;
    COFD[1719] = 1.42364115E-02;
    COFD[1720] = -1.71754154E+01;
    COFD[1721] = 4.13131681E+00;
    COFD[1722] = -3.23897559E-01;
    COFD[1723] = 1.41438222E-02;
    COFD[1724] = -1.60074211E+01;
    COFD[1725] = 3.63723937E+00;
    COFD[1726] = -2.60754222E-01;
    COFD[1727] = 1.14428814E-02;
    COFD[1728] = -1.93483692E+01;
    COFD[1729] = 4.79506290E+00;
    COFD[1730] = -4.04621659E-01;
    COFD[1731] = 1.74244230E-02;
    COFD[1732] = -2.14255087E+01;
    COFD[1733] = 5.52240865E+00;
    COFD[1734] = -4.84699537E-01;
    COFD[1735] = 2.03247833E-02;
    COFD[1736] = -2.14353267E+01;
    COFD[1737] = 5.52240865E+00;
    COFD[1738] = -4.84699537E-01;
    COFD[1739] = 2.03247833E-02;
    COFD[1740] = -2.10026861E+01;
    COFD[1741] = 5.38326647E+00;
    COFD[1742] = -4.71201048E-01;
    COFD[1743] = 1.99207516E-02;
    COFD[1744] = -1.98040322E+01;
    COFD[1745] = 4.97569695E+00;
    COFD[1746] = -4.26123307E-01;
    COFD[1747] = 1.82788664E-02;
    COFD[1748] = -1.94928815E+01;
    COFD[1749] = 4.83189721E+00;
    COFD[1750] = -4.08932249E-01;
    COFD[1751] = 1.75924650E-02;
    COFD[1752] = -1.95026789E+01;
    COFD[1753] = 4.83189721E+00;
    COFD[1754] = -4.08932249E-01;
    COFD[1755] = 1.75924650E-02;
    COFD[1756] = -1.59877782E+01;
    COFD[1757] = 3.63340763E+00;
    COFD[1758] = -2.60307961E-01;
    COFD[1759] = 1.14256954E-02;
    COFD[1760] = -1.72273911E+01;
    COFD[1761] = 4.09361913E+00;
    COFD[1762] = -3.19258125E-01;
    COFD[1763] = 1.39526981E-02;
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
    COFTD[4] = -1.44152190E-01;
    COFTD[5] = -7.99993584E-05;
    COFTD[6] = 4.89707442E-08;
    COFTD[7] = -9.14277269E-12;
    COFTD[8] = 4.06682492E-01;
    COFTD[9] = 3.84705248E-05;
    COFTD[10] = -2.54846868E-08;
    COFTD[11] = 5.86302354E-12;
    COFTD[12] = 4.26579943E-01;
    COFTD[13] = 1.20407274E-04;
    COFTD[14] = -7.67298757E-08;
    COFTD[15] = 1.52090336E-11;
    COFTD[16] = 4.12895615E-01;
    COFTD[17] = 3.90582612E-05;
    COFTD[18] = -2.58740310E-08;
    COFTD[19] = 5.95259633E-12;
    COFTD[20] = 2.27469146E-02;
    COFTD[21] = 6.73078907E-04;
    COFTD[22] = -3.40935843E-07;
    COFTD[23] = 5.48499211E-11;
    COFTD[24] = 4.28230888E-01;
    COFTD[25] = 1.20873273E-04;
    COFTD[26] = -7.70268349E-08;
    COFTD[27] = 1.52678954E-11;
    COFTD[28] = 3.24747031E-01;
    COFTD[29] = 1.77798548E-04;
    COFTD[30] = -1.08934732E-07;
    COFTD[31] = 2.03595881E-11;
    COFTD[32] = 3.24747031E-01;
    COFTD[33] = 1.77798548E-04;
    COFTD[34] = -1.08934732E-07;
    COFTD[35] = 2.03595881E-11;
    COFTD[36] = 3.31191185E-01;
    COFTD[37] = 1.81326714E-04;
    COFTD[38] = -1.11096391E-07;
    COFTD[39] = 2.07635959E-11;
    COFTD[40] = 3.39557243E-01;
    COFTD[41] = 1.79335036E-04;
    COFTD[42] = -1.10135705E-07;
    COFTD[43] = 2.06427239E-11;
    COFTD[44] = 4.30605547E-01;
    COFTD[45] = 9.35961902E-05;
    COFTD[46] = -6.03983623E-08;
    COFTD[47] = 1.23115170E-11;
    COFTD[48] = 2.93191523E-01;
    COFTD[49] = 4.01430006E-04;
    COFTD[50] = -2.30705763E-07;
    COFTD[51] = 4.05176586E-11;
    COFTD[52] = 1.22119780E-01;
    COFTD[53] = 6.18373616E-04;
    COFTD[54] = -3.28422593E-07;
    COFTD[55] = 5.44603522E-11;
    COFTD[56] = 1.22693382E-01;
    COFTD[57] = 6.21278143E-04;
    COFTD[58] = -3.29965208E-07;
    COFTD[59] = 5.47161548E-11;
    COFTD[60] = 1.40314191E-01;
    COFTD[61] = 6.01266129E-04;
    COFTD[62] = -3.21915137E-07;
    COFTD[63] = 5.36679068E-11;
    COFTD[64] = 2.49017478E-01;
    COFTD[65] = 4.29036573E-04;
    COFTD[66] = -2.42668617E-07;
    COFTD[67] = 4.20801371E-11;
    COFTD[68] = 2.72759599E-01;
    COFTD[69] = 3.94402719E-04;
    COFTD[70] = -2.25800520E-07;
    COFTD[71] = 3.95325634E-11;
    COFTD[72] = 2.74036956E-01;
    COFTD[73] = 3.96249742E-04;
    COFTD[74] = -2.26857964E-07;
    COFTD[75] = 3.97176979E-11;
    COFTD[76] = 4.31331269E-01;
    COFTD[77] = 9.20536800E-05;
    COFTD[78] = -5.94509616E-08;
    COFTD[79] = 1.21437993E-11;
    COFTD[80] = 4.01012808E-01;
    COFTD[81] = 1.97252826E-04;
    COFTD[82] = -1.21698146E-07;
    COFTD[83] = 2.29408847E-11;
    COFTD[84] = 1.44152190E-01;
    COFTD[85] = 7.99993584E-05;
    COFTD[86] = -4.89707442E-08;
    COFTD[87] = 9.14277269E-12;
    COFTD[88] = 0.00000000E+00;
    COFTD[89] = 0.00000000E+00;
    COFTD[90] = 0.00000000E+00;
    COFTD[91] = 0.00000000E+00;
    COFTD[92] = 2.35283119E-01;
    COFTD[93] = 4.65670599E-04;
    COFTD[94] = -2.60939824E-07;
    COFTD[95] = 4.49271822E-11;
    COFTD[96] = 1.79840299E-01;
    COFTD[97] = 6.01722902E-04;
    COFTD[98] = -3.26433894E-07;
    COFTD[99] = 5.49112302E-11;
    COFTD[100] = 2.37053352E-01;
    COFTD[101] = 4.69174231E-04;
    COFTD[102] = -2.62903094E-07;
    COFTD[103] = 4.52652072E-11;
    COFTD[104] = -1.74352698E-01;
    COFTD[105] = 8.62246873E-04;
    COFTD[106] = -3.79545489E-07;
    COFTD[107] = 5.60262093E-11;
    COFTD[108] = 1.80186965E-01;
    COFTD[109] = 6.02882805E-04;
    COFTD[110] = -3.27063140E-07;
    COFTD[111] = 5.50170790E-11;
    COFTD[112] = 9.90752318E-02;
    COFTD[113] = 6.44201384E-04;
    COFTD[114] = -3.38485953E-07;
    COFTD[115] = 5.57356746E-11;
    COFTD[116] = 9.90752318E-02;
    COFTD[117] = 6.44201384E-04;
    COFTD[118] = -3.38485953E-07;
    COFTD[119] = 5.57356746E-11;
    COFTD[120] = 1.00039110E-01;
    COFTD[121] = 6.50468660E-04;
    COFTD[122] = -3.41778999E-07;
    COFTD[123] = 5.62779132E-11;
    COFTD[124] = 1.05124122E-01;
    COFTD[125] = 6.50665957E-04;
    COFTD[126] = -3.42564538E-07;
    COFTD[127] = 5.64804120E-11;
    COFTD[128] = 2.00119897E-01;
    COFTD[129] = 5.64793704E-04;
    COFTD[130] = -3.09445484E-07;
    COFTD[131] = 5.24139335E-11;
    COFTD[132] = -2.00309448E-02;
    COFTD[133] = 8.50440115E-04;
    COFTD[134] = -4.21064468E-07;
    COFTD[135] = 6.67959710E-11;
    COFTD[136] = -1.60981264E-01;
    COFTD[137] = 9.03807572E-04;
    COFTD[138] = -4.06927941E-07;
    COFTD[139] = 6.09202254E-11;
    COFTD[140] = -1.61357564E-01;
    COFTD[141] = 9.05920260E-04;
    COFTD[142] = -4.07879153E-07;
    COFTD[143] = 6.10626290E-11;
    COFTD[144] = -1.31244519E-01;
    COFTD[145] = 9.03901384E-04;
    COFTD[146] = -4.17831507E-07;
    COFTD[147] = 6.35725667E-11;
    COFTD[148] = -5.08744745E-02;
    COFTD[149] = 8.54342586E-04;
    COFTD[150] = -4.15926453E-07;
    COFTD[151] = 6.53063261E-11;
    COFTD[152] = -2.71690558E-02;
    COFTD[153] = 8.37233133E-04;
    COFTD[154] = -4.12887636E-07;
    COFTD[155] = 6.53405197E-11;
    COFTD[156] = -2.72323768E-02;
    COFTD[157] = 8.39184413E-04;
    COFTD[158] = -4.13849924E-07;
    COFTD[159] = 6.54928043E-11;
    COFTD[160] = 2.01521643E-01;
    COFTD[161] = 5.62744089E-04;
    COFTD[162] = -3.08519239E-07;
    COFTD[163] = 5.22805986E-11;
    COFTD[164] = 1.22193921E-01;
    COFTD[165] = 6.90321128E-04;
    COFTD[166] = -3.64844875E-07;
    COFTD[167] = 6.03054876E-11;
}

#ifndef AMREX_USE_GPU
/* Initializes parameter database */
void CKINIT()
{

    // (0):  O + H + M <=> OH + M
    // (0):  O + H + M <=> OH + M
    fwd_A[8]     = 5e+17;
    fwd_beta[8]  = -1;
    fwd_Ea[8]    = 0;
    prefactor_units[8]  = 1.0000000000000002e-12;
    activation_units[8] = 0.50321666580471969;
    phase_units[8]      = pow(10,-12.000000);
    is_PD[8] = 0;
    nTB[8] = 7;
    TB[8] = (amrex::Real *) malloc(7 * sizeof(amrex::Real));
    TBid[8] = (int *) malloc(7 * sizeof(int));
    TBid[8][0] = 0; TB[8][0] = 2; // H2
    TBid[8][1] = 5; TB[8][1] = 6; // H2O
    TBid[8][2] = 10; TB[8][2] = 2; // CH4
    TBid[8][3] = 11; TB[8][3] = 1.5; // CO
    TBid[8][4] = 12; TB[8][4] = 2; // CO2
    TBid[8][5] = 18; TB[8][5] = 3; // C2H6
    TBid[8][6] = 20; TB[8][6] = 0.69999999999999996; // AR

    // (1):  O + H2 <=> H + OH
    // (1):  O + H2 <=> H + OH
    fwd_A[14]     = 50000;
    fwd_beta[14]  = 2.6699999999999999;
    fwd_Ea[14]    = 6290;
    prefactor_units[14]  = 1.0000000000000002e-06;
    activation_units[14] = 0.50321666580471969;
    phase_units[14]      = pow(10,-12.000000);
    is_PD[14] = 0;
    nTB[14] = 0;

    // (2):  O + HO2 <=> OH + O2
    // (2):  O + HO2 <=> OH + O2
    fwd_A[15]     = 20000000000000;
    fwd_beta[15]  = 0;
    fwd_Ea[15]    = 0;
    prefactor_units[15]  = 1.0000000000000002e-06;
    activation_units[15] = 0.50321666580471969;
    phase_units[15]      = pow(10,-12.000000);
    is_PD[15] = 0;
    nTB[15] = 0;

    // (3):  O + CH2 <=> H + HCO
    // (3):  O + CH2 <=> H + HCO
    fwd_A[16]     = 80000000000000;
    fwd_beta[16]  = 0;
    fwd_Ea[16]    = 0;
    prefactor_units[16]  = 1.0000000000000002e-06;
    activation_units[16] = 0.50321666580471969;
    phase_units[16]      = pow(10,-12.000000);
    is_PD[16] = 0;
    nTB[16] = 0;

    // (4):  O + CH2(S) <=> H + HCO
    // (4):  O + CH2(S) <=> H + HCO
    fwd_A[17]     = 15000000000000;
    fwd_beta[17]  = 0;
    fwd_Ea[17]    = 0;
    prefactor_units[17]  = 1.0000000000000002e-06;
    activation_units[17] = 0.50321666580471969;
    phase_units[17]      = pow(10,-12.000000);
    is_PD[17] = 0;
    nTB[17] = 0;

    // (5):  O + CH3 <=> H + CH2O
    // (5):  O + CH3 <=> H + CH2O
    fwd_A[18]     = 84300000000000;
    fwd_beta[18]  = 0;
    fwd_Ea[18]    = 0;
    prefactor_units[18]  = 1.0000000000000002e-06;
    activation_units[18] = 0.50321666580471969;
    phase_units[18]      = pow(10,-12.000000);
    is_PD[18] = 0;
    nTB[18] = 0;

    // (6):  O + CH4 <=> OH + CH3
    // (6):  O + CH4 <=> OH + CH3
    fwd_A[19]     = 1020000000;
    fwd_beta[19]  = 1.5;
    fwd_Ea[19]    = 8600;
    prefactor_units[19]  = 1.0000000000000002e-06;
    activation_units[19] = 0.50321666580471969;
    phase_units[19]      = pow(10,-12.000000);
    is_PD[19] = 0;
    nTB[19] = 0;

    // (7):  O + CO + M <=> CO2 + M
    // (7):  O + CO + M <=> CO2 + M
    fwd_A[9]     = 602000000000000;
    fwd_beta[9]  = 0;
    fwd_Ea[9]    = 3000;
    prefactor_units[9]  = 1.0000000000000002e-12;
    activation_units[9] = 0.50321666580471969;
    phase_units[9]      = pow(10,-12.000000);
    is_PD[9] = 0;
    nTB[9] = 8;
    TB[9] = (amrex::Real *) malloc(8 * sizeof(amrex::Real));
    TBid[9] = (int *) malloc(8 * sizeof(int));
    TBid[9][0] = 0; TB[9][0] = 2; // H2
    TBid[9][1] = 3; TB[9][1] = 6; // O2
    TBid[9][2] = 5; TB[9][2] = 6; // H2O
    TBid[9][3] = 10; TB[9][3] = 2; // CH4
    TBid[9][4] = 11; TB[9][4] = 1.5; // CO
    TBid[9][5] = 12; TB[9][5] = 3.5; // CO2
    TBid[9][6] = 18; TB[9][6] = 3; // C2H6
    TBid[9][7] = 20; TB[9][7] = 0.5; // AR

    // (8):  O + HCO <=> OH + CO
    // (8):  O + HCO <=> OH + CO
    fwd_A[20]     = 30000000000000;
    fwd_beta[20]  = 0;
    fwd_Ea[20]    = 0;
    prefactor_units[20]  = 1.0000000000000002e-06;
    activation_units[20] = 0.50321666580471969;
    phase_units[20]      = pow(10,-12.000000);
    is_PD[20] = 0;
    nTB[20] = 0;

    // (9):  O + HCO <=> H + CO2
    // (9):  O + HCO <=> H + CO2
    fwd_A[21]     = 30000000000000;
    fwd_beta[21]  = 0;
    fwd_Ea[21]    = 0;
    prefactor_units[21]  = 1.0000000000000002e-06;
    activation_units[21] = 0.50321666580471969;
    phase_units[21]      = pow(10,-12.000000);
    is_PD[21] = 0;
    nTB[21] = 0;

    // (10):  O + CH2O <=> OH + HCO
    // (10):  O + CH2O <=> OH + HCO
    fwd_A[22]     = 39000000000000;
    fwd_beta[22]  = 0;
    fwd_Ea[22]    = 3540;
    prefactor_units[22]  = 1.0000000000000002e-06;
    activation_units[22] = 0.50321666580471969;
    phase_units[22]      = pow(10,-12.000000);
    is_PD[22] = 0;
    nTB[22] = 0;

    // (11):  O + C2H4 <=> CH3 + HCO
    // (11):  O + C2H4 <=> CH3 + HCO
    fwd_A[23]     = 19200000;
    fwd_beta[23]  = 1.8300000000000001;
    fwd_Ea[23]    = 220;
    prefactor_units[23]  = 1.0000000000000002e-06;
    activation_units[23] = 0.50321666580471969;
    phase_units[23]      = pow(10,-12.000000);
    is_PD[23] = 0;
    nTB[23] = 0;

    // (12):  O + C2H5 <=> CH3 + CH2O
    // (12):  O + C2H5 <=> CH3 + CH2O
    fwd_A[24]     = 132000000000000;
    fwd_beta[24]  = 0;
    fwd_Ea[24]    = 0;
    prefactor_units[24]  = 1.0000000000000002e-06;
    activation_units[24] = 0.50321666580471969;
    phase_units[24]      = pow(10,-12.000000);
    is_PD[24] = 0;
    nTB[24] = 0;

    // (13):  O + C2H6 <=> OH + C2H5
    // (13):  O + C2H6 <=> OH + C2H5
    fwd_A[25]     = 89800000;
    fwd_beta[25]  = 1.9199999999999999;
    fwd_Ea[25]    = 5690;
    prefactor_units[25]  = 1.0000000000000002e-06;
    activation_units[25] = 0.50321666580471969;
    phase_units[25]      = pow(10,-12.000000);
    is_PD[25] = 0;
    nTB[25] = 0;

    // (14):  O2 + CO <=> O + CO2
    // (14):  O2 + CO <=> O + CO2
    fwd_A[26]     = 2500000000000;
    fwd_beta[26]  = 0;
    fwd_Ea[26]    = 47800;
    prefactor_units[26]  = 1.0000000000000002e-06;
    activation_units[26] = 0.50321666580471969;
    phase_units[26]      = pow(10,-12.000000);
    is_PD[26] = 0;
    nTB[26] = 0;

    // (15):  O2 + CH2O <=> HO2 + HCO
    // (15):  O2 + CH2O <=> HO2 + HCO
    fwd_A[27]     = 100000000000000;
    fwd_beta[27]  = 0;
    fwd_Ea[27]    = 40000;
    prefactor_units[27]  = 1.0000000000000002e-06;
    activation_units[27] = 0.50321666580471969;
    phase_units[27]      = pow(10,-12.000000);
    is_PD[27] = 0;
    nTB[27] = 0;

    // (16):  H + O2 + M <=> HO2 + M
    // (16):  H + O2 + M <=> HO2 + M
    fwd_A[10]     = 2.8e+18;
    fwd_beta[10]  = -0.85999999999999999;
    fwd_Ea[10]    = 0;
    prefactor_units[10]  = 1.0000000000000002e-12;
    activation_units[10] = 0.50321666580471969;
    phase_units[10]      = pow(10,-12.000000);
    is_PD[10] = 0;
    nTB[10] = 7;
    TB[10] = (amrex::Real *) malloc(7 * sizeof(amrex::Real));
    TBid[10] = (int *) malloc(7 * sizeof(int));
    TBid[10][0] = 3; TB[10][0] = 0; // O2
    TBid[10][1] = 5; TB[10][1] = 0; // H2O
    TBid[10][2] = 11; TB[10][2] = 0.75; // CO
    TBid[10][3] = 12; TB[10][3] = 1.5; // CO2
    TBid[10][4] = 18; TB[10][4] = 1.5; // C2H6
    TBid[10][5] = 19; TB[10][5] = 0; // N2
    TBid[10][6] = 20; TB[10][6] = 0; // AR

    // (17):  H + 2.000000 O2 <=> HO2 + O2
    // (17):  H + 2.000000 O2 <=> HO2 + O2
    fwd_A[28]     = 3e+20;
    fwd_beta[28]  = -1.72;
    fwd_Ea[28]    = 0;
    prefactor_units[28]  = 1.0000000000000002e-12;
    activation_units[28] = 0.50321666580471969;
    phase_units[28]      = pow(10,-18.000000);
    is_PD[28] = 0;
    nTB[28] = 0;

    // (18):  H + O2 + H2O <=> HO2 + H2O
    // (18):  H + O2 + H2O <=> HO2 + H2O
    fwd_A[29]     = 9.38e+18;
    fwd_beta[29]  = -0.76000000000000001;
    fwd_Ea[29]    = 0;
    prefactor_units[29]  = 1.0000000000000002e-12;
    activation_units[29] = 0.50321666580471969;
    phase_units[29]      = pow(10,-18.000000);
    is_PD[29] = 0;
    nTB[29] = 0;

    // (19):  H + O2 + N2 <=> HO2 + N2
    // (19):  H + O2 + N2 <=> HO2 + N2
    fwd_A[30]     = 3.75e+20;
    fwd_beta[30]  = -1.72;
    fwd_Ea[30]    = 0;
    prefactor_units[30]  = 1.0000000000000002e-12;
    activation_units[30] = 0.50321666580471969;
    phase_units[30]      = pow(10,-18.000000);
    is_PD[30] = 0;
    nTB[30] = 0;

    // (20):  H + O2 + AR <=> HO2 + AR
    // (20):  H + O2 + AR <=> HO2 + AR
    fwd_A[31]     = 7e+17;
    fwd_beta[31]  = -0.80000000000000004;
    fwd_Ea[31]    = 0;
    prefactor_units[31]  = 1.0000000000000002e-12;
    activation_units[31] = 0.50321666580471969;
    phase_units[31]      = pow(10,-18.000000);
    is_PD[31] = 0;
    nTB[31] = 0;

    // (21):  H + O2 <=> O + OH
    // (21):  H + O2 <=> O + OH
    fwd_A[32]     = 83000000000000;
    fwd_beta[32]  = 0;
    fwd_Ea[32]    = 14413;
    prefactor_units[32]  = 1.0000000000000002e-06;
    activation_units[32] = 0.50321666580471969;
    phase_units[32]      = pow(10,-12.000000);
    is_PD[32] = 0;
    nTB[32] = 0;

    // (22):  2.000000 H + M <=> H2 + M
    // (22):  2.000000 H + M <=> H2 + M
    fwd_A[11]     = 1e+18;
    fwd_beta[11]  = -1;
    fwd_Ea[11]    = 0;
    prefactor_units[11]  = 1.0000000000000002e-12;
    activation_units[11] = 0.50321666580471969;
    phase_units[11]      = pow(10,-12.000000);
    is_PD[11] = 0;
    nTB[11] = 6;
    TB[11] = (amrex::Real *) malloc(6 * sizeof(amrex::Real));
    TBid[11] = (int *) malloc(6 * sizeof(int));
    TBid[11][0] = 0; TB[11][0] = 0; // H2
    TBid[11][1] = 5; TB[11][1] = 0; // H2O
    TBid[11][2] = 10; TB[11][2] = 2; // CH4
    TBid[11][3] = 12; TB[11][3] = 0; // CO2
    TBid[11][4] = 18; TB[11][4] = 3; // C2H6
    TBid[11][5] = 20; TB[11][5] = 0.63; // AR

    // (23):  2.000000 H + H2 <=> 2.000000 H2
    // (23):  2.000000 H + H2 <=> 2.000000 H2
    fwd_A[33]     = 90000000000000000;
    fwd_beta[33]  = -0.59999999999999998;
    fwd_Ea[33]    = 0;
    prefactor_units[33]  = 1.0000000000000002e-12;
    activation_units[33] = 0.50321666580471969;
    phase_units[33]      = pow(10,-18.000000);
    is_PD[33] = 0;
    nTB[33] = 0;

    // (24):  2.000000 H + H2O <=> H2 + H2O
    // (24):  2.000000 H + H2O <=> H2 + H2O
    fwd_A[34]     = 6e+19;
    fwd_beta[34]  = -1.25;
    fwd_Ea[34]    = 0;
    prefactor_units[34]  = 1.0000000000000002e-12;
    activation_units[34] = 0.50321666580471969;
    phase_units[34]      = pow(10,-18.000000);
    is_PD[34] = 0;
    nTB[34] = 0;

    // (25):  2.000000 H + CO2 <=> H2 + CO2
    // (25):  2.000000 H + CO2 <=> H2 + CO2
    fwd_A[35]     = 5.5e+20;
    fwd_beta[35]  = -2;
    fwd_Ea[35]    = 0;
    prefactor_units[35]  = 1.0000000000000002e-12;
    activation_units[35] = 0.50321666580471969;
    phase_units[35]      = pow(10,-18.000000);
    is_PD[35] = 0;
    nTB[35] = 0;

    // (26):  H + OH + M <=> H2O + M
    // (26):  H + OH + M <=> H2O + M
    fwd_A[12]     = 2.2e+22;
    fwd_beta[12]  = -2;
    fwd_Ea[12]    = 0;
    prefactor_units[12]  = 1.0000000000000002e-12;
    activation_units[12] = 0.50321666580471969;
    phase_units[12]      = pow(10,-12.000000);
    is_PD[12] = 0;
    nTB[12] = 5;
    TB[12] = (amrex::Real *) malloc(5 * sizeof(amrex::Real));
    TBid[12] = (int *) malloc(5 * sizeof(int));
    TBid[12][0] = 0; TB[12][0] = 0.72999999999999998; // H2
    TBid[12][1] = 5; TB[12][1] = 3.6499999999999999; // H2O
    TBid[12][2] = 10; TB[12][2] = 2; // CH4
    TBid[12][3] = 18; TB[12][3] = 3; // C2H6
    TBid[12][4] = 20; TB[12][4] = 0.38; // AR

    // (27):  H + HO2 <=> O2 + H2
    // (27):  H + HO2 <=> O2 + H2
    fwd_A[36]     = 28000000000000;
    fwd_beta[36]  = 0;
    fwd_Ea[36]    = 1068;
    prefactor_units[36]  = 1.0000000000000002e-06;
    activation_units[36] = 0.50321666580471969;
    phase_units[36]      = pow(10,-12.000000);
    is_PD[36] = 0;
    nTB[36] = 0;

    // (28):  H + HO2 <=> 2.000000 OH
    // (28):  H + HO2 <=> 2.000000 OH
    fwd_A[37]     = 134000000000000;
    fwd_beta[37]  = 0;
    fwd_Ea[37]    = 635;
    prefactor_units[37]  = 1.0000000000000002e-06;
    activation_units[37] = 0.50321666580471969;
    phase_units[37]      = pow(10,-12.000000);
    is_PD[37] = 0;
    nTB[37] = 0;

    // (29):  H + CH2 (+M) <=> CH3 (+M)
    // (29):  H + CH2 (+M) <=> CH3 (+M)
    fwd_A[0]     = 25000000000000000;
    fwd_beta[0]  = -0.80000000000000004;
    fwd_Ea[0]    = 0;
    low_A[0]     = 3.2000000000000002e+27;
    low_beta[0]  = -3.1400000000000001;
    low_Ea[0]    = 1230;
    troe_a[0]    = 0.68000000000000005;
    troe_Tsss[0] = 78;
    troe_Ts[0]   = 1995;
    troe_Tss[0]  = 5590;
    troe_len[0]  = 4;
    prefactor_units[0]  = 1.0000000000000002e-06;
    activation_units[0] = 0.50321666580471969;
    phase_units[0]      = pow(10,-12.000000);
    is_PD[0] = 1;
    nTB[0] = 7;
    TB[0] = (amrex::Real *) malloc(7 * sizeof(amrex::Real));
    TBid[0] = (int *) malloc(7 * sizeof(int));
    TBid[0][0] = 0; TB[0][0] = 2; // H2
    TBid[0][1] = 5; TB[0][1] = 6; // H2O
    TBid[0][2] = 10; TB[0][2] = 2; // CH4
    TBid[0][3] = 11; TB[0][3] = 1.5; // CO
    TBid[0][4] = 12; TB[0][4] = 2; // CO2
    TBid[0][5] = 18; TB[0][5] = 3; // C2H6
    TBid[0][6] = 20; TB[0][6] = 0.69999999999999996; // AR

    // (30):  H + CH3 (+M) <=> CH4 (+M)
    // (30):  H + CH3 (+M) <=> CH4 (+M)
    fwd_A[1]     = 12700000000000000;
    fwd_beta[1]  = -0.63;
    fwd_Ea[1]    = 383;
    low_A[1]     = 2.4769999999999999e+33;
    low_beta[1]  = -4.7599999999999998;
    low_Ea[1]    = 2440;
    troe_a[1]    = 0.78300000000000003;
    troe_Tsss[1] = 74;
    troe_Ts[1]   = 2941;
    troe_Tss[1]  = 6964;
    troe_len[1]  = 4;
    prefactor_units[1]  = 1.0000000000000002e-06;
    activation_units[1] = 0.50321666580471969;
    phase_units[1]      = pow(10,-12.000000);
    is_PD[1] = 1;
    nTB[1] = 7;
    TB[1] = (amrex::Real *) malloc(7 * sizeof(amrex::Real));
    TBid[1] = (int *) malloc(7 * sizeof(int));
    TBid[1][0] = 0; TB[1][0] = 2; // H2
    TBid[1][1] = 5; TB[1][1] = 6; // H2O
    TBid[1][2] = 10; TB[1][2] = 2; // CH4
    TBid[1][3] = 11; TB[1][3] = 1.5; // CO
    TBid[1][4] = 12; TB[1][4] = 2; // CO2
    TBid[1][5] = 18; TB[1][5] = 3; // C2H6
    TBid[1][6] = 20; TB[1][6] = 0.69999999999999996; // AR

    // (31):  H + CH4 <=> CH3 + H2
    // (31):  H + CH4 <=> CH3 + H2
    fwd_A[38]     = 660000000;
    fwd_beta[38]  = 1.6200000000000001;
    fwd_Ea[38]    = 10840;
    prefactor_units[38]  = 1.0000000000000002e-06;
    activation_units[38] = 0.50321666580471969;
    phase_units[38]      = pow(10,-12.000000);
    is_PD[38] = 0;
    nTB[38] = 0;

    // (32):  H + HCO (+M) <=> CH2O (+M)
    // (32):  H + HCO (+M) <=> CH2O (+M)
    fwd_A[2]     = 1090000000000;
    fwd_beta[2]  = 0.47999999999999998;
    fwd_Ea[2]    = -260;
    low_A[2]     = 1.35e+24;
    low_beta[2]  = -2.5699999999999998;
    low_Ea[2]    = 1425;
    troe_a[2]    = 0.78239999999999998;
    troe_Tsss[2] = 271;
    troe_Ts[2]   = 2755;
    troe_Tss[2]  = 6570;
    troe_len[2]  = 4;
    prefactor_units[2]  = 1.0000000000000002e-06;
    activation_units[2] = 0.50321666580471969;
    phase_units[2]      = pow(10,-12.000000);
    is_PD[2] = 1;
    nTB[2] = 7;
    TB[2] = (amrex::Real *) malloc(7 * sizeof(amrex::Real));
    TBid[2] = (int *) malloc(7 * sizeof(int));
    TBid[2][0] = 0; TB[2][0] = 2; // H2
    TBid[2][1] = 5; TB[2][1] = 6; // H2O
    TBid[2][2] = 10; TB[2][2] = 2; // CH4
    TBid[2][3] = 11; TB[2][3] = 1.5; // CO
    TBid[2][4] = 12; TB[2][4] = 2; // CO2
    TBid[2][5] = 18; TB[2][5] = 3; // C2H6
    TBid[2][6] = 20; TB[2][6] = 0.69999999999999996; // AR

    // (33):  H + HCO <=> H2 + CO
    // (33):  H + HCO <=> H2 + CO
    fwd_A[39]     = 73400000000000;
    fwd_beta[39]  = 0;
    fwd_Ea[39]    = 0;
    prefactor_units[39]  = 1.0000000000000002e-06;
    activation_units[39] = 0.50321666580471969;
    phase_units[39]      = pow(10,-12.000000);
    is_PD[39] = 0;
    nTB[39] = 0;

    // (34):  H + CH2O (+M) <=> CH3O (+M)
    // (34):  H + CH2O (+M) <=> CH3O (+M)
    fwd_A[3]     = 540000000000;
    fwd_beta[3]  = 0.45400000000000001;
    fwd_Ea[3]    = 2600;
    low_A[3]     = 2.2e+30;
    low_beta[3]  = -4.7999999999999998;
    low_Ea[3]    = 5560;
    troe_a[3]    = 0.75800000000000001;
    troe_Tsss[3] = 94;
    troe_Ts[3]   = 1555;
    troe_Tss[3]  = 4200;
    troe_len[3]  = 4;
    prefactor_units[3]  = 1.0000000000000002e-06;
    activation_units[3] = 0.50321666580471969;
    phase_units[3]      = pow(10,-12.000000);
    is_PD[3] = 1;
    nTB[3] = 6;
    TB[3] = (amrex::Real *) malloc(6 * sizeof(amrex::Real));
    TBid[3] = (int *) malloc(6 * sizeof(int));
    TBid[3][0] = 0; TB[3][0] = 2; // H2
    TBid[3][1] = 5; TB[3][1] = 6; // H2O
    TBid[3][2] = 10; TB[3][2] = 2; // CH4
    TBid[3][3] = 11; TB[3][3] = 1.5; // CO
    TBid[3][4] = 12; TB[3][4] = 2; // CO2
    TBid[3][5] = 18; TB[3][5] = 3; // C2H6

    // (35):  H + CH2O <=> HCO + H2
    // (35):  H + CH2O <=> HCO + H2
    fwd_A[40]     = 23000000000;
    fwd_beta[40]  = 1.05;
    fwd_Ea[40]    = 3275;
    prefactor_units[40]  = 1.0000000000000002e-06;
    activation_units[40] = 0.50321666580471969;
    phase_units[40]      = pow(10,-12.000000);
    is_PD[40] = 0;
    nTB[40] = 0;

    // (36):  H + CH3O <=> OH + CH3
    // (36):  H + CH3O <=> OH + CH3
    fwd_A[41]     = 32000000000000;
    fwd_beta[41]  = 0;
    fwd_Ea[41]    = 0;
    prefactor_units[41]  = 1.0000000000000002e-06;
    activation_units[41] = 0.50321666580471969;
    phase_units[41]      = pow(10,-12.000000);
    is_PD[41] = 0;
    nTB[41] = 0;

    // (37):  H + C2H4 (+M) <=> C2H5 (+M)
    // (37):  H + C2H4 (+M) <=> C2H5 (+M)
    fwd_A[4]     = 1080000000000;
    fwd_beta[4]  = 0.45400000000000001;
    fwd_Ea[4]    = 1820;
    low_A[4]     = 1.1999999999999999e+42;
    low_beta[4]  = -7.6200000000000001;
    low_Ea[4]    = 6970;
    troe_a[4]    = 0.97529999999999994;
    troe_Tsss[4] = 210;
    troe_Ts[4]   = 984;
    troe_Tss[4]  = 4374;
    troe_len[4]  = 4;
    prefactor_units[4]  = 1.0000000000000002e-06;
    activation_units[4] = 0.50321666580471969;
    phase_units[4]      = pow(10,-12.000000);
    is_PD[4] = 1;
    nTB[4] = 7;
    TB[4] = (amrex::Real *) malloc(7 * sizeof(amrex::Real));
    TBid[4] = (int *) malloc(7 * sizeof(int));
    TBid[4][0] = 0; TB[4][0] = 2; // H2
    TBid[4][1] = 5; TB[4][1] = 6; // H2O
    TBid[4][2] = 10; TB[4][2] = 2; // CH4
    TBid[4][3] = 11; TB[4][3] = 1.5; // CO
    TBid[4][4] = 12; TB[4][4] = 2; // CO2
    TBid[4][5] = 18; TB[4][5] = 3; // C2H6
    TBid[4][6] = 20; TB[4][6] = 0.69999999999999996; // AR

    // (38):  H + C2H5 (+M) <=> C2H6 (+M)
    // (38):  H + C2H5 (+M) <=> C2H6 (+M)
    fwd_A[5]     = 5.21e+17;
    fwd_beta[5]  = -0.98999999999999999;
    fwd_Ea[5]    = 1580;
    low_A[5]     = 1.9900000000000001e+41;
    low_beta[5]  = -7.0800000000000001;
    low_Ea[5]    = 6685;
    troe_a[5]    = 0.84219999999999995;
    troe_Tsss[5] = 125;
    troe_Ts[5]   = 2219;
    troe_Tss[5]  = 6882;
    troe_len[5]  = 4;
    prefactor_units[5]  = 1.0000000000000002e-06;
    activation_units[5] = 0.50321666580471969;
    phase_units[5]      = pow(10,-12.000000);
    is_PD[5] = 1;
    nTB[5] = 7;
    TB[5] = (amrex::Real *) malloc(7 * sizeof(amrex::Real));
    TBid[5] = (int *) malloc(7 * sizeof(int));
    TBid[5][0] = 0; TB[5][0] = 2; // H2
    TBid[5][1] = 5; TB[5][1] = 6; // H2O
    TBid[5][2] = 10; TB[5][2] = 2; // CH4
    TBid[5][3] = 11; TB[5][3] = 1.5; // CO
    TBid[5][4] = 12; TB[5][4] = 2; // CO2
    TBid[5][5] = 18; TB[5][5] = 3; // C2H6
    TBid[5][6] = 20; TB[5][6] = 0.69999999999999996; // AR

    // (39):  H + C2H6 <=> C2H5 + H2
    // (39):  H + C2H6 <=> C2H5 + H2
    fwd_A[42]     = 115000000;
    fwd_beta[42]  = 1.8999999999999999;
    fwd_Ea[42]    = 7530;
    prefactor_units[42]  = 1.0000000000000002e-06;
    activation_units[42] = 0.50321666580471969;
    phase_units[42]      = pow(10,-12.000000);
    is_PD[42] = 0;
    nTB[42] = 0;

    // (40):  H2 + CO (+M) <=> CH2O (+M)
    // (40):  H2 + CO (+M) <=> CH2O (+M)
    fwd_A[6]     = 43000000;
    fwd_beta[6]  = 1.5;
    fwd_Ea[6]    = 79600;
    low_A[6]     = 5.0699999999999998e+27;
    low_beta[6]  = -3.4199999999999999;
    low_Ea[6]    = 84350;
    troe_a[6]    = 0.93200000000000005;
    troe_Tsss[6] = 197;
    troe_Ts[6]   = 1540;
    troe_Tss[6]  = 10300;
    troe_len[6]  = 4;
    prefactor_units[6]  = 1.0000000000000002e-06;
    activation_units[6] = 0.50321666580471969;
    phase_units[6]      = pow(10,-12.000000);
    is_PD[6] = 1;
    nTB[6] = 7;
    TB[6] = (amrex::Real *) malloc(7 * sizeof(amrex::Real));
    TBid[6] = (int *) malloc(7 * sizeof(int));
    TBid[6][0] = 0; TB[6][0] = 2; // H2
    TBid[6][1] = 5; TB[6][1] = 6; // H2O
    TBid[6][2] = 10; TB[6][2] = 2; // CH4
    TBid[6][3] = 11; TB[6][3] = 1.5; // CO
    TBid[6][4] = 12; TB[6][4] = 2; // CO2
    TBid[6][5] = 18; TB[6][5] = 3; // C2H6
    TBid[6][6] = 20; TB[6][6] = 0.69999999999999996; // AR

    // (41):  OH + H2 <=> H + H2O
    // (41):  OH + H2 <=> H + H2O
    fwd_A[43]     = 216000000;
    fwd_beta[43]  = 1.51;
    fwd_Ea[43]    = 3430;
    prefactor_units[43]  = 1.0000000000000002e-06;
    activation_units[43] = 0.50321666580471969;
    phase_units[43]      = pow(10,-12.000000);
    is_PD[43] = 0;
    nTB[43] = 0;

    // (42):  2.000000 OH <=> O + H2O
    // (42):  2.000000 OH <=> O + H2O
    fwd_A[44]     = 35700;
    fwd_beta[44]  = 2.3999999999999999;
    fwd_Ea[44]    = -2110;
    prefactor_units[44]  = 1.0000000000000002e-06;
    activation_units[44] = 0.50321666580471969;
    phase_units[44]      = pow(10,-12.000000);
    is_PD[44] = 0;
    nTB[44] = 0;

    // (43):  OH + HO2 <=> O2 + H2O
    // (43):  OH + HO2 <=> O2 + H2O
    fwd_A[45]     = 29000000000000;
    fwd_beta[45]  = 0;
    fwd_Ea[45]    = -500;
    prefactor_units[45]  = 1.0000000000000002e-06;
    activation_units[45] = 0.50321666580471969;
    phase_units[45]      = pow(10,-12.000000);
    is_PD[45] = 0;
    nTB[45] = 0;

    // (44):  OH + CH2 <=> H + CH2O
    // (44):  OH + CH2 <=> H + CH2O
    fwd_A[46]     = 20000000000000;
    fwd_beta[46]  = 0;
    fwd_Ea[46]    = 0;
    prefactor_units[46]  = 1.0000000000000002e-06;
    activation_units[46] = 0.50321666580471969;
    phase_units[46]      = pow(10,-12.000000);
    is_PD[46] = 0;
    nTB[46] = 0;

    // (45):  OH + CH2(S) <=> H + CH2O
    // (45):  OH + CH2(S) <=> H + CH2O
    fwd_A[47]     = 30000000000000;
    fwd_beta[47]  = 0;
    fwd_Ea[47]    = 0;
    prefactor_units[47]  = 1.0000000000000002e-06;
    activation_units[47] = 0.50321666580471969;
    phase_units[47]      = pow(10,-12.000000);
    is_PD[47] = 0;
    nTB[47] = 0;

    // (46):  OH + CH3 <=> CH2 + H2O
    // (46):  OH + CH3 <=> CH2 + H2O
    fwd_A[48]     = 56000000;
    fwd_beta[48]  = 1.6000000000000001;
    fwd_Ea[48]    = 5420;
    prefactor_units[48]  = 1.0000000000000002e-06;
    activation_units[48] = 0.50321666580471969;
    phase_units[48]      = pow(10,-12.000000);
    is_PD[48] = 0;
    nTB[48] = 0;

    // (47):  OH + CH3 <=> CH2(S) + H2O
    // (47):  OH + CH3 <=> CH2(S) + H2O
    fwd_A[49]     = 25010000000000;
    fwd_beta[49]  = 0;
    fwd_Ea[49]    = 0;
    prefactor_units[49]  = 1.0000000000000002e-06;
    activation_units[49] = 0.50321666580471969;
    phase_units[49]      = pow(10,-12.000000);
    is_PD[49] = 0;
    nTB[49] = 0;

    // (48):  OH + CH4 <=> CH3 + H2O
    // (48):  OH + CH4 <=> CH3 + H2O
    fwd_A[50]     = 100000000;
    fwd_beta[50]  = 1.6000000000000001;
    fwd_Ea[50]    = 3120;
    prefactor_units[50]  = 1.0000000000000002e-06;
    activation_units[50] = 0.50321666580471969;
    phase_units[50]      = pow(10,-12.000000);
    is_PD[50] = 0;
    nTB[50] = 0;

    // (49):  OH + CO <=> H + CO2
    // (49):  OH + CO <=> H + CO2
    fwd_A[51]     = 47600000;
    fwd_beta[51]  = 1.228;
    fwd_Ea[51]    = 70;
    prefactor_units[51]  = 1.0000000000000002e-06;
    activation_units[51] = 0.50321666580471969;
    phase_units[51]      = pow(10,-12.000000);
    is_PD[51] = 0;
    nTB[51] = 0;

    // (50):  OH + HCO <=> H2O + CO
    // (50):  OH + HCO <=> H2O + CO
    fwd_A[52]     = 50000000000000;
    fwd_beta[52]  = 0;
    fwd_Ea[52]    = 0;
    prefactor_units[52]  = 1.0000000000000002e-06;
    activation_units[52] = 0.50321666580471969;
    phase_units[52]      = pow(10,-12.000000);
    is_PD[52] = 0;
    nTB[52] = 0;

    // (51):  OH + CH2O <=> HCO + H2O
    // (51):  OH + CH2O <=> HCO + H2O
    fwd_A[53]     = 3430000000;
    fwd_beta[53]  = 1.1799999999999999;
    fwd_Ea[53]    = -447;
    prefactor_units[53]  = 1.0000000000000002e-06;
    activation_units[53] = 0.50321666580471969;
    phase_units[53]      = pow(10,-12.000000);
    is_PD[53] = 0;
    nTB[53] = 0;

    // (52):  OH + C2H6 <=> C2H5 + H2O
    // (52):  OH + C2H6 <=> C2H5 + H2O
    fwd_A[54]     = 3540000;
    fwd_beta[54]  = 2.1200000000000001;
    fwd_Ea[54]    = 870;
    prefactor_units[54]  = 1.0000000000000002e-06;
    activation_units[54] = 0.50321666580471969;
    phase_units[54]      = pow(10,-12.000000);
    is_PD[54] = 0;
    nTB[54] = 0;

    // (53):  HO2 + CH2 <=> OH + CH2O
    // (53):  HO2 + CH2 <=> OH + CH2O
    fwd_A[55]     = 20000000000000;
    fwd_beta[55]  = 0;
    fwd_Ea[55]    = 0;
    prefactor_units[55]  = 1.0000000000000002e-06;
    activation_units[55] = 0.50321666580471969;
    phase_units[55]      = pow(10,-12.000000);
    is_PD[55] = 0;
    nTB[55] = 0;

    // (54):  HO2 + CH3 <=> O2 + CH4
    // (54):  HO2 + CH3 <=> O2 + CH4
    fwd_A[56]     = 1000000000000;
    fwd_beta[56]  = 0;
    fwd_Ea[56]    = 0;
    prefactor_units[56]  = 1.0000000000000002e-06;
    activation_units[56] = 0.50321666580471969;
    phase_units[56]      = pow(10,-12.000000);
    is_PD[56] = 0;
    nTB[56] = 0;

    // (55):  HO2 + CH3 <=> OH + CH3O
    // (55):  HO2 + CH3 <=> OH + CH3O
    fwd_A[57]     = 20000000000000;
    fwd_beta[57]  = 0;
    fwd_Ea[57]    = 0;
    prefactor_units[57]  = 1.0000000000000002e-06;
    activation_units[57] = 0.50321666580471969;
    phase_units[57]      = pow(10,-12.000000);
    is_PD[57] = 0;
    nTB[57] = 0;

    // (56):  HO2 + CO <=> OH + CO2
    // (56):  HO2 + CO <=> OH + CO2
    fwd_A[58]     = 150000000000000;
    fwd_beta[58]  = 0;
    fwd_Ea[58]    = 23600;
    prefactor_units[58]  = 1.0000000000000002e-06;
    activation_units[58] = 0.50321666580471969;
    phase_units[58]      = pow(10,-12.000000);
    is_PD[58] = 0;
    nTB[58] = 0;

    // (57):  CH2 + O2 <=> OH + HCO
    // (57):  CH2 + O2 <=> OH + HCO
    fwd_A[59]     = 13200000000000;
    fwd_beta[59]  = 0;
    fwd_Ea[59]    = 1500;
    prefactor_units[59]  = 1.0000000000000002e-06;
    activation_units[59] = 0.50321666580471969;
    phase_units[59]      = pow(10,-12.000000);
    is_PD[59] = 0;
    nTB[59] = 0;

    // (58):  CH2 + H2 <=> H + CH3
    // (58):  CH2 + H2 <=> H + CH3
    fwd_A[60]     = 500000;
    fwd_beta[60]  = 2;
    fwd_Ea[60]    = 7230;
    prefactor_units[60]  = 1.0000000000000002e-06;
    activation_units[60] = 0.50321666580471969;
    phase_units[60]      = pow(10,-12.000000);
    is_PD[60] = 0;
    nTB[60] = 0;

    // (59):  CH2 + CH3 <=> H + C2H4
    // (59):  CH2 + CH3 <=> H + C2H4
    fwd_A[61]     = 40000000000000;
    fwd_beta[61]  = 0;
    fwd_Ea[61]    = 0;
    prefactor_units[61]  = 1.0000000000000002e-06;
    activation_units[61] = 0.50321666580471969;
    phase_units[61]      = pow(10,-12.000000);
    is_PD[61] = 0;
    nTB[61] = 0;

    // (60):  CH2 + CH4 <=> 2.000000 CH3
    // (60):  CH2 + CH4 <=> 2.000000 CH3
    fwd_A[62]     = 2460000;
    fwd_beta[62]  = 2;
    fwd_Ea[62]    = 8270;
    prefactor_units[62]  = 1.0000000000000002e-06;
    activation_units[62] = 0.50321666580471969;
    phase_units[62]      = pow(10,-12.000000);
    is_PD[62] = 0;
    nTB[62] = 0;

    // (61):  CH2(S) + N2 <=> CH2 + N2
    // (61):  CH2(S) + N2 <=> CH2 + N2
    fwd_A[63]     = 15000000000000;
    fwd_beta[63]  = 0;
    fwd_Ea[63]    = 600;
    prefactor_units[63]  = 1.0000000000000002e-06;
    activation_units[63] = 0.50321666580471969;
    phase_units[63]      = pow(10,-12.000000);
    is_PD[63] = 0;
    nTB[63] = 0;

    // (62):  CH2(S) + AR <=> CH2 + AR
    // (62):  CH2(S) + AR <=> CH2 + AR
    fwd_A[64]     = 9000000000000;
    fwd_beta[64]  = 0;
    fwd_Ea[64]    = 600;
    prefactor_units[64]  = 1.0000000000000002e-06;
    activation_units[64] = 0.50321666580471969;
    phase_units[64]      = pow(10,-12.000000);
    is_PD[64] = 0;
    nTB[64] = 0;

    // (63):  CH2(S) + O2 <=> H + OH + CO
    // (63):  CH2(S) + O2 <=> H + OH + CO
    fwd_A[65]     = 28000000000000;
    fwd_beta[65]  = 0;
    fwd_Ea[65]    = 0;
    prefactor_units[65]  = 1.0000000000000002e-06;
    activation_units[65] = 0.50321666580471969;
    phase_units[65]      = pow(10,-12.000000);
    is_PD[65] = 0;
    nTB[65] = 0;

    // (64):  CH2(S) + O2 <=> CO + H2O
    // (64):  CH2(S) + O2 <=> CO + H2O
    fwd_A[66]     = 12000000000000;
    fwd_beta[66]  = 0;
    fwd_Ea[66]    = 0;
    prefactor_units[66]  = 1.0000000000000002e-06;
    activation_units[66] = 0.50321666580471969;
    phase_units[66]      = pow(10,-12.000000);
    is_PD[66] = 0;
    nTB[66] = 0;

    // (65):  CH2(S) + H2 <=> CH3 + H
    // (65):  CH2(S) + H2 <=> CH3 + H
    fwd_A[67]     = 70000000000000;
    fwd_beta[67]  = 0;
    fwd_Ea[67]    = 0;
    prefactor_units[67]  = 1.0000000000000002e-06;
    activation_units[67] = 0.50321666580471969;
    phase_units[67]      = pow(10,-12.000000);
    is_PD[67] = 0;
    nTB[67] = 0;

    // (66):  CH2(S) + H2O <=> CH2 + H2O
    // (66):  CH2(S) + H2O <=> CH2 + H2O
    fwd_A[68]     = 30000000000000;
    fwd_beta[68]  = 0;
    fwd_Ea[68]    = 0;
    prefactor_units[68]  = 1.0000000000000002e-06;
    activation_units[68] = 0.50321666580471969;
    phase_units[68]      = pow(10,-12.000000);
    is_PD[68] = 0;
    nTB[68] = 0;

    // (67):  CH2(S) + CH3 <=> H + C2H4
    // (67):  CH2(S) + CH3 <=> H + C2H4
    fwd_A[69]     = 12000000000000;
    fwd_beta[69]  = 0;
    fwd_Ea[69]    = -570;
    prefactor_units[69]  = 1.0000000000000002e-06;
    activation_units[69] = 0.50321666580471969;
    phase_units[69]      = pow(10,-12.000000);
    is_PD[69] = 0;
    nTB[69] = 0;

    // (68):  CH2(S) + CH4 <=> 2.000000 CH3
    // (68):  CH2(S) + CH4 <=> 2.000000 CH3
    fwd_A[70]     = 16000000000000;
    fwd_beta[70]  = 0;
    fwd_Ea[70]    = -570;
    prefactor_units[70]  = 1.0000000000000002e-06;
    activation_units[70] = 0.50321666580471969;
    phase_units[70]      = pow(10,-12.000000);
    is_PD[70] = 0;
    nTB[70] = 0;

    // (69):  CH2(S) + CO <=> CH2 + CO
    // (69):  CH2(S) + CO <=> CH2 + CO
    fwd_A[71]     = 9000000000000;
    fwd_beta[71]  = 0;
    fwd_Ea[71]    = 0;
    prefactor_units[71]  = 1.0000000000000002e-06;
    activation_units[71] = 0.50321666580471969;
    phase_units[71]      = pow(10,-12.000000);
    is_PD[71] = 0;
    nTB[71] = 0;

    // (70):  CH2(S) + CO2 <=> CH2 + CO2
    // (70):  CH2(S) + CO2 <=> CH2 + CO2
    fwd_A[72]     = 7000000000000;
    fwd_beta[72]  = 0;
    fwd_Ea[72]    = 0;
    prefactor_units[72]  = 1.0000000000000002e-06;
    activation_units[72] = 0.50321666580471969;
    phase_units[72]      = pow(10,-12.000000);
    is_PD[72] = 0;
    nTB[72] = 0;

    // (71):  CH2(S) + CO2 <=> CO + CH2O
    // (71):  CH2(S) + CO2 <=> CO + CH2O
    fwd_A[73]     = 14000000000000;
    fwd_beta[73]  = 0;
    fwd_Ea[73]    = 0;
    prefactor_units[73]  = 1.0000000000000002e-06;
    activation_units[73] = 0.50321666580471969;
    phase_units[73]      = pow(10,-12.000000);
    is_PD[73] = 0;
    nTB[73] = 0;

    // (72):  CH3 + O2 <=> O + CH3O
    // (72):  CH3 + O2 <=> O + CH3O
    fwd_A[74]     = 26750000000000;
    fwd_beta[74]  = 0;
    fwd_Ea[74]    = 28800;
    prefactor_units[74]  = 1.0000000000000002e-06;
    activation_units[74] = 0.50321666580471969;
    phase_units[74]      = pow(10,-12.000000);
    is_PD[74] = 0;
    nTB[74] = 0;

    // (73):  CH3 + O2 <=> OH + CH2O
    // (73):  CH3 + O2 <=> OH + CH2O
    fwd_A[75]     = 36000000000;
    fwd_beta[75]  = 0;
    fwd_Ea[75]    = 8940;
    prefactor_units[75]  = 1.0000000000000002e-06;
    activation_units[75] = 0.50321666580471969;
    phase_units[75]      = pow(10,-12.000000);
    is_PD[75] = 0;
    nTB[75] = 0;

    // (74):  2.000000 CH3 (+M) <=> C2H6 (+M)
    // (74):  2.000000 CH3 (+M) <=> C2H6 (+M)
    fwd_A[7]     = 21200000000000000;
    fwd_beta[7]  = -0.96999999999999997;
    fwd_Ea[7]    = 620;
    low_A[7]     = 1.7700000000000001e+50;
    low_beta[7]  = -9.6699999999999999;
    low_Ea[7]    = 6220;
    troe_a[7]    = 0.53249999999999997;
    troe_Tsss[7] = 151;
    troe_Ts[7]   = 1038;
    troe_Tss[7]  = 4970;
    troe_len[7]  = 4;
    prefactor_units[7]  = 1.0000000000000002e-06;
    activation_units[7] = 0.50321666580471969;
    phase_units[7]      = pow(10,-12.000000);
    is_PD[7] = 1;
    nTB[7] = 7;
    TB[7] = (amrex::Real *) malloc(7 * sizeof(amrex::Real));
    TBid[7] = (int *) malloc(7 * sizeof(int));
    TBid[7][0] = 0; TB[7][0] = 2; // H2
    TBid[7][1] = 5; TB[7][1] = 6; // H2O
    TBid[7][2] = 10; TB[7][2] = 2; // CH4
    TBid[7][3] = 11; TB[7][3] = 1.5; // CO
    TBid[7][4] = 12; TB[7][4] = 2; // CO2
    TBid[7][5] = 18; TB[7][5] = 3; // C2H6
    TBid[7][6] = 20; TB[7][6] = 0.69999999999999996; // AR

    // (75):  2.000000 CH3 <=> H + C2H5
    // (75):  2.000000 CH3 <=> H + C2H5
    fwd_A[76]     = 4990000000000;
    fwd_beta[76]  = 0.10000000000000001;
    fwd_Ea[76]    = 10600;
    prefactor_units[76]  = 1.0000000000000002e-06;
    activation_units[76] = 0.50321666580471969;
    phase_units[76]      = pow(10,-12.000000);
    is_PD[76] = 0;
    nTB[76] = 0;

    // (76):  CH3 + HCO <=> CH4 + CO
    // (76):  CH3 + HCO <=> CH4 + CO
    fwd_A[77]     = 26480000000000;
    fwd_beta[77]  = 0;
    fwd_Ea[77]    = 0;
    prefactor_units[77]  = 1.0000000000000002e-06;
    activation_units[77] = 0.50321666580471969;
    phase_units[77]      = pow(10,-12.000000);
    is_PD[77] = 0;
    nTB[77] = 0;

    // (77):  CH3 + CH2O <=> HCO + CH4
    // (77):  CH3 + CH2O <=> HCO + CH4
    fwd_A[78]     = 3320;
    fwd_beta[78]  = 2.8100000000000001;
    fwd_Ea[78]    = 5860;
    prefactor_units[78]  = 1.0000000000000002e-06;
    activation_units[78] = 0.50321666580471969;
    phase_units[78]      = pow(10,-12.000000);
    is_PD[78] = 0;
    nTB[78] = 0;

    // (78):  CH3 + C2H6 <=> C2H5 + CH4
    // (78):  CH3 + C2H6 <=> C2H5 + CH4
    fwd_A[79]     = 6140000;
    fwd_beta[79]  = 1.74;
    fwd_Ea[79]    = 10450;
    prefactor_units[79]  = 1.0000000000000002e-06;
    activation_units[79] = 0.50321666580471969;
    phase_units[79]      = pow(10,-12.000000);
    is_PD[79] = 0;
    nTB[79] = 0;

    // (79):  HCO + H2O <=> H + CO + H2O
    // (79):  HCO + H2O <=> H + CO + H2O
    fwd_A[80]     = 2.244e+18;
    fwd_beta[80]  = -1;
    fwd_Ea[80]    = 17000;
    prefactor_units[80]  = 1.0000000000000002e-06;
    activation_units[80] = 0.50321666580471969;
    phase_units[80]      = pow(10,-12.000000);
    is_PD[80] = 0;
    nTB[80] = 0;

    // (80):  HCO + M <=> H + CO + M
    // (80):  HCO + M <=> H + CO + M
    fwd_A[13]     = 1.87e+17;
    fwd_beta[13]  = -1;
    fwd_Ea[13]    = 17000;
    prefactor_units[13]  = 1.0000000000000002e-06;
    activation_units[13] = 0.50321666580471969;
    phase_units[13]      = pow(10,-6.000000);
    is_PD[13] = 0;
    nTB[13] = 6;
    TB[13] = (amrex::Real *) malloc(6 * sizeof(amrex::Real));
    TBid[13] = (int *) malloc(6 * sizeof(int));
    TBid[13][0] = 0; TB[13][0] = 2; // H2
    TBid[13][1] = 5; TB[13][1] = 0; // H2O
    TBid[13][2] = 10; TB[13][2] = 2; // CH4
    TBid[13][3] = 11; TB[13][3] = 1.5; // CO
    TBid[13][4] = 12; TB[13][4] = 2; // CO2
    TBid[13][5] = 18; TB[13][5] = 3; // C2H6

    // (81):  HCO + O2 <=> HO2 + CO
    // (81):  HCO + O2 <=> HO2 + CO
    fwd_A[81]     = 7600000000000;
    fwd_beta[81]  = 0;
    fwd_Ea[81]    = 400;
    prefactor_units[81]  = 1.0000000000000002e-06;
    activation_units[81] = 0.50321666580471969;
    phase_units[81]      = pow(10,-12.000000);
    is_PD[81] = 0;
    nTB[81] = 0;

    // (82):  CH3O + O2 <=> HO2 + CH2O
    // (82):  CH3O + O2 <=> HO2 + CH2O
    fwd_A[82]     = 4.2799999999999999e-13;
    fwd_beta[82]  = 7.5999999999999996;
    fwd_Ea[82]    = -3530;
    prefactor_units[82]  = 1.0000000000000002e-06;
    activation_units[82] = 0.50321666580471969;
    phase_units[82]      = pow(10,-12.000000);
    is_PD[82] = 0;
    nTB[82] = 0;

    // (83):  C2H5 + O2 <=> HO2 + C2H4
    // (83):  C2H5 + O2 <=> HO2 + C2H4
    fwd_A[83]     = 840000000000;
    fwd_beta[83]  = 0;
    fwd_Ea[83]    = 3875;
    prefactor_units[83]  = 1.0000000000000002e-06;
    activation_units[83] = 0.50321666580471969;
    phase_units[83]      = pow(10,-12.000000);
    is_PD[83] = 0;
    nTB[83] = 0;

}


/* Finalizes parameter database */
void CKFINALIZE()
{
  for (int i=0; i<84; ++i) {
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

    /*CH2 */
    ncf[ 7 * kd + 2 ] = 1; /*C */
    ncf[ 7 * kd + 1 ] = 2; /*H */

    /*CH2(S) */
    ncf[ 8 * kd + 2 ] = 1; /*C */
    ncf[ 8 * kd + 1 ] = 2; /*H */

    /*CH3 */
    ncf[ 9 * kd + 2 ] = 1; /*C */
    ncf[ 9 * kd + 1 ] = 3; /*H */

    /*CH4 */
    ncf[ 10 * kd + 2 ] = 1; /*C */
    ncf[ 10 * kd + 1 ] = 4; /*H */

    /*CO */
    ncf[ 11 * kd + 2 ] = 1; /*C */
    ncf[ 11 * kd + 0 ] = 1; /*O */

    /*CO2 */
    ncf[ 12 * kd + 2 ] = 1; /*C */
    ncf[ 12 * kd + 0 ] = 2; /*O */

    /*HCO */
    ncf[ 13 * kd + 1 ] = 1; /*H */
    ncf[ 13 * kd + 2 ] = 1; /*C */
    ncf[ 13 * kd + 0 ] = 1; /*O */

    /*CH2O */
    ncf[ 14 * kd + 1 ] = 2; /*H */
    ncf[ 14 * kd + 2 ] = 1; /*C */
    ncf[ 14 * kd + 0 ] = 1; /*O */

    /*CH3O */
    ncf[ 15 * kd + 2 ] = 1; /*C */
    ncf[ 15 * kd + 1 ] = 3; /*H */
    ncf[ 15 * kd + 0 ] = 1; /*O */

    /*C2H4 */
    ncf[ 16 * kd + 2 ] = 2; /*C */
    ncf[ 16 * kd + 1 ] = 4; /*H */

    /*C2H5 */
    ncf[ 17 * kd + 2 ] = 2; /*C */
    ncf[ 17 * kd + 1 ] = 5; /*H */

    /*C2H6 */
    ncf[ 18 * kd + 2 ] = 2; /*C */
    ncf[ 18 * kd + 1 ] = 6; /*H */

    /*N2 */
    ncf[ 19 * kd + 3 ] = 2; /*N */

    /*AR */
    ncf[ 20 * kd + 4 ] = 1; /*AR */

}


/* Returns the vector of strings of element names */
void CKSYME_STR(amrex::Vector<std::string>& ename)
{
    ename.resize(5);
    ename[0] = "O";
    ename[1] = "H";
    ename[2] = "C";
    ename[3] = "N";
    ename[4] = "AR";
}


/* Returns the vector of strings of species names */
void CKSYMS_STR(amrex::Vector<std::string>& kname)
{
    kname.resize(21);
    kname[0] = "H2";
    kname[1] = "H";
    kname[2] = "O";
    kname[3] = "O2";
    kname[4] = "OH";
    kname[5] = "H2O";
    kname[6] = "HO2";
    kname[7] = "CH2";
    kname[8] = "CH2(S)";
    kname[9] = "CH3";
    kname[10] = "CH4";
    kname[11] = "CO";
    kname[12] = "CO2";
    kname[13] = "HCO";
    kname[14] = "CH2O";
    kname[15] = "CH3O";
    kname[16] = "C2H4";
    kname[17] = "C2H5";
    kname[18] = "C2H6";
    kname[19] = "N2";
    kname[20] = "AR";
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
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
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
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
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
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
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
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
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
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
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
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
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
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
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
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
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

