#ifndef MECHANISM_CPP
#define MECHANISM_CPP

#include "chemistry_file.H"

/* PURE CPU stuff  */
#ifndef AMREX_USE_GPU
namespace thermo
{
    amrex::Real fwd_A[175], fwd_beta[175], fwd_Ea[175];
    amrex::Real low_A[175], low_beta[175], low_Ea[175];
    amrex::Real rev_A[175], rev_beta[175], rev_Ea[175];
    amrex::Real troe_a[175],troe_Ts[175], troe_Tss[175], troe_Tsss[175];
    amrex::Real sri_a[175], sri_b[175], sri_c[175], sri_d[175], sri_e[175];
    amrex::Real activation_units[175], prefactor_units[175], phase_units[175];
    int is_PD[175], troe_len[175], sri_len[175], nTB[175], *TBid[175];
    amrex::Real *TB[175];
    std::vector<std::vector<amrex::Real>> kiv(175); 
    std::vector<std::vector<amrex::Real>> nuv(175); 
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
        if (*i > 175) {
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
    amrex::Real imw[39];

    get_imw(imw);

    for (int i=0; i<(*np); i++) {
        YOW[i] = 0.0;
    }

    for (int n=0; n<39; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] = y[n*(*np)+i] * imw[n];
            YOW[i] += x[n*(*np)+i];
        }
    }

    for (int i=0; i<(*np); i++) {
        YOW[i] = 1.0/YOW[i];
    }

    for (int n=0; n<39; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] *=  YOW[i];
        }
    }
}

/*Returns enthalpy in mass units (Eq 27.) */
void VCKHMS(int *  np, amrex::Real *  T,  amrex::Real *  hms)
{
    amrex::Real tc[5], h[39];
    amrex::Real imw[39];

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
        hms[21*(*np)+i] = h[21];
        hms[22*(*np)+i] = h[22];
        hms[23*(*np)+i] = h[23];
        hms[24*(*np)+i] = h[24];
        hms[25*(*np)+i] = h[25];
        hms[26*(*np)+i] = h[26];
        hms[27*(*np)+i] = h[27];
        hms[28*(*np)+i] = h[28];
        hms[29*(*np)+i] = h[29];
        hms[30*(*np)+i] = h[30];
        hms[31*(*np)+i] = h[31];
        hms[32*(*np)+i] = h[32];
        hms[33*(*np)+i] = h[33];
        hms[34*(*np)+i] = h[34];
        hms[35*(*np)+i] = h[35];
        hms[36*(*np)+i] = h[36];
        hms[37*(*np)+i] = h[37];
        hms[38*(*np)+i] = h[38];
    }

    for (int n=0; n<39; n++) {
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
    amrex::Real c[39*(*np)]; /*temporary storage */
    amrex::Real imw[39];

    get_imw(imw);

    /*See Eq 8 with an extra 1e6 so c goes to SI */
    for (int n=0; n<39; n++) {
        for (int i=0; i<(*np); i++) {
            c[n*(*np)+i] = 1.0e6 * rho[i] * y[n*(*np)+i] * imw[n];
        }
    }

    /*call productionRate */
    vproductionRate(*np, wdot, c, T);

    /*convert to chemkin units */
    for (int i=0; i<39*(*np); i++) {
        wdot[i] *= 1.0e-6;
    }
}


/*Compute P = rhoRT/W(y) */
void VCKPY(int *  np, amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  P)
{
    amrex::Real YOW[*np];
    amrex::Real imw[39];

    get_imw(imw);

    for (int i=0; i<(*np); i++) {
        YOW[i] = 0.0;
    }

    for (int n=0; n<39; n++) {
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
    amrex::Real k_f_s[175*npt], Kc_s[175*npt], mixture[npt], g_RT[39*npt];
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

    for (int n=0; n<39; n++) {
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
    vcomp_wdot_101_150(npt, wdot, mixture, sc, k_f_s, Kc_s, tc, invT, T);
    vcomp_wdot_151_175(npt, wdot, mixture, sc, k_f_s, Kc_s, tc, invT, T);
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
        k_f_s[124*npt+i] = prefactor_units[124] * fwd_A[124] * exp(fwd_beta[124] * tc[i] - activation_units[124] * fwd_Ea[124] * invT[i]);
        k_f_s[125*npt+i] = prefactor_units[125] * fwd_A[125] * exp(fwd_beta[125] * tc[i] - activation_units[125] * fwd_Ea[125] * invT[i]);
        k_f_s[126*npt+i] = prefactor_units[126] * fwd_A[126] * exp(fwd_beta[126] * tc[i] - activation_units[126] * fwd_Ea[126] * invT[i]);
        k_f_s[127*npt+i] = prefactor_units[127] * fwd_A[127] * exp(fwd_beta[127] * tc[i] - activation_units[127] * fwd_Ea[127] * invT[i]);
        k_f_s[128*npt+i] = prefactor_units[128] * fwd_A[128] * exp(fwd_beta[128] * tc[i] - activation_units[128] * fwd_Ea[128] * invT[i]);
        k_f_s[129*npt+i] = prefactor_units[129] * fwd_A[129] * exp(fwd_beta[129] * tc[i] - activation_units[129] * fwd_Ea[129] * invT[i]);
        k_f_s[130*npt+i] = prefactor_units[130] * fwd_A[130] * exp(fwd_beta[130] * tc[i] - activation_units[130] * fwd_Ea[130] * invT[i]);
        k_f_s[131*npt+i] = prefactor_units[131] * fwd_A[131] * exp(fwd_beta[131] * tc[i] - activation_units[131] * fwd_Ea[131] * invT[i]);
        k_f_s[132*npt+i] = prefactor_units[132] * fwd_A[132] * exp(fwd_beta[132] * tc[i] - activation_units[132] * fwd_Ea[132] * invT[i]);
        k_f_s[133*npt+i] = prefactor_units[133] * fwd_A[133] * exp(fwd_beta[133] * tc[i] - activation_units[133] * fwd_Ea[133] * invT[i]);
        k_f_s[134*npt+i] = prefactor_units[134] * fwd_A[134] * exp(fwd_beta[134] * tc[i] - activation_units[134] * fwd_Ea[134] * invT[i]);
        k_f_s[135*npt+i] = prefactor_units[135] * fwd_A[135] * exp(fwd_beta[135] * tc[i] - activation_units[135] * fwd_Ea[135] * invT[i]);
        k_f_s[136*npt+i] = prefactor_units[136] * fwd_A[136] * exp(fwd_beta[136] * tc[i] - activation_units[136] * fwd_Ea[136] * invT[i]);
        k_f_s[137*npt+i] = prefactor_units[137] * fwd_A[137] * exp(fwd_beta[137] * tc[i] - activation_units[137] * fwd_Ea[137] * invT[i]);
        k_f_s[138*npt+i] = prefactor_units[138] * fwd_A[138] * exp(fwd_beta[138] * tc[i] - activation_units[138] * fwd_Ea[138] * invT[i]);
        k_f_s[139*npt+i] = prefactor_units[139] * fwd_A[139] * exp(fwd_beta[139] * tc[i] - activation_units[139] * fwd_Ea[139] * invT[i]);
        k_f_s[140*npt+i] = prefactor_units[140] * fwd_A[140] * exp(fwd_beta[140] * tc[i] - activation_units[140] * fwd_Ea[140] * invT[i]);
        k_f_s[141*npt+i] = prefactor_units[141] * fwd_A[141] * exp(fwd_beta[141] * tc[i] - activation_units[141] * fwd_Ea[141] * invT[i]);
        k_f_s[142*npt+i] = prefactor_units[142] * fwd_A[142] * exp(fwd_beta[142] * tc[i] - activation_units[142] * fwd_Ea[142] * invT[i]);
        k_f_s[143*npt+i] = prefactor_units[143] * fwd_A[143] * exp(fwd_beta[143] * tc[i] - activation_units[143] * fwd_Ea[143] * invT[i]);
        k_f_s[144*npt+i] = prefactor_units[144] * fwd_A[144] * exp(fwd_beta[144] * tc[i] - activation_units[144] * fwd_Ea[144] * invT[i]);
        k_f_s[145*npt+i] = prefactor_units[145] * fwd_A[145] * exp(fwd_beta[145] * tc[i] - activation_units[145] * fwd_Ea[145] * invT[i]);
        k_f_s[146*npt+i] = prefactor_units[146] * fwd_A[146] * exp(fwd_beta[146] * tc[i] - activation_units[146] * fwd_Ea[146] * invT[i]);
        k_f_s[147*npt+i] = prefactor_units[147] * fwd_A[147] * exp(fwd_beta[147] * tc[i] - activation_units[147] * fwd_Ea[147] * invT[i]);
        k_f_s[148*npt+i] = prefactor_units[148] * fwd_A[148] * exp(fwd_beta[148] * tc[i] - activation_units[148] * fwd_Ea[148] * invT[i]);
        k_f_s[149*npt+i] = prefactor_units[149] * fwd_A[149] * exp(fwd_beta[149] * tc[i] - activation_units[149] * fwd_Ea[149] * invT[i]);
        k_f_s[150*npt+i] = prefactor_units[150] * fwd_A[150] * exp(fwd_beta[150] * tc[i] - activation_units[150] * fwd_Ea[150] * invT[i]);
        k_f_s[151*npt+i] = prefactor_units[151] * fwd_A[151] * exp(fwd_beta[151] * tc[i] - activation_units[151] * fwd_Ea[151] * invT[i]);
        k_f_s[152*npt+i] = prefactor_units[152] * fwd_A[152] * exp(fwd_beta[152] * tc[i] - activation_units[152] * fwd_Ea[152] * invT[i]);
        k_f_s[153*npt+i] = prefactor_units[153] * fwd_A[153] * exp(fwd_beta[153] * tc[i] - activation_units[153] * fwd_Ea[153] * invT[i]);
        k_f_s[154*npt+i] = prefactor_units[154] * fwd_A[154] * exp(fwd_beta[154] * tc[i] - activation_units[154] * fwd_Ea[154] * invT[i]);
        k_f_s[155*npt+i] = prefactor_units[155] * fwd_A[155] * exp(fwd_beta[155] * tc[i] - activation_units[155] * fwd_Ea[155] * invT[i]);
        k_f_s[156*npt+i] = prefactor_units[156] * fwd_A[156] * exp(fwd_beta[156] * tc[i] - activation_units[156] * fwd_Ea[156] * invT[i]);
        k_f_s[157*npt+i] = prefactor_units[157] * fwd_A[157] * exp(fwd_beta[157] * tc[i] - activation_units[157] * fwd_Ea[157] * invT[i]);
        k_f_s[158*npt+i] = prefactor_units[158] * fwd_A[158] * exp(fwd_beta[158] * tc[i] - activation_units[158] * fwd_Ea[158] * invT[i]);
        k_f_s[159*npt+i] = prefactor_units[159] * fwd_A[159] * exp(fwd_beta[159] * tc[i] - activation_units[159] * fwd_Ea[159] * invT[i]);
        k_f_s[160*npt+i] = prefactor_units[160] * fwd_A[160] * exp(fwd_beta[160] * tc[i] - activation_units[160] * fwd_Ea[160] * invT[i]);
        k_f_s[161*npt+i] = prefactor_units[161] * fwd_A[161] * exp(fwd_beta[161] * tc[i] - activation_units[161] * fwd_Ea[161] * invT[i]);
        k_f_s[162*npt+i] = prefactor_units[162] * fwd_A[162] * exp(fwd_beta[162] * tc[i] - activation_units[162] * fwd_Ea[162] * invT[i]);
        k_f_s[163*npt+i] = prefactor_units[163] * fwd_A[163] * exp(fwd_beta[163] * tc[i] - activation_units[163] * fwd_Ea[163] * invT[i]);
        k_f_s[164*npt+i] = prefactor_units[164] * fwd_A[164] * exp(fwd_beta[164] * tc[i] - activation_units[164] * fwd_Ea[164] * invT[i]);
        k_f_s[165*npt+i] = prefactor_units[165] * fwd_A[165] * exp(fwd_beta[165] * tc[i] - activation_units[165] * fwd_Ea[165] * invT[i]);
        k_f_s[166*npt+i] = prefactor_units[166] * fwd_A[166] * exp(fwd_beta[166] * tc[i] - activation_units[166] * fwd_Ea[166] * invT[i]);
        k_f_s[167*npt+i] = prefactor_units[167] * fwd_A[167] * exp(fwd_beta[167] * tc[i] - activation_units[167] * fwd_Ea[167] * invT[i]);
        k_f_s[168*npt+i] = prefactor_units[168] * fwd_A[168] * exp(fwd_beta[168] * tc[i] - activation_units[168] * fwd_Ea[168] * invT[i]);
        k_f_s[169*npt+i] = prefactor_units[169] * fwd_A[169] * exp(fwd_beta[169] * tc[i] - activation_units[169] * fwd_Ea[169] * invT[i]);
        k_f_s[170*npt+i] = prefactor_units[170] * fwd_A[170] * exp(fwd_beta[170] * tc[i] - activation_units[170] * fwd_Ea[170] * invT[i]);
        k_f_s[171*npt+i] = prefactor_units[171] * fwd_A[171] * exp(fwd_beta[171] * tc[i] - activation_units[171] * fwd_Ea[171] * invT[i]);
        k_f_s[172*npt+i] = prefactor_units[172] * fwd_A[172] * exp(fwd_beta[172] * tc[i] - activation_units[172] * fwd_Ea[172] * invT[i]);
        k_f_s[173*npt+i] = prefactor_units[173] * fwd_A[173] * exp(fwd_beta[173] * tc[i] - activation_units[173] * fwd_Ea[173] * invT[i]);
        k_f_s[174*npt+i] = prefactor_units[174] * fwd_A[174] * exp(fwd_beta[174] * tc[i] - activation_units[174] * fwd_Ea[174] * invT[i]);
    }
}

void vcomp_gibbs(int npt, amrex::Real *  g_RT, amrex::Real *  tc)
{
    /*compute the Gibbs free energy */
    for (int i=0; i<npt; i++) {
        amrex::Real tg[5], g[39];
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
        g_RT[21*npt+i] = g[21];
        g_RT[22*npt+i] = g[22];
        g_RT[23*npt+i] = g[23];
        g_RT[24*npt+i] = g[24];
        g_RT[25*npt+i] = g[25];
        g_RT[26*npt+i] = g[26];
        g_RT[27*npt+i] = g[27];
        g_RT[28*npt+i] = g[28];
        g_RT[29*npt+i] = g[29];
        g_RT[30*npt+i] = g[30];
        g_RT[31*npt+i] = g[31];
        g_RT[32*npt+i] = g[32];
        g_RT[33*npt+i] = g[33];
        g_RT[34*npt+i] = g[34];
        g_RT[35*npt+i] = g[35];
        g_RT[36*npt+i] = g[36];
        g_RT[37*npt+i] = g[37];
        g_RT[38*npt+i] = g[38];
    }
}

void vcomp_Kc(int npt, amrex::Real *  Kc_s, amrex::Real *  g_RT, amrex::Real *  invT)
{
    for (int i=0; i<npt; i++) {
        /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
        amrex::Real refC = (101325. / 8.31451) * invT[i];
        amrex::Real refCinv = 1.0 / refC;

        Kc_s[0*npt+i] = refCinv * exp((g_RT[0*npt+i] + g_RT[19*npt+i]) - (g_RT[20*npt+i]));
        Kc_s[1*npt+i] = refC * exp((g_RT[21*npt+i]) - (g_RT[7*npt+i] + g_RT[7*npt+i]));
        Kc_s[2*npt+i] = refCinv * exp((g_RT[4*npt+i] + g_RT[4*npt+i]) - (g_RT[16*npt+i]));
        Kc_s[3*npt+i] = refCinv * exp((g_RT[0*npt+i] + g_RT[4*npt+i]) - (g_RT[6*npt+i]));
        Kc_s[4*npt+i] = refCinv * exp((g_RT[0*npt+i] + g_RT[14*npt+i]) - (g_RT[16*npt+i]));
        Kc_s[5*npt+i] = refC * exp((g_RT[12*npt+i]) - (g_RT[1*npt+i] + g_RT[9*npt+i]));
        Kc_s[6*npt+i] = refCinv * exp((g_RT[0*npt+i] + g_RT[12*npt+i]) - (g_RT[14*npt+i]));
        Kc_s[7*npt+i] = refCinv * exp((g_RT[0*npt+i] + g_RT[10*npt+i]) - (g_RT[12*npt+i]));
        Kc_s[8*npt+i] = refCinv * exp((g_RT[0*npt+i] + g_RT[9*npt+i]) - (g_RT[10*npt+i]));
        Kc_s[9*npt+i] = refCinv * exp((g_RT[0*npt+i] + g_RT[2*npt+i]) - (g_RT[4*npt+i]));
        Kc_s[10*npt+i] = refCinv * exp((g_RT[5*npt+i] + g_RT[11*npt+i]) - (g_RT[22*npt+i]));
        Kc_s[11*npt+i] = refC * exp((g_RT[1*npt+i]) - (g_RT[0*npt+i] + g_RT[0*npt+i]));
        Kc_s[12*npt+i] = refCinv * exp((g_RT[5*npt+i] + g_RT[5*npt+i]) - (g_RT[19*npt+i]));
        Kc_s[13*npt+i] = refCinv * exp((g_RT[0*npt+i] + g_RT[5*npt+i]) - (g_RT[7*npt+i]));
        Kc_s[14*npt+i] = refCinv * exp((g_RT[0*npt+i] + g_RT[7*npt+i]) - (g_RT[8*npt+i]));
        Kc_s[15*npt+i] = refC * exp((g_RT[13*npt+i]) - (g_RT[0*npt+i] + g_RT[11*npt+i]));
        Kc_s[16*npt+i] = refC * exp((g_RT[15*npt+i]) - (g_RT[0*npt+i] + g_RT[13*npt+i]));
        Kc_s[17*npt+i] = refC * exp((g_RT[15*npt+i]) - (g_RT[1*npt+i] + g_RT[11*npt+i]));
        Kc_s[18*npt+i] = refC * exp((g_RT[17*npt+i]) - (g_RT[0*npt+i] + g_RT[15*npt+i]));
        Kc_s[19*npt+i] = refC * exp((g_RT[18*npt+i]) - (g_RT[0*npt+i] + g_RT[15*npt+i]));
        Kc_s[20*npt+i] = exp((g_RT[3*npt+i]) - (g_RT[2*npt+i]));
        Kc_s[21*npt+i] = refC * exp((g_RT[25*npt+i]) - (g_RT[8*npt+i] + g_RT[11*npt+i]));
        Kc_s[22*npt+i] = refC * exp((g_RT[25*npt+i]) - (g_RT[1*npt+i] + g_RT[22*npt+i]));
        Kc_s[23*npt+i] = exp((g_RT[0*npt+i] + g_RT[19*npt+i]) - (g_RT[5*npt+i] + g_RT[7*npt+i]));
        Kc_s[24*npt+i] = exp((g_RT[1*npt+i] + g_RT[5*npt+i]) - (g_RT[0*npt+i] + g_RT[7*npt+i]));
        Kc_s[25*npt+i] = exp((g_RT[1*npt+i] + g_RT[7*npt+i]) - (g_RT[0*npt+i] + g_RT[8*npt+i]));
        Kc_s[26*npt+i] = exp((g_RT[5*npt+i] + g_RT[8*npt+i]) - (g_RT[7*npt+i] + g_RT[7*npt+i]));
        Kc_s[27*npt+i] = exp((g_RT[0*npt+i] + g_RT[20*npt+i]) - (g_RT[1*npt+i] + g_RT[19*npt+i]));
        Kc_s[28*npt+i] = exp((g_RT[0*npt+i] + g_RT[20*npt+i]) - (g_RT[7*npt+i] + g_RT[7*npt+i]));
        Kc_s[29*npt+i] = exp((g_RT[5*npt+i] + g_RT[20*npt+i]) - (g_RT[7*npt+i] + g_RT[19*npt+i]));
        Kc_s[30*npt+i] = exp((g_RT[7*npt+i] + g_RT[20*npt+i]) - (g_RT[8*npt+i] + g_RT[19*npt+i]));
        Kc_s[31*npt+i] = exp((g_RT[20*npt+i] + g_RT[20*npt+i]) - (g_RT[19*npt+i] + g_RT[21*npt+i]));
        Kc_s[32*npt+i] = exp((g_RT[20*npt+i] + g_RT[20*npt+i]) - (g_RT[19*npt+i] + g_RT[21*npt+i]));
        Kc_s[33*npt+i] = exp((g_RT[0*npt+i] + g_RT[21*npt+i]) - (g_RT[7*npt+i] + g_RT[8*npt+i]));
        Kc_s[34*npt+i] = exp((g_RT[0*npt+i] + g_RT[21*npt+i]) - (g_RT[1*npt+i] + g_RT[20*npt+i]));
        Kc_s[35*npt+i] = exp((g_RT[5*npt+i] + g_RT[21*npt+i]) - (g_RT[7*npt+i] + g_RT[20*npt+i]));
        Kc_s[36*npt+i] = exp((g_RT[7*npt+i] + g_RT[21*npt+i]) - (g_RT[8*npt+i] + g_RT[20*npt+i]));
        Kc_s[37*npt+i] = exp((g_RT[7*npt+i] + g_RT[21*npt+i]) - (g_RT[8*npt+i] + g_RT[20*npt+i]));
        Kc_s[38*npt+i] = exp((g_RT[11*npt+i] + g_RT[19*npt+i]) - (g_RT[5*npt+i] + g_RT[22*npt+i]));
        Kc_s[39*npt+i] = exp((g_RT[11*npt+i] + g_RT[20*npt+i]) - (g_RT[7*npt+i] + g_RT[22*npt+i]));
        Kc_s[40*npt+i] = exp((g_RT[7*npt+i] + g_RT[11*npt+i]) - (g_RT[0*npt+i] + g_RT[22*npt+i]));
        Kc_s[41*npt+i] = exp((g_RT[13*npt+i] + g_RT[19*npt+i]) - (g_RT[11*npt+i] + g_RT[20*npt+i]));
        Kc_s[42*npt+i] = exp((g_RT[0*npt+i] + g_RT[13*npt+i]) - (g_RT[1*npt+i] + g_RT[11*npt+i]));
        Kc_s[43*npt+i] = exp((g_RT[5*npt+i] + g_RT[13*npt+i]) - (g_RT[7*npt+i] + g_RT[11*npt+i]));
        Kc_s[44*npt+i] = exp((g_RT[7*npt+i] + g_RT[13*npt+i]) - (g_RT[8*npt+i] + g_RT[11*npt+i]));
        Kc_s[45*npt+i] = exp((g_RT[5*npt+i] + g_RT[13*npt+i]) - (g_RT[0*npt+i] + g_RT[22*npt+i]));
        Kc_s[46*npt+i] = refC * exp((g_RT[13*npt+i] + g_RT[20*npt+i]) - (g_RT[0*npt+i] + g_RT[7*npt+i] + g_RT[22*npt+i]));
        Kc_s[47*npt+i] = refC * exp((g_RT[13*npt+i] + g_RT[13*npt+i]) - (g_RT[1*npt+i] + g_RT[11*npt+i] + g_RT[11*npt+i]));
        Kc_s[48*npt+i] = exp((g_RT[4*npt+i] + g_RT[13*npt+i]) - (g_RT[6*npt+i] + g_RT[11*npt+i]));
        Kc_s[49*npt+i] = exp((g_RT[13*npt+i] + g_RT[13*npt+i]) - (g_RT[11*npt+i] + g_RT[15*npt+i]));
        Kc_s[50*npt+i] = exp((g_RT[0*npt+i] + g_RT[15*npt+i]) - (g_RT[1*npt+i] + g_RT[13*npt+i]));
        Kc_s[51*npt+i] = exp((g_RT[5*npt+i] + g_RT[15*npt+i]) - (g_RT[7*npt+i] + g_RT[13*npt+i]));
        Kc_s[52*npt+i] = exp((g_RT[7*npt+i] + g_RT[15*npt+i]) - (g_RT[8*npt+i] + g_RT[13*npt+i]));
        Kc_s[53*npt+i] = exp((g_RT[15*npt+i] + g_RT[19*npt+i]) - (g_RT[13*npt+i] + g_RT[20*npt+i]));
        Kc_s[54*npt+i] = exp((g_RT[15*npt+i] + g_RT[20*npt+i]) - (g_RT[13*npt+i] + g_RT[21*npt+i]));
        Kc_s[55*npt+i] = exp((g_RT[4*npt+i] + g_RT[15*npt+i]) - (g_RT[6*npt+i] + g_RT[13*npt+i]));
        Kc_s[56*npt+i] = exp((g_RT[4*npt+i] + g_RT[5*npt+i]) - (g_RT[0*npt+i] + g_RT[15*npt+i]));
        Kc_s[57*npt+i] = exp((g_RT[4*npt+i] + g_RT[19*npt+i]) - (g_RT[5*npt+i] + g_RT[18*npt+i]));
        Kc_s[58*npt+i] = exp((g_RT[4*npt+i] + g_RT[19*npt+i]) - (g_RT[7*npt+i] + g_RT[15*npt+i]));
        Kc_s[59*npt+i] = exp((g_RT[4*npt+i] + g_RT[20*npt+i]) - (g_RT[7*npt+i] + g_RT[18*npt+i]));
        Kc_s[60*npt+i] = exp((g_RT[0*npt+i] + g_RT[6*npt+i]) - (g_RT[1*npt+i] + g_RT[4*npt+i]));
        Kc_s[61*npt+i] = exp((g_RT[5*npt+i] + g_RT[6*npt+i]) - (g_RT[4*npt+i] + g_RT[7*npt+i]));
        Kc_s[62*npt+i] = exp((g_RT[6*npt+i] + g_RT[7*npt+i]) - (g_RT[4*npt+i] + g_RT[8*npt+i]));
        Kc_s[63*npt+i] = exp((g_RT[4*npt+i] + g_RT[20*npt+i]) - (g_RT[6*npt+i] + g_RT[19*npt+i]));
        Kc_s[64*npt+i] = exp((g_RT[6*npt+i] + g_RT[20*npt+i]) - (g_RT[4*npt+i] + g_RT[21*npt+i]));
        Kc_s[65*npt+i] = exp((g_RT[0*npt+i] + g_RT[17*npt+i]) - (g_RT[1*npt+i] + g_RT[15*npt+i]));
        Kc_s[66*npt+i] = exp((g_RT[0*npt+i] + g_RT[17*npt+i]) - (g_RT[4*npt+i] + g_RT[7*npt+i]));
        Kc_s[67*npt+i] = exp((g_RT[5*npt+i] + g_RT[17*npt+i]) - (g_RT[7*npt+i] + g_RT[15*npt+i]));
        Kc_s[68*npt+i] = exp((g_RT[7*npt+i] + g_RT[17*npt+i]) - (g_RT[8*npt+i] + g_RT[15*npt+i]));
        Kc_s[69*npt+i] = exp((g_RT[17*npt+i] + g_RT[19*npt+i]) - (g_RT[15*npt+i] + g_RT[20*npt+i]));
        Kc_s[70*npt+i] = exp((g_RT[17*npt+i] + g_RT[19*npt+i]) - (g_RT[15*npt+i] + g_RT[20*npt+i]));
        Kc_s[71*npt+i] = exp((g_RT[17*npt+i] + g_RT[20*npt+i]) - (g_RT[15*npt+i] + g_RT[21*npt+i]));
        Kc_s[72*npt+i] = exp((g_RT[13*npt+i] + g_RT[17*npt+i]) - (g_RT[15*npt+i] + g_RT[15*npt+i]));
        Kc_s[73*npt+i] = exp((g_RT[0*npt+i] + g_RT[18*npt+i]) - (g_RT[4*npt+i] + g_RT[7*npt+i]));
        Kc_s[74*npt+i] = exp((g_RT[5*npt+i] + g_RT[18*npt+i]) - (g_RT[7*npt+i] + g_RT[15*npt+i]));
        Kc_s[75*npt+i] = exp((g_RT[7*npt+i] + g_RT[18*npt+i]) - (g_RT[8*npt+i] + g_RT[15*npt+i]));
        Kc_s[76*npt+i] = exp((g_RT[18*npt+i] + g_RT[19*npt+i]) - (g_RT[15*npt+i] + g_RT[20*npt+i]));
        Kc_s[77*npt+i] = exp((g_RT[18*npt+i] + g_RT[19*npt+i]) - (g_RT[15*npt+i] + g_RT[20*npt+i]));
        Kc_s[78*npt+i] = exp((g_RT[18*npt+i] + g_RT[20*npt+i]) - (g_RT[15*npt+i] + g_RT[21*npt+i]));
        Kc_s[79*npt+i] = exp((g_RT[11*npt+i] + g_RT[18*npt+i]) - (g_RT[4*npt+i] + g_RT[22*npt+i]));
        Kc_s[80*npt+i] = exp((g_RT[4*npt+i] + g_RT[4*npt+i]) - (g_RT[0*npt+i] + g_RT[14*npt+i]));
        Kc_s[81*npt+i] = exp((g_RT[2*npt+i] + g_RT[6*npt+i]) - (g_RT[4*npt+i] + g_RT[4*npt+i]));
        Kc_s[82*npt+i] = exp((g_RT[3*npt+i] + g_RT[6*npt+i]) - (g_RT[4*npt+i] + g_RT[4*npt+i]));
        Kc_s[83*npt+i] = exp((g_RT[4*npt+i] + g_RT[7*npt+i]) - (g_RT[2*npt+i] + g_RT[8*npt+i]));
        Kc_s[84*npt+i] = exp((g_RT[4*npt+i] + g_RT[7*npt+i]) - (g_RT[3*npt+i] + g_RT[8*npt+i]));
        Kc_s[85*npt+i] = exp((g_RT[2*npt+i] + g_RT[4*npt+i]) - (g_RT[0*npt+i] + g_RT[12*npt+i]));
        Kc_s[86*npt+i] = exp((g_RT[3*npt+i] + g_RT[4*npt+i]) - (g_RT[0*npt+i] + g_RT[12*npt+i]));
        Kc_s[87*npt+i] = exp((g_RT[0*npt+i] + g_RT[18*npt+i]) - (g_RT[3*npt+i] + g_RT[8*npt+i]));
        Kc_s[88*npt+i] = exp((g_RT[0*npt+i] + g_RT[16*npt+i]) - (g_RT[1*npt+i] + g_RT[14*npt+i]));
        Kc_s[89*npt+i] = exp((g_RT[5*npt+i] + g_RT[16*npt+i]) - (g_RT[7*npt+i] + g_RT[14*npt+i]));
        Kc_s[90*npt+i] = exp((g_RT[7*npt+i] + g_RT[16*npt+i]) - (g_RT[8*npt+i] + g_RT[14*npt+i]));
        Kc_s[91*npt+i] = exp((g_RT[16*npt+i] + g_RT[19*npt+i]) - (g_RT[14*npt+i] + g_RT[20*npt+i]));
        Kc_s[92*npt+i] = exp((g_RT[16*npt+i] + g_RT[20*npt+i]) - (g_RT[14*npt+i] + g_RT[21*npt+i]));
        Kc_s[93*npt+i] = exp((g_RT[4*npt+i] + g_RT[16*npt+i]) - (g_RT[6*npt+i] + g_RT[14*npt+i]));
        Kc_s[94*npt+i] = exp((g_RT[0*npt+i] + g_RT[14*npt+i]) - (g_RT[1*npt+i] + g_RT[12*npt+i]));
        Kc_s[95*npt+i] = exp((g_RT[5*npt+i] + g_RT[14*npt+i]) - (g_RT[4*npt+i] + g_RT[15*npt+i]));
        Kc_s[96*npt+i] = exp((g_RT[14*npt+i] + g_RT[19*npt+i]) - (g_RT[12*npt+i] + g_RT[20*npt+i]));
        Kc_s[97*npt+i] = exp((g_RT[14*npt+i] + g_RT[14*npt+i]) - (g_RT[12*npt+i] + g_RT[16*npt+i]));
        Kc_s[98*npt+i] = exp((g_RT[13*npt+i] + g_RT[14*npt+i]) - (g_RT[11*npt+i] + g_RT[16*npt+i]));
        Kc_s[99*npt+i] = exp((g_RT[5*npt+i] + g_RT[14*npt+i]) - (g_RT[0*npt+i] + g_RT[23*npt+i]));
        Kc_s[100*npt+i] = exp((g_RT[0*npt+i] + g_RT[12*npt+i]) - (g_RT[1*npt+i] + g_RT[10*npt+i]));
        Kc_s[101*npt+i] = exp((g_RT[7*npt+i] + g_RT[12*npt+i]) - (g_RT[8*npt+i] + g_RT[10*npt+i]));
        Kc_s[102*npt+i] = exp((g_RT[4*npt+i] + g_RT[12*npt+i]) - (g_RT[6*npt+i] + g_RT[10*npt+i]));
        Kc_s[103*npt+i] = exp((g_RT[5*npt+i] + g_RT[12*npt+i]) - (g_RT[4*npt+i] + g_RT[13*npt+i]));
        Kc_s[104*npt+i] = exp((g_RT[7*npt+i] + g_RT[10*npt+i]) - (g_RT[8*npt+i] + g_RT[9*npt+i]));
        Kc_s[105*npt+i] = exp((g_RT[5*npt+i] + g_RT[12*npt+i]) - (g_RT[7*npt+i] + g_RT[10*npt+i]));
        Kc_s[106*npt+i] = exp((g_RT[12*npt+i] + g_RT[19*npt+i]) - (g_RT[10*npt+i] + g_RT[20*npt+i]));
        Kc_s[107*npt+i] = exp((g_RT[0*npt+i] + g_RT[10*npt+i]) - (g_RT[1*npt+i] + g_RT[9*npt+i]));
        Kc_s[108*npt+i] = exp((g_RT[10*npt+i] + g_RT[21*npt+i]) - (g_RT[12*npt+i] + g_RT[20*npt+i]));
        Kc_s[109*npt+i] = exp((g_RT[4*npt+i] + g_RT[10*npt+i]) - (g_RT[6*npt+i] + g_RT[9*npt+i]));
        Kc_s[110*npt+i] = exp((g_RT[10*npt+i] + g_RT[10*npt+i]) - (g_RT[9*npt+i] + g_RT[12*npt+i]));
        Kc_s[111*npt+i] = exp((g_RT[10*npt+i] + g_RT[19*npt+i]) - (g_RT[13*npt+i] + g_RT[15*npt+i]));
        Kc_s[112*npt+i] = exp((g_RT[10*npt+i] + g_RT[19*npt+i]) - (g_RT[9*npt+i] + g_RT[20*npt+i]));
        Kc_s[113*npt+i] = exp((g_RT[5*npt+i] + g_RT[9*npt+i]) - (g_RT[2*npt+i] + g_RT[11*npt+i]));
        Kc_s[114*npt+i] = exp((g_RT[7*npt+i] + g_RT[9*npt+i]) - (g_RT[4*npt+i] + g_RT[11*npt+i]));
        Kc_s[115*npt+i] = exp((g_RT[2*npt+i] + g_RT[5*npt+i]) - (g_RT[0*npt+i] + g_RT[13*npt+i]));
        Kc_s[116*npt+i] = exp((g_RT[2*npt+i] + g_RT[7*npt+i]) - (g_RT[0*npt+i] + g_RT[15*npt+i]));
        Kc_s[117*npt+i] = exp((g_RT[1*npt+i] + g_RT[2*npt+i]) - (g_RT[0*npt+i] + g_RT[4*npt+i]));
        Kc_s[118*npt+i] = exp((g_RT[2*npt+i] + g_RT[19*npt+i]) - (g_RT[7*npt+i] + g_RT[13*npt+i]));
        Kc_s[119*npt+i] = exp((g_RT[2*npt+i] + g_RT[20*npt+i]) - (g_RT[7*npt+i] + g_RT[15*npt+i]));
        Kc_s[120*npt+i] = exp((g_RT[2*npt+i] + g_RT[2*npt+i]) - (g_RT[1*npt+i] + g_RT[9*npt+i]));
        Kc_s[121*npt+i] = exp((g_RT[3*npt+i] + g_RT[8*npt+i]) - (g_RT[2*npt+i] + g_RT[8*npt+i]));
        Kc_s[122*npt+i] = exp((g_RT[3*npt+i] + g_RT[11*npt+i]) - (g_RT[2*npt+i] + g_RT[11*npt+i]));
        Kc_s[123*npt+i] = exp((g_RT[3*npt+i] + g_RT[22*npt+i]) - (g_RT[2*npt+i] + g_RT[22*npt+i]));
        Kc_s[124*npt+i] = exp((g_RT[3*npt+i] + g_RT[5*npt+i]) - (g_RT[1*npt+i] + g_RT[11*npt+i]));
        Kc_s[125*npt+i] = exp((g_RT[3*npt+i] + g_RT[5*npt+i]) - (g_RT[0*npt+i] + g_RT[13*npt+i]));
        Kc_s[126*npt+i] = exp((g_RT[3*npt+i] + g_RT[7*npt+i]) - (g_RT[0*npt+i] + g_RT[15*npt+i]));
        Kc_s[127*npt+i] = exp((g_RT[1*npt+i] + g_RT[3*npt+i]) - (g_RT[0*npt+i] + g_RT[4*npt+i]));
        Kc_s[128*npt+i] = refC * exp((g_RT[3*npt+i] + g_RT[19*npt+i]) - (g_RT[0*npt+i] + g_RT[7*npt+i] + g_RT[11*npt+i]));
        Kc_s[129*npt+i] = exp((g_RT[3*npt+i] + g_RT[19*npt+i]) - (g_RT[8*npt+i] + g_RT[11*npt+i]));
        Kc_s[130*npt+i] = exp((g_RT[3*npt+i] + g_RT[22*npt+i]) - (g_RT[11*npt+i] + g_RT[15*npt+i]));
        Kc_s[131*npt+i] = refC * exp((g_RT[23*npt+i]) - (g_RT[4*npt+i] + g_RT[13*npt+i]));
        Kc_s[132*npt+i] = refC * exp((g_RT[26*npt+i]) - (g_RT[4*npt+i] + g_RT[18*npt+i]));
        Kc_s[133*npt+i] = exp((g_RT[7*npt+i] + g_RT[26*npt+i]) - (g_RT[8*npt+i] + g_RT[24*npt+i]));
        Kc_s[134*npt+i] = exp((g_RT[0*npt+i] + g_RT[26*npt+i]) - (g_RT[1*npt+i] + g_RT[24*npt+i]));
        Kc_s[135*npt+i] = exp((g_RT[4*npt+i] + g_RT[26*npt+i]) - (g_RT[6*npt+i] + g_RT[24*npt+i]));
        Kc_s[136*npt+i] = exp((g_RT[5*npt+i] + g_RT[26*npt+i]) - (g_RT[7*npt+i] + g_RT[24*npt+i]));
        Kc_s[137*npt+i] = exp((g_RT[20*npt+i] + g_RT[26*npt+i]) - (g_RT[21*npt+i] + g_RT[24*npt+i]));
        Kc_s[138*npt+i] = exp((g_RT[19*npt+i] + g_RT[26*npt+i]) - (g_RT[20*npt+i] + g_RT[24*npt+i]));
        Kc_s[139*npt+i] = refC * exp((g_RT[24*npt+i]) - (g_RT[4*npt+i] + g_RT[15*npt+i]));
        Kc_s[140*npt+i] = exp((g_RT[18*npt+i] + g_RT[24*npt+i]) - (g_RT[15*npt+i] + g_RT[26*npt+i]));
        Kc_s[141*npt+i] = exp((g_RT[15*npt+i] + g_RT[24*npt+i]) - (g_RT[13*npt+i] + g_RT[26*npt+i]));
        Kc_s[142*npt+i] = exp((g_RT[20*npt+i] + g_RT[24*npt+i]) - (g_RT[7*npt+i] + g_RT[30*npt+i]));
        Kc_s[143*npt+i] = refC * exp((g_RT[30*npt+i]) - (g_RT[0*npt+i] + g_RT[29*npt+i]));
        Kc_s[144*npt+i] = exp((g_RT[19*npt+i] + g_RT[29*npt+i]) - (g_RT[20*npt+i] + g_RT[28*npt+i]));
        Kc_s[145*npt+i] = exp((g_RT[7*npt+i] + g_RT[29*npt+i]) - (g_RT[8*npt+i] + g_RT[28*npt+i]));
        Kc_s[146*npt+i] = exp((g_RT[20*npt+i] + g_RT[29*npt+i]) - (g_RT[21*npt+i] + g_RT[28*npt+i]));
        Kc_s[147*npt+i] = exp((g_RT[5*npt+i] + g_RT[29*npt+i]) - (g_RT[7*npt+i] + g_RT[28*npt+i]));
        Kc_s[148*npt+i] = exp((g_RT[0*npt+i] + g_RT[29*npt+i]) - (g_RT[1*npt+i] + g_RT[28*npt+i]));
        Kc_s[149*npt+i] = exp((g_RT[4*npt+i] + g_RT[29*npt+i]) - (g_RT[6*npt+i] + g_RT[28*npt+i]));
        Kc_s[150*npt+i] = refC * exp((g_RT[28*npt+i]) - (g_RT[11*npt+i] + g_RT[18*npt+i]));
        Kc_s[151*npt+i] = refC * exp((g_RT[28*npt+i]) - (g_RT[4*npt+i] + g_RT[22*npt+i]));
        Kc_s[152*npt+i] = refCinv * exp((g_RT[19*npt+i] + g_RT[24*npt+i]) - (g_RT[34*npt+i]));
        Kc_s[153*npt+i] = refC * exp((g_RT[34*npt+i] + g_RT[34*npt+i]) - (g_RT[19*npt+i] + g_RT[30*npt+i] + g_RT[30*npt+i]));
        Kc_s[154*npt+i] = refC * exp((g_RT[34*npt+i] + g_RT[34*npt+i]) - (g_RT[19*npt+i] + g_RT[29*npt+i] + g_RT[31*npt+i]));
        Kc_s[155*npt+i] = refC * exp((g_RT[30*npt+i]) - (g_RT[15*npt+i] + g_RT[18*npt+i]));
        Kc_s[156*npt+i] = exp((g_RT[19*npt+i] + g_RT[30*npt+i]) - (g_RT[20*npt+i] + g_RT[29*npt+i]));
        Kc_s[157*npt+i] = exp((g_RT[34*npt+i]) - (g_RT[35*npt+i]));
        Kc_s[158*npt+i] = pow(refC,2.000000) * exp((g_RT[35*npt+i]) - (g_RT[7*npt+i] + g_RT[15*npt+i] + g_RT[15*npt+i]));
        Kc_s[159*npt+i] = refCinv * exp((g_RT[19*npt+i] + g_RT[35*npt+i]) - (g_RT[37*npt+i]));
        Kc_s[160*npt+i] = refC * exp((g_RT[37*npt+i]) - (g_RT[7*npt+i] + g_RT[36*npt+i]));
        Kc_s[161*npt+i] = refC * exp((g_RT[36*npt+i]) - (g_RT[7*npt+i] + g_RT[32*npt+i]));
        Kc_s[162*npt+i] = exp((g_RT[32*npt+i]) - (g_RT[33*npt+i]));
        Kc_s[163*npt+i] = refC * exp((g_RT[33*npt+i]) - (g_RT[11*npt+i] + g_RT[27*npt+i]));
        Kc_s[164*npt+i] = refC * exp((g_RT[33*npt+i]) - (g_RT[17*npt+i] + g_RT[22*npt+i]));
        Kc_s[165*npt+i] = refC * exp((g_RT[27*npt+i]) - (g_RT[0*npt+i] + g_RT[25*npt+i]));
        Kc_s[166*npt+i] = refCinv * exp((g_RT[7*npt+i] + g_RT[15*npt+i]) - (g_RT[27*npt+i]));
        Kc_s[167*npt+i] = refC * exp((g_RT[25*npt+i]) - (g_RT[7*npt+i] + g_RT[13*npt+i]));
        Kc_s[168*npt+i] = refC * exp((g_RT[7*npt+i] + g_RT[25*npt+i]) - (g_RT[0*npt+i] + g_RT[8*npt+i] + g_RT[22*npt+i]));
        Kc_s[169*npt+i] = refC * exp((g_RT[7*npt+i] + g_RT[25*npt+i]) - (g_RT[7*npt+i] + g_RT[8*npt+i] + g_RT[11*npt+i]));
        Kc_s[170*npt+i] = refC * exp((g_RT[0*npt+i] + g_RT[25*npt+i]) - (g_RT[0*npt+i] + g_RT[1*npt+i] + g_RT[22*npt+i]));
        Kc_s[171*npt+i] = refC * exp((g_RT[0*npt+i] + g_RT[25*npt+i]) - (g_RT[1*npt+i] + g_RT[7*npt+i] + g_RT[11*npt+i]));
        Kc_s[172*npt+i] = refC * exp((g_RT[4*npt+i] + g_RT[25*npt+i]) - (g_RT[6*npt+i] + g_RT[7*npt+i] + g_RT[11*npt+i]));
        Kc_s[173*npt+i] = refC * exp((g_RT[20*npt+i] + g_RT[25*npt+i]) - (g_RT[7*npt+i] + g_RT[11*npt+i] + g_RT[21*npt+i]));
        Kc_s[174*npt+i] = refC * exp((g_RT[5*npt+i] + g_RT[25*npt+i]) - (g_RT[7*npt+i] + g_RT[7*npt+i] + g_RT[11*npt+i]));
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
        phi_f = sc[0*npt+i]*sc[19*npt+i];
        alpha = mixture[i] + (TB[0][0] - 1)*sc[1*npt+i] + (TB[0][1] - 1)*sc[8*npt+i] + (TB[0][2] - 1)*sc[19*npt+i] + (TB[0][3] - 1)*sc[11*npt+i] + (TB[0][4] - 1)*sc[22*npt+i];
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
        phi_r = sc[20*npt+i];
        Kc = Kc_s[0*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[19*npt+i] -= qdot;
        wdot[20*npt+i] += qdot;

        /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
        phi_f = sc[21*npt+i];
        alpha = mixture[i] + (TB[1][0] - 1)*sc[1*npt+i] + (TB[1][1] - 1)*sc[8*npt+i] + (TB[1][2] - 1)*sc[11*npt+i] + (TB[1][3] - 1)*sc[22*npt+i];
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
        phi_r = sc[7*npt+i]*sc[7*npt+i];
        Kc = Kc_s[1*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[7*npt+i] += qdot;
        wdot[7*npt+i] += qdot;
        wdot[21*npt+i] -= qdot;

        /*reaction 3: CH3 + CH3 (+M) <=> C2H6 (+M) */
        phi_f = sc[4*npt+i]*sc[4*npt+i];
        alpha = mixture[i] + (TB[2][0] - 1)*sc[8*npt+i] + (TB[2][1] - 1)*sc[11*npt+i] + (TB[2][2] - 1)*sc[22*npt+i];
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
        phi_r = sc[16*npt+i];
        Kc = Kc_s[2*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[4*npt+i] -= qdot;
        wdot[16*npt+i] += qdot;

        /*reaction 4: CH3 + H (+M) <=> CH4 (+M) */
        phi_f = sc[0*npt+i]*sc[4*npt+i];
        alpha = mixture[i] + (TB[3][0] - 1)*sc[1*npt+i] + (TB[3][1] - 1)*sc[8*npt+i] + (TB[3][2] - 1)*sc[6*npt+i] + (TB[3][3] - 1)*sc[11*npt+i] + (TB[3][4] - 1)*sc[22*npt+i] + (TB[3][5] - 1)*sc[16*npt+i];
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
        phi_r = sc[6*npt+i];
        Kc = Kc_s[3*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[4*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;

        /*reaction 5: C2H5 + H (+M) <=> C2H6 (+M) */
        phi_f = sc[0*npt+i]*sc[14*npt+i];
        alpha = mixture[i] + (TB[4][0] - 1)*sc[1*npt+i] + (TB[4][1] - 1)*sc[8*npt+i] + (TB[4][2] - 1)*sc[6*npt+i] + (TB[4][3] - 1)*sc[11*npt+i] + (TB[4][4] - 1)*sc[22*npt+i] + (TB[4][5] - 1)*sc[16*npt+i];
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
        wdot[0*npt+i] -= qdot;
        wdot[14*npt+i] -= qdot;
        wdot[16*npt+i] += qdot;

        /*reaction 6: C2H4 (+M) <=> H2 + C2H2 (+M) */
        phi_f = sc[12*npt+i];
        alpha = mixture[i] + (TB[5][0] - 1)*sc[1*npt+i] + (TB[5][1] - 1)*sc[8*npt+i] + (TB[5][2] - 1)*sc[6*npt+i] + (TB[5][3] - 1)*sc[11*npt+i] + (TB[5][4] - 1)*sc[22*npt+i] + (TB[5][5] - 1)*sc[16*npt+i];
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
        phi_r = sc[1*npt+i]*sc[9*npt+i];
        Kc = Kc_s[5*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[9*npt+i] += qdot;
        wdot[12*npt+i] -= qdot;

        /*reaction 7: C2H4 + H (+M) <=> C2H5 (+M) */
        phi_f = sc[0*npt+i]*sc[12*npt+i];
        alpha = mixture[i] + (TB[6][0] - 1)*sc[1*npt+i] + (TB[6][1] - 1)*sc[8*npt+i] + (TB[6][2] - 1)*sc[6*npt+i] + (TB[6][3] - 1)*sc[11*npt+i] + (TB[6][4] - 1)*sc[22*npt+i] + (TB[6][5] - 1)*sc[16*npt+i];
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
        wdot[12*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;

        /*reaction 8: C2H3 + H (+M) <=> C2H4 (+M) */
        phi_f = sc[0*npt+i]*sc[10*npt+i];
        alpha = mixture[i] + (TB[7][0] - 1)*sc[1*npt+i] + (TB[7][1] - 1)*sc[8*npt+i] + (TB[7][2] - 1)*sc[6*npt+i] + (TB[7][3] - 1)*sc[11*npt+i] + (TB[7][4] - 1)*sc[22*npt+i] + (TB[7][5] - 1)*sc[16*npt+i];
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
        phi_r = sc[12*npt+i];
        Kc = Kc_s[7*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[10*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;

        /*reaction 9: C2H2 + H (+M) <=> C2H3 (+M) */
        phi_f = sc[0*npt+i]*sc[9*npt+i];
        alpha = mixture[i] + (TB[8][0] - 1)*sc[1*npt+i] + (TB[8][1] - 1)*sc[8*npt+i] + (TB[8][2] - 1)*sc[6*npt+i] + (TB[8][3] - 1)*sc[11*npt+i] + (TB[8][4] - 1)*sc[22*npt+i] + (TB[8][5] - 1)*sc[16*npt+i];
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
        phi_r = sc[10*npt+i];
        Kc = Kc_s[8*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[9*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;

        /*reaction 10: CH2 + H (+M) <=> CH3 (+M) */
        phi_f = sc[0*npt+i]*sc[2*npt+i];
        alpha = mixture[i] + (TB[9][0] - 1)*sc[1*npt+i] + (TB[9][1] - 1)*sc[8*npt+i] + (TB[9][2] - 1)*sc[6*npt+i] + (TB[9][3] - 1)*sc[11*npt+i] + (TB[9][4] - 1)*sc[22*npt+i] + (TB[9][5] - 1)*sc[16*npt+i];
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
        phi_r = sc[4*npt+i];
        Kc = Kc_s[9*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;

        /*reaction 11: CO + O (+M) <=> CO2 (+M) */
        phi_f = sc[5*npt+i]*sc[11*npt+i];
        alpha = mixture[i] + (TB[10][0] - 1)*sc[1*npt+i] + (TB[10][1] - 1)*sc[8*npt+i] + (TB[10][2] - 1)*sc[11*npt+i] + (TB[10][3] - 1)*sc[22*npt+i];
        k_f = k_f_s[10*npt+i];
        redP = alpha / k_f * phase_units[10] * low_A[10] * exp(low_beta[10] * tc[i] - activation_units[10] * low_Ea[10] * invT[i]);
        F = redP / (1 + redP);
        k_f *= F;
        q_f = phi_f * k_f;
        phi_r = sc[22*npt+i];
        Kc = Kc_s[10*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[11*npt+i] -= qdot;
        wdot[22*npt+i] += qdot;

        /*reaction 12: H2 + M <=> H + H + M */
        phi_f = sc[1*npt+i];
        alpha = mixture[i] + (TB[11][0] - 1)*sc[1*npt+i] + (TB[11][1] - 1)*sc[8*npt+i] + (TB[11][2] - 1)*sc[11*npt+i] + (TB[11][3] - 1)*sc[22*npt+i];
        k_f = alpha * k_f_s[11*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[0*npt+i];
        Kc = Kc_s[11*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;

        /*reaction 13: O + O + M <=> O2 + M */
        phi_f = sc[5*npt+i]*sc[5*npt+i];
        alpha = mixture[i] + (TB[12][0] - 1)*sc[1*npt+i] + (TB[12][1] - 1)*sc[8*npt+i] + (TB[12][2] - 1)*sc[11*npt+i] + (TB[12][3] - 1)*sc[22*npt+i];
        k_f = alpha * k_f_s[12*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[19*npt+i];
        Kc = Kc_s[12*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[5*npt+i] -= qdot;
        wdot[19*npt+i] += qdot;

        /*reaction 14: O + H + M <=> OH + M */
        phi_f = sc[0*npt+i]*sc[5*npt+i];
        alpha = mixture[i] + (TB[13][0] - 1)*sc[1*npt+i] + (TB[13][1] - 1)*sc[8*npt+i] + (TB[13][2] - 1)*sc[11*npt+i] + (TB[13][3] - 1)*sc[22*npt+i];
        k_f = alpha * k_f_s[13*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i];
        Kc = Kc_s[13*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[5*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;

        /*reaction 15: H + OH + M <=> H2O + M */
        phi_f = sc[0*npt+i]*sc[7*npt+i];
        alpha = mixture[i] + (TB[14][0] - 1)*sc[1*npt+i] + (TB[14][1] - 1)*sc[8*npt+i] + (TB[14][2] - 1)*sc[11*npt+i] + (TB[14][3] - 1)*sc[22*npt+i];
        k_f = alpha * k_f_s[14*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i];
        Kc = Kc_s[14*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[7*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;

        /*reaction 16: HCO + M <=> H + CO + M */
        phi_f = sc[13*npt+i];
        alpha = mixture[i] + (TB[15][0] - 1)*sc[1*npt+i] + (TB[15][1] - 1)*sc[8*npt+i] + (TB[15][2] - 1)*sc[11*npt+i] + (TB[15][3] - 1)*sc[22*npt+i];
        k_f = alpha * k_f_s[15*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[11*npt+i];
        Kc = Kc_s[15*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[11*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;

        /*reaction 17: CH2O + M <=> HCO + H + M */
        phi_f = sc[15*npt+i];
        alpha = mixture[i] + (TB[16][0] - 1)*sc[1*npt+i] + (TB[16][1] - 1)*sc[8*npt+i] + (TB[16][2] - 1)*sc[11*npt+i] + (TB[16][3] - 1)*sc[22*npt+i];
        k_f = alpha * k_f_s[16*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[13*npt+i];
        Kc = Kc_s[16*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[13*npt+i] += qdot;
        wdot[15*npt+i] -= qdot;

        /*reaction 18: CH2O + M <=> CO + H2 + M */
        phi_f = sc[15*npt+i];
        alpha = mixture[i] + (TB[17][0] - 1)*sc[1*npt+i] + (TB[17][1] - 1)*sc[8*npt+i] + (TB[17][2] - 1)*sc[11*npt+i] + (TB[17][3] - 1)*sc[22*npt+i];
        k_f = alpha * k_f_s[17*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[11*npt+i];
        Kc = Kc_s[17*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[11*npt+i] += qdot;
        wdot[15*npt+i] -= qdot;

        /*reaction 19: CH2OH + M <=> CH2O + H + M */
        phi_f = sc[17*npt+i];
        alpha = mixture[i];
        k_f = alpha * k_f_s[18*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[15*npt+i];
        Kc = Kc_s[18*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[15*npt+i] += qdot;
        wdot[17*npt+i] -= qdot;

        /*reaction 20: CH3O + M <=> CH2O + H + M */
        phi_f = sc[18*npt+i];
        alpha = mixture[i];
        k_f = alpha * k_f_s[19*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[15*npt+i];
        Kc = Kc_s[19*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[15*npt+i] += qdot;
        wdot[18*npt+i] -= qdot;

        /*reaction 21: CH2(S) + M <=> CH2 + M */
        phi_f = sc[3*npt+i];
        alpha = mixture[i] + (TB[20][0] - 1)*sc[8*npt+i] + (TB[20][1] - 1)*sc[11*npt+i] + (TB[20][2] - 1)*sc[22*npt+i];
        k_f = alpha * k_f_s[20*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i];
        Kc = Kc_s[20*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;

        /*reaction 22: HCOOH + M <=> CO + H2O + M */
        phi_f = sc[25*npt+i];
        alpha = mixture[i];
        k_f = alpha * k_f_s[21*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i]*sc[11*npt+i];
        Kc = Kc_s[21*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[8*npt+i] += qdot;
        wdot[11*npt+i] += qdot;
        wdot[25*npt+i] -= qdot;

        /*reaction 23: HCOOH + M <=> CO2 + H2 + M */
        phi_f = sc[25*npt+i];
        alpha = mixture[i];
        k_f = alpha * k_f_s[22*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[22*npt+i];
        Kc = Kc_s[22*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[22*npt+i] += qdot;
        wdot[25*npt+i] -= qdot;

        /*reaction 24: H + O2 <=> O + OH */
        phi_f = sc[0*npt+i]*sc[19*npt+i];
        k_f = k_f_s[23*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[7*npt+i];
        Kc = Kc_s[23*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[7*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 25: O + H2 <=> H + OH */
        phi_f = sc[1*npt+i]*sc[5*npt+i];
        k_f = k_f_s[24*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[7*npt+i];
        Kc = Kc_s[24*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[5*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;

        /*reaction 26: H2 + OH <=> H2O + H */
        phi_f = sc[1*npt+i]*sc[7*npt+i];
        k_f = k_f_s[25*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[8*npt+i];
        Kc = Kc_s[25*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[7*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;

        /*reaction 27: O + H2O <=> OH + OH */
        phi_f = sc[5*npt+i]*sc[8*npt+i];
        k_f = k_f_s[26*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[7*npt+i];
        Kc = Kc_s[26*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[7*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;

        /*reaction 28: HO2 + H <=> H2 + O2 */
        phi_f = sc[0*npt+i]*sc[20*npt+i];
        k_f = k_f_s[27*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[19*npt+i];
        Kc = Kc_s[27*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[19*npt+i] += qdot;
        wdot[20*npt+i] -= qdot;

        /*reaction 29: HO2 + H <=> OH + OH */
        phi_f = sc[0*npt+i]*sc[20*npt+i];
        k_f = k_f_s[28*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[7*npt+i];
        Kc = Kc_s[28*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[7*npt+i] += qdot;
        wdot[20*npt+i] -= qdot;

        /*reaction 30: HO2 + O <=> O2 + OH */
        phi_f = sc[5*npt+i]*sc[20*npt+i];
        k_f = k_f_s[29*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[19*npt+i];
        Kc = Kc_s[29*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[19*npt+i] += qdot;
        wdot[20*npt+i] -= qdot;

        /*reaction 31: HO2 + OH <=> H2O + O2 */
        phi_f = sc[7*npt+i]*sc[20*npt+i];
        k_f = k_f_s[30*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i]*sc[19*npt+i];
        Kc = Kc_s[30*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[7*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[19*npt+i] += qdot;
        wdot[20*npt+i] -= qdot;

        /*reaction 32: HO2 + HO2 <=> H2O2 + O2 */
        phi_f = sc[20*npt+i]*sc[20*npt+i];
        k_f = k_f_s[31*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[19*npt+i]*sc[21*npt+i];
        Kc = Kc_s[31*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[19*npt+i] += qdot;
        wdot[20*npt+i] -= qdot;
        wdot[20*npt+i] -= qdot;
        wdot[21*npt+i] += qdot;

        /*reaction 33: HO2 + HO2 <=> H2O2 + O2 */
        phi_f = sc[20*npt+i]*sc[20*npt+i];
        k_f = k_f_s[32*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[19*npt+i]*sc[21*npt+i];
        Kc = Kc_s[32*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[19*npt+i] += qdot;
        wdot[20*npt+i] -= qdot;
        wdot[20*npt+i] -= qdot;
        wdot[21*npt+i] += qdot;

        /*reaction 34: H2O2 + H <=> H2O + OH */
        phi_f = sc[0*npt+i]*sc[21*npt+i];
        k_f = k_f_s[33*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[8*npt+i];
        Kc = Kc_s[33*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[8*npt+i] += qdot;
        wdot[21*npt+i] -= qdot;

        /*reaction 35: H2O2 + H <=> HO2 + H2 */
        phi_f = sc[0*npt+i]*sc[21*npt+i];
        k_f = k_f_s[34*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[20*npt+i];
        Kc = Kc_s[34*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[20*npt+i] += qdot;
        wdot[21*npt+i] -= qdot;

        /*reaction 36: H2O2 + O <=> OH + HO2 */
        phi_f = sc[5*npt+i]*sc[21*npt+i];
        k_f = k_f_s[35*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[20*npt+i];
        Kc = Kc_s[35*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[20*npt+i] += qdot;
        wdot[21*npt+i] -= qdot;

        /*reaction 37: H2O2 + OH <=> HO2 + H2O */
        phi_f = sc[7*npt+i]*sc[21*npt+i];
        k_f = k_f_s[36*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i]*sc[20*npt+i];
        Kc = Kc_s[36*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[7*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[20*npt+i] += qdot;
        wdot[21*npt+i] -= qdot;

        /*reaction 38: H2O2 + OH <=> HO2 + H2O */
        phi_f = sc[7*npt+i]*sc[21*npt+i];
        k_f = k_f_s[37*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i]*sc[20*npt+i];
        Kc = Kc_s[37*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[7*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[20*npt+i] += qdot;
        wdot[21*npt+i] -= qdot;

        /*reaction 39: CO + O2 <=> CO2 + O */
        phi_f = sc[11*npt+i]*sc[19*npt+i];
        k_f = k_f_s[38*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[22*npt+i];
        Kc = Kc_s[38*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;
        wdot[19*npt+i] -= qdot;
        wdot[22*npt+i] += qdot;

        /*reaction 40: CO + HO2 <=> CO2 + OH */
        phi_f = sc[11*npt+i]*sc[20*npt+i];
        k_f = k_f_s[39*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[22*npt+i];
        Kc = Kc_s[39*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[7*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;
        wdot[20*npt+i] -= qdot;
        wdot[22*npt+i] += qdot;

        /*reaction 41: CO + OH <=> CO2 + H */
        phi_f = sc[7*npt+i]*sc[11*npt+i];
        k_f = k_f_s[40*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[22*npt+i];
        Kc = Kc_s[40*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;
        wdot[11*npt+i] -= qdot;
        wdot[22*npt+i] += qdot;

        /*reaction 42: HCO + O2 <=> CO + HO2 */
        phi_f = sc[13*npt+i]*sc[19*npt+i];
        k_f = k_f_s[41*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[11*npt+i]*sc[20*npt+i];
        Kc = Kc_s[41*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[11*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;
        wdot[19*npt+i] -= qdot;
        wdot[20*npt+i] += qdot;

        /*reaction 43: HCO + H <=> CO + H2 */
        phi_f = sc[0*npt+i]*sc[13*npt+i];
        k_f = k_f_s[42*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[11*npt+i];
        Kc = Kc_s[42*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[11*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;

        /*reaction 44: HCO + O <=> CO + OH */
        phi_f = sc[5*npt+i]*sc[13*npt+i];
        k_f = k_f_s[43*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[11*npt+i];
        Kc = Kc_s[43*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[11*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;

        /*reaction 45: HCO + OH <=> CO + H2O */
        phi_f = sc[7*npt+i]*sc[13*npt+i];
        k_f = k_f_s[44*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i]*sc[11*npt+i];
        Kc = Kc_s[44*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[7*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[11*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;

        /*reaction 46: HCO + O <=> CO2 + H */
        phi_f = sc[5*npt+i]*sc[13*npt+i];
        k_f = k_f_s[45*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[22*npt+i];
        Kc = Kc_s[45*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[13*npt+i] -= qdot;
        wdot[22*npt+i] += qdot;

        /*reaction 47: HCO + HO2 <=> CO2 + OH + H */
        phi_f = sc[13*npt+i]*sc[20*npt+i];
        k_f = k_f_s[46*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[7*npt+i]*sc[22*npt+i];
        Kc = Kc_s[46*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[7*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;
        wdot[20*npt+i] -= qdot;
        wdot[22*npt+i] += qdot;

        /*reaction 48: HCO + HCO <=> H2 + CO + CO */
        phi_f = sc[13*npt+i]*sc[13*npt+i];
        k_f = k_f_s[47*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[11*npt+i]*sc[11*npt+i];
        Kc = Kc_s[47*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[11*npt+i] += qdot;
        wdot[11*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;
        wdot[13*npt+i] -= qdot;

        /*reaction 49: HCO + CH3 <=> CO + CH4 */
        phi_f = sc[4*npt+i]*sc[13*npt+i];
        k_f = k_f_s[48*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[11*npt+i];
        Kc = Kc_s[48*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[11*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;

        /*reaction 50: HCO + HCO <=> CH2O + CO */
        phi_f = sc[13*npt+i]*sc[13*npt+i];
        k_f = k_f_s[49*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[11*npt+i]*sc[15*npt+i];
        Kc = Kc_s[49*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[11*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;
        wdot[13*npt+i] -= qdot;
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

        /*reaction 51: CH2O + H <=> HCO + H2 */
        phi_f = sc[0*npt+i]*sc[15*npt+i];
        k_f = k_f_s[50*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[13*npt+i];
        Kc = Kc_s[50*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[13*npt+i] += qdot;
        wdot[15*npt+i] -= qdot;

        /*reaction 52: CH2O + O <=> HCO + OH */
        phi_f = sc[5*npt+i]*sc[15*npt+i];
        k_f = k_f_s[51*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[13*npt+i];
        Kc = Kc_s[51*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[13*npt+i] += qdot;
        wdot[15*npt+i] -= qdot;

        /*reaction 53: CH2O + OH <=> HCO + H2O */
        phi_f = sc[7*npt+i]*sc[15*npt+i];
        k_f = k_f_s[52*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i]*sc[13*npt+i];
        Kc = Kc_s[52*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[7*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[13*npt+i] += qdot;
        wdot[15*npt+i] -= qdot;

        /*reaction 54: CH2O + O2 <=> HCO + HO2 */
        phi_f = sc[15*npt+i]*sc[19*npt+i];
        k_f = k_f_s[53*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[13*npt+i]*sc[20*npt+i];
        Kc = Kc_s[53*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[13*npt+i] += qdot;
        wdot[15*npt+i] -= qdot;
        wdot[19*npt+i] -= qdot;
        wdot[20*npt+i] += qdot;

        /*reaction 55: CH2O + HO2 <=> HCO + H2O2 */
        phi_f = sc[15*npt+i]*sc[20*npt+i];
        k_f = k_f_s[54*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[13*npt+i]*sc[21*npt+i];
        Kc = Kc_s[54*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[13*npt+i] += qdot;
        wdot[15*npt+i] -= qdot;
        wdot[20*npt+i] -= qdot;
        wdot[21*npt+i] += qdot;

        /*reaction 56: CH2O + CH3 <=> HCO + CH4 */
        phi_f = sc[4*npt+i]*sc[15*npt+i];
        k_f = k_f_s[55*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[13*npt+i];
        Kc = Kc_s[55*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[13*npt+i] += qdot;
        wdot[15*npt+i] -= qdot;

        /*reaction 57: CH3 + O <=> CH2O + H */
        phi_f = sc[4*npt+i]*sc[5*npt+i];
        k_f = k_f_s[56*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[15*npt+i];
        Kc = Kc_s[56*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;

        /*reaction 58: CH3 + O2 <=> CH3O + O */
        phi_f = sc[4*npt+i]*sc[19*npt+i];
        k_f = k_f_s[57*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[18*npt+i];
        Kc = Kc_s[57*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[18*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 59: CH3 + O2 <=> CH2O + OH */
        phi_f = sc[4*npt+i]*sc[19*npt+i];
        k_f = k_f_s[58*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[15*npt+i];
        Kc = Kc_s[58*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[15*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 60: CH3 + HO2 <=> CH3O + OH */
        phi_f = sc[4*npt+i]*sc[20*npt+i];
        k_f = k_f_s[59*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[18*npt+i];
        Kc = Kc_s[59*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[18*npt+i] += qdot;
        wdot[20*npt+i] -= qdot;

        /*reaction 61: CH4 + H <=> CH3 + H2 */
        phi_f = sc[0*npt+i]*sc[6*npt+i];
        k_f = k_f_s[60*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[4*npt+i];
        Kc = Kc_s[60*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[4*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 62: CH4 + O <=> CH3 + OH */
        phi_f = sc[5*npt+i]*sc[6*npt+i];
        k_f = k_f_s[61*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[7*npt+i];
        Kc = Kc_s[61*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[6*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;

        /*reaction 63: CH4 + OH <=> CH3 + H2O */
        phi_f = sc[6*npt+i]*sc[7*npt+i];
        k_f = k_f_s[62*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[8*npt+i];
        Kc = Kc_s[62*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[7*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;

        /*reaction 64: CH3 + HO2 <=> CH4 + O2 */
        phi_f = sc[4*npt+i]*sc[20*npt+i];
        k_f = k_f_s[63*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[19*npt+i];
        Kc = Kc_s[63*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[19*npt+i] += qdot;
        wdot[20*npt+i] -= qdot;

        /*reaction 65: CH4 + HO2 <=> CH3 + H2O2 */
        phi_f = sc[6*npt+i]*sc[20*npt+i];
        k_f = k_f_s[64*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[21*npt+i];
        Kc = Kc_s[64*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[20*npt+i] -= qdot;
        wdot[21*npt+i] += qdot;

        /*reaction 66: CH2OH + H <=> CH2O + H2 */
        phi_f = sc[0*npt+i]*sc[17*npt+i];
        k_f = k_f_s[65*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[15*npt+i];
        Kc = Kc_s[65*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[15*npt+i] += qdot;
        wdot[17*npt+i] -= qdot;

        /*reaction 67: CH2OH + H <=> CH3 + OH */
        phi_f = sc[0*npt+i]*sc[17*npt+i];
        k_f = k_f_s[66*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[7*npt+i];
        Kc = Kc_s[66*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[7*npt+i] += qdot;
        wdot[17*npt+i] -= qdot;

        /*reaction 68: CH2OH + O <=> CH2O + OH */
        phi_f = sc[5*npt+i]*sc[17*npt+i];
        k_f = k_f_s[67*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[15*npt+i];
        Kc = Kc_s[67*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[15*npt+i] += qdot;
        wdot[17*npt+i] -= qdot;

        /*reaction 69: CH2OH + OH <=> CH2O + H2O */
        phi_f = sc[7*npt+i]*sc[17*npt+i];
        k_f = k_f_s[68*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i]*sc[15*npt+i];
        Kc = Kc_s[68*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[7*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[15*npt+i] += qdot;
        wdot[17*npt+i] -= qdot;

        /*reaction 70: CH2OH + O2 <=> CH2O + HO2 */
        phi_f = sc[17*npt+i]*sc[19*npt+i];
        k_f = k_f_s[69*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[15*npt+i]*sc[20*npt+i];
        Kc = Kc_s[69*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[15*npt+i] += qdot;
        wdot[17*npt+i] -= qdot;
        wdot[19*npt+i] -= qdot;
        wdot[20*npt+i] += qdot;

        /*reaction 71: CH2OH + O2 <=> CH2O + HO2 */
        phi_f = sc[17*npt+i]*sc[19*npt+i];
        k_f = k_f_s[70*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[15*npt+i]*sc[20*npt+i];
        Kc = Kc_s[70*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[15*npt+i] += qdot;
        wdot[17*npt+i] -= qdot;
        wdot[19*npt+i] -= qdot;
        wdot[20*npt+i] += qdot;

        /*reaction 72: CH2OH + HO2 <=> CH2O + H2O2 */
        phi_f = sc[17*npt+i]*sc[20*npt+i];
        k_f = k_f_s[71*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[15*npt+i]*sc[21*npt+i];
        Kc = Kc_s[71*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[15*npt+i] += qdot;
        wdot[17*npt+i] -= qdot;
        wdot[20*npt+i] -= qdot;
        wdot[21*npt+i] += qdot;

        /*reaction 73: CH2OH + HCO <=> CH2O + CH2O */
        phi_f = sc[13*npt+i]*sc[17*npt+i];
        k_f = k_f_s[72*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[15*npt+i]*sc[15*npt+i];
        Kc = Kc_s[72*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[13*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;
        wdot[15*npt+i] += qdot;
        wdot[17*npt+i] -= qdot;

        /*reaction 74: CH3O + H <=> CH3 + OH */
        phi_f = sc[0*npt+i]*sc[18*npt+i];
        k_f = k_f_s[73*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[7*npt+i];
        Kc = Kc_s[73*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[7*npt+i] += qdot;
        wdot[18*npt+i] -= qdot;

        /*reaction 75: CH3O + O <=> CH2O + OH */
        phi_f = sc[5*npt+i]*sc[18*npt+i];
        k_f = k_f_s[74*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[15*npt+i];
        Kc = Kc_s[74*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[15*npt+i] += qdot;
        wdot[18*npt+i] -= qdot;

        /*reaction 76: CH3O + OH <=> CH2O + H2O */
        phi_f = sc[7*npt+i]*sc[18*npt+i];
        k_f = k_f_s[75*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i]*sc[15*npt+i];
        Kc = Kc_s[75*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[7*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[15*npt+i] += qdot;
        wdot[18*npt+i] -= qdot;

        /*reaction 77: CH3O + O2 <=> CH2O + HO2 */
        phi_f = sc[18*npt+i]*sc[19*npt+i];
        k_f = k_f_s[76*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[15*npt+i]*sc[20*npt+i];
        Kc = Kc_s[76*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[15*npt+i] += qdot;
        wdot[18*npt+i] -= qdot;
        wdot[19*npt+i] -= qdot;
        wdot[20*npt+i] += qdot;

        /*reaction 78: CH3O + O2 <=> CH2O + HO2 */
        phi_f = sc[18*npt+i]*sc[19*npt+i];
        k_f = k_f_s[77*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[15*npt+i]*sc[20*npt+i];
        Kc = Kc_s[77*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[15*npt+i] += qdot;
        wdot[18*npt+i] -= qdot;
        wdot[19*npt+i] -= qdot;
        wdot[20*npt+i] += qdot;

        /*reaction 79: CH3O + HO2 <=> CH2O + H2O2 */
        phi_f = sc[18*npt+i]*sc[20*npt+i];
        k_f = k_f_s[78*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[15*npt+i]*sc[21*npt+i];
        Kc = Kc_s[78*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[15*npt+i] += qdot;
        wdot[18*npt+i] -= qdot;
        wdot[20*npt+i] -= qdot;
        wdot[21*npt+i] += qdot;

        /*reaction 80: CH3O + CO <=> CH3 + CO2 */
        phi_f = sc[11*npt+i]*sc[18*npt+i];
        k_f = k_f_s[79*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[22*npt+i];
        Kc = Kc_s[79*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;
        wdot[18*npt+i] -= qdot;
        wdot[22*npt+i] += qdot;

        /*reaction 81: CH3 + CH3 <=> H + C2H5 */
        phi_f = sc[4*npt+i]*sc[4*npt+i];
        k_f = k_f_s[80*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[14*npt+i];
        Kc = Kc_s[80*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[4*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;

        /*reaction 82: CH4 + CH2 <=> CH3 + CH3 */
        phi_f = sc[2*npt+i]*sc[6*npt+i];
        k_f = k_f_s[81*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[4*npt+i];
        Kc = Kc_s[81*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[4*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 83: CH4 + CH2(S) <=> CH3 + CH3 */
        phi_f = sc[3*npt+i]*sc[6*npt+i];
        k_f = k_f_s[82*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[4*npt+i];
        Kc = Kc_s[82*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[4*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 84: CH3 + OH <=> CH2 + H2O */
        phi_f = sc[4*npt+i]*sc[7*npt+i];
        k_f = k_f_s[83*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[8*npt+i];
        Kc = Kc_s[83*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[7*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;

        /*reaction 85: CH3 + OH <=> CH2(S) + H2O */
        phi_f = sc[4*npt+i]*sc[7*npt+i];
        k_f = k_f_s[84*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[8*npt+i];
        Kc = Kc_s[84*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[7*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;

        /*reaction 86: CH3 + CH2 <=> C2H4 + H */
        phi_f = sc[2*npt+i]*sc[4*npt+i];
        k_f = k_f_s[85*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[12*npt+i];
        Kc = Kc_s[85*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;

        /*reaction 87: CH3 + CH2(S) <=> C2H4 + H */
        phi_f = sc[3*npt+i]*sc[4*npt+i];
        k_f = k_f_s[86*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[12*npt+i];
        Kc = Kc_s[86*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;

        /*reaction 88: CH3O + H <=> CH2(S) + H2O */
        phi_f = sc[0*npt+i]*sc[18*npt+i];
        k_f = k_f_s[87*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[8*npt+i];
        Kc = Kc_s[87*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[8*npt+i] += qdot;
        wdot[18*npt+i] -= qdot;

        /*reaction 89: C2H6 + H <=> C2H5 + H2 */
        phi_f = sc[0*npt+i]*sc[16*npt+i];
        k_f = k_f_s[88*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[14*npt+i];
        Kc = Kc_s[88*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[14*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;

        /*reaction 90: C2H6 + O <=> C2H5 + OH */
        phi_f = sc[5*npt+i]*sc[16*npt+i];
        k_f = k_f_s[89*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[14*npt+i];
        Kc = Kc_s[89*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[14*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;

        /*reaction 91: C2H6 + OH <=> C2H5 + H2O */
        phi_f = sc[7*npt+i]*sc[16*npt+i];
        k_f = k_f_s[90*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i]*sc[14*npt+i];
        Kc = Kc_s[90*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[7*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[14*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;

        /*reaction 92: C2H6 + O2 <=> C2H5 + HO2 */
        phi_f = sc[16*npt+i]*sc[19*npt+i];
        k_f = k_f_s[91*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[14*npt+i]*sc[20*npt+i];
        Kc = Kc_s[91*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[14*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;
        wdot[19*npt+i] -= qdot;
        wdot[20*npt+i] += qdot;

        /*reaction 93: C2H6 + HO2 <=> C2H5 + H2O2 */
        phi_f = sc[16*npt+i]*sc[20*npt+i];
        k_f = k_f_s[92*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[14*npt+i]*sc[21*npt+i];
        Kc = Kc_s[92*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[14*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;
        wdot[20*npt+i] -= qdot;
        wdot[21*npt+i] += qdot;

        /*reaction 94: C2H6 + CH3 <=> C2H5 + CH4 */
        phi_f = sc[4*npt+i]*sc[16*npt+i];
        k_f = k_f_s[93*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[14*npt+i];
        Kc = Kc_s[93*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[14*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;

        /*reaction 95: C2H5 + H <=> C2H4 + H2 */
        phi_f = sc[0*npt+i]*sc[14*npt+i];
        k_f = k_f_s[94*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[12*npt+i];
        Kc = Kc_s[94*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[12*npt+i] += qdot;
        wdot[14*npt+i] -= qdot;

        /*reaction 96: C2H5 + O <=> CH3 + CH2O */
        phi_f = sc[5*npt+i]*sc[14*npt+i];
        k_f = k_f_s[95*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[15*npt+i];
        Kc = Kc_s[95*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[14*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;

        /*reaction 97: C2H5 + O2 <=> C2H4 + HO2 */
        phi_f = sc[14*npt+i]*sc[19*npt+i];
        k_f = k_f_s[96*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[12*npt+i]*sc[20*npt+i];
        Kc = Kc_s[96*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[12*npt+i] += qdot;
        wdot[14*npt+i] -= qdot;
        wdot[19*npt+i] -= qdot;
        wdot[20*npt+i] += qdot;

        /*reaction 98: C2H5 + C2H5 <=> C2H4 + C2H6 */
        phi_f = sc[14*npt+i]*sc[14*npt+i];
        k_f = k_f_s[97*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[12*npt+i]*sc[16*npt+i];
        Kc = Kc_s[97*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[12*npt+i] += qdot;
        wdot[14*npt+i] -= qdot;
        wdot[14*npt+i] -= qdot;
        wdot[16*npt+i] += qdot;

        /*reaction 99: C2H5 + HCO <=> C2H6 + CO */
        phi_f = sc[13*npt+i]*sc[14*npt+i];
        k_f = k_f_s[98*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[11*npt+i]*sc[16*npt+i];
        Kc = Kc_s[98*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[11*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;
        wdot[14*npt+i] -= qdot;
        wdot[16*npt+i] += qdot;

        /*reaction 100: C2H5 + O <=> CH3HCO + H */
        phi_f = sc[5*npt+i]*sc[14*npt+i];
        k_f = k_f_s[99*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[23*npt+i];
        Kc = Kc_s[99*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[14*npt+i] -= qdot;
        wdot[23*npt+i] += qdot;
    }
}

void vcomp_wdot_101_150(int npt, amrex::Real *  wdot, amrex::Real *  mixture, amrex::Real *  sc,
		amrex::Real *  k_f_s, amrex::Real *  Kc_s,
		amrex::Real *  tc, amrex::Real *  invT, amrex::Real *  T)
{
    for (int i=0; i<npt; i++) {
        amrex::Real qdot, q_f, q_r, phi_f, phi_r, k_f, k_r, Kc;
        amrex::Real alpha;

        /*reaction 101: C2H4 + H <=> C2H3 + H2 */
        phi_f = sc[0*npt+i]*sc[12*npt+i];
        k_f = k_f_s[100*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[10*npt+i];
        Kc = Kc_s[100*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[10*npt+i] += qdot;
        wdot[12*npt+i] -= qdot;

        /*reaction 102: C2H4 + OH <=> C2H3 + H2O */
        phi_f = sc[7*npt+i]*sc[12*npt+i];
        k_f = k_f_s[101*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i]*sc[10*npt+i];
        Kc = Kc_s[101*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[7*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[10*npt+i] += qdot;
        wdot[12*npt+i] -= qdot;

        /*reaction 103: C2H4 + CH3 <=> C2H3 + CH4 */
        phi_f = sc[4*npt+i]*sc[12*npt+i];
        k_f = k_f_s[102*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[10*npt+i];
        Kc = Kc_s[102*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[10*npt+i] += qdot;
        wdot[12*npt+i] -= qdot;

        /*reaction 104: C2H4 + O <=> CH3 + HCO */
        phi_f = sc[5*npt+i]*sc[12*npt+i];
        k_f = k_f_s[103*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[13*npt+i];
        Kc = Kc_s[103*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[12*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;

        /*reaction 105: C2H3 + OH <=> C2H2 + H2O */
        phi_f = sc[7*npt+i]*sc[10*npt+i];
        k_f = k_f_s[104*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i]*sc[9*npt+i];
        Kc = Kc_s[104*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[7*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[9*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;

        /*reaction 106: C2H4 + O <=> OH + C2H3 */
        phi_f = sc[5*npt+i]*sc[12*npt+i];
        k_f = k_f_s[105*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[10*npt+i];
        Kc = Kc_s[105*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[10*npt+i] += qdot;
        wdot[12*npt+i] -= qdot;

        /*reaction 107: C2H4 + O2 <=> C2H3 + HO2 */
        phi_f = sc[12*npt+i]*sc[19*npt+i];
        k_f = k_f_s[106*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[10*npt+i]*sc[20*npt+i];
        Kc = Kc_s[106*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[10*npt+i] += qdot;
        wdot[12*npt+i] -= qdot;
        wdot[19*npt+i] -= qdot;
        wdot[20*npt+i] += qdot;

        /*reaction 108: C2H3 + H <=> C2H2 + H2 */
        phi_f = sc[0*npt+i]*sc[10*npt+i];
        k_f = k_f_s[107*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[9*npt+i];
        Kc = Kc_s[107*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[9*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;

        /*reaction 109: C2H3 + H2O2 <=> C2H4 + HO2 */
        phi_f = sc[10*npt+i]*sc[21*npt+i];
        k_f = k_f_s[108*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[12*npt+i]*sc[20*npt+i];
        Kc = Kc_s[108*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[10*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;
        wdot[20*npt+i] += qdot;
        wdot[21*npt+i] -= qdot;

        /*reaction 110: C2H3 + CH3 <=> C2H2 + CH4 */
        phi_f = sc[4*npt+i]*sc[10*npt+i];
        k_f = k_f_s[109*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[9*npt+i];
        Kc = Kc_s[109*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[9*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;

        /*reaction 111: C2H3 + C2H3 <=> C2H4 + C2H2 */
        phi_f = sc[10*npt+i]*sc[10*npt+i];
        k_f = k_f_s[110*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[9*npt+i]*sc[12*npt+i];
        Kc = Kc_s[110*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[9*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;
        wdot[10*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;

        /*reaction 112: C2H3 + O2 <=> HCO + CH2O */
        phi_f = sc[10*npt+i]*sc[19*npt+i];
        k_f = k_f_s[111*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[13*npt+i]*sc[15*npt+i];
        Kc = Kc_s[111*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[10*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;
        wdot[15*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 113: C2H3 + O2 <=> HO2 + C2H2 */
        phi_f = sc[10*npt+i]*sc[19*npt+i];
        k_f = k_f_s[112*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[9*npt+i]*sc[20*npt+i];
        Kc = Kc_s[112*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[9*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;
        wdot[19*npt+i] -= qdot;
        wdot[20*npt+i] += qdot;

        /*reaction 114: C2H2 + O <=> CH2 + CO */
        phi_f = sc[5*npt+i]*sc[9*npt+i];
        k_f = k_f_s[113*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[11*npt+i];
        Kc = Kc_s[113*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[9*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;

        /*reaction 115: C2H2 + OH <=> CH3 + CO */
        phi_f = sc[7*npt+i]*sc[9*npt+i];
        k_f = k_f_s[114*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[11*npt+i];
        Kc = Kc_s[114*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;
        wdot[9*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;

        /*reaction 116: CH2 + O <=> HCO + H */
        phi_f = sc[2*npt+i]*sc[5*npt+i];
        k_f = k_f_s[115*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[13*npt+i];
        Kc = Kc_s[115*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[5*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;

        /*reaction 117: CH2 + OH <=> CH2O + H */
        phi_f = sc[2*npt+i]*sc[7*npt+i];
        k_f = k_f_s[116*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[15*npt+i];
        Kc = Kc_s[116*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[7*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;

        /*reaction 118: CH2 + H2 <=> H + CH3 */
        phi_f = sc[1*npt+i]*sc[2*npt+i];
        k_f = k_f_s[117*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[4*npt+i];
        Kc = Kc_s[117*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;

        /*reaction 119: CH2 + O2 <=> HCO + OH */
        phi_f = sc[2*npt+i]*sc[19*npt+i];
        k_f = k_f_s[118*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[13*npt+i];
        Kc = Kc_s[118*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[13*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 120: CH2 + HO2 <=> CH2O + OH */
        phi_f = sc[2*npt+i]*sc[20*npt+i];
        k_f = k_f_s[119*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[15*npt+i];
        Kc = Kc_s[119*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[15*npt+i] += qdot;
        wdot[20*npt+i] -= qdot;

        /*reaction 121: CH2 + CH2 <=> C2H2 + H2 */
        phi_f = sc[2*npt+i]*sc[2*npt+i];
        k_f = k_f_s[120*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[9*npt+i];
        Kc = Kc_s[120*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[2*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;

        /*reaction 122: CH2(S) + H2O <=> CH2 + H2O */
        phi_f = sc[3*npt+i]*sc[8*npt+i];
        k_f = k_f_s[121*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[8*npt+i];
        Kc = Kc_s[121*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[8*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;

        /*reaction 123: CH2(S) + CO <=> CH2 + CO */
        phi_f = sc[3*npt+i]*sc[11*npt+i];
        k_f = k_f_s[122*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[11*npt+i];
        Kc = Kc_s[122*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[11*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;

        /*reaction 124: CH2(S) + CO2 <=> CH2 + CO2 */
        phi_f = sc[3*npt+i]*sc[22*npt+i];
        k_f = k_f_s[123*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[22*npt+i];
        Kc = Kc_s[123*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[22*npt+i] -= qdot;
        wdot[22*npt+i] += qdot;

        /*reaction 125: CH2(S) + O <=> CO + H2 */
        phi_f = sc[3*npt+i]*sc[5*npt+i];
        k_f = k_f_s[124*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[11*npt+i];
        Kc = Kc_s[124*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[5*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;

        /*reaction 126: CH2(S) + O <=> HCO + H */
        phi_f = sc[3*npt+i]*sc[5*npt+i];
        k_f = k_f_s[125*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[13*npt+i];
        Kc = Kc_s[125*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[5*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;

        /*reaction 127: CH2(S) + OH <=> CH2O + H */
        phi_f = sc[3*npt+i]*sc[7*npt+i];
        k_f = k_f_s[126*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[15*npt+i];
        Kc = Kc_s[126*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[7*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;

        /*reaction 128: CH2(S) + H2 <=> CH3 + H */
        phi_f = sc[1*npt+i]*sc[3*npt+i];
        k_f = k_f_s[127*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[4*npt+i];
        Kc = Kc_s[127*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;

        /*reaction 129: CH2(S) + O2 <=> H + OH + CO */
        phi_f = sc[3*npt+i]*sc[19*npt+i];
        k_f = k_f_s[128*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[7*npt+i]*sc[11*npt+i];
        Kc = Kc_s[128*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[11*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 130: CH2(S) + O2 <=> CO + H2O */
        phi_f = sc[3*npt+i]*sc[19*npt+i];
        k_f = k_f_s[129*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i]*sc[11*npt+i];
        Kc = Kc_s[129*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[11*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 131: CH2(S) + CO2 <=> CH2O + CO */
        phi_f = sc[3*npt+i]*sc[22*npt+i];
        k_f = k_f_s[130*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[11*npt+i]*sc[15*npt+i];
        Kc = Kc_s[130*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;
        wdot[15*npt+i] += qdot;
        wdot[22*npt+i] -= qdot;

        /*reaction 132: CH3HCO <=> CH3 + HCO */
        phi_f = sc[23*npt+i];
        k_f = k_f_s[131*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[13*npt+i];
        Kc = Kc_s[131*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] += qdot;
        wdot[13*npt+i] += qdot;
        wdot[23*npt+i] -= qdot;

        /*reaction 133: CH3OCH3 <=> CH3 + CH3O */
        phi_f = sc[26*npt+i];
        k_f = k_f_s[132*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[18*npt+i];
        Kc = Kc_s[132*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] += qdot;
        wdot[18*npt+i] += qdot;
        wdot[26*npt+i] -= qdot;

        /*reaction 134: CH3OCH3 + OH <=> CH3OCH2 + H2O */
        phi_f = sc[7*npt+i]*sc[26*npt+i];
        k_f = k_f_s[133*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i]*sc[24*npt+i];
        Kc = Kc_s[133*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[7*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[24*npt+i] += qdot;
        wdot[26*npt+i] -= qdot;

        /*reaction 135: CH3OCH3 + H <=> CH3OCH2 + H2 */
        phi_f = sc[0*npt+i]*sc[26*npt+i];
        k_f = k_f_s[134*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[24*npt+i];
        Kc = Kc_s[134*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[24*npt+i] += qdot;
        wdot[26*npt+i] -= qdot;

        /*reaction 136: CH3OCH3 + CH3 <=> CH3OCH2 + CH4 */
        phi_f = sc[4*npt+i]*sc[26*npt+i];
        k_f = k_f_s[135*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[24*npt+i];
        Kc = Kc_s[135*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[24*npt+i] += qdot;
        wdot[26*npt+i] -= qdot;

        /*reaction 137: CH3OCH3 + O <=> CH3OCH2 + OH */
        phi_f = sc[5*npt+i]*sc[26*npt+i];
        k_f = k_f_s[136*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[24*npt+i];
        Kc = Kc_s[136*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[24*npt+i] += qdot;
        wdot[26*npt+i] -= qdot;

        /*reaction 138: CH3OCH3 + HO2 <=> CH3OCH2 + H2O2 */
        phi_f = sc[20*npt+i]*sc[26*npt+i];
        k_f = k_f_s[137*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[21*npt+i]*sc[24*npt+i];
        Kc = Kc_s[137*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[20*npt+i] -= qdot;
        wdot[21*npt+i] += qdot;
        wdot[24*npt+i] += qdot;
        wdot[26*npt+i] -= qdot;

        /*reaction 139: CH3OCH3 + O2 <=> CH3OCH2 + HO2 */
        phi_f = sc[19*npt+i]*sc[26*npt+i];
        k_f = k_f_s[138*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[20*npt+i]*sc[24*npt+i];
        Kc = Kc_s[138*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[19*npt+i] -= qdot;
        wdot[20*npt+i] += qdot;
        wdot[24*npt+i] += qdot;
        wdot[26*npt+i] -= qdot;

        /*reaction 140: CH3OCH2 <=> CH2O + CH3 */
        phi_f = sc[24*npt+i];
        k_f = k_f_s[139*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[15*npt+i];
        Kc = Kc_s[139*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] += qdot;
        wdot[15*npt+i] += qdot;
        wdot[24*npt+i] -= qdot;

        /*reaction 141: CH3OCH2 + CH3O <=> CH3OCH3 + CH2O */
        phi_f = sc[18*npt+i]*sc[24*npt+i];
        k_f = k_f_s[140*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[15*npt+i]*sc[26*npt+i];
        Kc = Kc_s[140*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[15*npt+i] += qdot;
        wdot[18*npt+i] -= qdot;
        wdot[24*npt+i] -= qdot;
        wdot[26*npt+i] += qdot;

        /*reaction 142: CH3OCH2 + CH2O <=> CH3OCH3 + HCO */
        phi_f = sc[15*npt+i]*sc[24*npt+i];
        k_f = k_f_s[141*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[13*npt+i]*sc[26*npt+i];
        Kc = Kc_s[141*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[13*npt+i] += qdot;
        wdot[15*npt+i] -= qdot;
        wdot[24*npt+i] -= qdot;
        wdot[26*npt+i] += qdot;

        /*reaction 143: CH3OCH2 + HO2 <=> CH3OCH2O + OH */
        phi_f = sc[20*npt+i]*sc[24*npt+i];
        k_f = k_f_s[142*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[30*npt+i];
        Kc = Kc_s[142*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[7*npt+i] += qdot;
        wdot[20*npt+i] -= qdot;
        wdot[24*npt+i] -= qdot;
        wdot[30*npt+i] += qdot;

        /*reaction 144: CH3OCH2O <=> CH3OCHO + H */
        phi_f = sc[30*npt+i];
        k_f = k_f_s[143*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[29*npt+i];
        Kc = Kc_s[143*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[29*npt+i] += qdot;
        wdot[30*npt+i] -= qdot;

        /*reaction 145: CH3OCHO + O2 <=> CH3OCO + HO2 */
        phi_f = sc[19*npt+i]*sc[29*npt+i];
        k_f = k_f_s[144*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[20*npt+i]*sc[28*npt+i];
        Kc = Kc_s[144*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[19*npt+i] -= qdot;
        wdot[20*npt+i] += qdot;
        wdot[28*npt+i] += qdot;
        wdot[29*npt+i] -= qdot;

        /*reaction 146: CH3OCHO + OH <=> CH3OCO + H2O */
        phi_f = sc[7*npt+i]*sc[29*npt+i];
        k_f = k_f_s[145*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i]*sc[28*npt+i];
        Kc = Kc_s[145*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[7*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[28*npt+i] += qdot;
        wdot[29*npt+i] -= qdot;

        /*reaction 147: CH3OCHO + HO2 <=> CH3OCO + H2O2 */
        phi_f = sc[20*npt+i]*sc[29*npt+i];
        k_f = k_f_s[146*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[21*npt+i]*sc[28*npt+i];
        Kc = Kc_s[146*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[20*npt+i] -= qdot;
        wdot[21*npt+i] += qdot;
        wdot[28*npt+i] += qdot;
        wdot[29*npt+i] -= qdot;

        /*reaction 148: CH3OCHO + O <=> CH3OCO + OH */
        phi_f = sc[5*npt+i]*sc[29*npt+i];
        k_f = k_f_s[147*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[28*npt+i];
        Kc = Kc_s[147*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[28*npt+i] += qdot;
        wdot[29*npt+i] -= qdot;

        /*reaction 149: CH3OCHO + H <=> CH3OCO + H2 */
        phi_f = sc[0*npt+i]*sc[29*npt+i];
        k_f = k_f_s[148*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[28*npt+i];
        Kc = Kc_s[148*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[28*npt+i] += qdot;
        wdot[29*npt+i] -= qdot;

        /*reaction 150: CH3OCHO + CH3 <=> CH3OCO + CH4 */
        phi_f = sc[4*npt+i]*sc[29*npt+i];
        k_f = k_f_s[149*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[28*npt+i];
        Kc = Kc_s[149*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[28*npt+i] += qdot;
        wdot[29*npt+i] -= qdot;
    }
}

void vcomp_wdot_151_175(int npt, amrex::Real *  wdot, amrex::Real *  mixture, amrex::Real *  sc,
		amrex::Real *  k_f_s, amrex::Real *  Kc_s,
		amrex::Real *  tc, amrex::Real *  invT, amrex::Real *  T)
{
    for (int i=0; i<npt; i++) {
        amrex::Real qdot, q_f, q_r, phi_f, phi_r, k_f, k_r, Kc;
        amrex::Real alpha;

        /*reaction 151: CH3OCO <=> CH3O + CO */
        phi_f = sc[28*npt+i];
        k_f = k_f_s[150*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[11*npt+i]*sc[18*npt+i];
        Kc = Kc_s[150*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[11*npt+i] += qdot;
        wdot[18*npt+i] += qdot;
        wdot[28*npt+i] -= qdot;

        /*reaction 152: CH3OCO <=> CH3 + CO2 */
        phi_f = sc[28*npt+i];
        k_f = k_f_s[151*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[22*npt+i];
        Kc = Kc_s[151*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] += qdot;
        wdot[22*npt+i] += qdot;
        wdot[28*npt+i] -= qdot;

        /*reaction 153: CH3OCH2 + O2 <=> CH3OCH2O2 */
        phi_f = sc[19*npt+i]*sc[24*npt+i];
        k_f = k_f_s[152*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[34*npt+i];
        Kc = Kc_s[152*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[19*npt+i] -= qdot;
        wdot[24*npt+i] -= qdot;
        wdot[34*npt+i] += qdot;

        /*reaction 154: CH3OCH2O2 + CH3OCH2O2 <=> O2 + CH3OCH2O + CH3OCH2O */
        phi_f = sc[34*npt+i]*sc[34*npt+i];
        k_f = k_f_s[153*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[19*npt+i]*sc[30*npt+i]*sc[30*npt+i];
        Kc = Kc_s[153*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[19*npt+i] += qdot;
        wdot[30*npt+i] += qdot;
        wdot[30*npt+i] += qdot;
        wdot[34*npt+i] -= qdot;
        wdot[34*npt+i] -= qdot;

        /*reaction 155: CH3OCH2O2 + CH3OCH2O2 <=> O2 + CH3OCHO + CH3OCH2OH */
        phi_f = sc[34*npt+i]*sc[34*npt+i];
        k_f = k_f_s[154*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[19*npt+i]*sc[29*npt+i]*sc[31*npt+i];
        Kc = Kc_s[154*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[19*npt+i] += qdot;
        wdot[29*npt+i] += qdot;
        wdot[31*npt+i] += qdot;
        wdot[34*npt+i] -= qdot;
        wdot[34*npt+i] -= qdot;

        /*reaction 156: CH3OCH2O <=> CH3O + CH2O */
        phi_f = sc[30*npt+i];
        k_f = k_f_s[155*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[15*npt+i]*sc[18*npt+i];
        Kc = Kc_s[155*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[15*npt+i] += qdot;
        wdot[18*npt+i] += qdot;
        wdot[30*npt+i] -= qdot;

        /*reaction 157: CH3OCH2O + O2 <=> CH3OCHO + HO2 */
        phi_f = sc[19*npt+i]*sc[30*npt+i];
        k_f = k_f_s[156*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[20*npt+i]*sc[29*npt+i];
        Kc = Kc_s[156*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[19*npt+i] -= qdot;
        wdot[20*npt+i] += qdot;
        wdot[29*npt+i] += qdot;
        wdot[30*npt+i] -= qdot;

        /*reaction 158: CH3OCH2O2 <=> CH2OCH2O2H */
        phi_f = sc[34*npt+i];
        k_f = k_f_s[157*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[35*npt+i];
        Kc = Kc_s[157*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[34*npt+i] -= qdot;
        wdot[35*npt+i] += qdot;

        /*reaction 159: CH2OCH2O2H <=> OH + CH2O + CH2O */
        phi_f = sc[35*npt+i];
        k_f = k_f_s[158*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[15*npt+i]*sc[15*npt+i];
        Kc = Kc_s[158*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[7*npt+i] += qdot;
        wdot[15*npt+i] += qdot;
        wdot[15*npt+i] += qdot;
        wdot[35*npt+i] -= qdot;

        /*reaction 160: CH2OCH2O2H + O2 <=> O2CH2OCH2O2H */
        phi_f = sc[19*npt+i]*sc[35*npt+i];
        k_f = k_f_s[159*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[37*npt+i];
        Kc = Kc_s[159*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[19*npt+i] -= qdot;
        wdot[35*npt+i] -= qdot;
        wdot[37*npt+i] += qdot;

        /*reaction 161: O2CH2OCH2O2H <=> HO2CH2OCHO + OH */
        phi_f = sc[37*npt+i];
        k_f = k_f_s[160*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[36*npt+i];
        Kc = Kc_s[160*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[7*npt+i] += qdot;
        wdot[36*npt+i] += qdot;
        wdot[37*npt+i] -= qdot;

        /*reaction 162: HO2CH2OCHO <=> OCH2OCHO + OH */
        phi_f = sc[36*npt+i];
        k_f = k_f_s[161*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[32*npt+i];
        Kc = Kc_s[161*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[7*npt+i] += qdot;
        wdot[32*npt+i] += qdot;
        wdot[36*npt+i] -= qdot;

        /*reaction 163: OCH2OCHO <=> HOCH2OCO */
        phi_f = sc[32*npt+i];
        k_f = k_f_s[162*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[33*npt+i];
        Kc = Kc_s[162*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[32*npt+i] -= qdot;
        wdot[33*npt+i] += qdot;

        /*reaction 164: HOCH2OCO <=> HOCH2O + CO */
        phi_f = sc[33*npt+i];
        k_f = k_f_s[163*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[11*npt+i]*sc[27*npt+i];
        Kc = Kc_s[163*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[11*npt+i] += qdot;
        wdot[27*npt+i] += qdot;
        wdot[33*npt+i] -= qdot;

        /*reaction 165: HOCH2OCO <=> CH2OH + CO2 */
        phi_f = sc[33*npt+i];
        k_f = k_f_s[164*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[17*npt+i]*sc[22*npt+i];
        Kc = Kc_s[164*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[17*npt+i] += qdot;
        wdot[22*npt+i] += qdot;
        wdot[33*npt+i] -= qdot;

        /*reaction 166: HOCH2O <=> HCOOH + H */
        phi_f = sc[27*npt+i];
        k_f = k_f_s[165*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[25*npt+i];
        Kc = Kc_s[165*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[25*npt+i] += qdot;
        wdot[27*npt+i] -= qdot;

        /*reaction 167: CH2O + OH <=> HOCH2O */
        phi_f = sc[7*npt+i]*sc[15*npt+i];
        k_f = k_f_s[166*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[27*npt+i];
        Kc = Kc_s[166*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[7*npt+i] -= qdot;
        wdot[15*npt+i] -= qdot;
        wdot[27*npt+i] += qdot;

        /*reaction 168: HCOOH <=> HCO + OH */
        phi_f = sc[25*npt+i];
        k_f = k_f_s[167*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[13*npt+i];
        Kc = Kc_s[167*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[7*npt+i] += qdot;
        wdot[13*npt+i] += qdot;
        wdot[25*npt+i] -= qdot;

        /*reaction 169: HCOOH + OH <=> H2O + CO2 + H */
        phi_f = sc[7*npt+i]*sc[25*npt+i];
        k_f = k_f_s[168*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[8*npt+i]*sc[22*npt+i];
        Kc = Kc_s[168*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[22*npt+i] += qdot;
        wdot[25*npt+i] -= qdot;

        /*reaction 170: HCOOH + OH <=> H2O + CO + OH */
        phi_f = sc[7*npt+i]*sc[25*npt+i];
        k_f = k_f_s[169*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[8*npt+i]*sc[11*npt+i];
        Kc = Kc_s[169*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[7*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[8*npt+i] += qdot;
        wdot[11*npt+i] += qdot;
        wdot[25*npt+i] -= qdot;

        /*reaction 171: HCOOH + H <=> H2 + CO2 + H */
        phi_f = sc[0*npt+i]*sc[25*npt+i];
        k_f = k_f_s[170*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[1*npt+i]*sc[22*npt+i];
        Kc = Kc_s[170*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] += qdot;
        wdot[22*npt+i] += qdot;
        wdot[25*npt+i] -= qdot;

        /*reaction 172: HCOOH + H <=> H2 + CO + OH */
        phi_f = sc[0*npt+i]*sc[25*npt+i];
        k_f = k_f_s[171*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[7*npt+i]*sc[11*npt+i];
        Kc = Kc_s[171*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[7*npt+i] += qdot;
        wdot[11*npt+i] += qdot;
        wdot[25*npt+i] -= qdot;

        /*reaction 173: HCOOH + CH3 <=> CH4 + CO + OH */
        phi_f = sc[4*npt+i]*sc[25*npt+i];
        k_f = k_f_s[172*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[7*npt+i]*sc[11*npt+i];
        Kc = Kc_s[172*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] += qdot;
        wdot[11*npt+i] += qdot;
        wdot[25*npt+i] -= qdot;

        /*reaction 174: HCOOH + HO2 <=> H2O2 + CO + OH */
        phi_f = sc[20*npt+i]*sc[25*npt+i];
        k_f = k_f_s[173*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[11*npt+i]*sc[21*npt+i];
        Kc = Kc_s[173*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[7*npt+i] += qdot;
        wdot[11*npt+i] += qdot;
        wdot[20*npt+i] -= qdot;
        wdot[21*npt+i] += qdot;
        wdot[25*npt+i] -= qdot;

        /*reaction 175: HCOOH + O <=> CO + OH + OH */
        phi_f = sc[5*npt+i]*sc[25*npt+i];
        k_f = k_f_s[174*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[7*npt+i]*sc[11*npt+i];
        Kc = Kc_s[174*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[7*npt+i] += qdot;
        wdot[11*npt+i] += qdot;
        wdot[25*npt+i] -= qdot;
    }
}


static amrex::Real T_save = -1;
#ifdef _OPENMP
#pragma omp threadprivate(T_save)
#endif

static amrex::Real k_f_save[175];
#ifdef _OPENMP
#pragma omp threadprivate(k_f_save)
#endif

static amrex::Real Kc_save[175];
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

    amrex::Real qdot, q_f[175], q_r[175];
    amrex::Real sc_qss[1];
    comp_qfqr(q_f, q_r, sc, sc_qss, tc, invT);

    for (int i = 0; i < 39; ++i) {
        wdot[i] = 0.0;
    }

    qdot = q_f[0]-q_r[0];
    wdot[0] -= qdot;
    wdot[19] -= qdot;
    wdot[20] += qdot;

    qdot = q_f[1]-q_r[1];
    wdot[7] += qdot;
    wdot[7] += qdot;
    wdot[21] -= qdot;

    qdot = q_f[2]-q_r[2];
    wdot[4] -= qdot;
    wdot[4] -= qdot;
    wdot[16] += qdot;

    qdot = q_f[3]-q_r[3];
    wdot[0] -= qdot;
    wdot[4] -= qdot;
    wdot[6] += qdot;

    qdot = q_f[4]-q_r[4];
    wdot[0] -= qdot;
    wdot[14] -= qdot;
    wdot[16] += qdot;

    qdot = q_f[5]-q_r[5];
    wdot[1] += qdot;
    wdot[9] += qdot;
    wdot[12] -= qdot;

    qdot = q_f[6]-q_r[6];
    wdot[0] -= qdot;
    wdot[12] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[7]-q_r[7];
    wdot[0] -= qdot;
    wdot[10] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[8]-q_r[8];
    wdot[0] -= qdot;
    wdot[9] -= qdot;
    wdot[10] += qdot;

    qdot = q_f[9]-q_r[9];
    wdot[0] -= qdot;
    wdot[2] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[10]-q_r[10];
    wdot[5] -= qdot;
    wdot[11] -= qdot;
    wdot[22] += qdot;

    qdot = q_f[11]-q_r[11];
    wdot[0] += qdot;
    wdot[0] += qdot;
    wdot[1] -= qdot;

    qdot = q_f[12]-q_r[12];
    wdot[5] -= qdot;
    wdot[5] -= qdot;
    wdot[19] += qdot;

    qdot = q_f[13]-q_r[13];
    wdot[0] -= qdot;
    wdot[5] -= qdot;
    wdot[7] += qdot;

    qdot = q_f[14]-q_r[14];
    wdot[0] -= qdot;
    wdot[7] -= qdot;
    wdot[8] += qdot;

    qdot = q_f[15]-q_r[15];
    wdot[0] += qdot;
    wdot[11] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[16]-q_r[16];
    wdot[0] += qdot;
    wdot[13] += qdot;
    wdot[15] -= qdot;

    qdot = q_f[17]-q_r[17];
    wdot[1] += qdot;
    wdot[11] += qdot;
    wdot[15] -= qdot;

    qdot = q_f[18]-q_r[18];
    wdot[0] += qdot;
    wdot[15] += qdot;
    wdot[17] -= qdot;

    qdot = q_f[19]-q_r[19];
    wdot[0] += qdot;
    wdot[15] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[20]-q_r[20];
    wdot[2] += qdot;
    wdot[3] -= qdot;

    qdot = q_f[21]-q_r[21];
    wdot[8] += qdot;
    wdot[11] += qdot;
    wdot[25] -= qdot;

    qdot = q_f[22]-q_r[22];
    wdot[1] += qdot;
    wdot[22] += qdot;
    wdot[25] -= qdot;

    qdot = q_f[23]-q_r[23];
    wdot[0] -= qdot;
    wdot[5] += qdot;
    wdot[7] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[24]-q_r[24];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[5] -= qdot;
    wdot[7] += qdot;

    qdot = q_f[25]-q_r[25];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[7] -= qdot;
    wdot[8] += qdot;

    qdot = q_f[26]-q_r[26];
    wdot[5] -= qdot;
    wdot[7] += qdot;
    wdot[7] += qdot;
    wdot[8] -= qdot;

    qdot = q_f[27]-q_r[27];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[19] += qdot;
    wdot[20] -= qdot;

    qdot = q_f[28]-q_r[28];
    wdot[0] -= qdot;
    wdot[7] += qdot;
    wdot[7] += qdot;
    wdot[20] -= qdot;

    qdot = q_f[29]-q_r[29];
    wdot[5] -= qdot;
    wdot[7] += qdot;
    wdot[19] += qdot;
    wdot[20] -= qdot;

    qdot = q_f[30]-q_r[30];
    wdot[7] -= qdot;
    wdot[8] += qdot;
    wdot[19] += qdot;
    wdot[20] -= qdot;

    qdot = q_f[31]-q_r[31];
    wdot[19] += qdot;
    wdot[20] -= qdot;
    wdot[20] -= qdot;
    wdot[21] += qdot;

    qdot = q_f[32]-q_r[32];
    wdot[19] += qdot;
    wdot[20] -= qdot;
    wdot[20] -= qdot;
    wdot[21] += qdot;

    qdot = q_f[33]-q_r[33];
    wdot[0] -= qdot;
    wdot[7] += qdot;
    wdot[8] += qdot;
    wdot[21] -= qdot;

    qdot = q_f[34]-q_r[34];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[20] += qdot;
    wdot[21] -= qdot;

    qdot = q_f[35]-q_r[35];
    wdot[5] -= qdot;
    wdot[7] += qdot;
    wdot[20] += qdot;
    wdot[21] -= qdot;

    qdot = q_f[36]-q_r[36];
    wdot[7] -= qdot;
    wdot[8] += qdot;
    wdot[20] += qdot;
    wdot[21] -= qdot;

    qdot = q_f[37]-q_r[37];
    wdot[7] -= qdot;
    wdot[8] += qdot;
    wdot[20] += qdot;
    wdot[21] -= qdot;

    qdot = q_f[38]-q_r[38];
    wdot[5] += qdot;
    wdot[11] -= qdot;
    wdot[19] -= qdot;
    wdot[22] += qdot;

    qdot = q_f[39]-q_r[39];
    wdot[7] += qdot;
    wdot[11] -= qdot;
    wdot[20] -= qdot;
    wdot[22] += qdot;

    qdot = q_f[40]-q_r[40];
    wdot[0] += qdot;
    wdot[7] -= qdot;
    wdot[11] -= qdot;
    wdot[22] += qdot;

    qdot = q_f[41]-q_r[41];
    wdot[11] += qdot;
    wdot[13] -= qdot;
    wdot[19] -= qdot;
    wdot[20] += qdot;

    qdot = q_f[42]-q_r[42];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[11] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[43]-q_r[43];
    wdot[5] -= qdot;
    wdot[7] += qdot;
    wdot[11] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[44]-q_r[44];
    wdot[7] -= qdot;
    wdot[8] += qdot;
    wdot[11] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[45]-q_r[45];
    wdot[0] += qdot;
    wdot[5] -= qdot;
    wdot[13] -= qdot;
    wdot[22] += qdot;

    qdot = q_f[46]-q_r[46];
    wdot[0] += qdot;
    wdot[7] += qdot;
    wdot[13] -= qdot;
    wdot[20] -= qdot;
    wdot[22] += qdot;

    qdot = q_f[47]-q_r[47];
    wdot[1] += qdot;
    wdot[11] += qdot;
    wdot[11] += qdot;
    wdot[13] -= qdot;
    wdot[13] -= qdot;

    qdot = q_f[48]-q_r[48];
    wdot[4] -= qdot;
    wdot[6] += qdot;
    wdot[11] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[49]-q_r[49];
    wdot[11] += qdot;
    wdot[13] -= qdot;
    wdot[13] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[50]-q_r[50];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[13] += qdot;
    wdot[15] -= qdot;

    qdot = q_f[51]-q_r[51];
    wdot[5] -= qdot;
    wdot[7] += qdot;
    wdot[13] += qdot;
    wdot[15] -= qdot;

    qdot = q_f[52]-q_r[52];
    wdot[7] -= qdot;
    wdot[8] += qdot;
    wdot[13] += qdot;
    wdot[15] -= qdot;

    qdot = q_f[53]-q_r[53];
    wdot[13] += qdot;
    wdot[15] -= qdot;
    wdot[19] -= qdot;
    wdot[20] += qdot;

    qdot = q_f[54]-q_r[54];
    wdot[13] += qdot;
    wdot[15] -= qdot;
    wdot[20] -= qdot;
    wdot[21] += qdot;

    qdot = q_f[55]-q_r[55];
    wdot[4] -= qdot;
    wdot[6] += qdot;
    wdot[13] += qdot;
    wdot[15] -= qdot;

    qdot = q_f[56]-q_r[56];
    wdot[0] += qdot;
    wdot[4] -= qdot;
    wdot[5] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[57]-q_r[57];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[18] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[58]-q_r[58];
    wdot[4] -= qdot;
    wdot[7] += qdot;
    wdot[15] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[59]-q_r[59];
    wdot[4] -= qdot;
    wdot[7] += qdot;
    wdot[18] += qdot;
    wdot[20] -= qdot;

    qdot = q_f[60]-q_r[60];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[4] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[61]-q_r[61];
    wdot[4] += qdot;
    wdot[5] -= qdot;
    wdot[6] -= qdot;
    wdot[7] += qdot;

    qdot = q_f[62]-q_r[62];
    wdot[4] += qdot;
    wdot[6] -= qdot;
    wdot[7] -= qdot;
    wdot[8] += qdot;

    qdot = q_f[63]-q_r[63];
    wdot[4] -= qdot;
    wdot[6] += qdot;
    wdot[19] += qdot;
    wdot[20] -= qdot;

    qdot = q_f[64]-q_r[64];
    wdot[4] += qdot;
    wdot[6] -= qdot;
    wdot[20] -= qdot;
    wdot[21] += qdot;

    qdot = q_f[65]-q_r[65];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[15] += qdot;
    wdot[17] -= qdot;

    qdot = q_f[66]-q_r[66];
    wdot[0] -= qdot;
    wdot[4] += qdot;
    wdot[7] += qdot;
    wdot[17] -= qdot;

    qdot = q_f[67]-q_r[67];
    wdot[5] -= qdot;
    wdot[7] += qdot;
    wdot[15] += qdot;
    wdot[17] -= qdot;

    qdot = q_f[68]-q_r[68];
    wdot[7] -= qdot;
    wdot[8] += qdot;
    wdot[15] += qdot;
    wdot[17] -= qdot;

    qdot = q_f[69]-q_r[69];
    wdot[15] += qdot;
    wdot[17] -= qdot;
    wdot[19] -= qdot;
    wdot[20] += qdot;

    qdot = q_f[70]-q_r[70];
    wdot[15] += qdot;
    wdot[17] -= qdot;
    wdot[19] -= qdot;
    wdot[20] += qdot;

    qdot = q_f[71]-q_r[71];
    wdot[15] += qdot;
    wdot[17] -= qdot;
    wdot[20] -= qdot;
    wdot[21] += qdot;

    qdot = q_f[72]-q_r[72];
    wdot[13] -= qdot;
    wdot[15] += qdot;
    wdot[15] += qdot;
    wdot[17] -= qdot;

    qdot = q_f[73]-q_r[73];
    wdot[0] -= qdot;
    wdot[4] += qdot;
    wdot[7] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[74]-q_r[74];
    wdot[5] -= qdot;
    wdot[7] += qdot;
    wdot[15] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[75]-q_r[75];
    wdot[7] -= qdot;
    wdot[8] += qdot;
    wdot[15] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[76]-q_r[76];
    wdot[15] += qdot;
    wdot[18] -= qdot;
    wdot[19] -= qdot;
    wdot[20] += qdot;

    qdot = q_f[77]-q_r[77];
    wdot[15] += qdot;
    wdot[18] -= qdot;
    wdot[19] -= qdot;
    wdot[20] += qdot;

    qdot = q_f[78]-q_r[78];
    wdot[15] += qdot;
    wdot[18] -= qdot;
    wdot[20] -= qdot;
    wdot[21] += qdot;

    qdot = q_f[79]-q_r[79];
    wdot[4] += qdot;
    wdot[11] -= qdot;
    wdot[18] -= qdot;
    wdot[22] += qdot;

    qdot = q_f[80]-q_r[80];
    wdot[0] += qdot;
    wdot[4] -= qdot;
    wdot[4] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[81]-q_r[81];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[4] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[82]-q_r[82];
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[4] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[83]-q_r[83];
    wdot[2] += qdot;
    wdot[4] -= qdot;
    wdot[7] -= qdot;
    wdot[8] += qdot;

    qdot = q_f[84]-q_r[84];
    wdot[3] += qdot;
    wdot[4] -= qdot;
    wdot[7] -= qdot;
    wdot[8] += qdot;

    qdot = q_f[85]-q_r[85];
    wdot[0] += qdot;
    wdot[2] -= qdot;
    wdot[4] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[86]-q_r[86];
    wdot[0] += qdot;
    wdot[3] -= qdot;
    wdot[4] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[87]-q_r[87];
    wdot[0] -= qdot;
    wdot[3] += qdot;
    wdot[8] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[88]-q_r[88];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[14] += qdot;
    wdot[16] -= qdot;

    qdot = q_f[89]-q_r[89];
    wdot[5] -= qdot;
    wdot[7] += qdot;
    wdot[14] += qdot;
    wdot[16] -= qdot;

    qdot = q_f[90]-q_r[90];
    wdot[7] -= qdot;
    wdot[8] += qdot;
    wdot[14] += qdot;
    wdot[16] -= qdot;

    qdot = q_f[91]-q_r[91];
    wdot[14] += qdot;
    wdot[16] -= qdot;
    wdot[19] -= qdot;
    wdot[20] += qdot;

    qdot = q_f[92]-q_r[92];
    wdot[14] += qdot;
    wdot[16] -= qdot;
    wdot[20] -= qdot;
    wdot[21] += qdot;

    qdot = q_f[93]-q_r[93];
    wdot[4] -= qdot;
    wdot[6] += qdot;
    wdot[14] += qdot;
    wdot[16] -= qdot;

    qdot = q_f[94]-q_r[94];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[12] += qdot;
    wdot[14] -= qdot;

    qdot = q_f[95]-q_r[95];
    wdot[4] += qdot;
    wdot[5] -= qdot;
    wdot[14] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[96]-q_r[96];
    wdot[12] += qdot;
    wdot[14] -= qdot;
    wdot[19] -= qdot;
    wdot[20] += qdot;

    qdot = q_f[97]-q_r[97];
    wdot[12] += qdot;
    wdot[14] -= qdot;
    wdot[14] -= qdot;
    wdot[16] += qdot;

    qdot = q_f[98]-q_r[98];
    wdot[11] += qdot;
    wdot[13] -= qdot;
    wdot[14] -= qdot;
    wdot[16] += qdot;

    qdot = q_f[99]-q_r[99];
    wdot[0] += qdot;
    wdot[5] -= qdot;
    wdot[14] -= qdot;
    wdot[23] += qdot;

    qdot = q_f[100]-q_r[100];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[10] += qdot;
    wdot[12] -= qdot;

    qdot = q_f[101]-q_r[101];
    wdot[7] -= qdot;
    wdot[8] += qdot;
    wdot[10] += qdot;
    wdot[12] -= qdot;

    qdot = q_f[102]-q_r[102];
    wdot[4] -= qdot;
    wdot[6] += qdot;
    wdot[10] += qdot;
    wdot[12] -= qdot;

    qdot = q_f[103]-q_r[103];
    wdot[4] += qdot;
    wdot[5] -= qdot;
    wdot[12] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[104]-q_r[104];
    wdot[7] -= qdot;
    wdot[8] += qdot;
    wdot[9] += qdot;
    wdot[10] -= qdot;

    qdot = q_f[105]-q_r[105];
    wdot[5] -= qdot;
    wdot[7] += qdot;
    wdot[10] += qdot;
    wdot[12] -= qdot;

    qdot = q_f[106]-q_r[106];
    wdot[10] += qdot;
    wdot[12] -= qdot;
    wdot[19] -= qdot;
    wdot[20] += qdot;

    qdot = q_f[107]-q_r[107];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[9] += qdot;
    wdot[10] -= qdot;

    qdot = q_f[108]-q_r[108];
    wdot[10] -= qdot;
    wdot[12] += qdot;
    wdot[20] += qdot;
    wdot[21] -= qdot;

    qdot = q_f[109]-q_r[109];
    wdot[4] -= qdot;
    wdot[6] += qdot;
    wdot[9] += qdot;
    wdot[10] -= qdot;

    qdot = q_f[110]-q_r[110];
    wdot[9] += qdot;
    wdot[10] -= qdot;
    wdot[10] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[111]-q_r[111];
    wdot[10] -= qdot;
    wdot[13] += qdot;
    wdot[15] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[112]-q_r[112];
    wdot[9] += qdot;
    wdot[10] -= qdot;
    wdot[19] -= qdot;
    wdot[20] += qdot;

    qdot = q_f[113]-q_r[113];
    wdot[2] += qdot;
    wdot[5] -= qdot;
    wdot[9] -= qdot;
    wdot[11] += qdot;

    qdot = q_f[114]-q_r[114];
    wdot[4] += qdot;
    wdot[7] -= qdot;
    wdot[9] -= qdot;
    wdot[11] += qdot;

    qdot = q_f[115]-q_r[115];
    wdot[0] += qdot;
    wdot[2] -= qdot;
    wdot[5] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[116]-q_r[116];
    wdot[0] += qdot;
    wdot[2] -= qdot;
    wdot[7] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[117]-q_r[117];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[2] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[118]-q_r[118];
    wdot[2] -= qdot;
    wdot[7] += qdot;
    wdot[13] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[119]-q_r[119];
    wdot[2] -= qdot;
    wdot[7] += qdot;
    wdot[15] += qdot;
    wdot[20] -= qdot;

    qdot = q_f[120]-q_r[120];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[2] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[121]-q_r[121];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[8] -= qdot;
    wdot[8] += qdot;

    qdot = q_f[122]-q_r[122];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[11] -= qdot;
    wdot[11] += qdot;

    qdot = q_f[123]-q_r[123];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[22] -= qdot;
    wdot[22] += qdot;

    qdot = q_f[124]-q_r[124];
    wdot[1] += qdot;
    wdot[3] -= qdot;
    wdot[5] -= qdot;
    wdot[11] += qdot;

    qdot = q_f[125]-q_r[125];
    wdot[0] += qdot;
    wdot[3] -= qdot;
    wdot[5] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[126]-q_r[126];
    wdot[0] += qdot;
    wdot[3] -= qdot;
    wdot[7] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[127]-q_r[127];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[128]-q_r[128];
    wdot[0] += qdot;
    wdot[3] -= qdot;
    wdot[7] += qdot;
    wdot[11] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[129]-q_r[129];
    wdot[3] -= qdot;
    wdot[8] += qdot;
    wdot[11] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[130]-q_r[130];
    wdot[3] -= qdot;
    wdot[11] += qdot;
    wdot[15] += qdot;
    wdot[22] -= qdot;

    qdot = q_f[131]-q_r[131];
    wdot[4] += qdot;
    wdot[13] += qdot;
    wdot[23] -= qdot;

    qdot = q_f[132]-q_r[132];
    wdot[4] += qdot;
    wdot[18] += qdot;
    wdot[26] -= qdot;

    qdot = q_f[133]-q_r[133];
    wdot[7] -= qdot;
    wdot[8] += qdot;
    wdot[24] += qdot;
    wdot[26] -= qdot;

    qdot = q_f[134]-q_r[134];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[24] += qdot;
    wdot[26] -= qdot;

    qdot = q_f[135]-q_r[135];
    wdot[4] -= qdot;
    wdot[6] += qdot;
    wdot[24] += qdot;
    wdot[26] -= qdot;

    qdot = q_f[136]-q_r[136];
    wdot[5] -= qdot;
    wdot[7] += qdot;
    wdot[24] += qdot;
    wdot[26] -= qdot;

    qdot = q_f[137]-q_r[137];
    wdot[20] -= qdot;
    wdot[21] += qdot;
    wdot[24] += qdot;
    wdot[26] -= qdot;

    qdot = q_f[138]-q_r[138];
    wdot[19] -= qdot;
    wdot[20] += qdot;
    wdot[24] += qdot;
    wdot[26] -= qdot;

    qdot = q_f[139]-q_r[139];
    wdot[4] += qdot;
    wdot[15] += qdot;
    wdot[24] -= qdot;

    qdot = q_f[140]-q_r[140];
    wdot[15] += qdot;
    wdot[18] -= qdot;
    wdot[24] -= qdot;
    wdot[26] += qdot;

    qdot = q_f[141]-q_r[141];
    wdot[13] += qdot;
    wdot[15] -= qdot;
    wdot[24] -= qdot;
    wdot[26] += qdot;

    qdot = q_f[142]-q_r[142];
    wdot[7] += qdot;
    wdot[20] -= qdot;
    wdot[24] -= qdot;
    wdot[30] += qdot;

    qdot = q_f[143]-q_r[143];
    wdot[0] += qdot;
    wdot[29] += qdot;
    wdot[30] -= qdot;

    qdot = q_f[144]-q_r[144];
    wdot[19] -= qdot;
    wdot[20] += qdot;
    wdot[28] += qdot;
    wdot[29] -= qdot;

    qdot = q_f[145]-q_r[145];
    wdot[7] -= qdot;
    wdot[8] += qdot;
    wdot[28] += qdot;
    wdot[29] -= qdot;

    qdot = q_f[146]-q_r[146];
    wdot[20] -= qdot;
    wdot[21] += qdot;
    wdot[28] += qdot;
    wdot[29] -= qdot;

    qdot = q_f[147]-q_r[147];
    wdot[5] -= qdot;
    wdot[7] += qdot;
    wdot[28] += qdot;
    wdot[29] -= qdot;

    qdot = q_f[148]-q_r[148];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[28] += qdot;
    wdot[29] -= qdot;

    qdot = q_f[149]-q_r[149];
    wdot[4] -= qdot;
    wdot[6] += qdot;
    wdot[28] += qdot;
    wdot[29] -= qdot;

    qdot = q_f[150]-q_r[150];
    wdot[11] += qdot;
    wdot[18] += qdot;
    wdot[28] -= qdot;

    qdot = q_f[151]-q_r[151];
    wdot[4] += qdot;
    wdot[22] += qdot;
    wdot[28] -= qdot;

    qdot = q_f[152]-q_r[152];
    wdot[19] -= qdot;
    wdot[24] -= qdot;
    wdot[34] += qdot;

    qdot = q_f[153]-q_r[153];
    wdot[19] += qdot;
    wdot[30] += qdot;
    wdot[30] += qdot;
    wdot[34] -= qdot;
    wdot[34] -= qdot;

    qdot = q_f[154]-q_r[154];
    wdot[19] += qdot;
    wdot[29] += qdot;
    wdot[31] += qdot;
    wdot[34] -= qdot;
    wdot[34] -= qdot;

    qdot = q_f[155]-q_r[155];
    wdot[15] += qdot;
    wdot[18] += qdot;
    wdot[30] -= qdot;

    qdot = q_f[156]-q_r[156];
    wdot[19] -= qdot;
    wdot[20] += qdot;
    wdot[29] += qdot;
    wdot[30] -= qdot;

    qdot = q_f[157]-q_r[157];
    wdot[34] -= qdot;
    wdot[35] += qdot;

    qdot = q_f[158]-q_r[158];
    wdot[7] += qdot;
    wdot[15] += qdot;
    wdot[15] += qdot;
    wdot[35] -= qdot;

    qdot = q_f[159]-q_r[159];
    wdot[19] -= qdot;
    wdot[35] -= qdot;
    wdot[37] += qdot;

    qdot = q_f[160]-q_r[160];
    wdot[7] += qdot;
    wdot[36] += qdot;
    wdot[37] -= qdot;

    qdot = q_f[161]-q_r[161];
    wdot[7] += qdot;
    wdot[32] += qdot;
    wdot[36] -= qdot;

    qdot = q_f[162]-q_r[162];
    wdot[32] -= qdot;
    wdot[33] += qdot;

    qdot = q_f[163]-q_r[163];
    wdot[11] += qdot;
    wdot[27] += qdot;
    wdot[33] -= qdot;

    qdot = q_f[164]-q_r[164];
    wdot[17] += qdot;
    wdot[22] += qdot;
    wdot[33] -= qdot;

    qdot = q_f[165]-q_r[165];
    wdot[0] += qdot;
    wdot[25] += qdot;
    wdot[27] -= qdot;

    qdot = q_f[166]-q_r[166];
    wdot[7] -= qdot;
    wdot[15] -= qdot;
    wdot[27] += qdot;

    qdot = q_f[167]-q_r[167];
    wdot[7] += qdot;
    wdot[13] += qdot;
    wdot[25] -= qdot;

    qdot = q_f[168]-q_r[168];
    wdot[0] += qdot;
    wdot[7] -= qdot;
    wdot[8] += qdot;
    wdot[22] += qdot;
    wdot[25] -= qdot;

    qdot = q_f[169]-q_r[169];
    wdot[7] -= qdot;
    wdot[7] += qdot;
    wdot[8] += qdot;
    wdot[11] += qdot;
    wdot[25] -= qdot;

    qdot = q_f[170]-q_r[170];
    wdot[0] -= qdot;
    wdot[0] += qdot;
    wdot[1] += qdot;
    wdot[22] += qdot;
    wdot[25] -= qdot;

    qdot = q_f[171]-q_r[171];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[7] += qdot;
    wdot[11] += qdot;
    wdot[25] -= qdot;

    qdot = q_f[172]-q_r[172];
    wdot[4] -= qdot;
    wdot[6] += qdot;
    wdot[7] += qdot;
    wdot[11] += qdot;
    wdot[25] -= qdot;

    qdot = q_f[173]-q_r[173];
    wdot[7] += qdot;
    wdot[11] += qdot;
    wdot[20] -= qdot;
    wdot[21] += qdot;
    wdot[25] -= qdot;

    qdot = q_f[174]-q_r[174];
    wdot[5] -= qdot;
    wdot[7] += qdot;
    wdot[7] += qdot;
    wdot[11] += qdot;
    wdot[25] -= qdot;

    return;
}

void comp_k_f(amrex::Real *  tc, amrex::Real invT, amrex::Real *  k_f)
{
    for (int i=0; i<175; ++i) {
        k_f[i] = prefactor_units[i] * fwd_A[i]
                    * exp(fwd_beta[i] * tc[0] - activation_units[i] * fwd_Ea[i] * invT);
    };
    return;
}

void comp_Kc(amrex::Real *  tc, amrex::Real invT, amrex::Real *  Kc)
{
    /*compute the Gibbs free energy */
    amrex::Real g_RT[39];
    gibbs(g_RT, tc);

    Kc[0] = g_RT[0] + g_RT[19] - g_RT[20];
    Kc[1] = -g_RT[7] - g_RT[7] + g_RT[21];
    Kc[2] = g_RT[4] + g_RT[4] - g_RT[16];
    Kc[3] = g_RT[0] + g_RT[4] - g_RT[6];
    Kc[4] = g_RT[0] + g_RT[14] - g_RT[16];
    Kc[5] = -g_RT[1] - g_RT[9] + g_RT[12];
    Kc[6] = g_RT[0] + g_RT[12] - g_RT[14];
    Kc[7] = g_RT[0] + g_RT[10] - g_RT[12];
    Kc[8] = g_RT[0] + g_RT[9] - g_RT[10];
    Kc[9] = g_RT[0] + g_RT[2] - g_RT[4];
    Kc[10] = g_RT[5] + g_RT[11] - g_RT[22];
    Kc[11] = -g_RT[0] - g_RT[0] + g_RT[1];
    Kc[12] = g_RT[5] + g_RT[5] - g_RT[19];
    Kc[13] = g_RT[0] + g_RT[5] - g_RT[7];
    Kc[14] = g_RT[0] + g_RT[7] - g_RT[8];
    Kc[15] = -g_RT[0] - g_RT[11] + g_RT[13];
    Kc[16] = -g_RT[0] - g_RT[13] + g_RT[15];
    Kc[17] = -g_RT[1] - g_RT[11] + g_RT[15];
    Kc[18] = -g_RT[0] - g_RT[15] + g_RT[17];
    Kc[19] = -g_RT[0] - g_RT[15] + g_RT[18];
    Kc[20] = -g_RT[2] + g_RT[3];
    Kc[21] = -g_RT[8] - g_RT[11] + g_RT[25];
    Kc[22] = -g_RT[1] - g_RT[22] + g_RT[25];
    Kc[23] = g_RT[0] - g_RT[5] - g_RT[7] + g_RT[19];
    Kc[24] = -g_RT[0] + g_RT[1] + g_RT[5] - g_RT[7];
    Kc[25] = -g_RT[0] + g_RT[1] + g_RT[7] - g_RT[8];
    Kc[26] = g_RT[5] - g_RT[7] - g_RT[7] + g_RT[8];
    Kc[27] = g_RT[0] - g_RT[1] - g_RT[19] + g_RT[20];
    Kc[28] = g_RT[0] - g_RT[7] - g_RT[7] + g_RT[20];
    Kc[29] = g_RT[5] - g_RT[7] - g_RT[19] + g_RT[20];
    Kc[30] = g_RT[7] - g_RT[8] - g_RT[19] + g_RT[20];
    Kc[31] = -g_RT[19] + g_RT[20] + g_RT[20] - g_RT[21];
    Kc[32] = -g_RT[19] + g_RT[20] + g_RT[20] - g_RT[21];
    Kc[33] = g_RT[0] - g_RT[7] - g_RT[8] + g_RT[21];
    Kc[34] = g_RT[0] - g_RT[1] - g_RT[20] + g_RT[21];
    Kc[35] = g_RT[5] - g_RT[7] - g_RT[20] + g_RT[21];
    Kc[36] = g_RT[7] - g_RT[8] - g_RT[20] + g_RT[21];
    Kc[37] = g_RT[7] - g_RT[8] - g_RT[20] + g_RT[21];
    Kc[38] = -g_RT[5] + g_RT[11] + g_RT[19] - g_RT[22];
    Kc[39] = -g_RT[7] + g_RT[11] + g_RT[20] - g_RT[22];
    Kc[40] = -g_RT[0] + g_RT[7] + g_RT[11] - g_RT[22];
    Kc[41] = -g_RT[11] + g_RT[13] + g_RT[19] - g_RT[20];
    Kc[42] = g_RT[0] - g_RT[1] - g_RT[11] + g_RT[13];
    Kc[43] = g_RT[5] - g_RT[7] - g_RT[11] + g_RT[13];
    Kc[44] = g_RT[7] - g_RT[8] - g_RT[11] + g_RT[13];
    Kc[45] = -g_RT[0] + g_RT[5] + g_RT[13] - g_RT[22];
    Kc[46] = -g_RT[0] - g_RT[7] + g_RT[13] + g_RT[20] - g_RT[22];
    Kc[47] = -g_RT[1] - g_RT[11] - g_RT[11] + g_RT[13] + g_RT[13];
    Kc[48] = g_RT[4] - g_RT[6] - g_RT[11] + g_RT[13];
    Kc[49] = -g_RT[11] + g_RT[13] + g_RT[13] - g_RT[15];
    Kc[50] = g_RT[0] - g_RT[1] - g_RT[13] + g_RT[15];
    Kc[51] = g_RT[5] - g_RT[7] - g_RT[13] + g_RT[15];
    Kc[52] = g_RT[7] - g_RT[8] - g_RT[13] + g_RT[15];
    Kc[53] = -g_RT[13] + g_RT[15] + g_RT[19] - g_RT[20];
    Kc[54] = -g_RT[13] + g_RT[15] + g_RT[20] - g_RT[21];
    Kc[55] = g_RT[4] - g_RT[6] - g_RT[13] + g_RT[15];
    Kc[56] = -g_RT[0] + g_RT[4] + g_RT[5] - g_RT[15];
    Kc[57] = g_RT[4] - g_RT[5] - g_RT[18] + g_RT[19];
    Kc[58] = g_RT[4] - g_RT[7] - g_RT[15] + g_RT[19];
    Kc[59] = g_RT[4] - g_RT[7] - g_RT[18] + g_RT[20];
    Kc[60] = g_RT[0] - g_RT[1] - g_RT[4] + g_RT[6];
    Kc[61] = -g_RT[4] + g_RT[5] + g_RT[6] - g_RT[7];
    Kc[62] = -g_RT[4] + g_RT[6] + g_RT[7] - g_RT[8];
    Kc[63] = g_RT[4] - g_RT[6] - g_RT[19] + g_RT[20];
    Kc[64] = -g_RT[4] + g_RT[6] + g_RT[20] - g_RT[21];
    Kc[65] = g_RT[0] - g_RT[1] - g_RT[15] + g_RT[17];
    Kc[66] = g_RT[0] - g_RT[4] - g_RT[7] + g_RT[17];
    Kc[67] = g_RT[5] - g_RT[7] - g_RT[15] + g_RT[17];
    Kc[68] = g_RT[7] - g_RT[8] - g_RT[15] + g_RT[17];
    Kc[69] = -g_RT[15] + g_RT[17] + g_RT[19] - g_RT[20];
    Kc[70] = -g_RT[15] + g_RT[17] + g_RT[19] - g_RT[20];
    Kc[71] = -g_RT[15] + g_RT[17] + g_RT[20] - g_RT[21];
    Kc[72] = g_RT[13] - g_RT[15] - g_RT[15] + g_RT[17];
    Kc[73] = g_RT[0] - g_RT[4] - g_RT[7] + g_RT[18];
    Kc[74] = g_RT[5] - g_RT[7] - g_RT[15] + g_RT[18];
    Kc[75] = g_RT[7] - g_RT[8] - g_RT[15] + g_RT[18];
    Kc[76] = -g_RT[15] + g_RT[18] + g_RT[19] - g_RT[20];
    Kc[77] = -g_RT[15] + g_RT[18] + g_RT[19] - g_RT[20];
    Kc[78] = -g_RT[15] + g_RT[18] + g_RT[20] - g_RT[21];
    Kc[79] = -g_RT[4] + g_RT[11] + g_RT[18] - g_RT[22];
    Kc[80] = -g_RT[0] + g_RT[4] + g_RT[4] - g_RT[14];
    Kc[81] = g_RT[2] - g_RT[4] - g_RT[4] + g_RT[6];
    Kc[82] = g_RT[3] - g_RT[4] - g_RT[4] + g_RT[6];
    Kc[83] = -g_RT[2] + g_RT[4] + g_RT[7] - g_RT[8];
    Kc[84] = -g_RT[3] + g_RT[4] + g_RT[7] - g_RT[8];
    Kc[85] = -g_RT[0] + g_RT[2] + g_RT[4] - g_RT[12];
    Kc[86] = -g_RT[0] + g_RT[3] + g_RT[4] - g_RT[12];
    Kc[87] = g_RT[0] - g_RT[3] - g_RT[8] + g_RT[18];
    Kc[88] = g_RT[0] - g_RT[1] - g_RT[14] + g_RT[16];
    Kc[89] = g_RT[5] - g_RT[7] - g_RT[14] + g_RT[16];
    Kc[90] = g_RT[7] - g_RT[8] - g_RT[14] + g_RT[16];
    Kc[91] = -g_RT[14] + g_RT[16] + g_RT[19] - g_RT[20];
    Kc[92] = -g_RT[14] + g_RT[16] + g_RT[20] - g_RT[21];
    Kc[93] = g_RT[4] - g_RT[6] - g_RT[14] + g_RT[16];
    Kc[94] = g_RT[0] - g_RT[1] - g_RT[12] + g_RT[14];
    Kc[95] = -g_RT[4] + g_RT[5] + g_RT[14] - g_RT[15];
    Kc[96] = -g_RT[12] + g_RT[14] + g_RT[19] - g_RT[20];
    Kc[97] = -g_RT[12] + g_RT[14] + g_RT[14] - g_RT[16];
    Kc[98] = -g_RT[11] + g_RT[13] + g_RT[14] - g_RT[16];
    Kc[99] = -g_RT[0] + g_RT[5] + g_RT[14] - g_RT[23];
    Kc[100] = g_RT[0] - g_RT[1] - g_RT[10] + g_RT[12];
    Kc[101] = g_RT[7] - g_RT[8] - g_RT[10] + g_RT[12];
    Kc[102] = g_RT[4] - g_RT[6] - g_RT[10] + g_RT[12];
    Kc[103] = -g_RT[4] + g_RT[5] + g_RT[12] - g_RT[13];
    Kc[104] = g_RT[7] - g_RT[8] - g_RT[9] + g_RT[10];
    Kc[105] = g_RT[5] - g_RT[7] - g_RT[10] + g_RT[12];
    Kc[106] = -g_RT[10] + g_RT[12] + g_RT[19] - g_RT[20];
    Kc[107] = g_RT[0] - g_RT[1] - g_RT[9] + g_RT[10];
    Kc[108] = g_RT[10] - g_RT[12] - g_RT[20] + g_RT[21];
    Kc[109] = g_RT[4] - g_RT[6] - g_RT[9] + g_RT[10];
    Kc[110] = -g_RT[9] + g_RT[10] + g_RT[10] - g_RT[12];
    Kc[111] = g_RT[10] - g_RT[13] - g_RT[15] + g_RT[19];
    Kc[112] = -g_RT[9] + g_RT[10] + g_RT[19] - g_RT[20];
    Kc[113] = -g_RT[2] + g_RT[5] + g_RT[9] - g_RT[11];
    Kc[114] = -g_RT[4] + g_RT[7] + g_RT[9] - g_RT[11];
    Kc[115] = -g_RT[0] + g_RT[2] + g_RT[5] - g_RT[13];
    Kc[116] = -g_RT[0] + g_RT[2] + g_RT[7] - g_RT[15];
    Kc[117] = -g_RT[0] + g_RT[1] + g_RT[2] - g_RT[4];
    Kc[118] = g_RT[2] - g_RT[7] - g_RT[13] + g_RT[19];
    Kc[119] = g_RT[2] - g_RT[7] - g_RT[15] + g_RT[20];
    Kc[120] = -g_RT[1] + g_RT[2] + g_RT[2] - g_RT[9];
    Kc[121] = -g_RT[2] + g_RT[3] + g_RT[8] - g_RT[8];
    Kc[122] = -g_RT[2] + g_RT[3] + g_RT[11] - g_RT[11];
    Kc[123] = -g_RT[2] + g_RT[3] + g_RT[22] - g_RT[22];
    Kc[124] = -g_RT[1] + g_RT[3] + g_RT[5] - g_RT[11];
    Kc[125] = -g_RT[0] + g_RT[3] + g_RT[5] - g_RT[13];
    Kc[126] = -g_RT[0] + g_RT[3] + g_RT[7] - g_RT[15];
    Kc[127] = -g_RT[0] + g_RT[1] + g_RT[3] - g_RT[4];
    Kc[128] = -g_RT[0] + g_RT[3] - g_RT[7] - g_RT[11] + g_RT[19];
    Kc[129] = g_RT[3] - g_RT[8] - g_RT[11] + g_RT[19];
    Kc[130] = g_RT[3] - g_RT[11] - g_RT[15] + g_RT[22];
    Kc[131] = -g_RT[4] - g_RT[13] + g_RT[23];
    Kc[132] = -g_RT[4] - g_RT[18] + g_RT[26];
    Kc[133] = g_RT[7] - g_RT[8] - g_RT[24] + g_RT[26];
    Kc[134] = g_RT[0] - g_RT[1] - g_RT[24] + g_RT[26];
    Kc[135] = g_RT[4] - g_RT[6] - g_RT[24] + g_RT[26];
    Kc[136] = g_RT[5] - g_RT[7] - g_RT[24] + g_RT[26];
    Kc[137] = g_RT[20] - g_RT[21] - g_RT[24] + g_RT[26];
    Kc[138] = g_RT[19] - g_RT[20] - g_RT[24] + g_RT[26];
    Kc[139] = -g_RT[4] - g_RT[15] + g_RT[24];
    Kc[140] = -g_RT[15] + g_RT[18] + g_RT[24] - g_RT[26];
    Kc[141] = -g_RT[13] + g_RT[15] + g_RT[24] - g_RT[26];
    Kc[142] = -g_RT[7] + g_RT[20] + g_RT[24] - g_RT[30];
    Kc[143] = -g_RT[0] - g_RT[29] + g_RT[30];
    Kc[144] = g_RT[19] - g_RT[20] - g_RT[28] + g_RT[29];
    Kc[145] = g_RT[7] - g_RT[8] - g_RT[28] + g_RT[29];
    Kc[146] = g_RT[20] - g_RT[21] - g_RT[28] + g_RT[29];
    Kc[147] = g_RT[5] - g_RT[7] - g_RT[28] + g_RT[29];
    Kc[148] = g_RT[0] - g_RT[1] - g_RT[28] + g_RT[29];
    Kc[149] = g_RT[4] - g_RT[6] - g_RT[28] + g_RT[29];
    Kc[150] = -g_RT[11] - g_RT[18] + g_RT[28];
    Kc[151] = -g_RT[4] - g_RT[22] + g_RT[28];
    Kc[152] = g_RT[19] + g_RT[24] - g_RT[34];
    Kc[153] = -g_RT[19] - g_RT[30] - g_RT[30] + g_RT[34] + g_RT[34];
    Kc[154] = -g_RT[19] - g_RT[29] - g_RT[31] + g_RT[34] + g_RT[34];
    Kc[155] = -g_RT[15] - g_RT[18] + g_RT[30];
    Kc[156] = g_RT[19] - g_RT[20] - g_RT[29] + g_RT[30];
    Kc[157] = g_RT[34] - g_RT[35];
    Kc[158] = -g_RT[7] - g_RT[15] - g_RT[15] + g_RT[35];
    Kc[159] = g_RT[19] + g_RT[35] - g_RT[37];
    Kc[160] = -g_RT[7] - g_RT[36] + g_RT[37];
    Kc[161] = -g_RT[7] - g_RT[32] + g_RT[36];
    Kc[162] = g_RT[32] - g_RT[33];
    Kc[163] = -g_RT[11] - g_RT[27] + g_RT[33];
    Kc[164] = -g_RT[17] - g_RT[22] + g_RT[33];
    Kc[165] = -g_RT[0] - g_RT[25] + g_RT[27];
    Kc[166] = g_RT[7] + g_RT[15] - g_RT[27];
    Kc[167] = -g_RT[7] - g_RT[13] + g_RT[25];
    Kc[168] = -g_RT[0] + g_RT[7] - g_RT[8] - g_RT[22] + g_RT[25];
    Kc[169] = g_RT[7] - g_RT[7] - g_RT[8] - g_RT[11] + g_RT[25];
    Kc[170] = g_RT[0] - g_RT[0] - g_RT[1] - g_RT[22] + g_RT[25];
    Kc[171] = g_RT[0] - g_RT[1] - g_RT[7] - g_RT[11] + g_RT[25];
    Kc[172] = g_RT[4] - g_RT[6] - g_RT[7] - g_RT[11] + g_RT[25];
    Kc[173] = -g_RT[7] - g_RT[11] + g_RT[20] - g_RT[21] + g_RT[25];
    Kc[174] = g_RT[5] - g_RT[7] - g_RT[7] - g_RT[11] + g_RT[25];

    for (int i=0; i<175; ++i) {
        Kc[i] = exp(Kc[i]);
    };

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    amrex::Real refC = 101325 / 8.31446 * invT;
    amrex::Real refCinv = 1 / refC;

    Kc[0] *= refCinv;
    Kc[1] *= refC;
    Kc[2] *= refCinv;
    Kc[3] *= refCinv;
    Kc[4] *= refCinv;
    Kc[5] *= refC;
    Kc[6] *= refCinv;
    Kc[7] *= refCinv;
    Kc[8] *= refCinv;
    Kc[9] *= refCinv;
    Kc[10] *= refCinv;
    Kc[11] *= refC;
    Kc[12] *= refCinv;
    Kc[13] *= refCinv;
    Kc[14] *= refCinv;
    Kc[15] *= refC;
    Kc[16] *= refC;
    Kc[17] *= refC;
    Kc[18] *= refC;
    Kc[19] *= refC;
    Kc[21] *= refC;
    Kc[22] *= refC;
    Kc[46] *= refC;
    Kc[47] *= refC;
    Kc[128] *= refC;
    Kc[131] *= refC;
    Kc[132] *= refC;
    Kc[139] *= refC;
    Kc[143] *= refC;
    Kc[150] *= refC;
    Kc[151] *= refC;
    Kc[152] *= refCinv;
    Kc[153] *= refC;
    Kc[154] *= refC;
    Kc[155] *= refC;
    Kc[158] *= pow(refC,2.000000);
    Kc[159] *= refCinv;
    Kc[160] *= refC;
    Kc[161] *= refC;
    Kc[163] *= refC;
    Kc[164] *= refC;
    Kc[165] *= refC;
    Kc[166] *= refCinv;
    Kc[167] *= refC;
    Kc[168] *= refC;
    Kc[169] *= refC;
    Kc[170] *= refC;
    Kc[171] *= refC;
    Kc[172] *= refC;
    Kc[173] *= refC;
    Kc[174] *= refC;

    return;
}

void comp_qfqr(amrex::Real *  qf, amrex::Real *  qr, amrex::Real *  sc, amrex::Real * qss_sc, amrex::Real *  tc, amrex::Real invT)
{

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    qf[0] = sc[0]*sc[19];
    qr[0] = sc[20];

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    qf[1] = sc[21];
    qr[1] = sc[7]*sc[7];

    /*reaction 3: CH3 + CH3 (+M) <=> C2H6 (+M) */
    qf[2] = sc[4]*sc[4];
    qr[2] = sc[16];

    /*reaction 4: CH3 + H (+M) <=> CH4 (+M) */
    qf[3] = sc[0]*sc[4];
    qr[3] = sc[6];

    /*reaction 5: C2H5 + H (+M) <=> C2H6 (+M) */
    qf[4] = sc[0]*sc[14];
    qr[4] = sc[16];

    /*reaction 6: C2H4 (+M) <=> H2 + C2H2 (+M) */
    qf[5] = sc[12];
    qr[5] = sc[1]*sc[9];

    /*reaction 7: C2H4 + H (+M) <=> C2H5 (+M) */
    qf[6] = sc[0]*sc[12];
    qr[6] = sc[14];

    /*reaction 8: C2H3 + H (+M) <=> C2H4 (+M) */
    qf[7] = sc[0]*sc[10];
    qr[7] = sc[12];

    /*reaction 9: C2H2 + H (+M) <=> C2H3 (+M) */
    qf[8] = sc[0]*sc[9];
    qr[8] = sc[10];

    /*reaction 10: CH2 + H (+M) <=> CH3 (+M) */
    qf[9] = sc[0]*sc[2];
    qr[9] = sc[4];

    /*reaction 11: CO + O (+M) <=> CO2 (+M) */
    qf[10] = sc[5]*sc[11];
    qr[10] = sc[22];

    /*reaction 12: H2 + M <=> H + H + M */
    qf[11] = sc[1];
    qr[11] = sc[0]*sc[0];

    /*reaction 13: O + O + M <=> O2 + M */
    qf[12] = sc[5]*sc[5];
    qr[12] = sc[19];

    /*reaction 14: O + H + M <=> OH + M */
    qf[13] = sc[0]*sc[5];
    qr[13] = sc[7];

    /*reaction 15: H + OH + M <=> H2O + M */
    qf[14] = sc[0]*sc[7];
    qr[14] = sc[8];

    /*reaction 16: HCO + M <=> H + CO + M */
    qf[15] = sc[13];
    qr[15] = sc[0]*sc[11];

    /*reaction 17: CH2O + M <=> HCO + H + M */
    qf[16] = sc[15];
    qr[16] = sc[0]*sc[13];

    /*reaction 18: CH2O + M <=> CO + H2 + M */
    qf[17] = sc[15];
    qr[17] = sc[1]*sc[11];

    /*reaction 19: CH2OH + M <=> CH2O + H + M */
    qf[18] = sc[17];
    qr[18] = sc[0]*sc[15];

    /*reaction 20: CH3O + M <=> CH2O + H + M */
    qf[19] = sc[18];
    qr[19] = sc[0]*sc[15];

    /*reaction 21: CH2(S) + M <=> CH2 + M */
    qf[20] = sc[3];
    qr[20] = sc[2];

    /*reaction 22: HCOOH + M <=> CO + H2O + M */
    qf[21] = sc[25];
    qr[21] = sc[8]*sc[11];

    /*reaction 23: HCOOH + M <=> CO2 + H2 + M */
    qf[22] = sc[25];
    qr[22] = sc[1]*sc[22];

    /*reaction 24: H + O2 <=> O + OH */
    qf[23] = sc[0]*sc[19];
    qr[23] = sc[5]*sc[7];

    /*reaction 25: O + H2 <=> H + OH */
    qf[24] = sc[1]*sc[5];
    qr[24] = sc[0]*sc[7];

    /*reaction 26: H2 + OH <=> H2O + H */
    qf[25] = sc[1]*sc[7];
    qr[25] = sc[0]*sc[8];

    /*reaction 27: O + H2O <=> OH + OH */
    qf[26] = sc[5]*sc[8];
    qr[26] = sc[7]*sc[7];

    /*reaction 28: HO2 + H <=> H2 + O2 */
    qf[27] = sc[0]*sc[20];
    qr[27] = sc[1]*sc[19];

    /*reaction 29: HO2 + H <=> OH + OH */
    qf[28] = sc[0]*sc[20];
    qr[28] = sc[7]*sc[7];

    /*reaction 30: HO2 + O <=> O2 + OH */
    qf[29] = sc[5]*sc[20];
    qr[29] = sc[7]*sc[19];

    /*reaction 31: HO2 + OH <=> H2O + O2 */
    qf[30] = sc[7]*sc[20];
    qr[30] = sc[8]*sc[19];

    /*reaction 32: HO2 + HO2 <=> H2O2 + O2 */
    qf[31] = sc[20]*sc[20];
    qr[31] = sc[19]*sc[21];

    /*reaction 33: HO2 + HO2 <=> H2O2 + O2 */
    qf[32] = sc[20]*sc[20];
    qr[32] = sc[19]*sc[21];

    /*reaction 34: H2O2 + H <=> H2O + OH */
    qf[33] = sc[0]*sc[21];
    qr[33] = sc[7]*sc[8];

    /*reaction 35: H2O2 + H <=> HO2 + H2 */
    qf[34] = sc[0]*sc[21];
    qr[34] = sc[1]*sc[20];

    /*reaction 36: H2O2 + O <=> OH + HO2 */
    qf[35] = sc[5]*sc[21];
    qr[35] = sc[7]*sc[20];

    /*reaction 37: H2O2 + OH <=> HO2 + H2O */
    qf[36] = sc[7]*sc[21];
    qr[36] = sc[8]*sc[20];

    /*reaction 38: H2O2 + OH <=> HO2 + H2O */
    qf[37] = sc[7]*sc[21];
    qr[37] = sc[8]*sc[20];

    /*reaction 39: CO + O2 <=> CO2 + O */
    qf[38] = sc[11]*sc[19];
    qr[38] = sc[5]*sc[22];

    /*reaction 40: CO + HO2 <=> CO2 + OH */
    qf[39] = sc[11]*sc[20];
    qr[39] = sc[7]*sc[22];

    /*reaction 41: CO + OH <=> CO2 + H */
    qf[40] = sc[7]*sc[11];
    qr[40] = sc[0]*sc[22];

    /*reaction 42: HCO + O2 <=> CO + HO2 */
    qf[41] = sc[13]*sc[19];
    qr[41] = sc[11]*sc[20];

    /*reaction 43: HCO + H <=> CO + H2 */
    qf[42] = sc[0]*sc[13];
    qr[42] = sc[1]*sc[11];

    /*reaction 44: HCO + O <=> CO + OH */
    qf[43] = sc[5]*sc[13];
    qr[43] = sc[7]*sc[11];

    /*reaction 45: HCO + OH <=> CO + H2O */
    qf[44] = sc[7]*sc[13];
    qr[44] = sc[8]*sc[11];

    /*reaction 46: HCO + O <=> CO2 + H */
    qf[45] = sc[5]*sc[13];
    qr[45] = sc[0]*sc[22];

    /*reaction 47: HCO + HO2 <=> CO2 + OH + H */
    qf[46] = sc[13]*sc[20];
    qr[46] = sc[0]*sc[7]*sc[22];

    /*reaction 48: HCO + HCO <=> H2 + CO + CO */
    qf[47] = sc[13]*sc[13];
    qr[47] = sc[1]*sc[11]*sc[11];

    /*reaction 49: HCO + CH3 <=> CO + CH4 */
    qf[48] = sc[4]*sc[13];
    qr[48] = sc[6]*sc[11];

    /*reaction 50: HCO + HCO <=> CH2O + CO */
    qf[49] = sc[13]*sc[13];
    qr[49] = sc[11]*sc[15];

    /*reaction 51: CH2O + H <=> HCO + H2 */
    qf[50] = sc[0]*sc[15];
    qr[50] = sc[1]*sc[13];

    /*reaction 52: CH2O + O <=> HCO + OH */
    qf[51] = sc[5]*sc[15];
    qr[51] = sc[7]*sc[13];

    /*reaction 53: CH2O + OH <=> HCO + H2O */
    qf[52] = sc[7]*sc[15];
    qr[52] = sc[8]*sc[13];

    /*reaction 54: CH2O + O2 <=> HCO + HO2 */
    qf[53] = sc[15]*sc[19];
    qr[53] = sc[13]*sc[20];

    /*reaction 55: CH2O + HO2 <=> HCO + H2O2 */
    qf[54] = sc[15]*sc[20];
    qr[54] = sc[13]*sc[21];

    /*reaction 56: CH2O + CH3 <=> HCO + CH4 */
    qf[55] = sc[4]*sc[15];
    qr[55] = sc[6]*sc[13];

    /*reaction 57: CH3 + O <=> CH2O + H */
    qf[56] = sc[4]*sc[5];
    qr[56] = sc[0]*sc[15];

    /*reaction 58: CH3 + O2 <=> CH3O + O */
    qf[57] = sc[4]*sc[19];
    qr[57] = sc[5]*sc[18];

    /*reaction 59: CH3 + O2 <=> CH2O + OH */
    qf[58] = sc[4]*sc[19];
    qr[58] = sc[7]*sc[15];

    /*reaction 60: CH3 + HO2 <=> CH3O + OH */
    qf[59] = sc[4]*sc[20];
    qr[59] = sc[7]*sc[18];

    /*reaction 61: CH4 + H <=> CH3 + H2 */
    qf[60] = sc[0]*sc[6];
    qr[60] = sc[1]*sc[4];

    /*reaction 62: CH4 + O <=> CH3 + OH */
    qf[61] = sc[5]*sc[6];
    qr[61] = sc[4]*sc[7];

    /*reaction 63: CH4 + OH <=> CH3 + H2O */
    qf[62] = sc[6]*sc[7];
    qr[62] = sc[4]*sc[8];

    /*reaction 64: CH3 + HO2 <=> CH4 + O2 */
    qf[63] = sc[4]*sc[20];
    qr[63] = sc[6]*sc[19];

    /*reaction 65: CH4 + HO2 <=> CH3 + H2O2 */
    qf[64] = sc[6]*sc[20];
    qr[64] = sc[4]*sc[21];

    /*reaction 66: CH2OH + H <=> CH2O + H2 */
    qf[65] = sc[0]*sc[17];
    qr[65] = sc[1]*sc[15];

    /*reaction 67: CH2OH + H <=> CH3 + OH */
    qf[66] = sc[0]*sc[17];
    qr[66] = sc[4]*sc[7];

    /*reaction 68: CH2OH + O <=> CH2O + OH */
    qf[67] = sc[5]*sc[17];
    qr[67] = sc[7]*sc[15];

    /*reaction 69: CH2OH + OH <=> CH2O + H2O */
    qf[68] = sc[7]*sc[17];
    qr[68] = sc[8]*sc[15];

    /*reaction 70: CH2OH + O2 <=> CH2O + HO2 */
    qf[69] = sc[17]*sc[19];
    qr[69] = sc[15]*sc[20];

    /*reaction 71: CH2OH + O2 <=> CH2O + HO2 */
    qf[70] = sc[17]*sc[19];
    qr[70] = sc[15]*sc[20];

    /*reaction 72: CH2OH + HO2 <=> CH2O + H2O2 */
    qf[71] = sc[17]*sc[20];
    qr[71] = sc[15]*sc[21];

    /*reaction 73: CH2OH + HCO <=> CH2O + CH2O */
    qf[72] = sc[13]*sc[17];
    qr[72] = sc[15]*sc[15];

    /*reaction 74: CH3O + H <=> CH3 + OH */
    qf[73] = sc[0]*sc[18];
    qr[73] = sc[4]*sc[7];

    /*reaction 75: CH3O + O <=> CH2O + OH */
    qf[74] = sc[5]*sc[18];
    qr[74] = sc[7]*sc[15];

    /*reaction 76: CH3O + OH <=> CH2O + H2O */
    qf[75] = sc[7]*sc[18];
    qr[75] = sc[8]*sc[15];

    /*reaction 77: CH3O + O2 <=> CH2O + HO2 */
    qf[76] = sc[18]*sc[19];
    qr[76] = sc[15]*sc[20];

    /*reaction 78: CH3O + O2 <=> CH2O + HO2 */
    qf[77] = sc[18]*sc[19];
    qr[77] = sc[15]*sc[20];

    /*reaction 79: CH3O + HO2 <=> CH2O + H2O2 */
    qf[78] = sc[18]*sc[20];
    qr[78] = sc[15]*sc[21];

    /*reaction 80: CH3O + CO <=> CH3 + CO2 */
    qf[79] = sc[11]*sc[18];
    qr[79] = sc[4]*sc[22];

    /*reaction 81: CH3 + CH3 <=> H + C2H5 */
    qf[80] = sc[4]*sc[4];
    qr[80] = sc[0]*sc[14];

    /*reaction 82: CH4 + CH2 <=> CH3 + CH3 */
    qf[81] = sc[2]*sc[6];
    qr[81] = sc[4]*sc[4];

    /*reaction 83: CH4 + CH2(S) <=> CH3 + CH3 */
    qf[82] = sc[3]*sc[6];
    qr[82] = sc[4]*sc[4];

    /*reaction 84: CH3 + OH <=> CH2 + H2O */
    qf[83] = sc[4]*sc[7];
    qr[83] = sc[2]*sc[8];

    /*reaction 85: CH3 + OH <=> CH2(S) + H2O */
    qf[84] = sc[4]*sc[7];
    qr[84] = sc[3]*sc[8];

    /*reaction 86: CH3 + CH2 <=> C2H4 + H */
    qf[85] = sc[2]*sc[4];
    qr[85] = sc[0]*sc[12];

    /*reaction 87: CH3 + CH2(S) <=> C2H4 + H */
    qf[86] = sc[3]*sc[4];
    qr[86] = sc[0]*sc[12];

    /*reaction 88: CH3O + H <=> CH2(S) + H2O */
    qf[87] = sc[0]*sc[18];
    qr[87] = sc[3]*sc[8];

    /*reaction 89: C2H6 + H <=> C2H5 + H2 */
    qf[88] = sc[0]*sc[16];
    qr[88] = sc[1]*sc[14];

    /*reaction 90: C2H6 + O <=> C2H5 + OH */
    qf[89] = sc[5]*sc[16];
    qr[89] = sc[7]*sc[14];

    /*reaction 91: C2H6 + OH <=> C2H5 + H2O */
    qf[90] = sc[7]*sc[16];
    qr[90] = sc[8]*sc[14];

    /*reaction 92: C2H6 + O2 <=> C2H5 + HO2 */
    qf[91] = sc[16]*sc[19];
    qr[91] = sc[14]*sc[20];

    /*reaction 93: C2H6 + HO2 <=> C2H5 + H2O2 */
    qf[92] = sc[16]*sc[20];
    qr[92] = sc[14]*sc[21];

    /*reaction 94: C2H6 + CH3 <=> C2H5 + CH4 */
    qf[93] = sc[4]*sc[16];
    qr[93] = sc[6]*sc[14];

    /*reaction 95: C2H5 + H <=> C2H4 + H2 */
    qf[94] = sc[0]*sc[14];
    qr[94] = sc[1]*sc[12];

    /*reaction 96: C2H5 + O <=> CH3 + CH2O */
    qf[95] = sc[5]*sc[14];
    qr[95] = sc[4]*sc[15];

    /*reaction 97: C2H5 + O2 <=> C2H4 + HO2 */
    qf[96] = sc[14]*sc[19];
    qr[96] = sc[12]*sc[20];

    /*reaction 98: C2H5 + C2H5 <=> C2H4 + C2H6 */
    qf[97] = sc[14]*sc[14];
    qr[97] = sc[12]*sc[16];

    /*reaction 99: C2H5 + HCO <=> C2H6 + CO */
    qf[98] = sc[13]*sc[14];
    qr[98] = sc[11]*sc[16];

    /*reaction 100: C2H5 + O <=> CH3HCO + H */
    qf[99] = sc[5]*sc[14];
    qr[99] = sc[0]*sc[23];

    /*reaction 101: C2H4 + H <=> C2H3 + H2 */
    qf[100] = sc[0]*sc[12];
    qr[100] = sc[1]*sc[10];

    /*reaction 102: C2H4 + OH <=> C2H3 + H2O */
    qf[101] = sc[7]*sc[12];
    qr[101] = sc[8]*sc[10];

    /*reaction 103: C2H4 + CH3 <=> C2H3 + CH4 */
    qf[102] = sc[4]*sc[12];
    qr[102] = sc[6]*sc[10];

    /*reaction 104: C2H4 + O <=> CH3 + HCO */
    qf[103] = sc[5]*sc[12];
    qr[103] = sc[4]*sc[13];

    /*reaction 105: C2H3 + OH <=> C2H2 + H2O */
    qf[104] = sc[7]*sc[10];
    qr[104] = sc[8]*sc[9];

    /*reaction 106: C2H4 + O <=> OH + C2H3 */
    qf[105] = sc[5]*sc[12];
    qr[105] = sc[7]*sc[10];

    /*reaction 107: C2H4 + O2 <=> C2H3 + HO2 */
    qf[106] = sc[12]*sc[19];
    qr[106] = sc[10]*sc[20];

    /*reaction 108: C2H3 + H <=> C2H2 + H2 */
    qf[107] = sc[0]*sc[10];
    qr[107] = sc[1]*sc[9];

    /*reaction 109: C2H3 + H2O2 <=> C2H4 + HO2 */
    qf[108] = sc[10]*sc[21];
    qr[108] = sc[12]*sc[20];

    /*reaction 110: C2H3 + CH3 <=> C2H2 + CH4 */
    qf[109] = sc[4]*sc[10];
    qr[109] = sc[6]*sc[9];

    /*reaction 111: C2H3 + C2H3 <=> C2H4 + C2H2 */
    qf[110] = sc[10]*sc[10];
    qr[110] = sc[9]*sc[12];

    /*reaction 112: C2H3 + O2 <=> HCO + CH2O */
    qf[111] = sc[10]*sc[19];
    qr[111] = sc[13]*sc[15];

    /*reaction 113: C2H3 + O2 <=> HO2 + C2H2 */
    qf[112] = sc[10]*sc[19];
    qr[112] = sc[9]*sc[20];

    /*reaction 114: C2H2 + O <=> CH2 + CO */
    qf[113] = sc[5]*sc[9];
    qr[113] = sc[2]*sc[11];

    /*reaction 115: C2H2 + OH <=> CH3 + CO */
    qf[114] = sc[7]*sc[9];
    qr[114] = sc[4]*sc[11];

    /*reaction 116: CH2 + O <=> HCO + H */
    qf[115] = sc[2]*sc[5];
    qr[115] = sc[0]*sc[13];

    /*reaction 117: CH2 + OH <=> CH2O + H */
    qf[116] = sc[2]*sc[7];
    qr[116] = sc[0]*sc[15];

    /*reaction 118: CH2 + H2 <=> H + CH3 */
    qf[117] = sc[1]*sc[2];
    qr[117] = sc[0]*sc[4];

    /*reaction 119: CH2 + O2 <=> HCO + OH */
    qf[118] = sc[2]*sc[19];
    qr[118] = sc[7]*sc[13];

    /*reaction 120: CH2 + HO2 <=> CH2O + OH */
    qf[119] = sc[2]*sc[20];
    qr[119] = sc[7]*sc[15];

    /*reaction 121: CH2 + CH2 <=> C2H2 + H2 */
    qf[120] = sc[2]*sc[2];
    qr[120] = sc[1]*sc[9];

    /*reaction 122: CH2(S) + H2O <=> CH2 + H2O */
    qf[121] = sc[3]*sc[8];
    qr[121] = sc[2]*sc[8];

    /*reaction 123: CH2(S) + CO <=> CH2 + CO */
    qf[122] = sc[3]*sc[11];
    qr[122] = sc[2]*sc[11];

    /*reaction 124: CH2(S) + CO2 <=> CH2 + CO2 */
    qf[123] = sc[3]*sc[22];
    qr[123] = sc[2]*sc[22];

    /*reaction 125: CH2(S) + O <=> CO + H2 */
    qf[124] = sc[3]*sc[5];
    qr[124] = sc[1]*sc[11];

    /*reaction 126: CH2(S) + O <=> HCO + H */
    qf[125] = sc[3]*sc[5];
    qr[125] = sc[0]*sc[13];

    /*reaction 127: CH2(S) + OH <=> CH2O + H */
    qf[126] = sc[3]*sc[7];
    qr[126] = sc[0]*sc[15];

    /*reaction 128: CH2(S) + H2 <=> CH3 + H */
    qf[127] = sc[1]*sc[3];
    qr[127] = sc[0]*sc[4];

    /*reaction 129: CH2(S) + O2 <=> H + OH + CO */
    qf[128] = sc[3]*sc[19];
    qr[128] = sc[0]*sc[7]*sc[11];

    /*reaction 130: CH2(S) + O2 <=> CO + H2O */
    qf[129] = sc[3]*sc[19];
    qr[129] = sc[8]*sc[11];

    /*reaction 131: CH2(S) + CO2 <=> CH2O + CO */
    qf[130] = sc[3]*sc[22];
    qr[130] = sc[11]*sc[15];

    /*reaction 132: CH3HCO <=> CH3 + HCO */
    qf[131] = sc[23];
    qr[131] = sc[4]*sc[13];

    /*reaction 133: CH3OCH3 <=> CH3 + CH3O */
    qf[132] = sc[26];
    qr[132] = sc[4]*sc[18];

    /*reaction 134: CH3OCH3 + OH <=> CH3OCH2 + H2O */
    qf[133] = sc[7]*sc[26];
    qr[133] = sc[8]*sc[24];

    /*reaction 135: CH3OCH3 + H <=> CH3OCH2 + H2 */
    qf[134] = sc[0]*sc[26];
    qr[134] = sc[1]*sc[24];

    /*reaction 136: CH3OCH3 + CH3 <=> CH3OCH2 + CH4 */
    qf[135] = sc[4]*sc[26];
    qr[135] = sc[6]*sc[24];

    /*reaction 137: CH3OCH3 + O <=> CH3OCH2 + OH */
    qf[136] = sc[5]*sc[26];
    qr[136] = sc[7]*sc[24];

    /*reaction 138: CH3OCH3 + HO2 <=> CH3OCH2 + H2O2 */
    qf[137] = sc[20]*sc[26];
    qr[137] = sc[21]*sc[24];

    /*reaction 139: CH3OCH3 + O2 <=> CH3OCH2 + HO2 */
    qf[138] = sc[19]*sc[26];
    qr[138] = sc[20]*sc[24];

    /*reaction 140: CH3OCH2 <=> CH2O + CH3 */
    qf[139] = sc[24];
    qr[139] = sc[4]*sc[15];

    /*reaction 141: CH3OCH2 + CH3O <=> CH3OCH3 + CH2O */
    qf[140] = sc[18]*sc[24];
    qr[140] = sc[15]*sc[26];

    /*reaction 142: CH3OCH2 + CH2O <=> CH3OCH3 + HCO */
    qf[141] = sc[15]*sc[24];
    qr[141] = sc[13]*sc[26];

    /*reaction 143: CH3OCH2 + HO2 <=> CH3OCH2O + OH */
    qf[142] = sc[20]*sc[24];
    qr[142] = sc[7]*sc[30];

    /*reaction 144: CH3OCH2O <=> CH3OCHO + H */
    qf[143] = sc[30];
    qr[143] = sc[0]*sc[29];

    /*reaction 145: CH3OCHO + O2 <=> CH3OCO + HO2 */
    qf[144] = sc[19]*sc[29];
    qr[144] = sc[20]*sc[28];

    /*reaction 146: CH3OCHO + OH <=> CH3OCO + H2O */
    qf[145] = sc[7]*sc[29];
    qr[145] = sc[8]*sc[28];

    /*reaction 147: CH3OCHO + HO2 <=> CH3OCO + H2O2 */
    qf[146] = sc[20]*sc[29];
    qr[146] = sc[21]*sc[28];

    /*reaction 148: CH3OCHO + O <=> CH3OCO + OH */
    qf[147] = sc[5]*sc[29];
    qr[147] = sc[7]*sc[28];

    /*reaction 149: CH3OCHO + H <=> CH3OCO + H2 */
    qf[148] = sc[0]*sc[29];
    qr[148] = sc[1]*sc[28];

    /*reaction 150: CH3OCHO + CH3 <=> CH3OCO + CH4 */
    qf[149] = sc[4]*sc[29];
    qr[149] = sc[6]*sc[28];

    /*reaction 151: CH3OCO <=> CH3O + CO */
    qf[150] = sc[28];
    qr[150] = sc[11]*sc[18];

    /*reaction 152: CH3OCO <=> CH3 + CO2 */
    qf[151] = sc[28];
    qr[151] = sc[4]*sc[22];

    /*reaction 153: CH3OCH2 + O2 <=> CH3OCH2O2 */
    qf[152] = sc[19]*sc[24];
    qr[152] = sc[34];

    /*reaction 154: CH3OCH2O2 + CH3OCH2O2 <=> O2 + CH3OCH2O + CH3OCH2O */
    qf[153] = sc[34]*sc[34];
    qr[153] = sc[19]*sc[30]*sc[30];

    /*reaction 155: CH3OCH2O2 + CH3OCH2O2 <=> O2 + CH3OCHO + CH3OCH2OH */
    qf[154] = sc[34]*sc[34];
    qr[154] = sc[19]*sc[29]*sc[31];

    /*reaction 156: CH3OCH2O <=> CH3O + CH2O */
    qf[155] = sc[30];
    qr[155] = sc[15]*sc[18];

    /*reaction 157: CH3OCH2O + O2 <=> CH3OCHO + HO2 */
    qf[156] = sc[19]*sc[30];
    qr[156] = sc[20]*sc[29];

    /*reaction 158: CH3OCH2O2 <=> CH2OCH2O2H */
    qf[157] = sc[34];
    qr[157] = sc[35];

    /*reaction 159: CH2OCH2O2H <=> OH + CH2O + CH2O */
    qf[158] = sc[35];
    qr[158] = sc[7]*sc[15]*sc[15];

    /*reaction 160: CH2OCH2O2H + O2 <=> O2CH2OCH2O2H */
    qf[159] = sc[19]*sc[35];
    qr[159] = sc[37];

    /*reaction 161: O2CH2OCH2O2H <=> HO2CH2OCHO + OH */
    qf[160] = sc[37];
    qr[160] = sc[7]*sc[36];

    /*reaction 162: HO2CH2OCHO <=> OCH2OCHO + OH */
    qf[161] = sc[36];
    qr[161] = sc[7]*sc[32];

    /*reaction 163: OCH2OCHO <=> HOCH2OCO */
    qf[162] = sc[32];
    qr[162] = sc[33];

    /*reaction 164: HOCH2OCO <=> HOCH2O + CO */
    qf[163] = sc[33];
    qr[163] = sc[11]*sc[27];

    /*reaction 165: HOCH2OCO <=> CH2OH + CO2 */
    qf[164] = sc[33];
    qr[164] = sc[17]*sc[22];

    /*reaction 166: HOCH2O <=> HCOOH + H */
    qf[165] = sc[27];
    qr[165] = sc[0]*sc[25];

    /*reaction 167: CH2O + OH <=> HOCH2O */
    qf[166] = sc[7]*sc[15];
    qr[166] = sc[27];

    /*reaction 168: HCOOH <=> HCO + OH */
    qf[167] = sc[25];
    qr[167] = sc[7]*sc[13];

    /*reaction 169: HCOOH + OH <=> H2O + CO2 + H */
    qf[168] = sc[7]*sc[25];
    qr[168] = sc[0]*sc[8]*sc[22];

    /*reaction 170: HCOOH + OH <=> H2O + CO + OH */
    qf[169] = sc[7]*sc[25];
    qr[169] = sc[7]*sc[8]*sc[11];

    /*reaction 171: HCOOH + H <=> H2 + CO2 + H */
    qf[170] = sc[0]*sc[25];
    qr[170] = sc[0]*sc[1]*sc[22];

    /*reaction 172: HCOOH + H <=> H2 + CO + OH */
    qf[171] = sc[0]*sc[25];
    qr[171] = sc[1]*sc[7]*sc[11];

    /*reaction 173: HCOOH + CH3 <=> CH4 + CO + OH */
    qf[172] = sc[4]*sc[25];
    qr[172] = sc[6]*sc[7]*sc[11];

    /*reaction 174: HCOOH + HO2 <=> H2O2 + CO + OH */
    qf[173] = sc[20]*sc[25];
    qr[173] = sc[7]*sc[11]*sc[21];

    /*reaction 175: HCOOH + O <=> CO + OH + OH */
    qf[174] = sc[5]*sc[25];
    qr[174] = sc[7]*sc[7]*sc[11];

    amrex::Real T = tc[1];

    /*compute the mixture concentration */
    amrex::Real mixture = 0.0;
    for (int i = 0; i < 39; ++i) {
        mixture += sc[i];
    }

    amrex::Real Corr[175];
    for (int i = 0; i < 175; ++i) {
        Corr[i] = 1.0;
    }

    /* troe */
    {
        amrex::Real alpha[10];
        alpha[0] = mixture + (TB[0][0] - 1)*sc[1] + (TB[0][1] - 1)*sc[8] + (TB[0][2] - 1)*sc[19] + (TB[0][3] - 1)*sc[11] + (TB[0][4] - 1)*sc[22];
        alpha[1] = mixture + (TB[1][0] - 1)*sc[1] + (TB[1][1] - 1)*sc[8] + (TB[1][2] - 1)*sc[11] + (TB[1][3] - 1)*sc[22];
        alpha[2] = mixture + (TB[2][0] - 1)*sc[8] + (TB[2][1] - 1)*sc[11] + (TB[2][2] - 1)*sc[22];
        alpha[3] = mixture + (TB[3][0] - 1)*sc[1] + (TB[3][1] - 1)*sc[8] + (TB[3][2] - 1)*sc[6] + (TB[3][3] - 1)*sc[11] + (TB[3][4] - 1)*sc[22] + (TB[3][5] - 1)*sc[16];
        alpha[4] = mixture + (TB[4][0] - 1)*sc[1] + (TB[4][1] - 1)*sc[8] + (TB[4][2] - 1)*sc[6] + (TB[4][3] - 1)*sc[11] + (TB[4][4] - 1)*sc[22] + (TB[4][5] - 1)*sc[16];
        alpha[5] = mixture + (TB[5][0] - 1)*sc[1] + (TB[5][1] - 1)*sc[8] + (TB[5][2] - 1)*sc[6] + (TB[5][3] - 1)*sc[11] + (TB[5][4] - 1)*sc[22] + (TB[5][5] - 1)*sc[16];
        alpha[6] = mixture + (TB[6][0] - 1)*sc[1] + (TB[6][1] - 1)*sc[8] + (TB[6][2] - 1)*sc[6] + (TB[6][3] - 1)*sc[11] + (TB[6][4] - 1)*sc[22] + (TB[6][5] - 1)*sc[16];
        alpha[7] = mixture + (TB[7][0] - 1)*sc[1] + (TB[7][1] - 1)*sc[8] + (TB[7][2] - 1)*sc[6] + (TB[7][3] - 1)*sc[11] + (TB[7][4] - 1)*sc[22] + (TB[7][5] - 1)*sc[16];
        alpha[8] = mixture + (TB[8][0] - 1)*sc[1] + (TB[8][1] - 1)*sc[8] + (TB[8][2] - 1)*sc[6] + (TB[8][3] - 1)*sc[11] + (TB[8][4] - 1)*sc[22] + (TB[8][5] - 1)*sc[16];
        alpha[9] = mixture + (TB[9][0] - 1)*sc[1] + (TB[9][1] - 1)*sc[8] + (TB[9][2] - 1)*sc[6] + (TB[9][3] - 1)*sc[11] + (TB[9][4] - 1)*sc[22] + (TB[9][5] - 1)*sc[16];
        for (int i=0; i<10; i++)
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
        alpha = mixture + (TB[10][0] - 1)*sc[1] + (TB[10][1] - 1)*sc[8] + (TB[10][2] - 1)*sc[11] + (TB[10][3] - 1)*sc[22];
        amrex::Real redP = alpha / k_f_save[10] * phase_units[10] * low_A[10] * exp(low_beta[10] * tc[0] - activation_units[10] * low_Ea[10] * invT);
        Corr[10] = redP / (1. + redP);
    }

    /* simple three-body correction */
    {
        amrex::Real alpha;
        alpha = mixture + (TB[11][0] - 1)*sc[1] + (TB[11][1] - 1)*sc[8] + (TB[11][2] - 1)*sc[11] + (TB[11][3] - 1)*sc[22];
        Corr[11] = alpha;
        alpha = mixture + (TB[12][0] - 1)*sc[1] + (TB[12][1] - 1)*sc[8] + (TB[12][2] - 1)*sc[11] + (TB[12][3] - 1)*sc[22];
        Corr[12] = alpha;
        alpha = mixture + (TB[13][0] - 1)*sc[1] + (TB[13][1] - 1)*sc[8] + (TB[13][2] - 1)*sc[11] + (TB[13][3] - 1)*sc[22];
        Corr[13] = alpha;
        alpha = mixture + (TB[14][0] - 1)*sc[1] + (TB[14][1] - 1)*sc[8] + (TB[14][2] - 1)*sc[11] + (TB[14][3] - 1)*sc[22];
        Corr[14] = alpha;
        alpha = mixture + (TB[15][0] - 1)*sc[1] + (TB[15][1] - 1)*sc[8] + (TB[15][2] - 1)*sc[11] + (TB[15][3] - 1)*sc[22];
        Corr[15] = alpha;
        alpha = mixture + (TB[16][0] - 1)*sc[1] + (TB[16][1] - 1)*sc[8] + (TB[16][2] - 1)*sc[11] + (TB[16][3] - 1)*sc[22];
        Corr[16] = alpha;
        alpha = mixture + (TB[17][0] - 1)*sc[1] + (TB[17][1] - 1)*sc[8] + (TB[17][2] - 1)*sc[11] + (TB[17][3] - 1)*sc[22];
        Corr[17] = alpha;
        alpha = mixture;
        Corr[18] = alpha;
        Corr[19] = alpha;
        alpha = mixture + (TB[20][0] - 1)*sc[8] + (TB[20][1] - 1)*sc[11] + (TB[20][2] - 1)*sc[22];
        Corr[20] = alpha;
        alpha = mixture;
        Corr[21] = alpha;
        Corr[22] = alpha;
    }

    for (int i=0; i<175; i++)
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

    amrex::Real q_f[175], q_r[175];
    amrex::Real sc_qss[1];
    comp_qfqr(q_f, q_r, sc, sc_qss, tc, invT);

    for (int i = 0; i < 175; ++i) {
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
    amrex::Real c[39]; /*temporary storage */
    amrex::Real PORT = 1e6 * (*P)/(8.31446261815324e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 39; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRateFR(q_f, q_r, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 175; ++id) {
        q_f[id] *= 1.0e-6;
        q_r[id] *= 1.0e-6;
    }
}


/*Returns the rate of progress for each reaction */
void CKQC(amrex::Real *  T, amrex::Real *  C, amrex::Real *  qdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 39; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    progressRate(qdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 39; ++id) {
        C[id] *= 1.0e-6;
    }

    for (id = 0; id < 175; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mass fractions */
void CKQYP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y, amrex::Real *  qdot)
{
    int id; /*loop counter */
    amrex::Real c[39]; /*temporary storage */
    amrex::Real YOW = 0; 
    amrex::Real PWORT; 
    amrex::Real imw[39];

    get_imw(imw);

    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H */
    YOW += y[1]*imw[1]; /*H2 */
    YOW += y[2]*imw[2]; /*CH2 */
    YOW += y[3]*imw[3]; /*CH2(S) */
    YOW += y[4]*imw[4]; /*CH3 */
    YOW += y[5]*imw[5]; /*O */
    YOW += y[6]*imw[6]; /*CH4 */
    YOW += y[7]*imw[7]; /*OH */
    YOW += y[8]*imw[8]; /*H2O */
    YOW += y[9]*imw[9]; /*C2H2 */
    YOW += y[10]*imw[10]; /*C2H3 */
    YOW += y[11]*imw[11]; /*CO */
    YOW += y[12]*imw[12]; /*C2H4 */
    YOW += y[13]*imw[13]; /*HCO */
    YOW += y[14]*imw[14]; /*C2H5 */
    YOW += y[15]*imw[15]; /*CH2O */
    YOW += y[16]*imw[16]; /*C2H6 */
    YOW += y[17]*imw[17]; /*CH2OH */
    YOW += y[18]*imw[18]; /*CH3O */
    YOW += y[19]*imw[19]; /*O2 */
    YOW += y[20]*imw[20]; /*HO2 */
    YOW += y[21]*imw[21]; /*H2O2 */
    YOW += y[22]*imw[22]; /*CO2 */
    YOW += y[23]*imw[23]; /*CH3HCO */
    YOW += y[24]*imw[24]; /*CH3OCH2 */
    YOW += y[25]*imw[25]; /*HCOOH */
    YOW += y[26]*imw[26]; /*CH3OCH3 */
    YOW += y[27]*imw[27]; /*HOCH2O */
    YOW += y[28]*imw[28]; /*CH3OCO */
    YOW += y[29]*imw[29]; /*CH3OCHO */
    YOW += y[30]*imw[30]; /*CH3OCH2O */
    YOW += y[31]*imw[31]; /*CH3OCH2OH */
    YOW += y[32]*imw[32]; /*OCH2OCHO */
    YOW += y[33]*imw[33]; /*HOCH2OCO */
    YOW += y[34]*imw[34]; /*CH3OCH2O2 */
    YOW += y[35]*imw[35]; /*CH2OCH2O2H */
    YOW += y[36]*imw[36]; /*HO2CH2OCHO */
    YOW += y[37]*imw[37]; /*O2CH2OCH2O2H */
    YOW += y[38]*imw[38]; /*N2 */
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
    c[21] = PWORT * y[21]*imw[21]; 
    c[22] = PWORT * y[22]*imw[22]; 
    c[23] = PWORT * y[23]*imw[23]; 
    c[24] = PWORT * y[24]*imw[24]; 
    c[25] = PWORT * y[25]*imw[25]; 
    c[26] = PWORT * y[26]*imw[26]; 
    c[27] = PWORT * y[27]*imw[27]; 
    c[28] = PWORT * y[28]*imw[28]; 
    c[29] = PWORT * y[29]*imw[29]; 
    c[30] = PWORT * y[30]*imw[30]; 
    c[31] = PWORT * y[31]*imw[31]; 
    c[32] = PWORT * y[32]*imw[32]; 
    c[33] = PWORT * y[33]*imw[33]; 
    c[34] = PWORT * y[34]*imw[34]; 
    c[35] = PWORT * y[35]*imw[35]; 
    c[36] = PWORT * y[36]*imw[36]; 
    c[37] = PWORT * y[37]*imw[37]; 
    c[38] = PWORT * y[38]*imw[38]; 

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 175; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKQXP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x, amrex::Real *  qdot)
{
    int id; /*loop counter */
    amrex::Real c[39]; /*temporary storage */
    amrex::Real PORT = 1e6 * (*P)/(8.31446261815324e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 39; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 175; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mass fractions */
void CKQYR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y, amrex::Real *  qdot)
{
    int id; /*loop counter */
    amrex::Real c[39]; /*temporary storage */
    amrex::Real imw[39];

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
    c[21] = 1e6 * (*rho) * y[21]*imw[21]; 
    c[22] = 1e6 * (*rho) * y[22]*imw[22]; 
    c[23] = 1e6 * (*rho) * y[23]*imw[23]; 
    c[24] = 1e6 * (*rho) * y[24]*imw[24]; 
    c[25] = 1e6 * (*rho) * y[25]*imw[25]; 
    c[26] = 1e6 * (*rho) * y[26]*imw[26]; 
    c[27] = 1e6 * (*rho) * y[27]*imw[27]; 
    c[28] = 1e6 * (*rho) * y[28]*imw[28]; 
    c[29] = 1e6 * (*rho) * y[29]*imw[29]; 
    c[30] = 1e6 * (*rho) * y[30]*imw[30]; 
    c[31] = 1e6 * (*rho) * y[31]*imw[31]; 
    c[32] = 1e6 * (*rho) * y[32]*imw[32]; 
    c[33] = 1e6 * (*rho) * y[33]*imw[33]; 
    c[34] = 1e6 * (*rho) * y[34]*imw[34]; 
    c[35] = 1e6 * (*rho) * y[35]*imw[35]; 
    c[36] = 1e6 * (*rho) * y[36]*imw[36]; 
    c[37] = 1e6 * (*rho) * y[37]*imw[37]; 
    c[38] = 1e6 * (*rho) * y[38]*imw[38]; 

    /*call progressRate */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 175; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mole fractions */
void CKQXR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  x, amrex::Real *  qdot)
{
    int id; /*loop counter */
    amrex::Real c[39]; /*temporary storage */
    amrex::Real XW = 0; /*See Eq 4, 11 in CK Manual */
    amrex::Real ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*2.015940; /*H2 */
    XW += x[2]*14.027090; /*CH2 */
    XW += x[3]*14.027090; /*CH2(S) */
    XW += x[4]*15.035060; /*CH3 */
    XW += x[5]*15.999400; /*O */
    XW += x[6]*16.043030; /*CH4 */
    XW += x[7]*17.007370; /*OH */
    XW += x[8]*18.015340; /*H2O */
    XW += x[9]*26.038240; /*C2H2 */
    XW += x[10]*27.046210; /*C2H3 */
    XW += x[11]*28.010550; /*CO */
    XW += x[12]*28.054180; /*C2H4 */
    XW += x[13]*29.018520; /*HCO */
    XW += x[14]*29.062150; /*C2H5 */
    XW += x[15]*30.026490; /*CH2O */
    XW += x[16]*30.070120; /*C2H6 */
    XW += x[17]*31.034460; /*CH2OH */
    XW += x[18]*31.034460; /*CH3O */
    XW += x[19]*31.998800; /*O2 */
    XW += x[20]*33.006770; /*HO2 */
    XW += x[21]*34.014740; /*H2O2 */
    XW += x[22]*44.009950; /*CO2 */
    XW += x[23]*44.053580; /*CH3HCO */
    XW += x[24]*45.061550; /*CH3OCH2 */
    XW += x[25]*46.025890; /*HCOOH */
    XW += x[26]*46.069520; /*CH3OCH3 */
    XW += x[27]*47.033860; /*HOCH2O */
    XW += x[28]*59.045010; /*CH3OCO */
    XW += x[29]*60.052980; /*CH3OCHO */
    XW += x[30]*61.060950; /*CH3OCH2O */
    XW += x[31]*62.068920; /*CH3OCH2OH */
    XW += x[32]*75.044410; /*OCH2OCHO */
    XW += x[33]*75.044410; /*HOCH2OCO */
    XW += x[34]*77.060350; /*CH3OCH2O2 */
    XW += x[35]*77.060350; /*CH2OCH2O2H */
    XW += x[36]*92.051780; /*HO2CH2OCHO */
    XW += x[37]*109.059150; /*O2CH2OCH2O2H */
    XW += x[38]*28.013400; /*N2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 39; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 175; ++id) {
        qdot[id] *= 1.0e-6;
    }
}

/*compute the reaction Jacobian on CPU */
void aJacobian(amrex::Real *  J, amrex::Real *  sc, amrex::Real T, int consP)
{
    for (int i=0; i<1600; i++) {
        J[i] = 0.0;
    }

    amrex::Real wdot[39];
    for (int k=0; k<39; k++) {
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
    for (int k = 0; k < 39; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    amrex::Real g_RT[39];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    amrex::Real h_RT[39];
    speciesEnthalpy(h_RT, tc);

    amrex::Real phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    amrex::Real dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    amrex::Real dqdci, dcdc_fac, dqdc[39];
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
    alpha = mixture + (TB[0][0] - 1)*sc[1] + (TB[0][1] - 1)*sc[8] + (TB[0][2] - 1)*sc[19] + (TB[0][3] - 1)*sc[11] + (TB[0][4] - 1)*sc[22];
    /* forward */
    phi_f = sc[0]*sc[19];
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
    phi_r = sc[20];
    Kc = refCinv * exp(g_RT[0] + g_RT[19] - g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[19]) + (h_RT[20]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[19] -= q; /* O2 */
    wdot[20] += q; /* HO2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[19];
        J[0] -= dqdci;                /* dwdot[H]/d[H] */
        J[19] -= dqdci;               /* dwdot[O2]/d[H] */
        J[20] += dqdci;               /* dwdot[HO2]/d[H] */
        /* d()/d[H2] */
        dqdci = (TB[0][0] - 1)*dcdc_fac;
        J[40] -= dqdci;               /* dwdot[H]/d[H2] */
        J[59] -= dqdci;               /* dwdot[O2]/d[H2] */
        J[60] += dqdci;               /* dwdot[HO2]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[0][1] - 1)*dcdc_fac;
        J[320] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[339] -= dqdci;              /* dwdot[O2]/d[H2O] */
        J[340] += dqdci;              /* dwdot[HO2]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[0][3] - 1)*dcdc_fac;
        J[440] -= dqdci;              /* dwdot[H]/d[CO] */
        J[459] -= dqdci;              /* dwdot[O2]/d[CO] */
        J[460] += dqdci;              /* dwdot[HO2]/d[CO] */
        /* d()/d[O2] */
        dqdci = (TB[0][2] - 1)*dcdc_fac + k_f*sc[0];
        J[760] -= dqdci;              /* dwdot[H]/d[O2] */
        J[779] -= dqdci;              /* dwdot[O2]/d[O2] */
        J[780] += dqdci;              /* dwdot[HO2]/d[O2] */
        /* d()/d[HO2] */
        dqdci =  - k_r;
        J[800] -= dqdci;              /* dwdot[H]/d[HO2] */
        J[819] -= dqdci;              /* dwdot[O2]/d[HO2] */
        J[820] += dqdci;              /* dwdot[HO2]/d[HO2] */
        /* d()/d[CO2] */
        dqdci = (TB[0][4] - 1)*dcdc_fac;
        J[880] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[899] -= dqdci;              /* dwdot[O2]/d[CO2] */
        J[900] += dqdci;              /* dwdot[HO2]/d[CO2] */
    }
    else {
        dqdc[0] = dcdc_fac + k_f*sc[19];
        dqdc[1] = TB[0][0]*dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = TB[0][1]*dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = TB[0][3]*dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = TB[0][2]*dcdc_fac + k_f*sc[0];
        dqdc[20] = dcdc_fac - k_r;
        dqdc[21] = dcdc_fac;
        dqdc[22] = TB[0][4]*dcdc_fac;
        dqdc[23] = dcdc_fac;
        dqdc[24] = dcdc_fac;
        dqdc[25] = dcdc_fac;
        dqdc[26] = dcdc_fac;
        dqdc[27] = dcdc_fac;
        dqdc[28] = dcdc_fac;
        dqdc[29] = dcdc_fac;
        dqdc[30] = dcdc_fac;
        dqdc[31] = dcdc_fac;
        dqdc[32] = dcdc_fac;
        dqdc[33] = dcdc_fac;
        dqdc[34] = dcdc_fac;
        dqdc[35] = dcdc_fac;
        dqdc[36] = dcdc_fac;
        dqdc[37] = dcdc_fac;
        dqdc[38] = dcdc_fac;
        for (int k=0; k<39; k++) {
            J[40*k+0] -= dqdc[k];
            J[40*k+19] -= dqdc[k];
            J[40*k+20] += dqdc[k];
        }
    }
    J[1560] -= dqdT; /* dwdot[H]/dT */
    J[1579] -= dqdT; /* dwdot[O2]/dT */
    J[1580] += dqdT; /* dwdot[HO2]/dT */

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[1][0] - 1)*sc[1] + (TB[1][1] - 1)*sc[8] + (TB[1][2] - 1)*sc[11] + (TB[1][3] - 1)*sc[22];
    /* forward */
    phi_f = sc[21];
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
    phi_r = pow(sc[7], 2.000000);
    Kc = refC * exp(-g_RT[7] - g_RT[7] + g_RT[21]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[21]) + (2.000000*h_RT[7]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[7] += 2 * q; /* OH */
    wdot[21] -= q; /* H2O2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[1][0] - 1)*dcdc_fac;
        J[47] += 2 * dqdci;           /* dwdot[OH]/d[H2] */
        J[61] -= dqdci;               /* dwdot[H2O2]/d[H2] */
        /* d()/d[OH] */
        dqdci =  - k_r*2.000000*sc[7];
        J[287] += 2 * dqdci;          /* dwdot[OH]/d[OH] */
        J[301] -= dqdci;              /* dwdot[H2O2]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (TB[1][1] - 1)*dcdc_fac;
        J[327] += 2 * dqdci;          /* dwdot[OH]/d[H2O] */
        J[341] -= dqdci;              /* dwdot[H2O2]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[1][2] - 1)*dcdc_fac;
        J[447] += 2 * dqdci;          /* dwdot[OH]/d[CO] */
        J[461] -= dqdci;              /* dwdot[H2O2]/d[CO] */
        /* d()/d[H2O2] */
        dqdci =  + k_f;
        J[847] += 2 * dqdci;          /* dwdot[OH]/d[H2O2] */
        J[861] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
        /* d()/d[CO2] */
        dqdci = (TB[1][3] - 1)*dcdc_fac;
        J[887] += 2 * dqdci;          /* dwdot[OH]/d[CO2] */
        J[901] -= dqdci;              /* dwdot[H2O2]/d[CO2] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = TB[1][0]*dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac - k_r*2.000000*sc[7];
        dqdc[8] = TB[1][1]*dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = TB[1][2]*dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac + k_f;
        dqdc[22] = TB[1][3]*dcdc_fac;
        dqdc[23] = dcdc_fac;
        dqdc[24] = dcdc_fac;
        dqdc[25] = dcdc_fac;
        dqdc[26] = dcdc_fac;
        dqdc[27] = dcdc_fac;
        dqdc[28] = dcdc_fac;
        dqdc[29] = dcdc_fac;
        dqdc[30] = dcdc_fac;
        dqdc[31] = dcdc_fac;
        dqdc[32] = dcdc_fac;
        dqdc[33] = dcdc_fac;
        dqdc[34] = dcdc_fac;
        dqdc[35] = dcdc_fac;
        dqdc[36] = dcdc_fac;
        dqdc[37] = dcdc_fac;
        dqdc[38] = dcdc_fac;
        for (int k=0; k<39; k++) {
            J[40*k+7] += 2 * dqdc[k];
            J[40*k+21] -= dqdc[k];
        }
    }
    J[1567] += 2 * dqdT; /* dwdot[OH]/dT */
    J[1581] -= dqdT; /* dwdot[H2O2]/dT */

    /*reaction 3: CH3 + CH3 (+M) <=> C2H6 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[2][0] - 1)*sc[8] + (TB[2][1] - 1)*sc[11] + (TB[2][2] - 1)*sc[22];
    /* forward */
    phi_f = pow(sc[4], 2.000000);
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
    phi_r = sc[16];
    Kc = refCinv * exp(g_RT[4] + g_RT[4] - g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[4]) + (h_RT[16]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[4] -= 2 * q; /* CH3 */
    wdot[16] += q; /* C2H6 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[CH3] */
        dqdci =  + k_f*2.000000*sc[4];
        J[164] += -2 * dqdci;         /* dwdot[CH3]/d[CH3] */
        J[176] += dqdci;              /* dwdot[C2H6]/d[CH3] */
        /* d()/d[H2O] */
        dqdci = (TB[2][0] - 1)*dcdc_fac;
        J[324] += -2 * dqdci;         /* dwdot[CH3]/d[H2O] */
        J[336] += dqdci;              /* dwdot[C2H6]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[2][1] - 1)*dcdc_fac;
        J[444] += -2 * dqdci;         /* dwdot[CH3]/d[CO] */
        J[456] += dqdci;              /* dwdot[C2H6]/d[CO] */
        /* d()/d[C2H6] */
        dqdci =  - k_r;
        J[644] += -2 * dqdci;         /* dwdot[CH3]/d[C2H6] */
        J[656] += dqdci;              /* dwdot[C2H6]/d[C2H6] */
        /* d()/d[CO2] */
        dqdci = (TB[2][2] - 1)*dcdc_fac;
        J[884] += -2 * dqdci;         /* dwdot[CH3]/d[CO2] */
        J[896] += dqdci;              /* dwdot[C2H6]/d[CO2] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac + k_f*2.000000*sc[4];
        dqdc[5] = dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = TB[2][0]*dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = TB[2][1]*dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac - k_r;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac;
        dqdc[22] = TB[2][2]*dcdc_fac;
        dqdc[23] = dcdc_fac;
        dqdc[24] = dcdc_fac;
        dqdc[25] = dcdc_fac;
        dqdc[26] = dcdc_fac;
        dqdc[27] = dcdc_fac;
        dqdc[28] = dcdc_fac;
        dqdc[29] = dcdc_fac;
        dqdc[30] = dcdc_fac;
        dqdc[31] = dcdc_fac;
        dqdc[32] = dcdc_fac;
        dqdc[33] = dcdc_fac;
        dqdc[34] = dcdc_fac;
        dqdc[35] = dcdc_fac;
        dqdc[36] = dcdc_fac;
        dqdc[37] = dcdc_fac;
        dqdc[38] = dcdc_fac;
        for (int k=0; k<39; k++) {
            J[40*k+4] += -2 * dqdc[k];
            J[40*k+16] += dqdc[k];
        }
    }
    J[1564] += -2 * dqdT; /* dwdot[CH3]/dT */
    J[1576] += dqdT; /* dwdot[C2H6]/dT */

    /*reaction 4: CH3 + H (+M) <=> CH4 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[3][0] - 1)*sc[1] + (TB[3][1] - 1)*sc[8] + (TB[3][2] - 1)*sc[6] + (TB[3][3] - 1)*sc[11] + (TB[3][4] - 1)*sc[22] + (TB[3][5] - 1)*sc[16];
    /* forward */
    phi_f = sc[0]*sc[4];
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
    phi_r = sc[6];
    Kc = refCinv * exp(g_RT[0] + g_RT[4] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[4]) + (h_RT[6]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[4] -= q; /* CH3 */
    wdot[6] += q; /* CH4 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[4];
        J[0] -= dqdci;                /* dwdot[H]/d[H] */
        J[4] -= dqdci;                /* dwdot[CH3]/d[H] */
        J[6] += dqdci;                /* dwdot[CH4]/d[H] */
        /* d()/d[H2] */
        dqdci = (TB[3][0] - 1)*dcdc_fac;
        J[40] -= dqdci;               /* dwdot[H]/d[H2] */
        J[44] -= dqdci;               /* dwdot[CH3]/d[H2] */
        J[46] += dqdci;               /* dwdot[CH4]/d[H2] */
        /* d()/d[CH3] */
        dqdci =  + k_f*sc[0];
        J[160] -= dqdci;              /* dwdot[H]/d[CH3] */
        J[164] -= dqdci;              /* dwdot[CH3]/d[CH3] */
        J[166] += dqdci;              /* dwdot[CH4]/d[CH3] */
        /* d()/d[CH4] */
        dqdci = (TB[3][2] - 1)*dcdc_fac - k_r;
        J[240] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[244] -= dqdci;              /* dwdot[CH3]/d[CH4] */
        J[246] += dqdci;              /* dwdot[CH4]/d[CH4] */
        /* d()/d[H2O] */
        dqdci = (TB[3][1] - 1)*dcdc_fac;
        J[320] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[324] -= dqdci;              /* dwdot[CH3]/d[H2O] */
        J[326] += dqdci;              /* dwdot[CH4]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[3][3] - 1)*dcdc_fac;
        J[440] -= dqdci;              /* dwdot[H]/d[CO] */
        J[444] -= dqdci;              /* dwdot[CH3]/d[CO] */
        J[446] += dqdci;              /* dwdot[CH4]/d[CO] */
        /* d()/d[C2H6] */
        dqdci = (TB[3][5] - 1)*dcdc_fac;
        J[640] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[644] -= dqdci;              /* dwdot[CH3]/d[C2H6] */
        J[646] += dqdci;              /* dwdot[CH4]/d[C2H6] */
        /* d()/d[CO2] */
        dqdci = (TB[3][4] - 1)*dcdc_fac;
        J[880] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[884] -= dqdci;              /* dwdot[CH3]/d[CO2] */
        J[886] += dqdci;              /* dwdot[CH4]/d[CO2] */
    }
    else {
        dqdc[0] = dcdc_fac + k_f*sc[4];
        dqdc[1] = TB[3][0]*dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac + k_f*sc[0];
        dqdc[5] = dcdc_fac;
        dqdc[6] = TB[3][2]*dcdc_fac - k_r;
        dqdc[7] = dcdc_fac;
        dqdc[8] = TB[3][1]*dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = TB[3][3]*dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = TB[3][5]*dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac;
        dqdc[22] = TB[3][4]*dcdc_fac;
        dqdc[23] = dcdc_fac;
        dqdc[24] = dcdc_fac;
        dqdc[25] = dcdc_fac;
        dqdc[26] = dcdc_fac;
        dqdc[27] = dcdc_fac;
        dqdc[28] = dcdc_fac;
        dqdc[29] = dcdc_fac;
        dqdc[30] = dcdc_fac;
        dqdc[31] = dcdc_fac;
        dqdc[32] = dcdc_fac;
        dqdc[33] = dcdc_fac;
        dqdc[34] = dcdc_fac;
        dqdc[35] = dcdc_fac;
        dqdc[36] = dcdc_fac;
        dqdc[37] = dcdc_fac;
        dqdc[38] = dcdc_fac;
        for (int k=0; k<39; k++) {
            J[40*k+0] -= dqdc[k];
            J[40*k+4] -= dqdc[k];
            J[40*k+6] += dqdc[k];
        }
    }
    J[1560] -= dqdT; /* dwdot[H]/dT */
    J[1564] -= dqdT; /* dwdot[CH3]/dT */
    J[1566] += dqdT; /* dwdot[CH4]/dT */

    /*reaction 5: C2H5 + H (+M) <=> C2H6 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[4][0] - 1)*sc[1] + (TB[4][1] - 1)*sc[8] + (TB[4][2] - 1)*sc[6] + (TB[4][3] - 1)*sc[11] + (TB[4][4] - 1)*sc[22] + (TB[4][5] - 1)*sc[16];
    /* forward */
    phi_f = sc[0]*sc[14];
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
    Kc = refCinv * exp(g_RT[0] + g_RT[14] - g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[14]) + (h_RT[16]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[14] -= q; /* C2H5 */
    wdot[16] += q; /* C2H6 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[14];
        J[0] -= dqdci;                /* dwdot[H]/d[H] */
        J[14] -= dqdci;               /* dwdot[C2H5]/d[H] */
        J[16] += dqdci;               /* dwdot[C2H6]/d[H] */
        /* d()/d[H2] */
        dqdci = (TB[4][0] - 1)*dcdc_fac;
        J[40] -= dqdci;               /* dwdot[H]/d[H2] */
        J[54] -= dqdci;               /* dwdot[C2H5]/d[H2] */
        J[56] += dqdci;               /* dwdot[C2H6]/d[H2] */
        /* d()/d[CH4] */
        dqdci = (TB[4][2] - 1)*dcdc_fac;
        J[240] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[254] -= dqdci;              /* dwdot[C2H5]/d[CH4] */
        J[256] += dqdci;              /* dwdot[C2H6]/d[CH4] */
        /* d()/d[H2O] */
        dqdci = (TB[4][1] - 1)*dcdc_fac;
        J[320] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[334] -= dqdci;              /* dwdot[C2H5]/d[H2O] */
        J[336] += dqdci;              /* dwdot[C2H6]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[4][3] - 1)*dcdc_fac;
        J[440] -= dqdci;              /* dwdot[H]/d[CO] */
        J[454] -= dqdci;              /* dwdot[C2H5]/d[CO] */
        J[456] += dqdci;              /* dwdot[C2H6]/d[CO] */
        /* d()/d[C2H5] */
        dqdci =  + k_f*sc[0];
        J[560] -= dqdci;              /* dwdot[H]/d[C2H5] */
        J[574] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
        J[576] += dqdci;              /* dwdot[C2H6]/d[C2H5] */
        /* d()/d[C2H6] */
        dqdci = (TB[4][5] - 1)*dcdc_fac - k_r;
        J[640] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[654] -= dqdci;              /* dwdot[C2H5]/d[C2H6] */
        J[656] += dqdci;              /* dwdot[C2H6]/d[C2H6] */
        /* d()/d[CO2] */
        dqdci = (TB[4][4] - 1)*dcdc_fac;
        J[880] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[894] -= dqdci;              /* dwdot[C2H5]/d[CO2] */
        J[896] += dqdci;              /* dwdot[C2H6]/d[CO2] */
    }
    else {
        dqdc[0] = dcdc_fac + k_f*sc[14];
        dqdc[1] = TB[4][0]*dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = TB[4][2]*dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = TB[4][1]*dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = TB[4][3]*dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac + k_f*sc[0];
        dqdc[15] = dcdc_fac;
        dqdc[16] = TB[4][5]*dcdc_fac - k_r;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac;
        dqdc[22] = TB[4][4]*dcdc_fac;
        dqdc[23] = dcdc_fac;
        dqdc[24] = dcdc_fac;
        dqdc[25] = dcdc_fac;
        dqdc[26] = dcdc_fac;
        dqdc[27] = dcdc_fac;
        dqdc[28] = dcdc_fac;
        dqdc[29] = dcdc_fac;
        dqdc[30] = dcdc_fac;
        dqdc[31] = dcdc_fac;
        dqdc[32] = dcdc_fac;
        dqdc[33] = dcdc_fac;
        dqdc[34] = dcdc_fac;
        dqdc[35] = dcdc_fac;
        dqdc[36] = dcdc_fac;
        dqdc[37] = dcdc_fac;
        dqdc[38] = dcdc_fac;
        for (int k=0; k<39; k++) {
            J[40*k+0] -= dqdc[k];
            J[40*k+14] -= dqdc[k];
            J[40*k+16] += dqdc[k];
        }
    }
    J[1560] -= dqdT; /* dwdot[H]/dT */
    J[1574] -= dqdT; /* dwdot[C2H5]/dT */
    J[1576] += dqdT; /* dwdot[C2H6]/dT */

    /*reaction 6: C2H4 (+M) <=> H2 + C2H2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[5][0] - 1)*sc[1] + (TB[5][1] - 1)*sc[8] + (TB[5][2] - 1)*sc[6] + (TB[5][3] - 1)*sc[11] + (TB[5][4] - 1)*sc[22] + (TB[5][5] - 1)*sc[16];
    /* forward */
    phi_f = sc[12];
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
    phi_r = sc[1]*sc[9];
    Kc = refC * exp(-g_RT[1] - g_RT[9] + g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[12]) + (h_RT[1] + h_RT[9]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] += q; /* H2 */
    wdot[9] += q; /* C2H2 */
    wdot[12] -= q; /* C2H4 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[5][0] - 1)*dcdc_fac - k_r*sc[9];
        J[41] += dqdci;               /* dwdot[H2]/d[H2] */
        J[49] += dqdci;               /* dwdot[C2H2]/d[H2] */
        J[52] -= dqdci;               /* dwdot[C2H4]/d[H2] */
        /* d()/d[CH4] */
        dqdci = (TB[5][2] - 1)*dcdc_fac;
        J[241] += dqdci;              /* dwdot[H2]/d[CH4] */
        J[249] += dqdci;              /* dwdot[C2H2]/d[CH4] */
        J[252] -= dqdci;              /* dwdot[C2H4]/d[CH4] */
        /* d()/d[H2O] */
        dqdci = (TB[5][1] - 1)*dcdc_fac;
        J[321] += dqdci;              /* dwdot[H2]/d[H2O] */
        J[329] += dqdci;              /* dwdot[C2H2]/d[H2O] */
        J[332] -= dqdci;              /* dwdot[C2H4]/d[H2O] */
        /* d()/d[C2H2] */
        dqdci =  - k_r*sc[1];
        J[361] += dqdci;              /* dwdot[H2]/d[C2H2] */
        J[369] += dqdci;              /* dwdot[C2H2]/d[C2H2] */
        J[372] -= dqdci;              /* dwdot[C2H4]/d[C2H2] */
        /* d()/d[CO] */
        dqdci = (TB[5][3] - 1)*dcdc_fac;
        J[441] += dqdci;              /* dwdot[H2]/d[CO] */
        J[449] += dqdci;              /* dwdot[C2H2]/d[CO] */
        J[452] -= dqdci;              /* dwdot[C2H4]/d[CO] */
        /* d()/d[C2H4] */
        dqdci =  + k_f;
        J[481] += dqdci;              /* dwdot[H2]/d[C2H4] */
        J[489] += dqdci;              /* dwdot[C2H2]/d[C2H4] */
        J[492] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
        /* d()/d[C2H6] */
        dqdci = (TB[5][5] - 1)*dcdc_fac;
        J[641] += dqdci;              /* dwdot[H2]/d[C2H6] */
        J[649] += dqdci;              /* dwdot[C2H2]/d[C2H6] */
        J[652] -= dqdci;              /* dwdot[C2H4]/d[C2H6] */
        /* d()/d[CO2] */
        dqdci = (TB[5][4] - 1)*dcdc_fac;
        J[881] += dqdci;              /* dwdot[H2]/d[CO2] */
        J[889] += dqdci;              /* dwdot[C2H2]/d[CO2] */
        J[892] -= dqdci;              /* dwdot[C2H4]/d[CO2] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = TB[5][0]*dcdc_fac - k_r*sc[9];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = TB[5][2]*dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = TB[5][1]*dcdc_fac;
        dqdc[9] = dcdc_fac - k_r*sc[1];
        dqdc[10] = dcdc_fac;
        dqdc[11] = TB[5][3]*dcdc_fac;
        dqdc[12] = dcdc_fac + k_f;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = TB[5][5]*dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac;
        dqdc[22] = TB[5][4]*dcdc_fac;
        dqdc[23] = dcdc_fac;
        dqdc[24] = dcdc_fac;
        dqdc[25] = dcdc_fac;
        dqdc[26] = dcdc_fac;
        dqdc[27] = dcdc_fac;
        dqdc[28] = dcdc_fac;
        dqdc[29] = dcdc_fac;
        dqdc[30] = dcdc_fac;
        dqdc[31] = dcdc_fac;
        dqdc[32] = dcdc_fac;
        dqdc[33] = dcdc_fac;
        dqdc[34] = dcdc_fac;
        dqdc[35] = dcdc_fac;
        dqdc[36] = dcdc_fac;
        dqdc[37] = dcdc_fac;
        dqdc[38] = dcdc_fac;
        for (int k=0; k<39; k++) {
            J[40*k+1] += dqdc[k];
            J[40*k+9] += dqdc[k];
            J[40*k+12] -= dqdc[k];
        }
    }
    J[1561] += dqdT; /* dwdot[H2]/dT */
    J[1569] += dqdT; /* dwdot[C2H2]/dT */
    J[1572] -= dqdT; /* dwdot[C2H4]/dT */

    /*reaction 7: C2H4 + H (+M) <=> C2H5 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[6][0] - 1)*sc[1] + (TB[6][1] - 1)*sc[8] + (TB[6][2] - 1)*sc[6] + (TB[6][3] - 1)*sc[11] + (TB[6][4] - 1)*sc[22] + (TB[6][5] - 1)*sc[16];
    /* forward */
    phi_f = sc[0]*sc[12];
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
    Kc = refCinv * exp(g_RT[0] + g_RT[12] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[12]) + (h_RT[14]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[12] -= q; /* C2H4 */
    wdot[14] += q; /* C2H5 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[12];
        J[0] -= dqdci;                /* dwdot[H]/d[H] */
        J[12] -= dqdci;               /* dwdot[C2H4]/d[H] */
        J[14] += dqdci;               /* dwdot[C2H5]/d[H] */
        /* d()/d[H2] */
        dqdci = (TB[6][0] - 1)*dcdc_fac;
        J[40] -= dqdci;               /* dwdot[H]/d[H2] */
        J[52] -= dqdci;               /* dwdot[C2H4]/d[H2] */
        J[54] += dqdci;               /* dwdot[C2H5]/d[H2] */
        /* d()/d[CH4] */
        dqdci = (TB[6][2] - 1)*dcdc_fac;
        J[240] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[252] -= dqdci;              /* dwdot[C2H4]/d[CH4] */
        J[254] += dqdci;              /* dwdot[C2H5]/d[CH4] */
        /* d()/d[H2O] */
        dqdci = (TB[6][1] - 1)*dcdc_fac;
        J[320] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[332] -= dqdci;              /* dwdot[C2H4]/d[H2O] */
        J[334] += dqdci;              /* dwdot[C2H5]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[6][3] - 1)*dcdc_fac;
        J[440] -= dqdci;              /* dwdot[H]/d[CO] */
        J[452] -= dqdci;              /* dwdot[C2H4]/d[CO] */
        J[454] += dqdci;              /* dwdot[C2H5]/d[CO] */
        /* d()/d[C2H4] */
        dqdci =  + k_f*sc[0];
        J[480] -= dqdci;              /* dwdot[H]/d[C2H4] */
        J[492] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
        J[494] += dqdci;              /* dwdot[C2H5]/d[C2H4] */
        /* d()/d[C2H5] */
        dqdci =  - k_r;
        J[560] -= dqdci;              /* dwdot[H]/d[C2H5] */
        J[572] -= dqdci;              /* dwdot[C2H4]/d[C2H5] */
        J[574] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
        /* d()/d[C2H6] */
        dqdci = (TB[6][5] - 1)*dcdc_fac;
        J[640] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[652] -= dqdci;              /* dwdot[C2H4]/d[C2H6] */
        J[654] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
        /* d()/d[CO2] */
        dqdci = (TB[6][4] - 1)*dcdc_fac;
        J[880] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[892] -= dqdci;              /* dwdot[C2H4]/d[CO2] */
        J[894] += dqdci;              /* dwdot[C2H5]/d[CO2] */
    }
    else {
        dqdc[0] = dcdc_fac + k_f*sc[12];
        dqdc[1] = TB[6][0]*dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = TB[6][2]*dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = TB[6][1]*dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = TB[6][3]*dcdc_fac;
        dqdc[12] = dcdc_fac + k_f*sc[0];
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac - k_r;
        dqdc[15] = dcdc_fac;
        dqdc[16] = TB[6][5]*dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac;
        dqdc[22] = TB[6][4]*dcdc_fac;
        dqdc[23] = dcdc_fac;
        dqdc[24] = dcdc_fac;
        dqdc[25] = dcdc_fac;
        dqdc[26] = dcdc_fac;
        dqdc[27] = dcdc_fac;
        dqdc[28] = dcdc_fac;
        dqdc[29] = dcdc_fac;
        dqdc[30] = dcdc_fac;
        dqdc[31] = dcdc_fac;
        dqdc[32] = dcdc_fac;
        dqdc[33] = dcdc_fac;
        dqdc[34] = dcdc_fac;
        dqdc[35] = dcdc_fac;
        dqdc[36] = dcdc_fac;
        dqdc[37] = dcdc_fac;
        dqdc[38] = dcdc_fac;
        for (int k=0; k<39; k++) {
            J[40*k+0] -= dqdc[k];
            J[40*k+12] -= dqdc[k];
            J[40*k+14] += dqdc[k];
        }
    }
    J[1560] -= dqdT; /* dwdot[H]/dT */
    J[1572] -= dqdT; /* dwdot[C2H4]/dT */
    J[1574] += dqdT; /* dwdot[C2H5]/dT */

    /*reaction 8: C2H3 + H (+M) <=> C2H4 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[7][0] - 1)*sc[1] + (TB[7][1] - 1)*sc[8] + (TB[7][2] - 1)*sc[6] + (TB[7][3] - 1)*sc[11] + (TB[7][4] - 1)*sc[22] + (TB[7][5] - 1)*sc[16];
    /* forward */
    phi_f = sc[0]*sc[10];
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
    phi_r = sc[12];
    Kc = refCinv * exp(g_RT[0] + g_RT[10] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[10]) + (h_RT[12]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[10] -= q; /* C2H3 */
    wdot[12] += q; /* C2H4 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[10];
        J[0] -= dqdci;                /* dwdot[H]/d[H] */
        J[10] -= dqdci;               /* dwdot[C2H3]/d[H] */
        J[12] += dqdci;               /* dwdot[C2H4]/d[H] */
        /* d()/d[H2] */
        dqdci = (TB[7][0] - 1)*dcdc_fac;
        J[40] -= dqdci;               /* dwdot[H]/d[H2] */
        J[50] -= dqdci;               /* dwdot[C2H3]/d[H2] */
        J[52] += dqdci;               /* dwdot[C2H4]/d[H2] */
        /* d()/d[CH4] */
        dqdci = (TB[7][2] - 1)*dcdc_fac;
        J[240] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[250] -= dqdci;              /* dwdot[C2H3]/d[CH4] */
        J[252] += dqdci;              /* dwdot[C2H4]/d[CH4] */
        /* d()/d[H2O] */
        dqdci = (TB[7][1] - 1)*dcdc_fac;
        J[320] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[330] -= dqdci;              /* dwdot[C2H3]/d[H2O] */
        J[332] += dqdci;              /* dwdot[C2H4]/d[H2O] */
        /* d()/d[C2H3] */
        dqdci =  + k_f*sc[0];
        J[400] -= dqdci;              /* dwdot[H]/d[C2H3] */
        J[410] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
        J[412] += dqdci;              /* dwdot[C2H4]/d[C2H3] */
        /* d()/d[CO] */
        dqdci = (TB[7][3] - 1)*dcdc_fac;
        J[440] -= dqdci;              /* dwdot[H]/d[CO] */
        J[450] -= dqdci;              /* dwdot[C2H3]/d[CO] */
        J[452] += dqdci;              /* dwdot[C2H4]/d[CO] */
        /* d()/d[C2H4] */
        dqdci =  - k_r;
        J[480] -= dqdci;              /* dwdot[H]/d[C2H4] */
        J[490] -= dqdci;              /* dwdot[C2H3]/d[C2H4] */
        J[492] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
        /* d()/d[C2H6] */
        dqdci = (TB[7][5] - 1)*dcdc_fac;
        J[640] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[650] -= dqdci;              /* dwdot[C2H3]/d[C2H6] */
        J[652] += dqdci;              /* dwdot[C2H4]/d[C2H6] */
        /* d()/d[CO2] */
        dqdci = (TB[7][4] - 1)*dcdc_fac;
        J[880] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[890] -= dqdci;              /* dwdot[C2H3]/d[CO2] */
        J[892] += dqdci;              /* dwdot[C2H4]/d[CO2] */
    }
    else {
        dqdc[0] = dcdc_fac + k_f*sc[10];
        dqdc[1] = TB[7][0]*dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = TB[7][2]*dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = TB[7][1]*dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac + k_f*sc[0];
        dqdc[11] = TB[7][3]*dcdc_fac;
        dqdc[12] = dcdc_fac - k_r;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = TB[7][5]*dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac;
        dqdc[22] = TB[7][4]*dcdc_fac;
        dqdc[23] = dcdc_fac;
        dqdc[24] = dcdc_fac;
        dqdc[25] = dcdc_fac;
        dqdc[26] = dcdc_fac;
        dqdc[27] = dcdc_fac;
        dqdc[28] = dcdc_fac;
        dqdc[29] = dcdc_fac;
        dqdc[30] = dcdc_fac;
        dqdc[31] = dcdc_fac;
        dqdc[32] = dcdc_fac;
        dqdc[33] = dcdc_fac;
        dqdc[34] = dcdc_fac;
        dqdc[35] = dcdc_fac;
        dqdc[36] = dcdc_fac;
        dqdc[37] = dcdc_fac;
        dqdc[38] = dcdc_fac;
        for (int k=0; k<39; k++) {
            J[40*k+0] -= dqdc[k];
            J[40*k+10] -= dqdc[k];
            J[40*k+12] += dqdc[k];
        }
    }
    J[1560] -= dqdT; /* dwdot[H]/dT */
    J[1570] -= dqdT; /* dwdot[C2H3]/dT */
    J[1572] += dqdT; /* dwdot[C2H4]/dT */

    /*reaction 9: C2H2 + H (+M) <=> C2H3 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[8][0] - 1)*sc[1] + (TB[8][1] - 1)*sc[8] + (TB[8][2] - 1)*sc[6] + (TB[8][3] - 1)*sc[11] + (TB[8][4] - 1)*sc[22] + (TB[8][5] - 1)*sc[16];
    /* forward */
    phi_f = sc[0]*sc[9];
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
    phi_r = sc[10];
    Kc = refCinv * exp(g_RT[0] + g_RT[9] - g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[9]) + (h_RT[10]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[9] -= q; /* C2H2 */
    wdot[10] += q; /* C2H3 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[9];
        J[0] -= dqdci;                /* dwdot[H]/d[H] */
        J[9] -= dqdci;                /* dwdot[C2H2]/d[H] */
        J[10] += dqdci;               /* dwdot[C2H3]/d[H] */
        /* d()/d[H2] */
        dqdci = (TB[8][0] - 1)*dcdc_fac;
        J[40] -= dqdci;               /* dwdot[H]/d[H2] */
        J[49] -= dqdci;               /* dwdot[C2H2]/d[H2] */
        J[50] += dqdci;               /* dwdot[C2H3]/d[H2] */
        /* d()/d[CH4] */
        dqdci = (TB[8][2] - 1)*dcdc_fac;
        J[240] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[249] -= dqdci;              /* dwdot[C2H2]/d[CH4] */
        J[250] += dqdci;              /* dwdot[C2H3]/d[CH4] */
        /* d()/d[H2O] */
        dqdci = (TB[8][1] - 1)*dcdc_fac;
        J[320] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[329] -= dqdci;              /* dwdot[C2H2]/d[H2O] */
        J[330] += dqdci;              /* dwdot[C2H3]/d[H2O] */
        /* d()/d[C2H2] */
        dqdci =  + k_f*sc[0];
        J[360] -= dqdci;              /* dwdot[H]/d[C2H2] */
        J[369] -= dqdci;              /* dwdot[C2H2]/d[C2H2] */
        J[370] += dqdci;              /* dwdot[C2H3]/d[C2H2] */
        /* d()/d[C2H3] */
        dqdci =  - k_r;
        J[400] -= dqdci;              /* dwdot[H]/d[C2H3] */
        J[409] -= dqdci;              /* dwdot[C2H2]/d[C2H3] */
        J[410] += dqdci;              /* dwdot[C2H3]/d[C2H3] */
        /* d()/d[CO] */
        dqdci = (TB[8][3] - 1)*dcdc_fac;
        J[440] -= dqdci;              /* dwdot[H]/d[CO] */
        J[449] -= dqdci;              /* dwdot[C2H2]/d[CO] */
        J[450] += dqdci;              /* dwdot[C2H3]/d[CO] */
        /* d()/d[C2H6] */
        dqdci = (TB[8][5] - 1)*dcdc_fac;
        J[640] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[649] -= dqdci;              /* dwdot[C2H2]/d[C2H6] */
        J[650] += dqdci;              /* dwdot[C2H3]/d[C2H6] */
        /* d()/d[CO2] */
        dqdci = (TB[8][4] - 1)*dcdc_fac;
        J[880] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[889] -= dqdci;              /* dwdot[C2H2]/d[CO2] */
        J[890] += dqdci;              /* dwdot[C2H3]/d[CO2] */
    }
    else {
        dqdc[0] = dcdc_fac + k_f*sc[9];
        dqdc[1] = TB[8][0]*dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = TB[8][2]*dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = TB[8][1]*dcdc_fac;
        dqdc[9] = dcdc_fac + k_f*sc[0];
        dqdc[10] = dcdc_fac - k_r;
        dqdc[11] = TB[8][3]*dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = TB[8][5]*dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac;
        dqdc[22] = TB[8][4]*dcdc_fac;
        dqdc[23] = dcdc_fac;
        dqdc[24] = dcdc_fac;
        dqdc[25] = dcdc_fac;
        dqdc[26] = dcdc_fac;
        dqdc[27] = dcdc_fac;
        dqdc[28] = dcdc_fac;
        dqdc[29] = dcdc_fac;
        dqdc[30] = dcdc_fac;
        dqdc[31] = dcdc_fac;
        dqdc[32] = dcdc_fac;
        dqdc[33] = dcdc_fac;
        dqdc[34] = dcdc_fac;
        dqdc[35] = dcdc_fac;
        dqdc[36] = dcdc_fac;
        dqdc[37] = dcdc_fac;
        dqdc[38] = dcdc_fac;
        for (int k=0; k<39; k++) {
            J[40*k+0] -= dqdc[k];
            J[40*k+9] -= dqdc[k];
            J[40*k+10] += dqdc[k];
        }
    }
    J[1560] -= dqdT; /* dwdot[H]/dT */
    J[1569] -= dqdT; /* dwdot[C2H2]/dT */
    J[1570] += dqdT; /* dwdot[C2H3]/dT */

    /*reaction 10: CH2 + H (+M) <=> CH3 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[9][0] - 1)*sc[1] + (TB[9][1] - 1)*sc[8] + (TB[9][2] - 1)*sc[6] + (TB[9][3] - 1)*sc[11] + (TB[9][4] - 1)*sc[22] + (TB[9][5] - 1)*sc[16];
    /* forward */
    phi_f = sc[0]*sc[2];
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
    phi_r = sc[4];
    Kc = refCinv * exp(g_RT[0] + g_RT[2] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[2]) + (h_RT[4]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[2] -= q; /* CH2 */
    wdot[4] += q; /* CH3 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[2];
        J[0] -= dqdci;                /* dwdot[H]/d[H] */
        J[2] -= dqdci;                /* dwdot[CH2]/d[H] */
        J[4] += dqdci;                /* dwdot[CH3]/d[H] */
        /* d()/d[H2] */
        dqdci = (TB[9][0] - 1)*dcdc_fac;
        J[40] -= dqdci;               /* dwdot[H]/d[H2] */
        J[42] -= dqdci;               /* dwdot[CH2]/d[H2] */
        J[44] += dqdci;               /* dwdot[CH3]/d[H2] */
        /* d()/d[CH2] */
        dqdci =  + k_f*sc[0];
        J[80] -= dqdci;               /* dwdot[H]/d[CH2] */
        J[82] -= dqdci;               /* dwdot[CH2]/d[CH2] */
        J[84] += dqdci;               /* dwdot[CH3]/d[CH2] */
        /* d()/d[CH3] */
        dqdci =  - k_r;
        J[160] -= dqdci;              /* dwdot[H]/d[CH3] */
        J[162] -= dqdci;              /* dwdot[CH2]/d[CH3] */
        J[164] += dqdci;              /* dwdot[CH3]/d[CH3] */
        /* d()/d[CH4] */
        dqdci = (TB[9][2] - 1)*dcdc_fac;
        J[240] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[242] -= dqdci;              /* dwdot[CH2]/d[CH4] */
        J[244] += dqdci;              /* dwdot[CH3]/d[CH4] */
        /* d()/d[H2O] */
        dqdci = (TB[9][1] - 1)*dcdc_fac;
        J[320] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[322] -= dqdci;              /* dwdot[CH2]/d[H2O] */
        J[324] += dqdci;              /* dwdot[CH3]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[9][3] - 1)*dcdc_fac;
        J[440] -= dqdci;              /* dwdot[H]/d[CO] */
        J[442] -= dqdci;              /* dwdot[CH2]/d[CO] */
        J[444] += dqdci;              /* dwdot[CH3]/d[CO] */
        /* d()/d[C2H6] */
        dqdci = (TB[9][5] - 1)*dcdc_fac;
        J[640] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[642] -= dqdci;              /* dwdot[CH2]/d[C2H6] */
        J[644] += dqdci;              /* dwdot[CH3]/d[C2H6] */
        /* d()/d[CO2] */
        dqdci = (TB[9][4] - 1)*dcdc_fac;
        J[880] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[882] -= dqdci;              /* dwdot[CH2]/d[CO2] */
        J[884] += dqdci;              /* dwdot[CH3]/d[CO2] */
    }
    else {
        dqdc[0] = dcdc_fac + k_f*sc[2];
        dqdc[1] = TB[9][0]*dcdc_fac;
        dqdc[2] = dcdc_fac + k_f*sc[0];
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac - k_r;
        dqdc[5] = dcdc_fac;
        dqdc[6] = TB[9][2]*dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = TB[9][1]*dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = TB[9][3]*dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = TB[9][5]*dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac;
        dqdc[22] = TB[9][4]*dcdc_fac;
        dqdc[23] = dcdc_fac;
        dqdc[24] = dcdc_fac;
        dqdc[25] = dcdc_fac;
        dqdc[26] = dcdc_fac;
        dqdc[27] = dcdc_fac;
        dqdc[28] = dcdc_fac;
        dqdc[29] = dcdc_fac;
        dqdc[30] = dcdc_fac;
        dqdc[31] = dcdc_fac;
        dqdc[32] = dcdc_fac;
        dqdc[33] = dcdc_fac;
        dqdc[34] = dcdc_fac;
        dqdc[35] = dcdc_fac;
        dqdc[36] = dcdc_fac;
        dqdc[37] = dcdc_fac;
        dqdc[38] = dcdc_fac;
        for (int k=0; k<39; k++) {
            J[40*k+0] -= dqdc[k];
            J[40*k+2] -= dqdc[k];
            J[40*k+4] += dqdc[k];
        }
    }
    J[1560] -= dqdT; /* dwdot[H]/dT */
    J[1562] -= dqdT; /* dwdot[CH2]/dT */
    J[1564] += dqdT; /* dwdot[CH3]/dT */

    /*reaction 11: CO + O (+M) <=> CO2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[10][0] - 1)*sc[1] + (TB[10][1] - 1)*sc[8] + (TB[10][2] - 1)*sc[11] + (TB[10][3] - 1)*sc[22];
    /* forward */
    phi_f = sc[5]*sc[11];
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
    /* Lindemann form */
    F = 1.0;
    dlogFdlogPr = 0.0;
    dlogFdT = 0.0;
    /* reverse */
    phi_r = sc[22];
    Kc = refCinv * exp(g_RT[5] + g_RT[11] - g_RT[22]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[11]) + (h_RT[22]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[5] -= q; /* O */
    wdot[11] -= q; /* CO */
    wdot[22] += q; /* CO2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[10][0] - 1)*dcdc_fac;
        J[45] -= dqdci;               /* dwdot[O]/d[H2] */
        J[51] -= dqdci;               /* dwdot[CO]/d[H2] */
        J[62] += dqdci;               /* dwdot[CO2]/d[H2] */
        /* d()/d[O] */
        dqdci =  + k_f*sc[11];
        J[205] -= dqdci;              /* dwdot[O]/d[O] */
        J[211] -= dqdci;              /* dwdot[CO]/d[O] */
        J[222] += dqdci;              /* dwdot[CO2]/d[O] */
        /* d()/d[H2O] */
        dqdci = (TB[10][1] - 1)*dcdc_fac;
        J[325] -= dqdci;              /* dwdot[O]/d[H2O] */
        J[331] -= dqdci;              /* dwdot[CO]/d[H2O] */
        J[342] += dqdci;              /* dwdot[CO2]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[10][2] - 1)*dcdc_fac + k_f*sc[5];
        J[445] -= dqdci;              /* dwdot[O]/d[CO] */
        J[451] -= dqdci;              /* dwdot[CO]/d[CO] */
        J[462] += dqdci;              /* dwdot[CO2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[10][3] - 1)*dcdc_fac - k_r;
        J[885] -= dqdci;              /* dwdot[O]/d[CO2] */
        J[891] -= dqdci;              /* dwdot[CO]/d[CO2] */
        J[902] += dqdci;              /* dwdot[CO2]/d[CO2] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = TB[10][0]*dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac + k_f*sc[11];
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = TB[10][1]*dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = TB[10][2]*dcdc_fac + k_f*sc[5];
        dqdc[12] = dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac;
        dqdc[22] = TB[10][3]*dcdc_fac - k_r;
        dqdc[23] = dcdc_fac;
        dqdc[24] = dcdc_fac;
        dqdc[25] = dcdc_fac;
        dqdc[26] = dcdc_fac;
        dqdc[27] = dcdc_fac;
        dqdc[28] = dcdc_fac;
        dqdc[29] = dcdc_fac;
        dqdc[30] = dcdc_fac;
        dqdc[31] = dcdc_fac;
        dqdc[32] = dcdc_fac;
        dqdc[33] = dcdc_fac;
        dqdc[34] = dcdc_fac;
        dqdc[35] = dcdc_fac;
        dqdc[36] = dcdc_fac;
        dqdc[37] = dcdc_fac;
        dqdc[38] = dcdc_fac;
        for (int k=0; k<39; k++) {
            J[40*k+5] -= dqdc[k];
            J[40*k+11] -= dqdc[k];
            J[40*k+22] += dqdc[k];
        }
    }
    J[1565] -= dqdT; /* dwdot[O]/dT */
    J[1571] -= dqdT; /* dwdot[CO]/dT */
    J[1582] += dqdT; /* dwdot[CO2]/dT */

    /*reaction 12: H2 + M <=> H + H + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[11][0] - 1)*sc[1] + (TB[11][1] - 1)*sc[8] + (TB[11][2] - 1)*sc[11] + (TB[11][3] - 1)*sc[22];
    /* forward */
    phi_f = sc[1];
    k_f = prefactor_units[11] * fwd_A[11]
                * exp(fwd_beta[11] * tc[0] - activation_units[11] * fwd_Ea[11] * invT);
    dlnkfdT = fwd_beta[11] * invT + activation_units[11] * fwd_Ea[11] * invT2;
    /* reverse */
    phi_r = pow(sc[0], 2.000000);
    Kc = refC * exp(-g_RT[0] - g_RT[0] + g_RT[1]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1]) + (2.000000*h_RT[0]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += 2 * q; /* H */
    wdot[1] -= q; /* H2 */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H] */
        dqdci =  - k_r*2.000000*sc[0];
        J[0] += 2 * dqdci;            /* dwdot[H]/d[H] */
        J[1] -= dqdci;                /* dwdot[H2]/d[H] */
        /* d()/d[H2] */
        dqdci = (TB[11][0] - 1)*q_nocor + k_f;
        J[40] += 2 * dqdci;           /* dwdot[H]/d[H2] */
        J[41] -= dqdci;               /* dwdot[H2]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[11][1] - 1)*q_nocor;
        J[320] += 2 * dqdci;          /* dwdot[H]/d[H2O] */
        J[321] -= dqdci;              /* dwdot[H2]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[11][2] - 1)*q_nocor;
        J[440] += 2 * dqdci;          /* dwdot[H]/d[CO] */
        J[441] -= dqdci;              /* dwdot[H2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[11][3] - 1)*q_nocor;
        J[880] += 2 * dqdci;          /* dwdot[H]/d[CO2] */
        J[881] -= dqdci;              /* dwdot[H2]/d[CO2] */
    }
    else {
        dqdc[0] = q_nocor - k_r*2.000000*sc[0];
        dqdc[1] = TB[11][0]*q_nocor + k_f;
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor;
        dqdc[5] = q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = TB[11][1]*q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = TB[11][2]*q_nocor;
        dqdc[12] = q_nocor;
        dqdc[13] = q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = q_nocor;
        dqdc[19] = q_nocor;
        dqdc[20] = q_nocor;
        dqdc[21] = q_nocor;
        dqdc[22] = TB[11][3]*q_nocor;
        dqdc[23] = q_nocor;
        dqdc[24] = q_nocor;
        dqdc[25] = q_nocor;
        dqdc[26] = q_nocor;
        dqdc[27] = q_nocor;
        dqdc[28] = q_nocor;
        dqdc[29] = q_nocor;
        dqdc[30] = q_nocor;
        dqdc[31] = q_nocor;
        dqdc[32] = q_nocor;
        dqdc[33] = q_nocor;
        dqdc[34] = q_nocor;
        dqdc[35] = q_nocor;
        dqdc[36] = q_nocor;
        dqdc[37] = q_nocor;
        dqdc[38] = q_nocor;
        for (int k=0; k<39; k++) {
            J[40*k+0] += 2 * dqdc[k];
            J[40*k+1] -= dqdc[k];
        }
    }
    J[1560] += 2 * dqdT; /* dwdot[H]/dT */
    J[1561] -= dqdT; /* dwdot[H2]/dT */

    /*reaction 13: O + O + M <=> O2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[12][0] - 1)*sc[1] + (TB[12][1] - 1)*sc[8] + (TB[12][2] - 1)*sc[11] + (TB[12][3] - 1)*sc[22];
    /* forward */
    phi_f = pow(sc[5], 2.000000);
    k_f = prefactor_units[12] * fwd_A[12]
                * exp(fwd_beta[12] * tc[0] - activation_units[12] * fwd_Ea[12] * invT);
    dlnkfdT = fwd_beta[12] * invT + activation_units[12] * fwd_Ea[12] * invT2;
    /* reverse */
    phi_r = sc[19];
    Kc = refCinv * exp(g_RT[5] + g_RT[5] - g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[5]) + (h_RT[19]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] -= 2 * q; /* O */
    wdot[19] += q; /* O2 */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[12][0] - 1)*q_nocor;
        J[45] += -2 * dqdci;          /* dwdot[O]/d[H2] */
        J[59] += dqdci;               /* dwdot[O2]/d[H2] */
        /* d()/d[O] */
        dqdci =  + k_f*2.000000*sc[5];
        J[205] += -2 * dqdci;         /* dwdot[O]/d[O] */
        J[219] += dqdci;              /* dwdot[O2]/d[O] */
        /* d()/d[H2O] */
        dqdci = (TB[12][1] - 1)*q_nocor;
        J[325] += -2 * dqdci;         /* dwdot[O]/d[H2O] */
        J[339] += dqdci;              /* dwdot[O2]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[12][2] - 1)*q_nocor;
        J[445] += -2 * dqdci;         /* dwdot[O]/d[CO] */
        J[459] += dqdci;              /* dwdot[O2]/d[CO] */
        /* d()/d[O2] */
        dqdci =  - k_r;
        J[765] += -2 * dqdci;         /* dwdot[O]/d[O2] */
        J[779] += dqdci;              /* dwdot[O2]/d[O2] */
        /* d()/d[CO2] */
        dqdci = (TB[12][3] - 1)*q_nocor;
        J[885] += -2 * dqdci;         /* dwdot[O]/d[CO2] */
        J[899] += dqdci;              /* dwdot[O2]/d[CO2] */
    }
    else {
        dqdc[0] = q_nocor;
        dqdc[1] = TB[12][0]*q_nocor;
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor;
        dqdc[5] = q_nocor + k_f*2.000000*sc[5];
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = TB[12][1]*q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = TB[12][2]*q_nocor;
        dqdc[12] = q_nocor;
        dqdc[13] = q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = q_nocor;
        dqdc[19] = q_nocor - k_r;
        dqdc[20] = q_nocor;
        dqdc[21] = q_nocor;
        dqdc[22] = TB[12][3]*q_nocor;
        dqdc[23] = q_nocor;
        dqdc[24] = q_nocor;
        dqdc[25] = q_nocor;
        dqdc[26] = q_nocor;
        dqdc[27] = q_nocor;
        dqdc[28] = q_nocor;
        dqdc[29] = q_nocor;
        dqdc[30] = q_nocor;
        dqdc[31] = q_nocor;
        dqdc[32] = q_nocor;
        dqdc[33] = q_nocor;
        dqdc[34] = q_nocor;
        dqdc[35] = q_nocor;
        dqdc[36] = q_nocor;
        dqdc[37] = q_nocor;
        dqdc[38] = q_nocor;
        for (int k=0; k<39; k++) {
            J[40*k+5] += -2 * dqdc[k];
            J[40*k+19] += dqdc[k];
        }
    }
    J[1565] += -2 * dqdT; /* dwdot[O]/dT */
    J[1579] += dqdT; /* dwdot[O2]/dT */

    /*reaction 14: O + H + M <=> OH + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[13][0] - 1)*sc[1] + (TB[13][1] - 1)*sc[8] + (TB[13][2] - 1)*sc[11] + (TB[13][3] - 1)*sc[22];
    /* forward */
    phi_f = sc[0]*sc[5];
    k_f = prefactor_units[13] * fwd_A[13]
                * exp(fwd_beta[13] * tc[0] - activation_units[13] * fwd_Ea[13] * invT);
    dlnkfdT = fwd_beta[13] * invT + activation_units[13] * fwd_Ea[13] * invT2;
    /* reverse */
    phi_r = sc[7];
    Kc = refCinv * exp(g_RT[0] + g_RT[5] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[5]) + (h_RT[7]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[5] -= q; /* O */
    wdot[7] += q; /* OH */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[5];
        J[0] -= dqdci;                /* dwdot[H]/d[H] */
        J[5] -= dqdci;                /* dwdot[O]/d[H] */
        J[7] += dqdci;                /* dwdot[OH]/d[H] */
        /* d()/d[H2] */
        dqdci = (TB[13][0] - 1)*q_nocor;
        J[40] -= dqdci;               /* dwdot[H]/d[H2] */
        J[45] -= dqdci;               /* dwdot[O]/d[H2] */
        J[47] += dqdci;               /* dwdot[OH]/d[H2] */
        /* d()/d[O] */
        dqdci =  + k_f*sc[0];
        J[200] -= dqdci;              /* dwdot[H]/d[O] */
        J[205] -= dqdci;              /* dwdot[O]/d[O] */
        J[207] += dqdci;              /* dwdot[OH]/d[O] */
        /* d()/d[OH] */
        dqdci =  - k_r;
        J[280] -= dqdci;              /* dwdot[H]/d[OH] */
        J[285] -= dqdci;              /* dwdot[O]/d[OH] */
        J[287] += dqdci;              /* dwdot[OH]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (TB[13][1] - 1)*q_nocor;
        J[320] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[325] -= dqdci;              /* dwdot[O]/d[H2O] */
        J[327] += dqdci;              /* dwdot[OH]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[13][2] - 1)*q_nocor;
        J[440] -= dqdci;              /* dwdot[H]/d[CO] */
        J[445] -= dqdci;              /* dwdot[O]/d[CO] */
        J[447] += dqdci;              /* dwdot[OH]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[13][3] - 1)*q_nocor;
        J[880] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[885] -= dqdci;              /* dwdot[O]/d[CO2] */
        J[887] += dqdci;              /* dwdot[OH]/d[CO2] */
    }
    else {
        dqdc[0] = q_nocor + k_f*sc[5];
        dqdc[1] = TB[13][0]*q_nocor;
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor;
        dqdc[5] = q_nocor + k_f*sc[0];
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor - k_r;
        dqdc[8] = TB[13][1]*q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = TB[13][2]*q_nocor;
        dqdc[12] = q_nocor;
        dqdc[13] = q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = q_nocor;
        dqdc[19] = q_nocor;
        dqdc[20] = q_nocor;
        dqdc[21] = q_nocor;
        dqdc[22] = TB[13][3]*q_nocor;
        dqdc[23] = q_nocor;
        dqdc[24] = q_nocor;
        dqdc[25] = q_nocor;
        dqdc[26] = q_nocor;
        dqdc[27] = q_nocor;
        dqdc[28] = q_nocor;
        dqdc[29] = q_nocor;
        dqdc[30] = q_nocor;
        dqdc[31] = q_nocor;
        dqdc[32] = q_nocor;
        dqdc[33] = q_nocor;
        dqdc[34] = q_nocor;
        dqdc[35] = q_nocor;
        dqdc[36] = q_nocor;
        dqdc[37] = q_nocor;
        dqdc[38] = q_nocor;
        for (int k=0; k<39; k++) {
            J[40*k+0] -= dqdc[k];
            J[40*k+5] -= dqdc[k];
            J[40*k+7] += dqdc[k];
        }
    }
    J[1560] -= dqdT; /* dwdot[H]/dT */
    J[1565] -= dqdT; /* dwdot[O]/dT */
    J[1567] += dqdT; /* dwdot[OH]/dT */

    /*reaction 15: H + OH + M <=> H2O + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[14][0] - 1)*sc[1] + (TB[14][1] - 1)*sc[8] + (TB[14][2] - 1)*sc[11] + (TB[14][3] - 1)*sc[22];
    /* forward */
    phi_f = sc[0]*sc[7];
    k_f = prefactor_units[14] * fwd_A[14]
                * exp(fwd_beta[14] * tc[0] - activation_units[14] * fwd_Ea[14] * invT);
    dlnkfdT = fwd_beta[14] * invT + activation_units[14] * fwd_Ea[14] * invT2;
    /* reverse */
    phi_r = sc[8];
    Kc = refCinv * exp(g_RT[0] + g_RT[7] - g_RT[8]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[7]) + (h_RT[8]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[7] -= q; /* OH */
    wdot[8] += q; /* H2O */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[7];
        J[0] -= dqdci;                /* dwdot[H]/d[H] */
        J[7] -= dqdci;                /* dwdot[OH]/d[H] */
        J[8] += dqdci;                /* dwdot[H2O]/d[H] */
        /* d()/d[H2] */
        dqdci = (TB[14][0] - 1)*q_nocor;
        J[40] -= dqdci;               /* dwdot[H]/d[H2] */
        J[47] -= dqdci;               /* dwdot[OH]/d[H2] */
        J[48] += dqdci;               /* dwdot[H2O]/d[H2] */
        /* d()/d[OH] */
        dqdci =  + k_f*sc[0];
        J[280] -= dqdci;              /* dwdot[H]/d[OH] */
        J[287] -= dqdci;              /* dwdot[OH]/d[OH] */
        J[288] += dqdci;              /* dwdot[H2O]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (TB[14][1] - 1)*q_nocor - k_r;
        J[320] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[327] -= dqdci;              /* dwdot[OH]/d[H2O] */
        J[328] += dqdci;              /* dwdot[H2O]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[14][2] - 1)*q_nocor;
        J[440] -= dqdci;              /* dwdot[H]/d[CO] */
        J[447] -= dqdci;              /* dwdot[OH]/d[CO] */
        J[448] += dqdci;              /* dwdot[H2O]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[14][3] - 1)*q_nocor;
        J[880] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[887] -= dqdci;              /* dwdot[OH]/d[CO2] */
        J[888] += dqdci;              /* dwdot[H2O]/d[CO2] */
    }
    else {
        dqdc[0] = q_nocor + k_f*sc[7];
        dqdc[1] = TB[14][0]*q_nocor;
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor;
        dqdc[5] = q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor + k_f*sc[0];
        dqdc[8] = TB[14][1]*q_nocor - k_r;
        dqdc[9] = q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = TB[14][2]*q_nocor;
        dqdc[12] = q_nocor;
        dqdc[13] = q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = q_nocor;
        dqdc[19] = q_nocor;
        dqdc[20] = q_nocor;
        dqdc[21] = q_nocor;
        dqdc[22] = TB[14][3]*q_nocor;
        dqdc[23] = q_nocor;
        dqdc[24] = q_nocor;
        dqdc[25] = q_nocor;
        dqdc[26] = q_nocor;
        dqdc[27] = q_nocor;
        dqdc[28] = q_nocor;
        dqdc[29] = q_nocor;
        dqdc[30] = q_nocor;
        dqdc[31] = q_nocor;
        dqdc[32] = q_nocor;
        dqdc[33] = q_nocor;
        dqdc[34] = q_nocor;
        dqdc[35] = q_nocor;
        dqdc[36] = q_nocor;
        dqdc[37] = q_nocor;
        dqdc[38] = q_nocor;
        for (int k=0; k<39; k++) {
            J[40*k+0] -= dqdc[k];
            J[40*k+7] -= dqdc[k];
            J[40*k+8] += dqdc[k];
        }
    }
    J[1560] -= dqdT; /* dwdot[H]/dT */
    J[1567] -= dqdT; /* dwdot[OH]/dT */
    J[1568] += dqdT; /* dwdot[H2O]/dT */

    /*reaction 16: HCO + M <=> H + CO + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[15][0] - 1)*sc[1] + (TB[15][1] - 1)*sc[8] + (TB[15][2] - 1)*sc[11] + (TB[15][3] - 1)*sc[22];
    /* forward */
    phi_f = sc[13];
    k_f = prefactor_units[15] * fwd_A[15]
                * exp(fwd_beta[15] * tc[0] - activation_units[15] * fwd_Ea[15] * invT);
    dlnkfdT = fwd_beta[15] * invT + activation_units[15] * fwd_Ea[15] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[11];
    Kc = refC * exp(-g_RT[0] - g_RT[11] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[13]) + (h_RT[0] + h_RT[11]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[11] += q; /* CO */
    wdot[13] -= q; /* HCO */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H] */
        dqdci =  - k_r*sc[11];
        J[0] += dqdci;                /* dwdot[H]/d[H] */
        J[11] += dqdci;               /* dwdot[CO]/d[H] */
        J[13] -= dqdci;               /* dwdot[HCO]/d[H] */
        /* d()/d[H2] */
        dqdci = (TB[15][0] - 1)*q_nocor;
        J[40] += dqdci;               /* dwdot[H]/d[H2] */
        J[51] += dqdci;               /* dwdot[CO]/d[H2] */
        J[53] -= dqdci;               /* dwdot[HCO]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[15][1] - 1)*q_nocor;
        J[320] += dqdci;              /* dwdot[H]/d[H2O] */
        J[331] += dqdci;              /* dwdot[CO]/d[H2O] */
        J[333] -= dqdci;              /* dwdot[HCO]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[15][2] - 1)*q_nocor - k_r*sc[0];
        J[440] += dqdci;              /* dwdot[H]/d[CO] */
        J[451] += dqdci;              /* dwdot[CO]/d[CO] */
        J[453] -= dqdci;              /* dwdot[HCO]/d[CO] */
        /* d()/d[HCO] */
        dqdci =  + k_f;
        J[520] += dqdci;              /* dwdot[H]/d[HCO] */
        J[531] += dqdci;              /* dwdot[CO]/d[HCO] */
        J[533] -= dqdci;              /* dwdot[HCO]/d[HCO] */
        /* d()/d[CO2] */
        dqdci = (TB[15][3] - 1)*q_nocor;
        J[880] += dqdci;              /* dwdot[H]/d[CO2] */
        J[891] += dqdci;              /* dwdot[CO]/d[CO2] */
        J[893] -= dqdci;              /* dwdot[HCO]/d[CO2] */
    }
    else {
        dqdc[0] = q_nocor - k_r*sc[11];
        dqdc[1] = TB[15][0]*q_nocor;
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor;
        dqdc[5] = q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = TB[15][1]*q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = TB[15][2]*q_nocor - k_r*sc[0];
        dqdc[12] = q_nocor;
        dqdc[13] = q_nocor + k_f;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = q_nocor;
        dqdc[19] = q_nocor;
        dqdc[20] = q_nocor;
        dqdc[21] = q_nocor;
        dqdc[22] = TB[15][3]*q_nocor;
        dqdc[23] = q_nocor;
        dqdc[24] = q_nocor;
        dqdc[25] = q_nocor;
        dqdc[26] = q_nocor;
        dqdc[27] = q_nocor;
        dqdc[28] = q_nocor;
        dqdc[29] = q_nocor;
        dqdc[30] = q_nocor;
        dqdc[31] = q_nocor;
        dqdc[32] = q_nocor;
        dqdc[33] = q_nocor;
        dqdc[34] = q_nocor;
        dqdc[35] = q_nocor;
        dqdc[36] = q_nocor;
        dqdc[37] = q_nocor;
        dqdc[38] = q_nocor;
        for (int k=0; k<39; k++) {
            J[40*k+0] += dqdc[k];
            J[40*k+11] += dqdc[k];
            J[40*k+13] -= dqdc[k];
        }
    }
    J[1560] += dqdT; /* dwdot[H]/dT */
    J[1571] += dqdT; /* dwdot[CO]/dT */
    J[1573] -= dqdT; /* dwdot[HCO]/dT */

    /*reaction 17: CH2O + M <=> HCO + H + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[16][0] - 1)*sc[1] + (TB[16][1] - 1)*sc[8] + (TB[16][2] - 1)*sc[11] + (TB[16][3] - 1)*sc[22];
    /* forward */
    phi_f = sc[15];
    k_f = prefactor_units[16] * fwd_A[16]
                * exp(fwd_beta[16] * tc[0] - activation_units[16] * fwd_Ea[16] * invT);
    dlnkfdT = fwd_beta[16] * invT + activation_units[16] * fwd_Ea[16] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[13];
    Kc = refC * exp(-g_RT[0] - g_RT[13] + g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[15]) + (h_RT[0] + h_RT[13]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[13] += q; /* HCO */
    wdot[15] -= q; /* CH2O */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H] */
        dqdci =  - k_r*sc[13];
        J[0] += dqdci;                /* dwdot[H]/d[H] */
        J[13] += dqdci;               /* dwdot[HCO]/d[H] */
        J[15] -= dqdci;               /* dwdot[CH2O]/d[H] */
        /* d()/d[H2] */
        dqdci = (TB[16][0] - 1)*q_nocor;
        J[40] += dqdci;               /* dwdot[H]/d[H2] */
        J[53] += dqdci;               /* dwdot[HCO]/d[H2] */
        J[55] -= dqdci;               /* dwdot[CH2O]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[16][1] - 1)*q_nocor;
        J[320] += dqdci;              /* dwdot[H]/d[H2O] */
        J[333] += dqdci;              /* dwdot[HCO]/d[H2O] */
        J[335] -= dqdci;              /* dwdot[CH2O]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[16][2] - 1)*q_nocor;
        J[440] += dqdci;              /* dwdot[H]/d[CO] */
        J[453] += dqdci;              /* dwdot[HCO]/d[CO] */
        J[455] -= dqdci;              /* dwdot[CH2O]/d[CO] */
        /* d()/d[HCO] */
        dqdci =  - k_r*sc[0];
        J[520] += dqdci;              /* dwdot[H]/d[HCO] */
        J[533] += dqdci;              /* dwdot[HCO]/d[HCO] */
        J[535] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
        /* d()/d[CH2O] */
        dqdci =  + k_f;
        J[600] += dqdci;              /* dwdot[H]/d[CH2O] */
        J[613] += dqdci;              /* dwdot[HCO]/d[CH2O] */
        J[615] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
        /* d()/d[CO2] */
        dqdci = (TB[16][3] - 1)*q_nocor;
        J[880] += dqdci;              /* dwdot[H]/d[CO2] */
        J[893] += dqdci;              /* dwdot[HCO]/d[CO2] */
        J[895] -= dqdci;              /* dwdot[CH2O]/d[CO2] */
    }
    else {
        dqdc[0] = q_nocor - k_r*sc[13];
        dqdc[1] = TB[16][0]*q_nocor;
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor;
        dqdc[5] = q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = TB[16][1]*q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = TB[16][2]*q_nocor;
        dqdc[12] = q_nocor;
        dqdc[13] = q_nocor - k_r*sc[0];
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor + k_f;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = q_nocor;
        dqdc[19] = q_nocor;
        dqdc[20] = q_nocor;
        dqdc[21] = q_nocor;
        dqdc[22] = TB[16][3]*q_nocor;
        dqdc[23] = q_nocor;
        dqdc[24] = q_nocor;
        dqdc[25] = q_nocor;
        dqdc[26] = q_nocor;
        dqdc[27] = q_nocor;
        dqdc[28] = q_nocor;
        dqdc[29] = q_nocor;
        dqdc[30] = q_nocor;
        dqdc[31] = q_nocor;
        dqdc[32] = q_nocor;
        dqdc[33] = q_nocor;
        dqdc[34] = q_nocor;
        dqdc[35] = q_nocor;
        dqdc[36] = q_nocor;
        dqdc[37] = q_nocor;
        dqdc[38] = q_nocor;
        for (int k=0; k<39; k++) {
            J[40*k+0] += dqdc[k];
            J[40*k+13] += dqdc[k];
            J[40*k+15] -= dqdc[k];
        }
    }
    J[1560] += dqdT; /* dwdot[H]/dT */
    J[1573] += dqdT; /* dwdot[HCO]/dT */
    J[1575] -= dqdT; /* dwdot[CH2O]/dT */

    /*reaction 18: CH2O + M <=> CO + H2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[17][0] - 1)*sc[1] + (TB[17][1] - 1)*sc[8] + (TB[17][2] - 1)*sc[11] + (TB[17][3] - 1)*sc[22];
    /* forward */
    phi_f = sc[15];
    k_f = prefactor_units[17] * fwd_A[17]
                * exp(fwd_beta[17] * tc[0] - activation_units[17] * fwd_Ea[17] * invT);
    dlnkfdT = fwd_beta[17] * invT + activation_units[17] * fwd_Ea[17] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[11];
    Kc = refC * exp(-g_RT[1] - g_RT[11] + g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[15]) + (h_RT[1] + h_RT[11]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H2 */
    wdot[11] += q; /* CO */
    wdot[15] -= q; /* CH2O */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[17][0] - 1)*q_nocor - k_r*sc[11];
        J[41] += dqdci;               /* dwdot[H2]/d[H2] */
        J[51] += dqdci;               /* dwdot[CO]/d[H2] */
        J[55] -= dqdci;               /* dwdot[CH2O]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[17][1] - 1)*q_nocor;
        J[321] += dqdci;              /* dwdot[H2]/d[H2O] */
        J[331] += dqdci;              /* dwdot[CO]/d[H2O] */
        J[335] -= dqdci;              /* dwdot[CH2O]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[17][2] - 1)*q_nocor - k_r*sc[1];
        J[441] += dqdci;              /* dwdot[H2]/d[CO] */
        J[451] += dqdci;              /* dwdot[CO]/d[CO] */
        J[455] -= dqdci;              /* dwdot[CH2O]/d[CO] */
        /* d()/d[CH2O] */
        dqdci =  + k_f;
        J[601] += dqdci;              /* dwdot[H2]/d[CH2O] */
        J[611] += dqdci;              /* dwdot[CO]/d[CH2O] */
        J[615] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
        /* d()/d[CO2] */
        dqdci = (TB[17][3] - 1)*q_nocor;
        J[881] += dqdci;              /* dwdot[H2]/d[CO2] */
        J[891] += dqdci;              /* dwdot[CO]/d[CO2] */
        J[895] -= dqdci;              /* dwdot[CH2O]/d[CO2] */
    }
    else {
        dqdc[0] = q_nocor;
        dqdc[1] = TB[17][0]*q_nocor - k_r*sc[11];
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor;
        dqdc[5] = q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = TB[17][1]*q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = TB[17][2]*q_nocor - k_r*sc[1];
        dqdc[12] = q_nocor;
        dqdc[13] = q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor + k_f;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = q_nocor;
        dqdc[19] = q_nocor;
        dqdc[20] = q_nocor;
        dqdc[21] = q_nocor;
        dqdc[22] = TB[17][3]*q_nocor;
        dqdc[23] = q_nocor;
        dqdc[24] = q_nocor;
        dqdc[25] = q_nocor;
        dqdc[26] = q_nocor;
        dqdc[27] = q_nocor;
        dqdc[28] = q_nocor;
        dqdc[29] = q_nocor;
        dqdc[30] = q_nocor;
        dqdc[31] = q_nocor;
        dqdc[32] = q_nocor;
        dqdc[33] = q_nocor;
        dqdc[34] = q_nocor;
        dqdc[35] = q_nocor;
        dqdc[36] = q_nocor;
        dqdc[37] = q_nocor;
        dqdc[38] = q_nocor;
        for (int k=0; k<39; k++) {
            J[40*k+1] += dqdc[k];
            J[40*k+11] += dqdc[k];
            J[40*k+15] -= dqdc[k];
        }
    }
    J[1561] += dqdT; /* dwdot[H2]/dT */
    J[1571] += dqdT; /* dwdot[CO]/dT */
    J[1575] -= dqdT; /* dwdot[CH2O]/dT */

    /*reaction 19: CH2OH + M <=> CH2O + H + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture;
    /* forward */
    phi_f = sc[17];
    k_f = prefactor_units[18] * fwd_A[18]
                * exp(fwd_beta[18] * tc[0] - activation_units[18] * fwd_Ea[18] * invT);
    dlnkfdT = fwd_beta[18] * invT + activation_units[18] * fwd_Ea[18] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[15];
    Kc = refC * exp(-g_RT[0] - g_RT[15] + g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[17]) + (h_RT[0] + h_RT[15]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[15] += q; /* CH2O */
    wdot[17] -= q; /* CH2OH */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H] */
        dqdci =  - k_r*sc[15];
        J[0] += dqdci;                /* dwdot[H]/d[H] */
        J[15] += dqdci;               /* dwdot[CH2O]/d[H] */
        J[17] -= dqdci;               /* dwdot[CH2OH]/d[H] */
        /* d()/d[CH2O] */
        dqdci =  - k_r*sc[0];
        J[600] += dqdci;              /* dwdot[H]/d[CH2O] */
        J[615] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
        J[617] -= dqdci;              /* dwdot[CH2OH]/d[CH2O] */
        /* d()/d[CH2OH] */
        dqdci =  + k_f;
        J[680] += dqdci;              /* dwdot[H]/d[CH2OH] */
        J[695] += dqdci;              /* dwdot[CH2O]/d[CH2OH] */
        J[697] -= dqdci;              /* dwdot[CH2OH]/d[CH2OH] */
    }
    else {
        dqdc[0] = q_nocor - k_r*sc[15];
        dqdc[1] = q_nocor;
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor;
        dqdc[5] = q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = q_nocor;
        dqdc[12] = q_nocor;
        dqdc[13] = q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor - k_r*sc[0];
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor + k_f;
        dqdc[18] = q_nocor;
        dqdc[19] = q_nocor;
        dqdc[20] = q_nocor;
        dqdc[21] = q_nocor;
        dqdc[22] = q_nocor;
        dqdc[23] = q_nocor;
        dqdc[24] = q_nocor;
        dqdc[25] = q_nocor;
        dqdc[26] = q_nocor;
        dqdc[27] = q_nocor;
        dqdc[28] = q_nocor;
        dqdc[29] = q_nocor;
        dqdc[30] = q_nocor;
        dqdc[31] = q_nocor;
        dqdc[32] = q_nocor;
        dqdc[33] = q_nocor;
        dqdc[34] = q_nocor;
        dqdc[35] = q_nocor;
        dqdc[36] = q_nocor;
        dqdc[37] = q_nocor;
        dqdc[38] = q_nocor;
        for (int k=0; k<39; k++) {
            J[40*k+0] += dqdc[k];
            J[40*k+15] += dqdc[k];
            J[40*k+17] -= dqdc[k];
        }
    }
    J[1560] += dqdT; /* dwdot[H]/dT */
    J[1575] += dqdT; /* dwdot[CH2O]/dT */
    J[1577] -= dqdT; /* dwdot[CH2OH]/dT */

    /*reaction 20: CH3O + M <=> CH2O + H + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture;
    /* forward */
    phi_f = sc[18];
    k_f = prefactor_units[19] * fwd_A[19]
                * exp(fwd_beta[19] * tc[0] - activation_units[19] * fwd_Ea[19] * invT);
    dlnkfdT = fwd_beta[19] * invT + activation_units[19] * fwd_Ea[19] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[15];
    Kc = refC * exp(-g_RT[0] - g_RT[15] + g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[18]) + (h_RT[0] + h_RT[15]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[15] += q; /* CH2O */
    wdot[18] -= q; /* CH3O */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H] */
        dqdci =  - k_r*sc[15];
        J[0] += dqdci;                /* dwdot[H]/d[H] */
        J[15] += dqdci;               /* dwdot[CH2O]/d[H] */
        J[18] -= dqdci;               /* dwdot[CH3O]/d[H] */
        /* d()/d[CH2O] */
        dqdci =  - k_r*sc[0];
        J[600] += dqdci;              /* dwdot[H]/d[CH2O] */
        J[615] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
        J[618] -= dqdci;              /* dwdot[CH3O]/d[CH2O] */
        /* d()/d[CH3O] */
        dqdci =  + k_f;
        J[720] += dqdci;              /* dwdot[H]/d[CH3O] */
        J[735] += dqdci;              /* dwdot[CH2O]/d[CH3O] */
        J[738] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    }
    else {
        dqdc[0] = q_nocor - k_r*sc[15];
        dqdc[1] = q_nocor;
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor;
        dqdc[5] = q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = q_nocor;
        dqdc[12] = q_nocor;
        dqdc[13] = q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor - k_r*sc[0];
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = q_nocor + k_f;
        dqdc[19] = q_nocor;
        dqdc[20] = q_nocor;
        dqdc[21] = q_nocor;
        dqdc[22] = q_nocor;
        dqdc[23] = q_nocor;
        dqdc[24] = q_nocor;
        dqdc[25] = q_nocor;
        dqdc[26] = q_nocor;
        dqdc[27] = q_nocor;
        dqdc[28] = q_nocor;
        dqdc[29] = q_nocor;
        dqdc[30] = q_nocor;
        dqdc[31] = q_nocor;
        dqdc[32] = q_nocor;
        dqdc[33] = q_nocor;
        dqdc[34] = q_nocor;
        dqdc[35] = q_nocor;
        dqdc[36] = q_nocor;
        dqdc[37] = q_nocor;
        dqdc[38] = q_nocor;
        for (int k=0; k<39; k++) {
            J[40*k+0] += dqdc[k];
            J[40*k+15] += dqdc[k];
            J[40*k+18] -= dqdc[k];
        }
    }
    J[1560] += dqdT; /* dwdot[H]/dT */
    J[1575] += dqdT; /* dwdot[CH2O]/dT */
    J[1578] -= dqdT; /* dwdot[CH3O]/dT */

    /*reaction 21: CH2(S) + M <=> CH2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[20][0] - 1)*sc[8] + (TB[20][1] - 1)*sc[11] + (TB[20][2] - 1)*sc[22];
    /* forward */
    phi_f = sc[3];
    k_f = prefactor_units[20] * fwd_A[20]
                * exp(fwd_beta[20] * tc[0] - activation_units[20] * fwd_Ea[20] * invT);
    dlnkfdT = fwd_beta[20] * invT + activation_units[20] * fwd_Ea[20] * invT2;
    /* reverse */
    phi_r = sc[2];
    Kc = exp(-g_RT[2] + g_RT[3]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3]) + (h_RT[2]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* CH2 */
    wdot[3] -= q; /* CH2(S) */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[CH2] */
        dqdci =  - k_r;
        J[82] += dqdci;               /* dwdot[CH2]/d[CH2] */
        J[83] -= dqdci;               /* dwdot[CH2(S)]/d[CH2] */
        /* d()/d[CH2(S)] */
        dqdci =  + k_f;
        J[122] += dqdci;              /* dwdot[CH2]/d[CH2(S)] */
        J[123] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
        /* d()/d[H2O] */
        dqdci = (TB[20][0] - 1)*q_nocor;
        J[322] += dqdci;              /* dwdot[CH2]/d[H2O] */
        J[323] -= dqdci;              /* dwdot[CH2(S)]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[20][1] - 1)*q_nocor;
        J[442] += dqdci;              /* dwdot[CH2]/d[CO] */
        J[443] -= dqdci;              /* dwdot[CH2(S)]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[20][2] - 1)*q_nocor;
        J[882] += dqdci;              /* dwdot[CH2]/d[CO2] */
        J[883] -= dqdci;              /* dwdot[CH2(S)]/d[CO2] */
    }
    else {
        dqdc[0] = q_nocor;
        dqdc[1] = q_nocor;
        dqdc[2] = q_nocor - k_r;
        dqdc[3] = q_nocor + k_f;
        dqdc[4] = q_nocor;
        dqdc[5] = q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = TB[20][0]*q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = TB[20][1]*q_nocor;
        dqdc[12] = q_nocor;
        dqdc[13] = q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = q_nocor;
        dqdc[19] = q_nocor;
        dqdc[20] = q_nocor;
        dqdc[21] = q_nocor;
        dqdc[22] = TB[20][2]*q_nocor;
        dqdc[23] = q_nocor;
        dqdc[24] = q_nocor;
        dqdc[25] = q_nocor;
        dqdc[26] = q_nocor;
        dqdc[27] = q_nocor;
        dqdc[28] = q_nocor;
        dqdc[29] = q_nocor;
        dqdc[30] = q_nocor;
        dqdc[31] = q_nocor;
        dqdc[32] = q_nocor;
        dqdc[33] = q_nocor;
        dqdc[34] = q_nocor;
        dqdc[35] = q_nocor;
        dqdc[36] = q_nocor;
        dqdc[37] = q_nocor;
        dqdc[38] = q_nocor;
        for (int k=0; k<39; k++) {
            J[40*k+2] += dqdc[k];
            J[40*k+3] -= dqdc[k];
        }
    }
    J[1562] += dqdT; /* dwdot[CH2]/dT */
    J[1563] -= dqdT; /* dwdot[CH2(S)]/dT */

    /*reaction 22: HCOOH + M <=> CO + H2O + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture;
    /* forward */
    phi_f = sc[25];
    k_f = prefactor_units[21] * fwd_A[21]
                * exp(fwd_beta[21] * tc[0] - activation_units[21] * fwd_Ea[21] * invT);
    dlnkfdT = fwd_beta[21] * invT + activation_units[21] * fwd_Ea[21] * invT2;
    /* reverse */
    phi_r = sc[8]*sc[11];
    Kc = refC * exp(-g_RT[8] - g_RT[11] + g_RT[25]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[25]) + (h_RT[8] + h_RT[11]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[8] += q; /* H2O */
    wdot[11] += q; /* CO */
    wdot[25] -= q; /* HCOOH */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2O] */
        dqdci =  - k_r*sc[11];
        J[328] += dqdci;              /* dwdot[H2O]/d[H2O] */
        J[331] += dqdci;              /* dwdot[CO]/d[H2O] */
        J[345] -= dqdci;              /* dwdot[HCOOH]/d[H2O] */
        /* d()/d[CO] */
        dqdci =  - k_r*sc[8];
        J[448] += dqdci;              /* dwdot[H2O]/d[CO] */
        J[451] += dqdci;              /* dwdot[CO]/d[CO] */
        J[465] -= dqdci;              /* dwdot[HCOOH]/d[CO] */
        /* d()/d[HCOOH] */
        dqdci =  + k_f;
        J[1008] += dqdci;             /* dwdot[H2O]/d[HCOOH] */
        J[1011] += dqdci;             /* dwdot[CO]/d[HCOOH] */
        J[1025] -= dqdci;             /* dwdot[HCOOH]/d[HCOOH] */
    }
    else {
        dqdc[0] = q_nocor;
        dqdc[1] = q_nocor;
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor;
        dqdc[5] = q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor - k_r*sc[11];
        dqdc[9] = q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = q_nocor - k_r*sc[8];
        dqdc[12] = q_nocor;
        dqdc[13] = q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = q_nocor;
        dqdc[19] = q_nocor;
        dqdc[20] = q_nocor;
        dqdc[21] = q_nocor;
        dqdc[22] = q_nocor;
        dqdc[23] = q_nocor;
        dqdc[24] = q_nocor;
        dqdc[25] = q_nocor + k_f;
        dqdc[26] = q_nocor;
        dqdc[27] = q_nocor;
        dqdc[28] = q_nocor;
        dqdc[29] = q_nocor;
        dqdc[30] = q_nocor;
        dqdc[31] = q_nocor;
        dqdc[32] = q_nocor;
        dqdc[33] = q_nocor;
        dqdc[34] = q_nocor;
        dqdc[35] = q_nocor;
        dqdc[36] = q_nocor;
        dqdc[37] = q_nocor;
        dqdc[38] = q_nocor;
        for (int k=0; k<39; k++) {
            J[40*k+8] += dqdc[k];
            J[40*k+11] += dqdc[k];
            J[40*k+25] -= dqdc[k];
        }
    }
    J[1568] += dqdT; /* dwdot[H2O]/dT */
    J[1571] += dqdT; /* dwdot[CO]/dT */
    J[1585] -= dqdT; /* dwdot[HCOOH]/dT */

    /*reaction 23: HCOOH + M <=> CO2 + H2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture;
    /* forward */
    phi_f = sc[25];
    k_f = prefactor_units[22] * fwd_A[22]
                * exp(fwd_beta[22] * tc[0] - activation_units[22] * fwd_Ea[22] * invT);
    dlnkfdT = fwd_beta[22] * invT + activation_units[22] * fwd_Ea[22] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[22];
    Kc = refC * exp(-g_RT[1] - g_RT[22] + g_RT[25]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[25]) + (h_RT[1] + h_RT[22]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H2 */
    wdot[22] += q; /* CO2 */
    wdot[25] -= q; /* HCOOH */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci =  - k_r*sc[22];
        J[41] += dqdci;               /* dwdot[H2]/d[H2] */
        J[62] += dqdci;               /* dwdot[CO2]/d[H2] */
        J[65] -= dqdci;               /* dwdot[HCOOH]/d[H2] */
        /* d()/d[CO2] */
        dqdci =  - k_r*sc[1];
        J[881] += dqdci;              /* dwdot[H2]/d[CO2] */
        J[902] += dqdci;              /* dwdot[CO2]/d[CO2] */
        J[905] -= dqdci;              /* dwdot[HCOOH]/d[CO2] */
        /* d()/d[HCOOH] */
        dqdci =  + k_f;
        J[1001] += dqdci;             /* dwdot[H2]/d[HCOOH] */
        J[1022] += dqdci;             /* dwdot[CO2]/d[HCOOH] */
        J[1025] -= dqdci;             /* dwdot[HCOOH]/d[HCOOH] */
    }
    else {
        dqdc[0] = q_nocor;
        dqdc[1] = q_nocor - k_r*sc[22];
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor;
        dqdc[5] = q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = q_nocor;
        dqdc[12] = q_nocor;
        dqdc[13] = q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = q_nocor;
        dqdc[19] = q_nocor;
        dqdc[20] = q_nocor;
        dqdc[21] = q_nocor;
        dqdc[22] = q_nocor - k_r*sc[1];
        dqdc[23] = q_nocor;
        dqdc[24] = q_nocor;
        dqdc[25] = q_nocor + k_f;
        dqdc[26] = q_nocor;
        dqdc[27] = q_nocor;
        dqdc[28] = q_nocor;
        dqdc[29] = q_nocor;
        dqdc[30] = q_nocor;
        dqdc[31] = q_nocor;
        dqdc[32] = q_nocor;
        dqdc[33] = q_nocor;
        dqdc[34] = q_nocor;
        dqdc[35] = q_nocor;
        dqdc[36] = q_nocor;
        dqdc[37] = q_nocor;
        dqdc[38] = q_nocor;
        for (int k=0; k<39; k++) {
            J[40*k+1] += dqdc[k];
            J[40*k+22] += dqdc[k];
            J[40*k+25] -= dqdc[k];
        }
    }
    J[1561] += dqdT; /* dwdot[H2]/dT */
    J[1582] += dqdT; /* dwdot[CO2]/dT */
    J[1585] -= dqdT; /* dwdot[HCOOH]/dT */

    /*reaction 24: H + O2 <=> O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[19];
    k_f = prefactor_units[23] * fwd_A[23]
                * exp(fwd_beta[23] * tc[0] - activation_units[23] * fwd_Ea[23] * invT);
    dlnkfdT = fwd_beta[23] * invT + activation_units[23] * fwd_Ea[23] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[7];
    Kc = exp(g_RT[0] - g_RT[5] - g_RT[7] + g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[19]) + (h_RT[5] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[5] += q; /* O */
    wdot[7] += q; /* OH */
    wdot[19] -= q; /* O2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[19];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[5] += dqdci;                /* dwdot[O]/d[H] */
    J[7] += dqdci;                /* dwdot[OH]/d[H] */
    J[19] -= dqdci;               /* dwdot[O2]/d[H] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[7];
    J[200] -= dqdci;              /* dwdot[H]/d[O] */
    J[205] += dqdci;              /* dwdot[O]/d[O] */
    J[207] += dqdci;              /* dwdot[OH]/d[O] */
    J[219] -= dqdci;              /* dwdot[O2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[5];
    J[280] -= dqdci;              /* dwdot[H]/d[OH] */
    J[285] += dqdci;              /* dwdot[O]/d[OH] */
    J[287] += dqdci;              /* dwdot[OH]/d[OH] */
    J[299] -= dqdci;              /* dwdot[O2]/d[OH] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[0];
    J[760] -= dqdci;              /* dwdot[H]/d[O2] */
    J[765] += dqdci;              /* dwdot[O]/d[O2] */
    J[767] += dqdci;              /* dwdot[OH]/d[O2] */
    J[779] -= dqdci;              /* dwdot[O2]/d[O2] */
    /* d()/dT */
    J[1560] -= dqdT;              /* dwdot[H]/dT */
    J[1565] += dqdT;              /* dwdot[O]/dT */
    J[1567] += dqdT;              /* dwdot[OH]/dT */
    J[1579] -= dqdT;              /* dwdot[O2]/dT */

    /*reaction 25: O + H2 <=> H + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[5];
    k_f = prefactor_units[24] * fwd_A[24]
                * exp(fwd_beta[24] * tc[0] - activation_units[24] * fwd_Ea[24] * invT);
    dlnkfdT = fwd_beta[24] * invT + activation_units[24] * fwd_Ea[24] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[7];
    Kc = exp(-g_RT[0] + g_RT[1] + g_RT[5] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[5]) + (h_RT[0] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[1] -= q; /* H2 */
    wdot[5] -= q; /* O */
    wdot[7] += q; /* OH */
    /* d()/d[H] */
    dqdci =  - k_r*sc[7];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[1] -= dqdci;                /* dwdot[H2]/d[H] */
    J[5] -= dqdci;                /* dwdot[O]/d[H] */
    J[7] += dqdci;                /* dwdot[OH]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[5];
    J[40] += dqdci;               /* dwdot[H]/d[H2] */
    J[41] -= dqdci;               /* dwdot[H2]/d[H2] */
    J[45] -= dqdci;               /* dwdot[O]/d[H2] */
    J[47] += dqdci;               /* dwdot[OH]/d[H2] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[1];
    J[200] += dqdci;              /* dwdot[H]/d[O] */
    J[201] -= dqdci;              /* dwdot[H2]/d[O] */
    J[205] -= dqdci;              /* dwdot[O]/d[O] */
    J[207] += dqdci;              /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[0];
    J[280] += dqdci;              /* dwdot[H]/d[OH] */
    J[281] -= dqdci;              /* dwdot[H2]/d[OH] */
    J[285] -= dqdci;              /* dwdot[O]/d[OH] */
    J[287] += dqdci;              /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[1560] += dqdT;              /* dwdot[H]/dT */
    J[1561] -= dqdT;              /* dwdot[H2]/dT */
    J[1565] -= dqdT;              /* dwdot[O]/dT */
    J[1567] += dqdT;              /* dwdot[OH]/dT */

    /*reaction 26: H2 + OH <=> H2O + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[7];
    k_f = prefactor_units[25] * fwd_A[25]
                * exp(fwd_beta[25] * tc[0] - activation_units[25] * fwd_Ea[25] * invT);
    dlnkfdT = fwd_beta[25] * invT + activation_units[25] * fwd_Ea[25] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[8];
    Kc = exp(-g_RT[0] + g_RT[1] + g_RT[7] - g_RT[8]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[7]) + (h_RT[0] + h_RT[8]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[1] -= q; /* H2 */
    wdot[7] -= q; /* OH */
    wdot[8] += q; /* H2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[8];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[1] -= dqdci;                /* dwdot[H2]/d[H] */
    J[7] -= dqdci;                /* dwdot[OH]/d[H] */
    J[8] += dqdci;                /* dwdot[H2O]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[7];
    J[40] += dqdci;               /* dwdot[H]/d[H2] */
    J[41] -= dqdci;               /* dwdot[H2]/d[H2] */
    J[47] -= dqdci;               /* dwdot[OH]/d[H2] */
    J[48] += dqdci;               /* dwdot[H2O]/d[H2] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[1];
    J[280] += dqdci;              /* dwdot[H]/d[OH] */
    J[281] -= dqdci;              /* dwdot[H2]/d[OH] */
    J[287] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[288] += dqdci;              /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[0];
    J[320] += dqdci;              /* dwdot[H]/d[H2O] */
    J[321] -= dqdci;              /* dwdot[H2]/d[H2O] */
    J[327] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[328] += dqdci;              /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[1560] += dqdT;              /* dwdot[H]/dT */
    J[1561] -= dqdT;              /* dwdot[H2]/dT */
    J[1567] -= dqdT;              /* dwdot[OH]/dT */
    J[1568] += dqdT;              /* dwdot[H2O]/dT */

    /*reaction 27: O + H2O <=> OH + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[8];
    k_f = prefactor_units[26] * fwd_A[26]
                * exp(fwd_beta[26] * tc[0] - activation_units[26] * fwd_Ea[26] * invT);
    dlnkfdT = fwd_beta[26] * invT + activation_units[26] * fwd_Ea[26] * invT2;
    /* reverse */
    phi_r = pow(sc[7], 2.000000);
    Kc = exp(g_RT[5] - g_RT[7] - g_RT[7] + g_RT[8]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[8]) + (2.000000*h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] -= q; /* O */
    wdot[7] += 2 * q; /* OH */
    wdot[8] -= q; /* H2O */
    /* d()/d[O] */
    dqdci =  + k_f*sc[8];
    J[205] -= dqdci;              /* dwdot[O]/d[O] */
    J[207] += 2 * dqdci;          /* dwdot[OH]/d[O] */
    J[208] -= dqdci;              /* dwdot[H2O]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*2.000000*sc[7];
    J[285] -= dqdci;              /* dwdot[O]/d[OH] */
    J[287] += 2 * dqdci;          /* dwdot[OH]/d[OH] */
    J[288] -= dqdci;              /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[5];
    J[325] -= dqdci;              /* dwdot[O]/d[H2O] */
    J[327] += 2 * dqdci;          /* dwdot[OH]/d[H2O] */
    J[328] -= dqdci;              /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[1565] -= dqdT;              /* dwdot[O]/dT */
    J[1567] += 2 * dqdT;          /* dwdot[OH]/dT */
    J[1568] -= dqdT;              /* dwdot[H2O]/dT */

    /*reaction 28: HO2 + H <=> H2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[20];
    k_f = prefactor_units[27] * fwd_A[27]
                * exp(fwd_beta[27] * tc[0] - activation_units[27] * fwd_Ea[27] * invT);
    dlnkfdT = fwd_beta[27] * invT + activation_units[27] * fwd_Ea[27] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[19];
    Kc = exp(g_RT[0] - g_RT[1] - g_RT[19] + g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[20]) + (h_RT[1] + h_RT[19]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[1] += q; /* H2 */
    wdot[19] += q; /* O2 */
    wdot[20] -= q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[20];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[1] += dqdci;                /* dwdot[H2]/d[H] */
    J[19] += dqdci;               /* dwdot[O2]/d[H] */
    J[20] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[19];
    J[40] -= dqdci;               /* dwdot[H]/d[H2] */
    J[41] += dqdci;               /* dwdot[H2]/d[H2] */
    J[59] += dqdci;               /* dwdot[O2]/d[H2] */
    J[60] -= dqdci;               /* dwdot[HO2]/d[H2] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[1];
    J[760] -= dqdci;              /* dwdot[H]/d[O2] */
    J[761] += dqdci;              /* dwdot[H2]/d[O2] */
    J[779] += dqdci;              /* dwdot[O2]/d[O2] */
    J[780] -= dqdci;              /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[0];
    J[800] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[801] += dqdci;              /* dwdot[H2]/d[HO2] */
    J[819] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[820] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1560] -= dqdT;              /* dwdot[H]/dT */
    J[1561] += dqdT;              /* dwdot[H2]/dT */
    J[1579] += dqdT;              /* dwdot[O2]/dT */
    J[1580] -= dqdT;              /* dwdot[HO2]/dT */

    /*reaction 29: HO2 + H <=> OH + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[20];
    k_f = prefactor_units[28] * fwd_A[28]
                * exp(fwd_beta[28] * tc[0] - activation_units[28] * fwd_Ea[28] * invT);
    dlnkfdT = fwd_beta[28] * invT + activation_units[28] * fwd_Ea[28] * invT2;
    /* reverse */
    phi_r = pow(sc[7], 2.000000);
    Kc = exp(g_RT[0] - g_RT[7] - g_RT[7] + g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[20]) + (2.000000*h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[7] += 2 * q; /* OH */
    wdot[20] -= q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[20];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[7] += 2 * dqdci;            /* dwdot[OH]/d[H] */
    J[20] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*2.000000*sc[7];
    J[280] -= dqdci;              /* dwdot[H]/d[OH] */
    J[287] += 2 * dqdci;          /* dwdot[OH]/d[OH] */
    J[300] -= dqdci;              /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[0];
    J[800] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[807] += 2 * dqdci;          /* dwdot[OH]/d[HO2] */
    J[820] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1560] -= dqdT;              /* dwdot[H]/dT */
    J[1567] += 2 * dqdT;          /* dwdot[OH]/dT */
    J[1580] -= dqdT;              /* dwdot[HO2]/dT */

    /*reaction 30: HO2 + O <=> O2 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[20];
    k_f = prefactor_units[29] * fwd_A[29]
                * exp(fwd_beta[29] * tc[0] - activation_units[29] * fwd_Ea[29] * invT);
    dlnkfdT = fwd_beta[29] * invT + activation_units[29] * fwd_Ea[29] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[19];
    Kc = exp(g_RT[5] - g_RT[7] - g_RT[19] + g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[20]) + (h_RT[7] + h_RT[19]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] -= q; /* O */
    wdot[7] += q; /* OH */
    wdot[19] += q; /* O2 */
    wdot[20] -= q; /* HO2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[20];
    J[205] -= dqdci;              /* dwdot[O]/d[O] */
    J[207] += dqdci;              /* dwdot[OH]/d[O] */
    J[219] += dqdci;              /* dwdot[O2]/d[O] */
    J[220] -= dqdci;              /* dwdot[HO2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[19];
    J[285] -= dqdci;              /* dwdot[O]/d[OH] */
    J[287] += dqdci;              /* dwdot[OH]/d[OH] */
    J[299] += dqdci;              /* dwdot[O2]/d[OH] */
    J[300] -= dqdci;              /* dwdot[HO2]/d[OH] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[7];
    J[765] -= dqdci;              /* dwdot[O]/d[O2] */
    J[767] += dqdci;              /* dwdot[OH]/d[O2] */
    J[779] += dqdci;              /* dwdot[O2]/d[O2] */
    J[780] -= dqdci;              /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[5];
    J[805] -= dqdci;              /* dwdot[O]/d[HO2] */
    J[807] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[819] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[820] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1565] -= dqdT;              /* dwdot[O]/dT */
    J[1567] += dqdT;              /* dwdot[OH]/dT */
    J[1579] += dqdT;              /* dwdot[O2]/dT */
    J[1580] -= dqdT;              /* dwdot[HO2]/dT */

    /*reaction 31: HO2 + OH <=> H2O + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[7]*sc[20];
    k_f = prefactor_units[30] * fwd_A[30]
                * exp(fwd_beta[30] * tc[0] - activation_units[30] * fwd_Ea[30] * invT);
    dlnkfdT = fwd_beta[30] * invT + activation_units[30] * fwd_Ea[30] * invT2;
    /* reverse */
    phi_r = sc[8]*sc[19];
    Kc = exp(g_RT[7] - g_RT[8] - g_RT[19] + g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7] + h_RT[20]) + (h_RT[8] + h_RT[19]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] -= q; /* OH */
    wdot[8] += q; /* H2O */
    wdot[19] += q; /* O2 */
    wdot[20] -= q; /* HO2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[20];
    J[287] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[288] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[299] += dqdci;              /* dwdot[O2]/d[OH] */
    J[300] -= dqdci;              /* dwdot[HO2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[19];
    J[327] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[328] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[339] += dqdci;              /* dwdot[O2]/d[H2O] */
    J[340] -= dqdci;              /* dwdot[HO2]/d[H2O] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[8];
    J[767] -= dqdci;              /* dwdot[OH]/d[O2] */
    J[768] += dqdci;              /* dwdot[H2O]/d[O2] */
    J[779] += dqdci;              /* dwdot[O2]/d[O2] */
    J[780] -= dqdci;              /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[7];
    J[807] -= dqdci;              /* dwdot[OH]/d[HO2] */
    J[808] += dqdci;              /* dwdot[H2O]/d[HO2] */
    J[819] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[820] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1567] -= dqdT;              /* dwdot[OH]/dT */
    J[1568] += dqdT;              /* dwdot[H2O]/dT */
    J[1579] += dqdT;              /* dwdot[O2]/dT */
    J[1580] -= dqdT;              /* dwdot[HO2]/dT */

    /*reaction 32: HO2 + HO2 <=> H2O2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[20], 2.000000);
    k_f = prefactor_units[31] * fwd_A[31]
                * exp(fwd_beta[31] * tc[0] - activation_units[31] * fwd_Ea[31] * invT);
    dlnkfdT = fwd_beta[31] * invT + activation_units[31] * fwd_Ea[31] * invT2;
    /* reverse */
    phi_r = sc[19]*sc[21];
    Kc = exp(-g_RT[19] + g_RT[20] + g_RT[20] - g_RT[21]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[20]) + (h_RT[19] + h_RT[21]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[19] += q; /* O2 */
    wdot[20] -= 2 * q; /* HO2 */
    wdot[21] += q; /* H2O2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[21];
    J[779] += dqdci;              /* dwdot[O2]/d[O2] */
    J[780] += -2 * dqdci;         /* dwdot[HO2]/d[O2] */
    J[781] += dqdci;              /* dwdot[H2O2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*2.000000*sc[20];
    J[819] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[820] += -2 * dqdci;         /* dwdot[HO2]/d[HO2] */
    J[821] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[19];
    J[859] += dqdci;              /* dwdot[O2]/d[H2O2] */
    J[860] += -2 * dqdci;         /* dwdot[HO2]/d[H2O2] */
    J[861] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[1579] += dqdT;              /* dwdot[O2]/dT */
    J[1580] += -2 * dqdT;         /* dwdot[HO2]/dT */
    J[1581] += dqdT;              /* dwdot[H2O2]/dT */

    /*reaction 33: HO2 + HO2 <=> H2O2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[20], 2.000000);
    k_f = prefactor_units[32] * fwd_A[32]
                * exp(fwd_beta[32] * tc[0] - activation_units[32] * fwd_Ea[32] * invT);
    dlnkfdT = fwd_beta[32] * invT + activation_units[32] * fwd_Ea[32] * invT2;
    /* reverse */
    phi_r = sc[19]*sc[21];
    Kc = exp(-g_RT[19] + g_RT[20] + g_RT[20] - g_RT[21]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[20]) + (h_RT[19] + h_RT[21]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[19] += q; /* O2 */
    wdot[20] -= 2 * q; /* HO2 */
    wdot[21] += q; /* H2O2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[21];
    J[779] += dqdci;              /* dwdot[O2]/d[O2] */
    J[780] += -2 * dqdci;         /* dwdot[HO2]/d[O2] */
    J[781] += dqdci;              /* dwdot[H2O2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*2.000000*sc[20];
    J[819] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[820] += -2 * dqdci;         /* dwdot[HO2]/d[HO2] */
    J[821] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[19];
    J[859] += dqdci;              /* dwdot[O2]/d[H2O2] */
    J[860] += -2 * dqdci;         /* dwdot[HO2]/d[H2O2] */
    J[861] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[1579] += dqdT;              /* dwdot[O2]/dT */
    J[1580] += -2 * dqdT;         /* dwdot[HO2]/dT */
    J[1581] += dqdT;              /* dwdot[H2O2]/dT */

    /*reaction 34: H2O2 + H <=> H2O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[21];
    k_f = prefactor_units[33] * fwd_A[33]
                * exp(fwd_beta[33] * tc[0] - activation_units[33] * fwd_Ea[33] * invT);
    dlnkfdT = fwd_beta[33] * invT + activation_units[33] * fwd_Ea[33] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[8];
    Kc = exp(g_RT[0] - g_RT[7] - g_RT[8] + g_RT[21]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[21]) + (h_RT[7] + h_RT[8]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[7] += q; /* OH */
    wdot[8] += q; /* H2O */
    wdot[21] -= q; /* H2O2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[21];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[7] += dqdci;                /* dwdot[OH]/d[H] */
    J[8] += dqdci;                /* dwdot[H2O]/d[H] */
    J[21] -= dqdci;               /* dwdot[H2O2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[8];
    J[280] -= dqdci;              /* dwdot[H]/d[OH] */
    J[287] += dqdci;              /* dwdot[OH]/d[OH] */
    J[288] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[301] -= dqdci;              /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[7];
    J[320] -= dqdci;              /* dwdot[H]/d[H2O] */
    J[327] += dqdci;              /* dwdot[OH]/d[H2O] */
    J[328] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[341] -= dqdci;              /* dwdot[H2O2]/d[H2O] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[0];
    J[840] -= dqdci;              /* dwdot[H]/d[H2O2] */
    J[847] += dqdci;              /* dwdot[OH]/d[H2O2] */
    J[848] += dqdci;              /* dwdot[H2O]/d[H2O2] */
    J[861] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[1560] -= dqdT;              /* dwdot[H]/dT */
    J[1567] += dqdT;              /* dwdot[OH]/dT */
    J[1568] += dqdT;              /* dwdot[H2O]/dT */
    J[1581] -= dqdT;              /* dwdot[H2O2]/dT */

    /*reaction 35: H2O2 + H <=> HO2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[21];
    k_f = prefactor_units[34] * fwd_A[34]
                * exp(fwd_beta[34] * tc[0] - activation_units[34] * fwd_Ea[34] * invT);
    dlnkfdT = fwd_beta[34] * invT + activation_units[34] * fwd_Ea[34] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[20];
    Kc = exp(g_RT[0] - g_RT[1] - g_RT[20] + g_RT[21]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[21]) + (h_RT[1] + h_RT[20]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[1] += q; /* H2 */
    wdot[20] += q; /* HO2 */
    wdot[21] -= q; /* H2O2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[21];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[1] += dqdci;                /* dwdot[H2]/d[H] */
    J[20] += dqdci;               /* dwdot[HO2]/d[H] */
    J[21] -= dqdci;               /* dwdot[H2O2]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[20];
    J[40] -= dqdci;               /* dwdot[H]/d[H2] */
    J[41] += dqdci;               /* dwdot[H2]/d[H2] */
    J[60] += dqdci;               /* dwdot[HO2]/d[H2] */
    J[61] -= dqdci;               /* dwdot[H2O2]/d[H2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[1];
    J[800] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[801] += dqdci;              /* dwdot[H2]/d[HO2] */
    J[820] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[821] -= dqdci;              /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[0];
    J[840] -= dqdci;              /* dwdot[H]/d[H2O2] */
    J[841] += dqdci;              /* dwdot[H2]/d[H2O2] */
    J[860] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[861] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[1560] -= dqdT;              /* dwdot[H]/dT */
    J[1561] += dqdT;              /* dwdot[H2]/dT */
    J[1580] += dqdT;              /* dwdot[HO2]/dT */
    J[1581] -= dqdT;              /* dwdot[H2O2]/dT */

    /*reaction 36: H2O2 + O <=> OH + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[21];
    k_f = prefactor_units[35] * fwd_A[35]
                * exp(fwd_beta[35] * tc[0] - activation_units[35] * fwd_Ea[35] * invT);
    dlnkfdT = fwd_beta[35] * invT + activation_units[35] * fwd_Ea[35] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[20];
    Kc = exp(g_RT[5] - g_RT[7] - g_RT[20] + g_RT[21]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[21]) + (h_RT[7] + h_RT[20]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] -= q; /* O */
    wdot[7] += q; /* OH */
    wdot[20] += q; /* HO2 */
    wdot[21] -= q; /* H2O2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[21];
    J[205] -= dqdci;              /* dwdot[O]/d[O] */
    J[207] += dqdci;              /* dwdot[OH]/d[O] */
    J[220] += dqdci;              /* dwdot[HO2]/d[O] */
    J[221] -= dqdci;              /* dwdot[H2O2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[20];
    J[285] -= dqdci;              /* dwdot[O]/d[OH] */
    J[287] += dqdci;              /* dwdot[OH]/d[OH] */
    J[300] += dqdci;              /* dwdot[HO2]/d[OH] */
    J[301] -= dqdci;              /* dwdot[H2O2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[7];
    J[805] -= dqdci;              /* dwdot[O]/d[HO2] */
    J[807] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[820] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[821] -= dqdci;              /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[5];
    J[845] -= dqdci;              /* dwdot[O]/d[H2O2] */
    J[847] += dqdci;              /* dwdot[OH]/d[H2O2] */
    J[860] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[861] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[1565] -= dqdT;              /* dwdot[O]/dT */
    J[1567] += dqdT;              /* dwdot[OH]/dT */
    J[1580] += dqdT;              /* dwdot[HO2]/dT */
    J[1581] -= dqdT;              /* dwdot[H2O2]/dT */

    /*reaction 37: H2O2 + OH <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[7]*sc[21];
    k_f = prefactor_units[36] * fwd_A[36]
                * exp(fwd_beta[36] * tc[0] - activation_units[36] * fwd_Ea[36] * invT);
    dlnkfdT = fwd_beta[36] * invT + activation_units[36] * fwd_Ea[36] * invT2;
    /* reverse */
    phi_r = sc[8]*sc[20];
    Kc = exp(g_RT[7] - g_RT[8] - g_RT[20] + g_RT[21]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7] + h_RT[21]) + (h_RT[8] + h_RT[20]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] -= q; /* OH */
    wdot[8] += q; /* H2O */
    wdot[20] += q; /* HO2 */
    wdot[21] -= q; /* H2O2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[21];
    J[287] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[288] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[300] += dqdci;              /* dwdot[HO2]/d[OH] */
    J[301] -= dqdci;              /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[20];
    J[327] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[328] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[340] += dqdci;              /* dwdot[HO2]/d[H2O] */
    J[341] -= dqdci;              /* dwdot[H2O2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[8];
    J[807] -= dqdci;              /* dwdot[OH]/d[HO2] */
    J[808] += dqdci;              /* dwdot[H2O]/d[HO2] */
    J[820] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[821] -= dqdci;              /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[7];
    J[847] -= dqdci;              /* dwdot[OH]/d[H2O2] */
    J[848] += dqdci;              /* dwdot[H2O]/d[H2O2] */
    J[860] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[861] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[1567] -= dqdT;              /* dwdot[OH]/dT */
    J[1568] += dqdT;              /* dwdot[H2O]/dT */
    J[1580] += dqdT;              /* dwdot[HO2]/dT */
    J[1581] -= dqdT;              /* dwdot[H2O2]/dT */

    /*reaction 38: H2O2 + OH <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[7]*sc[21];
    k_f = prefactor_units[37] * fwd_A[37]
                * exp(fwd_beta[37] * tc[0] - activation_units[37] * fwd_Ea[37] * invT);
    dlnkfdT = fwd_beta[37] * invT + activation_units[37] * fwd_Ea[37] * invT2;
    /* reverse */
    phi_r = sc[8]*sc[20];
    Kc = exp(g_RT[7] - g_RT[8] - g_RT[20] + g_RT[21]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7] + h_RT[21]) + (h_RT[8] + h_RT[20]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] -= q; /* OH */
    wdot[8] += q; /* H2O */
    wdot[20] += q; /* HO2 */
    wdot[21] -= q; /* H2O2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[21];
    J[287] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[288] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[300] += dqdci;              /* dwdot[HO2]/d[OH] */
    J[301] -= dqdci;              /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[20];
    J[327] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[328] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[340] += dqdci;              /* dwdot[HO2]/d[H2O] */
    J[341] -= dqdci;              /* dwdot[H2O2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[8];
    J[807] -= dqdci;              /* dwdot[OH]/d[HO2] */
    J[808] += dqdci;              /* dwdot[H2O]/d[HO2] */
    J[820] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[821] -= dqdci;              /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[7];
    J[847] -= dqdci;              /* dwdot[OH]/d[H2O2] */
    J[848] += dqdci;              /* dwdot[H2O]/d[H2O2] */
    J[860] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[861] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[1567] -= dqdT;              /* dwdot[OH]/dT */
    J[1568] += dqdT;              /* dwdot[H2O]/dT */
    J[1580] += dqdT;              /* dwdot[HO2]/dT */
    J[1581] -= dqdT;              /* dwdot[H2O2]/dT */

    /*reaction 39: CO + O2 <=> CO2 + O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[19];
    k_f = prefactor_units[38] * fwd_A[38]
                * exp(fwd_beta[38] * tc[0] - activation_units[38] * fwd_Ea[38] * invT);
    dlnkfdT = fwd_beta[38] * invT + activation_units[38] * fwd_Ea[38] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[22];
    Kc = exp(-g_RT[5] + g_RT[11] + g_RT[19] - g_RT[22]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[11] + h_RT[19]) + (h_RT[5] + h_RT[22]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] += q; /* O */
    wdot[11] -= q; /* CO */
    wdot[19] -= q; /* O2 */
    wdot[22] += q; /* CO2 */
    /* d()/d[O] */
    dqdci =  - k_r*sc[22];
    J[205] += dqdci;              /* dwdot[O]/d[O] */
    J[211] -= dqdci;              /* dwdot[CO]/d[O] */
    J[219] -= dqdci;              /* dwdot[O2]/d[O] */
    J[222] += dqdci;              /* dwdot[CO2]/d[O] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[19];
    J[445] += dqdci;              /* dwdot[O]/d[CO] */
    J[451] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[459] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[462] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[11];
    J[765] += dqdci;              /* dwdot[O]/d[O2] */
    J[771] -= dqdci;              /* dwdot[CO]/d[O2] */
    J[779] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[782] += dqdci;              /* dwdot[CO2]/d[O2] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[5];
    J[885] += dqdci;              /* dwdot[O]/d[CO2] */
    J[891] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[899] -= dqdci;              /* dwdot[O2]/d[CO2] */
    J[902] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[1565] += dqdT;              /* dwdot[O]/dT */
    J[1571] -= dqdT;              /* dwdot[CO]/dT */
    J[1579] -= dqdT;              /* dwdot[O2]/dT */
    J[1582] += dqdT;              /* dwdot[CO2]/dT */

    /*reaction 40: CO + HO2 <=> CO2 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[20];
    k_f = prefactor_units[39] * fwd_A[39]
                * exp(fwd_beta[39] * tc[0] - activation_units[39] * fwd_Ea[39] * invT);
    dlnkfdT = fwd_beta[39] * invT + activation_units[39] * fwd_Ea[39] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[22];
    Kc = exp(-g_RT[7] + g_RT[11] + g_RT[20] - g_RT[22]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[11] + h_RT[20]) + (h_RT[7] + h_RT[22]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] += q; /* OH */
    wdot[11] -= q; /* CO */
    wdot[20] -= q; /* HO2 */
    wdot[22] += q; /* CO2 */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[22];
    J[287] += dqdci;              /* dwdot[OH]/d[OH] */
    J[291] -= dqdci;              /* dwdot[CO]/d[OH] */
    J[300] -= dqdci;              /* dwdot[HO2]/d[OH] */
    J[302] += dqdci;              /* dwdot[CO2]/d[OH] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[20];
    J[447] += dqdci;              /* dwdot[OH]/d[CO] */
    J[451] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[460] -= dqdci;              /* dwdot[HO2]/d[CO] */
    J[462] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[11];
    J[807] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[811] -= dqdci;              /* dwdot[CO]/d[HO2] */
    J[820] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[822] += dqdci;              /* dwdot[CO2]/d[HO2] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[7];
    J[887] += dqdci;              /* dwdot[OH]/d[CO2] */
    J[891] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[900] -= dqdci;              /* dwdot[HO2]/d[CO2] */
    J[902] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[1567] += dqdT;              /* dwdot[OH]/dT */
    J[1571] -= dqdT;              /* dwdot[CO]/dT */
    J[1580] -= dqdT;              /* dwdot[HO2]/dT */
    J[1582] += dqdT;              /* dwdot[CO2]/dT */

    /*reaction 41: CO + OH <=> CO2 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[7]*sc[11];
    k_f = prefactor_units[40] * fwd_A[40]
                * exp(fwd_beta[40] * tc[0] - activation_units[40] * fwd_Ea[40] * invT);
    dlnkfdT = fwd_beta[40] * invT + activation_units[40] * fwd_Ea[40] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[22];
    Kc = exp(-g_RT[0] + g_RT[7] + g_RT[11] - g_RT[22]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7] + h_RT[11]) + (h_RT[0] + h_RT[22]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[7] -= q; /* OH */
    wdot[11] -= q; /* CO */
    wdot[22] += q; /* CO2 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[22];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[7] -= dqdci;                /* dwdot[OH]/d[H] */
    J[11] -= dqdci;               /* dwdot[CO]/d[H] */
    J[22] += dqdci;               /* dwdot[CO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[11];
    J[280] += dqdci;              /* dwdot[H]/d[OH] */
    J[287] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[291] -= dqdci;              /* dwdot[CO]/d[OH] */
    J[302] += dqdci;              /* dwdot[CO2]/d[OH] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[7];
    J[440] += dqdci;              /* dwdot[H]/d[CO] */
    J[447] -= dqdci;              /* dwdot[OH]/d[CO] */
    J[451] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[462] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[0];
    J[880] += dqdci;              /* dwdot[H]/d[CO2] */
    J[887] -= dqdci;              /* dwdot[OH]/d[CO2] */
    J[891] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[902] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[1560] += dqdT;              /* dwdot[H]/dT */
    J[1567] -= dqdT;              /* dwdot[OH]/dT */
    J[1571] -= dqdT;              /* dwdot[CO]/dT */
    J[1582] += dqdT;              /* dwdot[CO2]/dT */

    /*reaction 42: HCO + O2 <=> CO + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[13]*sc[19];
    k_f = prefactor_units[41] * fwd_A[41]
                * exp(fwd_beta[41] * tc[0] - activation_units[41] * fwd_Ea[41] * invT);
    dlnkfdT = fwd_beta[41] * invT + activation_units[41] * fwd_Ea[41] * invT2;
    /* reverse */
    phi_r = sc[11]*sc[20];
    Kc = exp(-g_RT[11] + g_RT[13] + g_RT[19] - g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[13] + h_RT[19]) + (h_RT[11] + h_RT[20]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[11] += q; /* CO */
    wdot[13] -= q; /* HCO */
    wdot[19] -= q; /* O2 */
    wdot[20] += q; /* HO2 */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[20];
    J[451] += dqdci;              /* dwdot[CO]/d[CO] */
    J[453] -= dqdci;              /* dwdot[HCO]/d[CO] */
    J[459] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[460] += dqdci;              /* dwdot[HO2]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[19];
    J[531] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[533] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    J[539] -= dqdci;              /* dwdot[O2]/d[HCO] */
    J[540] += dqdci;              /* dwdot[HO2]/d[HCO] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[13];
    J[771] += dqdci;              /* dwdot[CO]/d[O2] */
    J[773] -= dqdci;              /* dwdot[HCO]/d[O2] */
    J[779] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[780] += dqdci;              /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[11];
    J[811] += dqdci;              /* dwdot[CO]/d[HO2] */
    J[813] -= dqdci;              /* dwdot[HCO]/d[HO2] */
    J[819] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[820] += dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1571] += dqdT;              /* dwdot[CO]/dT */
    J[1573] -= dqdT;              /* dwdot[HCO]/dT */
    J[1579] -= dqdT;              /* dwdot[O2]/dT */
    J[1580] += dqdT;              /* dwdot[HO2]/dT */

    /*reaction 43: HCO + H <=> CO + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[13];
    k_f = prefactor_units[42] * fwd_A[42]
                * exp(fwd_beta[42] * tc[0] - activation_units[42] * fwd_Ea[42] * invT);
    dlnkfdT = fwd_beta[42] * invT + activation_units[42] * fwd_Ea[42] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[11];
    Kc = exp(g_RT[0] - g_RT[1] - g_RT[11] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[13]) + (h_RT[1] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[1] += q; /* H2 */
    wdot[11] += q; /* CO */
    wdot[13] -= q; /* HCO */
    /* d()/d[H] */
    dqdci =  + k_f*sc[13];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[1] += dqdci;                /* dwdot[H2]/d[H] */
    J[11] += dqdci;               /* dwdot[CO]/d[H] */
    J[13] -= dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[11];
    J[40] -= dqdci;               /* dwdot[H]/d[H2] */
    J[41] += dqdci;               /* dwdot[H2]/d[H2] */
    J[51] += dqdci;               /* dwdot[CO]/d[H2] */
    J[53] -= dqdci;               /* dwdot[HCO]/d[H2] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[1];
    J[440] -= dqdci;              /* dwdot[H]/d[CO] */
    J[441] += dqdci;              /* dwdot[H2]/d[CO] */
    J[451] += dqdci;              /* dwdot[CO]/d[CO] */
    J[453] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[0];
    J[520] -= dqdci;              /* dwdot[H]/d[HCO] */
    J[521] += dqdci;              /* dwdot[H2]/d[HCO] */
    J[531] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[533] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[1560] -= dqdT;              /* dwdot[H]/dT */
    J[1561] += dqdT;              /* dwdot[H2]/dT */
    J[1571] += dqdT;              /* dwdot[CO]/dT */
    J[1573] -= dqdT;              /* dwdot[HCO]/dT */

    /*reaction 44: HCO + O <=> CO + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[13];
    k_f = prefactor_units[43] * fwd_A[43]
                * exp(fwd_beta[43] * tc[0] - activation_units[43] * fwd_Ea[43] * invT);
    dlnkfdT = fwd_beta[43] * invT + activation_units[43] * fwd_Ea[43] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[11];
    Kc = exp(g_RT[5] - g_RT[7] - g_RT[11] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[13]) + (h_RT[7] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] -= q; /* O */
    wdot[7] += q; /* OH */
    wdot[11] += q; /* CO */
    wdot[13] -= q; /* HCO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[13];
    J[205] -= dqdci;              /* dwdot[O]/d[O] */
    J[207] += dqdci;              /* dwdot[OH]/d[O] */
    J[211] += dqdci;              /* dwdot[CO]/d[O] */
    J[213] -= dqdci;              /* dwdot[HCO]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[11];
    J[285] -= dqdci;              /* dwdot[O]/d[OH] */
    J[287] += dqdci;              /* dwdot[OH]/d[OH] */
    J[291] += dqdci;              /* dwdot[CO]/d[OH] */
    J[293] -= dqdci;              /* dwdot[HCO]/d[OH] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[7];
    J[445] -= dqdci;              /* dwdot[O]/d[CO] */
    J[447] += dqdci;              /* dwdot[OH]/d[CO] */
    J[451] += dqdci;              /* dwdot[CO]/d[CO] */
    J[453] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[5];
    J[525] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[527] += dqdci;              /* dwdot[OH]/d[HCO] */
    J[531] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[533] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[1565] -= dqdT;              /* dwdot[O]/dT */
    J[1567] += dqdT;              /* dwdot[OH]/dT */
    J[1571] += dqdT;              /* dwdot[CO]/dT */
    J[1573] -= dqdT;              /* dwdot[HCO]/dT */

    /*reaction 45: HCO + OH <=> CO + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[7]*sc[13];
    k_f = prefactor_units[44] * fwd_A[44]
                * exp(fwd_beta[44] * tc[0] - activation_units[44] * fwd_Ea[44] * invT);
    dlnkfdT = fwd_beta[44] * invT + activation_units[44] * fwd_Ea[44] * invT2;
    /* reverse */
    phi_r = sc[8]*sc[11];
    Kc = exp(g_RT[7] - g_RT[8] - g_RT[11] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7] + h_RT[13]) + (h_RT[8] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] -= q; /* OH */
    wdot[8] += q; /* H2O */
    wdot[11] += q; /* CO */
    wdot[13] -= q; /* HCO */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[13];
    J[287] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[288] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[291] += dqdci;              /* dwdot[CO]/d[OH] */
    J[293] -= dqdci;              /* dwdot[HCO]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[11];
    J[327] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[328] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[331] += dqdci;              /* dwdot[CO]/d[H2O] */
    J[333] -= dqdci;              /* dwdot[HCO]/d[H2O] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[8];
    J[447] -= dqdci;              /* dwdot[OH]/d[CO] */
    J[448] += dqdci;              /* dwdot[H2O]/d[CO] */
    J[451] += dqdci;              /* dwdot[CO]/d[CO] */
    J[453] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[7];
    J[527] -= dqdci;              /* dwdot[OH]/d[HCO] */
    J[528] += dqdci;              /* dwdot[H2O]/d[HCO] */
    J[531] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[533] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[1567] -= dqdT;              /* dwdot[OH]/dT */
    J[1568] += dqdT;              /* dwdot[H2O]/dT */
    J[1571] += dqdT;              /* dwdot[CO]/dT */
    J[1573] -= dqdT;              /* dwdot[HCO]/dT */

    /*reaction 46: HCO + O <=> CO2 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[13];
    k_f = prefactor_units[45] * fwd_A[45]
                * exp(fwd_beta[45] * tc[0] - activation_units[45] * fwd_Ea[45] * invT);
    dlnkfdT = fwd_beta[45] * invT + activation_units[45] * fwd_Ea[45] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[22];
    Kc = exp(-g_RT[0] + g_RT[5] + g_RT[13] - g_RT[22]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[13]) + (h_RT[0] + h_RT[22]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[5] -= q; /* O */
    wdot[13] -= q; /* HCO */
    wdot[22] += q; /* CO2 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[22];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[5] -= dqdci;                /* dwdot[O]/d[H] */
    J[13] -= dqdci;               /* dwdot[HCO]/d[H] */
    J[22] += dqdci;               /* dwdot[CO2]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[13];
    J[200] += dqdci;              /* dwdot[H]/d[O] */
    J[205] -= dqdci;              /* dwdot[O]/d[O] */
    J[213] -= dqdci;              /* dwdot[HCO]/d[O] */
    J[222] += dqdci;              /* dwdot[CO2]/d[O] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[5];
    J[520] += dqdci;              /* dwdot[H]/d[HCO] */
    J[525] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[533] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    J[542] += dqdci;              /* dwdot[CO2]/d[HCO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[0];
    J[880] += dqdci;              /* dwdot[H]/d[CO2] */
    J[885] -= dqdci;              /* dwdot[O]/d[CO2] */
    J[893] -= dqdci;              /* dwdot[HCO]/d[CO2] */
    J[902] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[1560] += dqdT;              /* dwdot[H]/dT */
    J[1565] -= dqdT;              /* dwdot[O]/dT */
    J[1573] -= dqdT;              /* dwdot[HCO]/dT */
    J[1582] += dqdT;              /* dwdot[CO2]/dT */

    /*reaction 47: HCO + HO2 <=> CO2 + OH + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[13]*sc[20];
    k_f = prefactor_units[46] * fwd_A[46]
                * exp(fwd_beta[46] * tc[0] - activation_units[46] * fwd_Ea[46] * invT);
    dlnkfdT = fwd_beta[46] * invT + activation_units[46] * fwd_Ea[46] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[7]*sc[22];
    Kc = refC * exp(-g_RT[0] - g_RT[7] + g_RT[13] + g_RT[20] - g_RT[22]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[13] + h_RT[20]) + (h_RT[0] + h_RT[7] + h_RT[22]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[7] += q; /* OH */
    wdot[13] -= q; /* HCO */
    wdot[20] -= q; /* HO2 */
    wdot[22] += q; /* CO2 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[7]*sc[22];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[7] += dqdci;                /* dwdot[OH]/d[H] */
    J[13] -= dqdci;               /* dwdot[HCO]/d[H] */
    J[20] -= dqdci;               /* dwdot[HO2]/d[H] */
    J[22] += dqdci;               /* dwdot[CO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[0]*sc[22];
    J[280] += dqdci;              /* dwdot[H]/d[OH] */
    J[287] += dqdci;              /* dwdot[OH]/d[OH] */
    J[293] -= dqdci;              /* dwdot[HCO]/d[OH] */
    J[300] -= dqdci;              /* dwdot[HO2]/d[OH] */
    J[302] += dqdci;              /* dwdot[CO2]/d[OH] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[20];
    J[520] += dqdci;              /* dwdot[H]/d[HCO] */
    J[527] += dqdci;              /* dwdot[OH]/d[HCO] */
    J[533] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    J[540] -= dqdci;              /* dwdot[HO2]/d[HCO] */
    J[542] += dqdci;              /* dwdot[CO2]/d[HCO] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[13];
    J[800] += dqdci;              /* dwdot[H]/d[HO2] */
    J[807] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[813] -= dqdci;              /* dwdot[HCO]/d[HO2] */
    J[820] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[822] += dqdci;              /* dwdot[CO2]/d[HO2] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[0]*sc[7];
    J[880] += dqdci;              /* dwdot[H]/d[CO2] */
    J[887] += dqdci;              /* dwdot[OH]/d[CO2] */
    J[893] -= dqdci;              /* dwdot[HCO]/d[CO2] */
    J[900] -= dqdci;              /* dwdot[HO2]/d[CO2] */
    J[902] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[1560] += dqdT;              /* dwdot[H]/dT */
    J[1567] += dqdT;              /* dwdot[OH]/dT */
    J[1573] -= dqdT;              /* dwdot[HCO]/dT */
    J[1580] -= dqdT;              /* dwdot[HO2]/dT */
    J[1582] += dqdT;              /* dwdot[CO2]/dT */

    /*reaction 48: HCO + HCO <=> H2 + CO + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[13], 2.000000);
    k_f = prefactor_units[47] * fwd_A[47]
                * exp(fwd_beta[47] * tc[0] - activation_units[47] * fwd_Ea[47] * invT);
    dlnkfdT = fwd_beta[47] * invT + activation_units[47] * fwd_Ea[47] * invT2;
    /* reverse */
    phi_r = sc[1]*pow(sc[11], 2.000000);
    Kc = refC * exp(-g_RT[1] - g_RT[11] - g_RT[11] + g_RT[13] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[13]) + (h_RT[1] + 2.000000*h_RT[11]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H2 */
    wdot[11] += 2 * q; /* CO */
    wdot[13] -= 2 * q; /* HCO */
    /* d()/d[H2] */
    dqdci =  - k_r*pow(sc[11], 2.000000);
    J[41] += dqdci;               /* dwdot[H2]/d[H2] */
    J[51] += 2 * dqdci;           /* dwdot[CO]/d[H2] */
    J[53] += -2 * dqdci;          /* dwdot[HCO]/d[H2] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[1]*2.000000*sc[11];
    J[441] += dqdci;              /* dwdot[H2]/d[CO] */
    J[451] += 2 * dqdci;          /* dwdot[CO]/d[CO] */
    J[453] += -2 * dqdci;         /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*2.000000*sc[13];
    J[521] += dqdci;              /* dwdot[H2]/d[HCO] */
    J[531] += 2 * dqdci;          /* dwdot[CO]/d[HCO] */
    J[533] += -2 * dqdci;         /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[1561] += dqdT;              /* dwdot[H2]/dT */
    J[1571] += 2 * dqdT;          /* dwdot[CO]/dT */
    J[1573] += -2 * dqdT;         /* dwdot[HCO]/dT */

    /*reaction 49: HCO + CH3 <=> CO + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[13];
    k_f = prefactor_units[48] * fwd_A[48]
                * exp(fwd_beta[48] * tc[0] - activation_units[48] * fwd_Ea[48] * invT);
    dlnkfdT = fwd_beta[48] * invT + activation_units[48] * fwd_Ea[48] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[11];
    Kc = exp(g_RT[4] - g_RT[6] - g_RT[11] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[13]) + (h_RT[6] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* CH3 */
    wdot[6] += q; /* CH4 */
    wdot[11] += q; /* CO */
    wdot[13] -= q; /* HCO */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[13];
    J[164] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[166] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[171] += dqdci;              /* dwdot[CO]/d[CH3] */
    J[173] -= dqdci;              /* dwdot[HCO]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[11];
    J[244] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[246] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[251] += dqdci;              /* dwdot[CO]/d[CH4] */
    J[253] -= dqdci;              /* dwdot[HCO]/d[CH4] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[6];
    J[444] -= dqdci;              /* dwdot[CH3]/d[CO] */
    J[446] += dqdci;              /* dwdot[CH4]/d[CO] */
    J[451] += dqdci;              /* dwdot[CO]/d[CO] */
    J[453] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[4];
    J[524] -= dqdci;              /* dwdot[CH3]/d[HCO] */
    J[526] += dqdci;              /* dwdot[CH4]/d[HCO] */
    J[531] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[533] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[1564] -= dqdT;              /* dwdot[CH3]/dT */
    J[1566] += dqdT;              /* dwdot[CH4]/dT */
    J[1571] += dqdT;              /* dwdot[CO]/dT */
    J[1573] -= dqdT;              /* dwdot[HCO]/dT */

    /*reaction 50: HCO + HCO <=> CH2O + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[13], 2.000000);
    k_f = prefactor_units[49] * fwd_A[49]
                * exp(fwd_beta[49] * tc[0] - activation_units[49] * fwd_Ea[49] * invT);
    dlnkfdT = fwd_beta[49] * invT + activation_units[49] * fwd_Ea[49] * invT2;
    /* reverse */
    phi_r = sc[11]*sc[15];
    Kc = exp(-g_RT[11] + g_RT[13] + g_RT[13] - g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[13]) + (h_RT[11] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[11] += q; /* CO */
    wdot[13] -= 2 * q; /* HCO */
    wdot[15] += q; /* CH2O */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[15];
    J[451] += dqdci;              /* dwdot[CO]/d[CO] */
    J[453] += -2 * dqdci;         /* dwdot[HCO]/d[CO] */
    J[455] += dqdci;              /* dwdot[CH2O]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*2.000000*sc[13];
    J[531] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[533] += -2 * dqdci;         /* dwdot[HCO]/d[HCO] */
    J[535] += dqdci;              /* dwdot[CH2O]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[11];
    J[611] += dqdci;              /* dwdot[CO]/d[CH2O] */
    J[613] += -2 * dqdci;         /* dwdot[HCO]/d[CH2O] */
    J[615] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[1571] += dqdT;              /* dwdot[CO]/dT */
    J[1573] += -2 * dqdT;         /* dwdot[HCO]/dT */
    J[1575] += dqdT;              /* dwdot[CH2O]/dT */

    /*reaction 51: CH2O + H <=> HCO + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[15];
    k_f = prefactor_units[50] * fwd_A[50]
                * exp(fwd_beta[50] * tc[0] - activation_units[50] * fwd_Ea[50] * invT);
    dlnkfdT = fwd_beta[50] * invT + activation_units[50] * fwd_Ea[50] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[13];
    Kc = exp(g_RT[0] - g_RT[1] - g_RT[13] + g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[15]) + (h_RT[1] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[1] += q; /* H2 */
    wdot[13] += q; /* HCO */
    wdot[15] -= q; /* CH2O */
    /* d()/d[H] */
    dqdci =  + k_f*sc[15];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[1] += dqdci;                /* dwdot[H2]/d[H] */
    J[13] += dqdci;               /* dwdot[HCO]/d[H] */
    J[15] -= dqdci;               /* dwdot[CH2O]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[13];
    J[40] -= dqdci;               /* dwdot[H]/d[H2] */
    J[41] += dqdci;               /* dwdot[H2]/d[H2] */
    J[53] += dqdci;               /* dwdot[HCO]/d[H2] */
    J[55] -= dqdci;               /* dwdot[CH2O]/d[H2] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[1];
    J[520] -= dqdci;              /* dwdot[H]/d[HCO] */
    J[521] += dqdci;              /* dwdot[H2]/d[HCO] */
    J[533] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[535] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[0];
    J[600] -= dqdci;              /* dwdot[H]/d[CH2O] */
    J[601] += dqdci;              /* dwdot[H2]/d[CH2O] */
    J[613] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[615] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[1560] -= dqdT;              /* dwdot[H]/dT */
    J[1561] += dqdT;              /* dwdot[H2]/dT */
    J[1573] += dqdT;              /* dwdot[HCO]/dT */
    J[1575] -= dqdT;              /* dwdot[CH2O]/dT */

    /*reaction 52: CH2O + O <=> HCO + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[15];
    k_f = prefactor_units[51] * fwd_A[51]
                * exp(fwd_beta[51] * tc[0] - activation_units[51] * fwd_Ea[51] * invT);
    dlnkfdT = fwd_beta[51] * invT + activation_units[51] * fwd_Ea[51] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[13];
    Kc = exp(g_RT[5] - g_RT[7] - g_RT[13] + g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[15]) + (h_RT[7] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] -= q; /* O */
    wdot[7] += q; /* OH */
    wdot[13] += q; /* HCO */
    wdot[15] -= q; /* CH2O */
    /* d()/d[O] */
    dqdci =  + k_f*sc[15];
    J[205] -= dqdci;              /* dwdot[O]/d[O] */
    J[207] += dqdci;              /* dwdot[OH]/d[O] */
    J[213] += dqdci;              /* dwdot[HCO]/d[O] */
    J[215] -= dqdci;              /* dwdot[CH2O]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[13];
    J[285] -= dqdci;              /* dwdot[O]/d[OH] */
    J[287] += dqdci;              /* dwdot[OH]/d[OH] */
    J[293] += dqdci;              /* dwdot[HCO]/d[OH] */
    J[295] -= dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[7];
    J[525] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[527] += dqdci;              /* dwdot[OH]/d[HCO] */
    J[533] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[535] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[5];
    J[605] -= dqdci;              /* dwdot[O]/d[CH2O] */
    J[607] += dqdci;              /* dwdot[OH]/d[CH2O] */
    J[613] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[615] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[1565] -= dqdT;              /* dwdot[O]/dT */
    J[1567] += dqdT;              /* dwdot[OH]/dT */
    J[1573] += dqdT;              /* dwdot[HCO]/dT */
    J[1575] -= dqdT;              /* dwdot[CH2O]/dT */

    /*reaction 53: CH2O + OH <=> HCO + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[7]*sc[15];
    k_f = prefactor_units[52] * fwd_A[52]
                * exp(fwd_beta[52] * tc[0] - activation_units[52] * fwd_Ea[52] * invT);
    dlnkfdT = fwd_beta[52] * invT + activation_units[52] * fwd_Ea[52] * invT2;
    /* reverse */
    phi_r = sc[8]*sc[13];
    Kc = exp(g_RT[7] - g_RT[8] - g_RT[13] + g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7] + h_RT[15]) + (h_RT[8] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] -= q; /* OH */
    wdot[8] += q; /* H2O */
    wdot[13] += q; /* HCO */
    wdot[15] -= q; /* CH2O */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[15];
    J[287] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[288] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[293] += dqdci;              /* dwdot[HCO]/d[OH] */
    J[295] -= dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[13];
    J[327] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[328] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[333] += dqdci;              /* dwdot[HCO]/d[H2O] */
    J[335] -= dqdci;              /* dwdot[CH2O]/d[H2O] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[8];
    J[527] -= dqdci;              /* dwdot[OH]/d[HCO] */
    J[528] += dqdci;              /* dwdot[H2O]/d[HCO] */
    J[533] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[535] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[7];
    J[607] -= dqdci;              /* dwdot[OH]/d[CH2O] */
    J[608] += dqdci;              /* dwdot[H2O]/d[CH2O] */
    J[613] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[615] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[1567] -= dqdT;              /* dwdot[OH]/dT */
    J[1568] += dqdT;              /* dwdot[H2O]/dT */
    J[1573] += dqdT;              /* dwdot[HCO]/dT */
    J[1575] -= dqdT;              /* dwdot[CH2O]/dT */

    /*reaction 54: CH2O + O2 <=> HCO + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[15]*sc[19];
    k_f = prefactor_units[53] * fwd_A[53]
                * exp(fwd_beta[53] * tc[0] - activation_units[53] * fwd_Ea[53] * invT);
    dlnkfdT = fwd_beta[53] * invT + activation_units[53] * fwd_Ea[53] * invT2;
    /* reverse */
    phi_r = sc[13]*sc[20];
    Kc = exp(-g_RT[13] + g_RT[15] + g_RT[19] - g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[15] + h_RT[19]) + (h_RT[13] + h_RT[20]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[13] += q; /* HCO */
    wdot[15] -= q; /* CH2O */
    wdot[19] -= q; /* O2 */
    wdot[20] += q; /* HO2 */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[20];
    J[533] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[535] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    J[539] -= dqdci;              /* dwdot[O2]/d[HCO] */
    J[540] += dqdci;              /* dwdot[HO2]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[19];
    J[613] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[615] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[619] -= dqdci;              /* dwdot[O2]/d[CH2O] */
    J[620] += dqdci;              /* dwdot[HO2]/d[CH2O] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[15];
    J[773] += dqdci;              /* dwdot[HCO]/d[O2] */
    J[775] -= dqdci;              /* dwdot[CH2O]/d[O2] */
    J[779] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[780] += dqdci;              /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[13];
    J[813] += dqdci;              /* dwdot[HCO]/d[HO2] */
    J[815] -= dqdci;              /* dwdot[CH2O]/d[HO2] */
    J[819] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[820] += dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1573] += dqdT;              /* dwdot[HCO]/dT */
    J[1575] -= dqdT;              /* dwdot[CH2O]/dT */
    J[1579] -= dqdT;              /* dwdot[O2]/dT */
    J[1580] += dqdT;              /* dwdot[HO2]/dT */

    /*reaction 55: CH2O + HO2 <=> HCO + H2O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[15]*sc[20];
    k_f = prefactor_units[54] * fwd_A[54]
                * exp(fwd_beta[54] * tc[0] - activation_units[54] * fwd_Ea[54] * invT);
    dlnkfdT = fwd_beta[54] * invT + activation_units[54] * fwd_Ea[54] * invT2;
    /* reverse */
    phi_r = sc[13]*sc[21];
    Kc = exp(-g_RT[13] + g_RT[15] + g_RT[20] - g_RT[21]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[15] + h_RT[20]) + (h_RT[13] + h_RT[21]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[13] += q; /* HCO */
    wdot[15] -= q; /* CH2O */
    wdot[20] -= q; /* HO2 */
    wdot[21] += q; /* H2O2 */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[21];
    J[533] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[535] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    J[540] -= dqdci;              /* dwdot[HO2]/d[HCO] */
    J[541] += dqdci;              /* dwdot[H2O2]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[20];
    J[613] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[615] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[620] -= dqdci;              /* dwdot[HO2]/d[CH2O] */
    J[621] += dqdci;              /* dwdot[H2O2]/d[CH2O] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[15];
    J[813] += dqdci;              /* dwdot[HCO]/d[HO2] */
    J[815] -= dqdci;              /* dwdot[CH2O]/d[HO2] */
    J[820] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[821] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[13];
    J[853] += dqdci;              /* dwdot[HCO]/d[H2O2] */
    J[855] -= dqdci;              /* dwdot[CH2O]/d[H2O2] */
    J[860] -= dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[861] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[1573] += dqdT;              /* dwdot[HCO]/dT */
    J[1575] -= dqdT;              /* dwdot[CH2O]/dT */
    J[1580] -= dqdT;              /* dwdot[HO2]/dT */
    J[1581] += dqdT;              /* dwdot[H2O2]/dT */

    /*reaction 56: CH2O + CH3 <=> HCO + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[15];
    k_f = prefactor_units[55] * fwd_A[55]
                * exp(fwd_beta[55] * tc[0] - activation_units[55] * fwd_Ea[55] * invT);
    dlnkfdT = fwd_beta[55] * invT + activation_units[55] * fwd_Ea[55] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[13];
    Kc = exp(g_RT[4] - g_RT[6] - g_RT[13] + g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[15]) + (h_RT[6] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* CH3 */
    wdot[6] += q; /* CH4 */
    wdot[13] += q; /* HCO */
    wdot[15] -= q; /* CH2O */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[15];
    J[164] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[166] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[173] += dqdci;              /* dwdot[HCO]/d[CH3] */
    J[175] -= dqdci;              /* dwdot[CH2O]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[13];
    J[244] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[246] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[253] += dqdci;              /* dwdot[HCO]/d[CH4] */
    J[255] -= dqdci;              /* dwdot[CH2O]/d[CH4] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[6];
    J[524] -= dqdci;              /* dwdot[CH3]/d[HCO] */
    J[526] += dqdci;              /* dwdot[CH4]/d[HCO] */
    J[533] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[535] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[4];
    J[604] -= dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[606] += dqdci;              /* dwdot[CH4]/d[CH2O] */
    J[613] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[615] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[1564] -= dqdT;              /* dwdot[CH3]/dT */
    J[1566] += dqdT;              /* dwdot[CH4]/dT */
    J[1573] += dqdT;              /* dwdot[HCO]/dT */
    J[1575] -= dqdT;              /* dwdot[CH2O]/dT */

    /*reaction 57: CH3 + O <=> CH2O + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[5];
    k_f = prefactor_units[56] * fwd_A[56]
                * exp(fwd_beta[56] * tc[0] - activation_units[56] * fwd_Ea[56] * invT);
    dlnkfdT = fwd_beta[56] * invT + activation_units[56] * fwd_Ea[56] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[15];
    Kc = exp(-g_RT[0] + g_RT[4] + g_RT[5] - g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[5]) + (h_RT[0] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[4] -= q; /* CH3 */
    wdot[5] -= q; /* O */
    wdot[15] += q; /* CH2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[15];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[4] -= dqdci;                /* dwdot[CH3]/d[H] */
    J[5] -= dqdci;                /* dwdot[O]/d[H] */
    J[15] += dqdci;               /* dwdot[CH2O]/d[H] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[5];
    J[160] += dqdci;              /* dwdot[H]/d[CH3] */
    J[164] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[165] -= dqdci;              /* dwdot[O]/d[CH3] */
    J[175] += dqdci;              /* dwdot[CH2O]/d[CH3] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[4];
    J[200] += dqdci;              /* dwdot[H]/d[O] */
    J[204] -= dqdci;              /* dwdot[CH3]/d[O] */
    J[205] -= dqdci;              /* dwdot[O]/d[O] */
    J[215] += dqdci;              /* dwdot[CH2O]/d[O] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[0];
    J[600] += dqdci;              /* dwdot[H]/d[CH2O] */
    J[604] -= dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[605] -= dqdci;              /* dwdot[O]/d[CH2O] */
    J[615] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[1560] += dqdT;              /* dwdot[H]/dT */
    J[1564] -= dqdT;              /* dwdot[CH3]/dT */
    J[1565] -= dqdT;              /* dwdot[O]/dT */
    J[1575] += dqdT;              /* dwdot[CH2O]/dT */

    /*reaction 58: CH3 + O2 <=> CH3O + O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[19];
    k_f = prefactor_units[57] * fwd_A[57]
                * exp(fwd_beta[57] * tc[0] - activation_units[57] * fwd_Ea[57] * invT);
    dlnkfdT = fwd_beta[57] * invT + activation_units[57] * fwd_Ea[57] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[18];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[18] + g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[19]) + (h_RT[5] + h_RT[18]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* CH3 */
    wdot[5] += q; /* O */
    wdot[18] += q; /* CH3O */
    wdot[19] -= q; /* O2 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[19];
    J[164] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[165] += dqdci;              /* dwdot[O]/d[CH3] */
    J[178] += dqdci;              /* dwdot[CH3O]/d[CH3] */
    J[179] -= dqdci;              /* dwdot[O2]/d[CH3] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[18];
    J[204] -= dqdci;              /* dwdot[CH3]/d[O] */
    J[205] += dqdci;              /* dwdot[O]/d[O] */
    J[218] += dqdci;              /* dwdot[CH3O]/d[O] */
    J[219] -= dqdci;              /* dwdot[O2]/d[O] */
    /* d()/d[CH3O] */
    dqdci =  - k_r*sc[5];
    J[724] -= dqdci;              /* dwdot[CH3]/d[CH3O] */
    J[725] += dqdci;              /* dwdot[O]/d[CH3O] */
    J[738] += dqdci;              /* dwdot[CH3O]/d[CH3O] */
    J[739] -= dqdci;              /* dwdot[O2]/d[CH3O] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[4];
    J[764] -= dqdci;              /* dwdot[CH3]/d[O2] */
    J[765] += dqdci;              /* dwdot[O]/d[O2] */
    J[778] += dqdci;              /* dwdot[CH3O]/d[O2] */
    J[779] -= dqdci;              /* dwdot[O2]/d[O2] */
    /* d()/dT */
    J[1564] -= dqdT;              /* dwdot[CH3]/dT */
    J[1565] += dqdT;              /* dwdot[O]/dT */
    J[1578] += dqdT;              /* dwdot[CH3O]/dT */
    J[1579] -= dqdT;              /* dwdot[O2]/dT */

    /*reaction 59: CH3 + O2 <=> CH2O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[19];
    k_f = prefactor_units[58] * fwd_A[58]
                * exp(fwd_beta[58] * tc[0] - activation_units[58] * fwd_Ea[58] * invT);
    dlnkfdT = fwd_beta[58] * invT + activation_units[58] * fwd_Ea[58] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[15];
    Kc = exp(g_RT[4] - g_RT[7] - g_RT[15] + g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[19]) + (h_RT[7] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* CH3 */
    wdot[7] += q; /* OH */
    wdot[15] += q; /* CH2O */
    wdot[19] -= q; /* O2 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[19];
    J[164] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[167] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[175] += dqdci;              /* dwdot[CH2O]/d[CH3] */
    J[179] -= dqdci;              /* dwdot[O2]/d[CH3] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[15];
    J[284] -= dqdci;              /* dwdot[CH3]/d[OH] */
    J[287] += dqdci;              /* dwdot[OH]/d[OH] */
    J[295] += dqdci;              /* dwdot[CH2O]/d[OH] */
    J[299] -= dqdci;              /* dwdot[O2]/d[OH] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[7];
    J[604] -= dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[607] += dqdci;              /* dwdot[OH]/d[CH2O] */
    J[615] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[619] -= dqdci;              /* dwdot[O2]/d[CH2O] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[4];
    J[764] -= dqdci;              /* dwdot[CH3]/d[O2] */
    J[767] += dqdci;              /* dwdot[OH]/d[O2] */
    J[775] += dqdci;              /* dwdot[CH2O]/d[O2] */
    J[779] -= dqdci;              /* dwdot[O2]/d[O2] */
    /* d()/dT */
    J[1564] -= dqdT;              /* dwdot[CH3]/dT */
    J[1567] += dqdT;              /* dwdot[OH]/dT */
    J[1575] += dqdT;              /* dwdot[CH2O]/dT */
    J[1579] -= dqdT;              /* dwdot[O2]/dT */

    /*reaction 60: CH3 + HO2 <=> CH3O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[20];
    k_f = prefactor_units[59] * fwd_A[59]
                * exp(fwd_beta[59] * tc[0] - activation_units[59] * fwd_Ea[59] * invT);
    dlnkfdT = fwd_beta[59] * invT + activation_units[59] * fwd_Ea[59] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[18];
    Kc = exp(g_RT[4] - g_RT[7] - g_RT[18] + g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[20]) + (h_RT[7] + h_RT[18]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* CH3 */
    wdot[7] += q; /* OH */
    wdot[18] += q; /* CH3O */
    wdot[20] -= q; /* HO2 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[20];
    J[164] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[167] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[178] += dqdci;              /* dwdot[CH3O]/d[CH3] */
    J[180] -= dqdci;              /* dwdot[HO2]/d[CH3] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[18];
    J[284] -= dqdci;              /* dwdot[CH3]/d[OH] */
    J[287] += dqdci;              /* dwdot[OH]/d[OH] */
    J[298] += dqdci;              /* dwdot[CH3O]/d[OH] */
    J[300] -= dqdci;              /* dwdot[HO2]/d[OH] */
    /* d()/d[CH3O] */
    dqdci =  - k_r*sc[7];
    J[724] -= dqdci;              /* dwdot[CH3]/d[CH3O] */
    J[727] += dqdci;              /* dwdot[OH]/d[CH3O] */
    J[738] += dqdci;              /* dwdot[CH3O]/d[CH3O] */
    J[740] -= dqdci;              /* dwdot[HO2]/d[CH3O] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[4];
    J[804] -= dqdci;              /* dwdot[CH3]/d[HO2] */
    J[807] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[818] += dqdci;              /* dwdot[CH3O]/d[HO2] */
    J[820] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1564] -= dqdT;              /* dwdot[CH3]/dT */
    J[1567] += dqdT;              /* dwdot[OH]/dT */
    J[1578] += dqdT;              /* dwdot[CH3O]/dT */
    J[1580] -= dqdT;              /* dwdot[HO2]/dT */

    /*reaction 61: CH4 + H <=> CH3 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[6];
    k_f = prefactor_units[60] * fwd_A[60]
                * exp(fwd_beta[60] * tc[0] - activation_units[60] * fwd_Ea[60] * invT);
    dlnkfdT = fwd_beta[60] * invT + activation_units[60] * fwd_Ea[60] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[4];
    Kc = exp(g_RT[0] - g_RT[1] - g_RT[4] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[6]) + (h_RT[1] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[1] += q; /* H2 */
    wdot[4] += q; /* CH3 */
    wdot[6] -= q; /* CH4 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[1] += dqdci;                /* dwdot[H2]/d[H] */
    J[4] += dqdci;                /* dwdot[CH3]/d[H] */
    J[6] -= dqdci;                /* dwdot[CH4]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[4];
    J[40] -= dqdci;               /* dwdot[H]/d[H2] */
    J[41] += dqdci;               /* dwdot[H2]/d[H2] */
    J[44] += dqdci;               /* dwdot[CH3]/d[H2] */
    J[46] -= dqdci;               /* dwdot[CH4]/d[H2] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[1];
    J[160] -= dqdci;              /* dwdot[H]/d[CH3] */
    J[161] += dqdci;              /* dwdot[H2]/d[CH3] */
    J[164] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[166] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[0];
    J[240] -= dqdci;              /* dwdot[H]/d[CH4] */
    J[241] += dqdci;              /* dwdot[H2]/d[CH4] */
    J[244] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[246] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[1560] -= dqdT;              /* dwdot[H]/dT */
    J[1561] += dqdT;              /* dwdot[H2]/dT */
    J[1564] += dqdT;              /* dwdot[CH3]/dT */
    J[1566] -= dqdT;              /* dwdot[CH4]/dT */

    /*reaction 62: CH4 + O <=> CH3 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[6];
    k_f = prefactor_units[61] * fwd_A[61]
                * exp(fwd_beta[61] * tc[0] - activation_units[61] * fwd_Ea[61] * invT);
    dlnkfdT = fwd_beta[61] * invT + activation_units[61] * fwd_Ea[61] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[7];
    Kc = exp(-g_RT[4] + g_RT[5] + g_RT[6] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[6]) + (h_RT[4] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] += q; /* CH3 */
    wdot[5] -= q; /* O */
    wdot[6] -= q; /* CH4 */
    wdot[7] += q; /* OH */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[7];
    J[164] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[165] -= dqdci;              /* dwdot[O]/d[CH3] */
    J[166] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    J[167] += dqdci;              /* dwdot[OH]/d[CH3] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[6];
    J[204] += dqdci;              /* dwdot[CH3]/d[O] */
    J[205] -= dqdci;              /* dwdot[O]/d[O] */
    J[206] -= dqdci;              /* dwdot[CH4]/d[O] */
    J[207] += dqdci;              /* dwdot[OH]/d[O] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[5];
    J[244] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[245] -= dqdci;              /* dwdot[O]/d[CH4] */
    J[246] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    J[247] += dqdci;              /* dwdot[OH]/d[CH4] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[4];
    J[284] += dqdci;              /* dwdot[CH3]/d[OH] */
    J[285] -= dqdci;              /* dwdot[O]/d[OH] */
    J[286] -= dqdci;              /* dwdot[CH4]/d[OH] */
    J[287] += dqdci;              /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[1564] += dqdT;              /* dwdot[CH3]/dT */
    J[1565] -= dqdT;              /* dwdot[O]/dT */
    J[1566] -= dqdT;              /* dwdot[CH4]/dT */
    J[1567] += dqdT;              /* dwdot[OH]/dT */

    /*reaction 63: CH4 + OH <=> CH3 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[7];
    k_f = prefactor_units[62] * fwd_A[62]
                * exp(fwd_beta[62] * tc[0] - activation_units[62] * fwd_Ea[62] * invT);
    dlnkfdT = fwd_beta[62] * invT + activation_units[62] * fwd_Ea[62] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[8];
    Kc = exp(-g_RT[4] + g_RT[6] + g_RT[7] - g_RT[8]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[7]) + (h_RT[4] + h_RT[8]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] += q; /* CH3 */
    wdot[6] -= q; /* CH4 */
    wdot[7] -= q; /* OH */
    wdot[8] += q; /* H2O */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[8];
    J[164] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[166] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    J[167] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[168] += dqdci;              /* dwdot[H2O]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[7];
    J[244] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[246] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    J[247] -= dqdci;              /* dwdot[OH]/d[CH4] */
    J[248] += dqdci;              /* dwdot[H2O]/d[CH4] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[6];
    J[284] += dqdci;              /* dwdot[CH3]/d[OH] */
    J[286] -= dqdci;              /* dwdot[CH4]/d[OH] */
    J[287] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[288] += dqdci;              /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[4];
    J[324] += dqdci;              /* dwdot[CH3]/d[H2O] */
    J[326] -= dqdci;              /* dwdot[CH4]/d[H2O] */
    J[327] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[328] += dqdci;              /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[1564] += dqdT;              /* dwdot[CH3]/dT */
    J[1566] -= dqdT;              /* dwdot[CH4]/dT */
    J[1567] -= dqdT;              /* dwdot[OH]/dT */
    J[1568] += dqdT;              /* dwdot[H2O]/dT */

    /*reaction 64: CH3 + HO2 <=> CH4 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[20];
    k_f = prefactor_units[63] * fwd_A[63]
                * exp(fwd_beta[63] * tc[0] - activation_units[63] * fwd_Ea[63] * invT);
    dlnkfdT = fwd_beta[63] * invT + activation_units[63] * fwd_Ea[63] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[19];
    Kc = exp(g_RT[4] - g_RT[6] - g_RT[19] + g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[20]) + (h_RT[6] + h_RT[19]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* CH3 */
    wdot[6] += q; /* CH4 */
    wdot[19] += q; /* O2 */
    wdot[20] -= q; /* HO2 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[20];
    J[164] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[166] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[179] += dqdci;              /* dwdot[O2]/d[CH3] */
    J[180] -= dqdci;              /* dwdot[HO2]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[19];
    J[244] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[246] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[259] += dqdci;              /* dwdot[O2]/d[CH4] */
    J[260] -= dqdci;              /* dwdot[HO2]/d[CH4] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[6];
    J[764] -= dqdci;              /* dwdot[CH3]/d[O2] */
    J[766] += dqdci;              /* dwdot[CH4]/d[O2] */
    J[779] += dqdci;              /* dwdot[O2]/d[O2] */
    J[780] -= dqdci;              /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[4];
    J[804] -= dqdci;              /* dwdot[CH3]/d[HO2] */
    J[806] += dqdci;              /* dwdot[CH4]/d[HO2] */
    J[819] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[820] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1564] -= dqdT;              /* dwdot[CH3]/dT */
    J[1566] += dqdT;              /* dwdot[CH4]/dT */
    J[1579] += dqdT;              /* dwdot[O2]/dT */
    J[1580] -= dqdT;              /* dwdot[HO2]/dT */

    /*reaction 65: CH4 + HO2 <=> CH3 + H2O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[20];
    k_f = prefactor_units[64] * fwd_A[64]
                * exp(fwd_beta[64] * tc[0] - activation_units[64] * fwd_Ea[64] * invT);
    dlnkfdT = fwd_beta[64] * invT + activation_units[64] * fwd_Ea[64] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[21];
    Kc = exp(-g_RT[4] + g_RT[6] + g_RT[20] - g_RT[21]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[20]) + (h_RT[4] + h_RT[21]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] += q; /* CH3 */
    wdot[6] -= q; /* CH4 */
    wdot[20] -= q; /* HO2 */
    wdot[21] += q; /* H2O2 */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[21];
    J[164] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[166] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    J[180] -= dqdci;              /* dwdot[HO2]/d[CH3] */
    J[181] += dqdci;              /* dwdot[H2O2]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[20];
    J[244] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[246] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    J[260] -= dqdci;              /* dwdot[HO2]/d[CH4] */
    J[261] += dqdci;              /* dwdot[H2O2]/d[CH4] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[6];
    J[804] += dqdci;              /* dwdot[CH3]/d[HO2] */
    J[806] -= dqdci;              /* dwdot[CH4]/d[HO2] */
    J[820] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[821] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[4];
    J[844] += dqdci;              /* dwdot[CH3]/d[H2O2] */
    J[846] -= dqdci;              /* dwdot[CH4]/d[H2O2] */
    J[860] -= dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[861] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[1564] += dqdT;              /* dwdot[CH3]/dT */
    J[1566] -= dqdT;              /* dwdot[CH4]/dT */
    J[1580] -= dqdT;              /* dwdot[HO2]/dT */
    J[1581] += dqdT;              /* dwdot[H2O2]/dT */

    /*reaction 66: CH2OH + H <=> CH2O + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[17];
    k_f = prefactor_units[65] * fwd_A[65]
                * exp(fwd_beta[65] * tc[0] - activation_units[65] * fwd_Ea[65] * invT);
    dlnkfdT = fwd_beta[65] * invT + activation_units[65] * fwd_Ea[65] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[15];
    Kc = exp(g_RT[0] - g_RT[1] - g_RT[15] + g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[17]) + (h_RT[1] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[1] += q; /* H2 */
    wdot[15] += q; /* CH2O */
    wdot[17] -= q; /* CH2OH */
    /* d()/d[H] */
    dqdci =  + k_f*sc[17];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[1] += dqdci;                /* dwdot[H2]/d[H] */
    J[15] += dqdci;               /* dwdot[CH2O]/d[H] */
    J[17] -= dqdci;               /* dwdot[CH2OH]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[15];
    J[40] -= dqdci;               /* dwdot[H]/d[H2] */
    J[41] += dqdci;               /* dwdot[H2]/d[H2] */
    J[55] += dqdci;               /* dwdot[CH2O]/d[H2] */
    J[57] -= dqdci;               /* dwdot[CH2OH]/d[H2] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[1];
    J[600] -= dqdci;              /* dwdot[H]/d[CH2O] */
    J[601] += dqdci;              /* dwdot[H2]/d[CH2O] */
    J[615] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[617] -= dqdci;              /* dwdot[CH2OH]/d[CH2O] */
    /* d()/d[CH2OH] */
    dqdci =  + k_f*sc[0];
    J[680] -= dqdci;              /* dwdot[H]/d[CH2OH] */
    J[681] += dqdci;              /* dwdot[H2]/d[CH2OH] */
    J[695] += dqdci;              /* dwdot[CH2O]/d[CH2OH] */
    J[697] -= dqdci;              /* dwdot[CH2OH]/d[CH2OH] */
    /* d()/dT */
    J[1560] -= dqdT;              /* dwdot[H]/dT */
    J[1561] += dqdT;              /* dwdot[H2]/dT */
    J[1575] += dqdT;              /* dwdot[CH2O]/dT */
    J[1577] -= dqdT;              /* dwdot[CH2OH]/dT */

    /*reaction 67: CH2OH + H <=> CH3 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[17];
    k_f = prefactor_units[66] * fwd_A[66]
                * exp(fwd_beta[66] * tc[0] - activation_units[66] * fwd_Ea[66] * invT);
    dlnkfdT = fwd_beta[66] * invT + activation_units[66] * fwd_Ea[66] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[7];
    Kc = exp(g_RT[0] - g_RT[4] - g_RT[7] + g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[17]) + (h_RT[4] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[4] += q; /* CH3 */
    wdot[7] += q; /* OH */
    wdot[17] -= q; /* CH2OH */
    /* d()/d[H] */
    dqdci =  + k_f*sc[17];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[4] += dqdci;                /* dwdot[CH3]/d[H] */
    J[7] += dqdci;                /* dwdot[OH]/d[H] */
    J[17] -= dqdci;               /* dwdot[CH2OH]/d[H] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[7];
    J[160] -= dqdci;              /* dwdot[H]/d[CH3] */
    J[164] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[167] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[177] -= dqdci;              /* dwdot[CH2OH]/d[CH3] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[4];
    J[280] -= dqdci;              /* dwdot[H]/d[OH] */
    J[284] += dqdci;              /* dwdot[CH3]/d[OH] */
    J[287] += dqdci;              /* dwdot[OH]/d[OH] */
    J[297] -= dqdci;              /* dwdot[CH2OH]/d[OH] */
    /* d()/d[CH2OH] */
    dqdci =  + k_f*sc[0];
    J[680] -= dqdci;              /* dwdot[H]/d[CH2OH] */
    J[684] += dqdci;              /* dwdot[CH3]/d[CH2OH] */
    J[687] += dqdci;              /* dwdot[OH]/d[CH2OH] */
    J[697] -= dqdci;              /* dwdot[CH2OH]/d[CH2OH] */
    /* d()/dT */
    J[1560] -= dqdT;              /* dwdot[H]/dT */
    J[1564] += dqdT;              /* dwdot[CH3]/dT */
    J[1567] += dqdT;              /* dwdot[OH]/dT */
    J[1577] -= dqdT;              /* dwdot[CH2OH]/dT */

    /*reaction 68: CH2OH + O <=> CH2O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[17];
    k_f = prefactor_units[67] * fwd_A[67]
                * exp(fwd_beta[67] * tc[0] - activation_units[67] * fwd_Ea[67] * invT);
    dlnkfdT = fwd_beta[67] * invT + activation_units[67] * fwd_Ea[67] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[15];
    Kc = exp(g_RT[5] - g_RT[7] - g_RT[15] + g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[17]) + (h_RT[7] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] -= q; /* O */
    wdot[7] += q; /* OH */
    wdot[15] += q; /* CH2O */
    wdot[17] -= q; /* CH2OH */
    /* d()/d[O] */
    dqdci =  + k_f*sc[17];
    J[205] -= dqdci;              /* dwdot[O]/d[O] */
    J[207] += dqdci;              /* dwdot[OH]/d[O] */
    J[215] += dqdci;              /* dwdot[CH2O]/d[O] */
    J[217] -= dqdci;              /* dwdot[CH2OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[15];
    J[285] -= dqdci;              /* dwdot[O]/d[OH] */
    J[287] += dqdci;              /* dwdot[OH]/d[OH] */
    J[295] += dqdci;              /* dwdot[CH2O]/d[OH] */
    J[297] -= dqdci;              /* dwdot[CH2OH]/d[OH] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[7];
    J[605] -= dqdci;              /* dwdot[O]/d[CH2O] */
    J[607] += dqdci;              /* dwdot[OH]/d[CH2O] */
    J[615] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[617] -= dqdci;              /* dwdot[CH2OH]/d[CH2O] */
    /* d()/d[CH2OH] */
    dqdci =  + k_f*sc[5];
    J[685] -= dqdci;              /* dwdot[O]/d[CH2OH] */
    J[687] += dqdci;              /* dwdot[OH]/d[CH2OH] */
    J[695] += dqdci;              /* dwdot[CH2O]/d[CH2OH] */
    J[697] -= dqdci;              /* dwdot[CH2OH]/d[CH2OH] */
    /* d()/dT */
    J[1565] -= dqdT;              /* dwdot[O]/dT */
    J[1567] += dqdT;              /* dwdot[OH]/dT */
    J[1575] += dqdT;              /* dwdot[CH2O]/dT */
    J[1577] -= dqdT;              /* dwdot[CH2OH]/dT */

    /*reaction 69: CH2OH + OH <=> CH2O + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[7]*sc[17];
    k_f = prefactor_units[68] * fwd_A[68]
                * exp(fwd_beta[68] * tc[0] - activation_units[68] * fwd_Ea[68] * invT);
    dlnkfdT = fwd_beta[68] * invT + activation_units[68] * fwd_Ea[68] * invT2;
    /* reverse */
    phi_r = sc[8]*sc[15];
    Kc = exp(g_RT[7] - g_RT[8] - g_RT[15] + g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7] + h_RT[17]) + (h_RT[8] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] -= q; /* OH */
    wdot[8] += q; /* H2O */
    wdot[15] += q; /* CH2O */
    wdot[17] -= q; /* CH2OH */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[17];
    J[287] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[288] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[295] += dqdci;              /* dwdot[CH2O]/d[OH] */
    J[297] -= dqdci;              /* dwdot[CH2OH]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[15];
    J[327] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[328] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[335] += dqdci;              /* dwdot[CH2O]/d[H2O] */
    J[337] -= dqdci;              /* dwdot[CH2OH]/d[H2O] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[8];
    J[607] -= dqdci;              /* dwdot[OH]/d[CH2O] */
    J[608] += dqdci;              /* dwdot[H2O]/d[CH2O] */
    J[615] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[617] -= dqdci;              /* dwdot[CH2OH]/d[CH2O] */
    /* d()/d[CH2OH] */
    dqdci =  + k_f*sc[7];
    J[687] -= dqdci;              /* dwdot[OH]/d[CH2OH] */
    J[688] += dqdci;              /* dwdot[H2O]/d[CH2OH] */
    J[695] += dqdci;              /* dwdot[CH2O]/d[CH2OH] */
    J[697] -= dqdci;              /* dwdot[CH2OH]/d[CH2OH] */
    /* d()/dT */
    J[1567] -= dqdT;              /* dwdot[OH]/dT */
    J[1568] += dqdT;              /* dwdot[H2O]/dT */
    J[1575] += dqdT;              /* dwdot[CH2O]/dT */
    J[1577] -= dqdT;              /* dwdot[CH2OH]/dT */

    /*reaction 70: CH2OH + O2 <=> CH2O + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[17]*sc[19];
    k_f = prefactor_units[69] * fwd_A[69]
                * exp(fwd_beta[69] * tc[0] - activation_units[69] * fwd_Ea[69] * invT);
    dlnkfdT = fwd_beta[69] * invT + activation_units[69] * fwd_Ea[69] * invT2;
    /* reverse */
    phi_r = sc[15]*sc[20];
    Kc = exp(-g_RT[15] + g_RT[17] + g_RT[19] - g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[17] + h_RT[19]) + (h_RT[15] + h_RT[20]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[15] += q; /* CH2O */
    wdot[17] -= q; /* CH2OH */
    wdot[19] -= q; /* O2 */
    wdot[20] += q; /* HO2 */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[20];
    J[615] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[617] -= dqdci;              /* dwdot[CH2OH]/d[CH2O] */
    J[619] -= dqdci;              /* dwdot[O2]/d[CH2O] */
    J[620] += dqdci;              /* dwdot[HO2]/d[CH2O] */
    /* d()/d[CH2OH] */
    dqdci =  + k_f*sc[19];
    J[695] += dqdci;              /* dwdot[CH2O]/d[CH2OH] */
    J[697] -= dqdci;              /* dwdot[CH2OH]/d[CH2OH] */
    J[699] -= dqdci;              /* dwdot[O2]/d[CH2OH] */
    J[700] += dqdci;              /* dwdot[HO2]/d[CH2OH] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[17];
    J[775] += dqdci;              /* dwdot[CH2O]/d[O2] */
    J[777] -= dqdci;              /* dwdot[CH2OH]/d[O2] */
    J[779] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[780] += dqdci;              /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[15];
    J[815] += dqdci;              /* dwdot[CH2O]/d[HO2] */
    J[817] -= dqdci;              /* dwdot[CH2OH]/d[HO2] */
    J[819] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[820] += dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1575] += dqdT;              /* dwdot[CH2O]/dT */
    J[1577] -= dqdT;              /* dwdot[CH2OH]/dT */
    J[1579] -= dqdT;              /* dwdot[O2]/dT */
    J[1580] += dqdT;              /* dwdot[HO2]/dT */

    /*reaction 71: CH2OH + O2 <=> CH2O + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[17]*sc[19];
    k_f = prefactor_units[70] * fwd_A[70]
                * exp(fwd_beta[70] * tc[0] - activation_units[70] * fwd_Ea[70] * invT);
    dlnkfdT = fwd_beta[70] * invT + activation_units[70] * fwd_Ea[70] * invT2;
    /* reverse */
    phi_r = sc[15]*sc[20];
    Kc = exp(-g_RT[15] + g_RT[17] + g_RT[19] - g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[17] + h_RT[19]) + (h_RT[15] + h_RT[20]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[15] += q; /* CH2O */
    wdot[17] -= q; /* CH2OH */
    wdot[19] -= q; /* O2 */
    wdot[20] += q; /* HO2 */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[20];
    J[615] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[617] -= dqdci;              /* dwdot[CH2OH]/d[CH2O] */
    J[619] -= dqdci;              /* dwdot[O2]/d[CH2O] */
    J[620] += dqdci;              /* dwdot[HO2]/d[CH2O] */
    /* d()/d[CH2OH] */
    dqdci =  + k_f*sc[19];
    J[695] += dqdci;              /* dwdot[CH2O]/d[CH2OH] */
    J[697] -= dqdci;              /* dwdot[CH2OH]/d[CH2OH] */
    J[699] -= dqdci;              /* dwdot[O2]/d[CH2OH] */
    J[700] += dqdci;              /* dwdot[HO2]/d[CH2OH] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[17];
    J[775] += dqdci;              /* dwdot[CH2O]/d[O2] */
    J[777] -= dqdci;              /* dwdot[CH2OH]/d[O2] */
    J[779] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[780] += dqdci;              /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[15];
    J[815] += dqdci;              /* dwdot[CH2O]/d[HO2] */
    J[817] -= dqdci;              /* dwdot[CH2OH]/d[HO2] */
    J[819] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[820] += dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1575] += dqdT;              /* dwdot[CH2O]/dT */
    J[1577] -= dqdT;              /* dwdot[CH2OH]/dT */
    J[1579] -= dqdT;              /* dwdot[O2]/dT */
    J[1580] += dqdT;              /* dwdot[HO2]/dT */

    /*reaction 72: CH2OH + HO2 <=> CH2O + H2O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[17]*sc[20];
    k_f = prefactor_units[71] * fwd_A[71]
                * exp(fwd_beta[71] * tc[0] - activation_units[71] * fwd_Ea[71] * invT);
    dlnkfdT = fwd_beta[71] * invT + activation_units[71] * fwd_Ea[71] * invT2;
    /* reverse */
    phi_r = sc[15]*sc[21];
    Kc = exp(-g_RT[15] + g_RT[17] + g_RT[20] - g_RT[21]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[17] + h_RT[20]) + (h_RT[15] + h_RT[21]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[15] += q; /* CH2O */
    wdot[17] -= q; /* CH2OH */
    wdot[20] -= q; /* HO2 */
    wdot[21] += q; /* H2O2 */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[21];
    J[615] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[617] -= dqdci;              /* dwdot[CH2OH]/d[CH2O] */
    J[620] -= dqdci;              /* dwdot[HO2]/d[CH2O] */
    J[621] += dqdci;              /* dwdot[H2O2]/d[CH2O] */
    /* d()/d[CH2OH] */
    dqdci =  + k_f*sc[20];
    J[695] += dqdci;              /* dwdot[CH2O]/d[CH2OH] */
    J[697] -= dqdci;              /* dwdot[CH2OH]/d[CH2OH] */
    J[700] -= dqdci;              /* dwdot[HO2]/d[CH2OH] */
    J[701] += dqdci;              /* dwdot[H2O2]/d[CH2OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[17];
    J[815] += dqdci;              /* dwdot[CH2O]/d[HO2] */
    J[817] -= dqdci;              /* dwdot[CH2OH]/d[HO2] */
    J[820] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[821] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[15];
    J[855] += dqdci;              /* dwdot[CH2O]/d[H2O2] */
    J[857] -= dqdci;              /* dwdot[CH2OH]/d[H2O2] */
    J[860] -= dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[861] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[1575] += dqdT;              /* dwdot[CH2O]/dT */
    J[1577] -= dqdT;              /* dwdot[CH2OH]/dT */
    J[1580] -= dqdT;              /* dwdot[HO2]/dT */
    J[1581] += dqdT;              /* dwdot[H2O2]/dT */

    /*reaction 73: CH2OH + HCO <=> CH2O + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[13]*sc[17];
    k_f = prefactor_units[72] * fwd_A[72]
                * exp(fwd_beta[72] * tc[0] - activation_units[72] * fwd_Ea[72] * invT);
    dlnkfdT = fwd_beta[72] * invT + activation_units[72] * fwd_Ea[72] * invT2;
    /* reverse */
    phi_r = pow(sc[15], 2.000000);
    Kc = exp(g_RT[13] - g_RT[15] - g_RT[15] + g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[13] + h_RT[17]) + (2.000000*h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[13] -= q; /* HCO */
    wdot[15] += 2 * q; /* CH2O */
    wdot[17] -= q; /* CH2OH */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[17];
    J[533] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    J[535] += 2 * dqdci;          /* dwdot[CH2O]/d[HCO] */
    J[537] -= dqdci;              /* dwdot[CH2OH]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*2.000000*sc[15];
    J[613] -= dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[615] += 2 * dqdci;          /* dwdot[CH2O]/d[CH2O] */
    J[617] -= dqdci;              /* dwdot[CH2OH]/d[CH2O] */
    /* d()/d[CH2OH] */
    dqdci =  + k_f*sc[13];
    J[693] -= dqdci;              /* dwdot[HCO]/d[CH2OH] */
    J[695] += 2 * dqdci;          /* dwdot[CH2O]/d[CH2OH] */
    J[697] -= dqdci;              /* dwdot[CH2OH]/d[CH2OH] */
    /* d()/dT */
    J[1573] -= dqdT;              /* dwdot[HCO]/dT */
    J[1575] += 2 * dqdT;          /* dwdot[CH2O]/dT */
    J[1577] -= dqdT;              /* dwdot[CH2OH]/dT */

    /*reaction 74: CH3O + H <=> CH3 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[18];
    k_f = prefactor_units[73] * fwd_A[73]
                * exp(fwd_beta[73] * tc[0] - activation_units[73] * fwd_Ea[73] * invT);
    dlnkfdT = fwd_beta[73] * invT + activation_units[73] * fwd_Ea[73] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[7];
    Kc = exp(g_RT[0] - g_RT[4] - g_RT[7] + g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[18]) + (h_RT[4] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[4] += q; /* CH3 */
    wdot[7] += q; /* OH */
    wdot[18] -= q; /* CH3O */
    /* d()/d[H] */
    dqdci =  + k_f*sc[18];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[4] += dqdci;                /* dwdot[CH3]/d[H] */
    J[7] += dqdci;                /* dwdot[OH]/d[H] */
    J[18] -= dqdci;               /* dwdot[CH3O]/d[H] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[7];
    J[160] -= dqdci;              /* dwdot[H]/d[CH3] */
    J[164] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[167] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[178] -= dqdci;              /* dwdot[CH3O]/d[CH3] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[4];
    J[280] -= dqdci;              /* dwdot[H]/d[OH] */
    J[284] += dqdci;              /* dwdot[CH3]/d[OH] */
    J[287] += dqdci;              /* dwdot[OH]/d[OH] */
    J[298] -= dqdci;              /* dwdot[CH3O]/d[OH] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[0];
    J[720] -= dqdci;              /* dwdot[H]/d[CH3O] */
    J[724] += dqdci;              /* dwdot[CH3]/d[CH3O] */
    J[727] += dqdci;              /* dwdot[OH]/d[CH3O] */
    J[738] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[1560] -= dqdT;              /* dwdot[H]/dT */
    J[1564] += dqdT;              /* dwdot[CH3]/dT */
    J[1567] += dqdT;              /* dwdot[OH]/dT */
    J[1578] -= dqdT;              /* dwdot[CH3O]/dT */

    /*reaction 75: CH3O + O <=> CH2O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[18];
    k_f = prefactor_units[74] * fwd_A[74]
                * exp(fwd_beta[74] * tc[0] - activation_units[74] * fwd_Ea[74] * invT);
    dlnkfdT = fwd_beta[74] * invT + activation_units[74] * fwd_Ea[74] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[15];
    Kc = exp(g_RT[5] - g_RT[7] - g_RT[15] + g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[18]) + (h_RT[7] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] -= q; /* O */
    wdot[7] += q; /* OH */
    wdot[15] += q; /* CH2O */
    wdot[18] -= q; /* CH3O */
    /* d()/d[O] */
    dqdci =  + k_f*sc[18];
    J[205] -= dqdci;              /* dwdot[O]/d[O] */
    J[207] += dqdci;              /* dwdot[OH]/d[O] */
    J[215] += dqdci;              /* dwdot[CH2O]/d[O] */
    J[218] -= dqdci;              /* dwdot[CH3O]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[15];
    J[285] -= dqdci;              /* dwdot[O]/d[OH] */
    J[287] += dqdci;              /* dwdot[OH]/d[OH] */
    J[295] += dqdci;              /* dwdot[CH2O]/d[OH] */
    J[298] -= dqdci;              /* dwdot[CH3O]/d[OH] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[7];
    J[605] -= dqdci;              /* dwdot[O]/d[CH2O] */
    J[607] += dqdci;              /* dwdot[OH]/d[CH2O] */
    J[615] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[618] -= dqdci;              /* dwdot[CH3O]/d[CH2O] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[5];
    J[725] -= dqdci;              /* dwdot[O]/d[CH3O] */
    J[727] += dqdci;              /* dwdot[OH]/d[CH3O] */
    J[735] += dqdci;              /* dwdot[CH2O]/d[CH3O] */
    J[738] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[1565] -= dqdT;              /* dwdot[O]/dT */
    J[1567] += dqdT;              /* dwdot[OH]/dT */
    J[1575] += dqdT;              /* dwdot[CH2O]/dT */
    J[1578] -= dqdT;              /* dwdot[CH3O]/dT */

    /*reaction 76: CH3O + OH <=> CH2O + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[7]*sc[18];
    k_f = prefactor_units[75] * fwd_A[75]
                * exp(fwd_beta[75] * tc[0] - activation_units[75] * fwd_Ea[75] * invT);
    dlnkfdT = fwd_beta[75] * invT + activation_units[75] * fwd_Ea[75] * invT2;
    /* reverse */
    phi_r = sc[8]*sc[15];
    Kc = exp(g_RT[7] - g_RT[8] - g_RT[15] + g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7] + h_RT[18]) + (h_RT[8] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] -= q; /* OH */
    wdot[8] += q; /* H2O */
    wdot[15] += q; /* CH2O */
    wdot[18] -= q; /* CH3O */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[18];
    J[287] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[288] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[295] += dqdci;              /* dwdot[CH2O]/d[OH] */
    J[298] -= dqdci;              /* dwdot[CH3O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[15];
    J[327] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[328] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[335] += dqdci;              /* dwdot[CH2O]/d[H2O] */
    J[338] -= dqdci;              /* dwdot[CH3O]/d[H2O] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[8];
    J[607] -= dqdci;              /* dwdot[OH]/d[CH2O] */
    J[608] += dqdci;              /* dwdot[H2O]/d[CH2O] */
    J[615] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[618] -= dqdci;              /* dwdot[CH3O]/d[CH2O] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[7];
    J[727] -= dqdci;              /* dwdot[OH]/d[CH3O] */
    J[728] += dqdci;              /* dwdot[H2O]/d[CH3O] */
    J[735] += dqdci;              /* dwdot[CH2O]/d[CH3O] */
    J[738] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[1567] -= dqdT;              /* dwdot[OH]/dT */
    J[1568] += dqdT;              /* dwdot[H2O]/dT */
    J[1575] += dqdT;              /* dwdot[CH2O]/dT */
    J[1578] -= dqdT;              /* dwdot[CH3O]/dT */

    /*reaction 77: CH3O + O2 <=> CH2O + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[18]*sc[19];
    k_f = prefactor_units[76] * fwd_A[76]
                * exp(fwd_beta[76] * tc[0] - activation_units[76] * fwd_Ea[76] * invT);
    dlnkfdT = fwd_beta[76] * invT + activation_units[76] * fwd_Ea[76] * invT2;
    /* reverse */
    phi_r = sc[15]*sc[20];
    Kc = exp(-g_RT[15] + g_RT[18] + g_RT[19] - g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[18] + h_RT[19]) + (h_RT[15] + h_RT[20]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[15] += q; /* CH2O */
    wdot[18] -= q; /* CH3O */
    wdot[19] -= q; /* O2 */
    wdot[20] += q; /* HO2 */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[20];
    J[615] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[618] -= dqdci;              /* dwdot[CH3O]/d[CH2O] */
    J[619] -= dqdci;              /* dwdot[O2]/d[CH2O] */
    J[620] += dqdci;              /* dwdot[HO2]/d[CH2O] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[19];
    J[735] += dqdci;              /* dwdot[CH2O]/d[CH3O] */
    J[738] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    J[739] -= dqdci;              /* dwdot[O2]/d[CH3O] */
    J[740] += dqdci;              /* dwdot[HO2]/d[CH3O] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[18];
    J[775] += dqdci;              /* dwdot[CH2O]/d[O2] */
    J[778] -= dqdci;              /* dwdot[CH3O]/d[O2] */
    J[779] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[780] += dqdci;              /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[15];
    J[815] += dqdci;              /* dwdot[CH2O]/d[HO2] */
    J[818] -= dqdci;              /* dwdot[CH3O]/d[HO2] */
    J[819] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[820] += dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1575] += dqdT;              /* dwdot[CH2O]/dT */
    J[1578] -= dqdT;              /* dwdot[CH3O]/dT */
    J[1579] -= dqdT;              /* dwdot[O2]/dT */
    J[1580] += dqdT;              /* dwdot[HO2]/dT */

    /*reaction 78: CH3O + O2 <=> CH2O + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[18]*sc[19];
    k_f = prefactor_units[77] * fwd_A[77]
                * exp(fwd_beta[77] * tc[0] - activation_units[77] * fwd_Ea[77] * invT);
    dlnkfdT = fwd_beta[77] * invT + activation_units[77] * fwd_Ea[77] * invT2;
    /* reverse */
    phi_r = sc[15]*sc[20];
    Kc = exp(-g_RT[15] + g_RT[18] + g_RT[19] - g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[18] + h_RT[19]) + (h_RT[15] + h_RT[20]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[15] += q; /* CH2O */
    wdot[18] -= q; /* CH3O */
    wdot[19] -= q; /* O2 */
    wdot[20] += q; /* HO2 */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[20];
    J[615] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[618] -= dqdci;              /* dwdot[CH3O]/d[CH2O] */
    J[619] -= dqdci;              /* dwdot[O2]/d[CH2O] */
    J[620] += dqdci;              /* dwdot[HO2]/d[CH2O] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[19];
    J[735] += dqdci;              /* dwdot[CH2O]/d[CH3O] */
    J[738] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    J[739] -= dqdci;              /* dwdot[O2]/d[CH3O] */
    J[740] += dqdci;              /* dwdot[HO2]/d[CH3O] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[18];
    J[775] += dqdci;              /* dwdot[CH2O]/d[O2] */
    J[778] -= dqdci;              /* dwdot[CH3O]/d[O2] */
    J[779] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[780] += dqdci;              /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[15];
    J[815] += dqdci;              /* dwdot[CH2O]/d[HO2] */
    J[818] -= dqdci;              /* dwdot[CH3O]/d[HO2] */
    J[819] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[820] += dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1575] += dqdT;              /* dwdot[CH2O]/dT */
    J[1578] -= dqdT;              /* dwdot[CH3O]/dT */
    J[1579] -= dqdT;              /* dwdot[O2]/dT */
    J[1580] += dqdT;              /* dwdot[HO2]/dT */

    /*reaction 79: CH3O + HO2 <=> CH2O + H2O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[18]*sc[20];
    k_f = prefactor_units[78] * fwd_A[78]
                * exp(fwd_beta[78] * tc[0] - activation_units[78] * fwd_Ea[78] * invT);
    dlnkfdT = fwd_beta[78] * invT + activation_units[78] * fwd_Ea[78] * invT2;
    /* reverse */
    phi_r = sc[15]*sc[21];
    Kc = exp(-g_RT[15] + g_RT[18] + g_RT[20] - g_RT[21]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[18] + h_RT[20]) + (h_RT[15] + h_RT[21]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[15] += q; /* CH2O */
    wdot[18] -= q; /* CH3O */
    wdot[20] -= q; /* HO2 */
    wdot[21] += q; /* H2O2 */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[21];
    J[615] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[618] -= dqdci;              /* dwdot[CH3O]/d[CH2O] */
    J[620] -= dqdci;              /* dwdot[HO2]/d[CH2O] */
    J[621] += dqdci;              /* dwdot[H2O2]/d[CH2O] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[20];
    J[735] += dqdci;              /* dwdot[CH2O]/d[CH3O] */
    J[738] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    J[740] -= dqdci;              /* dwdot[HO2]/d[CH3O] */
    J[741] += dqdci;              /* dwdot[H2O2]/d[CH3O] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[18];
    J[815] += dqdci;              /* dwdot[CH2O]/d[HO2] */
    J[818] -= dqdci;              /* dwdot[CH3O]/d[HO2] */
    J[820] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[821] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[15];
    J[855] += dqdci;              /* dwdot[CH2O]/d[H2O2] */
    J[858] -= dqdci;              /* dwdot[CH3O]/d[H2O2] */
    J[860] -= dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[861] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[1575] += dqdT;              /* dwdot[CH2O]/dT */
    J[1578] -= dqdT;              /* dwdot[CH3O]/dT */
    J[1580] -= dqdT;              /* dwdot[HO2]/dT */
    J[1581] += dqdT;              /* dwdot[H2O2]/dT */

    /*reaction 80: CH3O + CO <=> CH3 + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[18];
    k_f = prefactor_units[79] * fwd_A[79]
                * exp(fwd_beta[79] * tc[0] - activation_units[79] * fwd_Ea[79] * invT);
    dlnkfdT = fwd_beta[79] * invT + activation_units[79] * fwd_Ea[79] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[22];
    Kc = exp(-g_RT[4] + g_RT[11] + g_RT[18] - g_RT[22]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[11] + h_RT[18]) + (h_RT[4] + h_RT[22]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] += q; /* CH3 */
    wdot[11] -= q; /* CO */
    wdot[18] -= q; /* CH3O */
    wdot[22] += q; /* CO2 */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[22];
    J[164] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[171] -= dqdci;              /* dwdot[CO]/d[CH3] */
    J[178] -= dqdci;              /* dwdot[CH3O]/d[CH3] */
    J[182] += dqdci;              /* dwdot[CO2]/d[CH3] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[18];
    J[444] += dqdci;              /* dwdot[CH3]/d[CO] */
    J[451] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[458] -= dqdci;              /* dwdot[CH3O]/d[CO] */
    J[462] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[11];
    J[724] += dqdci;              /* dwdot[CH3]/d[CH3O] */
    J[731] -= dqdci;              /* dwdot[CO]/d[CH3O] */
    J[738] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    J[742] += dqdci;              /* dwdot[CO2]/d[CH3O] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[4];
    J[884] += dqdci;              /* dwdot[CH3]/d[CO2] */
    J[891] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[898] -= dqdci;              /* dwdot[CH3O]/d[CO2] */
    J[902] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[1564] += dqdT;              /* dwdot[CH3]/dT */
    J[1571] -= dqdT;              /* dwdot[CO]/dT */
    J[1578] -= dqdT;              /* dwdot[CH3O]/dT */
    J[1582] += dqdT;              /* dwdot[CO2]/dT */

    /*reaction 81: CH3 + CH3 <=> H + C2H5 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[4], 2.000000);
    k_f = prefactor_units[80] * fwd_A[80]
                * exp(fwd_beta[80] * tc[0] - activation_units[80] * fwd_Ea[80] * invT);
    dlnkfdT = fwd_beta[80] * invT + activation_units[80] * fwd_Ea[80] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[14];
    Kc = exp(-g_RT[0] + g_RT[4] + g_RT[4] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[4]) + (h_RT[0] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[4] -= 2 * q; /* CH3 */
    wdot[14] += q; /* C2H5 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[14];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[4] += -2 * dqdci;           /* dwdot[CH3]/d[H] */
    J[14] += dqdci;               /* dwdot[C2H5]/d[H] */
    /* d()/d[CH3] */
    dqdci =  + k_f*2.000000*sc[4];
    J[160] += dqdci;              /* dwdot[H]/d[CH3] */
    J[164] += -2 * dqdci;         /* dwdot[CH3]/d[CH3] */
    J[174] += dqdci;              /* dwdot[C2H5]/d[CH3] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[0];
    J[560] += dqdci;              /* dwdot[H]/d[C2H5] */
    J[564] += -2 * dqdci;         /* dwdot[CH3]/d[C2H5] */
    J[574] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    /* d()/dT */
    J[1560] += dqdT;              /* dwdot[H]/dT */
    J[1564] += -2 * dqdT;         /* dwdot[CH3]/dT */
    J[1574] += dqdT;              /* dwdot[C2H5]/dT */

    /*reaction 82: CH4 + CH2 <=> CH3 + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[6];
    k_f = prefactor_units[81] * fwd_A[81]
                * exp(fwd_beta[81] * tc[0] - activation_units[81] * fwd_Ea[81] * invT);
    dlnkfdT = fwd_beta[81] * invT + activation_units[81] * fwd_Ea[81] * invT2;
    /* reverse */
    phi_r = pow(sc[4], 2.000000);
    Kc = exp(g_RT[2] - g_RT[4] - g_RT[4] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[6]) + (2.000000*h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* CH2 */
    wdot[4] += 2 * q; /* CH3 */
    wdot[6] -= q; /* CH4 */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[6];
    J[82] -= dqdci;               /* dwdot[CH2]/d[CH2] */
    J[84] += 2 * dqdci;           /* dwdot[CH3]/d[CH2] */
    J[86] -= dqdci;               /* dwdot[CH4]/d[CH2] */
    /* d()/d[CH3] */
    dqdci =  - k_r*2.000000*sc[4];
    J[162] -= dqdci;              /* dwdot[CH2]/d[CH3] */
    J[164] += 2 * dqdci;          /* dwdot[CH3]/d[CH3] */
    J[166] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[2];
    J[242] -= dqdci;              /* dwdot[CH2]/d[CH4] */
    J[244] += 2 * dqdci;          /* dwdot[CH3]/d[CH4] */
    J[246] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[1562] -= dqdT;              /* dwdot[CH2]/dT */
    J[1564] += 2 * dqdT;          /* dwdot[CH3]/dT */
    J[1566] -= dqdT;              /* dwdot[CH4]/dT */

    /*reaction 83: CH4 + CH2(S) <=> CH3 + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[6];
    k_f = prefactor_units[82] * fwd_A[82]
                * exp(fwd_beta[82] * tc[0] - activation_units[82] * fwd_Ea[82] * invT);
    dlnkfdT = fwd_beta[82] * invT + activation_units[82] * fwd_Ea[82] * invT2;
    /* reverse */
    phi_r = pow(sc[4], 2.000000);
    Kc = exp(g_RT[3] - g_RT[4] - g_RT[4] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[6]) + (2.000000*h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* CH2(S) */
    wdot[4] += 2 * q; /* CH3 */
    wdot[6] -= q; /* CH4 */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[6];
    J[123] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[124] += 2 * dqdci;          /* dwdot[CH3]/d[CH2(S)] */
    J[126] -= dqdci;              /* dwdot[CH4]/d[CH2(S)] */
    /* d()/d[CH3] */
    dqdci =  - k_r*2.000000*sc[4];
    J[163] -= dqdci;              /* dwdot[CH2(S)]/d[CH3] */
    J[164] += 2 * dqdci;          /* dwdot[CH3]/d[CH3] */
    J[166] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[3];
    J[243] -= dqdci;              /* dwdot[CH2(S)]/d[CH4] */
    J[244] += 2 * dqdci;          /* dwdot[CH3]/d[CH4] */
    J[246] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[1563] -= dqdT;              /* dwdot[CH2(S)]/dT */
    J[1564] += 2 * dqdT;          /* dwdot[CH3]/dT */
    J[1566] -= dqdT;              /* dwdot[CH4]/dT */

    /*reaction 84: CH3 + OH <=> CH2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[7];
    k_f = prefactor_units[83] * fwd_A[83]
                * exp(fwd_beta[83] * tc[0] - activation_units[83] * fwd_Ea[83] * invT);
    dlnkfdT = fwd_beta[83] * invT + activation_units[83] * fwd_Ea[83] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[8];
    Kc = exp(-g_RT[2] + g_RT[4] + g_RT[7] - g_RT[8]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[7]) + (h_RT[2] + h_RT[8]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* CH2 */
    wdot[4] -= q; /* CH3 */
    wdot[7] -= q; /* OH */
    wdot[8] += q; /* H2O */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[8];
    J[82] += dqdci;               /* dwdot[CH2]/d[CH2] */
    J[84] -= dqdci;               /* dwdot[CH3]/d[CH2] */
    J[87] -= dqdci;               /* dwdot[OH]/d[CH2] */
    J[88] += dqdci;               /* dwdot[H2O]/d[CH2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[7];
    J[162] += dqdci;              /* dwdot[CH2]/d[CH3] */
    J[164] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[167] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[168] += dqdci;              /* dwdot[H2O]/d[CH3] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[4];
    J[282] += dqdci;              /* dwdot[CH2]/d[OH] */
    J[284] -= dqdci;              /* dwdot[CH3]/d[OH] */
    J[287] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[288] += dqdci;              /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[2];
    J[322] += dqdci;              /* dwdot[CH2]/d[H2O] */
    J[324] -= dqdci;              /* dwdot[CH3]/d[H2O] */
    J[327] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[328] += dqdci;              /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[1562] += dqdT;              /* dwdot[CH2]/dT */
    J[1564] -= dqdT;              /* dwdot[CH3]/dT */
    J[1567] -= dqdT;              /* dwdot[OH]/dT */
    J[1568] += dqdT;              /* dwdot[H2O]/dT */

    /*reaction 85: CH3 + OH <=> CH2(S) + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[7];
    k_f = prefactor_units[84] * fwd_A[84]
                * exp(fwd_beta[84] * tc[0] - activation_units[84] * fwd_Ea[84] * invT);
    dlnkfdT = fwd_beta[84] * invT + activation_units[84] * fwd_Ea[84] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[8];
    Kc = exp(-g_RT[3] + g_RT[4] + g_RT[7] - g_RT[8]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[7]) + (h_RT[3] + h_RT[8]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += q; /* CH2(S) */
    wdot[4] -= q; /* CH3 */
    wdot[7] -= q; /* OH */
    wdot[8] += q; /* H2O */
    /* d()/d[CH2(S)] */
    dqdci =  - k_r*sc[8];
    J[123] += dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[124] -= dqdci;              /* dwdot[CH3]/d[CH2(S)] */
    J[127] -= dqdci;              /* dwdot[OH]/d[CH2(S)] */
    J[128] += dqdci;              /* dwdot[H2O]/d[CH2(S)] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[7];
    J[163] += dqdci;              /* dwdot[CH2(S)]/d[CH3] */
    J[164] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[167] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[168] += dqdci;              /* dwdot[H2O]/d[CH3] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[4];
    J[283] += dqdci;              /* dwdot[CH2(S)]/d[OH] */
    J[284] -= dqdci;              /* dwdot[CH3]/d[OH] */
    J[287] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[288] += dqdci;              /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[3];
    J[323] += dqdci;              /* dwdot[CH2(S)]/d[H2O] */
    J[324] -= dqdci;              /* dwdot[CH3]/d[H2O] */
    J[327] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[328] += dqdci;              /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[1563] += dqdT;              /* dwdot[CH2(S)]/dT */
    J[1564] -= dqdT;              /* dwdot[CH3]/dT */
    J[1567] -= dqdT;              /* dwdot[OH]/dT */
    J[1568] += dqdT;              /* dwdot[H2O]/dT */

    /*reaction 86: CH3 + CH2 <=> C2H4 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[4];
    k_f = prefactor_units[85] * fwd_A[85]
                * exp(fwd_beta[85] * tc[0] - activation_units[85] * fwd_Ea[85] * invT);
    dlnkfdT = fwd_beta[85] * invT + activation_units[85] * fwd_Ea[85] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[12];
    Kc = exp(-g_RT[0] + g_RT[2] + g_RT[4] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[4]) + (h_RT[0] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[2] -= q; /* CH2 */
    wdot[4] -= q; /* CH3 */
    wdot[12] += q; /* C2H4 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[12];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[2] -= dqdci;                /* dwdot[CH2]/d[H] */
    J[4] -= dqdci;                /* dwdot[CH3]/d[H] */
    J[12] += dqdci;               /* dwdot[C2H4]/d[H] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[4];
    J[80] += dqdci;               /* dwdot[H]/d[CH2] */
    J[82] -= dqdci;               /* dwdot[CH2]/d[CH2] */
    J[84] -= dqdci;               /* dwdot[CH3]/d[CH2] */
    J[92] += dqdci;               /* dwdot[C2H4]/d[CH2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[2];
    J[160] += dqdci;              /* dwdot[H]/d[CH3] */
    J[162] -= dqdci;              /* dwdot[CH2]/d[CH3] */
    J[164] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[172] += dqdci;              /* dwdot[C2H4]/d[CH3] */
    /* d()/d[C2H4] */
    dqdci =  - k_r*sc[0];
    J[480] += dqdci;              /* dwdot[H]/d[C2H4] */
    J[482] -= dqdci;              /* dwdot[CH2]/d[C2H4] */
    J[484] -= dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[492] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
    /* d()/dT */
    J[1560] += dqdT;              /* dwdot[H]/dT */
    J[1562] -= dqdT;              /* dwdot[CH2]/dT */
    J[1564] -= dqdT;              /* dwdot[CH3]/dT */
    J[1572] += dqdT;              /* dwdot[C2H4]/dT */

    /*reaction 87: CH3 + CH2(S) <=> C2H4 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[4];
    k_f = prefactor_units[86] * fwd_A[86]
                * exp(fwd_beta[86] * tc[0] - activation_units[86] * fwd_Ea[86] * invT);
    dlnkfdT = fwd_beta[86] * invT + activation_units[86] * fwd_Ea[86] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[12];
    Kc = exp(-g_RT[0] + g_RT[3] + g_RT[4] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[4]) + (h_RT[0] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[3] -= q; /* CH2(S) */
    wdot[4] -= q; /* CH3 */
    wdot[12] += q; /* C2H4 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[12];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[3] -= dqdci;                /* dwdot[CH2(S)]/d[H] */
    J[4] -= dqdci;                /* dwdot[CH3]/d[H] */
    J[12] += dqdci;               /* dwdot[C2H4]/d[H] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[4];
    J[120] += dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[123] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[124] -= dqdci;              /* dwdot[CH3]/d[CH2(S)] */
    J[132] += dqdci;              /* dwdot[C2H4]/d[CH2(S)] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[3];
    J[160] += dqdci;              /* dwdot[H]/d[CH3] */
    J[163] -= dqdci;              /* dwdot[CH2(S)]/d[CH3] */
    J[164] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[172] += dqdci;              /* dwdot[C2H4]/d[CH3] */
    /* d()/d[C2H4] */
    dqdci =  - k_r*sc[0];
    J[480] += dqdci;              /* dwdot[H]/d[C2H4] */
    J[483] -= dqdci;              /* dwdot[CH2(S)]/d[C2H4] */
    J[484] -= dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[492] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
    /* d()/dT */
    J[1560] += dqdT;              /* dwdot[H]/dT */
    J[1563] -= dqdT;              /* dwdot[CH2(S)]/dT */
    J[1564] -= dqdT;              /* dwdot[CH3]/dT */
    J[1572] += dqdT;              /* dwdot[C2H4]/dT */

    /*reaction 88: CH3O + H <=> CH2(S) + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[18];
    k_f = prefactor_units[87] * fwd_A[87]
                * exp(fwd_beta[87] * tc[0] - activation_units[87] * fwd_Ea[87] * invT);
    dlnkfdT = fwd_beta[87] * invT + activation_units[87] * fwd_Ea[87] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[8];
    Kc = exp(g_RT[0] - g_RT[3] - g_RT[8] + g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[18]) + (h_RT[3] + h_RT[8]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[3] += q; /* CH2(S) */
    wdot[8] += q; /* H2O */
    wdot[18] -= q; /* CH3O */
    /* d()/d[H] */
    dqdci =  + k_f*sc[18];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[3] += dqdci;                /* dwdot[CH2(S)]/d[H] */
    J[8] += dqdci;                /* dwdot[H2O]/d[H] */
    J[18] -= dqdci;               /* dwdot[CH3O]/d[H] */
    /* d()/d[CH2(S)] */
    dqdci =  - k_r*sc[8];
    J[120] -= dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[123] += dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[128] += dqdci;              /* dwdot[H2O]/d[CH2(S)] */
    J[138] -= dqdci;              /* dwdot[CH3O]/d[CH2(S)] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[3];
    J[320] -= dqdci;              /* dwdot[H]/d[H2O] */
    J[323] += dqdci;              /* dwdot[CH2(S)]/d[H2O] */
    J[328] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[338] -= dqdci;              /* dwdot[CH3O]/d[H2O] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[0];
    J[720] -= dqdci;              /* dwdot[H]/d[CH3O] */
    J[723] += dqdci;              /* dwdot[CH2(S)]/d[CH3O] */
    J[728] += dqdci;              /* dwdot[H2O]/d[CH3O] */
    J[738] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[1560] -= dqdT;              /* dwdot[H]/dT */
    J[1563] += dqdT;              /* dwdot[CH2(S)]/dT */
    J[1568] += dqdT;              /* dwdot[H2O]/dT */
    J[1578] -= dqdT;              /* dwdot[CH3O]/dT */

    /*reaction 89: C2H6 + H <=> C2H5 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[16];
    k_f = prefactor_units[88] * fwd_A[88]
                * exp(fwd_beta[88] * tc[0] - activation_units[88] * fwd_Ea[88] * invT);
    dlnkfdT = fwd_beta[88] * invT + activation_units[88] * fwd_Ea[88] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[14];
    Kc = exp(g_RT[0] - g_RT[1] - g_RT[14] + g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[16]) + (h_RT[1] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[1] += q; /* H2 */
    wdot[14] += q; /* C2H5 */
    wdot[16] -= q; /* C2H6 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[16];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[1] += dqdci;                /* dwdot[H2]/d[H] */
    J[14] += dqdci;               /* dwdot[C2H5]/d[H] */
    J[16] -= dqdci;               /* dwdot[C2H6]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[14];
    J[40] -= dqdci;               /* dwdot[H]/d[H2] */
    J[41] += dqdci;               /* dwdot[H2]/d[H2] */
    J[54] += dqdci;               /* dwdot[C2H5]/d[H2] */
    J[56] -= dqdci;               /* dwdot[C2H6]/d[H2] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[1];
    J[560] -= dqdci;              /* dwdot[H]/d[C2H5] */
    J[561] += dqdci;              /* dwdot[H2]/d[C2H5] */
    J[574] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[576] -= dqdci;              /* dwdot[C2H6]/d[C2H5] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[0];
    J[640] -= dqdci;              /* dwdot[H]/d[C2H6] */
    J[641] += dqdci;              /* dwdot[H2]/d[C2H6] */
    J[654] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[656] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[1560] -= dqdT;              /* dwdot[H]/dT */
    J[1561] += dqdT;              /* dwdot[H2]/dT */
    J[1574] += dqdT;              /* dwdot[C2H5]/dT */
    J[1576] -= dqdT;              /* dwdot[C2H6]/dT */

    /*reaction 90: C2H6 + O <=> C2H5 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[16];
    k_f = prefactor_units[89] * fwd_A[89]
                * exp(fwd_beta[89] * tc[0] - activation_units[89] * fwd_Ea[89] * invT);
    dlnkfdT = fwd_beta[89] * invT + activation_units[89] * fwd_Ea[89] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[14];
    Kc = exp(g_RT[5] - g_RT[7] - g_RT[14] + g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[16]) + (h_RT[7] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] -= q; /* O */
    wdot[7] += q; /* OH */
    wdot[14] += q; /* C2H5 */
    wdot[16] -= q; /* C2H6 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[16];
    J[205] -= dqdci;              /* dwdot[O]/d[O] */
    J[207] += dqdci;              /* dwdot[OH]/d[O] */
    J[214] += dqdci;              /* dwdot[C2H5]/d[O] */
    J[216] -= dqdci;              /* dwdot[C2H6]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[14];
    J[285] -= dqdci;              /* dwdot[O]/d[OH] */
    J[287] += dqdci;              /* dwdot[OH]/d[OH] */
    J[294] += dqdci;              /* dwdot[C2H5]/d[OH] */
    J[296] -= dqdci;              /* dwdot[C2H6]/d[OH] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[7];
    J[565] -= dqdci;              /* dwdot[O]/d[C2H5] */
    J[567] += dqdci;              /* dwdot[OH]/d[C2H5] */
    J[574] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[576] -= dqdci;              /* dwdot[C2H6]/d[C2H5] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[5];
    J[645] -= dqdci;              /* dwdot[O]/d[C2H6] */
    J[647] += dqdci;              /* dwdot[OH]/d[C2H6] */
    J[654] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[656] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[1565] -= dqdT;              /* dwdot[O]/dT */
    J[1567] += dqdT;              /* dwdot[OH]/dT */
    J[1574] += dqdT;              /* dwdot[C2H5]/dT */
    J[1576] -= dqdT;              /* dwdot[C2H6]/dT */

    /*reaction 91: C2H6 + OH <=> C2H5 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[7]*sc[16];
    k_f = prefactor_units[90] * fwd_A[90]
                * exp(fwd_beta[90] * tc[0] - activation_units[90] * fwd_Ea[90] * invT);
    dlnkfdT = fwd_beta[90] * invT + activation_units[90] * fwd_Ea[90] * invT2;
    /* reverse */
    phi_r = sc[8]*sc[14];
    Kc = exp(g_RT[7] - g_RT[8] - g_RT[14] + g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7] + h_RT[16]) + (h_RT[8] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] -= q; /* OH */
    wdot[8] += q; /* H2O */
    wdot[14] += q; /* C2H5 */
    wdot[16] -= q; /* C2H6 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[16];
    J[287] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[288] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[294] += dqdci;              /* dwdot[C2H5]/d[OH] */
    J[296] -= dqdci;              /* dwdot[C2H6]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[14];
    J[327] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[328] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[334] += dqdci;              /* dwdot[C2H5]/d[H2O] */
    J[336] -= dqdci;              /* dwdot[C2H6]/d[H2O] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[8];
    J[567] -= dqdci;              /* dwdot[OH]/d[C2H5] */
    J[568] += dqdci;              /* dwdot[H2O]/d[C2H5] */
    J[574] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[576] -= dqdci;              /* dwdot[C2H6]/d[C2H5] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[7];
    J[647] -= dqdci;              /* dwdot[OH]/d[C2H6] */
    J[648] += dqdci;              /* dwdot[H2O]/d[C2H6] */
    J[654] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[656] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[1567] -= dqdT;              /* dwdot[OH]/dT */
    J[1568] += dqdT;              /* dwdot[H2O]/dT */
    J[1574] += dqdT;              /* dwdot[C2H5]/dT */
    J[1576] -= dqdT;              /* dwdot[C2H6]/dT */

    /*reaction 92: C2H6 + O2 <=> C2H5 + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[16]*sc[19];
    k_f = prefactor_units[91] * fwd_A[91]
                * exp(fwd_beta[91] * tc[0] - activation_units[91] * fwd_Ea[91] * invT);
    dlnkfdT = fwd_beta[91] * invT + activation_units[91] * fwd_Ea[91] * invT2;
    /* reverse */
    phi_r = sc[14]*sc[20];
    Kc = exp(-g_RT[14] + g_RT[16] + g_RT[19] - g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[16] + h_RT[19]) + (h_RT[14] + h_RT[20]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[14] += q; /* C2H5 */
    wdot[16] -= q; /* C2H6 */
    wdot[19] -= q; /* O2 */
    wdot[20] += q; /* HO2 */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[20];
    J[574] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[576] -= dqdci;              /* dwdot[C2H6]/d[C2H5] */
    J[579] -= dqdci;              /* dwdot[O2]/d[C2H5] */
    J[580] += dqdci;              /* dwdot[HO2]/d[C2H5] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[19];
    J[654] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[656] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    J[659] -= dqdci;              /* dwdot[O2]/d[C2H6] */
    J[660] += dqdci;              /* dwdot[HO2]/d[C2H6] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[16];
    J[774] += dqdci;              /* dwdot[C2H5]/d[O2] */
    J[776] -= dqdci;              /* dwdot[C2H6]/d[O2] */
    J[779] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[780] += dqdci;              /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[14];
    J[814] += dqdci;              /* dwdot[C2H5]/d[HO2] */
    J[816] -= dqdci;              /* dwdot[C2H6]/d[HO2] */
    J[819] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[820] += dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1574] += dqdT;              /* dwdot[C2H5]/dT */
    J[1576] -= dqdT;              /* dwdot[C2H6]/dT */
    J[1579] -= dqdT;              /* dwdot[O2]/dT */
    J[1580] += dqdT;              /* dwdot[HO2]/dT */

    /*reaction 93: C2H6 + HO2 <=> C2H5 + H2O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[16]*sc[20];
    k_f = prefactor_units[92] * fwd_A[92]
                * exp(fwd_beta[92] * tc[0] - activation_units[92] * fwd_Ea[92] * invT);
    dlnkfdT = fwd_beta[92] * invT + activation_units[92] * fwd_Ea[92] * invT2;
    /* reverse */
    phi_r = sc[14]*sc[21];
    Kc = exp(-g_RT[14] + g_RT[16] + g_RT[20] - g_RT[21]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[16] + h_RT[20]) + (h_RT[14] + h_RT[21]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[14] += q; /* C2H5 */
    wdot[16] -= q; /* C2H6 */
    wdot[20] -= q; /* HO2 */
    wdot[21] += q; /* H2O2 */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[21];
    J[574] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[576] -= dqdci;              /* dwdot[C2H6]/d[C2H5] */
    J[580] -= dqdci;              /* dwdot[HO2]/d[C2H5] */
    J[581] += dqdci;              /* dwdot[H2O2]/d[C2H5] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[20];
    J[654] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[656] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    J[660] -= dqdci;              /* dwdot[HO2]/d[C2H6] */
    J[661] += dqdci;              /* dwdot[H2O2]/d[C2H6] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[16];
    J[814] += dqdci;              /* dwdot[C2H5]/d[HO2] */
    J[816] -= dqdci;              /* dwdot[C2H6]/d[HO2] */
    J[820] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[821] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[14];
    J[854] += dqdci;              /* dwdot[C2H5]/d[H2O2] */
    J[856] -= dqdci;              /* dwdot[C2H6]/d[H2O2] */
    J[860] -= dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[861] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[1574] += dqdT;              /* dwdot[C2H5]/dT */
    J[1576] -= dqdT;              /* dwdot[C2H6]/dT */
    J[1580] -= dqdT;              /* dwdot[HO2]/dT */
    J[1581] += dqdT;              /* dwdot[H2O2]/dT */

    /*reaction 94: C2H6 + CH3 <=> C2H5 + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[16];
    k_f = prefactor_units[93] * fwd_A[93]
                * exp(fwd_beta[93] * tc[0] - activation_units[93] * fwd_Ea[93] * invT);
    dlnkfdT = fwd_beta[93] * invT + activation_units[93] * fwd_Ea[93] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[14];
    Kc = exp(g_RT[4] - g_RT[6] - g_RT[14] + g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[16]) + (h_RT[6] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* CH3 */
    wdot[6] += q; /* CH4 */
    wdot[14] += q; /* C2H5 */
    wdot[16] -= q; /* C2H6 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[16];
    J[164] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[166] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[174] += dqdci;              /* dwdot[C2H5]/d[CH3] */
    J[176] -= dqdci;              /* dwdot[C2H6]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[14];
    J[244] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[246] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[254] += dqdci;              /* dwdot[C2H5]/d[CH4] */
    J[256] -= dqdci;              /* dwdot[C2H6]/d[CH4] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[6];
    J[564] -= dqdci;              /* dwdot[CH3]/d[C2H5] */
    J[566] += dqdci;              /* dwdot[CH4]/d[C2H5] */
    J[574] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[576] -= dqdci;              /* dwdot[C2H6]/d[C2H5] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[4];
    J[644] -= dqdci;              /* dwdot[CH3]/d[C2H6] */
    J[646] += dqdci;              /* dwdot[CH4]/d[C2H6] */
    J[654] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[656] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[1564] -= dqdT;              /* dwdot[CH3]/dT */
    J[1566] += dqdT;              /* dwdot[CH4]/dT */
    J[1574] += dqdT;              /* dwdot[C2H5]/dT */
    J[1576] -= dqdT;              /* dwdot[C2H6]/dT */

    /*reaction 95: C2H5 + H <=> C2H4 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[14];
    k_f = prefactor_units[94] * fwd_A[94]
                * exp(fwd_beta[94] * tc[0] - activation_units[94] * fwd_Ea[94] * invT);
    dlnkfdT = fwd_beta[94] * invT + activation_units[94] * fwd_Ea[94] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[12];
    Kc = exp(g_RT[0] - g_RT[1] - g_RT[12] + g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[14]) + (h_RT[1] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[1] += q; /* H2 */
    wdot[12] += q; /* C2H4 */
    wdot[14] -= q; /* C2H5 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[14];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[1] += dqdci;                /* dwdot[H2]/d[H] */
    J[12] += dqdci;               /* dwdot[C2H4]/d[H] */
    J[14] -= dqdci;               /* dwdot[C2H5]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[12];
    J[40] -= dqdci;               /* dwdot[H]/d[H2] */
    J[41] += dqdci;               /* dwdot[H2]/d[H2] */
    J[52] += dqdci;               /* dwdot[C2H4]/d[H2] */
    J[54] -= dqdci;               /* dwdot[C2H5]/d[H2] */
    /* d()/d[C2H4] */
    dqdci =  - k_r*sc[1];
    J[480] -= dqdci;              /* dwdot[H]/d[C2H4] */
    J[481] += dqdci;              /* dwdot[H2]/d[C2H4] */
    J[492] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[494] -= dqdci;              /* dwdot[C2H5]/d[C2H4] */
    /* d()/d[C2H5] */
    dqdci =  + k_f*sc[0];
    J[560] -= dqdci;              /* dwdot[H]/d[C2H5] */
    J[561] += dqdci;              /* dwdot[H2]/d[C2H5] */
    J[572] += dqdci;              /* dwdot[C2H4]/d[C2H5] */
    J[574] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
    /* d()/dT */
    J[1560] -= dqdT;              /* dwdot[H]/dT */
    J[1561] += dqdT;              /* dwdot[H2]/dT */
    J[1572] += dqdT;              /* dwdot[C2H4]/dT */
    J[1574] -= dqdT;              /* dwdot[C2H5]/dT */

    /*reaction 96: C2H5 + O <=> CH3 + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[14];
    k_f = prefactor_units[95] * fwd_A[95]
                * exp(fwd_beta[95] * tc[0] - activation_units[95] * fwd_Ea[95] * invT);
    dlnkfdT = fwd_beta[95] * invT + activation_units[95] * fwd_Ea[95] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[15];
    Kc = exp(-g_RT[4] + g_RT[5] + g_RT[14] - g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[14]) + (h_RT[4] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] += q; /* CH3 */
    wdot[5] -= q; /* O */
    wdot[14] -= q; /* C2H5 */
    wdot[15] += q; /* CH2O */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[15];
    J[164] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[165] -= dqdci;              /* dwdot[O]/d[CH3] */
    J[174] -= dqdci;              /* dwdot[C2H5]/d[CH3] */
    J[175] += dqdci;              /* dwdot[CH2O]/d[CH3] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[14];
    J[204] += dqdci;              /* dwdot[CH3]/d[O] */
    J[205] -= dqdci;              /* dwdot[O]/d[O] */
    J[214] -= dqdci;              /* dwdot[C2H5]/d[O] */
    J[215] += dqdci;              /* dwdot[CH2O]/d[O] */
    /* d()/d[C2H5] */
    dqdci =  + k_f*sc[5];
    J[564] += dqdci;              /* dwdot[CH3]/d[C2H5] */
    J[565] -= dqdci;              /* dwdot[O]/d[C2H5] */
    J[574] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[575] += dqdci;              /* dwdot[CH2O]/d[C2H5] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[4];
    J[604] += dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[605] -= dqdci;              /* dwdot[O]/d[CH2O] */
    J[614] -= dqdci;              /* dwdot[C2H5]/d[CH2O] */
    J[615] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[1564] += dqdT;              /* dwdot[CH3]/dT */
    J[1565] -= dqdT;              /* dwdot[O]/dT */
    J[1574] -= dqdT;              /* dwdot[C2H5]/dT */
    J[1575] += dqdT;              /* dwdot[CH2O]/dT */

    /*reaction 97: C2H5 + O2 <=> C2H4 + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[14]*sc[19];
    k_f = prefactor_units[96] * fwd_A[96]
                * exp(fwd_beta[96] * tc[0] - activation_units[96] * fwd_Ea[96] * invT);
    dlnkfdT = fwd_beta[96] * invT + activation_units[96] * fwd_Ea[96] * invT2;
    /* reverse */
    phi_r = sc[12]*sc[20];
    Kc = exp(-g_RT[12] + g_RT[14] + g_RT[19] - g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[14] + h_RT[19]) + (h_RT[12] + h_RT[20]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[12] += q; /* C2H4 */
    wdot[14] -= q; /* C2H5 */
    wdot[19] -= q; /* O2 */
    wdot[20] += q; /* HO2 */
    /* d()/d[C2H4] */
    dqdci =  - k_r*sc[20];
    J[492] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[494] -= dqdci;              /* dwdot[C2H5]/d[C2H4] */
    J[499] -= dqdci;              /* dwdot[O2]/d[C2H4] */
    J[500] += dqdci;              /* dwdot[HO2]/d[C2H4] */
    /* d()/d[C2H5] */
    dqdci =  + k_f*sc[19];
    J[572] += dqdci;              /* dwdot[C2H4]/d[C2H5] */
    J[574] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[579] -= dqdci;              /* dwdot[O2]/d[C2H5] */
    J[580] += dqdci;              /* dwdot[HO2]/d[C2H5] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[14];
    J[772] += dqdci;              /* dwdot[C2H4]/d[O2] */
    J[774] -= dqdci;              /* dwdot[C2H5]/d[O2] */
    J[779] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[780] += dqdci;              /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[12];
    J[812] += dqdci;              /* dwdot[C2H4]/d[HO2] */
    J[814] -= dqdci;              /* dwdot[C2H5]/d[HO2] */
    J[819] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[820] += dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1572] += dqdT;              /* dwdot[C2H4]/dT */
    J[1574] -= dqdT;              /* dwdot[C2H5]/dT */
    J[1579] -= dqdT;              /* dwdot[O2]/dT */
    J[1580] += dqdT;              /* dwdot[HO2]/dT */

    /*reaction 98: C2H5 + C2H5 <=> C2H4 + C2H6 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[14], 2.000000);
    k_f = prefactor_units[97] * fwd_A[97]
                * exp(fwd_beta[97] * tc[0] - activation_units[97] * fwd_Ea[97] * invT);
    dlnkfdT = fwd_beta[97] * invT + activation_units[97] * fwd_Ea[97] * invT2;
    /* reverse */
    phi_r = sc[12]*sc[16];
    Kc = exp(-g_RT[12] + g_RT[14] + g_RT[14] - g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[14]) + (h_RT[12] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[12] += q; /* C2H4 */
    wdot[14] -= 2 * q; /* C2H5 */
    wdot[16] += q; /* C2H6 */
    /* d()/d[C2H4] */
    dqdci =  - k_r*sc[16];
    J[492] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[494] += -2 * dqdci;         /* dwdot[C2H5]/d[C2H4] */
    J[496] += dqdci;              /* dwdot[C2H6]/d[C2H4] */
    /* d()/d[C2H5] */
    dqdci =  + k_f*2.000000*sc[14];
    J[572] += dqdci;              /* dwdot[C2H4]/d[C2H5] */
    J[574] += -2 * dqdci;         /* dwdot[C2H5]/d[C2H5] */
    J[576] += dqdci;              /* dwdot[C2H6]/d[C2H5] */
    /* d()/d[C2H6] */
    dqdci =  - k_r*sc[12];
    J[652] += dqdci;              /* dwdot[C2H4]/d[C2H6] */
    J[654] += -2 * dqdci;         /* dwdot[C2H5]/d[C2H6] */
    J[656] += dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[1572] += dqdT;              /* dwdot[C2H4]/dT */
    J[1574] += -2 * dqdT;         /* dwdot[C2H5]/dT */
    J[1576] += dqdT;              /* dwdot[C2H6]/dT */

    /*reaction 99: C2H5 + HCO <=> C2H6 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[13]*sc[14];
    k_f = prefactor_units[98] * fwd_A[98]
                * exp(fwd_beta[98] * tc[0] - activation_units[98] * fwd_Ea[98] * invT);
    dlnkfdT = fwd_beta[98] * invT + activation_units[98] * fwd_Ea[98] * invT2;
    /* reverse */
    phi_r = sc[11]*sc[16];
    Kc = exp(-g_RT[11] + g_RT[13] + g_RT[14] - g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[13] + h_RT[14]) + (h_RT[11] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[11] += q; /* CO */
    wdot[13] -= q; /* HCO */
    wdot[14] -= q; /* C2H5 */
    wdot[16] += q; /* C2H6 */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[16];
    J[451] += dqdci;              /* dwdot[CO]/d[CO] */
    J[453] -= dqdci;              /* dwdot[HCO]/d[CO] */
    J[454] -= dqdci;              /* dwdot[C2H5]/d[CO] */
    J[456] += dqdci;              /* dwdot[C2H6]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[14];
    J[531] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[533] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    J[534] -= dqdci;              /* dwdot[C2H5]/d[HCO] */
    J[536] += dqdci;              /* dwdot[C2H6]/d[HCO] */
    /* d()/d[C2H5] */
    dqdci =  + k_f*sc[13];
    J[571] += dqdci;              /* dwdot[CO]/d[C2H5] */
    J[573] -= dqdci;              /* dwdot[HCO]/d[C2H5] */
    J[574] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[576] += dqdci;              /* dwdot[C2H6]/d[C2H5] */
    /* d()/d[C2H6] */
    dqdci =  - k_r*sc[11];
    J[651] += dqdci;              /* dwdot[CO]/d[C2H6] */
    J[653] -= dqdci;              /* dwdot[HCO]/d[C2H6] */
    J[654] -= dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[656] += dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[1571] += dqdT;              /* dwdot[CO]/dT */
    J[1573] -= dqdT;              /* dwdot[HCO]/dT */
    J[1574] -= dqdT;              /* dwdot[C2H5]/dT */
    J[1576] += dqdT;              /* dwdot[C2H6]/dT */

    /*reaction 100: C2H5 + O <=> CH3HCO + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[14];
    k_f = prefactor_units[99] * fwd_A[99]
                * exp(fwd_beta[99] * tc[0] - activation_units[99] * fwd_Ea[99] * invT);
    dlnkfdT = fwd_beta[99] * invT + activation_units[99] * fwd_Ea[99] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[23];
    Kc = exp(-g_RT[0] + g_RT[5] + g_RT[14] - g_RT[23]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[14]) + (h_RT[0] + h_RT[23]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[5] -= q; /* O */
    wdot[14] -= q; /* C2H5 */
    wdot[23] += q; /* CH3HCO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[23];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[5] -= dqdci;                /* dwdot[O]/d[H] */
    J[14] -= dqdci;               /* dwdot[C2H5]/d[H] */
    J[23] += dqdci;               /* dwdot[CH3HCO]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[14];
    J[200] += dqdci;              /* dwdot[H]/d[O] */
    J[205] -= dqdci;              /* dwdot[O]/d[O] */
    J[214] -= dqdci;              /* dwdot[C2H5]/d[O] */
    J[223] += dqdci;              /* dwdot[CH3HCO]/d[O] */
    /* d()/d[C2H5] */
    dqdci =  + k_f*sc[5];
    J[560] += dqdci;              /* dwdot[H]/d[C2H5] */
    J[565] -= dqdci;              /* dwdot[O]/d[C2H5] */
    J[574] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[583] += dqdci;              /* dwdot[CH3HCO]/d[C2H5] */
    /* d()/d[CH3HCO] */
    dqdci =  - k_r*sc[0];
    J[920] += dqdci;              /* dwdot[H]/d[CH3HCO] */
    J[925] -= dqdci;              /* dwdot[O]/d[CH3HCO] */
    J[934] -= dqdci;              /* dwdot[C2H5]/d[CH3HCO] */
    J[943] += dqdci;              /* dwdot[CH3HCO]/d[CH3HCO] */
    /* d()/dT */
    J[1560] += dqdT;              /* dwdot[H]/dT */
    J[1565] -= dqdT;              /* dwdot[O]/dT */
    J[1574] -= dqdT;              /* dwdot[C2H5]/dT */
    J[1583] += dqdT;              /* dwdot[CH3HCO]/dT */

    /*reaction 101: C2H4 + H <=> C2H3 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[12];
    k_f = prefactor_units[100] * fwd_A[100]
                * exp(fwd_beta[100] * tc[0] - activation_units[100] * fwd_Ea[100] * invT);
    dlnkfdT = fwd_beta[100] * invT + activation_units[100] * fwd_Ea[100] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[10];
    Kc = exp(g_RT[0] - g_RT[1] - g_RT[10] + g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[12]) + (h_RT[1] + h_RT[10]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[1] += q; /* H2 */
    wdot[10] += q; /* C2H3 */
    wdot[12] -= q; /* C2H4 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[12];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[1] += dqdci;                /* dwdot[H2]/d[H] */
    J[10] += dqdci;               /* dwdot[C2H3]/d[H] */
    J[12] -= dqdci;               /* dwdot[C2H4]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[10];
    J[40] -= dqdci;               /* dwdot[H]/d[H2] */
    J[41] += dqdci;               /* dwdot[H2]/d[H2] */
    J[50] += dqdci;               /* dwdot[C2H3]/d[H2] */
    J[52] -= dqdci;               /* dwdot[C2H4]/d[H2] */
    /* d()/d[C2H3] */
    dqdci =  - k_r*sc[1];
    J[400] -= dqdci;              /* dwdot[H]/d[C2H3] */
    J[401] += dqdci;              /* dwdot[H2]/d[C2H3] */
    J[410] += dqdci;              /* dwdot[C2H3]/d[C2H3] */
    J[412] -= dqdci;              /* dwdot[C2H4]/d[C2H3] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[0];
    J[480] -= dqdci;              /* dwdot[H]/d[C2H4] */
    J[481] += dqdci;              /* dwdot[H2]/d[C2H4] */
    J[490] += dqdci;              /* dwdot[C2H3]/d[C2H4] */
    J[492] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    /* d()/dT */
    J[1560] -= dqdT;              /* dwdot[H]/dT */
    J[1561] += dqdT;              /* dwdot[H2]/dT */
    J[1570] += dqdT;              /* dwdot[C2H3]/dT */
    J[1572] -= dqdT;              /* dwdot[C2H4]/dT */

    /*reaction 102: C2H4 + OH <=> C2H3 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[7]*sc[12];
    k_f = prefactor_units[101] * fwd_A[101]
                * exp(fwd_beta[101] * tc[0] - activation_units[101] * fwd_Ea[101] * invT);
    dlnkfdT = fwd_beta[101] * invT + activation_units[101] * fwd_Ea[101] * invT2;
    /* reverse */
    phi_r = sc[8]*sc[10];
    Kc = exp(g_RT[7] - g_RT[8] - g_RT[10] + g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7] + h_RT[12]) + (h_RT[8] + h_RT[10]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] -= q; /* OH */
    wdot[8] += q; /* H2O */
    wdot[10] += q; /* C2H3 */
    wdot[12] -= q; /* C2H4 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[12];
    J[287] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[288] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[290] += dqdci;              /* dwdot[C2H3]/d[OH] */
    J[292] -= dqdci;              /* dwdot[C2H4]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[10];
    J[327] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[328] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[330] += dqdci;              /* dwdot[C2H3]/d[H2O] */
    J[332] -= dqdci;              /* dwdot[C2H4]/d[H2O] */
    /* d()/d[C2H3] */
    dqdci =  - k_r*sc[8];
    J[407] -= dqdci;              /* dwdot[OH]/d[C2H3] */
    J[408] += dqdci;              /* dwdot[H2O]/d[C2H3] */
    J[410] += dqdci;              /* dwdot[C2H3]/d[C2H3] */
    J[412] -= dqdci;              /* dwdot[C2H4]/d[C2H3] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[7];
    J[487] -= dqdci;              /* dwdot[OH]/d[C2H4] */
    J[488] += dqdci;              /* dwdot[H2O]/d[C2H4] */
    J[490] += dqdci;              /* dwdot[C2H3]/d[C2H4] */
    J[492] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    /* d()/dT */
    J[1567] -= dqdT;              /* dwdot[OH]/dT */
    J[1568] += dqdT;              /* dwdot[H2O]/dT */
    J[1570] += dqdT;              /* dwdot[C2H3]/dT */
    J[1572] -= dqdT;              /* dwdot[C2H4]/dT */

    /*reaction 103: C2H4 + CH3 <=> C2H3 + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[12];
    k_f = prefactor_units[102] * fwd_A[102]
                * exp(fwd_beta[102] * tc[0] - activation_units[102] * fwd_Ea[102] * invT);
    dlnkfdT = fwd_beta[102] * invT + activation_units[102] * fwd_Ea[102] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[10];
    Kc = exp(g_RT[4] - g_RT[6] - g_RT[10] + g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[12]) + (h_RT[6] + h_RT[10]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* CH3 */
    wdot[6] += q; /* CH4 */
    wdot[10] += q; /* C2H3 */
    wdot[12] -= q; /* C2H4 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[12];
    J[164] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[166] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[170] += dqdci;              /* dwdot[C2H3]/d[CH3] */
    J[172] -= dqdci;              /* dwdot[C2H4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[10];
    J[244] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[246] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[250] += dqdci;              /* dwdot[C2H3]/d[CH4] */
    J[252] -= dqdci;              /* dwdot[C2H4]/d[CH4] */
    /* d()/d[C2H3] */
    dqdci =  - k_r*sc[6];
    J[404] -= dqdci;              /* dwdot[CH3]/d[C2H3] */
    J[406] += dqdci;              /* dwdot[CH4]/d[C2H3] */
    J[410] += dqdci;              /* dwdot[C2H3]/d[C2H3] */
    J[412] -= dqdci;              /* dwdot[C2H4]/d[C2H3] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[4];
    J[484] -= dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[486] += dqdci;              /* dwdot[CH4]/d[C2H4] */
    J[490] += dqdci;              /* dwdot[C2H3]/d[C2H4] */
    J[492] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    /* d()/dT */
    J[1564] -= dqdT;              /* dwdot[CH3]/dT */
    J[1566] += dqdT;              /* dwdot[CH4]/dT */
    J[1570] += dqdT;              /* dwdot[C2H3]/dT */
    J[1572] -= dqdT;              /* dwdot[C2H4]/dT */

    /*reaction 104: C2H4 + O <=> CH3 + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[12];
    k_f = prefactor_units[103] * fwd_A[103]
                * exp(fwd_beta[103] * tc[0] - activation_units[103] * fwd_Ea[103] * invT);
    dlnkfdT = fwd_beta[103] * invT + activation_units[103] * fwd_Ea[103] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[13];
    Kc = exp(-g_RT[4] + g_RT[5] + g_RT[12] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[12]) + (h_RT[4] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] += q; /* CH3 */
    wdot[5] -= q; /* O */
    wdot[12] -= q; /* C2H4 */
    wdot[13] += q; /* HCO */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[13];
    J[164] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[165] -= dqdci;              /* dwdot[O]/d[CH3] */
    J[172] -= dqdci;              /* dwdot[C2H4]/d[CH3] */
    J[173] += dqdci;              /* dwdot[HCO]/d[CH3] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[12];
    J[204] += dqdci;              /* dwdot[CH3]/d[O] */
    J[205] -= dqdci;              /* dwdot[O]/d[O] */
    J[212] -= dqdci;              /* dwdot[C2H4]/d[O] */
    J[213] += dqdci;              /* dwdot[HCO]/d[O] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[5];
    J[484] += dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[485] -= dqdci;              /* dwdot[O]/d[C2H4] */
    J[492] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[493] += dqdci;              /* dwdot[HCO]/d[C2H4] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[4];
    J[524] += dqdci;              /* dwdot[CH3]/d[HCO] */
    J[525] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[532] -= dqdci;              /* dwdot[C2H4]/d[HCO] */
    J[533] += dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[1564] += dqdT;              /* dwdot[CH3]/dT */
    J[1565] -= dqdT;              /* dwdot[O]/dT */
    J[1572] -= dqdT;              /* dwdot[C2H4]/dT */
    J[1573] += dqdT;              /* dwdot[HCO]/dT */

    /*reaction 105: C2H3 + OH <=> C2H2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[7]*sc[10];
    k_f = prefactor_units[104] * fwd_A[104]
                * exp(fwd_beta[104] * tc[0] - activation_units[104] * fwd_Ea[104] * invT);
    dlnkfdT = fwd_beta[104] * invT + activation_units[104] * fwd_Ea[104] * invT2;
    /* reverse */
    phi_r = sc[8]*sc[9];
    Kc = exp(g_RT[7] - g_RT[8] - g_RT[9] + g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7] + h_RT[10]) + (h_RT[8] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] -= q; /* OH */
    wdot[8] += q; /* H2O */
    wdot[9] += q; /* C2H2 */
    wdot[10] -= q; /* C2H3 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[10];
    J[287] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[288] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[289] += dqdci;              /* dwdot[C2H2]/d[OH] */
    J[290] -= dqdci;              /* dwdot[C2H3]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[9];
    J[327] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[328] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[329] += dqdci;              /* dwdot[C2H2]/d[H2O] */
    J[330] -= dqdci;              /* dwdot[C2H3]/d[H2O] */
    /* d()/d[C2H2] */
    dqdci =  - k_r*sc[8];
    J[367] -= dqdci;              /* dwdot[OH]/d[C2H2] */
    J[368] += dqdci;              /* dwdot[H2O]/d[C2H2] */
    J[369] += dqdci;              /* dwdot[C2H2]/d[C2H2] */
    J[370] -= dqdci;              /* dwdot[C2H3]/d[C2H2] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[7];
    J[407] -= dqdci;              /* dwdot[OH]/d[C2H3] */
    J[408] += dqdci;              /* dwdot[H2O]/d[C2H3] */
    J[409] += dqdci;              /* dwdot[C2H2]/d[C2H3] */
    J[410] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    /* d()/dT */
    J[1567] -= dqdT;              /* dwdot[OH]/dT */
    J[1568] += dqdT;              /* dwdot[H2O]/dT */
    J[1569] += dqdT;              /* dwdot[C2H2]/dT */
    J[1570] -= dqdT;              /* dwdot[C2H3]/dT */

    /*reaction 106: C2H4 + O <=> OH + C2H3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[12];
    k_f = prefactor_units[105] * fwd_A[105]
                * exp(fwd_beta[105] * tc[0] - activation_units[105] * fwd_Ea[105] * invT);
    dlnkfdT = fwd_beta[105] * invT + activation_units[105] * fwd_Ea[105] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[10];
    Kc = exp(g_RT[5] - g_RT[7] - g_RT[10] + g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[12]) + (h_RT[7] + h_RT[10]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] -= q; /* O */
    wdot[7] += q; /* OH */
    wdot[10] += q; /* C2H3 */
    wdot[12] -= q; /* C2H4 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[12];
    J[205] -= dqdci;              /* dwdot[O]/d[O] */
    J[207] += dqdci;              /* dwdot[OH]/d[O] */
    J[210] += dqdci;              /* dwdot[C2H3]/d[O] */
    J[212] -= dqdci;              /* dwdot[C2H4]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[10];
    J[285] -= dqdci;              /* dwdot[O]/d[OH] */
    J[287] += dqdci;              /* dwdot[OH]/d[OH] */
    J[290] += dqdci;              /* dwdot[C2H3]/d[OH] */
    J[292] -= dqdci;              /* dwdot[C2H4]/d[OH] */
    /* d()/d[C2H3] */
    dqdci =  - k_r*sc[7];
    J[405] -= dqdci;              /* dwdot[O]/d[C2H3] */
    J[407] += dqdci;              /* dwdot[OH]/d[C2H3] */
    J[410] += dqdci;              /* dwdot[C2H3]/d[C2H3] */
    J[412] -= dqdci;              /* dwdot[C2H4]/d[C2H3] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[5];
    J[485] -= dqdci;              /* dwdot[O]/d[C2H4] */
    J[487] += dqdci;              /* dwdot[OH]/d[C2H4] */
    J[490] += dqdci;              /* dwdot[C2H3]/d[C2H4] */
    J[492] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    /* d()/dT */
    J[1565] -= dqdT;              /* dwdot[O]/dT */
    J[1567] += dqdT;              /* dwdot[OH]/dT */
    J[1570] += dqdT;              /* dwdot[C2H3]/dT */
    J[1572] -= dqdT;              /* dwdot[C2H4]/dT */

    /*reaction 107: C2H4 + O2 <=> C2H3 + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[12]*sc[19];
    k_f = prefactor_units[106] * fwd_A[106]
                * exp(fwd_beta[106] * tc[0] - activation_units[106] * fwd_Ea[106] * invT);
    dlnkfdT = fwd_beta[106] * invT + activation_units[106] * fwd_Ea[106] * invT2;
    /* reverse */
    phi_r = sc[10]*sc[20];
    Kc = exp(-g_RT[10] + g_RT[12] + g_RT[19] - g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[12] + h_RT[19]) + (h_RT[10] + h_RT[20]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[10] += q; /* C2H3 */
    wdot[12] -= q; /* C2H4 */
    wdot[19] -= q; /* O2 */
    wdot[20] += q; /* HO2 */
    /* d()/d[C2H3] */
    dqdci =  - k_r*sc[20];
    J[410] += dqdci;              /* dwdot[C2H3]/d[C2H3] */
    J[412] -= dqdci;              /* dwdot[C2H4]/d[C2H3] */
    J[419] -= dqdci;              /* dwdot[O2]/d[C2H3] */
    J[420] += dqdci;              /* dwdot[HO2]/d[C2H3] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[19];
    J[490] += dqdci;              /* dwdot[C2H3]/d[C2H4] */
    J[492] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[499] -= dqdci;              /* dwdot[O2]/d[C2H4] */
    J[500] += dqdci;              /* dwdot[HO2]/d[C2H4] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[12];
    J[770] += dqdci;              /* dwdot[C2H3]/d[O2] */
    J[772] -= dqdci;              /* dwdot[C2H4]/d[O2] */
    J[779] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[780] += dqdci;              /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[10];
    J[810] += dqdci;              /* dwdot[C2H3]/d[HO2] */
    J[812] -= dqdci;              /* dwdot[C2H4]/d[HO2] */
    J[819] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[820] += dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1570] += dqdT;              /* dwdot[C2H3]/dT */
    J[1572] -= dqdT;              /* dwdot[C2H4]/dT */
    J[1579] -= dqdT;              /* dwdot[O2]/dT */
    J[1580] += dqdT;              /* dwdot[HO2]/dT */

    /*reaction 108: C2H3 + H <=> C2H2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[10];
    k_f = prefactor_units[107] * fwd_A[107]
                * exp(fwd_beta[107] * tc[0] - activation_units[107] * fwd_Ea[107] * invT);
    dlnkfdT = fwd_beta[107] * invT + activation_units[107] * fwd_Ea[107] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[9];
    Kc = exp(g_RT[0] - g_RT[1] - g_RT[9] + g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[10]) + (h_RT[1] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[1] += q; /* H2 */
    wdot[9] += q; /* C2H2 */
    wdot[10] -= q; /* C2H3 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[10];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[1] += dqdci;                /* dwdot[H2]/d[H] */
    J[9] += dqdci;                /* dwdot[C2H2]/d[H] */
    J[10] -= dqdci;               /* dwdot[C2H3]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[9];
    J[40] -= dqdci;               /* dwdot[H]/d[H2] */
    J[41] += dqdci;               /* dwdot[H2]/d[H2] */
    J[49] += dqdci;               /* dwdot[C2H2]/d[H2] */
    J[50] -= dqdci;               /* dwdot[C2H3]/d[H2] */
    /* d()/d[C2H2] */
    dqdci =  - k_r*sc[1];
    J[360] -= dqdci;              /* dwdot[H]/d[C2H2] */
    J[361] += dqdci;              /* dwdot[H2]/d[C2H2] */
    J[369] += dqdci;              /* dwdot[C2H2]/d[C2H2] */
    J[370] -= dqdci;              /* dwdot[C2H3]/d[C2H2] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[0];
    J[400] -= dqdci;              /* dwdot[H]/d[C2H3] */
    J[401] += dqdci;              /* dwdot[H2]/d[C2H3] */
    J[409] += dqdci;              /* dwdot[C2H2]/d[C2H3] */
    J[410] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    /* d()/dT */
    J[1560] -= dqdT;              /* dwdot[H]/dT */
    J[1561] += dqdT;              /* dwdot[H2]/dT */
    J[1569] += dqdT;              /* dwdot[C2H2]/dT */
    J[1570] -= dqdT;              /* dwdot[C2H3]/dT */

    /*reaction 109: C2H3 + H2O2 <=> C2H4 + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[10]*sc[21];
    k_f = prefactor_units[108] * fwd_A[108]
                * exp(fwd_beta[108] * tc[0] - activation_units[108] * fwd_Ea[108] * invT);
    dlnkfdT = fwd_beta[108] * invT + activation_units[108] * fwd_Ea[108] * invT2;
    /* reverse */
    phi_r = sc[12]*sc[20];
    Kc = exp(g_RT[10] - g_RT[12] - g_RT[20] + g_RT[21]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[10] + h_RT[21]) + (h_RT[12] + h_RT[20]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[10] -= q; /* C2H3 */
    wdot[12] += q; /* C2H4 */
    wdot[20] += q; /* HO2 */
    wdot[21] -= q; /* H2O2 */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[21];
    J[410] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    J[412] += dqdci;              /* dwdot[C2H4]/d[C2H3] */
    J[420] += dqdci;              /* dwdot[HO2]/d[C2H3] */
    J[421] -= dqdci;              /* dwdot[H2O2]/d[C2H3] */
    /* d()/d[C2H4] */
    dqdci =  - k_r*sc[20];
    J[490] -= dqdci;              /* dwdot[C2H3]/d[C2H4] */
    J[492] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[500] += dqdci;              /* dwdot[HO2]/d[C2H4] */
    J[501] -= dqdci;              /* dwdot[H2O2]/d[C2H4] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[12];
    J[810] -= dqdci;              /* dwdot[C2H3]/d[HO2] */
    J[812] += dqdci;              /* dwdot[C2H4]/d[HO2] */
    J[820] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[821] -= dqdci;              /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[10];
    J[850] -= dqdci;              /* dwdot[C2H3]/d[H2O2] */
    J[852] += dqdci;              /* dwdot[C2H4]/d[H2O2] */
    J[860] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[861] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[1570] -= dqdT;              /* dwdot[C2H3]/dT */
    J[1572] += dqdT;              /* dwdot[C2H4]/dT */
    J[1580] += dqdT;              /* dwdot[HO2]/dT */
    J[1581] -= dqdT;              /* dwdot[H2O2]/dT */

    /*reaction 110: C2H3 + CH3 <=> C2H2 + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[10];
    k_f = prefactor_units[109] * fwd_A[109]
                * exp(fwd_beta[109] * tc[0] - activation_units[109] * fwd_Ea[109] * invT);
    dlnkfdT = fwd_beta[109] * invT + activation_units[109] * fwd_Ea[109] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[9];
    Kc = exp(g_RT[4] - g_RT[6] - g_RT[9] + g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[10]) + (h_RT[6] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* CH3 */
    wdot[6] += q; /* CH4 */
    wdot[9] += q; /* C2H2 */
    wdot[10] -= q; /* C2H3 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[10];
    J[164] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[166] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[169] += dqdci;              /* dwdot[C2H2]/d[CH3] */
    J[170] -= dqdci;              /* dwdot[C2H3]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[9];
    J[244] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[246] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[249] += dqdci;              /* dwdot[C2H2]/d[CH4] */
    J[250] -= dqdci;              /* dwdot[C2H3]/d[CH4] */
    /* d()/d[C2H2] */
    dqdci =  - k_r*sc[6];
    J[364] -= dqdci;              /* dwdot[CH3]/d[C2H2] */
    J[366] += dqdci;              /* dwdot[CH4]/d[C2H2] */
    J[369] += dqdci;              /* dwdot[C2H2]/d[C2H2] */
    J[370] -= dqdci;              /* dwdot[C2H3]/d[C2H2] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[4];
    J[404] -= dqdci;              /* dwdot[CH3]/d[C2H3] */
    J[406] += dqdci;              /* dwdot[CH4]/d[C2H3] */
    J[409] += dqdci;              /* dwdot[C2H2]/d[C2H3] */
    J[410] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    /* d()/dT */
    J[1564] -= dqdT;              /* dwdot[CH3]/dT */
    J[1566] += dqdT;              /* dwdot[CH4]/dT */
    J[1569] += dqdT;              /* dwdot[C2H2]/dT */
    J[1570] -= dqdT;              /* dwdot[C2H3]/dT */

    /*reaction 111: C2H3 + C2H3 <=> C2H4 + C2H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[10], 2.000000);
    k_f = prefactor_units[110] * fwd_A[110]
                * exp(fwd_beta[110] * tc[0] - activation_units[110] * fwd_Ea[110] * invT);
    dlnkfdT = fwd_beta[110] * invT + activation_units[110] * fwd_Ea[110] * invT2;
    /* reverse */
    phi_r = sc[9]*sc[12];
    Kc = exp(-g_RT[9] + g_RT[10] + g_RT[10] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[10]) + (h_RT[9] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[9] += q; /* C2H2 */
    wdot[10] -= 2 * q; /* C2H3 */
    wdot[12] += q; /* C2H4 */
    /* d()/d[C2H2] */
    dqdci =  - k_r*sc[12];
    J[369] += dqdci;              /* dwdot[C2H2]/d[C2H2] */
    J[370] += -2 * dqdci;         /* dwdot[C2H3]/d[C2H2] */
    J[372] += dqdci;              /* dwdot[C2H4]/d[C2H2] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*2.000000*sc[10];
    J[409] += dqdci;              /* dwdot[C2H2]/d[C2H3] */
    J[410] += -2 * dqdci;         /* dwdot[C2H3]/d[C2H3] */
    J[412] += dqdci;              /* dwdot[C2H4]/d[C2H3] */
    /* d()/d[C2H4] */
    dqdci =  - k_r*sc[9];
    J[489] += dqdci;              /* dwdot[C2H2]/d[C2H4] */
    J[490] += -2 * dqdci;         /* dwdot[C2H3]/d[C2H4] */
    J[492] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
    /* d()/dT */
    J[1569] += dqdT;              /* dwdot[C2H2]/dT */
    J[1570] += -2 * dqdT;         /* dwdot[C2H3]/dT */
    J[1572] += dqdT;              /* dwdot[C2H4]/dT */

    /*reaction 112: C2H3 + O2 <=> HCO + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[10]*sc[19];
    k_f = prefactor_units[111] * fwd_A[111]
                * exp(fwd_beta[111] * tc[0] - activation_units[111] * fwd_Ea[111] * invT);
    dlnkfdT = fwd_beta[111] * invT + activation_units[111] * fwd_Ea[111] * invT2;
    /* reverse */
    phi_r = sc[13]*sc[15];
    Kc = exp(g_RT[10] - g_RT[13] - g_RT[15] + g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[10] + h_RT[19]) + (h_RT[13] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[10] -= q; /* C2H3 */
    wdot[13] += q; /* HCO */
    wdot[15] += q; /* CH2O */
    wdot[19] -= q; /* O2 */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[19];
    J[410] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    J[413] += dqdci;              /* dwdot[HCO]/d[C2H3] */
    J[415] += dqdci;              /* dwdot[CH2O]/d[C2H3] */
    J[419] -= dqdci;              /* dwdot[O2]/d[C2H3] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[15];
    J[530] -= dqdci;              /* dwdot[C2H3]/d[HCO] */
    J[533] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[535] += dqdci;              /* dwdot[CH2O]/d[HCO] */
    J[539] -= dqdci;              /* dwdot[O2]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[13];
    J[610] -= dqdci;              /* dwdot[C2H3]/d[CH2O] */
    J[613] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[615] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[619] -= dqdci;              /* dwdot[O2]/d[CH2O] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[10];
    J[770] -= dqdci;              /* dwdot[C2H3]/d[O2] */
    J[773] += dqdci;              /* dwdot[HCO]/d[O2] */
    J[775] += dqdci;              /* dwdot[CH2O]/d[O2] */
    J[779] -= dqdci;              /* dwdot[O2]/d[O2] */
    /* d()/dT */
    J[1570] -= dqdT;              /* dwdot[C2H3]/dT */
    J[1573] += dqdT;              /* dwdot[HCO]/dT */
    J[1575] += dqdT;              /* dwdot[CH2O]/dT */
    J[1579] -= dqdT;              /* dwdot[O2]/dT */

    /*reaction 113: C2H3 + O2 <=> HO2 + C2H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[10]*sc[19];
    k_f = prefactor_units[112] * fwd_A[112]
                * exp(fwd_beta[112] * tc[0] - activation_units[112] * fwd_Ea[112] * invT);
    dlnkfdT = fwd_beta[112] * invT + activation_units[112] * fwd_Ea[112] * invT2;
    /* reverse */
    phi_r = sc[9]*sc[20];
    Kc = exp(-g_RT[9] + g_RT[10] + g_RT[19] - g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[10] + h_RT[19]) + (h_RT[9] + h_RT[20]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[9] += q; /* C2H2 */
    wdot[10] -= q; /* C2H3 */
    wdot[19] -= q; /* O2 */
    wdot[20] += q; /* HO2 */
    /* d()/d[C2H2] */
    dqdci =  - k_r*sc[20];
    J[369] += dqdci;              /* dwdot[C2H2]/d[C2H2] */
    J[370] -= dqdci;              /* dwdot[C2H3]/d[C2H2] */
    J[379] -= dqdci;              /* dwdot[O2]/d[C2H2] */
    J[380] += dqdci;              /* dwdot[HO2]/d[C2H2] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[19];
    J[409] += dqdci;              /* dwdot[C2H2]/d[C2H3] */
    J[410] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    J[419] -= dqdci;              /* dwdot[O2]/d[C2H3] */
    J[420] += dqdci;              /* dwdot[HO2]/d[C2H3] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[10];
    J[769] += dqdci;              /* dwdot[C2H2]/d[O2] */
    J[770] -= dqdci;              /* dwdot[C2H3]/d[O2] */
    J[779] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[780] += dqdci;              /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[9];
    J[809] += dqdci;              /* dwdot[C2H2]/d[HO2] */
    J[810] -= dqdci;              /* dwdot[C2H3]/d[HO2] */
    J[819] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[820] += dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1569] += dqdT;              /* dwdot[C2H2]/dT */
    J[1570] -= dqdT;              /* dwdot[C2H3]/dT */
    J[1579] -= dqdT;              /* dwdot[O2]/dT */
    J[1580] += dqdT;              /* dwdot[HO2]/dT */

    /*reaction 114: C2H2 + O <=> CH2 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[9];
    k_f = prefactor_units[113] * fwd_A[113]
                * exp(fwd_beta[113] * tc[0] - activation_units[113] * fwd_Ea[113] * invT);
    dlnkfdT = fwd_beta[113] * invT + activation_units[113] * fwd_Ea[113] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[11];
    Kc = exp(-g_RT[2] + g_RT[5] + g_RT[9] - g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[9]) + (h_RT[2] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* CH2 */
    wdot[5] -= q; /* O */
    wdot[9] -= q; /* C2H2 */
    wdot[11] += q; /* CO */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[11];
    J[82] += dqdci;               /* dwdot[CH2]/d[CH2] */
    J[85] -= dqdci;               /* dwdot[O]/d[CH2] */
    J[89] -= dqdci;               /* dwdot[C2H2]/d[CH2] */
    J[91] += dqdci;               /* dwdot[CO]/d[CH2] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[9];
    J[202] += dqdci;              /* dwdot[CH2]/d[O] */
    J[205] -= dqdci;              /* dwdot[O]/d[O] */
    J[209] -= dqdci;              /* dwdot[C2H2]/d[O] */
    J[211] += dqdci;              /* dwdot[CO]/d[O] */
    /* d()/d[C2H2] */
    dqdci =  + k_f*sc[5];
    J[362] += dqdci;              /* dwdot[CH2]/d[C2H2] */
    J[365] -= dqdci;              /* dwdot[O]/d[C2H2] */
    J[369] -= dqdci;              /* dwdot[C2H2]/d[C2H2] */
    J[371] += dqdci;              /* dwdot[CO]/d[C2H2] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[2];
    J[442] += dqdci;              /* dwdot[CH2]/d[CO] */
    J[445] -= dqdci;              /* dwdot[O]/d[CO] */
    J[449] -= dqdci;              /* dwdot[C2H2]/d[CO] */
    J[451] += dqdci;              /* dwdot[CO]/d[CO] */
    /* d()/dT */
    J[1562] += dqdT;              /* dwdot[CH2]/dT */
    J[1565] -= dqdT;              /* dwdot[O]/dT */
    J[1569] -= dqdT;              /* dwdot[C2H2]/dT */
    J[1571] += dqdT;              /* dwdot[CO]/dT */

    /*reaction 115: C2H2 + OH <=> CH3 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[7]*sc[9];
    k_f = prefactor_units[114] * fwd_A[114]
                * exp(fwd_beta[114] * tc[0] - activation_units[114] * fwd_Ea[114] * invT);
    dlnkfdT = fwd_beta[114] * invT + activation_units[114] * fwd_Ea[114] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[11];
    Kc = exp(-g_RT[4] + g_RT[7] + g_RT[9] - g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7] + h_RT[9]) + (h_RT[4] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] += q; /* CH3 */
    wdot[7] -= q; /* OH */
    wdot[9] -= q; /* C2H2 */
    wdot[11] += q; /* CO */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[11];
    J[164] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[167] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[169] -= dqdci;              /* dwdot[C2H2]/d[CH3] */
    J[171] += dqdci;              /* dwdot[CO]/d[CH3] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[9];
    J[284] += dqdci;              /* dwdot[CH3]/d[OH] */
    J[287] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[289] -= dqdci;              /* dwdot[C2H2]/d[OH] */
    J[291] += dqdci;              /* dwdot[CO]/d[OH] */
    /* d()/d[C2H2] */
    dqdci =  + k_f*sc[7];
    J[364] += dqdci;              /* dwdot[CH3]/d[C2H2] */
    J[367] -= dqdci;              /* dwdot[OH]/d[C2H2] */
    J[369] -= dqdci;              /* dwdot[C2H2]/d[C2H2] */
    J[371] += dqdci;              /* dwdot[CO]/d[C2H2] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[4];
    J[444] += dqdci;              /* dwdot[CH3]/d[CO] */
    J[447] -= dqdci;              /* dwdot[OH]/d[CO] */
    J[449] -= dqdci;              /* dwdot[C2H2]/d[CO] */
    J[451] += dqdci;              /* dwdot[CO]/d[CO] */
    /* d()/dT */
    J[1564] += dqdT;              /* dwdot[CH3]/dT */
    J[1567] -= dqdT;              /* dwdot[OH]/dT */
    J[1569] -= dqdT;              /* dwdot[C2H2]/dT */
    J[1571] += dqdT;              /* dwdot[CO]/dT */

    /*reaction 116: CH2 + O <=> HCO + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[5];
    k_f = prefactor_units[115] * fwd_A[115]
                * exp(fwd_beta[115] * tc[0] - activation_units[115] * fwd_Ea[115] * invT);
    dlnkfdT = fwd_beta[115] * invT + activation_units[115] * fwd_Ea[115] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[13];
    Kc = exp(-g_RT[0] + g_RT[2] + g_RT[5] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[5]) + (h_RT[0] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[2] -= q; /* CH2 */
    wdot[5] -= q; /* O */
    wdot[13] += q; /* HCO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[13];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[2] -= dqdci;                /* dwdot[CH2]/d[H] */
    J[5] -= dqdci;                /* dwdot[O]/d[H] */
    J[13] += dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[5];
    J[80] += dqdci;               /* dwdot[H]/d[CH2] */
    J[82] -= dqdci;               /* dwdot[CH2]/d[CH2] */
    J[85] -= dqdci;               /* dwdot[O]/d[CH2] */
    J[93] += dqdci;               /* dwdot[HCO]/d[CH2] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[2];
    J[200] += dqdci;              /* dwdot[H]/d[O] */
    J[202] -= dqdci;              /* dwdot[CH2]/d[O] */
    J[205] -= dqdci;              /* dwdot[O]/d[O] */
    J[213] += dqdci;              /* dwdot[HCO]/d[O] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[0];
    J[520] += dqdci;              /* dwdot[H]/d[HCO] */
    J[522] -= dqdci;              /* dwdot[CH2]/d[HCO] */
    J[525] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[533] += dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[1560] += dqdT;              /* dwdot[H]/dT */
    J[1562] -= dqdT;              /* dwdot[CH2]/dT */
    J[1565] -= dqdT;              /* dwdot[O]/dT */
    J[1573] += dqdT;              /* dwdot[HCO]/dT */

    /*reaction 117: CH2 + OH <=> CH2O + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[7];
    k_f = prefactor_units[116] * fwd_A[116]
                * exp(fwd_beta[116] * tc[0] - activation_units[116] * fwd_Ea[116] * invT);
    dlnkfdT = fwd_beta[116] * invT + activation_units[116] * fwd_Ea[116] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[15];
    Kc = exp(-g_RT[0] + g_RT[2] + g_RT[7] - g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[7]) + (h_RT[0] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[2] -= q; /* CH2 */
    wdot[7] -= q; /* OH */
    wdot[15] += q; /* CH2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[15];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[2] -= dqdci;                /* dwdot[CH2]/d[H] */
    J[7] -= dqdci;                /* dwdot[OH]/d[H] */
    J[15] += dqdci;               /* dwdot[CH2O]/d[H] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[7];
    J[80] += dqdci;               /* dwdot[H]/d[CH2] */
    J[82] -= dqdci;               /* dwdot[CH2]/d[CH2] */
    J[87] -= dqdci;               /* dwdot[OH]/d[CH2] */
    J[95] += dqdci;               /* dwdot[CH2O]/d[CH2] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[2];
    J[280] += dqdci;              /* dwdot[H]/d[OH] */
    J[282] -= dqdci;              /* dwdot[CH2]/d[OH] */
    J[287] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[295] += dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[0];
    J[600] += dqdci;              /* dwdot[H]/d[CH2O] */
    J[602] -= dqdci;              /* dwdot[CH2]/d[CH2O] */
    J[607] -= dqdci;              /* dwdot[OH]/d[CH2O] */
    J[615] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[1560] += dqdT;              /* dwdot[H]/dT */
    J[1562] -= dqdT;              /* dwdot[CH2]/dT */
    J[1567] -= dqdT;              /* dwdot[OH]/dT */
    J[1575] += dqdT;              /* dwdot[CH2O]/dT */

    /*reaction 118: CH2 + H2 <=> H + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[2];
    k_f = prefactor_units[117] * fwd_A[117]
                * exp(fwd_beta[117] * tc[0] - activation_units[117] * fwd_Ea[117] * invT);
    dlnkfdT = fwd_beta[117] * invT + activation_units[117] * fwd_Ea[117] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[4];
    Kc = exp(-g_RT[0] + g_RT[1] + g_RT[2] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[2]) + (h_RT[0] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[1] -= q; /* H2 */
    wdot[2] -= q; /* CH2 */
    wdot[4] += q; /* CH3 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[4];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[1] -= dqdci;                /* dwdot[H2]/d[H] */
    J[2] -= dqdci;                /* dwdot[CH2]/d[H] */
    J[4] += dqdci;                /* dwdot[CH3]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[2];
    J[40] += dqdci;               /* dwdot[H]/d[H2] */
    J[41] -= dqdci;               /* dwdot[H2]/d[H2] */
    J[42] -= dqdci;               /* dwdot[CH2]/d[H2] */
    J[44] += dqdci;               /* dwdot[CH3]/d[H2] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[1];
    J[80] += dqdci;               /* dwdot[H]/d[CH2] */
    J[81] -= dqdci;               /* dwdot[H2]/d[CH2] */
    J[82] -= dqdci;               /* dwdot[CH2]/d[CH2] */
    J[84] += dqdci;               /* dwdot[CH3]/d[CH2] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[0];
    J[160] += dqdci;              /* dwdot[H]/d[CH3] */
    J[161] -= dqdci;              /* dwdot[H2]/d[CH3] */
    J[162] -= dqdci;              /* dwdot[CH2]/d[CH3] */
    J[164] += dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[1560] += dqdT;              /* dwdot[H]/dT */
    J[1561] -= dqdT;              /* dwdot[H2]/dT */
    J[1562] -= dqdT;              /* dwdot[CH2]/dT */
    J[1564] += dqdT;              /* dwdot[CH3]/dT */

    /*reaction 119: CH2 + O2 <=> HCO + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[19];
    k_f = prefactor_units[118] * fwd_A[118]
                * exp(fwd_beta[118] * tc[0] - activation_units[118] * fwd_Ea[118] * invT);
    dlnkfdT = fwd_beta[118] * invT + activation_units[118] * fwd_Ea[118] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[13];
    Kc = exp(g_RT[2] - g_RT[7] - g_RT[13] + g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[19]) + (h_RT[7] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* CH2 */
    wdot[7] += q; /* OH */
    wdot[13] += q; /* HCO */
    wdot[19] -= q; /* O2 */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[19];
    J[82] -= dqdci;               /* dwdot[CH2]/d[CH2] */
    J[87] += dqdci;               /* dwdot[OH]/d[CH2] */
    J[93] += dqdci;               /* dwdot[HCO]/d[CH2] */
    J[99] -= dqdci;               /* dwdot[O2]/d[CH2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[13];
    J[282] -= dqdci;              /* dwdot[CH2]/d[OH] */
    J[287] += dqdci;              /* dwdot[OH]/d[OH] */
    J[293] += dqdci;              /* dwdot[HCO]/d[OH] */
    J[299] -= dqdci;              /* dwdot[O2]/d[OH] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[7];
    J[522] -= dqdci;              /* dwdot[CH2]/d[HCO] */
    J[527] += dqdci;              /* dwdot[OH]/d[HCO] */
    J[533] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[539] -= dqdci;              /* dwdot[O2]/d[HCO] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[2];
    J[762] -= dqdci;              /* dwdot[CH2]/d[O2] */
    J[767] += dqdci;              /* dwdot[OH]/d[O2] */
    J[773] += dqdci;              /* dwdot[HCO]/d[O2] */
    J[779] -= dqdci;              /* dwdot[O2]/d[O2] */
    /* d()/dT */
    J[1562] -= dqdT;              /* dwdot[CH2]/dT */
    J[1567] += dqdT;              /* dwdot[OH]/dT */
    J[1573] += dqdT;              /* dwdot[HCO]/dT */
    J[1579] -= dqdT;              /* dwdot[O2]/dT */

    /*reaction 120: CH2 + HO2 <=> CH2O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[20];
    k_f = prefactor_units[119] * fwd_A[119]
                * exp(fwd_beta[119] * tc[0] - activation_units[119] * fwd_Ea[119] * invT);
    dlnkfdT = fwd_beta[119] * invT + activation_units[119] * fwd_Ea[119] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[15];
    Kc = exp(g_RT[2] - g_RT[7] - g_RT[15] + g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[20]) + (h_RT[7] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* CH2 */
    wdot[7] += q; /* OH */
    wdot[15] += q; /* CH2O */
    wdot[20] -= q; /* HO2 */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[20];
    J[82] -= dqdci;               /* dwdot[CH2]/d[CH2] */
    J[87] += dqdci;               /* dwdot[OH]/d[CH2] */
    J[95] += dqdci;               /* dwdot[CH2O]/d[CH2] */
    J[100] -= dqdci;              /* dwdot[HO2]/d[CH2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[15];
    J[282] -= dqdci;              /* dwdot[CH2]/d[OH] */
    J[287] += dqdci;              /* dwdot[OH]/d[OH] */
    J[295] += dqdci;              /* dwdot[CH2O]/d[OH] */
    J[300] -= dqdci;              /* dwdot[HO2]/d[OH] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[7];
    J[602] -= dqdci;              /* dwdot[CH2]/d[CH2O] */
    J[607] += dqdci;              /* dwdot[OH]/d[CH2O] */
    J[615] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[620] -= dqdci;              /* dwdot[HO2]/d[CH2O] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[2];
    J[802] -= dqdci;              /* dwdot[CH2]/d[HO2] */
    J[807] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[815] += dqdci;              /* dwdot[CH2O]/d[HO2] */
    J[820] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1562] -= dqdT;              /* dwdot[CH2]/dT */
    J[1567] += dqdT;              /* dwdot[OH]/dT */
    J[1575] += dqdT;              /* dwdot[CH2O]/dT */
    J[1580] -= dqdT;              /* dwdot[HO2]/dT */

    /*reaction 121: CH2 + CH2 <=> C2H2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[2], 2.000000);
    k_f = prefactor_units[120] * fwd_A[120]
                * exp(fwd_beta[120] * tc[0] - activation_units[120] * fwd_Ea[120] * invT);
    dlnkfdT = fwd_beta[120] * invT + activation_units[120] * fwd_Ea[120] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[9];
    Kc = exp(-g_RT[1] + g_RT[2] + g_RT[2] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[2]) + (h_RT[1] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H2 */
    wdot[2] -= 2 * q; /* CH2 */
    wdot[9] += q; /* C2H2 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[9];
    J[41] += dqdci;               /* dwdot[H2]/d[H2] */
    J[42] += -2 * dqdci;          /* dwdot[CH2]/d[H2] */
    J[49] += dqdci;               /* dwdot[C2H2]/d[H2] */
    /* d()/d[CH2] */
    dqdci =  + k_f*2.000000*sc[2];
    J[81] += dqdci;               /* dwdot[H2]/d[CH2] */
    J[82] += -2 * dqdci;          /* dwdot[CH2]/d[CH2] */
    J[89] += dqdci;               /* dwdot[C2H2]/d[CH2] */
    /* d()/d[C2H2] */
    dqdci =  - k_r*sc[1];
    J[361] += dqdci;              /* dwdot[H2]/d[C2H2] */
    J[362] += -2 * dqdci;         /* dwdot[CH2]/d[C2H2] */
    J[369] += dqdci;              /* dwdot[C2H2]/d[C2H2] */
    /* d()/dT */
    J[1561] += dqdT;              /* dwdot[H2]/dT */
    J[1562] += -2 * dqdT;         /* dwdot[CH2]/dT */
    J[1569] += dqdT;              /* dwdot[C2H2]/dT */

    /*reaction 122: CH2(S) + H2O <=> CH2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[8];
    k_f = prefactor_units[121] * fwd_A[121]
                * exp(fwd_beta[121] * tc[0] - activation_units[121] * fwd_Ea[121] * invT);
    dlnkfdT = fwd_beta[121] * invT + activation_units[121] * fwd_Ea[121] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[8];
    Kc = exp(-g_RT[2] + g_RT[3] + g_RT[8] - g_RT[8]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[8]) + (h_RT[2] + h_RT[8]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* CH2 */
    wdot[3] -= q; /* CH2(S) */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[8];
    J[82] += dqdci;               /* dwdot[CH2]/d[CH2] */
    J[83] -= dqdci;               /* dwdot[CH2(S)]/d[CH2] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[8];
    J[122] += dqdci;              /* dwdot[CH2]/d[CH2(S)] */
    J[123] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[3] - k_r*sc[2];
    J[322] += dqdci;              /* dwdot[CH2]/d[H2O] */
    J[323] -= dqdci;              /* dwdot[CH2(S)]/d[H2O] */
    /* d()/dT */
    J[1562] += dqdT;              /* dwdot[CH2]/dT */
    J[1563] -= dqdT;              /* dwdot[CH2(S)]/dT */

    /*reaction 123: CH2(S) + CO <=> CH2 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[11];
    k_f = prefactor_units[122] * fwd_A[122]
                * exp(fwd_beta[122] * tc[0] - activation_units[122] * fwd_Ea[122] * invT);
    dlnkfdT = fwd_beta[122] * invT + activation_units[122] * fwd_Ea[122] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[11];
    Kc = exp(-g_RT[2] + g_RT[3] + g_RT[11] - g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[11]) + (h_RT[2] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* CH2 */
    wdot[3] -= q; /* CH2(S) */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[11];
    J[82] += dqdci;               /* dwdot[CH2]/d[CH2] */
    J[83] -= dqdci;               /* dwdot[CH2(S)]/d[CH2] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[11];
    J[122] += dqdci;              /* dwdot[CH2]/d[CH2(S)] */
    J[123] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[3] - k_r*sc[2];
    J[442] += dqdci;              /* dwdot[CH2]/d[CO] */
    J[443] -= dqdci;              /* dwdot[CH2(S)]/d[CO] */
    /* d()/dT */
    J[1562] += dqdT;              /* dwdot[CH2]/dT */
    J[1563] -= dqdT;              /* dwdot[CH2(S)]/dT */

    /*reaction 124: CH2(S) + CO2 <=> CH2 + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[22];
    k_f = prefactor_units[123] * fwd_A[123]
                * exp(fwd_beta[123] * tc[0] - activation_units[123] * fwd_Ea[123] * invT);
    dlnkfdT = fwd_beta[123] * invT + activation_units[123] * fwd_Ea[123] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[22];
    Kc = exp(-g_RT[2] + g_RT[3] + g_RT[22] - g_RT[22]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[22]) + (h_RT[2] + h_RT[22]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* CH2 */
    wdot[3] -= q; /* CH2(S) */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[22];
    J[82] += dqdci;               /* dwdot[CH2]/d[CH2] */
    J[83] -= dqdci;               /* dwdot[CH2(S)]/d[CH2] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[22];
    J[122] += dqdci;              /* dwdot[CH2]/d[CH2(S)] */
    J[123] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    /* d()/d[CO2] */
    dqdci =  + k_f*sc[3] - k_r*sc[2];
    J[882] += dqdci;              /* dwdot[CH2]/d[CO2] */
    J[883] -= dqdci;              /* dwdot[CH2(S)]/d[CO2] */
    /* d()/dT */
    J[1562] += dqdT;              /* dwdot[CH2]/dT */
    J[1563] -= dqdT;              /* dwdot[CH2(S)]/dT */

    /*reaction 125: CH2(S) + O <=> CO + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[5];
    k_f = prefactor_units[124] * fwd_A[124]
                * exp(fwd_beta[124] * tc[0] - activation_units[124] * fwd_Ea[124] * invT);
    dlnkfdT = fwd_beta[124] * invT + activation_units[124] * fwd_Ea[124] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[11];
    Kc = exp(-g_RT[1] + g_RT[3] + g_RT[5] - g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[5]) + (h_RT[1] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H2 */
    wdot[3] -= q; /* CH2(S) */
    wdot[5] -= q; /* O */
    wdot[11] += q; /* CO */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[11];
    J[41] += dqdci;               /* dwdot[H2]/d[H2] */
    J[43] -= dqdci;               /* dwdot[CH2(S)]/d[H2] */
    J[45] -= dqdci;               /* dwdot[O]/d[H2] */
    J[51] += dqdci;               /* dwdot[CO]/d[H2] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[5];
    J[121] += dqdci;              /* dwdot[H2]/d[CH2(S)] */
    J[123] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[125] -= dqdci;              /* dwdot[O]/d[CH2(S)] */
    J[131] += dqdci;              /* dwdot[CO]/d[CH2(S)] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[3];
    J[201] += dqdci;              /* dwdot[H2]/d[O] */
    J[203] -= dqdci;              /* dwdot[CH2(S)]/d[O] */
    J[205] -= dqdci;              /* dwdot[O]/d[O] */
    J[211] += dqdci;              /* dwdot[CO]/d[O] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[1];
    J[441] += dqdci;              /* dwdot[H2]/d[CO] */
    J[443] -= dqdci;              /* dwdot[CH2(S)]/d[CO] */
    J[445] -= dqdci;              /* dwdot[O]/d[CO] */
    J[451] += dqdci;              /* dwdot[CO]/d[CO] */
    /* d()/dT */
    J[1561] += dqdT;              /* dwdot[H2]/dT */
    J[1563] -= dqdT;              /* dwdot[CH2(S)]/dT */
    J[1565] -= dqdT;              /* dwdot[O]/dT */
    J[1571] += dqdT;              /* dwdot[CO]/dT */

    /*reaction 126: CH2(S) + O <=> HCO + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[5];
    k_f = prefactor_units[125] * fwd_A[125]
                * exp(fwd_beta[125] * tc[0] - activation_units[125] * fwd_Ea[125] * invT);
    dlnkfdT = fwd_beta[125] * invT + activation_units[125] * fwd_Ea[125] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[13];
    Kc = exp(-g_RT[0] + g_RT[3] + g_RT[5] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[5]) + (h_RT[0] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[3] -= q; /* CH2(S) */
    wdot[5] -= q; /* O */
    wdot[13] += q; /* HCO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[13];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[3] -= dqdci;                /* dwdot[CH2(S)]/d[H] */
    J[5] -= dqdci;                /* dwdot[O]/d[H] */
    J[13] += dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[5];
    J[120] += dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[123] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[125] -= dqdci;              /* dwdot[O]/d[CH2(S)] */
    J[133] += dqdci;              /* dwdot[HCO]/d[CH2(S)] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[3];
    J[200] += dqdci;              /* dwdot[H]/d[O] */
    J[203] -= dqdci;              /* dwdot[CH2(S)]/d[O] */
    J[205] -= dqdci;              /* dwdot[O]/d[O] */
    J[213] += dqdci;              /* dwdot[HCO]/d[O] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[0];
    J[520] += dqdci;              /* dwdot[H]/d[HCO] */
    J[523] -= dqdci;              /* dwdot[CH2(S)]/d[HCO] */
    J[525] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[533] += dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[1560] += dqdT;              /* dwdot[H]/dT */
    J[1563] -= dqdT;              /* dwdot[CH2(S)]/dT */
    J[1565] -= dqdT;              /* dwdot[O]/dT */
    J[1573] += dqdT;              /* dwdot[HCO]/dT */

    /*reaction 127: CH2(S) + OH <=> CH2O + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[7];
    k_f = prefactor_units[126] * fwd_A[126]
                * exp(fwd_beta[126] * tc[0] - activation_units[126] * fwd_Ea[126] * invT);
    dlnkfdT = fwd_beta[126] * invT + activation_units[126] * fwd_Ea[126] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[15];
    Kc = exp(-g_RT[0] + g_RT[3] + g_RT[7] - g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[7]) + (h_RT[0] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[3] -= q; /* CH2(S) */
    wdot[7] -= q; /* OH */
    wdot[15] += q; /* CH2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[15];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[3] -= dqdci;                /* dwdot[CH2(S)]/d[H] */
    J[7] -= dqdci;                /* dwdot[OH]/d[H] */
    J[15] += dqdci;               /* dwdot[CH2O]/d[H] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[7];
    J[120] += dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[123] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[127] -= dqdci;              /* dwdot[OH]/d[CH2(S)] */
    J[135] += dqdci;              /* dwdot[CH2O]/d[CH2(S)] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[3];
    J[280] += dqdci;              /* dwdot[H]/d[OH] */
    J[283] -= dqdci;              /* dwdot[CH2(S)]/d[OH] */
    J[287] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[295] += dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[0];
    J[600] += dqdci;              /* dwdot[H]/d[CH2O] */
    J[603] -= dqdci;              /* dwdot[CH2(S)]/d[CH2O] */
    J[607] -= dqdci;              /* dwdot[OH]/d[CH2O] */
    J[615] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[1560] += dqdT;              /* dwdot[H]/dT */
    J[1563] -= dqdT;              /* dwdot[CH2(S)]/dT */
    J[1567] -= dqdT;              /* dwdot[OH]/dT */
    J[1575] += dqdT;              /* dwdot[CH2O]/dT */

    /*reaction 128: CH2(S) + H2 <=> CH3 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = prefactor_units[127] * fwd_A[127]
                * exp(fwd_beta[127] * tc[0] - activation_units[127] * fwd_Ea[127] * invT);
    dlnkfdT = fwd_beta[127] * invT + activation_units[127] * fwd_Ea[127] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[4];
    Kc = exp(-g_RT[0] + g_RT[1] + g_RT[3] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3]) + (h_RT[0] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[1] -= q; /* H2 */
    wdot[3] -= q; /* CH2(S) */
    wdot[4] += q; /* CH3 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[4];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[1] -= dqdci;                /* dwdot[H2]/d[H] */
    J[3] -= dqdci;                /* dwdot[CH2(S)]/d[H] */
    J[4] += dqdci;                /* dwdot[CH3]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[3];
    J[40] += dqdci;               /* dwdot[H]/d[H2] */
    J[41] -= dqdci;               /* dwdot[H2]/d[H2] */
    J[43] -= dqdci;               /* dwdot[CH2(S)]/d[H2] */
    J[44] += dqdci;               /* dwdot[CH3]/d[H2] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[1];
    J[120] += dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[121] -= dqdci;              /* dwdot[H2]/d[CH2(S)] */
    J[123] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[124] += dqdci;              /* dwdot[CH3]/d[CH2(S)] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[0];
    J[160] += dqdci;              /* dwdot[H]/d[CH3] */
    J[161] -= dqdci;              /* dwdot[H2]/d[CH3] */
    J[163] -= dqdci;              /* dwdot[CH2(S)]/d[CH3] */
    J[164] += dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[1560] += dqdT;              /* dwdot[H]/dT */
    J[1561] -= dqdT;              /* dwdot[H2]/dT */
    J[1563] -= dqdT;              /* dwdot[CH2(S)]/dT */
    J[1564] += dqdT;              /* dwdot[CH3]/dT */

    /*reaction 129: CH2(S) + O2 <=> H + OH + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[19];
    k_f = prefactor_units[128] * fwd_A[128]
                * exp(fwd_beta[128] * tc[0] - activation_units[128] * fwd_Ea[128] * invT);
    dlnkfdT = fwd_beta[128] * invT + activation_units[128] * fwd_Ea[128] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[7]*sc[11];
    Kc = refC * exp(-g_RT[0] + g_RT[3] - g_RT[7] - g_RT[11] + g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[19]) + (h_RT[0] + h_RT[7] + h_RT[11]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[3] -= q; /* CH2(S) */
    wdot[7] += q; /* OH */
    wdot[11] += q; /* CO */
    wdot[19] -= q; /* O2 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[7]*sc[11];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[3] -= dqdci;                /* dwdot[CH2(S)]/d[H] */
    J[7] += dqdci;                /* dwdot[OH]/d[H] */
    J[11] += dqdci;               /* dwdot[CO]/d[H] */
    J[19] -= dqdci;               /* dwdot[O2]/d[H] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[19];
    J[120] += dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[123] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[127] += dqdci;              /* dwdot[OH]/d[CH2(S)] */
    J[131] += dqdci;              /* dwdot[CO]/d[CH2(S)] */
    J[139] -= dqdci;              /* dwdot[O2]/d[CH2(S)] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[0]*sc[11];
    J[280] += dqdci;              /* dwdot[H]/d[OH] */
    J[283] -= dqdci;              /* dwdot[CH2(S)]/d[OH] */
    J[287] += dqdci;              /* dwdot[OH]/d[OH] */
    J[291] += dqdci;              /* dwdot[CO]/d[OH] */
    J[299] -= dqdci;              /* dwdot[O2]/d[OH] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[0]*sc[7];
    J[440] += dqdci;              /* dwdot[H]/d[CO] */
    J[443] -= dqdci;              /* dwdot[CH2(S)]/d[CO] */
    J[447] += dqdci;              /* dwdot[OH]/d[CO] */
    J[451] += dqdci;              /* dwdot[CO]/d[CO] */
    J[459] -= dqdci;              /* dwdot[O2]/d[CO] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[3];
    J[760] += dqdci;              /* dwdot[H]/d[O2] */
    J[763] -= dqdci;              /* dwdot[CH2(S)]/d[O2] */
    J[767] += dqdci;              /* dwdot[OH]/d[O2] */
    J[771] += dqdci;              /* dwdot[CO]/d[O2] */
    J[779] -= dqdci;              /* dwdot[O2]/d[O2] */
    /* d()/dT */
    J[1560] += dqdT;              /* dwdot[H]/dT */
    J[1563] -= dqdT;              /* dwdot[CH2(S)]/dT */
    J[1567] += dqdT;              /* dwdot[OH]/dT */
    J[1571] += dqdT;              /* dwdot[CO]/dT */
    J[1579] -= dqdT;              /* dwdot[O2]/dT */

    /*reaction 130: CH2(S) + O2 <=> CO + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[19];
    k_f = prefactor_units[129] * fwd_A[129]
                * exp(fwd_beta[129] * tc[0] - activation_units[129] * fwd_Ea[129] * invT);
    dlnkfdT = fwd_beta[129] * invT + activation_units[129] * fwd_Ea[129] * invT2;
    /* reverse */
    phi_r = sc[8]*sc[11];
    Kc = exp(g_RT[3] - g_RT[8] - g_RT[11] + g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[19]) + (h_RT[8] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* CH2(S) */
    wdot[8] += q; /* H2O */
    wdot[11] += q; /* CO */
    wdot[19] -= q; /* O2 */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[19];
    J[123] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[128] += dqdci;              /* dwdot[H2O]/d[CH2(S)] */
    J[131] += dqdci;              /* dwdot[CO]/d[CH2(S)] */
    J[139] -= dqdci;              /* dwdot[O2]/d[CH2(S)] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[11];
    J[323] -= dqdci;              /* dwdot[CH2(S)]/d[H2O] */
    J[328] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[331] += dqdci;              /* dwdot[CO]/d[H2O] */
    J[339] -= dqdci;              /* dwdot[O2]/d[H2O] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[8];
    J[443] -= dqdci;              /* dwdot[CH2(S)]/d[CO] */
    J[448] += dqdci;              /* dwdot[H2O]/d[CO] */
    J[451] += dqdci;              /* dwdot[CO]/d[CO] */
    J[459] -= dqdci;              /* dwdot[O2]/d[CO] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[3];
    J[763] -= dqdci;              /* dwdot[CH2(S)]/d[O2] */
    J[768] += dqdci;              /* dwdot[H2O]/d[O2] */
    J[771] += dqdci;              /* dwdot[CO]/d[O2] */
    J[779] -= dqdci;              /* dwdot[O2]/d[O2] */
    /* d()/dT */
    J[1563] -= dqdT;              /* dwdot[CH2(S)]/dT */
    J[1568] += dqdT;              /* dwdot[H2O]/dT */
    J[1571] += dqdT;              /* dwdot[CO]/dT */
    J[1579] -= dqdT;              /* dwdot[O2]/dT */

    /*reaction 131: CH2(S) + CO2 <=> CH2O + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[22];
    k_f = prefactor_units[130] * fwd_A[130]
                * exp(fwd_beta[130] * tc[0] - activation_units[130] * fwd_Ea[130] * invT);
    dlnkfdT = fwd_beta[130] * invT + activation_units[130] * fwd_Ea[130] * invT2;
    /* reverse */
    phi_r = sc[11]*sc[15];
    Kc = exp(g_RT[3] - g_RT[11] - g_RT[15] + g_RT[22]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[22]) + (h_RT[11] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* CH2(S) */
    wdot[11] += q; /* CO */
    wdot[15] += q; /* CH2O */
    wdot[22] -= q; /* CO2 */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[22];
    J[123] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[131] += dqdci;              /* dwdot[CO]/d[CH2(S)] */
    J[135] += dqdci;              /* dwdot[CH2O]/d[CH2(S)] */
    J[142] -= dqdci;              /* dwdot[CO2]/d[CH2(S)] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[15];
    J[443] -= dqdci;              /* dwdot[CH2(S)]/d[CO] */
    J[451] += dqdci;              /* dwdot[CO]/d[CO] */
    J[455] += dqdci;              /* dwdot[CH2O]/d[CO] */
    J[462] -= dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[11];
    J[603] -= dqdci;              /* dwdot[CH2(S)]/d[CH2O] */
    J[611] += dqdci;              /* dwdot[CO]/d[CH2O] */
    J[615] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[622] -= dqdci;              /* dwdot[CO2]/d[CH2O] */
    /* d()/d[CO2] */
    dqdci =  + k_f*sc[3];
    J[883] -= dqdci;              /* dwdot[CH2(S)]/d[CO2] */
    J[891] += dqdci;              /* dwdot[CO]/d[CO2] */
    J[895] += dqdci;              /* dwdot[CH2O]/d[CO2] */
    J[902] -= dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[1563] -= dqdT;              /* dwdot[CH2(S)]/dT */
    J[1571] += dqdT;              /* dwdot[CO]/dT */
    J[1575] += dqdT;              /* dwdot[CH2O]/dT */
    J[1582] -= dqdT;              /* dwdot[CO2]/dT */

    /*reaction 132: CH3HCO <=> CH3 + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[23];
    k_f = prefactor_units[131] * fwd_A[131]
                * exp(fwd_beta[131] * tc[0] - activation_units[131] * fwd_Ea[131] * invT);
    dlnkfdT = fwd_beta[131] * invT + activation_units[131] * fwd_Ea[131] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[13];
    Kc = refC * exp(-g_RT[4] - g_RT[13] + g_RT[23]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[23]) + (h_RT[4] + h_RT[13]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] += q; /* CH3 */
    wdot[13] += q; /* HCO */
    wdot[23] -= q; /* CH3HCO */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[13];
    J[164] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[173] += dqdci;              /* dwdot[HCO]/d[CH3] */
    J[183] -= dqdci;              /* dwdot[CH3HCO]/d[CH3] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[4];
    J[524] += dqdci;              /* dwdot[CH3]/d[HCO] */
    J[533] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[543] -= dqdci;              /* dwdot[CH3HCO]/d[HCO] */
    /* d()/d[CH3HCO] */
    dqdci =  + k_f;
    J[924] += dqdci;              /* dwdot[CH3]/d[CH3HCO] */
    J[933] += dqdci;              /* dwdot[HCO]/d[CH3HCO] */
    J[943] -= dqdci;              /* dwdot[CH3HCO]/d[CH3HCO] */
    /* d()/dT */
    J[1564] += dqdT;              /* dwdot[CH3]/dT */
    J[1573] += dqdT;              /* dwdot[HCO]/dT */
    J[1583] -= dqdT;              /* dwdot[CH3HCO]/dT */

    /*reaction 133: CH3OCH3 <=> CH3 + CH3O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[26];
    k_f = prefactor_units[132] * fwd_A[132]
                * exp(fwd_beta[132] * tc[0] - activation_units[132] * fwd_Ea[132] * invT);
    dlnkfdT = fwd_beta[132] * invT + activation_units[132] * fwd_Ea[132] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[18];
    Kc = refC * exp(-g_RT[4] - g_RT[18] + g_RT[26]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[26]) + (h_RT[4] + h_RT[18]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] += q; /* CH3 */
    wdot[18] += q; /* CH3O */
    wdot[26] -= q; /* CH3OCH3 */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[18];
    J[164] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[178] += dqdci;              /* dwdot[CH3O]/d[CH3] */
    J[186] -= dqdci;              /* dwdot[CH3OCH3]/d[CH3] */
    /* d()/d[CH3O] */
    dqdci =  - k_r*sc[4];
    J[724] += dqdci;              /* dwdot[CH3]/d[CH3O] */
    J[738] += dqdci;              /* dwdot[CH3O]/d[CH3O] */
    J[746] -= dqdci;              /* dwdot[CH3OCH3]/d[CH3O] */
    /* d()/d[CH3OCH3] */
    dqdci =  + k_f;
    J[1044] += dqdci;             /* dwdot[CH3]/d[CH3OCH3] */
    J[1058] += dqdci;             /* dwdot[CH3O]/d[CH3OCH3] */
    J[1066] -= dqdci;             /* dwdot[CH3OCH3]/d[CH3OCH3] */
    /* d()/dT */
    J[1564] += dqdT;              /* dwdot[CH3]/dT */
    J[1578] += dqdT;              /* dwdot[CH3O]/dT */
    J[1586] -= dqdT;              /* dwdot[CH3OCH3]/dT */

    /*reaction 134: CH3OCH3 + OH <=> CH3OCH2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[7]*sc[26];
    k_f = prefactor_units[133] * fwd_A[133]
                * exp(fwd_beta[133] * tc[0] - activation_units[133] * fwd_Ea[133] * invT);
    dlnkfdT = fwd_beta[133] * invT + activation_units[133] * fwd_Ea[133] * invT2;
    /* reverse */
    phi_r = sc[8]*sc[24];
    Kc = exp(g_RT[7] - g_RT[8] - g_RT[24] + g_RT[26]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7] + h_RT[26]) + (h_RT[8] + h_RT[24]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] -= q; /* OH */
    wdot[8] += q; /* H2O */
    wdot[24] += q; /* CH3OCH2 */
    wdot[26] -= q; /* CH3OCH3 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[26];
    J[287] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[288] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[304] += dqdci;              /* dwdot[CH3OCH2]/d[OH] */
    J[306] -= dqdci;              /* dwdot[CH3OCH3]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[24];
    J[327] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[328] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[344] += dqdci;              /* dwdot[CH3OCH2]/d[H2O] */
    J[346] -= dqdci;              /* dwdot[CH3OCH3]/d[H2O] */
    /* d()/d[CH3OCH2] */
    dqdci =  - k_r*sc[8];
    J[967] -= dqdci;              /* dwdot[OH]/d[CH3OCH2] */
    J[968] += dqdci;              /* dwdot[H2O]/d[CH3OCH2] */
    J[984] += dqdci;              /* dwdot[CH3OCH2]/d[CH3OCH2] */
    J[986] -= dqdci;              /* dwdot[CH3OCH3]/d[CH3OCH2] */
    /* d()/d[CH3OCH3] */
    dqdci =  + k_f*sc[7];
    J[1047] -= dqdci;             /* dwdot[OH]/d[CH3OCH3] */
    J[1048] += dqdci;             /* dwdot[H2O]/d[CH3OCH3] */
    J[1064] += dqdci;             /* dwdot[CH3OCH2]/d[CH3OCH3] */
    J[1066] -= dqdci;             /* dwdot[CH3OCH3]/d[CH3OCH3] */
    /* d()/dT */
    J[1567] -= dqdT;              /* dwdot[OH]/dT */
    J[1568] += dqdT;              /* dwdot[H2O]/dT */
    J[1584] += dqdT;              /* dwdot[CH3OCH2]/dT */
    J[1586] -= dqdT;              /* dwdot[CH3OCH3]/dT */

    /*reaction 135: CH3OCH3 + H <=> CH3OCH2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[26];
    k_f = prefactor_units[134] * fwd_A[134]
                * exp(fwd_beta[134] * tc[0] - activation_units[134] * fwd_Ea[134] * invT);
    dlnkfdT = fwd_beta[134] * invT + activation_units[134] * fwd_Ea[134] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[24];
    Kc = exp(g_RT[0] - g_RT[1] - g_RT[24] + g_RT[26]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[26]) + (h_RT[1] + h_RT[24]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[1] += q; /* H2 */
    wdot[24] += q; /* CH3OCH2 */
    wdot[26] -= q; /* CH3OCH3 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[26];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[1] += dqdci;                /* dwdot[H2]/d[H] */
    J[24] += dqdci;               /* dwdot[CH3OCH2]/d[H] */
    J[26] -= dqdci;               /* dwdot[CH3OCH3]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[24];
    J[40] -= dqdci;               /* dwdot[H]/d[H2] */
    J[41] += dqdci;               /* dwdot[H2]/d[H2] */
    J[64] += dqdci;               /* dwdot[CH3OCH2]/d[H2] */
    J[66] -= dqdci;               /* dwdot[CH3OCH3]/d[H2] */
    /* d()/d[CH3OCH2] */
    dqdci =  - k_r*sc[1];
    J[960] -= dqdci;              /* dwdot[H]/d[CH3OCH2] */
    J[961] += dqdci;              /* dwdot[H2]/d[CH3OCH2] */
    J[984] += dqdci;              /* dwdot[CH3OCH2]/d[CH3OCH2] */
    J[986] -= dqdci;              /* dwdot[CH3OCH3]/d[CH3OCH2] */
    /* d()/d[CH3OCH3] */
    dqdci =  + k_f*sc[0];
    J[1040] -= dqdci;             /* dwdot[H]/d[CH3OCH3] */
    J[1041] += dqdci;             /* dwdot[H2]/d[CH3OCH3] */
    J[1064] += dqdci;             /* dwdot[CH3OCH2]/d[CH3OCH3] */
    J[1066] -= dqdci;             /* dwdot[CH3OCH3]/d[CH3OCH3] */
    /* d()/dT */
    J[1560] -= dqdT;              /* dwdot[H]/dT */
    J[1561] += dqdT;              /* dwdot[H2]/dT */
    J[1584] += dqdT;              /* dwdot[CH3OCH2]/dT */
    J[1586] -= dqdT;              /* dwdot[CH3OCH3]/dT */

    /*reaction 136: CH3OCH3 + CH3 <=> CH3OCH2 + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[26];
    k_f = prefactor_units[135] * fwd_A[135]
                * exp(fwd_beta[135] * tc[0] - activation_units[135] * fwd_Ea[135] * invT);
    dlnkfdT = fwd_beta[135] * invT + activation_units[135] * fwd_Ea[135] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[24];
    Kc = exp(g_RT[4] - g_RT[6] - g_RT[24] + g_RT[26]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[26]) + (h_RT[6] + h_RT[24]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* CH3 */
    wdot[6] += q; /* CH4 */
    wdot[24] += q; /* CH3OCH2 */
    wdot[26] -= q; /* CH3OCH3 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[26];
    J[164] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[166] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[184] += dqdci;              /* dwdot[CH3OCH2]/d[CH3] */
    J[186] -= dqdci;              /* dwdot[CH3OCH3]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[24];
    J[244] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[246] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[264] += dqdci;              /* dwdot[CH3OCH2]/d[CH4] */
    J[266] -= dqdci;              /* dwdot[CH3OCH3]/d[CH4] */
    /* d()/d[CH3OCH2] */
    dqdci =  - k_r*sc[6];
    J[964] -= dqdci;              /* dwdot[CH3]/d[CH3OCH2] */
    J[966] += dqdci;              /* dwdot[CH4]/d[CH3OCH2] */
    J[984] += dqdci;              /* dwdot[CH3OCH2]/d[CH3OCH2] */
    J[986] -= dqdci;              /* dwdot[CH3OCH3]/d[CH3OCH2] */
    /* d()/d[CH3OCH3] */
    dqdci =  + k_f*sc[4];
    J[1044] -= dqdci;             /* dwdot[CH3]/d[CH3OCH3] */
    J[1046] += dqdci;             /* dwdot[CH4]/d[CH3OCH3] */
    J[1064] += dqdci;             /* dwdot[CH3OCH2]/d[CH3OCH3] */
    J[1066] -= dqdci;             /* dwdot[CH3OCH3]/d[CH3OCH3] */
    /* d()/dT */
    J[1564] -= dqdT;              /* dwdot[CH3]/dT */
    J[1566] += dqdT;              /* dwdot[CH4]/dT */
    J[1584] += dqdT;              /* dwdot[CH3OCH2]/dT */
    J[1586] -= dqdT;              /* dwdot[CH3OCH3]/dT */

    /*reaction 137: CH3OCH3 + O <=> CH3OCH2 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[26];
    k_f = prefactor_units[136] * fwd_A[136]
                * exp(fwd_beta[136] * tc[0] - activation_units[136] * fwd_Ea[136] * invT);
    dlnkfdT = fwd_beta[136] * invT + activation_units[136] * fwd_Ea[136] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[24];
    Kc = exp(g_RT[5] - g_RT[7] - g_RT[24] + g_RT[26]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[26]) + (h_RT[7] + h_RT[24]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] -= q; /* O */
    wdot[7] += q; /* OH */
    wdot[24] += q; /* CH3OCH2 */
    wdot[26] -= q; /* CH3OCH3 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[26];
    J[205] -= dqdci;              /* dwdot[O]/d[O] */
    J[207] += dqdci;              /* dwdot[OH]/d[O] */
    J[224] += dqdci;              /* dwdot[CH3OCH2]/d[O] */
    J[226] -= dqdci;              /* dwdot[CH3OCH3]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[24];
    J[285] -= dqdci;              /* dwdot[O]/d[OH] */
    J[287] += dqdci;              /* dwdot[OH]/d[OH] */
    J[304] += dqdci;              /* dwdot[CH3OCH2]/d[OH] */
    J[306] -= dqdci;              /* dwdot[CH3OCH3]/d[OH] */
    /* d()/d[CH3OCH2] */
    dqdci =  - k_r*sc[7];
    J[965] -= dqdci;              /* dwdot[O]/d[CH3OCH2] */
    J[967] += dqdci;              /* dwdot[OH]/d[CH3OCH2] */
    J[984] += dqdci;              /* dwdot[CH3OCH2]/d[CH3OCH2] */
    J[986] -= dqdci;              /* dwdot[CH3OCH3]/d[CH3OCH2] */
    /* d()/d[CH3OCH3] */
    dqdci =  + k_f*sc[5];
    J[1045] -= dqdci;             /* dwdot[O]/d[CH3OCH3] */
    J[1047] += dqdci;             /* dwdot[OH]/d[CH3OCH3] */
    J[1064] += dqdci;             /* dwdot[CH3OCH2]/d[CH3OCH3] */
    J[1066] -= dqdci;             /* dwdot[CH3OCH3]/d[CH3OCH3] */
    /* d()/dT */
    J[1565] -= dqdT;              /* dwdot[O]/dT */
    J[1567] += dqdT;              /* dwdot[OH]/dT */
    J[1584] += dqdT;              /* dwdot[CH3OCH2]/dT */
    J[1586] -= dqdT;              /* dwdot[CH3OCH3]/dT */

    /*reaction 138: CH3OCH3 + HO2 <=> CH3OCH2 + H2O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[20]*sc[26];
    k_f = prefactor_units[137] * fwd_A[137]
                * exp(fwd_beta[137] * tc[0] - activation_units[137] * fwd_Ea[137] * invT);
    dlnkfdT = fwd_beta[137] * invT + activation_units[137] * fwd_Ea[137] * invT2;
    /* reverse */
    phi_r = sc[21]*sc[24];
    Kc = exp(g_RT[20] - g_RT[21] - g_RT[24] + g_RT[26]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[20] + h_RT[26]) + (h_RT[21] + h_RT[24]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[20] -= q; /* HO2 */
    wdot[21] += q; /* H2O2 */
    wdot[24] += q; /* CH3OCH2 */
    wdot[26] -= q; /* CH3OCH3 */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[26];
    J[820] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[821] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    J[824] += dqdci;              /* dwdot[CH3OCH2]/d[HO2] */
    J[826] -= dqdci;              /* dwdot[CH3OCH3]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[24];
    J[860] -= dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[861] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
    J[864] += dqdci;              /* dwdot[CH3OCH2]/d[H2O2] */
    J[866] -= dqdci;              /* dwdot[CH3OCH3]/d[H2O2] */
    /* d()/d[CH3OCH2] */
    dqdci =  - k_r*sc[21];
    J[980] -= dqdci;              /* dwdot[HO2]/d[CH3OCH2] */
    J[981] += dqdci;              /* dwdot[H2O2]/d[CH3OCH2] */
    J[984] += dqdci;              /* dwdot[CH3OCH2]/d[CH3OCH2] */
    J[986] -= dqdci;              /* dwdot[CH3OCH3]/d[CH3OCH2] */
    /* d()/d[CH3OCH3] */
    dqdci =  + k_f*sc[20];
    J[1060] -= dqdci;             /* dwdot[HO2]/d[CH3OCH3] */
    J[1061] += dqdci;             /* dwdot[H2O2]/d[CH3OCH3] */
    J[1064] += dqdci;             /* dwdot[CH3OCH2]/d[CH3OCH3] */
    J[1066] -= dqdci;             /* dwdot[CH3OCH3]/d[CH3OCH3] */
    /* d()/dT */
    J[1580] -= dqdT;              /* dwdot[HO2]/dT */
    J[1581] += dqdT;              /* dwdot[H2O2]/dT */
    J[1584] += dqdT;              /* dwdot[CH3OCH2]/dT */
    J[1586] -= dqdT;              /* dwdot[CH3OCH3]/dT */

    /*reaction 139: CH3OCH3 + O2 <=> CH3OCH2 + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[19]*sc[26];
    k_f = prefactor_units[138] * fwd_A[138]
                * exp(fwd_beta[138] * tc[0] - activation_units[138] * fwd_Ea[138] * invT);
    dlnkfdT = fwd_beta[138] * invT + activation_units[138] * fwd_Ea[138] * invT2;
    /* reverse */
    phi_r = sc[20]*sc[24];
    Kc = exp(g_RT[19] - g_RT[20] - g_RT[24] + g_RT[26]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[19] + h_RT[26]) + (h_RT[20] + h_RT[24]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[19] -= q; /* O2 */
    wdot[20] += q; /* HO2 */
    wdot[24] += q; /* CH3OCH2 */
    wdot[26] -= q; /* CH3OCH3 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[26];
    J[779] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[780] += dqdci;              /* dwdot[HO2]/d[O2] */
    J[784] += dqdci;              /* dwdot[CH3OCH2]/d[O2] */
    J[786] -= dqdci;              /* dwdot[CH3OCH3]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[24];
    J[819] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[820] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[824] += dqdci;              /* dwdot[CH3OCH2]/d[HO2] */
    J[826] -= dqdci;              /* dwdot[CH3OCH3]/d[HO2] */
    /* d()/d[CH3OCH2] */
    dqdci =  - k_r*sc[20];
    J[979] -= dqdci;              /* dwdot[O2]/d[CH3OCH2] */
    J[980] += dqdci;              /* dwdot[HO2]/d[CH3OCH2] */
    J[984] += dqdci;              /* dwdot[CH3OCH2]/d[CH3OCH2] */
    J[986] -= dqdci;              /* dwdot[CH3OCH3]/d[CH3OCH2] */
    /* d()/d[CH3OCH3] */
    dqdci =  + k_f*sc[19];
    J[1059] -= dqdci;             /* dwdot[O2]/d[CH3OCH3] */
    J[1060] += dqdci;             /* dwdot[HO2]/d[CH3OCH3] */
    J[1064] += dqdci;             /* dwdot[CH3OCH2]/d[CH3OCH3] */
    J[1066] -= dqdci;             /* dwdot[CH3OCH3]/d[CH3OCH3] */
    /* d()/dT */
    J[1579] -= dqdT;              /* dwdot[O2]/dT */
    J[1580] += dqdT;              /* dwdot[HO2]/dT */
    J[1584] += dqdT;              /* dwdot[CH3OCH2]/dT */
    J[1586] -= dqdT;              /* dwdot[CH3OCH3]/dT */

    /*reaction 140: CH3OCH2 <=> CH2O + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[24];
    k_f = prefactor_units[139] * fwd_A[139]
                * exp(fwd_beta[139] * tc[0] - activation_units[139] * fwd_Ea[139] * invT);
    dlnkfdT = fwd_beta[139] * invT + activation_units[139] * fwd_Ea[139] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[15];
    Kc = refC * exp(-g_RT[4] - g_RT[15] + g_RT[24]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[24]) + (h_RT[4] + h_RT[15]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] += q; /* CH3 */
    wdot[15] += q; /* CH2O */
    wdot[24] -= q; /* CH3OCH2 */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[15];
    J[164] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[175] += dqdci;              /* dwdot[CH2O]/d[CH3] */
    J[184] -= dqdci;              /* dwdot[CH3OCH2]/d[CH3] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[4];
    J[604] += dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[615] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[624] -= dqdci;              /* dwdot[CH3OCH2]/d[CH2O] */
    /* d()/d[CH3OCH2] */
    dqdci =  + k_f;
    J[964] += dqdci;              /* dwdot[CH3]/d[CH3OCH2] */
    J[975] += dqdci;              /* dwdot[CH2O]/d[CH3OCH2] */
    J[984] -= dqdci;              /* dwdot[CH3OCH2]/d[CH3OCH2] */
    /* d()/dT */
    J[1564] += dqdT;              /* dwdot[CH3]/dT */
    J[1575] += dqdT;              /* dwdot[CH2O]/dT */
    J[1584] -= dqdT;              /* dwdot[CH3OCH2]/dT */

    /*reaction 141: CH3OCH2 + CH3O <=> CH3OCH3 + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[18]*sc[24];
    k_f = prefactor_units[140] * fwd_A[140]
                * exp(fwd_beta[140] * tc[0] - activation_units[140] * fwd_Ea[140] * invT);
    dlnkfdT = fwd_beta[140] * invT + activation_units[140] * fwd_Ea[140] * invT2;
    /* reverse */
    phi_r = sc[15]*sc[26];
    Kc = exp(-g_RT[15] + g_RT[18] + g_RT[24] - g_RT[26]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[18] + h_RT[24]) + (h_RT[15] + h_RT[26]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[15] += q; /* CH2O */
    wdot[18] -= q; /* CH3O */
    wdot[24] -= q; /* CH3OCH2 */
    wdot[26] += q; /* CH3OCH3 */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[26];
    J[615] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[618] -= dqdci;              /* dwdot[CH3O]/d[CH2O] */
    J[624] -= dqdci;              /* dwdot[CH3OCH2]/d[CH2O] */
    J[626] += dqdci;              /* dwdot[CH3OCH3]/d[CH2O] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[24];
    J[735] += dqdci;              /* dwdot[CH2O]/d[CH3O] */
    J[738] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    J[744] -= dqdci;              /* dwdot[CH3OCH2]/d[CH3O] */
    J[746] += dqdci;              /* dwdot[CH3OCH3]/d[CH3O] */
    /* d()/d[CH3OCH2] */
    dqdci =  + k_f*sc[18];
    J[975] += dqdci;              /* dwdot[CH2O]/d[CH3OCH2] */
    J[978] -= dqdci;              /* dwdot[CH3O]/d[CH3OCH2] */
    J[984] -= dqdci;              /* dwdot[CH3OCH2]/d[CH3OCH2] */
    J[986] += dqdci;              /* dwdot[CH3OCH3]/d[CH3OCH2] */
    /* d()/d[CH3OCH3] */
    dqdci =  - k_r*sc[15];
    J[1055] += dqdci;             /* dwdot[CH2O]/d[CH3OCH3] */
    J[1058] -= dqdci;             /* dwdot[CH3O]/d[CH3OCH3] */
    J[1064] -= dqdci;             /* dwdot[CH3OCH2]/d[CH3OCH3] */
    J[1066] += dqdci;             /* dwdot[CH3OCH3]/d[CH3OCH3] */
    /* d()/dT */
    J[1575] += dqdT;              /* dwdot[CH2O]/dT */
    J[1578] -= dqdT;              /* dwdot[CH3O]/dT */
    J[1584] -= dqdT;              /* dwdot[CH3OCH2]/dT */
    J[1586] += dqdT;              /* dwdot[CH3OCH3]/dT */

    /*reaction 142: CH3OCH2 + CH2O <=> CH3OCH3 + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[15]*sc[24];
    k_f = prefactor_units[141] * fwd_A[141]
                * exp(fwd_beta[141] * tc[0] - activation_units[141] * fwd_Ea[141] * invT);
    dlnkfdT = fwd_beta[141] * invT + activation_units[141] * fwd_Ea[141] * invT2;
    /* reverse */
    phi_r = sc[13]*sc[26];
    Kc = exp(-g_RT[13] + g_RT[15] + g_RT[24] - g_RT[26]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[15] + h_RT[24]) + (h_RT[13] + h_RT[26]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[13] += q; /* HCO */
    wdot[15] -= q; /* CH2O */
    wdot[24] -= q; /* CH3OCH2 */
    wdot[26] += q; /* CH3OCH3 */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[26];
    J[533] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[535] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    J[544] -= dqdci;              /* dwdot[CH3OCH2]/d[HCO] */
    J[546] += dqdci;              /* dwdot[CH3OCH3]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[24];
    J[613] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[615] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[624] -= dqdci;              /* dwdot[CH3OCH2]/d[CH2O] */
    J[626] += dqdci;              /* dwdot[CH3OCH3]/d[CH2O] */
    /* d()/d[CH3OCH2] */
    dqdci =  + k_f*sc[15];
    J[973] += dqdci;              /* dwdot[HCO]/d[CH3OCH2] */
    J[975] -= dqdci;              /* dwdot[CH2O]/d[CH3OCH2] */
    J[984] -= dqdci;              /* dwdot[CH3OCH2]/d[CH3OCH2] */
    J[986] += dqdci;              /* dwdot[CH3OCH3]/d[CH3OCH2] */
    /* d()/d[CH3OCH3] */
    dqdci =  - k_r*sc[13];
    J[1053] += dqdci;             /* dwdot[HCO]/d[CH3OCH3] */
    J[1055] -= dqdci;             /* dwdot[CH2O]/d[CH3OCH3] */
    J[1064] -= dqdci;             /* dwdot[CH3OCH2]/d[CH3OCH3] */
    J[1066] += dqdci;             /* dwdot[CH3OCH3]/d[CH3OCH3] */
    /* d()/dT */
    J[1573] += dqdT;              /* dwdot[HCO]/dT */
    J[1575] -= dqdT;              /* dwdot[CH2O]/dT */
    J[1584] -= dqdT;              /* dwdot[CH3OCH2]/dT */
    J[1586] += dqdT;              /* dwdot[CH3OCH3]/dT */

    /*reaction 143: CH3OCH2 + HO2 <=> CH3OCH2O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[20]*sc[24];
    k_f = prefactor_units[142] * fwd_A[142]
                * exp(fwd_beta[142] * tc[0] - activation_units[142] * fwd_Ea[142] * invT);
    dlnkfdT = fwd_beta[142] * invT + activation_units[142] * fwd_Ea[142] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[30];
    Kc = exp(-g_RT[7] + g_RT[20] + g_RT[24] - g_RT[30]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[20] + h_RT[24]) + (h_RT[7] + h_RT[30]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] += q; /* OH */
    wdot[20] -= q; /* HO2 */
    wdot[24] -= q; /* CH3OCH2 */
    wdot[30] += q; /* CH3OCH2O */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[30];
    J[287] += dqdci;              /* dwdot[OH]/d[OH] */
    J[300] -= dqdci;              /* dwdot[HO2]/d[OH] */
    J[304] -= dqdci;              /* dwdot[CH3OCH2]/d[OH] */
    J[310] += dqdci;              /* dwdot[CH3OCH2O]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[24];
    J[807] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[820] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[824] -= dqdci;              /* dwdot[CH3OCH2]/d[HO2] */
    J[830] += dqdci;              /* dwdot[CH3OCH2O]/d[HO2] */
    /* d()/d[CH3OCH2] */
    dqdci =  + k_f*sc[20];
    J[967] += dqdci;              /* dwdot[OH]/d[CH3OCH2] */
    J[980] -= dqdci;              /* dwdot[HO2]/d[CH3OCH2] */
    J[984] -= dqdci;              /* dwdot[CH3OCH2]/d[CH3OCH2] */
    J[990] += dqdci;              /* dwdot[CH3OCH2O]/d[CH3OCH2] */
    /* d()/d[CH3OCH2O] */
    dqdci =  - k_r*sc[7];
    J[1207] += dqdci;             /* dwdot[OH]/d[CH3OCH2O] */
    J[1220] -= dqdci;             /* dwdot[HO2]/d[CH3OCH2O] */
    J[1224] -= dqdci;             /* dwdot[CH3OCH2]/d[CH3OCH2O] */
    J[1230] += dqdci;             /* dwdot[CH3OCH2O]/d[CH3OCH2O] */
    /* d()/dT */
    J[1567] += dqdT;              /* dwdot[OH]/dT */
    J[1580] -= dqdT;              /* dwdot[HO2]/dT */
    J[1584] -= dqdT;              /* dwdot[CH3OCH2]/dT */
    J[1590] += dqdT;              /* dwdot[CH3OCH2O]/dT */

    /*reaction 144: CH3OCH2O <=> CH3OCHO + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[30];
    k_f = prefactor_units[143] * fwd_A[143]
                * exp(fwd_beta[143] * tc[0] - activation_units[143] * fwd_Ea[143] * invT);
    dlnkfdT = fwd_beta[143] * invT + activation_units[143] * fwd_Ea[143] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[29];
    Kc = refC * exp(-g_RT[0] - g_RT[29] + g_RT[30]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[30]) + (h_RT[0] + h_RT[29]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[29] += q; /* CH3OCHO */
    wdot[30] -= q; /* CH3OCH2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[29];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[29] += dqdci;               /* dwdot[CH3OCHO]/d[H] */
    J[30] -= dqdci;               /* dwdot[CH3OCH2O]/d[H] */
    /* d()/d[CH3OCHO] */
    dqdci =  - k_r*sc[0];
    J[1160] += dqdci;             /* dwdot[H]/d[CH3OCHO] */
    J[1189] += dqdci;             /* dwdot[CH3OCHO]/d[CH3OCHO] */
    J[1190] -= dqdci;             /* dwdot[CH3OCH2O]/d[CH3OCHO] */
    /* d()/d[CH3OCH2O] */
    dqdci =  + k_f;
    J[1200] += dqdci;             /* dwdot[H]/d[CH3OCH2O] */
    J[1229] += dqdci;             /* dwdot[CH3OCHO]/d[CH3OCH2O] */
    J[1230] -= dqdci;             /* dwdot[CH3OCH2O]/d[CH3OCH2O] */
    /* d()/dT */
    J[1560] += dqdT;              /* dwdot[H]/dT */
    J[1589] += dqdT;              /* dwdot[CH3OCHO]/dT */
    J[1590] -= dqdT;              /* dwdot[CH3OCH2O]/dT */

    /*reaction 145: CH3OCHO + O2 <=> CH3OCO + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[19]*sc[29];
    k_f = prefactor_units[144] * fwd_A[144]
                * exp(fwd_beta[144] * tc[0] - activation_units[144] * fwd_Ea[144] * invT);
    dlnkfdT = fwd_beta[144] * invT + activation_units[144] * fwd_Ea[144] * invT2;
    /* reverse */
    phi_r = sc[20]*sc[28];
    Kc = exp(g_RT[19] - g_RT[20] - g_RT[28] + g_RT[29]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[19] + h_RT[29]) + (h_RT[20] + h_RT[28]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[19] -= q; /* O2 */
    wdot[20] += q; /* HO2 */
    wdot[28] += q; /* CH3OCO */
    wdot[29] -= q; /* CH3OCHO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[29];
    J[779] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[780] += dqdci;              /* dwdot[HO2]/d[O2] */
    J[788] += dqdci;              /* dwdot[CH3OCO]/d[O2] */
    J[789] -= dqdci;              /* dwdot[CH3OCHO]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[28];
    J[819] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[820] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[828] += dqdci;              /* dwdot[CH3OCO]/d[HO2] */
    J[829] -= dqdci;              /* dwdot[CH3OCHO]/d[HO2] */
    /* d()/d[CH3OCO] */
    dqdci =  - k_r*sc[20];
    J[1139] -= dqdci;             /* dwdot[O2]/d[CH3OCO] */
    J[1140] += dqdci;             /* dwdot[HO2]/d[CH3OCO] */
    J[1148] += dqdci;             /* dwdot[CH3OCO]/d[CH3OCO] */
    J[1149] -= dqdci;             /* dwdot[CH3OCHO]/d[CH3OCO] */
    /* d()/d[CH3OCHO] */
    dqdci =  + k_f*sc[19];
    J[1179] -= dqdci;             /* dwdot[O2]/d[CH3OCHO] */
    J[1180] += dqdci;             /* dwdot[HO2]/d[CH3OCHO] */
    J[1188] += dqdci;             /* dwdot[CH3OCO]/d[CH3OCHO] */
    J[1189] -= dqdci;             /* dwdot[CH3OCHO]/d[CH3OCHO] */
    /* d()/dT */
    J[1579] -= dqdT;              /* dwdot[O2]/dT */
    J[1580] += dqdT;              /* dwdot[HO2]/dT */
    J[1588] += dqdT;              /* dwdot[CH3OCO]/dT */
    J[1589] -= dqdT;              /* dwdot[CH3OCHO]/dT */

    /*reaction 146: CH3OCHO + OH <=> CH3OCO + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[7]*sc[29];
    k_f = prefactor_units[145] * fwd_A[145]
                * exp(fwd_beta[145] * tc[0] - activation_units[145] * fwd_Ea[145] * invT);
    dlnkfdT = fwd_beta[145] * invT + activation_units[145] * fwd_Ea[145] * invT2;
    /* reverse */
    phi_r = sc[8]*sc[28];
    Kc = exp(g_RT[7] - g_RT[8] - g_RT[28] + g_RT[29]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7] + h_RT[29]) + (h_RT[8] + h_RT[28]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] -= q; /* OH */
    wdot[8] += q; /* H2O */
    wdot[28] += q; /* CH3OCO */
    wdot[29] -= q; /* CH3OCHO */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[29];
    J[287] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[288] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[308] += dqdci;              /* dwdot[CH3OCO]/d[OH] */
    J[309] -= dqdci;              /* dwdot[CH3OCHO]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[28];
    J[327] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[328] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[348] += dqdci;              /* dwdot[CH3OCO]/d[H2O] */
    J[349] -= dqdci;              /* dwdot[CH3OCHO]/d[H2O] */
    /* d()/d[CH3OCO] */
    dqdci =  - k_r*sc[8];
    J[1127] -= dqdci;             /* dwdot[OH]/d[CH3OCO] */
    J[1128] += dqdci;             /* dwdot[H2O]/d[CH3OCO] */
    J[1148] += dqdci;             /* dwdot[CH3OCO]/d[CH3OCO] */
    J[1149] -= dqdci;             /* dwdot[CH3OCHO]/d[CH3OCO] */
    /* d()/d[CH3OCHO] */
    dqdci =  + k_f*sc[7];
    J[1167] -= dqdci;             /* dwdot[OH]/d[CH3OCHO] */
    J[1168] += dqdci;             /* dwdot[H2O]/d[CH3OCHO] */
    J[1188] += dqdci;             /* dwdot[CH3OCO]/d[CH3OCHO] */
    J[1189] -= dqdci;             /* dwdot[CH3OCHO]/d[CH3OCHO] */
    /* d()/dT */
    J[1567] -= dqdT;              /* dwdot[OH]/dT */
    J[1568] += dqdT;              /* dwdot[H2O]/dT */
    J[1588] += dqdT;              /* dwdot[CH3OCO]/dT */
    J[1589] -= dqdT;              /* dwdot[CH3OCHO]/dT */

    /*reaction 147: CH3OCHO + HO2 <=> CH3OCO + H2O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[20]*sc[29];
    k_f = prefactor_units[146] * fwd_A[146]
                * exp(fwd_beta[146] * tc[0] - activation_units[146] * fwd_Ea[146] * invT);
    dlnkfdT = fwd_beta[146] * invT + activation_units[146] * fwd_Ea[146] * invT2;
    /* reverse */
    phi_r = sc[21]*sc[28];
    Kc = exp(g_RT[20] - g_RT[21] - g_RT[28] + g_RT[29]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[20] + h_RT[29]) + (h_RT[21] + h_RT[28]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[20] -= q; /* HO2 */
    wdot[21] += q; /* H2O2 */
    wdot[28] += q; /* CH3OCO */
    wdot[29] -= q; /* CH3OCHO */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[29];
    J[820] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[821] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    J[828] += dqdci;              /* dwdot[CH3OCO]/d[HO2] */
    J[829] -= dqdci;              /* dwdot[CH3OCHO]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[28];
    J[860] -= dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[861] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
    J[868] += dqdci;              /* dwdot[CH3OCO]/d[H2O2] */
    J[869] -= dqdci;              /* dwdot[CH3OCHO]/d[H2O2] */
    /* d()/d[CH3OCO] */
    dqdci =  - k_r*sc[21];
    J[1140] -= dqdci;             /* dwdot[HO2]/d[CH3OCO] */
    J[1141] += dqdci;             /* dwdot[H2O2]/d[CH3OCO] */
    J[1148] += dqdci;             /* dwdot[CH3OCO]/d[CH3OCO] */
    J[1149] -= dqdci;             /* dwdot[CH3OCHO]/d[CH3OCO] */
    /* d()/d[CH3OCHO] */
    dqdci =  + k_f*sc[20];
    J[1180] -= dqdci;             /* dwdot[HO2]/d[CH3OCHO] */
    J[1181] += dqdci;             /* dwdot[H2O2]/d[CH3OCHO] */
    J[1188] += dqdci;             /* dwdot[CH3OCO]/d[CH3OCHO] */
    J[1189] -= dqdci;             /* dwdot[CH3OCHO]/d[CH3OCHO] */
    /* d()/dT */
    J[1580] -= dqdT;              /* dwdot[HO2]/dT */
    J[1581] += dqdT;              /* dwdot[H2O2]/dT */
    J[1588] += dqdT;              /* dwdot[CH3OCO]/dT */
    J[1589] -= dqdT;              /* dwdot[CH3OCHO]/dT */

    /*reaction 148: CH3OCHO + O <=> CH3OCO + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[29];
    k_f = prefactor_units[147] * fwd_A[147]
                * exp(fwd_beta[147] * tc[0] - activation_units[147] * fwd_Ea[147] * invT);
    dlnkfdT = fwd_beta[147] * invT + activation_units[147] * fwd_Ea[147] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[28];
    Kc = exp(g_RT[5] - g_RT[7] - g_RT[28] + g_RT[29]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[29]) + (h_RT[7] + h_RT[28]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] -= q; /* O */
    wdot[7] += q; /* OH */
    wdot[28] += q; /* CH3OCO */
    wdot[29] -= q; /* CH3OCHO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[29];
    J[205] -= dqdci;              /* dwdot[O]/d[O] */
    J[207] += dqdci;              /* dwdot[OH]/d[O] */
    J[228] += dqdci;              /* dwdot[CH3OCO]/d[O] */
    J[229] -= dqdci;              /* dwdot[CH3OCHO]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[28];
    J[285] -= dqdci;              /* dwdot[O]/d[OH] */
    J[287] += dqdci;              /* dwdot[OH]/d[OH] */
    J[308] += dqdci;              /* dwdot[CH3OCO]/d[OH] */
    J[309] -= dqdci;              /* dwdot[CH3OCHO]/d[OH] */
    /* d()/d[CH3OCO] */
    dqdci =  - k_r*sc[7];
    J[1125] -= dqdci;             /* dwdot[O]/d[CH3OCO] */
    J[1127] += dqdci;             /* dwdot[OH]/d[CH3OCO] */
    J[1148] += dqdci;             /* dwdot[CH3OCO]/d[CH3OCO] */
    J[1149] -= dqdci;             /* dwdot[CH3OCHO]/d[CH3OCO] */
    /* d()/d[CH3OCHO] */
    dqdci =  + k_f*sc[5];
    J[1165] -= dqdci;             /* dwdot[O]/d[CH3OCHO] */
    J[1167] += dqdci;             /* dwdot[OH]/d[CH3OCHO] */
    J[1188] += dqdci;             /* dwdot[CH3OCO]/d[CH3OCHO] */
    J[1189] -= dqdci;             /* dwdot[CH3OCHO]/d[CH3OCHO] */
    /* d()/dT */
    J[1565] -= dqdT;              /* dwdot[O]/dT */
    J[1567] += dqdT;              /* dwdot[OH]/dT */
    J[1588] += dqdT;              /* dwdot[CH3OCO]/dT */
    J[1589] -= dqdT;              /* dwdot[CH3OCHO]/dT */

    /*reaction 149: CH3OCHO + H <=> CH3OCO + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[29];
    k_f = prefactor_units[148] * fwd_A[148]
                * exp(fwd_beta[148] * tc[0] - activation_units[148] * fwd_Ea[148] * invT);
    dlnkfdT = fwd_beta[148] * invT + activation_units[148] * fwd_Ea[148] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[28];
    Kc = exp(g_RT[0] - g_RT[1] - g_RT[28] + g_RT[29]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[29]) + (h_RT[1] + h_RT[28]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[1] += q; /* H2 */
    wdot[28] += q; /* CH3OCO */
    wdot[29] -= q; /* CH3OCHO */
    /* d()/d[H] */
    dqdci =  + k_f*sc[29];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[1] += dqdci;                /* dwdot[H2]/d[H] */
    J[28] += dqdci;               /* dwdot[CH3OCO]/d[H] */
    J[29] -= dqdci;               /* dwdot[CH3OCHO]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[28];
    J[40] -= dqdci;               /* dwdot[H]/d[H2] */
    J[41] += dqdci;               /* dwdot[H2]/d[H2] */
    J[68] += dqdci;               /* dwdot[CH3OCO]/d[H2] */
    J[69] -= dqdci;               /* dwdot[CH3OCHO]/d[H2] */
    /* d()/d[CH3OCO] */
    dqdci =  - k_r*sc[1];
    J[1120] -= dqdci;             /* dwdot[H]/d[CH3OCO] */
    J[1121] += dqdci;             /* dwdot[H2]/d[CH3OCO] */
    J[1148] += dqdci;             /* dwdot[CH3OCO]/d[CH3OCO] */
    J[1149] -= dqdci;             /* dwdot[CH3OCHO]/d[CH3OCO] */
    /* d()/d[CH3OCHO] */
    dqdci =  + k_f*sc[0];
    J[1160] -= dqdci;             /* dwdot[H]/d[CH3OCHO] */
    J[1161] += dqdci;             /* dwdot[H2]/d[CH3OCHO] */
    J[1188] += dqdci;             /* dwdot[CH3OCO]/d[CH3OCHO] */
    J[1189] -= dqdci;             /* dwdot[CH3OCHO]/d[CH3OCHO] */
    /* d()/dT */
    J[1560] -= dqdT;              /* dwdot[H]/dT */
    J[1561] += dqdT;              /* dwdot[H2]/dT */
    J[1588] += dqdT;              /* dwdot[CH3OCO]/dT */
    J[1589] -= dqdT;              /* dwdot[CH3OCHO]/dT */

    /*reaction 150: CH3OCHO + CH3 <=> CH3OCO + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[29];
    k_f = prefactor_units[149] * fwd_A[149]
                * exp(fwd_beta[149] * tc[0] - activation_units[149] * fwd_Ea[149] * invT);
    dlnkfdT = fwd_beta[149] * invT + activation_units[149] * fwd_Ea[149] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[28];
    Kc = exp(g_RT[4] - g_RT[6] - g_RT[28] + g_RT[29]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[29]) + (h_RT[6] + h_RT[28]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* CH3 */
    wdot[6] += q; /* CH4 */
    wdot[28] += q; /* CH3OCO */
    wdot[29] -= q; /* CH3OCHO */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[29];
    J[164] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[166] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[188] += dqdci;              /* dwdot[CH3OCO]/d[CH3] */
    J[189] -= dqdci;              /* dwdot[CH3OCHO]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[28];
    J[244] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[246] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[268] += dqdci;              /* dwdot[CH3OCO]/d[CH4] */
    J[269] -= dqdci;              /* dwdot[CH3OCHO]/d[CH4] */
    /* d()/d[CH3OCO] */
    dqdci =  - k_r*sc[6];
    J[1124] -= dqdci;             /* dwdot[CH3]/d[CH3OCO] */
    J[1126] += dqdci;             /* dwdot[CH4]/d[CH3OCO] */
    J[1148] += dqdci;             /* dwdot[CH3OCO]/d[CH3OCO] */
    J[1149] -= dqdci;             /* dwdot[CH3OCHO]/d[CH3OCO] */
    /* d()/d[CH3OCHO] */
    dqdci =  + k_f*sc[4];
    J[1164] -= dqdci;             /* dwdot[CH3]/d[CH3OCHO] */
    J[1166] += dqdci;             /* dwdot[CH4]/d[CH3OCHO] */
    J[1188] += dqdci;             /* dwdot[CH3OCO]/d[CH3OCHO] */
    J[1189] -= dqdci;             /* dwdot[CH3OCHO]/d[CH3OCHO] */
    /* d()/dT */
    J[1564] -= dqdT;              /* dwdot[CH3]/dT */
    J[1566] += dqdT;              /* dwdot[CH4]/dT */
    J[1588] += dqdT;              /* dwdot[CH3OCO]/dT */
    J[1589] -= dqdT;              /* dwdot[CH3OCHO]/dT */

    /*reaction 151: CH3OCO <=> CH3O + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[28];
    k_f = prefactor_units[150] * fwd_A[150]
                * exp(fwd_beta[150] * tc[0] - activation_units[150] * fwd_Ea[150] * invT);
    dlnkfdT = fwd_beta[150] * invT + activation_units[150] * fwd_Ea[150] * invT2;
    /* reverse */
    phi_r = sc[11]*sc[18];
    Kc = refC * exp(-g_RT[11] - g_RT[18] + g_RT[28]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[28]) + (h_RT[11] + h_RT[18]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[11] += q; /* CO */
    wdot[18] += q; /* CH3O */
    wdot[28] -= q; /* CH3OCO */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[18];
    J[451] += dqdci;              /* dwdot[CO]/d[CO] */
    J[458] += dqdci;              /* dwdot[CH3O]/d[CO] */
    J[468] -= dqdci;              /* dwdot[CH3OCO]/d[CO] */
    /* d()/d[CH3O] */
    dqdci =  - k_r*sc[11];
    J[731] += dqdci;              /* dwdot[CO]/d[CH3O] */
    J[738] += dqdci;              /* dwdot[CH3O]/d[CH3O] */
    J[748] -= dqdci;              /* dwdot[CH3OCO]/d[CH3O] */
    /* d()/d[CH3OCO] */
    dqdci =  + k_f;
    J[1131] += dqdci;             /* dwdot[CO]/d[CH3OCO] */
    J[1138] += dqdci;             /* dwdot[CH3O]/d[CH3OCO] */
    J[1148] -= dqdci;             /* dwdot[CH3OCO]/d[CH3OCO] */
    /* d()/dT */
    J[1571] += dqdT;              /* dwdot[CO]/dT */
    J[1578] += dqdT;              /* dwdot[CH3O]/dT */
    J[1588] -= dqdT;              /* dwdot[CH3OCO]/dT */

    /*reaction 152: CH3OCO <=> CH3 + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[28];
    k_f = prefactor_units[151] * fwd_A[151]
                * exp(fwd_beta[151] * tc[0] - activation_units[151] * fwd_Ea[151] * invT);
    dlnkfdT = fwd_beta[151] * invT + activation_units[151] * fwd_Ea[151] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[22];
    Kc = refC * exp(-g_RT[4] - g_RT[22] + g_RT[28]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[28]) + (h_RT[4] + h_RT[22]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] += q; /* CH3 */
    wdot[22] += q; /* CO2 */
    wdot[28] -= q; /* CH3OCO */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[22];
    J[164] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[182] += dqdci;              /* dwdot[CO2]/d[CH3] */
    J[188] -= dqdci;              /* dwdot[CH3OCO]/d[CH3] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[4];
    J[884] += dqdci;              /* dwdot[CH3]/d[CO2] */
    J[902] += dqdci;              /* dwdot[CO2]/d[CO2] */
    J[908] -= dqdci;              /* dwdot[CH3OCO]/d[CO2] */
    /* d()/d[CH3OCO] */
    dqdci =  + k_f;
    J[1124] += dqdci;             /* dwdot[CH3]/d[CH3OCO] */
    J[1142] += dqdci;             /* dwdot[CO2]/d[CH3OCO] */
    J[1148] -= dqdci;             /* dwdot[CH3OCO]/d[CH3OCO] */
    /* d()/dT */
    J[1564] += dqdT;              /* dwdot[CH3]/dT */
    J[1582] += dqdT;              /* dwdot[CO2]/dT */
    J[1588] -= dqdT;              /* dwdot[CH3OCO]/dT */

    /*reaction 153: CH3OCH2 + O2 <=> CH3OCH2O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[19]*sc[24];
    k_f = prefactor_units[152] * fwd_A[152]
                * exp(fwd_beta[152] * tc[0] - activation_units[152] * fwd_Ea[152] * invT);
    dlnkfdT = fwd_beta[152] * invT + activation_units[152] * fwd_Ea[152] * invT2;
    /* reverse */
    phi_r = sc[34];
    Kc = refCinv * exp(g_RT[19] + g_RT[24] - g_RT[34]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[19] + h_RT[24]) + (h_RT[34]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[19] -= q; /* O2 */
    wdot[24] -= q; /* CH3OCH2 */
    wdot[34] += q; /* CH3OCH2O2 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[24];
    J[779] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[784] -= dqdci;              /* dwdot[CH3OCH2]/d[O2] */
    J[794] += dqdci;              /* dwdot[CH3OCH2O2]/d[O2] */
    /* d()/d[CH3OCH2] */
    dqdci =  + k_f*sc[19];
    J[979] -= dqdci;              /* dwdot[O2]/d[CH3OCH2] */
    J[984] -= dqdci;              /* dwdot[CH3OCH2]/d[CH3OCH2] */
    J[994] += dqdci;              /* dwdot[CH3OCH2O2]/d[CH3OCH2] */
    /* d()/d[CH3OCH2O2] */
    dqdci =  - k_r;
    J[1379] -= dqdci;             /* dwdot[O2]/d[CH3OCH2O2] */
    J[1384] -= dqdci;             /* dwdot[CH3OCH2]/d[CH3OCH2O2] */
    J[1394] += dqdci;             /* dwdot[CH3OCH2O2]/d[CH3OCH2O2] */
    /* d()/dT */
    J[1579] -= dqdT;              /* dwdot[O2]/dT */
    J[1584] -= dqdT;              /* dwdot[CH3OCH2]/dT */
    J[1594] += dqdT;              /* dwdot[CH3OCH2O2]/dT */

    /*reaction 154: CH3OCH2O2 + CH3OCH2O2 <=> O2 + CH3OCH2O + CH3OCH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[34], 2.000000);
    k_f = prefactor_units[153] * fwd_A[153]
                * exp(fwd_beta[153] * tc[0] - activation_units[153] * fwd_Ea[153] * invT);
    dlnkfdT = fwd_beta[153] * invT + activation_units[153] * fwd_Ea[153] * invT2;
    /* reverse */
    phi_r = sc[19]*pow(sc[30], 2.000000);
    Kc = refC * exp(-g_RT[19] - g_RT[30] - g_RT[30] + g_RT[34] + g_RT[34]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[34]) + (h_RT[19] + 2.000000*h_RT[30]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[19] += q; /* O2 */
    wdot[30] += 2 * q; /* CH3OCH2O */
    wdot[34] -= 2 * q; /* CH3OCH2O2 */
    /* d()/d[O2] */
    dqdci =  - k_r*pow(sc[30], 2.000000);
    J[779] += dqdci;              /* dwdot[O2]/d[O2] */
    J[790] += 2 * dqdci;          /* dwdot[CH3OCH2O]/d[O2] */
    J[794] += -2 * dqdci;         /* dwdot[CH3OCH2O2]/d[O2] */
    /* d()/d[CH3OCH2O] */
    dqdci =  - k_r*sc[19]*2.000000*sc[30];
    J[1219] += dqdci;             /* dwdot[O2]/d[CH3OCH2O] */
    J[1230] += 2 * dqdci;         /* dwdot[CH3OCH2O]/d[CH3OCH2O] */
    J[1234] += -2 * dqdci;        /* dwdot[CH3OCH2O2]/d[CH3OCH2O] */
    /* d()/d[CH3OCH2O2] */
    dqdci =  + k_f*2.000000*sc[34];
    J[1379] += dqdci;             /* dwdot[O2]/d[CH3OCH2O2] */
    J[1390] += 2 * dqdci;         /* dwdot[CH3OCH2O]/d[CH3OCH2O2] */
    J[1394] += -2 * dqdci;        /* dwdot[CH3OCH2O2]/d[CH3OCH2O2] */
    /* d()/dT */
    J[1579] += dqdT;              /* dwdot[O2]/dT */
    J[1590] += 2 * dqdT;          /* dwdot[CH3OCH2O]/dT */
    J[1594] += -2 * dqdT;         /* dwdot[CH3OCH2O2]/dT */

    /*reaction 155: CH3OCH2O2 + CH3OCH2O2 <=> O2 + CH3OCHO + CH3OCH2OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[34], 2.000000);
    k_f = prefactor_units[154] * fwd_A[154]
                * exp(fwd_beta[154] * tc[0] - activation_units[154] * fwd_Ea[154] * invT);
    dlnkfdT = fwd_beta[154] * invT + activation_units[154] * fwd_Ea[154] * invT2;
    /* reverse */
    phi_r = sc[19]*sc[29]*sc[31];
    Kc = refC * exp(-g_RT[19] - g_RT[29] - g_RT[31] + g_RT[34] + g_RT[34]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[34]) + (h_RT[19] + h_RT[29] + h_RT[31]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[19] += q; /* O2 */
    wdot[29] += q; /* CH3OCHO */
    wdot[31] += q; /* CH3OCH2OH */
    wdot[34] -= 2 * q; /* CH3OCH2O2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[29]*sc[31];
    J[779] += dqdci;              /* dwdot[O2]/d[O2] */
    J[789] += dqdci;              /* dwdot[CH3OCHO]/d[O2] */
    J[791] += dqdci;              /* dwdot[CH3OCH2OH]/d[O2] */
    J[794] += -2 * dqdci;         /* dwdot[CH3OCH2O2]/d[O2] */
    /* d()/d[CH3OCHO] */
    dqdci =  - k_r*sc[19]*sc[31];
    J[1179] += dqdci;             /* dwdot[O2]/d[CH3OCHO] */
    J[1189] += dqdci;             /* dwdot[CH3OCHO]/d[CH3OCHO] */
    J[1191] += dqdci;             /* dwdot[CH3OCH2OH]/d[CH3OCHO] */
    J[1194] += -2 * dqdci;        /* dwdot[CH3OCH2O2]/d[CH3OCHO] */
    /* d()/d[CH3OCH2OH] */
    dqdci =  - k_r*sc[19]*sc[29];
    J[1259] += dqdci;             /* dwdot[O2]/d[CH3OCH2OH] */
    J[1269] += dqdci;             /* dwdot[CH3OCHO]/d[CH3OCH2OH] */
    J[1271] += dqdci;             /* dwdot[CH3OCH2OH]/d[CH3OCH2OH] */
    J[1274] += -2 * dqdci;        /* dwdot[CH3OCH2O2]/d[CH3OCH2OH] */
    /* d()/d[CH3OCH2O2] */
    dqdci =  + k_f*2.000000*sc[34];
    J[1379] += dqdci;             /* dwdot[O2]/d[CH3OCH2O2] */
    J[1389] += dqdci;             /* dwdot[CH3OCHO]/d[CH3OCH2O2] */
    J[1391] += dqdci;             /* dwdot[CH3OCH2OH]/d[CH3OCH2O2] */
    J[1394] += -2 * dqdci;        /* dwdot[CH3OCH2O2]/d[CH3OCH2O2] */
    /* d()/dT */
    J[1579] += dqdT;              /* dwdot[O2]/dT */
    J[1589] += dqdT;              /* dwdot[CH3OCHO]/dT */
    J[1591] += dqdT;              /* dwdot[CH3OCH2OH]/dT */
    J[1594] += -2 * dqdT;         /* dwdot[CH3OCH2O2]/dT */

    /*reaction 156: CH3OCH2O <=> CH3O + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[30];
    k_f = prefactor_units[155] * fwd_A[155]
                * exp(fwd_beta[155] * tc[0] - activation_units[155] * fwd_Ea[155] * invT);
    dlnkfdT = fwd_beta[155] * invT + activation_units[155] * fwd_Ea[155] * invT2;
    /* reverse */
    phi_r = sc[15]*sc[18];
    Kc = refC * exp(-g_RT[15] - g_RT[18] + g_RT[30]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[30]) + (h_RT[15] + h_RT[18]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[15] += q; /* CH2O */
    wdot[18] += q; /* CH3O */
    wdot[30] -= q; /* CH3OCH2O */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[18];
    J[615] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[618] += dqdci;              /* dwdot[CH3O]/d[CH2O] */
    J[630] -= dqdci;              /* dwdot[CH3OCH2O]/d[CH2O] */
    /* d()/d[CH3O] */
    dqdci =  - k_r*sc[15];
    J[735] += dqdci;              /* dwdot[CH2O]/d[CH3O] */
    J[738] += dqdci;              /* dwdot[CH3O]/d[CH3O] */
    J[750] -= dqdci;              /* dwdot[CH3OCH2O]/d[CH3O] */
    /* d()/d[CH3OCH2O] */
    dqdci =  + k_f;
    J[1215] += dqdci;             /* dwdot[CH2O]/d[CH3OCH2O] */
    J[1218] += dqdci;             /* dwdot[CH3O]/d[CH3OCH2O] */
    J[1230] -= dqdci;             /* dwdot[CH3OCH2O]/d[CH3OCH2O] */
    /* d()/dT */
    J[1575] += dqdT;              /* dwdot[CH2O]/dT */
    J[1578] += dqdT;              /* dwdot[CH3O]/dT */
    J[1590] -= dqdT;              /* dwdot[CH3OCH2O]/dT */

    /*reaction 157: CH3OCH2O + O2 <=> CH3OCHO + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[19]*sc[30];
    k_f = prefactor_units[156] * fwd_A[156]
                * exp(fwd_beta[156] * tc[0] - activation_units[156] * fwd_Ea[156] * invT);
    dlnkfdT = fwd_beta[156] * invT + activation_units[156] * fwd_Ea[156] * invT2;
    /* reverse */
    phi_r = sc[20]*sc[29];
    Kc = exp(g_RT[19] - g_RT[20] - g_RT[29] + g_RT[30]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[19] + h_RT[30]) + (h_RT[20] + h_RT[29]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[19] -= q; /* O2 */
    wdot[20] += q; /* HO2 */
    wdot[29] += q; /* CH3OCHO */
    wdot[30] -= q; /* CH3OCH2O */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[30];
    J[779] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[780] += dqdci;              /* dwdot[HO2]/d[O2] */
    J[789] += dqdci;              /* dwdot[CH3OCHO]/d[O2] */
    J[790] -= dqdci;              /* dwdot[CH3OCH2O]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[29];
    J[819] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[820] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[829] += dqdci;              /* dwdot[CH3OCHO]/d[HO2] */
    J[830] -= dqdci;              /* dwdot[CH3OCH2O]/d[HO2] */
    /* d()/d[CH3OCHO] */
    dqdci =  - k_r*sc[20];
    J[1179] -= dqdci;             /* dwdot[O2]/d[CH3OCHO] */
    J[1180] += dqdci;             /* dwdot[HO2]/d[CH3OCHO] */
    J[1189] += dqdci;             /* dwdot[CH3OCHO]/d[CH3OCHO] */
    J[1190] -= dqdci;             /* dwdot[CH3OCH2O]/d[CH3OCHO] */
    /* d()/d[CH3OCH2O] */
    dqdci =  + k_f*sc[19];
    J[1219] -= dqdci;             /* dwdot[O2]/d[CH3OCH2O] */
    J[1220] += dqdci;             /* dwdot[HO2]/d[CH3OCH2O] */
    J[1229] += dqdci;             /* dwdot[CH3OCHO]/d[CH3OCH2O] */
    J[1230] -= dqdci;             /* dwdot[CH3OCH2O]/d[CH3OCH2O] */
    /* d()/dT */
    J[1579] -= dqdT;              /* dwdot[O2]/dT */
    J[1580] += dqdT;              /* dwdot[HO2]/dT */
    J[1589] += dqdT;              /* dwdot[CH3OCHO]/dT */
    J[1590] -= dqdT;              /* dwdot[CH3OCH2O]/dT */

    /*reaction 158: CH3OCH2O2 <=> CH2OCH2O2H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[34];
    k_f = prefactor_units[157] * fwd_A[157]
                * exp(fwd_beta[157] * tc[0] - activation_units[157] * fwd_Ea[157] * invT);
    dlnkfdT = fwd_beta[157] * invT + activation_units[157] * fwd_Ea[157] * invT2;
    /* reverse */
    phi_r = sc[35];
    Kc = exp(g_RT[34] - g_RT[35]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[34]) + (h_RT[35]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[34] -= q; /* CH3OCH2O2 */
    wdot[35] += q; /* CH2OCH2O2H */
    /* d()/d[CH3OCH2O2] */
    dqdci =  + k_f;
    J[1394] -= dqdci;             /* dwdot[CH3OCH2O2]/d[CH3OCH2O2] */
    J[1395] += dqdci;             /* dwdot[CH2OCH2O2H]/d[CH3OCH2O2] */
    /* d()/d[CH2OCH2O2H] */
    dqdci =  - k_r;
    J[1434] -= dqdci;             /* dwdot[CH3OCH2O2]/d[CH2OCH2O2H] */
    J[1435] += dqdci;             /* dwdot[CH2OCH2O2H]/d[CH2OCH2O2H] */
    /* d()/dT */
    J[1594] -= dqdT;              /* dwdot[CH3OCH2O2]/dT */
    J[1595] += dqdT;              /* dwdot[CH2OCH2O2H]/dT */

    /*reaction 159: CH2OCH2O2H <=> OH + CH2O + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[35];
    k_f = prefactor_units[158] * fwd_A[158]
                * exp(fwd_beta[158] * tc[0] - activation_units[158] * fwd_Ea[158] * invT);
    dlnkfdT = fwd_beta[158] * invT + activation_units[158] * fwd_Ea[158] * invT2;
    /* reverse */
    phi_r = sc[7]*pow(sc[15], 2.000000);
    Kc = pow(refC,2.000000) * exp(-g_RT[7] - g_RT[15] - g_RT[15] + g_RT[35]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[35]) + (h_RT[7] + 2.000000*h_RT[15]) - 2.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] += q; /* OH */
    wdot[15] += 2 * q; /* CH2O */
    wdot[35] -= q; /* CH2OCH2O2H */
    /* d()/d[OH] */
    dqdci =  - k_r*pow(sc[15], 2.000000);
    J[287] += dqdci;              /* dwdot[OH]/d[OH] */
    J[295] += 2 * dqdci;          /* dwdot[CH2O]/d[OH] */
    J[315] -= dqdci;              /* dwdot[CH2OCH2O2H]/d[OH] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[7]*2.000000*sc[15];
    J[607] += dqdci;              /* dwdot[OH]/d[CH2O] */
    J[615] += 2 * dqdci;          /* dwdot[CH2O]/d[CH2O] */
    J[635] -= dqdci;              /* dwdot[CH2OCH2O2H]/d[CH2O] */
    /* d()/d[CH2OCH2O2H] */
    dqdci =  + k_f;
    J[1407] += dqdci;             /* dwdot[OH]/d[CH2OCH2O2H] */
    J[1415] += 2 * dqdci;         /* dwdot[CH2O]/d[CH2OCH2O2H] */
    J[1435] -= dqdci;             /* dwdot[CH2OCH2O2H]/d[CH2OCH2O2H] */
    /* d()/dT */
    J[1567] += dqdT;              /* dwdot[OH]/dT */
    J[1575] += 2 * dqdT;          /* dwdot[CH2O]/dT */
    J[1595] -= dqdT;              /* dwdot[CH2OCH2O2H]/dT */

    /*reaction 160: CH2OCH2O2H + O2 <=> O2CH2OCH2O2H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[19]*sc[35];
    k_f = prefactor_units[159] * fwd_A[159]
                * exp(fwd_beta[159] * tc[0] - activation_units[159] * fwd_Ea[159] * invT);
    dlnkfdT = fwd_beta[159] * invT + activation_units[159] * fwd_Ea[159] * invT2;
    /* reverse */
    phi_r = sc[37];
    Kc = refCinv * exp(g_RT[19] + g_RT[35] - g_RT[37]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[19] + h_RT[35]) + (h_RT[37]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[19] -= q; /* O2 */
    wdot[35] -= q; /* CH2OCH2O2H */
    wdot[37] += q; /* O2CH2OCH2O2H */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[35];
    J[779] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[795] -= dqdci;              /* dwdot[CH2OCH2O2H]/d[O2] */
    J[797] += dqdci;              /* dwdot[O2CH2OCH2O2H]/d[O2] */
    /* d()/d[CH2OCH2O2H] */
    dqdci =  + k_f*sc[19];
    J[1419] -= dqdci;             /* dwdot[O2]/d[CH2OCH2O2H] */
    J[1435] -= dqdci;             /* dwdot[CH2OCH2O2H]/d[CH2OCH2O2H] */
    J[1437] += dqdci;             /* dwdot[O2CH2OCH2O2H]/d[CH2OCH2O2H] */
    /* d()/d[O2CH2OCH2O2H] */
    dqdci =  - k_r;
    J[1499] -= dqdci;             /* dwdot[O2]/d[O2CH2OCH2O2H] */
    J[1515] -= dqdci;             /* dwdot[CH2OCH2O2H]/d[O2CH2OCH2O2H] */
    J[1517] += dqdci;             /* dwdot[O2CH2OCH2O2H]/d[O2CH2OCH2O2H] */
    /* d()/dT */
    J[1579] -= dqdT;              /* dwdot[O2]/dT */
    J[1595] -= dqdT;              /* dwdot[CH2OCH2O2H]/dT */
    J[1597] += dqdT;              /* dwdot[O2CH2OCH2O2H]/dT */

    /*reaction 161: O2CH2OCH2O2H <=> HO2CH2OCHO + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[37];
    k_f = prefactor_units[160] * fwd_A[160]
                * exp(fwd_beta[160] * tc[0] - activation_units[160] * fwd_Ea[160] * invT);
    dlnkfdT = fwd_beta[160] * invT + activation_units[160] * fwd_Ea[160] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[36];
    Kc = refC * exp(-g_RT[7] - g_RT[36] + g_RT[37]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[37]) + (h_RT[7] + h_RT[36]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] += q; /* OH */
    wdot[36] += q; /* HO2CH2OCHO */
    wdot[37] -= q; /* O2CH2OCH2O2H */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[36];
    J[287] += dqdci;              /* dwdot[OH]/d[OH] */
    J[316] += dqdci;              /* dwdot[HO2CH2OCHO]/d[OH] */
    J[317] -= dqdci;              /* dwdot[O2CH2OCH2O2H]/d[OH] */
    /* d()/d[HO2CH2OCHO] */
    dqdci =  - k_r*sc[7];
    J[1447] += dqdci;             /* dwdot[OH]/d[HO2CH2OCHO] */
    J[1476] += dqdci;             /* dwdot[HO2CH2OCHO]/d[HO2CH2OCHO] */
    J[1477] -= dqdci;             /* dwdot[O2CH2OCH2O2H]/d[HO2CH2OCHO] */
    /* d()/d[O2CH2OCH2O2H] */
    dqdci =  + k_f;
    J[1487] += dqdci;             /* dwdot[OH]/d[O2CH2OCH2O2H] */
    J[1516] += dqdci;             /* dwdot[HO2CH2OCHO]/d[O2CH2OCH2O2H] */
    J[1517] -= dqdci;             /* dwdot[O2CH2OCH2O2H]/d[O2CH2OCH2O2H] */
    /* d()/dT */
    J[1567] += dqdT;              /* dwdot[OH]/dT */
    J[1596] += dqdT;              /* dwdot[HO2CH2OCHO]/dT */
    J[1597] -= dqdT;              /* dwdot[O2CH2OCH2O2H]/dT */

    /*reaction 162: HO2CH2OCHO <=> OCH2OCHO + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[36];
    k_f = prefactor_units[161] * fwd_A[161]
                * exp(fwd_beta[161] * tc[0] - activation_units[161] * fwd_Ea[161] * invT);
    dlnkfdT = fwd_beta[161] * invT + activation_units[161] * fwd_Ea[161] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[32];
    Kc = refC * exp(-g_RT[7] - g_RT[32] + g_RT[36]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[36]) + (h_RT[7] + h_RT[32]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] += q; /* OH */
    wdot[32] += q; /* OCH2OCHO */
    wdot[36] -= q; /* HO2CH2OCHO */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[32];
    J[287] += dqdci;              /* dwdot[OH]/d[OH] */
    J[312] += dqdci;              /* dwdot[OCH2OCHO]/d[OH] */
    J[316] -= dqdci;              /* dwdot[HO2CH2OCHO]/d[OH] */
    /* d()/d[OCH2OCHO] */
    dqdci =  - k_r*sc[7];
    J[1287] += dqdci;             /* dwdot[OH]/d[OCH2OCHO] */
    J[1312] += dqdci;             /* dwdot[OCH2OCHO]/d[OCH2OCHO] */
    J[1316] -= dqdci;             /* dwdot[HO2CH2OCHO]/d[OCH2OCHO] */
    /* d()/d[HO2CH2OCHO] */
    dqdci =  + k_f;
    J[1447] += dqdci;             /* dwdot[OH]/d[HO2CH2OCHO] */
    J[1472] += dqdci;             /* dwdot[OCH2OCHO]/d[HO2CH2OCHO] */
    J[1476] -= dqdci;             /* dwdot[HO2CH2OCHO]/d[HO2CH2OCHO] */
    /* d()/dT */
    J[1567] += dqdT;              /* dwdot[OH]/dT */
    J[1592] += dqdT;              /* dwdot[OCH2OCHO]/dT */
    J[1596] -= dqdT;              /* dwdot[HO2CH2OCHO]/dT */

    /*reaction 163: OCH2OCHO <=> HOCH2OCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[32];
    k_f = prefactor_units[162] * fwd_A[162]
                * exp(fwd_beta[162] * tc[0] - activation_units[162] * fwd_Ea[162] * invT);
    dlnkfdT = fwd_beta[162] * invT + activation_units[162] * fwd_Ea[162] * invT2;
    /* reverse */
    phi_r = sc[33];
    Kc = exp(g_RT[32] - g_RT[33]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[32]) + (h_RT[33]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[32] -= q; /* OCH2OCHO */
    wdot[33] += q; /* HOCH2OCO */
    /* d()/d[OCH2OCHO] */
    dqdci =  + k_f;
    J[1312] -= dqdci;             /* dwdot[OCH2OCHO]/d[OCH2OCHO] */
    J[1313] += dqdci;             /* dwdot[HOCH2OCO]/d[OCH2OCHO] */
    /* d()/d[HOCH2OCO] */
    dqdci =  - k_r;
    J[1352] -= dqdci;             /* dwdot[OCH2OCHO]/d[HOCH2OCO] */
    J[1353] += dqdci;             /* dwdot[HOCH2OCO]/d[HOCH2OCO] */
    /* d()/dT */
    J[1592] -= dqdT;              /* dwdot[OCH2OCHO]/dT */
    J[1593] += dqdT;              /* dwdot[HOCH2OCO]/dT */

    /*reaction 164: HOCH2OCO <=> HOCH2O + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[33];
    k_f = prefactor_units[163] * fwd_A[163]
                * exp(fwd_beta[163] * tc[0] - activation_units[163] * fwd_Ea[163] * invT);
    dlnkfdT = fwd_beta[163] * invT + activation_units[163] * fwd_Ea[163] * invT2;
    /* reverse */
    phi_r = sc[11]*sc[27];
    Kc = refC * exp(-g_RT[11] - g_RT[27] + g_RT[33]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[33]) + (h_RT[11] + h_RT[27]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[11] += q; /* CO */
    wdot[27] += q; /* HOCH2O */
    wdot[33] -= q; /* HOCH2OCO */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[27];
    J[451] += dqdci;              /* dwdot[CO]/d[CO] */
    J[467] += dqdci;              /* dwdot[HOCH2O]/d[CO] */
    J[473] -= dqdci;              /* dwdot[HOCH2OCO]/d[CO] */
    /* d()/d[HOCH2O] */
    dqdci =  - k_r*sc[11];
    J[1091] += dqdci;             /* dwdot[CO]/d[HOCH2O] */
    J[1107] += dqdci;             /* dwdot[HOCH2O]/d[HOCH2O] */
    J[1113] -= dqdci;             /* dwdot[HOCH2OCO]/d[HOCH2O] */
    /* d()/d[HOCH2OCO] */
    dqdci =  + k_f;
    J[1331] += dqdci;             /* dwdot[CO]/d[HOCH2OCO] */
    J[1347] += dqdci;             /* dwdot[HOCH2O]/d[HOCH2OCO] */
    J[1353] -= dqdci;             /* dwdot[HOCH2OCO]/d[HOCH2OCO] */
    /* d()/dT */
    J[1571] += dqdT;              /* dwdot[CO]/dT */
    J[1587] += dqdT;              /* dwdot[HOCH2O]/dT */
    J[1593] -= dqdT;              /* dwdot[HOCH2OCO]/dT */

    /*reaction 165: HOCH2OCO <=> CH2OH + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[33];
    k_f = prefactor_units[164] * fwd_A[164]
                * exp(fwd_beta[164] * tc[0] - activation_units[164] * fwd_Ea[164] * invT);
    dlnkfdT = fwd_beta[164] * invT + activation_units[164] * fwd_Ea[164] * invT2;
    /* reverse */
    phi_r = sc[17]*sc[22];
    Kc = refC * exp(-g_RT[17] - g_RT[22] + g_RT[33]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[33]) + (h_RT[17] + h_RT[22]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[17] += q; /* CH2OH */
    wdot[22] += q; /* CO2 */
    wdot[33] -= q; /* HOCH2OCO */
    /* d()/d[CH2OH] */
    dqdci =  - k_r*sc[22];
    J[697] += dqdci;              /* dwdot[CH2OH]/d[CH2OH] */
    J[702] += dqdci;              /* dwdot[CO2]/d[CH2OH] */
    J[713] -= dqdci;              /* dwdot[HOCH2OCO]/d[CH2OH] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[17];
    J[897] += dqdci;              /* dwdot[CH2OH]/d[CO2] */
    J[902] += dqdci;              /* dwdot[CO2]/d[CO2] */
    J[913] -= dqdci;              /* dwdot[HOCH2OCO]/d[CO2] */
    /* d()/d[HOCH2OCO] */
    dqdci =  + k_f;
    J[1337] += dqdci;             /* dwdot[CH2OH]/d[HOCH2OCO] */
    J[1342] += dqdci;             /* dwdot[CO2]/d[HOCH2OCO] */
    J[1353] -= dqdci;             /* dwdot[HOCH2OCO]/d[HOCH2OCO] */
    /* d()/dT */
    J[1577] += dqdT;              /* dwdot[CH2OH]/dT */
    J[1582] += dqdT;              /* dwdot[CO2]/dT */
    J[1593] -= dqdT;              /* dwdot[HOCH2OCO]/dT */

    /*reaction 166: HOCH2O <=> HCOOH + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[27];
    k_f = prefactor_units[165] * fwd_A[165]
                * exp(fwd_beta[165] * tc[0] - activation_units[165] * fwd_Ea[165] * invT);
    dlnkfdT = fwd_beta[165] * invT + activation_units[165] * fwd_Ea[165] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[25];
    Kc = refC * exp(-g_RT[0] - g_RT[25] + g_RT[27]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[27]) + (h_RT[0] + h_RT[25]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[25] += q; /* HCOOH */
    wdot[27] -= q; /* HOCH2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[25];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[25] += dqdci;               /* dwdot[HCOOH]/d[H] */
    J[27] -= dqdci;               /* dwdot[HOCH2O]/d[H] */
    /* d()/d[HCOOH] */
    dqdci =  - k_r*sc[0];
    J[1000] += dqdci;             /* dwdot[H]/d[HCOOH] */
    J[1025] += dqdci;             /* dwdot[HCOOH]/d[HCOOH] */
    J[1027] -= dqdci;             /* dwdot[HOCH2O]/d[HCOOH] */
    /* d()/d[HOCH2O] */
    dqdci =  + k_f;
    J[1080] += dqdci;             /* dwdot[H]/d[HOCH2O] */
    J[1105] += dqdci;             /* dwdot[HCOOH]/d[HOCH2O] */
    J[1107] -= dqdci;             /* dwdot[HOCH2O]/d[HOCH2O] */
    /* d()/dT */
    J[1560] += dqdT;              /* dwdot[H]/dT */
    J[1585] += dqdT;              /* dwdot[HCOOH]/dT */
    J[1587] -= dqdT;              /* dwdot[HOCH2O]/dT */

    /*reaction 167: CH2O + OH <=> HOCH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[7]*sc[15];
    k_f = prefactor_units[166] * fwd_A[166]
                * exp(fwd_beta[166] * tc[0] - activation_units[166] * fwd_Ea[166] * invT);
    dlnkfdT = fwd_beta[166] * invT + activation_units[166] * fwd_Ea[166] * invT2;
    /* reverse */
    phi_r = sc[27];
    Kc = refCinv * exp(g_RT[7] + g_RT[15] - g_RT[27]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7] + h_RT[15]) + (h_RT[27]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] -= q; /* OH */
    wdot[15] -= q; /* CH2O */
    wdot[27] += q; /* HOCH2O */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[15];
    J[287] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[295] -= dqdci;              /* dwdot[CH2O]/d[OH] */
    J[307] += dqdci;              /* dwdot[HOCH2O]/d[OH] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[7];
    J[607] -= dqdci;              /* dwdot[OH]/d[CH2O] */
    J[615] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[627] += dqdci;              /* dwdot[HOCH2O]/d[CH2O] */
    /* d()/d[HOCH2O] */
    dqdci =  - k_r;
    J[1087] -= dqdci;             /* dwdot[OH]/d[HOCH2O] */
    J[1095] -= dqdci;             /* dwdot[CH2O]/d[HOCH2O] */
    J[1107] += dqdci;             /* dwdot[HOCH2O]/d[HOCH2O] */
    /* d()/dT */
    J[1567] -= dqdT;              /* dwdot[OH]/dT */
    J[1575] -= dqdT;              /* dwdot[CH2O]/dT */
    J[1587] += dqdT;              /* dwdot[HOCH2O]/dT */

    /*reaction 168: HCOOH <=> HCO + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[25];
    k_f = prefactor_units[167] * fwd_A[167]
                * exp(fwd_beta[167] * tc[0] - activation_units[167] * fwd_Ea[167] * invT);
    dlnkfdT = fwd_beta[167] * invT + activation_units[167] * fwd_Ea[167] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[13];
    Kc = refC * exp(-g_RT[7] - g_RT[13] + g_RT[25]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[25]) + (h_RT[7] + h_RT[13]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] += q; /* OH */
    wdot[13] += q; /* HCO */
    wdot[25] -= q; /* HCOOH */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[13];
    J[287] += dqdci;              /* dwdot[OH]/d[OH] */
    J[293] += dqdci;              /* dwdot[HCO]/d[OH] */
    J[305] -= dqdci;              /* dwdot[HCOOH]/d[OH] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[7];
    J[527] += dqdci;              /* dwdot[OH]/d[HCO] */
    J[533] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[545] -= dqdci;              /* dwdot[HCOOH]/d[HCO] */
    /* d()/d[HCOOH] */
    dqdci =  + k_f;
    J[1007] += dqdci;             /* dwdot[OH]/d[HCOOH] */
    J[1013] += dqdci;             /* dwdot[HCO]/d[HCOOH] */
    J[1025] -= dqdci;             /* dwdot[HCOOH]/d[HCOOH] */
    /* d()/dT */
    J[1567] += dqdT;              /* dwdot[OH]/dT */
    J[1573] += dqdT;              /* dwdot[HCO]/dT */
    J[1585] -= dqdT;              /* dwdot[HCOOH]/dT */

    /*reaction 169: HCOOH + OH <=> H2O + CO2 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[7]*sc[25];
    k_f = prefactor_units[168] * fwd_A[168]
                * exp(fwd_beta[168] * tc[0] - activation_units[168] * fwd_Ea[168] * invT);
    dlnkfdT = fwd_beta[168] * invT + activation_units[168] * fwd_Ea[168] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[8]*sc[22];
    Kc = refC * exp(-g_RT[0] + g_RT[7] - g_RT[8] - g_RT[22] + g_RT[25]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7] + h_RT[25]) + (h_RT[0] + h_RT[8] + h_RT[22]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[7] -= q; /* OH */
    wdot[8] += q; /* H2O */
    wdot[22] += q; /* CO2 */
    wdot[25] -= q; /* HCOOH */
    /* d()/d[H] */
    dqdci =  - k_r*sc[8]*sc[22];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[7] -= dqdci;                /* dwdot[OH]/d[H] */
    J[8] += dqdci;                /* dwdot[H2O]/d[H] */
    J[22] += dqdci;               /* dwdot[CO2]/d[H] */
    J[25] -= dqdci;               /* dwdot[HCOOH]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[25];
    J[280] += dqdci;              /* dwdot[H]/d[OH] */
    J[287] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[288] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[302] += dqdci;              /* dwdot[CO2]/d[OH] */
    J[305] -= dqdci;              /* dwdot[HCOOH]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[0]*sc[22];
    J[320] += dqdci;              /* dwdot[H]/d[H2O] */
    J[327] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[328] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[342] += dqdci;              /* dwdot[CO2]/d[H2O] */
    J[345] -= dqdci;              /* dwdot[HCOOH]/d[H2O] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[0]*sc[8];
    J[880] += dqdci;              /* dwdot[H]/d[CO2] */
    J[887] -= dqdci;              /* dwdot[OH]/d[CO2] */
    J[888] += dqdci;              /* dwdot[H2O]/d[CO2] */
    J[902] += dqdci;              /* dwdot[CO2]/d[CO2] */
    J[905] -= dqdci;              /* dwdot[HCOOH]/d[CO2] */
    /* d()/d[HCOOH] */
    dqdci =  + k_f*sc[7];
    J[1000] += dqdci;             /* dwdot[H]/d[HCOOH] */
    J[1007] -= dqdci;             /* dwdot[OH]/d[HCOOH] */
    J[1008] += dqdci;             /* dwdot[H2O]/d[HCOOH] */
    J[1022] += dqdci;             /* dwdot[CO2]/d[HCOOH] */
    J[1025] -= dqdci;             /* dwdot[HCOOH]/d[HCOOH] */
    /* d()/dT */
    J[1560] += dqdT;              /* dwdot[H]/dT */
    J[1567] -= dqdT;              /* dwdot[OH]/dT */
    J[1568] += dqdT;              /* dwdot[H2O]/dT */
    J[1582] += dqdT;              /* dwdot[CO2]/dT */
    J[1585] -= dqdT;              /* dwdot[HCOOH]/dT */

    /*reaction 170: HCOOH + OH <=> H2O + CO + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[7]*sc[25];
    k_f = prefactor_units[169] * fwd_A[169]
                * exp(fwd_beta[169] * tc[0] - activation_units[169] * fwd_Ea[169] * invT);
    dlnkfdT = fwd_beta[169] * invT + activation_units[169] * fwd_Ea[169] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[8]*sc[11];
    Kc = refC * exp(g_RT[7] - g_RT[7] - g_RT[8] - g_RT[11] + g_RT[25]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7] + h_RT[25]) + (h_RT[7] + h_RT[8] + h_RT[11]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[8] += q; /* H2O */
    wdot[11] += q; /* CO */
    wdot[25] -= q; /* HCOOH */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[25] - k_r*sc[8]*sc[11];
    J[288] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[291] += dqdci;              /* dwdot[CO]/d[OH] */
    J[305] -= dqdci;              /* dwdot[HCOOH]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[7]*sc[11];
    J[328] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[331] += dqdci;              /* dwdot[CO]/d[H2O] */
    J[345] -= dqdci;              /* dwdot[HCOOH]/d[H2O] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[7]*sc[8];
    J[448] += dqdci;              /* dwdot[H2O]/d[CO] */
    J[451] += dqdci;              /* dwdot[CO]/d[CO] */
    J[465] -= dqdci;              /* dwdot[HCOOH]/d[CO] */
    /* d()/d[HCOOH] */
    dqdci =  + k_f*sc[7];
    J[1008] += dqdci;             /* dwdot[H2O]/d[HCOOH] */
    J[1011] += dqdci;             /* dwdot[CO]/d[HCOOH] */
    J[1025] -= dqdci;             /* dwdot[HCOOH]/d[HCOOH] */
    /* d()/dT */
    J[1568] += dqdT;              /* dwdot[H2O]/dT */
    J[1571] += dqdT;              /* dwdot[CO]/dT */
    J[1585] -= dqdT;              /* dwdot[HCOOH]/dT */

    /*reaction 171: HCOOH + H <=> H2 + CO2 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[25];
    k_f = prefactor_units[170] * fwd_A[170]
                * exp(fwd_beta[170] * tc[0] - activation_units[170] * fwd_Ea[170] * invT);
    dlnkfdT = fwd_beta[170] * invT + activation_units[170] * fwd_Ea[170] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[1]*sc[22];
    Kc = refC * exp(g_RT[0] - g_RT[0] - g_RT[1] - g_RT[22] + g_RT[25]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[25]) + (h_RT[0] + h_RT[1] + h_RT[22]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H2 */
    wdot[22] += q; /* CO2 */
    wdot[25] -= q; /* HCOOH */
    /* d()/d[H] */
    dqdci =  + k_f*sc[25] - k_r*sc[1]*sc[22];
    J[1] += dqdci;                /* dwdot[H2]/d[H] */
    J[22] += dqdci;               /* dwdot[CO2]/d[H] */
    J[25] -= dqdci;               /* dwdot[HCOOH]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[0]*sc[22];
    J[41] += dqdci;               /* dwdot[H2]/d[H2] */
    J[62] += dqdci;               /* dwdot[CO2]/d[H2] */
    J[65] -= dqdci;               /* dwdot[HCOOH]/d[H2] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[0]*sc[1];
    J[881] += dqdci;              /* dwdot[H2]/d[CO2] */
    J[902] += dqdci;              /* dwdot[CO2]/d[CO2] */
    J[905] -= dqdci;              /* dwdot[HCOOH]/d[CO2] */
    /* d()/d[HCOOH] */
    dqdci =  + k_f*sc[0];
    J[1001] += dqdci;             /* dwdot[H2]/d[HCOOH] */
    J[1022] += dqdci;             /* dwdot[CO2]/d[HCOOH] */
    J[1025] -= dqdci;             /* dwdot[HCOOH]/d[HCOOH] */
    /* d()/dT */
    J[1561] += dqdT;              /* dwdot[H2]/dT */
    J[1582] += dqdT;              /* dwdot[CO2]/dT */
    J[1585] -= dqdT;              /* dwdot[HCOOH]/dT */

    /*reaction 172: HCOOH + H <=> H2 + CO + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[25];
    k_f = prefactor_units[171] * fwd_A[171]
                * exp(fwd_beta[171] * tc[0] - activation_units[171] * fwd_Ea[171] * invT);
    dlnkfdT = fwd_beta[171] * invT + activation_units[171] * fwd_Ea[171] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[7]*sc[11];
    Kc = refC * exp(g_RT[0] - g_RT[1] - g_RT[7] - g_RT[11] + g_RT[25]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[25]) + (h_RT[1] + h_RT[7] + h_RT[11]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[1] += q; /* H2 */
    wdot[7] += q; /* OH */
    wdot[11] += q; /* CO */
    wdot[25] -= q; /* HCOOH */
    /* d()/d[H] */
    dqdci =  + k_f*sc[25];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[1] += dqdci;                /* dwdot[H2]/d[H] */
    J[7] += dqdci;                /* dwdot[OH]/d[H] */
    J[11] += dqdci;               /* dwdot[CO]/d[H] */
    J[25] -= dqdci;               /* dwdot[HCOOH]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[7]*sc[11];
    J[40] -= dqdci;               /* dwdot[H]/d[H2] */
    J[41] += dqdci;               /* dwdot[H2]/d[H2] */
    J[47] += dqdci;               /* dwdot[OH]/d[H2] */
    J[51] += dqdci;               /* dwdot[CO]/d[H2] */
    J[65] -= dqdci;               /* dwdot[HCOOH]/d[H2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[1]*sc[11];
    J[280] -= dqdci;              /* dwdot[H]/d[OH] */
    J[281] += dqdci;              /* dwdot[H2]/d[OH] */
    J[287] += dqdci;              /* dwdot[OH]/d[OH] */
    J[291] += dqdci;              /* dwdot[CO]/d[OH] */
    J[305] -= dqdci;              /* dwdot[HCOOH]/d[OH] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[1]*sc[7];
    J[440] -= dqdci;              /* dwdot[H]/d[CO] */
    J[441] += dqdci;              /* dwdot[H2]/d[CO] */
    J[447] += dqdci;              /* dwdot[OH]/d[CO] */
    J[451] += dqdci;              /* dwdot[CO]/d[CO] */
    J[465] -= dqdci;              /* dwdot[HCOOH]/d[CO] */
    /* d()/d[HCOOH] */
    dqdci =  + k_f*sc[0];
    J[1000] -= dqdci;             /* dwdot[H]/d[HCOOH] */
    J[1001] += dqdci;             /* dwdot[H2]/d[HCOOH] */
    J[1007] += dqdci;             /* dwdot[OH]/d[HCOOH] */
    J[1011] += dqdci;             /* dwdot[CO]/d[HCOOH] */
    J[1025] -= dqdci;             /* dwdot[HCOOH]/d[HCOOH] */
    /* d()/dT */
    J[1560] -= dqdT;              /* dwdot[H]/dT */
    J[1561] += dqdT;              /* dwdot[H2]/dT */
    J[1567] += dqdT;              /* dwdot[OH]/dT */
    J[1571] += dqdT;              /* dwdot[CO]/dT */
    J[1585] -= dqdT;              /* dwdot[HCOOH]/dT */

    /*reaction 173: HCOOH + CH3 <=> CH4 + CO + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[25];
    k_f = prefactor_units[172] * fwd_A[172]
                * exp(fwd_beta[172] * tc[0] - activation_units[172] * fwd_Ea[172] * invT);
    dlnkfdT = fwd_beta[172] * invT + activation_units[172] * fwd_Ea[172] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[7]*sc[11];
    Kc = refC * exp(g_RT[4] - g_RT[6] - g_RT[7] - g_RT[11] + g_RT[25]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[25]) + (h_RT[6] + h_RT[7] + h_RT[11]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* CH3 */
    wdot[6] += q; /* CH4 */
    wdot[7] += q; /* OH */
    wdot[11] += q; /* CO */
    wdot[25] -= q; /* HCOOH */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[25];
    J[164] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[166] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[167] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[171] += dqdci;              /* dwdot[CO]/d[CH3] */
    J[185] -= dqdci;              /* dwdot[HCOOH]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[7]*sc[11];
    J[244] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[246] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[247] += dqdci;              /* dwdot[OH]/d[CH4] */
    J[251] += dqdci;              /* dwdot[CO]/d[CH4] */
    J[265] -= dqdci;              /* dwdot[HCOOH]/d[CH4] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[6]*sc[11];
    J[284] -= dqdci;              /* dwdot[CH3]/d[OH] */
    J[286] += dqdci;              /* dwdot[CH4]/d[OH] */
    J[287] += dqdci;              /* dwdot[OH]/d[OH] */
    J[291] += dqdci;              /* dwdot[CO]/d[OH] */
    J[305] -= dqdci;              /* dwdot[HCOOH]/d[OH] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[6]*sc[7];
    J[444] -= dqdci;              /* dwdot[CH3]/d[CO] */
    J[446] += dqdci;              /* dwdot[CH4]/d[CO] */
    J[447] += dqdci;              /* dwdot[OH]/d[CO] */
    J[451] += dqdci;              /* dwdot[CO]/d[CO] */
    J[465] -= dqdci;              /* dwdot[HCOOH]/d[CO] */
    /* d()/d[HCOOH] */
    dqdci =  + k_f*sc[4];
    J[1004] -= dqdci;             /* dwdot[CH3]/d[HCOOH] */
    J[1006] += dqdci;             /* dwdot[CH4]/d[HCOOH] */
    J[1007] += dqdci;             /* dwdot[OH]/d[HCOOH] */
    J[1011] += dqdci;             /* dwdot[CO]/d[HCOOH] */
    J[1025] -= dqdci;             /* dwdot[HCOOH]/d[HCOOH] */
    /* d()/dT */
    J[1564] -= dqdT;              /* dwdot[CH3]/dT */
    J[1566] += dqdT;              /* dwdot[CH4]/dT */
    J[1567] += dqdT;              /* dwdot[OH]/dT */
    J[1571] += dqdT;              /* dwdot[CO]/dT */
    J[1585] -= dqdT;              /* dwdot[HCOOH]/dT */

    /*reaction 174: HCOOH + HO2 <=> H2O2 + CO + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[20]*sc[25];
    k_f = prefactor_units[173] * fwd_A[173]
                * exp(fwd_beta[173] * tc[0] - activation_units[173] * fwd_Ea[173] * invT);
    dlnkfdT = fwd_beta[173] * invT + activation_units[173] * fwd_Ea[173] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[11]*sc[21];
    Kc = refC * exp(-g_RT[7] - g_RT[11] + g_RT[20] - g_RT[21] + g_RT[25]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[20] + h_RT[25]) + (h_RT[7] + h_RT[11] + h_RT[21]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[7] += q; /* OH */
    wdot[11] += q; /* CO */
    wdot[20] -= q; /* HO2 */
    wdot[21] += q; /* H2O2 */
    wdot[25] -= q; /* HCOOH */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[11]*sc[21];
    J[287] += dqdci;              /* dwdot[OH]/d[OH] */
    J[291] += dqdci;              /* dwdot[CO]/d[OH] */
    J[300] -= dqdci;              /* dwdot[HO2]/d[OH] */
    J[301] += dqdci;              /* dwdot[H2O2]/d[OH] */
    J[305] -= dqdci;              /* dwdot[HCOOH]/d[OH] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[7]*sc[21];
    J[447] += dqdci;              /* dwdot[OH]/d[CO] */
    J[451] += dqdci;              /* dwdot[CO]/d[CO] */
    J[460] -= dqdci;              /* dwdot[HO2]/d[CO] */
    J[461] += dqdci;              /* dwdot[H2O2]/d[CO] */
    J[465] -= dqdci;              /* dwdot[HCOOH]/d[CO] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[25];
    J[807] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[811] += dqdci;              /* dwdot[CO]/d[HO2] */
    J[820] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[821] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    J[825] -= dqdci;              /* dwdot[HCOOH]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[7]*sc[11];
    J[847] += dqdci;              /* dwdot[OH]/d[H2O2] */
    J[851] += dqdci;              /* dwdot[CO]/d[H2O2] */
    J[860] -= dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[861] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
    J[865] -= dqdci;              /* dwdot[HCOOH]/d[H2O2] */
    /* d()/d[HCOOH] */
    dqdci =  + k_f*sc[20];
    J[1007] += dqdci;             /* dwdot[OH]/d[HCOOH] */
    J[1011] += dqdci;             /* dwdot[CO]/d[HCOOH] */
    J[1020] -= dqdci;             /* dwdot[HO2]/d[HCOOH] */
    J[1021] += dqdci;             /* dwdot[H2O2]/d[HCOOH] */
    J[1025] -= dqdci;             /* dwdot[HCOOH]/d[HCOOH] */
    /* d()/dT */
    J[1567] += dqdT;              /* dwdot[OH]/dT */
    J[1571] += dqdT;              /* dwdot[CO]/dT */
    J[1580] -= dqdT;              /* dwdot[HO2]/dT */
    J[1581] += dqdT;              /* dwdot[H2O2]/dT */
    J[1585] -= dqdT;              /* dwdot[HCOOH]/dT */

    /*reaction 175: HCOOH + O <=> CO + OH + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[25];
    k_f = prefactor_units[174] * fwd_A[174]
                * exp(fwd_beta[174] * tc[0] - activation_units[174] * fwd_Ea[174] * invT);
    dlnkfdT = fwd_beta[174] * invT + activation_units[174] * fwd_Ea[174] * invT2;
    /* reverse */
    phi_r = pow(sc[7], 2.000000)*sc[11];
    Kc = refC * exp(g_RT[5] - g_RT[7] - g_RT[7] - g_RT[11] + g_RT[25]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[25]) + (2.000000*h_RT[7] + h_RT[11]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] -= q; /* O */
    wdot[7] += 2 * q; /* OH */
    wdot[11] += q; /* CO */
    wdot[25] -= q; /* HCOOH */
    /* d()/d[O] */
    dqdci =  + k_f*sc[25];
    J[205] -= dqdci;              /* dwdot[O]/d[O] */
    J[207] += 2 * dqdci;          /* dwdot[OH]/d[O] */
    J[211] += dqdci;              /* dwdot[CO]/d[O] */
    J[225] -= dqdci;              /* dwdot[HCOOH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*2.000000*sc[7]*sc[11];
    J[285] -= dqdci;              /* dwdot[O]/d[OH] */
    J[287] += 2 * dqdci;          /* dwdot[OH]/d[OH] */
    J[291] += dqdci;              /* dwdot[CO]/d[OH] */
    J[305] -= dqdci;              /* dwdot[HCOOH]/d[OH] */
    /* d()/d[CO] */
    dqdci =  - k_r*pow(sc[7], 2.000000);
    J[445] -= dqdci;              /* dwdot[O]/d[CO] */
    J[447] += 2 * dqdci;          /* dwdot[OH]/d[CO] */
    J[451] += dqdci;              /* dwdot[CO]/d[CO] */
    J[465] -= dqdci;              /* dwdot[HCOOH]/d[CO] */
    /* d()/d[HCOOH] */
    dqdci =  + k_f*sc[5];
    J[1005] -= dqdci;             /* dwdot[O]/d[HCOOH] */
    J[1007] += 2 * dqdci;         /* dwdot[OH]/d[HCOOH] */
    J[1011] += dqdci;             /* dwdot[CO]/d[HCOOH] */
    J[1025] -= dqdci;             /* dwdot[HCOOH]/d[HCOOH] */
    /* d()/dT */
    J[1565] -= dqdT;              /* dwdot[O]/dT */
    J[1567] += 2 * dqdT;          /* dwdot[OH]/dT */
    J[1571] += dqdT;              /* dwdot[CO]/dT */
    J[1585] -= dqdT;              /* dwdot[HCOOH]/dT */

    amrex::Real c_R[39], dcRdT[39], e_RT[39];
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
    for (int k = 0; k < 39; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[1560+k];
    }

    amrex::Real cmixinv = 1.0/cmix;
    amrex::Real tmp1 = ehmix*cmixinv;
    amrex::Real tmp3 = cmixinv*T;
    amrex::Real tmp2 = tmp1*tmp3;
    amrex::Real dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 39; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 39; ++m) {
            dehmixdc += eh_RT[m]*J[k*40+m];
        }
        J[k*40+39] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[1599] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;
}

#endif
/* Transport function declarations  */


void egtransetLENIMC(int* LENIMC ) {
    *LENIMC = 158;}


void egtransetLENRMC(int* LENRMC ) {
    *LENRMC = 30264;}


void egtransetNO(int* NO ) {
    *NO = 4;}


void egtransetKK(int* KK ) {
    *KK = 39;}


void egtransetNLITE(int* NLITE ) {
    *NLITE = 2;}


/*Patm in ergs/cm3 */
void egtransetPATM(amrex::Real* PATM) {
    *PATM =   0.1013250000000000E+07;}


/*the molecular weights in g/mol */
void egtransetWT(amrex::Real* WT ) {
    WT[0] = 1.00797000E+00;
    WT[1] = 2.01594000E+00;
    WT[2] = 1.40270900E+01;
    WT[3] = 1.40270900E+01;
    WT[4] = 1.50350600E+01;
    WT[5] = 1.59994000E+01;
    WT[6] = 1.60430300E+01;
    WT[7] = 1.70073700E+01;
    WT[8] = 1.80153400E+01;
    WT[9] = 2.60382400E+01;
    WT[10] = 2.70462100E+01;
    WT[11] = 2.80105500E+01;
    WT[12] = 2.80541800E+01;
    WT[13] = 2.90185200E+01;
    WT[14] = 2.90621500E+01;
    WT[15] = 3.00264900E+01;
    WT[16] = 3.00701200E+01;
    WT[17] = 3.10344600E+01;
    WT[18] = 3.10344600E+01;
    WT[19] = 3.19988000E+01;
    WT[20] = 3.30067700E+01;
    WT[21] = 3.40147400E+01;
    WT[22] = 4.40099500E+01;
    WT[23] = 4.40535800E+01;
    WT[24] = 4.50615500E+01;
    WT[25] = 4.60258900E+01;
    WT[26] = 4.60695200E+01;
    WT[27] = 4.70338600E+01;
    WT[28] = 5.90450100E+01;
    WT[29] = 6.00529800E+01;
    WT[30] = 6.10609500E+01;
    WT[31] = 6.20689200E+01;
    WT[32] = 7.50444100E+01;
    WT[33] = 7.50444100E+01;
    WT[34] = 7.70603500E+01;
    WT[35] = 7.70603500E+01;
    WT[36] = 9.20517800E+01;
    WT[37] = 1.09059150E+02;
    WT[38] = 2.80134000E+01;
}


/*the lennard-jones potential well depth eps/kb in K */
void egtransetEPS(amrex::Real* EPS ) {
    EPS[0] = 1.45000000E+02;
    EPS[1] = 3.80000000E+01;
    EPS[2] = 1.44000000E+02;
    EPS[3] = 1.44000000E+02;
    EPS[4] = 1.44000000E+02;
    EPS[5] = 8.00000000E+01;
    EPS[6] = 1.41400000E+02;
    EPS[7] = 8.00000000E+01;
    EPS[8] = 5.72400000E+02;
    EPS[9] = 2.65300000E+02;
    EPS[10] = 2.65300000E+02;
    EPS[11] = 9.81000000E+01;
    EPS[12] = 2.38400000E+02;
    EPS[13] = 4.98000000E+02;
    EPS[14] = 2.47500000E+02;
    EPS[15] = 4.98000000E+02;
    EPS[16] = 2.47500000E+02;
    EPS[17] = 4.17000000E+02;
    EPS[18] = 4.17000000E+02;
    EPS[19] = 1.07400000E+02;
    EPS[20] = 1.07400000E+02;
    EPS[21] = 1.07400000E+02;
    EPS[22] = 2.44000000E+02;
    EPS[23] = 4.36000000E+02;
    EPS[24] = 3.29400000E+02;
    EPS[25] = 4.70600000E+02;
    EPS[26] = 3.29400000E+02;
    EPS[27] = 4.81800000E+02;
    EPS[28] = 4.06500000E+02;
    EPS[29] = 4.06500000E+02;
    EPS[30] = 4.70900000E+02;
    EPS[31] = 3.29400000E+02;
    EPS[32] = 3.29400000E+02;
    EPS[33] = 3.29400000E+02;
    EPS[34] = 3.29400000E+02;
    EPS[35] = 3.29400000E+02;
    EPS[36] = 3.29400000E+02;
    EPS[37] = 3.29400000E+02;
    EPS[38] = 9.75300000E+01;
}


/*the lennard-jones collision diameter in Angstroms */
void egtransetSIG(amrex::Real* SIG ) {
    SIG[0] = 2.05000000E+00;
    SIG[1] = 2.92000000E+00;
    SIG[2] = 3.80000000E+00;
    SIG[3] = 3.80000000E+00;
    SIG[4] = 3.80000000E+00;
    SIG[5] = 2.75000000E+00;
    SIG[6] = 3.74600000E+00;
    SIG[7] = 2.75000000E+00;
    SIG[8] = 2.60500000E+00;
    SIG[9] = 3.72100000E+00;
    SIG[10] = 3.72100000E+00;
    SIG[11] = 3.65000000E+00;
    SIG[12] = 3.49600000E+00;
    SIG[13] = 3.59000000E+00;
    SIG[14] = 4.35000000E+00;
    SIG[15] = 3.59000000E+00;
    SIG[16] = 4.35000000E+00;
    SIG[17] = 3.69000000E+00;
    SIG[18] = 3.69000000E+00;
    SIG[19] = 3.45800000E+00;
    SIG[20] = 3.45800000E+00;
    SIG[21] = 3.45800000E+00;
    SIG[22] = 3.76300000E+00;
    SIG[23] = 3.97000000E+00;
    SIG[24] = 4.62400000E+00;
    SIG[25] = 4.41000000E+00;
    SIG[26] = 4.62400000E+00;
    SIG[27] = 3.62600000E+00;
    SIG[28] = 4.70900000E+00;
    SIG[29] = 4.70900000E+00;
    SIG[30] = 4.86200000E+00;
    SIG[31] = 4.62400000E+00;
    SIG[32] = 4.62400000E+00;
    SIG[33] = 4.62400000E+00;
    SIG[34] = 4.62400000E+00;
    SIG[35] = 4.62400000E+00;
    SIG[36] = 4.62400000E+00;
    SIG[37] = 4.62400000E+00;
    SIG[38] = 3.62100000E+00;
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
    DIP[8] = 1.84400000E+00;
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
    DIP[21] = 0.00000000E+00;
    DIP[22] = 0.00000000E+00;
    DIP[23] = 0.00000000E+00;
    DIP[24] = 0.00000000E+00;
    DIP[25] = 0.00000000E+00;
    DIP[26] = 0.00000000E+00;
    DIP[27] = 0.00000000E+00;
    DIP[28] = 0.00000000E+00;
    DIP[29] = 0.00000000E+00;
    DIP[30] = 0.00000000E+00;
    DIP[31] = 0.00000000E+00;
    DIP[32] = 0.00000000E+00;
    DIP[33] = 0.00000000E+00;
    DIP[34] = 0.00000000E+00;
    DIP[35] = 0.00000000E+00;
    DIP[36] = 0.00000000E+00;
    DIP[37] = 0.00000000E+00;
    DIP[38] = 0.00000000E+00;
}


/*the polarizability in cubic Angstroms */
void egtransetPOL(amrex::Real* POL ) {
    POL[0] = 0.00000000E+00;
    POL[1] = 7.90000000E-01;
    POL[2] = 0.00000000E+00;
    POL[3] = 0.00000000E+00;
    POL[4] = 0.00000000E+00;
    POL[5] = 0.00000000E+00;
    POL[6] = 2.60000000E+00;
    POL[7] = 0.00000000E+00;
    POL[8] = 0.00000000E+00;
    POL[9] = 0.00000000E+00;
    POL[10] = 0.00000000E+00;
    POL[11] = 1.95000000E+00;
    POL[12] = 0.00000000E+00;
    POL[13] = 0.00000000E+00;
    POL[14] = 0.00000000E+00;
    POL[15] = 0.00000000E+00;
    POL[16] = 0.00000000E+00;
    POL[17] = 0.00000000E+00;
    POL[18] = 0.00000000E+00;
    POL[19] = 1.60000000E+00;
    POL[20] = 0.00000000E+00;
    POL[21] = 0.00000000E+00;
    POL[22] = 2.65000000E+00;
    POL[23] = 0.00000000E+00;
    POL[24] = 0.00000000E+00;
    POL[25] = 0.00000000E+00;
    POL[26] = 0.00000000E+00;
    POL[27] = 0.00000000E+00;
    POL[28] = 0.00000000E+00;
    POL[29] = 0.00000000E+00;
    POL[30] = 0.00000000E+00;
    POL[31] = 0.00000000E+00;
    POL[32] = 0.00000000E+00;
    POL[33] = 0.00000000E+00;
    POL[34] = 0.00000000E+00;
    POL[35] = 0.00000000E+00;
    POL[36] = 0.00000000E+00;
    POL[37] = 0.00000000E+00;
    POL[38] = 1.76000000E+00;
}


/*the rotational relaxation collision number at 298 K */
void egtransetZROT(amrex::Real* ZROT ) {
    ZROT[0] = 0.00000000E+00;
    ZROT[1] = 2.80000000E+02;
    ZROT[2] = 0.00000000E+00;
    ZROT[3] = 0.00000000E+00;
    ZROT[4] = 0.00000000E+00;
    ZROT[5] = 0.00000000E+00;
    ZROT[6] = 1.30000000E+01;
    ZROT[7] = 0.00000000E+00;
    ZROT[8] = 4.00000000E+00;
    ZROT[9] = 2.50000000E+00;
    ZROT[10] = 1.00000000E+00;
    ZROT[11] = 1.80000000E+00;
    ZROT[12] = 1.50000000E+00;
    ZROT[13] = 0.00000000E+00;
    ZROT[14] = 1.50000000E+00;
    ZROT[15] = 2.00000000E+00;
    ZROT[16] = 1.50000000E+00;
    ZROT[17] = 2.00000000E+00;
    ZROT[18] = 2.00000000E+00;
    ZROT[19] = 3.80000000E+00;
    ZROT[20] = 1.00000000E+00;
    ZROT[21] = 3.80000000E+00;
    ZROT[22] = 2.10000000E+00;
    ZROT[23] = 2.00000000E+00;
    ZROT[24] = 1.00000000E+00;
    ZROT[25] = 1.50000000E+00;
    ZROT[26] = 1.00000000E+00;
    ZROT[27] = 1.00000000E+00;
    ZROT[28] = 1.00000000E+00;
    ZROT[29] = 1.00000000E+00;
    ZROT[30] = 1.00000000E+00;
    ZROT[31] = 1.00000000E+00;
    ZROT[32] = 1.00000000E+00;
    ZROT[33] = 1.00000000E+00;
    ZROT[34] = 1.00000000E+00;
    ZROT[35] = 1.00000000E+00;
    ZROT[36] = 1.00000000E+00;
    ZROT[37] = 1.00000000E+00;
    ZROT[38] = 4.00000000E+00;
}


/*0: monoatomic, 1: linear, 2: nonlinear */
void egtransetNLIN(int* NLIN) {
    NLIN[0] = 0;
    NLIN[1] = 1;
    NLIN[2] = 1;
    NLIN[3] = 1;
    NLIN[4] = 1;
    NLIN[5] = 0;
    NLIN[6] = 2;
    NLIN[7] = 1;
    NLIN[8] = 2;
    NLIN[9] = 1;
    NLIN[10] = 2;
    NLIN[11] = 1;
    NLIN[12] = 2;
    NLIN[13] = 2;
    NLIN[14] = 2;
    NLIN[15] = 2;
    NLIN[16] = 2;
    NLIN[17] = 2;
    NLIN[18] = 2;
    NLIN[19] = 1;
    NLIN[20] = 2;
    NLIN[21] = 2;
    NLIN[22] = 1;
    NLIN[23] = 2;
    NLIN[24] = 2;
    NLIN[25] = 2;
    NLIN[26] = 2;
    NLIN[27] = 2;
    NLIN[28] = 2;
    NLIN[29] = 2;
    NLIN[30] = 2;
    NLIN[31] = 2;
    NLIN[32] = 2;
    NLIN[33] = 2;
    NLIN[34] = 2;
    NLIN[35] = 2;
    NLIN[36] = 2;
    NLIN[37] = 2;
    NLIN[38] = 1;
}


/*Poly fits for the viscosities, dim NO*KK */
void egtransetCOFETA(amrex::Real* COFETA) {
    COFETA[0] = -2.04078397E+01;
    COFETA[1] = 3.65436395E+00;
    COFETA[2] = -3.98339635E-01;
    COFETA[3] = 1.75883009E-02;
    COFETA[4] = -1.38347699E+01;
    COFETA[5] = 1.00106621E+00;
    COFETA[6] = -4.98105694E-02;
    COFETA[7] = 2.31450475E-03;
    COFETA[8] = -2.02663469E+01;
    COFETA[9] = 3.63241793E+00;
    COFETA[10] = -3.95581049E-01;
    COFETA[11] = 1.74725495E-02;
    COFETA[12] = -2.02663469E+01;
    COFETA[13] = 3.63241793E+00;
    COFETA[14] = -3.95581049E-01;
    COFETA[15] = 1.74725495E-02;
    COFETA[16] = -2.02316497E+01;
    COFETA[17] = 3.63241793E+00;
    COFETA[18] = -3.95581049E-01;
    COFETA[19] = 1.74725495E-02;
    COFETA[20] = -1.50926240E+01;
    COFETA[21] = 1.92606504E+00;
    COFETA[22] = -1.73487476E-01;
    COFETA[23] = 7.82572931E-03;
    COFETA[24] = -2.00094664E+01;
    COFETA[25] = 3.57220167E+00;
    COFETA[26] = -3.87936446E-01;
    COFETA[27] = 1.71483254E-02;
    COFETA[28] = -1.50620763E+01;
    COFETA[29] = 1.92606504E+00;
    COFETA[30] = -1.73487476E-01;
    COFETA[31] = 7.82572931E-03;
    COFETA[32] = -1.05420863E+01;
    COFETA[33] = -1.37777096E+00;
    COFETA[34] = 4.20502308E-01;
    COFETA[35] = -2.40627230E-02;
    COFETA[36] = -2.47697856E+01;
    COFETA[37] = 5.30039568E+00;
    COFETA[38] = -5.89273639E-01;
    COFETA[39] = 2.49261407E-02;
    COFETA[40] = -2.47507953E+01;
    COFETA[41] = 5.30039568E+00;
    COFETA[42] = -5.89273639E-01;
    COFETA[43] = 2.49261407E-02;
    COFETA[44] = -1.66188336E+01;
    COFETA[45] = 2.40307799E+00;
    COFETA[46] = -2.36167638E-01;
    COFETA[47] = 1.05714061E-02;
    COFETA[48] = -2.39690472E+01;
    COFETA[49] = 5.11436059E+00;
    COFETA[50] = -5.71999954E-01;
    COFETA[51] = 2.44581334E-02;
    COFETA[52] = -1.98501306E+01;
    COFETA[53] = 2.69480162E+00;
    COFETA[54] = -1.65880845E-01;
    COFETA[55] = 3.14504769E-03;
    COFETA[56] = -2.45602637E+01;
    COFETA[57] = 5.15878990E+00;
    COFETA[58] = -5.75274341E-01;
    COFETA[59] = 2.44975136E-02;
    COFETA[60] = -1.98330577E+01;
    COFETA[61] = 2.69480162E+00;
    COFETA[62] = -1.65880845E-01;
    COFETA[63] = 3.14504769E-03;
    COFETA[64] = -2.45432160E+01;
    COFETA[65] = 5.15878990E+00;
    COFETA[66] = -5.75274341E-01;
    COFETA[67] = 2.44975136E-02;
    COFETA[68] = -1.99945919E+01;
    COFETA[69] = 2.86923313E+00;
    COFETA[70] = -2.03325661E-01;
    COFETA[71] = 5.39056989E-03;
    COFETA[72] = -1.99945919E+01;
    COFETA[73] = 2.86923313E+00;
    COFETA[74] = -2.03325661E-01;
    COFETA[75] = 5.39056989E-03;
    COFETA[76] = -1.71618309E+01;
    COFETA[77] = 2.68036374E+00;
    COFETA[78] = -2.72570227E-01;
    COFETA[79] = 1.21650964E-02;
    COFETA[80] = -1.71463238E+01;
    COFETA[81] = 2.68036374E+00;
    COFETA[82] = -2.72570227E-01;
    COFETA[83] = 1.21650964E-02;
    COFETA[84] = -1.71312832E+01;
    COFETA[85] = 2.68036374E+00;
    COFETA[86] = -2.72570227E-01;
    COFETA[87] = 1.21650964E-02;
    COFETA[88] = -2.40014975E+01;
    COFETA[89] = 5.14359547E+00;
    COFETA[90] = -5.74269731E-01;
    COFETA[91] = 2.44937679E-02;
    COFETA[92] = -2.23161441E+01;
    COFETA[93] = 3.86433912E+00;
    COFETA[94] = -3.41553983E-01;
    COFETA[95] = 1.17083447E-02;
    COFETA[96] = -2.50662656E+01;
    COFETA[97] = 5.17000432E+00;
    COFETA[98] = -5.51481603E-01;
    COFETA[99] = 2.24387705E-02;
    COFETA[100] = -2.12394850E+01;
    COFETA[101] = 3.25726898E+00;
    COFETA[102] = -2.49519605E-01;
    COFETA[103] = 7.19215196E-03;
    COFETA[104] = -2.50552045E+01;
    COFETA[105] = 5.17000432E+00;
    COFETA[106] = -5.51481603E-01;
    COFETA[107] = 2.24387705E-02;
    COFETA[108] = -2.03725491E+01;
    COFETA[109] = 3.03946431E+00;
    COFETA[110] = -2.16994867E-01;
    COFETA[111] = 5.61394012E-03;
    COFETA[112] = -2.34833807E+01;
    COFETA[113] = 4.34232030E+00;
    COFETA[114] = -4.15133257E-01;
    COFETA[115] = 1.53602178E-02;
    COFETA[116] = -2.34749172E+01;
    COFETA[117] = 4.34232030E+00;
    COFETA[118] = -4.15133257E-01;
    COFETA[119] = 1.53602178E-02;
    COFETA[120] = -2.12812276E+01;
    COFETA[121] = 3.25157556E+00;
    COFETA[122] = -2.48665998E-01;
    COFETA[123] = 7.15060391E-03;
    COFETA[124] = -2.49061576E+01;
    COFETA[125] = 5.17000432E+00;
    COFETA[126] = -5.51481603E-01;
    COFETA[127] = 2.24387705E-02;
    COFETA[128] = -2.48112403E+01;
    COFETA[129] = 5.17000432E+00;
    COFETA[130] = -5.51481603E-01;
    COFETA[131] = 2.24387705E-02;
    COFETA[132] = -2.48112403E+01;
    COFETA[133] = 5.17000432E+00;
    COFETA[134] = -5.51481603E-01;
    COFETA[135] = 2.24387705E-02;
    COFETA[136] = -2.47979859E+01;
    COFETA[137] = 5.17000432E+00;
    COFETA[138] = -5.51481603E-01;
    COFETA[139] = 2.24387705E-02;
    COFETA[140] = -2.47979859E+01;
    COFETA[141] = 5.17000432E+00;
    COFETA[142] = -5.51481603E-01;
    COFETA[143] = 2.24387705E-02;
    COFETA[144] = -2.47091047E+01;
    COFETA[145] = 5.17000432E+00;
    COFETA[146] = -5.51481603E-01;
    COFETA[147] = 2.24387705E-02;
    COFETA[148] = -2.46243351E+01;
    COFETA[149] = 5.17000432E+00;
    COFETA[150] = -5.51481603E-01;
    COFETA[151] = 2.24387705E-02;
    COFETA[152] = -1.65695594E+01;
    COFETA[153] = 2.39056562E+00;
    COFETA[154] = -2.34558144E-01;
    COFETA[155] = 1.05024037E-02;
}


/*Poly fits for the conductivities, dim NO*KK */
void egtransetCOFLAM(amrex::Real* COFLAM) {
    COFLAM[0] = -8.57929284E-01;
    COFLAM[1] = 3.65436395E+00;
    COFLAM[2] = -3.98339635E-01;
    COFLAM[3] = 1.75883009E-02;
    COFLAM[4] = 9.13734407E+00;
    COFLAM[5] = -4.36833375E-01;
    COFLAM[6] = 1.12981765E-01;
    COFLAM[7] = -2.54610412E-03;
    COFLAM[8] = 1.29177902E+01;
    COFLAM[9] = -3.73745535E+00;
    COFLAM[10] = 7.15831021E-01;
    COFLAM[11] = -3.63846910E-02;
    COFLAM[12] = 1.89383266E+01;
    COFLAM[13] = -6.51018128E+00;
    COFLAM[14] = 1.13292060E+00;
    COFLAM[15] = -5.69603226E-02;
    COFLAM[16] = 1.32973935E+01;
    COFLAM[17] = -4.31770522E+00;
    COFLAM[18] = 8.57562813E-01;
    COFLAM[19] = -4.51644164E-02;
    COFLAM[20] = 1.69267361E+00;
    COFLAM[21] = 1.92606504E+00;
    COFLAM[22] = -1.73487476E-01;
    COFLAM[23] = 7.82572931E-03;
    COFLAM[24] = 1.29622480E+01;
    COFLAM[25] = -4.85747192E+00;
    COFLAM[26] = 1.02918185E+00;
    COFLAM[27] = -5.69931976E-02;
    COFLAM[28] = 1.50119731E+01;
    COFLAM[29] = -3.63267854E+00;
    COFLAM[30] = 5.92839101E-01;
    COFLAM[31] = -2.62920439E-02;
    COFLAM[32] = 2.35086340E+01;
    COFLAM[33] = -9.05997330E+00;
    COFLAM[34] = 1.54816025E+00;
    COFLAM[35] = -7.71639384E-02;
    COFLAM[36] = -9.07948800E+00;
    COFLAM[37] = 5.06981141E+00;
    COFLAM[38] = -4.58219178E-01;
    COFLAM[39] = 1.59367685E-02;
    COFLAM[40] = -1.04615545E+01;
    COFLAM[41] = 5.04178845E+00;
    COFLAM[42] = -3.70776164E-01;
    COFLAM[43] = 8.42040125E-03;
    COFLAM[44] = 1.15794127E+01;
    COFLAM[45] = -3.02088602E+00;
    COFLAM[46] = 5.82091962E-01;
    COFLAM[47] = -2.93406692E-02;
    COFLAM[48] = -1.34346703E+01;
    COFLAM[49] = 6.08825042E+00;
    COFLAM[50] = -4.77261530E-01;
    COFLAM[51] = 1.18191214E-02;
    COFLAM[52] = 5.62029187E+00;
    COFLAM[53] = -1.91574800E+00;
    COFLAM[54] = 5.90225152E-01;
    COFLAM[55] = -3.57616080E-02;
    COFLAM[56] = -9.76065737E+00;
    COFLAM[57] = 4.41719885E+00;
    COFLAM[58] = -2.47339046E-01;
    COFLAM[59] = 1.40174222E-03;
    COFLAM[60] = 3.57128442E+00;
    COFLAM[61] = -1.58527059E+00;
    COFLAM[62] = 6.20991255E-01;
    COFLAM[63] = -4.00762909E-02;
    COFLAM[64] = -1.39341971E+01;
    COFLAM[65] = 6.08925205E+00;
    COFLAM[66] = -4.67702719E-01;
    COFLAM[67] = 1.12481374E-02;
    COFLAM[68] = -5.56720763E+00;
    COFLAM[69] = 2.89052547E+00;
    COFLAM[70] = -7.30059663E-02;
    COFLAM[71] = -5.12272887E-03;
    COFLAM[72] = -6.14587484E+00;
    COFLAM[73] = 2.47428533E+00;
    COFLAM[74] = 6.44004931E-02;
    COFLAM[75] = -1.45368612E-02;
    COFLAM[76] = -1.93888086E+00;
    COFLAM[77] = 2.89244157E+00;
    COFLAM[78] = -2.71258557E-01;
    COFLAM[79] = 1.15340544E-02;
    COFLAM[80] = -1.12960913E+00;
    COFLAM[81] = 2.34014305E+00;
    COFLAM[82] = -1.63245030E-01;
    COFLAM[83] = 5.80319600E-03;
    COFLAM[84] = 8.99254403E-01;
    COFLAM[85] = 1.32506478E+00;
    COFLAM[86] = 1.81955930E-02;
    COFLAM[87] = -4.46691285E-03;
    COFLAM[88] = -1.15552013E+01;
    COFLAM[89] = 5.97444378E+00;
    COFLAM[90] = -5.83493959E-01;
    COFLAM[91] = 2.11390997E-02;
    COFLAM[92] = -9.88397056E+00;
    COFLAM[93] = 4.13891036E+00;
    COFLAM[94] = -1.74918892E-01;
    COFLAM[95] = -3.28392643E-03;
    COFLAM[96] = -5.98455417E+00;
    COFLAM[97] = 2.74597482E+00;
    COFLAM[98] = -1.33426244E-02;
    COFLAM[99] = -9.61096027E-03;
    COFLAM[100] = -1.36375211E+01;
    COFLAM[101] = 5.70113661E+00;
    COFLAM[102] = -4.14213135E-01;
    COFLAM[103] = 8.34373882E-03;
    COFLAM[104] = -8.31367313E+00;
    COFLAM[105] = 3.45763526E+00;
    COFLAM[106] = -7.71003858E-02;
    COFLAM[107] = -8.00162224E-03;
    COFLAM[108] = 1.68422913E+00;
    COFLAM[109] = -6.01950082E-01;
    COFLAM[110] = 4.68956133E-01;
    COFLAM[111] = -3.23944583E-02;
    COFLAM[112] = -1.79264269E+01;
    COFLAM[113] = 8.09034865E+00;
    COFLAM[114] = -8.01396311E-01;
    COFLAM[115] = 2.80655203E-02;
    COFLAM[116] = -8.64004520E+00;
    COFLAM[117] = 3.57815725E+00;
    COFLAM[118] = -1.01159505E-01;
    COFLAM[119] = -6.85522605E-03;
    COFLAM[120] = -6.07331497E+00;
    COFLAM[121] = 2.38074449E+00;
    COFLAM[122] = 7.42085915E-02;
    COFLAM[123] = -1.50407494E-02;
    COFLAM[124] = -8.20670492E+00;
    COFLAM[125] = 3.61610027E+00;
    COFLAM[126] = -1.24905892E-01;
    COFLAM[127] = -4.79916006E-03;
    COFLAM[128] = -9.03345078E+00;
    COFLAM[129] = 4.13214827E+00;
    COFLAM[130] = -2.21155748E-01;
    COFLAM[131] = 1.92816462E-04;
    COFLAM[132] = -4.75221469E+00;
    COFLAM[133] = 2.39466705E+00;
    COFLAM[134] = 1.10520831E-02;
    COFLAM[135] = -1.00968213E-02;
    COFLAM[136] = -1.49512219E+01;
    COFLAM[137] = 6.75380005E+00;
    COFLAM[138] = -6.00471360E-01;
    COFLAM[139] = 1.86094829E-02;
    COFLAM[140] = -1.90602373E+01;
    COFLAM[141] = 8.74377154E+00;
    COFLAM[142] = -9.05452607E-01;
    COFLAM[143] = 3.36467064E-02;
    COFLAM[144] = -1.84034518E+01;
    COFLAM[145] = 8.43437025E+00;
    COFLAM[146] = -8.58775267E-01;
    COFLAM[147] = 3.12042072E-02;
    COFLAM[148] = -2.24994186E+01;
    COFLAM[149] = 1.02734334E+01;
    COFLAM[150] = -1.12700167E+00;
    COFLAM[151] = 4.41962434E-02;
    COFLAM[152] = 1.29306625E+01;
    COFLAM[153] = -3.52819393E+00;
    COFLAM[154] = 6.45501923E-01;
    COFLAM[155] = -3.19376675E-02;
}


/*Poly fits for the diffusion coefficients, dim NO*KK*KK */
void egtransetCOFD(amrex::Real* COFD) {
    COFD[0] = -1.47968712E+01;
    COFD[1] = 4.23027636E+00;
    COFD[2] = -3.36139991E-01;
    COFD[3] = 1.46507621E-02;
    COFD[4] = -1.14366381E+01;
    COFD[5] = 2.78323501E+00;
    COFD[6] = -1.51214064E-01;
    COFD[7] = 6.75150012E-03;
    COFD[8] = -1.57972369E+01;
    COFD[9] = 4.22225052E+00;
    COFD[10] = -3.35156428E-01;
    COFD[11] = 1.46104855E-02;
    COFD[12] = -1.57972369E+01;
    COFD[13] = 4.22225052E+00;
    COFD[14] = -3.35156428E-01;
    COFD[15] = 1.46104855E-02;
    COFD[16] = -1.57994893E+01;
    COFD[17] = 4.22225052E+00;
    COFD[18] = -3.35156428E-01;
    COFD[19] = 1.46104855E-02;
    COFD[20] = -1.34230272E+01;
    COFD[21] = 3.48624238E+00;
    COFD[22] = -2.41554467E-01;
    COFD[23] = 1.06263545E-02;
    COFD[24] = -1.57199037E+01;
    COFD[25] = 4.19936335E+00;
    COFD[26] = -3.32311009E-01;
    COFD[27] = 1.44921003E-02;
    COFD[28] = -1.34247866E+01;
    COFD[29] = 3.48624238E+00;
    COFD[30] = -2.41554467E-01;
    COFD[31] = 1.06263545E-02;
    COFD[32] = -1.95739570E+01;
    COFD[33] = 5.61113230E+00;
    COFD[34] = -4.90190187E-01;
    COFD[35] = 2.03260675E-02;
    COFD[36] = -1.79310765E+01;
    COFD[37] = 4.98037650E+00;
    COFD[38] = -4.26676911E-01;
    COFD[39] = 1.83007231E-02;
    COFD[40] = -1.79317714E+01;
    COFD[41] = 4.98037650E+00;
    COFD[42] = -4.26676911E-01;
    COFD[43] = 1.83007231E-02;
    COFD[44] = -1.43151174E+01;
    COFD[45] = 3.68038508E+00;
    COFD[46] = -2.65779346E-01;
    COFD[47] = 1.16360771E-02;
    COFD[48] = -1.74407963E+01;
    COFD[49] = 4.83580036E+00;
    COFD[50] = -4.09383573E-01;
    COFD[51] = 1.76098175E-02;
    COFD[52] = -1.97544450E+01;
    COFD[53] = 5.56931926E+00;
    COFD[54] = -4.89105511E-01;
    COFD[55] = 2.04493129E-02;
    COFD[56] = -1.78631557E+01;
    COFD[57] = 4.88268692E+00;
    COFD[58] = -4.14917638E-01;
    COFD[59] = 1.78274298E-02;
    COFD[60] = -1.97550088E+01;
    COFD[61] = 5.56931926E+00;
    COFD[62] = -4.89105511E-01;
    COFD[63] = 2.04493129E-02;
    COFD[64] = -1.78637178E+01;
    COFD[65] = 4.88268692E+00;
    COFD[66] = -4.14917638E-01;
    COFD[67] = 1.78274298E-02;
    COFD[68] = -1.92718582E+01;
    COFD[69] = 5.41172124E+00;
    COFD[70] = -4.73213887E-01;
    COFD[71] = 1.99405473E-02;
    COFD[72] = -1.92718582E+01;
    COFD[73] = 5.41172124E+00;
    COFD[74] = -4.73213887E-01;
    COFD[75] = 1.99405473E-02;
    COFD[76] = -1.46550083E+01;
    COFD[77] = 3.83606243E+00;
    COFD[78] = -2.86076532E-01;
    COFD[79] = 1.25205829E-02;
    COFD[80] = -1.46554748E+01;
    COFD[81] = 3.83606243E+00;
    COFD[82] = -2.86076532E-01;
    COFD[83] = 1.25205829E-02;
    COFD[84] = -1.46559141E+01;
    COFD[85] = 3.83606243E+00;
    COFD[86] = -2.86076532E-01;
    COFD[87] = 1.25205829E-02;
    COFD[88] = -1.76147026E+01;
    COFD[89] = 4.86049500E+00;
    COFD[90] = -4.12200578E-01;
    COFD[91] = 1.77160971E-02;
    COFD[92] = -1.94694048E+01;
    COFD[93] = 5.43830787E+00;
    COFD[94] = -4.75472880E-01;
    COFD[95] = 1.99909996E-02;
    COFD[96] = -1.89222647E+01;
    COFD[97] = 5.20981862E+00;
    COFD[98] = -4.52489825E-01;
    COFD[99] = 1.92609226E-02;
    COFD[100] = -1.98808664E+01;
    COFD[101] = 5.52555673E+00;
    COFD[102] = -4.84999851E-01;
    COFD[103] = 2.03334931E-02;
    COFD[104] = -1.89225041E+01;
    COFD[105] = 5.20981862E+00;
    COFD[106] = -4.52489825E-01;
    COFD[107] = 1.92609226E-02;
    COFD[108] = -1.96914944E+01;
    COFD[109] = 5.54637286E+00;
    COFD[110] = -4.87070324E-01;
    COFD[111] = 2.03983467E-02;
    COFD[112] = -1.95591756E+01;
    COFD[113] = 5.40012228E+00;
    COFD[114] = -4.72416152E-01;
    COFD[115] = 1.99342987E-02;
    COFD[116] = -1.95593165E+01;
    COFD[117] = 5.40012228E+00;
    COFD[118] = -4.72416152E-01;
    COFD[119] = 1.99342987E-02;
    COFD[120] = -2.00206731E+01;
    COFD[121] = 5.52613176E+00;
    COFD[122] = -4.85057412E-01;
    COFD[123] = 2.03353154E-02;
    COFD[124] = -1.89252713E+01;
    COFD[125] = 5.20981862E+00;
    COFD[126] = -4.52489825E-01;
    COFD[127] = 1.92609226E-02;
    COFD[128] = -1.89266547E+01;
    COFD[129] = 5.20981862E+00;
    COFD[130] = -4.52489825E-01;
    COFD[131] = 1.92609226E-02;
    COFD[132] = -1.89266547E+01;
    COFD[133] = 5.20981862E+00;
    COFD[134] = -4.52489825E-01;
    COFD[135] = 1.92609226E-02;
    COFD[136] = -1.89268281E+01;
    COFD[137] = 5.20981862E+00;
    COFD[138] = -4.52489825E-01;
    COFD[139] = 1.92609226E-02;
    COFD[140] = -1.89268281E+01;
    COFD[141] = 5.20981862E+00;
    COFD[142] = -4.52489825E-01;
    COFD[143] = 1.92609226E-02;
    COFD[144] = -1.89278805E+01;
    COFD[145] = 5.20981862E+00;
    COFD[146] = -4.52489825E-01;
    COFD[147] = 1.92609226E-02;
    COFD[148] = -1.89287258E+01;
    COFD[149] = 5.20981862E+00;
    COFD[150] = -4.52489825E-01;
    COFD[151] = 1.92609226E-02;
    COFD[152] = -1.42894441E+01;
    COFD[153] = 3.67490723E+00;
    COFD[154] = -2.65114792E-01;
    COFD[155] = 1.16092671E-02;
    COFD[156] = -1.14366381E+01;
    COFD[157] = 2.78323501E+00;
    COFD[158] = -1.51214064E-01;
    COFD[159] = 6.75150012E-03;
    COFD[160] = -1.03270606E+01;
    COFD[161] = 2.19285409E+00;
    COFD[162] = -7.54492786E-02;
    COFD[163] = 3.51398213E-03;
    COFD[164] = -1.25098960E+01;
    COFD[165] = 2.77873601E+00;
    COFD[166] = -1.50637360E-01;
    COFD[167] = 6.72684281E-03;
    COFD[168] = -1.25098960E+01;
    COFD[169] = 2.77873601E+00;
    COFD[170] = -1.50637360E-01;
    COFD[171] = 6.72684281E-03;
    COFD[172] = -1.25141260E+01;
    COFD[173] = 2.77873601E+00;
    COFD[174] = -1.50637360E-01;
    COFD[175] = 6.72684281E-03;
    COFD[176] = -1.09595712E+01;
    COFD[177] = 2.30836460E+00;
    COFD[178] = -8.76339315E-02;
    COFD[179] = 3.90878445E-03;
    COFD[180] = -1.24693568E+01;
    COFD[181] = 2.76686648E+00;
    COFD[182] = -1.49120141E-01;
    COFD[183] = 6.66220432E-03;
    COFD[184] = -1.09628982E+01;
    COFD[185] = 2.30836460E+00;
    COFD[186] = -8.76339315E-02;
    COFD[187] = 3.90878445E-03;
    COFD[188] = -1.71982995E+01;
    COFD[189] = 4.63881404E+00;
    COFD[190] = -3.86139633E-01;
    COFD[191] = 1.66955081E-02;
    COFD[192] = -1.39315266E+01;
    COFD[193] = 3.30394764E+00;
    COFD[194] = -2.17920112E-01;
    COFD[195] = 9.60284243E-03;
    COFD[196] = -1.39328674E+01;
    COFD[197] = 3.30394764E+00;
    COFD[198] = -2.17920112E-01;
    COFD[199] = 9.60284243E-03;
    COFD[200] = -1.17159737E+01;
    COFD[201] = 2.48123210E+00;
    COFD[202] = -1.11322604E-01;
    COFD[203] = 4.99282389E-03;
    COFD[204] = -1.36336373E+01;
    COFD[205] = 3.22088176E+00;
    COFD[206] = -2.07623790E-01;
    COFD[207] = 9.17771542E-03;
    COFD[208] = -1.60517370E+01;
    COFD[209] = 4.11188603E+00;
    COFD[210] = -3.21540884E-01;
    COFD[211] = 1.40482564E-02;
    COFD[212] = -1.39648112E+01;
    COFD[213] = 3.24966086E+00;
    COFD[214] = -2.11199992E-01;
    COFD[215] = 9.32580661E-03;
    COFD[216] = -1.60528285E+01;
    COFD[217] = 4.11188603E+00;
    COFD[218] = -3.21540884E-01;
    COFD[219] = 1.40482564E-02;
    COFD[220] = -1.39658996E+01;
    COFD[221] = 3.24966086E+00;
    COFD[222] = -2.11199992E-01;
    COFD[223] = 9.32580661E-03;
    COFD[224] = -1.58456300E+01;
    COFD[225] = 4.02074783E+00;
    COFD[226] = -3.10018522E-01;
    COFD[227] = 1.35599552E-02;
    COFD[228] = -1.58456300E+01;
    COFD[229] = 4.02074783E+00;
    COFD[230] = -3.10018522E-01;
    COFD[231] = 1.35599552E-02;
    COFD[232] = -1.18988955E+01;
    COFD[233] = 2.57507000E+00;
    COFD[234] = -1.24033737E-01;
    COFD[235] = 5.56694959E-03;
    COFD[236] = -1.18998012E+01;
    COFD[237] = 2.57507000E+00;
    COFD[238] = -1.24033737E-01;
    COFD[239] = 5.56694959E-03;
    COFD[240] = -1.19006548E+01;
    COFD[241] = 2.57507000E+00;
    COFD[242] = -1.24033737E-01;
    COFD[243] = 5.56694959E-03;
    COFD[244] = -1.37794315E+01;
    COFD[245] = 3.23973858E+00;
    COFD[246] = -2.09989036E-01;
    COFD[247] = 9.27667906E-03;
    COFD[248] = -1.57045332E+01;
    COFD[249] = 3.93614244E+00;
    COFD[250] = -2.99111497E-01;
    COFD[251] = 1.30888229E-02;
    COFD[252] = -1.49070584E+01;
    COFD[253] = 3.56974825E+00;
    COFD[254] = -2.52221138E-01;
    COFD[255] = 1.10819767E-02;
    COFD[256] = -1.61121174E+01;
    COFD[257] = 4.04227735E+00;
    COFD[258] = -3.12745253E-01;
    COFD[259] = 1.36756977E-02;
    COFD[260] = -1.49075271E+01;
    COFD[261] = 3.56974825E+00;
    COFD[262] = -2.52221138E-01;
    COFD[263] = 1.10819767E-02;
    COFD[264] = -1.59632479E+01;
    COFD[265] = 4.07051484E+00;
    COFD[266] = -3.16303109E-01;
    COFD[267] = 1.38259377E-02;
    COFD[268] = -1.56146754E+01;
    COFD[269] = 3.82176300E+00;
    COFD[270] = -2.84202467E-01;
    COFD[271] = 1.24384328E-02;
    COFD[272] = -1.56149526E+01;
    COFD[273] = 3.82176300E+00;
    COFD[274] = -2.84202467E-01;
    COFD[275] = 1.24384328E-02;
    COFD[276] = -1.62390700E+01;
    COFD[277] = 4.04304103E+00;
    COFD[278] = -3.12841265E-01;
    COFD[279] = 1.36797409E-02;
    COFD[280] = -1.49129599E+01;
    COFD[281] = 3.56974825E+00;
    COFD[282] = -2.52221138E-01;
    COFD[283] = 1.10819767E-02;
    COFD[284] = -1.49156868E+01;
    COFD[285] = 3.56974825E+00;
    COFD[286] = -2.52221138E-01;
    COFD[287] = 1.10819767E-02;
    COFD[288] = -1.49156868E+01;
    COFD[289] = 3.56974825E+00;
    COFD[290] = -2.52221138E-01;
    COFD[291] = 1.10819767E-02;
    COFD[292] = -1.49160291E+01;
    COFD[293] = 3.56974825E+00;
    COFD[294] = -2.52221138E-01;
    COFD[295] = 1.10819767E-02;
    COFD[296] = -1.49160291E+01;
    COFD[297] = 3.56974825E+00;
    COFD[298] = -2.52221138E-01;
    COFD[299] = 1.10819767E-02;
    COFD[300] = -1.49181094E+01;
    COFD[301] = 3.56974825E+00;
    COFD[302] = -2.52221138E-01;
    COFD[303] = 1.10819767E-02;
    COFD[304] = -1.49197832E+01;
    COFD[305] = 3.56974825E+00;
    COFD[306] = -2.52221138E-01;
    COFD[307] = 1.10819767E-02;
    COFD[308] = -1.16906297E+01;
    COFD[309] = 2.47469981E+00;
    COFD[310] = -1.10436257E-01;
    COFD[311] = 4.95273813E-03;
    COFD[312] = -1.57972369E+01;
    COFD[313] = 4.22225052E+00;
    COFD[314] = -3.35156428E-01;
    COFD[315] = 1.46104855E-02;
    COFD[316] = -1.25098960E+01;
    COFD[317] = 2.77873601E+00;
    COFD[318] = -1.50637360E-01;
    COFD[319] = 6.72684281E-03;
    COFD[320] = -1.73027557E+01;
    COFD[321] = 4.21416723E+00;
    COFD[322] = -3.34163932E-01;
    COFD[323] = 1.45697432E-02;
    COFD[324] = -1.73027557E+01;
    COFD[325] = 4.21416723E+00;
    COFD[326] = -3.34163932E-01;
    COFD[327] = 1.45697432E-02;
    COFD[328] = -1.73198034E+01;
    COFD[329] = 4.21416723E+00;
    COFD[330] = -3.34163932E-01;
    COFD[331] = 1.45697432E-02;
    COFD[332] = -1.50584249E+01;
    COFD[333] = 3.47945612E+00;
    COFD[334] = -2.40703722E-01;
    COFD[335] = 1.05907441E-02;
    COFD[336] = -1.72556729E+01;
    COFD[337] = 4.19029808E+00;
    COFD[338] = -3.31177076E-01;
    COFD[339] = 1.44446234E-02;
    COFD[340] = -1.50724636E+01;
    COFD[341] = 3.47945612E+00;
    COFD[342] = -2.40703722E-01;
    COFD[343] = 1.05907441E-02;
    COFD[344] = -2.12639214E+01;
    COFD[345] = 5.61184117E+00;
    COFD[346] = -4.90532156E-01;
    COFD[347] = 2.03507922E-02;
    COFD[348] = -1.95548230E+01;
    COFD[349] = 4.97133070E+00;
    COFD[350] = -4.25604177E-01;
    COFD[351] = 1.82582594E-02;
    COFD[352] = -1.95613899E+01;
    COFD[353] = 4.97133070E+00;
    COFD[354] = -4.25604177E-01;
    COFD[355] = 1.82582594E-02;
    COFD[356] = -1.59634533E+01;
    COFD[357] = 3.67388294E+00;
    COFD[358] = -2.64990709E-01;
    COFD[359] = 1.16042706E-02;
    COFD[360] = -1.90996795E+01;
    COFD[361] = 4.82869066E+00;
    COFD[362] = -4.08564514E-01;
    COFD[363] = 1.75784675E-02;
    COFD[364] = -2.14160703E+01;
    COFD[365] = 5.56531152E+00;
    COFD[366] = -4.88789821E-01;
    COFD[367] = 2.04437116E-02;
    COFD[368] = -1.94530250E+01;
    COFD[369] = 4.87180830E+00;
    COFD[370] = -4.13582958E-01;
    COFD[371] = 1.77726094E-02;
    COFD[372] = -2.14215700E+01;
    COFD[373] = 5.56531152E+00;
    COFD[374] = -4.88789821E-01;
    COFD[375] = 2.04437116E-02;
    COFD[376] = -1.94585111E+01;
    COFD[377] = 4.87180830E+00;
    COFD[378] = -4.13582958E-01;
    COFD[379] = 1.77726094E-02;
    COFD[380] = -2.09376196E+01;
    COFD[381] = 5.40870099E+00;
    COFD[382] = -4.73017610E-01;
    COFD[383] = 1.99399066E-02;
    COFD[384] = -2.09376196E+01;
    COFD[385] = 5.40870099E+00;
    COFD[386] = -4.73017610E-01;
    COFD[387] = 1.99399066E-02;
    COFD[388] = -1.63254691E+01;
    COFD[389] = 3.82388595E+00;
    COFD[390] = -2.84480724E-01;
    COFD[391] = 1.24506311E-02;
    COFD[392] = -1.63301444E+01;
    COFD[393] = 3.82388595E+00;
    COFD[394] = -2.84480724E-01;
    COFD[395] = 1.24506311E-02;
    COFD[396] = -1.63345829E+01;
    COFD[397] = 3.82388595E+00;
    COFD[398] = -2.84480724E-01;
    COFD[399] = 1.24506311E-02;
    COFD[400] = -1.93015555E+01;
    COFD[401] = 4.85015581E+00;
    COFD[402] = -4.10945109E-01;
    COFD[403] = 1.76651398E-02;
    COFD[404] = -2.11406662E+01;
    COFD[405] = 5.42846112E+00;
    COFD[406] = -4.74321870E-01;
    COFD[407] = 1.99459749E-02;
    COFD[408] = -2.05539642E+01;
    COFD[409] = 5.20087227E+00;
    COFD[410] = -4.51444972E-01;
    COFD[411] = 1.92201898E-02;
    COFD[412] = -2.15351292E+01;
    COFD[413] = 5.51982454E+00;
    COFD[414] = -4.84452039E-01;
    COFD[415] = 2.03175522E-02;
    COFD[416] = -2.05565679E+01;
    COFD[417] = 5.20087227E+00;
    COFD[418] = -4.51444972E-01;
    COFD[419] = 1.92201898E-02;
    COFD[420] = -2.14048982E+01;
    COFD[421] = 5.54007827E+00;
    COFD[422] = -4.86434511E-01;
    COFD[423] = 2.03779006E-02;
    COFD[424] = -2.12255802E+01;
    COFD[425] = 5.39710574E+00;
    COFD[426] = -4.72221942E-01;
    COFD[427] = 1.99338345E-02;
    COFD[428] = -2.12271939E+01;
    COFD[429] = 5.39710574E+00;
    COFD[430] = -4.72221942E-01;
    COFD[431] = 1.99338345E-02;
    COFD[432] = -2.16737420E+01;
    COFD[433] = 5.52035712E+00;
    COFD[434] = -4.84503119E-01;
    COFD[435] = 2.03190484E-02;
    COFD[436] = -2.05875936E+01;
    COFD[437] = 5.20087227E+00;
    COFD[438] = -4.51444972E-01;
    COFD[439] = 1.92201898E-02;
    COFD[440] = -2.06037892E+01;
    COFD[441] = 5.20087227E+00;
    COFD[442] = -4.51444972E-01;
    COFD[443] = 1.92201898E-02;
    COFD[444] = -2.06037892E+01;
    COFD[445] = 5.20087227E+00;
    COFD[446] = -4.51444972E-01;
    COFD[447] = 1.92201898E-02;
    COFD[448] = -2.06058533E+01;
    COFD[449] = 5.20087227E+00;
    COFD[450] = -4.51444972E-01;
    COFD[451] = 1.92201898E-02;
    COFD[452] = -2.06058533E+01;
    COFD[453] = 5.20087227E+00;
    COFD[454] = -4.51444972E-01;
    COFD[455] = 1.92201898E-02;
    COFD[456] = -2.06185530E+01;
    COFD[457] = 5.20087227E+00;
    COFD[458] = -4.51444972E-01;
    COFD[459] = 1.92201898E-02;
    COFD[460] = -2.06289714E+01;
    COFD[461] = 5.20087227E+00;
    COFD[462] = -4.51444972E-01;
    COFD[463] = 1.92201898E-02;
    COFD[464] = -1.59404882E+01;
    COFD[465] = 3.66853818E+00;
    COFD[466] = -2.64346221E-01;
    COFD[467] = 1.15784613E-02;
    COFD[468] = -1.57972369E+01;
    COFD[469] = 4.22225052E+00;
    COFD[470] = -3.35156428E-01;
    COFD[471] = 1.46104855E-02;
    COFD[472] = -1.25098960E+01;
    COFD[473] = 2.77873601E+00;
    COFD[474] = -1.50637360E-01;
    COFD[475] = 6.72684281E-03;
    COFD[476] = -1.73027557E+01;
    COFD[477] = 4.21416723E+00;
    COFD[478] = -3.34163932E-01;
    COFD[479] = 1.45697432E-02;
    COFD[480] = -1.73027557E+01;
    COFD[481] = 4.21416723E+00;
    COFD[482] = -3.34163932E-01;
    COFD[483] = 1.45697432E-02;
    COFD[484] = -1.73198034E+01;
    COFD[485] = 4.21416723E+00;
    COFD[486] = -3.34163932E-01;
    COFD[487] = 1.45697432E-02;
    COFD[488] = -1.50584249E+01;
    COFD[489] = 3.47945612E+00;
    COFD[490] = -2.40703722E-01;
    COFD[491] = 1.05907441E-02;
    COFD[492] = -1.72556729E+01;
    COFD[493] = 4.19029808E+00;
    COFD[494] = -3.31177076E-01;
    COFD[495] = 1.44446234E-02;
    COFD[496] = -1.50724636E+01;
    COFD[497] = 3.47945612E+00;
    COFD[498] = -2.40703722E-01;
    COFD[499] = 1.05907441E-02;
    COFD[500] = -2.12639214E+01;
    COFD[501] = 5.61184117E+00;
    COFD[502] = -4.90532156E-01;
    COFD[503] = 2.03507922E-02;
    COFD[504] = -1.95548230E+01;
    COFD[505] = 4.97133070E+00;
    COFD[506] = -4.25604177E-01;
    COFD[507] = 1.82582594E-02;
    COFD[508] = -1.95613899E+01;
    COFD[509] = 4.97133070E+00;
    COFD[510] = -4.25604177E-01;
    COFD[511] = 1.82582594E-02;
    COFD[512] = -1.59634533E+01;
    COFD[513] = 3.67388294E+00;
    COFD[514] = -2.64990709E-01;
    COFD[515] = 1.16042706E-02;
    COFD[516] = -1.90996795E+01;
    COFD[517] = 4.82869066E+00;
    COFD[518] = -4.08564514E-01;
    COFD[519] = 1.75784675E-02;
    COFD[520] = -2.14160703E+01;
    COFD[521] = 5.56531152E+00;
    COFD[522] = -4.88789821E-01;
    COFD[523] = 2.04437116E-02;
    COFD[524] = -1.94530250E+01;
    COFD[525] = 4.87180830E+00;
    COFD[526] = -4.13582958E-01;
    COFD[527] = 1.77726094E-02;
    COFD[528] = -2.14215700E+01;
    COFD[529] = 5.56531152E+00;
    COFD[530] = -4.88789821E-01;
    COFD[531] = 2.04437116E-02;
    COFD[532] = -1.94585111E+01;
    COFD[533] = 4.87180830E+00;
    COFD[534] = -4.13582958E-01;
    COFD[535] = 1.77726094E-02;
    COFD[536] = -2.09376196E+01;
    COFD[537] = 5.40870099E+00;
    COFD[538] = -4.73017610E-01;
    COFD[539] = 1.99399066E-02;
    COFD[540] = -2.09376196E+01;
    COFD[541] = 5.40870099E+00;
    COFD[542] = -4.73017610E-01;
    COFD[543] = 1.99399066E-02;
    COFD[544] = -1.63254691E+01;
    COFD[545] = 3.82388595E+00;
    COFD[546] = -2.84480724E-01;
    COFD[547] = 1.24506311E-02;
    COFD[548] = -1.63301444E+01;
    COFD[549] = 3.82388595E+00;
    COFD[550] = -2.84480724E-01;
    COFD[551] = 1.24506311E-02;
    COFD[552] = -1.63345829E+01;
    COFD[553] = 3.82388595E+00;
    COFD[554] = -2.84480724E-01;
    COFD[555] = 1.24506311E-02;
    COFD[556] = -1.93015555E+01;
    COFD[557] = 4.85015581E+00;
    COFD[558] = -4.10945109E-01;
    COFD[559] = 1.76651398E-02;
    COFD[560] = -2.11406662E+01;
    COFD[561] = 5.42846112E+00;
    COFD[562] = -4.74321870E-01;
    COFD[563] = 1.99459749E-02;
    COFD[564] = -2.05539642E+01;
    COFD[565] = 5.20087227E+00;
    COFD[566] = -4.51444972E-01;
    COFD[567] = 1.92201898E-02;
    COFD[568] = -2.15351292E+01;
    COFD[569] = 5.51982454E+00;
    COFD[570] = -4.84452039E-01;
    COFD[571] = 2.03175522E-02;
    COFD[572] = -2.05565679E+01;
    COFD[573] = 5.20087227E+00;
    COFD[574] = -4.51444972E-01;
    COFD[575] = 1.92201898E-02;
    COFD[576] = -2.14048982E+01;
    COFD[577] = 5.54007827E+00;
    COFD[578] = -4.86434511E-01;
    COFD[579] = 2.03779006E-02;
    COFD[580] = -2.12255802E+01;
    COFD[581] = 5.39710574E+00;
    COFD[582] = -4.72221942E-01;
    COFD[583] = 1.99338345E-02;
    COFD[584] = -2.12271939E+01;
    COFD[585] = 5.39710574E+00;
    COFD[586] = -4.72221942E-01;
    COFD[587] = 1.99338345E-02;
    COFD[588] = -2.16737420E+01;
    COFD[589] = 5.52035712E+00;
    COFD[590] = -4.84503119E-01;
    COFD[591] = 2.03190484E-02;
    COFD[592] = -2.05875936E+01;
    COFD[593] = 5.20087227E+00;
    COFD[594] = -4.51444972E-01;
    COFD[595] = 1.92201898E-02;
    COFD[596] = -2.06037892E+01;
    COFD[597] = 5.20087227E+00;
    COFD[598] = -4.51444972E-01;
    COFD[599] = 1.92201898E-02;
    COFD[600] = -2.06037892E+01;
    COFD[601] = 5.20087227E+00;
    COFD[602] = -4.51444972E-01;
    COFD[603] = 1.92201898E-02;
    COFD[604] = -2.06058533E+01;
    COFD[605] = 5.20087227E+00;
    COFD[606] = -4.51444972E-01;
    COFD[607] = 1.92201898E-02;
    COFD[608] = -2.06058533E+01;
    COFD[609] = 5.20087227E+00;
    COFD[610] = -4.51444972E-01;
    COFD[611] = 1.92201898E-02;
    COFD[612] = -2.06185530E+01;
    COFD[613] = 5.20087227E+00;
    COFD[614] = -4.51444972E-01;
    COFD[615] = 1.92201898E-02;
    COFD[616] = -2.06289714E+01;
    COFD[617] = 5.20087227E+00;
    COFD[618] = -4.51444972E-01;
    COFD[619] = 1.92201898E-02;
    COFD[620] = -1.59404882E+01;
    COFD[621] = 3.66853818E+00;
    COFD[622] = -2.64346221E-01;
    COFD[623] = 1.15784613E-02;
    COFD[624] = -1.57994893E+01;
    COFD[625] = 4.22225052E+00;
    COFD[626] = -3.35156428E-01;
    COFD[627] = 1.46104855E-02;
    COFD[628] = -1.25141260E+01;
    COFD[629] = 2.77873601E+00;
    COFD[630] = -1.50637360E-01;
    COFD[631] = 6.72684281E-03;
    COFD[632] = -1.73198034E+01;
    COFD[633] = 4.21416723E+00;
    COFD[634] = -3.34163932E-01;
    COFD[635] = 1.45697432E-02;
    COFD[636] = -1.73198034E+01;
    COFD[637] = 4.21416723E+00;
    COFD[638] = -3.34163932E-01;
    COFD[639] = 1.45697432E-02;
    COFD[640] = -1.73374529E+01;
    COFD[641] = 4.21416723E+00;
    COFD[642] = -3.34163932E-01;
    COFD[643] = 1.45697432E-02;
    COFD[644] = -1.50766130E+01;
    COFD[645] = 3.47945612E+00;
    COFD[646] = -2.40703722E-01;
    COFD[647] = 1.05907441E-02;
    COFD[648] = -1.72738845E+01;
    COFD[649] = 4.19029808E+00;
    COFD[650] = -3.31177076E-01;
    COFD[651] = 1.44446234E-02;
    COFD[652] = -1.50911794E+01;
    COFD[653] = 3.47945612E+00;
    COFD[654] = -2.40703722E-01;
    COFD[655] = 1.05907441E-02;
    COFD[656] = -2.12831323E+01;
    COFD[657] = 5.61184117E+00;
    COFD[658] = -4.90532156E-01;
    COFD[659] = 2.03507922E-02;
    COFD[660] = -1.95770968E+01;
    COFD[661] = 4.97133070E+00;
    COFD[662] = -4.25604177E-01;
    COFD[663] = 1.82582594E-02;
    COFD[664] = -1.95839648E+01;
    COFD[665] = 4.97133070E+00;
    COFD[666] = -4.25604177E-01;
    COFD[667] = 1.82582594E-02;
    COFD[668] = -1.59863030E+01;
    COFD[669] = 3.67388294E+00;
    COFD[670] = -2.64990709E-01;
    COFD[671] = 1.16042706E-02;
    COFD[672] = -1.91225414E+01;
    COFD[673] = 4.82869066E+00;
    COFD[674] = -4.08564514E-01;
    COFD[675] = 1.75784675E-02;
    COFD[676] = -2.14391943E+01;
    COFD[677] = 5.56531152E+00;
    COFD[678] = -4.88789821E-01;
    COFD[679] = 2.04437116E-02;
    COFD[680] = -1.94761606E+01;
    COFD[681] = 4.87180830E+00;
    COFD[682] = -4.13582958E-01;
    COFD[683] = 1.77726094E-02;
    COFD[684] = -2.14449559E+01;
    COFD[685] = 5.56531152E+00;
    COFD[686] = -4.88789821E-01;
    COFD[687] = 2.04437116E-02;
    COFD[688] = -1.94819080E+01;
    COFD[689] = 4.87180830E+00;
    COFD[690] = -4.13582958E-01;
    COFD[691] = 1.77726094E-02;
    COFD[692] = -2.09612557E+01;
    COFD[693] = 5.40870099E+00;
    COFD[694] = -4.73017610E-01;
    COFD[695] = 1.99399066E-02;
    COFD[696] = -2.09612557E+01;
    COFD[697] = 5.40870099E+00;
    COFD[698] = -4.73017610E-01;
    COFD[699] = 1.99399066E-02;
    COFD[700] = -1.63493345E+01;
    COFD[701] = 3.82388595E+00;
    COFD[702] = -2.84480724E-01;
    COFD[703] = 1.24506311E-02;
    COFD[704] = -1.63542394E+01;
    COFD[705] = 3.82388595E+00;
    COFD[706] = -2.84480724E-01;
    COFD[707] = 1.24506311E-02;
    COFD[708] = -1.63588981E+01;
    COFD[709] = 3.82388595E+00;
    COFD[710] = -2.84480724E-01;
    COFD[711] = 1.24506311E-02;
    COFD[712] = -1.93276434E+01;
    COFD[713] = 4.85015581E+00;
    COFD[714] = -4.10945109E-01;
    COFD[715] = 1.76651398E-02;
    COFD[716] = -2.11667605E+01;
    COFD[717] = 5.42846112E+00;
    COFD[718] = -4.74321870E-01;
    COFD[719] = 1.99459749E-02;
    COFD[720] = -2.05802040E+01;
    COFD[721] = 5.20087227E+00;
    COFD[722] = -4.51444972E-01;
    COFD[723] = 1.92201898E-02;
    COFD[724] = -2.15615037E+01;
    COFD[725] = 5.51982454E+00;
    COFD[726] = -4.84452039E-01;
    COFD[727] = 2.03175522E-02;
    COFD[728] = -2.05829484E+01;
    COFD[729] = 5.20087227E+00;
    COFD[730] = -4.51444972E-01;
    COFD[731] = 1.92201898E-02;
    COFD[732] = -2.14314090E+01;
    COFD[733] = 5.54007827E+00;
    COFD[734] = -4.86434511E-01;
    COFD[735] = 2.03779006E-02;
    COFD[736] = -2.12534275E+01;
    COFD[737] = 5.39710574E+00;
    COFD[738] = -4.72221942E-01;
    COFD[739] = 1.99338345E-02;
    COFD[740] = -2.12551337E+01;
    COFD[741] = 5.39710574E+00;
    COFD[742] = -4.72221942E-01;
    COFD[743] = 1.99338345E-02;
    COFD[744] = -2.17017719E+01;
    COFD[745] = 5.52035712E+00;
    COFD[746] = -4.84503119E-01;
    COFD[747] = 2.03190484E-02;
    COFD[748] = -2.06157113E+01;
    COFD[749] = 5.20087227E+00;
    COFD[750] = -4.51444972E-01;
    COFD[751] = 1.92201898E-02;
    COFD[752] = -2.06328599E+01;
    COFD[753] = 5.20087227E+00;
    COFD[754] = -4.51444972E-01;
    COFD[755] = 1.92201898E-02;
    COFD[756] = -2.06328599E+01;
    COFD[757] = 5.20087227E+00;
    COFD[758] = -4.51444972E-01;
    COFD[759] = 1.92201898E-02;
    COFD[760] = -2.06350479E+01;
    COFD[761] = 5.20087227E+00;
    COFD[762] = -4.51444972E-01;
    COFD[763] = 1.92201898E-02;
    COFD[764] = -2.06350479E+01;
    COFD[765] = 5.20087227E+00;
    COFD[766] = -4.51444972E-01;
    COFD[767] = 1.92201898E-02;
    COFD[768] = -2.06485216E+01;
    COFD[769] = 5.20087227E+00;
    COFD[770] = -4.51444972E-01;
    COFD[771] = 1.92201898E-02;
    COFD[772] = -2.06595907E+01;
    COFD[773] = 5.20087227E+00;
    COFD[774] = -4.51444972E-01;
    COFD[775] = 1.92201898E-02;
    COFD[776] = -1.59633387E+01;
    COFD[777] = 3.66853818E+00;
    COFD[778] = -2.64346221E-01;
    COFD[779] = 1.15784613E-02;
    COFD[780] = -1.34230272E+01;
    COFD[781] = 3.48624238E+00;
    COFD[782] = -2.41554467E-01;
    COFD[783] = 1.06263545E-02;
    COFD[784] = -1.09595712E+01;
    COFD[785] = 2.30836460E+00;
    COFD[786] = -8.76339315E-02;
    COFD[787] = 3.90878445E-03;
    COFD[788] = -1.50584249E+01;
    COFD[789] = 3.47945612E+00;
    COFD[790] = -2.40703722E-01;
    COFD[791] = 1.05907441E-02;
    COFD[792] = -1.50584249E+01;
    COFD[793] = 3.47945612E+00;
    COFD[794] = -2.40703722E-01;
    COFD[795] = 1.05907441E-02;
    COFD[796] = -1.50766130E+01;
    COFD[797] = 3.47945612E+00;
    COFD[798] = -2.40703722E-01;
    COFD[799] = 1.05907441E-02;
    COFD[800] = -1.32093628E+01;
    COFD[801] = 2.90778936E+00;
    COFD[802] = -1.67388544E-01;
    COFD[803] = 7.45220609E-03;
    COFD[804] = -1.50270339E+01;
    COFD[805] = 3.46140064E+00;
    COFD[806] = -2.38440092E-01;
    COFD[807] = 1.04960087E-02;
    COFD[808] = -1.32244035E+01;
    COFD[809] = 2.90778936E+00;
    COFD[810] = -1.67388544E-01;
    COFD[811] = 7.45220609E-03;
    COFD[812] = -1.94093572E+01;
    COFD[813] = 5.16013126E+00;
    COFD[814] = -4.46824543E-01;
    COFD[815] = 1.90464887E-02;
    COFD[816] = -1.72286007E+01;
    COFD[817] = 4.24084025E+00;
    COFD[818] = -3.37428619E-01;
    COFD[819] = 1.47032793E-02;
    COFD[820] = -1.72357436E+01;
    COFD[821] = 4.24084025E+00;
    COFD[822] = -3.37428619E-01;
    COFD[823] = 1.47032793E-02;
    COFD[824] = -1.40999008E+01;
    COFD[825] = 3.08120012E+00;
    COFD[826] = -1.89629903E-01;
    COFD[827] = 8.40361952E-03;
    COFD[828] = -1.68343393E+01;
    COFD[829] = 4.11954900E+00;
    COFD[830] = -3.22470391E-01;
    COFD[831] = 1.40859564E-02;
    COFD[832] = -1.94313116E+01;
    COFD[833] = 5.02567894E+00;
    COFD[834] = -4.32045169E-01;
    COFD[835] = 1.85132214E-02;
    COFD[836] = -1.72053106E+01;
    COFD[837] = 4.15807461E+00;
    COFD[838] = -3.27178539E-01;
    COFD[839] = 1.42784349E-02;
    COFD[840] = -1.94373127E+01;
    COFD[841] = 5.02567894E+00;
    COFD[842] = -4.32045169E-01;
    COFD[843] = 1.85132214E-02;
    COFD[844] = -1.72112971E+01;
    COFD[845] = 4.15807461E+00;
    COFD[846] = -3.27178539E-01;
    COFD[847] = 1.42784349E-02;
    COFD[848] = -1.88179418E+01;
    COFD[849] = 4.79683898E+00;
    COFD[850] = -4.04829719E-01;
    COFD[851] = 1.74325475E-02;
    COFD[852] = -1.88179418E+01;
    COFD[853] = 4.79683898E+00;
    COFD[854] = -4.04829719E-01;
    COFD[855] = 1.74325475E-02;
    COFD[856] = -1.43139231E+01;
    COFD[857] = 3.17651319E+00;
    COFD[858] = -2.02028974E-01;
    COFD[859] = 8.94232502E-03;
    COFD[860] = -1.43190389E+01;
    COFD[861] = 3.17651319E+00;
    COFD[862] = -2.02028974E-01;
    COFD[863] = 8.94232502E-03;
    COFD[864] = -1.43238998E+01;
    COFD[865] = 3.17651319E+00;
    COFD[866] = -2.02028974E-01;
    COFD[867] = 8.94232502E-03;
    COFD[868] = -1.70534856E+01;
    COFD[869] = 4.14240922E+00;
    COFD[870] = -3.25239774E-01;
    COFD[871] = 1.41980687E-02;
    COFD[872] = -1.90946745E+01;
    COFD[873] = 4.84384483E+00;
    COFD[874] = -4.10265575E-01;
    COFD[875] = 1.76414287E-02;
    COFD[876] = -1.83038349E+01;
    COFD[877] = 4.49977587E+00;
    COFD[878] = -3.68989022E-01;
    COFD[879] = 1.59879891E-02;
    COFD[880] = -1.95342215E+01;
    COFD[881] = 4.95249173E+00;
    COFD[882] = -4.23376552E-01;
    COFD[883] = 1.81703714E-02;
    COFD[884] = -1.83067096E+01;
    COFD[885] = 4.49977587E+00;
    COFD[886] = -3.68989022E-01;
    COFD[887] = 1.59879891E-02;
    COFD[888] = -1.93925667E+01;
    COFD[889] = 4.98286777E+00;
    COFD[890] = -4.26970814E-01;
    COFD[891] = 1.83122917E-02;
    COFD[892] = -1.91187832E+01;
    COFD[893] = 4.76920246E+00;
    COFD[894] = -4.01609656E-01;
    COFD[895] = 1.73077246E-02;
    COFD[896] = -1.91205757E+01;
    COFD[897] = 4.76920246E+00;
    COFD[898] = -4.01609656E-01;
    COFD[899] = 1.73077246E-02;
    COFD[900] = -1.96918321E+01;
    COFD[901] = 4.95331445E+00;
    COFD[902] = -4.23474055E-01;
    COFD[903] = 1.81742301E-02;
    COFD[904] = -1.83410870E+01;
    COFD[905] = 4.49977587E+00;
    COFD[906] = -3.68989022E-01;
    COFD[907] = 1.59879891E-02;
    COFD[908] = -1.83591261E+01;
    COFD[909] = 4.49977587E+00;
    COFD[910] = -3.68989022E-01;
    COFD[911] = 1.59879891E-02;
    COFD[912] = -1.83591261E+01;
    COFD[913] = 4.49977587E+00;
    COFD[914] = -3.68989022E-01;
    COFD[915] = 1.59879891E-02;
    COFD[916] = -1.83614301E+01;
    COFD[917] = 4.49977587E+00;
    COFD[918] = -3.68989022E-01;
    COFD[919] = 1.59879891E-02;
    COFD[920] = -1.83614301E+01;
    COFD[921] = 4.49977587E+00;
    COFD[922] = -3.68989022E-01;
    COFD[923] = 1.59879891E-02;
    COFD[924] = -1.83756296E+01;
    COFD[925] = 4.49977587E+00;
    COFD[926] = -3.68989022E-01;
    COFD[927] = 1.59879891E-02;
    COFD[928] = -1.83873107E+01;
    COFD[929] = 4.49977587E+00;
    COFD[930] = -3.68989022E-01;
    COFD[931] = 1.59879891E-02;
    COFD[932] = -1.40756935E+01;
    COFD[933] = 3.07549274E+00;
    COFD[934] = -1.88889344E-01;
    COFD[935] = 8.37152866E-03;
    COFD[936] = -1.57199037E+01;
    COFD[937] = 4.19936335E+00;
    COFD[938] = -3.32311009E-01;
    COFD[939] = 1.44921003E-02;
    COFD[940] = -1.24693568E+01;
    COFD[941] = 2.76686648E+00;
    COFD[942] = -1.49120141E-01;
    COFD[943] = 6.66220432E-03;
    COFD[944] = -1.72556729E+01;
    COFD[945] = 4.19029808E+00;
    COFD[946] = -3.31177076E-01;
    COFD[947] = 1.44446234E-02;
    COFD[948] = -1.72556729E+01;
    COFD[949] = 4.19029808E+00;
    COFD[950] = -3.31177076E-01;
    COFD[951] = 1.44446234E-02;
    COFD[952] = -1.72738845E+01;
    COFD[953] = 4.19029808E+00;
    COFD[954] = -3.31177076E-01;
    COFD[955] = 1.44446234E-02;
    COFD[956] = -1.50270339E+01;
    COFD[957] = 3.46140064E+00;
    COFD[958] = -2.38440092E-01;
    COFD[959] = 1.04960087E-02;
    COFD[960] = -1.72167708E+01;
    COFD[961] = 4.16886779E+00;
    COFD[962] = -3.28518156E-01;
    COFD[963] = 1.43341626E-02;
    COFD[964] = -1.50420953E+01;
    COFD[965] = 3.46140064E+00;
    COFD[966] = -2.38440092E-01;
    COFD[967] = 1.04960087E-02;
    COFD[968] = -2.14087397E+01;
    COFD[969] = 5.57282008E+00;
    COFD[970] = -4.76690890E-01;
    COFD[971] = 1.94000719E-02;
    COFD[972] = -1.95154079E+01;
    COFD[973] = 4.94787350E+00;
    COFD[974] = -4.22829292E-01;
    COFD[975] = 1.81487163E-02;
    COFD[976] = -1.95225629E+01;
    COFD[977] = 4.94787350E+00;
    COFD[978] = -4.22829292E-01;
    COFD[979] = 1.81487163E-02;
    COFD[980] = -1.59525102E+01;
    COFD[981] = 3.66023858E+00;
    COFD[982] = -2.63401043E-01;
    COFD[983] = 1.15432000E-02;
    COFD[984] = -1.90692595E+01;
    COFD[985] = 4.80830699E+00;
    COFD[986] = -4.06171933E-01;
    COFD[987] = 1.74848791E-02;
    COFD[988] = -2.14022336E+01;
    COFD[989] = 5.55346617E+00;
    COFD[990] = -4.87783156E-01;
    COFD[991] = 2.04210886E-02;
    COFD[992] = -1.94126575E+01;
    COFD[993] = 4.84669430E+00;
    COFD[994] = -4.10571455E-01;
    COFD[995] = 1.76520543E-02;
    COFD[996] = -2.14082453E+01;
    COFD[997] = 5.55346617E+00;
    COFD[998] = -4.87783156E-01;
    COFD[999] = 2.04210886E-02;
    COFD[1000] = -1.94186547E+01;
    COFD[1001] = 4.84669430E+00;
    COFD[1002] = -4.10571455E-01;
    COFD[1003] = 1.76520543E-02;
    COFD[1004] = -2.11381508E+01;
    COFD[1005] = 5.45574440E+00;
    COFD[1006] = -4.77436155E-01;
    COFD[1007] = 2.00644596E-02;
    COFD[1008] = -2.11381508E+01;
    COFD[1009] = 5.45574440E+00;
    COFD[1010] = -4.77436155E-01;
    COFD[1011] = 2.00644596E-02;
    COFD[1012] = -1.62724462E+01;
    COFD[1013] = 3.79163564E+00;
    COFD[1014] = -2.80257365E-01;
    COFD[1015] = 1.22656902E-02;
    COFD[1016] = -1.62775714E+01;
    COFD[1017] = 3.79163564E+00;
    COFD[1018] = -2.80257365E-01;
    COFD[1019] = 1.22656902E-02;
    COFD[1020] = -1.62824412E+01;
    COFD[1021] = 3.79163564E+00;
    COFD[1022] = -2.80257365E-01;
    COFD[1023] = 1.22656902E-02;
    COFD[1024] = -1.92867554E+01;
    COFD[1025] = 4.83375900E+00;
    COFD[1026] = -4.09146560E-01;
    COFD[1027] = 1.76006599E-02;
    COFD[1028] = -2.11372811E+01;
    COFD[1029] = 5.41773516E+00;
    COFD[1030] = -4.73414338E-01;
    COFD[1031] = 1.99258685E-02;
    COFD[1032] = -2.05236383E+01;
    COFD[1033] = 5.17771473E+00;
    COFD[1034] = -4.48750201E-01;
    COFD[1035] = 1.91155567E-02;
    COFD[1036] = -2.15095161E+01;
    COFD[1037] = 5.49964831E+00;
    COFD[1038] = -4.82275380E-01;
    COFD[1039] = 2.02405072E-02;
    COFD[1040] = -2.05265188E+01;
    COFD[1041] = 5.17771473E+00;
    COFD[1042] = -4.48750201E-01;
    COFD[1043] = 1.91155567E-02;
    COFD[1044] = -2.13881945E+01;
    COFD[1045] = 5.52422470E+00;
    COFD[1046] = -4.84872944E-01;
    COFD[1047] = 2.03298213E-02;
    COFD[1048] = -2.12314042E+01;
    COFD[1049] = 5.38826297E+00;
    COFD[1050] = -4.71579030E-01;
    COFD[1051] = 1.99262449E-02;
    COFD[1052] = -2.12332005E+01;
    COFD[1053] = 5.38826297E+00;
    COFD[1054] = -4.71579030E-01;
    COFD[1055] = 1.99262449E-02;
    COFD[1056] = -2.16525572E+01;
    COFD[1057] = 5.50038516E+00;
    COFD[1058] = -4.82355440E-01;
    COFD[1059] = 2.02433624E-02;
    COFD[1060] = -2.05609682E+01;
    COFD[1061] = 5.17771473E+00;
    COFD[1062] = -4.48750201E-01;
    COFD[1063] = 1.91155567E-02;
    COFD[1064] = -2.05790471E+01;
    COFD[1065] = 5.17771473E+00;
    COFD[1066] = -4.48750201E-01;
    COFD[1067] = 1.91155567E-02;
    COFD[1068] = -2.05790471E+01;
    COFD[1069] = 5.17771473E+00;
    COFD[1070] = -4.48750201E-01;
    COFD[1071] = 1.91155567E-02;
    COFD[1072] = -2.05813563E+01;
    COFD[1073] = 5.17771473E+00;
    COFD[1074] = -4.48750201E-01;
    COFD[1075] = 1.91155567E-02;
    COFD[1076] = -2.05813563E+01;
    COFD[1077] = 5.17771473E+00;
    COFD[1078] = -4.48750201E-01;
    COFD[1079] = 1.91155567E-02;
    COFD[1080] = -2.05955883E+01;
    COFD[1081] = 5.17771473E+00;
    COFD[1082] = -4.48750201E-01;
    COFD[1083] = 1.91155567E-02;
    COFD[1084] = -2.06072968E+01;
    COFD[1085] = 5.17771473E+00;
    COFD[1086] = -4.48750201E-01;
    COFD[1087] = 1.91155567E-02;
    COFD[1088] = -1.59327297E+01;
    COFD[1089] = 3.65620899E+00;
    COFD[1090] = -2.62933804E-01;
    COFD[1091] = 1.15253223E-02;
    COFD[1092] = -1.34247866E+01;
    COFD[1093] = 3.48624238E+00;
    COFD[1094] = -2.41554467E-01;
    COFD[1095] = 1.06263545E-02;
    COFD[1096] = -1.09628982E+01;
    COFD[1097] = 2.30836460E+00;
    COFD[1098] = -8.76339315E-02;
    COFD[1099] = 3.90878445E-03;
    COFD[1100] = -1.50724636E+01;
    COFD[1101] = 3.47945612E+00;
    COFD[1102] = -2.40703722E-01;
    COFD[1103] = 1.05907441E-02;
    COFD[1104] = -1.50724636E+01;
    COFD[1105] = 3.47945612E+00;
    COFD[1106] = -2.40703722E-01;
    COFD[1107] = 1.05907441E-02;
    COFD[1108] = -1.50911794E+01;
    COFD[1109] = 3.47945612E+00;
    COFD[1110] = -2.40703722E-01;
    COFD[1111] = 1.05907441E-02;
    COFD[1112] = -1.32244035E+01;
    COFD[1113] = 2.90778936E+00;
    COFD[1114] = -1.67388544E-01;
    COFD[1115] = 7.45220609E-03;
    COFD[1116] = -1.50420953E+01;
    COFD[1117] = 3.46140064E+00;
    COFD[1118] = -2.38440092E-01;
    COFD[1119] = 1.04960087E-02;
    COFD[1120] = -1.32399106E+01;
    COFD[1121] = 2.90778936E+00;
    COFD[1122] = -1.67388544E-01;
    COFD[1123] = 7.45220609E-03;
    COFD[1124] = -1.94253036E+01;
    COFD[1125] = 5.16013126E+00;
    COFD[1126] = -4.46824543E-01;
    COFD[1127] = 1.90464887E-02;
    COFD[1128] = -1.72473011E+01;
    COFD[1129] = 4.24084025E+00;
    COFD[1130] = -3.37428619E-01;
    COFD[1131] = 1.47032793E-02;
    COFD[1132] = -1.72547182E+01;
    COFD[1133] = 4.24084025E+00;
    COFD[1134] = -3.37428619E-01;
    COFD[1135] = 1.47032793E-02;
    COFD[1136] = -1.41191261E+01;
    COFD[1137] = 3.08120012E+00;
    COFD[1138] = -1.89629903E-01;
    COFD[1139] = 8.40361952E-03;
    COFD[1140] = -1.68535757E+01;
    COFD[1141] = 4.11954900E+00;
    COFD[1142] = -3.22470391E-01;
    COFD[1143] = 1.40859564E-02;
    COFD[1144] = -1.94507876E+01;
    COFD[1145] = 5.02567894E+00;
    COFD[1146] = -4.32045169E-01;
    COFD[1147] = 1.85132214E-02;
    COFD[1148] = -1.72247972E+01;
    COFD[1149] = 4.15807461E+00;
    COFD[1150] = -3.27178539E-01;
    COFD[1151] = 1.42784349E-02;
    COFD[1152] = -1.94570287E+01;
    COFD[1153] = 5.02567894E+00;
    COFD[1154] = -4.32045169E-01;
    COFD[1155] = 1.85132214E-02;
    COFD[1156] = -1.72310232E+01;
    COFD[1157] = 4.15807461E+00;
    COFD[1158] = -3.27178539E-01;
    COFD[1159] = 1.42784349E-02;
    COFD[1160] = -1.88378874E+01;
    COFD[1161] = 4.79683898E+00;
    COFD[1162] = -4.04829719E-01;
    COFD[1163] = 1.74325475E-02;
    COFD[1164] = -1.88378874E+01;
    COFD[1165] = 4.79683898E+00;
    COFD[1166] = -4.04829719E-01;
    COFD[1167] = 1.74325475E-02;
    COFD[1168] = -1.43340796E+01;
    COFD[1169] = 3.17651319E+00;
    COFD[1170] = -2.02028974E-01;
    COFD[1171] = 8.94232502E-03;
    COFD[1172] = -1.43394069E+01;
    COFD[1173] = 3.17651319E+00;
    COFD[1174] = -2.02028974E-01;
    COFD[1175] = 8.94232502E-03;
    COFD[1176] = -1.43444709E+01;
    COFD[1177] = 3.17651319E+00;
    COFD[1178] = -2.02028974E-01;
    COFD[1179] = 8.94232502E-03;
    COFD[1180] = -1.70757047E+01;
    COFD[1181] = 4.14240922E+00;
    COFD[1182] = -3.25239774E-01;
    COFD[1183] = 1.41980687E-02;
    COFD[1184] = -1.91168996E+01;
    COFD[1185] = 4.84384483E+00;
    COFD[1186] = -4.10265575E-01;
    COFD[1187] = 1.76414287E-02;
    COFD[1188] = -1.83261963E+01;
    COFD[1189] = 4.49977587E+00;
    COFD[1190] = -3.68989022E-01;
    COFD[1191] = 1.59879891E-02;
    COFD[1192] = -1.95567091E+01;
    COFD[1193] = 4.95249173E+00;
    COFD[1194] = -4.23376552E-01;
    COFD[1195] = 1.81703714E-02;
    COFD[1196] = -1.83292028E+01;
    COFD[1197] = 4.49977587E+00;
    COFD[1198] = -3.68989022E-01;
    COFD[1199] = 1.59879891E-02;
    COFD[1200] = -1.94151822E+01;
    COFD[1201] = 4.98286777E+00;
    COFD[1202] = -4.26970814E-01;
    COFD[1203] = 1.83122917E-02;
    COFD[1204] = -1.91426599E+01;
    COFD[1205] = 4.76920246E+00;
    COFD[1206] = -4.01609656E-01;
    COFD[1207] = 1.73077246E-02;
    COFD[1208] = -1.91445402E+01;
    COFD[1209] = 4.76920246E+00;
    COFD[1210] = -4.01609656E-01;
    COFD[1211] = 1.73077246E-02;
    COFD[1212] = -1.97158821E+01;
    COFD[1213] = 4.95331445E+00;
    COFD[1214] = -4.23474055E-01;
    COFD[1215] = 1.81742301E-02;
    COFD[1216] = -1.83652204E+01;
    COFD[1217] = 4.49977587E+00;
    COFD[1218] = -3.68989022E-01;
    COFD[1219] = 1.59879891E-02;
    COFD[1220] = -1.83841687E+01;
    COFD[1221] = 4.49977587E+00;
    COFD[1222] = -3.68989022E-01;
    COFD[1223] = 1.59879891E-02;
    COFD[1224] = -1.83841687E+01;
    COFD[1225] = 4.49977587E+00;
    COFD[1226] = -3.68989022E-01;
    COFD[1227] = 1.59879891E-02;
    COFD[1228] = -1.83865912E+01;
    COFD[1229] = 4.49977587E+00;
    COFD[1230] = -3.68989022E-01;
    COFD[1231] = 1.59879891E-02;
    COFD[1232] = -1.83865912E+01;
    COFD[1233] = 4.49977587E+00;
    COFD[1234] = -3.68989022E-01;
    COFD[1235] = 1.59879891E-02;
    COFD[1236] = -1.84015347E+01;
    COFD[1237] = 4.49977587E+00;
    COFD[1238] = -3.68989022E-01;
    COFD[1239] = 1.59879891E-02;
    COFD[1240] = -1.84138446E+01;
    COFD[1241] = 4.49977587E+00;
    COFD[1242] = -3.68989022E-01;
    COFD[1243] = 1.59879891E-02;
    COFD[1244] = -1.40949196E+01;
    COFD[1245] = 3.07549274E+00;
    COFD[1246] = -1.88889344E-01;
    COFD[1247] = 8.37152866E-03;
    COFD[1248] = -1.95739570E+01;
    COFD[1249] = 5.61113230E+00;
    COFD[1250] = -4.90190187E-01;
    COFD[1251] = 2.03260675E-02;
    COFD[1252] = -1.71982995E+01;
    COFD[1253] = 4.63881404E+00;
    COFD[1254] = -3.86139633E-01;
    COFD[1255] = 1.66955081E-02;
    COFD[1256] = -2.12639214E+01;
    COFD[1257] = 5.61184117E+00;
    COFD[1258] = -4.90532156E-01;
    COFD[1259] = 2.03507922E-02;
    COFD[1260] = -2.12639214E+01;
    COFD[1261] = 5.61184117E+00;
    COFD[1262] = -4.90532156E-01;
    COFD[1263] = 2.03507922E-02;
    COFD[1264] = -2.12831323E+01;
    COFD[1265] = 5.61184117E+00;
    COFD[1266] = -4.90532156E-01;
    COFD[1267] = 2.03507922E-02;
    COFD[1268] = -1.94093572E+01;
    COFD[1269] = 5.16013126E+00;
    COFD[1270] = -4.46824543E-01;
    COFD[1271] = 1.90464887E-02;
    COFD[1272] = -2.14087397E+01;
    COFD[1273] = 5.57282008E+00;
    COFD[1274] = -4.76690890E-01;
    COFD[1275] = 1.94000719E-02;
    COFD[1276] = -1.94253036E+01;
    COFD[1277] = 5.16013126E+00;
    COFD[1278] = -4.46824543E-01;
    COFD[1279] = 1.90464887E-02;
    COFD[1280] = -1.19157919E+01;
    COFD[1281] = 9.28955130E-01;
    COFD[1282] = 2.42107090E-01;
    COFD[1283] = -1.59823963E-02;
    COFD[1284] = -2.09565916E+01;
    COFD[1285] = 5.18380539E+00;
    COFD[1286] = -4.06234719E-01;
    COFD[1287] = 1.55515345E-02;
    COFD[1288] = -2.09642705E+01;
    COFD[1289] = 5.18380539E+00;
    COFD[1290] = -4.06234719E-01;
    COFD[1291] = 1.55515345E-02;
    COFD[1292] = -2.11388331E+01;
    COFD[1293] = 5.55529675E+00;
    COFD[1294] = -4.87942518E-01;
    COFD[1295] = 2.04249054E-02;
    COFD[1296] = -2.11309197E+01;
    COFD[1297] = 5.32644193E+00;
    COFD[1298] = -4.30581064E-01;
    COFD[1299] = 1.68379725E-02;
    COFD[1300] = -1.77498543E+01;
    COFD[1301] = 3.57475686E+00;
    COFD[1302] = -1.56396297E-01;
    COFD[1303] = 3.12157721E-03;
    COFD[1304] = -2.13084334E+01;
    COFD[1305] = 5.27210469E+00;
    COFD[1306] = -4.21419216E-01;
    COFD[1307] = 1.63567178E-02;
    COFD[1308] = -1.77563250E+01;
    COFD[1309] = 3.57475686E+00;
    COFD[1310] = -1.56396297E-01;
    COFD[1311] = 3.12157721E-03;
    COFD[1312] = -2.13148887E+01;
    COFD[1313] = 5.27210469E+00;
    COFD[1314] = -4.21419216E-01;
    COFD[1315] = 1.63567178E-02;
    COFD[1316] = -1.65295288E+01;
    COFD[1317] = 2.97569206E+00;
    COFD[1318] = -6.75652842E-02;
    COFD[1319] = -1.08648422E-03;
    COFD[1320] = -1.65295288E+01;
    COFD[1321] = 2.97569206E+00;
    COFD[1322] = -6.75652842E-02;
    COFD[1323] = -1.08648422E-03;
    COFD[1324] = -2.12652533E+01;
    COFD[1325] = 5.59961818E+00;
    COFD[1326] = -4.91624858E-01;
    COFD[1327] = 2.05035550E-02;
    COFD[1328] = -2.06463744E+01;
    COFD[1329] = 5.41688482E+00;
    COFD[1330] = -4.73387188E-01;
    COFD[1331] = 1.99280175E-02;
    COFD[1332] = -2.06516336E+01;
    COFD[1333] = 5.41688482E+00;
    COFD[1334] = -4.73387188E-01;
    COFD[1335] = 1.99280175E-02;
    COFD[1336] = -2.07653719E+01;
    COFD[1337] = 5.01092022E+00;
    COFD[1338] = -3.77985635E-01;
    COFD[1339] = 1.40968645E-02;
    COFD[1340] = -1.87453067E+01;
    COFD[1341] = 3.96926341E+00;
    COFD[1342] = -2.16412264E-01;
    COFD[1343] = 6.06012078E-03;
    COFD[1344] = -2.05013151E+01;
    COFD[1345] = 4.74677479E+00;
    COFD[1346] = -3.36160335E-01;
    COFD[1347] = 1.19858600E-02;
    COFD[1348] = -1.84568379E+01;
    COFD[1349] = 3.75912079E+00;
    COFD[1350] = -1.84235105E-01;
    COFD[1351] = 4.47800951E-03;
    COFD[1352] = -2.05044494E+01;
    COFD[1353] = 4.74677479E+00;
    COFD[1354] = -3.36160335E-01;
    COFD[1355] = 1.19858600E-02;
    COFD[1356] = -1.80862867E+01;
    COFD[1357] = 3.69199168E+00;
    COFD[1358] = -1.74005516E-01;
    COFD[1359] = 3.97694372E-03;
    COFD[1360] = -1.95040271E+01;
    COFD[1361] = 4.21243642E+00;
    COFD[1362] = -2.53087979E-01;
    COFD[1363] = 7.84955719E-03;
    COFD[1364] = -1.95059929E+01;
    COFD[1365] = 4.21243642E+00;
    COFD[1366] = -2.53087979E-01;
    COFD[1367] = 7.84955719E-03;
    COFD[1368] = -1.86141216E+01;
    COFD[1369] = 3.75739003E+00;
    COFD[1370] = -1.83970343E-01;
    COFD[1371] = 4.46501245E-03;
    COFD[1372] = -2.05420607E+01;
    COFD[1373] = 4.74677479E+00;
    COFD[1374] = -3.36160335E-01;
    COFD[1375] = 1.19858600E-02;
    COFD[1376] = -2.05618968E+01;
    COFD[1377] = 4.74677479E+00;
    COFD[1378] = -3.36160335E-01;
    COFD[1379] = 1.19858600E-02;
    COFD[1380] = -2.05618968E+01;
    COFD[1381] = 4.74677479E+00;
    COFD[1382] = -3.36160335E-01;
    COFD[1383] = 1.19858600E-02;
    COFD[1384] = -2.05644355E+01;
    COFD[1385] = 4.74677479E+00;
    COFD[1386] = -3.36160335E-01;
    COFD[1387] = 1.19858600E-02;
    COFD[1388] = -2.05644355E+01;
    COFD[1389] = 4.74677479E+00;
    COFD[1390] = -3.36160335E-01;
    COFD[1391] = 1.19858600E-02;
    COFD[1392] = -2.05801081E+01;
    COFD[1393] = 4.74677479E+00;
    COFD[1394] = -3.36160335E-01;
    COFD[1395] = 1.19858600E-02;
    COFD[1396] = -2.05930361E+01;
    COFD[1397] = 4.74677479E+00;
    COFD[1398] = -3.36160335E-01;
    COFD[1399] = 1.19858600E-02;
    COFD[1400] = -2.10643259E+01;
    COFD[1401] = 5.53614847E+00;
    COFD[1402] = -4.86046736E-01;
    COFD[1403] = 2.03659188E-02;
    COFD[1404] = -1.79310765E+01;
    COFD[1405] = 4.98037650E+00;
    COFD[1406] = -4.26676911E-01;
    COFD[1407] = 1.83007231E-02;
    COFD[1408] = -1.39315266E+01;
    COFD[1409] = 3.30394764E+00;
    COFD[1410] = -2.17920112E-01;
    COFD[1411] = 9.60284243E-03;
    COFD[1412] = -1.95548230E+01;
    COFD[1413] = 4.97133070E+00;
    COFD[1414] = -4.25604177E-01;
    COFD[1415] = 1.82582594E-02;
    COFD[1416] = -1.95548230E+01;
    COFD[1417] = 4.97133070E+00;
    COFD[1418] = -4.25604177E-01;
    COFD[1419] = 1.82582594E-02;
    COFD[1420] = -1.95770968E+01;
    COFD[1421] = 4.97133070E+00;
    COFD[1422] = -4.25604177E-01;
    COFD[1423] = 1.82582594E-02;
    COFD[1424] = -1.72286007E+01;
    COFD[1425] = 4.24084025E+00;
    COFD[1426] = -3.37428619E-01;
    COFD[1427] = 1.47032793E-02;
    COFD[1428] = -1.95154079E+01;
    COFD[1429] = 4.94787350E+00;
    COFD[1430] = -4.22829292E-01;
    COFD[1431] = 1.81487163E-02;
    COFD[1432] = -1.72473011E+01;
    COFD[1433] = 4.24084025E+00;
    COFD[1434] = -3.37428619E-01;
    COFD[1435] = 1.47032793E-02;
    COFD[1436] = -2.09565916E+01;
    COFD[1437] = 5.18380539E+00;
    COFD[1438] = -4.06234719E-01;
    COFD[1439] = 1.55515345E-02;
    COFD[1440] = -2.15453676E+01;
    COFD[1441] = 5.55313619E+00;
    COFD[1442] = -4.87753729E-01;
    COFD[1443] = 2.04203421E-02;
    COFD[1444] = -2.15547727E+01;
    COFD[1445] = 5.55313619E+00;
    COFD[1446] = -4.87753729E-01;
    COFD[1447] = 2.04203421E-02;
    COFD[1448] = -1.83296965E+01;
    COFD[1449] = 4.48570999E+00;
    COFD[1450] = -3.67301524E-01;
    COFD[1451] = 1.59204254E-02;
    COFD[1452] = -2.11427744E+01;
    COFD[1453] = 5.43893233E+00;
    COFD[1454] = -4.75546039E-01;
    COFD[1455] = 1.99938690E-02;
    COFD[1456] = -2.16718247E+01;
    COFD[1457] = 5.36811769E+00;
    COFD[1458] = -4.37727086E-01;
    COFD[1459] = 1.72167686E-02;
    COFD[1460] = -2.15126310E+01;
    COFD[1461] = 5.48426911E+00;
    COFD[1462] = -4.80606512E-01;
    COFD[1463] = 2.01811046E-02;
    COFD[1464] = -2.16798265E+01;
    COFD[1465] = 5.36811769E+00;
    COFD[1466] = -4.37727086E-01;
    COFD[1467] = 1.72167686E-02;
    COFD[1468] = -2.15206146E+01;
    COFD[1469] = 5.48426911E+00;
    COFD[1470] = -4.80606512E-01;
    COFD[1471] = 2.01811046E-02;
    COFD[1472] = -2.18848136E+01;
    COFD[1473] = 5.51302074E+00;
    COFD[1474] = -4.65263979E-01;
    COFD[1475] = 1.87580679E-02;
    COFD[1476] = -2.18848136E+01;
    COFD[1477] = 5.51302074E+00;
    COFD[1478] = -4.65263979E-01;
    COFD[1479] = 1.87580679E-02;
    COFD[1480] = -1.86507213E+01;
    COFD[1481] = 4.60874797E+00;
    COFD[1482] = -3.82368716E-01;
    COFD[1483] = 1.65370164E-02;
    COFD[1484] = -1.86576191E+01;
    COFD[1485] = 4.60874797E+00;
    COFD[1486] = -3.82368716E-01;
    COFD[1487] = 1.65370164E-02;
    COFD[1488] = -1.86641962E+01;
    COFD[1489] = 4.60874797E+00;
    COFD[1490] = -3.82368716E-01;
    COFD[1491] = 1.65370164E-02;
    COFD[1492] = -2.13961414E+01;
    COFD[1493] = 5.46685775E+00;
    COFD[1494] = -4.78665416E-01;
    COFD[1495] = 2.01093915E-02;
    COFD[1496] = -2.20351084E+01;
    COFD[1497] = 5.49663315E+00;
    COFD[1498] = -4.61182837E-01;
    COFD[1499] = 1.85035558E-02;
    COFD[1500] = -2.21632697E+01;
    COFD[1501] = 5.59047891E+00;
    COFD[1502] = -4.85359237E-01;
    COFD[1503] = 2.00302091E-02;
    COFD[1504] = -2.20585906E+01;
    COFD[1505] = 5.42445100E+00;
    COFD[1506] = -4.47918761E-01;
    COFD[1507] = 1.77729995E-02;
    COFD[1508] = -2.21672921E+01;
    COFD[1509] = 5.59047891E+00;
    COFD[1510] = -4.85359237E-01;
    COFD[1511] = 2.00302091E-02;
    COFD[1512] = -2.18318278E+01;
    COFD[1513] = 5.40298848E+00;
    COFD[1514] = -4.43954594E-01;
    COFD[1515] = 1.75542998E-02;
    COFD[1516] = -2.22849472E+01;
    COFD[1517] = 5.53187091E+00;
    COFD[1518] = -4.68918850E-01;
    COFD[1519] = 1.89648616E-02;
    COFD[1520] = -2.22875222E+01;
    COFD[1521] = 5.53187091E+00;
    COFD[1522] = -4.68918850E-01;
    COFD[1523] = 1.89648616E-02;
    COFD[1524] = -2.22125735E+01;
    COFD[1525] = 5.42385157E+00;
    COFD[1526] = -4.47809271E-01;
    COFD[1527] = 1.77669962E-02;
    COFD[1528] = -2.22161430E+01;
    COFD[1529] = 5.59047891E+00;
    COFD[1530] = -4.85359237E-01;
    COFD[1531] = 2.00302091E-02;
    COFD[1532] = -2.22423680E+01;
    COFD[1533] = 5.59047891E+00;
    COFD[1534] = -4.85359237E-01;
    COFD[1535] = 2.00302091E-02;
    COFD[1536] = -2.22423680E+01;
    COFD[1537] = 5.59047891E+00;
    COFD[1538] = -4.85359237E-01;
    COFD[1539] = 2.00302091E-02;
    COFD[1540] = -2.22457488E+01;
    COFD[1541] = 5.59047891E+00;
    COFD[1542] = -4.85359237E-01;
    COFD[1543] = 2.00302091E-02;
    COFD[1544] = -2.22457488E+01;
    COFD[1545] = 5.59047891E+00;
    COFD[1546] = -4.85359237E-01;
    COFD[1547] = 2.00302091E-02;
    COFD[1548] = -2.22667492E+01;
    COFD[1549] = 5.59047891E+00;
    COFD[1550] = -4.85359237E-01;
    COFD[1551] = 2.00302091E-02;
    COFD[1552] = -2.22842444E+01;
    COFD[1553] = 5.59047891E+00;
    COFD[1554] = -4.85359237E-01;
    COFD[1555] = 2.00302091E-02;
    COFD[1556] = -1.83039618E+01;
    COFD[1557] = 4.47952077E+00;
    COFD[1558] = -3.66569471E-01;
    COFD[1559] = 1.58916129E-02;
    COFD[1560] = -1.79317714E+01;
    COFD[1561] = 4.98037650E+00;
    COFD[1562] = -4.26676911E-01;
    COFD[1563] = 1.83007231E-02;
    COFD[1564] = -1.39328674E+01;
    COFD[1565] = 3.30394764E+00;
    COFD[1566] = -2.17920112E-01;
    COFD[1567] = 9.60284243E-03;
    COFD[1568] = -1.95613899E+01;
    COFD[1569] = 4.97133070E+00;
    COFD[1570] = -4.25604177E-01;
    COFD[1571] = 1.82582594E-02;
    COFD[1572] = -1.95613899E+01;
    COFD[1573] = 4.97133070E+00;
    COFD[1574] = -4.25604177E-01;
    COFD[1575] = 1.82582594E-02;
    COFD[1576] = -1.95839648E+01;
    COFD[1577] = 4.97133070E+00;
    COFD[1578] = -4.25604177E-01;
    COFD[1579] = 1.82582594E-02;
    COFD[1580] = -1.72357436E+01;
    COFD[1581] = 4.24084025E+00;
    COFD[1582] = -3.37428619E-01;
    COFD[1583] = 1.47032793E-02;
    COFD[1584] = -1.95225629E+01;
    COFD[1585] = 4.94787350E+00;
    COFD[1586] = -4.22829292E-01;
    COFD[1587] = 1.81487163E-02;
    COFD[1588] = -1.72547182E+01;
    COFD[1589] = 4.24084025E+00;
    COFD[1590] = -3.37428619E-01;
    COFD[1591] = 1.47032793E-02;
    COFD[1592] = -2.09642705E+01;
    COFD[1593] = 5.18380539E+00;
    COFD[1594] = -4.06234719E-01;
    COFD[1595] = 1.55515345E-02;
    COFD[1596] = -2.15547727E+01;
    COFD[1597] = 5.55313619E+00;
    COFD[1598] = -4.87753729E-01;
    COFD[1599] = 2.04203421E-02;
    COFD[1600] = -2.15643580E+01;
    COFD[1601] = 5.55313619E+00;
    COFD[1602] = -4.87753729E-01;
    COFD[1603] = 2.04203421E-02;
    COFD[1604] = -1.83394481E+01;
    COFD[1605] = 4.48570999E+00;
    COFD[1606] = -3.67301524E-01;
    COFD[1607] = 1.59204254E-02;
    COFD[1608] = -2.11525334E+01;
    COFD[1609] = 5.43893233E+00;
    COFD[1610] = -4.75546039E-01;
    COFD[1611] = 1.99938690E-02;
    COFD[1612] = -2.16817439E+01;
    COFD[1613] = 5.36811769E+00;
    COFD[1614] = -4.37727086E-01;
    COFD[1615] = 1.72167686E-02;
    COFD[1616] = -2.15225573E+01;
    COFD[1617] = 5.48426911E+00;
    COFD[1618] = -4.80606512E-01;
    COFD[1619] = 2.01811046E-02;
    COFD[1620] = -2.16899073E+01;
    COFD[1621] = 5.36811769E+00;
    COFD[1622] = -4.37727086E-01;
    COFD[1623] = 1.72167686E-02;
    COFD[1624] = -2.15307023E+01;
    COFD[1625] = 5.48426911E+00;
    COFD[1626] = -4.80606512E-01;
    COFD[1627] = 2.01811046E-02;
    COFD[1628] = -2.18950505E+01;
    COFD[1629] = 5.51302074E+00;
    COFD[1630] = -4.65263979E-01;
    COFD[1631] = 1.87580679E-02;
    COFD[1632] = -2.18950505E+01;
    COFD[1633] = 5.51302074E+00;
    COFD[1634] = -4.65263979E-01;
    COFD[1635] = 1.87580679E-02;
    COFD[1636] = -1.86611023E+01;
    COFD[1637] = 4.60874797E+00;
    COFD[1638] = -3.82368716E-01;
    COFD[1639] = 1.65370164E-02;
    COFD[1640] = -1.86681459E+01;
    COFD[1641] = 4.60874797E+00;
    COFD[1642] = -3.82368716E-01;
    COFD[1643] = 1.65370164E-02;
    COFD[1644] = -1.86748638E+01;
    COFD[1645] = 4.60874797E+00;
    COFD[1646] = -3.82368716E-01;
    COFD[1647] = 1.65370164E-02;
    COFD[1648] = -2.14079882E+01;
    COFD[1649] = 5.46685775E+00;
    COFD[1650] = -4.78665416E-01;
    COFD[1651] = 2.01093915E-02;
    COFD[1652] = -2.20469595E+01;
    COFD[1653] = 5.49663315E+00;
    COFD[1654] = -4.61182837E-01;
    COFD[1655] = 1.85035558E-02;
    COFD[1656] = -2.21752213E+01;
    COFD[1657] = 5.59047891E+00;
    COFD[1658] = -4.85359237E-01;
    COFD[1659] = 2.00302091E-02;
    COFD[1660] = -2.20706358E+01;
    COFD[1661] = 5.42445100E+00;
    COFD[1662] = -4.47918761E-01;
    COFD[1663] = 1.77729995E-02;
    COFD[1664] = -2.21793415E+01;
    COFD[1665] = 5.59047891E+00;
    COFD[1666] = -4.85359237E-01;
    COFD[1667] = 2.00302091E-02;
    COFD[1668] = -2.18439681E+01;
    COFD[1669] = 5.40298848E+00;
    COFD[1670] = -4.43954594E-01;
    COFD[1671] = 1.75542998E-02;
    COFD[1672] = -2.22980489E+01;
    COFD[1673] = 5.53187091E+00;
    COFD[1674] = -4.68918850E-01;
    COFD[1675] = 1.89648616E-02;
    COFD[1676] = -2.23006924E+01;
    COFD[1677] = 5.53187091E+00;
    COFD[1678] = -4.68918850E-01;
    COFD[1679] = 1.89648616E-02;
    COFD[1680] = -2.22258107E+01;
    COFD[1681] = 5.42385157E+00;
    COFD[1682] = -4.47809271E-01;
    COFD[1683] = 1.77669962E-02;
    COFD[1684] = -2.22294456E+01;
    COFD[1685] = 5.59047891E+00;
    COFD[1686] = -4.85359237E-01;
    COFD[1687] = 2.00302091E-02;
    COFD[1688] = -2.22563971E+01;
    COFD[1689] = 5.59047891E+00;
    COFD[1690] = -4.85359237E-01;
    COFD[1691] = 2.00302091E-02;
    COFD[1692] = -2.22563971E+01;
    COFD[1693] = 5.59047891E+00;
    COFD[1694] = -4.85359237E-01;
    COFD[1695] = 2.00302091E-02;
    COFD[1696] = -2.22598745E+01;
    COFD[1697] = 5.59047891E+00;
    COFD[1698] = -4.85359237E-01;
    COFD[1699] = 2.00302091E-02;
    COFD[1700] = -2.22598745E+01;
    COFD[1701] = 5.59047891E+00;
    COFD[1702] = -4.85359237E-01;
    COFD[1703] = 2.00302091E-02;
    COFD[1704] = -2.22814898E+01;
    COFD[1705] = 5.59047891E+00;
    COFD[1706] = -4.85359237E-01;
    COFD[1707] = 2.00302091E-02;
    COFD[1708] = -2.22995181E+01;
    COFD[1709] = 5.59047891E+00;
    COFD[1710] = -4.85359237E-01;
    COFD[1711] = 2.00302091E-02;
    COFD[1712] = -1.83137139E+01;
    COFD[1713] = 4.47952077E+00;
    COFD[1714] = -3.66569471E-01;
    COFD[1715] = 1.58916129E-02;
    COFD[1716] = -1.43151174E+01;
    COFD[1717] = 3.68038508E+00;
    COFD[1718] = -2.65779346E-01;
    COFD[1719] = 1.16360771E-02;
    COFD[1720] = -1.17159737E+01;
    COFD[1721] = 2.48123210E+00;
    COFD[1722] = -1.11322604E-01;
    COFD[1723] = 4.99282389E-03;
    COFD[1724] = -1.59634533E+01;
    COFD[1725] = 3.67388294E+00;
    COFD[1726] = -2.64990709E-01;
    COFD[1727] = 1.16042706E-02;
    COFD[1728] = -1.59634533E+01;
    COFD[1729] = 3.67388294E+00;
    COFD[1730] = -2.64990709E-01;
    COFD[1731] = 1.16042706E-02;
    COFD[1732] = -1.59863030E+01;
    COFD[1733] = 3.67388294E+00;
    COFD[1734] = -2.64990709E-01;
    COFD[1735] = 1.16042706E-02;
    COFD[1736] = -1.40999008E+01;
    COFD[1737] = 3.08120012E+00;
    COFD[1738] = -1.89629903E-01;
    COFD[1739] = 8.40361952E-03;
    COFD[1740] = -1.59525102E+01;
    COFD[1741] = 3.66023858E+00;
    COFD[1742] = -2.63401043E-01;
    COFD[1743] = 1.15432000E-02;
    COFD[1744] = -1.41191261E+01;
    COFD[1745] = 3.08120012E+00;
    COFD[1746] = -1.89629903E-01;
    COFD[1747] = 8.40361952E-03;
    COFD[1748] = -2.11388331E+01;
    COFD[1749] = 5.55529675E+00;
    COFD[1750] = -4.87942518E-01;
    COFD[1751] = 2.04249054E-02;
    COFD[1752] = -1.83296965E+01;
    COFD[1753] = 4.48570999E+00;
    COFD[1754] = -3.67301524E-01;
    COFD[1755] = 1.59204254E-02;
    COFD[1756] = -1.83394481E+01;
    COFD[1757] = 4.48570999E+00;
    COFD[1758] = -3.67301524E-01;
    COFD[1759] = 1.59204254E-02;
    COFD[1760] = -1.50233475E+01;
    COFD[1761] = 3.26660767E+00;
    COFD[1762] = -2.13287177E-01;
    COFD[1763] = 9.41137857E-03;
    COFD[1764] = -1.79116531E+01;
    COFD[1765] = 4.35148286E+00;
    COFD[1766] = -3.50886647E-01;
    COFD[1767] = 1.52498573E-02;
    COFD[1768] = -2.05045578E+01;
    COFD[1769] = 5.23843909E+00;
    COFD[1770] = -4.55815614E-01;
    COFD[1771] = 1.93898040E-02;
    COFD[1772] = -1.82872310E+01;
    COFD[1773] = 4.40289649E+00;
    COFD[1774] = -3.57289765E-01;
    COFD[1775] = 1.55166804E-02;
    COFD[1776] = -2.05128705E+01;
    COFD[1777] = 5.23843909E+00;
    COFD[1778] = -4.55815614E-01;
    COFD[1779] = 1.93898040E-02;
    COFD[1780] = -1.82955252E+01;
    COFD[1781] = 4.40289649E+00;
    COFD[1782] = -3.57289765E-01;
    COFD[1783] = 1.55166804E-02;
    COFD[1784] = -2.02642227E+01;
    COFD[1785] = 5.14499740E+00;
    COFD[1786] = -4.45694430E-01;
    COFD[1787] = 1.90318646E-02;
    COFD[1788] = -2.02642227E+01;
    COFD[1789] = 5.14499740E+00;
    COFD[1790] = -4.45694430E-01;
    COFD[1791] = 1.90318646E-02;
    COFD[1792] = -1.52721107E+01;
    COFD[1793] = 3.36790500E+00;
    COFD[1794] = -2.26321740E-01;
    COFD[1795] = 9.97135055E-03;
    COFD[1796] = -1.52792891E+01;
    COFD[1797] = 3.36790500E+00;
    COFD[1798] = -2.26321740E-01;
    COFD[1799] = 9.97135055E-03;
    COFD[1800] = -1.52861376E+01;
    COFD[1801] = 3.36790500E+00;
    COFD[1802] = -2.26321740E-01;
    COFD[1803] = 9.97135055E-03;
    COFD[1804] = -1.81735763E+01;
    COFD[1805] = 4.38391495E+00;
    COFD[1806] = -3.54941287E-01;
    COFD[1807] = 1.54195107E-02;
    COFD[1808] = -2.03015042E+01;
    COFD[1809] = 5.11106992E+00;
    COFD[1810] = -4.42047129E-01;
    COFD[1811] = 1.89042990E-02;
    COFD[1812] = -1.94622154E+01;
    COFD[1813] = 4.76145602E+00;
    COFD[1814] = -4.00684587E-01;
    COFD[1815] = 1.72708322E-02;
    COFD[1816] = -2.06106760E+01;
    COFD[1817] = 5.16748146E+00;
    COFD[1818] = -4.47594939E-01;
    COFD[1819] = 1.90724110E-02;
    COFD[1820] = -1.94664266E+01;
    COFD[1821] = 4.76145602E+00;
    COFD[1822] = -4.00684587E-01;
    COFD[1823] = 1.72708322E-02;
    COFD[1824] = -2.04949373E+01;
    COFD[1825] = 5.19614628E+00;
    COFD[1826] = -4.50889164E-01;
    COFD[1827] = 1.91983328E-02;
    COFD[1828] = -2.02929645E+01;
    COFD[1829] = 5.02679716E+00;
    COFD[1830] = -4.32175153E-01;
    COFD[1831] = 1.85182518E-02;
    COFD[1832] = -2.02956721E+01;
    COFD[1833] = 5.02679716E+00;
    COFD[1834] = -4.32175153E-01;
    COFD[1835] = 1.85182518E-02;
    COFD[1836] = -2.07708766E+01;
    COFD[1837] = 5.16820334E+00;
    COFD[1838] = -4.47675828E-01;
    COFD[1839] = 1.90754011E-02;
    COFD[1840] = -1.95177006E+01;
    COFD[1841] = 4.76145602E+00;
    COFD[1842] = -4.00684587E-01;
    COFD[1843] = 1.72708322E-02;
    COFD[1844] = -1.95453329E+01;
    COFD[1845] = 4.76145602E+00;
    COFD[1846] = -4.00684587E-01;
    COFD[1847] = 1.72708322E-02;
    COFD[1848] = -1.95453329E+01;
    COFD[1849] = 4.76145602E+00;
    COFD[1850] = -4.00684587E-01;
    COFD[1851] = 1.72708322E-02;
    COFD[1852] = -1.95489008E+01;
    COFD[1853] = 4.76145602E+00;
    COFD[1854] = -4.00684587E-01;
    COFD[1855] = 1.72708322E-02;
    COFD[1856] = -1.95489008E+01;
    COFD[1857] = 4.76145602E+00;
    COFD[1858] = -4.00684587E-01;
    COFD[1859] = 1.72708322E-02;
    COFD[1860] = -1.95710942E+01;
    COFD[1861] = 4.76145602E+00;
    COFD[1862] = -4.00684587E-01;
    COFD[1863] = 1.72708322E-02;
    COFD[1864] = -1.95896245E+01;
    COFD[1865] = 4.76145602E+00;
    COFD[1866] = -4.00684587E-01;
    COFD[1867] = 1.72708322E-02;
    COFD[1868] = -1.50031687E+01;
    COFD[1869] = 3.26223357E+00;
    COFD[1870] = -2.12746642E-01;
    COFD[1871] = 9.38912883E-03;
    COFD[1872] = -1.74407963E+01;
    COFD[1873] = 4.83580036E+00;
    COFD[1874] = -4.09383573E-01;
    COFD[1875] = 1.76098175E-02;
    COFD[1876] = -1.36336373E+01;
    COFD[1877] = 3.22088176E+00;
    COFD[1878] = -2.07623790E-01;
    COFD[1879] = 9.17771542E-03;
    COFD[1880] = -1.90996795E+01;
    COFD[1881] = 4.82869066E+00;
    COFD[1882] = -4.08564514E-01;
    COFD[1883] = 1.75784675E-02;
    COFD[1884] = -1.90996795E+01;
    COFD[1885] = 4.82869066E+00;
    COFD[1886] = -4.08564514E-01;
    COFD[1887] = 1.75784675E-02;
    COFD[1888] = -1.91225414E+01;
    COFD[1889] = 4.82869066E+00;
    COFD[1890] = -4.08564514E-01;
    COFD[1891] = 1.75784675E-02;
    COFD[1892] = -1.68343393E+01;
    COFD[1893] = 4.11954900E+00;
    COFD[1894] = -3.22470391E-01;
    COFD[1895] = 1.40859564E-02;
    COFD[1896] = -1.90692595E+01;
    COFD[1897] = 4.80830699E+00;
    COFD[1898] = -4.06171933E-01;
    COFD[1899] = 1.74848791E-02;
    COFD[1900] = -1.68535757E+01;
    COFD[1901] = 4.11954900E+00;
    COFD[1902] = -3.22470391E-01;
    COFD[1903] = 1.40859564E-02;
    COFD[1904] = -2.11309197E+01;
    COFD[1905] = 5.32644193E+00;
    COFD[1906] = -4.30581064E-01;
    COFD[1907] = 1.68379725E-02;
    COFD[1908] = -2.11427744E+01;
    COFD[1909] = 5.43893233E+00;
    COFD[1910] = -4.75546039E-01;
    COFD[1911] = 1.99938690E-02;
    COFD[1912] = -2.11525334E+01;
    COFD[1913] = 5.43893233E+00;
    COFD[1914] = -4.75546039E-01;
    COFD[1915] = 1.99938690E-02;
    COFD[1916] = -1.79116531E+01;
    COFD[1917] = 4.35148286E+00;
    COFD[1918] = -3.50886647E-01;
    COFD[1919] = 1.52498573E-02;
    COFD[1920] = -2.08820897E+01;
    COFD[1921] = 5.38250415E+00;
    COFD[1922] = -4.71144140E-01;
    COFD[1923] = 1.99199779E-02;
    COFD[1924] = -2.17771954E+01;
    COFD[1925] = 5.47519298E+00;
    COFD[1926] = -4.57113040E-01;
    COFD[1927] = 1.82758312E-02;
    COFD[1928] = -2.11931178E+01;
    COFD[1929] = 5.40060531E+00;
    COFD[1930] = -4.72449699E-01;
    COFD[1931] = 1.99345817E-02;
    COFD[1932] = -2.17855148E+01;
    COFD[1933] = 5.47519298E+00;
    COFD[1934] = -4.57113040E-01;
    COFD[1935] = 1.82758312E-02;
    COFD[1936] = -2.12014186E+01;
    COFD[1937] = 5.40060531E+00;
    COFD[1938] = -4.72449699E-01;
    COFD[1939] = 1.99345817E-02;
    COFD[1940] = -2.19136842E+01;
    COFD[1941] = 5.58503445E+00;
    COFD[1942] = -4.79552117E-01;
    COFD[1943] = 1.95750393E-02;
    COFD[1944] = -2.19136842E+01;
    COFD[1945] = 5.58503445E+00;
    COFD[1946] = -4.79552117E-01;
    COFD[1947] = 1.95750393E-02;
    COFD[1948] = -1.82145353E+01;
    COFD[1949] = 4.46848269E+00;
    COFD[1950] = -3.65269718E-01;
    COFD[1951] = 1.58407652E-02;
    COFD[1952] = -1.82217198E+01;
    COFD[1953] = 4.46848269E+00;
    COFD[1954] = -3.65269718E-01;
    COFD[1955] = 1.58407652E-02;
    COFD[1956] = -1.82285740E+01;
    COFD[1957] = 4.46848269E+00;
    COFD[1958] = -3.65269718E-01;
    COFD[1959] = 1.58407652E-02;
    COFD[1960] = -2.11031143E+01;
    COFD[1961] = 5.39439999E+00;
    COFD[1962] = -4.72050184E-01;
    COFD[1963] = 1.99336257E-02;
    COFD[1964] = -2.20490757E+01;
    COFD[1965] = 5.56049839E+00;
    COFD[1966] = -4.74367872E-01;
    COFD[1967] = 1.92702787E-02;
    COFD[1968] = -2.20752743E+01;
    COFD[1969] = 5.60509076E+00;
    COFD[1970] = -4.91228646E-01;
    COFD[1971] = 2.04425499E-02;
    COFD[1972] = -2.21129692E+01;
    COFD[1973] = 5.50506115E+00;
    COFD[1974] = -4.63563533E-01;
    COFD[1975] = 1.86575247E-02;
    COFD[1976] = -2.20794895E+01;
    COFD[1977] = 5.60509076E+00;
    COFD[1978] = -4.91228646E-01;
    COFD[1979] = 2.04425499E-02;
    COFD[1980] = -2.19162360E+01;
    COFD[1981] = 5.49906960E+00;
    COFD[1982] = -4.61793001E-01;
    COFD[1983] = 1.85415189E-02;
    COFD[1984] = -2.22838940E+01;
    COFD[1985] = 5.58518325E+00;
    COFD[1986] = -4.80534710E-01;
    COFD[1987] = 1.96556785E-02;
    COFD[1988] = -2.22866045E+01;
    COFD[1989] = 5.58518325E+00;
    COFD[1990] = -4.80534710E-01;
    COFD[1991] = 1.96556785E-02;
    COFD[1992] = -2.22731582E+01;
    COFD[1993] = 5.50482046E+00;
    COFD[1994] = -4.63503431E-01;
    COFD[1995] = 1.86537631E-02;
    COFD[1996] = -2.21308159E+01;
    COFD[1997] = 5.60509076E+00;
    COFD[1998] = -4.91228646E-01;
    COFD[1999] = 2.04425499E-02;
    COFD[2000] = -2.21584786E+01;
    COFD[2001] = 5.60509076E+00;
    COFD[2002] = -4.91228646E-01;
    COFD[2003] = 2.04425499E-02;
    COFD[2004] = -2.21584786E+01;
    COFD[2005] = 5.60509076E+00;
    COFD[2006] = -4.91228646E-01;
    COFD[2007] = 2.04425499E-02;
    COFD[2008] = -2.21620506E+01;
    COFD[2009] = 5.60509076E+00;
    COFD[2010] = -4.91228646E-01;
    COFD[2011] = 2.04425499E-02;
    COFD[2012] = -2.21620506E+01;
    COFD[2013] = 5.60509076E+00;
    COFD[2014] = -4.91228646E-01;
    COFD[2015] = 2.04425499E-02;
    COFD[2016] = -2.21842699E+01;
    COFD[2017] = 5.60509076E+00;
    COFD[2018] = -4.91228646E-01;
    COFD[2019] = 2.04425499E-02;
    COFD[2020] = -2.22028227E+01;
    COFD[2021] = 5.60509076E+00;
    COFD[2022] = -4.91228646E-01;
    COFD[2023] = 2.04425499E-02;
    COFD[2024] = -1.78815889E+01;
    COFD[2025] = 4.34347890E+00;
    COFD[2026] = -3.49890003E-01;
    COFD[2027] = 1.52083459E-02;
    COFD[2028] = -1.97544450E+01;
    COFD[2029] = 5.56931926E+00;
    COFD[2030] = -4.89105511E-01;
    COFD[2031] = 2.04493129E-02;
    COFD[2032] = -1.60517370E+01;
    COFD[2033] = 4.11188603E+00;
    COFD[2034] = -3.21540884E-01;
    COFD[2035] = 1.40482564E-02;
    COFD[2036] = -2.14160703E+01;
    COFD[2037] = 5.56531152E+00;
    COFD[2038] = -4.88789821E-01;
    COFD[2039] = 2.04437116E-02;
    COFD[2040] = -2.14160703E+01;
    COFD[2041] = 5.56531152E+00;
    COFD[2042] = -4.88789821E-01;
    COFD[2043] = 2.04437116E-02;
    COFD[2044] = -2.14391943E+01;
    COFD[2045] = 5.56531152E+00;
    COFD[2046] = -4.88789821E-01;
    COFD[2047] = 2.04437116E-02;
    COFD[2048] = -1.94313116E+01;
    COFD[2049] = 5.02567894E+00;
    COFD[2050] = -4.32045169E-01;
    COFD[2051] = 1.85132214E-02;
    COFD[2052] = -2.14022336E+01;
    COFD[2053] = 5.55346617E+00;
    COFD[2054] = -4.87783156E-01;
    COFD[2055] = 2.04210886E-02;
    COFD[2056] = -1.94507876E+01;
    COFD[2057] = 5.02567894E+00;
    COFD[2058] = -4.32045169E-01;
    COFD[2059] = 1.85132214E-02;
    COFD[2060] = -1.77498543E+01;
    COFD[2061] = 3.57475686E+00;
    COFD[2062] = -1.56396297E-01;
    COFD[2063] = 3.12157721E-03;
    COFD[2064] = -2.16718247E+01;
    COFD[2065] = 5.36811769E+00;
    COFD[2066] = -4.37727086E-01;
    COFD[2067] = 1.72167686E-02;
    COFD[2068] = -2.16817439E+01;
    COFD[2069] = 5.36811769E+00;
    COFD[2070] = -4.37727086E-01;
    COFD[2071] = 1.72167686E-02;
    COFD[2072] = -2.05045578E+01;
    COFD[2073] = 5.23843909E+00;
    COFD[2074] = -4.55815614E-01;
    COFD[2075] = 1.93898040E-02;
    COFD[2076] = -2.17771954E+01;
    COFD[2077] = 5.47519298E+00;
    COFD[2078] = -4.57113040E-01;
    COFD[2079] = 1.82758312E-02;
    COFD[2080] = -1.90328712E+01;
    COFD[2081] = 3.99221757E+00;
    COFD[2082] = -2.19854880E-01;
    COFD[2083] = 6.22736279E-03;
    COFD[2084] = -2.19615570E+01;
    COFD[2085] = 5.43750833E+00;
    COFD[2086] = -4.50273329E-01;
    COFD[2087] = 1.79013718E-02;
    COFD[2088] = -1.90413348E+01;
    COFD[2089] = 3.99221757E+00;
    COFD[2090] = -2.19854880E-01;
    COFD[2091] = 6.22736279E-03;
    COFD[2092] = -2.19700018E+01;
    COFD[2093] = 5.43750833E+00;
    COFD[2094] = -4.50273329E-01;
    COFD[2095] = 1.79013718E-02;
    COFD[2096] = -2.01801667E+01;
    COFD[2097] = 4.53183330E+00;
    COFD[2098] = -3.02186760E-01;
    COFD[2099] = 1.02756490E-02;
    COFD[2100] = -2.01801667E+01;
    COFD[2101] = 4.53183330E+00;
    COFD[2102] = -3.02186760E-01;
    COFD[2103] = 1.02756490E-02;
    COFD[2104] = -2.08204449E+01;
    COFD[2105] = 5.35267674E+00;
    COFD[2106] = -4.69010505E-01;
    COFD[2107] = 1.98979152E-02;
    COFD[2108] = -2.08277598E+01;
    COFD[2109] = 5.35267674E+00;
    COFD[2110] = -4.69010505E-01;
    COFD[2111] = 1.98979152E-02;
    COFD[2112] = -2.08347403E+01;
    COFD[2113] = 5.35267674E+00;
    COFD[2114] = -4.69010505E-01;
    COFD[2115] = 1.98979152E-02;
    COFD[2116] = -2.19215555E+01;
    COFD[2117] = 5.45216133E+00;
    COFD[2118] = -4.52916925E-01;
    COFD[2119] = 1.80456400E-02;
    COFD[2120] = -2.01009366E+01;
    COFD[2121] = 4.41511629E+00;
    COFD[2122] = -2.84086963E-01;
    COFD[2123] = 9.37586971E-03;
    COFD[2124] = -2.14954618E+01;
    COFD[2125] = 5.05104142E+00;
    COFD[2126] = -3.84428177E-01;
    COFD[2127] = 1.44249285E-02;
    COFD[2128] = -1.97587116E+01;
    COFD[2129] = 4.18758010E+00;
    COFD[2130] = -2.49327776E-01;
    COFD[2131] = 7.66559103E-03;
    COFD[2132] = -2.14997655E+01;
    COFD[2133] = 5.05104142E+00;
    COFD[2134] = -3.84428177E-01;
    COFD[2135] = 1.44249285E-02;
    COFD[2136] = -1.93946947E+01;
    COFD[2137] = 4.10954793E+00;
    COFD[2138] = -2.37523329E-01;
    COFD[2139] = 7.08858141E-03;
    COFD[2140] = -2.06900561E+01;
    COFD[2141] = 4.59219167E+00;
    COFD[2142] = -3.11627566E-01;
    COFD[2143] = 1.07472587E-02;
    COFD[2144] = -2.06928292E+01;
    COFD[2145] = 4.59219167E+00;
    COFD[2146] = -3.11627566E-01;
    COFD[2147] = 1.07472587E-02;
    COFD[2148] = -1.99140853E+01;
    COFD[2149] = 4.18537863E+00;
    COFD[2150] = -2.48994776E-01;
    COFD[2151] = 7.64930102E-03;
    COFD[2152] = -2.15522331E+01;
    COFD[2153] = 5.05104142E+00;
    COFD[2154] = -3.84428177E-01;
    COFD[2155] = 1.44249285E-02;
    COFD[2156] = -2.15805625E+01;
    COFD[2157] = 5.05104142E+00;
    COFD[2158] = -3.84428177E-01;
    COFD[2159] = 1.44249285E-02;
    COFD[2160] = -2.15805625E+01;
    COFD[2161] = 5.05104142E+00;
    COFD[2162] = -3.84428177E-01;
    COFD[2163] = 1.44249285E-02;
    COFD[2164] = -2.15842234E+01;
    COFD[2165] = 5.05104142E+00;
    COFD[2166] = -3.84428177E-01;
    COFD[2167] = 1.44249285E-02;
    COFD[2168] = -2.15842234E+01;
    COFD[2169] = 5.05104142E+00;
    COFD[2170] = -3.84428177E-01;
    COFD[2171] = 1.44249285E-02;
    COFD[2172] = -2.16070103E+01;
    COFD[2173] = 5.05104142E+00;
    COFD[2174] = -3.84428177E-01;
    COFD[2175] = 1.44249285E-02;
    COFD[2176] = -2.16260574E+01;
    COFD[2177] = 5.05104142E+00;
    COFD[2178] = -3.84428177E-01;
    COFD[2179] = 1.44249285E-02;
    COFD[2180] = -2.04750581E+01;
    COFD[2181] = 5.23112374E+00;
    COFD[2182] = -4.54967682E-01;
    COFD[2183] = 1.93570423E-02;
    COFD[2184] = -1.78631557E+01;
    COFD[2185] = 4.88268692E+00;
    COFD[2186] = -4.14917638E-01;
    COFD[2187] = 1.78274298E-02;
    COFD[2188] = -1.39648112E+01;
    COFD[2189] = 3.24966086E+00;
    COFD[2190] = -2.11199992E-01;
    COFD[2191] = 9.32580661E-03;
    COFD[2192] = -1.94530250E+01;
    COFD[2193] = 4.87180830E+00;
    COFD[2194] = -4.13582958E-01;
    COFD[2195] = 1.77726094E-02;
    COFD[2196] = -1.94530250E+01;
    COFD[2197] = 4.87180830E+00;
    COFD[2198] = -4.13582958E-01;
    COFD[2199] = 1.77726094E-02;
    COFD[2200] = -1.94761606E+01;
    COFD[2201] = 4.87180830E+00;
    COFD[2202] = -4.13582958E-01;
    COFD[2203] = 1.77726094E-02;
    COFD[2204] = -1.72053106E+01;
    COFD[2205] = 4.15807461E+00;
    COFD[2206] = -3.27178539E-01;
    COFD[2207] = 1.42784349E-02;
    COFD[2208] = -1.94126575E+01;
    COFD[2209] = 4.84669430E+00;
    COFD[2210] = -4.10571455E-01;
    COFD[2211] = 1.76520543E-02;
    COFD[2212] = -1.72247972E+01;
    COFD[2213] = 4.15807461E+00;
    COFD[2214] = -3.27178539E-01;
    COFD[2215] = 1.42784349E-02;
    COFD[2216] = -2.13084334E+01;
    COFD[2217] = 5.27210469E+00;
    COFD[2218] = -4.21419216E-01;
    COFD[2219] = 1.63567178E-02;
    COFD[2220] = -2.15126310E+01;
    COFD[2221] = 5.48426911E+00;
    COFD[2222] = -4.80606512E-01;
    COFD[2223] = 2.01811046E-02;
    COFD[2224] = -2.15225573E+01;
    COFD[2225] = 5.48426911E+00;
    COFD[2226] = -4.80606512E-01;
    COFD[2227] = 2.01811046E-02;
    COFD[2228] = -1.82872310E+01;
    COFD[2229] = 4.40289649E+00;
    COFD[2230] = -3.57289765E-01;
    COFD[2231] = 1.55166804E-02;
    COFD[2232] = -2.11931178E+01;
    COFD[2233] = 5.40060531E+00;
    COFD[2234] = -4.72449699E-01;
    COFD[2235] = 1.99345817E-02;
    COFD[2236] = -2.19615570E+01;
    COFD[2237] = 5.43750833E+00;
    COFD[2238] = -4.50273329E-01;
    COFD[2239] = 1.79013718E-02;
    COFD[2240] = -2.14737305E+01;
    COFD[2241] = 5.41585806E+00;
    COFD[2242] = -4.73359323E-01;
    COFD[2243] = 1.99310239E-02;
    COFD[2244] = -2.19700270E+01;
    COFD[2245] = 5.43750833E+00;
    COFD[2246] = -4.50273329E-01;
    COFD[2247] = 1.79013718E-02;
    COFD[2248] = -2.14821817E+01;
    COFD[2249] = 5.41585806E+00;
    COFD[2250] = -4.73359323E-01;
    COFD[2251] = 1.99310239E-02;
    COFD[2252] = -2.21384805E+01;
    COFD[2253] = 5.56656297E+00;
    COFD[2254] = -4.75500048E-01;
    COFD[2255] = 1.93332291E-02;
    COFD[2256] = -2.21384805E+01;
    COFD[2257] = 5.56656297E+00;
    COFD[2258] = -4.75500048E-01;
    COFD[2259] = 1.93332291E-02;
    COFD[2260] = -1.85756076E+01;
    COFD[2261] = 4.51052425E+00;
    COFD[2262] = -3.70301627E-01;
    COFD[2263] = 1.60416153E-02;
    COFD[2264] = -1.85829283E+01;
    COFD[2265] = 4.51052425E+00;
    COFD[2266] = -3.70301627E-01;
    COFD[2267] = 1.60416153E-02;
    COFD[2268] = -1.85899144E+01;
    COFD[2269] = 4.51052425E+00;
    COFD[2270] = -3.70301627E-01;
    COFD[2271] = 1.60416153E-02;
    COFD[2272] = -2.14049543E+01;
    COFD[2273] = 5.41122754E+00;
    COFD[2274] = -4.73185889E-01;
    COFD[2275] = 1.99407905E-02;
    COFD[2276] = -2.22424134E+01;
    COFD[2277] = 5.53139819E+00;
    COFD[2278] = -4.68828555E-01;
    COFD[2279] = 1.89597887E-02;
    COFD[2280] = -2.23333008E+01;
    COFD[2281] = 5.61213963E+00;
    COFD[2282] = -4.90953451E-01;
    COFD[2283] = 2.03841518E-02;
    COFD[2284] = -2.23342576E+01;
    COFD[2285] = 5.49239750E+00;
    COFD[2286] = -4.60320987E-01;
    COFD[2287] = 1.84538922E-02;
    COFD[2288] = -2.23376085E+01;
    COFD[2289] = 5.61213963E+00;
    COFD[2290] = -4.90953451E-01;
    COFD[2291] = 2.03841518E-02;
    COFD[2292] = -2.21229141E+01;
    COFD[2293] = 5.47072190E+00;
    COFD[2294] = -4.56301261E-01;
    COFD[2295] = 1.82313566E-02;
    COFD[2296] = -2.25191362E+01;
    COFD[2297] = 5.58153976E+00;
    COFD[2298] = -4.78600113E-01;
    COFD[2299] = 1.95139120E-02;
    COFD[2300] = -2.25219121E+01;
    COFD[2301] = 5.58153976E+00;
    COFD[2302] = -4.78600113E-01;
    COFD[2303] = 1.95139120E-02;
    COFD[2304] = -2.24844780E+01;
    COFD[2305] = 5.49194178E+00;
    COFD[2306] = -4.60232233E-01;
    COFD[2307] = 1.84488745E-02;
    COFD[2308] = -2.23901271E+01;
    COFD[2309] = 5.61213963E+00;
    COFD[2310] = -4.90953451E-01;
    COFD[2311] = 2.03841518E-02;
    COFD[2312] = -2.24184864E+01;
    COFD[2313] = 5.61213963E+00;
    COFD[2314] = -4.90953451E-01;
    COFD[2315] = 2.03841518E-02;
    COFD[2316] = -2.24184864E+01;
    COFD[2317] = 5.61213963E+00;
    COFD[2318] = -4.90953451E-01;
    COFD[2319] = 2.03841518E-02;
    COFD[2320] = -2.24221512E+01;
    COFD[2321] = 5.61213963E+00;
    COFD[2322] = -4.90953451E-01;
    COFD[2323] = 2.03841518E-02;
    COFD[2324] = -2.24221512E+01;
    COFD[2325] = 5.61213963E+00;
    COFD[2326] = -4.90953451E-01;
    COFD[2327] = 2.03841518E-02;
    COFD[2328] = -2.24449636E+01;
    COFD[2329] = 5.61213963E+00;
    COFD[2330] = -4.90953451E-01;
    COFD[2331] = 2.03841518E-02;
    COFD[2332] = -2.24640329E+01;
    COFD[2333] = 5.61213963E+00;
    COFD[2334] = -4.90953451E-01;
    COFD[2335] = 2.03841518E-02;
    COFD[2336] = -1.82590824E+01;
    COFD[2337] = 4.39538102E+00;
    COFD[2338] = -3.56367230E-01;
    COFD[2339] = 1.54788461E-02;
    COFD[2340] = -1.97550088E+01;
    COFD[2341] = 5.56931926E+00;
    COFD[2342] = -4.89105511E-01;
    COFD[2343] = 2.04493129E-02;
    COFD[2344] = -1.60528285E+01;
    COFD[2345] = 4.11188603E+00;
    COFD[2346] = -3.21540884E-01;
    COFD[2347] = 1.40482564E-02;
    COFD[2348] = -2.14215700E+01;
    COFD[2349] = 5.56531152E+00;
    COFD[2350] = -4.88789821E-01;
    COFD[2351] = 2.04437116E-02;
    COFD[2352] = -2.14215700E+01;
    COFD[2353] = 5.56531152E+00;
    COFD[2354] = -4.88789821E-01;
    COFD[2355] = 2.04437116E-02;
    COFD[2356] = -2.14449559E+01;
    COFD[2357] = 5.56531152E+00;
    COFD[2358] = -4.88789821E-01;
    COFD[2359] = 2.04437116E-02;
    COFD[2360] = -1.94373127E+01;
    COFD[2361] = 5.02567894E+00;
    COFD[2362] = -4.32045169E-01;
    COFD[2363] = 1.85132214E-02;
    COFD[2364] = -2.14082453E+01;
    COFD[2365] = 5.55346617E+00;
    COFD[2366] = -4.87783156E-01;
    COFD[2367] = 2.04210886E-02;
    COFD[2368] = -1.94570287E+01;
    COFD[2369] = 5.02567894E+00;
    COFD[2370] = -4.32045169E-01;
    COFD[2371] = 1.85132214E-02;
    COFD[2372] = -1.77563250E+01;
    COFD[2373] = 3.57475686E+00;
    COFD[2374] = -1.56396297E-01;
    COFD[2375] = 3.12157721E-03;
    COFD[2376] = -2.16798265E+01;
    COFD[2377] = 5.36811769E+00;
    COFD[2378] = -4.37727086E-01;
    COFD[2379] = 1.72167686E-02;
    COFD[2380] = -2.16899073E+01;
    COFD[2381] = 5.36811769E+00;
    COFD[2382] = -4.37727086E-01;
    COFD[2383] = 1.72167686E-02;
    COFD[2384] = -2.05128705E+01;
    COFD[2385] = 5.23843909E+00;
    COFD[2386] = -4.55815614E-01;
    COFD[2387] = 1.93898040E-02;
    COFD[2388] = -2.17855148E+01;
    COFD[2389] = 5.47519298E+00;
    COFD[2390] = -4.57113040E-01;
    COFD[2391] = 1.82758312E-02;
    COFD[2392] = -1.90413348E+01;
    COFD[2393] = 3.99221757E+00;
    COFD[2394] = -2.19854880E-01;
    COFD[2395] = 6.22736279E-03;
    COFD[2396] = -2.19700270E+01;
    COFD[2397] = 5.43750833E+00;
    COFD[2398] = -4.50273329E-01;
    COFD[2399] = 1.79013718E-02;
    COFD[2400] = -1.90499441E+01;
    COFD[2401] = 3.99221757E+00;
    COFD[2402] = -2.19854880E-01;
    COFD[2403] = 6.22736279E-03;
    COFD[2404] = -2.19786173E+01;
    COFD[2405] = 5.43750833E+00;
    COFD[2406] = -4.50273329E-01;
    COFD[2407] = 1.79013718E-02;
    COFD[2408] = -2.01889168E+01;
    COFD[2409] = 4.53183330E+00;
    COFD[2410] = -3.02186760E-01;
    COFD[2411] = 1.02756490E-02;
    COFD[2412] = -2.01889168E+01;
    COFD[2413] = 4.53183330E+00;
    COFD[2414] = -3.02186760E-01;
    COFD[2415] = 1.02756490E-02;
    COFD[2416] = -2.08293255E+01;
    COFD[2417] = 5.35267674E+00;
    COFD[2418] = -4.69010505E-01;
    COFD[2419] = 1.98979152E-02;
    COFD[2420] = -2.08367725E+01;
    COFD[2421] = 5.35267674E+00;
    COFD[2422] = -4.69010505E-01;
    COFD[2423] = 1.98979152E-02;
    COFD[2424] = -2.08438809E+01;
    COFD[2425] = 5.35267674E+00;
    COFD[2426] = -4.69010505E-01;
    COFD[2427] = 1.98979152E-02;
    COFD[2428] = -2.19317743E+01;
    COFD[2429] = 5.45216133E+00;
    COFD[2430] = -4.52916925E-01;
    COFD[2431] = 1.80456400E-02;
    COFD[2432] = -2.01111595E+01;
    COFD[2433] = 4.41511629E+00;
    COFD[2434] = -2.84086963E-01;
    COFD[2435] = 9.37586971E-03;
    COFD[2436] = -2.15057773E+01;
    COFD[2437] = 5.05104142E+00;
    COFD[2438] = -3.84428177E-01;
    COFD[2439] = 1.44249285E-02;
    COFD[2440] = -1.97691133E+01;
    COFD[2441] = 4.18758010E+00;
    COFD[2442] = -2.49327776E-01;
    COFD[2443] = 7.66559103E-03;
    COFD[2444] = -2.15101711E+01;
    COFD[2445] = 5.05104142E+00;
    COFD[2446] = -3.84428177E-01;
    COFD[2447] = 1.44249285E-02;
    COFD[2448] = -1.94051843E+01;
    COFD[2449] = 4.10954793E+00;
    COFD[2450] = -2.37523329E-01;
    COFD[2451] = 7.08858141E-03;
    COFD[2452] = -2.07014385E+01;
    COFD[2453] = 4.59219167E+00;
    COFD[2454] = -3.11627566E-01;
    COFD[2455] = 1.07472587E-02;
    COFD[2456] = -2.07042756E+01;
    COFD[2457] = 4.59219167E+00;
    COFD[2458] = -3.11627566E-01;
    COFD[2459] = 1.07472587E-02;
    COFD[2460] = -1.99255944E+01;
    COFD[2461] = 4.18537863E+00;
    COFD[2462] = -2.48994776E-01;
    COFD[2463] = 7.64930102E-03;
    COFD[2464] = -2.15638034E+01;
    COFD[2465] = 5.05104142E+00;
    COFD[2466] = -3.84428177E-01;
    COFD[2467] = 1.44249285E-02;
    COFD[2468] = -2.15928156E+01;
    COFD[2469] = 5.05104142E+00;
    COFD[2470] = -3.84428177E-01;
    COFD[2471] = 1.44249285E-02;
    COFD[2472] = -2.15928156E+01;
    COFD[2473] = 5.05104142E+00;
    COFD[2474] = -3.84428177E-01;
    COFD[2475] = 1.44249285E-02;
    COFD[2476] = -2.15965677E+01;
    COFD[2477] = 5.05104142E+00;
    COFD[2478] = -3.84428177E-01;
    COFD[2479] = 1.44249285E-02;
    COFD[2480] = -2.15965677E+01;
    COFD[2481] = 5.05104142E+00;
    COFD[2482] = -3.84428177E-01;
    COFD[2483] = 1.44249285E-02;
    COFD[2484] = -2.16199377E+01;
    COFD[2485] = 5.05104142E+00;
    COFD[2486] = -3.84428177E-01;
    COFD[2487] = 1.44249285E-02;
    COFD[2488] = -2.16394936E+01;
    COFD[2489] = 5.05104142E+00;
    COFD[2490] = -3.84428177E-01;
    COFD[2491] = 1.44249285E-02;
    COFD[2492] = -2.04833713E+01;
    COFD[2493] = 5.23112374E+00;
    COFD[2494] = -4.54967682E-01;
    COFD[2495] = 1.93570423E-02;
    COFD[2496] = -1.78637178E+01;
    COFD[2497] = 4.88268692E+00;
    COFD[2498] = -4.14917638E-01;
    COFD[2499] = 1.78274298E-02;
    COFD[2500] = -1.39658996E+01;
    COFD[2501] = 3.24966086E+00;
    COFD[2502] = -2.11199992E-01;
    COFD[2503] = 9.32580661E-03;
    COFD[2504] = -1.94585111E+01;
    COFD[2505] = 4.87180830E+00;
    COFD[2506] = -4.13582958E-01;
    COFD[2507] = 1.77726094E-02;
    COFD[2508] = -1.94585111E+01;
    COFD[2509] = 4.87180830E+00;
    COFD[2510] = -4.13582958E-01;
    COFD[2511] = 1.77726094E-02;
    COFD[2512] = -1.94819080E+01;
    COFD[2513] = 4.87180830E+00;
    COFD[2514] = -4.13582958E-01;
    COFD[2515] = 1.77726094E-02;
    COFD[2516] = -1.72112971E+01;
    COFD[2517] = 4.15807461E+00;
    COFD[2518] = -3.27178539E-01;
    COFD[2519] = 1.42784349E-02;
    COFD[2520] = -1.94186547E+01;
    COFD[2521] = 4.84669430E+00;
    COFD[2522] = -4.10571455E-01;
    COFD[2523] = 1.76520543E-02;
    COFD[2524] = -1.72310232E+01;
    COFD[2525] = 4.15807461E+00;
    COFD[2526] = -3.27178539E-01;
    COFD[2527] = 1.42784349E-02;
    COFD[2528] = -2.13148887E+01;
    COFD[2529] = 5.27210469E+00;
    COFD[2530] = -4.21419216E-01;
    COFD[2531] = 1.63567178E-02;
    COFD[2532] = -2.15206146E+01;
    COFD[2533] = 5.48426911E+00;
    COFD[2534] = -4.80606512E-01;
    COFD[2535] = 2.01811046E-02;
    COFD[2536] = -2.15307023E+01;
    COFD[2537] = 5.48426911E+00;
    COFD[2538] = -4.80606512E-01;
    COFD[2539] = 2.01811046E-02;
    COFD[2540] = -1.82955252E+01;
    COFD[2541] = 4.40289649E+00;
    COFD[2542] = -3.57289765E-01;
    COFD[2543] = 1.55166804E-02;
    COFD[2544] = -2.12014186E+01;
    COFD[2545] = 5.40060531E+00;
    COFD[2546] = -4.72449699E-01;
    COFD[2547] = 1.99345817E-02;
    COFD[2548] = -2.19700018E+01;
    COFD[2549] = 5.43750833E+00;
    COFD[2550] = -4.50273329E-01;
    COFD[2551] = 1.79013718E-02;
    COFD[2552] = -2.14821817E+01;
    COFD[2553] = 5.41585806E+00;
    COFD[2554] = -4.73359323E-01;
    COFD[2555] = 1.99310239E-02;
    COFD[2556] = -2.19786173E+01;
    COFD[2557] = 5.43750833E+00;
    COFD[2558] = -4.50273329E-01;
    COFD[2559] = 1.79013718E-02;
    COFD[2560] = -2.14907782E+01;
    COFD[2561] = 5.41585806E+00;
    COFD[2562] = -4.73359323E-01;
    COFD[2563] = 1.99310239E-02;
    COFD[2564] = -2.21472114E+01;
    COFD[2565] = 5.56656297E+00;
    COFD[2566] = -4.75500048E-01;
    COFD[2567] = 1.93332291E-02;
    COFD[2568] = -2.21472114E+01;
    COFD[2569] = 5.56656297E+00;
    COFD[2570] = -4.75500048E-01;
    COFD[2571] = 1.93332291E-02;
    COFD[2572] = -1.85844688E+01;
    COFD[2573] = 4.51052425E+00;
    COFD[2574] = -3.70301627E-01;
    COFD[2575] = 1.60416153E-02;
    COFD[2576] = -1.85919214E+01;
    COFD[2577] = 4.51052425E+00;
    COFD[2578] = -3.70301627E-01;
    COFD[2579] = 1.60416153E-02;
    COFD[2580] = -1.85990352E+01;
    COFD[2581] = 4.51052425E+00;
    COFD[2582] = -3.70301627E-01;
    COFD[2583] = 1.60416153E-02;
    COFD[2584] = -2.14151520E+01;
    COFD[2585] = 5.41122754E+00;
    COFD[2586] = -4.73185889E-01;
    COFD[2587] = 1.99407905E-02;
    COFD[2588] = -2.22526152E+01;
    COFD[2589] = 5.53139819E+00;
    COFD[2590] = -4.68828555E-01;
    COFD[2591] = 1.89597887E-02;
    COFD[2592] = -2.23435951E+01;
    COFD[2593] = 5.61213963E+00;
    COFD[2594] = -4.90953451E-01;
    COFD[2595] = 2.03841518E-02;
    COFD[2596] = -2.23446381E+01;
    COFD[2597] = 5.49239750E+00;
    COFD[2598] = -4.60320987E-01;
    COFD[2599] = 1.84538922E-02;
    COFD[2600] = -2.23479928E+01;
    COFD[2601] = 5.61213963E+00;
    COFD[2602] = -4.90953451E-01;
    COFD[2603] = 2.03841518E-02;
    COFD[2604] = -2.21333822E+01;
    COFD[2605] = 5.47072190E+00;
    COFD[2606] = -4.56301261E-01;
    COFD[2607] = 1.82313566E-02;
    COFD[2608] = -2.25304962E+01;
    COFD[2609] = 5.58153976E+00;
    COFD[2610] = -4.78600113E-01;
    COFD[2611] = 1.95139120E-02;
    COFD[2612] = -2.25333361E+01;
    COFD[2613] = 5.58153976E+00;
    COFD[2614] = -4.78600113E-01;
    COFD[2615] = 1.95139120E-02;
    COFD[2616] = -2.24959645E+01;
    COFD[2617] = 5.49194178E+00;
    COFD[2618] = -4.60232233E-01;
    COFD[2619] = 1.84488745E-02;
    COFD[2620] = -2.24016748E+01;
    COFD[2621] = 5.61213963E+00;
    COFD[2622] = -4.90953451E-01;
    COFD[2623] = 2.03841518E-02;
    COFD[2624] = -2.24307163E+01;
    COFD[2625] = 5.61213963E+00;
    COFD[2626] = -4.90953451E-01;
    COFD[2627] = 2.03841518E-02;
    COFD[2628] = -2.24307163E+01;
    COFD[2629] = 5.61213963E+00;
    COFD[2630] = -4.90953451E-01;
    COFD[2631] = 2.03841518E-02;
    COFD[2632] = -2.24344722E+01;
    COFD[2633] = 5.61213963E+00;
    COFD[2634] = -4.90953451E-01;
    COFD[2635] = 2.03841518E-02;
    COFD[2636] = -2.24344722E+01;
    COFD[2637] = 5.61213963E+00;
    COFD[2638] = -4.90953451E-01;
    COFD[2639] = 2.03841518E-02;
    COFD[2640] = -2.24578673E+01;
    COFD[2641] = 5.61213963E+00;
    COFD[2642] = -4.90953451E-01;
    COFD[2643] = 2.03841518E-02;
    COFD[2644] = -2.24774450E+01;
    COFD[2645] = 5.61213963E+00;
    COFD[2646] = -4.90953451E-01;
    COFD[2647] = 2.03841518E-02;
    COFD[2648] = -1.82673770E+01;
    COFD[2649] = 4.39538102E+00;
    COFD[2650] = -3.56367230E-01;
    COFD[2651] = 1.54788461E-02;
    COFD[2652] = -1.92718582E+01;
    COFD[2653] = 5.41172124E+00;
    COFD[2654] = -4.73213887E-01;
    COFD[2655] = 1.99405473E-02;
    COFD[2656] = -1.58456300E+01;
    COFD[2657] = 4.02074783E+00;
    COFD[2658] = -3.10018522E-01;
    COFD[2659] = 1.35599552E-02;
    COFD[2660] = -2.09376196E+01;
    COFD[2661] = 5.40870099E+00;
    COFD[2662] = -4.73017610E-01;
    COFD[2663] = 1.99399066E-02;
    COFD[2664] = -2.09376196E+01;
    COFD[2665] = 5.40870099E+00;
    COFD[2666] = -4.73017610E-01;
    COFD[2667] = 1.99399066E-02;
    COFD[2668] = -2.09612557E+01;
    COFD[2669] = 5.40870099E+00;
    COFD[2670] = -4.73017610E-01;
    COFD[2671] = 1.99399066E-02;
    COFD[2672] = -1.88179418E+01;
    COFD[2673] = 4.79683898E+00;
    COFD[2674] = -4.04829719E-01;
    COFD[2675] = 1.74325475E-02;
    COFD[2676] = -2.11381508E+01;
    COFD[2677] = 5.45574440E+00;
    COFD[2678] = -4.77436155E-01;
    COFD[2679] = 2.00644596E-02;
    COFD[2680] = -1.88378874E+01;
    COFD[2681] = 4.79683898E+00;
    COFD[2682] = -4.04829719E-01;
    COFD[2683] = 1.74325475E-02;
    COFD[2684] = -1.65295288E+01;
    COFD[2685] = 2.97569206E+00;
    COFD[2686] = -6.75652842E-02;
    COFD[2687] = -1.08648422E-03;
    COFD[2688] = -2.18848136E+01;
    COFD[2689] = 5.51302074E+00;
    COFD[2690] = -4.65263979E-01;
    COFD[2691] = 1.87580679E-02;
    COFD[2692] = -2.18950505E+01;
    COFD[2693] = 5.51302074E+00;
    COFD[2694] = -4.65263979E-01;
    COFD[2695] = 1.87580679E-02;
    COFD[2696] = -2.02642227E+01;
    COFD[2697] = 5.14499740E+00;
    COFD[2698] = -4.45694430E-01;
    COFD[2699] = 1.90318646E-02;
    COFD[2700] = -2.19136842E+01;
    COFD[2701] = 5.58503445E+00;
    COFD[2702] = -4.79552117E-01;
    COFD[2703] = 1.95750393E-02;
    COFD[2704] = -2.01801667E+01;
    COFD[2705] = 4.53183330E+00;
    COFD[2706] = -3.02186760E-01;
    COFD[2707] = 1.02756490E-02;
    COFD[2708] = -2.21384805E+01;
    COFD[2709] = 5.56656297E+00;
    COFD[2710] = -4.75500048E-01;
    COFD[2711] = 1.93332291E-02;
    COFD[2712] = -2.01889168E+01;
    COFD[2713] = 4.53183330E+00;
    COFD[2714] = -3.02186760E-01;
    COFD[2715] = 1.02756490E-02;
    COFD[2716] = -2.21472114E+01;
    COFD[2717] = 5.56656297E+00;
    COFD[2718] = -4.75500048E-01;
    COFD[2719] = 1.93332291E-02;
    COFD[2720] = -1.95877017E+01;
    COFD[2721] = 4.27643051E+00;
    COFD[2722] = -2.68040901E-01;
    COFD[2723] = 8.77650113E-03;
    COFD[2724] = -1.95877017E+01;
    COFD[2725] = 4.27643051E+00;
    COFD[2726] = -2.68040901E-01;
    COFD[2727] = 8.77650113E-03;
    COFD[2728] = -2.04928958E+01;
    COFD[2729] = 5.22397933E+00;
    COFD[2730] = -4.54138171E-01;
    COFD[2731] = 1.93249285E-02;
    COFD[2732] = -2.02637994E+01;
    COFD[2733] = 5.14984081E+00;
    COFD[2734] = -4.46093018E-01;
    COFD[2735] = 1.90396647E-02;
    COFD[2736] = -2.02710316E+01;
    COFD[2737] = 5.14984081E+00;
    COFD[2738] = -4.46093018E-01;
    COFD[2739] = 1.90396647E-02;
    COFD[2740] = -2.20421041E+01;
    COFD[2741] = 5.52708332E+00;
    COFD[2742] = -4.68000808E-01;
    COFD[2743] = 1.89131908E-02;
    COFD[2744] = -2.09320587E+01;
    COFD[2745] = 4.82184721E+00;
    COFD[2746] = -3.48128875E-01;
    COFD[2747] = 1.25918978E-02;
    COFD[2748] = -2.19952270E+01;
    COFD[2749] = 5.31843386E+00;
    COFD[2750] = -4.29198961E-01;
    COFD[2751] = 1.67645078E-02;
    COFD[2752] = -2.07399014E+01;
    COFD[2753] = 4.65728078E+00;
    COFD[2754] = -3.22002062E-01;
    COFD[2755] = 1.12723316E-02;
    COFD[2756] = -2.19997085E+01;
    COFD[2757] = 5.31843386E+00;
    COFD[2758] = -4.29198961E-01;
    COFD[2759] = 1.67645078E-02;
    COFD[2760] = -2.04451935E+01;
    COFD[2761] = 4.60682543E+00;
    COFD[2762] = -3.13971634E-01;
    COFD[2763] = 1.08661011E-02;
    COFD[2764] = -2.14849242E+01;
    COFD[2765] = 4.98394310E+00;
    COFD[2766] = -3.73694768E-01;
    COFD[2767] = 1.38795565E-02;
    COFD[2768] = -2.14878240E+01;
    COFD[2769] = 4.98394310E+00;
    COFD[2770] = -3.73694768E-01;
    COFD[2771] = 1.38795565E-02;
    COFD[2772] = -2.08978461E+01;
    COFD[2773] = 4.65579398E+00;
    COFD[2774] = -3.21767942E-01;
    COFD[2775] = 1.12605667E-02;
    COFD[2776] = -2.20544777E+01;
    COFD[2777] = 5.31843386E+00;
    COFD[2778] = -4.29198961E-01;
    COFD[2779] = 1.67645078E-02;
    COFD[2780] = -2.20841588E+01;
    COFD[2781] = 5.31843386E+00;
    COFD[2782] = -4.29198961E-01;
    COFD[2783] = 1.67645078E-02;
    COFD[2784] = -2.20841588E+01;
    COFD[2785] = 5.31843386E+00;
    COFD[2786] = -4.29198961E-01;
    COFD[2787] = 1.67645078E-02;
    COFD[2788] = -2.20880003E+01;
    COFD[2789] = 5.31843386E+00;
    COFD[2790] = -4.29198961E-01;
    COFD[2791] = 1.67645078E-02;
    COFD[2792] = -2.20880003E+01;
    COFD[2793] = 5.31843386E+00;
    COFD[2794] = -4.29198961E-01;
    COFD[2795] = 1.67645078E-02;
    COFD[2796] = -2.21119432E+01;
    COFD[2797] = 5.31843386E+00;
    COFD[2798] = -4.29198961E-01;
    COFD[2799] = 1.67645078E-02;
    COFD[2800] = -2.21320000E+01;
    COFD[2801] = 5.31843386E+00;
    COFD[2802] = -4.29198961E-01;
    COFD[2803] = 1.67645078E-02;
    COFD[2804] = -2.02268902E+01;
    COFD[2805] = 5.13632093E+00;
    COFD[2806] = -4.44839124E-01;
    COFD[2807] = 1.90058354E-02;
    COFD[2808] = -1.92718582E+01;
    COFD[2809] = 5.41172124E+00;
    COFD[2810] = -4.73213887E-01;
    COFD[2811] = 1.99405473E-02;
    COFD[2812] = -1.58456300E+01;
    COFD[2813] = 4.02074783E+00;
    COFD[2814] = -3.10018522E-01;
    COFD[2815] = 1.35599552E-02;
    COFD[2816] = -2.09376196E+01;
    COFD[2817] = 5.40870099E+00;
    COFD[2818] = -4.73017610E-01;
    COFD[2819] = 1.99399066E-02;
    COFD[2820] = -2.09376196E+01;
    COFD[2821] = 5.40870099E+00;
    COFD[2822] = -4.73017610E-01;
    COFD[2823] = 1.99399066E-02;
    COFD[2824] = -2.09612557E+01;
    COFD[2825] = 5.40870099E+00;
    COFD[2826] = -4.73017610E-01;
    COFD[2827] = 1.99399066E-02;
    COFD[2828] = -1.88179418E+01;
    COFD[2829] = 4.79683898E+00;
    COFD[2830] = -4.04829719E-01;
    COFD[2831] = 1.74325475E-02;
    COFD[2832] = -2.11381508E+01;
    COFD[2833] = 5.45574440E+00;
    COFD[2834] = -4.77436155E-01;
    COFD[2835] = 2.00644596E-02;
    COFD[2836] = -1.88378874E+01;
    COFD[2837] = 4.79683898E+00;
    COFD[2838] = -4.04829719E-01;
    COFD[2839] = 1.74325475E-02;
    COFD[2840] = -1.65295288E+01;
    COFD[2841] = 2.97569206E+00;
    COFD[2842] = -6.75652842E-02;
    COFD[2843] = -1.08648422E-03;
    COFD[2844] = -2.18848136E+01;
    COFD[2845] = 5.51302074E+00;
    COFD[2846] = -4.65263979E-01;
    COFD[2847] = 1.87580679E-02;
    COFD[2848] = -2.18950505E+01;
    COFD[2849] = 5.51302074E+00;
    COFD[2850] = -4.65263979E-01;
    COFD[2851] = 1.87580679E-02;
    COFD[2852] = -2.02642227E+01;
    COFD[2853] = 5.14499740E+00;
    COFD[2854] = -4.45694430E-01;
    COFD[2855] = 1.90318646E-02;
    COFD[2856] = -2.19136842E+01;
    COFD[2857] = 5.58503445E+00;
    COFD[2858] = -4.79552117E-01;
    COFD[2859] = 1.95750393E-02;
    COFD[2860] = -2.01801667E+01;
    COFD[2861] = 4.53183330E+00;
    COFD[2862] = -3.02186760E-01;
    COFD[2863] = 1.02756490E-02;
    COFD[2864] = -2.21384805E+01;
    COFD[2865] = 5.56656297E+00;
    COFD[2866] = -4.75500048E-01;
    COFD[2867] = 1.93332291E-02;
    COFD[2868] = -2.01889168E+01;
    COFD[2869] = 4.53183330E+00;
    COFD[2870] = -3.02186760E-01;
    COFD[2871] = 1.02756490E-02;
    COFD[2872] = -2.21472114E+01;
    COFD[2873] = 5.56656297E+00;
    COFD[2874] = -4.75500048E-01;
    COFD[2875] = 1.93332291E-02;
    COFD[2876] = -1.95877017E+01;
    COFD[2877] = 4.27643051E+00;
    COFD[2878] = -2.68040901E-01;
    COFD[2879] = 8.77650113E-03;
    COFD[2880] = -1.95877017E+01;
    COFD[2881] = 4.27643051E+00;
    COFD[2882] = -2.68040901E-01;
    COFD[2883] = 8.77650113E-03;
    COFD[2884] = -2.04928958E+01;
    COFD[2885] = 5.22397933E+00;
    COFD[2886] = -4.54138171E-01;
    COFD[2887] = 1.93249285E-02;
    COFD[2888] = -2.02637994E+01;
    COFD[2889] = 5.14984081E+00;
    COFD[2890] = -4.46093018E-01;
    COFD[2891] = 1.90396647E-02;
    COFD[2892] = -2.02710316E+01;
    COFD[2893] = 5.14984081E+00;
    COFD[2894] = -4.46093018E-01;
    COFD[2895] = 1.90396647E-02;
    COFD[2896] = -2.20421041E+01;
    COFD[2897] = 5.52708332E+00;
    COFD[2898] = -4.68000808E-01;
    COFD[2899] = 1.89131908E-02;
    COFD[2900] = -2.09320587E+01;
    COFD[2901] = 4.82184721E+00;
    COFD[2902] = -3.48128875E-01;
    COFD[2903] = 1.25918978E-02;
    COFD[2904] = -2.19952270E+01;
    COFD[2905] = 5.31843386E+00;
    COFD[2906] = -4.29198961E-01;
    COFD[2907] = 1.67645078E-02;
    COFD[2908] = -2.07399014E+01;
    COFD[2909] = 4.65728078E+00;
    COFD[2910] = -3.22002062E-01;
    COFD[2911] = 1.12723316E-02;
    COFD[2912] = -2.19997085E+01;
    COFD[2913] = 5.31843386E+00;
    COFD[2914] = -4.29198961E-01;
    COFD[2915] = 1.67645078E-02;
    COFD[2916] = -2.04451935E+01;
    COFD[2917] = 4.60682543E+00;
    COFD[2918] = -3.13971634E-01;
    COFD[2919] = 1.08661011E-02;
    COFD[2920] = -2.14849242E+01;
    COFD[2921] = 4.98394310E+00;
    COFD[2922] = -3.73694768E-01;
    COFD[2923] = 1.38795565E-02;
    COFD[2924] = -2.14878240E+01;
    COFD[2925] = 4.98394310E+00;
    COFD[2926] = -3.73694768E-01;
    COFD[2927] = 1.38795565E-02;
    COFD[2928] = -2.08978461E+01;
    COFD[2929] = 4.65579398E+00;
    COFD[2930] = -3.21767942E-01;
    COFD[2931] = 1.12605667E-02;
    COFD[2932] = -2.20544777E+01;
    COFD[2933] = 5.31843386E+00;
    COFD[2934] = -4.29198961E-01;
    COFD[2935] = 1.67645078E-02;
    COFD[2936] = -2.20841588E+01;
    COFD[2937] = 5.31843386E+00;
    COFD[2938] = -4.29198961E-01;
    COFD[2939] = 1.67645078E-02;
    COFD[2940] = -2.20841588E+01;
    COFD[2941] = 5.31843386E+00;
    COFD[2942] = -4.29198961E-01;
    COFD[2943] = 1.67645078E-02;
    COFD[2944] = -2.20880003E+01;
    COFD[2945] = 5.31843386E+00;
    COFD[2946] = -4.29198961E-01;
    COFD[2947] = 1.67645078E-02;
    COFD[2948] = -2.20880003E+01;
    COFD[2949] = 5.31843386E+00;
    COFD[2950] = -4.29198961E-01;
    COFD[2951] = 1.67645078E-02;
    COFD[2952] = -2.21119432E+01;
    COFD[2953] = 5.31843386E+00;
    COFD[2954] = -4.29198961E-01;
    COFD[2955] = 1.67645078E-02;
    COFD[2956] = -2.21320000E+01;
    COFD[2957] = 5.31843386E+00;
    COFD[2958] = -4.29198961E-01;
    COFD[2959] = 1.67645078E-02;
    COFD[2960] = -2.02268902E+01;
    COFD[2961] = 5.13632093E+00;
    COFD[2962] = -4.44839124E-01;
    COFD[2963] = 1.90058354E-02;
    COFD[2964] = -1.46550083E+01;
    COFD[2965] = 3.83606243E+00;
    COFD[2966] = -2.86076532E-01;
    COFD[2967] = 1.25205829E-02;
    COFD[2968] = -1.18988955E+01;
    COFD[2969] = 2.57507000E+00;
    COFD[2970] = -1.24033737E-01;
    COFD[2971] = 5.56694959E-03;
    COFD[2972] = -1.63254691E+01;
    COFD[2973] = 3.82388595E+00;
    COFD[2974] = -2.84480724E-01;
    COFD[2975] = 1.24506311E-02;
    COFD[2976] = -1.63254691E+01;
    COFD[2977] = 3.82388595E+00;
    COFD[2978] = -2.84480724E-01;
    COFD[2979] = 1.24506311E-02;
    COFD[2980] = -1.63493345E+01;
    COFD[2981] = 3.82388595E+00;
    COFD[2982] = -2.84480724E-01;
    COFD[2983] = 1.24506311E-02;
    COFD[2984] = -1.43139231E+01;
    COFD[2985] = 3.17651319E+00;
    COFD[2986] = -2.02028974E-01;
    COFD[2987] = 8.94232502E-03;
    COFD[2988] = -1.62724462E+01;
    COFD[2989] = 3.79163564E+00;
    COFD[2990] = -2.80257365E-01;
    COFD[2991] = 1.22656902E-02;
    COFD[2992] = -1.43340796E+01;
    COFD[2993] = 3.17651319E+00;
    COFD[2994] = -2.02028974E-01;
    COFD[2995] = 8.94232502E-03;
    COFD[2996] = -2.12652533E+01;
    COFD[2997] = 5.59961818E+00;
    COFD[2998] = -4.91624858E-01;
    COFD[2999] = 2.05035550E-02;
    COFD[3000] = -1.86507213E+01;
    COFD[3001] = 4.60874797E+00;
    COFD[3002] = -3.82368716E-01;
    COFD[3003] = 1.65370164E-02;
    COFD[3004] = -1.86611023E+01;
    COFD[3005] = 4.60874797E+00;
    COFD[3006] = -3.82368716E-01;
    COFD[3007] = 1.65370164E-02;
    COFD[3008] = -1.52721107E+01;
    COFD[3009] = 3.36790500E+00;
    COFD[3010] = -2.26321740E-01;
    COFD[3011] = 9.97135055E-03;
    COFD[3012] = -1.82145353E+01;
    COFD[3013] = 4.46848269E+00;
    COFD[3014] = -3.65269718E-01;
    COFD[3015] = 1.58407652E-02;
    COFD[3016] = -2.08204449E+01;
    COFD[3017] = 5.35267674E+00;
    COFD[3018] = -4.69010505E-01;
    COFD[3019] = 1.98979152E-02;
    COFD[3020] = -1.85756076E+01;
    COFD[3021] = 4.51052425E+00;
    COFD[3022] = -3.70301627E-01;
    COFD[3023] = 1.60416153E-02;
    COFD[3024] = -2.08293255E+01;
    COFD[3025] = 5.35267674E+00;
    COFD[3026] = -4.69010505E-01;
    COFD[3027] = 1.98979152E-02;
    COFD[3028] = -1.85844688E+01;
    COFD[3029] = 4.51052425E+00;
    COFD[3030] = -3.70301627E-01;
    COFD[3031] = 1.60416153E-02;
    COFD[3032] = -2.04928958E+01;
    COFD[3033] = 5.22397933E+00;
    COFD[3034] = -4.54138171E-01;
    COFD[3035] = 1.93249285E-02;
    COFD[3036] = -2.04928958E+01;
    COFD[3037] = 5.22397933E+00;
    COFD[3038] = -4.54138171E-01;
    COFD[3039] = 1.93249285E-02;
    COFD[3040] = -1.55511344E+01;
    COFD[3041] = 3.48070094E+00;
    COFD[3042] = -2.40859499E-01;
    COFD[3043] = 1.05972514E-02;
    COFD[3044] = -1.55588279E+01;
    COFD[3045] = 3.48070094E+00;
    COFD[3046] = -2.40859499E-01;
    COFD[3047] = 1.05972514E-02;
    COFD[3048] = -1.55661750E+01;
    COFD[3049] = 3.48070094E+00;
    COFD[3050] = -2.40859499E-01;
    COFD[3051] = 1.05972514E-02;
    COFD[3052] = -1.84688406E+01;
    COFD[3053] = 4.49330851E+00;
    COFD[3054] = -3.68208715E-01;
    COFD[3055] = 1.59565402E-02;
    COFD[3056] = -2.05284752E+01;
    COFD[3057] = 5.18417470E+00;
    COFD[3058] = -4.49491573E-01;
    COFD[3059] = 1.91438508E-02;
    COFD[3060] = -1.97454130E+01;
    COFD[3061] = 4.86038833E+00;
    COFD[3062] = -4.12187543E-01;
    COFD[3063] = 1.77155642E-02;
    COFD[3064] = -2.09524743E+01;
    COFD[3065] = 5.28755355E+00;
    COFD[3066] = -4.61641920E-01;
    COFD[3067] = 1.96208961E-02;
    COFD[3068] = -1.97499763E+01;
    COFD[3069] = 4.86038833E+00;
    COFD[3070] = -4.12187543E-01;
    COFD[3071] = 1.77155642E-02;
    COFD[3072] = -2.08463209E+01;
    COFD[3073] = 5.32244593E+00;
    COFD[3074] = -4.65829403E-01;
    COFD[3075] = 1.97895274E-02;
    COFD[3076] = -2.06003483E+01;
    COFD[3077] = 5.13220015E+00;
    COFD[3078] = -4.44412321E-01;
    COFD[3079] = 1.89917261E-02;
    COFD[3080] = -2.06033066E+01;
    COFD[3081] = 5.13220015E+00;
    COFD[3082] = -4.44412321E-01;
    COFD[3083] = 1.89917261E-02;
    COFD[3084] = -2.11202704E+01;
    COFD[3085] = 5.28856610E+00;
    COFD[3086] = -4.61764591E-01;
    COFD[3087] = 1.96258820E-02;
    COFD[3088] = -1.98058079E+01;
    COFD[3089] = 4.86038833E+00;
    COFD[3090] = -4.12187543E-01;
    COFD[3091] = 1.77155642E-02;
    COFD[3092] = -1.98361165E+01;
    COFD[3093] = 4.86038833E+00;
    COFD[3094] = -4.12187543E-01;
    COFD[3095] = 1.77155642E-02;
    COFD[3096] = -1.98361165E+01;
    COFD[3097] = 4.86038833E+00;
    COFD[3098] = -4.12187543E-01;
    COFD[3099] = 1.77155642E-02;
    COFD[3100] = -1.98400420E+01;
    COFD[3101] = 4.86038833E+00;
    COFD[3102] = -4.12187543E-01;
    COFD[3103] = 1.77155642E-02;
    COFD[3104] = -1.98400420E+01;
    COFD[3105] = 4.86038833E+00;
    COFD[3106] = -4.12187543E-01;
    COFD[3107] = 1.77155642E-02;
    COFD[3108] = -1.98645237E+01;
    COFD[3109] = 4.86038833E+00;
    COFD[3110] = -4.12187543E-01;
    COFD[3111] = 1.77155642E-02;
    COFD[3112] = -1.98850525E+01;
    COFD[3113] = 4.86038833E+00;
    COFD[3114] = -4.12187543E-01;
    COFD[3115] = 1.77155642E-02;
    COFD[3116] = -1.52414485E+01;
    COFD[3117] = 3.35922578E+00;
    COFD[3118] = -2.25181399E-01;
    COFD[3119] = 9.92132878E-03;
    COFD[3120] = -1.46554748E+01;
    COFD[3121] = 3.83606243E+00;
    COFD[3122] = -2.86076532E-01;
    COFD[3123] = 1.25205829E-02;
    COFD[3124] = -1.18998012E+01;
    COFD[3125] = 2.57507000E+00;
    COFD[3126] = -1.24033737E-01;
    COFD[3127] = 5.56694959E-03;
    COFD[3128] = -1.63301444E+01;
    COFD[3129] = 3.82388595E+00;
    COFD[3130] = -2.84480724E-01;
    COFD[3131] = 1.24506311E-02;
    COFD[3132] = -1.63301444E+01;
    COFD[3133] = 3.82388595E+00;
    COFD[3134] = -2.84480724E-01;
    COFD[3135] = 1.24506311E-02;
    COFD[3136] = -1.63542394E+01;
    COFD[3137] = 3.82388595E+00;
    COFD[3138] = -2.84480724E-01;
    COFD[3139] = 1.24506311E-02;
    COFD[3140] = -1.43190389E+01;
    COFD[3141] = 3.17651319E+00;
    COFD[3142] = -2.02028974E-01;
    COFD[3143] = 8.94232502E-03;
    COFD[3144] = -1.62775714E+01;
    COFD[3145] = 3.79163564E+00;
    COFD[3146] = -2.80257365E-01;
    COFD[3147] = 1.22656902E-02;
    COFD[3148] = -1.43394069E+01;
    COFD[3149] = 3.17651319E+00;
    COFD[3150] = -2.02028974E-01;
    COFD[3151] = 8.94232502E-03;
    COFD[3152] = -2.06463744E+01;
    COFD[3153] = 5.41688482E+00;
    COFD[3154] = -4.73387188E-01;
    COFD[3155] = 1.99280175E-02;
    COFD[3156] = -1.86576191E+01;
    COFD[3157] = 4.60874797E+00;
    COFD[3158] = -3.82368716E-01;
    COFD[3159] = 1.65370164E-02;
    COFD[3160] = -1.86681459E+01;
    COFD[3161] = 4.60874797E+00;
    COFD[3162] = -3.82368716E-01;
    COFD[3163] = 1.65370164E-02;
    COFD[3164] = -1.52792891E+01;
    COFD[3165] = 3.36790500E+00;
    COFD[3166] = -2.26321740E-01;
    COFD[3167] = 9.97135055E-03;
    COFD[3168] = -1.82217198E+01;
    COFD[3169] = 4.46848269E+00;
    COFD[3170] = -3.65269718E-01;
    COFD[3171] = 1.58407652E-02;
    COFD[3172] = -2.08277598E+01;
    COFD[3173] = 5.35267674E+00;
    COFD[3174] = -4.69010505E-01;
    COFD[3175] = 1.98979152E-02;
    COFD[3176] = -1.85829283E+01;
    COFD[3177] = 4.51052425E+00;
    COFD[3178] = -3.70301627E-01;
    COFD[3179] = 1.60416153E-02;
    COFD[3180] = -2.08367725E+01;
    COFD[3181] = 5.35267674E+00;
    COFD[3182] = -4.69010505E-01;
    COFD[3183] = 1.98979152E-02;
    COFD[3184] = -1.85919214E+01;
    COFD[3185] = 4.51052425E+00;
    COFD[3186] = -3.70301627E-01;
    COFD[3187] = 1.60416153E-02;
    COFD[3188] = -2.02637994E+01;
    COFD[3189] = 5.14984081E+00;
    COFD[3190] = -4.46093018E-01;
    COFD[3191] = 1.90396647E-02;
    COFD[3192] = -2.02637994E+01;
    COFD[3193] = 5.14984081E+00;
    COFD[3194] = -4.46093018E-01;
    COFD[3195] = 1.90396647E-02;
    COFD[3196] = -1.55588279E+01;
    COFD[3197] = 3.48070094E+00;
    COFD[3198] = -2.40859499E-01;
    COFD[3199] = 1.05972514E-02;
    COFD[3200] = -1.55666415E+01;
    COFD[3201] = 3.48070094E+00;
    COFD[3202] = -2.40859499E-01;
    COFD[3203] = 1.05972514E-02;
    COFD[3204] = -1.55741053E+01;
    COFD[3205] = 3.48070094E+00;
    COFD[3206] = -2.40859499E-01;
    COFD[3207] = 1.05972514E-02;
    COFD[3208] = -1.84777607E+01;
    COFD[3209] = 4.49330851E+00;
    COFD[3210] = -3.68208715E-01;
    COFD[3211] = 1.59565402E-02;
    COFD[3212] = -2.05373990E+01;
    COFD[3213] = 5.18417470E+00;
    COFD[3214] = -4.49491573E-01;
    COFD[3215] = 1.91438508E-02;
    COFD[3216] = -1.97544224E+01;
    COFD[3217] = 4.86038833E+00;
    COFD[3218] = -4.12187543E-01;
    COFD[3219] = 1.77155642E-02;
    COFD[3220] = -2.09615635E+01;
    COFD[3221] = 5.28755355E+00;
    COFD[3222] = -4.61641920E-01;
    COFD[3223] = 1.96208961E-02;
    COFD[3224] = -1.97590691E+01;
    COFD[3225] = 4.86038833E+00;
    COFD[3226] = -4.12187543E-01;
    COFD[3227] = 1.77155642E-02;
    COFD[3228] = -2.08554914E+01;
    COFD[3229] = 5.32244593E+00;
    COFD[3230] = -4.65829403E-01;
    COFD[3231] = 1.97895274E-02;
    COFD[3232] = -2.06103502E+01;
    COFD[3233] = 5.13220015E+00;
    COFD[3234] = -4.44412321E-01;
    COFD[3235] = 1.89917261E-02;
    COFD[3236] = -2.06133685E+01;
    COFD[3237] = 5.13220015E+00;
    COFD[3238] = -4.44412321E-01;
    COFD[3239] = 1.89917261E-02;
    COFD[3240] = -2.11303909E+01;
    COFD[3241] = 5.28856610E+00;
    COFD[3242] = -4.61764591E-01;
    COFD[3243] = 1.96258820E-02;
    COFD[3244] = -1.98159859E+01;
    COFD[3245] = 4.86038833E+00;
    COFD[3246] = -4.12187543E-01;
    COFD[3247] = 1.77155642E-02;
    COFD[3248] = -1.98469374E+01;
    COFD[3249] = 4.86038833E+00;
    COFD[3250] = -4.12187543E-01;
    COFD[3251] = 1.77155642E-02;
    COFD[3252] = -1.98469374E+01;
    COFD[3253] = 4.86038833E+00;
    COFD[3254] = -4.12187543E-01;
    COFD[3255] = 1.77155642E-02;
    COFD[3256] = -1.98509491E+01;
    COFD[3257] = 4.86038833E+00;
    COFD[3258] = -4.12187543E-01;
    COFD[3259] = 1.77155642E-02;
    COFD[3260] = -1.98509491E+01;
    COFD[3261] = 4.86038833E+00;
    COFD[3262] = -4.12187543E-01;
    COFD[3263] = 1.77155642E-02;
    COFD[3264] = -1.98759845E+01;
    COFD[3265] = 4.86038833E+00;
    COFD[3266] = -4.12187543E-01;
    COFD[3267] = 1.77155642E-02;
    COFD[3268] = -1.98969995E+01;
    COFD[3269] = 4.86038833E+00;
    COFD[3270] = -4.12187543E-01;
    COFD[3271] = 1.77155642E-02;
    COFD[3272] = -1.52486273E+01;
    COFD[3273] = 3.35922578E+00;
    COFD[3274] = -2.25181399E-01;
    COFD[3275] = 9.92132878E-03;
    COFD[3276] = -1.46559141E+01;
    COFD[3277] = 3.83606243E+00;
    COFD[3278] = -2.86076532E-01;
    COFD[3279] = 1.25205829E-02;
    COFD[3280] = -1.19006548E+01;
    COFD[3281] = 2.57507000E+00;
    COFD[3282] = -1.24033737E-01;
    COFD[3283] = 5.56694959E-03;
    COFD[3284] = -1.63345829E+01;
    COFD[3285] = 3.82388595E+00;
    COFD[3286] = -2.84480724E-01;
    COFD[3287] = 1.24506311E-02;
    COFD[3288] = -1.63345829E+01;
    COFD[3289] = 3.82388595E+00;
    COFD[3290] = -2.84480724E-01;
    COFD[3291] = 1.24506311E-02;
    COFD[3292] = -1.63588981E+01;
    COFD[3293] = 3.82388595E+00;
    COFD[3294] = -2.84480724E-01;
    COFD[3295] = 1.24506311E-02;
    COFD[3296] = -1.43238998E+01;
    COFD[3297] = 3.17651319E+00;
    COFD[3298] = -2.02028974E-01;
    COFD[3299] = 8.94232502E-03;
    COFD[3300] = -1.62824412E+01;
    COFD[3301] = 3.79163564E+00;
    COFD[3302] = -2.80257365E-01;
    COFD[3303] = 1.22656902E-02;
    COFD[3304] = -1.43444709E+01;
    COFD[3305] = 3.17651319E+00;
    COFD[3306] = -2.02028974E-01;
    COFD[3307] = 8.94232502E-03;
    COFD[3308] = -2.06516336E+01;
    COFD[3309] = 5.41688482E+00;
    COFD[3310] = -4.73387188E-01;
    COFD[3311] = 1.99280175E-02;
    COFD[3312] = -1.86641962E+01;
    COFD[3313] = 4.60874797E+00;
    COFD[3314] = -3.82368716E-01;
    COFD[3315] = 1.65370164E-02;
    COFD[3316] = -1.86748638E+01;
    COFD[3317] = 4.60874797E+00;
    COFD[3318] = -3.82368716E-01;
    COFD[3319] = 1.65370164E-02;
    COFD[3320] = -1.52861376E+01;
    COFD[3321] = 3.36790500E+00;
    COFD[3322] = -2.26321740E-01;
    COFD[3323] = 9.97135055E-03;
    COFD[3324] = -1.82285740E+01;
    COFD[3325] = 4.46848269E+00;
    COFD[3326] = -3.65269718E-01;
    COFD[3327] = 1.58407652E-02;
    COFD[3328] = -2.08347403E+01;
    COFD[3329] = 5.35267674E+00;
    COFD[3330] = -4.69010505E-01;
    COFD[3331] = 1.98979152E-02;
    COFD[3332] = -1.85899144E+01;
    COFD[3333] = 4.51052425E+00;
    COFD[3334] = -3.70301627E-01;
    COFD[3335] = 1.60416153E-02;
    COFD[3336] = -2.08438809E+01;
    COFD[3337] = 5.35267674E+00;
    COFD[3338] = -4.69010505E-01;
    COFD[3339] = 1.98979152E-02;
    COFD[3340] = -1.85990352E+01;
    COFD[3341] = 4.51052425E+00;
    COFD[3342] = -3.70301627E-01;
    COFD[3343] = 1.60416153E-02;
    COFD[3344] = -2.02710316E+01;
    COFD[3345] = 5.14984081E+00;
    COFD[3346] = -4.46093018E-01;
    COFD[3347] = 1.90396647E-02;
    COFD[3348] = -2.02710316E+01;
    COFD[3349] = 5.14984081E+00;
    COFD[3350] = -4.46093018E-01;
    COFD[3351] = 1.90396647E-02;
    COFD[3352] = -1.55661750E+01;
    COFD[3353] = 3.48070094E+00;
    COFD[3354] = -2.40859499E-01;
    COFD[3355] = 1.05972514E-02;
    COFD[3356] = -1.55741053E+01;
    COFD[3357] = 3.48070094E+00;
    COFD[3358] = -2.40859499E-01;
    COFD[3359] = 1.05972514E-02;
    COFD[3360] = -1.55816822E+01;
    COFD[3361] = 3.48070094E+00;
    COFD[3362] = -2.40859499E-01;
    COFD[3363] = 1.05972514E-02;
    COFD[3364] = -1.84863000E+01;
    COFD[3365] = 4.49330851E+00;
    COFD[3366] = -3.68208715E-01;
    COFD[3367] = 1.59565402E-02;
    COFD[3368] = -2.05459419E+01;
    COFD[3369] = 5.18417470E+00;
    COFD[3370] = -4.49491573E-01;
    COFD[3371] = 1.91438508E-02;
    COFD[3372] = -1.97630486E+01;
    COFD[3373] = 4.86038833E+00;
    COFD[3374] = -4.12187543E-01;
    COFD[3375] = 1.77155642E-02;
    COFD[3376] = -2.09702675E+01;
    COFD[3377] = 5.28755355E+00;
    COFD[3378] = -4.61641920E-01;
    COFD[3379] = 1.96208961E-02;
    COFD[3380] = -1.97677766E+01;
    COFD[3381] = 4.86038833E+00;
    COFD[3382] = -4.12187543E-01;
    COFD[3383] = 1.77155642E-02;
    COFD[3384] = -2.08642748E+01;
    COFD[3385] = 5.32244593E+00;
    COFD[3386] = -4.65829403E-01;
    COFD[3387] = 1.97895274E-02;
    COFD[3388] = -2.06199456E+01;
    COFD[3389] = 5.13220015E+00;
    COFD[3390] = -4.44412321E-01;
    COFD[3391] = 1.89917261E-02;
    COFD[3392] = -2.06230226E+01;
    COFD[3393] = 5.13220015E+00;
    COFD[3394] = -4.44412321E-01;
    COFD[3395] = 1.89917261E-02;
    COFD[3396] = -2.11401024E+01;
    COFD[3397] = 5.28856610E+00;
    COFD[3398] = -4.61764591E-01;
    COFD[3399] = 1.96258820E-02;
    COFD[3400] = -1.98257536E+01;
    COFD[3401] = 4.86038833E+00;
    COFD[3402] = -4.12187543E-01;
    COFD[3403] = 1.77155642E-02;
    COFD[3404] = -1.98573353E+01;
    COFD[3405] = 4.86038833E+00;
    COFD[3406] = -4.12187543E-01;
    COFD[3407] = 1.77155642E-02;
    COFD[3408] = -1.98573353E+01;
    COFD[3409] = 4.86038833E+00;
    COFD[3410] = -4.12187543E-01;
    COFD[3411] = 1.77155642E-02;
    COFD[3412] = -1.98614317E+01;
    COFD[3413] = 4.86038833E+00;
    COFD[3414] = -4.12187543E-01;
    COFD[3415] = 1.77155642E-02;
    COFD[3416] = -1.98614317E+01;
    COFD[3417] = 4.86038833E+00;
    COFD[3418] = -4.12187543E-01;
    COFD[3419] = 1.77155642E-02;
    COFD[3420] = -1.98870113E+01;
    COFD[3421] = 4.86038833E+00;
    COFD[3422] = -4.12187543E-01;
    COFD[3423] = 1.77155642E-02;
    COFD[3424] = -1.99085051E+01;
    COFD[3425] = 4.86038833E+00;
    COFD[3426] = -4.12187543E-01;
    COFD[3427] = 1.77155642E-02;
    COFD[3428] = -1.52554761E+01;
    COFD[3429] = 3.35922578E+00;
    COFD[3430] = -2.25181399E-01;
    COFD[3431] = 9.92132878E-03;
    COFD[3432] = -1.76147026E+01;
    COFD[3433] = 4.86049500E+00;
    COFD[3434] = -4.12200578E-01;
    COFD[3435] = 1.77160971E-02;
    COFD[3436] = -1.37794315E+01;
    COFD[3437] = 3.23973858E+00;
    COFD[3438] = -2.09989036E-01;
    COFD[3439] = 9.27667906E-03;
    COFD[3440] = -1.93015555E+01;
    COFD[3441] = 4.85015581E+00;
    COFD[3442] = -4.10945109E-01;
    COFD[3443] = 1.76651398E-02;
    COFD[3444] = -1.93015555E+01;
    COFD[3445] = 4.85015581E+00;
    COFD[3446] = -4.10945109E-01;
    COFD[3447] = 1.76651398E-02;
    COFD[3448] = -1.93276434E+01;
    COFD[3449] = 4.85015581E+00;
    COFD[3450] = -4.10945109E-01;
    COFD[3451] = 1.76651398E-02;
    COFD[3452] = -1.70534856E+01;
    COFD[3453] = 4.14240922E+00;
    COFD[3454] = -3.25239774E-01;
    COFD[3455] = 1.41980687E-02;
    COFD[3456] = -1.92867554E+01;
    COFD[3457] = 4.83375900E+00;
    COFD[3458] = -4.09146560E-01;
    COFD[3459] = 1.76006599E-02;
    COFD[3460] = -1.70757047E+01;
    COFD[3461] = 4.14240922E+00;
    COFD[3462] = -3.25239774E-01;
    COFD[3463] = 1.41980687E-02;
    COFD[3464] = -2.07653719E+01;
    COFD[3465] = 5.01092022E+00;
    COFD[3466] = -3.77985635E-01;
    COFD[3467] = 1.40968645E-02;
    COFD[3468] = -2.13961414E+01;
    COFD[3469] = 5.46685775E+00;
    COFD[3470] = -4.78665416E-01;
    COFD[3471] = 2.01093915E-02;
    COFD[3472] = -2.14079882E+01;
    COFD[3473] = 5.46685775E+00;
    COFD[3474] = -4.78665416E-01;
    COFD[3475] = 2.01093915E-02;
    COFD[3476] = -1.81735763E+01;
    COFD[3477] = 4.38391495E+00;
    COFD[3478] = -3.54941287E-01;
    COFD[3479] = 1.54195107E-02;
    COFD[3480] = -2.11031143E+01;
    COFD[3481] = 5.39439999E+00;
    COFD[3482] = -4.72050184E-01;
    COFD[3483] = 1.99336257E-02;
    COFD[3484] = -2.19215555E+01;
    COFD[3485] = 5.45216133E+00;
    COFD[3486] = -4.52916925E-01;
    COFD[3487] = 1.80456400E-02;
    COFD[3488] = -2.14049543E+01;
    COFD[3489] = 5.41122754E+00;
    COFD[3490] = -4.73185889E-01;
    COFD[3491] = 1.99407905E-02;
    COFD[3492] = -2.19317743E+01;
    COFD[3493] = 5.45216133E+00;
    COFD[3494] = -4.52916925E-01;
    COFD[3495] = 1.80456400E-02;
    COFD[3496] = -2.14151520E+01;
    COFD[3497] = 5.41122754E+00;
    COFD[3498] = -4.73185889E-01;
    COFD[3499] = 1.99407905E-02;
    COFD[3500] = -2.20421041E+01;
    COFD[3501] = 5.52708332E+00;
    COFD[3502] = -4.68000808E-01;
    COFD[3503] = 1.89131908E-02;
    COFD[3504] = -2.20421041E+01;
    COFD[3505] = 5.52708332E+00;
    COFD[3506] = -4.68000808E-01;
    COFD[3507] = 1.89131908E-02;
    COFD[3508] = -1.84688406E+01;
    COFD[3509] = 4.49330851E+00;
    COFD[3510] = -3.68208715E-01;
    COFD[3511] = 1.59565402E-02;
    COFD[3512] = -1.84777607E+01;
    COFD[3513] = 4.49330851E+00;
    COFD[3514] = -3.68208715E-01;
    COFD[3515] = 1.59565402E-02;
    COFD[3516] = -1.84863000E+01;
    COFD[3517] = 4.49330851E+00;
    COFD[3518] = -3.68208715E-01;
    COFD[3519] = 1.59565402E-02;
    COFD[3520] = -2.13425698E+01;
    COFD[3521] = 5.40460130E+00;
    COFD[3522] = -4.72718910E-01;
    COFD[3523] = 1.99362717E-02;
    COFD[3524] = -2.22235122E+01;
    COFD[3525] = 5.54251230E+00;
    COFD[3526] = -4.70946314E-01;
    COFD[3527] = 1.90785869E-02;
    COFD[3528] = -2.23017378E+01;
    COFD[3529] = 5.61208013E+00;
    COFD[3530] = -4.91433206E-01;
    COFD[3531] = 2.04241259E-02;
    COFD[3532] = -2.23150204E+01;
    COFD[3533] = 5.49916900E+00;
    COFD[3534] = -4.61818485E-01;
    COFD[3535] = 1.85431163E-02;
    COFD[3536] = -2.23071725E+01;
    COFD[3537] = 5.61208013E+00;
    COFD[3538] = -4.91433206E-01;
    COFD[3539] = 2.04241259E-02;
    COFD[3540] = -2.21083035E+01;
    COFD[3541] = 5.48540187E+00;
    COFD[3542] = -4.58962148E-01;
    COFD[3543] = 1.83770355E-02;
    COFD[3544] = -2.25096031E+01;
    COFD[3545] = 5.58506701E+00;
    COFD[3546] = -4.79642457E-01;
    COFD[3547] = 1.95823546E-02;
    COFD[3548] = -2.25131999E+01;
    COFD[3549] = 5.58506701E+00;
    COFD[3550] = -4.79642457E-01;
    COFD[3551] = 1.95823546E-02;
    COFD[3552] = -2.24868912E+01;
    COFD[3553] = 5.49895896E+00;
    COFD[3554] = -4.61764691E-01;
    COFD[3555] = 1.85397456E-02;
    COFD[3556] = -2.23744741E+01;
    COFD[3557] = 5.61208013E+00;
    COFD[3558] = -4.91433206E-01;
    COFD[3559] = 2.04241259E-02;
    COFD[3560] = -2.24116928E+01;
    COFD[3561] = 5.61208013E+00;
    COFD[3562] = -4.91433206E-01;
    COFD[3563] = 2.04241259E-02;
    COFD[3564] = -2.24116928E+01;
    COFD[3565] = 5.61208013E+00;
    COFD[3566] = -4.91433206E-01;
    COFD[3567] = 2.04241259E-02;
    COFD[3568] = -2.24165516E+01;
    COFD[3569] = 5.61208013E+00;
    COFD[3570] = -4.91433206E-01;
    COFD[3571] = 2.04241259E-02;
    COFD[3572] = -2.24165516E+01;
    COFD[3573] = 5.61208013E+00;
    COFD[3574] = -4.91433206E-01;
    COFD[3575] = 2.04241259E-02;
    COFD[3576] = -2.24470642E+01;
    COFD[3577] = 5.61208013E+00;
    COFD[3578] = -4.91433206E-01;
    COFD[3579] = 2.04241259E-02;
    COFD[3580] = -2.24729433E+01;
    COFD[3581] = 5.61208013E+00;
    COFD[3582] = -4.91433206E-01;
    COFD[3583] = 2.04241259E-02;
    COFD[3584] = -1.81432461E+01;
    COFD[3585] = 4.37565431E+00;
    COFD[3586] = -3.53906025E-01;
    COFD[3587] = 1.53760786E-02;
    COFD[3588] = -1.94694048E+01;
    COFD[3589] = 5.43830787E+00;
    COFD[3590] = -4.75472880E-01;
    COFD[3591] = 1.99909996E-02;
    COFD[3592] = -1.57045332E+01;
    COFD[3593] = 3.93614244E+00;
    COFD[3594] = -2.99111497E-01;
    COFD[3595] = 1.30888229E-02;
    COFD[3596] = -2.11406662E+01;
    COFD[3597] = 5.42846112E+00;
    COFD[3598] = -4.74321870E-01;
    COFD[3599] = 1.99459749E-02;
    COFD[3600] = -2.11406662E+01;
    COFD[3601] = 5.42846112E+00;
    COFD[3602] = -4.74321870E-01;
    COFD[3603] = 1.99459749E-02;
    COFD[3604] = -2.11667605E+01;
    COFD[3605] = 5.42846112E+00;
    COFD[3606] = -4.74321870E-01;
    COFD[3607] = 1.99459749E-02;
    COFD[3608] = -1.90946745E+01;
    COFD[3609] = 4.84384483E+00;
    COFD[3610] = -4.10265575E-01;
    COFD[3611] = 1.76414287E-02;
    COFD[3612] = -2.11372811E+01;
    COFD[3613] = 5.41773516E+00;
    COFD[3614] = -4.73414338E-01;
    COFD[3615] = 1.99258685E-02;
    COFD[3616] = -1.91168996E+01;
    COFD[3617] = 4.84384483E+00;
    COFD[3618] = -4.10265575E-01;
    COFD[3619] = 1.76414287E-02;
    COFD[3620] = -1.87453067E+01;
    COFD[3621] = 3.96926341E+00;
    COFD[3622] = -2.16412264E-01;
    COFD[3623] = 6.06012078E-03;
    COFD[3624] = -2.20351084E+01;
    COFD[3625] = 5.49663315E+00;
    COFD[3626] = -4.61182837E-01;
    COFD[3627] = 1.85035558E-02;
    COFD[3628] = -2.20469595E+01;
    COFD[3629] = 5.49663315E+00;
    COFD[3630] = -4.61182837E-01;
    COFD[3631] = 1.85035558E-02;
    COFD[3632] = -2.03015042E+01;
    COFD[3633] = 5.11106992E+00;
    COFD[3634] = -4.42047129E-01;
    COFD[3635] = 1.89042990E-02;
    COFD[3636] = -2.20490757E+01;
    COFD[3637] = 5.56049839E+00;
    COFD[3638] = -4.74367872E-01;
    COFD[3639] = 1.92702787E-02;
    COFD[3640] = -2.01009366E+01;
    COFD[3641] = 4.41511629E+00;
    COFD[3642] = -2.84086963E-01;
    COFD[3643] = 9.37586971E-03;
    COFD[3644] = -2.22424134E+01;
    COFD[3645] = 5.53139819E+00;
    COFD[3646] = -4.68828555E-01;
    COFD[3647] = 1.89597887E-02;
    COFD[3648] = -2.01111595E+01;
    COFD[3649] = 4.41511629E+00;
    COFD[3650] = -2.84086963E-01;
    COFD[3651] = 9.37586971E-03;
    COFD[3652] = -2.22526152E+01;
    COFD[3653] = 5.53139819E+00;
    COFD[3654] = -4.68828555E-01;
    COFD[3655] = 1.89597887E-02;
    COFD[3656] = -2.09320587E+01;
    COFD[3657] = 4.82184721E+00;
    COFD[3658] = -3.48128875E-01;
    COFD[3659] = 1.25918978E-02;
    COFD[3660] = -2.09320587E+01;
    COFD[3661] = 4.82184721E+00;
    COFD[3662] = -3.48128875E-01;
    COFD[3663] = 1.25918978E-02;
    COFD[3664] = -2.05284752E+01;
    COFD[3665] = 5.18417470E+00;
    COFD[3666] = -4.49491573E-01;
    COFD[3667] = 1.91438508E-02;
    COFD[3668] = -2.05373990E+01;
    COFD[3669] = 5.18417470E+00;
    COFD[3670] = -4.49491573E-01;
    COFD[3671] = 1.91438508E-02;
    COFD[3672] = -2.05459419E+01;
    COFD[3673] = 5.18417470E+00;
    COFD[3674] = -4.49491573E-01;
    COFD[3675] = 1.91438508E-02;
    COFD[3676] = -2.22235122E+01;
    COFD[3677] = 5.54251230E+00;
    COFD[3678] = -4.70946314E-01;
    COFD[3679] = 1.90785869E-02;
    COFD[3680] = -2.09236948E+01;
    COFD[3681] = 4.72895031E+00;
    COFD[3682] = -3.33332771E-01;
    COFD[3683] = 1.18431478E-02;
    COFD[3684] = -2.20535434E+01;
    COFD[3685] = 5.25559319E+00;
    COFD[3686] = -4.18539684E-01;
    COFD[3687] = 1.62026965E-02;
    COFD[3688] = -2.07244002E+01;
    COFD[3689] = 4.56211059E+00;
    COFD[3690] = -3.06895158E-01;
    COFD[3691] = 1.05100393E-02;
    COFD[3692] = -2.20589808E+01;
    COFD[3693] = 5.25559319E+00;
    COFD[3694] = -4.18539684E-01;
    COFD[3695] = 1.62026965E-02;
    COFD[3696] = -2.04156275E+01;
    COFD[3697] = 4.50250781E+00;
    COFD[3698] = -2.97622106E-01;
    COFD[3699] = 1.00481473E-02;
    COFD[3700] = -2.14564561E+01;
    COFD[3701] = 4.87969239E+00;
    COFD[3702] = -3.57259098E-01;
    COFD[3703] = 1.30518671E-02;
    COFD[3704] = -2.14600550E+01;
    COFD[3705] = 4.87969239E+00;
    COFD[3706] = -3.57259098E-01;
    COFD[3707] = 1.30518671E-02;
    COFD[3708] = -2.08906722E+01;
    COFD[3709] = 4.56059926E+00;
    COFD[3710] = -3.06658953E-01;
    COFD[3711] = 1.04982481E-02;
    COFD[3712] = -2.21263189E+01;
    COFD[3713] = 5.25559319E+00;
    COFD[3714] = -4.18539684E-01;
    COFD[3715] = 1.62026965E-02;
    COFD[3716] = -2.21635600E+01;
    COFD[3717] = 5.25559319E+00;
    COFD[3718] = -4.18539684E-01;
    COFD[3719] = 1.62026965E-02;
    COFD[3720] = -2.21635600E+01;
    COFD[3721] = 5.25559319E+00;
    COFD[3722] = -4.18539684E-01;
    COFD[3723] = 1.62026965E-02;
    COFD[3724] = -2.21684219E+01;
    COFD[3725] = 5.25559319E+00;
    COFD[3726] = -4.18539684E-01;
    COFD[3727] = 1.62026965E-02;
    COFD[3728] = -2.21684219E+01;
    COFD[3729] = 5.25559319E+00;
    COFD[3730] = -4.18539684E-01;
    COFD[3731] = 1.62026965E-02;
    COFD[3732] = -2.21989542E+01;
    COFD[3733] = 5.25559319E+00;
    COFD[3734] = -4.18539684E-01;
    COFD[3735] = 1.62026965E-02;
    COFD[3736] = -2.22248512E+01;
    COFD[3737] = 5.25559319E+00;
    COFD[3738] = -4.18539684E-01;
    COFD[3739] = 1.62026965E-02;
    COFD[3740] = -2.02738958E+01;
    COFD[3741] = 5.10426133E+00;
    COFD[3742] = -4.41256919E-01;
    COFD[3743] = 1.88737290E-02;
    COFD[3744] = -1.89222647E+01;
    COFD[3745] = 5.20981862E+00;
    COFD[3746] = -4.52489825E-01;
    COFD[3747] = 1.92609226E-02;
    COFD[3748] = -1.49070584E+01;
    COFD[3749] = 3.56974825E+00;
    COFD[3750] = -2.52221138E-01;
    COFD[3751] = 1.10819767E-02;
    COFD[3752] = -2.05539642E+01;
    COFD[3753] = 5.20087227E+00;
    COFD[3754] = -4.51444972E-01;
    COFD[3755] = 1.92201898E-02;
    COFD[3756] = -2.05539642E+01;
    COFD[3757] = 5.20087227E+00;
    COFD[3758] = -4.51444972E-01;
    COFD[3759] = 1.92201898E-02;
    COFD[3760] = -2.05802040E+01;
    COFD[3761] = 5.20087227E+00;
    COFD[3762] = -4.51444972E-01;
    COFD[3763] = 1.92201898E-02;
    COFD[3764] = -1.83038349E+01;
    COFD[3765] = 4.49977587E+00;
    COFD[3766] = -3.68989022E-01;
    COFD[3767] = 1.59879891E-02;
    COFD[3768] = -2.05236383E+01;
    COFD[3769] = 5.17771473E+00;
    COFD[3770] = -4.48750201E-01;
    COFD[3771] = 1.91155567E-02;
    COFD[3772] = -1.83261963E+01;
    COFD[3773] = 4.49977587E+00;
    COFD[3774] = -3.68989022E-01;
    COFD[3775] = 1.59879891E-02;
    COFD[3776] = -2.05013151E+01;
    COFD[3777] = 4.74677479E+00;
    COFD[3778] = -3.36160335E-01;
    COFD[3779] = 1.19858600E-02;
    COFD[3780] = -2.21632697E+01;
    COFD[3781] = 5.59047891E+00;
    COFD[3782] = -4.85359237E-01;
    COFD[3783] = 2.00302091E-02;
    COFD[3784] = -2.21752213E+01;
    COFD[3785] = 5.59047891E+00;
    COFD[3786] = -4.85359237E-01;
    COFD[3787] = 2.00302091E-02;
    COFD[3788] = -1.94622154E+01;
    COFD[3789] = 4.76145602E+00;
    COFD[3790] = -4.00684587E-01;
    COFD[3791] = 1.72708322E-02;
    COFD[3792] = -2.20752743E+01;
    COFD[3793] = 5.60509076E+00;
    COFD[3794] = -4.91228646E-01;
    COFD[3795] = 2.04425499E-02;
    COFD[3796] = -2.14954618E+01;
    COFD[3797] = 5.05104142E+00;
    COFD[3798] = -3.84428177E-01;
    COFD[3799] = 1.44249285E-02;
    COFD[3800] = -2.23333008E+01;
    COFD[3801] = 5.61213963E+00;
    COFD[3802] = -4.90953451E-01;
    COFD[3803] = 2.03841518E-02;
    COFD[3804] = -2.15057773E+01;
    COFD[3805] = 5.05104142E+00;
    COFD[3806] = -3.84428177E-01;
    COFD[3807] = 1.44249285E-02;
    COFD[3808] = -2.23435951E+01;
    COFD[3809] = 5.61213963E+00;
    COFD[3810] = -4.90953451E-01;
    COFD[3811] = 2.03841518E-02;
    COFD[3812] = -2.19952270E+01;
    COFD[3813] = 5.31843386E+00;
    COFD[3814] = -4.29198961E-01;
    COFD[3815] = 1.67645078E-02;
    COFD[3816] = -2.19952270E+01;
    COFD[3817] = 5.31843386E+00;
    COFD[3818] = -4.29198961E-01;
    COFD[3819] = 1.67645078E-02;
    COFD[3820] = -1.97454130E+01;
    COFD[3821] = 4.86038833E+00;
    COFD[3822] = -4.12187543E-01;
    COFD[3823] = 1.77155642E-02;
    COFD[3824] = -1.97544224E+01;
    COFD[3825] = 4.86038833E+00;
    COFD[3826] = -4.12187543E-01;
    COFD[3827] = 1.77155642E-02;
    COFD[3828] = -1.97630486E+01;
    COFD[3829] = 4.86038833E+00;
    COFD[3830] = -4.12187543E-01;
    COFD[3831] = 1.77155642E-02;
    COFD[3832] = -2.23017378E+01;
    COFD[3833] = 5.61208013E+00;
    COFD[3834] = -4.91433206E-01;
    COFD[3835] = 2.04241259E-02;
    COFD[3836] = -2.20535434E+01;
    COFD[3837] = 5.25559319E+00;
    COFD[3838] = -4.18539684E-01;
    COFD[3839] = 1.62026965E-02;
    COFD[3840] = -2.25751604E+01;
    COFD[3841] = 5.52719085E+00;
    COFD[3842] = -4.68021516E-01;
    COFD[3843] = 1.89143588E-02;
    COFD[3844] = -2.19952997E+01;
    COFD[3845] = 5.15225414E+00;
    COFD[3846] = -4.00960404E-01;
    COFD[3847] = 1.52761358E-02;
    COFD[3848] = -2.25806604E+01;
    COFD[3849] = 5.52719085E+00;
    COFD[3850] = -4.68021516E-01;
    COFD[3851] = 1.89143588E-02;
    COFD[3852] = -2.17504079E+01;
    COFD[3853] = 5.11317705E+00;
    COFD[3854] = -3.94506858E-01;
    COFD[3855] = 1.49416998E-02;
    COFD[3856] = -2.24427771E+01;
    COFD[3857] = 5.35063033E+00;
    COFD[3858] = -4.34734059E-01;
    COFD[3859] = 1.70582710E-02;
    COFD[3860] = -2.24464229E+01;
    COFD[3861] = 5.35063033E+00;
    COFD[3862] = -4.34734059E-01;
    COFD[3863] = 1.70582710E-02;
    COFD[3864] = -2.21563247E+01;
    COFD[3865] = 5.15132481E+00;
    COFD[3866] = -4.00804233E-01;
    COFD[3867] = 1.52679625E-02;
    COFD[3868] = -2.26488330E+01;
    COFD[3869] = 5.52719085E+00;
    COFD[3870] = -4.68021516E-01;
    COFD[3871] = 1.89143588E-02;
    COFD[3872] = -2.26865869E+01;
    COFD[3873] = 5.52719085E+00;
    COFD[3874] = -4.68021516E-01;
    COFD[3875] = 1.89143588E-02;
    COFD[3876] = -2.26865869E+01;
    COFD[3877] = 5.52719085E+00;
    COFD[3878] = -4.68021516E-01;
    COFD[3879] = 1.89143588E-02;
    COFD[3880] = -2.26915186E+01;
    COFD[3881] = 5.52719085E+00;
    COFD[3882] = -4.68021516E-01;
    COFD[3883] = 1.89143588E-02;
    COFD[3884] = -2.26915186E+01;
    COFD[3885] = 5.52719085E+00;
    COFD[3886] = -4.68021516E-01;
    COFD[3887] = 1.89143588E-02;
    COFD[3888] = -2.27225058E+01;
    COFD[3889] = 5.52719085E+00;
    COFD[3890] = -4.68021516E-01;
    COFD[3891] = 1.89143588E-02;
    COFD[3892] = -2.27488112E+01;
    COFD[3893] = 5.52719085E+00;
    COFD[3894] = -4.68021516E-01;
    COFD[3895] = 1.89143588E-02;
    COFD[3896] = -1.94344389E+01;
    COFD[3897] = 4.75414932E+00;
    COFD[3898] = -3.99807084E-01;
    COFD[3899] = 1.72356137E-02;
    COFD[3900] = -1.98808664E+01;
    COFD[3901] = 5.52555673E+00;
    COFD[3902] = -4.84999851E-01;
    COFD[3903] = 2.03334931E-02;
    COFD[3904] = -1.61121174E+01;
    COFD[3905] = 4.04227735E+00;
    COFD[3906] = -3.12745253E-01;
    COFD[3907] = 1.36756977E-02;
    COFD[3908] = -2.15351292E+01;
    COFD[3909] = 5.51982454E+00;
    COFD[3910] = -4.84452039E-01;
    COFD[3911] = 2.03175522E-02;
    COFD[3912] = -2.15351292E+01;
    COFD[3913] = 5.51982454E+00;
    COFD[3914] = -4.84452039E-01;
    COFD[3915] = 2.03175522E-02;
    COFD[3916] = -2.15615037E+01;
    COFD[3917] = 5.51982454E+00;
    COFD[3918] = -4.84452039E-01;
    COFD[3919] = 2.03175522E-02;
    COFD[3920] = -1.95342215E+01;
    COFD[3921] = 4.95249173E+00;
    COFD[3922] = -4.23376552E-01;
    COFD[3923] = 1.81703714E-02;
    COFD[3924] = -2.15095161E+01;
    COFD[3925] = 5.49964831E+00;
    COFD[3926] = -4.82275380E-01;
    COFD[3927] = 2.02405072E-02;
    COFD[3928] = -1.95567091E+01;
    COFD[3929] = 4.95249173E+00;
    COFD[3930] = -4.23376552E-01;
    COFD[3931] = 1.81703714E-02;
    COFD[3932] = -1.84568379E+01;
    COFD[3933] = 3.75912079E+00;
    COFD[3934] = -1.84235105E-01;
    COFD[3935] = 4.47800951E-03;
    COFD[3936] = -2.20585906E+01;
    COFD[3937] = 5.42445100E+00;
    COFD[3938] = -4.47918761E-01;
    COFD[3939] = 1.77729995E-02;
    COFD[3940] = -2.20706358E+01;
    COFD[3941] = 5.42445100E+00;
    COFD[3942] = -4.47918761E-01;
    COFD[3943] = 1.77729995E-02;
    COFD[3944] = -2.06106760E+01;
    COFD[3945] = 5.16748146E+00;
    COFD[3946] = -4.47594939E-01;
    COFD[3947] = 1.90724110E-02;
    COFD[3948] = -2.21129692E+01;
    COFD[3949] = 5.50506115E+00;
    COFD[3950] = -4.63563533E-01;
    COFD[3951] = 1.86575247E-02;
    COFD[3952] = -1.97587116E+01;
    COFD[3953] = 4.18758010E+00;
    COFD[3954] = -2.49327776E-01;
    COFD[3955] = 7.66559103E-03;
    COFD[3956] = -2.23342576E+01;
    COFD[3957] = 5.49239750E+00;
    COFD[3958] = -4.60320987E-01;
    COFD[3959] = 1.84538922E-02;
    COFD[3960] = -1.97691133E+01;
    COFD[3961] = 4.18758010E+00;
    COFD[3962] = -2.49327776E-01;
    COFD[3963] = 7.66559103E-03;
    COFD[3964] = -2.23446381E+01;
    COFD[3965] = 5.49239750E+00;
    COFD[3966] = -4.60320987E-01;
    COFD[3967] = 1.84538922E-02;
    COFD[3968] = -2.07399014E+01;
    COFD[3969] = 4.65728078E+00;
    COFD[3970] = -3.22002062E-01;
    COFD[3971] = 1.12723316E-02;
    COFD[3972] = -2.07399014E+01;
    COFD[3973] = 4.65728078E+00;
    COFD[3974] = -3.22002062E-01;
    COFD[3975] = 1.12723316E-02;
    COFD[3976] = -2.09524743E+01;
    COFD[3977] = 5.28755355E+00;
    COFD[3978] = -4.61641920E-01;
    COFD[3979] = 1.96208961E-02;
    COFD[3980] = -2.09615635E+01;
    COFD[3981] = 5.28755355E+00;
    COFD[3982] = -4.61641920E-01;
    COFD[3983] = 1.96208961E-02;
    COFD[3984] = -2.09702675E+01;
    COFD[3985] = 5.28755355E+00;
    COFD[3986] = -4.61641920E-01;
    COFD[3987] = 1.96208961E-02;
    COFD[3988] = -2.23150204E+01;
    COFD[3989] = 5.49916900E+00;
    COFD[3990] = -4.61818485E-01;
    COFD[3991] = 1.85431163E-02;
    COFD[3992] = -2.07244002E+01;
    COFD[3993] = 4.56211059E+00;
    COFD[3994] = -3.06895158E-01;
    COFD[3995] = 1.05100393E-02;
    COFD[3996] = -2.19952997E+01;
    COFD[3997] = 5.15225414E+00;
    COFD[3998] = -4.00960404E-01;
    COFD[3999] = 1.52761358E-02;
    COFD[4000] = -2.04327449E+01;
    COFD[4001] = 4.35883159E+00;
    COFD[4002] = -2.75434484E-01;
    COFD[4003] = 8.94819804E-03;
    COFD[4004] = -2.20008583E+01;
    COFD[4005] = 5.15225414E+00;
    COFD[4006] = -4.00960404E-01;
    COFD[4007] = 1.52761358E-02;
    COFD[4008] = -2.01151169E+01;
    COFD[4009] = 4.29138907E+00;
    COFD[4010] = -2.65108149E-01;
    COFD[4011] = 8.43949637E-03;
    COFD[4012] = -2.12543792E+01;
    COFD[4013] = 4.71504706E+00;
    COFD[4014] = -3.31128704E-01;
    COFD[4015] = 1.17319185E-02;
    COFD[4016] = -2.12580690E+01;
    COFD[4017] = 4.71504706E+00;
    COFD[4018] = -3.31128704E-01;
    COFD[4019] = 1.17319185E-02;
    COFD[4020] = -2.05947670E+01;
    COFD[4021] = 4.35703694E+00;
    COFD[4022] = -2.75158780E-01;
    COFD[4023] = 8.93458153E-03;
    COFD[4024] = -2.20698134E+01;
    COFD[4025] = 5.15225414E+00;
    COFD[4026] = -4.00960404E-01;
    COFD[4027] = 1.52761358E-02;
    COFD[4028] = -2.21080494E+01;
    COFD[4029] = 5.15225414E+00;
    COFD[4030] = -4.00960404E-01;
    COFD[4031] = 1.52761358E-02;
    COFD[4032] = -2.21080494E+01;
    COFD[4033] = 5.15225414E+00;
    COFD[4034] = -4.00960404E-01;
    COFD[4035] = 1.52761358E-02;
    COFD[4036] = -2.21130468E+01;
    COFD[4037] = 5.15225414E+00;
    COFD[4038] = -4.00960404E-01;
    COFD[4039] = 1.52761358E-02;
    COFD[4040] = -2.21130468E+01;
    COFD[4041] = 5.15225414E+00;
    COFD[4042] = -4.00960404E-01;
    COFD[4043] = 1.52761358E-02;
    COFD[4044] = -2.21444625E+01;
    COFD[4045] = 5.15225414E+00;
    COFD[4046] = -4.00960404E-01;
    COFD[4047] = 1.52761358E-02;
    COFD[4048] = -2.21711534E+01;
    COFD[4049] = 5.15225414E+00;
    COFD[4050] = -4.00960404E-01;
    COFD[4051] = 1.52761358E-02;
    COFD[4052] = -2.05842618E+01;
    COFD[4053] = 5.16117916E+00;
    COFD[4054] = -4.46897404E-01;
    COFD[4055] = 1.90470443E-02;
    COFD[4056] = -1.89225041E+01;
    COFD[4057] = 5.20981862E+00;
    COFD[4058] = -4.52489825E-01;
    COFD[4059] = 1.92609226E-02;
    COFD[4060] = -1.49075271E+01;
    COFD[4061] = 3.56974825E+00;
    COFD[4062] = -2.52221138E-01;
    COFD[4063] = 1.10819767E-02;
    COFD[4064] = -2.05565679E+01;
    COFD[4065] = 5.20087227E+00;
    COFD[4066] = -4.51444972E-01;
    COFD[4067] = 1.92201898E-02;
    COFD[4068] = -2.05565679E+01;
    COFD[4069] = 5.20087227E+00;
    COFD[4070] = -4.51444972E-01;
    COFD[4071] = 1.92201898E-02;
    COFD[4072] = -2.05829484E+01;
    COFD[4073] = 5.20087227E+00;
    COFD[4074] = -4.51444972E-01;
    COFD[4075] = 1.92201898E-02;
    COFD[4076] = -1.83067096E+01;
    COFD[4077] = 4.49977587E+00;
    COFD[4078] = -3.68989022E-01;
    COFD[4079] = 1.59879891E-02;
    COFD[4080] = -2.05265188E+01;
    COFD[4081] = 5.17771473E+00;
    COFD[4082] = -4.48750201E-01;
    COFD[4083] = 1.91155567E-02;
    COFD[4084] = -1.83292028E+01;
    COFD[4085] = 4.49977587E+00;
    COFD[4086] = -3.68989022E-01;
    COFD[4087] = 1.59879891E-02;
    COFD[4088] = -2.05044494E+01;
    COFD[4089] = 4.74677479E+00;
    COFD[4090] = -3.36160335E-01;
    COFD[4091] = 1.19858600E-02;
    COFD[4092] = -2.21672921E+01;
    COFD[4093] = 5.59047891E+00;
    COFD[4094] = -4.85359237E-01;
    COFD[4095] = 2.00302091E-02;
    COFD[4096] = -2.21793415E+01;
    COFD[4097] = 5.59047891E+00;
    COFD[4098] = -4.85359237E-01;
    COFD[4099] = 2.00302091E-02;
    COFD[4100] = -1.94664266E+01;
    COFD[4101] = 4.76145602E+00;
    COFD[4102] = -4.00684587E-01;
    COFD[4103] = 1.72708322E-02;
    COFD[4104] = -2.20794895E+01;
    COFD[4105] = 5.60509076E+00;
    COFD[4106] = -4.91228646E-01;
    COFD[4107] = 2.04425499E-02;
    COFD[4108] = -2.14997655E+01;
    COFD[4109] = 5.05104142E+00;
    COFD[4110] = -3.84428177E-01;
    COFD[4111] = 1.44249285E-02;
    COFD[4112] = -2.23376085E+01;
    COFD[4113] = 5.61213963E+00;
    COFD[4114] = -4.90953451E-01;
    COFD[4115] = 2.03841518E-02;
    COFD[4116] = -2.15101711E+01;
    COFD[4117] = 5.05104142E+00;
    COFD[4118] = -3.84428177E-01;
    COFD[4119] = 1.44249285E-02;
    COFD[4120] = -2.23479928E+01;
    COFD[4121] = 5.61213963E+00;
    COFD[4122] = -4.90953451E-01;
    COFD[4123] = 2.03841518E-02;
    COFD[4124] = -2.19997085E+01;
    COFD[4125] = 5.31843386E+00;
    COFD[4126] = -4.29198961E-01;
    COFD[4127] = 1.67645078E-02;
    COFD[4128] = -2.19997085E+01;
    COFD[4129] = 5.31843386E+00;
    COFD[4130] = -4.29198961E-01;
    COFD[4131] = 1.67645078E-02;
    COFD[4132] = -1.97499763E+01;
    COFD[4133] = 4.86038833E+00;
    COFD[4134] = -4.12187543E-01;
    COFD[4135] = 1.77155642E-02;
    COFD[4136] = -1.97590691E+01;
    COFD[4137] = 4.86038833E+00;
    COFD[4138] = -4.12187543E-01;
    COFD[4139] = 1.77155642E-02;
    COFD[4140] = -1.97677766E+01;
    COFD[4141] = 4.86038833E+00;
    COFD[4142] = -4.12187543E-01;
    COFD[4143] = 1.77155642E-02;
    COFD[4144] = -2.23071725E+01;
    COFD[4145] = 5.61208013E+00;
    COFD[4146] = -4.91433206E-01;
    COFD[4147] = 2.04241259E-02;
    COFD[4148] = -2.20589808E+01;
    COFD[4149] = 5.25559319E+00;
    COFD[4150] = -4.18539684E-01;
    COFD[4151] = 1.62026965E-02;
    COFD[4152] = -2.25806604E+01;
    COFD[4153] = 5.52719085E+00;
    COFD[4154] = -4.68021516E-01;
    COFD[4155] = 1.89143588E-02;
    COFD[4156] = -2.20008583E+01;
    COFD[4157] = 5.15225414E+00;
    COFD[4158] = -4.00960404E-01;
    COFD[4159] = 1.52761358E-02;
    COFD[4160] = -2.25862216E+01;
    COFD[4161] = 5.52719085E+00;
    COFD[4162] = -4.68021516E-01;
    COFD[4163] = 1.89143588E-02;
    COFD[4164] = -2.17560264E+01;
    COFD[4165] = 5.11317705E+00;
    COFD[4166] = -3.94506858E-01;
    COFD[4167] = 1.49416998E-02;
    COFD[4168] = -2.24490205E+01;
    COFD[4169] = 5.35063033E+00;
    COFD[4170] = -4.34734059E-01;
    COFD[4171] = 1.70582710E-02;
    COFD[4172] = -2.24527123E+01;
    COFD[4173] = 5.35063033E+00;
    COFD[4174] = -4.34734059E-01;
    COFD[4175] = 1.70582710E-02;
    COFD[4176] = -2.21626591E+01;
    COFD[4177] = 5.15132481E+00;
    COFD[4178] = -4.00804233E-01;
    COFD[4179] = 1.52679625E-02;
    COFD[4180] = -2.26552117E+01;
    COFD[4181] = 5.52719085E+00;
    COFD[4182] = -4.68021516E-01;
    COFD[4183] = 1.89143588E-02;
    COFD[4184] = -2.26934694E+01;
    COFD[4185] = 5.52719085E+00;
    COFD[4186] = -4.68021516E-01;
    COFD[4187] = 1.89143588E-02;
    COFD[4188] = -2.26934694E+01;
    COFD[4189] = 5.52719085E+00;
    COFD[4190] = -4.68021516E-01;
    COFD[4191] = 1.89143588E-02;
    COFD[4192] = -2.26984698E+01;
    COFD[4193] = 5.52719085E+00;
    COFD[4194] = -4.68021516E-01;
    COFD[4195] = 1.89143588E-02;
    COFD[4196] = -2.26984698E+01;
    COFD[4197] = 5.52719085E+00;
    COFD[4198] = -4.68021516E-01;
    COFD[4199] = 1.89143588E-02;
    COFD[4200] = -2.27299046E+01;
    COFD[4201] = 5.52719085E+00;
    COFD[4202] = -4.68021516E-01;
    COFD[4203] = 1.89143588E-02;
    COFD[4204] = -2.27566129E+01;
    COFD[4205] = 5.52719085E+00;
    COFD[4206] = -4.68021516E-01;
    COFD[4207] = 1.89143588E-02;
    COFD[4208] = -1.94386503E+01;
    COFD[4209] = 4.75414932E+00;
    COFD[4210] = -3.99807084E-01;
    COFD[4211] = 1.72356137E-02;
    COFD[4212] = -1.96914944E+01;
    COFD[4213] = 5.54637286E+00;
    COFD[4214] = -4.87070324E-01;
    COFD[4215] = 2.03983467E-02;
    COFD[4216] = -1.59632479E+01;
    COFD[4217] = 4.07051484E+00;
    COFD[4218] = -3.16303109E-01;
    COFD[4219] = 1.38259377E-02;
    COFD[4220] = -2.14048982E+01;
    COFD[4221] = 5.54007827E+00;
    COFD[4222] = -4.86434511E-01;
    COFD[4223] = 2.03779006E-02;
    COFD[4224] = -2.14048982E+01;
    COFD[4225] = 5.54007827E+00;
    COFD[4226] = -4.86434511E-01;
    COFD[4227] = 2.03779006E-02;
    COFD[4228] = -2.14314090E+01;
    COFD[4229] = 5.54007827E+00;
    COFD[4230] = -4.86434511E-01;
    COFD[4231] = 2.03779006E-02;
    COFD[4232] = -1.93925667E+01;
    COFD[4233] = 4.98286777E+00;
    COFD[4234] = -4.26970814E-01;
    COFD[4235] = 1.83122917E-02;
    COFD[4236] = -2.13881945E+01;
    COFD[4237] = 5.52422470E+00;
    COFD[4238] = -4.84872944E-01;
    COFD[4239] = 2.03298213E-02;
    COFD[4240] = -1.94151822E+01;
    COFD[4241] = 4.98286777E+00;
    COFD[4242] = -4.26970814E-01;
    COFD[4243] = 1.83122917E-02;
    COFD[4244] = -1.80862867E+01;
    COFD[4245] = 3.69199168E+00;
    COFD[4246] = -1.74005516E-01;
    COFD[4247] = 3.97694372E-03;
    COFD[4248] = -2.18318278E+01;
    COFD[4249] = 5.40298848E+00;
    COFD[4250] = -4.43954594E-01;
    COFD[4251] = 1.75542998E-02;
    COFD[4252] = -2.18439681E+01;
    COFD[4253] = 5.40298848E+00;
    COFD[4254] = -4.43954594E-01;
    COFD[4255] = 1.75542998E-02;
    COFD[4256] = -2.04949373E+01;
    COFD[4257] = 5.19614628E+00;
    COFD[4258] = -4.50889164E-01;
    COFD[4259] = 1.91983328E-02;
    COFD[4260] = -2.19162360E+01;
    COFD[4261] = 5.49906960E+00;
    COFD[4262] = -4.61793001E-01;
    COFD[4263] = 1.85415189E-02;
    COFD[4264] = -1.93946947E+01;
    COFD[4265] = 4.10954793E+00;
    COFD[4266] = -2.37523329E-01;
    COFD[4267] = 7.08858141E-03;
    COFD[4268] = -2.21229141E+01;
    COFD[4269] = 5.47072190E+00;
    COFD[4270] = -4.56301261E-01;
    COFD[4271] = 1.82313566E-02;
    COFD[4272] = -1.94051843E+01;
    COFD[4273] = 4.10954793E+00;
    COFD[4274] = -2.37523329E-01;
    COFD[4275] = 7.08858141E-03;
    COFD[4276] = -2.21333822E+01;
    COFD[4277] = 5.47072190E+00;
    COFD[4278] = -4.56301261E-01;
    COFD[4279] = 1.82313566E-02;
    COFD[4280] = -2.04451935E+01;
    COFD[4281] = 4.60682543E+00;
    COFD[4282] = -3.13971634E-01;
    COFD[4283] = 1.08661011E-02;
    COFD[4284] = -2.04451935E+01;
    COFD[4285] = 4.60682543E+00;
    COFD[4286] = -3.13971634E-01;
    COFD[4287] = 1.08661011E-02;
    COFD[4288] = -2.08463209E+01;
    COFD[4289] = 5.32244593E+00;
    COFD[4290] = -4.65829403E-01;
    COFD[4291] = 1.97895274E-02;
    COFD[4292] = -2.08554914E+01;
    COFD[4293] = 5.32244593E+00;
    COFD[4294] = -4.65829403E-01;
    COFD[4295] = 1.97895274E-02;
    COFD[4296] = -2.08642748E+01;
    COFD[4297] = 5.32244593E+00;
    COFD[4298] = -4.65829403E-01;
    COFD[4299] = 1.97895274E-02;
    COFD[4300] = -2.21083035E+01;
    COFD[4301] = 5.48540187E+00;
    COFD[4302] = -4.58962148E-01;
    COFD[4303] = 1.83770355E-02;
    COFD[4304] = -2.04156275E+01;
    COFD[4305] = 4.50250781E+00;
    COFD[4306] = -2.97622106E-01;
    COFD[4307] = 1.00481473E-02;
    COFD[4308] = -2.17504079E+01;
    COFD[4309] = 5.11317705E+00;
    COFD[4310] = -3.94506858E-01;
    COFD[4311] = 1.49416998E-02;
    COFD[4312] = -2.01151169E+01;
    COFD[4313] = 4.29138907E+00;
    COFD[4314] = -2.65108149E-01;
    COFD[4315] = 8.43949637E-03;
    COFD[4316] = -2.17560264E+01;
    COFD[4317] = 5.11317705E+00;
    COFD[4318] = -3.94506858E-01;
    COFD[4319] = 1.49416998E-02;
    COFD[4320] = -1.97704178E+01;
    COFD[4321] = 4.22062499E+00;
    COFD[4322] = -2.54326872E-01;
    COFD[4323] = 7.91017784E-03;
    COFD[4324] = -2.09786696E+01;
    COFD[4325] = 4.66190324E+00;
    COFD[4326] = -3.22729985E-01;
    COFD[4327] = 1.13089121E-02;
    COFD[4328] = -2.09824046E+01;
    COFD[4329] = 4.66190324E+00;
    COFD[4330] = -3.22729985E-01;
    COFD[4331] = 1.13089121E-02;
    COFD[4332] = -2.02873532E+01;
    COFD[4333] = 4.28958705E+00;
    COFD[4334] = -2.64832562E-01;
    COFD[4335] = 8.42593887E-03;
    COFD[4336] = -2.18257833E+01;
    COFD[4337] = 5.11317705E+00;
    COFD[4338] = -3.94506858E-01;
    COFD[4339] = 1.49416998E-02;
    COFD[4340] = -2.18645147E+01;
    COFD[4341] = 5.11317705E+00;
    COFD[4342] = -3.94506858E-01;
    COFD[4343] = 1.49416998E-02;
    COFD[4344] = -2.18645147E+01;
    COFD[4345] = 5.11317705E+00;
    COFD[4346] = -3.94506858E-01;
    COFD[4347] = 1.49416998E-02;
    COFD[4348] = -2.18695798E+01;
    COFD[4349] = 5.11317705E+00;
    COFD[4350] = -3.94506858E-01;
    COFD[4351] = 1.49416998E-02;
    COFD[4352] = -2.18695798E+01;
    COFD[4353] = 5.11317705E+00;
    COFD[4354] = -3.94506858E-01;
    COFD[4355] = 1.49416998E-02;
    COFD[4356] = -2.19014365E+01;
    COFD[4357] = 5.11317705E+00;
    COFD[4358] = -3.94506858E-01;
    COFD[4359] = 1.49416998E-02;
    COFD[4360] = -2.19285250E+01;
    COFD[4361] = 5.11317705E+00;
    COFD[4362] = -3.94506858E-01;
    COFD[4363] = 1.49416998E-02;
    COFD[4364] = -2.04649069E+01;
    COFD[4365] = 5.18856872E+00;
    COFD[4366] = -4.50001829E-01;
    COFD[4367] = 1.91636142E-02;
    COFD[4368] = -1.95591756E+01;
    COFD[4369] = 5.40012228E+00;
    COFD[4370] = -4.72416152E-01;
    COFD[4371] = 1.99342987E-02;
    COFD[4372] = -1.56146754E+01;
    COFD[4373] = 3.82176300E+00;
    COFD[4374] = -2.84202467E-01;
    COFD[4375] = 1.24384328E-02;
    COFD[4376] = -2.12255802E+01;
    COFD[4377] = 5.39710574E+00;
    COFD[4378] = -4.72221942E-01;
    COFD[4379] = 1.99338345E-02;
    COFD[4380] = -2.12255802E+01;
    COFD[4381] = 5.39710574E+00;
    COFD[4382] = -4.72221942E-01;
    COFD[4383] = 1.99338345E-02;
    COFD[4384] = -2.12534275E+01;
    COFD[4385] = 5.39710574E+00;
    COFD[4386] = -4.72221942E-01;
    COFD[4387] = 1.99338345E-02;
    COFD[4388] = -1.91187832E+01;
    COFD[4389] = 4.76920246E+00;
    COFD[4390] = -4.01609656E-01;
    COFD[4391] = 1.73077246E-02;
    COFD[4392] = -2.12314042E+01;
    COFD[4393] = 5.38826297E+00;
    COFD[4394] = -4.71579030E-01;
    COFD[4395] = 1.99262449E-02;
    COFD[4396] = -1.91426599E+01;
    COFD[4397] = 4.76920246E+00;
    COFD[4398] = -4.01609656E-01;
    COFD[4399] = 1.73077246E-02;
    COFD[4400] = -1.95040271E+01;
    COFD[4401] = 4.21243642E+00;
    COFD[4402] = -2.53087979E-01;
    COFD[4403] = 7.84955719E-03;
    COFD[4404] = -2.22849472E+01;
    COFD[4405] = 5.53187091E+00;
    COFD[4406] = -4.68918850E-01;
    COFD[4407] = 1.89648616E-02;
    COFD[4408] = -2.22980489E+01;
    COFD[4409] = 5.53187091E+00;
    COFD[4410] = -4.68918850E-01;
    COFD[4411] = 1.89648616E-02;
    COFD[4412] = -2.02929645E+01;
    COFD[4413] = 5.02679716E+00;
    COFD[4414] = -4.32175153E-01;
    COFD[4415] = 1.85182518E-02;
    COFD[4416] = -2.22838940E+01;
    COFD[4417] = 5.58518325E+00;
    COFD[4418] = -4.80534710E-01;
    COFD[4419] = 1.96556785E-02;
    COFD[4420] = -2.06900561E+01;
    COFD[4421] = 4.59219167E+00;
    COFD[4422] = -3.11627566E-01;
    COFD[4423] = 1.07472587E-02;
    COFD[4424] = -2.25191362E+01;
    COFD[4425] = 5.58153976E+00;
    COFD[4426] = -4.78600113E-01;
    COFD[4427] = 1.95139120E-02;
    COFD[4428] = -2.07014385E+01;
    COFD[4429] = 4.59219167E+00;
    COFD[4430] = -3.11627566E-01;
    COFD[4431] = 1.07472587E-02;
    COFD[4432] = -2.25304962E+01;
    COFD[4433] = 5.58153976E+00;
    COFD[4434] = -4.78600113E-01;
    COFD[4435] = 1.95139120E-02;
    COFD[4436] = -2.14849242E+01;
    COFD[4437] = 4.98394310E+00;
    COFD[4438] = -3.73694768E-01;
    COFD[4439] = 1.38795565E-02;
    COFD[4440] = -2.14849242E+01;
    COFD[4441] = 4.98394310E+00;
    COFD[4442] = -3.73694768E-01;
    COFD[4443] = 1.38795565E-02;
    COFD[4444] = -2.06003483E+01;
    COFD[4445] = 5.13220015E+00;
    COFD[4446] = -4.44412321E-01;
    COFD[4447] = 1.89917261E-02;
    COFD[4448] = -2.06103502E+01;
    COFD[4449] = 5.13220015E+00;
    COFD[4450] = -4.44412321E-01;
    COFD[4451] = 1.89917261E-02;
    COFD[4452] = -2.06199456E+01;
    COFD[4453] = 5.13220015E+00;
    COFD[4454] = -4.44412321E-01;
    COFD[4455] = 1.89917261E-02;
    COFD[4456] = -2.25096031E+01;
    COFD[4457] = 5.58506701E+00;
    COFD[4458] = -4.79642457E-01;
    COFD[4459] = 1.95823546E-02;
    COFD[4460] = -2.14564561E+01;
    COFD[4461] = 4.87969239E+00;
    COFD[4462] = -3.57259098E-01;
    COFD[4463] = 1.30518671E-02;
    COFD[4464] = -2.24427771E+01;
    COFD[4465] = 5.35063033E+00;
    COFD[4466] = -4.34734059E-01;
    COFD[4467] = 1.70582710E-02;
    COFD[4468] = -2.12543792E+01;
    COFD[4469] = 4.71504706E+00;
    COFD[4470] = -3.31128704E-01;
    COFD[4471] = 1.17319185E-02;
    COFD[4472] = -2.24490205E+01;
    COFD[4473] = 5.35063033E+00;
    COFD[4474] = -4.34734059E-01;
    COFD[4475] = 1.70582710E-02;
    COFD[4476] = -2.09786696E+01;
    COFD[4477] = 4.66190324E+00;
    COFD[4478] = -3.22729985E-01;
    COFD[4479] = 1.13089121E-02;
    COFD[4480] = -2.20007672E+01;
    COFD[4481] = 5.03720506E+00;
    COFD[4482] = -3.82190496E-01;
    COFD[4483] = 1.43104209E-02;
    COFD[4484] = -2.20049811E+01;
    COFD[4485] = 5.03720506E+00;
    COFD[4486] = -3.82190496E-01;
    COFD[4487] = 1.43104209E-02;
    COFD[4488] = -2.14229024E+01;
    COFD[4489] = 4.71363498E+00;
    COFD[4490] = -3.30905003E-01;
    COFD[4491] = 1.17206344E-02;
    COFD[4492] = -2.25272268E+01;
    COFD[4493] = 5.35063033E+00;
    COFD[4494] = -4.34734059E-01;
    COFD[4495] = 1.70582710E-02;
    COFD[4496] = -2.25712566E+01;
    COFD[4497] = 5.35063033E+00;
    COFD[4498] = -4.34734059E-01;
    COFD[4499] = 1.70582710E-02;
    COFD[4500] = -2.25712566E+01;
    COFD[4501] = 5.35063033E+00;
    COFD[4502] = -4.34734059E-01;
    COFD[4503] = 1.70582710E-02;
    COFD[4504] = -2.25770498E+01;
    COFD[4505] = 5.35063033E+00;
    COFD[4506] = -4.34734059E-01;
    COFD[4507] = 1.70582710E-02;
    COFD[4508] = -2.25770498E+01;
    COFD[4509] = 5.35063033E+00;
    COFD[4510] = -4.34734059E-01;
    COFD[4511] = 1.70582710E-02;
    COFD[4512] = -2.26136853E+01;
    COFD[4513] = 5.35063033E+00;
    COFD[4514] = -4.34734059E-01;
    COFD[4515] = 1.70582710E-02;
    COFD[4516] = -2.26451233E+01;
    COFD[4517] = 5.35063033E+00;
    COFD[4518] = -4.34734059E-01;
    COFD[4519] = 1.70582710E-02;
    COFD[4520] = -2.02659445E+01;
    COFD[4521] = 5.01994886E+00;
    COFD[4522] = -4.31380254E-01;
    COFD[4523] = 1.84875367E-02;
    COFD[4524] = -1.95593165E+01;
    COFD[4525] = 5.40012228E+00;
    COFD[4526] = -4.72416152E-01;
    COFD[4527] = 1.99342987E-02;
    COFD[4528] = -1.56149526E+01;
    COFD[4529] = 3.82176300E+00;
    COFD[4530] = -2.84202467E-01;
    COFD[4531] = 1.24384328E-02;
    COFD[4532] = -2.12271939E+01;
    COFD[4533] = 5.39710574E+00;
    COFD[4534] = -4.72221942E-01;
    COFD[4535] = 1.99338345E-02;
    COFD[4536] = -2.12271939E+01;
    COFD[4537] = 5.39710574E+00;
    COFD[4538] = -4.72221942E-01;
    COFD[4539] = 1.99338345E-02;
    COFD[4540] = -2.12551337E+01;
    COFD[4541] = 5.39710574E+00;
    COFD[4542] = -4.72221942E-01;
    COFD[4543] = 1.99338345E-02;
    COFD[4544] = -1.91205757E+01;
    COFD[4545] = 4.76920246E+00;
    COFD[4546] = -4.01609656E-01;
    COFD[4547] = 1.73077246E-02;
    COFD[4548] = -2.12332005E+01;
    COFD[4549] = 5.38826297E+00;
    COFD[4550] = -4.71579030E-01;
    COFD[4551] = 1.99262449E-02;
    COFD[4552] = -1.91445402E+01;
    COFD[4553] = 4.76920246E+00;
    COFD[4554] = -4.01609656E-01;
    COFD[4555] = 1.73077246E-02;
    COFD[4556] = -1.95059929E+01;
    COFD[4557] = 4.21243642E+00;
    COFD[4558] = -2.53087979E-01;
    COFD[4559] = 7.84955719E-03;
    COFD[4560] = -2.22875222E+01;
    COFD[4561] = 5.53187091E+00;
    COFD[4562] = -4.68918850E-01;
    COFD[4563] = 1.89648616E-02;
    COFD[4564] = -2.23006924E+01;
    COFD[4565] = 5.53187091E+00;
    COFD[4566] = -4.68918850E-01;
    COFD[4567] = 1.89648616E-02;
    COFD[4568] = -2.02956721E+01;
    COFD[4569] = 5.02679716E+00;
    COFD[4570] = -4.32175153E-01;
    COFD[4571] = 1.85182518E-02;
    COFD[4572] = -2.22866045E+01;
    COFD[4573] = 5.58518325E+00;
    COFD[4574] = -4.80534710E-01;
    COFD[4575] = 1.96556785E-02;
    COFD[4576] = -2.06928292E+01;
    COFD[4577] = 4.59219167E+00;
    COFD[4578] = -3.11627566E-01;
    COFD[4579] = 1.07472587E-02;
    COFD[4580] = -2.25219121E+01;
    COFD[4581] = 5.58153976E+00;
    COFD[4582] = -4.78600113E-01;
    COFD[4583] = 1.95139120E-02;
    COFD[4584] = -2.07042756E+01;
    COFD[4585] = 4.59219167E+00;
    COFD[4586] = -3.11627566E-01;
    COFD[4587] = 1.07472587E-02;
    COFD[4588] = -2.25333361E+01;
    COFD[4589] = 5.58153976E+00;
    COFD[4590] = -4.78600113E-01;
    COFD[4591] = 1.95139120E-02;
    COFD[4592] = -2.14878240E+01;
    COFD[4593] = 4.98394310E+00;
    COFD[4594] = -3.73694768E-01;
    COFD[4595] = 1.38795565E-02;
    COFD[4596] = -2.14878240E+01;
    COFD[4597] = 4.98394310E+00;
    COFD[4598] = -3.73694768E-01;
    COFD[4599] = 1.38795565E-02;
    COFD[4600] = -2.06033066E+01;
    COFD[4601] = 5.13220015E+00;
    COFD[4602] = -4.44412321E-01;
    COFD[4603] = 1.89917261E-02;
    COFD[4604] = -2.06133685E+01;
    COFD[4605] = 5.13220015E+00;
    COFD[4606] = -4.44412321E-01;
    COFD[4607] = 1.89917261E-02;
    COFD[4608] = -2.06230226E+01;
    COFD[4609] = 5.13220015E+00;
    COFD[4610] = -4.44412321E-01;
    COFD[4611] = 1.89917261E-02;
    COFD[4612] = -2.25131999E+01;
    COFD[4613] = 5.58506701E+00;
    COFD[4614] = -4.79642457E-01;
    COFD[4615] = 1.95823546E-02;
    COFD[4616] = -2.14600550E+01;
    COFD[4617] = 4.87969239E+00;
    COFD[4618] = -3.57259098E-01;
    COFD[4619] = 1.30518671E-02;
    COFD[4620] = -2.24464229E+01;
    COFD[4621] = 5.35063033E+00;
    COFD[4622] = -4.34734059E-01;
    COFD[4623] = 1.70582710E-02;
    COFD[4624] = -2.12580690E+01;
    COFD[4625] = 4.71504706E+00;
    COFD[4626] = -3.31128704E-01;
    COFD[4627] = 1.17319185E-02;
    COFD[4628] = -2.24527123E+01;
    COFD[4629] = 5.35063033E+00;
    COFD[4630] = -4.34734059E-01;
    COFD[4631] = 1.70582710E-02;
    COFD[4632] = -2.09824046E+01;
    COFD[4633] = 4.66190324E+00;
    COFD[4634] = -3.22729985E-01;
    COFD[4635] = 1.13089121E-02;
    COFD[4636] = -2.20049811E+01;
    COFD[4637] = 5.03720506E+00;
    COFD[4638] = -3.82190496E-01;
    COFD[4639] = 1.43104209E-02;
    COFD[4640] = -2.20092308E+01;
    COFD[4641] = 5.03720506E+00;
    COFD[4642] = -3.82190496E-01;
    COFD[4643] = 1.43104209E-02;
    COFD[4644] = -2.14271873E+01;
    COFD[4645] = 4.71363498E+00;
    COFD[4646] = -3.30905003E-01;
    COFD[4647] = 1.17206344E-02;
    COFD[4648] = -2.25315464E+01;
    COFD[4649] = 5.35063033E+00;
    COFD[4650] = -4.34734059E-01;
    COFD[4651] = 1.70582710E-02;
    COFD[4652] = -2.25759756E+01;
    COFD[4653] = 5.35063033E+00;
    COFD[4654] = -4.34734059E-01;
    COFD[4655] = 1.70582710E-02;
    COFD[4656] = -2.25759756E+01;
    COFD[4657] = 5.35063033E+00;
    COFD[4658] = -4.34734059E-01;
    COFD[4659] = 1.70582710E-02;
    COFD[4660] = -2.25818241E+01;
    COFD[4661] = 5.35063033E+00;
    COFD[4662] = -4.34734059E-01;
    COFD[4663] = 1.70582710E-02;
    COFD[4664] = -2.25818241E+01;
    COFD[4665] = 5.35063033E+00;
    COFD[4666] = -4.34734059E-01;
    COFD[4667] = 1.70582710E-02;
    COFD[4668] = -2.26188244E+01;
    COFD[4669] = 5.35063033E+00;
    COFD[4670] = -4.34734059E-01;
    COFD[4671] = 1.70582710E-02;
    COFD[4672] = -2.26505977E+01;
    COFD[4673] = 5.35063033E+00;
    COFD[4674] = -4.34734059E-01;
    COFD[4675] = 1.70582710E-02;
    COFD[4676] = -2.02686523E+01;
    COFD[4677] = 5.01994886E+00;
    COFD[4678] = -4.31380254E-01;
    COFD[4679] = 1.84875367E-02;
    COFD[4680] = -2.00206731E+01;
    COFD[4681] = 5.52613176E+00;
    COFD[4682] = -4.85057412E-01;
    COFD[4683] = 2.03353154E-02;
    COFD[4684] = -1.62390700E+01;
    COFD[4685] = 4.04304103E+00;
    COFD[4686] = -3.12841265E-01;
    COFD[4687] = 1.36797409E-02;
    COFD[4688] = -2.16737420E+01;
    COFD[4689] = 5.52035712E+00;
    COFD[4690] = -4.84503119E-01;
    COFD[4691] = 2.03190484E-02;
    COFD[4692] = -2.16737420E+01;
    COFD[4693] = 5.52035712E+00;
    COFD[4694] = -4.84503119E-01;
    COFD[4695] = 2.03190484E-02;
    COFD[4696] = -2.17017719E+01;
    COFD[4697] = 5.52035712E+00;
    COFD[4698] = -4.84503119E-01;
    COFD[4699] = 2.03190484E-02;
    COFD[4700] = -1.96918321E+01;
    COFD[4701] = 4.95331445E+00;
    COFD[4702] = -4.23474055E-01;
    COFD[4703] = 1.81742301E-02;
    COFD[4704] = -2.16525572E+01;
    COFD[4705] = 5.50038516E+00;
    COFD[4706] = -4.82355440E-01;
    COFD[4707] = 2.02433624E-02;
    COFD[4708] = -1.97158821E+01;
    COFD[4709] = 4.95331445E+00;
    COFD[4710] = -4.23474055E-01;
    COFD[4711] = 1.81742301E-02;
    COFD[4712] = -1.86141216E+01;
    COFD[4713] = 3.75739003E+00;
    COFD[4714] = -1.83970343E-01;
    COFD[4715] = 4.46501245E-03;
    COFD[4716] = -2.22125735E+01;
    COFD[4717] = 5.42385157E+00;
    COFD[4718] = -4.47809271E-01;
    COFD[4719] = 1.77669962E-02;
    COFD[4720] = -2.22258107E+01;
    COFD[4721] = 5.42385157E+00;
    COFD[4722] = -4.47809271E-01;
    COFD[4723] = 1.77669962E-02;
    COFD[4724] = -2.07708766E+01;
    COFD[4725] = 5.16820334E+00;
    COFD[4726] = -4.47675828E-01;
    COFD[4727] = 1.90754011E-02;
    COFD[4728] = -2.22731582E+01;
    COFD[4729] = 5.50482046E+00;
    COFD[4730] = -4.63503431E-01;
    COFD[4731] = 1.86537631E-02;
    COFD[4732] = -1.99140853E+01;
    COFD[4733] = 4.18537863E+00;
    COFD[4734] = -2.48994776E-01;
    COFD[4735] = 7.64930102E-03;
    COFD[4736] = -2.24844780E+01;
    COFD[4737] = 5.49194178E+00;
    COFD[4738] = -4.60232233E-01;
    COFD[4739] = 1.84488745E-02;
    COFD[4740] = -1.99255944E+01;
    COFD[4741] = 4.18537863E+00;
    COFD[4742] = -2.48994776E-01;
    COFD[4743] = 7.64930102E-03;
    COFD[4744] = -2.24959645E+01;
    COFD[4745] = 5.49194178E+00;
    COFD[4746] = -4.60232233E-01;
    COFD[4747] = 1.84488745E-02;
    COFD[4748] = -2.08978461E+01;
    COFD[4749] = 4.65579398E+00;
    COFD[4750] = -3.21767942E-01;
    COFD[4751] = 1.12605667E-02;
    COFD[4752] = -2.08978461E+01;
    COFD[4753] = 4.65579398E+00;
    COFD[4754] = -3.21767942E-01;
    COFD[4755] = 1.12605667E-02;
    COFD[4756] = -2.11202704E+01;
    COFD[4757] = 5.28856610E+00;
    COFD[4758] = -4.61764591E-01;
    COFD[4759] = 1.96258820E-02;
    COFD[4760] = -2.11303909E+01;
    COFD[4761] = 5.28856610E+00;
    COFD[4762] = -4.61764591E-01;
    COFD[4763] = 1.96258820E-02;
    COFD[4764] = -2.11401024E+01;
    COFD[4765] = 5.28856610E+00;
    COFD[4766] = -4.61764591E-01;
    COFD[4767] = 1.96258820E-02;
    COFD[4768] = -2.24868912E+01;
    COFD[4769] = 5.49895896E+00;
    COFD[4770] = -4.61764691E-01;
    COFD[4771] = 1.85397456E-02;
    COFD[4772] = -2.08906722E+01;
    COFD[4773] = 4.56059926E+00;
    COFD[4774] = -3.06658953E-01;
    COFD[4775] = 1.04982481E-02;
    COFD[4776] = -2.21563247E+01;
    COFD[4777] = 5.15132481E+00;
    COFD[4778] = -4.00804233E-01;
    COFD[4779] = 1.52679625E-02;
    COFD[4780] = -2.05947670E+01;
    COFD[4781] = 4.35703694E+00;
    COFD[4782] = -2.75158780E-01;
    COFD[4783] = 8.93458153E-03;
    COFD[4784] = -2.21626591E+01;
    COFD[4785] = 5.15132481E+00;
    COFD[4786] = -4.00804233E-01;
    COFD[4787] = 1.52679625E-02;
    COFD[4788] = -2.02873532E+01;
    COFD[4789] = 4.28958705E+00;
    COFD[4790] = -2.64832562E-01;
    COFD[4791] = 8.42593887E-03;
    COFD[4792] = -2.14229024E+01;
    COFD[4793] = 4.71363498E+00;
    COFD[4794] = -3.30905003E-01;
    COFD[4795] = 1.17206344E-02;
    COFD[4796] = -2.14271873E+01;
    COFD[4797] = 4.71363498E+00;
    COFD[4798] = -3.30905003E-01;
    COFD[4799] = 1.17206344E-02;
    COFD[4800] = -2.07619819E+01;
    COFD[4801] = 4.35524100E+00;
    COFD[4802] = -2.74882904E-01;
    COFD[4803] = 8.92095745E-03;
    COFD[4804] = -2.22421099E+01;
    COFD[4805] = 5.15132481E+00;
    COFD[4806] = -4.00804233E-01;
    COFD[4807] = 1.52679625E-02;
    COFD[4808] = -2.22869324E+01;
    COFD[4809] = 5.15132481E+00;
    COFD[4810] = -4.00804233E-01;
    COFD[4811] = 1.52679625E-02;
    COFD[4812] = -2.22869324E+01;
    COFD[4813] = 5.15132481E+00;
    COFD[4814] = -4.00804233E-01;
    COFD[4815] = 1.52679625E-02;
    COFD[4816] = -2.22928354E+01;
    COFD[4817] = 5.15132481E+00;
    COFD[4818] = -4.00804233E-01;
    COFD[4819] = 1.52679625E-02;
    COFD[4820] = -2.22928354E+01;
    COFD[4821] = 5.15132481E+00;
    COFD[4822] = -4.00804233E-01;
    COFD[4823] = 1.52679625E-02;
    COFD[4824] = -2.23301955E+01;
    COFD[4825] = 5.15132481E+00;
    COFD[4826] = -4.00804233E-01;
    COFD[4827] = 1.52679625E-02;
    COFD[4828] = -2.23622999E+01;
    COFD[4829] = 5.15132481E+00;
    COFD[4830] = -4.00804233E-01;
    COFD[4831] = 1.52679625E-02;
    COFD[4832] = -2.07441338E+01;
    COFD[4833] = 5.16159824E+00;
    COFD[4834] = -4.46935718E-01;
    COFD[4835] = 1.90480458E-02;
    COFD[4836] = -1.89252713E+01;
    COFD[4837] = 5.20981862E+00;
    COFD[4838] = -4.52489825E-01;
    COFD[4839] = 1.92609226E-02;
    COFD[4840] = -1.49129599E+01;
    COFD[4841] = 3.56974825E+00;
    COFD[4842] = -2.52221138E-01;
    COFD[4843] = 1.10819767E-02;
    COFD[4844] = -2.05875936E+01;
    COFD[4845] = 5.20087227E+00;
    COFD[4846] = -4.51444972E-01;
    COFD[4847] = 1.92201898E-02;
    COFD[4848] = -2.05875936E+01;
    COFD[4849] = 5.20087227E+00;
    COFD[4850] = -4.51444972E-01;
    COFD[4851] = 1.92201898E-02;
    COFD[4852] = -2.06157113E+01;
    COFD[4853] = 5.20087227E+00;
    COFD[4854] = -4.51444972E-01;
    COFD[4855] = 1.92201898E-02;
    COFD[4856] = -1.83410870E+01;
    COFD[4857] = 4.49977587E+00;
    COFD[4858] = -3.68989022E-01;
    COFD[4859] = 1.59879891E-02;
    COFD[4860] = -2.05609682E+01;
    COFD[4861] = 5.17771473E+00;
    COFD[4862] = -4.48750201E-01;
    COFD[4863] = 1.91155567E-02;
    COFD[4864] = -1.83652204E+01;
    COFD[4865] = 4.49977587E+00;
    COFD[4866] = -3.68989022E-01;
    COFD[4867] = 1.59879891E-02;
    COFD[4868] = -2.05420607E+01;
    COFD[4869] = 4.74677479E+00;
    COFD[4870] = -3.36160335E-01;
    COFD[4871] = 1.19858600E-02;
    COFD[4872] = -2.22161430E+01;
    COFD[4873] = 5.59047891E+00;
    COFD[4874] = -4.85359237E-01;
    COFD[4875] = 2.00302091E-02;
    COFD[4876] = -2.22294456E+01;
    COFD[4877] = 5.59047891E+00;
    COFD[4878] = -4.85359237E-01;
    COFD[4879] = 2.00302091E-02;
    COFD[4880] = -1.95177006E+01;
    COFD[4881] = 4.76145602E+00;
    COFD[4882] = -4.00684587E-01;
    COFD[4883] = 1.72708322E-02;
    COFD[4884] = -2.21308159E+01;
    COFD[4885] = 5.60509076E+00;
    COFD[4886] = -4.91228646E-01;
    COFD[4887] = 2.04425499E-02;
    COFD[4888] = -2.15522331E+01;
    COFD[4889] = 5.05104142E+00;
    COFD[4890] = -3.84428177E-01;
    COFD[4891] = 1.44249285E-02;
    COFD[4892] = -2.23901271E+01;
    COFD[4893] = 5.61213963E+00;
    COFD[4894] = -4.90953451E-01;
    COFD[4895] = 2.03841518E-02;
    COFD[4896] = -2.15638034E+01;
    COFD[4897] = 5.05104142E+00;
    COFD[4898] = -3.84428177E-01;
    COFD[4899] = 1.44249285E-02;
    COFD[4900] = -2.24016748E+01;
    COFD[4901] = 5.61213963E+00;
    COFD[4902] = -4.90953451E-01;
    COFD[4903] = 2.03841518E-02;
    COFD[4904] = -2.20544777E+01;
    COFD[4905] = 5.31843386E+00;
    COFD[4906] = -4.29198961E-01;
    COFD[4907] = 1.67645078E-02;
    COFD[4908] = -2.20544777E+01;
    COFD[4909] = 5.31843386E+00;
    COFD[4910] = -4.29198961E-01;
    COFD[4911] = 1.67645078E-02;
    COFD[4912] = -1.98058079E+01;
    COFD[4913] = 4.86038833E+00;
    COFD[4914] = -4.12187543E-01;
    COFD[4915] = 1.77155642E-02;
    COFD[4916] = -1.98159859E+01;
    COFD[4917] = 4.86038833E+00;
    COFD[4918] = -4.12187543E-01;
    COFD[4919] = 1.77155642E-02;
    COFD[4920] = -1.98257536E+01;
    COFD[4921] = 4.86038833E+00;
    COFD[4922] = -4.12187543E-01;
    COFD[4923] = 1.77155642E-02;
    COFD[4924] = -2.23744741E+01;
    COFD[4925] = 5.61208013E+00;
    COFD[4926] = -4.91433206E-01;
    COFD[4927] = 2.04241259E-02;
    COFD[4928] = -2.21263189E+01;
    COFD[4929] = 5.25559319E+00;
    COFD[4930] = -4.18539684E-01;
    COFD[4931] = 1.62026965E-02;
    COFD[4932] = -2.26488330E+01;
    COFD[4933] = 5.52719085E+00;
    COFD[4934] = -4.68021516E-01;
    COFD[4935] = 1.89143588E-02;
    COFD[4936] = -2.20698134E+01;
    COFD[4937] = 5.15225414E+00;
    COFD[4938] = -4.00960404E-01;
    COFD[4939] = 1.52761358E-02;
    COFD[4940] = -2.26552117E+01;
    COFD[4941] = 5.52719085E+00;
    COFD[4942] = -4.68021516E-01;
    COFD[4943] = 1.89143588E-02;
    COFD[4944] = -2.18257833E+01;
    COFD[4945] = 5.11317705E+00;
    COFD[4946] = -3.94506858E-01;
    COFD[4947] = 1.49416998E-02;
    COFD[4948] = -2.25272268E+01;
    COFD[4949] = 5.35063033E+00;
    COFD[4950] = -4.34734059E-01;
    COFD[4951] = 1.70582710E-02;
    COFD[4952] = -2.25315464E+01;
    COFD[4953] = 5.35063033E+00;
    COFD[4954] = -4.34734059E-01;
    COFD[4955] = 1.70582710E-02;
    COFD[4956] = -2.22421099E+01;
    COFD[4957] = 5.15132481E+00;
    COFD[4958] = -4.00804233E-01;
    COFD[4959] = 1.52679625E-02;
    COFD[4960] = -2.27352685E+01;
    COFD[4961] = 5.52719085E+00;
    COFD[4962] = -4.68021516E-01;
    COFD[4963] = 1.89143588E-02;
    COFD[4964] = -2.27804782E+01;
    COFD[4965] = 5.52719085E+00;
    COFD[4966] = -4.68021516E-01;
    COFD[4967] = 1.89143588E-02;
    COFD[4968] = -2.27804782E+01;
    COFD[4969] = 5.52719085E+00;
    COFD[4970] = -4.68021516E-01;
    COFD[4971] = 1.89143588E-02;
    COFD[4972] = -2.27864348E+01;
    COFD[4973] = 5.52719085E+00;
    COFD[4974] = -4.68021516E-01;
    COFD[4975] = 1.89143588E-02;
    COFD[4976] = -2.27864348E+01;
    COFD[4977] = 5.52719085E+00;
    COFD[4978] = -4.68021516E-01;
    COFD[4979] = 1.89143588E-02;
    COFD[4980] = -2.28241497E+01;
    COFD[4981] = 5.52719085E+00;
    COFD[4982] = -4.68021516E-01;
    COFD[4983] = 1.89143588E-02;
    COFD[4984] = -2.28565811E+01;
    COFD[4985] = 5.52719085E+00;
    COFD[4986] = -4.68021516E-01;
    COFD[4987] = 1.89143588E-02;
    COFD[4988] = -1.94899278E+01;
    COFD[4989] = 4.75414932E+00;
    COFD[4990] = -3.99807084E-01;
    COFD[4991] = 1.72356137E-02;
    COFD[4992] = -1.89266547E+01;
    COFD[4993] = 5.20981862E+00;
    COFD[4994] = -4.52489825E-01;
    COFD[4995] = 1.92609226E-02;
    COFD[4996] = -1.49156868E+01;
    COFD[4997] = 3.56974825E+00;
    COFD[4998] = -2.52221138E-01;
    COFD[4999] = 1.10819767E-02;
    COFD[5000] = -2.06037892E+01;
    COFD[5001] = 5.20087227E+00;
    COFD[5002] = -4.51444972E-01;
    COFD[5003] = 1.92201898E-02;
    COFD[5004] = -2.06037892E+01;
    COFD[5005] = 5.20087227E+00;
    COFD[5006] = -4.51444972E-01;
    COFD[5007] = 1.92201898E-02;
    COFD[5008] = -2.06328599E+01;
    COFD[5009] = 5.20087227E+00;
    COFD[5010] = -4.51444972E-01;
    COFD[5011] = 1.92201898E-02;
    COFD[5012] = -1.83591261E+01;
    COFD[5013] = 4.49977587E+00;
    COFD[5014] = -3.68989022E-01;
    COFD[5015] = 1.59879891E-02;
    COFD[5016] = -2.05790471E+01;
    COFD[5017] = 5.17771473E+00;
    COFD[5018] = -4.48750201E-01;
    COFD[5019] = 1.91155567E-02;
    COFD[5020] = -1.83841687E+01;
    COFD[5021] = 4.49977587E+00;
    COFD[5022] = -3.68989022E-01;
    COFD[5023] = 1.59879891E-02;
    COFD[5024] = -2.05618968E+01;
    COFD[5025] = 4.74677479E+00;
    COFD[5026] = -3.36160335E-01;
    COFD[5027] = 1.19858600E-02;
    COFD[5028] = -2.22423680E+01;
    COFD[5029] = 5.59047891E+00;
    COFD[5030] = -4.85359237E-01;
    COFD[5031] = 2.00302091E-02;
    COFD[5032] = -2.22563971E+01;
    COFD[5033] = 5.59047891E+00;
    COFD[5034] = -4.85359237E-01;
    COFD[5035] = 2.00302091E-02;
    COFD[5036] = -1.95453329E+01;
    COFD[5037] = 4.76145602E+00;
    COFD[5038] = -4.00684587E-01;
    COFD[5039] = 1.72708322E-02;
    COFD[5040] = -2.21584786E+01;
    COFD[5041] = 5.60509076E+00;
    COFD[5042] = -4.91228646E-01;
    COFD[5043] = 2.04425499E-02;
    COFD[5044] = -2.15805625E+01;
    COFD[5045] = 5.05104142E+00;
    COFD[5046] = -3.84428177E-01;
    COFD[5047] = 1.44249285E-02;
    COFD[5048] = -2.24184864E+01;
    COFD[5049] = 5.61213963E+00;
    COFD[5050] = -4.90953451E-01;
    COFD[5051] = 2.03841518E-02;
    COFD[5052] = -2.15928156E+01;
    COFD[5053] = 5.05104142E+00;
    COFD[5054] = -3.84428177E-01;
    COFD[5055] = 1.44249285E-02;
    COFD[5056] = -2.24307163E+01;
    COFD[5057] = 5.61213963E+00;
    COFD[5058] = -4.90953451E-01;
    COFD[5059] = 2.03841518E-02;
    COFD[5060] = -2.20841588E+01;
    COFD[5061] = 5.31843386E+00;
    COFD[5062] = -4.29198961E-01;
    COFD[5063] = 1.67645078E-02;
    COFD[5064] = -2.20841588E+01;
    COFD[5065] = 5.31843386E+00;
    COFD[5066] = -4.29198961E-01;
    COFD[5067] = 1.67645078E-02;
    COFD[5068] = -1.98361165E+01;
    COFD[5069] = 4.86038833E+00;
    COFD[5070] = -4.12187543E-01;
    COFD[5071] = 1.77155642E-02;
    COFD[5072] = -1.98469374E+01;
    COFD[5073] = 4.86038833E+00;
    COFD[5074] = -4.12187543E-01;
    COFD[5075] = 1.77155642E-02;
    COFD[5076] = -1.98573353E+01;
    COFD[5077] = 4.86038833E+00;
    COFD[5078] = -4.12187543E-01;
    COFD[5079] = 1.77155642E-02;
    COFD[5080] = -2.24116928E+01;
    COFD[5081] = 5.61208013E+00;
    COFD[5082] = -4.91433206E-01;
    COFD[5083] = 2.04241259E-02;
    COFD[5084] = -2.21635600E+01;
    COFD[5085] = 5.25559319E+00;
    COFD[5086] = -4.18539684E-01;
    COFD[5087] = 1.62026965E-02;
    COFD[5088] = -2.26865869E+01;
    COFD[5089] = 5.52719085E+00;
    COFD[5090] = -4.68021516E-01;
    COFD[5091] = 1.89143588E-02;
    COFD[5092] = -2.21080494E+01;
    COFD[5093] = 5.15225414E+00;
    COFD[5094] = -4.00960404E-01;
    COFD[5095] = 1.52761358E-02;
    COFD[5096] = -2.26934694E+01;
    COFD[5097] = 5.52719085E+00;
    COFD[5098] = -4.68021516E-01;
    COFD[5099] = 1.89143588E-02;
    COFD[5100] = -2.18645147E+01;
    COFD[5101] = 5.11317705E+00;
    COFD[5102] = -3.94506858E-01;
    COFD[5103] = 1.49416998E-02;
    COFD[5104] = -2.25712566E+01;
    COFD[5105] = 5.35063033E+00;
    COFD[5106] = -4.34734059E-01;
    COFD[5107] = 1.70582710E-02;
    COFD[5108] = -2.25759756E+01;
    COFD[5109] = 5.35063033E+00;
    COFD[5110] = -4.34734059E-01;
    COFD[5111] = 1.70582710E-02;
    COFD[5112] = -2.22869324E+01;
    COFD[5113] = 5.15132481E+00;
    COFD[5114] = -4.00804233E-01;
    COFD[5115] = 1.52679625E-02;
    COFD[5116] = -2.27804782E+01;
    COFD[5117] = 5.52719085E+00;
    COFD[5118] = -4.68021516E-01;
    COFD[5119] = 1.89143588E-02;
    COFD[5120] = -2.28301858E+01;
    COFD[5121] = 5.52719085E+00;
    COFD[5122] = -4.68021516E-01;
    COFD[5123] = 1.89143588E-02;
    COFD[5124] = -2.28301858E+01;
    COFD[5125] = 5.52719085E+00;
    COFD[5126] = -4.68021516E-01;
    COFD[5127] = 1.89143588E-02;
    COFD[5128] = -2.28367691E+01;
    COFD[5129] = 5.52719085E+00;
    COFD[5130] = -4.68021516E-01;
    COFD[5131] = 1.89143588E-02;
    COFD[5132] = -2.28367691E+01;
    COFD[5133] = 5.52719085E+00;
    COFD[5134] = -4.68021516E-01;
    COFD[5135] = 1.89143588E-02;
    COFD[5136] = -2.28786502E+01;
    COFD[5137] = 5.52719085E+00;
    COFD[5138] = -4.68021516E-01;
    COFD[5139] = 1.89143588E-02;
    COFD[5140] = -2.29149554E+01;
    COFD[5141] = 5.52719085E+00;
    COFD[5142] = -4.68021516E-01;
    COFD[5143] = 1.89143588E-02;
    COFD[5144] = -1.95175620E+01;
    COFD[5145] = 4.75414932E+00;
    COFD[5146] = -3.99807084E-01;
    COFD[5147] = 1.72356137E-02;
    COFD[5148] = -1.89266547E+01;
    COFD[5149] = 5.20981862E+00;
    COFD[5150] = -4.52489825E-01;
    COFD[5151] = 1.92609226E-02;
    COFD[5152] = -1.49156868E+01;
    COFD[5153] = 3.56974825E+00;
    COFD[5154] = -2.52221138E-01;
    COFD[5155] = 1.10819767E-02;
    COFD[5156] = -2.06037892E+01;
    COFD[5157] = 5.20087227E+00;
    COFD[5158] = -4.51444972E-01;
    COFD[5159] = 1.92201898E-02;
    COFD[5160] = -2.06037892E+01;
    COFD[5161] = 5.20087227E+00;
    COFD[5162] = -4.51444972E-01;
    COFD[5163] = 1.92201898E-02;
    COFD[5164] = -2.06328599E+01;
    COFD[5165] = 5.20087227E+00;
    COFD[5166] = -4.51444972E-01;
    COFD[5167] = 1.92201898E-02;
    COFD[5168] = -1.83591261E+01;
    COFD[5169] = 4.49977587E+00;
    COFD[5170] = -3.68989022E-01;
    COFD[5171] = 1.59879891E-02;
    COFD[5172] = -2.05790471E+01;
    COFD[5173] = 5.17771473E+00;
    COFD[5174] = -4.48750201E-01;
    COFD[5175] = 1.91155567E-02;
    COFD[5176] = -1.83841687E+01;
    COFD[5177] = 4.49977587E+00;
    COFD[5178] = -3.68989022E-01;
    COFD[5179] = 1.59879891E-02;
    COFD[5180] = -2.05618968E+01;
    COFD[5181] = 4.74677479E+00;
    COFD[5182] = -3.36160335E-01;
    COFD[5183] = 1.19858600E-02;
    COFD[5184] = -2.22423680E+01;
    COFD[5185] = 5.59047891E+00;
    COFD[5186] = -4.85359237E-01;
    COFD[5187] = 2.00302091E-02;
    COFD[5188] = -2.22563971E+01;
    COFD[5189] = 5.59047891E+00;
    COFD[5190] = -4.85359237E-01;
    COFD[5191] = 2.00302091E-02;
    COFD[5192] = -1.95453329E+01;
    COFD[5193] = 4.76145602E+00;
    COFD[5194] = -4.00684587E-01;
    COFD[5195] = 1.72708322E-02;
    COFD[5196] = -2.21584786E+01;
    COFD[5197] = 5.60509076E+00;
    COFD[5198] = -4.91228646E-01;
    COFD[5199] = 2.04425499E-02;
    COFD[5200] = -2.15805625E+01;
    COFD[5201] = 5.05104142E+00;
    COFD[5202] = -3.84428177E-01;
    COFD[5203] = 1.44249285E-02;
    COFD[5204] = -2.24184864E+01;
    COFD[5205] = 5.61213963E+00;
    COFD[5206] = -4.90953451E-01;
    COFD[5207] = 2.03841518E-02;
    COFD[5208] = -2.15928156E+01;
    COFD[5209] = 5.05104142E+00;
    COFD[5210] = -3.84428177E-01;
    COFD[5211] = 1.44249285E-02;
    COFD[5212] = -2.24307163E+01;
    COFD[5213] = 5.61213963E+00;
    COFD[5214] = -4.90953451E-01;
    COFD[5215] = 2.03841518E-02;
    COFD[5216] = -2.20841588E+01;
    COFD[5217] = 5.31843386E+00;
    COFD[5218] = -4.29198961E-01;
    COFD[5219] = 1.67645078E-02;
    COFD[5220] = -2.20841588E+01;
    COFD[5221] = 5.31843386E+00;
    COFD[5222] = -4.29198961E-01;
    COFD[5223] = 1.67645078E-02;
    COFD[5224] = -1.98361165E+01;
    COFD[5225] = 4.86038833E+00;
    COFD[5226] = -4.12187543E-01;
    COFD[5227] = 1.77155642E-02;
    COFD[5228] = -1.98469374E+01;
    COFD[5229] = 4.86038833E+00;
    COFD[5230] = -4.12187543E-01;
    COFD[5231] = 1.77155642E-02;
    COFD[5232] = -1.98573353E+01;
    COFD[5233] = 4.86038833E+00;
    COFD[5234] = -4.12187543E-01;
    COFD[5235] = 1.77155642E-02;
    COFD[5236] = -2.24116928E+01;
    COFD[5237] = 5.61208013E+00;
    COFD[5238] = -4.91433206E-01;
    COFD[5239] = 2.04241259E-02;
    COFD[5240] = -2.21635600E+01;
    COFD[5241] = 5.25559319E+00;
    COFD[5242] = -4.18539684E-01;
    COFD[5243] = 1.62026965E-02;
    COFD[5244] = -2.26865869E+01;
    COFD[5245] = 5.52719085E+00;
    COFD[5246] = -4.68021516E-01;
    COFD[5247] = 1.89143588E-02;
    COFD[5248] = -2.21080494E+01;
    COFD[5249] = 5.15225414E+00;
    COFD[5250] = -4.00960404E-01;
    COFD[5251] = 1.52761358E-02;
    COFD[5252] = -2.26934694E+01;
    COFD[5253] = 5.52719085E+00;
    COFD[5254] = -4.68021516E-01;
    COFD[5255] = 1.89143588E-02;
    COFD[5256] = -2.18645147E+01;
    COFD[5257] = 5.11317705E+00;
    COFD[5258] = -3.94506858E-01;
    COFD[5259] = 1.49416998E-02;
    COFD[5260] = -2.25712566E+01;
    COFD[5261] = 5.35063033E+00;
    COFD[5262] = -4.34734059E-01;
    COFD[5263] = 1.70582710E-02;
    COFD[5264] = -2.25759756E+01;
    COFD[5265] = 5.35063033E+00;
    COFD[5266] = -4.34734059E-01;
    COFD[5267] = 1.70582710E-02;
    COFD[5268] = -2.22869324E+01;
    COFD[5269] = 5.15132481E+00;
    COFD[5270] = -4.00804233E-01;
    COFD[5271] = 1.52679625E-02;
    COFD[5272] = -2.27804782E+01;
    COFD[5273] = 5.52719085E+00;
    COFD[5274] = -4.68021516E-01;
    COFD[5275] = 1.89143588E-02;
    COFD[5276] = -2.28301858E+01;
    COFD[5277] = 5.52719085E+00;
    COFD[5278] = -4.68021516E-01;
    COFD[5279] = 1.89143588E-02;
    COFD[5280] = -2.28301858E+01;
    COFD[5281] = 5.52719085E+00;
    COFD[5282] = -4.68021516E-01;
    COFD[5283] = 1.89143588E-02;
    COFD[5284] = -2.28367691E+01;
    COFD[5285] = 5.52719085E+00;
    COFD[5286] = -4.68021516E-01;
    COFD[5287] = 1.89143588E-02;
    COFD[5288] = -2.28367691E+01;
    COFD[5289] = 5.52719085E+00;
    COFD[5290] = -4.68021516E-01;
    COFD[5291] = 1.89143588E-02;
    COFD[5292] = -2.28786502E+01;
    COFD[5293] = 5.52719085E+00;
    COFD[5294] = -4.68021516E-01;
    COFD[5295] = 1.89143588E-02;
    COFD[5296] = -2.29149554E+01;
    COFD[5297] = 5.52719085E+00;
    COFD[5298] = -4.68021516E-01;
    COFD[5299] = 1.89143588E-02;
    COFD[5300] = -1.95175620E+01;
    COFD[5301] = 4.75414932E+00;
    COFD[5302] = -3.99807084E-01;
    COFD[5303] = 1.72356137E-02;
    COFD[5304] = -1.89268281E+01;
    COFD[5305] = 5.20981862E+00;
    COFD[5306] = -4.52489825E-01;
    COFD[5307] = 1.92609226E-02;
    COFD[5308] = -1.49160291E+01;
    COFD[5309] = 3.56974825E+00;
    COFD[5310] = -2.52221138E-01;
    COFD[5311] = 1.10819767E-02;
    COFD[5312] = -2.06058533E+01;
    COFD[5313] = 5.20087227E+00;
    COFD[5314] = -4.51444972E-01;
    COFD[5315] = 1.92201898E-02;
    COFD[5316] = -2.06058533E+01;
    COFD[5317] = 5.20087227E+00;
    COFD[5318] = -4.51444972E-01;
    COFD[5319] = 1.92201898E-02;
    COFD[5320] = -2.06350479E+01;
    COFD[5321] = 5.20087227E+00;
    COFD[5322] = -4.51444972E-01;
    COFD[5323] = 1.92201898E-02;
    COFD[5324] = -1.83614301E+01;
    COFD[5325] = 4.49977587E+00;
    COFD[5326] = -3.68989022E-01;
    COFD[5327] = 1.59879891E-02;
    COFD[5328] = -2.05813563E+01;
    COFD[5329] = 5.17771473E+00;
    COFD[5330] = -4.48750201E-01;
    COFD[5331] = 1.91155567E-02;
    COFD[5332] = -1.83865912E+01;
    COFD[5333] = 4.49977587E+00;
    COFD[5334] = -3.68989022E-01;
    COFD[5335] = 1.59879891E-02;
    COFD[5336] = -2.05644355E+01;
    COFD[5337] = 4.74677479E+00;
    COFD[5338] = -3.36160335E-01;
    COFD[5339] = 1.19858600E-02;
    COFD[5340] = -2.22457488E+01;
    COFD[5341] = 5.59047891E+00;
    COFD[5342] = -4.85359237E-01;
    COFD[5343] = 2.00302091E-02;
    COFD[5344] = -2.22598745E+01;
    COFD[5345] = 5.59047891E+00;
    COFD[5346] = -4.85359237E-01;
    COFD[5347] = 2.00302091E-02;
    COFD[5348] = -1.95489008E+01;
    COFD[5349] = 4.76145602E+00;
    COFD[5350] = -4.00684587E-01;
    COFD[5351] = 1.72708322E-02;
    COFD[5352] = -2.21620506E+01;
    COFD[5353] = 5.60509076E+00;
    COFD[5354] = -4.91228646E-01;
    COFD[5355] = 2.04425499E-02;
    COFD[5356] = -2.15842234E+01;
    COFD[5357] = 5.05104142E+00;
    COFD[5358] = -3.84428177E-01;
    COFD[5359] = 1.44249285E-02;
    COFD[5360] = -2.24221512E+01;
    COFD[5361] = 5.61213963E+00;
    COFD[5362] = -4.90953451E-01;
    COFD[5363] = 2.03841518E-02;
    COFD[5364] = -2.15965677E+01;
    COFD[5365] = 5.05104142E+00;
    COFD[5366] = -3.84428177E-01;
    COFD[5367] = 1.44249285E-02;
    COFD[5368] = -2.24344722E+01;
    COFD[5369] = 5.61213963E+00;
    COFD[5370] = -4.90953451E-01;
    COFD[5371] = 2.03841518E-02;
    COFD[5372] = -2.20880003E+01;
    COFD[5373] = 5.31843386E+00;
    COFD[5374] = -4.29198961E-01;
    COFD[5375] = 1.67645078E-02;
    COFD[5376] = -2.20880003E+01;
    COFD[5377] = 5.31843386E+00;
    COFD[5378] = -4.29198961E-01;
    COFD[5379] = 1.67645078E-02;
    COFD[5380] = -1.98400420E+01;
    COFD[5381] = 4.86038833E+00;
    COFD[5382] = -4.12187543E-01;
    COFD[5383] = 1.77155642E-02;
    COFD[5384] = -1.98509491E+01;
    COFD[5385] = 4.86038833E+00;
    COFD[5386] = -4.12187543E-01;
    COFD[5387] = 1.77155642E-02;
    COFD[5388] = -1.98614317E+01;
    COFD[5389] = 4.86038833E+00;
    COFD[5390] = -4.12187543E-01;
    COFD[5391] = 1.77155642E-02;
    COFD[5392] = -2.24165516E+01;
    COFD[5393] = 5.61208013E+00;
    COFD[5394] = -4.91433206E-01;
    COFD[5395] = 2.04241259E-02;
    COFD[5396] = -2.21684219E+01;
    COFD[5397] = 5.25559319E+00;
    COFD[5398] = -4.18539684E-01;
    COFD[5399] = 1.62026965E-02;
    COFD[5400] = -2.26915186E+01;
    COFD[5401] = 5.52719085E+00;
    COFD[5402] = -4.68021516E-01;
    COFD[5403] = 1.89143588E-02;
    COFD[5404] = -2.21130468E+01;
    COFD[5405] = 5.15225414E+00;
    COFD[5406] = -4.00960404E-01;
    COFD[5407] = 1.52761358E-02;
    COFD[5408] = -2.26984698E+01;
    COFD[5409] = 5.52719085E+00;
    COFD[5410] = -4.68021516E-01;
    COFD[5411] = 1.89143588E-02;
    COFD[5412] = -2.18695798E+01;
    COFD[5413] = 5.11317705E+00;
    COFD[5414] = -3.94506858E-01;
    COFD[5415] = 1.49416998E-02;
    COFD[5416] = -2.25770498E+01;
    COFD[5417] = 5.35063033E+00;
    COFD[5418] = -4.34734059E-01;
    COFD[5419] = 1.70582710E-02;
    COFD[5420] = -2.25818241E+01;
    COFD[5421] = 5.35063033E+00;
    COFD[5422] = -4.34734059E-01;
    COFD[5423] = 1.70582710E-02;
    COFD[5424] = -2.22928354E+01;
    COFD[5425] = 5.15132481E+00;
    COFD[5426] = -4.00804233E-01;
    COFD[5427] = 1.52679625E-02;
    COFD[5428] = -2.27864348E+01;
    COFD[5429] = 5.52719085E+00;
    COFD[5430] = -4.68021516E-01;
    COFD[5431] = 1.89143588E-02;
    COFD[5432] = -2.28367691E+01;
    COFD[5433] = 5.52719085E+00;
    COFD[5434] = -4.68021516E-01;
    COFD[5435] = 1.89143588E-02;
    COFD[5436] = -2.28367691E+01;
    COFD[5437] = 5.52719085E+00;
    COFD[5438] = -4.68021516E-01;
    COFD[5439] = 1.89143588E-02;
    COFD[5440] = -2.28434402E+01;
    COFD[5441] = 5.52719085E+00;
    COFD[5442] = -4.68021516E-01;
    COFD[5443] = 1.89143588E-02;
    COFD[5444] = -2.28434402E+01;
    COFD[5445] = 5.52719085E+00;
    COFD[5446] = -4.68021516E-01;
    COFD[5447] = 1.89143588E-02;
    COFD[5448] = -2.28859084E+01;
    COFD[5449] = 5.52719085E+00;
    COFD[5450] = -4.68021516E-01;
    COFD[5451] = 1.89143588E-02;
    COFD[5452] = -2.29227645E+01;
    COFD[5453] = 5.52719085E+00;
    COFD[5454] = -4.68021516E-01;
    COFD[5455] = 1.89143588E-02;
    COFD[5456] = -1.95211303E+01;
    COFD[5457] = 4.75414932E+00;
    COFD[5458] = -3.99807084E-01;
    COFD[5459] = 1.72356137E-02;
    COFD[5460] = -1.89268281E+01;
    COFD[5461] = 5.20981862E+00;
    COFD[5462] = -4.52489825E-01;
    COFD[5463] = 1.92609226E-02;
    COFD[5464] = -1.49160291E+01;
    COFD[5465] = 3.56974825E+00;
    COFD[5466] = -2.52221138E-01;
    COFD[5467] = 1.10819767E-02;
    COFD[5468] = -2.06058533E+01;
    COFD[5469] = 5.20087227E+00;
    COFD[5470] = -4.51444972E-01;
    COFD[5471] = 1.92201898E-02;
    COFD[5472] = -2.06058533E+01;
    COFD[5473] = 5.20087227E+00;
    COFD[5474] = -4.51444972E-01;
    COFD[5475] = 1.92201898E-02;
    COFD[5476] = -2.06350479E+01;
    COFD[5477] = 5.20087227E+00;
    COFD[5478] = -4.51444972E-01;
    COFD[5479] = 1.92201898E-02;
    COFD[5480] = -1.83614301E+01;
    COFD[5481] = 4.49977587E+00;
    COFD[5482] = -3.68989022E-01;
    COFD[5483] = 1.59879891E-02;
    COFD[5484] = -2.05813563E+01;
    COFD[5485] = 5.17771473E+00;
    COFD[5486] = -4.48750201E-01;
    COFD[5487] = 1.91155567E-02;
    COFD[5488] = -1.83865912E+01;
    COFD[5489] = 4.49977587E+00;
    COFD[5490] = -3.68989022E-01;
    COFD[5491] = 1.59879891E-02;
    COFD[5492] = -2.05644355E+01;
    COFD[5493] = 4.74677479E+00;
    COFD[5494] = -3.36160335E-01;
    COFD[5495] = 1.19858600E-02;
    COFD[5496] = -2.22457488E+01;
    COFD[5497] = 5.59047891E+00;
    COFD[5498] = -4.85359237E-01;
    COFD[5499] = 2.00302091E-02;
    COFD[5500] = -2.22598745E+01;
    COFD[5501] = 5.59047891E+00;
    COFD[5502] = -4.85359237E-01;
    COFD[5503] = 2.00302091E-02;
    COFD[5504] = -1.95489008E+01;
    COFD[5505] = 4.76145602E+00;
    COFD[5506] = -4.00684587E-01;
    COFD[5507] = 1.72708322E-02;
    COFD[5508] = -2.21620506E+01;
    COFD[5509] = 5.60509076E+00;
    COFD[5510] = -4.91228646E-01;
    COFD[5511] = 2.04425499E-02;
    COFD[5512] = -2.15842234E+01;
    COFD[5513] = 5.05104142E+00;
    COFD[5514] = -3.84428177E-01;
    COFD[5515] = 1.44249285E-02;
    COFD[5516] = -2.24221512E+01;
    COFD[5517] = 5.61213963E+00;
    COFD[5518] = -4.90953451E-01;
    COFD[5519] = 2.03841518E-02;
    COFD[5520] = -2.15965677E+01;
    COFD[5521] = 5.05104142E+00;
    COFD[5522] = -3.84428177E-01;
    COFD[5523] = 1.44249285E-02;
    COFD[5524] = -2.24344722E+01;
    COFD[5525] = 5.61213963E+00;
    COFD[5526] = -4.90953451E-01;
    COFD[5527] = 2.03841518E-02;
    COFD[5528] = -2.20880003E+01;
    COFD[5529] = 5.31843386E+00;
    COFD[5530] = -4.29198961E-01;
    COFD[5531] = 1.67645078E-02;
    COFD[5532] = -2.20880003E+01;
    COFD[5533] = 5.31843386E+00;
    COFD[5534] = -4.29198961E-01;
    COFD[5535] = 1.67645078E-02;
    COFD[5536] = -1.98400420E+01;
    COFD[5537] = 4.86038833E+00;
    COFD[5538] = -4.12187543E-01;
    COFD[5539] = 1.77155642E-02;
    COFD[5540] = -1.98509491E+01;
    COFD[5541] = 4.86038833E+00;
    COFD[5542] = -4.12187543E-01;
    COFD[5543] = 1.77155642E-02;
    COFD[5544] = -1.98614317E+01;
    COFD[5545] = 4.86038833E+00;
    COFD[5546] = -4.12187543E-01;
    COFD[5547] = 1.77155642E-02;
    COFD[5548] = -2.24165516E+01;
    COFD[5549] = 5.61208013E+00;
    COFD[5550] = -4.91433206E-01;
    COFD[5551] = 2.04241259E-02;
    COFD[5552] = -2.21684219E+01;
    COFD[5553] = 5.25559319E+00;
    COFD[5554] = -4.18539684E-01;
    COFD[5555] = 1.62026965E-02;
    COFD[5556] = -2.26915186E+01;
    COFD[5557] = 5.52719085E+00;
    COFD[5558] = -4.68021516E-01;
    COFD[5559] = 1.89143588E-02;
    COFD[5560] = -2.21130468E+01;
    COFD[5561] = 5.15225414E+00;
    COFD[5562] = -4.00960404E-01;
    COFD[5563] = 1.52761358E-02;
    COFD[5564] = -2.26984698E+01;
    COFD[5565] = 5.52719085E+00;
    COFD[5566] = -4.68021516E-01;
    COFD[5567] = 1.89143588E-02;
    COFD[5568] = -2.18695798E+01;
    COFD[5569] = 5.11317705E+00;
    COFD[5570] = -3.94506858E-01;
    COFD[5571] = 1.49416998E-02;
    COFD[5572] = -2.25770498E+01;
    COFD[5573] = 5.35063033E+00;
    COFD[5574] = -4.34734059E-01;
    COFD[5575] = 1.70582710E-02;
    COFD[5576] = -2.25818241E+01;
    COFD[5577] = 5.35063033E+00;
    COFD[5578] = -4.34734059E-01;
    COFD[5579] = 1.70582710E-02;
    COFD[5580] = -2.22928354E+01;
    COFD[5581] = 5.15132481E+00;
    COFD[5582] = -4.00804233E-01;
    COFD[5583] = 1.52679625E-02;
    COFD[5584] = -2.27864348E+01;
    COFD[5585] = 5.52719085E+00;
    COFD[5586] = -4.68021516E-01;
    COFD[5587] = 1.89143588E-02;
    COFD[5588] = -2.28367691E+01;
    COFD[5589] = 5.52719085E+00;
    COFD[5590] = -4.68021516E-01;
    COFD[5591] = 1.89143588E-02;
    COFD[5592] = -2.28367691E+01;
    COFD[5593] = 5.52719085E+00;
    COFD[5594] = -4.68021516E-01;
    COFD[5595] = 1.89143588E-02;
    COFD[5596] = -2.28434402E+01;
    COFD[5597] = 5.52719085E+00;
    COFD[5598] = -4.68021516E-01;
    COFD[5599] = 1.89143588E-02;
    COFD[5600] = -2.28434402E+01;
    COFD[5601] = 5.52719085E+00;
    COFD[5602] = -4.68021516E-01;
    COFD[5603] = 1.89143588E-02;
    COFD[5604] = -2.28859084E+01;
    COFD[5605] = 5.52719085E+00;
    COFD[5606] = -4.68021516E-01;
    COFD[5607] = 1.89143588E-02;
    COFD[5608] = -2.29227645E+01;
    COFD[5609] = 5.52719085E+00;
    COFD[5610] = -4.68021516E-01;
    COFD[5611] = 1.89143588E-02;
    COFD[5612] = -1.95211303E+01;
    COFD[5613] = 4.75414932E+00;
    COFD[5614] = -3.99807084E-01;
    COFD[5615] = 1.72356137E-02;
    COFD[5616] = -1.89278805E+01;
    COFD[5617] = 5.20981862E+00;
    COFD[5618] = -4.52489825E-01;
    COFD[5619] = 1.92609226E-02;
    COFD[5620] = -1.49181094E+01;
    COFD[5621] = 3.56974825E+00;
    COFD[5622] = -2.52221138E-01;
    COFD[5623] = 1.10819767E-02;
    COFD[5624] = -2.06185530E+01;
    COFD[5625] = 5.20087227E+00;
    COFD[5626] = -4.51444972E-01;
    COFD[5627] = 1.92201898E-02;
    COFD[5628] = -2.06185530E+01;
    COFD[5629] = 5.20087227E+00;
    COFD[5630] = -4.51444972E-01;
    COFD[5631] = 1.92201898E-02;
    COFD[5632] = -2.06485216E+01;
    COFD[5633] = 5.20087227E+00;
    COFD[5634] = -4.51444972E-01;
    COFD[5635] = 1.92201898E-02;
    COFD[5636] = -1.83756296E+01;
    COFD[5637] = 4.49977587E+00;
    COFD[5638] = -3.68989022E-01;
    COFD[5639] = 1.59879891E-02;
    COFD[5640] = -2.05955883E+01;
    COFD[5641] = 5.17771473E+00;
    COFD[5642] = -4.48750201E-01;
    COFD[5643] = 1.91155567E-02;
    COFD[5644] = -1.84015347E+01;
    COFD[5645] = 4.49977587E+00;
    COFD[5646] = -3.68989022E-01;
    COFD[5647] = 1.59879891E-02;
    COFD[5648] = -2.05801081E+01;
    COFD[5649] = 4.74677479E+00;
    COFD[5650] = -3.36160335E-01;
    COFD[5651] = 1.19858600E-02;
    COFD[5652] = -2.22667492E+01;
    COFD[5653] = 5.59047891E+00;
    COFD[5654] = -4.85359237E-01;
    COFD[5655] = 2.00302091E-02;
    COFD[5656] = -2.22814898E+01;
    COFD[5657] = 5.59047891E+00;
    COFD[5658] = -4.85359237E-01;
    COFD[5659] = 2.00302091E-02;
    COFD[5660] = -1.95710942E+01;
    COFD[5661] = 4.76145602E+00;
    COFD[5662] = -4.00684587E-01;
    COFD[5663] = 1.72708322E-02;
    COFD[5664] = -2.21842699E+01;
    COFD[5665] = 5.60509076E+00;
    COFD[5666] = -4.91228646E-01;
    COFD[5667] = 2.04425499E-02;
    COFD[5668] = -2.16070103E+01;
    COFD[5669] = 5.05104142E+00;
    COFD[5670] = -3.84428177E-01;
    COFD[5671] = 1.44249285E-02;
    COFD[5672] = -2.24449636E+01;
    COFD[5673] = 5.61213963E+00;
    COFD[5674] = -4.90953451E-01;
    COFD[5675] = 2.03841518E-02;
    COFD[5676] = -2.16199377E+01;
    COFD[5677] = 5.05104142E+00;
    COFD[5678] = -3.84428177E-01;
    COFD[5679] = 1.44249285E-02;
    COFD[5680] = -2.24578673E+01;
    COFD[5681] = 5.61213963E+00;
    COFD[5682] = -4.90953451E-01;
    COFD[5683] = 2.03841518E-02;
    COFD[5684] = -2.21119432E+01;
    COFD[5685] = 5.31843386E+00;
    COFD[5686] = -4.29198961E-01;
    COFD[5687] = 1.67645078E-02;
    COFD[5688] = -2.21119432E+01;
    COFD[5689] = 5.31843386E+00;
    COFD[5690] = -4.29198961E-01;
    COFD[5691] = 1.67645078E-02;
    COFD[5692] = -1.98645237E+01;
    COFD[5693] = 4.86038833E+00;
    COFD[5694] = -4.12187543E-01;
    COFD[5695] = 1.77155642E-02;
    COFD[5696] = -1.98759845E+01;
    COFD[5697] = 4.86038833E+00;
    COFD[5698] = -4.12187543E-01;
    COFD[5699] = 1.77155642E-02;
    COFD[5700] = -1.98870113E+01;
    COFD[5701] = 4.86038833E+00;
    COFD[5702] = -4.12187543E-01;
    COFD[5703] = 1.77155642E-02;
    COFD[5704] = -2.24470642E+01;
    COFD[5705] = 5.61208013E+00;
    COFD[5706] = -4.91433206E-01;
    COFD[5707] = 2.04241259E-02;
    COFD[5708] = -2.21989542E+01;
    COFD[5709] = 5.25559319E+00;
    COFD[5710] = -4.18539684E-01;
    COFD[5711] = 1.62026965E-02;
    COFD[5712] = -2.27225058E+01;
    COFD[5713] = 5.52719085E+00;
    COFD[5714] = -4.68021516E-01;
    COFD[5715] = 1.89143588E-02;
    COFD[5716] = -2.21444625E+01;
    COFD[5717] = 5.15225414E+00;
    COFD[5718] = -4.00960404E-01;
    COFD[5719] = 1.52761358E-02;
    COFD[5720] = -2.27299046E+01;
    COFD[5721] = 5.52719085E+00;
    COFD[5722] = -4.68021516E-01;
    COFD[5723] = 1.89143588E-02;
    COFD[5724] = -2.19014365E+01;
    COFD[5725] = 5.11317705E+00;
    COFD[5726] = -3.94506858E-01;
    COFD[5727] = 1.49416998E-02;
    COFD[5728] = -2.26136853E+01;
    COFD[5729] = 5.35063033E+00;
    COFD[5730] = -4.34734059E-01;
    COFD[5731] = 1.70582710E-02;
    COFD[5732] = -2.26188244E+01;
    COFD[5733] = 5.35063033E+00;
    COFD[5734] = -4.34734059E-01;
    COFD[5735] = 1.70582710E-02;
    COFD[5736] = -2.23301955E+01;
    COFD[5737] = 5.15132481E+00;
    COFD[5738] = -4.00804233E-01;
    COFD[5739] = 1.52679625E-02;
    COFD[5740] = -2.28241497E+01;
    COFD[5741] = 5.52719085E+00;
    COFD[5742] = -4.68021516E-01;
    COFD[5743] = 1.89143588E-02;
    COFD[5744] = -2.28786502E+01;
    COFD[5745] = 5.52719085E+00;
    COFD[5746] = -4.68021516E-01;
    COFD[5747] = 1.89143588E-02;
    COFD[5748] = -2.28786502E+01;
    COFD[5749] = 5.52719085E+00;
    COFD[5750] = -4.68021516E-01;
    COFD[5751] = 1.89143588E-02;
    COFD[5752] = -2.28859084E+01;
    COFD[5753] = 5.52719085E+00;
    COFD[5754] = -4.68021516E-01;
    COFD[5755] = 1.89143588E-02;
    COFD[5756] = -2.28859084E+01;
    COFD[5757] = 5.52719085E+00;
    COFD[5758] = -4.68021516E-01;
    COFD[5759] = 1.89143588E-02;
    COFD[5760] = -2.29323214E+01;
    COFD[5761] = 5.52719085E+00;
    COFD[5762] = -4.68021516E-01;
    COFD[5763] = 1.89143588E-02;
    COFD[5764] = -2.29729119E+01;
    COFD[5765] = 5.52719085E+00;
    COFD[5766] = -4.68021516E-01;
    COFD[5767] = 1.89143588E-02;
    COFD[5768] = -1.95433253E+01;
    COFD[5769] = 4.75414932E+00;
    COFD[5770] = -3.99807084E-01;
    COFD[5771] = 1.72356137E-02;
    COFD[5772] = -1.89287258E+01;
    COFD[5773] = 5.20981862E+00;
    COFD[5774] = -4.52489825E-01;
    COFD[5775] = 1.92609226E-02;
    COFD[5776] = -1.49197832E+01;
    COFD[5777] = 3.56974825E+00;
    COFD[5778] = -2.52221138E-01;
    COFD[5779] = 1.10819767E-02;
    COFD[5780] = -2.06289714E+01;
    COFD[5781] = 5.20087227E+00;
    COFD[5782] = -4.51444972E-01;
    COFD[5783] = 1.92201898E-02;
    COFD[5784] = -2.06289714E+01;
    COFD[5785] = 5.20087227E+00;
    COFD[5786] = -4.51444972E-01;
    COFD[5787] = 1.92201898E-02;
    COFD[5788] = -2.06595907E+01;
    COFD[5789] = 5.20087227E+00;
    COFD[5790] = -4.51444972E-01;
    COFD[5791] = 1.92201898E-02;
    COFD[5792] = -1.83873107E+01;
    COFD[5793] = 4.49977587E+00;
    COFD[5794] = -3.68989022E-01;
    COFD[5795] = 1.59879891E-02;
    COFD[5796] = -2.06072968E+01;
    COFD[5797] = 5.17771473E+00;
    COFD[5798] = -4.48750201E-01;
    COFD[5799] = 1.91155567E-02;
    COFD[5800] = -1.84138446E+01;
    COFD[5801] = 4.49977587E+00;
    COFD[5802] = -3.68989022E-01;
    COFD[5803] = 1.59879891E-02;
    COFD[5804] = -2.05930361E+01;
    COFD[5805] = 4.74677479E+00;
    COFD[5806] = -3.36160335E-01;
    COFD[5807] = 1.19858600E-02;
    COFD[5808] = -2.22842444E+01;
    COFD[5809] = 5.59047891E+00;
    COFD[5810] = -4.85359237E-01;
    COFD[5811] = 2.00302091E-02;
    COFD[5812] = -2.22995181E+01;
    COFD[5813] = 5.59047891E+00;
    COFD[5814] = -4.85359237E-01;
    COFD[5815] = 2.00302091E-02;
    COFD[5816] = -1.95896245E+01;
    COFD[5817] = 4.76145602E+00;
    COFD[5818] = -4.00684587E-01;
    COFD[5819] = 1.72708322E-02;
    COFD[5820] = -2.22028227E+01;
    COFD[5821] = 5.60509076E+00;
    COFD[5822] = -4.91228646E-01;
    COFD[5823] = 2.04425499E-02;
    COFD[5824] = -2.16260574E+01;
    COFD[5825] = 5.05104142E+00;
    COFD[5826] = -3.84428177E-01;
    COFD[5827] = 1.44249285E-02;
    COFD[5828] = -2.24640329E+01;
    COFD[5829] = 5.61213963E+00;
    COFD[5830] = -4.90953451E-01;
    COFD[5831] = 2.03841518E-02;
    COFD[5832] = -2.16394936E+01;
    COFD[5833] = 5.05104142E+00;
    COFD[5834] = -3.84428177E-01;
    COFD[5835] = 1.44249285E-02;
    COFD[5836] = -2.24774450E+01;
    COFD[5837] = 5.61213963E+00;
    COFD[5838] = -4.90953451E-01;
    COFD[5839] = 2.03841518E-02;
    COFD[5840] = -2.21320000E+01;
    COFD[5841] = 5.31843386E+00;
    COFD[5842] = -4.29198961E-01;
    COFD[5843] = 1.67645078E-02;
    COFD[5844] = -2.21320000E+01;
    COFD[5845] = 5.31843386E+00;
    COFD[5846] = -4.29198961E-01;
    COFD[5847] = 1.67645078E-02;
    COFD[5848] = -1.98850525E+01;
    COFD[5849] = 4.86038833E+00;
    COFD[5850] = -4.12187543E-01;
    COFD[5851] = 1.77155642E-02;
    COFD[5852] = -1.98969995E+01;
    COFD[5853] = 4.86038833E+00;
    COFD[5854] = -4.12187543E-01;
    COFD[5855] = 1.77155642E-02;
    COFD[5856] = -1.99085051E+01;
    COFD[5857] = 4.86038833E+00;
    COFD[5858] = -4.12187543E-01;
    COFD[5859] = 1.77155642E-02;
    COFD[5860] = -2.24729433E+01;
    COFD[5861] = 5.61208013E+00;
    COFD[5862] = -4.91433206E-01;
    COFD[5863] = 2.04241259E-02;
    COFD[5864] = -2.22248512E+01;
    COFD[5865] = 5.25559319E+00;
    COFD[5866] = -4.18539684E-01;
    COFD[5867] = 1.62026965E-02;
    COFD[5868] = -2.27488112E+01;
    COFD[5869] = 5.52719085E+00;
    COFD[5870] = -4.68021516E-01;
    COFD[5871] = 1.89143588E-02;
    COFD[5872] = -2.21711534E+01;
    COFD[5873] = 5.15225414E+00;
    COFD[5874] = -4.00960404E-01;
    COFD[5875] = 1.52761358E-02;
    COFD[5876] = -2.27566129E+01;
    COFD[5877] = 5.52719085E+00;
    COFD[5878] = -4.68021516E-01;
    COFD[5879] = 1.89143588E-02;
    COFD[5880] = -2.19285250E+01;
    COFD[5881] = 5.11317705E+00;
    COFD[5882] = -3.94506858E-01;
    COFD[5883] = 1.49416998E-02;
    COFD[5884] = -2.26451233E+01;
    COFD[5885] = 5.35063033E+00;
    COFD[5886] = -4.34734059E-01;
    COFD[5887] = 1.70582710E-02;
    COFD[5888] = -2.26505977E+01;
    COFD[5889] = 5.35063033E+00;
    COFD[5890] = -4.34734059E-01;
    COFD[5891] = 1.70582710E-02;
    COFD[5892] = -2.23622999E+01;
    COFD[5893] = 5.15132481E+00;
    COFD[5894] = -4.00804233E-01;
    COFD[5895] = 1.52679625E-02;
    COFD[5896] = -2.28565811E+01;
    COFD[5897] = 5.52719085E+00;
    COFD[5898] = -4.68021516E-01;
    COFD[5899] = 1.89143588E-02;
    COFD[5900] = -2.29149554E+01;
    COFD[5901] = 5.52719085E+00;
    COFD[5902] = -4.68021516E-01;
    COFD[5903] = 1.89143588E-02;
    COFD[5904] = -2.29149554E+01;
    COFD[5905] = 5.52719085E+00;
    COFD[5906] = -4.68021516E-01;
    COFD[5907] = 1.89143588E-02;
    COFD[5908] = -2.29227645E+01;
    COFD[5909] = 5.52719085E+00;
    COFD[5910] = -4.68021516E-01;
    COFD[5911] = 1.89143588E-02;
    COFD[5912] = -2.29227645E+01;
    COFD[5913] = 5.52719085E+00;
    COFD[5914] = -4.68021516E-01;
    COFD[5915] = 1.89143588E-02;
    COFD[5916] = -2.29729119E+01;
    COFD[5917] = 5.52719085E+00;
    COFD[5918] = -4.68021516E-01;
    COFD[5919] = 1.89143588E-02;
    COFD[5920] = -2.30170910E+01;
    COFD[5921] = 5.52719085E+00;
    COFD[5922] = -4.68021516E-01;
    COFD[5923] = 1.89143588E-02;
    COFD[5924] = -1.95618571E+01;
    COFD[5925] = 4.75414932E+00;
    COFD[5926] = -3.99807084E-01;
    COFD[5927] = 1.72356137E-02;
    COFD[5928] = -1.42894441E+01;
    COFD[5929] = 3.67490723E+00;
    COFD[5930] = -2.65114792E-01;
    COFD[5931] = 1.16092671E-02;
    COFD[5932] = -1.16906297E+01;
    COFD[5933] = 2.47469981E+00;
    COFD[5934] = -1.10436257E-01;
    COFD[5935] = 4.95273813E-03;
    COFD[5936] = -1.59404882E+01;
    COFD[5937] = 3.66853818E+00;
    COFD[5938] = -2.64346221E-01;
    COFD[5939] = 1.15784613E-02;
    COFD[5940] = -1.59404882E+01;
    COFD[5941] = 3.66853818E+00;
    COFD[5942] = -2.64346221E-01;
    COFD[5943] = 1.15784613E-02;
    COFD[5944] = -1.59633387E+01;
    COFD[5945] = 3.66853818E+00;
    COFD[5946] = -2.64346221E-01;
    COFD[5947] = 1.15784613E-02;
    COFD[5948] = -1.40756935E+01;
    COFD[5949] = 3.07549274E+00;
    COFD[5950] = -1.88889344E-01;
    COFD[5951] = 8.37152866E-03;
    COFD[5952] = -1.59327297E+01;
    COFD[5953] = 3.65620899E+00;
    COFD[5954] = -2.62933804E-01;
    COFD[5955] = 1.15253223E-02;
    COFD[5956] = -1.40949196E+01;
    COFD[5957] = 3.07549274E+00;
    COFD[5958] = -1.88889344E-01;
    COFD[5959] = 8.37152866E-03;
    COFD[5960] = -2.10643259E+01;
    COFD[5961] = 5.53614847E+00;
    COFD[5962] = -4.86046736E-01;
    COFD[5963] = 2.03659188E-02;
    COFD[5964] = -1.83039618E+01;
    COFD[5965] = 4.47952077E+00;
    COFD[5966] = -3.66569471E-01;
    COFD[5967] = 1.58916129E-02;
    COFD[5968] = -1.83137139E+01;
    COFD[5969] = 4.47952077E+00;
    COFD[5970] = -3.66569471E-01;
    COFD[5971] = 1.58916129E-02;
    COFD[5972] = -1.50031687E+01;
    COFD[5973] = 3.26223357E+00;
    COFD[5974] = -2.12746642E-01;
    COFD[5975] = 9.38912883E-03;
    COFD[5976] = -1.78815889E+01;
    COFD[5977] = 4.34347890E+00;
    COFD[5978] = -3.49890003E-01;
    COFD[5979] = 1.52083459E-02;
    COFD[5980] = -2.04750581E+01;
    COFD[5981] = 5.23112374E+00;
    COFD[5982] = -4.54967682E-01;
    COFD[5983] = 1.93570423E-02;
    COFD[5984] = -1.82590824E+01;
    COFD[5985] = 4.39538102E+00;
    COFD[5986] = -3.56367230E-01;
    COFD[5987] = 1.54788461E-02;
    COFD[5988] = -2.04833713E+01;
    COFD[5989] = 5.23112374E+00;
    COFD[5990] = -4.54967682E-01;
    COFD[5991] = 1.93570423E-02;
    COFD[5992] = -1.82673770E+01;
    COFD[5993] = 4.39538102E+00;
    COFD[5994] = -3.56367230E-01;
    COFD[5995] = 1.54788461E-02;
    COFD[5996] = -2.02268902E+01;
    COFD[5997] = 5.13632093E+00;
    COFD[5998] = -4.44839124E-01;
    COFD[5999] = 1.90058354E-02;
    COFD[6000] = -2.02268902E+01;
    COFD[6001] = 5.13632093E+00;
    COFD[6002] = -4.44839124E-01;
    COFD[6003] = 1.90058354E-02;
    COFD[6004] = -1.52414485E+01;
    COFD[6005] = 3.35922578E+00;
    COFD[6006] = -2.25181399E-01;
    COFD[6007] = 9.92132878E-03;
    COFD[6008] = -1.52486273E+01;
    COFD[6009] = 3.35922578E+00;
    COFD[6010] = -2.25181399E-01;
    COFD[6011] = 9.92132878E-03;
    COFD[6012] = -1.52554761E+01;
    COFD[6013] = 3.35922578E+00;
    COFD[6014] = -2.25181399E-01;
    COFD[6015] = 9.92132878E-03;
    COFD[6016] = -1.81432461E+01;
    COFD[6017] = 4.37565431E+00;
    COFD[6018] = -3.53906025E-01;
    COFD[6019] = 1.53760786E-02;
    COFD[6020] = -2.02738958E+01;
    COFD[6021] = 5.10426133E+00;
    COFD[6022] = -4.41256919E-01;
    COFD[6023] = 1.88737290E-02;
    COFD[6024] = -1.94344389E+01;
    COFD[6025] = 4.75414932E+00;
    COFD[6026] = -3.99807084E-01;
    COFD[6027] = 1.72356137E-02;
    COFD[6028] = -2.05842618E+01;
    COFD[6029] = 5.16117916E+00;
    COFD[6030] = -4.46897404E-01;
    COFD[6031] = 1.90470443E-02;
    COFD[6032] = -1.94386503E+01;
    COFD[6033] = 4.75414932E+00;
    COFD[6034] = -3.99807084E-01;
    COFD[6035] = 1.72356137E-02;
    COFD[6036] = -2.04649069E+01;
    COFD[6037] = 5.18856872E+00;
    COFD[6038] = -4.50001829E-01;
    COFD[6039] = 1.91636142E-02;
    COFD[6040] = -2.02659445E+01;
    COFD[6041] = 5.01994886E+00;
    COFD[6042] = -4.31380254E-01;
    COFD[6043] = 1.84875367E-02;
    COFD[6044] = -2.02686523E+01;
    COFD[6045] = 5.01994886E+00;
    COFD[6046] = -4.31380254E-01;
    COFD[6047] = 1.84875367E-02;
    COFD[6048] = -2.07441338E+01;
    COFD[6049] = 5.16159824E+00;
    COFD[6050] = -4.46935718E-01;
    COFD[6051] = 1.90480458E-02;
    COFD[6052] = -1.94899278E+01;
    COFD[6053] = 4.75414932E+00;
    COFD[6054] = -3.99807084E-01;
    COFD[6055] = 1.72356137E-02;
    COFD[6056] = -1.95175620E+01;
    COFD[6057] = 4.75414932E+00;
    COFD[6058] = -3.99807084E-01;
    COFD[6059] = 1.72356137E-02;
    COFD[6060] = -1.95175620E+01;
    COFD[6061] = 4.75414932E+00;
    COFD[6062] = -3.99807084E-01;
    COFD[6063] = 1.72356137E-02;
    COFD[6064] = -1.95211303E+01;
    COFD[6065] = 4.75414932E+00;
    COFD[6066] = -3.99807084E-01;
    COFD[6067] = 1.72356137E-02;
    COFD[6068] = -1.95211303E+01;
    COFD[6069] = 4.75414932E+00;
    COFD[6070] = -3.99807084E-01;
    COFD[6071] = 1.72356137E-02;
    COFD[6072] = -1.95433253E+01;
    COFD[6073] = 4.75414932E+00;
    COFD[6074] = -3.99807084E-01;
    COFD[6075] = 1.72356137E-02;
    COFD[6076] = -1.95618571E+01;
    COFD[6077] = 4.75414932E+00;
    COFD[6078] = -3.99807084E-01;
    COFD[6079] = 1.72356137E-02;
    COFD[6080] = -1.49828430E+01;
    COFD[6081] = 3.25781069E+00;
    COFD[6082] = -2.12199367E-01;
    COFD[6083] = 9.36657283E-03;
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
    COFTD[4] = 1.44152190E-01;
    COFTD[5] = 7.99993584E-05;
    COFTD[6] = -4.89707442E-08;
    COFTD[7] = 9.14277269E-12;
    COFTD[8] = 9.90752318E-02;
    COFTD[9] = 6.44201384E-04;
    COFTD[10] = -3.38485953E-07;
    COFTD[11] = 5.57356746E-11;
    COFTD[12] = 9.90752318E-02;
    COFTD[13] = 6.44201384E-04;
    COFTD[14] = -3.38485953E-07;
    COFTD[15] = 5.57356746E-11;
    COFTD[16] = 1.00039110E-01;
    COFTD[17] = 6.50468660E-04;
    COFTD[18] = -3.41778999E-07;
    COFTD[19] = 5.62779132E-11;
    COFTD[20] = 2.35283119E-01;
    COFTD[21] = 4.65670599E-04;
    COFTD[22] = -2.60939824E-07;
    COFTD[23] = 4.49271822E-11;
    COFTD[24] = 1.05124122E-01;
    COFTD[25] = 6.50665957E-04;
    COFTD[26] = -3.42564538E-07;
    COFTD[27] = 5.64804120E-11;
    COFTD[28] = 2.37053352E-01;
    COFTD[29] = 4.69174231E-04;
    COFTD[30] = -2.62903094E-07;
    COFTD[31] = 4.52652072E-11;
    COFTD[32] = -1.74352698E-01;
    COFTD[33] = 8.62246873E-04;
    COFTD[34] = -3.79545489E-07;
    COFTD[35] = 5.60262093E-11;
    COFTD[36] = -3.81470765E-02;
    COFTD[37] = 8.39833490E-04;
    COFTD[38] = -4.11688915E-07;
    COFTD[39] = 6.49124952E-11;
    COFTD[40] = -3.82574649E-02;
    COFTD[41] = 8.42263764E-04;
    COFTD[42] = -4.12880242E-07;
    COFTD[43] = 6.51003362E-11;
    COFTD[44] = 2.00119897E-01;
    COFTD[45] = 5.64793704E-04;
    COFTD[46] = -3.09445484E-07;
    COFTD[47] = 5.24139335E-11;
    COFTD[48] = -1.42100396E-02;
    COFTD[49] = 8.23812102E-04;
    COFTD[50] = -4.08995515E-07;
    COFTD[51] = 6.49899310E-11;
    COFTD[52] = -1.60981264E-01;
    COFTD[53] = 9.03807572E-04;
    COFTD[54] = -4.06927941E-07;
    COFTD[55] = 6.09202254E-11;
    COFTD[56] = -2.28105944E-02;
    COFTD[57] = 8.33470403E-04;
    COFTD[58] = -4.11969112E-07;
    COFTD[59] = 6.52859371E-11;
    COFTD[60] = -1.61357564E-01;
    COFTD[61] = 9.05920260E-04;
    COFTD[62] = -4.07879153E-07;
    COFTD[63] = 6.10626290E-11;
    COFTD[64] = -2.28637575E-02;
    COFTD[65] = 8.35412914E-04;
    COFTD[66] = -4.12929260E-07;
    COFTD[67] = 6.54380945E-11;
    COFTD[68] = -1.31244519E-01;
    COFTD[69] = 9.03901384E-04;
    COFTD[70] = -4.17831507E-07;
    COFTD[71] = 6.35725667E-11;
    COFTD[72] = -1.31244519E-01;
    COFTD[73] = 9.03901384E-04;
    COFTD[74] = -4.17831507E-07;
    COFTD[75] = 6.35725667E-11;
    COFTD[76] = 1.79840299E-01;
    COFTD[77] = 6.01722902E-04;
    COFTD[78] = -3.26433894E-07;
    COFTD[79] = 5.49112302E-11;
    COFTD[80] = 1.80186965E-01;
    COFTD[81] = 6.02882805E-04;
    COFTD[82] = -3.27063140E-07;
    COFTD[83] = 5.50170790E-11;
    COFTD[84] = 1.80513677E-01;
    COFTD[85] = 6.03975942E-04;
    COFTD[86] = -3.27656165E-07;
    COFTD[87] = 5.51168351E-11;
    COFTD[88] = -2.00309448E-02;
    COFTD[89] = 8.50440115E-04;
    COFTD[90] = -4.21064468E-07;
    COFTD[91] = 6.67959710E-11;
    COFTD[92] = -1.41951848E-01;
    COFTD[93] = 9.23429679E-04;
    COFTD[94] = -4.24140376E-07;
    COFTD[95] = 6.42810196E-11;
    COFTD[96] = -8.70374736E-02;
    COFTD[97] = 9.01418310E-04;
    COFTD[98] = -4.30150021E-07;
    COFTD[99] = 6.67243061E-11;
    COFTD[100] = -1.55682554E-01;
    COFTD[101] = 9.27159178E-04;
    COFTD[102] = -4.21074432E-07;
    COFTD[103] = 6.33759629E-11;
    COFTD[104] = -8.71227516E-02;
    COFTD[105] = 9.02301506E-04;
    COFTD[106] = -4.30571476E-07;
    COFTD[107] = 6.67896816E-11;
    COFTD[108] = -1.59826932E-01;
    COFTD[109] = 9.28231324E-04;
    COFTD[110] = -4.20059750E-07;
    COFTD[111] = 6.30844146E-11;
    COFTD[112] = -1.30522623E-01;
    COFTD[113] = 9.30761849E-04;
    COFTD[114] = -4.31797867E-07;
    COFTD[115] = 6.58417879E-11;
    COFTD[116] = -1.30597464E-01;
    COFTD[117] = 9.31295546E-04;
    COFTD[118] = -4.32045459E-07;
    COFTD[119] = 6.58795415E-11;
    COFTD[120] = -1.57481646E-01;
    COFTD[121] = 9.37224347E-04;
    COFTD[122] = -4.25604687E-07;
    COFTD[123] = 6.40540242E-11;
    COFTD[124] = -8.81113877E-02;
    COFTD[125] = 9.12540483E-04;
    COFTD[126] = -4.35457439E-07;
    COFTD[127] = 6.75475857E-11;
    COFTD[128] = -8.86077023E-02;
    COFTD[129] = 9.17680649E-04;
    COFTD[130] = -4.37910287E-07;
    COFTD[131] = 6.79280683E-11;
    COFTD[132] = -8.86077023E-02;
    COFTD[133] = 9.17680649E-04;
    COFTD[134] = -4.37910287E-07;
    COFTD[135] = 6.79280683E-11;
    COFTD[136] = -8.86700049E-02;
    COFTD[137] = 9.18325897E-04;
    COFTD[138] = -4.38218194E-07;
    COFTD[139] = 6.79758305E-11;
    COFTD[140] = -8.86700049E-02;
    COFTD[141] = 9.18325897E-04;
    COFTD[142] = -4.38218194E-07;
    COFTD[143] = 6.79758305E-11;
    COFTD[144] = -8.90486416E-02;
    COFTD[145] = 9.22247312E-04;
    COFTD[146] = -4.40089464E-07;
    COFTD[147] = 6.82660994E-11;
    COFTD[148] = -8.93533147E-02;
    COFTD[149] = 9.25402710E-04;
    COFTD[150] = -4.41595196E-07;
    COFTD[151] = 6.84996666E-11;
    COFTD[152] = 2.01521643E-01;
    COFTD[153] = 5.62744089E-04;
    COFTD[154] = -3.08519239E-07;
    COFTD[155] = 5.22805986E-11;
    COFTD[156] = -1.44152190E-01;
    COFTD[157] = -7.99993584E-05;
    COFTD[158] = 4.89707442E-08;
    COFTD[159] = -9.14277269E-12;
    COFTD[160] = 0.00000000E+00;
    COFTD[161] = 0.00000000E+00;
    COFTD[162] = 0.00000000E+00;
    COFTD[163] = 0.00000000E+00;
    COFTD[164] = 3.24747031E-01;
    COFTD[165] = 1.77798548E-04;
    COFTD[166] = -1.08934732E-07;
    COFTD[167] = 2.03595881E-11;
    COFTD[168] = 3.24747031E-01;
    COFTD[169] = 1.77798548E-04;
    COFTD[170] = -1.08934732E-07;
    COFTD[171] = 2.03595881E-11;
    COFTD[172] = 3.31191185E-01;
    COFTD[173] = 1.81326714E-04;
    COFTD[174] = -1.11096391E-07;
    COFTD[175] = 2.07635959E-11;
    COFTD[176] = 4.06682492E-01;
    COFTD[177] = 3.84705248E-05;
    COFTD[178] = -2.54846868E-08;
    COFTD[179] = 5.86302354E-12;
    COFTD[180] = 3.39557243E-01;
    COFTD[181] = 1.79335036E-04;
    COFTD[182] = -1.10135705E-07;
    COFTD[183] = 2.06427239E-11;
    COFTD[184] = 4.12895615E-01;
    COFTD[185] = 3.90582612E-05;
    COFTD[186] = -2.58740310E-08;
    COFTD[187] = 5.95259633E-12;
    COFTD[188] = 2.27469146E-02;
    COFTD[189] = 6.73078907E-04;
    COFTD[190] = -3.40935843E-07;
    COFTD[191] = 5.48499211E-11;
    COFTD[192] = 2.58066832E-01;
    COFTD[193] = 4.05072593E-04;
    COFTD[194] = -2.30587443E-07;
    COFTD[195] = 4.01863841E-11;
    COFTD[196] = 2.59569092E-01;
    COFTD[197] = 4.07430603E-04;
    COFTD[198] = -2.31929740E-07;
    COFTD[199] = 4.04203173E-11;
    COFTD[200] = 4.30605547E-01;
    COFTD[201] = 9.35961902E-05;
    COFTD[202] = -6.03983623E-08;
    COFTD[203] = 1.23115170E-11;
    COFTD[204] = 2.82974392E-01;
    COFTD[205] = 3.73032949E-04;
    COFTD[206] = -2.14959161E-07;
    COFTD[207] = 3.78355155E-11;
    COFTD[208] = 1.22119780E-01;
    COFTD[209] = 6.18373616E-04;
    COFTD[210] = -3.28422593E-07;
    COFTD[211] = 5.44603522E-11;
    COFTD[212] = 2.76725963E-01;
    COFTD[213] = 3.87792818E-04;
    COFTD[214] = -2.22504581E-07;
    COFTD[215] = 3.90251143E-11;
    COFTD[216] = 1.22693382E-01;
    COFTD[217] = 6.21278143E-04;
    COFTD[218] = -3.29965208E-07;
    COFTD[219] = 5.47161548E-11;
    COFTD[220] = 2.78021896E-01;
    COFTD[221] = 3.89608886E-04;
    COFTD[222] = -2.23546590E-07;
    COFTD[223] = 3.92078724E-11;
    COFTD[224] = 1.40314191E-01;
    COFTD[225] = 6.01266129E-04;
    COFTD[226] = -3.21915137E-07;
    COFTD[227] = 5.36679068E-11;
    COFTD[228] = 1.40314191E-01;
    COFTD[229] = 6.01266129E-04;
    COFTD[230] = -3.21915137E-07;
    COFTD[231] = 5.36679068E-11;
    COFTD[232] = 4.26579943E-01;
    COFTD[233] = 1.20407274E-04;
    COFTD[234] = -7.67298757E-08;
    COFTD[235] = 1.52090336E-11;
    COFTD[236] = 4.28230888E-01;
    COFTD[237] = 1.20873273E-04;
    COFTD[238] = -7.70268349E-08;
    COFTD[239] = 1.52678954E-11;
    COFTD[240] = 4.29789463E-01;
    COFTD[241] = 1.21313199E-04;
    COFTD[242] = -7.73071792E-08;
    COFTD[243] = 1.53234639E-11;
    COFTD[244] = 2.93191523E-01;
    COFTD[245] = 4.01430006E-04;
    COFTD[246] = -2.30705763E-07;
    COFTD[247] = 4.05176586E-11;
    COFTD[248] = 1.59991186E-01;
    COFTD[249] = 6.05491303E-04;
    COFTD[250] = -3.26269573E-07;
    COFTD[251] = 5.46306751E-11;
    COFTD[252] = 2.26564429E-01;
    COFTD[253] = 5.10157718E-04;
    COFTD[254] = -2.83466602E-07;
    COFTD[255] = 4.85011808E-11;
    COFTD[256] = 1.42235861E-01;
    COFTD[257] = 6.32941361E-04;
    COFTD[258] = -3.38237792E-07;
    COFTD[259] = 5.63182114E-11;
    COFTD[260] = 2.27009269E-01;
    COFTD[261] = 5.11159370E-04;
    COFTD[262] = -2.84023165E-07;
    COFTD[263] = 4.85964088E-11;
    COFTD[264] = 1.36817715E-01;
    COFTD[265] = 6.41727473E-04;
    COFTD[266] = -3.42055963E-07;
    COFTD[267] = 5.68567648E-11;
    COFTD[268] = 1.80870457E-01;
    COFTD[269] = 5.95738417E-04;
    COFTD[270] = -3.23473700E-07;
    COFTD[271] = 5.44465614E-11;
    COFTD[272] = 1.81078117E-01;
    COFTD[273] = 5.96422392E-04;
    COFTD[274] = -3.23845085E-07;
    COFTD[275] = 5.45090723E-11;
    COFTD[276] = 1.45184792E-01;
    COFTD[277] = 6.46973651E-04;
    COFTD[278] = -3.45712529E-07;
    COFTD[279] = 5.75601214E-11;
    COFTD[280] = 2.32196293E-01;
    COFTD[281] = 5.22839051E-04;
    COFTD[282] = -2.90512922E-07;
    COFTD[283] = 4.97068072E-11;
    COFTD[284] = 2.34821243E-01;
    COFTD[285] = 5.28749682E-04;
    COFTD[286] = -2.93797135E-07;
    COFTD[287] = 5.02687366E-11;
    COFTD[288] = 2.34821243E-01;
    COFTD[289] = 5.28749682E-04;
    COFTD[290] = -2.93797135E-07;
    COFTD[291] = 5.02687366E-11;
    COFTD[292] = 2.35151753E-01;
    COFTD[293] = 5.29493894E-04;
    COFTD[294] = -2.94210652E-07;
    COFTD[295] = 5.03394895E-11;
    COFTD[296] = 2.35151753E-01;
    COFTD[297] = 5.29493894E-04;
    COFTD[298] = -2.94210652E-07;
    COFTD[299] = 5.03394895E-11;
    COFTD[300] = 2.37165199E-01;
    COFTD[301] = 5.34027593E-04;
    COFTD[302] = -2.96729780E-07;
    COFTD[303] = 5.07705126E-11;
    COFTD[304] = 2.38791360E-01;
    COFTD[305] = 5.37689239E-04;
    COFTD[306] = -2.98764355E-07;
    COFTD[307] = 5.11186288E-11;
    COFTD[308] = 4.31331269E-01;
    COFTD[309] = 9.20536800E-05;
    COFTD[310] = -5.94509616E-08;
    COFTD[311] = 1.21437993E-11;
}

#ifndef AMREX_USE_GPU
/* Initializes parameter database */
void CKINIT()
{

    // (0):  H + O2 <=> O + OH
    kiv[23] = {0,19,5,7};
    nuv[23] = {-1,-1,1,1};
    // (0):  H + O2 <=> O + OH
    fwd_A[23]     = 3547000000000000;
    fwd_beta[23]  = -0.40600000000000003;
    fwd_Ea[23]    = 16599;
    prefactor_units[23]  = 1.0000000000000002e-06;
    activation_units[23] = 0.50321666580471969;
    phase_units[23]      = pow(10,-12.000000);
    is_PD[23] = 0;
    nTB[23] = 0;

    // (1):  O + H2 <=> H + OH
    kiv[24] = {5,1,0,7};
    nuv[24] = {-1,-1,1,1};
    // (1):  O + H2 <=> H + OH
    fwd_A[24]     = 50800;
    fwd_beta[24]  = 2.6699999999999999;
    fwd_Ea[24]    = 6290;
    prefactor_units[24]  = 1.0000000000000002e-06;
    activation_units[24] = 0.50321666580471969;
    phase_units[24]      = pow(10,-12.000000);
    is_PD[24] = 0;
    nTB[24] = 0;

    // (2):  H2 + OH <=> H2O + H
    kiv[25] = {1,7,8,0};
    nuv[25] = {-1,-1,1,1};
    // (2):  H2 + OH <=> H2O + H
    fwd_A[25]     = 216000000;
    fwd_beta[25]  = 1.51;
    fwd_Ea[25]    = 3430;
    prefactor_units[25]  = 1.0000000000000002e-06;
    activation_units[25] = 0.50321666580471969;
    phase_units[25]      = pow(10,-12.000000);
    is_PD[25] = 0;
    nTB[25] = 0;

    // (3):  O + H2O <=> OH + OH
    kiv[26] = {5,8,7,7};
    nuv[26] = {-1,-1,1,1};
    // (3):  O + H2O <=> OH + OH
    fwd_A[26]     = 2970000;
    fwd_beta[26]  = 2.02;
    fwd_Ea[26]    = 13400;
    prefactor_units[26]  = 1.0000000000000002e-06;
    activation_units[26] = 0.50321666580471969;
    phase_units[26]      = pow(10,-12.000000);
    is_PD[26] = 0;
    nTB[26] = 0;

    // (4):  H2 + M <=> H + H + M
    kiv[11] = {1,0,0};
    nuv[11] = {-1,1,1};
    // (4):  H2 + M <=> H + H + M
    fwd_A[11]     = 4.577e+19;
    fwd_beta[11]  = -1.3999999999999999;
    fwd_Ea[11]    = 104380;
    prefactor_units[11]  = 1.0000000000000002e-06;
    activation_units[11] = 0.50321666580471969;
    phase_units[11]      = pow(10,-6.000000);
    is_PD[11] = 0;
    nTB[11] = 4;
    TB[11] = (amrex::Real *) malloc(4 * sizeof(amrex::Real));
    TBid[11] = (int *) malloc(4 * sizeof(int));
    TBid[11][0] = 1; TB[11][0] = 2.5; // H2
    TBid[11][1] = 8; TB[11][1] = 12; // H2O
    TBid[11][2] = 11; TB[11][2] = 1.8999999999999999; // CO
    TBid[11][3] = 22; TB[11][3] = 3.7999999999999998; // CO2

    // (5):  O + O + M <=> O2 + M
    kiv[12] = {5,5,19};
    nuv[12] = {-1,-1,1};
    // (5):  O + O + M <=> O2 + M
    fwd_A[12]     = 6165000000000000;
    fwd_beta[12]  = -0.5;
    fwd_Ea[12]    = 0;
    prefactor_units[12]  = 1.0000000000000002e-12;
    activation_units[12] = 0.50321666580471969;
    phase_units[12]      = pow(10,-12.000000);
    is_PD[12] = 0;
    nTB[12] = 4;
    TB[12] = (amrex::Real *) malloc(4 * sizeof(amrex::Real));
    TBid[12] = (int *) malloc(4 * sizeof(int));
    TBid[12][0] = 1; TB[12][0] = 2.5; // H2
    TBid[12][1] = 8; TB[12][1] = 12; // H2O
    TBid[12][2] = 11; TB[12][2] = 1.8999999999999999; // CO
    TBid[12][3] = 22; TB[12][3] = 3.7999999999999998; // CO2

    // (6):  O + H + M <=> OH + M
    kiv[13] = {5,0,7};
    nuv[13] = {-1,-1,1};
    // (6):  O + H + M <=> OH + M
    fwd_A[13]     = 4.714e+18;
    fwd_beta[13]  = -1;
    fwd_Ea[13]    = 0;
    prefactor_units[13]  = 1.0000000000000002e-12;
    activation_units[13] = 0.50321666580471969;
    phase_units[13]      = pow(10,-12.000000);
    is_PD[13] = 0;
    nTB[13] = 4;
    TB[13] = (amrex::Real *) malloc(4 * sizeof(amrex::Real));
    TBid[13] = (int *) malloc(4 * sizeof(int));
    TBid[13][0] = 1; TB[13][0] = 2.5; // H2
    TBid[13][1] = 8; TB[13][1] = 12; // H2O
    TBid[13][2] = 11; TB[13][2] = 1.8999999999999999; // CO
    TBid[13][3] = 22; TB[13][3] = 3.7999999999999998; // CO2

    // (7):  H + OH + M <=> H2O + M
    kiv[14] = {0,7,8};
    nuv[14] = {-1,-1,1};
    // (7):  H + OH + M <=> H2O + M
    fwd_A[14]     = 3.8000000000000004e+22;
    fwd_beta[14]  = -2;
    fwd_Ea[14]    = 0;
    prefactor_units[14]  = 1.0000000000000002e-12;
    activation_units[14] = 0.50321666580471969;
    phase_units[14]      = pow(10,-12.000000);
    is_PD[14] = 0;
    nTB[14] = 4;
    TB[14] = (amrex::Real *) malloc(4 * sizeof(amrex::Real));
    TBid[14] = (int *) malloc(4 * sizeof(int));
    TBid[14][0] = 1; TB[14][0] = 2.5; // H2
    TBid[14][1] = 8; TB[14][1] = 12; // H2O
    TBid[14][2] = 11; TB[14][2] = 1.8999999999999999; // CO
    TBid[14][3] = 22; TB[14][3] = 3.7999999999999998; // CO2

    // (8):  H + O2 (+M) <=> HO2 (+M)
    kiv[0] = {0,19,20};
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
    nTB[0] = 5;
    TB[0] = (amrex::Real *) malloc(5 * sizeof(amrex::Real));
    TBid[0] = (int *) malloc(5 * sizeof(int));
    TBid[0][0] = 1; TB[0][0] = 2; // H2
    TBid[0][1] = 8; TB[0][1] = 11; // H2O
    TBid[0][2] = 19; TB[0][2] = 0.78000000000000003; // O2
    TBid[0][3] = 11; TB[0][3] = 1.8999999999999999; // CO
    TBid[0][4] = 22; TB[0][4] = 3.7999999999999998; // CO2

    // (9):  HO2 + H <=> H2 + O2
    kiv[27] = {20,0,1,19};
    nuv[27] = {-1,-1,1,1};
    // (9):  HO2 + H <=> H2 + O2
    fwd_A[27]     = 16600000000000;
    fwd_beta[27]  = 0;
    fwd_Ea[27]    = 823;
    prefactor_units[27]  = 1.0000000000000002e-06;
    activation_units[27] = 0.50321666580471969;
    phase_units[27]      = pow(10,-12.000000);
    is_PD[27] = 0;
    nTB[27] = 0;

    // (10):  HO2 + H <=> OH + OH
    kiv[28] = {20,0,7,7};
    nuv[28] = {-1,-1,1,1};
    // (10):  HO2 + H <=> OH + OH
    fwd_A[28]     = 70790000000000;
    fwd_beta[28]  = 0;
    fwd_Ea[28]    = 295;
    prefactor_units[28]  = 1.0000000000000002e-06;
    activation_units[28] = 0.50321666580471969;
    phase_units[28]      = pow(10,-12.000000);
    is_PD[28] = 0;
    nTB[28] = 0;

    // (11):  HO2 + O <=> O2 + OH
    kiv[29] = {20,5,19,7};
    nuv[29] = {-1,-1,1,1};
    // (11):  HO2 + O <=> O2 + OH
    fwd_A[29]     = 32500000000000;
    fwd_beta[29]  = 0;
    fwd_Ea[29]    = 0;
    prefactor_units[29]  = 1.0000000000000002e-06;
    activation_units[29] = 0.50321666580471969;
    phase_units[29]      = pow(10,-12.000000);
    is_PD[29] = 0;
    nTB[29] = 0;

    // (12):  HO2 + OH <=> H2O + O2
    kiv[30] = {20,7,8,19};
    nuv[30] = {-1,-1,1,1};
    // (12):  HO2 + OH <=> H2O + O2
    fwd_A[30]     = 28900000000000;
    fwd_beta[30]  = 0;
    fwd_Ea[30]    = -497;
    prefactor_units[30]  = 1.0000000000000002e-06;
    activation_units[30] = 0.50321666580471969;
    phase_units[30]      = pow(10,-12.000000);
    is_PD[30] = 0;
    nTB[30] = 0;

    // (13):  HO2 + HO2 <=> H2O2 + O2
    kiv[31] = {20,20,21,19};
    nuv[31] = {-1,-1,1,1};
    // (13):  HO2 + HO2 <=> H2O2 + O2
    fwd_A[31]     = 420000000000000;
    fwd_beta[31]  = 0;
    fwd_Ea[31]    = 11982;
    prefactor_units[31]  = 1.0000000000000002e-06;
    activation_units[31] = 0.50321666580471969;
    phase_units[31]      = pow(10,-12.000000);
    is_PD[31] = 0;
    nTB[31] = 0;

    // (14):  HO2 + HO2 <=> H2O2 + O2
    kiv[32] = {20,20,21,19};
    nuv[32] = {-1,-1,1,1};
    // (14):  HO2 + HO2 <=> H2O2 + O2
    fwd_A[32]     = 130000000000;
    fwd_beta[32]  = 0;
    fwd_Ea[32]    = -1629.3;
    prefactor_units[32]  = 1.0000000000000002e-06;
    activation_units[32] = 0.50321666580471969;
    phase_units[32]      = pow(10,-12.000000);
    is_PD[32] = 0;
    nTB[32] = 0;

    // (15):  H2O2 (+M) <=> OH + OH (+M)
    kiv[1] = {21,7,7};
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
    nTB[1] = 4;
    TB[1] = (amrex::Real *) malloc(4 * sizeof(amrex::Real));
    TBid[1] = (int *) malloc(4 * sizeof(int));
    TBid[1][0] = 1; TB[1][0] = 2.5; // H2
    TBid[1][1] = 8; TB[1][1] = 12; // H2O
    TBid[1][2] = 11; TB[1][2] = 1.8999999999999999; // CO
    TBid[1][3] = 22; TB[1][3] = 3.7999999999999998; // CO2

    // (16):  H2O2 + H <=> H2O + OH
    kiv[33] = {21,0,8,7};
    nuv[33] = {-1,-1,1,1};
    // (16):  H2O2 + H <=> H2O + OH
    fwd_A[33]     = 24100000000000;
    fwd_beta[33]  = 0;
    fwd_Ea[33]    = 3970;
    prefactor_units[33]  = 1.0000000000000002e-06;
    activation_units[33] = 0.50321666580471969;
    phase_units[33]      = pow(10,-12.000000);
    is_PD[33] = 0;
    nTB[33] = 0;

    // (17):  H2O2 + H <=> HO2 + H2
    kiv[34] = {21,0,20,1};
    nuv[34] = {-1,-1,1,1};
    // (17):  H2O2 + H <=> HO2 + H2
    fwd_A[34]     = 48200000000000;
    fwd_beta[34]  = 0;
    fwd_Ea[34]    = 7950;
    prefactor_units[34]  = 1.0000000000000002e-06;
    activation_units[34] = 0.50321666580471969;
    phase_units[34]      = pow(10,-12.000000);
    is_PD[34] = 0;
    nTB[34] = 0;

    // (18):  H2O2 + O <=> OH + HO2
    kiv[35] = {21,5,7,20};
    nuv[35] = {-1,-1,1,1};
    // (18):  H2O2 + O <=> OH + HO2
    fwd_A[35]     = 9550000;
    fwd_beta[35]  = 2;
    fwd_Ea[35]    = 3970;
    prefactor_units[35]  = 1.0000000000000002e-06;
    activation_units[35] = 0.50321666580471969;
    phase_units[35]      = pow(10,-12.000000);
    is_PD[35] = 0;
    nTB[35] = 0;

    // (19):  H2O2 + OH <=> HO2 + H2O
    kiv[36] = {21,7,20,8};
    nuv[36] = {-1,-1,1,1};
    // (19):  H2O2 + OH <=> HO2 + H2O
    fwd_A[36]     = 1000000000000;
    fwd_beta[36]  = 0;
    fwd_Ea[36]    = 0;
    prefactor_units[36]  = 1.0000000000000002e-06;
    activation_units[36] = 0.50321666580471969;
    phase_units[36]      = pow(10,-12.000000);
    is_PD[36] = 0;
    nTB[36] = 0;

    // (20):  H2O2 + OH <=> HO2 + H2O
    kiv[37] = {21,7,20,8};
    nuv[37] = {-1,-1,1,1};
    // (20):  H2O2 + OH <=> HO2 + H2O
    fwd_A[37]     = 580000000000000;
    fwd_beta[37]  = 0;
    fwd_Ea[37]    = 9557;
    prefactor_units[37]  = 1.0000000000000002e-06;
    activation_units[37] = 0.50321666580471969;
    phase_units[37]      = pow(10,-12.000000);
    is_PD[37] = 0;
    nTB[37] = 0;

    // (21):  CO + O (+M) <=> CO2 (+M)
    kiv[10] = {11,5,22};
    nuv[10] = {-1,-1,1};
    // (21):  CO + O (+M) <=> CO2 (+M)
    fwd_A[10]     = 18000000000;
    fwd_beta[10]  = 0;
    fwd_Ea[10]    = 2384;
    low_A[10]     = 1.5500000000000001e+24;
    low_beta[10]  = -2.79;
    low_Ea[10]    = 4191;
    prefactor_units[10]  = 1.0000000000000002e-06;
    activation_units[10] = 0.50321666580471969;
    phase_units[10]      = pow(10,-12.000000);
    is_PD[10] = 1;
    nTB[10] = 4;
    TB[10] = (amrex::Real *) malloc(4 * sizeof(amrex::Real));
    TBid[10] = (int *) malloc(4 * sizeof(int));
    TBid[10][0] = 1; TB[10][0] = 2.5; // H2
    TBid[10][1] = 8; TB[10][1] = 12; // H2O
    TBid[10][2] = 11; TB[10][2] = 1.8999999999999999; // CO
    TBid[10][3] = 22; TB[10][3] = 3.7999999999999998; // CO2

    // (22):  CO + O2 <=> CO2 + O
    kiv[38] = {11,19,22,5};
    nuv[38] = {-1,-1,1,1};
    // (22):  CO + O2 <=> CO2 + O
    fwd_A[38]     = 2530000000000;
    fwd_beta[38]  = 0;
    fwd_Ea[38]    = 47700;
    prefactor_units[38]  = 1.0000000000000002e-06;
    activation_units[38] = 0.50321666580471969;
    phase_units[38]      = pow(10,-12.000000);
    is_PD[38] = 0;
    nTB[38] = 0;

    // (23):  CO + HO2 <=> CO2 + OH
    kiv[39] = {11,20,22,7};
    nuv[39] = {-1,-1,1,1};
    // (23):  CO + HO2 <=> CO2 + OH
    fwd_A[39]     = 30100000000000;
    fwd_beta[39]  = 0;
    fwd_Ea[39]    = 23000;
    prefactor_units[39]  = 1.0000000000000002e-06;
    activation_units[39] = 0.50321666580471969;
    phase_units[39]      = pow(10,-12.000000);
    is_PD[39] = 0;
    nTB[39] = 0;

    // (24):  CO + OH <=> CO2 + H
    kiv[40] = {11,7,22,0};
    nuv[40] = {-1,-1,1,1};
    // (24):  CO + OH <=> CO2 + H
    fwd_A[40]     = 222900;
    fwd_beta[40]  = 1.8899999999999999;
    fwd_Ea[40]    = -1158.7;
    prefactor_units[40]  = 1.0000000000000002e-06;
    activation_units[40] = 0.50321666580471969;
    phase_units[40]      = pow(10,-12.000000);
    is_PD[40] = 0;
    nTB[40] = 0;

    // (25):  HCO + M <=> H + CO + M
    kiv[15] = {13,0,11};
    nuv[15] = {-1,1,1};
    // (25):  HCO + M <=> H + CO + M
    fwd_A[15]     = 474850000000;
    fwd_beta[15]  = 0.65900000000000003;
    fwd_Ea[15]    = 14874;
    prefactor_units[15]  = 1.0000000000000002e-06;
    activation_units[15] = 0.50321666580471969;
    phase_units[15]      = pow(10,-6.000000);
    is_PD[15] = 0;
    nTB[15] = 4;
    TB[15] = (amrex::Real *) malloc(4 * sizeof(amrex::Real));
    TBid[15] = (int *) malloc(4 * sizeof(int));
    TBid[15][0] = 1; TB[15][0] = 2.5; // H2
    TBid[15][1] = 8; TB[15][1] = 6; // H2O
    TBid[15][2] = 11; TB[15][2] = 1.8999999999999999; // CO
    TBid[15][3] = 22; TB[15][3] = 3.7999999999999998; // CO2

    // (26):  HCO + O2 <=> CO + HO2
    kiv[41] = {13,19,11,20};
    nuv[41] = {-1,-1,1,1};
    // (26):  HCO + O2 <=> CO + HO2
    fwd_A[41]     = 7580000000000;
    fwd_beta[41]  = 0;
    fwd_Ea[41]    = 410;
    prefactor_units[41]  = 1.0000000000000002e-06;
    activation_units[41] = 0.50321666580471969;
    phase_units[41]      = pow(10,-12.000000);
    is_PD[41] = 0;
    nTB[41] = 0;

    // (27):  HCO + H <=> CO + H2
    kiv[42] = {13,0,11,1};
    nuv[42] = {-1,-1,1,1};
    // (27):  HCO + H <=> CO + H2
    fwd_A[42]     = 72300000000000;
    fwd_beta[42]  = 0;
    fwd_Ea[42]    = 0;
    prefactor_units[42]  = 1.0000000000000002e-06;
    activation_units[42] = 0.50321666580471969;
    phase_units[42]      = pow(10,-12.000000);
    is_PD[42] = 0;
    nTB[42] = 0;

    // (28):  HCO + O <=> CO + OH
    kiv[43] = {13,5,11,7};
    nuv[43] = {-1,-1,1,1};
    // (28):  HCO + O <=> CO + OH
    fwd_A[43]     = 30200000000000;
    fwd_beta[43]  = 0;
    fwd_Ea[43]    = 0;
    prefactor_units[43]  = 1.0000000000000002e-06;
    activation_units[43] = 0.50321666580471969;
    phase_units[43]      = pow(10,-12.000000);
    is_PD[43] = 0;
    nTB[43] = 0;

    // (29):  HCO + OH <=> CO + H2O
    kiv[44] = {13,7,11,8};
    nuv[44] = {-1,-1,1,1};
    // (29):  HCO + OH <=> CO + H2O
    fwd_A[44]     = 30200000000000;
    fwd_beta[44]  = 0;
    fwd_Ea[44]    = 0;
    prefactor_units[44]  = 1.0000000000000002e-06;
    activation_units[44] = 0.50321666580471969;
    phase_units[44]      = pow(10,-12.000000);
    is_PD[44] = 0;
    nTB[44] = 0;

    // (30):  HCO + O <=> CO2 + H
    kiv[45] = {13,5,22,0};
    nuv[45] = {-1,-1,1,1};
    // (30):  HCO + O <=> CO2 + H
    fwd_A[45]     = 30000000000000;
    fwd_beta[45]  = 0;
    fwd_Ea[45]    = 0;
    prefactor_units[45]  = 1.0000000000000002e-06;
    activation_units[45] = 0.50321666580471969;
    phase_units[45]      = pow(10,-12.000000);
    is_PD[45] = 0;
    nTB[45] = 0;

    // (31):  HCO + HO2 <=> CO2 + OH + H
    kiv[46] = {13,20,22,7,0};
    nuv[46] = {-1,-1,1,1,1};
    // (31):  HCO + HO2 <=> CO2 + OH + H
    fwd_A[46]     = 30000000000000;
    fwd_beta[46]  = 0;
    fwd_Ea[46]    = 0;
    prefactor_units[46]  = 1.0000000000000002e-06;
    activation_units[46] = 0.50321666580471969;
    phase_units[46]      = pow(10,-12.000000);
    is_PD[46] = 0;
    nTB[46] = 0;

    // (32):  HCO + HCO <=> H2 + CO + CO
    kiv[47] = {13,13,1,11,11};
    nuv[47] = {-1,-1,1,1,1};
    // (32):  HCO + HCO <=> H2 + CO + CO
    fwd_A[47]     = 3000000000000;
    fwd_beta[47]  = 0;
    fwd_Ea[47]    = 0;
    prefactor_units[47]  = 1.0000000000000002e-06;
    activation_units[47] = 0.50321666580471969;
    phase_units[47]      = pow(10,-12.000000);
    is_PD[47] = 0;
    nTB[47] = 0;

    // (33):  HCO + CH3 <=> CO + CH4
    kiv[48] = {13,4,11,6};
    nuv[48] = {-1,-1,1,1};
    // (33):  HCO + CH3 <=> CO + CH4
    fwd_A[48]     = 26500000000000;
    fwd_beta[48]  = 0;
    fwd_Ea[48]    = 0;
    prefactor_units[48]  = 1.0000000000000002e-06;
    activation_units[48] = 0.50321666580471969;
    phase_units[48]      = pow(10,-12.000000);
    is_PD[48] = 0;
    nTB[48] = 0;

    // (34):  HCO + HCO <=> CH2O + CO
    kiv[49] = {13,13,15,11};
    nuv[49] = {-1,-1,1,1};
    // (34):  HCO + HCO <=> CH2O + CO
    fwd_A[49]     = 30000000000000;
    fwd_beta[49]  = 0;
    fwd_Ea[49]    = 0;
    prefactor_units[49]  = 1.0000000000000002e-06;
    activation_units[49] = 0.50321666580471969;
    phase_units[49]      = pow(10,-12.000000);
    is_PD[49] = 0;
    nTB[49] = 0;

    // (35):  CH2O + M <=> HCO + H + M
    kiv[16] = {15,13,0};
    nuv[16] = {-1,1,1};
    // (35):  CH2O + M <=> HCO + H + M
    fwd_A[16]     = 3.3000000000000002e+39;
    fwd_beta[16]  = -6.2999999999999998;
    fwd_Ea[16]    = 99900;
    prefactor_units[16]  = 1.0000000000000002e-06;
    activation_units[16] = 0.50321666580471969;
    phase_units[16]      = pow(10,-6.000000);
    is_PD[16] = 0;
    nTB[16] = 4;
    TB[16] = (amrex::Real *) malloc(4 * sizeof(amrex::Real));
    TBid[16] = (int *) malloc(4 * sizeof(int));
    TBid[16][0] = 1; TB[16][0] = 2.5; // H2
    TBid[16][1] = 8; TB[16][1] = 12; // H2O
    TBid[16][2] = 11; TB[16][2] = 1.8999999999999999; // CO
    TBid[16][3] = 22; TB[16][3] = 3.7999999999999998; // CO2

    // (36):  CH2O + M <=> CO + H2 + M
    kiv[17] = {15,11,1};
    nuv[17] = {-1,1,1};
    // (36):  CH2O + M <=> CO + H2 + M
    fwd_A[17]     = 3.0999999999999999e+45;
    fwd_beta[17]  = -8;
    fwd_Ea[17]    = 97510;
    prefactor_units[17]  = 1.0000000000000002e-06;
    activation_units[17] = 0.50321666580471969;
    phase_units[17]      = pow(10,-6.000000);
    is_PD[17] = 0;
    nTB[17] = 4;
    TB[17] = (amrex::Real *) malloc(4 * sizeof(amrex::Real));
    TBid[17] = (int *) malloc(4 * sizeof(int));
    TBid[17][0] = 1; TB[17][0] = 2.5; // H2
    TBid[17][1] = 8; TB[17][1] = 12; // H2O
    TBid[17][2] = 11; TB[17][2] = 1.8999999999999999; // CO
    TBid[17][3] = 22; TB[17][3] = 3.7999999999999998; // CO2

    // (37):  CH2O + H <=> HCO + H2
    kiv[50] = {15,0,13,1};
    nuv[50] = {-1,-1,1,1};
    // (37):  CH2O + H <=> HCO + H2
    fwd_A[50]     = 57400000;
    fwd_beta[50]  = 1.8999999999999999;
    fwd_Ea[50]    = 2748.5999999999999;
    prefactor_units[50]  = 1.0000000000000002e-06;
    activation_units[50] = 0.50321666580471969;
    phase_units[50]      = pow(10,-12.000000);
    is_PD[50] = 0;
    nTB[50] = 0;

    // (38):  CH2O + O <=> HCO + OH
    kiv[51] = {15,5,13,7};
    nuv[51] = {-1,-1,1,1};
    // (38):  CH2O + O <=> HCO + OH
    fwd_A[51]     = 18100000000000;
    fwd_beta[51]  = 0;
    fwd_Ea[51]    = 3080;
    prefactor_units[51]  = 1.0000000000000002e-06;
    activation_units[51] = 0.50321666580471969;
    phase_units[51]      = pow(10,-12.000000);
    is_PD[51] = 0;
    nTB[51] = 0;

    // (39):  CH2O + OH <=> HCO + H2O
    kiv[52] = {15,7,13,8};
    nuv[52] = {-1,-1,1,1};
    // (39):  CH2O + OH <=> HCO + H2O
    fwd_A[52]     = 3430000000;
    fwd_beta[52]  = 1.1799999999999999;
    fwd_Ea[52]    = -447;
    prefactor_units[52]  = 1.0000000000000002e-06;
    activation_units[52] = 0.50321666580471969;
    phase_units[52]      = pow(10,-12.000000);
    is_PD[52] = 0;
    nTB[52] = 0;

    // (40):  CH2O + O2 <=> HCO + HO2
    kiv[53] = {15,19,13,20};
    nuv[53] = {-1,-1,1,1};
    // (40):  CH2O + O2 <=> HCO + HO2
    fwd_A[53]     = 1230000;
    fwd_beta[53]  = 3;
    fwd_Ea[53]    = 52000;
    prefactor_units[53]  = 1.0000000000000002e-06;
    activation_units[53] = 0.50321666580471969;
    phase_units[53]      = pow(10,-12.000000);
    is_PD[53] = 0;
    nTB[53] = 0;

    // (41):  CH2O + HO2 <=> HCO + H2O2
    kiv[54] = {15,20,13,21};
    nuv[54] = {-1,-1,1,1};
    // (41):  CH2O + HO2 <=> HCO + H2O2
    fwd_A[54]     = 41100;
    fwd_beta[54]  = 2.5;
    fwd_Ea[54]    = 10210;
    prefactor_units[54]  = 1.0000000000000002e-06;
    activation_units[54] = 0.50321666580471969;
    phase_units[54]      = pow(10,-12.000000);
    is_PD[54] = 0;
    nTB[54] = 0;

    // (42):  CH2O + CH3 <=> HCO + CH4
    kiv[55] = {15,4,13,6};
    nuv[55] = {-1,-1,1,1};
    // (42):  CH2O + CH3 <=> HCO + CH4
    fwd_A[55]     = 3.636e-06;
    fwd_beta[55]  = 5.4199999999999999;
    fwd_Ea[55]    = 998;
    prefactor_units[55]  = 1.0000000000000002e-06;
    activation_units[55] = 0.50321666580471969;
    phase_units[55]      = pow(10,-12.000000);
    is_PD[55] = 0;
    nTB[55] = 0;

    // (43):  CH3 + O <=> CH2O + H
    kiv[56] = {4,5,15,0};
    nuv[56] = {-1,-1,1,1};
    // (43):  CH3 + O <=> CH2O + H
    fwd_A[56]     = 84300000000000;
    fwd_beta[56]  = 0;
    fwd_Ea[56]    = 0;
    prefactor_units[56]  = 1.0000000000000002e-06;
    activation_units[56] = 0.50321666580471969;
    phase_units[56]      = pow(10,-12.000000);
    is_PD[56] = 0;
    nTB[56] = 0;

    // (44):  CH3 + O2 <=> CH3O + O
    kiv[57] = {4,19,18,5};
    nuv[57] = {-1,-1,1,1};
    // (44):  CH3 + O2 <=> CH3O + O
    fwd_A[57]     = 1.99e+18;
    fwd_beta[57]  = -1.5700000000000001;
    fwd_Ea[57]    = 29230;
    prefactor_units[57]  = 1.0000000000000002e-06;
    activation_units[57] = 0.50321666580471969;
    phase_units[57]      = pow(10,-12.000000);
    is_PD[57] = 0;
    nTB[57] = 0;

    // (45):  CH3 + O2 <=> CH2O + OH
    kiv[58] = {4,19,15,7};
    nuv[58] = {-1,-1,1,1};
    // (45):  CH3 + O2 <=> CH2O + OH
    fwd_A[58]     = 374000000000;
    fwd_beta[58]  = 0;
    fwd_Ea[58]    = 14640;
    prefactor_units[58]  = 1.0000000000000002e-06;
    activation_units[58] = 0.50321666580471969;
    phase_units[58]      = pow(10,-12.000000);
    is_PD[58] = 0;
    nTB[58] = 0;

    // (46):  CH3 + HO2 <=> CH3O + OH
    kiv[59] = {4,20,18,7};
    nuv[59] = {-1,-1,1,1};
    // (46):  CH3 + HO2 <=> CH3O + OH
    fwd_A[59]     = 24100000000;
    fwd_beta[59]  = 0.76000000000000001;
    fwd_Ea[59]    = -2325;
    prefactor_units[59]  = 1.0000000000000002e-06;
    activation_units[59] = 0.50321666580471969;
    phase_units[59]      = pow(10,-12.000000);
    is_PD[59] = 0;
    nTB[59] = 0;

    // (47):  CH3 + CH3 (+M) <=> C2H6 (+M)
    kiv[2] = {4,4,16};
    nuv[2] = {-1,-1,1};
    // (47):  CH3 + CH3 (+M) <=> C2H6 (+M)
    fwd_A[2]     = 2277000000000000;
    fwd_beta[2]  = -0.68999999999999995;
    fwd_Ea[2]    = 174.86000000000001;
    low_A[2]     = 8.0540000000000002e+31;
    low_beta[2]  = -3.75;
    low_Ea[2]    = 981.60000000000002;
    troe_a[2]    = 0;
    troe_Tsss[2] = 570;
    troe_Ts[2]   = 0;
    troe_Tss[2]  = 1e+30;
    troe_len[2]  = 4;
    prefactor_units[2]  = 1.0000000000000002e-06;
    activation_units[2] = 0.50321666580471969;
    phase_units[2]      = pow(10,-12.000000);
    is_PD[2] = 1;
    nTB[2] = 3;
    TB[2] = (amrex::Real *) malloc(3 * sizeof(amrex::Real));
    TBid[2] = (int *) malloc(3 * sizeof(int));
    TBid[2][0] = 8; TB[2][0] = 5; // H2O
    TBid[2][1] = 11; TB[2][1] = 2; // CO
    TBid[2][2] = 22; TB[2][2] = 3; // CO2

    // (48):  CH3 + H (+M) <=> CH4 (+M)
    kiv[3] = {4,0,6};
    nuv[3] = {-1,-1,1};
    // (48):  CH3 + H (+M) <=> CH4 (+M)
    fwd_A[3]     = 12700000000000000;
    fwd_beta[3]  = -0.63;
    fwd_Ea[3]    = 383;
    low_A[3]     = 2.4769999999999999e+33;
    low_beta[3]  = -4.7599999999999998;
    low_Ea[3]    = 2440;
    troe_a[3]    = 0.78300000000000003;
    troe_Tsss[3] = 74;
    troe_Ts[3]   = 2941;
    troe_Tss[3]  = 6964;
    troe_len[3]  = 4;
    prefactor_units[3]  = 1.0000000000000002e-06;
    activation_units[3] = 0.50321666580471969;
    phase_units[3]      = pow(10,-12.000000);
    is_PD[3] = 1;
    nTB[3] = 6;
    TB[3] = (amrex::Real *) malloc(6 * sizeof(amrex::Real));
    TBid[3] = (int *) malloc(6 * sizeof(int));
    TBid[3][0] = 1; TB[3][0] = 2; // H2
    TBid[3][1] = 8; TB[3][1] = 6; // H2O
    TBid[3][2] = 6; TB[3][2] = 2; // CH4
    TBid[3][3] = 11; TB[3][3] = 1.5; // CO
    TBid[3][4] = 22; TB[3][4] = 2; // CO2
    TBid[3][5] = 16; TB[3][5] = 3; // C2H6

    // (49):  CH4 + H <=> CH3 + H2
    kiv[60] = {6,0,4,1};
    nuv[60] = {-1,-1,1,1};
    // (49):  CH4 + H <=> CH3 + H2
    fwd_A[60]     = 54700000;
    fwd_beta[60]  = 1.97;
    fwd_Ea[60]    = 11210;
    prefactor_units[60]  = 1.0000000000000002e-06;
    activation_units[60] = 0.50321666580471969;
    phase_units[60]      = pow(10,-12.000000);
    is_PD[60] = 0;
    nTB[60] = 0;

    // (50):  CH4 + O <=> CH3 + OH
    kiv[61] = {6,5,4,7};
    nuv[61] = {-1,-1,1,1};
    // (50):  CH4 + O <=> CH3 + OH
    fwd_A[61]     = 3150000000000;
    fwd_beta[61]  = 0.5;
    fwd_Ea[61]    = 10290;
    prefactor_units[61]  = 1.0000000000000002e-06;
    activation_units[61] = 0.50321666580471969;
    phase_units[61]      = pow(10,-12.000000);
    is_PD[61] = 0;
    nTB[61] = 0;

    // (51):  CH4 + OH <=> CH3 + H2O
    kiv[62] = {6,7,4,8};
    nuv[62] = {-1,-1,1,1};
    // (51):  CH4 + OH <=> CH3 + H2O
    fwd_A[62]     = 5720000;
    fwd_beta[62]  = 1.96;
    fwd_Ea[62]    = 2639;
    prefactor_units[62]  = 1.0000000000000002e-06;
    activation_units[62] = 0.50321666580471969;
    phase_units[62]      = pow(10,-12.000000);
    is_PD[62] = 0;
    nTB[62] = 0;

    // (52):  CH3 + HO2 <=> CH4 + O2
    kiv[63] = {4,20,6,19};
    nuv[63] = {-1,-1,1,1};
    // (52):  CH3 + HO2 <=> CH4 + O2
    fwd_A[63]     = 3160000000000;
    fwd_beta[63]  = 0;
    fwd_Ea[63]    = 0;
    prefactor_units[63]  = 1.0000000000000002e-06;
    activation_units[63] = 0.50321666580471969;
    phase_units[63]      = pow(10,-12.000000);
    is_PD[63] = 0;
    nTB[63] = 0;

    // (53):  CH4 + HO2 <=> CH3 + H2O2
    kiv[64] = {6,20,4,21};
    nuv[64] = {-1,-1,1,1};
    // (53):  CH4 + HO2 <=> CH3 + H2O2
    fwd_A[64]     = 181000000000;
    fwd_beta[64]  = 0;
    fwd_Ea[64]    = 18580;
    prefactor_units[64]  = 1.0000000000000002e-06;
    activation_units[64] = 0.50321666580471969;
    phase_units[64]      = pow(10,-12.000000);
    is_PD[64] = 0;
    nTB[64] = 0;

    // (54):  CH2OH + M <=> CH2O + H + M
    kiv[18] = {17,15,0};
    nuv[18] = {-1,1,1};
    // (54):  CH2OH + M <=> CH2O + H + M
    fwd_A[18]     = 100000000000000;
    fwd_beta[18]  = 0;
    fwd_Ea[18]    = 25100;
    prefactor_units[18]  = 1.0000000000000002e-06;
    activation_units[18] = 0.50321666580471969;
    phase_units[18]      = pow(10,-6.000000);
    is_PD[18] = 0;
    nTB[18] = 0;

    // (55):  CH2OH + H <=> CH2O + H2
    kiv[65] = {17,0,15,1};
    nuv[65] = {-1,-1,1,1};
    // (55):  CH2OH + H <=> CH2O + H2
    fwd_A[65]     = 6000000000000;
    fwd_beta[65]  = 0;
    fwd_Ea[65]    = 0;
    prefactor_units[65]  = 1.0000000000000002e-06;
    activation_units[65] = 0.50321666580471969;
    phase_units[65]      = pow(10,-12.000000);
    is_PD[65] = 0;
    nTB[65] = 0;

    // (56):  CH2OH + H <=> CH3 + OH
    kiv[66] = {17,0,4,7};
    nuv[66] = {-1,-1,1,1};
    // (56):  CH2OH + H <=> CH3 + OH
    fwd_A[66]     = 96350000000000;
    fwd_beta[66]  = 0;
    fwd_Ea[66]    = 0;
    prefactor_units[66]  = 1.0000000000000002e-06;
    activation_units[66] = 0.50321666580471969;
    phase_units[66]      = pow(10,-12.000000);
    is_PD[66] = 0;
    nTB[66] = 0;

    // (57):  CH2OH + O <=> CH2O + OH
    kiv[67] = {17,5,15,7};
    nuv[67] = {-1,-1,1,1};
    // (57):  CH2OH + O <=> CH2O + OH
    fwd_A[67]     = 42000000000000;
    fwd_beta[67]  = 0;
    fwd_Ea[67]    = 0;
    prefactor_units[67]  = 1.0000000000000002e-06;
    activation_units[67] = 0.50321666580471969;
    phase_units[67]      = pow(10,-12.000000);
    is_PD[67] = 0;
    nTB[67] = 0;

    // (58):  CH2OH + OH <=> CH2O + H2O
    kiv[68] = {17,7,15,8};
    nuv[68] = {-1,-1,1,1};
    // (58):  CH2OH + OH <=> CH2O + H2O
    fwd_A[68]     = 24000000000000;
    fwd_beta[68]  = 0;
    fwd_Ea[68]    = 0;
    prefactor_units[68]  = 1.0000000000000002e-06;
    activation_units[68] = 0.50321666580471969;
    phase_units[68]      = pow(10,-12.000000);
    is_PD[68] = 0;
    nTB[68] = 0;

    // (59):  CH2OH + O2 <=> CH2O + HO2
    kiv[69] = {17,19,15,20};
    nuv[69] = {-1,-1,1,1};
    // (59):  CH2OH + O2 <=> CH2O + HO2
    fwd_A[69]     = 241000000000000;
    fwd_beta[69]  = 0;
    fwd_Ea[69]    = 5017;
    prefactor_units[69]  = 1.0000000000000002e-06;
    activation_units[69] = 0.50321666580471969;
    phase_units[69]      = pow(10,-12.000000);
    is_PD[69] = 0;
    nTB[69] = 0;

    // (60):  CH2OH + O2 <=> CH2O + HO2
    kiv[70] = {17,19,15,20};
    nuv[70] = {-1,-1,1,1};
    // (60):  CH2OH + O2 <=> CH2O + HO2
    fwd_A[70]     = 1510000000000000;
    fwd_beta[70]  = -1;
    fwd_Ea[70]    = 0;
    prefactor_units[70]  = 1.0000000000000002e-06;
    activation_units[70] = 0.50321666580471969;
    phase_units[70]      = pow(10,-12.000000);
    is_PD[70] = 0;
    nTB[70] = 0;

    // (61):  CH2OH + HO2 <=> CH2O + H2O2
    kiv[71] = {17,20,15,21};
    nuv[71] = {-1,-1,1,1};
    // (61):  CH2OH + HO2 <=> CH2O + H2O2
    fwd_A[71]     = 12000000000000;
    fwd_beta[71]  = 0;
    fwd_Ea[71]    = 0;
    prefactor_units[71]  = 1.0000000000000002e-06;
    activation_units[71] = 0.50321666580471969;
    phase_units[71]      = pow(10,-12.000000);
    is_PD[71] = 0;
    nTB[71] = 0;

    // (62):  CH2OH + HCO <=> CH2O + CH2O
    kiv[72] = {17,13,15,15};
    nuv[72] = {-1,-1,1,1};
    // (62):  CH2OH + HCO <=> CH2O + CH2O
    fwd_A[72]     = 15000000000000;
    fwd_beta[72]  = 0;
    fwd_Ea[72]    = 0;
    prefactor_units[72]  = 1.0000000000000002e-06;
    activation_units[72] = 0.50321666580471969;
    phase_units[72]      = pow(10,-12.000000);
    is_PD[72] = 0;
    nTB[72] = 0;

    // (63):  CH3O + M <=> CH2O + H + M
    kiv[19] = {18,15,0};
    nuv[19] = {-1,1,1};
    // (63):  CH3O + M <=> CH2O + H + M
    fwd_A[19]     = 8.3e+17;
    fwd_beta[19]  = -1.2;
    fwd_Ea[19]    = 15500;
    prefactor_units[19]  = 1.0000000000000002e-06;
    activation_units[19] = 0.50321666580471969;
    phase_units[19]      = pow(10,-6.000000);
    is_PD[19] = 0;
    nTB[19] = 0;

    // (64):  CH3O + H <=> CH3 + OH
    kiv[73] = {18,0,4,7};
    nuv[73] = {-1,-1,1,1};
    // (64):  CH3O + H <=> CH3 + OH
    fwd_A[73]     = 32000000000000;
    fwd_beta[73]  = 0;
    fwd_Ea[73]    = 0;
    prefactor_units[73]  = 1.0000000000000002e-06;
    activation_units[73] = 0.50321666580471969;
    phase_units[73]      = pow(10,-12.000000);
    is_PD[73] = 0;
    nTB[73] = 0;

    // (65):  CH3O + O <=> CH2O + OH
    kiv[74] = {18,5,15,7};
    nuv[74] = {-1,-1,1,1};
    // (65):  CH3O + O <=> CH2O + OH
    fwd_A[74]     = 6000000000000;
    fwd_beta[74]  = 0;
    fwd_Ea[74]    = 0;
    prefactor_units[74]  = 1.0000000000000002e-06;
    activation_units[74] = 0.50321666580471969;
    phase_units[74]      = pow(10,-12.000000);
    is_PD[74] = 0;
    nTB[74] = 0;

    // (66):  CH3O + OH <=> CH2O + H2O
    kiv[75] = {18,7,15,8};
    nuv[75] = {-1,-1,1,1};
    // (66):  CH3O + OH <=> CH2O + H2O
    fwd_A[75]     = 18000000000000;
    fwd_beta[75]  = 0;
    fwd_Ea[75]    = 0;
    prefactor_units[75]  = 1.0000000000000002e-06;
    activation_units[75] = 0.50321666580471969;
    phase_units[75]      = pow(10,-12.000000);
    is_PD[75] = 0;
    nTB[75] = 0;

    // (67):  CH3O + O2 <=> CH2O + HO2
    kiv[76] = {18,19,15,20};
    nuv[76] = {-1,-1,1,1};
    // (67):  CH3O + O2 <=> CH2O + HO2
    fwd_A[76]     = 90330000000000;
    fwd_beta[76]  = 0;
    fwd_Ea[76]    = 11980;
    prefactor_units[76]  = 1.0000000000000002e-06;
    activation_units[76] = 0.50321666580471969;
    phase_units[76]      = pow(10,-12.000000);
    is_PD[76] = 0;
    nTB[76] = 0;

    // (68):  CH3O + O2 <=> CH2O + HO2
    kiv[77] = {18,19,15,20};
    nuv[77] = {-1,-1,1,1};
    // (68):  CH3O + O2 <=> CH2O + HO2
    fwd_A[77]     = 22000000000;
    fwd_beta[77]  = 0;
    fwd_Ea[77]    = 1748;
    prefactor_units[77]  = 1.0000000000000002e-06;
    activation_units[77] = 0.50321666580471969;
    phase_units[77]      = pow(10,-12.000000);
    is_PD[77] = 0;
    nTB[77] = 0;

    // (69):  CH3O + HO2 <=> CH2O + H2O2
    kiv[78] = {18,20,15,21};
    nuv[78] = {-1,-1,1,1};
    // (69):  CH3O + HO2 <=> CH2O + H2O2
    fwd_A[78]     = 300000000000;
    fwd_beta[78]  = 0;
    fwd_Ea[78]    = 0;
    prefactor_units[78]  = 1.0000000000000002e-06;
    activation_units[78] = 0.50321666580471969;
    phase_units[78]      = pow(10,-12.000000);
    is_PD[78] = 0;
    nTB[78] = 0;

    // (70):  CH3O + CO <=> CH3 + CO2
    kiv[79] = {18,11,4,22};
    nuv[79] = {-1,-1,1,1};
    // (70):  CH3O + CO <=> CH3 + CO2
    fwd_A[79]     = 16000000000000;
    fwd_beta[79]  = 0;
    fwd_Ea[79]    = 11800;
    prefactor_units[79]  = 1.0000000000000002e-06;
    activation_units[79] = 0.50321666580471969;
    phase_units[79]      = pow(10,-12.000000);
    is_PD[79] = 0;
    nTB[79] = 0;

    // (71):  CH3 + CH3 <=> H + C2H5
    kiv[80] = {4,4,0,14};
    nuv[80] = {-1,-1,1,1};
    // (71):  CH3 + CH3 <=> H + C2H5
    fwd_A[80]     = 4990000000000;
    fwd_beta[80]  = 0.10000000000000001;
    fwd_Ea[80]    = 10600;
    prefactor_units[80]  = 1.0000000000000002e-06;
    activation_units[80] = 0.50321666580471969;
    phase_units[80]      = pow(10,-12.000000);
    is_PD[80] = 0;
    nTB[80] = 0;

    // (72):  CH4 + CH2 <=> CH3 + CH3
    kiv[81] = {6,2,4,4};
    nuv[81] = {-1,-1,1,1};
    // (72):  CH4 + CH2 <=> CH3 + CH3
    fwd_A[81]     = 2460000;
    fwd_beta[81]  = 2;
    fwd_Ea[81]    = 8270;
    prefactor_units[81]  = 1.0000000000000002e-06;
    activation_units[81] = 0.50321666580471969;
    phase_units[81]      = pow(10,-12.000000);
    is_PD[81] = 0;
    nTB[81] = 0;

    // (73):  CH4 + CH2(S) <=> CH3 + CH3
    kiv[82] = {6,3,4,4};
    nuv[82] = {-1,-1,1,1};
    // (73):  CH4 + CH2(S) <=> CH3 + CH3
    fwd_A[82]     = 16000000000000;
    fwd_beta[82]  = 0;
    fwd_Ea[82]    = -570;
    prefactor_units[82]  = 1.0000000000000002e-06;
    activation_units[82] = 0.50321666580471969;
    phase_units[82]      = pow(10,-12.000000);
    is_PD[82] = 0;
    nTB[82] = 0;

    // (74):  CH3 + OH <=> CH2 + H2O
    kiv[83] = {4,7,2,8};
    nuv[83] = {-1,-1,1,1};
    // (74):  CH3 + OH <=> CH2 + H2O
    fwd_A[83]     = 56000000;
    fwd_beta[83]  = 1.6000000000000001;
    fwd_Ea[83]    = 5420;
    prefactor_units[83]  = 1.0000000000000002e-06;
    activation_units[83] = 0.50321666580471969;
    phase_units[83]      = pow(10,-12.000000);
    is_PD[83] = 0;
    nTB[83] = 0;

    // (75):  CH3 + OH <=> CH2(S) + H2O
    kiv[84] = {4,7,3,8};
    nuv[84] = {-1,-1,1,1};
    // (75):  CH3 + OH <=> CH2(S) + H2O
    fwd_A[84]     = 25010000000000;
    fwd_beta[84]  = 0;
    fwd_Ea[84]    = 0;
    prefactor_units[84]  = 1.0000000000000002e-06;
    activation_units[84] = 0.50321666580471969;
    phase_units[84]      = pow(10,-12.000000);
    is_PD[84] = 0;
    nTB[84] = 0;

    // (76):  CH3 + CH2 <=> C2H4 + H
    kiv[85] = {4,2,12,0};
    nuv[85] = {-1,-1,1,1};
    // (76):  CH3 + CH2 <=> C2H4 + H
    fwd_A[85]     = 40000000000000;
    fwd_beta[85]  = 0;
    fwd_Ea[85]    = 0;
    prefactor_units[85]  = 1.0000000000000002e-06;
    activation_units[85] = 0.50321666580471969;
    phase_units[85]      = pow(10,-12.000000);
    is_PD[85] = 0;
    nTB[85] = 0;

    // (77):  CH3 + CH2(S) <=> C2H4 + H
    kiv[86] = {4,3,12,0};
    nuv[86] = {-1,-1,1,1};
    // (77):  CH3 + CH2(S) <=> C2H4 + H
    fwd_A[86]     = 12000000000000;
    fwd_beta[86]  = 0;
    fwd_Ea[86]    = -570;
    prefactor_units[86]  = 1.0000000000000002e-06;
    activation_units[86] = 0.50321666580471969;
    phase_units[86]      = pow(10,-12.000000);
    is_PD[86] = 0;
    nTB[86] = 0;

    // (78):  CH3O + H <=> CH2(S) + H2O
    kiv[87] = {18,0,3,8};
    nuv[87] = {-1,-1,1,1};
    // (78):  CH3O + H <=> CH2(S) + H2O
    fwd_A[87]     = 16000000000000;
    fwd_beta[87]  = 0;
    fwd_Ea[87]    = 0;
    prefactor_units[87]  = 1.0000000000000002e-06;
    activation_units[87] = 0.50321666580471969;
    phase_units[87]      = pow(10,-12.000000);
    is_PD[87] = 0;
    nTB[87] = 0;

    // (79):  C2H6 + H <=> C2H5 + H2
    kiv[88] = {16,0,14,1};
    nuv[88] = {-1,-1,1,1};
    // (79):  C2H6 + H <=> C2H5 + H2
    fwd_A[88]     = 115000000;
    fwd_beta[88]  = 1.8999999999999999;
    fwd_Ea[88]    = 7530;
    prefactor_units[88]  = 1.0000000000000002e-06;
    activation_units[88] = 0.50321666580471969;
    phase_units[88]      = pow(10,-12.000000);
    is_PD[88] = 0;
    nTB[88] = 0;

    // (80):  C2H6 + O <=> C2H5 + OH
    kiv[89] = {16,5,14,7};
    nuv[89] = {-1,-1,1,1};
    // (80):  C2H6 + O <=> C2H5 + OH
    fwd_A[89]     = 89800000;
    fwd_beta[89]  = 1.9199999999999999;
    fwd_Ea[89]    = 5690;
    prefactor_units[89]  = 1.0000000000000002e-06;
    activation_units[89] = 0.50321666580471969;
    phase_units[89]      = pow(10,-12.000000);
    is_PD[89] = 0;
    nTB[89] = 0;

    // (81):  C2H6 + OH <=> C2H5 + H2O
    kiv[90] = {16,7,14,8};
    nuv[90] = {-1,-1,1,1};
    // (81):  C2H6 + OH <=> C2H5 + H2O
    fwd_A[90]     = 3540000;
    fwd_beta[90]  = 2.1200000000000001;
    fwd_Ea[90]    = 870;
    prefactor_units[90]  = 1.0000000000000002e-06;
    activation_units[90] = 0.50321666580471969;
    phase_units[90]      = pow(10,-12.000000);
    is_PD[90] = 0;
    nTB[90] = 0;

    // (82):  C2H6 + O2 <=> C2H5 + HO2
    kiv[91] = {16,19,14,20};
    nuv[91] = {-1,-1,1,1};
    // (82):  C2H6 + O2 <=> C2H5 + HO2
    fwd_A[91]     = 40000000000000;
    fwd_beta[91]  = 0;
    fwd_Ea[91]    = 50900;
    prefactor_units[91]  = 1.0000000000000002e-06;
    activation_units[91] = 0.50321666580471969;
    phase_units[91]      = pow(10,-12.000000);
    is_PD[91] = 0;
    nTB[91] = 0;

    // (83):  C2H6 + HO2 <=> C2H5 + H2O2
    kiv[92] = {16,20,14,21};
    nuv[92] = {-1,-1,1,1};
    // (83):  C2H6 + HO2 <=> C2H5 + H2O2
    fwd_A[92]     = 294000000000;
    fwd_beta[92]  = 0;
    fwd_Ea[92]    = 14940;
    prefactor_units[92]  = 1.0000000000000002e-06;
    activation_units[92] = 0.50321666580471969;
    phase_units[92]      = pow(10,-12.000000);
    is_PD[92] = 0;
    nTB[92] = 0;

    // (84):  C2H6 + CH3 <=> C2H5 + CH4
    kiv[93] = {16,4,14,6};
    nuv[93] = {-1,-1,1,1};
    // (84):  C2H6 + CH3 <=> C2H5 + CH4
    fwd_A[93]     = 6140000;
    fwd_beta[93]  = 1.74;
    fwd_Ea[93]    = 10450;
    prefactor_units[93]  = 1.0000000000000002e-06;
    activation_units[93] = 0.50321666580471969;
    phase_units[93]      = pow(10,-12.000000);
    is_PD[93] = 0;
    nTB[93] = 0;

    // (85):  C2H5 + H (+M) <=> C2H6 (+M)
    kiv[4] = {14,0,16};
    nuv[4] = {-1,-1,1};
    // (85):  C2H5 + H (+M) <=> C2H6 (+M)
    fwd_A[4]     = 5.21e+17;
    fwd_beta[4]  = -0.98999999999999999;
    fwd_Ea[4]    = 1580;
    low_A[4]     = 1.9900000000000001e+41;
    low_beta[4]  = -7.0800000000000001;
    low_Ea[4]    = 6685;
    troe_a[4]    = 0.84219999999999995;
    troe_Tsss[4] = 125;
    troe_Ts[4]   = 2219;
    troe_Tss[4]  = 6882;
    troe_len[4]  = 4;
    prefactor_units[4]  = 1.0000000000000002e-06;
    activation_units[4] = 0.50321666580471969;
    phase_units[4]      = pow(10,-12.000000);
    is_PD[4] = 1;
    nTB[4] = 6;
    TB[4] = (amrex::Real *) malloc(6 * sizeof(amrex::Real));
    TBid[4] = (int *) malloc(6 * sizeof(int));
    TBid[4][0] = 1; TB[4][0] = 2; // H2
    TBid[4][1] = 8; TB[4][1] = 6; // H2O
    TBid[4][2] = 6; TB[4][2] = 2; // CH4
    TBid[4][3] = 11; TB[4][3] = 1.5; // CO
    TBid[4][4] = 22; TB[4][4] = 2; // CO2
    TBid[4][5] = 16; TB[4][5] = 3; // C2H6

    // (86):  C2H5 + H <=> C2H4 + H2
    kiv[94] = {14,0,12,1};
    nuv[94] = {-1,-1,1,1};
    // (86):  C2H5 + H <=> C2H4 + H2
    fwd_A[94]     = 2000000000000;
    fwd_beta[94]  = 0;
    fwd_Ea[94]    = 0;
    prefactor_units[94]  = 1.0000000000000002e-06;
    activation_units[94] = 0.50321666580471969;
    phase_units[94]      = pow(10,-12.000000);
    is_PD[94] = 0;
    nTB[94] = 0;

    // (87):  C2H5 + O <=> CH3 + CH2O
    kiv[95] = {14,5,4,15};
    nuv[95] = {-1,-1,1,1};
    // (87):  C2H5 + O <=> CH3 + CH2O
    fwd_A[95]     = 132000000000000;
    fwd_beta[95]  = 0;
    fwd_Ea[95]    = 0;
    prefactor_units[95]  = 1.0000000000000002e-06;
    activation_units[95] = 0.50321666580471969;
    phase_units[95]      = pow(10,-12.000000);
    is_PD[95] = 0;
    nTB[95] = 0;

    // (88):  C2H5 + O2 <=> C2H4 + HO2
    kiv[96] = {14,19,12,20};
    nuv[96] = {-1,-1,1,1};
    // (88):  C2H5 + O2 <=> C2H4 + HO2
    fwd_A[96]     = 20000000000;
    fwd_beta[96]  = 0;
    fwd_Ea[96]    = 0;
    prefactor_units[96]  = 1.0000000000000002e-06;
    activation_units[96] = 0.50321666580471969;
    phase_units[96]      = pow(10,-12.000000);
    is_PD[96] = 0;
    nTB[96] = 0;

    // (89):  C2H5 + C2H5 <=> C2H4 + C2H6
    kiv[97] = {14,14,12,16};
    nuv[97] = {-1,-1,1,1};
    // (89):  C2H5 + C2H5 <=> C2H4 + C2H6
    fwd_A[97]     = 1400000000000;
    fwd_beta[97]  = 0;
    fwd_Ea[97]    = 0;
    prefactor_units[97]  = 1.0000000000000002e-06;
    activation_units[97] = 0.50321666580471969;
    phase_units[97]      = pow(10,-12.000000);
    is_PD[97] = 0;
    nTB[97] = 0;

    // (90):  C2H5 + HCO <=> C2H6 + CO
    kiv[98] = {14,13,16,11};
    nuv[98] = {-1,-1,1,1};
    // (90):  C2H5 + HCO <=> C2H6 + CO
    fwd_A[98]     = 120000000000000;
    fwd_beta[98]  = 0;
    fwd_Ea[98]    = 0;
    prefactor_units[98]  = 1.0000000000000002e-06;
    activation_units[98] = 0.50321666580471969;
    phase_units[98]      = pow(10,-12.000000);
    is_PD[98] = 0;
    nTB[98] = 0;

    // (91):  C2H5 + O <=> CH3HCO + H
    kiv[99] = {14,5,23,0};
    nuv[99] = {-1,-1,1,1};
    // (91):  C2H5 + O <=> CH3HCO + H
    fwd_A[99]     = 80200000000000;
    fwd_beta[99]  = 0;
    fwd_Ea[99]    = 0;
    prefactor_units[99]  = 1.0000000000000002e-06;
    activation_units[99] = 0.50321666580471969;
    phase_units[99]      = pow(10,-12.000000);
    is_PD[99] = 0;
    nTB[99] = 0;

    // (92):  C2H4 (+M) <=> H2 + C2H2 (+M)
    kiv[5] = {12,1,9};
    nuv[5] = {-1,1,1};
    // (92):  C2H4 (+M) <=> H2 + C2H2 (+M)
    fwd_A[5]     = 8000000000000;
    fwd_beta[5]  = 0.44;
    fwd_Ea[5]    = 88770;
    low_A[5]     = 7.0000000000000001e+50;
    low_beta[5]  = -9.3100000000000005;
    low_Ea[5]    = 99860;
    troe_a[5]    = 0.73450000000000004;
    troe_Tsss[5] = 180;
    troe_Ts[5]   = 1035;
    troe_Tss[5]  = 5417;
    troe_len[5]  = 4;
    prefactor_units[5]  = 1;
    activation_units[5] = 0.50321666580471969;
    phase_units[5]      = pow(10,-6.000000);
    is_PD[5] = 1;
    nTB[5] = 6;
    TB[5] = (amrex::Real *) malloc(6 * sizeof(amrex::Real));
    TBid[5] = (int *) malloc(6 * sizeof(int));
    TBid[5][0] = 1; TB[5][0] = 2; // H2
    TBid[5][1] = 8; TB[5][1] = 6; // H2O
    TBid[5][2] = 6; TB[5][2] = 2; // CH4
    TBid[5][3] = 11; TB[5][3] = 1.5; // CO
    TBid[5][4] = 22; TB[5][4] = 2; // CO2
    TBid[5][5] = 16; TB[5][5] = 3; // C2H6

    // (93):  C2H4 + H (+M) <=> C2H5 (+M)
    kiv[6] = {12,0,14};
    nuv[6] = {-1,-1,1};
    // (93):  C2H4 + H (+M) <=> C2H5 (+M)
    fwd_A[6]     = 1080000000000;
    fwd_beta[6]  = 0.45400000000000001;
    fwd_Ea[6]    = 1820;
    low_A[6]     = 1.1999999999999999e+42;
    low_beta[6]  = -7.6200000000000001;
    low_Ea[6]    = 6970;
    troe_a[6]    = 0.97529999999999994;
    troe_Tsss[6] = 210;
    troe_Ts[6]   = 984;
    troe_Tss[6]  = 4374;
    troe_len[6]  = 4;
    prefactor_units[6]  = 1.0000000000000002e-06;
    activation_units[6] = 0.50321666580471969;
    phase_units[6]      = pow(10,-12.000000);
    is_PD[6] = 1;
    nTB[6] = 6;
    TB[6] = (amrex::Real *) malloc(6 * sizeof(amrex::Real));
    TBid[6] = (int *) malloc(6 * sizeof(int));
    TBid[6][0] = 1; TB[6][0] = 2; // H2
    TBid[6][1] = 8; TB[6][1] = 6; // H2O
    TBid[6][2] = 6; TB[6][2] = 2; // CH4
    TBid[6][3] = 11; TB[6][3] = 1.5; // CO
    TBid[6][4] = 22; TB[6][4] = 2; // CO2
    TBid[6][5] = 16; TB[6][5] = 3; // C2H6

    // (94):  C2H3 + H (+M) <=> C2H4 (+M)
    kiv[7] = {10,0,12};
    nuv[7] = {-1,-1,1};
    // (94):  C2H3 + H (+M) <=> C2H4 (+M)
    fwd_A[7]     = 6080000000000;
    fwd_beta[7]  = 0.27000000000000002;
    fwd_Ea[7]    = 280;
    low_A[7]     = 1.3999999999999999e+30;
    low_beta[7]  = -3.8599999999999999;
    low_Ea[7]    = 3320;
    troe_a[7]    = 0.78200000000000003;
    troe_Tsss[7] = 207.5;
    troe_Ts[7]   = 2663;
    troe_Tss[7]  = 6095;
    troe_len[7]  = 4;
    prefactor_units[7]  = 1.0000000000000002e-06;
    activation_units[7] = 0.50321666580471969;
    phase_units[7]      = pow(10,-12.000000);
    is_PD[7] = 1;
    nTB[7] = 6;
    TB[7] = (amrex::Real *) malloc(6 * sizeof(amrex::Real));
    TBid[7] = (int *) malloc(6 * sizeof(int));
    TBid[7][0] = 1; TB[7][0] = 2; // H2
    TBid[7][1] = 8; TB[7][1] = 6; // H2O
    TBid[7][2] = 6; TB[7][2] = 2; // CH4
    TBid[7][3] = 11; TB[7][3] = 1.5; // CO
    TBid[7][4] = 22; TB[7][4] = 2; // CO2
    TBid[7][5] = 16; TB[7][5] = 3; // C2H6

    // (95):  C2H4 + H <=> C2H3 + H2
    kiv[100] = {12,0,10,1};
    nuv[100] = {-1,-1,1,1};
    // (95):  C2H4 + H <=> C2H3 + H2
    fwd_A[100]     = 1325000;
    fwd_beta[100]  = 2.5299999999999998;
    fwd_Ea[100]    = 12240;
    prefactor_units[100]  = 1.0000000000000002e-06;
    activation_units[100] = 0.50321666580471969;
    phase_units[100]      = pow(10,-12.000000);
    is_PD[100] = 0;
    nTB[100] = 0;

    // (96):  C2H4 + OH <=> C2H3 + H2O
    kiv[101] = {12,7,10,8};
    nuv[101] = {-1,-1,1,1};
    // (96):  C2H4 + OH <=> C2H3 + H2O
    fwd_A[101]     = 1800000;
    fwd_beta[101]  = 2;
    fwd_Ea[101]    = 2500;
    prefactor_units[101]  = 1.0000000000000002e-06;
    activation_units[101] = 0.50321666580471969;
    phase_units[101]      = pow(10,-12.000000);
    is_PD[101] = 0;
    nTB[101] = 0;

    // (97):  C2H4 + CH3 <=> C2H3 + CH4
    kiv[102] = {12,4,10,6};
    nuv[102] = {-1,-1,1,1};
    // (97):  C2H4 + CH3 <=> C2H3 + CH4
    fwd_A[102]     = 227000;
    fwd_beta[102]  = 2;
    fwd_Ea[102]    = 9200;
    prefactor_units[102]  = 1.0000000000000002e-06;
    activation_units[102] = 0.50321666580471969;
    phase_units[102]      = pow(10,-12.000000);
    is_PD[102] = 0;
    nTB[102] = 0;

    // (98):  C2H4 + O <=> CH3 + HCO
    kiv[103] = {12,5,4,13};
    nuv[103] = {-1,-1,1,1};
    // (98):  C2H4 + O <=> CH3 + HCO
    fwd_A[103]     = 19200000;
    fwd_beta[103]  = 1.8300000000000001;
    fwd_Ea[103]    = 220;
    prefactor_units[103]  = 1.0000000000000002e-06;
    activation_units[103] = 0.50321666580471969;
    phase_units[103]      = pow(10,-12.000000);
    is_PD[103] = 0;
    nTB[103] = 0;

    // (99):  C2H3 + OH <=> C2H2 + H2O
    kiv[104] = {10,7,9,8};
    nuv[104] = {-1,-1,1,1};
    // (99):  C2H3 + OH <=> C2H2 + H2O
    fwd_A[104]     = 5000000000000;
    fwd_beta[104]  = 0;
    fwd_Ea[104]    = 0;
    prefactor_units[104]  = 1.0000000000000002e-06;
    activation_units[104] = 0.50321666580471969;
    phase_units[104]      = pow(10,-12.000000);
    is_PD[104] = 0;
    nTB[104] = 0;

    // (100):  C2H4 + O <=> OH + C2H3
    kiv[105] = {12,5,7,10};
    nuv[105] = {-1,-1,1,1};
    // (100):  C2H4 + O <=> OH + C2H3
    fwd_A[105]     = 15100000;
    fwd_beta[105]  = 1.9099999999999999;
    fwd_Ea[105]    = 3740;
    prefactor_units[105]  = 1.0000000000000002e-06;
    activation_units[105] = 0.50321666580471969;
    phase_units[105]      = pow(10,-12.000000);
    is_PD[105] = 0;
    nTB[105] = 0;

    // (101):  C2H4 + O2 <=> C2H3 + HO2
    kiv[106] = {12,19,10,20};
    nuv[106] = {-1,-1,1,1};
    // (101):  C2H4 + O2 <=> C2H3 + HO2
    fwd_A[106]     = 42150000000000;
    fwd_beta[106]  = 0;
    fwd_Ea[106]    = 57600;
    prefactor_units[106]  = 1.0000000000000002e-06;
    activation_units[106] = 0.50321666580471969;
    phase_units[106]      = pow(10,-12.000000);
    is_PD[106] = 0;
    nTB[106] = 0;

    // (102):  C2H3 + H <=> C2H2 + H2
    kiv[107] = {10,0,9,1};
    nuv[107] = {-1,-1,1,1};
    // (102):  C2H3 + H <=> C2H2 + H2
    fwd_A[107]     = 96400000000000;
    fwd_beta[107]  = 0;
    fwd_Ea[107]    = 0;
    prefactor_units[107]  = 1.0000000000000002e-06;
    activation_units[107] = 0.50321666580471969;
    phase_units[107]      = pow(10,-12.000000);
    is_PD[107] = 0;
    nTB[107] = 0;

    // (103):  C2H3 + H2O2 <=> C2H4 + HO2
    kiv[108] = {10,21,12,20};
    nuv[108] = {-1,-1,1,1};
    // (103):  C2H3 + H2O2 <=> C2H4 + HO2
    fwd_A[108]     = 12100000000;
    fwd_beta[108]  = 0;
    fwd_Ea[108]    = -596;
    prefactor_units[108]  = 1.0000000000000002e-06;
    activation_units[108] = 0.50321666580471969;
    phase_units[108]      = pow(10,-12.000000);
    is_PD[108] = 0;
    nTB[108] = 0;

    // (104):  C2H3 + CH3 <=> C2H2 + CH4
    kiv[109] = {10,4,9,6};
    nuv[109] = {-1,-1,1,1};
    // (104):  C2H3 + CH3 <=> C2H2 + CH4
    fwd_A[109]     = 390000000000;
    fwd_beta[109]  = 0;
    fwd_Ea[109]    = 0;
    prefactor_units[109]  = 1.0000000000000002e-06;
    activation_units[109] = 0.50321666580471969;
    phase_units[109]      = pow(10,-12.000000);
    is_PD[109] = 0;
    nTB[109] = 0;

    // (105):  C2H3 + C2H3 <=> C2H4 + C2H2
    kiv[110] = {10,10,12,9};
    nuv[110] = {-1,-1,1,1};
    // (105):  C2H3 + C2H3 <=> C2H4 + C2H2
    fwd_A[110]     = 960000000000;
    fwd_beta[110]  = 0;
    fwd_Ea[110]    = 0;
    prefactor_units[110]  = 1.0000000000000002e-06;
    activation_units[110] = 0.50321666580471969;
    phase_units[110]      = pow(10,-12.000000);
    is_PD[110] = 0;
    nTB[110] = 0;

    // (106):  C2H3 + O2 <=> HCO + CH2O
    kiv[111] = {10,19,13,15};
    nuv[111] = {-1,-1,1,1};
    // (106):  C2H3 + O2 <=> HCO + CH2O
    fwd_A[111]     = 45800000000000000;
    fwd_beta[111]  = -1.3899999999999999;
    fwd_Ea[111]    = 1015;
    prefactor_units[111]  = 1.0000000000000002e-06;
    activation_units[111] = 0.50321666580471969;
    phase_units[111]      = pow(10,-12.000000);
    is_PD[111] = 0;
    nTB[111] = 0;

    // (107):  C2H3 + O2 <=> HO2 + C2H2
    kiv[112] = {10,19,20,9};
    nuv[112] = {-1,-1,1,1};
    // (107):  C2H3 + O2 <=> HO2 + C2H2
    fwd_A[112]     = 1337000;
    fwd_beta[112]  = 1.6100000000000001;
    fwd_Ea[112]    = -384;
    prefactor_units[112]  = 1.0000000000000002e-06;
    activation_units[112] = 0.50321666580471969;
    phase_units[112]      = pow(10,-12.000000);
    is_PD[112] = 0;
    nTB[112] = 0;

    // (108):  C2H2 + H (+M) <=> C2H3 (+M)
    kiv[8] = {9,0,10};
    nuv[8] = {-1,-1,1};
    // (108):  C2H2 + H (+M) <=> C2H3 (+M)
    fwd_A[8]     = 5600000000000;
    fwd_beta[8]  = 0;
    fwd_Ea[8]    = 2400;
    low_A[8]     = 3.8e+40;
    low_beta[8]  = -7.2699999999999996;
    low_Ea[8]    = 7220;
    troe_a[8]    = 0.75070000000000003;
    troe_Tsss[8] = 98.5;
    troe_Ts[8]   = 1302;
    troe_Tss[8]  = 4167;
    troe_len[8]  = 4;
    prefactor_units[8]  = 1.0000000000000002e-06;
    activation_units[8] = 0.50321666580471969;
    phase_units[8]      = pow(10,-12.000000);
    is_PD[8] = 1;
    nTB[8] = 6;
    TB[8] = (amrex::Real *) malloc(6 * sizeof(amrex::Real));
    TBid[8] = (int *) malloc(6 * sizeof(int));
    TBid[8][0] = 1; TB[8][0] = 2; // H2
    TBid[8][1] = 8; TB[8][1] = 6; // H2O
    TBid[8][2] = 6; TB[8][2] = 2; // CH4
    TBid[8][3] = 11; TB[8][3] = 1.5; // CO
    TBid[8][4] = 22; TB[8][4] = 2; // CO2
    TBid[8][5] = 16; TB[8][5] = 3; // C2H6

    // (109):  C2H2 + O <=> CH2 + CO
    kiv[113] = {9,5,2,11};
    nuv[113] = {-1,-1,1,1};
    // (109):  C2H2 + O <=> CH2 + CO
    fwd_A[113]     = 4080000;
    fwd_beta[113]  = 2;
    fwd_Ea[113]    = 1900;
    prefactor_units[113]  = 1.0000000000000002e-06;
    activation_units[113] = 0.50321666580471969;
    phase_units[113]      = pow(10,-12.000000);
    is_PD[113] = 0;
    nTB[113] = 0;

    // (110):  C2H2 + OH <=> CH3 + CO
    kiv[114] = {9,7,4,11};
    nuv[114] = {-1,-1,1,1};
    // (110):  C2H2 + OH <=> CH3 + CO
    fwd_A[114]     = 0.00048299999999999998;
    fwd_beta[114]  = 4;
    fwd_Ea[114]    = -2000;
    prefactor_units[114]  = 1.0000000000000002e-06;
    activation_units[114] = 0.50321666580471969;
    phase_units[114]      = pow(10,-12.000000);
    is_PD[114] = 0;
    nTB[114] = 0;

    // (111):  CH2 + H (+M) <=> CH3 (+M)
    kiv[9] = {2,0,4};
    nuv[9] = {-1,-1,1};
    // (111):  CH2 + H (+M) <=> CH3 (+M)
    fwd_A[9]     = 25000000000000000;
    fwd_beta[9]  = -0.80000000000000004;
    fwd_Ea[9]    = 0;
    low_A[9]     = 3.2000000000000002e+27;
    low_beta[9]  = -3.1400000000000001;
    low_Ea[9]    = 1230;
    troe_a[9]    = 0.68000000000000005;
    troe_Tsss[9] = 78;
    troe_Ts[9]   = 1995;
    troe_Tss[9]  = 5590;
    troe_len[9]  = 4;
    prefactor_units[9]  = 1.0000000000000002e-06;
    activation_units[9] = 0.50321666580471969;
    phase_units[9]      = pow(10,-12.000000);
    is_PD[9] = 1;
    nTB[9] = 6;
    TB[9] = (amrex::Real *) malloc(6 * sizeof(amrex::Real));
    TBid[9] = (int *) malloc(6 * sizeof(int));
    TBid[9][0] = 1; TB[9][0] = 2; // H2
    TBid[9][1] = 8; TB[9][1] = 6; // H2O
    TBid[9][2] = 6; TB[9][2] = 2; // CH4
    TBid[9][3] = 11; TB[9][3] = 1.5; // CO
    TBid[9][4] = 22; TB[9][4] = 2; // CO2
    TBid[9][5] = 16; TB[9][5] = 3; // C2H6

    // (112):  CH2 + O <=> HCO + H
    kiv[115] = {2,5,13,0};
    nuv[115] = {-1,-1,1,1};
    // (112):  CH2 + O <=> HCO + H
    fwd_A[115]     = 80000000000000;
    fwd_beta[115]  = 0;
    fwd_Ea[115]    = 0;
    prefactor_units[115]  = 1.0000000000000002e-06;
    activation_units[115] = 0.50321666580471969;
    phase_units[115]      = pow(10,-12.000000);
    is_PD[115] = 0;
    nTB[115] = 0;

    // (113):  CH2 + OH <=> CH2O + H
    kiv[116] = {2,7,15,0};
    nuv[116] = {-1,-1,1,1};
    // (113):  CH2 + OH <=> CH2O + H
    fwd_A[116]     = 20000000000000;
    fwd_beta[116]  = 0;
    fwd_Ea[116]    = 0;
    prefactor_units[116]  = 1.0000000000000002e-06;
    activation_units[116] = 0.50321666580471969;
    phase_units[116]      = pow(10,-12.000000);
    is_PD[116] = 0;
    nTB[116] = 0;

    // (114):  CH2 + H2 <=> H + CH3
    kiv[117] = {2,1,0,4};
    nuv[117] = {-1,-1,1,1};
    // (114):  CH2 + H2 <=> H + CH3
    fwd_A[117]     = 500000;
    fwd_beta[117]  = 2;
    fwd_Ea[117]    = 7230;
    prefactor_units[117]  = 1.0000000000000002e-06;
    activation_units[117] = 0.50321666580471969;
    phase_units[117]      = pow(10,-12.000000);
    is_PD[117] = 0;
    nTB[117] = 0;

    // (115):  CH2 + O2 <=> HCO + OH
    kiv[118] = {2,19,13,7};
    nuv[118] = {-1,-1,1,1};
    // (115):  CH2 + O2 <=> HCO + OH
    fwd_A[118]     = 13200000000000;
    fwd_beta[118]  = 0;
    fwd_Ea[118]    = 1500;
    prefactor_units[118]  = 1.0000000000000002e-06;
    activation_units[118] = 0.50321666580471969;
    phase_units[118]      = pow(10,-12.000000);
    is_PD[118] = 0;
    nTB[118] = 0;

    // (116):  CH2 + HO2 <=> CH2O + OH
    kiv[119] = {2,20,15,7};
    nuv[119] = {-1,-1,1,1};
    // (116):  CH2 + HO2 <=> CH2O + OH
    fwd_A[119]     = 20000000000000;
    fwd_beta[119]  = 0;
    fwd_Ea[119]    = 0;
    prefactor_units[119]  = 1.0000000000000002e-06;
    activation_units[119] = 0.50321666580471969;
    phase_units[119]      = pow(10,-12.000000);
    is_PD[119] = 0;
    nTB[119] = 0;

    // (117):  CH2 + CH2 <=> C2H2 + H2
    kiv[120] = {2,2,9,1};
    nuv[120] = {-1,-1,1,1};
    // (117):  CH2 + CH2 <=> C2H2 + H2
    fwd_A[120]     = 32000000000000;
    fwd_beta[120]  = 0;
    fwd_Ea[120]    = 0;
    prefactor_units[120]  = 1.0000000000000002e-06;
    activation_units[120] = 0.50321666580471969;
    phase_units[120]      = pow(10,-12.000000);
    is_PD[120] = 0;
    nTB[120] = 0;

    // (118):  CH2(S) + M <=> CH2 + M
    kiv[20] = {3,2};
    nuv[20] = {-1,1};
    // (118):  CH2(S) + M <=> CH2 + M
    fwd_A[20]     = 9000000000000;
    fwd_beta[20]  = 0;
    fwd_Ea[20]    = 600;
    prefactor_units[20]  = 1.0000000000000002e-06;
    activation_units[20] = 0.50321666580471969;
    phase_units[20]      = pow(10,-6.000000);
    is_PD[20] = 0;
    nTB[20] = 3;
    TB[20] = (amrex::Real *) malloc(3 * sizeof(amrex::Real));
    TBid[20] = (int *) malloc(3 * sizeof(int));
    TBid[20][0] = 8; TB[20][0] = 0; // H2O
    TBid[20][1] = 11; TB[20][1] = 0; // CO
    TBid[20][2] = 22; TB[20][2] = 0; // CO2

    // (119):  CH2(S) + H2O <=> CH2 + H2O
    kiv[121] = {3,8,2,8};
    nuv[121] = {-1,-1,1,1};
    // (119):  CH2(S) + H2O <=> CH2 + H2O
    fwd_A[121]     = 30000000000000;
    fwd_beta[121]  = 0;
    fwd_Ea[121]    = 0;
    prefactor_units[121]  = 1.0000000000000002e-06;
    activation_units[121] = 0.50321666580471969;
    phase_units[121]      = pow(10,-12.000000);
    is_PD[121] = 0;
    nTB[121] = 0;

    // (120):  CH2(S) + CO <=> CH2 + CO
    kiv[122] = {3,11,2,11};
    nuv[122] = {-1,-1,1,1};
    // (120):  CH2(S) + CO <=> CH2 + CO
    fwd_A[122]     = 9000000000000;
    fwd_beta[122]  = 0;
    fwd_Ea[122]    = 0;
    prefactor_units[122]  = 1.0000000000000002e-06;
    activation_units[122] = 0.50321666580471969;
    phase_units[122]      = pow(10,-12.000000);
    is_PD[122] = 0;
    nTB[122] = 0;

    // (121):  CH2(S) + CO2 <=> CH2 + CO2
    kiv[123] = {3,22,2,22};
    nuv[123] = {-1,-1,1,1};
    // (121):  CH2(S) + CO2 <=> CH2 + CO2
    fwd_A[123]     = 7000000000000;
    fwd_beta[123]  = 0;
    fwd_Ea[123]    = 0;
    prefactor_units[123]  = 1.0000000000000002e-06;
    activation_units[123] = 0.50321666580471969;
    phase_units[123]      = pow(10,-12.000000);
    is_PD[123] = 0;
    nTB[123] = 0;

    // (122):  CH2(S) + O <=> CO + H2
    kiv[124] = {3,5,11,1};
    nuv[124] = {-1,-1,1,1};
    // (122):  CH2(S) + O <=> CO + H2
    fwd_A[124]     = 15000000000000;
    fwd_beta[124]  = 0;
    fwd_Ea[124]    = 0;
    prefactor_units[124]  = 1.0000000000000002e-06;
    activation_units[124] = 0.50321666580471969;
    phase_units[124]      = pow(10,-12.000000);
    is_PD[124] = 0;
    nTB[124] = 0;

    // (123):  CH2(S) + O <=> HCO + H
    kiv[125] = {3,5,13,0};
    nuv[125] = {-1,-1,1,1};
    // (123):  CH2(S) + O <=> HCO + H
    fwd_A[125]     = 15000000000000;
    fwd_beta[125]  = 0;
    fwd_Ea[125]    = 0;
    prefactor_units[125]  = 1.0000000000000002e-06;
    activation_units[125] = 0.50321666580471969;
    phase_units[125]      = pow(10,-12.000000);
    is_PD[125] = 0;
    nTB[125] = 0;

    // (124):  CH2(S) + OH <=> CH2O + H
    kiv[126] = {3,7,15,0};
    nuv[126] = {-1,-1,1,1};
    // (124):  CH2(S) + OH <=> CH2O + H
    fwd_A[126]     = 30000000000000;
    fwd_beta[126]  = 0;
    fwd_Ea[126]    = 0;
    prefactor_units[126]  = 1.0000000000000002e-06;
    activation_units[126] = 0.50321666580471969;
    phase_units[126]      = pow(10,-12.000000);
    is_PD[126] = 0;
    nTB[126] = 0;

    // (125):  CH2(S) + H2 <=> CH3 + H
    kiv[127] = {3,1,4,0};
    nuv[127] = {-1,-1,1,1};
    // (125):  CH2(S) + H2 <=> CH3 + H
    fwd_A[127]     = 70000000000000;
    fwd_beta[127]  = 0;
    fwd_Ea[127]    = 0;
    prefactor_units[127]  = 1.0000000000000002e-06;
    activation_units[127] = 0.50321666580471969;
    phase_units[127]      = pow(10,-12.000000);
    is_PD[127] = 0;
    nTB[127] = 0;

    // (126):  CH2(S) + O2 <=> H + OH + CO
    kiv[128] = {3,19,0,7,11};
    nuv[128] = {-1,-1,1,1,1};
    // (126):  CH2(S) + O2 <=> H + OH + CO
    fwd_A[128]     = 28000000000000;
    fwd_beta[128]  = 0;
    fwd_Ea[128]    = 0;
    prefactor_units[128]  = 1.0000000000000002e-06;
    activation_units[128] = 0.50321666580471969;
    phase_units[128]      = pow(10,-12.000000);
    is_PD[128] = 0;
    nTB[128] = 0;

    // (127):  CH2(S) + O2 <=> CO + H2O
    kiv[129] = {3,19,11,8};
    nuv[129] = {-1,-1,1,1};
    // (127):  CH2(S) + O2 <=> CO + H2O
    fwd_A[129]     = 12000000000000;
    fwd_beta[129]  = 0;
    fwd_Ea[129]    = 0;
    prefactor_units[129]  = 1.0000000000000002e-06;
    activation_units[129] = 0.50321666580471969;
    phase_units[129]      = pow(10,-12.000000);
    is_PD[129] = 0;
    nTB[129] = 0;

    // (128):  CH2(S) + CO2 <=> CH2O + CO
    kiv[130] = {3,22,15,11};
    nuv[130] = {-1,-1,1,1};
    // (128):  CH2(S) + CO2 <=> CH2O + CO
    fwd_A[130]     = 14000000000000;
    fwd_beta[130]  = 0;
    fwd_Ea[130]    = 0;
    prefactor_units[130]  = 1.0000000000000002e-06;
    activation_units[130] = 0.50321666580471969;
    phase_units[130]      = pow(10,-12.000000);
    is_PD[130] = 0;
    nTB[130] = 0;

    // (129):  CH3HCO <=> CH3 + HCO
    kiv[131] = {23,4,13};
    nuv[131] = {-1,1,1};
    // (129):  CH3HCO <=> CH3 + HCO
    fwd_A[131]     = 7000000000000000;
    fwd_beta[131]  = 0;
    fwd_Ea[131]    = 81674;
    prefactor_units[131]  = 1;
    activation_units[131] = 0.50321666580471969;
    phase_units[131]      = pow(10,-6.000000);
    is_PD[131] = 0;
    nTB[131] = 0;

    // (130):  CH3OCH3 <=> CH3 + CH3O
    kiv[132] = {26,4,18};
    nuv[132] = {-1,1,1};
    // (130):  CH3OCH3 <=> CH3 + CH3O
    fwd_A[132]     = 1.6989700000000001e+42;
    fwd_beta[132]  = -7.9535900000000002;
    fwd_Ea[132]    = 91806.600000000006;
    prefactor_units[132]  = 1;
    activation_units[132] = 0.50321666580471969;
    phase_units[132]      = pow(10,-6.000000);
    is_PD[132] = 0;
    nTB[132] = 0;

    // (131):  CH3OCH3 + OH <=> CH3OCH2 + H2O
    kiv[133] = {26,7,24,8};
    nuv[133] = {-1,-1,1,1};
    // (131):  CH3OCH3 + OH <=> CH3OCH2 + H2O
    fwd_A[133]     = 6710000;
    fwd_beta[133]  = 2;
    fwd_Ea[133]    = -629.88;
    prefactor_units[133]  = 1.0000000000000002e-06;
    activation_units[133] = 0.50321666580471969;
    phase_units[133]      = pow(10,-12.000000);
    is_PD[133] = 0;
    nTB[133] = 0;

    // (132):  CH3OCH3 + H <=> CH3OCH2 + H2
    kiv[134] = {26,0,24,1};
    nuv[134] = {-1,-1,1,1};
    // (132):  CH3OCH3 + H <=> CH3OCH2 + H2
    fwd_A[134]     = 29700000;
    fwd_beta[134]  = 2;
    fwd_Ea[134]    = 4033.6100000000001;
    prefactor_units[134]  = 1.0000000000000002e-06;
    activation_units[134] = 0.50321666580471969;
    phase_units[134]      = pow(10,-12.000000);
    is_PD[134] = 0;
    nTB[134] = 0;

    // (133):  CH3OCH3 + CH3 <=> CH3OCH2 + CH4
    kiv[135] = {26,4,24,6};
    nuv[135] = {-1,-1,1,1};
    // (133):  CH3OCH3 + CH3 <=> CH3OCH2 + CH4
    fwd_A[135]     = 26.800000000000001;
    fwd_beta[135]  = 3.7778999999999998;
    fwd_Ea[135]    = 9631.2999999999993;
    prefactor_units[135]  = 1.0000000000000002e-06;
    activation_units[135] = 0.50321666580471969;
    phase_units[135]      = pow(10,-12.000000);
    is_PD[135] = 0;
    nTB[135] = 0;

    // (134):  CH3OCH3 + O <=> CH3OCH2 + OH
    kiv[136] = {26,5,24,7};
    nuv[136] = {-1,-1,1,1};
    // (134):  CH3OCH3 + O <=> CH3OCH2 + OH
    fwd_A[136]     = 0.0018550000000000001;
    fwd_beta[136]  = 5.29;
    fwd_Ea[136]    = -109;
    prefactor_units[136]  = 1.0000000000000002e-06;
    activation_units[136] = 0.50321666580471969;
    phase_units[136]      = pow(10,-12.000000);
    is_PD[136] = 0;
    nTB[136] = 0;

    // (135):  CH3OCH3 + HO2 <=> CH3OCH2 + H2O2
    kiv[137] = {26,20,24,21};
    nuv[137] = {-1,-1,1,1};
    // (135):  CH3OCH3 + HO2 <=> CH3OCH2 + H2O2
    fwd_A[137]     = 20000000000000;
    fwd_beta[137]  = 0;
    fwd_Ea[137]    = 16500;
    prefactor_units[137]  = 1.0000000000000002e-06;
    activation_units[137] = 0.50321666580471969;
    phase_units[137]      = pow(10,-12.000000);
    is_PD[137] = 0;
    nTB[137] = 0;

    // (136):  CH3OCH3 + O2 <=> CH3OCH2 + HO2
    kiv[138] = {26,19,24,20};
    nuv[138] = {-1,-1,1,1};
    // (136):  CH3OCH3 + O2 <=> CH3OCH2 + HO2
    fwd_A[138]     = 41000000000000;
    fwd_beta[138]  = 0;
    fwd_Ea[138]    = 44910;
    prefactor_units[138]  = 1.0000000000000002e-06;
    activation_units[138] = 0.50321666580471969;
    phase_units[138]      = pow(10,-12.000000);
    is_PD[138] = 0;
    nTB[138] = 0;

    // (137):  CH3OCH2 <=> CH2O + CH3
    kiv[139] = {24,15,4};
    nuv[139] = {-1,1,1};
    // (137):  CH3OCH2 <=> CH2O + CH3
    fwd_A[139]     = 12000000000000;
    fwd_beta[139]  = 0;
    fwd_Ea[139]    = 25750;
    prefactor_units[139]  = 1;
    activation_units[139] = 0.50321666580471969;
    phase_units[139]      = pow(10,-6.000000);
    is_PD[139] = 0;
    nTB[139] = 0;

    // (138):  CH3OCH2 + CH3O <=> CH3OCH3 + CH2O
    kiv[140] = {24,18,26,15};
    nuv[140] = {-1,-1,1,1};
    // (138):  CH3OCH2 + CH3O <=> CH3OCH3 + CH2O
    fwd_A[140]     = 24100000000000;
    fwd_beta[140]  = 0;
    fwd_Ea[140]    = 0;
    prefactor_units[140]  = 1.0000000000000002e-06;
    activation_units[140] = 0.50321666580471969;
    phase_units[140]      = pow(10,-12.000000);
    is_PD[140] = 0;
    nTB[140] = 0;

    // (139):  CH3OCH2 + CH2O <=> CH3OCH3 + HCO
    kiv[141] = {24,15,26,13};
    nuv[141] = {-1,-1,1,1};
    // (139):  CH3OCH2 + CH2O <=> CH3OCH3 + HCO
    fwd_A[141]     = 5490;
    fwd_beta[141]  = 2.7999999999999998;
    fwd_Ea[141]    = 5862;
    prefactor_units[141]  = 1.0000000000000002e-06;
    activation_units[141] = 0.50321666580471969;
    phase_units[141]      = pow(10,-12.000000);
    is_PD[141] = 0;
    nTB[141] = 0;

    // (140):  CH3OCH2 + HO2 <=> CH3OCH2O + OH
    kiv[142] = {24,20,30,7};
    nuv[142] = {-1,-1,1,1};
    // (140):  CH3OCH2 + HO2 <=> CH3OCH2O + OH
    fwd_A[142]     = 9000000000000;
    fwd_beta[142]  = 0;
    fwd_Ea[142]    = 0;
    prefactor_units[142]  = 1.0000000000000002e-06;
    activation_units[142] = 0.50321666580471969;
    phase_units[142]      = pow(10,-12.000000);
    is_PD[142] = 0;
    nTB[142] = 0;

    // (141):  CH3OCH2O <=> CH3OCHO + H
    kiv[143] = {30,29,0};
    nuv[143] = {-1,1,1};
    // (141):  CH3OCH2O <=> CH3OCHO + H
    fwd_A[143]     = 17450000000000000;
    fwd_beta[143]  = -0.66000000000000003;
    fwd_Ea[143]    = 11720;
    prefactor_units[143]  = 1;
    activation_units[143] = 0.50321666580471969;
    phase_units[143]      = pow(10,-6.000000);
    is_PD[143] = 0;
    nTB[143] = 0;

    // (142):  CH3OCHO + O2 <=> CH3OCO + HO2
    kiv[144] = {29,19,28,20};
    nuv[144] = {-1,-1,1,1};
    // (142):  CH3OCHO + O2 <=> CH3OCO + HO2
    fwd_A[144]     = 10000000000000;
    fwd_beta[144]  = 0;
    fwd_Ea[144]    = 49700;
    prefactor_units[144]  = 1.0000000000000002e-06;
    activation_units[144] = 0.50321666580471969;
    phase_units[144]      = pow(10,-12.000000);
    is_PD[144] = 0;
    nTB[144] = 0;

    // (143):  CH3OCHO + OH <=> CH3OCO + H2O
    kiv[145] = {29,7,28,8};
    nuv[145] = {-1,-1,1,1};
    // (143):  CH3OCHO + OH <=> CH3OCO + H2O
    fwd_A[145]     = 23400000;
    fwd_beta[145]  = 1.6100000000000001;
    fwd_Ea[145]    = -35;
    prefactor_units[145]  = 1.0000000000000002e-06;
    activation_units[145] = 0.50321666580471969;
    phase_units[145]      = pow(10,-12.000000);
    is_PD[145] = 0;
    nTB[145] = 0;

    // (144):  CH3OCHO + HO2 <=> CH3OCO + H2O2
    kiv[146] = {29,20,28,21};
    nuv[146] = {-1,-1,1,1};
    // (144):  CH3OCHO + HO2 <=> CH3OCO + H2O2
    fwd_A[146]     = 1220000000000;
    fwd_beta[146]  = 0;
    fwd_Ea[146]    = 17000;
    prefactor_units[146]  = 1.0000000000000002e-06;
    activation_units[146] = 0.50321666580471969;
    phase_units[146]      = pow(10,-12.000000);
    is_PD[146] = 0;
    nTB[146] = 0;

    // (145):  CH3OCHO + O <=> CH3OCO + OH
    kiv[147] = {29,5,28,7};
    nuv[147] = {-1,-1,1,1};
    // (145):  CH3OCHO + O <=> CH3OCO + OH
    fwd_A[147]     = 235000;
    fwd_beta[147]  = 2.5;
    fwd_Ea[147]    = 2230;
    prefactor_units[147]  = 1.0000000000000002e-06;
    activation_units[147] = 0.50321666580471969;
    phase_units[147]      = pow(10,-12.000000);
    is_PD[147] = 0;
    nTB[147] = 0;

    // (146):  CH3OCHO + H <=> CH3OCO + H2
    kiv[148] = {29,0,28,1};
    nuv[148] = {-1,-1,1,1};
    // (146):  CH3OCHO + H <=> CH3OCO + H2
    fwd_A[148]     = 4550000;
    fwd_beta[148]  = 2;
    fwd_Ea[148]    = 5000;
    prefactor_units[148]  = 1.0000000000000002e-06;
    activation_units[148] = 0.50321666580471969;
    phase_units[148]      = pow(10,-12.000000);
    is_PD[148] = 0;
    nTB[148] = 0;

    // (147):  CH3OCHO + CH3 <=> CH3OCO + CH4
    kiv[149] = {29,4,28,6};
    nuv[149] = {-1,-1,1,1};
    // (147):  CH3OCHO + CH3 <=> CH3OCO + CH4
    fwd_A[149]     = 0.755;
    fwd_beta[149]  = 3.46;
    fwd_Ea[149]    = 5481;
    prefactor_units[149]  = 1.0000000000000002e-06;
    activation_units[149] = 0.50321666580471969;
    phase_units[149]      = pow(10,-12.000000);
    is_PD[149] = 0;
    nTB[149] = 0;

    // (148):  CH3OCO <=> CH3O + CO
    kiv[150] = {28,18,11};
    nuv[150] = {-1,1,1};
    // (148):  CH3OCO <=> CH3O + CO
    fwd_A[150]     = 7451000000000;
    fwd_beta[150]  = -1.76;
    fwd_Ea[150]    = 17150;
    prefactor_units[150]  = 1;
    activation_units[150] = 0.50321666580471969;
    phase_units[150]      = pow(10,-6.000000);
    is_PD[150] = 0;
    nTB[150] = 0;

    // (149):  CH3OCO <=> CH3 + CO2
    kiv[151] = {28,4,22};
    nuv[151] = {-1,1,1};
    // (149):  CH3OCO <=> CH3 + CO2
    fwd_A[151]     = 1514000000000;
    fwd_beta[151]  = -1.78;
    fwd_Ea[151]    = 13820;
    prefactor_units[151]  = 1;
    activation_units[151] = 0.50321666580471969;
    phase_units[151]      = pow(10,-6.000000);
    is_PD[151] = 0;
    nTB[151] = 0;

    // (150):  CH3OCH2 + O2 <=> CH3OCH2O2
    kiv[152] = {24,19,34};
    nuv[152] = {-1,-1,1};
    // (150):  CH3OCH2 + O2 <=> CH3OCH2O2
    fwd_A[152]     = 2000000000000;
    fwd_beta[152]  = 0;
    fwd_Ea[152]    = 0;
    prefactor_units[152]  = 1.0000000000000002e-06;
    activation_units[152] = 0.50321666580471969;
    phase_units[152]      = pow(10,-12.000000);
    is_PD[152] = 0;
    nTB[152] = 0;

    // (151):  CH3OCH2O2 + CH3OCH2O2 <=> O2 + CH3OCH2O + CH3OCH2O
    kiv[153] = {34,34,19,30,30};
    nuv[153] = {-1,-1,1,1,1};
    // (151):  CH3OCH2O2 + CH3OCH2O2 <=> O2 + CH3OCH2O + CH3OCH2O
    fwd_A[153]     = 1.5969999999999999e+23;
    fwd_beta[153]  = -4.5;
    fwd_Ea[153]    = 0;
    prefactor_units[153]  = 1.0000000000000002e-06;
    activation_units[153] = 0.50321666580471969;
    phase_units[153]      = pow(10,-12.000000);
    is_PD[153] = 0;
    nTB[153] = 0;

    // (152):  CH3OCH2O2 + CH3OCH2O2 <=> O2 + CH3OCHO + CH3OCH2OH
    kiv[154] = {34,34,19,29,31};
    nuv[154] = {-1,-1,1,1,1};
    // (152):  CH3OCH2O2 + CH3OCH2O2 <=> O2 + CH3OCHO + CH3OCH2OH
    fwd_A[154]     = 6.8440000000000002e+22;
    fwd_beta[154]  = -4.5;
    fwd_Ea[154]    = 0;
    prefactor_units[154]  = 1.0000000000000002e-06;
    activation_units[154] = 0.50321666580471969;
    phase_units[154]      = pow(10,-12.000000);
    is_PD[154] = 0;
    nTB[154] = 0;

    // (153):  CH3OCH2O <=> CH3O + CH2O
    kiv[155] = {30,18,15};
    nuv[155] = {-1,1,1};
    // (153):  CH3OCH2O <=> CH3O + CH2O
    fwd_A[155]     = 9722000000000000;
    fwd_beta[155]  = -1.1000000000000001;
    fwd_Ea[155]    = 20640;
    prefactor_units[155]  = 1;
    activation_units[155] = 0.50321666580471969;
    phase_units[155]      = pow(10,-6.000000);
    is_PD[155] = 0;
    nTB[155] = 0;

    // (154):  CH3OCH2O + O2 <=> CH3OCHO + HO2
    kiv[156] = {30,19,29,20};
    nuv[156] = {-1,-1,1,1};
    // (154):  CH3OCH2O + O2 <=> CH3OCHO + HO2
    fwd_A[156]     = 50000000000;
    fwd_beta[156]  = 0;
    fwd_Ea[156]    = 500;
    prefactor_units[156]  = 1.0000000000000002e-06;
    activation_units[156] = 0.50321666580471969;
    phase_units[156]      = pow(10,-12.000000);
    is_PD[156] = 0;
    nTB[156] = 0;

    // (155):  CH3OCH2O2 <=> CH2OCH2O2H
    kiv[157] = {34,35};
    nuv[157] = {-1,1};
    // (155):  CH3OCH2O2 <=> CH2OCH2O2H
    fwd_A[157]     = 60000000000;
    fwd_beta[157]  = 0;
    fwd_Ea[157]    = 21500;
    prefactor_units[157]  = 1;
    activation_units[157] = 0.50321666580471969;
    phase_units[157]      = pow(10,-6.000000);
    is_PD[157] = 0;
    nTB[157] = 0;

    // (156):  CH2OCH2O2H <=> OH + CH2O + CH2O
    kiv[158] = {35,7,15,15};
    nuv[158] = {-1,1,1,1};
    // (156):  CH2OCH2O2H <=> OH + CH2O + CH2O
    fwd_A[158]     = 15000000000000;
    fwd_beta[158]  = 0;
    fwd_Ea[158]    = 20500;
    prefactor_units[158]  = 1;
    activation_units[158] = 0.50321666580471969;
    phase_units[158]      = pow(10,-6.000000);
    is_PD[158] = 0;
    nTB[158] = 0;

    // (157):  CH2OCH2O2H + O2 <=> O2CH2OCH2O2H
    kiv[159] = {35,19,37};
    nuv[159] = {-1,-1,1};
    // (157):  CH2OCH2O2H + O2 <=> O2CH2OCH2O2H
    fwd_A[159]     = 700000000000;
    fwd_beta[159]  = 0;
    fwd_Ea[159]    = 0;
    prefactor_units[159]  = 1.0000000000000002e-06;
    activation_units[159] = 0.50321666580471969;
    phase_units[159]      = pow(10,-12.000000);
    is_PD[159] = 0;
    nTB[159] = 0;

    // (158):  O2CH2OCH2O2H <=> HO2CH2OCHO + OH
    kiv[160] = {37,36,7};
    nuv[160] = {-1,1,1};
    // (158):  O2CH2OCH2O2H <=> HO2CH2OCHO + OH
    fwd_A[160]     = 40000000000;
    fwd_beta[160]  = 0;
    fwd_Ea[160]    = 18500;
    prefactor_units[160]  = 1;
    activation_units[160] = 0.50321666580471969;
    phase_units[160]      = pow(10,-6.000000);
    is_PD[160] = 0;
    nTB[160] = 0;

    // (159):  HO2CH2OCHO <=> OCH2OCHO + OH
    kiv[161] = {36,32,7};
    nuv[161] = {-1,1,1};
    // (159):  HO2CH2OCHO <=> OCH2OCHO + OH
    fwd_A[161]     = 30000000000000000;
    fwd_beta[161]  = 0;
    fwd_Ea[161]    = 40000;
    prefactor_units[161]  = 1;
    activation_units[161] = 0.50321666580471969;
    phase_units[161]      = pow(10,-6.000000);
    is_PD[161] = 0;
    nTB[161] = 0;

    // (160):  OCH2OCHO <=> HOCH2OCO
    kiv[162] = {32,33};
    nuv[162] = {-1,1};
    // (160):  OCH2OCHO <=> HOCH2OCO
    fwd_A[162]     = 100000000000;
    fwd_beta[162]  = 0;
    fwd_Ea[162]    = 14000;
    prefactor_units[162]  = 1;
    activation_units[162] = 0.50321666580471969;
    phase_units[162]      = pow(10,-6.000000);
    is_PD[162] = 0;
    nTB[162] = 0;

    // (161):  HOCH2OCO <=> HOCH2O + CO
    kiv[163] = {33,27,11};
    nuv[163] = {-1,1,1};
    // (161):  HOCH2OCO <=> HOCH2O + CO
    fwd_A[163]     = 21770000000000000;
    fwd_beta[163]  = -2.6899999999999999;
    fwd_Ea[163]    = 17200;
    prefactor_units[163]  = 1;
    activation_units[163] = 0.50321666580471969;
    phase_units[163]      = pow(10,-6.000000);
    is_PD[163] = 0;
    nTB[163] = 0;

    // (162):  HOCH2OCO <=> CH2OH + CO2
    kiv[164] = {33,17,22};
    nuv[164] = {-1,1,1};
    // (162):  HOCH2OCO <=> CH2OH + CO2
    fwd_A[164]     = 5311000000000000;
    fwd_beta[164]  = -2.6099999999999999;
    fwd_Ea[164]    = 20810;
    prefactor_units[164]  = 1;
    activation_units[164] = 0.50321666580471969;
    phase_units[164]      = pow(10,-6.000000);
    is_PD[164] = 0;
    nTB[164] = 0;

    // (163):  HOCH2O <=> HCOOH + H
    kiv[165] = {27,25,0};
    nuv[165] = {-1,1,1};
    // (163):  HOCH2O <=> HCOOH + H
    fwd_A[165]     = 100000000000000;
    fwd_beta[165]  = 0;
    fwd_Ea[165]    = 14900;
    prefactor_units[165]  = 1;
    activation_units[165] = 0.50321666580471969;
    phase_units[165]      = pow(10,-6.000000);
    is_PD[165] = 0;
    nTB[165] = 0;

    // (164):  CH2O + OH <=> HOCH2O
    kiv[166] = {15,7,27};
    nuv[166] = {-1,-1,1};
    // (164):  CH2O + OH <=> HOCH2O
    fwd_A[166]     = 4500000000000000;
    fwd_beta[166]  = -1.1100000000000001;
    fwd_Ea[166]    = 0;
    prefactor_units[166]  = 1.0000000000000002e-06;
    activation_units[166] = 0.50321666580471969;
    phase_units[166]      = pow(10,-12.000000);
    is_PD[166] = 0;
    nTB[166] = 0;

    // (165):  HCOOH + M <=> CO + H2O + M
    kiv[21] = {25,11,8};
    nuv[21] = {-1,1,1};
    // (165):  HCOOH + M <=> CO + H2O + M
    fwd_A[21]     = 23000000000000;
    fwd_beta[21]  = 0;
    fwd_Ea[21]    = 50000;
    prefactor_units[21]  = 1.0000000000000002e-06;
    activation_units[21] = 0.50321666580471969;
    phase_units[21]      = pow(10,-6.000000);
    is_PD[21] = 0;
    nTB[21] = 0;

    // (166):  HCOOH + M <=> CO2 + H2 + M
    kiv[22] = {25,22,1};
    nuv[22] = {-1,1,1};
    // (166):  HCOOH + M <=> CO2 + H2 + M
    fwd_A[22]     = 15000000000000000;
    fwd_beta[22]  = 0;
    fwd_Ea[22]    = 57000;
    prefactor_units[22]  = 1.0000000000000002e-06;
    activation_units[22] = 0.50321666580471969;
    phase_units[22]      = pow(10,-6.000000);
    is_PD[22] = 0;
    nTB[22] = 0;

    // (167):  HCOOH <=> HCO + OH
    kiv[167] = {25,13,7};
    nuv[167] = {-1,1,1};
    // (167):  HCOOH <=> HCO + OH
    fwd_A[167]     = 4.593e+18;
    fwd_beta[167]  = -0.46000000000000002;
    fwd_Ea[167]    = 108300;
    prefactor_units[167]  = 1;
    activation_units[167] = 0.50321666580471969;
    phase_units[167]      = pow(10,-6.000000);
    is_PD[167] = 0;
    nTB[167] = 0;

    // (168):  HCOOH + OH <=> H2O + CO2 + H
    kiv[168] = {25,7,8,22,0};
    nuv[168] = {-1,-1,1,1,1};
    // (168):  HCOOH + OH <=> H2O + CO2 + H
    fwd_A[168]     = 2620000;
    fwd_beta[168]  = 2.0600000000000001;
    fwd_Ea[168]    = 916;
    prefactor_units[168]  = 1.0000000000000002e-06;
    activation_units[168] = 0.50321666580471969;
    phase_units[168]      = pow(10,-12.000000);
    is_PD[168] = 0;
    nTB[168] = 0;

    // (169):  HCOOH + OH <=> H2O + CO + OH
    kiv[169] = {25,7,8,11,7};
    nuv[169] = {-1,-1,1,1,1};
    // (169):  HCOOH + OH <=> H2O + CO + OH
    fwd_A[169]     = 18500000;
    fwd_beta[169]  = 1.51;
    fwd_Ea[169]    = -962;
    prefactor_units[169]  = 1.0000000000000002e-06;
    activation_units[169] = 0.50321666580471969;
    phase_units[169]      = pow(10,-12.000000);
    is_PD[169] = 0;
    nTB[169] = 0;

    // (170):  HCOOH + H <=> H2 + CO2 + H
    kiv[170] = {25,0,1,22,0};
    nuv[170] = {-1,-1,1,1,1};
    // (170):  HCOOH + H <=> H2 + CO2 + H
    fwd_A[170]     = 4240000;
    fwd_beta[170]  = 2.1000000000000001;
    fwd_Ea[170]    = 4868;
    prefactor_units[170]  = 1.0000000000000002e-06;
    activation_units[170] = 0.50321666580471969;
    phase_units[170]      = pow(10,-12.000000);
    is_PD[170] = 0;
    nTB[170] = 0;

    // (171):  HCOOH + H <=> H2 + CO + OH
    kiv[171] = {25,0,1,11,7};
    nuv[171] = {-1,-1,1,1,1};
    // (171):  HCOOH + H <=> H2 + CO + OH
    fwd_A[171]     = 60300000000000;
    fwd_beta[171]  = -0.34999999999999998;
    fwd_Ea[171]    = 2988;
    prefactor_units[171]  = 1.0000000000000002e-06;
    activation_units[171] = 0.50321666580471969;
    phase_units[171]      = pow(10,-12.000000);
    is_PD[171] = 0;
    nTB[171] = 0;

    // (172):  HCOOH + CH3 <=> CH4 + CO + OH
    kiv[172] = {25,4,6,11,7};
    nuv[172] = {-1,-1,1,1,1};
    // (172):  HCOOH + CH3 <=> CH4 + CO + OH
    fwd_A[172]     = 3.9000000000000002e-07;
    fwd_beta[172]  = 5.7999999999999998;
    fwd_Ea[172]    = 2200;
    prefactor_units[172]  = 1.0000000000000002e-06;
    activation_units[172] = 0.50321666580471969;
    phase_units[172]      = pow(10,-12.000000);
    is_PD[172] = 0;
    nTB[172] = 0;

    // (173):  HCOOH + HO2 <=> H2O2 + CO + OH
    kiv[173] = {25,20,21,11,7};
    nuv[173] = {-1,-1,1,1,1};
    // (173):  HCOOH + HO2 <=> H2O2 + CO + OH
    fwd_A[173]     = 1000000000000;
    fwd_beta[173]  = 0;
    fwd_Ea[173]    = 11920;
    prefactor_units[173]  = 1.0000000000000002e-06;
    activation_units[173] = 0.50321666580471969;
    phase_units[173]      = pow(10,-12.000000);
    is_PD[173] = 0;
    nTB[173] = 0;

    // (174):  HCOOH + O <=> CO + OH + OH
    kiv[174] = {25,5,11,7,7};
    nuv[174] = {-1,-1,1,1,1};
    // (174):  HCOOH + O <=> CO + OH + OH
    fwd_A[174]     = 1.77e+18;
    fwd_beta[174]  = -1.8999999999999999;
    fwd_Ea[174]    = 2975;
    prefactor_units[174]  = 1.0000000000000002e-06;
    activation_units[174] = 0.50321666580471969;
    phase_units[174]      = pow(10,-12.000000);
    is_PD[174] = 0;
    nTB[174] = 0;

}


/* Finalizes parameter database */
void CKFINALIZE()
{
  for (int i=0; i<175; ++i) {
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
    for (id = 0; id < kd * 39; ++ id) {
         ncf[id] = 0; 
    }

    /*H */
    ncf[ 0 * kd + 1 ] = 1; /*H */

    /*H2 */
    ncf[ 1 * kd + 1 ] = 2; /*H */

    /*CH2 */
    ncf[ 2 * kd + 0 ] = 1; /*C */
    ncf[ 2 * kd + 1 ] = 2; /*H */

    /*CH2(S) */
    ncf[ 3 * kd + 0 ] = 1; /*C */
    ncf[ 3 * kd + 1 ] = 2; /*H */

    /*CH3 */
    ncf[ 4 * kd + 0 ] = 1; /*C */
    ncf[ 4 * kd + 1 ] = 3; /*H */

    /*O */
    ncf[ 5 * kd + 2 ] = 1; /*O */

    /*CH4 */
    ncf[ 6 * kd + 0 ] = 1; /*C */
    ncf[ 6 * kd + 1 ] = 4; /*H */

    /*OH */
    ncf[ 7 * kd + 2 ] = 1; /*O */
    ncf[ 7 * kd + 1 ] = 1; /*H */

    /*H2O */
    ncf[ 8 * kd + 1 ] = 2; /*H */
    ncf[ 8 * kd + 2 ] = 1; /*O */

    /*C2H2 */
    ncf[ 9 * kd + 0 ] = 2; /*C */
    ncf[ 9 * kd + 1 ] = 2; /*H */

    /*C2H3 */
    ncf[ 10 * kd + 0 ] = 2; /*C */
    ncf[ 10 * kd + 1 ] = 3; /*H */

    /*CO */
    ncf[ 11 * kd + 0 ] = 1; /*C */
    ncf[ 11 * kd + 2 ] = 1; /*O */

    /*C2H4 */
    ncf[ 12 * kd + 0 ] = 2; /*C */
    ncf[ 12 * kd + 1 ] = 4; /*H */

    /*HCO */
    ncf[ 13 * kd + 1 ] = 1; /*H */
    ncf[ 13 * kd + 0 ] = 1; /*C */
    ncf[ 13 * kd + 2 ] = 1; /*O */

    /*C2H5 */
    ncf[ 14 * kd + 0 ] = 2; /*C */
    ncf[ 14 * kd + 1 ] = 5; /*H */

    /*CH2O */
    ncf[ 15 * kd + 0 ] = 1; /*C */
    ncf[ 15 * kd + 1 ] = 2; /*H */
    ncf[ 15 * kd + 2 ] = 1; /*O */

    /*C2H6 */
    ncf[ 16 * kd + 0 ] = 2; /*C */
    ncf[ 16 * kd + 1 ] = 6; /*H */

    /*CH2OH */
    ncf[ 17 * kd + 0 ] = 1; /*C */
    ncf[ 17 * kd + 1 ] = 3; /*H */
    ncf[ 17 * kd + 2 ] = 1; /*O */

    /*CH3O */
    ncf[ 18 * kd + 0 ] = 1; /*C */
    ncf[ 18 * kd + 1 ] = 3; /*H */
    ncf[ 18 * kd + 2 ] = 1; /*O */

    /*O2 */
    ncf[ 19 * kd + 2 ] = 2; /*O */

    /*HO2 */
    ncf[ 20 * kd + 1 ] = 1; /*H */
    ncf[ 20 * kd + 2 ] = 2; /*O */

    /*H2O2 */
    ncf[ 21 * kd + 1 ] = 2; /*H */
    ncf[ 21 * kd + 2 ] = 2; /*O */

    /*CO2 */
    ncf[ 22 * kd + 0 ] = 1; /*C */
    ncf[ 22 * kd + 2 ] = 2; /*O */

    /*CH3HCO */
    ncf[ 23 * kd + 0 ] = 2; /*C */
    ncf[ 23 * kd + 1 ] = 4; /*H */
    ncf[ 23 * kd + 2 ] = 1; /*O */

    /*CH3OCH2 */
    ncf[ 24 * kd + 0 ] = 2; /*C */
    ncf[ 24 * kd + 1 ] = 5; /*H */
    ncf[ 24 * kd + 2 ] = 1; /*O */

    /*HCOOH */
    ncf[ 25 * kd + 0 ] = 1; /*C */
    ncf[ 25 * kd + 1 ] = 2; /*H */
    ncf[ 25 * kd + 2 ] = 2; /*O */

    /*CH3OCH3 */
    ncf[ 26 * kd + 0 ] = 2; /*C */
    ncf[ 26 * kd + 1 ] = 6; /*H */
    ncf[ 26 * kd + 2 ] = 1; /*O */

    /*HOCH2O */
    ncf[ 27 * kd + 0 ] = 1; /*C */
    ncf[ 27 * kd + 1 ] = 3; /*H */
    ncf[ 27 * kd + 2 ] = 2; /*O */

    /*CH3OCO */
    ncf[ 28 * kd + 0 ] = 2; /*C */
    ncf[ 28 * kd + 1 ] = 3; /*H */
    ncf[ 28 * kd + 2 ] = 2; /*O */

    /*CH3OCHO */
    ncf[ 29 * kd + 0 ] = 2; /*C */
    ncf[ 29 * kd + 1 ] = 4; /*H */
    ncf[ 29 * kd + 2 ] = 2; /*O */

    /*CH3OCH2O */
    ncf[ 30 * kd + 0 ] = 2; /*C */
    ncf[ 30 * kd + 1 ] = 5; /*H */
    ncf[ 30 * kd + 2 ] = 2; /*O */

    /*CH3OCH2OH */
    ncf[ 31 * kd + 0 ] = 2; /*C */
    ncf[ 31 * kd + 1 ] = 6; /*H */
    ncf[ 31 * kd + 2 ] = 2; /*O */

    /*OCH2OCHO */
    ncf[ 32 * kd + 0 ] = 2; /*C */
    ncf[ 32 * kd + 1 ] = 3; /*H */
    ncf[ 32 * kd + 2 ] = 3; /*O */

    /*HOCH2OCO */
    ncf[ 33 * kd + 0 ] = 2; /*C */
    ncf[ 33 * kd + 1 ] = 3; /*H */
    ncf[ 33 * kd + 2 ] = 3; /*O */

    /*CH3OCH2O2 */
    ncf[ 34 * kd + 0 ] = 2; /*C */
    ncf[ 34 * kd + 1 ] = 5; /*H */
    ncf[ 34 * kd + 2 ] = 3; /*O */

    /*CH2OCH2O2H */
    ncf[ 35 * kd + 0 ] = 2; /*C */
    ncf[ 35 * kd + 1 ] = 5; /*H */
    ncf[ 35 * kd + 2 ] = 3; /*O */

    /*HO2CH2OCHO */
    ncf[ 36 * kd + 0 ] = 2; /*C */
    ncf[ 36 * kd + 1 ] = 4; /*H */
    ncf[ 36 * kd + 2 ] = 4; /*O */

    /*O2CH2OCH2O2H */
    ncf[ 37 * kd + 0 ] = 2; /*C */
    ncf[ 37 * kd + 1 ] = 5; /*H */
    ncf[ 37 * kd + 2 ] = 5; /*O */

    /*N2 */
    ncf[ 38 * kd + 3 ] = 2; /*N */

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
    kname.resize(39);
    kname[0] = "H";
    kname[1] = "H2";
    kname[2] = "CH2";
    kname[3] = "CH2(S)";
    kname[4] = "CH3";
    kname[5] = "O";
    kname[6] = "CH4";
    kname[7] = "OH";
    kname[8] = "H2O";
    kname[9] = "C2H2";
    kname[10] = "C2H3";
    kname[11] = "CO";
    kname[12] = "C2H4";
    kname[13] = "HCO";
    kname[14] = "C2H5";
    kname[15] = "CH2O";
    kname[16] = "C2H6";
    kname[17] = "CH2OH";
    kname[18] = "CH3O";
    kname[19] = "O2";
    kname[20] = "HO2";
    kname[21] = "H2O2";
    kname[22] = "CO2";
    kname[23] = "CH3HCO";
    kname[24] = "CH3OCH2";
    kname[25] = "HCOOH";
    kname[26] = "CH3OCH3";
    kname[27] = "HOCH2O";
    kname[28] = "CH3OCO";
    kname[29] = "CH3OCHO";
    kname[30] = "CH3OCH2O";
    kname[31] = "CH3OCH2OH";
    kname[32] = "OCH2OCHO";
    kname[33] = "HOCH2OCO";
    kname[34] = "CH3OCH2O2";
    kname[35] = "CH2OCH2O2H";
    kname[36] = "HO2CH2OCHO";
    kname[37] = "O2CH2OCH2O2H";
    kname[38] = "N2";
}

/*compute the sparsity pattern of the chemistry Jacobian */
void SPARSITY_INFO( int * nJdata, int * consP, int NCELLS)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(1600);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(39);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[1600];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            for (int k=0; k<39; k++) {
                c_d[k] = 1.0/ 39.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<40; k++) {
        for (int l=0; l<40; l++) {
            if(J_h[ 40 * k + l] != 0.0){
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(1600);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(39);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[1600];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            for (int k=0; k<39; k++) {
                c_d[k] = 1.0/ 39.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<40; k++) {
        for (int l=0; l<40; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J_h[ 40 * k + l] != 0.0){
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(1600);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(39);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[1600];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            for (int k=0; k<39; k++) {
                c_d[k] = 1.0/ 39.000000 ;
            }
            aJacobian_precond(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<40; k++) {
        for (int l=0; l<40; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J_h[ 40 * k + l] != 0.0){
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

    amrex::Gpu::DeviceVector<amrex::Real> J_v(1600);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(39);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[1600];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            for (int k=0; k<39; k++) {
                c_d[k] = 1.0/ 39.000000 ;
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
        offset_row = nc * 40;
        offset_col = nc * 40;
        for (int k=0; k<40; k++) {
            for (int l=0; l<40; l++) {
                if(J_h[40*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(1600);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(39);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[1600];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            for (int k=0; k<39; k++) {
                c_d[k] = 1.0/ 39.000000 ;
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
            offset = nc * 40;
            for (int l=0; l<40; l++) {
                for (int k=0; k<40; k++) {
                    if(J_h[40*k + l] != 0.0) {
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
            offset = nc * 40;
            for (int l=0; l<40; l++) {
                for (int k=0; k<40; k++) {
                    if(J_h[40*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(1600);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(39);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[1600];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            for (int k=0; k<39; k++) {
                c_d[k] = 1.0/ 39.000000 ;
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
            offset = nc * 40;
            for (int l=0; l<40; l++) {
                for (int k=0; k<40; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp-1] = l+1 + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(J_h[40*k + l] != 0.0) {
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
            offset = nc * 40;
            for (int l=0; l<40; l++) {
                for (int k=0; k<40; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp] = l + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(J_h[40*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(1600);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(39);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[1600];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            for (int k=0; k<39; k++) {
                c_d[k] = 1.0/ 39.000000 ;
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
    for (int k=0; k<40; k++) {
        for (int l=0; l<40; l++) {
            if (k == l) {
                rowVals[nJdata_tmp] = l; 
                indx[nJdata_tmp] = 40*k + l;
                nJdata_tmp = nJdata_tmp + 1; 
            } else {
                if(J_h[40*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l; 
                    indx[nJdata_tmp] = 40*k + l;
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(1600);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(39);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[1600];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            for (int k=0; k<39; k++) {
                c_d[k] = 1.0/ 39.000000 ;
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
        for (int l=0; l<40; l++) {
            for (int k=0; k<40; k++) {
                if (k == l) {
                    colVals[nJdata_tmp-1] = l+1; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J_h[40*k + l] != 0.0) {
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
        for (int l=0; l<40; l++) {
            for (int k=0; k<40; k++) {
                if (k == l) {
                    colVals[nJdata_tmp] = l; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J_h[40*k + l] != 0.0) {
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

