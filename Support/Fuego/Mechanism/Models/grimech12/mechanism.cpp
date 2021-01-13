#ifndef MECHANISM_CPP
#define MECHANISM_CPP

#include "chemistry_file.H"

/* PURE CPU stuff  */
#ifndef AMREX_USE_GPU
namespace thermo
{
    amrex::Real fwd_A[177], fwd_beta[177], fwd_Ea[177];
    amrex::Real low_A[177], low_beta[177], low_Ea[177];
    amrex::Real rev_A[177], rev_beta[177], rev_Ea[177];
    amrex::Real troe_a[177],troe_Ts[177], troe_Tss[177], troe_Tsss[177];
    amrex::Real sri_a[177], sri_b[177], sri_c[177], sri_d[177], sri_e[177];
    amrex::Real activation_units[177], prefactor_units[177], phase_units[177];
    int is_PD[177], troe_len[177], sri_len[177], nTB[177], *TBid[177];
    amrex::Real *TB[177];
    std::vector<std::vector<amrex::Real>> kiv(177); 
    std::vector<std::vector<amrex::Real>> nuv(177); 
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
        if (*i > 177) {
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
    amrex::Real imw[32];

    get_imw(imw);

    for (int i=0; i<(*np); i++) {
        YOW[i] = 0.0;
    }

    for (int n=0; n<32; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] = y[n*(*np)+i] * imw[n];
            YOW[i] += x[n*(*np)+i];
        }
    }

    for (int i=0; i<(*np); i++) {
        YOW[i] = 1.0/YOW[i];
    }

    for (int n=0; n<32; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] *=  YOW[i];
        }
    }
}

/*Returns enthalpy in mass units (Eq 27.) */
void VCKHMS(int *  np, amrex::Real *  T,  amrex::Real *  hms)
{
    amrex::Real tc[5], h[32];
    amrex::Real imw[32];

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
    }

    for (int n=0; n<32; n++) {
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
    amrex::Real c[32*(*np)]; /*temporary storage */
    amrex::Real imw[32];

    get_imw(imw);

    /*See Eq 8 with an extra 1e6 so c goes to SI */
    for (int n=0; n<32; n++) {
        for (int i=0; i<(*np); i++) {
            c[n*(*np)+i] = 1.0e6 * rho[i] * y[n*(*np)+i] * imw[n];
        }
    }

    /*call productionRate */
    vproductionRate(*np, wdot, c, T);

    /*convert to chemkin units */
    for (int i=0; i<32*(*np); i++) {
        wdot[i] *= 1.0e-6;
    }
}


/*Compute P = rhoRT/W(y) */
void VCKPY(int *  np, amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  P)
{
    amrex::Real YOW[*np];
    amrex::Real imw[32];

    get_imw(imw);

    for (int i=0; i<(*np); i++) {
        YOW[i] = 0.0;
    }

    for (int n=0; n<32; n++) {
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
    amrex::Real k_f_s[177*npt], Kc_s[177*npt], mixture[npt], g_RT[32*npt];
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

    for (int n=0; n<32; n++) {
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
    vcomp_wdot_151_177(npt, wdot, mixture, sc, k_f_s, Kc_s, tc, invT, T);
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
        k_f_s[175*npt+i] = prefactor_units[175] * fwd_A[175] * exp(fwd_beta[175] * tc[i] - activation_units[175] * fwd_Ea[175] * invT[i]);
        k_f_s[176*npt+i] = prefactor_units[176] * fwd_A[176] * exp(fwd_beta[176] * tc[i] - activation_units[176] * fwd_Ea[176] * invT[i]);
    }
}

void vcomp_gibbs(int npt, amrex::Real *  g_RT, amrex::Real *  tc)
{
    /*compute the Gibbs free energy */
    for (int i=0; i<npt; i++) {
        amrex::Real tg[5], g[32];
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
    }
}

void vcomp_Kc(int npt, amrex::Real *  Kc_s, amrex::Real *  g_RT, amrex::Real *  invT)
{
    for (int i=0; i<npt; i++) {
        /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
        amrex::Real refC = (101325. / 8.31451) * invT[i];
        amrex::Real refCinv = 1.0 / refC;

        Kc_s[0*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[10*npt+i]) - (g_RT[12*npt+i]));
        Kc_s[1*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[12*npt+i]) - (g_RT[13*npt+i]));
        Kc_s[2*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[16*npt+i]) - (g_RT[17*npt+i]));
        Kc_s[3*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[17*npt+i]) - (g_RT[18*npt+i]));
        Kc_s[4*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[17*npt+i]) - (g_RT[19*npt+i]));
        Kc_s[5*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[18*npt+i]) - (g_RT[20*npt+i]));
        Kc_s[6*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[19*npt+i]) - (g_RT[20*npt+i]));
        Kc_s[7*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[21*npt+i]) - (g_RT[22*npt+i]));
        Kc_s[8*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[22*npt+i]) - (g_RT[23*npt+i]));
        Kc_s[9*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[23*npt+i]) - (g_RT[24*npt+i]));
        Kc_s[10*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[24*npt+i]) - (g_RT[25*npt+i]));
        Kc_s[11*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[25*npt+i]) - (g_RT[26*npt+i]));
        Kc_s[12*npt+i] = refCinv * exp((g_RT[0*npt+i] + g_RT[14*npt+i]) - (g_RT[17*npt+i]));
        Kc_s[13*npt+i] = refCinv * exp((2.000000 * g_RT[4*npt+i]) - (g_RT[7*npt+i]));
        Kc_s[14*npt+i] = refCinv * exp((g_RT[4*npt+i] + g_RT[12*npt+i]) - (g_RT[20*npt+i]));
        Kc_s[15*npt+i] = refCinv * exp((g_RT[9*npt+i] + g_RT[14*npt+i]) - (g_RT[27*npt+i]));
        Kc_s[16*npt+i] = refCinv * exp((g_RT[10*npt+i] + g_RT[14*npt+i]) - (g_RT[28*npt+i]));
        Kc_s[17*npt+i] = refCinv * exp((g_RT[5*npt+i] + g_RT[11*npt+i]) - (g_RT[20*npt+i]));
        Kc_s[18*npt+i] = refCinv * exp((2.000000 * g_RT[12*npt+i]) - (g_RT[26*npt+i]));
        Kc_s[19*npt+i] = refC * exp((g_RT[24*npt+i]) - (g_RT[0*npt+i] + g_RT[22*npt+i]));
        Kc_s[20*npt+i] = refCinv * exp((2.000000 * g_RT[2*npt+i]) - (g_RT[3*npt+i]));
        Kc_s[21*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[2*npt+i]) - (g_RT[4*npt+i]));
        Kc_s[22*npt+i] = refCinv * exp((g_RT[2*npt+i] + g_RT[14*npt+i]) - (g_RT[15*npt+i]));
        Kc_s[23*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[3*npt+i]) - (g_RT[6*npt+i]));
        Kc_s[24*npt+i] = refCinv * exp((2.000000 * g_RT[1*npt+i]) - (g_RT[0*npt+i]));
        Kc_s[25*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[4*npt+i]) - (g_RT[5*npt+i]));
        Kc_s[26*npt+i] = refC * exp((g_RT[16*npt+i]) - (g_RT[1*npt+i] + g_RT[14*npt+i]));
        Kc_s[27*npt+i] = exp((g_RT[0*npt+i] + g_RT[2*npt+i]) - (g_RT[1*npt+i] + g_RT[4*npt+i]));
        Kc_s[28*npt+i] = exp((g_RT[2*npt+i] + g_RT[6*npt+i]) - (g_RT[3*npt+i] + g_RT[4*npt+i]));
        Kc_s[29*npt+i] = exp((g_RT[2*npt+i] + g_RT[7*npt+i]) - (g_RT[4*npt+i] + g_RT[6*npt+i]));
        Kc_s[30*npt+i] = exp((g_RT[2*npt+i] + g_RT[9*npt+i]) - (g_RT[1*npt+i] + g_RT[14*npt+i]));
        Kc_s[31*npt+i] = exp((g_RT[2*npt+i] + g_RT[10*npt+i]) - (g_RT[1*npt+i] + g_RT[16*npt+i]));
        Kc_s[32*npt+i] = exp((g_RT[2*npt+i] + g_RT[11*npt+i]) - (g_RT[0*npt+i] + g_RT[14*npt+i]));
        Kc_s[33*npt+i] = exp((g_RT[2*npt+i] + g_RT[11*npt+i]) - (g_RT[1*npt+i] + g_RT[16*npt+i]));
        Kc_s[34*npt+i] = exp((g_RT[2*npt+i] + g_RT[12*npt+i]) - (g_RT[1*npt+i] + g_RT[17*npt+i]));
        Kc_s[35*npt+i] = exp((g_RT[2*npt+i] + g_RT[13*npt+i]) - (g_RT[4*npt+i] + g_RT[12*npt+i]));
        Kc_s[36*npt+i] = exp((g_RT[2*npt+i] + g_RT[16*npt+i]) - (g_RT[4*npt+i] + g_RT[14*npt+i]));
        Kc_s[37*npt+i] = exp((g_RT[2*npt+i] + g_RT[16*npt+i]) - (g_RT[1*npt+i] + g_RT[15*npt+i]));
        Kc_s[38*npt+i] = exp((g_RT[2*npt+i] + g_RT[17*npt+i]) - (g_RT[4*npt+i] + g_RT[16*npt+i]));
        Kc_s[39*npt+i] = exp((g_RT[2*npt+i] + g_RT[18*npt+i]) - (g_RT[4*npt+i] + g_RT[17*npt+i]));
        Kc_s[40*npt+i] = exp((g_RT[2*npt+i] + g_RT[19*npt+i]) - (g_RT[4*npt+i] + g_RT[17*npt+i]));
        Kc_s[41*npt+i] = exp((g_RT[2*npt+i] + g_RT[20*npt+i]) - (g_RT[4*npt+i] + g_RT[18*npt+i]));
        Kc_s[42*npt+i] = exp((g_RT[2*npt+i] + g_RT[20*npt+i]) - (g_RT[4*npt+i] + g_RT[19*npt+i]));
        Kc_s[43*npt+i] = exp((g_RT[2*npt+i] + g_RT[21*npt+i]) - (g_RT[9*npt+i] + g_RT[14*npt+i]));
        Kc_s[44*npt+i] = exp((g_RT[2*npt+i] + g_RT[22*npt+i]) - (g_RT[1*npt+i] + g_RT[27*npt+i]));
        Kc_s[45*npt+i] = exp((g_RT[2*npt+i] + g_RT[22*npt+i]) - (g_RT[4*npt+i] + g_RT[21*npt+i]));
        Kc_s[46*npt+i] = exp((g_RT[2*npt+i] + g_RT[22*npt+i]) - (g_RT[10*npt+i] + g_RT[14*npt+i]));
        Kc_s[47*npt+i] = exp((g_RT[2*npt+i] + g_RT[23*npt+i]) - (g_RT[1*npt+i] + g_RT[28*npt+i]));
        Kc_s[48*npt+i] = exp((g_RT[2*npt+i] + g_RT[24*npt+i]) - (g_RT[12*npt+i] + g_RT[16*npt+i]));
        Kc_s[49*npt+i] = exp((g_RT[2*npt+i] + g_RT[25*npt+i]) - (g_RT[12*npt+i] + g_RT[17*npt+i]));
        Kc_s[50*npt+i] = exp((g_RT[2*npt+i] + g_RT[26*npt+i]) - (g_RT[4*npt+i] + g_RT[25*npt+i]));
        Kc_s[51*npt+i] = refC * exp((g_RT[2*npt+i] + g_RT[27*npt+i]) - (g_RT[1*npt+i] + 2.000000 * g_RT[14*npt+i]));
        Kc_s[52*npt+i] = exp((g_RT[2*npt+i] + g_RT[28*npt+i]) - (g_RT[4*npt+i] + g_RT[27*npt+i]));
        Kc_s[53*npt+i] = exp((g_RT[2*npt+i] + g_RT[28*npt+i]) - (g_RT[10*npt+i] + g_RT[15*npt+i]));
        Kc_s[54*npt+i] = exp((g_RT[3*npt+i] + g_RT[14*npt+i]) - (g_RT[2*npt+i] + g_RT[15*npt+i]));
        Kc_s[55*npt+i] = exp((g_RT[3*npt+i] + g_RT[17*npt+i]) - (g_RT[6*npt+i] + g_RT[16*npt+i]));
        Kc_s[56*npt+i] = refCinv * exp((g_RT[1*npt+i] + 2.000000 * g_RT[3*npt+i]) - (g_RT[3*npt+i] + g_RT[6*npt+i]));
        Kc_s[57*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[3*npt+i] + g_RT[5*npt+i]) - (g_RT[5*npt+i] + g_RT[6*npt+i]));
        Kc_s[58*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[3*npt+i] + g_RT[30*npt+i]) - (g_RT[6*npt+i] + g_RT[30*npt+i]));
        Kc_s[59*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[3*npt+i] + g_RT[31*npt+i]) - (g_RT[6*npt+i] + g_RT[31*npt+i]));
        Kc_s[60*npt+i] = exp((g_RT[1*npt+i] + g_RT[3*npt+i]) - (g_RT[2*npt+i] + g_RT[4*npt+i]));
        Kc_s[61*npt+i] = refCinv * exp((g_RT[0*npt+i] + 2.000000 * g_RT[1*npt+i]) - (2.000000 * g_RT[0*npt+i]));
        Kc_s[62*npt+i] = refCinv * exp((2.000000 * g_RT[1*npt+i] + g_RT[5*npt+i]) - (g_RT[0*npt+i] + g_RT[5*npt+i]));
        Kc_s[63*npt+i] = refCinv * exp((2.000000 * g_RT[1*npt+i] + g_RT[15*npt+i]) - (g_RT[0*npt+i] + g_RT[15*npt+i]));
        Kc_s[64*npt+i] = exp((g_RT[1*npt+i] + g_RT[6*npt+i]) - (g_RT[2*npt+i] + g_RT[5*npt+i]));
        Kc_s[65*npt+i] = exp((g_RT[1*npt+i] + g_RT[6*npt+i]) - (g_RT[0*npt+i] + g_RT[3*npt+i]));
        Kc_s[66*npt+i] = exp((g_RT[1*npt+i] + g_RT[6*npt+i]) - (2.000000 * g_RT[4*npt+i]));
        Kc_s[67*npt+i] = exp((g_RT[1*npt+i] + g_RT[7*npt+i]) - (g_RT[0*npt+i] + g_RT[6*npt+i]));
        Kc_s[68*npt+i] = exp((g_RT[1*npt+i] + g_RT[7*npt+i]) - (g_RT[4*npt+i] + g_RT[5*npt+i]));
        Kc_s[69*npt+i] = exp((g_RT[1*npt+i] + g_RT[9*npt+i]) - (g_RT[0*npt+i] + g_RT[8*npt+i]));
        Kc_s[70*npt+i] = exp((g_RT[1*npt+i] + g_RT[11*npt+i]) - (g_RT[0*npt+i] + g_RT[9*npt+i]));
        Kc_s[71*npt+i] = exp((g_RT[1*npt+i] + g_RT[13*npt+i]) - (g_RT[0*npt+i] + g_RT[12*npt+i]));
        Kc_s[72*npt+i] = exp((g_RT[1*npt+i] + g_RT[16*npt+i]) - (g_RT[0*npt+i] + g_RT[14*npt+i]));
        Kc_s[73*npt+i] = exp((g_RT[1*npt+i] + g_RT[17*npt+i]) - (g_RT[0*npt+i] + g_RT[16*npt+i]));
        Kc_s[74*npt+i] = exp((g_RT[1*npt+i] + g_RT[18*npt+i]) - (g_RT[0*npt+i] + g_RT[17*npt+i]));
        Kc_s[75*npt+i] = exp((g_RT[1*npt+i] + g_RT[18*npt+i]) - (g_RT[4*npt+i] + g_RT[12*npt+i]));
        Kc_s[76*npt+i] = exp((g_RT[1*npt+i] + g_RT[18*npt+i]) - (g_RT[5*npt+i] + g_RT[11*npt+i]));
        Kc_s[77*npt+i] = exp((g_RT[1*npt+i] + g_RT[19*npt+i]) - (g_RT[1*npt+i] + g_RT[18*npt+i]));
        Kc_s[78*npt+i] = exp((g_RT[1*npt+i] + g_RT[19*npt+i]) - (g_RT[0*npt+i] + g_RT[17*npt+i]));
        Kc_s[79*npt+i] = exp((g_RT[1*npt+i] + g_RT[19*npt+i]) - (g_RT[4*npt+i] + g_RT[12*npt+i]));
        Kc_s[80*npt+i] = exp((g_RT[1*npt+i] + g_RT[19*npt+i]) - (g_RT[5*npt+i] + g_RT[11*npt+i]));
        Kc_s[81*npt+i] = exp((g_RT[1*npt+i] + g_RT[20*npt+i]) - (g_RT[0*npt+i] + g_RT[18*npt+i]));
        Kc_s[82*npt+i] = exp((g_RT[1*npt+i] + g_RT[20*npt+i]) - (g_RT[0*npt+i] + g_RT[19*npt+i]));
        Kc_s[83*npt+i] = exp((g_RT[1*npt+i] + g_RT[23*npt+i]) - (g_RT[0*npt+i] + g_RT[22*npt+i]));
        Kc_s[84*npt+i] = exp((g_RT[1*npt+i] + g_RT[24*npt+i]) - (g_RT[0*npt+i] + g_RT[23*npt+i]));
        Kc_s[85*npt+i] = exp((g_RT[1*npt+i] + g_RT[25*npt+i]) - (g_RT[0*npt+i] + g_RT[24*npt+i]));
        Kc_s[86*npt+i] = exp((g_RT[1*npt+i] + g_RT[26*npt+i]) - (g_RT[0*npt+i] + g_RT[25*npt+i]));
        Kc_s[87*npt+i] = exp((g_RT[1*npt+i] + g_RT[27*npt+i]) - (g_RT[11*npt+i] + g_RT[14*npt+i]));
        Kc_s[88*npt+i] = exp((g_RT[1*npt+i] + g_RT[28*npt+i]) - (g_RT[0*npt+i] + g_RT[27*npt+i]));
        Kc_s[89*npt+i] = exp((g_RT[1*npt+i] + g_RT[28*npt+i]) - (g_RT[12*npt+i] + g_RT[14*npt+i]));
        Kc_s[90*npt+i] = exp((g_RT[1*npt+i] + g_RT[29*npt+i]) - (g_RT[1*npt+i] + g_RT[28*npt+i]));
        Kc_s[91*npt+i] = exp((g_RT[0*npt+i] + g_RT[4*npt+i]) - (g_RT[1*npt+i] + g_RT[5*npt+i]));
        Kc_s[92*npt+i] = exp((2.000000 * g_RT[4*npt+i]) - (g_RT[2*npt+i] + g_RT[5*npt+i]));
        Kc_s[93*npt+i] = exp((g_RT[4*npt+i] + g_RT[6*npt+i]) - (g_RT[3*npt+i] + g_RT[5*npt+i]));
        Kc_s[94*npt+i] = exp((g_RT[4*npt+i] + g_RT[7*npt+i]) - (g_RT[5*npt+i] + g_RT[6*npt+i]));
        Kc_s[95*npt+i] = exp((g_RT[4*npt+i] + g_RT[7*npt+i]) - (g_RT[5*npt+i] + g_RT[6*npt+i]));
        Kc_s[96*npt+i] = exp((g_RT[4*npt+i] + g_RT[8*npt+i]) - (g_RT[1*npt+i] + g_RT[14*npt+i]));
        Kc_s[97*npt+i] = exp((g_RT[4*npt+i] + g_RT[9*npt+i]) - (g_RT[1*npt+i] + g_RT[16*npt+i]));
        Kc_s[98*npt+i] = exp((g_RT[4*npt+i] + g_RT[10*npt+i]) - (g_RT[1*npt+i] + g_RT[17*npt+i]));
        Kc_s[99*npt+i] = exp((g_RT[4*npt+i] + g_RT[10*npt+i]) - (g_RT[5*npt+i] + g_RT[9*npt+i]));
        Kc_s[100*npt+i] = exp((g_RT[4*npt+i] + g_RT[11*npt+i]) - (g_RT[1*npt+i] + g_RT[17*npt+i]));
        Kc_s[101*npt+i] = exp((g_RT[4*npt+i] + g_RT[12*npt+i]) - (g_RT[5*npt+i] + g_RT[10*npt+i]));
        Kc_s[102*npt+i] = exp((g_RT[4*npt+i] + g_RT[12*npt+i]) - (g_RT[5*npt+i] + g_RT[11*npt+i]));
        Kc_s[103*npt+i] = exp((g_RT[4*npt+i] + g_RT[13*npt+i]) - (g_RT[5*npt+i] + g_RT[12*npt+i]));
        Kc_s[104*npt+i] = exp((g_RT[4*npt+i] + g_RT[14*npt+i]) - (g_RT[1*npt+i] + g_RT[15*npt+i]));
        Kc_s[105*npt+i] = exp((g_RT[4*npt+i] + g_RT[16*npt+i]) - (g_RT[5*npt+i] + g_RT[14*npt+i]));
        Kc_s[106*npt+i] = exp((g_RT[4*npt+i] + g_RT[17*npt+i]) - (g_RT[5*npt+i] + g_RT[16*npt+i]));
        Kc_s[107*npt+i] = exp((g_RT[4*npt+i] + g_RT[18*npt+i]) - (g_RT[5*npt+i] + g_RT[17*npt+i]));
        Kc_s[108*npt+i] = exp((g_RT[4*npt+i] + g_RT[19*npt+i]) - (g_RT[5*npt+i] + g_RT[17*npt+i]));
        Kc_s[109*npt+i] = exp((g_RT[4*npt+i] + g_RT[20*npt+i]) - (g_RT[5*npt+i] + g_RT[18*npt+i]));
        Kc_s[110*npt+i] = exp((g_RT[4*npt+i] + g_RT[20*npt+i]) - (g_RT[5*npt+i] + g_RT[19*npt+i]));
        Kc_s[111*npt+i] = exp((g_RT[4*npt+i] + g_RT[21*npt+i]) - (g_RT[1*npt+i] + g_RT[27*npt+i]));
        Kc_s[112*npt+i] = exp((g_RT[4*npt+i] + g_RT[22*npt+i]) - (g_RT[1*npt+i] + g_RT[28*npt+i]));
        Kc_s[113*npt+i] = exp((g_RT[4*npt+i] + g_RT[22*npt+i]) - (g_RT[1*npt+i] + g_RT[29*npt+i]));
        Kc_s[114*npt+i] = exp((g_RT[4*npt+i] + g_RT[22*npt+i]) - (g_RT[5*npt+i] + g_RT[21*npt+i]));
        Kc_s[115*npt+i] = exp((g_RT[4*npt+i] + g_RT[22*npt+i]) - (g_RT[12*npt+i] + g_RT[14*npt+i]));
        Kc_s[116*npt+i] = exp((g_RT[4*npt+i] + g_RT[23*npt+i]) - (g_RT[5*npt+i] + g_RT[22*npt+i]));
        Kc_s[117*npt+i] = exp((g_RT[4*npt+i] + g_RT[24*npt+i]) - (g_RT[5*npt+i] + g_RT[23*npt+i]));
        Kc_s[118*npt+i] = exp((g_RT[4*npt+i] + g_RT[26*npt+i]) - (g_RT[5*npt+i] + g_RT[25*npt+i]));
        Kc_s[119*npt+i] = exp((g_RT[4*npt+i] + g_RT[28*npt+i]) - (g_RT[5*npt+i] + g_RT[27*npt+i]));
        Kc_s[120*npt+i] = exp((2.000000 * g_RT[6*npt+i]) - (g_RT[3*npt+i] + g_RT[7*npt+i]));
        Kc_s[121*npt+i] = exp((2.000000 * g_RT[6*npt+i]) - (g_RT[3*npt+i] + g_RT[7*npt+i]));
        Kc_s[122*npt+i] = exp((g_RT[6*npt+i] + g_RT[10*npt+i]) - (g_RT[4*npt+i] + g_RT[17*npt+i]));
        Kc_s[123*npt+i] = exp((g_RT[6*npt+i] + g_RT[12*npt+i]) - (g_RT[3*npt+i] + g_RT[13*npt+i]));
        Kc_s[124*npt+i] = exp((g_RT[6*npt+i] + g_RT[12*npt+i]) - (g_RT[4*npt+i] + g_RT[19*npt+i]));
        Kc_s[125*npt+i] = exp((g_RT[6*npt+i] + g_RT[14*npt+i]) - (g_RT[4*npt+i] + g_RT[15*npt+i]));
        Kc_s[126*npt+i] = exp((g_RT[6*npt+i] + g_RT[17*npt+i]) - (g_RT[7*npt+i] + g_RT[16*npt+i]));
        Kc_s[127*npt+i] = exp((g_RT[3*npt+i] + g_RT[8*npt+i]) - (g_RT[2*npt+i] + g_RT[14*npt+i]));
        Kc_s[128*npt+i] = exp((g_RT[8*npt+i] + g_RT[10*npt+i]) - (g_RT[1*npt+i] + g_RT[21*npt+i]));
        Kc_s[129*npt+i] = exp((g_RT[8*npt+i] + g_RT[12*npt+i]) - (g_RT[1*npt+i] + g_RT[22*npt+i]));
        Kc_s[130*npt+i] = exp((g_RT[3*npt+i] + g_RT[9*npt+i]) - (g_RT[2*npt+i] + g_RT[16*npt+i]));
        Kc_s[131*npt+i] = exp((g_RT[0*npt+i] + g_RT[9*npt+i]) - (g_RT[1*npt+i] + g_RT[10*npt+i]));
        Kc_s[132*npt+i] = exp((g_RT[5*npt+i] + g_RT[9*npt+i]) - (g_RT[1*npt+i] + g_RT[17*npt+i]));
        Kc_s[133*npt+i] = exp((g_RT[9*npt+i] + g_RT[10*npt+i]) - (g_RT[1*npt+i] + g_RT[22*npt+i]));
        Kc_s[134*npt+i] = exp((g_RT[9*npt+i] + g_RT[12*npt+i]) - (g_RT[1*npt+i] + g_RT[23*npt+i]));
        Kc_s[135*npt+i] = exp((g_RT[9*npt+i] + g_RT[13*npt+i]) - (g_RT[1*npt+i] + g_RT[24*npt+i]));
        Kc_s[136*npt+i] = exp((g_RT[9*npt+i] + g_RT[15*npt+i]) - (g_RT[14*npt+i] + g_RT[16*npt+i]));
        Kc_s[137*npt+i] = exp((g_RT[9*npt+i] + g_RT[17*npt+i]) - (g_RT[1*npt+i] + g_RT[28*npt+i]));
        Kc_s[138*npt+i] = exp((g_RT[9*npt+i] + g_RT[27*npt+i]) - (g_RT[14*npt+i] + g_RT[22*npt+i]));
        Kc_s[139*npt+i] = exp((g_RT[3*npt+i] + g_RT[10*npt+i]) - (g_RT[4*npt+i] + g_RT[16*npt+i]));
        Kc_s[140*npt+i] = exp((g_RT[0*npt+i] + g_RT[10*npt+i]) - (g_RT[1*npt+i] + g_RT[12*npt+i]));
        Kc_s[141*npt+i] = exp((2.000000 * g_RT[10*npt+i]) - (g_RT[0*npt+i] + g_RT[22*npt+i]));
        Kc_s[142*npt+i] = exp((g_RT[10*npt+i] + g_RT[12*npt+i]) - (g_RT[1*npt+i] + g_RT[24*npt+i]));
        Kc_s[143*npt+i] = exp((g_RT[10*npt+i] + g_RT[13*npt+i]) - (2.000000 * g_RT[12*npt+i]));
        Kc_s[144*npt+i] = exp((g_RT[10*npt+i] + g_RT[27*npt+i]) - (g_RT[14*npt+i] + g_RT[23*npt+i]));
        Kc_s[145*npt+i] = exp((g_RT[11*npt+i] + g_RT[30*npt+i]) - (g_RT[10*npt+i] + g_RT[30*npt+i]));
        Kc_s[146*npt+i] = exp((g_RT[11*npt+i] + g_RT[31*npt+i]) - (g_RT[10*npt+i] + g_RT[31*npt+i]));
        Kc_s[147*npt+i] = refC * exp((g_RT[3*npt+i] + g_RT[11*npt+i]) - (g_RT[1*npt+i] + g_RT[4*npt+i] + g_RT[14*npt+i]));
        Kc_s[148*npt+i] = exp((g_RT[3*npt+i] + g_RT[11*npt+i]) - (g_RT[5*npt+i] + g_RT[14*npt+i]));
        Kc_s[149*npt+i] = exp((g_RT[0*npt+i] + g_RT[11*npt+i]) - (g_RT[1*npt+i] + g_RT[12*npt+i]));
        Kc_s[150*npt+i] = exp((g_RT[5*npt+i] + g_RT[11*npt+i]) - (g_RT[5*npt+i] + g_RT[10*npt+i]));
        Kc_s[151*npt+i] = exp((g_RT[11*npt+i] + g_RT[12*npt+i]) - (g_RT[1*npt+i] + g_RT[24*npt+i]));
        Kc_s[152*npt+i] = exp((g_RT[11*npt+i] + g_RT[13*npt+i]) - (2.000000 * g_RT[12*npt+i]));
        Kc_s[153*npt+i] = exp((g_RT[11*npt+i] + g_RT[14*npt+i]) - (g_RT[10*npt+i] + g_RT[14*npt+i]));
        Kc_s[154*npt+i] = exp((g_RT[11*npt+i] + g_RT[15*npt+i]) - (g_RT[10*npt+i] + g_RT[15*npt+i]));
        Kc_s[155*npt+i] = exp((g_RT[11*npt+i] + g_RT[15*npt+i]) - (g_RT[14*npt+i] + g_RT[17*npt+i]));
        Kc_s[156*npt+i] = exp((g_RT[11*npt+i] + g_RT[26*npt+i]) - (g_RT[12*npt+i] + g_RT[25*npt+i]));
        Kc_s[157*npt+i] = exp((g_RT[3*npt+i] + g_RT[12*npt+i]) - (g_RT[2*npt+i] + g_RT[19*npt+i]));
        Kc_s[158*npt+i] = exp((g_RT[3*npt+i] + g_RT[12*npt+i]) - (g_RT[4*npt+i] + g_RT[17*npt+i]));
        Kc_s[159*npt+i] = exp((g_RT[7*npt+i] + g_RT[12*npt+i]) - (g_RT[6*npt+i] + g_RT[13*npt+i]));
        Kc_s[160*npt+i] = exp((2.000000 * g_RT[12*npt+i]) - (g_RT[1*npt+i] + g_RT[25*npt+i]));
        Kc_s[161*npt+i] = exp((g_RT[12*npt+i] + g_RT[16*npt+i]) - (g_RT[13*npt+i] + g_RT[14*npt+i]));
        Kc_s[162*npt+i] = exp((g_RT[12*npt+i] + g_RT[17*npt+i]) - (g_RT[13*npt+i] + g_RT[16*npt+i]));
        Kc_s[163*npt+i] = exp((g_RT[12*npt+i] + g_RT[20*npt+i]) - (g_RT[13*npt+i] + g_RT[18*npt+i]));
        Kc_s[164*npt+i] = exp((g_RT[12*npt+i] + g_RT[20*npt+i]) - (g_RT[13*npt+i] + g_RT[19*npt+i]));
        Kc_s[165*npt+i] = exp((g_RT[12*npt+i] + g_RT[24*npt+i]) - (g_RT[13*npt+i] + g_RT[23*npt+i]));
        Kc_s[166*npt+i] = exp((g_RT[12*npt+i] + g_RT[26*npt+i]) - (g_RT[13*npt+i] + g_RT[25*npt+i]));
        Kc_s[167*npt+i] = refC * exp((g_RT[5*npt+i] + g_RT[16*npt+i]) - (g_RT[1*npt+i] + g_RT[5*npt+i] + g_RT[14*npt+i]));
        Kc_s[168*npt+i] = exp((g_RT[3*npt+i] + g_RT[16*npt+i]) - (g_RT[6*npt+i] + g_RT[14*npt+i]));
        Kc_s[169*npt+i] = exp((g_RT[3*npt+i] + g_RT[18*npt+i]) - (g_RT[6*npt+i] + g_RT[17*npt+i]));
        Kc_s[170*npt+i] = exp((g_RT[3*npt+i] + g_RT[19*npt+i]) - (g_RT[6*npt+i] + g_RT[17*npt+i]));
        Kc_s[171*npt+i] = exp((g_RT[3*npt+i] + g_RT[21*npt+i]) - (g_RT[14*npt+i] + g_RT[16*npt+i]));
        Kc_s[172*npt+i] = exp((g_RT[0*npt+i] + g_RT[21*npt+i]) - (g_RT[1*npt+i] + g_RT[22*npt+i]));
        Kc_s[173*npt+i] = exp((g_RT[3*npt+i] + g_RT[23*npt+i]) - (g_RT[16*npt+i] + g_RT[17*npt+i]));
        Kc_s[174*npt+i] = exp((g_RT[3*npt+i] + g_RT[25*npt+i]) - (g_RT[6*npt+i] + g_RT[24*npt+i]));
        Kc_s[175*npt+i] = refC * exp((g_RT[3*npt+i] + g_RT[27*npt+i]) - (g_RT[4*npt+i] + 2.000000 * g_RT[14*npt+i]));
        Kc_s[176*npt+i] = refC * exp((2.000000 * g_RT[27*npt+i]) - (2.000000 * g_RT[14*npt+i] + g_RT[22*npt+i]));
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
        phi_f = sc[1*npt+i]*sc[10*npt+i];
        alpha = mixture[i] + (TB[0][0] - 1)*sc[0*npt+i] + (TB[0][1] - 1)*sc[5*npt+i] + (TB[0][2] - 1)*sc[13*npt+i] + (TB[0][3] - 1)*sc[14*npt+i] + (TB[0][4] - 1)*sc[15*npt+i] + (TB[0][5] - 1)*sc[26*npt+i] + (TB[0][6] - 1)*sc[31*npt+i];
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
        phi_r = sc[12*npt+i];
        Kc = Kc_s[0*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[10*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;

        /*reaction 2: H + CH3 (+M) <=> CH4 (+M) */
        phi_f = sc[1*npt+i]*sc[12*npt+i];
        alpha = mixture[i] + (TB[1][0] - 1)*sc[0*npt+i] + (TB[1][1] - 1)*sc[5*npt+i] + (TB[1][2] - 1)*sc[13*npt+i] + (TB[1][3] - 1)*sc[14*npt+i] + (TB[1][4] - 1)*sc[15*npt+i] + (TB[1][5] - 1)*sc[26*npt+i] + (TB[1][6] - 1)*sc[31*npt+i];
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
        wdot[12*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;

        /*reaction 3: H + HCO (+M) <=> CH2O (+M) */
        phi_f = sc[1*npt+i]*sc[16*npt+i];
        alpha = mixture[i] + (TB[2][0] - 1)*sc[0*npt+i] + (TB[2][1] - 1)*sc[5*npt+i] + (TB[2][2] - 1)*sc[13*npt+i] + (TB[2][3] - 1)*sc[14*npt+i] + (TB[2][4] - 1)*sc[15*npt+i] + (TB[2][5] - 1)*sc[26*npt+i] + (TB[2][6] - 1)*sc[31*npt+i];
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
        phi_r = sc[17*npt+i];
        Kc = Kc_s[2*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[16*npt+i] -= qdot;
        wdot[17*npt+i] += qdot;

        /*reaction 4: H + CH2O (+M) <=> CH2OH (+M) */
        phi_f = sc[1*npt+i]*sc[17*npt+i];
        alpha = mixture[i] + (TB[3][0] - 1)*sc[0*npt+i] + (TB[3][1] - 1)*sc[5*npt+i] + (TB[3][2] - 1)*sc[13*npt+i] + (TB[3][3] - 1)*sc[14*npt+i] + (TB[3][4] - 1)*sc[15*npt+i] + (TB[3][5] - 1)*sc[26*npt+i];
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
        phi_r = sc[18*npt+i];
        Kc = Kc_s[3*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[17*npt+i] -= qdot;
        wdot[18*npt+i] += qdot;

        /*reaction 5: H + CH2O (+M) <=> CH3O (+M) */
        phi_f = sc[1*npt+i]*sc[17*npt+i];
        alpha = mixture[i] + (TB[4][0] - 1)*sc[0*npt+i] + (TB[4][1] - 1)*sc[5*npt+i] + (TB[4][2] - 1)*sc[13*npt+i] + (TB[4][3] - 1)*sc[14*npt+i] + (TB[4][4] - 1)*sc[15*npt+i] + (TB[4][5] - 1)*sc[26*npt+i];
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
        phi_r = sc[19*npt+i];
        Kc = Kc_s[4*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[17*npt+i] -= qdot;
        wdot[19*npt+i] += qdot;

        /*reaction 6: H + CH2OH (+M) <=> CH3OH (+M) */
        phi_f = sc[1*npt+i]*sc[18*npt+i];
        alpha = mixture[i] + (TB[5][0] - 1)*sc[0*npt+i] + (TB[5][1] - 1)*sc[5*npt+i] + (TB[5][2] - 1)*sc[13*npt+i] + (TB[5][3] - 1)*sc[14*npt+i] + (TB[5][4] - 1)*sc[15*npt+i] + (TB[5][5] - 1)*sc[26*npt+i];
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
        phi_r = sc[20*npt+i];
        Kc = Kc_s[5*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[18*npt+i] -= qdot;
        wdot[20*npt+i] += qdot;

        /*reaction 7: H + CH3O (+M) <=> CH3OH (+M) */
        phi_f = sc[1*npt+i]*sc[19*npt+i];
        alpha = mixture[i] + (TB[6][0] - 1)*sc[0*npt+i] + (TB[6][1] - 1)*sc[5*npt+i] + (TB[6][2] - 1)*sc[13*npt+i] + (TB[6][3] - 1)*sc[14*npt+i] + (TB[6][4] - 1)*sc[15*npt+i] + (TB[6][5] - 1)*sc[26*npt+i];
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
        phi_r = sc[20*npt+i];
        Kc = Kc_s[6*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[19*npt+i] -= qdot;
        wdot[20*npt+i] += qdot;

        /*reaction 8: H + C2H (+M) <=> C2H2 (+M) */
        phi_f = sc[1*npt+i]*sc[21*npt+i];
        alpha = mixture[i] + (TB[7][0] - 1)*sc[0*npt+i] + (TB[7][1] - 1)*sc[5*npt+i] + (TB[7][2] - 1)*sc[13*npt+i] + (TB[7][3] - 1)*sc[14*npt+i] + (TB[7][4] - 1)*sc[15*npt+i] + (TB[7][5] - 1)*sc[26*npt+i] + (TB[7][6] - 1)*sc[31*npt+i];
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
        phi_r = sc[22*npt+i];
        Kc = Kc_s[7*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[21*npt+i] -= qdot;
        wdot[22*npt+i] += qdot;

        /*reaction 9: H + C2H2 (+M) <=> C2H3 (+M) */
        phi_f = sc[1*npt+i]*sc[22*npt+i];
        alpha = mixture[i] + (TB[8][0] - 1)*sc[0*npt+i] + (TB[8][1] - 1)*sc[5*npt+i] + (TB[8][2] - 1)*sc[13*npt+i] + (TB[8][3] - 1)*sc[14*npt+i] + (TB[8][4] - 1)*sc[15*npt+i] + (TB[8][5] - 1)*sc[26*npt+i] + (TB[8][6] - 1)*sc[31*npt+i];
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
        phi_r = sc[23*npt+i];
        Kc = Kc_s[8*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[22*npt+i] -= qdot;
        wdot[23*npt+i] += qdot;

        /*reaction 10: H + C2H3 (+M) <=> C2H4 (+M) */
        phi_f = sc[1*npt+i]*sc[23*npt+i];
        alpha = mixture[i] + (TB[9][0] - 1)*sc[0*npt+i] + (TB[9][1] - 1)*sc[5*npt+i] + (TB[9][2] - 1)*sc[13*npt+i] + (TB[9][3] - 1)*sc[14*npt+i] + (TB[9][4] - 1)*sc[15*npt+i] + (TB[9][5] - 1)*sc[26*npt+i] + (TB[9][6] - 1)*sc[31*npt+i];
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
        phi_r = sc[24*npt+i];
        Kc = Kc_s[9*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[23*npt+i] -= qdot;
        wdot[24*npt+i] += qdot;

        /*reaction 11: H + C2H4 (+M) <=> C2H5 (+M) */
        phi_f = sc[1*npt+i]*sc[24*npt+i];
        alpha = mixture[i] + (TB[10][0] - 1)*sc[0*npt+i] + (TB[10][1] - 1)*sc[5*npt+i] + (TB[10][2] - 1)*sc[13*npt+i] + (TB[10][3] - 1)*sc[14*npt+i] + (TB[10][4] - 1)*sc[15*npt+i] + (TB[10][5] - 1)*sc[26*npt+i] + (TB[10][6] - 1)*sc[31*npt+i];
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
        phi_r = sc[25*npt+i];
        Kc = Kc_s[10*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[24*npt+i] -= qdot;
        wdot[25*npt+i] += qdot;

        /*reaction 12: H + C2H5 (+M) <=> C2H6 (+M) */
        phi_f = sc[1*npt+i]*sc[25*npt+i];
        alpha = mixture[i] + (TB[11][0] - 1)*sc[0*npt+i] + (TB[11][1] - 1)*sc[5*npt+i] + (TB[11][2] - 1)*sc[13*npt+i] + (TB[11][3] - 1)*sc[14*npt+i] + (TB[11][4] - 1)*sc[15*npt+i] + (TB[11][5] - 1)*sc[26*npt+i] + (TB[11][6] - 1)*sc[31*npt+i];
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
        phi_r = sc[26*npt+i];
        Kc = Kc_s[11*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[25*npt+i] -= qdot;
        wdot[26*npt+i] += qdot;

        /*reaction 13: H2 + CO (+M) <=> CH2O (+M) */
        phi_f = sc[0*npt+i]*sc[14*npt+i];
        alpha = mixture[i] + (TB[12][0] - 1)*sc[0*npt+i] + (TB[12][1] - 1)*sc[5*npt+i] + (TB[12][2] - 1)*sc[13*npt+i] + (TB[12][3] - 1)*sc[14*npt+i] + (TB[12][4] - 1)*sc[15*npt+i] + (TB[12][5] - 1)*sc[26*npt+i] + (TB[12][6] - 1)*sc[31*npt+i];
        k_f = k_f_s[12*npt+i];
        redP = alpha / k_f * phase_units[12] * low_A[12] * exp(low_beta[12] * tc[i] - activation_units[12] * low_Ea[12] * invT[i]);
        F = redP / (1 + redP);
        logPred = log10(redP);
        logFcent = log10(
            (fabs(troe_Tsss[12]) > 1.e-100 ? (1.-troe_a[12])*exp(-T[i]/troe_Tsss[12]) : 0.) 
            + (fabs(troe_Ts[12]) > 1.e-100 ? troe_a[12] * exp(-T[i]/troe_Ts[12]) : 0.) 
            + (troe_len[12] == 4 ? exp(-troe_Tss[12] * invT[i]) : 0.) );
        troe_c = -.4 - .67 * logFcent;
        troe_n = .75 - 1.27 * logFcent;
        troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
        F_troe = pow(10., logFcent / (1.0 + troe*troe));
        F *= F_troe;
        k_f *= F;
        q_f = phi_f * k_f;
        phi_r = sc[17*npt+i];
        Kc = Kc_s[12*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[14*npt+i] -= qdot;
        wdot[17*npt+i] += qdot;

        /*reaction 14: 2.000000 OH (+M) <=> H2O2 (+M) */
        phi_f = pow(sc[4*npt+i], 2.000000);
        alpha = mixture[i] + (TB[13][0] - 1)*sc[0*npt+i] + (TB[13][1] - 1)*sc[5*npt+i] + (TB[13][2] - 1)*sc[13*npt+i] + (TB[13][3] - 1)*sc[14*npt+i] + (TB[13][4] - 1)*sc[15*npt+i] + (TB[13][5] - 1)*sc[26*npt+i] + (TB[13][6] - 1)*sc[31*npt+i];
        k_f = k_f_s[13*npt+i];
        redP = alpha / k_f * phase_units[13] * low_A[13] * exp(low_beta[13] * tc[i] - activation_units[13] * low_Ea[13] * invT[i]);
        F = redP / (1 + redP);
        logPred = log10(redP);
        logFcent = log10(
            (fabs(troe_Tsss[13]) > 1.e-100 ? (1.-troe_a[13])*exp(-T[i]/troe_Tsss[13]) : 0.) 
            + (fabs(troe_Ts[13]) > 1.e-100 ? troe_a[13] * exp(-T[i]/troe_Ts[13]) : 0.) 
            + (troe_len[13] == 4 ? exp(-troe_Tss[13] * invT[i]) : 0.) );
        troe_c = -.4 - .67 * logFcent;
        troe_n = .75 - 1.27 * logFcent;
        troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
        F_troe = pow(10., logFcent / (1.0 + troe*troe));
        F *= F_troe;
        k_f *= F;
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i];
        Kc = Kc_s[13*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= 2.000000 * qdot;
        wdot[7*npt+i] += qdot;

        /*reaction 15: OH + CH3 (+M) <=> CH3OH (+M) */
        phi_f = sc[4*npt+i]*sc[12*npt+i];
        alpha = mixture[i] + (TB[14][0] - 1)*sc[0*npt+i] + (TB[14][1] - 1)*sc[5*npt+i] + (TB[14][2] - 1)*sc[13*npt+i] + (TB[14][3] - 1)*sc[14*npt+i] + (TB[14][4] - 1)*sc[15*npt+i] + (TB[14][5] - 1)*sc[26*npt+i];
        k_f = k_f_s[14*npt+i];
        redP = alpha / k_f * phase_units[14] * low_A[14] * exp(low_beta[14] * tc[i] - activation_units[14] * low_Ea[14] * invT[i]);
        F = redP / (1 + redP);
        logPred = log10(redP);
        logFcent = log10(
            (fabs(troe_Tsss[14]) > 1.e-100 ? (1.-troe_a[14])*exp(-T[i]/troe_Tsss[14]) : 0.) 
            + (fabs(troe_Ts[14]) > 1.e-100 ? troe_a[14] * exp(-T[i]/troe_Ts[14]) : 0.) 
            + (troe_len[14] == 4 ? exp(-troe_Tss[14] * invT[i]) : 0.) );
        troe_c = -.4 - .67 * logFcent;
        troe_n = .75 - 1.27 * logFcent;
        troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
        F_troe = pow(10., logFcent / (1.0 + troe*troe));
        F *= F_troe;
        k_f *= F;
        q_f = phi_f * k_f;
        phi_r = sc[20*npt+i];
        Kc = Kc_s[14*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[12*npt+i] -= qdot;
        wdot[20*npt+i] += qdot;

        /*reaction 16: CH + CO (+M) <=> HCCO (+M) */
        phi_f = sc[9*npt+i]*sc[14*npt+i];
        alpha = mixture[i] + (TB[15][0] - 1)*sc[0*npt+i] + (TB[15][1] - 1)*sc[5*npt+i] + (TB[15][2] - 1)*sc[13*npt+i] + (TB[15][3] - 1)*sc[14*npt+i] + (TB[15][4] - 1)*sc[15*npt+i] + (TB[15][5] - 1)*sc[26*npt+i] + (TB[15][6] - 1)*sc[31*npt+i];
        k_f = k_f_s[15*npt+i];
        redP = alpha / k_f * phase_units[15] * low_A[15] * exp(low_beta[15] * tc[i] - activation_units[15] * low_Ea[15] * invT[i]);
        F = redP / (1 + redP);
        logPred = log10(redP);
        logFcent = log10(
            (fabs(troe_Tsss[15]) > 1.e-100 ? (1.-troe_a[15])*exp(-T[i]/troe_Tsss[15]) : 0.) 
            + (fabs(troe_Ts[15]) > 1.e-100 ? troe_a[15] * exp(-T[i]/troe_Ts[15]) : 0.) 
            + (troe_len[15] == 4 ? exp(-troe_Tss[15] * invT[i]) : 0.) );
        troe_c = -.4 - .67 * logFcent;
        troe_n = .75 - 1.27 * logFcent;
        troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
        F_troe = pow(10., logFcent / (1.0 + troe*troe));
        F *= F_troe;
        k_f *= F;
        q_f = phi_f * k_f;
        phi_r = sc[27*npt+i];
        Kc = Kc_s[15*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[9*npt+i] -= qdot;
        wdot[14*npt+i] -= qdot;
        wdot[27*npt+i] += qdot;

        /*reaction 17: CH2 + CO (+M) <=> CH2CO (+M) */
        phi_f = sc[10*npt+i]*sc[14*npt+i];
        alpha = mixture[i] + (TB[16][0] - 1)*sc[0*npt+i] + (TB[16][1] - 1)*sc[5*npt+i] + (TB[16][2] - 1)*sc[13*npt+i] + (TB[16][3] - 1)*sc[14*npt+i] + (TB[16][4] - 1)*sc[15*npt+i] + (TB[16][5] - 1)*sc[26*npt+i] + (TB[16][6] - 1)*sc[31*npt+i];
        k_f = k_f_s[16*npt+i];
        redP = alpha / k_f * phase_units[16] * low_A[16] * exp(low_beta[16] * tc[i] - activation_units[16] * low_Ea[16] * invT[i]);
        F = redP / (1 + redP);
        logPred = log10(redP);
        logFcent = log10(
            (fabs(troe_Tsss[16]) > 1.e-100 ? (1.-troe_a[16])*exp(-T[i]/troe_Tsss[16]) : 0.) 
            + (fabs(troe_Ts[16]) > 1.e-100 ? troe_a[16] * exp(-T[i]/troe_Ts[16]) : 0.) 
            + (troe_len[16] == 4 ? exp(-troe_Tss[16] * invT[i]) : 0.) );
        troe_c = -.4 - .67 * logFcent;
        troe_n = .75 - 1.27 * logFcent;
        troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
        F_troe = pow(10., logFcent / (1.0 + troe*troe));
        F *= F_troe;
        k_f *= F;
        q_f = phi_f * k_f;
        phi_r = sc[28*npt+i];
        Kc = Kc_s[16*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[10*npt+i] -= qdot;
        wdot[14*npt+i] -= qdot;
        wdot[28*npt+i] += qdot;

        /*reaction 18: CH2(S) + H2O (+M) <=> CH3OH (+M) */
        phi_f = sc[5*npt+i]*sc[11*npt+i];
        alpha = mixture[i] + (TB[17][0] - 1)*sc[0*npt+i] + (TB[17][1] - 1)*sc[5*npt+i] + (TB[17][2] - 1)*sc[13*npt+i] + (TB[17][3] - 1)*sc[14*npt+i] + (TB[17][4] - 1)*sc[15*npt+i] + (TB[17][5] - 1)*sc[26*npt+i];
        k_f = k_f_s[17*npt+i];
        redP = alpha / k_f * phase_units[17] * low_A[17] * exp(low_beta[17] * tc[i] - activation_units[17] * low_Ea[17] * invT[i]);
        F = redP / (1 + redP);
        logPred = log10(redP);
        logFcent = log10(
            (fabs(troe_Tsss[17]) > 1.e-100 ? (1.-troe_a[17])*exp(-T[i]/troe_Tsss[17]) : 0.) 
            + (fabs(troe_Ts[17]) > 1.e-100 ? troe_a[17] * exp(-T[i]/troe_Ts[17]) : 0.) 
            + (troe_len[17] == 4 ? exp(-troe_Tss[17] * invT[i]) : 0.) );
        troe_c = -.4 - .67 * logFcent;
        troe_n = .75 - 1.27 * logFcent;
        troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
        F_troe = pow(10., logFcent / (1.0 + troe*troe));
        F *= F_troe;
        k_f *= F;
        q_f = phi_f * k_f;
        phi_r = sc[20*npt+i];
        Kc = Kc_s[17*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[11*npt+i] -= qdot;
        wdot[20*npt+i] += qdot;

        /*reaction 19: 2.000000 CH3 (+M) <=> C2H6 (+M) */
        phi_f = pow(sc[12*npt+i], 2.000000);
        alpha = mixture[i] + (TB[18][0] - 1)*sc[0*npt+i] + (TB[18][1] - 1)*sc[5*npt+i] + (TB[18][2] - 1)*sc[13*npt+i] + (TB[18][3] - 1)*sc[14*npt+i] + (TB[18][4] - 1)*sc[15*npt+i] + (TB[18][5] - 1)*sc[26*npt+i] + (TB[18][6] - 1)*sc[31*npt+i];
        k_f = k_f_s[18*npt+i];
        redP = alpha / k_f * phase_units[18] * low_A[18] * exp(low_beta[18] * tc[i] - activation_units[18] * low_Ea[18] * invT[i]);
        F = redP / (1 + redP);
        logPred = log10(redP);
        logFcent = log10(
            (fabs(troe_Tsss[18]) > 1.e-100 ? (1.-troe_a[18])*exp(-T[i]/troe_Tsss[18]) : 0.) 
            + (fabs(troe_Ts[18]) > 1.e-100 ? troe_a[18] * exp(-T[i]/troe_Ts[18]) : 0.) 
            + (troe_len[18] == 4 ? exp(-troe_Tss[18] * invT[i]) : 0.) );
        troe_c = -.4 - .67 * logFcent;
        troe_n = .75 - 1.27 * logFcent;
        troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
        F_troe = pow(10., logFcent / (1.0 + troe*troe));
        F *= F_troe;
        k_f *= F;
        q_f = phi_f * k_f;
        phi_r = sc[26*npt+i];
        Kc = Kc_s[18*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[12*npt+i] -= 2.000000 * qdot;
        wdot[26*npt+i] += qdot;

        /*reaction 20: C2H4 (+M) <=> H2 + C2H2 (+M) */
        phi_f = sc[24*npt+i];
        alpha = mixture[i] + (TB[19][0] - 1)*sc[0*npt+i] + (TB[19][1] - 1)*sc[5*npt+i] + (TB[19][2] - 1)*sc[13*npt+i] + (TB[19][3] - 1)*sc[14*npt+i] + (TB[19][4] - 1)*sc[15*npt+i] + (TB[19][5] - 1)*sc[26*npt+i] + (TB[19][6] - 1)*sc[31*npt+i];
        k_f = k_f_s[19*npt+i];
        redP = alpha / k_f * phase_units[19] * low_A[19] * exp(low_beta[19] * tc[i] - activation_units[19] * low_Ea[19] * invT[i]);
        F = redP / (1 + redP);
        logPred = log10(redP);
        logFcent = log10(
            (fabs(troe_Tsss[19]) > 1.e-100 ? (1.-troe_a[19])*exp(-T[i]/troe_Tsss[19]) : 0.) 
            + (fabs(troe_Ts[19]) > 1.e-100 ? troe_a[19] * exp(-T[i]/troe_Ts[19]) : 0.) 
            + (troe_len[19] == 4 ? exp(-troe_Tss[19] * invT[i]) : 0.) );
        troe_c = -.4 - .67 * logFcent;
        troe_n = .75 - 1.27 * logFcent;
        troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
        F_troe = pow(10., logFcent / (1.0 + troe*troe));
        F *= F_troe;
        k_f *= F;
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[22*npt+i];
        Kc = Kc_s[19*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[22*npt+i] += qdot;
        wdot[24*npt+i] -= qdot;

        /*reaction 21: 2.000000 O + M <=> O2 + M */
        phi_f = pow(sc[2*npt+i], 2.000000);
        alpha = mixture[i] + (TB[20][0] - 1)*sc[0*npt+i] + (TB[20][1] - 1)*sc[5*npt+i] + (TB[20][2] - 1)*sc[13*npt+i] + (TB[20][3] - 1)*sc[14*npt+i] + (TB[20][4] - 1)*sc[15*npt+i] + (TB[20][5] - 1)*sc[26*npt+i] + (TB[20][6] - 1)*sc[31*npt+i];
        k_f = alpha * k_f_s[20*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i];
        Kc = Kc_s[20*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= 2.000000 * qdot;
        wdot[3*npt+i] += qdot;

        /*reaction 22: O + H + M <=> OH + M */
        phi_f = sc[1*npt+i]*sc[2*npt+i];
        alpha = mixture[i] + (TB[21][0] - 1)*sc[0*npt+i] + (TB[21][1] - 1)*sc[5*npt+i] + (TB[21][2] - 1)*sc[13*npt+i] + (TB[21][3] - 1)*sc[14*npt+i] + (TB[21][4] - 1)*sc[15*npt+i] + (TB[21][5] - 1)*sc[26*npt+i] + (TB[21][6] - 1)*sc[31*npt+i];
        k_f = alpha * k_f_s[21*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i];
        Kc = Kc_s[21*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;

        /*reaction 23: O + CO + M <=> CO2 + M */
        phi_f = sc[2*npt+i]*sc[14*npt+i];
        alpha = mixture[i] + (TB[22][0] - 1)*sc[0*npt+i] + (TB[22][1] - 1)*sc[3*npt+i] + (TB[22][2] - 1)*sc[5*npt+i] + (TB[22][3] - 1)*sc[13*npt+i] + (TB[22][4] - 1)*sc[14*npt+i] + (TB[22][5] - 1)*sc[15*npt+i] + (TB[22][6] - 1)*sc[26*npt+i] + (TB[22][7] - 1)*sc[31*npt+i];
        k_f = alpha * k_f_s[22*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[15*npt+i];
        Kc = Kc_s[22*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[14*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;

        /*reaction 24: H + O2 + M <=> HO2 + M */
        phi_f = sc[1*npt+i]*sc[3*npt+i];
        alpha = mixture[i] + (TB[23][0] - 1)*sc[3*npt+i] + (TB[23][1] - 1)*sc[5*npt+i] + (TB[23][2] - 1)*sc[14*npt+i] + (TB[23][3] - 1)*sc[15*npt+i] + (TB[23][4] - 1)*sc[26*npt+i] + (TB[23][5] - 1)*sc[30*npt+i] + (TB[23][6] - 1)*sc[31*npt+i];
        k_f = alpha * k_f_s[23*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i];
        Kc = Kc_s[23*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;

        /*reaction 25: 2.000000 H + M <=> H2 + M */
        phi_f = pow(sc[1*npt+i], 2.000000);
        alpha = mixture[i] + (TB[24][0] - 1)*sc[0*npt+i] + (TB[24][1] - 1)*sc[5*npt+i] + (TB[24][2] - 1)*sc[13*npt+i] + (TB[24][3] - 1)*sc[15*npt+i] + (TB[24][4] - 1)*sc[26*npt+i] + (TB[24][5] - 1)*sc[31*npt+i];
        k_f = alpha * k_f_s[24*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i];
        Kc = Kc_s[24*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= 2.000000 * qdot;

        /*reaction 26: H + OH + M <=> H2O + M */
        phi_f = sc[1*npt+i]*sc[4*npt+i];
        alpha = mixture[i] + (TB[25][0] - 1)*sc[0*npt+i] + (TB[25][1] - 1)*sc[5*npt+i] + (TB[25][2] - 1)*sc[13*npt+i] + (TB[25][3] - 1)*sc[26*npt+i] + (TB[25][4] - 1)*sc[31*npt+i];
        k_f = alpha * k_f_s[25*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i];
        Kc = Kc_s[25*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;

        /*reaction 27: HCO + M <=> H + CO + M */
        phi_f = sc[16*npt+i];
        alpha = mixture[i] + (TB[26][0] - 1)*sc[0*npt+i] + (TB[26][1] - 1)*sc[5*npt+i] + (TB[26][2] - 1)*sc[13*npt+i] + (TB[26][3] - 1)*sc[14*npt+i] + (TB[26][4] - 1)*sc[15*npt+i] + (TB[26][5] - 1)*sc[26*npt+i];
        k_f = alpha * k_f_s[26*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[14*npt+i];
        Kc = Kc_s[26*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[14*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;

        /*reaction 28: O + H2 <=> H + OH */
        phi_f = sc[0*npt+i]*sc[2*npt+i];
        k_f = k_f_s[27*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[4*npt+i];
        Kc = Kc_s[27*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;

        /*reaction 29: O + HO2 <=> OH + O2 */
        phi_f = sc[2*npt+i]*sc[6*npt+i];
        k_f = k_f_s[28*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[4*npt+i];
        Kc = Kc_s[28*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[4*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 30: O + H2O2 <=> OH + HO2 */
        phi_f = sc[2*npt+i]*sc[7*npt+i];
        k_f = k_f_s[29*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[6*npt+i];
        Kc = Kc_s[29*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 31: O + CH <=> H + CO */
        phi_f = sc[2*npt+i]*sc[9*npt+i];
        k_f = k_f_s[30*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[14*npt+i];
        Kc = Kc_s[30*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[9*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;

        /*reaction 32: O + CH2 <=> H + HCO */
        phi_f = sc[2*npt+i]*sc[10*npt+i];
        k_f = k_f_s[31*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[16*npt+i];
        Kc = Kc_s[31*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[10*npt+i] -= qdot;
        wdot[16*npt+i] += qdot;

        /*reaction 33: O + CH2(S) <=> H2 + CO */
        phi_f = sc[2*npt+i]*sc[11*npt+i];
        k_f = k_f_s[32*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[14*npt+i];
        Kc = Kc_s[32*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[11*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;

        /*reaction 34: O + CH2(S) <=> H + HCO */
        phi_f = sc[2*npt+i]*sc[11*npt+i];
        k_f = k_f_s[33*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[16*npt+i];
        Kc = Kc_s[33*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[11*npt+i] -= qdot;
        wdot[16*npt+i] += qdot;

        /*reaction 35: O + CH3 <=> H + CH2O */
        phi_f = sc[2*npt+i]*sc[12*npt+i];
        k_f = k_f_s[34*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[17*npt+i];
        Kc = Kc_s[34*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[12*npt+i] -= qdot;
        wdot[17*npt+i] += qdot;

        /*reaction 36: O + CH4 <=> OH + CH3 */
        phi_f = sc[2*npt+i]*sc[13*npt+i];
        k_f = k_f_s[35*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[12*npt+i];
        Kc = Kc_s[35*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[12*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;

        /*reaction 37: O + HCO <=> OH + CO */
        phi_f = sc[2*npt+i]*sc[16*npt+i];
        k_f = k_f_s[36*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[14*npt+i];
        Kc = Kc_s[36*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[14*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;

        /*reaction 38: O + HCO <=> H + CO2 */
        phi_f = sc[2*npt+i]*sc[16*npt+i];
        k_f = k_f_s[37*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[15*npt+i];
        Kc = Kc_s[37*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;

        /*reaction 39: O + CH2O <=> OH + HCO */
        phi_f = sc[2*npt+i]*sc[17*npt+i];
        k_f = k_f_s[38*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[16*npt+i];
        Kc = Kc_s[38*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[16*npt+i] += qdot;
        wdot[17*npt+i] -= qdot;

        /*reaction 40: O + CH2OH <=> OH + CH2O */
        phi_f = sc[2*npt+i]*sc[18*npt+i];
        k_f = k_f_s[39*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[17*npt+i];
        Kc = Kc_s[39*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[17*npt+i] += qdot;
        wdot[18*npt+i] -= qdot;

        /*reaction 41: O + CH3O <=> OH + CH2O */
        phi_f = sc[2*npt+i]*sc[19*npt+i];
        k_f = k_f_s[40*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[17*npt+i];
        Kc = Kc_s[40*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[17*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 42: O + CH3OH <=> OH + CH2OH */
        phi_f = sc[2*npt+i]*sc[20*npt+i];
        k_f = k_f_s[41*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[18*npt+i];
        Kc = Kc_s[41*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[18*npt+i] += qdot;
        wdot[20*npt+i] -= qdot;

        /*reaction 43: O + CH3OH <=> OH + CH3O */
        phi_f = sc[2*npt+i]*sc[20*npt+i];
        k_f = k_f_s[42*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[19*npt+i];
        Kc = Kc_s[42*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[19*npt+i] += qdot;
        wdot[20*npt+i] -= qdot;

        /*reaction 44: O + C2H <=> CH + CO */
        phi_f = sc[2*npt+i]*sc[21*npt+i];
        k_f = k_f_s[43*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[9*npt+i]*sc[14*npt+i];
        Kc = Kc_s[43*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;
        wdot[14*npt+i] += qdot;
        wdot[21*npt+i] -= qdot;

        /*reaction 45: O + C2H2 <=> H + HCCO */
        phi_f = sc[2*npt+i]*sc[22*npt+i];
        k_f = k_f_s[44*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[27*npt+i];
        Kc = Kc_s[44*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[22*npt+i] -= qdot;
        wdot[27*npt+i] += qdot;

        /*reaction 46: O + C2H2 <=> OH + C2H */
        phi_f = sc[2*npt+i]*sc[22*npt+i];
        k_f = k_f_s[45*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[21*npt+i];
        Kc = Kc_s[45*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[21*npt+i] += qdot;
        wdot[22*npt+i] -= qdot;

        /*reaction 47: O + C2H2 <=> CO + CH2 */
        phi_f = sc[2*npt+i]*sc[22*npt+i];
        k_f = k_f_s[46*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[10*npt+i]*sc[14*npt+i];
        Kc = Kc_s[46*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;
        wdot[14*npt+i] += qdot;
        wdot[22*npt+i] -= qdot;

        /*reaction 48: O + C2H3 <=> H + CH2CO */
        phi_f = sc[2*npt+i]*sc[23*npt+i];
        k_f = k_f_s[47*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[28*npt+i];
        Kc = Kc_s[47*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[23*npt+i] -= qdot;
        wdot[28*npt+i] += qdot;

        /*reaction 49: O + C2H4 <=> CH3 + HCO */
        phi_f = sc[2*npt+i]*sc[24*npt+i];
        k_f = k_f_s[48*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[12*npt+i]*sc[16*npt+i];
        Kc = Kc_s[48*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;
        wdot[16*npt+i] += qdot;
        wdot[24*npt+i] -= qdot;

        /*reaction 50: O + C2H5 <=> CH3 + CH2O */
        phi_f = sc[2*npt+i]*sc[25*npt+i];
        k_f = k_f_s[49*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[12*npt+i]*sc[17*npt+i];
        Kc = Kc_s[49*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;
        wdot[17*npt+i] += qdot;
        wdot[25*npt+i] -= qdot;
    }
}

void vcomp_wdot_51_100(int npt, amrex::Real *  wdot, amrex::Real *  mixture, amrex::Real *  sc,
		amrex::Real *  k_f_s, amrex::Real *  Kc_s,
		amrex::Real *  tc, amrex::Real *  invT, amrex::Real *  T)
{
    for (int i=0; i<npt; i++) {
        amrex::Real qdot, q_f, q_r, phi_f, phi_r, k_f, k_r, Kc;
        amrex::Real alpha;

        /*reaction 51: O + C2H6 <=> OH + C2H5 */
        phi_f = sc[2*npt+i]*sc[26*npt+i];
        k_f = k_f_s[50*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[25*npt+i];
        Kc = Kc_s[50*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[25*npt+i] += qdot;
        wdot[26*npt+i] -= qdot;

        /*reaction 52: O + HCCO <=> H + 2.000000 CO */
        phi_f = sc[2*npt+i]*sc[27*npt+i];
        k_f = k_f_s[51*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*pow(sc[14*npt+i], 2.000000);
        Kc = Kc_s[51*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[14*npt+i] += 2.000000 * qdot;
        wdot[27*npt+i] -= qdot;

        /*reaction 53: O + CH2CO <=> OH + HCCO */
        phi_f = sc[2*npt+i]*sc[28*npt+i];
        k_f = k_f_s[52*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[27*npt+i];
        Kc = Kc_s[52*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[27*npt+i] += qdot;
        wdot[28*npt+i] -= qdot;

        /*reaction 54: O + CH2CO <=> CH2 + CO2 */
        phi_f = sc[2*npt+i]*sc[28*npt+i];
        k_f = k_f_s[53*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[10*npt+i]*sc[15*npt+i];
        Kc = Kc_s[53*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;
        wdot[15*npt+i] += qdot;
        wdot[28*npt+i] -= qdot;

        /*reaction 55: O2 + CO <=> O + CO2 */
        phi_f = sc[3*npt+i]*sc[14*npt+i];
        k_f = k_f_s[54*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[15*npt+i];
        Kc = Kc_s[54*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[14*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;

        /*reaction 56: O2 + CH2O <=> HO2 + HCO */
        phi_f = sc[3*npt+i]*sc[17*npt+i];
        k_f = k_f_s[55*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[16*npt+i];
        Kc = Kc_s[55*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[16*npt+i] += qdot;
        wdot[17*npt+i] -= qdot;

        /*reaction 57: H + 2.000000 O2 <=> HO2 + O2 */
        phi_f = sc[1*npt+i]*pow(sc[3*npt+i], 2.000000);
        k_f = k_f_s[56*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[6*npt+i];
        Kc = Kc_s[56*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] -= 2.000000 * qdot;
        wdot[3*npt+i] += qdot;
        wdot[6*npt+i] += qdot;

        /*reaction 58: H + O2 + H2O <=> HO2 + H2O */
        phi_f = sc[1*npt+i]*sc[3*npt+i]*sc[5*npt+i];
        k_f = k_f_s[57*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[6*npt+i];
        Kc = Kc_s[57*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] -= qdot;
        wdot[5*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] += qdot;

        /*reaction 59: H + O2 + N2 <=> HO2 + N2 */
        phi_f = sc[1*npt+i]*sc[3*npt+i]*sc[30*npt+i];
        k_f = k_f_s[58*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[30*npt+i];
        Kc = Kc_s[58*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[30*npt+i] -= qdot;
        wdot[30*npt+i] += qdot;

        /*reaction 60: H + O2 + AR <=> HO2 + AR */
        phi_f = sc[1*npt+i]*sc[3*npt+i]*sc[31*npt+i];
        k_f = k_f_s[59*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[31*npt+i];
        Kc = Kc_s[59*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[31*npt+i] -= qdot;
        wdot[31*npt+i] += qdot;

        /*reaction 61: H + O2 <=> O + OH */
        phi_f = sc[1*npt+i]*sc[3*npt+i];
        k_f = k_f_s[60*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[4*npt+i];
        Kc = Kc_s[60*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;

        /*reaction 62: 2.000000 H + H2 <=> 2.000000 H2 */
        phi_f = sc[0*npt+i]*pow(sc[1*npt+i], 2.000000);
        k_f = k_f_s[61*npt+i];
        q_f = phi_f * k_f;
        phi_r = pow(sc[0*npt+i], 2.000000);
        Kc = Kc_s[61*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[0*npt+i] += 2.000000 * qdot;
        wdot[1*npt+i] -= 2.000000 * qdot;

        /*reaction 63: 2.000000 H + H2O <=> H2 + H2O */
        phi_f = pow(sc[1*npt+i], 2.000000)*sc[5*npt+i];
        k_f = k_f_s[62*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[5*npt+i];
        Kc = Kc_s[62*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= 2.000000 * qdot;
        wdot[5*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;

        /*reaction 64: 2.000000 H + CO2 <=> H2 + CO2 */
        phi_f = pow(sc[1*npt+i], 2.000000)*sc[15*npt+i];
        k_f = k_f_s[63*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[15*npt+i];
        Kc = Kc_s[63*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= 2.000000 * qdot;
        wdot[15*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;

        /*reaction 65: H + HO2 <=> O + H2O */
        phi_f = sc[1*npt+i]*sc[6*npt+i];
        k_f = k_f_s[64*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[5*npt+i];
        Kc = Kc_s[64*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 66: H + HO2 <=> O2 + H2 */
        phi_f = sc[1*npt+i]*sc[6*npt+i];
        k_f = k_f_s[65*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[3*npt+i];
        Kc = Kc_s[65*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 67: H + HO2 <=> 2.000000 OH */
        phi_f = sc[1*npt+i]*sc[6*npt+i];
        k_f = k_f_s[66*npt+i];
        q_f = phi_f * k_f;
        phi_r = pow(sc[4*npt+i], 2.000000);
        Kc = Kc_s[66*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[4*npt+i] += 2.000000 * qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 68: H + H2O2 <=> HO2 + H2 */
        phi_f = sc[1*npt+i]*sc[7*npt+i];
        k_f = k_f_s[67*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[6*npt+i];
        Kc = Kc_s[67*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 69: H + H2O2 <=> OH + H2O */
        phi_f = sc[1*npt+i]*sc[7*npt+i];
        k_f = k_f_s[68*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[5*npt+i];
        Kc = Kc_s[68*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[5*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 70: H + CH <=> C + H2 */
        phi_f = sc[1*npt+i]*sc[9*npt+i];
        k_f = k_f_s[69*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[8*npt+i];
        Kc = Kc_s[69*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[9*npt+i] -= qdot;

        /*reaction 71: H + CH2(S) <=> CH + H2 */
        phi_f = sc[1*npt+i]*sc[11*npt+i];
        k_f = k_f_s[70*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[9*npt+i];
        Kc = Kc_s[70*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;

        /*reaction 72: H + CH4 <=> CH3 + H2 */
        phi_f = sc[1*npt+i]*sc[13*npt+i];
        k_f = k_f_s[71*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[12*npt+i];
        Kc = Kc_s[71*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;

        /*reaction 73: H + HCO <=> H2 + CO */
        phi_f = sc[1*npt+i]*sc[16*npt+i];
        k_f = k_f_s[72*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[14*npt+i];
        Kc = Kc_s[72*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;

        /*reaction 74: H + CH2O <=> HCO + H2 */
        phi_f = sc[1*npt+i]*sc[17*npt+i];
        k_f = k_f_s[73*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[16*npt+i];
        Kc = Kc_s[73*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[16*npt+i] += qdot;
        wdot[17*npt+i] -= qdot;

        /*reaction 75: H + CH2OH <=> H2 + CH2O */
        phi_f = sc[1*npt+i]*sc[18*npt+i];
        k_f = k_f_s[74*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[17*npt+i];
        Kc = Kc_s[74*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[17*npt+i] += qdot;
        wdot[18*npt+i] -= qdot;

        /*reaction 76: H + CH2OH <=> OH + CH3 */
        phi_f = sc[1*npt+i]*sc[18*npt+i];
        k_f = k_f_s[75*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[12*npt+i];
        Kc = Kc_s[75*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[12*npt+i] += qdot;
        wdot[18*npt+i] -= qdot;

        /*reaction 77: H + CH2OH <=> CH2(S) + H2O */
        phi_f = sc[1*npt+i]*sc[18*npt+i];
        k_f = k_f_s[76*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[11*npt+i];
        Kc = Kc_s[76*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[11*npt+i] += qdot;
        wdot[18*npt+i] -= qdot;

        /*reaction 78: H + CH3O <=> H + CH2OH */
        phi_f = sc[1*npt+i]*sc[19*npt+i];
        k_f = k_f_s[77*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[18*npt+i];
        Kc = Kc_s[77*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[18*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 79: H + CH3O <=> H2 + CH2O */
        phi_f = sc[1*npt+i]*sc[19*npt+i];
        k_f = k_f_s[78*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[17*npt+i];
        Kc = Kc_s[78*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[17*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 80: H + CH3O <=> OH + CH3 */
        phi_f = sc[1*npt+i]*sc[19*npt+i];
        k_f = k_f_s[79*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[12*npt+i];
        Kc = Kc_s[79*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[12*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 81: H + CH3O <=> CH2(S) + H2O */
        phi_f = sc[1*npt+i]*sc[19*npt+i];
        k_f = k_f_s[80*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[11*npt+i];
        Kc = Kc_s[80*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[11*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 82: H + CH3OH <=> CH2OH + H2 */
        phi_f = sc[1*npt+i]*sc[20*npt+i];
        k_f = k_f_s[81*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[18*npt+i];
        Kc = Kc_s[81*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[18*npt+i] += qdot;
        wdot[20*npt+i] -= qdot;

        /*reaction 83: H + CH3OH <=> CH3O + H2 */
        phi_f = sc[1*npt+i]*sc[20*npt+i];
        k_f = k_f_s[82*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[19*npt+i];
        Kc = Kc_s[82*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[19*npt+i] += qdot;
        wdot[20*npt+i] -= qdot;

        /*reaction 84: H + C2H3 <=> H2 + C2H2 */
        phi_f = sc[1*npt+i]*sc[23*npt+i];
        k_f = k_f_s[83*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[22*npt+i];
        Kc = Kc_s[83*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[22*npt+i] += qdot;
        wdot[23*npt+i] -= qdot;

        /*reaction 85: H + C2H4 <=> C2H3 + H2 */
        phi_f = sc[1*npt+i]*sc[24*npt+i];
        k_f = k_f_s[84*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[23*npt+i];
        Kc = Kc_s[84*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[23*npt+i] += qdot;
        wdot[24*npt+i] -= qdot;

        /*reaction 86: H + C2H5 <=> H2 + C2H4 */
        phi_f = sc[1*npt+i]*sc[25*npt+i];
        k_f = k_f_s[85*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[24*npt+i];
        Kc = Kc_s[85*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[24*npt+i] += qdot;
        wdot[25*npt+i] -= qdot;

        /*reaction 87: H + C2H6 <=> C2H5 + H2 */
        phi_f = sc[1*npt+i]*sc[26*npt+i];
        k_f = k_f_s[86*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[25*npt+i];
        Kc = Kc_s[86*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[25*npt+i] += qdot;
        wdot[26*npt+i] -= qdot;

        /*reaction 88: H + HCCO <=> CH2(S) + CO */
        phi_f = sc[1*npt+i]*sc[27*npt+i];
        k_f = k_f_s[87*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[11*npt+i]*sc[14*npt+i];
        Kc = Kc_s[87*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;
        wdot[14*npt+i] += qdot;
        wdot[27*npt+i] -= qdot;

        /*reaction 89: H + CH2CO <=> HCCO + H2 */
        phi_f = sc[1*npt+i]*sc[28*npt+i];
        k_f = k_f_s[88*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[27*npt+i];
        Kc = Kc_s[88*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[27*npt+i] += qdot;
        wdot[28*npt+i] -= qdot;

        /*reaction 90: H + CH2CO <=> CH3 + CO */
        phi_f = sc[1*npt+i]*sc[28*npt+i];
        k_f = k_f_s[89*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[12*npt+i]*sc[14*npt+i];
        Kc = Kc_s[89*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;
        wdot[14*npt+i] += qdot;
        wdot[28*npt+i] -= qdot;

        /*reaction 91: H + HCCOH <=> H + CH2CO */
        phi_f = sc[1*npt+i]*sc[29*npt+i];
        k_f = k_f_s[90*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[28*npt+i];
        Kc = Kc_s[90*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[28*npt+i] += qdot;
        wdot[29*npt+i] -= qdot;

        /*reaction 92: OH + H2 <=> H + H2O */
        phi_f = sc[0*npt+i]*sc[4*npt+i];
        k_f = k_f_s[91*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[5*npt+i];
        Kc = Kc_s[91*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;

        /*reaction 93: 2.000000 OH <=> O + H2O */
        phi_f = pow(sc[4*npt+i], 2.000000);
        k_f = k_f_s[92*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[5*npt+i];
        Kc = Kc_s[92*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[4*npt+i] -= 2.000000 * qdot;
        wdot[5*npt+i] += qdot;

        /*reaction 94: OH + HO2 <=> O2 + H2O */
        phi_f = sc[4*npt+i]*sc[6*npt+i];
        k_f = k_f_s[93*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[5*npt+i];
        Kc = Kc_s[93*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 95: OH + H2O2 <=> HO2 + H2O */
        phi_f = sc[4*npt+i]*sc[7*npt+i];
        k_f = k_f_s[94*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[6*npt+i];
        Kc = Kc_s[94*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 96: OH + H2O2 <=> HO2 + H2O */
        phi_f = sc[4*npt+i]*sc[7*npt+i];
        k_f = k_f_s[95*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[6*npt+i];
        Kc = Kc_s[95*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 97: OH + C <=> H + CO */
        phi_f = sc[4*npt+i]*sc[8*npt+i];
        k_f = k_f_s[96*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[14*npt+i];
        Kc = Kc_s[96*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[8*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;

        /*reaction 98: OH + CH <=> H + HCO */
        phi_f = sc[4*npt+i]*sc[9*npt+i];
        k_f = k_f_s[97*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[16*npt+i];
        Kc = Kc_s[97*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[9*npt+i] -= qdot;
        wdot[16*npt+i] += qdot;

        /*reaction 99: OH + CH2 <=> H + CH2O */
        phi_f = sc[4*npt+i]*sc[10*npt+i];
        k_f = k_f_s[98*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[17*npt+i];
        Kc = Kc_s[98*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[10*npt+i] -= qdot;
        wdot[17*npt+i] += qdot;

        /*reaction 100: OH + CH2 <=> CH + H2O */
        phi_f = sc[4*npt+i]*sc[10*npt+i];
        k_f = k_f_s[99*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[9*npt+i];
        Kc = Kc_s[99*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[9*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;
    }
}

void vcomp_wdot_101_150(int npt, amrex::Real *  wdot, amrex::Real *  mixture, amrex::Real *  sc,
		amrex::Real *  k_f_s, amrex::Real *  Kc_s,
		amrex::Real *  tc, amrex::Real *  invT, amrex::Real *  T)
{
    for (int i=0; i<npt; i++) {
        amrex::Real qdot, q_f, q_r, phi_f, phi_r, k_f, k_r, Kc;
        amrex::Real alpha;

        /*reaction 101: OH + CH2(S) <=> H + CH2O */
        phi_f = sc[4*npt+i]*sc[11*npt+i];
        k_f = k_f_s[100*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[17*npt+i];
        Kc = Kc_s[100*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[11*npt+i] -= qdot;
        wdot[17*npt+i] += qdot;

        /*reaction 102: OH + CH3 <=> CH2 + H2O */
        phi_f = sc[4*npt+i]*sc[12*npt+i];
        k_f = k_f_s[101*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[10*npt+i];
        Kc = Kc_s[101*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[10*npt+i] += qdot;
        wdot[12*npt+i] -= qdot;

        /*reaction 103: OH + CH3 <=> CH2(S) + H2O */
        phi_f = sc[4*npt+i]*sc[12*npt+i];
        k_f = k_f_s[102*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[11*npt+i];
        Kc = Kc_s[102*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[11*npt+i] += qdot;
        wdot[12*npt+i] -= qdot;

        /*reaction 104: OH + CH4 <=> CH3 + H2O */
        phi_f = sc[4*npt+i]*sc[13*npt+i];
        k_f = k_f_s[103*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[12*npt+i];
        Kc = Kc_s[103*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[12*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;

        /*reaction 105: OH + CO <=> H + CO2 */
        phi_f = sc[4*npt+i]*sc[14*npt+i];
        k_f = k_f_s[104*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[15*npt+i];
        Kc = Kc_s[104*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[14*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;

        /*reaction 106: OH + HCO <=> H2O + CO */
        phi_f = sc[4*npt+i]*sc[16*npt+i];
        k_f = k_f_s[105*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[14*npt+i];
        Kc = Kc_s[105*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[14*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;

        /*reaction 107: OH + CH2O <=> HCO + H2O */
        phi_f = sc[4*npt+i]*sc[17*npt+i];
        k_f = k_f_s[106*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[16*npt+i];
        Kc = Kc_s[106*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[16*npt+i] += qdot;
        wdot[17*npt+i] -= qdot;

        /*reaction 108: OH + CH2OH <=> H2O + CH2O */
        phi_f = sc[4*npt+i]*sc[18*npt+i];
        k_f = k_f_s[107*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[17*npt+i];
        Kc = Kc_s[107*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[17*npt+i] += qdot;
        wdot[18*npt+i] -= qdot;

        /*reaction 109: OH + CH3O <=> H2O + CH2O */
        phi_f = sc[4*npt+i]*sc[19*npt+i];
        k_f = k_f_s[108*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[17*npt+i];
        Kc = Kc_s[108*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[17*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 110: OH + CH3OH <=> CH2OH + H2O */
        phi_f = sc[4*npt+i]*sc[20*npt+i];
        k_f = k_f_s[109*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[18*npt+i];
        Kc = Kc_s[109*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[18*npt+i] += qdot;
        wdot[20*npt+i] -= qdot;

        /*reaction 111: OH + CH3OH <=> CH3O + H2O */
        phi_f = sc[4*npt+i]*sc[20*npt+i];
        k_f = k_f_s[110*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[19*npt+i];
        Kc = Kc_s[110*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[19*npt+i] += qdot;
        wdot[20*npt+i] -= qdot;

        /*reaction 112: OH + C2H <=> H + HCCO */
        phi_f = sc[4*npt+i]*sc[21*npt+i];
        k_f = k_f_s[111*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[27*npt+i];
        Kc = Kc_s[111*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[21*npt+i] -= qdot;
        wdot[27*npt+i] += qdot;

        /*reaction 113: OH + C2H2 <=> H + CH2CO */
        phi_f = sc[4*npt+i]*sc[22*npt+i];
        k_f = k_f_s[112*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[28*npt+i];
        Kc = Kc_s[112*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[22*npt+i] -= qdot;
        wdot[28*npt+i] += qdot;

        /*reaction 114: OH + C2H2 <=> H + HCCOH */
        phi_f = sc[4*npt+i]*sc[22*npt+i];
        k_f = k_f_s[113*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[29*npt+i];
        Kc = Kc_s[113*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[22*npt+i] -= qdot;
        wdot[29*npt+i] += qdot;

        /*reaction 115: OH + C2H2 <=> C2H + H2O */
        phi_f = sc[4*npt+i]*sc[22*npt+i];
        k_f = k_f_s[114*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[21*npt+i];
        Kc = Kc_s[114*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[21*npt+i] += qdot;
        wdot[22*npt+i] -= qdot;

        /*reaction 116: OH + C2H2 <=> CH3 + CO */
        phi_f = sc[4*npt+i]*sc[22*npt+i];
        k_f = k_f_s[115*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[12*npt+i]*sc[14*npt+i];
        Kc = Kc_s[115*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;
        wdot[14*npt+i] += qdot;
        wdot[22*npt+i] -= qdot;

        /*reaction 117: OH + C2H3 <=> H2O + C2H2 */
        phi_f = sc[4*npt+i]*sc[23*npt+i];
        k_f = k_f_s[116*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[22*npt+i];
        Kc = Kc_s[116*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[22*npt+i] += qdot;
        wdot[23*npt+i] -= qdot;

        /*reaction 118: OH + C2H4 <=> C2H3 + H2O */
        phi_f = sc[4*npt+i]*sc[24*npt+i];
        k_f = k_f_s[117*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[23*npt+i];
        Kc = Kc_s[117*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[23*npt+i] += qdot;
        wdot[24*npt+i] -= qdot;

        /*reaction 119: OH + C2H6 <=> C2H5 + H2O */
        phi_f = sc[4*npt+i]*sc[26*npt+i];
        k_f = k_f_s[118*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[25*npt+i];
        Kc = Kc_s[118*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[25*npt+i] += qdot;
        wdot[26*npt+i] -= qdot;

        /*reaction 120: OH + CH2CO <=> HCCO + H2O */
        phi_f = sc[4*npt+i]*sc[28*npt+i];
        k_f = k_f_s[119*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[27*npt+i];
        Kc = Kc_s[119*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[27*npt+i] += qdot;
        wdot[28*npt+i] -= qdot;

        /*reaction 121: 2.000000 HO2 <=> O2 + H2O2 */
        phi_f = pow(sc[6*npt+i], 2.000000);
        k_f = k_f_s[120*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[7*npt+i];
        Kc = Kc_s[120*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[6*npt+i] -= 2.000000 * qdot;
        wdot[7*npt+i] += qdot;

        /*reaction 122: 2.000000 HO2 <=> O2 + H2O2 */
        phi_f = pow(sc[6*npt+i], 2.000000);
        k_f = k_f_s[121*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[7*npt+i];
        Kc = Kc_s[121*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[6*npt+i] -= 2.000000 * qdot;
        wdot[7*npt+i] += qdot;

        /*reaction 123: HO2 + CH2 <=> OH + CH2O */
        phi_f = sc[6*npt+i]*sc[10*npt+i];
        k_f = k_f_s[122*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[17*npt+i];
        Kc = Kc_s[122*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[10*npt+i] -= qdot;
        wdot[17*npt+i] += qdot;

        /*reaction 124: HO2 + CH3 <=> O2 + CH4 */
        phi_f = sc[6*npt+i]*sc[12*npt+i];
        k_f = k_f_s[123*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[13*npt+i];
        Kc = Kc_s[123*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[12*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;

        /*reaction 125: HO2 + CH3 <=> OH + CH3O */
        phi_f = sc[6*npt+i]*sc[12*npt+i];
        k_f = k_f_s[124*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[19*npt+i];
        Kc = Kc_s[124*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[12*npt+i] -= qdot;
        wdot[19*npt+i] += qdot;

        /*reaction 126: HO2 + CO <=> OH + CO2 */
        phi_f = sc[6*npt+i]*sc[14*npt+i];
        k_f = k_f_s[125*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[15*npt+i];
        Kc = Kc_s[125*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[14*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;

        /*reaction 127: HO2 + CH2O <=> HCO + H2O2 */
        phi_f = sc[6*npt+i]*sc[17*npt+i];
        k_f = k_f_s[126*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[16*npt+i];
        Kc = Kc_s[126*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[16*npt+i] += qdot;
        wdot[17*npt+i] -= qdot;

        /*reaction 128: C + O2 <=> O + CO */
        phi_f = sc[3*npt+i]*sc[8*npt+i];
        k_f = k_f_s[127*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[14*npt+i];
        Kc = Kc_s[127*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[8*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;

        /*reaction 129: C + CH2 <=> H + C2H */
        phi_f = sc[8*npt+i]*sc[10*npt+i];
        k_f = k_f_s[128*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[21*npt+i];
        Kc = Kc_s[128*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[10*npt+i] -= qdot;
        wdot[21*npt+i] += qdot;

        /*reaction 130: C + CH3 <=> H + C2H2 */
        phi_f = sc[8*npt+i]*sc[12*npt+i];
        k_f = k_f_s[129*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[22*npt+i];
        Kc = Kc_s[129*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[12*npt+i] -= qdot;
        wdot[22*npt+i] += qdot;

        /*reaction 131: CH + O2 <=> O + HCO */
        phi_f = sc[3*npt+i]*sc[9*npt+i];
        k_f = k_f_s[130*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[16*npt+i];
        Kc = Kc_s[130*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[9*npt+i] -= qdot;
        wdot[16*npt+i] += qdot;

        /*reaction 132: CH + H2 <=> H + CH2 */
        phi_f = sc[0*npt+i]*sc[9*npt+i];
        k_f = k_f_s[131*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[10*npt+i];
        Kc = Kc_s[131*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[9*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;

        /*reaction 133: CH + H2O <=> H + CH2O */
        phi_f = sc[5*npt+i]*sc[9*npt+i];
        k_f = k_f_s[132*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[17*npt+i];
        Kc = Kc_s[132*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[9*npt+i] -= qdot;
        wdot[17*npt+i] += qdot;

        /*reaction 134: CH + CH2 <=> H + C2H2 */
        phi_f = sc[9*npt+i]*sc[10*npt+i];
        k_f = k_f_s[133*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[22*npt+i];
        Kc = Kc_s[133*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[9*npt+i] -= qdot;
        wdot[10*npt+i] -= qdot;
        wdot[22*npt+i] += qdot;

        /*reaction 135: CH + CH3 <=> H + C2H3 */
        phi_f = sc[9*npt+i]*sc[12*npt+i];
        k_f = k_f_s[134*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[23*npt+i];
        Kc = Kc_s[134*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[9*npt+i] -= qdot;
        wdot[12*npt+i] -= qdot;
        wdot[23*npt+i] += qdot;

        /*reaction 136: CH + CH4 <=> H + C2H4 */
        phi_f = sc[9*npt+i]*sc[13*npt+i];
        k_f = k_f_s[135*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[24*npt+i];
        Kc = Kc_s[135*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[9*npt+i] -= qdot;
        wdot[13*npt+i] -= qdot;
        wdot[24*npt+i] += qdot;

        /*reaction 137: CH + CO2 <=> HCO + CO */
        phi_f = sc[9*npt+i]*sc[15*npt+i];
        k_f = k_f_s[136*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[14*npt+i]*sc[16*npt+i];
        Kc = Kc_s[136*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[9*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;
        wdot[15*npt+i] -= qdot;
        wdot[16*npt+i] += qdot;

        /*reaction 138: CH + CH2O <=> H + CH2CO */
        phi_f = sc[9*npt+i]*sc[17*npt+i];
        k_f = k_f_s[137*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[28*npt+i];
        Kc = Kc_s[137*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[9*npt+i] -= qdot;
        wdot[17*npt+i] -= qdot;
        wdot[28*npt+i] += qdot;

        /*reaction 139: CH + HCCO <=> CO + C2H2 */
        phi_f = sc[9*npt+i]*sc[27*npt+i];
        k_f = k_f_s[138*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[14*npt+i]*sc[22*npt+i];
        Kc = Kc_s[138*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[9*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;
        wdot[22*npt+i] += qdot;
        wdot[27*npt+i] -= qdot;

        /*reaction 140: CH2 + O2 <=> OH + HCO */
        phi_f = sc[3*npt+i]*sc[10*npt+i];
        k_f = k_f_s[139*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[16*npt+i];
        Kc = Kc_s[139*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;
        wdot[16*npt+i] += qdot;

        /*reaction 141: CH2 + H2 <=> H + CH3 */
        phi_f = sc[0*npt+i]*sc[10*npt+i];
        k_f = k_f_s[140*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[12*npt+i];
        Kc = Kc_s[140*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;

        /*reaction 142: 2.000000 CH2 <=> H2 + C2H2 */
        phi_f = pow(sc[10*npt+i], 2.000000);
        k_f = k_f_s[141*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[22*npt+i];
        Kc = Kc_s[141*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[10*npt+i] -= 2.000000 * qdot;
        wdot[22*npt+i] += qdot;

        /*reaction 143: CH2 + CH3 <=> H + C2H4 */
        phi_f = sc[10*npt+i]*sc[12*npt+i];
        k_f = k_f_s[142*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[24*npt+i];
        Kc = Kc_s[142*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;
        wdot[12*npt+i] -= qdot;
        wdot[24*npt+i] += qdot;

        /*reaction 144: CH2 + CH4 <=> 2.000000 CH3 */
        phi_f = sc[10*npt+i]*sc[13*npt+i];
        k_f = k_f_s[143*npt+i];
        q_f = phi_f * k_f;
        phi_r = pow(sc[12*npt+i], 2.000000);
        Kc = Kc_s[143*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[10*npt+i] -= qdot;
        wdot[12*npt+i] += 2.000000 * qdot;
        wdot[13*npt+i] -= qdot;

        /*reaction 145: CH2 + HCCO <=> C2H3 + CO */
        phi_f = sc[10*npt+i]*sc[27*npt+i];
        k_f = k_f_s[144*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[14*npt+i]*sc[23*npt+i];
        Kc = Kc_s[144*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[10*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;
        wdot[23*npt+i] += qdot;
        wdot[27*npt+i] -= qdot;

        /*reaction 146: CH2(S) + N2 <=> CH2 + N2 */
        phi_f = sc[11*npt+i]*sc[30*npt+i];
        k_f = k_f_s[145*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[10*npt+i]*sc[30*npt+i];
        Kc = Kc_s[145*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[10*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;
        wdot[30*npt+i] -= qdot;
        wdot[30*npt+i] += qdot;

        /*reaction 147: CH2(S) + AR <=> CH2 + AR */
        phi_f = sc[11*npt+i]*sc[31*npt+i];
        k_f = k_f_s[146*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[10*npt+i]*sc[31*npt+i];
        Kc = Kc_s[146*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[10*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;
        wdot[31*npt+i] -= qdot;
        wdot[31*npt+i] += qdot;

        /*reaction 148: CH2(S) + O2 <=> H + OH + CO */
        phi_f = sc[3*npt+i]*sc[11*npt+i];
        k_f = k_f_s[147*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[4*npt+i]*sc[14*npt+i];
        Kc = Kc_s[147*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;

        /*reaction 149: CH2(S) + O2 <=> CO + H2O */
        phi_f = sc[3*npt+i]*sc[11*npt+i];
        k_f = k_f_s[148*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[14*npt+i];
        Kc = Kc_s[148*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;

        /*reaction 150: CH2(S) + H2 <=> CH3 + H */
        phi_f = sc[0*npt+i]*sc[11*npt+i];
        k_f = k_f_s[149*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[12*npt+i];
        Kc = Kc_s[149*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;
    }
}

void vcomp_wdot_151_177(int npt, amrex::Real *  wdot, amrex::Real *  mixture, amrex::Real *  sc,
		amrex::Real *  k_f_s, amrex::Real *  Kc_s,
		amrex::Real *  tc, amrex::Real *  invT, amrex::Real *  T)
{
    for (int i=0; i<npt; i++) {
        amrex::Real qdot, q_f, q_r, phi_f, phi_r, k_f, k_r, Kc;
        amrex::Real alpha;

        /*reaction 151: CH2(S) + H2O <=> CH2 + H2O */
        phi_f = sc[5*npt+i]*sc[11*npt+i];
        k_f = k_f_s[150*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[10*npt+i];
        Kc = Kc_s[150*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[10*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;

        /*reaction 152: CH2(S) + CH3 <=> H + C2H4 */
        phi_f = sc[11*npt+i]*sc[12*npt+i];
        k_f = k_f_s[151*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[24*npt+i];
        Kc = Kc_s[151*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;
        wdot[12*npt+i] -= qdot;
        wdot[24*npt+i] += qdot;

        /*reaction 153: CH2(S) + CH4 <=> 2.000000 CH3 */
        phi_f = sc[11*npt+i]*sc[13*npt+i];
        k_f = k_f_s[152*npt+i];
        q_f = phi_f * k_f;
        phi_r = pow(sc[12*npt+i], 2.000000);
        Kc = Kc_s[152*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[11*npt+i] -= qdot;
        wdot[12*npt+i] += 2.000000 * qdot;
        wdot[13*npt+i] -= qdot;

        /*reaction 154: CH2(S) + CO <=> CH2 + CO */
        phi_f = sc[11*npt+i]*sc[14*npt+i];
        k_f = k_f_s[153*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[10*npt+i]*sc[14*npt+i];
        Kc = Kc_s[153*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[10*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;
        wdot[14*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;

        /*reaction 155: CH2(S) + CO2 <=> CH2 + CO2 */
        phi_f = sc[11*npt+i]*sc[15*npt+i];
        k_f = k_f_s[154*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[10*npt+i]*sc[15*npt+i];
        Kc = Kc_s[154*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[10*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;
        wdot[15*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;

        /*reaction 156: CH2(S) + CO2 <=> CO + CH2O */
        phi_f = sc[11*npt+i]*sc[15*npt+i];
        k_f = k_f_s[155*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[14*npt+i]*sc[17*npt+i];
        Kc = Kc_s[155*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[11*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;
        wdot[15*npt+i] -= qdot;
        wdot[17*npt+i] += qdot;

        /*reaction 157: CH2(S) + C2H6 <=> CH3 + C2H5 */
        phi_f = sc[11*npt+i]*sc[26*npt+i];
        k_f = k_f_s[156*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[12*npt+i]*sc[25*npt+i];
        Kc = Kc_s[156*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[11*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;
        wdot[25*npt+i] += qdot;
        wdot[26*npt+i] -= qdot;

        /*reaction 158: CH3 + O2 <=> O + CH3O */
        phi_f = sc[3*npt+i]*sc[12*npt+i];
        k_f = k_f_s[157*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[19*npt+i];
        Kc = Kc_s[157*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[12*npt+i] -= qdot;
        wdot[19*npt+i] += qdot;

        /*reaction 159: CH3 + O2 <=> OH + CH2O */
        phi_f = sc[3*npt+i]*sc[12*npt+i];
        k_f = k_f_s[158*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[17*npt+i];
        Kc = Kc_s[158*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[12*npt+i] -= qdot;
        wdot[17*npt+i] += qdot;

        /*reaction 160: CH3 + H2O2 <=> HO2 + CH4 */
        phi_f = sc[7*npt+i]*sc[12*npt+i];
        k_f = k_f_s[159*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[13*npt+i];
        Kc = Kc_s[159*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;
        wdot[12*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;

        /*reaction 161: 2.000000 CH3 <=> H + C2H5 */
        phi_f = pow(sc[12*npt+i], 2.000000);
        k_f = k_f_s[160*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[25*npt+i];
        Kc = Kc_s[160*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[12*npt+i] -= 2.000000 * qdot;
        wdot[25*npt+i] += qdot;

        /*reaction 162: CH3 + HCO <=> CH4 + CO */
        phi_f = sc[12*npt+i]*sc[16*npt+i];
        k_f = k_f_s[161*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[13*npt+i]*sc[14*npt+i];
        Kc = Kc_s[161*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[12*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;
        wdot[14*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;

        /*reaction 163: CH3 + CH2O <=> HCO + CH4 */
        phi_f = sc[12*npt+i]*sc[17*npt+i];
        k_f = k_f_s[162*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[13*npt+i]*sc[16*npt+i];
        Kc = Kc_s[162*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[12*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;
        wdot[16*npt+i] += qdot;
        wdot[17*npt+i] -= qdot;

        /*reaction 164: CH3 + CH3OH <=> CH2OH + CH4 */
        phi_f = sc[12*npt+i]*sc[20*npt+i];
        k_f = k_f_s[163*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[13*npt+i]*sc[18*npt+i];
        Kc = Kc_s[163*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[12*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;
        wdot[18*npt+i] += qdot;
        wdot[20*npt+i] -= qdot;

        /*reaction 165: CH3 + CH3OH <=> CH3O + CH4 */
        phi_f = sc[12*npt+i]*sc[20*npt+i];
        k_f = k_f_s[164*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[13*npt+i]*sc[19*npt+i];
        Kc = Kc_s[164*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[12*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;
        wdot[19*npt+i] += qdot;
        wdot[20*npt+i] -= qdot;

        /*reaction 166: CH3 + C2H4 <=> C2H3 + CH4 */
        phi_f = sc[12*npt+i]*sc[24*npt+i];
        k_f = k_f_s[165*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[13*npt+i]*sc[23*npt+i];
        Kc = Kc_s[165*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[12*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;
        wdot[23*npt+i] += qdot;
        wdot[24*npt+i] -= qdot;

        /*reaction 167: CH3 + C2H6 <=> C2H5 + CH4 */
        phi_f = sc[12*npt+i]*sc[26*npt+i];
        k_f = k_f_s[166*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[13*npt+i]*sc[25*npt+i];
        Kc = Kc_s[166*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[12*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;
        wdot[25*npt+i] += qdot;
        wdot[26*npt+i] -= qdot;

        /*reaction 168: HCO + H2O <=> H + CO + H2O */
        phi_f = sc[5*npt+i]*sc[16*npt+i];
        k_f = k_f_s[167*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[5*npt+i]*sc[14*npt+i];
        Kc = Kc_s[167*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[14*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;

        /*reaction 169: HCO + O2 <=> HO2 + CO */
        phi_f = sc[3*npt+i]*sc[16*npt+i];
        k_f = k_f_s[168*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[14*npt+i];
        Kc = Kc_s[168*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[14*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;

        /*reaction 170: CH2OH + O2 <=> HO2 + CH2O */
        phi_f = sc[3*npt+i]*sc[18*npt+i];
        k_f = k_f_s[169*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[17*npt+i];
        Kc = Kc_s[169*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[17*npt+i] += qdot;
        wdot[18*npt+i] -= qdot;

        /*reaction 171: CH3O + O2 <=> HO2 + CH2O */
        phi_f = sc[3*npt+i]*sc[19*npt+i];
        k_f = k_f_s[170*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[17*npt+i];
        Kc = Kc_s[170*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[17*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 172: C2H + O2 <=> HCO + CO */
        phi_f = sc[3*npt+i]*sc[21*npt+i];
        k_f = k_f_s[171*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[14*npt+i]*sc[16*npt+i];
        Kc = Kc_s[171*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;
        wdot[16*npt+i] += qdot;
        wdot[21*npt+i] -= qdot;

        /*reaction 173: C2H + H2 <=> H + C2H2 */
        phi_f = sc[0*npt+i]*sc[21*npt+i];
        k_f = k_f_s[172*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[22*npt+i];
        Kc = Kc_s[172*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[21*npt+i] -= qdot;
        wdot[22*npt+i] += qdot;

        /*reaction 174: C2H3 + O2 <=> HCO + CH2O */
        phi_f = sc[3*npt+i]*sc[23*npt+i];
        k_f = k_f_s[173*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[16*npt+i]*sc[17*npt+i];
        Kc = Kc_s[173*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[16*npt+i] += qdot;
        wdot[17*npt+i] += qdot;
        wdot[23*npt+i] -= qdot;

        /*reaction 175: C2H5 + O2 <=> HO2 + C2H4 */
        phi_f = sc[3*npt+i]*sc[25*npt+i];
        k_f = k_f_s[174*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[24*npt+i];
        Kc = Kc_s[174*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[24*npt+i] += qdot;
        wdot[25*npt+i] -= qdot;

        /*reaction 176: HCCO + O2 <=> OH + 2.000000 CO */
        phi_f = sc[3*npt+i]*sc[27*npt+i];
        k_f = k_f_s[175*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*pow(sc[14*npt+i], 2.000000);
        Kc = Kc_s[175*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[14*npt+i] += 2.000000 * qdot;
        wdot[27*npt+i] -= qdot;

        /*reaction 177: 2.000000 HCCO <=> 2.000000 CO + C2H2 */
        phi_f = pow(sc[27*npt+i], 2.000000);
        k_f = k_f_s[176*npt+i];
        q_f = phi_f * k_f;
        phi_r = pow(sc[14*npt+i], 2.000000)*sc[22*npt+i];
        Kc = Kc_s[176*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[14*npt+i] += 2.000000 * qdot;
        wdot[22*npt+i] += qdot;
        wdot[27*npt+i] -= 2.000000 * qdot;
    }
}


static amrex::Real T_save = -1;
#ifdef _OPENMP
#pragma omp threadprivate(T_save)
#endif

static amrex::Real k_f_save[177];
#ifdef _OPENMP
#pragma omp threadprivate(k_f_save)
#endif

static amrex::Real Kc_save[177];
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

    amrex::Real qdot, q_f[177], q_r[177];
    amrex::Real sc_qss[1];
    comp_qfqr(q_f, q_r, sc, sc_qss, tc, invT);

    for (int i = 0; i < 32; ++i) {
        wdot[i] = 0.0;
    }

    qdot = q_f[0]-q_r[0];
    wdot[1] -= qdot;
    wdot[10] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[1]-q_r[1];
    wdot[1] -= qdot;
    wdot[12] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[2]-q_r[2];
    wdot[1] -= qdot;
    wdot[16] -= qdot;
    wdot[17] += qdot;

    qdot = q_f[3]-q_r[3];
    wdot[1] -= qdot;
    wdot[17] -= qdot;
    wdot[18] += qdot;

    qdot = q_f[4]-q_r[4];
    wdot[1] -= qdot;
    wdot[17] -= qdot;
    wdot[19] += qdot;

    qdot = q_f[5]-q_r[5];
    wdot[1] -= qdot;
    wdot[18] -= qdot;
    wdot[20] += qdot;

    qdot = q_f[6]-q_r[6];
    wdot[1] -= qdot;
    wdot[19] -= qdot;
    wdot[20] += qdot;

    qdot = q_f[7]-q_r[7];
    wdot[1] -= qdot;
    wdot[21] -= qdot;
    wdot[22] += qdot;

    qdot = q_f[8]-q_r[8];
    wdot[1] -= qdot;
    wdot[22] -= qdot;
    wdot[23] += qdot;

    qdot = q_f[9]-q_r[9];
    wdot[1] -= qdot;
    wdot[23] -= qdot;
    wdot[24] += qdot;

    qdot = q_f[10]-q_r[10];
    wdot[1] -= qdot;
    wdot[24] -= qdot;
    wdot[25] += qdot;

    qdot = q_f[11]-q_r[11];
    wdot[1] -= qdot;
    wdot[25] -= qdot;
    wdot[26] += qdot;

    qdot = q_f[12]-q_r[12];
    wdot[0] -= qdot;
    wdot[14] -= qdot;
    wdot[17] += qdot;

    qdot = q_f[13]-q_r[13];
    wdot[4] -= 2.000000 * qdot;
    wdot[7] += qdot;

    qdot = q_f[14]-q_r[14];
    wdot[4] -= qdot;
    wdot[12] -= qdot;
    wdot[20] += qdot;

    qdot = q_f[15]-q_r[15];
    wdot[9] -= qdot;
    wdot[14] -= qdot;
    wdot[27] += qdot;

    qdot = q_f[16]-q_r[16];
    wdot[10] -= qdot;
    wdot[14] -= qdot;
    wdot[28] += qdot;

    qdot = q_f[17]-q_r[17];
    wdot[5] -= qdot;
    wdot[11] -= qdot;
    wdot[20] += qdot;

    qdot = q_f[18]-q_r[18];
    wdot[12] -= 2.000000 * qdot;
    wdot[26] += qdot;

    qdot = q_f[19]-q_r[19];
    wdot[0] += qdot;
    wdot[22] += qdot;
    wdot[24] -= qdot;

    qdot = q_f[20]-q_r[20];
    wdot[2] -= 2.000000 * qdot;
    wdot[3] += qdot;

    qdot = q_f[21]-q_r[21];
    wdot[1] -= qdot;
    wdot[2] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[22]-q_r[22];
    wdot[2] -= qdot;
    wdot[14] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[23]-q_r[23];
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[6] += qdot;

    qdot = q_f[24]-q_r[24];
    wdot[0] += qdot;
    wdot[1] -= 2.000000 * qdot;

    qdot = q_f[25]-q_r[25];
    wdot[1] -= qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[26]-q_r[26];
    wdot[1] += qdot;
    wdot[14] += qdot;
    wdot[16] -= qdot;

    qdot = q_f[27]-q_r[27];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[28]-q_r[28];
    wdot[2] -= qdot;
    wdot[3] += qdot;
    wdot[4] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[29]-q_r[29];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[30]-q_r[30];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[9] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[31]-q_r[31];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[10] -= qdot;
    wdot[16] += qdot;

    qdot = q_f[32]-q_r[32];
    wdot[0] += qdot;
    wdot[2] -= qdot;
    wdot[11] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[33]-q_r[33];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[11] -= qdot;
    wdot[16] += qdot;

    qdot = q_f[34]-q_r[34];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[12] -= qdot;
    wdot[17] += qdot;

    qdot = q_f[35]-q_r[35];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[12] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[36]-q_r[36];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[14] += qdot;
    wdot[16] -= qdot;

    qdot = q_f[37]-q_r[37];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[15] += qdot;
    wdot[16] -= qdot;

    qdot = q_f[38]-q_r[38];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[16] += qdot;
    wdot[17] -= qdot;

    qdot = q_f[39]-q_r[39];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[17] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[40]-q_r[40];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[17] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[41]-q_r[41];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[18] += qdot;
    wdot[20] -= qdot;

    qdot = q_f[42]-q_r[42];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[19] += qdot;
    wdot[20] -= qdot;

    qdot = q_f[43]-q_r[43];
    wdot[2] -= qdot;
    wdot[9] += qdot;
    wdot[14] += qdot;
    wdot[21] -= qdot;

    qdot = q_f[44]-q_r[44];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[22] -= qdot;
    wdot[27] += qdot;

    qdot = q_f[45]-q_r[45];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[21] += qdot;
    wdot[22] -= qdot;

    qdot = q_f[46]-q_r[46];
    wdot[2] -= qdot;
    wdot[10] += qdot;
    wdot[14] += qdot;
    wdot[22] -= qdot;

    qdot = q_f[47]-q_r[47];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[23] -= qdot;
    wdot[28] += qdot;

    qdot = q_f[48]-q_r[48];
    wdot[2] -= qdot;
    wdot[12] += qdot;
    wdot[16] += qdot;
    wdot[24] -= qdot;

    qdot = q_f[49]-q_r[49];
    wdot[2] -= qdot;
    wdot[12] += qdot;
    wdot[17] += qdot;
    wdot[25] -= qdot;

    qdot = q_f[50]-q_r[50];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[25] += qdot;
    wdot[26] -= qdot;

    qdot = q_f[51]-q_r[51];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[14] += 2.000000 * qdot;
    wdot[27] -= qdot;

    qdot = q_f[52]-q_r[52];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[27] += qdot;
    wdot[28] -= qdot;

    qdot = q_f[53]-q_r[53];
    wdot[2] -= qdot;
    wdot[10] += qdot;
    wdot[15] += qdot;
    wdot[28] -= qdot;

    qdot = q_f[54]-q_r[54];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[14] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[55]-q_r[55];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[16] += qdot;
    wdot[17] -= qdot;

    qdot = q_f[56]-q_r[56];
    wdot[1] -= qdot;
    wdot[3] -= 2.000000 * qdot;
    wdot[3] += qdot;
    wdot[6] += qdot;

    qdot = q_f[57]-q_r[57];
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[5] -= qdot;
    wdot[5] += qdot;
    wdot[6] += qdot;

    qdot = q_f[58]-q_r[58];
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[30] -= qdot;
    wdot[30] += qdot;

    qdot = q_f[59]-q_r[59];
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[31] -= qdot;
    wdot[31] += qdot;

    qdot = q_f[60]-q_r[60];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[61]-q_r[61];
    wdot[0] -= qdot;
    wdot[0] += 2.000000 * qdot;
    wdot[1] -= 2.000000 * qdot;

    qdot = q_f[62]-q_r[62];
    wdot[0] += qdot;
    wdot[1] -= 2.000000 * qdot;
    wdot[5] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[63]-q_r[63];
    wdot[0] += qdot;
    wdot[1] -= 2.000000 * qdot;
    wdot[15] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[64]-q_r[64];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[5] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[65]-q_r[65];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[3] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[66]-q_r[66];
    wdot[1] -= qdot;
    wdot[4] += 2.000000 * qdot;
    wdot[6] -= qdot;

    qdot = q_f[67]-q_r[67];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[68]-q_r[68];
    wdot[1] -= qdot;
    wdot[4] += qdot;
    wdot[5] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[69]-q_r[69];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[8] += qdot;
    wdot[9] -= qdot;

    qdot = q_f[70]-q_r[70];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[9] += qdot;
    wdot[11] -= qdot;

    qdot = q_f[71]-q_r[71];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[12] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[72]-q_r[72];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[14] += qdot;
    wdot[16] -= qdot;

    qdot = q_f[73]-q_r[73];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[16] += qdot;
    wdot[17] -= qdot;

    qdot = q_f[74]-q_r[74];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[17] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[75]-q_r[75];
    wdot[1] -= qdot;
    wdot[4] += qdot;
    wdot[12] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[76]-q_r[76];
    wdot[1] -= qdot;
    wdot[5] += qdot;
    wdot[11] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[77]-q_r[77];
    wdot[1] -= qdot;
    wdot[1] += qdot;
    wdot[18] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[78]-q_r[78];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[17] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[79]-q_r[79];
    wdot[1] -= qdot;
    wdot[4] += qdot;
    wdot[12] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[80]-q_r[80];
    wdot[1] -= qdot;
    wdot[5] += qdot;
    wdot[11] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[81]-q_r[81];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[18] += qdot;
    wdot[20] -= qdot;

    qdot = q_f[82]-q_r[82];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[19] += qdot;
    wdot[20] -= qdot;

    qdot = q_f[83]-q_r[83];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[22] += qdot;
    wdot[23] -= qdot;

    qdot = q_f[84]-q_r[84];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[23] += qdot;
    wdot[24] -= qdot;

    qdot = q_f[85]-q_r[85];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[24] += qdot;
    wdot[25] -= qdot;

    qdot = q_f[86]-q_r[86];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[25] += qdot;
    wdot[26] -= qdot;

    qdot = q_f[87]-q_r[87];
    wdot[1] -= qdot;
    wdot[11] += qdot;
    wdot[14] += qdot;
    wdot[27] -= qdot;

    qdot = q_f[88]-q_r[88];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[27] += qdot;
    wdot[28] -= qdot;

    qdot = q_f[89]-q_r[89];
    wdot[1] -= qdot;
    wdot[12] += qdot;
    wdot[14] += qdot;
    wdot[28] -= qdot;

    qdot = q_f[90]-q_r[90];
    wdot[1] -= qdot;
    wdot[1] += qdot;
    wdot[28] += qdot;
    wdot[29] -= qdot;

    qdot = q_f[91]-q_r[91];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[92]-q_r[92];
    wdot[2] += qdot;
    wdot[4] -= 2.000000 * qdot;
    wdot[5] += qdot;

    qdot = q_f[93]-q_r[93];
    wdot[3] += qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[94]-q_r[94];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[95]-q_r[95];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[96]-q_r[96];
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[8] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[97]-q_r[97];
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[9] -= qdot;
    wdot[16] += qdot;

    qdot = q_f[98]-q_r[98];
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[10] -= qdot;
    wdot[17] += qdot;

    qdot = q_f[99]-q_r[99];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[9] += qdot;
    wdot[10] -= qdot;

    qdot = q_f[100]-q_r[100];
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[11] -= qdot;
    wdot[17] += qdot;

    qdot = q_f[101]-q_r[101];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[10] += qdot;
    wdot[12] -= qdot;

    qdot = q_f[102]-q_r[102];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[11] += qdot;
    wdot[12] -= qdot;

    qdot = q_f[103]-q_r[103];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[12] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[104]-q_r[104];
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[14] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[105]-q_r[105];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[14] += qdot;
    wdot[16] -= qdot;

    qdot = q_f[106]-q_r[106];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[16] += qdot;
    wdot[17] -= qdot;

    qdot = q_f[107]-q_r[107];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[17] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[108]-q_r[108];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[17] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[109]-q_r[109];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[18] += qdot;
    wdot[20] -= qdot;

    qdot = q_f[110]-q_r[110];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[19] += qdot;
    wdot[20] -= qdot;

    qdot = q_f[111]-q_r[111];
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[21] -= qdot;
    wdot[27] += qdot;

    qdot = q_f[112]-q_r[112];
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[22] -= qdot;
    wdot[28] += qdot;

    qdot = q_f[113]-q_r[113];
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[22] -= qdot;
    wdot[29] += qdot;

    qdot = q_f[114]-q_r[114];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[21] += qdot;
    wdot[22] -= qdot;

    qdot = q_f[115]-q_r[115];
    wdot[4] -= qdot;
    wdot[12] += qdot;
    wdot[14] += qdot;
    wdot[22] -= qdot;

    qdot = q_f[116]-q_r[116];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[22] += qdot;
    wdot[23] -= qdot;

    qdot = q_f[117]-q_r[117];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[23] += qdot;
    wdot[24] -= qdot;

    qdot = q_f[118]-q_r[118];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[25] += qdot;
    wdot[26] -= qdot;

    qdot = q_f[119]-q_r[119];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[27] += qdot;
    wdot[28] -= qdot;

    qdot = q_f[120]-q_r[120];
    wdot[3] += qdot;
    wdot[6] -= 2.000000 * qdot;
    wdot[7] += qdot;

    qdot = q_f[121]-q_r[121];
    wdot[3] += qdot;
    wdot[6] -= 2.000000 * qdot;
    wdot[7] += qdot;

    qdot = q_f[122]-q_r[122];
    wdot[4] += qdot;
    wdot[6] -= qdot;
    wdot[10] -= qdot;
    wdot[17] += qdot;

    qdot = q_f[123]-q_r[123];
    wdot[3] += qdot;
    wdot[6] -= qdot;
    wdot[12] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[124]-q_r[124];
    wdot[4] += qdot;
    wdot[6] -= qdot;
    wdot[12] -= qdot;
    wdot[19] += qdot;

    qdot = q_f[125]-q_r[125];
    wdot[4] += qdot;
    wdot[6] -= qdot;
    wdot[14] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[126]-q_r[126];
    wdot[6] -= qdot;
    wdot[7] += qdot;
    wdot[16] += qdot;
    wdot[17] -= qdot;

    qdot = q_f[127]-q_r[127];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[8] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[128]-q_r[128];
    wdot[1] += qdot;
    wdot[8] -= qdot;
    wdot[10] -= qdot;
    wdot[21] += qdot;

    qdot = q_f[129]-q_r[129];
    wdot[1] += qdot;
    wdot[8] -= qdot;
    wdot[12] -= qdot;
    wdot[22] += qdot;

    qdot = q_f[130]-q_r[130];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[9] -= qdot;
    wdot[16] += qdot;

    qdot = q_f[131]-q_r[131];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[9] -= qdot;
    wdot[10] += qdot;

    qdot = q_f[132]-q_r[132];
    wdot[1] += qdot;
    wdot[5] -= qdot;
    wdot[9] -= qdot;
    wdot[17] += qdot;

    qdot = q_f[133]-q_r[133];
    wdot[1] += qdot;
    wdot[9] -= qdot;
    wdot[10] -= qdot;
    wdot[22] += qdot;

    qdot = q_f[134]-q_r[134];
    wdot[1] += qdot;
    wdot[9] -= qdot;
    wdot[12] -= qdot;
    wdot[23] += qdot;

    qdot = q_f[135]-q_r[135];
    wdot[1] += qdot;
    wdot[9] -= qdot;
    wdot[13] -= qdot;
    wdot[24] += qdot;

    qdot = q_f[136]-q_r[136];
    wdot[9] -= qdot;
    wdot[14] += qdot;
    wdot[15] -= qdot;
    wdot[16] += qdot;

    qdot = q_f[137]-q_r[137];
    wdot[1] += qdot;
    wdot[9] -= qdot;
    wdot[17] -= qdot;
    wdot[28] += qdot;

    qdot = q_f[138]-q_r[138];
    wdot[9] -= qdot;
    wdot[14] += qdot;
    wdot[22] += qdot;
    wdot[27] -= qdot;

    qdot = q_f[139]-q_r[139];
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[10] -= qdot;
    wdot[16] += qdot;

    qdot = q_f[140]-q_r[140];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[10] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[141]-q_r[141];
    wdot[0] += qdot;
    wdot[10] -= 2.000000 * qdot;
    wdot[22] += qdot;

    qdot = q_f[142]-q_r[142];
    wdot[1] += qdot;
    wdot[10] -= qdot;
    wdot[12] -= qdot;
    wdot[24] += qdot;

    qdot = q_f[143]-q_r[143];
    wdot[10] -= qdot;
    wdot[12] += 2.000000 * qdot;
    wdot[13] -= qdot;

    qdot = q_f[144]-q_r[144];
    wdot[10] -= qdot;
    wdot[14] += qdot;
    wdot[23] += qdot;
    wdot[27] -= qdot;

    qdot = q_f[145]-q_r[145];
    wdot[10] += qdot;
    wdot[11] -= qdot;
    wdot[30] -= qdot;
    wdot[30] += qdot;

    qdot = q_f[146]-q_r[146];
    wdot[10] += qdot;
    wdot[11] -= qdot;
    wdot[31] -= qdot;
    wdot[31] += qdot;

    qdot = q_f[147]-q_r[147];
    wdot[1] += qdot;
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[11] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[148]-q_r[148];
    wdot[3] -= qdot;
    wdot[5] += qdot;
    wdot[11] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[149]-q_r[149];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[11] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[150]-q_r[150];
    wdot[5] -= qdot;
    wdot[5] += qdot;
    wdot[10] += qdot;
    wdot[11] -= qdot;

    qdot = q_f[151]-q_r[151];
    wdot[1] += qdot;
    wdot[11] -= qdot;
    wdot[12] -= qdot;
    wdot[24] += qdot;

    qdot = q_f[152]-q_r[152];
    wdot[11] -= qdot;
    wdot[12] += 2.000000 * qdot;
    wdot[13] -= qdot;

    qdot = q_f[153]-q_r[153];
    wdot[10] += qdot;
    wdot[11] -= qdot;
    wdot[14] -= qdot;
    wdot[14] += qdot;

    qdot = q_f[154]-q_r[154];
    wdot[10] += qdot;
    wdot[11] -= qdot;
    wdot[15] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[155]-q_r[155];
    wdot[11] -= qdot;
    wdot[14] += qdot;
    wdot[15] -= qdot;
    wdot[17] += qdot;

    qdot = q_f[156]-q_r[156];
    wdot[11] -= qdot;
    wdot[12] += qdot;
    wdot[25] += qdot;
    wdot[26] -= qdot;

    qdot = q_f[157]-q_r[157];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[12] -= qdot;
    wdot[19] += qdot;

    qdot = q_f[158]-q_r[158];
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[12] -= qdot;
    wdot[17] += qdot;

    qdot = q_f[159]-q_r[159];
    wdot[6] += qdot;
    wdot[7] -= qdot;
    wdot[12] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[160]-q_r[160];
    wdot[1] += qdot;
    wdot[12] -= 2.000000 * qdot;
    wdot[25] += qdot;

    qdot = q_f[161]-q_r[161];
    wdot[12] -= qdot;
    wdot[13] += qdot;
    wdot[14] += qdot;
    wdot[16] -= qdot;

    qdot = q_f[162]-q_r[162];
    wdot[12] -= qdot;
    wdot[13] += qdot;
    wdot[16] += qdot;
    wdot[17] -= qdot;

    qdot = q_f[163]-q_r[163];
    wdot[12] -= qdot;
    wdot[13] += qdot;
    wdot[18] += qdot;
    wdot[20] -= qdot;

    qdot = q_f[164]-q_r[164];
    wdot[12] -= qdot;
    wdot[13] += qdot;
    wdot[19] += qdot;
    wdot[20] -= qdot;

    qdot = q_f[165]-q_r[165];
    wdot[12] -= qdot;
    wdot[13] += qdot;
    wdot[23] += qdot;
    wdot[24] -= qdot;

    qdot = q_f[166]-q_r[166];
    wdot[12] -= qdot;
    wdot[13] += qdot;
    wdot[25] += qdot;
    wdot[26] -= qdot;

    qdot = q_f[167]-q_r[167];
    wdot[1] += qdot;
    wdot[5] -= qdot;
    wdot[5] += qdot;
    wdot[14] += qdot;
    wdot[16] -= qdot;

    qdot = q_f[168]-q_r[168];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[14] += qdot;
    wdot[16] -= qdot;

    qdot = q_f[169]-q_r[169];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[17] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[170]-q_r[170];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[17] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[171]-q_r[171];
    wdot[3] -= qdot;
    wdot[14] += qdot;
    wdot[16] += qdot;
    wdot[21] -= qdot;

    qdot = q_f[172]-q_r[172];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[21] -= qdot;
    wdot[22] += qdot;

    qdot = q_f[173]-q_r[173];
    wdot[3] -= qdot;
    wdot[16] += qdot;
    wdot[17] += qdot;
    wdot[23] -= qdot;

    qdot = q_f[174]-q_r[174];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[24] += qdot;
    wdot[25] -= qdot;

    qdot = q_f[175]-q_r[175];
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[14] += 2.000000 * qdot;
    wdot[27] -= qdot;

    qdot = q_f[176]-q_r[176];
    wdot[14] += 2.000000 * qdot;
    wdot[22] += qdot;
    wdot[27] -= 2.000000 * qdot;

    return;
}

void comp_k_f(amrex::Real *  tc, amrex::Real invT, amrex::Real *  k_f)
{
    for (int i=0; i<177; ++i) {
        k_f[i] = prefactor_units[i] * fwd_A[i]
                    * exp(fwd_beta[i] * tc[0] - activation_units[i] * fwd_Ea[i] * invT);
    };
    return;
}

void comp_Kc(amrex::Real *  tc, amrex::Real invT, amrex::Real *  Kc)
{
    /*compute the Gibbs free energy */
    amrex::Real g_RT[32];
    gibbs(g_RT, tc);

    Kc[0] = g_RT[1] + g_RT[10] - g_RT[12];
    Kc[1] = g_RT[1] + g_RT[12] - g_RT[13];
    Kc[2] = g_RT[1] + g_RT[16] - g_RT[17];
    Kc[3] = g_RT[1] + g_RT[17] - g_RT[18];
    Kc[4] = g_RT[1] + g_RT[17] - g_RT[19];
    Kc[5] = g_RT[1] + g_RT[18] - g_RT[20];
    Kc[6] = g_RT[1] + g_RT[19] - g_RT[20];
    Kc[7] = g_RT[1] + g_RT[21] - g_RT[22];
    Kc[8] = g_RT[1] + g_RT[22] - g_RT[23];
    Kc[9] = g_RT[1] + g_RT[23] - g_RT[24];
    Kc[10] = g_RT[1] + g_RT[24] - g_RT[25];
    Kc[11] = g_RT[1] + g_RT[25] - g_RT[26];
    Kc[12] = g_RT[0] + g_RT[14] - g_RT[17];
    Kc[13] = 2.000000*g_RT[4] - g_RT[7];
    Kc[14] = g_RT[4] + g_RT[12] - g_RT[20];
    Kc[15] = g_RT[9] + g_RT[14] - g_RT[27];
    Kc[16] = g_RT[10] + g_RT[14] - g_RT[28];
    Kc[17] = g_RT[5] + g_RT[11] - g_RT[20];
    Kc[18] = 2.000000*g_RT[12] - g_RT[26];
    Kc[19] = -g_RT[0] - g_RT[22] + g_RT[24];
    Kc[20] = 2.000000*g_RT[2] - g_RT[3];
    Kc[21] = g_RT[1] + g_RT[2] - g_RT[4];
    Kc[22] = g_RT[2] + g_RT[14] - g_RT[15];
    Kc[23] = g_RT[1] + g_RT[3] - g_RT[6];
    Kc[24] = -g_RT[0] + 2.000000*g_RT[1];
    Kc[25] = g_RT[1] + g_RT[4] - g_RT[5];
    Kc[26] = -g_RT[1] - g_RT[14] + g_RT[16];
    Kc[27] = g_RT[0] - g_RT[1] + g_RT[2] - g_RT[4];
    Kc[28] = g_RT[2] - g_RT[3] - g_RT[4] + g_RT[6];
    Kc[29] = g_RT[2] - g_RT[4] - g_RT[6] + g_RT[7];
    Kc[30] = -g_RT[1] + g_RT[2] + g_RT[9] - g_RT[14];
    Kc[31] = -g_RT[1] + g_RT[2] + g_RT[10] - g_RT[16];
    Kc[32] = -g_RT[0] + g_RT[2] + g_RT[11] - g_RT[14];
    Kc[33] = -g_RT[1] + g_RT[2] + g_RT[11] - g_RT[16];
    Kc[34] = -g_RT[1] + g_RT[2] + g_RT[12] - g_RT[17];
    Kc[35] = g_RT[2] - g_RT[4] - g_RT[12] + g_RT[13];
    Kc[36] = g_RT[2] - g_RT[4] - g_RT[14] + g_RT[16];
    Kc[37] = -g_RT[1] + g_RT[2] - g_RT[15] + g_RT[16];
    Kc[38] = g_RT[2] - g_RT[4] - g_RT[16] + g_RT[17];
    Kc[39] = g_RT[2] - g_RT[4] - g_RT[17] + g_RT[18];
    Kc[40] = g_RT[2] - g_RT[4] - g_RT[17] + g_RT[19];
    Kc[41] = g_RT[2] - g_RT[4] - g_RT[18] + g_RT[20];
    Kc[42] = g_RT[2] - g_RT[4] - g_RT[19] + g_RT[20];
    Kc[43] = g_RT[2] - g_RT[9] - g_RT[14] + g_RT[21];
    Kc[44] = -g_RT[1] + g_RT[2] + g_RT[22] - g_RT[27];
    Kc[45] = g_RT[2] - g_RT[4] - g_RT[21] + g_RT[22];
    Kc[46] = g_RT[2] - g_RT[10] - g_RT[14] + g_RT[22];
    Kc[47] = -g_RT[1] + g_RT[2] + g_RT[23] - g_RT[28];
    Kc[48] = g_RT[2] - g_RT[12] - g_RT[16] + g_RT[24];
    Kc[49] = g_RT[2] - g_RT[12] - g_RT[17] + g_RT[25];
    Kc[50] = g_RT[2] - g_RT[4] - g_RT[25] + g_RT[26];
    Kc[51] = -g_RT[1] + g_RT[2] - 2.000000*g_RT[14] + g_RT[27];
    Kc[52] = g_RT[2] - g_RT[4] - g_RT[27] + g_RT[28];
    Kc[53] = g_RT[2] - g_RT[10] - g_RT[15] + g_RT[28];
    Kc[54] = -g_RT[2] + g_RT[3] + g_RT[14] - g_RT[15];
    Kc[55] = g_RT[3] - g_RT[6] - g_RT[16] + g_RT[17];
    Kc[56] = g_RT[1] + 2.000000*g_RT[3] - g_RT[3] - g_RT[6];
    Kc[57] = g_RT[1] + g_RT[3] + g_RT[5] - g_RT[5] - g_RT[6];
    Kc[58] = g_RT[1] + g_RT[3] - g_RT[6] + g_RT[30] - g_RT[30];
    Kc[59] = g_RT[1] + g_RT[3] - g_RT[6] + g_RT[31] - g_RT[31];
    Kc[60] = g_RT[1] - g_RT[2] + g_RT[3] - g_RT[4];
    Kc[61] = g_RT[0] - 2.000000*g_RT[0] + 2.000000*g_RT[1];
    Kc[62] = -g_RT[0] + 2.000000*g_RT[1] + g_RT[5] - g_RT[5];
    Kc[63] = -g_RT[0] + 2.000000*g_RT[1] + g_RT[15] - g_RT[15];
    Kc[64] = g_RT[1] - g_RT[2] - g_RT[5] + g_RT[6];
    Kc[65] = -g_RT[0] + g_RT[1] - g_RT[3] + g_RT[6];
    Kc[66] = g_RT[1] - 2.000000*g_RT[4] + g_RT[6];
    Kc[67] = -g_RT[0] + g_RT[1] - g_RT[6] + g_RT[7];
    Kc[68] = g_RT[1] - g_RT[4] - g_RT[5] + g_RT[7];
    Kc[69] = -g_RT[0] + g_RT[1] - g_RT[8] + g_RT[9];
    Kc[70] = -g_RT[0] + g_RT[1] - g_RT[9] + g_RT[11];
    Kc[71] = -g_RT[0] + g_RT[1] - g_RT[12] + g_RT[13];
    Kc[72] = -g_RT[0] + g_RT[1] - g_RT[14] + g_RT[16];
    Kc[73] = -g_RT[0] + g_RT[1] - g_RT[16] + g_RT[17];
    Kc[74] = -g_RT[0] + g_RT[1] - g_RT[17] + g_RT[18];
    Kc[75] = g_RT[1] - g_RT[4] - g_RT[12] + g_RT[18];
    Kc[76] = g_RT[1] - g_RT[5] - g_RT[11] + g_RT[18];
    Kc[77] = g_RT[1] - g_RT[1] - g_RT[18] + g_RT[19];
    Kc[78] = -g_RT[0] + g_RT[1] - g_RT[17] + g_RT[19];
    Kc[79] = g_RT[1] - g_RT[4] - g_RT[12] + g_RT[19];
    Kc[80] = g_RT[1] - g_RT[5] - g_RT[11] + g_RT[19];
    Kc[81] = -g_RT[0] + g_RT[1] - g_RT[18] + g_RT[20];
    Kc[82] = -g_RT[0] + g_RT[1] - g_RT[19] + g_RT[20];
    Kc[83] = -g_RT[0] + g_RT[1] - g_RT[22] + g_RT[23];
    Kc[84] = -g_RT[0] + g_RT[1] - g_RT[23] + g_RT[24];
    Kc[85] = -g_RT[0] + g_RT[1] - g_RT[24] + g_RT[25];
    Kc[86] = -g_RT[0] + g_RT[1] - g_RT[25] + g_RT[26];
    Kc[87] = g_RT[1] - g_RT[11] - g_RT[14] + g_RT[27];
    Kc[88] = -g_RT[0] + g_RT[1] - g_RT[27] + g_RT[28];
    Kc[89] = g_RT[1] - g_RT[12] - g_RT[14] + g_RT[28];
    Kc[90] = g_RT[1] - g_RT[1] - g_RT[28] + g_RT[29];
    Kc[91] = g_RT[0] - g_RT[1] + g_RT[4] - g_RT[5];
    Kc[92] = -g_RT[2] + 2.000000*g_RT[4] - g_RT[5];
    Kc[93] = -g_RT[3] + g_RT[4] - g_RT[5] + g_RT[6];
    Kc[94] = g_RT[4] - g_RT[5] - g_RT[6] + g_RT[7];
    Kc[95] = g_RT[4] - g_RT[5] - g_RT[6] + g_RT[7];
    Kc[96] = -g_RT[1] + g_RT[4] + g_RT[8] - g_RT[14];
    Kc[97] = -g_RT[1] + g_RT[4] + g_RT[9] - g_RT[16];
    Kc[98] = -g_RT[1] + g_RT[4] + g_RT[10] - g_RT[17];
    Kc[99] = g_RT[4] - g_RT[5] - g_RT[9] + g_RT[10];
    Kc[100] = -g_RT[1] + g_RT[4] + g_RT[11] - g_RT[17];
    Kc[101] = g_RT[4] - g_RT[5] - g_RT[10] + g_RT[12];
    Kc[102] = g_RT[4] - g_RT[5] - g_RT[11] + g_RT[12];
    Kc[103] = g_RT[4] - g_RT[5] - g_RT[12] + g_RT[13];
    Kc[104] = -g_RT[1] + g_RT[4] + g_RT[14] - g_RT[15];
    Kc[105] = g_RT[4] - g_RT[5] - g_RT[14] + g_RT[16];
    Kc[106] = g_RT[4] - g_RT[5] - g_RT[16] + g_RT[17];
    Kc[107] = g_RT[4] - g_RT[5] - g_RT[17] + g_RT[18];
    Kc[108] = g_RT[4] - g_RT[5] - g_RT[17] + g_RT[19];
    Kc[109] = g_RT[4] - g_RT[5] - g_RT[18] + g_RT[20];
    Kc[110] = g_RT[4] - g_RT[5] - g_RT[19] + g_RT[20];
    Kc[111] = -g_RT[1] + g_RT[4] + g_RT[21] - g_RT[27];
    Kc[112] = -g_RT[1] + g_RT[4] + g_RT[22] - g_RT[28];
    Kc[113] = -g_RT[1] + g_RT[4] + g_RT[22] - g_RT[29];
    Kc[114] = g_RT[4] - g_RT[5] - g_RT[21] + g_RT[22];
    Kc[115] = g_RT[4] - g_RT[12] - g_RT[14] + g_RT[22];
    Kc[116] = g_RT[4] - g_RT[5] - g_RT[22] + g_RT[23];
    Kc[117] = g_RT[4] - g_RT[5] - g_RT[23] + g_RT[24];
    Kc[118] = g_RT[4] - g_RT[5] - g_RT[25] + g_RT[26];
    Kc[119] = g_RT[4] - g_RT[5] - g_RT[27] + g_RT[28];
    Kc[120] = -g_RT[3] + 2.000000*g_RT[6] - g_RT[7];
    Kc[121] = -g_RT[3] + 2.000000*g_RT[6] - g_RT[7];
    Kc[122] = -g_RT[4] + g_RT[6] + g_RT[10] - g_RT[17];
    Kc[123] = -g_RT[3] + g_RT[6] + g_RT[12] - g_RT[13];
    Kc[124] = -g_RT[4] + g_RT[6] + g_RT[12] - g_RT[19];
    Kc[125] = -g_RT[4] + g_RT[6] + g_RT[14] - g_RT[15];
    Kc[126] = g_RT[6] - g_RT[7] - g_RT[16] + g_RT[17];
    Kc[127] = -g_RT[2] + g_RT[3] + g_RT[8] - g_RT[14];
    Kc[128] = -g_RT[1] + g_RT[8] + g_RT[10] - g_RT[21];
    Kc[129] = -g_RT[1] + g_RT[8] + g_RT[12] - g_RT[22];
    Kc[130] = -g_RT[2] + g_RT[3] + g_RT[9] - g_RT[16];
    Kc[131] = g_RT[0] - g_RT[1] + g_RT[9] - g_RT[10];
    Kc[132] = -g_RT[1] + g_RT[5] + g_RT[9] - g_RT[17];
    Kc[133] = -g_RT[1] + g_RT[9] + g_RT[10] - g_RT[22];
    Kc[134] = -g_RT[1] + g_RT[9] + g_RT[12] - g_RT[23];
    Kc[135] = -g_RT[1] + g_RT[9] + g_RT[13] - g_RT[24];
    Kc[136] = g_RT[9] - g_RT[14] + g_RT[15] - g_RT[16];
    Kc[137] = -g_RT[1] + g_RT[9] + g_RT[17] - g_RT[28];
    Kc[138] = g_RT[9] - g_RT[14] - g_RT[22] + g_RT[27];
    Kc[139] = g_RT[3] - g_RT[4] + g_RT[10] - g_RT[16];
    Kc[140] = g_RT[0] - g_RT[1] + g_RT[10] - g_RT[12];
    Kc[141] = -g_RT[0] + 2.000000*g_RT[10] - g_RT[22];
    Kc[142] = -g_RT[1] + g_RT[10] + g_RT[12] - g_RT[24];
    Kc[143] = g_RT[10] - 2.000000*g_RT[12] + g_RT[13];
    Kc[144] = g_RT[10] - g_RT[14] - g_RT[23] + g_RT[27];
    Kc[145] = -g_RT[10] + g_RT[11] + g_RT[30] - g_RT[30];
    Kc[146] = -g_RT[10] + g_RT[11] + g_RT[31] - g_RT[31];
    Kc[147] = -g_RT[1] + g_RT[3] - g_RT[4] + g_RT[11] - g_RT[14];
    Kc[148] = g_RT[3] - g_RT[5] + g_RT[11] - g_RT[14];
    Kc[149] = g_RT[0] - g_RT[1] + g_RT[11] - g_RT[12];
    Kc[150] = g_RT[5] - g_RT[5] - g_RT[10] + g_RT[11];
    Kc[151] = -g_RT[1] + g_RT[11] + g_RT[12] - g_RT[24];
    Kc[152] = g_RT[11] - 2.000000*g_RT[12] + g_RT[13];
    Kc[153] = -g_RT[10] + g_RT[11] + g_RT[14] - g_RT[14];
    Kc[154] = -g_RT[10] + g_RT[11] + g_RT[15] - g_RT[15];
    Kc[155] = g_RT[11] - g_RT[14] + g_RT[15] - g_RT[17];
    Kc[156] = g_RT[11] - g_RT[12] - g_RT[25] + g_RT[26];
    Kc[157] = -g_RT[2] + g_RT[3] + g_RT[12] - g_RT[19];
    Kc[158] = g_RT[3] - g_RT[4] + g_RT[12] - g_RT[17];
    Kc[159] = -g_RT[6] + g_RT[7] + g_RT[12] - g_RT[13];
    Kc[160] = -g_RT[1] + 2.000000*g_RT[12] - g_RT[25];
    Kc[161] = g_RT[12] - g_RT[13] - g_RT[14] + g_RT[16];
    Kc[162] = g_RT[12] - g_RT[13] - g_RT[16] + g_RT[17];
    Kc[163] = g_RT[12] - g_RT[13] - g_RT[18] + g_RT[20];
    Kc[164] = g_RT[12] - g_RT[13] - g_RT[19] + g_RT[20];
    Kc[165] = g_RT[12] - g_RT[13] - g_RT[23] + g_RT[24];
    Kc[166] = g_RT[12] - g_RT[13] - g_RT[25] + g_RT[26];
    Kc[167] = -g_RT[1] + g_RT[5] - g_RT[5] - g_RT[14] + g_RT[16];
    Kc[168] = g_RT[3] - g_RT[6] - g_RT[14] + g_RT[16];
    Kc[169] = g_RT[3] - g_RT[6] - g_RT[17] + g_RT[18];
    Kc[170] = g_RT[3] - g_RT[6] - g_RT[17] + g_RT[19];
    Kc[171] = g_RT[3] - g_RT[14] - g_RT[16] + g_RT[21];
    Kc[172] = g_RT[0] - g_RT[1] + g_RT[21] - g_RT[22];
    Kc[173] = g_RT[3] - g_RT[16] - g_RT[17] + g_RT[23];
    Kc[174] = g_RT[3] - g_RT[6] - g_RT[24] + g_RT[25];
    Kc[175] = g_RT[3] - g_RT[4] - 2.000000*g_RT[14] + g_RT[27];
    Kc[176] = -2.000000*g_RT[14] - g_RT[22] + 2.000000*g_RT[27];

    for (int i=0; i<177; ++i) {
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
    Kc[13] *= refCinv;
    Kc[14] *= refCinv;
    Kc[15] *= refCinv;
    Kc[16] *= refCinv;
    Kc[17] *= refCinv;
    Kc[18] *= refCinv;
    Kc[19] *= refC;
    Kc[20] *= refCinv;
    Kc[21] *= refCinv;
    Kc[22] *= refCinv;
    Kc[23] *= refCinv;
    Kc[24] *= refCinv;
    Kc[25] *= refCinv;
    Kc[26] *= refC;
    Kc[51] *= refC;
    Kc[56] *= refCinv;
    Kc[57] *= refCinv;
    Kc[58] *= refCinv;
    Kc[59] *= refCinv;
    Kc[61] *= refCinv;
    Kc[62] *= refCinv;
    Kc[63] *= refCinv;
    Kc[147] *= refC;
    Kc[167] *= refC;
    Kc[175] *= refC;
    Kc[176] *= refC;

    return;
}

void comp_qfqr(amrex::Real *  qf, amrex::Real *  qr, amrex::Real *  sc, amrex::Real * qss_sc, amrex::Real *  tc, amrex::Real invT)
{

    /*reaction 1: H + CH2 (+M) <=> CH3 (+M) */
    qf[0] = sc[1]*sc[10];
    qr[0] = sc[12];

    /*reaction 2: H + CH3 (+M) <=> CH4 (+M) */
    qf[1] = sc[1]*sc[12];
    qr[1] = sc[13];

    /*reaction 3: H + HCO (+M) <=> CH2O (+M) */
    qf[2] = sc[1]*sc[16];
    qr[2] = sc[17];

    /*reaction 4: H + CH2O (+M) <=> CH2OH (+M) */
    qf[3] = sc[1]*sc[17];
    qr[3] = sc[18];

    /*reaction 5: H + CH2O (+M) <=> CH3O (+M) */
    qf[4] = sc[1]*sc[17];
    qr[4] = sc[19];

    /*reaction 6: H + CH2OH (+M) <=> CH3OH (+M) */
    qf[5] = sc[1]*sc[18];
    qr[5] = sc[20];

    /*reaction 7: H + CH3O (+M) <=> CH3OH (+M) */
    qf[6] = sc[1]*sc[19];
    qr[6] = sc[20];

    /*reaction 8: H + C2H (+M) <=> C2H2 (+M) */
    qf[7] = sc[1]*sc[21];
    qr[7] = sc[22];

    /*reaction 9: H + C2H2 (+M) <=> C2H3 (+M) */
    qf[8] = sc[1]*sc[22];
    qr[8] = sc[23];

    /*reaction 10: H + C2H3 (+M) <=> C2H4 (+M) */
    qf[9] = sc[1]*sc[23];
    qr[9] = sc[24];

    /*reaction 11: H + C2H4 (+M) <=> C2H5 (+M) */
    qf[10] = sc[1]*sc[24];
    qr[10] = sc[25];

    /*reaction 12: H + C2H5 (+M) <=> C2H6 (+M) */
    qf[11] = sc[1]*sc[25];
    qr[11] = sc[26];

    /*reaction 13: H2 + CO (+M) <=> CH2O (+M) */
    qf[12] = sc[0]*sc[14];
    qr[12] = sc[17];

    /*reaction 14: 2.000000 OH (+M) <=> H2O2 (+M) */
    qf[13] = pow(sc[4], 2.000000);
    qr[13] = sc[7];

    /*reaction 15: OH + CH3 (+M) <=> CH3OH (+M) */
    qf[14] = sc[4]*sc[12];
    qr[14] = sc[20];

    /*reaction 16: CH + CO (+M) <=> HCCO (+M) */
    qf[15] = sc[9]*sc[14];
    qr[15] = sc[27];

    /*reaction 17: CH2 + CO (+M) <=> CH2CO (+M) */
    qf[16] = sc[10]*sc[14];
    qr[16] = sc[28];

    /*reaction 18: CH2(S) + H2O (+M) <=> CH3OH (+M) */
    qf[17] = sc[5]*sc[11];
    qr[17] = sc[20];

    /*reaction 19: 2.000000 CH3 (+M) <=> C2H6 (+M) */
    qf[18] = pow(sc[12], 2.000000);
    qr[18] = sc[26];

    /*reaction 20: C2H4 (+M) <=> H2 + C2H2 (+M) */
    qf[19] = sc[24];
    qr[19] = sc[0]*sc[22];

    /*reaction 21: 2.000000 O + M <=> O2 + M */
    qf[20] = pow(sc[2], 2.000000);
    qr[20] = sc[3];

    /*reaction 22: O + H + M <=> OH + M */
    qf[21] = sc[1]*sc[2];
    qr[21] = sc[4];

    /*reaction 23: O + CO + M <=> CO2 + M */
    qf[22] = sc[2]*sc[14];
    qr[22] = sc[15];

    /*reaction 24: H + O2 + M <=> HO2 + M */
    qf[23] = sc[1]*sc[3];
    qr[23] = sc[6];

    /*reaction 25: 2.000000 H + M <=> H2 + M */
    qf[24] = pow(sc[1], 2.000000);
    qr[24] = sc[0];

    /*reaction 26: H + OH + M <=> H2O + M */
    qf[25] = sc[1]*sc[4];
    qr[25] = sc[5];

    /*reaction 27: HCO + M <=> H + CO + M */
    qf[26] = sc[16];
    qr[26] = sc[1]*sc[14];

    /*reaction 28: O + H2 <=> H + OH */
    qf[27] = sc[0]*sc[2];
    qr[27] = sc[1]*sc[4];

    /*reaction 29: O + HO2 <=> OH + O2 */
    qf[28] = sc[2]*sc[6];
    qr[28] = sc[3]*sc[4];

    /*reaction 30: O + H2O2 <=> OH + HO2 */
    qf[29] = sc[2]*sc[7];
    qr[29] = sc[4]*sc[6];

    /*reaction 31: O + CH <=> H + CO */
    qf[30] = sc[2]*sc[9];
    qr[30] = sc[1]*sc[14];

    /*reaction 32: O + CH2 <=> H + HCO */
    qf[31] = sc[2]*sc[10];
    qr[31] = sc[1]*sc[16];

    /*reaction 33: O + CH2(S) <=> H2 + CO */
    qf[32] = sc[2]*sc[11];
    qr[32] = sc[0]*sc[14];

    /*reaction 34: O + CH2(S) <=> H + HCO */
    qf[33] = sc[2]*sc[11];
    qr[33] = sc[1]*sc[16];

    /*reaction 35: O + CH3 <=> H + CH2O */
    qf[34] = sc[2]*sc[12];
    qr[34] = sc[1]*sc[17];

    /*reaction 36: O + CH4 <=> OH + CH3 */
    qf[35] = sc[2]*sc[13];
    qr[35] = sc[4]*sc[12];

    /*reaction 37: O + HCO <=> OH + CO */
    qf[36] = sc[2]*sc[16];
    qr[36] = sc[4]*sc[14];

    /*reaction 38: O + HCO <=> H + CO2 */
    qf[37] = sc[2]*sc[16];
    qr[37] = sc[1]*sc[15];

    /*reaction 39: O + CH2O <=> OH + HCO */
    qf[38] = sc[2]*sc[17];
    qr[38] = sc[4]*sc[16];

    /*reaction 40: O + CH2OH <=> OH + CH2O */
    qf[39] = sc[2]*sc[18];
    qr[39] = sc[4]*sc[17];

    /*reaction 41: O + CH3O <=> OH + CH2O */
    qf[40] = sc[2]*sc[19];
    qr[40] = sc[4]*sc[17];

    /*reaction 42: O + CH3OH <=> OH + CH2OH */
    qf[41] = sc[2]*sc[20];
    qr[41] = sc[4]*sc[18];

    /*reaction 43: O + CH3OH <=> OH + CH3O */
    qf[42] = sc[2]*sc[20];
    qr[42] = sc[4]*sc[19];

    /*reaction 44: O + C2H <=> CH + CO */
    qf[43] = sc[2]*sc[21];
    qr[43] = sc[9]*sc[14];

    /*reaction 45: O + C2H2 <=> H + HCCO */
    qf[44] = sc[2]*sc[22];
    qr[44] = sc[1]*sc[27];

    /*reaction 46: O + C2H2 <=> OH + C2H */
    qf[45] = sc[2]*sc[22];
    qr[45] = sc[4]*sc[21];

    /*reaction 47: O + C2H2 <=> CO + CH2 */
    qf[46] = sc[2]*sc[22];
    qr[46] = sc[10]*sc[14];

    /*reaction 48: O + C2H3 <=> H + CH2CO */
    qf[47] = sc[2]*sc[23];
    qr[47] = sc[1]*sc[28];

    /*reaction 49: O + C2H4 <=> CH3 + HCO */
    qf[48] = sc[2]*sc[24];
    qr[48] = sc[12]*sc[16];

    /*reaction 50: O + C2H5 <=> CH3 + CH2O */
    qf[49] = sc[2]*sc[25];
    qr[49] = sc[12]*sc[17];

    /*reaction 51: O + C2H6 <=> OH + C2H5 */
    qf[50] = sc[2]*sc[26];
    qr[50] = sc[4]*sc[25];

    /*reaction 52: O + HCCO <=> H + 2.000000 CO */
    qf[51] = sc[2]*sc[27];
    qr[51] = sc[1]*pow(sc[14], 2.000000);

    /*reaction 53: O + CH2CO <=> OH + HCCO */
    qf[52] = sc[2]*sc[28];
    qr[52] = sc[4]*sc[27];

    /*reaction 54: O + CH2CO <=> CH2 + CO2 */
    qf[53] = sc[2]*sc[28];
    qr[53] = sc[10]*sc[15];

    /*reaction 55: O2 + CO <=> O + CO2 */
    qf[54] = sc[3]*sc[14];
    qr[54] = sc[2]*sc[15];

    /*reaction 56: O2 + CH2O <=> HO2 + HCO */
    qf[55] = sc[3]*sc[17];
    qr[55] = sc[6]*sc[16];

    /*reaction 57: H + 2.000000 O2 <=> HO2 + O2 */
    qf[56] = sc[1]*pow(sc[3], 2.000000);
    qr[56] = sc[3]*sc[6];

    /*reaction 58: H + O2 + H2O <=> HO2 + H2O */
    qf[57] = sc[1]*sc[3]*sc[5];
    qr[57] = sc[5]*sc[6];

    /*reaction 59: H + O2 + N2 <=> HO2 + N2 */
    qf[58] = sc[1]*sc[3]*sc[30];
    qr[58] = sc[6]*sc[30];

    /*reaction 60: H + O2 + AR <=> HO2 + AR */
    qf[59] = sc[1]*sc[3]*sc[31];
    qr[59] = sc[6]*sc[31];

    /*reaction 61: H + O2 <=> O + OH */
    qf[60] = sc[1]*sc[3];
    qr[60] = sc[2]*sc[4];

    /*reaction 62: 2.000000 H + H2 <=> 2.000000 H2 */
    qf[61] = sc[0]*pow(sc[1], 2.000000);
    qr[61] = pow(sc[0], 2.000000);

    /*reaction 63: 2.000000 H + H2O <=> H2 + H2O */
    qf[62] = pow(sc[1], 2.000000)*sc[5];
    qr[62] = sc[0]*sc[5];

    /*reaction 64: 2.000000 H + CO2 <=> H2 + CO2 */
    qf[63] = pow(sc[1], 2.000000)*sc[15];
    qr[63] = sc[0]*sc[15];

    /*reaction 65: H + HO2 <=> O + H2O */
    qf[64] = sc[1]*sc[6];
    qr[64] = sc[2]*sc[5];

    /*reaction 66: H + HO2 <=> O2 + H2 */
    qf[65] = sc[1]*sc[6];
    qr[65] = sc[0]*sc[3];

    /*reaction 67: H + HO2 <=> 2.000000 OH */
    qf[66] = sc[1]*sc[6];
    qr[66] = pow(sc[4], 2.000000);

    /*reaction 68: H + H2O2 <=> HO2 + H2 */
    qf[67] = sc[1]*sc[7];
    qr[67] = sc[0]*sc[6];

    /*reaction 69: H + H2O2 <=> OH + H2O */
    qf[68] = sc[1]*sc[7];
    qr[68] = sc[4]*sc[5];

    /*reaction 70: H + CH <=> C + H2 */
    qf[69] = sc[1]*sc[9];
    qr[69] = sc[0]*sc[8];

    /*reaction 71: H + CH2(S) <=> CH + H2 */
    qf[70] = sc[1]*sc[11];
    qr[70] = sc[0]*sc[9];

    /*reaction 72: H + CH4 <=> CH3 + H2 */
    qf[71] = sc[1]*sc[13];
    qr[71] = sc[0]*sc[12];

    /*reaction 73: H + HCO <=> H2 + CO */
    qf[72] = sc[1]*sc[16];
    qr[72] = sc[0]*sc[14];

    /*reaction 74: H + CH2O <=> HCO + H2 */
    qf[73] = sc[1]*sc[17];
    qr[73] = sc[0]*sc[16];

    /*reaction 75: H + CH2OH <=> H2 + CH2O */
    qf[74] = sc[1]*sc[18];
    qr[74] = sc[0]*sc[17];

    /*reaction 76: H + CH2OH <=> OH + CH3 */
    qf[75] = sc[1]*sc[18];
    qr[75] = sc[4]*sc[12];

    /*reaction 77: H + CH2OH <=> CH2(S) + H2O */
    qf[76] = sc[1]*sc[18];
    qr[76] = sc[5]*sc[11];

    /*reaction 78: H + CH3O <=> H + CH2OH */
    qf[77] = sc[1]*sc[19];
    qr[77] = sc[1]*sc[18];

    /*reaction 79: H + CH3O <=> H2 + CH2O */
    qf[78] = sc[1]*sc[19];
    qr[78] = sc[0]*sc[17];

    /*reaction 80: H + CH3O <=> OH + CH3 */
    qf[79] = sc[1]*sc[19];
    qr[79] = sc[4]*sc[12];

    /*reaction 81: H + CH3O <=> CH2(S) + H2O */
    qf[80] = sc[1]*sc[19];
    qr[80] = sc[5]*sc[11];

    /*reaction 82: H + CH3OH <=> CH2OH + H2 */
    qf[81] = sc[1]*sc[20];
    qr[81] = sc[0]*sc[18];

    /*reaction 83: H + CH3OH <=> CH3O + H2 */
    qf[82] = sc[1]*sc[20];
    qr[82] = sc[0]*sc[19];

    /*reaction 84: H + C2H3 <=> H2 + C2H2 */
    qf[83] = sc[1]*sc[23];
    qr[83] = sc[0]*sc[22];

    /*reaction 85: H + C2H4 <=> C2H3 + H2 */
    qf[84] = sc[1]*sc[24];
    qr[84] = sc[0]*sc[23];

    /*reaction 86: H + C2H5 <=> H2 + C2H4 */
    qf[85] = sc[1]*sc[25];
    qr[85] = sc[0]*sc[24];

    /*reaction 87: H + C2H6 <=> C2H5 + H2 */
    qf[86] = sc[1]*sc[26];
    qr[86] = sc[0]*sc[25];

    /*reaction 88: H + HCCO <=> CH2(S) + CO */
    qf[87] = sc[1]*sc[27];
    qr[87] = sc[11]*sc[14];

    /*reaction 89: H + CH2CO <=> HCCO + H2 */
    qf[88] = sc[1]*sc[28];
    qr[88] = sc[0]*sc[27];

    /*reaction 90: H + CH2CO <=> CH3 + CO */
    qf[89] = sc[1]*sc[28];
    qr[89] = sc[12]*sc[14];

    /*reaction 91: H + HCCOH <=> H + CH2CO */
    qf[90] = sc[1]*sc[29];
    qr[90] = sc[1]*sc[28];

    /*reaction 92: OH + H2 <=> H + H2O */
    qf[91] = sc[0]*sc[4];
    qr[91] = sc[1]*sc[5];

    /*reaction 93: 2.000000 OH <=> O + H2O */
    qf[92] = pow(sc[4], 2.000000);
    qr[92] = sc[2]*sc[5];

    /*reaction 94: OH + HO2 <=> O2 + H2O */
    qf[93] = sc[4]*sc[6];
    qr[93] = sc[3]*sc[5];

    /*reaction 95: OH + H2O2 <=> HO2 + H2O */
    qf[94] = sc[4]*sc[7];
    qr[94] = sc[5]*sc[6];

    /*reaction 96: OH + H2O2 <=> HO2 + H2O */
    qf[95] = sc[4]*sc[7];
    qr[95] = sc[5]*sc[6];

    /*reaction 97: OH + C <=> H + CO */
    qf[96] = sc[4]*sc[8];
    qr[96] = sc[1]*sc[14];

    /*reaction 98: OH + CH <=> H + HCO */
    qf[97] = sc[4]*sc[9];
    qr[97] = sc[1]*sc[16];

    /*reaction 99: OH + CH2 <=> H + CH2O */
    qf[98] = sc[4]*sc[10];
    qr[98] = sc[1]*sc[17];

    /*reaction 100: OH + CH2 <=> CH + H2O */
    qf[99] = sc[4]*sc[10];
    qr[99] = sc[5]*sc[9];

    /*reaction 101: OH + CH2(S) <=> H + CH2O */
    qf[100] = sc[4]*sc[11];
    qr[100] = sc[1]*sc[17];

    /*reaction 102: OH + CH3 <=> CH2 + H2O */
    qf[101] = sc[4]*sc[12];
    qr[101] = sc[5]*sc[10];

    /*reaction 103: OH + CH3 <=> CH2(S) + H2O */
    qf[102] = sc[4]*sc[12];
    qr[102] = sc[5]*sc[11];

    /*reaction 104: OH + CH4 <=> CH3 + H2O */
    qf[103] = sc[4]*sc[13];
    qr[103] = sc[5]*sc[12];

    /*reaction 105: OH + CO <=> H + CO2 */
    qf[104] = sc[4]*sc[14];
    qr[104] = sc[1]*sc[15];

    /*reaction 106: OH + HCO <=> H2O + CO */
    qf[105] = sc[4]*sc[16];
    qr[105] = sc[5]*sc[14];

    /*reaction 107: OH + CH2O <=> HCO + H2O */
    qf[106] = sc[4]*sc[17];
    qr[106] = sc[5]*sc[16];

    /*reaction 108: OH + CH2OH <=> H2O + CH2O */
    qf[107] = sc[4]*sc[18];
    qr[107] = sc[5]*sc[17];

    /*reaction 109: OH + CH3O <=> H2O + CH2O */
    qf[108] = sc[4]*sc[19];
    qr[108] = sc[5]*sc[17];

    /*reaction 110: OH + CH3OH <=> CH2OH + H2O */
    qf[109] = sc[4]*sc[20];
    qr[109] = sc[5]*sc[18];

    /*reaction 111: OH + CH3OH <=> CH3O + H2O */
    qf[110] = sc[4]*sc[20];
    qr[110] = sc[5]*sc[19];

    /*reaction 112: OH + C2H <=> H + HCCO */
    qf[111] = sc[4]*sc[21];
    qr[111] = sc[1]*sc[27];

    /*reaction 113: OH + C2H2 <=> H + CH2CO */
    qf[112] = sc[4]*sc[22];
    qr[112] = sc[1]*sc[28];

    /*reaction 114: OH + C2H2 <=> H + HCCOH */
    qf[113] = sc[4]*sc[22];
    qr[113] = sc[1]*sc[29];

    /*reaction 115: OH + C2H2 <=> C2H + H2O */
    qf[114] = sc[4]*sc[22];
    qr[114] = sc[5]*sc[21];

    /*reaction 116: OH + C2H2 <=> CH3 + CO */
    qf[115] = sc[4]*sc[22];
    qr[115] = sc[12]*sc[14];

    /*reaction 117: OH + C2H3 <=> H2O + C2H2 */
    qf[116] = sc[4]*sc[23];
    qr[116] = sc[5]*sc[22];

    /*reaction 118: OH + C2H4 <=> C2H3 + H2O */
    qf[117] = sc[4]*sc[24];
    qr[117] = sc[5]*sc[23];

    /*reaction 119: OH + C2H6 <=> C2H5 + H2O */
    qf[118] = sc[4]*sc[26];
    qr[118] = sc[5]*sc[25];

    /*reaction 120: OH + CH2CO <=> HCCO + H2O */
    qf[119] = sc[4]*sc[28];
    qr[119] = sc[5]*sc[27];

    /*reaction 121: 2.000000 HO2 <=> O2 + H2O2 */
    qf[120] = pow(sc[6], 2.000000);
    qr[120] = sc[3]*sc[7];

    /*reaction 122: 2.000000 HO2 <=> O2 + H2O2 */
    qf[121] = pow(sc[6], 2.000000);
    qr[121] = sc[3]*sc[7];

    /*reaction 123: HO2 + CH2 <=> OH + CH2O */
    qf[122] = sc[6]*sc[10];
    qr[122] = sc[4]*sc[17];

    /*reaction 124: HO2 + CH3 <=> O2 + CH4 */
    qf[123] = sc[6]*sc[12];
    qr[123] = sc[3]*sc[13];

    /*reaction 125: HO2 + CH3 <=> OH + CH3O */
    qf[124] = sc[6]*sc[12];
    qr[124] = sc[4]*sc[19];

    /*reaction 126: HO2 + CO <=> OH + CO2 */
    qf[125] = sc[6]*sc[14];
    qr[125] = sc[4]*sc[15];

    /*reaction 127: HO2 + CH2O <=> HCO + H2O2 */
    qf[126] = sc[6]*sc[17];
    qr[126] = sc[7]*sc[16];

    /*reaction 128: C + O2 <=> O + CO */
    qf[127] = sc[3]*sc[8];
    qr[127] = sc[2]*sc[14];

    /*reaction 129: C + CH2 <=> H + C2H */
    qf[128] = sc[8]*sc[10];
    qr[128] = sc[1]*sc[21];

    /*reaction 130: C + CH3 <=> H + C2H2 */
    qf[129] = sc[8]*sc[12];
    qr[129] = sc[1]*sc[22];

    /*reaction 131: CH + O2 <=> O + HCO */
    qf[130] = sc[3]*sc[9];
    qr[130] = sc[2]*sc[16];

    /*reaction 132: CH + H2 <=> H + CH2 */
    qf[131] = sc[0]*sc[9];
    qr[131] = sc[1]*sc[10];

    /*reaction 133: CH + H2O <=> H + CH2O */
    qf[132] = sc[5]*sc[9];
    qr[132] = sc[1]*sc[17];

    /*reaction 134: CH + CH2 <=> H + C2H2 */
    qf[133] = sc[9]*sc[10];
    qr[133] = sc[1]*sc[22];

    /*reaction 135: CH + CH3 <=> H + C2H3 */
    qf[134] = sc[9]*sc[12];
    qr[134] = sc[1]*sc[23];

    /*reaction 136: CH + CH4 <=> H + C2H4 */
    qf[135] = sc[9]*sc[13];
    qr[135] = sc[1]*sc[24];

    /*reaction 137: CH + CO2 <=> HCO + CO */
    qf[136] = sc[9]*sc[15];
    qr[136] = sc[14]*sc[16];

    /*reaction 138: CH + CH2O <=> H + CH2CO */
    qf[137] = sc[9]*sc[17];
    qr[137] = sc[1]*sc[28];

    /*reaction 139: CH + HCCO <=> CO + C2H2 */
    qf[138] = sc[9]*sc[27];
    qr[138] = sc[14]*sc[22];

    /*reaction 140: CH2 + O2 <=> OH + HCO */
    qf[139] = sc[3]*sc[10];
    qr[139] = sc[4]*sc[16];

    /*reaction 141: CH2 + H2 <=> H + CH3 */
    qf[140] = sc[0]*sc[10];
    qr[140] = sc[1]*sc[12];

    /*reaction 142: 2.000000 CH2 <=> H2 + C2H2 */
    qf[141] = pow(sc[10], 2.000000);
    qr[141] = sc[0]*sc[22];

    /*reaction 143: CH2 + CH3 <=> H + C2H4 */
    qf[142] = sc[10]*sc[12];
    qr[142] = sc[1]*sc[24];

    /*reaction 144: CH2 + CH4 <=> 2.000000 CH3 */
    qf[143] = sc[10]*sc[13];
    qr[143] = pow(sc[12], 2.000000);

    /*reaction 145: CH2 + HCCO <=> C2H3 + CO */
    qf[144] = sc[10]*sc[27];
    qr[144] = sc[14]*sc[23];

    /*reaction 146: CH2(S) + N2 <=> CH2 + N2 */
    qf[145] = sc[11]*sc[30];
    qr[145] = sc[10]*sc[30];

    /*reaction 147: CH2(S) + AR <=> CH2 + AR */
    qf[146] = sc[11]*sc[31];
    qr[146] = sc[10]*sc[31];

    /*reaction 148: CH2(S) + O2 <=> H + OH + CO */
    qf[147] = sc[3]*sc[11];
    qr[147] = sc[1]*sc[4]*sc[14];

    /*reaction 149: CH2(S) + O2 <=> CO + H2O */
    qf[148] = sc[3]*sc[11];
    qr[148] = sc[5]*sc[14];

    /*reaction 150: CH2(S) + H2 <=> CH3 + H */
    qf[149] = sc[0]*sc[11];
    qr[149] = sc[1]*sc[12];

    /*reaction 151: CH2(S) + H2O <=> CH2 + H2O */
    qf[150] = sc[5]*sc[11];
    qr[150] = sc[5]*sc[10];

    /*reaction 152: CH2(S) + CH3 <=> H + C2H4 */
    qf[151] = sc[11]*sc[12];
    qr[151] = sc[1]*sc[24];

    /*reaction 153: CH2(S) + CH4 <=> 2.000000 CH3 */
    qf[152] = sc[11]*sc[13];
    qr[152] = pow(sc[12], 2.000000);

    /*reaction 154: CH2(S) + CO <=> CH2 + CO */
    qf[153] = sc[11]*sc[14];
    qr[153] = sc[10]*sc[14];

    /*reaction 155: CH2(S) + CO2 <=> CH2 + CO2 */
    qf[154] = sc[11]*sc[15];
    qr[154] = sc[10]*sc[15];

    /*reaction 156: CH2(S) + CO2 <=> CO + CH2O */
    qf[155] = sc[11]*sc[15];
    qr[155] = sc[14]*sc[17];

    /*reaction 157: CH2(S) + C2H6 <=> CH3 + C2H5 */
    qf[156] = sc[11]*sc[26];
    qr[156] = sc[12]*sc[25];

    /*reaction 158: CH3 + O2 <=> O + CH3O */
    qf[157] = sc[3]*sc[12];
    qr[157] = sc[2]*sc[19];

    /*reaction 159: CH3 + O2 <=> OH + CH2O */
    qf[158] = sc[3]*sc[12];
    qr[158] = sc[4]*sc[17];

    /*reaction 160: CH3 + H2O2 <=> HO2 + CH4 */
    qf[159] = sc[7]*sc[12];
    qr[159] = sc[6]*sc[13];

    /*reaction 161: 2.000000 CH3 <=> H + C2H5 */
    qf[160] = pow(sc[12], 2.000000);
    qr[160] = sc[1]*sc[25];

    /*reaction 162: CH3 + HCO <=> CH4 + CO */
    qf[161] = sc[12]*sc[16];
    qr[161] = sc[13]*sc[14];

    /*reaction 163: CH3 + CH2O <=> HCO + CH4 */
    qf[162] = sc[12]*sc[17];
    qr[162] = sc[13]*sc[16];

    /*reaction 164: CH3 + CH3OH <=> CH2OH + CH4 */
    qf[163] = sc[12]*sc[20];
    qr[163] = sc[13]*sc[18];

    /*reaction 165: CH3 + CH3OH <=> CH3O + CH4 */
    qf[164] = sc[12]*sc[20];
    qr[164] = sc[13]*sc[19];

    /*reaction 166: CH3 + C2H4 <=> C2H3 + CH4 */
    qf[165] = sc[12]*sc[24];
    qr[165] = sc[13]*sc[23];

    /*reaction 167: CH3 + C2H6 <=> C2H5 + CH4 */
    qf[166] = sc[12]*sc[26];
    qr[166] = sc[13]*sc[25];

    /*reaction 168: HCO + H2O <=> H + CO + H2O */
    qf[167] = sc[5]*sc[16];
    qr[167] = sc[1]*sc[5]*sc[14];

    /*reaction 169: HCO + O2 <=> HO2 + CO */
    qf[168] = sc[3]*sc[16];
    qr[168] = sc[6]*sc[14];

    /*reaction 170: CH2OH + O2 <=> HO2 + CH2O */
    qf[169] = sc[3]*sc[18];
    qr[169] = sc[6]*sc[17];

    /*reaction 171: CH3O + O2 <=> HO2 + CH2O */
    qf[170] = sc[3]*sc[19];
    qr[170] = sc[6]*sc[17];

    /*reaction 172: C2H + O2 <=> HCO + CO */
    qf[171] = sc[3]*sc[21];
    qr[171] = sc[14]*sc[16];

    /*reaction 173: C2H + H2 <=> H + C2H2 */
    qf[172] = sc[0]*sc[21];
    qr[172] = sc[1]*sc[22];

    /*reaction 174: C2H3 + O2 <=> HCO + CH2O */
    qf[173] = sc[3]*sc[23];
    qr[173] = sc[16]*sc[17];

    /*reaction 175: C2H5 + O2 <=> HO2 + C2H4 */
    qf[174] = sc[3]*sc[25];
    qr[174] = sc[6]*sc[24];

    /*reaction 176: HCCO + O2 <=> OH + 2.000000 CO */
    qf[175] = sc[3]*sc[27];
    qr[175] = sc[4]*pow(sc[14], 2.000000);

    /*reaction 177: 2.000000 HCCO <=> 2.000000 CO + C2H2 */
    qf[176] = pow(sc[27], 2.000000);
    qr[176] = pow(sc[14], 2.000000)*sc[22];

    amrex::Real T = tc[1];

    /*compute the mixture concentration */
    amrex::Real mixture = 0.0;
    for (int i = 0; i < 32; ++i) {
        mixture += sc[i];
    }

    amrex::Real Corr[177];
    for (int i = 0; i < 177; ++i) {
        Corr[i] = 1.0;
    }

    /* troe */
    {
        amrex::Real alpha[20];
        alpha[0] = mixture + (TB[0][0] - 1)*sc[0] + (TB[0][1] - 1)*sc[5] + (TB[0][2] - 1)*sc[13] + (TB[0][3] - 1)*sc[14] + (TB[0][4] - 1)*sc[15] + (TB[0][5] - 1)*sc[26] + (TB[0][6] - 1)*sc[31];
        alpha[1] = mixture + (TB[1][0] - 1)*sc[0] + (TB[1][1] - 1)*sc[5] + (TB[1][2] - 1)*sc[13] + (TB[1][3] - 1)*sc[14] + (TB[1][4] - 1)*sc[15] + (TB[1][5] - 1)*sc[26] + (TB[1][6] - 1)*sc[31];
        alpha[2] = mixture + (TB[2][0] - 1)*sc[0] + (TB[2][1] - 1)*sc[5] + (TB[2][2] - 1)*sc[13] + (TB[2][3] - 1)*sc[14] + (TB[2][4] - 1)*sc[15] + (TB[2][5] - 1)*sc[26] + (TB[2][6] - 1)*sc[31];
        alpha[3] = mixture + (TB[3][0] - 1)*sc[0] + (TB[3][1] - 1)*sc[5] + (TB[3][2] - 1)*sc[13] + (TB[3][3] - 1)*sc[14] + (TB[3][4] - 1)*sc[15] + (TB[3][5] - 1)*sc[26];
        alpha[4] = mixture + (TB[4][0] - 1)*sc[0] + (TB[4][1] - 1)*sc[5] + (TB[4][2] - 1)*sc[13] + (TB[4][3] - 1)*sc[14] + (TB[4][4] - 1)*sc[15] + (TB[4][5] - 1)*sc[26];
        alpha[5] = mixture + (TB[5][0] - 1)*sc[0] + (TB[5][1] - 1)*sc[5] + (TB[5][2] - 1)*sc[13] + (TB[5][3] - 1)*sc[14] + (TB[5][4] - 1)*sc[15] + (TB[5][5] - 1)*sc[26];
        alpha[6] = mixture + (TB[6][0] - 1)*sc[0] + (TB[6][1] - 1)*sc[5] + (TB[6][2] - 1)*sc[13] + (TB[6][3] - 1)*sc[14] + (TB[6][4] - 1)*sc[15] + (TB[6][5] - 1)*sc[26];
        alpha[7] = mixture + (TB[7][0] - 1)*sc[0] + (TB[7][1] - 1)*sc[5] + (TB[7][2] - 1)*sc[13] + (TB[7][3] - 1)*sc[14] + (TB[7][4] - 1)*sc[15] + (TB[7][5] - 1)*sc[26] + (TB[7][6] - 1)*sc[31];
        alpha[8] = mixture + (TB[8][0] - 1)*sc[0] + (TB[8][1] - 1)*sc[5] + (TB[8][2] - 1)*sc[13] + (TB[8][3] - 1)*sc[14] + (TB[8][4] - 1)*sc[15] + (TB[8][5] - 1)*sc[26] + (TB[8][6] - 1)*sc[31];
        alpha[9] = mixture + (TB[9][0] - 1)*sc[0] + (TB[9][1] - 1)*sc[5] + (TB[9][2] - 1)*sc[13] + (TB[9][3] - 1)*sc[14] + (TB[9][4] - 1)*sc[15] + (TB[9][5] - 1)*sc[26] + (TB[9][6] - 1)*sc[31];
        alpha[10] = mixture + (TB[10][0] - 1)*sc[0] + (TB[10][1] - 1)*sc[5] + (TB[10][2] - 1)*sc[13] + (TB[10][3] - 1)*sc[14] + (TB[10][4] - 1)*sc[15] + (TB[10][5] - 1)*sc[26] + (TB[10][6] - 1)*sc[31];
        alpha[11] = mixture + (TB[11][0] - 1)*sc[0] + (TB[11][1] - 1)*sc[5] + (TB[11][2] - 1)*sc[13] + (TB[11][3] - 1)*sc[14] + (TB[11][4] - 1)*sc[15] + (TB[11][5] - 1)*sc[26] + (TB[11][6] - 1)*sc[31];
        alpha[12] = mixture + (TB[12][0] - 1)*sc[0] + (TB[12][1] - 1)*sc[5] + (TB[12][2] - 1)*sc[13] + (TB[12][3] - 1)*sc[14] + (TB[12][4] - 1)*sc[15] + (TB[12][5] - 1)*sc[26] + (TB[12][6] - 1)*sc[31];
        alpha[13] = mixture + (TB[13][0] - 1)*sc[0] + (TB[13][1] - 1)*sc[5] + (TB[13][2] - 1)*sc[13] + (TB[13][3] - 1)*sc[14] + (TB[13][4] - 1)*sc[15] + (TB[13][5] - 1)*sc[26] + (TB[13][6] - 1)*sc[31];
        alpha[14] = mixture + (TB[14][0] - 1)*sc[0] + (TB[14][1] - 1)*sc[5] + (TB[14][2] - 1)*sc[13] + (TB[14][3] - 1)*sc[14] + (TB[14][4] - 1)*sc[15] + (TB[14][5] - 1)*sc[26];
        alpha[15] = mixture + (TB[15][0] - 1)*sc[0] + (TB[15][1] - 1)*sc[5] + (TB[15][2] - 1)*sc[13] + (TB[15][3] - 1)*sc[14] + (TB[15][4] - 1)*sc[15] + (TB[15][5] - 1)*sc[26] + (TB[15][6] - 1)*sc[31];
        alpha[16] = mixture + (TB[16][0] - 1)*sc[0] + (TB[16][1] - 1)*sc[5] + (TB[16][2] - 1)*sc[13] + (TB[16][3] - 1)*sc[14] + (TB[16][4] - 1)*sc[15] + (TB[16][5] - 1)*sc[26] + (TB[16][6] - 1)*sc[31];
        alpha[17] = mixture + (TB[17][0] - 1)*sc[0] + (TB[17][1] - 1)*sc[5] + (TB[17][2] - 1)*sc[13] + (TB[17][3] - 1)*sc[14] + (TB[17][4] - 1)*sc[15] + (TB[17][5] - 1)*sc[26];
        alpha[18] = mixture + (TB[18][0] - 1)*sc[0] + (TB[18][1] - 1)*sc[5] + (TB[18][2] - 1)*sc[13] + (TB[18][3] - 1)*sc[14] + (TB[18][4] - 1)*sc[15] + (TB[18][5] - 1)*sc[26] + (TB[18][6] - 1)*sc[31];
        alpha[19] = mixture + (TB[19][0] - 1)*sc[0] + (TB[19][1] - 1)*sc[5] + (TB[19][2] - 1)*sc[13] + (TB[19][3] - 1)*sc[14] + (TB[19][4] - 1)*sc[15] + (TB[19][5] - 1)*sc[26] + (TB[19][6] - 1)*sc[31];
        for (int i=0; i<20; i++)
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
        alpha = mixture + (TB[20][0] - 1)*sc[0] + (TB[20][1] - 1)*sc[5] + (TB[20][2] - 1)*sc[13] + (TB[20][3] - 1)*sc[14] + (TB[20][4] - 1)*sc[15] + (TB[20][5] - 1)*sc[26] + (TB[20][6] - 1)*sc[31];
        Corr[20] = alpha;
        alpha = mixture + (TB[21][0] - 1)*sc[0] + (TB[21][1] - 1)*sc[5] + (TB[21][2] - 1)*sc[13] + (TB[21][3] - 1)*sc[14] + (TB[21][4] - 1)*sc[15] + (TB[21][5] - 1)*sc[26] + (TB[21][6] - 1)*sc[31];
        Corr[21] = alpha;
        alpha = mixture + (TB[22][0] - 1)*sc[0] + (TB[22][1] - 1)*sc[3] + (TB[22][2] - 1)*sc[5] + (TB[22][3] - 1)*sc[13] + (TB[22][4] - 1)*sc[14] + (TB[22][5] - 1)*sc[15] + (TB[22][6] - 1)*sc[26] + (TB[22][7] - 1)*sc[31];
        Corr[22] = alpha;
        alpha = mixture + (TB[23][0] - 1)*sc[3] + (TB[23][1] - 1)*sc[5] + (TB[23][2] - 1)*sc[14] + (TB[23][3] - 1)*sc[15] + (TB[23][4] - 1)*sc[26] + (TB[23][5] - 1)*sc[30] + (TB[23][6] - 1)*sc[31];
        Corr[23] = alpha;
        alpha = mixture + (TB[24][0] - 1)*sc[0] + (TB[24][1] - 1)*sc[5] + (TB[24][2] - 1)*sc[13] + (TB[24][3] - 1)*sc[15] + (TB[24][4] - 1)*sc[26] + (TB[24][5] - 1)*sc[31];
        Corr[24] = alpha;
        alpha = mixture + (TB[25][0] - 1)*sc[0] + (TB[25][1] - 1)*sc[5] + (TB[25][2] - 1)*sc[13] + (TB[25][3] - 1)*sc[26] + (TB[25][4] - 1)*sc[31];
        Corr[25] = alpha;
        alpha = mixture + (TB[26][0] - 1)*sc[0] + (TB[26][1] - 1)*sc[5] + (TB[26][2] - 1)*sc[13] + (TB[26][3] - 1)*sc[14] + (TB[26][4] - 1)*sc[15] + (TB[26][5] - 1)*sc[26];
        Corr[26] = alpha;
    }

    for (int i=0; i<177; i++)
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

    amrex::Real q_f[177], q_r[177];
    amrex::Real sc_qss[1];
    comp_qfqr(q_f, q_r, sc, sc_qss, tc, invT);

    for (int i = 0; i < 177; ++i) {
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
    amrex::Real c[32]; /*temporary storage */
    amrex::Real PORT = 1e6 * (*P)/(8.31446261815324e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 32; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRateFR(q_f, q_r, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 177; ++id) {
        q_f[id] *= 1.0e-6;
        q_r[id] *= 1.0e-6;
    }
}


/*Returns the rate of progress for each reaction */
void CKQC(amrex::Real *  T, amrex::Real *  C, amrex::Real *  qdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 32; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    progressRate(qdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 32; ++id) {
        C[id] *= 1.0e-6;
    }

    for (id = 0; id < 177; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mass fractions */
void CKQYP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y, amrex::Real *  qdot)
{
    int id; /*loop counter */
    amrex::Real c[32]; /*temporary storage */
    amrex::Real YOW = 0; 
    amrex::Real PWORT; 
    amrex::Real imw[32];

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
    YOW += y[8]*imw[8]; /*C */
    YOW += y[9]*imw[9]; /*CH */
    YOW += y[10]*imw[10]; /*CH2 */
    YOW += y[11]*imw[11]; /*CH2(S) */
    YOW += y[12]*imw[12]; /*CH3 */
    YOW += y[13]*imw[13]; /*CH4 */
    YOW += y[14]*imw[14]; /*CO */
    YOW += y[15]*imw[15]; /*CO2 */
    YOW += y[16]*imw[16]; /*HCO */
    YOW += y[17]*imw[17]; /*CH2O */
    YOW += y[18]*imw[18]; /*CH2OH */
    YOW += y[19]*imw[19]; /*CH3O */
    YOW += y[20]*imw[20]; /*CH3OH */
    YOW += y[21]*imw[21]; /*C2H */
    YOW += y[22]*imw[22]; /*C2H2 */
    YOW += y[23]*imw[23]; /*C2H3 */
    YOW += y[24]*imw[24]; /*C2H4 */
    YOW += y[25]*imw[25]; /*C2H5 */
    YOW += y[26]*imw[26]; /*C2H6 */
    YOW += y[27]*imw[27]; /*HCCO */
    YOW += y[28]*imw[28]; /*CH2CO */
    YOW += y[29]*imw[29]; /*HCCOH */
    YOW += y[30]*imw[30]; /*N2 */
    YOW += y[31]*imw[31]; /*AR */
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

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 177; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKQXP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x, amrex::Real *  qdot)
{
    int id; /*loop counter */
    amrex::Real c[32]; /*temporary storage */
    amrex::Real PORT = 1e6 * (*P)/(8.31446261815324e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 32; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 177; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mass fractions */
void CKQYR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y, amrex::Real *  qdot)
{
    int id; /*loop counter */
    amrex::Real c[32]; /*temporary storage */
    amrex::Real imw[32];

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

    /*call progressRate */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 177; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mole fractions */
void CKQXR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  x, amrex::Real *  qdot)
{
    int id; /*loop counter */
    amrex::Real c[32]; /*temporary storage */
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
    XW += x[8]*12.011150; /*C */
    XW += x[9]*13.019120; /*CH */
    XW += x[10]*14.027090; /*CH2 */
    XW += x[11]*14.027090; /*CH2(S) */
    XW += x[12]*15.035060; /*CH3 */
    XW += x[13]*16.043030; /*CH4 */
    XW += x[14]*28.010550; /*CO */
    XW += x[15]*44.009950; /*CO2 */
    XW += x[16]*29.018520; /*HCO */
    XW += x[17]*30.026490; /*CH2O */
    XW += x[18]*31.034460; /*CH2OH */
    XW += x[19]*31.034460; /*CH3O */
    XW += x[20]*32.042430; /*CH3OH */
    XW += x[21]*25.030270; /*C2H */
    XW += x[22]*26.038240; /*C2H2 */
    XW += x[23]*27.046210; /*C2H3 */
    XW += x[24]*28.054180; /*C2H4 */
    XW += x[25]*29.062150; /*C2H5 */
    XW += x[26]*30.070120; /*C2H6 */
    XW += x[27]*41.029670; /*HCCO */
    XW += x[28]*42.037640; /*CH2CO */
    XW += x[29]*42.037640; /*HCCOH */
    XW += x[30]*28.013400; /*N2 */
    XW += x[31]*39.948000; /*AR */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 32; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 177; ++id) {
        qdot[id] *= 1.0e-6;
    }
}

/*compute the reaction Jacobian on CPU */
void aJacobian(amrex::Real *  J, amrex::Real *  sc, amrex::Real T, int consP)
{
    for (int i=0; i<1089; i++) {
        J[i] = 0.0;
    }

    amrex::Real wdot[32];
    for (int k=0; k<32; k++) {
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
    for (int k = 0; k < 32; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    amrex::Real g_RT[32];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    amrex::Real h_RT[32];
    speciesEnthalpy(h_RT, tc);

    amrex::Real phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    amrex::Real dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    amrex::Real dqdci, dcdc_fac, dqdc[32];
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
    alpha = mixture + (TB[0][0] - 1)*sc[0] + (TB[0][1] - 1)*sc[5] + (TB[0][2] - 1)*sc[13] + (TB[0][3] - 1)*sc[14] + (TB[0][4] - 1)*sc[15] + (TB[0][5] - 1)*sc[26] + (TB[0][6] - 1)*sc[31];
    /* forward */
    phi_f = sc[1]*sc[10];
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
    phi_r = sc[12];
    Kc = refCinv * exp(g_RT[1] + g_RT[10] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[10]) + (h_RT[12]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[10] -= q; /* CH2 */
    wdot[12] += q; /* CH3 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[0][0] - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[10] -= dqdci;               /* dwdot[CH2]/d[H2] */
        J[12] += dqdci;               /* dwdot[CH3]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[10];
        J[34] -= dqdci;               /* dwdot[H]/d[H] */
        J[43] -= dqdci;               /* dwdot[CH2]/d[H] */
        J[45] += dqdci;               /* dwdot[CH3]/d[H] */
        /* d()/d[H2O] */
        dqdci = (TB[0][1] - 1)*dcdc_fac;
        J[166] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[175] -= dqdci;              /* dwdot[CH2]/d[H2O] */
        J[177] += dqdci;              /* dwdot[CH3]/d[H2O] */
        /* d()/d[CH2] */
        dqdci =  + k_f*sc[1];
        J[331] -= dqdci;              /* dwdot[H]/d[CH2] */
        J[340] -= dqdci;              /* dwdot[CH2]/d[CH2] */
        J[342] += dqdci;              /* dwdot[CH3]/d[CH2] */
        /* d()/d[CH3] */
        dqdci =  - k_r;
        J[397] -= dqdci;              /* dwdot[H]/d[CH3] */
        J[406] -= dqdci;              /* dwdot[CH2]/d[CH3] */
        J[408] += dqdci;              /* dwdot[CH3]/d[CH3] */
        /* d()/d[CH4] */
        dqdci = (TB[0][2] - 1)*dcdc_fac;
        J[430] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[439] -= dqdci;              /* dwdot[CH2]/d[CH4] */
        J[441] += dqdci;              /* dwdot[CH3]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[0][3] - 1)*dcdc_fac;
        J[463] -= dqdci;              /* dwdot[H]/d[CO] */
        J[472] -= dqdci;              /* dwdot[CH2]/d[CO] */
        J[474] += dqdci;              /* dwdot[CH3]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[0][4] - 1)*dcdc_fac;
        J[496] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[505] -= dqdci;              /* dwdot[CH2]/d[CO2] */
        J[507] += dqdci;              /* dwdot[CH3]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (TB[0][5] - 1)*dcdc_fac;
        J[859] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[868] -= dqdci;              /* dwdot[CH2]/d[C2H6] */
        J[870] += dqdci;              /* dwdot[CH3]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (TB[0][6] - 1)*dcdc_fac;
        J[1024] -= dqdci;             /* dwdot[H]/d[AR] */
        J[1033] -= dqdci;             /* dwdot[CH2]/d[AR] */
        J[1035] += dqdci;             /* dwdot[CH3]/d[AR] */
    }
    else {
        dqdc[0] = TB[0][0]*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[10];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = TB[0][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac + k_f*sc[1];
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac - k_r;
        dqdc[13] = TB[0][2]*dcdc_fac;
        dqdc[14] = TB[0][3]*dcdc_fac;
        dqdc[15] = TB[0][4]*dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac;
        dqdc[22] = dcdc_fac;
        dqdc[23] = dcdc_fac;
        dqdc[24] = dcdc_fac;
        dqdc[25] = dcdc_fac;
        dqdc[26] = TB[0][5]*dcdc_fac;
        dqdc[27] = dcdc_fac;
        dqdc[28] = dcdc_fac;
        dqdc[29] = dcdc_fac;
        dqdc[30] = dcdc_fac;
        dqdc[31] = TB[0][6]*dcdc_fac;
        for (int k=0; k<32; k++) {
            J[33*k+1] -= dqdc[k];
            J[33*k+10] -= dqdc[k];
            J[33*k+12] += dqdc[k];
        }
    }
    J[1057] -= dqdT; /* dwdot[H]/dT */
    J[1066] -= dqdT; /* dwdot[CH2]/dT */
    J[1068] += dqdT; /* dwdot[CH3]/dT */

    /*reaction 2: H + CH3 (+M) <=> CH4 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[1][0] - 1)*sc[0] + (TB[1][1] - 1)*sc[5] + (TB[1][2] - 1)*sc[13] + (TB[1][3] - 1)*sc[14] + (TB[1][4] - 1)*sc[15] + (TB[1][5] - 1)*sc[26] + (TB[1][6] - 1)*sc[31];
    /* forward */
    phi_f = sc[1]*sc[12];
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
    Kc = refCinv * exp(g_RT[1] + g_RT[12] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[12]) + (h_RT[13]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[12] -= q; /* CH3 */
    wdot[13] += q; /* CH4 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[1][0] - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[12] -= dqdci;               /* dwdot[CH3]/d[H2] */
        J[13] += dqdci;               /* dwdot[CH4]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[12];
        J[34] -= dqdci;               /* dwdot[H]/d[H] */
        J[45] -= dqdci;               /* dwdot[CH3]/d[H] */
        J[46] += dqdci;               /* dwdot[CH4]/d[H] */
        /* d()/d[H2O] */
        dqdci = (TB[1][1] - 1)*dcdc_fac;
        J[166] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[177] -= dqdci;              /* dwdot[CH3]/d[H2O] */
        J[178] += dqdci;              /* dwdot[CH4]/d[H2O] */
        /* d()/d[CH3] */
        dqdci =  + k_f*sc[1];
        J[397] -= dqdci;              /* dwdot[H]/d[CH3] */
        J[408] -= dqdci;              /* dwdot[CH3]/d[CH3] */
        J[409] += dqdci;              /* dwdot[CH4]/d[CH3] */
        /* d()/d[CH4] */
        dqdci = (TB[1][2] - 1)*dcdc_fac - k_r;
        J[430] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[441] -= dqdci;              /* dwdot[CH3]/d[CH4] */
        J[442] += dqdci;              /* dwdot[CH4]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[1][3] - 1)*dcdc_fac;
        J[463] -= dqdci;              /* dwdot[H]/d[CO] */
        J[474] -= dqdci;              /* dwdot[CH3]/d[CO] */
        J[475] += dqdci;              /* dwdot[CH4]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[1][4] - 1)*dcdc_fac;
        J[496] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[507] -= dqdci;              /* dwdot[CH3]/d[CO2] */
        J[508] += dqdci;              /* dwdot[CH4]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (TB[1][5] - 1)*dcdc_fac;
        J[859] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[870] -= dqdci;              /* dwdot[CH3]/d[C2H6] */
        J[871] += dqdci;              /* dwdot[CH4]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (TB[1][6] - 1)*dcdc_fac;
        J[1024] -= dqdci;             /* dwdot[H]/d[AR] */
        J[1035] -= dqdci;             /* dwdot[CH3]/d[AR] */
        J[1036] += dqdci;             /* dwdot[CH4]/d[AR] */
    }
    else {
        dqdc[0] = TB[1][0]*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[12];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = TB[1][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac + k_f*sc[1];
        dqdc[13] = TB[1][2]*dcdc_fac - k_r;
        dqdc[14] = TB[1][3]*dcdc_fac;
        dqdc[15] = TB[1][4]*dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac;
        dqdc[22] = dcdc_fac;
        dqdc[23] = dcdc_fac;
        dqdc[24] = dcdc_fac;
        dqdc[25] = dcdc_fac;
        dqdc[26] = TB[1][5]*dcdc_fac;
        dqdc[27] = dcdc_fac;
        dqdc[28] = dcdc_fac;
        dqdc[29] = dcdc_fac;
        dqdc[30] = dcdc_fac;
        dqdc[31] = TB[1][6]*dcdc_fac;
        for (int k=0; k<32; k++) {
            J[33*k+1] -= dqdc[k];
            J[33*k+12] -= dqdc[k];
            J[33*k+13] += dqdc[k];
        }
    }
    J[1057] -= dqdT; /* dwdot[H]/dT */
    J[1068] -= dqdT; /* dwdot[CH3]/dT */
    J[1069] += dqdT; /* dwdot[CH4]/dT */

    /*reaction 3: H + HCO (+M) <=> CH2O (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[2][0] - 1)*sc[0] + (TB[2][1] - 1)*sc[5] + (TB[2][2] - 1)*sc[13] + (TB[2][3] - 1)*sc[14] + (TB[2][4] - 1)*sc[15] + (TB[2][5] - 1)*sc[26] + (TB[2][6] - 1)*sc[31];
    /* forward */
    phi_f = sc[1]*sc[16];
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
    wdot[16] -= q; /* HCO */
    wdot[17] += q; /* CH2O */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[2][0] - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[16] -= dqdci;               /* dwdot[HCO]/d[H2] */
        J[17] += dqdci;               /* dwdot[CH2O]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[16];
        J[34] -= dqdci;               /* dwdot[H]/d[H] */
        J[49] -= dqdci;               /* dwdot[HCO]/d[H] */
        J[50] += dqdci;               /* dwdot[CH2O]/d[H] */
        /* d()/d[H2O] */
        dqdci = (TB[2][1] - 1)*dcdc_fac;
        J[166] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[181] -= dqdci;              /* dwdot[HCO]/d[H2O] */
        J[182] += dqdci;              /* dwdot[CH2O]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (TB[2][2] - 1)*dcdc_fac;
        J[430] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[445] -= dqdci;              /* dwdot[HCO]/d[CH4] */
        J[446] += dqdci;              /* dwdot[CH2O]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[2][3] - 1)*dcdc_fac;
        J[463] -= dqdci;              /* dwdot[H]/d[CO] */
        J[478] -= dqdci;              /* dwdot[HCO]/d[CO] */
        J[479] += dqdci;              /* dwdot[CH2O]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[2][4] - 1)*dcdc_fac;
        J[496] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[511] -= dqdci;              /* dwdot[HCO]/d[CO2] */
        J[512] += dqdci;              /* dwdot[CH2O]/d[CO2] */
        /* d()/d[HCO] */
        dqdci =  + k_f*sc[1];
        J[529] -= dqdci;              /* dwdot[H]/d[HCO] */
        J[544] -= dqdci;              /* dwdot[HCO]/d[HCO] */
        J[545] += dqdci;              /* dwdot[CH2O]/d[HCO] */
        /* d()/d[CH2O] */
        dqdci =  - k_r;
        J[562] -= dqdci;              /* dwdot[H]/d[CH2O] */
        J[577] -= dqdci;              /* dwdot[HCO]/d[CH2O] */
        J[578] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
        /* d()/d[C2H6] */
        dqdci = (TB[2][5] - 1)*dcdc_fac;
        J[859] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[874] -= dqdci;              /* dwdot[HCO]/d[C2H6] */
        J[875] += dqdci;              /* dwdot[CH2O]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (TB[2][6] - 1)*dcdc_fac;
        J[1024] -= dqdci;             /* dwdot[H]/d[AR] */
        J[1039] -= dqdci;             /* dwdot[HCO]/d[AR] */
        J[1040] += dqdci;             /* dwdot[CH2O]/d[AR] */
    }
    else {
        dqdc[0] = TB[2][0]*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[16];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = TB[2][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = TB[2][2]*dcdc_fac;
        dqdc[14] = TB[2][3]*dcdc_fac;
        dqdc[15] = TB[2][4]*dcdc_fac;
        dqdc[16] = dcdc_fac + k_f*sc[1];
        dqdc[17] = dcdc_fac - k_r;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac;
        dqdc[22] = dcdc_fac;
        dqdc[23] = dcdc_fac;
        dqdc[24] = dcdc_fac;
        dqdc[25] = dcdc_fac;
        dqdc[26] = TB[2][5]*dcdc_fac;
        dqdc[27] = dcdc_fac;
        dqdc[28] = dcdc_fac;
        dqdc[29] = dcdc_fac;
        dqdc[30] = dcdc_fac;
        dqdc[31] = TB[2][6]*dcdc_fac;
        for (int k=0; k<32; k++) {
            J[33*k+1] -= dqdc[k];
            J[33*k+16] -= dqdc[k];
            J[33*k+17] += dqdc[k];
        }
    }
    J[1057] -= dqdT; /* dwdot[H]/dT */
    J[1072] -= dqdT; /* dwdot[HCO]/dT */
    J[1073] += dqdT; /* dwdot[CH2O]/dT */

    /*reaction 4: H + CH2O (+M) <=> CH2OH (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[3][0] - 1)*sc[0] + (TB[3][1] - 1)*sc[5] + (TB[3][2] - 1)*sc[13] + (TB[3][3] - 1)*sc[14] + (TB[3][4] - 1)*sc[15] + (TB[3][5] - 1)*sc[26];
    /* forward */
    phi_f = sc[1]*sc[17];
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
    wdot[17] -= q; /* CH2O */
    wdot[18] += q; /* CH2OH */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[3][0] - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[17] -= dqdci;               /* dwdot[CH2O]/d[H2] */
        J[18] += dqdci;               /* dwdot[CH2OH]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[17];
        J[34] -= dqdci;               /* dwdot[H]/d[H] */
        J[50] -= dqdci;               /* dwdot[CH2O]/d[H] */
        J[51] += dqdci;               /* dwdot[CH2OH]/d[H] */
        /* d()/d[H2O] */
        dqdci = (TB[3][1] - 1)*dcdc_fac;
        J[166] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[182] -= dqdci;              /* dwdot[CH2O]/d[H2O] */
        J[183] += dqdci;              /* dwdot[CH2OH]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (TB[3][2] - 1)*dcdc_fac;
        J[430] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[446] -= dqdci;              /* dwdot[CH2O]/d[CH4] */
        J[447] += dqdci;              /* dwdot[CH2OH]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[3][3] - 1)*dcdc_fac;
        J[463] -= dqdci;              /* dwdot[H]/d[CO] */
        J[479] -= dqdci;              /* dwdot[CH2O]/d[CO] */
        J[480] += dqdci;              /* dwdot[CH2OH]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[3][4] - 1)*dcdc_fac;
        J[496] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[512] -= dqdci;              /* dwdot[CH2O]/d[CO2] */
        J[513] += dqdci;              /* dwdot[CH2OH]/d[CO2] */
        /* d()/d[CH2O] */
        dqdci =  + k_f*sc[1];
        J[562] -= dqdci;              /* dwdot[H]/d[CH2O] */
        J[578] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
        J[579] += dqdci;              /* dwdot[CH2OH]/d[CH2O] */
        /* d()/d[CH2OH] */
        dqdci =  - k_r;
        J[595] -= dqdci;              /* dwdot[H]/d[CH2OH] */
        J[611] -= dqdci;              /* dwdot[CH2O]/d[CH2OH] */
        J[612] += dqdci;              /* dwdot[CH2OH]/d[CH2OH] */
        /* d()/d[C2H6] */
        dqdci = (TB[3][5] - 1)*dcdc_fac;
        J[859] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[875] -= dqdci;              /* dwdot[CH2O]/d[C2H6] */
        J[876] += dqdci;              /* dwdot[CH2OH]/d[C2H6] */
    }
    else {
        dqdc[0] = TB[3][0]*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[17];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = TB[3][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = TB[3][2]*dcdc_fac;
        dqdc[14] = TB[3][3]*dcdc_fac;
        dqdc[15] = TB[3][4]*dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac + k_f*sc[1];
        dqdc[18] = dcdc_fac - k_r;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac;
        dqdc[22] = dcdc_fac;
        dqdc[23] = dcdc_fac;
        dqdc[24] = dcdc_fac;
        dqdc[25] = dcdc_fac;
        dqdc[26] = TB[3][5]*dcdc_fac;
        dqdc[27] = dcdc_fac;
        dqdc[28] = dcdc_fac;
        dqdc[29] = dcdc_fac;
        dqdc[30] = dcdc_fac;
        dqdc[31] = dcdc_fac;
        for (int k=0; k<32; k++) {
            J[33*k+1] -= dqdc[k];
            J[33*k+17] -= dqdc[k];
            J[33*k+18] += dqdc[k];
        }
    }
    J[1057] -= dqdT; /* dwdot[H]/dT */
    J[1073] -= dqdT; /* dwdot[CH2O]/dT */
    J[1074] += dqdT; /* dwdot[CH2OH]/dT */

    /*reaction 5: H + CH2O (+M) <=> CH3O (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[4][0] - 1)*sc[0] + (TB[4][1] - 1)*sc[5] + (TB[4][2] - 1)*sc[13] + (TB[4][3] - 1)*sc[14] + (TB[4][4] - 1)*sc[15] + (TB[4][5] - 1)*sc[26];
    /* forward */
    phi_f = sc[1]*sc[17];
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
    phi_r = sc[19];
    Kc = refCinv * exp(g_RT[1] + g_RT[17] - g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[17]) + (h_RT[19]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[17] -= q; /* CH2O */
    wdot[19] += q; /* CH3O */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[4][0] - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[17] -= dqdci;               /* dwdot[CH2O]/d[H2] */
        J[19] += dqdci;               /* dwdot[CH3O]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[17];
        J[34] -= dqdci;               /* dwdot[H]/d[H] */
        J[50] -= dqdci;               /* dwdot[CH2O]/d[H] */
        J[52] += dqdci;               /* dwdot[CH3O]/d[H] */
        /* d()/d[H2O] */
        dqdci = (TB[4][1] - 1)*dcdc_fac;
        J[166] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[182] -= dqdci;              /* dwdot[CH2O]/d[H2O] */
        J[184] += dqdci;              /* dwdot[CH3O]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (TB[4][2] - 1)*dcdc_fac;
        J[430] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[446] -= dqdci;              /* dwdot[CH2O]/d[CH4] */
        J[448] += dqdci;              /* dwdot[CH3O]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[4][3] - 1)*dcdc_fac;
        J[463] -= dqdci;              /* dwdot[H]/d[CO] */
        J[479] -= dqdci;              /* dwdot[CH2O]/d[CO] */
        J[481] += dqdci;              /* dwdot[CH3O]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[4][4] - 1)*dcdc_fac;
        J[496] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[512] -= dqdci;              /* dwdot[CH2O]/d[CO2] */
        J[514] += dqdci;              /* dwdot[CH3O]/d[CO2] */
        /* d()/d[CH2O] */
        dqdci =  + k_f*sc[1];
        J[562] -= dqdci;              /* dwdot[H]/d[CH2O] */
        J[578] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
        J[580] += dqdci;              /* dwdot[CH3O]/d[CH2O] */
        /* d()/d[CH3O] */
        dqdci =  - k_r;
        J[628] -= dqdci;              /* dwdot[H]/d[CH3O] */
        J[644] -= dqdci;              /* dwdot[CH2O]/d[CH3O] */
        J[646] += dqdci;              /* dwdot[CH3O]/d[CH3O] */
        /* d()/d[C2H6] */
        dqdci = (TB[4][5] - 1)*dcdc_fac;
        J[859] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[875] -= dqdci;              /* dwdot[CH2O]/d[C2H6] */
        J[877] += dqdci;              /* dwdot[CH3O]/d[C2H6] */
    }
    else {
        dqdc[0] = TB[4][0]*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[17];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = TB[4][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = TB[4][2]*dcdc_fac;
        dqdc[14] = TB[4][3]*dcdc_fac;
        dqdc[15] = TB[4][4]*dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac + k_f*sc[1];
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac - k_r;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac;
        dqdc[22] = dcdc_fac;
        dqdc[23] = dcdc_fac;
        dqdc[24] = dcdc_fac;
        dqdc[25] = dcdc_fac;
        dqdc[26] = TB[4][5]*dcdc_fac;
        dqdc[27] = dcdc_fac;
        dqdc[28] = dcdc_fac;
        dqdc[29] = dcdc_fac;
        dqdc[30] = dcdc_fac;
        dqdc[31] = dcdc_fac;
        for (int k=0; k<32; k++) {
            J[33*k+1] -= dqdc[k];
            J[33*k+17] -= dqdc[k];
            J[33*k+19] += dqdc[k];
        }
    }
    J[1057] -= dqdT; /* dwdot[H]/dT */
    J[1073] -= dqdT; /* dwdot[CH2O]/dT */
    J[1075] += dqdT; /* dwdot[CH3O]/dT */

    /*reaction 6: H + CH2OH (+M) <=> CH3OH (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[5][0] - 1)*sc[0] + (TB[5][1] - 1)*sc[5] + (TB[5][2] - 1)*sc[13] + (TB[5][3] - 1)*sc[14] + (TB[5][4] - 1)*sc[15] + (TB[5][5] - 1)*sc[26];
    /* forward */
    phi_f = sc[1]*sc[18];
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
    phi_r = sc[20];
    Kc = refCinv * exp(g_RT[1] + g_RT[18] - g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[18]) + (h_RT[20]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[18] -= q; /* CH2OH */
    wdot[20] += q; /* CH3OH */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[5][0] - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[18] -= dqdci;               /* dwdot[CH2OH]/d[H2] */
        J[20] += dqdci;               /* dwdot[CH3OH]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[18];
        J[34] -= dqdci;               /* dwdot[H]/d[H] */
        J[51] -= dqdci;               /* dwdot[CH2OH]/d[H] */
        J[53] += dqdci;               /* dwdot[CH3OH]/d[H] */
        /* d()/d[H2O] */
        dqdci = (TB[5][1] - 1)*dcdc_fac;
        J[166] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[183] -= dqdci;              /* dwdot[CH2OH]/d[H2O] */
        J[185] += dqdci;              /* dwdot[CH3OH]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (TB[5][2] - 1)*dcdc_fac;
        J[430] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[447] -= dqdci;              /* dwdot[CH2OH]/d[CH4] */
        J[449] += dqdci;              /* dwdot[CH3OH]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[5][3] - 1)*dcdc_fac;
        J[463] -= dqdci;              /* dwdot[H]/d[CO] */
        J[480] -= dqdci;              /* dwdot[CH2OH]/d[CO] */
        J[482] += dqdci;              /* dwdot[CH3OH]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[5][4] - 1)*dcdc_fac;
        J[496] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[513] -= dqdci;              /* dwdot[CH2OH]/d[CO2] */
        J[515] += dqdci;              /* dwdot[CH3OH]/d[CO2] */
        /* d()/d[CH2OH] */
        dqdci =  + k_f*sc[1];
        J[595] -= dqdci;              /* dwdot[H]/d[CH2OH] */
        J[612] -= dqdci;              /* dwdot[CH2OH]/d[CH2OH] */
        J[614] += dqdci;              /* dwdot[CH3OH]/d[CH2OH] */
        /* d()/d[CH3OH] */
        dqdci =  - k_r;
        J[661] -= dqdci;              /* dwdot[H]/d[CH3OH] */
        J[678] -= dqdci;              /* dwdot[CH2OH]/d[CH3OH] */
        J[680] += dqdci;              /* dwdot[CH3OH]/d[CH3OH] */
        /* d()/d[C2H6] */
        dqdci = (TB[5][5] - 1)*dcdc_fac;
        J[859] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[876] -= dqdci;              /* dwdot[CH2OH]/d[C2H6] */
        J[878] += dqdci;              /* dwdot[CH3OH]/d[C2H6] */
    }
    else {
        dqdc[0] = TB[5][0]*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[18];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = TB[5][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = TB[5][2]*dcdc_fac;
        dqdc[14] = TB[5][3]*dcdc_fac;
        dqdc[15] = TB[5][4]*dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac + k_f*sc[1];
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac - k_r;
        dqdc[21] = dcdc_fac;
        dqdc[22] = dcdc_fac;
        dqdc[23] = dcdc_fac;
        dqdc[24] = dcdc_fac;
        dqdc[25] = dcdc_fac;
        dqdc[26] = TB[5][5]*dcdc_fac;
        dqdc[27] = dcdc_fac;
        dqdc[28] = dcdc_fac;
        dqdc[29] = dcdc_fac;
        dqdc[30] = dcdc_fac;
        dqdc[31] = dcdc_fac;
        for (int k=0; k<32; k++) {
            J[33*k+1] -= dqdc[k];
            J[33*k+18] -= dqdc[k];
            J[33*k+20] += dqdc[k];
        }
    }
    J[1057] -= dqdT; /* dwdot[H]/dT */
    J[1074] -= dqdT; /* dwdot[CH2OH]/dT */
    J[1076] += dqdT; /* dwdot[CH3OH]/dT */

    /*reaction 7: H + CH3O (+M) <=> CH3OH (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[6][0] - 1)*sc[0] + (TB[6][1] - 1)*sc[5] + (TB[6][2] - 1)*sc[13] + (TB[6][3] - 1)*sc[14] + (TB[6][4] - 1)*sc[15] + (TB[6][5] - 1)*sc[26];
    /* forward */
    phi_f = sc[1]*sc[19];
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
    phi_r = sc[20];
    Kc = refCinv * exp(g_RT[1] + g_RT[19] - g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[19]) + (h_RT[20]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[19] -= q; /* CH3O */
    wdot[20] += q; /* CH3OH */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[6][0] - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[19] -= dqdci;               /* dwdot[CH3O]/d[H2] */
        J[20] += dqdci;               /* dwdot[CH3OH]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[19];
        J[34] -= dqdci;               /* dwdot[H]/d[H] */
        J[52] -= dqdci;               /* dwdot[CH3O]/d[H] */
        J[53] += dqdci;               /* dwdot[CH3OH]/d[H] */
        /* d()/d[H2O] */
        dqdci = (TB[6][1] - 1)*dcdc_fac;
        J[166] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[184] -= dqdci;              /* dwdot[CH3O]/d[H2O] */
        J[185] += dqdci;              /* dwdot[CH3OH]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (TB[6][2] - 1)*dcdc_fac;
        J[430] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[448] -= dqdci;              /* dwdot[CH3O]/d[CH4] */
        J[449] += dqdci;              /* dwdot[CH3OH]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[6][3] - 1)*dcdc_fac;
        J[463] -= dqdci;              /* dwdot[H]/d[CO] */
        J[481] -= dqdci;              /* dwdot[CH3O]/d[CO] */
        J[482] += dqdci;              /* dwdot[CH3OH]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[6][4] - 1)*dcdc_fac;
        J[496] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[514] -= dqdci;              /* dwdot[CH3O]/d[CO2] */
        J[515] += dqdci;              /* dwdot[CH3OH]/d[CO2] */
        /* d()/d[CH3O] */
        dqdci =  + k_f*sc[1];
        J[628] -= dqdci;              /* dwdot[H]/d[CH3O] */
        J[646] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
        J[647] += dqdci;              /* dwdot[CH3OH]/d[CH3O] */
        /* d()/d[CH3OH] */
        dqdci =  - k_r;
        J[661] -= dqdci;              /* dwdot[H]/d[CH3OH] */
        J[679] -= dqdci;              /* dwdot[CH3O]/d[CH3OH] */
        J[680] += dqdci;              /* dwdot[CH3OH]/d[CH3OH] */
        /* d()/d[C2H6] */
        dqdci = (TB[6][5] - 1)*dcdc_fac;
        J[859] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[877] -= dqdci;              /* dwdot[CH3O]/d[C2H6] */
        J[878] += dqdci;              /* dwdot[CH3OH]/d[C2H6] */
    }
    else {
        dqdc[0] = TB[6][0]*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[19];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = TB[6][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = TB[6][2]*dcdc_fac;
        dqdc[14] = TB[6][3]*dcdc_fac;
        dqdc[15] = TB[6][4]*dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac + k_f*sc[1];
        dqdc[20] = dcdc_fac - k_r;
        dqdc[21] = dcdc_fac;
        dqdc[22] = dcdc_fac;
        dqdc[23] = dcdc_fac;
        dqdc[24] = dcdc_fac;
        dqdc[25] = dcdc_fac;
        dqdc[26] = TB[6][5]*dcdc_fac;
        dqdc[27] = dcdc_fac;
        dqdc[28] = dcdc_fac;
        dqdc[29] = dcdc_fac;
        dqdc[30] = dcdc_fac;
        dqdc[31] = dcdc_fac;
        for (int k=0; k<32; k++) {
            J[33*k+1] -= dqdc[k];
            J[33*k+19] -= dqdc[k];
            J[33*k+20] += dqdc[k];
        }
    }
    J[1057] -= dqdT; /* dwdot[H]/dT */
    J[1075] -= dqdT; /* dwdot[CH3O]/dT */
    J[1076] += dqdT; /* dwdot[CH3OH]/dT */

    /*reaction 8: H + C2H (+M) <=> C2H2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[7][0] - 1)*sc[0] + (TB[7][1] - 1)*sc[5] + (TB[7][2] - 1)*sc[13] + (TB[7][3] - 1)*sc[14] + (TB[7][4] - 1)*sc[15] + (TB[7][5] - 1)*sc[26] + (TB[7][6] - 1)*sc[31];
    /* forward */
    phi_f = sc[1]*sc[21];
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
    phi_r = sc[22];
    Kc = refCinv * exp(g_RT[1] + g_RT[21] - g_RT[22]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[21]) + (h_RT[22]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[21] -= q; /* C2H */
    wdot[22] += q; /* C2H2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[7][0] - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[21] -= dqdci;               /* dwdot[C2H]/d[H2] */
        J[22] += dqdci;               /* dwdot[C2H2]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[21];
        J[34] -= dqdci;               /* dwdot[H]/d[H] */
        J[54] -= dqdci;               /* dwdot[C2H]/d[H] */
        J[55] += dqdci;               /* dwdot[C2H2]/d[H] */
        /* d()/d[H2O] */
        dqdci = (TB[7][1] - 1)*dcdc_fac;
        J[166] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[186] -= dqdci;              /* dwdot[C2H]/d[H2O] */
        J[187] += dqdci;              /* dwdot[C2H2]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (TB[7][2] - 1)*dcdc_fac;
        J[430] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[450] -= dqdci;              /* dwdot[C2H]/d[CH4] */
        J[451] += dqdci;              /* dwdot[C2H2]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[7][3] - 1)*dcdc_fac;
        J[463] -= dqdci;              /* dwdot[H]/d[CO] */
        J[483] -= dqdci;              /* dwdot[C2H]/d[CO] */
        J[484] += dqdci;              /* dwdot[C2H2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[7][4] - 1)*dcdc_fac;
        J[496] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[516] -= dqdci;              /* dwdot[C2H]/d[CO2] */
        J[517] += dqdci;              /* dwdot[C2H2]/d[CO2] */
        /* d()/d[C2H] */
        dqdci =  + k_f*sc[1];
        J[694] -= dqdci;              /* dwdot[H]/d[C2H] */
        J[714] -= dqdci;              /* dwdot[C2H]/d[C2H] */
        J[715] += dqdci;              /* dwdot[C2H2]/d[C2H] */
        /* d()/d[C2H2] */
        dqdci =  - k_r;
        J[727] -= dqdci;              /* dwdot[H]/d[C2H2] */
        J[747] -= dqdci;              /* dwdot[C2H]/d[C2H2] */
        J[748] += dqdci;              /* dwdot[C2H2]/d[C2H2] */
        /* d()/d[C2H6] */
        dqdci = (TB[7][5] - 1)*dcdc_fac;
        J[859] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[879] -= dqdci;              /* dwdot[C2H]/d[C2H6] */
        J[880] += dqdci;              /* dwdot[C2H2]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (TB[7][6] - 1)*dcdc_fac;
        J[1024] -= dqdci;             /* dwdot[H]/d[AR] */
        J[1044] -= dqdci;             /* dwdot[C2H]/d[AR] */
        J[1045] += dqdci;             /* dwdot[C2H2]/d[AR] */
    }
    else {
        dqdc[0] = TB[7][0]*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[21];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = TB[7][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = TB[7][2]*dcdc_fac;
        dqdc[14] = TB[7][3]*dcdc_fac;
        dqdc[15] = TB[7][4]*dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac + k_f*sc[1];
        dqdc[22] = dcdc_fac - k_r;
        dqdc[23] = dcdc_fac;
        dqdc[24] = dcdc_fac;
        dqdc[25] = dcdc_fac;
        dqdc[26] = TB[7][5]*dcdc_fac;
        dqdc[27] = dcdc_fac;
        dqdc[28] = dcdc_fac;
        dqdc[29] = dcdc_fac;
        dqdc[30] = dcdc_fac;
        dqdc[31] = TB[7][6]*dcdc_fac;
        for (int k=0; k<32; k++) {
            J[33*k+1] -= dqdc[k];
            J[33*k+21] -= dqdc[k];
            J[33*k+22] += dqdc[k];
        }
    }
    J[1057] -= dqdT; /* dwdot[H]/dT */
    J[1077] -= dqdT; /* dwdot[C2H]/dT */
    J[1078] += dqdT; /* dwdot[C2H2]/dT */

    /*reaction 9: H + C2H2 (+M) <=> C2H3 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[8][0] - 1)*sc[0] + (TB[8][1] - 1)*sc[5] + (TB[8][2] - 1)*sc[13] + (TB[8][3] - 1)*sc[14] + (TB[8][4] - 1)*sc[15] + (TB[8][5] - 1)*sc[26] + (TB[8][6] - 1)*sc[31];
    /* forward */
    phi_f = sc[1]*sc[22];
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
    phi_r = sc[23];
    Kc = refCinv * exp(g_RT[1] + g_RT[22] - g_RT[23]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[22]) + (h_RT[23]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[22] -= q; /* C2H2 */
    wdot[23] += q; /* C2H3 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[8][0] - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[22] -= dqdci;               /* dwdot[C2H2]/d[H2] */
        J[23] += dqdci;               /* dwdot[C2H3]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[22];
        J[34] -= dqdci;               /* dwdot[H]/d[H] */
        J[55] -= dqdci;               /* dwdot[C2H2]/d[H] */
        J[56] += dqdci;               /* dwdot[C2H3]/d[H] */
        /* d()/d[H2O] */
        dqdci = (TB[8][1] - 1)*dcdc_fac;
        J[166] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[187] -= dqdci;              /* dwdot[C2H2]/d[H2O] */
        J[188] += dqdci;              /* dwdot[C2H3]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (TB[8][2] - 1)*dcdc_fac;
        J[430] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[451] -= dqdci;              /* dwdot[C2H2]/d[CH4] */
        J[452] += dqdci;              /* dwdot[C2H3]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[8][3] - 1)*dcdc_fac;
        J[463] -= dqdci;              /* dwdot[H]/d[CO] */
        J[484] -= dqdci;              /* dwdot[C2H2]/d[CO] */
        J[485] += dqdci;              /* dwdot[C2H3]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[8][4] - 1)*dcdc_fac;
        J[496] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[517] -= dqdci;              /* dwdot[C2H2]/d[CO2] */
        J[518] += dqdci;              /* dwdot[C2H3]/d[CO2] */
        /* d()/d[C2H2] */
        dqdci =  + k_f*sc[1];
        J[727] -= dqdci;              /* dwdot[H]/d[C2H2] */
        J[748] -= dqdci;              /* dwdot[C2H2]/d[C2H2] */
        J[749] += dqdci;              /* dwdot[C2H3]/d[C2H2] */
        /* d()/d[C2H3] */
        dqdci =  - k_r;
        J[760] -= dqdci;              /* dwdot[H]/d[C2H3] */
        J[781] -= dqdci;              /* dwdot[C2H2]/d[C2H3] */
        J[782] += dqdci;              /* dwdot[C2H3]/d[C2H3] */
        /* d()/d[C2H6] */
        dqdci = (TB[8][5] - 1)*dcdc_fac;
        J[859] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[880] -= dqdci;              /* dwdot[C2H2]/d[C2H6] */
        J[881] += dqdci;              /* dwdot[C2H3]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (TB[8][6] - 1)*dcdc_fac;
        J[1024] -= dqdci;             /* dwdot[H]/d[AR] */
        J[1045] -= dqdci;             /* dwdot[C2H2]/d[AR] */
        J[1046] += dqdci;             /* dwdot[C2H3]/d[AR] */
    }
    else {
        dqdc[0] = TB[8][0]*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[22];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = TB[8][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = TB[8][2]*dcdc_fac;
        dqdc[14] = TB[8][3]*dcdc_fac;
        dqdc[15] = TB[8][4]*dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac;
        dqdc[22] = dcdc_fac + k_f*sc[1];
        dqdc[23] = dcdc_fac - k_r;
        dqdc[24] = dcdc_fac;
        dqdc[25] = dcdc_fac;
        dqdc[26] = TB[8][5]*dcdc_fac;
        dqdc[27] = dcdc_fac;
        dqdc[28] = dcdc_fac;
        dqdc[29] = dcdc_fac;
        dqdc[30] = dcdc_fac;
        dqdc[31] = TB[8][6]*dcdc_fac;
        for (int k=0; k<32; k++) {
            J[33*k+1] -= dqdc[k];
            J[33*k+22] -= dqdc[k];
            J[33*k+23] += dqdc[k];
        }
    }
    J[1057] -= dqdT; /* dwdot[H]/dT */
    J[1078] -= dqdT; /* dwdot[C2H2]/dT */
    J[1079] += dqdT; /* dwdot[C2H3]/dT */

    /*reaction 10: H + C2H3 (+M) <=> C2H4 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[9][0] - 1)*sc[0] + (TB[9][1] - 1)*sc[5] + (TB[9][2] - 1)*sc[13] + (TB[9][3] - 1)*sc[14] + (TB[9][4] - 1)*sc[15] + (TB[9][5] - 1)*sc[26] + (TB[9][6] - 1)*sc[31];
    /* forward */
    phi_f = sc[1]*sc[23];
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
    phi_r = sc[24];
    Kc = refCinv * exp(g_RT[1] + g_RT[23] - g_RT[24]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[23]) + (h_RT[24]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[23] -= q; /* C2H3 */
    wdot[24] += q; /* C2H4 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[9][0] - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[23] -= dqdci;               /* dwdot[C2H3]/d[H2] */
        J[24] += dqdci;               /* dwdot[C2H4]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[23];
        J[34] -= dqdci;               /* dwdot[H]/d[H] */
        J[56] -= dqdci;               /* dwdot[C2H3]/d[H] */
        J[57] += dqdci;               /* dwdot[C2H4]/d[H] */
        /* d()/d[H2O] */
        dqdci = (TB[9][1] - 1)*dcdc_fac;
        J[166] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[188] -= dqdci;              /* dwdot[C2H3]/d[H2O] */
        J[189] += dqdci;              /* dwdot[C2H4]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (TB[9][2] - 1)*dcdc_fac;
        J[430] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[452] -= dqdci;              /* dwdot[C2H3]/d[CH4] */
        J[453] += dqdci;              /* dwdot[C2H4]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[9][3] - 1)*dcdc_fac;
        J[463] -= dqdci;              /* dwdot[H]/d[CO] */
        J[485] -= dqdci;              /* dwdot[C2H3]/d[CO] */
        J[486] += dqdci;              /* dwdot[C2H4]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[9][4] - 1)*dcdc_fac;
        J[496] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[518] -= dqdci;              /* dwdot[C2H3]/d[CO2] */
        J[519] += dqdci;              /* dwdot[C2H4]/d[CO2] */
        /* d()/d[C2H3] */
        dqdci =  + k_f*sc[1];
        J[760] -= dqdci;              /* dwdot[H]/d[C2H3] */
        J[782] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
        J[783] += dqdci;              /* dwdot[C2H4]/d[C2H3] */
        /* d()/d[C2H4] */
        dqdci =  - k_r;
        J[793] -= dqdci;              /* dwdot[H]/d[C2H4] */
        J[815] -= dqdci;              /* dwdot[C2H3]/d[C2H4] */
        J[816] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
        /* d()/d[C2H6] */
        dqdci = (TB[9][5] - 1)*dcdc_fac;
        J[859] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[881] -= dqdci;              /* dwdot[C2H3]/d[C2H6] */
        J[882] += dqdci;              /* dwdot[C2H4]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (TB[9][6] - 1)*dcdc_fac;
        J[1024] -= dqdci;             /* dwdot[H]/d[AR] */
        J[1046] -= dqdci;             /* dwdot[C2H3]/d[AR] */
        J[1047] += dqdci;             /* dwdot[C2H4]/d[AR] */
    }
    else {
        dqdc[0] = TB[9][0]*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[23];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = TB[9][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = TB[9][2]*dcdc_fac;
        dqdc[14] = TB[9][3]*dcdc_fac;
        dqdc[15] = TB[9][4]*dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac;
        dqdc[22] = dcdc_fac;
        dqdc[23] = dcdc_fac + k_f*sc[1];
        dqdc[24] = dcdc_fac - k_r;
        dqdc[25] = dcdc_fac;
        dqdc[26] = TB[9][5]*dcdc_fac;
        dqdc[27] = dcdc_fac;
        dqdc[28] = dcdc_fac;
        dqdc[29] = dcdc_fac;
        dqdc[30] = dcdc_fac;
        dqdc[31] = TB[9][6]*dcdc_fac;
        for (int k=0; k<32; k++) {
            J[33*k+1] -= dqdc[k];
            J[33*k+23] -= dqdc[k];
            J[33*k+24] += dqdc[k];
        }
    }
    J[1057] -= dqdT; /* dwdot[H]/dT */
    J[1079] -= dqdT; /* dwdot[C2H3]/dT */
    J[1080] += dqdT; /* dwdot[C2H4]/dT */

    /*reaction 11: H + C2H4 (+M) <=> C2H5 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[10][0] - 1)*sc[0] + (TB[10][1] - 1)*sc[5] + (TB[10][2] - 1)*sc[13] + (TB[10][3] - 1)*sc[14] + (TB[10][4] - 1)*sc[15] + (TB[10][5] - 1)*sc[26] + (TB[10][6] - 1)*sc[31];
    /* forward */
    phi_f = sc[1]*sc[24];
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
    phi_r = sc[25];
    Kc = refCinv * exp(g_RT[1] + g_RT[24] - g_RT[25]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[24]) + (h_RT[25]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[24] -= q; /* C2H4 */
    wdot[25] += q; /* C2H5 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[10][0] - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[24] -= dqdci;               /* dwdot[C2H4]/d[H2] */
        J[25] += dqdci;               /* dwdot[C2H5]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[24];
        J[34] -= dqdci;               /* dwdot[H]/d[H] */
        J[57] -= dqdci;               /* dwdot[C2H4]/d[H] */
        J[58] += dqdci;               /* dwdot[C2H5]/d[H] */
        /* d()/d[H2O] */
        dqdci = (TB[10][1] - 1)*dcdc_fac;
        J[166] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[189] -= dqdci;              /* dwdot[C2H4]/d[H2O] */
        J[190] += dqdci;              /* dwdot[C2H5]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (TB[10][2] - 1)*dcdc_fac;
        J[430] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[453] -= dqdci;              /* dwdot[C2H4]/d[CH4] */
        J[454] += dqdci;              /* dwdot[C2H5]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[10][3] - 1)*dcdc_fac;
        J[463] -= dqdci;              /* dwdot[H]/d[CO] */
        J[486] -= dqdci;              /* dwdot[C2H4]/d[CO] */
        J[487] += dqdci;              /* dwdot[C2H5]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[10][4] - 1)*dcdc_fac;
        J[496] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[519] -= dqdci;              /* dwdot[C2H4]/d[CO2] */
        J[520] += dqdci;              /* dwdot[C2H5]/d[CO2] */
        /* d()/d[C2H4] */
        dqdci =  + k_f*sc[1];
        J[793] -= dqdci;              /* dwdot[H]/d[C2H4] */
        J[816] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
        J[817] += dqdci;              /* dwdot[C2H5]/d[C2H4] */
        /* d()/d[C2H5] */
        dqdci =  - k_r;
        J[826] -= dqdci;              /* dwdot[H]/d[C2H5] */
        J[849] -= dqdci;              /* dwdot[C2H4]/d[C2H5] */
        J[850] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
        /* d()/d[C2H6] */
        dqdci = (TB[10][5] - 1)*dcdc_fac;
        J[859] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[882] -= dqdci;              /* dwdot[C2H4]/d[C2H6] */
        J[883] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (TB[10][6] - 1)*dcdc_fac;
        J[1024] -= dqdci;             /* dwdot[H]/d[AR] */
        J[1047] -= dqdci;             /* dwdot[C2H4]/d[AR] */
        J[1048] += dqdci;             /* dwdot[C2H5]/d[AR] */
    }
    else {
        dqdc[0] = TB[10][0]*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[24];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = TB[10][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = TB[10][2]*dcdc_fac;
        dqdc[14] = TB[10][3]*dcdc_fac;
        dqdc[15] = TB[10][4]*dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac;
        dqdc[22] = dcdc_fac;
        dqdc[23] = dcdc_fac;
        dqdc[24] = dcdc_fac + k_f*sc[1];
        dqdc[25] = dcdc_fac - k_r;
        dqdc[26] = TB[10][5]*dcdc_fac;
        dqdc[27] = dcdc_fac;
        dqdc[28] = dcdc_fac;
        dqdc[29] = dcdc_fac;
        dqdc[30] = dcdc_fac;
        dqdc[31] = TB[10][6]*dcdc_fac;
        for (int k=0; k<32; k++) {
            J[33*k+1] -= dqdc[k];
            J[33*k+24] -= dqdc[k];
            J[33*k+25] += dqdc[k];
        }
    }
    J[1057] -= dqdT; /* dwdot[H]/dT */
    J[1080] -= dqdT; /* dwdot[C2H4]/dT */
    J[1081] += dqdT; /* dwdot[C2H5]/dT */

    /*reaction 12: H + C2H5 (+M) <=> C2H6 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[11][0] - 1)*sc[0] + (TB[11][1] - 1)*sc[5] + (TB[11][2] - 1)*sc[13] + (TB[11][3] - 1)*sc[14] + (TB[11][4] - 1)*sc[15] + (TB[11][5] - 1)*sc[26] + (TB[11][6] - 1)*sc[31];
    /* forward */
    phi_f = sc[1]*sc[25];
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
    phi_r = sc[26];
    Kc = refCinv * exp(g_RT[1] + g_RT[25] - g_RT[26]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[25]) + (h_RT[26]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[25] -= q; /* C2H5 */
    wdot[26] += q; /* C2H6 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[11][0] - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[25] -= dqdci;               /* dwdot[C2H5]/d[H2] */
        J[26] += dqdci;               /* dwdot[C2H6]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[25];
        J[34] -= dqdci;               /* dwdot[H]/d[H] */
        J[58] -= dqdci;               /* dwdot[C2H5]/d[H] */
        J[59] += dqdci;               /* dwdot[C2H6]/d[H] */
        /* d()/d[H2O] */
        dqdci = (TB[11][1] - 1)*dcdc_fac;
        J[166] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[190] -= dqdci;              /* dwdot[C2H5]/d[H2O] */
        J[191] += dqdci;              /* dwdot[C2H6]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (TB[11][2] - 1)*dcdc_fac;
        J[430] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[454] -= dqdci;              /* dwdot[C2H5]/d[CH4] */
        J[455] += dqdci;              /* dwdot[C2H6]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[11][3] - 1)*dcdc_fac;
        J[463] -= dqdci;              /* dwdot[H]/d[CO] */
        J[487] -= dqdci;              /* dwdot[C2H5]/d[CO] */
        J[488] += dqdci;              /* dwdot[C2H6]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[11][4] - 1)*dcdc_fac;
        J[496] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[520] -= dqdci;              /* dwdot[C2H5]/d[CO2] */
        J[521] += dqdci;              /* dwdot[C2H6]/d[CO2] */
        /* d()/d[C2H5] */
        dqdci =  + k_f*sc[1];
        J[826] -= dqdci;              /* dwdot[H]/d[C2H5] */
        J[850] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
        J[851] += dqdci;              /* dwdot[C2H6]/d[C2H5] */
        /* d()/d[C2H6] */
        dqdci = (TB[11][5] - 1)*dcdc_fac - k_r;
        J[859] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[883] -= dqdci;              /* dwdot[C2H5]/d[C2H6] */
        J[884] += dqdci;              /* dwdot[C2H6]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (TB[11][6] - 1)*dcdc_fac;
        J[1024] -= dqdci;             /* dwdot[H]/d[AR] */
        J[1048] -= dqdci;             /* dwdot[C2H5]/d[AR] */
        J[1049] += dqdci;             /* dwdot[C2H6]/d[AR] */
    }
    else {
        dqdc[0] = TB[11][0]*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[25];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = TB[11][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = TB[11][2]*dcdc_fac;
        dqdc[14] = TB[11][3]*dcdc_fac;
        dqdc[15] = TB[11][4]*dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac;
        dqdc[22] = dcdc_fac;
        dqdc[23] = dcdc_fac;
        dqdc[24] = dcdc_fac;
        dqdc[25] = dcdc_fac + k_f*sc[1];
        dqdc[26] = TB[11][5]*dcdc_fac - k_r;
        dqdc[27] = dcdc_fac;
        dqdc[28] = dcdc_fac;
        dqdc[29] = dcdc_fac;
        dqdc[30] = dcdc_fac;
        dqdc[31] = TB[11][6]*dcdc_fac;
        for (int k=0; k<32; k++) {
            J[33*k+1] -= dqdc[k];
            J[33*k+25] -= dqdc[k];
            J[33*k+26] += dqdc[k];
        }
    }
    J[1057] -= dqdT; /* dwdot[H]/dT */
    J[1081] -= dqdT; /* dwdot[C2H5]/dT */
    J[1082] += dqdT; /* dwdot[C2H6]/dT */

    /*reaction 13: H2 + CO (+M) <=> CH2O (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[12][0] - 1)*sc[0] + (TB[12][1] - 1)*sc[5] + (TB[12][2] - 1)*sc[13] + (TB[12][3] - 1)*sc[14] + (TB[12][4] - 1)*sc[15] + (TB[12][5] - 1)*sc[26] + (TB[12][6] - 1)*sc[31];
    /* forward */
    phi_f = sc[0]*sc[14];
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
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(troe_Tsss[12]) > 1.e-100 ? (1.-troe_a[12])*exp(-T/troe_Tsss[12]) : 0.);
    Fcent2 = (fabs(troe_Ts[12]) > 1.e-100 ? troe_a[12] * exp(-T/troe_Ts[12]) : 0.);
    Fcent3 = (troe_len[12] == 4 ? exp(-troe_Tss[12] * invT) : 0.);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        (fabs(troe_Tsss[12]) > 1.e-100 ? -Fcent1/troe_Tsss[12] : 0.)
      + (fabs(troe_Ts[12]) > 1.e-100 ? -Fcent2/troe_Ts[12] : 0.)
      + (troe_len[12] == 4 ? Fcent3*troe_Tss[12]*invT2 : 0.) );
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[17];
    Kc = refCinv * exp(g_RT[0] + g_RT[14] - g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[14]) + (h_RT[17]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[14] -= q; /* CO */
    wdot[17] += q; /* CH2O */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[12][0] - 1)*dcdc_fac + k_f*sc[14];
        J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
        J[14] -= dqdci;               /* dwdot[CO]/d[H2] */
        J[17] += dqdci;               /* dwdot[CH2O]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[12][1] - 1)*dcdc_fac;
        J[165] -= dqdci;              /* dwdot[H2]/d[H2O] */
        J[179] -= dqdci;              /* dwdot[CO]/d[H2O] */
        J[182] += dqdci;              /* dwdot[CH2O]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (TB[12][2] - 1)*dcdc_fac;
        J[429] -= dqdci;              /* dwdot[H2]/d[CH4] */
        J[443] -= dqdci;              /* dwdot[CO]/d[CH4] */
        J[446] += dqdci;              /* dwdot[CH2O]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[12][3] - 1)*dcdc_fac + k_f*sc[0];
        J[462] -= dqdci;              /* dwdot[H2]/d[CO] */
        J[476] -= dqdci;              /* dwdot[CO]/d[CO] */
        J[479] += dqdci;              /* dwdot[CH2O]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[12][4] - 1)*dcdc_fac;
        J[495] -= dqdci;              /* dwdot[H2]/d[CO2] */
        J[509] -= dqdci;              /* dwdot[CO]/d[CO2] */
        J[512] += dqdci;              /* dwdot[CH2O]/d[CO2] */
        /* d()/d[CH2O] */
        dqdci =  - k_r;
        J[561] -= dqdci;              /* dwdot[H2]/d[CH2O] */
        J[575] -= dqdci;              /* dwdot[CO]/d[CH2O] */
        J[578] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
        /* d()/d[C2H6] */
        dqdci = (TB[12][5] - 1)*dcdc_fac;
        J[858] -= dqdci;              /* dwdot[H2]/d[C2H6] */
        J[872] -= dqdci;              /* dwdot[CO]/d[C2H6] */
        J[875] += dqdci;              /* dwdot[CH2O]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (TB[12][6] - 1)*dcdc_fac;
        J[1023] -= dqdci;             /* dwdot[H2]/d[AR] */
        J[1037] -= dqdci;             /* dwdot[CO]/d[AR] */
        J[1040] += dqdci;             /* dwdot[CH2O]/d[AR] */
    }
    else {
        dqdc[0] = TB[12][0]*dcdc_fac + k_f*sc[14];
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = TB[12][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = TB[12][2]*dcdc_fac;
        dqdc[14] = TB[12][3]*dcdc_fac + k_f*sc[0];
        dqdc[15] = TB[12][4]*dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac - k_r;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac;
        dqdc[22] = dcdc_fac;
        dqdc[23] = dcdc_fac;
        dqdc[24] = dcdc_fac;
        dqdc[25] = dcdc_fac;
        dqdc[26] = TB[12][5]*dcdc_fac;
        dqdc[27] = dcdc_fac;
        dqdc[28] = dcdc_fac;
        dqdc[29] = dcdc_fac;
        dqdc[30] = dcdc_fac;
        dqdc[31] = TB[12][6]*dcdc_fac;
        for (int k=0; k<32; k++) {
            J[33*k+0] -= dqdc[k];
            J[33*k+14] -= dqdc[k];
            J[33*k+17] += dqdc[k];
        }
    }
    J[1056] -= dqdT; /* dwdot[H2]/dT */
    J[1070] -= dqdT; /* dwdot[CO]/dT */
    J[1073] += dqdT; /* dwdot[CH2O]/dT */

    /*reaction 14: 2.000000 OH (+M) <=> H2O2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[13][0] - 1)*sc[0] + (TB[13][1] - 1)*sc[5] + (TB[13][2] - 1)*sc[13] + (TB[13][3] - 1)*sc[14] + (TB[13][4] - 1)*sc[15] + (TB[13][5] - 1)*sc[26] + (TB[13][6] - 1)*sc[31];
    /* forward */
    phi_f = pow(sc[4], 2.000000);
    k_f = prefactor_units[13] * fwd_A[13]
                * exp(fwd_beta[13] * tc[0] - activation_units[13] * fwd_Ea[13] * invT);
    dlnkfdT = fwd_beta[13] * invT + activation_units[13] * fwd_Ea[13] * invT2;
    /* pressure-fall-off */
    k_0 = low_A[13] * exp(low_beta[13] * tc[0] - activation_units[13] * low_Ea[13] * invT);
    Pr = phase_units[13] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = low_beta[13] * invT + activation_units[13] * low_Ea[13] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(troe_Tsss[13]) > 1.e-100 ? (1.-troe_a[13])*exp(-T/troe_Tsss[13]) : 0.);
    Fcent2 = (fabs(troe_Ts[13]) > 1.e-100 ? troe_a[13] * exp(-T/troe_Ts[13]) : 0.);
    Fcent3 = (troe_len[13] == 4 ? exp(-troe_Tss[13] * invT) : 0.);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        (fabs(troe_Tsss[13]) > 1.e-100 ? -Fcent1/troe_Tsss[13] : 0.)
      + (fabs(troe_Ts[13]) > 1.e-100 ? -Fcent2/troe_Ts[13] : 0.)
      + (troe_len[13] == 4 ? Fcent3*troe_Tss[13]*invT2 : 0.) );
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
        dqdci = (TB[13][0] - 1)*dcdc_fac;
        J[4] += -2 * dqdci;           /* dwdot[OH]/d[H2] */
        J[7] += dqdci;                /* dwdot[H2O2]/d[H2] */
        /* d()/d[OH] */
        dqdci =  + k_f*2.000000*sc[4];
        J[136] += -2 * dqdci;         /* dwdot[OH]/d[OH] */
        J[139] += dqdci;              /* dwdot[H2O2]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (TB[13][1] - 1)*dcdc_fac;
        J[169] += -2 * dqdci;         /* dwdot[OH]/d[H2O] */
        J[172] += dqdci;              /* dwdot[H2O2]/d[H2O] */
        /* d()/d[H2O2] */
        dqdci =  - k_r;
        J[235] += -2 * dqdci;         /* dwdot[OH]/d[H2O2] */
        J[238] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
        /* d()/d[CH4] */
        dqdci = (TB[13][2] - 1)*dcdc_fac;
        J[433] += -2 * dqdci;         /* dwdot[OH]/d[CH4] */
        J[436] += dqdci;              /* dwdot[H2O2]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[13][3] - 1)*dcdc_fac;
        J[466] += -2 * dqdci;         /* dwdot[OH]/d[CO] */
        J[469] += dqdci;              /* dwdot[H2O2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[13][4] - 1)*dcdc_fac;
        J[499] += -2 * dqdci;         /* dwdot[OH]/d[CO2] */
        J[502] += dqdci;              /* dwdot[H2O2]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (TB[13][5] - 1)*dcdc_fac;
        J[862] += -2 * dqdci;         /* dwdot[OH]/d[C2H6] */
        J[865] += dqdci;              /* dwdot[H2O2]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (TB[13][6] - 1)*dcdc_fac;
        J[1027] += -2 * dqdci;        /* dwdot[OH]/d[AR] */
        J[1030] += dqdci;             /* dwdot[H2O2]/d[AR] */
    }
    else {
        dqdc[0] = TB[13][0]*dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac + k_f*2.000000*sc[4];
        dqdc[5] = TB[13][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac - k_r;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = TB[13][2]*dcdc_fac;
        dqdc[14] = TB[13][3]*dcdc_fac;
        dqdc[15] = TB[13][4]*dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac;
        dqdc[22] = dcdc_fac;
        dqdc[23] = dcdc_fac;
        dqdc[24] = dcdc_fac;
        dqdc[25] = dcdc_fac;
        dqdc[26] = TB[13][5]*dcdc_fac;
        dqdc[27] = dcdc_fac;
        dqdc[28] = dcdc_fac;
        dqdc[29] = dcdc_fac;
        dqdc[30] = dcdc_fac;
        dqdc[31] = TB[13][6]*dcdc_fac;
        for (int k=0; k<32; k++) {
            J[33*k+4] += -2 * dqdc[k];
            J[33*k+7] += dqdc[k];
        }
    }
    J[1060] += -2 * dqdT; /* dwdot[OH]/dT */
    J[1063] += dqdT; /* dwdot[H2O2]/dT */

    /*reaction 15: OH + CH3 (+M) <=> CH3OH (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[14][0] - 1)*sc[0] + (TB[14][1] - 1)*sc[5] + (TB[14][2] - 1)*sc[13] + (TB[14][3] - 1)*sc[14] + (TB[14][4] - 1)*sc[15] + (TB[14][5] - 1)*sc[26];
    /* forward */
    phi_f = sc[4]*sc[12];
    k_f = prefactor_units[14] * fwd_A[14]
                * exp(fwd_beta[14] * tc[0] - activation_units[14] * fwd_Ea[14] * invT);
    dlnkfdT = fwd_beta[14] * invT + activation_units[14] * fwd_Ea[14] * invT2;
    /* pressure-fall-off */
    k_0 = low_A[14] * exp(low_beta[14] * tc[0] - activation_units[14] * low_Ea[14] * invT);
    Pr = phase_units[14] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = low_beta[14] * invT + activation_units[14] * low_Ea[14] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(troe_Tsss[14]) > 1.e-100 ? (1.-troe_a[14])*exp(-T/troe_Tsss[14]) : 0.);
    Fcent2 = (fabs(troe_Ts[14]) > 1.e-100 ? troe_a[14] * exp(-T/troe_Ts[14]) : 0.);
    Fcent3 = (troe_len[14] == 4 ? exp(-troe_Tss[14] * invT) : 0.);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        (fabs(troe_Tsss[14]) > 1.e-100 ? -Fcent1/troe_Tsss[14] : 0.)
      + (fabs(troe_Ts[14]) > 1.e-100 ? -Fcent2/troe_Ts[14] : 0.)
      + (troe_len[14] == 4 ? Fcent3*troe_Tss[14]*invT2 : 0.) );
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[20];
    Kc = refCinv * exp(g_RT[4] + g_RT[12] - g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[12]) + (h_RT[20]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[12] -= q; /* CH3 */
    wdot[20] += q; /* CH3OH */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[14][0] - 1)*dcdc_fac;
        J[4] -= dqdci;                /* dwdot[OH]/d[H2] */
        J[12] -= dqdci;               /* dwdot[CH3]/d[H2] */
        J[20] += dqdci;               /* dwdot[CH3OH]/d[H2] */
        /* d()/d[OH] */
        dqdci =  + k_f*sc[12];
        J[136] -= dqdci;              /* dwdot[OH]/d[OH] */
        J[144] -= dqdci;              /* dwdot[CH3]/d[OH] */
        J[152] += dqdci;              /* dwdot[CH3OH]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (TB[14][1] - 1)*dcdc_fac;
        J[169] -= dqdci;              /* dwdot[OH]/d[H2O] */
        J[177] -= dqdci;              /* dwdot[CH3]/d[H2O] */
        J[185] += dqdci;              /* dwdot[CH3OH]/d[H2O] */
        /* d()/d[CH3] */
        dqdci =  + k_f*sc[4];
        J[400] -= dqdci;              /* dwdot[OH]/d[CH3] */
        J[408] -= dqdci;              /* dwdot[CH3]/d[CH3] */
        J[416] += dqdci;              /* dwdot[CH3OH]/d[CH3] */
        /* d()/d[CH4] */
        dqdci = (TB[14][2] - 1)*dcdc_fac;
        J[433] -= dqdci;              /* dwdot[OH]/d[CH4] */
        J[441] -= dqdci;              /* dwdot[CH3]/d[CH4] */
        J[449] += dqdci;              /* dwdot[CH3OH]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[14][3] - 1)*dcdc_fac;
        J[466] -= dqdci;              /* dwdot[OH]/d[CO] */
        J[474] -= dqdci;              /* dwdot[CH3]/d[CO] */
        J[482] += dqdci;              /* dwdot[CH3OH]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[14][4] - 1)*dcdc_fac;
        J[499] -= dqdci;              /* dwdot[OH]/d[CO2] */
        J[507] -= dqdci;              /* dwdot[CH3]/d[CO2] */
        J[515] += dqdci;              /* dwdot[CH3OH]/d[CO2] */
        /* d()/d[CH3OH] */
        dqdci =  - k_r;
        J[664] -= dqdci;              /* dwdot[OH]/d[CH3OH] */
        J[672] -= dqdci;              /* dwdot[CH3]/d[CH3OH] */
        J[680] += dqdci;              /* dwdot[CH3OH]/d[CH3OH] */
        /* d()/d[C2H6] */
        dqdci = (TB[14][5] - 1)*dcdc_fac;
        J[862] -= dqdci;              /* dwdot[OH]/d[C2H6] */
        J[870] -= dqdci;              /* dwdot[CH3]/d[C2H6] */
        J[878] += dqdci;              /* dwdot[CH3OH]/d[C2H6] */
    }
    else {
        dqdc[0] = TB[14][0]*dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac + k_f*sc[12];
        dqdc[5] = TB[14][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac + k_f*sc[4];
        dqdc[13] = TB[14][2]*dcdc_fac;
        dqdc[14] = TB[14][3]*dcdc_fac;
        dqdc[15] = TB[14][4]*dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac - k_r;
        dqdc[21] = dcdc_fac;
        dqdc[22] = dcdc_fac;
        dqdc[23] = dcdc_fac;
        dqdc[24] = dcdc_fac;
        dqdc[25] = dcdc_fac;
        dqdc[26] = TB[14][5]*dcdc_fac;
        dqdc[27] = dcdc_fac;
        dqdc[28] = dcdc_fac;
        dqdc[29] = dcdc_fac;
        dqdc[30] = dcdc_fac;
        dqdc[31] = dcdc_fac;
        for (int k=0; k<32; k++) {
            J[33*k+4] -= dqdc[k];
            J[33*k+12] -= dqdc[k];
            J[33*k+20] += dqdc[k];
        }
    }
    J[1060] -= dqdT; /* dwdot[OH]/dT */
    J[1068] -= dqdT; /* dwdot[CH3]/dT */
    J[1076] += dqdT; /* dwdot[CH3OH]/dT */

    /*reaction 16: CH + CO (+M) <=> HCCO (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[15][0] - 1)*sc[0] + (TB[15][1] - 1)*sc[5] + (TB[15][2] - 1)*sc[13] + (TB[15][3] - 1)*sc[14] + (TB[15][4] - 1)*sc[15] + (TB[15][5] - 1)*sc[26] + (TB[15][6] - 1)*sc[31];
    /* forward */
    phi_f = sc[9]*sc[14];
    k_f = prefactor_units[15] * fwd_A[15]
                * exp(fwd_beta[15] * tc[0] - activation_units[15] * fwd_Ea[15] * invT);
    dlnkfdT = fwd_beta[15] * invT + activation_units[15] * fwd_Ea[15] * invT2;
    /* pressure-fall-off */
    k_0 = low_A[15] * exp(low_beta[15] * tc[0] - activation_units[15] * low_Ea[15] * invT);
    Pr = phase_units[15] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = low_beta[15] * invT + activation_units[15] * low_Ea[15] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(troe_Tsss[15]) > 1.e-100 ? (1.-troe_a[15])*exp(-T/troe_Tsss[15]) : 0.);
    Fcent2 = (fabs(troe_Ts[15]) > 1.e-100 ? troe_a[15] * exp(-T/troe_Ts[15]) : 0.);
    Fcent3 = (troe_len[15] == 4 ? exp(-troe_Tss[15] * invT) : 0.);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        (fabs(troe_Tsss[15]) > 1.e-100 ? -Fcent1/troe_Tsss[15] : 0.)
      + (fabs(troe_Ts[15]) > 1.e-100 ? -Fcent2/troe_Ts[15] : 0.)
      + (troe_len[15] == 4 ? Fcent3*troe_Tss[15]*invT2 : 0.) );
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[27];
    Kc = refCinv * exp(g_RT[9] + g_RT[14] - g_RT[27]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[9] + h_RT[14]) + (h_RT[27]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[9] -= q; /* CH */
    wdot[14] -= q; /* CO */
    wdot[27] += q; /* HCCO */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[15][0] - 1)*dcdc_fac;
        J[9] -= dqdci;                /* dwdot[CH]/d[H2] */
        J[14] -= dqdci;               /* dwdot[CO]/d[H2] */
        J[27] += dqdci;               /* dwdot[HCCO]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[15][1] - 1)*dcdc_fac;
        J[174] -= dqdci;              /* dwdot[CH]/d[H2O] */
        J[179] -= dqdci;              /* dwdot[CO]/d[H2O] */
        J[192] += dqdci;              /* dwdot[HCCO]/d[H2O] */
        /* d()/d[CH] */
        dqdci =  + k_f*sc[14];
        J[306] -= dqdci;              /* dwdot[CH]/d[CH] */
        J[311] -= dqdci;              /* dwdot[CO]/d[CH] */
        J[324] += dqdci;              /* dwdot[HCCO]/d[CH] */
        /* d()/d[CH4] */
        dqdci = (TB[15][2] - 1)*dcdc_fac;
        J[438] -= dqdci;              /* dwdot[CH]/d[CH4] */
        J[443] -= dqdci;              /* dwdot[CO]/d[CH4] */
        J[456] += dqdci;              /* dwdot[HCCO]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[15][3] - 1)*dcdc_fac + k_f*sc[9];
        J[471] -= dqdci;              /* dwdot[CH]/d[CO] */
        J[476] -= dqdci;              /* dwdot[CO]/d[CO] */
        J[489] += dqdci;              /* dwdot[HCCO]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[15][4] - 1)*dcdc_fac;
        J[504] -= dqdci;              /* dwdot[CH]/d[CO2] */
        J[509] -= dqdci;              /* dwdot[CO]/d[CO2] */
        J[522] += dqdci;              /* dwdot[HCCO]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (TB[15][5] - 1)*dcdc_fac;
        J[867] -= dqdci;              /* dwdot[CH]/d[C2H6] */
        J[872] -= dqdci;              /* dwdot[CO]/d[C2H6] */
        J[885] += dqdci;              /* dwdot[HCCO]/d[C2H6] */
        /* d()/d[HCCO] */
        dqdci =  - k_r;
        J[900] -= dqdci;              /* dwdot[CH]/d[HCCO] */
        J[905] -= dqdci;              /* dwdot[CO]/d[HCCO] */
        J[918] += dqdci;              /* dwdot[HCCO]/d[HCCO] */
        /* d()/d[AR] */
        dqdci = (TB[15][6] - 1)*dcdc_fac;
        J[1032] -= dqdci;             /* dwdot[CH]/d[AR] */
        J[1037] -= dqdci;             /* dwdot[CO]/d[AR] */
        J[1050] += dqdci;             /* dwdot[HCCO]/d[AR] */
    }
    else {
        dqdc[0] = TB[15][0]*dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = TB[15][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac + k_f*sc[14];
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = TB[15][2]*dcdc_fac;
        dqdc[14] = TB[15][3]*dcdc_fac + k_f*sc[9];
        dqdc[15] = TB[15][4]*dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac;
        dqdc[22] = dcdc_fac;
        dqdc[23] = dcdc_fac;
        dqdc[24] = dcdc_fac;
        dqdc[25] = dcdc_fac;
        dqdc[26] = TB[15][5]*dcdc_fac;
        dqdc[27] = dcdc_fac - k_r;
        dqdc[28] = dcdc_fac;
        dqdc[29] = dcdc_fac;
        dqdc[30] = dcdc_fac;
        dqdc[31] = TB[15][6]*dcdc_fac;
        for (int k=0; k<32; k++) {
            J[33*k+9] -= dqdc[k];
            J[33*k+14] -= dqdc[k];
            J[33*k+27] += dqdc[k];
        }
    }
    J[1065] -= dqdT; /* dwdot[CH]/dT */
    J[1070] -= dqdT; /* dwdot[CO]/dT */
    J[1083] += dqdT; /* dwdot[HCCO]/dT */

    /*reaction 17: CH2 + CO (+M) <=> CH2CO (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[16][0] - 1)*sc[0] + (TB[16][1] - 1)*sc[5] + (TB[16][2] - 1)*sc[13] + (TB[16][3] - 1)*sc[14] + (TB[16][4] - 1)*sc[15] + (TB[16][5] - 1)*sc[26] + (TB[16][6] - 1)*sc[31];
    /* forward */
    phi_f = sc[10]*sc[14];
    k_f = prefactor_units[16] * fwd_A[16]
                * exp(fwd_beta[16] * tc[0] - activation_units[16] * fwd_Ea[16] * invT);
    dlnkfdT = fwd_beta[16] * invT + activation_units[16] * fwd_Ea[16] * invT2;
    /* pressure-fall-off */
    k_0 = low_A[16] * exp(low_beta[16] * tc[0] - activation_units[16] * low_Ea[16] * invT);
    Pr = phase_units[16] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = low_beta[16] * invT + activation_units[16] * low_Ea[16] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(troe_Tsss[16]) > 1.e-100 ? (1.-troe_a[16])*exp(-T/troe_Tsss[16]) : 0.);
    Fcent2 = (fabs(troe_Ts[16]) > 1.e-100 ? troe_a[16] * exp(-T/troe_Ts[16]) : 0.);
    Fcent3 = (troe_len[16] == 4 ? exp(-troe_Tss[16] * invT) : 0.);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        (fabs(troe_Tsss[16]) > 1.e-100 ? -Fcent1/troe_Tsss[16] : 0.)
      + (fabs(troe_Ts[16]) > 1.e-100 ? -Fcent2/troe_Ts[16] : 0.)
      + (troe_len[16] == 4 ? Fcent3*troe_Tss[16]*invT2 : 0.) );
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[28];
    Kc = refCinv * exp(g_RT[10] + g_RT[14] - g_RT[28]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[10] + h_RT[14]) + (h_RT[28]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[10] -= q; /* CH2 */
    wdot[14] -= q; /* CO */
    wdot[28] += q; /* CH2CO */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[16][0] - 1)*dcdc_fac;
        J[10] -= dqdci;               /* dwdot[CH2]/d[H2] */
        J[14] -= dqdci;               /* dwdot[CO]/d[H2] */
        J[28] += dqdci;               /* dwdot[CH2CO]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[16][1] - 1)*dcdc_fac;
        J[175] -= dqdci;              /* dwdot[CH2]/d[H2O] */
        J[179] -= dqdci;              /* dwdot[CO]/d[H2O] */
        J[193] += dqdci;              /* dwdot[CH2CO]/d[H2O] */
        /* d()/d[CH2] */
        dqdci =  + k_f*sc[14];
        J[340] -= dqdci;              /* dwdot[CH2]/d[CH2] */
        J[344] -= dqdci;              /* dwdot[CO]/d[CH2] */
        J[358] += dqdci;              /* dwdot[CH2CO]/d[CH2] */
        /* d()/d[CH4] */
        dqdci = (TB[16][2] - 1)*dcdc_fac;
        J[439] -= dqdci;              /* dwdot[CH2]/d[CH4] */
        J[443] -= dqdci;              /* dwdot[CO]/d[CH4] */
        J[457] += dqdci;              /* dwdot[CH2CO]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[16][3] - 1)*dcdc_fac + k_f*sc[10];
        J[472] -= dqdci;              /* dwdot[CH2]/d[CO] */
        J[476] -= dqdci;              /* dwdot[CO]/d[CO] */
        J[490] += dqdci;              /* dwdot[CH2CO]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[16][4] - 1)*dcdc_fac;
        J[505] -= dqdci;              /* dwdot[CH2]/d[CO2] */
        J[509] -= dqdci;              /* dwdot[CO]/d[CO2] */
        J[523] += dqdci;              /* dwdot[CH2CO]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (TB[16][5] - 1)*dcdc_fac;
        J[868] -= dqdci;              /* dwdot[CH2]/d[C2H6] */
        J[872] -= dqdci;              /* dwdot[CO]/d[C2H6] */
        J[886] += dqdci;              /* dwdot[CH2CO]/d[C2H6] */
        /* d()/d[CH2CO] */
        dqdci =  - k_r;
        J[934] -= dqdci;              /* dwdot[CH2]/d[CH2CO] */
        J[938] -= dqdci;              /* dwdot[CO]/d[CH2CO] */
        J[952] += dqdci;              /* dwdot[CH2CO]/d[CH2CO] */
        /* d()/d[AR] */
        dqdci = (TB[16][6] - 1)*dcdc_fac;
        J[1033] -= dqdci;             /* dwdot[CH2]/d[AR] */
        J[1037] -= dqdci;             /* dwdot[CO]/d[AR] */
        J[1051] += dqdci;             /* dwdot[CH2CO]/d[AR] */
    }
    else {
        dqdc[0] = TB[16][0]*dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = TB[16][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac + k_f*sc[14];
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = TB[16][2]*dcdc_fac;
        dqdc[14] = TB[16][3]*dcdc_fac + k_f*sc[10];
        dqdc[15] = TB[16][4]*dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac;
        dqdc[22] = dcdc_fac;
        dqdc[23] = dcdc_fac;
        dqdc[24] = dcdc_fac;
        dqdc[25] = dcdc_fac;
        dqdc[26] = TB[16][5]*dcdc_fac;
        dqdc[27] = dcdc_fac;
        dqdc[28] = dcdc_fac - k_r;
        dqdc[29] = dcdc_fac;
        dqdc[30] = dcdc_fac;
        dqdc[31] = TB[16][6]*dcdc_fac;
        for (int k=0; k<32; k++) {
            J[33*k+10] -= dqdc[k];
            J[33*k+14] -= dqdc[k];
            J[33*k+28] += dqdc[k];
        }
    }
    J[1066] -= dqdT; /* dwdot[CH2]/dT */
    J[1070] -= dqdT; /* dwdot[CO]/dT */
    J[1084] += dqdT; /* dwdot[CH2CO]/dT */

    /*reaction 18: CH2(S) + H2O (+M) <=> CH3OH (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[17][0] - 1)*sc[0] + (TB[17][1] - 1)*sc[5] + (TB[17][2] - 1)*sc[13] + (TB[17][3] - 1)*sc[14] + (TB[17][4] - 1)*sc[15] + (TB[17][5] - 1)*sc[26];
    /* forward */
    phi_f = sc[5]*sc[11];
    k_f = prefactor_units[17] * fwd_A[17]
                * exp(fwd_beta[17] * tc[0] - activation_units[17] * fwd_Ea[17] * invT);
    dlnkfdT = fwd_beta[17] * invT + activation_units[17] * fwd_Ea[17] * invT2;
    /* pressure-fall-off */
    k_0 = low_A[17] * exp(low_beta[17] * tc[0] - activation_units[17] * low_Ea[17] * invT);
    Pr = phase_units[17] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = low_beta[17] * invT + activation_units[17] * low_Ea[17] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(troe_Tsss[17]) > 1.e-100 ? (1.-troe_a[17])*exp(-T/troe_Tsss[17]) : 0.);
    Fcent2 = (fabs(troe_Ts[17]) > 1.e-100 ? troe_a[17] * exp(-T/troe_Ts[17]) : 0.);
    Fcent3 = (troe_len[17] == 4 ? exp(-troe_Tss[17] * invT) : 0.);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        (fabs(troe_Tsss[17]) > 1.e-100 ? -Fcent1/troe_Tsss[17] : 0.)
      + (fabs(troe_Ts[17]) > 1.e-100 ? -Fcent2/troe_Ts[17] : 0.)
      + (troe_len[17] == 4 ? Fcent3*troe_Tss[17]*invT2 : 0.) );
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[20];
    Kc = refCinv * exp(g_RT[5] + g_RT[11] - g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[11]) + (h_RT[20]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[5] -= q; /* H2O */
    wdot[11] -= q; /* CH2(S) */
    wdot[20] += q; /* CH3OH */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[17][0] - 1)*dcdc_fac;
        J[5] -= dqdci;                /* dwdot[H2O]/d[H2] */
        J[11] -= dqdci;               /* dwdot[CH2(S)]/d[H2] */
        J[20] += dqdci;               /* dwdot[CH3OH]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[17][1] - 1)*dcdc_fac + k_f*sc[11];
        J[170] -= dqdci;              /* dwdot[H2O]/d[H2O] */
        J[176] -= dqdci;              /* dwdot[CH2(S)]/d[H2O] */
        J[185] += dqdci;              /* dwdot[CH3OH]/d[H2O] */
        /* d()/d[CH2(S)] */
        dqdci =  + k_f*sc[5];
        J[368] -= dqdci;              /* dwdot[H2O]/d[CH2(S)] */
        J[374] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
        J[383] += dqdci;              /* dwdot[CH3OH]/d[CH2(S)] */
        /* d()/d[CH4] */
        dqdci = (TB[17][2] - 1)*dcdc_fac;
        J[434] -= dqdci;              /* dwdot[H2O]/d[CH4] */
        J[440] -= dqdci;              /* dwdot[CH2(S)]/d[CH4] */
        J[449] += dqdci;              /* dwdot[CH3OH]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[17][3] - 1)*dcdc_fac;
        J[467] -= dqdci;              /* dwdot[H2O]/d[CO] */
        J[473] -= dqdci;              /* dwdot[CH2(S)]/d[CO] */
        J[482] += dqdci;              /* dwdot[CH3OH]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[17][4] - 1)*dcdc_fac;
        J[500] -= dqdci;              /* dwdot[H2O]/d[CO2] */
        J[506] -= dqdci;              /* dwdot[CH2(S)]/d[CO2] */
        J[515] += dqdci;              /* dwdot[CH3OH]/d[CO2] */
        /* d()/d[CH3OH] */
        dqdci =  - k_r;
        J[665] -= dqdci;              /* dwdot[H2O]/d[CH3OH] */
        J[671] -= dqdci;              /* dwdot[CH2(S)]/d[CH3OH] */
        J[680] += dqdci;              /* dwdot[CH3OH]/d[CH3OH] */
        /* d()/d[C2H6] */
        dqdci = (TB[17][5] - 1)*dcdc_fac;
        J[863] -= dqdci;              /* dwdot[H2O]/d[C2H6] */
        J[869] -= dqdci;              /* dwdot[CH2(S)]/d[C2H6] */
        J[878] += dqdci;              /* dwdot[CH3OH]/d[C2H6] */
    }
    else {
        dqdc[0] = TB[17][0]*dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = TB[17][1]*dcdc_fac + k_f*sc[11];
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac + k_f*sc[5];
        dqdc[12] = dcdc_fac;
        dqdc[13] = TB[17][2]*dcdc_fac;
        dqdc[14] = TB[17][3]*dcdc_fac;
        dqdc[15] = TB[17][4]*dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac - k_r;
        dqdc[21] = dcdc_fac;
        dqdc[22] = dcdc_fac;
        dqdc[23] = dcdc_fac;
        dqdc[24] = dcdc_fac;
        dqdc[25] = dcdc_fac;
        dqdc[26] = TB[17][5]*dcdc_fac;
        dqdc[27] = dcdc_fac;
        dqdc[28] = dcdc_fac;
        dqdc[29] = dcdc_fac;
        dqdc[30] = dcdc_fac;
        dqdc[31] = dcdc_fac;
        for (int k=0; k<32; k++) {
            J[33*k+5] -= dqdc[k];
            J[33*k+11] -= dqdc[k];
            J[33*k+20] += dqdc[k];
        }
    }
    J[1061] -= dqdT; /* dwdot[H2O]/dT */
    J[1067] -= dqdT; /* dwdot[CH2(S)]/dT */
    J[1076] += dqdT; /* dwdot[CH3OH]/dT */

    /*reaction 19: 2.000000 CH3 (+M) <=> C2H6 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[18][0] - 1)*sc[0] + (TB[18][1] - 1)*sc[5] + (TB[18][2] - 1)*sc[13] + (TB[18][3] - 1)*sc[14] + (TB[18][4] - 1)*sc[15] + (TB[18][5] - 1)*sc[26] + (TB[18][6] - 1)*sc[31];
    /* forward */
    phi_f = pow(sc[12], 2.000000);
    k_f = prefactor_units[18] * fwd_A[18]
                * exp(fwd_beta[18] * tc[0] - activation_units[18] * fwd_Ea[18] * invT);
    dlnkfdT = fwd_beta[18] * invT + activation_units[18] * fwd_Ea[18] * invT2;
    /* pressure-fall-off */
    k_0 = low_A[18] * exp(low_beta[18] * tc[0] - activation_units[18] * low_Ea[18] * invT);
    Pr = phase_units[18] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = low_beta[18] * invT + activation_units[18] * low_Ea[18] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(troe_Tsss[18]) > 1.e-100 ? (1.-troe_a[18])*exp(-T/troe_Tsss[18]) : 0.);
    Fcent2 = (fabs(troe_Ts[18]) > 1.e-100 ? troe_a[18] * exp(-T/troe_Ts[18]) : 0.);
    Fcent3 = (troe_len[18] == 4 ? exp(-troe_Tss[18] * invT) : 0.);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        (fabs(troe_Tsss[18]) > 1.e-100 ? -Fcent1/troe_Tsss[18] : 0.)
      + (fabs(troe_Ts[18]) > 1.e-100 ? -Fcent2/troe_Ts[18] : 0.)
      + (troe_len[18] == 4 ? Fcent3*troe_Tss[18]*invT2 : 0.) );
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[26];
    Kc = refCinv * exp(2.000000*g_RT[12] - g_RT[26]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[12]) + (h_RT[26]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[12] -= 2 * q; /* CH3 */
    wdot[26] += q; /* C2H6 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[18][0] - 1)*dcdc_fac;
        J[12] += -2 * dqdci;          /* dwdot[CH3]/d[H2] */
        J[26] += dqdci;               /* dwdot[C2H6]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[18][1] - 1)*dcdc_fac;
        J[177] += -2 * dqdci;         /* dwdot[CH3]/d[H2O] */
        J[191] += dqdci;              /* dwdot[C2H6]/d[H2O] */
        /* d()/d[CH3] */
        dqdci =  + k_f*2.000000*sc[12];
        J[408] += -2 * dqdci;         /* dwdot[CH3]/d[CH3] */
        J[422] += dqdci;              /* dwdot[C2H6]/d[CH3] */
        /* d()/d[CH4] */
        dqdci = (TB[18][2] - 1)*dcdc_fac;
        J[441] += -2 * dqdci;         /* dwdot[CH3]/d[CH4] */
        J[455] += dqdci;              /* dwdot[C2H6]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[18][3] - 1)*dcdc_fac;
        J[474] += -2 * dqdci;         /* dwdot[CH3]/d[CO] */
        J[488] += dqdci;              /* dwdot[C2H6]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[18][4] - 1)*dcdc_fac;
        J[507] += -2 * dqdci;         /* dwdot[CH3]/d[CO2] */
        J[521] += dqdci;              /* dwdot[C2H6]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (TB[18][5] - 1)*dcdc_fac - k_r;
        J[870] += -2 * dqdci;         /* dwdot[CH3]/d[C2H6] */
        J[884] += dqdci;              /* dwdot[C2H6]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (TB[18][6] - 1)*dcdc_fac;
        J[1035] += -2 * dqdci;        /* dwdot[CH3]/d[AR] */
        J[1049] += dqdci;             /* dwdot[C2H6]/d[AR] */
    }
    else {
        dqdc[0] = TB[18][0]*dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = TB[18][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac + k_f*2.000000*sc[12];
        dqdc[13] = TB[18][2]*dcdc_fac;
        dqdc[14] = TB[18][3]*dcdc_fac;
        dqdc[15] = TB[18][4]*dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac;
        dqdc[22] = dcdc_fac;
        dqdc[23] = dcdc_fac;
        dqdc[24] = dcdc_fac;
        dqdc[25] = dcdc_fac;
        dqdc[26] = TB[18][5]*dcdc_fac - k_r;
        dqdc[27] = dcdc_fac;
        dqdc[28] = dcdc_fac;
        dqdc[29] = dcdc_fac;
        dqdc[30] = dcdc_fac;
        dqdc[31] = TB[18][6]*dcdc_fac;
        for (int k=0; k<32; k++) {
            J[33*k+12] += -2 * dqdc[k];
            J[33*k+26] += dqdc[k];
        }
    }
    J[1068] += -2 * dqdT; /* dwdot[CH3]/dT */
    J[1082] += dqdT; /* dwdot[C2H6]/dT */

    /*reaction 20: C2H4 (+M) <=> H2 + C2H2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[19][0] - 1)*sc[0] + (TB[19][1] - 1)*sc[5] + (TB[19][2] - 1)*sc[13] + (TB[19][3] - 1)*sc[14] + (TB[19][4] - 1)*sc[15] + (TB[19][5] - 1)*sc[26] + (TB[19][6] - 1)*sc[31];
    /* forward */
    phi_f = sc[24];
    k_f = prefactor_units[19] * fwd_A[19]
                * exp(fwd_beta[19] * tc[0] - activation_units[19] * fwd_Ea[19] * invT);
    dlnkfdT = fwd_beta[19] * invT + activation_units[19] * fwd_Ea[19] * invT2;
    /* pressure-fall-off */
    k_0 = low_A[19] * exp(low_beta[19] * tc[0] - activation_units[19] * low_Ea[19] * invT);
    Pr = phase_units[19] * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = low_beta[19] * invT + activation_units[19] * low_Ea[19] * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (fabs(troe_Tsss[19]) > 1.e-100 ? (1.-troe_a[19])*exp(-T/troe_Tsss[19]) : 0.);
    Fcent2 = (fabs(troe_Ts[19]) > 1.e-100 ? troe_a[19] * exp(-T/troe_Ts[19]) : 0.);
    Fcent3 = (troe_len[19] == 4 ? exp(-troe_Tss[19] * invT) : 0.);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        (fabs(troe_Tsss[19]) > 1.e-100 ? -Fcent1/troe_Tsss[19] : 0.)
      + (fabs(troe_Ts[19]) > 1.e-100 ? -Fcent2/troe_Ts[19] : 0.)
      + (troe_len[19] == 4 ? Fcent3*troe_Tss[19]*invT2 : 0.) );
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[0]*sc[22];
    Kc = refC * exp(-g_RT[0] - g_RT[22] + g_RT[24]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[24]) + (h_RT[0] + h_RT[22]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[22] += q; /* C2H2 */
    wdot[24] -= q; /* C2H4 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[19][0] - 1)*dcdc_fac - k_r*sc[22];
        J[0] += dqdci;                /* dwdot[H2]/d[H2] */
        J[22] += dqdci;               /* dwdot[C2H2]/d[H2] */
        J[24] -= dqdci;               /* dwdot[C2H4]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[19][1] - 1)*dcdc_fac;
        J[165] += dqdci;              /* dwdot[H2]/d[H2O] */
        J[187] += dqdci;              /* dwdot[C2H2]/d[H2O] */
        J[189] -= dqdci;              /* dwdot[C2H4]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (TB[19][2] - 1)*dcdc_fac;
        J[429] += dqdci;              /* dwdot[H2]/d[CH4] */
        J[451] += dqdci;              /* dwdot[C2H2]/d[CH4] */
        J[453] -= dqdci;              /* dwdot[C2H4]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[19][3] - 1)*dcdc_fac;
        J[462] += dqdci;              /* dwdot[H2]/d[CO] */
        J[484] += dqdci;              /* dwdot[C2H2]/d[CO] */
        J[486] -= dqdci;              /* dwdot[C2H4]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[19][4] - 1)*dcdc_fac;
        J[495] += dqdci;              /* dwdot[H2]/d[CO2] */
        J[517] += dqdci;              /* dwdot[C2H2]/d[CO2] */
        J[519] -= dqdci;              /* dwdot[C2H4]/d[CO2] */
        /* d()/d[C2H2] */
        dqdci =  - k_r*sc[0];
        J[726] += dqdci;              /* dwdot[H2]/d[C2H2] */
        J[748] += dqdci;              /* dwdot[C2H2]/d[C2H2] */
        J[750] -= dqdci;              /* dwdot[C2H4]/d[C2H2] */
        /* d()/d[C2H4] */
        dqdci =  + k_f;
        J[792] += dqdci;              /* dwdot[H2]/d[C2H4] */
        J[814] += dqdci;              /* dwdot[C2H2]/d[C2H4] */
        J[816] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
        /* d()/d[C2H6] */
        dqdci = (TB[19][5] - 1)*dcdc_fac;
        J[858] += dqdci;              /* dwdot[H2]/d[C2H6] */
        J[880] += dqdci;              /* dwdot[C2H2]/d[C2H6] */
        J[882] -= dqdci;              /* dwdot[C2H4]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (TB[19][6] - 1)*dcdc_fac;
        J[1023] += dqdci;             /* dwdot[H2]/d[AR] */
        J[1045] += dqdci;             /* dwdot[C2H2]/d[AR] */
        J[1047] -= dqdci;             /* dwdot[C2H4]/d[AR] */
    }
    else {
        dqdc[0] = TB[19][0]*dcdc_fac - k_r*sc[22];
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = TB[19][1]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = TB[19][2]*dcdc_fac;
        dqdc[14] = TB[19][3]*dcdc_fac;
        dqdc[15] = TB[19][4]*dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac;
        dqdc[22] = dcdc_fac - k_r*sc[0];
        dqdc[23] = dcdc_fac;
        dqdc[24] = dcdc_fac + k_f;
        dqdc[25] = dcdc_fac;
        dqdc[26] = TB[19][5]*dcdc_fac;
        dqdc[27] = dcdc_fac;
        dqdc[28] = dcdc_fac;
        dqdc[29] = dcdc_fac;
        dqdc[30] = dcdc_fac;
        dqdc[31] = TB[19][6]*dcdc_fac;
        for (int k=0; k<32; k++) {
            J[33*k+0] += dqdc[k];
            J[33*k+22] += dqdc[k];
            J[33*k+24] -= dqdc[k];
        }
    }
    J[1056] += dqdT; /* dwdot[H2]/dT */
    J[1078] += dqdT; /* dwdot[C2H2]/dT */
    J[1080] -= dqdT; /* dwdot[C2H4]/dT */

    /*reaction 21: 2.000000 O + M <=> O2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[20][0] - 1)*sc[0] + (TB[20][1] - 1)*sc[5] + (TB[20][2] - 1)*sc[13] + (TB[20][3] - 1)*sc[14] + (TB[20][4] - 1)*sc[15] + (TB[20][5] - 1)*sc[26] + (TB[20][6] - 1)*sc[31];
    /* forward */
    phi_f = pow(sc[2], 2.000000);
    k_f = prefactor_units[20] * fwd_A[20]
                * exp(fwd_beta[20] * tc[0] - activation_units[20] * fwd_Ea[20] * invT);
    dlnkfdT = fwd_beta[20] * invT + activation_units[20] * fwd_Ea[20] * invT2;
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
        dqdci = (TB[20][0] - 1)*q_nocor;
        J[2] += -2 * dqdci;           /* dwdot[O]/d[H2] */
        J[3] += dqdci;                /* dwdot[O2]/d[H2] */
        /* d()/d[O] */
        dqdci =  + k_f*2.000000*sc[2];
        J[68] += -2 * dqdci;          /* dwdot[O]/d[O] */
        J[69] += dqdci;               /* dwdot[O2]/d[O] */
        /* d()/d[O2] */
        dqdci =  - k_r;
        J[101] += -2 * dqdci;         /* dwdot[O]/d[O2] */
        J[102] += dqdci;              /* dwdot[O2]/d[O2] */
        /* d()/d[H2O] */
        dqdci = (TB[20][1] - 1)*q_nocor;
        J[167] += -2 * dqdci;         /* dwdot[O]/d[H2O] */
        J[168] += dqdci;              /* dwdot[O2]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (TB[20][2] - 1)*q_nocor;
        J[431] += -2 * dqdci;         /* dwdot[O]/d[CH4] */
        J[432] += dqdci;              /* dwdot[O2]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[20][3] - 1)*q_nocor;
        J[464] += -2 * dqdci;         /* dwdot[O]/d[CO] */
        J[465] += dqdci;              /* dwdot[O2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[20][4] - 1)*q_nocor;
        J[497] += -2 * dqdci;         /* dwdot[O]/d[CO2] */
        J[498] += dqdci;              /* dwdot[O2]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (TB[20][5] - 1)*q_nocor;
        J[860] += -2 * dqdci;         /* dwdot[O]/d[C2H6] */
        J[861] += dqdci;              /* dwdot[O2]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (TB[20][6] - 1)*q_nocor;
        J[1025] += -2 * dqdci;        /* dwdot[O]/d[AR] */
        J[1026] += dqdci;             /* dwdot[O2]/d[AR] */
    }
    else {
        dqdc[0] = TB[20][0]*q_nocor;
        dqdc[1] = q_nocor;
        dqdc[2] = q_nocor + k_f*2.000000*sc[2];
        dqdc[3] = q_nocor - k_r;
        dqdc[4] = q_nocor;
        dqdc[5] = TB[20][1]*q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = q_nocor;
        dqdc[12] = q_nocor;
        dqdc[13] = TB[20][2]*q_nocor;
        dqdc[14] = TB[20][3]*q_nocor;
        dqdc[15] = TB[20][4]*q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = q_nocor;
        dqdc[19] = q_nocor;
        dqdc[20] = q_nocor;
        dqdc[21] = q_nocor;
        dqdc[22] = q_nocor;
        dqdc[23] = q_nocor;
        dqdc[24] = q_nocor;
        dqdc[25] = q_nocor;
        dqdc[26] = TB[20][5]*q_nocor;
        dqdc[27] = q_nocor;
        dqdc[28] = q_nocor;
        dqdc[29] = q_nocor;
        dqdc[30] = q_nocor;
        dqdc[31] = TB[20][6]*q_nocor;
        for (int k=0; k<32; k++) {
            J[33*k+2] += -2 * dqdc[k];
            J[33*k+3] += dqdc[k];
        }
    }
    J[1058] += -2 * dqdT; /* dwdot[O]/dT */
    J[1059] += dqdT; /* dwdot[O2]/dT */

    /*reaction 22: O + H + M <=> OH + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[21][0] - 1)*sc[0] + (TB[21][1] - 1)*sc[5] + (TB[21][2] - 1)*sc[13] + (TB[21][3] - 1)*sc[14] + (TB[21][4] - 1)*sc[15] + (TB[21][5] - 1)*sc[26] + (TB[21][6] - 1)*sc[31];
    /* forward */
    phi_f = sc[1]*sc[2];
    k_f = prefactor_units[21] * fwd_A[21]
                * exp(fwd_beta[21] * tc[0] - activation_units[21] * fwd_Ea[21] * invT);
    dlnkfdT = fwd_beta[21] * invT + activation_units[21] * fwd_Ea[21] * invT2;
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
        dqdci = (TB[21][0] - 1)*q_nocor;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[2] -= dqdci;                /* dwdot[O]/d[H2] */
        J[4] += dqdci;                /* dwdot[OH]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[2];
        J[34] -= dqdci;               /* dwdot[H]/d[H] */
        J[35] -= dqdci;               /* dwdot[O]/d[H] */
        J[37] += dqdci;               /* dwdot[OH]/d[H] */
        /* d()/d[O] */
        dqdci =  + k_f*sc[1];
        J[67] -= dqdci;               /* dwdot[H]/d[O] */
        J[68] -= dqdci;               /* dwdot[O]/d[O] */
        J[70] += dqdci;               /* dwdot[OH]/d[O] */
        /* d()/d[OH] */
        dqdci =  - k_r;
        J[133] -= dqdci;              /* dwdot[H]/d[OH] */
        J[134] -= dqdci;              /* dwdot[O]/d[OH] */
        J[136] += dqdci;              /* dwdot[OH]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (TB[21][1] - 1)*q_nocor;
        J[166] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[167] -= dqdci;              /* dwdot[O]/d[H2O] */
        J[169] += dqdci;              /* dwdot[OH]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (TB[21][2] - 1)*q_nocor;
        J[430] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[431] -= dqdci;              /* dwdot[O]/d[CH4] */
        J[433] += dqdci;              /* dwdot[OH]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[21][3] - 1)*q_nocor;
        J[463] -= dqdci;              /* dwdot[H]/d[CO] */
        J[464] -= dqdci;              /* dwdot[O]/d[CO] */
        J[466] += dqdci;              /* dwdot[OH]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[21][4] - 1)*q_nocor;
        J[496] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[497] -= dqdci;              /* dwdot[O]/d[CO2] */
        J[499] += dqdci;              /* dwdot[OH]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (TB[21][5] - 1)*q_nocor;
        J[859] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[860] -= dqdci;              /* dwdot[O]/d[C2H6] */
        J[862] += dqdci;              /* dwdot[OH]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (TB[21][6] - 1)*q_nocor;
        J[1024] -= dqdci;             /* dwdot[H]/d[AR] */
        J[1025] -= dqdci;             /* dwdot[O]/d[AR] */
        J[1027] += dqdci;             /* dwdot[OH]/d[AR] */
    }
    else {
        dqdc[0] = TB[21][0]*q_nocor;
        dqdc[1] = q_nocor + k_f*sc[2];
        dqdc[2] = q_nocor + k_f*sc[1];
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor - k_r;
        dqdc[5] = TB[21][1]*q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = q_nocor;
        dqdc[12] = q_nocor;
        dqdc[13] = TB[21][2]*q_nocor;
        dqdc[14] = TB[21][3]*q_nocor;
        dqdc[15] = TB[21][4]*q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = q_nocor;
        dqdc[19] = q_nocor;
        dqdc[20] = q_nocor;
        dqdc[21] = q_nocor;
        dqdc[22] = q_nocor;
        dqdc[23] = q_nocor;
        dqdc[24] = q_nocor;
        dqdc[25] = q_nocor;
        dqdc[26] = TB[21][5]*q_nocor;
        dqdc[27] = q_nocor;
        dqdc[28] = q_nocor;
        dqdc[29] = q_nocor;
        dqdc[30] = q_nocor;
        dqdc[31] = TB[21][6]*q_nocor;
        for (int k=0; k<32; k++) {
            J[33*k+1] -= dqdc[k];
            J[33*k+2] -= dqdc[k];
            J[33*k+4] += dqdc[k];
        }
    }
    J[1057] -= dqdT; /* dwdot[H]/dT */
    J[1058] -= dqdT; /* dwdot[O]/dT */
    J[1060] += dqdT; /* dwdot[OH]/dT */

    /*reaction 23: O + CO + M <=> CO2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[22][0] - 1)*sc[0] + (TB[22][1] - 1)*sc[3] + (TB[22][2] - 1)*sc[5] + (TB[22][3] - 1)*sc[13] + (TB[22][4] - 1)*sc[14] + (TB[22][5] - 1)*sc[15] + (TB[22][6] - 1)*sc[26] + (TB[22][7] - 1)*sc[31];
    /* forward */
    phi_f = sc[2]*sc[14];
    k_f = prefactor_units[22] * fwd_A[22]
                * exp(fwd_beta[22] * tc[0] - activation_units[22] * fwd_Ea[22] * invT);
    dlnkfdT = fwd_beta[22] * invT + activation_units[22] * fwd_Ea[22] * invT2;
    /* reverse */
    phi_r = sc[15];
    Kc = refCinv * exp(g_RT[2] + g_RT[14] - g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[14]) + (h_RT[15]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[14] -= q; /* CO */
    wdot[15] += q; /* CO2 */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[22][0] - 1)*q_nocor;
        J[2] -= dqdci;                /* dwdot[O]/d[H2] */
        J[14] -= dqdci;               /* dwdot[CO]/d[H2] */
        J[15] += dqdci;               /* dwdot[CO2]/d[H2] */
        /* d()/d[O] */
        dqdci =  + k_f*sc[14];
        J[68] -= dqdci;               /* dwdot[O]/d[O] */
        J[80] -= dqdci;               /* dwdot[CO]/d[O] */
        J[81] += dqdci;               /* dwdot[CO2]/d[O] */
        /* d()/d[O2] */
        dqdci = (TB[22][1] - 1)*q_nocor;
        J[101] -= dqdci;              /* dwdot[O]/d[O2] */
        J[113] -= dqdci;              /* dwdot[CO]/d[O2] */
        J[114] += dqdci;              /* dwdot[CO2]/d[O2] */
        /* d()/d[H2O] */
        dqdci = (TB[22][2] - 1)*q_nocor;
        J[167] -= dqdci;              /* dwdot[O]/d[H2O] */
        J[179] -= dqdci;              /* dwdot[CO]/d[H2O] */
        J[180] += dqdci;              /* dwdot[CO2]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (TB[22][3] - 1)*q_nocor;
        J[431] -= dqdci;              /* dwdot[O]/d[CH4] */
        J[443] -= dqdci;              /* dwdot[CO]/d[CH4] */
        J[444] += dqdci;              /* dwdot[CO2]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[22][4] - 1)*q_nocor + k_f*sc[2];
        J[464] -= dqdci;              /* dwdot[O]/d[CO] */
        J[476] -= dqdci;              /* dwdot[CO]/d[CO] */
        J[477] += dqdci;              /* dwdot[CO2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[22][5] - 1)*q_nocor - k_r;
        J[497] -= dqdci;              /* dwdot[O]/d[CO2] */
        J[509] -= dqdci;              /* dwdot[CO]/d[CO2] */
        J[510] += dqdci;              /* dwdot[CO2]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (TB[22][6] - 1)*q_nocor;
        J[860] -= dqdci;              /* dwdot[O]/d[C2H6] */
        J[872] -= dqdci;              /* dwdot[CO]/d[C2H6] */
        J[873] += dqdci;              /* dwdot[CO2]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (TB[22][7] - 1)*q_nocor;
        J[1025] -= dqdci;             /* dwdot[O]/d[AR] */
        J[1037] -= dqdci;             /* dwdot[CO]/d[AR] */
        J[1038] += dqdci;             /* dwdot[CO2]/d[AR] */
    }
    else {
        dqdc[0] = TB[22][0]*q_nocor;
        dqdc[1] = q_nocor;
        dqdc[2] = q_nocor + k_f*sc[14];
        dqdc[3] = TB[22][1]*q_nocor;
        dqdc[4] = q_nocor;
        dqdc[5] = TB[22][2]*q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = q_nocor;
        dqdc[12] = q_nocor;
        dqdc[13] = TB[22][3]*q_nocor;
        dqdc[14] = TB[22][4]*q_nocor + k_f*sc[2];
        dqdc[15] = TB[22][5]*q_nocor - k_r;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = q_nocor;
        dqdc[19] = q_nocor;
        dqdc[20] = q_nocor;
        dqdc[21] = q_nocor;
        dqdc[22] = q_nocor;
        dqdc[23] = q_nocor;
        dqdc[24] = q_nocor;
        dqdc[25] = q_nocor;
        dqdc[26] = TB[22][6]*q_nocor;
        dqdc[27] = q_nocor;
        dqdc[28] = q_nocor;
        dqdc[29] = q_nocor;
        dqdc[30] = q_nocor;
        dqdc[31] = TB[22][7]*q_nocor;
        for (int k=0; k<32; k++) {
            J[33*k+2] -= dqdc[k];
            J[33*k+14] -= dqdc[k];
            J[33*k+15] += dqdc[k];
        }
    }
    J[1058] -= dqdT; /* dwdot[O]/dT */
    J[1070] -= dqdT; /* dwdot[CO]/dT */
    J[1071] += dqdT; /* dwdot[CO2]/dT */

    /*reaction 24: H + O2 + M <=> HO2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[23][0] - 1)*sc[3] + (TB[23][1] - 1)*sc[5] + (TB[23][2] - 1)*sc[14] + (TB[23][3] - 1)*sc[15] + (TB[23][4] - 1)*sc[26] + (TB[23][5] - 1)*sc[30] + (TB[23][6] - 1)*sc[31];
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = prefactor_units[23] * fwd_A[23]
                * exp(fwd_beta[23] * tc[0] - activation_units[23] * fwd_Ea[23] * invT);
    dlnkfdT = fwd_beta[23] * invT + activation_units[23] * fwd_Ea[23] * invT2;
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
        J[34] -= dqdci;               /* dwdot[H]/d[H] */
        J[36] -= dqdci;               /* dwdot[O2]/d[H] */
        J[39] += dqdci;               /* dwdot[HO2]/d[H] */
        /* d()/d[O2] */
        dqdci = (TB[23][0] - 1)*q_nocor + k_f*sc[1];
        J[100] -= dqdci;              /* dwdot[H]/d[O2] */
        J[102] -= dqdci;              /* dwdot[O2]/d[O2] */
        J[105] += dqdci;              /* dwdot[HO2]/d[O2] */
        /* d()/d[H2O] */
        dqdci = (TB[23][1] - 1)*q_nocor;
        J[166] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[168] -= dqdci;              /* dwdot[O2]/d[H2O] */
        J[171] += dqdci;              /* dwdot[HO2]/d[H2O] */
        /* d()/d[HO2] */
        dqdci =  - k_r;
        J[199] -= dqdci;              /* dwdot[H]/d[HO2] */
        J[201] -= dqdci;              /* dwdot[O2]/d[HO2] */
        J[204] += dqdci;              /* dwdot[HO2]/d[HO2] */
        /* d()/d[CO] */
        dqdci = (TB[23][2] - 1)*q_nocor;
        J[463] -= dqdci;              /* dwdot[H]/d[CO] */
        J[465] -= dqdci;              /* dwdot[O2]/d[CO] */
        J[468] += dqdci;              /* dwdot[HO2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[23][3] - 1)*q_nocor;
        J[496] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[498] -= dqdci;              /* dwdot[O2]/d[CO2] */
        J[501] += dqdci;              /* dwdot[HO2]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (TB[23][4] - 1)*q_nocor;
        J[859] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[861] -= dqdci;              /* dwdot[O2]/d[C2H6] */
        J[864] += dqdci;              /* dwdot[HO2]/d[C2H6] */
        /* d()/d[N2] */
        dqdci = (TB[23][5] - 1)*q_nocor;
        J[991] -= dqdci;              /* dwdot[H]/d[N2] */
        J[993] -= dqdci;              /* dwdot[O2]/d[N2] */
        J[996] += dqdci;              /* dwdot[HO2]/d[N2] */
        /* d()/d[AR] */
        dqdci = (TB[23][6] - 1)*q_nocor;
        J[1024] -= dqdci;             /* dwdot[H]/d[AR] */
        J[1026] -= dqdci;             /* dwdot[O2]/d[AR] */
        J[1029] += dqdci;             /* dwdot[HO2]/d[AR] */
    }
    else {
        dqdc[0] = q_nocor;
        dqdc[1] = q_nocor + k_f*sc[3];
        dqdc[2] = q_nocor;
        dqdc[3] = TB[23][0]*q_nocor + k_f*sc[1];
        dqdc[4] = q_nocor;
        dqdc[5] = TB[23][1]*q_nocor;
        dqdc[6] = q_nocor - k_r;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = q_nocor;
        dqdc[12] = q_nocor;
        dqdc[13] = q_nocor;
        dqdc[14] = TB[23][2]*q_nocor;
        dqdc[15] = TB[23][3]*q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = q_nocor;
        dqdc[19] = q_nocor;
        dqdc[20] = q_nocor;
        dqdc[21] = q_nocor;
        dqdc[22] = q_nocor;
        dqdc[23] = q_nocor;
        dqdc[24] = q_nocor;
        dqdc[25] = q_nocor;
        dqdc[26] = TB[23][4]*q_nocor;
        dqdc[27] = q_nocor;
        dqdc[28] = q_nocor;
        dqdc[29] = q_nocor;
        dqdc[30] = TB[23][5]*q_nocor;
        dqdc[31] = TB[23][6]*q_nocor;
        for (int k=0; k<32; k++) {
            J[33*k+1] -= dqdc[k];
            J[33*k+3] -= dqdc[k];
            J[33*k+6] += dqdc[k];
        }
    }
    J[1057] -= dqdT; /* dwdot[H]/dT */
    J[1059] -= dqdT; /* dwdot[O2]/dT */
    J[1062] += dqdT; /* dwdot[HO2]/dT */

    /*reaction 25: 2.000000 H + M <=> H2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[24][0] - 1)*sc[0] + (TB[24][1] - 1)*sc[5] + (TB[24][2] - 1)*sc[13] + (TB[24][3] - 1)*sc[15] + (TB[24][4] - 1)*sc[26] + (TB[24][5] - 1)*sc[31];
    /* forward */
    phi_f = pow(sc[1], 2.000000);
    k_f = prefactor_units[24] * fwd_A[24]
                * exp(fwd_beta[24] * tc[0] - activation_units[24] * fwd_Ea[24] * invT);
    dlnkfdT = fwd_beta[24] * invT + activation_units[24] * fwd_Ea[24] * invT2;
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
        dqdci = (TB[24][0] - 1)*q_nocor - k_r;
        J[0] += dqdci;                /* dwdot[H2]/d[H2] */
        J[1] += -2 * dqdci;           /* dwdot[H]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*2.000000*sc[1];
        J[33] += dqdci;               /* dwdot[H2]/d[H] */
        J[34] += -2 * dqdci;          /* dwdot[H]/d[H] */
        /* d()/d[H2O] */
        dqdci = (TB[24][1] - 1)*q_nocor;
        J[165] += dqdci;              /* dwdot[H2]/d[H2O] */
        J[166] += -2 * dqdci;         /* dwdot[H]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (TB[24][2] - 1)*q_nocor;
        J[429] += dqdci;              /* dwdot[H2]/d[CH4] */
        J[430] += -2 * dqdci;         /* dwdot[H]/d[CH4] */
        /* d()/d[CO2] */
        dqdci = (TB[24][3] - 1)*q_nocor;
        J[495] += dqdci;              /* dwdot[H2]/d[CO2] */
        J[496] += -2 * dqdci;         /* dwdot[H]/d[CO2] */
        /* d()/d[C2H6] */
        dqdci = (TB[24][4] - 1)*q_nocor;
        J[858] += dqdci;              /* dwdot[H2]/d[C2H6] */
        J[859] += -2 * dqdci;         /* dwdot[H]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (TB[24][5] - 1)*q_nocor;
        J[1023] += dqdci;             /* dwdot[H2]/d[AR] */
        J[1024] += -2 * dqdci;        /* dwdot[H]/d[AR] */
    }
    else {
        dqdc[0] = TB[24][0]*q_nocor - k_r;
        dqdc[1] = q_nocor + k_f*2.000000*sc[1];
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor;
        dqdc[5] = TB[24][1]*q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = q_nocor;
        dqdc[12] = q_nocor;
        dqdc[13] = TB[24][2]*q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = TB[24][3]*q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = q_nocor;
        dqdc[19] = q_nocor;
        dqdc[20] = q_nocor;
        dqdc[21] = q_nocor;
        dqdc[22] = q_nocor;
        dqdc[23] = q_nocor;
        dqdc[24] = q_nocor;
        dqdc[25] = q_nocor;
        dqdc[26] = TB[24][4]*q_nocor;
        dqdc[27] = q_nocor;
        dqdc[28] = q_nocor;
        dqdc[29] = q_nocor;
        dqdc[30] = q_nocor;
        dqdc[31] = TB[24][5]*q_nocor;
        for (int k=0; k<32; k++) {
            J[33*k+0] += dqdc[k];
            J[33*k+1] += -2 * dqdc[k];
        }
    }
    J[1056] += dqdT; /* dwdot[H2]/dT */
    J[1057] += -2 * dqdT; /* dwdot[H]/dT */

    /*reaction 26: H + OH + M <=> H2O + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[25][0] - 1)*sc[0] + (TB[25][1] - 1)*sc[5] + (TB[25][2] - 1)*sc[13] + (TB[25][3] - 1)*sc[26] + (TB[25][4] - 1)*sc[31];
    /* forward */
    phi_f = sc[1]*sc[4];
    k_f = prefactor_units[25] * fwd_A[25]
                * exp(fwd_beta[25] * tc[0] - activation_units[25] * fwd_Ea[25] * invT);
    dlnkfdT = fwd_beta[25] * invT + activation_units[25] * fwd_Ea[25] * invT2;
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
        dqdci = (TB[25][0] - 1)*q_nocor;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[4] -= dqdci;                /* dwdot[OH]/d[H2] */
        J[5] += dqdci;                /* dwdot[H2O]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[4];
        J[34] -= dqdci;               /* dwdot[H]/d[H] */
        J[37] -= dqdci;               /* dwdot[OH]/d[H] */
        J[38] += dqdci;               /* dwdot[H2O]/d[H] */
        /* d()/d[OH] */
        dqdci =  + k_f*sc[1];
        J[133] -= dqdci;              /* dwdot[H]/d[OH] */
        J[136] -= dqdci;              /* dwdot[OH]/d[OH] */
        J[137] += dqdci;              /* dwdot[H2O]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (TB[25][1] - 1)*q_nocor - k_r;
        J[166] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[169] -= dqdci;              /* dwdot[OH]/d[H2O] */
        J[170] += dqdci;              /* dwdot[H2O]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (TB[25][2] - 1)*q_nocor;
        J[430] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[433] -= dqdci;              /* dwdot[OH]/d[CH4] */
        J[434] += dqdci;              /* dwdot[H2O]/d[CH4] */
        /* d()/d[C2H6] */
        dqdci = (TB[25][3] - 1)*q_nocor;
        J[859] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[862] -= dqdci;              /* dwdot[OH]/d[C2H6] */
        J[863] += dqdci;              /* dwdot[H2O]/d[C2H6] */
        /* d()/d[AR] */
        dqdci = (TB[25][4] - 1)*q_nocor;
        J[1024] -= dqdci;             /* dwdot[H]/d[AR] */
        J[1027] -= dqdci;             /* dwdot[OH]/d[AR] */
        J[1028] += dqdci;             /* dwdot[H2O]/d[AR] */
    }
    else {
        dqdc[0] = TB[25][0]*q_nocor;
        dqdc[1] = q_nocor + k_f*sc[4];
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor + k_f*sc[1];
        dqdc[5] = TB[25][1]*q_nocor - k_r;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = q_nocor;
        dqdc[12] = q_nocor;
        dqdc[13] = TB[25][2]*q_nocor;
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
        dqdc[25] = q_nocor;
        dqdc[26] = TB[25][3]*q_nocor;
        dqdc[27] = q_nocor;
        dqdc[28] = q_nocor;
        dqdc[29] = q_nocor;
        dqdc[30] = q_nocor;
        dqdc[31] = TB[25][4]*q_nocor;
        for (int k=0; k<32; k++) {
            J[33*k+1] -= dqdc[k];
            J[33*k+4] -= dqdc[k];
            J[33*k+5] += dqdc[k];
        }
    }
    J[1057] -= dqdT; /* dwdot[H]/dT */
    J[1060] -= dqdT; /* dwdot[OH]/dT */
    J[1061] += dqdT; /* dwdot[H2O]/dT */

    /*reaction 27: HCO + M <=> H + CO + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[26][0] - 1)*sc[0] + (TB[26][1] - 1)*sc[5] + (TB[26][2] - 1)*sc[13] + (TB[26][3] - 1)*sc[14] + (TB[26][4] - 1)*sc[15] + (TB[26][5] - 1)*sc[26];
    /* forward */
    phi_f = sc[16];
    k_f = prefactor_units[26] * fwd_A[26]
                * exp(fwd_beta[26] * tc[0] - activation_units[26] * fwd_Ea[26] * invT);
    dlnkfdT = fwd_beta[26] * invT + activation_units[26] * fwd_Ea[26] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[14];
    Kc = refC * exp(-g_RT[1] - g_RT[14] + g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[16]) + (h_RT[1] + h_RT[14]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[14] += q; /* CO */
    wdot[16] -= q; /* HCO */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[26][0] - 1)*q_nocor;
        J[1] += dqdci;                /* dwdot[H]/d[H2] */
        J[14] += dqdci;               /* dwdot[CO]/d[H2] */
        J[16] -= dqdci;               /* dwdot[HCO]/d[H2] */
        /* d()/d[H] */
        dqdci =  - k_r*sc[14];
        J[34] += dqdci;               /* dwdot[H]/d[H] */
        J[47] += dqdci;               /* dwdot[CO]/d[H] */
        J[49] -= dqdci;               /* dwdot[HCO]/d[H] */
        /* d()/d[H2O] */
        dqdci = (TB[26][1] - 1)*q_nocor;
        J[166] += dqdci;              /* dwdot[H]/d[H2O] */
        J[179] += dqdci;              /* dwdot[CO]/d[H2O] */
        J[181] -= dqdci;              /* dwdot[HCO]/d[H2O] */
        /* d()/d[CH4] */
        dqdci = (TB[26][2] - 1)*q_nocor;
        J[430] += dqdci;              /* dwdot[H]/d[CH4] */
        J[443] += dqdci;              /* dwdot[CO]/d[CH4] */
        J[445] -= dqdci;              /* dwdot[HCO]/d[CH4] */
        /* d()/d[CO] */
        dqdci = (TB[26][3] - 1)*q_nocor - k_r*sc[1];
        J[463] += dqdci;              /* dwdot[H]/d[CO] */
        J[476] += dqdci;              /* dwdot[CO]/d[CO] */
        J[478] -= dqdci;              /* dwdot[HCO]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[26][4] - 1)*q_nocor;
        J[496] += dqdci;              /* dwdot[H]/d[CO2] */
        J[509] += dqdci;              /* dwdot[CO]/d[CO2] */
        J[511] -= dqdci;              /* dwdot[HCO]/d[CO2] */
        /* d()/d[HCO] */
        dqdci =  + k_f;
        J[529] += dqdci;              /* dwdot[H]/d[HCO] */
        J[542] += dqdci;              /* dwdot[CO]/d[HCO] */
        J[544] -= dqdci;              /* dwdot[HCO]/d[HCO] */
        /* d()/d[C2H6] */
        dqdci = (TB[26][5] - 1)*q_nocor;
        J[859] += dqdci;              /* dwdot[H]/d[C2H6] */
        J[872] += dqdci;              /* dwdot[CO]/d[C2H6] */
        J[874] -= dqdci;              /* dwdot[HCO]/d[C2H6] */
    }
    else {
        dqdc[0] = TB[26][0]*q_nocor;
        dqdc[1] = q_nocor - k_r*sc[14];
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor;
        dqdc[5] = TB[26][1]*q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = q_nocor;
        dqdc[12] = q_nocor;
        dqdc[13] = TB[26][2]*q_nocor;
        dqdc[14] = TB[26][3]*q_nocor - k_r*sc[1];
        dqdc[15] = TB[26][4]*q_nocor;
        dqdc[16] = q_nocor + k_f;
        dqdc[17] = q_nocor;
        dqdc[18] = q_nocor;
        dqdc[19] = q_nocor;
        dqdc[20] = q_nocor;
        dqdc[21] = q_nocor;
        dqdc[22] = q_nocor;
        dqdc[23] = q_nocor;
        dqdc[24] = q_nocor;
        dqdc[25] = q_nocor;
        dqdc[26] = TB[26][5]*q_nocor;
        dqdc[27] = q_nocor;
        dqdc[28] = q_nocor;
        dqdc[29] = q_nocor;
        dqdc[30] = q_nocor;
        dqdc[31] = q_nocor;
        for (int k=0; k<32; k++) {
            J[33*k+1] += dqdc[k];
            J[33*k+14] += dqdc[k];
            J[33*k+16] -= dqdc[k];
        }
    }
    J[1057] += dqdT; /* dwdot[H]/dT */
    J[1070] += dqdT; /* dwdot[CO]/dT */
    J[1072] -= dqdT; /* dwdot[HCO]/dT */

    /*reaction 28: O + H2 <=> H + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[2];
    k_f = prefactor_units[27] * fwd_A[27]
                * exp(fwd_beta[27] * tc[0] - activation_units[27] * fwd_Ea[27] * invT);
    dlnkfdT = fwd_beta[27] * invT + activation_units[27] * fwd_Ea[27] * invT2;
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
    J[33] -= dqdci;               /* dwdot[H2]/d[H] */
    J[34] += dqdci;               /* dwdot[H]/d[H] */
    J[35] -= dqdci;               /* dwdot[O]/d[H] */
    J[37] += dqdci;               /* dwdot[OH]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[0];
    J[66] -= dqdci;               /* dwdot[H2]/d[O] */
    J[67] += dqdci;               /* dwdot[H]/d[O] */
    J[68] -= dqdci;               /* dwdot[O]/d[O] */
    J[70] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[1];
    J[132] -= dqdci;              /* dwdot[H2]/d[OH] */
    J[133] += dqdci;              /* dwdot[H]/d[OH] */
    J[134] -= dqdci;              /* dwdot[O]/d[OH] */
    J[136] += dqdci;              /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[1056] -= dqdT;              /* dwdot[H2]/dT */
    J[1057] += dqdT;              /* dwdot[H]/dT */
    J[1058] -= dqdT;              /* dwdot[O]/dT */
    J[1060] += dqdT;              /* dwdot[OH]/dT */

    /*reaction 29: O + HO2 <=> OH + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[6];
    k_f = prefactor_units[28] * fwd_A[28]
                * exp(fwd_beta[28] * tc[0] - activation_units[28] * fwd_Ea[28] * invT);
    dlnkfdT = fwd_beta[28] * invT + activation_units[28] * fwd_Ea[28] * invT2;
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
    J[68] -= dqdci;               /* dwdot[O]/d[O] */
    J[69] += dqdci;               /* dwdot[O2]/d[O] */
    J[70] += dqdci;               /* dwdot[OH]/d[O] */
    J[72] -= dqdci;               /* dwdot[HO2]/d[O] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[4];
    J[101] -= dqdci;              /* dwdot[O]/d[O2] */
    J[102] += dqdci;              /* dwdot[O2]/d[O2] */
    J[103] += dqdci;              /* dwdot[OH]/d[O2] */
    J[105] -= dqdci;              /* dwdot[HO2]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[3];
    J[134] -= dqdci;              /* dwdot[O]/d[OH] */
    J[135] += dqdci;              /* dwdot[O2]/d[OH] */
    J[136] += dqdci;              /* dwdot[OH]/d[OH] */
    J[138] -= dqdci;              /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[2];
    J[200] -= dqdci;              /* dwdot[O]/d[HO2] */
    J[201] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[202] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[204] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1058] -= dqdT;              /* dwdot[O]/dT */
    J[1059] += dqdT;              /* dwdot[O2]/dT */
    J[1060] += dqdT;              /* dwdot[OH]/dT */
    J[1062] -= dqdT;              /* dwdot[HO2]/dT */

    /*reaction 30: O + H2O2 <=> OH + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[7];
    k_f = prefactor_units[29] * fwd_A[29]
                * exp(fwd_beta[29] * tc[0] - activation_units[29] * fwd_Ea[29] * invT);
    dlnkfdT = fwd_beta[29] * invT + activation_units[29] * fwd_Ea[29] * invT2;
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
    J[68] -= dqdci;               /* dwdot[O]/d[O] */
    J[70] += dqdci;               /* dwdot[OH]/d[O] */
    J[72] += dqdci;               /* dwdot[HO2]/d[O] */
    J[73] -= dqdci;               /* dwdot[H2O2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[6];
    J[134] -= dqdci;              /* dwdot[O]/d[OH] */
    J[136] += dqdci;              /* dwdot[OH]/d[OH] */
    J[138] += dqdci;              /* dwdot[HO2]/d[OH] */
    J[139] -= dqdci;              /* dwdot[H2O2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[4];
    J[200] -= dqdci;              /* dwdot[O]/d[HO2] */
    J[202] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[204] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[205] -= dqdci;              /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[2];
    J[233] -= dqdci;              /* dwdot[O]/d[H2O2] */
    J[235] += dqdci;              /* dwdot[OH]/d[H2O2] */
    J[237] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[238] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[1058] -= dqdT;              /* dwdot[O]/dT */
    J[1060] += dqdT;              /* dwdot[OH]/dT */
    J[1062] += dqdT;              /* dwdot[HO2]/dT */
    J[1063] -= dqdT;              /* dwdot[H2O2]/dT */

    /*reaction 31: O + CH <=> H + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[9];
    k_f = prefactor_units[30] * fwd_A[30]
                * exp(fwd_beta[30] * tc[0] - activation_units[30] * fwd_Ea[30] * invT);
    dlnkfdT = fwd_beta[30] * invT + activation_units[30] * fwd_Ea[30] * invT2;
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
    wdot[9] -= q; /* CH */
    wdot[14] += q; /* CO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[14];
    J[34] += dqdci;               /* dwdot[H]/d[H] */
    J[35] -= dqdci;               /* dwdot[O]/d[H] */
    J[42] -= dqdci;               /* dwdot[CH]/d[H] */
    J[47] += dqdci;               /* dwdot[CO]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[9];
    J[67] += dqdci;               /* dwdot[H]/d[O] */
    J[68] -= dqdci;               /* dwdot[O]/d[O] */
    J[75] -= dqdci;               /* dwdot[CH]/d[O] */
    J[80] += dqdci;               /* dwdot[CO]/d[O] */
    /* d()/d[CH] */
    dqdci =  + k_f*sc[2];
    J[298] += dqdci;              /* dwdot[H]/d[CH] */
    J[299] -= dqdci;              /* dwdot[O]/d[CH] */
    J[306] -= dqdci;              /* dwdot[CH]/d[CH] */
    J[311] += dqdci;              /* dwdot[CO]/d[CH] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[1];
    J[463] += dqdci;              /* dwdot[H]/d[CO] */
    J[464] -= dqdci;              /* dwdot[O]/d[CO] */
    J[471] -= dqdci;              /* dwdot[CH]/d[CO] */
    J[476] += dqdci;              /* dwdot[CO]/d[CO] */
    /* d()/dT */
    J[1057] += dqdT;              /* dwdot[H]/dT */
    J[1058] -= dqdT;              /* dwdot[O]/dT */
    J[1065] -= dqdT;              /* dwdot[CH]/dT */
    J[1070] += dqdT;              /* dwdot[CO]/dT */

    /*reaction 32: O + CH2 <=> H + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[10];
    k_f = prefactor_units[31] * fwd_A[31]
                * exp(fwd_beta[31] * tc[0] - activation_units[31] * fwd_Ea[31] * invT);
    dlnkfdT = fwd_beta[31] * invT + activation_units[31] * fwd_Ea[31] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[16];
    Kc = exp(-g_RT[1] + g_RT[2] + g_RT[10] - g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[10]) + (h_RT[1] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[2] -= q; /* O */
    wdot[10] -= q; /* CH2 */
    wdot[16] += q; /* HCO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[16];
    J[34] += dqdci;               /* dwdot[H]/d[H] */
    J[35] -= dqdci;               /* dwdot[O]/d[H] */
    J[43] -= dqdci;               /* dwdot[CH2]/d[H] */
    J[49] += dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[10];
    J[67] += dqdci;               /* dwdot[H]/d[O] */
    J[68] -= dqdci;               /* dwdot[O]/d[O] */
    J[76] -= dqdci;               /* dwdot[CH2]/d[O] */
    J[82] += dqdci;               /* dwdot[HCO]/d[O] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[2];
    J[331] += dqdci;              /* dwdot[H]/d[CH2] */
    J[332] -= dqdci;              /* dwdot[O]/d[CH2] */
    J[340] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[346] += dqdci;              /* dwdot[HCO]/d[CH2] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[1];
    J[529] += dqdci;              /* dwdot[H]/d[HCO] */
    J[530] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[538] -= dqdci;              /* dwdot[CH2]/d[HCO] */
    J[544] += dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[1057] += dqdT;              /* dwdot[H]/dT */
    J[1058] -= dqdT;              /* dwdot[O]/dT */
    J[1066] -= dqdT;              /* dwdot[CH2]/dT */
    J[1072] += dqdT;              /* dwdot[HCO]/dT */

    /*reaction 33: O + CH2(S) <=> H2 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[11];
    k_f = prefactor_units[32] * fwd_A[32]
                * exp(fwd_beta[32] * tc[0] - activation_units[32] * fwd_Ea[32] * invT);
    dlnkfdT = fwd_beta[32] * invT + activation_units[32] * fwd_Ea[32] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[14];
    Kc = exp(-g_RT[0] + g_RT[2] + g_RT[11] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[11]) + (h_RT[0] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[2] -= q; /* O */
    wdot[11] -= q; /* CH2(S) */
    wdot[14] += q; /* CO */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[14];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[2] -= dqdci;                /* dwdot[O]/d[H2] */
    J[11] -= dqdci;               /* dwdot[CH2(S)]/d[H2] */
    J[14] += dqdci;               /* dwdot[CO]/d[H2] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[11];
    J[66] += dqdci;               /* dwdot[H2]/d[O] */
    J[68] -= dqdci;               /* dwdot[O]/d[O] */
    J[77] -= dqdci;               /* dwdot[CH2(S)]/d[O] */
    J[80] += dqdci;               /* dwdot[CO]/d[O] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[2];
    J[363] += dqdci;              /* dwdot[H2]/d[CH2(S)] */
    J[365] -= dqdci;              /* dwdot[O]/d[CH2(S)] */
    J[374] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[377] += dqdci;              /* dwdot[CO]/d[CH2(S)] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[0];
    J[462] += dqdci;              /* dwdot[H2]/d[CO] */
    J[464] -= dqdci;              /* dwdot[O]/d[CO] */
    J[473] -= dqdci;              /* dwdot[CH2(S)]/d[CO] */
    J[476] += dqdci;              /* dwdot[CO]/d[CO] */
    /* d()/dT */
    J[1056] += dqdT;              /* dwdot[H2]/dT */
    J[1058] -= dqdT;              /* dwdot[O]/dT */
    J[1067] -= dqdT;              /* dwdot[CH2(S)]/dT */
    J[1070] += dqdT;              /* dwdot[CO]/dT */

    /*reaction 34: O + CH2(S) <=> H + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[11];
    k_f = prefactor_units[33] * fwd_A[33]
                * exp(fwd_beta[33] * tc[0] - activation_units[33] * fwd_Ea[33] * invT);
    dlnkfdT = fwd_beta[33] * invT + activation_units[33] * fwd_Ea[33] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[16];
    Kc = exp(-g_RT[1] + g_RT[2] + g_RT[11] - g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[11]) + (h_RT[1] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[2] -= q; /* O */
    wdot[11] -= q; /* CH2(S) */
    wdot[16] += q; /* HCO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[16];
    J[34] += dqdci;               /* dwdot[H]/d[H] */
    J[35] -= dqdci;               /* dwdot[O]/d[H] */
    J[44] -= dqdci;               /* dwdot[CH2(S)]/d[H] */
    J[49] += dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[11];
    J[67] += dqdci;               /* dwdot[H]/d[O] */
    J[68] -= dqdci;               /* dwdot[O]/d[O] */
    J[77] -= dqdci;               /* dwdot[CH2(S)]/d[O] */
    J[82] += dqdci;               /* dwdot[HCO]/d[O] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[2];
    J[364] += dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[365] -= dqdci;              /* dwdot[O]/d[CH2(S)] */
    J[374] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[379] += dqdci;              /* dwdot[HCO]/d[CH2(S)] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[1];
    J[529] += dqdci;              /* dwdot[H]/d[HCO] */
    J[530] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[539] -= dqdci;              /* dwdot[CH2(S)]/d[HCO] */
    J[544] += dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[1057] += dqdT;              /* dwdot[H]/dT */
    J[1058] -= dqdT;              /* dwdot[O]/dT */
    J[1067] -= dqdT;              /* dwdot[CH2(S)]/dT */
    J[1072] += dqdT;              /* dwdot[HCO]/dT */

    /*reaction 35: O + CH3 <=> H + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[12];
    k_f = prefactor_units[34] * fwd_A[34]
                * exp(fwd_beta[34] * tc[0] - activation_units[34] * fwd_Ea[34] * invT);
    dlnkfdT = fwd_beta[34] * invT + activation_units[34] * fwd_Ea[34] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[17];
    Kc = exp(-g_RT[1] + g_RT[2] + g_RT[12] - g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[12]) + (h_RT[1] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[2] -= q; /* O */
    wdot[12] -= q; /* CH3 */
    wdot[17] += q; /* CH2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[17];
    J[34] += dqdci;               /* dwdot[H]/d[H] */
    J[35] -= dqdci;               /* dwdot[O]/d[H] */
    J[45] -= dqdci;               /* dwdot[CH3]/d[H] */
    J[50] += dqdci;               /* dwdot[CH2O]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[12];
    J[67] += dqdci;               /* dwdot[H]/d[O] */
    J[68] -= dqdci;               /* dwdot[O]/d[O] */
    J[78] -= dqdci;               /* dwdot[CH3]/d[O] */
    J[83] += dqdci;               /* dwdot[CH2O]/d[O] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[2];
    J[397] += dqdci;              /* dwdot[H]/d[CH3] */
    J[398] -= dqdci;              /* dwdot[O]/d[CH3] */
    J[408] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[413] += dqdci;              /* dwdot[CH2O]/d[CH3] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[1];
    J[562] += dqdci;              /* dwdot[H]/d[CH2O] */
    J[563] -= dqdci;              /* dwdot[O]/d[CH2O] */
    J[573] -= dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[578] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[1057] += dqdT;              /* dwdot[H]/dT */
    J[1058] -= dqdT;              /* dwdot[O]/dT */
    J[1068] -= dqdT;              /* dwdot[CH3]/dT */
    J[1073] += dqdT;              /* dwdot[CH2O]/dT */

    /*reaction 36: O + CH4 <=> OH + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[13];
    k_f = prefactor_units[35] * fwd_A[35]
                * exp(fwd_beta[35] * tc[0] - activation_units[35] * fwd_Ea[35] * invT);
    dlnkfdT = fwd_beta[35] * invT + activation_units[35] * fwd_Ea[35] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[12];
    Kc = exp(g_RT[2] - g_RT[4] - g_RT[12] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[13]) + (h_RT[4] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    wdot[12] += q; /* CH3 */
    wdot[13] -= q; /* CH4 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[13];
    J[68] -= dqdci;               /* dwdot[O]/d[O] */
    J[70] += dqdci;               /* dwdot[OH]/d[O] */
    J[78] += dqdci;               /* dwdot[CH3]/d[O] */
    J[79] -= dqdci;               /* dwdot[CH4]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[12];
    J[134] -= dqdci;              /* dwdot[O]/d[OH] */
    J[136] += dqdci;              /* dwdot[OH]/d[OH] */
    J[144] += dqdci;              /* dwdot[CH3]/d[OH] */
    J[145] -= dqdci;              /* dwdot[CH4]/d[OH] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[4];
    J[398] -= dqdci;              /* dwdot[O]/d[CH3] */
    J[400] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[408] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[409] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[2];
    J[431] -= dqdci;              /* dwdot[O]/d[CH4] */
    J[433] += dqdci;              /* dwdot[OH]/d[CH4] */
    J[441] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[442] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[1058] -= dqdT;              /* dwdot[O]/dT */
    J[1060] += dqdT;              /* dwdot[OH]/dT */
    J[1068] += dqdT;              /* dwdot[CH3]/dT */
    J[1069] -= dqdT;              /* dwdot[CH4]/dT */

    /*reaction 37: O + HCO <=> OH + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[16];
    k_f = prefactor_units[36] * fwd_A[36]
                * exp(fwd_beta[36] * tc[0] - activation_units[36] * fwd_Ea[36] * invT);
    dlnkfdT = fwd_beta[36] * invT + activation_units[36] * fwd_Ea[36] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[14];
    Kc = exp(g_RT[2] - g_RT[4] - g_RT[14] + g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[16]) + (h_RT[4] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    wdot[14] += q; /* CO */
    wdot[16] -= q; /* HCO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[16];
    J[68] -= dqdci;               /* dwdot[O]/d[O] */
    J[70] += dqdci;               /* dwdot[OH]/d[O] */
    J[80] += dqdci;               /* dwdot[CO]/d[O] */
    J[82] -= dqdci;               /* dwdot[HCO]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[14];
    J[134] -= dqdci;              /* dwdot[O]/d[OH] */
    J[136] += dqdci;              /* dwdot[OH]/d[OH] */
    J[146] += dqdci;              /* dwdot[CO]/d[OH] */
    J[148] -= dqdci;              /* dwdot[HCO]/d[OH] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[4];
    J[464] -= dqdci;              /* dwdot[O]/d[CO] */
    J[466] += dqdci;              /* dwdot[OH]/d[CO] */
    J[476] += dqdci;              /* dwdot[CO]/d[CO] */
    J[478] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[2];
    J[530] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[532] += dqdci;              /* dwdot[OH]/d[HCO] */
    J[542] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[544] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[1058] -= dqdT;              /* dwdot[O]/dT */
    J[1060] += dqdT;              /* dwdot[OH]/dT */
    J[1070] += dqdT;              /* dwdot[CO]/dT */
    J[1072] -= dqdT;              /* dwdot[HCO]/dT */

    /*reaction 38: O + HCO <=> H + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[16];
    k_f = prefactor_units[37] * fwd_A[37]
                * exp(fwd_beta[37] * tc[0] - activation_units[37] * fwd_Ea[37] * invT);
    dlnkfdT = fwd_beta[37] * invT + activation_units[37] * fwd_Ea[37] * invT2;
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
    wdot[1] += q; /* H */
    wdot[2] -= q; /* O */
    wdot[15] += q; /* CO2 */
    wdot[16] -= q; /* HCO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[15];
    J[34] += dqdci;               /* dwdot[H]/d[H] */
    J[35] -= dqdci;               /* dwdot[O]/d[H] */
    J[48] += dqdci;               /* dwdot[CO2]/d[H] */
    J[49] -= dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[16];
    J[67] += dqdci;               /* dwdot[H]/d[O] */
    J[68] -= dqdci;               /* dwdot[O]/d[O] */
    J[81] += dqdci;               /* dwdot[CO2]/d[O] */
    J[82] -= dqdci;               /* dwdot[HCO]/d[O] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[1];
    J[496] += dqdci;              /* dwdot[H]/d[CO2] */
    J[497] -= dqdci;              /* dwdot[O]/d[CO2] */
    J[510] += dqdci;              /* dwdot[CO2]/d[CO2] */
    J[511] -= dqdci;              /* dwdot[HCO]/d[CO2] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[2];
    J[529] += dqdci;              /* dwdot[H]/d[HCO] */
    J[530] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[543] += dqdci;              /* dwdot[CO2]/d[HCO] */
    J[544] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[1057] += dqdT;              /* dwdot[H]/dT */
    J[1058] -= dqdT;              /* dwdot[O]/dT */
    J[1071] += dqdT;              /* dwdot[CO2]/dT */
    J[1072] -= dqdT;              /* dwdot[HCO]/dT */

    /*reaction 39: O + CH2O <=> OH + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[17];
    k_f = prefactor_units[38] * fwd_A[38]
                * exp(fwd_beta[38] * tc[0] - activation_units[38] * fwd_Ea[38] * invT);
    dlnkfdT = fwd_beta[38] * invT + activation_units[38] * fwd_Ea[38] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[16];
    Kc = exp(g_RT[2] - g_RT[4] - g_RT[16] + g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[17]) + (h_RT[4] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    wdot[16] += q; /* HCO */
    wdot[17] -= q; /* CH2O */
    /* d()/d[O] */
    dqdci =  + k_f*sc[17];
    J[68] -= dqdci;               /* dwdot[O]/d[O] */
    J[70] += dqdci;               /* dwdot[OH]/d[O] */
    J[82] += dqdci;               /* dwdot[HCO]/d[O] */
    J[83] -= dqdci;               /* dwdot[CH2O]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[16];
    J[134] -= dqdci;              /* dwdot[O]/d[OH] */
    J[136] += dqdci;              /* dwdot[OH]/d[OH] */
    J[148] += dqdci;              /* dwdot[HCO]/d[OH] */
    J[149] -= dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[4];
    J[530] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[532] += dqdci;              /* dwdot[OH]/d[HCO] */
    J[544] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[545] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[2];
    J[563] -= dqdci;              /* dwdot[O]/d[CH2O] */
    J[565] += dqdci;              /* dwdot[OH]/d[CH2O] */
    J[577] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[578] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[1058] -= dqdT;              /* dwdot[O]/dT */
    J[1060] += dqdT;              /* dwdot[OH]/dT */
    J[1072] += dqdT;              /* dwdot[HCO]/dT */
    J[1073] -= dqdT;              /* dwdot[CH2O]/dT */

    /*reaction 40: O + CH2OH <=> OH + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[18];
    k_f = prefactor_units[39] * fwd_A[39]
                * exp(fwd_beta[39] * tc[0] - activation_units[39] * fwd_Ea[39] * invT);
    dlnkfdT = fwd_beta[39] * invT + activation_units[39] * fwd_Ea[39] * invT2;
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
    wdot[17] += q; /* CH2O */
    wdot[18] -= q; /* CH2OH */
    /* d()/d[O] */
    dqdci =  + k_f*sc[18];
    J[68] -= dqdci;               /* dwdot[O]/d[O] */
    J[70] += dqdci;               /* dwdot[OH]/d[O] */
    J[83] += dqdci;               /* dwdot[CH2O]/d[O] */
    J[84] -= dqdci;               /* dwdot[CH2OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[17];
    J[134] -= dqdci;              /* dwdot[O]/d[OH] */
    J[136] += dqdci;              /* dwdot[OH]/d[OH] */
    J[149] += dqdci;              /* dwdot[CH2O]/d[OH] */
    J[150] -= dqdci;              /* dwdot[CH2OH]/d[OH] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[4];
    J[563] -= dqdci;              /* dwdot[O]/d[CH2O] */
    J[565] += dqdci;              /* dwdot[OH]/d[CH2O] */
    J[578] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[579] -= dqdci;              /* dwdot[CH2OH]/d[CH2O] */
    /* d()/d[CH2OH] */
    dqdci =  + k_f*sc[2];
    J[596] -= dqdci;              /* dwdot[O]/d[CH2OH] */
    J[598] += dqdci;              /* dwdot[OH]/d[CH2OH] */
    J[611] += dqdci;              /* dwdot[CH2O]/d[CH2OH] */
    J[612] -= dqdci;              /* dwdot[CH2OH]/d[CH2OH] */
    /* d()/dT */
    J[1058] -= dqdT;              /* dwdot[O]/dT */
    J[1060] += dqdT;              /* dwdot[OH]/dT */
    J[1073] += dqdT;              /* dwdot[CH2O]/dT */
    J[1074] -= dqdT;              /* dwdot[CH2OH]/dT */

    /*reaction 41: O + CH3O <=> OH + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[19];
    k_f = prefactor_units[40] * fwd_A[40]
                * exp(fwd_beta[40] * tc[0] - activation_units[40] * fwd_Ea[40] * invT);
    dlnkfdT = fwd_beta[40] * invT + activation_units[40] * fwd_Ea[40] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[17];
    Kc = exp(g_RT[2] - g_RT[4] - g_RT[17] + g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[19]) + (h_RT[4] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    wdot[17] += q; /* CH2O */
    wdot[19] -= q; /* CH3O */
    /* d()/d[O] */
    dqdci =  + k_f*sc[19];
    J[68] -= dqdci;               /* dwdot[O]/d[O] */
    J[70] += dqdci;               /* dwdot[OH]/d[O] */
    J[83] += dqdci;               /* dwdot[CH2O]/d[O] */
    J[85] -= dqdci;               /* dwdot[CH3O]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[17];
    J[134] -= dqdci;              /* dwdot[O]/d[OH] */
    J[136] += dqdci;              /* dwdot[OH]/d[OH] */
    J[149] += dqdci;              /* dwdot[CH2O]/d[OH] */
    J[151] -= dqdci;              /* dwdot[CH3O]/d[OH] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[4];
    J[563] -= dqdci;              /* dwdot[O]/d[CH2O] */
    J[565] += dqdci;              /* dwdot[OH]/d[CH2O] */
    J[578] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[580] -= dqdci;              /* dwdot[CH3O]/d[CH2O] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[2];
    J[629] -= dqdci;              /* dwdot[O]/d[CH3O] */
    J[631] += dqdci;              /* dwdot[OH]/d[CH3O] */
    J[644] += dqdci;              /* dwdot[CH2O]/d[CH3O] */
    J[646] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[1058] -= dqdT;              /* dwdot[O]/dT */
    J[1060] += dqdT;              /* dwdot[OH]/dT */
    J[1073] += dqdT;              /* dwdot[CH2O]/dT */
    J[1075] -= dqdT;              /* dwdot[CH3O]/dT */

    /*reaction 42: O + CH3OH <=> OH + CH2OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[20];
    k_f = prefactor_units[41] * fwd_A[41]
                * exp(fwd_beta[41] * tc[0] - activation_units[41] * fwd_Ea[41] * invT);
    dlnkfdT = fwd_beta[41] * invT + activation_units[41] * fwd_Ea[41] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[18];
    Kc = exp(g_RT[2] - g_RT[4] - g_RT[18] + g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[20]) + (h_RT[4] + h_RT[18]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    wdot[18] += q; /* CH2OH */
    wdot[20] -= q; /* CH3OH */
    /* d()/d[O] */
    dqdci =  + k_f*sc[20];
    J[68] -= dqdci;               /* dwdot[O]/d[O] */
    J[70] += dqdci;               /* dwdot[OH]/d[O] */
    J[84] += dqdci;               /* dwdot[CH2OH]/d[O] */
    J[86] -= dqdci;               /* dwdot[CH3OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[18];
    J[134] -= dqdci;              /* dwdot[O]/d[OH] */
    J[136] += dqdci;              /* dwdot[OH]/d[OH] */
    J[150] += dqdci;              /* dwdot[CH2OH]/d[OH] */
    J[152] -= dqdci;              /* dwdot[CH3OH]/d[OH] */
    /* d()/d[CH2OH] */
    dqdci =  - k_r*sc[4];
    J[596] -= dqdci;              /* dwdot[O]/d[CH2OH] */
    J[598] += dqdci;              /* dwdot[OH]/d[CH2OH] */
    J[612] += dqdci;              /* dwdot[CH2OH]/d[CH2OH] */
    J[614] -= dqdci;              /* dwdot[CH3OH]/d[CH2OH] */
    /* d()/d[CH3OH] */
    dqdci =  + k_f*sc[2];
    J[662] -= dqdci;              /* dwdot[O]/d[CH3OH] */
    J[664] += dqdci;              /* dwdot[OH]/d[CH3OH] */
    J[678] += dqdci;              /* dwdot[CH2OH]/d[CH3OH] */
    J[680] -= dqdci;              /* dwdot[CH3OH]/d[CH3OH] */
    /* d()/dT */
    J[1058] -= dqdT;              /* dwdot[O]/dT */
    J[1060] += dqdT;              /* dwdot[OH]/dT */
    J[1074] += dqdT;              /* dwdot[CH2OH]/dT */
    J[1076] -= dqdT;              /* dwdot[CH3OH]/dT */

    /*reaction 43: O + CH3OH <=> OH + CH3O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[20];
    k_f = prefactor_units[42] * fwd_A[42]
                * exp(fwd_beta[42] * tc[0] - activation_units[42] * fwd_Ea[42] * invT);
    dlnkfdT = fwd_beta[42] * invT + activation_units[42] * fwd_Ea[42] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[19];
    Kc = exp(g_RT[2] - g_RT[4] - g_RT[19] + g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[20]) + (h_RT[4] + h_RT[19]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    wdot[19] += q; /* CH3O */
    wdot[20] -= q; /* CH3OH */
    /* d()/d[O] */
    dqdci =  + k_f*sc[20];
    J[68] -= dqdci;               /* dwdot[O]/d[O] */
    J[70] += dqdci;               /* dwdot[OH]/d[O] */
    J[85] += dqdci;               /* dwdot[CH3O]/d[O] */
    J[86] -= dqdci;               /* dwdot[CH3OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[19];
    J[134] -= dqdci;              /* dwdot[O]/d[OH] */
    J[136] += dqdci;              /* dwdot[OH]/d[OH] */
    J[151] += dqdci;              /* dwdot[CH3O]/d[OH] */
    J[152] -= dqdci;              /* dwdot[CH3OH]/d[OH] */
    /* d()/d[CH3O] */
    dqdci =  - k_r*sc[4];
    J[629] -= dqdci;              /* dwdot[O]/d[CH3O] */
    J[631] += dqdci;              /* dwdot[OH]/d[CH3O] */
    J[646] += dqdci;              /* dwdot[CH3O]/d[CH3O] */
    J[647] -= dqdci;              /* dwdot[CH3OH]/d[CH3O] */
    /* d()/d[CH3OH] */
    dqdci =  + k_f*sc[2];
    J[662] -= dqdci;              /* dwdot[O]/d[CH3OH] */
    J[664] += dqdci;              /* dwdot[OH]/d[CH3OH] */
    J[679] += dqdci;              /* dwdot[CH3O]/d[CH3OH] */
    J[680] -= dqdci;              /* dwdot[CH3OH]/d[CH3OH] */
    /* d()/dT */
    J[1058] -= dqdT;              /* dwdot[O]/dT */
    J[1060] += dqdT;              /* dwdot[OH]/dT */
    J[1075] += dqdT;              /* dwdot[CH3O]/dT */
    J[1076] -= dqdT;              /* dwdot[CH3OH]/dT */

    /*reaction 44: O + C2H <=> CH + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[21];
    k_f = prefactor_units[43] * fwd_A[43]
                * exp(fwd_beta[43] * tc[0] - activation_units[43] * fwd_Ea[43] * invT);
    dlnkfdT = fwd_beta[43] * invT + activation_units[43] * fwd_Ea[43] * invT2;
    /* reverse */
    phi_r = sc[9]*sc[14];
    Kc = exp(g_RT[2] - g_RT[9] - g_RT[14] + g_RT[21]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[21]) + (h_RT[9] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[9] += q; /* CH */
    wdot[14] += q; /* CO */
    wdot[21] -= q; /* C2H */
    /* d()/d[O] */
    dqdci =  + k_f*sc[21];
    J[68] -= dqdci;               /* dwdot[O]/d[O] */
    J[75] += dqdci;               /* dwdot[CH]/d[O] */
    J[80] += dqdci;               /* dwdot[CO]/d[O] */
    J[87] -= dqdci;               /* dwdot[C2H]/d[O] */
    /* d()/d[CH] */
    dqdci =  - k_r*sc[14];
    J[299] -= dqdci;              /* dwdot[O]/d[CH] */
    J[306] += dqdci;              /* dwdot[CH]/d[CH] */
    J[311] += dqdci;              /* dwdot[CO]/d[CH] */
    J[318] -= dqdci;              /* dwdot[C2H]/d[CH] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[9];
    J[464] -= dqdci;              /* dwdot[O]/d[CO] */
    J[471] += dqdci;              /* dwdot[CH]/d[CO] */
    J[476] += dqdci;              /* dwdot[CO]/d[CO] */
    J[483] -= dqdci;              /* dwdot[C2H]/d[CO] */
    /* d()/d[C2H] */
    dqdci =  + k_f*sc[2];
    J[695] -= dqdci;              /* dwdot[O]/d[C2H] */
    J[702] += dqdci;              /* dwdot[CH]/d[C2H] */
    J[707] += dqdci;              /* dwdot[CO]/d[C2H] */
    J[714] -= dqdci;              /* dwdot[C2H]/d[C2H] */
    /* d()/dT */
    J[1058] -= dqdT;              /* dwdot[O]/dT */
    J[1065] += dqdT;              /* dwdot[CH]/dT */
    J[1070] += dqdT;              /* dwdot[CO]/dT */
    J[1077] -= dqdT;              /* dwdot[C2H]/dT */

    /*reaction 45: O + C2H2 <=> H + HCCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[22];
    k_f = prefactor_units[44] * fwd_A[44]
                * exp(fwd_beta[44] * tc[0] - activation_units[44] * fwd_Ea[44] * invT);
    dlnkfdT = fwd_beta[44] * invT + activation_units[44] * fwd_Ea[44] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[27];
    Kc = exp(-g_RT[1] + g_RT[2] + g_RT[22] - g_RT[27]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[22]) + (h_RT[1] + h_RT[27]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[2] -= q; /* O */
    wdot[22] -= q; /* C2H2 */
    wdot[27] += q; /* HCCO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[27];
    J[34] += dqdci;               /* dwdot[H]/d[H] */
    J[35] -= dqdci;               /* dwdot[O]/d[H] */
    J[55] -= dqdci;               /* dwdot[C2H2]/d[H] */
    J[60] += dqdci;               /* dwdot[HCCO]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[22];
    J[67] += dqdci;               /* dwdot[H]/d[O] */
    J[68] -= dqdci;               /* dwdot[O]/d[O] */
    J[88] -= dqdci;               /* dwdot[C2H2]/d[O] */
    J[93] += dqdci;               /* dwdot[HCCO]/d[O] */
    /* d()/d[C2H2] */
    dqdci =  + k_f*sc[2];
    J[727] += dqdci;              /* dwdot[H]/d[C2H2] */
    J[728] -= dqdci;              /* dwdot[O]/d[C2H2] */
    J[748] -= dqdci;              /* dwdot[C2H2]/d[C2H2] */
    J[753] += dqdci;              /* dwdot[HCCO]/d[C2H2] */
    /* d()/d[HCCO] */
    dqdci =  - k_r*sc[1];
    J[892] += dqdci;              /* dwdot[H]/d[HCCO] */
    J[893] -= dqdci;              /* dwdot[O]/d[HCCO] */
    J[913] -= dqdci;              /* dwdot[C2H2]/d[HCCO] */
    J[918] += dqdci;              /* dwdot[HCCO]/d[HCCO] */
    /* d()/dT */
    J[1057] += dqdT;              /* dwdot[H]/dT */
    J[1058] -= dqdT;              /* dwdot[O]/dT */
    J[1078] -= dqdT;              /* dwdot[C2H2]/dT */
    J[1083] += dqdT;              /* dwdot[HCCO]/dT */

    /*reaction 46: O + C2H2 <=> OH + C2H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[22];
    k_f = prefactor_units[45] * fwd_A[45]
                * exp(fwd_beta[45] * tc[0] - activation_units[45] * fwd_Ea[45] * invT);
    dlnkfdT = fwd_beta[45] * invT + activation_units[45] * fwd_Ea[45] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[21];
    Kc = exp(g_RT[2] - g_RT[4] - g_RT[21] + g_RT[22]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[22]) + (h_RT[4] + h_RT[21]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    wdot[21] += q; /* C2H */
    wdot[22] -= q; /* C2H2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[22];
    J[68] -= dqdci;               /* dwdot[O]/d[O] */
    J[70] += dqdci;               /* dwdot[OH]/d[O] */
    J[87] += dqdci;               /* dwdot[C2H]/d[O] */
    J[88] -= dqdci;               /* dwdot[C2H2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[21];
    J[134] -= dqdci;              /* dwdot[O]/d[OH] */
    J[136] += dqdci;              /* dwdot[OH]/d[OH] */
    J[153] += dqdci;              /* dwdot[C2H]/d[OH] */
    J[154] -= dqdci;              /* dwdot[C2H2]/d[OH] */
    /* d()/d[C2H] */
    dqdci =  - k_r*sc[4];
    J[695] -= dqdci;              /* dwdot[O]/d[C2H] */
    J[697] += dqdci;              /* dwdot[OH]/d[C2H] */
    J[714] += dqdci;              /* dwdot[C2H]/d[C2H] */
    J[715] -= dqdci;              /* dwdot[C2H2]/d[C2H] */
    /* d()/d[C2H2] */
    dqdci =  + k_f*sc[2];
    J[728] -= dqdci;              /* dwdot[O]/d[C2H2] */
    J[730] += dqdci;              /* dwdot[OH]/d[C2H2] */
    J[747] += dqdci;              /* dwdot[C2H]/d[C2H2] */
    J[748] -= dqdci;              /* dwdot[C2H2]/d[C2H2] */
    /* d()/dT */
    J[1058] -= dqdT;              /* dwdot[O]/dT */
    J[1060] += dqdT;              /* dwdot[OH]/dT */
    J[1077] += dqdT;              /* dwdot[C2H]/dT */
    J[1078] -= dqdT;              /* dwdot[C2H2]/dT */

    /*reaction 47: O + C2H2 <=> CO + CH2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[22];
    k_f = prefactor_units[46] * fwd_A[46]
                * exp(fwd_beta[46] * tc[0] - activation_units[46] * fwd_Ea[46] * invT);
    dlnkfdT = fwd_beta[46] * invT + activation_units[46] * fwd_Ea[46] * invT2;
    /* reverse */
    phi_r = sc[10]*sc[14];
    Kc = exp(g_RT[2] - g_RT[10] - g_RT[14] + g_RT[22]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[22]) + (h_RT[10] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[10] += q; /* CH2 */
    wdot[14] += q; /* CO */
    wdot[22] -= q; /* C2H2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[22];
    J[68] -= dqdci;               /* dwdot[O]/d[O] */
    J[76] += dqdci;               /* dwdot[CH2]/d[O] */
    J[80] += dqdci;               /* dwdot[CO]/d[O] */
    J[88] -= dqdci;               /* dwdot[C2H2]/d[O] */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[14];
    J[332] -= dqdci;              /* dwdot[O]/d[CH2] */
    J[340] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[344] += dqdci;              /* dwdot[CO]/d[CH2] */
    J[352] -= dqdci;              /* dwdot[C2H2]/d[CH2] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[10];
    J[464] -= dqdci;              /* dwdot[O]/d[CO] */
    J[472] += dqdci;              /* dwdot[CH2]/d[CO] */
    J[476] += dqdci;              /* dwdot[CO]/d[CO] */
    J[484] -= dqdci;              /* dwdot[C2H2]/d[CO] */
    /* d()/d[C2H2] */
    dqdci =  + k_f*sc[2];
    J[728] -= dqdci;              /* dwdot[O]/d[C2H2] */
    J[736] += dqdci;              /* dwdot[CH2]/d[C2H2] */
    J[740] += dqdci;              /* dwdot[CO]/d[C2H2] */
    J[748] -= dqdci;              /* dwdot[C2H2]/d[C2H2] */
    /* d()/dT */
    J[1058] -= dqdT;              /* dwdot[O]/dT */
    J[1066] += dqdT;              /* dwdot[CH2]/dT */
    J[1070] += dqdT;              /* dwdot[CO]/dT */
    J[1078] -= dqdT;              /* dwdot[C2H2]/dT */

    /*reaction 48: O + C2H3 <=> H + CH2CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[23];
    k_f = prefactor_units[47] * fwd_A[47]
                * exp(fwd_beta[47] * tc[0] - activation_units[47] * fwd_Ea[47] * invT);
    dlnkfdT = fwd_beta[47] * invT + activation_units[47] * fwd_Ea[47] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[28];
    Kc = exp(-g_RT[1] + g_RT[2] + g_RT[23] - g_RT[28]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[23]) + (h_RT[1] + h_RT[28]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[2] -= q; /* O */
    wdot[23] -= q; /* C2H3 */
    wdot[28] += q; /* CH2CO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[28];
    J[34] += dqdci;               /* dwdot[H]/d[H] */
    J[35] -= dqdci;               /* dwdot[O]/d[H] */
    J[56] -= dqdci;               /* dwdot[C2H3]/d[H] */
    J[61] += dqdci;               /* dwdot[CH2CO]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[23];
    J[67] += dqdci;               /* dwdot[H]/d[O] */
    J[68] -= dqdci;               /* dwdot[O]/d[O] */
    J[89] -= dqdci;               /* dwdot[C2H3]/d[O] */
    J[94] += dqdci;               /* dwdot[CH2CO]/d[O] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[2];
    J[760] += dqdci;              /* dwdot[H]/d[C2H3] */
    J[761] -= dqdci;              /* dwdot[O]/d[C2H3] */
    J[782] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    J[787] += dqdci;              /* dwdot[CH2CO]/d[C2H3] */
    /* d()/d[CH2CO] */
    dqdci =  - k_r*sc[1];
    J[925] += dqdci;              /* dwdot[H]/d[CH2CO] */
    J[926] -= dqdci;              /* dwdot[O]/d[CH2CO] */
    J[947] -= dqdci;              /* dwdot[C2H3]/d[CH2CO] */
    J[952] += dqdci;              /* dwdot[CH2CO]/d[CH2CO] */
    /* d()/dT */
    J[1057] += dqdT;              /* dwdot[H]/dT */
    J[1058] -= dqdT;              /* dwdot[O]/dT */
    J[1079] -= dqdT;              /* dwdot[C2H3]/dT */
    J[1084] += dqdT;              /* dwdot[CH2CO]/dT */

    /*reaction 49: O + C2H4 <=> CH3 + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[24];
    k_f = prefactor_units[48] * fwd_A[48]
                * exp(fwd_beta[48] * tc[0] - activation_units[48] * fwd_Ea[48] * invT);
    dlnkfdT = fwd_beta[48] * invT + activation_units[48] * fwd_Ea[48] * invT2;
    /* reverse */
    phi_r = sc[12]*sc[16];
    Kc = exp(g_RT[2] - g_RT[12] - g_RT[16] + g_RT[24]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[24]) + (h_RT[12] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[12] += q; /* CH3 */
    wdot[16] += q; /* HCO */
    wdot[24] -= q; /* C2H4 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[24];
    J[68] -= dqdci;               /* dwdot[O]/d[O] */
    J[78] += dqdci;               /* dwdot[CH3]/d[O] */
    J[82] += dqdci;               /* dwdot[HCO]/d[O] */
    J[90] -= dqdci;               /* dwdot[C2H4]/d[O] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[16];
    J[398] -= dqdci;              /* dwdot[O]/d[CH3] */
    J[408] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[412] += dqdci;              /* dwdot[HCO]/d[CH3] */
    J[420] -= dqdci;              /* dwdot[C2H4]/d[CH3] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[12];
    J[530] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[540] += dqdci;              /* dwdot[CH3]/d[HCO] */
    J[544] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[552] -= dqdci;              /* dwdot[C2H4]/d[HCO] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[2];
    J[794] -= dqdci;              /* dwdot[O]/d[C2H4] */
    J[804] += dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[808] += dqdci;              /* dwdot[HCO]/d[C2H4] */
    J[816] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    /* d()/dT */
    J[1058] -= dqdT;              /* dwdot[O]/dT */
    J[1068] += dqdT;              /* dwdot[CH3]/dT */
    J[1072] += dqdT;              /* dwdot[HCO]/dT */
    J[1080] -= dqdT;              /* dwdot[C2H4]/dT */

    /*reaction 50: O + C2H5 <=> CH3 + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[25];
    k_f = prefactor_units[49] * fwd_A[49]
                * exp(fwd_beta[49] * tc[0] - activation_units[49] * fwd_Ea[49] * invT);
    dlnkfdT = fwd_beta[49] * invT + activation_units[49] * fwd_Ea[49] * invT2;
    /* reverse */
    phi_r = sc[12]*sc[17];
    Kc = exp(g_RT[2] - g_RT[12] - g_RT[17] + g_RT[25]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[25]) + (h_RT[12] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[12] += q; /* CH3 */
    wdot[17] += q; /* CH2O */
    wdot[25] -= q; /* C2H5 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[25];
    J[68] -= dqdci;               /* dwdot[O]/d[O] */
    J[78] += dqdci;               /* dwdot[CH3]/d[O] */
    J[83] += dqdci;               /* dwdot[CH2O]/d[O] */
    J[91] -= dqdci;               /* dwdot[C2H5]/d[O] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[17];
    J[398] -= dqdci;              /* dwdot[O]/d[CH3] */
    J[408] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[413] += dqdci;              /* dwdot[CH2O]/d[CH3] */
    J[421] -= dqdci;              /* dwdot[C2H5]/d[CH3] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[12];
    J[563] -= dqdci;              /* dwdot[O]/d[CH2O] */
    J[573] += dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[578] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[586] -= dqdci;              /* dwdot[C2H5]/d[CH2O] */
    /* d()/d[C2H5] */
    dqdci =  + k_f*sc[2];
    J[827] -= dqdci;              /* dwdot[O]/d[C2H5] */
    J[837] += dqdci;              /* dwdot[CH3]/d[C2H5] */
    J[842] += dqdci;              /* dwdot[CH2O]/d[C2H5] */
    J[850] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
    /* d()/dT */
    J[1058] -= dqdT;              /* dwdot[O]/dT */
    J[1068] += dqdT;              /* dwdot[CH3]/dT */
    J[1073] += dqdT;              /* dwdot[CH2O]/dT */
    J[1081] -= dqdT;              /* dwdot[C2H5]/dT */

    /*reaction 51: O + C2H6 <=> OH + C2H5 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[26];
    k_f = prefactor_units[50] * fwd_A[50]
                * exp(fwd_beta[50] * tc[0] - activation_units[50] * fwd_Ea[50] * invT);
    dlnkfdT = fwd_beta[50] * invT + activation_units[50] * fwd_Ea[50] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[25];
    Kc = exp(g_RT[2] - g_RT[4] - g_RT[25] + g_RT[26]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[26]) + (h_RT[4] + h_RT[25]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    wdot[25] += q; /* C2H5 */
    wdot[26] -= q; /* C2H6 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[26];
    J[68] -= dqdci;               /* dwdot[O]/d[O] */
    J[70] += dqdci;               /* dwdot[OH]/d[O] */
    J[91] += dqdci;               /* dwdot[C2H5]/d[O] */
    J[92] -= dqdci;               /* dwdot[C2H6]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[25];
    J[134] -= dqdci;              /* dwdot[O]/d[OH] */
    J[136] += dqdci;              /* dwdot[OH]/d[OH] */
    J[157] += dqdci;              /* dwdot[C2H5]/d[OH] */
    J[158] -= dqdci;              /* dwdot[C2H6]/d[OH] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[4];
    J[827] -= dqdci;              /* dwdot[O]/d[C2H5] */
    J[829] += dqdci;              /* dwdot[OH]/d[C2H5] */
    J[850] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[851] -= dqdci;              /* dwdot[C2H6]/d[C2H5] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[2];
    J[860] -= dqdci;              /* dwdot[O]/d[C2H6] */
    J[862] += dqdci;              /* dwdot[OH]/d[C2H6] */
    J[883] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[884] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[1058] -= dqdT;              /* dwdot[O]/dT */
    J[1060] += dqdT;              /* dwdot[OH]/dT */
    J[1081] += dqdT;              /* dwdot[C2H5]/dT */
    J[1082] -= dqdT;              /* dwdot[C2H6]/dT */

    /*reaction 52: O + HCCO <=> H + 2.000000 CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[27];
    k_f = prefactor_units[51] * fwd_A[51]
                * exp(fwd_beta[51] * tc[0] - activation_units[51] * fwd_Ea[51] * invT);
    dlnkfdT = fwd_beta[51] * invT + activation_units[51] * fwd_Ea[51] * invT2;
    /* reverse */
    phi_r = sc[1]*pow(sc[14], 2.000000);
    Kc = refC * exp(-g_RT[1] + g_RT[2] - 2.000000*g_RT[14] + g_RT[27]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[27]) + (h_RT[1] + 2.000000*h_RT[14]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[2] -= q; /* O */
    wdot[14] += 2 * q; /* CO */
    wdot[27] -= q; /* HCCO */
    /* d()/d[H] */
    dqdci =  - k_r*pow(sc[14], 2.000000);
    J[34] += dqdci;               /* dwdot[H]/d[H] */
    J[35] -= dqdci;               /* dwdot[O]/d[H] */
    J[47] += 2 * dqdci;           /* dwdot[CO]/d[H] */
    J[60] -= dqdci;               /* dwdot[HCCO]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[27];
    J[67] += dqdci;               /* dwdot[H]/d[O] */
    J[68] -= dqdci;               /* dwdot[O]/d[O] */
    J[80] += 2 * dqdci;           /* dwdot[CO]/d[O] */
    J[93] -= dqdci;               /* dwdot[HCCO]/d[O] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[1]*2.000000*sc[14];
    J[463] += dqdci;              /* dwdot[H]/d[CO] */
    J[464] -= dqdci;              /* dwdot[O]/d[CO] */
    J[476] += 2 * dqdci;          /* dwdot[CO]/d[CO] */
    J[489] -= dqdci;              /* dwdot[HCCO]/d[CO] */
    /* d()/d[HCCO] */
    dqdci =  + k_f*sc[2];
    J[892] += dqdci;              /* dwdot[H]/d[HCCO] */
    J[893] -= dqdci;              /* dwdot[O]/d[HCCO] */
    J[905] += 2 * dqdci;          /* dwdot[CO]/d[HCCO] */
    J[918] -= dqdci;              /* dwdot[HCCO]/d[HCCO] */
    /* d()/dT */
    J[1057] += dqdT;              /* dwdot[H]/dT */
    J[1058] -= dqdT;              /* dwdot[O]/dT */
    J[1070] += 2 * dqdT;          /* dwdot[CO]/dT */
    J[1083] -= dqdT;              /* dwdot[HCCO]/dT */

    /*reaction 53: O + CH2CO <=> OH + HCCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[28];
    k_f = prefactor_units[52] * fwd_A[52]
                * exp(fwd_beta[52] * tc[0] - activation_units[52] * fwd_Ea[52] * invT);
    dlnkfdT = fwd_beta[52] * invT + activation_units[52] * fwd_Ea[52] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[27];
    Kc = exp(g_RT[2] - g_RT[4] - g_RT[27] + g_RT[28]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[28]) + (h_RT[4] + h_RT[27]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    wdot[27] += q; /* HCCO */
    wdot[28] -= q; /* CH2CO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[28];
    J[68] -= dqdci;               /* dwdot[O]/d[O] */
    J[70] += dqdci;               /* dwdot[OH]/d[O] */
    J[93] += dqdci;               /* dwdot[HCCO]/d[O] */
    J[94] -= dqdci;               /* dwdot[CH2CO]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[27];
    J[134] -= dqdci;              /* dwdot[O]/d[OH] */
    J[136] += dqdci;              /* dwdot[OH]/d[OH] */
    J[159] += dqdci;              /* dwdot[HCCO]/d[OH] */
    J[160] -= dqdci;              /* dwdot[CH2CO]/d[OH] */
    /* d()/d[HCCO] */
    dqdci =  - k_r*sc[4];
    J[893] -= dqdci;              /* dwdot[O]/d[HCCO] */
    J[895] += dqdci;              /* dwdot[OH]/d[HCCO] */
    J[918] += dqdci;              /* dwdot[HCCO]/d[HCCO] */
    J[919] -= dqdci;              /* dwdot[CH2CO]/d[HCCO] */
    /* d()/d[CH2CO] */
    dqdci =  + k_f*sc[2];
    J[926] -= dqdci;              /* dwdot[O]/d[CH2CO] */
    J[928] += dqdci;              /* dwdot[OH]/d[CH2CO] */
    J[951] += dqdci;              /* dwdot[HCCO]/d[CH2CO] */
    J[952] -= dqdci;              /* dwdot[CH2CO]/d[CH2CO] */
    /* d()/dT */
    J[1058] -= dqdT;              /* dwdot[O]/dT */
    J[1060] += dqdT;              /* dwdot[OH]/dT */
    J[1083] += dqdT;              /* dwdot[HCCO]/dT */
    J[1084] -= dqdT;              /* dwdot[CH2CO]/dT */

    /*reaction 54: O + CH2CO <=> CH2 + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[28];
    k_f = prefactor_units[53] * fwd_A[53]
                * exp(fwd_beta[53] * tc[0] - activation_units[53] * fwd_Ea[53] * invT);
    dlnkfdT = fwd_beta[53] * invT + activation_units[53] * fwd_Ea[53] * invT2;
    /* reverse */
    phi_r = sc[10]*sc[15];
    Kc = exp(g_RT[2] - g_RT[10] - g_RT[15] + g_RT[28]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[28]) + (h_RT[10] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[10] += q; /* CH2 */
    wdot[15] += q; /* CO2 */
    wdot[28] -= q; /* CH2CO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[28];
    J[68] -= dqdci;               /* dwdot[O]/d[O] */
    J[76] += dqdci;               /* dwdot[CH2]/d[O] */
    J[81] += dqdci;               /* dwdot[CO2]/d[O] */
    J[94] -= dqdci;               /* dwdot[CH2CO]/d[O] */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[15];
    J[332] -= dqdci;              /* dwdot[O]/d[CH2] */
    J[340] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[345] += dqdci;              /* dwdot[CO2]/d[CH2] */
    J[358] -= dqdci;              /* dwdot[CH2CO]/d[CH2] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[10];
    J[497] -= dqdci;              /* dwdot[O]/d[CO2] */
    J[505] += dqdci;              /* dwdot[CH2]/d[CO2] */
    J[510] += dqdci;              /* dwdot[CO2]/d[CO2] */
    J[523] -= dqdci;              /* dwdot[CH2CO]/d[CO2] */
    /* d()/d[CH2CO] */
    dqdci =  + k_f*sc[2];
    J[926] -= dqdci;              /* dwdot[O]/d[CH2CO] */
    J[934] += dqdci;              /* dwdot[CH2]/d[CH2CO] */
    J[939] += dqdci;              /* dwdot[CO2]/d[CH2CO] */
    J[952] -= dqdci;              /* dwdot[CH2CO]/d[CH2CO] */
    /* d()/dT */
    J[1058] -= dqdT;              /* dwdot[O]/dT */
    J[1066] += dqdT;              /* dwdot[CH2]/dT */
    J[1071] += dqdT;              /* dwdot[CO2]/dT */
    J[1084] -= dqdT;              /* dwdot[CH2CO]/dT */

    /*reaction 55: O2 + CO <=> O + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[14];
    k_f = prefactor_units[54] * fwd_A[54]
                * exp(fwd_beta[54] * tc[0] - activation_units[54] * fwd_Ea[54] * invT);
    dlnkfdT = fwd_beta[54] * invT + activation_units[54] * fwd_Ea[54] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[15];
    Kc = exp(-g_RT[2] + g_RT[3] + g_RT[14] - g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[14]) + (h_RT[2] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* O */
    wdot[3] -= q; /* O2 */
    wdot[14] -= q; /* CO */
    wdot[15] += q; /* CO2 */
    /* d()/d[O] */
    dqdci =  - k_r*sc[15];
    J[68] += dqdci;               /* dwdot[O]/d[O] */
    J[69] -= dqdci;               /* dwdot[O2]/d[O] */
    J[80] -= dqdci;               /* dwdot[CO]/d[O] */
    J[81] += dqdci;               /* dwdot[CO2]/d[O] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[14];
    J[101] += dqdci;              /* dwdot[O]/d[O2] */
    J[102] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[113] -= dqdci;              /* dwdot[CO]/d[O2] */
    J[114] += dqdci;              /* dwdot[CO2]/d[O2] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[3];
    J[464] += dqdci;              /* dwdot[O]/d[CO] */
    J[465] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[476] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[477] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[2];
    J[497] += dqdci;              /* dwdot[O]/d[CO2] */
    J[498] -= dqdci;              /* dwdot[O2]/d[CO2] */
    J[509] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[510] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[1058] += dqdT;              /* dwdot[O]/dT */
    J[1059] -= dqdT;              /* dwdot[O2]/dT */
    J[1070] -= dqdT;              /* dwdot[CO]/dT */
    J[1071] += dqdT;              /* dwdot[CO2]/dT */

    /*reaction 56: O2 + CH2O <=> HO2 + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[17];
    k_f = prefactor_units[55] * fwd_A[55]
                * exp(fwd_beta[55] * tc[0] - activation_units[55] * fwd_Ea[55] * invT);
    dlnkfdT = fwd_beta[55] * invT + activation_units[55] * fwd_Ea[55] * invT2;
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
    wdot[16] += q; /* HCO */
    wdot[17] -= q; /* CH2O */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[17];
    J[102] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[105] += dqdci;              /* dwdot[HO2]/d[O2] */
    J[115] += dqdci;              /* dwdot[HCO]/d[O2] */
    J[116] -= dqdci;              /* dwdot[CH2O]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[16];
    J[201] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[204] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[214] += dqdci;              /* dwdot[HCO]/d[HO2] */
    J[215] -= dqdci;              /* dwdot[CH2O]/d[HO2] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[6];
    J[531] -= dqdci;              /* dwdot[O2]/d[HCO] */
    J[534] += dqdci;              /* dwdot[HO2]/d[HCO] */
    J[544] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[545] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[3];
    J[564] -= dqdci;              /* dwdot[O2]/d[CH2O] */
    J[567] += dqdci;              /* dwdot[HO2]/d[CH2O] */
    J[577] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[578] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[1059] -= dqdT;              /* dwdot[O2]/dT */
    J[1062] += dqdT;              /* dwdot[HO2]/dT */
    J[1072] += dqdT;              /* dwdot[HCO]/dT */
    J[1073] -= dqdT;              /* dwdot[CH2O]/dT */

    /*reaction 57: H + 2.000000 O2 <=> HO2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*pow(sc[3], 2.000000);
    k_f = prefactor_units[56] * fwd_A[56]
                * exp(fwd_beta[56] * tc[0] - activation_units[56] * fwd_Ea[56] * invT);
    dlnkfdT = fwd_beta[56] * invT + activation_units[56] * fwd_Ea[56] * invT2;
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
    J[34] -= dqdci;               /* dwdot[H]/d[H] */
    J[36] -= dqdci;               /* dwdot[O2]/d[H] */
    J[39] += dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1]*2.000000*sc[3] - k_r*sc[6];
    J[100] -= dqdci;              /* dwdot[H]/d[O2] */
    J[102] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[105] += dqdci;              /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[3];
    J[199] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[201] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[204] += dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1057] -= dqdT;              /* dwdot[H]/dT */
    J[1059] -= dqdT;              /* dwdot[O2]/dT */
    J[1062] += dqdT;              /* dwdot[HO2]/dT */

    /*reaction 58: H + O2 + H2O <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3]*sc[5];
    k_f = prefactor_units[57] * fwd_A[57]
                * exp(fwd_beta[57] * tc[0] - activation_units[57] * fwd_Ea[57] * invT);
    dlnkfdT = fwd_beta[57] * invT + activation_units[57] * fwd_Ea[57] * invT2;
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
    J[34] -= dqdci;               /* dwdot[H]/d[H] */
    J[36] -= dqdci;               /* dwdot[O2]/d[H] */
    J[39] += dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1]*sc[5];
    J[100] -= dqdci;              /* dwdot[H]/d[O2] */
    J[102] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[105] += dqdci;              /* dwdot[HO2]/d[O2] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[1]*sc[3] - k_r*sc[6];
    J[166] -= dqdci;              /* dwdot[H]/d[H2O] */
    J[168] -= dqdci;              /* dwdot[O2]/d[H2O] */
    J[171] += dqdci;              /* dwdot[HO2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[5];
    J[199] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[201] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[204] += dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1057] -= dqdT;              /* dwdot[H]/dT */
    J[1059] -= dqdT;              /* dwdot[O2]/dT */
    J[1062] += dqdT;              /* dwdot[HO2]/dT */

    /*reaction 59: H + O2 + N2 <=> HO2 + N2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3]*sc[30];
    k_f = prefactor_units[58] * fwd_A[58]
                * exp(fwd_beta[58] * tc[0] - activation_units[58] * fwd_Ea[58] * invT);
    dlnkfdT = fwd_beta[58] * invT + activation_units[58] * fwd_Ea[58] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[30];
    Kc = refCinv * exp(g_RT[1] + g_RT[3] - g_RT[6] + g_RT[30] - g_RT[30]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3] + h_RT[30]) + (h_RT[6] + h_RT[30]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[3]*sc[30];
    J[34] -= dqdci;               /* dwdot[H]/d[H] */
    J[36] -= dqdci;               /* dwdot[O2]/d[H] */
    J[39] += dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1]*sc[30];
    J[100] -= dqdci;              /* dwdot[H]/d[O2] */
    J[102] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[105] += dqdci;              /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[30];
    J[199] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[201] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[204] += dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/d[N2] */
    dqdci =  + k_f*sc[1]*sc[3] - k_r*sc[6];
    J[991] -= dqdci;              /* dwdot[H]/d[N2] */
    J[993] -= dqdci;              /* dwdot[O2]/d[N2] */
    J[996] += dqdci;              /* dwdot[HO2]/d[N2] */
    /* d()/dT */
    J[1057] -= dqdT;              /* dwdot[H]/dT */
    J[1059] -= dqdT;              /* dwdot[O2]/dT */
    J[1062] += dqdT;              /* dwdot[HO2]/dT */

    /*reaction 60: H + O2 + AR <=> HO2 + AR */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3]*sc[31];
    k_f = prefactor_units[59] * fwd_A[59]
                * exp(fwd_beta[59] * tc[0] - activation_units[59] * fwd_Ea[59] * invT);
    dlnkfdT = fwd_beta[59] * invT + activation_units[59] * fwd_Ea[59] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[31];
    Kc = refCinv * exp(g_RT[1] + g_RT[3] - g_RT[6] + g_RT[31] - g_RT[31]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3] + h_RT[31]) + (h_RT[6] + h_RT[31]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[3]*sc[31];
    J[34] -= dqdci;               /* dwdot[H]/d[H] */
    J[36] -= dqdci;               /* dwdot[O2]/d[H] */
    J[39] += dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1]*sc[31];
    J[100] -= dqdci;              /* dwdot[H]/d[O2] */
    J[102] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[105] += dqdci;              /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[31];
    J[199] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[201] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[204] += dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/d[AR] */
    dqdci =  + k_f*sc[1]*sc[3] - k_r*sc[6];
    J[1024] -= dqdci;             /* dwdot[H]/d[AR] */
    J[1026] -= dqdci;             /* dwdot[O2]/d[AR] */
    J[1029] += dqdci;             /* dwdot[HO2]/d[AR] */
    /* d()/dT */
    J[1057] -= dqdT;              /* dwdot[H]/dT */
    J[1059] -= dqdT;              /* dwdot[O2]/dT */
    J[1062] += dqdT;              /* dwdot[HO2]/dT */

    /*reaction 61: H + O2 <=> O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = prefactor_units[60] * fwd_A[60]
                * exp(fwd_beta[60] * tc[0] - activation_units[60] * fwd_Ea[60] * invT);
    dlnkfdT = fwd_beta[60] * invT + activation_units[60] * fwd_Ea[60] * invT2;
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
    J[34] -= dqdci;               /* dwdot[H]/d[H] */
    J[35] += dqdci;               /* dwdot[O]/d[H] */
    J[36] -= dqdci;               /* dwdot[O2]/d[H] */
    J[37] += dqdci;               /* dwdot[OH]/d[H] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[4];
    J[67] -= dqdci;               /* dwdot[H]/d[O] */
    J[68] += dqdci;               /* dwdot[O]/d[O] */
    J[69] -= dqdci;               /* dwdot[O2]/d[O] */
    J[70] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1];
    J[100] -= dqdci;              /* dwdot[H]/d[O2] */
    J[101] += dqdci;              /* dwdot[O]/d[O2] */
    J[102] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[103] += dqdci;              /* dwdot[OH]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[2];
    J[133] -= dqdci;              /* dwdot[H]/d[OH] */
    J[134] += dqdci;              /* dwdot[O]/d[OH] */
    J[135] -= dqdci;              /* dwdot[O2]/d[OH] */
    J[136] += dqdci;              /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[1057] -= dqdT;              /* dwdot[H]/dT */
    J[1058] += dqdT;              /* dwdot[O]/dT */
    J[1059] -= dqdT;              /* dwdot[O2]/dT */
    J[1060] += dqdT;              /* dwdot[OH]/dT */

    /*reaction 62: 2.000000 H + H2 <=> 2.000000 H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*pow(sc[1], 2.000000);
    k_f = prefactor_units[61] * fwd_A[61]
                * exp(fwd_beta[61] * tc[0] - activation_units[61] * fwd_Ea[61] * invT);
    dlnkfdT = fwd_beta[61] * invT + activation_units[61] * fwd_Ea[61] * invT2;
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
    J[33] += dqdci;               /* dwdot[H2]/d[H] */
    J[34] += -2 * dqdci;          /* dwdot[H]/d[H] */
    /* d()/dT */
    J[1056] += dqdT;              /* dwdot[H2]/dT */
    J[1057] += -2 * dqdT;         /* dwdot[H]/dT */

    /*reaction 63: 2.000000 H + H2O <=> H2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[1], 2.000000)*sc[5];
    k_f = prefactor_units[62] * fwd_A[62]
                * exp(fwd_beta[62] * tc[0] - activation_units[62] * fwd_Ea[62] * invT);
    dlnkfdT = fwd_beta[62] * invT + activation_units[62] * fwd_Ea[62] * invT2;
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
    J[33] += dqdci;               /* dwdot[H2]/d[H] */
    J[34] += -2 * dqdci;          /* dwdot[H]/d[H] */
    /* d()/d[H2O] */
    dqdci =  + k_f*pow(sc[1], 2.000000) - k_r*sc[0];
    J[165] += dqdci;              /* dwdot[H2]/d[H2O] */
    J[166] += -2 * dqdci;         /* dwdot[H]/d[H2O] */
    /* d()/dT */
    J[1056] += dqdT;              /* dwdot[H2]/dT */
    J[1057] += -2 * dqdT;         /* dwdot[H]/dT */

    /*reaction 64: 2.000000 H + CO2 <=> H2 + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[1], 2.000000)*sc[15];
    k_f = prefactor_units[63] * fwd_A[63]
                * exp(fwd_beta[63] * tc[0] - activation_units[63] * fwd_Ea[63] * invT);
    dlnkfdT = fwd_beta[63] * invT + activation_units[63] * fwd_Ea[63] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[15];
    Kc = refCinv * exp(-g_RT[0] + 2.000000*g_RT[1] + g_RT[15] - g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[1] + h_RT[15]) + (h_RT[0] + h_RT[15]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= 2 * q; /* H */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[15];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += -2 * dqdci;           /* dwdot[H]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*2.000000*sc[1]*sc[15];
    J[33] += dqdci;               /* dwdot[H2]/d[H] */
    J[34] += -2 * dqdci;          /* dwdot[H]/d[H] */
    /* d()/d[CO2] */
    dqdci =  + k_f*pow(sc[1], 2.000000) - k_r*sc[0];
    J[495] += dqdci;              /* dwdot[H2]/d[CO2] */
    J[496] += -2 * dqdci;         /* dwdot[H]/d[CO2] */
    /* d()/dT */
    J[1056] += dqdT;              /* dwdot[H2]/dT */
    J[1057] += -2 * dqdT;         /* dwdot[H]/dT */

    /*reaction 65: H + HO2 <=> O + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[6];
    k_f = prefactor_units[64] * fwd_A[64]
                * exp(fwd_beta[64] * tc[0] - activation_units[64] * fwd_Ea[64] * invT);
    dlnkfdT = fwd_beta[64] * invT + activation_units[64] * fwd_Ea[64] * invT2;
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
    J[34] -= dqdci;               /* dwdot[H]/d[H] */
    J[35] += dqdci;               /* dwdot[O]/d[H] */
    J[38] += dqdci;               /* dwdot[H2O]/d[H] */
    J[39] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[5];
    J[67] -= dqdci;               /* dwdot[H]/d[O] */
    J[68] += dqdci;               /* dwdot[O]/d[O] */
    J[71] += dqdci;               /* dwdot[H2O]/d[O] */
    J[72] -= dqdci;               /* dwdot[HO2]/d[O] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[2];
    J[166] -= dqdci;              /* dwdot[H]/d[H2O] */
    J[167] += dqdci;              /* dwdot[O]/d[H2O] */
    J[170] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[171] -= dqdci;              /* dwdot[HO2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[1];
    J[199] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[200] += dqdci;              /* dwdot[O]/d[HO2] */
    J[203] += dqdci;              /* dwdot[H2O]/d[HO2] */
    J[204] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1057] -= dqdT;              /* dwdot[H]/dT */
    J[1058] += dqdT;              /* dwdot[O]/dT */
    J[1061] += dqdT;              /* dwdot[H2O]/dT */
    J[1062] -= dqdT;              /* dwdot[HO2]/dT */

    /*reaction 66: H + HO2 <=> O2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[6];
    k_f = prefactor_units[65] * fwd_A[65]
                * exp(fwd_beta[65] * tc[0] - activation_units[65] * fwd_Ea[65] * invT);
    dlnkfdT = fwd_beta[65] * invT + activation_units[65] * fwd_Ea[65] * invT2;
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
    J[33] += dqdci;               /* dwdot[H2]/d[H] */
    J[34] -= dqdci;               /* dwdot[H]/d[H] */
    J[36] += dqdci;               /* dwdot[O2]/d[H] */
    J[39] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[0];
    J[99] += dqdci;               /* dwdot[H2]/d[O2] */
    J[100] -= dqdci;              /* dwdot[H]/d[O2] */
    J[102] += dqdci;              /* dwdot[O2]/d[O2] */
    J[105] -= dqdci;              /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[1];
    J[198] += dqdci;              /* dwdot[H2]/d[HO2] */
    J[199] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[201] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[204] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1056] += dqdT;              /* dwdot[H2]/dT */
    J[1057] -= dqdT;              /* dwdot[H]/dT */
    J[1059] += dqdT;              /* dwdot[O2]/dT */
    J[1062] -= dqdT;              /* dwdot[HO2]/dT */

    /*reaction 67: H + HO2 <=> 2.000000 OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[6];
    k_f = prefactor_units[66] * fwd_A[66]
                * exp(fwd_beta[66] * tc[0] - activation_units[66] * fwd_Ea[66] * invT);
    dlnkfdT = fwd_beta[66] * invT + activation_units[66] * fwd_Ea[66] * invT2;
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
    J[34] -= dqdci;               /* dwdot[H]/d[H] */
    J[37] += 2 * dqdci;           /* dwdot[OH]/d[H] */
    J[39] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*2.000000*sc[4];
    J[133] -= dqdci;              /* dwdot[H]/d[OH] */
    J[136] += 2 * dqdci;          /* dwdot[OH]/d[OH] */
    J[138] -= dqdci;              /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[1];
    J[199] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[202] += 2 * dqdci;          /* dwdot[OH]/d[HO2] */
    J[204] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1057] -= dqdT;              /* dwdot[H]/dT */
    J[1060] += 2 * dqdT;          /* dwdot[OH]/dT */
    J[1062] -= dqdT;              /* dwdot[HO2]/dT */

    /*reaction 68: H + H2O2 <=> HO2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[7];
    k_f = prefactor_units[67] * fwd_A[67]
                * exp(fwd_beta[67] * tc[0] - activation_units[67] * fwd_Ea[67] * invT);
    dlnkfdT = fwd_beta[67] * invT + activation_units[67] * fwd_Ea[67] * invT2;
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
    J[33] += dqdci;               /* dwdot[H2]/d[H] */
    J[34] -= dqdci;               /* dwdot[H]/d[H] */
    J[39] += dqdci;               /* dwdot[HO2]/d[H] */
    J[40] -= dqdci;               /* dwdot[H2O2]/d[H] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[0];
    J[198] += dqdci;              /* dwdot[H2]/d[HO2] */
    J[199] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[204] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[205] -= dqdci;              /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[1];
    J[231] += dqdci;              /* dwdot[H2]/d[H2O2] */
    J[232] -= dqdci;              /* dwdot[H]/d[H2O2] */
    J[237] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[238] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[1056] += dqdT;              /* dwdot[H2]/dT */
    J[1057] -= dqdT;              /* dwdot[H]/dT */
    J[1062] += dqdT;              /* dwdot[HO2]/dT */
    J[1063] -= dqdT;              /* dwdot[H2O2]/dT */

    /*reaction 69: H + H2O2 <=> OH + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[7];
    k_f = prefactor_units[68] * fwd_A[68]
                * exp(fwd_beta[68] * tc[0] - activation_units[68] * fwd_Ea[68] * invT);
    dlnkfdT = fwd_beta[68] * invT + activation_units[68] * fwd_Ea[68] * invT2;
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
    J[34] -= dqdci;               /* dwdot[H]/d[H] */
    J[37] += dqdci;               /* dwdot[OH]/d[H] */
    J[38] += dqdci;               /* dwdot[H2O]/d[H] */
    J[40] -= dqdci;               /* dwdot[H2O2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[5];
    J[133] -= dqdci;              /* dwdot[H]/d[OH] */
    J[136] += dqdci;              /* dwdot[OH]/d[OH] */
    J[137] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[139] -= dqdci;              /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[4];
    J[166] -= dqdci;              /* dwdot[H]/d[H2O] */
    J[169] += dqdci;              /* dwdot[OH]/d[H2O] */
    J[170] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[172] -= dqdci;              /* dwdot[H2O2]/d[H2O] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[1];
    J[232] -= dqdci;              /* dwdot[H]/d[H2O2] */
    J[235] += dqdci;              /* dwdot[OH]/d[H2O2] */
    J[236] += dqdci;              /* dwdot[H2O]/d[H2O2] */
    J[238] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[1057] -= dqdT;              /* dwdot[H]/dT */
    J[1060] += dqdT;              /* dwdot[OH]/dT */
    J[1061] += dqdT;              /* dwdot[H2O]/dT */
    J[1063] -= dqdT;              /* dwdot[H2O2]/dT */

    /*reaction 70: H + CH <=> C + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[9];
    k_f = prefactor_units[69] * fwd_A[69]
                * exp(fwd_beta[69] * tc[0] - activation_units[69] * fwd_Ea[69] * invT);
    dlnkfdT = fwd_beta[69] * invT + activation_units[69] * fwd_Ea[69] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[8];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[8] + g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[9]) + (h_RT[0] + h_RT[8]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[8] += q; /* C */
    wdot[9] -= q; /* CH */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[8];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[8] += dqdci;                /* dwdot[C]/d[H2] */
    J[9] -= dqdci;                /* dwdot[CH]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[9];
    J[33] += dqdci;               /* dwdot[H2]/d[H] */
    J[34] -= dqdci;               /* dwdot[H]/d[H] */
    J[41] += dqdci;               /* dwdot[C]/d[H] */
    J[42] -= dqdci;               /* dwdot[CH]/d[H] */
    /* d()/d[C] */
    dqdci =  - k_r*sc[0];
    J[264] += dqdci;              /* dwdot[H2]/d[C] */
    J[265] -= dqdci;              /* dwdot[H]/d[C] */
    J[272] += dqdci;              /* dwdot[C]/d[C] */
    J[273] -= dqdci;              /* dwdot[CH]/d[C] */
    /* d()/d[CH] */
    dqdci =  + k_f*sc[1];
    J[297] += dqdci;              /* dwdot[H2]/d[CH] */
    J[298] -= dqdci;              /* dwdot[H]/d[CH] */
    J[305] += dqdci;              /* dwdot[C]/d[CH] */
    J[306] -= dqdci;              /* dwdot[CH]/d[CH] */
    /* d()/dT */
    J[1056] += dqdT;              /* dwdot[H2]/dT */
    J[1057] -= dqdT;              /* dwdot[H]/dT */
    J[1064] += dqdT;              /* dwdot[C]/dT */
    J[1065] -= dqdT;              /* dwdot[CH]/dT */

    /*reaction 71: H + CH2(S) <=> CH + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[11];
    k_f = prefactor_units[70] * fwd_A[70]
                * exp(fwd_beta[70] * tc[0] - activation_units[70] * fwd_Ea[70] * invT);
    dlnkfdT = fwd_beta[70] * invT + activation_units[70] * fwd_Ea[70] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[9];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[9] + g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[11]) + (h_RT[0] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[9] += q; /* CH */
    wdot[11] -= q; /* CH2(S) */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[9];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[9] += dqdci;                /* dwdot[CH]/d[H2] */
    J[11] -= dqdci;               /* dwdot[CH2(S)]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[11];
    J[33] += dqdci;               /* dwdot[H2]/d[H] */
    J[34] -= dqdci;               /* dwdot[H]/d[H] */
    J[42] += dqdci;               /* dwdot[CH]/d[H] */
    J[44] -= dqdci;               /* dwdot[CH2(S)]/d[H] */
    /* d()/d[CH] */
    dqdci =  - k_r*sc[0];
    J[297] += dqdci;              /* dwdot[H2]/d[CH] */
    J[298] -= dqdci;              /* dwdot[H]/d[CH] */
    J[306] += dqdci;              /* dwdot[CH]/d[CH] */
    J[308] -= dqdci;              /* dwdot[CH2(S)]/d[CH] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[1];
    J[363] += dqdci;              /* dwdot[H2]/d[CH2(S)] */
    J[364] -= dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[372] += dqdci;              /* dwdot[CH]/d[CH2(S)] */
    J[374] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    /* d()/dT */
    J[1056] += dqdT;              /* dwdot[H2]/dT */
    J[1057] -= dqdT;              /* dwdot[H]/dT */
    J[1065] += dqdT;              /* dwdot[CH]/dT */
    J[1067] -= dqdT;              /* dwdot[CH2(S)]/dT */

    /*reaction 72: H + CH4 <=> CH3 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[13];
    k_f = prefactor_units[71] * fwd_A[71]
                * exp(fwd_beta[71] * tc[0] - activation_units[71] * fwd_Ea[71] * invT);
    dlnkfdT = fwd_beta[71] * invT + activation_units[71] * fwd_Ea[71] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[12];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[12] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[13]) + (h_RT[0] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[12] += q; /* CH3 */
    wdot[13] -= q; /* CH4 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[12];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[12] += dqdci;               /* dwdot[CH3]/d[H2] */
    J[13] -= dqdci;               /* dwdot[CH4]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[13];
    J[33] += dqdci;               /* dwdot[H2]/d[H] */
    J[34] -= dqdci;               /* dwdot[H]/d[H] */
    J[45] += dqdci;               /* dwdot[CH3]/d[H] */
    J[46] -= dqdci;               /* dwdot[CH4]/d[H] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[0];
    J[396] += dqdci;              /* dwdot[H2]/d[CH3] */
    J[397] -= dqdci;              /* dwdot[H]/d[CH3] */
    J[408] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[409] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[1];
    J[429] += dqdci;              /* dwdot[H2]/d[CH4] */
    J[430] -= dqdci;              /* dwdot[H]/d[CH4] */
    J[441] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[442] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[1056] += dqdT;              /* dwdot[H2]/dT */
    J[1057] -= dqdT;              /* dwdot[H]/dT */
    J[1068] += dqdT;              /* dwdot[CH3]/dT */
    J[1069] -= dqdT;              /* dwdot[CH4]/dT */

    /*reaction 73: H + HCO <=> H2 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[16];
    k_f = prefactor_units[72] * fwd_A[72]
                * exp(fwd_beta[72] * tc[0] - activation_units[72] * fwd_Ea[72] * invT);
    dlnkfdT = fwd_beta[72] * invT + activation_units[72] * fwd_Ea[72] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[14];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[14] + g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[16]) + (h_RT[0] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[14] += q; /* CO */
    wdot[16] -= q; /* HCO */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[14];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[14] += dqdci;               /* dwdot[CO]/d[H2] */
    J[16] -= dqdci;               /* dwdot[HCO]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[16];
    J[33] += dqdci;               /* dwdot[H2]/d[H] */
    J[34] -= dqdci;               /* dwdot[H]/d[H] */
    J[47] += dqdci;               /* dwdot[CO]/d[H] */
    J[49] -= dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[0];
    J[462] += dqdci;              /* dwdot[H2]/d[CO] */
    J[463] -= dqdci;              /* dwdot[H]/d[CO] */
    J[476] += dqdci;              /* dwdot[CO]/d[CO] */
    J[478] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[1];
    J[528] += dqdci;              /* dwdot[H2]/d[HCO] */
    J[529] -= dqdci;              /* dwdot[H]/d[HCO] */
    J[542] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[544] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[1056] += dqdT;              /* dwdot[H2]/dT */
    J[1057] -= dqdT;              /* dwdot[H]/dT */
    J[1070] += dqdT;              /* dwdot[CO]/dT */
    J[1072] -= dqdT;              /* dwdot[HCO]/dT */

    /*reaction 74: H + CH2O <=> HCO + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[17];
    k_f = prefactor_units[73] * fwd_A[73]
                * exp(fwd_beta[73] * tc[0] - activation_units[73] * fwd_Ea[73] * invT);
    dlnkfdT = fwd_beta[73] * invT + activation_units[73] * fwd_Ea[73] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[16];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[16] + g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[17]) + (h_RT[0] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[16] += q; /* HCO */
    wdot[17] -= q; /* CH2O */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[16];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[16] += dqdci;               /* dwdot[HCO]/d[H2] */
    J[17] -= dqdci;               /* dwdot[CH2O]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[17];
    J[33] += dqdci;               /* dwdot[H2]/d[H] */
    J[34] -= dqdci;               /* dwdot[H]/d[H] */
    J[49] += dqdci;               /* dwdot[HCO]/d[H] */
    J[50] -= dqdci;               /* dwdot[CH2O]/d[H] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[0];
    J[528] += dqdci;              /* dwdot[H2]/d[HCO] */
    J[529] -= dqdci;              /* dwdot[H]/d[HCO] */
    J[544] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[545] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[1];
    J[561] += dqdci;              /* dwdot[H2]/d[CH2O] */
    J[562] -= dqdci;              /* dwdot[H]/d[CH2O] */
    J[577] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[578] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[1056] += dqdT;              /* dwdot[H2]/dT */
    J[1057] -= dqdT;              /* dwdot[H]/dT */
    J[1072] += dqdT;              /* dwdot[HCO]/dT */
    J[1073] -= dqdT;              /* dwdot[CH2O]/dT */

    /*reaction 75: H + CH2OH <=> H2 + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[18];
    k_f = prefactor_units[74] * fwd_A[74]
                * exp(fwd_beta[74] * tc[0] - activation_units[74] * fwd_Ea[74] * invT);
    dlnkfdT = fwd_beta[74] * invT + activation_units[74] * fwd_Ea[74] * invT2;
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
    wdot[17] += q; /* CH2O */
    wdot[18] -= q; /* CH2OH */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[17];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[17] += dqdci;               /* dwdot[CH2O]/d[H2] */
    J[18] -= dqdci;               /* dwdot[CH2OH]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[18];
    J[33] += dqdci;               /* dwdot[H2]/d[H] */
    J[34] -= dqdci;               /* dwdot[H]/d[H] */
    J[50] += dqdci;               /* dwdot[CH2O]/d[H] */
    J[51] -= dqdci;               /* dwdot[CH2OH]/d[H] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[0];
    J[561] += dqdci;              /* dwdot[H2]/d[CH2O] */
    J[562] -= dqdci;              /* dwdot[H]/d[CH2O] */
    J[578] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[579] -= dqdci;              /* dwdot[CH2OH]/d[CH2O] */
    /* d()/d[CH2OH] */
    dqdci =  + k_f*sc[1];
    J[594] += dqdci;              /* dwdot[H2]/d[CH2OH] */
    J[595] -= dqdci;              /* dwdot[H]/d[CH2OH] */
    J[611] += dqdci;              /* dwdot[CH2O]/d[CH2OH] */
    J[612] -= dqdci;              /* dwdot[CH2OH]/d[CH2OH] */
    /* d()/dT */
    J[1056] += dqdT;              /* dwdot[H2]/dT */
    J[1057] -= dqdT;              /* dwdot[H]/dT */
    J[1073] += dqdT;              /* dwdot[CH2O]/dT */
    J[1074] -= dqdT;              /* dwdot[CH2OH]/dT */

    /*reaction 76: H + CH2OH <=> OH + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[18];
    k_f = prefactor_units[75] * fwd_A[75]
                * exp(fwd_beta[75] * tc[0] - activation_units[75] * fwd_Ea[75] * invT);
    dlnkfdT = fwd_beta[75] * invT + activation_units[75] * fwd_Ea[75] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[12];
    Kc = exp(g_RT[1] - g_RT[4] - g_RT[12] + g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[18]) + (h_RT[4] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[4] += q; /* OH */
    wdot[12] += q; /* CH3 */
    wdot[18] -= q; /* CH2OH */
    /* d()/d[H] */
    dqdci =  + k_f*sc[18];
    J[34] -= dqdci;               /* dwdot[H]/d[H] */
    J[37] += dqdci;               /* dwdot[OH]/d[H] */
    J[45] += dqdci;               /* dwdot[CH3]/d[H] */
    J[51] -= dqdci;               /* dwdot[CH2OH]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[12];
    J[133] -= dqdci;              /* dwdot[H]/d[OH] */
    J[136] += dqdci;              /* dwdot[OH]/d[OH] */
    J[144] += dqdci;              /* dwdot[CH3]/d[OH] */
    J[150] -= dqdci;              /* dwdot[CH2OH]/d[OH] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[4];
    J[397] -= dqdci;              /* dwdot[H]/d[CH3] */
    J[400] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[408] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[414] -= dqdci;              /* dwdot[CH2OH]/d[CH3] */
    /* d()/d[CH2OH] */
    dqdci =  + k_f*sc[1];
    J[595] -= dqdci;              /* dwdot[H]/d[CH2OH] */
    J[598] += dqdci;              /* dwdot[OH]/d[CH2OH] */
    J[606] += dqdci;              /* dwdot[CH3]/d[CH2OH] */
    J[612] -= dqdci;              /* dwdot[CH2OH]/d[CH2OH] */
    /* d()/dT */
    J[1057] -= dqdT;              /* dwdot[H]/dT */
    J[1060] += dqdT;              /* dwdot[OH]/dT */
    J[1068] += dqdT;              /* dwdot[CH3]/dT */
    J[1074] -= dqdT;              /* dwdot[CH2OH]/dT */

    /*reaction 77: H + CH2OH <=> CH2(S) + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[18];
    k_f = prefactor_units[76] * fwd_A[76]
                * exp(fwd_beta[76] * tc[0] - activation_units[76] * fwd_Ea[76] * invT);
    dlnkfdT = fwd_beta[76] * invT + activation_units[76] * fwd_Ea[76] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[11];
    Kc = exp(g_RT[1] - g_RT[5] - g_RT[11] + g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[18]) + (h_RT[5] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[5] += q; /* H2O */
    wdot[11] += q; /* CH2(S) */
    wdot[18] -= q; /* CH2OH */
    /* d()/d[H] */
    dqdci =  + k_f*sc[18];
    J[34] -= dqdci;               /* dwdot[H]/d[H] */
    J[38] += dqdci;               /* dwdot[H2O]/d[H] */
    J[44] += dqdci;               /* dwdot[CH2(S)]/d[H] */
    J[51] -= dqdci;               /* dwdot[CH2OH]/d[H] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[11];
    J[166] -= dqdci;              /* dwdot[H]/d[H2O] */
    J[170] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[176] += dqdci;              /* dwdot[CH2(S)]/d[H2O] */
    J[183] -= dqdci;              /* dwdot[CH2OH]/d[H2O] */
    /* d()/d[CH2(S)] */
    dqdci =  - k_r*sc[5];
    J[364] -= dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[368] += dqdci;              /* dwdot[H2O]/d[CH2(S)] */
    J[374] += dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[381] -= dqdci;              /* dwdot[CH2OH]/d[CH2(S)] */
    /* d()/d[CH2OH] */
    dqdci =  + k_f*sc[1];
    J[595] -= dqdci;              /* dwdot[H]/d[CH2OH] */
    J[599] += dqdci;              /* dwdot[H2O]/d[CH2OH] */
    J[605] += dqdci;              /* dwdot[CH2(S)]/d[CH2OH] */
    J[612] -= dqdci;              /* dwdot[CH2OH]/d[CH2OH] */
    /* d()/dT */
    J[1057] -= dqdT;              /* dwdot[H]/dT */
    J[1061] += dqdT;              /* dwdot[H2O]/dT */
    J[1067] += dqdT;              /* dwdot[CH2(S)]/dT */
    J[1074] -= dqdT;              /* dwdot[CH2OH]/dT */

    /*reaction 78: H + CH3O <=> H + CH2OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[19];
    k_f = prefactor_units[77] * fwd_A[77]
                * exp(fwd_beta[77] * tc[0] - activation_units[77] * fwd_Ea[77] * invT);
    dlnkfdT = fwd_beta[77] * invT + activation_units[77] * fwd_Ea[77] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[18];
    Kc = exp(g_RT[1] - g_RT[1] - g_RT[18] + g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[19]) + (h_RT[1] + h_RT[18]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[18] += q; /* CH2OH */
    wdot[19] -= q; /* CH3O */
    /* d()/d[H] */
    dqdci =  + k_f*sc[19] - k_r*sc[18];
    J[51] += dqdci;               /* dwdot[CH2OH]/d[H] */
    J[52] -= dqdci;               /* dwdot[CH3O]/d[H] */
    /* d()/d[CH2OH] */
    dqdci =  - k_r*sc[1];
    J[612] += dqdci;              /* dwdot[CH2OH]/d[CH2OH] */
    J[613] -= dqdci;              /* dwdot[CH3O]/d[CH2OH] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[1];
    J[645] += dqdci;              /* dwdot[CH2OH]/d[CH3O] */
    J[646] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[1074] += dqdT;              /* dwdot[CH2OH]/dT */
    J[1075] -= dqdT;              /* dwdot[CH3O]/dT */

    /*reaction 79: H + CH3O <=> H2 + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[19];
    k_f = prefactor_units[78] * fwd_A[78]
                * exp(fwd_beta[78] * tc[0] - activation_units[78] * fwd_Ea[78] * invT);
    dlnkfdT = fwd_beta[78] * invT + activation_units[78] * fwd_Ea[78] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[17];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[17] + g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[19]) + (h_RT[0] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[17] += q; /* CH2O */
    wdot[19] -= q; /* CH3O */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[17];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[17] += dqdci;               /* dwdot[CH2O]/d[H2] */
    J[19] -= dqdci;               /* dwdot[CH3O]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[19];
    J[33] += dqdci;               /* dwdot[H2]/d[H] */
    J[34] -= dqdci;               /* dwdot[H]/d[H] */
    J[50] += dqdci;               /* dwdot[CH2O]/d[H] */
    J[52] -= dqdci;               /* dwdot[CH3O]/d[H] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[0];
    J[561] += dqdci;              /* dwdot[H2]/d[CH2O] */
    J[562] -= dqdci;              /* dwdot[H]/d[CH2O] */
    J[578] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[580] -= dqdci;              /* dwdot[CH3O]/d[CH2O] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[1];
    J[627] += dqdci;              /* dwdot[H2]/d[CH3O] */
    J[628] -= dqdci;              /* dwdot[H]/d[CH3O] */
    J[644] += dqdci;              /* dwdot[CH2O]/d[CH3O] */
    J[646] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[1056] += dqdT;              /* dwdot[H2]/dT */
    J[1057] -= dqdT;              /* dwdot[H]/dT */
    J[1073] += dqdT;              /* dwdot[CH2O]/dT */
    J[1075] -= dqdT;              /* dwdot[CH3O]/dT */

    /*reaction 80: H + CH3O <=> OH + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[19];
    k_f = prefactor_units[79] * fwd_A[79]
                * exp(fwd_beta[79] * tc[0] - activation_units[79] * fwd_Ea[79] * invT);
    dlnkfdT = fwd_beta[79] * invT + activation_units[79] * fwd_Ea[79] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[12];
    Kc = exp(g_RT[1] - g_RT[4] - g_RT[12] + g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[19]) + (h_RT[4] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[4] += q; /* OH */
    wdot[12] += q; /* CH3 */
    wdot[19] -= q; /* CH3O */
    /* d()/d[H] */
    dqdci =  + k_f*sc[19];
    J[34] -= dqdci;               /* dwdot[H]/d[H] */
    J[37] += dqdci;               /* dwdot[OH]/d[H] */
    J[45] += dqdci;               /* dwdot[CH3]/d[H] */
    J[52] -= dqdci;               /* dwdot[CH3O]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[12];
    J[133] -= dqdci;              /* dwdot[H]/d[OH] */
    J[136] += dqdci;              /* dwdot[OH]/d[OH] */
    J[144] += dqdci;              /* dwdot[CH3]/d[OH] */
    J[151] -= dqdci;              /* dwdot[CH3O]/d[OH] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[4];
    J[397] -= dqdci;              /* dwdot[H]/d[CH3] */
    J[400] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[408] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[415] -= dqdci;              /* dwdot[CH3O]/d[CH3] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[1];
    J[628] -= dqdci;              /* dwdot[H]/d[CH3O] */
    J[631] += dqdci;              /* dwdot[OH]/d[CH3O] */
    J[639] += dqdci;              /* dwdot[CH3]/d[CH3O] */
    J[646] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[1057] -= dqdT;              /* dwdot[H]/dT */
    J[1060] += dqdT;              /* dwdot[OH]/dT */
    J[1068] += dqdT;              /* dwdot[CH3]/dT */
    J[1075] -= dqdT;              /* dwdot[CH3O]/dT */

    /*reaction 81: H + CH3O <=> CH2(S) + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[19];
    k_f = prefactor_units[80] * fwd_A[80]
                * exp(fwd_beta[80] * tc[0] - activation_units[80] * fwd_Ea[80] * invT);
    dlnkfdT = fwd_beta[80] * invT + activation_units[80] * fwd_Ea[80] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[11];
    Kc = exp(g_RT[1] - g_RT[5] - g_RT[11] + g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[19]) + (h_RT[5] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[5] += q; /* H2O */
    wdot[11] += q; /* CH2(S) */
    wdot[19] -= q; /* CH3O */
    /* d()/d[H] */
    dqdci =  + k_f*sc[19];
    J[34] -= dqdci;               /* dwdot[H]/d[H] */
    J[38] += dqdci;               /* dwdot[H2O]/d[H] */
    J[44] += dqdci;               /* dwdot[CH2(S)]/d[H] */
    J[52] -= dqdci;               /* dwdot[CH3O]/d[H] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[11];
    J[166] -= dqdci;              /* dwdot[H]/d[H2O] */
    J[170] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[176] += dqdci;              /* dwdot[CH2(S)]/d[H2O] */
    J[184] -= dqdci;              /* dwdot[CH3O]/d[H2O] */
    /* d()/d[CH2(S)] */
    dqdci =  - k_r*sc[5];
    J[364] -= dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[368] += dqdci;              /* dwdot[H2O]/d[CH2(S)] */
    J[374] += dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[382] -= dqdci;              /* dwdot[CH3O]/d[CH2(S)] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[1];
    J[628] -= dqdci;              /* dwdot[H]/d[CH3O] */
    J[632] += dqdci;              /* dwdot[H2O]/d[CH3O] */
    J[638] += dqdci;              /* dwdot[CH2(S)]/d[CH3O] */
    J[646] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[1057] -= dqdT;              /* dwdot[H]/dT */
    J[1061] += dqdT;              /* dwdot[H2O]/dT */
    J[1067] += dqdT;              /* dwdot[CH2(S)]/dT */
    J[1075] -= dqdT;              /* dwdot[CH3O]/dT */

    /*reaction 82: H + CH3OH <=> CH2OH + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[20];
    k_f = prefactor_units[81] * fwd_A[81]
                * exp(fwd_beta[81] * tc[0] - activation_units[81] * fwd_Ea[81] * invT);
    dlnkfdT = fwd_beta[81] * invT + activation_units[81] * fwd_Ea[81] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[18];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[18] + g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[20]) + (h_RT[0] + h_RT[18]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[18] += q; /* CH2OH */
    wdot[20] -= q; /* CH3OH */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[18];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[18] += dqdci;               /* dwdot[CH2OH]/d[H2] */
    J[20] -= dqdci;               /* dwdot[CH3OH]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[20];
    J[33] += dqdci;               /* dwdot[H2]/d[H] */
    J[34] -= dqdci;               /* dwdot[H]/d[H] */
    J[51] += dqdci;               /* dwdot[CH2OH]/d[H] */
    J[53] -= dqdci;               /* dwdot[CH3OH]/d[H] */
    /* d()/d[CH2OH] */
    dqdci =  - k_r*sc[0];
    J[594] += dqdci;              /* dwdot[H2]/d[CH2OH] */
    J[595] -= dqdci;              /* dwdot[H]/d[CH2OH] */
    J[612] += dqdci;              /* dwdot[CH2OH]/d[CH2OH] */
    J[614] -= dqdci;              /* dwdot[CH3OH]/d[CH2OH] */
    /* d()/d[CH3OH] */
    dqdci =  + k_f*sc[1];
    J[660] += dqdci;              /* dwdot[H2]/d[CH3OH] */
    J[661] -= dqdci;              /* dwdot[H]/d[CH3OH] */
    J[678] += dqdci;              /* dwdot[CH2OH]/d[CH3OH] */
    J[680] -= dqdci;              /* dwdot[CH3OH]/d[CH3OH] */
    /* d()/dT */
    J[1056] += dqdT;              /* dwdot[H2]/dT */
    J[1057] -= dqdT;              /* dwdot[H]/dT */
    J[1074] += dqdT;              /* dwdot[CH2OH]/dT */
    J[1076] -= dqdT;              /* dwdot[CH3OH]/dT */

    /*reaction 83: H + CH3OH <=> CH3O + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[20];
    k_f = prefactor_units[82] * fwd_A[82]
                * exp(fwd_beta[82] * tc[0] - activation_units[82] * fwd_Ea[82] * invT);
    dlnkfdT = fwd_beta[82] * invT + activation_units[82] * fwd_Ea[82] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[19];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[19] + g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[20]) + (h_RT[0] + h_RT[19]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[19] += q; /* CH3O */
    wdot[20] -= q; /* CH3OH */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[19];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[19] += dqdci;               /* dwdot[CH3O]/d[H2] */
    J[20] -= dqdci;               /* dwdot[CH3OH]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[20];
    J[33] += dqdci;               /* dwdot[H2]/d[H] */
    J[34] -= dqdci;               /* dwdot[H]/d[H] */
    J[52] += dqdci;               /* dwdot[CH3O]/d[H] */
    J[53] -= dqdci;               /* dwdot[CH3OH]/d[H] */
    /* d()/d[CH3O] */
    dqdci =  - k_r*sc[0];
    J[627] += dqdci;              /* dwdot[H2]/d[CH3O] */
    J[628] -= dqdci;              /* dwdot[H]/d[CH3O] */
    J[646] += dqdci;              /* dwdot[CH3O]/d[CH3O] */
    J[647] -= dqdci;              /* dwdot[CH3OH]/d[CH3O] */
    /* d()/d[CH3OH] */
    dqdci =  + k_f*sc[1];
    J[660] += dqdci;              /* dwdot[H2]/d[CH3OH] */
    J[661] -= dqdci;              /* dwdot[H]/d[CH3OH] */
    J[679] += dqdci;              /* dwdot[CH3O]/d[CH3OH] */
    J[680] -= dqdci;              /* dwdot[CH3OH]/d[CH3OH] */
    /* d()/dT */
    J[1056] += dqdT;              /* dwdot[H2]/dT */
    J[1057] -= dqdT;              /* dwdot[H]/dT */
    J[1075] += dqdT;              /* dwdot[CH3O]/dT */
    J[1076] -= dqdT;              /* dwdot[CH3OH]/dT */

    /*reaction 84: H + C2H3 <=> H2 + C2H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[23];
    k_f = prefactor_units[83] * fwd_A[83]
                * exp(fwd_beta[83] * tc[0] - activation_units[83] * fwd_Ea[83] * invT);
    dlnkfdT = fwd_beta[83] * invT + activation_units[83] * fwd_Ea[83] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[22];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[22] + g_RT[23]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[23]) + (h_RT[0] + h_RT[22]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[22] += q; /* C2H2 */
    wdot[23] -= q; /* C2H3 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[22];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[22] += dqdci;               /* dwdot[C2H2]/d[H2] */
    J[23] -= dqdci;               /* dwdot[C2H3]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[23];
    J[33] += dqdci;               /* dwdot[H2]/d[H] */
    J[34] -= dqdci;               /* dwdot[H]/d[H] */
    J[55] += dqdci;               /* dwdot[C2H2]/d[H] */
    J[56] -= dqdci;               /* dwdot[C2H3]/d[H] */
    /* d()/d[C2H2] */
    dqdci =  - k_r*sc[0];
    J[726] += dqdci;              /* dwdot[H2]/d[C2H2] */
    J[727] -= dqdci;              /* dwdot[H]/d[C2H2] */
    J[748] += dqdci;              /* dwdot[C2H2]/d[C2H2] */
    J[749] -= dqdci;              /* dwdot[C2H3]/d[C2H2] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[1];
    J[759] += dqdci;              /* dwdot[H2]/d[C2H3] */
    J[760] -= dqdci;              /* dwdot[H]/d[C2H3] */
    J[781] += dqdci;              /* dwdot[C2H2]/d[C2H3] */
    J[782] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    /* d()/dT */
    J[1056] += dqdT;              /* dwdot[H2]/dT */
    J[1057] -= dqdT;              /* dwdot[H]/dT */
    J[1078] += dqdT;              /* dwdot[C2H2]/dT */
    J[1079] -= dqdT;              /* dwdot[C2H3]/dT */

    /*reaction 85: H + C2H4 <=> C2H3 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[24];
    k_f = prefactor_units[84] * fwd_A[84]
                * exp(fwd_beta[84] * tc[0] - activation_units[84] * fwd_Ea[84] * invT);
    dlnkfdT = fwd_beta[84] * invT + activation_units[84] * fwd_Ea[84] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[23];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[23] + g_RT[24]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[24]) + (h_RT[0] + h_RT[23]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[23] += q; /* C2H3 */
    wdot[24] -= q; /* C2H4 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[23];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[23] += dqdci;               /* dwdot[C2H3]/d[H2] */
    J[24] -= dqdci;               /* dwdot[C2H4]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[24];
    J[33] += dqdci;               /* dwdot[H2]/d[H] */
    J[34] -= dqdci;               /* dwdot[H]/d[H] */
    J[56] += dqdci;               /* dwdot[C2H3]/d[H] */
    J[57] -= dqdci;               /* dwdot[C2H4]/d[H] */
    /* d()/d[C2H3] */
    dqdci =  - k_r*sc[0];
    J[759] += dqdci;              /* dwdot[H2]/d[C2H3] */
    J[760] -= dqdci;              /* dwdot[H]/d[C2H3] */
    J[782] += dqdci;              /* dwdot[C2H3]/d[C2H3] */
    J[783] -= dqdci;              /* dwdot[C2H4]/d[C2H3] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[1];
    J[792] += dqdci;              /* dwdot[H2]/d[C2H4] */
    J[793] -= dqdci;              /* dwdot[H]/d[C2H4] */
    J[815] += dqdci;              /* dwdot[C2H3]/d[C2H4] */
    J[816] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    /* d()/dT */
    J[1056] += dqdT;              /* dwdot[H2]/dT */
    J[1057] -= dqdT;              /* dwdot[H]/dT */
    J[1079] += dqdT;              /* dwdot[C2H3]/dT */
    J[1080] -= dqdT;              /* dwdot[C2H4]/dT */

    /*reaction 86: H + C2H5 <=> H2 + C2H4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[25];
    k_f = prefactor_units[85] * fwd_A[85]
                * exp(fwd_beta[85] * tc[0] - activation_units[85] * fwd_Ea[85] * invT);
    dlnkfdT = fwd_beta[85] * invT + activation_units[85] * fwd_Ea[85] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[24];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[24] + g_RT[25]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[25]) + (h_RT[0] + h_RT[24]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[24] += q; /* C2H4 */
    wdot[25] -= q; /* C2H5 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[24];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[24] += dqdci;               /* dwdot[C2H4]/d[H2] */
    J[25] -= dqdci;               /* dwdot[C2H5]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[25];
    J[33] += dqdci;               /* dwdot[H2]/d[H] */
    J[34] -= dqdci;               /* dwdot[H]/d[H] */
    J[57] += dqdci;               /* dwdot[C2H4]/d[H] */
    J[58] -= dqdci;               /* dwdot[C2H5]/d[H] */
    /* d()/d[C2H4] */
    dqdci =  - k_r*sc[0];
    J[792] += dqdci;              /* dwdot[H2]/d[C2H4] */
    J[793] -= dqdci;              /* dwdot[H]/d[C2H4] */
    J[816] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[817] -= dqdci;              /* dwdot[C2H5]/d[C2H4] */
    /* d()/d[C2H5] */
    dqdci =  + k_f*sc[1];
    J[825] += dqdci;              /* dwdot[H2]/d[C2H5] */
    J[826] -= dqdci;              /* dwdot[H]/d[C2H5] */
    J[849] += dqdci;              /* dwdot[C2H4]/d[C2H5] */
    J[850] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
    /* d()/dT */
    J[1056] += dqdT;              /* dwdot[H2]/dT */
    J[1057] -= dqdT;              /* dwdot[H]/dT */
    J[1080] += dqdT;              /* dwdot[C2H4]/dT */
    J[1081] -= dqdT;              /* dwdot[C2H5]/dT */

    /*reaction 87: H + C2H6 <=> C2H5 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[26];
    k_f = prefactor_units[86] * fwd_A[86]
                * exp(fwd_beta[86] * tc[0] - activation_units[86] * fwd_Ea[86] * invT);
    dlnkfdT = fwd_beta[86] * invT + activation_units[86] * fwd_Ea[86] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[25];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[25] + g_RT[26]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[26]) + (h_RT[0] + h_RT[25]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[25] += q; /* C2H5 */
    wdot[26] -= q; /* C2H6 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[25];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[25] += dqdci;               /* dwdot[C2H5]/d[H2] */
    J[26] -= dqdci;               /* dwdot[C2H6]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[26];
    J[33] += dqdci;               /* dwdot[H2]/d[H] */
    J[34] -= dqdci;               /* dwdot[H]/d[H] */
    J[58] += dqdci;               /* dwdot[C2H5]/d[H] */
    J[59] -= dqdci;               /* dwdot[C2H6]/d[H] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[0];
    J[825] += dqdci;              /* dwdot[H2]/d[C2H5] */
    J[826] -= dqdci;              /* dwdot[H]/d[C2H5] */
    J[850] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[851] -= dqdci;              /* dwdot[C2H6]/d[C2H5] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[1];
    J[858] += dqdci;              /* dwdot[H2]/d[C2H6] */
    J[859] -= dqdci;              /* dwdot[H]/d[C2H6] */
    J[883] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[884] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[1056] += dqdT;              /* dwdot[H2]/dT */
    J[1057] -= dqdT;              /* dwdot[H]/dT */
    J[1081] += dqdT;              /* dwdot[C2H5]/dT */
    J[1082] -= dqdT;              /* dwdot[C2H6]/dT */

    /*reaction 88: H + HCCO <=> CH2(S) + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[27];
    k_f = prefactor_units[87] * fwd_A[87]
                * exp(fwd_beta[87] * tc[0] - activation_units[87] * fwd_Ea[87] * invT);
    dlnkfdT = fwd_beta[87] * invT + activation_units[87] * fwd_Ea[87] * invT2;
    /* reverse */
    phi_r = sc[11]*sc[14];
    Kc = exp(g_RT[1] - g_RT[11] - g_RT[14] + g_RT[27]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[27]) + (h_RT[11] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[11] += q; /* CH2(S) */
    wdot[14] += q; /* CO */
    wdot[27] -= q; /* HCCO */
    /* d()/d[H] */
    dqdci =  + k_f*sc[27];
    J[34] -= dqdci;               /* dwdot[H]/d[H] */
    J[44] += dqdci;               /* dwdot[CH2(S)]/d[H] */
    J[47] += dqdci;               /* dwdot[CO]/d[H] */
    J[60] -= dqdci;               /* dwdot[HCCO]/d[H] */
    /* d()/d[CH2(S)] */
    dqdci =  - k_r*sc[14];
    J[364] -= dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[374] += dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[377] += dqdci;              /* dwdot[CO]/d[CH2(S)] */
    J[390] -= dqdci;              /* dwdot[HCCO]/d[CH2(S)] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[11];
    J[463] -= dqdci;              /* dwdot[H]/d[CO] */
    J[473] += dqdci;              /* dwdot[CH2(S)]/d[CO] */
    J[476] += dqdci;              /* dwdot[CO]/d[CO] */
    J[489] -= dqdci;              /* dwdot[HCCO]/d[CO] */
    /* d()/d[HCCO] */
    dqdci =  + k_f*sc[1];
    J[892] -= dqdci;              /* dwdot[H]/d[HCCO] */
    J[902] += dqdci;              /* dwdot[CH2(S)]/d[HCCO] */
    J[905] += dqdci;              /* dwdot[CO]/d[HCCO] */
    J[918] -= dqdci;              /* dwdot[HCCO]/d[HCCO] */
    /* d()/dT */
    J[1057] -= dqdT;              /* dwdot[H]/dT */
    J[1067] += dqdT;              /* dwdot[CH2(S)]/dT */
    J[1070] += dqdT;              /* dwdot[CO]/dT */
    J[1083] -= dqdT;              /* dwdot[HCCO]/dT */

    /*reaction 89: H + CH2CO <=> HCCO + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[28];
    k_f = prefactor_units[88] * fwd_A[88]
                * exp(fwd_beta[88] * tc[0] - activation_units[88] * fwd_Ea[88] * invT);
    dlnkfdT = fwd_beta[88] * invT + activation_units[88] * fwd_Ea[88] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[27];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[27] + g_RT[28]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[28]) + (h_RT[0] + h_RT[27]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[27] += q; /* HCCO */
    wdot[28] -= q; /* CH2CO */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[27];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[27] += dqdci;               /* dwdot[HCCO]/d[H2] */
    J[28] -= dqdci;               /* dwdot[CH2CO]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[28];
    J[33] += dqdci;               /* dwdot[H2]/d[H] */
    J[34] -= dqdci;               /* dwdot[H]/d[H] */
    J[60] += dqdci;               /* dwdot[HCCO]/d[H] */
    J[61] -= dqdci;               /* dwdot[CH2CO]/d[H] */
    /* d()/d[HCCO] */
    dqdci =  - k_r*sc[0];
    J[891] += dqdci;              /* dwdot[H2]/d[HCCO] */
    J[892] -= dqdci;              /* dwdot[H]/d[HCCO] */
    J[918] += dqdci;              /* dwdot[HCCO]/d[HCCO] */
    J[919] -= dqdci;              /* dwdot[CH2CO]/d[HCCO] */
    /* d()/d[CH2CO] */
    dqdci =  + k_f*sc[1];
    J[924] += dqdci;              /* dwdot[H2]/d[CH2CO] */
    J[925] -= dqdci;              /* dwdot[H]/d[CH2CO] */
    J[951] += dqdci;              /* dwdot[HCCO]/d[CH2CO] */
    J[952] -= dqdci;              /* dwdot[CH2CO]/d[CH2CO] */
    /* d()/dT */
    J[1056] += dqdT;              /* dwdot[H2]/dT */
    J[1057] -= dqdT;              /* dwdot[H]/dT */
    J[1083] += dqdT;              /* dwdot[HCCO]/dT */
    J[1084] -= dqdT;              /* dwdot[CH2CO]/dT */

    /*reaction 90: H + CH2CO <=> CH3 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[28];
    k_f = prefactor_units[89] * fwd_A[89]
                * exp(fwd_beta[89] * tc[0] - activation_units[89] * fwd_Ea[89] * invT);
    dlnkfdT = fwd_beta[89] * invT + activation_units[89] * fwd_Ea[89] * invT2;
    /* reverse */
    phi_r = sc[12]*sc[14];
    Kc = exp(g_RT[1] - g_RT[12] - g_RT[14] + g_RT[28]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[28]) + (h_RT[12] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[12] += q; /* CH3 */
    wdot[14] += q; /* CO */
    wdot[28] -= q; /* CH2CO */
    /* d()/d[H] */
    dqdci =  + k_f*sc[28];
    J[34] -= dqdci;               /* dwdot[H]/d[H] */
    J[45] += dqdci;               /* dwdot[CH3]/d[H] */
    J[47] += dqdci;               /* dwdot[CO]/d[H] */
    J[61] -= dqdci;               /* dwdot[CH2CO]/d[H] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[14];
    J[397] -= dqdci;              /* dwdot[H]/d[CH3] */
    J[408] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[410] += dqdci;              /* dwdot[CO]/d[CH3] */
    J[424] -= dqdci;              /* dwdot[CH2CO]/d[CH3] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[12];
    J[463] -= dqdci;              /* dwdot[H]/d[CO] */
    J[474] += dqdci;              /* dwdot[CH3]/d[CO] */
    J[476] += dqdci;              /* dwdot[CO]/d[CO] */
    J[490] -= dqdci;              /* dwdot[CH2CO]/d[CO] */
    /* d()/d[CH2CO] */
    dqdci =  + k_f*sc[1];
    J[925] -= dqdci;              /* dwdot[H]/d[CH2CO] */
    J[936] += dqdci;              /* dwdot[CH3]/d[CH2CO] */
    J[938] += dqdci;              /* dwdot[CO]/d[CH2CO] */
    J[952] -= dqdci;              /* dwdot[CH2CO]/d[CH2CO] */
    /* d()/dT */
    J[1057] -= dqdT;              /* dwdot[H]/dT */
    J[1068] += dqdT;              /* dwdot[CH3]/dT */
    J[1070] += dqdT;              /* dwdot[CO]/dT */
    J[1084] -= dqdT;              /* dwdot[CH2CO]/dT */

    /*reaction 91: H + HCCOH <=> H + CH2CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[29];
    k_f = prefactor_units[90] * fwd_A[90]
                * exp(fwd_beta[90] * tc[0] - activation_units[90] * fwd_Ea[90] * invT);
    dlnkfdT = fwd_beta[90] * invT + activation_units[90] * fwd_Ea[90] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[28];
    Kc = exp(g_RT[1] - g_RT[1] - g_RT[28] + g_RT[29]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[29]) + (h_RT[1] + h_RT[28]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[28] += q; /* CH2CO */
    wdot[29] -= q; /* HCCOH */
    /* d()/d[H] */
    dqdci =  + k_f*sc[29] - k_r*sc[28];
    J[61] += dqdci;               /* dwdot[CH2CO]/d[H] */
    J[62] -= dqdci;               /* dwdot[HCCOH]/d[H] */
    /* d()/d[CH2CO] */
    dqdci =  - k_r*sc[1];
    J[952] += dqdci;              /* dwdot[CH2CO]/d[CH2CO] */
    J[953] -= dqdci;              /* dwdot[HCCOH]/d[CH2CO] */
    /* d()/d[HCCOH] */
    dqdci =  + k_f*sc[1];
    J[985] += dqdci;              /* dwdot[CH2CO]/d[HCCOH] */
    J[986] -= dqdci;              /* dwdot[HCCOH]/d[HCCOH] */
    /* d()/dT */
    J[1084] += dqdT;              /* dwdot[CH2CO]/dT */
    J[1085] -= dqdT;              /* dwdot[HCCOH]/dT */

    /*reaction 92: OH + H2 <=> H + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[4];
    k_f = prefactor_units[91] * fwd_A[91]
                * exp(fwd_beta[91] * tc[0] - activation_units[91] * fwd_Ea[91] * invT);
    dlnkfdT = fwd_beta[91] * invT + activation_units[91] * fwd_Ea[91] * invT2;
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
    J[33] -= dqdci;               /* dwdot[H2]/d[H] */
    J[34] += dqdci;               /* dwdot[H]/d[H] */
    J[37] -= dqdci;               /* dwdot[OH]/d[H] */
    J[38] += dqdci;               /* dwdot[H2O]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[0];
    J[132] -= dqdci;              /* dwdot[H2]/d[OH] */
    J[133] += dqdci;              /* dwdot[H]/d[OH] */
    J[136] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[137] += dqdci;              /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[1];
    J[165] -= dqdci;              /* dwdot[H2]/d[H2O] */
    J[166] += dqdci;              /* dwdot[H]/d[H2O] */
    J[169] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[170] += dqdci;              /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[1056] -= dqdT;              /* dwdot[H2]/dT */
    J[1057] += dqdT;              /* dwdot[H]/dT */
    J[1060] -= dqdT;              /* dwdot[OH]/dT */
    J[1061] += dqdT;              /* dwdot[H2O]/dT */

    /*reaction 93: 2.000000 OH <=> O + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[4], 2.000000);
    k_f = prefactor_units[92] * fwd_A[92]
                * exp(fwd_beta[92] * tc[0] - activation_units[92] * fwd_Ea[92] * invT);
    dlnkfdT = fwd_beta[92] * invT + activation_units[92] * fwd_Ea[92] * invT2;
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
    J[68] += dqdci;               /* dwdot[O]/d[O] */
    J[70] += -2 * dqdci;          /* dwdot[OH]/d[O] */
    J[71] += dqdci;               /* dwdot[H2O]/d[O] */
    /* d()/d[OH] */
    dqdci =  + k_f*2.000000*sc[4];
    J[134] += dqdci;              /* dwdot[O]/d[OH] */
    J[136] += -2 * dqdci;         /* dwdot[OH]/d[OH] */
    J[137] += dqdci;              /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[2];
    J[167] += dqdci;              /* dwdot[O]/d[H2O] */
    J[169] += -2 * dqdci;         /* dwdot[OH]/d[H2O] */
    J[170] += dqdci;              /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[1058] += dqdT;              /* dwdot[O]/dT */
    J[1060] += -2 * dqdT;         /* dwdot[OH]/dT */
    J[1061] += dqdT;              /* dwdot[H2O]/dT */

    /*reaction 94: OH + HO2 <=> O2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[6];
    k_f = prefactor_units[93] * fwd_A[93]
                * exp(fwd_beta[93] * tc[0] - activation_units[93] * fwd_Ea[93] * invT);
    dlnkfdT = fwd_beta[93] * invT + activation_units[93] * fwd_Ea[93] * invT2;
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
    J[102] += dqdci;              /* dwdot[O2]/d[O2] */
    J[103] -= dqdci;              /* dwdot[OH]/d[O2] */
    J[104] += dqdci;              /* dwdot[H2O]/d[O2] */
    J[105] -= dqdci;              /* dwdot[HO2]/d[O2] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[6];
    J[135] += dqdci;              /* dwdot[O2]/d[OH] */
    J[136] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[137] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[138] -= dqdci;              /* dwdot[HO2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[3];
    J[168] += dqdci;              /* dwdot[O2]/d[H2O] */
    J[169] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[170] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[171] -= dqdci;              /* dwdot[HO2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[4];
    J[201] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[202] -= dqdci;              /* dwdot[OH]/d[HO2] */
    J[203] += dqdci;              /* dwdot[H2O]/d[HO2] */
    J[204] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1059] += dqdT;              /* dwdot[O2]/dT */
    J[1060] -= dqdT;              /* dwdot[OH]/dT */
    J[1061] += dqdT;              /* dwdot[H2O]/dT */
    J[1062] -= dqdT;              /* dwdot[HO2]/dT */

    /*reaction 95: OH + H2O2 <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[7];
    k_f = prefactor_units[94] * fwd_A[94]
                * exp(fwd_beta[94] * tc[0] - activation_units[94] * fwd_Ea[94] * invT);
    dlnkfdT = fwd_beta[94] * invT + activation_units[94] * fwd_Ea[94] * invT2;
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
    J[136] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[137] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[138] += dqdci;              /* dwdot[HO2]/d[OH] */
    J[139] -= dqdci;              /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[6];
    J[169] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[170] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[171] += dqdci;              /* dwdot[HO2]/d[H2O] */
    J[172] -= dqdci;              /* dwdot[H2O2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[5];
    J[202] -= dqdci;              /* dwdot[OH]/d[HO2] */
    J[203] += dqdci;              /* dwdot[H2O]/d[HO2] */
    J[204] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[205] -= dqdci;              /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[4];
    J[235] -= dqdci;              /* dwdot[OH]/d[H2O2] */
    J[236] += dqdci;              /* dwdot[H2O]/d[H2O2] */
    J[237] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[238] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[1060] -= dqdT;              /* dwdot[OH]/dT */
    J[1061] += dqdT;              /* dwdot[H2O]/dT */
    J[1062] += dqdT;              /* dwdot[HO2]/dT */
    J[1063] -= dqdT;              /* dwdot[H2O2]/dT */

    /*reaction 96: OH + H2O2 <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[7];
    k_f = prefactor_units[95] * fwd_A[95]
                * exp(fwd_beta[95] * tc[0] - activation_units[95] * fwd_Ea[95] * invT);
    dlnkfdT = fwd_beta[95] * invT + activation_units[95] * fwd_Ea[95] * invT2;
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
    J[136] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[137] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[138] += dqdci;              /* dwdot[HO2]/d[OH] */
    J[139] -= dqdci;              /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[6];
    J[169] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[170] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[171] += dqdci;              /* dwdot[HO2]/d[H2O] */
    J[172] -= dqdci;              /* dwdot[H2O2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[5];
    J[202] -= dqdci;              /* dwdot[OH]/d[HO2] */
    J[203] += dqdci;              /* dwdot[H2O]/d[HO2] */
    J[204] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[205] -= dqdci;              /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[4];
    J[235] -= dqdci;              /* dwdot[OH]/d[H2O2] */
    J[236] += dqdci;              /* dwdot[H2O]/d[H2O2] */
    J[237] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[238] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[1060] -= dqdT;              /* dwdot[OH]/dT */
    J[1061] += dqdT;              /* dwdot[H2O]/dT */
    J[1062] += dqdT;              /* dwdot[HO2]/dT */
    J[1063] -= dqdT;              /* dwdot[H2O2]/dT */

    /*reaction 97: OH + C <=> H + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[8];
    k_f = prefactor_units[96] * fwd_A[96]
                * exp(fwd_beta[96] * tc[0] - activation_units[96] * fwd_Ea[96] * invT);
    dlnkfdT = fwd_beta[96] * invT + activation_units[96] * fwd_Ea[96] * invT2;
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
    wdot[8] -= q; /* C */
    wdot[14] += q; /* CO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[14];
    J[34] += dqdci;               /* dwdot[H]/d[H] */
    J[37] -= dqdci;               /* dwdot[OH]/d[H] */
    J[41] -= dqdci;               /* dwdot[C]/d[H] */
    J[47] += dqdci;               /* dwdot[CO]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[8];
    J[133] += dqdci;              /* dwdot[H]/d[OH] */
    J[136] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[140] -= dqdci;              /* dwdot[C]/d[OH] */
    J[146] += dqdci;              /* dwdot[CO]/d[OH] */
    /* d()/d[C] */
    dqdci =  + k_f*sc[4];
    J[265] += dqdci;              /* dwdot[H]/d[C] */
    J[268] -= dqdci;              /* dwdot[OH]/d[C] */
    J[272] -= dqdci;              /* dwdot[C]/d[C] */
    J[278] += dqdci;              /* dwdot[CO]/d[C] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[1];
    J[463] += dqdci;              /* dwdot[H]/d[CO] */
    J[466] -= dqdci;              /* dwdot[OH]/d[CO] */
    J[470] -= dqdci;              /* dwdot[C]/d[CO] */
    J[476] += dqdci;              /* dwdot[CO]/d[CO] */
    /* d()/dT */
    J[1057] += dqdT;              /* dwdot[H]/dT */
    J[1060] -= dqdT;              /* dwdot[OH]/dT */
    J[1064] -= dqdT;              /* dwdot[C]/dT */
    J[1070] += dqdT;              /* dwdot[CO]/dT */

    /*reaction 98: OH + CH <=> H + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[9];
    k_f = prefactor_units[97] * fwd_A[97]
                * exp(fwd_beta[97] * tc[0] - activation_units[97] * fwd_Ea[97] * invT);
    dlnkfdT = fwd_beta[97] * invT + activation_units[97] * fwd_Ea[97] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[16];
    Kc = exp(-g_RT[1] + g_RT[4] + g_RT[9] - g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[9]) + (h_RT[1] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[4] -= q; /* OH */
    wdot[9] -= q; /* CH */
    wdot[16] += q; /* HCO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[16];
    J[34] += dqdci;               /* dwdot[H]/d[H] */
    J[37] -= dqdci;               /* dwdot[OH]/d[H] */
    J[42] -= dqdci;               /* dwdot[CH]/d[H] */
    J[49] += dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[9];
    J[133] += dqdci;              /* dwdot[H]/d[OH] */
    J[136] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[141] -= dqdci;              /* dwdot[CH]/d[OH] */
    J[148] += dqdci;              /* dwdot[HCO]/d[OH] */
    /* d()/d[CH] */
    dqdci =  + k_f*sc[4];
    J[298] += dqdci;              /* dwdot[H]/d[CH] */
    J[301] -= dqdci;              /* dwdot[OH]/d[CH] */
    J[306] -= dqdci;              /* dwdot[CH]/d[CH] */
    J[313] += dqdci;              /* dwdot[HCO]/d[CH] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[1];
    J[529] += dqdci;              /* dwdot[H]/d[HCO] */
    J[532] -= dqdci;              /* dwdot[OH]/d[HCO] */
    J[537] -= dqdci;              /* dwdot[CH]/d[HCO] */
    J[544] += dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[1057] += dqdT;              /* dwdot[H]/dT */
    J[1060] -= dqdT;              /* dwdot[OH]/dT */
    J[1065] -= dqdT;              /* dwdot[CH]/dT */
    J[1072] += dqdT;              /* dwdot[HCO]/dT */

    /*reaction 99: OH + CH2 <=> H + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[10];
    k_f = prefactor_units[98] * fwd_A[98]
                * exp(fwd_beta[98] * tc[0] - activation_units[98] * fwd_Ea[98] * invT);
    dlnkfdT = fwd_beta[98] * invT + activation_units[98] * fwd_Ea[98] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[17];
    Kc = exp(-g_RT[1] + g_RT[4] + g_RT[10] - g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[10]) + (h_RT[1] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[4] -= q; /* OH */
    wdot[10] -= q; /* CH2 */
    wdot[17] += q; /* CH2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[17];
    J[34] += dqdci;               /* dwdot[H]/d[H] */
    J[37] -= dqdci;               /* dwdot[OH]/d[H] */
    J[43] -= dqdci;               /* dwdot[CH2]/d[H] */
    J[50] += dqdci;               /* dwdot[CH2O]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[10];
    J[133] += dqdci;              /* dwdot[H]/d[OH] */
    J[136] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[142] -= dqdci;              /* dwdot[CH2]/d[OH] */
    J[149] += dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[4];
    J[331] += dqdci;              /* dwdot[H]/d[CH2] */
    J[334] -= dqdci;              /* dwdot[OH]/d[CH2] */
    J[340] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[347] += dqdci;              /* dwdot[CH2O]/d[CH2] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[1];
    J[562] += dqdci;              /* dwdot[H]/d[CH2O] */
    J[565] -= dqdci;              /* dwdot[OH]/d[CH2O] */
    J[571] -= dqdci;              /* dwdot[CH2]/d[CH2O] */
    J[578] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[1057] += dqdT;              /* dwdot[H]/dT */
    J[1060] -= dqdT;              /* dwdot[OH]/dT */
    J[1066] -= dqdT;              /* dwdot[CH2]/dT */
    J[1073] += dqdT;              /* dwdot[CH2O]/dT */

    /*reaction 100: OH + CH2 <=> CH + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[10];
    k_f = prefactor_units[99] * fwd_A[99]
                * exp(fwd_beta[99] * tc[0] - activation_units[99] * fwd_Ea[99] * invT);
    dlnkfdT = fwd_beta[99] * invT + activation_units[99] * fwd_Ea[99] * invT2;
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
    wdot[9] += q; /* CH */
    wdot[10] -= q; /* CH2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[10];
    J[136] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[137] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[141] += dqdci;              /* dwdot[CH]/d[OH] */
    J[142] -= dqdci;              /* dwdot[CH2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[9];
    J[169] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[170] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[174] += dqdci;              /* dwdot[CH]/d[H2O] */
    J[175] -= dqdci;              /* dwdot[CH2]/d[H2O] */
    /* d()/d[CH] */
    dqdci =  - k_r*sc[5];
    J[301] -= dqdci;              /* dwdot[OH]/d[CH] */
    J[302] += dqdci;              /* dwdot[H2O]/d[CH] */
    J[306] += dqdci;              /* dwdot[CH]/d[CH] */
    J[307] -= dqdci;              /* dwdot[CH2]/d[CH] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[4];
    J[334] -= dqdci;              /* dwdot[OH]/d[CH2] */
    J[335] += dqdci;              /* dwdot[H2O]/d[CH2] */
    J[339] += dqdci;              /* dwdot[CH]/d[CH2] */
    J[340] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    /* d()/dT */
    J[1060] -= dqdT;              /* dwdot[OH]/dT */
    J[1061] += dqdT;              /* dwdot[H2O]/dT */
    J[1065] += dqdT;              /* dwdot[CH]/dT */
    J[1066] -= dqdT;              /* dwdot[CH2]/dT */

    /*reaction 101: OH + CH2(S) <=> H + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[11];
    k_f = prefactor_units[100] * fwd_A[100]
                * exp(fwd_beta[100] * tc[0] - activation_units[100] * fwd_Ea[100] * invT);
    dlnkfdT = fwd_beta[100] * invT + activation_units[100] * fwd_Ea[100] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[17];
    Kc = exp(-g_RT[1] + g_RT[4] + g_RT[11] - g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[11]) + (h_RT[1] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[4] -= q; /* OH */
    wdot[11] -= q; /* CH2(S) */
    wdot[17] += q; /* CH2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[17];
    J[34] += dqdci;               /* dwdot[H]/d[H] */
    J[37] -= dqdci;               /* dwdot[OH]/d[H] */
    J[44] -= dqdci;               /* dwdot[CH2(S)]/d[H] */
    J[50] += dqdci;               /* dwdot[CH2O]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[11];
    J[133] += dqdci;              /* dwdot[H]/d[OH] */
    J[136] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[143] -= dqdci;              /* dwdot[CH2(S)]/d[OH] */
    J[149] += dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[4];
    J[364] += dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[367] -= dqdci;              /* dwdot[OH]/d[CH2(S)] */
    J[374] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[380] += dqdci;              /* dwdot[CH2O]/d[CH2(S)] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[1];
    J[562] += dqdci;              /* dwdot[H]/d[CH2O] */
    J[565] -= dqdci;              /* dwdot[OH]/d[CH2O] */
    J[572] -= dqdci;              /* dwdot[CH2(S)]/d[CH2O] */
    J[578] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[1057] += dqdT;              /* dwdot[H]/dT */
    J[1060] -= dqdT;              /* dwdot[OH]/dT */
    J[1067] -= dqdT;              /* dwdot[CH2(S)]/dT */
    J[1073] += dqdT;              /* dwdot[CH2O]/dT */

    /*reaction 102: OH + CH3 <=> CH2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[12];
    k_f = prefactor_units[101] * fwd_A[101]
                * exp(fwd_beta[101] * tc[0] - activation_units[101] * fwd_Ea[101] * invT);
    dlnkfdT = fwd_beta[101] * invT + activation_units[101] * fwd_Ea[101] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[10];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[10] + g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[12]) + (h_RT[5] + h_RT[10]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[10] += q; /* CH2 */
    wdot[12] -= q; /* CH3 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[12];
    J[136] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[137] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[142] += dqdci;              /* dwdot[CH2]/d[OH] */
    J[144] -= dqdci;              /* dwdot[CH3]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[10];
    J[169] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[170] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[175] += dqdci;              /* dwdot[CH2]/d[H2O] */
    J[177] -= dqdci;              /* dwdot[CH3]/d[H2O] */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[5];
    J[334] -= dqdci;              /* dwdot[OH]/d[CH2] */
    J[335] += dqdci;              /* dwdot[H2O]/d[CH2] */
    J[340] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[342] -= dqdci;              /* dwdot[CH3]/d[CH2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[4];
    J[400] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[401] += dqdci;              /* dwdot[H2O]/d[CH3] */
    J[406] += dqdci;              /* dwdot[CH2]/d[CH3] */
    J[408] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[1060] -= dqdT;              /* dwdot[OH]/dT */
    J[1061] += dqdT;              /* dwdot[H2O]/dT */
    J[1066] += dqdT;              /* dwdot[CH2]/dT */
    J[1068] -= dqdT;              /* dwdot[CH3]/dT */

    /*reaction 103: OH + CH3 <=> CH2(S) + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[12];
    k_f = prefactor_units[102] * fwd_A[102]
                * exp(fwd_beta[102] * tc[0] - activation_units[102] * fwd_Ea[102] * invT);
    dlnkfdT = fwd_beta[102] * invT + activation_units[102] * fwd_Ea[102] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[11];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[11] + g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[12]) + (h_RT[5] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[11] += q; /* CH2(S) */
    wdot[12] -= q; /* CH3 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[12];
    J[136] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[137] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[143] += dqdci;              /* dwdot[CH2(S)]/d[OH] */
    J[144] -= dqdci;              /* dwdot[CH3]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[11];
    J[169] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[170] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[176] += dqdci;              /* dwdot[CH2(S)]/d[H2O] */
    J[177] -= dqdci;              /* dwdot[CH3]/d[H2O] */
    /* d()/d[CH2(S)] */
    dqdci =  - k_r*sc[5];
    J[367] -= dqdci;              /* dwdot[OH]/d[CH2(S)] */
    J[368] += dqdci;              /* dwdot[H2O]/d[CH2(S)] */
    J[374] += dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[375] -= dqdci;              /* dwdot[CH3]/d[CH2(S)] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[4];
    J[400] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[401] += dqdci;              /* dwdot[H2O]/d[CH3] */
    J[407] += dqdci;              /* dwdot[CH2(S)]/d[CH3] */
    J[408] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[1060] -= dqdT;              /* dwdot[OH]/dT */
    J[1061] += dqdT;              /* dwdot[H2O]/dT */
    J[1067] += dqdT;              /* dwdot[CH2(S)]/dT */
    J[1068] -= dqdT;              /* dwdot[CH3]/dT */

    /*reaction 104: OH + CH4 <=> CH3 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[13];
    k_f = prefactor_units[103] * fwd_A[103]
                * exp(fwd_beta[103] * tc[0] - activation_units[103] * fwd_Ea[103] * invT);
    dlnkfdT = fwd_beta[103] * invT + activation_units[103] * fwd_Ea[103] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[12];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[12] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[13]) + (h_RT[5] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[12] += q; /* CH3 */
    wdot[13] -= q; /* CH4 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[13];
    J[136] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[137] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[144] += dqdci;              /* dwdot[CH3]/d[OH] */
    J[145] -= dqdci;              /* dwdot[CH4]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[12];
    J[169] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[170] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[177] += dqdci;              /* dwdot[CH3]/d[H2O] */
    J[178] -= dqdci;              /* dwdot[CH4]/d[H2O] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[5];
    J[400] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[401] += dqdci;              /* dwdot[H2O]/d[CH3] */
    J[408] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[409] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[4];
    J[433] -= dqdci;              /* dwdot[OH]/d[CH4] */
    J[434] += dqdci;              /* dwdot[H2O]/d[CH4] */
    J[441] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[442] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[1060] -= dqdT;              /* dwdot[OH]/dT */
    J[1061] += dqdT;              /* dwdot[H2O]/dT */
    J[1068] += dqdT;              /* dwdot[CH3]/dT */
    J[1069] -= dqdT;              /* dwdot[CH4]/dT */

    /*reaction 105: OH + CO <=> H + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[14];
    k_f = prefactor_units[104] * fwd_A[104]
                * exp(fwd_beta[104] * tc[0] - activation_units[104] * fwd_Ea[104] * invT);
    dlnkfdT = fwd_beta[104] * invT + activation_units[104] * fwd_Ea[104] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[15];
    Kc = exp(-g_RT[1] + g_RT[4] + g_RT[14] - g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[14]) + (h_RT[1] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[4] -= q; /* OH */
    wdot[14] -= q; /* CO */
    wdot[15] += q; /* CO2 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[15];
    J[34] += dqdci;               /* dwdot[H]/d[H] */
    J[37] -= dqdci;               /* dwdot[OH]/d[H] */
    J[47] -= dqdci;               /* dwdot[CO]/d[H] */
    J[48] += dqdci;               /* dwdot[CO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[14];
    J[133] += dqdci;              /* dwdot[H]/d[OH] */
    J[136] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[146] -= dqdci;              /* dwdot[CO]/d[OH] */
    J[147] += dqdci;              /* dwdot[CO2]/d[OH] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[4];
    J[463] += dqdci;              /* dwdot[H]/d[CO] */
    J[466] -= dqdci;              /* dwdot[OH]/d[CO] */
    J[476] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[477] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[1];
    J[496] += dqdci;              /* dwdot[H]/d[CO2] */
    J[499] -= dqdci;              /* dwdot[OH]/d[CO2] */
    J[509] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[510] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[1057] += dqdT;              /* dwdot[H]/dT */
    J[1060] -= dqdT;              /* dwdot[OH]/dT */
    J[1070] -= dqdT;              /* dwdot[CO]/dT */
    J[1071] += dqdT;              /* dwdot[CO2]/dT */

    /*reaction 106: OH + HCO <=> H2O + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[16];
    k_f = prefactor_units[105] * fwd_A[105]
                * exp(fwd_beta[105] * tc[0] - activation_units[105] * fwd_Ea[105] * invT);
    dlnkfdT = fwd_beta[105] * invT + activation_units[105] * fwd_Ea[105] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[14];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[14] + g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[16]) + (h_RT[5] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[14] += q; /* CO */
    wdot[16] -= q; /* HCO */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[16];
    J[136] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[137] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[146] += dqdci;              /* dwdot[CO]/d[OH] */
    J[148] -= dqdci;              /* dwdot[HCO]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[14];
    J[169] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[170] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[179] += dqdci;              /* dwdot[CO]/d[H2O] */
    J[181] -= dqdci;              /* dwdot[HCO]/d[H2O] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[5];
    J[466] -= dqdci;              /* dwdot[OH]/d[CO] */
    J[467] += dqdci;              /* dwdot[H2O]/d[CO] */
    J[476] += dqdci;              /* dwdot[CO]/d[CO] */
    J[478] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[4];
    J[532] -= dqdci;              /* dwdot[OH]/d[HCO] */
    J[533] += dqdci;              /* dwdot[H2O]/d[HCO] */
    J[542] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[544] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[1060] -= dqdT;              /* dwdot[OH]/dT */
    J[1061] += dqdT;              /* dwdot[H2O]/dT */
    J[1070] += dqdT;              /* dwdot[CO]/dT */
    J[1072] -= dqdT;              /* dwdot[HCO]/dT */

    /*reaction 107: OH + CH2O <=> HCO + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[17];
    k_f = prefactor_units[106] * fwd_A[106]
                * exp(fwd_beta[106] * tc[0] - activation_units[106] * fwd_Ea[106] * invT);
    dlnkfdT = fwd_beta[106] * invT + activation_units[106] * fwd_Ea[106] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[16];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[16] + g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[17]) + (h_RT[5] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[16] += q; /* HCO */
    wdot[17] -= q; /* CH2O */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[17];
    J[136] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[137] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[148] += dqdci;              /* dwdot[HCO]/d[OH] */
    J[149] -= dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[16];
    J[169] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[170] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[181] += dqdci;              /* dwdot[HCO]/d[H2O] */
    J[182] -= dqdci;              /* dwdot[CH2O]/d[H2O] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[5];
    J[532] -= dqdci;              /* dwdot[OH]/d[HCO] */
    J[533] += dqdci;              /* dwdot[H2O]/d[HCO] */
    J[544] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[545] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[4];
    J[565] -= dqdci;              /* dwdot[OH]/d[CH2O] */
    J[566] += dqdci;              /* dwdot[H2O]/d[CH2O] */
    J[577] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[578] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[1060] -= dqdT;              /* dwdot[OH]/dT */
    J[1061] += dqdT;              /* dwdot[H2O]/dT */
    J[1072] += dqdT;              /* dwdot[HCO]/dT */
    J[1073] -= dqdT;              /* dwdot[CH2O]/dT */

    /*reaction 108: OH + CH2OH <=> H2O + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[18];
    k_f = prefactor_units[107] * fwd_A[107]
                * exp(fwd_beta[107] * tc[0] - activation_units[107] * fwd_Ea[107] * invT);
    dlnkfdT = fwd_beta[107] * invT + activation_units[107] * fwd_Ea[107] * invT2;
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
    wdot[17] += q; /* CH2O */
    wdot[18] -= q; /* CH2OH */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[18];
    J[136] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[137] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[149] += dqdci;              /* dwdot[CH2O]/d[OH] */
    J[150] -= dqdci;              /* dwdot[CH2OH]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[17];
    J[169] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[170] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[182] += dqdci;              /* dwdot[CH2O]/d[H2O] */
    J[183] -= dqdci;              /* dwdot[CH2OH]/d[H2O] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[5];
    J[565] -= dqdci;              /* dwdot[OH]/d[CH2O] */
    J[566] += dqdci;              /* dwdot[H2O]/d[CH2O] */
    J[578] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[579] -= dqdci;              /* dwdot[CH2OH]/d[CH2O] */
    /* d()/d[CH2OH] */
    dqdci =  + k_f*sc[4];
    J[598] -= dqdci;              /* dwdot[OH]/d[CH2OH] */
    J[599] += dqdci;              /* dwdot[H2O]/d[CH2OH] */
    J[611] += dqdci;              /* dwdot[CH2O]/d[CH2OH] */
    J[612] -= dqdci;              /* dwdot[CH2OH]/d[CH2OH] */
    /* d()/dT */
    J[1060] -= dqdT;              /* dwdot[OH]/dT */
    J[1061] += dqdT;              /* dwdot[H2O]/dT */
    J[1073] += dqdT;              /* dwdot[CH2O]/dT */
    J[1074] -= dqdT;              /* dwdot[CH2OH]/dT */

    /*reaction 109: OH + CH3O <=> H2O + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[19];
    k_f = prefactor_units[108] * fwd_A[108]
                * exp(fwd_beta[108] * tc[0] - activation_units[108] * fwd_Ea[108] * invT);
    dlnkfdT = fwd_beta[108] * invT + activation_units[108] * fwd_Ea[108] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[17];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[17] + g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[19]) + (h_RT[5] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[17] += q; /* CH2O */
    wdot[19] -= q; /* CH3O */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[19];
    J[136] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[137] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[149] += dqdci;              /* dwdot[CH2O]/d[OH] */
    J[151] -= dqdci;              /* dwdot[CH3O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[17];
    J[169] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[170] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[182] += dqdci;              /* dwdot[CH2O]/d[H2O] */
    J[184] -= dqdci;              /* dwdot[CH3O]/d[H2O] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[5];
    J[565] -= dqdci;              /* dwdot[OH]/d[CH2O] */
    J[566] += dqdci;              /* dwdot[H2O]/d[CH2O] */
    J[578] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[580] -= dqdci;              /* dwdot[CH3O]/d[CH2O] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[4];
    J[631] -= dqdci;              /* dwdot[OH]/d[CH3O] */
    J[632] += dqdci;              /* dwdot[H2O]/d[CH3O] */
    J[644] += dqdci;              /* dwdot[CH2O]/d[CH3O] */
    J[646] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[1060] -= dqdT;              /* dwdot[OH]/dT */
    J[1061] += dqdT;              /* dwdot[H2O]/dT */
    J[1073] += dqdT;              /* dwdot[CH2O]/dT */
    J[1075] -= dqdT;              /* dwdot[CH3O]/dT */

    /*reaction 110: OH + CH3OH <=> CH2OH + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[20];
    k_f = prefactor_units[109] * fwd_A[109]
                * exp(fwd_beta[109] * tc[0] - activation_units[109] * fwd_Ea[109] * invT);
    dlnkfdT = fwd_beta[109] * invT + activation_units[109] * fwd_Ea[109] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[18];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[18] + g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[20]) + (h_RT[5] + h_RT[18]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[18] += q; /* CH2OH */
    wdot[20] -= q; /* CH3OH */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[20];
    J[136] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[137] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[150] += dqdci;              /* dwdot[CH2OH]/d[OH] */
    J[152] -= dqdci;              /* dwdot[CH3OH]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[18];
    J[169] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[170] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[183] += dqdci;              /* dwdot[CH2OH]/d[H2O] */
    J[185] -= dqdci;              /* dwdot[CH3OH]/d[H2O] */
    /* d()/d[CH2OH] */
    dqdci =  - k_r*sc[5];
    J[598] -= dqdci;              /* dwdot[OH]/d[CH2OH] */
    J[599] += dqdci;              /* dwdot[H2O]/d[CH2OH] */
    J[612] += dqdci;              /* dwdot[CH2OH]/d[CH2OH] */
    J[614] -= dqdci;              /* dwdot[CH3OH]/d[CH2OH] */
    /* d()/d[CH3OH] */
    dqdci =  + k_f*sc[4];
    J[664] -= dqdci;              /* dwdot[OH]/d[CH3OH] */
    J[665] += dqdci;              /* dwdot[H2O]/d[CH3OH] */
    J[678] += dqdci;              /* dwdot[CH2OH]/d[CH3OH] */
    J[680] -= dqdci;              /* dwdot[CH3OH]/d[CH3OH] */
    /* d()/dT */
    J[1060] -= dqdT;              /* dwdot[OH]/dT */
    J[1061] += dqdT;              /* dwdot[H2O]/dT */
    J[1074] += dqdT;              /* dwdot[CH2OH]/dT */
    J[1076] -= dqdT;              /* dwdot[CH3OH]/dT */

    /*reaction 111: OH + CH3OH <=> CH3O + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[20];
    k_f = prefactor_units[110] * fwd_A[110]
                * exp(fwd_beta[110] * tc[0] - activation_units[110] * fwd_Ea[110] * invT);
    dlnkfdT = fwd_beta[110] * invT + activation_units[110] * fwd_Ea[110] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[19];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[19] + g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[20]) + (h_RT[5] + h_RT[19]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[19] += q; /* CH3O */
    wdot[20] -= q; /* CH3OH */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[20];
    J[136] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[137] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[151] += dqdci;              /* dwdot[CH3O]/d[OH] */
    J[152] -= dqdci;              /* dwdot[CH3OH]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[19];
    J[169] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[170] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[184] += dqdci;              /* dwdot[CH3O]/d[H2O] */
    J[185] -= dqdci;              /* dwdot[CH3OH]/d[H2O] */
    /* d()/d[CH3O] */
    dqdci =  - k_r*sc[5];
    J[631] -= dqdci;              /* dwdot[OH]/d[CH3O] */
    J[632] += dqdci;              /* dwdot[H2O]/d[CH3O] */
    J[646] += dqdci;              /* dwdot[CH3O]/d[CH3O] */
    J[647] -= dqdci;              /* dwdot[CH3OH]/d[CH3O] */
    /* d()/d[CH3OH] */
    dqdci =  + k_f*sc[4];
    J[664] -= dqdci;              /* dwdot[OH]/d[CH3OH] */
    J[665] += dqdci;              /* dwdot[H2O]/d[CH3OH] */
    J[679] += dqdci;              /* dwdot[CH3O]/d[CH3OH] */
    J[680] -= dqdci;              /* dwdot[CH3OH]/d[CH3OH] */
    /* d()/dT */
    J[1060] -= dqdT;              /* dwdot[OH]/dT */
    J[1061] += dqdT;              /* dwdot[H2O]/dT */
    J[1075] += dqdT;              /* dwdot[CH3O]/dT */
    J[1076] -= dqdT;              /* dwdot[CH3OH]/dT */

    /*reaction 112: OH + C2H <=> H + HCCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[21];
    k_f = prefactor_units[111] * fwd_A[111]
                * exp(fwd_beta[111] * tc[0] - activation_units[111] * fwd_Ea[111] * invT);
    dlnkfdT = fwd_beta[111] * invT + activation_units[111] * fwd_Ea[111] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[27];
    Kc = exp(-g_RT[1] + g_RT[4] + g_RT[21] - g_RT[27]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[21]) + (h_RT[1] + h_RT[27]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[4] -= q; /* OH */
    wdot[21] -= q; /* C2H */
    wdot[27] += q; /* HCCO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[27];
    J[34] += dqdci;               /* dwdot[H]/d[H] */
    J[37] -= dqdci;               /* dwdot[OH]/d[H] */
    J[54] -= dqdci;               /* dwdot[C2H]/d[H] */
    J[60] += dqdci;               /* dwdot[HCCO]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[21];
    J[133] += dqdci;              /* dwdot[H]/d[OH] */
    J[136] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[153] -= dqdci;              /* dwdot[C2H]/d[OH] */
    J[159] += dqdci;              /* dwdot[HCCO]/d[OH] */
    /* d()/d[C2H] */
    dqdci =  + k_f*sc[4];
    J[694] += dqdci;              /* dwdot[H]/d[C2H] */
    J[697] -= dqdci;              /* dwdot[OH]/d[C2H] */
    J[714] -= dqdci;              /* dwdot[C2H]/d[C2H] */
    J[720] += dqdci;              /* dwdot[HCCO]/d[C2H] */
    /* d()/d[HCCO] */
    dqdci =  - k_r*sc[1];
    J[892] += dqdci;              /* dwdot[H]/d[HCCO] */
    J[895] -= dqdci;              /* dwdot[OH]/d[HCCO] */
    J[912] -= dqdci;              /* dwdot[C2H]/d[HCCO] */
    J[918] += dqdci;              /* dwdot[HCCO]/d[HCCO] */
    /* d()/dT */
    J[1057] += dqdT;              /* dwdot[H]/dT */
    J[1060] -= dqdT;              /* dwdot[OH]/dT */
    J[1077] -= dqdT;              /* dwdot[C2H]/dT */
    J[1083] += dqdT;              /* dwdot[HCCO]/dT */

    /*reaction 113: OH + C2H2 <=> H + CH2CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[22];
    k_f = prefactor_units[112] * fwd_A[112]
                * exp(fwd_beta[112] * tc[0] - activation_units[112] * fwd_Ea[112] * invT);
    dlnkfdT = fwd_beta[112] * invT + activation_units[112] * fwd_Ea[112] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[28];
    Kc = exp(-g_RT[1] + g_RT[4] + g_RT[22] - g_RT[28]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[22]) + (h_RT[1] + h_RT[28]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[4] -= q; /* OH */
    wdot[22] -= q; /* C2H2 */
    wdot[28] += q; /* CH2CO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[28];
    J[34] += dqdci;               /* dwdot[H]/d[H] */
    J[37] -= dqdci;               /* dwdot[OH]/d[H] */
    J[55] -= dqdci;               /* dwdot[C2H2]/d[H] */
    J[61] += dqdci;               /* dwdot[CH2CO]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[22];
    J[133] += dqdci;              /* dwdot[H]/d[OH] */
    J[136] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[154] -= dqdci;              /* dwdot[C2H2]/d[OH] */
    J[160] += dqdci;              /* dwdot[CH2CO]/d[OH] */
    /* d()/d[C2H2] */
    dqdci =  + k_f*sc[4];
    J[727] += dqdci;              /* dwdot[H]/d[C2H2] */
    J[730] -= dqdci;              /* dwdot[OH]/d[C2H2] */
    J[748] -= dqdci;              /* dwdot[C2H2]/d[C2H2] */
    J[754] += dqdci;              /* dwdot[CH2CO]/d[C2H2] */
    /* d()/d[CH2CO] */
    dqdci =  - k_r*sc[1];
    J[925] += dqdci;              /* dwdot[H]/d[CH2CO] */
    J[928] -= dqdci;              /* dwdot[OH]/d[CH2CO] */
    J[946] -= dqdci;              /* dwdot[C2H2]/d[CH2CO] */
    J[952] += dqdci;              /* dwdot[CH2CO]/d[CH2CO] */
    /* d()/dT */
    J[1057] += dqdT;              /* dwdot[H]/dT */
    J[1060] -= dqdT;              /* dwdot[OH]/dT */
    J[1078] -= dqdT;              /* dwdot[C2H2]/dT */
    J[1084] += dqdT;              /* dwdot[CH2CO]/dT */

    /*reaction 114: OH + C2H2 <=> H + HCCOH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[22];
    k_f = prefactor_units[113] * fwd_A[113]
                * exp(fwd_beta[113] * tc[0] - activation_units[113] * fwd_Ea[113] * invT);
    dlnkfdT = fwd_beta[113] * invT + activation_units[113] * fwd_Ea[113] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[29];
    Kc = exp(-g_RT[1] + g_RT[4] + g_RT[22] - g_RT[29]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[22]) + (h_RT[1] + h_RT[29]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[4] -= q; /* OH */
    wdot[22] -= q; /* C2H2 */
    wdot[29] += q; /* HCCOH */
    /* d()/d[H] */
    dqdci =  - k_r*sc[29];
    J[34] += dqdci;               /* dwdot[H]/d[H] */
    J[37] -= dqdci;               /* dwdot[OH]/d[H] */
    J[55] -= dqdci;               /* dwdot[C2H2]/d[H] */
    J[62] += dqdci;               /* dwdot[HCCOH]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[22];
    J[133] += dqdci;              /* dwdot[H]/d[OH] */
    J[136] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[154] -= dqdci;              /* dwdot[C2H2]/d[OH] */
    J[161] += dqdci;              /* dwdot[HCCOH]/d[OH] */
    /* d()/d[C2H2] */
    dqdci =  + k_f*sc[4];
    J[727] += dqdci;              /* dwdot[H]/d[C2H2] */
    J[730] -= dqdci;              /* dwdot[OH]/d[C2H2] */
    J[748] -= dqdci;              /* dwdot[C2H2]/d[C2H2] */
    J[755] += dqdci;              /* dwdot[HCCOH]/d[C2H2] */
    /* d()/d[HCCOH] */
    dqdci =  - k_r*sc[1];
    J[958] += dqdci;              /* dwdot[H]/d[HCCOH] */
    J[961] -= dqdci;              /* dwdot[OH]/d[HCCOH] */
    J[979] -= dqdci;              /* dwdot[C2H2]/d[HCCOH] */
    J[986] += dqdci;              /* dwdot[HCCOH]/d[HCCOH] */
    /* d()/dT */
    J[1057] += dqdT;              /* dwdot[H]/dT */
    J[1060] -= dqdT;              /* dwdot[OH]/dT */
    J[1078] -= dqdT;              /* dwdot[C2H2]/dT */
    J[1085] += dqdT;              /* dwdot[HCCOH]/dT */

    /*reaction 115: OH + C2H2 <=> C2H + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[22];
    k_f = prefactor_units[114] * fwd_A[114]
                * exp(fwd_beta[114] * tc[0] - activation_units[114] * fwd_Ea[114] * invT);
    dlnkfdT = fwd_beta[114] * invT + activation_units[114] * fwd_Ea[114] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[21];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[21] + g_RT[22]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[22]) + (h_RT[5] + h_RT[21]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[21] += q; /* C2H */
    wdot[22] -= q; /* C2H2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[22];
    J[136] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[137] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[153] += dqdci;              /* dwdot[C2H]/d[OH] */
    J[154] -= dqdci;              /* dwdot[C2H2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[21];
    J[169] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[170] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[186] += dqdci;              /* dwdot[C2H]/d[H2O] */
    J[187] -= dqdci;              /* dwdot[C2H2]/d[H2O] */
    /* d()/d[C2H] */
    dqdci =  - k_r*sc[5];
    J[697] -= dqdci;              /* dwdot[OH]/d[C2H] */
    J[698] += dqdci;              /* dwdot[H2O]/d[C2H] */
    J[714] += dqdci;              /* dwdot[C2H]/d[C2H] */
    J[715] -= dqdci;              /* dwdot[C2H2]/d[C2H] */
    /* d()/d[C2H2] */
    dqdci =  + k_f*sc[4];
    J[730] -= dqdci;              /* dwdot[OH]/d[C2H2] */
    J[731] += dqdci;              /* dwdot[H2O]/d[C2H2] */
    J[747] += dqdci;              /* dwdot[C2H]/d[C2H2] */
    J[748] -= dqdci;              /* dwdot[C2H2]/d[C2H2] */
    /* d()/dT */
    J[1060] -= dqdT;              /* dwdot[OH]/dT */
    J[1061] += dqdT;              /* dwdot[H2O]/dT */
    J[1077] += dqdT;              /* dwdot[C2H]/dT */
    J[1078] -= dqdT;              /* dwdot[C2H2]/dT */

    /*reaction 116: OH + C2H2 <=> CH3 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[22];
    k_f = prefactor_units[115] * fwd_A[115]
                * exp(fwd_beta[115] * tc[0] - activation_units[115] * fwd_Ea[115] * invT);
    dlnkfdT = fwd_beta[115] * invT + activation_units[115] * fwd_Ea[115] * invT2;
    /* reverse */
    phi_r = sc[12]*sc[14];
    Kc = exp(g_RT[4] - g_RT[12] - g_RT[14] + g_RT[22]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[22]) + (h_RT[12] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[12] += q; /* CH3 */
    wdot[14] += q; /* CO */
    wdot[22] -= q; /* C2H2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[22];
    J[136] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[144] += dqdci;              /* dwdot[CH3]/d[OH] */
    J[146] += dqdci;              /* dwdot[CO]/d[OH] */
    J[154] -= dqdci;              /* dwdot[C2H2]/d[OH] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[14];
    J[400] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[408] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[410] += dqdci;              /* dwdot[CO]/d[CH3] */
    J[418] -= dqdci;              /* dwdot[C2H2]/d[CH3] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[12];
    J[466] -= dqdci;              /* dwdot[OH]/d[CO] */
    J[474] += dqdci;              /* dwdot[CH3]/d[CO] */
    J[476] += dqdci;              /* dwdot[CO]/d[CO] */
    J[484] -= dqdci;              /* dwdot[C2H2]/d[CO] */
    /* d()/d[C2H2] */
    dqdci =  + k_f*sc[4];
    J[730] -= dqdci;              /* dwdot[OH]/d[C2H2] */
    J[738] += dqdci;              /* dwdot[CH3]/d[C2H2] */
    J[740] += dqdci;              /* dwdot[CO]/d[C2H2] */
    J[748] -= dqdci;              /* dwdot[C2H2]/d[C2H2] */
    /* d()/dT */
    J[1060] -= dqdT;              /* dwdot[OH]/dT */
    J[1068] += dqdT;              /* dwdot[CH3]/dT */
    J[1070] += dqdT;              /* dwdot[CO]/dT */
    J[1078] -= dqdT;              /* dwdot[C2H2]/dT */

    /*reaction 117: OH + C2H3 <=> H2O + C2H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[23];
    k_f = prefactor_units[116] * fwd_A[116]
                * exp(fwd_beta[116] * tc[0] - activation_units[116] * fwd_Ea[116] * invT);
    dlnkfdT = fwd_beta[116] * invT + activation_units[116] * fwd_Ea[116] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[22];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[22] + g_RT[23]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[23]) + (h_RT[5] + h_RT[22]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[22] += q; /* C2H2 */
    wdot[23] -= q; /* C2H3 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[23];
    J[136] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[137] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[154] += dqdci;              /* dwdot[C2H2]/d[OH] */
    J[155] -= dqdci;              /* dwdot[C2H3]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[22];
    J[169] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[170] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[187] += dqdci;              /* dwdot[C2H2]/d[H2O] */
    J[188] -= dqdci;              /* dwdot[C2H3]/d[H2O] */
    /* d()/d[C2H2] */
    dqdci =  - k_r*sc[5];
    J[730] -= dqdci;              /* dwdot[OH]/d[C2H2] */
    J[731] += dqdci;              /* dwdot[H2O]/d[C2H2] */
    J[748] += dqdci;              /* dwdot[C2H2]/d[C2H2] */
    J[749] -= dqdci;              /* dwdot[C2H3]/d[C2H2] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[4];
    J[763] -= dqdci;              /* dwdot[OH]/d[C2H3] */
    J[764] += dqdci;              /* dwdot[H2O]/d[C2H3] */
    J[781] += dqdci;              /* dwdot[C2H2]/d[C2H3] */
    J[782] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    /* d()/dT */
    J[1060] -= dqdT;              /* dwdot[OH]/dT */
    J[1061] += dqdT;              /* dwdot[H2O]/dT */
    J[1078] += dqdT;              /* dwdot[C2H2]/dT */
    J[1079] -= dqdT;              /* dwdot[C2H3]/dT */

    /*reaction 118: OH + C2H4 <=> C2H3 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[24];
    k_f = prefactor_units[117] * fwd_A[117]
                * exp(fwd_beta[117] * tc[0] - activation_units[117] * fwd_Ea[117] * invT);
    dlnkfdT = fwd_beta[117] * invT + activation_units[117] * fwd_Ea[117] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[23];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[23] + g_RT[24]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[24]) + (h_RT[5] + h_RT[23]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[23] += q; /* C2H3 */
    wdot[24] -= q; /* C2H4 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[24];
    J[136] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[137] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[155] += dqdci;              /* dwdot[C2H3]/d[OH] */
    J[156] -= dqdci;              /* dwdot[C2H4]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[23];
    J[169] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[170] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[188] += dqdci;              /* dwdot[C2H3]/d[H2O] */
    J[189] -= dqdci;              /* dwdot[C2H4]/d[H2O] */
    /* d()/d[C2H3] */
    dqdci =  - k_r*sc[5];
    J[763] -= dqdci;              /* dwdot[OH]/d[C2H3] */
    J[764] += dqdci;              /* dwdot[H2O]/d[C2H3] */
    J[782] += dqdci;              /* dwdot[C2H3]/d[C2H3] */
    J[783] -= dqdci;              /* dwdot[C2H4]/d[C2H3] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[4];
    J[796] -= dqdci;              /* dwdot[OH]/d[C2H4] */
    J[797] += dqdci;              /* dwdot[H2O]/d[C2H4] */
    J[815] += dqdci;              /* dwdot[C2H3]/d[C2H4] */
    J[816] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    /* d()/dT */
    J[1060] -= dqdT;              /* dwdot[OH]/dT */
    J[1061] += dqdT;              /* dwdot[H2O]/dT */
    J[1079] += dqdT;              /* dwdot[C2H3]/dT */
    J[1080] -= dqdT;              /* dwdot[C2H4]/dT */

    /*reaction 119: OH + C2H6 <=> C2H5 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[26];
    k_f = prefactor_units[118] * fwd_A[118]
                * exp(fwd_beta[118] * tc[0] - activation_units[118] * fwd_Ea[118] * invT);
    dlnkfdT = fwd_beta[118] * invT + activation_units[118] * fwd_Ea[118] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[25];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[25] + g_RT[26]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[26]) + (h_RT[5] + h_RT[25]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[25] += q; /* C2H5 */
    wdot[26] -= q; /* C2H6 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[26];
    J[136] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[137] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[157] += dqdci;              /* dwdot[C2H5]/d[OH] */
    J[158] -= dqdci;              /* dwdot[C2H6]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[25];
    J[169] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[170] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[190] += dqdci;              /* dwdot[C2H5]/d[H2O] */
    J[191] -= dqdci;              /* dwdot[C2H6]/d[H2O] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[5];
    J[829] -= dqdci;              /* dwdot[OH]/d[C2H5] */
    J[830] += dqdci;              /* dwdot[H2O]/d[C2H5] */
    J[850] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[851] -= dqdci;              /* dwdot[C2H6]/d[C2H5] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[4];
    J[862] -= dqdci;              /* dwdot[OH]/d[C2H6] */
    J[863] += dqdci;              /* dwdot[H2O]/d[C2H6] */
    J[883] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[884] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[1060] -= dqdT;              /* dwdot[OH]/dT */
    J[1061] += dqdT;              /* dwdot[H2O]/dT */
    J[1081] += dqdT;              /* dwdot[C2H5]/dT */
    J[1082] -= dqdT;              /* dwdot[C2H6]/dT */

    /*reaction 120: OH + CH2CO <=> HCCO + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[28];
    k_f = prefactor_units[119] * fwd_A[119]
                * exp(fwd_beta[119] * tc[0] - activation_units[119] * fwd_Ea[119] * invT);
    dlnkfdT = fwd_beta[119] * invT + activation_units[119] * fwd_Ea[119] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[27];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[27] + g_RT[28]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[28]) + (h_RT[5] + h_RT[27]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[27] += q; /* HCCO */
    wdot[28] -= q; /* CH2CO */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[28];
    J[136] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[137] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[159] += dqdci;              /* dwdot[HCCO]/d[OH] */
    J[160] -= dqdci;              /* dwdot[CH2CO]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[27];
    J[169] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[170] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[192] += dqdci;              /* dwdot[HCCO]/d[H2O] */
    J[193] -= dqdci;              /* dwdot[CH2CO]/d[H2O] */
    /* d()/d[HCCO] */
    dqdci =  - k_r*sc[5];
    J[895] -= dqdci;              /* dwdot[OH]/d[HCCO] */
    J[896] += dqdci;              /* dwdot[H2O]/d[HCCO] */
    J[918] += dqdci;              /* dwdot[HCCO]/d[HCCO] */
    J[919] -= dqdci;              /* dwdot[CH2CO]/d[HCCO] */
    /* d()/d[CH2CO] */
    dqdci =  + k_f*sc[4];
    J[928] -= dqdci;              /* dwdot[OH]/d[CH2CO] */
    J[929] += dqdci;              /* dwdot[H2O]/d[CH2CO] */
    J[951] += dqdci;              /* dwdot[HCCO]/d[CH2CO] */
    J[952] -= dqdci;              /* dwdot[CH2CO]/d[CH2CO] */
    /* d()/dT */
    J[1060] -= dqdT;              /* dwdot[OH]/dT */
    J[1061] += dqdT;              /* dwdot[H2O]/dT */
    J[1083] += dqdT;              /* dwdot[HCCO]/dT */
    J[1084] -= dqdT;              /* dwdot[CH2CO]/dT */

    /*reaction 121: 2.000000 HO2 <=> O2 + H2O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[6], 2.000000);
    k_f = prefactor_units[120] * fwd_A[120]
                * exp(fwd_beta[120] * tc[0] - activation_units[120] * fwd_Ea[120] * invT);
    dlnkfdT = fwd_beta[120] * invT + activation_units[120] * fwd_Ea[120] * invT2;
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
    J[102] += dqdci;              /* dwdot[O2]/d[O2] */
    J[105] += -2 * dqdci;         /* dwdot[HO2]/d[O2] */
    J[106] += dqdci;              /* dwdot[H2O2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*2.000000*sc[6];
    J[201] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[204] += -2 * dqdci;         /* dwdot[HO2]/d[HO2] */
    J[205] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[3];
    J[234] += dqdci;              /* dwdot[O2]/d[H2O2] */
    J[237] += -2 * dqdci;         /* dwdot[HO2]/d[H2O2] */
    J[238] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[1059] += dqdT;              /* dwdot[O2]/dT */
    J[1062] += -2 * dqdT;         /* dwdot[HO2]/dT */
    J[1063] += dqdT;              /* dwdot[H2O2]/dT */

    /*reaction 122: 2.000000 HO2 <=> O2 + H2O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[6], 2.000000);
    k_f = prefactor_units[121] * fwd_A[121]
                * exp(fwd_beta[121] * tc[0] - activation_units[121] * fwd_Ea[121] * invT);
    dlnkfdT = fwd_beta[121] * invT + activation_units[121] * fwd_Ea[121] * invT2;
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
    J[102] += dqdci;              /* dwdot[O2]/d[O2] */
    J[105] += -2 * dqdci;         /* dwdot[HO2]/d[O2] */
    J[106] += dqdci;              /* dwdot[H2O2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*2.000000*sc[6];
    J[201] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[204] += -2 * dqdci;         /* dwdot[HO2]/d[HO2] */
    J[205] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[3];
    J[234] += dqdci;              /* dwdot[O2]/d[H2O2] */
    J[237] += -2 * dqdci;         /* dwdot[HO2]/d[H2O2] */
    J[238] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[1059] += dqdT;              /* dwdot[O2]/dT */
    J[1062] += -2 * dqdT;         /* dwdot[HO2]/dT */
    J[1063] += dqdT;              /* dwdot[H2O2]/dT */

    /*reaction 123: HO2 + CH2 <=> OH + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[10];
    k_f = prefactor_units[122] * fwd_A[122]
                * exp(fwd_beta[122] * tc[0] - activation_units[122] * fwd_Ea[122] * invT);
    dlnkfdT = fwd_beta[122] * invT + activation_units[122] * fwd_Ea[122] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[17];
    Kc = exp(-g_RT[4] + g_RT[6] + g_RT[10] - g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[10]) + (h_RT[4] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] += q; /* OH */
    wdot[6] -= q; /* HO2 */
    wdot[10] -= q; /* CH2 */
    wdot[17] += q; /* CH2O */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[17];
    J[136] += dqdci;              /* dwdot[OH]/d[OH] */
    J[138] -= dqdci;              /* dwdot[HO2]/d[OH] */
    J[142] -= dqdci;              /* dwdot[CH2]/d[OH] */
    J[149] += dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[10];
    J[202] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[204] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[208] -= dqdci;              /* dwdot[CH2]/d[HO2] */
    J[215] += dqdci;              /* dwdot[CH2O]/d[HO2] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[6];
    J[334] += dqdci;              /* dwdot[OH]/d[CH2] */
    J[336] -= dqdci;              /* dwdot[HO2]/d[CH2] */
    J[340] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[347] += dqdci;              /* dwdot[CH2O]/d[CH2] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[4];
    J[565] += dqdci;              /* dwdot[OH]/d[CH2O] */
    J[567] -= dqdci;              /* dwdot[HO2]/d[CH2O] */
    J[571] -= dqdci;              /* dwdot[CH2]/d[CH2O] */
    J[578] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[1060] += dqdT;              /* dwdot[OH]/dT */
    J[1062] -= dqdT;              /* dwdot[HO2]/dT */
    J[1066] -= dqdT;              /* dwdot[CH2]/dT */
    J[1073] += dqdT;              /* dwdot[CH2O]/dT */

    /*reaction 124: HO2 + CH3 <=> O2 + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[12];
    k_f = prefactor_units[123] * fwd_A[123]
                * exp(fwd_beta[123] * tc[0] - activation_units[123] * fwd_Ea[123] * invT);
    dlnkfdT = fwd_beta[123] * invT + activation_units[123] * fwd_Ea[123] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[13];
    Kc = exp(-g_RT[3] + g_RT[6] + g_RT[12] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[12]) + (h_RT[3] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += q; /* O2 */
    wdot[6] -= q; /* HO2 */
    wdot[12] -= q; /* CH3 */
    wdot[13] += q; /* CH4 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[13];
    J[102] += dqdci;              /* dwdot[O2]/d[O2] */
    J[105] -= dqdci;              /* dwdot[HO2]/d[O2] */
    J[111] -= dqdci;              /* dwdot[CH3]/d[O2] */
    J[112] += dqdci;              /* dwdot[CH4]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[12];
    J[201] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[204] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[210] -= dqdci;              /* dwdot[CH3]/d[HO2] */
    J[211] += dqdci;              /* dwdot[CH4]/d[HO2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[6];
    J[399] += dqdci;              /* dwdot[O2]/d[CH3] */
    J[402] -= dqdci;              /* dwdot[HO2]/d[CH3] */
    J[408] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[409] += dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[3];
    J[432] += dqdci;              /* dwdot[O2]/d[CH4] */
    J[435] -= dqdci;              /* dwdot[HO2]/d[CH4] */
    J[441] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[442] += dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[1059] += dqdT;              /* dwdot[O2]/dT */
    J[1062] -= dqdT;              /* dwdot[HO2]/dT */
    J[1068] -= dqdT;              /* dwdot[CH3]/dT */
    J[1069] += dqdT;              /* dwdot[CH4]/dT */

    /*reaction 125: HO2 + CH3 <=> OH + CH3O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[12];
    k_f = prefactor_units[124] * fwd_A[124]
                * exp(fwd_beta[124] * tc[0] - activation_units[124] * fwd_Ea[124] * invT);
    dlnkfdT = fwd_beta[124] * invT + activation_units[124] * fwd_Ea[124] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[19];
    Kc = exp(-g_RT[4] + g_RT[6] + g_RT[12] - g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[12]) + (h_RT[4] + h_RT[19]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] += q; /* OH */
    wdot[6] -= q; /* HO2 */
    wdot[12] -= q; /* CH3 */
    wdot[19] += q; /* CH3O */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[19];
    J[136] += dqdci;              /* dwdot[OH]/d[OH] */
    J[138] -= dqdci;              /* dwdot[HO2]/d[OH] */
    J[144] -= dqdci;              /* dwdot[CH3]/d[OH] */
    J[151] += dqdci;              /* dwdot[CH3O]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[12];
    J[202] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[204] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[210] -= dqdci;              /* dwdot[CH3]/d[HO2] */
    J[217] += dqdci;              /* dwdot[CH3O]/d[HO2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[6];
    J[400] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[402] -= dqdci;              /* dwdot[HO2]/d[CH3] */
    J[408] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[415] += dqdci;              /* dwdot[CH3O]/d[CH3] */
    /* d()/d[CH3O] */
    dqdci =  - k_r*sc[4];
    J[631] += dqdci;              /* dwdot[OH]/d[CH3O] */
    J[633] -= dqdci;              /* dwdot[HO2]/d[CH3O] */
    J[639] -= dqdci;              /* dwdot[CH3]/d[CH3O] */
    J[646] += dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[1060] += dqdT;              /* dwdot[OH]/dT */
    J[1062] -= dqdT;              /* dwdot[HO2]/dT */
    J[1068] -= dqdT;              /* dwdot[CH3]/dT */
    J[1075] += dqdT;              /* dwdot[CH3O]/dT */

    /*reaction 126: HO2 + CO <=> OH + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[14];
    k_f = prefactor_units[125] * fwd_A[125]
                * exp(fwd_beta[125] * tc[0] - activation_units[125] * fwd_Ea[125] * invT);
    dlnkfdT = fwd_beta[125] * invT + activation_units[125] * fwd_Ea[125] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[15];
    Kc = exp(-g_RT[4] + g_RT[6] + g_RT[14] - g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[14]) + (h_RT[4] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] += q; /* OH */
    wdot[6] -= q; /* HO2 */
    wdot[14] -= q; /* CO */
    wdot[15] += q; /* CO2 */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[15];
    J[136] += dqdci;              /* dwdot[OH]/d[OH] */
    J[138] -= dqdci;              /* dwdot[HO2]/d[OH] */
    J[146] -= dqdci;              /* dwdot[CO]/d[OH] */
    J[147] += dqdci;              /* dwdot[CO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[14];
    J[202] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[204] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[212] -= dqdci;              /* dwdot[CO]/d[HO2] */
    J[213] += dqdci;              /* dwdot[CO2]/d[HO2] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[6];
    J[466] += dqdci;              /* dwdot[OH]/d[CO] */
    J[468] -= dqdci;              /* dwdot[HO2]/d[CO] */
    J[476] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[477] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[4];
    J[499] += dqdci;              /* dwdot[OH]/d[CO2] */
    J[501] -= dqdci;              /* dwdot[HO2]/d[CO2] */
    J[509] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[510] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[1060] += dqdT;              /* dwdot[OH]/dT */
    J[1062] -= dqdT;              /* dwdot[HO2]/dT */
    J[1070] -= dqdT;              /* dwdot[CO]/dT */
    J[1071] += dqdT;              /* dwdot[CO2]/dT */

    /*reaction 127: HO2 + CH2O <=> HCO + H2O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[17];
    k_f = prefactor_units[126] * fwd_A[126]
                * exp(fwd_beta[126] * tc[0] - activation_units[126] * fwd_Ea[126] * invT);
    dlnkfdT = fwd_beta[126] * invT + activation_units[126] * fwd_Ea[126] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[16];
    Kc = exp(g_RT[6] - g_RT[7] - g_RT[16] + g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[17]) + (h_RT[7] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[6] -= q; /* HO2 */
    wdot[7] += q; /* H2O2 */
    wdot[16] += q; /* HCO */
    wdot[17] -= q; /* CH2O */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[17];
    J[204] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[205] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    J[214] += dqdci;              /* dwdot[HCO]/d[HO2] */
    J[215] -= dqdci;              /* dwdot[CH2O]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[16];
    J[237] -= dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[238] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
    J[247] += dqdci;              /* dwdot[HCO]/d[H2O2] */
    J[248] -= dqdci;              /* dwdot[CH2O]/d[H2O2] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[7];
    J[534] -= dqdci;              /* dwdot[HO2]/d[HCO] */
    J[535] += dqdci;              /* dwdot[H2O2]/d[HCO] */
    J[544] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[545] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[6];
    J[567] -= dqdci;              /* dwdot[HO2]/d[CH2O] */
    J[568] += dqdci;              /* dwdot[H2O2]/d[CH2O] */
    J[577] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[578] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[1062] -= dqdT;              /* dwdot[HO2]/dT */
    J[1063] += dqdT;              /* dwdot[H2O2]/dT */
    J[1072] += dqdT;              /* dwdot[HCO]/dT */
    J[1073] -= dqdT;              /* dwdot[CH2O]/dT */

    /*reaction 128: C + O2 <=> O + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[8];
    k_f = prefactor_units[127] * fwd_A[127]
                * exp(fwd_beta[127] * tc[0] - activation_units[127] * fwd_Ea[127] * invT);
    dlnkfdT = fwd_beta[127] * invT + activation_units[127] * fwd_Ea[127] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[14];
    Kc = exp(-g_RT[2] + g_RT[3] + g_RT[8] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[8]) + (h_RT[2] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* O */
    wdot[3] -= q; /* O2 */
    wdot[8] -= q; /* C */
    wdot[14] += q; /* CO */
    /* d()/d[O] */
    dqdci =  - k_r*sc[14];
    J[68] += dqdci;               /* dwdot[O]/d[O] */
    J[69] -= dqdci;               /* dwdot[O2]/d[O] */
    J[74] -= dqdci;               /* dwdot[C]/d[O] */
    J[80] += dqdci;               /* dwdot[CO]/d[O] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[8];
    J[101] += dqdci;              /* dwdot[O]/d[O2] */
    J[102] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[107] -= dqdci;              /* dwdot[C]/d[O2] */
    J[113] += dqdci;              /* dwdot[CO]/d[O2] */
    /* d()/d[C] */
    dqdci =  + k_f*sc[3];
    J[266] += dqdci;              /* dwdot[O]/d[C] */
    J[267] -= dqdci;              /* dwdot[O2]/d[C] */
    J[272] -= dqdci;              /* dwdot[C]/d[C] */
    J[278] += dqdci;              /* dwdot[CO]/d[C] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[2];
    J[464] += dqdci;              /* dwdot[O]/d[CO] */
    J[465] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[470] -= dqdci;              /* dwdot[C]/d[CO] */
    J[476] += dqdci;              /* dwdot[CO]/d[CO] */
    /* d()/dT */
    J[1058] += dqdT;              /* dwdot[O]/dT */
    J[1059] -= dqdT;              /* dwdot[O2]/dT */
    J[1064] -= dqdT;              /* dwdot[C]/dT */
    J[1070] += dqdT;              /* dwdot[CO]/dT */

    /*reaction 129: C + CH2 <=> H + C2H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[10];
    k_f = prefactor_units[128] * fwd_A[128]
                * exp(fwd_beta[128] * tc[0] - activation_units[128] * fwd_Ea[128] * invT);
    dlnkfdT = fwd_beta[128] * invT + activation_units[128] * fwd_Ea[128] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[21];
    Kc = exp(-g_RT[1] + g_RT[8] + g_RT[10] - g_RT[21]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[10]) + (h_RT[1] + h_RT[21]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[8] -= q; /* C */
    wdot[10] -= q; /* CH2 */
    wdot[21] += q; /* C2H */
    /* d()/d[H] */
    dqdci =  - k_r*sc[21];
    J[34] += dqdci;               /* dwdot[H]/d[H] */
    J[41] -= dqdci;               /* dwdot[C]/d[H] */
    J[43] -= dqdci;               /* dwdot[CH2]/d[H] */
    J[54] += dqdci;               /* dwdot[C2H]/d[H] */
    /* d()/d[C] */
    dqdci =  + k_f*sc[10];
    J[265] += dqdci;              /* dwdot[H]/d[C] */
    J[272] -= dqdci;              /* dwdot[C]/d[C] */
    J[274] -= dqdci;              /* dwdot[CH2]/d[C] */
    J[285] += dqdci;              /* dwdot[C2H]/d[C] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[8];
    J[331] += dqdci;              /* dwdot[H]/d[CH2] */
    J[338] -= dqdci;              /* dwdot[C]/d[CH2] */
    J[340] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[351] += dqdci;              /* dwdot[C2H]/d[CH2] */
    /* d()/d[C2H] */
    dqdci =  - k_r*sc[1];
    J[694] += dqdci;              /* dwdot[H]/d[C2H] */
    J[701] -= dqdci;              /* dwdot[C]/d[C2H] */
    J[703] -= dqdci;              /* dwdot[CH2]/d[C2H] */
    J[714] += dqdci;              /* dwdot[C2H]/d[C2H] */
    /* d()/dT */
    J[1057] += dqdT;              /* dwdot[H]/dT */
    J[1064] -= dqdT;              /* dwdot[C]/dT */
    J[1066] -= dqdT;              /* dwdot[CH2]/dT */
    J[1077] += dqdT;              /* dwdot[C2H]/dT */

    /*reaction 130: C + CH3 <=> H + C2H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[12];
    k_f = prefactor_units[129] * fwd_A[129]
                * exp(fwd_beta[129] * tc[0] - activation_units[129] * fwd_Ea[129] * invT);
    dlnkfdT = fwd_beta[129] * invT + activation_units[129] * fwd_Ea[129] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[22];
    Kc = exp(-g_RT[1] + g_RT[8] + g_RT[12] - g_RT[22]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[12]) + (h_RT[1] + h_RT[22]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[8] -= q; /* C */
    wdot[12] -= q; /* CH3 */
    wdot[22] += q; /* C2H2 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[22];
    J[34] += dqdci;               /* dwdot[H]/d[H] */
    J[41] -= dqdci;               /* dwdot[C]/d[H] */
    J[45] -= dqdci;               /* dwdot[CH3]/d[H] */
    J[55] += dqdci;               /* dwdot[C2H2]/d[H] */
    /* d()/d[C] */
    dqdci =  + k_f*sc[12];
    J[265] += dqdci;              /* dwdot[H]/d[C] */
    J[272] -= dqdci;              /* dwdot[C]/d[C] */
    J[276] -= dqdci;              /* dwdot[CH3]/d[C] */
    J[286] += dqdci;              /* dwdot[C2H2]/d[C] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[8];
    J[397] += dqdci;              /* dwdot[H]/d[CH3] */
    J[404] -= dqdci;              /* dwdot[C]/d[CH3] */
    J[408] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[418] += dqdci;              /* dwdot[C2H2]/d[CH3] */
    /* d()/d[C2H2] */
    dqdci =  - k_r*sc[1];
    J[727] += dqdci;              /* dwdot[H]/d[C2H2] */
    J[734] -= dqdci;              /* dwdot[C]/d[C2H2] */
    J[738] -= dqdci;              /* dwdot[CH3]/d[C2H2] */
    J[748] += dqdci;              /* dwdot[C2H2]/d[C2H2] */
    /* d()/dT */
    J[1057] += dqdT;              /* dwdot[H]/dT */
    J[1064] -= dqdT;              /* dwdot[C]/dT */
    J[1068] -= dqdT;              /* dwdot[CH3]/dT */
    J[1078] += dqdT;              /* dwdot[C2H2]/dT */

    /*reaction 131: CH + O2 <=> O + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[9];
    k_f = prefactor_units[130] * fwd_A[130]
                * exp(fwd_beta[130] * tc[0] - activation_units[130] * fwd_Ea[130] * invT);
    dlnkfdT = fwd_beta[130] * invT + activation_units[130] * fwd_Ea[130] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[16];
    Kc = exp(-g_RT[2] + g_RT[3] + g_RT[9] - g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[9]) + (h_RT[2] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* O */
    wdot[3] -= q; /* O2 */
    wdot[9] -= q; /* CH */
    wdot[16] += q; /* HCO */
    /* d()/d[O] */
    dqdci =  - k_r*sc[16];
    J[68] += dqdci;               /* dwdot[O]/d[O] */
    J[69] -= dqdci;               /* dwdot[O2]/d[O] */
    J[75] -= dqdci;               /* dwdot[CH]/d[O] */
    J[82] += dqdci;               /* dwdot[HCO]/d[O] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[9];
    J[101] += dqdci;              /* dwdot[O]/d[O2] */
    J[102] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[108] -= dqdci;              /* dwdot[CH]/d[O2] */
    J[115] += dqdci;              /* dwdot[HCO]/d[O2] */
    /* d()/d[CH] */
    dqdci =  + k_f*sc[3];
    J[299] += dqdci;              /* dwdot[O]/d[CH] */
    J[300] -= dqdci;              /* dwdot[O2]/d[CH] */
    J[306] -= dqdci;              /* dwdot[CH]/d[CH] */
    J[313] += dqdci;              /* dwdot[HCO]/d[CH] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[2];
    J[530] += dqdci;              /* dwdot[O]/d[HCO] */
    J[531] -= dqdci;              /* dwdot[O2]/d[HCO] */
    J[537] -= dqdci;              /* dwdot[CH]/d[HCO] */
    J[544] += dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[1058] += dqdT;              /* dwdot[O]/dT */
    J[1059] -= dqdT;              /* dwdot[O2]/dT */
    J[1065] -= dqdT;              /* dwdot[CH]/dT */
    J[1072] += dqdT;              /* dwdot[HCO]/dT */

    /*reaction 132: CH + H2 <=> H + CH2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[9];
    k_f = prefactor_units[131] * fwd_A[131]
                * exp(fwd_beta[131] * tc[0] - activation_units[131] * fwd_Ea[131] * invT);
    dlnkfdT = fwd_beta[131] * invT + activation_units[131] * fwd_Ea[131] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[10];
    Kc = exp(g_RT[0] - g_RT[1] + g_RT[9] - g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[9]) + (h_RT[1] + h_RT[10]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[1] += q; /* H */
    wdot[9] -= q; /* CH */
    wdot[10] += q; /* CH2 */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[9];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += dqdci;                /* dwdot[H]/d[H2] */
    J[9] -= dqdci;                /* dwdot[CH]/d[H2] */
    J[10] += dqdci;               /* dwdot[CH2]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[10];
    J[33] -= dqdci;               /* dwdot[H2]/d[H] */
    J[34] += dqdci;               /* dwdot[H]/d[H] */
    J[42] -= dqdci;               /* dwdot[CH]/d[H] */
    J[43] += dqdci;               /* dwdot[CH2]/d[H] */
    /* d()/d[CH] */
    dqdci =  + k_f*sc[0];
    J[297] -= dqdci;              /* dwdot[H2]/d[CH] */
    J[298] += dqdci;              /* dwdot[H]/d[CH] */
    J[306] -= dqdci;              /* dwdot[CH]/d[CH] */
    J[307] += dqdci;              /* dwdot[CH2]/d[CH] */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[1];
    J[330] -= dqdci;              /* dwdot[H2]/d[CH2] */
    J[331] += dqdci;              /* dwdot[H]/d[CH2] */
    J[339] -= dqdci;              /* dwdot[CH]/d[CH2] */
    J[340] += dqdci;              /* dwdot[CH2]/d[CH2] */
    /* d()/dT */
    J[1056] -= dqdT;              /* dwdot[H2]/dT */
    J[1057] += dqdT;              /* dwdot[H]/dT */
    J[1065] -= dqdT;              /* dwdot[CH]/dT */
    J[1066] += dqdT;              /* dwdot[CH2]/dT */

    /*reaction 133: CH + H2O <=> H + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[9];
    k_f = prefactor_units[132] * fwd_A[132]
                * exp(fwd_beta[132] * tc[0] - activation_units[132] * fwd_Ea[132] * invT);
    dlnkfdT = fwd_beta[132] * invT + activation_units[132] * fwd_Ea[132] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[17];
    Kc = exp(-g_RT[1] + g_RT[5] + g_RT[9] - g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[9]) + (h_RT[1] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[5] -= q; /* H2O */
    wdot[9] -= q; /* CH */
    wdot[17] += q; /* CH2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[17];
    J[34] += dqdci;               /* dwdot[H]/d[H] */
    J[38] -= dqdci;               /* dwdot[H2O]/d[H] */
    J[42] -= dqdci;               /* dwdot[CH]/d[H] */
    J[50] += dqdci;               /* dwdot[CH2O]/d[H] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[9];
    J[166] += dqdci;              /* dwdot[H]/d[H2O] */
    J[170] -= dqdci;              /* dwdot[H2O]/d[H2O] */
    J[174] -= dqdci;              /* dwdot[CH]/d[H2O] */
    J[182] += dqdci;              /* dwdot[CH2O]/d[H2O] */
    /* d()/d[CH] */
    dqdci =  + k_f*sc[5];
    J[298] += dqdci;              /* dwdot[H]/d[CH] */
    J[302] -= dqdci;              /* dwdot[H2O]/d[CH] */
    J[306] -= dqdci;              /* dwdot[CH]/d[CH] */
    J[314] += dqdci;              /* dwdot[CH2O]/d[CH] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[1];
    J[562] += dqdci;              /* dwdot[H]/d[CH2O] */
    J[566] -= dqdci;              /* dwdot[H2O]/d[CH2O] */
    J[570] -= dqdci;              /* dwdot[CH]/d[CH2O] */
    J[578] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[1057] += dqdT;              /* dwdot[H]/dT */
    J[1061] -= dqdT;              /* dwdot[H2O]/dT */
    J[1065] -= dqdT;              /* dwdot[CH]/dT */
    J[1073] += dqdT;              /* dwdot[CH2O]/dT */

    /*reaction 134: CH + CH2 <=> H + C2H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[10];
    k_f = prefactor_units[133] * fwd_A[133]
                * exp(fwd_beta[133] * tc[0] - activation_units[133] * fwd_Ea[133] * invT);
    dlnkfdT = fwd_beta[133] * invT + activation_units[133] * fwd_Ea[133] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[22];
    Kc = exp(-g_RT[1] + g_RT[9] + g_RT[10] - g_RT[22]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[9] + h_RT[10]) + (h_RT[1] + h_RT[22]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[9] -= q; /* CH */
    wdot[10] -= q; /* CH2 */
    wdot[22] += q; /* C2H2 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[22];
    J[34] += dqdci;               /* dwdot[H]/d[H] */
    J[42] -= dqdci;               /* dwdot[CH]/d[H] */
    J[43] -= dqdci;               /* dwdot[CH2]/d[H] */
    J[55] += dqdci;               /* dwdot[C2H2]/d[H] */
    /* d()/d[CH] */
    dqdci =  + k_f*sc[10];
    J[298] += dqdci;              /* dwdot[H]/d[CH] */
    J[306] -= dqdci;              /* dwdot[CH]/d[CH] */
    J[307] -= dqdci;              /* dwdot[CH2]/d[CH] */
    J[319] += dqdci;              /* dwdot[C2H2]/d[CH] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[9];
    J[331] += dqdci;              /* dwdot[H]/d[CH2] */
    J[339] -= dqdci;              /* dwdot[CH]/d[CH2] */
    J[340] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[352] += dqdci;              /* dwdot[C2H2]/d[CH2] */
    /* d()/d[C2H2] */
    dqdci =  - k_r*sc[1];
    J[727] += dqdci;              /* dwdot[H]/d[C2H2] */
    J[735] -= dqdci;              /* dwdot[CH]/d[C2H2] */
    J[736] -= dqdci;              /* dwdot[CH2]/d[C2H2] */
    J[748] += dqdci;              /* dwdot[C2H2]/d[C2H2] */
    /* d()/dT */
    J[1057] += dqdT;              /* dwdot[H]/dT */
    J[1065] -= dqdT;              /* dwdot[CH]/dT */
    J[1066] -= dqdT;              /* dwdot[CH2]/dT */
    J[1078] += dqdT;              /* dwdot[C2H2]/dT */

    /*reaction 135: CH + CH3 <=> H + C2H3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[12];
    k_f = prefactor_units[134] * fwd_A[134]
                * exp(fwd_beta[134] * tc[0] - activation_units[134] * fwd_Ea[134] * invT);
    dlnkfdT = fwd_beta[134] * invT + activation_units[134] * fwd_Ea[134] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[23];
    Kc = exp(-g_RT[1] + g_RT[9] + g_RT[12] - g_RT[23]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[9] + h_RT[12]) + (h_RT[1] + h_RT[23]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[9] -= q; /* CH */
    wdot[12] -= q; /* CH3 */
    wdot[23] += q; /* C2H3 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[23];
    J[34] += dqdci;               /* dwdot[H]/d[H] */
    J[42] -= dqdci;               /* dwdot[CH]/d[H] */
    J[45] -= dqdci;               /* dwdot[CH3]/d[H] */
    J[56] += dqdci;               /* dwdot[C2H3]/d[H] */
    /* d()/d[CH] */
    dqdci =  + k_f*sc[12];
    J[298] += dqdci;              /* dwdot[H]/d[CH] */
    J[306] -= dqdci;              /* dwdot[CH]/d[CH] */
    J[309] -= dqdci;              /* dwdot[CH3]/d[CH] */
    J[320] += dqdci;              /* dwdot[C2H3]/d[CH] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[9];
    J[397] += dqdci;              /* dwdot[H]/d[CH3] */
    J[405] -= dqdci;              /* dwdot[CH]/d[CH3] */
    J[408] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[419] += dqdci;              /* dwdot[C2H3]/d[CH3] */
    /* d()/d[C2H3] */
    dqdci =  - k_r*sc[1];
    J[760] += dqdci;              /* dwdot[H]/d[C2H3] */
    J[768] -= dqdci;              /* dwdot[CH]/d[C2H3] */
    J[771] -= dqdci;              /* dwdot[CH3]/d[C2H3] */
    J[782] += dqdci;              /* dwdot[C2H3]/d[C2H3] */
    /* d()/dT */
    J[1057] += dqdT;              /* dwdot[H]/dT */
    J[1065] -= dqdT;              /* dwdot[CH]/dT */
    J[1068] -= dqdT;              /* dwdot[CH3]/dT */
    J[1079] += dqdT;              /* dwdot[C2H3]/dT */

    /*reaction 136: CH + CH4 <=> H + C2H4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[13];
    k_f = prefactor_units[135] * fwd_A[135]
                * exp(fwd_beta[135] * tc[0] - activation_units[135] * fwd_Ea[135] * invT);
    dlnkfdT = fwd_beta[135] * invT + activation_units[135] * fwd_Ea[135] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[24];
    Kc = exp(-g_RT[1] + g_RT[9] + g_RT[13] - g_RT[24]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[9] + h_RT[13]) + (h_RT[1] + h_RT[24]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[9] -= q; /* CH */
    wdot[13] -= q; /* CH4 */
    wdot[24] += q; /* C2H4 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[24];
    J[34] += dqdci;               /* dwdot[H]/d[H] */
    J[42] -= dqdci;               /* dwdot[CH]/d[H] */
    J[46] -= dqdci;               /* dwdot[CH4]/d[H] */
    J[57] += dqdci;               /* dwdot[C2H4]/d[H] */
    /* d()/d[CH] */
    dqdci =  + k_f*sc[13];
    J[298] += dqdci;              /* dwdot[H]/d[CH] */
    J[306] -= dqdci;              /* dwdot[CH]/d[CH] */
    J[310] -= dqdci;              /* dwdot[CH4]/d[CH] */
    J[321] += dqdci;              /* dwdot[C2H4]/d[CH] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[9];
    J[430] += dqdci;              /* dwdot[H]/d[CH4] */
    J[438] -= dqdci;              /* dwdot[CH]/d[CH4] */
    J[442] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    J[453] += dqdci;              /* dwdot[C2H4]/d[CH4] */
    /* d()/d[C2H4] */
    dqdci =  - k_r*sc[1];
    J[793] += dqdci;              /* dwdot[H]/d[C2H4] */
    J[801] -= dqdci;              /* dwdot[CH]/d[C2H4] */
    J[805] -= dqdci;              /* dwdot[CH4]/d[C2H4] */
    J[816] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
    /* d()/dT */
    J[1057] += dqdT;              /* dwdot[H]/dT */
    J[1065] -= dqdT;              /* dwdot[CH]/dT */
    J[1069] -= dqdT;              /* dwdot[CH4]/dT */
    J[1080] += dqdT;              /* dwdot[C2H4]/dT */

    /*reaction 137: CH + CO2 <=> HCO + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[15];
    k_f = prefactor_units[136] * fwd_A[136]
                * exp(fwd_beta[136] * tc[0] - activation_units[136] * fwd_Ea[136] * invT);
    dlnkfdT = fwd_beta[136] * invT + activation_units[136] * fwd_Ea[136] * invT2;
    /* reverse */
    phi_r = sc[14]*sc[16];
    Kc = exp(g_RT[9] - g_RT[14] + g_RT[15] - g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[9] + h_RT[15]) + (h_RT[14] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[9] -= q; /* CH */
    wdot[14] += q; /* CO */
    wdot[15] -= q; /* CO2 */
    wdot[16] += q; /* HCO */
    /* d()/d[CH] */
    dqdci =  + k_f*sc[15];
    J[306] -= dqdci;              /* dwdot[CH]/d[CH] */
    J[311] += dqdci;              /* dwdot[CO]/d[CH] */
    J[312] -= dqdci;              /* dwdot[CO2]/d[CH] */
    J[313] += dqdci;              /* dwdot[HCO]/d[CH] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[16];
    J[471] -= dqdci;              /* dwdot[CH]/d[CO] */
    J[476] += dqdci;              /* dwdot[CO]/d[CO] */
    J[477] -= dqdci;              /* dwdot[CO2]/d[CO] */
    J[478] += dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  + k_f*sc[9];
    J[504] -= dqdci;              /* dwdot[CH]/d[CO2] */
    J[509] += dqdci;              /* dwdot[CO]/d[CO2] */
    J[510] -= dqdci;              /* dwdot[CO2]/d[CO2] */
    J[511] += dqdci;              /* dwdot[HCO]/d[CO2] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[14];
    J[537] -= dqdci;              /* dwdot[CH]/d[HCO] */
    J[542] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[543] -= dqdci;              /* dwdot[CO2]/d[HCO] */
    J[544] += dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[1065] -= dqdT;              /* dwdot[CH]/dT */
    J[1070] += dqdT;              /* dwdot[CO]/dT */
    J[1071] -= dqdT;              /* dwdot[CO2]/dT */
    J[1072] += dqdT;              /* dwdot[HCO]/dT */

    /*reaction 138: CH + CH2O <=> H + CH2CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[17];
    k_f = prefactor_units[137] * fwd_A[137]
                * exp(fwd_beta[137] * tc[0] - activation_units[137] * fwd_Ea[137] * invT);
    dlnkfdT = fwd_beta[137] * invT + activation_units[137] * fwd_Ea[137] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[28];
    Kc = exp(-g_RT[1] + g_RT[9] + g_RT[17] - g_RT[28]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[9] + h_RT[17]) + (h_RT[1] + h_RT[28]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[9] -= q; /* CH */
    wdot[17] -= q; /* CH2O */
    wdot[28] += q; /* CH2CO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[28];
    J[34] += dqdci;               /* dwdot[H]/d[H] */
    J[42] -= dqdci;               /* dwdot[CH]/d[H] */
    J[50] -= dqdci;               /* dwdot[CH2O]/d[H] */
    J[61] += dqdci;               /* dwdot[CH2CO]/d[H] */
    /* d()/d[CH] */
    dqdci =  + k_f*sc[17];
    J[298] += dqdci;              /* dwdot[H]/d[CH] */
    J[306] -= dqdci;              /* dwdot[CH]/d[CH] */
    J[314] -= dqdci;              /* dwdot[CH2O]/d[CH] */
    J[325] += dqdci;              /* dwdot[CH2CO]/d[CH] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[9];
    J[562] += dqdci;              /* dwdot[H]/d[CH2O] */
    J[570] -= dqdci;              /* dwdot[CH]/d[CH2O] */
    J[578] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[589] += dqdci;              /* dwdot[CH2CO]/d[CH2O] */
    /* d()/d[CH2CO] */
    dqdci =  - k_r*sc[1];
    J[925] += dqdci;              /* dwdot[H]/d[CH2CO] */
    J[933] -= dqdci;              /* dwdot[CH]/d[CH2CO] */
    J[941] -= dqdci;              /* dwdot[CH2O]/d[CH2CO] */
    J[952] += dqdci;              /* dwdot[CH2CO]/d[CH2CO] */
    /* d()/dT */
    J[1057] += dqdT;              /* dwdot[H]/dT */
    J[1065] -= dqdT;              /* dwdot[CH]/dT */
    J[1073] -= dqdT;              /* dwdot[CH2O]/dT */
    J[1084] += dqdT;              /* dwdot[CH2CO]/dT */

    /*reaction 139: CH + HCCO <=> CO + C2H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[27];
    k_f = prefactor_units[138] * fwd_A[138]
                * exp(fwd_beta[138] * tc[0] - activation_units[138] * fwd_Ea[138] * invT);
    dlnkfdT = fwd_beta[138] * invT + activation_units[138] * fwd_Ea[138] * invT2;
    /* reverse */
    phi_r = sc[14]*sc[22];
    Kc = exp(g_RT[9] - g_RT[14] - g_RT[22] + g_RT[27]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[9] + h_RT[27]) + (h_RT[14] + h_RT[22]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[9] -= q; /* CH */
    wdot[14] += q; /* CO */
    wdot[22] += q; /* C2H2 */
    wdot[27] -= q; /* HCCO */
    /* d()/d[CH] */
    dqdci =  + k_f*sc[27];
    J[306] -= dqdci;              /* dwdot[CH]/d[CH] */
    J[311] += dqdci;              /* dwdot[CO]/d[CH] */
    J[319] += dqdci;              /* dwdot[C2H2]/d[CH] */
    J[324] -= dqdci;              /* dwdot[HCCO]/d[CH] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[22];
    J[471] -= dqdci;              /* dwdot[CH]/d[CO] */
    J[476] += dqdci;              /* dwdot[CO]/d[CO] */
    J[484] += dqdci;              /* dwdot[C2H2]/d[CO] */
    J[489] -= dqdci;              /* dwdot[HCCO]/d[CO] */
    /* d()/d[C2H2] */
    dqdci =  - k_r*sc[14];
    J[735] -= dqdci;              /* dwdot[CH]/d[C2H2] */
    J[740] += dqdci;              /* dwdot[CO]/d[C2H2] */
    J[748] += dqdci;              /* dwdot[C2H2]/d[C2H2] */
    J[753] -= dqdci;              /* dwdot[HCCO]/d[C2H2] */
    /* d()/d[HCCO] */
    dqdci =  + k_f*sc[9];
    J[900] -= dqdci;              /* dwdot[CH]/d[HCCO] */
    J[905] += dqdci;              /* dwdot[CO]/d[HCCO] */
    J[913] += dqdci;              /* dwdot[C2H2]/d[HCCO] */
    J[918] -= dqdci;              /* dwdot[HCCO]/d[HCCO] */
    /* d()/dT */
    J[1065] -= dqdT;              /* dwdot[CH]/dT */
    J[1070] += dqdT;              /* dwdot[CO]/dT */
    J[1078] += dqdT;              /* dwdot[C2H2]/dT */
    J[1083] -= dqdT;              /* dwdot[HCCO]/dT */

    /*reaction 140: CH2 + O2 <=> OH + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[10];
    k_f = prefactor_units[139] * fwd_A[139]
                * exp(fwd_beta[139] * tc[0] - activation_units[139] * fwd_Ea[139] * invT);
    dlnkfdT = fwd_beta[139] * invT + activation_units[139] * fwd_Ea[139] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[16];
    Kc = exp(g_RT[3] - g_RT[4] + g_RT[10] - g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[10]) + (h_RT[4] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[4] += q; /* OH */
    wdot[10] -= q; /* CH2 */
    wdot[16] += q; /* HCO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[10];
    J[102] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[103] += dqdci;              /* dwdot[OH]/d[O2] */
    J[109] -= dqdci;              /* dwdot[CH2]/d[O2] */
    J[115] += dqdci;              /* dwdot[HCO]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[16];
    J[135] -= dqdci;              /* dwdot[O2]/d[OH] */
    J[136] += dqdci;              /* dwdot[OH]/d[OH] */
    J[142] -= dqdci;              /* dwdot[CH2]/d[OH] */
    J[148] += dqdci;              /* dwdot[HCO]/d[OH] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[3];
    J[333] -= dqdci;              /* dwdot[O2]/d[CH2] */
    J[334] += dqdci;              /* dwdot[OH]/d[CH2] */
    J[340] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[346] += dqdci;              /* dwdot[HCO]/d[CH2] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[4];
    J[531] -= dqdci;              /* dwdot[O2]/d[HCO] */
    J[532] += dqdci;              /* dwdot[OH]/d[HCO] */
    J[538] -= dqdci;              /* dwdot[CH2]/d[HCO] */
    J[544] += dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[1059] -= dqdT;              /* dwdot[O2]/dT */
    J[1060] += dqdT;              /* dwdot[OH]/dT */
    J[1066] -= dqdT;              /* dwdot[CH2]/dT */
    J[1072] += dqdT;              /* dwdot[HCO]/dT */

    /*reaction 141: CH2 + H2 <=> H + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[10];
    k_f = prefactor_units[140] * fwd_A[140]
                * exp(fwd_beta[140] * tc[0] - activation_units[140] * fwd_Ea[140] * invT);
    dlnkfdT = fwd_beta[140] * invT + activation_units[140] * fwd_Ea[140] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[12];
    Kc = exp(g_RT[0] - g_RT[1] + g_RT[10] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[10]) + (h_RT[1] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[1] += q; /* H */
    wdot[10] -= q; /* CH2 */
    wdot[12] += q; /* CH3 */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[10];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += dqdci;                /* dwdot[H]/d[H2] */
    J[10] -= dqdci;               /* dwdot[CH2]/d[H2] */
    J[12] += dqdci;               /* dwdot[CH3]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[12];
    J[33] -= dqdci;               /* dwdot[H2]/d[H] */
    J[34] += dqdci;               /* dwdot[H]/d[H] */
    J[43] -= dqdci;               /* dwdot[CH2]/d[H] */
    J[45] += dqdci;               /* dwdot[CH3]/d[H] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[0];
    J[330] -= dqdci;              /* dwdot[H2]/d[CH2] */
    J[331] += dqdci;              /* dwdot[H]/d[CH2] */
    J[340] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[342] += dqdci;              /* dwdot[CH3]/d[CH2] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[1];
    J[396] -= dqdci;              /* dwdot[H2]/d[CH3] */
    J[397] += dqdci;              /* dwdot[H]/d[CH3] */
    J[406] -= dqdci;              /* dwdot[CH2]/d[CH3] */
    J[408] += dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[1056] -= dqdT;              /* dwdot[H2]/dT */
    J[1057] += dqdT;              /* dwdot[H]/dT */
    J[1066] -= dqdT;              /* dwdot[CH2]/dT */
    J[1068] += dqdT;              /* dwdot[CH3]/dT */

    /*reaction 142: 2.000000 CH2 <=> H2 + C2H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[10], 2.000000);
    k_f = prefactor_units[141] * fwd_A[141]
                * exp(fwd_beta[141] * tc[0] - activation_units[141] * fwd_Ea[141] * invT);
    dlnkfdT = fwd_beta[141] * invT + activation_units[141] * fwd_Ea[141] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[22];
    Kc = exp(-g_RT[0] + 2.000000*g_RT[10] - g_RT[22]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[10]) + (h_RT[0] + h_RT[22]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[10] -= 2 * q; /* CH2 */
    wdot[22] += q; /* C2H2 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[22];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[10] += -2 * dqdci;          /* dwdot[CH2]/d[H2] */
    J[22] += dqdci;               /* dwdot[C2H2]/d[H2] */
    /* d()/d[CH2] */
    dqdci =  + k_f*2.000000*sc[10];
    J[330] += dqdci;              /* dwdot[H2]/d[CH2] */
    J[340] += -2 * dqdci;         /* dwdot[CH2]/d[CH2] */
    J[352] += dqdci;              /* dwdot[C2H2]/d[CH2] */
    /* d()/d[C2H2] */
    dqdci =  - k_r*sc[0];
    J[726] += dqdci;              /* dwdot[H2]/d[C2H2] */
    J[736] += -2 * dqdci;         /* dwdot[CH2]/d[C2H2] */
    J[748] += dqdci;              /* dwdot[C2H2]/d[C2H2] */
    /* d()/dT */
    J[1056] += dqdT;              /* dwdot[H2]/dT */
    J[1066] += -2 * dqdT;         /* dwdot[CH2]/dT */
    J[1078] += dqdT;              /* dwdot[C2H2]/dT */

    /*reaction 143: CH2 + CH3 <=> H + C2H4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[10]*sc[12];
    k_f = prefactor_units[142] * fwd_A[142]
                * exp(fwd_beta[142] * tc[0] - activation_units[142] * fwd_Ea[142] * invT);
    dlnkfdT = fwd_beta[142] * invT + activation_units[142] * fwd_Ea[142] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[24];
    Kc = exp(-g_RT[1] + g_RT[10] + g_RT[12] - g_RT[24]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[10] + h_RT[12]) + (h_RT[1] + h_RT[24]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[10] -= q; /* CH2 */
    wdot[12] -= q; /* CH3 */
    wdot[24] += q; /* C2H4 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[24];
    J[34] += dqdci;               /* dwdot[H]/d[H] */
    J[43] -= dqdci;               /* dwdot[CH2]/d[H] */
    J[45] -= dqdci;               /* dwdot[CH3]/d[H] */
    J[57] += dqdci;               /* dwdot[C2H4]/d[H] */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[12];
    J[331] += dqdci;              /* dwdot[H]/d[CH2] */
    J[340] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[342] -= dqdci;              /* dwdot[CH3]/d[CH2] */
    J[354] += dqdci;              /* dwdot[C2H4]/d[CH2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[10];
    J[397] += dqdci;              /* dwdot[H]/d[CH3] */
    J[406] -= dqdci;              /* dwdot[CH2]/d[CH3] */
    J[408] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[420] += dqdci;              /* dwdot[C2H4]/d[CH3] */
    /* d()/d[C2H4] */
    dqdci =  - k_r*sc[1];
    J[793] += dqdci;              /* dwdot[H]/d[C2H4] */
    J[802] -= dqdci;              /* dwdot[CH2]/d[C2H4] */
    J[804] -= dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[816] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
    /* d()/dT */
    J[1057] += dqdT;              /* dwdot[H]/dT */
    J[1066] -= dqdT;              /* dwdot[CH2]/dT */
    J[1068] -= dqdT;              /* dwdot[CH3]/dT */
    J[1080] += dqdT;              /* dwdot[C2H4]/dT */

    /*reaction 144: CH2 + CH4 <=> 2.000000 CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[10]*sc[13];
    k_f = prefactor_units[143] * fwd_A[143]
                * exp(fwd_beta[143] * tc[0] - activation_units[143] * fwd_Ea[143] * invT);
    dlnkfdT = fwd_beta[143] * invT + activation_units[143] * fwd_Ea[143] * invT2;
    /* reverse */
    phi_r = pow(sc[12], 2.000000);
    Kc = exp(g_RT[10] - 2.000000*g_RT[12] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[10] + h_RT[13]) + (2.000000*h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[10] -= q; /* CH2 */
    wdot[12] += 2 * q; /* CH3 */
    wdot[13] -= q; /* CH4 */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[13];
    J[340] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[342] += 2 * dqdci;          /* dwdot[CH3]/d[CH2] */
    J[343] -= dqdci;              /* dwdot[CH4]/d[CH2] */
    /* d()/d[CH3] */
    dqdci =  - k_r*2.000000*sc[12];
    J[406] -= dqdci;              /* dwdot[CH2]/d[CH3] */
    J[408] += 2 * dqdci;          /* dwdot[CH3]/d[CH3] */
    J[409] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[10];
    J[439] -= dqdci;              /* dwdot[CH2]/d[CH4] */
    J[441] += 2 * dqdci;          /* dwdot[CH3]/d[CH4] */
    J[442] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[1066] -= dqdT;              /* dwdot[CH2]/dT */
    J[1068] += 2 * dqdT;          /* dwdot[CH3]/dT */
    J[1069] -= dqdT;              /* dwdot[CH4]/dT */

    /*reaction 145: CH2 + HCCO <=> C2H3 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[10]*sc[27];
    k_f = prefactor_units[144] * fwd_A[144]
                * exp(fwd_beta[144] * tc[0] - activation_units[144] * fwd_Ea[144] * invT);
    dlnkfdT = fwd_beta[144] * invT + activation_units[144] * fwd_Ea[144] * invT2;
    /* reverse */
    phi_r = sc[14]*sc[23];
    Kc = exp(g_RT[10] - g_RT[14] - g_RT[23] + g_RT[27]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[10] + h_RT[27]) + (h_RT[14] + h_RT[23]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[10] -= q; /* CH2 */
    wdot[14] += q; /* CO */
    wdot[23] += q; /* C2H3 */
    wdot[27] -= q; /* HCCO */
    /* d()/d[CH2] */
    dqdci =  + k_f*sc[27];
    J[340] -= dqdci;              /* dwdot[CH2]/d[CH2] */
    J[344] += dqdci;              /* dwdot[CO]/d[CH2] */
    J[353] += dqdci;              /* dwdot[C2H3]/d[CH2] */
    J[357] -= dqdci;              /* dwdot[HCCO]/d[CH2] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[23];
    J[472] -= dqdci;              /* dwdot[CH2]/d[CO] */
    J[476] += dqdci;              /* dwdot[CO]/d[CO] */
    J[485] += dqdci;              /* dwdot[C2H3]/d[CO] */
    J[489] -= dqdci;              /* dwdot[HCCO]/d[CO] */
    /* d()/d[C2H3] */
    dqdci =  - k_r*sc[14];
    J[769] -= dqdci;              /* dwdot[CH2]/d[C2H3] */
    J[773] += dqdci;              /* dwdot[CO]/d[C2H3] */
    J[782] += dqdci;              /* dwdot[C2H3]/d[C2H3] */
    J[786] -= dqdci;              /* dwdot[HCCO]/d[C2H3] */
    /* d()/d[HCCO] */
    dqdci =  + k_f*sc[10];
    J[901] -= dqdci;              /* dwdot[CH2]/d[HCCO] */
    J[905] += dqdci;              /* dwdot[CO]/d[HCCO] */
    J[914] += dqdci;              /* dwdot[C2H3]/d[HCCO] */
    J[918] -= dqdci;              /* dwdot[HCCO]/d[HCCO] */
    /* d()/dT */
    J[1066] -= dqdT;              /* dwdot[CH2]/dT */
    J[1070] += dqdT;              /* dwdot[CO]/dT */
    J[1079] += dqdT;              /* dwdot[C2H3]/dT */
    J[1083] -= dqdT;              /* dwdot[HCCO]/dT */

    /*reaction 146: CH2(S) + N2 <=> CH2 + N2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[30];
    k_f = prefactor_units[145] * fwd_A[145]
                * exp(fwd_beta[145] * tc[0] - activation_units[145] * fwd_Ea[145] * invT);
    dlnkfdT = fwd_beta[145] * invT + activation_units[145] * fwd_Ea[145] * invT2;
    /* reverse */
    phi_r = sc[10]*sc[30];
    Kc = exp(-g_RT[10] + g_RT[11] + g_RT[30] - g_RT[30]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[11] + h_RT[30]) + (h_RT[10] + h_RT[30]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[10] += q; /* CH2 */
    wdot[11] -= q; /* CH2(S) */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[30];
    J[340] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[341] -= dqdci;              /* dwdot[CH2(S)]/d[CH2] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[30];
    J[373] += dqdci;              /* dwdot[CH2]/d[CH2(S)] */
    J[374] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    /* d()/d[N2] */
    dqdci =  + k_f*sc[11] - k_r*sc[10];
    J[1000] += dqdci;             /* dwdot[CH2]/d[N2] */
    J[1001] -= dqdci;             /* dwdot[CH2(S)]/d[N2] */
    /* d()/dT */
    J[1066] += dqdT;              /* dwdot[CH2]/dT */
    J[1067] -= dqdT;              /* dwdot[CH2(S)]/dT */

    /*reaction 147: CH2(S) + AR <=> CH2 + AR */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[31];
    k_f = prefactor_units[146] * fwd_A[146]
                * exp(fwd_beta[146] * tc[0] - activation_units[146] * fwd_Ea[146] * invT);
    dlnkfdT = fwd_beta[146] * invT + activation_units[146] * fwd_Ea[146] * invT2;
    /* reverse */
    phi_r = sc[10]*sc[31];
    Kc = exp(-g_RT[10] + g_RT[11] + g_RT[31] - g_RT[31]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[11] + h_RT[31]) + (h_RT[10] + h_RT[31]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[10] += q; /* CH2 */
    wdot[11] -= q; /* CH2(S) */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[31];
    J[340] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[341] -= dqdci;              /* dwdot[CH2(S)]/d[CH2] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[31];
    J[373] += dqdci;              /* dwdot[CH2]/d[CH2(S)] */
    J[374] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    /* d()/d[AR] */
    dqdci =  + k_f*sc[11] - k_r*sc[10];
    J[1033] += dqdci;             /* dwdot[CH2]/d[AR] */
    J[1034] -= dqdci;             /* dwdot[CH2(S)]/d[AR] */
    /* d()/dT */
    J[1066] += dqdT;              /* dwdot[CH2]/dT */
    J[1067] -= dqdT;              /* dwdot[CH2(S)]/dT */

    /*reaction 148: CH2(S) + O2 <=> H + OH + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[11];
    k_f = prefactor_units[147] * fwd_A[147]
                * exp(fwd_beta[147] * tc[0] - activation_units[147] * fwd_Ea[147] * invT);
    dlnkfdT = fwd_beta[147] * invT + activation_units[147] * fwd_Ea[147] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[4]*sc[14];
    Kc = refC * exp(-g_RT[1] + g_RT[3] - g_RT[4] + g_RT[11] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[11]) + (h_RT[1] + h_RT[4] + h_RT[14]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[3] -= q; /* O2 */
    wdot[4] += q; /* OH */
    wdot[11] -= q; /* CH2(S) */
    wdot[14] += q; /* CO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[4]*sc[14];
    J[34] += dqdci;               /* dwdot[H]/d[H] */
    J[36] -= dqdci;               /* dwdot[O2]/d[H] */
    J[37] += dqdci;               /* dwdot[OH]/d[H] */
    J[44] -= dqdci;               /* dwdot[CH2(S)]/d[H] */
    J[47] += dqdci;               /* dwdot[CO]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[11];
    J[100] += dqdci;              /* dwdot[H]/d[O2] */
    J[102] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[103] += dqdci;              /* dwdot[OH]/d[O2] */
    J[110] -= dqdci;              /* dwdot[CH2(S)]/d[O2] */
    J[113] += dqdci;              /* dwdot[CO]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[1]*sc[14];
    J[133] += dqdci;              /* dwdot[H]/d[OH] */
    J[135] -= dqdci;              /* dwdot[O2]/d[OH] */
    J[136] += dqdci;              /* dwdot[OH]/d[OH] */
    J[143] -= dqdci;              /* dwdot[CH2(S)]/d[OH] */
    J[146] += dqdci;              /* dwdot[CO]/d[OH] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[3];
    J[364] += dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[366] -= dqdci;              /* dwdot[O2]/d[CH2(S)] */
    J[367] += dqdci;              /* dwdot[OH]/d[CH2(S)] */
    J[374] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[377] += dqdci;              /* dwdot[CO]/d[CH2(S)] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[1]*sc[4];
    J[463] += dqdci;              /* dwdot[H]/d[CO] */
    J[465] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[466] += dqdci;              /* dwdot[OH]/d[CO] */
    J[473] -= dqdci;              /* dwdot[CH2(S)]/d[CO] */
    J[476] += dqdci;              /* dwdot[CO]/d[CO] */
    /* d()/dT */
    J[1057] += dqdT;              /* dwdot[H]/dT */
    J[1059] -= dqdT;              /* dwdot[O2]/dT */
    J[1060] += dqdT;              /* dwdot[OH]/dT */
    J[1067] -= dqdT;              /* dwdot[CH2(S)]/dT */
    J[1070] += dqdT;              /* dwdot[CO]/dT */

    /*reaction 149: CH2(S) + O2 <=> CO + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[11];
    k_f = prefactor_units[148] * fwd_A[148]
                * exp(fwd_beta[148] * tc[0] - activation_units[148] * fwd_Ea[148] * invT);
    dlnkfdT = fwd_beta[148] * invT + activation_units[148] * fwd_Ea[148] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[14];
    Kc = exp(g_RT[3] - g_RT[5] + g_RT[11] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[11]) + (h_RT[5] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[5] += q; /* H2O */
    wdot[11] -= q; /* CH2(S) */
    wdot[14] += q; /* CO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[11];
    J[102] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[104] += dqdci;              /* dwdot[H2O]/d[O2] */
    J[110] -= dqdci;              /* dwdot[CH2(S)]/d[O2] */
    J[113] += dqdci;              /* dwdot[CO]/d[O2] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[14];
    J[168] -= dqdci;              /* dwdot[O2]/d[H2O] */
    J[170] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[176] -= dqdci;              /* dwdot[CH2(S)]/d[H2O] */
    J[179] += dqdci;              /* dwdot[CO]/d[H2O] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[3];
    J[366] -= dqdci;              /* dwdot[O2]/d[CH2(S)] */
    J[368] += dqdci;              /* dwdot[H2O]/d[CH2(S)] */
    J[374] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[377] += dqdci;              /* dwdot[CO]/d[CH2(S)] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[5];
    J[465] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[467] += dqdci;              /* dwdot[H2O]/d[CO] */
    J[473] -= dqdci;              /* dwdot[CH2(S)]/d[CO] */
    J[476] += dqdci;              /* dwdot[CO]/d[CO] */
    /* d()/dT */
    J[1059] -= dqdT;              /* dwdot[O2]/dT */
    J[1061] += dqdT;              /* dwdot[H2O]/dT */
    J[1067] -= dqdT;              /* dwdot[CH2(S)]/dT */
    J[1070] += dqdT;              /* dwdot[CO]/dT */

    /*reaction 150: CH2(S) + H2 <=> CH3 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[11];
    k_f = prefactor_units[149] * fwd_A[149]
                * exp(fwd_beta[149] * tc[0] - activation_units[149] * fwd_Ea[149] * invT);
    dlnkfdT = fwd_beta[149] * invT + activation_units[149] * fwd_Ea[149] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[12];
    Kc = exp(g_RT[0] - g_RT[1] + g_RT[11] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[11]) + (h_RT[1] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[1] += q; /* H */
    wdot[11] -= q; /* CH2(S) */
    wdot[12] += q; /* CH3 */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[11];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += dqdci;                /* dwdot[H]/d[H2] */
    J[11] -= dqdci;               /* dwdot[CH2(S)]/d[H2] */
    J[12] += dqdci;               /* dwdot[CH3]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[12];
    J[33] -= dqdci;               /* dwdot[H2]/d[H] */
    J[34] += dqdci;               /* dwdot[H]/d[H] */
    J[44] -= dqdci;               /* dwdot[CH2(S)]/d[H] */
    J[45] += dqdci;               /* dwdot[CH3]/d[H] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[0];
    J[363] -= dqdci;              /* dwdot[H2]/d[CH2(S)] */
    J[364] += dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[374] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[375] += dqdci;              /* dwdot[CH3]/d[CH2(S)] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[1];
    J[396] -= dqdci;              /* dwdot[H2]/d[CH3] */
    J[397] += dqdci;              /* dwdot[H]/d[CH3] */
    J[407] -= dqdci;              /* dwdot[CH2(S)]/d[CH3] */
    J[408] += dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[1056] -= dqdT;              /* dwdot[H2]/dT */
    J[1057] += dqdT;              /* dwdot[H]/dT */
    J[1067] -= dqdT;              /* dwdot[CH2(S)]/dT */
    J[1068] += dqdT;              /* dwdot[CH3]/dT */

    /*reaction 151: CH2(S) + H2O <=> CH2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[11];
    k_f = prefactor_units[150] * fwd_A[150]
                * exp(fwd_beta[150] * tc[0] - activation_units[150] * fwd_Ea[150] * invT);
    dlnkfdT = fwd_beta[150] * invT + activation_units[150] * fwd_Ea[150] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[10];
    Kc = exp(g_RT[5] - g_RT[5] - g_RT[10] + g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[11]) + (h_RT[5] + h_RT[10]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[10] += q; /* CH2 */
    wdot[11] -= q; /* CH2(S) */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[11] - k_r*sc[10];
    J[175] += dqdci;              /* dwdot[CH2]/d[H2O] */
    J[176] -= dqdci;              /* dwdot[CH2(S)]/d[H2O] */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[5];
    J[340] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[341] -= dqdci;              /* dwdot[CH2(S)]/d[CH2] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[5];
    J[373] += dqdci;              /* dwdot[CH2]/d[CH2(S)] */
    J[374] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    /* d()/dT */
    J[1066] += dqdT;              /* dwdot[CH2]/dT */
    J[1067] -= dqdT;              /* dwdot[CH2(S)]/dT */

    /*reaction 152: CH2(S) + CH3 <=> H + C2H4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[12];
    k_f = prefactor_units[151] * fwd_A[151]
                * exp(fwd_beta[151] * tc[0] - activation_units[151] * fwd_Ea[151] * invT);
    dlnkfdT = fwd_beta[151] * invT + activation_units[151] * fwd_Ea[151] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[24];
    Kc = exp(-g_RT[1] + g_RT[11] + g_RT[12] - g_RT[24]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[11] + h_RT[12]) + (h_RT[1] + h_RT[24]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[11] -= q; /* CH2(S) */
    wdot[12] -= q; /* CH3 */
    wdot[24] += q; /* C2H4 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[24];
    J[34] += dqdci;               /* dwdot[H]/d[H] */
    J[44] -= dqdci;               /* dwdot[CH2(S)]/d[H] */
    J[45] -= dqdci;               /* dwdot[CH3]/d[H] */
    J[57] += dqdci;               /* dwdot[C2H4]/d[H] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[12];
    J[364] += dqdci;              /* dwdot[H]/d[CH2(S)] */
    J[374] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[375] -= dqdci;              /* dwdot[CH3]/d[CH2(S)] */
    J[387] += dqdci;              /* dwdot[C2H4]/d[CH2(S)] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[11];
    J[397] += dqdci;              /* dwdot[H]/d[CH3] */
    J[407] -= dqdci;              /* dwdot[CH2(S)]/d[CH3] */
    J[408] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[420] += dqdci;              /* dwdot[C2H4]/d[CH3] */
    /* d()/d[C2H4] */
    dqdci =  - k_r*sc[1];
    J[793] += dqdci;              /* dwdot[H]/d[C2H4] */
    J[803] -= dqdci;              /* dwdot[CH2(S)]/d[C2H4] */
    J[804] -= dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[816] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
    /* d()/dT */
    J[1057] += dqdT;              /* dwdot[H]/dT */
    J[1067] -= dqdT;              /* dwdot[CH2(S)]/dT */
    J[1068] -= dqdT;              /* dwdot[CH3]/dT */
    J[1080] += dqdT;              /* dwdot[C2H4]/dT */

    /*reaction 153: CH2(S) + CH4 <=> 2.000000 CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[13];
    k_f = prefactor_units[152] * fwd_A[152]
                * exp(fwd_beta[152] * tc[0] - activation_units[152] * fwd_Ea[152] * invT);
    dlnkfdT = fwd_beta[152] * invT + activation_units[152] * fwd_Ea[152] * invT2;
    /* reverse */
    phi_r = pow(sc[12], 2.000000);
    Kc = exp(g_RT[11] - 2.000000*g_RT[12] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[11] + h_RT[13]) + (2.000000*h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[11] -= q; /* CH2(S) */
    wdot[12] += 2 * q; /* CH3 */
    wdot[13] -= q; /* CH4 */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[13];
    J[374] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[375] += 2 * dqdci;          /* dwdot[CH3]/d[CH2(S)] */
    J[376] -= dqdci;              /* dwdot[CH4]/d[CH2(S)] */
    /* d()/d[CH3] */
    dqdci =  - k_r*2.000000*sc[12];
    J[407] -= dqdci;              /* dwdot[CH2(S)]/d[CH3] */
    J[408] += 2 * dqdci;          /* dwdot[CH3]/d[CH3] */
    J[409] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[11];
    J[440] -= dqdci;              /* dwdot[CH2(S)]/d[CH4] */
    J[441] += 2 * dqdci;          /* dwdot[CH3]/d[CH4] */
    J[442] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[1067] -= dqdT;              /* dwdot[CH2(S)]/dT */
    J[1068] += 2 * dqdT;          /* dwdot[CH3]/dT */
    J[1069] -= dqdT;              /* dwdot[CH4]/dT */

    /*reaction 154: CH2(S) + CO <=> CH2 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[14];
    k_f = prefactor_units[153] * fwd_A[153]
                * exp(fwd_beta[153] * tc[0] - activation_units[153] * fwd_Ea[153] * invT);
    dlnkfdT = fwd_beta[153] * invT + activation_units[153] * fwd_Ea[153] * invT2;
    /* reverse */
    phi_r = sc[10]*sc[14];
    Kc = exp(-g_RT[10] + g_RT[11] + g_RT[14] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[11] + h_RT[14]) + (h_RT[10] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[10] += q; /* CH2 */
    wdot[11] -= q; /* CH2(S) */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[14];
    J[340] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[341] -= dqdci;              /* dwdot[CH2(S)]/d[CH2] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[14];
    J[373] += dqdci;              /* dwdot[CH2]/d[CH2(S)] */
    J[374] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[11] - k_r*sc[10];
    J[472] += dqdci;              /* dwdot[CH2]/d[CO] */
    J[473] -= dqdci;              /* dwdot[CH2(S)]/d[CO] */
    /* d()/dT */
    J[1066] += dqdT;              /* dwdot[CH2]/dT */
    J[1067] -= dqdT;              /* dwdot[CH2(S)]/dT */

    /*reaction 155: CH2(S) + CO2 <=> CH2 + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[15];
    k_f = prefactor_units[154] * fwd_A[154]
                * exp(fwd_beta[154] * tc[0] - activation_units[154] * fwd_Ea[154] * invT);
    dlnkfdT = fwd_beta[154] * invT + activation_units[154] * fwd_Ea[154] * invT2;
    /* reverse */
    phi_r = sc[10]*sc[15];
    Kc = exp(-g_RT[10] + g_RT[11] + g_RT[15] - g_RT[15]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[11] + h_RT[15]) + (h_RT[10] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[10] += q; /* CH2 */
    wdot[11] -= q; /* CH2(S) */
    /* d()/d[CH2] */
    dqdci =  - k_r*sc[15];
    J[340] += dqdci;              /* dwdot[CH2]/d[CH2] */
    J[341] -= dqdci;              /* dwdot[CH2(S)]/d[CH2] */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[15];
    J[373] += dqdci;              /* dwdot[CH2]/d[CH2(S)] */
    J[374] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    /* d()/d[CO2] */
    dqdci =  + k_f*sc[11] - k_r*sc[10];
    J[505] += dqdci;              /* dwdot[CH2]/d[CO2] */
    J[506] -= dqdci;              /* dwdot[CH2(S)]/d[CO2] */
    /* d()/dT */
    J[1066] += dqdT;              /* dwdot[CH2]/dT */
    J[1067] -= dqdT;              /* dwdot[CH2(S)]/dT */

    /*reaction 156: CH2(S) + CO2 <=> CO + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[15];
    k_f = prefactor_units[155] * fwd_A[155]
                * exp(fwd_beta[155] * tc[0] - activation_units[155] * fwd_Ea[155] * invT);
    dlnkfdT = fwd_beta[155] * invT + activation_units[155] * fwd_Ea[155] * invT2;
    /* reverse */
    phi_r = sc[14]*sc[17];
    Kc = exp(g_RT[11] - g_RT[14] + g_RT[15] - g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[11] + h_RT[15]) + (h_RT[14] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[11] -= q; /* CH2(S) */
    wdot[14] += q; /* CO */
    wdot[15] -= q; /* CO2 */
    wdot[17] += q; /* CH2O */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[15];
    J[374] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[377] += dqdci;              /* dwdot[CO]/d[CH2(S)] */
    J[378] -= dqdci;              /* dwdot[CO2]/d[CH2(S)] */
    J[380] += dqdci;              /* dwdot[CH2O]/d[CH2(S)] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[17];
    J[473] -= dqdci;              /* dwdot[CH2(S)]/d[CO] */
    J[476] += dqdci;              /* dwdot[CO]/d[CO] */
    J[477] -= dqdci;              /* dwdot[CO2]/d[CO] */
    J[479] += dqdci;              /* dwdot[CH2O]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  + k_f*sc[11];
    J[506] -= dqdci;              /* dwdot[CH2(S)]/d[CO2] */
    J[509] += dqdci;              /* dwdot[CO]/d[CO2] */
    J[510] -= dqdci;              /* dwdot[CO2]/d[CO2] */
    J[512] += dqdci;              /* dwdot[CH2O]/d[CO2] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[14];
    J[572] -= dqdci;              /* dwdot[CH2(S)]/d[CH2O] */
    J[575] += dqdci;              /* dwdot[CO]/d[CH2O] */
    J[576] -= dqdci;              /* dwdot[CO2]/d[CH2O] */
    J[578] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[1067] -= dqdT;              /* dwdot[CH2(S)]/dT */
    J[1070] += dqdT;              /* dwdot[CO]/dT */
    J[1071] -= dqdT;              /* dwdot[CO2]/dT */
    J[1073] += dqdT;              /* dwdot[CH2O]/dT */

    /*reaction 157: CH2(S) + C2H6 <=> CH3 + C2H5 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[26];
    k_f = prefactor_units[156] * fwd_A[156]
                * exp(fwd_beta[156] * tc[0] - activation_units[156] * fwd_Ea[156] * invT);
    dlnkfdT = fwd_beta[156] * invT + activation_units[156] * fwd_Ea[156] * invT2;
    /* reverse */
    phi_r = sc[12]*sc[25];
    Kc = exp(g_RT[11] - g_RT[12] - g_RT[25] + g_RT[26]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[11] + h_RT[26]) + (h_RT[12] + h_RT[25]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[11] -= q; /* CH2(S) */
    wdot[12] += q; /* CH3 */
    wdot[25] += q; /* C2H5 */
    wdot[26] -= q; /* C2H6 */
    /* d()/d[CH2(S)] */
    dqdci =  + k_f*sc[26];
    J[374] -= dqdci;              /* dwdot[CH2(S)]/d[CH2(S)] */
    J[375] += dqdci;              /* dwdot[CH3]/d[CH2(S)] */
    J[388] += dqdci;              /* dwdot[C2H5]/d[CH2(S)] */
    J[389] -= dqdci;              /* dwdot[C2H6]/d[CH2(S)] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[25];
    J[407] -= dqdci;              /* dwdot[CH2(S)]/d[CH3] */
    J[408] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[421] += dqdci;              /* dwdot[C2H5]/d[CH3] */
    J[422] -= dqdci;              /* dwdot[C2H6]/d[CH3] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[12];
    J[836] -= dqdci;              /* dwdot[CH2(S)]/d[C2H5] */
    J[837] += dqdci;              /* dwdot[CH3]/d[C2H5] */
    J[850] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[851] -= dqdci;              /* dwdot[C2H6]/d[C2H5] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[11];
    J[869] -= dqdci;              /* dwdot[CH2(S)]/d[C2H6] */
    J[870] += dqdci;              /* dwdot[CH3]/d[C2H6] */
    J[883] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[884] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[1067] -= dqdT;              /* dwdot[CH2(S)]/dT */
    J[1068] += dqdT;              /* dwdot[CH3]/dT */
    J[1081] += dqdT;              /* dwdot[C2H5]/dT */
    J[1082] -= dqdT;              /* dwdot[C2H6]/dT */

    /*reaction 158: CH3 + O2 <=> O + CH3O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[12];
    k_f = prefactor_units[157] * fwd_A[157]
                * exp(fwd_beta[157] * tc[0] - activation_units[157] * fwd_Ea[157] * invT);
    dlnkfdT = fwd_beta[157] * invT + activation_units[157] * fwd_Ea[157] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[19];
    Kc = exp(-g_RT[2] + g_RT[3] + g_RT[12] - g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[12]) + (h_RT[2] + h_RT[19]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* O */
    wdot[3] -= q; /* O2 */
    wdot[12] -= q; /* CH3 */
    wdot[19] += q; /* CH3O */
    /* d()/d[O] */
    dqdci =  - k_r*sc[19];
    J[68] += dqdci;               /* dwdot[O]/d[O] */
    J[69] -= dqdci;               /* dwdot[O2]/d[O] */
    J[78] -= dqdci;               /* dwdot[CH3]/d[O] */
    J[85] += dqdci;               /* dwdot[CH3O]/d[O] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[12];
    J[101] += dqdci;              /* dwdot[O]/d[O2] */
    J[102] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[111] -= dqdci;              /* dwdot[CH3]/d[O2] */
    J[118] += dqdci;              /* dwdot[CH3O]/d[O2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[3];
    J[398] += dqdci;              /* dwdot[O]/d[CH3] */
    J[399] -= dqdci;              /* dwdot[O2]/d[CH3] */
    J[408] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[415] += dqdci;              /* dwdot[CH3O]/d[CH3] */
    /* d()/d[CH3O] */
    dqdci =  - k_r*sc[2];
    J[629] += dqdci;              /* dwdot[O]/d[CH3O] */
    J[630] -= dqdci;              /* dwdot[O2]/d[CH3O] */
    J[639] -= dqdci;              /* dwdot[CH3]/d[CH3O] */
    J[646] += dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[1058] += dqdT;              /* dwdot[O]/dT */
    J[1059] -= dqdT;              /* dwdot[O2]/dT */
    J[1068] -= dqdT;              /* dwdot[CH3]/dT */
    J[1075] += dqdT;              /* dwdot[CH3O]/dT */

    /*reaction 159: CH3 + O2 <=> OH + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[12];
    k_f = prefactor_units[158] * fwd_A[158]
                * exp(fwd_beta[158] * tc[0] - activation_units[158] * fwd_Ea[158] * invT);
    dlnkfdT = fwd_beta[158] * invT + activation_units[158] * fwd_Ea[158] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[17];
    Kc = exp(g_RT[3] - g_RT[4] + g_RT[12] - g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[12]) + (h_RT[4] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[4] += q; /* OH */
    wdot[12] -= q; /* CH3 */
    wdot[17] += q; /* CH2O */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[12];
    J[102] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[103] += dqdci;              /* dwdot[OH]/d[O2] */
    J[111] -= dqdci;              /* dwdot[CH3]/d[O2] */
    J[116] += dqdci;              /* dwdot[CH2O]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[17];
    J[135] -= dqdci;              /* dwdot[O2]/d[OH] */
    J[136] += dqdci;              /* dwdot[OH]/d[OH] */
    J[144] -= dqdci;              /* dwdot[CH3]/d[OH] */
    J[149] += dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[3];
    J[399] -= dqdci;              /* dwdot[O2]/d[CH3] */
    J[400] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[408] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[413] += dqdci;              /* dwdot[CH2O]/d[CH3] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[4];
    J[564] -= dqdci;              /* dwdot[O2]/d[CH2O] */
    J[565] += dqdci;              /* dwdot[OH]/d[CH2O] */
    J[573] -= dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[578] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[1059] -= dqdT;              /* dwdot[O2]/dT */
    J[1060] += dqdT;              /* dwdot[OH]/dT */
    J[1068] -= dqdT;              /* dwdot[CH3]/dT */
    J[1073] += dqdT;              /* dwdot[CH2O]/dT */

    /*reaction 160: CH3 + H2O2 <=> HO2 + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[7]*sc[12];
    k_f = prefactor_units[159] * fwd_A[159]
                * exp(fwd_beta[159] * tc[0] - activation_units[159] * fwd_Ea[159] * invT);
    dlnkfdT = fwd_beta[159] * invT + activation_units[159] * fwd_Ea[159] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[13];
    Kc = exp(-g_RT[6] + g_RT[7] + g_RT[12] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7] + h_RT[12]) + (h_RT[6] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    wdot[12] -= q; /* CH3 */
    wdot[13] += q; /* CH4 */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[13];
    J[204] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[205] -= dqdci;              /* dwdot[H2O2]/d[HO2] */
    J[210] -= dqdci;              /* dwdot[CH3]/d[HO2] */
    J[211] += dqdci;              /* dwdot[CH4]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[12];
    J[237] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[238] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    J[243] -= dqdci;              /* dwdot[CH3]/d[H2O2] */
    J[244] += dqdci;              /* dwdot[CH4]/d[H2O2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[7];
    J[402] += dqdci;              /* dwdot[HO2]/d[CH3] */
    J[403] -= dqdci;              /* dwdot[H2O2]/d[CH3] */
    J[408] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[409] += dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[6];
    J[435] += dqdci;              /* dwdot[HO2]/d[CH4] */
    J[436] -= dqdci;              /* dwdot[H2O2]/d[CH4] */
    J[441] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[442] += dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[1062] += dqdT;              /* dwdot[HO2]/dT */
    J[1063] -= dqdT;              /* dwdot[H2O2]/dT */
    J[1068] -= dqdT;              /* dwdot[CH3]/dT */
    J[1069] += dqdT;              /* dwdot[CH4]/dT */

    /*reaction 161: 2.000000 CH3 <=> H + C2H5 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[12], 2.000000);
    k_f = prefactor_units[160] * fwd_A[160]
                * exp(fwd_beta[160] * tc[0] - activation_units[160] * fwd_Ea[160] * invT);
    dlnkfdT = fwd_beta[160] * invT + activation_units[160] * fwd_Ea[160] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[25];
    Kc = exp(-g_RT[1] + 2.000000*g_RT[12] - g_RT[25]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[12]) + (h_RT[1] + h_RT[25]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[12] -= 2 * q; /* CH3 */
    wdot[25] += q; /* C2H5 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[25];
    J[34] += dqdci;               /* dwdot[H]/d[H] */
    J[45] += -2 * dqdci;          /* dwdot[CH3]/d[H] */
    J[58] += dqdci;               /* dwdot[C2H5]/d[H] */
    /* d()/d[CH3] */
    dqdci =  + k_f*2.000000*sc[12];
    J[397] += dqdci;              /* dwdot[H]/d[CH3] */
    J[408] += -2 * dqdci;         /* dwdot[CH3]/d[CH3] */
    J[421] += dqdci;              /* dwdot[C2H5]/d[CH3] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[1];
    J[826] += dqdci;              /* dwdot[H]/d[C2H5] */
    J[837] += -2 * dqdci;         /* dwdot[CH3]/d[C2H5] */
    J[850] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    /* d()/dT */
    J[1057] += dqdT;              /* dwdot[H]/dT */
    J[1068] += -2 * dqdT;         /* dwdot[CH3]/dT */
    J[1081] += dqdT;              /* dwdot[C2H5]/dT */

    /*reaction 162: CH3 + HCO <=> CH4 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[12]*sc[16];
    k_f = prefactor_units[161] * fwd_A[161]
                * exp(fwd_beta[161] * tc[0] - activation_units[161] * fwd_Ea[161] * invT);
    dlnkfdT = fwd_beta[161] * invT + activation_units[161] * fwd_Ea[161] * invT2;
    /* reverse */
    phi_r = sc[13]*sc[14];
    Kc = exp(g_RT[12] - g_RT[13] - g_RT[14] + g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[12] + h_RT[16]) + (h_RT[13] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[12] -= q; /* CH3 */
    wdot[13] += q; /* CH4 */
    wdot[14] += q; /* CO */
    wdot[16] -= q; /* HCO */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[16];
    J[408] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[409] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[410] += dqdci;              /* dwdot[CO]/d[CH3] */
    J[412] -= dqdci;              /* dwdot[HCO]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[14];
    J[441] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[442] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[443] += dqdci;              /* dwdot[CO]/d[CH4] */
    J[445] -= dqdci;              /* dwdot[HCO]/d[CH4] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[13];
    J[474] -= dqdci;              /* dwdot[CH3]/d[CO] */
    J[475] += dqdci;              /* dwdot[CH4]/d[CO] */
    J[476] += dqdci;              /* dwdot[CO]/d[CO] */
    J[478] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[12];
    J[540] -= dqdci;              /* dwdot[CH3]/d[HCO] */
    J[541] += dqdci;              /* dwdot[CH4]/d[HCO] */
    J[542] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[544] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[1068] -= dqdT;              /* dwdot[CH3]/dT */
    J[1069] += dqdT;              /* dwdot[CH4]/dT */
    J[1070] += dqdT;              /* dwdot[CO]/dT */
    J[1072] -= dqdT;              /* dwdot[HCO]/dT */

    /*reaction 163: CH3 + CH2O <=> HCO + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[12]*sc[17];
    k_f = prefactor_units[162] * fwd_A[162]
                * exp(fwd_beta[162] * tc[0] - activation_units[162] * fwd_Ea[162] * invT);
    dlnkfdT = fwd_beta[162] * invT + activation_units[162] * fwd_Ea[162] * invT2;
    /* reverse */
    phi_r = sc[13]*sc[16];
    Kc = exp(g_RT[12] - g_RT[13] - g_RT[16] + g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[12] + h_RT[17]) + (h_RT[13] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[12] -= q; /* CH3 */
    wdot[13] += q; /* CH4 */
    wdot[16] += q; /* HCO */
    wdot[17] -= q; /* CH2O */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[17];
    J[408] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[409] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[412] += dqdci;              /* dwdot[HCO]/d[CH3] */
    J[413] -= dqdci;              /* dwdot[CH2O]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[16];
    J[441] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[442] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[445] += dqdci;              /* dwdot[HCO]/d[CH4] */
    J[446] -= dqdci;              /* dwdot[CH2O]/d[CH4] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[13];
    J[540] -= dqdci;              /* dwdot[CH3]/d[HCO] */
    J[541] += dqdci;              /* dwdot[CH4]/d[HCO] */
    J[544] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[545] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[12];
    J[573] -= dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[574] += dqdci;              /* dwdot[CH4]/d[CH2O] */
    J[577] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[578] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[1068] -= dqdT;              /* dwdot[CH3]/dT */
    J[1069] += dqdT;              /* dwdot[CH4]/dT */
    J[1072] += dqdT;              /* dwdot[HCO]/dT */
    J[1073] -= dqdT;              /* dwdot[CH2O]/dT */

    /*reaction 164: CH3 + CH3OH <=> CH2OH + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[12]*sc[20];
    k_f = prefactor_units[163] * fwd_A[163]
                * exp(fwd_beta[163] * tc[0] - activation_units[163] * fwd_Ea[163] * invT);
    dlnkfdT = fwd_beta[163] * invT + activation_units[163] * fwd_Ea[163] * invT2;
    /* reverse */
    phi_r = sc[13]*sc[18];
    Kc = exp(g_RT[12] - g_RT[13] - g_RT[18] + g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[12] + h_RT[20]) + (h_RT[13] + h_RT[18]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[12] -= q; /* CH3 */
    wdot[13] += q; /* CH4 */
    wdot[18] += q; /* CH2OH */
    wdot[20] -= q; /* CH3OH */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[20];
    J[408] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[409] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[414] += dqdci;              /* dwdot[CH2OH]/d[CH3] */
    J[416] -= dqdci;              /* dwdot[CH3OH]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[18];
    J[441] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[442] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[447] += dqdci;              /* dwdot[CH2OH]/d[CH4] */
    J[449] -= dqdci;              /* dwdot[CH3OH]/d[CH4] */
    /* d()/d[CH2OH] */
    dqdci =  - k_r*sc[13];
    J[606] -= dqdci;              /* dwdot[CH3]/d[CH2OH] */
    J[607] += dqdci;              /* dwdot[CH4]/d[CH2OH] */
    J[612] += dqdci;              /* dwdot[CH2OH]/d[CH2OH] */
    J[614] -= dqdci;              /* dwdot[CH3OH]/d[CH2OH] */
    /* d()/d[CH3OH] */
    dqdci =  + k_f*sc[12];
    J[672] -= dqdci;              /* dwdot[CH3]/d[CH3OH] */
    J[673] += dqdci;              /* dwdot[CH4]/d[CH3OH] */
    J[678] += dqdci;              /* dwdot[CH2OH]/d[CH3OH] */
    J[680] -= dqdci;              /* dwdot[CH3OH]/d[CH3OH] */
    /* d()/dT */
    J[1068] -= dqdT;              /* dwdot[CH3]/dT */
    J[1069] += dqdT;              /* dwdot[CH4]/dT */
    J[1074] += dqdT;              /* dwdot[CH2OH]/dT */
    J[1076] -= dqdT;              /* dwdot[CH3OH]/dT */

    /*reaction 165: CH3 + CH3OH <=> CH3O + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[12]*sc[20];
    k_f = prefactor_units[164] * fwd_A[164]
                * exp(fwd_beta[164] * tc[0] - activation_units[164] * fwd_Ea[164] * invT);
    dlnkfdT = fwd_beta[164] * invT + activation_units[164] * fwd_Ea[164] * invT2;
    /* reverse */
    phi_r = sc[13]*sc[19];
    Kc = exp(g_RT[12] - g_RT[13] - g_RT[19] + g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[12] + h_RT[20]) + (h_RT[13] + h_RT[19]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[12] -= q; /* CH3 */
    wdot[13] += q; /* CH4 */
    wdot[19] += q; /* CH3O */
    wdot[20] -= q; /* CH3OH */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[20];
    J[408] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[409] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[415] += dqdci;              /* dwdot[CH3O]/d[CH3] */
    J[416] -= dqdci;              /* dwdot[CH3OH]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[19];
    J[441] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[442] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[448] += dqdci;              /* dwdot[CH3O]/d[CH4] */
    J[449] -= dqdci;              /* dwdot[CH3OH]/d[CH4] */
    /* d()/d[CH3O] */
    dqdci =  - k_r*sc[13];
    J[639] -= dqdci;              /* dwdot[CH3]/d[CH3O] */
    J[640] += dqdci;              /* dwdot[CH4]/d[CH3O] */
    J[646] += dqdci;              /* dwdot[CH3O]/d[CH3O] */
    J[647] -= dqdci;              /* dwdot[CH3OH]/d[CH3O] */
    /* d()/d[CH3OH] */
    dqdci =  + k_f*sc[12];
    J[672] -= dqdci;              /* dwdot[CH3]/d[CH3OH] */
    J[673] += dqdci;              /* dwdot[CH4]/d[CH3OH] */
    J[679] += dqdci;              /* dwdot[CH3O]/d[CH3OH] */
    J[680] -= dqdci;              /* dwdot[CH3OH]/d[CH3OH] */
    /* d()/dT */
    J[1068] -= dqdT;              /* dwdot[CH3]/dT */
    J[1069] += dqdT;              /* dwdot[CH4]/dT */
    J[1075] += dqdT;              /* dwdot[CH3O]/dT */
    J[1076] -= dqdT;              /* dwdot[CH3OH]/dT */

    /*reaction 166: CH3 + C2H4 <=> C2H3 + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[12]*sc[24];
    k_f = prefactor_units[165] * fwd_A[165]
                * exp(fwd_beta[165] * tc[0] - activation_units[165] * fwd_Ea[165] * invT);
    dlnkfdT = fwd_beta[165] * invT + activation_units[165] * fwd_Ea[165] * invT2;
    /* reverse */
    phi_r = sc[13]*sc[23];
    Kc = exp(g_RT[12] - g_RT[13] - g_RT[23] + g_RT[24]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[12] + h_RT[24]) + (h_RT[13] + h_RT[23]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[12] -= q; /* CH3 */
    wdot[13] += q; /* CH4 */
    wdot[23] += q; /* C2H3 */
    wdot[24] -= q; /* C2H4 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[24];
    J[408] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[409] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[419] += dqdci;              /* dwdot[C2H3]/d[CH3] */
    J[420] -= dqdci;              /* dwdot[C2H4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[23];
    J[441] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[442] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[452] += dqdci;              /* dwdot[C2H3]/d[CH4] */
    J[453] -= dqdci;              /* dwdot[C2H4]/d[CH4] */
    /* d()/d[C2H3] */
    dqdci =  - k_r*sc[13];
    J[771] -= dqdci;              /* dwdot[CH3]/d[C2H3] */
    J[772] += dqdci;              /* dwdot[CH4]/d[C2H3] */
    J[782] += dqdci;              /* dwdot[C2H3]/d[C2H3] */
    J[783] -= dqdci;              /* dwdot[C2H4]/d[C2H3] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[12];
    J[804] -= dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[805] += dqdci;              /* dwdot[CH4]/d[C2H4] */
    J[815] += dqdci;              /* dwdot[C2H3]/d[C2H4] */
    J[816] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    /* d()/dT */
    J[1068] -= dqdT;              /* dwdot[CH3]/dT */
    J[1069] += dqdT;              /* dwdot[CH4]/dT */
    J[1079] += dqdT;              /* dwdot[C2H3]/dT */
    J[1080] -= dqdT;              /* dwdot[C2H4]/dT */

    /*reaction 167: CH3 + C2H6 <=> C2H5 + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[12]*sc[26];
    k_f = prefactor_units[166] * fwd_A[166]
                * exp(fwd_beta[166] * tc[0] - activation_units[166] * fwd_Ea[166] * invT);
    dlnkfdT = fwd_beta[166] * invT + activation_units[166] * fwd_Ea[166] * invT2;
    /* reverse */
    phi_r = sc[13]*sc[25];
    Kc = exp(g_RT[12] - g_RT[13] - g_RT[25] + g_RT[26]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[12] + h_RT[26]) + (h_RT[13] + h_RT[25]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[12] -= q; /* CH3 */
    wdot[13] += q; /* CH4 */
    wdot[25] += q; /* C2H5 */
    wdot[26] -= q; /* C2H6 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[26];
    J[408] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[409] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[421] += dqdci;              /* dwdot[C2H5]/d[CH3] */
    J[422] -= dqdci;              /* dwdot[C2H6]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[25];
    J[441] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[442] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[454] += dqdci;              /* dwdot[C2H5]/d[CH4] */
    J[455] -= dqdci;              /* dwdot[C2H6]/d[CH4] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[13];
    J[837] -= dqdci;              /* dwdot[CH3]/d[C2H5] */
    J[838] += dqdci;              /* dwdot[CH4]/d[C2H5] */
    J[850] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[851] -= dqdci;              /* dwdot[C2H6]/d[C2H5] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[12];
    J[870] -= dqdci;              /* dwdot[CH3]/d[C2H6] */
    J[871] += dqdci;              /* dwdot[CH4]/d[C2H6] */
    J[883] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[884] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[1068] -= dqdT;              /* dwdot[CH3]/dT */
    J[1069] += dqdT;              /* dwdot[CH4]/dT */
    J[1081] += dqdT;              /* dwdot[C2H5]/dT */
    J[1082] -= dqdT;              /* dwdot[C2H6]/dT */

    /*reaction 168: HCO + H2O <=> H + CO + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[16];
    k_f = prefactor_units[167] * fwd_A[167]
                * exp(fwd_beta[167] * tc[0] - activation_units[167] * fwd_Ea[167] * invT);
    dlnkfdT = fwd_beta[167] * invT + activation_units[167] * fwd_Ea[167] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[5]*sc[14];
    Kc = refC * exp(-g_RT[1] + g_RT[5] - g_RT[5] - g_RT[14] + g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[16]) + (h_RT[1] + h_RT[5] + h_RT[14]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[14] += q; /* CO */
    wdot[16] -= q; /* HCO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[5]*sc[14];
    J[34] += dqdci;               /* dwdot[H]/d[H] */
    J[47] += dqdci;               /* dwdot[CO]/d[H] */
    J[49] -= dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[16] - k_r*sc[1]*sc[14];
    J[166] += dqdci;              /* dwdot[H]/d[H2O] */
    J[179] += dqdci;              /* dwdot[CO]/d[H2O] */
    J[181] -= dqdci;              /* dwdot[HCO]/d[H2O] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[1]*sc[5];
    J[463] += dqdci;              /* dwdot[H]/d[CO] */
    J[476] += dqdci;              /* dwdot[CO]/d[CO] */
    J[478] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[5];
    J[529] += dqdci;              /* dwdot[H]/d[HCO] */
    J[542] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[544] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[1057] += dqdT;              /* dwdot[H]/dT */
    J[1070] += dqdT;              /* dwdot[CO]/dT */
    J[1072] -= dqdT;              /* dwdot[HCO]/dT */

    /*reaction 169: HCO + O2 <=> HO2 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[16];
    k_f = prefactor_units[168] * fwd_A[168]
                * exp(fwd_beta[168] * tc[0] - activation_units[168] * fwd_Ea[168] * invT);
    dlnkfdT = fwd_beta[168] * invT + activation_units[168] * fwd_Ea[168] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[14];
    Kc = exp(g_RT[3] - g_RT[6] - g_RT[14] + g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[16]) + (h_RT[6] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    wdot[14] += q; /* CO */
    wdot[16] -= q; /* HCO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[16];
    J[102] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[105] += dqdci;              /* dwdot[HO2]/d[O2] */
    J[113] += dqdci;              /* dwdot[CO]/d[O2] */
    J[115] -= dqdci;              /* dwdot[HCO]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[14];
    J[201] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[204] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[212] += dqdci;              /* dwdot[CO]/d[HO2] */
    J[214] -= dqdci;              /* dwdot[HCO]/d[HO2] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[6];
    J[465] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[468] += dqdci;              /* dwdot[HO2]/d[CO] */
    J[476] += dqdci;              /* dwdot[CO]/d[CO] */
    J[478] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[3];
    J[531] -= dqdci;              /* dwdot[O2]/d[HCO] */
    J[534] += dqdci;              /* dwdot[HO2]/d[HCO] */
    J[542] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[544] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[1059] -= dqdT;              /* dwdot[O2]/dT */
    J[1062] += dqdT;              /* dwdot[HO2]/dT */
    J[1070] += dqdT;              /* dwdot[CO]/dT */
    J[1072] -= dqdT;              /* dwdot[HCO]/dT */

    /*reaction 170: CH2OH + O2 <=> HO2 + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[18];
    k_f = prefactor_units[169] * fwd_A[169]
                * exp(fwd_beta[169] * tc[0] - activation_units[169] * fwd_Ea[169] * invT);
    dlnkfdT = fwd_beta[169] * invT + activation_units[169] * fwd_Ea[169] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[17];
    Kc = exp(g_RT[3] - g_RT[6] - g_RT[17] + g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[18]) + (h_RT[6] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    wdot[17] += q; /* CH2O */
    wdot[18] -= q; /* CH2OH */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[18];
    J[102] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[105] += dqdci;              /* dwdot[HO2]/d[O2] */
    J[116] += dqdci;              /* dwdot[CH2O]/d[O2] */
    J[117] -= dqdci;              /* dwdot[CH2OH]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[17];
    J[201] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[204] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[215] += dqdci;              /* dwdot[CH2O]/d[HO2] */
    J[216] -= dqdci;              /* dwdot[CH2OH]/d[HO2] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[6];
    J[564] -= dqdci;              /* dwdot[O2]/d[CH2O] */
    J[567] += dqdci;              /* dwdot[HO2]/d[CH2O] */
    J[578] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[579] -= dqdci;              /* dwdot[CH2OH]/d[CH2O] */
    /* d()/d[CH2OH] */
    dqdci =  + k_f*sc[3];
    J[597] -= dqdci;              /* dwdot[O2]/d[CH2OH] */
    J[600] += dqdci;              /* dwdot[HO2]/d[CH2OH] */
    J[611] += dqdci;              /* dwdot[CH2O]/d[CH2OH] */
    J[612] -= dqdci;              /* dwdot[CH2OH]/d[CH2OH] */
    /* d()/dT */
    J[1059] -= dqdT;              /* dwdot[O2]/dT */
    J[1062] += dqdT;              /* dwdot[HO2]/dT */
    J[1073] += dqdT;              /* dwdot[CH2O]/dT */
    J[1074] -= dqdT;              /* dwdot[CH2OH]/dT */

    /*reaction 171: CH3O + O2 <=> HO2 + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[19];
    k_f = prefactor_units[170] * fwd_A[170]
                * exp(fwd_beta[170] * tc[0] - activation_units[170] * fwd_Ea[170] * invT);
    dlnkfdT = fwd_beta[170] * invT + activation_units[170] * fwd_Ea[170] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[17];
    Kc = exp(g_RT[3] - g_RT[6] - g_RT[17] + g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[19]) + (h_RT[6] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    wdot[17] += q; /* CH2O */
    wdot[19] -= q; /* CH3O */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[19];
    J[102] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[105] += dqdci;              /* dwdot[HO2]/d[O2] */
    J[116] += dqdci;              /* dwdot[CH2O]/d[O2] */
    J[118] -= dqdci;              /* dwdot[CH3O]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[17];
    J[201] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[204] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[215] += dqdci;              /* dwdot[CH2O]/d[HO2] */
    J[217] -= dqdci;              /* dwdot[CH3O]/d[HO2] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[6];
    J[564] -= dqdci;              /* dwdot[O2]/d[CH2O] */
    J[567] += dqdci;              /* dwdot[HO2]/d[CH2O] */
    J[578] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[580] -= dqdci;              /* dwdot[CH3O]/d[CH2O] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[3];
    J[630] -= dqdci;              /* dwdot[O2]/d[CH3O] */
    J[633] += dqdci;              /* dwdot[HO2]/d[CH3O] */
    J[644] += dqdci;              /* dwdot[CH2O]/d[CH3O] */
    J[646] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[1059] -= dqdT;              /* dwdot[O2]/dT */
    J[1062] += dqdT;              /* dwdot[HO2]/dT */
    J[1073] += dqdT;              /* dwdot[CH2O]/dT */
    J[1075] -= dqdT;              /* dwdot[CH3O]/dT */

    /*reaction 172: C2H + O2 <=> HCO + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[21];
    k_f = prefactor_units[171] * fwd_A[171]
                * exp(fwd_beta[171] * tc[0] - activation_units[171] * fwd_Ea[171] * invT);
    dlnkfdT = fwd_beta[171] * invT + activation_units[171] * fwd_Ea[171] * invT2;
    /* reverse */
    phi_r = sc[14]*sc[16];
    Kc = exp(g_RT[3] - g_RT[14] - g_RT[16] + g_RT[21]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[21]) + (h_RT[14] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[14] += q; /* CO */
    wdot[16] += q; /* HCO */
    wdot[21] -= q; /* C2H */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[21];
    J[102] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[113] += dqdci;              /* dwdot[CO]/d[O2] */
    J[115] += dqdci;              /* dwdot[HCO]/d[O2] */
    J[120] -= dqdci;              /* dwdot[C2H]/d[O2] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[16];
    J[465] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[476] += dqdci;              /* dwdot[CO]/d[CO] */
    J[478] += dqdci;              /* dwdot[HCO]/d[CO] */
    J[483] -= dqdci;              /* dwdot[C2H]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[14];
    J[531] -= dqdci;              /* dwdot[O2]/d[HCO] */
    J[542] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[544] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[549] -= dqdci;              /* dwdot[C2H]/d[HCO] */
    /* d()/d[C2H] */
    dqdci =  + k_f*sc[3];
    J[696] -= dqdci;              /* dwdot[O2]/d[C2H] */
    J[707] += dqdci;              /* dwdot[CO]/d[C2H] */
    J[709] += dqdci;              /* dwdot[HCO]/d[C2H] */
    J[714] -= dqdci;              /* dwdot[C2H]/d[C2H] */
    /* d()/dT */
    J[1059] -= dqdT;              /* dwdot[O2]/dT */
    J[1070] += dqdT;              /* dwdot[CO]/dT */
    J[1072] += dqdT;              /* dwdot[HCO]/dT */
    J[1077] -= dqdT;              /* dwdot[C2H]/dT */

    /*reaction 173: C2H + H2 <=> H + C2H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[21];
    k_f = prefactor_units[172] * fwd_A[172]
                * exp(fwd_beta[172] * tc[0] - activation_units[172] * fwd_Ea[172] * invT);
    dlnkfdT = fwd_beta[172] * invT + activation_units[172] * fwd_Ea[172] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[22];
    Kc = exp(g_RT[0] - g_RT[1] + g_RT[21] - g_RT[22]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[21]) + (h_RT[1] + h_RT[22]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[1] += q; /* H */
    wdot[21] -= q; /* C2H */
    wdot[22] += q; /* C2H2 */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[21];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += dqdci;                /* dwdot[H]/d[H2] */
    J[21] -= dqdci;               /* dwdot[C2H]/d[H2] */
    J[22] += dqdci;               /* dwdot[C2H2]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[22];
    J[33] -= dqdci;               /* dwdot[H2]/d[H] */
    J[34] += dqdci;               /* dwdot[H]/d[H] */
    J[54] -= dqdci;               /* dwdot[C2H]/d[H] */
    J[55] += dqdci;               /* dwdot[C2H2]/d[H] */
    /* d()/d[C2H] */
    dqdci =  + k_f*sc[0];
    J[693] -= dqdci;              /* dwdot[H2]/d[C2H] */
    J[694] += dqdci;              /* dwdot[H]/d[C2H] */
    J[714] -= dqdci;              /* dwdot[C2H]/d[C2H] */
    J[715] += dqdci;              /* dwdot[C2H2]/d[C2H] */
    /* d()/d[C2H2] */
    dqdci =  - k_r*sc[1];
    J[726] -= dqdci;              /* dwdot[H2]/d[C2H2] */
    J[727] += dqdci;              /* dwdot[H]/d[C2H2] */
    J[747] -= dqdci;              /* dwdot[C2H]/d[C2H2] */
    J[748] += dqdci;              /* dwdot[C2H2]/d[C2H2] */
    /* d()/dT */
    J[1056] -= dqdT;              /* dwdot[H2]/dT */
    J[1057] += dqdT;              /* dwdot[H]/dT */
    J[1077] -= dqdT;              /* dwdot[C2H]/dT */
    J[1078] += dqdT;              /* dwdot[C2H2]/dT */

    /*reaction 174: C2H3 + O2 <=> HCO + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[23];
    k_f = prefactor_units[173] * fwd_A[173]
                * exp(fwd_beta[173] * tc[0] - activation_units[173] * fwd_Ea[173] * invT);
    dlnkfdT = fwd_beta[173] * invT + activation_units[173] * fwd_Ea[173] * invT2;
    /* reverse */
    phi_r = sc[16]*sc[17];
    Kc = exp(g_RT[3] - g_RT[16] - g_RT[17] + g_RT[23]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[23]) + (h_RT[16] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[16] += q; /* HCO */
    wdot[17] += q; /* CH2O */
    wdot[23] -= q; /* C2H3 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[23];
    J[102] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[115] += dqdci;              /* dwdot[HCO]/d[O2] */
    J[116] += dqdci;              /* dwdot[CH2O]/d[O2] */
    J[122] -= dqdci;              /* dwdot[C2H3]/d[O2] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[17];
    J[531] -= dqdci;              /* dwdot[O2]/d[HCO] */
    J[544] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[545] += dqdci;              /* dwdot[CH2O]/d[HCO] */
    J[551] -= dqdci;              /* dwdot[C2H3]/d[HCO] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[16];
    J[564] -= dqdci;              /* dwdot[O2]/d[CH2O] */
    J[577] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[578] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[584] -= dqdci;              /* dwdot[C2H3]/d[CH2O] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[3];
    J[762] -= dqdci;              /* dwdot[O2]/d[C2H3] */
    J[775] += dqdci;              /* dwdot[HCO]/d[C2H3] */
    J[776] += dqdci;              /* dwdot[CH2O]/d[C2H3] */
    J[782] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    /* d()/dT */
    J[1059] -= dqdT;              /* dwdot[O2]/dT */
    J[1072] += dqdT;              /* dwdot[HCO]/dT */
    J[1073] += dqdT;              /* dwdot[CH2O]/dT */
    J[1079] -= dqdT;              /* dwdot[C2H3]/dT */

    /*reaction 175: C2H5 + O2 <=> HO2 + C2H4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[25];
    k_f = prefactor_units[174] * fwd_A[174]
                * exp(fwd_beta[174] * tc[0] - activation_units[174] * fwd_Ea[174] * invT);
    dlnkfdT = fwd_beta[174] * invT + activation_units[174] * fwd_Ea[174] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[24];
    Kc = exp(g_RT[3] - g_RT[6] - g_RT[24] + g_RT[25]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[25]) + (h_RT[6] + h_RT[24]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    wdot[24] += q; /* C2H4 */
    wdot[25] -= q; /* C2H5 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[25];
    J[102] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[105] += dqdci;              /* dwdot[HO2]/d[O2] */
    J[123] += dqdci;              /* dwdot[C2H4]/d[O2] */
    J[124] -= dqdci;              /* dwdot[C2H5]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[24];
    J[201] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[204] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[222] += dqdci;              /* dwdot[C2H4]/d[HO2] */
    J[223] -= dqdci;              /* dwdot[C2H5]/d[HO2] */
    /* d()/d[C2H4] */
    dqdci =  - k_r*sc[6];
    J[795] -= dqdci;              /* dwdot[O2]/d[C2H4] */
    J[798] += dqdci;              /* dwdot[HO2]/d[C2H4] */
    J[816] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[817] -= dqdci;              /* dwdot[C2H5]/d[C2H4] */
    /* d()/d[C2H5] */
    dqdci =  + k_f*sc[3];
    J[828] -= dqdci;              /* dwdot[O2]/d[C2H5] */
    J[831] += dqdci;              /* dwdot[HO2]/d[C2H5] */
    J[849] += dqdci;              /* dwdot[C2H4]/d[C2H5] */
    J[850] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
    /* d()/dT */
    J[1059] -= dqdT;              /* dwdot[O2]/dT */
    J[1062] += dqdT;              /* dwdot[HO2]/dT */
    J[1080] += dqdT;              /* dwdot[C2H4]/dT */
    J[1081] -= dqdT;              /* dwdot[C2H5]/dT */

    /*reaction 176: HCCO + O2 <=> OH + 2.000000 CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[27];
    k_f = prefactor_units[175] * fwd_A[175]
                * exp(fwd_beta[175] * tc[0] - activation_units[175] * fwd_Ea[175] * invT);
    dlnkfdT = fwd_beta[175] * invT + activation_units[175] * fwd_Ea[175] * invT2;
    /* reverse */
    phi_r = sc[4]*pow(sc[14], 2.000000);
    Kc = refC * exp(g_RT[3] - g_RT[4] - 2.000000*g_RT[14] + g_RT[27]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[27]) + (h_RT[4] + 2.000000*h_RT[14]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O2 */
    wdot[4] += q; /* OH */
    wdot[14] += 2 * q; /* CO */
    wdot[27] -= q; /* HCCO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[27];
    J[102] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[103] += dqdci;              /* dwdot[OH]/d[O2] */
    J[113] += 2 * dqdci;          /* dwdot[CO]/d[O2] */
    J[126] -= dqdci;              /* dwdot[HCCO]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*pow(sc[14], 2.000000);
    J[135] -= dqdci;              /* dwdot[O2]/d[OH] */
    J[136] += dqdci;              /* dwdot[OH]/d[OH] */
    J[146] += 2 * dqdci;          /* dwdot[CO]/d[OH] */
    J[159] -= dqdci;              /* dwdot[HCCO]/d[OH] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[4]*2.000000*sc[14];
    J[465] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[466] += dqdci;              /* dwdot[OH]/d[CO] */
    J[476] += 2 * dqdci;          /* dwdot[CO]/d[CO] */
    J[489] -= dqdci;              /* dwdot[HCCO]/d[CO] */
    /* d()/d[HCCO] */
    dqdci =  + k_f*sc[3];
    J[894] -= dqdci;              /* dwdot[O2]/d[HCCO] */
    J[895] += dqdci;              /* dwdot[OH]/d[HCCO] */
    J[905] += 2 * dqdci;          /* dwdot[CO]/d[HCCO] */
    J[918] -= dqdci;              /* dwdot[HCCO]/d[HCCO] */
    /* d()/dT */
    J[1059] -= dqdT;              /* dwdot[O2]/dT */
    J[1060] += dqdT;              /* dwdot[OH]/dT */
    J[1070] += 2 * dqdT;          /* dwdot[CO]/dT */
    J[1083] -= dqdT;              /* dwdot[HCCO]/dT */

    /*reaction 177: 2.000000 HCCO <=> 2.000000 CO + C2H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[27], 2.000000);
    k_f = prefactor_units[176] * fwd_A[176]
                * exp(fwd_beta[176] * tc[0] - activation_units[176] * fwd_Ea[176] * invT);
    dlnkfdT = fwd_beta[176] * invT + activation_units[176] * fwd_Ea[176] * invT2;
    /* reverse */
    phi_r = pow(sc[14], 2.000000)*sc[22];
    Kc = refC * exp(-2.000000*g_RT[14] - g_RT[22] + 2.000000*g_RT[27]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[27]) + (2.000000*h_RT[14] + h_RT[22]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[14] += 2 * q; /* CO */
    wdot[22] += q; /* C2H2 */
    wdot[27] -= 2 * q; /* HCCO */
    /* d()/d[CO] */
    dqdci =  - k_r*2.000000*sc[14]*sc[22];
    J[476] += 2 * dqdci;          /* dwdot[CO]/d[CO] */
    J[484] += dqdci;              /* dwdot[C2H2]/d[CO] */
    J[489] += -2 * dqdci;         /* dwdot[HCCO]/d[CO] */
    /* d()/d[C2H2] */
    dqdci =  - k_r*pow(sc[14], 2.000000);
    J[740] += 2 * dqdci;          /* dwdot[CO]/d[C2H2] */
    J[748] += dqdci;              /* dwdot[C2H2]/d[C2H2] */
    J[753] += -2 * dqdci;         /* dwdot[HCCO]/d[C2H2] */
    /* d()/d[HCCO] */
    dqdci =  + k_f*2.000000*sc[27];
    J[905] += 2 * dqdci;          /* dwdot[CO]/d[HCCO] */
    J[913] += dqdci;              /* dwdot[C2H2]/d[HCCO] */
    J[918] += -2 * dqdci;         /* dwdot[HCCO]/d[HCCO] */
    /* d()/dT */
    J[1070] += 2 * dqdT;          /* dwdot[CO]/dT */
    J[1078] += dqdT;              /* dwdot[C2H2]/dT */
    J[1083] += -2 * dqdT;         /* dwdot[HCCO]/dT */

    amrex::Real c_R[32], dcRdT[32], e_RT[32];
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
    for (int k = 0; k < 32; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[1056+k];
    }

    amrex::Real cmixinv = 1.0/cmix;
    amrex::Real tmp1 = ehmix*cmixinv;
    amrex::Real tmp3 = cmixinv*T;
    amrex::Real tmp2 = tmp1*tmp3;
    amrex::Real dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 32; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 32; ++m) {
            dehmixdc += eh_RT[m]*J[k*33+m];
        }
        J[k*33+32] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[1088] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;
}

#endif
/* Transport function declarations  */


void egtransetLENIMC(int* LENIMC ) {
    *LENIMC = 130;}


void egtransetLENRMC(int* LENRMC ) {
    *LENRMC = 20576;}


void egtransetNO(int* NO ) {
    *NO = 4;}


void egtransetKK(int* KK ) {
    *KK = 32;}


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
    WT[8] = 1.20111500E+01;
    WT[9] = 1.30191200E+01;
    WT[10] = 1.40270900E+01;
    WT[11] = 1.40270900E+01;
    WT[12] = 1.50350600E+01;
    WT[13] = 1.60430300E+01;
    WT[14] = 2.80105500E+01;
    WT[15] = 4.40099500E+01;
    WT[16] = 2.90185200E+01;
    WT[17] = 3.00264900E+01;
    WT[18] = 3.10344600E+01;
    WT[19] = 3.10344600E+01;
    WT[20] = 3.20424300E+01;
    WT[21] = 2.50302700E+01;
    WT[22] = 2.60382400E+01;
    WT[23] = 2.70462100E+01;
    WT[24] = 2.80541800E+01;
    WT[25] = 2.90621500E+01;
    WT[26] = 3.00701200E+01;
    WT[27] = 4.10296700E+01;
    WT[28] = 4.20376400E+01;
    WT[29] = 4.20376400E+01;
    WT[30] = 2.80134000E+01;
    WT[31] = 3.99480000E+01;
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
    EPS[8] = 7.14000000E+01;
    EPS[9] = 8.00000000E+01;
    EPS[10] = 1.44000000E+02;
    EPS[11] = 1.44000000E+02;
    EPS[12] = 1.44000000E+02;
    EPS[13] = 1.41400000E+02;
    EPS[14] = 9.81000000E+01;
    EPS[15] = 2.44000000E+02;
    EPS[16] = 4.98000000E+02;
    EPS[17] = 4.98000000E+02;
    EPS[18] = 4.17000000E+02;
    EPS[19] = 4.17000000E+02;
    EPS[20] = 4.81800000E+02;
    EPS[21] = 2.09000000E+02;
    EPS[22] = 2.09000000E+02;
    EPS[23] = 2.09000000E+02;
    EPS[24] = 2.80800000E+02;
    EPS[25] = 2.52300000E+02;
    EPS[26] = 2.52300000E+02;
    EPS[27] = 1.50000000E+02;
    EPS[28] = 4.36000000E+02;
    EPS[29] = 4.36000000E+02;
    EPS[30] = 9.75300000E+01;
    EPS[31] = 1.36500000E+02;
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
    SIG[8] = 3.29800000E+00;
    SIG[9] = 2.75000000E+00;
    SIG[10] = 3.80000000E+00;
    SIG[11] = 3.80000000E+00;
    SIG[12] = 3.80000000E+00;
    SIG[13] = 3.74600000E+00;
    SIG[14] = 3.65000000E+00;
    SIG[15] = 3.76300000E+00;
    SIG[16] = 3.59000000E+00;
    SIG[17] = 3.59000000E+00;
    SIG[18] = 3.69000000E+00;
    SIG[19] = 3.69000000E+00;
    SIG[20] = 3.62600000E+00;
    SIG[21] = 4.10000000E+00;
    SIG[22] = 4.10000000E+00;
    SIG[23] = 4.10000000E+00;
    SIG[24] = 3.97100000E+00;
    SIG[25] = 4.30200000E+00;
    SIG[26] = 4.30200000E+00;
    SIG[27] = 2.50000000E+00;
    SIG[28] = 3.97000000E+00;
    SIG[29] = 3.97000000E+00;
    SIG[30] = 3.62100000E+00;
    SIG[31] = 3.33000000E+00;
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
    DIP[15] = 0.00000000E+00;
    DIP[16] = 0.00000000E+00;
    DIP[17] = 0.00000000E+00;
    DIP[18] = 1.70000000E+00;
    DIP[19] = 1.70000000E+00;
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
    POL[10] = 0.00000000E+00;
    POL[11] = 0.00000000E+00;
    POL[12] = 0.00000000E+00;
    POL[13] = 2.60000000E+00;
    POL[14] = 1.95000000E+00;
    POL[15] = 2.65000000E+00;
    POL[16] = 0.00000000E+00;
    POL[17] = 0.00000000E+00;
    POL[18] = 0.00000000E+00;
    POL[19] = 0.00000000E+00;
    POL[20] = 0.00000000E+00;
    POL[21] = 0.00000000E+00;
    POL[22] = 0.00000000E+00;
    POL[23] = 0.00000000E+00;
    POL[24] = 0.00000000E+00;
    POL[25] = 0.00000000E+00;
    POL[26] = 0.00000000E+00;
    POL[27] = 0.00000000E+00;
    POL[28] = 0.00000000E+00;
    POL[29] = 0.00000000E+00;
    POL[30] = 1.76000000E+00;
    POL[31] = 0.00000000E+00;
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
    ZROT[8] = 0.00000000E+00;
    ZROT[9] = 0.00000000E+00;
    ZROT[10] = 0.00000000E+00;
    ZROT[11] = 0.00000000E+00;
    ZROT[12] = 0.00000000E+00;
    ZROT[13] = 1.30000000E+01;
    ZROT[14] = 1.80000000E+00;
    ZROT[15] = 2.10000000E+00;
    ZROT[16] = 0.00000000E+00;
    ZROT[17] = 2.00000000E+00;
    ZROT[18] = 2.00000000E+00;
    ZROT[19] = 2.00000000E+00;
    ZROT[20] = 1.00000000E+00;
    ZROT[21] = 2.50000000E+00;
    ZROT[22] = 2.50000000E+00;
    ZROT[23] = 1.00000000E+00;
    ZROT[24] = 1.50000000E+00;
    ZROT[25] = 1.50000000E+00;
    ZROT[26] = 1.50000000E+00;
    ZROT[27] = 1.00000000E+00;
    ZROT[28] = 2.00000000E+00;
    ZROT[29] = 2.00000000E+00;
    ZROT[30] = 4.00000000E+00;
    ZROT[31] = 0.00000000E+00;
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
    NLIN[8] = 0;
    NLIN[9] = 1;
    NLIN[10] = 1;
    NLIN[11] = 1;
    NLIN[12] = 1;
    NLIN[13] = 2;
    NLIN[14] = 1;
    NLIN[15] = 1;
    NLIN[16] = 2;
    NLIN[17] = 2;
    NLIN[18] = 2;
    NLIN[19] = 2;
    NLIN[20] = 2;
    NLIN[21] = 1;
    NLIN[22] = 1;
    NLIN[23] = 2;
    NLIN[24] = 2;
    NLIN[25] = 2;
    NLIN[26] = 2;
    NLIN[27] = 2;
    NLIN[28] = 2;
    NLIN[29] = 2;
    NLIN[30] = 1;
    NLIN[31] = 0;
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
    COFETA[28] = -1.71312832E+01;
    COFETA[29] = 2.68036374E+00;
    COFETA[30] = -2.72570227E-01;
    COFETA[31] = 1.21650964E-02;
    COFETA[32] = -1.50067648E+01;
    COFETA[33] = 1.69625877E+00;
    COFETA[34] = -1.42936462E-01;
    COFETA[35] = 6.47223426E-03;
    COFETA[36] = -1.51956901E+01;
    COFETA[37] = 1.92606504E+00;
    COFETA[38] = -1.73487476E-01;
    COFETA[39] = 7.82572931E-03;
    COFETA[40] = -2.02663469E+01;
    COFETA[41] = 3.63241793E+00;
    COFETA[42] = -3.95581049E-01;
    COFETA[43] = 1.74725495E-02;
    COFETA[44] = -2.02663469E+01;
    COFETA[45] = 3.63241793E+00;
    COFETA[46] = -3.95581049E-01;
    COFETA[47] = 1.74725495E-02;
    COFETA[48] = -2.02316497E+01;
    COFETA[49] = 3.63241793E+00;
    COFETA[50] = -3.95581049E-01;
    COFETA[51] = 1.74725495E-02;
    COFETA[52] = -2.00094664E+01;
    COFETA[53] = 3.57220167E+00;
    COFETA[54] = -3.87936446E-01;
    COFETA[55] = 1.71483254E-02;
    COFETA[56] = -1.66188336E+01;
    COFETA[57] = 2.40307799E+00;
    COFETA[58] = -2.36167638E-01;
    COFETA[59] = 1.05714061E-02;
    COFETA[60] = -2.40014975E+01;
    COFETA[61] = 5.14359547E+00;
    COFETA[62] = -5.74269731E-01;
    COFETA[63] = 2.44937679E-02;
    COFETA[64] = -1.98501306E+01;
    COFETA[65] = 2.69480162E+00;
    COFETA[66] = -1.65880845E-01;
    COFETA[67] = 3.14504769E-03;
    COFETA[68] = -1.98330577E+01;
    COFETA[69] = 2.69480162E+00;
    COFETA[70] = -1.65880845E-01;
    COFETA[71] = 3.14504769E-03;
    COFETA[72] = -1.99945919E+01;
    COFETA[73] = 2.86923313E+00;
    COFETA[74] = -2.03325661E-01;
    COFETA[75] = 5.39056989E-03;
    COFETA[76] = -1.99945919E+01;
    COFETA[77] = 2.86923313E+00;
    COFETA[78] = -2.03325661E-01;
    COFETA[79] = 5.39056989E-03;
    COFETA[80] = -2.05644525E+01;
    COFETA[81] = 3.03946431E+00;
    COFETA[82] = -2.16994867E-01;
    COFETA[83] = 5.61394012E-03;
    COFETA[84] = -2.33863848E+01;
    COFETA[85] = 4.80350223E+00;
    COFETA[86] = -5.38341336E-01;
    COFETA[87] = 2.32747213E-02;
    COFETA[88] = -2.33666446E+01;
    COFETA[89] = 4.80350223E+00;
    COFETA[90] = -5.38341336E-01;
    COFETA[91] = 2.32747213E-02;
    COFETA[92] = -2.33476543E+01;
    COFETA[93] = 4.80350223E+00;
    COFETA[94] = -5.38341336E-01;
    COFETA[95] = 2.32747213E-02;
    COFETA[96] = -2.50655444E+01;
    COFETA[97] = 5.33982977E+00;
    COFETA[98] = -5.89982992E-01;
    COFETA[99] = 2.47780650E-02;
    COFETA[100] = -2.46581414E+01;
    COFETA[101] = 5.19497183E+00;
    COFETA[102] = -5.78827339E-01;
    COFETA[103] = 2.46050921E-02;
    COFETA[104] = -2.46410937E+01;
    COFETA[105] = 5.19497183E+00;
    COFETA[106] = -5.78827339E-01;
    COFETA[107] = 2.46050921E-02;
    COFETA[108] = -1.92183831E+01;
    COFETA[109] = 3.75164499E+00;
    COFETA[110] = -4.10390993E-01;
    COFETA[111] = 1.80861665E-02;
    COFETA[112] = -2.23395647E+01;
    COFETA[113] = 3.86433912E+00;
    COFETA[114] = -3.41553983E-01;
    COFETA[115] = 1.17083447E-02;
    COFETA[116] = -2.23395647E+01;
    COFETA[117] = 3.86433912E+00;
    COFETA[118] = -3.41553983E-01;
    COFETA[119] = 1.17083447E-02;
    COFETA[120] = -1.65695594E+01;
    COFETA[121] = 2.39056562E+00;
    COFETA[122] = -2.34558144E-01;
    COFETA[123] = 1.05024037E-02;
    COFETA[124] = -1.90422907E+01;
    COFETA[125] = 3.47025711E+00;
    COFETA[126] = -3.75102111E-01;
    COFETA[127] = 1.66086076E-02;
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
    COFLAM[28] = 8.83996545E-01;
    COFLAM[29] = 1.31525428E+00;
    COFLAM[30] = 1.91774576E-02;
    COFLAM[31] = -4.41642722E-03;
    COFLAM[32] = 2.06524872E+00;
    COFLAM[33] = 1.69625877E+00;
    COFLAM[34] = -1.42936462E-01;
    COFLAM[35] = 6.47223426E-03;
    COFLAM[36] = 2.08093824E+01;
    COFLAM[37] = -6.24180163E+00;
    COFLAM[38] = 9.82387427E-01;
    COFLAM[39] = -4.50360664E-02;
    COFLAM[40] = 1.29177902E+01;
    COFLAM[41] = -3.73745535E+00;
    COFLAM[42] = 7.15831021E-01;
    COFLAM[43] = -3.63846910E-02;
    COFLAM[44] = 1.89383266E+01;
    COFLAM[45] = -6.51018128E+00;
    COFLAM[46] = 1.13292060E+00;
    COFLAM[47] = -5.69603226E-02;
    COFLAM[48] = 1.39937901E+01;
    COFLAM[49] = -4.64256494E+00;
    COFLAM[50] = 9.07728674E-01;
    COFLAM[51] = -4.77274469E-02;
    COFLAM[52] = 1.33091602E+01;
    COFLAM[53] = -4.96140261E+00;
    COFLAM[54] = 1.03295505E+00;
    COFLAM[55] = -5.63420090E-02;
    COFLAM[56] = 1.18777264E+01;
    COFLAM[57] = -3.15463949E+00;
    COFLAM[58] = 6.01973268E-01;
    COFLAM[59] = -3.03211261E-02;
    COFLAM[60] = -1.13649314E+01;
    COFLAM[61] = 5.88177395E+00;
    COFLAM[62] = -5.68651819E-01;
    COFLAM[63] = 2.03561485E-02;
    COFLAM[64] = 6.30243508E+00;
    COFLAM[65] = -2.22810801E+00;
    COFLAM[66] = 6.37340514E-01;
    COFLAM[67] = -3.81056018E-02;
    COFLAM[68] = 5.39305086E+00;
    COFLAM[69] = -2.39312375E+00;
    COFLAM[70] = 7.39585221E-01;
    COFLAM[71] = -4.58435589E-02;
    COFLAM[72] = 1.05517076E+00;
    COFLAM[73] = 4.29805638E-02;
    COFLAM[74] = 3.31460301E-01;
    COFLAM[75] = -2.40803935E-02;
    COFLAM[76] = -6.14588187E+00;
    COFLAM[77] = 2.47428873E+00;
    COFLAM[78] = 6.43999571E-02;
    COFLAM[79] = -1.45368336E-02;
    COFLAM[80] = -1.83491607E+00;
    COFLAM[81] = 5.55823162E-01;
    COFLAM[82] = 3.51309919E-01;
    COFLAM[83] = -2.85527489E-02;
    COFLAM[84] = 1.33513522E+01;
    COFLAM[85] = -4.32126196E+00;
    COFLAM[86] = 8.36087620E-01;
    COFLAM[87] = -4.37693086E-02;
    COFLAM[88] = -7.70164502E+00;
    COFLAM[89] = 4.56884453E+00;
    COFLAM[90] = -4.04747583E-01;
    COFLAM[91] = 1.40841060E-02;
    COFLAM[92] = -9.10384516E+00;
    COFLAM[93] = 4.54798869E+00;
    COFLAM[94] = -3.18114458E-01;
    COFLAM[95] = 6.59577576E-03;
    COFLAM[96] = -1.46152839E+01;
    COFLAM[97] = 6.36251406E+00;
    COFLAM[98] = -5.03832130E-01;
    COFLAM[99] = 1.26121050E-02;
    COFLAM[100] = -8.95009705E+00;
    COFLAM[101] = 4.02515080E+00;
    COFLAM[102] = -1.84063946E-01;
    COFLAM[103] = -1.94054752E-03;
    COFLAM[104] = -1.09902209E+01;
    COFLAM[105] = 4.70647707E+00;
    COFLAM[106] = -2.52272495E-01;
    COFLAM[107] = 1.75193258E-04;
    COFLAM[108] = -5.96382813E+00;
    COFLAM[109] = 4.39356020E+00;
    COFLAM[110] = -3.97598260E-01;
    COFLAM[111] = 1.39749624E-02;
    COFLAM[112] = -8.32871231E+00;
    COFLAM[113] = 3.97067262E+00;
    COFLAM[114] = -2.21252287E-01;
    COFLAM[115] = 1.47870329E-03;
    COFLAM[116] = -5.36683560E+00;
    COFLAM[117] = 2.98311144E+00;
    COFLAM[118] = -1.14643230E-01;
    COFLAM[119] = -2.12830349E-03;
    COFLAM[120] = 1.29306158E+01;
    COFLAM[121] = -3.52817362E+00;
    COFLAM[122] = 6.45499013E-01;
    COFLAM[123] = -3.19375299E-02;
    COFLAM[124] = -3.17202048E+00;
    COFLAM[125] = 3.47025711E+00;
    COFLAM[126] = -3.75102111E-01;
    COFLAM[127] = 1.66086076E-02;
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
    COFD[28] = -1.19006548E+01;
    COFD[29] = 2.57507000E+00;
    COFD[30] = -1.24033737E-01;
    COFD[31] = 5.56694959E-03;
    COFD[32] = -1.08369481E+01;
    COFD[33] = 2.19094415E+00;
    COFD[34] = -7.11992510E-02;
    COFD[35] = 3.14105973E-03;
    COFD[36] = -1.09469245E+01;
    COFD[37] = 2.30836460E+00;
    COFD[38] = -8.76339315E-02;
    COFD[39] = 3.90878445E-03;
    COFD[40] = -1.25098960E+01;
    COFD[41] = 2.77873601E+00;
    COFD[42] = -1.50637360E-01;
    COFD[43] = 6.72684281E-03;
    COFD[44] = -1.25098960E+01;
    COFD[45] = 2.77873601E+00;
    COFD[46] = -1.50637360E-01;
    COFD[47] = 6.72684281E-03;
    COFD[48] = -1.25141260E+01;
    COFD[49] = 2.77873601E+00;
    COFD[50] = -1.50637360E-01;
    COFD[51] = 6.72684281E-03;
    COFD[52] = -1.24693568E+01;
    COFD[53] = 2.76686648E+00;
    COFD[54] = -1.49120141E-01;
    COFD[55] = 6.66220432E-03;
    COFD[56] = -1.17159737E+01;
    COFD[57] = 2.48123210E+00;
    COFD[58] = -1.11322604E-01;
    COFD[59] = 4.99282389E-03;
    COFD[60] = -1.37794315E+01;
    COFD[61] = 3.23973858E+00;
    COFD[62] = -2.09989036E-01;
    COFD[63] = 9.27667906E-03;
    COFD[64] = -1.60517370E+01;
    COFD[65] = 4.11188603E+00;
    COFD[66] = -3.21540884E-01;
    COFD[67] = 1.40482564E-02;
    COFD[68] = -1.60528285E+01;
    COFD[69] = 4.11188603E+00;
    COFD[70] = -3.21540884E-01;
    COFD[71] = 1.40482564E-02;
    COFD[72] = -1.58456300E+01;
    COFD[73] = 4.02074783E+00;
    COFD[74] = -3.10018522E-01;
    COFD[75] = 1.35599552E-02;
    COFD[76] = -1.58456300E+01;
    COFD[77] = 4.02074783E+00;
    COFD[78] = -3.10018522E-01;
    COFD[79] = 1.35599552E-02;
    COFD[80] = -1.59537247E+01;
    COFD[81] = 4.07051484E+00;
    COFD[82] = -3.16303109E-01;
    COFD[83] = 1.38259377E-02;
    COFD[84] = -1.34695359E+01;
    COFD[85] = 3.09379603E+00;
    COFD[86] = -1.91268635E-01;
    COFD[87] = 8.47480224E-03;
    COFD[88] = -1.34709807E+01;
    COFD[89] = 3.09379603E+00;
    COFD[90] = -1.91268635E-01;
    COFD[91] = 8.47480224E-03;
    COFD[92] = -1.34723215E+01;
    COFD[93] = 3.09379603E+00;
    COFD[94] = -1.91268635E-01;
    COFD[95] = 8.47480224E-03;
    COFD[96] = -1.42229194E+01;
    COFD[97] = 3.38669384E+00;
    COFD[98] = -2.28784122E-01;
    COFD[99] = 1.00790953E-02;
    COFD[100] = -1.39913897E+01;
    COFD[101] = 3.26384506E+00;
    COFD[102] = -2.12947087E-01;
    COFD[103] = 9.39743888E-03;
    COFD[104] = -1.39924781E+01;
    COFD[105] = 3.26384506E+00;
    COFD[106] = -2.12947087E-01;
    COFD[107] = 9.39743888E-03;
    COFD[108] = -1.22004324E+01;
    COFD[109] = 2.80725489E+00;
    COFD[110] = -1.54291406E-01;
    COFD[111] = 6.88290911E-03;
    COFD[112] = -1.57034851E+01;
    COFD[113] = 3.93614244E+00;
    COFD[114] = -2.99111497E-01;
    COFD[115] = 1.30888229E-02;
    COFD[116] = -1.57034851E+01;
    COFD[117] = 3.93614244E+00;
    COFD[118] = -2.99111497E-01;
    COFD[119] = 1.30888229E-02;
    COFD[120] = -1.16906297E+01;
    COFD[121] = 2.47469981E+00;
    COFD[122] = -1.10436257E-01;
    COFD[123] = 4.95273813E-03;
    COFD[124] = -1.23130152E+01;
    COFD[125] = 2.74418790E+00;
    COFD[126] = -1.46230156E-01;
    COFD[127] = 6.53948886E-03;
    COFD[128] = -1.14366381E+01;
    COFD[129] = 2.78323501E+00;
    COFD[130] = -1.51214064E-01;
    COFD[131] = 6.75150012E-03;
    COFD[132] = -1.47968712E+01;
    COFD[133] = 4.23027636E+00;
    COFD[134] = -3.36139991E-01;
    COFD[135] = 1.46507621E-02;
    COFD[136] = -1.34230272E+01;
    COFD[137] = 3.48624238E+00;
    COFD[138] = -2.41554467E-01;
    COFD[139] = 1.06263545E-02;
    COFD[140] = -1.46550083E+01;
    COFD[141] = 3.83606243E+00;
    COFD[142] = -2.86076532E-01;
    COFD[143] = 1.25205829E-02;
    COFD[144] = -1.34247866E+01;
    COFD[145] = 3.48624238E+00;
    COFD[146] = -2.41554467E-01;
    COFD[147] = 1.06263545E-02;
    COFD[148] = -1.95739570E+01;
    COFD[149] = 5.61113230E+00;
    COFD[150] = -4.90190187E-01;
    COFD[151] = 2.03260675E-02;
    COFD[152] = -1.46554748E+01;
    COFD[153] = 3.83606243E+00;
    COFD[154] = -2.86076532E-01;
    COFD[155] = 1.25205829E-02;
    COFD[156] = -1.46559141E+01;
    COFD[157] = 3.83606243E+00;
    COFD[158] = -2.86076532E-01;
    COFD[159] = 1.25205829E-02;
    COFD[160] = -1.32467319E+01;
    COFD[161] = 3.34156587E+00;
    COFD[162] = -2.22853306E-01;
    COFD[163] = 9.81883417E-03;
    COFD[164] = -1.34162893E+01;
    COFD[165] = 3.48624238E+00;
    COFD[166] = -2.41554467E-01;
    COFD[167] = 1.06263545E-02;
    COFD[168] = -1.57972369E+01;
    COFD[169] = 4.22225052E+00;
    COFD[170] = -3.35156428E-01;
    COFD[171] = 1.46104855E-02;
    COFD[172] = -1.57972369E+01;
    COFD[173] = 4.22225052E+00;
    COFD[174] = -3.35156428E-01;
    COFD[175] = 1.46104855E-02;
    COFD[176] = -1.57994893E+01;
    COFD[177] = 4.22225052E+00;
    COFD[178] = -3.35156428E-01;
    COFD[179] = 1.46104855E-02;
    COFD[180] = -1.57199037E+01;
    COFD[181] = 4.19936335E+00;
    COFD[182] = -3.32311009E-01;
    COFD[183] = 1.44921003E-02;
    COFD[184] = -1.43151174E+01;
    COFD[185] = 3.68038508E+00;
    COFD[186] = -2.65779346E-01;
    COFD[187] = 1.16360771E-02;
    COFD[188] = -1.76147026E+01;
    COFD[189] = 4.86049500E+00;
    COFD[190] = -4.12200578E-01;
    COFD[191] = 1.77160971E-02;
    COFD[192] = -1.97544450E+01;
    COFD[193] = 5.56931926E+00;
    COFD[194] = -4.89105511E-01;
    COFD[195] = 2.04493129E-02;
    COFD[196] = -1.97550088E+01;
    COFD[197] = 5.56931926E+00;
    COFD[198] = -4.89105511E-01;
    COFD[199] = 2.04493129E-02;
    COFD[200] = -1.92718582E+01;
    COFD[201] = 5.41172124E+00;
    COFD[202] = -4.73213887E-01;
    COFD[203] = 1.99405473E-02;
    COFD[204] = -1.92718582E+01;
    COFD[205] = 5.41172124E+00;
    COFD[206] = -4.73213887E-01;
    COFD[207] = 1.99405473E-02;
    COFD[208] = -1.96866103E+01;
    COFD[209] = 5.54637286E+00;
    COFD[210] = -4.87070324E-01;
    COFD[211] = 2.03983467E-02;
    COFD[212] = -1.72224724E+01;
    COFD[213] = 4.69060745E+00;
    COFD[214] = -3.92369888E-01;
    COFD[215] = 1.69459661E-02;
    COFD[216] = -1.72232223E+01;
    COFD[217] = 4.69060745E+00;
    COFD[218] = -3.92369888E-01;
    COFD[219] = 1.69459661E-02;
    COFD[220] = -1.72239172E+01;
    COFD[221] = 4.69060745E+00;
    COFD[222] = -3.92369888E-01;
    COFD[223] = 1.69459661E-02;
    COFD[224] = -1.82251914E+01;
    COFD[225] = 5.05237312E+00;
    COFD[226] = -4.35182396E-01;
    COFD[227] = 1.86363074E-02;
    COFD[228] = -1.79339327E+01;
    COFD[229] = 4.91373893E+00;
    COFD[230] = -4.18747629E-01;
    COFD[231] = 1.79856610E-02;
    COFD[232] = -1.79344949E+01;
    COFD[233] = 4.91373893E+00;
    COFD[234] = -4.18747629E-01;
    COFD[235] = 1.79856610E-02;
    COFD[236] = -1.54460820E+01;
    COFD[237] = 4.26819983E+00;
    COFD[238] = -3.40766379E-01;
    COFD[239] = 1.48393361E-02;
    COFD[240] = -1.94688688E+01;
    COFD[241] = 5.43830787E+00;
    COFD[242] = -4.75472880E-01;
    COFD[243] = 1.99909996E-02;
    COFD[244] = -1.94688688E+01;
    COFD[245] = 5.43830787E+00;
    COFD[246] = -4.75472880E-01;
    COFD[247] = 1.99909996E-02;
    COFD[248] = -1.42894441E+01;
    COFD[249] = 3.67490723E+00;
    COFD[250] = -2.65114792E-01;
    COFD[251] = 1.16092671E-02;
    COFD[252] = -1.54738604E+01;
    COFD[253] = 4.15765300E+00;
    COFD[254] = -3.27126237E-01;
    COFD[255] = 1.42762611E-02;
    COFD[256] = -1.09595712E+01;
    COFD[257] = 2.30836460E+00;
    COFD[258] = -8.76339315E-02;
    COFD[259] = 3.90878445E-03;
    COFD[260] = -1.34230272E+01;
    COFD[261] = 3.48624238E+00;
    COFD[262] = -2.41554467E-01;
    COFD[263] = 1.06263545E-02;
    COFD[264] = -1.32093628E+01;
    COFD[265] = 2.90778936E+00;
    COFD[266] = -1.67388544E-01;
    COFD[267] = 7.45220609E-03;
    COFD[268] = -1.43139231E+01;
    COFD[269] = 3.17651319E+00;
    COFD[270] = -2.02028974E-01;
    COFD[271] = 8.94232502E-03;
    COFD[272] = -1.32244035E+01;
    COFD[273] = 2.90778936E+00;
    COFD[274] = -1.67388544E-01;
    COFD[275] = 7.45220609E-03;
    COFD[276] = -1.94093572E+01;
    COFD[277] = 5.16013126E+00;
    COFD[278] = -4.46824543E-01;
    COFD[279] = 1.90464887E-02;
    COFD[280] = -1.43190389E+01;
    COFD[281] = 3.17651319E+00;
    COFD[282] = -2.02028974E-01;
    COFD[283] = 8.94232502E-03;
    COFD[284] = -1.43238998E+01;
    COFD[285] = 3.17651319E+00;
    COFD[286] = -2.02028974E-01;
    COFD[287] = 8.94232502E-03;
    COFD[288] = -1.30610083E+01;
    COFD[289] = 2.80913567E+00;
    COFD[290] = -1.54536855E-01;
    COFD[291] = 6.89359313E-03;
    COFD[292] = -1.31551788E+01;
    COFD[293] = 2.90778936E+00;
    COFD[294] = -1.67388544E-01;
    COFD[295] = 7.45220609E-03;
    COFD[296] = -1.50584249E+01;
    COFD[297] = 3.47945612E+00;
    COFD[298] = -2.40703722E-01;
    COFD[299] = 1.05907441E-02;
    COFD[300] = -1.50584249E+01;
    COFD[301] = 3.47945612E+00;
    COFD[302] = -2.40703722E-01;
    COFD[303] = 1.05907441E-02;
    COFD[304] = -1.50766130E+01;
    COFD[305] = 3.47945612E+00;
    COFD[306] = -2.40703722E-01;
    COFD[307] = 1.05907441E-02;
    COFD[308] = -1.50270339E+01;
    COFD[309] = 3.46140064E+00;
    COFD[310] = -2.38440092E-01;
    COFD[311] = 1.04960087E-02;
    COFD[312] = -1.40999008E+01;
    COFD[313] = 3.08120012E+00;
    COFD[314] = -1.89629903E-01;
    COFD[315] = 8.40361952E-03;
    COFD[316] = -1.70534856E+01;
    COFD[317] = 4.14240922E+00;
    COFD[318] = -3.25239774E-01;
    COFD[319] = 1.41980687E-02;
    COFD[320] = -1.94313116E+01;
    COFD[321] = 5.02567894E+00;
    COFD[322] = -4.32045169E-01;
    COFD[323] = 1.85132214E-02;
    COFD[324] = -1.94373127E+01;
    COFD[325] = 5.02567894E+00;
    COFD[326] = -4.32045169E-01;
    COFD[327] = 1.85132214E-02;
    COFD[328] = -1.88179418E+01;
    COFD[329] = 4.79683898E+00;
    COFD[330] = -4.04829719E-01;
    COFD[331] = 1.74325475E-02;
    COFD[332] = -1.88179418E+01;
    COFD[333] = 4.79683898E+00;
    COFD[334] = -4.04829719E-01;
    COFD[335] = 1.74325475E-02;
    COFD[336] = -1.93364585E+01;
    COFD[337] = 4.98286777E+00;
    COFD[338] = -4.26970814E-01;
    COFD[339] = 1.83122917E-02;
    COFD[340] = -1.65412306E+01;
    COFD[341] = 3.95035840E+00;
    COFD[342] = -3.00959418E-01;
    COFD[343] = 1.31692593E-02;
    COFD[344] = -1.65488358E+01;
    COFD[345] = 3.95035840E+00;
    COFD[346] = -3.00959418E-01;
    COFD[347] = 1.31692593E-02;
    COFD[348] = -1.65559787E+01;
    COFD[349] = 3.95035840E+00;
    COFD[350] = -3.00959418E-01;
    COFD[351] = 1.31692593E-02;
    COFD[352] = -1.74792112E+01;
    COFD[353] = 4.29676909E+00;
    COFD[354] = -3.44085306E-01;
    COFD[355] = 1.49671135E-02;
    COFD[356] = -1.72496634E+01;
    COFD[357] = 4.17889917E+00;
    COFD[358] = -3.29752510E-01;
    COFD[359] = 1.43850275E-02;
    COFD[360] = -1.72556499E+01;
    COFD[361] = 4.17889917E+00;
    COFD[362] = -3.29752510E-01;
    COFD[363] = 1.43850275E-02;
    COFD[364] = -1.49500357E+01;
    COFD[365] = 3.52327209E+00;
    COFD[366] = -2.46286208E-01;
    COFD[367] = 1.08285963E-02;
    COFD[368] = -1.90883268E+01;
    COFD[369] = 4.84384483E+00;
    COFD[370] = -4.10265575E-01;
    COFD[371] = 1.76414287E-02;
    COFD[372] = -1.90883268E+01;
    COFD[373] = 4.84384483E+00;
    COFD[374] = -4.10265575E-01;
    COFD[375] = 1.76414287E-02;
    COFD[376] = -1.40756935E+01;
    COFD[377] = 3.07549274E+00;
    COFD[378] = -1.88889344E-01;
    COFD[379] = 8.37152866E-03;
    COFD[380] = -1.49610527E+01;
    COFD[381] = 3.41988961E+00;
    COFD[382] = -2.33128386E-01;
    COFD[383] = 1.02689994E-02;
    COFD[384] = -1.18988955E+01;
    COFD[385] = 2.57507000E+00;
    COFD[386] = -1.24033737E-01;
    COFD[387] = 5.56694959E-03;
    COFD[388] = -1.46550083E+01;
    COFD[389] = 3.83606243E+00;
    COFD[390] = -2.86076532E-01;
    COFD[391] = 1.25205829E-02;
    COFD[392] = -1.43139231E+01;
    COFD[393] = 3.17651319E+00;
    COFD[394] = -2.02028974E-01;
    COFD[395] = 8.94232502E-03;
    COFD[396] = -1.55511344E+01;
    COFD[397] = 3.48070094E+00;
    COFD[398] = -2.40859499E-01;
    COFD[399] = 1.05972514E-02;
    COFD[400] = -1.43340796E+01;
    COFD[401] = 3.17651319E+00;
    COFD[402] = -2.02028974E-01;
    COFD[403] = 8.94232502E-03;
    COFD[404] = -2.12652533E+01;
    COFD[405] = 5.59961818E+00;
    COFD[406] = -4.91624858E-01;
    COFD[407] = 2.05035550E-02;
    COFD[408] = -1.55588279E+01;
    COFD[409] = 3.48070094E+00;
    COFD[410] = -2.40859499E-01;
    COFD[411] = 1.05972514E-02;
    COFD[412] = -1.55661750E+01;
    COFD[413] = 3.48070094E+00;
    COFD[414] = -2.40859499E-01;
    COFD[415] = 1.05972514E-02;
    COFD[416] = -1.40707473E+01;
    COFD[417] = 3.05837263E+00;
    COFD[418] = -1.86672802E-01;
    COFD[419] = 8.27575734E-03;
    COFD[420] = -1.42429085E+01;
    COFD[421] = 3.17651319E+00;
    COFD[422] = -2.02028974E-01;
    COFD[423] = 8.94232502E-03;
    COFD[424] = -1.63254691E+01;
    COFD[425] = 3.82388595E+00;
    COFD[426] = -2.84480724E-01;
    COFD[427] = 1.24506311E-02;
    COFD[428] = -1.63254691E+01;
    COFD[429] = 3.82388595E+00;
    COFD[430] = -2.84480724E-01;
    COFD[431] = 1.24506311E-02;
    COFD[432] = -1.63493345E+01;
    COFD[433] = 3.82388595E+00;
    COFD[434] = -2.84480724E-01;
    COFD[435] = 1.24506311E-02;
    COFD[436] = -1.62724462E+01;
    COFD[437] = 3.79163564E+00;
    COFD[438] = -2.80257365E-01;
    COFD[439] = 1.22656902E-02;
    COFD[440] = -1.52721107E+01;
    COFD[441] = 3.36790500E+00;
    COFD[442] = -2.26321740E-01;
    COFD[443] = 9.97135055E-03;
    COFD[444] = -1.84688406E+01;
    COFD[445] = 4.49330851E+00;
    COFD[446] = -3.68208715E-01;
    COFD[447] = 1.59565402E-02;
    COFD[448] = -2.08204449E+01;
    COFD[449] = 5.35267674E+00;
    COFD[450] = -4.69010505E-01;
    COFD[451] = 1.98979152E-02;
    COFD[452] = -2.08293255E+01;
    COFD[453] = 5.35267674E+00;
    COFD[454] = -4.69010505E-01;
    COFD[455] = 1.98979152E-02;
    COFD[456] = -2.04928958E+01;
    COFD[457] = 5.22397933E+00;
    COFD[458] = -4.54138171E-01;
    COFD[459] = 1.93249285E-02;
    COFD[460] = -2.04928958E+01;
    COFD[461] = 5.22397933E+00;
    COFD[462] = -4.54138171E-01;
    COFD[463] = 1.93249285E-02;
    COFD[464] = -2.07595845E+01;
    COFD[465] = 5.32244593E+00;
    COFD[466] = -4.65829403E-01;
    COFD[467] = 1.97895274E-02;
    COFD[468] = -1.78725135E+01;
    COFD[469] = 4.29613154E+00;
    COFD[470] = -3.44012526E-01;
    COFD[471] = 1.49643715E-02;
    COFD[472] = -1.78834935E+01;
    COFD[473] = 4.29613154E+00;
    COFD[474] = -3.44012526E-01;
    COFD[475] = 1.49643715E-02;
    COFD[476] = -1.78938745E+01;
    COFD[477] = 4.29613154E+00;
    COFD[478] = -3.44012526E-01;
    COFD[479] = 1.49643715E-02;
    COFD[480] = -1.89544778E+01;
    COFD[481] = 4.68595732E+00;
    COFD[482] = -3.91842840E-01;
    COFD[483] = 1.69262542E-02;
    COFD[484] = -1.86335932E+01;
    COFD[485] = 4.53572533E+00;
    COFD[486] = -3.73386925E-01;
    COFD[487] = 1.61678881E-02;
    COFD[488] = -1.86424545E+01;
    COFD[489] = 4.53572533E+00;
    COFD[490] = -3.73386925E-01;
    COFD[491] = 1.61678881E-02;
    COFD[492] = -1.64169433E+01;
    COFD[493] = 3.89309916E+00;
    COFD[494] = -2.93528188E-01;
    COFD[495] = 1.28463177E-02;
    COFD[496] = -2.05184870E+01;
    COFD[497] = 5.18417470E+00;
    COFD[498] = -4.49491573E-01;
    COFD[499] = 1.91438508E-02;
    COFD[500] = -2.05184870E+01;
    COFD[501] = 5.18417470E+00;
    COFD[502] = -4.49491573E-01;
    COFD[503] = 1.91438508E-02;
    COFD[504] = -1.52414485E+01;
    COFD[505] = 3.35922578E+00;
    COFD[506] = -2.25181399E-01;
    COFD[507] = 9.92132878E-03;
    COFD[508] = -1.62380075E+01;
    COFD[509] = 3.72612300E+00;
    COFD[510] = -2.71663673E-01;
    COFD[511] = 1.18889643E-02;
    COFD[512] = -1.09628982E+01;
    COFD[513] = 2.30836460E+00;
    COFD[514] = -8.76339315E-02;
    COFD[515] = 3.90878445E-03;
    COFD[516] = -1.34247866E+01;
    COFD[517] = 3.48624238E+00;
    COFD[518] = -2.41554467E-01;
    COFD[519] = 1.06263545E-02;
    COFD[520] = -1.32244035E+01;
    COFD[521] = 2.90778936E+00;
    COFD[522] = -1.67388544E-01;
    COFD[523] = 7.45220609E-03;
    COFD[524] = -1.43340796E+01;
    COFD[525] = 3.17651319E+00;
    COFD[526] = -2.02028974E-01;
    COFD[527] = 8.94232502E-03;
    COFD[528] = -1.32399106E+01;
    COFD[529] = 2.90778936E+00;
    COFD[530] = -1.67388544E-01;
    COFD[531] = 7.45220609E-03;
    COFD[532] = -1.94253036E+01;
    COFD[533] = 5.16013126E+00;
    COFD[534] = -4.46824543E-01;
    COFD[535] = 1.90464887E-02;
    COFD[536] = -1.43394069E+01;
    COFD[537] = 3.17651319E+00;
    COFD[538] = -2.02028974E-01;
    COFD[539] = 8.94232502E-03;
    COFD[540] = -1.43444709E+01;
    COFD[541] = 3.17651319E+00;
    COFD[542] = -2.02028974E-01;
    COFD[543] = 8.94232502E-03;
    COFD[544] = -1.30738796E+01;
    COFD[545] = 2.80913567E+00;
    COFD[546] = -1.54536855E-01;
    COFD[547] = 6.89359313E-03;
    COFD[548] = -1.31686537E+01;
    COFD[549] = 2.90778936E+00;
    COFD[550] = -1.67388544E-01;
    COFD[551] = 7.45220609E-03;
    COFD[552] = -1.50724636E+01;
    COFD[553] = 3.47945612E+00;
    COFD[554] = -2.40703722E-01;
    COFD[555] = 1.05907441E-02;
    COFD[556] = -1.50724636E+01;
    COFD[557] = 3.47945612E+00;
    COFD[558] = -2.40703722E-01;
    COFD[559] = 1.05907441E-02;
    COFD[560] = -1.50911794E+01;
    COFD[561] = 3.47945612E+00;
    COFD[562] = -2.40703722E-01;
    COFD[563] = 1.05907441E-02;
    COFD[564] = -1.50420953E+01;
    COFD[565] = 3.46140064E+00;
    COFD[566] = -2.38440092E-01;
    COFD[567] = 1.04960087E-02;
    COFD[568] = -1.41191261E+01;
    COFD[569] = 3.08120012E+00;
    COFD[570] = -1.89629903E-01;
    COFD[571] = 8.40361952E-03;
    COFD[572] = -1.70757047E+01;
    COFD[573] = 4.14240922E+00;
    COFD[574] = -3.25239774E-01;
    COFD[575] = 1.41980687E-02;
    COFD[576] = -1.94507876E+01;
    COFD[577] = 5.02567894E+00;
    COFD[578] = -4.32045169E-01;
    COFD[579] = 1.85132214E-02;
    COFD[580] = -1.94570287E+01;
    COFD[581] = 5.02567894E+00;
    COFD[582] = -4.32045169E-01;
    COFD[583] = 1.85132214E-02;
    COFD[584] = -1.88378874E+01;
    COFD[585] = 4.79683898E+00;
    COFD[586] = -4.04829719E-01;
    COFD[587] = 1.74325475E-02;
    COFD[588] = -1.88378874E+01;
    COFD[589] = 4.79683898E+00;
    COFD[590] = -4.04829719E-01;
    COFD[591] = 1.74325475E-02;
    COFD[592] = -1.93566243E+01;
    COFD[593] = 4.98286777E+00;
    COFD[594] = -4.26970814E-01;
    COFD[595] = 1.83122917E-02;
    COFD[596] = -1.65596434E+01;
    COFD[597] = 3.95035840E+00;
    COFD[598] = -3.00959418E-01;
    COFD[599] = 1.31692593E-02;
    COFD[600] = -1.65675362E+01;
    COFD[601] = 3.95035840E+00;
    COFD[602] = -3.00959418E-01;
    COFD[603] = 1.31692593E-02;
    COFD[604] = -1.65749533E+01;
    COFD[605] = 3.95035840E+00;
    COFD[606] = -3.00959418E-01;
    COFD[607] = 1.31692593E-02;
    COFD[608] = -1.74984476E+01;
    COFD[609] = 4.29676909E+00;
    COFD[610] = -3.44085306E-01;
    COFD[611] = 1.49671135E-02;
    COFD[612] = -1.72691500E+01;
    COFD[613] = 4.17889917E+00;
    COFD[614] = -3.29752510E-01;
    COFD[615] = 1.43850275E-02;
    COFD[616] = -1.72753760E+01;
    COFD[617] = 4.17889917E+00;
    COFD[618] = -3.29752510E-01;
    COFD[619] = 1.43850275E-02;
    COFD[620] = -1.49718233E+01;
    COFD[621] = 3.52327209E+00;
    COFD[622] = -2.46286208E-01;
    COFD[623] = 1.08285963E-02;
    COFD[624] = -1.91102652E+01;
    COFD[625] = 4.84384483E+00;
    COFD[626] = -4.10265575E-01;
    COFD[627] = 1.76414287E-02;
    COFD[628] = -1.91102652E+01;
    COFD[629] = 4.84384483E+00;
    COFD[630] = -4.10265575E-01;
    COFD[631] = 1.76414287E-02;
    COFD[632] = -1.40949196E+01;
    COFD[633] = 3.07549274E+00;
    COFD[634] = -1.88889344E-01;
    COFD[635] = 8.37152866E-03;
    COFD[636] = -1.49826725E+01;
    COFD[637] = 3.41988961E+00;
    COFD[638] = -2.33128386E-01;
    COFD[639] = 1.02689994E-02;
    COFD[640] = -1.71982995E+01;
    COFD[641] = 4.63881404E+00;
    COFD[642] = -3.86139633E-01;
    COFD[643] = 1.66955081E-02;
    COFD[644] = -1.95739570E+01;
    COFD[645] = 5.61113230E+00;
    COFD[646] = -4.90190187E-01;
    COFD[647] = 2.03260675E-02;
    COFD[648] = -1.94093572E+01;
    COFD[649] = 5.16013126E+00;
    COFD[650] = -4.46824543E-01;
    COFD[651] = 1.90464887E-02;
    COFD[652] = -2.12652533E+01;
    COFD[653] = 5.59961818E+00;
    COFD[654] = -4.91624858E-01;
    COFD[655] = 2.05035550E-02;
    COFD[656] = -1.94253036E+01;
    COFD[657] = 5.16013126E+00;
    COFD[658] = -4.46824543E-01;
    COFD[659] = 1.90464887E-02;
    COFD[660] = -1.19157919E+01;
    COFD[661] = 9.28955130E-01;
    COFD[662] = 2.42107090E-01;
    COFD[663] = -1.59823963E-02;
    COFD[664] = -2.06463744E+01;
    COFD[665] = 5.41688482E+00;
    COFD[666] = -4.73387188E-01;
    COFD[667] = 1.99280175E-02;
    COFD[668] = -2.06516336E+01;
    COFD[669] = 5.41688482E+00;
    COFD[670] = -4.73387188E-01;
    COFD[671] = 1.99280175E-02;
    COFD[672] = -1.92003817E+01;
    COFD[673] = 5.05708283E+00;
    COFD[674] = -4.35739290E-01;
    COFD[675] = 1.86583205E-02;
    COFD[676] = -1.93521390E+01;
    COFD[677] = 5.16013126E+00;
    COFD[678] = -4.46824543E-01;
    COFD[679] = 1.90464887E-02;
    COFD[680] = -2.12639214E+01;
    COFD[681] = 5.61184117E+00;
    COFD[682] = -4.90532156E-01;
    COFD[683] = 2.03507922E-02;
    COFD[684] = -2.12639214E+01;
    COFD[685] = 5.61184117E+00;
    COFD[686] = -4.90532156E-01;
    COFD[687] = 2.03507922E-02;
    COFD[688] = -2.12831323E+01;
    COFD[689] = 5.61184117E+00;
    COFD[690] = -4.90532156E-01;
    COFD[691] = 2.03507922E-02;
    COFD[692] = -2.14087397E+01;
    COFD[693] = 5.57282008E+00;
    COFD[694] = -4.76690890E-01;
    COFD[695] = 1.94000719E-02;
    COFD[696] = -2.11388331E+01;
    COFD[697] = 5.55529675E+00;
    COFD[698] = -4.87942518E-01;
    COFD[699] = 2.04249054E-02;
    COFD[700] = -2.07653719E+01;
    COFD[701] = 5.01092022E+00;
    COFD[702] = -3.77985635E-01;
    COFD[703] = 1.40968645E-02;
    COFD[704] = -1.77498543E+01;
    COFD[705] = 3.57475686E+00;
    COFD[706] = -1.56396297E-01;
    COFD[707] = 3.12157721E-03;
    COFD[708] = -1.77563250E+01;
    COFD[709] = 3.57475686E+00;
    COFD[710] = -1.56396297E-01;
    COFD[711] = 3.12157721E-03;
    COFD[712] = -1.65295288E+01;
    COFD[713] = 2.97569206E+00;
    COFD[714] = -6.75652842E-02;
    COFD[715] = -1.08648422E-03;
    COFD[716] = -1.65295288E+01;
    COFD[717] = 2.97569206E+00;
    COFD[718] = -6.75652842E-02;
    COFD[719] = -1.08648422E-03;
    COFD[720] = -1.80253664E+01;
    COFD[721] = 3.69199168E+00;
    COFD[722] = -1.74005516E-01;
    COFD[723] = 3.97694372E-03;
    COFD[724] = -2.15014310E+01;
    COFD[725] = 5.46737673E+00;
    COFD[726] = -4.55696085E-01;
    COFD[727] = 1.81982625E-02;
    COFD[728] = -2.15095980E+01;
    COFD[729] = 5.46737673E+00;
    COFD[730] = -4.55696085E-01;
    COFD[731] = 1.81982625E-02;
    COFD[732] = -2.15172770E+01;
    COFD[733] = 5.46737673E+00;
    COFD[734] = -4.55696085E-01;
    COFD[735] = 1.81982625E-02;
    COFD[736] = -2.08812333E+01;
    COFD[737] = 5.08859217E+00;
    COFD[738] = -3.90525428E-01;
    COFD[739] = 1.47376395E-02;
    COFD[740] = -2.12597312E+01;
    COFD[741] = 5.24930667E+00;
    COFD[742] = -4.17435088E-01;
    COFD[743] = 1.61434424E-02;
    COFD[744] = -2.12661865E+01;
    COFD[745] = 5.24930667E+00;
    COFD[746] = -4.17435088E-01;
    COFD[747] = 1.61434424E-02;
    COFD[748] = -2.10440675E+01;
    COFD[749] = 5.59806282E+00;
    COFD[750] = -4.87109535E-01;
    COFD[751] = 2.01370226E-02;
    COFD[752] = -1.87383952E+01;
    COFD[753] = 3.96926341E+00;
    COFD[754] = -2.16412264E-01;
    COFD[755] = 6.06012078E-03;
    COFD[756] = -1.87383952E+01;
    COFD[757] = 3.96926341E+00;
    COFD[758] = -2.16412264E-01;
    COFD[759] = 6.06012078E-03;
    COFD[760] = -2.10643259E+01;
    COFD[761] = 5.53614847E+00;
    COFD[762] = -4.86046736E-01;
    COFD[763] = 2.03659188E-02;
    COFD[764] = -2.12755888E+01;
    COFD[765] = 5.60381989E+00;
    COFD[766] = -4.91225459E-01;
    COFD[767] = 2.04487844E-02;
    COFD[768] = -1.18998012E+01;
    COFD[769] = 2.57507000E+00;
    COFD[770] = -1.24033737E-01;
    COFD[771] = 5.56694959E-03;
    COFD[772] = -1.46554748E+01;
    COFD[773] = 3.83606243E+00;
    COFD[774] = -2.86076532E-01;
    COFD[775] = 1.25205829E-02;
    COFD[776] = -1.43190389E+01;
    COFD[777] = 3.17651319E+00;
    COFD[778] = -2.02028974E-01;
    COFD[779] = 8.94232502E-03;
    COFD[780] = -1.55588279E+01;
    COFD[781] = 3.48070094E+00;
    COFD[782] = -2.40859499E-01;
    COFD[783] = 1.05972514E-02;
    COFD[784] = -1.43394069E+01;
    COFD[785] = 3.17651319E+00;
    COFD[786] = -2.02028974E-01;
    COFD[787] = 8.94232502E-03;
    COFD[788] = -2.06463744E+01;
    COFD[789] = 5.41688482E+00;
    COFD[790] = -4.73387188E-01;
    COFD[791] = 1.99280175E-02;
    COFD[792] = -1.55666415E+01;
    COFD[793] = 3.48070094E+00;
    COFD[794] = -2.40859499E-01;
    COFD[795] = 1.05972514E-02;
    COFD[796] = -1.55741053E+01;
    COFD[797] = 3.48070094E+00;
    COFD[798] = -2.40859499E-01;
    COFD[799] = 1.05972514E-02;
    COFD[800] = -1.40749320E+01;
    COFD[801] = 3.05837263E+00;
    COFD[802] = -1.86672802E-01;
    COFD[803] = 8.27575734E-03;
    COFD[804] = -1.42473439E+01;
    COFD[805] = 3.17651319E+00;
    COFD[806] = -2.02028974E-01;
    COFD[807] = 8.94232502E-03;
    COFD[808] = -1.63301444E+01;
    COFD[809] = 3.82388595E+00;
    COFD[810] = -2.84480724E-01;
    COFD[811] = 1.24506311E-02;
    COFD[812] = -1.63301444E+01;
    COFD[813] = 3.82388595E+00;
    COFD[814] = -2.84480724E-01;
    COFD[815] = 1.24506311E-02;
    COFD[816] = -1.63542394E+01;
    COFD[817] = 3.82388595E+00;
    COFD[818] = -2.84480724E-01;
    COFD[819] = 1.24506311E-02;
    COFD[820] = -1.62775714E+01;
    COFD[821] = 3.79163564E+00;
    COFD[822] = -2.80257365E-01;
    COFD[823] = 1.22656902E-02;
    COFD[824] = -1.52792891E+01;
    COFD[825] = 3.36790500E+00;
    COFD[826] = -2.26321740E-01;
    COFD[827] = 9.97135055E-03;
    COFD[828] = -1.84777607E+01;
    COFD[829] = 4.49330851E+00;
    COFD[830] = -3.68208715E-01;
    COFD[831] = 1.59565402E-02;
    COFD[832] = -2.08277598E+01;
    COFD[833] = 5.35267674E+00;
    COFD[834] = -4.69010505E-01;
    COFD[835] = 1.98979152E-02;
    COFD[836] = -2.08367725E+01;
    COFD[837] = 5.35267674E+00;
    COFD[838] = -4.69010505E-01;
    COFD[839] = 1.98979152E-02;
    COFD[840] = -2.02637994E+01;
    COFD[841] = 5.14984081E+00;
    COFD[842] = -4.46093018E-01;
    COFD[843] = 1.90396647E-02;
    COFD[844] = -2.02637994E+01;
    COFD[845] = 5.14984081E+00;
    COFD[846] = -4.46093018E-01;
    COFD[847] = 1.90396647E-02;
    COFD[848] = -2.07672833E+01;
    COFD[849] = 5.32244593E+00;
    COFD[850] = -4.65829403E-01;
    COFD[851] = 1.97895274E-02;
    COFD[852] = -1.78792605E+01;
    COFD[853] = 4.29613154E+00;
    COFD[854] = -3.44012526E-01;
    COFD[855] = 1.49643715E-02;
    COFD[856] = -1.78903913E+01;
    COFD[857] = 4.29613154E+00;
    COFD[858] = -3.44012526E-01;
    COFD[859] = 1.49643715E-02;
    COFD[860] = -1.79009181E+01;
    COFD[861] = 4.29613154E+00;
    COFD[862] = -3.44012526E-01;
    COFD[863] = 1.49643715E-02;
    COFD[864] = -1.89616623E+01;
    COFD[865] = 4.68595732E+00;
    COFD[866] = -3.91842840E-01;
    COFD[867] = 1.69262542E-02;
    COFD[868] = -1.86409139E+01;
    COFD[869] = 4.53572533E+00;
    COFD[870] = -3.73386925E-01;
    COFD[871] = 1.61678881E-02;
    COFD[872] = -1.86499071E+01;
    COFD[873] = 4.53572533E+00;
    COFD[874] = -3.73386925E-01;
    COFD[875] = 1.61678881E-02;
    COFD[876] = -1.64255964E+01;
    COFD[877] = 3.89309916E+00;
    COFD[878] = -2.93528188E-01;
    COFD[879] = 1.28463177E-02;
    COFD[880] = -2.05272328E+01;
    COFD[881] = 5.18417470E+00;
    COFD[882] = -4.49491573E-01;
    COFD[883] = 1.91438508E-02;
    COFD[884] = -2.05272328E+01;
    COFD[885] = 5.18417470E+00;
    COFD[886] = -4.49491573E-01;
    COFD[887] = 1.91438508E-02;
    COFD[888] = -1.52486273E+01;
    COFD[889] = 3.35922578E+00;
    COFD[890] = -2.25181399E-01;
    COFD[891] = 9.92132878E-03;
    COFD[892] = -1.62465583E+01;
    COFD[893] = 3.72612300E+00;
    COFD[894] = -2.71663673E-01;
    COFD[895] = 1.18889643E-02;
    COFD[896] = -1.19006548E+01;
    COFD[897] = 2.57507000E+00;
    COFD[898] = -1.24033737E-01;
    COFD[899] = 5.56694959E-03;
    COFD[900] = -1.46559141E+01;
    COFD[901] = 3.83606243E+00;
    COFD[902] = -2.86076532E-01;
    COFD[903] = 1.25205829E-02;
    COFD[904] = -1.43238998E+01;
    COFD[905] = 3.17651319E+00;
    COFD[906] = -2.02028974E-01;
    COFD[907] = 8.94232502E-03;
    COFD[908] = -1.55661750E+01;
    COFD[909] = 3.48070094E+00;
    COFD[910] = -2.40859499E-01;
    COFD[911] = 1.05972514E-02;
    COFD[912] = -1.43444709E+01;
    COFD[913] = 3.17651319E+00;
    COFD[914] = -2.02028974E-01;
    COFD[915] = 8.94232502E-03;
    COFD[916] = -2.06516336E+01;
    COFD[917] = 5.41688482E+00;
    COFD[918] = -4.73387188E-01;
    COFD[919] = 1.99280175E-02;
    COFD[920] = -1.55741053E+01;
    COFD[921] = 3.48070094E+00;
    COFD[922] = -2.40859499E-01;
    COFD[923] = 1.05972514E-02;
    COFD[924] = -1.55816822E+01;
    COFD[925] = 3.48070094E+00;
    COFD[926] = -2.40859499E-01;
    COFD[927] = 1.05972514E-02;
    COFD[928] = -1.40789009E+01;
    COFD[929] = 3.05837263E+00;
    COFD[930] = -1.86672802E-01;
    COFD[931] = 8.27575734E-03;
    COFD[932] = -1.42515527E+01;
    COFD[933] = 3.17651319E+00;
    COFD[934] = -2.02028974E-01;
    COFD[935] = 8.94232502E-03;
    COFD[936] = -1.63345829E+01;
    COFD[937] = 3.82388595E+00;
    COFD[938] = -2.84480724E-01;
    COFD[939] = 1.24506311E-02;
    COFD[940] = -1.63345829E+01;
    COFD[941] = 3.82388595E+00;
    COFD[942] = -2.84480724E-01;
    COFD[943] = 1.24506311E-02;
    COFD[944] = -1.63588981E+01;
    COFD[945] = 3.82388595E+00;
    COFD[946] = -2.84480724E-01;
    COFD[947] = 1.24506311E-02;
    COFD[948] = -1.62824412E+01;
    COFD[949] = 3.79163564E+00;
    COFD[950] = -2.80257365E-01;
    COFD[951] = 1.22656902E-02;
    COFD[952] = -1.52861376E+01;
    COFD[953] = 3.36790500E+00;
    COFD[954] = -2.26321740E-01;
    COFD[955] = 9.97135055E-03;
    COFD[956] = -1.84863000E+01;
    COFD[957] = 4.49330851E+00;
    COFD[958] = -3.68208715E-01;
    COFD[959] = 1.59565402E-02;
    COFD[960] = -2.08347403E+01;
    COFD[961] = 5.35267674E+00;
    COFD[962] = -4.69010505E-01;
    COFD[963] = 1.98979152E-02;
    COFD[964] = -2.08438809E+01;
    COFD[965] = 5.35267674E+00;
    COFD[966] = -4.69010505E-01;
    COFD[967] = 1.98979152E-02;
    COFD[968] = -2.02710316E+01;
    COFD[969] = 5.14984081E+00;
    COFD[970] = -4.46093018E-01;
    COFD[971] = 1.90396647E-02;
    COFD[972] = -2.02710316E+01;
    COFD[973] = 5.14984081E+00;
    COFD[974] = -4.46093018E-01;
    COFD[975] = 1.90396647E-02;
    COFD[976] = -2.07746356E+01;
    COFD[977] = 5.32244593E+00;
    COFD[978] = -4.65829403E-01;
    COFD[979] = 1.97895274E-02;
    COFD[980] = -1.78856918E+01;
    COFD[981] = 4.29613154E+00;
    COFD[982] = -3.44012526E-01;
    COFD[983] = 1.49643715E-02;
    COFD[984] = -1.78969684E+01;
    COFD[985] = 4.29613154E+00;
    COFD[986] = -3.44012526E-01;
    COFD[987] = 1.49643715E-02;
    COFD[988] = -1.79076361E+01;
    COFD[989] = 4.29613154E+00;
    COFD[990] = -3.44012526E-01;
    COFD[991] = 1.49643715E-02;
    COFD[992] = -1.89685165E+01;
    COFD[993] = 4.68595732E+00;
    COFD[994] = -3.91842840E-01;
    COFD[995] = 1.69262542E-02;
    COFD[996] = -1.86479000E+01;
    COFD[997] = 4.53572533E+00;
    COFD[998] = -3.73386925E-01;
    COFD[999] = 1.61678881E-02;
    COFD[1000] = -1.86570209E+01;
    COFD[1001] = 4.53572533E+00;
    COFD[1002] = -3.73386925E-01;
    COFD[1003] = 1.61678881E-02;
    COFD[1004] = -1.64338757E+01;
    COFD[1005] = 3.89309916E+00;
    COFD[1006] = -2.93528188E-01;
    COFD[1007] = 1.28463177E-02;
    COFD[1008] = -2.05356023E+01;
    COFD[1009] = 5.18417470E+00;
    COFD[1010] = -4.49491573E-01;
    COFD[1011] = 1.91438508E-02;
    COFD[1012] = -2.05356023E+01;
    COFD[1013] = 5.18417470E+00;
    COFD[1014] = -4.49491573E-01;
    COFD[1015] = 1.91438508E-02;
    COFD[1016] = -1.52554761E+01;
    COFD[1017] = 3.35922578E+00;
    COFD[1018] = -2.25181399E-01;
    COFD[1019] = 9.92132878E-03;
    COFD[1020] = -1.62547381E+01;
    COFD[1021] = 3.72612300E+00;
    COFD[1022] = -2.71663673E-01;
    COFD[1023] = 1.18889643E-02;
    COFD[1024] = -1.08369481E+01;
    COFD[1025] = 2.19094415E+00;
    COFD[1026] = -7.11992510E-02;
    COFD[1027] = 3.14105973E-03;
    COFD[1028] = -1.32467319E+01;
    COFD[1029] = 3.34156587E+00;
    COFD[1030] = -2.22853306E-01;
    COFD[1031] = 9.81883417E-03;
    COFD[1032] = -1.30610083E+01;
    COFD[1033] = 2.80913567E+00;
    COFD[1034] = -1.54536855E-01;
    COFD[1035] = 6.89359313E-03;
    COFD[1036] = -1.40707473E+01;
    COFD[1037] = 3.05837263E+00;
    COFD[1038] = -1.86672802E-01;
    COFD[1039] = 8.27575734E-03;
    COFD[1040] = -1.30738796E+01;
    COFD[1041] = 2.80913567E+00;
    COFD[1042] = -1.54536855E-01;
    COFD[1043] = 6.89359313E-03;
    COFD[1044] = -1.92003817E+01;
    COFD[1045] = 5.05708283E+00;
    COFD[1046] = -4.35739290E-01;
    COFD[1047] = 1.86583205E-02;
    COFD[1048] = -1.40749320E+01;
    COFD[1049] = 3.05837263E+00;
    COFD[1050] = -1.86672802E-01;
    COFD[1051] = 8.27575734E-03;
    COFD[1052] = -1.40789009E+01;
    COFD[1053] = 3.05837263E+00;
    COFD[1054] = -1.86672802E-01;
    COFD[1055] = 8.27575734E-03;
    COFD[1056] = -1.29573076E+01;
    COFD[1057] = 2.73155251E+00;
    COFD[1058] = -1.44594082E-01;
    COFD[1059] = 6.46883252E-03;
    COFD[1060] = -1.30141899E+01;
    COFD[1061] = 2.80913567E+00;
    COFD[1062] = -1.54536855E-01;
    COFD[1063] = 6.89359313E-03;
    COFD[1064] = -1.47558850E+01;
    COFD[1065] = 3.33113524E+00;
    COFD[1066] = -2.21479057E-01;
    COFD[1067] = 9.75837737E-03;
    COFD[1068] = -1.47558850E+01;
    COFD[1069] = 3.33113524E+00;
    COFD[1070] = -2.21479057E-01;
    COFD[1071] = 9.75837737E-03;
    COFD[1072] = -1.47715918E+01;
    COFD[1073] = 3.33113524E+00;
    COFD[1074] = -2.21479057E-01;
    COFD[1075] = 9.75837737E-03;
    COFD[1076] = -1.47048024E+01;
    COFD[1077] = 3.30594991E+00;
    COFD[1078] = -2.18182207E-01;
    COFD[1079] = 9.61429447E-03;
    COFD[1080] = -1.38862839E+01;
    COFD[1081] = 2.97564184E+00;
    COFD[1082] = -1.76025309E-01;
    COFD[1083] = 7.81869993E-03;
    COFD[1084] = -1.67390425E+01;
    COFD[1085] = 4.00828594E+00;
    COFD[1086] = -3.08414344E-01;
    COFD[1087] = 1.34907430E-02;
    COFD[1088] = -1.90526941E+01;
    COFD[1089] = 4.86821670E+00;
    COFD[1090] = -4.13144121E-01;
    COFD[1091] = 1.77546701E-02;
    COFD[1092] = -1.90576320E+01;
    COFD[1093] = 4.86821670E+00;
    COFD[1094] = -4.13144121E-01;
    COFD[1095] = 1.77546701E-02;
    COFD[1096] = -1.85157414E+01;
    COFD[1097] = 4.67076124E+00;
    COFD[1098] = -3.90022427E-01;
    COFD[1099] = 1.68533953E-02;
    COFD[1100] = -1.85157414E+01;
    COFD[1101] = 4.67076124E+00;
    COFD[1102] = -3.90022427E-01;
    COFD[1103] = 1.68533953E-02;
    COFD[1104] = -1.89671752E+01;
    COFD[1105] = 4.83076737E+00;
    COFD[1106] = -4.08802573E-01;
    COFD[1107] = 1.75875241E-02;
    COFD[1108] = -1.61038806E+01;
    COFD[1109] = 3.75910622E+00;
    COFD[1110] = -2.75986578E-01;
    COFD[1111] = 1.20782843E-02;
    COFD[1112] = -1.61101966E+01;
    COFD[1113] = 3.75910622E+00;
    COFD[1114] = -2.75986578E-01;
    COFD[1115] = 1.20782843E-02;
    COFD[1116] = -1.61161138E+01;
    COFD[1117] = 3.75910622E+00;
    COFD[1118] = -2.75986578E-01;
    COFD[1119] = 1.20782843E-02;
    COFD[1120] = -1.71884218E+01;
    COFD[1121] = 4.17190426E+00;
    COFD[1122] = -3.28894681E-01;
    COFD[1123] = 1.43498101E-02;
    COFD[1124] = -1.69491115E+01;
    COFD[1125] = 4.05099737E+00;
    COFD[1126] = -3.13841660E-01;
    COFD[1127] = 1.37218854E-02;
    COFD[1128] = -1.69540369E+01;
    COFD[1129] = 4.05099737E+00;
    COFD[1130] = -3.13841660E-01;
    COFD[1131] = 1.37218854E-02;
    COFD[1132] = -1.46907028E+01;
    COFD[1133] = 3.39229020E+00;
    COFD[1134] = -2.29520232E-01;
    COFD[1135] = 1.01114311E-02;
    COFD[1136] = -1.87688110E+01;
    COFD[1137] = 4.71729964E+00;
    COFD[1138] = -3.95432573E-01;
    COFD[1139] = 1.70623691E-02;
    COFD[1140] = -1.87688110E+01;
    COFD[1141] = 4.71729964E+00;
    COFD[1142] = -3.95432573E-01;
    COFD[1143] = 1.70623691E-02;
    COFD[1144] = -1.38661480E+01;
    COFD[1145] = 2.97137588E+00;
    COFD[1146] = -1.75491257E-01;
    COFD[1147] = 7.79646773E-03;
    COFD[1148] = -1.46461881E+01;
    COFD[1149] = 3.27505697E+00;
    COFD[1150] = -2.14306851E-01;
    COFD[1151] = 9.45219335E-03;
    COFD[1152] = -1.09469245E+01;
    COFD[1153] = 2.30836460E+00;
    COFD[1154] = -8.76339315E-02;
    COFD[1155] = 3.90878445E-03;
    COFD[1156] = -1.34162893E+01;
    COFD[1157] = 3.48624238E+00;
    COFD[1158] = -2.41554467E-01;
    COFD[1159] = 1.06263545E-02;
    COFD[1160] = -1.31551788E+01;
    COFD[1161] = 2.90778936E+00;
    COFD[1162] = -1.67388544E-01;
    COFD[1163] = 7.45220609E-03;
    COFD[1164] = -1.42429085E+01;
    COFD[1165] = 3.17651319E+00;
    COFD[1166] = -2.02028974E-01;
    COFD[1167] = 8.94232502E-03;
    COFD[1168] = -1.31686537E+01;
    COFD[1169] = 2.90778936E+00;
    COFD[1170] = -1.67388544E-01;
    COFD[1171] = 7.45220609E-03;
    COFD[1172] = -1.93521390E+01;
    COFD[1173] = 5.16013126E+00;
    COFD[1174] = -4.46824543E-01;
    COFD[1175] = 1.90464887E-02;
    COFD[1176] = -1.42473439E+01;
    COFD[1177] = 3.17651319E+00;
    COFD[1178] = -2.02028974E-01;
    COFD[1179] = 8.94232502E-03;
    COFD[1180] = -1.42515527E+01;
    COFD[1181] = 3.17651319E+00;
    COFD[1182] = -2.02028974E-01;
    COFD[1183] = 8.94232502E-03;
    COFD[1184] = -1.30141899E+01;
    COFD[1185] = 2.80913567E+00;
    COFD[1186] = -1.54536855E-01;
    COFD[1187] = 6.89359313E-03;
    COFD[1188] = -1.31062967E+01;
    COFD[1189] = 2.90778936E+00;
    COFD[1190] = -1.67388544E-01;
    COFD[1191] = 7.45220609E-03;
    COFD[1192] = -1.50076254E+01;
    COFD[1193] = 3.47945612E+00;
    COFD[1194] = -2.40703722E-01;
    COFD[1195] = 1.05907441E-02;
    COFD[1196] = -1.50076254E+01;
    COFD[1197] = 3.47945612E+00;
    COFD[1198] = -2.40703722E-01;
    COFD[1199] = 1.05907441E-02;
    COFD[1200] = -1.50240272E+01;
    COFD[1201] = 3.47945612E+00;
    COFD[1202] = -2.40703722E-01;
    COFD[1203] = 1.05907441E-02;
    COFD[1204] = -1.49727799E+01;
    COFD[1205] = 3.46140064E+00;
    COFD[1206] = -2.38440092E-01;
    COFD[1207] = 1.04960087E-02;
    COFD[1208] = -1.40318948E+01;
    COFD[1209] = 3.08120012E+00;
    COFD[1210] = -1.89629903E-01;
    COFD[1211] = 8.40361952E-03;
    COFD[1212] = -1.69758891E+01;
    COFD[1213] = 4.14240922E+00;
    COFD[1214] = -3.25239774E-01;
    COFD[1215] = 1.41980687E-02;
    COFD[1216] = -1.93624931E+01;
    COFD[1217] = 5.02567894E+00;
    COFD[1218] = -4.32045169E-01;
    COFD[1219] = 1.85132214E-02;
    COFD[1220] = -1.93677186E+01;
    COFD[1221] = 5.02567894E+00;
    COFD[1222] = -4.32045169E-01;
    COFD[1223] = 1.85132214E-02;
    COFD[1224] = -1.87476063E+01;
    COFD[1225] = 4.79683898E+00;
    COFD[1226] = -4.04829719E-01;
    COFD[1227] = 1.74325475E-02;
    COFD[1228] = -1.87476063E+01;
    COFD[1229] = 4.79683898E+00;
    COFD[1230] = -4.04829719E-01;
    COFD[1231] = 1.74325475E-02;
    COFD[1232] = -1.92654138E+01;
    COFD[1233] = 4.98286777E+00;
    COFD[1234] = -4.26970814E-01;
    COFD[1235] = 1.83122917E-02;
    COFD[1236] = -1.64758697E+01;
    COFD[1237] = 3.95035840E+00;
    COFD[1238] = -3.00959418E-01;
    COFD[1239] = 1.31692593E-02;
    COFD[1240] = -1.64825368E+01;
    COFD[1241] = 3.95035840E+00;
    COFD[1242] = -3.00959418E-01;
    COFD[1243] = 1.31692593E-02;
    COFD[1244] = -1.64887871E+01;
    COFD[1245] = 3.95035840E+00;
    COFD[1246] = -3.00959418E-01;
    COFD[1247] = 1.31692593E-02;
    COFD[1248] = -1.74111692E+01;
    COFD[1249] = 4.29676909E+00;
    COFD[1250] = -3.44085306E-01;
    COFD[1251] = 1.49671135E-02;
    COFD[1252] = -1.71808106E+01;
    COFD[1253] = 4.17889917E+00;
    COFD[1254] = -3.29752510E-01;
    COFD[1255] = 1.43850275E-02;
    COFD[1256] = -1.71860230E+01;
    COFD[1257] = 4.17889917E+00;
    COFD[1258] = -3.29752510E-01;
    COFD[1259] = 1.43850275E-02;
    COFD[1260] = -1.48738066E+01;
    COFD[1261] = 3.52327209E+00;
    COFD[1262] = -2.46286208E-01;
    COFD[1263] = 1.08285963E-02;
    COFD[1264] = -1.90116191E+01;
    COFD[1265] = 4.84384483E+00;
    COFD[1266] = -4.10265575E-01;
    COFD[1267] = 1.76414287E-02;
    COFD[1268] = -1.90116191E+01;
    COFD[1269] = 4.84384483E+00;
    COFD[1270] = -4.10265575E-01;
    COFD[1271] = 1.76414287E-02;
    COFD[1272] = -1.40076852E+01;
    COFD[1273] = 3.07549274E+00;
    COFD[1274] = -1.88889344E-01;
    COFD[1275] = 8.37152866E-03;
    COFD[1276] = -1.48853569E+01;
    COFD[1277] = 3.41988961E+00;
    COFD[1278] = -2.33128386E-01;
    COFD[1279] = 1.02689994E-02;
    COFD[1280] = -1.25098960E+01;
    COFD[1281] = 2.77873601E+00;
    COFD[1282] = -1.50637360E-01;
    COFD[1283] = 6.72684281E-03;
    COFD[1284] = -1.57972369E+01;
    COFD[1285] = 4.22225052E+00;
    COFD[1286] = -3.35156428E-01;
    COFD[1287] = 1.46104855E-02;
    COFD[1288] = -1.50584249E+01;
    COFD[1289] = 3.47945612E+00;
    COFD[1290] = -2.40703722E-01;
    COFD[1291] = 1.05907441E-02;
    COFD[1292] = -1.63254691E+01;
    COFD[1293] = 3.82388595E+00;
    COFD[1294] = -2.84480724E-01;
    COFD[1295] = 1.24506311E-02;
    COFD[1296] = -1.50724636E+01;
    COFD[1297] = 3.47945612E+00;
    COFD[1298] = -2.40703722E-01;
    COFD[1299] = 1.05907441E-02;
    COFD[1300] = -2.12639214E+01;
    COFD[1301] = 5.61184117E+00;
    COFD[1302] = -4.90532156E-01;
    COFD[1303] = 2.03507922E-02;
    COFD[1304] = -1.63301444E+01;
    COFD[1305] = 3.82388595E+00;
    COFD[1306] = -2.84480724E-01;
    COFD[1307] = 1.24506311E-02;
    COFD[1308] = -1.63345829E+01;
    COFD[1309] = 3.82388595E+00;
    COFD[1310] = -2.84480724E-01;
    COFD[1311] = 1.24506311E-02;
    COFD[1312] = -1.47558850E+01;
    COFD[1313] = 3.33113524E+00;
    COFD[1314] = -2.21479057E-01;
    COFD[1315] = 9.75837737E-03;
    COFD[1316] = -1.50076254E+01;
    COFD[1317] = 3.47945612E+00;
    COFD[1318] = -2.40703722E-01;
    COFD[1319] = 1.05907441E-02;
    COFD[1320] = -1.73027557E+01;
    COFD[1321] = 4.21416723E+00;
    COFD[1322] = -3.34163932E-01;
    COFD[1323] = 1.45697432E-02;
    COFD[1324] = -1.73027557E+01;
    COFD[1325] = 4.21416723E+00;
    COFD[1326] = -3.34163932E-01;
    COFD[1327] = 1.45697432E-02;
    COFD[1328] = -1.73198034E+01;
    COFD[1329] = 4.21416723E+00;
    COFD[1330] = -3.34163932E-01;
    COFD[1331] = 1.45697432E-02;
    COFD[1332] = -1.72556729E+01;
    COFD[1333] = 4.19029808E+00;
    COFD[1334] = -3.31177076E-01;
    COFD[1335] = 1.44446234E-02;
    COFD[1336] = -1.59634533E+01;
    COFD[1337] = 3.67388294E+00;
    COFD[1338] = -2.64990709E-01;
    COFD[1339] = 1.16042706E-02;
    COFD[1340] = -1.93015555E+01;
    COFD[1341] = 4.85015581E+00;
    COFD[1342] = -4.10945109E-01;
    COFD[1343] = 1.76651398E-02;
    COFD[1344] = -2.14160703E+01;
    COFD[1345] = 5.56531152E+00;
    COFD[1346] = -4.88789821E-01;
    COFD[1347] = 2.04437116E-02;
    COFD[1348] = -2.14215700E+01;
    COFD[1349] = 5.56531152E+00;
    COFD[1350] = -4.88789821E-01;
    COFD[1351] = 2.04437116E-02;
    COFD[1352] = -2.09376196E+01;
    COFD[1353] = 5.40870099E+00;
    COFD[1354] = -4.73017610E-01;
    COFD[1355] = 1.99399066E-02;
    COFD[1356] = -2.09376196E+01;
    COFD[1357] = 5.40870099E+00;
    COFD[1358] = -4.73017610E-01;
    COFD[1359] = 1.99399066E-02;
    COFD[1360] = -2.13538553E+01;
    COFD[1361] = 5.54007827E+00;
    COFD[1362] = -4.86434511E-01;
    COFD[1363] = 2.03779006E-02;
    COFD[1364] = -1.88171077E+01;
    COFD[1365] = 4.68393046E+00;
    COFD[1366] = -3.91610863E-01;
    COFD[1367] = 1.69174645E-02;
    COFD[1368] = -1.88241079E+01;
    COFD[1369] = 4.68393046E+00;
    COFD[1370] = -3.91610863E-01;
    COFD[1371] = 1.69174645E-02;
    COFD[1372] = -1.88306747E+01;
    COFD[1373] = 4.68393046E+00;
    COFD[1374] = -3.91610863E-01;
    COFD[1375] = 1.69174645E-02;
    COFD[1376] = -1.98418115E+01;
    COFD[1377] = 5.04367502E+00;
    COFD[1378] = -4.34153325E-01;
    COFD[1379] = 1.85956055E-02;
    COFD[1380] = -1.95263312E+01;
    COFD[1381] = 4.90255048E+00;
    COFD[1382] = -4.17368501E-01;
    COFD[1383] = 1.79287358E-02;
    COFD[1384] = -1.95318173E+01;
    COFD[1385] = 4.90255048E+00;
    COFD[1386] = -4.17368501E-01;
    COFD[1387] = 1.79287358E-02;
    COFD[1388] = -1.72572042E+01;
    COFD[1389] = 4.26063341E+00;
    COFD[1390] = -3.39848064E-01;
    COFD[1391] = 1.48021313E-02;
    COFD[1392] = -2.11349086E+01;
    COFD[1393] = 5.42846112E+00;
    COFD[1394] = -4.74321870E-01;
    COFD[1395] = 1.99459749E-02;
    COFD[1396] = -2.11349086E+01;
    COFD[1397] = 5.42846112E+00;
    COFD[1398] = -4.74321870E-01;
    COFD[1399] = 1.99459749E-02;
    COFD[1400] = -1.59404882E+01;
    COFD[1401] = 3.66853818E+00;
    COFD[1402] = -2.64346221E-01;
    COFD[1403] = 1.15784613E-02;
    COFD[1404] = -1.71942502E+01;
    COFD[1405] = 4.14993355E+00;
    COFD[1406] = -3.26168062E-01;
    COFD[1407] = 1.42364115E-02;
    COFD[1408] = -1.25098960E+01;
    COFD[1409] = 2.77873601E+00;
    COFD[1410] = -1.50637360E-01;
    COFD[1411] = 6.72684281E-03;
    COFD[1412] = -1.57972369E+01;
    COFD[1413] = 4.22225052E+00;
    COFD[1414] = -3.35156428E-01;
    COFD[1415] = 1.46104855E-02;
    COFD[1416] = -1.50584249E+01;
    COFD[1417] = 3.47945612E+00;
    COFD[1418] = -2.40703722E-01;
    COFD[1419] = 1.05907441E-02;
    COFD[1420] = -1.63254691E+01;
    COFD[1421] = 3.82388595E+00;
    COFD[1422] = -2.84480724E-01;
    COFD[1423] = 1.24506311E-02;
    COFD[1424] = -1.50724636E+01;
    COFD[1425] = 3.47945612E+00;
    COFD[1426] = -2.40703722E-01;
    COFD[1427] = 1.05907441E-02;
    COFD[1428] = -2.12639214E+01;
    COFD[1429] = 5.61184117E+00;
    COFD[1430] = -4.90532156E-01;
    COFD[1431] = 2.03507922E-02;
    COFD[1432] = -1.63301444E+01;
    COFD[1433] = 3.82388595E+00;
    COFD[1434] = -2.84480724E-01;
    COFD[1435] = 1.24506311E-02;
    COFD[1436] = -1.63345829E+01;
    COFD[1437] = 3.82388595E+00;
    COFD[1438] = -2.84480724E-01;
    COFD[1439] = 1.24506311E-02;
    COFD[1440] = -1.47558850E+01;
    COFD[1441] = 3.33113524E+00;
    COFD[1442] = -2.21479057E-01;
    COFD[1443] = 9.75837737E-03;
    COFD[1444] = -1.50076254E+01;
    COFD[1445] = 3.47945612E+00;
    COFD[1446] = -2.40703722E-01;
    COFD[1447] = 1.05907441E-02;
    COFD[1448] = -1.73027557E+01;
    COFD[1449] = 4.21416723E+00;
    COFD[1450] = -3.34163932E-01;
    COFD[1451] = 1.45697432E-02;
    COFD[1452] = -1.73027557E+01;
    COFD[1453] = 4.21416723E+00;
    COFD[1454] = -3.34163932E-01;
    COFD[1455] = 1.45697432E-02;
    COFD[1456] = -1.73198034E+01;
    COFD[1457] = 4.21416723E+00;
    COFD[1458] = -3.34163932E-01;
    COFD[1459] = 1.45697432E-02;
    COFD[1460] = -1.72556729E+01;
    COFD[1461] = 4.19029808E+00;
    COFD[1462] = -3.31177076E-01;
    COFD[1463] = 1.44446234E-02;
    COFD[1464] = -1.59634533E+01;
    COFD[1465] = 3.67388294E+00;
    COFD[1466] = -2.64990709E-01;
    COFD[1467] = 1.16042706E-02;
    COFD[1468] = -1.93015555E+01;
    COFD[1469] = 4.85015581E+00;
    COFD[1470] = -4.10945109E-01;
    COFD[1471] = 1.76651398E-02;
    COFD[1472] = -2.14160703E+01;
    COFD[1473] = 5.56531152E+00;
    COFD[1474] = -4.88789821E-01;
    COFD[1475] = 2.04437116E-02;
    COFD[1476] = -2.14215700E+01;
    COFD[1477] = 5.56531152E+00;
    COFD[1478] = -4.88789821E-01;
    COFD[1479] = 2.04437116E-02;
    COFD[1480] = -2.09376196E+01;
    COFD[1481] = 5.40870099E+00;
    COFD[1482] = -4.73017610E-01;
    COFD[1483] = 1.99399066E-02;
    COFD[1484] = -2.09376196E+01;
    COFD[1485] = 5.40870099E+00;
    COFD[1486] = -4.73017610E-01;
    COFD[1487] = 1.99399066E-02;
    COFD[1488] = -2.13538553E+01;
    COFD[1489] = 5.54007827E+00;
    COFD[1490] = -4.86434511E-01;
    COFD[1491] = 2.03779006E-02;
    COFD[1492] = -1.88171077E+01;
    COFD[1493] = 4.68393046E+00;
    COFD[1494] = -3.91610863E-01;
    COFD[1495] = 1.69174645E-02;
    COFD[1496] = -1.88241079E+01;
    COFD[1497] = 4.68393046E+00;
    COFD[1498] = -3.91610863E-01;
    COFD[1499] = 1.69174645E-02;
    COFD[1500] = -1.88306747E+01;
    COFD[1501] = 4.68393046E+00;
    COFD[1502] = -3.91610863E-01;
    COFD[1503] = 1.69174645E-02;
    COFD[1504] = -1.98418115E+01;
    COFD[1505] = 5.04367502E+00;
    COFD[1506] = -4.34153325E-01;
    COFD[1507] = 1.85956055E-02;
    COFD[1508] = -1.95263312E+01;
    COFD[1509] = 4.90255048E+00;
    COFD[1510] = -4.17368501E-01;
    COFD[1511] = 1.79287358E-02;
    COFD[1512] = -1.95318173E+01;
    COFD[1513] = 4.90255048E+00;
    COFD[1514] = -4.17368501E-01;
    COFD[1515] = 1.79287358E-02;
    COFD[1516] = -1.72572042E+01;
    COFD[1517] = 4.26063341E+00;
    COFD[1518] = -3.39848064E-01;
    COFD[1519] = 1.48021313E-02;
    COFD[1520] = -2.11349086E+01;
    COFD[1521] = 5.42846112E+00;
    COFD[1522] = -4.74321870E-01;
    COFD[1523] = 1.99459749E-02;
    COFD[1524] = -2.11349086E+01;
    COFD[1525] = 5.42846112E+00;
    COFD[1526] = -4.74321870E-01;
    COFD[1527] = 1.99459749E-02;
    COFD[1528] = -1.59404882E+01;
    COFD[1529] = 3.66853818E+00;
    COFD[1530] = -2.64346221E-01;
    COFD[1531] = 1.15784613E-02;
    COFD[1532] = -1.71942502E+01;
    COFD[1533] = 4.14993355E+00;
    COFD[1534] = -3.26168062E-01;
    COFD[1535] = 1.42364115E-02;
    COFD[1536] = -1.25141260E+01;
    COFD[1537] = 2.77873601E+00;
    COFD[1538] = -1.50637360E-01;
    COFD[1539] = 6.72684281E-03;
    COFD[1540] = -1.57994893E+01;
    COFD[1541] = 4.22225052E+00;
    COFD[1542] = -3.35156428E-01;
    COFD[1543] = 1.46104855E-02;
    COFD[1544] = -1.50766130E+01;
    COFD[1545] = 3.47945612E+00;
    COFD[1546] = -2.40703722E-01;
    COFD[1547] = 1.05907441E-02;
    COFD[1548] = -1.63493345E+01;
    COFD[1549] = 3.82388595E+00;
    COFD[1550] = -2.84480724E-01;
    COFD[1551] = 1.24506311E-02;
    COFD[1552] = -1.50911794E+01;
    COFD[1553] = 3.47945612E+00;
    COFD[1554] = -2.40703722E-01;
    COFD[1555] = 1.05907441E-02;
    COFD[1556] = -2.12831323E+01;
    COFD[1557] = 5.61184117E+00;
    COFD[1558] = -4.90532156E-01;
    COFD[1559] = 2.03507922E-02;
    COFD[1560] = -1.63542394E+01;
    COFD[1561] = 3.82388595E+00;
    COFD[1562] = -2.84480724E-01;
    COFD[1563] = 1.24506311E-02;
    COFD[1564] = -1.63588981E+01;
    COFD[1565] = 3.82388595E+00;
    COFD[1566] = -2.84480724E-01;
    COFD[1567] = 1.24506311E-02;
    COFD[1568] = -1.47715918E+01;
    COFD[1569] = 3.33113524E+00;
    COFD[1570] = -2.21479057E-01;
    COFD[1571] = 9.75837737E-03;
    COFD[1572] = -1.50240272E+01;
    COFD[1573] = 3.47945612E+00;
    COFD[1574] = -2.40703722E-01;
    COFD[1575] = 1.05907441E-02;
    COFD[1576] = -1.73198034E+01;
    COFD[1577] = 4.21416723E+00;
    COFD[1578] = -3.34163932E-01;
    COFD[1579] = 1.45697432E-02;
    COFD[1580] = -1.73198034E+01;
    COFD[1581] = 4.21416723E+00;
    COFD[1582] = -3.34163932E-01;
    COFD[1583] = 1.45697432E-02;
    COFD[1584] = -1.73374529E+01;
    COFD[1585] = 4.21416723E+00;
    COFD[1586] = -3.34163932E-01;
    COFD[1587] = 1.45697432E-02;
    COFD[1588] = -1.72738845E+01;
    COFD[1589] = 4.19029808E+00;
    COFD[1590] = -3.31177076E-01;
    COFD[1591] = 1.44446234E-02;
    COFD[1592] = -1.59863030E+01;
    COFD[1593] = 3.67388294E+00;
    COFD[1594] = -2.64990709E-01;
    COFD[1595] = 1.16042706E-02;
    COFD[1596] = -1.93276434E+01;
    COFD[1597] = 4.85015581E+00;
    COFD[1598] = -4.10945109E-01;
    COFD[1599] = 1.76651398E-02;
    COFD[1600] = -2.14391943E+01;
    COFD[1601] = 5.56531152E+00;
    COFD[1602] = -4.88789821E-01;
    COFD[1603] = 2.04437116E-02;
    COFD[1604] = -2.14449559E+01;
    COFD[1605] = 5.56531152E+00;
    COFD[1606] = -4.88789821E-01;
    COFD[1607] = 2.04437116E-02;
    COFD[1608] = -2.09612557E+01;
    COFD[1609] = 5.40870099E+00;
    COFD[1610] = -4.73017610E-01;
    COFD[1611] = 1.99399066E-02;
    COFD[1612] = -2.09612557E+01;
    COFD[1613] = 5.40870099E+00;
    COFD[1614] = -4.73017610E-01;
    COFD[1615] = 1.99399066E-02;
    COFD[1616] = -2.13777308E+01;
    COFD[1617] = 5.54007827E+00;
    COFD[1618] = -4.86434511E-01;
    COFD[1619] = 2.03779006E-02;
    COFD[1620] = -1.88390649E+01;
    COFD[1621] = 4.68393046E+00;
    COFD[1622] = -3.91610863E-01;
    COFD[1623] = 1.69174645E-02;
    COFD[1624] = -1.88463816E+01;
    COFD[1625] = 4.68393046E+00;
    COFD[1626] = -3.91610863E-01;
    COFD[1627] = 1.69174645E-02;
    COFD[1628] = -1.88532497E+01;
    COFD[1629] = 4.68393046E+00;
    COFD[1630] = -3.91610863E-01;
    COFD[1631] = 1.69174645E-02;
    COFD[1632] = -1.98646734E+01;
    COFD[1633] = 5.04367502E+00;
    COFD[1634] = -4.34153325E-01;
    COFD[1635] = 1.85956055E-02;
    COFD[1636] = -1.95494668E+01;
    COFD[1637] = 4.90255048E+00;
    COFD[1638] = -4.17368501E-01;
    COFD[1639] = 1.79287358E-02;
    COFD[1640] = -1.95552142E+01;
    COFD[1641] = 4.90255048E+00;
    COFD[1642] = -4.17368501E-01;
    COFD[1643] = 1.79287358E-02;
    COFD[1644] = -1.72828302E+01;
    COFD[1645] = 4.26063341E+00;
    COFD[1646] = -3.39848064E-01;
    COFD[1647] = 1.48021313E-02;
    COFD[1648] = -2.11606963E+01;
    COFD[1649] = 5.42846112E+00;
    COFD[1650] = -4.74321870E-01;
    COFD[1651] = 1.99459749E-02;
    COFD[1652] = -2.11606963E+01;
    COFD[1653] = 5.42846112E+00;
    COFD[1654] = -4.74321870E-01;
    COFD[1655] = 1.99459749E-02;
    COFD[1656] = -1.59633387E+01;
    COFD[1657] = 3.66853818E+00;
    COFD[1658] = -2.64346221E-01;
    COFD[1659] = 1.15784613E-02;
    COFD[1660] = -1.72196961E+01;
    COFD[1661] = 4.14993355E+00;
    COFD[1662] = -3.26168062E-01;
    COFD[1663] = 1.42364115E-02;
    COFD[1664] = -1.24693568E+01;
    COFD[1665] = 2.76686648E+00;
    COFD[1666] = -1.49120141E-01;
    COFD[1667] = 6.66220432E-03;
    COFD[1668] = -1.57199037E+01;
    COFD[1669] = 4.19936335E+00;
    COFD[1670] = -3.32311009E-01;
    COFD[1671] = 1.44921003E-02;
    COFD[1672] = -1.50270339E+01;
    COFD[1673] = 3.46140064E+00;
    COFD[1674] = -2.38440092E-01;
    COFD[1675] = 1.04960087E-02;
    COFD[1676] = -1.62724462E+01;
    COFD[1677] = 3.79163564E+00;
    COFD[1678] = -2.80257365E-01;
    COFD[1679] = 1.22656902E-02;
    COFD[1680] = -1.50420953E+01;
    COFD[1681] = 3.46140064E+00;
    COFD[1682] = -2.38440092E-01;
    COFD[1683] = 1.04960087E-02;
    COFD[1684] = -2.14087397E+01;
    COFD[1685] = 5.57282008E+00;
    COFD[1686] = -4.76690890E-01;
    COFD[1687] = 1.94000719E-02;
    COFD[1688] = -1.62775714E+01;
    COFD[1689] = 3.79163564E+00;
    COFD[1690] = -2.80257365E-01;
    COFD[1691] = 1.22656902E-02;
    COFD[1692] = -1.62824412E+01;
    COFD[1693] = 3.79163564E+00;
    COFD[1694] = -2.80257365E-01;
    COFD[1695] = 1.22656902E-02;
    COFD[1696] = -1.47048024E+01;
    COFD[1697] = 3.30594991E+00;
    COFD[1698] = -2.18182207E-01;
    COFD[1699] = 9.61429447E-03;
    COFD[1700] = -1.49727799E+01;
    COFD[1701] = 3.46140064E+00;
    COFD[1702] = -2.38440092E-01;
    COFD[1703] = 1.04960087E-02;
    COFD[1704] = -1.72556729E+01;
    COFD[1705] = 4.19029808E+00;
    COFD[1706] = -3.31177076E-01;
    COFD[1707] = 1.44446234E-02;
    COFD[1708] = -1.72556729E+01;
    COFD[1709] = 4.19029808E+00;
    COFD[1710] = -3.31177076E-01;
    COFD[1711] = 1.44446234E-02;
    COFD[1712] = -1.72738845E+01;
    COFD[1713] = 4.19029808E+00;
    COFD[1714] = -3.31177076E-01;
    COFD[1715] = 1.44446234E-02;
    COFD[1716] = -1.72167708E+01;
    COFD[1717] = 4.16886779E+00;
    COFD[1718] = -3.28518156E-01;
    COFD[1719] = 1.43341626E-02;
    COFD[1720] = -1.59525102E+01;
    COFD[1721] = 3.66023858E+00;
    COFD[1722] = -2.63401043E-01;
    COFD[1723] = 1.15432000E-02;
    COFD[1724] = -1.92867554E+01;
    COFD[1725] = 4.83375900E+00;
    COFD[1726] = -4.09146560E-01;
    COFD[1727] = 1.76006599E-02;
    COFD[1728] = -2.14022336E+01;
    COFD[1729] = 5.55346617E+00;
    COFD[1730] = -4.87783156E-01;
    COFD[1731] = 2.04210886E-02;
    COFD[1732] = -2.14082453E+01;
    COFD[1733] = 5.55346617E+00;
    COFD[1734] = -4.87783156E-01;
    COFD[1735] = 2.04210886E-02;
    COFD[1736] = -2.11381508E+01;
    COFD[1737] = 5.45574440E+00;
    COFD[1738] = -4.77436155E-01;
    COFD[1739] = 2.00644596E-02;
    COFD[1740] = -2.11381508E+01;
    COFD[1741] = 5.45574440E+00;
    COFD[1742] = -4.77436155E-01;
    COFD[1743] = 2.00644596E-02;
    COFD[1744] = -2.13319784E+01;
    COFD[1745] = 5.52422470E+00;
    COFD[1746] = -4.84872944E-01;
    COFD[1747] = 2.03298213E-02;
    COFD[1748] = -1.87821119E+01;
    COFD[1749] = 4.66162351E+00;
    COFD[1750] = -3.88920477E-01;
    COFD[1751] = 1.68089648E-02;
    COFD[1752] = -1.87897298E+01;
    COFD[1753] = 4.66162351E+00;
    COFD[1754] = -3.88920477E-01;
    COFD[1755] = 1.68089648E-02;
    COFD[1756] = -1.87968848E+01;
    COFD[1757] = 4.66162351E+00;
    COFD[1758] = -3.88920477E-01;
    COFD[1759] = 1.68089648E-02;
    COFD[1760] = -1.98075055E+01;
    COFD[1761] = 5.02169524E+00;
    COFD[1762] = -4.31582804E-01;
    COFD[1763] = 1.84953568E-02;
    COFD[1764] = -1.94763688E+01;
    COFD[1765] = 4.87333294E+00;
    COFD[1766] = -4.13769241E-01;
    COFD[1767] = 1.77802244E-02;
    COFD[1768] = -1.94823660E+01;
    COFD[1769] = 4.87333294E+00;
    COFD[1770] = -4.13769241E-01;
    COFD[1771] = 1.77802244E-02;
    COFD[1772] = -1.72316148E+01;
    COFD[1773] = 4.24011069E+00;
    COFD[1774] = -3.37339810E-01;
    COFD[1775] = 1.46996679E-02;
    COFD[1776] = -2.11309207E+01;
    COFD[1777] = 5.41773516E+00;
    COFD[1778] = -4.73414338E-01;
    COFD[1779] = 1.99258685E-02;
    COFD[1780] = -2.11309207E+01;
    COFD[1781] = 5.41773516E+00;
    COFD[1782] = -4.73414338E-01;
    COFD[1783] = 1.99258685E-02;
    COFD[1784] = -1.59327297E+01;
    COFD[1785] = 3.65620899E+00;
    COFD[1786] = -2.62933804E-01;
    COFD[1787] = 1.15253223E-02;
    COFD[1788] = -1.71754154E+01;
    COFD[1789] = 4.13131681E+00;
    COFD[1790] = -3.23897559E-01;
    COFD[1791] = 1.41438222E-02;
    COFD[1792] = -1.17159737E+01;
    COFD[1793] = 2.48123210E+00;
    COFD[1794] = -1.11322604E-01;
    COFD[1795] = 4.99282389E-03;
    COFD[1796] = -1.43151174E+01;
    COFD[1797] = 3.68038508E+00;
    COFD[1798] = -2.65779346E-01;
    COFD[1799] = 1.16360771E-02;
    COFD[1800] = -1.40999008E+01;
    COFD[1801] = 3.08120012E+00;
    COFD[1802] = -1.89629903E-01;
    COFD[1803] = 8.40361952E-03;
    COFD[1804] = -1.52721107E+01;
    COFD[1805] = 3.36790500E+00;
    COFD[1806] = -2.26321740E-01;
    COFD[1807] = 9.97135055E-03;
    COFD[1808] = -1.41191261E+01;
    COFD[1809] = 3.08120012E+00;
    COFD[1810] = -1.89629903E-01;
    COFD[1811] = 8.40361952E-03;
    COFD[1812] = -2.11388331E+01;
    COFD[1813] = 5.55529675E+00;
    COFD[1814] = -4.87942518E-01;
    COFD[1815] = 2.04249054E-02;
    COFD[1816] = -1.52792891E+01;
    COFD[1817] = 3.36790500E+00;
    COFD[1818] = -2.26321740E-01;
    COFD[1819] = 9.97135055E-03;
    COFD[1820] = -1.52861376E+01;
    COFD[1821] = 3.36790500E+00;
    COFD[1822] = -2.26321740E-01;
    COFD[1823] = 9.97135055E-03;
    COFD[1824] = -1.38862839E+01;
    COFD[1825] = 2.97564184E+00;
    COFD[1826] = -1.76025309E-01;
    COFD[1827] = 7.81869993E-03;
    COFD[1828] = -1.40318948E+01;
    COFD[1829] = 3.08120012E+00;
    COFD[1830] = -1.89629903E-01;
    COFD[1831] = 8.40361952E-03;
    COFD[1832] = -1.59634533E+01;
    COFD[1833] = 3.67388294E+00;
    COFD[1834] = -2.64990709E-01;
    COFD[1835] = 1.16042706E-02;
    COFD[1836] = -1.59634533E+01;
    COFD[1837] = 3.67388294E+00;
    COFD[1838] = -2.64990709E-01;
    COFD[1839] = 1.16042706E-02;
    COFD[1840] = -1.59863030E+01;
    COFD[1841] = 3.67388294E+00;
    COFD[1842] = -2.64990709E-01;
    COFD[1843] = 1.16042706E-02;
    COFD[1844] = -1.59525102E+01;
    COFD[1845] = 3.66023858E+00;
    COFD[1846] = -2.63401043E-01;
    COFD[1847] = 1.15432000E-02;
    COFD[1848] = -1.50233475E+01;
    COFD[1849] = 3.26660767E+00;
    COFD[1850] = -2.13287177E-01;
    COFD[1851] = 9.41137857E-03;
    COFD[1852] = -1.81735763E+01;
    COFD[1853] = 4.38391495E+00;
    COFD[1854] = -3.54941287E-01;
    COFD[1855] = 1.54195107E-02;
    COFD[1856] = -2.05045578E+01;
    COFD[1857] = 5.23843909E+00;
    COFD[1858] = -4.55815614E-01;
    COFD[1859] = 1.93898040E-02;
    COFD[1860] = -2.05128705E+01;
    COFD[1861] = 5.23843909E+00;
    COFD[1862] = -4.55815614E-01;
    COFD[1863] = 1.93898040E-02;
    COFD[1864] = -2.02642227E+01;
    COFD[1865] = 5.14499740E+00;
    COFD[1866] = -4.45694430E-01;
    COFD[1867] = 1.90318646E-02;
    COFD[1868] = -2.02642227E+01;
    COFD[1869] = 5.14499740E+00;
    COFD[1870] = -4.45694430E-01;
    COFD[1871] = 1.90318646E-02;
    COFD[1872] = -2.04144604E+01;
    COFD[1873] = 5.19614628E+00;
    COFD[1874] = -4.50889164E-01;
    COFD[1875] = 1.91983328E-02;
    COFD[1876] = -1.76182365E+01;
    COFD[1877] = 4.19935698E+00;
    COFD[1878] = -3.32310212E-01;
    COFD[1879] = 1.44920670E-02;
    COFD[1880] = -1.76285640E+01;
    COFD[1881] = 4.19935698E+00;
    COFD[1882] = -3.32310212E-01;
    COFD[1883] = 1.44920670E-02;
    COFD[1884] = -1.76383156E+01;
    COFD[1885] = 4.19935698E+00;
    COFD[1886] = -3.32310212E-01;
    COFD[1887] = 1.44920670E-02;
    COFD[1888] = -1.86157761E+01;
    COFD[1889] = 4.55689508E+00;
    COFD[1890] = -3.75937921E-01;
    COFD[1891] = 1.62703488E-02;
    COFD[1892] = -1.83455435E+01;
    COFD[1893] = 4.42828044E+00;
    COFD[1894] = -3.60417833E-01;
    COFD[1895] = 1.56455103E-02;
    COFD[1896] = -1.83538377E+01;
    COFD[1897] = 4.42828044E+00;
    COFD[1898] = -3.60417833E-01;
    COFD[1899] = 1.56455103E-02;
    COFD[1900] = -1.60261675E+01;
    COFD[1901] = 3.73312045E+00;
    COFD[1902] = -2.72579779E-01;
    COFD[1903] = 1.19290272E-02;
    COFD[1904] = -2.02922701E+01;
    COFD[1905] = 5.11106992E+00;
    COFD[1906] = -4.42047129E-01;
    COFD[1907] = 1.89042990E-02;
    COFD[1908] = -2.02922701E+01;
    COFD[1909] = 5.11106992E+00;
    COFD[1910] = -4.42047129E-01;
    COFD[1911] = 1.89042990E-02;
    COFD[1912] = -1.50031687E+01;
    COFD[1913] = 3.26223357E+00;
    COFD[1914] = -2.12746642E-01;
    COFD[1915] = 9.38912883E-03;
    COFD[1916] = -1.60074211E+01;
    COFD[1917] = 3.63723937E+00;
    COFD[1918] = -2.60754222E-01;
    COFD[1919] = 1.14428814E-02;
    COFD[1920] = -1.37794315E+01;
    COFD[1921] = 3.23973858E+00;
    COFD[1922] = -2.09989036E-01;
    COFD[1923] = 9.27667906E-03;
    COFD[1924] = -1.76147026E+01;
    COFD[1925] = 4.86049500E+00;
    COFD[1926] = -4.12200578E-01;
    COFD[1927] = 1.77160971E-02;
    COFD[1928] = -1.70534856E+01;
    COFD[1929] = 4.14240922E+00;
    COFD[1930] = -3.25239774E-01;
    COFD[1931] = 1.41980687E-02;
    COFD[1932] = -1.84688406E+01;
    COFD[1933] = 4.49330851E+00;
    COFD[1934] = -3.68208715E-01;
    COFD[1935] = 1.59565402E-02;
    COFD[1936] = -1.70757047E+01;
    COFD[1937] = 4.14240922E+00;
    COFD[1938] = -3.25239774E-01;
    COFD[1939] = 1.41980687E-02;
    COFD[1940] = -2.07653719E+01;
    COFD[1941] = 5.01092022E+00;
    COFD[1942] = -3.77985635E-01;
    COFD[1943] = 1.40968645E-02;
    COFD[1944] = -1.84777607E+01;
    COFD[1945] = 4.49330851E+00;
    COFD[1946] = -3.68208715E-01;
    COFD[1947] = 1.59565402E-02;
    COFD[1948] = -1.84863000E+01;
    COFD[1949] = 4.49330851E+00;
    COFD[1950] = -3.68208715E-01;
    COFD[1951] = 1.59565402E-02;
    COFD[1952] = -1.67390425E+01;
    COFD[1953] = 4.00828594E+00;
    COFD[1954] = -3.08414344E-01;
    COFD[1955] = 1.34907430E-02;
    COFD[1956] = -1.69758891E+01;
    COFD[1957] = 4.14240922E+00;
    COFD[1958] = -3.25239774E-01;
    COFD[1959] = 1.41980687E-02;
    COFD[1960] = -1.93015555E+01;
    COFD[1961] = 4.85015581E+00;
    COFD[1962] = -4.10945109E-01;
    COFD[1963] = 1.76651398E-02;
    COFD[1964] = -1.93015555E+01;
    COFD[1965] = 4.85015581E+00;
    COFD[1966] = -4.10945109E-01;
    COFD[1967] = 1.76651398E-02;
    COFD[1968] = -1.93276434E+01;
    COFD[1969] = 4.85015581E+00;
    COFD[1970] = -4.10945109E-01;
    COFD[1971] = 1.76651398E-02;
    COFD[1972] = -1.92867554E+01;
    COFD[1973] = 4.83375900E+00;
    COFD[1974] = -4.09146560E-01;
    COFD[1975] = 1.76006599E-02;
    COFD[1976] = -1.81735763E+01;
    COFD[1977] = 4.38391495E+00;
    COFD[1978] = -3.54941287E-01;
    COFD[1979] = 1.54195107E-02;
    COFD[1980] = -2.13425698E+01;
    COFD[1981] = 5.40460130E+00;
    COFD[1982] = -4.72718910E-01;
    COFD[1983] = 1.99362717E-02;
    COFD[1984] = -2.19215555E+01;
    COFD[1985] = 5.45216133E+00;
    COFD[1986] = -4.52916925E-01;
    COFD[1987] = 1.80456400E-02;
    COFD[1988] = -2.19317743E+01;
    COFD[1989] = 5.45216133E+00;
    COFD[1990] = -4.52916925E-01;
    COFD[1991] = 1.80456400E-02;
    COFD[1992] = -2.20421041E+01;
    COFD[1993] = 5.52708332E+00;
    COFD[1994] = -4.68000808E-01;
    COFD[1995] = 1.89131908E-02;
    COFD[1996] = -2.20421041E+01;
    COFD[1997] = 5.52708332E+00;
    COFD[1998] = -4.68000808E-01;
    COFD[1999] = 1.89131908E-02;
    COFD[2000] = -2.20063594E+01;
    COFD[2001] = 5.48540187E+00;
    COFD[2002] = -4.58962148E-01;
    COFD[2003] = 1.83770355E-02;
    COFD[2004] = -2.09066354E+01;
    COFD[2005] = 5.30153901E+00;
    COFD[2006] = -4.63335119E-01;
    COFD[2007] = 1.96897053E-02;
    COFD[2008] = -2.09191285E+01;
    COFD[2009] = 5.30153901E+00;
    COFD[2010] = -4.63335119E-01;
    COFD[2011] = 1.96897053E-02;
    COFD[2012] = -2.09309753E+01;
    COFD[2013] = 5.30153901E+00;
    COFD[2014] = -4.63335119E-01;
    COFD[2015] = 1.96897053E-02;
    COFD[2016] = -2.16802612E+01;
    COFD[2017] = 5.52918296E+00;
    COFD[2018] = -4.85360709E-01;
    COFD[2019] = 2.03448006E-02;
    COFD[2020] = -2.14224484E+01;
    COFD[2021] = 5.41729961E+00;
    COFD[2022] = -4.73400281E-01;
    COFD[2023] = 1.99269567E-02;
    COFD[2024] = -2.14326461E+01;
    COFD[2025] = 5.41729961E+00;
    COFD[2026] = -4.73400281E-01;
    COFD[2027] = 1.99269567E-02;
    COFD[2028] = -1.94485982E+01;
    COFD[2029] = 4.91446566E+00;
    COFD[2030] = -4.18837152E-01;
    COFD[2031] = 1.79893537E-02;
    COFD[2032] = -2.22116706E+01;
    COFD[2033] = 5.54251230E+00;
    COFD[2034] = -4.70946314E-01;
    COFD[2035] = 1.90785869E-02;
    COFD[2036] = -2.22116706E+01;
    COFD[2037] = 5.54251230E+00;
    COFD[2038] = -4.70946314E-01;
    COFD[2039] = 1.90785869E-02;
    COFD[2040] = -1.81432461E+01;
    COFD[2041] = 4.37565431E+00;
    COFD[2042] = -3.53906025E-01;
    COFD[2043] = 1.53760786E-02;
    COFD[2044] = -1.93483692E+01;
    COFD[2045] = 4.79506290E+00;
    COFD[2046] = -4.04621659E-01;
    COFD[2047] = 1.74244230E-02;
    COFD[2048] = -1.60517370E+01;
    COFD[2049] = 4.11188603E+00;
    COFD[2050] = -3.21540884E-01;
    COFD[2051] = 1.40482564E-02;
    COFD[2052] = -1.97544450E+01;
    COFD[2053] = 5.56931926E+00;
    COFD[2054] = -4.89105511E-01;
    COFD[2055] = 2.04493129E-02;
    COFD[2056] = -1.94313116E+01;
    COFD[2057] = 5.02567894E+00;
    COFD[2058] = -4.32045169E-01;
    COFD[2059] = 1.85132214E-02;
    COFD[2060] = -2.08204449E+01;
    COFD[2061] = 5.35267674E+00;
    COFD[2062] = -4.69010505E-01;
    COFD[2063] = 1.98979152E-02;
    COFD[2064] = -1.94507876E+01;
    COFD[2065] = 5.02567894E+00;
    COFD[2066] = -4.32045169E-01;
    COFD[2067] = 1.85132214E-02;
    COFD[2068] = -1.77498543E+01;
    COFD[2069] = 3.57475686E+00;
    COFD[2070] = -1.56396297E-01;
    COFD[2071] = 3.12157721E-03;
    COFD[2072] = -2.08277598E+01;
    COFD[2073] = 5.35267674E+00;
    COFD[2074] = -4.69010505E-01;
    COFD[2075] = 1.98979152E-02;
    COFD[2076] = -2.08347403E+01;
    COFD[2077] = 5.35267674E+00;
    COFD[2078] = -4.69010505E-01;
    COFD[2079] = 1.98979152E-02;
    COFD[2080] = -1.90526941E+01;
    COFD[2081] = 4.86821670E+00;
    COFD[2082] = -4.13144121E-01;
    COFD[2083] = 1.77546701E-02;
    COFD[2084] = -1.93624931E+01;
    COFD[2085] = 5.02567894E+00;
    COFD[2086] = -4.32045169E-01;
    COFD[2087] = 1.85132214E-02;
    COFD[2088] = -2.14160703E+01;
    COFD[2089] = 5.56531152E+00;
    COFD[2090] = -4.88789821E-01;
    COFD[2091] = 2.04437116E-02;
    COFD[2092] = -2.14160703E+01;
    COFD[2093] = 5.56531152E+00;
    COFD[2094] = -4.88789821E-01;
    COFD[2095] = 2.04437116E-02;
    COFD[2096] = -2.14391943E+01;
    COFD[2097] = 5.56531152E+00;
    COFD[2098] = -4.88789821E-01;
    COFD[2099] = 2.04437116E-02;
    COFD[2100] = -2.14022336E+01;
    COFD[2101] = 5.55346617E+00;
    COFD[2102] = -4.87783156E-01;
    COFD[2103] = 2.04210886E-02;
    COFD[2104] = -2.05045578E+01;
    COFD[2105] = 5.23843909E+00;
    COFD[2106] = -4.55815614E-01;
    COFD[2107] = 1.93898040E-02;
    COFD[2108] = -2.19215555E+01;
    COFD[2109] = 5.45216133E+00;
    COFD[2110] = -4.52916925E-01;
    COFD[2111] = 1.80456400E-02;
    COFD[2112] = -1.90328712E+01;
    COFD[2113] = 3.99221757E+00;
    COFD[2114] = -2.19854880E-01;
    COFD[2115] = 6.22736279E-03;
    COFD[2116] = -1.90413348E+01;
    COFD[2117] = 3.99221757E+00;
    COFD[2118] = -2.19854880E-01;
    COFD[2119] = 6.22736279E-03;
    COFD[2120] = -2.01801667E+01;
    COFD[2121] = 4.53183330E+00;
    COFD[2122] = -3.02186760E-01;
    COFD[2123] = 1.02756490E-02;
    COFD[2124] = -2.01801667E+01;
    COFD[2125] = 4.53183330E+00;
    COFD[2126] = -3.02186760E-01;
    COFD[2127] = 1.02756490E-02;
    COFD[2128] = -1.93125662E+01;
    COFD[2129] = 4.10954793E+00;
    COFD[2130] = -2.37523329E-01;
    COFD[2131] = 7.08858141E-03;
    COFD[2132] = -2.19851338E+01;
    COFD[2133] = 5.55935694E+00;
    COFD[2134] = -4.74154740E-01;
    COFD[2135] = 1.92584304E-02;
    COFD[2136] = -2.19956352E+01;
    COFD[2137] = 5.55935694E+00;
    COFD[2138] = -4.74154740E-01;
    COFD[2139] = 1.92584304E-02;
    COFD[2140] = -2.20055544E+01;
    COFD[2141] = 5.55935694E+00;
    COFD[2142] = -4.74154740E-01;
    COFD[2143] = 1.92584304E-02;
    COFD[2144] = -2.16296373E+01;
    COFD[2145] = 5.29019717E+00;
    COFD[2146] = -4.24502606E-01;
    COFD[2147] = 1.65197343E-02;
    COFD[2148] = -2.19229190E+01;
    COFD[2149] = 5.41841631E+00;
    COFD[2150] = -4.46818971E-01;
    COFD[2151] = 1.77127652E-02;
    COFD[2152] = -2.19313638E+01;
    COFD[2153] = 5.41841631E+00;
    COFD[2154] = -4.46818971E-01;
    COFD[2155] = 1.77127652E-02;
    COFD[2156] = -2.14204185E+01;
    COFD[2157] = 5.59268435E+00;
    COFD[2158] = -4.91232974E-01;
    COFD[2159] = 2.05064746E-02;
    COFD[2160] = -2.00915040E+01;
    COFD[2161] = 4.41511629E+00;
    COFD[2162] = -2.84086963E-01;
    COFD[2163] = 9.37586971E-03;
    COFD[2164] = -2.00915040E+01;
    COFD[2165] = 4.41511629E+00;
    COFD[2166] = -2.84086963E-01;
    COFD[2167] = 9.37586971E-03;
    COFD[2168] = -2.04750581E+01;
    COFD[2169] = 5.23112374E+00;
    COFD[2170] = -4.54967682E-01;
    COFD[2171] = 1.93570423E-02;
    COFD[2172] = -2.14255087E+01;
    COFD[2173] = 5.52240865E+00;
    COFD[2174] = -4.84699537E-01;
    COFD[2175] = 2.03247833E-02;
    COFD[2176] = -1.60528285E+01;
    COFD[2177] = 4.11188603E+00;
    COFD[2178] = -3.21540884E-01;
    COFD[2179] = 1.40482564E-02;
    COFD[2180] = -1.97550088E+01;
    COFD[2181] = 5.56931926E+00;
    COFD[2182] = -4.89105511E-01;
    COFD[2183] = 2.04493129E-02;
    COFD[2184] = -1.94373127E+01;
    COFD[2185] = 5.02567894E+00;
    COFD[2186] = -4.32045169E-01;
    COFD[2187] = 1.85132214E-02;
    COFD[2188] = -2.08293255E+01;
    COFD[2189] = 5.35267674E+00;
    COFD[2190] = -4.69010505E-01;
    COFD[2191] = 1.98979152E-02;
    COFD[2192] = -1.94570287E+01;
    COFD[2193] = 5.02567894E+00;
    COFD[2194] = -4.32045169E-01;
    COFD[2195] = 1.85132214E-02;
    COFD[2196] = -1.77563250E+01;
    COFD[2197] = 3.57475686E+00;
    COFD[2198] = -1.56396297E-01;
    COFD[2199] = 3.12157721E-03;
    COFD[2200] = -2.08367725E+01;
    COFD[2201] = 5.35267674E+00;
    COFD[2202] = -4.69010505E-01;
    COFD[2203] = 1.98979152E-02;
    COFD[2204] = -2.08438809E+01;
    COFD[2205] = 5.35267674E+00;
    COFD[2206] = -4.69010505E-01;
    COFD[2207] = 1.98979152E-02;
    COFD[2208] = -1.90576320E+01;
    COFD[2209] = 4.86821670E+00;
    COFD[2210] = -4.13144121E-01;
    COFD[2211] = 1.77546701E-02;
    COFD[2212] = -1.93677186E+01;
    COFD[2213] = 5.02567894E+00;
    COFD[2214] = -4.32045169E-01;
    COFD[2215] = 1.85132214E-02;
    COFD[2216] = -2.14215700E+01;
    COFD[2217] = 5.56531152E+00;
    COFD[2218] = -4.88789821E-01;
    COFD[2219] = 2.04437116E-02;
    COFD[2220] = -2.14215700E+01;
    COFD[2221] = 5.56531152E+00;
    COFD[2222] = -4.88789821E-01;
    COFD[2223] = 2.04437116E-02;
    COFD[2224] = -2.14449559E+01;
    COFD[2225] = 5.56531152E+00;
    COFD[2226] = -4.88789821E-01;
    COFD[2227] = 2.04437116E-02;
    COFD[2228] = -2.14082453E+01;
    COFD[2229] = 5.55346617E+00;
    COFD[2230] = -4.87783156E-01;
    COFD[2231] = 2.04210886E-02;
    COFD[2232] = -2.05128705E+01;
    COFD[2233] = 5.23843909E+00;
    COFD[2234] = -4.55815614E-01;
    COFD[2235] = 1.93898040E-02;
    COFD[2236] = -2.19317743E+01;
    COFD[2237] = 5.45216133E+00;
    COFD[2238] = -4.52916925E-01;
    COFD[2239] = 1.80456400E-02;
    COFD[2240] = -1.90413348E+01;
    COFD[2241] = 3.99221757E+00;
    COFD[2242] = -2.19854880E-01;
    COFD[2243] = 6.22736279E-03;
    COFD[2244] = -1.90499441E+01;
    COFD[2245] = 3.99221757E+00;
    COFD[2246] = -2.19854880E-01;
    COFD[2247] = 6.22736279E-03;
    COFD[2248] = -2.01889168E+01;
    COFD[2249] = 4.53183330E+00;
    COFD[2250] = -3.02186760E-01;
    COFD[2251] = 1.02756490E-02;
    COFD[2252] = -2.01889168E+01;
    COFD[2253] = 4.53183330E+00;
    COFD[2254] = -3.02186760E-01;
    COFD[2255] = 1.02756490E-02;
    COFD[2256] = -1.93214527E+01;
    COFD[2257] = 4.10954793E+00;
    COFD[2258] = -2.37523329E-01;
    COFD[2259] = 7.08858141E-03;
    COFD[2260] = -2.19929679E+01;
    COFD[2261] = 5.55935694E+00;
    COFD[2262] = -4.74154740E-01;
    COFD[2263] = 1.92584304E-02;
    COFD[2264] = -2.20036369E+01;
    COFD[2265] = 5.55935694E+00;
    COFD[2266] = -4.74154740E-01;
    COFD[2267] = 1.92584304E-02;
    COFD[2268] = -2.20137178E+01;
    COFD[2269] = 5.55935694E+00;
    COFD[2270] = -4.74154740E-01;
    COFD[2271] = 1.92584304E-02;
    COFD[2272] = -2.16379567E+01;
    COFD[2273] = 5.29019717E+00;
    COFD[2274] = -4.24502606E-01;
    COFD[2275] = 1.65197343E-02;
    COFD[2276] = -2.19313890E+01;
    COFD[2277] = 5.41841631E+00;
    COFD[2278] = -4.46818971E-01;
    COFD[2279] = 1.77127652E-02;
    COFD[2280] = -2.19399793E+01;
    COFD[2281] = 5.41841631E+00;
    COFD[2282] = -4.46818971E-01;
    COFD[2283] = 1.77127652E-02;
    COFD[2284] = -2.14303479E+01;
    COFD[2285] = 5.59268435E+00;
    COFD[2286] = -4.91232974E-01;
    COFD[2287] = 2.05064746E-02;
    COFD[2288] = -2.01015340E+01;
    COFD[2289] = 4.41511629E+00;
    COFD[2290] = -2.84086963E-01;
    COFD[2291] = 9.37586971E-03;
    COFD[2292] = -2.01015340E+01;
    COFD[2293] = 4.41511629E+00;
    COFD[2294] = -2.84086963E-01;
    COFD[2295] = 9.37586971E-03;
    COFD[2296] = -2.04833713E+01;
    COFD[2297] = 5.23112374E+00;
    COFD[2298] = -4.54967682E-01;
    COFD[2299] = 1.93570423E-02;
    COFD[2300] = -2.14353267E+01;
    COFD[2301] = 5.52240865E+00;
    COFD[2302] = -4.84699537E-01;
    COFD[2303] = 2.03247833E-02;
    COFD[2304] = -1.58456300E+01;
    COFD[2305] = 4.02074783E+00;
    COFD[2306] = -3.10018522E-01;
    COFD[2307] = 1.35599552E-02;
    COFD[2308] = -1.92718582E+01;
    COFD[2309] = 5.41172124E+00;
    COFD[2310] = -4.73213887E-01;
    COFD[2311] = 1.99405473E-02;
    COFD[2312] = -1.88179418E+01;
    COFD[2313] = 4.79683898E+00;
    COFD[2314] = -4.04829719E-01;
    COFD[2315] = 1.74325475E-02;
    COFD[2316] = -2.04928958E+01;
    COFD[2317] = 5.22397933E+00;
    COFD[2318] = -4.54138171E-01;
    COFD[2319] = 1.93249285E-02;
    COFD[2320] = -1.88378874E+01;
    COFD[2321] = 4.79683898E+00;
    COFD[2322] = -4.04829719E-01;
    COFD[2323] = 1.74325475E-02;
    COFD[2324] = -1.65295288E+01;
    COFD[2325] = 2.97569206E+00;
    COFD[2326] = -6.75652842E-02;
    COFD[2327] = -1.08648422E-03;
    COFD[2328] = -2.02637994E+01;
    COFD[2329] = 5.14984081E+00;
    COFD[2330] = -4.46093018E-01;
    COFD[2331] = 1.90396647E-02;
    COFD[2332] = -2.02710316E+01;
    COFD[2333] = 5.14984081E+00;
    COFD[2334] = -4.46093018E-01;
    COFD[2335] = 1.90396647E-02;
    COFD[2336] = -1.85157414E+01;
    COFD[2337] = 4.67076124E+00;
    COFD[2338] = -3.90022427E-01;
    COFD[2339] = 1.68533953E-02;
    COFD[2340] = -1.87476063E+01;
    COFD[2341] = 4.79683898E+00;
    COFD[2342] = -4.04829719E-01;
    COFD[2343] = 1.74325475E-02;
    COFD[2344] = -2.09376196E+01;
    COFD[2345] = 5.40870099E+00;
    COFD[2346] = -4.73017610E-01;
    COFD[2347] = 1.99399066E-02;
    COFD[2348] = -2.09376196E+01;
    COFD[2349] = 5.40870099E+00;
    COFD[2350] = -4.73017610E-01;
    COFD[2351] = 1.99399066E-02;
    COFD[2352] = -2.09612557E+01;
    COFD[2353] = 5.40870099E+00;
    COFD[2354] = -4.73017610E-01;
    COFD[2355] = 1.99399066E-02;
    COFD[2356] = -2.11381508E+01;
    COFD[2357] = 5.45574440E+00;
    COFD[2358] = -4.77436155E-01;
    COFD[2359] = 2.00644596E-02;
    COFD[2360] = -2.02642227E+01;
    COFD[2361] = 5.14499740E+00;
    COFD[2362] = -4.45694430E-01;
    COFD[2363] = 1.90318646E-02;
    COFD[2364] = -2.20421041E+01;
    COFD[2365] = 5.52708332E+00;
    COFD[2366] = -4.68000808E-01;
    COFD[2367] = 1.89131908E-02;
    COFD[2368] = -2.01801667E+01;
    COFD[2369] = 4.53183330E+00;
    COFD[2370] = -3.02186760E-01;
    COFD[2371] = 1.02756490E-02;
    COFD[2372] = -2.01889168E+01;
    COFD[2373] = 4.53183330E+00;
    COFD[2374] = -3.02186760E-01;
    COFD[2375] = 1.02756490E-02;
    COFD[2376] = -1.95877017E+01;
    COFD[2377] = 4.27643051E+00;
    COFD[2378] = -2.68040901E-01;
    COFD[2379] = 8.77650113E-03;
    COFD[2380] = -1.95877017E+01;
    COFD[2381] = 4.27643051E+00;
    COFD[2382] = -2.68040901E-01;
    COFD[2383] = 8.77650113E-03;
    COFD[2384] = -2.03599050E+01;
    COFD[2385] = 4.60682543E+00;
    COFD[2386] = -3.13971634E-01;
    COFD[2387] = 1.08661011E-02;
    COFD[2388] = -2.19383403E+01;
    COFD[2389] = 5.59157589E+00;
    COFD[2390] = -4.85617912E-01;
    COFD[2391] = 2.00461138E-02;
    COFD[2392] = -2.19491710E+01;
    COFD[2393] = 5.59157589E+00;
    COFD[2394] = -4.85617912E-01;
    COFD[2395] = 2.00461138E-02;
    COFD[2396] = -2.19594078E+01;
    COFD[2397] = 5.59157589E+00;
    COFD[2398] = -4.85617912E-01;
    COFD[2399] = 2.00461138E-02;
    COFD[2400] = -2.19670848E+01;
    COFD[2401] = 5.48847873E+00;
    COFD[2402] = -4.59558930E-01;
    COFD[2403] = 1.84107961E-02;
    COFD[2404] = -2.21070030E+01;
    COFD[2405] = 5.55072945E+00;
    COFD[2406] = -4.72525345E-01;
    COFD[2407] = 1.91674202E-02;
    COFD[2408] = -2.21157340E+01;
    COFD[2409] = 5.55072945E+00;
    COFD[2410] = -4.72525345E-01;
    COFD[2411] = 1.91674202E-02;
    COFD[2412] = -2.09241647E+01;
    COFD[2413] = 5.42316225E+00;
    COFD[2414] = -4.73702801E-01;
    COFD[2415] = 1.99217718E-02;
    COFD[2416] = -2.09222454E+01;
    COFD[2417] = 4.82184721E+00;
    COFD[2418] = -3.48128875E-01;
    COFD[2419] = 1.25918978E-02;
    COFD[2420] = -2.09222454E+01;
    COFD[2421] = 4.82184721E+00;
    COFD[2422] = -3.48128875E-01;
    COFD[2423] = 1.25918978E-02;
    COFD[2424] = -2.02268902E+01;
    COFD[2425] = 5.13632093E+00;
    COFD[2426] = -4.44839124E-01;
    COFD[2427] = 1.90058354E-02;
    COFD[2428] = -2.10026861E+01;
    COFD[2429] = 5.38326647E+00;
    COFD[2430] = -4.71201048E-01;
    COFD[2431] = 1.99207516E-02;
    COFD[2432] = -1.58456300E+01;
    COFD[2433] = 4.02074783E+00;
    COFD[2434] = -3.10018522E-01;
    COFD[2435] = 1.35599552E-02;
    COFD[2436] = -1.92718582E+01;
    COFD[2437] = 5.41172124E+00;
    COFD[2438] = -4.73213887E-01;
    COFD[2439] = 1.99405473E-02;
    COFD[2440] = -1.88179418E+01;
    COFD[2441] = 4.79683898E+00;
    COFD[2442] = -4.04829719E-01;
    COFD[2443] = 1.74325475E-02;
    COFD[2444] = -2.04928958E+01;
    COFD[2445] = 5.22397933E+00;
    COFD[2446] = -4.54138171E-01;
    COFD[2447] = 1.93249285E-02;
    COFD[2448] = -1.88378874E+01;
    COFD[2449] = 4.79683898E+00;
    COFD[2450] = -4.04829719E-01;
    COFD[2451] = 1.74325475E-02;
    COFD[2452] = -1.65295288E+01;
    COFD[2453] = 2.97569206E+00;
    COFD[2454] = -6.75652842E-02;
    COFD[2455] = -1.08648422E-03;
    COFD[2456] = -2.02637994E+01;
    COFD[2457] = 5.14984081E+00;
    COFD[2458] = -4.46093018E-01;
    COFD[2459] = 1.90396647E-02;
    COFD[2460] = -2.02710316E+01;
    COFD[2461] = 5.14984081E+00;
    COFD[2462] = -4.46093018E-01;
    COFD[2463] = 1.90396647E-02;
    COFD[2464] = -1.85157414E+01;
    COFD[2465] = 4.67076124E+00;
    COFD[2466] = -3.90022427E-01;
    COFD[2467] = 1.68533953E-02;
    COFD[2468] = -1.87476063E+01;
    COFD[2469] = 4.79683898E+00;
    COFD[2470] = -4.04829719E-01;
    COFD[2471] = 1.74325475E-02;
    COFD[2472] = -2.09376196E+01;
    COFD[2473] = 5.40870099E+00;
    COFD[2474] = -4.73017610E-01;
    COFD[2475] = 1.99399066E-02;
    COFD[2476] = -2.09376196E+01;
    COFD[2477] = 5.40870099E+00;
    COFD[2478] = -4.73017610E-01;
    COFD[2479] = 1.99399066E-02;
    COFD[2480] = -2.09612557E+01;
    COFD[2481] = 5.40870099E+00;
    COFD[2482] = -4.73017610E-01;
    COFD[2483] = 1.99399066E-02;
    COFD[2484] = -2.11381508E+01;
    COFD[2485] = 5.45574440E+00;
    COFD[2486] = -4.77436155E-01;
    COFD[2487] = 2.00644596E-02;
    COFD[2488] = -2.02642227E+01;
    COFD[2489] = 5.14499740E+00;
    COFD[2490] = -4.45694430E-01;
    COFD[2491] = 1.90318646E-02;
    COFD[2492] = -2.20421041E+01;
    COFD[2493] = 5.52708332E+00;
    COFD[2494] = -4.68000808E-01;
    COFD[2495] = 1.89131908E-02;
    COFD[2496] = -2.01801667E+01;
    COFD[2497] = 4.53183330E+00;
    COFD[2498] = -3.02186760E-01;
    COFD[2499] = 1.02756490E-02;
    COFD[2500] = -2.01889168E+01;
    COFD[2501] = 4.53183330E+00;
    COFD[2502] = -3.02186760E-01;
    COFD[2503] = 1.02756490E-02;
    COFD[2504] = -1.95877017E+01;
    COFD[2505] = 4.27643051E+00;
    COFD[2506] = -2.68040901E-01;
    COFD[2507] = 8.77650113E-03;
    COFD[2508] = -1.95877017E+01;
    COFD[2509] = 4.27643051E+00;
    COFD[2510] = -2.68040901E-01;
    COFD[2511] = 8.77650113E-03;
    COFD[2512] = -2.03599050E+01;
    COFD[2513] = 4.60682543E+00;
    COFD[2514] = -3.13971634E-01;
    COFD[2515] = 1.08661011E-02;
    COFD[2516] = -2.19383403E+01;
    COFD[2517] = 5.59157589E+00;
    COFD[2518] = -4.85617912E-01;
    COFD[2519] = 2.00461138E-02;
    COFD[2520] = -2.19491710E+01;
    COFD[2521] = 5.59157589E+00;
    COFD[2522] = -4.85617912E-01;
    COFD[2523] = 2.00461138E-02;
    COFD[2524] = -2.19594078E+01;
    COFD[2525] = 5.59157589E+00;
    COFD[2526] = -4.85617912E-01;
    COFD[2527] = 2.00461138E-02;
    COFD[2528] = -2.19670848E+01;
    COFD[2529] = 5.48847873E+00;
    COFD[2530] = -4.59558930E-01;
    COFD[2531] = 1.84107961E-02;
    COFD[2532] = -2.21070030E+01;
    COFD[2533] = 5.55072945E+00;
    COFD[2534] = -4.72525345E-01;
    COFD[2535] = 1.91674202E-02;
    COFD[2536] = -2.21157340E+01;
    COFD[2537] = 5.55072945E+00;
    COFD[2538] = -4.72525345E-01;
    COFD[2539] = 1.91674202E-02;
    COFD[2540] = -2.09241647E+01;
    COFD[2541] = 5.42316225E+00;
    COFD[2542] = -4.73702801E-01;
    COFD[2543] = 1.99217718E-02;
    COFD[2544] = -2.09222454E+01;
    COFD[2545] = 4.82184721E+00;
    COFD[2546] = -3.48128875E-01;
    COFD[2547] = 1.25918978E-02;
    COFD[2548] = -2.09222454E+01;
    COFD[2549] = 4.82184721E+00;
    COFD[2550] = -3.48128875E-01;
    COFD[2551] = 1.25918978E-02;
    COFD[2552] = -2.02268902E+01;
    COFD[2553] = 5.13632093E+00;
    COFD[2554] = -4.44839124E-01;
    COFD[2555] = 1.90058354E-02;
    COFD[2556] = -2.10026861E+01;
    COFD[2557] = 5.38326647E+00;
    COFD[2558] = -4.71201048E-01;
    COFD[2559] = 1.99207516E-02;
    COFD[2560] = -1.59537247E+01;
    COFD[2561] = 4.07051484E+00;
    COFD[2562] = -3.16303109E-01;
    COFD[2563] = 1.38259377E-02;
    COFD[2564] = -1.96866103E+01;
    COFD[2565] = 5.54637286E+00;
    COFD[2566] = -4.87070324E-01;
    COFD[2567] = 2.03983467E-02;
    COFD[2568] = -1.93364585E+01;
    COFD[2569] = 4.98286777E+00;
    COFD[2570] = -4.26970814E-01;
    COFD[2571] = 1.83122917E-02;
    COFD[2572] = -2.07595845E+01;
    COFD[2573] = 5.32244593E+00;
    COFD[2574] = -4.65829403E-01;
    COFD[2575] = 1.97895274E-02;
    COFD[2576] = -1.93566243E+01;
    COFD[2577] = 4.98286777E+00;
    COFD[2578] = -4.26970814E-01;
    COFD[2579] = 1.83122917E-02;
    COFD[2580] = -1.80253664E+01;
    COFD[2581] = 3.69199168E+00;
    COFD[2582] = -1.74005516E-01;
    COFD[2583] = 3.97694372E-03;
    COFD[2584] = -2.07672833E+01;
    COFD[2585] = 5.32244593E+00;
    COFD[2586] = -4.65829403E-01;
    COFD[2587] = 1.97895274E-02;
    COFD[2588] = -2.07746356E+01;
    COFD[2589] = 5.32244593E+00;
    COFD[2590] = -4.65829403E-01;
    COFD[2591] = 1.97895274E-02;
    COFD[2592] = -1.89671752E+01;
    COFD[2593] = 4.83076737E+00;
    COFD[2594] = -4.08802573E-01;
    COFD[2595] = 1.75875241E-02;
    COFD[2596] = -1.92654138E+01;
    COFD[2597] = 4.98286777E+00;
    COFD[2598] = -4.26970814E-01;
    COFD[2599] = 1.83122917E-02;
    COFD[2600] = -2.13538553E+01;
    COFD[2601] = 5.54007827E+00;
    COFD[2602] = -4.86434511E-01;
    COFD[2603] = 2.03779006E-02;
    COFD[2604] = -2.13538553E+01;
    COFD[2605] = 5.54007827E+00;
    COFD[2606] = -4.86434511E-01;
    COFD[2607] = 2.03779006E-02;
    COFD[2608] = -2.13777308E+01;
    COFD[2609] = 5.54007827E+00;
    COFD[2610] = -4.86434511E-01;
    COFD[2611] = 2.03779006E-02;
    COFD[2612] = -2.13319784E+01;
    COFD[2613] = 5.52422470E+00;
    COFD[2614] = -4.84872944E-01;
    COFD[2615] = 2.03298213E-02;
    COFD[2616] = -2.04144604E+01;
    COFD[2617] = 5.19614628E+00;
    COFD[2618] = -4.50889164E-01;
    COFD[2619] = 1.91983328E-02;
    COFD[2620] = -2.20063594E+01;
    COFD[2621] = 5.48540187E+00;
    COFD[2622] = -4.58962148E-01;
    COFD[2623] = 1.83770355E-02;
    COFD[2624] = -1.93125662E+01;
    COFD[2625] = 4.10954793E+00;
    COFD[2626] = -2.37523329E-01;
    COFD[2627] = 7.08858141E-03;
    COFD[2628] = -1.93214527E+01;
    COFD[2629] = 4.10954793E+00;
    COFD[2630] = -2.37523329E-01;
    COFD[2631] = 7.08858141E-03;
    COFD[2632] = -2.03599050E+01;
    COFD[2633] = 4.60682543E+00;
    COFD[2634] = -3.13971634E-01;
    COFD[2635] = 1.08661011E-02;
    COFD[2636] = -2.03599050E+01;
    COFD[2637] = 4.60682543E+00;
    COFD[2638] = -3.13971634E-01;
    COFD[2639] = 1.08661011E-02;
    COFD[2640] = -1.95785144E+01;
    COFD[2641] = 4.22062499E+00;
    COFD[2642] = -2.54326872E-01;
    COFD[2643] = 7.91017784E-03;
    COFD[2644] = -2.20378456E+01;
    COFD[2645] = 5.58129885E+00;
    COFD[2646] = -4.78532921E-01;
    COFD[2647] = 1.95095699E-02;
    COFD[2648] = -2.20488323E+01;
    COFD[2649] = 5.58129885E+00;
    COFD[2650] = -4.78532921E-01;
    COFD[2651] = 1.95095699E-02;
    COFD[2652] = -2.20592197E+01;
    COFD[2653] = 5.58129885E+00;
    COFD[2654] = -4.78532921E-01;
    COFD[2655] = 1.95095699E-02;
    COFD[2656] = -2.17419207E+01;
    COFD[2657] = 5.33732875E+00;
    COFD[2658] = -4.32453425E-01;
    COFD[2659] = 1.69373833E-02;
    COFD[2660] = -2.20028184E+01;
    COFD[2661] = 5.45178028E+00;
    COFD[2662] = -4.52847771E-01;
    COFD[2663] = 1.80418544E-02;
    COFD[2664] = -2.20116855E+01;
    COFD[2665] = 5.45178028E+00;
    COFD[2666] = -4.52847771E-01;
    COFD[2667] = 1.80418544E-02;
    COFD[2668] = -2.13796303E+01;
    COFD[2669] = 5.56978987E+00;
    COFD[2670] = -4.89141980E-01;
    COFD[2671] = 2.04499210E-02;
    COFD[2672] = -2.03036402E+01;
    COFD[2673] = 4.50250781E+00;
    COFD[2674] = -2.97622106E-01;
    COFD[2675] = 1.00481473E-02;
    COFD[2676] = -2.03036402E+01;
    COFD[2677] = 4.50250781E+00;
    COFD[2678] = -2.97622106E-01;
    COFD[2679] = 1.00481473E-02;
    COFD[2680] = -2.03844252E+01;
    COFD[2681] = 5.18856872E+00;
    COFD[2682] = -4.50001829E-01;
    COFD[2683] = 1.91636142E-02;
    COFD[2684] = -2.13502344E+01;
    COFD[2685] = 5.48617067E+00;
    COFD[2686] = -4.80816776E-01;
    COFD[2687] = 2.01887868E-02;
    COFD[2688] = -1.34695359E+01;
    COFD[2689] = 3.09379603E+00;
    COFD[2690] = -1.91268635E-01;
    COFD[2691] = 8.47480224E-03;
    COFD[2692] = -1.72224724E+01;
    COFD[2693] = 4.69060745E+00;
    COFD[2694] = -3.92369888E-01;
    COFD[2695] = 1.69459661E-02;
    COFD[2696] = -1.65412306E+01;
    COFD[2697] = 3.95035840E+00;
    COFD[2698] = -3.00959418E-01;
    COFD[2699] = 1.31692593E-02;
    COFD[2700] = -1.78725135E+01;
    COFD[2701] = 4.29613154E+00;
    COFD[2702] = -3.44012526E-01;
    COFD[2703] = 1.49643715E-02;
    COFD[2704] = -1.65596434E+01;
    COFD[2705] = 3.95035840E+00;
    COFD[2706] = -3.00959418E-01;
    COFD[2707] = 1.31692593E-02;
    COFD[2708] = -2.15014310E+01;
    COFD[2709] = 5.46737673E+00;
    COFD[2710] = -4.55696085E-01;
    COFD[2711] = 1.81982625E-02;
    COFD[2712] = -1.78792605E+01;
    COFD[2713] = 4.29613154E+00;
    COFD[2714] = -3.44012526E-01;
    COFD[2715] = 1.49643715E-02;
    COFD[2716] = -1.78856918E+01;
    COFD[2717] = 4.29613154E+00;
    COFD[2718] = -3.44012526E-01;
    COFD[2719] = 1.49643715E-02;
    COFD[2720] = -1.61038806E+01;
    COFD[2721] = 3.75910622E+00;
    COFD[2722] = -2.75986578E-01;
    COFD[2723] = 1.20782843E-02;
    COFD[2724] = -1.64758697E+01;
    COFD[2725] = 3.95035840E+00;
    COFD[2726] = -3.00959418E-01;
    COFD[2727] = 1.31692593E-02;
    COFD[2728] = -1.88171077E+01;
    COFD[2729] = 4.68393046E+00;
    COFD[2730] = -3.91610863E-01;
    COFD[2731] = 1.69174645E-02;
    COFD[2732] = -1.88171077E+01;
    COFD[2733] = 4.68393046E+00;
    COFD[2734] = -3.91610863E-01;
    COFD[2735] = 1.69174645E-02;
    COFD[2736] = -1.88390649E+01;
    COFD[2737] = 4.68393046E+00;
    COFD[2738] = -3.91610863E-01;
    COFD[2739] = 1.69174645E-02;
    COFD[2740] = -1.87821119E+01;
    COFD[2741] = 4.66162351E+00;
    COFD[2742] = -3.88920477E-01;
    COFD[2743] = 1.68089648E-02;
    COFD[2744] = -1.76182365E+01;
    COFD[2745] = 4.19935698E+00;
    COFD[2746] = -3.32310212E-01;
    COFD[2747] = 1.44920670E-02;
    COFD[2748] = -2.09066354E+01;
    COFD[2749] = 5.30153901E+00;
    COFD[2750] = -4.63335119E-01;
    COFD[2751] = 1.96897053E-02;
    COFD[2752] = -2.19851338E+01;
    COFD[2753] = 5.55935694E+00;
    COFD[2754] = -4.74154740E-01;
    COFD[2755] = 1.92584304E-02;
    COFD[2756] = -2.19929679E+01;
    COFD[2757] = 5.55935694E+00;
    COFD[2758] = -4.74154740E-01;
    COFD[2759] = 1.92584304E-02;
    COFD[2760] = -2.19383403E+01;
    COFD[2761] = 5.59157589E+00;
    COFD[2762] = -4.85617912E-01;
    COFD[2763] = 2.00461138E-02;
    COFD[2764] = -2.19383403E+01;
    COFD[2765] = 5.59157589E+00;
    COFD[2766] = -4.85617912E-01;
    COFD[2767] = 2.00461138E-02;
    COFD[2768] = -2.20378456E+01;
    COFD[2769] = 5.58129885E+00;
    COFD[2770] = -4.78532921E-01;
    COFD[2771] = 1.95095699E-02;
    COFD[2772] = -2.03569548E+01;
    COFD[2773] = 5.13263469E+00;
    COFD[2774] = -4.44457285E-01;
    COFD[2775] = 1.89932102E-02;
    COFD[2776] = -2.03667275E+01;
    COFD[2777] = 5.13263469E+00;
    COFD[2778] = -4.44457285E-01;
    COFD[2779] = 1.89932102E-02;
    COFD[2780] = -2.03759451E+01;
    COFD[2781] = 5.13263469E+00;
    COFD[2782] = -4.44457285E-01;
    COFD[2783] = 1.89932102E-02;
    COFD[2784] = -2.12018018E+01;
    COFD[2785] = 5.39823225E+00;
    COFD[2786] = -4.72294645E-01;
    COFD[2787] = 1.99340225E-02;
    COFD[2788] = -2.10759163E+01;
    COFD[2789] = 5.34286099E+00;
    COFD[2790] = -4.68100992E-01;
    COFD[2791] = 1.98731399E-02;
    COFD[2792] = -2.10837327E+01;
    COFD[2793] = 5.34286099E+00;
    COFD[2794] = -4.68100992E-01;
    COFD[2795] = 1.98731399E-02;
    COFD[2796] = -1.88524385E+01;
    COFD[2797] = 4.72476764E+00;
    COFD[2798] = -3.96306836E-01;
    COFD[2799] = 1.70964541E-02;
    COFD[2800] = -2.20942491E+01;
    COFD[2801] = 5.58360799E+00;
    COFD[2802] = -4.82701436E-01;
    COFD[2803] = 1.98437922E-02;
    COFD[2804] = -2.20942491E+01;
    COFD[2805] = 5.58360799E+00;
    COFD[2806] = -4.82701436E-01;
    COFD[2807] = 1.98437922E-02;
    COFD[2808] = -1.75898751E+01;
    COFD[2809] = 4.19171952E+00;
    COFD[2810] = -3.31354810E-01;
    COFD[2811] = 1.44520623E-02;
    COFD[2812] = -1.87596895E+01;
    COFD[2813] = 4.61078776E+00;
    COFD[2814] = -3.82625667E-01;
    COFD[2815] = 1.65478601E-02;
    COFD[2816] = -1.34709807E+01;
    COFD[2817] = 3.09379603E+00;
    COFD[2818] = -1.91268635E-01;
    COFD[2819] = 8.47480224E-03;
    COFD[2820] = -1.72232223E+01;
    COFD[2821] = 4.69060745E+00;
    COFD[2822] = -3.92369888E-01;
    COFD[2823] = 1.69459661E-02;
    COFD[2824] = -1.65488358E+01;
    COFD[2825] = 3.95035840E+00;
    COFD[2826] = -3.00959418E-01;
    COFD[2827] = 1.31692593E-02;
    COFD[2828] = -1.78834935E+01;
    COFD[2829] = 4.29613154E+00;
    COFD[2830] = -3.44012526E-01;
    COFD[2831] = 1.49643715E-02;
    COFD[2832] = -1.65675362E+01;
    COFD[2833] = 3.95035840E+00;
    COFD[2834] = -3.00959418E-01;
    COFD[2835] = 1.31692593E-02;
    COFD[2836] = -2.15095980E+01;
    COFD[2837] = 5.46737673E+00;
    COFD[2838] = -4.55696085E-01;
    COFD[2839] = 1.81982625E-02;
    COFD[2840] = -1.78903913E+01;
    COFD[2841] = 4.29613154E+00;
    COFD[2842] = -3.44012526E-01;
    COFD[2843] = 1.49643715E-02;
    COFD[2844] = -1.78969684E+01;
    COFD[2845] = 4.29613154E+00;
    COFD[2846] = -3.44012526E-01;
    COFD[2847] = 1.49643715E-02;
    COFD[2848] = -1.61101966E+01;
    COFD[2849] = 3.75910622E+00;
    COFD[2850] = -2.75986578E-01;
    COFD[2851] = 1.20782843E-02;
    COFD[2852] = -1.64825368E+01;
    COFD[2853] = 3.95035840E+00;
    COFD[2854] = -3.00959418E-01;
    COFD[2855] = 1.31692593E-02;
    COFD[2856] = -1.88241079E+01;
    COFD[2857] = 4.68393046E+00;
    COFD[2858] = -3.91610863E-01;
    COFD[2859] = 1.69174645E-02;
    COFD[2860] = -1.88241079E+01;
    COFD[2861] = 4.68393046E+00;
    COFD[2862] = -3.91610863E-01;
    COFD[2863] = 1.69174645E-02;
    COFD[2864] = -1.88463816E+01;
    COFD[2865] = 4.68393046E+00;
    COFD[2866] = -3.91610863E-01;
    COFD[2867] = 1.69174645E-02;
    COFD[2868] = -1.87897298E+01;
    COFD[2869] = 4.66162351E+00;
    COFD[2870] = -3.88920477E-01;
    COFD[2871] = 1.68089648E-02;
    COFD[2872] = -1.76285640E+01;
    COFD[2873] = 4.19935698E+00;
    COFD[2874] = -3.32310212E-01;
    COFD[2875] = 1.44920670E-02;
    COFD[2876] = -2.09191285E+01;
    COFD[2877] = 5.30153901E+00;
    COFD[2878] = -4.63335119E-01;
    COFD[2879] = 1.96897053E-02;
    COFD[2880] = -2.19956352E+01;
    COFD[2881] = 5.55935694E+00;
    COFD[2882] = -4.74154740E-01;
    COFD[2883] = 1.92584304E-02;
    COFD[2884] = -2.20036369E+01;
    COFD[2885] = 5.55935694E+00;
    COFD[2886] = -4.74154740E-01;
    COFD[2887] = 1.92584304E-02;
    COFD[2888] = -2.19491710E+01;
    COFD[2889] = 5.59157589E+00;
    COFD[2890] = -4.85617912E-01;
    COFD[2891] = 2.00461138E-02;
    COFD[2892] = -2.19491710E+01;
    COFD[2893] = 5.59157589E+00;
    COFD[2894] = -4.85617912E-01;
    COFD[2895] = 2.00461138E-02;
    COFD[2896] = -2.20488323E+01;
    COFD[2897] = 5.58129885E+00;
    COFD[2898] = -4.78532921E-01;
    COFD[2899] = 1.95095699E-02;
    COFD[2900] = -2.03667275E+01;
    COFD[2901] = 5.13263469E+00;
    COFD[2902] = -4.44457285E-01;
    COFD[2903] = 1.89932102E-02;
    COFD[2904] = -2.03766950E+01;
    COFD[2905] = 5.13263469E+00;
    COFD[2906] = -4.44457285E-01;
    COFD[2907] = 1.89932102E-02;
    COFD[2908] = -2.03861000E+01;
    COFD[2909] = 5.13263469E+00;
    COFD[2910] = -4.44457285E-01;
    COFD[2911] = 1.89932102E-02;
    COFD[2912] = -2.12121370E+01;
    COFD[2913] = 5.39823225E+00;
    COFD[2914] = -4.72294645E-01;
    COFD[2915] = 1.99340225E-02;
    COFD[2916] = -2.10864251E+01;
    COFD[2917] = 5.34286099E+00;
    COFD[2918] = -4.68100992E-01;
    COFD[2919] = 1.98731399E-02;
    COFD[2920] = -2.10944088E+01;
    COFD[2921] = 5.34286099E+00;
    COFD[2922] = -4.68100992E-01;
    COFD[2923] = 1.98731399E-02;
    COFD[2924] = -1.88646070E+01;
    COFD[2925] = 4.72476764E+00;
    COFD[2926] = -3.96306836E-01;
    COFD[2927] = 1.70964541E-02;
    COFD[2928] = -2.21065306E+01;
    COFD[2929] = 5.58360799E+00;
    COFD[2930] = -4.82701436E-01;
    COFD[2931] = 1.98437922E-02;
    COFD[2932] = -2.21065306E+01;
    COFD[2933] = 5.58360799E+00;
    COFD[2934] = -4.82701436E-01;
    COFD[2935] = 1.98437922E-02;
    COFD[2936] = -1.76002031E+01;
    COFD[2937] = 4.19171952E+00;
    COFD[2938] = -3.31354810E-01;
    COFD[2939] = 1.44520623E-02;
    COFD[2940] = -1.87717330E+01;
    COFD[2941] = 4.61078776E+00;
    COFD[2942] = -3.82625667E-01;
    COFD[2943] = 1.65478601E-02;
    COFD[2944] = -1.34723215E+01;
    COFD[2945] = 3.09379603E+00;
    COFD[2946] = -1.91268635E-01;
    COFD[2947] = 8.47480224E-03;
    COFD[2948] = -1.72239172E+01;
    COFD[2949] = 4.69060745E+00;
    COFD[2950] = -3.92369888E-01;
    COFD[2951] = 1.69459661E-02;
    COFD[2952] = -1.65559787E+01;
    COFD[2953] = 3.95035840E+00;
    COFD[2954] = -3.00959418E-01;
    COFD[2955] = 1.31692593E-02;
    COFD[2956] = -1.78938745E+01;
    COFD[2957] = 4.29613154E+00;
    COFD[2958] = -3.44012526E-01;
    COFD[2959] = 1.49643715E-02;
    COFD[2960] = -1.65749533E+01;
    COFD[2961] = 3.95035840E+00;
    COFD[2962] = -3.00959418E-01;
    COFD[2963] = 1.31692593E-02;
    COFD[2964] = -2.15172770E+01;
    COFD[2965] = 5.46737673E+00;
    COFD[2966] = -4.55696085E-01;
    COFD[2967] = 1.81982625E-02;
    COFD[2968] = -1.79009181E+01;
    COFD[2969] = 4.29613154E+00;
    COFD[2970] = -3.44012526E-01;
    COFD[2971] = 1.49643715E-02;
    COFD[2972] = -1.79076361E+01;
    COFD[2973] = 4.29613154E+00;
    COFD[2974] = -3.44012526E-01;
    COFD[2975] = 1.49643715E-02;
    COFD[2976] = -1.61161138E+01;
    COFD[2977] = 3.75910622E+00;
    COFD[2978] = -2.75986578E-01;
    COFD[2979] = 1.20782843E-02;
    COFD[2980] = -1.64887871E+01;
    COFD[2981] = 3.95035840E+00;
    COFD[2982] = -3.00959418E-01;
    COFD[2983] = 1.31692593E-02;
    COFD[2984] = -1.88306747E+01;
    COFD[2985] = 4.68393046E+00;
    COFD[2986] = -3.91610863E-01;
    COFD[2987] = 1.69174645E-02;
    COFD[2988] = -1.88306747E+01;
    COFD[2989] = 4.68393046E+00;
    COFD[2990] = -3.91610863E-01;
    COFD[2991] = 1.69174645E-02;
    COFD[2992] = -1.88532497E+01;
    COFD[2993] = 4.68393046E+00;
    COFD[2994] = -3.91610863E-01;
    COFD[2995] = 1.69174645E-02;
    COFD[2996] = -1.87968848E+01;
    COFD[2997] = 4.66162351E+00;
    COFD[2998] = -3.88920477E-01;
    COFD[2999] = 1.68089648E-02;
    COFD[3000] = -1.76383156E+01;
    COFD[3001] = 4.19935698E+00;
    COFD[3002] = -3.32310212E-01;
    COFD[3003] = 1.44920670E-02;
    COFD[3004] = -2.09309753E+01;
    COFD[3005] = 5.30153901E+00;
    COFD[3006] = -4.63335119E-01;
    COFD[3007] = 1.96897053E-02;
    COFD[3008] = -2.20055544E+01;
    COFD[3009] = 5.55935694E+00;
    COFD[3010] = -4.74154740E-01;
    COFD[3011] = 1.92584304E-02;
    COFD[3012] = -2.20137178E+01;
    COFD[3013] = 5.55935694E+00;
    COFD[3014] = -4.74154740E-01;
    COFD[3015] = 1.92584304E-02;
    COFD[3016] = -2.19594078E+01;
    COFD[3017] = 5.59157589E+00;
    COFD[3018] = -4.85617912E-01;
    COFD[3019] = 2.00461138E-02;
    COFD[3020] = -2.19594078E+01;
    COFD[3021] = 5.59157589E+00;
    COFD[3022] = -4.85617912E-01;
    COFD[3023] = 2.00461138E-02;
    COFD[3024] = -2.20592197E+01;
    COFD[3025] = 5.58129885E+00;
    COFD[3026] = -4.78532921E-01;
    COFD[3027] = 1.95095699E-02;
    COFD[3028] = -2.03759451E+01;
    COFD[3029] = 5.13263469E+00;
    COFD[3030] = -4.44457285E-01;
    COFD[3031] = 1.89932102E-02;
    COFD[3032] = -2.03861000E+01;
    COFD[3033] = 5.13263469E+00;
    COFD[3034] = -4.44457285E-01;
    COFD[3035] = 1.89932102E-02;
    COFD[3036] = -2.03956853E+01;
    COFD[3037] = 5.13263469E+00;
    COFD[3038] = -4.44457285E-01;
    COFD[3039] = 1.89932102E-02;
    COFD[3040] = -2.12218960E+01;
    COFD[3041] = 5.39823225E+00;
    COFD[3042] = -4.72294645E-01;
    COFD[3043] = 1.99340225E-02;
    COFD[3044] = -2.10963515E+01;
    COFD[3045] = 5.34286099E+00;
    COFD[3046] = -4.68100992E-01;
    COFD[3047] = 1.98731399E-02;
    COFD[3048] = -2.11044965E+01;
    COFD[3049] = 5.34286099E+00;
    COFD[3050] = -4.68100992E-01;
    COFD[3051] = 1.98731399E-02;
    COFD[3052] = -1.88761387E+01;
    COFD[3053] = 4.72476764E+00;
    COFD[3054] = -3.96306836E-01;
    COFD[3055] = 1.70964541E-02;
    COFD[3056] = -2.21181720E+01;
    COFD[3057] = 5.58360799E+00;
    COFD[3058] = -4.82701436E-01;
    COFD[3059] = 1.98437922E-02;
    COFD[3060] = -2.21181720E+01;
    COFD[3061] = 5.58360799E+00;
    COFD[3062] = -4.82701436E-01;
    COFD[3063] = 1.98437922E-02;
    COFD[3064] = -1.76099552E+01;
    COFD[3065] = 4.19171952E+00;
    COFD[3066] = -3.31354810E-01;
    COFD[3067] = 1.44520623E-02;
    COFD[3068] = -1.87831434E+01;
    COFD[3069] = 4.61078776E+00;
    COFD[3070] = -3.82625667E-01;
    COFD[3071] = 1.65478601E-02;
    COFD[3072] = -1.42229194E+01;
    COFD[3073] = 3.38669384E+00;
    COFD[3074] = -2.28784122E-01;
    COFD[3075] = 1.00790953E-02;
    COFD[3076] = -1.82251914E+01;
    COFD[3077] = 5.05237312E+00;
    COFD[3078] = -4.35182396E-01;
    COFD[3079] = 1.86363074E-02;
    COFD[3080] = -1.74792112E+01;
    COFD[3081] = 4.29676909E+00;
    COFD[3082] = -3.44085306E-01;
    COFD[3083] = 1.49671135E-02;
    COFD[3084] = -1.89544778E+01;
    COFD[3085] = 4.68595732E+00;
    COFD[3086] = -3.91842840E-01;
    COFD[3087] = 1.69262542E-02;
    COFD[3088] = -1.74984476E+01;
    COFD[3089] = 4.29676909E+00;
    COFD[3090] = -3.44085306E-01;
    COFD[3091] = 1.49671135E-02;
    COFD[3092] = -2.08812333E+01;
    COFD[3093] = 5.08859217E+00;
    COFD[3094] = -3.90525428E-01;
    COFD[3095] = 1.47376395E-02;
    COFD[3096] = -1.89616623E+01;
    COFD[3097] = 4.68595732E+00;
    COFD[3098] = -3.91842840E-01;
    COFD[3099] = 1.69262542E-02;
    COFD[3100] = -1.89685165E+01;
    COFD[3101] = 4.68595732E+00;
    COFD[3102] = -3.91842840E-01;
    COFD[3103] = 1.69262542E-02;
    COFD[3104] = -1.71884218E+01;
    COFD[3105] = 4.17190426E+00;
    COFD[3106] = -3.28894681E-01;
    COFD[3107] = 1.43498101E-02;
    COFD[3108] = -1.74111692E+01;
    COFD[3109] = 4.29676909E+00;
    COFD[3110] = -3.44085306E-01;
    COFD[3111] = 1.49671135E-02;
    COFD[3112] = -1.98418115E+01;
    COFD[3113] = 5.04367502E+00;
    COFD[3114] = -4.34153325E-01;
    COFD[3115] = 1.85956055E-02;
    COFD[3116] = -1.98418115E+01;
    COFD[3117] = 5.04367502E+00;
    COFD[3118] = -4.34153325E-01;
    COFD[3119] = 1.85956055E-02;
    COFD[3120] = -1.98646734E+01;
    COFD[3121] = 5.04367502E+00;
    COFD[3122] = -4.34153325E-01;
    COFD[3123] = 1.85956055E-02;
    COFD[3124] = -1.98075055E+01;
    COFD[3125] = 5.02169524E+00;
    COFD[3126] = -4.31582804E-01;
    COFD[3127] = 1.84953568E-02;
    COFD[3128] = -1.86157761E+01;
    COFD[3129] = 4.55689508E+00;
    COFD[3130] = -3.75937921E-01;
    COFD[3131] = 1.62703488E-02;
    COFD[3132] = -2.16802612E+01;
    COFD[3133] = 5.52918296E+00;
    COFD[3134] = -4.85360709E-01;
    COFD[3135] = 2.03448006E-02;
    COFD[3136] = -2.16296373E+01;
    COFD[3137] = 5.29019717E+00;
    COFD[3138] = -4.24502606E-01;
    COFD[3139] = 1.65197343E-02;
    COFD[3140] = -2.16379567E+01;
    COFD[3141] = 5.29019717E+00;
    COFD[3142] = -4.24502606E-01;
    COFD[3143] = 1.65197343E-02;
    COFD[3144] = -2.19670848E+01;
    COFD[3145] = 5.48847873E+00;
    COFD[3146] = -4.59558930E-01;
    COFD[3147] = 1.84107961E-02;
    COFD[3148] = -2.19670848E+01;
    COFD[3149] = 5.48847873E+00;
    COFD[3150] = -4.59558930E-01;
    COFD[3151] = 1.84107961E-02;
    COFD[3152] = -2.17419207E+01;
    COFD[3153] = 5.33732875E+00;
    COFD[3154] = -4.32453425E-01;
    COFD[3155] = 1.69373833E-02;
    COFD[3156] = -2.12018018E+01;
    COFD[3157] = 5.39823225E+00;
    COFD[3158] = -4.72294645E-01;
    COFD[3159] = 1.99340225E-02;
    COFD[3160] = -2.12121370E+01;
    COFD[3161] = 5.39823225E+00;
    COFD[3162] = -4.72294645E-01;
    COFD[3163] = 1.99340225E-02;
    COFD[3164] = -2.12218960E+01;
    COFD[3165] = 5.39823225E+00;
    COFD[3166] = -4.72294645E-01;
    COFD[3167] = 1.99340225E-02;
    COFD[3168] = -2.19327397E+01;
    COFD[3169] = 5.60638188E+00;
    COFD[3170] = -4.91272522E-01;
    COFD[3171] = 2.04396264E-02;
    COFD[3172] = -2.18190539E+01;
    COFD[3173] = 5.55753905E+00;
    COFD[3174] = -4.88136714E-01;
    COFD[3175] = 2.04294957E-02;
    COFD[3176] = -2.18273547E+01;
    COFD[3177] = 5.55753905E+00;
    COFD[3178] = -4.88136714E-01;
    COFD[3179] = 2.04294957E-02;
    COFD[3180] = -1.99081556E+01;
    COFD[3181] = 5.09311649E+00;
    COFD[3182] = -4.39965178E-01;
    COFD[3183] = 1.88238537E-02;
    COFD[3184] = -2.20453723E+01;
    COFD[3185] = 5.44448440E+00;
    COFD[3186] = -4.51529024E-01;
    COFD[3187] = 1.79698119E-02;
    COFD[3188] = -2.20453723E+01;
    COFD[3189] = 5.44448440E+00;
    COFD[3190] = -4.51529024E-01;
    COFD[3191] = 1.79698119E-02;
    COFD[3192] = -1.85864144E+01;
    COFD[3193] = 4.54915847E+00;
    COFD[3194] = -3.75000738E-01;
    COFD[3195] = 1.62324821E-02;
    COFD[3196] = -1.98040322E+01;
    COFD[3197] = 4.97569695E+00;
    COFD[3198] = -4.26123307E-01;
    COFD[3199] = 1.82788664E-02;
    COFD[3200] = -1.39913897E+01;
    COFD[3201] = 3.26384506E+00;
    COFD[3202] = -2.12947087E-01;
    COFD[3203] = 9.39743888E-03;
    COFD[3204] = -1.79339327E+01;
    COFD[3205] = 4.91373893E+00;
    COFD[3206] = -4.18747629E-01;
    COFD[3207] = 1.79856610E-02;
    COFD[3208] = -1.72496634E+01;
    COFD[3209] = 4.17889917E+00;
    COFD[3210] = -3.29752510E-01;
    COFD[3211] = 1.43850275E-02;
    COFD[3212] = -1.86335932E+01;
    COFD[3213] = 4.53572533E+00;
    COFD[3214] = -3.73386925E-01;
    COFD[3215] = 1.61678881E-02;
    COFD[3216] = -1.72691500E+01;
    COFD[3217] = 4.17889917E+00;
    COFD[3218] = -3.29752510E-01;
    COFD[3219] = 1.43850275E-02;
    COFD[3220] = -2.12597312E+01;
    COFD[3221] = 5.24930667E+00;
    COFD[3222] = -4.17435088E-01;
    COFD[3223] = 1.61434424E-02;
    COFD[3224] = -1.86409139E+01;
    COFD[3225] = 4.53572533E+00;
    COFD[3226] = -3.73386925E-01;
    COFD[3227] = 1.61678881E-02;
    COFD[3228] = -1.86479000E+01;
    COFD[3229] = 4.53572533E+00;
    COFD[3230] = -3.73386925E-01;
    COFD[3231] = 1.61678881E-02;
    COFD[3232] = -1.69491115E+01;
    COFD[3233] = 4.05099737E+00;
    COFD[3234] = -3.13841660E-01;
    COFD[3235] = 1.37218854E-02;
    COFD[3236] = -1.71808106E+01;
    COFD[3237] = 4.17889917E+00;
    COFD[3238] = -3.29752510E-01;
    COFD[3239] = 1.43850275E-02;
    COFD[3240] = -1.95263312E+01;
    COFD[3241] = 4.90255048E+00;
    COFD[3242] = -4.17368501E-01;
    COFD[3243] = 1.79287358E-02;
    COFD[3244] = -1.95263312E+01;
    COFD[3245] = 4.90255048E+00;
    COFD[3246] = -4.17368501E-01;
    COFD[3247] = 1.79287358E-02;
    COFD[3248] = -1.95494668E+01;
    COFD[3249] = 4.90255048E+00;
    COFD[3250] = -4.17368501E-01;
    COFD[3251] = 1.79287358E-02;
    COFD[3252] = -1.94763688E+01;
    COFD[3253] = 4.87333294E+00;
    COFD[3254] = -4.13769241E-01;
    COFD[3255] = 1.77802244E-02;
    COFD[3256] = -1.83455435E+01;
    COFD[3257] = 4.42828044E+00;
    COFD[3258] = -3.60417833E-01;
    COFD[3259] = 1.56455103E-02;
    COFD[3260] = -2.14224484E+01;
    COFD[3261] = 5.41729961E+00;
    COFD[3262] = -4.73400281E-01;
    COFD[3263] = 1.99269567E-02;
    COFD[3264] = -2.19229190E+01;
    COFD[3265] = 5.41841631E+00;
    COFD[3266] = -4.46818971E-01;
    COFD[3267] = 1.77127652E-02;
    COFD[3268] = -2.19313890E+01;
    COFD[3269] = 5.41841631E+00;
    COFD[3270] = -4.46818971E-01;
    COFD[3271] = 1.77127652E-02;
    COFD[3272] = -2.21070030E+01;
    COFD[3273] = 5.55072945E+00;
    COFD[3274] = -4.72525345E-01;
    COFD[3275] = 1.91674202E-02;
    COFD[3276] = -2.21070030E+01;
    COFD[3277] = 5.55072945E+00;
    COFD[3278] = -4.72525345E-01;
    COFD[3279] = 1.91674202E-02;
    COFD[3280] = -2.20028184E+01;
    COFD[3281] = 5.45178028E+00;
    COFD[3282] = -4.52847771E-01;
    COFD[3283] = 1.80418544E-02;
    COFD[3284] = -2.10759163E+01;
    COFD[3285] = 5.34286099E+00;
    COFD[3286] = -4.68100992E-01;
    COFD[3287] = 1.98731399E-02;
    COFD[3288] = -2.10864251E+01;
    COFD[3289] = 5.34286099E+00;
    COFD[3290] = -4.68100992E-01;
    COFD[3291] = 1.98731399E-02;
    COFD[3292] = -2.10963515E+01;
    COFD[3293] = 5.34286099E+00;
    COFD[3294] = -4.68100992E-01;
    COFD[3295] = 1.98731399E-02;
    COFD[3296] = -2.18190539E+01;
    COFD[3297] = 5.55753905E+00;
    COFD[3298] = -4.88136714E-01;
    COFD[3299] = 2.04294957E-02;
    COFD[3300] = -2.15575659E+01;
    COFD[3301] = 5.44803850E+00;
    COFD[3302] = -4.76610560E-01;
    COFD[3303] = 2.00355294E-02;
    COFD[3304] = -2.15660171E+01;
    COFD[3305] = 5.44803850E+00;
    COFD[3306] = -4.76610560E-01;
    COFD[3307] = 2.00355294E-02;
    COFD[3308] = -1.96309042E+01;
    COFD[3309] = 4.95923807E+00;
    COFD[3310] = -4.24176182E-01;
    COFD[3311] = 1.82020215E-02;
    COFD[3312] = -2.22059540E+01;
    COFD[3313] = 5.51722375E+00;
    COFD[3314] = -4.66081431E-01;
    COFD[3315] = 1.88044011E-02;
    COFD[3316] = -2.22059540E+01;
    COFD[3317] = 5.51722375E+00;
    COFD[3318] = -4.66081431E-01;
    COFD[3319] = 1.88044011E-02;
    COFD[3320] = -1.83166353E+01;
    COFD[3321] = 4.42045763E+00;
    COFD[3322] = -3.59451578E-01;
    COFD[3323] = 1.56056164E-02;
    COFD[3324] = -1.94928815E+01;
    COFD[3325] = 4.83189721E+00;
    COFD[3326] = -4.08932249E-01;
    COFD[3327] = 1.75924650E-02;
    COFD[3328] = -1.39924781E+01;
    COFD[3329] = 3.26384506E+00;
    COFD[3330] = -2.12947087E-01;
    COFD[3331] = 9.39743888E-03;
    COFD[3332] = -1.79344949E+01;
    COFD[3333] = 4.91373893E+00;
    COFD[3334] = -4.18747629E-01;
    COFD[3335] = 1.79856610E-02;
    COFD[3336] = -1.72556499E+01;
    COFD[3337] = 4.17889917E+00;
    COFD[3338] = -3.29752510E-01;
    COFD[3339] = 1.43850275E-02;
    COFD[3340] = -1.86424545E+01;
    COFD[3341] = 4.53572533E+00;
    COFD[3342] = -3.73386925E-01;
    COFD[3343] = 1.61678881E-02;
    COFD[3344] = -1.72753760E+01;
    COFD[3345] = 4.17889917E+00;
    COFD[3346] = -3.29752510E-01;
    COFD[3347] = 1.43850275E-02;
    COFD[3348] = -2.12661865E+01;
    COFD[3349] = 5.24930667E+00;
    COFD[3350] = -4.17435088E-01;
    COFD[3351] = 1.61434424E-02;
    COFD[3352] = -1.86499071E+01;
    COFD[3353] = 4.53572533E+00;
    COFD[3354] = -3.73386925E-01;
    COFD[3355] = 1.61678881E-02;
    COFD[3356] = -1.86570209E+01;
    COFD[3357] = 4.53572533E+00;
    COFD[3358] = -3.73386925E-01;
    COFD[3359] = 1.61678881E-02;
    COFD[3360] = -1.69540369E+01;
    COFD[3361] = 4.05099737E+00;
    COFD[3362] = -3.13841660E-01;
    COFD[3363] = 1.37218854E-02;
    COFD[3364] = -1.71860230E+01;
    COFD[3365] = 4.17889917E+00;
    COFD[3366] = -3.29752510E-01;
    COFD[3367] = 1.43850275E-02;
    COFD[3368] = -1.95318173E+01;
    COFD[3369] = 4.90255048E+00;
    COFD[3370] = -4.17368501E-01;
    COFD[3371] = 1.79287358E-02;
    COFD[3372] = -1.95318173E+01;
    COFD[3373] = 4.90255048E+00;
    COFD[3374] = -4.17368501E-01;
    COFD[3375] = 1.79287358E-02;
    COFD[3376] = -1.95552142E+01;
    COFD[3377] = 4.90255048E+00;
    COFD[3378] = -4.17368501E-01;
    COFD[3379] = 1.79287358E-02;
    COFD[3380] = -1.94823660E+01;
    COFD[3381] = 4.87333294E+00;
    COFD[3382] = -4.13769241E-01;
    COFD[3383] = 1.77802244E-02;
    COFD[3384] = -1.83538377E+01;
    COFD[3385] = 4.42828044E+00;
    COFD[3386] = -3.60417833E-01;
    COFD[3387] = 1.56455103E-02;
    COFD[3388] = -2.14326461E+01;
    COFD[3389] = 5.41729961E+00;
    COFD[3390] = -4.73400281E-01;
    COFD[3391] = 1.99269567E-02;
    COFD[3392] = -2.19313638E+01;
    COFD[3393] = 5.41841631E+00;
    COFD[3394] = -4.46818971E-01;
    COFD[3395] = 1.77127652E-02;
    COFD[3396] = -2.19399793E+01;
    COFD[3397] = 5.41841631E+00;
    COFD[3398] = -4.46818971E-01;
    COFD[3399] = 1.77127652E-02;
    COFD[3400] = -2.21157340E+01;
    COFD[3401] = 5.55072945E+00;
    COFD[3402] = -4.72525345E-01;
    COFD[3403] = 1.91674202E-02;
    COFD[3404] = -2.21157340E+01;
    COFD[3405] = 5.55072945E+00;
    COFD[3406] = -4.72525345E-01;
    COFD[3407] = 1.91674202E-02;
    COFD[3408] = -2.20116855E+01;
    COFD[3409] = 5.45178028E+00;
    COFD[3410] = -4.52847771E-01;
    COFD[3411] = 1.80418544E-02;
    COFD[3412] = -2.10837327E+01;
    COFD[3413] = 5.34286099E+00;
    COFD[3414] = -4.68100992E-01;
    COFD[3415] = 1.98731399E-02;
    COFD[3416] = -2.10944088E+01;
    COFD[3417] = 5.34286099E+00;
    COFD[3418] = -4.68100992E-01;
    COFD[3419] = 1.98731399E-02;
    COFD[3420] = -2.11044965E+01;
    COFD[3421] = 5.34286099E+00;
    COFD[3422] = -4.68100992E-01;
    COFD[3423] = 1.98731399E-02;
    COFD[3424] = -2.18273547E+01;
    COFD[3425] = 5.55753905E+00;
    COFD[3426] = -4.88136714E-01;
    COFD[3427] = 2.04294957E-02;
    COFD[3428] = -2.15660171E+01;
    COFD[3429] = 5.44803850E+00;
    COFD[3430] = -4.76610560E-01;
    COFD[3431] = 2.00355294E-02;
    COFD[3432] = -2.15746136E+01;
    COFD[3433] = 5.44803850E+00;
    COFD[3434] = -4.76610560E-01;
    COFD[3435] = 2.00355294E-02;
    COFD[3436] = -1.96408127E+01;
    COFD[3437] = 4.95923807E+00;
    COFD[3438] = -4.24176182E-01;
    COFD[3439] = 1.82020215E-02;
    COFD[3440] = -2.22159630E+01;
    COFD[3441] = 5.51722375E+00;
    COFD[3442] = -4.66081431E-01;
    COFD[3443] = 1.88044011E-02;
    COFD[3444] = -2.22159630E+01;
    COFD[3445] = 5.51722375E+00;
    COFD[3446] = -4.66081431E-01;
    COFD[3447] = 1.88044011E-02;
    COFD[3448] = -1.83249299E+01;
    COFD[3449] = 4.42045763E+00;
    COFD[3450] = -3.59451578E-01;
    COFD[3451] = 1.56056164E-02;
    COFD[3452] = -1.95026789E+01;
    COFD[3453] = 4.83189721E+00;
    COFD[3454] = -4.08932249E-01;
    COFD[3455] = 1.75924650E-02;
    COFD[3456] = -1.22004324E+01;
    COFD[3457] = 2.80725489E+00;
    COFD[3458] = -1.54291406E-01;
    COFD[3459] = 6.88290911E-03;
    COFD[3460] = -1.54460820E+01;
    COFD[3461] = 4.26819983E+00;
    COFD[3462] = -3.40766379E-01;
    COFD[3463] = 1.48393361E-02;
    COFD[3464] = -1.49500357E+01;
    COFD[3465] = 3.52327209E+00;
    COFD[3466] = -2.46286208E-01;
    COFD[3467] = 1.08285963E-02;
    COFD[3468] = -1.64169433E+01;
    COFD[3469] = 3.89309916E+00;
    COFD[3470] = -2.93528188E-01;
    COFD[3471] = 1.28463177E-02;
    COFD[3472] = -1.49718233E+01;
    COFD[3473] = 3.52327209E+00;
    COFD[3474] = -2.46286208E-01;
    COFD[3475] = 1.08285963E-02;
    COFD[3476] = -2.10440675E+01;
    COFD[3477] = 5.59806282E+00;
    COFD[3478] = -4.87109535E-01;
    COFD[3479] = 2.01370226E-02;
    COFD[3480] = -1.64255964E+01;
    COFD[3481] = 3.89309916E+00;
    COFD[3482] = -2.93528188E-01;
    COFD[3483] = 1.28463177E-02;
    COFD[3484] = -1.64338757E+01;
    COFD[3485] = 3.89309916E+00;
    COFD[3486] = -2.93528188E-01;
    COFD[3487] = 1.28463177E-02;
    COFD[3488] = -1.46907028E+01;
    COFD[3489] = 3.39229020E+00;
    COFD[3490] = -2.29520232E-01;
    COFD[3491] = 1.01114311E-02;
    COFD[3492] = -1.48738066E+01;
    COFD[3493] = 3.52327209E+00;
    COFD[3494] = -2.46286208E-01;
    COFD[3495] = 1.08285963E-02;
    COFD[3496] = -1.72572042E+01;
    COFD[3497] = 4.26063341E+00;
    COFD[3498] = -3.39848064E-01;
    COFD[3499] = 1.48021313E-02;
    COFD[3500] = -1.72572042E+01;
    COFD[3501] = 4.26063341E+00;
    COFD[3502] = -3.39848064E-01;
    COFD[3503] = 1.48021313E-02;
    COFD[3504] = -1.72828302E+01;
    COFD[3505] = 4.26063341E+00;
    COFD[3506] = -3.39848064E-01;
    COFD[3507] = 1.48021313E-02;
    COFD[3508] = -1.72316148E+01;
    COFD[3509] = 4.24011069E+00;
    COFD[3510] = -3.37339810E-01;
    COFD[3511] = 1.46996679E-02;
    COFD[3512] = -1.60261675E+01;
    COFD[3513] = 3.73312045E+00;
    COFD[3514] = -2.72579779E-01;
    COFD[3515] = 1.19290272E-02;
    COFD[3516] = -1.94485982E+01;
    COFD[3517] = 4.91446566E+00;
    COFD[3518] = -4.18837152E-01;
    COFD[3519] = 1.79893537E-02;
    COFD[3520] = -2.14204185E+01;
    COFD[3521] = 5.59268435E+00;
    COFD[3522] = -4.91232974E-01;
    COFD[3523] = 2.05064746E-02;
    COFD[3524] = -2.14303479E+01;
    COFD[3525] = 5.59268435E+00;
    COFD[3526] = -4.91232974E-01;
    COFD[3527] = 2.05064746E-02;
    COFD[3528] = -2.09241647E+01;
    COFD[3529] = 5.42316225E+00;
    COFD[3530] = -4.73702801E-01;
    COFD[3531] = 1.99217718E-02;
    COFD[3532] = -2.09241647E+01;
    COFD[3533] = 5.42316225E+00;
    COFD[3534] = -4.73702801E-01;
    COFD[3535] = 1.99217718E-02;
    COFD[3536] = -2.13796303E+01;
    COFD[3537] = 5.56978987E+00;
    COFD[3538] = -4.89141980E-01;
    COFD[3539] = 2.04499210E-02;
    COFD[3540] = -1.88524385E+01;
    COFD[3541] = 4.72476764E+00;
    COFD[3542] = -3.96306836E-01;
    COFD[3543] = 1.70964541E-02;
    COFD[3544] = -1.88646070E+01;
    COFD[3545] = 4.72476764E+00;
    COFD[3546] = -3.96306836E-01;
    COFD[3547] = 1.70964541E-02;
    COFD[3548] = -1.88761387E+01;
    COFD[3549] = 4.72476764E+00;
    COFD[3550] = -3.96306836E-01;
    COFD[3551] = 1.70964541E-02;
    COFD[3552] = -1.99081556E+01;
    COFD[3553] = 5.09311649E+00;
    COFD[3554] = -4.39965178E-01;
    COFD[3555] = 1.88238537E-02;
    COFD[3556] = -1.96309042E+01;
    COFD[3557] = 4.95923807E+00;
    COFD[3558] = -4.24176182E-01;
    COFD[3559] = 1.82020215E-02;
    COFD[3560] = -1.96408127E+01;
    COFD[3561] = 4.95923807E+00;
    COFD[3562] = -4.24176182E-01;
    COFD[3563] = 1.82020215E-02;
    COFD[3564] = -1.72414862E+01;
    COFD[3565] = 4.29808578E+00;
    COFD[3566] = -3.44235570E-01;
    COFD[3567] = 1.49727727E-02;
    COFD[3568] = -2.12621914E+01;
    COFD[3569] = 5.47935225E+00;
    COFD[3570] = -4.80056796E-01;
    COFD[3571] = 2.01607180E-02;
    COFD[3572] = -2.12621914E+01;
    COFD[3573] = 5.47935225E+00;
    COFD[3574] = -4.80056796E-01;
    COFD[3575] = 2.01607180E-02;
    COFD[3576] = -1.59884305E+01;
    COFD[3577] = 3.72220402E+00;
    COFD[3578] = -2.71150591E-01;
    COFD[3579] = 1.18665265E-02;
    COFD[3580] = -1.72570607E+01;
    COFD[3581] = 4.19757624E+00;
    COFD[3582] = -3.32087529E-01;
    COFD[3583] = 1.44827462E-02;
    COFD[3584] = -1.57034851E+01;
    COFD[3585] = 3.93614244E+00;
    COFD[3586] = -2.99111497E-01;
    COFD[3587] = 1.30888229E-02;
    COFD[3588] = -1.94688688E+01;
    COFD[3589] = 5.43830787E+00;
    COFD[3590] = -4.75472880E-01;
    COFD[3591] = 1.99909996E-02;
    COFD[3592] = -1.90883268E+01;
    COFD[3593] = 4.84384483E+00;
    COFD[3594] = -4.10265575E-01;
    COFD[3595] = 1.76414287E-02;
    COFD[3596] = -2.05184870E+01;
    COFD[3597] = 5.18417470E+00;
    COFD[3598] = -4.49491573E-01;
    COFD[3599] = 1.91438508E-02;
    COFD[3600] = -1.91102652E+01;
    COFD[3601] = 4.84384483E+00;
    COFD[3602] = -4.10265575E-01;
    COFD[3603] = 1.76414287E-02;
    COFD[3604] = -1.87383952E+01;
    COFD[3605] = 3.96926341E+00;
    COFD[3606] = -2.16412264E-01;
    COFD[3607] = 6.06012078E-03;
    COFD[3608] = -2.05272328E+01;
    COFD[3609] = 5.18417470E+00;
    COFD[3610] = -4.49491573E-01;
    COFD[3611] = 1.91438508E-02;
    COFD[3612] = -2.05356023E+01;
    COFD[3613] = 5.18417470E+00;
    COFD[3614] = -4.49491573E-01;
    COFD[3615] = 1.91438508E-02;
    COFD[3616] = -1.87688110E+01;
    COFD[3617] = 4.71729964E+00;
    COFD[3618] = -3.95432573E-01;
    COFD[3619] = 1.70623691E-02;
    COFD[3620] = -1.90116191E+01;
    COFD[3621] = 4.84384483E+00;
    COFD[3622] = -4.10265575E-01;
    COFD[3623] = 1.76414287E-02;
    COFD[3624] = -2.11349086E+01;
    COFD[3625] = 5.42846112E+00;
    COFD[3626] = -4.74321870E-01;
    COFD[3627] = 1.99459749E-02;
    COFD[3628] = -2.11349086E+01;
    COFD[3629] = 5.42846112E+00;
    COFD[3630] = -4.74321870E-01;
    COFD[3631] = 1.99459749E-02;
    COFD[3632] = -2.11606963E+01;
    COFD[3633] = 5.42846112E+00;
    COFD[3634] = -4.74321870E-01;
    COFD[3635] = 1.99459749E-02;
    COFD[3636] = -2.11309207E+01;
    COFD[3637] = 5.41773516E+00;
    COFD[3638] = -4.73414338E-01;
    COFD[3639] = 1.99258685E-02;
    COFD[3640] = -2.02922701E+01;
    COFD[3641] = 5.11106992E+00;
    COFD[3642] = -4.42047129E-01;
    COFD[3643] = 1.89042990E-02;
    COFD[3644] = -2.22116706E+01;
    COFD[3645] = 5.54251230E+00;
    COFD[3646] = -4.70946314E-01;
    COFD[3647] = 1.90785869E-02;
    COFD[3648] = -2.00915040E+01;
    COFD[3649] = 4.41511629E+00;
    COFD[3650] = -2.84086963E-01;
    COFD[3651] = 9.37586971E-03;
    COFD[3652] = -2.01015340E+01;
    COFD[3653] = 4.41511629E+00;
    COFD[3654] = -2.84086963E-01;
    COFD[3655] = 9.37586971E-03;
    COFD[3656] = -2.09222454E+01;
    COFD[3657] = 4.82184721E+00;
    COFD[3658] = -3.48128875E-01;
    COFD[3659] = 1.25918978E-02;
    COFD[3660] = -2.09222454E+01;
    COFD[3661] = 4.82184721E+00;
    COFD[3662] = -3.48128875E-01;
    COFD[3663] = 1.25918978E-02;
    COFD[3664] = -2.03036402E+01;
    COFD[3665] = 4.50250781E+00;
    COFD[3666] = -2.97622106E-01;
    COFD[3667] = 1.00481473E-02;
    COFD[3668] = -2.20942491E+01;
    COFD[3669] = 5.58360799E+00;
    COFD[3670] = -4.82701436E-01;
    COFD[3671] = 1.98437922E-02;
    COFD[3672] = -2.21065306E+01;
    COFD[3673] = 5.58360799E+00;
    COFD[3674] = -4.82701436E-01;
    COFD[3675] = 1.98437922E-02;
    COFD[3676] = -2.21181720E+01;
    COFD[3677] = 5.58360799E+00;
    COFD[3678] = -4.82701436E-01;
    COFD[3679] = 1.98437922E-02;
    COFD[3680] = -2.20453723E+01;
    COFD[3681] = 5.44448440E+00;
    COFD[3682] = -4.51529024E-01;
    COFD[3683] = 1.79698119E-02;
    COFD[3684] = -2.22059540E+01;
    COFD[3685] = 5.51722375E+00;
    COFD[3686] = -4.66081431E-01;
    COFD[3687] = 1.88044011E-02;
    COFD[3688] = -2.22159630E+01;
    COFD[3689] = 5.51722375E+00;
    COFD[3690] = -4.66081431E-01;
    COFD[3691] = 1.88044011E-02;
    COFD[3692] = -2.12621914E+01;
    COFD[3693] = 5.47935225E+00;
    COFD[3694] = -4.80056796E-01;
    COFD[3695] = 2.01607180E-02;
    COFD[3696] = -2.09002742E+01;
    COFD[3697] = 4.72895031E+00;
    COFD[3698] = -3.33332771E-01;
    COFD[3699] = 1.18431478E-02;
    COFD[3700] = -2.09002742E+01;
    COFD[3701] = 4.72895031E+00;
    COFD[3702] = -3.33332771E-01;
    COFD[3703] = 1.18431478E-02;
    COFD[3704] = -2.02646611E+01;
    COFD[3705] = 5.10426133E+00;
    COFD[3706] = -4.41256919E-01;
    COFD[3707] = 1.88737290E-02;
    COFD[3708] = -2.12450965E+01;
    COFD[3709] = 5.40444222E+00;
    COFD[3710] = -4.72708609E-01;
    COFD[3711] = 1.99362392E-02;
    COFD[3712] = -1.57034851E+01;
    COFD[3713] = 3.93614244E+00;
    COFD[3714] = -2.99111497E-01;
    COFD[3715] = 1.30888229E-02;
    COFD[3716] = -1.94688688E+01;
    COFD[3717] = 5.43830787E+00;
    COFD[3718] = -4.75472880E-01;
    COFD[3719] = 1.99909996E-02;
    COFD[3720] = -1.90883268E+01;
    COFD[3721] = 4.84384483E+00;
    COFD[3722] = -4.10265575E-01;
    COFD[3723] = 1.76414287E-02;
    COFD[3724] = -2.05184870E+01;
    COFD[3725] = 5.18417470E+00;
    COFD[3726] = -4.49491573E-01;
    COFD[3727] = 1.91438508E-02;
    COFD[3728] = -1.91102652E+01;
    COFD[3729] = 4.84384483E+00;
    COFD[3730] = -4.10265575E-01;
    COFD[3731] = 1.76414287E-02;
    COFD[3732] = -1.87383952E+01;
    COFD[3733] = 3.96926341E+00;
    COFD[3734] = -2.16412264E-01;
    COFD[3735] = 6.06012078E-03;
    COFD[3736] = -2.05272328E+01;
    COFD[3737] = 5.18417470E+00;
    COFD[3738] = -4.49491573E-01;
    COFD[3739] = 1.91438508E-02;
    COFD[3740] = -2.05356023E+01;
    COFD[3741] = 5.18417470E+00;
    COFD[3742] = -4.49491573E-01;
    COFD[3743] = 1.91438508E-02;
    COFD[3744] = -1.87688110E+01;
    COFD[3745] = 4.71729964E+00;
    COFD[3746] = -3.95432573E-01;
    COFD[3747] = 1.70623691E-02;
    COFD[3748] = -1.90116191E+01;
    COFD[3749] = 4.84384483E+00;
    COFD[3750] = -4.10265575E-01;
    COFD[3751] = 1.76414287E-02;
    COFD[3752] = -2.11349086E+01;
    COFD[3753] = 5.42846112E+00;
    COFD[3754] = -4.74321870E-01;
    COFD[3755] = 1.99459749E-02;
    COFD[3756] = -2.11349086E+01;
    COFD[3757] = 5.42846112E+00;
    COFD[3758] = -4.74321870E-01;
    COFD[3759] = 1.99459749E-02;
    COFD[3760] = -2.11606963E+01;
    COFD[3761] = 5.42846112E+00;
    COFD[3762] = -4.74321870E-01;
    COFD[3763] = 1.99459749E-02;
    COFD[3764] = -2.11309207E+01;
    COFD[3765] = 5.41773516E+00;
    COFD[3766] = -4.73414338E-01;
    COFD[3767] = 1.99258685E-02;
    COFD[3768] = -2.02922701E+01;
    COFD[3769] = 5.11106992E+00;
    COFD[3770] = -4.42047129E-01;
    COFD[3771] = 1.89042990E-02;
    COFD[3772] = -2.22116706E+01;
    COFD[3773] = 5.54251230E+00;
    COFD[3774] = -4.70946314E-01;
    COFD[3775] = 1.90785869E-02;
    COFD[3776] = -2.00915040E+01;
    COFD[3777] = 4.41511629E+00;
    COFD[3778] = -2.84086963E-01;
    COFD[3779] = 9.37586971E-03;
    COFD[3780] = -2.01015340E+01;
    COFD[3781] = 4.41511629E+00;
    COFD[3782] = -2.84086963E-01;
    COFD[3783] = 9.37586971E-03;
    COFD[3784] = -2.09222454E+01;
    COFD[3785] = 4.82184721E+00;
    COFD[3786] = -3.48128875E-01;
    COFD[3787] = 1.25918978E-02;
    COFD[3788] = -2.09222454E+01;
    COFD[3789] = 4.82184721E+00;
    COFD[3790] = -3.48128875E-01;
    COFD[3791] = 1.25918978E-02;
    COFD[3792] = -2.03036402E+01;
    COFD[3793] = 4.50250781E+00;
    COFD[3794] = -2.97622106E-01;
    COFD[3795] = 1.00481473E-02;
    COFD[3796] = -2.20942491E+01;
    COFD[3797] = 5.58360799E+00;
    COFD[3798] = -4.82701436E-01;
    COFD[3799] = 1.98437922E-02;
    COFD[3800] = -2.21065306E+01;
    COFD[3801] = 5.58360799E+00;
    COFD[3802] = -4.82701436E-01;
    COFD[3803] = 1.98437922E-02;
    COFD[3804] = -2.21181720E+01;
    COFD[3805] = 5.58360799E+00;
    COFD[3806] = -4.82701436E-01;
    COFD[3807] = 1.98437922E-02;
    COFD[3808] = -2.20453723E+01;
    COFD[3809] = 5.44448440E+00;
    COFD[3810] = -4.51529024E-01;
    COFD[3811] = 1.79698119E-02;
    COFD[3812] = -2.22059540E+01;
    COFD[3813] = 5.51722375E+00;
    COFD[3814] = -4.66081431E-01;
    COFD[3815] = 1.88044011E-02;
    COFD[3816] = -2.22159630E+01;
    COFD[3817] = 5.51722375E+00;
    COFD[3818] = -4.66081431E-01;
    COFD[3819] = 1.88044011E-02;
    COFD[3820] = -2.12621914E+01;
    COFD[3821] = 5.47935225E+00;
    COFD[3822] = -4.80056796E-01;
    COFD[3823] = 2.01607180E-02;
    COFD[3824] = -2.09002742E+01;
    COFD[3825] = 4.72895031E+00;
    COFD[3826] = -3.33332771E-01;
    COFD[3827] = 1.18431478E-02;
    COFD[3828] = -2.09002742E+01;
    COFD[3829] = 4.72895031E+00;
    COFD[3830] = -3.33332771E-01;
    COFD[3831] = 1.18431478E-02;
    COFD[3832] = -2.02646611E+01;
    COFD[3833] = 5.10426133E+00;
    COFD[3834] = -4.41256919E-01;
    COFD[3835] = 1.88737290E-02;
    COFD[3836] = -2.12450965E+01;
    COFD[3837] = 5.40444222E+00;
    COFD[3838] = -4.72708609E-01;
    COFD[3839] = 1.99362392E-02;
    COFD[3840] = -1.16906297E+01;
    COFD[3841] = 2.47469981E+00;
    COFD[3842] = -1.10436257E-01;
    COFD[3843] = 4.95273813E-03;
    COFD[3844] = -1.42894441E+01;
    COFD[3845] = 3.67490723E+00;
    COFD[3846] = -2.65114792E-01;
    COFD[3847] = 1.16092671E-02;
    COFD[3848] = -1.40756935E+01;
    COFD[3849] = 3.07549274E+00;
    COFD[3850] = -1.88889344E-01;
    COFD[3851] = 8.37152866E-03;
    COFD[3852] = -1.52414485E+01;
    COFD[3853] = 3.35922578E+00;
    COFD[3854] = -2.25181399E-01;
    COFD[3855] = 9.92132878E-03;
    COFD[3856] = -1.40949196E+01;
    COFD[3857] = 3.07549274E+00;
    COFD[3858] = -1.88889344E-01;
    COFD[3859] = 8.37152866E-03;
    COFD[3860] = -2.10643259E+01;
    COFD[3861] = 5.53614847E+00;
    COFD[3862] = -4.86046736E-01;
    COFD[3863] = 2.03659188E-02;
    COFD[3864] = -1.52486273E+01;
    COFD[3865] = 3.35922578E+00;
    COFD[3866] = -2.25181399E-01;
    COFD[3867] = 9.92132878E-03;
    COFD[3868] = -1.52554761E+01;
    COFD[3869] = 3.35922578E+00;
    COFD[3870] = -2.25181399E-01;
    COFD[3871] = 9.92132878E-03;
    COFD[3872] = -1.38661480E+01;
    COFD[3873] = 2.97137588E+00;
    COFD[3874] = -1.75491257E-01;
    COFD[3875] = 7.79646773E-03;
    COFD[3876] = -1.40076852E+01;
    COFD[3877] = 3.07549274E+00;
    COFD[3878] = -1.88889344E-01;
    COFD[3879] = 8.37152866E-03;
    COFD[3880] = -1.59404882E+01;
    COFD[3881] = 3.66853818E+00;
    COFD[3882] = -2.64346221E-01;
    COFD[3883] = 1.15784613E-02;
    COFD[3884] = -1.59404882E+01;
    COFD[3885] = 3.66853818E+00;
    COFD[3886] = -2.64346221E-01;
    COFD[3887] = 1.15784613E-02;
    COFD[3888] = -1.59633387E+01;
    COFD[3889] = 3.66853818E+00;
    COFD[3890] = -2.64346221E-01;
    COFD[3891] = 1.15784613E-02;
    COFD[3892] = -1.59327297E+01;
    COFD[3893] = 3.65620899E+00;
    COFD[3894] = -2.62933804E-01;
    COFD[3895] = 1.15253223E-02;
    COFD[3896] = -1.50031687E+01;
    COFD[3897] = 3.26223357E+00;
    COFD[3898] = -2.12746642E-01;
    COFD[3899] = 9.38912883E-03;
    COFD[3900] = -1.81432461E+01;
    COFD[3901] = 4.37565431E+00;
    COFD[3902] = -3.53906025E-01;
    COFD[3903] = 1.53760786E-02;
    COFD[3904] = -2.04750581E+01;
    COFD[3905] = 5.23112374E+00;
    COFD[3906] = -4.54967682E-01;
    COFD[3907] = 1.93570423E-02;
    COFD[3908] = -2.04833713E+01;
    COFD[3909] = 5.23112374E+00;
    COFD[3910] = -4.54967682E-01;
    COFD[3911] = 1.93570423E-02;
    COFD[3912] = -2.02268902E+01;
    COFD[3913] = 5.13632093E+00;
    COFD[3914] = -4.44839124E-01;
    COFD[3915] = 1.90058354E-02;
    COFD[3916] = -2.02268902E+01;
    COFD[3917] = 5.13632093E+00;
    COFD[3918] = -4.44839124E-01;
    COFD[3919] = 1.90058354E-02;
    COFD[3920] = -2.03844252E+01;
    COFD[3921] = 5.18856872E+00;
    COFD[3922] = -4.50001829E-01;
    COFD[3923] = 1.91636142E-02;
    COFD[3924] = -1.75898751E+01;
    COFD[3925] = 4.19171952E+00;
    COFD[3926] = -3.31354810E-01;
    COFD[3927] = 1.44520623E-02;
    COFD[3928] = -1.76002031E+01;
    COFD[3929] = 4.19171952E+00;
    COFD[3930] = -3.31354810E-01;
    COFD[3931] = 1.44520623E-02;
    COFD[3932] = -1.76099552E+01;
    COFD[3933] = 4.19171952E+00;
    COFD[3934] = -3.31354810E-01;
    COFD[3935] = 1.44520623E-02;
    COFD[3936] = -1.85864144E+01;
    COFD[3937] = 4.54915847E+00;
    COFD[3938] = -3.75000738E-01;
    COFD[3939] = 1.62324821E-02;
    COFD[3940] = -1.83166353E+01;
    COFD[3941] = 4.42045763E+00;
    COFD[3942] = -3.59451578E-01;
    COFD[3943] = 1.56056164E-02;
    COFD[3944] = -1.83249299E+01;
    COFD[3945] = 4.42045763E+00;
    COFD[3946] = -3.59451578E-01;
    COFD[3947] = 1.56056164E-02;
    COFD[3948] = -1.59884305E+01;
    COFD[3949] = 3.72220402E+00;
    COFD[3950] = -2.71150591E-01;
    COFD[3951] = 1.18665265E-02;
    COFD[3952] = -2.02646611E+01;
    COFD[3953] = 5.10426133E+00;
    COFD[3954] = -4.41256919E-01;
    COFD[3955] = 1.88737290E-02;
    COFD[3956] = -2.02646611E+01;
    COFD[3957] = 5.10426133E+00;
    COFD[3958] = -4.41256919E-01;
    COFD[3959] = 1.88737290E-02;
    COFD[3960] = -1.49828430E+01;
    COFD[3961] = 3.25781069E+00;
    COFD[3962] = -2.12199367E-01;
    COFD[3963] = 9.36657283E-03;
    COFD[3964] = -1.59877782E+01;
    COFD[3965] = 3.63340763E+00;
    COFD[3966] = -2.60307961E-01;
    COFD[3967] = 1.14256954E-02;
    COFD[3968] = -1.23130152E+01;
    COFD[3969] = 2.74418790E+00;
    COFD[3970] = -1.46230156E-01;
    COFD[3971] = 6.53948886E-03;
    COFD[3972] = -1.54738604E+01;
    COFD[3973] = 4.15765300E+00;
    COFD[3974] = -3.27126237E-01;
    COFD[3975] = 1.42762611E-02;
    COFD[3976] = -1.49610527E+01;
    COFD[3977] = 3.41988961E+00;
    COFD[3978] = -2.33128386E-01;
    COFD[3979] = 1.02689994E-02;
    COFD[3980] = -1.62380075E+01;
    COFD[3981] = 3.72612300E+00;
    COFD[3982] = -2.71663673E-01;
    COFD[3983] = 1.18889643E-02;
    COFD[3984] = -1.49826725E+01;
    COFD[3985] = 3.41988961E+00;
    COFD[3986] = -2.33128386E-01;
    COFD[3987] = 1.02689994E-02;
    COFD[3988] = -2.12755888E+01;
    COFD[3989] = 5.60381989E+00;
    COFD[3990] = -4.91225459E-01;
    COFD[3991] = 2.04487844E-02;
    COFD[3992] = -1.62465583E+01;
    COFD[3993] = 3.72612300E+00;
    COFD[3994] = -2.71663673E-01;
    COFD[3995] = 1.18889643E-02;
    COFD[3996] = -1.62547381E+01;
    COFD[3997] = 3.72612300E+00;
    COFD[3998] = -2.71663673E-01;
    COFD[3999] = 1.18889643E-02;
    COFD[4000] = -1.46461881E+01;
    COFD[4001] = 3.27505697E+00;
    COFD[4002] = -2.14306851E-01;
    COFD[4003] = 9.45219335E-03;
    COFD[4004] = -1.48853569E+01;
    COFD[4005] = 3.41988961E+00;
    COFD[4006] = -2.33128386E-01;
    COFD[4007] = 1.02689994E-02;
    COFD[4008] = -1.71942502E+01;
    COFD[4009] = 4.14993355E+00;
    COFD[4010] = -3.26168062E-01;
    COFD[4011] = 1.42364115E-02;
    COFD[4012] = -1.71942502E+01;
    COFD[4013] = 4.14993355E+00;
    COFD[4014] = -3.26168062E-01;
    COFD[4015] = 1.42364115E-02;
    COFD[4016] = -1.72196961E+01;
    COFD[4017] = 4.14993355E+00;
    COFD[4018] = -3.26168062E-01;
    COFD[4019] = 1.42364115E-02;
    COFD[4020] = -1.71754154E+01;
    COFD[4021] = 4.13131681E+00;
    COFD[4022] = -3.23897559E-01;
    COFD[4023] = 1.41438222E-02;
    COFD[4024] = -1.60074211E+01;
    COFD[4025] = 3.63723937E+00;
    COFD[4026] = -2.60754222E-01;
    COFD[4027] = 1.14428814E-02;
    COFD[4028] = -1.93483692E+01;
    COFD[4029] = 4.79506290E+00;
    COFD[4030] = -4.04621659E-01;
    COFD[4031] = 1.74244230E-02;
    COFD[4032] = -2.14255087E+01;
    COFD[4033] = 5.52240865E+00;
    COFD[4034] = -4.84699537E-01;
    COFD[4035] = 2.03247833E-02;
    COFD[4036] = -2.14353267E+01;
    COFD[4037] = 5.52240865E+00;
    COFD[4038] = -4.84699537E-01;
    COFD[4039] = 2.03247833E-02;
    COFD[4040] = -2.10026861E+01;
    COFD[4041] = 5.38326647E+00;
    COFD[4042] = -4.71201048E-01;
    COFD[4043] = 1.99207516E-02;
    COFD[4044] = -2.10026861E+01;
    COFD[4045] = 5.38326647E+00;
    COFD[4046] = -4.71201048E-01;
    COFD[4047] = 1.99207516E-02;
    COFD[4048] = -2.13502344E+01;
    COFD[4049] = 5.48617067E+00;
    COFD[4050] = -4.80816776E-01;
    COFD[4051] = 2.01887868E-02;
    COFD[4052] = -1.87596895E+01;
    COFD[4053] = 4.61078776E+00;
    COFD[4054] = -3.82625667E-01;
    COFD[4055] = 1.65478601E-02;
    COFD[4056] = -1.87717330E+01;
    COFD[4057] = 4.61078776E+00;
    COFD[4058] = -3.82625667E-01;
    COFD[4059] = 1.65478601E-02;
    COFD[4060] = -1.87831434E+01;
    COFD[4061] = 4.61078776E+00;
    COFD[4062] = -3.82625667E-01;
    COFD[4063] = 1.65478601E-02;
    COFD[4064] = -1.98040322E+01;
    COFD[4065] = 4.97569695E+00;
    COFD[4066] = -4.26123307E-01;
    COFD[4067] = 1.82788664E-02;
    COFD[4068] = -1.94928815E+01;
    COFD[4069] = 4.83189721E+00;
    COFD[4070] = -4.08932249E-01;
    COFD[4071] = 1.75924650E-02;
    COFD[4072] = -1.95026789E+01;
    COFD[4073] = 4.83189721E+00;
    COFD[4074] = -4.08932249E-01;
    COFD[4075] = 1.75924650E-02;
    COFD[4076] = -1.72570607E+01;
    COFD[4077] = 4.19757624E+00;
    COFD[4078] = -3.32087529E-01;
    COFD[4079] = 1.44827462E-02;
    COFD[4080] = -2.12450965E+01;
    COFD[4081] = 5.40444222E+00;
    COFD[4082] = -4.72708609E-01;
    COFD[4083] = 1.99362392E-02;
    COFD[4084] = -2.12450965E+01;
    COFD[4085] = 5.40444222E+00;
    COFD[4086] = -4.72708609E-01;
    COFD[4087] = 1.99362392E-02;
    COFD[4088] = -1.59877782E+01;
    COFD[4089] = 3.63340763E+00;
    COFD[4090] = -2.60307961E-01;
    COFD[4091] = 1.14256954E-02;
    COFD[4092] = -1.72273911E+01;
    COFD[4093] = 4.09361913E+00;
    COFD[4094] = -3.19258125E-01;
    COFD[4095] = 1.39526981E-02;
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
    COFTD[28] = 4.29789463E-01;
    COFTD[29] = 1.21313199E-04;
    COFTD[30] = -7.73071792E-08;
    COFTD[31] = 1.53234639E-11;
    COFTD[32] = 3.82245475E-01;
    COFTD[33] = 1.47167068E-05;
    COFTD[34] = -9.75995257E-09;
    COFTD[35] = 2.83217152E-12;
    COFTD[36] = 3.83439056E-01;
    COFTD[37] = 3.62717894E-05;
    COFTD[38] = -2.40281409E-08;
    COFTD[39] = 5.52792966E-12;
    COFTD[40] = 3.24747031E-01;
    COFTD[41] = 1.77798548E-04;
    COFTD[42] = -1.08934732E-07;
    COFTD[43] = 2.03595881E-11;
    COFTD[44] = 3.24747031E-01;
    COFTD[45] = 1.77798548E-04;
    COFTD[46] = -1.08934732E-07;
    COFTD[47] = 2.03595881E-11;
    COFTD[48] = 3.31191185E-01;
    COFTD[49] = 1.81326714E-04;
    COFTD[50] = -1.11096391E-07;
    COFTD[51] = 2.07635959E-11;
    COFTD[52] = 3.39557243E-01;
    COFTD[53] = 1.79335036E-04;
    COFTD[54] = -1.10135705E-07;
    COFTD[55] = 2.06427239E-11;
    COFTD[56] = 4.30605547E-01;
    COFTD[57] = 9.35961902E-05;
    COFTD[58] = -6.03983623E-08;
    COFTD[59] = 1.23115170E-11;
    COFTD[60] = 2.93191523E-01;
    COFTD[61] = 4.01430006E-04;
    COFTD[62] = -2.30705763E-07;
    COFTD[63] = 4.05176586E-11;
    COFTD[64] = 1.22119780E-01;
    COFTD[65] = 6.18373616E-04;
    COFTD[66] = -3.28422593E-07;
    COFTD[67] = 5.44603522E-11;
    COFTD[68] = 1.22693382E-01;
    COFTD[69] = 6.21278143E-04;
    COFTD[70] = -3.29965208E-07;
    COFTD[71] = 5.47161548E-11;
    COFTD[72] = 1.40314191E-01;
    COFTD[73] = 6.01266129E-04;
    COFTD[74] = -3.21915137E-07;
    COFTD[75] = 5.36679068E-11;
    COFTD[76] = 1.40314191E-01;
    COFTD[77] = 6.01266129E-04;
    COFTD[78] = -3.21915137E-07;
    COFTD[79] = 5.36679068E-11;
    COFTD[80] = 1.31424053E-01;
    COFTD[81] = 6.16429134E-04;
    COFTD[82] = -3.28571348E-07;
    COFTD[83] = 5.46153434E-11;
    COFTD[84] = 3.03701584E-01;
    COFTD[85] = 3.22476070E-04;
    COFTD[86] = -1.88701794E-07;
    COFTD[87] = 3.36545091E-11;
    COFTD[88] = 3.05613225E-01;
    COFTD[89] = 3.24505886E-04;
    COFTD[90] = -1.89889572E-07;
    COFTD[91] = 3.38663465E-11;
    COFTD[92] = 3.07392263E-01;
    COFTD[93] = 3.26394901E-04;
    COFTD[94] = -1.90994958E-07;
    COFTD[95] = 3.40634894E-11;
    COFTD[96] = 2.49017478E-01;
    COFTD[97] = 4.29036573E-04;
    COFTD[98] = -2.42668617E-07;
    COFTD[99] = 4.20801371E-11;
    COFTD[100] = 2.72759599E-01;
    COFTD[101] = 3.94402719E-04;
    COFTD[102] = -2.25800520E-07;
    COFTD[103] = 3.95325634E-11;
    COFTD[104] = 2.74036956E-01;
    COFTD[105] = 3.96249742E-04;
    COFTD[106] = -2.26857964E-07;
    COFTD[107] = 3.97176979E-11;
    COFTD[108] = 3.86107464E-01;
    COFTD[109] = 2.28760446E-04;
    COFTD[110] = -1.39425040E-07;
    COFTD[111] = 2.58989754E-11;
    COFTD[112] = 1.59288984E-01;
    COFTD[113] = 6.02833801E-04;
    COFTD[114] = -3.24837576E-07;
    COFTD[115] = 5.43909010E-11;
    COFTD[116] = 1.59288984E-01;
    COFTD[117] = 6.02833801E-04;
    COFTD[118] = -3.24837576E-07;
    COFTD[119] = 5.43909010E-11;
    COFTD[120] = 4.31331269E-01;
    COFTD[121] = 9.20536800E-05;
    COFTD[122] = -5.94509616E-08;
    COFTD[123] = 1.21437993E-11;
    COFTD[124] = 4.01012808E-01;
    COFTD[125] = 1.97252826E-04;
    COFTD[126] = -1.21698146E-07;
    COFTD[127] = 2.29408847E-11;
    COFTD[128] = 1.44152190E-01;
    COFTD[129] = 7.99993584E-05;
    COFTD[130] = -4.89707442E-08;
    COFTD[131] = 9.14277269E-12;
    COFTD[132] = 0.00000000E+00;
    COFTD[133] = 0.00000000E+00;
    COFTD[134] = 0.00000000E+00;
    COFTD[135] = 0.00000000E+00;
    COFTD[136] = 2.35283119E-01;
    COFTD[137] = 4.65670599E-04;
    COFTD[138] = -2.60939824E-07;
    COFTD[139] = 4.49271822E-11;
    COFTD[140] = 1.79840299E-01;
    COFTD[141] = 6.01722902E-04;
    COFTD[142] = -3.26433894E-07;
    COFTD[143] = 5.49112302E-11;
    COFTD[144] = 2.37053352E-01;
    COFTD[145] = 4.69174231E-04;
    COFTD[146] = -2.62903094E-07;
    COFTD[147] = 4.52652072E-11;
    COFTD[148] = -1.74352698E-01;
    COFTD[149] = 8.62246873E-04;
    COFTD[150] = -3.79545489E-07;
    COFTD[151] = 5.60262093E-11;
    COFTD[152] = 1.80186965E-01;
    COFTD[153] = 6.02882805E-04;
    COFTD[154] = -3.27063140E-07;
    COFTD[155] = 5.50170790E-11;
    COFTD[156] = 1.80513677E-01;
    COFTD[157] = 6.03975942E-04;
    COFTD[158] = -3.27656165E-07;
    COFTD[159] = 5.51168351E-11;
    COFTD[160] = 2.49272491E-01;
    COFTD[161] = 4.08682510E-04;
    COFTD[162] = -2.31943878E-07;
    COFTD[163] = 4.03271405E-11;
    COFTD[164] = 2.28560867E-01;
    COFTD[165] = 4.52365967E-04;
    COFTD[166] = -2.53484536E-07;
    COFTD[167] = 4.36435719E-11;
    COFTD[168] = 9.90752318E-02;
    COFTD[169] = 6.44201384E-04;
    COFTD[170] = -3.38485953E-07;
    COFTD[171] = 5.57356746E-11;
    COFTD[172] = 9.90752318E-02;
    COFTD[173] = 6.44201384E-04;
    COFTD[174] = -3.38485953E-07;
    COFTD[175] = 5.57356746E-11;
    COFTD[176] = 1.00039110E-01;
    COFTD[177] = 6.50468660E-04;
    COFTD[178] = -3.41778999E-07;
    COFTD[179] = 5.62779132E-11;
    COFTD[180] = 1.05124122E-01;
    COFTD[181] = 6.50665957E-04;
    COFTD[182] = -3.42564538E-07;
    COFTD[183] = 5.64804120E-11;
    COFTD[184] = 2.00119897E-01;
    COFTD[185] = 5.64793704E-04;
    COFTD[186] = -3.09445484E-07;
    COFTD[187] = 5.24139335E-11;
    COFTD[188] = -2.00309448E-02;
    COFTD[189] = 8.50440115E-04;
    COFTD[190] = -4.21064468E-07;
    COFTD[191] = 6.67959710E-11;
    COFTD[192] = -1.60981264E-01;
    COFTD[193] = 9.03807572E-04;
    COFTD[194] = -4.06927941E-07;
    COFTD[195] = 6.09202254E-11;
    COFTD[196] = -1.61357564E-01;
    COFTD[197] = 9.05920260E-04;
    COFTD[198] = -4.07879153E-07;
    COFTD[199] = 6.10626290E-11;
    COFTD[200] = -1.31244519E-01;
    COFTD[201] = 9.03901384E-04;
    COFTD[202] = -4.17831507E-07;
    COFTD[203] = 6.35725667E-11;
    COFTD[204] = -1.31244519E-01;
    COFTD[205] = 9.03901384E-04;
    COFTD[206] = -4.17831507E-07;
    COFTD[207] = 6.35725667E-11;
    COFTD[208] = -1.56651581E-01;
    COFTD[209] = 9.09789751E-04;
    COFTD[210] = -4.11714242E-07;
    COFTD[211] = 6.18310893E-11;
    COFTD[212] = 1.62736132E-02;
    COFTD[213] = 7.87669911E-04;
    COFTD[214] = -3.97050662E-07;
    COFTD[215] = 6.36859622E-11;
    COFTD[216] = 1.63245097E-02;
    COFTD[217] = 7.90133388E-04;
    COFTD[218] = -3.98292458E-07;
    COFTD[219] = 6.38851432E-11;
    COFTD[220] = 1.63717489E-02;
    COFTD[221] = 7.92419842E-04;
    COFTD[222] = -3.99445020E-07;
    COFTD[223] = 6.40700113E-11;
    COFTD[224] = -5.08744745E-02;
    COFTD[225] = 8.54342586E-04;
    COFTD[226] = -4.15926453E-07;
    COFTD[227] = 6.53063261E-11;
    COFTD[228] = -2.71690558E-02;
    COFTD[229] = 8.37233133E-04;
    COFTD[230] = -4.12887636E-07;
    COFTD[231] = 6.53405197E-11;
    COFTD[232] = -2.72323768E-02;
    COFTD[233] = 8.39184413E-04;
    COFTD[234] = -4.13849924E-07;
    COFTD[235] = 6.54928043E-11;
    COFTD[236] = 9.86934401E-02;
    COFTD[237] = 7.20974863E-04;
    COFTD[238] = -3.77135221E-07;
    COFTD[239] = 6.19202579E-11;
    COFTD[240] = -1.41640506E-01;
    COFTD[241] = 9.21404324E-04;
    COFTD[242] = -4.23210110E-07;
    COFTD[243] = 6.41400322E-11;
    COFTD[244] = -1.41640506E-01;
    COFTD[245] = 9.21404324E-04;
    COFTD[246] = -4.23210110E-07;
    COFTD[247] = 6.41400322E-11;
    COFTD[248] = 2.01521643E-01;
    COFTD[249] = 5.62744089E-04;
    COFTD[250] = -3.08519239E-07;
    COFTD[251] = 5.22805986E-11;
    COFTD[252] = 1.22193921E-01;
    COFTD[253] = 6.90321128E-04;
    COFTD[254] = -3.64844875E-07;
    COFTD[255] = 6.03054876E-11;
}

#ifndef AMREX_USE_GPU
/* Initializes parameter database */
void CKINIT()
{

    // (0):  2.000000 O + M <=> O2 + M
    kiv[20] = {2,3};
    nuv[20] = {-2.0,1};
    // (0):  2.000000 O + M <=> O2 + M
    fwd_A[20]     = 1.2e+17;
    fwd_beta[20]  = -1;
    fwd_Ea[20]    = 0;
    prefactor_units[20]  = 1.0000000000000002e-12;
    activation_units[20] = 0.50321666580471969;
    phase_units[20]      = pow(10,-12.000000);
    is_PD[20] = 0;
    nTB[20] = 7;
    TB[20] = (amrex::Real *) malloc(7 * sizeof(amrex::Real));
    TBid[20] = (int *) malloc(7 * sizeof(int));
    TBid[20][0] = 0; TB[20][0] = 2.3999999999999999; // H2
    TBid[20][1] = 5; TB[20][1] = 15.4; // H2O
    TBid[20][2] = 13; TB[20][2] = 2; // CH4
    TBid[20][3] = 14; TB[20][3] = 1.75; // CO
    TBid[20][4] = 15; TB[20][4] = 3.6000000000000001; // CO2
    TBid[20][5] = 26; TB[20][5] = 3; // C2H6
    TBid[20][6] = 31; TB[20][6] = 0.82999999999999996; // AR

    // (1):  O + H + M <=> OH + M
    kiv[21] = {2,1,4};
    nuv[21] = {-1,-1,1};
    // (1):  O + H + M <=> OH + M
    fwd_A[21]     = 5e+17;
    fwd_beta[21]  = -1;
    fwd_Ea[21]    = 0;
    prefactor_units[21]  = 1.0000000000000002e-12;
    activation_units[21] = 0.50321666580471969;
    phase_units[21]      = pow(10,-12.000000);
    is_PD[21] = 0;
    nTB[21] = 7;
    TB[21] = (amrex::Real *) malloc(7 * sizeof(amrex::Real));
    TBid[21] = (int *) malloc(7 * sizeof(int));
    TBid[21][0] = 0; TB[21][0] = 2; // H2
    TBid[21][1] = 5; TB[21][1] = 6; // H2O
    TBid[21][2] = 13; TB[21][2] = 2; // CH4
    TBid[21][3] = 14; TB[21][3] = 1.5; // CO
    TBid[21][4] = 15; TB[21][4] = 2; // CO2
    TBid[21][5] = 26; TB[21][5] = 3; // C2H6
    TBid[21][6] = 31; TB[21][6] = 0.69999999999999996; // AR

    // (2):  O + H2 <=> H + OH
    kiv[27] = {2,0,1,4};
    nuv[27] = {-1,-1,1,1};
    // (2):  O + H2 <=> H + OH
    fwd_A[27]     = 50000;
    fwd_beta[27]  = 2.6699999999999999;
    fwd_Ea[27]    = 6290;
    prefactor_units[27]  = 1.0000000000000002e-06;
    activation_units[27] = 0.50321666580471969;
    phase_units[27]      = pow(10,-12.000000);
    is_PD[27] = 0;
    nTB[27] = 0;

    // (3):  O + HO2 <=> OH + O2
    kiv[28] = {2,6,4,3};
    nuv[28] = {-1,-1,1,1};
    // (3):  O + HO2 <=> OH + O2
    fwd_A[28]     = 20000000000000;
    fwd_beta[28]  = 0;
    fwd_Ea[28]    = 0;
    prefactor_units[28]  = 1.0000000000000002e-06;
    activation_units[28] = 0.50321666580471969;
    phase_units[28]      = pow(10,-12.000000);
    is_PD[28] = 0;
    nTB[28] = 0;

    // (4):  O + H2O2 <=> OH + HO2
    kiv[29] = {2,7,4,6};
    nuv[29] = {-1,-1,1,1};
    // (4):  O + H2O2 <=> OH + HO2
    fwd_A[29]     = 9630000;
    fwd_beta[29]  = 2;
    fwd_Ea[29]    = 4000;
    prefactor_units[29]  = 1.0000000000000002e-06;
    activation_units[29] = 0.50321666580471969;
    phase_units[29]      = pow(10,-12.000000);
    is_PD[29] = 0;
    nTB[29] = 0;

    // (5):  O + CH <=> H + CO
    kiv[30] = {2,9,1,14};
    nuv[30] = {-1,-1,1,1};
    // (5):  O + CH <=> H + CO
    fwd_A[30]     = 57000000000000;
    fwd_beta[30]  = 0;
    fwd_Ea[30]    = 0;
    prefactor_units[30]  = 1.0000000000000002e-06;
    activation_units[30] = 0.50321666580471969;
    phase_units[30]      = pow(10,-12.000000);
    is_PD[30] = 0;
    nTB[30] = 0;

    // (6):  O + CH2 <=> H + HCO
    kiv[31] = {2,10,1,16};
    nuv[31] = {-1,-1,1,1};
    // (6):  O + CH2 <=> H + HCO
    fwd_A[31]     = 80000000000000;
    fwd_beta[31]  = 0;
    fwd_Ea[31]    = 0;
    prefactor_units[31]  = 1.0000000000000002e-06;
    activation_units[31] = 0.50321666580471969;
    phase_units[31]      = pow(10,-12.000000);
    is_PD[31] = 0;
    nTB[31] = 0;

    // (7):  O + CH2(S) <=> H2 + CO
    kiv[32] = {2,11,0,14};
    nuv[32] = {-1,-1,1,1};
    // (7):  O + CH2(S) <=> H2 + CO
    fwd_A[32]     = 15000000000000;
    fwd_beta[32]  = 0;
    fwd_Ea[32]    = 0;
    prefactor_units[32]  = 1.0000000000000002e-06;
    activation_units[32] = 0.50321666580471969;
    phase_units[32]      = pow(10,-12.000000);
    is_PD[32] = 0;
    nTB[32] = 0;

    // (8):  O + CH2(S) <=> H + HCO
    kiv[33] = {2,11,1,16};
    nuv[33] = {-1,-1,1,1};
    // (8):  O + CH2(S) <=> H + HCO
    fwd_A[33]     = 15000000000000;
    fwd_beta[33]  = 0;
    fwd_Ea[33]    = 0;
    prefactor_units[33]  = 1.0000000000000002e-06;
    activation_units[33] = 0.50321666580471969;
    phase_units[33]      = pow(10,-12.000000);
    is_PD[33] = 0;
    nTB[33] = 0;

    // (9):  O + CH3 <=> H + CH2O
    kiv[34] = {2,12,1,17};
    nuv[34] = {-1,-1,1,1};
    // (9):  O + CH3 <=> H + CH2O
    fwd_A[34]     = 84300000000000;
    fwd_beta[34]  = 0;
    fwd_Ea[34]    = 0;
    prefactor_units[34]  = 1.0000000000000002e-06;
    activation_units[34] = 0.50321666580471969;
    phase_units[34]      = pow(10,-12.000000);
    is_PD[34] = 0;
    nTB[34] = 0;

    // (10):  O + CH4 <=> OH + CH3
    kiv[35] = {2,13,4,12};
    nuv[35] = {-1,-1,1,1};
    // (10):  O + CH4 <=> OH + CH3
    fwd_A[35]     = 1020000000;
    fwd_beta[35]  = 1.5;
    fwd_Ea[35]    = 8600;
    prefactor_units[35]  = 1.0000000000000002e-06;
    activation_units[35] = 0.50321666580471969;
    phase_units[35]      = pow(10,-12.000000);
    is_PD[35] = 0;
    nTB[35] = 0;

    // (11):  O + CO + M <=> CO2 + M
    kiv[22] = {2,14,15};
    nuv[22] = {-1,-1,1};
    // (11):  O + CO + M <=> CO2 + M
    fwd_A[22]     = 602000000000000;
    fwd_beta[22]  = 0;
    fwd_Ea[22]    = 3000;
    prefactor_units[22]  = 1.0000000000000002e-12;
    activation_units[22] = 0.50321666580471969;
    phase_units[22]      = pow(10,-12.000000);
    is_PD[22] = 0;
    nTB[22] = 8;
    TB[22] = (amrex::Real *) malloc(8 * sizeof(amrex::Real));
    TBid[22] = (int *) malloc(8 * sizeof(int));
    TBid[22][0] = 0; TB[22][0] = 2; // H2
    TBid[22][1] = 3; TB[22][1] = 6; // O2
    TBid[22][2] = 5; TB[22][2] = 6; // H2O
    TBid[22][3] = 13; TB[22][3] = 2; // CH4
    TBid[22][4] = 14; TB[22][4] = 1.5; // CO
    TBid[22][5] = 15; TB[22][5] = 3.5; // CO2
    TBid[22][6] = 26; TB[22][6] = 3; // C2H6
    TBid[22][7] = 31; TB[22][7] = 0.5; // AR

    // (12):  O + HCO <=> OH + CO
    kiv[36] = {2,16,4,14};
    nuv[36] = {-1,-1,1,1};
    // (12):  O + HCO <=> OH + CO
    fwd_A[36]     = 30000000000000;
    fwd_beta[36]  = 0;
    fwd_Ea[36]    = 0;
    prefactor_units[36]  = 1.0000000000000002e-06;
    activation_units[36] = 0.50321666580471969;
    phase_units[36]      = pow(10,-12.000000);
    is_PD[36] = 0;
    nTB[36] = 0;

    // (13):  O + HCO <=> H + CO2
    kiv[37] = {2,16,1,15};
    nuv[37] = {-1,-1,1,1};
    // (13):  O + HCO <=> H + CO2
    fwd_A[37]     = 30000000000000;
    fwd_beta[37]  = 0;
    fwd_Ea[37]    = 0;
    prefactor_units[37]  = 1.0000000000000002e-06;
    activation_units[37] = 0.50321666580471969;
    phase_units[37]      = pow(10,-12.000000);
    is_PD[37] = 0;
    nTB[37] = 0;

    // (14):  O + CH2O <=> OH + HCO
    kiv[38] = {2,17,4,16};
    nuv[38] = {-1,-1,1,1};
    // (14):  O + CH2O <=> OH + HCO
    fwd_A[38]     = 39000000000000;
    fwd_beta[38]  = 0;
    fwd_Ea[38]    = 3540;
    prefactor_units[38]  = 1.0000000000000002e-06;
    activation_units[38] = 0.50321666580471969;
    phase_units[38]      = pow(10,-12.000000);
    is_PD[38] = 0;
    nTB[38] = 0;

    // (15):  O + CH2OH <=> OH + CH2O
    kiv[39] = {2,18,4,17};
    nuv[39] = {-1,-1,1,1};
    // (15):  O + CH2OH <=> OH + CH2O
    fwd_A[39]     = 10000000000000;
    fwd_beta[39]  = 0;
    fwd_Ea[39]    = 0;
    prefactor_units[39]  = 1.0000000000000002e-06;
    activation_units[39] = 0.50321666580471969;
    phase_units[39]      = pow(10,-12.000000);
    is_PD[39] = 0;
    nTB[39] = 0;

    // (16):  O + CH3O <=> OH + CH2O
    kiv[40] = {2,19,4,17};
    nuv[40] = {-1,-1,1,1};
    // (16):  O + CH3O <=> OH + CH2O
    fwd_A[40]     = 10000000000000;
    fwd_beta[40]  = 0;
    fwd_Ea[40]    = 0;
    prefactor_units[40]  = 1.0000000000000002e-06;
    activation_units[40] = 0.50321666580471969;
    phase_units[40]      = pow(10,-12.000000);
    is_PD[40] = 0;
    nTB[40] = 0;

    // (17):  O + CH3OH <=> OH + CH2OH
    kiv[41] = {2,20,4,18};
    nuv[41] = {-1,-1,1,1};
    // (17):  O + CH3OH <=> OH + CH2OH
    fwd_A[41]     = 388000;
    fwd_beta[41]  = 2.5;
    fwd_Ea[41]    = 3100;
    prefactor_units[41]  = 1.0000000000000002e-06;
    activation_units[41] = 0.50321666580471969;
    phase_units[41]      = pow(10,-12.000000);
    is_PD[41] = 0;
    nTB[41] = 0;

    // (18):  O + CH3OH <=> OH + CH3O
    kiv[42] = {2,20,4,19};
    nuv[42] = {-1,-1,1,1};
    // (18):  O + CH3OH <=> OH + CH3O
    fwd_A[42]     = 130000;
    fwd_beta[42]  = 2.5;
    fwd_Ea[42]    = 5000;
    prefactor_units[42]  = 1.0000000000000002e-06;
    activation_units[42] = 0.50321666580471969;
    phase_units[42]      = pow(10,-12.000000);
    is_PD[42] = 0;
    nTB[42] = 0;

    // (19):  O + C2H <=> CH + CO
    kiv[43] = {2,21,9,14};
    nuv[43] = {-1,-1,1,1};
    // (19):  O + C2H <=> CH + CO
    fwd_A[43]     = 50000000000000;
    fwd_beta[43]  = 0;
    fwd_Ea[43]    = 0;
    prefactor_units[43]  = 1.0000000000000002e-06;
    activation_units[43] = 0.50321666580471969;
    phase_units[43]      = pow(10,-12.000000);
    is_PD[43] = 0;
    nTB[43] = 0;

    // (20):  O + C2H2 <=> H + HCCO
    kiv[44] = {2,22,1,27};
    nuv[44] = {-1,-1,1,1};
    // (20):  O + C2H2 <=> H + HCCO
    fwd_A[44]     = 10200000;
    fwd_beta[44]  = 2;
    fwd_Ea[44]    = 1900;
    prefactor_units[44]  = 1.0000000000000002e-06;
    activation_units[44] = 0.50321666580471969;
    phase_units[44]      = pow(10,-12.000000);
    is_PD[44] = 0;
    nTB[44] = 0;

    // (21):  O + C2H2 <=> OH + C2H
    kiv[45] = {2,22,4,21};
    nuv[45] = {-1,-1,1,1};
    // (21):  O + C2H2 <=> OH + C2H
    fwd_A[45]     = 4.6e+19;
    fwd_beta[45]  = -1.4099999999999999;
    fwd_Ea[45]    = 28950;
    prefactor_units[45]  = 1.0000000000000002e-06;
    activation_units[45] = 0.50321666580471969;
    phase_units[45]      = pow(10,-12.000000);
    is_PD[45] = 0;
    nTB[45] = 0;

    // (22):  O + C2H2 <=> CO + CH2
    kiv[46] = {2,22,14,10};
    nuv[46] = {-1,-1,1,1};
    // (22):  O + C2H2 <=> CO + CH2
    fwd_A[46]     = 10200000;
    fwd_beta[46]  = 2;
    fwd_Ea[46]    = 1900;
    prefactor_units[46]  = 1.0000000000000002e-06;
    activation_units[46] = 0.50321666580471969;
    phase_units[46]      = pow(10,-12.000000);
    is_PD[46] = 0;
    nTB[46] = 0;

    // (23):  O + C2H3 <=> H + CH2CO
    kiv[47] = {2,23,1,28};
    nuv[47] = {-1,-1,1,1};
    // (23):  O + C2H3 <=> H + CH2CO
    fwd_A[47]     = 30000000000000;
    fwd_beta[47]  = 0;
    fwd_Ea[47]    = 0;
    prefactor_units[47]  = 1.0000000000000002e-06;
    activation_units[47] = 0.50321666580471969;
    phase_units[47]      = pow(10,-12.000000);
    is_PD[47] = 0;
    nTB[47] = 0;

    // (24):  O + C2H4 <=> CH3 + HCO
    kiv[48] = {2,24,12,16};
    nuv[48] = {-1,-1,1,1};
    // (24):  O + C2H4 <=> CH3 + HCO
    fwd_A[48]     = 19200000;
    fwd_beta[48]  = 1.8300000000000001;
    fwd_Ea[48]    = 220;
    prefactor_units[48]  = 1.0000000000000002e-06;
    activation_units[48] = 0.50321666580471969;
    phase_units[48]      = pow(10,-12.000000);
    is_PD[48] = 0;
    nTB[48] = 0;

    // (25):  O + C2H5 <=> CH3 + CH2O
    kiv[49] = {2,25,12,17};
    nuv[49] = {-1,-1,1,1};
    // (25):  O + C2H5 <=> CH3 + CH2O
    fwd_A[49]     = 132000000000000;
    fwd_beta[49]  = 0;
    fwd_Ea[49]    = 0;
    prefactor_units[49]  = 1.0000000000000002e-06;
    activation_units[49] = 0.50321666580471969;
    phase_units[49]      = pow(10,-12.000000);
    is_PD[49] = 0;
    nTB[49] = 0;

    // (26):  O + C2H6 <=> OH + C2H5
    kiv[50] = {2,26,4,25};
    nuv[50] = {-1,-1,1,1};
    // (26):  O + C2H6 <=> OH + C2H5
    fwd_A[50]     = 89800000;
    fwd_beta[50]  = 1.9199999999999999;
    fwd_Ea[50]    = 5690;
    prefactor_units[50]  = 1.0000000000000002e-06;
    activation_units[50] = 0.50321666580471969;
    phase_units[50]      = pow(10,-12.000000);
    is_PD[50] = 0;
    nTB[50] = 0;

    // (27):  O + HCCO <=> H + 2.000000 CO
    kiv[51] = {2,27,1,14};
    nuv[51] = {-1,-1,1,2.0};
    // (27):  O + HCCO <=> H + 2.000000 CO
    fwd_A[51]     = 100000000000000;
    fwd_beta[51]  = 0;
    fwd_Ea[51]    = 0;
    prefactor_units[51]  = 1.0000000000000002e-06;
    activation_units[51] = 0.50321666580471969;
    phase_units[51]      = pow(10,-12.000000);
    is_PD[51] = 0;
    nTB[51] = 0;

    // (28):  O + CH2CO <=> OH + HCCO
    kiv[52] = {2,28,4,27};
    nuv[52] = {-1,-1,1,1};
    // (28):  O + CH2CO <=> OH + HCCO
    fwd_A[52]     = 10000000000000;
    fwd_beta[52]  = 0;
    fwd_Ea[52]    = 8000;
    prefactor_units[52]  = 1.0000000000000002e-06;
    activation_units[52] = 0.50321666580471969;
    phase_units[52]      = pow(10,-12.000000);
    is_PD[52] = 0;
    nTB[52] = 0;

    // (29):  O + CH2CO <=> CH2 + CO2
    kiv[53] = {2,28,10,15};
    nuv[53] = {-1,-1,1,1};
    // (29):  O + CH2CO <=> CH2 + CO2
    fwd_A[53]     = 1750000000000;
    fwd_beta[53]  = 0;
    fwd_Ea[53]    = 1350;
    prefactor_units[53]  = 1.0000000000000002e-06;
    activation_units[53] = 0.50321666580471969;
    phase_units[53]      = pow(10,-12.000000);
    is_PD[53] = 0;
    nTB[53] = 0;

    // (30):  O2 + CO <=> O + CO2
    kiv[54] = {3,14,2,15};
    nuv[54] = {-1,-1,1,1};
    // (30):  O2 + CO <=> O + CO2
    fwd_A[54]     = 2500000000000;
    fwd_beta[54]  = 0;
    fwd_Ea[54]    = 47800;
    prefactor_units[54]  = 1.0000000000000002e-06;
    activation_units[54] = 0.50321666580471969;
    phase_units[54]      = pow(10,-12.000000);
    is_PD[54] = 0;
    nTB[54] = 0;

    // (31):  O2 + CH2O <=> HO2 + HCO
    kiv[55] = {3,17,6,16};
    nuv[55] = {-1,-1,1,1};
    // (31):  O2 + CH2O <=> HO2 + HCO
    fwd_A[55]     = 100000000000000;
    fwd_beta[55]  = 0;
    fwd_Ea[55]    = 40000;
    prefactor_units[55]  = 1.0000000000000002e-06;
    activation_units[55] = 0.50321666580471969;
    phase_units[55]      = pow(10,-12.000000);
    is_PD[55] = 0;
    nTB[55] = 0;

    // (32):  H + O2 + M <=> HO2 + M
    kiv[23] = {1,3,6};
    nuv[23] = {-1,-1,1};
    // (32):  H + O2 + M <=> HO2 + M
    fwd_A[23]     = 2.8e+18;
    fwd_beta[23]  = -0.85999999999999999;
    fwd_Ea[23]    = 0;
    prefactor_units[23]  = 1.0000000000000002e-12;
    activation_units[23] = 0.50321666580471969;
    phase_units[23]      = pow(10,-12.000000);
    is_PD[23] = 0;
    nTB[23] = 7;
    TB[23] = (amrex::Real *) malloc(7 * sizeof(amrex::Real));
    TBid[23] = (int *) malloc(7 * sizeof(int));
    TBid[23][0] = 3; TB[23][0] = 0; // O2
    TBid[23][1] = 5; TB[23][1] = 0; // H2O
    TBid[23][2] = 14; TB[23][2] = 0.75; // CO
    TBid[23][3] = 15; TB[23][3] = 1.5; // CO2
    TBid[23][4] = 26; TB[23][4] = 1.5; // C2H6
    TBid[23][5] = 30; TB[23][5] = 0; // N2
    TBid[23][6] = 31; TB[23][6] = 0; // AR

    // (33):  H + 2.000000 O2 <=> HO2 + O2
    kiv[56] = {1,3,6,3};
    nuv[56] = {-1,-2.0,1,1};
    // (33):  H + 2.000000 O2 <=> HO2 + O2
    fwd_A[56]     = 3e+20;
    fwd_beta[56]  = -1.72;
    fwd_Ea[56]    = 0;
    prefactor_units[56]  = 1.0000000000000002e-12;
    activation_units[56] = 0.50321666580471969;
    phase_units[56]      = pow(10,-18.000000);
    is_PD[56] = 0;
    nTB[56] = 0;

    // (34):  H + O2 + H2O <=> HO2 + H2O
    kiv[57] = {1,3,5,6,5};
    nuv[57] = {-1,-1,-1,1,1};
    // (34):  H + O2 + H2O <=> HO2 + H2O
    fwd_A[57]     = 9.38e+18;
    fwd_beta[57]  = -0.76000000000000001;
    fwd_Ea[57]    = 0;
    prefactor_units[57]  = 1.0000000000000002e-12;
    activation_units[57] = 0.50321666580471969;
    phase_units[57]      = pow(10,-18.000000);
    is_PD[57] = 0;
    nTB[57] = 0;

    // (35):  H + O2 + N2 <=> HO2 + N2
    kiv[58] = {1,3,30,6,30};
    nuv[58] = {-1,-1,-1,1,1};
    // (35):  H + O2 + N2 <=> HO2 + N2
    fwd_A[58]     = 3.75e+20;
    fwd_beta[58]  = -1.72;
    fwd_Ea[58]    = 0;
    prefactor_units[58]  = 1.0000000000000002e-12;
    activation_units[58] = 0.50321666580471969;
    phase_units[58]      = pow(10,-18.000000);
    is_PD[58] = 0;
    nTB[58] = 0;

    // (36):  H + O2 + AR <=> HO2 + AR
    kiv[59] = {1,3,31,6,31};
    nuv[59] = {-1,-1,-1,1,1};
    // (36):  H + O2 + AR <=> HO2 + AR
    fwd_A[59]     = 7e+17;
    fwd_beta[59]  = -0.80000000000000004;
    fwd_Ea[59]    = 0;
    prefactor_units[59]  = 1.0000000000000002e-12;
    activation_units[59] = 0.50321666580471969;
    phase_units[59]      = pow(10,-18.000000);
    is_PD[59] = 0;
    nTB[59] = 0;

    // (37):  H + O2 <=> O + OH
    kiv[60] = {1,3,2,4};
    nuv[60] = {-1,-1,1,1};
    // (37):  H + O2 <=> O + OH
    fwd_A[60]     = 83000000000000;
    fwd_beta[60]  = 0;
    fwd_Ea[60]    = 14413;
    prefactor_units[60]  = 1.0000000000000002e-06;
    activation_units[60] = 0.50321666580471969;
    phase_units[60]      = pow(10,-12.000000);
    is_PD[60] = 0;
    nTB[60] = 0;

    // (38):  2.000000 H + M <=> H2 + M
    kiv[24] = {1,0};
    nuv[24] = {-2.0,1};
    // (38):  2.000000 H + M <=> H2 + M
    fwd_A[24]     = 1e+18;
    fwd_beta[24]  = -1;
    fwd_Ea[24]    = 0;
    prefactor_units[24]  = 1.0000000000000002e-12;
    activation_units[24] = 0.50321666580471969;
    phase_units[24]      = pow(10,-12.000000);
    is_PD[24] = 0;
    nTB[24] = 6;
    TB[24] = (amrex::Real *) malloc(6 * sizeof(amrex::Real));
    TBid[24] = (int *) malloc(6 * sizeof(int));
    TBid[24][0] = 0; TB[24][0] = 0; // H2
    TBid[24][1] = 5; TB[24][1] = 0; // H2O
    TBid[24][2] = 13; TB[24][2] = 2; // CH4
    TBid[24][3] = 15; TB[24][3] = 0; // CO2
    TBid[24][4] = 26; TB[24][4] = 3; // C2H6
    TBid[24][5] = 31; TB[24][5] = 0.63; // AR

    // (39):  2.000000 H + H2 <=> 2.000000 H2
    kiv[61] = {1,0,0};
    nuv[61] = {-2.0,-1,2.0};
    // (39):  2.000000 H + H2 <=> 2.000000 H2
    fwd_A[61]     = 90000000000000000;
    fwd_beta[61]  = -0.59999999999999998;
    fwd_Ea[61]    = 0;
    prefactor_units[61]  = 1.0000000000000002e-12;
    activation_units[61] = 0.50321666580471969;
    phase_units[61]      = pow(10,-18.000000);
    is_PD[61] = 0;
    nTB[61] = 0;

    // (40):  2.000000 H + H2O <=> H2 + H2O
    kiv[62] = {1,5,0,5};
    nuv[62] = {-2.0,-1,1,1};
    // (40):  2.000000 H + H2O <=> H2 + H2O
    fwd_A[62]     = 6e+19;
    fwd_beta[62]  = -1.25;
    fwd_Ea[62]    = 0;
    prefactor_units[62]  = 1.0000000000000002e-12;
    activation_units[62] = 0.50321666580471969;
    phase_units[62]      = pow(10,-18.000000);
    is_PD[62] = 0;
    nTB[62] = 0;

    // (41):  2.000000 H + CO2 <=> H2 + CO2
    kiv[63] = {1,15,0,15};
    nuv[63] = {-2.0,-1,1,1};
    // (41):  2.000000 H + CO2 <=> H2 + CO2
    fwd_A[63]     = 5.5e+20;
    fwd_beta[63]  = -2;
    fwd_Ea[63]    = 0;
    prefactor_units[63]  = 1.0000000000000002e-12;
    activation_units[63] = 0.50321666580471969;
    phase_units[63]      = pow(10,-18.000000);
    is_PD[63] = 0;
    nTB[63] = 0;

    // (42):  H + OH + M <=> H2O + M
    kiv[25] = {1,4,5};
    nuv[25] = {-1,-1,1};
    // (42):  H + OH + M <=> H2O + M
    fwd_A[25]     = 2.2e+22;
    fwd_beta[25]  = -2;
    fwd_Ea[25]    = 0;
    prefactor_units[25]  = 1.0000000000000002e-12;
    activation_units[25] = 0.50321666580471969;
    phase_units[25]      = pow(10,-12.000000);
    is_PD[25] = 0;
    nTB[25] = 5;
    TB[25] = (amrex::Real *) malloc(5 * sizeof(amrex::Real));
    TBid[25] = (int *) malloc(5 * sizeof(int));
    TBid[25][0] = 0; TB[25][0] = 0.72999999999999998; // H2
    TBid[25][1] = 5; TB[25][1] = 3.6499999999999999; // H2O
    TBid[25][2] = 13; TB[25][2] = 2; // CH4
    TBid[25][3] = 26; TB[25][3] = 3; // C2H6
    TBid[25][4] = 31; TB[25][4] = 0.38; // AR

    // (43):  H + HO2 <=> O + H2O
    kiv[64] = {1,6,2,5};
    nuv[64] = {-1,-1,1,1};
    // (43):  H + HO2 <=> O + H2O
    fwd_A[64]     = 3970000000000;
    fwd_beta[64]  = 0;
    fwd_Ea[64]    = 671;
    prefactor_units[64]  = 1.0000000000000002e-06;
    activation_units[64] = 0.50321666580471969;
    phase_units[64]      = pow(10,-12.000000);
    is_PD[64] = 0;
    nTB[64] = 0;

    // (44):  H + HO2 <=> O2 + H2
    kiv[65] = {1,6,3,0};
    nuv[65] = {-1,-1,1,1};
    // (44):  H + HO2 <=> O2 + H2
    fwd_A[65]     = 28000000000000;
    fwd_beta[65]  = 0;
    fwd_Ea[65]    = 1068;
    prefactor_units[65]  = 1.0000000000000002e-06;
    activation_units[65] = 0.50321666580471969;
    phase_units[65]      = pow(10,-12.000000);
    is_PD[65] = 0;
    nTB[65] = 0;

    // (45):  H + HO2 <=> 2.000000 OH
    kiv[66] = {1,6,4};
    nuv[66] = {-1,-1,2.0};
    // (45):  H + HO2 <=> 2.000000 OH
    fwd_A[66]     = 134000000000000;
    fwd_beta[66]  = 0;
    fwd_Ea[66]    = 635;
    prefactor_units[66]  = 1.0000000000000002e-06;
    activation_units[66] = 0.50321666580471969;
    phase_units[66]      = pow(10,-12.000000);
    is_PD[66] = 0;
    nTB[66] = 0;

    // (46):  H + H2O2 <=> HO2 + H2
    kiv[67] = {1,7,6,0};
    nuv[67] = {-1,-1,1,1};
    // (46):  H + H2O2 <=> HO2 + H2
    fwd_A[67]     = 12100000;
    fwd_beta[67]  = 2;
    fwd_Ea[67]    = 5200;
    prefactor_units[67]  = 1.0000000000000002e-06;
    activation_units[67] = 0.50321666580471969;
    phase_units[67]      = pow(10,-12.000000);
    is_PD[67] = 0;
    nTB[67] = 0;

    // (47):  H + H2O2 <=> OH + H2O
    kiv[68] = {1,7,4,5};
    nuv[68] = {-1,-1,1,1};
    // (47):  H + H2O2 <=> OH + H2O
    fwd_A[68]     = 10000000000000;
    fwd_beta[68]  = 0;
    fwd_Ea[68]    = 3600;
    prefactor_units[68]  = 1.0000000000000002e-06;
    activation_units[68] = 0.50321666580471969;
    phase_units[68]      = pow(10,-12.000000);
    is_PD[68] = 0;
    nTB[68] = 0;

    // (48):  H + CH <=> C + H2
    kiv[69] = {1,9,8,0};
    nuv[69] = {-1,-1,1,1};
    // (48):  H + CH <=> C + H2
    fwd_A[69]     = 110000000000000;
    fwd_beta[69]  = 0;
    fwd_Ea[69]    = 0;
    prefactor_units[69]  = 1.0000000000000002e-06;
    activation_units[69] = 0.50321666580471969;
    phase_units[69]      = pow(10,-12.000000);
    is_PD[69] = 0;
    nTB[69] = 0;

    // (49):  H + CH2 (+M) <=> CH3 (+M)
    kiv[0] = {1,10,12};
    nuv[0] = {-1,-1,1};
    // (49):  H + CH2 (+M) <=> CH3 (+M)
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
    TBid[0][2] = 13; TB[0][2] = 2; // CH4
    TBid[0][3] = 14; TB[0][3] = 1.5; // CO
    TBid[0][4] = 15; TB[0][4] = 2; // CO2
    TBid[0][5] = 26; TB[0][5] = 3; // C2H6
    TBid[0][6] = 31; TB[0][6] = 0.69999999999999996; // AR

    // (50):  H + CH2(S) <=> CH + H2
    kiv[70] = {1,11,9,0};
    nuv[70] = {-1,-1,1,1};
    // (50):  H + CH2(S) <=> CH + H2
    fwd_A[70]     = 30000000000000;
    fwd_beta[70]  = 0;
    fwd_Ea[70]    = 0;
    prefactor_units[70]  = 1.0000000000000002e-06;
    activation_units[70] = 0.50321666580471969;
    phase_units[70]      = pow(10,-12.000000);
    is_PD[70] = 0;
    nTB[70] = 0;

    // (51):  H + CH3 (+M) <=> CH4 (+M)
    kiv[1] = {1,12,13};
    nuv[1] = {-1,-1,1};
    // (51):  H + CH3 (+M) <=> CH4 (+M)
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
    TBid[1][2] = 13; TB[1][2] = 2; // CH4
    TBid[1][3] = 14; TB[1][3] = 1.5; // CO
    TBid[1][4] = 15; TB[1][4] = 2; // CO2
    TBid[1][5] = 26; TB[1][5] = 3; // C2H6
    TBid[1][6] = 31; TB[1][6] = 0.69999999999999996; // AR

    // (52):  H + CH4 <=> CH3 + H2
    kiv[71] = {1,13,12,0};
    nuv[71] = {-1,-1,1,1};
    // (52):  H + CH4 <=> CH3 + H2
    fwd_A[71]     = 660000000;
    fwd_beta[71]  = 1.6200000000000001;
    fwd_Ea[71]    = 10840;
    prefactor_units[71]  = 1.0000000000000002e-06;
    activation_units[71] = 0.50321666580471969;
    phase_units[71]      = pow(10,-12.000000);
    is_PD[71] = 0;
    nTB[71] = 0;

    // (53):  H + HCO (+M) <=> CH2O (+M)
    kiv[2] = {1,16,17};
    nuv[2] = {-1,-1,1};
    // (53):  H + HCO (+M) <=> CH2O (+M)
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
    TBid[2][2] = 13; TB[2][2] = 2; // CH4
    TBid[2][3] = 14; TB[2][3] = 1.5; // CO
    TBid[2][4] = 15; TB[2][4] = 2; // CO2
    TBid[2][5] = 26; TB[2][5] = 3; // C2H6
    TBid[2][6] = 31; TB[2][6] = 0.69999999999999996; // AR

    // (54):  H + HCO <=> H2 + CO
    kiv[72] = {1,16,0,14};
    nuv[72] = {-1,-1,1,1};
    // (54):  H + HCO <=> H2 + CO
    fwd_A[72]     = 73400000000000;
    fwd_beta[72]  = 0;
    fwd_Ea[72]    = 0;
    prefactor_units[72]  = 1.0000000000000002e-06;
    activation_units[72] = 0.50321666580471969;
    phase_units[72]      = pow(10,-12.000000);
    is_PD[72] = 0;
    nTB[72] = 0;

    // (55):  H + CH2O (+M) <=> CH2OH (+M)
    kiv[3] = {1,17,18};
    nuv[3] = {-1,-1,1};
    // (55):  H + CH2O (+M) <=> CH2OH (+M)
    fwd_A[3]     = 540000000000;
    fwd_beta[3]  = 0.45400000000000001;
    fwd_Ea[3]    = 3600;
    low_A[3]     = 1.27e+32;
    low_beta[3]  = -4.8200000000000003;
    low_Ea[3]    = 6530;
    troe_a[3]    = 0.71870000000000001;
    troe_Tsss[3] = 103;
    troe_Ts[3]   = 1291;
    troe_Tss[3]  = 4160;
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
    TBid[3][2] = 13; TB[3][2] = 2; // CH4
    TBid[3][3] = 14; TB[3][3] = 1.5; // CO
    TBid[3][4] = 15; TB[3][4] = 2; // CO2
    TBid[3][5] = 26; TB[3][5] = 3; // C2H6

    // (56):  H + CH2O (+M) <=> CH3O (+M)
    kiv[4] = {1,17,19};
    nuv[4] = {-1,-1,1};
    // (56):  H + CH2O (+M) <=> CH3O (+M)
    fwd_A[4]     = 540000000000;
    fwd_beta[4]  = 0.45400000000000001;
    fwd_Ea[4]    = 2600;
    low_A[4]     = 2.2e+30;
    low_beta[4]  = -4.7999999999999998;
    low_Ea[4]    = 5560;
    troe_a[4]    = 0.75800000000000001;
    troe_Tsss[4] = 94;
    troe_Ts[4]   = 1555;
    troe_Tss[4]  = 4200;
    troe_len[4]  = 4;
    prefactor_units[4]  = 1.0000000000000002e-06;
    activation_units[4] = 0.50321666580471969;
    phase_units[4]      = pow(10,-12.000000);
    is_PD[4] = 1;
    nTB[4] = 6;
    TB[4] = (amrex::Real *) malloc(6 * sizeof(amrex::Real));
    TBid[4] = (int *) malloc(6 * sizeof(int));
    TBid[4][0] = 0; TB[4][0] = 2; // H2
    TBid[4][1] = 5; TB[4][1] = 6; // H2O
    TBid[4][2] = 13; TB[4][2] = 2; // CH4
    TBid[4][3] = 14; TB[4][3] = 1.5; // CO
    TBid[4][4] = 15; TB[4][4] = 2; // CO2
    TBid[4][5] = 26; TB[4][5] = 3; // C2H6

    // (57):  H + CH2O <=> HCO + H2
    kiv[73] = {1,17,16,0};
    nuv[73] = {-1,-1,1,1};
    // (57):  H + CH2O <=> HCO + H2
    fwd_A[73]     = 23000000000;
    fwd_beta[73]  = 1.05;
    fwd_Ea[73]    = 3275;
    prefactor_units[73]  = 1.0000000000000002e-06;
    activation_units[73] = 0.50321666580471969;
    phase_units[73]      = pow(10,-12.000000);
    is_PD[73] = 0;
    nTB[73] = 0;

    // (58):  H + CH2OH (+M) <=> CH3OH (+M)
    kiv[5] = {1,18,20};
    nuv[5] = {-1,-1,1};
    // (58):  H + CH2OH (+M) <=> CH3OH (+M)
    fwd_A[5]     = 18000000000000;
    fwd_beta[5]  = 0;
    fwd_Ea[5]    = 0;
    low_A[5]     = 2.9999999999999999e+31;
    low_beta[5]  = -4.7999999999999998;
    low_Ea[5]    = 3300;
    troe_a[5]    = 0.76790000000000003;
    troe_Tsss[5] = 338;
    troe_Ts[5]   = 1812;
    troe_Tss[5]  = 5081;
    troe_len[5]  = 4;
    prefactor_units[5]  = 1.0000000000000002e-06;
    activation_units[5] = 0.50321666580471969;
    phase_units[5]      = pow(10,-12.000000);
    is_PD[5] = 1;
    nTB[5] = 6;
    TB[5] = (amrex::Real *) malloc(6 * sizeof(amrex::Real));
    TBid[5] = (int *) malloc(6 * sizeof(int));
    TBid[5][0] = 0; TB[5][0] = 2; // H2
    TBid[5][1] = 5; TB[5][1] = 6; // H2O
    TBid[5][2] = 13; TB[5][2] = 2; // CH4
    TBid[5][3] = 14; TB[5][3] = 1.5; // CO
    TBid[5][4] = 15; TB[5][4] = 2; // CO2
    TBid[5][5] = 26; TB[5][5] = 3; // C2H6

    // (59):  H + CH2OH <=> H2 + CH2O
    kiv[74] = {1,18,0,17};
    nuv[74] = {-1,-1,1,1};
    // (59):  H + CH2OH <=> H2 + CH2O
    fwd_A[74]     = 20000000000000;
    fwd_beta[74]  = 0;
    fwd_Ea[74]    = 0;
    prefactor_units[74]  = 1.0000000000000002e-06;
    activation_units[74] = 0.50321666580471969;
    phase_units[74]      = pow(10,-12.000000);
    is_PD[74] = 0;
    nTB[74] = 0;

    // (60):  H + CH2OH <=> OH + CH3
    kiv[75] = {1,18,4,12};
    nuv[75] = {-1,-1,1,1};
    // (60):  H + CH2OH <=> OH + CH3
    fwd_A[75]     = 12000000000000;
    fwd_beta[75]  = 0;
    fwd_Ea[75]    = 0;
    prefactor_units[75]  = 1.0000000000000002e-06;
    activation_units[75] = 0.50321666580471969;
    phase_units[75]      = pow(10,-12.000000);
    is_PD[75] = 0;
    nTB[75] = 0;

    // (61):  H + CH2OH <=> CH2(S) + H2O
    kiv[76] = {1,18,11,5};
    nuv[76] = {-1,-1,1,1};
    // (61):  H + CH2OH <=> CH2(S) + H2O
    fwd_A[76]     = 6000000000000;
    fwd_beta[76]  = 0;
    fwd_Ea[76]    = 0;
    prefactor_units[76]  = 1.0000000000000002e-06;
    activation_units[76] = 0.50321666580471969;
    phase_units[76]      = pow(10,-12.000000);
    is_PD[76] = 0;
    nTB[76] = 0;

    // (62):  H + CH3O (+M) <=> CH3OH (+M)
    kiv[6] = {1,19,20};
    nuv[6] = {-1,-1,1};
    // (62):  H + CH3O (+M) <=> CH3OH (+M)
    fwd_A[6]     = 50000000000000;
    fwd_beta[6]  = 0;
    fwd_Ea[6]    = 0;
    low_A[6]     = 8.5999999999999995e+28;
    low_beta[6]  = -4;
    low_Ea[6]    = 3025;
    troe_a[6]    = 0.89019999999999999;
    troe_Tsss[6] = 144;
    troe_Ts[6]   = 2838;
    troe_Tss[6]  = 45569;
    troe_len[6]  = 4;
    prefactor_units[6]  = 1.0000000000000002e-06;
    activation_units[6] = 0.50321666580471969;
    phase_units[6]      = pow(10,-12.000000);
    is_PD[6] = 1;
    nTB[6] = 6;
    TB[6] = (amrex::Real *) malloc(6 * sizeof(amrex::Real));
    TBid[6] = (int *) malloc(6 * sizeof(int));
    TBid[6][0] = 0; TB[6][0] = 2; // H2
    TBid[6][1] = 5; TB[6][1] = 6; // H2O
    TBid[6][2] = 13; TB[6][2] = 2; // CH4
    TBid[6][3] = 14; TB[6][3] = 1.5; // CO
    TBid[6][4] = 15; TB[6][4] = 2; // CO2
    TBid[6][5] = 26; TB[6][5] = 3; // C2H6

    // (63):  H + CH3O <=> H + CH2OH
    kiv[77] = {1,19,1,18};
    nuv[77] = {-1,-1,1,1};
    // (63):  H + CH3O <=> H + CH2OH
    fwd_A[77]     = 3400000;
    fwd_beta[77]  = 1.6000000000000001;
    fwd_Ea[77]    = 0;
    prefactor_units[77]  = 1.0000000000000002e-06;
    activation_units[77] = 0.50321666580471969;
    phase_units[77]      = pow(10,-12.000000);
    is_PD[77] = 0;
    nTB[77] = 0;

    // (64):  H + CH3O <=> H2 + CH2O
    kiv[78] = {1,19,0,17};
    nuv[78] = {-1,-1,1,1};
    // (64):  H + CH3O <=> H2 + CH2O
    fwd_A[78]     = 20000000000000;
    fwd_beta[78]  = 0;
    fwd_Ea[78]    = 0;
    prefactor_units[78]  = 1.0000000000000002e-06;
    activation_units[78] = 0.50321666580471969;
    phase_units[78]      = pow(10,-12.000000);
    is_PD[78] = 0;
    nTB[78] = 0;

    // (65):  H + CH3O <=> OH + CH3
    kiv[79] = {1,19,4,12};
    nuv[79] = {-1,-1,1,1};
    // (65):  H + CH3O <=> OH + CH3
    fwd_A[79]     = 32000000000000;
    fwd_beta[79]  = 0;
    fwd_Ea[79]    = 0;
    prefactor_units[79]  = 1.0000000000000002e-06;
    activation_units[79] = 0.50321666580471969;
    phase_units[79]      = pow(10,-12.000000);
    is_PD[79] = 0;
    nTB[79] = 0;

    // (66):  H + CH3O <=> CH2(S) + H2O
    kiv[80] = {1,19,11,5};
    nuv[80] = {-1,-1,1,1};
    // (66):  H + CH3O <=> CH2(S) + H2O
    fwd_A[80]     = 16000000000000;
    fwd_beta[80]  = 0;
    fwd_Ea[80]    = 0;
    prefactor_units[80]  = 1.0000000000000002e-06;
    activation_units[80] = 0.50321666580471969;
    phase_units[80]      = pow(10,-12.000000);
    is_PD[80] = 0;
    nTB[80] = 0;

    // (67):  H + CH3OH <=> CH2OH + H2
    kiv[81] = {1,20,18,0};
    nuv[81] = {-1,-1,1,1};
    // (67):  H + CH3OH <=> CH2OH + H2
    fwd_A[81]     = 17000000;
    fwd_beta[81]  = 2.1000000000000001;
    fwd_Ea[81]    = 4870;
    prefactor_units[81]  = 1.0000000000000002e-06;
    activation_units[81] = 0.50321666580471969;
    phase_units[81]      = pow(10,-12.000000);
    is_PD[81] = 0;
    nTB[81] = 0;

    // (68):  H + CH3OH <=> CH3O + H2
    kiv[82] = {1,20,19,0};
    nuv[82] = {-1,-1,1,1};
    // (68):  H + CH3OH <=> CH3O + H2
    fwd_A[82]     = 4200000;
    fwd_beta[82]  = 2.1000000000000001;
    fwd_Ea[82]    = 4870;
    prefactor_units[82]  = 1.0000000000000002e-06;
    activation_units[82] = 0.50321666580471969;
    phase_units[82]      = pow(10,-12.000000);
    is_PD[82] = 0;
    nTB[82] = 0;

    // (69):  H + C2H (+M) <=> C2H2 (+M)
    kiv[7] = {1,21,22};
    nuv[7] = {-1,-1,1};
    // (69):  H + C2H (+M) <=> C2H2 (+M)
    fwd_A[7]     = 1e+17;
    fwd_beta[7]  = -1;
    fwd_Ea[7]    = 0;
    low_A[7]     = 3.7500000000000002e+33;
    low_beta[7]  = -4.7999999999999998;
    low_Ea[7]    = 1900;
    troe_a[7]    = 0.64639999999999997;
    troe_Tsss[7] = 132;
    troe_Ts[7]   = 1315;
    troe_Tss[7]  = 5566;
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
    TBid[7][2] = 13; TB[7][2] = 2; // CH4
    TBid[7][3] = 14; TB[7][3] = 1.5; // CO
    TBid[7][4] = 15; TB[7][4] = 2; // CO2
    TBid[7][5] = 26; TB[7][5] = 3; // C2H6
    TBid[7][6] = 31; TB[7][6] = 0.69999999999999996; // AR

    // (70):  H + C2H2 (+M) <=> C2H3 (+M)
    kiv[8] = {1,22,23};
    nuv[8] = {-1,-1,1};
    // (70):  H + C2H2 (+M) <=> C2H3 (+M)
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
    nTB[8] = 7;
    TB[8] = (amrex::Real *) malloc(7 * sizeof(amrex::Real));
    TBid[8] = (int *) malloc(7 * sizeof(int));
    TBid[8][0] = 0; TB[8][0] = 2; // H2
    TBid[8][1] = 5; TB[8][1] = 6; // H2O
    TBid[8][2] = 13; TB[8][2] = 2; // CH4
    TBid[8][3] = 14; TB[8][3] = 1.5; // CO
    TBid[8][4] = 15; TB[8][4] = 2; // CO2
    TBid[8][5] = 26; TB[8][5] = 3; // C2H6
    TBid[8][6] = 31; TB[8][6] = 0.69999999999999996; // AR

    // (71):  H + C2H3 (+M) <=> C2H4 (+M)
    kiv[9] = {1,23,24};
    nuv[9] = {-1,-1,1};
    // (71):  H + C2H3 (+M) <=> C2H4 (+M)
    fwd_A[9]     = 6080000000000;
    fwd_beta[9]  = 0.27000000000000002;
    fwd_Ea[9]    = 280;
    low_A[9]     = 1.3999999999999999e+30;
    low_beta[9]  = -3.8599999999999999;
    low_Ea[9]    = 3320;
    troe_a[9]    = 0.78200000000000003;
    troe_Tsss[9] = 207.5;
    troe_Ts[9]   = 2663;
    troe_Tss[9]  = 6095;
    troe_len[9]  = 4;
    prefactor_units[9]  = 1.0000000000000002e-06;
    activation_units[9] = 0.50321666580471969;
    phase_units[9]      = pow(10,-12.000000);
    is_PD[9] = 1;
    nTB[9] = 7;
    TB[9] = (amrex::Real *) malloc(7 * sizeof(amrex::Real));
    TBid[9] = (int *) malloc(7 * sizeof(int));
    TBid[9][0] = 0; TB[9][0] = 2; // H2
    TBid[9][1] = 5; TB[9][1] = 6; // H2O
    TBid[9][2] = 13; TB[9][2] = 2; // CH4
    TBid[9][3] = 14; TB[9][3] = 1.5; // CO
    TBid[9][4] = 15; TB[9][4] = 2; // CO2
    TBid[9][5] = 26; TB[9][5] = 3; // C2H6
    TBid[9][6] = 31; TB[9][6] = 0.69999999999999996; // AR

    // (72):  H + C2H3 <=> H2 + C2H2
    kiv[83] = {1,23,0,22};
    nuv[83] = {-1,-1,1,1};
    // (72):  H + C2H3 <=> H2 + C2H2
    fwd_A[83]     = 30000000000000;
    fwd_beta[83]  = 0;
    fwd_Ea[83]    = 0;
    prefactor_units[83]  = 1.0000000000000002e-06;
    activation_units[83] = 0.50321666580471969;
    phase_units[83]      = pow(10,-12.000000);
    is_PD[83] = 0;
    nTB[83] = 0;

    // (73):  H + C2H4 (+M) <=> C2H5 (+M)
    kiv[10] = {1,24,25};
    nuv[10] = {-1,-1,1};
    // (73):  H + C2H4 (+M) <=> C2H5 (+M)
    fwd_A[10]     = 1080000000000;
    fwd_beta[10]  = 0.45400000000000001;
    fwd_Ea[10]    = 1820;
    low_A[10]     = 1.1999999999999999e+42;
    low_beta[10]  = -7.6200000000000001;
    low_Ea[10]    = 6970;
    troe_a[10]    = 0.97529999999999994;
    troe_Tsss[10] = 210;
    troe_Ts[10]   = 984;
    troe_Tss[10]  = 4374;
    troe_len[10]  = 4;
    prefactor_units[10]  = 1.0000000000000002e-06;
    activation_units[10] = 0.50321666580471969;
    phase_units[10]      = pow(10,-12.000000);
    is_PD[10] = 1;
    nTB[10] = 7;
    TB[10] = (amrex::Real *) malloc(7 * sizeof(amrex::Real));
    TBid[10] = (int *) malloc(7 * sizeof(int));
    TBid[10][0] = 0; TB[10][0] = 2; // H2
    TBid[10][1] = 5; TB[10][1] = 6; // H2O
    TBid[10][2] = 13; TB[10][2] = 2; // CH4
    TBid[10][3] = 14; TB[10][3] = 1.5; // CO
    TBid[10][4] = 15; TB[10][4] = 2; // CO2
    TBid[10][5] = 26; TB[10][5] = 3; // C2H6
    TBid[10][6] = 31; TB[10][6] = 0.69999999999999996; // AR

    // (74):  H + C2H4 <=> C2H3 + H2
    kiv[84] = {1,24,23,0};
    nuv[84] = {-1,-1,1,1};
    // (74):  H + C2H4 <=> C2H3 + H2
    fwd_A[84]     = 1325000;
    fwd_beta[84]  = 2.5299999999999998;
    fwd_Ea[84]    = 12240;
    prefactor_units[84]  = 1.0000000000000002e-06;
    activation_units[84] = 0.50321666580471969;
    phase_units[84]      = pow(10,-12.000000);
    is_PD[84] = 0;
    nTB[84] = 0;

    // (75):  H + C2H5 (+M) <=> C2H6 (+M)
    kiv[11] = {1,25,26};
    nuv[11] = {-1,-1,1};
    // (75):  H + C2H5 (+M) <=> C2H6 (+M)
    fwd_A[11]     = 5.21e+17;
    fwd_beta[11]  = -0.98999999999999999;
    fwd_Ea[11]    = 1580;
    low_A[11]     = 1.9900000000000001e+41;
    low_beta[11]  = -7.0800000000000001;
    low_Ea[11]    = 6685;
    troe_a[11]    = 0.84219999999999995;
    troe_Tsss[11] = 125;
    troe_Ts[11]   = 2219;
    troe_Tss[11]  = 6882;
    troe_len[11]  = 4;
    prefactor_units[11]  = 1.0000000000000002e-06;
    activation_units[11] = 0.50321666580471969;
    phase_units[11]      = pow(10,-12.000000);
    is_PD[11] = 1;
    nTB[11] = 7;
    TB[11] = (amrex::Real *) malloc(7 * sizeof(amrex::Real));
    TBid[11] = (int *) malloc(7 * sizeof(int));
    TBid[11][0] = 0; TB[11][0] = 2; // H2
    TBid[11][1] = 5; TB[11][1] = 6; // H2O
    TBid[11][2] = 13; TB[11][2] = 2; // CH4
    TBid[11][3] = 14; TB[11][3] = 1.5; // CO
    TBid[11][4] = 15; TB[11][4] = 2; // CO2
    TBid[11][5] = 26; TB[11][5] = 3; // C2H6
    TBid[11][6] = 31; TB[11][6] = 0.69999999999999996; // AR

    // (76):  H + C2H5 <=> H2 + C2H4
    kiv[85] = {1,25,0,24};
    nuv[85] = {-1,-1,1,1};
    // (76):  H + C2H5 <=> H2 + C2H4
    fwd_A[85]     = 2000000000000;
    fwd_beta[85]  = 0;
    fwd_Ea[85]    = 0;
    prefactor_units[85]  = 1.0000000000000002e-06;
    activation_units[85] = 0.50321666580471969;
    phase_units[85]      = pow(10,-12.000000);
    is_PD[85] = 0;
    nTB[85] = 0;

    // (77):  H + C2H6 <=> C2H5 + H2
    kiv[86] = {1,26,25,0};
    nuv[86] = {-1,-1,1,1};
    // (77):  H + C2H6 <=> C2H5 + H2
    fwd_A[86]     = 115000000;
    fwd_beta[86]  = 1.8999999999999999;
    fwd_Ea[86]    = 7530;
    prefactor_units[86]  = 1.0000000000000002e-06;
    activation_units[86] = 0.50321666580471969;
    phase_units[86]      = pow(10,-12.000000);
    is_PD[86] = 0;
    nTB[86] = 0;

    // (78):  H + HCCO <=> CH2(S) + CO
    kiv[87] = {1,27,11,14};
    nuv[87] = {-1,-1,1,1};
    // (78):  H + HCCO <=> CH2(S) + CO
    fwd_A[87]     = 100000000000000;
    fwd_beta[87]  = 0;
    fwd_Ea[87]    = 0;
    prefactor_units[87]  = 1.0000000000000002e-06;
    activation_units[87] = 0.50321666580471969;
    phase_units[87]      = pow(10,-12.000000);
    is_PD[87] = 0;
    nTB[87] = 0;

    // (79):  H + CH2CO <=> HCCO + H2
    kiv[88] = {1,28,27,0};
    nuv[88] = {-1,-1,1,1};
    // (79):  H + CH2CO <=> HCCO + H2
    fwd_A[88]     = 50000000000000;
    fwd_beta[88]  = 0;
    fwd_Ea[88]    = 8000;
    prefactor_units[88]  = 1.0000000000000002e-06;
    activation_units[88] = 0.50321666580471969;
    phase_units[88]      = pow(10,-12.000000);
    is_PD[88] = 0;
    nTB[88] = 0;

    // (80):  H + CH2CO <=> CH3 + CO
    kiv[89] = {1,28,12,14};
    nuv[89] = {-1,-1,1,1};
    // (80):  H + CH2CO <=> CH3 + CO
    fwd_A[89]     = 11300000000000;
    fwd_beta[89]  = 0;
    fwd_Ea[89]    = 3428;
    prefactor_units[89]  = 1.0000000000000002e-06;
    activation_units[89] = 0.50321666580471969;
    phase_units[89]      = pow(10,-12.000000);
    is_PD[89] = 0;
    nTB[89] = 0;

    // (81):  H + HCCOH <=> H + CH2CO
    kiv[90] = {1,29,1,28};
    nuv[90] = {-1,-1,1,1};
    // (81):  H + HCCOH <=> H + CH2CO
    fwd_A[90]     = 10000000000000;
    fwd_beta[90]  = 0;
    fwd_Ea[90]    = 0;
    prefactor_units[90]  = 1.0000000000000002e-06;
    activation_units[90] = 0.50321666580471969;
    phase_units[90]      = pow(10,-12.000000);
    is_PD[90] = 0;
    nTB[90] = 0;

    // (82):  H2 + CO (+M) <=> CH2O (+M)
    kiv[12] = {0,14,17};
    nuv[12] = {-1,-1,1};
    // (82):  H2 + CO (+M) <=> CH2O (+M)
    fwd_A[12]     = 43000000;
    fwd_beta[12]  = 1.5;
    fwd_Ea[12]    = 79600;
    low_A[12]     = 5.0699999999999998e+27;
    low_beta[12]  = -3.4199999999999999;
    low_Ea[12]    = 84350;
    troe_a[12]    = 0.93200000000000005;
    troe_Tsss[12] = 197;
    troe_Ts[12]   = 1540;
    troe_Tss[12]  = 10300;
    troe_len[12]  = 4;
    prefactor_units[12]  = 1.0000000000000002e-06;
    activation_units[12] = 0.50321666580471969;
    phase_units[12]      = pow(10,-12.000000);
    is_PD[12] = 1;
    nTB[12] = 7;
    TB[12] = (amrex::Real *) malloc(7 * sizeof(amrex::Real));
    TBid[12] = (int *) malloc(7 * sizeof(int));
    TBid[12][0] = 0; TB[12][0] = 2; // H2
    TBid[12][1] = 5; TB[12][1] = 6; // H2O
    TBid[12][2] = 13; TB[12][2] = 2; // CH4
    TBid[12][3] = 14; TB[12][3] = 1.5; // CO
    TBid[12][4] = 15; TB[12][4] = 2; // CO2
    TBid[12][5] = 26; TB[12][5] = 3; // C2H6
    TBid[12][6] = 31; TB[12][6] = 0.69999999999999996; // AR

    // (83):  OH + H2 <=> H + H2O
    kiv[91] = {4,0,1,5};
    nuv[91] = {-1,-1,1,1};
    // (83):  OH + H2 <=> H + H2O
    fwd_A[91]     = 216000000;
    fwd_beta[91]  = 1.51;
    fwd_Ea[91]    = 3430;
    prefactor_units[91]  = 1.0000000000000002e-06;
    activation_units[91] = 0.50321666580471969;
    phase_units[91]      = pow(10,-12.000000);
    is_PD[91] = 0;
    nTB[91] = 0;

    // (84):  2.000000 OH (+M) <=> H2O2 (+M)
    kiv[13] = {4,7};
    nuv[13] = {-2.0,1};
    // (84):  2.000000 OH (+M) <=> H2O2 (+M)
    fwd_A[13]     = 74000000000000;
    fwd_beta[13]  = -0.37;
    fwd_Ea[13]    = 0;
    low_A[13]     = 2.3e+18;
    low_beta[13]  = -0.90000000000000002;
    low_Ea[13]    = -1700;
    troe_a[13]    = 0.73460000000000003;
    troe_Tsss[13] = 94;
    troe_Ts[13]   = 1756;
    troe_Tss[13]  = 5182;
    troe_len[13]  = 4;
    prefactor_units[13]  = 1.0000000000000002e-06;
    activation_units[13] = 0.50321666580471969;
    phase_units[13]      = pow(10,-12.000000);
    is_PD[13] = 1;
    nTB[13] = 7;
    TB[13] = (amrex::Real *) malloc(7 * sizeof(amrex::Real));
    TBid[13] = (int *) malloc(7 * sizeof(int));
    TBid[13][0] = 0; TB[13][0] = 2; // H2
    TBid[13][1] = 5; TB[13][1] = 6; // H2O
    TBid[13][2] = 13; TB[13][2] = 2; // CH4
    TBid[13][3] = 14; TB[13][3] = 1.5; // CO
    TBid[13][4] = 15; TB[13][4] = 2; // CO2
    TBid[13][5] = 26; TB[13][5] = 3; // C2H6
    TBid[13][6] = 31; TB[13][6] = 0.69999999999999996; // AR

    // (85):  2.000000 OH <=> O + H2O
    kiv[92] = {4,2,5};
    nuv[92] = {-2.0,1,1};
    // (85):  2.000000 OH <=> O + H2O
    fwd_A[92]     = 35700;
    fwd_beta[92]  = 2.3999999999999999;
    fwd_Ea[92]    = -2110;
    prefactor_units[92]  = 1.0000000000000002e-06;
    activation_units[92] = 0.50321666580471969;
    phase_units[92]      = pow(10,-12.000000);
    is_PD[92] = 0;
    nTB[92] = 0;

    // (86):  OH + HO2 <=> O2 + H2O
    kiv[93] = {4,6,3,5};
    nuv[93] = {-1,-1,1,1};
    // (86):  OH + HO2 <=> O2 + H2O
    fwd_A[93]     = 29000000000000;
    fwd_beta[93]  = 0;
    fwd_Ea[93]    = -500;
    prefactor_units[93]  = 1.0000000000000002e-06;
    activation_units[93] = 0.50321666580471969;
    phase_units[93]      = pow(10,-12.000000);
    is_PD[93] = 0;
    nTB[93] = 0;

    // (87):  OH + H2O2 <=> HO2 + H2O
    kiv[94] = {4,7,6,5};
    nuv[94] = {-1,-1,1,1};
    // (87):  OH + H2O2 <=> HO2 + H2O
    fwd_A[94]     = 1750000000000;
    fwd_beta[94]  = 0;
    fwd_Ea[94]    = 320;
    prefactor_units[94]  = 1.0000000000000002e-06;
    activation_units[94] = 0.50321666580471969;
    phase_units[94]      = pow(10,-12.000000);
    is_PD[94] = 0;
    nTB[94] = 0;

    // (88):  OH + H2O2 <=> HO2 + H2O
    kiv[95] = {4,7,6,5};
    nuv[95] = {-1,-1,1,1};
    // (88):  OH + H2O2 <=> HO2 + H2O
    fwd_A[95]     = 580000000000000;
    fwd_beta[95]  = 0;
    fwd_Ea[95]    = 9560;
    prefactor_units[95]  = 1.0000000000000002e-06;
    activation_units[95] = 0.50321666580471969;
    phase_units[95]      = pow(10,-12.000000);
    is_PD[95] = 0;
    nTB[95] = 0;

    // (89):  OH + C <=> H + CO
    kiv[96] = {4,8,1,14};
    nuv[96] = {-1,-1,1,1};
    // (89):  OH + C <=> H + CO
    fwd_A[96]     = 50000000000000;
    fwd_beta[96]  = 0;
    fwd_Ea[96]    = 0;
    prefactor_units[96]  = 1.0000000000000002e-06;
    activation_units[96] = 0.50321666580471969;
    phase_units[96]      = pow(10,-12.000000);
    is_PD[96] = 0;
    nTB[96] = 0;

    // (90):  OH + CH <=> H + HCO
    kiv[97] = {4,9,1,16};
    nuv[97] = {-1,-1,1,1};
    // (90):  OH + CH <=> H + HCO
    fwd_A[97]     = 30000000000000;
    fwd_beta[97]  = 0;
    fwd_Ea[97]    = 0;
    prefactor_units[97]  = 1.0000000000000002e-06;
    activation_units[97] = 0.50321666580471969;
    phase_units[97]      = pow(10,-12.000000);
    is_PD[97] = 0;
    nTB[97] = 0;

    // (91):  OH + CH2 <=> H + CH2O
    kiv[98] = {4,10,1,17};
    nuv[98] = {-1,-1,1,1};
    // (91):  OH + CH2 <=> H + CH2O
    fwd_A[98]     = 20000000000000;
    fwd_beta[98]  = 0;
    fwd_Ea[98]    = 0;
    prefactor_units[98]  = 1.0000000000000002e-06;
    activation_units[98] = 0.50321666580471969;
    phase_units[98]      = pow(10,-12.000000);
    is_PD[98] = 0;
    nTB[98] = 0;

    // (92):  OH + CH2 <=> CH + H2O
    kiv[99] = {4,10,9,5};
    nuv[99] = {-1,-1,1,1};
    // (92):  OH + CH2 <=> CH + H2O
    fwd_A[99]     = 11300000;
    fwd_beta[99]  = 2;
    fwd_Ea[99]    = 3000;
    prefactor_units[99]  = 1.0000000000000002e-06;
    activation_units[99] = 0.50321666580471969;
    phase_units[99]      = pow(10,-12.000000);
    is_PD[99] = 0;
    nTB[99] = 0;

    // (93):  OH + CH2(S) <=> H + CH2O
    kiv[100] = {4,11,1,17};
    nuv[100] = {-1,-1,1,1};
    // (93):  OH + CH2(S) <=> H + CH2O
    fwd_A[100]     = 30000000000000;
    fwd_beta[100]  = 0;
    fwd_Ea[100]    = 0;
    prefactor_units[100]  = 1.0000000000000002e-06;
    activation_units[100] = 0.50321666580471969;
    phase_units[100]      = pow(10,-12.000000);
    is_PD[100] = 0;
    nTB[100] = 0;

    // (94):  OH + CH3 (+M) <=> CH3OH (+M)
    kiv[14] = {4,12,20};
    nuv[14] = {-1,-1,1};
    // (94):  OH + CH3 (+M) <=> CH3OH (+M)
    fwd_A[14]     = 63000000000000;
    fwd_beta[14]  = 0;
    fwd_Ea[14]    = 0;
    low_A[14]     = 2.7e+38;
    low_beta[14]  = -6.2999999999999998;
    low_Ea[14]    = 3100;
    troe_a[14]    = 0.21049999999999999;
    troe_Tsss[14] = 83.5;
    troe_Ts[14]   = 5398;
    troe_Tss[14]  = 8370;
    troe_len[14]  = 4;
    prefactor_units[14]  = 1.0000000000000002e-06;
    activation_units[14] = 0.50321666580471969;
    phase_units[14]      = pow(10,-12.000000);
    is_PD[14] = 1;
    nTB[14] = 6;
    TB[14] = (amrex::Real *) malloc(6 * sizeof(amrex::Real));
    TBid[14] = (int *) malloc(6 * sizeof(int));
    TBid[14][0] = 0; TB[14][0] = 2; // H2
    TBid[14][1] = 5; TB[14][1] = 6; // H2O
    TBid[14][2] = 13; TB[14][2] = 2; // CH4
    TBid[14][3] = 14; TB[14][3] = 1.5; // CO
    TBid[14][4] = 15; TB[14][4] = 2; // CO2
    TBid[14][5] = 26; TB[14][5] = 3; // C2H6

    // (95):  OH + CH3 <=> CH2 + H2O
    kiv[101] = {4,12,10,5};
    nuv[101] = {-1,-1,1,1};
    // (95):  OH + CH3 <=> CH2 + H2O
    fwd_A[101]     = 56000000;
    fwd_beta[101]  = 1.6000000000000001;
    fwd_Ea[101]    = 5420;
    prefactor_units[101]  = 1.0000000000000002e-06;
    activation_units[101] = 0.50321666580471969;
    phase_units[101]      = pow(10,-12.000000);
    is_PD[101] = 0;
    nTB[101] = 0;

    // (96):  OH + CH3 <=> CH2(S) + H2O
    kiv[102] = {4,12,11,5};
    nuv[102] = {-1,-1,1,1};
    // (96):  OH + CH3 <=> CH2(S) + H2O
    fwd_A[102]     = 25010000000000;
    fwd_beta[102]  = 0;
    fwd_Ea[102]    = 0;
    prefactor_units[102]  = 1.0000000000000002e-06;
    activation_units[102] = 0.50321666580471969;
    phase_units[102]      = pow(10,-12.000000);
    is_PD[102] = 0;
    nTB[102] = 0;

    // (97):  OH + CH4 <=> CH3 + H2O
    kiv[103] = {4,13,12,5};
    nuv[103] = {-1,-1,1,1};
    // (97):  OH + CH4 <=> CH3 + H2O
    fwd_A[103]     = 100000000;
    fwd_beta[103]  = 1.6000000000000001;
    fwd_Ea[103]    = 3120;
    prefactor_units[103]  = 1.0000000000000002e-06;
    activation_units[103] = 0.50321666580471969;
    phase_units[103]      = pow(10,-12.000000);
    is_PD[103] = 0;
    nTB[103] = 0;

    // (98):  OH + CO <=> H + CO2
    kiv[104] = {4,14,1,15};
    nuv[104] = {-1,-1,1,1};
    // (98):  OH + CO <=> H + CO2
    fwd_A[104]     = 47600000;
    fwd_beta[104]  = 1.228;
    fwd_Ea[104]    = 70;
    prefactor_units[104]  = 1.0000000000000002e-06;
    activation_units[104] = 0.50321666580471969;
    phase_units[104]      = pow(10,-12.000000);
    is_PD[104] = 0;
    nTB[104] = 0;

    // (99):  OH + HCO <=> H2O + CO
    kiv[105] = {4,16,5,14};
    nuv[105] = {-1,-1,1,1};
    // (99):  OH + HCO <=> H2O + CO
    fwd_A[105]     = 50000000000000;
    fwd_beta[105]  = 0;
    fwd_Ea[105]    = 0;
    prefactor_units[105]  = 1.0000000000000002e-06;
    activation_units[105] = 0.50321666580471969;
    phase_units[105]      = pow(10,-12.000000);
    is_PD[105] = 0;
    nTB[105] = 0;

    // (100):  OH + CH2O <=> HCO + H2O
    kiv[106] = {4,17,16,5};
    nuv[106] = {-1,-1,1,1};
    // (100):  OH + CH2O <=> HCO + H2O
    fwd_A[106]     = 3430000000;
    fwd_beta[106]  = 1.1799999999999999;
    fwd_Ea[106]    = -447;
    prefactor_units[106]  = 1.0000000000000002e-06;
    activation_units[106] = 0.50321666580471969;
    phase_units[106]      = pow(10,-12.000000);
    is_PD[106] = 0;
    nTB[106] = 0;

    // (101):  OH + CH2OH <=> H2O + CH2O
    kiv[107] = {4,18,5,17};
    nuv[107] = {-1,-1,1,1};
    // (101):  OH + CH2OH <=> H2O + CH2O
    fwd_A[107]     = 5000000000000;
    fwd_beta[107]  = 0;
    fwd_Ea[107]    = 0;
    prefactor_units[107]  = 1.0000000000000002e-06;
    activation_units[107] = 0.50321666580471969;
    phase_units[107]      = pow(10,-12.000000);
    is_PD[107] = 0;
    nTB[107] = 0;

    // (102):  OH + CH3O <=> H2O + CH2O
    kiv[108] = {4,19,5,17};
    nuv[108] = {-1,-1,1,1};
    // (102):  OH + CH3O <=> H2O + CH2O
    fwd_A[108]     = 5000000000000;
    fwd_beta[108]  = 0;
    fwd_Ea[108]    = 0;
    prefactor_units[108]  = 1.0000000000000002e-06;
    activation_units[108] = 0.50321666580471969;
    phase_units[108]      = pow(10,-12.000000);
    is_PD[108] = 0;
    nTB[108] = 0;

    // (103):  OH + CH3OH <=> CH2OH + H2O
    kiv[109] = {4,20,18,5};
    nuv[109] = {-1,-1,1,1};
    // (103):  OH + CH3OH <=> CH2OH + H2O
    fwd_A[109]     = 1440000;
    fwd_beta[109]  = 2;
    fwd_Ea[109]    = -840;
    prefactor_units[109]  = 1.0000000000000002e-06;
    activation_units[109] = 0.50321666580471969;
    phase_units[109]      = pow(10,-12.000000);
    is_PD[109] = 0;
    nTB[109] = 0;

    // (104):  OH + CH3OH <=> CH3O + H2O
    kiv[110] = {4,20,19,5};
    nuv[110] = {-1,-1,1,1};
    // (104):  OH + CH3OH <=> CH3O + H2O
    fwd_A[110]     = 6300000;
    fwd_beta[110]  = 2;
    fwd_Ea[110]    = 1500;
    prefactor_units[110]  = 1.0000000000000002e-06;
    activation_units[110] = 0.50321666580471969;
    phase_units[110]      = pow(10,-12.000000);
    is_PD[110] = 0;
    nTB[110] = 0;

    // (105):  OH + C2H <=> H + HCCO
    kiv[111] = {4,21,1,27};
    nuv[111] = {-1,-1,1,1};
    // (105):  OH + C2H <=> H + HCCO
    fwd_A[111]     = 20000000000000;
    fwd_beta[111]  = 0;
    fwd_Ea[111]    = 0;
    prefactor_units[111]  = 1.0000000000000002e-06;
    activation_units[111] = 0.50321666580471969;
    phase_units[111]      = pow(10,-12.000000);
    is_PD[111] = 0;
    nTB[111] = 0;

    // (106):  OH + C2H2 <=> H + CH2CO
    kiv[112] = {4,22,1,28};
    nuv[112] = {-1,-1,1,1};
    // (106):  OH + C2H2 <=> H + CH2CO
    fwd_A[112]     = 0.00021800000000000001;
    fwd_beta[112]  = 4.5;
    fwd_Ea[112]    = -1000;
    prefactor_units[112]  = 1.0000000000000002e-06;
    activation_units[112] = 0.50321666580471969;
    phase_units[112]      = pow(10,-12.000000);
    is_PD[112] = 0;
    nTB[112] = 0;

    // (107):  OH + C2H2 <=> H + HCCOH
    kiv[113] = {4,22,1,29};
    nuv[113] = {-1,-1,1,1};
    // (107):  OH + C2H2 <=> H + HCCOH
    fwd_A[113]     = 504000;
    fwd_beta[113]  = 2.2999999999999998;
    fwd_Ea[113]    = 13500;
    prefactor_units[113]  = 1.0000000000000002e-06;
    activation_units[113] = 0.50321666580471969;
    phase_units[113]      = pow(10,-12.000000);
    is_PD[113] = 0;
    nTB[113] = 0;

    // (108):  OH + C2H2 <=> C2H + H2O
    kiv[114] = {4,22,21,5};
    nuv[114] = {-1,-1,1,1};
    // (108):  OH + C2H2 <=> C2H + H2O
    fwd_A[114]     = 33700000;
    fwd_beta[114]  = 2;
    fwd_Ea[114]    = 14000;
    prefactor_units[114]  = 1.0000000000000002e-06;
    activation_units[114] = 0.50321666580471969;
    phase_units[114]      = pow(10,-12.000000);
    is_PD[114] = 0;
    nTB[114] = 0;

    // (109):  OH + C2H2 <=> CH3 + CO
    kiv[115] = {4,22,12,14};
    nuv[115] = {-1,-1,1,1};
    // (109):  OH + C2H2 <=> CH3 + CO
    fwd_A[115]     = 0.00048299999999999998;
    fwd_beta[115]  = 4;
    fwd_Ea[115]    = -2000;
    prefactor_units[115]  = 1.0000000000000002e-06;
    activation_units[115] = 0.50321666580471969;
    phase_units[115]      = pow(10,-12.000000);
    is_PD[115] = 0;
    nTB[115] = 0;

    // (110):  OH + C2H3 <=> H2O + C2H2
    kiv[116] = {4,23,5,22};
    nuv[116] = {-1,-1,1,1};
    // (110):  OH + C2H3 <=> H2O + C2H2
    fwd_A[116]     = 5000000000000;
    fwd_beta[116]  = 0;
    fwd_Ea[116]    = 0;
    prefactor_units[116]  = 1.0000000000000002e-06;
    activation_units[116] = 0.50321666580471969;
    phase_units[116]      = pow(10,-12.000000);
    is_PD[116] = 0;
    nTB[116] = 0;

    // (111):  OH + C2H4 <=> C2H3 + H2O
    kiv[117] = {4,24,23,5};
    nuv[117] = {-1,-1,1,1};
    // (111):  OH + C2H4 <=> C2H3 + H2O
    fwd_A[117]     = 3600000;
    fwd_beta[117]  = 2;
    fwd_Ea[117]    = 2500;
    prefactor_units[117]  = 1.0000000000000002e-06;
    activation_units[117] = 0.50321666580471969;
    phase_units[117]      = pow(10,-12.000000);
    is_PD[117] = 0;
    nTB[117] = 0;

    // (112):  OH + C2H6 <=> C2H5 + H2O
    kiv[118] = {4,26,25,5};
    nuv[118] = {-1,-1,1,1};
    // (112):  OH + C2H6 <=> C2H5 + H2O
    fwd_A[118]     = 3540000;
    fwd_beta[118]  = 2.1200000000000001;
    fwd_Ea[118]    = 870;
    prefactor_units[118]  = 1.0000000000000002e-06;
    activation_units[118] = 0.50321666580471969;
    phase_units[118]      = pow(10,-12.000000);
    is_PD[118] = 0;
    nTB[118] = 0;

    // (113):  OH + CH2CO <=> HCCO + H2O
    kiv[119] = {4,28,27,5};
    nuv[119] = {-1,-1,1,1};
    // (113):  OH + CH2CO <=> HCCO + H2O
    fwd_A[119]     = 7500000000000;
    fwd_beta[119]  = 0;
    fwd_Ea[119]    = 2000;
    prefactor_units[119]  = 1.0000000000000002e-06;
    activation_units[119] = 0.50321666580471969;
    phase_units[119]      = pow(10,-12.000000);
    is_PD[119] = 0;
    nTB[119] = 0;

    // (114):  2.000000 HO2 <=> O2 + H2O2
    kiv[120] = {6,3,7};
    nuv[120] = {-2.0,1,1};
    // (114):  2.000000 HO2 <=> O2 + H2O2
    fwd_A[120]     = 130000000000;
    fwd_beta[120]  = 0;
    fwd_Ea[120]    = -1630;
    prefactor_units[120]  = 1.0000000000000002e-06;
    activation_units[120] = 0.50321666580471969;
    phase_units[120]      = pow(10,-12.000000);
    is_PD[120] = 0;
    nTB[120] = 0;

    // (115):  2.000000 HO2 <=> O2 + H2O2
    kiv[121] = {6,3,7};
    nuv[121] = {-2.0,1,1};
    // (115):  2.000000 HO2 <=> O2 + H2O2
    fwd_A[121]     = 420000000000000;
    fwd_beta[121]  = 0;
    fwd_Ea[121]    = 12000;
    prefactor_units[121]  = 1.0000000000000002e-06;
    activation_units[121] = 0.50321666580471969;
    phase_units[121]      = pow(10,-12.000000);
    is_PD[121] = 0;
    nTB[121] = 0;

    // (116):  HO2 + CH2 <=> OH + CH2O
    kiv[122] = {6,10,4,17};
    nuv[122] = {-1,-1,1,1};
    // (116):  HO2 + CH2 <=> OH + CH2O
    fwd_A[122]     = 20000000000000;
    fwd_beta[122]  = 0;
    fwd_Ea[122]    = 0;
    prefactor_units[122]  = 1.0000000000000002e-06;
    activation_units[122] = 0.50321666580471969;
    phase_units[122]      = pow(10,-12.000000);
    is_PD[122] = 0;
    nTB[122] = 0;

    // (117):  HO2 + CH3 <=> O2 + CH4
    kiv[123] = {6,12,3,13};
    nuv[123] = {-1,-1,1,1};
    // (117):  HO2 + CH3 <=> O2 + CH4
    fwd_A[123]     = 1000000000000;
    fwd_beta[123]  = 0;
    fwd_Ea[123]    = 0;
    prefactor_units[123]  = 1.0000000000000002e-06;
    activation_units[123] = 0.50321666580471969;
    phase_units[123]      = pow(10,-12.000000);
    is_PD[123] = 0;
    nTB[123] = 0;

    // (118):  HO2 + CH3 <=> OH + CH3O
    kiv[124] = {6,12,4,19};
    nuv[124] = {-1,-1,1,1};
    // (118):  HO2 + CH3 <=> OH + CH3O
    fwd_A[124]     = 20000000000000;
    fwd_beta[124]  = 0;
    fwd_Ea[124]    = 0;
    prefactor_units[124]  = 1.0000000000000002e-06;
    activation_units[124] = 0.50321666580471969;
    phase_units[124]      = pow(10,-12.000000);
    is_PD[124] = 0;
    nTB[124] = 0;

    // (119):  HO2 + CO <=> OH + CO2
    kiv[125] = {6,14,4,15};
    nuv[125] = {-1,-1,1,1};
    // (119):  HO2 + CO <=> OH + CO2
    fwd_A[125]     = 150000000000000;
    fwd_beta[125]  = 0;
    fwd_Ea[125]    = 23600;
    prefactor_units[125]  = 1.0000000000000002e-06;
    activation_units[125] = 0.50321666580471969;
    phase_units[125]      = pow(10,-12.000000);
    is_PD[125] = 0;
    nTB[125] = 0;

    // (120):  HO2 + CH2O <=> HCO + H2O2
    kiv[126] = {6,17,16,7};
    nuv[126] = {-1,-1,1,1};
    // (120):  HO2 + CH2O <=> HCO + H2O2
    fwd_A[126]     = 1000000000000;
    fwd_beta[126]  = 0;
    fwd_Ea[126]    = 8000;
    prefactor_units[126]  = 1.0000000000000002e-06;
    activation_units[126] = 0.50321666580471969;
    phase_units[126]      = pow(10,-12.000000);
    is_PD[126] = 0;
    nTB[126] = 0;

    // (121):  C + O2 <=> O + CO
    kiv[127] = {8,3,2,14};
    nuv[127] = {-1,-1,1,1};
    // (121):  C + O2 <=> O + CO
    fwd_A[127]     = 58000000000000;
    fwd_beta[127]  = 0;
    fwd_Ea[127]    = 576;
    prefactor_units[127]  = 1.0000000000000002e-06;
    activation_units[127] = 0.50321666580471969;
    phase_units[127]      = pow(10,-12.000000);
    is_PD[127] = 0;
    nTB[127] = 0;

    // (122):  C + CH2 <=> H + C2H
    kiv[128] = {8,10,1,21};
    nuv[128] = {-1,-1,1,1};
    // (122):  C + CH2 <=> H + C2H
    fwd_A[128]     = 50000000000000;
    fwd_beta[128]  = 0;
    fwd_Ea[128]    = 0;
    prefactor_units[128]  = 1.0000000000000002e-06;
    activation_units[128] = 0.50321666580471969;
    phase_units[128]      = pow(10,-12.000000);
    is_PD[128] = 0;
    nTB[128] = 0;

    // (123):  C + CH3 <=> H + C2H2
    kiv[129] = {8,12,1,22};
    nuv[129] = {-1,-1,1,1};
    // (123):  C + CH3 <=> H + C2H2
    fwd_A[129]     = 50000000000000;
    fwd_beta[129]  = 0;
    fwd_Ea[129]    = 0;
    prefactor_units[129]  = 1.0000000000000002e-06;
    activation_units[129] = 0.50321666580471969;
    phase_units[129]      = pow(10,-12.000000);
    is_PD[129] = 0;
    nTB[129] = 0;

    // (124):  CH + O2 <=> O + HCO
    kiv[130] = {9,3,2,16};
    nuv[130] = {-1,-1,1,1};
    // (124):  CH + O2 <=> O + HCO
    fwd_A[130]     = 33000000000000;
    fwd_beta[130]  = 0;
    fwd_Ea[130]    = 0;
    prefactor_units[130]  = 1.0000000000000002e-06;
    activation_units[130] = 0.50321666580471969;
    phase_units[130]      = pow(10,-12.000000);
    is_PD[130] = 0;
    nTB[130] = 0;

    // (125):  CH + H2 <=> H + CH2
    kiv[131] = {9,0,1,10};
    nuv[131] = {-1,-1,1,1};
    // (125):  CH + H2 <=> H + CH2
    fwd_A[131]     = 110700000;
    fwd_beta[131]  = 1.79;
    fwd_Ea[131]    = 1670;
    prefactor_units[131]  = 1.0000000000000002e-06;
    activation_units[131] = 0.50321666580471969;
    phase_units[131]      = pow(10,-12.000000);
    is_PD[131] = 0;
    nTB[131] = 0;

    // (126):  CH + H2O <=> H + CH2O
    kiv[132] = {9,5,1,17};
    nuv[132] = {-1,-1,1,1};
    // (126):  CH + H2O <=> H + CH2O
    fwd_A[132]     = 5710000000000;
    fwd_beta[132]  = 0;
    fwd_Ea[132]    = -755;
    prefactor_units[132]  = 1.0000000000000002e-06;
    activation_units[132] = 0.50321666580471969;
    phase_units[132]      = pow(10,-12.000000);
    is_PD[132] = 0;
    nTB[132] = 0;

    // (127):  CH + CH2 <=> H + C2H2
    kiv[133] = {9,10,1,22};
    nuv[133] = {-1,-1,1,1};
    // (127):  CH + CH2 <=> H + C2H2
    fwd_A[133]     = 40000000000000;
    fwd_beta[133]  = 0;
    fwd_Ea[133]    = 0;
    prefactor_units[133]  = 1.0000000000000002e-06;
    activation_units[133] = 0.50321666580471969;
    phase_units[133]      = pow(10,-12.000000);
    is_PD[133] = 0;
    nTB[133] = 0;

    // (128):  CH + CH3 <=> H + C2H3
    kiv[134] = {9,12,1,23};
    nuv[134] = {-1,-1,1,1};
    // (128):  CH + CH3 <=> H + C2H3
    fwd_A[134]     = 30000000000000;
    fwd_beta[134]  = 0;
    fwd_Ea[134]    = 0;
    prefactor_units[134]  = 1.0000000000000002e-06;
    activation_units[134] = 0.50321666580471969;
    phase_units[134]      = pow(10,-12.000000);
    is_PD[134] = 0;
    nTB[134] = 0;

    // (129):  CH + CH4 <=> H + C2H4
    kiv[135] = {9,13,1,24};
    nuv[135] = {-1,-1,1,1};
    // (129):  CH + CH4 <=> H + C2H4
    fwd_A[135]     = 60000000000000;
    fwd_beta[135]  = 0;
    fwd_Ea[135]    = 0;
    prefactor_units[135]  = 1.0000000000000002e-06;
    activation_units[135] = 0.50321666580471969;
    phase_units[135]      = pow(10,-12.000000);
    is_PD[135] = 0;
    nTB[135] = 0;

    // (130):  CH + CO (+M) <=> HCCO (+M)
    kiv[15] = {9,14,27};
    nuv[15] = {-1,-1,1};
    // (130):  CH + CO (+M) <=> HCCO (+M)
    fwd_A[15]     = 50000000000000;
    fwd_beta[15]  = 0;
    fwd_Ea[15]    = 0;
    low_A[15]     = 2.6899999999999998e+28;
    low_beta[15]  = -3.7400000000000002;
    low_Ea[15]    = 1936;
    troe_a[15]    = 0.57569999999999999;
    troe_Tsss[15] = 237;
    troe_Ts[15]   = 1652;
    troe_Tss[15]  = 5069;
    troe_len[15]  = 4;
    prefactor_units[15]  = 1.0000000000000002e-06;
    activation_units[15] = 0.50321666580471969;
    phase_units[15]      = pow(10,-12.000000);
    is_PD[15] = 1;
    nTB[15] = 7;
    TB[15] = (amrex::Real *) malloc(7 * sizeof(amrex::Real));
    TBid[15] = (int *) malloc(7 * sizeof(int));
    TBid[15][0] = 0; TB[15][0] = 2; // H2
    TBid[15][1] = 5; TB[15][1] = 6; // H2O
    TBid[15][2] = 13; TB[15][2] = 2; // CH4
    TBid[15][3] = 14; TB[15][3] = 1.5; // CO
    TBid[15][4] = 15; TB[15][4] = 2; // CO2
    TBid[15][5] = 26; TB[15][5] = 3; // C2H6
    TBid[15][6] = 31; TB[15][6] = 0.69999999999999996; // AR

    // (131):  CH + CO2 <=> HCO + CO
    kiv[136] = {9,15,16,14};
    nuv[136] = {-1,-1,1,1};
    // (131):  CH + CO2 <=> HCO + CO
    fwd_A[136]     = 3400000000000;
    fwd_beta[136]  = 0;
    fwd_Ea[136]    = 690;
    prefactor_units[136]  = 1.0000000000000002e-06;
    activation_units[136] = 0.50321666580471969;
    phase_units[136]      = pow(10,-12.000000);
    is_PD[136] = 0;
    nTB[136] = 0;

    // (132):  CH + CH2O <=> H + CH2CO
    kiv[137] = {9,17,1,28};
    nuv[137] = {-1,-1,1,1};
    // (132):  CH + CH2O <=> H + CH2CO
    fwd_A[137]     = 94600000000000;
    fwd_beta[137]  = 0;
    fwd_Ea[137]    = -515;
    prefactor_units[137]  = 1.0000000000000002e-06;
    activation_units[137] = 0.50321666580471969;
    phase_units[137]      = pow(10,-12.000000);
    is_PD[137] = 0;
    nTB[137] = 0;

    // (133):  CH + HCCO <=> CO + C2H2
    kiv[138] = {9,27,14,22};
    nuv[138] = {-1,-1,1,1};
    // (133):  CH + HCCO <=> CO + C2H2
    fwd_A[138]     = 50000000000000;
    fwd_beta[138]  = 0;
    fwd_Ea[138]    = 0;
    prefactor_units[138]  = 1.0000000000000002e-06;
    activation_units[138] = 0.50321666580471969;
    phase_units[138]      = pow(10,-12.000000);
    is_PD[138] = 0;
    nTB[138] = 0;

    // (134):  CH2 + O2 <=> OH + HCO
    kiv[139] = {10,3,4,16};
    nuv[139] = {-1,-1,1,1};
    // (134):  CH2 + O2 <=> OH + HCO
    fwd_A[139]     = 13200000000000;
    fwd_beta[139]  = 0;
    fwd_Ea[139]    = 1500;
    prefactor_units[139]  = 1.0000000000000002e-06;
    activation_units[139] = 0.50321666580471969;
    phase_units[139]      = pow(10,-12.000000);
    is_PD[139] = 0;
    nTB[139] = 0;

    // (135):  CH2 + H2 <=> H + CH3
    kiv[140] = {10,0,1,12};
    nuv[140] = {-1,-1,1,1};
    // (135):  CH2 + H2 <=> H + CH3
    fwd_A[140]     = 500000;
    fwd_beta[140]  = 2;
    fwd_Ea[140]    = 7230;
    prefactor_units[140]  = 1.0000000000000002e-06;
    activation_units[140] = 0.50321666580471969;
    phase_units[140]      = pow(10,-12.000000);
    is_PD[140] = 0;
    nTB[140] = 0;

    // (136):  2.000000 CH2 <=> H2 + C2H2
    kiv[141] = {10,0,22};
    nuv[141] = {-2.0,1,1};
    // (136):  2.000000 CH2 <=> H2 + C2H2
    fwd_A[141]     = 32000000000000;
    fwd_beta[141]  = 0;
    fwd_Ea[141]    = 0;
    prefactor_units[141]  = 1.0000000000000002e-06;
    activation_units[141] = 0.50321666580471969;
    phase_units[141]      = pow(10,-12.000000);
    is_PD[141] = 0;
    nTB[141] = 0;

    // (137):  CH2 + CH3 <=> H + C2H4
    kiv[142] = {10,12,1,24};
    nuv[142] = {-1,-1,1,1};
    // (137):  CH2 + CH3 <=> H + C2H4
    fwd_A[142]     = 40000000000000;
    fwd_beta[142]  = 0;
    fwd_Ea[142]    = 0;
    prefactor_units[142]  = 1.0000000000000002e-06;
    activation_units[142] = 0.50321666580471969;
    phase_units[142]      = pow(10,-12.000000);
    is_PD[142] = 0;
    nTB[142] = 0;

    // (138):  CH2 + CH4 <=> 2.000000 CH3
    kiv[143] = {10,13,12};
    nuv[143] = {-1,-1,2.0};
    // (138):  CH2 + CH4 <=> 2.000000 CH3
    fwd_A[143]     = 2460000;
    fwd_beta[143]  = 2;
    fwd_Ea[143]    = 8270;
    prefactor_units[143]  = 1.0000000000000002e-06;
    activation_units[143] = 0.50321666580471969;
    phase_units[143]      = pow(10,-12.000000);
    is_PD[143] = 0;
    nTB[143] = 0;

    // (139):  CH2 + CO (+M) <=> CH2CO (+M)
    kiv[16] = {10,14,28};
    nuv[16] = {-1,-1,1};
    // (139):  CH2 + CO (+M) <=> CH2CO (+M)
    fwd_A[16]     = 810000000000;
    fwd_beta[16]  = 0.5;
    fwd_Ea[16]    = 4510;
    low_A[16]     = 2.69e+33;
    low_beta[16]  = -5.1100000000000003;
    low_Ea[16]    = 7095;
    troe_a[16]    = 0.5907;
    troe_Tsss[16] = 275;
    troe_Ts[16]   = 1226;
    troe_Tss[16]  = 5185;
    troe_len[16]  = 4;
    prefactor_units[16]  = 1.0000000000000002e-06;
    activation_units[16] = 0.50321666580471969;
    phase_units[16]      = pow(10,-12.000000);
    is_PD[16] = 1;
    nTB[16] = 7;
    TB[16] = (amrex::Real *) malloc(7 * sizeof(amrex::Real));
    TBid[16] = (int *) malloc(7 * sizeof(int));
    TBid[16][0] = 0; TB[16][0] = 2; // H2
    TBid[16][1] = 5; TB[16][1] = 6; // H2O
    TBid[16][2] = 13; TB[16][2] = 2; // CH4
    TBid[16][3] = 14; TB[16][3] = 1.5; // CO
    TBid[16][4] = 15; TB[16][4] = 2; // CO2
    TBid[16][5] = 26; TB[16][5] = 3; // C2H6
    TBid[16][6] = 31; TB[16][6] = 0.69999999999999996; // AR

    // (140):  CH2 + HCCO <=> C2H3 + CO
    kiv[144] = {10,27,23,14};
    nuv[144] = {-1,-1,1,1};
    // (140):  CH2 + HCCO <=> C2H3 + CO
    fwd_A[144]     = 30000000000000;
    fwd_beta[144]  = 0;
    fwd_Ea[144]    = 0;
    prefactor_units[144]  = 1.0000000000000002e-06;
    activation_units[144] = 0.50321666580471969;
    phase_units[144]      = pow(10,-12.000000);
    is_PD[144] = 0;
    nTB[144] = 0;

    // (141):  CH2(S) + N2 <=> CH2 + N2
    kiv[145] = {11,30,10,30};
    nuv[145] = {-1,-1,1,1};
    // (141):  CH2(S) + N2 <=> CH2 + N2
    fwd_A[145]     = 15000000000000;
    fwd_beta[145]  = 0;
    fwd_Ea[145]    = 600;
    prefactor_units[145]  = 1.0000000000000002e-06;
    activation_units[145] = 0.50321666580471969;
    phase_units[145]      = pow(10,-12.000000);
    is_PD[145] = 0;
    nTB[145] = 0;

    // (142):  CH2(S) + AR <=> CH2 + AR
    kiv[146] = {11,31,10,31};
    nuv[146] = {-1,-1,1,1};
    // (142):  CH2(S) + AR <=> CH2 + AR
    fwd_A[146]     = 9000000000000;
    fwd_beta[146]  = 0;
    fwd_Ea[146]    = 600;
    prefactor_units[146]  = 1.0000000000000002e-06;
    activation_units[146] = 0.50321666580471969;
    phase_units[146]      = pow(10,-12.000000);
    is_PD[146] = 0;
    nTB[146] = 0;

    // (143):  CH2(S) + O2 <=> H + OH + CO
    kiv[147] = {11,3,1,4,14};
    nuv[147] = {-1,-1,1,1,1};
    // (143):  CH2(S) + O2 <=> H + OH + CO
    fwd_A[147]     = 28000000000000;
    fwd_beta[147]  = 0;
    fwd_Ea[147]    = 0;
    prefactor_units[147]  = 1.0000000000000002e-06;
    activation_units[147] = 0.50321666580471969;
    phase_units[147]      = pow(10,-12.000000);
    is_PD[147] = 0;
    nTB[147] = 0;

    // (144):  CH2(S) + O2 <=> CO + H2O
    kiv[148] = {11,3,14,5};
    nuv[148] = {-1,-1,1,1};
    // (144):  CH2(S) + O2 <=> CO + H2O
    fwd_A[148]     = 12000000000000;
    fwd_beta[148]  = 0;
    fwd_Ea[148]    = 0;
    prefactor_units[148]  = 1.0000000000000002e-06;
    activation_units[148] = 0.50321666580471969;
    phase_units[148]      = pow(10,-12.000000);
    is_PD[148] = 0;
    nTB[148] = 0;

    // (145):  CH2(S) + H2 <=> CH3 + H
    kiv[149] = {11,0,12,1};
    nuv[149] = {-1,-1,1,1};
    // (145):  CH2(S) + H2 <=> CH3 + H
    fwd_A[149]     = 70000000000000;
    fwd_beta[149]  = 0;
    fwd_Ea[149]    = 0;
    prefactor_units[149]  = 1.0000000000000002e-06;
    activation_units[149] = 0.50321666580471969;
    phase_units[149]      = pow(10,-12.000000);
    is_PD[149] = 0;
    nTB[149] = 0;

    // (146):  CH2(S) + H2O (+M) <=> CH3OH (+M)
    kiv[17] = {11,5,20};
    nuv[17] = {-1,-1,1};
    // (146):  CH2(S) + H2O (+M) <=> CH3OH (+M)
    fwd_A[17]     = 20000000000000;
    fwd_beta[17]  = 0;
    fwd_Ea[17]    = 0;
    low_A[17]     = 2.7e+38;
    low_beta[17]  = -6.2999999999999998;
    low_Ea[17]    = 3100;
    troe_a[17]    = 0.1507;
    troe_Tsss[17] = 134;
    troe_Ts[17]   = 2383;
    troe_Tss[17]  = 7265;
    troe_len[17]  = 4;
    prefactor_units[17]  = 1.0000000000000002e-06;
    activation_units[17] = 0.50321666580471969;
    phase_units[17]      = pow(10,-12.000000);
    is_PD[17] = 1;
    nTB[17] = 6;
    TB[17] = (amrex::Real *) malloc(6 * sizeof(amrex::Real));
    TBid[17] = (int *) malloc(6 * sizeof(int));
    TBid[17][0] = 0; TB[17][0] = 2; // H2
    TBid[17][1] = 5; TB[17][1] = 6; // H2O
    TBid[17][2] = 13; TB[17][2] = 2; // CH4
    TBid[17][3] = 14; TB[17][3] = 1.5; // CO
    TBid[17][4] = 15; TB[17][4] = 2; // CO2
    TBid[17][5] = 26; TB[17][5] = 3; // C2H6

    // (147):  CH2(S) + H2O <=> CH2 + H2O
    kiv[150] = {11,5,10,5};
    nuv[150] = {-1,-1,1,1};
    // (147):  CH2(S) + H2O <=> CH2 + H2O
    fwd_A[150]     = 30000000000000;
    fwd_beta[150]  = 0;
    fwd_Ea[150]    = 0;
    prefactor_units[150]  = 1.0000000000000002e-06;
    activation_units[150] = 0.50321666580471969;
    phase_units[150]      = pow(10,-12.000000);
    is_PD[150] = 0;
    nTB[150] = 0;

    // (148):  CH2(S) + CH3 <=> H + C2H4
    kiv[151] = {11,12,1,24};
    nuv[151] = {-1,-1,1,1};
    // (148):  CH2(S) + CH3 <=> H + C2H4
    fwd_A[151]     = 12000000000000;
    fwd_beta[151]  = 0;
    fwd_Ea[151]    = -570;
    prefactor_units[151]  = 1.0000000000000002e-06;
    activation_units[151] = 0.50321666580471969;
    phase_units[151]      = pow(10,-12.000000);
    is_PD[151] = 0;
    nTB[151] = 0;

    // (149):  CH2(S) + CH4 <=> 2.000000 CH3
    kiv[152] = {11,13,12};
    nuv[152] = {-1,-1,2.0};
    // (149):  CH2(S) + CH4 <=> 2.000000 CH3
    fwd_A[152]     = 16000000000000;
    fwd_beta[152]  = 0;
    fwd_Ea[152]    = -570;
    prefactor_units[152]  = 1.0000000000000002e-06;
    activation_units[152] = 0.50321666580471969;
    phase_units[152]      = pow(10,-12.000000);
    is_PD[152] = 0;
    nTB[152] = 0;

    // (150):  CH2(S) + CO <=> CH2 + CO
    kiv[153] = {11,14,10,14};
    nuv[153] = {-1,-1,1,1};
    // (150):  CH2(S) + CO <=> CH2 + CO
    fwd_A[153]     = 9000000000000;
    fwd_beta[153]  = 0;
    fwd_Ea[153]    = 0;
    prefactor_units[153]  = 1.0000000000000002e-06;
    activation_units[153] = 0.50321666580471969;
    phase_units[153]      = pow(10,-12.000000);
    is_PD[153] = 0;
    nTB[153] = 0;

    // (151):  CH2(S) + CO2 <=> CH2 + CO2
    kiv[154] = {11,15,10,15};
    nuv[154] = {-1,-1,1,1};
    // (151):  CH2(S) + CO2 <=> CH2 + CO2
    fwd_A[154]     = 7000000000000;
    fwd_beta[154]  = 0;
    fwd_Ea[154]    = 0;
    prefactor_units[154]  = 1.0000000000000002e-06;
    activation_units[154] = 0.50321666580471969;
    phase_units[154]      = pow(10,-12.000000);
    is_PD[154] = 0;
    nTB[154] = 0;

    // (152):  CH2(S) + CO2 <=> CO + CH2O
    kiv[155] = {11,15,14,17};
    nuv[155] = {-1,-1,1,1};
    // (152):  CH2(S) + CO2 <=> CO + CH2O
    fwd_A[155]     = 14000000000000;
    fwd_beta[155]  = 0;
    fwd_Ea[155]    = 0;
    prefactor_units[155]  = 1.0000000000000002e-06;
    activation_units[155] = 0.50321666580471969;
    phase_units[155]      = pow(10,-12.000000);
    is_PD[155] = 0;
    nTB[155] = 0;

    // (153):  CH2(S) + C2H6 <=> CH3 + C2H5
    kiv[156] = {11,26,12,25};
    nuv[156] = {-1,-1,1,1};
    // (153):  CH2(S) + C2H6 <=> CH3 + C2H5
    fwd_A[156]     = 40000000000000;
    fwd_beta[156]  = 0;
    fwd_Ea[156]    = -550;
    prefactor_units[156]  = 1.0000000000000002e-06;
    activation_units[156] = 0.50321666580471969;
    phase_units[156]      = pow(10,-12.000000);
    is_PD[156] = 0;
    nTB[156] = 0;

    // (154):  CH3 + O2 <=> O + CH3O
    kiv[157] = {12,3,2,19};
    nuv[157] = {-1,-1,1,1};
    // (154):  CH3 + O2 <=> O + CH3O
    fwd_A[157]     = 26750000000000;
    fwd_beta[157]  = 0;
    fwd_Ea[157]    = 28800;
    prefactor_units[157]  = 1.0000000000000002e-06;
    activation_units[157] = 0.50321666580471969;
    phase_units[157]      = pow(10,-12.000000);
    is_PD[157] = 0;
    nTB[157] = 0;

    // (155):  CH3 + O2 <=> OH + CH2O
    kiv[158] = {12,3,4,17};
    nuv[158] = {-1,-1,1,1};
    // (155):  CH3 + O2 <=> OH + CH2O
    fwd_A[158]     = 36000000000;
    fwd_beta[158]  = 0;
    fwd_Ea[158]    = 8940;
    prefactor_units[158]  = 1.0000000000000002e-06;
    activation_units[158] = 0.50321666580471969;
    phase_units[158]      = pow(10,-12.000000);
    is_PD[158] = 0;
    nTB[158] = 0;

    // (156):  CH3 + H2O2 <=> HO2 + CH4
    kiv[159] = {12,7,6,13};
    nuv[159] = {-1,-1,1,1};
    // (156):  CH3 + H2O2 <=> HO2 + CH4
    fwd_A[159]     = 24500;
    fwd_beta[159]  = 2.4700000000000002;
    fwd_Ea[159]    = 5180;
    prefactor_units[159]  = 1.0000000000000002e-06;
    activation_units[159] = 0.50321666580471969;
    phase_units[159]      = pow(10,-12.000000);
    is_PD[159] = 0;
    nTB[159] = 0;

    // (157):  2.000000 CH3 (+M) <=> C2H6 (+M)
    kiv[18] = {12,26};
    nuv[18] = {-2.0,1};
    // (157):  2.000000 CH3 (+M) <=> C2H6 (+M)
    fwd_A[18]     = 21200000000000000;
    fwd_beta[18]  = -0.96999999999999997;
    fwd_Ea[18]    = 620;
    low_A[18]     = 1.7700000000000001e+50;
    low_beta[18]  = -9.6699999999999999;
    low_Ea[18]    = 6220;
    troe_a[18]    = 0.53249999999999997;
    troe_Tsss[18] = 151;
    troe_Ts[18]   = 1038;
    troe_Tss[18]  = 4970;
    troe_len[18]  = 4;
    prefactor_units[18]  = 1.0000000000000002e-06;
    activation_units[18] = 0.50321666580471969;
    phase_units[18]      = pow(10,-12.000000);
    is_PD[18] = 1;
    nTB[18] = 7;
    TB[18] = (amrex::Real *) malloc(7 * sizeof(amrex::Real));
    TBid[18] = (int *) malloc(7 * sizeof(int));
    TBid[18][0] = 0; TB[18][0] = 2; // H2
    TBid[18][1] = 5; TB[18][1] = 6; // H2O
    TBid[18][2] = 13; TB[18][2] = 2; // CH4
    TBid[18][3] = 14; TB[18][3] = 1.5; // CO
    TBid[18][4] = 15; TB[18][4] = 2; // CO2
    TBid[18][5] = 26; TB[18][5] = 3; // C2H6
    TBid[18][6] = 31; TB[18][6] = 0.69999999999999996; // AR

    // (158):  2.000000 CH3 <=> H + C2H5
    kiv[160] = {12,1,25};
    nuv[160] = {-2.0,1,1};
    // (158):  2.000000 CH3 <=> H + C2H5
    fwd_A[160]     = 4990000000000;
    fwd_beta[160]  = 0.10000000000000001;
    fwd_Ea[160]    = 10600;
    prefactor_units[160]  = 1.0000000000000002e-06;
    activation_units[160] = 0.50321666580471969;
    phase_units[160]      = pow(10,-12.000000);
    is_PD[160] = 0;
    nTB[160] = 0;

    // (159):  CH3 + HCO <=> CH4 + CO
    kiv[161] = {12,16,13,14};
    nuv[161] = {-1,-1,1,1};
    // (159):  CH3 + HCO <=> CH4 + CO
    fwd_A[161]     = 26480000000000;
    fwd_beta[161]  = 0;
    fwd_Ea[161]    = 0;
    prefactor_units[161]  = 1.0000000000000002e-06;
    activation_units[161] = 0.50321666580471969;
    phase_units[161]      = pow(10,-12.000000);
    is_PD[161] = 0;
    nTB[161] = 0;

    // (160):  CH3 + CH2O <=> HCO + CH4
    kiv[162] = {12,17,16,13};
    nuv[162] = {-1,-1,1,1};
    // (160):  CH3 + CH2O <=> HCO + CH4
    fwd_A[162]     = 3320;
    fwd_beta[162]  = 2.8100000000000001;
    fwd_Ea[162]    = 5860;
    prefactor_units[162]  = 1.0000000000000002e-06;
    activation_units[162] = 0.50321666580471969;
    phase_units[162]      = pow(10,-12.000000);
    is_PD[162] = 0;
    nTB[162] = 0;

    // (161):  CH3 + CH3OH <=> CH2OH + CH4
    kiv[163] = {12,20,18,13};
    nuv[163] = {-1,-1,1,1};
    // (161):  CH3 + CH3OH <=> CH2OH + CH4
    fwd_A[163]     = 30000000;
    fwd_beta[163]  = 1.5;
    fwd_Ea[163]    = 9940;
    prefactor_units[163]  = 1.0000000000000002e-06;
    activation_units[163] = 0.50321666580471969;
    phase_units[163]      = pow(10,-12.000000);
    is_PD[163] = 0;
    nTB[163] = 0;

    // (162):  CH3 + CH3OH <=> CH3O + CH4
    kiv[164] = {12,20,19,13};
    nuv[164] = {-1,-1,1,1};
    // (162):  CH3 + CH3OH <=> CH3O + CH4
    fwd_A[164]     = 10000000;
    fwd_beta[164]  = 1.5;
    fwd_Ea[164]    = 9940;
    prefactor_units[164]  = 1.0000000000000002e-06;
    activation_units[164] = 0.50321666580471969;
    phase_units[164]      = pow(10,-12.000000);
    is_PD[164] = 0;
    nTB[164] = 0;

    // (163):  CH3 + C2H4 <=> C2H3 + CH4
    kiv[165] = {12,24,23,13};
    nuv[165] = {-1,-1,1,1};
    // (163):  CH3 + C2H4 <=> C2H3 + CH4
    fwd_A[165]     = 227000;
    fwd_beta[165]  = 2;
    fwd_Ea[165]    = 9200;
    prefactor_units[165]  = 1.0000000000000002e-06;
    activation_units[165] = 0.50321666580471969;
    phase_units[165]      = pow(10,-12.000000);
    is_PD[165] = 0;
    nTB[165] = 0;

    // (164):  CH3 + C2H6 <=> C2H5 + CH4
    kiv[166] = {12,26,25,13};
    nuv[166] = {-1,-1,1,1};
    // (164):  CH3 + C2H6 <=> C2H5 + CH4
    fwd_A[166]     = 6140000;
    fwd_beta[166]  = 1.74;
    fwd_Ea[166]    = 10450;
    prefactor_units[166]  = 1.0000000000000002e-06;
    activation_units[166] = 0.50321666580471969;
    phase_units[166]      = pow(10,-12.000000);
    is_PD[166] = 0;
    nTB[166] = 0;

    // (165):  HCO + H2O <=> H + CO + H2O
    kiv[167] = {16,5,1,14,5};
    nuv[167] = {-1,-1,1,1,1};
    // (165):  HCO + H2O <=> H + CO + H2O
    fwd_A[167]     = 2.244e+18;
    fwd_beta[167]  = -1;
    fwd_Ea[167]    = 17000;
    prefactor_units[167]  = 1.0000000000000002e-06;
    activation_units[167] = 0.50321666580471969;
    phase_units[167]      = pow(10,-12.000000);
    is_PD[167] = 0;
    nTB[167] = 0;

    // (166):  HCO + M <=> H + CO + M
    kiv[26] = {16,1,14};
    nuv[26] = {-1,1,1};
    // (166):  HCO + M <=> H + CO + M
    fwd_A[26]     = 1.87e+17;
    fwd_beta[26]  = -1;
    fwd_Ea[26]    = 17000;
    prefactor_units[26]  = 1.0000000000000002e-06;
    activation_units[26] = 0.50321666580471969;
    phase_units[26]      = pow(10,-6.000000);
    is_PD[26] = 0;
    nTB[26] = 6;
    TB[26] = (amrex::Real *) malloc(6 * sizeof(amrex::Real));
    TBid[26] = (int *) malloc(6 * sizeof(int));
    TBid[26][0] = 0; TB[26][0] = 2; // H2
    TBid[26][1] = 5; TB[26][1] = 0; // H2O
    TBid[26][2] = 13; TB[26][2] = 2; // CH4
    TBid[26][3] = 14; TB[26][3] = 1.5; // CO
    TBid[26][4] = 15; TB[26][4] = 2; // CO2
    TBid[26][5] = 26; TB[26][5] = 3; // C2H6

    // (167):  HCO + O2 <=> HO2 + CO
    kiv[168] = {16,3,6,14};
    nuv[168] = {-1,-1,1,1};
    // (167):  HCO + O2 <=> HO2 + CO
    fwd_A[168]     = 7600000000000;
    fwd_beta[168]  = 0;
    fwd_Ea[168]    = 400;
    prefactor_units[168]  = 1.0000000000000002e-06;
    activation_units[168] = 0.50321666580471969;
    phase_units[168]      = pow(10,-12.000000);
    is_PD[168] = 0;
    nTB[168] = 0;

    // (168):  CH2OH + O2 <=> HO2 + CH2O
    kiv[169] = {18,3,6,17};
    nuv[169] = {-1,-1,1,1};
    // (168):  CH2OH + O2 <=> HO2 + CH2O
    fwd_A[169]     = 18000000000000;
    fwd_beta[169]  = 0;
    fwd_Ea[169]    = 900;
    prefactor_units[169]  = 1.0000000000000002e-06;
    activation_units[169] = 0.50321666580471969;
    phase_units[169]      = pow(10,-12.000000);
    is_PD[169] = 0;
    nTB[169] = 0;

    // (169):  CH3O + O2 <=> HO2 + CH2O
    kiv[170] = {19,3,6,17};
    nuv[170] = {-1,-1,1,1};
    // (169):  CH3O + O2 <=> HO2 + CH2O
    fwd_A[170]     = 4.2799999999999999e-13;
    fwd_beta[170]  = 7.5999999999999996;
    fwd_Ea[170]    = -3530;
    prefactor_units[170]  = 1.0000000000000002e-06;
    activation_units[170] = 0.50321666580471969;
    phase_units[170]      = pow(10,-12.000000);
    is_PD[170] = 0;
    nTB[170] = 0;

    // (170):  C2H + O2 <=> HCO + CO
    kiv[171] = {21,3,16,14};
    nuv[171] = {-1,-1,1,1};
    // (170):  C2H + O2 <=> HCO + CO
    fwd_A[171]     = 50000000000000;
    fwd_beta[171]  = 0;
    fwd_Ea[171]    = 1500;
    prefactor_units[171]  = 1.0000000000000002e-06;
    activation_units[171] = 0.50321666580471969;
    phase_units[171]      = pow(10,-12.000000);
    is_PD[171] = 0;
    nTB[171] = 0;

    // (171):  C2H + H2 <=> H + C2H2
    kiv[172] = {21,0,1,22};
    nuv[172] = {-1,-1,1,1};
    // (171):  C2H + H2 <=> H + C2H2
    fwd_A[172]     = 407000;
    fwd_beta[172]  = 2.3999999999999999;
    fwd_Ea[172]    = 200;
    prefactor_units[172]  = 1.0000000000000002e-06;
    activation_units[172] = 0.50321666580471969;
    phase_units[172]      = pow(10,-12.000000);
    is_PD[172] = 0;
    nTB[172] = 0;

    // (172):  C2H3 + O2 <=> HCO + CH2O
    kiv[173] = {23,3,16,17};
    nuv[173] = {-1,-1,1,1};
    // (172):  C2H3 + O2 <=> HCO + CH2O
    fwd_A[173]     = 3980000000000;
    fwd_beta[173]  = 0;
    fwd_Ea[173]    = -240;
    prefactor_units[173]  = 1.0000000000000002e-06;
    activation_units[173] = 0.50321666580471969;
    phase_units[173]      = pow(10,-12.000000);
    is_PD[173] = 0;
    nTB[173] = 0;

    // (173):  C2H4 (+M) <=> H2 + C2H2 (+M)
    kiv[19] = {24,0,22};
    nuv[19] = {-1,1,1};
    // (173):  C2H4 (+M) <=> H2 + C2H2 (+M)
    fwd_A[19]     = 8000000000000;
    fwd_beta[19]  = 0.44;
    fwd_Ea[19]    = 88770;
    low_A[19]     = 7.0000000000000001e+50;
    low_beta[19]  = -9.3100000000000005;
    low_Ea[19]    = 99860;
    troe_a[19]    = 0.73450000000000004;
    troe_Tsss[19] = 180;
    troe_Ts[19]   = 1035;
    troe_Tss[19]  = 5417;
    troe_len[19]  = 4;
    prefactor_units[19]  = 1;
    activation_units[19] = 0.50321666580471969;
    phase_units[19]      = pow(10,-6.000000);
    is_PD[19] = 1;
    nTB[19] = 7;
    TB[19] = (amrex::Real *) malloc(7 * sizeof(amrex::Real));
    TBid[19] = (int *) malloc(7 * sizeof(int));
    TBid[19][0] = 0; TB[19][0] = 2; // H2
    TBid[19][1] = 5; TB[19][1] = 6; // H2O
    TBid[19][2] = 13; TB[19][2] = 2; // CH4
    TBid[19][3] = 14; TB[19][3] = 1.5; // CO
    TBid[19][4] = 15; TB[19][4] = 2; // CO2
    TBid[19][5] = 26; TB[19][5] = 3; // C2H6
    TBid[19][6] = 31; TB[19][6] = 0.69999999999999996; // AR

    // (174):  C2H5 + O2 <=> HO2 + C2H4
    kiv[174] = {25,3,6,24};
    nuv[174] = {-1,-1,1,1};
    // (174):  C2H5 + O2 <=> HO2 + C2H4
    fwd_A[174]     = 840000000000;
    fwd_beta[174]  = 0;
    fwd_Ea[174]    = 3875;
    prefactor_units[174]  = 1.0000000000000002e-06;
    activation_units[174] = 0.50321666580471969;
    phase_units[174]      = pow(10,-12.000000);
    is_PD[174] = 0;
    nTB[174] = 0;

    // (175):  HCCO + O2 <=> OH + 2.000000 CO
    kiv[175] = {27,3,4,14};
    nuv[175] = {-1,-1,1,2.0};
    // (175):  HCCO + O2 <=> OH + 2.000000 CO
    fwd_A[175]     = 1600000000000;
    fwd_beta[175]  = 0;
    fwd_Ea[175]    = 854;
    prefactor_units[175]  = 1.0000000000000002e-06;
    activation_units[175] = 0.50321666580471969;
    phase_units[175]      = pow(10,-12.000000);
    is_PD[175] = 0;
    nTB[175] = 0;

    // (176):  2.000000 HCCO <=> 2.000000 CO + C2H2
    kiv[176] = {27,14,22};
    nuv[176] = {-2.0,2.0,1};
    // (176):  2.000000 HCCO <=> 2.000000 CO + C2H2
    fwd_A[176]     = 10000000000000;
    fwd_beta[176]  = 0;
    fwd_Ea[176]    = 0;
    prefactor_units[176]  = 1.0000000000000002e-06;
    activation_units[176] = 0.50321666580471969;
    phase_units[176]      = pow(10,-12.000000);
    is_PD[176] = 0;
    nTB[176] = 0;

}


/* Finalizes parameter database */
void CKFINALIZE()
{
  for (int i=0; i<177; ++i) {
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
    for (id = 0; id < kd * 32; ++ id) {
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

    /*C */
    ncf[ 8 * kd + 2 ] = 1; /*C */

    /*CH */
    ncf[ 9 * kd + 2 ] = 1; /*C */
    ncf[ 9 * kd + 1 ] = 1; /*H */

    /*CH2 */
    ncf[ 10 * kd + 2 ] = 1; /*C */
    ncf[ 10 * kd + 1 ] = 2; /*H */

    /*CH2(S) */
    ncf[ 11 * kd + 2 ] = 1; /*C */
    ncf[ 11 * kd + 1 ] = 2; /*H */

    /*CH3 */
    ncf[ 12 * kd + 2 ] = 1; /*C */
    ncf[ 12 * kd + 1 ] = 3; /*H */

    /*CH4 */
    ncf[ 13 * kd + 2 ] = 1; /*C */
    ncf[ 13 * kd + 1 ] = 4; /*H */

    /*CO */
    ncf[ 14 * kd + 2 ] = 1; /*C */
    ncf[ 14 * kd + 0 ] = 1; /*O */

    /*CO2 */
    ncf[ 15 * kd + 2 ] = 1; /*C */
    ncf[ 15 * kd + 0 ] = 2; /*O */

    /*HCO */
    ncf[ 16 * kd + 1 ] = 1; /*H */
    ncf[ 16 * kd + 2 ] = 1; /*C */
    ncf[ 16 * kd + 0 ] = 1; /*O */

    /*CH2O */
    ncf[ 17 * kd + 1 ] = 2; /*H */
    ncf[ 17 * kd + 2 ] = 1; /*C */
    ncf[ 17 * kd + 0 ] = 1; /*O */

    /*CH2OH */
    ncf[ 18 * kd + 2 ] = 1; /*C */
    ncf[ 18 * kd + 1 ] = 3; /*H */
    ncf[ 18 * kd + 0 ] = 1; /*O */

    /*CH3O */
    ncf[ 19 * kd + 2 ] = 1; /*C */
    ncf[ 19 * kd + 1 ] = 3; /*H */
    ncf[ 19 * kd + 0 ] = 1; /*O */

    /*CH3OH */
    ncf[ 20 * kd + 2 ] = 1; /*C */
    ncf[ 20 * kd + 1 ] = 4; /*H */
    ncf[ 20 * kd + 0 ] = 1; /*O */

    /*C2H */
    ncf[ 21 * kd + 2 ] = 2; /*C */
    ncf[ 21 * kd + 1 ] = 1; /*H */

    /*C2H2 */
    ncf[ 22 * kd + 2 ] = 2; /*C */
    ncf[ 22 * kd + 1 ] = 2; /*H */

    /*C2H3 */
    ncf[ 23 * kd + 2 ] = 2; /*C */
    ncf[ 23 * kd + 1 ] = 3; /*H */

    /*C2H4 */
    ncf[ 24 * kd + 2 ] = 2; /*C */
    ncf[ 24 * kd + 1 ] = 4; /*H */

    /*C2H5 */
    ncf[ 25 * kd + 2 ] = 2; /*C */
    ncf[ 25 * kd + 1 ] = 5; /*H */

    /*C2H6 */
    ncf[ 26 * kd + 2 ] = 2; /*C */
    ncf[ 26 * kd + 1 ] = 6; /*H */

    /*HCCO */
    ncf[ 27 * kd + 1 ] = 1; /*H */
    ncf[ 27 * kd + 2 ] = 2; /*C */
    ncf[ 27 * kd + 0 ] = 1; /*O */

    /*CH2CO */
    ncf[ 28 * kd + 2 ] = 2; /*C */
    ncf[ 28 * kd + 1 ] = 2; /*H */
    ncf[ 28 * kd + 0 ] = 1; /*O */

    /*HCCOH */
    ncf[ 29 * kd + 2 ] = 2; /*C */
    ncf[ 29 * kd + 0 ] = 1; /*O */
    ncf[ 29 * kd + 1 ] = 2; /*H */

    /*N2 */
    ncf[ 30 * kd + 3 ] = 2; /*N */

    /*AR */
    ncf[ 31 * kd + 4 ] = 1; /*AR */

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
    kname.resize(32);
    kname[0] = "H2";
    kname[1] = "H";
    kname[2] = "O";
    kname[3] = "O2";
    kname[4] = "OH";
    kname[5] = "H2O";
    kname[6] = "HO2";
    kname[7] = "H2O2";
    kname[8] = "C";
    kname[9] = "CH";
    kname[10] = "CH2";
    kname[11] = "CH2(S)";
    kname[12] = "CH3";
    kname[13] = "CH4";
    kname[14] = "CO";
    kname[15] = "CO2";
    kname[16] = "HCO";
    kname[17] = "CH2O";
    kname[18] = "CH2OH";
    kname[19] = "CH3O";
    kname[20] = "CH3OH";
    kname[21] = "C2H";
    kname[22] = "C2H2";
    kname[23] = "C2H3";
    kname[24] = "C2H4";
    kname[25] = "C2H5";
    kname[26] = "C2H6";
    kname[27] = "HCCO";
    kname[28] = "CH2CO";
    kname[29] = "HCCOH";
    kname[30] = "N2";
    kname[31] = "AR";
}

/*compute the sparsity pattern of the chemistry Jacobian */
void SPARSITY_INFO( int * nJdata, int * consP, int NCELLS)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(1089);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(32);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[1089];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<32; k++) {
                c_d[k] = 1.0/ 32.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<33; k++) {
        for (int l=0; l<33; l++) {
            if(J_h[ 33 * k + l] != 0.0){
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(1089);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(32);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[1089];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<32; k++) {
                c_d[k] = 1.0/ 32.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<33; k++) {
        for (int l=0; l<33; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J_h[ 33 * k + l] != 0.0){
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(1089);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(32);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[1089];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<32; k++) {
                c_d[k] = 1.0/ 32.000000 ;
            }
            aJacobian_precond(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<33; k++) {
        for (int l=0; l<33; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J_h[ 33 * k + l] != 0.0){
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

    amrex::Gpu::DeviceVector<amrex::Real> J_v(1089);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(32);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[1089];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<32; k++) {
                c_d[k] = 1.0/ 32.000000 ;
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
        offset_row = nc * 33;
        offset_col = nc * 33;
        for (int k=0; k<33; k++) {
            for (int l=0; l<33; l++) {
                if(J_h[33*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(1089);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(32);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[1089];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<32; k++) {
                c_d[k] = 1.0/ 32.000000 ;
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
            offset = nc * 33;
            for (int l=0; l<33; l++) {
                for (int k=0; k<33; k++) {
                    if(J_h[33*k + l] != 0.0) {
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
            offset = nc * 33;
            for (int l=0; l<33; l++) {
                for (int k=0; k<33; k++) {
                    if(J_h[33*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(1089);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(32);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[1089];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<32; k++) {
                c_d[k] = 1.0/ 32.000000 ;
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
            offset = nc * 33;
            for (int l=0; l<33; l++) {
                for (int k=0; k<33; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp-1] = l+1 + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(J_h[33*k + l] != 0.0) {
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
            offset = nc * 33;
            for (int l=0; l<33; l++) {
                for (int k=0; k<33; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp] = l + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(J_h[33*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(1089);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(32);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[1089];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<32; k++) {
                c_d[k] = 1.0/ 32.000000 ;
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
    for (int k=0; k<33; k++) {
        for (int l=0; l<33; l++) {
            if (k == l) {
                rowVals[nJdata_tmp] = l; 
                indx[nJdata_tmp] = 33*k + l;
                nJdata_tmp = nJdata_tmp + 1; 
            } else {
                if(J_h[33*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l; 
                    indx[nJdata_tmp] = 33*k + l;
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(1089);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(32);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[1089];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<32; k++) {
                c_d[k] = 1.0/ 32.000000 ;
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
        for (int l=0; l<33; l++) {
            for (int k=0; k<33; k++) {
                if (k == l) {
                    colVals[nJdata_tmp-1] = l+1; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J_h[33*k + l] != 0.0) {
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
        for (int l=0; l<33; l++) {
            for (int k=0; k<33; k++) {
                if (k == l) {
                    colVals[nJdata_tmp] = l; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J_h[33*k + l] != 0.0) {
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

