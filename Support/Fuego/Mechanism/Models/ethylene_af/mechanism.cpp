#ifndef MECHANISM_CPP
#define MECHANISM_CPP

#include "chemistry_file.H"

/* PURE CPU stuff  */
#ifndef AMREX_USE_GPU
namespace thermo
{
    amrex::Real fwd_A[172], fwd_beta[172], fwd_Ea[172];
    amrex::Real low_A[172], low_beta[172], low_Ea[172];
    amrex::Real rev_A[172], rev_beta[172], rev_Ea[172];
    amrex::Real troe_a[172],troe_Ts[172], troe_Tss[172], troe_Tsss[172];
    amrex::Real sri_a[172], sri_b[172], sri_c[172], sri_d[172], sri_e[172];
    amrex::Real activation_units[172], prefactor_units[172], phase_units[172];
    int is_PD[172], troe_len[172], sri_len[172], nTB[172], *TBid[172];
    amrex::Real *TB[172];
    std::vector<std::vector<amrex::Real>> kiv(172); 
    std::vector<std::vector<amrex::Real>> nuv(172); 
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
        if (*i > 172) {
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
    amrex::Real imw[29];

    get_imw(imw);

    for (int i=0; i<(*np); i++) {
        YOW[i] = 0.0;
    }

    for (int n=0; n<29; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] = y[n*(*np)+i] * imw[n];
            YOW[i] += x[n*(*np)+i];
        }
    }

    for (int i=0; i<(*np); i++) {
        YOW[i] = 1.0/YOW[i];
    }

    for (int n=0; n<29; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] *=  YOW[i];
        }
    }
}

/*Returns enthalpy in mass units (Eq 27.) */
void VCKHMS(int *  np, amrex::Real *  T,  amrex::Real *  hms)
{
    amrex::Real tc[5], h[29];
    amrex::Real imw[29];

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
    }

    for (int n=0; n<29; n++) {
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
    amrex::Real c[29*(*np)]; /*temporary storage */
    amrex::Real imw[29];

    get_imw(imw);

    /*See Eq 8 with an extra 1e6 so c goes to SI */
    for (int n=0; n<29; n++) {
        for (int i=0; i<(*np); i++) {
            c[n*(*np)+i] = 1.0e6 * rho[i] * y[n*(*np)+i] * imw[n];
        }
    }

    /*call productionRate */
    vproductionRate(*np, wdot, c, T);

    /*convert to chemkin units */
    for (int i=0; i<29*(*np); i++) {
        wdot[i] *= 1.0e-6;
    }
}


/*Compute P = rhoRT/W(y) */
void VCKPY(int *  np, amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  P)
{
    amrex::Real YOW[*np];
    amrex::Real imw[29];

    get_imw(imw);

    for (int i=0; i<(*np); i++) {
        YOW[i] = 0.0;
    }

    for (int n=0; n<29; n++) {
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
    amrex::Real k_f_s[172*npt], Kc_s[172*npt], mixture[npt], g_RT[29*npt];
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

    for (int n=0; n<29; n++) {
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
    vcomp_wdot_151_172(npt, wdot, mixture, sc, k_f_s, Kc_s, tc, invT, T);
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
    }
}

void vcomp_gibbs(int npt, amrex::Real *  g_RT, amrex::Real *  tc)
{
    /*compute the Gibbs free energy */
    for (int i=0; i<npt; i++) {
        amrex::Real tg[5], g[29];
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
    }
}

void vcomp_Kc(int npt, amrex::Real *  Kc_s, amrex::Real *  g_RT, amrex::Real *  invT)
{
    for (int i=0; i<npt; i++) {
        /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
        amrex::Real refC = (101325. / 8.31451) * invT[i];
        amrex::Real refCinv = 1.0 / refC;

        Kc_s[0*npt+i] = refCinv * exp((2.000000 * g_RT[4*npt+i]) - (g_RT[6*npt+i]));
        Kc_s[1*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[5*npt+i]) - (g_RT[8*npt+i]));
        Kc_s[2*npt+i] = refCinv * exp((g_RT[2*npt+i] + g_RT[19*npt+i]) - (g_RT[10*npt+i]));
        Kc_s[3*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[21*npt+i]) - (g_RT[10*npt+i]));
        Kc_s[4*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[10*npt+i]) - (g_RT[13*npt+i]));
        Kc_s[5*npt+i] = refCinv * exp((g_RT[9*npt+i] + g_RT[21*npt+i]) - (g_RT[16*npt+i]));
        Kc_s[6*npt+i] = refCinv * exp((g_RT[2*npt+i] + g_RT[9*npt+i]) - (g_RT[11*npt+i]));
        Kc_s[7*npt+i] = refCinv * exp((g_RT[9*npt+i] + g_RT[19*npt+i]) - (g_RT[25*npt+i]));
        Kc_s[8*npt+i] = refCinv * exp((g_RT[3*npt+i] + g_RT[9*npt+i]) - (g_RT[12*npt+i]));
        Kc_s[9*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[20*npt+i]) - (g_RT[11*npt+i]));
        Kc_s[10*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[14*npt+i]) - (g_RT[23*npt+i]));
        Kc_s[11*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[23*npt+i]) - (g_RT[15*npt+i]));
        Kc_s[12*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[15*npt+i]) - (g_RT[24*npt+i]));
        Kc_s[13*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[24*npt+i]) - (g_RT[17*npt+i]));
        Kc_s[14*npt+i] = refC * exp((g_RT[17*npt+i]) - (2.000000 * g_RT[10*npt+i]));
        Kc_s[15*npt+i] = refCinv * exp((2.000000 * g_RT[1*npt+i]) - (g_RT[2*npt+i]));
        Kc_s[16*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[3*npt+i]) - (g_RT[4*npt+i]));
        Kc_s[17*npt+i] = refCinv * exp((2.000000 * g_RT[3*npt+i]) - (g_RT[5*npt+i]));
        Kc_s[18*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[4*npt+i]) - (g_RT[7*npt+i]));
        Kc_s[19*npt+i] = refC * exp((g_RT[20*npt+i]) - (g_RT[1*npt+i] + g_RT[9*npt+i]));
        Kc_s[20*npt+i] = refCinv * exp((2.000000 * g_RT[1*npt+i] + g_RT[2*npt+i]) - (2.000000 * g_RT[2*npt+i]));
        Kc_s[21*npt+i] = exp((g_RT[2*npt+i] + g_RT[3*npt+i]) - (g_RT[1*npt+i] + g_RT[4*npt+i]));
        Kc_s[22*npt+i] = exp((2.000000 * g_RT[4*npt+i]) - (g_RT[3*npt+i] + g_RT[7*npt+i]));
        Kc_s[23*npt+i] = exp((g_RT[2*npt+i] + g_RT[4*npt+i]) - (g_RT[1*npt+i] + g_RT[7*npt+i]));
        Kc_s[24*npt+i] = refCinv * exp((2.000000 * g_RT[1*npt+i] + g_RT[7*npt+i]) - (g_RT[2*npt+i] + g_RT[7*npt+i]));
        Kc_s[25*npt+i] = exp((g_RT[1*npt+i] + g_RT[5*npt+i]) - (g_RT[3*npt+i] + g_RT[4*npt+i]));
        Kc_s[26*npt+i] = exp((g_RT[2*npt+i] + g_RT[5*npt+i]) - (g_RT[1*npt+i] + g_RT[8*npt+i]));
        Kc_s[27*npt+i] = exp((g_RT[4*npt+i] + g_RT[8*npt+i]) - (g_RT[5*npt+i] + g_RT[7*npt+i]));
        Kc_s[28*npt+i] = exp((g_RT[1*npt+i] + g_RT[8*npt+i]) - (2.000000 * g_RT[4*npt+i]));
        Kc_s[29*npt+i] = exp((g_RT[3*npt+i] + g_RT[8*npt+i]) - (g_RT[4*npt+i] + g_RT[5*npt+i]));
        Kc_s[30*npt+i] = exp((g_RT[1*npt+i] + g_RT[8*npt+i]) - (g_RT[3*npt+i] + g_RT[7*npt+i]));
        Kc_s[31*npt+i] = exp((g_RT[4*npt+i] + g_RT[8*npt+i]) - (g_RT[5*npt+i] + g_RT[7*npt+i]));
        Kc_s[32*npt+i] = exp((g_RT[3*npt+i] + g_RT[6*npt+i]) - (g_RT[4*npt+i] + g_RT[8*npt+i]));
        Kc_s[33*npt+i] = exp((g_RT[1*npt+i] + g_RT[6*npt+i]) - (g_RT[4*npt+i] + g_RT[7*npt+i]));
        Kc_s[34*npt+i] = exp((g_RT[1*npt+i] + g_RT[6*npt+i]) - (g_RT[2*npt+i] + g_RT[8*npt+i]));
        Kc_s[35*npt+i] = exp((g_RT[4*npt+i] + g_RT[6*npt+i]) - (g_RT[7*npt+i] + g_RT[8*npt+i]));
        Kc_s[36*npt+i] = exp((g_RT[4*npt+i] + g_RT[6*npt+i]) - (g_RT[7*npt+i] + g_RT[8*npt+i]));
        Kc_s[37*npt+i] = exp((g_RT[5*npt+i] + g_RT[18*npt+i]) - (g_RT[3*npt+i] + g_RT[9*npt+i]));
        Kc_s[38*npt+i] = exp((g_RT[4*npt+i] + g_RT[18*npt+i]) - (g_RT[1*npt+i] + g_RT[9*npt+i]));
        Kc_s[39*npt+i] = exp((g_RT[4*npt+i] + g_RT[19*npt+i]) - (g_RT[1*npt+i] + g_RT[20*npt+i]));
        Kc_s[40*npt+i] = exp((g_RT[2*npt+i] + g_RT[19*npt+i]) - (g_RT[1*npt+i] + g_RT[21*npt+i]));
        Kc_s[41*npt+i] = exp((g_RT[3*npt+i] + g_RT[19*npt+i]) - (g_RT[1*npt+i] + g_RT[9*npt+i]));
        Kc_s[42*npt+i] = exp((g_RT[5*npt+i] + g_RT[19*npt+i]) - (g_RT[3*npt+i] + g_RT[20*npt+i]));
        Kc_s[43*npt+i] = exp((g_RT[1*npt+i] + g_RT[19*npt+i]) - (g_RT[2*npt+i] + g_RT[18*npt+i]));
        Kc_s[44*npt+i] = exp((g_RT[7*npt+i] + g_RT[19*npt+i]) - (g_RT[1*npt+i] + g_RT[11*npt+i]));
        Kc_s[45*npt+i] = refC * exp((g_RT[5*npt+i] + g_RT[21*npt+i]) - (g_RT[1*npt+i] + g_RT[4*npt+i] + g_RT[9*npt+i]));
        Kc_s[46*npt+i] = exp((g_RT[5*npt+i] + g_RT[21*npt+i]) - (g_RT[3*npt+i] + g_RT[11*npt+i]));
        Kc_s[47*npt+i] = exp((g_RT[4*npt+i] + g_RT[21*npt+i]) - (g_RT[1*npt+i] + g_RT[11*npt+i]));
        Kc_s[48*npt+i] = exp((g_RT[8*npt+i] + g_RT[21*npt+i]) - (g_RT[4*npt+i] + g_RT[11*npt+i]));
        Kc_s[49*npt+i] = refC * exp((g_RT[5*npt+i] + g_RT[21*npt+i]) - (2.000000 * g_RT[1*npt+i] + g_RT[12*npt+i]));
        Kc_s[50*npt+i] = exp((g_RT[4*npt+i] + g_RT[21*npt+i]) - (g_RT[7*npt+i] + g_RT[19*npt+i]));
        Kc_s[51*npt+i] = exp((g_RT[3*npt+i] + g_RT[21*npt+i]) - (g_RT[1*npt+i] + g_RT[20*npt+i]));
        Kc_s[52*npt+i] = exp((g_RT[2*npt+i] + g_RT[21*npt+i]) - (g_RT[1*npt+i] + g_RT[10*npt+i]));
        Kc_s[53*npt+i] = exp((g_RT[7*npt+i] + g_RT[22*npt+i]) - (g_RT[7*npt+i] + g_RT[21*npt+i]));
        Kc_s[54*npt+i] = exp((g_RT[1*npt+i] + g_RT[22*npt+i]) - (g_RT[2*npt+i] + g_RT[19*npt+i]));
        Kc_s[55*npt+i] = refC * exp((g_RT[5*npt+i] + g_RT[22*npt+i]) - (g_RT[1*npt+i] + g_RT[4*npt+i] + g_RT[9*npt+i]));
        Kc_s[56*npt+i] = exp((g_RT[3*npt+i] + g_RT[22*npt+i]) - (g_RT[2*npt+i] + g_RT[9*npt+i]));
        Kc_s[57*npt+i] = exp((g_RT[5*npt+i] + g_RT[22*npt+i]) - (g_RT[7*npt+i] + g_RT[9*npt+i]));
        Kc_s[58*npt+i] = exp((g_RT[2*npt+i] + g_RT[22*npt+i]) - (g_RT[1*npt+i] + g_RT[10*npt+i]));
        Kc_s[59*npt+i] = exp((g_RT[3*npt+i] + g_RT[22*npt+i]) - (g_RT[1*npt+i] + g_RT[20*npt+i]));
        Kc_s[60*npt+i] = exp((g_RT[7*npt+i] + g_RT[22*npt+i]) - (g_RT[2*npt+i] + g_RT[11*npt+i]));
        Kc_s[61*npt+i] = exp((g_RT[4*npt+i] + g_RT[22*npt+i]) - (g_RT[1*npt+i] + g_RT[11*npt+i]));
        Kc_s[62*npt+i] = exp((g_RT[4*npt+i] + g_RT[10*npt+i]) - (g_RT[2*npt+i] + g_RT[11*npt+i]));
        Kc_s[63*npt+i] = exp((g_RT[6*npt+i] + g_RT[10*npt+i]) - (g_RT[8*npt+i] + g_RT[13*npt+i]));
        Kc_s[64*npt+i] = exp((g_RT[5*npt+i] + g_RT[10*npt+i]) - (g_RT[4*npt+i] + g_RT[11*npt+i]));
        Kc_s[65*npt+i] = exp((g_RT[10*npt+i] + g_RT[19*npt+i]) - (g_RT[1*npt+i] + g_RT[23*npt+i]));
        Kc_s[66*npt+i] = exp((g_RT[3*npt+i] + g_RT[10*npt+i]) - (g_RT[1*npt+i] + g_RT[11*npt+i]));
        Kc_s[67*npt+i] = exp((g_RT[10*npt+i] + g_RT[18*npt+i]) - (g_RT[1*npt+i] + g_RT[14*npt+i]));
        Kc_s[68*npt+i] = exp((g_RT[4*npt+i] + g_RT[10*npt+i]) - (g_RT[7*npt+i] + g_RT[21*npt+i]));
        Kc_s[69*npt+i] = exp((g_RT[10*npt+i] + g_RT[22*npt+i]) - (g_RT[1*npt+i] + g_RT[15*npt+i]));
        Kc_s[70*npt+i] = exp((g_RT[4*npt+i] + g_RT[10*npt+i]) - (g_RT[7*npt+i] + g_RT[22*npt+i]));
        Kc_s[71*npt+i] = exp((2.000000 * g_RT[10*npt+i]) - (g_RT[1*npt+i] + g_RT[24*npt+i]));
        Kc_s[72*npt+i] = exp((g_RT[8*npt+i] + g_RT[10*npt+i]) - (g_RT[5*npt+i] + g_RT[13*npt+i]));
        Kc_s[73*npt+i] = exp((g_RT[10*npt+i] + g_RT[21*npt+i]) - (g_RT[1*npt+i] + g_RT[15*npt+i]));
        Kc_s[74*npt+i] = refC * exp((g_RT[3*npt+i] + g_RT[10*npt+i]) - (g_RT[1*npt+i] + g_RT[2*npt+i] + g_RT[9*npt+i]));
        Kc_s[75*npt+i] = exp((g_RT[13*npt+i] + g_RT[19*npt+i]) - (g_RT[1*npt+i] + g_RT[15*npt+i]));
        Kc_s[76*npt+i] = exp((g_RT[13*npt+i] + g_RT[22*npt+i]) - (2.000000 * g_RT[10*npt+i]));
        Kc_s[77*npt+i] = exp((g_RT[3*npt+i] + g_RT[13*npt+i]) - (g_RT[4*npt+i] + g_RT[10*npt+i]));
        Kc_s[78*npt+i] = exp((g_RT[4*npt+i] + g_RT[13*npt+i]) - (g_RT[7*npt+i] + g_RT[10*npt+i]));
        Kc_s[79*npt+i] = exp((g_RT[13*npt+i] + g_RT[21*npt+i]) - (2.000000 * g_RT[10*npt+i]));
        Kc_s[80*npt+i] = exp((g_RT[1*npt+i] + g_RT[13*npt+i]) - (g_RT[2*npt+i] + g_RT[10*npt+i]));
        Kc_s[81*npt+i] = exp((g_RT[9*npt+i] + g_RT[22*npt+i]) - (g_RT[9*npt+i] + g_RT[21*npt+i]));
        Kc_s[82*npt+i] = exp((g_RT[5*npt+i] + g_RT[9*npt+i]) - (g_RT[3*npt+i] + g_RT[12*npt+i]));
        Kc_s[83*npt+i] = exp((g_RT[4*npt+i] + g_RT[9*npt+i]) - (g_RT[1*npt+i] + g_RT[12*npt+i]));
        Kc_s[84*npt+i] = exp((g_RT[4*npt+i] + g_RT[9*npt+i]) - (g_RT[1*npt+i] + g_RT[12*npt+i]));
        Kc_s[85*npt+i] = exp((g_RT[8*npt+i] + g_RT[9*npt+i]) - (g_RT[4*npt+i] + g_RT[12*npt+i]));
        Kc_s[86*npt+i] = exp((g_RT[1*npt+i] + g_RT[20*npt+i]) - (g_RT[2*npt+i] + g_RT[9*npt+i]));
        Kc_s[87*npt+i] = refCinv * exp((g_RT[10*npt+i] + g_RT[20*npt+i]) - (g_RT[26*npt+i]));
        Kc_s[88*npt+i] = refC * exp((g_RT[7*npt+i] + g_RT[20*npt+i]) - (g_RT[1*npt+i] + g_RT[7*npt+i] + g_RT[9*npt+i]));
        Kc_s[89*npt+i] = exp((g_RT[3*npt+i] + g_RT[20*npt+i]) - (g_RT[4*npt+i] + g_RT[9*npt+i]));
        Kc_s[90*npt+i] = exp((g_RT[4*npt+i] + g_RT[20*npt+i]) - (g_RT[7*npt+i] + g_RT[9*npt+i]));
        Kc_s[91*npt+i] = exp((g_RT[10*npt+i] + g_RT[20*npt+i]) - (g_RT[9*npt+i] + g_RT[13*npt+i]));
        Kc_s[92*npt+i] = exp((g_RT[3*npt+i] + g_RT[20*npt+i]) - (g_RT[1*npt+i] + g_RT[12*npt+i]));
        Kc_s[93*npt+i] = exp((g_RT[5*npt+i] + g_RT[20*npt+i]) - (g_RT[8*npt+i] + g_RT[9*npt+i]));
        Kc_s[94*npt+i] = exp((g_RT[1*npt+i] + g_RT[11*npt+i]) - (g_RT[2*npt+i] + g_RT[20*npt+i]));
        Kc_s[95*npt+i] = exp((g_RT[3*npt+i] + g_RT[11*npt+i]) - (g_RT[4*npt+i] + g_RT[20*npt+i]));
        Kc_s[96*npt+i] = exp((g_RT[10*npt+i] + g_RT[11*npt+i]) - (g_RT[13*npt+i] + g_RT[20*npt+i]));
        Kc_s[97*npt+i] = exp((g_RT[4*npt+i] + g_RT[11*npt+i]) - (g_RT[7*npt+i] + g_RT[20*npt+i]));
        Kc_s[98*npt+i] = exp((g_RT[11*npt+i] + g_RT[19*npt+i]) - (g_RT[1*npt+i] + g_RT[16*npt+i]));
        Kc_s[99*npt+i] = exp((g_RT[5*npt+i] + g_RT[11*npt+i]) - (g_RT[8*npt+i] + g_RT[20*npt+i]));
        Kc_s[100*npt+i] = exp((g_RT[8*npt+i] + g_RT[11*npt+i]) - (g_RT[6*npt+i] + g_RT[20*npt+i]));
        Kc_s[101*npt+i] = refCinv * exp((2.000000 * g_RT[1*npt+i] + g_RT[12*npt+i]) - (g_RT[2*npt+i] + g_RT[12*npt+i]));
        Kc_s[102*npt+i] = exp((g_RT[12*npt+i] + g_RT[22*npt+i]) - (g_RT[12*npt+i] + g_RT[21*npt+i]));
        Kc_s[103*npt+i] = exp((g_RT[12*npt+i] + g_RT[22*npt+i]) - (g_RT[9*npt+i] + g_RT[11*npt+i]));
        Kc_s[104*npt+i] = exp((g_RT[12*npt+i] + g_RT[19*npt+i]) - (g_RT[9*npt+i] + g_RT[20*npt+i]));
        Kc_s[105*npt+i] = exp((g_RT[3*npt+i] + g_RT[14*npt+i]) - (g_RT[9*npt+i] + g_RT[21*npt+i]));
        Kc_s[106*npt+i] = exp((g_RT[4*npt+i] + g_RT[14*npt+i]) - (g_RT[9*npt+i] + g_RT[10*npt+i]));
        Kc_s[107*npt+i] = exp((g_RT[4*npt+i] + g_RT[14*npt+i]) - (g_RT[1*npt+i] + g_RT[16*npt+i]));
        Kc_s[108*npt+i] = exp((g_RT[3*npt+i] + g_RT[14*npt+i]) - (g_RT[1*npt+i] + g_RT[25*npt+i]));
        Kc_s[109*npt+i] = exp((g_RT[4*npt+i] + g_RT[23*npt+i]) - (g_RT[7*npt+i] + g_RT[14*npt+i]));
        Kc_s[110*npt+i] = exp((g_RT[5*npt+i] + g_RT[23*npt+i]) - (g_RT[3*npt+i] + g_RT[27*npt+i]));
        Kc_s[111*npt+i] = refCinv * exp((g_RT[3*npt+i] + g_RT[23*npt+i]) - (g_RT[27*npt+i]));
        Kc_s[112*npt+i] = exp((g_RT[1*npt+i] + g_RT[23*npt+i]) - (g_RT[2*npt+i] + g_RT[14*npt+i]));
        Kc_s[113*npt+i] = exp((g_RT[10*npt+i] + g_RT[23*npt+i]) - (g_RT[13*npt+i] + g_RT[14*npt+i]));
        Kc_s[114*npt+i] = exp((g_RT[5*npt+i] + g_RT[23*npt+i]) - (g_RT[11*npt+i] + g_RT[20*npt+i]));
        Kc_s[115*npt+i] = exp((g_RT[6*npt+i] + g_RT[23*npt+i]) - (g_RT[8*npt+i] + g_RT[15*npt+i]));
        Kc_s[116*npt+i] = exp((g_RT[5*npt+i] + g_RT[23*npt+i]) - (g_RT[8*npt+i] + g_RT[14*npt+i]));
        Kc_s[117*npt+i] = exp((g_RT[10*npt+i] + g_RT[15*npt+i]) - (g_RT[13*npt+i] + g_RT[23*npt+i]));
        Kc_s[118*npt+i] = refC * exp((g_RT[5*npt+i] + g_RT[15*npt+i]) - (g_RT[1*npt+i] + g_RT[10*npt+i] + g_RT[12*npt+i]));
        Kc_s[119*npt+i] = exp((g_RT[4*npt+i] + g_RT[15*npt+i]) - (g_RT[7*npt+i] + g_RT[23*npt+i]));
        Kc_s[120*npt+i] = refCinv * exp((g_RT[4*npt+i] + g_RT[15*npt+i]) - (g_RT[28*npt+i]));
        Kc_s[121*npt+i] = exp((g_RT[3*npt+i] + g_RT[15*npt+i]) - (g_RT[1*npt+i] + g_RT[27*npt+i]));
        Kc_s[122*npt+i] = exp((g_RT[3*npt+i] + g_RT[15*npt+i]) - (g_RT[10*npt+i] + g_RT[20*npt+i]));
        Kc_s[123*npt+i] = exp((g_RT[5*npt+i] + g_RT[15*npt+i]) - (g_RT[8*npt+i] + g_RT[23*npt+i]));
        Kc_s[124*npt+i] = exp((g_RT[1*npt+i] + g_RT[15*npt+i]) - (g_RT[2*npt+i] + g_RT[23*npt+i]));
        Kc_s[125*npt+i] = exp((g_RT[3*npt+i] + g_RT[15*npt+i]) - (g_RT[11*npt+i] + g_RT[21*npt+i]));
        Kc_s[126*npt+i] = exp((g_RT[8*npt+i] + g_RT[24*npt+i]) - (g_RT[6*npt+i] + g_RT[15*npt+i]));
        Kc_s[127*npt+i] = exp((g_RT[8*npt+i] + g_RT[24*npt+i]) - (g_RT[4*npt+i] + g_RT[28*npt+i]));
        Kc_s[128*npt+i] = refCinv * exp((g_RT[3*npt+i] + g_RT[24*npt+i]) - (g_RT[28*npt+i]));
        Kc_s[129*npt+i] = exp((g_RT[1*npt+i] + g_RT[24*npt+i]) - (g_RT[2*npt+i] + g_RT[15*npt+i]));
        Kc_s[130*npt+i] = exp((g_RT[5*npt+i] + g_RT[24*npt+i]) - (g_RT[8*npt+i] + g_RT[15*npt+i]));
        Kc_s[131*npt+i] = exp((g_RT[8*npt+i] + g_RT[24*npt+i]) - (g_RT[5*npt+i] + g_RT[17*npt+i]));
        Kc_s[132*npt+i] = exp((g_RT[10*npt+i] + g_RT[24*npt+i]) - (g_RT[13*npt+i] + g_RT[15*npt+i]));
        Kc_s[133*npt+i] = exp((g_RT[17*npt+i] + g_RT[22*npt+i]) - (g_RT[10*npt+i] + g_RT[24*npt+i]));
        Kc_s[134*npt+i] = exp((g_RT[10*npt+i] + g_RT[17*npt+i]) - (g_RT[13*npt+i] + g_RT[24*npt+i]));
        Kc_s[135*npt+i] = exp((g_RT[3*npt+i] + g_RT[17*npt+i]) - (g_RT[4*npt+i] + g_RT[24*npt+i]));
        Kc_s[136*npt+i] = exp((g_RT[8*npt+i] + g_RT[17*npt+i]) - (g_RT[6*npt+i] + g_RT[24*npt+i]));
        Kc_s[137*npt+i] = exp((g_RT[1*npt+i] + g_RT[17*npt+i]) - (g_RT[2*npt+i] + g_RT[24*npt+i]));
        Kc_s[138*npt+i] = exp((g_RT[4*npt+i] + g_RT[17*npt+i]) - (g_RT[7*npt+i] + g_RT[24*npt+i]));
        Kc_s[139*npt+i] = refC * exp((g_RT[5*npt+i] + g_RT[25*npt+i]) - (g_RT[4*npt+i] + 2.000000 * g_RT[9*npt+i]));
        Kc_s[140*npt+i] = refC * exp((g_RT[3*npt+i] + g_RT[25*npt+i]) - (g_RT[1*npt+i] + 2.000000 * g_RT[9*npt+i]));
        Kc_s[141*npt+i] = exp((g_RT[10*npt+i] + g_RT[25*npt+i]) - (g_RT[9*npt+i] + g_RT[15*npt+i]));
        Kc_s[142*npt+i] = exp((g_RT[1*npt+i] + g_RT[25*npt+i]) - (g_RT[9*npt+i] + g_RT[22*npt+i]));
        Kc_s[143*npt+i] = exp((g_RT[1*npt+i] + g_RT[16*npt+i]) - (g_RT[9*npt+i] + g_RT[10*npt+i]));
        Kc_s[144*npt+i] = exp((g_RT[16*npt+i] + g_RT[21*npt+i]) - (g_RT[9*npt+i] + g_RT[15*npt+i]));
        Kc_s[145*npt+i] = exp((g_RT[3*npt+i] + g_RT[16*npt+i]) - (g_RT[4*npt+i] + g_RT[25*npt+i]));
        Kc_s[146*npt+i] = exp((g_RT[10*npt+i] + g_RT[16*npt+i]) - (g_RT[13*npt+i] + g_RT[25*npt+i]));
        Kc_s[147*npt+i] = exp((g_RT[3*npt+i] + g_RT[16*npt+i]) - (g_RT[12*npt+i] + g_RT[21*npt+i]));
        Kc_s[148*npt+i] = exp((g_RT[10*npt+i] + g_RT[16*npt+i]) - (g_RT[9*npt+i] + g_RT[24*npt+i]));
        Kc_s[149*npt+i] = exp((g_RT[4*npt+i] + g_RT[16*npt+i]) - (g_RT[7*npt+i] + g_RT[25*npt+i]));
        Kc_s[150*npt+i] = exp((g_RT[1*npt+i] + g_RT[16*npt+i]) - (g_RT[2*npt+i] + g_RT[25*npt+i]));
        Kc_s[151*npt+i] = exp((g_RT[16*npt+i] + g_RT[21*npt+i]) - (g_RT[10*npt+i] + g_RT[25*npt+i]));
        Kc_s[152*npt+i] = exp((g_RT[3*npt+i] + g_RT[27*npt+i]) - (g_RT[11*npt+i] + g_RT[20*npt+i]));
        Kc_s[153*npt+i] = refC * exp((g_RT[27*npt+i]) - (g_RT[1*npt+i] + g_RT[16*npt+i]));
        Kc_s[154*npt+i] = exp((g_RT[4*npt+i] + g_RT[27*npt+i]) - (g_RT[7*npt+i] + g_RT[16*npt+i]));
        Kc_s[155*npt+i] = exp((g_RT[1*npt+i] + g_RT[27*npt+i]) - (g_RT[2*npt+i] + g_RT[16*npt+i]));
        Kc_s[156*npt+i] = refC * exp((g_RT[5*npt+i] + g_RT[27*npt+i]) - (g_RT[4*npt+i] + g_RT[9*npt+i] + g_RT[11*npt+i]));
        Kc_s[157*npt+i] = refC * exp((g_RT[27*npt+i]) - (g_RT[9*npt+i] + g_RT[10*npt+i]));
        Kc_s[158*npt+i] = refC * exp((g_RT[5*npt+i] + g_RT[27*npt+i]) - (g_RT[4*npt+i] + 2.000000 * g_RT[20*npt+i]));
        Kc_s[159*npt+i] = exp((g_RT[1*npt+i] + g_RT[27*npt+i]) - (g_RT[10*npt+i] + g_RT[20*npt+i]));
        Kc_s[160*npt+i] = refC * exp((g_RT[3*npt+i] + g_RT[26*npt+i]) - (g_RT[4*npt+i] + g_RT[9*npt+i] + g_RT[10*npt+i]));
        Kc_s[161*npt+i] = refC * exp((g_RT[5*npt+i] + g_RT[26*npt+i]) - (g_RT[8*npt+i] + g_RT[9*npt+i] + g_RT[10*npt+i]));
        Kc_s[162*npt+i] = refC * exp((g_RT[4*npt+i] + g_RT[26*npt+i]) - (g_RT[7*npt+i] + g_RT[9*npt+i] + g_RT[10*npt+i]));
        Kc_s[163*npt+i] = exp((g_RT[1*npt+i] + g_RT[26*npt+i]) - (g_RT[2*npt+i] + g_RT[27*npt+i]));
        Kc_s[164*npt+i] = refC * exp((g_RT[1*npt+i] + g_RT[26*npt+i]) - (g_RT[2*npt+i] + g_RT[9*npt+i] + g_RT[10*npt+i]));
        Kc_s[165*npt+i] = exp((g_RT[3*npt+i] + g_RT[26*npt+i]) - (g_RT[4*npt+i] + g_RT[27*npt+i]));
        Kc_s[166*npt+i] = refC * exp((g_RT[10*npt+i] + g_RT[26*npt+i]) - (g_RT[9*npt+i] + g_RT[10*npt+i] + g_RT[13*npt+i]));
        Kc_s[167*npt+i] = refC * exp((g_RT[8*npt+i] + g_RT[26*npt+i]) - (g_RT[6*npt+i] + g_RT[9*npt+i] + g_RT[10*npt+i]));
        Kc_s[168*npt+i] = refC * exp((g_RT[28*npt+i]) - (g_RT[10*npt+i] + g_RT[11*npt+i]));
        Kc_s[169*npt+i] = refC * exp((g_RT[28*npt+i]) - (g_RT[1*npt+i] + g_RT[26*npt+i]));
        Kc_s[170*npt+i] = exp((g_RT[5*npt+i] + g_RT[28*npt+i]) - (g_RT[8*npt+i] + g_RT[26*npt+i]));
        Kc_s[171*npt+i] = exp((g_RT[0*npt+i] + g_RT[22*npt+i]) - (g_RT[0*npt+i] + g_RT[21*npt+i]));
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

        /*reaction 1: 2.000000 OH (+M) <=> H2O2 (+M) */
        phi_f = pow(sc[4*npt+i], 2.000000);
        alpha = mixture[i] + (TB[0][0] - 1)*sc[2*npt+i] + (TB[0][1] - 1)*sc[7*npt+i] + (TB[0][2] - 1)*sc[18*npt+i] + (TB[0][3] - 1)*sc[9*npt+i] + (TB[0][4] - 1)*sc[19*npt+i] + (TB[0][5] - 1)*sc[20*npt+i] + (TB[0][6] - 1)*sc[21*npt+i] + (TB[0][7] - 1)*sc[12*npt+i] + (TB[0][8] - 1)*sc[22*npt+i] + (TB[0][9] - 1)*sc[13*npt+i] + (TB[0][10] - 1)*sc[23*npt+i] + (TB[0][11] - 1)*sc[24*npt+i] + (TB[0][12] - 1)*sc[25*npt+i] + (TB[0][13] - 1)*sc[26*npt+i] + (TB[0][14] - 1)*sc[27*npt+i] + (TB[0][15] - 1)*sc[28*npt+i] + (TB[0][16] - 1)*sc[17*npt+i];
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
        wdot[4*npt+i] -= 2.000000 * qdot;
        wdot[6*npt+i] += qdot;

        /*reaction 2: H + O2 (+M) <=> HO2 (+M) */
        phi_f = sc[1*npt+i]*sc[5*npt+i];
        alpha = mixture[i] + (TB[1][0] - 1)*sc[2*npt+i] + (TB[1][1] - 1)*sc[5*npt+i] + (TB[1][2] - 1)*sc[7*npt+i] + (TB[1][3] - 1)*sc[18*npt+i] + (TB[1][4] - 1)*sc[9*npt+i] + (TB[1][5] - 1)*sc[19*npt+i] + (TB[1][6] - 1)*sc[20*npt+i] + (TB[1][7] - 1)*sc[21*npt+i] + (TB[1][8] - 1)*sc[12*npt+i] + (TB[1][9] - 1)*sc[22*npt+i] + (TB[1][10] - 1)*sc[23*npt+i] + (TB[1][11] - 1)*sc[24*npt+i] + (TB[1][12] - 1)*sc[25*npt+i] + (TB[1][13] - 1)*sc[26*npt+i] + (TB[1][14] - 1)*sc[27*npt+i] + (TB[1][15] - 1)*sc[28*npt+i];
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
        phi_r = sc[8*npt+i];
        Kc = Kc_s[1*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[5*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;

        /*reaction 3: CH + H2 (+M) <=> CH3 (+M) */
        phi_f = sc[2*npt+i]*sc[19*npt+i];
        alpha = mixture[i] + (TB[2][0] - 1)*sc[2*npt+i] + (TB[2][1] - 1)*sc[7*npt+i] + (TB[2][2] - 1)*sc[18*npt+i] + (TB[2][3] - 1)*sc[9*npt+i] + (TB[2][4] - 1)*sc[19*npt+i] + (TB[2][5] - 1)*sc[20*npt+i] + (TB[2][6] - 1)*sc[21*npt+i] + (TB[2][7] - 1)*sc[12*npt+i] + (TB[2][8] - 1)*sc[22*npt+i] + (TB[2][9] - 1)*sc[13*npt+i] + (TB[2][10] - 1)*sc[23*npt+i] + (TB[2][11] - 1)*sc[24*npt+i] + (TB[2][12] - 1)*sc[25*npt+i] + (TB[2][13] - 1)*sc[26*npt+i] + (TB[2][14] - 1)*sc[27*npt+i] + (TB[2][15] - 1)*sc[28*npt+i] + (TB[2][16] - 1)*sc[17*npt+i];
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
        phi_r = sc[10*npt+i];
        Kc = Kc_s[2*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 4: TXCH2 + H (+M) <=> CH3 (+M) */
        phi_f = sc[1*npt+i]*sc[21*npt+i];
        alpha = mixture[i] + (TB[3][0] - 1)*sc[2*npt+i] + (TB[3][1] - 1)*sc[7*npt+i] + (TB[3][2] - 1)*sc[18*npt+i] + (TB[3][3] - 1)*sc[9*npt+i] + (TB[3][4] - 1)*sc[19*npt+i] + (TB[3][5] - 1)*sc[20*npt+i] + (TB[3][6] - 1)*sc[21*npt+i] + (TB[3][7] - 1)*sc[12*npt+i] + (TB[3][8] - 1)*sc[22*npt+i] + (TB[3][9] - 1)*sc[13*npt+i] + (TB[3][10] - 1)*sc[23*npt+i] + (TB[3][11] - 1)*sc[24*npt+i] + (TB[3][12] - 1)*sc[25*npt+i] + (TB[3][13] - 1)*sc[26*npt+i] + (TB[3][14] - 1)*sc[27*npt+i] + (TB[3][15] - 1)*sc[28*npt+i] + (TB[3][16] - 1)*sc[17*npt+i];
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
        phi_r = sc[10*npt+i];
        Kc = Kc_s[3*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;
        wdot[21*npt+i] -= qdot;

        /*reaction 5: CH3 + H (+M) <=> CH4 (+M) */
        phi_f = sc[1*npt+i]*sc[10*npt+i];
        alpha = mixture[i] + (TB[4][0] - 1)*sc[2*npt+i] + (TB[4][1] - 1)*sc[7*npt+i] + (TB[4][2] - 1)*sc[18*npt+i] + (TB[4][3] - 1)*sc[9*npt+i] + (TB[4][4] - 1)*sc[19*npt+i] + (TB[4][5] - 1)*sc[20*npt+i] + (TB[4][6] - 1)*sc[21*npt+i] + (TB[4][7] - 1)*sc[12*npt+i] + (TB[4][8] - 1)*sc[22*npt+i] + (TB[4][9] - 1)*sc[13*npt+i] + (TB[4][10] - 1)*sc[23*npt+i] + (TB[4][11] - 1)*sc[24*npt+i] + (TB[4][12] - 1)*sc[25*npt+i] + (TB[4][13] - 1)*sc[26*npt+i] + (TB[4][14] - 1)*sc[27*npt+i] + (TB[4][15] - 1)*sc[28*npt+i] + (TB[4][16] - 1)*sc[17*npt+i];
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
        phi_r = sc[13*npt+i];
        Kc = Kc_s[4*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[10*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;

        /*reaction 6: TXCH2 + CO (+M) <=> CH2CO (+M) */
        phi_f = sc[9*npt+i]*sc[21*npt+i];
        alpha = mixture[i] + (TB[5][0] - 1)*sc[2*npt+i] + (TB[5][1] - 1)*sc[7*npt+i] + (TB[5][2] - 1)*sc[18*npt+i] + (TB[5][3] - 1)*sc[9*npt+i] + (TB[5][4] - 1)*sc[19*npt+i] + (TB[5][5] - 1)*sc[20*npt+i] + (TB[5][6] - 1)*sc[21*npt+i] + (TB[5][7] - 1)*sc[12*npt+i] + (TB[5][8] - 1)*sc[22*npt+i] + (TB[5][9] - 1)*sc[13*npt+i] + (TB[5][10] - 1)*sc[23*npt+i] + (TB[5][11] - 1)*sc[24*npt+i] + (TB[5][12] - 1)*sc[25*npt+i] + (TB[5][13] - 1)*sc[26*npt+i] + (TB[5][14] - 1)*sc[27*npt+i] + (TB[5][15] - 1)*sc[28*npt+i] + (TB[5][16] - 1)*sc[17*npt+i];
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
        phi_r = sc[16*npt+i];
        Kc = Kc_s[5*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[9*npt+i] -= qdot;
        wdot[16*npt+i] += qdot;
        wdot[21*npt+i] -= qdot;

        /*reaction 7: CO + H2 (+M) <=> CH2O (+M) */
        phi_f = sc[2*npt+i]*sc[9*npt+i];
        alpha = mixture[i] + (TB[6][0] - 1)*sc[2*npt+i] + (TB[6][1] - 1)*sc[7*npt+i] + (TB[6][2] - 1)*sc[18*npt+i] + (TB[6][3] - 1)*sc[9*npt+i] + (TB[6][4] - 1)*sc[19*npt+i] + (TB[6][5] - 1)*sc[20*npt+i] + (TB[6][6] - 1)*sc[21*npt+i] + (TB[6][7] - 1)*sc[12*npt+i] + (TB[6][8] - 1)*sc[22*npt+i] + (TB[6][9] - 1)*sc[13*npt+i] + (TB[6][10] - 1)*sc[23*npt+i] + (TB[6][11] - 1)*sc[24*npt+i] + (TB[6][12] - 1)*sc[25*npt+i] + (TB[6][13] - 1)*sc[26*npt+i] + (TB[6][14] - 1)*sc[27*npt+i] + (TB[6][15] - 1)*sc[28*npt+i] + (TB[6][16] - 1)*sc[17*npt+i];
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
        phi_r = sc[11*npt+i];
        Kc = Kc_s[6*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[9*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;

        /*reaction 8: CH + CO (+M) <=> HCCO (+M) */
        phi_f = sc[9*npt+i]*sc[19*npt+i];
        alpha = mixture[i] + (TB[7][0] - 1)*sc[2*npt+i] + (TB[7][1] - 1)*sc[7*npt+i] + (TB[7][2] - 1)*sc[18*npt+i] + (TB[7][3] - 1)*sc[9*npt+i] + (TB[7][4] - 1)*sc[19*npt+i] + (TB[7][5] - 1)*sc[20*npt+i] + (TB[7][6] - 1)*sc[21*npt+i] + (TB[7][7] - 1)*sc[12*npt+i] + (TB[7][8] - 1)*sc[22*npt+i] + (TB[7][9] - 1)*sc[13*npt+i] + (TB[7][10] - 1)*sc[23*npt+i] + (TB[7][11] - 1)*sc[24*npt+i] + (TB[7][12] - 1)*sc[25*npt+i] + (TB[7][13] - 1)*sc[26*npt+i] + (TB[7][14] - 1)*sc[27*npt+i] + (TB[7][15] - 1)*sc[28*npt+i] + (TB[7][16] - 1)*sc[17*npt+i];
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
        phi_r = sc[25*npt+i];
        Kc = Kc_s[7*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[9*npt+i] -= qdot;
        wdot[19*npt+i] -= qdot;
        wdot[25*npt+i] += qdot;

        /*reaction 9: CO + O (+M) <=> CO2 (+M) */
        phi_f = sc[3*npt+i]*sc[9*npt+i];
        alpha = mixture[i] + (TB[8][0] - 1)*sc[2*npt+i] + (TB[8][1] - 1)*sc[7*npt+i] + (TB[8][2] - 1)*sc[18*npt+i] + (TB[8][3] - 1)*sc[9*npt+i] + (TB[8][4] - 1)*sc[19*npt+i] + (TB[8][5] - 1)*sc[20*npt+i] + (TB[8][6] - 1)*sc[21*npt+i] + (TB[8][7] - 1)*sc[12*npt+i] + (TB[8][8] - 1)*sc[22*npt+i] + (TB[8][9] - 1)*sc[13*npt+i] + (TB[8][10] - 1)*sc[23*npt+i] + (TB[8][11] - 1)*sc[24*npt+i] + (TB[8][12] - 1)*sc[25*npt+i] + (TB[8][13] - 1)*sc[26*npt+i] + (TB[8][14] - 1)*sc[27*npt+i] + (TB[8][15] - 1)*sc[28*npt+i] + (TB[8][16] - 1)*sc[17*npt+i];
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
        phi_r = sc[12*npt+i];
        Kc = Kc_s[8*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[9*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;

        /*reaction 10: HCO + H (+M) <=> CH2O (+M) */
        phi_f = sc[1*npt+i]*sc[20*npt+i];
        alpha = mixture[i] + (TB[9][0] - 1)*sc[2*npt+i] + (TB[9][1] - 1)*sc[7*npt+i] + (TB[9][2] - 1)*sc[18*npt+i] + (TB[9][3] - 1)*sc[9*npt+i] + (TB[9][4] - 1)*sc[19*npt+i] + (TB[9][5] - 1)*sc[20*npt+i] + (TB[9][6] - 1)*sc[21*npt+i] + (TB[9][7] - 1)*sc[12*npt+i] + (TB[9][8] - 1)*sc[22*npt+i] + (TB[9][9] - 1)*sc[13*npt+i] + (TB[9][10] - 1)*sc[23*npt+i] + (TB[9][11] - 1)*sc[24*npt+i] + (TB[9][12] - 1)*sc[25*npt+i] + (TB[9][13] - 1)*sc[26*npt+i] + (TB[9][14] - 1)*sc[27*npt+i] + (TB[9][15] - 1)*sc[28*npt+i] + (TB[9][16] - 1)*sc[17*npt+i];
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
        phi_r = sc[11*npt+i];
        Kc = Kc_s[9*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;
        wdot[20*npt+i] -= qdot;

        /*reaction 11: C2H2 + H (+M) <=> C2H3 (+M) */
        phi_f = sc[1*npt+i]*sc[14*npt+i];
        alpha = mixture[i] + (TB[10][0] - 1)*sc[2*npt+i] + (TB[10][1] - 1)*sc[7*npt+i] + (TB[10][2] - 1)*sc[18*npt+i] + (TB[10][3] - 1)*sc[9*npt+i] + (TB[10][4] - 1)*sc[19*npt+i] + (TB[10][5] - 1)*sc[20*npt+i] + (TB[10][6] - 1)*sc[21*npt+i] + (TB[10][7] - 1)*sc[12*npt+i] + (TB[10][8] - 1)*sc[22*npt+i] + (TB[10][9] - 1)*sc[13*npt+i] + (TB[10][10] - 1)*sc[23*npt+i] + (TB[10][11] - 1)*sc[24*npt+i] + (TB[10][12] - 1)*sc[25*npt+i] + (TB[10][13] - 1)*sc[26*npt+i] + (TB[10][14] - 1)*sc[27*npt+i] + (TB[10][15] - 1)*sc[28*npt+i] + (TB[10][16] - 1)*sc[17*npt+i];
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
        phi_r = sc[23*npt+i];
        Kc = Kc_s[10*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[14*npt+i] -= qdot;
        wdot[23*npt+i] += qdot;

        /*reaction 12: C2H3 + H (+M) <=> C2H4 (+M) */
        phi_f = sc[1*npt+i]*sc[23*npt+i];
        alpha = mixture[i] + (TB[11][0] - 1)*sc[2*npt+i] + (TB[11][1] - 1)*sc[7*npt+i] + (TB[11][2] - 1)*sc[18*npt+i] + (TB[11][3] - 1)*sc[9*npt+i] + (TB[11][4] - 1)*sc[19*npt+i] + (TB[11][5] - 1)*sc[20*npt+i] + (TB[11][6] - 1)*sc[21*npt+i] + (TB[11][7] - 1)*sc[12*npt+i] + (TB[11][8] - 1)*sc[22*npt+i] + (TB[11][9] - 1)*sc[13*npt+i] + (TB[11][10] - 1)*sc[23*npt+i] + (TB[11][11] - 1)*sc[24*npt+i] + (TB[11][12] - 1)*sc[25*npt+i] + (TB[11][13] - 1)*sc[26*npt+i] + (TB[11][14] - 1)*sc[27*npt+i] + (TB[11][15] - 1)*sc[28*npt+i] + (TB[11][16] - 1)*sc[17*npt+i];
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
        phi_r = sc[15*npt+i];
        Kc = Kc_s[11*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;
        wdot[23*npt+i] -= qdot;

        /*reaction 13: C2H4 + H (+M) <=> C2H5 (+M) */
        phi_f = sc[1*npt+i]*sc[15*npt+i];
        alpha = mixture[i] + (TB[12][0] - 1)*sc[2*npt+i] + (TB[12][1] - 1)*sc[7*npt+i] + (TB[12][2] - 1)*sc[18*npt+i] + (TB[12][3] - 1)*sc[9*npt+i] + (TB[12][4] - 1)*sc[19*npt+i] + (TB[12][5] - 1)*sc[20*npt+i] + (TB[12][6] - 1)*sc[21*npt+i] + (TB[12][7] - 1)*sc[12*npt+i] + (TB[12][8] - 1)*sc[22*npt+i] + (TB[12][9] - 1)*sc[13*npt+i] + (TB[12][10] - 1)*sc[23*npt+i] + (TB[12][11] - 1)*sc[24*npt+i] + (TB[12][12] - 1)*sc[25*npt+i] + (TB[12][13] - 1)*sc[26*npt+i] + (TB[12][14] - 1)*sc[27*npt+i] + (TB[12][15] - 1)*sc[28*npt+i] + (TB[12][16] - 1)*sc[17*npt+i];
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
        phi_r = sc[24*npt+i];
        Kc = Kc_s[12*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[15*npt+i] -= qdot;
        wdot[24*npt+i] += qdot;

        /*reaction 14: C2H5 + H (+M) <=> C2H6 (+M) */
        phi_f = sc[1*npt+i]*sc[24*npt+i];
        alpha = mixture[i] + (TB[13][0] - 1)*sc[2*npt+i] + (TB[13][1] - 1)*sc[7*npt+i] + (TB[13][2] - 1)*sc[18*npt+i] + (TB[13][3] - 1)*sc[9*npt+i] + (TB[13][4] - 1)*sc[19*npt+i] + (TB[13][5] - 1)*sc[20*npt+i] + (TB[13][6] - 1)*sc[21*npt+i] + (TB[13][7] - 1)*sc[12*npt+i] + (TB[13][8] - 1)*sc[22*npt+i] + (TB[13][9] - 1)*sc[13*npt+i] + (TB[13][10] - 1)*sc[23*npt+i] + (TB[13][11] - 1)*sc[24*npt+i] + (TB[13][12] - 1)*sc[25*npt+i] + (TB[13][13] - 1)*sc[26*npt+i] + (TB[13][14] - 1)*sc[27*npt+i] + (TB[13][15] - 1)*sc[28*npt+i] + (TB[13][16] - 1)*sc[17*npt+i];
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
        phi_r = sc[17*npt+i];
        Kc = Kc_s[13*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[17*npt+i] += qdot;
        wdot[24*npt+i] -= qdot;

        /*reaction 15: C2H6 (+M) <=> 2.000000 CH3 (+M) */
        phi_f = sc[17*npt+i];
        alpha = mixture[i] + (TB[14][0] - 1)*sc[2*npt+i] + (TB[14][1] - 1)*sc[7*npt+i] + (TB[14][2] - 1)*sc[18*npt+i] + (TB[14][3] - 1)*sc[9*npt+i] + (TB[14][4] - 1)*sc[19*npt+i] + (TB[14][5] - 1)*sc[20*npt+i] + (TB[14][6] - 1)*sc[21*npt+i] + (TB[14][7] - 1)*sc[12*npt+i] + (TB[14][8] - 1)*sc[22*npt+i] + (TB[14][9] - 1)*sc[13*npt+i] + (TB[14][10] - 1)*sc[23*npt+i] + (TB[14][11] - 1)*sc[24*npt+i] + (TB[14][12] - 1)*sc[25*npt+i] + (TB[14][13] - 1)*sc[26*npt+i] + (TB[14][14] - 1)*sc[27*npt+i] + (TB[14][15] - 1)*sc[28*npt+i] + (TB[14][16] - 1)*sc[17*npt+i];
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
        phi_r = pow(sc[10*npt+i], 2.000000);
        Kc = Kc_s[14*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[10*npt+i] += 2.000000 * qdot;
        wdot[17*npt+i] -= qdot;

        /*reaction 16: 2.000000 H + M <=> H2 + M */
        phi_f = pow(sc[1*npt+i], 2.000000);
        alpha = mixture[i] + (TB[15][0] - 1)*sc[2*npt+i] + (TB[15][1] - 1)*sc[7*npt+i] + (TB[15][2] - 1)*sc[18*npt+i] + (TB[15][3] - 1)*sc[19*npt+i] + (TB[15][4] - 1)*sc[20*npt+i] + (TB[15][5] - 1)*sc[21*npt+i] + (TB[15][6] - 1)*sc[12*npt+i] + (TB[15][7] - 1)*sc[22*npt+i] + (TB[15][8] - 1)*sc[13*npt+i] + (TB[15][9] - 1)*sc[23*npt+i] + (TB[15][10] - 1)*sc[24*npt+i] + (TB[15][11] - 1)*sc[25*npt+i] + (TB[15][12] - 1)*sc[26*npt+i] + (TB[15][13] - 1)*sc[27*npt+i] + (TB[15][14] - 1)*sc[28*npt+i] + (TB[15][15] - 1)*sc[17*npt+i];
        k_f = alpha * k_f_s[15*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i];
        Kc = Kc_s[15*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= 2.000000 * qdot;
        wdot[2*npt+i] += qdot;

        /*reaction 17: H + O + M <=> OH + M */
        phi_f = sc[1*npt+i]*sc[3*npt+i];
        alpha = mixture[i] + (TB[16][0] - 1)*sc[2*npt+i] + (TB[16][1] - 1)*sc[7*npt+i] + (TB[16][2] - 1)*sc[18*npt+i] + (TB[16][3] - 1)*sc[9*npt+i] + (TB[16][4] - 1)*sc[19*npt+i] + (TB[16][5] - 1)*sc[20*npt+i] + (TB[16][6] - 1)*sc[21*npt+i] + (TB[16][7] - 1)*sc[12*npt+i] + (TB[16][8] - 1)*sc[22*npt+i] + (TB[16][9] - 1)*sc[13*npt+i] + (TB[16][10] - 1)*sc[23*npt+i] + (TB[16][11] - 1)*sc[24*npt+i] + (TB[16][12] - 1)*sc[25*npt+i] + (TB[16][13] - 1)*sc[26*npt+i] + (TB[16][14] - 1)*sc[27*npt+i] + (TB[16][15] - 1)*sc[28*npt+i] + (TB[16][16] - 1)*sc[17*npt+i];
        k_f = alpha * k_f_s[16*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i];
        Kc = Kc_s[16*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;

        /*reaction 18: 2.000000 O + M <=> O2 + M */
        phi_f = pow(sc[3*npt+i], 2.000000);
        alpha = mixture[i] + (TB[17][0] - 1)*sc[2*npt+i] + (TB[17][1] - 1)*sc[7*npt+i] + (TB[17][2] - 1)*sc[18*npt+i] + (TB[17][3] - 1)*sc[9*npt+i] + (TB[17][4] - 1)*sc[19*npt+i] + (TB[17][5] - 1)*sc[20*npt+i] + (TB[17][6] - 1)*sc[21*npt+i] + (TB[17][7] - 1)*sc[12*npt+i] + (TB[17][8] - 1)*sc[22*npt+i] + (TB[17][9] - 1)*sc[13*npt+i] + (TB[17][10] - 1)*sc[23*npt+i] + (TB[17][11] - 1)*sc[24*npt+i] + (TB[17][12] - 1)*sc[25*npt+i] + (TB[17][13] - 1)*sc[26*npt+i] + (TB[17][14] - 1)*sc[27*npt+i] + (TB[17][15] - 1)*sc[28*npt+i] + (TB[17][16] - 1)*sc[17*npt+i];
        k_f = alpha * k_f_s[17*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i];
        Kc = Kc_s[17*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= 2.000000 * qdot;
        wdot[5*npt+i] += qdot;

        /*reaction 19: H + OH + M <=> H2O + M */
        phi_f = sc[1*npt+i]*sc[4*npt+i];
        alpha = mixture[i] + (TB[18][0] - 1)*sc[2*npt+i] + (TB[18][1] - 1)*sc[7*npt+i] + (TB[18][2] - 1)*sc[18*npt+i] + (TB[18][3] - 1)*sc[9*npt+i] + (TB[18][4] - 1)*sc[19*npt+i] + (TB[18][5] - 1)*sc[20*npt+i] + (TB[18][6] - 1)*sc[21*npt+i] + (TB[18][7] - 1)*sc[12*npt+i] + (TB[18][8] - 1)*sc[22*npt+i] + (TB[18][9] - 1)*sc[13*npt+i] + (TB[18][10] - 1)*sc[23*npt+i] + (TB[18][11] - 1)*sc[24*npt+i] + (TB[18][12] - 1)*sc[25*npt+i] + (TB[18][13] - 1)*sc[26*npt+i] + (TB[18][14] - 1)*sc[27*npt+i] + (TB[18][15] - 1)*sc[28*npt+i] + (TB[18][16] - 1)*sc[17*npt+i];
        k_f = alpha * k_f_s[18*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i];
        Kc = Kc_s[18*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[4*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;

        /*reaction 20: HCO + M <=> CO + H + M */
        phi_f = sc[20*npt+i];
        alpha = mixture[i] + (TB[19][0] - 1)*sc[2*npt+i] + (TB[19][1] - 1)*sc[7*npt+i] + (TB[19][2] - 1)*sc[18*npt+i] + (TB[19][3] - 1)*sc[9*npt+i] + (TB[19][4] - 1)*sc[19*npt+i] + (TB[19][5] - 1)*sc[20*npt+i] + (TB[19][6] - 1)*sc[21*npt+i] + (TB[19][7] - 1)*sc[12*npt+i] + (TB[19][8] - 1)*sc[22*npt+i] + (TB[19][9] - 1)*sc[13*npt+i] + (TB[19][10] - 1)*sc[23*npt+i] + (TB[19][11] - 1)*sc[24*npt+i] + (TB[19][12] - 1)*sc[25*npt+i] + (TB[19][13] - 1)*sc[26*npt+i] + (TB[19][14] - 1)*sc[27*npt+i] + (TB[19][15] - 1)*sc[28*npt+i] + (TB[19][16] - 1)*sc[17*npt+i];
        k_f = alpha * k_f_s[19*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[9*npt+i];
        Kc = Kc_s[19*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[9*npt+i] += qdot;
        wdot[20*npt+i] -= qdot;

        /*reaction 21: 2.000000 H + H2 <=> 2.000000 H2 */
        phi_f = pow(sc[1*npt+i], 2.000000)*sc[2*npt+i];
        k_f = k_f_s[20*npt+i];
        q_f = phi_f * k_f;
        phi_r = pow(sc[2*npt+i], 2.000000);
        Kc = Kc_s[20*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= 2.000000 * qdot;
        wdot[2*npt+i] -= qdot;
        wdot[2*npt+i] += 2.000000 * qdot;

        /*reaction 22: O + H2 <=> H + OH */
        phi_f = sc[2*npt+i]*sc[3*npt+i];
        k_f = k_f_s[21*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[4*npt+i];
        Kc = Kc_s[21*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;

        /*reaction 23: 2.000000 OH <=> O + H2O */
        phi_f = pow(sc[4*npt+i], 2.000000);
        k_f = k_f_s[22*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[7*npt+i];
        Kc = Kc_s[22*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[4*npt+i] -= 2.000000 * qdot;
        wdot[7*npt+i] += qdot;

        /*reaction 24: OH + H2 <=> H + H2O */
        phi_f = sc[2*npt+i]*sc[4*npt+i];
        k_f = k_f_s[23*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[7*npt+i];
        Kc = Kc_s[23*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;

        /*reaction 25: 2.000000 H + H2O <=> H2 + H2O */
        phi_f = pow(sc[1*npt+i], 2.000000)*sc[7*npt+i];
        k_f = k_f_s[24*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[7*npt+i];
        Kc = Kc_s[24*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= 2.000000 * qdot;
        wdot[2*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;

        /*reaction 26: H + O2 <=> O + OH */
        phi_f = sc[1*npt+i]*sc[5*npt+i];
        k_f = k_f_s[25*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[4*npt+i];
        Kc = Kc_s[25*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[4*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;

        /*reaction 27: H2 + O2 <=> HO2 + H */
        phi_f = sc[2*npt+i]*sc[5*npt+i];
        k_f = k_f_s[26*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[8*npt+i];
        Kc = Kc_s[26*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[5*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;

        /*reaction 28: HO2 + OH <=> H2O + O2 */
        phi_f = sc[4*npt+i]*sc[8*npt+i];
        k_f = k_f_s[27*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[7*npt+i];
        Kc = Kc_s[27*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[7*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;

        /*reaction 29: HO2 + H <=> 2.000000 OH */
        phi_f = sc[1*npt+i]*sc[8*npt+i];
        k_f = k_f_s[28*npt+i];
        q_f = phi_f * k_f;
        phi_r = pow(sc[4*npt+i], 2.000000);
        Kc = Kc_s[28*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[4*npt+i] += 2.000000 * qdot;
        wdot[8*npt+i] -= qdot;

        /*reaction 30: HO2 + O <=> OH + O2 */
        phi_f = sc[3*npt+i]*sc[8*npt+i];
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
        wdot[8*npt+i] -= qdot;

        /*reaction 31: HO2 + H <=> O + H2O */
        phi_f = sc[1*npt+i]*sc[8*npt+i];
        k_f = k_f_s[30*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[7*npt+i];
        Kc = Kc_s[30*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[7*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;

        /*reaction 32: HO2 + OH <=> H2O + O2 */
        phi_f = sc[4*npt+i]*sc[8*npt+i];
        k_f = k_f_s[31*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[7*npt+i];
        Kc = Kc_s[31*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[7*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;

        /*reaction 33: H2O2 + O <=> HO2 + OH */
        phi_f = sc[3*npt+i]*sc[6*npt+i];
        k_f = k_f_s[32*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[8*npt+i];
        Kc = Kc_s[32*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;

        /*reaction 34: H2O2 + H <=> H2O + OH */
        phi_f = sc[1*npt+i]*sc[6*npt+i];
        k_f = k_f_s[33*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[7*npt+i];
        Kc = Kc_s[33*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;

        /*reaction 35: H2O2 + H <=> HO2 + H2 */
        phi_f = sc[1*npt+i]*sc[6*npt+i];
        k_f = k_f_s[34*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[8*npt+i];
        Kc = Kc_s[34*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;

        /*reaction 36: H2O2 + OH <=> HO2 + H2O */
        phi_f = sc[4*npt+i]*sc[6*npt+i];
        k_f = k_f_s[35*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[8*npt+i];
        Kc = Kc_s[35*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[6*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[8*npt+i] += qdot;

        /*reaction 37: H2O2 + OH <=> HO2 + H2O */
        phi_f = sc[4*npt+i]*sc[6*npt+i];
        k_f = k_f_s[36*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[8*npt+i];
        Kc = Kc_s[36*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[6*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[8*npt+i] += qdot;

        /*reaction 38: C + O2 <=> CO + O */
        phi_f = sc[5*npt+i]*sc[18*npt+i];
        k_f = k_f_s[37*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[9*npt+i];
        Kc = Kc_s[37*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;
        wdot[18*npt+i] -= qdot;

        /*reaction 39: C + OH <=> CO + H */
        phi_f = sc[4*npt+i]*sc[18*npt+i];
        k_f = k_f_s[38*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[9*npt+i];
        Kc = Kc_s[38*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;
        wdot[18*npt+i] -= qdot;

        /*reaction 40: CH + OH <=> HCO + H */
        phi_f = sc[4*npt+i]*sc[19*npt+i];
        k_f = k_f_s[39*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[20*npt+i];
        Kc = Kc_s[39*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[19*npt+i] -= qdot;
        wdot[20*npt+i] += qdot;

        /*reaction 41: CH + H2 <=> TXCH2 + H */
        phi_f = sc[2*npt+i]*sc[19*npt+i];
        k_f = k_f_s[40*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[21*npt+i];
        Kc = Kc_s[40*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[19*npt+i] -= qdot;
        wdot[21*npt+i] += qdot;

        /*reaction 42: CH + O <=> CO + H */
        phi_f = sc[3*npt+i]*sc[19*npt+i];
        k_f = k_f_s[41*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[9*npt+i];
        Kc = Kc_s[41*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 43: CH + O2 <=> HCO + O */
        phi_f = sc[5*npt+i]*sc[19*npt+i];
        k_f = k_f_s[42*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[20*npt+i];
        Kc = Kc_s[42*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[19*npt+i] -= qdot;
        wdot[20*npt+i] += qdot;

        /*reaction 44: CH + H <=> C + H2 */
        phi_f = sc[1*npt+i]*sc[19*npt+i];
        k_f = k_f_s[43*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[18*npt+i];
        Kc = Kc_s[43*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[18*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 45: CH + H2O <=> CH2O + H */
        phi_f = sc[7*npt+i]*sc[19*npt+i];
        k_f = k_f_s[44*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[11*npt+i];
        Kc = Kc_s[44*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 46: TXCH2 + O2 => OH + H + CO */
        phi_f = sc[5*npt+i]*sc[21*npt+i];
        k_f = k_f_s[45*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[4*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;
        wdot[21*npt+i] -= qdot;

        /*reaction 47: TXCH2 + O2 <=> CH2O + O */
        phi_f = sc[5*npt+i]*sc[21*npt+i];
        k_f = k_f_s[46*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[11*npt+i];
        Kc = Kc_s[46*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;
        wdot[21*npt+i] -= qdot;

        /*reaction 48: TXCH2 + OH <=> CH2O + H */
        phi_f = sc[4*npt+i]*sc[21*npt+i];
        k_f = k_f_s[47*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[11*npt+i];
        Kc = Kc_s[47*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;
        wdot[21*npt+i] -= qdot;

        /*reaction 49: TXCH2 + HO2 <=> CH2O + OH */
        phi_f = sc[8*npt+i]*sc[21*npt+i];
        k_f = k_f_s[48*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[11*npt+i];
        Kc = Kc_s[48*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;
        wdot[21*npt+i] -= qdot;

        /*reaction 50: TXCH2 + O2 => CO2 + 2.000000 H */
        phi_f = sc[5*npt+i]*sc[21*npt+i];
        k_f = k_f_s[49*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[1*npt+i] += 2.000000 * qdot;
        wdot[5*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;
        wdot[21*npt+i] -= qdot;
    }
}

void vcomp_wdot_51_100(int npt, amrex::Real *  wdot, amrex::Real *  mixture, amrex::Real *  sc,
		amrex::Real *  k_f_s, amrex::Real *  Kc_s,
		amrex::Real *  tc, amrex::Real *  invT, amrex::Real *  T)
{
    for (int i=0; i<npt; i++) {
        amrex::Real qdot, q_f, q_r, phi_f, phi_r, k_f, k_r, Kc;
        amrex::Real alpha;

        /*reaction 51: TXCH2 + OH <=> CH + H2O */
        phi_f = sc[4*npt+i]*sc[21*npt+i];
        k_f = k_f_s[50*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[19*npt+i];
        Kc = Kc_s[50*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[19*npt+i] += qdot;
        wdot[21*npt+i] -= qdot;

        /*reaction 52: TXCH2 + O <=> HCO + H */
        phi_f = sc[3*npt+i]*sc[21*npt+i];
        k_f = k_f_s[51*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[20*npt+i];
        Kc = Kc_s[51*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[20*npt+i] += qdot;
        wdot[21*npt+i] -= qdot;

        /*reaction 53: TXCH2 + H2 <=> H + CH3 */
        phi_f = sc[2*npt+i]*sc[21*npt+i];
        k_f = k_f_s[52*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[10*npt+i];
        Kc = Kc_s[52*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;
        wdot[21*npt+i] -= qdot;

        /*reaction 54: SXCH2 + H2O <=> TXCH2 + H2O */
        phi_f = sc[7*npt+i]*sc[22*npt+i];
        k_f = k_f_s[53*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[21*npt+i];
        Kc = Kc_s[53*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[7*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[21*npt+i] += qdot;
        wdot[22*npt+i] -= qdot;

        /*reaction 55: SXCH2 + H <=> CH + H2 */
        phi_f = sc[1*npt+i]*sc[22*npt+i];
        k_f = k_f_s[54*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[19*npt+i];
        Kc = Kc_s[54*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[19*npt+i] += qdot;
        wdot[22*npt+i] -= qdot;

        /*reaction 56: SXCH2 + O2 <=> H + OH + CO */
        phi_f = sc[5*npt+i]*sc[22*npt+i];
        k_f = k_f_s[55*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[4*npt+i]*sc[9*npt+i];
        Kc = Kc_s[55*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[4*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;
        wdot[22*npt+i] -= qdot;

        /*reaction 57: SXCH2 + O <=> CO + H2 */
        phi_f = sc[3*npt+i]*sc[22*npt+i];
        k_f = k_f_s[56*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[9*npt+i];
        Kc = Kc_s[56*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;
        wdot[22*npt+i] -= qdot;

        /*reaction 58: SXCH2 + O2 <=> CO + H2O */
        phi_f = sc[5*npt+i]*sc[22*npt+i];
        k_f = k_f_s[57*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[9*npt+i];
        Kc = Kc_s[57*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[9*npt+i] += qdot;
        wdot[22*npt+i] -= qdot;

        /*reaction 59: SXCH2 + H2 <=> CH3 + H */
        phi_f = sc[2*npt+i]*sc[22*npt+i];
        k_f = k_f_s[58*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[10*npt+i];
        Kc = Kc_s[58*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;
        wdot[22*npt+i] -= qdot;

        /*reaction 60: SXCH2 + O <=> HCO + H */
        phi_f = sc[3*npt+i]*sc[22*npt+i];
        k_f = k_f_s[59*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[20*npt+i];
        Kc = Kc_s[59*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[20*npt+i] += qdot;
        wdot[22*npt+i] -= qdot;

        /*reaction 61: SXCH2 + H2O => H2 + CH2O */
        phi_f = sc[7*npt+i]*sc[22*npt+i];
        k_f = k_f_s[60*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;
        wdot[22*npt+i] -= qdot;

        /*reaction 62: SXCH2 + OH <=> CH2O + H */
        phi_f = sc[4*npt+i]*sc[22*npt+i];
        k_f = k_f_s[61*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[11*npt+i];
        Kc = Kc_s[61*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;
        wdot[22*npt+i] -= qdot;

        /*reaction 63: CH3 + OH => H2 + CH2O */
        phi_f = sc[4*npt+i]*sc[10*npt+i];
        k_f = k_f_s[62*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[10*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;

        /*reaction 64: CH3 + H2O2 <=> CH4 + HO2 */
        phi_f = sc[6*npt+i]*sc[10*npt+i];
        k_f = k_f_s[63*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i]*sc[13*npt+i];
        Kc = Kc_s[63*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;

        /*reaction 65: CH3 + O2 <=> CH2O + OH */
        phi_f = sc[5*npt+i]*sc[10*npt+i];
        k_f = k_f_s[64*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[11*npt+i];
        Kc = Kc_s[64*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[10*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;

        /*reaction 66: CH3 + CH <=> C2H3 + H */
        phi_f = sc[10*npt+i]*sc[19*npt+i];
        k_f = k_f_s[65*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[23*npt+i];
        Kc = Kc_s[65*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;
        wdot[19*npt+i] -= qdot;
        wdot[23*npt+i] += qdot;

        /*reaction 67: CH3 + O <=> CH2O + H */
        phi_f = sc[3*npt+i]*sc[10*npt+i];
        k_f = k_f_s[66*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[11*npt+i];
        Kc = Kc_s[66*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[10*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;

        /*reaction 68: CH3 + C <=> C2H2 + H */
        phi_f = sc[10*npt+i]*sc[18*npt+i];
        k_f = k_f_s[67*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[14*npt+i];
        Kc = Kc_s[67*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;
        wdot[18*npt+i] -= qdot;

        /*reaction 69: CH3 + OH <=> TXCH2 + H2O */
        phi_f = sc[4*npt+i]*sc[10*npt+i];
        k_f = k_f_s[68*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[21*npt+i];
        Kc = Kc_s[68*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;
        wdot[21*npt+i] += qdot;

        /*reaction 70: CH3 + SXCH2 <=> C2H4 + H */
        phi_f = sc[10*npt+i]*sc[22*npt+i];
        k_f = k_f_s[69*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[15*npt+i];
        Kc = Kc_s[69*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;
        wdot[22*npt+i] -= qdot;

        /*reaction 71: CH3 + OH <=> SXCH2 + H2O */
        phi_f = sc[4*npt+i]*sc[10*npt+i];
        k_f = k_f_s[70*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[22*npt+i];
        Kc = Kc_s[70*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;
        wdot[22*npt+i] += qdot;

        /*reaction 72: 2.000000 CH3 <=> C2H5 + H */
        phi_f = pow(sc[10*npt+i], 2.000000);
        k_f = k_f_s[71*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[24*npt+i];
        Kc = Kc_s[71*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[10*npt+i] -= 2.000000 * qdot;
        wdot[24*npt+i] += qdot;

        /*reaction 73: CH3 + HO2 <=> CH4 + O2 */
        phi_f = sc[8*npt+i]*sc[10*npt+i];
        k_f = k_f_s[72*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[13*npt+i];
        Kc = Kc_s[72*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[10*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;

        /*reaction 74: CH3 + TXCH2 <=> C2H4 + H */
        phi_f = sc[10*npt+i]*sc[21*npt+i];
        k_f = k_f_s[73*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[15*npt+i];
        Kc = Kc_s[73*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;
        wdot[21*npt+i] -= qdot;

        /*reaction 75: CH3 + O => H + H2 + CO */
        phi_f = sc[3*npt+i]*sc[10*npt+i];
        k_f = k_f_s[74*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;

        /*reaction 76: CH4 + CH <=> C2H4 + H */
        phi_f = sc[13*npt+i]*sc[19*npt+i];
        k_f = k_f_s[75*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[15*npt+i];
        Kc = Kc_s[75*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 77: CH4 + SXCH2 <=> 2.000000 CH3 */
        phi_f = sc[13*npt+i]*sc[22*npt+i];
        k_f = k_f_s[76*npt+i];
        q_f = phi_f * k_f;
        phi_r = pow(sc[10*npt+i], 2.000000);
        Kc = Kc_s[76*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[10*npt+i] += 2.000000 * qdot;
        wdot[13*npt+i] -= qdot;
        wdot[22*npt+i] -= qdot;

        /*reaction 78: CH4 + O <=> CH3 + OH */
        phi_f = sc[3*npt+i]*sc[13*npt+i];
        k_f = k_f_s[77*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[10*npt+i];
        Kc = Kc_s[77*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[10*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;

        /*reaction 79: CH4 + OH <=> CH3 + H2O */
        phi_f = sc[4*npt+i]*sc[13*npt+i];
        k_f = k_f_s[78*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[10*npt+i];
        Kc = Kc_s[78*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[10*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;

        /*reaction 80: CH4 + TXCH2 <=> 2.000000 CH3 */
        phi_f = sc[13*npt+i]*sc[21*npt+i];
        k_f = k_f_s[79*npt+i];
        q_f = phi_f * k_f;
        phi_r = pow(sc[10*npt+i], 2.000000);
        Kc = Kc_s[79*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[10*npt+i] += 2.000000 * qdot;
        wdot[13*npt+i] -= qdot;
        wdot[21*npt+i] -= qdot;

        /*reaction 81: CH4 + H <=> CH3 + H2 */
        phi_f = sc[1*npt+i]*sc[13*npt+i];
        k_f = k_f_s[80*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[10*npt+i];
        Kc = Kc_s[80*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[10*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;

        /*reaction 82: SXCH2 + CO <=> TXCH2 + CO */
        phi_f = sc[9*npt+i]*sc[22*npt+i];
        k_f = k_f_s[81*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[9*npt+i]*sc[21*npt+i];
        Kc = Kc_s[81*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[9*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;
        wdot[21*npt+i] += qdot;
        wdot[22*npt+i] -= qdot;

        /*reaction 83: CO + O2 <=> CO2 + O */
        phi_f = sc[5*npt+i]*sc[9*npt+i];
        k_f = k_f_s[82*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[12*npt+i];
        Kc = Kc_s[82*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[9*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;

        /*reaction 84: CO + OH <=> CO2 + H */
        phi_f = sc[4*npt+i]*sc[9*npt+i];
        k_f = k_f_s[83*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[12*npt+i];
        Kc = Kc_s[83*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[9*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;

        /*reaction 85: CO + OH <=> CO2 + H */
        phi_f = sc[4*npt+i]*sc[9*npt+i];
        k_f = k_f_s[84*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[12*npt+i];
        Kc = Kc_s[84*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[9*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;

        /*reaction 86: CO + HO2 <=> CO2 + OH */
        phi_f = sc[8*npt+i]*sc[9*npt+i];
        k_f = k_f_s[85*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[12*npt+i];
        Kc = Kc_s[85*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[9*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;

        /*reaction 87: HCO + H <=> CO + H2 */
        phi_f = sc[1*npt+i]*sc[20*npt+i];
        k_f = k_f_s[86*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[9*npt+i];
        Kc = Kc_s[86*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[9*npt+i] += qdot;
        wdot[20*npt+i] -= qdot;

        /*reaction 88: CH3 + HCO <=> CH3CHO */
        phi_f = sc[10*npt+i]*sc[20*npt+i];
        k_f = k_f_s[87*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[26*npt+i];
        Kc = Kc_s[87*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[10*npt+i] -= qdot;
        wdot[20*npt+i] -= qdot;
        wdot[26*npt+i] += qdot;

        /*reaction 89: HCO + H2O <=> CO + H + H2O */
        phi_f = sc[7*npt+i]*sc[20*npt+i];
        k_f = k_f_s[88*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[7*npt+i]*sc[9*npt+i];
        Kc = Kc_s[88*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[9*npt+i] += qdot;
        wdot[20*npt+i] -= qdot;

        /*reaction 90: HCO + O <=> CO + OH */
        phi_f = sc[3*npt+i]*sc[20*npt+i];
        k_f = k_f_s[89*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[9*npt+i];
        Kc = Kc_s[89*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[9*npt+i] += qdot;
        wdot[20*npt+i] -= qdot;

        /*reaction 91: HCO + OH <=> CO + H2O */
        phi_f = sc[4*npt+i]*sc[20*npt+i];
        k_f = k_f_s[90*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[9*npt+i];
        Kc = Kc_s[90*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[9*npt+i] += qdot;
        wdot[20*npt+i] -= qdot;

        /*reaction 92: CH3 + HCO <=> CH4 + CO */
        phi_f = sc[10*npt+i]*sc[20*npt+i];
        k_f = k_f_s[91*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[9*npt+i]*sc[13*npt+i];
        Kc = Kc_s[91*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[9*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;
        wdot[20*npt+i] -= qdot;

        /*reaction 93: HCO + O <=> CO2 + H */
        phi_f = sc[3*npt+i]*sc[20*npt+i];
        k_f = k_f_s[92*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[12*npt+i];
        Kc = Kc_s[92*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;
        wdot[20*npt+i] -= qdot;

        /*reaction 94: HCO + O2 <=> CO + HO2 */
        phi_f = sc[5*npt+i]*sc[20*npt+i];
        k_f = k_f_s[93*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i]*sc[9*npt+i];
        Kc = Kc_s[93*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[9*npt+i] += qdot;
        wdot[20*npt+i] -= qdot;

        /*reaction 95: CH2O + H <=> HCO + H2 */
        phi_f = sc[1*npt+i]*sc[11*npt+i];
        k_f = k_f_s[94*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[20*npt+i];
        Kc = Kc_s[94*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;
        wdot[20*npt+i] += qdot;

        /*reaction 96: CH2O + O <=> HCO + OH */
        phi_f = sc[3*npt+i]*sc[11*npt+i];
        k_f = k_f_s[95*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[20*npt+i];
        Kc = Kc_s[95*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;
        wdot[20*npt+i] += qdot;

        /*reaction 97: CH3 + CH2O <=> CH4 + HCO */
        phi_f = sc[10*npt+i]*sc[11*npt+i];
        k_f = k_f_s[96*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[13*npt+i]*sc[20*npt+i];
        Kc = Kc_s[96*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[10*npt+i] -= qdot;
        wdot[11*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;
        wdot[20*npt+i] += qdot;

        /*reaction 98: CH2O + OH <=> HCO + H2O */
        phi_f = sc[4*npt+i]*sc[11*npt+i];
        k_f = k_f_s[97*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[20*npt+i];
        Kc = Kc_s[97*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;
        wdot[20*npt+i] += qdot;

        /*reaction 99: CH2O + CH <=> CH2CO + H */
        phi_f = sc[11*npt+i]*sc[19*npt+i];
        k_f = k_f_s[98*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[16*npt+i];
        Kc = Kc_s[98*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;
        wdot[16*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 100: CH2O + O2 <=> HCO + HO2 */
        phi_f = sc[5*npt+i]*sc[11*npt+i];
        k_f = k_f_s[99*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i]*sc[20*npt+i];
        Kc = Kc_s[99*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;
        wdot[20*npt+i] += qdot;
    }
}

void vcomp_wdot_101_150(int npt, amrex::Real *  wdot, amrex::Real *  mixture, amrex::Real *  sc,
		amrex::Real *  k_f_s, amrex::Real *  Kc_s,
		amrex::Real *  tc, amrex::Real *  invT, amrex::Real *  T)
{
    for (int i=0; i<npt; i++) {
        amrex::Real qdot, q_f, q_r, phi_f, phi_r, k_f, k_r, Kc;
        amrex::Real alpha;

        /*reaction 101: CH2O + HO2 <=> HCO + H2O2 */
        phi_f = sc[8*npt+i]*sc[11*npt+i];
        k_f = k_f_s[100*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[20*npt+i];
        Kc = Kc_s[100*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[11*npt+i] -= qdot;
        wdot[20*npt+i] += qdot;

        /*reaction 102: 2.000000 H + CO2 <=> H2 + CO2 */
        phi_f = pow(sc[1*npt+i], 2.000000)*sc[12*npt+i];
        k_f = k_f_s[101*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[12*npt+i];
        Kc = Kc_s[101*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= 2.000000 * qdot;
        wdot[2*npt+i] += qdot;
        wdot[12*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;

        /*reaction 103: SXCH2 + CO2 <=> TXCH2 + CO2 */
        phi_f = sc[12*npt+i]*sc[22*npt+i];
        k_f = k_f_s[102*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[12*npt+i]*sc[21*npt+i];
        Kc = Kc_s[102*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[12*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;
        wdot[21*npt+i] += qdot;
        wdot[22*npt+i] -= qdot;

        /*reaction 104: SXCH2 + CO2 <=> CH2O + CO */
        phi_f = sc[12*npt+i]*sc[22*npt+i];
        k_f = k_f_s[103*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[9*npt+i]*sc[11*npt+i];
        Kc = Kc_s[103*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[9*npt+i] += qdot;
        wdot[11*npt+i] += qdot;
        wdot[12*npt+i] -= qdot;
        wdot[22*npt+i] -= qdot;

        /*reaction 105: CH + CO2 <=> HCO + CO */
        phi_f = sc[12*npt+i]*sc[19*npt+i];
        k_f = k_f_s[104*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[9*npt+i]*sc[20*npt+i];
        Kc = Kc_s[104*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[9*npt+i] += qdot;
        wdot[12*npt+i] -= qdot;
        wdot[19*npt+i] -= qdot;
        wdot[20*npt+i] += qdot;

        /*reaction 106: C2H2 + O <=> TXCH2 + CO */
        phi_f = sc[3*npt+i]*sc[14*npt+i];
        k_f = k_f_s[105*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[9*npt+i]*sc[21*npt+i];
        Kc = Kc_s[105*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;
        wdot[14*npt+i] -= qdot;
        wdot[21*npt+i] += qdot;

        /*reaction 107: C2H2 + OH <=> CH3 + CO */
        phi_f = sc[4*npt+i]*sc[14*npt+i];
        k_f = k_f_s[106*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[9*npt+i]*sc[10*npt+i];
        Kc = Kc_s[106*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;
        wdot[10*npt+i] += qdot;
        wdot[14*npt+i] -= qdot;

        /*reaction 108: C2H2 + OH <=> CH2CO + H */
        phi_f = sc[4*npt+i]*sc[14*npt+i];
        k_f = k_f_s[107*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[16*npt+i];
        Kc = Kc_s[107*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[14*npt+i] -= qdot;
        wdot[16*npt+i] += qdot;

        /*reaction 109: C2H2 + O <=> HCCO + H */
        phi_f = sc[3*npt+i]*sc[14*npt+i];
        k_f = k_f_s[108*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[25*npt+i];
        Kc = Kc_s[108*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[14*npt+i] -= qdot;
        wdot[25*npt+i] += qdot;

        /*reaction 110: C2H3 + OH <=> C2H2 + H2O */
        phi_f = sc[4*npt+i]*sc[23*npt+i];
        k_f = k_f_s[109*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[14*npt+i];
        Kc = Kc_s[109*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[14*npt+i] += qdot;
        wdot[23*npt+i] -= qdot;

        /*reaction 111: C2H3 + O2 <=> CH2CHO + O */
        phi_f = sc[5*npt+i]*sc[23*npt+i];
        k_f = k_f_s[110*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[27*npt+i];
        Kc = Kc_s[110*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[23*npt+i] -= qdot;
        wdot[27*npt+i] += qdot;

        /*reaction 112: C2H3 + O <=> CH2CHO */
        phi_f = sc[3*npt+i]*sc[23*npt+i];
        k_f = k_f_s[111*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[27*npt+i];
        Kc = Kc_s[111*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[23*npt+i] -= qdot;
        wdot[27*npt+i] += qdot;

        /*reaction 113: C2H3 + H <=> C2H2 + H2 */
        phi_f = sc[1*npt+i]*sc[23*npt+i];
        k_f = k_f_s[112*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[14*npt+i];
        Kc = Kc_s[112*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[14*npt+i] += qdot;
        wdot[23*npt+i] -= qdot;

        /*reaction 114: C2H3 + CH3 <=> C2H2 + CH4 */
        phi_f = sc[10*npt+i]*sc[23*npt+i];
        k_f = k_f_s[113*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[13*npt+i]*sc[14*npt+i];
        Kc = Kc_s[113*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[10*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;
        wdot[14*npt+i] += qdot;
        wdot[23*npt+i] -= qdot;

        /*reaction 115: C2H3 + O2 <=> HCO + CH2O */
        phi_f = sc[5*npt+i]*sc[23*npt+i];
        k_f = k_f_s[114*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[11*npt+i]*sc[20*npt+i];
        Kc = Kc_s[114*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;
        wdot[20*npt+i] += qdot;
        wdot[23*npt+i] -= qdot;

        /*reaction 116: C2H3 + H2O2 <=> C2H4 + HO2 */
        phi_f = sc[6*npt+i]*sc[23*npt+i];
        k_f = k_f_s[115*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i]*sc[15*npt+i];
        Kc = Kc_s[115*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[15*npt+i] += qdot;
        wdot[23*npt+i] -= qdot;

        /*reaction 117: C2H3 + O2 <=> C2H2 + HO2 */
        phi_f = sc[5*npt+i]*sc[23*npt+i];
        k_f = k_f_s[116*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i]*sc[14*npt+i];
        Kc = Kc_s[116*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[14*npt+i] += qdot;
        wdot[23*npt+i] -= qdot;

        /*reaction 118: C2H4 + CH3 <=> C2H3 + CH4 */
        phi_f = sc[10*npt+i]*sc[15*npt+i];
        k_f = k_f_s[117*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[13*npt+i]*sc[23*npt+i];
        Kc = Kc_s[117*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[10*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;
        wdot[15*npt+i] -= qdot;
        wdot[23*npt+i] += qdot;

        /*reaction 119: C2H4 + O2 => CH3 + CO2 + H */
        phi_f = sc[5*npt+i]*sc[15*npt+i];
        k_f = k_f_s[118*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;
        wdot[12*npt+i] += qdot;
        wdot[15*npt+i] -= qdot;

        /*reaction 120: C2H4 + OH <=> C2H3 + H2O */
        phi_f = sc[4*npt+i]*sc[15*npt+i];
        k_f = k_f_s[119*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[23*npt+i];
        Kc = Kc_s[119*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[15*npt+i] -= qdot;
        wdot[23*npt+i] += qdot;

        /*reaction 121: C2H4 + OH <=> C2H5O */
        phi_f = sc[4*npt+i]*sc[15*npt+i];
        k_f = k_f_s[120*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[28*npt+i];
        Kc = Kc_s[120*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[15*npt+i] -= qdot;
        wdot[28*npt+i] += qdot;

        /*reaction 122: C2H4 + O <=> CH2CHO + H */
        phi_f = sc[3*npt+i]*sc[15*npt+i];
        k_f = k_f_s[121*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[27*npt+i];
        Kc = Kc_s[121*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[15*npt+i] -= qdot;
        wdot[27*npt+i] += qdot;

        /*reaction 123: C2H4 + O <=> CH3 + HCO */
        phi_f = sc[3*npt+i]*sc[15*npt+i];
        k_f = k_f_s[122*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[10*npt+i]*sc[20*npt+i];
        Kc = Kc_s[122*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;
        wdot[15*npt+i] -= qdot;
        wdot[20*npt+i] += qdot;

        /*reaction 124: C2H4 + O2 <=> C2H3 + HO2 */
        phi_f = sc[5*npt+i]*sc[15*npt+i];
        k_f = k_f_s[123*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i]*sc[23*npt+i];
        Kc = Kc_s[123*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[15*npt+i] -= qdot;
        wdot[23*npt+i] += qdot;

        /*reaction 125: C2H4 + H <=> C2H3 + H2 */
        phi_f = sc[1*npt+i]*sc[15*npt+i];
        k_f = k_f_s[124*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[23*npt+i];
        Kc = Kc_s[124*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[15*npt+i] -= qdot;
        wdot[23*npt+i] += qdot;

        /*reaction 126: C2H4 + O <=> TXCH2 + CH2O */
        phi_f = sc[3*npt+i]*sc[15*npt+i];
        k_f = k_f_s[125*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[11*npt+i]*sc[21*npt+i];
        Kc = Kc_s[125*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;
        wdot[15*npt+i] -= qdot;
        wdot[21*npt+i] += qdot;

        /*reaction 127: C2H5 + HO2 <=> C2H4 + H2O2 */
        phi_f = sc[8*npt+i]*sc[24*npt+i];
        k_f = k_f_s[126*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[15*npt+i];
        Kc = Kc_s[126*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;
        wdot[24*npt+i] -= qdot;

        /*reaction 128: C2H5 + HO2 <=> C2H5O + OH */
        phi_f = sc[8*npt+i]*sc[24*npt+i];
        k_f = k_f_s[127*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[28*npt+i];
        Kc = Kc_s[127*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[24*npt+i] -= qdot;
        wdot[28*npt+i] += qdot;

        /*reaction 129: C2H5 + O <=> C2H5O */
        phi_f = sc[3*npt+i]*sc[24*npt+i];
        k_f = k_f_s[128*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[28*npt+i];
        Kc = Kc_s[128*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[24*npt+i] -= qdot;
        wdot[28*npt+i] += qdot;

        /*reaction 130: C2H5 + H <=> C2H4 + H2 */
        phi_f = sc[1*npt+i]*sc[24*npt+i];
        k_f = k_f_s[129*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[15*npt+i];
        Kc = Kc_s[129*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[15*npt+i] += qdot;
        wdot[24*npt+i] -= qdot;

        /*reaction 131: C2H5 + O2 <=> C2H4 + HO2 */
        phi_f = sc[5*npt+i]*sc[24*npt+i];
        k_f = k_f_s[130*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i]*sc[15*npt+i];
        Kc = Kc_s[130*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[15*npt+i] += qdot;
        wdot[24*npt+i] -= qdot;

        /*reaction 132: C2H5 + HO2 <=> C2H6 + O2 */
        phi_f = sc[8*npt+i]*sc[24*npt+i];
        k_f = k_f_s[131*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[17*npt+i];
        Kc = Kc_s[131*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[17*npt+i] += qdot;
        wdot[24*npt+i] -= qdot;

        /*reaction 133: C2H5 + CH3 <=> C2H4 + CH4 */
        phi_f = sc[10*npt+i]*sc[24*npt+i];
        k_f = k_f_s[132*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[13*npt+i]*sc[15*npt+i];
        Kc = Kc_s[132*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[10*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;
        wdot[15*npt+i] += qdot;
        wdot[24*npt+i] -= qdot;

        /*reaction 134: C2H6 + SXCH2 <=> C2H5 + CH3 */
        phi_f = sc[17*npt+i]*sc[22*npt+i];
        k_f = k_f_s[133*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[10*npt+i]*sc[24*npt+i];
        Kc = Kc_s[133*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[10*npt+i] += qdot;
        wdot[17*npt+i] -= qdot;
        wdot[22*npt+i] -= qdot;
        wdot[24*npt+i] += qdot;

        /*reaction 135: C2H6 + CH3 <=> C2H5 + CH4 */
        phi_f = sc[10*npt+i]*sc[17*npt+i];
        k_f = k_f_s[134*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[13*npt+i]*sc[24*npt+i];
        Kc = Kc_s[134*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[10*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;
        wdot[17*npt+i] -= qdot;
        wdot[24*npt+i] += qdot;

        /*reaction 136: C2H6 + O <=> C2H5 + OH */
        phi_f = sc[3*npt+i]*sc[17*npt+i];
        k_f = k_f_s[135*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[24*npt+i];
        Kc = Kc_s[135*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[17*npt+i] -= qdot;
        wdot[24*npt+i] += qdot;

        /*reaction 137: C2H6 + HO2 <=> C2H5 + H2O2 */
        phi_f = sc[8*npt+i]*sc[17*npt+i];
        k_f = k_f_s[136*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[6*npt+i]*sc[24*npt+i];
        Kc = Kc_s[136*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[6*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[17*npt+i] -= qdot;
        wdot[24*npt+i] += qdot;

        /*reaction 138: C2H6 + H <=> C2H5 + H2 */
        phi_f = sc[1*npt+i]*sc[17*npt+i];
        k_f = k_f_s[137*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[24*npt+i];
        Kc = Kc_s[137*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[17*npt+i] -= qdot;
        wdot[24*npt+i] += qdot;

        /*reaction 139: C2H6 + OH <=> C2H5 + H2O */
        phi_f = sc[4*npt+i]*sc[17*npt+i];
        k_f = k_f_s[138*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[24*npt+i];
        Kc = Kc_s[138*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[17*npt+i] -= qdot;
        wdot[24*npt+i] += qdot;

        /*reaction 140: HCCO + O2 <=> OH + 2.000000 CO */
        phi_f = sc[5*npt+i]*sc[25*npt+i];
        k_f = k_f_s[139*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*pow(sc[9*npt+i], 2.000000);
        Kc = Kc_s[139*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[9*npt+i] += 2.000000 * qdot;
        wdot[25*npt+i] -= qdot;

        /*reaction 141: HCCO + O <=> H + 2.000000 CO */
        phi_f = sc[3*npt+i]*sc[25*npt+i];
        k_f = k_f_s[140*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*pow(sc[9*npt+i], 2.000000);
        Kc = Kc_s[140*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[9*npt+i] += 2.000000 * qdot;
        wdot[25*npt+i] -= qdot;

        /*reaction 142: HCCO + CH3 <=> C2H4 + CO */
        phi_f = sc[10*npt+i]*sc[25*npt+i];
        k_f = k_f_s[141*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[9*npt+i]*sc[15*npt+i];
        Kc = Kc_s[141*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[9*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;
        wdot[25*npt+i] -= qdot;

        /*reaction 143: HCCO + H <=> SXCH2 + CO */
        phi_f = sc[1*npt+i]*sc[25*npt+i];
        k_f = k_f_s[142*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[9*npt+i]*sc[22*npt+i];
        Kc = Kc_s[142*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;
        wdot[22*npt+i] += qdot;
        wdot[25*npt+i] -= qdot;

        /*reaction 144: CH2CO + H <=> CH3 + CO */
        phi_f = sc[1*npt+i]*sc[16*npt+i];
        k_f = k_f_s[143*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[9*npt+i]*sc[10*npt+i];
        Kc = Kc_s[143*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;
        wdot[10*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;

        /*reaction 145: CH2CO + TXCH2 <=> C2H4 + CO */
        phi_f = sc[16*npt+i]*sc[21*npt+i];
        k_f = k_f_s[144*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[9*npt+i]*sc[15*npt+i];
        Kc = Kc_s[144*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[9*npt+i] += qdot;
        wdot[15*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;
        wdot[21*npt+i] -= qdot;

        /*reaction 146: CH2CO + O <=> HCCO + OH */
        phi_f = sc[3*npt+i]*sc[16*npt+i];
        k_f = k_f_s[145*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[25*npt+i];
        Kc = Kc_s[145*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;
        wdot[25*npt+i] += qdot;

        /*reaction 147: CH2CO + CH3 <=> HCCO + CH4 */
        phi_f = sc[10*npt+i]*sc[16*npt+i];
        k_f = k_f_s[146*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[13*npt+i]*sc[25*npt+i];
        Kc = Kc_s[146*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[10*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;
        wdot[25*npt+i] += qdot;

        /*reaction 148: CH2CO + O <=> TXCH2 + CO2 */
        phi_f = sc[3*npt+i]*sc[16*npt+i];
        k_f = k_f_s[147*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[12*npt+i]*sc[21*npt+i];
        Kc = Kc_s[147*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;
        wdot[21*npt+i] += qdot;

        /*reaction 149: CH2CO + CH3 <=> C2H5 + CO */
        phi_f = sc[10*npt+i]*sc[16*npt+i];
        k_f = k_f_s[148*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[9*npt+i]*sc[24*npt+i];
        Kc = Kc_s[148*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[9*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;
        wdot[16*npt+i] -= qdot;
        wdot[24*npt+i] += qdot;

        /*reaction 150: CH2CO + OH <=> HCCO + H2O */
        phi_f = sc[4*npt+i]*sc[16*npt+i];
        k_f = k_f_s[149*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[25*npt+i];
        Kc = Kc_s[149*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;
        wdot[25*npt+i] += qdot;
    }
}

void vcomp_wdot_151_172(int npt, amrex::Real *  wdot, amrex::Real *  mixture, amrex::Real *  sc,
		amrex::Real *  k_f_s, amrex::Real *  Kc_s,
		amrex::Real *  tc, amrex::Real *  invT, amrex::Real *  T)
{
    for (int i=0; i<npt; i++) {
        amrex::Real qdot, q_f, q_r, phi_f, phi_r, k_f, k_r, Kc;
        amrex::Real alpha;

        /*reaction 151: CH2CO + H <=> HCCO + H2 */
        phi_f = sc[1*npt+i]*sc[16*npt+i];
        k_f = k_f_s[150*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[25*npt+i];
        Kc = Kc_s[150*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;
        wdot[25*npt+i] += qdot;

        /*reaction 152: CH2CO + TXCH2 <=> HCCO + CH3 */
        phi_f = sc[16*npt+i]*sc[21*npt+i];
        k_f = k_f_s[151*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[10*npt+i]*sc[25*npt+i];
        Kc = Kc_s[151*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[10*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;
        wdot[21*npt+i] -= qdot;
        wdot[25*npt+i] += qdot;

        /*reaction 153: CH2CHO + O <=> CH2O + HCO */
        phi_f = sc[3*npt+i]*sc[27*npt+i];
        k_f = k_f_s[152*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[11*npt+i]*sc[20*npt+i];
        Kc = Kc_s[152*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;
        wdot[20*npt+i] += qdot;
        wdot[27*npt+i] -= qdot;

        /*reaction 154: CH2CHO <=> CH2CO + H */
        phi_f = sc[27*npt+i];
        k_f = k_f_s[153*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[16*npt+i];
        Kc = Kc_s[153*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[16*npt+i] += qdot;
        wdot[27*npt+i] -= qdot;

        /*reaction 155: CH2CHO + OH <=> H2O + CH2CO */
        phi_f = sc[4*npt+i]*sc[27*npt+i];
        k_f = k_f_s[154*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[7*npt+i]*sc[16*npt+i];
        Kc = Kc_s[154*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[16*npt+i] += qdot;
        wdot[27*npt+i] -= qdot;

        /*reaction 156: CH2CHO + H <=> CH2CO + H2 */
        phi_f = sc[1*npt+i]*sc[27*npt+i];
        k_f = k_f_s[155*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[16*npt+i];
        Kc = Kc_s[155*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[16*npt+i] += qdot;
        wdot[27*npt+i] -= qdot;

        /*reaction 157: CH2CHO + O2 => OH + CO + CH2O */
        phi_f = sc[5*npt+i]*sc[27*npt+i];
        k_f = k_f_s[156*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[4*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;
        wdot[11*npt+i] += qdot;
        wdot[27*npt+i] -= qdot;

        /*reaction 158: CH2CHO <=> CH3 + CO */
        phi_f = sc[27*npt+i];
        k_f = k_f_s[157*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[9*npt+i]*sc[10*npt+i];
        Kc = Kc_s[157*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[9*npt+i] += qdot;
        wdot[10*npt+i] += qdot;
        wdot[27*npt+i] -= qdot;

        /*reaction 159: CH2CHO + O2 => OH + 2.000000 HCO */
        phi_f = sc[5*npt+i]*sc[27*npt+i];
        k_f = k_f_s[158*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[4*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[20*npt+i] += 2.000000 * qdot;
        wdot[27*npt+i] -= qdot;

        /*reaction 160: CH2CHO + H <=> CH3 + HCO */
        phi_f = sc[1*npt+i]*sc[27*npt+i];
        k_f = k_f_s[159*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[10*npt+i]*sc[20*npt+i];
        Kc = Kc_s[159*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;
        wdot[20*npt+i] += qdot;
        wdot[27*npt+i] -= qdot;

        /*reaction 161: CH3CHO + O => CH3 + CO + OH */
        phi_f = sc[3*npt+i]*sc[26*npt+i];
        k_f = k_f_s[160*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[9*npt+i] += qdot;
        wdot[10*npt+i] += qdot;
        wdot[26*npt+i] -= qdot;

        /*reaction 162: CH3CHO + O2 => CH3 + CO + HO2 */
        phi_f = sc[5*npt+i]*sc[26*npt+i];
        k_f = k_f_s[161*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[9*npt+i] += qdot;
        wdot[10*npt+i] += qdot;
        wdot[26*npt+i] -= qdot;

        /*reaction 163: CH3CHO + OH => CH3 + CO + H2O */
        phi_f = sc[4*npt+i]*sc[26*npt+i];
        k_f = k_f_s[162*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[4*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[9*npt+i] += qdot;
        wdot[10*npt+i] += qdot;
        wdot[26*npt+i] -= qdot;

        /*reaction 164: CH3CHO + H <=> CH2CHO + H2 */
        phi_f = sc[1*npt+i]*sc[26*npt+i];
        k_f = k_f_s[163*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[27*npt+i];
        Kc = Kc_s[163*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[26*npt+i] -= qdot;
        wdot[27*npt+i] += qdot;

        /*reaction 165: CH3CHO + H => CH3 + CO + H2 */
        phi_f = sc[1*npt+i]*sc[26*npt+i];
        k_f = k_f_s[164*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[9*npt+i] += qdot;
        wdot[10*npt+i] += qdot;
        wdot[26*npt+i] -= qdot;

        /*reaction 166: CH3CHO + O <=> CH2CHO + OH */
        phi_f = sc[3*npt+i]*sc[26*npt+i];
        k_f = k_f_s[165*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[27*npt+i];
        Kc = Kc_s[165*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[26*npt+i] -= qdot;
        wdot[27*npt+i] += qdot;

        /*reaction 167: CH3CHO + CH3 => CH3 + CO + CH4 */
        phi_f = sc[10*npt+i]*sc[26*npt+i];
        k_f = k_f_s[166*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[9*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;
        wdot[13*npt+i] += qdot;
        wdot[26*npt+i] -= qdot;

        /*reaction 168: CH3CHO + HO2 => CH3 + CO + H2O2 */
        phi_f = sc[8*npt+i]*sc[26*npt+i];
        k_f = k_f_s[167*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[6*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;
        wdot[10*npt+i] += qdot;
        wdot[26*npt+i] -= qdot;

        /*reaction 169: C2H5O <=> CH3 + CH2O */
        phi_f = sc[28*npt+i];
        k_f = k_f_s[168*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[10*npt+i]*sc[11*npt+i];
        Kc = Kc_s[168*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[10*npt+i] += qdot;
        wdot[11*npt+i] += qdot;
        wdot[28*npt+i] -= qdot;

        /*reaction 170: C2H5O <=> CH3CHO + H */
        phi_f = sc[28*npt+i];
        k_f = k_f_s[169*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[26*npt+i];
        Kc = Kc_s[169*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[26*npt+i] += qdot;
        wdot[28*npt+i] -= qdot;

        /*reaction 171: C2H5O + O2 <=> CH3CHO + HO2 */
        phi_f = sc[5*npt+i]*sc[28*npt+i];
        k_f = k_f_s[170*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[8*npt+i]*sc[26*npt+i];
        Kc = Kc_s[170*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[26*npt+i] += qdot;
        wdot[28*npt+i] -= qdot;

        /*reaction 172: SXCH2 + N2 <=> TXCH2 + N2 */
        phi_f = sc[0*npt+i]*sc[22*npt+i];
        k_f = k_f_s[171*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[21*npt+i];
        Kc = Kc_s[171*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[0*npt+i] += qdot;
        wdot[21*npt+i] += qdot;
        wdot[22*npt+i] -= qdot;
    }
}


static amrex::Real T_save = -1;
#ifdef _OPENMP
#pragma omp threadprivate(T_save)
#endif

static amrex::Real k_f_save[172];
#ifdef _OPENMP
#pragma omp threadprivate(k_f_save)
#endif

static amrex::Real Kc_save[172];
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

    amrex::Real qdot, q_f[172], q_r[172];
    amrex::Real sc_qss[1];
    comp_qfqr(q_f, q_r, sc, sc_qss, tc, invT);

    for (int i = 0; i < 29; ++i) {
        wdot[i] = 0.0;
    }

    qdot = q_f[0]-q_r[0];
    wdot[4] -= 2.000000 * qdot;
    wdot[6] += qdot;

    qdot = q_f[1]-q_r[1];
    wdot[1] -= qdot;
    wdot[5] -= qdot;
    wdot[8] += qdot;

    qdot = q_f[2]-q_r[2];
    wdot[2] -= qdot;
    wdot[10] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[3]-q_r[3];
    wdot[1] -= qdot;
    wdot[10] += qdot;
    wdot[21] -= qdot;

    qdot = q_f[4]-q_r[4];
    wdot[1] -= qdot;
    wdot[10] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[5]-q_r[5];
    wdot[9] -= qdot;
    wdot[16] += qdot;
    wdot[21] -= qdot;

    qdot = q_f[6]-q_r[6];
    wdot[2] -= qdot;
    wdot[9] -= qdot;
    wdot[11] += qdot;

    qdot = q_f[7]-q_r[7];
    wdot[9] -= qdot;
    wdot[19] -= qdot;
    wdot[25] += qdot;

    qdot = q_f[8]-q_r[8];
    wdot[3] -= qdot;
    wdot[9] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[9]-q_r[9];
    wdot[1] -= qdot;
    wdot[11] += qdot;
    wdot[20] -= qdot;

    qdot = q_f[10]-q_r[10];
    wdot[1] -= qdot;
    wdot[14] -= qdot;
    wdot[23] += qdot;

    qdot = q_f[11]-q_r[11];
    wdot[1] -= qdot;
    wdot[15] += qdot;
    wdot[23] -= qdot;

    qdot = q_f[12]-q_r[12];
    wdot[1] -= qdot;
    wdot[15] -= qdot;
    wdot[24] += qdot;

    qdot = q_f[13]-q_r[13];
    wdot[1] -= qdot;
    wdot[17] += qdot;
    wdot[24] -= qdot;

    qdot = q_f[14]-q_r[14];
    wdot[10] += 2.000000 * qdot;
    wdot[17] -= qdot;

    qdot = q_f[15]-q_r[15];
    wdot[1] -= 2.000000 * qdot;
    wdot[2] += qdot;

    qdot = q_f[16]-q_r[16];
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[17]-q_r[17];
    wdot[3] -= 2.000000 * qdot;
    wdot[5] += qdot;

    qdot = q_f[18]-q_r[18];
    wdot[1] -= qdot;
    wdot[4] -= qdot;
    wdot[7] += qdot;

    qdot = q_f[19]-q_r[19];
    wdot[1] += qdot;
    wdot[9] += qdot;
    wdot[20] -= qdot;

    qdot = q_f[20]-q_r[20];
    wdot[1] -= 2.000000 * qdot;
    wdot[2] -= qdot;
    wdot[2] += 2.000000 * qdot;

    qdot = q_f[21]-q_r[21];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[3] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[22]-q_r[22];
    wdot[3] += qdot;
    wdot[4] -= 2.000000 * qdot;
    wdot[7] += qdot;

    qdot = q_f[23]-q_r[23];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[4] -= qdot;
    wdot[7] += qdot;

    qdot = q_f[24]-q_r[24];
    wdot[1] -= 2.000000 * qdot;
    wdot[2] += qdot;
    wdot[7] -= qdot;
    wdot[7] += qdot;

    qdot = q_f[25]-q_r[25];
    wdot[1] -= qdot;
    wdot[3] += qdot;
    wdot[4] += qdot;
    wdot[5] -= qdot;

    qdot = q_f[26]-q_r[26];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[5] -= qdot;
    wdot[8] += qdot;

    qdot = q_f[27]-q_r[27];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[7] += qdot;
    wdot[8] -= qdot;

    qdot = q_f[28]-q_r[28];
    wdot[1] -= qdot;
    wdot[4] += 2.000000 * qdot;
    wdot[8] -= qdot;

    qdot = q_f[29]-q_r[29];
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[5] += qdot;
    wdot[8] -= qdot;

    qdot = q_f[30]-q_r[30];
    wdot[1] -= qdot;
    wdot[3] += qdot;
    wdot[7] += qdot;
    wdot[8] -= qdot;

    qdot = q_f[31]-q_r[31];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[7] += qdot;
    wdot[8] -= qdot;

    qdot = q_f[32]-q_r[32];
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[6] -= qdot;
    wdot[8] += qdot;

    qdot = q_f[33]-q_r[33];
    wdot[1] -= qdot;
    wdot[4] += qdot;
    wdot[6] -= qdot;
    wdot[7] += qdot;

    qdot = q_f[34]-q_r[34];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[6] -= qdot;
    wdot[8] += qdot;

    qdot = q_f[35]-q_r[35];
    wdot[4] -= qdot;
    wdot[6] -= qdot;
    wdot[7] += qdot;
    wdot[8] += qdot;

    qdot = q_f[36]-q_r[36];
    wdot[4] -= qdot;
    wdot[6] -= qdot;
    wdot[7] += qdot;
    wdot[8] += qdot;

    qdot = q_f[37]-q_r[37];
    wdot[3] += qdot;
    wdot[5] -= qdot;
    wdot[9] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[38]-q_r[38];
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[9] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[39]-q_r[39];
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[19] -= qdot;
    wdot[20] += qdot;

    qdot = q_f[40]-q_r[40];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[19] -= qdot;
    wdot[21] += qdot;

    qdot = q_f[41]-q_r[41];
    wdot[1] += qdot;
    wdot[3] -= qdot;
    wdot[9] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[42]-q_r[42];
    wdot[3] += qdot;
    wdot[5] -= qdot;
    wdot[19] -= qdot;
    wdot[20] += qdot;

    qdot = q_f[43]-q_r[43];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[18] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[44]-q_r[44];
    wdot[1] += qdot;
    wdot[7] -= qdot;
    wdot[11] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[45]-q_r[45];
    wdot[1] += qdot;
    wdot[4] += qdot;
    wdot[5] -= qdot;
    wdot[9] += qdot;
    wdot[21] -= qdot;

    qdot = q_f[46]-q_r[46];
    wdot[3] += qdot;
    wdot[5] -= qdot;
    wdot[11] += qdot;
    wdot[21] -= qdot;

    qdot = q_f[47]-q_r[47];
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[11] += qdot;
    wdot[21] -= qdot;

    qdot = q_f[48]-q_r[48];
    wdot[4] += qdot;
    wdot[8] -= qdot;
    wdot[11] += qdot;
    wdot[21] -= qdot;

    qdot = q_f[49]-q_r[49];
    wdot[1] += 2.000000 * qdot;
    wdot[5] -= qdot;
    wdot[12] += qdot;
    wdot[21] -= qdot;

    qdot = q_f[50]-q_r[50];
    wdot[4] -= qdot;
    wdot[7] += qdot;
    wdot[19] += qdot;
    wdot[21] -= qdot;

    qdot = q_f[51]-q_r[51];
    wdot[1] += qdot;
    wdot[3] -= qdot;
    wdot[20] += qdot;
    wdot[21] -= qdot;

    qdot = q_f[52]-q_r[52];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[10] += qdot;
    wdot[21] -= qdot;

    qdot = q_f[53]-q_r[53];
    wdot[7] -= qdot;
    wdot[7] += qdot;
    wdot[21] += qdot;
    wdot[22] -= qdot;

    qdot = q_f[54]-q_r[54];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[19] += qdot;
    wdot[22] -= qdot;

    qdot = q_f[55]-q_r[55];
    wdot[1] += qdot;
    wdot[4] += qdot;
    wdot[5] -= qdot;
    wdot[9] += qdot;
    wdot[22] -= qdot;

    qdot = q_f[56]-q_r[56];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[9] += qdot;
    wdot[22] -= qdot;

    qdot = q_f[57]-q_r[57];
    wdot[5] -= qdot;
    wdot[7] += qdot;
    wdot[9] += qdot;
    wdot[22] -= qdot;

    qdot = q_f[58]-q_r[58];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[10] += qdot;
    wdot[22] -= qdot;

    qdot = q_f[59]-q_r[59];
    wdot[1] += qdot;
    wdot[3] -= qdot;
    wdot[20] += qdot;
    wdot[22] -= qdot;

    qdot = q_f[60]-q_r[60];
    wdot[2] += qdot;
    wdot[7] -= qdot;
    wdot[11] += qdot;
    wdot[22] -= qdot;

    qdot = q_f[61]-q_r[61];
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[11] += qdot;
    wdot[22] -= qdot;

    qdot = q_f[62]-q_r[62];
    wdot[2] += qdot;
    wdot[4] -= qdot;
    wdot[10] -= qdot;
    wdot[11] += qdot;

    qdot = q_f[63]-q_r[63];
    wdot[6] -= qdot;
    wdot[8] += qdot;
    wdot[10] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[64]-q_r[64];
    wdot[4] += qdot;
    wdot[5] -= qdot;
    wdot[10] -= qdot;
    wdot[11] += qdot;

    qdot = q_f[65]-q_r[65];
    wdot[1] += qdot;
    wdot[10] -= qdot;
    wdot[19] -= qdot;
    wdot[23] += qdot;

    qdot = q_f[66]-q_r[66];
    wdot[1] += qdot;
    wdot[3] -= qdot;
    wdot[10] -= qdot;
    wdot[11] += qdot;

    qdot = q_f[67]-q_r[67];
    wdot[1] += qdot;
    wdot[10] -= qdot;
    wdot[14] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[68]-q_r[68];
    wdot[4] -= qdot;
    wdot[7] += qdot;
    wdot[10] -= qdot;
    wdot[21] += qdot;

    qdot = q_f[69]-q_r[69];
    wdot[1] += qdot;
    wdot[10] -= qdot;
    wdot[15] += qdot;
    wdot[22] -= qdot;

    qdot = q_f[70]-q_r[70];
    wdot[4] -= qdot;
    wdot[7] += qdot;
    wdot[10] -= qdot;
    wdot[22] += qdot;

    qdot = q_f[71]-q_r[71];
    wdot[1] += qdot;
    wdot[10] -= 2.000000 * qdot;
    wdot[24] += qdot;

    qdot = q_f[72]-q_r[72];
    wdot[5] += qdot;
    wdot[8] -= qdot;
    wdot[10] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[73]-q_r[73];
    wdot[1] += qdot;
    wdot[10] -= qdot;
    wdot[15] += qdot;
    wdot[21] -= qdot;

    qdot = q_f[74]-q_r[74];
    wdot[1] += qdot;
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[9] += qdot;
    wdot[10] -= qdot;

    qdot = q_f[75]-q_r[75];
    wdot[1] += qdot;
    wdot[13] -= qdot;
    wdot[15] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[76]-q_r[76];
    wdot[10] += 2.000000 * qdot;
    wdot[13] -= qdot;
    wdot[22] -= qdot;

    qdot = q_f[77]-q_r[77];
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[10] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[78]-q_r[78];
    wdot[4] -= qdot;
    wdot[7] += qdot;
    wdot[10] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[79]-q_r[79];
    wdot[10] += 2.000000 * qdot;
    wdot[13] -= qdot;
    wdot[21] -= qdot;

    qdot = q_f[80]-q_r[80];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[10] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[81]-q_r[81];
    wdot[9] -= qdot;
    wdot[9] += qdot;
    wdot[21] += qdot;
    wdot[22] -= qdot;

    qdot = q_f[82]-q_r[82];
    wdot[3] += qdot;
    wdot[5] -= qdot;
    wdot[9] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[83]-q_r[83];
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[9] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[84]-q_r[84];
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[9] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[85]-q_r[85];
    wdot[4] += qdot;
    wdot[8] -= qdot;
    wdot[9] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[86]-q_r[86];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[9] += qdot;
    wdot[20] -= qdot;

    qdot = q_f[87]-q_r[87];
    wdot[10] -= qdot;
    wdot[20] -= qdot;
    wdot[26] += qdot;

    qdot = q_f[88]-q_r[88];
    wdot[1] += qdot;
    wdot[7] -= qdot;
    wdot[7] += qdot;
    wdot[9] += qdot;
    wdot[20] -= qdot;

    qdot = q_f[89]-q_r[89];
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[9] += qdot;
    wdot[20] -= qdot;

    qdot = q_f[90]-q_r[90];
    wdot[4] -= qdot;
    wdot[7] += qdot;
    wdot[9] += qdot;
    wdot[20] -= qdot;

    qdot = q_f[91]-q_r[91];
    wdot[9] += qdot;
    wdot[10] -= qdot;
    wdot[13] += qdot;
    wdot[20] -= qdot;

    qdot = q_f[92]-q_r[92];
    wdot[1] += qdot;
    wdot[3] -= qdot;
    wdot[12] += qdot;
    wdot[20] -= qdot;

    qdot = q_f[93]-q_r[93];
    wdot[5] -= qdot;
    wdot[8] += qdot;
    wdot[9] += qdot;
    wdot[20] -= qdot;

    qdot = q_f[94]-q_r[94];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[11] -= qdot;
    wdot[20] += qdot;

    qdot = q_f[95]-q_r[95];
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[11] -= qdot;
    wdot[20] += qdot;

    qdot = q_f[96]-q_r[96];
    wdot[10] -= qdot;
    wdot[11] -= qdot;
    wdot[13] += qdot;
    wdot[20] += qdot;

    qdot = q_f[97]-q_r[97];
    wdot[4] -= qdot;
    wdot[7] += qdot;
    wdot[11] -= qdot;
    wdot[20] += qdot;

    qdot = q_f[98]-q_r[98];
    wdot[1] += qdot;
    wdot[11] -= qdot;
    wdot[16] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[99]-q_r[99];
    wdot[5] -= qdot;
    wdot[8] += qdot;
    wdot[11] -= qdot;
    wdot[20] += qdot;

    qdot = q_f[100]-q_r[100];
    wdot[6] += qdot;
    wdot[8] -= qdot;
    wdot[11] -= qdot;
    wdot[20] += qdot;

    qdot = q_f[101]-q_r[101];
    wdot[1] -= 2.000000 * qdot;
    wdot[2] += qdot;
    wdot[12] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[102]-q_r[102];
    wdot[12] -= qdot;
    wdot[12] += qdot;
    wdot[21] += qdot;
    wdot[22] -= qdot;

    qdot = q_f[103]-q_r[103];
    wdot[9] += qdot;
    wdot[11] += qdot;
    wdot[12] -= qdot;
    wdot[22] -= qdot;

    qdot = q_f[104]-q_r[104];
    wdot[9] += qdot;
    wdot[12] -= qdot;
    wdot[19] -= qdot;
    wdot[20] += qdot;

    qdot = q_f[105]-q_r[105];
    wdot[3] -= qdot;
    wdot[9] += qdot;
    wdot[14] -= qdot;
    wdot[21] += qdot;

    qdot = q_f[106]-q_r[106];
    wdot[4] -= qdot;
    wdot[9] += qdot;
    wdot[10] += qdot;
    wdot[14] -= qdot;

    qdot = q_f[107]-q_r[107];
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[14] -= qdot;
    wdot[16] += qdot;

    qdot = q_f[108]-q_r[108];
    wdot[1] += qdot;
    wdot[3] -= qdot;
    wdot[14] -= qdot;
    wdot[25] += qdot;

    qdot = q_f[109]-q_r[109];
    wdot[4] -= qdot;
    wdot[7] += qdot;
    wdot[14] += qdot;
    wdot[23] -= qdot;

    qdot = q_f[110]-q_r[110];
    wdot[3] += qdot;
    wdot[5] -= qdot;
    wdot[23] -= qdot;
    wdot[27] += qdot;

    qdot = q_f[111]-q_r[111];
    wdot[3] -= qdot;
    wdot[23] -= qdot;
    wdot[27] += qdot;

    qdot = q_f[112]-q_r[112];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[14] += qdot;
    wdot[23] -= qdot;

    qdot = q_f[113]-q_r[113];
    wdot[10] -= qdot;
    wdot[13] += qdot;
    wdot[14] += qdot;
    wdot[23] -= qdot;

    qdot = q_f[114]-q_r[114];
    wdot[5] -= qdot;
    wdot[11] += qdot;
    wdot[20] += qdot;
    wdot[23] -= qdot;

    qdot = q_f[115]-q_r[115];
    wdot[6] -= qdot;
    wdot[8] += qdot;
    wdot[15] += qdot;
    wdot[23] -= qdot;

    qdot = q_f[116]-q_r[116];
    wdot[5] -= qdot;
    wdot[8] += qdot;
    wdot[14] += qdot;
    wdot[23] -= qdot;

    qdot = q_f[117]-q_r[117];
    wdot[10] -= qdot;
    wdot[13] += qdot;
    wdot[15] -= qdot;
    wdot[23] += qdot;

    qdot = q_f[118]-q_r[118];
    wdot[1] += qdot;
    wdot[5] -= qdot;
    wdot[10] += qdot;
    wdot[12] += qdot;
    wdot[15] -= qdot;

    qdot = q_f[119]-q_r[119];
    wdot[4] -= qdot;
    wdot[7] += qdot;
    wdot[15] -= qdot;
    wdot[23] += qdot;

    qdot = q_f[120]-q_r[120];
    wdot[4] -= qdot;
    wdot[15] -= qdot;
    wdot[28] += qdot;

    qdot = q_f[121]-q_r[121];
    wdot[1] += qdot;
    wdot[3] -= qdot;
    wdot[15] -= qdot;
    wdot[27] += qdot;

    qdot = q_f[122]-q_r[122];
    wdot[3] -= qdot;
    wdot[10] += qdot;
    wdot[15] -= qdot;
    wdot[20] += qdot;

    qdot = q_f[123]-q_r[123];
    wdot[5] -= qdot;
    wdot[8] += qdot;
    wdot[15] -= qdot;
    wdot[23] += qdot;

    qdot = q_f[124]-q_r[124];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[15] -= qdot;
    wdot[23] += qdot;

    qdot = q_f[125]-q_r[125];
    wdot[3] -= qdot;
    wdot[11] += qdot;
    wdot[15] -= qdot;
    wdot[21] += qdot;

    qdot = q_f[126]-q_r[126];
    wdot[6] += qdot;
    wdot[8] -= qdot;
    wdot[15] += qdot;
    wdot[24] -= qdot;

    qdot = q_f[127]-q_r[127];
    wdot[4] += qdot;
    wdot[8] -= qdot;
    wdot[24] -= qdot;
    wdot[28] += qdot;

    qdot = q_f[128]-q_r[128];
    wdot[3] -= qdot;
    wdot[24] -= qdot;
    wdot[28] += qdot;

    qdot = q_f[129]-q_r[129];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[15] += qdot;
    wdot[24] -= qdot;

    qdot = q_f[130]-q_r[130];
    wdot[5] -= qdot;
    wdot[8] += qdot;
    wdot[15] += qdot;
    wdot[24] -= qdot;

    qdot = q_f[131]-q_r[131];
    wdot[5] += qdot;
    wdot[8] -= qdot;
    wdot[17] += qdot;
    wdot[24] -= qdot;

    qdot = q_f[132]-q_r[132];
    wdot[10] -= qdot;
    wdot[13] += qdot;
    wdot[15] += qdot;
    wdot[24] -= qdot;

    qdot = q_f[133]-q_r[133];
    wdot[10] += qdot;
    wdot[17] -= qdot;
    wdot[22] -= qdot;
    wdot[24] += qdot;

    qdot = q_f[134]-q_r[134];
    wdot[10] -= qdot;
    wdot[13] += qdot;
    wdot[17] -= qdot;
    wdot[24] += qdot;

    qdot = q_f[135]-q_r[135];
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[17] -= qdot;
    wdot[24] += qdot;

    qdot = q_f[136]-q_r[136];
    wdot[6] += qdot;
    wdot[8] -= qdot;
    wdot[17] -= qdot;
    wdot[24] += qdot;

    qdot = q_f[137]-q_r[137];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[17] -= qdot;
    wdot[24] += qdot;

    qdot = q_f[138]-q_r[138];
    wdot[4] -= qdot;
    wdot[7] += qdot;
    wdot[17] -= qdot;
    wdot[24] += qdot;

    qdot = q_f[139]-q_r[139];
    wdot[4] += qdot;
    wdot[5] -= qdot;
    wdot[9] += 2.000000 * qdot;
    wdot[25] -= qdot;

    qdot = q_f[140]-q_r[140];
    wdot[1] += qdot;
    wdot[3] -= qdot;
    wdot[9] += 2.000000 * qdot;
    wdot[25] -= qdot;

    qdot = q_f[141]-q_r[141];
    wdot[9] += qdot;
    wdot[10] -= qdot;
    wdot[15] += qdot;
    wdot[25] -= qdot;

    qdot = q_f[142]-q_r[142];
    wdot[1] -= qdot;
    wdot[9] += qdot;
    wdot[22] += qdot;
    wdot[25] -= qdot;

    qdot = q_f[143]-q_r[143];
    wdot[1] -= qdot;
    wdot[9] += qdot;
    wdot[10] += qdot;
    wdot[16] -= qdot;

    qdot = q_f[144]-q_r[144];
    wdot[9] += qdot;
    wdot[15] += qdot;
    wdot[16] -= qdot;
    wdot[21] -= qdot;

    qdot = q_f[145]-q_r[145];
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[16] -= qdot;
    wdot[25] += qdot;

    qdot = q_f[146]-q_r[146];
    wdot[10] -= qdot;
    wdot[13] += qdot;
    wdot[16] -= qdot;
    wdot[25] += qdot;

    qdot = q_f[147]-q_r[147];
    wdot[3] -= qdot;
    wdot[12] += qdot;
    wdot[16] -= qdot;
    wdot[21] += qdot;

    qdot = q_f[148]-q_r[148];
    wdot[9] += qdot;
    wdot[10] -= qdot;
    wdot[16] -= qdot;
    wdot[24] += qdot;

    qdot = q_f[149]-q_r[149];
    wdot[4] -= qdot;
    wdot[7] += qdot;
    wdot[16] -= qdot;
    wdot[25] += qdot;

    qdot = q_f[150]-q_r[150];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[16] -= qdot;
    wdot[25] += qdot;

    qdot = q_f[151]-q_r[151];
    wdot[10] += qdot;
    wdot[16] -= qdot;
    wdot[21] -= qdot;
    wdot[25] += qdot;

    qdot = q_f[152]-q_r[152];
    wdot[3] -= qdot;
    wdot[11] += qdot;
    wdot[20] += qdot;
    wdot[27] -= qdot;

    qdot = q_f[153]-q_r[153];
    wdot[1] += qdot;
    wdot[16] += qdot;
    wdot[27] -= qdot;

    qdot = q_f[154]-q_r[154];
    wdot[4] -= qdot;
    wdot[7] += qdot;
    wdot[16] += qdot;
    wdot[27] -= qdot;

    qdot = q_f[155]-q_r[155];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[16] += qdot;
    wdot[27] -= qdot;

    qdot = q_f[156]-q_r[156];
    wdot[4] += qdot;
    wdot[5] -= qdot;
    wdot[9] += qdot;
    wdot[11] += qdot;
    wdot[27] -= qdot;

    qdot = q_f[157]-q_r[157];
    wdot[9] += qdot;
    wdot[10] += qdot;
    wdot[27] -= qdot;

    qdot = q_f[158]-q_r[158];
    wdot[4] += qdot;
    wdot[5] -= qdot;
    wdot[20] += 2.000000 * qdot;
    wdot[27] -= qdot;

    qdot = q_f[159]-q_r[159];
    wdot[1] -= qdot;
    wdot[10] += qdot;
    wdot[20] += qdot;
    wdot[27] -= qdot;

    qdot = q_f[160]-q_r[160];
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[9] += qdot;
    wdot[10] += qdot;
    wdot[26] -= qdot;

    qdot = q_f[161]-q_r[161];
    wdot[5] -= qdot;
    wdot[8] += qdot;
    wdot[9] += qdot;
    wdot[10] += qdot;
    wdot[26] -= qdot;

    qdot = q_f[162]-q_r[162];
    wdot[4] -= qdot;
    wdot[7] += qdot;
    wdot[9] += qdot;
    wdot[10] += qdot;
    wdot[26] -= qdot;

    qdot = q_f[163]-q_r[163];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[26] -= qdot;
    wdot[27] += qdot;

    qdot = q_f[164]-q_r[164];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[9] += qdot;
    wdot[10] += qdot;
    wdot[26] -= qdot;

    qdot = q_f[165]-q_r[165];
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[26] -= qdot;
    wdot[27] += qdot;

    qdot = q_f[166]-q_r[166];
    wdot[9] += qdot;
    wdot[10] -= qdot;
    wdot[10] += qdot;
    wdot[13] += qdot;
    wdot[26] -= qdot;

    qdot = q_f[167]-q_r[167];
    wdot[6] += qdot;
    wdot[8] -= qdot;
    wdot[9] += qdot;
    wdot[10] += qdot;
    wdot[26] -= qdot;

    qdot = q_f[168]-q_r[168];
    wdot[10] += qdot;
    wdot[11] += qdot;
    wdot[28] -= qdot;

    qdot = q_f[169]-q_r[169];
    wdot[1] += qdot;
    wdot[26] += qdot;
    wdot[28] -= qdot;

    qdot = q_f[170]-q_r[170];
    wdot[5] -= qdot;
    wdot[8] += qdot;
    wdot[26] += qdot;
    wdot[28] -= qdot;

    qdot = q_f[171]-q_r[171];
    wdot[0] -= qdot;
    wdot[0] += qdot;
    wdot[21] += qdot;
    wdot[22] -= qdot;

    return;
}

void comp_k_f(amrex::Real *  tc, amrex::Real invT, amrex::Real *  k_f)
{
    for (int i=0; i<172; ++i) {
        k_f[i] = prefactor_units[i] * fwd_A[i]
                    * exp(fwd_beta[i] * tc[0] - activation_units[i] * fwd_Ea[i] * invT);
    };
    return;
}

void comp_Kc(amrex::Real *  tc, amrex::Real invT, amrex::Real *  Kc)
{
    /*compute the Gibbs free energy */
    amrex::Real g_RT[29];
    gibbs(g_RT, tc);

    Kc[0] = 2.000000*g_RT[4] - g_RT[6];
    Kc[1] = g_RT[1] + g_RT[5] - g_RT[8];
    Kc[2] = g_RT[2] - g_RT[10] + g_RT[19];
    Kc[3] = g_RT[1] - g_RT[10] + g_RT[21];
    Kc[4] = g_RT[1] + g_RT[10] - g_RT[13];
    Kc[5] = g_RT[9] - g_RT[16] + g_RT[21];
    Kc[6] = g_RT[2] + g_RT[9] - g_RT[11];
    Kc[7] = g_RT[9] + g_RT[19] - g_RT[25];
    Kc[8] = g_RT[3] + g_RT[9] - g_RT[12];
    Kc[9] = g_RT[1] - g_RT[11] + g_RT[20];
    Kc[10] = g_RT[1] + g_RT[14] - g_RT[23];
    Kc[11] = g_RT[1] - g_RT[15] + g_RT[23];
    Kc[12] = g_RT[1] + g_RT[15] - g_RT[24];
    Kc[13] = g_RT[1] - g_RT[17] + g_RT[24];
    Kc[14] = -2.000000*g_RT[10] + g_RT[17];
    Kc[15] = 2.000000*g_RT[1] - g_RT[2];
    Kc[16] = g_RT[1] + g_RT[3] - g_RT[4];
    Kc[17] = 2.000000*g_RT[3] - g_RT[5];
    Kc[18] = g_RT[1] + g_RT[4] - g_RT[7];
    Kc[19] = -g_RT[1] - g_RT[9] + g_RT[20];
    Kc[20] = 2.000000*g_RT[1] + g_RT[2] - 2.000000*g_RT[2];
    Kc[21] = -g_RT[1] + g_RT[2] + g_RT[3] - g_RT[4];
    Kc[22] = -g_RT[3] + 2.000000*g_RT[4] - g_RT[7];
    Kc[23] = -g_RT[1] + g_RT[2] + g_RT[4] - g_RT[7];
    Kc[24] = 2.000000*g_RT[1] - g_RT[2] + g_RT[7] - g_RT[7];
    Kc[25] = g_RT[1] - g_RT[3] - g_RT[4] + g_RT[5];
    Kc[26] = -g_RT[1] + g_RT[2] + g_RT[5] - g_RT[8];
    Kc[27] = g_RT[4] - g_RT[5] - g_RT[7] + g_RT[8];
    Kc[28] = g_RT[1] - 2.000000*g_RT[4] + g_RT[8];
    Kc[29] = g_RT[3] - g_RT[4] - g_RT[5] + g_RT[8];
    Kc[30] = g_RT[1] - g_RT[3] - g_RT[7] + g_RT[8];
    Kc[31] = g_RT[4] - g_RT[5] - g_RT[7] + g_RT[8];
    Kc[32] = g_RT[3] - g_RT[4] + g_RT[6] - g_RT[8];
    Kc[33] = g_RT[1] - g_RT[4] + g_RT[6] - g_RT[7];
    Kc[34] = g_RT[1] - g_RT[2] + g_RT[6] - g_RT[8];
    Kc[35] = g_RT[4] + g_RT[6] - g_RT[7] - g_RT[8];
    Kc[36] = g_RT[4] + g_RT[6] - g_RT[7] - g_RT[8];
    Kc[37] = -g_RT[3] + g_RT[5] - g_RT[9] + g_RT[18];
    Kc[38] = -g_RT[1] + g_RT[4] - g_RT[9] + g_RT[18];
    Kc[39] = -g_RT[1] + g_RT[4] + g_RT[19] - g_RT[20];
    Kc[40] = -g_RT[1] + g_RT[2] + g_RT[19] - g_RT[21];
    Kc[41] = -g_RT[1] + g_RT[3] - g_RT[9] + g_RT[19];
    Kc[42] = -g_RT[3] + g_RT[5] + g_RT[19] - g_RT[20];
    Kc[43] = g_RT[1] - g_RT[2] - g_RT[18] + g_RT[19];
    Kc[44] = -g_RT[1] + g_RT[7] - g_RT[11] + g_RT[19];
    Kc[45] = -g_RT[1] - g_RT[4] + g_RT[5] - g_RT[9] + g_RT[21];
    Kc[46] = -g_RT[3] + g_RT[5] - g_RT[11] + g_RT[21];
    Kc[47] = -g_RT[1] + g_RT[4] - g_RT[11] + g_RT[21];
    Kc[48] = -g_RT[4] + g_RT[8] - g_RT[11] + g_RT[21];
    Kc[49] = -2.000000*g_RT[1] + g_RT[5] - g_RT[12] + g_RT[21];
    Kc[50] = g_RT[4] - g_RT[7] - g_RT[19] + g_RT[21];
    Kc[51] = -g_RT[1] + g_RT[3] - g_RT[20] + g_RT[21];
    Kc[52] = -g_RT[1] + g_RT[2] - g_RT[10] + g_RT[21];
    Kc[53] = g_RT[7] - g_RT[7] - g_RT[21] + g_RT[22];
    Kc[54] = g_RT[1] - g_RT[2] - g_RT[19] + g_RT[22];
    Kc[55] = -g_RT[1] - g_RT[4] + g_RT[5] - g_RT[9] + g_RT[22];
    Kc[56] = -g_RT[2] + g_RT[3] - g_RT[9] + g_RT[22];
    Kc[57] = g_RT[5] - g_RT[7] - g_RT[9] + g_RT[22];
    Kc[58] = -g_RT[1] + g_RT[2] - g_RT[10] + g_RT[22];
    Kc[59] = -g_RT[1] + g_RT[3] - g_RT[20] + g_RT[22];
    Kc[60] = -g_RT[2] + g_RT[7] - g_RT[11] + g_RT[22];
    Kc[61] = -g_RT[1] + g_RT[4] - g_RT[11] + g_RT[22];
    Kc[62] = -g_RT[2] + g_RT[4] + g_RT[10] - g_RT[11];
    Kc[63] = g_RT[6] - g_RT[8] + g_RT[10] - g_RT[13];
    Kc[64] = -g_RT[4] + g_RT[5] + g_RT[10] - g_RT[11];
    Kc[65] = -g_RT[1] + g_RT[10] + g_RT[19] - g_RT[23];
    Kc[66] = -g_RT[1] + g_RT[3] + g_RT[10] - g_RT[11];
    Kc[67] = -g_RT[1] + g_RT[10] - g_RT[14] + g_RT[18];
    Kc[68] = g_RT[4] - g_RT[7] + g_RT[10] - g_RT[21];
    Kc[69] = -g_RT[1] + g_RT[10] - g_RT[15] + g_RT[22];
    Kc[70] = g_RT[4] - g_RT[7] + g_RT[10] - g_RT[22];
    Kc[71] = -g_RT[1] + 2.000000*g_RT[10] - g_RT[24];
    Kc[72] = -g_RT[5] + g_RT[8] + g_RT[10] - g_RT[13];
    Kc[73] = -g_RT[1] + g_RT[10] - g_RT[15] + g_RT[21];
    Kc[74] = -g_RT[1] - g_RT[2] + g_RT[3] - g_RT[9] + g_RT[10];
    Kc[75] = -g_RT[1] + g_RT[13] - g_RT[15] + g_RT[19];
    Kc[76] = -2.000000*g_RT[10] + g_RT[13] + g_RT[22];
    Kc[77] = g_RT[3] - g_RT[4] - g_RT[10] + g_RT[13];
    Kc[78] = g_RT[4] - g_RT[7] - g_RT[10] + g_RT[13];
    Kc[79] = -2.000000*g_RT[10] + g_RT[13] + g_RT[21];
    Kc[80] = g_RT[1] - g_RT[2] - g_RT[10] + g_RT[13];
    Kc[81] = g_RT[9] - g_RT[9] - g_RT[21] + g_RT[22];
    Kc[82] = -g_RT[3] + g_RT[5] + g_RT[9] - g_RT[12];
    Kc[83] = -g_RT[1] + g_RT[4] + g_RT[9] - g_RT[12];
    Kc[84] = -g_RT[1] + g_RT[4] + g_RT[9] - g_RT[12];
    Kc[85] = -g_RT[4] + g_RT[8] + g_RT[9] - g_RT[12];
    Kc[86] = g_RT[1] - g_RT[2] - g_RT[9] + g_RT[20];
    Kc[87] = g_RT[10] + g_RT[20] - g_RT[26];
    Kc[88] = -g_RT[1] + g_RT[7] - g_RT[7] - g_RT[9] + g_RT[20];
    Kc[89] = g_RT[3] - g_RT[4] - g_RT[9] + g_RT[20];
    Kc[90] = g_RT[4] - g_RT[7] - g_RT[9] + g_RT[20];
    Kc[91] = -g_RT[9] + g_RT[10] - g_RT[13] + g_RT[20];
    Kc[92] = -g_RT[1] + g_RT[3] - g_RT[12] + g_RT[20];
    Kc[93] = g_RT[5] - g_RT[8] - g_RT[9] + g_RT[20];
    Kc[94] = g_RT[1] - g_RT[2] + g_RT[11] - g_RT[20];
    Kc[95] = g_RT[3] - g_RT[4] + g_RT[11] - g_RT[20];
    Kc[96] = g_RT[10] + g_RT[11] - g_RT[13] - g_RT[20];
    Kc[97] = g_RT[4] - g_RT[7] + g_RT[11] - g_RT[20];
    Kc[98] = -g_RT[1] + g_RT[11] - g_RT[16] + g_RT[19];
    Kc[99] = g_RT[5] - g_RT[8] + g_RT[11] - g_RT[20];
    Kc[100] = -g_RT[6] + g_RT[8] + g_RT[11] - g_RT[20];
    Kc[101] = 2.000000*g_RT[1] - g_RT[2] + g_RT[12] - g_RT[12];
    Kc[102] = g_RT[12] - g_RT[12] - g_RT[21] + g_RT[22];
    Kc[103] = -g_RT[9] - g_RT[11] + g_RT[12] + g_RT[22];
    Kc[104] = -g_RT[9] + g_RT[12] + g_RT[19] - g_RT[20];
    Kc[105] = g_RT[3] - g_RT[9] + g_RT[14] - g_RT[21];
    Kc[106] = g_RT[4] - g_RT[9] - g_RT[10] + g_RT[14];
    Kc[107] = -g_RT[1] + g_RT[4] + g_RT[14] - g_RT[16];
    Kc[108] = -g_RT[1] + g_RT[3] + g_RT[14] - g_RT[25];
    Kc[109] = g_RT[4] - g_RT[7] - g_RT[14] + g_RT[23];
    Kc[110] = -g_RT[3] + g_RT[5] + g_RT[23] - g_RT[27];
    Kc[111] = g_RT[3] + g_RT[23] - g_RT[27];
    Kc[112] = g_RT[1] - g_RT[2] - g_RT[14] + g_RT[23];
    Kc[113] = g_RT[10] - g_RT[13] - g_RT[14] + g_RT[23];
    Kc[114] = g_RT[5] - g_RT[11] - g_RT[20] + g_RT[23];
    Kc[115] = g_RT[6] - g_RT[8] - g_RT[15] + g_RT[23];
    Kc[116] = g_RT[5] - g_RT[8] - g_RT[14] + g_RT[23];
    Kc[117] = g_RT[10] - g_RT[13] + g_RT[15] - g_RT[23];
    Kc[118] = -g_RT[1] + g_RT[5] - g_RT[10] - g_RT[12] + g_RT[15];
    Kc[119] = g_RT[4] - g_RT[7] + g_RT[15] - g_RT[23];
    Kc[120] = g_RT[4] + g_RT[15] - g_RT[28];
    Kc[121] = -g_RT[1] + g_RT[3] + g_RT[15] - g_RT[27];
    Kc[122] = g_RT[3] - g_RT[10] + g_RT[15] - g_RT[20];
    Kc[123] = g_RT[5] - g_RT[8] + g_RT[15] - g_RT[23];
    Kc[124] = g_RT[1] - g_RT[2] + g_RT[15] - g_RT[23];
    Kc[125] = g_RT[3] - g_RT[11] + g_RT[15] - g_RT[21];
    Kc[126] = -g_RT[6] + g_RT[8] - g_RT[15] + g_RT[24];
    Kc[127] = -g_RT[4] + g_RT[8] + g_RT[24] - g_RT[28];
    Kc[128] = g_RT[3] + g_RT[24] - g_RT[28];
    Kc[129] = g_RT[1] - g_RT[2] - g_RT[15] + g_RT[24];
    Kc[130] = g_RT[5] - g_RT[8] - g_RT[15] + g_RT[24];
    Kc[131] = -g_RT[5] + g_RT[8] - g_RT[17] + g_RT[24];
    Kc[132] = g_RT[10] - g_RT[13] - g_RT[15] + g_RT[24];
    Kc[133] = -g_RT[10] + g_RT[17] + g_RT[22] - g_RT[24];
    Kc[134] = g_RT[10] - g_RT[13] + g_RT[17] - g_RT[24];
    Kc[135] = g_RT[3] - g_RT[4] + g_RT[17] - g_RT[24];
    Kc[136] = -g_RT[6] + g_RT[8] + g_RT[17] - g_RT[24];
    Kc[137] = g_RT[1] - g_RT[2] + g_RT[17] - g_RT[24];
    Kc[138] = g_RT[4] - g_RT[7] + g_RT[17] - g_RT[24];
    Kc[139] = -g_RT[4] + g_RT[5] - 2.000000*g_RT[9] + g_RT[25];
    Kc[140] = -g_RT[1] + g_RT[3] - 2.000000*g_RT[9] + g_RT[25];
    Kc[141] = -g_RT[9] + g_RT[10] - g_RT[15] + g_RT[25];
    Kc[142] = g_RT[1] - g_RT[9] - g_RT[22] + g_RT[25];
    Kc[143] = g_RT[1] - g_RT[9] - g_RT[10] + g_RT[16];
    Kc[144] = -g_RT[9] - g_RT[15] + g_RT[16] + g_RT[21];
    Kc[145] = g_RT[3] - g_RT[4] + g_RT[16] - g_RT[25];
    Kc[146] = g_RT[10] - g_RT[13] + g_RT[16] - g_RT[25];
    Kc[147] = g_RT[3] - g_RT[12] + g_RT[16] - g_RT[21];
    Kc[148] = -g_RT[9] + g_RT[10] + g_RT[16] - g_RT[24];
    Kc[149] = g_RT[4] - g_RT[7] + g_RT[16] - g_RT[25];
    Kc[150] = g_RT[1] - g_RT[2] + g_RT[16] - g_RT[25];
    Kc[151] = -g_RT[10] + g_RT[16] + g_RT[21] - g_RT[25];
    Kc[152] = g_RT[3] - g_RT[11] - g_RT[20] + g_RT[27];
    Kc[153] = -g_RT[1] - g_RT[16] + g_RT[27];
    Kc[154] = g_RT[4] - g_RT[7] - g_RT[16] + g_RT[27];
    Kc[155] = g_RT[1] - g_RT[2] - g_RT[16] + g_RT[27];
    Kc[156] = -g_RT[4] + g_RT[5] - g_RT[9] - g_RT[11] + g_RT[27];
    Kc[157] = -g_RT[9] - g_RT[10] + g_RT[27];
    Kc[158] = -g_RT[4] + g_RT[5] - 2.000000*g_RT[20] + g_RT[27];
    Kc[159] = g_RT[1] - g_RT[10] - g_RT[20] + g_RT[27];
    Kc[160] = g_RT[3] - g_RT[4] - g_RT[9] - g_RT[10] + g_RT[26];
    Kc[161] = g_RT[5] - g_RT[8] - g_RT[9] - g_RT[10] + g_RT[26];
    Kc[162] = g_RT[4] - g_RT[7] - g_RT[9] - g_RT[10] + g_RT[26];
    Kc[163] = g_RT[1] - g_RT[2] + g_RT[26] - g_RT[27];
    Kc[164] = g_RT[1] - g_RT[2] - g_RT[9] - g_RT[10] + g_RT[26];
    Kc[165] = g_RT[3] - g_RT[4] + g_RT[26] - g_RT[27];
    Kc[166] = -g_RT[9] + g_RT[10] - g_RT[10] - g_RT[13] + g_RT[26];
    Kc[167] = -g_RT[6] + g_RT[8] - g_RT[9] - g_RT[10] + g_RT[26];
    Kc[168] = -g_RT[10] - g_RT[11] + g_RT[28];
    Kc[169] = -g_RT[1] - g_RT[26] + g_RT[28];
    Kc[170] = g_RT[5] - g_RT[8] - g_RT[26] + g_RT[28];
    Kc[171] = g_RT[0] - g_RT[0] - g_RT[21] + g_RT[22];

    for (int i=0; i<172; ++i) {
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
    Kc[14] *= refC;
    Kc[15] *= refCinv;
    Kc[16] *= refCinv;
    Kc[17] *= refCinv;
    Kc[18] *= refCinv;
    Kc[19] *= refC;
    Kc[20] *= refCinv;
    Kc[24] *= refCinv;
    Kc[45] *= refC;
    Kc[49] *= refC;
    Kc[55] *= refC;
    Kc[74] *= refC;
    Kc[87] *= refCinv;
    Kc[88] *= refC;
    Kc[101] *= refCinv;
    Kc[111] *= refCinv;
    Kc[118] *= refC;
    Kc[120] *= refCinv;
    Kc[128] *= refCinv;
    Kc[139] *= refC;
    Kc[140] *= refC;
    Kc[153] *= refC;
    Kc[156] *= refC;
    Kc[157] *= refC;
    Kc[158] *= refC;
    Kc[160] *= refC;
    Kc[161] *= refC;
    Kc[162] *= refC;
    Kc[164] *= refC;
    Kc[166] *= refC;
    Kc[167] *= refC;
    Kc[168] *= refC;
    Kc[169] *= refC;

    return;
}

void comp_qfqr(amrex::Real *  qf, amrex::Real *  qr, amrex::Real *  sc, amrex::Real * qss_sc, amrex::Real *  tc, amrex::Real invT)
{

    /*reaction 1: 2.000000 OH (+M) <=> H2O2 (+M) */
    qf[0] = pow(sc[4], 2.000000);
    qr[0] = sc[6];

    /*reaction 2: H + O2 (+M) <=> HO2 (+M) */
    qf[1] = sc[1]*sc[5];
    qr[1] = sc[8];

    /*reaction 3: CH + H2 (+M) <=> CH3 (+M) */
    qf[2] = sc[2]*sc[19];
    qr[2] = sc[10];

    /*reaction 4: TXCH2 + H (+M) <=> CH3 (+M) */
    qf[3] = sc[1]*sc[21];
    qr[3] = sc[10];

    /*reaction 5: CH3 + H (+M) <=> CH4 (+M) */
    qf[4] = sc[1]*sc[10];
    qr[4] = sc[13];

    /*reaction 6: TXCH2 + CO (+M) <=> CH2CO (+M) */
    qf[5] = sc[9]*sc[21];
    qr[5] = sc[16];

    /*reaction 7: CO + H2 (+M) <=> CH2O (+M) */
    qf[6] = sc[2]*sc[9];
    qr[6] = sc[11];

    /*reaction 8: CH + CO (+M) <=> HCCO (+M) */
    qf[7] = sc[9]*sc[19];
    qr[7] = sc[25];

    /*reaction 9: CO + O (+M) <=> CO2 (+M) */
    qf[8] = sc[3]*sc[9];
    qr[8] = sc[12];

    /*reaction 10: HCO + H (+M) <=> CH2O (+M) */
    qf[9] = sc[1]*sc[20];
    qr[9] = sc[11];

    /*reaction 11: C2H2 + H (+M) <=> C2H3 (+M) */
    qf[10] = sc[1]*sc[14];
    qr[10] = sc[23];

    /*reaction 12: C2H3 + H (+M) <=> C2H4 (+M) */
    qf[11] = sc[1]*sc[23];
    qr[11] = sc[15];

    /*reaction 13: C2H4 + H (+M) <=> C2H5 (+M) */
    qf[12] = sc[1]*sc[15];
    qr[12] = sc[24];

    /*reaction 14: C2H5 + H (+M) <=> C2H6 (+M) */
    qf[13] = sc[1]*sc[24];
    qr[13] = sc[17];

    /*reaction 15: C2H6 (+M) <=> 2.000000 CH3 (+M) */
    qf[14] = sc[17];
    qr[14] = pow(sc[10], 2.000000);

    /*reaction 16: 2.000000 H + M <=> H2 + M */
    qf[15] = pow(sc[1], 2.000000);
    qr[15] = sc[2];

    /*reaction 17: H + O + M <=> OH + M */
    qf[16] = sc[1]*sc[3];
    qr[16] = sc[4];

    /*reaction 18: 2.000000 O + M <=> O2 + M */
    qf[17] = pow(sc[3], 2.000000);
    qr[17] = sc[5];

    /*reaction 19: H + OH + M <=> H2O + M */
    qf[18] = sc[1]*sc[4];
    qr[18] = sc[7];

    /*reaction 20: HCO + M <=> CO + H + M */
    qf[19] = sc[20];
    qr[19] = sc[1]*sc[9];

    /*reaction 21: 2.000000 H + H2 <=> 2.000000 H2 */
    qf[20] = pow(sc[1], 2.000000)*sc[2];
    qr[20] = pow(sc[2], 2.000000);

    /*reaction 22: O + H2 <=> H + OH */
    qf[21] = sc[2]*sc[3];
    qr[21] = sc[1]*sc[4];

    /*reaction 23: 2.000000 OH <=> O + H2O */
    qf[22] = pow(sc[4], 2.000000);
    qr[22] = sc[3]*sc[7];

    /*reaction 24: OH + H2 <=> H + H2O */
    qf[23] = sc[2]*sc[4];
    qr[23] = sc[1]*sc[7];

    /*reaction 25: 2.000000 H + H2O <=> H2 + H2O */
    qf[24] = pow(sc[1], 2.000000)*sc[7];
    qr[24] = sc[2]*sc[7];

    /*reaction 26: H + O2 <=> O + OH */
    qf[25] = sc[1]*sc[5];
    qr[25] = sc[3]*sc[4];

    /*reaction 27: H2 + O2 <=> HO2 + H */
    qf[26] = sc[2]*sc[5];
    qr[26] = sc[1]*sc[8];

    /*reaction 28: HO2 + OH <=> H2O + O2 */
    qf[27] = sc[4]*sc[8];
    qr[27] = sc[5]*sc[7];

    /*reaction 29: HO2 + H <=> 2.000000 OH */
    qf[28] = sc[1]*sc[8];
    qr[28] = pow(sc[4], 2.000000);

    /*reaction 30: HO2 + O <=> OH + O2 */
    qf[29] = sc[3]*sc[8];
    qr[29] = sc[4]*sc[5];

    /*reaction 31: HO2 + H <=> O + H2O */
    qf[30] = sc[1]*sc[8];
    qr[30] = sc[3]*sc[7];

    /*reaction 32: HO2 + OH <=> H2O + O2 */
    qf[31] = sc[4]*sc[8];
    qr[31] = sc[5]*sc[7];

    /*reaction 33: H2O2 + O <=> HO2 + OH */
    qf[32] = sc[3]*sc[6];
    qr[32] = sc[4]*sc[8];

    /*reaction 34: H2O2 + H <=> H2O + OH */
    qf[33] = sc[1]*sc[6];
    qr[33] = sc[4]*sc[7];

    /*reaction 35: H2O2 + H <=> HO2 + H2 */
    qf[34] = sc[1]*sc[6];
    qr[34] = sc[2]*sc[8];

    /*reaction 36: H2O2 + OH <=> HO2 + H2O */
    qf[35] = sc[4]*sc[6];
    qr[35] = sc[7]*sc[8];

    /*reaction 37: H2O2 + OH <=> HO2 + H2O */
    qf[36] = sc[4]*sc[6];
    qr[36] = sc[7]*sc[8];

    /*reaction 38: C + O2 <=> CO + O */
    qf[37] = sc[5]*sc[18];
    qr[37] = sc[3]*sc[9];

    /*reaction 39: C + OH <=> CO + H */
    qf[38] = sc[4]*sc[18];
    qr[38] = sc[1]*sc[9];

    /*reaction 40: CH + OH <=> HCO + H */
    qf[39] = sc[4]*sc[19];
    qr[39] = sc[1]*sc[20];

    /*reaction 41: CH + H2 <=> TXCH2 + H */
    qf[40] = sc[2]*sc[19];
    qr[40] = sc[1]*sc[21];

    /*reaction 42: CH + O <=> CO + H */
    qf[41] = sc[3]*sc[19];
    qr[41] = sc[1]*sc[9];

    /*reaction 43: CH + O2 <=> HCO + O */
    qf[42] = sc[5]*sc[19];
    qr[42] = sc[3]*sc[20];

    /*reaction 44: CH + H <=> C + H2 */
    qf[43] = sc[1]*sc[19];
    qr[43] = sc[2]*sc[18];

    /*reaction 45: CH + H2O <=> CH2O + H */
    qf[44] = sc[7]*sc[19];
    qr[44] = sc[1]*sc[11];

    /*reaction 46: TXCH2 + O2 => OH + H + CO */
    qf[45] = sc[5]*sc[21];
    qr[45] = 0.0;

    /*reaction 47: TXCH2 + O2 <=> CH2O + O */
    qf[46] = sc[5]*sc[21];
    qr[46] = sc[3]*sc[11];

    /*reaction 48: TXCH2 + OH <=> CH2O + H */
    qf[47] = sc[4]*sc[21];
    qr[47] = sc[1]*sc[11];

    /*reaction 49: TXCH2 + HO2 <=> CH2O + OH */
    qf[48] = sc[8]*sc[21];
    qr[48] = sc[4]*sc[11];

    /*reaction 50: TXCH2 + O2 => CO2 + 2.000000 H */
    qf[49] = sc[5]*sc[21];
    qr[49] = 0.0;

    /*reaction 51: TXCH2 + OH <=> CH + H2O */
    qf[50] = sc[4]*sc[21];
    qr[50] = sc[7]*sc[19];

    /*reaction 52: TXCH2 + O <=> HCO + H */
    qf[51] = sc[3]*sc[21];
    qr[51] = sc[1]*sc[20];

    /*reaction 53: TXCH2 + H2 <=> H + CH3 */
    qf[52] = sc[2]*sc[21];
    qr[52] = sc[1]*sc[10];

    /*reaction 54: SXCH2 + H2O <=> TXCH2 + H2O */
    qf[53] = sc[7]*sc[22];
    qr[53] = sc[7]*sc[21];

    /*reaction 55: SXCH2 + H <=> CH + H2 */
    qf[54] = sc[1]*sc[22];
    qr[54] = sc[2]*sc[19];

    /*reaction 56: SXCH2 + O2 <=> H + OH + CO */
    qf[55] = sc[5]*sc[22];
    qr[55] = sc[1]*sc[4]*sc[9];

    /*reaction 57: SXCH2 + O <=> CO + H2 */
    qf[56] = sc[3]*sc[22];
    qr[56] = sc[2]*sc[9];

    /*reaction 58: SXCH2 + O2 <=> CO + H2O */
    qf[57] = sc[5]*sc[22];
    qr[57] = sc[7]*sc[9];

    /*reaction 59: SXCH2 + H2 <=> CH3 + H */
    qf[58] = sc[2]*sc[22];
    qr[58] = sc[1]*sc[10];

    /*reaction 60: SXCH2 + O <=> HCO + H */
    qf[59] = sc[3]*sc[22];
    qr[59] = sc[1]*sc[20];

    /*reaction 61: SXCH2 + H2O => H2 + CH2O */
    qf[60] = sc[7]*sc[22];
    qr[60] = 0.0;

    /*reaction 62: SXCH2 + OH <=> CH2O + H */
    qf[61] = sc[4]*sc[22];
    qr[61] = sc[1]*sc[11];

    /*reaction 63: CH3 + OH => H2 + CH2O */
    qf[62] = sc[4]*sc[10];
    qr[62] = 0.0;

    /*reaction 64: CH3 + H2O2 <=> CH4 + HO2 */
    qf[63] = sc[6]*sc[10];
    qr[63] = sc[8]*sc[13];

    /*reaction 65: CH3 + O2 <=> CH2O + OH */
    qf[64] = sc[5]*sc[10];
    qr[64] = sc[4]*sc[11];

    /*reaction 66: CH3 + CH <=> C2H3 + H */
    qf[65] = sc[10]*sc[19];
    qr[65] = sc[1]*sc[23];

    /*reaction 67: CH3 + O <=> CH2O + H */
    qf[66] = sc[3]*sc[10];
    qr[66] = sc[1]*sc[11];

    /*reaction 68: CH3 + C <=> C2H2 + H */
    qf[67] = sc[10]*sc[18];
    qr[67] = sc[1]*sc[14];

    /*reaction 69: CH3 + OH <=> TXCH2 + H2O */
    qf[68] = sc[4]*sc[10];
    qr[68] = sc[7]*sc[21];

    /*reaction 70: CH3 + SXCH2 <=> C2H4 + H */
    qf[69] = sc[10]*sc[22];
    qr[69] = sc[1]*sc[15];

    /*reaction 71: CH3 + OH <=> SXCH2 + H2O */
    qf[70] = sc[4]*sc[10];
    qr[70] = sc[7]*sc[22];

    /*reaction 72: 2.000000 CH3 <=> C2H5 + H */
    qf[71] = pow(sc[10], 2.000000);
    qr[71] = sc[1]*sc[24];

    /*reaction 73: CH3 + HO2 <=> CH4 + O2 */
    qf[72] = sc[8]*sc[10];
    qr[72] = sc[5]*sc[13];

    /*reaction 74: CH3 + TXCH2 <=> C2H4 + H */
    qf[73] = sc[10]*sc[21];
    qr[73] = sc[1]*sc[15];

    /*reaction 75: CH3 + O => H + H2 + CO */
    qf[74] = sc[3]*sc[10];
    qr[74] = 0.0;

    /*reaction 76: CH4 + CH <=> C2H4 + H */
    qf[75] = sc[13]*sc[19];
    qr[75] = sc[1]*sc[15];

    /*reaction 77: CH4 + SXCH2 <=> 2.000000 CH3 */
    qf[76] = sc[13]*sc[22];
    qr[76] = pow(sc[10], 2.000000);

    /*reaction 78: CH4 + O <=> CH3 + OH */
    qf[77] = sc[3]*sc[13];
    qr[77] = sc[4]*sc[10];

    /*reaction 79: CH4 + OH <=> CH3 + H2O */
    qf[78] = sc[4]*sc[13];
    qr[78] = sc[7]*sc[10];

    /*reaction 80: CH4 + TXCH2 <=> 2.000000 CH3 */
    qf[79] = sc[13]*sc[21];
    qr[79] = pow(sc[10], 2.000000);

    /*reaction 81: CH4 + H <=> CH3 + H2 */
    qf[80] = sc[1]*sc[13];
    qr[80] = sc[2]*sc[10];

    /*reaction 82: SXCH2 + CO <=> TXCH2 + CO */
    qf[81] = sc[9]*sc[22];
    qr[81] = sc[9]*sc[21];

    /*reaction 83: CO + O2 <=> CO2 + O */
    qf[82] = sc[5]*sc[9];
    qr[82] = sc[3]*sc[12];

    /*reaction 84: CO + OH <=> CO2 + H */
    qf[83] = sc[4]*sc[9];
    qr[83] = sc[1]*sc[12];

    /*reaction 85: CO + OH <=> CO2 + H */
    qf[84] = sc[4]*sc[9];
    qr[84] = sc[1]*sc[12];

    /*reaction 86: CO + HO2 <=> CO2 + OH */
    qf[85] = sc[8]*sc[9];
    qr[85] = sc[4]*sc[12];

    /*reaction 87: HCO + H <=> CO + H2 */
    qf[86] = sc[1]*sc[20];
    qr[86] = sc[2]*sc[9];

    /*reaction 88: CH3 + HCO <=> CH3CHO */
    qf[87] = sc[10]*sc[20];
    qr[87] = sc[26];

    /*reaction 89: HCO + H2O <=> CO + H + H2O */
    qf[88] = sc[7]*sc[20];
    qr[88] = sc[1]*sc[7]*sc[9];

    /*reaction 90: HCO + O <=> CO + OH */
    qf[89] = sc[3]*sc[20];
    qr[89] = sc[4]*sc[9];

    /*reaction 91: HCO + OH <=> CO + H2O */
    qf[90] = sc[4]*sc[20];
    qr[90] = sc[7]*sc[9];

    /*reaction 92: CH3 + HCO <=> CH4 + CO */
    qf[91] = sc[10]*sc[20];
    qr[91] = sc[9]*sc[13];

    /*reaction 93: HCO + O <=> CO2 + H */
    qf[92] = sc[3]*sc[20];
    qr[92] = sc[1]*sc[12];

    /*reaction 94: HCO + O2 <=> CO + HO2 */
    qf[93] = sc[5]*sc[20];
    qr[93] = sc[8]*sc[9];

    /*reaction 95: CH2O + H <=> HCO + H2 */
    qf[94] = sc[1]*sc[11];
    qr[94] = sc[2]*sc[20];

    /*reaction 96: CH2O + O <=> HCO + OH */
    qf[95] = sc[3]*sc[11];
    qr[95] = sc[4]*sc[20];

    /*reaction 97: CH3 + CH2O <=> CH4 + HCO */
    qf[96] = sc[10]*sc[11];
    qr[96] = sc[13]*sc[20];

    /*reaction 98: CH2O + OH <=> HCO + H2O */
    qf[97] = sc[4]*sc[11];
    qr[97] = sc[7]*sc[20];

    /*reaction 99: CH2O + CH <=> CH2CO + H */
    qf[98] = sc[11]*sc[19];
    qr[98] = sc[1]*sc[16];

    /*reaction 100: CH2O + O2 <=> HCO + HO2 */
    qf[99] = sc[5]*sc[11];
    qr[99] = sc[8]*sc[20];

    /*reaction 101: CH2O + HO2 <=> HCO + H2O2 */
    qf[100] = sc[8]*sc[11];
    qr[100] = sc[6]*sc[20];

    /*reaction 102: 2.000000 H + CO2 <=> H2 + CO2 */
    qf[101] = pow(sc[1], 2.000000)*sc[12];
    qr[101] = sc[2]*sc[12];

    /*reaction 103: SXCH2 + CO2 <=> TXCH2 + CO2 */
    qf[102] = sc[12]*sc[22];
    qr[102] = sc[12]*sc[21];

    /*reaction 104: SXCH2 + CO2 <=> CH2O + CO */
    qf[103] = sc[12]*sc[22];
    qr[103] = sc[9]*sc[11];

    /*reaction 105: CH + CO2 <=> HCO + CO */
    qf[104] = sc[12]*sc[19];
    qr[104] = sc[9]*sc[20];

    /*reaction 106: C2H2 + O <=> TXCH2 + CO */
    qf[105] = sc[3]*sc[14];
    qr[105] = sc[9]*sc[21];

    /*reaction 107: C2H2 + OH <=> CH3 + CO */
    qf[106] = sc[4]*sc[14];
    qr[106] = sc[9]*sc[10];

    /*reaction 108: C2H2 + OH <=> CH2CO + H */
    qf[107] = sc[4]*sc[14];
    qr[107] = sc[1]*sc[16];

    /*reaction 109: C2H2 + O <=> HCCO + H */
    qf[108] = sc[3]*sc[14];
    qr[108] = sc[1]*sc[25];

    /*reaction 110: C2H3 + OH <=> C2H2 + H2O */
    qf[109] = sc[4]*sc[23];
    qr[109] = sc[7]*sc[14];

    /*reaction 111: C2H3 + O2 <=> CH2CHO + O */
    qf[110] = sc[5]*sc[23];
    qr[110] = sc[3]*sc[27];

    /*reaction 112: C2H3 + O <=> CH2CHO */
    qf[111] = sc[3]*sc[23];
    qr[111] = sc[27];

    /*reaction 113: C2H3 + H <=> C2H2 + H2 */
    qf[112] = sc[1]*sc[23];
    qr[112] = sc[2]*sc[14];

    /*reaction 114: C2H3 + CH3 <=> C2H2 + CH4 */
    qf[113] = sc[10]*sc[23];
    qr[113] = sc[13]*sc[14];

    /*reaction 115: C2H3 + O2 <=> HCO + CH2O */
    qf[114] = sc[5]*sc[23];
    qr[114] = sc[11]*sc[20];

    /*reaction 116: C2H3 + H2O2 <=> C2H4 + HO2 */
    qf[115] = sc[6]*sc[23];
    qr[115] = sc[8]*sc[15];

    /*reaction 117: C2H3 + O2 <=> C2H2 + HO2 */
    qf[116] = sc[5]*sc[23];
    qr[116] = sc[8]*sc[14];

    /*reaction 118: C2H4 + CH3 <=> C2H3 + CH4 */
    qf[117] = sc[10]*sc[15];
    qr[117] = sc[13]*sc[23];

    /*reaction 119: C2H4 + O2 => CH3 + CO2 + H */
    qf[118] = sc[5]*sc[15];
    qr[118] = 0.0;

    /*reaction 120: C2H4 + OH <=> C2H3 + H2O */
    qf[119] = sc[4]*sc[15];
    qr[119] = sc[7]*sc[23];

    /*reaction 121: C2H4 + OH <=> C2H5O */
    qf[120] = sc[4]*sc[15];
    qr[120] = sc[28];

    /*reaction 122: C2H4 + O <=> CH2CHO + H */
    qf[121] = sc[3]*sc[15];
    qr[121] = sc[1]*sc[27];

    /*reaction 123: C2H4 + O <=> CH3 + HCO */
    qf[122] = sc[3]*sc[15];
    qr[122] = sc[10]*sc[20];

    /*reaction 124: C2H4 + O2 <=> C2H3 + HO2 */
    qf[123] = sc[5]*sc[15];
    qr[123] = sc[8]*sc[23];

    /*reaction 125: C2H4 + H <=> C2H3 + H2 */
    qf[124] = sc[1]*sc[15];
    qr[124] = sc[2]*sc[23];

    /*reaction 126: C2H4 + O <=> TXCH2 + CH2O */
    qf[125] = sc[3]*sc[15];
    qr[125] = sc[11]*sc[21];

    /*reaction 127: C2H5 + HO2 <=> C2H4 + H2O2 */
    qf[126] = sc[8]*sc[24];
    qr[126] = sc[6]*sc[15];

    /*reaction 128: C2H5 + HO2 <=> C2H5O + OH */
    qf[127] = sc[8]*sc[24];
    qr[127] = sc[4]*sc[28];

    /*reaction 129: C2H5 + O <=> C2H5O */
    qf[128] = sc[3]*sc[24];
    qr[128] = sc[28];

    /*reaction 130: C2H5 + H <=> C2H4 + H2 */
    qf[129] = sc[1]*sc[24];
    qr[129] = sc[2]*sc[15];

    /*reaction 131: C2H5 + O2 <=> C2H4 + HO2 */
    qf[130] = sc[5]*sc[24];
    qr[130] = sc[8]*sc[15];

    /*reaction 132: C2H5 + HO2 <=> C2H6 + O2 */
    qf[131] = sc[8]*sc[24];
    qr[131] = sc[5]*sc[17];

    /*reaction 133: C2H5 + CH3 <=> C2H4 + CH4 */
    qf[132] = sc[10]*sc[24];
    qr[132] = sc[13]*sc[15];

    /*reaction 134: C2H6 + SXCH2 <=> C2H5 + CH3 */
    qf[133] = sc[17]*sc[22];
    qr[133] = sc[10]*sc[24];

    /*reaction 135: C2H6 + CH3 <=> C2H5 + CH4 */
    qf[134] = sc[10]*sc[17];
    qr[134] = sc[13]*sc[24];

    /*reaction 136: C2H6 + O <=> C2H5 + OH */
    qf[135] = sc[3]*sc[17];
    qr[135] = sc[4]*sc[24];

    /*reaction 137: C2H6 + HO2 <=> C2H5 + H2O2 */
    qf[136] = sc[8]*sc[17];
    qr[136] = sc[6]*sc[24];

    /*reaction 138: C2H6 + H <=> C2H5 + H2 */
    qf[137] = sc[1]*sc[17];
    qr[137] = sc[2]*sc[24];

    /*reaction 139: C2H6 + OH <=> C2H5 + H2O */
    qf[138] = sc[4]*sc[17];
    qr[138] = sc[7]*sc[24];

    /*reaction 140: HCCO + O2 <=> OH + 2.000000 CO */
    qf[139] = sc[5]*sc[25];
    qr[139] = sc[4]*pow(sc[9], 2.000000);

    /*reaction 141: HCCO + O <=> H + 2.000000 CO */
    qf[140] = sc[3]*sc[25];
    qr[140] = sc[1]*pow(sc[9], 2.000000);

    /*reaction 142: HCCO + CH3 <=> C2H4 + CO */
    qf[141] = sc[10]*sc[25];
    qr[141] = sc[9]*sc[15];

    /*reaction 143: HCCO + H <=> SXCH2 + CO */
    qf[142] = sc[1]*sc[25];
    qr[142] = sc[9]*sc[22];

    /*reaction 144: CH2CO + H <=> CH3 + CO */
    qf[143] = sc[1]*sc[16];
    qr[143] = sc[9]*sc[10];

    /*reaction 145: CH2CO + TXCH2 <=> C2H4 + CO */
    qf[144] = sc[16]*sc[21];
    qr[144] = sc[9]*sc[15];

    /*reaction 146: CH2CO + O <=> HCCO + OH */
    qf[145] = sc[3]*sc[16];
    qr[145] = sc[4]*sc[25];

    /*reaction 147: CH2CO + CH3 <=> HCCO + CH4 */
    qf[146] = sc[10]*sc[16];
    qr[146] = sc[13]*sc[25];

    /*reaction 148: CH2CO + O <=> TXCH2 + CO2 */
    qf[147] = sc[3]*sc[16];
    qr[147] = sc[12]*sc[21];

    /*reaction 149: CH2CO + CH3 <=> C2H5 + CO */
    qf[148] = sc[10]*sc[16];
    qr[148] = sc[9]*sc[24];

    /*reaction 150: CH2CO + OH <=> HCCO + H2O */
    qf[149] = sc[4]*sc[16];
    qr[149] = sc[7]*sc[25];

    /*reaction 151: CH2CO + H <=> HCCO + H2 */
    qf[150] = sc[1]*sc[16];
    qr[150] = sc[2]*sc[25];

    /*reaction 152: CH2CO + TXCH2 <=> HCCO + CH3 */
    qf[151] = sc[16]*sc[21];
    qr[151] = sc[10]*sc[25];

    /*reaction 153: CH2CHO + O <=> CH2O + HCO */
    qf[152] = sc[3]*sc[27];
    qr[152] = sc[11]*sc[20];

    /*reaction 154: CH2CHO <=> CH2CO + H */
    qf[153] = sc[27];
    qr[153] = sc[1]*sc[16];

    /*reaction 155: CH2CHO + OH <=> H2O + CH2CO */
    qf[154] = sc[4]*sc[27];
    qr[154] = sc[7]*sc[16];

    /*reaction 156: CH2CHO + H <=> CH2CO + H2 */
    qf[155] = sc[1]*sc[27];
    qr[155] = sc[2]*sc[16];

    /*reaction 157: CH2CHO + O2 => OH + CO + CH2O */
    qf[156] = sc[5]*sc[27];
    qr[156] = 0.0;

    /*reaction 158: CH2CHO <=> CH3 + CO */
    qf[157] = sc[27];
    qr[157] = sc[9]*sc[10];

    /*reaction 159: CH2CHO + O2 => OH + 2.000000 HCO */
    qf[158] = sc[5]*sc[27];
    qr[158] = 0.0;

    /*reaction 160: CH2CHO + H <=> CH3 + HCO */
    qf[159] = sc[1]*sc[27];
    qr[159] = sc[10]*sc[20];

    /*reaction 161: CH3CHO + O => CH3 + CO + OH */
    qf[160] = sc[3]*sc[26];
    qr[160] = 0.0;

    /*reaction 162: CH3CHO + O2 => CH3 + CO + HO2 */
    qf[161] = sc[5]*sc[26];
    qr[161] = 0.0;

    /*reaction 163: CH3CHO + OH => CH3 + CO + H2O */
    qf[162] = sc[4]*sc[26];
    qr[162] = 0.0;

    /*reaction 164: CH3CHO + H <=> CH2CHO + H2 */
    qf[163] = sc[1]*sc[26];
    qr[163] = sc[2]*sc[27];

    /*reaction 165: CH3CHO + H => CH3 + CO + H2 */
    qf[164] = sc[1]*sc[26];
    qr[164] = 0.0;

    /*reaction 166: CH3CHO + O <=> CH2CHO + OH */
    qf[165] = sc[3]*sc[26];
    qr[165] = sc[4]*sc[27];

    /*reaction 167: CH3CHO + CH3 => CH3 + CO + CH4 */
    qf[166] = sc[10]*sc[26];
    qr[166] = 0.0;

    /*reaction 168: CH3CHO + HO2 => CH3 + CO + H2O2 */
    qf[167] = sc[8]*sc[26];
    qr[167] = 0.0;

    /*reaction 169: C2H5O <=> CH3 + CH2O */
    qf[168] = sc[28];
    qr[168] = sc[10]*sc[11];

    /*reaction 170: C2H5O <=> CH3CHO + H */
    qf[169] = sc[28];
    qr[169] = sc[1]*sc[26];

    /*reaction 171: C2H5O + O2 <=> CH3CHO + HO2 */
    qf[170] = sc[5]*sc[28];
    qr[170] = sc[8]*sc[26];

    /*reaction 172: SXCH2 + N2 <=> TXCH2 + N2 */
    qf[171] = sc[0]*sc[22];
    qr[171] = sc[0]*sc[21];

    amrex::Real T = tc[1];

    /*compute the mixture concentration */
    amrex::Real mixture = 0.0;
    for (int i = 0; i < 29; ++i) {
        mixture += sc[i];
    }

    amrex::Real Corr[172];
    for (int i = 0; i < 172; ++i) {
        Corr[i] = 1.0;
    }

    /* troe */
    {
        amrex::Real alpha[15];
        alpha[0] = mixture + (TB[0][0] - 1)*sc[2] + (TB[0][1] - 1)*sc[7] + (TB[0][2] - 1)*sc[18] + (TB[0][3] - 1)*sc[9] + (TB[0][4] - 1)*sc[19] + (TB[0][5] - 1)*sc[20] + (TB[0][6] - 1)*sc[21] + (TB[0][7] - 1)*sc[12] + (TB[0][8] - 1)*sc[22] + (TB[0][9] - 1)*sc[13] + (TB[0][10] - 1)*sc[23] + (TB[0][11] - 1)*sc[24] + (TB[0][12] - 1)*sc[25] + (TB[0][13] - 1)*sc[26] + (TB[0][14] - 1)*sc[27] + (TB[0][15] - 1)*sc[28] + (TB[0][16] - 1)*sc[17];
        alpha[1] = mixture + (TB[1][0] - 1)*sc[2] + (TB[1][1] - 1)*sc[5] + (TB[1][2] - 1)*sc[7] + (TB[1][3] - 1)*sc[18] + (TB[1][4] - 1)*sc[9] + (TB[1][5] - 1)*sc[19] + (TB[1][6] - 1)*sc[20] + (TB[1][7] - 1)*sc[21] + (TB[1][8] - 1)*sc[12] + (TB[1][9] - 1)*sc[22] + (TB[1][10] - 1)*sc[23] + (TB[1][11] - 1)*sc[24] + (TB[1][12] - 1)*sc[25] + (TB[1][13] - 1)*sc[26] + (TB[1][14] - 1)*sc[27] + (TB[1][15] - 1)*sc[28];
        alpha[2] = mixture + (TB[2][0] - 1)*sc[2] + (TB[2][1] - 1)*sc[7] + (TB[2][2] - 1)*sc[18] + (TB[2][3] - 1)*sc[9] + (TB[2][4] - 1)*sc[19] + (TB[2][5] - 1)*sc[20] + (TB[2][6] - 1)*sc[21] + (TB[2][7] - 1)*sc[12] + (TB[2][8] - 1)*sc[22] + (TB[2][9] - 1)*sc[13] + (TB[2][10] - 1)*sc[23] + (TB[2][11] - 1)*sc[24] + (TB[2][12] - 1)*sc[25] + (TB[2][13] - 1)*sc[26] + (TB[2][14] - 1)*sc[27] + (TB[2][15] - 1)*sc[28] + (TB[2][16] - 1)*sc[17];
        alpha[3] = mixture + (TB[3][0] - 1)*sc[2] + (TB[3][1] - 1)*sc[7] + (TB[3][2] - 1)*sc[18] + (TB[3][3] - 1)*sc[9] + (TB[3][4] - 1)*sc[19] + (TB[3][5] - 1)*sc[20] + (TB[3][6] - 1)*sc[21] + (TB[3][7] - 1)*sc[12] + (TB[3][8] - 1)*sc[22] + (TB[3][9] - 1)*sc[13] + (TB[3][10] - 1)*sc[23] + (TB[3][11] - 1)*sc[24] + (TB[3][12] - 1)*sc[25] + (TB[3][13] - 1)*sc[26] + (TB[3][14] - 1)*sc[27] + (TB[3][15] - 1)*sc[28] + (TB[3][16] - 1)*sc[17];
        alpha[4] = mixture + (TB[4][0] - 1)*sc[2] + (TB[4][1] - 1)*sc[7] + (TB[4][2] - 1)*sc[18] + (TB[4][3] - 1)*sc[9] + (TB[4][4] - 1)*sc[19] + (TB[4][5] - 1)*sc[20] + (TB[4][6] - 1)*sc[21] + (TB[4][7] - 1)*sc[12] + (TB[4][8] - 1)*sc[22] + (TB[4][9] - 1)*sc[13] + (TB[4][10] - 1)*sc[23] + (TB[4][11] - 1)*sc[24] + (TB[4][12] - 1)*sc[25] + (TB[4][13] - 1)*sc[26] + (TB[4][14] - 1)*sc[27] + (TB[4][15] - 1)*sc[28] + (TB[4][16] - 1)*sc[17];
        alpha[5] = mixture + (TB[5][0] - 1)*sc[2] + (TB[5][1] - 1)*sc[7] + (TB[5][2] - 1)*sc[18] + (TB[5][3] - 1)*sc[9] + (TB[5][4] - 1)*sc[19] + (TB[5][5] - 1)*sc[20] + (TB[5][6] - 1)*sc[21] + (TB[5][7] - 1)*sc[12] + (TB[5][8] - 1)*sc[22] + (TB[5][9] - 1)*sc[13] + (TB[5][10] - 1)*sc[23] + (TB[5][11] - 1)*sc[24] + (TB[5][12] - 1)*sc[25] + (TB[5][13] - 1)*sc[26] + (TB[5][14] - 1)*sc[27] + (TB[5][15] - 1)*sc[28] + (TB[5][16] - 1)*sc[17];
        alpha[6] = mixture + (TB[6][0] - 1)*sc[2] + (TB[6][1] - 1)*sc[7] + (TB[6][2] - 1)*sc[18] + (TB[6][3] - 1)*sc[9] + (TB[6][4] - 1)*sc[19] + (TB[6][5] - 1)*sc[20] + (TB[6][6] - 1)*sc[21] + (TB[6][7] - 1)*sc[12] + (TB[6][8] - 1)*sc[22] + (TB[6][9] - 1)*sc[13] + (TB[6][10] - 1)*sc[23] + (TB[6][11] - 1)*sc[24] + (TB[6][12] - 1)*sc[25] + (TB[6][13] - 1)*sc[26] + (TB[6][14] - 1)*sc[27] + (TB[6][15] - 1)*sc[28] + (TB[6][16] - 1)*sc[17];
        alpha[7] = mixture + (TB[7][0] - 1)*sc[2] + (TB[7][1] - 1)*sc[7] + (TB[7][2] - 1)*sc[18] + (TB[7][3] - 1)*sc[9] + (TB[7][4] - 1)*sc[19] + (TB[7][5] - 1)*sc[20] + (TB[7][6] - 1)*sc[21] + (TB[7][7] - 1)*sc[12] + (TB[7][8] - 1)*sc[22] + (TB[7][9] - 1)*sc[13] + (TB[7][10] - 1)*sc[23] + (TB[7][11] - 1)*sc[24] + (TB[7][12] - 1)*sc[25] + (TB[7][13] - 1)*sc[26] + (TB[7][14] - 1)*sc[27] + (TB[7][15] - 1)*sc[28] + (TB[7][16] - 1)*sc[17];
        alpha[8] = mixture + (TB[8][0] - 1)*sc[2] + (TB[8][1] - 1)*sc[7] + (TB[8][2] - 1)*sc[18] + (TB[8][3] - 1)*sc[9] + (TB[8][4] - 1)*sc[19] + (TB[8][5] - 1)*sc[20] + (TB[8][6] - 1)*sc[21] + (TB[8][7] - 1)*sc[12] + (TB[8][8] - 1)*sc[22] + (TB[8][9] - 1)*sc[13] + (TB[8][10] - 1)*sc[23] + (TB[8][11] - 1)*sc[24] + (TB[8][12] - 1)*sc[25] + (TB[8][13] - 1)*sc[26] + (TB[8][14] - 1)*sc[27] + (TB[8][15] - 1)*sc[28] + (TB[8][16] - 1)*sc[17];
        alpha[9] = mixture + (TB[9][0] - 1)*sc[2] + (TB[9][1] - 1)*sc[7] + (TB[9][2] - 1)*sc[18] + (TB[9][3] - 1)*sc[9] + (TB[9][4] - 1)*sc[19] + (TB[9][5] - 1)*sc[20] + (TB[9][6] - 1)*sc[21] + (TB[9][7] - 1)*sc[12] + (TB[9][8] - 1)*sc[22] + (TB[9][9] - 1)*sc[13] + (TB[9][10] - 1)*sc[23] + (TB[9][11] - 1)*sc[24] + (TB[9][12] - 1)*sc[25] + (TB[9][13] - 1)*sc[26] + (TB[9][14] - 1)*sc[27] + (TB[9][15] - 1)*sc[28] + (TB[9][16] - 1)*sc[17];
        alpha[10] = mixture + (TB[10][0] - 1)*sc[2] + (TB[10][1] - 1)*sc[7] + (TB[10][2] - 1)*sc[18] + (TB[10][3] - 1)*sc[9] + (TB[10][4] - 1)*sc[19] + (TB[10][5] - 1)*sc[20] + (TB[10][6] - 1)*sc[21] + (TB[10][7] - 1)*sc[12] + (TB[10][8] - 1)*sc[22] + (TB[10][9] - 1)*sc[13] + (TB[10][10] - 1)*sc[23] + (TB[10][11] - 1)*sc[24] + (TB[10][12] - 1)*sc[25] + (TB[10][13] - 1)*sc[26] + (TB[10][14] - 1)*sc[27] + (TB[10][15] - 1)*sc[28] + (TB[10][16] - 1)*sc[17];
        alpha[11] = mixture + (TB[11][0] - 1)*sc[2] + (TB[11][1] - 1)*sc[7] + (TB[11][2] - 1)*sc[18] + (TB[11][3] - 1)*sc[9] + (TB[11][4] - 1)*sc[19] + (TB[11][5] - 1)*sc[20] + (TB[11][6] - 1)*sc[21] + (TB[11][7] - 1)*sc[12] + (TB[11][8] - 1)*sc[22] + (TB[11][9] - 1)*sc[13] + (TB[11][10] - 1)*sc[23] + (TB[11][11] - 1)*sc[24] + (TB[11][12] - 1)*sc[25] + (TB[11][13] - 1)*sc[26] + (TB[11][14] - 1)*sc[27] + (TB[11][15] - 1)*sc[28] + (TB[11][16] - 1)*sc[17];
        alpha[12] = mixture + (TB[12][0] - 1)*sc[2] + (TB[12][1] - 1)*sc[7] + (TB[12][2] - 1)*sc[18] + (TB[12][3] - 1)*sc[9] + (TB[12][4] - 1)*sc[19] + (TB[12][5] - 1)*sc[20] + (TB[12][6] - 1)*sc[21] + (TB[12][7] - 1)*sc[12] + (TB[12][8] - 1)*sc[22] + (TB[12][9] - 1)*sc[13] + (TB[12][10] - 1)*sc[23] + (TB[12][11] - 1)*sc[24] + (TB[12][12] - 1)*sc[25] + (TB[12][13] - 1)*sc[26] + (TB[12][14] - 1)*sc[27] + (TB[12][15] - 1)*sc[28] + (TB[12][16] - 1)*sc[17];
        alpha[13] = mixture + (TB[13][0] - 1)*sc[2] + (TB[13][1] - 1)*sc[7] + (TB[13][2] - 1)*sc[18] + (TB[13][3] - 1)*sc[9] + (TB[13][4] - 1)*sc[19] + (TB[13][5] - 1)*sc[20] + (TB[13][6] - 1)*sc[21] + (TB[13][7] - 1)*sc[12] + (TB[13][8] - 1)*sc[22] + (TB[13][9] - 1)*sc[13] + (TB[13][10] - 1)*sc[23] + (TB[13][11] - 1)*sc[24] + (TB[13][12] - 1)*sc[25] + (TB[13][13] - 1)*sc[26] + (TB[13][14] - 1)*sc[27] + (TB[13][15] - 1)*sc[28] + (TB[13][16] - 1)*sc[17];
        alpha[14] = mixture + (TB[14][0] - 1)*sc[2] + (TB[14][1] - 1)*sc[7] + (TB[14][2] - 1)*sc[18] + (TB[14][3] - 1)*sc[9] + (TB[14][4] - 1)*sc[19] + (TB[14][5] - 1)*sc[20] + (TB[14][6] - 1)*sc[21] + (TB[14][7] - 1)*sc[12] + (TB[14][8] - 1)*sc[22] + (TB[14][9] - 1)*sc[13] + (TB[14][10] - 1)*sc[23] + (TB[14][11] - 1)*sc[24] + (TB[14][12] - 1)*sc[25] + (TB[14][13] - 1)*sc[26] + (TB[14][14] - 1)*sc[27] + (TB[14][15] - 1)*sc[28] + (TB[14][16] - 1)*sc[17];
        for (int i=0; i<15; i++)
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
        alpha = mixture + (TB[15][0] - 1)*sc[2] + (TB[15][1] - 1)*sc[7] + (TB[15][2] - 1)*sc[18] + (TB[15][3] - 1)*sc[19] + (TB[15][4] - 1)*sc[20] + (TB[15][5] - 1)*sc[21] + (TB[15][6] - 1)*sc[12] + (TB[15][7] - 1)*sc[22] + (TB[15][8] - 1)*sc[13] + (TB[15][9] - 1)*sc[23] + (TB[15][10] - 1)*sc[24] + (TB[15][11] - 1)*sc[25] + (TB[15][12] - 1)*sc[26] + (TB[15][13] - 1)*sc[27] + (TB[15][14] - 1)*sc[28] + (TB[15][15] - 1)*sc[17];
        Corr[15] = alpha;
        alpha = mixture + (TB[16][0] - 1)*sc[2] + (TB[16][1] - 1)*sc[7] + (TB[16][2] - 1)*sc[18] + (TB[16][3] - 1)*sc[9] + (TB[16][4] - 1)*sc[19] + (TB[16][5] - 1)*sc[20] + (TB[16][6] - 1)*sc[21] + (TB[16][7] - 1)*sc[12] + (TB[16][8] - 1)*sc[22] + (TB[16][9] - 1)*sc[13] + (TB[16][10] - 1)*sc[23] + (TB[16][11] - 1)*sc[24] + (TB[16][12] - 1)*sc[25] + (TB[16][13] - 1)*sc[26] + (TB[16][14] - 1)*sc[27] + (TB[16][15] - 1)*sc[28] + (TB[16][16] - 1)*sc[17];
        Corr[16] = alpha;
        alpha = mixture + (TB[17][0] - 1)*sc[2] + (TB[17][1] - 1)*sc[7] + (TB[17][2] - 1)*sc[18] + (TB[17][3] - 1)*sc[9] + (TB[17][4] - 1)*sc[19] + (TB[17][5] - 1)*sc[20] + (TB[17][6] - 1)*sc[21] + (TB[17][7] - 1)*sc[12] + (TB[17][8] - 1)*sc[22] + (TB[17][9] - 1)*sc[13] + (TB[17][10] - 1)*sc[23] + (TB[17][11] - 1)*sc[24] + (TB[17][12] - 1)*sc[25] + (TB[17][13] - 1)*sc[26] + (TB[17][14] - 1)*sc[27] + (TB[17][15] - 1)*sc[28] + (TB[17][16] - 1)*sc[17];
        Corr[17] = alpha;
        alpha = mixture + (TB[18][0] - 1)*sc[2] + (TB[18][1] - 1)*sc[7] + (TB[18][2] - 1)*sc[18] + (TB[18][3] - 1)*sc[9] + (TB[18][4] - 1)*sc[19] + (TB[18][5] - 1)*sc[20] + (TB[18][6] - 1)*sc[21] + (TB[18][7] - 1)*sc[12] + (TB[18][8] - 1)*sc[22] + (TB[18][9] - 1)*sc[13] + (TB[18][10] - 1)*sc[23] + (TB[18][11] - 1)*sc[24] + (TB[18][12] - 1)*sc[25] + (TB[18][13] - 1)*sc[26] + (TB[18][14] - 1)*sc[27] + (TB[18][15] - 1)*sc[28] + (TB[18][16] - 1)*sc[17];
        Corr[18] = alpha;
        alpha = mixture + (TB[19][0] - 1)*sc[2] + (TB[19][1] - 1)*sc[7] + (TB[19][2] - 1)*sc[18] + (TB[19][3] - 1)*sc[9] + (TB[19][4] - 1)*sc[19] + (TB[19][5] - 1)*sc[20] + (TB[19][6] - 1)*sc[21] + (TB[19][7] - 1)*sc[12] + (TB[19][8] - 1)*sc[22] + (TB[19][9] - 1)*sc[13] + (TB[19][10] - 1)*sc[23] + (TB[19][11] - 1)*sc[24] + (TB[19][12] - 1)*sc[25] + (TB[19][13] - 1)*sc[26] + (TB[19][14] - 1)*sc[27] + (TB[19][15] - 1)*sc[28] + (TB[19][16] - 1)*sc[17];
        Corr[19] = alpha;
    }

    for (int i=0; i<172; i++)
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

    amrex::Real q_f[172], q_r[172];
    amrex::Real sc_qss[1];
    comp_qfqr(q_f, q_r, sc, sc_qss, tc, invT);

    for (int i = 0; i < 172; ++i) {
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
    amrex::Real c[29]; /*temporary storage */
    amrex::Real PORT = 1e6 * (*P)/(8.31446261815324e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 29; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRateFR(q_f, q_r, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 172; ++id) {
        q_f[id] *= 1.0e-6;
        q_r[id] *= 1.0e-6;
    }
}


/*Returns the rate of progress for each reaction */
void CKQC(amrex::Real *  T, amrex::Real *  C, amrex::Real *  qdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 29; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    progressRate(qdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 29; ++id) {
        C[id] *= 1.0e-6;
    }

    for (id = 0; id < 172; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mass fractions */
void CKQYP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y, amrex::Real *  qdot)
{
    int id; /*loop counter */
    amrex::Real c[29]; /*temporary storage */
    amrex::Real YOW = 0; 
    amrex::Real PWORT; 
    amrex::Real imw[29];

    get_imw(imw);

    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*N2 */
    YOW += y[1]*imw[1]; /*H */
    YOW += y[2]*imw[2]; /*H2 */
    YOW += y[3]*imw[3]; /*O */
    YOW += y[4]*imw[4]; /*OH */
    YOW += y[5]*imw[5]; /*O2 */
    YOW += y[6]*imw[6]; /*H2O2 */
    YOW += y[7]*imw[7]; /*H2O */
    YOW += y[8]*imw[8]; /*HO2 */
    YOW += y[9]*imw[9]; /*CO */
    YOW += y[10]*imw[10]; /*CH3 */
    YOW += y[11]*imw[11]; /*CH2O */
    YOW += y[12]*imw[12]; /*CO2 */
    YOW += y[13]*imw[13]; /*CH4 */
    YOW += y[14]*imw[14]; /*C2H2 */
    YOW += y[15]*imw[15]; /*C2H4 */
    YOW += y[16]*imw[16]; /*CH2CO */
    YOW += y[17]*imw[17]; /*C2H6 */
    YOW += y[18]*imw[18]; /*C */
    YOW += y[19]*imw[19]; /*CH */
    YOW += y[20]*imw[20]; /*HCO */
    YOW += y[21]*imw[21]; /*TXCH2 */
    YOW += y[22]*imw[22]; /*SXCH2 */
    YOW += y[23]*imw[23]; /*C2H3 */
    YOW += y[24]*imw[24]; /*C2H5 */
    YOW += y[25]*imw[25]; /*HCCO */
    YOW += y[26]*imw[26]; /*CH3CHO */
    YOW += y[27]*imw[27]; /*CH2CHO */
    YOW += y[28]*imw[28]; /*C2H5O */
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

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 172; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKQXP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x, amrex::Real *  qdot)
{
    int id; /*loop counter */
    amrex::Real c[29]; /*temporary storage */
    amrex::Real PORT = 1e6 * (*P)/(8.31446261815324e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 29; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 172; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mass fractions */
void CKQYR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y, amrex::Real *  qdot)
{
    int id; /*loop counter */
    amrex::Real c[29]; /*temporary storage */
    amrex::Real imw[29];

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

    /*call progressRate */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 172; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mole fractions */
void CKQXR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  x, amrex::Real *  qdot)
{
    int id; /*loop counter */
    amrex::Real c[29]; /*temporary storage */
    amrex::Real XW = 0; /*See Eq 4, 11 in CK Manual */
    amrex::Real ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*28.013400; /*N2 */
    XW += x[1]*1.007970; /*H */
    XW += x[2]*2.015940; /*H2 */
    XW += x[3]*15.999400; /*O */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*31.998800; /*O2 */
    XW += x[6]*34.014740; /*H2O2 */
    XW += x[7]*18.015340; /*H2O */
    XW += x[8]*33.006770; /*HO2 */
    XW += x[9]*28.010550; /*CO */
    XW += x[10]*15.035060; /*CH3 */
    XW += x[11]*30.026490; /*CH2O */
    XW += x[12]*44.009950; /*CO2 */
    XW += x[13]*16.043030; /*CH4 */
    XW += x[14]*26.038240; /*C2H2 */
    XW += x[15]*28.054180; /*C2H4 */
    XW += x[16]*42.037640; /*CH2CO */
    XW += x[17]*30.070120; /*C2H6 */
    XW += x[18]*12.011150; /*C */
    XW += x[19]*13.019120; /*CH */
    XW += x[20]*29.018520; /*HCO */
    XW += x[21]*14.027090; /*TXCH2 */
    XW += x[22]*14.027090; /*SXCH2 */
    XW += x[23]*27.046210; /*C2H3 */
    XW += x[24]*29.062150; /*C2H5 */
    XW += x[25]*41.029670; /*HCCO */
    XW += x[26]*44.053580; /*CH3CHO */
    XW += x[27]*43.045610; /*CH2CHO */
    XW += x[28]*45.061550; /*C2H5O */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 29; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 172; ++id) {
        qdot[id] *= 1.0e-6;
    }
}

/*compute the reaction Jacobian on CPU */
void aJacobian(amrex::Real *  J, amrex::Real *  sc, amrex::Real T, int consP)
{
    for (int i=0; i<900; i++) {
        J[i] = 0.0;
    }

    amrex::Real wdot[29];
    for (int k=0; k<29; k++) {
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
    for (int k = 0; k < 29; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    amrex::Real g_RT[29];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    amrex::Real h_RT[29];
    speciesEnthalpy(h_RT, tc);

    amrex::Real phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    amrex::Real dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    amrex::Real dqdci, dcdc_fac, dqdc[29];
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
    alpha = mixture + (TB[0][0] - 1)*sc[2] + (TB[0][1] - 1)*sc[7] + (TB[0][2] - 1)*sc[18] + (TB[0][3] - 1)*sc[9] + (TB[0][4] - 1)*sc[19] + (TB[0][5] - 1)*sc[20] + (TB[0][6] - 1)*sc[21] + (TB[0][7] - 1)*sc[12] + (TB[0][8] - 1)*sc[22] + (TB[0][9] - 1)*sc[13] + (TB[0][10] - 1)*sc[23] + (TB[0][11] - 1)*sc[24] + (TB[0][12] - 1)*sc[25] + (TB[0][13] - 1)*sc[26] + (TB[0][14] - 1)*sc[27] + (TB[0][15] - 1)*sc[28] + (TB[0][16] - 1)*sc[17];
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
    phi_r = sc[6];
    Kc = refCinv * exp(2.000000*g_RT[4] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[4]) + (h_RT[6]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[4] -= 2 * q; /* OH */
    wdot[6] += q; /* H2O2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[0][0] - 1)*dcdc_fac;
        J[64] += -2 * dqdci;          /* dwdot[OH]/d[H2] */
        J[66] += dqdci;               /* dwdot[H2O2]/d[H2] */
        /* d()/d[OH] */
        dqdci =  + k_f*2.000000*sc[4];
        J[124] += -2 * dqdci;         /* dwdot[OH]/d[OH] */
        J[126] += dqdci;              /* dwdot[H2O2]/d[OH] */
        /* d()/d[H2O2] */
        dqdci =  - k_r;
        J[184] += -2 * dqdci;         /* dwdot[OH]/d[H2O2] */
        J[186] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
        /* d()/d[H2O] */
        dqdci = (TB[0][1] - 1)*dcdc_fac;
        J[214] += -2 * dqdci;         /* dwdot[OH]/d[H2O] */
        J[216] += dqdci;              /* dwdot[H2O2]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[0][3] - 1)*dcdc_fac;
        J[274] += -2 * dqdci;         /* dwdot[OH]/d[CO] */
        J[276] += dqdci;              /* dwdot[H2O2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[0][7] - 1)*dcdc_fac;
        J[364] += -2 * dqdci;         /* dwdot[OH]/d[CO2] */
        J[366] += dqdci;              /* dwdot[H2O2]/d[CO2] */
        /* d()/d[CH4] */
        dqdci = (TB[0][9] - 1)*dcdc_fac;
        J[394] += -2 * dqdci;         /* dwdot[OH]/d[CH4] */
        J[396] += dqdci;              /* dwdot[H2O2]/d[CH4] */
        /* d()/d[C2H6] */
        dqdci = (TB[0][16] - 1)*dcdc_fac;
        J[514] += -2 * dqdci;         /* dwdot[OH]/d[C2H6] */
        J[516] += dqdci;              /* dwdot[H2O2]/d[C2H6] */
        /* d()/d[C] */
        dqdci = (TB[0][2] - 1)*dcdc_fac;
        J[544] += -2 * dqdci;         /* dwdot[OH]/d[C] */
        J[546] += dqdci;              /* dwdot[H2O2]/d[C] */
        /* d()/d[CH] */
        dqdci = (TB[0][4] - 1)*dcdc_fac;
        J[574] += -2 * dqdci;         /* dwdot[OH]/d[CH] */
        J[576] += dqdci;              /* dwdot[H2O2]/d[CH] */
        /* d()/d[HCO] */
        dqdci = (TB[0][5] - 1)*dcdc_fac;
        J[604] += -2 * dqdci;         /* dwdot[OH]/d[HCO] */
        J[606] += dqdci;              /* dwdot[H2O2]/d[HCO] */
        /* d()/d[TXCH2] */
        dqdci = (TB[0][6] - 1)*dcdc_fac;
        J[634] += -2 * dqdci;         /* dwdot[OH]/d[TXCH2] */
        J[636] += dqdci;              /* dwdot[H2O2]/d[TXCH2] */
        /* d()/d[SXCH2] */
        dqdci = (TB[0][8] - 1)*dcdc_fac;
        J[664] += -2 * dqdci;         /* dwdot[OH]/d[SXCH2] */
        J[666] += dqdci;              /* dwdot[H2O2]/d[SXCH2] */
        /* d()/d[C2H3] */
        dqdci = (TB[0][10] - 1)*dcdc_fac;
        J[694] += -2 * dqdci;         /* dwdot[OH]/d[C2H3] */
        J[696] += dqdci;              /* dwdot[H2O2]/d[C2H3] */
        /* d()/d[C2H5] */
        dqdci = (TB[0][11] - 1)*dcdc_fac;
        J[724] += -2 * dqdci;         /* dwdot[OH]/d[C2H5] */
        J[726] += dqdci;              /* dwdot[H2O2]/d[C2H5] */
        /* d()/d[HCCO] */
        dqdci = (TB[0][12] - 1)*dcdc_fac;
        J[754] += -2 * dqdci;         /* dwdot[OH]/d[HCCO] */
        J[756] += dqdci;              /* dwdot[H2O2]/d[HCCO] */
        /* d()/d[CH3CHO] */
        dqdci = (TB[0][13] - 1)*dcdc_fac;
        J[784] += -2 * dqdci;         /* dwdot[OH]/d[CH3CHO] */
        J[786] += dqdci;              /* dwdot[H2O2]/d[CH3CHO] */
        /* d()/d[CH2CHO] */
        dqdci = (TB[0][14] - 1)*dcdc_fac;
        J[814] += -2 * dqdci;         /* dwdot[OH]/d[CH2CHO] */
        J[816] += dqdci;              /* dwdot[H2O2]/d[CH2CHO] */
        /* d()/d[C2H5O] */
        dqdci = (TB[0][15] - 1)*dcdc_fac;
        J[844] += -2 * dqdci;         /* dwdot[OH]/d[C2H5O] */
        J[846] += dqdci;              /* dwdot[H2O2]/d[C2H5O] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = TB[0][0]*dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac + k_f*2.000000*sc[4];
        dqdc[5] = dcdc_fac;
        dqdc[6] = dcdc_fac - k_r;
        dqdc[7] = TB[0][1]*dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = TB[0][3]*dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = TB[0][7]*dcdc_fac;
        dqdc[13] = TB[0][9]*dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = TB[0][16]*dcdc_fac;
        dqdc[18] = TB[0][2]*dcdc_fac;
        dqdc[19] = TB[0][4]*dcdc_fac;
        dqdc[20] = TB[0][5]*dcdc_fac;
        dqdc[21] = TB[0][6]*dcdc_fac;
        dqdc[22] = TB[0][8]*dcdc_fac;
        dqdc[23] = TB[0][10]*dcdc_fac;
        dqdc[24] = TB[0][11]*dcdc_fac;
        dqdc[25] = TB[0][12]*dcdc_fac;
        dqdc[26] = TB[0][13]*dcdc_fac;
        dqdc[27] = TB[0][14]*dcdc_fac;
        dqdc[28] = TB[0][15]*dcdc_fac;
        for (int k=0; k<29; k++) {
            J[30*k+4] += -2 * dqdc[k];
            J[30*k+6] += dqdc[k];
        }
    }
    J[874] += -2 * dqdT; /* dwdot[OH]/dT */
    J[876] += dqdT; /* dwdot[H2O2]/dT */

    /*reaction 2: H + O2 (+M) <=> HO2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[1][0] - 1)*sc[2] + (TB[1][1] - 1)*sc[5] + (TB[1][2] - 1)*sc[7] + (TB[1][3] - 1)*sc[18] + (TB[1][4] - 1)*sc[9] + (TB[1][5] - 1)*sc[19] + (TB[1][6] - 1)*sc[20] + (TB[1][7] - 1)*sc[21] + (TB[1][8] - 1)*sc[12] + (TB[1][9] - 1)*sc[22] + (TB[1][10] - 1)*sc[23] + (TB[1][11] - 1)*sc[24] + (TB[1][12] - 1)*sc[25] + (TB[1][13] - 1)*sc[26] + (TB[1][14] - 1)*sc[27] + (TB[1][15] - 1)*sc[28];
    /* forward */
    phi_f = sc[1]*sc[5];
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
    phi_r = sc[8];
    Kc = refCinv * exp(g_RT[1] + g_RT[5] - g_RT[8]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[5]) + (h_RT[8]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[5] -= q; /* O2 */
    wdot[8] += q; /* HO2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[5];
        J[31] -= dqdci;               /* dwdot[H]/d[H] */
        J[35] -= dqdci;               /* dwdot[O2]/d[H] */
        J[38] += dqdci;               /* dwdot[HO2]/d[H] */
        /* d()/d[H2] */
        dqdci = (TB[1][0] - 1)*dcdc_fac;
        J[61] -= dqdci;               /* dwdot[H]/d[H2] */
        J[65] -= dqdci;               /* dwdot[O2]/d[H2] */
        J[68] += dqdci;               /* dwdot[HO2]/d[H2] */
        /* d()/d[O2] */
        dqdci = (TB[1][1] - 1)*dcdc_fac + k_f*sc[1];
        J[151] -= dqdci;              /* dwdot[H]/d[O2] */
        J[155] -= dqdci;              /* dwdot[O2]/d[O2] */
        J[158] += dqdci;              /* dwdot[HO2]/d[O2] */
        /* d()/d[H2O] */
        dqdci = (TB[1][2] - 1)*dcdc_fac;
        J[211] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[215] -= dqdci;              /* dwdot[O2]/d[H2O] */
        J[218] += dqdci;              /* dwdot[HO2]/d[H2O] */
        /* d()/d[HO2] */
        dqdci =  - k_r;
        J[241] -= dqdci;              /* dwdot[H]/d[HO2] */
        J[245] -= dqdci;              /* dwdot[O2]/d[HO2] */
        J[248] += dqdci;              /* dwdot[HO2]/d[HO2] */
        /* d()/d[CO] */
        dqdci = (TB[1][4] - 1)*dcdc_fac;
        J[271] -= dqdci;              /* dwdot[H]/d[CO] */
        J[275] -= dqdci;              /* dwdot[O2]/d[CO] */
        J[278] += dqdci;              /* dwdot[HO2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[1][8] - 1)*dcdc_fac;
        J[361] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[365] -= dqdci;              /* dwdot[O2]/d[CO2] */
        J[368] += dqdci;              /* dwdot[HO2]/d[CO2] */
        /* d()/d[C] */
        dqdci = (TB[1][3] - 1)*dcdc_fac;
        J[541] -= dqdci;              /* dwdot[H]/d[C] */
        J[545] -= dqdci;              /* dwdot[O2]/d[C] */
        J[548] += dqdci;              /* dwdot[HO2]/d[C] */
        /* d()/d[CH] */
        dqdci = (TB[1][5] - 1)*dcdc_fac;
        J[571] -= dqdci;              /* dwdot[H]/d[CH] */
        J[575] -= dqdci;              /* dwdot[O2]/d[CH] */
        J[578] += dqdci;              /* dwdot[HO2]/d[CH] */
        /* d()/d[HCO] */
        dqdci = (TB[1][6] - 1)*dcdc_fac;
        J[601] -= dqdci;              /* dwdot[H]/d[HCO] */
        J[605] -= dqdci;              /* dwdot[O2]/d[HCO] */
        J[608] += dqdci;              /* dwdot[HO2]/d[HCO] */
        /* d()/d[TXCH2] */
        dqdci = (TB[1][7] - 1)*dcdc_fac;
        J[631] -= dqdci;              /* dwdot[H]/d[TXCH2] */
        J[635] -= dqdci;              /* dwdot[O2]/d[TXCH2] */
        J[638] += dqdci;              /* dwdot[HO2]/d[TXCH2] */
        /* d()/d[SXCH2] */
        dqdci = (TB[1][9] - 1)*dcdc_fac;
        J[661] -= dqdci;              /* dwdot[H]/d[SXCH2] */
        J[665] -= dqdci;              /* dwdot[O2]/d[SXCH2] */
        J[668] += dqdci;              /* dwdot[HO2]/d[SXCH2] */
        /* d()/d[C2H3] */
        dqdci = (TB[1][10] - 1)*dcdc_fac;
        J[691] -= dqdci;              /* dwdot[H]/d[C2H3] */
        J[695] -= dqdci;              /* dwdot[O2]/d[C2H3] */
        J[698] += dqdci;              /* dwdot[HO2]/d[C2H3] */
        /* d()/d[C2H5] */
        dqdci = (TB[1][11] - 1)*dcdc_fac;
        J[721] -= dqdci;              /* dwdot[H]/d[C2H5] */
        J[725] -= dqdci;              /* dwdot[O2]/d[C2H5] */
        J[728] += dqdci;              /* dwdot[HO2]/d[C2H5] */
        /* d()/d[HCCO] */
        dqdci = (TB[1][12] - 1)*dcdc_fac;
        J[751] -= dqdci;              /* dwdot[H]/d[HCCO] */
        J[755] -= dqdci;              /* dwdot[O2]/d[HCCO] */
        J[758] += dqdci;              /* dwdot[HO2]/d[HCCO] */
        /* d()/d[CH3CHO] */
        dqdci = (TB[1][13] - 1)*dcdc_fac;
        J[781] -= dqdci;              /* dwdot[H]/d[CH3CHO] */
        J[785] -= dqdci;              /* dwdot[O2]/d[CH3CHO] */
        J[788] += dqdci;              /* dwdot[HO2]/d[CH3CHO] */
        /* d()/d[CH2CHO] */
        dqdci = (TB[1][14] - 1)*dcdc_fac;
        J[811] -= dqdci;              /* dwdot[H]/d[CH2CHO] */
        J[815] -= dqdci;              /* dwdot[O2]/d[CH2CHO] */
        J[818] += dqdci;              /* dwdot[HO2]/d[CH2CHO] */
        /* d()/d[C2H5O] */
        dqdci = (TB[1][15] - 1)*dcdc_fac;
        J[841] -= dqdci;              /* dwdot[H]/d[C2H5O] */
        J[845] -= dqdci;              /* dwdot[O2]/d[C2H5O] */
        J[848] += dqdci;              /* dwdot[HO2]/d[C2H5O] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[5];
        dqdc[2] = TB[1][0]*dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = TB[1][1]*dcdc_fac + k_f*sc[1];
        dqdc[6] = dcdc_fac;
        dqdc[7] = TB[1][2]*dcdc_fac;
        dqdc[8] = dcdc_fac - k_r;
        dqdc[9] = TB[1][4]*dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = TB[1][8]*dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = TB[1][3]*dcdc_fac;
        dqdc[19] = TB[1][5]*dcdc_fac;
        dqdc[20] = TB[1][6]*dcdc_fac;
        dqdc[21] = TB[1][7]*dcdc_fac;
        dqdc[22] = TB[1][9]*dcdc_fac;
        dqdc[23] = TB[1][10]*dcdc_fac;
        dqdc[24] = TB[1][11]*dcdc_fac;
        dqdc[25] = TB[1][12]*dcdc_fac;
        dqdc[26] = TB[1][13]*dcdc_fac;
        dqdc[27] = TB[1][14]*dcdc_fac;
        dqdc[28] = TB[1][15]*dcdc_fac;
        for (int k=0; k<29; k++) {
            J[30*k+1] -= dqdc[k];
            J[30*k+5] -= dqdc[k];
            J[30*k+8] += dqdc[k];
        }
    }
    J[871] -= dqdT; /* dwdot[H]/dT */
    J[875] -= dqdT; /* dwdot[O2]/dT */
    J[878] += dqdT; /* dwdot[HO2]/dT */

    /*reaction 3: CH + H2 (+M) <=> CH3 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[2][0] - 1)*sc[2] + (TB[2][1] - 1)*sc[7] + (TB[2][2] - 1)*sc[18] + (TB[2][3] - 1)*sc[9] + (TB[2][4] - 1)*sc[19] + (TB[2][5] - 1)*sc[20] + (TB[2][6] - 1)*sc[21] + (TB[2][7] - 1)*sc[12] + (TB[2][8] - 1)*sc[22] + (TB[2][9] - 1)*sc[13] + (TB[2][10] - 1)*sc[23] + (TB[2][11] - 1)*sc[24] + (TB[2][12] - 1)*sc[25] + (TB[2][13] - 1)*sc[26] + (TB[2][14] - 1)*sc[27] + (TB[2][15] - 1)*sc[28] + (TB[2][16] - 1)*sc[17];
    /* forward */
    phi_f = sc[2]*sc[19];
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
    phi_r = sc[10];
    Kc = refCinv * exp(g_RT[2] - g_RT[10] + g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[19]) + (h_RT[10]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[2] -= q; /* H2 */
    wdot[10] += q; /* CH3 */
    wdot[19] -= q; /* CH */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[2][0] - 1)*dcdc_fac + k_f*sc[19];
        J[62] -= dqdci;               /* dwdot[H2]/d[H2] */
        J[70] += dqdci;               /* dwdot[CH3]/d[H2] */
        J[79] -= dqdci;               /* dwdot[CH]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[2][1] - 1)*dcdc_fac;
        J[212] -= dqdci;              /* dwdot[H2]/d[H2O] */
        J[220] += dqdci;              /* dwdot[CH3]/d[H2O] */
        J[229] -= dqdci;              /* dwdot[CH]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[2][3] - 1)*dcdc_fac;
        J[272] -= dqdci;              /* dwdot[H2]/d[CO] */
        J[280] += dqdci;              /* dwdot[CH3]/d[CO] */
        J[289] -= dqdci;              /* dwdot[CH]/d[CO] */
        /* d()/d[CH3] */
        dqdci =  - k_r;
        J[302] -= dqdci;              /* dwdot[H2]/d[CH3] */
        J[310] += dqdci;              /* dwdot[CH3]/d[CH3] */
        J[319] -= dqdci;              /* dwdot[CH]/d[CH3] */
        /* d()/d[CO2] */
        dqdci = (TB[2][7] - 1)*dcdc_fac;
        J[362] -= dqdci;              /* dwdot[H2]/d[CO2] */
        J[370] += dqdci;              /* dwdot[CH3]/d[CO2] */
        J[379] -= dqdci;              /* dwdot[CH]/d[CO2] */
        /* d()/d[CH4] */
        dqdci = (TB[2][9] - 1)*dcdc_fac;
        J[392] -= dqdci;              /* dwdot[H2]/d[CH4] */
        J[400] += dqdci;              /* dwdot[CH3]/d[CH4] */
        J[409] -= dqdci;              /* dwdot[CH]/d[CH4] */
        /* d()/d[C2H6] */
        dqdci = (TB[2][16] - 1)*dcdc_fac;
        J[512] -= dqdci;              /* dwdot[H2]/d[C2H6] */
        J[520] += dqdci;              /* dwdot[CH3]/d[C2H6] */
        J[529] -= dqdci;              /* dwdot[CH]/d[C2H6] */
        /* d()/d[C] */
        dqdci = (TB[2][2] - 1)*dcdc_fac;
        J[542] -= dqdci;              /* dwdot[H2]/d[C] */
        J[550] += dqdci;              /* dwdot[CH3]/d[C] */
        J[559] -= dqdci;              /* dwdot[CH]/d[C] */
        /* d()/d[CH] */
        dqdci = (TB[2][4] - 1)*dcdc_fac + k_f*sc[2];
        J[572] -= dqdci;              /* dwdot[H2]/d[CH] */
        J[580] += dqdci;              /* dwdot[CH3]/d[CH] */
        J[589] -= dqdci;              /* dwdot[CH]/d[CH] */
        /* d()/d[HCO] */
        dqdci = (TB[2][5] - 1)*dcdc_fac;
        J[602] -= dqdci;              /* dwdot[H2]/d[HCO] */
        J[610] += dqdci;              /* dwdot[CH3]/d[HCO] */
        J[619] -= dqdci;              /* dwdot[CH]/d[HCO] */
        /* d()/d[TXCH2] */
        dqdci = (TB[2][6] - 1)*dcdc_fac;
        J[632] -= dqdci;              /* dwdot[H2]/d[TXCH2] */
        J[640] += dqdci;              /* dwdot[CH3]/d[TXCH2] */
        J[649] -= dqdci;              /* dwdot[CH]/d[TXCH2] */
        /* d()/d[SXCH2] */
        dqdci = (TB[2][8] - 1)*dcdc_fac;
        J[662] -= dqdci;              /* dwdot[H2]/d[SXCH2] */
        J[670] += dqdci;              /* dwdot[CH3]/d[SXCH2] */
        J[679] -= dqdci;              /* dwdot[CH]/d[SXCH2] */
        /* d()/d[C2H3] */
        dqdci = (TB[2][10] - 1)*dcdc_fac;
        J[692] -= dqdci;              /* dwdot[H2]/d[C2H3] */
        J[700] += dqdci;              /* dwdot[CH3]/d[C2H3] */
        J[709] -= dqdci;              /* dwdot[CH]/d[C2H3] */
        /* d()/d[C2H5] */
        dqdci = (TB[2][11] - 1)*dcdc_fac;
        J[722] -= dqdci;              /* dwdot[H2]/d[C2H5] */
        J[730] += dqdci;              /* dwdot[CH3]/d[C2H5] */
        J[739] -= dqdci;              /* dwdot[CH]/d[C2H5] */
        /* d()/d[HCCO] */
        dqdci = (TB[2][12] - 1)*dcdc_fac;
        J[752] -= dqdci;              /* dwdot[H2]/d[HCCO] */
        J[760] += dqdci;              /* dwdot[CH3]/d[HCCO] */
        J[769] -= dqdci;              /* dwdot[CH]/d[HCCO] */
        /* d()/d[CH3CHO] */
        dqdci = (TB[2][13] - 1)*dcdc_fac;
        J[782] -= dqdci;              /* dwdot[H2]/d[CH3CHO] */
        J[790] += dqdci;              /* dwdot[CH3]/d[CH3CHO] */
        J[799] -= dqdci;              /* dwdot[CH]/d[CH3CHO] */
        /* d()/d[CH2CHO] */
        dqdci = (TB[2][14] - 1)*dcdc_fac;
        J[812] -= dqdci;              /* dwdot[H2]/d[CH2CHO] */
        J[820] += dqdci;              /* dwdot[CH3]/d[CH2CHO] */
        J[829] -= dqdci;              /* dwdot[CH]/d[CH2CHO] */
        /* d()/d[C2H5O] */
        dqdci = (TB[2][15] - 1)*dcdc_fac;
        J[842] -= dqdci;              /* dwdot[H2]/d[C2H5O] */
        J[850] += dqdci;              /* dwdot[CH3]/d[C2H5O] */
        J[859] -= dqdci;              /* dwdot[CH]/d[C2H5O] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = TB[2][0]*dcdc_fac + k_f*sc[19];
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = TB[2][1]*dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = TB[2][3]*dcdc_fac;
        dqdc[10] = dcdc_fac - k_r;
        dqdc[11] = dcdc_fac;
        dqdc[12] = TB[2][7]*dcdc_fac;
        dqdc[13] = TB[2][9]*dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = TB[2][16]*dcdc_fac;
        dqdc[18] = TB[2][2]*dcdc_fac;
        dqdc[19] = TB[2][4]*dcdc_fac + k_f*sc[2];
        dqdc[20] = TB[2][5]*dcdc_fac;
        dqdc[21] = TB[2][6]*dcdc_fac;
        dqdc[22] = TB[2][8]*dcdc_fac;
        dqdc[23] = TB[2][10]*dcdc_fac;
        dqdc[24] = TB[2][11]*dcdc_fac;
        dqdc[25] = TB[2][12]*dcdc_fac;
        dqdc[26] = TB[2][13]*dcdc_fac;
        dqdc[27] = TB[2][14]*dcdc_fac;
        dqdc[28] = TB[2][15]*dcdc_fac;
        for (int k=0; k<29; k++) {
            J[30*k+2] -= dqdc[k];
            J[30*k+10] += dqdc[k];
            J[30*k+19] -= dqdc[k];
        }
    }
    J[872] -= dqdT; /* dwdot[H2]/dT */
    J[880] += dqdT; /* dwdot[CH3]/dT */
    J[889] -= dqdT; /* dwdot[CH]/dT */

    /*reaction 4: TXCH2 + H (+M) <=> CH3 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[3][0] - 1)*sc[2] + (TB[3][1] - 1)*sc[7] + (TB[3][2] - 1)*sc[18] + (TB[3][3] - 1)*sc[9] + (TB[3][4] - 1)*sc[19] + (TB[3][5] - 1)*sc[20] + (TB[3][6] - 1)*sc[21] + (TB[3][7] - 1)*sc[12] + (TB[3][8] - 1)*sc[22] + (TB[3][9] - 1)*sc[13] + (TB[3][10] - 1)*sc[23] + (TB[3][11] - 1)*sc[24] + (TB[3][12] - 1)*sc[25] + (TB[3][13] - 1)*sc[26] + (TB[3][14] - 1)*sc[27] + (TB[3][15] - 1)*sc[28] + (TB[3][16] - 1)*sc[17];
    /* forward */
    phi_f = sc[1]*sc[21];
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
    phi_r = sc[10];
    Kc = refCinv * exp(g_RT[1] - g_RT[10] + g_RT[21]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[21]) + (h_RT[10]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[10] += q; /* CH3 */
    wdot[21] -= q; /* TXCH2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[21];
        J[31] -= dqdci;               /* dwdot[H]/d[H] */
        J[40] += dqdci;               /* dwdot[CH3]/d[H] */
        J[51] -= dqdci;               /* dwdot[TXCH2]/d[H] */
        /* d()/d[H2] */
        dqdci = (TB[3][0] - 1)*dcdc_fac;
        J[61] -= dqdci;               /* dwdot[H]/d[H2] */
        J[70] += dqdci;               /* dwdot[CH3]/d[H2] */
        J[81] -= dqdci;               /* dwdot[TXCH2]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[3][1] - 1)*dcdc_fac;
        J[211] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[220] += dqdci;              /* dwdot[CH3]/d[H2O] */
        J[231] -= dqdci;              /* dwdot[TXCH2]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[3][3] - 1)*dcdc_fac;
        J[271] -= dqdci;              /* dwdot[H]/d[CO] */
        J[280] += dqdci;              /* dwdot[CH3]/d[CO] */
        J[291] -= dqdci;              /* dwdot[TXCH2]/d[CO] */
        /* d()/d[CH3] */
        dqdci =  - k_r;
        J[301] -= dqdci;              /* dwdot[H]/d[CH3] */
        J[310] += dqdci;              /* dwdot[CH3]/d[CH3] */
        J[321] -= dqdci;              /* dwdot[TXCH2]/d[CH3] */
        /* d()/d[CO2] */
        dqdci = (TB[3][7] - 1)*dcdc_fac;
        J[361] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[370] += dqdci;              /* dwdot[CH3]/d[CO2] */
        J[381] -= dqdci;              /* dwdot[TXCH2]/d[CO2] */
        /* d()/d[CH4] */
        dqdci = (TB[3][9] - 1)*dcdc_fac;
        J[391] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[400] += dqdci;              /* dwdot[CH3]/d[CH4] */
        J[411] -= dqdci;              /* dwdot[TXCH2]/d[CH4] */
        /* d()/d[C2H6] */
        dqdci = (TB[3][16] - 1)*dcdc_fac;
        J[511] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[520] += dqdci;              /* dwdot[CH3]/d[C2H6] */
        J[531] -= dqdci;              /* dwdot[TXCH2]/d[C2H6] */
        /* d()/d[C] */
        dqdci = (TB[3][2] - 1)*dcdc_fac;
        J[541] -= dqdci;              /* dwdot[H]/d[C] */
        J[550] += dqdci;              /* dwdot[CH3]/d[C] */
        J[561] -= dqdci;              /* dwdot[TXCH2]/d[C] */
        /* d()/d[CH] */
        dqdci = (TB[3][4] - 1)*dcdc_fac;
        J[571] -= dqdci;              /* dwdot[H]/d[CH] */
        J[580] += dqdci;              /* dwdot[CH3]/d[CH] */
        J[591] -= dqdci;              /* dwdot[TXCH2]/d[CH] */
        /* d()/d[HCO] */
        dqdci = (TB[3][5] - 1)*dcdc_fac;
        J[601] -= dqdci;              /* dwdot[H]/d[HCO] */
        J[610] += dqdci;              /* dwdot[CH3]/d[HCO] */
        J[621] -= dqdci;              /* dwdot[TXCH2]/d[HCO] */
        /* d()/d[TXCH2] */
        dqdci = (TB[3][6] - 1)*dcdc_fac + k_f*sc[1];
        J[631] -= dqdci;              /* dwdot[H]/d[TXCH2] */
        J[640] += dqdci;              /* dwdot[CH3]/d[TXCH2] */
        J[651] -= dqdci;              /* dwdot[TXCH2]/d[TXCH2] */
        /* d()/d[SXCH2] */
        dqdci = (TB[3][8] - 1)*dcdc_fac;
        J[661] -= dqdci;              /* dwdot[H]/d[SXCH2] */
        J[670] += dqdci;              /* dwdot[CH3]/d[SXCH2] */
        J[681] -= dqdci;              /* dwdot[TXCH2]/d[SXCH2] */
        /* d()/d[C2H3] */
        dqdci = (TB[3][10] - 1)*dcdc_fac;
        J[691] -= dqdci;              /* dwdot[H]/d[C2H3] */
        J[700] += dqdci;              /* dwdot[CH3]/d[C2H3] */
        J[711] -= dqdci;              /* dwdot[TXCH2]/d[C2H3] */
        /* d()/d[C2H5] */
        dqdci = (TB[3][11] - 1)*dcdc_fac;
        J[721] -= dqdci;              /* dwdot[H]/d[C2H5] */
        J[730] += dqdci;              /* dwdot[CH3]/d[C2H5] */
        J[741] -= dqdci;              /* dwdot[TXCH2]/d[C2H5] */
        /* d()/d[HCCO] */
        dqdci = (TB[3][12] - 1)*dcdc_fac;
        J[751] -= dqdci;              /* dwdot[H]/d[HCCO] */
        J[760] += dqdci;              /* dwdot[CH3]/d[HCCO] */
        J[771] -= dqdci;              /* dwdot[TXCH2]/d[HCCO] */
        /* d()/d[CH3CHO] */
        dqdci = (TB[3][13] - 1)*dcdc_fac;
        J[781] -= dqdci;              /* dwdot[H]/d[CH3CHO] */
        J[790] += dqdci;              /* dwdot[CH3]/d[CH3CHO] */
        J[801] -= dqdci;              /* dwdot[TXCH2]/d[CH3CHO] */
        /* d()/d[CH2CHO] */
        dqdci = (TB[3][14] - 1)*dcdc_fac;
        J[811] -= dqdci;              /* dwdot[H]/d[CH2CHO] */
        J[820] += dqdci;              /* dwdot[CH3]/d[CH2CHO] */
        J[831] -= dqdci;              /* dwdot[TXCH2]/d[CH2CHO] */
        /* d()/d[C2H5O] */
        dqdci = (TB[3][15] - 1)*dcdc_fac;
        J[841] -= dqdci;              /* dwdot[H]/d[C2H5O] */
        J[850] += dqdci;              /* dwdot[CH3]/d[C2H5O] */
        J[861] -= dqdci;              /* dwdot[TXCH2]/d[C2H5O] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[21];
        dqdc[2] = TB[3][0]*dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = TB[3][1]*dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = TB[3][3]*dcdc_fac;
        dqdc[10] = dcdc_fac - k_r;
        dqdc[11] = dcdc_fac;
        dqdc[12] = TB[3][7]*dcdc_fac;
        dqdc[13] = TB[3][9]*dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = TB[3][16]*dcdc_fac;
        dqdc[18] = TB[3][2]*dcdc_fac;
        dqdc[19] = TB[3][4]*dcdc_fac;
        dqdc[20] = TB[3][5]*dcdc_fac;
        dqdc[21] = TB[3][6]*dcdc_fac + k_f*sc[1];
        dqdc[22] = TB[3][8]*dcdc_fac;
        dqdc[23] = TB[3][10]*dcdc_fac;
        dqdc[24] = TB[3][11]*dcdc_fac;
        dqdc[25] = TB[3][12]*dcdc_fac;
        dqdc[26] = TB[3][13]*dcdc_fac;
        dqdc[27] = TB[3][14]*dcdc_fac;
        dqdc[28] = TB[3][15]*dcdc_fac;
        for (int k=0; k<29; k++) {
            J[30*k+1] -= dqdc[k];
            J[30*k+10] += dqdc[k];
            J[30*k+21] -= dqdc[k];
        }
    }
    J[871] -= dqdT; /* dwdot[H]/dT */
    J[880] += dqdT; /* dwdot[CH3]/dT */
    J[891] -= dqdT; /* dwdot[TXCH2]/dT */

    /*reaction 5: CH3 + H (+M) <=> CH4 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[4][0] - 1)*sc[2] + (TB[4][1] - 1)*sc[7] + (TB[4][2] - 1)*sc[18] + (TB[4][3] - 1)*sc[9] + (TB[4][4] - 1)*sc[19] + (TB[4][5] - 1)*sc[20] + (TB[4][6] - 1)*sc[21] + (TB[4][7] - 1)*sc[12] + (TB[4][8] - 1)*sc[22] + (TB[4][9] - 1)*sc[13] + (TB[4][10] - 1)*sc[23] + (TB[4][11] - 1)*sc[24] + (TB[4][12] - 1)*sc[25] + (TB[4][13] - 1)*sc[26] + (TB[4][14] - 1)*sc[27] + (TB[4][15] - 1)*sc[28] + (TB[4][16] - 1)*sc[17];
    /* forward */
    phi_f = sc[1]*sc[10];
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
    wdot[1] -= q; /* H */
    wdot[10] -= q; /* CH3 */
    wdot[13] += q; /* CH4 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[10];
        J[31] -= dqdci;               /* dwdot[H]/d[H] */
        J[40] -= dqdci;               /* dwdot[CH3]/d[H] */
        J[43] += dqdci;               /* dwdot[CH4]/d[H] */
        /* d()/d[H2] */
        dqdci = (TB[4][0] - 1)*dcdc_fac;
        J[61] -= dqdci;               /* dwdot[H]/d[H2] */
        J[70] -= dqdci;               /* dwdot[CH3]/d[H2] */
        J[73] += dqdci;               /* dwdot[CH4]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[4][1] - 1)*dcdc_fac;
        J[211] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[220] -= dqdci;              /* dwdot[CH3]/d[H2O] */
        J[223] += dqdci;              /* dwdot[CH4]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[4][3] - 1)*dcdc_fac;
        J[271] -= dqdci;              /* dwdot[H]/d[CO] */
        J[280] -= dqdci;              /* dwdot[CH3]/d[CO] */
        J[283] += dqdci;              /* dwdot[CH4]/d[CO] */
        /* d()/d[CH3] */
        dqdci =  + k_f*sc[1];
        J[301] -= dqdci;              /* dwdot[H]/d[CH3] */
        J[310] -= dqdci;              /* dwdot[CH3]/d[CH3] */
        J[313] += dqdci;              /* dwdot[CH4]/d[CH3] */
        /* d()/d[CO2] */
        dqdci = (TB[4][7] - 1)*dcdc_fac;
        J[361] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[370] -= dqdci;              /* dwdot[CH3]/d[CO2] */
        J[373] += dqdci;              /* dwdot[CH4]/d[CO2] */
        /* d()/d[CH4] */
        dqdci = (TB[4][9] - 1)*dcdc_fac - k_r;
        J[391] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[400] -= dqdci;              /* dwdot[CH3]/d[CH4] */
        J[403] += dqdci;              /* dwdot[CH4]/d[CH4] */
        /* d()/d[C2H6] */
        dqdci = (TB[4][16] - 1)*dcdc_fac;
        J[511] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[520] -= dqdci;              /* dwdot[CH3]/d[C2H6] */
        J[523] += dqdci;              /* dwdot[CH4]/d[C2H6] */
        /* d()/d[C] */
        dqdci = (TB[4][2] - 1)*dcdc_fac;
        J[541] -= dqdci;              /* dwdot[H]/d[C] */
        J[550] -= dqdci;              /* dwdot[CH3]/d[C] */
        J[553] += dqdci;              /* dwdot[CH4]/d[C] */
        /* d()/d[CH] */
        dqdci = (TB[4][4] - 1)*dcdc_fac;
        J[571] -= dqdci;              /* dwdot[H]/d[CH] */
        J[580] -= dqdci;              /* dwdot[CH3]/d[CH] */
        J[583] += dqdci;              /* dwdot[CH4]/d[CH] */
        /* d()/d[HCO] */
        dqdci = (TB[4][5] - 1)*dcdc_fac;
        J[601] -= dqdci;              /* dwdot[H]/d[HCO] */
        J[610] -= dqdci;              /* dwdot[CH3]/d[HCO] */
        J[613] += dqdci;              /* dwdot[CH4]/d[HCO] */
        /* d()/d[TXCH2] */
        dqdci = (TB[4][6] - 1)*dcdc_fac;
        J[631] -= dqdci;              /* dwdot[H]/d[TXCH2] */
        J[640] -= dqdci;              /* dwdot[CH3]/d[TXCH2] */
        J[643] += dqdci;              /* dwdot[CH4]/d[TXCH2] */
        /* d()/d[SXCH2] */
        dqdci = (TB[4][8] - 1)*dcdc_fac;
        J[661] -= dqdci;              /* dwdot[H]/d[SXCH2] */
        J[670] -= dqdci;              /* dwdot[CH3]/d[SXCH2] */
        J[673] += dqdci;              /* dwdot[CH4]/d[SXCH2] */
        /* d()/d[C2H3] */
        dqdci = (TB[4][10] - 1)*dcdc_fac;
        J[691] -= dqdci;              /* dwdot[H]/d[C2H3] */
        J[700] -= dqdci;              /* dwdot[CH3]/d[C2H3] */
        J[703] += dqdci;              /* dwdot[CH4]/d[C2H3] */
        /* d()/d[C2H5] */
        dqdci = (TB[4][11] - 1)*dcdc_fac;
        J[721] -= dqdci;              /* dwdot[H]/d[C2H5] */
        J[730] -= dqdci;              /* dwdot[CH3]/d[C2H5] */
        J[733] += dqdci;              /* dwdot[CH4]/d[C2H5] */
        /* d()/d[HCCO] */
        dqdci = (TB[4][12] - 1)*dcdc_fac;
        J[751] -= dqdci;              /* dwdot[H]/d[HCCO] */
        J[760] -= dqdci;              /* dwdot[CH3]/d[HCCO] */
        J[763] += dqdci;              /* dwdot[CH4]/d[HCCO] */
        /* d()/d[CH3CHO] */
        dqdci = (TB[4][13] - 1)*dcdc_fac;
        J[781] -= dqdci;              /* dwdot[H]/d[CH3CHO] */
        J[790] -= dqdci;              /* dwdot[CH3]/d[CH3CHO] */
        J[793] += dqdci;              /* dwdot[CH4]/d[CH3CHO] */
        /* d()/d[CH2CHO] */
        dqdci = (TB[4][14] - 1)*dcdc_fac;
        J[811] -= dqdci;              /* dwdot[H]/d[CH2CHO] */
        J[820] -= dqdci;              /* dwdot[CH3]/d[CH2CHO] */
        J[823] += dqdci;              /* dwdot[CH4]/d[CH2CHO] */
        /* d()/d[C2H5O] */
        dqdci = (TB[4][15] - 1)*dcdc_fac;
        J[841] -= dqdci;              /* dwdot[H]/d[C2H5O] */
        J[850] -= dqdci;              /* dwdot[CH3]/d[C2H5O] */
        J[853] += dqdci;              /* dwdot[CH4]/d[C2H5O] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[10];
        dqdc[2] = TB[4][0]*dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = TB[4][1]*dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = TB[4][3]*dcdc_fac;
        dqdc[10] = dcdc_fac + k_f*sc[1];
        dqdc[11] = dcdc_fac;
        dqdc[12] = TB[4][7]*dcdc_fac;
        dqdc[13] = TB[4][9]*dcdc_fac - k_r;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = TB[4][16]*dcdc_fac;
        dqdc[18] = TB[4][2]*dcdc_fac;
        dqdc[19] = TB[4][4]*dcdc_fac;
        dqdc[20] = TB[4][5]*dcdc_fac;
        dqdc[21] = TB[4][6]*dcdc_fac;
        dqdc[22] = TB[4][8]*dcdc_fac;
        dqdc[23] = TB[4][10]*dcdc_fac;
        dqdc[24] = TB[4][11]*dcdc_fac;
        dqdc[25] = TB[4][12]*dcdc_fac;
        dqdc[26] = TB[4][13]*dcdc_fac;
        dqdc[27] = TB[4][14]*dcdc_fac;
        dqdc[28] = TB[4][15]*dcdc_fac;
        for (int k=0; k<29; k++) {
            J[30*k+1] -= dqdc[k];
            J[30*k+10] -= dqdc[k];
            J[30*k+13] += dqdc[k];
        }
    }
    J[871] -= dqdT; /* dwdot[H]/dT */
    J[880] -= dqdT; /* dwdot[CH3]/dT */
    J[883] += dqdT; /* dwdot[CH4]/dT */

    /*reaction 6: TXCH2 + CO (+M) <=> CH2CO (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[5][0] - 1)*sc[2] + (TB[5][1] - 1)*sc[7] + (TB[5][2] - 1)*sc[18] + (TB[5][3] - 1)*sc[9] + (TB[5][4] - 1)*sc[19] + (TB[5][5] - 1)*sc[20] + (TB[5][6] - 1)*sc[21] + (TB[5][7] - 1)*sc[12] + (TB[5][8] - 1)*sc[22] + (TB[5][9] - 1)*sc[13] + (TB[5][10] - 1)*sc[23] + (TB[5][11] - 1)*sc[24] + (TB[5][12] - 1)*sc[25] + (TB[5][13] - 1)*sc[26] + (TB[5][14] - 1)*sc[27] + (TB[5][15] - 1)*sc[28] + (TB[5][16] - 1)*sc[17];
    /* forward */
    phi_f = sc[9]*sc[21];
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
    phi_r = sc[16];
    Kc = refCinv * exp(g_RT[9] - g_RT[16] + g_RT[21]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[9] + h_RT[21]) + (h_RT[16]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[9] -= q; /* CO */
    wdot[16] += q; /* CH2CO */
    wdot[21] -= q; /* TXCH2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[5][0] - 1)*dcdc_fac;
        J[69] -= dqdci;               /* dwdot[CO]/d[H2] */
        J[76] += dqdci;               /* dwdot[CH2CO]/d[H2] */
        J[81] -= dqdci;               /* dwdot[TXCH2]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[5][1] - 1)*dcdc_fac;
        J[219] -= dqdci;              /* dwdot[CO]/d[H2O] */
        J[226] += dqdci;              /* dwdot[CH2CO]/d[H2O] */
        J[231] -= dqdci;              /* dwdot[TXCH2]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[5][3] - 1)*dcdc_fac + k_f*sc[21];
        J[279] -= dqdci;              /* dwdot[CO]/d[CO] */
        J[286] += dqdci;              /* dwdot[CH2CO]/d[CO] */
        J[291] -= dqdci;              /* dwdot[TXCH2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[5][7] - 1)*dcdc_fac;
        J[369] -= dqdci;              /* dwdot[CO]/d[CO2] */
        J[376] += dqdci;              /* dwdot[CH2CO]/d[CO2] */
        J[381] -= dqdci;              /* dwdot[TXCH2]/d[CO2] */
        /* d()/d[CH4] */
        dqdci = (TB[5][9] - 1)*dcdc_fac;
        J[399] -= dqdci;              /* dwdot[CO]/d[CH4] */
        J[406] += dqdci;              /* dwdot[CH2CO]/d[CH4] */
        J[411] -= dqdci;              /* dwdot[TXCH2]/d[CH4] */
        /* d()/d[CH2CO] */
        dqdci =  - k_r;
        J[489] -= dqdci;              /* dwdot[CO]/d[CH2CO] */
        J[496] += dqdci;              /* dwdot[CH2CO]/d[CH2CO] */
        J[501] -= dqdci;              /* dwdot[TXCH2]/d[CH2CO] */
        /* d()/d[C2H6] */
        dqdci = (TB[5][16] - 1)*dcdc_fac;
        J[519] -= dqdci;              /* dwdot[CO]/d[C2H6] */
        J[526] += dqdci;              /* dwdot[CH2CO]/d[C2H6] */
        J[531] -= dqdci;              /* dwdot[TXCH2]/d[C2H6] */
        /* d()/d[C] */
        dqdci = (TB[5][2] - 1)*dcdc_fac;
        J[549] -= dqdci;              /* dwdot[CO]/d[C] */
        J[556] += dqdci;              /* dwdot[CH2CO]/d[C] */
        J[561] -= dqdci;              /* dwdot[TXCH2]/d[C] */
        /* d()/d[CH] */
        dqdci = (TB[5][4] - 1)*dcdc_fac;
        J[579] -= dqdci;              /* dwdot[CO]/d[CH] */
        J[586] += dqdci;              /* dwdot[CH2CO]/d[CH] */
        J[591] -= dqdci;              /* dwdot[TXCH2]/d[CH] */
        /* d()/d[HCO] */
        dqdci = (TB[5][5] - 1)*dcdc_fac;
        J[609] -= dqdci;              /* dwdot[CO]/d[HCO] */
        J[616] += dqdci;              /* dwdot[CH2CO]/d[HCO] */
        J[621] -= dqdci;              /* dwdot[TXCH2]/d[HCO] */
        /* d()/d[TXCH2] */
        dqdci = (TB[5][6] - 1)*dcdc_fac + k_f*sc[9];
        J[639] -= dqdci;              /* dwdot[CO]/d[TXCH2] */
        J[646] += dqdci;              /* dwdot[CH2CO]/d[TXCH2] */
        J[651] -= dqdci;              /* dwdot[TXCH2]/d[TXCH2] */
        /* d()/d[SXCH2] */
        dqdci = (TB[5][8] - 1)*dcdc_fac;
        J[669] -= dqdci;              /* dwdot[CO]/d[SXCH2] */
        J[676] += dqdci;              /* dwdot[CH2CO]/d[SXCH2] */
        J[681] -= dqdci;              /* dwdot[TXCH2]/d[SXCH2] */
        /* d()/d[C2H3] */
        dqdci = (TB[5][10] - 1)*dcdc_fac;
        J[699] -= dqdci;              /* dwdot[CO]/d[C2H3] */
        J[706] += dqdci;              /* dwdot[CH2CO]/d[C2H3] */
        J[711] -= dqdci;              /* dwdot[TXCH2]/d[C2H3] */
        /* d()/d[C2H5] */
        dqdci = (TB[5][11] - 1)*dcdc_fac;
        J[729] -= dqdci;              /* dwdot[CO]/d[C2H5] */
        J[736] += dqdci;              /* dwdot[CH2CO]/d[C2H5] */
        J[741] -= dqdci;              /* dwdot[TXCH2]/d[C2H5] */
        /* d()/d[HCCO] */
        dqdci = (TB[5][12] - 1)*dcdc_fac;
        J[759] -= dqdci;              /* dwdot[CO]/d[HCCO] */
        J[766] += dqdci;              /* dwdot[CH2CO]/d[HCCO] */
        J[771] -= dqdci;              /* dwdot[TXCH2]/d[HCCO] */
        /* d()/d[CH3CHO] */
        dqdci = (TB[5][13] - 1)*dcdc_fac;
        J[789] -= dqdci;              /* dwdot[CO]/d[CH3CHO] */
        J[796] += dqdci;              /* dwdot[CH2CO]/d[CH3CHO] */
        J[801] -= dqdci;              /* dwdot[TXCH2]/d[CH3CHO] */
        /* d()/d[CH2CHO] */
        dqdci = (TB[5][14] - 1)*dcdc_fac;
        J[819] -= dqdci;              /* dwdot[CO]/d[CH2CHO] */
        J[826] += dqdci;              /* dwdot[CH2CO]/d[CH2CHO] */
        J[831] -= dqdci;              /* dwdot[TXCH2]/d[CH2CHO] */
        /* d()/d[C2H5O] */
        dqdci = (TB[5][15] - 1)*dcdc_fac;
        J[849] -= dqdci;              /* dwdot[CO]/d[C2H5O] */
        J[856] += dqdci;              /* dwdot[CH2CO]/d[C2H5O] */
        J[861] -= dqdci;              /* dwdot[TXCH2]/d[C2H5O] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = TB[5][0]*dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = TB[5][1]*dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = TB[5][3]*dcdc_fac + k_f*sc[21];
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = TB[5][7]*dcdc_fac;
        dqdc[13] = TB[5][9]*dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac - k_r;
        dqdc[17] = TB[5][16]*dcdc_fac;
        dqdc[18] = TB[5][2]*dcdc_fac;
        dqdc[19] = TB[5][4]*dcdc_fac;
        dqdc[20] = TB[5][5]*dcdc_fac;
        dqdc[21] = TB[5][6]*dcdc_fac + k_f*sc[9];
        dqdc[22] = TB[5][8]*dcdc_fac;
        dqdc[23] = TB[5][10]*dcdc_fac;
        dqdc[24] = TB[5][11]*dcdc_fac;
        dqdc[25] = TB[5][12]*dcdc_fac;
        dqdc[26] = TB[5][13]*dcdc_fac;
        dqdc[27] = TB[5][14]*dcdc_fac;
        dqdc[28] = TB[5][15]*dcdc_fac;
        for (int k=0; k<29; k++) {
            J[30*k+9] -= dqdc[k];
            J[30*k+16] += dqdc[k];
            J[30*k+21] -= dqdc[k];
        }
    }
    J[879] -= dqdT; /* dwdot[CO]/dT */
    J[886] += dqdT; /* dwdot[CH2CO]/dT */
    J[891] -= dqdT; /* dwdot[TXCH2]/dT */

    /*reaction 7: CO + H2 (+M) <=> CH2O (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[6][0] - 1)*sc[2] + (TB[6][1] - 1)*sc[7] + (TB[6][2] - 1)*sc[18] + (TB[6][3] - 1)*sc[9] + (TB[6][4] - 1)*sc[19] + (TB[6][5] - 1)*sc[20] + (TB[6][6] - 1)*sc[21] + (TB[6][7] - 1)*sc[12] + (TB[6][8] - 1)*sc[22] + (TB[6][9] - 1)*sc[13] + (TB[6][10] - 1)*sc[23] + (TB[6][11] - 1)*sc[24] + (TB[6][12] - 1)*sc[25] + (TB[6][13] - 1)*sc[26] + (TB[6][14] - 1)*sc[27] + (TB[6][15] - 1)*sc[28] + (TB[6][16] - 1)*sc[17];
    /* forward */
    phi_f = sc[2]*sc[9];
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
    phi_r = sc[11];
    Kc = refCinv * exp(g_RT[2] + g_RT[9] - g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[9]) + (h_RT[11]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[2] -= q; /* H2 */
    wdot[9] -= q; /* CO */
    wdot[11] += q; /* CH2O */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[6][0] - 1)*dcdc_fac + k_f*sc[9];
        J[62] -= dqdci;               /* dwdot[H2]/d[H2] */
        J[69] -= dqdci;               /* dwdot[CO]/d[H2] */
        J[71] += dqdci;               /* dwdot[CH2O]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[6][1] - 1)*dcdc_fac;
        J[212] -= dqdci;              /* dwdot[H2]/d[H2O] */
        J[219] -= dqdci;              /* dwdot[CO]/d[H2O] */
        J[221] += dqdci;              /* dwdot[CH2O]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[6][3] - 1)*dcdc_fac + k_f*sc[2];
        J[272] -= dqdci;              /* dwdot[H2]/d[CO] */
        J[279] -= dqdci;              /* dwdot[CO]/d[CO] */
        J[281] += dqdci;              /* dwdot[CH2O]/d[CO] */
        /* d()/d[CH2O] */
        dqdci =  - k_r;
        J[332] -= dqdci;              /* dwdot[H2]/d[CH2O] */
        J[339] -= dqdci;              /* dwdot[CO]/d[CH2O] */
        J[341] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
        /* d()/d[CO2] */
        dqdci = (TB[6][7] - 1)*dcdc_fac;
        J[362] -= dqdci;              /* dwdot[H2]/d[CO2] */
        J[369] -= dqdci;              /* dwdot[CO]/d[CO2] */
        J[371] += dqdci;              /* dwdot[CH2O]/d[CO2] */
        /* d()/d[CH4] */
        dqdci = (TB[6][9] - 1)*dcdc_fac;
        J[392] -= dqdci;              /* dwdot[H2]/d[CH4] */
        J[399] -= dqdci;              /* dwdot[CO]/d[CH4] */
        J[401] += dqdci;              /* dwdot[CH2O]/d[CH4] */
        /* d()/d[C2H6] */
        dqdci = (TB[6][16] - 1)*dcdc_fac;
        J[512] -= dqdci;              /* dwdot[H2]/d[C2H6] */
        J[519] -= dqdci;              /* dwdot[CO]/d[C2H6] */
        J[521] += dqdci;              /* dwdot[CH2O]/d[C2H6] */
        /* d()/d[C] */
        dqdci = (TB[6][2] - 1)*dcdc_fac;
        J[542] -= dqdci;              /* dwdot[H2]/d[C] */
        J[549] -= dqdci;              /* dwdot[CO]/d[C] */
        J[551] += dqdci;              /* dwdot[CH2O]/d[C] */
        /* d()/d[CH] */
        dqdci = (TB[6][4] - 1)*dcdc_fac;
        J[572] -= dqdci;              /* dwdot[H2]/d[CH] */
        J[579] -= dqdci;              /* dwdot[CO]/d[CH] */
        J[581] += dqdci;              /* dwdot[CH2O]/d[CH] */
        /* d()/d[HCO] */
        dqdci = (TB[6][5] - 1)*dcdc_fac;
        J[602] -= dqdci;              /* dwdot[H2]/d[HCO] */
        J[609] -= dqdci;              /* dwdot[CO]/d[HCO] */
        J[611] += dqdci;              /* dwdot[CH2O]/d[HCO] */
        /* d()/d[TXCH2] */
        dqdci = (TB[6][6] - 1)*dcdc_fac;
        J[632] -= dqdci;              /* dwdot[H2]/d[TXCH2] */
        J[639] -= dqdci;              /* dwdot[CO]/d[TXCH2] */
        J[641] += dqdci;              /* dwdot[CH2O]/d[TXCH2] */
        /* d()/d[SXCH2] */
        dqdci = (TB[6][8] - 1)*dcdc_fac;
        J[662] -= dqdci;              /* dwdot[H2]/d[SXCH2] */
        J[669] -= dqdci;              /* dwdot[CO]/d[SXCH2] */
        J[671] += dqdci;              /* dwdot[CH2O]/d[SXCH2] */
        /* d()/d[C2H3] */
        dqdci = (TB[6][10] - 1)*dcdc_fac;
        J[692] -= dqdci;              /* dwdot[H2]/d[C2H3] */
        J[699] -= dqdci;              /* dwdot[CO]/d[C2H3] */
        J[701] += dqdci;              /* dwdot[CH2O]/d[C2H3] */
        /* d()/d[C2H5] */
        dqdci = (TB[6][11] - 1)*dcdc_fac;
        J[722] -= dqdci;              /* dwdot[H2]/d[C2H5] */
        J[729] -= dqdci;              /* dwdot[CO]/d[C2H5] */
        J[731] += dqdci;              /* dwdot[CH2O]/d[C2H5] */
        /* d()/d[HCCO] */
        dqdci = (TB[6][12] - 1)*dcdc_fac;
        J[752] -= dqdci;              /* dwdot[H2]/d[HCCO] */
        J[759] -= dqdci;              /* dwdot[CO]/d[HCCO] */
        J[761] += dqdci;              /* dwdot[CH2O]/d[HCCO] */
        /* d()/d[CH3CHO] */
        dqdci = (TB[6][13] - 1)*dcdc_fac;
        J[782] -= dqdci;              /* dwdot[H2]/d[CH3CHO] */
        J[789] -= dqdci;              /* dwdot[CO]/d[CH3CHO] */
        J[791] += dqdci;              /* dwdot[CH2O]/d[CH3CHO] */
        /* d()/d[CH2CHO] */
        dqdci = (TB[6][14] - 1)*dcdc_fac;
        J[812] -= dqdci;              /* dwdot[H2]/d[CH2CHO] */
        J[819] -= dqdci;              /* dwdot[CO]/d[CH2CHO] */
        J[821] += dqdci;              /* dwdot[CH2O]/d[CH2CHO] */
        /* d()/d[C2H5O] */
        dqdci = (TB[6][15] - 1)*dcdc_fac;
        J[842] -= dqdci;              /* dwdot[H2]/d[C2H5O] */
        J[849] -= dqdci;              /* dwdot[CO]/d[C2H5O] */
        J[851] += dqdci;              /* dwdot[CH2O]/d[C2H5O] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = TB[6][0]*dcdc_fac + k_f*sc[9];
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = TB[6][1]*dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = TB[6][3]*dcdc_fac + k_f*sc[2];
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac - k_r;
        dqdc[12] = TB[6][7]*dcdc_fac;
        dqdc[13] = TB[6][9]*dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = TB[6][16]*dcdc_fac;
        dqdc[18] = TB[6][2]*dcdc_fac;
        dqdc[19] = TB[6][4]*dcdc_fac;
        dqdc[20] = TB[6][5]*dcdc_fac;
        dqdc[21] = TB[6][6]*dcdc_fac;
        dqdc[22] = TB[6][8]*dcdc_fac;
        dqdc[23] = TB[6][10]*dcdc_fac;
        dqdc[24] = TB[6][11]*dcdc_fac;
        dqdc[25] = TB[6][12]*dcdc_fac;
        dqdc[26] = TB[6][13]*dcdc_fac;
        dqdc[27] = TB[6][14]*dcdc_fac;
        dqdc[28] = TB[6][15]*dcdc_fac;
        for (int k=0; k<29; k++) {
            J[30*k+2] -= dqdc[k];
            J[30*k+9] -= dqdc[k];
            J[30*k+11] += dqdc[k];
        }
    }
    J[872] -= dqdT; /* dwdot[H2]/dT */
    J[879] -= dqdT; /* dwdot[CO]/dT */
    J[881] += dqdT; /* dwdot[CH2O]/dT */

    /*reaction 8: CH + CO (+M) <=> HCCO (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[7][0] - 1)*sc[2] + (TB[7][1] - 1)*sc[7] + (TB[7][2] - 1)*sc[18] + (TB[7][3] - 1)*sc[9] + (TB[7][4] - 1)*sc[19] + (TB[7][5] - 1)*sc[20] + (TB[7][6] - 1)*sc[21] + (TB[7][7] - 1)*sc[12] + (TB[7][8] - 1)*sc[22] + (TB[7][9] - 1)*sc[13] + (TB[7][10] - 1)*sc[23] + (TB[7][11] - 1)*sc[24] + (TB[7][12] - 1)*sc[25] + (TB[7][13] - 1)*sc[26] + (TB[7][14] - 1)*sc[27] + (TB[7][15] - 1)*sc[28] + (TB[7][16] - 1)*sc[17];
    /* forward */
    phi_f = sc[9]*sc[19];
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
    phi_r = sc[25];
    Kc = refCinv * exp(g_RT[9] + g_RT[19] - g_RT[25]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[9] + h_RT[19]) + (h_RT[25]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[9] -= q; /* CO */
    wdot[19] -= q; /* CH */
    wdot[25] += q; /* HCCO */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[7][0] - 1)*dcdc_fac;
        J[69] -= dqdci;               /* dwdot[CO]/d[H2] */
        J[79] -= dqdci;               /* dwdot[CH]/d[H2] */
        J[85] += dqdci;               /* dwdot[HCCO]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[7][1] - 1)*dcdc_fac;
        J[219] -= dqdci;              /* dwdot[CO]/d[H2O] */
        J[229] -= dqdci;              /* dwdot[CH]/d[H2O] */
        J[235] += dqdci;              /* dwdot[HCCO]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[7][3] - 1)*dcdc_fac + k_f*sc[19];
        J[279] -= dqdci;              /* dwdot[CO]/d[CO] */
        J[289] -= dqdci;              /* dwdot[CH]/d[CO] */
        J[295] += dqdci;              /* dwdot[HCCO]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[7][7] - 1)*dcdc_fac;
        J[369] -= dqdci;              /* dwdot[CO]/d[CO2] */
        J[379] -= dqdci;              /* dwdot[CH]/d[CO2] */
        J[385] += dqdci;              /* dwdot[HCCO]/d[CO2] */
        /* d()/d[CH4] */
        dqdci = (TB[7][9] - 1)*dcdc_fac;
        J[399] -= dqdci;              /* dwdot[CO]/d[CH4] */
        J[409] -= dqdci;              /* dwdot[CH]/d[CH4] */
        J[415] += dqdci;              /* dwdot[HCCO]/d[CH4] */
        /* d()/d[C2H6] */
        dqdci = (TB[7][16] - 1)*dcdc_fac;
        J[519] -= dqdci;              /* dwdot[CO]/d[C2H6] */
        J[529] -= dqdci;              /* dwdot[CH]/d[C2H6] */
        J[535] += dqdci;              /* dwdot[HCCO]/d[C2H6] */
        /* d()/d[C] */
        dqdci = (TB[7][2] - 1)*dcdc_fac;
        J[549] -= dqdci;              /* dwdot[CO]/d[C] */
        J[559] -= dqdci;              /* dwdot[CH]/d[C] */
        J[565] += dqdci;              /* dwdot[HCCO]/d[C] */
        /* d()/d[CH] */
        dqdci = (TB[7][4] - 1)*dcdc_fac + k_f*sc[9];
        J[579] -= dqdci;              /* dwdot[CO]/d[CH] */
        J[589] -= dqdci;              /* dwdot[CH]/d[CH] */
        J[595] += dqdci;              /* dwdot[HCCO]/d[CH] */
        /* d()/d[HCO] */
        dqdci = (TB[7][5] - 1)*dcdc_fac;
        J[609] -= dqdci;              /* dwdot[CO]/d[HCO] */
        J[619] -= dqdci;              /* dwdot[CH]/d[HCO] */
        J[625] += dqdci;              /* dwdot[HCCO]/d[HCO] */
        /* d()/d[TXCH2] */
        dqdci = (TB[7][6] - 1)*dcdc_fac;
        J[639] -= dqdci;              /* dwdot[CO]/d[TXCH2] */
        J[649] -= dqdci;              /* dwdot[CH]/d[TXCH2] */
        J[655] += dqdci;              /* dwdot[HCCO]/d[TXCH2] */
        /* d()/d[SXCH2] */
        dqdci = (TB[7][8] - 1)*dcdc_fac;
        J[669] -= dqdci;              /* dwdot[CO]/d[SXCH2] */
        J[679] -= dqdci;              /* dwdot[CH]/d[SXCH2] */
        J[685] += dqdci;              /* dwdot[HCCO]/d[SXCH2] */
        /* d()/d[C2H3] */
        dqdci = (TB[7][10] - 1)*dcdc_fac;
        J[699] -= dqdci;              /* dwdot[CO]/d[C2H3] */
        J[709] -= dqdci;              /* dwdot[CH]/d[C2H3] */
        J[715] += dqdci;              /* dwdot[HCCO]/d[C2H3] */
        /* d()/d[C2H5] */
        dqdci = (TB[7][11] - 1)*dcdc_fac;
        J[729] -= dqdci;              /* dwdot[CO]/d[C2H5] */
        J[739] -= dqdci;              /* dwdot[CH]/d[C2H5] */
        J[745] += dqdci;              /* dwdot[HCCO]/d[C2H5] */
        /* d()/d[HCCO] */
        dqdci = (TB[7][12] - 1)*dcdc_fac - k_r;
        J[759] -= dqdci;              /* dwdot[CO]/d[HCCO] */
        J[769] -= dqdci;              /* dwdot[CH]/d[HCCO] */
        J[775] += dqdci;              /* dwdot[HCCO]/d[HCCO] */
        /* d()/d[CH3CHO] */
        dqdci = (TB[7][13] - 1)*dcdc_fac;
        J[789] -= dqdci;              /* dwdot[CO]/d[CH3CHO] */
        J[799] -= dqdci;              /* dwdot[CH]/d[CH3CHO] */
        J[805] += dqdci;              /* dwdot[HCCO]/d[CH3CHO] */
        /* d()/d[CH2CHO] */
        dqdci = (TB[7][14] - 1)*dcdc_fac;
        J[819] -= dqdci;              /* dwdot[CO]/d[CH2CHO] */
        J[829] -= dqdci;              /* dwdot[CH]/d[CH2CHO] */
        J[835] += dqdci;              /* dwdot[HCCO]/d[CH2CHO] */
        /* d()/d[C2H5O] */
        dqdci = (TB[7][15] - 1)*dcdc_fac;
        J[849] -= dqdci;              /* dwdot[CO]/d[C2H5O] */
        J[859] -= dqdci;              /* dwdot[CH]/d[C2H5O] */
        J[865] += dqdci;              /* dwdot[HCCO]/d[C2H5O] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = TB[7][0]*dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = TB[7][1]*dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = TB[7][3]*dcdc_fac + k_f*sc[19];
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = TB[7][7]*dcdc_fac;
        dqdc[13] = TB[7][9]*dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = TB[7][16]*dcdc_fac;
        dqdc[18] = TB[7][2]*dcdc_fac;
        dqdc[19] = TB[7][4]*dcdc_fac + k_f*sc[9];
        dqdc[20] = TB[7][5]*dcdc_fac;
        dqdc[21] = TB[7][6]*dcdc_fac;
        dqdc[22] = TB[7][8]*dcdc_fac;
        dqdc[23] = TB[7][10]*dcdc_fac;
        dqdc[24] = TB[7][11]*dcdc_fac;
        dqdc[25] = TB[7][12]*dcdc_fac - k_r;
        dqdc[26] = TB[7][13]*dcdc_fac;
        dqdc[27] = TB[7][14]*dcdc_fac;
        dqdc[28] = TB[7][15]*dcdc_fac;
        for (int k=0; k<29; k++) {
            J[30*k+9] -= dqdc[k];
            J[30*k+19] -= dqdc[k];
            J[30*k+25] += dqdc[k];
        }
    }
    J[879] -= dqdT; /* dwdot[CO]/dT */
    J[889] -= dqdT; /* dwdot[CH]/dT */
    J[895] += dqdT; /* dwdot[HCCO]/dT */

    /*reaction 9: CO + O (+M) <=> CO2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[8][0] - 1)*sc[2] + (TB[8][1] - 1)*sc[7] + (TB[8][2] - 1)*sc[18] + (TB[8][3] - 1)*sc[9] + (TB[8][4] - 1)*sc[19] + (TB[8][5] - 1)*sc[20] + (TB[8][6] - 1)*sc[21] + (TB[8][7] - 1)*sc[12] + (TB[8][8] - 1)*sc[22] + (TB[8][9] - 1)*sc[13] + (TB[8][10] - 1)*sc[23] + (TB[8][11] - 1)*sc[24] + (TB[8][12] - 1)*sc[25] + (TB[8][13] - 1)*sc[26] + (TB[8][14] - 1)*sc[27] + (TB[8][15] - 1)*sc[28] + (TB[8][16] - 1)*sc[17];
    /* forward */
    phi_f = sc[3]*sc[9];
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
    phi_r = sc[12];
    Kc = refCinv * exp(g_RT[3] + g_RT[9] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[9]) + (h_RT[12]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[3] -= q; /* O */
    wdot[9] -= q; /* CO */
    wdot[12] += q; /* CO2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[8][0] - 1)*dcdc_fac;
        J[63] -= dqdci;               /* dwdot[O]/d[H2] */
        J[69] -= dqdci;               /* dwdot[CO]/d[H2] */
        J[72] += dqdci;               /* dwdot[CO2]/d[H2] */
        /* d()/d[O] */
        dqdci =  + k_f*sc[9];
        J[93] -= dqdci;               /* dwdot[O]/d[O] */
        J[99] -= dqdci;               /* dwdot[CO]/d[O] */
        J[102] += dqdci;              /* dwdot[CO2]/d[O] */
        /* d()/d[H2O] */
        dqdci = (TB[8][1] - 1)*dcdc_fac;
        J[213] -= dqdci;              /* dwdot[O]/d[H2O] */
        J[219] -= dqdci;              /* dwdot[CO]/d[H2O] */
        J[222] += dqdci;              /* dwdot[CO2]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[8][3] - 1)*dcdc_fac + k_f*sc[3];
        J[273] -= dqdci;              /* dwdot[O]/d[CO] */
        J[279] -= dqdci;              /* dwdot[CO]/d[CO] */
        J[282] += dqdci;              /* dwdot[CO2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[8][7] - 1)*dcdc_fac - k_r;
        J[363] -= dqdci;              /* dwdot[O]/d[CO2] */
        J[369] -= dqdci;              /* dwdot[CO]/d[CO2] */
        J[372] += dqdci;              /* dwdot[CO2]/d[CO2] */
        /* d()/d[CH4] */
        dqdci = (TB[8][9] - 1)*dcdc_fac;
        J[393] -= dqdci;              /* dwdot[O]/d[CH4] */
        J[399] -= dqdci;              /* dwdot[CO]/d[CH4] */
        J[402] += dqdci;              /* dwdot[CO2]/d[CH4] */
        /* d()/d[C2H6] */
        dqdci = (TB[8][16] - 1)*dcdc_fac;
        J[513] -= dqdci;              /* dwdot[O]/d[C2H6] */
        J[519] -= dqdci;              /* dwdot[CO]/d[C2H6] */
        J[522] += dqdci;              /* dwdot[CO2]/d[C2H6] */
        /* d()/d[C] */
        dqdci = (TB[8][2] - 1)*dcdc_fac;
        J[543] -= dqdci;              /* dwdot[O]/d[C] */
        J[549] -= dqdci;              /* dwdot[CO]/d[C] */
        J[552] += dqdci;              /* dwdot[CO2]/d[C] */
        /* d()/d[CH] */
        dqdci = (TB[8][4] - 1)*dcdc_fac;
        J[573] -= dqdci;              /* dwdot[O]/d[CH] */
        J[579] -= dqdci;              /* dwdot[CO]/d[CH] */
        J[582] += dqdci;              /* dwdot[CO2]/d[CH] */
        /* d()/d[HCO] */
        dqdci = (TB[8][5] - 1)*dcdc_fac;
        J[603] -= dqdci;              /* dwdot[O]/d[HCO] */
        J[609] -= dqdci;              /* dwdot[CO]/d[HCO] */
        J[612] += dqdci;              /* dwdot[CO2]/d[HCO] */
        /* d()/d[TXCH2] */
        dqdci = (TB[8][6] - 1)*dcdc_fac;
        J[633] -= dqdci;              /* dwdot[O]/d[TXCH2] */
        J[639] -= dqdci;              /* dwdot[CO]/d[TXCH2] */
        J[642] += dqdci;              /* dwdot[CO2]/d[TXCH2] */
        /* d()/d[SXCH2] */
        dqdci = (TB[8][8] - 1)*dcdc_fac;
        J[663] -= dqdci;              /* dwdot[O]/d[SXCH2] */
        J[669] -= dqdci;              /* dwdot[CO]/d[SXCH2] */
        J[672] += dqdci;              /* dwdot[CO2]/d[SXCH2] */
        /* d()/d[C2H3] */
        dqdci = (TB[8][10] - 1)*dcdc_fac;
        J[693] -= dqdci;              /* dwdot[O]/d[C2H3] */
        J[699] -= dqdci;              /* dwdot[CO]/d[C2H3] */
        J[702] += dqdci;              /* dwdot[CO2]/d[C2H3] */
        /* d()/d[C2H5] */
        dqdci = (TB[8][11] - 1)*dcdc_fac;
        J[723] -= dqdci;              /* dwdot[O]/d[C2H5] */
        J[729] -= dqdci;              /* dwdot[CO]/d[C2H5] */
        J[732] += dqdci;              /* dwdot[CO2]/d[C2H5] */
        /* d()/d[HCCO] */
        dqdci = (TB[8][12] - 1)*dcdc_fac;
        J[753] -= dqdci;              /* dwdot[O]/d[HCCO] */
        J[759] -= dqdci;              /* dwdot[CO]/d[HCCO] */
        J[762] += dqdci;              /* dwdot[CO2]/d[HCCO] */
        /* d()/d[CH3CHO] */
        dqdci = (TB[8][13] - 1)*dcdc_fac;
        J[783] -= dqdci;              /* dwdot[O]/d[CH3CHO] */
        J[789] -= dqdci;              /* dwdot[CO]/d[CH3CHO] */
        J[792] += dqdci;              /* dwdot[CO2]/d[CH3CHO] */
        /* d()/d[CH2CHO] */
        dqdci = (TB[8][14] - 1)*dcdc_fac;
        J[813] -= dqdci;              /* dwdot[O]/d[CH2CHO] */
        J[819] -= dqdci;              /* dwdot[CO]/d[CH2CHO] */
        J[822] += dqdci;              /* dwdot[CO2]/d[CH2CHO] */
        /* d()/d[C2H5O] */
        dqdci = (TB[8][15] - 1)*dcdc_fac;
        J[843] -= dqdci;              /* dwdot[O]/d[C2H5O] */
        J[849] -= dqdci;              /* dwdot[CO]/d[C2H5O] */
        J[852] += dqdci;              /* dwdot[CO2]/d[C2H5O] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = TB[8][0]*dcdc_fac;
        dqdc[3] = dcdc_fac + k_f*sc[9];
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = TB[8][1]*dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = TB[8][3]*dcdc_fac + k_f*sc[3];
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = TB[8][7]*dcdc_fac - k_r;
        dqdc[13] = TB[8][9]*dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = TB[8][16]*dcdc_fac;
        dqdc[18] = TB[8][2]*dcdc_fac;
        dqdc[19] = TB[8][4]*dcdc_fac;
        dqdc[20] = TB[8][5]*dcdc_fac;
        dqdc[21] = TB[8][6]*dcdc_fac;
        dqdc[22] = TB[8][8]*dcdc_fac;
        dqdc[23] = TB[8][10]*dcdc_fac;
        dqdc[24] = TB[8][11]*dcdc_fac;
        dqdc[25] = TB[8][12]*dcdc_fac;
        dqdc[26] = TB[8][13]*dcdc_fac;
        dqdc[27] = TB[8][14]*dcdc_fac;
        dqdc[28] = TB[8][15]*dcdc_fac;
        for (int k=0; k<29; k++) {
            J[30*k+3] -= dqdc[k];
            J[30*k+9] -= dqdc[k];
            J[30*k+12] += dqdc[k];
        }
    }
    J[873] -= dqdT; /* dwdot[O]/dT */
    J[879] -= dqdT; /* dwdot[CO]/dT */
    J[882] += dqdT; /* dwdot[CO2]/dT */

    /*reaction 10: HCO + H (+M) <=> CH2O (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[9][0] - 1)*sc[2] + (TB[9][1] - 1)*sc[7] + (TB[9][2] - 1)*sc[18] + (TB[9][3] - 1)*sc[9] + (TB[9][4] - 1)*sc[19] + (TB[9][5] - 1)*sc[20] + (TB[9][6] - 1)*sc[21] + (TB[9][7] - 1)*sc[12] + (TB[9][8] - 1)*sc[22] + (TB[9][9] - 1)*sc[13] + (TB[9][10] - 1)*sc[23] + (TB[9][11] - 1)*sc[24] + (TB[9][12] - 1)*sc[25] + (TB[9][13] - 1)*sc[26] + (TB[9][14] - 1)*sc[27] + (TB[9][15] - 1)*sc[28] + (TB[9][16] - 1)*sc[17];
    /* forward */
    phi_f = sc[1]*sc[20];
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
    phi_r = sc[11];
    Kc = refCinv * exp(g_RT[1] - g_RT[11] + g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[20]) + (h_RT[11]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[11] += q; /* CH2O */
    wdot[20] -= q; /* HCO */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[20];
        J[31] -= dqdci;               /* dwdot[H]/d[H] */
        J[41] += dqdci;               /* dwdot[CH2O]/d[H] */
        J[50] -= dqdci;               /* dwdot[HCO]/d[H] */
        /* d()/d[H2] */
        dqdci = (TB[9][0] - 1)*dcdc_fac;
        J[61] -= dqdci;               /* dwdot[H]/d[H2] */
        J[71] += dqdci;               /* dwdot[CH2O]/d[H2] */
        J[80] -= dqdci;               /* dwdot[HCO]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[9][1] - 1)*dcdc_fac;
        J[211] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[221] += dqdci;              /* dwdot[CH2O]/d[H2O] */
        J[230] -= dqdci;              /* dwdot[HCO]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[9][3] - 1)*dcdc_fac;
        J[271] -= dqdci;              /* dwdot[H]/d[CO] */
        J[281] += dqdci;              /* dwdot[CH2O]/d[CO] */
        J[290] -= dqdci;              /* dwdot[HCO]/d[CO] */
        /* d()/d[CH2O] */
        dqdci =  - k_r;
        J[331] -= dqdci;              /* dwdot[H]/d[CH2O] */
        J[341] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
        J[350] -= dqdci;              /* dwdot[HCO]/d[CH2O] */
        /* d()/d[CO2] */
        dqdci = (TB[9][7] - 1)*dcdc_fac;
        J[361] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[371] += dqdci;              /* dwdot[CH2O]/d[CO2] */
        J[380] -= dqdci;              /* dwdot[HCO]/d[CO2] */
        /* d()/d[CH4] */
        dqdci = (TB[9][9] - 1)*dcdc_fac;
        J[391] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[401] += dqdci;              /* dwdot[CH2O]/d[CH4] */
        J[410] -= dqdci;              /* dwdot[HCO]/d[CH4] */
        /* d()/d[C2H6] */
        dqdci = (TB[9][16] - 1)*dcdc_fac;
        J[511] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[521] += dqdci;              /* dwdot[CH2O]/d[C2H6] */
        J[530] -= dqdci;              /* dwdot[HCO]/d[C2H6] */
        /* d()/d[C] */
        dqdci = (TB[9][2] - 1)*dcdc_fac;
        J[541] -= dqdci;              /* dwdot[H]/d[C] */
        J[551] += dqdci;              /* dwdot[CH2O]/d[C] */
        J[560] -= dqdci;              /* dwdot[HCO]/d[C] */
        /* d()/d[CH] */
        dqdci = (TB[9][4] - 1)*dcdc_fac;
        J[571] -= dqdci;              /* dwdot[H]/d[CH] */
        J[581] += dqdci;              /* dwdot[CH2O]/d[CH] */
        J[590] -= dqdci;              /* dwdot[HCO]/d[CH] */
        /* d()/d[HCO] */
        dqdci = (TB[9][5] - 1)*dcdc_fac + k_f*sc[1];
        J[601] -= dqdci;              /* dwdot[H]/d[HCO] */
        J[611] += dqdci;              /* dwdot[CH2O]/d[HCO] */
        J[620] -= dqdci;              /* dwdot[HCO]/d[HCO] */
        /* d()/d[TXCH2] */
        dqdci = (TB[9][6] - 1)*dcdc_fac;
        J[631] -= dqdci;              /* dwdot[H]/d[TXCH2] */
        J[641] += dqdci;              /* dwdot[CH2O]/d[TXCH2] */
        J[650] -= dqdci;              /* dwdot[HCO]/d[TXCH2] */
        /* d()/d[SXCH2] */
        dqdci = (TB[9][8] - 1)*dcdc_fac;
        J[661] -= dqdci;              /* dwdot[H]/d[SXCH2] */
        J[671] += dqdci;              /* dwdot[CH2O]/d[SXCH2] */
        J[680] -= dqdci;              /* dwdot[HCO]/d[SXCH2] */
        /* d()/d[C2H3] */
        dqdci = (TB[9][10] - 1)*dcdc_fac;
        J[691] -= dqdci;              /* dwdot[H]/d[C2H3] */
        J[701] += dqdci;              /* dwdot[CH2O]/d[C2H3] */
        J[710] -= dqdci;              /* dwdot[HCO]/d[C2H3] */
        /* d()/d[C2H5] */
        dqdci = (TB[9][11] - 1)*dcdc_fac;
        J[721] -= dqdci;              /* dwdot[H]/d[C2H5] */
        J[731] += dqdci;              /* dwdot[CH2O]/d[C2H5] */
        J[740] -= dqdci;              /* dwdot[HCO]/d[C2H5] */
        /* d()/d[HCCO] */
        dqdci = (TB[9][12] - 1)*dcdc_fac;
        J[751] -= dqdci;              /* dwdot[H]/d[HCCO] */
        J[761] += dqdci;              /* dwdot[CH2O]/d[HCCO] */
        J[770] -= dqdci;              /* dwdot[HCO]/d[HCCO] */
        /* d()/d[CH3CHO] */
        dqdci = (TB[9][13] - 1)*dcdc_fac;
        J[781] -= dqdci;              /* dwdot[H]/d[CH3CHO] */
        J[791] += dqdci;              /* dwdot[CH2O]/d[CH3CHO] */
        J[800] -= dqdci;              /* dwdot[HCO]/d[CH3CHO] */
        /* d()/d[CH2CHO] */
        dqdci = (TB[9][14] - 1)*dcdc_fac;
        J[811] -= dqdci;              /* dwdot[H]/d[CH2CHO] */
        J[821] += dqdci;              /* dwdot[CH2O]/d[CH2CHO] */
        J[830] -= dqdci;              /* dwdot[HCO]/d[CH2CHO] */
        /* d()/d[C2H5O] */
        dqdci = (TB[9][15] - 1)*dcdc_fac;
        J[841] -= dqdci;              /* dwdot[H]/d[C2H5O] */
        J[851] += dqdci;              /* dwdot[CH2O]/d[C2H5O] */
        J[860] -= dqdci;              /* dwdot[HCO]/d[C2H5O] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[20];
        dqdc[2] = TB[9][0]*dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = TB[9][1]*dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = TB[9][3]*dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac - k_r;
        dqdc[12] = TB[9][7]*dcdc_fac;
        dqdc[13] = TB[9][9]*dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = TB[9][16]*dcdc_fac;
        dqdc[18] = TB[9][2]*dcdc_fac;
        dqdc[19] = TB[9][4]*dcdc_fac;
        dqdc[20] = TB[9][5]*dcdc_fac + k_f*sc[1];
        dqdc[21] = TB[9][6]*dcdc_fac;
        dqdc[22] = TB[9][8]*dcdc_fac;
        dqdc[23] = TB[9][10]*dcdc_fac;
        dqdc[24] = TB[9][11]*dcdc_fac;
        dqdc[25] = TB[9][12]*dcdc_fac;
        dqdc[26] = TB[9][13]*dcdc_fac;
        dqdc[27] = TB[9][14]*dcdc_fac;
        dqdc[28] = TB[9][15]*dcdc_fac;
        for (int k=0; k<29; k++) {
            J[30*k+1] -= dqdc[k];
            J[30*k+11] += dqdc[k];
            J[30*k+20] -= dqdc[k];
        }
    }
    J[871] -= dqdT; /* dwdot[H]/dT */
    J[881] += dqdT; /* dwdot[CH2O]/dT */
    J[890] -= dqdT; /* dwdot[HCO]/dT */

    /*reaction 11: C2H2 + H (+M) <=> C2H3 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[10][0] - 1)*sc[2] + (TB[10][1] - 1)*sc[7] + (TB[10][2] - 1)*sc[18] + (TB[10][3] - 1)*sc[9] + (TB[10][4] - 1)*sc[19] + (TB[10][5] - 1)*sc[20] + (TB[10][6] - 1)*sc[21] + (TB[10][7] - 1)*sc[12] + (TB[10][8] - 1)*sc[22] + (TB[10][9] - 1)*sc[13] + (TB[10][10] - 1)*sc[23] + (TB[10][11] - 1)*sc[24] + (TB[10][12] - 1)*sc[25] + (TB[10][13] - 1)*sc[26] + (TB[10][14] - 1)*sc[27] + (TB[10][15] - 1)*sc[28] + (TB[10][16] - 1)*sc[17];
    /* forward */
    phi_f = sc[1]*sc[14];
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
    phi_r = sc[23];
    Kc = refCinv * exp(g_RT[1] + g_RT[14] - g_RT[23]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[14]) + (h_RT[23]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[14] -= q; /* C2H2 */
    wdot[23] += q; /* C2H3 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[14];
        J[31] -= dqdci;               /* dwdot[H]/d[H] */
        J[44] -= dqdci;               /* dwdot[C2H2]/d[H] */
        J[53] += dqdci;               /* dwdot[C2H3]/d[H] */
        /* d()/d[H2] */
        dqdci = (TB[10][0] - 1)*dcdc_fac;
        J[61] -= dqdci;               /* dwdot[H]/d[H2] */
        J[74] -= dqdci;               /* dwdot[C2H2]/d[H2] */
        J[83] += dqdci;               /* dwdot[C2H3]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[10][1] - 1)*dcdc_fac;
        J[211] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[224] -= dqdci;              /* dwdot[C2H2]/d[H2O] */
        J[233] += dqdci;              /* dwdot[C2H3]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[10][3] - 1)*dcdc_fac;
        J[271] -= dqdci;              /* dwdot[H]/d[CO] */
        J[284] -= dqdci;              /* dwdot[C2H2]/d[CO] */
        J[293] += dqdci;              /* dwdot[C2H3]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[10][7] - 1)*dcdc_fac;
        J[361] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[374] -= dqdci;              /* dwdot[C2H2]/d[CO2] */
        J[383] += dqdci;              /* dwdot[C2H3]/d[CO2] */
        /* d()/d[CH4] */
        dqdci = (TB[10][9] - 1)*dcdc_fac;
        J[391] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[404] -= dqdci;              /* dwdot[C2H2]/d[CH4] */
        J[413] += dqdci;              /* dwdot[C2H3]/d[CH4] */
        /* d()/d[C2H2] */
        dqdci =  + k_f*sc[1];
        J[421] -= dqdci;              /* dwdot[H]/d[C2H2] */
        J[434] -= dqdci;              /* dwdot[C2H2]/d[C2H2] */
        J[443] += dqdci;              /* dwdot[C2H3]/d[C2H2] */
        /* d()/d[C2H6] */
        dqdci = (TB[10][16] - 1)*dcdc_fac;
        J[511] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[524] -= dqdci;              /* dwdot[C2H2]/d[C2H6] */
        J[533] += dqdci;              /* dwdot[C2H3]/d[C2H6] */
        /* d()/d[C] */
        dqdci = (TB[10][2] - 1)*dcdc_fac;
        J[541] -= dqdci;              /* dwdot[H]/d[C] */
        J[554] -= dqdci;              /* dwdot[C2H2]/d[C] */
        J[563] += dqdci;              /* dwdot[C2H3]/d[C] */
        /* d()/d[CH] */
        dqdci = (TB[10][4] - 1)*dcdc_fac;
        J[571] -= dqdci;              /* dwdot[H]/d[CH] */
        J[584] -= dqdci;              /* dwdot[C2H2]/d[CH] */
        J[593] += dqdci;              /* dwdot[C2H3]/d[CH] */
        /* d()/d[HCO] */
        dqdci = (TB[10][5] - 1)*dcdc_fac;
        J[601] -= dqdci;              /* dwdot[H]/d[HCO] */
        J[614] -= dqdci;              /* dwdot[C2H2]/d[HCO] */
        J[623] += dqdci;              /* dwdot[C2H3]/d[HCO] */
        /* d()/d[TXCH2] */
        dqdci = (TB[10][6] - 1)*dcdc_fac;
        J[631] -= dqdci;              /* dwdot[H]/d[TXCH2] */
        J[644] -= dqdci;              /* dwdot[C2H2]/d[TXCH2] */
        J[653] += dqdci;              /* dwdot[C2H3]/d[TXCH2] */
        /* d()/d[SXCH2] */
        dqdci = (TB[10][8] - 1)*dcdc_fac;
        J[661] -= dqdci;              /* dwdot[H]/d[SXCH2] */
        J[674] -= dqdci;              /* dwdot[C2H2]/d[SXCH2] */
        J[683] += dqdci;              /* dwdot[C2H3]/d[SXCH2] */
        /* d()/d[C2H3] */
        dqdci = (TB[10][10] - 1)*dcdc_fac - k_r;
        J[691] -= dqdci;              /* dwdot[H]/d[C2H3] */
        J[704] -= dqdci;              /* dwdot[C2H2]/d[C2H3] */
        J[713] += dqdci;              /* dwdot[C2H3]/d[C2H3] */
        /* d()/d[C2H5] */
        dqdci = (TB[10][11] - 1)*dcdc_fac;
        J[721] -= dqdci;              /* dwdot[H]/d[C2H5] */
        J[734] -= dqdci;              /* dwdot[C2H2]/d[C2H5] */
        J[743] += dqdci;              /* dwdot[C2H3]/d[C2H5] */
        /* d()/d[HCCO] */
        dqdci = (TB[10][12] - 1)*dcdc_fac;
        J[751] -= dqdci;              /* dwdot[H]/d[HCCO] */
        J[764] -= dqdci;              /* dwdot[C2H2]/d[HCCO] */
        J[773] += dqdci;              /* dwdot[C2H3]/d[HCCO] */
        /* d()/d[CH3CHO] */
        dqdci = (TB[10][13] - 1)*dcdc_fac;
        J[781] -= dqdci;              /* dwdot[H]/d[CH3CHO] */
        J[794] -= dqdci;              /* dwdot[C2H2]/d[CH3CHO] */
        J[803] += dqdci;              /* dwdot[C2H3]/d[CH3CHO] */
        /* d()/d[CH2CHO] */
        dqdci = (TB[10][14] - 1)*dcdc_fac;
        J[811] -= dqdci;              /* dwdot[H]/d[CH2CHO] */
        J[824] -= dqdci;              /* dwdot[C2H2]/d[CH2CHO] */
        J[833] += dqdci;              /* dwdot[C2H3]/d[CH2CHO] */
        /* d()/d[C2H5O] */
        dqdci = (TB[10][15] - 1)*dcdc_fac;
        J[841] -= dqdci;              /* dwdot[H]/d[C2H5O] */
        J[854] -= dqdci;              /* dwdot[C2H2]/d[C2H5O] */
        J[863] += dqdci;              /* dwdot[C2H3]/d[C2H5O] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[14];
        dqdc[2] = TB[10][0]*dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = TB[10][1]*dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = TB[10][3]*dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = TB[10][7]*dcdc_fac;
        dqdc[13] = TB[10][9]*dcdc_fac;
        dqdc[14] = dcdc_fac + k_f*sc[1];
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = TB[10][16]*dcdc_fac;
        dqdc[18] = TB[10][2]*dcdc_fac;
        dqdc[19] = TB[10][4]*dcdc_fac;
        dqdc[20] = TB[10][5]*dcdc_fac;
        dqdc[21] = TB[10][6]*dcdc_fac;
        dqdc[22] = TB[10][8]*dcdc_fac;
        dqdc[23] = TB[10][10]*dcdc_fac - k_r;
        dqdc[24] = TB[10][11]*dcdc_fac;
        dqdc[25] = TB[10][12]*dcdc_fac;
        dqdc[26] = TB[10][13]*dcdc_fac;
        dqdc[27] = TB[10][14]*dcdc_fac;
        dqdc[28] = TB[10][15]*dcdc_fac;
        for (int k=0; k<29; k++) {
            J[30*k+1] -= dqdc[k];
            J[30*k+14] -= dqdc[k];
            J[30*k+23] += dqdc[k];
        }
    }
    J[871] -= dqdT; /* dwdot[H]/dT */
    J[884] -= dqdT; /* dwdot[C2H2]/dT */
    J[893] += dqdT; /* dwdot[C2H3]/dT */

    /*reaction 12: C2H3 + H (+M) <=> C2H4 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[11][0] - 1)*sc[2] + (TB[11][1] - 1)*sc[7] + (TB[11][2] - 1)*sc[18] + (TB[11][3] - 1)*sc[9] + (TB[11][4] - 1)*sc[19] + (TB[11][5] - 1)*sc[20] + (TB[11][6] - 1)*sc[21] + (TB[11][7] - 1)*sc[12] + (TB[11][8] - 1)*sc[22] + (TB[11][9] - 1)*sc[13] + (TB[11][10] - 1)*sc[23] + (TB[11][11] - 1)*sc[24] + (TB[11][12] - 1)*sc[25] + (TB[11][13] - 1)*sc[26] + (TB[11][14] - 1)*sc[27] + (TB[11][15] - 1)*sc[28] + (TB[11][16] - 1)*sc[17];
    /* forward */
    phi_f = sc[1]*sc[23];
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
    phi_r = sc[15];
    Kc = refCinv * exp(g_RT[1] - g_RT[15] + g_RT[23]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[23]) + (h_RT[15]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[15] += q; /* C2H4 */
    wdot[23] -= q; /* C2H3 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[23];
        J[31] -= dqdci;               /* dwdot[H]/d[H] */
        J[45] += dqdci;               /* dwdot[C2H4]/d[H] */
        J[53] -= dqdci;               /* dwdot[C2H3]/d[H] */
        /* d()/d[H2] */
        dqdci = (TB[11][0] - 1)*dcdc_fac;
        J[61] -= dqdci;               /* dwdot[H]/d[H2] */
        J[75] += dqdci;               /* dwdot[C2H4]/d[H2] */
        J[83] -= dqdci;               /* dwdot[C2H3]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[11][1] - 1)*dcdc_fac;
        J[211] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[225] += dqdci;              /* dwdot[C2H4]/d[H2O] */
        J[233] -= dqdci;              /* dwdot[C2H3]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[11][3] - 1)*dcdc_fac;
        J[271] -= dqdci;              /* dwdot[H]/d[CO] */
        J[285] += dqdci;              /* dwdot[C2H4]/d[CO] */
        J[293] -= dqdci;              /* dwdot[C2H3]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[11][7] - 1)*dcdc_fac;
        J[361] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[375] += dqdci;              /* dwdot[C2H4]/d[CO2] */
        J[383] -= dqdci;              /* dwdot[C2H3]/d[CO2] */
        /* d()/d[CH4] */
        dqdci = (TB[11][9] - 1)*dcdc_fac;
        J[391] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[405] += dqdci;              /* dwdot[C2H4]/d[CH4] */
        J[413] -= dqdci;              /* dwdot[C2H3]/d[CH4] */
        /* d()/d[C2H4] */
        dqdci =  - k_r;
        J[451] -= dqdci;              /* dwdot[H]/d[C2H4] */
        J[465] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
        J[473] -= dqdci;              /* dwdot[C2H3]/d[C2H4] */
        /* d()/d[C2H6] */
        dqdci = (TB[11][16] - 1)*dcdc_fac;
        J[511] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[525] += dqdci;              /* dwdot[C2H4]/d[C2H6] */
        J[533] -= dqdci;              /* dwdot[C2H3]/d[C2H6] */
        /* d()/d[C] */
        dqdci = (TB[11][2] - 1)*dcdc_fac;
        J[541] -= dqdci;              /* dwdot[H]/d[C] */
        J[555] += dqdci;              /* dwdot[C2H4]/d[C] */
        J[563] -= dqdci;              /* dwdot[C2H3]/d[C] */
        /* d()/d[CH] */
        dqdci = (TB[11][4] - 1)*dcdc_fac;
        J[571] -= dqdci;              /* dwdot[H]/d[CH] */
        J[585] += dqdci;              /* dwdot[C2H4]/d[CH] */
        J[593] -= dqdci;              /* dwdot[C2H3]/d[CH] */
        /* d()/d[HCO] */
        dqdci = (TB[11][5] - 1)*dcdc_fac;
        J[601] -= dqdci;              /* dwdot[H]/d[HCO] */
        J[615] += dqdci;              /* dwdot[C2H4]/d[HCO] */
        J[623] -= dqdci;              /* dwdot[C2H3]/d[HCO] */
        /* d()/d[TXCH2] */
        dqdci = (TB[11][6] - 1)*dcdc_fac;
        J[631] -= dqdci;              /* dwdot[H]/d[TXCH2] */
        J[645] += dqdci;              /* dwdot[C2H4]/d[TXCH2] */
        J[653] -= dqdci;              /* dwdot[C2H3]/d[TXCH2] */
        /* d()/d[SXCH2] */
        dqdci = (TB[11][8] - 1)*dcdc_fac;
        J[661] -= dqdci;              /* dwdot[H]/d[SXCH2] */
        J[675] += dqdci;              /* dwdot[C2H4]/d[SXCH2] */
        J[683] -= dqdci;              /* dwdot[C2H3]/d[SXCH2] */
        /* d()/d[C2H3] */
        dqdci = (TB[11][10] - 1)*dcdc_fac + k_f*sc[1];
        J[691] -= dqdci;              /* dwdot[H]/d[C2H3] */
        J[705] += dqdci;              /* dwdot[C2H4]/d[C2H3] */
        J[713] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
        /* d()/d[C2H5] */
        dqdci = (TB[11][11] - 1)*dcdc_fac;
        J[721] -= dqdci;              /* dwdot[H]/d[C2H5] */
        J[735] += dqdci;              /* dwdot[C2H4]/d[C2H5] */
        J[743] -= dqdci;              /* dwdot[C2H3]/d[C2H5] */
        /* d()/d[HCCO] */
        dqdci = (TB[11][12] - 1)*dcdc_fac;
        J[751] -= dqdci;              /* dwdot[H]/d[HCCO] */
        J[765] += dqdci;              /* dwdot[C2H4]/d[HCCO] */
        J[773] -= dqdci;              /* dwdot[C2H3]/d[HCCO] */
        /* d()/d[CH3CHO] */
        dqdci = (TB[11][13] - 1)*dcdc_fac;
        J[781] -= dqdci;              /* dwdot[H]/d[CH3CHO] */
        J[795] += dqdci;              /* dwdot[C2H4]/d[CH3CHO] */
        J[803] -= dqdci;              /* dwdot[C2H3]/d[CH3CHO] */
        /* d()/d[CH2CHO] */
        dqdci = (TB[11][14] - 1)*dcdc_fac;
        J[811] -= dqdci;              /* dwdot[H]/d[CH2CHO] */
        J[825] += dqdci;              /* dwdot[C2H4]/d[CH2CHO] */
        J[833] -= dqdci;              /* dwdot[C2H3]/d[CH2CHO] */
        /* d()/d[C2H5O] */
        dqdci = (TB[11][15] - 1)*dcdc_fac;
        J[841] -= dqdci;              /* dwdot[H]/d[C2H5O] */
        J[855] += dqdci;              /* dwdot[C2H4]/d[C2H5O] */
        J[863] -= dqdci;              /* dwdot[C2H3]/d[C2H5O] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[23];
        dqdc[2] = TB[11][0]*dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = TB[11][1]*dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = TB[11][3]*dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = TB[11][7]*dcdc_fac;
        dqdc[13] = TB[11][9]*dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac - k_r;
        dqdc[16] = dcdc_fac;
        dqdc[17] = TB[11][16]*dcdc_fac;
        dqdc[18] = TB[11][2]*dcdc_fac;
        dqdc[19] = TB[11][4]*dcdc_fac;
        dqdc[20] = TB[11][5]*dcdc_fac;
        dqdc[21] = TB[11][6]*dcdc_fac;
        dqdc[22] = TB[11][8]*dcdc_fac;
        dqdc[23] = TB[11][10]*dcdc_fac + k_f*sc[1];
        dqdc[24] = TB[11][11]*dcdc_fac;
        dqdc[25] = TB[11][12]*dcdc_fac;
        dqdc[26] = TB[11][13]*dcdc_fac;
        dqdc[27] = TB[11][14]*dcdc_fac;
        dqdc[28] = TB[11][15]*dcdc_fac;
        for (int k=0; k<29; k++) {
            J[30*k+1] -= dqdc[k];
            J[30*k+15] += dqdc[k];
            J[30*k+23] -= dqdc[k];
        }
    }
    J[871] -= dqdT; /* dwdot[H]/dT */
    J[885] += dqdT; /* dwdot[C2H4]/dT */
    J[893] -= dqdT; /* dwdot[C2H3]/dT */

    /*reaction 13: C2H4 + H (+M) <=> C2H5 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[12][0] - 1)*sc[2] + (TB[12][1] - 1)*sc[7] + (TB[12][2] - 1)*sc[18] + (TB[12][3] - 1)*sc[9] + (TB[12][4] - 1)*sc[19] + (TB[12][5] - 1)*sc[20] + (TB[12][6] - 1)*sc[21] + (TB[12][7] - 1)*sc[12] + (TB[12][8] - 1)*sc[22] + (TB[12][9] - 1)*sc[13] + (TB[12][10] - 1)*sc[23] + (TB[12][11] - 1)*sc[24] + (TB[12][12] - 1)*sc[25] + (TB[12][13] - 1)*sc[26] + (TB[12][14] - 1)*sc[27] + (TB[12][15] - 1)*sc[28] + (TB[12][16] - 1)*sc[17];
    /* forward */
    phi_f = sc[1]*sc[15];
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
    phi_r = sc[24];
    Kc = refCinv * exp(g_RT[1] + g_RT[15] - g_RT[24]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[15]) + (h_RT[24]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[15] -= q; /* C2H4 */
    wdot[24] += q; /* C2H5 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[15];
        J[31] -= dqdci;               /* dwdot[H]/d[H] */
        J[45] -= dqdci;               /* dwdot[C2H4]/d[H] */
        J[54] += dqdci;               /* dwdot[C2H5]/d[H] */
        /* d()/d[H2] */
        dqdci = (TB[12][0] - 1)*dcdc_fac;
        J[61] -= dqdci;               /* dwdot[H]/d[H2] */
        J[75] -= dqdci;               /* dwdot[C2H4]/d[H2] */
        J[84] += dqdci;               /* dwdot[C2H5]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[12][1] - 1)*dcdc_fac;
        J[211] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[225] -= dqdci;              /* dwdot[C2H4]/d[H2O] */
        J[234] += dqdci;              /* dwdot[C2H5]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[12][3] - 1)*dcdc_fac;
        J[271] -= dqdci;              /* dwdot[H]/d[CO] */
        J[285] -= dqdci;              /* dwdot[C2H4]/d[CO] */
        J[294] += dqdci;              /* dwdot[C2H5]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[12][7] - 1)*dcdc_fac;
        J[361] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[375] -= dqdci;              /* dwdot[C2H4]/d[CO2] */
        J[384] += dqdci;              /* dwdot[C2H5]/d[CO2] */
        /* d()/d[CH4] */
        dqdci = (TB[12][9] - 1)*dcdc_fac;
        J[391] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[405] -= dqdci;              /* dwdot[C2H4]/d[CH4] */
        J[414] += dqdci;              /* dwdot[C2H5]/d[CH4] */
        /* d()/d[C2H4] */
        dqdci =  + k_f*sc[1];
        J[451] -= dqdci;              /* dwdot[H]/d[C2H4] */
        J[465] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
        J[474] += dqdci;              /* dwdot[C2H5]/d[C2H4] */
        /* d()/d[C2H6] */
        dqdci = (TB[12][16] - 1)*dcdc_fac;
        J[511] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[525] -= dqdci;              /* dwdot[C2H4]/d[C2H6] */
        J[534] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
        /* d()/d[C] */
        dqdci = (TB[12][2] - 1)*dcdc_fac;
        J[541] -= dqdci;              /* dwdot[H]/d[C] */
        J[555] -= dqdci;              /* dwdot[C2H4]/d[C] */
        J[564] += dqdci;              /* dwdot[C2H5]/d[C] */
        /* d()/d[CH] */
        dqdci = (TB[12][4] - 1)*dcdc_fac;
        J[571] -= dqdci;              /* dwdot[H]/d[CH] */
        J[585] -= dqdci;              /* dwdot[C2H4]/d[CH] */
        J[594] += dqdci;              /* dwdot[C2H5]/d[CH] */
        /* d()/d[HCO] */
        dqdci = (TB[12][5] - 1)*dcdc_fac;
        J[601] -= dqdci;              /* dwdot[H]/d[HCO] */
        J[615] -= dqdci;              /* dwdot[C2H4]/d[HCO] */
        J[624] += dqdci;              /* dwdot[C2H5]/d[HCO] */
        /* d()/d[TXCH2] */
        dqdci = (TB[12][6] - 1)*dcdc_fac;
        J[631] -= dqdci;              /* dwdot[H]/d[TXCH2] */
        J[645] -= dqdci;              /* dwdot[C2H4]/d[TXCH2] */
        J[654] += dqdci;              /* dwdot[C2H5]/d[TXCH2] */
        /* d()/d[SXCH2] */
        dqdci = (TB[12][8] - 1)*dcdc_fac;
        J[661] -= dqdci;              /* dwdot[H]/d[SXCH2] */
        J[675] -= dqdci;              /* dwdot[C2H4]/d[SXCH2] */
        J[684] += dqdci;              /* dwdot[C2H5]/d[SXCH2] */
        /* d()/d[C2H3] */
        dqdci = (TB[12][10] - 1)*dcdc_fac;
        J[691] -= dqdci;              /* dwdot[H]/d[C2H3] */
        J[705] -= dqdci;              /* dwdot[C2H4]/d[C2H3] */
        J[714] += dqdci;              /* dwdot[C2H5]/d[C2H3] */
        /* d()/d[C2H5] */
        dqdci = (TB[12][11] - 1)*dcdc_fac - k_r;
        J[721] -= dqdci;              /* dwdot[H]/d[C2H5] */
        J[735] -= dqdci;              /* dwdot[C2H4]/d[C2H5] */
        J[744] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
        /* d()/d[HCCO] */
        dqdci = (TB[12][12] - 1)*dcdc_fac;
        J[751] -= dqdci;              /* dwdot[H]/d[HCCO] */
        J[765] -= dqdci;              /* dwdot[C2H4]/d[HCCO] */
        J[774] += dqdci;              /* dwdot[C2H5]/d[HCCO] */
        /* d()/d[CH3CHO] */
        dqdci = (TB[12][13] - 1)*dcdc_fac;
        J[781] -= dqdci;              /* dwdot[H]/d[CH3CHO] */
        J[795] -= dqdci;              /* dwdot[C2H4]/d[CH3CHO] */
        J[804] += dqdci;              /* dwdot[C2H5]/d[CH3CHO] */
        /* d()/d[CH2CHO] */
        dqdci = (TB[12][14] - 1)*dcdc_fac;
        J[811] -= dqdci;              /* dwdot[H]/d[CH2CHO] */
        J[825] -= dqdci;              /* dwdot[C2H4]/d[CH2CHO] */
        J[834] += dqdci;              /* dwdot[C2H5]/d[CH2CHO] */
        /* d()/d[C2H5O] */
        dqdci = (TB[12][15] - 1)*dcdc_fac;
        J[841] -= dqdci;              /* dwdot[H]/d[C2H5O] */
        J[855] -= dqdci;              /* dwdot[C2H4]/d[C2H5O] */
        J[864] += dqdci;              /* dwdot[C2H5]/d[C2H5O] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[15];
        dqdc[2] = TB[12][0]*dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = TB[12][1]*dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = TB[12][3]*dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = TB[12][7]*dcdc_fac;
        dqdc[13] = TB[12][9]*dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac + k_f*sc[1];
        dqdc[16] = dcdc_fac;
        dqdc[17] = TB[12][16]*dcdc_fac;
        dqdc[18] = TB[12][2]*dcdc_fac;
        dqdc[19] = TB[12][4]*dcdc_fac;
        dqdc[20] = TB[12][5]*dcdc_fac;
        dqdc[21] = TB[12][6]*dcdc_fac;
        dqdc[22] = TB[12][8]*dcdc_fac;
        dqdc[23] = TB[12][10]*dcdc_fac;
        dqdc[24] = TB[12][11]*dcdc_fac - k_r;
        dqdc[25] = TB[12][12]*dcdc_fac;
        dqdc[26] = TB[12][13]*dcdc_fac;
        dqdc[27] = TB[12][14]*dcdc_fac;
        dqdc[28] = TB[12][15]*dcdc_fac;
        for (int k=0; k<29; k++) {
            J[30*k+1] -= dqdc[k];
            J[30*k+15] -= dqdc[k];
            J[30*k+24] += dqdc[k];
        }
    }
    J[871] -= dqdT; /* dwdot[H]/dT */
    J[885] -= dqdT; /* dwdot[C2H4]/dT */
    J[894] += dqdT; /* dwdot[C2H5]/dT */

    /*reaction 14: C2H5 + H (+M) <=> C2H6 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[13][0] - 1)*sc[2] + (TB[13][1] - 1)*sc[7] + (TB[13][2] - 1)*sc[18] + (TB[13][3] - 1)*sc[9] + (TB[13][4] - 1)*sc[19] + (TB[13][5] - 1)*sc[20] + (TB[13][6] - 1)*sc[21] + (TB[13][7] - 1)*sc[12] + (TB[13][8] - 1)*sc[22] + (TB[13][9] - 1)*sc[13] + (TB[13][10] - 1)*sc[23] + (TB[13][11] - 1)*sc[24] + (TB[13][12] - 1)*sc[25] + (TB[13][13] - 1)*sc[26] + (TB[13][14] - 1)*sc[27] + (TB[13][15] - 1)*sc[28] + (TB[13][16] - 1)*sc[17];
    /* forward */
    phi_f = sc[1]*sc[24];
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
    phi_r = sc[17];
    Kc = refCinv * exp(g_RT[1] - g_RT[17] + g_RT[24]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[24]) + (h_RT[17]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[17] += q; /* C2H6 */
    wdot[24] -= q; /* C2H5 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[24];
        J[31] -= dqdci;               /* dwdot[H]/d[H] */
        J[47] += dqdci;               /* dwdot[C2H6]/d[H] */
        J[54] -= dqdci;               /* dwdot[C2H5]/d[H] */
        /* d()/d[H2] */
        dqdci = (TB[13][0] - 1)*dcdc_fac;
        J[61] -= dqdci;               /* dwdot[H]/d[H2] */
        J[77] += dqdci;               /* dwdot[C2H6]/d[H2] */
        J[84] -= dqdci;               /* dwdot[C2H5]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[13][1] - 1)*dcdc_fac;
        J[211] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[227] += dqdci;              /* dwdot[C2H6]/d[H2O] */
        J[234] -= dqdci;              /* dwdot[C2H5]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[13][3] - 1)*dcdc_fac;
        J[271] -= dqdci;              /* dwdot[H]/d[CO] */
        J[287] += dqdci;              /* dwdot[C2H6]/d[CO] */
        J[294] -= dqdci;              /* dwdot[C2H5]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[13][7] - 1)*dcdc_fac;
        J[361] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[377] += dqdci;              /* dwdot[C2H6]/d[CO2] */
        J[384] -= dqdci;              /* dwdot[C2H5]/d[CO2] */
        /* d()/d[CH4] */
        dqdci = (TB[13][9] - 1)*dcdc_fac;
        J[391] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[407] += dqdci;              /* dwdot[C2H6]/d[CH4] */
        J[414] -= dqdci;              /* dwdot[C2H5]/d[CH4] */
        /* d()/d[C2H6] */
        dqdci = (TB[13][16] - 1)*dcdc_fac - k_r;
        J[511] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[527] += dqdci;              /* dwdot[C2H6]/d[C2H6] */
        J[534] -= dqdci;              /* dwdot[C2H5]/d[C2H6] */
        /* d()/d[C] */
        dqdci = (TB[13][2] - 1)*dcdc_fac;
        J[541] -= dqdci;              /* dwdot[H]/d[C] */
        J[557] += dqdci;              /* dwdot[C2H6]/d[C] */
        J[564] -= dqdci;              /* dwdot[C2H5]/d[C] */
        /* d()/d[CH] */
        dqdci = (TB[13][4] - 1)*dcdc_fac;
        J[571] -= dqdci;              /* dwdot[H]/d[CH] */
        J[587] += dqdci;              /* dwdot[C2H6]/d[CH] */
        J[594] -= dqdci;              /* dwdot[C2H5]/d[CH] */
        /* d()/d[HCO] */
        dqdci = (TB[13][5] - 1)*dcdc_fac;
        J[601] -= dqdci;              /* dwdot[H]/d[HCO] */
        J[617] += dqdci;              /* dwdot[C2H6]/d[HCO] */
        J[624] -= dqdci;              /* dwdot[C2H5]/d[HCO] */
        /* d()/d[TXCH2] */
        dqdci = (TB[13][6] - 1)*dcdc_fac;
        J[631] -= dqdci;              /* dwdot[H]/d[TXCH2] */
        J[647] += dqdci;              /* dwdot[C2H6]/d[TXCH2] */
        J[654] -= dqdci;              /* dwdot[C2H5]/d[TXCH2] */
        /* d()/d[SXCH2] */
        dqdci = (TB[13][8] - 1)*dcdc_fac;
        J[661] -= dqdci;              /* dwdot[H]/d[SXCH2] */
        J[677] += dqdci;              /* dwdot[C2H6]/d[SXCH2] */
        J[684] -= dqdci;              /* dwdot[C2H5]/d[SXCH2] */
        /* d()/d[C2H3] */
        dqdci = (TB[13][10] - 1)*dcdc_fac;
        J[691] -= dqdci;              /* dwdot[H]/d[C2H3] */
        J[707] += dqdci;              /* dwdot[C2H6]/d[C2H3] */
        J[714] -= dqdci;              /* dwdot[C2H5]/d[C2H3] */
        /* d()/d[C2H5] */
        dqdci = (TB[13][11] - 1)*dcdc_fac + k_f*sc[1];
        J[721] -= dqdci;              /* dwdot[H]/d[C2H5] */
        J[737] += dqdci;              /* dwdot[C2H6]/d[C2H5] */
        J[744] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
        /* d()/d[HCCO] */
        dqdci = (TB[13][12] - 1)*dcdc_fac;
        J[751] -= dqdci;              /* dwdot[H]/d[HCCO] */
        J[767] += dqdci;              /* dwdot[C2H6]/d[HCCO] */
        J[774] -= dqdci;              /* dwdot[C2H5]/d[HCCO] */
        /* d()/d[CH3CHO] */
        dqdci = (TB[13][13] - 1)*dcdc_fac;
        J[781] -= dqdci;              /* dwdot[H]/d[CH3CHO] */
        J[797] += dqdci;              /* dwdot[C2H6]/d[CH3CHO] */
        J[804] -= dqdci;              /* dwdot[C2H5]/d[CH3CHO] */
        /* d()/d[CH2CHO] */
        dqdci = (TB[13][14] - 1)*dcdc_fac;
        J[811] -= dqdci;              /* dwdot[H]/d[CH2CHO] */
        J[827] += dqdci;              /* dwdot[C2H6]/d[CH2CHO] */
        J[834] -= dqdci;              /* dwdot[C2H5]/d[CH2CHO] */
        /* d()/d[C2H5O] */
        dqdci = (TB[13][15] - 1)*dcdc_fac;
        J[841] -= dqdci;              /* dwdot[H]/d[C2H5O] */
        J[857] += dqdci;              /* dwdot[C2H6]/d[C2H5O] */
        J[864] -= dqdci;              /* dwdot[C2H5]/d[C2H5O] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[24];
        dqdc[2] = TB[13][0]*dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = TB[13][1]*dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = TB[13][3]*dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = TB[13][7]*dcdc_fac;
        dqdc[13] = TB[13][9]*dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = TB[13][16]*dcdc_fac - k_r;
        dqdc[18] = TB[13][2]*dcdc_fac;
        dqdc[19] = TB[13][4]*dcdc_fac;
        dqdc[20] = TB[13][5]*dcdc_fac;
        dqdc[21] = TB[13][6]*dcdc_fac;
        dqdc[22] = TB[13][8]*dcdc_fac;
        dqdc[23] = TB[13][10]*dcdc_fac;
        dqdc[24] = TB[13][11]*dcdc_fac + k_f*sc[1];
        dqdc[25] = TB[13][12]*dcdc_fac;
        dqdc[26] = TB[13][13]*dcdc_fac;
        dqdc[27] = TB[13][14]*dcdc_fac;
        dqdc[28] = TB[13][15]*dcdc_fac;
        for (int k=0; k<29; k++) {
            J[30*k+1] -= dqdc[k];
            J[30*k+17] += dqdc[k];
            J[30*k+24] -= dqdc[k];
        }
    }
    J[871] -= dqdT; /* dwdot[H]/dT */
    J[887] += dqdT; /* dwdot[C2H6]/dT */
    J[894] -= dqdT; /* dwdot[C2H5]/dT */

    /*reaction 15: C2H6 (+M) <=> 2.000000 CH3 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[14][0] - 1)*sc[2] + (TB[14][1] - 1)*sc[7] + (TB[14][2] - 1)*sc[18] + (TB[14][3] - 1)*sc[9] + (TB[14][4] - 1)*sc[19] + (TB[14][5] - 1)*sc[20] + (TB[14][6] - 1)*sc[21] + (TB[14][7] - 1)*sc[12] + (TB[14][8] - 1)*sc[22] + (TB[14][9] - 1)*sc[13] + (TB[14][10] - 1)*sc[23] + (TB[14][11] - 1)*sc[24] + (TB[14][12] - 1)*sc[25] + (TB[14][13] - 1)*sc[26] + (TB[14][14] - 1)*sc[27] + (TB[14][15] - 1)*sc[28] + (TB[14][16] - 1)*sc[17];
    /* forward */
    phi_f = sc[17];
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
    phi_r = pow(sc[10], 2.000000);
    Kc = refC * exp(-2.000000*g_RT[10] + g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[17]) + (2.000000*h_RT[10]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[10] += 2 * q; /* CH3 */
    wdot[17] -= q; /* C2H6 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[14][0] - 1)*dcdc_fac;
        J[70] += 2 * dqdci;           /* dwdot[CH3]/d[H2] */
        J[77] -= dqdci;               /* dwdot[C2H6]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[14][1] - 1)*dcdc_fac;
        J[220] += 2 * dqdci;          /* dwdot[CH3]/d[H2O] */
        J[227] -= dqdci;              /* dwdot[C2H6]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[14][3] - 1)*dcdc_fac;
        J[280] += 2 * dqdci;          /* dwdot[CH3]/d[CO] */
        J[287] -= dqdci;              /* dwdot[C2H6]/d[CO] */
        /* d()/d[CH3] */
        dqdci =  - k_r*2.000000*sc[10];
        J[310] += 2 * dqdci;          /* dwdot[CH3]/d[CH3] */
        J[317] -= dqdci;              /* dwdot[C2H6]/d[CH3] */
        /* d()/d[CO2] */
        dqdci = (TB[14][7] - 1)*dcdc_fac;
        J[370] += 2 * dqdci;          /* dwdot[CH3]/d[CO2] */
        J[377] -= dqdci;              /* dwdot[C2H6]/d[CO2] */
        /* d()/d[CH4] */
        dqdci = (TB[14][9] - 1)*dcdc_fac;
        J[400] += 2 * dqdci;          /* dwdot[CH3]/d[CH4] */
        J[407] -= dqdci;              /* dwdot[C2H6]/d[CH4] */
        /* d()/d[C2H6] */
        dqdci = (TB[14][16] - 1)*dcdc_fac + k_f;
        J[520] += 2 * dqdci;          /* dwdot[CH3]/d[C2H6] */
        J[527] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
        /* d()/d[C] */
        dqdci = (TB[14][2] - 1)*dcdc_fac;
        J[550] += 2 * dqdci;          /* dwdot[CH3]/d[C] */
        J[557] -= dqdci;              /* dwdot[C2H6]/d[C] */
        /* d()/d[CH] */
        dqdci = (TB[14][4] - 1)*dcdc_fac;
        J[580] += 2 * dqdci;          /* dwdot[CH3]/d[CH] */
        J[587] -= dqdci;              /* dwdot[C2H6]/d[CH] */
        /* d()/d[HCO] */
        dqdci = (TB[14][5] - 1)*dcdc_fac;
        J[610] += 2 * dqdci;          /* dwdot[CH3]/d[HCO] */
        J[617] -= dqdci;              /* dwdot[C2H6]/d[HCO] */
        /* d()/d[TXCH2] */
        dqdci = (TB[14][6] - 1)*dcdc_fac;
        J[640] += 2 * dqdci;          /* dwdot[CH3]/d[TXCH2] */
        J[647] -= dqdci;              /* dwdot[C2H6]/d[TXCH2] */
        /* d()/d[SXCH2] */
        dqdci = (TB[14][8] - 1)*dcdc_fac;
        J[670] += 2 * dqdci;          /* dwdot[CH3]/d[SXCH2] */
        J[677] -= dqdci;              /* dwdot[C2H6]/d[SXCH2] */
        /* d()/d[C2H3] */
        dqdci = (TB[14][10] - 1)*dcdc_fac;
        J[700] += 2 * dqdci;          /* dwdot[CH3]/d[C2H3] */
        J[707] -= dqdci;              /* dwdot[C2H6]/d[C2H3] */
        /* d()/d[C2H5] */
        dqdci = (TB[14][11] - 1)*dcdc_fac;
        J[730] += 2 * dqdci;          /* dwdot[CH3]/d[C2H5] */
        J[737] -= dqdci;              /* dwdot[C2H6]/d[C2H5] */
        /* d()/d[HCCO] */
        dqdci = (TB[14][12] - 1)*dcdc_fac;
        J[760] += 2 * dqdci;          /* dwdot[CH3]/d[HCCO] */
        J[767] -= dqdci;              /* dwdot[C2H6]/d[HCCO] */
        /* d()/d[CH3CHO] */
        dqdci = (TB[14][13] - 1)*dcdc_fac;
        J[790] += 2 * dqdci;          /* dwdot[CH3]/d[CH3CHO] */
        J[797] -= dqdci;              /* dwdot[C2H6]/d[CH3CHO] */
        /* d()/d[CH2CHO] */
        dqdci = (TB[14][14] - 1)*dcdc_fac;
        J[820] += 2 * dqdci;          /* dwdot[CH3]/d[CH2CHO] */
        J[827] -= dqdci;              /* dwdot[C2H6]/d[CH2CHO] */
        /* d()/d[C2H5O] */
        dqdci = (TB[14][15] - 1)*dcdc_fac;
        J[850] += 2 * dqdci;          /* dwdot[CH3]/d[C2H5O] */
        J[857] -= dqdci;              /* dwdot[C2H6]/d[C2H5O] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = TB[14][0]*dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = TB[14][1]*dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = TB[14][3]*dcdc_fac;
        dqdc[10] = dcdc_fac - k_r*2.000000*sc[10];
        dqdc[11] = dcdc_fac;
        dqdc[12] = TB[14][7]*dcdc_fac;
        dqdc[13] = TB[14][9]*dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = TB[14][16]*dcdc_fac + k_f;
        dqdc[18] = TB[14][2]*dcdc_fac;
        dqdc[19] = TB[14][4]*dcdc_fac;
        dqdc[20] = TB[14][5]*dcdc_fac;
        dqdc[21] = TB[14][6]*dcdc_fac;
        dqdc[22] = TB[14][8]*dcdc_fac;
        dqdc[23] = TB[14][10]*dcdc_fac;
        dqdc[24] = TB[14][11]*dcdc_fac;
        dqdc[25] = TB[14][12]*dcdc_fac;
        dqdc[26] = TB[14][13]*dcdc_fac;
        dqdc[27] = TB[14][14]*dcdc_fac;
        dqdc[28] = TB[14][15]*dcdc_fac;
        for (int k=0; k<29; k++) {
            J[30*k+10] += 2 * dqdc[k];
            J[30*k+17] -= dqdc[k];
        }
    }
    J[880] += 2 * dqdT; /* dwdot[CH3]/dT */
    J[887] -= dqdT; /* dwdot[C2H6]/dT */

    /*reaction 16: 2.000000 H + M <=> H2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[15][0] - 1)*sc[2] + (TB[15][1] - 1)*sc[7] + (TB[15][2] - 1)*sc[18] + (TB[15][3] - 1)*sc[19] + (TB[15][4] - 1)*sc[20] + (TB[15][5] - 1)*sc[21] + (TB[15][6] - 1)*sc[12] + (TB[15][7] - 1)*sc[22] + (TB[15][8] - 1)*sc[13] + (TB[15][9] - 1)*sc[23] + (TB[15][10] - 1)*sc[24] + (TB[15][11] - 1)*sc[25] + (TB[15][12] - 1)*sc[26] + (TB[15][13] - 1)*sc[27] + (TB[15][14] - 1)*sc[28] + (TB[15][15] - 1)*sc[17];
    /* forward */
    phi_f = pow(sc[1], 2.000000);
    k_f = prefactor_units[15] * fwd_A[15]
                * exp(fwd_beta[15] * tc[0] - activation_units[15] * fwd_Ea[15] * invT);
    dlnkfdT = fwd_beta[15] * invT + activation_units[15] * fwd_Ea[15] * invT2;
    /* reverse */
    phi_r = sc[2];
    Kc = refCinv * exp(2.000000*g_RT[1] - g_RT[2]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[1]) + (h_RT[2]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= 2 * q; /* H */
    wdot[2] += q; /* H2 */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*2.000000*sc[1];
        J[31] += -2 * dqdci;          /* dwdot[H]/d[H] */
        J[32] += dqdci;               /* dwdot[H2]/d[H] */
        /* d()/d[H2] */
        dqdci = (TB[15][0] - 1)*q_nocor - k_r;
        J[61] += -2 * dqdci;          /* dwdot[H]/d[H2] */
        J[62] += dqdci;               /* dwdot[H2]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[15][1] - 1)*q_nocor;
        J[211] += -2 * dqdci;         /* dwdot[H]/d[H2O] */
        J[212] += dqdci;              /* dwdot[H2]/d[H2O] */
        /* d()/d[CO2] */
        dqdci = (TB[15][6] - 1)*q_nocor;
        J[361] += -2 * dqdci;         /* dwdot[H]/d[CO2] */
        J[362] += dqdci;              /* dwdot[H2]/d[CO2] */
        /* d()/d[CH4] */
        dqdci = (TB[15][8] - 1)*q_nocor;
        J[391] += -2 * dqdci;         /* dwdot[H]/d[CH4] */
        J[392] += dqdci;              /* dwdot[H2]/d[CH4] */
        /* d()/d[C2H6] */
        dqdci = (TB[15][15] - 1)*q_nocor;
        J[511] += -2 * dqdci;         /* dwdot[H]/d[C2H6] */
        J[512] += dqdci;              /* dwdot[H2]/d[C2H6] */
        /* d()/d[C] */
        dqdci = (TB[15][2] - 1)*q_nocor;
        J[541] += -2 * dqdci;         /* dwdot[H]/d[C] */
        J[542] += dqdci;              /* dwdot[H2]/d[C] */
        /* d()/d[CH] */
        dqdci = (TB[15][3] - 1)*q_nocor;
        J[571] += -2 * dqdci;         /* dwdot[H]/d[CH] */
        J[572] += dqdci;              /* dwdot[H2]/d[CH] */
        /* d()/d[HCO] */
        dqdci = (TB[15][4] - 1)*q_nocor;
        J[601] += -2 * dqdci;         /* dwdot[H]/d[HCO] */
        J[602] += dqdci;              /* dwdot[H2]/d[HCO] */
        /* d()/d[TXCH2] */
        dqdci = (TB[15][5] - 1)*q_nocor;
        J[631] += -2 * dqdci;         /* dwdot[H]/d[TXCH2] */
        J[632] += dqdci;              /* dwdot[H2]/d[TXCH2] */
        /* d()/d[SXCH2] */
        dqdci = (TB[15][7] - 1)*q_nocor;
        J[661] += -2 * dqdci;         /* dwdot[H]/d[SXCH2] */
        J[662] += dqdci;              /* dwdot[H2]/d[SXCH2] */
        /* d()/d[C2H3] */
        dqdci = (TB[15][9] - 1)*q_nocor;
        J[691] += -2 * dqdci;         /* dwdot[H]/d[C2H3] */
        J[692] += dqdci;              /* dwdot[H2]/d[C2H3] */
        /* d()/d[C2H5] */
        dqdci = (TB[15][10] - 1)*q_nocor;
        J[721] += -2 * dqdci;         /* dwdot[H]/d[C2H5] */
        J[722] += dqdci;              /* dwdot[H2]/d[C2H5] */
        /* d()/d[HCCO] */
        dqdci = (TB[15][11] - 1)*q_nocor;
        J[751] += -2 * dqdci;         /* dwdot[H]/d[HCCO] */
        J[752] += dqdci;              /* dwdot[H2]/d[HCCO] */
        /* d()/d[CH3CHO] */
        dqdci = (TB[15][12] - 1)*q_nocor;
        J[781] += -2 * dqdci;         /* dwdot[H]/d[CH3CHO] */
        J[782] += dqdci;              /* dwdot[H2]/d[CH3CHO] */
        /* d()/d[CH2CHO] */
        dqdci = (TB[15][13] - 1)*q_nocor;
        J[811] += -2 * dqdci;         /* dwdot[H]/d[CH2CHO] */
        J[812] += dqdci;              /* dwdot[H2]/d[CH2CHO] */
        /* d()/d[C2H5O] */
        dqdci = (TB[15][14] - 1)*q_nocor;
        J[841] += -2 * dqdci;         /* dwdot[H]/d[C2H5O] */
        J[842] += dqdci;              /* dwdot[H2]/d[C2H5O] */
    }
    else {
        dqdc[0] = q_nocor;
        dqdc[1] = q_nocor + k_f*2.000000*sc[1];
        dqdc[2] = TB[15][0]*q_nocor - k_r;
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor;
        dqdc[5] = q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = TB[15][1]*q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = q_nocor;
        dqdc[12] = TB[15][6]*q_nocor;
        dqdc[13] = TB[15][8]*q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = TB[15][15]*q_nocor;
        dqdc[18] = TB[15][2]*q_nocor;
        dqdc[19] = TB[15][3]*q_nocor;
        dqdc[20] = TB[15][4]*q_nocor;
        dqdc[21] = TB[15][5]*q_nocor;
        dqdc[22] = TB[15][7]*q_nocor;
        dqdc[23] = TB[15][9]*q_nocor;
        dqdc[24] = TB[15][10]*q_nocor;
        dqdc[25] = TB[15][11]*q_nocor;
        dqdc[26] = TB[15][12]*q_nocor;
        dqdc[27] = TB[15][13]*q_nocor;
        dqdc[28] = TB[15][14]*q_nocor;
        for (int k=0; k<29; k++) {
            J[30*k+1] += -2 * dqdc[k];
            J[30*k+2] += dqdc[k];
        }
    }
    J[871] += -2 * dqdT; /* dwdot[H]/dT */
    J[872] += dqdT; /* dwdot[H2]/dT */

    /*reaction 17: H + O + M <=> OH + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[16][0] - 1)*sc[2] + (TB[16][1] - 1)*sc[7] + (TB[16][2] - 1)*sc[18] + (TB[16][3] - 1)*sc[9] + (TB[16][4] - 1)*sc[19] + (TB[16][5] - 1)*sc[20] + (TB[16][6] - 1)*sc[21] + (TB[16][7] - 1)*sc[12] + (TB[16][8] - 1)*sc[22] + (TB[16][9] - 1)*sc[13] + (TB[16][10] - 1)*sc[23] + (TB[16][11] - 1)*sc[24] + (TB[16][12] - 1)*sc[25] + (TB[16][13] - 1)*sc[26] + (TB[16][14] - 1)*sc[27] + (TB[16][15] - 1)*sc[28] + (TB[16][16] - 1)*sc[17];
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = prefactor_units[16] * fwd_A[16]
                * exp(fwd_beta[16] * tc[0] - activation_units[16] * fwd_Ea[16] * invT);
    dlnkfdT = fwd_beta[16] * invT + activation_units[16] * fwd_Ea[16] * invT2;
    /* reverse */
    phi_r = sc[4];
    Kc = refCinv * exp(g_RT[1] + g_RT[3] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3]) + (h_RT[4]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[3] -= q; /* O */
    wdot[4] += q; /* OH */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[3];
        J[31] -= dqdci;               /* dwdot[H]/d[H] */
        J[33] -= dqdci;               /* dwdot[O]/d[H] */
        J[34] += dqdci;               /* dwdot[OH]/d[H] */
        /* d()/d[H2] */
        dqdci = (TB[16][0] - 1)*q_nocor;
        J[61] -= dqdci;               /* dwdot[H]/d[H2] */
        J[63] -= dqdci;               /* dwdot[O]/d[H2] */
        J[64] += dqdci;               /* dwdot[OH]/d[H2] */
        /* d()/d[O] */
        dqdci =  + k_f*sc[1];
        J[91] -= dqdci;               /* dwdot[H]/d[O] */
        J[93] -= dqdci;               /* dwdot[O]/d[O] */
        J[94] += dqdci;               /* dwdot[OH]/d[O] */
        /* d()/d[OH] */
        dqdci =  - k_r;
        J[121] -= dqdci;              /* dwdot[H]/d[OH] */
        J[123] -= dqdci;              /* dwdot[O]/d[OH] */
        J[124] += dqdci;              /* dwdot[OH]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (TB[16][1] - 1)*q_nocor;
        J[211] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[213] -= dqdci;              /* dwdot[O]/d[H2O] */
        J[214] += dqdci;              /* dwdot[OH]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[16][3] - 1)*q_nocor;
        J[271] -= dqdci;              /* dwdot[H]/d[CO] */
        J[273] -= dqdci;              /* dwdot[O]/d[CO] */
        J[274] += dqdci;              /* dwdot[OH]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[16][7] - 1)*q_nocor;
        J[361] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[363] -= dqdci;              /* dwdot[O]/d[CO2] */
        J[364] += dqdci;              /* dwdot[OH]/d[CO2] */
        /* d()/d[CH4] */
        dqdci = (TB[16][9] - 1)*q_nocor;
        J[391] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[393] -= dqdci;              /* dwdot[O]/d[CH4] */
        J[394] += dqdci;              /* dwdot[OH]/d[CH4] */
        /* d()/d[C2H6] */
        dqdci = (TB[16][16] - 1)*q_nocor;
        J[511] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[513] -= dqdci;              /* dwdot[O]/d[C2H6] */
        J[514] += dqdci;              /* dwdot[OH]/d[C2H6] */
        /* d()/d[C] */
        dqdci = (TB[16][2] - 1)*q_nocor;
        J[541] -= dqdci;              /* dwdot[H]/d[C] */
        J[543] -= dqdci;              /* dwdot[O]/d[C] */
        J[544] += dqdci;              /* dwdot[OH]/d[C] */
        /* d()/d[CH] */
        dqdci = (TB[16][4] - 1)*q_nocor;
        J[571] -= dqdci;              /* dwdot[H]/d[CH] */
        J[573] -= dqdci;              /* dwdot[O]/d[CH] */
        J[574] += dqdci;              /* dwdot[OH]/d[CH] */
        /* d()/d[HCO] */
        dqdci = (TB[16][5] - 1)*q_nocor;
        J[601] -= dqdci;              /* dwdot[H]/d[HCO] */
        J[603] -= dqdci;              /* dwdot[O]/d[HCO] */
        J[604] += dqdci;              /* dwdot[OH]/d[HCO] */
        /* d()/d[TXCH2] */
        dqdci = (TB[16][6] - 1)*q_nocor;
        J[631] -= dqdci;              /* dwdot[H]/d[TXCH2] */
        J[633] -= dqdci;              /* dwdot[O]/d[TXCH2] */
        J[634] += dqdci;              /* dwdot[OH]/d[TXCH2] */
        /* d()/d[SXCH2] */
        dqdci = (TB[16][8] - 1)*q_nocor;
        J[661] -= dqdci;              /* dwdot[H]/d[SXCH2] */
        J[663] -= dqdci;              /* dwdot[O]/d[SXCH2] */
        J[664] += dqdci;              /* dwdot[OH]/d[SXCH2] */
        /* d()/d[C2H3] */
        dqdci = (TB[16][10] - 1)*q_nocor;
        J[691] -= dqdci;              /* dwdot[H]/d[C2H3] */
        J[693] -= dqdci;              /* dwdot[O]/d[C2H3] */
        J[694] += dqdci;              /* dwdot[OH]/d[C2H3] */
        /* d()/d[C2H5] */
        dqdci = (TB[16][11] - 1)*q_nocor;
        J[721] -= dqdci;              /* dwdot[H]/d[C2H5] */
        J[723] -= dqdci;              /* dwdot[O]/d[C2H5] */
        J[724] += dqdci;              /* dwdot[OH]/d[C2H5] */
        /* d()/d[HCCO] */
        dqdci = (TB[16][12] - 1)*q_nocor;
        J[751] -= dqdci;              /* dwdot[H]/d[HCCO] */
        J[753] -= dqdci;              /* dwdot[O]/d[HCCO] */
        J[754] += dqdci;              /* dwdot[OH]/d[HCCO] */
        /* d()/d[CH3CHO] */
        dqdci = (TB[16][13] - 1)*q_nocor;
        J[781] -= dqdci;              /* dwdot[H]/d[CH3CHO] */
        J[783] -= dqdci;              /* dwdot[O]/d[CH3CHO] */
        J[784] += dqdci;              /* dwdot[OH]/d[CH3CHO] */
        /* d()/d[CH2CHO] */
        dqdci = (TB[16][14] - 1)*q_nocor;
        J[811] -= dqdci;              /* dwdot[H]/d[CH2CHO] */
        J[813] -= dqdci;              /* dwdot[O]/d[CH2CHO] */
        J[814] += dqdci;              /* dwdot[OH]/d[CH2CHO] */
        /* d()/d[C2H5O] */
        dqdci = (TB[16][15] - 1)*q_nocor;
        J[841] -= dqdci;              /* dwdot[H]/d[C2H5O] */
        J[843] -= dqdci;              /* dwdot[O]/d[C2H5O] */
        J[844] += dqdci;              /* dwdot[OH]/d[C2H5O] */
    }
    else {
        dqdc[0] = q_nocor;
        dqdc[1] = q_nocor + k_f*sc[3];
        dqdc[2] = TB[16][0]*q_nocor;
        dqdc[3] = q_nocor + k_f*sc[1];
        dqdc[4] = q_nocor - k_r;
        dqdc[5] = q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = TB[16][1]*q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = TB[16][3]*q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = q_nocor;
        dqdc[12] = TB[16][7]*q_nocor;
        dqdc[13] = TB[16][9]*q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = TB[16][16]*q_nocor;
        dqdc[18] = TB[16][2]*q_nocor;
        dqdc[19] = TB[16][4]*q_nocor;
        dqdc[20] = TB[16][5]*q_nocor;
        dqdc[21] = TB[16][6]*q_nocor;
        dqdc[22] = TB[16][8]*q_nocor;
        dqdc[23] = TB[16][10]*q_nocor;
        dqdc[24] = TB[16][11]*q_nocor;
        dqdc[25] = TB[16][12]*q_nocor;
        dqdc[26] = TB[16][13]*q_nocor;
        dqdc[27] = TB[16][14]*q_nocor;
        dqdc[28] = TB[16][15]*q_nocor;
        for (int k=0; k<29; k++) {
            J[30*k+1] -= dqdc[k];
            J[30*k+3] -= dqdc[k];
            J[30*k+4] += dqdc[k];
        }
    }
    J[871] -= dqdT; /* dwdot[H]/dT */
    J[873] -= dqdT; /* dwdot[O]/dT */
    J[874] += dqdT; /* dwdot[OH]/dT */

    /*reaction 18: 2.000000 O + M <=> O2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[17][0] - 1)*sc[2] + (TB[17][1] - 1)*sc[7] + (TB[17][2] - 1)*sc[18] + (TB[17][3] - 1)*sc[9] + (TB[17][4] - 1)*sc[19] + (TB[17][5] - 1)*sc[20] + (TB[17][6] - 1)*sc[21] + (TB[17][7] - 1)*sc[12] + (TB[17][8] - 1)*sc[22] + (TB[17][9] - 1)*sc[13] + (TB[17][10] - 1)*sc[23] + (TB[17][11] - 1)*sc[24] + (TB[17][12] - 1)*sc[25] + (TB[17][13] - 1)*sc[26] + (TB[17][14] - 1)*sc[27] + (TB[17][15] - 1)*sc[28] + (TB[17][16] - 1)*sc[17];
    /* forward */
    phi_f = pow(sc[3], 2.000000);
    k_f = prefactor_units[17] * fwd_A[17]
                * exp(fwd_beta[17] * tc[0] - activation_units[17] * fwd_Ea[17] * invT);
    dlnkfdT = fwd_beta[17] * invT + activation_units[17] * fwd_Ea[17] * invT2;
    /* reverse */
    phi_r = sc[5];
    Kc = refCinv * exp(2.000000*g_RT[3] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[3]) + (h_RT[5]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= 2 * q; /* O */
    wdot[5] += q; /* O2 */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[17][0] - 1)*q_nocor;
        J[63] += -2 * dqdci;          /* dwdot[O]/d[H2] */
        J[65] += dqdci;               /* dwdot[O2]/d[H2] */
        /* d()/d[O] */
        dqdci =  + k_f*2.000000*sc[3];
        J[93] += -2 * dqdci;          /* dwdot[O]/d[O] */
        J[95] += dqdci;               /* dwdot[O2]/d[O] */
        /* d()/d[O2] */
        dqdci =  - k_r;
        J[153] += -2 * dqdci;         /* dwdot[O]/d[O2] */
        J[155] += dqdci;              /* dwdot[O2]/d[O2] */
        /* d()/d[H2O] */
        dqdci = (TB[17][1] - 1)*q_nocor;
        J[213] += -2 * dqdci;         /* dwdot[O]/d[H2O] */
        J[215] += dqdci;              /* dwdot[O2]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[17][3] - 1)*q_nocor;
        J[273] += -2 * dqdci;         /* dwdot[O]/d[CO] */
        J[275] += dqdci;              /* dwdot[O2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[17][7] - 1)*q_nocor;
        J[363] += -2 * dqdci;         /* dwdot[O]/d[CO2] */
        J[365] += dqdci;              /* dwdot[O2]/d[CO2] */
        /* d()/d[CH4] */
        dqdci = (TB[17][9] - 1)*q_nocor;
        J[393] += -2 * dqdci;         /* dwdot[O]/d[CH4] */
        J[395] += dqdci;              /* dwdot[O2]/d[CH4] */
        /* d()/d[C2H6] */
        dqdci = (TB[17][16] - 1)*q_nocor;
        J[513] += -2 * dqdci;         /* dwdot[O]/d[C2H6] */
        J[515] += dqdci;              /* dwdot[O2]/d[C2H6] */
        /* d()/d[C] */
        dqdci = (TB[17][2] - 1)*q_nocor;
        J[543] += -2 * dqdci;         /* dwdot[O]/d[C] */
        J[545] += dqdci;              /* dwdot[O2]/d[C] */
        /* d()/d[CH] */
        dqdci = (TB[17][4] - 1)*q_nocor;
        J[573] += -2 * dqdci;         /* dwdot[O]/d[CH] */
        J[575] += dqdci;              /* dwdot[O2]/d[CH] */
        /* d()/d[HCO] */
        dqdci = (TB[17][5] - 1)*q_nocor;
        J[603] += -2 * dqdci;         /* dwdot[O]/d[HCO] */
        J[605] += dqdci;              /* dwdot[O2]/d[HCO] */
        /* d()/d[TXCH2] */
        dqdci = (TB[17][6] - 1)*q_nocor;
        J[633] += -2 * dqdci;         /* dwdot[O]/d[TXCH2] */
        J[635] += dqdci;              /* dwdot[O2]/d[TXCH2] */
        /* d()/d[SXCH2] */
        dqdci = (TB[17][8] - 1)*q_nocor;
        J[663] += -2 * dqdci;         /* dwdot[O]/d[SXCH2] */
        J[665] += dqdci;              /* dwdot[O2]/d[SXCH2] */
        /* d()/d[C2H3] */
        dqdci = (TB[17][10] - 1)*q_nocor;
        J[693] += -2 * dqdci;         /* dwdot[O]/d[C2H3] */
        J[695] += dqdci;              /* dwdot[O2]/d[C2H3] */
        /* d()/d[C2H5] */
        dqdci = (TB[17][11] - 1)*q_nocor;
        J[723] += -2 * dqdci;         /* dwdot[O]/d[C2H5] */
        J[725] += dqdci;              /* dwdot[O2]/d[C2H5] */
        /* d()/d[HCCO] */
        dqdci = (TB[17][12] - 1)*q_nocor;
        J[753] += -2 * dqdci;         /* dwdot[O]/d[HCCO] */
        J[755] += dqdci;              /* dwdot[O2]/d[HCCO] */
        /* d()/d[CH3CHO] */
        dqdci = (TB[17][13] - 1)*q_nocor;
        J[783] += -2 * dqdci;         /* dwdot[O]/d[CH3CHO] */
        J[785] += dqdci;              /* dwdot[O2]/d[CH3CHO] */
        /* d()/d[CH2CHO] */
        dqdci = (TB[17][14] - 1)*q_nocor;
        J[813] += -2 * dqdci;         /* dwdot[O]/d[CH2CHO] */
        J[815] += dqdci;              /* dwdot[O2]/d[CH2CHO] */
        /* d()/d[C2H5O] */
        dqdci = (TB[17][15] - 1)*q_nocor;
        J[843] += -2 * dqdci;         /* dwdot[O]/d[C2H5O] */
        J[845] += dqdci;              /* dwdot[O2]/d[C2H5O] */
    }
    else {
        dqdc[0] = q_nocor;
        dqdc[1] = q_nocor;
        dqdc[2] = TB[17][0]*q_nocor;
        dqdc[3] = q_nocor + k_f*2.000000*sc[3];
        dqdc[4] = q_nocor;
        dqdc[5] = q_nocor - k_r;
        dqdc[6] = q_nocor;
        dqdc[7] = TB[17][1]*q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = TB[17][3]*q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = q_nocor;
        dqdc[12] = TB[17][7]*q_nocor;
        dqdc[13] = TB[17][9]*q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = TB[17][16]*q_nocor;
        dqdc[18] = TB[17][2]*q_nocor;
        dqdc[19] = TB[17][4]*q_nocor;
        dqdc[20] = TB[17][5]*q_nocor;
        dqdc[21] = TB[17][6]*q_nocor;
        dqdc[22] = TB[17][8]*q_nocor;
        dqdc[23] = TB[17][10]*q_nocor;
        dqdc[24] = TB[17][11]*q_nocor;
        dqdc[25] = TB[17][12]*q_nocor;
        dqdc[26] = TB[17][13]*q_nocor;
        dqdc[27] = TB[17][14]*q_nocor;
        dqdc[28] = TB[17][15]*q_nocor;
        for (int k=0; k<29; k++) {
            J[30*k+3] += -2 * dqdc[k];
            J[30*k+5] += dqdc[k];
        }
    }
    J[873] += -2 * dqdT; /* dwdot[O]/dT */
    J[875] += dqdT; /* dwdot[O2]/dT */

    /*reaction 19: H + OH + M <=> H2O + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[18][0] - 1)*sc[2] + (TB[18][1] - 1)*sc[7] + (TB[18][2] - 1)*sc[18] + (TB[18][3] - 1)*sc[9] + (TB[18][4] - 1)*sc[19] + (TB[18][5] - 1)*sc[20] + (TB[18][6] - 1)*sc[21] + (TB[18][7] - 1)*sc[12] + (TB[18][8] - 1)*sc[22] + (TB[18][9] - 1)*sc[13] + (TB[18][10] - 1)*sc[23] + (TB[18][11] - 1)*sc[24] + (TB[18][12] - 1)*sc[25] + (TB[18][13] - 1)*sc[26] + (TB[18][14] - 1)*sc[27] + (TB[18][15] - 1)*sc[28] + (TB[18][16] - 1)*sc[17];
    /* forward */
    phi_f = sc[1]*sc[4];
    k_f = prefactor_units[18] * fwd_A[18]
                * exp(fwd_beta[18] * tc[0] - activation_units[18] * fwd_Ea[18] * invT);
    dlnkfdT = fwd_beta[18] * invT + activation_units[18] * fwd_Ea[18] * invT2;
    /* reverse */
    phi_r = sc[7];
    Kc = refCinv * exp(g_RT[1] + g_RT[4] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[4]) + (h_RT[7]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[4] -= q; /* OH */
    wdot[7] += q; /* H2O */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[4];
        J[31] -= dqdci;               /* dwdot[H]/d[H] */
        J[34] -= dqdci;               /* dwdot[OH]/d[H] */
        J[37] += dqdci;               /* dwdot[H2O]/d[H] */
        /* d()/d[H2] */
        dqdci = (TB[18][0] - 1)*q_nocor;
        J[61] -= dqdci;               /* dwdot[H]/d[H2] */
        J[64] -= dqdci;               /* dwdot[OH]/d[H2] */
        J[67] += dqdci;               /* dwdot[H2O]/d[H2] */
        /* d()/d[OH] */
        dqdci =  + k_f*sc[1];
        J[121] -= dqdci;              /* dwdot[H]/d[OH] */
        J[124] -= dqdci;              /* dwdot[OH]/d[OH] */
        J[127] += dqdci;              /* dwdot[H2O]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (TB[18][1] - 1)*q_nocor - k_r;
        J[211] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[214] -= dqdci;              /* dwdot[OH]/d[H2O] */
        J[217] += dqdci;              /* dwdot[H2O]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[18][3] - 1)*q_nocor;
        J[271] -= dqdci;              /* dwdot[H]/d[CO] */
        J[274] -= dqdci;              /* dwdot[OH]/d[CO] */
        J[277] += dqdci;              /* dwdot[H2O]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[18][7] - 1)*q_nocor;
        J[361] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[364] -= dqdci;              /* dwdot[OH]/d[CO2] */
        J[367] += dqdci;              /* dwdot[H2O]/d[CO2] */
        /* d()/d[CH4] */
        dqdci = (TB[18][9] - 1)*q_nocor;
        J[391] -= dqdci;              /* dwdot[H]/d[CH4] */
        J[394] -= dqdci;              /* dwdot[OH]/d[CH4] */
        J[397] += dqdci;              /* dwdot[H2O]/d[CH4] */
        /* d()/d[C2H6] */
        dqdci = (TB[18][16] - 1)*q_nocor;
        J[511] -= dqdci;              /* dwdot[H]/d[C2H6] */
        J[514] -= dqdci;              /* dwdot[OH]/d[C2H6] */
        J[517] += dqdci;              /* dwdot[H2O]/d[C2H6] */
        /* d()/d[C] */
        dqdci = (TB[18][2] - 1)*q_nocor;
        J[541] -= dqdci;              /* dwdot[H]/d[C] */
        J[544] -= dqdci;              /* dwdot[OH]/d[C] */
        J[547] += dqdci;              /* dwdot[H2O]/d[C] */
        /* d()/d[CH] */
        dqdci = (TB[18][4] - 1)*q_nocor;
        J[571] -= dqdci;              /* dwdot[H]/d[CH] */
        J[574] -= dqdci;              /* dwdot[OH]/d[CH] */
        J[577] += dqdci;              /* dwdot[H2O]/d[CH] */
        /* d()/d[HCO] */
        dqdci = (TB[18][5] - 1)*q_nocor;
        J[601] -= dqdci;              /* dwdot[H]/d[HCO] */
        J[604] -= dqdci;              /* dwdot[OH]/d[HCO] */
        J[607] += dqdci;              /* dwdot[H2O]/d[HCO] */
        /* d()/d[TXCH2] */
        dqdci = (TB[18][6] - 1)*q_nocor;
        J[631] -= dqdci;              /* dwdot[H]/d[TXCH2] */
        J[634] -= dqdci;              /* dwdot[OH]/d[TXCH2] */
        J[637] += dqdci;              /* dwdot[H2O]/d[TXCH2] */
        /* d()/d[SXCH2] */
        dqdci = (TB[18][8] - 1)*q_nocor;
        J[661] -= dqdci;              /* dwdot[H]/d[SXCH2] */
        J[664] -= dqdci;              /* dwdot[OH]/d[SXCH2] */
        J[667] += dqdci;              /* dwdot[H2O]/d[SXCH2] */
        /* d()/d[C2H3] */
        dqdci = (TB[18][10] - 1)*q_nocor;
        J[691] -= dqdci;              /* dwdot[H]/d[C2H3] */
        J[694] -= dqdci;              /* dwdot[OH]/d[C2H3] */
        J[697] += dqdci;              /* dwdot[H2O]/d[C2H3] */
        /* d()/d[C2H5] */
        dqdci = (TB[18][11] - 1)*q_nocor;
        J[721] -= dqdci;              /* dwdot[H]/d[C2H5] */
        J[724] -= dqdci;              /* dwdot[OH]/d[C2H5] */
        J[727] += dqdci;              /* dwdot[H2O]/d[C2H5] */
        /* d()/d[HCCO] */
        dqdci = (TB[18][12] - 1)*q_nocor;
        J[751] -= dqdci;              /* dwdot[H]/d[HCCO] */
        J[754] -= dqdci;              /* dwdot[OH]/d[HCCO] */
        J[757] += dqdci;              /* dwdot[H2O]/d[HCCO] */
        /* d()/d[CH3CHO] */
        dqdci = (TB[18][13] - 1)*q_nocor;
        J[781] -= dqdci;              /* dwdot[H]/d[CH3CHO] */
        J[784] -= dqdci;              /* dwdot[OH]/d[CH3CHO] */
        J[787] += dqdci;              /* dwdot[H2O]/d[CH3CHO] */
        /* d()/d[CH2CHO] */
        dqdci = (TB[18][14] - 1)*q_nocor;
        J[811] -= dqdci;              /* dwdot[H]/d[CH2CHO] */
        J[814] -= dqdci;              /* dwdot[OH]/d[CH2CHO] */
        J[817] += dqdci;              /* dwdot[H2O]/d[CH2CHO] */
        /* d()/d[C2H5O] */
        dqdci = (TB[18][15] - 1)*q_nocor;
        J[841] -= dqdci;              /* dwdot[H]/d[C2H5O] */
        J[844] -= dqdci;              /* dwdot[OH]/d[C2H5O] */
        J[847] += dqdci;              /* dwdot[H2O]/d[C2H5O] */
    }
    else {
        dqdc[0] = q_nocor;
        dqdc[1] = q_nocor + k_f*sc[4];
        dqdc[2] = TB[18][0]*q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor + k_f*sc[1];
        dqdc[5] = q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = TB[18][1]*q_nocor - k_r;
        dqdc[8] = q_nocor;
        dqdc[9] = TB[18][3]*q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = q_nocor;
        dqdc[12] = TB[18][7]*q_nocor;
        dqdc[13] = TB[18][9]*q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = TB[18][16]*q_nocor;
        dqdc[18] = TB[18][2]*q_nocor;
        dqdc[19] = TB[18][4]*q_nocor;
        dqdc[20] = TB[18][5]*q_nocor;
        dqdc[21] = TB[18][6]*q_nocor;
        dqdc[22] = TB[18][8]*q_nocor;
        dqdc[23] = TB[18][10]*q_nocor;
        dqdc[24] = TB[18][11]*q_nocor;
        dqdc[25] = TB[18][12]*q_nocor;
        dqdc[26] = TB[18][13]*q_nocor;
        dqdc[27] = TB[18][14]*q_nocor;
        dqdc[28] = TB[18][15]*q_nocor;
        for (int k=0; k<29; k++) {
            J[30*k+1] -= dqdc[k];
            J[30*k+4] -= dqdc[k];
            J[30*k+7] += dqdc[k];
        }
    }
    J[871] -= dqdT; /* dwdot[H]/dT */
    J[874] -= dqdT; /* dwdot[OH]/dT */
    J[877] += dqdT; /* dwdot[H2O]/dT */

    /*reaction 20: HCO + M <=> CO + H + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[19][0] - 1)*sc[2] + (TB[19][1] - 1)*sc[7] + (TB[19][2] - 1)*sc[18] + (TB[19][3] - 1)*sc[9] + (TB[19][4] - 1)*sc[19] + (TB[19][5] - 1)*sc[20] + (TB[19][6] - 1)*sc[21] + (TB[19][7] - 1)*sc[12] + (TB[19][8] - 1)*sc[22] + (TB[19][9] - 1)*sc[13] + (TB[19][10] - 1)*sc[23] + (TB[19][11] - 1)*sc[24] + (TB[19][12] - 1)*sc[25] + (TB[19][13] - 1)*sc[26] + (TB[19][14] - 1)*sc[27] + (TB[19][15] - 1)*sc[28] + (TB[19][16] - 1)*sc[17];
    /* forward */
    phi_f = sc[20];
    k_f = prefactor_units[19] * fwd_A[19]
                * exp(fwd_beta[19] * tc[0] - activation_units[19] * fwd_Ea[19] * invT);
    dlnkfdT = fwd_beta[19] * invT + activation_units[19] * fwd_Ea[19] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[9];
    Kc = refC * exp(-g_RT[1] - g_RT[9] + g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[20]) + (h_RT[1] + h_RT[9]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[9] += q; /* CO */
    wdot[20] -= q; /* HCO */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H] */
        dqdci =  - k_r*sc[9];
        J[31] += dqdci;               /* dwdot[H]/d[H] */
        J[39] += dqdci;               /* dwdot[CO]/d[H] */
        J[50] -= dqdci;               /* dwdot[HCO]/d[H] */
        /* d()/d[H2] */
        dqdci = (TB[19][0] - 1)*q_nocor;
        J[61] += dqdci;               /* dwdot[H]/d[H2] */
        J[69] += dqdci;               /* dwdot[CO]/d[H2] */
        J[80] -= dqdci;               /* dwdot[HCO]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[19][1] - 1)*q_nocor;
        J[211] += dqdci;              /* dwdot[H]/d[H2O] */
        J[219] += dqdci;              /* dwdot[CO]/d[H2O] */
        J[230] -= dqdci;              /* dwdot[HCO]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[19][3] - 1)*q_nocor - k_r*sc[1];
        J[271] += dqdci;              /* dwdot[H]/d[CO] */
        J[279] += dqdci;              /* dwdot[CO]/d[CO] */
        J[290] -= dqdci;              /* dwdot[HCO]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[19][7] - 1)*q_nocor;
        J[361] += dqdci;              /* dwdot[H]/d[CO2] */
        J[369] += dqdci;              /* dwdot[CO]/d[CO2] */
        J[380] -= dqdci;              /* dwdot[HCO]/d[CO2] */
        /* d()/d[CH4] */
        dqdci = (TB[19][9] - 1)*q_nocor;
        J[391] += dqdci;              /* dwdot[H]/d[CH4] */
        J[399] += dqdci;              /* dwdot[CO]/d[CH4] */
        J[410] -= dqdci;              /* dwdot[HCO]/d[CH4] */
        /* d()/d[C2H6] */
        dqdci = (TB[19][16] - 1)*q_nocor;
        J[511] += dqdci;              /* dwdot[H]/d[C2H6] */
        J[519] += dqdci;              /* dwdot[CO]/d[C2H6] */
        J[530] -= dqdci;              /* dwdot[HCO]/d[C2H6] */
        /* d()/d[C] */
        dqdci = (TB[19][2] - 1)*q_nocor;
        J[541] += dqdci;              /* dwdot[H]/d[C] */
        J[549] += dqdci;              /* dwdot[CO]/d[C] */
        J[560] -= dqdci;              /* dwdot[HCO]/d[C] */
        /* d()/d[CH] */
        dqdci = (TB[19][4] - 1)*q_nocor;
        J[571] += dqdci;              /* dwdot[H]/d[CH] */
        J[579] += dqdci;              /* dwdot[CO]/d[CH] */
        J[590] -= dqdci;              /* dwdot[HCO]/d[CH] */
        /* d()/d[HCO] */
        dqdci = (TB[19][5] - 1)*q_nocor + k_f;
        J[601] += dqdci;              /* dwdot[H]/d[HCO] */
        J[609] += dqdci;              /* dwdot[CO]/d[HCO] */
        J[620] -= dqdci;              /* dwdot[HCO]/d[HCO] */
        /* d()/d[TXCH2] */
        dqdci = (TB[19][6] - 1)*q_nocor;
        J[631] += dqdci;              /* dwdot[H]/d[TXCH2] */
        J[639] += dqdci;              /* dwdot[CO]/d[TXCH2] */
        J[650] -= dqdci;              /* dwdot[HCO]/d[TXCH2] */
        /* d()/d[SXCH2] */
        dqdci = (TB[19][8] - 1)*q_nocor;
        J[661] += dqdci;              /* dwdot[H]/d[SXCH2] */
        J[669] += dqdci;              /* dwdot[CO]/d[SXCH2] */
        J[680] -= dqdci;              /* dwdot[HCO]/d[SXCH2] */
        /* d()/d[C2H3] */
        dqdci = (TB[19][10] - 1)*q_nocor;
        J[691] += dqdci;              /* dwdot[H]/d[C2H3] */
        J[699] += dqdci;              /* dwdot[CO]/d[C2H3] */
        J[710] -= dqdci;              /* dwdot[HCO]/d[C2H3] */
        /* d()/d[C2H5] */
        dqdci = (TB[19][11] - 1)*q_nocor;
        J[721] += dqdci;              /* dwdot[H]/d[C2H5] */
        J[729] += dqdci;              /* dwdot[CO]/d[C2H5] */
        J[740] -= dqdci;              /* dwdot[HCO]/d[C2H5] */
        /* d()/d[HCCO] */
        dqdci = (TB[19][12] - 1)*q_nocor;
        J[751] += dqdci;              /* dwdot[H]/d[HCCO] */
        J[759] += dqdci;              /* dwdot[CO]/d[HCCO] */
        J[770] -= dqdci;              /* dwdot[HCO]/d[HCCO] */
        /* d()/d[CH3CHO] */
        dqdci = (TB[19][13] - 1)*q_nocor;
        J[781] += dqdci;              /* dwdot[H]/d[CH3CHO] */
        J[789] += dqdci;              /* dwdot[CO]/d[CH3CHO] */
        J[800] -= dqdci;              /* dwdot[HCO]/d[CH3CHO] */
        /* d()/d[CH2CHO] */
        dqdci = (TB[19][14] - 1)*q_nocor;
        J[811] += dqdci;              /* dwdot[H]/d[CH2CHO] */
        J[819] += dqdci;              /* dwdot[CO]/d[CH2CHO] */
        J[830] -= dqdci;              /* dwdot[HCO]/d[CH2CHO] */
        /* d()/d[C2H5O] */
        dqdci = (TB[19][15] - 1)*q_nocor;
        J[841] += dqdci;              /* dwdot[H]/d[C2H5O] */
        J[849] += dqdci;              /* dwdot[CO]/d[C2H5O] */
        J[860] -= dqdci;              /* dwdot[HCO]/d[C2H5O] */
    }
    else {
        dqdc[0] = q_nocor;
        dqdc[1] = q_nocor - k_r*sc[9];
        dqdc[2] = TB[19][0]*q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor;
        dqdc[5] = q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = TB[19][1]*q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = TB[19][3]*q_nocor - k_r*sc[1];
        dqdc[10] = q_nocor;
        dqdc[11] = q_nocor;
        dqdc[12] = TB[19][7]*q_nocor;
        dqdc[13] = TB[19][9]*q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = TB[19][16]*q_nocor;
        dqdc[18] = TB[19][2]*q_nocor;
        dqdc[19] = TB[19][4]*q_nocor;
        dqdc[20] = TB[19][5]*q_nocor + k_f;
        dqdc[21] = TB[19][6]*q_nocor;
        dqdc[22] = TB[19][8]*q_nocor;
        dqdc[23] = TB[19][10]*q_nocor;
        dqdc[24] = TB[19][11]*q_nocor;
        dqdc[25] = TB[19][12]*q_nocor;
        dqdc[26] = TB[19][13]*q_nocor;
        dqdc[27] = TB[19][14]*q_nocor;
        dqdc[28] = TB[19][15]*q_nocor;
        for (int k=0; k<29; k++) {
            J[30*k+1] += dqdc[k];
            J[30*k+9] += dqdc[k];
            J[30*k+20] -= dqdc[k];
        }
    }
    J[871] += dqdT; /* dwdot[H]/dT */
    J[879] += dqdT; /* dwdot[CO]/dT */
    J[890] -= dqdT; /* dwdot[HCO]/dT */

    /*reaction 21: 2.000000 H + H2 <=> 2.000000 H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[1], 2.000000)*sc[2];
    k_f = prefactor_units[20] * fwd_A[20]
                * exp(fwd_beta[20] * tc[0] - activation_units[20] * fwd_Ea[20] * invT);
    dlnkfdT = fwd_beta[20] * invT + activation_units[20] * fwd_Ea[20] * invT2;
    /* reverse */
    phi_r = pow(sc[2], 2.000000);
    Kc = refCinv * exp(2.000000*g_RT[1] + g_RT[2] - 2.000000*g_RT[2]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[1] + h_RT[2]) + (2.000000*h_RT[2]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= 2 * q; /* H */
    wdot[2] += q; /* H2 */
    /* d()/d[H] */
    dqdci =  + k_f*2.000000*sc[1]*sc[2];
    J[31] += -2 * dqdci;          /* dwdot[H]/d[H] */
    J[32] += dqdci;               /* dwdot[H2]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*pow(sc[1], 2.000000) - k_r*2.000000*sc[2];
    J[61] += -2 * dqdci;          /* dwdot[H]/d[H2] */
    J[62] += dqdci;               /* dwdot[H2]/d[H2] */
    /* d()/dT */
    J[871] += -2 * dqdT;          /* dwdot[H]/dT */
    J[872] += dqdT;               /* dwdot[H2]/dT */

    /*reaction 22: O + H2 <=> H + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[3];
    k_f = prefactor_units[21] * fwd_A[21]
                * exp(fwd_beta[21] * tc[0] - activation_units[21] * fwd_Ea[21] * invT);
    dlnkfdT = fwd_beta[21] * invT + activation_units[21] * fwd_Ea[21] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[4];
    Kc = exp(-g_RT[1] + g_RT[2] + g_RT[3] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[3]) + (h_RT[1] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[2] -= q; /* H2 */
    wdot[3] -= q; /* O */
    wdot[4] += q; /* OH */
    /* d()/d[H] */
    dqdci =  - k_r*sc[4];
    J[31] += dqdci;               /* dwdot[H]/d[H] */
    J[32] -= dqdci;               /* dwdot[H2]/d[H] */
    J[33] -= dqdci;               /* dwdot[O]/d[H] */
    J[34] += dqdci;               /* dwdot[OH]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[3];
    J[61] += dqdci;               /* dwdot[H]/d[H2] */
    J[62] -= dqdci;               /* dwdot[H2]/d[H2] */
    J[63] -= dqdci;               /* dwdot[O]/d[H2] */
    J[64] += dqdci;               /* dwdot[OH]/d[H2] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[2];
    J[91] += dqdci;               /* dwdot[H]/d[O] */
    J[92] -= dqdci;               /* dwdot[H2]/d[O] */
    J[93] -= dqdci;               /* dwdot[O]/d[O] */
    J[94] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[1];
    J[121] += dqdci;              /* dwdot[H]/d[OH] */
    J[122] -= dqdci;              /* dwdot[H2]/d[OH] */
    J[123] -= dqdci;              /* dwdot[O]/d[OH] */
    J[124] += dqdci;              /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[871] += dqdT;               /* dwdot[H]/dT */
    J[872] -= dqdT;               /* dwdot[H2]/dT */
    J[873] -= dqdT;               /* dwdot[O]/dT */
    J[874] += dqdT;               /* dwdot[OH]/dT */

    /*reaction 23: 2.000000 OH <=> O + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[4], 2.000000);
    k_f = prefactor_units[22] * fwd_A[22]
                * exp(fwd_beta[22] * tc[0] - activation_units[22] * fwd_Ea[22] * invT);
    dlnkfdT = fwd_beta[22] * invT + activation_units[22] * fwd_Ea[22] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[7];
    Kc = exp(-g_RT[3] + 2.000000*g_RT[4] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[4]) + (h_RT[3] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += q; /* O */
    wdot[4] -= 2 * q; /* OH */
    wdot[7] += q; /* H2O */
    /* d()/d[O] */
    dqdci =  - k_r*sc[7];
    J[93] += dqdci;               /* dwdot[O]/d[O] */
    J[94] += -2 * dqdci;          /* dwdot[OH]/d[O] */
    J[97] += dqdci;               /* dwdot[H2O]/d[O] */
    /* d()/d[OH] */
    dqdci =  + k_f*2.000000*sc[4];
    J[123] += dqdci;              /* dwdot[O]/d[OH] */
    J[124] += -2 * dqdci;         /* dwdot[OH]/d[OH] */
    J[127] += dqdci;              /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[3];
    J[213] += dqdci;              /* dwdot[O]/d[H2O] */
    J[214] += -2 * dqdci;         /* dwdot[OH]/d[H2O] */
    J[217] += dqdci;              /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[873] += dqdT;               /* dwdot[O]/dT */
    J[874] += -2 * dqdT;          /* dwdot[OH]/dT */
    J[877] += dqdT;               /* dwdot[H2O]/dT */

    /*reaction 24: OH + H2 <=> H + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[4];
    k_f = prefactor_units[23] * fwd_A[23]
                * exp(fwd_beta[23] * tc[0] - activation_units[23] * fwd_Ea[23] * invT);
    dlnkfdT = fwd_beta[23] * invT + activation_units[23] * fwd_Ea[23] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[7];
    Kc = exp(-g_RT[1] + g_RT[2] + g_RT[4] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[4]) + (h_RT[1] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[2] -= q; /* H2 */
    wdot[4] -= q; /* OH */
    wdot[7] += q; /* H2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[7];
    J[31] += dqdci;               /* dwdot[H]/d[H] */
    J[32] -= dqdci;               /* dwdot[H2]/d[H] */
    J[34] -= dqdci;               /* dwdot[OH]/d[H] */
    J[37] += dqdci;               /* dwdot[H2O]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[4];
    J[61] += dqdci;               /* dwdot[H]/d[H2] */
    J[62] -= dqdci;               /* dwdot[H2]/d[H2] */
    J[64] -= dqdci;               /* dwdot[OH]/d[H2] */
    J[67] += dqdci;               /* dwdot[H2O]/d[H2] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[2];
    J[121] += dqdci;              /* dwdot[H]/d[OH] */
    J[122] -= dqdci;              /* dwdot[H2]/d[OH] */
    J[124] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[127] += dqdci;              /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[1];
    J[211] += dqdci;              /* dwdot[H]/d[H2O] */
    J[212] -= dqdci;              /* dwdot[H2]/d[H2O] */
    J[214] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[217] += dqdci;              /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[871] += dqdT;               /* dwdot[H]/dT */
    J[872] -= dqdT;               /* dwdot[H2]/dT */
    J[874] -= dqdT;               /* dwdot[OH]/dT */
    J[877] += dqdT;               /* dwdot[H2O]/dT */

    /*reaction 25: 2.000000 H + H2O <=> H2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[1], 2.000000)*sc[7];
    k_f = prefactor_units[24] * fwd_A[24]
                * exp(fwd_beta[24] * tc[0] - activation_units[24] * fwd_Ea[24] * invT);
    dlnkfdT = fwd_beta[24] * invT + activation_units[24] * fwd_Ea[24] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[7];
    Kc = refCinv * exp(2.000000*g_RT[1] - g_RT[2] + g_RT[7] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[1] + h_RT[7]) + (h_RT[2] + h_RT[7]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= 2 * q; /* H */
    wdot[2] += q; /* H2 */
    /* d()/d[H] */
    dqdci =  + k_f*2.000000*sc[1]*sc[7];
    J[31] += -2 * dqdci;          /* dwdot[H]/d[H] */
    J[32] += dqdci;               /* dwdot[H2]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[7];
    J[61] += -2 * dqdci;          /* dwdot[H]/d[H2] */
    J[62] += dqdci;               /* dwdot[H2]/d[H2] */
    /* d()/d[H2O] */
    dqdci =  + k_f*pow(sc[1], 2.000000) - k_r*sc[2];
    J[211] += -2 * dqdci;         /* dwdot[H]/d[H2O] */
    J[212] += dqdci;              /* dwdot[H2]/d[H2O] */
    /* d()/dT */
    J[871] += -2 * dqdT;          /* dwdot[H]/dT */
    J[872] += dqdT;               /* dwdot[H2]/dT */

    /*reaction 26: H + O2 <=> O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[5];
    k_f = prefactor_units[25] * fwd_A[25]
                * exp(fwd_beta[25] * tc[0] - activation_units[25] * fwd_Ea[25] * invT);
    dlnkfdT = fwd_beta[25] * invT + activation_units[25] * fwd_Ea[25] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[4];
    Kc = exp(g_RT[1] - g_RT[3] - g_RT[4] + g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[5]) + (h_RT[3] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[3] += q; /* O */
    wdot[4] += q; /* OH */
    wdot[5] -= q; /* O2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[5];
    J[31] -= dqdci;               /* dwdot[H]/d[H] */
    J[33] += dqdci;               /* dwdot[O]/d[H] */
    J[34] += dqdci;               /* dwdot[OH]/d[H] */
    J[35] -= dqdci;               /* dwdot[O2]/d[H] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[4];
    J[91] -= dqdci;               /* dwdot[H]/d[O] */
    J[93] += dqdci;               /* dwdot[O]/d[O] */
    J[94] += dqdci;               /* dwdot[OH]/d[O] */
    J[95] -= dqdci;               /* dwdot[O2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[3];
    J[121] -= dqdci;              /* dwdot[H]/d[OH] */
    J[123] += dqdci;              /* dwdot[O]/d[OH] */
    J[124] += dqdci;              /* dwdot[OH]/d[OH] */
    J[125] -= dqdci;              /* dwdot[O2]/d[OH] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1];
    J[151] -= dqdci;              /* dwdot[H]/d[O2] */
    J[153] += dqdci;              /* dwdot[O]/d[O2] */
    J[154] += dqdci;              /* dwdot[OH]/d[O2] */
    J[155] -= dqdci;              /* dwdot[O2]/d[O2] */
    /* d()/dT */
    J[871] -= dqdT;               /* dwdot[H]/dT */
    J[873] += dqdT;               /* dwdot[O]/dT */
    J[874] += dqdT;               /* dwdot[OH]/dT */
    J[875] -= dqdT;               /* dwdot[O2]/dT */

    /*reaction 27: H2 + O2 <=> HO2 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[5];
    k_f = prefactor_units[26] * fwd_A[26]
                * exp(fwd_beta[26] * tc[0] - activation_units[26] * fwd_Ea[26] * invT);
    dlnkfdT = fwd_beta[26] * invT + activation_units[26] * fwd_Ea[26] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[8];
    Kc = exp(-g_RT[1] + g_RT[2] + g_RT[5] - g_RT[8]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[5]) + (h_RT[1] + h_RT[8]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[2] -= q; /* H2 */
    wdot[5] -= q; /* O2 */
    wdot[8] += q; /* HO2 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[8];
    J[31] += dqdci;               /* dwdot[H]/d[H] */
    J[32] -= dqdci;               /* dwdot[H2]/d[H] */
    J[35] -= dqdci;               /* dwdot[O2]/d[H] */
    J[38] += dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[5];
    J[61] += dqdci;               /* dwdot[H]/d[H2] */
    J[62] -= dqdci;               /* dwdot[H2]/d[H2] */
    J[65] -= dqdci;               /* dwdot[O2]/d[H2] */
    J[68] += dqdci;               /* dwdot[HO2]/d[H2] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[2];
    J[151] += dqdci;              /* dwdot[H]/d[O2] */
    J[152] -= dqdci;              /* dwdot[H2]/d[O2] */
    J[155] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[158] += dqdci;              /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[1];
    J[241] += dqdci;              /* dwdot[H]/d[HO2] */
    J[242] -= dqdci;              /* dwdot[H2]/d[HO2] */
    J[245] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[248] += dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[871] += dqdT;               /* dwdot[H]/dT */
    J[872] -= dqdT;               /* dwdot[H2]/dT */
    J[875] -= dqdT;               /* dwdot[O2]/dT */
    J[878] += dqdT;               /* dwdot[HO2]/dT */

    /*reaction 28: HO2 + OH <=> H2O + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[8];
    k_f = prefactor_units[27] * fwd_A[27]
                * exp(fwd_beta[27] * tc[0] - activation_units[27] * fwd_Ea[27] * invT);
    dlnkfdT = fwd_beta[27] * invT + activation_units[27] * fwd_Ea[27] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[7];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[7] + g_RT[8]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[8]) + (h_RT[5] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* O2 */
    wdot[7] += q; /* H2O */
    wdot[8] -= q; /* HO2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[8];
    J[124] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[125] += dqdci;              /* dwdot[O2]/d[OH] */
    J[127] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[128] -= dqdci;              /* dwdot[HO2]/d[OH] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[7];
    J[154] -= dqdci;              /* dwdot[OH]/d[O2] */
    J[155] += dqdci;              /* dwdot[O2]/d[O2] */
    J[157] += dqdci;              /* dwdot[H2O]/d[O2] */
    J[158] -= dqdci;              /* dwdot[HO2]/d[O2] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[5];
    J[214] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[215] += dqdci;              /* dwdot[O2]/d[H2O] */
    J[217] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[218] -= dqdci;              /* dwdot[HO2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[4];
    J[244] -= dqdci;              /* dwdot[OH]/d[HO2] */
    J[245] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[247] += dqdci;              /* dwdot[H2O]/d[HO2] */
    J[248] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[874] -= dqdT;               /* dwdot[OH]/dT */
    J[875] += dqdT;               /* dwdot[O2]/dT */
    J[877] += dqdT;               /* dwdot[H2O]/dT */
    J[878] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 29: HO2 + H <=> 2.000000 OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[8];
    k_f = prefactor_units[28] * fwd_A[28]
                * exp(fwd_beta[28] * tc[0] - activation_units[28] * fwd_Ea[28] * invT);
    dlnkfdT = fwd_beta[28] * invT + activation_units[28] * fwd_Ea[28] * invT2;
    /* reverse */
    phi_r = pow(sc[4], 2.000000);
    Kc = exp(g_RT[1] - 2.000000*g_RT[4] + g_RT[8]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[8]) + (2.000000*h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[4] += 2 * q; /* OH */
    wdot[8] -= q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[8];
    J[31] -= dqdci;               /* dwdot[H]/d[H] */
    J[34] += 2 * dqdci;           /* dwdot[OH]/d[H] */
    J[38] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*2.000000*sc[4];
    J[121] -= dqdci;              /* dwdot[H]/d[OH] */
    J[124] += 2 * dqdci;          /* dwdot[OH]/d[OH] */
    J[128] -= dqdci;              /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[1];
    J[241] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[244] += 2 * dqdci;          /* dwdot[OH]/d[HO2] */
    J[248] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[871] -= dqdT;               /* dwdot[H]/dT */
    J[874] += 2 * dqdT;           /* dwdot[OH]/dT */
    J[878] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 30: HO2 + O <=> OH + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[8];
    k_f = prefactor_units[29] * fwd_A[29]
                * exp(fwd_beta[29] * tc[0] - activation_units[29] * fwd_Ea[29] * invT);
    dlnkfdT = fwd_beta[29] * invT + activation_units[29] * fwd_Ea[29] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[5];
    Kc = exp(g_RT[3] - g_RT[4] - g_RT[5] + g_RT[8]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[8]) + (h_RT[4] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O */
    wdot[4] += q; /* OH */
    wdot[5] += q; /* O2 */
    wdot[8] -= q; /* HO2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[8];
    J[93] -= dqdci;               /* dwdot[O]/d[O] */
    J[94] += dqdci;               /* dwdot[OH]/d[O] */
    J[95] += dqdci;               /* dwdot[O2]/d[O] */
    J[98] -= dqdci;               /* dwdot[HO2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[5];
    J[123] -= dqdci;              /* dwdot[O]/d[OH] */
    J[124] += dqdci;              /* dwdot[OH]/d[OH] */
    J[125] += dqdci;              /* dwdot[O2]/d[OH] */
    J[128] -= dqdci;              /* dwdot[HO2]/d[OH] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[4];
    J[153] -= dqdci;              /* dwdot[O]/d[O2] */
    J[154] += dqdci;              /* dwdot[OH]/d[O2] */
    J[155] += dqdci;              /* dwdot[O2]/d[O2] */
    J[158] -= dqdci;              /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[3];
    J[243] -= dqdci;              /* dwdot[O]/d[HO2] */
    J[244] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[245] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[248] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[873] -= dqdT;               /* dwdot[O]/dT */
    J[874] += dqdT;               /* dwdot[OH]/dT */
    J[875] += dqdT;               /* dwdot[O2]/dT */
    J[878] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 31: HO2 + H <=> O + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[8];
    k_f = prefactor_units[30] * fwd_A[30]
                * exp(fwd_beta[30] * tc[0] - activation_units[30] * fwd_Ea[30] * invT);
    dlnkfdT = fwd_beta[30] * invT + activation_units[30] * fwd_Ea[30] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[7];
    Kc = exp(g_RT[1] - g_RT[3] - g_RT[7] + g_RT[8]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[8]) + (h_RT[3] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[3] += q; /* O */
    wdot[7] += q; /* H2O */
    wdot[8] -= q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[8];
    J[31] -= dqdci;               /* dwdot[H]/d[H] */
    J[33] += dqdci;               /* dwdot[O]/d[H] */
    J[37] += dqdci;               /* dwdot[H2O]/d[H] */
    J[38] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[7];
    J[91] -= dqdci;               /* dwdot[H]/d[O] */
    J[93] += dqdci;               /* dwdot[O]/d[O] */
    J[97] += dqdci;               /* dwdot[H2O]/d[O] */
    J[98] -= dqdci;               /* dwdot[HO2]/d[O] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[3];
    J[211] -= dqdci;              /* dwdot[H]/d[H2O] */
    J[213] += dqdci;              /* dwdot[O]/d[H2O] */
    J[217] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[218] -= dqdci;              /* dwdot[HO2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[1];
    J[241] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[243] += dqdci;              /* dwdot[O]/d[HO2] */
    J[247] += dqdci;              /* dwdot[H2O]/d[HO2] */
    J[248] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[871] -= dqdT;               /* dwdot[H]/dT */
    J[873] += dqdT;               /* dwdot[O]/dT */
    J[877] += dqdT;               /* dwdot[H2O]/dT */
    J[878] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 32: HO2 + OH <=> H2O + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[8];
    k_f = prefactor_units[31] * fwd_A[31]
                * exp(fwd_beta[31] * tc[0] - activation_units[31] * fwd_Ea[31] * invT);
    dlnkfdT = fwd_beta[31] * invT + activation_units[31] * fwd_Ea[31] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[7];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[7] + g_RT[8]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[8]) + (h_RT[5] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* O2 */
    wdot[7] += q; /* H2O */
    wdot[8] -= q; /* HO2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[8];
    J[124] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[125] += dqdci;              /* dwdot[O2]/d[OH] */
    J[127] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[128] -= dqdci;              /* dwdot[HO2]/d[OH] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[7];
    J[154] -= dqdci;              /* dwdot[OH]/d[O2] */
    J[155] += dqdci;              /* dwdot[O2]/d[O2] */
    J[157] += dqdci;              /* dwdot[H2O]/d[O2] */
    J[158] -= dqdci;              /* dwdot[HO2]/d[O2] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[5];
    J[214] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[215] += dqdci;              /* dwdot[O2]/d[H2O] */
    J[217] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[218] -= dqdci;              /* dwdot[HO2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[4];
    J[244] -= dqdci;              /* dwdot[OH]/d[HO2] */
    J[245] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[247] += dqdci;              /* dwdot[H2O]/d[HO2] */
    J[248] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[874] -= dqdT;               /* dwdot[OH]/dT */
    J[875] += dqdT;               /* dwdot[O2]/dT */
    J[877] += dqdT;               /* dwdot[H2O]/dT */
    J[878] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 33: H2O2 + O <=> HO2 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[6];
    k_f = prefactor_units[32] * fwd_A[32]
                * exp(fwd_beta[32] * tc[0] - activation_units[32] * fwd_Ea[32] * invT);
    dlnkfdT = fwd_beta[32] * invT + activation_units[32] * fwd_Ea[32] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[8];
    Kc = exp(g_RT[3] - g_RT[4] + g_RT[6] - g_RT[8]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[6]) + (h_RT[4] + h_RT[8]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O */
    wdot[4] += q; /* OH */
    wdot[6] -= q; /* H2O2 */
    wdot[8] += q; /* HO2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[6];
    J[93] -= dqdci;               /* dwdot[O]/d[O] */
    J[94] += dqdci;               /* dwdot[OH]/d[O] */
    J[96] -= dqdci;               /* dwdot[H2O2]/d[O] */
    J[98] += dqdci;               /* dwdot[HO2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[8];
    J[123] -= dqdci;              /* dwdot[O]/d[OH] */
    J[124] += dqdci;              /* dwdot[OH]/d[OH] */
    J[126] -= dqdci;              /* dwdot[H2O2]/d[OH] */
    J[128] += dqdci;              /* dwdot[HO2]/d[OH] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[3];
    J[183] -= dqdci;              /* dwdot[O]/d[H2O2] */
    J[184] += dqdci;              /* dwdot[OH]/d[H2O2] */
    J[186] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    J[188] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[4];
    J[243] -= dqdci;              /* dwdot[O]/d[HO2] */
    J[244] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[246] -= dqdci;              /* dwdot[H2O2]/d[HO2] */
    J[248] += dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[873] -= dqdT;               /* dwdot[O]/dT */
    J[874] += dqdT;               /* dwdot[OH]/dT */
    J[876] -= dqdT;               /* dwdot[H2O2]/dT */
    J[878] += dqdT;               /* dwdot[HO2]/dT */

    /*reaction 34: H2O2 + H <=> H2O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[6];
    k_f = prefactor_units[33] * fwd_A[33]
                * exp(fwd_beta[33] * tc[0] - activation_units[33] * fwd_Ea[33] * invT);
    dlnkfdT = fwd_beta[33] * invT + activation_units[33] * fwd_Ea[33] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[7];
    Kc = exp(g_RT[1] - g_RT[4] + g_RT[6] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[6]) + (h_RT[4] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[4] += q; /* OH */
    wdot[6] -= q; /* H2O2 */
    wdot[7] += q; /* H2O */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[31] -= dqdci;               /* dwdot[H]/d[H] */
    J[34] += dqdci;               /* dwdot[OH]/d[H] */
    J[36] -= dqdci;               /* dwdot[H2O2]/d[H] */
    J[37] += dqdci;               /* dwdot[H2O]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[7];
    J[121] -= dqdci;              /* dwdot[H]/d[OH] */
    J[124] += dqdci;              /* dwdot[OH]/d[OH] */
    J[126] -= dqdci;              /* dwdot[H2O2]/d[OH] */
    J[127] += dqdci;              /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[1];
    J[181] -= dqdci;              /* dwdot[H]/d[H2O2] */
    J[184] += dqdci;              /* dwdot[OH]/d[H2O2] */
    J[186] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    J[187] += dqdci;              /* dwdot[H2O]/d[H2O2] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[4];
    J[211] -= dqdci;              /* dwdot[H]/d[H2O] */
    J[214] += dqdci;              /* dwdot[OH]/d[H2O] */
    J[216] -= dqdci;              /* dwdot[H2O2]/d[H2O] */
    J[217] += dqdci;              /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[871] -= dqdT;               /* dwdot[H]/dT */
    J[874] += dqdT;               /* dwdot[OH]/dT */
    J[876] -= dqdT;               /* dwdot[H2O2]/dT */
    J[877] += dqdT;               /* dwdot[H2O]/dT */

    /*reaction 35: H2O2 + H <=> HO2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[6];
    k_f = prefactor_units[34] * fwd_A[34]
                * exp(fwd_beta[34] * tc[0] - activation_units[34] * fwd_Ea[34] * invT);
    dlnkfdT = fwd_beta[34] * invT + activation_units[34] * fwd_Ea[34] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[8];
    Kc = exp(g_RT[1] - g_RT[2] + g_RT[6] - g_RT[8]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[6]) + (h_RT[2] + h_RT[8]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[2] += q; /* H2 */
    wdot[6] -= q; /* H2O2 */
    wdot[8] += q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[31] -= dqdci;               /* dwdot[H]/d[H] */
    J[32] += dqdci;               /* dwdot[H2]/d[H] */
    J[36] -= dqdci;               /* dwdot[H2O2]/d[H] */
    J[38] += dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[8];
    J[61] -= dqdci;               /* dwdot[H]/d[H2] */
    J[62] += dqdci;               /* dwdot[H2]/d[H2] */
    J[66] -= dqdci;               /* dwdot[H2O2]/d[H2] */
    J[68] += dqdci;               /* dwdot[HO2]/d[H2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[1];
    J[181] -= dqdci;              /* dwdot[H]/d[H2O2] */
    J[182] += dqdci;              /* dwdot[H2]/d[H2O2] */
    J[186] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    J[188] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[2];
    J[241] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[242] += dqdci;              /* dwdot[H2]/d[HO2] */
    J[246] -= dqdci;              /* dwdot[H2O2]/d[HO2] */
    J[248] += dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[871] -= dqdT;               /* dwdot[H]/dT */
    J[872] += dqdT;               /* dwdot[H2]/dT */
    J[876] -= dqdT;               /* dwdot[H2O2]/dT */
    J[878] += dqdT;               /* dwdot[HO2]/dT */

    /*reaction 36: H2O2 + OH <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[6];
    k_f = prefactor_units[35] * fwd_A[35]
                * exp(fwd_beta[35] * tc[0] - activation_units[35] * fwd_Ea[35] * invT);
    dlnkfdT = fwd_beta[35] * invT + activation_units[35] * fwd_Ea[35] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[8];
    Kc = exp(g_RT[4] + g_RT[6] - g_RT[7] - g_RT[8]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[6]) + (h_RT[7] + h_RT[8]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[6] -= q; /* H2O2 */
    wdot[7] += q; /* H2O */
    wdot[8] += q; /* HO2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[6];
    J[124] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[126] -= dqdci;              /* dwdot[H2O2]/d[OH] */
    J[127] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[128] += dqdci;              /* dwdot[HO2]/d[OH] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[4];
    J[184] -= dqdci;              /* dwdot[OH]/d[H2O2] */
    J[186] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    J[187] += dqdci;              /* dwdot[H2O]/d[H2O2] */
    J[188] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[8];
    J[214] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[216] -= dqdci;              /* dwdot[H2O2]/d[H2O] */
    J[217] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[218] += dqdci;              /* dwdot[HO2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[7];
    J[244] -= dqdci;              /* dwdot[OH]/d[HO2] */
    J[246] -= dqdci;              /* dwdot[H2O2]/d[HO2] */
    J[247] += dqdci;              /* dwdot[H2O]/d[HO2] */
    J[248] += dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[874] -= dqdT;               /* dwdot[OH]/dT */
    J[876] -= dqdT;               /* dwdot[H2O2]/dT */
    J[877] += dqdT;               /* dwdot[H2O]/dT */
    J[878] += dqdT;               /* dwdot[HO2]/dT */

    /*reaction 37: H2O2 + OH <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[6];
    k_f = prefactor_units[36] * fwd_A[36]
                * exp(fwd_beta[36] * tc[0] - activation_units[36] * fwd_Ea[36] * invT);
    dlnkfdT = fwd_beta[36] * invT + activation_units[36] * fwd_Ea[36] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[8];
    Kc = exp(g_RT[4] + g_RT[6] - g_RT[7] - g_RT[8]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[6]) + (h_RT[7] + h_RT[8]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[6] -= q; /* H2O2 */
    wdot[7] += q; /* H2O */
    wdot[8] += q; /* HO2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[6];
    J[124] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[126] -= dqdci;              /* dwdot[H2O2]/d[OH] */
    J[127] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[128] += dqdci;              /* dwdot[HO2]/d[OH] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[4];
    J[184] -= dqdci;              /* dwdot[OH]/d[H2O2] */
    J[186] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    J[187] += dqdci;              /* dwdot[H2O]/d[H2O2] */
    J[188] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[8];
    J[214] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[216] -= dqdci;              /* dwdot[H2O2]/d[H2O] */
    J[217] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[218] += dqdci;              /* dwdot[HO2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[7];
    J[244] -= dqdci;              /* dwdot[OH]/d[HO2] */
    J[246] -= dqdci;              /* dwdot[H2O2]/d[HO2] */
    J[247] += dqdci;              /* dwdot[H2O]/d[HO2] */
    J[248] += dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[874] -= dqdT;               /* dwdot[OH]/dT */
    J[876] -= dqdT;               /* dwdot[H2O2]/dT */
    J[877] += dqdT;               /* dwdot[H2O]/dT */
    J[878] += dqdT;               /* dwdot[HO2]/dT */

    /*reaction 38: C + O2 <=> CO + O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[18];
    k_f = prefactor_units[37] * fwd_A[37]
                * exp(fwd_beta[37] * tc[0] - activation_units[37] * fwd_Ea[37] * invT);
    dlnkfdT = fwd_beta[37] * invT + activation_units[37] * fwd_Ea[37] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[9];
    Kc = exp(-g_RT[3] + g_RT[5] - g_RT[9] + g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[18]) + (h_RT[3] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += q; /* O */
    wdot[5] -= q; /* O2 */
    wdot[9] += q; /* CO */
    wdot[18] -= q; /* C */
    /* d()/d[O] */
    dqdci =  - k_r*sc[9];
    J[93] += dqdci;               /* dwdot[O]/d[O] */
    J[95] -= dqdci;               /* dwdot[O2]/d[O] */
    J[99] += dqdci;               /* dwdot[CO]/d[O] */
    J[108] -= dqdci;              /* dwdot[C]/d[O] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[18];
    J[153] += dqdci;              /* dwdot[O]/d[O2] */
    J[155] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[159] += dqdci;              /* dwdot[CO]/d[O2] */
    J[168] -= dqdci;              /* dwdot[C]/d[O2] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[3];
    J[273] += dqdci;              /* dwdot[O]/d[CO] */
    J[275] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[279] += dqdci;              /* dwdot[CO]/d[CO] */
    J[288] -= dqdci;              /* dwdot[C]/d[CO] */
    /* d()/d[C] */
    dqdci =  + k_f*sc[5];
    J[543] += dqdci;              /* dwdot[O]/d[C] */
    J[545] -= dqdci;              /* dwdot[O2]/d[C] */
    J[549] += dqdci;              /* dwdot[CO]/d[C] */
    J[558] -= dqdci;              /* dwdot[C]/d[C] */
    /* d()/dT */
    J[873] += dqdT;               /* dwdot[O]/dT */
    J[875] -= dqdT;               /* dwdot[O2]/dT */
    J[879] += dqdT;               /* dwdot[CO]/dT */
    J[888] -= dqdT;               /* dwdot[C]/dT */

    /*reaction 39: C + OH <=> CO + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[18];
    k_f = prefactor_units[38] * fwd_A[38]
                * exp(fwd_beta[38] * tc[0] - activation_units[38] * fwd_Ea[38] * invT);
    dlnkfdT = fwd_beta[38] * invT + activation_units[38] * fwd_Ea[38] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[9];
    Kc = exp(-g_RT[1] + g_RT[4] - g_RT[9] + g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[18]) + (h_RT[1] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[4] -= q; /* OH */
    wdot[9] += q; /* CO */
    wdot[18] -= q; /* C */
    /* d()/d[H] */
    dqdci =  - k_r*sc[9];
    J[31] += dqdci;               /* dwdot[H]/d[H] */
    J[34] -= dqdci;               /* dwdot[OH]/d[H] */
    J[39] += dqdci;               /* dwdot[CO]/d[H] */
    J[48] -= dqdci;               /* dwdot[C]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[18];
    J[121] += dqdci;              /* dwdot[H]/d[OH] */
    J[124] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[129] += dqdci;              /* dwdot[CO]/d[OH] */
    J[138] -= dqdci;              /* dwdot[C]/d[OH] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[1];
    J[271] += dqdci;              /* dwdot[H]/d[CO] */
    J[274] -= dqdci;              /* dwdot[OH]/d[CO] */
    J[279] += dqdci;              /* dwdot[CO]/d[CO] */
    J[288] -= dqdci;              /* dwdot[C]/d[CO] */
    /* d()/d[C] */
    dqdci =  + k_f*sc[4];
    J[541] += dqdci;              /* dwdot[H]/d[C] */
    J[544] -= dqdci;              /* dwdot[OH]/d[C] */
    J[549] += dqdci;              /* dwdot[CO]/d[C] */
    J[558] -= dqdci;              /* dwdot[C]/d[C] */
    /* d()/dT */
    J[871] += dqdT;               /* dwdot[H]/dT */
    J[874] -= dqdT;               /* dwdot[OH]/dT */
    J[879] += dqdT;               /* dwdot[CO]/dT */
    J[888] -= dqdT;               /* dwdot[C]/dT */

    /*reaction 40: CH + OH <=> HCO + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[19];
    k_f = prefactor_units[39] * fwd_A[39]
                * exp(fwd_beta[39] * tc[0] - activation_units[39] * fwd_Ea[39] * invT);
    dlnkfdT = fwd_beta[39] * invT + activation_units[39] * fwd_Ea[39] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[20];
    Kc = exp(-g_RT[1] + g_RT[4] + g_RT[19] - g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[19]) + (h_RT[1] + h_RT[20]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[4] -= q; /* OH */
    wdot[19] -= q; /* CH */
    wdot[20] += q; /* HCO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[20];
    J[31] += dqdci;               /* dwdot[H]/d[H] */
    J[34] -= dqdci;               /* dwdot[OH]/d[H] */
    J[49] -= dqdci;               /* dwdot[CH]/d[H] */
    J[50] += dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[19];
    J[121] += dqdci;              /* dwdot[H]/d[OH] */
    J[124] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[139] -= dqdci;              /* dwdot[CH]/d[OH] */
    J[140] += dqdci;              /* dwdot[HCO]/d[OH] */
    /* d()/d[CH] */
    dqdci =  + k_f*sc[4];
    J[571] += dqdci;              /* dwdot[H]/d[CH] */
    J[574] -= dqdci;              /* dwdot[OH]/d[CH] */
    J[589] -= dqdci;              /* dwdot[CH]/d[CH] */
    J[590] += dqdci;              /* dwdot[HCO]/d[CH] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[1];
    J[601] += dqdci;              /* dwdot[H]/d[HCO] */
    J[604] -= dqdci;              /* dwdot[OH]/d[HCO] */
    J[619] -= dqdci;              /* dwdot[CH]/d[HCO] */
    J[620] += dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[871] += dqdT;               /* dwdot[H]/dT */
    J[874] -= dqdT;               /* dwdot[OH]/dT */
    J[889] -= dqdT;               /* dwdot[CH]/dT */
    J[890] += dqdT;               /* dwdot[HCO]/dT */

    /*reaction 41: CH + H2 <=> TXCH2 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[19];
    k_f = prefactor_units[40] * fwd_A[40]
                * exp(fwd_beta[40] * tc[0] - activation_units[40] * fwd_Ea[40] * invT);
    dlnkfdT = fwd_beta[40] * invT + activation_units[40] * fwd_Ea[40] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[21];
    Kc = exp(-g_RT[1] + g_RT[2] + g_RT[19] - g_RT[21]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[19]) + (h_RT[1] + h_RT[21]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[2] -= q; /* H2 */
    wdot[19] -= q; /* CH */
    wdot[21] += q; /* TXCH2 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[21];
    J[31] += dqdci;               /* dwdot[H]/d[H] */
    J[32] -= dqdci;               /* dwdot[H2]/d[H] */
    J[49] -= dqdci;               /* dwdot[CH]/d[H] */
    J[51] += dqdci;               /* dwdot[TXCH2]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[19];
    J[61] += dqdci;               /* dwdot[H]/d[H2] */
    J[62] -= dqdci;               /* dwdot[H2]/d[H2] */
    J[79] -= dqdci;               /* dwdot[CH]/d[H2] */
    J[81] += dqdci;               /* dwdot[TXCH2]/d[H2] */
    /* d()/d[CH] */
    dqdci =  + k_f*sc[2];
    J[571] += dqdci;              /* dwdot[H]/d[CH] */
    J[572] -= dqdci;              /* dwdot[H2]/d[CH] */
    J[589] -= dqdci;              /* dwdot[CH]/d[CH] */
    J[591] += dqdci;              /* dwdot[TXCH2]/d[CH] */
    /* d()/d[TXCH2] */
    dqdci =  - k_r*sc[1];
    J[631] += dqdci;              /* dwdot[H]/d[TXCH2] */
    J[632] -= dqdci;              /* dwdot[H2]/d[TXCH2] */
    J[649] -= dqdci;              /* dwdot[CH]/d[TXCH2] */
    J[651] += dqdci;              /* dwdot[TXCH2]/d[TXCH2] */
    /* d()/dT */
    J[871] += dqdT;               /* dwdot[H]/dT */
    J[872] -= dqdT;               /* dwdot[H2]/dT */
    J[889] -= dqdT;               /* dwdot[CH]/dT */
    J[891] += dqdT;               /* dwdot[TXCH2]/dT */

    /*reaction 42: CH + O <=> CO + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[19];
    k_f = prefactor_units[41] * fwd_A[41]
                * exp(fwd_beta[41] * tc[0] - activation_units[41] * fwd_Ea[41] * invT);
    dlnkfdT = fwd_beta[41] * invT + activation_units[41] * fwd_Ea[41] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[9];
    Kc = exp(-g_RT[1] + g_RT[3] - g_RT[9] + g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[19]) + (h_RT[1] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[3] -= q; /* O */
    wdot[9] += q; /* CO */
    wdot[19] -= q; /* CH */
    /* d()/d[H] */
    dqdci =  - k_r*sc[9];
    J[31] += dqdci;               /* dwdot[H]/d[H] */
    J[33] -= dqdci;               /* dwdot[O]/d[H] */
    J[39] += dqdci;               /* dwdot[CO]/d[H] */
    J[49] -= dqdci;               /* dwdot[CH]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[19];
    J[91] += dqdci;               /* dwdot[H]/d[O] */
    J[93] -= dqdci;               /* dwdot[O]/d[O] */
    J[99] += dqdci;               /* dwdot[CO]/d[O] */
    J[109] -= dqdci;              /* dwdot[CH]/d[O] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[1];
    J[271] += dqdci;              /* dwdot[H]/d[CO] */
    J[273] -= dqdci;              /* dwdot[O]/d[CO] */
    J[279] += dqdci;              /* dwdot[CO]/d[CO] */
    J[289] -= dqdci;              /* dwdot[CH]/d[CO] */
    /* d()/d[CH] */
    dqdci =  + k_f*sc[3];
    J[571] += dqdci;              /* dwdot[H]/d[CH] */
    J[573] -= dqdci;              /* dwdot[O]/d[CH] */
    J[579] += dqdci;              /* dwdot[CO]/d[CH] */
    J[589] -= dqdci;              /* dwdot[CH]/d[CH] */
    /* d()/dT */
    J[871] += dqdT;               /* dwdot[H]/dT */
    J[873] -= dqdT;               /* dwdot[O]/dT */
    J[879] += dqdT;               /* dwdot[CO]/dT */
    J[889] -= dqdT;               /* dwdot[CH]/dT */

    /*reaction 43: CH + O2 <=> HCO + O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[19];
    k_f = prefactor_units[42] * fwd_A[42]
                * exp(fwd_beta[42] * tc[0] - activation_units[42] * fwd_Ea[42] * invT);
    dlnkfdT = fwd_beta[42] * invT + activation_units[42] * fwd_Ea[42] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[20];
    Kc = exp(-g_RT[3] + g_RT[5] + g_RT[19] - g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[19]) + (h_RT[3] + h_RT[20]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += q; /* O */
    wdot[5] -= q; /* O2 */
    wdot[19] -= q; /* CH */
    wdot[20] += q; /* HCO */
    /* d()/d[O] */
    dqdci =  - k_r*sc[20];
    J[93] += dqdci;               /* dwdot[O]/d[O] */
    J[95] -= dqdci;               /* dwdot[O2]/d[O] */
    J[109] -= dqdci;              /* dwdot[CH]/d[O] */
    J[110] += dqdci;              /* dwdot[HCO]/d[O] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[19];
    J[153] += dqdci;              /* dwdot[O]/d[O2] */
    J[155] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[169] -= dqdci;              /* dwdot[CH]/d[O2] */
    J[170] += dqdci;              /* dwdot[HCO]/d[O2] */
    /* d()/d[CH] */
    dqdci =  + k_f*sc[5];
    J[573] += dqdci;              /* dwdot[O]/d[CH] */
    J[575] -= dqdci;              /* dwdot[O2]/d[CH] */
    J[589] -= dqdci;              /* dwdot[CH]/d[CH] */
    J[590] += dqdci;              /* dwdot[HCO]/d[CH] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[3];
    J[603] += dqdci;              /* dwdot[O]/d[HCO] */
    J[605] -= dqdci;              /* dwdot[O2]/d[HCO] */
    J[619] -= dqdci;              /* dwdot[CH]/d[HCO] */
    J[620] += dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[873] += dqdT;               /* dwdot[O]/dT */
    J[875] -= dqdT;               /* dwdot[O2]/dT */
    J[889] -= dqdT;               /* dwdot[CH]/dT */
    J[890] += dqdT;               /* dwdot[HCO]/dT */

    /*reaction 44: CH + H <=> C + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[19];
    k_f = prefactor_units[43] * fwd_A[43]
                * exp(fwd_beta[43] * tc[0] - activation_units[43] * fwd_Ea[43] * invT);
    dlnkfdT = fwd_beta[43] * invT + activation_units[43] * fwd_Ea[43] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[18];
    Kc = exp(g_RT[1] - g_RT[2] - g_RT[18] + g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[19]) + (h_RT[2] + h_RT[18]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[2] += q; /* H2 */
    wdot[18] += q; /* C */
    wdot[19] -= q; /* CH */
    /* d()/d[H] */
    dqdci =  + k_f*sc[19];
    J[31] -= dqdci;               /* dwdot[H]/d[H] */
    J[32] += dqdci;               /* dwdot[H2]/d[H] */
    J[48] += dqdci;               /* dwdot[C]/d[H] */
    J[49] -= dqdci;               /* dwdot[CH]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[18];
    J[61] -= dqdci;               /* dwdot[H]/d[H2] */
    J[62] += dqdci;               /* dwdot[H2]/d[H2] */
    J[78] += dqdci;               /* dwdot[C]/d[H2] */
    J[79] -= dqdci;               /* dwdot[CH]/d[H2] */
    /* d()/d[C] */
    dqdci =  - k_r*sc[2];
    J[541] -= dqdci;              /* dwdot[H]/d[C] */
    J[542] += dqdci;              /* dwdot[H2]/d[C] */
    J[558] += dqdci;              /* dwdot[C]/d[C] */
    J[559] -= dqdci;              /* dwdot[CH]/d[C] */
    /* d()/d[CH] */
    dqdci =  + k_f*sc[1];
    J[571] -= dqdci;              /* dwdot[H]/d[CH] */
    J[572] += dqdci;              /* dwdot[H2]/d[CH] */
    J[588] += dqdci;              /* dwdot[C]/d[CH] */
    J[589] -= dqdci;              /* dwdot[CH]/d[CH] */
    /* d()/dT */
    J[871] -= dqdT;               /* dwdot[H]/dT */
    J[872] += dqdT;               /* dwdot[H2]/dT */
    J[888] += dqdT;               /* dwdot[C]/dT */
    J[889] -= dqdT;               /* dwdot[CH]/dT */

    /*reaction 45: CH + H2O <=> CH2O + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[7]*sc[19];
    k_f = prefactor_units[44] * fwd_A[44]
                * exp(fwd_beta[44] * tc[0] - activation_units[44] * fwd_Ea[44] * invT);
    dlnkfdT = fwd_beta[44] * invT + activation_units[44] * fwd_Ea[44] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[11];
    Kc = exp(-g_RT[1] + g_RT[7] - g_RT[11] + g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7] + h_RT[19]) + (h_RT[1] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[7] -= q; /* H2O */
    wdot[11] += q; /* CH2O */
    wdot[19] -= q; /* CH */
    /* d()/d[H] */
    dqdci =  - k_r*sc[11];
    J[31] += dqdci;               /* dwdot[H]/d[H] */
    J[37] -= dqdci;               /* dwdot[H2O]/d[H] */
    J[41] += dqdci;               /* dwdot[CH2O]/d[H] */
    J[49] -= dqdci;               /* dwdot[CH]/d[H] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[19];
    J[211] += dqdci;              /* dwdot[H]/d[H2O] */
    J[217] -= dqdci;              /* dwdot[H2O]/d[H2O] */
    J[221] += dqdci;              /* dwdot[CH2O]/d[H2O] */
    J[229] -= dqdci;              /* dwdot[CH]/d[H2O] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[1];
    J[331] += dqdci;              /* dwdot[H]/d[CH2O] */
    J[337] -= dqdci;              /* dwdot[H2O]/d[CH2O] */
    J[341] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[349] -= dqdci;              /* dwdot[CH]/d[CH2O] */
    /* d()/d[CH] */
    dqdci =  + k_f*sc[7];
    J[571] += dqdci;              /* dwdot[H]/d[CH] */
    J[577] -= dqdci;              /* dwdot[H2O]/d[CH] */
    J[581] += dqdci;              /* dwdot[CH2O]/d[CH] */
    J[589] -= dqdci;              /* dwdot[CH]/d[CH] */
    /* d()/dT */
    J[871] += dqdT;               /* dwdot[H]/dT */
    J[877] -= dqdT;               /* dwdot[H2O]/dT */
    J[881] += dqdT;               /* dwdot[CH2O]/dT */
    J[889] -= dqdT;               /* dwdot[CH]/dT */

    /*reaction 46: TXCH2 + O2 => OH + H + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[21];
    k_f = prefactor_units[45] * fwd_A[45]
                * exp(fwd_beta[45] * tc[0] - activation_units[45] * fwd_Ea[45] * invT);
    dlnkfdT = fwd_beta[45] * invT + activation_units[45] * fwd_Ea[45] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[4] += q; /* OH */
    wdot[5] -= q; /* O2 */
    wdot[9] += q; /* CO */
    wdot[21] -= q; /* TXCH2 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[21];
    J[151] += dqdci;              /* dwdot[H]/d[O2] */
    J[154] += dqdci;              /* dwdot[OH]/d[O2] */
    J[155] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[159] += dqdci;              /* dwdot[CO]/d[O2] */
    J[171] -= dqdci;              /* dwdot[TXCH2]/d[O2] */
    /* d()/d[TXCH2] */
    dqdci =  + k_f*sc[5];
    J[631] += dqdci;              /* dwdot[H]/d[TXCH2] */
    J[634] += dqdci;              /* dwdot[OH]/d[TXCH2] */
    J[635] -= dqdci;              /* dwdot[O2]/d[TXCH2] */
    J[639] += dqdci;              /* dwdot[CO]/d[TXCH2] */
    J[651] -= dqdci;              /* dwdot[TXCH2]/d[TXCH2] */
    /* d()/dT */
    J[871] += dqdT;               /* dwdot[H]/dT */
    J[874] += dqdT;               /* dwdot[OH]/dT */
    J[875] -= dqdT;               /* dwdot[O2]/dT */
    J[879] += dqdT;               /* dwdot[CO]/dT */
    J[891] -= dqdT;               /* dwdot[TXCH2]/dT */

    /*reaction 47: TXCH2 + O2 <=> CH2O + O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[21];
    k_f = prefactor_units[46] * fwd_A[46]
                * exp(fwd_beta[46] * tc[0] - activation_units[46] * fwd_Ea[46] * invT);
    dlnkfdT = fwd_beta[46] * invT + activation_units[46] * fwd_Ea[46] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[11];
    Kc = exp(-g_RT[3] + g_RT[5] - g_RT[11] + g_RT[21]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[21]) + (h_RT[3] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += q; /* O */
    wdot[5] -= q; /* O2 */
    wdot[11] += q; /* CH2O */
    wdot[21] -= q; /* TXCH2 */
    /* d()/d[O] */
    dqdci =  - k_r*sc[11];
    J[93] += dqdci;               /* dwdot[O]/d[O] */
    J[95] -= dqdci;               /* dwdot[O2]/d[O] */
    J[101] += dqdci;              /* dwdot[CH2O]/d[O] */
    J[111] -= dqdci;              /* dwdot[TXCH2]/d[O] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[21];
    J[153] += dqdci;              /* dwdot[O]/d[O2] */
    J[155] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[161] += dqdci;              /* dwdot[CH2O]/d[O2] */
    J[171] -= dqdci;              /* dwdot[TXCH2]/d[O2] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[3];
    J[333] += dqdci;              /* dwdot[O]/d[CH2O] */
    J[335] -= dqdci;              /* dwdot[O2]/d[CH2O] */
    J[341] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[351] -= dqdci;              /* dwdot[TXCH2]/d[CH2O] */
    /* d()/d[TXCH2] */
    dqdci =  + k_f*sc[5];
    J[633] += dqdci;              /* dwdot[O]/d[TXCH2] */
    J[635] -= dqdci;              /* dwdot[O2]/d[TXCH2] */
    J[641] += dqdci;              /* dwdot[CH2O]/d[TXCH2] */
    J[651] -= dqdci;              /* dwdot[TXCH2]/d[TXCH2] */
    /* d()/dT */
    J[873] += dqdT;               /* dwdot[O]/dT */
    J[875] -= dqdT;               /* dwdot[O2]/dT */
    J[881] += dqdT;               /* dwdot[CH2O]/dT */
    J[891] -= dqdT;               /* dwdot[TXCH2]/dT */

    /*reaction 48: TXCH2 + OH <=> CH2O + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[21];
    k_f = prefactor_units[47] * fwd_A[47]
                * exp(fwd_beta[47] * tc[0] - activation_units[47] * fwd_Ea[47] * invT);
    dlnkfdT = fwd_beta[47] * invT + activation_units[47] * fwd_Ea[47] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[11];
    Kc = exp(-g_RT[1] + g_RT[4] - g_RT[11] + g_RT[21]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[21]) + (h_RT[1] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[4] -= q; /* OH */
    wdot[11] += q; /* CH2O */
    wdot[21] -= q; /* TXCH2 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[11];
    J[31] += dqdci;               /* dwdot[H]/d[H] */
    J[34] -= dqdci;               /* dwdot[OH]/d[H] */
    J[41] += dqdci;               /* dwdot[CH2O]/d[H] */
    J[51] -= dqdci;               /* dwdot[TXCH2]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[21];
    J[121] += dqdci;              /* dwdot[H]/d[OH] */
    J[124] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[131] += dqdci;              /* dwdot[CH2O]/d[OH] */
    J[141] -= dqdci;              /* dwdot[TXCH2]/d[OH] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[1];
    J[331] += dqdci;              /* dwdot[H]/d[CH2O] */
    J[334] -= dqdci;              /* dwdot[OH]/d[CH2O] */
    J[341] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[351] -= dqdci;              /* dwdot[TXCH2]/d[CH2O] */
    /* d()/d[TXCH2] */
    dqdci =  + k_f*sc[4];
    J[631] += dqdci;              /* dwdot[H]/d[TXCH2] */
    J[634] -= dqdci;              /* dwdot[OH]/d[TXCH2] */
    J[641] += dqdci;              /* dwdot[CH2O]/d[TXCH2] */
    J[651] -= dqdci;              /* dwdot[TXCH2]/d[TXCH2] */
    /* d()/dT */
    J[871] += dqdT;               /* dwdot[H]/dT */
    J[874] -= dqdT;               /* dwdot[OH]/dT */
    J[881] += dqdT;               /* dwdot[CH2O]/dT */
    J[891] -= dqdT;               /* dwdot[TXCH2]/dT */

    /*reaction 49: TXCH2 + HO2 <=> CH2O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[21];
    k_f = prefactor_units[48] * fwd_A[48]
                * exp(fwd_beta[48] * tc[0] - activation_units[48] * fwd_Ea[48] * invT);
    dlnkfdT = fwd_beta[48] * invT + activation_units[48] * fwd_Ea[48] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[11];
    Kc = exp(-g_RT[4] + g_RT[8] - g_RT[11] + g_RT[21]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[21]) + (h_RT[4] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] += q; /* OH */
    wdot[8] -= q; /* HO2 */
    wdot[11] += q; /* CH2O */
    wdot[21] -= q; /* TXCH2 */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[11];
    J[124] += dqdci;              /* dwdot[OH]/d[OH] */
    J[128] -= dqdci;              /* dwdot[HO2]/d[OH] */
    J[131] += dqdci;              /* dwdot[CH2O]/d[OH] */
    J[141] -= dqdci;              /* dwdot[TXCH2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[21];
    J[244] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[248] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[251] += dqdci;              /* dwdot[CH2O]/d[HO2] */
    J[261] -= dqdci;              /* dwdot[TXCH2]/d[HO2] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[4];
    J[334] += dqdci;              /* dwdot[OH]/d[CH2O] */
    J[338] -= dqdci;              /* dwdot[HO2]/d[CH2O] */
    J[341] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[351] -= dqdci;              /* dwdot[TXCH2]/d[CH2O] */
    /* d()/d[TXCH2] */
    dqdci =  + k_f*sc[8];
    J[634] += dqdci;              /* dwdot[OH]/d[TXCH2] */
    J[638] -= dqdci;              /* dwdot[HO2]/d[TXCH2] */
    J[641] += dqdci;              /* dwdot[CH2O]/d[TXCH2] */
    J[651] -= dqdci;              /* dwdot[TXCH2]/d[TXCH2] */
    /* d()/dT */
    J[874] += dqdT;               /* dwdot[OH]/dT */
    J[878] -= dqdT;               /* dwdot[HO2]/dT */
    J[881] += dqdT;               /* dwdot[CH2O]/dT */
    J[891] -= dqdT;               /* dwdot[TXCH2]/dT */

    /*reaction 50: TXCH2 + O2 => CO2 + 2.000000 H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[21];
    k_f = prefactor_units[49] * fwd_A[49]
                * exp(fwd_beta[49] * tc[0] - activation_units[49] * fwd_Ea[49] * invT);
    dlnkfdT = fwd_beta[49] * invT + activation_units[49] * fwd_Ea[49] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] += 2 * q; /* H */
    wdot[5] -= q; /* O2 */
    wdot[12] += q; /* CO2 */
    wdot[21] -= q; /* TXCH2 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[21];
    J[151] += 2 * dqdci;          /* dwdot[H]/d[O2] */
    J[155] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[162] += dqdci;              /* dwdot[CO2]/d[O2] */
    J[171] -= dqdci;              /* dwdot[TXCH2]/d[O2] */
    /* d()/d[TXCH2] */
    dqdci =  + k_f*sc[5];
    J[631] += 2 * dqdci;          /* dwdot[H]/d[TXCH2] */
    J[635] -= dqdci;              /* dwdot[O2]/d[TXCH2] */
    J[642] += dqdci;              /* dwdot[CO2]/d[TXCH2] */
    J[651] -= dqdci;              /* dwdot[TXCH2]/d[TXCH2] */
    /* d()/dT */
    J[871] += 2 * dqdT;           /* dwdot[H]/dT */
    J[875] -= dqdT;               /* dwdot[O2]/dT */
    J[882] += dqdT;               /* dwdot[CO2]/dT */
    J[891] -= dqdT;               /* dwdot[TXCH2]/dT */

    /*reaction 51: TXCH2 + OH <=> CH + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[21];
    k_f = prefactor_units[50] * fwd_A[50]
                * exp(fwd_beta[50] * tc[0] - activation_units[50] * fwd_Ea[50] * invT);
    dlnkfdT = fwd_beta[50] * invT + activation_units[50] * fwd_Ea[50] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[19];
    Kc = exp(g_RT[4] - g_RT[7] - g_RT[19] + g_RT[21]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[21]) + (h_RT[7] + h_RT[19]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[7] += q; /* H2O */
    wdot[19] += q; /* CH */
    wdot[21] -= q; /* TXCH2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[21];
    J[124] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[127] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[139] += dqdci;              /* dwdot[CH]/d[OH] */
    J[141] -= dqdci;              /* dwdot[TXCH2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[19];
    J[214] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[217] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[229] += dqdci;              /* dwdot[CH]/d[H2O] */
    J[231] -= dqdci;              /* dwdot[TXCH2]/d[H2O] */
    /* d()/d[CH] */
    dqdci =  - k_r*sc[7];
    J[574] -= dqdci;              /* dwdot[OH]/d[CH] */
    J[577] += dqdci;              /* dwdot[H2O]/d[CH] */
    J[589] += dqdci;              /* dwdot[CH]/d[CH] */
    J[591] -= dqdci;              /* dwdot[TXCH2]/d[CH] */
    /* d()/d[TXCH2] */
    dqdci =  + k_f*sc[4];
    J[634] -= dqdci;              /* dwdot[OH]/d[TXCH2] */
    J[637] += dqdci;              /* dwdot[H2O]/d[TXCH2] */
    J[649] += dqdci;              /* dwdot[CH]/d[TXCH2] */
    J[651] -= dqdci;              /* dwdot[TXCH2]/d[TXCH2] */
    /* d()/dT */
    J[874] -= dqdT;               /* dwdot[OH]/dT */
    J[877] += dqdT;               /* dwdot[H2O]/dT */
    J[889] += dqdT;               /* dwdot[CH]/dT */
    J[891] -= dqdT;               /* dwdot[TXCH2]/dT */

    /*reaction 52: TXCH2 + O <=> HCO + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[21];
    k_f = prefactor_units[51] * fwd_A[51]
                * exp(fwd_beta[51] * tc[0] - activation_units[51] * fwd_Ea[51] * invT);
    dlnkfdT = fwd_beta[51] * invT + activation_units[51] * fwd_Ea[51] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[20];
    Kc = exp(-g_RT[1] + g_RT[3] - g_RT[20] + g_RT[21]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[21]) + (h_RT[1] + h_RT[20]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[3] -= q; /* O */
    wdot[20] += q; /* HCO */
    wdot[21] -= q; /* TXCH2 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[20];
    J[31] += dqdci;               /* dwdot[H]/d[H] */
    J[33] -= dqdci;               /* dwdot[O]/d[H] */
    J[50] += dqdci;               /* dwdot[HCO]/d[H] */
    J[51] -= dqdci;               /* dwdot[TXCH2]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[21];
    J[91] += dqdci;               /* dwdot[H]/d[O] */
    J[93] -= dqdci;               /* dwdot[O]/d[O] */
    J[110] += dqdci;              /* dwdot[HCO]/d[O] */
    J[111] -= dqdci;              /* dwdot[TXCH2]/d[O] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[1];
    J[601] += dqdci;              /* dwdot[H]/d[HCO] */
    J[603] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[620] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[621] -= dqdci;              /* dwdot[TXCH2]/d[HCO] */
    /* d()/d[TXCH2] */
    dqdci =  + k_f*sc[3];
    J[631] += dqdci;              /* dwdot[H]/d[TXCH2] */
    J[633] -= dqdci;              /* dwdot[O]/d[TXCH2] */
    J[650] += dqdci;              /* dwdot[HCO]/d[TXCH2] */
    J[651] -= dqdci;              /* dwdot[TXCH2]/d[TXCH2] */
    /* d()/dT */
    J[871] += dqdT;               /* dwdot[H]/dT */
    J[873] -= dqdT;               /* dwdot[O]/dT */
    J[890] += dqdT;               /* dwdot[HCO]/dT */
    J[891] -= dqdT;               /* dwdot[TXCH2]/dT */

    /*reaction 53: TXCH2 + H2 <=> H + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[21];
    k_f = prefactor_units[52] * fwd_A[52]
                * exp(fwd_beta[52] * tc[0] - activation_units[52] * fwd_Ea[52] * invT);
    dlnkfdT = fwd_beta[52] * invT + activation_units[52] * fwd_Ea[52] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[10];
    Kc = exp(-g_RT[1] + g_RT[2] - g_RT[10] + g_RT[21]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[21]) + (h_RT[1] + h_RT[10]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[2] -= q; /* H2 */
    wdot[10] += q; /* CH3 */
    wdot[21] -= q; /* TXCH2 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[10];
    J[31] += dqdci;               /* dwdot[H]/d[H] */
    J[32] -= dqdci;               /* dwdot[H2]/d[H] */
    J[40] += dqdci;               /* dwdot[CH3]/d[H] */
    J[51] -= dqdci;               /* dwdot[TXCH2]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[21];
    J[61] += dqdci;               /* dwdot[H]/d[H2] */
    J[62] -= dqdci;               /* dwdot[H2]/d[H2] */
    J[70] += dqdci;               /* dwdot[CH3]/d[H2] */
    J[81] -= dqdci;               /* dwdot[TXCH2]/d[H2] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[1];
    J[301] += dqdci;              /* dwdot[H]/d[CH3] */
    J[302] -= dqdci;              /* dwdot[H2]/d[CH3] */
    J[310] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[321] -= dqdci;              /* dwdot[TXCH2]/d[CH3] */
    /* d()/d[TXCH2] */
    dqdci =  + k_f*sc[2];
    J[631] += dqdci;              /* dwdot[H]/d[TXCH2] */
    J[632] -= dqdci;              /* dwdot[H2]/d[TXCH2] */
    J[640] += dqdci;              /* dwdot[CH3]/d[TXCH2] */
    J[651] -= dqdci;              /* dwdot[TXCH2]/d[TXCH2] */
    /* d()/dT */
    J[871] += dqdT;               /* dwdot[H]/dT */
    J[872] -= dqdT;               /* dwdot[H2]/dT */
    J[880] += dqdT;               /* dwdot[CH3]/dT */
    J[891] -= dqdT;               /* dwdot[TXCH2]/dT */

    /*reaction 54: SXCH2 + H2O <=> TXCH2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[7]*sc[22];
    k_f = prefactor_units[53] * fwd_A[53]
                * exp(fwd_beta[53] * tc[0] - activation_units[53] * fwd_Ea[53] * invT);
    dlnkfdT = fwd_beta[53] * invT + activation_units[53] * fwd_Ea[53] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[21];
    Kc = exp(g_RT[7] - g_RT[7] - g_RT[21] + g_RT[22]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7] + h_RT[22]) + (h_RT[7] + h_RT[21]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[21] += q; /* TXCH2 */
    wdot[22] -= q; /* SXCH2 */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[22] - k_r*sc[21];
    J[231] += dqdci;              /* dwdot[TXCH2]/d[H2O] */
    J[232] -= dqdci;              /* dwdot[SXCH2]/d[H2O] */
    /* d()/d[TXCH2] */
    dqdci =  - k_r*sc[7];
    J[651] += dqdci;              /* dwdot[TXCH2]/d[TXCH2] */
    J[652] -= dqdci;              /* dwdot[SXCH2]/d[TXCH2] */
    /* d()/d[SXCH2] */
    dqdci =  + k_f*sc[7];
    J[681] += dqdci;              /* dwdot[TXCH2]/d[SXCH2] */
    J[682] -= dqdci;              /* dwdot[SXCH2]/d[SXCH2] */
    /* d()/dT */
    J[891] += dqdT;               /* dwdot[TXCH2]/dT */
    J[892] -= dqdT;               /* dwdot[SXCH2]/dT */

    /*reaction 55: SXCH2 + H <=> CH + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[22];
    k_f = prefactor_units[54] * fwd_A[54]
                * exp(fwd_beta[54] * tc[0] - activation_units[54] * fwd_Ea[54] * invT);
    dlnkfdT = fwd_beta[54] * invT + activation_units[54] * fwd_Ea[54] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[19];
    Kc = exp(g_RT[1] - g_RT[2] - g_RT[19] + g_RT[22]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[22]) + (h_RT[2] + h_RT[19]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[2] += q; /* H2 */
    wdot[19] += q; /* CH */
    wdot[22] -= q; /* SXCH2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[22];
    J[31] -= dqdci;               /* dwdot[H]/d[H] */
    J[32] += dqdci;               /* dwdot[H2]/d[H] */
    J[49] += dqdci;               /* dwdot[CH]/d[H] */
    J[52] -= dqdci;               /* dwdot[SXCH2]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[19];
    J[61] -= dqdci;               /* dwdot[H]/d[H2] */
    J[62] += dqdci;               /* dwdot[H2]/d[H2] */
    J[79] += dqdci;               /* dwdot[CH]/d[H2] */
    J[82] -= dqdci;               /* dwdot[SXCH2]/d[H2] */
    /* d()/d[CH] */
    dqdci =  - k_r*sc[2];
    J[571] -= dqdci;              /* dwdot[H]/d[CH] */
    J[572] += dqdci;              /* dwdot[H2]/d[CH] */
    J[589] += dqdci;              /* dwdot[CH]/d[CH] */
    J[592] -= dqdci;              /* dwdot[SXCH2]/d[CH] */
    /* d()/d[SXCH2] */
    dqdci =  + k_f*sc[1];
    J[661] -= dqdci;              /* dwdot[H]/d[SXCH2] */
    J[662] += dqdci;              /* dwdot[H2]/d[SXCH2] */
    J[679] += dqdci;              /* dwdot[CH]/d[SXCH2] */
    J[682] -= dqdci;              /* dwdot[SXCH2]/d[SXCH2] */
    /* d()/dT */
    J[871] -= dqdT;               /* dwdot[H]/dT */
    J[872] += dqdT;               /* dwdot[H2]/dT */
    J[889] += dqdT;               /* dwdot[CH]/dT */
    J[892] -= dqdT;               /* dwdot[SXCH2]/dT */

    /*reaction 56: SXCH2 + O2 <=> H + OH + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[22];
    k_f = prefactor_units[55] * fwd_A[55]
                * exp(fwd_beta[55] * tc[0] - activation_units[55] * fwd_Ea[55] * invT);
    dlnkfdT = fwd_beta[55] * invT + activation_units[55] * fwd_Ea[55] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[4]*sc[9];
    Kc = refC * exp(-g_RT[1] - g_RT[4] + g_RT[5] - g_RT[9] + g_RT[22]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[22]) + (h_RT[1] + h_RT[4] + h_RT[9]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[4] += q; /* OH */
    wdot[5] -= q; /* O2 */
    wdot[9] += q; /* CO */
    wdot[22] -= q; /* SXCH2 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[4]*sc[9];
    J[31] += dqdci;               /* dwdot[H]/d[H] */
    J[34] += dqdci;               /* dwdot[OH]/d[H] */
    J[35] -= dqdci;               /* dwdot[O2]/d[H] */
    J[39] += dqdci;               /* dwdot[CO]/d[H] */
    J[52] -= dqdci;               /* dwdot[SXCH2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[1]*sc[9];
    J[121] += dqdci;              /* dwdot[H]/d[OH] */
    J[124] += dqdci;              /* dwdot[OH]/d[OH] */
    J[125] -= dqdci;              /* dwdot[O2]/d[OH] */
    J[129] += dqdci;              /* dwdot[CO]/d[OH] */
    J[142] -= dqdci;              /* dwdot[SXCH2]/d[OH] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[22];
    J[151] += dqdci;              /* dwdot[H]/d[O2] */
    J[154] += dqdci;              /* dwdot[OH]/d[O2] */
    J[155] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[159] += dqdci;              /* dwdot[CO]/d[O2] */
    J[172] -= dqdci;              /* dwdot[SXCH2]/d[O2] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[1]*sc[4];
    J[271] += dqdci;              /* dwdot[H]/d[CO] */
    J[274] += dqdci;              /* dwdot[OH]/d[CO] */
    J[275] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[279] += dqdci;              /* dwdot[CO]/d[CO] */
    J[292] -= dqdci;              /* dwdot[SXCH2]/d[CO] */
    /* d()/d[SXCH2] */
    dqdci =  + k_f*sc[5];
    J[661] += dqdci;              /* dwdot[H]/d[SXCH2] */
    J[664] += dqdci;              /* dwdot[OH]/d[SXCH2] */
    J[665] -= dqdci;              /* dwdot[O2]/d[SXCH2] */
    J[669] += dqdci;              /* dwdot[CO]/d[SXCH2] */
    J[682] -= dqdci;              /* dwdot[SXCH2]/d[SXCH2] */
    /* d()/dT */
    J[871] += dqdT;               /* dwdot[H]/dT */
    J[874] += dqdT;               /* dwdot[OH]/dT */
    J[875] -= dqdT;               /* dwdot[O2]/dT */
    J[879] += dqdT;               /* dwdot[CO]/dT */
    J[892] -= dqdT;               /* dwdot[SXCH2]/dT */

    /*reaction 57: SXCH2 + O <=> CO + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[22];
    k_f = prefactor_units[56] * fwd_A[56]
                * exp(fwd_beta[56] * tc[0] - activation_units[56] * fwd_Ea[56] * invT);
    dlnkfdT = fwd_beta[56] * invT + activation_units[56] * fwd_Ea[56] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[9];
    Kc = exp(-g_RT[2] + g_RT[3] - g_RT[9] + g_RT[22]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[22]) + (h_RT[2] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* H2 */
    wdot[3] -= q; /* O */
    wdot[9] += q; /* CO */
    wdot[22] -= q; /* SXCH2 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[9];
    J[62] += dqdci;               /* dwdot[H2]/d[H2] */
    J[63] -= dqdci;               /* dwdot[O]/d[H2] */
    J[69] += dqdci;               /* dwdot[CO]/d[H2] */
    J[82] -= dqdci;               /* dwdot[SXCH2]/d[H2] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[22];
    J[92] += dqdci;               /* dwdot[H2]/d[O] */
    J[93] -= dqdci;               /* dwdot[O]/d[O] */
    J[99] += dqdci;               /* dwdot[CO]/d[O] */
    J[112] -= dqdci;              /* dwdot[SXCH2]/d[O] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[2];
    J[272] += dqdci;              /* dwdot[H2]/d[CO] */
    J[273] -= dqdci;              /* dwdot[O]/d[CO] */
    J[279] += dqdci;              /* dwdot[CO]/d[CO] */
    J[292] -= dqdci;              /* dwdot[SXCH2]/d[CO] */
    /* d()/d[SXCH2] */
    dqdci =  + k_f*sc[3];
    J[662] += dqdci;              /* dwdot[H2]/d[SXCH2] */
    J[663] -= dqdci;              /* dwdot[O]/d[SXCH2] */
    J[669] += dqdci;              /* dwdot[CO]/d[SXCH2] */
    J[682] -= dqdci;              /* dwdot[SXCH2]/d[SXCH2] */
    /* d()/dT */
    J[872] += dqdT;               /* dwdot[H2]/dT */
    J[873] -= dqdT;               /* dwdot[O]/dT */
    J[879] += dqdT;               /* dwdot[CO]/dT */
    J[892] -= dqdT;               /* dwdot[SXCH2]/dT */

    /*reaction 58: SXCH2 + O2 <=> CO + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[22];
    k_f = prefactor_units[57] * fwd_A[57]
                * exp(fwd_beta[57] * tc[0] - activation_units[57] * fwd_Ea[57] * invT);
    dlnkfdT = fwd_beta[57] * invT + activation_units[57] * fwd_Ea[57] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[9];
    Kc = exp(g_RT[5] - g_RT[7] - g_RT[9] + g_RT[22]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[22]) + (h_RT[7] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[7] += q; /* H2O */
    wdot[9] += q; /* CO */
    wdot[22] -= q; /* SXCH2 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[22];
    J[155] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[157] += dqdci;              /* dwdot[H2O]/d[O2] */
    J[159] += dqdci;              /* dwdot[CO]/d[O2] */
    J[172] -= dqdci;              /* dwdot[SXCH2]/d[O2] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[9];
    J[215] -= dqdci;              /* dwdot[O2]/d[H2O] */
    J[217] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[219] += dqdci;              /* dwdot[CO]/d[H2O] */
    J[232] -= dqdci;              /* dwdot[SXCH2]/d[H2O] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[7];
    J[275] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[277] += dqdci;              /* dwdot[H2O]/d[CO] */
    J[279] += dqdci;              /* dwdot[CO]/d[CO] */
    J[292] -= dqdci;              /* dwdot[SXCH2]/d[CO] */
    /* d()/d[SXCH2] */
    dqdci =  + k_f*sc[5];
    J[665] -= dqdci;              /* dwdot[O2]/d[SXCH2] */
    J[667] += dqdci;              /* dwdot[H2O]/d[SXCH2] */
    J[669] += dqdci;              /* dwdot[CO]/d[SXCH2] */
    J[682] -= dqdci;              /* dwdot[SXCH2]/d[SXCH2] */
    /* d()/dT */
    J[875] -= dqdT;               /* dwdot[O2]/dT */
    J[877] += dqdT;               /* dwdot[H2O]/dT */
    J[879] += dqdT;               /* dwdot[CO]/dT */
    J[892] -= dqdT;               /* dwdot[SXCH2]/dT */

    /*reaction 59: SXCH2 + H2 <=> CH3 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[22];
    k_f = prefactor_units[58] * fwd_A[58]
                * exp(fwd_beta[58] * tc[0] - activation_units[58] * fwd_Ea[58] * invT);
    dlnkfdT = fwd_beta[58] * invT + activation_units[58] * fwd_Ea[58] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[10];
    Kc = exp(-g_RT[1] + g_RT[2] - g_RT[10] + g_RT[22]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[22]) + (h_RT[1] + h_RT[10]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[2] -= q; /* H2 */
    wdot[10] += q; /* CH3 */
    wdot[22] -= q; /* SXCH2 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[10];
    J[31] += dqdci;               /* dwdot[H]/d[H] */
    J[32] -= dqdci;               /* dwdot[H2]/d[H] */
    J[40] += dqdci;               /* dwdot[CH3]/d[H] */
    J[52] -= dqdci;               /* dwdot[SXCH2]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[22];
    J[61] += dqdci;               /* dwdot[H]/d[H2] */
    J[62] -= dqdci;               /* dwdot[H2]/d[H2] */
    J[70] += dqdci;               /* dwdot[CH3]/d[H2] */
    J[82] -= dqdci;               /* dwdot[SXCH2]/d[H2] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[1];
    J[301] += dqdci;              /* dwdot[H]/d[CH3] */
    J[302] -= dqdci;              /* dwdot[H2]/d[CH3] */
    J[310] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[322] -= dqdci;              /* dwdot[SXCH2]/d[CH3] */
    /* d()/d[SXCH2] */
    dqdci =  + k_f*sc[2];
    J[661] += dqdci;              /* dwdot[H]/d[SXCH2] */
    J[662] -= dqdci;              /* dwdot[H2]/d[SXCH2] */
    J[670] += dqdci;              /* dwdot[CH3]/d[SXCH2] */
    J[682] -= dqdci;              /* dwdot[SXCH2]/d[SXCH2] */
    /* d()/dT */
    J[871] += dqdT;               /* dwdot[H]/dT */
    J[872] -= dqdT;               /* dwdot[H2]/dT */
    J[880] += dqdT;               /* dwdot[CH3]/dT */
    J[892] -= dqdT;               /* dwdot[SXCH2]/dT */

    /*reaction 60: SXCH2 + O <=> HCO + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[22];
    k_f = prefactor_units[59] * fwd_A[59]
                * exp(fwd_beta[59] * tc[0] - activation_units[59] * fwd_Ea[59] * invT);
    dlnkfdT = fwd_beta[59] * invT + activation_units[59] * fwd_Ea[59] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[20];
    Kc = exp(-g_RT[1] + g_RT[3] - g_RT[20] + g_RT[22]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[22]) + (h_RT[1] + h_RT[20]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[3] -= q; /* O */
    wdot[20] += q; /* HCO */
    wdot[22] -= q; /* SXCH2 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[20];
    J[31] += dqdci;               /* dwdot[H]/d[H] */
    J[33] -= dqdci;               /* dwdot[O]/d[H] */
    J[50] += dqdci;               /* dwdot[HCO]/d[H] */
    J[52] -= dqdci;               /* dwdot[SXCH2]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[22];
    J[91] += dqdci;               /* dwdot[H]/d[O] */
    J[93] -= dqdci;               /* dwdot[O]/d[O] */
    J[110] += dqdci;              /* dwdot[HCO]/d[O] */
    J[112] -= dqdci;              /* dwdot[SXCH2]/d[O] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[1];
    J[601] += dqdci;              /* dwdot[H]/d[HCO] */
    J[603] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[620] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[622] -= dqdci;              /* dwdot[SXCH2]/d[HCO] */
    /* d()/d[SXCH2] */
    dqdci =  + k_f*sc[3];
    J[661] += dqdci;              /* dwdot[H]/d[SXCH2] */
    J[663] -= dqdci;              /* dwdot[O]/d[SXCH2] */
    J[680] += dqdci;              /* dwdot[HCO]/d[SXCH2] */
    J[682] -= dqdci;              /* dwdot[SXCH2]/d[SXCH2] */
    /* d()/dT */
    J[871] += dqdT;               /* dwdot[H]/dT */
    J[873] -= dqdT;               /* dwdot[O]/dT */
    J[890] += dqdT;               /* dwdot[HCO]/dT */
    J[892] -= dqdT;               /* dwdot[SXCH2]/dT */

    /*reaction 61: SXCH2 + H2O => H2 + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[7]*sc[22];
    k_f = prefactor_units[60] * fwd_A[60]
                * exp(fwd_beta[60] * tc[0] - activation_units[60] * fwd_Ea[60] * invT);
    dlnkfdT = fwd_beta[60] * invT + activation_units[60] * fwd_Ea[60] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H2 */
    wdot[7] -= q; /* H2O */
    wdot[11] += q; /* CH2O */
    wdot[22] -= q; /* SXCH2 */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[22];
    J[212] += dqdci;              /* dwdot[H2]/d[H2O] */
    J[217] -= dqdci;              /* dwdot[H2O]/d[H2O] */
    J[221] += dqdci;              /* dwdot[CH2O]/d[H2O] */
    J[232] -= dqdci;              /* dwdot[SXCH2]/d[H2O] */
    /* d()/d[SXCH2] */
    dqdci =  + k_f*sc[7];
    J[662] += dqdci;              /* dwdot[H2]/d[SXCH2] */
    J[667] -= dqdci;              /* dwdot[H2O]/d[SXCH2] */
    J[671] += dqdci;              /* dwdot[CH2O]/d[SXCH2] */
    J[682] -= dqdci;              /* dwdot[SXCH2]/d[SXCH2] */
    /* d()/dT */
    J[872] += dqdT;               /* dwdot[H2]/dT */
    J[877] -= dqdT;               /* dwdot[H2O]/dT */
    J[881] += dqdT;               /* dwdot[CH2O]/dT */
    J[892] -= dqdT;               /* dwdot[SXCH2]/dT */

    /*reaction 62: SXCH2 + OH <=> CH2O + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[22];
    k_f = prefactor_units[61] * fwd_A[61]
                * exp(fwd_beta[61] * tc[0] - activation_units[61] * fwd_Ea[61] * invT);
    dlnkfdT = fwd_beta[61] * invT + activation_units[61] * fwd_Ea[61] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[11];
    Kc = exp(-g_RT[1] + g_RT[4] - g_RT[11] + g_RT[22]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[22]) + (h_RT[1] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[4] -= q; /* OH */
    wdot[11] += q; /* CH2O */
    wdot[22] -= q; /* SXCH2 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[11];
    J[31] += dqdci;               /* dwdot[H]/d[H] */
    J[34] -= dqdci;               /* dwdot[OH]/d[H] */
    J[41] += dqdci;               /* dwdot[CH2O]/d[H] */
    J[52] -= dqdci;               /* dwdot[SXCH2]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[22];
    J[121] += dqdci;              /* dwdot[H]/d[OH] */
    J[124] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[131] += dqdci;              /* dwdot[CH2O]/d[OH] */
    J[142] -= dqdci;              /* dwdot[SXCH2]/d[OH] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[1];
    J[331] += dqdci;              /* dwdot[H]/d[CH2O] */
    J[334] -= dqdci;              /* dwdot[OH]/d[CH2O] */
    J[341] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[352] -= dqdci;              /* dwdot[SXCH2]/d[CH2O] */
    /* d()/d[SXCH2] */
    dqdci =  + k_f*sc[4];
    J[661] += dqdci;              /* dwdot[H]/d[SXCH2] */
    J[664] -= dqdci;              /* dwdot[OH]/d[SXCH2] */
    J[671] += dqdci;              /* dwdot[CH2O]/d[SXCH2] */
    J[682] -= dqdci;              /* dwdot[SXCH2]/d[SXCH2] */
    /* d()/dT */
    J[871] += dqdT;               /* dwdot[H]/dT */
    J[874] -= dqdT;               /* dwdot[OH]/dT */
    J[881] += dqdT;               /* dwdot[CH2O]/dT */
    J[892] -= dqdT;               /* dwdot[SXCH2]/dT */

    /*reaction 63: CH3 + OH => H2 + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[10];
    k_f = prefactor_units[62] * fwd_A[62]
                * exp(fwd_beta[62] * tc[0] - activation_units[62] * fwd_Ea[62] * invT);
    dlnkfdT = fwd_beta[62] * invT + activation_units[62] * fwd_Ea[62] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H2 */
    wdot[4] -= q; /* OH */
    wdot[10] -= q; /* CH3 */
    wdot[11] += q; /* CH2O */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[10];
    J[122] += dqdci;              /* dwdot[H2]/d[OH] */
    J[124] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[130] -= dqdci;              /* dwdot[CH3]/d[OH] */
    J[131] += dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[4];
    J[302] += dqdci;              /* dwdot[H2]/d[CH3] */
    J[304] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[310] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[311] += dqdci;              /* dwdot[CH2O]/d[CH3] */
    /* d()/dT */
    J[872] += dqdT;               /* dwdot[H2]/dT */
    J[874] -= dqdT;               /* dwdot[OH]/dT */
    J[880] -= dqdT;               /* dwdot[CH3]/dT */
    J[881] += dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 64: CH3 + H2O2 <=> CH4 + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[10];
    k_f = prefactor_units[63] * fwd_A[63]
                * exp(fwd_beta[63] * tc[0] - activation_units[63] * fwd_Ea[63] * invT);
    dlnkfdT = fwd_beta[63] * invT + activation_units[63] * fwd_Ea[63] * invT2;
    /* reverse */
    phi_r = sc[8]*sc[13];
    Kc = exp(g_RT[6] - g_RT[8] + g_RT[10] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[10]) + (h_RT[8] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[6] -= q; /* H2O2 */
    wdot[8] += q; /* HO2 */
    wdot[10] -= q; /* CH3 */
    wdot[13] += q; /* CH4 */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[10];
    J[186] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    J[188] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[190] -= dqdci;              /* dwdot[CH3]/d[H2O2] */
    J[193] += dqdci;              /* dwdot[CH4]/d[H2O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[13];
    J[246] -= dqdci;              /* dwdot[H2O2]/d[HO2] */
    J[248] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[250] -= dqdci;              /* dwdot[CH3]/d[HO2] */
    J[253] += dqdci;              /* dwdot[CH4]/d[HO2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[6];
    J[306] -= dqdci;              /* dwdot[H2O2]/d[CH3] */
    J[308] += dqdci;              /* dwdot[HO2]/d[CH3] */
    J[310] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[313] += dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[8];
    J[396] -= dqdci;              /* dwdot[H2O2]/d[CH4] */
    J[398] += dqdci;              /* dwdot[HO2]/d[CH4] */
    J[400] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[403] += dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[876] -= dqdT;               /* dwdot[H2O2]/dT */
    J[878] += dqdT;               /* dwdot[HO2]/dT */
    J[880] -= dqdT;               /* dwdot[CH3]/dT */
    J[883] += dqdT;               /* dwdot[CH4]/dT */

    /*reaction 65: CH3 + O2 <=> CH2O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[10];
    k_f = prefactor_units[64] * fwd_A[64]
                * exp(fwd_beta[64] * tc[0] - activation_units[64] * fwd_Ea[64] * invT);
    dlnkfdT = fwd_beta[64] * invT + activation_units[64] * fwd_Ea[64] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[11];
    Kc = exp(-g_RT[4] + g_RT[5] + g_RT[10] - g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[10]) + (h_RT[4] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] += q; /* OH */
    wdot[5] -= q; /* O2 */
    wdot[10] -= q; /* CH3 */
    wdot[11] += q; /* CH2O */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[11];
    J[124] += dqdci;              /* dwdot[OH]/d[OH] */
    J[125] -= dqdci;              /* dwdot[O2]/d[OH] */
    J[130] -= dqdci;              /* dwdot[CH3]/d[OH] */
    J[131] += dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[10];
    J[154] += dqdci;              /* dwdot[OH]/d[O2] */
    J[155] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[160] -= dqdci;              /* dwdot[CH3]/d[O2] */
    J[161] += dqdci;              /* dwdot[CH2O]/d[O2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[5];
    J[304] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[305] -= dqdci;              /* dwdot[O2]/d[CH3] */
    J[310] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[311] += dqdci;              /* dwdot[CH2O]/d[CH3] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[4];
    J[334] += dqdci;              /* dwdot[OH]/d[CH2O] */
    J[335] -= dqdci;              /* dwdot[O2]/d[CH2O] */
    J[340] -= dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[341] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[874] += dqdT;               /* dwdot[OH]/dT */
    J[875] -= dqdT;               /* dwdot[O2]/dT */
    J[880] -= dqdT;               /* dwdot[CH3]/dT */
    J[881] += dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 66: CH3 + CH <=> C2H3 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[10]*sc[19];
    k_f = prefactor_units[65] * fwd_A[65]
                * exp(fwd_beta[65] * tc[0] - activation_units[65] * fwd_Ea[65] * invT);
    dlnkfdT = fwd_beta[65] * invT + activation_units[65] * fwd_Ea[65] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[23];
    Kc = exp(-g_RT[1] + g_RT[10] + g_RT[19] - g_RT[23]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[10] + h_RT[19]) + (h_RT[1] + h_RT[23]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[10] -= q; /* CH3 */
    wdot[19] -= q; /* CH */
    wdot[23] += q; /* C2H3 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[23];
    J[31] += dqdci;               /* dwdot[H]/d[H] */
    J[40] -= dqdci;               /* dwdot[CH3]/d[H] */
    J[49] -= dqdci;               /* dwdot[CH]/d[H] */
    J[53] += dqdci;               /* dwdot[C2H3]/d[H] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[19];
    J[301] += dqdci;              /* dwdot[H]/d[CH3] */
    J[310] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[319] -= dqdci;              /* dwdot[CH]/d[CH3] */
    J[323] += dqdci;              /* dwdot[C2H3]/d[CH3] */
    /* d()/d[CH] */
    dqdci =  + k_f*sc[10];
    J[571] += dqdci;              /* dwdot[H]/d[CH] */
    J[580] -= dqdci;              /* dwdot[CH3]/d[CH] */
    J[589] -= dqdci;              /* dwdot[CH]/d[CH] */
    J[593] += dqdci;              /* dwdot[C2H3]/d[CH] */
    /* d()/d[C2H3] */
    dqdci =  - k_r*sc[1];
    J[691] += dqdci;              /* dwdot[H]/d[C2H3] */
    J[700] -= dqdci;              /* dwdot[CH3]/d[C2H3] */
    J[709] -= dqdci;              /* dwdot[CH]/d[C2H3] */
    J[713] += dqdci;              /* dwdot[C2H3]/d[C2H3] */
    /* d()/dT */
    J[871] += dqdT;               /* dwdot[H]/dT */
    J[880] -= dqdT;               /* dwdot[CH3]/dT */
    J[889] -= dqdT;               /* dwdot[CH]/dT */
    J[893] += dqdT;               /* dwdot[C2H3]/dT */

    /*reaction 67: CH3 + O <=> CH2O + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[10];
    k_f = prefactor_units[66] * fwd_A[66]
                * exp(fwd_beta[66] * tc[0] - activation_units[66] * fwd_Ea[66] * invT);
    dlnkfdT = fwd_beta[66] * invT + activation_units[66] * fwd_Ea[66] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[11];
    Kc = exp(-g_RT[1] + g_RT[3] + g_RT[10] - g_RT[11]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[10]) + (h_RT[1] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[3] -= q; /* O */
    wdot[10] -= q; /* CH3 */
    wdot[11] += q; /* CH2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[11];
    J[31] += dqdci;               /* dwdot[H]/d[H] */
    J[33] -= dqdci;               /* dwdot[O]/d[H] */
    J[40] -= dqdci;               /* dwdot[CH3]/d[H] */
    J[41] += dqdci;               /* dwdot[CH2O]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[10];
    J[91] += dqdci;               /* dwdot[H]/d[O] */
    J[93] -= dqdci;               /* dwdot[O]/d[O] */
    J[100] -= dqdci;              /* dwdot[CH3]/d[O] */
    J[101] += dqdci;              /* dwdot[CH2O]/d[O] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[3];
    J[301] += dqdci;              /* dwdot[H]/d[CH3] */
    J[303] -= dqdci;              /* dwdot[O]/d[CH3] */
    J[310] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[311] += dqdci;              /* dwdot[CH2O]/d[CH3] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[1];
    J[331] += dqdci;              /* dwdot[H]/d[CH2O] */
    J[333] -= dqdci;              /* dwdot[O]/d[CH2O] */
    J[340] -= dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[341] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    /* d()/dT */
    J[871] += dqdT;               /* dwdot[H]/dT */
    J[873] -= dqdT;               /* dwdot[O]/dT */
    J[880] -= dqdT;               /* dwdot[CH3]/dT */
    J[881] += dqdT;               /* dwdot[CH2O]/dT */

    /*reaction 68: CH3 + C <=> C2H2 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[10]*sc[18];
    k_f = prefactor_units[67] * fwd_A[67]
                * exp(fwd_beta[67] * tc[0] - activation_units[67] * fwd_Ea[67] * invT);
    dlnkfdT = fwd_beta[67] * invT + activation_units[67] * fwd_Ea[67] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[14];
    Kc = exp(-g_RT[1] + g_RT[10] - g_RT[14] + g_RT[18]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[10] + h_RT[18]) + (h_RT[1] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[10] -= q; /* CH3 */
    wdot[14] += q; /* C2H2 */
    wdot[18] -= q; /* C */
    /* d()/d[H] */
    dqdci =  - k_r*sc[14];
    J[31] += dqdci;               /* dwdot[H]/d[H] */
    J[40] -= dqdci;               /* dwdot[CH3]/d[H] */
    J[44] += dqdci;               /* dwdot[C2H2]/d[H] */
    J[48] -= dqdci;               /* dwdot[C]/d[H] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[18];
    J[301] += dqdci;              /* dwdot[H]/d[CH3] */
    J[310] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[314] += dqdci;              /* dwdot[C2H2]/d[CH3] */
    J[318] -= dqdci;              /* dwdot[C]/d[CH3] */
    /* d()/d[C2H2] */
    dqdci =  - k_r*sc[1];
    J[421] += dqdci;              /* dwdot[H]/d[C2H2] */
    J[430] -= dqdci;              /* dwdot[CH3]/d[C2H2] */
    J[434] += dqdci;              /* dwdot[C2H2]/d[C2H2] */
    J[438] -= dqdci;              /* dwdot[C]/d[C2H2] */
    /* d()/d[C] */
    dqdci =  + k_f*sc[10];
    J[541] += dqdci;              /* dwdot[H]/d[C] */
    J[550] -= dqdci;              /* dwdot[CH3]/d[C] */
    J[554] += dqdci;              /* dwdot[C2H2]/d[C] */
    J[558] -= dqdci;              /* dwdot[C]/d[C] */
    /* d()/dT */
    J[871] += dqdT;               /* dwdot[H]/dT */
    J[880] -= dqdT;               /* dwdot[CH3]/dT */
    J[884] += dqdT;               /* dwdot[C2H2]/dT */
    J[888] -= dqdT;               /* dwdot[C]/dT */

    /*reaction 69: CH3 + OH <=> TXCH2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[10];
    k_f = prefactor_units[68] * fwd_A[68]
                * exp(fwd_beta[68] * tc[0] - activation_units[68] * fwd_Ea[68] * invT);
    dlnkfdT = fwd_beta[68] * invT + activation_units[68] * fwd_Ea[68] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[21];
    Kc = exp(g_RT[4] - g_RT[7] + g_RT[10] - g_RT[21]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[10]) + (h_RT[7] + h_RT[21]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[7] += q; /* H2O */
    wdot[10] -= q; /* CH3 */
    wdot[21] += q; /* TXCH2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[10];
    J[124] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[127] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[130] -= dqdci;              /* dwdot[CH3]/d[OH] */
    J[141] += dqdci;              /* dwdot[TXCH2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[21];
    J[214] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[217] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[220] -= dqdci;              /* dwdot[CH3]/d[H2O] */
    J[231] += dqdci;              /* dwdot[TXCH2]/d[H2O] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[4];
    J[304] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[307] += dqdci;              /* dwdot[H2O]/d[CH3] */
    J[310] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[321] += dqdci;              /* dwdot[TXCH2]/d[CH3] */
    /* d()/d[TXCH2] */
    dqdci =  - k_r*sc[7];
    J[634] -= dqdci;              /* dwdot[OH]/d[TXCH2] */
    J[637] += dqdci;              /* dwdot[H2O]/d[TXCH2] */
    J[640] -= dqdci;              /* dwdot[CH3]/d[TXCH2] */
    J[651] += dqdci;              /* dwdot[TXCH2]/d[TXCH2] */
    /* d()/dT */
    J[874] -= dqdT;               /* dwdot[OH]/dT */
    J[877] += dqdT;               /* dwdot[H2O]/dT */
    J[880] -= dqdT;               /* dwdot[CH3]/dT */
    J[891] += dqdT;               /* dwdot[TXCH2]/dT */

    /*reaction 70: CH3 + SXCH2 <=> C2H4 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[10]*sc[22];
    k_f = prefactor_units[69] * fwd_A[69]
                * exp(fwd_beta[69] * tc[0] - activation_units[69] * fwd_Ea[69] * invT);
    dlnkfdT = fwd_beta[69] * invT + activation_units[69] * fwd_Ea[69] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[15];
    Kc = exp(-g_RT[1] + g_RT[10] - g_RT[15] + g_RT[22]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[10] + h_RT[22]) + (h_RT[1] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[10] -= q; /* CH3 */
    wdot[15] += q; /* C2H4 */
    wdot[22] -= q; /* SXCH2 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[15];
    J[31] += dqdci;               /* dwdot[H]/d[H] */
    J[40] -= dqdci;               /* dwdot[CH3]/d[H] */
    J[45] += dqdci;               /* dwdot[C2H4]/d[H] */
    J[52] -= dqdci;               /* dwdot[SXCH2]/d[H] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[22];
    J[301] += dqdci;              /* dwdot[H]/d[CH3] */
    J[310] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[315] += dqdci;              /* dwdot[C2H4]/d[CH3] */
    J[322] -= dqdci;              /* dwdot[SXCH2]/d[CH3] */
    /* d()/d[C2H4] */
    dqdci =  - k_r*sc[1];
    J[451] += dqdci;              /* dwdot[H]/d[C2H4] */
    J[460] -= dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[465] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[472] -= dqdci;              /* dwdot[SXCH2]/d[C2H4] */
    /* d()/d[SXCH2] */
    dqdci =  + k_f*sc[10];
    J[661] += dqdci;              /* dwdot[H]/d[SXCH2] */
    J[670] -= dqdci;              /* dwdot[CH3]/d[SXCH2] */
    J[675] += dqdci;              /* dwdot[C2H4]/d[SXCH2] */
    J[682] -= dqdci;              /* dwdot[SXCH2]/d[SXCH2] */
    /* d()/dT */
    J[871] += dqdT;               /* dwdot[H]/dT */
    J[880] -= dqdT;               /* dwdot[CH3]/dT */
    J[885] += dqdT;               /* dwdot[C2H4]/dT */
    J[892] -= dqdT;               /* dwdot[SXCH2]/dT */

    /*reaction 71: CH3 + OH <=> SXCH2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[10];
    k_f = prefactor_units[70] * fwd_A[70]
                * exp(fwd_beta[70] * tc[0] - activation_units[70] * fwd_Ea[70] * invT);
    dlnkfdT = fwd_beta[70] * invT + activation_units[70] * fwd_Ea[70] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[22];
    Kc = exp(g_RT[4] - g_RT[7] + g_RT[10] - g_RT[22]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[10]) + (h_RT[7] + h_RT[22]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[7] += q; /* H2O */
    wdot[10] -= q; /* CH3 */
    wdot[22] += q; /* SXCH2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[10];
    J[124] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[127] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[130] -= dqdci;              /* dwdot[CH3]/d[OH] */
    J[142] += dqdci;              /* dwdot[SXCH2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[22];
    J[214] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[217] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[220] -= dqdci;              /* dwdot[CH3]/d[H2O] */
    J[232] += dqdci;              /* dwdot[SXCH2]/d[H2O] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[4];
    J[304] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[307] += dqdci;              /* dwdot[H2O]/d[CH3] */
    J[310] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[322] += dqdci;              /* dwdot[SXCH2]/d[CH3] */
    /* d()/d[SXCH2] */
    dqdci =  - k_r*sc[7];
    J[664] -= dqdci;              /* dwdot[OH]/d[SXCH2] */
    J[667] += dqdci;              /* dwdot[H2O]/d[SXCH2] */
    J[670] -= dqdci;              /* dwdot[CH3]/d[SXCH2] */
    J[682] += dqdci;              /* dwdot[SXCH2]/d[SXCH2] */
    /* d()/dT */
    J[874] -= dqdT;               /* dwdot[OH]/dT */
    J[877] += dqdT;               /* dwdot[H2O]/dT */
    J[880] -= dqdT;               /* dwdot[CH3]/dT */
    J[892] += dqdT;               /* dwdot[SXCH2]/dT */

    /*reaction 72: 2.000000 CH3 <=> C2H5 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[10], 2.000000);
    k_f = prefactor_units[71] * fwd_A[71]
                * exp(fwd_beta[71] * tc[0] - activation_units[71] * fwd_Ea[71] * invT);
    dlnkfdT = fwd_beta[71] * invT + activation_units[71] * fwd_Ea[71] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[24];
    Kc = exp(-g_RT[1] + 2.000000*g_RT[10] - g_RT[24]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[10]) + (h_RT[1] + h_RT[24]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[10] -= 2 * q; /* CH3 */
    wdot[24] += q; /* C2H5 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[24];
    J[31] += dqdci;               /* dwdot[H]/d[H] */
    J[40] += -2 * dqdci;          /* dwdot[CH3]/d[H] */
    J[54] += dqdci;               /* dwdot[C2H5]/d[H] */
    /* d()/d[CH3] */
    dqdci =  + k_f*2.000000*sc[10];
    J[301] += dqdci;              /* dwdot[H]/d[CH3] */
    J[310] += -2 * dqdci;         /* dwdot[CH3]/d[CH3] */
    J[324] += dqdci;              /* dwdot[C2H5]/d[CH3] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[1];
    J[721] += dqdci;              /* dwdot[H]/d[C2H5] */
    J[730] += -2 * dqdci;         /* dwdot[CH3]/d[C2H5] */
    J[744] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    /* d()/dT */
    J[871] += dqdT;               /* dwdot[H]/dT */
    J[880] += -2 * dqdT;          /* dwdot[CH3]/dT */
    J[894] += dqdT;               /* dwdot[C2H5]/dT */

    /*reaction 73: CH3 + HO2 <=> CH4 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[10];
    k_f = prefactor_units[72] * fwd_A[72]
                * exp(fwd_beta[72] * tc[0] - activation_units[72] * fwd_Ea[72] * invT);
    dlnkfdT = fwd_beta[72] * invT + activation_units[72] * fwd_Ea[72] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[13];
    Kc = exp(-g_RT[5] + g_RT[8] + g_RT[10] - g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[10]) + (h_RT[5] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[8] -= q; /* HO2 */
    wdot[10] -= q; /* CH3 */
    wdot[13] += q; /* CH4 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[13];
    J[155] += dqdci;              /* dwdot[O2]/d[O2] */
    J[158] -= dqdci;              /* dwdot[HO2]/d[O2] */
    J[160] -= dqdci;              /* dwdot[CH3]/d[O2] */
    J[163] += dqdci;              /* dwdot[CH4]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[10];
    J[245] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[248] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[250] -= dqdci;              /* dwdot[CH3]/d[HO2] */
    J[253] += dqdci;              /* dwdot[CH4]/d[HO2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[8];
    J[305] += dqdci;              /* dwdot[O2]/d[CH3] */
    J[308] -= dqdci;              /* dwdot[HO2]/d[CH3] */
    J[310] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[313] += dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[5];
    J[395] += dqdci;              /* dwdot[O2]/d[CH4] */
    J[398] -= dqdci;              /* dwdot[HO2]/d[CH4] */
    J[400] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[403] += dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[875] += dqdT;               /* dwdot[O2]/dT */
    J[878] -= dqdT;               /* dwdot[HO2]/dT */
    J[880] -= dqdT;               /* dwdot[CH3]/dT */
    J[883] += dqdT;               /* dwdot[CH4]/dT */

    /*reaction 74: CH3 + TXCH2 <=> C2H4 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[10]*sc[21];
    k_f = prefactor_units[73] * fwd_A[73]
                * exp(fwd_beta[73] * tc[0] - activation_units[73] * fwd_Ea[73] * invT);
    dlnkfdT = fwd_beta[73] * invT + activation_units[73] * fwd_Ea[73] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[15];
    Kc = exp(-g_RT[1] + g_RT[10] - g_RT[15] + g_RT[21]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[10] + h_RT[21]) + (h_RT[1] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[10] -= q; /* CH3 */
    wdot[15] += q; /* C2H4 */
    wdot[21] -= q; /* TXCH2 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[15];
    J[31] += dqdci;               /* dwdot[H]/d[H] */
    J[40] -= dqdci;               /* dwdot[CH3]/d[H] */
    J[45] += dqdci;               /* dwdot[C2H4]/d[H] */
    J[51] -= dqdci;               /* dwdot[TXCH2]/d[H] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[21];
    J[301] += dqdci;              /* dwdot[H]/d[CH3] */
    J[310] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[315] += dqdci;              /* dwdot[C2H4]/d[CH3] */
    J[321] -= dqdci;              /* dwdot[TXCH2]/d[CH3] */
    /* d()/d[C2H4] */
    dqdci =  - k_r*sc[1];
    J[451] += dqdci;              /* dwdot[H]/d[C2H4] */
    J[460] -= dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[465] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[471] -= dqdci;              /* dwdot[TXCH2]/d[C2H4] */
    /* d()/d[TXCH2] */
    dqdci =  + k_f*sc[10];
    J[631] += dqdci;              /* dwdot[H]/d[TXCH2] */
    J[640] -= dqdci;              /* dwdot[CH3]/d[TXCH2] */
    J[645] += dqdci;              /* dwdot[C2H4]/d[TXCH2] */
    J[651] -= dqdci;              /* dwdot[TXCH2]/d[TXCH2] */
    /* d()/dT */
    J[871] += dqdT;               /* dwdot[H]/dT */
    J[880] -= dqdT;               /* dwdot[CH3]/dT */
    J[885] += dqdT;               /* dwdot[C2H4]/dT */
    J[891] -= dqdT;               /* dwdot[TXCH2]/dT */

    /*reaction 75: CH3 + O => H + H2 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[10];
    k_f = prefactor_units[74] * fwd_A[74]
                * exp(fwd_beta[74] * tc[0] - activation_units[74] * fwd_Ea[74] * invT);
    dlnkfdT = fwd_beta[74] * invT + activation_units[74] * fwd_Ea[74] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[2] += q; /* H2 */
    wdot[3] -= q; /* O */
    wdot[9] += q; /* CO */
    wdot[10] -= q; /* CH3 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[10];
    J[91] += dqdci;               /* dwdot[H]/d[O] */
    J[92] += dqdci;               /* dwdot[H2]/d[O] */
    J[93] -= dqdci;               /* dwdot[O]/d[O] */
    J[99] += dqdci;               /* dwdot[CO]/d[O] */
    J[100] -= dqdci;              /* dwdot[CH3]/d[O] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[3];
    J[301] += dqdci;              /* dwdot[H]/d[CH3] */
    J[302] += dqdci;              /* dwdot[H2]/d[CH3] */
    J[303] -= dqdci;              /* dwdot[O]/d[CH3] */
    J[309] += dqdci;              /* dwdot[CO]/d[CH3] */
    J[310] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[871] += dqdT;               /* dwdot[H]/dT */
    J[872] += dqdT;               /* dwdot[H2]/dT */
    J[873] -= dqdT;               /* dwdot[O]/dT */
    J[879] += dqdT;               /* dwdot[CO]/dT */
    J[880] -= dqdT;               /* dwdot[CH3]/dT */

    /*reaction 76: CH4 + CH <=> C2H4 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[13]*sc[19];
    k_f = prefactor_units[75] * fwd_A[75]
                * exp(fwd_beta[75] * tc[0] - activation_units[75] * fwd_Ea[75] * invT);
    dlnkfdT = fwd_beta[75] * invT + activation_units[75] * fwd_Ea[75] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[15];
    Kc = exp(-g_RT[1] + g_RT[13] - g_RT[15] + g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[13] + h_RT[19]) + (h_RT[1] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[13] -= q; /* CH4 */
    wdot[15] += q; /* C2H4 */
    wdot[19] -= q; /* CH */
    /* d()/d[H] */
    dqdci =  - k_r*sc[15];
    J[31] += dqdci;               /* dwdot[H]/d[H] */
    J[43] -= dqdci;               /* dwdot[CH4]/d[H] */
    J[45] += dqdci;               /* dwdot[C2H4]/d[H] */
    J[49] -= dqdci;               /* dwdot[CH]/d[H] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[19];
    J[391] += dqdci;              /* dwdot[H]/d[CH4] */
    J[403] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    J[405] += dqdci;              /* dwdot[C2H4]/d[CH4] */
    J[409] -= dqdci;              /* dwdot[CH]/d[CH4] */
    /* d()/d[C2H4] */
    dqdci =  - k_r*sc[1];
    J[451] += dqdci;              /* dwdot[H]/d[C2H4] */
    J[463] -= dqdci;              /* dwdot[CH4]/d[C2H4] */
    J[465] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[469] -= dqdci;              /* dwdot[CH]/d[C2H4] */
    /* d()/d[CH] */
    dqdci =  + k_f*sc[13];
    J[571] += dqdci;              /* dwdot[H]/d[CH] */
    J[583] -= dqdci;              /* dwdot[CH4]/d[CH] */
    J[585] += dqdci;              /* dwdot[C2H4]/d[CH] */
    J[589] -= dqdci;              /* dwdot[CH]/d[CH] */
    /* d()/dT */
    J[871] += dqdT;               /* dwdot[H]/dT */
    J[883] -= dqdT;               /* dwdot[CH4]/dT */
    J[885] += dqdT;               /* dwdot[C2H4]/dT */
    J[889] -= dqdT;               /* dwdot[CH]/dT */

    /*reaction 77: CH4 + SXCH2 <=> 2.000000 CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[13]*sc[22];
    k_f = prefactor_units[76] * fwd_A[76]
                * exp(fwd_beta[76] * tc[0] - activation_units[76] * fwd_Ea[76] * invT);
    dlnkfdT = fwd_beta[76] * invT + activation_units[76] * fwd_Ea[76] * invT2;
    /* reverse */
    phi_r = pow(sc[10], 2.000000);
    Kc = exp(-2.000000*g_RT[10] + g_RT[13] + g_RT[22]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[13] + h_RT[22]) + (2.000000*h_RT[10]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[10] += 2 * q; /* CH3 */
    wdot[13] -= q; /* CH4 */
    wdot[22] -= q; /* SXCH2 */
    /* d()/d[CH3] */
    dqdci =  - k_r*2.000000*sc[10];
    J[310] += 2 * dqdci;          /* dwdot[CH3]/d[CH3] */
    J[313] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    J[322] -= dqdci;              /* dwdot[SXCH2]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[22];
    J[400] += 2 * dqdci;          /* dwdot[CH3]/d[CH4] */
    J[403] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    J[412] -= dqdci;              /* dwdot[SXCH2]/d[CH4] */
    /* d()/d[SXCH2] */
    dqdci =  + k_f*sc[13];
    J[670] += 2 * dqdci;          /* dwdot[CH3]/d[SXCH2] */
    J[673] -= dqdci;              /* dwdot[CH4]/d[SXCH2] */
    J[682] -= dqdci;              /* dwdot[SXCH2]/d[SXCH2] */
    /* d()/dT */
    J[880] += 2 * dqdT;           /* dwdot[CH3]/dT */
    J[883] -= dqdT;               /* dwdot[CH4]/dT */
    J[892] -= dqdT;               /* dwdot[SXCH2]/dT */

    /*reaction 78: CH4 + O <=> CH3 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[13];
    k_f = prefactor_units[77] * fwd_A[77]
                * exp(fwd_beta[77] * tc[0] - activation_units[77] * fwd_Ea[77] * invT);
    dlnkfdT = fwd_beta[77] * invT + activation_units[77] * fwd_Ea[77] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[10];
    Kc = exp(g_RT[3] - g_RT[4] - g_RT[10] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[13]) + (h_RT[4] + h_RT[10]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O */
    wdot[4] += q; /* OH */
    wdot[10] += q; /* CH3 */
    wdot[13] -= q; /* CH4 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[13];
    J[93] -= dqdci;               /* dwdot[O]/d[O] */
    J[94] += dqdci;               /* dwdot[OH]/d[O] */
    J[100] += dqdci;              /* dwdot[CH3]/d[O] */
    J[103] -= dqdci;              /* dwdot[CH4]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[10];
    J[123] -= dqdci;              /* dwdot[O]/d[OH] */
    J[124] += dqdci;              /* dwdot[OH]/d[OH] */
    J[130] += dqdci;              /* dwdot[CH3]/d[OH] */
    J[133] -= dqdci;              /* dwdot[CH4]/d[OH] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[4];
    J[303] -= dqdci;              /* dwdot[O]/d[CH3] */
    J[304] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[310] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[313] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[3];
    J[393] -= dqdci;              /* dwdot[O]/d[CH4] */
    J[394] += dqdci;              /* dwdot[OH]/d[CH4] */
    J[400] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[403] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[873] -= dqdT;               /* dwdot[O]/dT */
    J[874] += dqdT;               /* dwdot[OH]/dT */
    J[880] += dqdT;               /* dwdot[CH3]/dT */
    J[883] -= dqdT;               /* dwdot[CH4]/dT */

    /*reaction 79: CH4 + OH <=> CH3 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[13];
    k_f = prefactor_units[78] * fwd_A[78]
                * exp(fwd_beta[78] * tc[0] - activation_units[78] * fwd_Ea[78] * invT);
    dlnkfdT = fwd_beta[78] * invT + activation_units[78] * fwd_Ea[78] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[10];
    Kc = exp(g_RT[4] - g_RT[7] - g_RT[10] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[13]) + (h_RT[7] + h_RT[10]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[7] += q; /* H2O */
    wdot[10] += q; /* CH3 */
    wdot[13] -= q; /* CH4 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[13];
    J[124] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[127] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[130] += dqdci;              /* dwdot[CH3]/d[OH] */
    J[133] -= dqdci;              /* dwdot[CH4]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[10];
    J[214] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[217] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[220] += dqdci;              /* dwdot[CH3]/d[H2O] */
    J[223] -= dqdci;              /* dwdot[CH4]/d[H2O] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[7];
    J[304] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[307] += dqdci;              /* dwdot[H2O]/d[CH3] */
    J[310] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[313] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[4];
    J[394] -= dqdci;              /* dwdot[OH]/d[CH4] */
    J[397] += dqdci;              /* dwdot[H2O]/d[CH4] */
    J[400] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[403] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[874] -= dqdT;               /* dwdot[OH]/dT */
    J[877] += dqdT;               /* dwdot[H2O]/dT */
    J[880] += dqdT;               /* dwdot[CH3]/dT */
    J[883] -= dqdT;               /* dwdot[CH4]/dT */

    /*reaction 80: CH4 + TXCH2 <=> 2.000000 CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[13]*sc[21];
    k_f = prefactor_units[79] * fwd_A[79]
                * exp(fwd_beta[79] * tc[0] - activation_units[79] * fwd_Ea[79] * invT);
    dlnkfdT = fwd_beta[79] * invT + activation_units[79] * fwd_Ea[79] * invT2;
    /* reverse */
    phi_r = pow(sc[10], 2.000000);
    Kc = exp(-2.000000*g_RT[10] + g_RT[13] + g_RT[21]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[13] + h_RT[21]) + (2.000000*h_RT[10]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[10] += 2 * q; /* CH3 */
    wdot[13] -= q; /* CH4 */
    wdot[21] -= q; /* TXCH2 */
    /* d()/d[CH3] */
    dqdci =  - k_r*2.000000*sc[10];
    J[310] += 2 * dqdci;          /* dwdot[CH3]/d[CH3] */
    J[313] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    J[321] -= dqdci;              /* dwdot[TXCH2]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[21];
    J[400] += 2 * dqdci;          /* dwdot[CH3]/d[CH4] */
    J[403] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    J[411] -= dqdci;              /* dwdot[TXCH2]/d[CH4] */
    /* d()/d[TXCH2] */
    dqdci =  + k_f*sc[13];
    J[640] += 2 * dqdci;          /* dwdot[CH3]/d[TXCH2] */
    J[643] -= dqdci;              /* dwdot[CH4]/d[TXCH2] */
    J[651] -= dqdci;              /* dwdot[TXCH2]/d[TXCH2] */
    /* d()/dT */
    J[880] += 2 * dqdT;           /* dwdot[CH3]/dT */
    J[883] -= dqdT;               /* dwdot[CH4]/dT */
    J[891] -= dqdT;               /* dwdot[TXCH2]/dT */

    /*reaction 81: CH4 + H <=> CH3 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[13];
    k_f = prefactor_units[80] * fwd_A[80]
                * exp(fwd_beta[80] * tc[0] - activation_units[80] * fwd_Ea[80] * invT);
    dlnkfdT = fwd_beta[80] * invT + activation_units[80] * fwd_Ea[80] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[10];
    Kc = exp(g_RT[1] - g_RT[2] - g_RT[10] + g_RT[13]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[13]) + (h_RT[2] + h_RT[10]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[2] += q; /* H2 */
    wdot[10] += q; /* CH3 */
    wdot[13] -= q; /* CH4 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[13];
    J[31] -= dqdci;               /* dwdot[H]/d[H] */
    J[32] += dqdci;               /* dwdot[H2]/d[H] */
    J[40] += dqdci;               /* dwdot[CH3]/d[H] */
    J[43] -= dqdci;               /* dwdot[CH4]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[10];
    J[61] -= dqdci;               /* dwdot[H]/d[H2] */
    J[62] += dqdci;               /* dwdot[H2]/d[H2] */
    J[70] += dqdci;               /* dwdot[CH3]/d[H2] */
    J[73] -= dqdci;               /* dwdot[CH4]/d[H2] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[2];
    J[301] -= dqdci;              /* dwdot[H]/d[CH3] */
    J[302] += dqdci;              /* dwdot[H2]/d[CH3] */
    J[310] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[313] -= dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[1];
    J[391] -= dqdci;              /* dwdot[H]/d[CH4] */
    J[392] += dqdci;              /* dwdot[H2]/d[CH4] */
    J[400] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[403] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[871] -= dqdT;               /* dwdot[H]/dT */
    J[872] += dqdT;               /* dwdot[H2]/dT */
    J[880] += dqdT;               /* dwdot[CH3]/dT */
    J[883] -= dqdT;               /* dwdot[CH4]/dT */

    /*reaction 82: SXCH2 + CO <=> TXCH2 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[22];
    k_f = prefactor_units[81] * fwd_A[81]
                * exp(fwd_beta[81] * tc[0] - activation_units[81] * fwd_Ea[81] * invT);
    dlnkfdT = fwd_beta[81] * invT + activation_units[81] * fwd_Ea[81] * invT2;
    /* reverse */
    phi_r = sc[9]*sc[21];
    Kc = exp(g_RT[9] - g_RT[9] - g_RT[21] + g_RT[22]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[9] + h_RT[22]) + (h_RT[9] + h_RT[21]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[21] += q; /* TXCH2 */
    wdot[22] -= q; /* SXCH2 */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[22] - k_r*sc[21];
    J[291] += dqdci;              /* dwdot[TXCH2]/d[CO] */
    J[292] -= dqdci;              /* dwdot[SXCH2]/d[CO] */
    /* d()/d[TXCH2] */
    dqdci =  - k_r*sc[9];
    J[651] += dqdci;              /* dwdot[TXCH2]/d[TXCH2] */
    J[652] -= dqdci;              /* dwdot[SXCH2]/d[TXCH2] */
    /* d()/d[SXCH2] */
    dqdci =  + k_f*sc[9];
    J[681] += dqdci;              /* dwdot[TXCH2]/d[SXCH2] */
    J[682] -= dqdci;              /* dwdot[SXCH2]/d[SXCH2] */
    /* d()/dT */
    J[891] += dqdT;               /* dwdot[TXCH2]/dT */
    J[892] -= dqdT;               /* dwdot[SXCH2]/dT */

    /*reaction 83: CO + O2 <=> CO2 + O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[9];
    k_f = prefactor_units[82] * fwd_A[82]
                * exp(fwd_beta[82] * tc[0] - activation_units[82] * fwd_Ea[82] * invT);
    dlnkfdT = fwd_beta[82] * invT + activation_units[82] * fwd_Ea[82] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[12];
    Kc = exp(-g_RT[3] + g_RT[5] + g_RT[9] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[9]) + (h_RT[3] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += q; /* O */
    wdot[5] -= q; /* O2 */
    wdot[9] -= q; /* CO */
    wdot[12] += q; /* CO2 */
    /* d()/d[O] */
    dqdci =  - k_r*sc[12];
    J[93] += dqdci;               /* dwdot[O]/d[O] */
    J[95] -= dqdci;               /* dwdot[O2]/d[O] */
    J[99] -= dqdci;               /* dwdot[CO]/d[O] */
    J[102] += dqdci;              /* dwdot[CO2]/d[O] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[9];
    J[153] += dqdci;              /* dwdot[O]/d[O2] */
    J[155] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[159] -= dqdci;              /* dwdot[CO]/d[O2] */
    J[162] += dqdci;              /* dwdot[CO2]/d[O2] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[5];
    J[273] += dqdci;              /* dwdot[O]/d[CO] */
    J[275] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[279] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[282] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[3];
    J[363] += dqdci;              /* dwdot[O]/d[CO2] */
    J[365] -= dqdci;              /* dwdot[O2]/d[CO2] */
    J[369] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[372] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[873] += dqdT;               /* dwdot[O]/dT */
    J[875] -= dqdT;               /* dwdot[O2]/dT */
    J[879] -= dqdT;               /* dwdot[CO]/dT */
    J[882] += dqdT;               /* dwdot[CO2]/dT */

    /*reaction 84: CO + OH <=> CO2 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[9];
    k_f = prefactor_units[83] * fwd_A[83]
                * exp(fwd_beta[83] * tc[0] - activation_units[83] * fwd_Ea[83] * invT);
    dlnkfdT = fwd_beta[83] * invT + activation_units[83] * fwd_Ea[83] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[12];
    Kc = exp(-g_RT[1] + g_RT[4] + g_RT[9] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[9]) + (h_RT[1] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[4] -= q; /* OH */
    wdot[9] -= q; /* CO */
    wdot[12] += q; /* CO2 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[12];
    J[31] += dqdci;               /* dwdot[H]/d[H] */
    J[34] -= dqdci;               /* dwdot[OH]/d[H] */
    J[39] -= dqdci;               /* dwdot[CO]/d[H] */
    J[42] += dqdci;               /* dwdot[CO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[9];
    J[121] += dqdci;              /* dwdot[H]/d[OH] */
    J[124] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[129] -= dqdci;              /* dwdot[CO]/d[OH] */
    J[132] += dqdci;              /* dwdot[CO2]/d[OH] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[4];
    J[271] += dqdci;              /* dwdot[H]/d[CO] */
    J[274] -= dqdci;              /* dwdot[OH]/d[CO] */
    J[279] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[282] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[1];
    J[361] += dqdci;              /* dwdot[H]/d[CO2] */
    J[364] -= dqdci;              /* dwdot[OH]/d[CO2] */
    J[369] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[372] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[871] += dqdT;               /* dwdot[H]/dT */
    J[874] -= dqdT;               /* dwdot[OH]/dT */
    J[879] -= dqdT;               /* dwdot[CO]/dT */
    J[882] += dqdT;               /* dwdot[CO2]/dT */

    /*reaction 85: CO + OH <=> CO2 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[9];
    k_f = prefactor_units[84] * fwd_A[84]
                * exp(fwd_beta[84] * tc[0] - activation_units[84] * fwd_Ea[84] * invT);
    dlnkfdT = fwd_beta[84] * invT + activation_units[84] * fwd_Ea[84] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[12];
    Kc = exp(-g_RT[1] + g_RT[4] + g_RT[9] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[9]) + (h_RT[1] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[4] -= q; /* OH */
    wdot[9] -= q; /* CO */
    wdot[12] += q; /* CO2 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[12];
    J[31] += dqdci;               /* dwdot[H]/d[H] */
    J[34] -= dqdci;               /* dwdot[OH]/d[H] */
    J[39] -= dqdci;               /* dwdot[CO]/d[H] */
    J[42] += dqdci;               /* dwdot[CO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[9];
    J[121] += dqdci;              /* dwdot[H]/d[OH] */
    J[124] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[129] -= dqdci;              /* dwdot[CO]/d[OH] */
    J[132] += dqdci;              /* dwdot[CO2]/d[OH] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[4];
    J[271] += dqdci;              /* dwdot[H]/d[CO] */
    J[274] -= dqdci;              /* dwdot[OH]/d[CO] */
    J[279] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[282] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[1];
    J[361] += dqdci;              /* dwdot[H]/d[CO2] */
    J[364] -= dqdci;              /* dwdot[OH]/d[CO2] */
    J[369] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[372] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[871] += dqdT;               /* dwdot[H]/dT */
    J[874] -= dqdT;               /* dwdot[OH]/dT */
    J[879] -= dqdT;               /* dwdot[CO]/dT */
    J[882] += dqdT;               /* dwdot[CO2]/dT */

    /*reaction 86: CO + HO2 <=> CO2 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[9];
    k_f = prefactor_units[85] * fwd_A[85]
                * exp(fwd_beta[85] * tc[0] - activation_units[85] * fwd_Ea[85] * invT);
    dlnkfdT = fwd_beta[85] * invT + activation_units[85] * fwd_Ea[85] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[12];
    Kc = exp(-g_RT[4] + g_RT[8] + g_RT[9] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[9]) + (h_RT[4] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] += q; /* OH */
    wdot[8] -= q; /* HO2 */
    wdot[9] -= q; /* CO */
    wdot[12] += q; /* CO2 */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[12];
    J[124] += dqdci;              /* dwdot[OH]/d[OH] */
    J[128] -= dqdci;              /* dwdot[HO2]/d[OH] */
    J[129] -= dqdci;              /* dwdot[CO]/d[OH] */
    J[132] += dqdci;              /* dwdot[CO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[9];
    J[244] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[248] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[249] -= dqdci;              /* dwdot[CO]/d[HO2] */
    J[252] += dqdci;              /* dwdot[CO2]/d[HO2] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[8];
    J[274] += dqdci;              /* dwdot[OH]/d[CO] */
    J[278] -= dqdci;              /* dwdot[HO2]/d[CO] */
    J[279] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[282] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[4];
    J[364] += dqdci;              /* dwdot[OH]/d[CO2] */
    J[368] -= dqdci;              /* dwdot[HO2]/d[CO2] */
    J[369] -= dqdci;              /* dwdot[CO]/d[CO2] */
    J[372] += dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[874] += dqdT;               /* dwdot[OH]/dT */
    J[878] -= dqdT;               /* dwdot[HO2]/dT */
    J[879] -= dqdT;               /* dwdot[CO]/dT */
    J[882] += dqdT;               /* dwdot[CO2]/dT */

    /*reaction 87: HCO + H <=> CO + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[20];
    k_f = prefactor_units[86] * fwd_A[86]
                * exp(fwd_beta[86] * tc[0] - activation_units[86] * fwd_Ea[86] * invT);
    dlnkfdT = fwd_beta[86] * invT + activation_units[86] * fwd_Ea[86] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[9];
    Kc = exp(g_RT[1] - g_RT[2] - g_RT[9] + g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[20]) + (h_RT[2] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[2] += q; /* H2 */
    wdot[9] += q; /* CO */
    wdot[20] -= q; /* HCO */
    /* d()/d[H] */
    dqdci =  + k_f*sc[20];
    J[31] -= dqdci;               /* dwdot[H]/d[H] */
    J[32] += dqdci;               /* dwdot[H2]/d[H] */
    J[39] += dqdci;               /* dwdot[CO]/d[H] */
    J[50] -= dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[9];
    J[61] -= dqdci;               /* dwdot[H]/d[H2] */
    J[62] += dqdci;               /* dwdot[H2]/d[H2] */
    J[69] += dqdci;               /* dwdot[CO]/d[H2] */
    J[80] -= dqdci;               /* dwdot[HCO]/d[H2] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[2];
    J[271] -= dqdci;              /* dwdot[H]/d[CO] */
    J[272] += dqdci;              /* dwdot[H2]/d[CO] */
    J[279] += dqdci;              /* dwdot[CO]/d[CO] */
    J[290] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[1];
    J[601] -= dqdci;              /* dwdot[H]/d[HCO] */
    J[602] += dqdci;              /* dwdot[H2]/d[HCO] */
    J[609] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[620] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[871] -= dqdT;               /* dwdot[H]/dT */
    J[872] += dqdT;               /* dwdot[H2]/dT */
    J[879] += dqdT;               /* dwdot[CO]/dT */
    J[890] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 88: CH3 + HCO <=> CH3CHO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[10]*sc[20];
    k_f = prefactor_units[87] * fwd_A[87]
                * exp(fwd_beta[87] * tc[0] - activation_units[87] * fwd_Ea[87] * invT);
    dlnkfdT = fwd_beta[87] * invT + activation_units[87] * fwd_Ea[87] * invT2;
    /* reverse */
    phi_r = sc[26];
    Kc = refCinv * exp(g_RT[10] + g_RT[20] - g_RT[26]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[10] + h_RT[20]) + (h_RT[26]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[10] -= q; /* CH3 */
    wdot[20] -= q; /* HCO */
    wdot[26] += q; /* CH3CHO */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[20];
    J[310] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[320] -= dqdci;              /* dwdot[HCO]/d[CH3] */
    J[326] += dqdci;              /* dwdot[CH3CHO]/d[CH3] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[10];
    J[610] -= dqdci;              /* dwdot[CH3]/d[HCO] */
    J[620] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    J[626] += dqdci;              /* dwdot[CH3CHO]/d[HCO] */
    /* d()/d[CH3CHO] */
    dqdci =  - k_r;
    J[790] -= dqdci;              /* dwdot[CH3]/d[CH3CHO] */
    J[800] -= dqdci;              /* dwdot[HCO]/d[CH3CHO] */
    J[806] += dqdci;              /* dwdot[CH3CHO]/d[CH3CHO] */
    /* d()/dT */
    J[880] -= dqdT;               /* dwdot[CH3]/dT */
    J[890] -= dqdT;               /* dwdot[HCO]/dT */
    J[896] += dqdT;               /* dwdot[CH3CHO]/dT */

    /*reaction 89: HCO + H2O <=> CO + H + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[7]*sc[20];
    k_f = prefactor_units[88] * fwd_A[88]
                * exp(fwd_beta[88] * tc[0] - activation_units[88] * fwd_Ea[88] * invT);
    dlnkfdT = fwd_beta[88] * invT + activation_units[88] * fwd_Ea[88] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[7]*sc[9];
    Kc = refC * exp(-g_RT[1] + g_RT[7] - g_RT[7] - g_RT[9] + g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7] + h_RT[20]) + (h_RT[1] + h_RT[7] + h_RT[9]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[9] += q; /* CO */
    wdot[20] -= q; /* HCO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[7]*sc[9];
    J[31] += dqdci;               /* dwdot[H]/d[H] */
    J[39] += dqdci;               /* dwdot[CO]/d[H] */
    J[50] -= dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[20] - k_r*sc[1]*sc[9];
    J[211] += dqdci;              /* dwdot[H]/d[H2O] */
    J[219] += dqdci;              /* dwdot[CO]/d[H2O] */
    J[230] -= dqdci;              /* dwdot[HCO]/d[H2O] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[1]*sc[7];
    J[271] += dqdci;              /* dwdot[H]/d[CO] */
    J[279] += dqdci;              /* dwdot[CO]/d[CO] */
    J[290] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[7];
    J[601] += dqdci;              /* dwdot[H]/d[HCO] */
    J[609] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[620] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[871] += dqdT;               /* dwdot[H]/dT */
    J[879] += dqdT;               /* dwdot[CO]/dT */
    J[890] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 90: HCO + O <=> CO + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[20];
    k_f = prefactor_units[89] * fwd_A[89]
                * exp(fwd_beta[89] * tc[0] - activation_units[89] * fwd_Ea[89] * invT);
    dlnkfdT = fwd_beta[89] * invT + activation_units[89] * fwd_Ea[89] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[9];
    Kc = exp(g_RT[3] - g_RT[4] - g_RT[9] + g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[20]) + (h_RT[4] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O */
    wdot[4] += q; /* OH */
    wdot[9] += q; /* CO */
    wdot[20] -= q; /* HCO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[20];
    J[93] -= dqdci;               /* dwdot[O]/d[O] */
    J[94] += dqdci;               /* dwdot[OH]/d[O] */
    J[99] += dqdci;               /* dwdot[CO]/d[O] */
    J[110] -= dqdci;              /* dwdot[HCO]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[9];
    J[123] -= dqdci;              /* dwdot[O]/d[OH] */
    J[124] += dqdci;              /* dwdot[OH]/d[OH] */
    J[129] += dqdci;              /* dwdot[CO]/d[OH] */
    J[140] -= dqdci;              /* dwdot[HCO]/d[OH] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[4];
    J[273] -= dqdci;              /* dwdot[O]/d[CO] */
    J[274] += dqdci;              /* dwdot[OH]/d[CO] */
    J[279] += dqdci;              /* dwdot[CO]/d[CO] */
    J[290] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[3];
    J[603] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[604] += dqdci;              /* dwdot[OH]/d[HCO] */
    J[609] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[620] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[873] -= dqdT;               /* dwdot[O]/dT */
    J[874] += dqdT;               /* dwdot[OH]/dT */
    J[879] += dqdT;               /* dwdot[CO]/dT */
    J[890] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 91: HCO + OH <=> CO + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[20];
    k_f = prefactor_units[90] * fwd_A[90]
                * exp(fwd_beta[90] * tc[0] - activation_units[90] * fwd_Ea[90] * invT);
    dlnkfdT = fwd_beta[90] * invT + activation_units[90] * fwd_Ea[90] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[9];
    Kc = exp(g_RT[4] - g_RT[7] - g_RT[9] + g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[20]) + (h_RT[7] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[7] += q; /* H2O */
    wdot[9] += q; /* CO */
    wdot[20] -= q; /* HCO */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[20];
    J[124] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[127] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[129] += dqdci;              /* dwdot[CO]/d[OH] */
    J[140] -= dqdci;              /* dwdot[HCO]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[9];
    J[214] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[217] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[219] += dqdci;              /* dwdot[CO]/d[H2O] */
    J[230] -= dqdci;              /* dwdot[HCO]/d[H2O] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[7];
    J[274] -= dqdci;              /* dwdot[OH]/d[CO] */
    J[277] += dqdci;              /* dwdot[H2O]/d[CO] */
    J[279] += dqdci;              /* dwdot[CO]/d[CO] */
    J[290] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[4];
    J[604] -= dqdci;              /* dwdot[OH]/d[HCO] */
    J[607] += dqdci;              /* dwdot[H2O]/d[HCO] */
    J[609] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[620] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[874] -= dqdT;               /* dwdot[OH]/dT */
    J[877] += dqdT;               /* dwdot[H2O]/dT */
    J[879] += dqdT;               /* dwdot[CO]/dT */
    J[890] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 92: CH3 + HCO <=> CH4 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[10]*sc[20];
    k_f = prefactor_units[91] * fwd_A[91]
                * exp(fwd_beta[91] * tc[0] - activation_units[91] * fwd_Ea[91] * invT);
    dlnkfdT = fwd_beta[91] * invT + activation_units[91] * fwd_Ea[91] * invT2;
    /* reverse */
    phi_r = sc[9]*sc[13];
    Kc = exp(-g_RT[9] + g_RT[10] - g_RT[13] + g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[10] + h_RT[20]) + (h_RT[9] + h_RT[13]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[9] += q; /* CO */
    wdot[10] -= q; /* CH3 */
    wdot[13] += q; /* CH4 */
    wdot[20] -= q; /* HCO */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[13];
    J[279] += dqdci;              /* dwdot[CO]/d[CO] */
    J[280] -= dqdci;              /* dwdot[CH3]/d[CO] */
    J[283] += dqdci;              /* dwdot[CH4]/d[CO] */
    J[290] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[20];
    J[309] += dqdci;              /* dwdot[CO]/d[CH3] */
    J[310] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[313] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[320] -= dqdci;              /* dwdot[HCO]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[9];
    J[399] += dqdci;              /* dwdot[CO]/d[CH4] */
    J[400] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[403] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[410] -= dqdci;              /* dwdot[HCO]/d[CH4] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[10];
    J[609] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[610] -= dqdci;              /* dwdot[CH3]/d[HCO] */
    J[613] += dqdci;              /* dwdot[CH4]/d[HCO] */
    J[620] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[879] += dqdT;               /* dwdot[CO]/dT */
    J[880] -= dqdT;               /* dwdot[CH3]/dT */
    J[883] += dqdT;               /* dwdot[CH4]/dT */
    J[890] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 93: HCO + O <=> CO2 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[20];
    k_f = prefactor_units[92] * fwd_A[92]
                * exp(fwd_beta[92] * tc[0] - activation_units[92] * fwd_Ea[92] * invT);
    dlnkfdT = fwd_beta[92] * invT + activation_units[92] * fwd_Ea[92] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[12];
    Kc = exp(-g_RT[1] + g_RT[3] - g_RT[12] + g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[20]) + (h_RT[1] + h_RT[12]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[3] -= q; /* O */
    wdot[12] += q; /* CO2 */
    wdot[20] -= q; /* HCO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[12];
    J[31] += dqdci;               /* dwdot[H]/d[H] */
    J[33] -= dqdci;               /* dwdot[O]/d[H] */
    J[42] += dqdci;               /* dwdot[CO2]/d[H] */
    J[50] -= dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[20];
    J[91] += dqdci;               /* dwdot[H]/d[O] */
    J[93] -= dqdci;               /* dwdot[O]/d[O] */
    J[102] += dqdci;              /* dwdot[CO2]/d[O] */
    J[110] -= dqdci;              /* dwdot[HCO]/d[O] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[1];
    J[361] += dqdci;              /* dwdot[H]/d[CO2] */
    J[363] -= dqdci;              /* dwdot[O]/d[CO2] */
    J[372] += dqdci;              /* dwdot[CO2]/d[CO2] */
    J[380] -= dqdci;              /* dwdot[HCO]/d[CO2] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[3];
    J[601] += dqdci;              /* dwdot[H]/d[HCO] */
    J[603] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[612] += dqdci;              /* dwdot[CO2]/d[HCO] */
    J[620] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[871] += dqdT;               /* dwdot[H]/dT */
    J[873] -= dqdT;               /* dwdot[O]/dT */
    J[882] += dqdT;               /* dwdot[CO2]/dT */
    J[890] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 94: HCO + O2 <=> CO + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[20];
    k_f = prefactor_units[93] * fwd_A[93]
                * exp(fwd_beta[93] * tc[0] - activation_units[93] * fwd_Ea[93] * invT);
    dlnkfdT = fwd_beta[93] * invT + activation_units[93] * fwd_Ea[93] * invT2;
    /* reverse */
    phi_r = sc[8]*sc[9];
    Kc = exp(g_RT[5] - g_RT[8] - g_RT[9] + g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[20]) + (h_RT[8] + h_RT[9]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[8] += q; /* HO2 */
    wdot[9] += q; /* CO */
    wdot[20] -= q; /* HCO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[20];
    J[155] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[158] += dqdci;              /* dwdot[HO2]/d[O2] */
    J[159] += dqdci;              /* dwdot[CO]/d[O2] */
    J[170] -= dqdci;              /* dwdot[HCO]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[9];
    J[245] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[248] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[249] += dqdci;              /* dwdot[CO]/d[HO2] */
    J[260] -= dqdci;              /* dwdot[HCO]/d[HO2] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[8];
    J[275] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[278] += dqdci;              /* dwdot[HO2]/d[CO] */
    J[279] += dqdci;              /* dwdot[CO]/d[CO] */
    J[290] -= dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[5];
    J[605] -= dqdci;              /* dwdot[O2]/d[HCO] */
    J[608] += dqdci;              /* dwdot[HO2]/d[HCO] */
    J[609] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[620] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[875] -= dqdT;               /* dwdot[O2]/dT */
    J[878] += dqdT;               /* dwdot[HO2]/dT */
    J[879] += dqdT;               /* dwdot[CO]/dT */
    J[890] -= dqdT;               /* dwdot[HCO]/dT */

    /*reaction 95: CH2O + H <=> HCO + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[11];
    k_f = prefactor_units[94] * fwd_A[94]
                * exp(fwd_beta[94] * tc[0] - activation_units[94] * fwd_Ea[94] * invT);
    dlnkfdT = fwd_beta[94] * invT + activation_units[94] * fwd_Ea[94] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[20];
    Kc = exp(g_RT[1] - g_RT[2] + g_RT[11] - g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[11]) + (h_RT[2] + h_RT[20]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[2] += q; /* H2 */
    wdot[11] -= q; /* CH2O */
    wdot[20] += q; /* HCO */
    /* d()/d[H] */
    dqdci =  + k_f*sc[11];
    J[31] -= dqdci;               /* dwdot[H]/d[H] */
    J[32] += dqdci;               /* dwdot[H2]/d[H] */
    J[41] -= dqdci;               /* dwdot[CH2O]/d[H] */
    J[50] += dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[20];
    J[61] -= dqdci;               /* dwdot[H]/d[H2] */
    J[62] += dqdci;               /* dwdot[H2]/d[H2] */
    J[71] -= dqdci;               /* dwdot[CH2O]/d[H2] */
    J[80] += dqdci;               /* dwdot[HCO]/d[H2] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[1];
    J[331] -= dqdci;              /* dwdot[H]/d[CH2O] */
    J[332] += dqdci;              /* dwdot[H2]/d[CH2O] */
    J[341] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[350] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[2];
    J[601] -= dqdci;              /* dwdot[H]/d[HCO] */
    J[602] += dqdci;              /* dwdot[H2]/d[HCO] */
    J[611] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    J[620] += dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[871] -= dqdT;               /* dwdot[H]/dT */
    J[872] += dqdT;               /* dwdot[H2]/dT */
    J[881] -= dqdT;               /* dwdot[CH2O]/dT */
    J[890] += dqdT;               /* dwdot[HCO]/dT */

    /*reaction 96: CH2O + O <=> HCO + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[11];
    k_f = prefactor_units[95] * fwd_A[95]
                * exp(fwd_beta[95] * tc[0] - activation_units[95] * fwd_Ea[95] * invT);
    dlnkfdT = fwd_beta[95] * invT + activation_units[95] * fwd_Ea[95] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[20];
    Kc = exp(g_RT[3] - g_RT[4] + g_RT[11] - g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[11]) + (h_RT[4] + h_RT[20]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O */
    wdot[4] += q; /* OH */
    wdot[11] -= q; /* CH2O */
    wdot[20] += q; /* HCO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[11];
    J[93] -= dqdci;               /* dwdot[O]/d[O] */
    J[94] += dqdci;               /* dwdot[OH]/d[O] */
    J[101] -= dqdci;              /* dwdot[CH2O]/d[O] */
    J[110] += dqdci;              /* dwdot[HCO]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[20];
    J[123] -= dqdci;              /* dwdot[O]/d[OH] */
    J[124] += dqdci;              /* dwdot[OH]/d[OH] */
    J[131] -= dqdci;              /* dwdot[CH2O]/d[OH] */
    J[140] += dqdci;              /* dwdot[HCO]/d[OH] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[3];
    J[333] -= dqdci;              /* dwdot[O]/d[CH2O] */
    J[334] += dqdci;              /* dwdot[OH]/d[CH2O] */
    J[341] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[350] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[4];
    J[603] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[604] += dqdci;              /* dwdot[OH]/d[HCO] */
    J[611] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    J[620] += dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[873] -= dqdT;               /* dwdot[O]/dT */
    J[874] += dqdT;               /* dwdot[OH]/dT */
    J[881] -= dqdT;               /* dwdot[CH2O]/dT */
    J[890] += dqdT;               /* dwdot[HCO]/dT */

    /*reaction 97: CH3 + CH2O <=> CH4 + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[10]*sc[11];
    k_f = prefactor_units[96] * fwd_A[96]
                * exp(fwd_beta[96] * tc[0] - activation_units[96] * fwd_Ea[96] * invT);
    dlnkfdT = fwd_beta[96] * invT + activation_units[96] * fwd_Ea[96] * invT2;
    /* reverse */
    phi_r = sc[13]*sc[20];
    Kc = exp(g_RT[10] + g_RT[11] - g_RT[13] - g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[10] + h_RT[11]) + (h_RT[13] + h_RT[20]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[10] -= q; /* CH3 */
    wdot[11] -= q; /* CH2O */
    wdot[13] += q; /* CH4 */
    wdot[20] += q; /* HCO */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[11];
    J[310] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[311] -= dqdci;              /* dwdot[CH2O]/d[CH3] */
    J[313] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[320] += dqdci;              /* dwdot[HCO]/d[CH3] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[10];
    J[340] -= dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[341] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[343] += dqdci;              /* dwdot[CH4]/d[CH2O] */
    J[350] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[20];
    J[400] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[401] -= dqdci;              /* dwdot[CH2O]/d[CH4] */
    J[403] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[410] += dqdci;              /* dwdot[HCO]/d[CH4] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[13];
    J[610] -= dqdci;              /* dwdot[CH3]/d[HCO] */
    J[611] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    J[613] += dqdci;              /* dwdot[CH4]/d[HCO] */
    J[620] += dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[880] -= dqdT;               /* dwdot[CH3]/dT */
    J[881] -= dqdT;               /* dwdot[CH2O]/dT */
    J[883] += dqdT;               /* dwdot[CH4]/dT */
    J[890] += dqdT;               /* dwdot[HCO]/dT */

    /*reaction 98: CH2O + OH <=> HCO + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[11];
    k_f = prefactor_units[97] * fwd_A[97]
                * exp(fwd_beta[97] * tc[0] - activation_units[97] * fwd_Ea[97] * invT);
    dlnkfdT = fwd_beta[97] * invT + activation_units[97] * fwd_Ea[97] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[20];
    Kc = exp(g_RT[4] - g_RT[7] + g_RT[11] - g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[11]) + (h_RT[7] + h_RT[20]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[7] += q; /* H2O */
    wdot[11] -= q; /* CH2O */
    wdot[20] += q; /* HCO */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[11];
    J[124] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[127] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[131] -= dqdci;              /* dwdot[CH2O]/d[OH] */
    J[140] += dqdci;              /* dwdot[HCO]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[20];
    J[214] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[217] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[221] -= dqdci;              /* dwdot[CH2O]/d[H2O] */
    J[230] += dqdci;              /* dwdot[HCO]/d[H2O] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[4];
    J[334] -= dqdci;              /* dwdot[OH]/d[CH2O] */
    J[337] += dqdci;              /* dwdot[H2O]/d[CH2O] */
    J[341] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[350] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[7];
    J[604] -= dqdci;              /* dwdot[OH]/d[HCO] */
    J[607] += dqdci;              /* dwdot[H2O]/d[HCO] */
    J[611] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    J[620] += dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[874] -= dqdT;               /* dwdot[OH]/dT */
    J[877] += dqdT;               /* dwdot[H2O]/dT */
    J[881] -= dqdT;               /* dwdot[CH2O]/dT */
    J[890] += dqdT;               /* dwdot[HCO]/dT */

    /*reaction 99: CH2O + CH <=> CH2CO + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[19];
    k_f = prefactor_units[98] * fwd_A[98]
                * exp(fwd_beta[98] * tc[0] - activation_units[98] * fwd_Ea[98] * invT);
    dlnkfdT = fwd_beta[98] * invT + activation_units[98] * fwd_Ea[98] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[16];
    Kc = exp(-g_RT[1] + g_RT[11] - g_RT[16] + g_RT[19]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[11] + h_RT[19]) + (h_RT[1] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[11] -= q; /* CH2O */
    wdot[16] += q; /* CH2CO */
    wdot[19] -= q; /* CH */
    /* d()/d[H] */
    dqdci =  - k_r*sc[16];
    J[31] += dqdci;               /* dwdot[H]/d[H] */
    J[41] -= dqdci;               /* dwdot[CH2O]/d[H] */
    J[46] += dqdci;               /* dwdot[CH2CO]/d[H] */
    J[49] -= dqdci;               /* dwdot[CH]/d[H] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[19];
    J[331] += dqdci;              /* dwdot[H]/d[CH2O] */
    J[341] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[346] += dqdci;              /* dwdot[CH2CO]/d[CH2O] */
    J[349] -= dqdci;              /* dwdot[CH]/d[CH2O] */
    /* d()/d[CH2CO] */
    dqdci =  - k_r*sc[1];
    J[481] += dqdci;              /* dwdot[H]/d[CH2CO] */
    J[491] -= dqdci;              /* dwdot[CH2O]/d[CH2CO] */
    J[496] += dqdci;              /* dwdot[CH2CO]/d[CH2CO] */
    J[499] -= dqdci;              /* dwdot[CH]/d[CH2CO] */
    /* d()/d[CH] */
    dqdci =  + k_f*sc[11];
    J[571] += dqdci;              /* dwdot[H]/d[CH] */
    J[581] -= dqdci;              /* dwdot[CH2O]/d[CH] */
    J[586] += dqdci;              /* dwdot[CH2CO]/d[CH] */
    J[589] -= dqdci;              /* dwdot[CH]/d[CH] */
    /* d()/dT */
    J[871] += dqdT;               /* dwdot[H]/dT */
    J[881] -= dqdT;               /* dwdot[CH2O]/dT */
    J[886] += dqdT;               /* dwdot[CH2CO]/dT */
    J[889] -= dqdT;               /* dwdot[CH]/dT */

    /*reaction 100: CH2O + O2 <=> HCO + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[11];
    k_f = prefactor_units[99] * fwd_A[99]
                * exp(fwd_beta[99] * tc[0] - activation_units[99] * fwd_Ea[99] * invT);
    dlnkfdT = fwd_beta[99] * invT + activation_units[99] * fwd_Ea[99] * invT2;
    /* reverse */
    phi_r = sc[8]*sc[20];
    Kc = exp(g_RT[5] - g_RT[8] + g_RT[11] - g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[11]) + (h_RT[8] + h_RT[20]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[8] += q; /* HO2 */
    wdot[11] -= q; /* CH2O */
    wdot[20] += q; /* HCO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[11];
    J[155] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[158] += dqdci;              /* dwdot[HO2]/d[O2] */
    J[161] -= dqdci;              /* dwdot[CH2O]/d[O2] */
    J[170] += dqdci;              /* dwdot[HCO]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[20];
    J[245] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[248] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[251] -= dqdci;              /* dwdot[CH2O]/d[HO2] */
    J[260] += dqdci;              /* dwdot[HCO]/d[HO2] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[5];
    J[335] -= dqdci;              /* dwdot[O2]/d[CH2O] */
    J[338] += dqdci;              /* dwdot[HO2]/d[CH2O] */
    J[341] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[350] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[8];
    J[605] -= dqdci;              /* dwdot[O2]/d[HCO] */
    J[608] += dqdci;              /* dwdot[HO2]/d[HCO] */
    J[611] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    J[620] += dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[875] -= dqdT;               /* dwdot[O2]/dT */
    J[878] += dqdT;               /* dwdot[HO2]/dT */
    J[881] -= dqdT;               /* dwdot[CH2O]/dT */
    J[890] += dqdT;               /* dwdot[HCO]/dT */

    /*reaction 101: CH2O + HO2 <=> HCO + H2O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[11];
    k_f = prefactor_units[100] * fwd_A[100]
                * exp(fwd_beta[100] * tc[0] - activation_units[100] * fwd_Ea[100] * invT);
    dlnkfdT = fwd_beta[100] * invT + activation_units[100] * fwd_Ea[100] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[20];
    Kc = exp(-g_RT[6] + g_RT[8] + g_RT[11] - g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[11]) + (h_RT[6] + h_RT[20]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[6] += q; /* H2O2 */
    wdot[8] -= q; /* HO2 */
    wdot[11] -= q; /* CH2O */
    wdot[20] += q; /* HCO */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[20];
    J[186] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
    J[188] -= dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[191] -= dqdci;              /* dwdot[CH2O]/d[H2O2] */
    J[200] += dqdci;              /* dwdot[HCO]/d[H2O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[11];
    J[246] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    J[248] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[251] -= dqdci;              /* dwdot[CH2O]/d[HO2] */
    J[260] += dqdci;              /* dwdot[HCO]/d[HO2] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[8];
    J[336] += dqdci;              /* dwdot[H2O2]/d[CH2O] */
    J[338] -= dqdci;              /* dwdot[HO2]/d[CH2O] */
    J[341] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[350] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[6];
    J[606] += dqdci;              /* dwdot[H2O2]/d[HCO] */
    J[608] -= dqdci;              /* dwdot[HO2]/d[HCO] */
    J[611] -= dqdci;              /* dwdot[CH2O]/d[HCO] */
    J[620] += dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[876] += dqdT;               /* dwdot[H2O2]/dT */
    J[878] -= dqdT;               /* dwdot[HO2]/dT */
    J[881] -= dqdT;               /* dwdot[CH2O]/dT */
    J[890] += dqdT;               /* dwdot[HCO]/dT */

    /*reaction 102: 2.000000 H + CO2 <=> H2 + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[1], 2.000000)*sc[12];
    k_f = prefactor_units[101] * fwd_A[101]
                * exp(fwd_beta[101] * tc[0] - activation_units[101] * fwd_Ea[101] * invT);
    dlnkfdT = fwd_beta[101] * invT + activation_units[101] * fwd_Ea[101] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[12];
    Kc = refCinv * exp(2.000000*g_RT[1] - g_RT[2] + g_RT[12] - g_RT[12]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[1] + h_RT[12]) + (h_RT[2] + h_RT[12]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= 2 * q; /* H */
    wdot[2] += q; /* H2 */
    /* d()/d[H] */
    dqdci =  + k_f*2.000000*sc[1]*sc[12];
    J[31] += -2 * dqdci;          /* dwdot[H]/d[H] */
    J[32] += dqdci;               /* dwdot[H2]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[12];
    J[61] += -2 * dqdci;          /* dwdot[H]/d[H2] */
    J[62] += dqdci;               /* dwdot[H2]/d[H2] */
    /* d()/d[CO2] */
    dqdci =  + k_f*pow(sc[1], 2.000000) - k_r*sc[2];
    J[361] += -2 * dqdci;         /* dwdot[H]/d[CO2] */
    J[362] += dqdci;              /* dwdot[H2]/d[CO2] */
    /* d()/dT */
    J[871] += -2 * dqdT;          /* dwdot[H]/dT */
    J[872] += dqdT;               /* dwdot[H2]/dT */

    /*reaction 103: SXCH2 + CO2 <=> TXCH2 + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[12]*sc[22];
    k_f = prefactor_units[102] * fwd_A[102]
                * exp(fwd_beta[102] * tc[0] - activation_units[102] * fwd_Ea[102] * invT);
    dlnkfdT = fwd_beta[102] * invT + activation_units[102] * fwd_Ea[102] * invT2;
    /* reverse */
    phi_r = sc[12]*sc[21];
    Kc = exp(g_RT[12] - g_RT[12] - g_RT[21] + g_RT[22]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[12] + h_RT[22]) + (h_RT[12] + h_RT[21]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[21] += q; /* TXCH2 */
    wdot[22] -= q; /* SXCH2 */
    /* d()/d[CO2] */
    dqdci =  + k_f*sc[22] - k_r*sc[21];
    J[381] += dqdci;              /* dwdot[TXCH2]/d[CO2] */
    J[382] -= dqdci;              /* dwdot[SXCH2]/d[CO2] */
    /* d()/d[TXCH2] */
    dqdci =  - k_r*sc[12];
    J[651] += dqdci;              /* dwdot[TXCH2]/d[TXCH2] */
    J[652] -= dqdci;              /* dwdot[SXCH2]/d[TXCH2] */
    /* d()/d[SXCH2] */
    dqdci =  + k_f*sc[12];
    J[681] += dqdci;              /* dwdot[TXCH2]/d[SXCH2] */
    J[682] -= dqdci;              /* dwdot[SXCH2]/d[SXCH2] */
    /* d()/dT */
    J[891] += dqdT;               /* dwdot[TXCH2]/dT */
    J[892] -= dqdT;               /* dwdot[SXCH2]/dT */

    /*reaction 104: SXCH2 + CO2 <=> CH2O + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[12]*sc[22];
    k_f = prefactor_units[103] * fwd_A[103]
                * exp(fwd_beta[103] * tc[0] - activation_units[103] * fwd_Ea[103] * invT);
    dlnkfdT = fwd_beta[103] * invT + activation_units[103] * fwd_Ea[103] * invT2;
    /* reverse */
    phi_r = sc[9]*sc[11];
    Kc = exp(-g_RT[9] - g_RT[11] + g_RT[12] + g_RT[22]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[12] + h_RT[22]) + (h_RT[9] + h_RT[11]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[9] += q; /* CO */
    wdot[11] += q; /* CH2O */
    wdot[12] -= q; /* CO2 */
    wdot[22] -= q; /* SXCH2 */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[11];
    J[279] += dqdci;              /* dwdot[CO]/d[CO] */
    J[281] += dqdci;              /* dwdot[CH2O]/d[CO] */
    J[282] -= dqdci;              /* dwdot[CO2]/d[CO] */
    J[292] -= dqdci;              /* dwdot[SXCH2]/d[CO] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[9];
    J[339] += dqdci;              /* dwdot[CO]/d[CH2O] */
    J[341] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[342] -= dqdci;              /* dwdot[CO2]/d[CH2O] */
    J[352] -= dqdci;              /* dwdot[SXCH2]/d[CH2O] */
    /* d()/d[CO2] */
    dqdci =  + k_f*sc[22];
    J[369] += dqdci;              /* dwdot[CO]/d[CO2] */
    J[371] += dqdci;              /* dwdot[CH2O]/d[CO2] */
    J[372] -= dqdci;              /* dwdot[CO2]/d[CO2] */
    J[382] -= dqdci;              /* dwdot[SXCH2]/d[CO2] */
    /* d()/d[SXCH2] */
    dqdci =  + k_f*sc[12];
    J[669] += dqdci;              /* dwdot[CO]/d[SXCH2] */
    J[671] += dqdci;              /* dwdot[CH2O]/d[SXCH2] */
    J[672] -= dqdci;              /* dwdot[CO2]/d[SXCH2] */
    J[682] -= dqdci;              /* dwdot[SXCH2]/d[SXCH2] */
    /* d()/dT */
    J[879] += dqdT;               /* dwdot[CO]/dT */
    J[881] += dqdT;               /* dwdot[CH2O]/dT */
    J[882] -= dqdT;               /* dwdot[CO2]/dT */
    J[892] -= dqdT;               /* dwdot[SXCH2]/dT */

    /*reaction 105: CH + CO2 <=> HCO + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[12]*sc[19];
    k_f = prefactor_units[104] * fwd_A[104]
                * exp(fwd_beta[104] * tc[0] - activation_units[104] * fwd_Ea[104] * invT);
    dlnkfdT = fwd_beta[104] * invT + activation_units[104] * fwd_Ea[104] * invT2;
    /* reverse */
    phi_r = sc[9]*sc[20];
    Kc = exp(-g_RT[9] + g_RT[12] + g_RT[19] - g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[12] + h_RT[19]) + (h_RT[9] + h_RT[20]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[9] += q; /* CO */
    wdot[12] -= q; /* CO2 */
    wdot[19] -= q; /* CH */
    wdot[20] += q; /* HCO */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[20];
    J[279] += dqdci;              /* dwdot[CO]/d[CO] */
    J[282] -= dqdci;              /* dwdot[CO2]/d[CO] */
    J[289] -= dqdci;              /* dwdot[CH]/d[CO] */
    J[290] += dqdci;              /* dwdot[HCO]/d[CO] */
    /* d()/d[CO2] */
    dqdci =  + k_f*sc[19];
    J[369] += dqdci;              /* dwdot[CO]/d[CO2] */
    J[372] -= dqdci;              /* dwdot[CO2]/d[CO2] */
    J[379] -= dqdci;              /* dwdot[CH]/d[CO2] */
    J[380] += dqdci;              /* dwdot[HCO]/d[CO2] */
    /* d()/d[CH] */
    dqdci =  + k_f*sc[12];
    J[579] += dqdci;              /* dwdot[CO]/d[CH] */
    J[582] -= dqdci;              /* dwdot[CO2]/d[CH] */
    J[589] -= dqdci;              /* dwdot[CH]/d[CH] */
    J[590] += dqdci;              /* dwdot[HCO]/d[CH] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[9];
    J[609] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[612] -= dqdci;              /* dwdot[CO2]/d[HCO] */
    J[619] -= dqdci;              /* dwdot[CH]/d[HCO] */
    J[620] += dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[879] += dqdT;               /* dwdot[CO]/dT */
    J[882] -= dqdT;               /* dwdot[CO2]/dT */
    J[889] -= dqdT;               /* dwdot[CH]/dT */
    J[890] += dqdT;               /* dwdot[HCO]/dT */

    /*reaction 106: C2H2 + O <=> TXCH2 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[14];
    k_f = prefactor_units[105] * fwd_A[105]
                * exp(fwd_beta[105] * tc[0] - activation_units[105] * fwd_Ea[105] * invT);
    dlnkfdT = fwd_beta[105] * invT + activation_units[105] * fwd_Ea[105] * invT2;
    /* reverse */
    phi_r = sc[9]*sc[21];
    Kc = exp(g_RT[3] - g_RT[9] + g_RT[14] - g_RT[21]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[14]) + (h_RT[9] + h_RT[21]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O */
    wdot[9] += q; /* CO */
    wdot[14] -= q; /* C2H2 */
    wdot[21] += q; /* TXCH2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[14];
    J[93] -= dqdci;               /* dwdot[O]/d[O] */
    J[99] += dqdci;               /* dwdot[CO]/d[O] */
    J[104] -= dqdci;              /* dwdot[C2H2]/d[O] */
    J[111] += dqdci;              /* dwdot[TXCH2]/d[O] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[21];
    J[273] -= dqdci;              /* dwdot[O]/d[CO] */
    J[279] += dqdci;              /* dwdot[CO]/d[CO] */
    J[284] -= dqdci;              /* dwdot[C2H2]/d[CO] */
    J[291] += dqdci;              /* dwdot[TXCH2]/d[CO] */
    /* d()/d[C2H2] */
    dqdci =  + k_f*sc[3];
    J[423] -= dqdci;              /* dwdot[O]/d[C2H2] */
    J[429] += dqdci;              /* dwdot[CO]/d[C2H2] */
    J[434] -= dqdci;              /* dwdot[C2H2]/d[C2H2] */
    J[441] += dqdci;              /* dwdot[TXCH2]/d[C2H2] */
    /* d()/d[TXCH2] */
    dqdci =  - k_r*sc[9];
    J[633] -= dqdci;              /* dwdot[O]/d[TXCH2] */
    J[639] += dqdci;              /* dwdot[CO]/d[TXCH2] */
    J[644] -= dqdci;              /* dwdot[C2H2]/d[TXCH2] */
    J[651] += dqdci;              /* dwdot[TXCH2]/d[TXCH2] */
    /* d()/dT */
    J[873] -= dqdT;               /* dwdot[O]/dT */
    J[879] += dqdT;               /* dwdot[CO]/dT */
    J[884] -= dqdT;               /* dwdot[C2H2]/dT */
    J[891] += dqdT;               /* dwdot[TXCH2]/dT */

    /*reaction 107: C2H2 + OH <=> CH3 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[14];
    k_f = prefactor_units[106] * fwd_A[106]
                * exp(fwd_beta[106] * tc[0] - activation_units[106] * fwd_Ea[106] * invT);
    dlnkfdT = fwd_beta[106] * invT + activation_units[106] * fwd_Ea[106] * invT2;
    /* reverse */
    phi_r = sc[9]*sc[10];
    Kc = exp(g_RT[4] - g_RT[9] - g_RT[10] + g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[14]) + (h_RT[9] + h_RT[10]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[9] += q; /* CO */
    wdot[10] += q; /* CH3 */
    wdot[14] -= q; /* C2H2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[14];
    J[124] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[129] += dqdci;              /* dwdot[CO]/d[OH] */
    J[130] += dqdci;              /* dwdot[CH3]/d[OH] */
    J[134] -= dqdci;              /* dwdot[C2H2]/d[OH] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[10];
    J[274] -= dqdci;              /* dwdot[OH]/d[CO] */
    J[279] += dqdci;              /* dwdot[CO]/d[CO] */
    J[280] += dqdci;              /* dwdot[CH3]/d[CO] */
    J[284] -= dqdci;              /* dwdot[C2H2]/d[CO] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[9];
    J[304] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[309] += dqdci;              /* dwdot[CO]/d[CH3] */
    J[310] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[314] -= dqdci;              /* dwdot[C2H2]/d[CH3] */
    /* d()/d[C2H2] */
    dqdci =  + k_f*sc[4];
    J[424] -= dqdci;              /* dwdot[OH]/d[C2H2] */
    J[429] += dqdci;              /* dwdot[CO]/d[C2H2] */
    J[430] += dqdci;              /* dwdot[CH3]/d[C2H2] */
    J[434] -= dqdci;              /* dwdot[C2H2]/d[C2H2] */
    /* d()/dT */
    J[874] -= dqdT;               /* dwdot[OH]/dT */
    J[879] += dqdT;               /* dwdot[CO]/dT */
    J[880] += dqdT;               /* dwdot[CH3]/dT */
    J[884] -= dqdT;               /* dwdot[C2H2]/dT */

    /*reaction 108: C2H2 + OH <=> CH2CO + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[14];
    k_f = prefactor_units[107] * fwd_A[107]
                * exp(fwd_beta[107] * tc[0] - activation_units[107] * fwd_Ea[107] * invT);
    dlnkfdT = fwd_beta[107] * invT + activation_units[107] * fwd_Ea[107] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[16];
    Kc = exp(-g_RT[1] + g_RT[4] + g_RT[14] - g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[14]) + (h_RT[1] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[4] -= q; /* OH */
    wdot[14] -= q; /* C2H2 */
    wdot[16] += q; /* CH2CO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[16];
    J[31] += dqdci;               /* dwdot[H]/d[H] */
    J[34] -= dqdci;               /* dwdot[OH]/d[H] */
    J[44] -= dqdci;               /* dwdot[C2H2]/d[H] */
    J[46] += dqdci;               /* dwdot[CH2CO]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[14];
    J[121] += dqdci;              /* dwdot[H]/d[OH] */
    J[124] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[134] -= dqdci;              /* dwdot[C2H2]/d[OH] */
    J[136] += dqdci;              /* dwdot[CH2CO]/d[OH] */
    /* d()/d[C2H2] */
    dqdci =  + k_f*sc[4];
    J[421] += dqdci;              /* dwdot[H]/d[C2H2] */
    J[424] -= dqdci;              /* dwdot[OH]/d[C2H2] */
    J[434] -= dqdci;              /* dwdot[C2H2]/d[C2H2] */
    J[436] += dqdci;              /* dwdot[CH2CO]/d[C2H2] */
    /* d()/d[CH2CO] */
    dqdci =  - k_r*sc[1];
    J[481] += dqdci;              /* dwdot[H]/d[CH2CO] */
    J[484] -= dqdci;              /* dwdot[OH]/d[CH2CO] */
    J[494] -= dqdci;              /* dwdot[C2H2]/d[CH2CO] */
    J[496] += dqdci;              /* dwdot[CH2CO]/d[CH2CO] */
    /* d()/dT */
    J[871] += dqdT;               /* dwdot[H]/dT */
    J[874] -= dqdT;               /* dwdot[OH]/dT */
    J[884] -= dqdT;               /* dwdot[C2H2]/dT */
    J[886] += dqdT;               /* dwdot[CH2CO]/dT */

    /*reaction 109: C2H2 + O <=> HCCO + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[14];
    k_f = prefactor_units[108] * fwd_A[108]
                * exp(fwd_beta[108] * tc[0] - activation_units[108] * fwd_Ea[108] * invT);
    dlnkfdT = fwd_beta[108] * invT + activation_units[108] * fwd_Ea[108] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[25];
    Kc = exp(-g_RT[1] + g_RT[3] + g_RT[14] - g_RT[25]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[14]) + (h_RT[1] + h_RT[25]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[3] -= q; /* O */
    wdot[14] -= q; /* C2H2 */
    wdot[25] += q; /* HCCO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[25];
    J[31] += dqdci;               /* dwdot[H]/d[H] */
    J[33] -= dqdci;               /* dwdot[O]/d[H] */
    J[44] -= dqdci;               /* dwdot[C2H2]/d[H] */
    J[55] += dqdci;               /* dwdot[HCCO]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[14];
    J[91] += dqdci;               /* dwdot[H]/d[O] */
    J[93] -= dqdci;               /* dwdot[O]/d[O] */
    J[104] -= dqdci;              /* dwdot[C2H2]/d[O] */
    J[115] += dqdci;              /* dwdot[HCCO]/d[O] */
    /* d()/d[C2H2] */
    dqdci =  + k_f*sc[3];
    J[421] += dqdci;              /* dwdot[H]/d[C2H2] */
    J[423] -= dqdci;              /* dwdot[O]/d[C2H2] */
    J[434] -= dqdci;              /* dwdot[C2H2]/d[C2H2] */
    J[445] += dqdci;              /* dwdot[HCCO]/d[C2H2] */
    /* d()/d[HCCO] */
    dqdci =  - k_r*sc[1];
    J[751] += dqdci;              /* dwdot[H]/d[HCCO] */
    J[753] -= dqdci;              /* dwdot[O]/d[HCCO] */
    J[764] -= dqdci;              /* dwdot[C2H2]/d[HCCO] */
    J[775] += dqdci;              /* dwdot[HCCO]/d[HCCO] */
    /* d()/dT */
    J[871] += dqdT;               /* dwdot[H]/dT */
    J[873] -= dqdT;               /* dwdot[O]/dT */
    J[884] -= dqdT;               /* dwdot[C2H2]/dT */
    J[895] += dqdT;               /* dwdot[HCCO]/dT */

    /*reaction 110: C2H3 + OH <=> C2H2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[23];
    k_f = prefactor_units[109] * fwd_A[109]
                * exp(fwd_beta[109] * tc[0] - activation_units[109] * fwd_Ea[109] * invT);
    dlnkfdT = fwd_beta[109] * invT + activation_units[109] * fwd_Ea[109] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[14];
    Kc = exp(g_RT[4] - g_RT[7] - g_RT[14] + g_RT[23]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[23]) + (h_RT[7] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[7] += q; /* H2O */
    wdot[14] += q; /* C2H2 */
    wdot[23] -= q; /* C2H3 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[23];
    J[124] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[127] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[134] += dqdci;              /* dwdot[C2H2]/d[OH] */
    J[143] -= dqdci;              /* dwdot[C2H3]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[14];
    J[214] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[217] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[224] += dqdci;              /* dwdot[C2H2]/d[H2O] */
    J[233] -= dqdci;              /* dwdot[C2H3]/d[H2O] */
    /* d()/d[C2H2] */
    dqdci =  - k_r*sc[7];
    J[424] -= dqdci;              /* dwdot[OH]/d[C2H2] */
    J[427] += dqdci;              /* dwdot[H2O]/d[C2H2] */
    J[434] += dqdci;              /* dwdot[C2H2]/d[C2H2] */
    J[443] -= dqdci;              /* dwdot[C2H3]/d[C2H2] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[4];
    J[694] -= dqdci;              /* dwdot[OH]/d[C2H3] */
    J[697] += dqdci;              /* dwdot[H2O]/d[C2H3] */
    J[704] += dqdci;              /* dwdot[C2H2]/d[C2H3] */
    J[713] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    /* d()/dT */
    J[874] -= dqdT;               /* dwdot[OH]/dT */
    J[877] += dqdT;               /* dwdot[H2O]/dT */
    J[884] += dqdT;               /* dwdot[C2H2]/dT */
    J[893] -= dqdT;               /* dwdot[C2H3]/dT */

    /*reaction 111: C2H3 + O2 <=> CH2CHO + O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[23];
    k_f = prefactor_units[110] * fwd_A[110]
                * exp(fwd_beta[110] * tc[0] - activation_units[110] * fwd_Ea[110] * invT);
    dlnkfdT = fwd_beta[110] * invT + activation_units[110] * fwd_Ea[110] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[27];
    Kc = exp(-g_RT[3] + g_RT[5] + g_RT[23] - g_RT[27]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[23]) + (h_RT[3] + h_RT[27]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += q; /* O */
    wdot[5] -= q; /* O2 */
    wdot[23] -= q; /* C2H3 */
    wdot[27] += q; /* CH2CHO */
    /* d()/d[O] */
    dqdci =  - k_r*sc[27];
    J[93] += dqdci;               /* dwdot[O]/d[O] */
    J[95] -= dqdci;               /* dwdot[O2]/d[O] */
    J[113] -= dqdci;              /* dwdot[C2H3]/d[O] */
    J[117] += dqdci;              /* dwdot[CH2CHO]/d[O] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[23];
    J[153] += dqdci;              /* dwdot[O]/d[O2] */
    J[155] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[173] -= dqdci;              /* dwdot[C2H3]/d[O2] */
    J[177] += dqdci;              /* dwdot[CH2CHO]/d[O2] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[5];
    J[693] += dqdci;              /* dwdot[O]/d[C2H3] */
    J[695] -= dqdci;              /* dwdot[O2]/d[C2H3] */
    J[713] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    J[717] += dqdci;              /* dwdot[CH2CHO]/d[C2H3] */
    /* d()/d[CH2CHO] */
    dqdci =  - k_r*sc[3];
    J[813] += dqdci;              /* dwdot[O]/d[CH2CHO] */
    J[815] -= dqdci;              /* dwdot[O2]/d[CH2CHO] */
    J[833] -= dqdci;              /* dwdot[C2H3]/d[CH2CHO] */
    J[837] += dqdci;              /* dwdot[CH2CHO]/d[CH2CHO] */
    /* d()/dT */
    J[873] += dqdT;               /* dwdot[O]/dT */
    J[875] -= dqdT;               /* dwdot[O2]/dT */
    J[893] -= dqdT;               /* dwdot[C2H3]/dT */
    J[897] += dqdT;               /* dwdot[CH2CHO]/dT */

    /*reaction 112: C2H3 + O <=> CH2CHO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[23];
    k_f = prefactor_units[111] * fwd_A[111]
                * exp(fwd_beta[111] * tc[0] - activation_units[111] * fwd_Ea[111] * invT);
    dlnkfdT = fwd_beta[111] * invT + activation_units[111] * fwd_Ea[111] * invT2;
    /* reverse */
    phi_r = sc[27];
    Kc = refCinv * exp(g_RT[3] + g_RT[23] - g_RT[27]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[23]) + (h_RT[27]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O */
    wdot[23] -= q; /* C2H3 */
    wdot[27] += q; /* CH2CHO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[23];
    J[93] -= dqdci;               /* dwdot[O]/d[O] */
    J[113] -= dqdci;              /* dwdot[C2H3]/d[O] */
    J[117] += dqdci;              /* dwdot[CH2CHO]/d[O] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[3];
    J[693] -= dqdci;              /* dwdot[O]/d[C2H3] */
    J[713] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    J[717] += dqdci;              /* dwdot[CH2CHO]/d[C2H3] */
    /* d()/d[CH2CHO] */
    dqdci =  - k_r;
    J[813] -= dqdci;              /* dwdot[O]/d[CH2CHO] */
    J[833] -= dqdci;              /* dwdot[C2H3]/d[CH2CHO] */
    J[837] += dqdci;              /* dwdot[CH2CHO]/d[CH2CHO] */
    /* d()/dT */
    J[873] -= dqdT;               /* dwdot[O]/dT */
    J[893] -= dqdT;               /* dwdot[C2H3]/dT */
    J[897] += dqdT;               /* dwdot[CH2CHO]/dT */

    /*reaction 113: C2H3 + H <=> C2H2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[23];
    k_f = prefactor_units[112] * fwd_A[112]
                * exp(fwd_beta[112] * tc[0] - activation_units[112] * fwd_Ea[112] * invT);
    dlnkfdT = fwd_beta[112] * invT + activation_units[112] * fwd_Ea[112] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[14];
    Kc = exp(g_RT[1] - g_RT[2] - g_RT[14] + g_RT[23]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[23]) + (h_RT[2] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[2] += q; /* H2 */
    wdot[14] += q; /* C2H2 */
    wdot[23] -= q; /* C2H3 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[23];
    J[31] -= dqdci;               /* dwdot[H]/d[H] */
    J[32] += dqdci;               /* dwdot[H2]/d[H] */
    J[44] += dqdci;               /* dwdot[C2H2]/d[H] */
    J[53] -= dqdci;               /* dwdot[C2H3]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[14];
    J[61] -= dqdci;               /* dwdot[H]/d[H2] */
    J[62] += dqdci;               /* dwdot[H2]/d[H2] */
    J[74] += dqdci;               /* dwdot[C2H2]/d[H2] */
    J[83] -= dqdci;               /* dwdot[C2H3]/d[H2] */
    /* d()/d[C2H2] */
    dqdci =  - k_r*sc[2];
    J[421] -= dqdci;              /* dwdot[H]/d[C2H2] */
    J[422] += dqdci;              /* dwdot[H2]/d[C2H2] */
    J[434] += dqdci;              /* dwdot[C2H2]/d[C2H2] */
    J[443] -= dqdci;              /* dwdot[C2H3]/d[C2H2] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[1];
    J[691] -= dqdci;              /* dwdot[H]/d[C2H3] */
    J[692] += dqdci;              /* dwdot[H2]/d[C2H3] */
    J[704] += dqdci;              /* dwdot[C2H2]/d[C2H3] */
    J[713] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    /* d()/dT */
    J[871] -= dqdT;               /* dwdot[H]/dT */
    J[872] += dqdT;               /* dwdot[H2]/dT */
    J[884] += dqdT;               /* dwdot[C2H2]/dT */
    J[893] -= dqdT;               /* dwdot[C2H3]/dT */

    /*reaction 114: C2H3 + CH3 <=> C2H2 + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[10]*sc[23];
    k_f = prefactor_units[113] * fwd_A[113]
                * exp(fwd_beta[113] * tc[0] - activation_units[113] * fwd_Ea[113] * invT);
    dlnkfdT = fwd_beta[113] * invT + activation_units[113] * fwd_Ea[113] * invT2;
    /* reverse */
    phi_r = sc[13]*sc[14];
    Kc = exp(g_RT[10] - g_RT[13] - g_RT[14] + g_RT[23]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[10] + h_RT[23]) + (h_RT[13] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[10] -= q; /* CH3 */
    wdot[13] += q; /* CH4 */
    wdot[14] += q; /* C2H2 */
    wdot[23] -= q; /* C2H3 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[23];
    J[310] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[313] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[314] += dqdci;              /* dwdot[C2H2]/d[CH3] */
    J[323] -= dqdci;              /* dwdot[C2H3]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[14];
    J[400] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[403] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[404] += dqdci;              /* dwdot[C2H2]/d[CH4] */
    J[413] -= dqdci;              /* dwdot[C2H3]/d[CH4] */
    /* d()/d[C2H2] */
    dqdci =  - k_r*sc[13];
    J[430] -= dqdci;              /* dwdot[CH3]/d[C2H2] */
    J[433] += dqdci;              /* dwdot[CH4]/d[C2H2] */
    J[434] += dqdci;              /* dwdot[C2H2]/d[C2H2] */
    J[443] -= dqdci;              /* dwdot[C2H3]/d[C2H2] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[10];
    J[700] -= dqdci;              /* dwdot[CH3]/d[C2H3] */
    J[703] += dqdci;              /* dwdot[CH4]/d[C2H3] */
    J[704] += dqdci;              /* dwdot[C2H2]/d[C2H3] */
    J[713] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    /* d()/dT */
    J[880] -= dqdT;               /* dwdot[CH3]/dT */
    J[883] += dqdT;               /* dwdot[CH4]/dT */
    J[884] += dqdT;               /* dwdot[C2H2]/dT */
    J[893] -= dqdT;               /* dwdot[C2H3]/dT */

    /*reaction 115: C2H3 + O2 <=> HCO + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[23];
    k_f = prefactor_units[114] * fwd_A[114]
                * exp(fwd_beta[114] * tc[0] - activation_units[114] * fwd_Ea[114] * invT);
    dlnkfdT = fwd_beta[114] * invT + activation_units[114] * fwd_Ea[114] * invT2;
    /* reverse */
    phi_r = sc[11]*sc[20];
    Kc = exp(g_RT[5] - g_RT[11] - g_RT[20] + g_RT[23]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[23]) + (h_RT[11] + h_RT[20]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[11] += q; /* CH2O */
    wdot[20] += q; /* HCO */
    wdot[23] -= q; /* C2H3 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[23];
    J[155] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[161] += dqdci;              /* dwdot[CH2O]/d[O2] */
    J[170] += dqdci;              /* dwdot[HCO]/d[O2] */
    J[173] -= dqdci;              /* dwdot[C2H3]/d[O2] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[20];
    J[335] -= dqdci;              /* dwdot[O2]/d[CH2O] */
    J[341] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[350] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[353] -= dqdci;              /* dwdot[C2H3]/d[CH2O] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[11];
    J[605] -= dqdci;              /* dwdot[O2]/d[HCO] */
    J[611] += dqdci;              /* dwdot[CH2O]/d[HCO] */
    J[620] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[623] -= dqdci;              /* dwdot[C2H3]/d[HCO] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[5];
    J[695] -= dqdci;              /* dwdot[O2]/d[C2H3] */
    J[701] += dqdci;              /* dwdot[CH2O]/d[C2H3] */
    J[710] += dqdci;              /* dwdot[HCO]/d[C2H3] */
    J[713] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    /* d()/dT */
    J[875] -= dqdT;               /* dwdot[O2]/dT */
    J[881] += dqdT;               /* dwdot[CH2O]/dT */
    J[890] += dqdT;               /* dwdot[HCO]/dT */
    J[893] -= dqdT;               /* dwdot[C2H3]/dT */

    /*reaction 116: C2H3 + H2O2 <=> C2H4 + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[23];
    k_f = prefactor_units[115] * fwd_A[115]
                * exp(fwd_beta[115] * tc[0] - activation_units[115] * fwd_Ea[115] * invT);
    dlnkfdT = fwd_beta[115] * invT + activation_units[115] * fwd_Ea[115] * invT2;
    /* reverse */
    phi_r = sc[8]*sc[15];
    Kc = exp(g_RT[6] - g_RT[8] - g_RT[15] + g_RT[23]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[6] + h_RT[23]) + (h_RT[8] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[6] -= q; /* H2O2 */
    wdot[8] += q; /* HO2 */
    wdot[15] += q; /* C2H4 */
    wdot[23] -= q; /* C2H3 */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[23];
    J[186] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    J[188] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[195] += dqdci;              /* dwdot[C2H4]/d[H2O2] */
    J[203] -= dqdci;              /* dwdot[C2H3]/d[H2O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[15];
    J[246] -= dqdci;              /* dwdot[H2O2]/d[HO2] */
    J[248] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[255] += dqdci;              /* dwdot[C2H4]/d[HO2] */
    J[263] -= dqdci;              /* dwdot[C2H3]/d[HO2] */
    /* d()/d[C2H4] */
    dqdci =  - k_r*sc[8];
    J[456] -= dqdci;              /* dwdot[H2O2]/d[C2H4] */
    J[458] += dqdci;              /* dwdot[HO2]/d[C2H4] */
    J[465] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[473] -= dqdci;              /* dwdot[C2H3]/d[C2H4] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[6];
    J[696] -= dqdci;              /* dwdot[H2O2]/d[C2H3] */
    J[698] += dqdci;              /* dwdot[HO2]/d[C2H3] */
    J[705] += dqdci;              /* dwdot[C2H4]/d[C2H3] */
    J[713] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    /* d()/dT */
    J[876] -= dqdT;               /* dwdot[H2O2]/dT */
    J[878] += dqdT;               /* dwdot[HO2]/dT */
    J[885] += dqdT;               /* dwdot[C2H4]/dT */
    J[893] -= dqdT;               /* dwdot[C2H3]/dT */

    /*reaction 117: C2H3 + O2 <=> C2H2 + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[23];
    k_f = prefactor_units[116] * fwd_A[116]
                * exp(fwd_beta[116] * tc[0] - activation_units[116] * fwd_Ea[116] * invT);
    dlnkfdT = fwd_beta[116] * invT + activation_units[116] * fwd_Ea[116] * invT2;
    /* reverse */
    phi_r = sc[8]*sc[14];
    Kc = exp(g_RT[5] - g_RT[8] - g_RT[14] + g_RT[23]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[23]) + (h_RT[8] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[8] += q; /* HO2 */
    wdot[14] += q; /* C2H2 */
    wdot[23] -= q; /* C2H3 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[23];
    J[155] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[158] += dqdci;              /* dwdot[HO2]/d[O2] */
    J[164] += dqdci;              /* dwdot[C2H2]/d[O2] */
    J[173] -= dqdci;              /* dwdot[C2H3]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[14];
    J[245] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[248] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[254] += dqdci;              /* dwdot[C2H2]/d[HO2] */
    J[263] -= dqdci;              /* dwdot[C2H3]/d[HO2] */
    /* d()/d[C2H2] */
    dqdci =  - k_r*sc[8];
    J[425] -= dqdci;              /* dwdot[O2]/d[C2H2] */
    J[428] += dqdci;              /* dwdot[HO2]/d[C2H2] */
    J[434] += dqdci;              /* dwdot[C2H2]/d[C2H2] */
    J[443] -= dqdci;              /* dwdot[C2H3]/d[C2H2] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[5];
    J[695] -= dqdci;              /* dwdot[O2]/d[C2H3] */
    J[698] += dqdci;              /* dwdot[HO2]/d[C2H3] */
    J[704] += dqdci;              /* dwdot[C2H2]/d[C2H3] */
    J[713] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    /* d()/dT */
    J[875] -= dqdT;               /* dwdot[O2]/dT */
    J[878] += dqdT;               /* dwdot[HO2]/dT */
    J[884] += dqdT;               /* dwdot[C2H2]/dT */
    J[893] -= dqdT;               /* dwdot[C2H3]/dT */

    /*reaction 118: C2H4 + CH3 <=> C2H3 + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[10]*sc[15];
    k_f = prefactor_units[117] * fwd_A[117]
                * exp(fwd_beta[117] * tc[0] - activation_units[117] * fwd_Ea[117] * invT);
    dlnkfdT = fwd_beta[117] * invT + activation_units[117] * fwd_Ea[117] * invT2;
    /* reverse */
    phi_r = sc[13]*sc[23];
    Kc = exp(g_RT[10] - g_RT[13] + g_RT[15] - g_RT[23]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[10] + h_RT[15]) + (h_RT[13] + h_RT[23]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[10] -= q; /* CH3 */
    wdot[13] += q; /* CH4 */
    wdot[15] -= q; /* C2H4 */
    wdot[23] += q; /* C2H3 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[15];
    J[310] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[313] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[315] -= dqdci;              /* dwdot[C2H4]/d[CH3] */
    J[323] += dqdci;              /* dwdot[C2H3]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[23];
    J[400] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[403] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[405] -= dqdci;              /* dwdot[C2H4]/d[CH4] */
    J[413] += dqdci;              /* dwdot[C2H3]/d[CH4] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[10];
    J[460] -= dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[463] += dqdci;              /* dwdot[CH4]/d[C2H4] */
    J[465] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[473] += dqdci;              /* dwdot[C2H3]/d[C2H4] */
    /* d()/d[C2H3] */
    dqdci =  - k_r*sc[13];
    J[700] -= dqdci;              /* dwdot[CH3]/d[C2H3] */
    J[703] += dqdci;              /* dwdot[CH4]/d[C2H3] */
    J[705] -= dqdci;              /* dwdot[C2H4]/d[C2H3] */
    J[713] += dqdci;              /* dwdot[C2H3]/d[C2H3] */
    /* d()/dT */
    J[880] -= dqdT;               /* dwdot[CH3]/dT */
    J[883] += dqdT;               /* dwdot[CH4]/dT */
    J[885] -= dqdT;               /* dwdot[C2H4]/dT */
    J[893] += dqdT;               /* dwdot[C2H3]/dT */

    /*reaction 119: C2H4 + O2 => CH3 + CO2 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[15];
    k_f = prefactor_units[118] * fwd_A[118]
                * exp(fwd_beta[118] * tc[0] - activation_units[118] * fwd_Ea[118] * invT);
    dlnkfdT = fwd_beta[118] * invT + activation_units[118] * fwd_Ea[118] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[5] -= q; /* O2 */
    wdot[10] += q; /* CH3 */
    wdot[12] += q; /* CO2 */
    wdot[15] -= q; /* C2H4 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[15];
    J[151] += dqdci;              /* dwdot[H]/d[O2] */
    J[155] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[160] += dqdci;              /* dwdot[CH3]/d[O2] */
    J[162] += dqdci;              /* dwdot[CO2]/d[O2] */
    J[165] -= dqdci;              /* dwdot[C2H4]/d[O2] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[5];
    J[451] += dqdci;              /* dwdot[H]/d[C2H4] */
    J[455] -= dqdci;              /* dwdot[O2]/d[C2H4] */
    J[460] += dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[462] += dqdci;              /* dwdot[CO2]/d[C2H4] */
    J[465] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    /* d()/dT */
    J[871] += dqdT;               /* dwdot[H]/dT */
    J[875] -= dqdT;               /* dwdot[O2]/dT */
    J[880] += dqdT;               /* dwdot[CH3]/dT */
    J[882] += dqdT;               /* dwdot[CO2]/dT */
    J[885] -= dqdT;               /* dwdot[C2H4]/dT */

    /*reaction 120: C2H4 + OH <=> C2H3 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[15];
    k_f = prefactor_units[119] * fwd_A[119]
                * exp(fwd_beta[119] * tc[0] - activation_units[119] * fwd_Ea[119] * invT);
    dlnkfdT = fwd_beta[119] * invT + activation_units[119] * fwd_Ea[119] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[23];
    Kc = exp(g_RT[4] - g_RT[7] + g_RT[15] - g_RT[23]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[15]) + (h_RT[7] + h_RT[23]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[7] += q; /* H2O */
    wdot[15] -= q; /* C2H4 */
    wdot[23] += q; /* C2H3 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[15];
    J[124] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[127] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[135] -= dqdci;              /* dwdot[C2H4]/d[OH] */
    J[143] += dqdci;              /* dwdot[C2H3]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[23];
    J[214] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[217] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[225] -= dqdci;              /* dwdot[C2H4]/d[H2O] */
    J[233] += dqdci;              /* dwdot[C2H3]/d[H2O] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[4];
    J[454] -= dqdci;              /* dwdot[OH]/d[C2H4] */
    J[457] += dqdci;              /* dwdot[H2O]/d[C2H4] */
    J[465] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[473] += dqdci;              /* dwdot[C2H3]/d[C2H4] */
    /* d()/d[C2H3] */
    dqdci =  - k_r*sc[7];
    J[694] -= dqdci;              /* dwdot[OH]/d[C2H3] */
    J[697] += dqdci;              /* dwdot[H2O]/d[C2H3] */
    J[705] -= dqdci;              /* dwdot[C2H4]/d[C2H3] */
    J[713] += dqdci;              /* dwdot[C2H3]/d[C2H3] */
    /* d()/dT */
    J[874] -= dqdT;               /* dwdot[OH]/dT */
    J[877] += dqdT;               /* dwdot[H2O]/dT */
    J[885] -= dqdT;               /* dwdot[C2H4]/dT */
    J[893] += dqdT;               /* dwdot[C2H3]/dT */

    /*reaction 121: C2H4 + OH <=> C2H5O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[15];
    k_f = prefactor_units[120] * fwd_A[120]
                * exp(fwd_beta[120] * tc[0] - activation_units[120] * fwd_Ea[120] * invT);
    dlnkfdT = fwd_beta[120] * invT + activation_units[120] * fwd_Ea[120] * invT2;
    /* reverse */
    phi_r = sc[28];
    Kc = refCinv * exp(g_RT[4] + g_RT[15] - g_RT[28]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[15]) + (h_RT[28]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[15] -= q; /* C2H4 */
    wdot[28] += q; /* C2H5O */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[15];
    J[124] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[135] -= dqdci;              /* dwdot[C2H4]/d[OH] */
    J[148] += dqdci;              /* dwdot[C2H5O]/d[OH] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[4];
    J[454] -= dqdci;              /* dwdot[OH]/d[C2H4] */
    J[465] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[478] += dqdci;              /* dwdot[C2H5O]/d[C2H4] */
    /* d()/d[C2H5O] */
    dqdci =  - k_r;
    J[844] -= dqdci;              /* dwdot[OH]/d[C2H5O] */
    J[855] -= dqdci;              /* dwdot[C2H4]/d[C2H5O] */
    J[868] += dqdci;              /* dwdot[C2H5O]/d[C2H5O] */
    /* d()/dT */
    J[874] -= dqdT;               /* dwdot[OH]/dT */
    J[885] -= dqdT;               /* dwdot[C2H4]/dT */
    J[898] += dqdT;               /* dwdot[C2H5O]/dT */

    /*reaction 122: C2H4 + O <=> CH2CHO + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[15];
    k_f = prefactor_units[121] * fwd_A[121]
                * exp(fwd_beta[121] * tc[0] - activation_units[121] * fwd_Ea[121] * invT);
    dlnkfdT = fwd_beta[121] * invT + activation_units[121] * fwd_Ea[121] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[27];
    Kc = exp(-g_RT[1] + g_RT[3] + g_RT[15] - g_RT[27]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[15]) + (h_RT[1] + h_RT[27]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[3] -= q; /* O */
    wdot[15] -= q; /* C2H4 */
    wdot[27] += q; /* CH2CHO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[27];
    J[31] += dqdci;               /* dwdot[H]/d[H] */
    J[33] -= dqdci;               /* dwdot[O]/d[H] */
    J[45] -= dqdci;               /* dwdot[C2H4]/d[H] */
    J[57] += dqdci;               /* dwdot[CH2CHO]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[15];
    J[91] += dqdci;               /* dwdot[H]/d[O] */
    J[93] -= dqdci;               /* dwdot[O]/d[O] */
    J[105] -= dqdci;              /* dwdot[C2H4]/d[O] */
    J[117] += dqdci;              /* dwdot[CH2CHO]/d[O] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[3];
    J[451] += dqdci;              /* dwdot[H]/d[C2H4] */
    J[453] -= dqdci;              /* dwdot[O]/d[C2H4] */
    J[465] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[477] += dqdci;              /* dwdot[CH2CHO]/d[C2H4] */
    /* d()/d[CH2CHO] */
    dqdci =  - k_r*sc[1];
    J[811] += dqdci;              /* dwdot[H]/d[CH2CHO] */
    J[813] -= dqdci;              /* dwdot[O]/d[CH2CHO] */
    J[825] -= dqdci;              /* dwdot[C2H4]/d[CH2CHO] */
    J[837] += dqdci;              /* dwdot[CH2CHO]/d[CH2CHO] */
    /* d()/dT */
    J[871] += dqdT;               /* dwdot[H]/dT */
    J[873] -= dqdT;               /* dwdot[O]/dT */
    J[885] -= dqdT;               /* dwdot[C2H4]/dT */
    J[897] += dqdT;               /* dwdot[CH2CHO]/dT */

    /*reaction 123: C2H4 + O <=> CH3 + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[15];
    k_f = prefactor_units[122] * fwd_A[122]
                * exp(fwd_beta[122] * tc[0] - activation_units[122] * fwd_Ea[122] * invT);
    dlnkfdT = fwd_beta[122] * invT + activation_units[122] * fwd_Ea[122] * invT2;
    /* reverse */
    phi_r = sc[10]*sc[20];
    Kc = exp(g_RT[3] - g_RT[10] + g_RT[15] - g_RT[20]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[15]) + (h_RT[10] + h_RT[20]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O */
    wdot[10] += q; /* CH3 */
    wdot[15] -= q; /* C2H4 */
    wdot[20] += q; /* HCO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[15];
    J[93] -= dqdci;               /* dwdot[O]/d[O] */
    J[100] += dqdci;              /* dwdot[CH3]/d[O] */
    J[105] -= dqdci;              /* dwdot[C2H4]/d[O] */
    J[110] += dqdci;              /* dwdot[HCO]/d[O] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[20];
    J[303] -= dqdci;              /* dwdot[O]/d[CH3] */
    J[310] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[315] -= dqdci;              /* dwdot[C2H4]/d[CH3] */
    J[320] += dqdci;              /* dwdot[HCO]/d[CH3] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[3];
    J[453] -= dqdci;              /* dwdot[O]/d[C2H4] */
    J[460] += dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[465] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[470] += dqdci;              /* dwdot[HCO]/d[C2H4] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[10];
    J[603] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[610] += dqdci;              /* dwdot[CH3]/d[HCO] */
    J[615] -= dqdci;              /* dwdot[C2H4]/d[HCO] */
    J[620] += dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[873] -= dqdT;               /* dwdot[O]/dT */
    J[880] += dqdT;               /* dwdot[CH3]/dT */
    J[885] -= dqdT;               /* dwdot[C2H4]/dT */
    J[890] += dqdT;               /* dwdot[HCO]/dT */

    /*reaction 124: C2H4 + O2 <=> C2H3 + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[15];
    k_f = prefactor_units[123] * fwd_A[123]
                * exp(fwd_beta[123] * tc[0] - activation_units[123] * fwd_Ea[123] * invT);
    dlnkfdT = fwd_beta[123] * invT + activation_units[123] * fwd_Ea[123] * invT2;
    /* reverse */
    phi_r = sc[8]*sc[23];
    Kc = exp(g_RT[5] - g_RT[8] + g_RT[15] - g_RT[23]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[15]) + (h_RT[8] + h_RT[23]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[8] += q; /* HO2 */
    wdot[15] -= q; /* C2H4 */
    wdot[23] += q; /* C2H3 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[15];
    J[155] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[158] += dqdci;              /* dwdot[HO2]/d[O2] */
    J[165] -= dqdci;              /* dwdot[C2H4]/d[O2] */
    J[173] += dqdci;              /* dwdot[C2H3]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[23];
    J[245] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[248] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[255] -= dqdci;              /* dwdot[C2H4]/d[HO2] */
    J[263] += dqdci;              /* dwdot[C2H3]/d[HO2] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[5];
    J[455] -= dqdci;              /* dwdot[O2]/d[C2H4] */
    J[458] += dqdci;              /* dwdot[HO2]/d[C2H4] */
    J[465] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[473] += dqdci;              /* dwdot[C2H3]/d[C2H4] */
    /* d()/d[C2H3] */
    dqdci =  - k_r*sc[8];
    J[695] -= dqdci;              /* dwdot[O2]/d[C2H3] */
    J[698] += dqdci;              /* dwdot[HO2]/d[C2H3] */
    J[705] -= dqdci;              /* dwdot[C2H4]/d[C2H3] */
    J[713] += dqdci;              /* dwdot[C2H3]/d[C2H3] */
    /* d()/dT */
    J[875] -= dqdT;               /* dwdot[O2]/dT */
    J[878] += dqdT;               /* dwdot[HO2]/dT */
    J[885] -= dqdT;               /* dwdot[C2H4]/dT */
    J[893] += dqdT;               /* dwdot[C2H3]/dT */

    /*reaction 125: C2H4 + H <=> C2H3 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[15];
    k_f = prefactor_units[124] * fwd_A[124]
                * exp(fwd_beta[124] * tc[0] - activation_units[124] * fwd_Ea[124] * invT);
    dlnkfdT = fwd_beta[124] * invT + activation_units[124] * fwd_Ea[124] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[23];
    Kc = exp(g_RT[1] - g_RT[2] + g_RT[15] - g_RT[23]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[15]) + (h_RT[2] + h_RT[23]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[2] += q; /* H2 */
    wdot[15] -= q; /* C2H4 */
    wdot[23] += q; /* C2H3 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[15];
    J[31] -= dqdci;               /* dwdot[H]/d[H] */
    J[32] += dqdci;               /* dwdot[H2]/d[H] */
    J[45] -= dqdci;               /* dwdot[C2H4]/d[H] */
    J[53] += dqdci;               /* dwdot[C2H3]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[23];
    J[61] -= dqdci;               /* dwdot[H]/d[H2] */
    J[62] += dqdci;               /* dwdot[H2]/d[H2] */
    J[75] -= dqdci;               /* dwdot[C2H4]/d[H2] */
    J[83] += dqdci;               /* dwdot[C2H3]/d[H2] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[1];
    J[451] -= dqdci;              /* dwdot[H]/d[C2H4] */
    J[452] += dqdci;              /* dwdot[H2]/d[C2H4] */
    J[465] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[473] += dqdci;              /* dwdot[C2H3]/d[C2H4] */
    /* d()/d[C2H3] */
    dqdci =  - k_r*sc[2];
    J[691] -= dqdci;              /* dwdot[H]/d[C2H3] */
    J[692] += dqdci;              /* dwdot[H2]/d[C2H3] */
    J[705] -= dqdci;              /* dwdot[C2H4]/d[C2H3] */
    J[713] += dqdci;              /* dwdot[C2H3]/d[C2H3] */
    /* d()/dT */
    J[871] -= dqdT;               /* dwdot[H]/dT */
    J[872] += dqdT;               /* dwdot[H2]/dT */
    J[885] -= dqdT;               /* dwdot[C2H4]/dT */
    J[893] += dqdT;               /* dwdot[C2H3]/dT */

    /*reaction 126: C2H4 + O <=> TXCH2 + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[15];
    k_f = prefactor_units[125] * fwd_A[125]
                * exp(fwd_beta[125] * tc[0] - activation_units[125] * fwd_Ea[125] * invT);
    dlnkfdT = fwd_beta[125] * invT + activation_units[125] * fwd_Ea[125] * invT2;
    /* reverse */
    phi_r = sc[11]*sc[21];
    Kc = exp(g_RT[3] - g_RT[11] + g_RT[15] - g_RT[21]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[15]) + (h_RT[11] + h_RT[21]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O */
    wdot[11] += q; /* CH2O */
    wdot[15] -= q; /* C2H4 */
    wdot[21] += q; /* TXCH2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[15];
    J[93] -= dqdci;               /* dwdot[O]/d[O] */
    J[101] += dqdci;              /* dwdot[CH2O]/d[O] */
    J[105] -= dqdci;              /* dwdot[C2H4]/d[O] */
    J[111] += dqdci;              /* dwdot[TXCH2]/d[O] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[21];
    J[333] -= dqdci;              /* dwdot[O]/d[CH2O] */
    J[341] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[345] -= dqdci;              /* dwdot[C2H4]/d[CH2O] */
    J[351] += dqdci;              /* dwdot[TXCH2]/d[CH2O] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[3];
    J[453] -= dqdci;              /* dwdot[O]/d[C2H4] */
    J[461] += dqdci;              /* dwdot[CH2O]/d[C2H4] */
    J[465] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[471] += dqdci;              /* dwdot[TXCH2]/d[C2H4] */
    /* d()/d[TXCH2] */
    dqdci =  - k_r*sc[11];
    J[633] -= dqdci;              /* dwdot[O]/d[TXCH2] */
    J[641] += dqdci;              /* dwdot[CH2O]/d[TXCH2] */
    J[645] -= dqdci;              /* dwdot[C2H4]/d[TXCH2] */
    J[651] += dqdci;              /* dwdot[TXCH2]/d[TXCH2] */
    /* d()/dT */
    J[873] -= dqdT;               /* dwdot[O]/dT */
    J[881] += dqdT;               /* dwdot[CH2O]/dT */
    J[885] -= dqdT;               /* dwdot[C2H4]/dT */
    J[891] += dqdT;               /* dwdot[TXCH2]/dT */

    /*reaction 127: C2H5 + HO2 <=> C2H4 + H2O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[24];
    k_f = prefactor_units[126] * fwd_A[126]
                * exp(fwd_beta[126] * tc[0] - activation_units[126] * fwd_Ea[126] * invT);
    dlnkfdT = fwd_beta[126] * invT + activation_units[126] * fwd_Ea[126] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[15];
    Kc = exp(-g_RT[6] + g_RT[8] - g_RT[15] + g_RT[24]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[24]) + (h_RT[6] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[6] += q; /* H2O2 */
    wdot[8] -= q; /* HO2 */
    wdot[15] += q; /* C2H4 */
    wdot[24] -= q; /* C2H5 */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[15];
    J[186] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
    J[188] -= dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[195] += dqdci;              /* dwdot[C2H4]/d[H2O2] */
    J[204] -= dqdci;              /* dwdot[C2H5]/d[H2O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[24];
    J[246] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    J[248] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[255] += dqdci;              /* dwdot[C2H4]/d[HO2] */
    J[264] -= dqdci;              /* dwdot[C2H5]/d[HO2] */
    /* d()/d[C2H4] */
    dqdci =  - k_r*sc[6];
    J[456] += dqdci;              /* dwdot[H2O2]/d[C2H4] */
    J[458] -= dqdci;              /* dwdot[HO2]/d[C2H4] */
    J[465] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[474] -= dqdci;              /* dwdot[C2H5]/d[C2H4] */
    /* d()/d[C2H5] */
    dqdci =  + k_f*sc[8];
    J[726] += dqdci;              /* dwdot[H2O2]/d[C2H5] */
    J[728] -= dqdci;              /* dwdot[HO2]/d[C2H5] */
    J[735] += dqdci;              /* dwdot[C2H4]/d[C2H5] */
    J[744] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
    /* d()/dT */
    J[876] += dqdT;               /* dwdot[H2O2]/dT */
    J[878] -= dqdT;               /* dwdot[HO2]/dT */
    J[885] += dqdT;               /* dwdot[C2H4]/dT */
    J[894] -= dqdT;               /* dwdot[C2H5]/dT */

    /*reaction 128: C2H5 + HO2 <=> C2H5O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[24];
    k_f = prefactor_units[127] * fwd_A[127]
                * exp(fwd_beta[127] * tc[0] - activation_units[127] * fwd_Ea[127] * invT);
    dlnkfdT = fwd_beta[127] * invT + activation_units[127] * fwd_Ea[127] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[28];
    Kc = exp(-g_RT[4] + g_RT[8] + g_RT[24] - g_RT[28]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[24]) + (h_RT[4] + h_RT[28]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] += q; /* OH */
    wdot[8] -= q; /* HO2 */
    wdot[24] -= q; /* C2H5 */
    wdot[28] += q; /* C2H5O */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[28];
    J[124] += dqdci;              /* dwdot[OH]/d[OH] */
    J[128] -= dqdci;              /* dwdot[HO2]/d[OH] */
    J[144] -= dqdci;              /* dwdot[C2H5]/d[OH] */
    J[148] += dqdci;              /* dwdot[C2H5O]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[24];
    J[244] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[248] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[264] -= dqdci;              /* dwdot[C2H5]/d[HO2] */
    J[268] += dqdci;              /* dwdot[C2H5O]/d[HO2] */
    /* d()/d[C2H5] */
    dqdci =  + k_f*sc[8];
    J[724] += dqdci;              /* dwdot[OH]/d[C2H5] */
    J[728] -= dqdci;              /* dwdot[HO2]/d[C2H5] */
    J[744] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[748] += dqdci;              /* dwdot[C2H5O]/d[C2H5] */
    /* d()/d[C2H5O] */
    dqdci =  - k_r*sc[4];
    J[844] += dqdci;              /* dwdot[OH]/d[C2H5O] */
    J[848] -= dqdci;              /* dwdot[HO2]/d[C2H5O] */
    J[864] -= dqdci;              /* dwdot[C2H5]/d[C2H5O] */
    J[868] += dqdci;              /* dwdot[C2H5O]/d[C2H5O] */
    /* d()/dT */
    J[874] += dqdT;               /* dwdot[OH]/dT */
    J[878] -= dqdT;               /* dwdot[HO2]/dT */
    J[894] -= dqdT;               /* dwdot[C2H5]/dT */
    J[898] += dqdT;               /* dwdot[C2H5O]/dT */

    /*reaction 129: C2H5 + O <=> C2H5O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[24];
    k_f = prefactor_units[128] * fwd_A[128]
                * exp(fwd_beta[128] * tc[0] - activation_units[128] * fwd_Ea[128] * invT);
    dlnkfdT = fwd_beta[128] * invT + activation_units[128] * fwd_Ea[128] * invT2;
    /* reverse */
    phi_r = sc[28];
    Kc = refCinv * exp(g_RT[3] + g_RT[24] - g_RT[28]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[24]) + (h_RT[28]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O */
    wdot[24] -= q; /* C2H5 */
    wdot[28] += q; /* C2H5O */
    /* d()/d[O] */
    dqdci =  + k_f*sc[24];
    J[93] -= dqdci;               /* dwdot[O]/d[O] */
    J[114] -= dqdci;              /* dwdot[C2H5]/d[O] */
    J[118] += dqdci;              /* dwdot[C2H5O]/d[O] */
    /* d()/d[C2H5] */
    dqdci =  + k_f*sc[3];
    J[723] -= dqdci;              /* dwdot[O]/d[C2H5] */
    J[744] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[748] += dqdci;              /* dwdot[C2H5O]/d[C2H5] */
    /* d()/d[C2H5O] */
    dqdci =  - k_r;
    J[843] -= dqdci;              /* dwdot[O]/d[C2H5O] */
    J[864] -= dqdci;              /* dwdot[C2H5]/d[C2H5O] */
    J[868] += dqdci;              /* dwdot[C2H5O]/d[C2H5O] */
    /* d()/dT */
    J[873] -= dqdT;               /* dwdot[O]/dT */
    J[894] -= dqdT;               /* dwdot[C2H5]/dT */
    J[898] += dqdT;               /* dwdot[C2H5O]/dT */

    /*reaction 130: C2H5 + H <=> C2H4 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[24];
    k_f = prefactor_units[129] * fwd_A[129]
                * exp(fwd_beta[129] * tc[0] - activation_units[129] * fwd_Ea[129] * invT);
    dlnkfdT = fwd_beta[129] * invT + activation_units[129] * fwd_Ea[129] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[15];
    Kc = exp(g_RT[1] - g_RT[2] - g_RT[15] + g_RT[24]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[24]) + (h_RT[2] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[2] += q; /* H2 */
    wdot[15] += q; /* C2H4 */
    wdot[24] -= q; /* C2H5 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[24];
    J[31] -= dqdci;               /* dwdot[H]/d[H] */
    J[32] += dqdci;               /* dwdot[H2]/d[H] */
    J[45] += dqdci;               /* dwdot[C2H4]/d[H] */
    J[54] -= dqdci;               /* dwdot[C2H5]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[15];
    J[61] -= dqdci;               /* dwdot[H]/d[H2] */
    J[62] += dqdci;               /* dwdot[H2]/d[H2] */
    J[75] += dqdci;               /* dwdot[C2H4]/d[H2] */
    J[84] -= dqdci;               /* dwdot[C2H5]/d[H2] */
    /* d()/d[C2H4] */
    dqdci =  - k_r*sc[2];
    J[451] -= dqdci;              /* dwdot[H]/d[C2H4] */
    J[452] += dqdci;              /* dwdot[H2]/d[C2H4] */
    J[465] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[474] -= dqdci;              /* dwdot[C2H5]/d[C2H4] */
    /* d()/d[C2H5] */
    dqdci =  + k_f*sc[1];
    J[721] -= dqdci;              /* dwdot[H]/d[C2H5] */
    J[722] += dqdci;              /* dwdot[H2]/d[C2H5] */
    J[735] += dqdci;              /* dwdot[C2H4]/d[C2H5] */
    J[744] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
    /* d()/dT */
    J[871] -= dqdT;               /* dwdot[H]/dT */
    J[872] += dqdT;               /* dwdot[H2]/dT */
    J[885] += dqdT;               /* dwdot[C2H4]/dT */
    J[894] -= dqdT;               /* dwdot[C2H5]/dT */

    /*reaction 131: C2H5 + O2 <=> C2H4 + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[24];
    k_f = prefactor_units[130] * fwd_A[130]
                * exp(fwd_beta[130] * tc[0] - activation_units[130] * fwd_Ea[130] * invT);
    dlnkfdT = fwd_beta[130] * invT + activation_units[130] * fwd_Ea[130] * invT2;
    /* reverse */
    phi_r = sc[8]*sc[15];
    Kc = exp(g_RT[5] - g_RT[8] - g_RT[15] + g_RT[24]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[24]) + (h_RT[8] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[8] += q; /* HO2 */
    wdot[15] += q; /* C2H4 */
    wdot[24] -= q; /* C2H5 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[24];
    J[155] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[158] += dqdci;              /* dwdot[HO2]/d[O2] */
    J[165] += dqdci;              /* dwdot[C2H4]/d[O2] */
    J[174] -= dqdci;              /* dwdot[C2H5]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[15];
    J[245] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[248] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[255] += dqdci;              /* dwdot[C2H4]/d[HO2] */
    J[264] -= dqdci;              /* dwdot[C2H5]/d[HO2] */
    /* d()/d[C2H4] */
    dqdci =  - k_r*sc[8];
    J[455] -= dqdci;              /* dwdot[O2]/d[C2H4] */
    J[458] += dqdci;              /* dwdot[HO2]/d[C2H4] */
    J[465] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[474] -= dqdci;              /* dwdot[C2H5]/d[C2H4] */
    /* d()/d[C2H5] */
    dqdci =  + k_f*sc[5];
    J[725] -= dqdci;              /* dwdot[O2]/d[C2H5] */
    J[728] += dqdci;              /* dwdot[HO2]/d[C2H5] */
    J[735] += dqdci;              /* dwdot[C2H4]/d[C2H5] */
    J[744] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
    /* d()/dT */
    J[875] -= dqdT;               /* dwdot[O2]/dT */
    J[878] += dqdT;               /* dwdot[HO2]/dT */
    J[885] += dqdT;               /* dwdot[C2H4]/dT */
    J[894] -= dqdT;               /* dwdot[C2H5]/dT */

    /*reaction 132: C2H5 + HO2 <=> C2H6 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[24];
    k_f = prefactor_units[131] * fwd_A[131]
                * exp(fwd_beta[131] * tc[0] - activation_units[131] * fwd_Ea[131] * invT);
    dlnkfdT = fwd_beta[131] * invT + activation_units[131] * fwd_Ea[131] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[17];
    Kc = exp(-g_RT[5] + g_RT[8] - g_RT[17] + g_RT[24]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[24]) + (h_RT[5] + h_RT[17]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[8] -= q; /* HO2 */
    wdot[17] += q; /* C2H6 */
    wdot[24] -= q; /* C2H5 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[17];
    J[155] += dqdci;              /* dwdot[O2]/d[O2] */
    J[158] -= dqdci;              /* dwdot[HO2]/d[O2] */
    J[167] += dqdci;              /* dwdot[C2H6]/d[O2] */
    J[174] -= dqdci;              /* dwdot[C2H5]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[24];
    J[245] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[248] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[257] += dqdci;              /* dwdot[C2H6]/d[HO2] */
    J[264] -= dqdci;              /* dwdot[C2H5]/d[HO2] */
    /* d()/d[C2H6] */
    dqdci =  - k_r*sc[5];
    J[515] += dqdci;              /* dwdot[O2]/d[C2H6] */
    J[518] -= dqdci;              /* dwdot[HO2]/d[C2H6] */
    J[527] += dqdci;              /* dwdot[C2H6]/d[C2H6] */
    J[534] -= dqdci;              /* dwdot[C2H5]/d[C2H6] */
    /* d()/d[C2H5] */
    dqdci =  + k_f*sc[8];
    J[725] += dqdci;              /* dwdot[O2]/d[C2H5] */
    J[728] -= dqdci;              /* dwdot[HO2]/d[C2H5] */
    J[737] += dqdci;              /* dwdot[C2H6]/d[C2H5] */
    J[744] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
    /* d()/dT */
    J[875] += dqdT;               /* dwdot[O2]/dT */
    J[878] -= dqdT;               /* dwdot[HO2]/dT */
    J[887] += dqdT;               /* dwdot[C2H6]/dT */
    J[894] -= dqdT;               /* dwdot[C2H5]/dT */

    /*reaction 133: C2H5 + CH3 <=> C2H4 + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[10]*sc[24];
    k_f = prefactor_units[132] * fwd_A[132]
                * exp(fwd_beta[132] * tc[0] - activation_units[132] * fwd_Ea[132] * invT);
    dlnkfdT = fwd_beta[132] * invT + activation_units[132] * fwd_Ea[132] * invT2;
    /* reverse */
    phi_r = sc[13]*sc[15];
    Kc = exp(g_RT[10] - g_RT[13] - g_RT[15] + g_RT[24]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[10] + h_RT[24]) + (h_RT[13] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[10] -= q; /* CH3 */
    wdot[13] += q; /* CH4 */
    wdot[15] += q; /* C2H4 */
    wdot[24] -= q; /* C2H5 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[24];
    J[310] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[313] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[315] += dqdci;              /* dwdot[C2H4]/d[CH3] */
    J[324] -= dqdci;              /* dwdot[C2H5]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[15];
    J[400] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[403] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[405] += dqdci;              /* dwdot[C2H4]/d[CH4] */
    J[414] -= dqdci;              /* dwdot[C2H5]/d[CH4] */
    /* d()/d[C2H4] */
    dqdci =  - k_r*sc[13];
    J[460] -= dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[463] += dqdci;              /* dwdot[CH4]/d[C2H4] */
    J[465] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[474] -= dqdci;              /* dwdot[C2H5]/d[C2H4] */
    /* d()/d[C2H5] */
    dqdci =  + k_f*sc[10];
    J[730] -= dqdci;              /* dwdot[CH3]/d[C2H5] */
    J[733] += dqdci;              /* dwdot[CH4]/d[C2H5] */
    J[735] += dqdci;              /* dwdot[C2H4]/d[C2H5] */
    J[744] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
    /* d()/dT */
    J[880] -= dqdT;               /* dwdot[CH3]/dT */
    J[883] += dqdT;               /* dwdot[CH4]/dT */
    J[885] += dqdT;               /* dwdot[C2H4]/dT */
    J[894] -= dqdT;               /* dwdot[C2H5]/dT */

    /*reaction 134: C2H6 + SXCH2 <=> C2H5 + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[17]*sc[22];
    k_f = prefactor_units[133] * fwd_A[133]
                * exp(fwd_beta[133] * tc[0] - activation_units[133] * fwd_Ea[133] * invT);
    dlnkfdT = fwd_beta[133] * invT + activation_units[133] * fwd_Ea[133] * invT2;
    /* reverse */
    phi_r = sc[10]*sc[24];
    Kc = exp(-g_RT[10] + g_RT[17] + g_RT[22] - g_RT[24]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[17] + h_RT[22]) + (h_RT[10] + h_RT[24]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[10] += q; /* CH3 */
    wdot[17] -= q; /* C2H6 */
    wdot[22] -= q; /* SXCH2 */
    wdot[24] += q; /* C2H5 */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[24];
    J[310] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[317] -= dqdci;              /* dwdot[C2H6]/d[CH3] */
    J[322] -= dqdci;              /* dwdot[SXCH2]/d[CH3] */
    J[324] += dqdci;              /* dwdot[C2H5]/d[CH3] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[22];
    J[520] += dqdci;              /* dwdot[CH3]/d[C2H6] */
    J[527] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    J[532] -= dqdci;              /* dwdot[SXCH2]/d[C2H6] */
    J[534] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    /* d()/d[SXCH2] */
    dqdci =  + k_f*sc[17];
    J[670] += dqdci;              /* dwdot[CH3]/d[SXCH2] */
    J[677] -= dqdci;              /* dwdot[C2H6]/d[SXCH2] */
    J[682] -= dqdci;              /* dwdot[SXCH2]/d[SXCH2] */
    J[684] += dqdci;              /* dwdot[C2H5]/d[SXCH2] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[10];
    J[730] += dqdci;              /* dwdot[CH3]/d[C2H5] */
    J[737] -= dqdci;              /* dwdot[C2H6]/d[C2H5] */
    J[742] -= dqdci;              /* dwdot[SXCH2]/d[C2H5] */
    J[744] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    /* d()/dT */
    J[880] += dqdT;               /* dwdot[CH3]/dT */
    J[887] -= dqdT;               /* dwdot[C2H6]/dT */
    J[892] -= dqdT;               /* dwdot[SXCH2]/dT */
    J[894] += dqdT;               /* dwdot[C2H5]/dT */

    /*reaction 135: C2H6 + CH3 <=> C2H5 + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[10]*sc[17];
    k_f = prefactor_units[134] * fwd_A[134]
                * exp(fwd_beta[134] * tc[0] - activation_units[134] * fwd_Ea[134] * invT);
    dlnkfdT = fwd_beta[134] * invT + activation_units[134] * fwd_Ea[134] * invT2;
    /* reverse */
    phi_r = sc[13]*sc[24];
    Kc = exp(g_RT[10] - g_RT[13] + g_RT[17] - g_RT[24]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[10] + h_RT[17]) + (h_RT[13] + h_RT[24]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[10] -= q; /* CH3 */
    wdot[13] += q; /* CH4 */
    wdot[17] -= q; /* C2H6 */
    wdot[24] += q; /* C2H5 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[17];
    J[310] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[313] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[317] -= dqdci;              /* dwdot[C2H6]/d[CH3] */
    J[324] += dqdci;              /* dwdot[C2H5]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[24];
    J[400] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[403] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[407] -= dqdci;              /* dwdot[C2H6]/d[CH4] */
    J[414] += dqdci;              /* dwdot[C2H5]/d[CH4] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[10];
    J[520] -= dqdci;              /* dwdot[CH3]/d[C2H6] */
    J[523] += dqdci;              /* dwdot[CH4]/d[C2H6] */
    J[527] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    J[534] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[13];
    J[730] -= dqdci;              /* dwdot[CH3]/d[C2H5] */
    J[733] += dqdci;              /* dwdot[CH4]/d[C2H5] */
    J[737] -= dqdci;              /* dwdot[C2H6]/d[C2H5] */
    J[744] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    /* d()/dT */
    J[880] -= dqdT;               /* dwdot[CH3]/dT */
    J[883] += dqdT;               /* dwdot[CH4]/dT */
    J[887] -= dqdT;               /* dwdot[C2H6]/dT */
    J[894] += dqdT;               /* dwdot[C2H5]/dT */

    /*reaction 136: C2H6 + O <=> C2H5 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[17];
    k_f = prefactor_units[135] * fwd_A[135]
                * exp(fwd_beta[135] * tc[0] - activation_units[135] * fwd_Ea[135] * invT);
    dlnkfdT = fwd_beta[135] * invT + activation_units[135] * fwd_Ea[135] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[24];
    Kc = exp(g_RT[3] - g_RT[4] + g_RT[17] - g_RT[24]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[17]) + (h_RT[4] + h_RT[24]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O */
    wdot[4] += q; /* OH */
    wdot[17] -= q; /* C2H6 */
    wdot[24] += q; /* C2H5 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[17];
    J[93] -= dqdci;               /* dwdot[O]/d[O] */
    J[94] += dqdci;               /* dwdot[OH]/d[O] */
    J[107] -= dqdci;              /* dwdot[C2H6]/d[O] */
    J[114] += dqdci;              /* dwdot[C2H5]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[24];
    J[123] -= dqdci;              /* dwdot[O]/d[OH] */
    J[124] += dqdci;              /* dwdot[OH]/d[OH] */
    J[137] -= dqdci;              /* dwdot[C2H6]/d[OH] */
    J[144] += dqdci;              /* dwdot[C2H5]/d[OH] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[3];
    J[513] -= dqdci;              /* dwdot[O]/d[C2H6] */
    J[514] += dqdci;              /* dwdot[OH]/d[C2H6] */
    J[527] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    J[534] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[4];
    J[723] -= dqdci;              /* dwdot[O]/d[C2H5] */
    J[724] += dqdci;              /* dwdot[OH]/d[C2H5] */
    J[737] -= dqdci;              /* dwdot[C2H6]/d[C2H5] */
    J[744] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    /* d()/dT */
    J[873] -= dqdT;               /* dwdot[O]/dT */
    J[874] += dqdT;               /* dwdot[OH]/dT */
    J[887] -= dqdT;               /* dwdot[C2H6]/dT */
    J[894] += dqdT;               /* dwdot[C2H5]/dT */

    /*reaction 137: C2H6 + HO2 <=> C2H5 + H2O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[17];
    k_f = prefactor_units[136] * fwd_A[136]
                * exp(fwd_beta[136] * tc[0] - activation_units[136] * fwd_Ea[136] * invT);
    dlnkfdT = fwd_beta[136] * invT + activation_units[136] * fwd_Ea[136] * invT2;
    /* reverse */
    phi_r = sc[6]*sc[24];
    Kc = exp(-g_RT[6] + g_RT[8] + g_RT[17] - g_RT[24]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[8] + h_RT[17]) + (h_RT[6] + h_RT[24]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[6] += q; /* H2O2 */
    wdot[8] -= q; /* HO2 */
    wdot[17] -= q; /* C2H6 */
    wdot[24] += q; /* C2H5 */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[24];
    J[186] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
    J[188] -= dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[197] -= dqdci;              /* dwdot[C2H6]/d[H2O2] */
    J[204] += dqdci;              /* dwdot[C2H5]/d[H2O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[17];
    J[246] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    J[248] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[257] -= dqdci;              /* dwdot[C2H6]/d[HO2] */
    J[264] += dqdci;              /* dwdot[C2H5]/d[HO2] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[8];
    J[516] += dqdci;              /* dwdot[H2O2]/d[C2H6] */
    J[518] -= dqdci;              /* dwdot[HO2]/d[C2H6] */
    J[527] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    J[534] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[6];
    J[726] += dqdci;              /* dwdot[H2O2]/d[C2H5] */
    J[728] -= dqdci;              /* dwdot[HO2]/d[C2H5] */
    J[737] -= dqdci;              /* dwdot[C2H6]/d[C2H5] */
    J[744] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    /* d()/dT */
    J[876] += dqdT;               /* dwdot[H2O2]/dT */
    J[878] -= dqdT;               /* dwdot[HO2]/dT */
    J[887] -= dqdT;               /* dwdot[C2H6]/dT */
    J[894] += dqdT;               /* dwdot[C2H5]/dT */

    /*reaction 138: C2H6 + H <=> C2H5 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[17];
    k_f = prefactor_units[137] * fwd_A[137]
                * exp(fwd_beta[137] * tc[0] - activation_units[137] * fwd_Ea[137] * invT);
    dlnkfdT = fwd_beta[137] * invT + activation_units[137] * fwd_Ea[137] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[24];
    Kc = exp(g_RT[1] - g_RT[2] + g_RT[17] - g_RT[24]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[17]) + (h_RT[2] + h_RT[24]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[2] += q; /* H2 */
    wdot[17] -= q; /* C2H6 */
    wdot[24] += q; /* C2H5 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[17];
    J[31] -= dqdci;               /* dwdot[H]/d[H] */
    J[32] += dqdci;               /* dwdot[H2]/d[H] */
    J[47] -= dqdci;               /* dwdot[C2H6]/d[H] */
    J[54] += dqdci;               /* dwdot[C2H5]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[24];
    J[61] -= dqdci;               /* dwdot[H]/d[H2] */
    J[62] += dqdci;               /* dwdot[H2]/d[H2] */
    J[77] -= dqdci;               /* dwdot[C2H6]/d[H2] */
    J[84] += dqdci;               /* dwdot[C2H5]/d[H2] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[1];
    J[511] -= dqdci;              /* dwdot[H]/d[C2H6] */
    J[512] += dqdci;              /* dwdot[H2]/d[C2H6] */
    J[527] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    J[534] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[2];
    J[721] -= dqdci;              /* dwdot[H]/d[C2H5] */
    J[722] += dqdci;              /* dwdot[H2]/d[C2H5] */
    J[737] -= dqdci;              /* dwdot[C2H6]/d[C2H5] */
    J[744] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    /* d()/dT */
    J[871] -= dqdT;               /* dwdot[H]/dT */
    J[872] += dqdT;               /* dwdot[H2]/dT */
    J[887] -= dqdT;               /* dwdot[C2H6]/dT */
    J[894] += dqdT;               /* dwdot[C2H5]/dT */

    /*reaction 139: C2H6 + OH <=> C2H5 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[17];
    k_f = prefactor_units[138] * fwd_A[138]
                * exp(fwd_beta[138] * tc[0] - activation_units[138] * fwd_Ea[138] * invT);
    dlnkfdT = fwd_beta[138] * invT + activation_units[138] * fwd_Ea[138] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[24];
    Kc = exp(g_RT[4] - g_RT[7] + g_RT[17] - g_RT[24]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[17]) + (h_RT[7] + h_RT[24]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[7] += q; /* H2O */
    wdot[17] -= q; /* C2H6 */
    wdot[24] += q; /* C2H5 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[17];
    J[124] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[127] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[137] -= dqdci;              /* dwdot[C2H6]/d[OH] */
    J[144] += dqdci;              /* dwdot[C2H5]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[24];
    J[214] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[217] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[227] -= dqdci;              /* dwdot[C2H6]/d[H2O] */
    J[234] += dqdci;              /* dwdot[C2H5]/d[H2O] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[4];
    J[514] -= dqdci;              /* dwdot[OH]/d[C2H6] */
    J[517] += dqdci;              /* dwdot[H2O]/d[C2H6] */
    J[527] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    J[534] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[7];
    J[724] -= dqdci;              /* dwdot[OH]/d[C2H5] */
    J[727] += dqdci;              /* dwdot[H2O]/d[C2H5] */
    J[737] -= dqdci;              /* dwdot[C2H6]/d[C2H5] */
    J[744] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    /* d()/dT */
    J[874] -= dqdT;               /* dwdot[OH]/dT */
    J[877] += dqdT;               /* dwdot[H2O]/dT */
    J[887] -= dqdT;               /* dwdot[C2H6]/dT */
    J[894] += dqdT;               /* dwdot[C2H5]/dT */

    /*reaction 140: HCCO + O2 <=> OH + 2.000000 CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[25];
    k_f = prefactor_units[139] * fwd_A[139]
                * exp(fwd_beta[139] * tc[0] - activation_units[139] * fwd_Ea[139] * invT);
    dlnkfdT = fwd_beta[139] * invT + activation_units[139] * fwd_Ea[139] * invT2;
    /* reverse */
    phi_r = sc[4]*pow(sc[9], 2.000000);
    Kc = refC * exp(-g_RT[4] + g_RT[5] - 2.000000*g_RT[9] + g_RT[25]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[25]) + (h_RT[4] + 2.000000*h_RT[9]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] += q; /* OH */
    wdot[5] -= q; /* O2 */
    wdot[9] += 2 * q; /* CO */
    wdot[25] -= q; /* HCCO */
    /* d()/d[OH] */
    dqdci =  - k_r*pow(sc[9], 2.000000);
    J[124] += dqdci;              /* dwdot[OH]/d[OH] */
    J[125] -= dqdci;              /* dwdot[O2]/d[OH] */
    J[129] += 2 * dqdci;          /* dwdot[CO]/d[OH] */
    J[145] -= dqdci;              /* dwdot[HCCO]/d[OH] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[25];
    J[154] += dqdci;              /* dwdot[OH]/d[O2] */
    J[155] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[159] += 2 * dqdci;          /* dwdot[CO]/d[O2] */
    J[175] -= dqdci;              /* dwdot[HCCO]/d[O2] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[4]*2.000000*sc[9];
    J[274] += dqdci;              /* dwdot[OH]/d[CO] */
    J[275] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[279] += 2 * dqdci;          /* dwdot[CO]/d[CO] */
    J[295] -= dqdci;              /* dwdot[HCCO]/d[CO] */
    /* d()/d[HCCO] */
    dqdci =  + k_f*sc[5];
    J[754] += dqdci;              /* dwdot[OH]/d[HCCO] */
    J[755] -= dqdci;              /* dwdot[O2]/d[HCCO] */
    J[759] += 2 * dqdci;          /* dwdot[CO]/d[HCCO] */
    J[775] -= dqdci;              /* dwdot[HCCO]/d[HCCO] */
    /* d()/dT */
    J[874] += dqdT;               /* dwdot[OH]/dT */
    J[875] -= dqdT;               /* dwdot[O2]/dT */
    J[879] += 2 * dqdT;           /* dwdot[CO]/dT */
    J[895] -= dqdT;               /* dwdot[HCCO]/dT */

    /*reaction 141: HCCO + O <=> H + 2.000000 CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[25];
    k_f = prefactor_units[140] * fwd_A[140]
                * exp(fwd_beta[140] * tc[0] - activation_units[140] * fwd_Ea[140] * invT);
    dlnkfdT = fwd_beta[140] * invT + activation_units[140] * fwd_Ea[140] * invT2;
    /* reverse */
    phi_r = sc[1]*pow(sc[9], 2.000000);
    Kc = refC * exp(-g_RT[1] + g_RT[3] - 2.000000*g_RT[9] + g_RT[25]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[25]) + (h_RT[1] + 2.000000*h_RT[9]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[3] -= q; /* O */
    wdot[9] += 2 * q; /* CO */
    wdot[25] -= q; /* HCCO */
    /* d()/d[H] */
    dqdci =  - k_r*pow(sc[9], 2.000000);
    J[31] += dqdci;               /* dwdot[H]/d[H] */
    J[33] -= dqdci;               /* dwdot[O]/d[H] */
    J[39] += 2 * dqdci;           /* dwdot[CO]/d[H] */
    J[55] -= dqdci;               /* dwdot[HCCO]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[25];
    J[91] += dqdci;               /* dwdot[H]/d[O] */
    J[93] -= dqdci;               /* dwdot[O]/d[O] */
    J[99] += 2 * dqdci;           /* dwdot[CO]/d[O] */
    J[115] -= dqdci;              /* dwdot[HCCO]/d[O] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[1]*2.000000*sc[9];
    J[271] += dqdci;              /* dwdot[H]/d[CO] */
    J[273] -= dqdci;              /* dwdot[O]/d[CO] */
    J[279] += 2 * dqdci;          /* dwdot[CO]/d[CO] */
    J[295] -= dqdci;              /* dwdot[HCCO]/d[CO] */
    /* d()/d[HCCO] */
    dqdci =  + k_f*sc[3];
    J[751] += dqdci;              /* dwdot[H]/d[HCCO] */
    J[753] -= dqdci;              /* dwdot[O]/d[HCCO] */
    J[759] += 2 * dqdci;          /* dwdot[CO]/d[HCCO] */
    J[775] -= dqdci;              /* dwdot[HCCO]/d[HCCO] */
    /* d()/dT */
    J[871] += dqdT;               /* dwdot[H]/dT */
    J[873] -= dqdT;               /* dwdot[O]/dT */
    J[879] += 2 * dqdT;           /* dwdot[CO]/dT */
    J[895] -= dqdT;               /* dwdot[HCCO]/dT */

    /*reaction 142: HCCO + CH3 <=> C2H4 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[10]*sc[25];
    k_f = prefactor_units[141] * fwd_A[141]
                * exp(fwd_beta[141] * tc[0] - activation_units[141] * fwd_Ea[141] * invT);
    dlnkfdT = fwd_beta[141] * invT + activation_units[141] * fwd_Ea[141] * invT2;
    /* reverse */
    phi_r = sc[9]*sc[15];
    Kc = exp(-g_RT[9] + g_RT[10] - g_RT[15] + g_RT[25]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[10] + h_RT[25]) + (h_RT[9] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[9] += q; /* CO */
    wdot[10] -= q; /* CH3 */
    wdot[15] += q; /* C2H4 */
    wdot[25] -= q; /* HCCO */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[15];
    J[279] += dqdci;              /* dwdot[CO]/d[CO] */
    J[280] -= dqdci;              /* dwdot[CH3]/d[CO] */
    J[285] += dqdci;              /* dwdot[C2H4]/d[CO] */
    J[295] -= dqdci;              /* dwdot[HCCO]/d[CO] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[25];
    J[309] += dqdci;              /* dwdot[CO]/d[CH3] */
    J[310] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[315] += dqdci;              /* dwdot[C2H4]/d[CH3] */
    J[325] -= dqdci;              /* dwdot[HCCO]/d[CH3] */
    /* d()/d[C2H4] */
    dqdci =  - k_r*sc[9];
    J[459] += dqdci;              /* dwdot[CO]/d[C2H4] */
    J[460] -= dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[465] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[475] -= dqdci;              /* dwdot[HCCO]/d[C2H4] */
    /* d()/d[HCCO] */
    dqdci =  + k_f*sc[10];
    J[759] += dqdci;              /* dwdot[CO]/d[HCCO] */
    J[760] -= dqdci;              /* dwdot[CH3]/d[HCCO] */
    J[765] += dqdci;              /* dwdot[C2H4]/d[HCCO] */
    J[775] -= dqdci;              /* dwdot[HCCO]/d[HCCO] */
    /* d()/dT */
    J[879] += dqdT;               /* dwdot[CO]/dT */
    J[880] -= dqdT;               /* dwdot[CH3]/dT */
    J[885] += dqdT;               /* dwdot[C2H4]/dT */
    J[895] -= dqdT;               /* dwdot[HCCO]/dT */

    /*reaction 143: HCCO + H <=> SXCH2 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[25];
    k_f = prefactor_units[142] * fwd_A[142]
                * exp(fwd_beta[142] * tc[0] - activation_units[142] * fwd_Ea[142] * invT);
    dlnkfdT = fwd_beta[142] * invT + activation_units[142] * fwd_Ea[142] * invT2;
    /* reverse */
    phi_r = sc[9]*sc[22];
    Kc = exp(g_RT[1] - g_RT[9] - g_RT[22] + g_RT[25]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[25]) + (h_RT[9] + h_RT[22]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[9] += q; /* CO */
    wdot[22] += q; /* SXCH2 */
    wdot[25] -= q; /* HCCO */
    /* d()/d[H] */
    dqdci =  + k_f*sc[25];
    J[31] -= dqdci;               /* dwdot[H]/d[H] */
    J[39] += dqdci;               /* dwdot[CO]/d[H] */
    J[52] += dqdci;               /* dwdot[SXCH2]/d[H] */
    J[55] -= dqdci;               /* dwdot[HCCO]/d[H] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[22];
    J[271] -= dqdci;              /* dwdot[H]/d[CO] */
    J[279] += dqdci;              /* dwdot[CO]/d[CO] */
    J[292] += dqdci;              /* dwdot[SXCH2]/d[CO] */
    J[295] -= dqdci;              /* dwdot[HCCO]/d[CO] */
    /* d()/d[SXCH2] */
    dqdci =  - k_r*sc[9];
    J[661] -= dqdci;              /* dwdot[H]/d[SXCH2] */
    J[669] += dqdci;              /* dwdot[CO]/d[SXCH2] */
    J[682] += dqdci;              /* dwdot[SXCH2]/d[SXCH2] */
    J[685] -= dqdci;              /* dwdot[HCCO]/d[SXCH2] */
    /* d()/d[HCCO] */
    dqdci =  + k_f*sc[1];
    J[751] -= dqdci;              /* dwdot[H]/d[HCCO] */
    J[759] += dqdci;              /* dwdot[CO]/d[HCCO] */
    J[772] += dqdci;              /* dwdot[SXCH2]/d[HCCO] */
    J[775] -= dqdci;              /* dwdot[HCCO]/d[HCCO] */
    /* d()/dT */
    J[871] -= dqdT;               /* dwdot[H]/dT */
    J[879] += dqdT;               /* dwdot[CO]/dT */
    J[892] += dqdT;               /* dwdot[SXCH2]/dT */
    J[895] -= dqdT;               /* dwdot[HCCO]/dT */

    /*reaction 144: CH2CO + H <=> CH3 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[16];
    k_f = prefactor_units[143] * fwd_A[143]
                * exp(fwd_beta[143] * tc[0] - activation_units[143] * fwd_Ea[143] * invT);
    dlnkfdT = fwd_beta[143] * invT + activation_units[143] * fwd_Ea[143] * invT2;
    /* reverse */
    phi_r = sc[9]*sc[10];
    Kc = exp(g_RT[1] - g_RT[9] - g_RT[10] + g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[16]) + (h_RT[9] + h_RT[10]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[9] += q; /* CO */
    wdot[10] += q; /* CH3 */
    wdot[16] -= q; /* CH2CO */
    /* d()/d[H] */
    dqdci =  + k_f*sc[16];
    J[31] -= dqdci;               /* dwdot[H]/d[H] */
    J[39] += dqdci;               /* dwdot[CO]/d[H] */
    J[40] += dqdci;               /* dwdot[CH3]/d[H] */
    J[46] -= dqdci;               /* dwdot[CH2CO]/d[H] */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[10];
    J[271] -= dqdci;              /* dwdot[H]/d[CO] */
    J[279] += dqdci;              /* dwdot[CO]/d[CO] */
    J[280] += dqdci;              /* dwdot[CH3]/d[CO] */
    J[286] -= dqdci;              /* dwdot[CH2CO]/d[CO] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[9];
    J[301] -= dqdci;              /* dwdot[H]/d[CH3] */
    J[309] += dqdci;              /* dwdot[CO]/d[CH3] */
    J[310] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[316] -= dqdci;              /* dwdot[CH2CO]/d[CH3] */
    /* d()/d[CH2CO] */
    dqdci =  + k_f*sc[1];
    J[481] -= dqdci;              /* dwdot[H]/d[CH2CO] */
    J[489] += dqdci;              /* dwdot[CO]/d[CH2CO] */
    J[490] += dqdci;              /* dwdot[CH3]/d[CH2CO] */
    J[496] -= dqdci;              /* dwdot[CH2CO]/d[CH2CO] */
    /* d()/dT */
    J[871] -= dqdT;               /* dwdot[H]/dT */
    J[879] += dqdT;               /* dwdot[CO]/dT */
    J[880] += dqdT;               /* dwdot[CH3]/dT */
    J[886] -= dqdT;               /* dwdot[CH2CO]/dT */

    /*reaction 145: CH2CO + TXCH2 <=> C2H4 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[16]*sc[21];
    k_f = prefactor_units[144] * fwd_A[144]
                * exp(fwd_beta[144] * tc[0] - activation_units[144] * fwd_Ea[144] * invT);
    dlnkfdT = fwd_beta[144] * invT + activation_units[144] * fwd_Ea[144] * invT2;
    /* reverse */
    phi_r = sc[9]*sc[15];
    Kc = exp(-g_RT[9] - g_RT[15] + g_RT[16] + g_RT[21]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[16] + h_RT[21]) + (h_RT[9] + h_RT[15]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[9] += q; /* CO */
    wdot[15] += q; /* C2H4 */
    wdot[16] -= q; /* CH2CO */
    wdot[21] -= q; /* TXCH2 */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[15];
    J[279] += dqdci;              /* dwdot[CO]/d[CO] */
    J[285] += dqdci;              /* dwdot[C2H4]/d[CO] */
    J[286] -= dqdci;              /* dwdot[CH2CO]/d[CO] */
    J[291] -= dqdci;              /* dwdot[TXCH2]/d[CO] */
    /* d()/d[C2H4] */
    dqdci =  - k_r*sc[9];
    J[459] += dqdci;              /* dwdot[CO]/d[C2H4] */
    J[465] += dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[466] -= dqdci;              /* dwdot[CH2CO]/d[C2H4] */
    J[471] -= dqdci;              /* dwdot[TXCH2]/d[C2H4] */
    /* d()/d[CH2CO] */
    dqdci =  + k_f*sc[21];
    J[489] += dqdci;              /* dwdot[CO]/d[CH2CO] */
    J[495] += dqdci;              /* dwdot[C2H4]/d[CH2CO] */
    J[496] -= dqdci;              /* dwdot[CH2CO]/d[CH2CO] */
    J[501] -= dqdci;              /* dwdot[TXCH2]/d[CH2CO] */
    /* d()/d[TXCH2] */
    dqdci =  + k_f*sc[16];
    J[639] += dqdci;              /* dwdot[CO]/d[TXCH2] */
    J[645] += dqdci;              /* dwdot[C2H4]/d[TXCH2] */
    J[646] -= dqdci;              /* dwdot[CH2CO]/d[TXCH2] */
    J[651] -= dqdci;              /* dwdot[TXCH2]/d[TXCH2] */
    /* d()/dT */
    J[879] += dqdT;               /* dwdot[CO]/dT */
    J[885] += dqdT;               /* dwdot[C2H4]/dT */
    J[886] -= dqdT;               /* dwdot[CH2CO]/dT */
    J[891] -= dqdT;               /* dwdot[TXCH2]/dT */

    /*reaction 146: CH2CO + O <=> HCCO + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[16];
    k_f = prefactor_units[145] * fwd_A[145]
                * exp(fwd_beta[145] * tc[0] - activation_units[145] * fwd_Ea[145] * invT);
    dlnkfdT = fwd_beta[145] * invT + activation_units[145] * fwd_Ea[145] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[25];
    Kc = exp(g_RT[3] - g_RT[4] + g_RT[16] - g_RT[25]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[16]) + (h_RT[4] + h_RT[25]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O */
    wdot[4] += q; /* OH */
    wdot[16] -= q; /* CH2CO */
    wdot[25] += q; /* HCCO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[16];
    J[93] -= dqdci;               /* dwdot[O]/d[O] */
    J[94] += dqdci;               /* dwdot[OH]/d[O] */
    J[106] -= dqdci;              /* dwdot[CH2CO]/d[O] */
    J[115] += dqdci;              /* dwdot[HCCO]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[25];
    J[123] -= dqdci;              /* dwdot[O]/d[OH] */
    J[124] += dqdci;              /* dwdot[OH]/d[OH] */
    J[136] -= dqdci;              /* dwdot[CH2CO]/d[OH] */
    J[145] += dqdci;              /* dwdot[HCCO]/d[OH] */
    /* d()/d[CH2CO] */
    dqdci =  + k_f*sc[3];
    J[483] -= dqdci;              /* dwdot[O]/d[CH2CO] */
    J[484] += dqdci;              /* dwdot[OH]/d[CH2CO] */
    J[496] -= dqdci;              /* dwdot[CH2CO]/d[CH2CO] */
    J[505] += dqdci;              /* dwdot[HCCO]/d[CH2CO] */
    /* d()/d[HCCO] */
    dqdci =  - k_r*sc[4];
    J[753] -= dqdci;              /* dwdot[O]/d[HCCO] */
    J[754] += dqdci;              /* dwdot[OH]/d[HCCO] */
    J[766] -= dqdci;              /* dwdot[CH2CO]/d[HCCO] */
    J[775] += dqdci;              /* dwdot[HCCO]/d[HCCO] */
    /* d()/dT */
    J[873] -= dqdT;               /* dwdot[O]/dT */
    J[874] += dqdT;               /* dwdot[OH]/dT */
    J[886] -= dqdT;               /* dwdot[CH2CO]/dT */
    J[895] += dqdT;               /* dwdot[HCCO]/dT */

    /*reaction 147: CH2CO + CH3 <=> HCCO + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[10]*sc[16];
    k_f = prefactor_units[146] * fwd_A[146]
                * exp(fwd_beta[146] * tc[0] - activation_units[146] * fwd_Ea[146] * invT);
    dlnkfdT = fwd_beta[146] * invT + activation_units[146] * fwd_Ea[146] * invT2;
    /* reverse */
    phi_r = sc[13]*sc[25];
    Kc = exp(g_RT[10] - g_RT[13] + g_RT[16] - g_RT[25]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[10] + h_RT[16]) + (h_RT[13] + h_RT[25]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[10] -= q; /* CH3 */
    wdot[13] += q; /* CH4 */
    wdot[16] -= q; /* CH2CO */
    wdot[25] += q; /* HCCO */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[16];
    J[310] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[313] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[316] -= dqdci;              /* dwdot[CH2CO]/d[CH3] */
    J[325] += dqdci;              /* dwdot[HCCO]/d[CH3] */
    /* d()/d[CH4] */
    dqdci =  - k_r*sc[25];
    J[400] -= dqdci;              /* dwdot[CH3]/d[CH4] */
    J[403] += dqdci;              /* dwdot[CH4]/d[CH4] */
    J[406] -= dqdci;              /* dwdot[CH2CO]/d[CH4] */
    J[415] += dqdci;              /* dwdot[HCCO]/d[CH4] */
    /* d()/d[CH2CO] */
    dqdci =  + k_f*sc[10];
    J[490] -= dqdci;              /* dwdot[CH3]/d[CH2CO] */
    J[493] += dqdci;              /* dwdot[CH4]/d[CH2CO] */
    J[496] -= dqdci;              /* dwdot[CH2CO]/d[CH2CO] */
    J[505] += dqdci;              /* dwdot[HCCO]/d[CH2CO] */
    /* d()/d[HCCO] */
    dqdci =  - k_r*sc[13];
    J[760] -= dqdci;              /* dwdot[CH3]/d[HCCO] */
    J[763] += dqdci;              /* dwdot[CH4]/d[HCCO] */
    J[766] -= dqdci;              /* dwdot[CH2CO]/d[HCCO] */
    J[775] += dqdci;              /* dwdot[HCCO]/d[HCCO] */
    /* d()/dT */
    J[880] -= dqdT;               /* dwdot[CH3]/dT */
    J[883] += dqdT;               /* dwdot[CH4]/dT */
    J[886] -= dqdT;               /* dwdot[CH2CO]/dT */
    J[895] += dqdT;               /* dwdot[HCCO]/dT */

    /*reaction 148: CH2CO + O <=> TXCH2 + CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[16];
    k_f = prefactor_units[147] * fwd_A[147]
                * exp(fwd_beta[147] * tc[0] - activation_units[147] * fwd_Ea[147] * invT);
    dlnkfdT = fwd_beta[147] * invT + activation_units[147] * fwd_Ea[147] * invT2;
    /* reverse */
    phi_r = sc[12]*sc[21];
    Kc = exp(g_RT[3] - g_RT[12] + g_RT[16] - g_RT[21]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[16]) + (h_RT[12] + h_RT[21]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O */
    wdot[12] += q; /* CO2 */
    wdot[16] -= q; /* CH2CO */
    wdot[21] += q; /* TXCH2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[16];
    J[93] -= dqdci;               /* dwdot[O]/d[O] */
    J[102] += dqdci;              /* dwdot[CO2]/d[O] */
    J[106] -= dqdci;              /* dwdot[CH2CO]/d[O] */
    J[111] += dqdci;              /* dwdot[TXCH2]/d[O] */
    /* d()/d[CO2] */
    dqdci =  - k_r*sc[21];
    J[363] -= dqdci;              /* dwdot[O]/d[CO2] */
    J[372] += dqdci;              /* dwdot[CO2]/d[CO2] */
    J[376] -= dqdci;              /* dwdot[CH2CO]/d[CO2] */
    J[381] += dqdci;              /* dwdot[TXCH2]/d[CO2] */
    /* d()/d[CH2CO] */
    dqdci =  + k_f*sc[3];
    J[483] -= dqdci;              /* dwdot[O]/d[CH2CO] */
    J[492] += dqdci;              /* dwdot[CO2]/d[CH2CO] */
    J[496] -= dqdci;              /* dwdot[CH2CO]/d[CH2CO] */
    J[501] += dqdci;              /* dwdot[TXCH2]/d[CH2CO] */
    /* d()/d[TXCH2] */
    dqdci =  - k_r*sc[12];
    J[633] -= dqdci;              /* dwdot[O]/d[TXCH2] */
    J[642] += dqdci;              /* dwdot[CO2]/d[TXCH2] */
    J[646] -= dqdci;              /* dwdot[CH2CO]/d[TXCH2] */
    J[651] += dqdci;              /* dwdot[TXCH2]/d[TXCH2] */
    /* d()/dT */
    J[873] -= dqdT;               /* dwdot[O]/dT */
    J[882] += dqdT;               /* dwdot[CO2]/dT */
    J[886] -= dqdT;               /* dwdot[CH2CO]/dT */
    J[891] += dqdT;               /* dwdot[TXCH2]/dT */

    /*reaction 149: CH2CO + CH3 <=> C2H5 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[10]*sc[16];
    k_f = prefactor_units[148] * fwd_A[148]
                * exp(fwd_beta[148] * tc[0] - activation_units[148] * fwd_Ea[148] * invT);
    dlnkfdT = fwd_beta[148] * invT + activation_units[148] * fwd_Ea[148] * invT2;
    /* reverse */
    phi_r = sc[9]*sc[24];
    Kc = exp(-g_RT[9] + g_RT[10] + g_RT[16] - g_RT[24]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[10] + h_RT[16]) + (h_RT[9] + h_RT[24]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[9] += q; /* CO */
    wdot[10] -= q; /* CH3 */
    wdot[16] -= q; /* CH2CO */
    wdot[24] += q; /* C2H5 */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[24];
    J[279] += dqdci;              /* dwdot[CO]/d[CO] */
    J[280] -= dqdci;              /* dwdot[CH3]/d[CO] */
    J[286] -= dqdci;              /* dwdot[CH2CO]/d[CO] */
    J[294] += dqdci;              /* dwdot[C2H5]/d[CO] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[16];
    J[309] += dqdci;              /* dwdot[CO]/d[CH3] */
    J[310] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[316] -= dqdci;              /* dwdot[CH2CO]/d[CH3] */
    J[324] += dqdci;              /* dwdot[C2H5]/d[CH3] */
    /* d()/d[CH2CO] */
    dqdci =  + k_f*sc[10];
    J[489] += dqdci;              /* dwdot[CO]/d[CH2CO] */
    J[490] -= dqdci;              /* dwdot[CH3]/d[CH2CO] */
    J[496] -= dqdci;              /* dwdot[CH2CO]/d[CH2CO] */
    J[504] += dqdci;              /* dwdot[C2H5]/d[CH2CO] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[9];
    J[729] += dqdci;              /* dwdot[CO]/d[C2H5] */
    J[730] -= dqdci;              /* dwdot[CH3]/d[C2H5] */
    J[736] -= dqdci;              /* dwdot[CH2CO]/d[C2H5] */
    J[744] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    /* d()/dT */
    J[879] += dqdT;               /* dwdot[CO]/dT */
    J[880] -= dqdT;               /* dwdot[CH3]/dT */
    J[886] -= dqdT;               /* dwdot[CH2CO]/dT */
    J[894] += dqdT;               /* dwdot[C2H5]/dT */

    /*reaction 150: CH2CO + OH <=> HCCO + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[16];
    k_f = prefactor_units[149] * fwd_A[149]
                * exp(fwd_beta[149] * tc[0] - activation_units[149] * fwd_Ea[149] * invT);
    dlnkfdT = fwd_beta[149] * invT + activation_units[149] * fwd_Ea[149] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[25];
    Kc = exp(g_RT[4] - g_RT[7] + g_RT[16] - g_RT[25]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[16]) + (h_RT[7] + h_RT[25]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[7] += q; /* H2O */
    wdot[16] -= q; /* CH2CO */
    wdot[25] += q; /* HCCO */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[16];
    J[124] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[127] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[136] -= dqdci;              /* dwdot[CH2CO]/d[OH] */
    J[145] += dqdci;              /* dwdot[HCCO]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[25];
    J[214] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[217] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[226] -= dqdci;              /* dwdot[CH2CO]/d[H2O] */
    J[235] += dqdci;              /* dwdot[HCCO]/d[H2O] */
    /* d()/d[CH2CO] */
    dqdci =  + k_f*sc[4];
    J[484] -= dqdci;              /* dwdot[OH]/d[CH2CO] */
    J[487] += dqdci;              /* dwdot[H2O]/d[CH2CO] */
    J[496] -= dqdci;              /* dwdot[CH2CO]/d[CH2CO] */
    J[505] += dqdci;              /* dwdot[HCCO]/d[CH2CO] */
    /* d()/d[HCCO] */
    dqdci =  - k_r*sc[7];
    J[754] -= dqdci;              /* dwdot[OH]/d[HCCO] */
    J[757] += dqdci;              /* dwdot[H2O]/d[HCCO] */
    J[766] -= dqdci;              /* dwdot[CH2CO]/d[HCCO] */
    J[775] += dqdci;              /* dwdot[HCCO]/d[HCCO] */
    /* d()/dT */
    J[874] -= dqdT;               /* dwdot[OH]/dT */
    J[877] += dqdT;               /* dwdot[H2O]/dT */
    J[886] -= dqdT;               /* dwdot[CH2CO]/dT */
    J[895] += dqdT;               /* dwdot[HCCO]/dT */

    /*reaction 151: CH2CO + H <=> HCCO + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[16];
    k_f = prefactor_units[150] * fwd_A[150]
                * exp(fwd_beta[150] * tc[0] - activation_units[150] * fwd_Ea[150] * invT);
    dlnkfdT = fwd_beta[150] * invT + activation_units[150] * fwd_Ea[150] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[25];
    Kc = exp(g_RT[1] - g_RT[2] + g_RT[16] - g_RT[25]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[16]) + (h_RT[2] + h_RT[25]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[2] += q; /* H2 */
    wdot[16] -= q; /* CH2CO */
    wdot[25] += q; /* HCCO */
    /* d()/d[H] */
    dqdci =  + k_f*sc[16];
    J[31] -= dqdci;               /* dwdot[H]/d[H] */
    J[32] += dqdci;               /* dwdot[H2]/d[H] */
    J[46] -= dqdci;               /* dwdot[CH2CO]/d[H] */
    J[55] += dqdci;               /* dwdot[HCCO]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[25];
    J[61] -= dqdci;               /* dwdot[H]/d[H2] */
    J[62] += dqdci;               /* dwdot[H2]/d[H2] */
    J[76] -= dqdci;               /* dwdot[CH2CO]/d[H2] */
    J[85] += dqdci;               /* dwdot[HCCO]/d[H2] */
    /* d()/d[CH2CO] */
    dqdci =  + k_f*sc[1];
    J[481] -= dqdci;              /* dwdot[H]/d[CH2CO] */
    J[482] += dqdci;              /* dwdot[H2]/d[CH2CO] */
    J[496] -= dqdci;              /* dwdot[CH2CO]/d[CH2CO] */
    J[505] += dqdci;              /* dwdot[HCCO]/d[CH2CO] */
    /* d()/d[HCCO] */
    dqdci =  - k_r*sc[2];
    J[751] -= dqdci;              /* dwdot[H]/d[HCCO] */
    J[752] += dqdci;              /* dwdot[H2]/d[HCCO] */
    J[766] -= dqdci;              /* dwdot[CH2CO]/d[HCCO] */
    J[775] += dqdci;              /* dwdot[HCCO]/d[HCCO] */
    /* d()/dT */
    J[871] -= dqdT;               /* dwdot[H]/dT */
    J[872] += dqdT;               /* dwdot[H2]/dT */
    J[886] -= dqdT;               /* dwdot[CH2CO]/dT */
    J[895] += dqdT;               /* dwdot[HCCO]/dT */

    /*reaction 152: CH2CO + TXCH2 <=> HCCO + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[16]*sc[21];
    k_f = prefactor_units[151] * fwd_A[151]
                * exp(fwd_beta[151] * tc[0] - activation_units[151] * fwd_Ea[151] * invT);
    dlnkfdT = fwd_beta[151] * invT + activation_units[151] * fwd_Ea[151] * invT2;
    /* reverse */
    phi_r = sc[10]*sc[25];
    Kc = exp(-g_RT[10] + g_RT[16] + g_RT[21] - g_RT[25]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[16] + h_RT[21]) + (h_RT[10] + h_RT[25]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[10] += q; /* CH3 */
    wdot[16] -= q; /* CH2CO */
    wdot[21] -= q; /* TXCH2 */
    wdot[25] += q; /* HCCO */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[25];
    J[310] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[316] -= dqdci;              /* dwdot[CH2CO]/d[CH3] */
    J[321] -= dqdci;              /* dwdot[TXCH2]/d[CH3] */
    J[325] += dqdci;              /* dwdot[HCCO]/d[CH3] */
    /* d()/d[CH2CO] */
    dqdci =  + k_f*sc[21];
    J[490] += dqdci;              /* dwdot[CH3]/d[CH2CO] */
    J[496] -= dqdci;              /* dwdot[CH2CO]/d[CH2CO] */
    J[501] -= dqdci;              /* dwdot[TXCH2]/d[CH2CO] */
    J[505] += dqdci;              /* dwdot[HCCO]/d[CH2CO] */
    /* d()/d[TXCH2] */
    dqdci =  + k_f*sc[16];
    J[640] += dqdci;              /* dwdot[CH3]/d[TXCH2] */
    J[646] -= dqdci;              /* dwdot[CH2CO]/d[TXCH2] */
    J[651] -= dqdci;              /* dwdot[TXCH2]/d[TXCH2] */
    J[655] += dqdci;              /* dwdot[HCCO]/d[TXCH2] */
    /* d()/d[HCCO] */
    dqdci =  - k_r*sc[10];
    J[760] += dqdci;              /* dwdot[CH3]/d[HCCO] */
    J[766] -= dqdci;              /* dwdot[CH2CO]/d[HCCO] */
    J[771] -= dqdci;              /* dwdot[TXCH2]/d[HCCO] */
    J[775] += dqdci;              /* dwdot[HCCO]/d[HCCO] */
    /* d()/dT */
    J[880] += dqdT;               /* dwdot[CH3]/dT */
    J[886] -= dqdT;               /* dwdot[CH2CO]/dT */
    J[891] -= dqdT;               /* dwdot[TXCH2]/dT */
    J[895] += dqdT;               /* dwdot[HCCO]/dT */

    /*reaction 153: CH2CHO + O <=> CH2O + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[27];
    k_f = prefactor_units[152] * fwd_A[152]
                * exp(fwd_beta[152] * tc[0] - activation_units[152] * fwd_Ea[152] * invT);
    dlnkfdT = fwd_beta[152] * invT + activation_units[152] * fwd_Ea[152] * invT2;
    /* reverse */
    phi_r = sc[11]*sc[20];
    Kc = exp(g_RT[3] - g_RT[11] - g_RT[20] + g_RT[27]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[27]) + (h_RT[11] + h_RT[20]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O */
    wdot[11] += q; /* CH2O */
    wdot[20] += q; /* HCO */
    wdot[27] -= q; /* CH2CHO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[27];
    J[93] -= dqdci;               /* dwdot[O]/d[O] */
    J[101] += dqdci;              /* dwdot[CH2O]/d[O] */
    J[110] += dqdci;              /* dwdot[HCO]/d[O] */
    J[117] -= dqdci;              /* dwdot[CH2CHO]/d[O] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[20];
    J[333] -= dqdci;              /* dwdot[O]/d[CH2O] */
    J[341] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[350] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[357] -= dqdci;              /* dwdot[CH2CHO]/d[CH2O] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[11];
    J[603] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[611] += dqdci;              /* dwdot[CH2O]/d[HCO] */
    J[620] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[627] -= dqdci;              /* dwdot[CH2CHO]/d[HCO] */
    /* d()/d[CH2CHO] */
    dqdci =  + k_f*sc[3];
    J[813] -= dqdci;              /* dwdot[O]/d[CH2CHO] */
    J[821] += dqdci;              /* dwdot[CH2O]/d[CH2CHO] */
    J[830] += dqdci;              /* dwdot[HCO]/d[CH2CHO] */
    J[837] -= dqdci;              /* dwdot[CH2CHO]/d[CH2CHO] */
    /* d()/dT */
    J[873] -= dqdT;               /* dwdot[O]/dT */
    J[881] += dqdT;               /* dwdot[CH2O]/dT */
    J[890] += dqdT;               /* dwdot[HCO]/dT */
    J[897] -= dqdT;               /* dwdot[CH2CHO]/dT */

    /*reaction 154: CH2CHO <=> CH2CO + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[27];
    k_f = prefactor_units[153] * fwd_A[153]
                * exp(fwd_beta[153] * tc[0] - activation_units[153] * fwd_Ea[153] * invT);
    dlnkfdT = fwd_beta[153] * invT + activation_units[153] * fwd_Ea[153] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[16];
    Kc = refC * exp(-g_RT[1] - g_RT[16] + g_RT[27]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[27]) + (h_RT[1] + h_RT[16]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[16] += q; /* CH2CO */
    wdot[27] -= q; /* CH2CHO */
    /* d()/d[H] */
    dqdci =  - k_r*sc[16];
    J[31] += dqdci;               /* dwdot[H]/d[H] */
    J[46] += dqdci;               /* dwdot[CH2CO]/d[H] */
    J[57] -= dqdci;               /* dwdot[CH2CHO]/d[H] */
    /* d()/d[CH2CO] */
    dqdci =  - k_r*sc[1];
    J[481] += dqdci;              /* dwdot[H]/d[CH2CO] */
    J[496] += dqdci;              /* dwdot[CH2CO]/d[CH2CO] */
    J[507] -= dqdci;              /* dwdot[CH2CHO]/d[CH2CO] */
    /* d()/d[CH2CHO] */
    dqdci =  + k_f;
    J[811] += dqdci;              /* dwdot[H]/d[CH2CHO] */
    J[826] += dqdci;              /* dwdot[CH2CO]/d[CH2CHO] */
    J[837] -= dqdci;              /* dwdot[CH2CHO]/d[CH2CHO] */
    /* d()/dT */
    J[871] += dqdT;               /* dwdot[H]/dT */
    J[886] += dqdT;               /* dwdot[CH2CO]/dT */
    J[897] -= dqdT;               /* dwdot[CH2CHO]/dT */

    /*reaction 155: CH2CHO + OH <=> H2O + CH2CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[27];
    k_f = prefactor_units[154] * fwd_A[154]
                * exp(fwd_beta[154] * tc[0] - activation_units[154] * fwd_Ea[154] * invT);
    dlnkfdT = fwd_beta[154] * invT + activation_units[154] * fwd_Ea[154] * invT2;
    /* reverse */
    phi_r = sc[7]*sc[16];
    Kc = exp(g_RT[4] - g_RT[7] - g_RT[16] + g_RT[27]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[27]) + (h_RT[7] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[7] += q; /* H2O */
    wdot[16] += q; /* CH2CO */
    wdot[27] -= q; /* CH2CHO */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[27];
    J[124] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[127] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[136] += dqdci;              /* dwdot[CH2CO]/d[OH] */
    J[147] -= dqdci;              /* dwdot[CH2CHO]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[16];
    J[214] -= dqdci;              /* dwdot[OH]/d[H2O] */
    J[217] += dqdci;              /* dwdot[H2O]/d[H2O] */
    J[226] += dqdci;              /* dwdot[CH2CO]/d[H2O] */
    J[237] -= dqdci;              /* dwdot[CH2CHO]/d[H2O] */
    /* d()/d[CH2CO] */
    dqdci =  - k_r*sc[7];
    J[484] -= dqdci;              /* dwdot[OH]/d[CH2CO] */
    J[487] += dqdci;              /* dwdot[H2O]/d[CH2CO] */
    J[496] += dqdci;              /* dwdot[CH2CO]/d[CH2CO] */
    J[507] -= dqdci;              /* dwdot[CH2CHO]/d[CH2CO] */
    /* d()/d[CH2CHO] */
    dqdci =  + k_f*sc[4];
    J[814] -= dqdci;              /* dwdot[OH]/d[CH2CHO] */
    J[817] += dqdci;              /* dwdot[H2O]/d[CH2CHO] */
    J[826] += dqdci;              /* dwdot[CH2CO]/d[CH2CHO] */
    J[837] -= dqdci;              /* dwdot[CH2CHO]/d[CH2CHO] */
    /* d()/dT */
    J[874] -= dqdT;               /* dwdot[OH]/dT */
    J[877] += dqdT;               /* dwdot[H2O]/dT */
    J[886] += dqdT;               /* dwdot[CH2CO]/dT */
    J[897] -= dqdT;               /* dwdot[CH2CHO]/dT */

    /*reaction 156: CH2CHO + H <=> CH2CO + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[27];
    k_f = prefactor_units[155] * fwd_A[155]
                * exp(fwd_beta[155] * tc[0] - activation_units[155] * fwd_Ea[155] * invT);
    dlnkfdT = fwd_beta[155] * invT + activation_units[155] * fwd_Ea[155] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[16];
    Kc = exp(g_RT[1] - g_RT[2] - g_RT[16] + g_RT[27]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[27]) + (h_RT[2] + h_RT[16]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[2] += q; /* H2 */
    wdot[16] += q; /* CH2CO */
    wdot[27] -= q; /* CH2CHO */
    /* d()/d[H] */
    dqdci =  + k_f*sc[27];
    J[31] -= dqdci;               /* dwdot[H]/d[H] */
    J[32] += dqdci;               /* dwdot[H2]/d[H] */
    J[46] += dqdci;               /* dwdot[CH2CO]/d[H] */
    J[57] -= dqdci;               /* dwdot[CH2CHO]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[16];
    J[61] -= dqdci;               /* dwdot[H]/d[H2] */
    J[62] += dqdci;               /* dwdot[H2]/d[H2] */
    J[76] += dqdci;               /* dwdot[CH2CO]/d[H2] */
    J[87] -= dqdci;               /* dwdot[CH2CHO]/d[H2] */
    /* d()/d[CH2CO] */
    dqdci =  - k_r*sc[2];
    J[481] -= dqdci;              /* dwdot[H]/d[CH2CO] */
    J[482] += dqdci;              /* dwdot[H2]/d[CH2CO] */
    J[496] += dqdci;              /* dwdot[CH2CO]/d[CH2CO] */
    J[507] -= dqdci;              /* dwdot[CH2CHO]/d[CH2CO] */
    /* d()/d[CH2CHO] */
    dqdci =  + k_f*sc[1];
    J[811] -= dqdci;              /* dwdot[H]/d[CH2CHO] */
    J[812] += dqdci;              /* dwdot[H2]/d[CH2CHO] */
    J[826] += dqdci;              /* dwdot[CH2CO]/d[CH2CHO] */
    J[837] -= dqdci;              /* dwdot[CH2CHO]/d[CH2CHO] */
    /* d()/dT */
    J[871] -= dqdT;               /* dwdot[H]/dT */
    J[872] += dqdT;               /* dwdot[H2]/dT */
    J[886] += dqdT;               /* dwdot[CH2CO]/dT */
    J[897] -= dqdT;               /* dwdot[CH2CHO]/dT */

    /*reaction 157: CH2CHO + O2 => OH + CO + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[27];
    k_f = prefactor_units[156] * fwd_A[156]
                * exp(fwd_beta[156] * tc[0] - activation_units[156] * fwd_Ea[156] * invT);
    dlnkfdT = fwd_beta[156] * invT + activation_units[156] * fwd_Ea[156] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[4] += q; /* OH */
    wdot[5] -= q; /* O2 */
    wdot[9] += q; /* CO */
    wdot[11] += q; /* CH2O */
    wdot[27] -= q; /* CH2CHO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[27];
    J[154] += dqdci;              /* dwdot[OH]/d[O2] */
    J[155] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[159] += dqdci;              /* dwdot[CO]/d[O2] */
    J[161] += dqdci;              /* dwdot[CH2O]/d[O2] */
    J[177] -= dqdci;              /* dwdot[CH2CHO]/d[O2] */
    /* d()/d[CH2CHO] */
    dqdci =  + k_f*sc[5];
    J[814] += dqdci;              /* dwdot[OH]/d[CH2CHO] */
    J[815] -= dqdci;              /* dwdot[O2]/d[CH2CHO] */
    J[819] += dqdci;              /* dwdot[CO]/d[CH2CHO] */
    J[821] += dqdci;              /* dwdot[CH2O]/d[CH2CHO] */
    J[837] -= dqdci;              /* dwdot[CH2CHO]/d[CH2CHO] */
    /* d()/dT */
    J[874] += dqdT;               /* dwdot[OH]/dT */
    J[875] -= dqdT;               /* dwdot[O2]/dT */
    J[879] += dqdT;               /* dwdot[CO]/dT */
    J[881] += dqdT;               /* dwdot[CH2O]/dT */
    J[897] -= dqdT;               /* dwdot[CH2CHO]/dT */

    /*reaction 158: CH2CHO <=> CH3 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[27];
    k_f = prefactor_units[157] * fwd_A[157]
                * exp(fwd_beta[157] * tc[0] - activation_units[157] * fwd_Ea[157] * invT);
    dlnkfdT = fwd_beta[157] * invT + activation_units[157] * fwd_Ea[157] * invT2;
    /* reverse */
    phi_r = sc[9]*sc[10];
    Kc = refC * exp(-g_RT[9] - g_RT[10] + g_RT[27]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[27]) + (h_RT[9] + h_RT[10]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[9] += q; /* CO */
    wdot[10] += q; /* CH3 */
    wdot[27] -= q; /* CH2CHO */
    /* d()/d[CO] */
    dqdci =  - k_r*sc[10];
    J[279] += dqdci;              /* dwdot[CO]/d[CO] */
    J[280] += dqdci;              /* dwdot[CH3]/d[CO] */
    J[297] -= dqdci;              /* dwdot[CH2CHO]/d[CO] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[9];
    J[309] += dqdci;              /* dwdot[CO]/d[CH3] */
    J[310] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[327] -= dqdci;              /* dwdot[CH2CHO]/d[CH3] */
    /* d()/d[CH2CHO] */
    dqdci =  + k_f;
    J[819] += dqdci;              /* dwdot[CO]/d[CH2CHO] */
    J[820] += dqdci;              /* dwdot[CH3]/d[CH2CHO] */
    J[837] -= dqdci;              /* dwdot[CH2CHO]/d[CH2CHO] */
    /* d()/dT */
    J[879] += dqdT;               /* dwdot[CO]/dT */
    J[880] += dqdT;               /* dwdot[CH3]/dT */
    J[897] -= dqdT;               /* dwdot[CH2CHO]/dT */

    /*reaction 159: CH2CHO + O2 => OH + 2.000000 HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[27];
    k_f = prefactor_units[158] * fwd_A[158]
                * exp(fwd_beta[158] * tc[0] - activation_units[158] * fwd_Ea[158] * invT);
    dlnkfdT = fwd_beta[158] * invT + activation_units[158] * fwd_Ea[158] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[4] += q; /* OH */
    wdot[5] -= q; /* O2 */
    wdot[20] += 2 * q; /* HCO */
    wdot[27] -= q; /* CH2CHO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[27];
    J[154] += dqdci;              /* dwdot[OH]/d[O2] */
    J[155] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[170] += 2 * dqdci;          /* dwdot[HCO]/d[O2] */
    J[177] -= dqdci;              /* dwdot[CH2CHO]/d[O2] */
    /* d()/d[CH2CHO] */
    dqdci =  + k_f*sc[5];
    J[814] += dqdci;              /* dwdot[OH]/d[CH2CHO] */
    J[815] -= dqdci;              /* dwdot[O2]/d[CH2CHO] */
    J[830] += 2 * dqdci;          /* dwdot[HCO]/d[CH2CHO] */
    J[837] -= dqdci;              /* dwdot[CH2CHO]/d[CH2CHO] */
    /* d()/dT */
    J[874] += dqdT;               /* dwdot[OH]/dT */
    J[875] -= dqdT;               /* dwdot[O2]/dT */
    J[890] += 2 * dqdT;           /* dwdot[HCO]/dT */
    J[897] -= dqdT;               /* dwdot[CH2CHO]/dT */

    /*reaction 160: CH2CHO + H <=> CH3 + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[27];
    k_f = prefactor_units[159] * fwd_A[159]
                * exp(fwd_beta[159] * tc[0] - activation_units[159] * fwd_Ea[159] * invT);
    dlnkfdT = fwd_beta[159] * invT + activation_units[159] * fwd_Ea[159] * invT2;
    /* reverse */
    phi_r = sc[10]*sc[20];
    Kc = exp(g_RT[1] - g_RT[10] - g_RT[20] + g_RT[27]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[27]) + (h_RT[10] + h_RT[20]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[10] += q; /* CH3 */
    wdot[20] += q; /* HCO */
    wdot[27] -= q; /* CH2CHO */
    /* d()/d[H] */
    dqdci =  + k_f*sc[27];
    J[31] -= dqdci;               /* dwdot[H]/d[H] */
    J[40] += dqdci;               /* dwdot[CH3]/d[H] */
    J[50] += dqdci;               /* dwdot[HCO]/d[H] */
    J[57] -= dqdci;               /* dwdot[CH2CHO]/d[H] */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[20];
    J[301] -= dqdci;              /* dwdot[H]/d[CH3] */
    J[310] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[320] += dqdci;              /* dwdot[HCO]/d[CH3] */
    J[327] -= dqdci;              /* dwdot[CH2CHO]/d[CH3] */
    /* d()/d[HCO] */
    dqdci =  - k_r*sc[10];
    J[601] -= dqdci;              /* dwdot[H]/d[HCO] */
    J[610] += dqdci;              /* dwdot[CH3]/d[HCO] */
    J[620] += dqdci;              /* dwdot[HCO]/d[HCO] */
    J[627] -= dqdci;              /* dwdot[CH2CHO]/d[HCO] */
    /* d()/d[CH2CHO] */
    dqdci =  + k_f*sc[1];
    J[811] -= dqdci;              /* dwdot[H]/d[CH2CHO] */
    J[820] += dqdci;              /* dwdot[CH3]/d[CH2CHO] */
    J[830] += dqdci;              /* dwdot[HCO]/d[CH2CHO] */
    J[837] -= dqdci;              /* dwdot[CH2CHO]/d[CH2CHO] */
    /* d()/dT */
    J[871] -= dqdT;               /* dwdot[H]/dT */
    J[880] += dqdT;               /* dwdot[CH3]/dT */
    J[890] += dqdT;               /* dwdot[HCO]/dT */
    J[897] -= dqdT;               /* dwdot[CH2CHO]/dT */

    /*reaction 161: CH3CHO + O => CH3 + CO + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[26];
    k_f = prefactor_units[160] * fwd_A[160]
                * exp(fwd_beta[160] * tc[0] - activation_units[160] * fwd_Ea[160] * invT);
    dlnkfdT = fwd_beta[160] * invT + activation_units[160] * fwd_Ea[160] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* O */
    wdot[4] += q; /* OH */
    wdot[9] += q; /* CO */
    wdot[10] += q; /* CH3 */
    wdot[26] -= q; /* CH3CHO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[26];
    J[93] -= dqdci;               /* dwdot[O]/d[O] */
    J[94] += dqdci;               /* dwdot[OH]/d[O] */
    J[99] += dqdci;               /* dwdot[CO]/d[O] */
    J[100] += dqdci;              /* dwdot[CH3]/d[O] */
    J[116] -= dqdci;              /* dwdot[CH3CHO]/d[O] */
    /* d()/d[CH3CHO] */
    dqdci =  + k_f*sc[3];
    J[783] -= dqdci;              /* dwdot[O]/d[CH3CHO] */
    J[784] += dqdci;              /* dwdot[OH]/d[CH3CHO] */
    J[789] += dqdci;              /* dwdot[CO]/d[CH3CHO] */
    J[790] += dqdci;              /* dwdot[CH3]/d[CH3CHO] */
    J[806] -= dqdci;              /* dwdot[CH3CHO]/d[CH3CHO] */
    /* d()/dT */
    J[873] -= dqdT;               /* dwdot[O]/dT */
    J[874] += dqdT;               /* dwdot[OH]/dT */
    J[879] += dqdT;               /* dwdot[CO]/dT */
    J[880] += dqdT;               /* dwdot[CH3]/dT */
    J[896] -= dqdT;               /* dwdot[CH3CHO]/dT */

    /*reaction 162: CH3CHO + O2 => CH3 + CO + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[26];
    k_f = prefactor_units[161] * fwd_A[161]
                * exp(fwd_beta[161] * tc[0] - activation_units[161] * fwd_Ea[161] * invT);
    dlnkfdT = fwd_beta[161] * invT + activation_units[161] * fwd_Ea[161] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[8] += q; /* HO2 */
    wdot[9] += q; /* CO */
    wdot[10] += q; /* CH3 */
    wdot[26] -= q; /* CH3CHO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[26];
    J[155] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[158] += dqdci;              /* dwdot[HO2]/d[O2] */
    J[159] += dqdci;              /* dwdot[CO]/d[O2] */
    J[160] += dqdci;              /* dwdot[CH3]/d[O2] */
    J[176] -= dqdci;              /* dwdot[CH3CHO]/d[O2] */
    /* d()/d[CH3CHO] */
    dqdci =  + k_f*sc[5];
    J[785] -= dqdci;              /* dwdot[O2]/d[CH3CHO] */
    J[788] += dqdci;              /* dwdot[HO2]/d[CH3CHO] */
    J[789] += dqdci;              /* dwdot[CO]/d[CH3CHO] */
    J[790] += dqdci;              /* dwdot[CH3]/d[CH3CHO] */
    J[806] -= dqdci;              /* dwdot[CH3CHO]/d[CH3CHO] */
    /* d()/dT */
    J[875] -= dqdT;               /* dwdot[O2]/dT */
    J[878] += dqdT;               /* dwdot[HO2]/dT */
    J[879] += dqdT;               /* dwdot[CO]/dT */
    J[880] += dqdT;               /* dwdot[CH3]/dT */
    J[896] -= dqdT;               /* dwdot[CH3CHO]/dT */

    /*reaction 163: CH3CHO + OH => CH3 + CO + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[26];
    k_f = prefactor_units[162] * fwd_A[162]
                * exp(fwd_beta[162] * tc[0] - activation_units[162] * fwd_Ea[162] * invT);
    dlnkfdT = fwd_beta[162] * invT + activation_units[162] * fwd_Ea[162] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[7] += q; /* H2O */
    wdot[9] += q; /* CO */
    wdot[10] += q; /* CH3 */
    wdot[26] -= q; /* CH3CHO */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[26];
    J[124] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[127] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[129] += dqdci;              /* dwdot[CO]/d[OH] */
    J[130] += dqdci;              /* dwdot[CH3]/d[OH] */
    J[146] -= dqdci;              /* dwdot[CH3CHO]/d[OH] */
    /* d()/d[CH3CHO] */
    dqdci =  + k_f*sc[4];
    J[784] -= dqdci;              /* dwdot[OH]/d[CH3CHO] */
    J[787] += dqdci;              /* dwdot[H2O]/d[CH3CHO] */
    J[789] += dqdci;              /* dwdot[CO]/d[CH3CHO] */
    J[790] += dqdci;              /* dwdot[CH3]/d[CH3CHO] */
    J[806] -= dqdci;              /* dwdot[CH3CHO]/d[CH3CHO] */
    /* d()/dT */
    J[874] -= dqdT;               /* dwdot[OH]/dT */
    J[877] += dqdT;               /* dwdot[H2O]/dT */
    J[879] += dqdT;               /* dwdot[CO]/dT */
    J[880] += dqdT;               /* dwdot[CH3]/dT */
    J[896] -= dqdT;               /* dwdot[CH3CHO]/dT */

    /*reaction 164: CH3CHO + H <=> CH2CHO + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[26];
    k_f = prefactor_units[163] * fwd_A[163]
                * exp(fwd_beta[163] * tc[0] - activation_units[163] * fwd_Ea[163] * invT);
    dlnkfdT = fwd_beta[163] * invT + activation_units[163] * fwd_Ea[163] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[27];
    Kc = exp(g_RT[1] - g_RT[2] + g_RT[26] - g_RT[27]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[26]) + (h_RT[2] + h_RT[27]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[2] += q; /* H2 */
    wdot[26] -= q; /* CH3CHO */
    wdot[27] += q; /* CH2CHO */
    /* d()/d[H] */
    dqdci =  + k_f*sc[26];
    J[31] -= dqdci;               /* dwdot[H]/d[H] */
    J[32] += dqdci;               /* dwdot[H2]/d[H] */
    J[56] -= dqdci;               /* dwdot[CH3CHO]/d[H] */
    J[57] += dqdci;               /* dwdot[CH2CHO]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[27];
    J[61] -= dqdci;               /* dwdot[H]/d[H2] */
    J[62] += dqdci;               /* dwdot[H2]/d[H2] */
    J[86] -= dqdci;               /* dwdot[CH3CHO]/d[H2] */
    J[87] += dqdci;               /* dwdot[CH2CHO]/d[H2] */
    /* d()/d[CH3CHO] */
    dqdci =  + k_f*sc[1];
    J[781] -= dqdci;              /* dwdot[H]/d[CH3CHO] */
    J[782] += dqdci;              /* dwdot[H2]/d[CH3CHO] */
    J[806] -= dqdci;              /* dwdot[CH3CHO]/d[CH3CHO] */
    J[807] += dqdci;              /* dwdot[CH2CHO]/d[CH3CHO] */
    /* d()/d[CH2CHO] */
    dqdci =  - k_r*sc[2];
    J[811] -= dqdci;              /* dwdot[H]/d[CH2CHO] */
    J[812] += dqdci;              /* dwdot[H2]/d[CH2CHO] */
    J[836] -= dqdci;              /* dwdot[CH3CHO]/d[CH2CHO] */
    J[837] += dqdci;              /* dwdot[CH2CHO]/d[CH2CHO] */
    /* d()/dT */
    J[871] -= dqdT;               /* dwdot[H]/dT */
    J[872] += dqdT;               /* dwdot[H2]/dT */
    J[896] -= dqdT;               /* dwdot[CH3CHO]/dT */
    J[897] += dqdT;               /* dwdot[CH2CHO]/dT */

    /*reaction 165: CH3CHO + H => CH3 + CO + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[26];
    k_f = prefactor_units[164] * fwd_A[164]
                * exp(fwd_beta[164] * tc[0] - activation_units[164] * fwd_Ea[164] * invT);
    dlnkfdT = fwd_beta[164] * invT + activation_units[164] * fwd_Ea[164] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[2] += q; /* H2 */
    wdot[9] += q; /* CO */
    wdot[10] += q; /* CH3 */
    wdot[26] -= q; /* CH3CHO */
    /* d()/d[H] */
    dqdci =  + k_f*sc[26];
    J[31] -= dqdci;               /* dwdot[H]/d[H] */
    J[32] += dqdci;               /* dwdot[H2]/d[H] */
    J[39] += dqdci;               /* dwdot[CO]/d[H] */
    J[40] += dqdci;               /* dwdot[CH3]/d[H] */
    J[56] -= dqdci;               /* dwdot[CH3CHO]/d[H] */
    /* d()/d[CH3CHO] */
    dqdci =  + k_f*sc[1];
    J[781] -= dqdci;              /* dwdot[H]/d[CH3CHO] */
    J[782] += dqdci;              /* dwdot[H2]/d[CH3CHO] */
    J[789] += dqdci;              /* dwdot[CO]/d[CH3CHO] */
    J[790] += dqdci;              /* dwdot[CH3]/d[CH3CHO] */
    J[806] -= dqdci;              /* dwdot[CH3CHO]/d[CH3CHO] */
    /* d()/dT */
    J[871] -= dqdT;               /* dwdot[H]/dT */
    J[872] += dqdT;               /* dwdot[H2]/dT */
    J[879] += dqdT;               /* dwdot[CO]/dT */
    J[880] += dqdT;               /* dwdot[CH3]/dT */
    J[896] -= dqdT;               /* dwdot[CH3CHO]/dT */

    /*reaction 166: CH3CHO + O <=> CH2CHO + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[26];
    k_f = prefactor_units[165] * fwd_A[165]
                * exp(fwd_beta[165] * tc[0] - activation_units[165] * fwd_Ea[165] * invT);
    dlnkfdT = fwd_beta[165] * invT + activation_units[165] * fwd_Ea[165] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[27];
    Kc = exp(g_RT[3] - g_RT[4] + g_RT[26] - g_RT[27]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[26]) + (h_RT[4] + h_RT[27]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* O */
    wdot[4] += q; /* OH */
    wdot[26] -= q; /* CH3CHO */
    wdot[27] += q; /* CH2CHO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[26];
    J[93] -= dqdci;               /* dwdot[O]/d[O] */
    J[94] += dqdci;               /* dwdot[OH]/d[O] */
    J[116] -= dqdci;              /* dwdot[CH3CHO]/d[O] */
    J[117] += dqdci;              /* dwdot[CH2CHO]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[27];
    J[123] -= dqdci;              /* dwdot[O]/d[OH] */
    J[124] += dqdci;              /* dwdot[OH]/d[OH] */
    J[146] -= dqdci;              /* dwdot[CH3CHO]/d[OH] */
    J[147] += dqdci;              /* dwdot[CH2CHO]/d[OH] */
    /* d()/d[CH3CHO] */
    dqdci =  + k_f*sc[3];
    J[783] -= dqdci;              /* dwdot[O]/d[CH3CHO] */
    J[784] += dqdci;              /* dwdot[OH]/d[CH3CHO] */
    J[806] -= dqdci;              /* dwdot[CH3CHO]/d[CH3CHO] */
    J[807] += dqdci;              /* dwdot[CH2CHO]/d[CH3CHO] */
    /* d()/d[CH2CHO] */
    dqdci =  - k_r*sc[4];
    J[813] -= dqdci;              /* dwdot[O]/d[CH2CHO] */
    J[814] += dqdci;              /* dwdot[OH]/d[CH2CHO] */
    J[836] -= dqdci;              /* dwdot[CH3CHO]/d[CH2CHO] */
    J[837] += dqdci;              /* dwdot[CH2CHO]/d[CH2CHO] */
    /* d()/dT */
    J[873] -= dqdT;               /* dwdot[O]/dT */
    J[874] += dqdT;               /* dwdot[OH]/dT */
    J[896] -= dqdT;               /* dwdot[CH3CHO]/dT */
    J[897] += dqdT;               /* dwdot[CH2CHO]/dT */

    /*reaction 167: CH3CHO + CH3 => CH3 + CO + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[10]*sc[26];
    k_f = prefactor_units[166] * fwd_A[166]
                * exp(fwd_beta[166] * tc[0] - activation_units[166] * fwd_Ea[166] * invT);
    dlnkfdT = fwd_beta[166] * invT + activation_units[166] * fwd_Ea[166] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[9] += q; /* CO */
    wdot[13] += q; /* CH4 */
    wdot[26] -= q; /* CH3CHO */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[26];
    J[309] += dqdci;              /* dwdot[CO]/d[CH3] */
    J[313] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[326] -= dqdci;              /* dwdot[CH3CHO]/d[CH3] */
    /* d()/d[CH3CHO] */
    dqdci =  + k_f*sc[10];
    J[789] += dqdci;              /* dwdot[CO]/d[CH3CHO] */
    J[793] += dqdci;              /* dwdot[CH4]/d[CH3CHO] */
    J[806] -= dqdci;              /* dwdot[CH3CHO]/d[CH3CHO] */
    /* d()/dT */
    J[879] += dqdT;               /* dwdot[CO]/dT */
    J[883] += dqdT;               /* dwdot[CH4]/dT */
    J[896] -= dqdT;               /* dwdot[CH3CHO]/dT */

    /*reaction 168: CH3CHO + HO2 => CH3 + CO + H2O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[26];
    k_f = prefactor_units[167] * fwd_A[167]
                * exp(fwd_beta[167] * tc[0] - activation_units[167] * fwd_Ea[167] * invT);
    dlnkfdT = fwd_beta[167] * invT + activation_units[167] * fwd_Ea[167] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[6] += q; /* H2O2 */
    wdot[8] -= q; /* HO2 */
    wdot[9] += q; /* CO */
    wdot[10] += q; /* CH3 */
    wdot[26] -= q; /* CH3CHO */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[26];
    J[246] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    J[248] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[249] += dqdci;              /* dwdot[CO]/d[HO2] */
    J[250] += dqdci;              /* dwdot[CH3]/d[HO2] */
    J[266] -= dqdci;              /* dwdot[CH3CHO]/d[HO2] */
    /* d()/d[CH3CHO] */
    dqdci =  + k_f*sc[8];
    J[786] += dqdci;              /* dwdot[H2O2]/d[CH3CHO] */
    J[788] -= dqdci;              /* dwdot[HO2]/d[CH3CHO] */
    J[789] += dqdci;              /* dwdot[CO]/d[CH3CHO] */
    J[790] += dqdci;              /* dwdot[CH3]/d[CH3CHO] */
    J[806] -= dqdci;              /* dwdot[CH3CHO]/d[CH3CHO] */
    /* d()/dT */
    J[876] += dqdT;               /* dwdot[H2O2]/dT */
    J[878] -= dqdT;               /* dwdot[HO2]/dT */
    J[879] += dqdT;               /* dwdot[CO]/dT */
    J[880] += dqdT;               /* dwdot[CH3]/dT */
    J[896] -= dqdT;               /* dwdot[CH3CHO]/dT */

    /*reaction 169: C2H5O <=> CH3 + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[28];
    k_f = prefactor_units[168] * fwd_A[168]
                * exp(fwd_beta[168] * tc[0] - activation_units[168] * fwd_Ea[168] * invT);
    dlnkfdT = fwd_beta[168] * invT + activation_units[168] * fwd_Ea[168] * invT2;
    /* reverse */
    phi_r = sc[10]*sc[11];
    Kc = refC * exp(-g_RT[10] - g_RT[11] + g_RT[28]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[28]) + (h_RT[10] + h_RT[11]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[10] += q; /* CH3 */
    wdot[11] += q; /* CH2O */
    wdot[28] -= q; /* C2H5O */
    /* d()/d[CH3] */
    dqdci =  - k_r*sc[11];
    J[310] += dqdci;              /* dwdot[CH3]/d[CH3] */
    J[311] += dqdci;              /* dwdot[CH2O]/d[CH3] */
    J[328] -= dqdci;              /* dwdot[C2H5O]/d[CH3] */
    /* d()/d[CH2O] */
    dqdci =  - k_r*sc[10];
    J[340] += dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[341] += dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[358] -= dqdci;              /* dwdot[C2H5O]/d[CH2O] */
    /* d()/d[C2H5O] */
    dqdci =  + k_f;
    J[850] += dqdci;              /* dwdot[CH3]/d[C2H5O] */
    J[851] += dqdci;              /* dwdot[CH2O]/d[C2H5O] */
    J[868] -= dqdci;              /* dwdot[C2H5O]/d[C2H5O] */
    /* d()/dT */
    J[880] += dqdT;               /* dwdot[CH3]/dT */
    J[881] += dqdT;               /* dwdot[CH2O]/dT */
    J[898] -= dqdT;               /* dwdot[C2H5O]/dT */

    /*reaction 170: C2H5O <=> CH3CHO + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[28];
    k_f = prefactor_units[169] * fwd_A[169]
                * exp(fwd_beta[169] * tc[0] - activation_units[169] * fwd_Ea[169] * invT);
    dlnkfdT = fwd_beta[169] * invT + activation_units[169] * fwd_Ea[169] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[26];
    Kc = refC * exp(-g_RT[1] - g_RT[26] + g_RT[28]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[28]) + (h_RT[1] + h_RT[26]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] += q; /* H */
    wdot[26] += q; /* CH3CHO */
    wdot[28] -= q; /* C2H5O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[26];
    J[31] += dqdci;               /* dwdot[H]/d[H] */
    J[56] += dqdci;               /* dwdot[CH3CHO]/d[H] */
    J[58] -= dqdci;               /* dwdot[C2H5O]/d[H] */
    /* d()/d[CH3CHO] */
    dqdci =  - k_r*sc[1];
    J[781] += dqdci;              /* dwdot[H]/d[CH3CHO] */
    J[806] += dqdci;              /* dwdot[CH3CHO]/d[CH3CHO] */
    J[808] -= dqdci;              /* dwdot[C2H5O]/d[CH3CHO] */
    /* d()/d[C2H5O] */
    dqdci =  + k_f;
    J[841] += dqdci;              /* dwdot[H]/d[C2H5O] */
    J[866] += dqdci;              /* dwdot[CH3CHO]/d[C2H5O] */
    J[868] -= dqdci;              /* dwdot[C2H5O]/d[C2H5O] */
    /* d()/dT */
    J[871] += dqdT;               /* dwdot[H]/dT */
    J[896] += dqdT;               /* dwdot[CH3CHO]/dT */
    J[898] -= dqdT;               /* dwdot[C2H5O]/dT */

    /*reaction 171: C2H5O + O2 <=> CH3CHO + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[28];
    k_f = prefactor_units[170] * fwd_A[170]
                * exp(fwd_beta[170] * tc[0] - activation_units[170] * fwd_Ea[170] * invT);
    dlnkfdT = fwd_beta[170] * invT + activation_units[170] * fwd_Ea[170] * invT2;
    /* reverse */
    phi_r = sc[8]*sc[26];
    Kc = exp(g_RT[5] - g_RT[8] - g_RT[26] + g_RT[28]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[28]) + (h_RT[8] + h_RT[26]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[8] += q; /* HO2 */
    wdot[26] += q; /* CH3CHO */
    wdot[28] -= q; /* C2H5O */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[28];
    J[155] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[158] += dqdci;              /* dwdot[HO2]/d[O2] */
    J[176] += dqdci;              /* dwdot[CH3CHO]/d[O2] */
    J[178] -= dqdci;              /* dwdot[C2H5O]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[26];
    J[245] -= dqdci;              /* dwdot[O2]/d[HO2] */
    J[248] += dqdci;              /* dwdot[HO2]/d[HO2] */
    J[266] += dqdci;              /* dwdot[CH3CHO]/d[HO2] */
    J[268] -= dqdci;              /* dwdot[C2H5O]/d[HO2] */
    /* d()/d[CH3CHO] */
    dqdci =  - k_r*sc[8];
    J[785] -= dqdci;              /* dwdot[O2]/d[CH3CHO] */
    J[788] += dqdci;              /* dwdot[HO2]/d[CH3CHO] */
    J[806] += dqdci;              /* dwdot[CH3CHO]/d[CH3CHO] */
    J[808] -= dqdci;              /* dwdot[C2H5O]/d[CH3CHO] */
    /* d()/d[C2H5O] */
    dqdci =  + k_f*sc[5];
    J[845] -= dqdci;              /* dwdot[O2]/d[C2H5O] */
    J[848] += dqdci;              /* dwdot[HO2]/d[C2H5O] */
    J[866] += dqdci;              /* dwdot[CH3CHO]/d[C2H5O] */
    J[868] -= dqdci;              /* dwdot[C2H5O]/d[C2H5O] */
    /* d()/dT */
    J[875] -= dqdT;               /* dwdot[O2]/dT */
    J[878] += dqdT;               /* dwdot[HO2]/dT */
    J[896] += dqdT;               /* dwdot[CH3CHO]/dT */
    J[898] -= dqdT;               /* dwdot[C2H5O]/dT */

    /*reaction 172: SXCH2 + N2 <=> TXCH2 + N2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[22];
    k_f = prefactor_units[171] * fwd_A[171]
                * exp(fwd_beta[171] * tc[0] - activation_units[171] * fwd_Ea[171] * invT);
    dlnkfdT = fwd_beta[171] * invT + activation_units[171] * fwd_Ea[171] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[21];
    Kc = exp(g_RT[0] - g_RT[0] - g_RT[21] + g_RT[22]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[22]) + (h_RT[0] + h_RT[21]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[21] += q; /* TXCH2 */
    wdot[22] -= q; /* SXCH2 */
    /* d()/d[N2] */
    dqdci =  + k_f*sc[22] - k_r*sc[21];
    J[21] += dqdci;               /* dwdot[TXCH2]/d[N2] */
    J[22] -= dqdci;               /* dwdot[SXCH2]/d[N2] */
    /* d()/d[TXCH2] */
    dqdci =  - k_r*sc[0];
    J[651] += dqdci;              /* dwdot[TXCH2]/d[TXCH2] */
    J[652] -= dqdci;              /* dwdot[SXCH2]/d[TXCH2] */
    /* d()/d[SXCH2] */
    dqdci =  + k_f*sc[0];
    J[681] += dqdci;              /* dwdot[TXCH2]/d[SXCH2] */
    J[682] -= dqdci;              /* dwdot[SXCH2]/d[SXCH2] */
    /* d()/dT */
    J[891] += dqdT;               /* dwdot[TXCH2]/dT */
    J[892] -= dqdT;               /* dwdot[SXCH2]/dT */

    amrex::Real c_R[29], dcRdT[29], e_RT[29];
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
    for (int k = 0; k < 29; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[870+k];
    }

    amrex::Real cmixinv = 1.0/cmix;
    amrex::Real tmp1 = ehmix*cmixinv;
    amrex::Real tmp3 = cmixinv*T;
    amrex::Real tmp2 = tmp1*tmp3;
    amrex::Real dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 29; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 29; ++m) {
            dehmixdc += eh_RT[m]*J[k*30+m];
        }
        J[k*30+29] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[899] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;
}

#endif
/* Transport function declarations  */


void egtransetLENIMC(int* LENIMC ) {
    *LENIMC = 118;}


void egtransetLENRMC(int* LENRMC ) {
    *LENRMC = 16994;}


void egtransetNO(int* NO ) {
    *NO = 4;}


void egtransetKK(int* KK ) {
    *KK = 29;}


void egtransetNLITE(int* NLITE ) {
    *NLITE = 2;}


/*Patm in ergs/cm3 */
void egtransetPATM(amrex::Real* PATM) {
    *PATM =   0.1013250000000000E+07;}


/*the molecular weights in g/mol */
void egtransetWT(amrex::Real* WT ) {
    WT[0] = 2.80134000E+01;
    WT[1] = 1.00797000E+00;
    WT[2] = 2.01594000E+00;
    WT[3] = 1.59994000E+01;
    WT[4] = 1.70073700E+01;
    WT[5] = 3.19988000E+01;
    WT[6] = 3.40147400E+01;
    WT[7] = 1.80153400E+01;
    WT[8] = 3.30067700E+01;
    WT[9] = 2.80105500E+01;
    WT[10] = 1.50350600E+01;
    WT[11] = 3.00264900E+01;
    WT[12] = 4.40099500E+01;
    WT[13] = 1.60430300E+01;
    WT[14] = 2.60382400E+01;
    WT[15] = 2.80541800E+01;
    WT[16] = 4.20376400E+01;
    WT[17] = 3.00701200E+01;
    WT[18] = 1.20111500E+01;
    WT[19] = 1.30191200E+01;
    WT[20] = 2.90185200E+01;
    WT[21] = 1.40270900E+01;
    WT[22] = 1.40270900E+01;
    WT[23] = 2.70462100E+01;
    WT[24] = 2.90621500E+01;
    WT[25] = 4.10296700E+01;
    WT[26] = 4.40535800E+01;
    WT[27] = 4.30456100E+01;
    WT[28] = 4.50615500E+01;
}


/*the lennard-jones potential well depth eps/kb in K */
void egtransetEPS(amrex::Real* EPS ) {
    EPS[0] = 9.75300000E+01;
    EPS[1] = 1.45000000E+02;
    EPS[2] = 3.80000000E+01;
    EPS[3] = 8.00000000E+01;
    EPS[4] = 8.00000000E+01;
    EPS[5] = 1.07400000E+02;
    EPS[6] = 1.07400000E+02;
    EPS[7] = 5.72400000E+02;
    EPS[8] = 1.07400000E+02;
    EPS[9] = 9.81000000E+01;
    EPS[10] = 1.44000000E+02;
    EPS[11] = 4.98000000E+02;
    EPS[12] = 2.44000000E+02;
    EPS[13] = 1.41400000E+02;
    EPS[14] = 2.09000000E+02;
    EPS[15] = 2.80800000E+02;
    EPS[16] = 4.36000000E+02;
    EPS[17] = 2.52300000E+02;
    EPS[18] = 7.14000000E+01;
    EPS[19] = 8.00000000E+01;
    EPS[20] = 4.98000000E+02;
    EPS[21] = 1.44000000E+02;
    EPS[22] = 1.44000000E+02;
    EPS[23] = 2.09000000E+02;
    EPS[24] = 2.52300000E+02;
    EPS[25] = 1.50000000E+02;
    EPS[26] = 4.36000000E+02;
    EPS[27] = 4.36000000E+02;
    EPS[28] = 4.70600000E+02;
}


/*the lennard-jones collision diameter in Angstroms */
void egtransetSIG(amrex::Real* SIG ) {
    SIG[0] = 3.62100000E+00;
    SIG[1] = 2.05000000E+00;
    SIG[2] = 2.92000000E+00;
    SIG[3] = 2.75000000E+00;
    SIG[4] = 2.75000000E+00;
    SIG[5] = 3.45800000E+00;
    SIG[6] = 3.45800000E+00;
    SIG[7] = 2.60500000E+00;
    SIG[8] = 3.45800000E+00;
    SIG[9] = 3.65000000E+00;
    SIG[10] = 3.80000000E+00;
    SIG[11] = 3.59000000E+00;
    SIG[12] = 3.76300000E+00;
    SIG[13] = 3.74600000E+00;
    SIG[14] = 4.10000000E+00;
    SIG[15] = 3.97100000E+00;
    SIG[16] = 3.97000000E+00;
    SIG[17] = 4.30200000E+00;
    SIG[18] = 3.29800000E+00;
    SIG[19] = 2.75000000E+00;
    SIG[20] = 3.59000000E+00;
    SIG[21] = 3.80000000E+00;
    SIG[22] = 3.80000000E+00;
    SIG[23] = 4.10000000E+00;
    SIG[24] = 4.30200000E+00;
    SIG[25] = 2.50000000E+00;
    SIG[26] = 3.97000000E+00;
    SIG[27] = 3.97000000E+00;
    SIG[28] = 4.41000000E+00;
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
    DIP[7] = 1.84400000E+00;
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
    DIP[18] = 0.00000000E+00;
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
}


/*the polarizability in cubic Angstroms */
void egtransetPOL(amrex::Real* POL ) {
    POL[0] = 1.76000000E+00;
    POL[1] = 0.00000000E+00;
    POL[2] = 7.90000000E-01;
    POL[3] = 0.00000000E+00;
    POL[4] = 0.00000000E+00;
    POL[5] = 1.60000000E+00;
    POL[6] = 0.00000000E+00;
    POL[7] = 0.00000000E+00;
    POL[8] = 0.00000000E+00;
    POL[9] = 1.95000000E+00;
    POL[10] = 0.00000000E+00;
    POL[11] = 0.00000000E+00;
    POL[12] = 2.65000000E+00;
    POL[13] = 2.60000000E+00;
    POL[14] = 0.00000000E+00;
    POL[15] = 0.00000000E+00;
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
}


/*the rotational relaxation collision number at 298 K */
void egtransetZROT(amrex::Real* ZROT ) {
    ZROT[0] = 4.00000000E+00;
    ZROT[1] = 0.00000000E+00;
    ZROT[2] = 2.80000000E+02;
    ZROT[3] = 0.00000000E+00;
    ZROT[4] = 0.00000000E+00;
    ZROT[5] = 3.80000000E+00;
    ZROT[6] = 3.80000000E+00;
    ZROT[7] = 4.00000000E+00;
    ZROT[8] = 1.00000000E+00;
    ZROT[9] = 1.80000000E+00;
    ZROT[10] = 0.00000000E+00;
    ZROT[11] = 2.00000000E+00;
    ZROT[12] = 2.10000000E+00;
    ZROT[13] = 1.30000000E+01;
    ZROT[14] = 2.50000000E+00;
    ZROT[15] = 1.50000000E+00;
    ZROT[16] = 2.00000000E+00;
    ZROT[17] = 1.50000000E+00;
    ZROT[18] = 0.00000000E+00;
    ZROT[19] = 0.00000000E+00;
    ZROT[20] = 0.00000000E+00;
    ZROT[21] = 0.00000000E+00;
    ZROT[22] = 0.00000000E+00;
    ZROT[23] = 1.00000000E+00;
    ZROT[24] = 1.50000000E+00;
    ZROT[25] = 1.00000000E+00;
    ZROT[26] = 2.00000000E+00;
    ZROT[27] = 2.00000000E+00;
    ZROT[28] = 1.50000000E+00;
}


/*0: monoatomic, 1: linear, 2: nonlinear */
void egtransetNLIN(int* NLIN) {
    NLIN[0] = 1;
    NLIN[1] = 0;
    NLIN[2] = 1;
    NLIN[3] = 0;
    NLIN[4] = 1;
    NLIN[5] = 1;
    NLIN[6] = 2;
    NLIN[7] = 2;
    NLIN[8] = 2;
    NLIN[9] = 1;
    NLIN[10] = 1;
    NLIN[11] = 2;
    NLIN[12] = 1;
    NLIN[13] = 2;
    NLIN[14] = 1;
    NLIN[15] = 2;
    NLIN[16] = 2;
    NLIN[17] = 2;
    NLIN[18] = 0;
    NLIN[19] = 1;
    NLIN[20] = 2;
    NLIN[21] = 1;
    NLIN[22] = 1;
    NLIN[23] = 2;
    NLIN[24] = 2;
    NLIN[25] = 2;
    NLIN[26] = 2;
    NLIN[27] = 2;
    NLIN[28] = 2;
}


/*Poly fits for the viscosities, dim NO*KK */
void egtransetCOFETA(amrex::Real* COFETA) {
    COFETA[0] = -1.65695594E+01;
    COFETA[1] = 2.39056562E+00;
    COFETA[2] = -2.34558144E-01;
    COFETA[3] = 1.05024037E-02;
    COFETA[4] = -2.04078397E+01;
    COFETA[5] = 3.65436395E+00;
    COFETA[6] = -3.98339635E-01;
    COFETA[7] = 1.75883009E-02;
    COFETA[8] = -1.38347699E+01;
    COFETA[9] = 1.00106621E+00;
    COFETA[10] = -4.98105694E-02;
    COFETA[11] = 2.31450475E-03;
    COFETA[12] = -1.50926240E+01;
    COFETA[13] = 1.92606504E+00;
    COFETA[14] = -1.73487476E-01;
    COFETA[15] = 7.82572931E-03;
    COFETA[16] = -1.50620763E+01;
    COFETA[17] = 1.92606504E+00;
    COFETA[18] = -1.73487476E-01;
    COFETA[19] = 7.82572931E-03;
    COFETA[20] = -1.71618309E+01;
    COFETA[21] = 2.68036374E+00;
    COFETA[22] = -2.72570227E-01;
    COFETA[23] = 1.21650964E-02;
    COFETA[24] = -1.71312832E+01;
    COFETA[25] = 2.68036374E+00;
    COFETA[26] = -2.72570227E-01;
    COFETA[27] = 1.21650964E-02;
    COFETA[28] = -1.05420863E+01;
    COFETA[29] = -1.37777096E+00;
    COFETA[30] = 4.20502308E-01;
    COFETA[31] = -2.40627230E-02;
    COFETA[32] = -1.71463238E+01;
    COFETA[33] = 2.68036374E+00;
    COFETA[34] = -2.72570227E-01;
    COFETA[35] = 1.21650964E-02;
    COFETA[36] = -1.66188336E+01;
    COFETA[37] = 2.40307799E+00;
    COFETA[38] = -2.36167638E-01;
    COFETA[39] = 1.05714061E-02;
    COFETA[40] = -2.02316497E+01;
    COFETA[41] = 3.63241793E+00;
    COFETA[42] = -3.95581049E-01;
    COFETA[43] = 1.74725495E-02;
    COFETA[44] = -1.98330577E+01;
    COFETA[45] = 2.69480162E+00;
    COFETA[46] = -1.65880845E-01;
    COFETA[47] = 3.14504769E-03;
    COFETA[48] = -2.40014975E+01;
    COFETA[49] = 5.14359547E+00;
    COFETA[50] = -5.74269731E-01;
    COFETA[51] = 2.44937679E-02;
    COFETA[52] = -2.00094664E+01;
    COFETA[53] = 3.57220167E+00;
    COFETA[54] = -3.87936446E-01;
    COFETA[55] = 1.71483254E-02;
    COFETA[56] = -2.33666446E+01;
    COFETA[57] = 4.80350223E+00;
    COFETA[58] = -5.38341336E-01;
    COFETA[59] = 2.32747213E-02;
    COFETA[60] = -2.50655444E+01;
    COFETA[61] = 5.33982977E+00;
    COFETA[62] = -5.89982992E-01;
    COFETA[63] = 2.47780650E-02;
    COFETA[64] = -2.23395647E+01;
    COFETA[65] = 3.86433912E+00;
    COFETA[66] = -3.41553983E-01;
    COFETA[67] = 1.17083447E-02;
    COFETA[68] = -2.46410937E+01;
    COFETA[69] = 5.19497183E+00;
    COFETA[70] = -5.78827339E-01;
    COFETA[71] = 2.46050921E-02;
    COFETA[72] = -1.50067648E+01;
    COFETA[73] = 1.69625877E+00;
    COFETA[74] = -1.42936462E-01;
    COFETA[75] = 6.47223426E-03;
    COFETA[76] = -1.51956901E+01;
    COFETA[77] = 1.92606504E+00;
    COFETA[78] = -1.73487476E-01;
    COFETA[79] = 7.82572931E-03;
    COFETA[80] = -1.98501306E+01;
    COFETA[81] = 2.69480162E+00;
    COFETA[82] = -1.65880845E-01;
    COFETA[83] = 3.14504769E-03;
    COFETA[84] = -2.02663469E+01;
    COFETA[85] = 3.63241793E+00;
    COFETA[86] = -3.95581049E-01;
    COFETA[87] = 1.74725495E-02;
    COFETA[88] = -2.02663469E+01;
    COFETA[89] = 3.63241793E+00;
    COFETA[90] = -3.95581049E-01;
    COFETA[91] = 1.74725495E-02;
    COFETA[92] = -2.33476543E+01;
    COFETA[93] = 4.80350223E+00;
    COFETA[94] = -5.38341336E-01;
    COFETA[95] = 2.32747213E-02;
    COFETA[96] = -2.46581414E+01;
    COFETA[97] = 5.19497183E+00;
    COFETA[98] = -5.78827339E-01;
    COFETA[99] = 2.46050921E-02;
    COFETA[100] = -1.92183831E+01;
    COFETA[101] = 3.75164499E+00;
    COFETA[102] = -4.10390993E-01;
    COFETA[103] = 1.80861665E-02;
    COFETA[104] = -2.23161441E+01;
    COFETA[105] = 3.86433912E+00;
    COFETA[106] = -3.41553983E-01;
    COFETA[107] = 1.17083447E-02;
    COFETA[108] = -2.23277173E+01;
    COFETA[109] = 3.86433912E+00;
    COFETA[110] = -3.41553983E-01;
    COFETA[111] = 1.17083447E-02;
    COFETA[112] = -2.12500723E+01;
    COFETA[113] = 3.25726898E+00;
    COFETA[114] = -2.49519605E-01;
    COFETA[115] = 7.19215196E-03;
}


/*Poly fits for the conductivities, dim NO*KK */
void egtransetCOFLAM(amrex::Real* COFLAM) {
    COFLAM[0] = 1.29306158E+01;
    COFLAM[1] = -3.52817362E+00;
    COFLAM[2] = 6.45499013E-01;
    COFLAM[3] = -3.19375299E-02;
    COFLAM[4] = -8.57929284E-01;
    COFLAM[5] = 3.65436395E+00;
    COFLAM[6] = -3.98339635E-01;
    COFLAM[7] = 1.75883009E-02;
    COFLAM[8] = 9.24084392E+00;
    COFLAM[9] = -4.69568931E-01;
    COFLAM[10] = 1.15980279E-01;
    COFLAM[11] = -2.61033830E-03;
    COFLAM[12] = 1.69267361E+00;
    COFLAM[13] = 1.92606504E+00;
    COFLAM[14] = -1.73487476E-01;
    COFLAM[15] = 7.82572931E-03;
    COFLAM[16] = 1.50119731E+01;
    COFLAM[17] = -3.63267854E+00;
    COFLAM[18] = 5.92839101E-01;
    COFLAM[19] = -2.62920439E-02;
    COFLAM[20] = -1.93718739E+00;
    COFLAM[21] = 2.89110219E+00;
    COFLAM[22] = -2.71096923E-01;
    COFLAM[23] = 1.15344907E-02;
    COFLAM[24] = 8.83996545E-01;
    COFLAM[25] = 1.31525428E+00;
    COFLAM[26] = 1.91774576E-02;
    COFLAM[27] = -4.41642722E-03;
    COFLAM[28] = 2.33729817E+01;
    COFLAM[29] = -8.96536433E+00;
    COFLAM[30] = 1.52828665E+00;
    COFLAM[31] = -7.58551778E-02;
    COFLAM[32] = -1.12960913E+00;
    COFLAM[33] = 2.34014305E+00;
    COFLAM[34] = -1.63245030E-01;
    COFLAM[35] = 5.80319600E-03;
    COFLAM[36] = 1.18777264E+01;
    COFLAM[37] = -3.15463949E+00;
    COFLAM[38] = 6.01973268E-01;
    COFLAM[39] = -3.03211261E-02;
    COFLAM[40] = 1.32973935E+01;
    COFLAM[41] = -4.31770522E+00;
    COFLAM[42] = 8.57562813E-01;
    COFLAM[43] = -4.51644164E-02;
    COFLAM[44] = 5.39305086E+00;
    COFLAM[45] = -2.39312375E+00;
    COFLAM[46] = 7.39585221E-01;
    COFLAM[47] = -4.58435589E-02;
    COFLAM[48] = -1.13649314E+01;
    COFLAM[49] = 5.88177395E+00;
    COFLAM[50] = -5.68651819E-01;
    COFLAM[51] = 2.03561485E-02;
    COFLAM[52] = 1.07413889E+01;
    COFLAM[53] = -3.75495104E+00;
    COFLAM[54] = 8.44678023E-01;
    COFLAM[55] = -4.65824234E-02;
    COFLAM[56] = -7.70164502E+00;
    COFLAM[57] = 4.56884453E+00;
    COFLAM[58] = -4.04747583E-01;
    COFLAM[59] = 1.40841060E-02;
    COFLAM[60] = -1.46152839E+01;
    COFLAM[61] = 6.36251406E+00;
    COFLAM[62] = -5.03832130E-01;
    COFLAM[63] = 1.26121050E-02;
    COFLAM[64] = -8.32871231E+00;
    COFLAM[65] = 3.97067262E+00;
    COFLAM[66] = -2.21252287E-01;
    COFLAM[67] = 1.47870329E-03;
    COFLAM[68] = -1.09902209E+01;
    COFLAM[69] = 4.70647707E+00;
    COFLAM[70] = -2.52272495E-01;
    COFLAM[71] = 1.75193258E-04;
    COFLAM[72] = 2.06524872E+00;
    COFLAM[73] = 1.69625877E+00;
    COFLAM[74] = -1.42936462E-01;
    COFLAM[75] = 6.47223426E-03;
    COFLAM[76] = 2.08093824E+01;
    COFLAM[77] = -6.24180163E+00;
    COFLAM[78] = 9.82387427E-01;
    COFLAM[79] = -4.50360664E-02;
    COFLAM[80] = 6.30243508E+00;
    COFLAM[81] = -2.22810801E+00;
    COFLAM[82] = 6.37340514E-01;
    COFLAM[83] = -3.81056018E-02;
    COFLAM[84] = 1.29177902E+01;
    COFLAM[85] = -3.73745535E+00;
    COFLAM[86] = 7.15831021E-01;
    COFLAM[87] = -3.63846910E-02;
    COFLAM[88] = 1.89383266E+01;
    COFLAM[89] = -6.51018128E+00;
    COFLAM[90] = 1.13292060E+00;
    COFLAM[91] = -5.69603226E-02;
    COFLAM[92] = -9.10384516E+00;
    COFLAM[93] = 4.54798869E+00;
    COFLAM[94] = -3.18114458E-01;
    COFLAM[95] = 6.59577576E-03;
    COFLAM[96] = -8.95009705E+00;
    COFLAM[97] = 4.02515080E+00;
    COFLAM[98] = -1.84063946E-01;
    COFLAM[99] = -1.94054752E-03;
    COFLAM[100] = -5.96382813E+00;
    COFLAM[101] = 4.39356020E+00;
    COFLAM[102] = -3.97598260E-01;
    COFLAM[103] = 1.39749624E-02;
    COFLAM[104] = -8.81913142E+00;
    COFLAM[105] = 3.70554072E+00;
    COFLAM[106] = -1.14942694E-01;
    COFLAM[107] = -6.09987938E-03;
    COFLAM[108] = -1.51505495E+01;
    COFLAM[109] = 6.60055580E+00;
    COFLAM[110] = -5.50266319E-01;
    COFLAM[111] = 1.52072675E-02;
    COFLAM[112] = -3.15963092E+01;
    COFLAM[113] = 1.36690300E+01;
    COFLAM[114] = -1.57444544E+00;
    COFLAM[115] = 6.51266833E-02;
}


/*Poly fits for the diffusion coefficients, dim NO*KK*KK */
void egtransetCOFD(amrex::Real* COFD) {
    COFD[0] = -1.49828430E+01;
    COFD[1] = 3.25781069E+00;
    COFD[2] = -2.12199367E-01;
    COFD[3] = 9.36657283E-03;
    COFD[4] = -1.42894441E+01;
    COFD[5] = 3.67490723E+00;
    COFD[6] = -2.65114792E-01;
    COFD[7] = 1.16092671E-02;
    COFD[8] = -1.16906297E+01;
    COFD[9] = 2.47469981E+00;
    COFD[10] = -1.10436257E-01;
    COFD[11] = 4.95273813E-03;
    COFD[12] = -1.40756935E+01;
    COFD[13] = 3.07549274E+00;
    COFD[14] = -1.88889344E-01;
    COFD[15] = 8.37152866E-03;
    COFD[16] = -1.40949196E+01;
    COFD[17] = 3.07549274E+00;
    COFD[18] = -1.88889344E-01;
    COFD[19] = 8.37152866E-03;
    COFD[20] = -1.52414485E+01;
    COFD[21] = 3.35922578E+00;
    COFD[22] = -2.25181399E-01;
    COFD[23] = 9.92132878E-03;
    COFD[24] = -1.52554761E+01;
    COFD[25] = 3.35922578E+00;
    COFD[26] = -2.25181399E-01;
    COFD[27] = 9.92132878E-03;
    COFD[28] = -2.10643259E+01;
    COFD[29] = 5.53614847E+00;
    COFD[30] = -4.86046736E-01;
    COFD[31] = 2.03659188E-02;
    COFD[32] = -1.52486273E+01;
    COFD[33] = 3.35922578E+00;
    COFD[34] = -2.25181399E-01;
    COFD[35] = 9.92132878E-03;
    COFD[36] = -1.50031687E+01;
    COFD[37] = 3.26223357E+00;
    COFD[38] = -2.12746642E-01;
    COFD[39] = 9.38912883E-03;
    COFD[40] = -1.59633387E+01;
    COFD[41] = 3.66853818E+00;
    COFD[42] = -2.64346221E-01;
    COFD[43] = 1.15784613E-02;
    COFD[44] = -2.04833713E+01;
    COFD[45] = 5.23112374E+00;
    COFD[46] = -4.54967682E-01;
    COFD[47] = 1.93570423E-02;
    COFD[48] = -1.81432461E+01;
    COFD[49] = 4.37565431E+00;
    COFD[50] = -3.53906025E-01;
    COFD[51] = 1.53760786E-02;
    COFD[52] = -1.59327297E+01;
    COFD[53] = 3.65620899E+00;
    COFD[54] = -2.62933804E-01;
    COFD[55] = 1.15253223E-02;
    COFD[56] = -1.76002031E+01;
    COFD[57] = 4.19171952E+00;
    COFD[58] = -3.31354810E-01;
    COFD[59] = 1.44520623E-02;
    COFD[60] = -1.85864144E+01;
    COFD[61] = 4.54915847E+00;
    COFD[62] = -3.75000738E-01;
    COFD[63] = 1.62324821E-02;
    COFD[64] = -2.02646611E+01;
    COFD[65] = 5.10426133E+00;
    COFD[66] = -4.41256919E-01;
    COFD[67] = 1.88737290E-02;
    COFD[68] = -1.83249299E+01;
    COFD[69] = 4.42045763E+00;
    COFD[70] = -3.59451578E-01;
    COFD[71] = 1.56056164E-02;
    COFD[72] = -1.38661480E+01;
    COFD[73] = 2.97137588E+00;
    COFD[74] = -1.75491257E-01;
    COFD[75] = 7.79646773E-03;
    COFD[76] = -1.40076852E+01;
    COFD[77] = 3.07549274E+00;
    COFD[78] = -1.88889344E-01;
    COFD[79] = 8.37152866E-03;
    COFD[80] = -2.04750581E+01;
    COFD[81] = 5.23112374E+00;
    COFD[82] = -4.54967682E-01;
    COFD[83] = 1.93570423E-02;
    COFD[84] = -1.59404882E+01;
    COFD[85] = 3.66853818E+00;
    COFD[86] = -2.64346221E-01;
    COFD[87] = 1.15784613E-02;
    COFD[88] = -1.59404882E+01;
    COFD[89] = 3.66853818E+00;
    COFD[90] = -2.64346221E-01;
    COFD[91] = 1.15784613E-02;
    COFD[92] = -1.76099552E+01;
    COFD[93] = 4.19171952E+00;
    COFD[94] = -3.31354810E-01;
    COFD[95] = 1.44520623E-02;
    COFD[96] = -1.83166353E+01;
    COFD[97] = 4.42045763E+00;
    COFD[98] = -3.59451578E-01;
    COFD[99] = 1.56056164E-02;
    COFD[100] = -1.59884305E+01;
    COFD[101] = 3.72220402E+00;
    COFD[102] = -2.71150591E-01;
    COFD[103] = 1.18665265E-02;
    COFD[104] = -2.02738958E+01;
    COFD[105] = 5.10426133E+00;
    COFD[106] = -4.41256919E-01;
    COFD[107] = 1.88737290E-02;
    COFD[108] = -2.02693653E+01;
    COFD[109] = 5.10426133E+00;
    COFD[110] = -4.41256919E-01;
    COFD[111] = 1.88737290E-02;
    COFD[112] = -2.05802296E+01;
    COFD[113] = 5.16117916E+00;
    COFD[114] = -4.46897404E-01;
    COFD[115] = 1.90470443E-02;
    COFD[116] = -1.42894441E+01;
    COFD[117] = 3.67490723E+00;
    COFD[118] = -2.65114792E-01;
    COFD[119] = 1.16092671E-02;
    COFD[120] = -1.47968712E+01;
    COFD[121] = 4.23027636E+00;
    COFD[122] = -3.36139991E-01;
    COFD[123] = 1.46507621E-02;
    COFD[124] = -1.14366381E+01;
    COFD[125] = 2.78323501E+00;
    COFD[126] = -1.51214064E-01;
    COFD[127] = 6.75150012E-03;
    COFD[128] = -1.34230272E+01;
    COFD[129] = 3.48624238E+00;
    COFD[130] = -2.41554467E-01;
    COFD[131] = 1.06263545E-02;
    COFD[132] = -1.34247866E+01;
    COFD[133] = 3.48624238E+00;
    COFD[134] = -2.41554467E-01;
    COFD[135] = 1.06263545E-02;
    COFD[136] = -1.46550083E+01;
    COFD[137] = 3.83606243E+00;
    COFD[138] = -2.86076532E-01;
    COFD[139] = 1.25205829E-02;
    COFD[140] = -1.46559141E+01;
    COFD[141] = 3.83606243E+00;
    COFD[142] = -2.86076532E-01;
    COFD[143] = 1.25205829E-02;
    COFD[144] = -1.95739570E+01;
    COFD[145] = 5.61113230E+00;
    COFD[146] = -4.90190187E-01;
    COFD[147] = 2.03260675E-02;
    COFD[148] = -1.46554748E+01;
    COFD[149] = 3.83606243E+00;
    COFD[150] = -2.86076532E-01;
    COFD[151] = 1.25205829E-02;
    COFD[152] = -1.43151174E+01;
    COFD[153] = 3.68038508E+00;
    COFD[154] = -2.65779346E-01;
    COFD[155] = 1.16360771E-02;
    COFD[156] = -1.57994893E+01;
    COFD[157] = 4.22225052E+00;
    COFD[158] = -3.35156428E-01;
    COFD[159] = 1.46104855E-02;
    COFD[160] = -1.97550088E+01;
    COFD[161] = 5.56931926E+00;
    COFD[162] = -4.89105511E-01;
    COFD[163] = 2.04493129E-02;
    COFD[164] = -1.76147026E+01;
    COFD[165] = 4.86049500E+00;
    COFD[166] = -4.12200578E-01;
    COFD[167] = 1.77160971E-02;
    COFD[168] = -1.57199037E+01;
    COFD[169] = 4.19936335E+00;
    COFD[170] = -3.32311009E-01;
    COFD[171] = 1.44921003E-02;
    COFD[172] = -1.72232223E+01;
    COFD[173] = 4.69060745E+00;
    COFD[174] = -3.92369888E-01;
    COFD[175] = 1.69459661E-02;
    COFD[176] = -1.82251914E+01;
    COFD[177] = 5.05237312E+00;
    COFD[178] = -4.35182396E-01;
    COFD[179] = 1.86363074E-02;
    COFD[180] = -1.94688688E+01;
    COFD[181] = 5.43830787E+00;
    COFD[182] = -4.75472880E-01;
    COFD[183] = 1.99909996E-02;
    COFD[184] = -1.79344949E+01;
    COFD[185] = 4.91373893E+00;
    COFD[186] = -4.18747629E-01;
    COFD[187] = 1.79856610E-02;
    COFD[188] = -1.32467319E+01;
    COFD[189] = 3.34156587E+00;
    COFD[190] = -2.22853306E-01;
    COFD[191] = 9.81883417E-03;
    COFD[192] = -1.34162893E+01;
    COFD[193] = 3.48624238E+00;
    COFD[194] = -2.41554467E-01;
    COFD[195] = 1.06263545E-02;
    COFD[196] = -1.97544450E+01;
    COFD[197] = 5.56931926E+00;
    COFD[198] = -4.89105511E-01;
    COFD[199] = 2.04493129E-02;
    COFD[200] = -1.57972369E+01;
    COFD[201] = 4.22225052E+00;
    COFD[202] = -3.35156428E-01;
    COFD[203] = 1.46104855E-02;
    COFD[204] = -1.57972369E+01;
    COFD[205] = 4.22225052E+00;
    COFD[206] = -3.35156428E-01;
    COFD[207] = 1.46104855E-02;
    COFD[208] = -1.72239172E+01;
    COFD[209] = 4.69060745E+00;
    COFD[210] = -3.92369888E-01;
    COFD[211] = 1.69459661E-02;
    COFD[212] = -1.79339327E+01;
    COFD[213] = 4.91373893E+00;
    COFD[214] = -4.18747629E-01;
    COFD[215] = 1.79856610E-02;
    COFD[216] = -1.54460820E+01;
    COFD[217] = 4.26819983E+00;
    COFD[218] = -3.40766379E-01;
    COFD[219] = 1.48393361E-02;
    COFD[220] = -1.94694048E+01;
    COFD[221] = 5.43830787E+00;
    COFD[222] = -4.75472880E-01;
    COFD[223] = 1.99909996E-02;
    COFD[224] = -1.94691430E+01;
    COFD[225] = 5.43830787E+00;
    COFD[226] = -4.75472880E-01;
    COFD[227] = 1.99909996E-02;
    COFD[228] = -1.98806372E+01;
    COFD[229] = 5.52555673E+00;
    COFD[230] = -4.84999851E-01;
    COFD[231] = 2.03334931E-02;
    COFD[232] = -1.16906297E+01;
    COFD[233] = 2.47469981E+00;
    COFD[234] = -1.10436257E-01;
    COFD[235] = 4.95273813E-03;
    COFD[236] = -1.14366381E+01;
    COFD[237] = 2.78323501E+00;
    COFD[238] = -1.51214064E-01;
    COFD[239] = 6.75150012E-03;
    COFD[240] = -1.03270606E+01;
    COFD[241] = 2.19285409E+00;
    COFD[242] = -7.54492786E-02;
    COFD[243] = 3.51398213E-03;
    COFD[244] = -1.09595712E+01;
    COFD[245] = 2.30836460E+00;
    COFD[246] = -8.76339315E-02;
    COFD[247] = 3.90878445E-03;
    COFD[248] = -1.09628982E+01;
    COFD[249] = 2.30836460E+00;
    COFD[250] = -8.76339315E-02;
    COFD[251] = 3.90878445E-03;
    COFD[252] = -1.18988955E+01;
    COFD[253] = 2.57507000E+00;
    COFD[254] = -1.24033737E-01;
    COFD[255] = 5.56694959E-03;
    COFD[256] = -1.19006548E+01;
    COFD[257] = 2.57507000E+00;
    COFD[258] = -1.24033737E-01;
    COFD[259] = 5.56694959E-03;
    COFD[260] = -1.71982995E+01;
    COFD[261] = 4.63881404E+00;
    COFD[262] = -3.86139633E-01;
    COFD[263] = 1.66955081E-02;
    COFD[264] = -1.18998012E+01;
    COFD[265] = 2.57507000E+00;
    COFD[266] = -1.24033737E-01;
    COFD[267] = 5.56694959E-03;
    COFD[268] = -1.17159737E+01;
    COFD[269] = 2.48123210E+00;
    COFD[270] = -1.11322604E-01;
    COFD[271] = 4.99282389E-03;
    COFD[272] = -1.25141260E+01;
    COFD[273] = 2.77873601E+00;
    COFD[274] = -1.50637360E-01;
    COFD[275] = 6.72684281E-03;
    COFD[276] = -1.60528285E+01;
    COFD[277] = 4.11188603E+00;
    COFD[278] = -3.21540884E-01;
    COFD[279] = 1.40482564E-02;
    COFD[280] = -1.37794315E+01;
    COFD[281] = 3.23973858E+00;
    COFD[282] = -2.09989036E-01;
    COFD[283] = 9.27667906E-03;
    COFD[284] = -1.24693568E+01;
    COFD[285] = 2.76686648E+00;
    COFD[286] = -1.49120141E-01;
    COFD[287] = 6.66220432E-03;
    COFD[288] = -1.34709807E+01;
    COFD[289] = 3.09379603E+00;
    COFD[290] = -1.91268635E-01;
    COFD[291] = 8.47480224E-03;
    COFD[292] = -1.42229194E+01;
    COFD[293] = 3.38669384E+00;
    COFD[294] = -2.28784122E-01;
    COFD[295] = 1.00790953E-02;
    COFD[296] = -1.57034851E+01;
    COFD[297] = 3.93614244E+00;
    COFD[298] = -2.99111497E-01;
    COFD[299] = 1.30888229E-02;
    COFD[300] = -1.39924781E+01;
    COFD[301] = 3.26384506E+00;
    COFD[302] = -2.12947087E-01;
    COFD[303] = 9.39743888E-03;
    COFD[304] = -1.08369481E+01;
    COFD[305] = 2.19094415E+00;
    COFD[306] = -7.11992510E-02;
    COFD[307] = 3.14105973E-03;
    COFD[308] = -1.09469245E+01;
    COFD[309] = 2.30836460E+00;
    COFD[310] = -8.76339315E-02;
    COFD[311] = 3.90878445E-03;
    COFD[312] = -1.60517370E+01;
    COFD[313] = 4.11188603E+00;
    COFD[314] = -3.21540884E-01;
    COFD[315] = 1.40482564E-02;
    COFD[316] = -1.25098960E+01;
    COFD[317] = 2.77873601E+00;
    COFD[318] = -1.50637360E-01;
    COFD[319] = 6.72684281E-03;
    COFD[320] = -1.25098960E+01;
    COFD[321] = 2.77873601E+00;
    COFD[322] = -1.50637360E-01;
    COFD[323] = 6.72684281E-03;
    COFD[324] = -1.34723215E+01;
    COFD[325] = 3.09379603E+00;
    COFD[326] = -1.91268635E-01;
    COFD[327] = 8.47480224E-03;
    COFD[328] = -1.39913897E+01;
    COFD[329] = 3.26384506E+00;
    COFD[330] = -2.12947087E-01;
    COFD[331] = 9.39743888E-03;
    COFD[332] = -1.22004324E+01;
    COFD[333] = 2.80725489E+00;
    COFD[334] = -1.54291406E-01;
    COFD[335] = 6.88290911E-03;
    COFD[336] = -1.57045332E+01;
    COFD[337] = 3.93614244E+00;
    COFD[338] = -2.99111497E-01;
    COFD[339] = 1.30888229E-02;
    COFD[340] = -1.57040212E+01;
    COFD[341] = 3.93614244E+00;
    COFD[342] = -2.99111497E-01;
    COFD[343] = 1.30888229E-02;
    COFD[344] = -1.61116686E+01;
    COFD[345] = 4.04227735E+00;
    COFD[346] = -3.12745253E-01;
    COFD[347] = 1.36756977E-02;
    COFD[348] = -1.40756935E+01;
    COFD[349] = 3.07549274E+00;
    COFD[350] = -1.88889344E-01;
    COFD[351] = 8.37152866E-03;
    COFD[352] = -1.34230272E+01;
    COFD[353] = 3.48624238E+00;
    COFD[354] = -2.41554467E-01;
    COFD[355] = 1.06263545E-02;
    COFD[356] = -1.09595712E+01;
    COFD[357] = 2.30836460E+00;
    COFD[358] = -8.76339315E-02;
    COFD[359] = 3.90878445E-03;
    COFD[360] = -1.32093628E+01;
    COFD[361] = 2.90778936E+00;
    COFD[362] = -1.67388544E-01;
    COFD[363] = 7.45220609E-03;
    COFD[364] = -1.32244035E+01;
    COFD[365] = 2.90778936E+00;
    COFD[366] = -1.67388544E-01;
    COFD[367] = 7.45220609E-03;
    COFD[368] = -1.43139231E+01;
    COFD[369] = 3.17651319E+00;
    COFD[370] = -2.02028974E-01;
    COFD[371] = 8.94232502E-03;
    COFD[372] = -1.43238998E+01;
    COFD[373] = 3.17651319E+00;
    COFD[374] = -2.02028974E-01;
    COFD[375] = 8.94232502E-03;
    COFD[376] = -1.94093572E+01;
    COFD[377] = 5.16013126E+00;
    COFD[378] = -4.46824543E-01;
    COFD[379] = 1.90464887E-02;
    COFD[380] = -1.43190389E+01;
    COFD[381] = 3.17651319E+00;
    COFD[382] = -2.02028974E-01;
    COFD[383] = 8.94232502E-03;
    COFD[384] = -1.40999008E+01;
    COFD[385] = 3.08120012E+00;
    COFD[386] = -1.89629903E-01;
    COFD[387] = 8.40361952E-03;
    COFD[388] = -1.50766130E+01;
    COFD[389] = 3.47945612E+00;
    COFD[390] = -2.40703722E-01;
    COFD[391] = 1.05907441E-02;
    COFD[392] = -1.94373127E+01;
    COFD[393] = 5.02567894E+00;
    COFD[394] = -4.32045169E-01;
    COFD[395] = 1.85132214E-02;
    COFD[396] = -1.70534856E+01;
    COFD[397] = 4.14240922E+00;
    COFD[398] = -3.25239774E-01;
    COFD[399] = 1.41980687E-02;
    COFD[400] = -1.50270339E+01;
    COFD[401] = 3.46140064E+00;
    COFD[402] = -2.38440092E-01;
    COFD[403] = 1.04960087E-02;
    COFD[404] = -1.65488358E+01;
    COFD[405] = 3.95035840E+00;
    COFD[406] = -3.00959418E-01;
    COFD[407] = 1.31692593E-02;
    COFD[408] = -1.74792112E+01;
    COFD[409] = 4.29676909E+00;
    COFD[410] = -3.44085306E-01;
    COFD[411] = 1.49671135E-02;
    COFD[412] = -1.90883268E+01;
    COFD[413] = 4.84384483E+00;
    COFD[414] = -4.10265575E-01;
    COFD[415] = 1.76414287E-02;
    COFD[416] = -1.72556499E+01;
    COFD[417] = 4.17889917E+00;
    COFD[418] = -3.29752510E-01;
    COFD[419] = 1.43850275E-02;
    COFD[420] = -1.30610083E+01;
    COFD[421] = 2.80913567E+00;
    COFD[422] = -1.54536855E-01;
    COFD[423] = 6.89359313E-03;
    COFD[424] = -1.31551788E+01;
    COFD[425] = 2.90778936E+00;
    COFD[426] = -1.67388544E-01;
    COFD[427] = 7.45220609E-03;
    COFD[428] = -1.94313116E+01;
    COFD[429] = 5.02567894E+00;
    COFD[430] = -4.32045169E-01;
    COFD[431] = 1.85132214E-02;
    COFD[432] = -1.50584249E+01;
    COFD[433] = 3.47945612E+00;
    COFD[434] = -2.40703722E-01;
    COFD[435] = 1.05907441E-02;
    COFD[436] = -1.50584249E+01;
    COFD[437] = 3.47945612E+00;
    COFD[438] = -2.40703722E-01;
    COFD[439] = 1.05907441E-02;
    COFD[440] = -1.65559787E+01;
    COFD[441] = 3.95035840E+00;
    COFD[442] = -3.00959418E-01;
    COFD[443] = 1.31692593E-02;
    COFD[444] = -1.72496634E+01;
    COFD[445] = 4.17889917E+00;
    COFD[446] = -3.29752510E-01;
    COFD[447] = 1.43850275E-02;
    COFD[448] = -1.49500357E+01;
    COFD[449] = 3.52327209E+00;
    COFD[450] = -2.46286208E-01;
    COFD[451] = 1.08285963E-02;
    COFD[452] = -1.90946745E+01;
    COFD[453] = 4.84384483E+00;
    COFD[454] = -4.10265575E-01;
    COFD[455] = 1.76414287E-02;
    COFD[456] = -1.90915649E+01;
    COFD[457] = 4.84384483E+00;
    COFD[458] = -4.10265575E-01;
    COFD[459] = 1.76414287E-02;
    COFD[460] = -1.95314689E+01;
    COFD[461] = 4.95249173E+00;
    COFD[462] = -4.23376552E-01;
    COFD[463] = 1.81703714E-02;
    COFD[464] = -1.40949196E+01;
    COFD[465] = 3.07549274E+00;
    COFD[466] = -1.88889344E-01;
    COFD[467] = 8.37152866E-03;
    COFD[468] = -1.34247866E+01;
    COFD[469] = 3.48624238E+00;
    COFD[470] = -2.41554467E-01;
    COFD[471] = 1.06263545E-02;
    COFD[472] = -1.09628982E+01;
    COFD[473] = 2.30836460E+00;
    COFD[474] = -8.76339315E-02;
    COFD[475] = 3.90878445E-03;
    COFD[476] = -1.32244035E+01;
    COFD[477] = 2.90778936E+00;
    COFD[478] = -1.67388544E-01;
    COFD[479] = 7.45220609E-03;
    COFD[480] = -1.32399106E+01;
    COFD[481] = 2.90778936E+00;
    COFD[482] = -1.67388544E-01;
    COFD[483] = 7.45220609E-03;
    COFD[484] = -1.43340796E+01;
    COFD[485] = 3.17651319E+00;
    COFD[486] = -2.02028974E-01;
    COFD[487] = 8.94232502E-03;
    COFD[488] = -1.43444709E+01;
    COFD[489] = 3.17651319E+00;
    COFD[490] = -2.02028974E-01;
    COFD[491] = 8.94232502E-03;
    COFD[492] = -1.94253036E+01;
    COFD[493] = 5.16013126E+00;
    COFD[494] = -4.46824543E-01;
    COFD[495] = 1.90464887E-02;
    COFD[496] = -1.43394069E+01;
    COFD[497] = 3.17651319E+00;
    COFD[498] = -2.02028974E-01;
    COFD[499] = 8.94232502E-03;
    COFD[500] = -1.41191261E+01;
    COFD[501] = 3.08120012E+00;
    COFD[502] = -1.89629903E-01;
    COFD[503] = 8.40361952E-03;
    COFD[504] = -1.50911794E+01;
    COFD[505] = 3.47945612E+00;
    COFD[506] = -2.40703722E-01;
    COFD[507] = 1.05907441E-02;
    COFD[508] = -1.94570287E+01;
    COFD[509] = 5.02567894E+00;
    COFD[510] = -4.32045169E-01;
    COFD[511] = 1.85132214E-02;
    COFD[512] = -1.70757047E+01;
    COFD[513] = 4.14240922E+00;
    COFD[514] = -3.25239774E-01;
    COFD[515] = 1.41980687E-02;
    COFD[516] = -1.50420953E+01;
    COFD[517] = 3.46140064E+00;
    COFD[518] = -2.38440092E-01;
    COFD[519] = 1.04960087E-02;
    COFD[520] = -1.65675362E+01;
    COFD[521] = 3.95035840E+00;
    COFD[522] = -3.00959418E-01;
    COFD[523] = 1.31692593E-02;
    COFD[524] = -1.74984476E+01;
    COFD[525] = 4.29676909E+00;
    COFD[526] = -3.44085306E-01;
    COFD[527] = 1.49671135E-02;
    COFD[528] = -1.91102652E+01;
    COFD[529] = 4.84384483E+00;
    COFD[530] = -4.10265575E-01;
    COFD[531] = 1.76414287E-02;
    COFD[532] = -1.72753760E+01;
    COFD[533] = 4.17889917E+00;
    COFD[534] = -3.29752510E-01;
    COFD[535] = 1.43850275E-02;
    COFD[536] = -1.30738796E+01;
    COFD[537] = 2.80913567E+00;
    COFD[538] = -1.54536855E-01;
    COFD[539] = 6.89359313E-03;
    COFD[540] = -1.31686537E+01;
    COFD[541] = 2.90778936E+00;
    COFD[542] = -1.67388544E-01;
    COFD[543] = 7.45220609E-03;
    COFD[544] = -1.94507876E+01;
    COFD[545] = 5.02567894E+00;
    COFD[546] = -4.32045169E-01;
    COFD[547] = 1.85132214E-02;
    COFD[548] = -1.50724636E+01;
    COFD[549] = 3.47945612E+00;
    COFD[550] = -2.40703722E-01;
    COFD[551] = 1.05907441E-02;
    COFD[552] = -1.50724636E+01;
    COFD[553] = 3.47945612E+00;
    COFD[554] = -2.40703722E-01;
    COFD[555] = 1.05907441E-02;
    COFD[556] = -1.65749533E+01;
    COFD[557] = 3.95035840E+00;
    COFD[558] = -3.00959418E-01;
    COFD[559] = 1.31692593E-02;
    COFD[560] = -1.72691500E+01;
    COFD[561] = 4.17889917E+00;
    COFD[562] = -3.29752510E-01;
    COFD[563] = 1.43850275E-02;
    COFD[564] = -1.49718233E+01;
    COFD[565] = 3.52327209E+00;
    COFD[566] = -2.46286208E-01;
    COFD[567] = 1.08285963E-02;
    COFD[568] = -1.91168996E+01;
    COFD[569] = 4.84384483E+00;
    COFD[570] = -4.10265575E-01;
    COFD[571] = 1.76414287E-02;
    COFD[572] = -1.91136491E+01;
    COFD[573] = 4.84384483E+00;
    COFD[574] = -4.10265575E-01;
    COFD[575] = 1.76414287E-02;
    COFD[576] = -1.95538303E+01;
    COFD[577] = 4.95249173E+00;
    COFD[578] = -4.23376552E-01;
    COFD[579] = 1.81703714E-02;
    COFD[580] = -1.52414485E+01;
    COFD[581] = 3.35922578E+00;
    COFD[582] = -2.25181399E-01;
    COFD[583] = 9.92132878E-03;
    COFD[584] = -1.46550083E+01;
    COFD[585] = 3.83606243E+00;
    COFD[586] = -2.86076532E-01;
    COFD[587] = 1.25205829E-02;
    COFD[588] = -1.18988955E+01;
    COFD[589] = 2.57507000E+00;
    COFD[590] = -1.24033737E-01;
    COFD[591] = 5.56694959E-03;
    COFD[592] = -1.43139231E+01;
    COFD[593] = 3.17651319E+00;
    COFD[594] = -2.02028974E-01;
    COFD[595] = 8.94232502E-03;
    COFD[596] = -1.43340796E+01;
    COFD[597] = 3.17651319E+00;
    COFD[598] = -2.02028974E-01;
    COFD[599] = 8.94232502E-03;
    COFD[600] = -1.55511344E+01;
    COFD[601] = 3.48070094E+00;
    COFD[602] = -2.40859499E-01;
    COFD[603] = 1.05972514E-02;
    COFD[604] = -1.55661750E+01;
    COFD[605] = 3.48070094E+00;
    COFD[606] = -2.40859499E-01;
    COFD[607] = 1.05972514E-02;
    COFD[608] = -2.12652533E+01;
    COFD[609] = 5.59961818E+00;
    COFD[610] = -4.91624858E-01;
    COFD[611] = 2.05035550E-02;
    COFD[612] = -1.55588279E+01;
    COFD[613] = 3.48070094E+00;
    COFD[614] = -2.40859499E-01;
    COFD[615] = 1.05972514E-02;
    COFD[616] = -1.52721107E+01;
    COFD[617] = 3.36790500E+00;
    COFD[618] = -2.26321740E-01;
    COFD[619] = 9.97135055E-03;
    COFD[620] = -1.63493345E+01;
    COFD[621] = 3.82388595E+00;
    COFD[622] = -2.84480724E-01;
    COFD[623] = 1.24506311E-02;
    COFD[624] = -2.08293255E+01;
    COFD[625] = 5.35267674E+00;
    COFD[626] = -4.69010505E-01;
    COFD[627] = 1.98979152E-02;
    COFD[628] = -1.84688406E+01;
    COFD[629] = 4.49330851E+00;
    COFD[630] = -3.68208715E-01;
    COFD[631] = 1.59565402E-02;
    COFD[632] = -1.62724462E+01;
    COFD[633] = 3.79163564E+00;
    COFD[634] = -2.80257365E-01;
    COFD[635] = 1.22656902E-02;
    COFD[636] = -1.78834935E+01;
    COFD[637] = 4.29613154E+00;
    COFD[638] = -3.44012526E-01;
    COFD[639] = 1.49643715E-02;
    COFD[640] = -1.89544778E+01;
    COFD[641] = 4.68595732E+00;
    COFD[642] = -3.91842840E-01;
    COFD[643] = 1.69262542E-02;
    COFD[644] = -2.05184870E+01;
    COFD[645] = 5.18417470E+00;
    COFD[646] = -4.49491573E-01;
    COFD[647] = 1.91438508E-02;
    COFD[648] = -1.86424545E+01;
    COFD[649] = 4.53572533E+00;
    COFD[650] = -3.73386925E-01;
    COFD[651] = 1.61678881E-02;
    COFD[652] = -1.40707473E+01;
    COFD[653] = 3.05837263E+00;
    COFD[654] = -1.86672802E-01;
    COFD[655] = 8.27575734E-03;
    COFD[656] = -1.42429085E+01;
    COFD[657] = 3.17651319E+00;
    COFD[658] = -2.02028974E-01;
    COFD[659] = 8.94232502E-03;
    COFD[660] = -2.08204449E+01;
    COFD[661] = 5.35267674E+00;
    COFD[662] = -4.69010505E-01;
    COFD[663] = 1.98979152E-02;
    COFD[664] = -1.63254691E+01;
    COFD[665] = 3.82388595E+00;
    COFD[666] = -2.84480724E-01;
    COFD[667] = 1.24506311E-02;
    COFD[668] = -1.63254691E+01;
    COFD[669] = 3.82388595E+00;
    COFD[670] = -2.84480724E-01;
    COFD[671] = 1.24506311E-02;
    COFD[672] = -1.78938745E+01;
    COFD[673] = 4.29613154E+00;
    COFD[674] = -3.44012526E-01;
    COFD[675] = 1.49643715E-02;
    COFD[676] = -1.86335932E+01;
    COFD[677] = 4.53572533E+00;
    COFD[678] = -3.73386925E-01;
    COFD[679] = 1.61678881E-02;
    COFD[680] = -1.64169433E+01;
    COFD[681] = 3.89309916E+00;
    COFD[682] = -2.93528188E-01;
    COFD[683] = 1.28463177E-02;
    COFD[684] = -2.05284752E+01;
    COFD[685] = 5.18417470E+00;
    COFD[686] = -4.49491573E-01;
    COFD[687] = 1.91438508E-02;
    COFD[688] = -2.05235731E+01;
    COFD[689] = 5.18417470E+00;
    COFD[690] = -4.49491573E-01;
    COFD[691] = 1.91438508E-02;
    COFD[692] = -2.09481051E+01;
    COFD[693] = 5.28755355E+00;
    COFD[694] = -4.61641920E-01;
    COFD[695] = 1.96208961E-02;
    COFD[696] = -1.52554761E+01;
    COFD[697] = 3.35922578E+00;
    COFD[698] = -2.25181399E-01;
    COFD[699] = 9.92132878E-03;
    COFD[700] = -1.46559141E+01;
    COFD[701] = 3.83606243E+00;
    COFD[702] = -2.86076532E-01;
    COFD[703] = 1.25205829E-02;
    COFD[704] = -1.19006548E+01;
    COFD[705] = 2.57507000E+00;
    COFD[706] = -1.24033737E-01;
    COFD[707] = 5.56694959E-03;
    COFD[708] = -1.43238998E+01;
    COFD[709] = 3.17651319E+00;
    COFD[710] = -2.02028974E-01;
    COFD[711] = 8.94232502E-03;
    COFD[712] = -1.43444709E+01;
    COFD[713] = 3.17651319E+00;
    COFD[714] = -2.02028974E-01;
    COFD[715] = 8.94232502E-03;
    COFD[716] = -1.55661750E+01;
    COFD[717] = 3.48070094E+00;
    COFD[718] = -2.40859499E-01;
    COFD[719] = 1.05972514E-02;
    COFD[720] = -1.55816822E+01;
    COFD[721] = 3.48070094E+00;
    COFD[722] = -2.40859499E-01;
    COFD[723] = 1.05972514E-02;
    COFD[724] = -2.06516336E+01;
    COFD[725] = 5.41688482E+00;
    COFD[726] = -4.73387188E-01;
    COFD[727] = 1.99280175E-02;
    COFD[728] = -1.55741053E+01;
    COFD[729] = 3.48070094E+00;
    COFD[730] = -2.40859499E-01;
    COFD[731] = 1.05972514E-02;
    COFD[732] = -1.52861376E+01;
    COFD[733] = 3.36790500E+00;
    COFD[734] = -2.26321740E-01;
    COFD[735] = 9.97135055E-03;
    COFD[736] = -1.63588981E+01;
    COFD[737] = 3.82388595E+00;
    COFD[738] = -2.84480724E-01;
    COFD[739] = 1.24506311E-02;
    COFD[740] = -2.08438809E+01;
    COFD[741] = 5.35267674E+00;
    COFD[742] = -4.69010505E-01;
    COFD[743] = 1.98979152E-02;
    COFD[744] = -1.84863000E+01;
    COFD[745] = 4.49330851E+00;
    COFD[746] = -3.68208715E-01;
    COFD[747] = 1.59565402E-02;
    COFD[748] = -1.62824412E+01;
    COFD[749] = 3.79163564E+00;
    COFD[750] = -2.80257365E-01;
    COFD[751] = 1.22656902E-02;
    COFD[752] = -1.78969684E+01;
    COFD[753] = 4.29613154E+00;
    COFD[754] = -3.44012526E-01;
    COFD[755] = 1.49643715E-02;
    COFD[756] = -1.89685165E+01;
    COFD[757] = 4.68595732E+00;
    COFD[758] = -3.91842840E-01;
    COFD[759] = 1.69262542E-02;
    COFD[760] = -2.05356023E+01;
    COFD[761] = 5.18417470E+00;
    COFD[762] = -4.49491573E-01;
    COFD[763] = 1.91438508E-02;
    COFD[764] = -1.86570209E+01;
    COFD[765] = 4.53572533E+00;
    COFD[766] = -3.73386925E-01;
    COFD[767] = 1.61678881E-02;
    COFD[768] = -1.40789009E+01;
    COFD[769] = 3.05837263E+00;
    COFD[770] = -1.86672802E-01;
    COFD[771] = 8.27575734E-03;
    COFD[772] = -1.42515527E+01;
    COFD[773] = 3.17651319E+00;
    COFD[774] = -2.02028974E-01;
    COFD[775] = 8.94232502E-03;
    COFD[776] = -2.08347403E+01;
    COFD[777] = 5.35267674E+00;
    COFD[778] = -4.69010505E-01;
    COFD[779] = 1.98979152E-02;
    COFD[780] = -1.63345829E+01;
    COFD[781] = 3.82388595E+00;
    COFD[782] = -2.84480724E-01;
    COFD[783] = 1.24506311E-02;
    COFD[784] = -1.63345829E+01;
    COFD[785] = 3.82388595E+00;
    COFD[786] = -2.84480724E-01;
    COFD[787] = 1.24506311E-02;
    COFD[788] = -1.79076361E+01;
    COFD[789] = 4.29613154E+00;
    COFD[790] = -3.44012526E-01;
    COFD[791] = 1.49643715E-02;
    COFD[792] = -1.86479000E+01;
    COFD[793] = 4.53572533E+00;
    COFD[794] = -3.73386925E-01;
    COFD[795] = 1.61678881E-02;
    COFD[796] = -1.64338757E+01;
    COFD[797] = 3.89309916E+00;
    COFD[798] = -2.93528188E-01;
    COFD[799] = 1.28463177E-02;
    COFD[800] = -2.05459419E+01;
    COFD[801] = 5.18417470E+00;
    COFD[802] = -4.49491573E-01;
    COFD[803] = 1.91438508E-02;
    COFD[804] = -2.05408665E+01;
    COFD[805] = 5.18417470E+00;
    COFD[806] = -4.49491573E-01;
    COFD[807] = 1.91438508E-02;
    COFD[808] = -2.09657408E+01;
    COFD[809] = 5.28755355E+00;
    COFD[810] = -4.61641920E-01;
    COFD[811] = 1.96208961E-02;
    COFD[812] = -2.10643259E+01;
    COFD[813] = 5.53614847E+00;
    COFD[814] = -4.86046736E-01;
    COFD[815] = 2.03659188E-02;
    COFD[816] = -1.95739570E+01;
    COFD[817] = 5.61113230E+00;
    COFD[818] = -4.90190187E-01;
    COFD[819] = 2.03260675E-02;
    COFD[820] = -1.71982995E+01;
    COFD[821] = 4.63881404E+00;
    COFD[822] = -3.86139633E-01;
    COFD[823] = 1.66955081E-02;
    COFD[824] = -1.94093572E+01;
    COFD[825] = 5.16013126E+00;
    COFD[826] = -4.46824543E-01;
    COFD[827] = 1.90464887E-02;
    COFD[828] = -1.94253036E+01;
    COFD[829] = 5.16013126E+00;
    COFD[830] = -4.46824543E-01;
    COFD[831] = 1.90464887E-02;
    COFD[832] = -2.12652533E+01;
    COFD[833] = 5.59961818E+00;
    COFD[834] = -4.91624858E-01;
    COFD[835] = 2.05035550E-02;
    COFD[836] = -2.06516336E+01;
    COFD[837] = 5.41688482E+00;
    COFD[838] = -4.73387188E-01;
    COFD[839] = 1.99280175E-02;
    COFD[840] = -1.19157919E+01;
    COFD[841] = 9.28955130E-01;
    COFD[842] = 2.42107090E-01;
    COFD[843] = -1.59823963E-02;
    COFD[844] = -2.06463744E+01;
    COFD[845] = 5.41688482E+00;
    COFD[846] = -4.73387188E-01;
    COFD[847] = 1.99280175E-02;
    COFD[848] = -2.11388331E+01;
    COFD[849] = 5.55529675E+00;
    COFD[850] = -4.87942518E-01;
    COFD[851] = 2.04249054E-02;
    COFD[852] = -2.12831323E+01;
    COFD[853] = 5.61184117E+00;
    COFD[854] = -4.90532156E-01;
    COFD[855] = 2.03507922E-02;
    COFD[856] = -1.77563250E+01;
    COFD[857] = 3.57475686E+00;
    COFD[858] = -1.56396297E-01;
    COFD[859] = 3.12157721E-03;
    COFD[860] = -2.07653719E+01;
    COFD[861] = 5.01092022E+00;
    COFD[862] = -3.77985635E-01;
    COFD[863] = 1.40968645E-02;
    COFD[864] = -2.14087397E+01;
    COFD[865] = 5.57282008E+00;
    COFD[866] = -4.76690890E-01;
    COFD[867] = 1.94000719E-02;
    COFD[868] = -2.15095980E+01;
    COFD[869] = 5.46737673E+00;
    COFD[870] = -4.55696085E-01;
    COFD[871] = 1.81982625E-02;
    COFD[872] = -2.08812333E+01;
    COFD[873] = 5.08859217E+00;
    COFD[874] = -3.90525428E-01;
    COFD[875] = 1.47376395E-02;
    COFD[876] = -1.87383952E+01;
    COFD[877] = 3.96926341E+00;
    COFD[878] = -2.16412264E-01;
    COFD[879] = 6.06012078E-03;
    COFD[880] = -2.12661865E+01;
    COFD[881] = 5.24930667E+00;
    COFD[882] = -4.17435088E-01;
    COFD[883] = 1.61434424E-02;
    COFD[884] = -1.92003817E+01;
    COFD[885] = 5.05708283E+00;
    COFD[886] = -4.35739290E-01;
    COFD[887] = 1.86583205E-02;
    COFD[888] = -1.93521390E+01;
    COFD[889] = 5.16013126E+00;
    COFD[890] = -4.46824543E-01;
    COFD[891] = 1.90464887E-02;
    COFD[892] = -1.77498543E+01;
    COFD[893] = 3.57475686E+00;
    COFD[894] = -1.56396297E-01;
    COFD[895] = 3.12157721E-03;
    COFD[896] = -2.12639214E+01;
    COFD[897] = 5.61184117E+00;
    COFD[898] = -4.90532156E-01;
    COFD[899] = 2.03507922E-02;
    COFD[900] = -2.12639214E+01;
    COFD[901] = 5.61184117E+00;
    COFD[902] = -4.90532156E-01;
    COFD[903] = 2.03507922E-02;
    COFD[904] = -2.15172770E+01;
    COFD[905] = 5.46737673E+00;
    COFD[906] = -4.55696085E-01;
    COFD[907] = 1.81982625E-02;
    COFD[908] = -2.12597312E+01;
    COFD[909] = 5.24930667E+00;
    COFD[910] = -4.17435088E-01;
    COFD[911] = 1.61434424E-02;
    COFD[912] = -2.10440675E+01;
    COFD[913] = 5.59806282E+00;
    COFD[914] = -4.87109535E-01;
    COFD[915] = 2.01370226E-02;
    COFD[916] = -1.87453067E+01;
    COFD[917] = 3.96926341E+00;
    COFD[918] = -2.16412264E-01;
    COFD[919] = 6.06012078E-03;
    COFD[920] = -1.87419199E+01;
    COFD[921] = 3.96926341E+00;
    COFD[922] = -2.16412264E-01;
    COFD[923] = 6.06012078E-03;
    COFD[924] = -1.84538368E+01;
    COFD[925] = 3.75912079E+00;
    COFD[926] = -1.84235105E-01;
    COFD[927] = 4.47800951E-03;
    COFD[928] = -1.52486273E+01;
    COFD[929] = 3.35922578E+00;
    COFD[930] = -2.25181399E-01;
    COFD[931] = 9.92132878E-03;
    COFD[932] = -1.46554748E+01;
    COFD[933] = 3.83606243E+00;
    COFD[934] = -2.86076532E-01;
    COFD[935] = 1.25205829E-02;
    COFD[936] = -1.18998012E+01;
    COFD[937] = 2.57507000E+00;
    COFD[938] = -1.24033737E-01;
    COFD[939] = 5.56694959E-03;
    COFD[940] = -1.43190389E+01;
    COFD[941] = 3.17651319E+00;
    COFD[942] = -2.02028974E-01;
    COFD[943] = 8.94232502E-03;
    COFD[944] = -1.43394069E+01;
    COFD[945] = 3.17651319E+00;
    COFD[946] = -2.02028974E-01;
    COFD[947] = 8.94232502E-03;
    COFD[948] = -1.55588279E+01;
    COFD[949] = 3.48070094E+00;
    COFD[950] = -2.40859499E-01;
    COFD[951] = 1.05972514E-02;
    COFD[952] = -1.55741053E+01;
    COFD[953] = 3.48070094E+00;
    COFD[954] = -2.40859499E-01;
    COFD[955] = 1.05972514E-02;
    COFD[956] = -2.06463744E+01;
    COFD[957] = 5.41688482E+00;
    COFD[958] = -4.73387188E-01;
    COFD[959] = 1.99280175E-02;
    COFD[960] = -1.55666415E+01;
    COFD[961] = 3.48070094E+00;
    COFD[962] = -2.40859499E-01;
    COFD[963] = 1.05972514E-02;
    COFD[964] = -1.52792891E+01;
    COFD[965] = 3.36790500E+00;
    COFD[966] = -2.26321740E-01;
    COFD[967] = 9.97135055E-03;
    COFD[968] = -1.63542394E+01;
    COFD[969] = 3.82388595E+00;
    COFD[970] = -2.84480724E-01;
    COFD[971] = 1.24506311E-02;
    COFD[972] = -2.08367725E+01;
    COFD[973] = 5.35267674E+00;
    COFD[974] = -4.69010505E-01;
    COFD[975] = 1.98979152E-02;
    COFD[976] = -1.84777607E+01;
    COFD[977] = 4.49330851E+00;
    COFD[978] = -3.68208715E-01;
    COFD[979] = 1.59565402E-02;
    COFD[980] = -1.62775714E+01;
    COFD[981] = 3.79163564E+00;
    COFD[982] = -2.80257365E-01;
    COFD[983] = 1.22656902E-02;
    COFD[984] = -1.78903913E+01;
    COFD[985] = 4.29613154E+00;
    COFD[986] = -3.44012526E-01;
    COFD[987] = 1.49643715E-02;
    COFD[988] = -1.89616623E+01;
    COFD[989] = 4.68595732E+00;
    COFD[990] = -3.91842840E-01;
    COFD[991] = 1.69262542E-02;
    COFD[992] = -2.05272328E+01;
    COFD[993] = 5.18417470E+00;
    COFD[994] = -4.49491573E-01;
    COFD[995] = 1.91438508E-02;
    COFD[996] = -1.86499071E+01;
    COFD[997] = 4.53572533E+00;
    COFD[998] = -3.73386925E-01;
    COFD[999] = 1.61678881E-02;
    COFD[1000] = -1.40749320E+01;
    COFD[1001] = 3.05837263E+00;
    COFD[1002] = -1.86672802E-01;
    COFD[1003] = 8.27575734E-03;
    COFD[1004] = -1.42473439E+01;
    COFD[1005] = 3.17651319E+00;
    COFD[1006] = -2.02028974E-01;
    COFD[1007] = 8.94232502E-03;
    COFD[1008] = -2.08277598E+01;
    COFD[1009] = 5.35267674E+00;
    COFD[1010] = -4.69010505E-01;
    COFD[1011] = 1.98979152E-02;
    COFD[1012] = -1.63301444E+01;
    COFD[1013] = 3.82388595E+00;
    COFD[1014] = -2.84480724E-01;
    COFD[1015] = 1.24506311E-02;
    COFD[1016] = -1.63301444E+01;
    COFD[1017] = 3.82388595E+00;
    COFD[1018] = -2.84480724E-01;
    COFD[1019] = 1.24506311E-02;
    COFD[1020] = -1.79009181E+01;
    COFD[1021] = 4.29613154E+00;
    COFD[1022] = -3.44012526E-01;
    COFD[1023] = 1.49643715E-02;
    COFD[1024] = -1.86409139E+01;
    COFD[1025] = 4.53572533E+00;
    COFD[1026] = -3.73386925E-01;
    COFD[1027] = 1.61678881E-02;
    COFD[1028] = -1.64255964E+01;
    COFD[1029] = 3.89309916E+00;
    COFD[1030] = -2.93528188E-01;
    COFD[1031] = 1.28463177E-02;
    COFD[1032] = -2.05373990E+01;
    COFD[1033] = 5.18417470E+00;
    COFD[1034] = -4.49491573E-01;
    COFD[1035] = 1.91438508E-02;
    COFD[1036] = -2.05324091E+01;
    COFD[1037] = 5.18417470E+00;
    COFD[1038] = -4.49491573E-01;
    COFD[1039] = 1.91438508E-02;
    COFD[1040] = -2.09571146E+01;
    COFD[1041] = 5.28755355E+00;
    COFD[1042] = -4.61641920E-01;
    COFD[1043] = 1.96208961E-02;
    COFD[1044] = -1.50031687E+01;
    COFD[1045] = 3.26223357E+00;
    COFD[1046] = -2.12746642E-01;
    COFD[1047] = 9.38912883E-03;
    COFD[1048] = -1.43151174E+01;
    COFD[1049] = 3.68038508E+00;
    COFD[1050] = -2.65779346E-01;
    COFD[1051] = 1.16360771E-02;
    COFD[1052] = -1.17159737E+01;
    COFD[1053] = 2.48123210E+00;
    COFD[1054] = -1.11322604E-01;
    COFD[1055] = 4.99282389E-03;
    COFD[1056] = -1.40999008E+01;
    COFD[1057] = 3.08120012E+00;
    COFD[1058] = -1.89629903E-01;
    COFD[1059] = 8.40361952E-03;
    COFD[1060] = -1.41191261E+01;
    COFD[1061] = 3.08120012E+00;
    COFD[1062] = -1.89629903E-01;
    COFD[1063] = 8.40361952E-03;
    COFD[1064] = -1.52721107E+01;
    COFD[1065] = 3.36790500E+00;
    COFD[1066] = -2.26321740E-01;
    COFD[1067] = 9.97135055E-03;
    COFD[1068] = -1.52861376E+01;
    COFD[1069] = 3.36790500E+00;
    COFD[1070] = -2.26321740E-01;
    COFD[1071] = 9.97135055E-03;
    COFD[1072] = -2.11388331E+01;
    COFD[1073] = 5.55529675E+00;
    COFD[1074] = -4.87942518E-01;
    COFD[1075] = 2.04249054E-02;
    COFD[1076] = -1.52792891E+01;
    COFD[1077] = 3.36790500E+00;
    COFD[1078] = -2.26321740E-01;
    COFD[1079] = 9.97135055E-03;
    COFD[1080] = -1.50233475E+01;
    COFD[1081] = 3.26660767E+00;
    COFD[1082] = -2.13287177E-01;
    COFD[1083] = 9.41137857E-03;
    COFD[1084] = -1.59863030E+01;
    COFD[1085] = 3.67388294E+00;
    COFD[1086] = -2.64990709E-01;
    COFD[1087] = 1.16042706E-02;
    COFD[1088] = -2.05128705E+01;
    COFD[1089] = 5.23843909E+00;
    COFD[1090] = -4.55815614E-01;
    COFD[1091] = 1.93898040E-02;
    COFD[1092] = -1.81735763E+01;
    COFD[1093] = 4.38391495E+00;
    COFD[1094] = -3.54941287E-01;
    COFD[1095] = 1.54195107E-02;
    COFD[1096] = -1.59525102E+01;
    COFD[1097] = 3.66023858E+00;
    COFD[1098] = -2.63401043E-01;
    COFD[1099] = 1.15432000E-02;
    COFD[1100] = -1.76285640E+01;
    COFD[1101] = 4.19935698E+00;
    COFD[1102] = -3.32310212E-01;
    COFD[1103] = 1.44920670E-02;
    COFD[1104] = -1.86157761E+01;
    COFD[1105] = 4.55689508E+00;
    COFD[1106] = -3.75937921E-01;
    COFD[1107] = 1.62703488E-02;
    COFD[1108] = -2.02922701E+01;
    COFD[1109] = 5.11106992E+00;
    COFD[1110] = -4.42047129E-01;
    COFD[1111] = 1.89042990E-02;
    COFD[1112] = -1.83538377E+01;
    COFD[1113] = 4.42828044E+00;
    COFD[1114] = -3.60417833E-01;
    COFD[1115] = 1.56455103E-02;
    COFD[1116] = -1.38862839E+01;
    COFD[1117] = 2.97564184E+00;
    COFD[1118] = -1.76025309E-01;
    COFD[1119] = 7.81869993E-03;
    COFD[1120] = -1.40318948E+01;
    COFD[1121] = 3.08120012E+00;
    COFD[1122] = -1.89629903E-01;
    COFD[1123] = 8.40361952E-03;
    COFD[1124] = -2.05045578E+01;
    COFD[1125] = 5.23843909E+00;
    COFD[1126] = -4.55815614E-01;
    COFD[1127] = 1.93898040E-02;
    COFD[1128] = -1.59634533E+01;
    COFD[1129] = 3.67388294E+00;
    COFD[1130] = -2.64990709E-01;
    COFD[1131] = 1.16042706E-02;
    COFD[1132] = -1.59634533E+01;
    COFD[1133] = 3.67388294E+00;
    COFD[1134] = -2.64990709E-01;
    COFD[1135] = 1.16042706E-02;
    COFD[1136] = -1.76383156E+01;
    COFD[1137] = 4.19935698E+00;
    COFD[1138] = -3.32310212E-01;
    COFD[1139] = 1.44920670E-02;
    COFD[1140] = -1.83455435E+01;
    COFD[1141] = 4.42828044E+00;
    COFD[1142] = -3.60417833E-01;
    COFD[1143] = 1.56455103E-02;
    COFD[1144] = -1.60261675E+01;
    COFD[1145] = 3.73312045E+00;
    COFD[1146] = -2.72579779E-01;
    COFD[1147] = 1.19290272E-02;
    COFD[1148] = -2.03015042E+01;
    COFD[1149] = 5.11106992E+00;
    COFD[1150] = -4.42047129E-01;
    COFD[1151] = 1.89042990E-02;
    COFD[1152] = -2.02969740E+01;
    COFD[1153] = 5.11106992E+00;
    COFD[1154] = -4.42047129E-01;
    COFD[1155] = 1.89042990E-02;
    COFD[1156] = -2.06066440E+01;
    COFD[1157] = 5.16748146E+00;
    COFD[1158] = -4.47594939E-01;
    COFD[1159] = 1.90724110E-02;
    COFD[1160] = -1.59633387E+01;
    COFD[1161] = 3.66853818E+00;
    COFD[1162] = -2.64346221E-01;
    COFD[1163] = 1.15784613E-02;
    COFD[1164] = -1.57994893E+01;
    COFD[1165] = 4.22225052E+00;
    COFD[1166] = -3.35156428E-01;
    COFD[1167] = 1.46104855E-02;
    COFD[1168] = -1.25141260E+01;
    COFD[1169] = 2.77873601E+00;
    COFD[1170] = -1.50637360E-01;
    COFD[1171] = 6.72684281E-03;
    COFD[1172] = -1.50766130E+01;
    COFD[1173] = 3.47945612E+00;
    COFD[1174] = -2.40703722E-01;
    COFD[1175] = 1.05907441E-02;
    COFD[1176] = -1.50911794E+01;
    COFD[1177] = 3.47945612E+00;
    COFD[1178] = -2.40703722E-01;
    COFD[1179] = 1.05907441E-02;
    COFD[1180] = -1.63493345E+01;
    COFD[1181] = 3.82388595E+00;
    COFD[1182] = -2.84480724E-01;
    COFD[1183] = 1.24506311E-02;
    COFD[1184] = -1.63588981E+01;
    COFD[1185] = 3.82388595E+00;
    COFD[1186] = -2.84480724E-01;
    COFD[1187] = 1.24506311E-02;
    COFD[1188] = -2.12831323E+01;
    COFD[1189] = 5.61184117E+00;
    COFD[1190] = -4.90532156E-01;
    COFD[1191] = 2.03507922E-02;
    COFD[1192] = -1.63542394E+01;
    COFD[1193] = 3.82388595E+00;
    COFD[1194] = -2.84480724E-01;
    COFD[1195] = 1.24506311E-02;
    COFD[1196] = -1.59863030E+01;
    COFD[1197] = 3.67388294E+00;
    COFD[1198] = -2.64990709E-01;
    COFD[1199] = 1.16042706E-02;
    COFD[1200] = -1.73374529E+01;
    COFD[1201] = 4.21416723E+00;
    COFD[1202] = -3.34163932E-01;
    COFD[1203] = 1.45697432E-02;
    COFD[1204] = -2.14449559E+01;
    COFD[1205] = 5.56531152E+00;
    COFD[1206] = -4.88789821E-01;
    COFD[1207] = 2.04437116E-02;
    COFD[1208] = -1.93276434E+01;
    COFD[1209] = 4.85015581E+00;
    COFD[1210] = -4.10945109E-01;
    COFD[1211] = 1.76651398E-02;
    COFD[1212] = -1.72738845E+01;
    COFD[1213] = 4.19029808E+00;
    COFD[1214] = -3.31177076E-01;
    COFD[1215] = 1.44446234E-02;
    COFD[1216] = -1.88463816E+01;
    COFD[1217] = 4.68393046E+00;
    COFD[1218] = -3.91610863E-01;
    COFD[1219] = 1.69174645E-02;
    COFD[1220] = -1.98646734E+01;
    COFD[1221] = 5.04367502E+00;
    COFD[1222] = -4.34153325E-01;
    COFD[1223] = 1.85956055E-02;
    COFD[1224] = -2.11606963E+01;
    COFD[1225] = 5.42846112E+00;
    COFD[1226] = -4.74321870E-01;
    COFD[1227] = 1.99459749E-02;
    COFD[1228] = -1.95552142E+01;
    COFD[1229] = 4.90255048E+00;
    COFD[1230] = -4.17368501E-01;
    COFD[1231] = 1.79287358E-02;
    COFD[1232] = -1.47715918E+01;
    COFD[1233] = 3.33113524E+00;
    COFD[1234] = -2.21479057E-01;
    COFD[1235] = 9.75837737E-03;
    COFD[1236] = -1.50240272E+01;
    COFD[1237] = 3.47945612E+00;
    COFD[1238] = -2.40703722E-01;
    COFD[1239] = 1.05907441E-02;
    COFD[1240] = -2.14391943E+01;
    COFD[1241] = 5.56531152E+00;
    COFD[1242] = -4.88789821E-01;
    COFD[1243] = 2.04437116E-02;
    COFD[1244] = -1.73198034E+01;
    COFD[1245] = 4.21416723E+00;
    COFD[1246] = -3.34163932E-01;
    COFD[1247] = 1.45697432E-02;
    COFD[1248] = -1.73198034E+01;
    COFD[1249] = 4.21416723E+00;
    COFD[1250] = -3.34163932E-01;
    COFD[1251] = 1.45697432E-02;
    COFD[1252] = -1.88532497E+01;
    COFD[1253] = 4.68393046E+00;
    COFD[1254] = -3.91610863E-01;
    COFD[1255] = 1.69174645E-02;
    COFD[1256] = -1.95494668E+01;
    COFD[1257] = 4.90255048E+00;
    COFD[1258] = -4.17368501E-01;
    COFD[1259] = 1.79287358E-02;
    COFD[1260] = -1.72828302E+01;
    COFD[1261] = 4.26063341E+00;
    COFD[1262] = -3.39848064E-01;
    COFD[1263] = 1.48021313E-02;
    COFD[1264] = -2.11667605E+01;
    COFD[1265] = 5.42846112E+00;
    COFD[1266] = -4.74321870E-01;
    COFD[1267] = 1.99459749E-02;
    COFD[1268] = -2.11637902E+01;
    COFD[1269] = 5.42846112E+00;
    COFD[1270] = -4.74321870E-01;
    COFD[1271] = 1.99459749E-02;
    COFD[1272] = -2.15588759E+01;
    COFD[1273] = 5.51982454E+00;
    COFD[1274] = -4.84452039E-01;
    COFD[1275] = 2.03175522E-02;
    COFD[1276] = -2.04833713E+01;
    COFD[1277] = 5.23112374E+00;
    COFD[1278] = -4.54967682E-01;
    COFD[1279] = 1.93570423E-02;
    COFD[1280] = -1.97550088E+01;
    COFD[1281] = 5.56931926E+00;
    COFD[1282] = -4.89105511E-01;
    COFD[1283] = 2.04493129E-02;
    COFD[1284] = -1.60528285E+01;
    COFD[1285] = 4.11188603E+00;
    COFD[1286] = -3.21540884E-01;
    COFD[1287] = 1.40482564E-02;
    COFD[1288] = -1.94373127E+01;
    COFD[1289] = 5.02567894E+00;
    COFD[1290] = -4.32045169E-01;
    COFD[1291] = 1.85132214E-02;
    COFD[1292] = -1.94570287E+01;
    COFD[1293] = 5.02567894E+00;
    COFD[1294] = -4.32045169E-01;
    COFD[1295] = 1.85132214E-02;
    COFD[1296] = -2.08293255E+01;
    COFD[1297] = 5.35267674E+00;
    COFD[1298] = -4.69010505E-01;
    COFD[1299] = 1.98979152E-02;
    COFD[1300] = -2.08438809E+01;
    COFD[1301] = 5.35267674E+00;
    COFD[1302] = -4.69010505E-01;
    COFD[1303] = 1.98979152E-02;
    COFD[1304] = -1.77563250E+01;
    COFD[1305] = 3.57475686E+00;
    COFD[1306] = -1.56396297E-01;
    COFD[1307] = 3.12157721E-03;
    COFD[1308] = -2.08367725E+01;
    COFD[1309] = 5.35267674E+00;
    COFD[1310] = -4.69010505E-01;
    COFD[1311] = 1.98979152E-02;
    COFD[1312] = -2.05128705E+01;
    COFD[1313] = 5.23843909E+00;
    COFD[1314] = -4.55815614E-01;
    COFD[1315] = 1.93898040E-02;
    COFD[1316] = -2.14449559E+01;
    COFD[1317] = 5.56531152E+00;
    COFD[1318] = -4.88789821E-01;
    COFD[1319] = 2.04437116E-02;
    COFD[1320] = -1.90499441E+01;
    COFD[1321] = 3.99221757E+00;
    COFD[1322] = -2.19854880E-01;
    COFD[1323] = 6.22736279E-03;
    COFD[1324] = -2.19317743E+01;
    COFD[1325] = 5.45216133E+00;
    COFD[1326] = -4.52916925E-01;
    COFD[1327] = 1.80456400E-02;
    COFD[1328] = -2.14082453E+01;
    COFD[1329] = 5.55346617E+00;
    COFD[1330] = -4.87783156E-01;
    COFD[1331] = 2.04210886E-02;
    COFD[1332] = -2.20036369E+01;
    COFD[1333] = 5.55935694E+00;
    COFD[1334] = -4.74154740E-01;
    COFD[1335] = 1.92584304E-02;
    COFD[1336] = -2.16379567E+01;
    COFD[1337] = 5.29019717E+00;
    COFD[1338] = -4.24502606E-01;
    COFD[1339] = 1.65197343E-02;
    COFD[1340] = -2.01015340E+01;
    COFD[1341] = 4.41511629E+00;
    COFD[1342] = -2.84086963E-01;
    COFD[1343] = 9.37586971E-03;
    COFD[1344] = -2.19399793E+01;
    COFD[1345] = 5.41841631E+00;
    COFD[1346] = -4.46818971E-01;
    COFD[1347] = 1.77127652E-02;
    COFD[1348] = -1.90576320E+01;
    COFD[1349] = 4.86821670E+00;
    COFD[1350] = -4.13144121E-01;
    COFD[1351] = 1.77546701E-02;
    COFD[1352] = -1.93677186E+01;
    COFD[1353] = 5.02567894E+00;
    COFD[1354] = -4.32045169E-01;
    COFD[1355] = 1.85132214E-02;
    COFD[1356] = -1.90413348E+01;
    COFD[1357] = 3.99221757E+00;
    COFD[1358] = -2.19854880E-01;
    COFD[1359] = 6.22736279E-03;
    COFD[1360] = -2.14215700E+01;
    COFD[1361] = 5.56531152E+00;
    COFD[1362] = -4.88789821E-01;
    COFD[1363] = 2.04437116E-02;
    COFD[1364] = -2.14215700E+01;
    COFD[1365] = 5.56531152E+00;
    COFD[1366] = -4.88789821E-01;
    COFD[1367] = 2.04437116E-02;
    COFD[1368] = -2.20137178E+01;
    COFD[1369] = 5.55935694E+00;
    COFD[1370] = -4.74154740E-01;
    COFD[1371] = 1.92584304E-02;
    COFD[1372] = -2.19313890E+01;
    COFD[1373] = 5.41841631E+00;
    COFD[1374] = -4.46818971E-01;
    COFD[1375] = 1.77127652E-02;
    COFD[1376] = -2.14303479E+01;
    COFD[1377] = 5.59268435E+00;
    COFD[1378] = -4.91232974E-01;
    COFD[1379] = 2.05064746E-02;
    COFD[1380] = -2.01111595E+01;
    COFD[1381] = 4.41511629E+00;
    COFD[1382] = -2.84086963E-01;
    COFD[1383] = 9.37586971E-03;
    COFD[1384] = -2.01064363E+01;
    COFD[1385] = 4.41511629E+00;
    COFD[1386] = -2.84086963E-01;
    COFD[1387] = 9.37586971E-03;
    COFD[1388] = -1.97649065E+01;
    COFD[1389] = 4.18758010E+00;
    COFD[1390] = -2.49327776E-01;
    COFD[1391] = 7.66559103E-03;
    COFD[1392] = -1.81432461E+01;
    COFD[1393] = 4.37565431E+00;
    COFD[1394] = -3.53906025E-01;
    COFD[1395] = 1.53760786E-02;
    COFD[1396] = -1.76147026E+01;
    COFD[1397] = 4.86049500E+00;
    COFD[1398] = -4.12200578E-01;
    COFD[1399] = 1.77160971E-02;
    COFD[1400] = -1.37794315E+01;
    COFD[1401] = 3.23973858E+00;
    COFD[1402] = -2.09989036E-01;
    COFD[1403] = 9.27667906E-03;
    COFD[1404] = -1.70534856E+01;
    COFD[1405] = 4.14240922E+00;
    COFD[1406] = -3.25239774E-01;
    COFD[1407] = 1.41980687E-02;
    COFD[1408] = -1.70757047E+01;
    COFD[1409] = 4.14240922E+00;
    COFD[1410] = -3.25239774E-01;
    COFD[1411] = 1.41980687E-02;
    COFD[1412] = -1.84688406E+01;
    COFD[1413] = 4.49330851E+00;
    COFD[1414] = -3.68208715E-01;
    COFD[1415] = 1.59565402E-02;
    COFD[1416] = -1.84863000E+01;
    COFD[1417] = 4.49330851E+00;
    COFD[1418] = -3.68208715E-01;
    COFD[1419] = 1.59565402E-02;
    COFD[1420] = -2.07653719E+01;
    COFD[1421] = 5.01092022E+00;
    COFD[1422] = -3.77985635E-01;
    COFD[1423] = 1.40968645E-02;
    COFD[1424] = -1.84777607E+01;
    COFD[1425] = 4.49330851E+00;
    COFD[1426] = -3.68208715E-01;
    COFD[1427] = 1.59565402E-02;
    COFD[1428] = -1.81735763E+01;
    COFD[1429] = 4.38391495E+00;
    COFD[1430] = -3.54941287E-01;
    COFD[1431] = 1.54195107E-02;
    COFD[1432] = -1.93276434E+01;
    COFD[1433] = 4.85015581E+00;
    COFD[1434] = -4.10945109E-01;
    COFD[1435] = 1.76651398E-02;
    COFD[1436] = -2.19317743E+01;
    COFD[1437] = 5.45216133E+00;
    COFD[1438] = -4.52916925E-01;
    COFD[1439] = 1.80456400E-02;
    COFD[1440] = -2.13425698E+01;
    COFD[1441] = 5.40460130E+00;
    COFD[1442] = -4.72718910E-01;
    COFD[1443] = 1.99362717E-02;
    COFD[1444] = -1.92867554E+01;
    COFD[1445] = 4.83375900E+00;
    COFD[1446] = -4.09146560E-01;
    COFD[1447] = 1.76006599E-02;
    COFD[1448] = -2.09191285E+01;
    COFD[1449] = 5.30153901E+00;
    COFD[1450] = -4.63335119E-01;
    COFD[1451] = 1.96897053E-02;
    COFD[1452] = -2.16802612E+01;
    COFD[1453] = 5.52918296E+00;
    COFD[1454] = -4.85360709E-01;
    COFD[1455] = 2.03448006E-02;
    COFD[1456] = -2.22116706E+01;
    COFD[1457] = 5.54251230E+00;
    COFD[1458] = -4.70946314E-01;
    COFD[1459] = 1.90785869E-02;
    COFD[1460] = -2.14326461E+01;
    COFD[1461] = 5.41729961E+00;
    COFD[1462] = -4.73400281E-01;
    COFD[1463] = 1.99269567E-02;
    COFD[1464] = -1.67390425E+01;
    COFD[1465] = 4.00828594E+00;
    COFD[1466] = -3.08414344E-01;
    COFD[1467] = 1.34907430E-02;
    COFD[1468] = -1.69758891E+01;
    COFD[1469] = 4.14240922E+00;
    COFD[1470] = -3.25239774E-01;
    COFD[1471] = 1.41980687E-02;
    COFD[1472] = -2.19215555E+01;
    COFD[1473] = 5.45216133E+00;
    COFD[1474] = -4.52916925E-01;
    COFD[1475] = 1.80456400E-02;
    COFD[1476] = -1.93015555E+01;
    COFD[1477] = 4.85015581E+00;
    COFD[1478] = -4.10945109E-01;
    COFD[1479] = 1.76651398E-02;
    COFD[1480] = -1.93015555E+01;
    COFD[1481] = 4.85015581E+00;
    COFD[1482] = -4.10945109E-01;
    COFD[1483] = 1.76651398E-02;
    COFD[1484] = -2.09309753E+01;
    COFD[1485] = 5.30153901E+00;
    COFD[1486] = -4.63335119E-01;
    COFD[1487] = 1.96897053E-02;
    COFD[1488] = -2.14224484E+01;
    COFD[1489] = 5.41729961E+00;
    COFD[1490] = -4.73400281E-01;
    COFD[1491] = 1.99269567E-02;
    COFD[1492] = -1.94485982E+01;
    COFD[1493] = 4.91446566E+00;
    COFD[1494] = -4.18837152E-01;
    COFD[1495] = 1.79893537E-02;
    COFD[1496] = -2.22235122E+01;
    COFD[1497] = 5.54251230E+00;
    COFD[1498] = -4.70946314E-01;
    COFD[1499] = 1.90785869E-02;
    COFD[1500] = -2.22176950E+01;
    COFD[1501] = 5.54251230E+00;
    COFD[1502] = -4.70946314E-01;
    COFD[1503] = 1.90785869E-02;
    COFD[1504] = -2.23098172E+01;
    COFD[1505] = 5.49916900E+00;
    COFD[1506] = -4.61818485E-01;
    COFD[1507] = 1.85431163E-02;
    COFD[1508] = -1.59327297E+01;
    COFD[1509] = 3.65620899E+00;
    COFD[1510] = -2.62933804E-01;
    COFD[1511] = 1.15253223E-02;
    COFD[1512] = -1.57199037E+01;
    COFD[1513] = 4.19936335E+00;
    COFD[1514] = -3.32311009E-01;
    COFD[1515] = 1.44921003E-02;
    COFD[1516] = -1.24693568E+01;
    COFD[1517] = 2.76686648E+00;
    COFD[1518] = -1.49120141E-01;
    COFD[1519] = 6.66220432E-03;
    COFD[1520] = -1.50270339E+01;
    COFD[1521] = 3.46140064E+00;
    COFD[1522] = -2.38440092E-01;
    COFD[1523] = 1.04960087E-02;
    COFD[1524] = -1.50420953E+01;
    COFD[1525] = 3.46140064E+00;
    COFD[1526] = -2.38440092E-01;
    COFD[1527] = 1.04960087E-02;
    COFD[1528] = -1.62724462E+01;
    COFD[1529] = 3.79163564E+00;
    COFD[1530] = -2.80257365E-01;
    COFD[1531] = 1.22656902E-02;
    COFD[1532] = -1.62824412E+01;
    COFD[1533] = 3.79163564E+00;
    COFD[1534] = -2.80257365E-01;
    COFD[1535] = 1.22656902E-02;
    COFD[1536] = -2.14087397E+01;
    COFD[1537] = 5.57282008E+00;
    COFD[1538] = -4.76690890E-01;
    COFD[1539] = 1.94000719E-02;
    COFD[1540] = -1.62775714E+01;
    COFD[1541] = 3.79163564E+00;
    COFD[1542] = -2.80257365E-01;
    COFD[1543] = 1.22656902E-02;
    COFD[1544] = -1.59525102E+01;
    COFD[1545] = 3.66023858E+00;
    COFD[1546] = -2.63401043E-01;
    COFD[1547] = 1.15432000E-02;
    COFD[1548] = -1.72738845E+01;
    COFD[1549] = 4.19029808E+00;
    COFD[1550] = -3.31177076E-01;
    COFD[1551] = 1.44446234E-02;
    COFD[1552] = -2.14082453E+01;
    COFD[1553] = 5.55346617E+00;
    COFD[1554] = -4.87783156E-01;
    COFD[1555] = 2.04210886E-02;
    COFD[1556] = -1.92867554E+01;
    COFD[1557] = 4.83375900E+00;
    COFD[1558] = -4.09146560E-01;
    COFD[1559] = 1.76006599E-02;
    COFD[1560] = -1.72167708E+01;
    COFD[1561] = 4.16886779E+00;
    COFD[1562] = -3.28518156E-01;
    COFD[1563] = 1.43341626E-02;
    COFD[1564] = -1.87897298E+01;
    COFD[1565] = 4.66162351E+00;
    COFD[1566] = -3.88920477E-01;
    COFD[1567] = 1.68089648E-02;
    COFD[1568] = -1.98075055E+01;
    COFD[1569] = 5.02169524E+00;
    COFD[1570] = -4.31582804E-01;
    COFD[1571] = 1.84953568E-02;
    COFD[1572] = -2.11309207E+01;
    COFD[1573] = 5.41773516E+00;
    COFD[1574] = -4.73414338E-01;
    COFD[1575] = 1.99258685E-02;
    COFD[1576] = -1.94823660E+01;
    COFD[1577] = 4.87333294E+00;
    COFD[1578] = -4.13769241E-01;
    COFD[1579] = 1.77802244E-02;
    COFD[1580] = -1.47048024E+01;
    COFD[1581] = 3.30594991E+00;
    COFD[1582] = -2.18182207E-01;
    COFD[1583] = 9.61429447E-03;
    COFD[1584] = -1.49727799E+01;
    COFD[1585] = 3.46140064E+00;
    COFD[1586] = -2.38440092E-01;
    COFD[1587] = 1.04960087E-02;
    COFD[1588] = -2.14022336E+01;
    COFD[1589] = 5.55346617E+00;
    COFD[1590] = -4.87783156E-01;
    COFD[1591] = 2.04210886E-02;
    COFD[1592] = -1.72556729E+01;
    COFD[1593] = 4.19029808E+00;
    COFD[1594] = -3.31177076E-01;
    COFD[1595] = 1.44446234E-02;
    COFD[1596] = -1.72556729E+01;
    COFD[1597] = 4.19029808E+00;
    COFD[1598] = -3.31177076E-01;
    COFD[1599] = 1.44446234E-02;
    COFD[1600] = -1.87968848E+01;
    COFD[1601] = 4.66162351E+00;
    COFD[1602] = -3.88920477E-01;
    COFD[1603] = 1.68089648E-02;
    COFD[1604] = -1.94763688E+01;
    COFD[1605] = 4.87333294E+00;
    COFD[1606] = -4.13769241E-01;
    COFD[1607] = 1.77802244E-02;
    COFD[1608] = -1.72316148E+01;
    COFD[1609] = 4.24011069E+00;
    COFD[1610] = -3.37339810E-01;
    COFD[1611] = 1.46996679E-02;
    COFD[1612] = -2.11372811E+01;
    COFD[1613] = 5.41773516E+00;
    COFD[1614] = -4.73414338E-01;
    COFD[1615] = 1.99258685E-02;
    COFD[1616] = -2.11341653E+01;
    COFD[1617] = 5.41773516E+00;
    COFD[1618] = -4.73414338E-01;
    COFD[1619] = 1.99258685E-02;
    COFD[1620] = -2.15067581E+01;
    COFD[1621] = 5.49964831E+00;
    COFD[1622] = -4.82275380E-01;
    COFD[1623] = 2.02405072E-02;
    COFD[1624] = -1.76002031E+01;
    COFD[1625] = 4.19171952E+00;
    COFD[1626] = -3.31354810E-01;
    COFD[1627] = 1.44520623E-02;
    COFD[1628] = -1.72232223E+01;
    COFD[1629] = 4.69060745E+00;
    COFD[1630] = -3.92369888E-01;
    COFD[1631] = 1.69459661E-02;
    COFD[1632] = -1.34709807E+01;
    COFD[1633] = 3.09379603E+00;
    COFD[1634] = -1.91268635E-01;
    COFD[1635] = 8.47480224E-03;
    COFD[1636] = -1.65488358E+01;
    COFD[1637] = 3.95035840E+00;
    COFD[1638] = -3.00959418E-01;
    COFD[1639] = 1.31692593E-02;
    COFD[1640] = -1.65675362E+01;
    COFD[1641] = 3.95035840E+00;
    COFD[1642] = -3.00959418E-01;
    COFD[1643] = 1.31692593E-02;
    COFD[1644] = -1.78834935E+01;
    COFD[1645] = 4.29613154E+00;
    COFD[1646] = -3.44012526E-01;
    COFD[1647] = 1.49643715E-02;
    COFD[1648] = -1.78969684E+01;
    COFD[1649] = 4.29613154E+00;
    COFD[1650] = -3.44012526E-01;
    COFD[1651] = 1.49643715E-02;
    COFD[1652] = -2.15095980E+01;
    COFD[1653] = 5.46737673E+00;
    COFD[1654] = -4.55696085E-01;
    COFD[1655] = 1.81982625E-02;
    COFD[1656] = -1.78903913E+01;
    COFD[1657] = 4.29613154E+00;
    COFD[1658] = -3.44012526E-01;
    COFD[1659] = 1.49643715E-02;
    COFD[1660] = -1.76285640E+01;
    COFD[1661] = 4.19935698E+00;
    COFD[1662] = -3.32310212E-01;
    COFD[1663] = 1.44920670E-02;
    COFD[1664] = -1.88463816E+01;
    COFD[1665] = 4.68393046E+00;
    COFD[1666] = -3.91610863E-01;
    COFD[1667] = 1.69174645E-02;
    COFD[1668] = -2.20036369E+01;
    COFD[1669] = 5.55935694E+00;
    COFD[1670] = -4.74154740E-01;
    COFD[1671] = 1.92584304E-02;
    COFD[1672] = -2.09191285E+01;
    COFD[1673] = 5.30153901E+00;
    COFD[1674] = -4.63335119E-01;
    COFD[1675] = 1.96897053E-02;
    COFD[1676] = -1.87897298E+01;
    COFD[1677] = 4.66162351E+00;
    COFD[1678] = -3.88920477E-01;
    COFD[1679] = 1.68089648E-02;
    COFD[1680] = -2.03766950E+01;
    COFD[1681] = 5.13263469E+00;
    COFD[1682] = -4.44457285E-01;
    COFD[1683] = 1.89932102E-02;
    COFD[1684] = -2.12121370E+01;
    COFD[1685] = 5.39823225E+00;
    COFD[1686] = -4.72294645E-01;
    COFD[1687] = 1.99340225E-02;
    COFD[1688] = -2.21065306E+01;
    COFD[1689] = 5.58360799E+00;
    COFD[1690] = -4.82701436E-01;
    COFD[1691] = 1.98437922E-02;
    COFD[1692] = -2.10944088E+01;
    COFD[1693] = 5.34286099E+00;
    COFD[1694] = -4.68100992E-01;
    COFD[1695] = 1.98731399E-02;
    COFD[1696] = -1.61101966E+01;
    COFD[1697] = 3.75910622E+00;
    COFD[1698] = -2.75986578E-01;
    COFD[1699] = 1.20782843E-02;
    COFD[1700] = -1.64825368E+01;
    COFD[1701] = 3.95035840E+00;
    COFD[1702] = -3.00959418E-01;
    COFD[1703] = 1.31692593E-02;
    COFD[1704] = -2.19956352E+01;
    COFD[1705] = 5.55935694E+00;
    COFD[1706] = -4.74154740E-01;
    COFD[1707] = 1.92584304E-02;
    COFD[1708] = -1.88241079E+01;
    COFD[1709] = 4.68393046E+00;
    COFD[1710] = -3.91610863E-01;
    COFD[1711] = 1.69174645E-02;
    COFD[1712] = -1.88241079E+01;
    COFD[1713] = 4.68393046E+00;
    COFD[1714] = -3.91610863E-01;
    COFD[1715] = 1.69174645E-02;
    COFD[1716] = -2.03861000E+01;
    COFD[1717] = 5.13263469E+00;
    COFD[1718] = -4.44457285E-01;
    COFD[1719] = 1.89932102E-02;
    COFD[1720] = -2.10864251E+01;
    COFD[1721] = 5.34286099E+00;
    COFD[1722] = -4.68100992E-01;
    COFD[1723] = 1.98731399E-02;
    COFD[1724] = -1.88646070E+01;
    COFD[1725] = 4.72476764E+00;
    COFD[1726] = -3.96306836E-01;
    COFD[1727] = 1.70964541E-02;
    COFD[1728] = -2.21153597E+01;
    COFD[1729] = 5.58360799E+00;
    COFD[1730] = -4.82701436E-01;
    COFD[1731] = 1.98437922E-02;
    COFD[1732] = -2.21110290E+01;
    COFD[1733] = 5.58360799E+00;
    COFD[1734] = -4.82701436E-01;
    COFD[1735] = 1.98437922E-02;
    COFD[1736] = -2.22982472E+01;
    COFD[1737] = 5.58471203E+00;
    COFD[1738] = -4.79905311E-01;
    COFD[1739] = 1.96058913E-02;
    COFD[1740] = -1.85864144E+01;
    COFD[1741] = 4.54915847E+00;
    COFD[1742] = -3.75000738E-01;
    COFD[1743] = 1.62324821E-02;
    COFD[1744] = -1.82251914E+01;
    COFD[1745] = 5.05237312E+00;
    COFD[1746] = -4.35182396E-01;
    COFD[1747] = 1.86363074E-02;
    COFD[1748] = -1.42229194E+01;
    COFD[1749] = 3.38669384E+00;
    COFD[1750] = -2.28784122E-01;
    COFD[1751] = 1.00790953E-02;
    COFD[1752] = -1.74792112E+01;
    COFD[1753] = 4.29676909E+00;
    COFD[1754] = -3.44085306E-01;
    COFD[1755] = 1.49671135E-02;
    COFD[1756] = -1.74984476E+01;
    COFD[1757] = 4.29676909E+00;
    COFD[1758] = -3.44085306E-01;
    COFD[1759] = 1.49671135E-02;
    COFD[1760] = -1.89544778E+01;
    COFD[1761] = 4.68595732E+00;
    COFD[1762] = -3.91842840E-01;
    COFD[1763] = 1.69262542E-02;
    COFD[1764] = -1.89685165E+01;
    COFD[1765] = 4.68595732E+00;
    COFD[1766] = -3.91842840E-01;
    COFD[1767] = 1.69262542E-02;
    COFD[1768] = -2.08812333E+01;
    COFD[1769] = 5.08859217E+00;
    COFD[1770] = -3.90525428E-01;
    COFD[1771] = 1.47376395E-02;
    COFD[1772] = -1.89616623E+01;
    COFD[1773] = 4.68595732E+00;
    COFD[1774] = -3.91842840E-01;
    COFD[1775] = 1.69262542E-02;
    COFD[1776] = -1.86157761E+01;
    COFD[1777] = 4.55689508E+00;
    COFD[1778] = -3.75937921E-01;
    COFD[1779] = 1.62703488E-02;
    COFD[1780] = -1.98646734E+01;
    COFD[1781] = 5.04367502E+00;
    COFD[1782] = -4.34153325E-01;
    COFD[1783] = 1.85956055E-02;
    COFD[1784] = -2.16379567E+01;
    COFD[1785] = 5.29019717E+00;
    COFD[1786] = -4.24502606E-01;
    COFD[1787] = 1.65197343E-02;
    COFD[1788] = -2.16802612E+01;
    COFD[1789] = 5.52918296E+00;
    COFD[1790] = -4.85360709E-01;
    COFD[1791] = 2.03448006E-02;
    COFD[1792] = -1.98075055E+01;
    COFD[1793] = 5.02169524E+00;
    COFD[1794] = -4.31582804E-01;
    COFD[1795] = 1.84953568E-02;
    COFD[1796] = -2.12121370E+01;
    COFD[1797] = 5.39823225E+00;
    COFD[1798] = -4.72294645E-01;
    COFD[1799] = 1.99340225E-02;
    COFD[1800] = -2.19327397E+01;
    COFD[1801] = 5.60638188E+00;
    COFD[1802] = -4.91272522E-01;
    COFD[1803] = 2.04396264E-02;
    COFD[1804] = -2.20453723E+01;
    COFD[1805] = 5.44448440E+00;
    COFD[1806] = -4.51529024E-01;
    COFD[1807] = 1.79698119E-02;
    COFD[1808] = -2.18273547E+01;
    COFD[1809] = 5.55753905E+00;
    COFD[1810] = -4.88136714E-01;
    COFD[1811] = 2.04294957E-02;
    COFD[1812] = -1.71884218E+01;
    COFD[1813] = 4.17190426E+00;
    COFD[1814] = -3.28894681E-01;
    COFD[1815] = 1.43498101E-02;
    COFD[1816] = -1.74111692E+01;
    COFD[1817] = 4.29676909E+00;
    COFD[1818] = -3.44085306E-01;
    COFD[1819] = 1.49671135E-02;
    COFD[1820] = -2.16296373E+01;
    COFD[1821] = 5.29019717E+00;
    COFD[1822] = -4.24502606E-01;
    COFD[1823] = 1.65197343E-02;
    COFD[1824] = -1.98418115E+01;
    COFD[1825] = 5.04367502E+00;
    COFD[1826] = -4.34153325E-01;
    COFD[1827] = 1.85956055E-02;
    COFD[1828] = -1.98418115E+01;
    COFD[1829] = 5.04367502E+00;
    COFD[1830] = -4.34153325E-01;
    COFD[1831] = 1.85956055E-02;
    COFD[1832] = -2.12218960E+01;
    COFD[1833] = 5.39823225E+00;
    COFD[1834] = -4.72294645E-01;
    COFD[1835] = 1.99340225E-02;
    COFD[1836] = -2.18190539E+01;
    COFD[1837] = 5.55753905E+00;
    COFD[1838] = -4.88136714E-01;
    COFD[1839] = 2.04294957E-02;
    COFD[1840] = -1.99081556E+01;
    COFD[1841] = 5.09311649E+00;
    COFD[1842] = -4.39965178E-01;
    COFD[1843] = 1.88238537E-02;
    COFD[1844] = -2.20546151E+01;
    COFD[1845] = 5.44448440E+00;
    COFD[1846] = -4.51529024E-01;
    COFD[1847] = 1.79698119E-02;
    COFD[1848] = -2.20500806E+01;
    COFD[1849] = 5.44448440E+00;
    COFD[1850] = -4.51529024E-01;
    COFD[1851] = 1.79698119E-02;
    COFD[1852] = -2.20600677E+01;
    COFD[1853] = 5.36785880E+00;
    COFD[1854] = -4.37683002E-01;
    COFD[1855] = 1.72144393E-02;
    COFD[1856] = -2.02646611E+01;
    COFD[1857] = 5.10426133E+00;
    COFD[1858] = -4.41256919E-01;
    COFD[1859] = 1.88737290E-02;
    COFD[1860] = -1.94688688E+01;
    COFD[1861] = 5.43830787E+00;
    COFD[1862] = -4.75472880E-01;
    COFD[1863] = 1.99909996E-02;
    COFD[1864] = -1.57034851E+01;
    COFD[1865] = 3.93614244E+00;
    COFD[1866] = -2.99111497E-01;
    COFD[1867] = 1.30888229E-02;
    COFD[1868] = -1.90883268E+01;
    COFD[1869] = 4.84384483E+00;
    COFD[1870] = -4.10265575E-01;
    COFD[1871] = 1.76414287E-02;
    COFD[1872] = -1.91102652E+01;
    COFD[1873] = 4.84384483E+00;
    COFD[1874] = -4.10265575E-01;
    COFD[1875] = 1.76414287E-02;
    COFD[1876] = -2.05184870E+01;
    COFD[1877] = 5.18417470E+00;
    COFD[1878] = -4.49491573E-01;
    COFD[1879] = 1.91438508E-02;
    COFD[1880] = -2.05356023E+01;
    COFD[1881] = 5.18417470E+00;
    COFD[1882] = -4.49491573E-01;
    COFD[1883] = 1.91438508E-02;
    COFD[1884] = -1.87383952E+01;
    COFD[1885] = 3.96926341E+00;
    COFD[1886] = -2.16412264E-01;
    COFD[1887] = 6.06012078E-03;
    COFD[1888] = -2.05272328E+01;
    COFD[1889] = 5.18417470E+00;
    COFD[1890] = -4.49491573E-01;
    COFD[1891] = 1.91438508E-02;
    COFD[1892] = -2.02922701E+01;
    COFD[1893] = 5.11106992E+00;
    COFD[1894] = -4.42047129E-01;
    COFD[1895] = 1.89042990E-02;
    COFD[1896] = -2.11606963E+01;
    COFD[1897] = 5.42846112E+00;
    COFD[1898] = -4.74321870E-01;
    COFD[1899] = 1.99459749E-02;
    COFD[1900] = -2.01015340E+01;
    COFD[1901] = 4.41511629E+00;
    COFD[1902] = -2.84086963E-01;
    COFD[1903] = 9.37586971E-03;
    COFD[1904] = -2.22116706E+01;
    COFD[1905] = 5.54251230E+00;
    COFD[1906] = -4.70946314E-01;
    COFD[1907] = 1.90785869E-02;
    COFD[1908] = -2.11309207E+01;
    COFD[1909] = 5.41773516E+00;
    COFD[1910] = -4.73414338E-01;
    COFD[1911] = 1.99258685E-02;
    COFD[1912] = -2.21065306E+01;
    COFD[1913] = 5.58360799E+00;
    COFD[1914] = -4.82701436E-01;
    COFD[1915] = 1.98437922E-02;
    COFD[1916] = -2.20453723E+01;
    COFD[1917] = 5.44448440E+00;
    COFD[1918] = -4.51529024E-01;
    COFD[1919] = 1.79698119E-02;
    COFD[1920] = -2.09002742E+01;
    COFD[1921] = 4.72895031E+00;
    COFD[1922] = -3.33332771E-01;
    COFD[1923] = 1.18431478E-02;
    COFD[1924] = -2.22159630E+01;
    COFD[1925] = 5.51722375E+00;
    COFD[1926] = -4.66081431E-01;
    COFD[1927] = 1.88044011E-02;
    COFD[1928] = -1.87688110E+01;
    COFD[1929] = 4.71729964E+00;
    COFD[1930] = -3.95432573E-01;
    COFD[1931] = 1.70623691E-02;
    COFD[1932] = -1.90116191E+01;
    COFD[1933] = 4.84384483E+00;
    COFD[1934] = -4.10265575E-01;
    COFD[1935] = 1.76414287E-02;
    COFD[1936] = -2.00915040E+01;
    COFD[1937] = 4.41511629E+00;
    COFD[1938] = -2.84086963E-01;
    COFD[1939] = 9.37586971E-03;
    COFD[1940] = -2.11349086E+01;
    COFD[1941] = 5.42846112E+00;
    COFD[1942] = -4.74321870E-01;
    COFD[1943] = 1.99459749E-02;
    COFD[1944] = -2.11349086E+01;
    COFD[1945] = 5.42846112E+00;
    COFD[1946] = -4.74321870E-01;
    COFD[1947] = 1.99459749E-02;
    COFD[1948] = -2.21181720E+01;
    COFD[1949] = 5.58360799E+00;
    COFD[1950] = -4.82701436E-01;
    COFD[1951] = 1.98437922E-02;
    COFD[1952] = -2.22059540E+01;
    COFD[1953] = 5.51722375E+00;
    COFD[1954] = -4.66081431E-01;
    COFD[1955] = 1.88044011E-02;
    COFD[1956] = -2.12621914E+01;
    COFD[1957] = 5.47935225E+00;
    COFD[1958] = -4.80056796E-01;
    COFD[1959] = 2.01607180E-02;
    COFD[1960] = -2.09118474E+01;
    COFD[1961] = 4.72895031E+00;
    COFD[1962] = -3.33332771E-01;
    COFD[1963] = 1.18431478E-02;
    COFD[1964] = -2.09061629E+01;
    COFD[1965] = 4.72895031E+00;
    COFD[1966] = -3.33332771E-01;
    COFD[1967] = 1.18431478E-02;
    COFD[1968] = -2.07072145E+01;
    COFD[1969] = 4.56211059E+00;
    COFD[1970] = -3.06895158E-01;
    COFD[1971] = 1.05100393E-02;
    COFD[1972] = -1.83249299E+01;
    COFD[1973] = 4.42045763E+00;
    COFD[1974] = -3.59451578E-01;
    COFD[1975] = 1.56056164E-02;
    COFD[1976] = -1.79344949E+01;
    COFD[1977] = 4.91373893E+00;
    COFD[1978] = -4.18747629E-01;
    COFD[1979] = 1.79856610E-02;
    COFD[1980] = -1.39924781E+01;
    COFD[1981] = 3.26384506E+00;
    COFD[1982] = -2.12947087E-01;
    COFD[1983] = 9.39743888E-03;
    COFD[1984] = -1.72556499E+01;
    COFD[1985] = 4.17889917E+00;
    COFD[1986] = -3.29752510E-01;
    COFD[1987] = 1.43850275E-02;
    COFD[1988] = -1.72753760E+01;
    COFD[1989] = 4.17889917E+00;
    COFD[1990] = -3.29752510E-01;
    COFD[1991] = 1.43850275E-02;
    COFD[1992] = -1.86424545E+01;
    COFD[1993] = 4.53572533E+00;
    COFD[1994] = -3.73386925E-01;
    COFD[1995] = 1.61678881E-02;
    COFD[1996] = -1.86570209E+01;
    COFD[1997] = 4.53572533E+00;
    COFD[1998] = -3.73386925E-01;
    COFD[1999] = 1.61678881E-02;
    COFD[2000] = -2.12661865E+01;
    COFD[2001] = 5.24930667E+00;
    COFD[2002] = -4.17435088E-01;
    COFD[2003] = 1.61434424E-02;
    COFD[2004] = -1.86499071E+01;
    COFD[2005] = 4.53572533E+00;
    COFD[2006] = -3.73386925E-01;
    COFD[2007] = 1.61678881E-02;
    COFD[2008] = -1.83538377E+01;
    COFD[2009] = 4.42828044E+00;
    COFD[2010] = -3.60417833E-01;
    COFD[2011] = 1.56455103E-02;
    COFD[2012] = -1.95552142E+01;
    COFD[2013] = 4.90255048E+00;
    COFD[2014] = -4.17368501E-01;
    COFD[2015] = 1.79287358E-02;
    COFD[2016] = -2.19399793E+01;
    COFD[2017] = 5.41841631E+00;
    COFD[2018] = -4.46818971E-01;
    COFD[2019] = 1.77127652E-02;
    COFD[2020] = -2.14326461E+01;
    COFD[2021] = 5.41729961E+00;
    COFD[2022] = -4.73400281E-01;
    COFD[2023] = 1.99269567E-02;
    COFD[2024] = -1.94823660E+01;
    COFD[2025] = 4.87333294E+00;
    COFD[2026] = -4.13769241E-01;
    COFD[2027] = 1.77802244E-02;
    COFD[2028] = -2.10944088E+01;
    COFD[2029] = 5.34286099E+00;
    COFD[2030] = -4.68100992E-01;
    COFD[2031] = 1.98731399E-02;
    COFD[2032] = -2.18273547E+01;
    COFD[2033] = 5.55753905E+00;
    COFD[2034] = -4.88136714E-01;
    COFD[2035] = 2.04294957E-02;
    COFD[2036] = -2.22159630E+01;
    COFD[2037] = 5.51722375E+00;
    COFD[2038] = -4.66081431E-01;
    COFD[2039] = 1.88044011E-02;
    COFD[2040] = -2.15746136E+01;
    COFD[2041] = 5.44803850E+00;
    COFD[2042] = -4.76610560E-01;
    COFD[2043] = 2.00355294E-02;
    COFD[2044] = -1.69540369E+01;
    COFD[2045] = 4.05099737E+00;
    COFD[2046] = -3.13841660E-01;
    COFD[2047] = 1.37218854E-02;
    COFD[2048] = -1.71860230E+01;
    COFD[2049] = 4.17889917E+00;
    COFD[2050] = -3.29752510E-01;
    COFD[2051] = 1.43850275E-02;
    COFD[2052] = -2.19313638E+01;
    COFD[2053] = 5.41841631E+00;
    COFD[2054] = -4.46818971E-01;
    COFD[2055] = 1.77127652E-02;
    COFD[2056] = -1.95318173E+01;
    COFD[2057] = 4.90255048E+00;
    COFD[2058] = -4.17368501E-01;
    COFD[2059] = 1.79287358E-02;
    COFD[2060] = -1.95318173E+01;
    COFD[2061] = 4.90255048E+00;
    COFD[2062] = -4.17368501E-01;
    COFD[2063] = 1.79287358E-02;
    COFD[2064] = -2.11044965E+01;
    COFD[2065] = 5.34286099E+00;
    COFD[2066] = -4.68100992E-01;
    COFD[2067] = 1.98731399E-02;
    COFD[2068] = -2.15660171E+01;
    COFD[2069] = 5.44803850E+00;
    COFD[2070] = -4.76610560E-01;
    COFD[2071] = 2.00355294E-02;
    COFD[2072] = -1.96408127E+01;
    COFD[2073] = 4.95923807E+00;
    COFD[2074] = -4.24176182E-01;
    COFD[2075] = 1.82020215E-02;
    COFD[2076] = -2.22255968E+01;
    COFD[2077] = 5.51722375E+00;
    COFD[2078] = -4.66081431E-01;
    COFD[2079] = 1.88044011E-02;
    COFD[2080] = -2.22208695E+01;
    COFD[2081] = 5.51722375E+00;
    COFD[2082] = -4.66081431E-01;
    COFD[2083] = 1.88044011E-02;
    COFD[2084] = -2.23075353E+01;
    COFD[2085] = 5.47511492E+00;
    COFD[2086] = -4.57098839E-01;
    COFD[2087] = 1.82750523E-02;
    COFD[2088] = -1.38661480E+01;
    COFD[2089] = 2.97137588E+00;
    COFD[2090] = -1.75491257E-01;
    COFD[2091] = 7.79646773E-03;
    COFD[2092] = -1.32467319E+01;
    COFD[2093] = 3.34156587E+00;
    COFD[2094] = -2.22853306E-01;
    COFD[2095] = 9.81883417E-03;
    COFD[2096] = -1.08369481E+01;
    COFD[2097] = 2.19094415E+00;
    COFD[2098] = -7.11992510E-02;
    COFD[2099] = 3.14105973E-03;
    COFD[2100] = -1.30610083E+01;
    COFD[2101] = 2.80913567E+00;
    COFD[2102] = -1.54536855E-01;
    COFD[2103] = 6.89359313E-03;
    COFD[2104] = -1.30738796E+01;
    COFD[2105] = 2.80913567E+00;
    COFD[2106] = -1.54536855E-01;
    COFD[2107] = 6.89359313E-03;
    COFD[2108] = -1.40707473E+01;
    COFD[2109] = 3.05837263E+00;
    COFD[2110] = -1.86672802E-01;
    COFD[2111] = 8.27575734E-03;
    COFD[2112] = -1.40789009E+01;
    COFD[2113] = 3.05837263E+00;
    COFD[2114] = -1.86672802E-01;
    COFD[2115] = 8.27575734E-03;
    COFD[2116] = -1.92003817E+01;
    COFD[2117] = 5.05708283E+00;
    COFD[2118] = -4.35739290E-01;
    COFD[2119] = 1.86583205E-02;
    COFD[2120] = -1.40749320E+01;
    COFD[2121] = 3.05837263E+00;
    COFD[2122] = -1.86672802E-01;
    COFD[2123] = 8.27575734E-03;
    COFD[2124] = -1.38862839E+01;
    COFD[2125] = 2.97564184E+00;
    COFD[2126] = -1.76025309E-01;
    COFD[2127] = 7.81869993E-03;
    COFD[2128] = -1.47715918E+01;
    COFD[2129] = 3.33113524E+00;
    COFD[2130] = -2.21479057E-01;
    COFD[2131] = 9.75837737E-03;
    COFD[2132] = -1.90576320E+01;
    COFD[2133] = 4.86821670E+00;
    COFD[2134] = -4.13144121E-01;
    COFD[2135] = 1.77546701E-02;
    COFD[2136] = -1.67390425E+01;
    COFD[2137] = 4.00828594E+00;
    COFD[2138] = -3.08414344E-01;
    COFD[2139] = 1.34907430E-02;
    COFD[2140] = -1.47048024E+01;
    COFD[2141] = 3.30594991E+00;
    COFD[2142] = -2.18182207E-01;
    COFD[2143] = 9.61429447E-03;
    COFD[2144] = -1.61101966E+01;
    COFD[2145] = 3.75910622E+00;
    COFD[2146] = -2.75986578E-01;
    COFD[2147] = 1.20782843E-02;
    COFD[2148] = -1.71884218E+01;
    COFD[2149] = 4.17190426E+00;
    COFD[2150] = -3.28894681E-01;
    COFD[2151] = 1.43498101E-02;
    COFD[2152] = -1.87688110E+01;
    COFD[2153] = 4.71729964E+00;
    COFD[2154] = -3.95432573E-01;
    COFD[2155] = 1.70623691E-02;
    COFD[2156] = -1.69540369E+01;
    COFD[2157] = 4.05099737E+00;
    COFD[2158] = -3.13841660E-01;
    COFD[2159] = 1.37218854E-02;
    COFD[2160] = -1.29573076E+01;
    COFD[2161] = 2.73155251E+00;
    COFD[2162] = -1.44594082E-01;
    COFD[2163] = 6.46883252E-03;
    COFD[2164] = -1.30141899E+01;
    COFD[2165] = 2.80913567E+00;
    COFD[2166] = -1.54536855E-01;
    COFD[2167] = 6.89359313E-03;
    COFD[2168] = -1.90526941E+01;
    COFD[2169] = 4.86821670E+00;
    COFD[2170] = -4.13144121E-01;
    COFD[2171] = 1.77546701E-02;
    COFD[2172] = -1.47558850E+01;
    COFD[2173] = 3.33113524E+00;
    COFD[2174] = -2.21479057E-01;
    COFD[2175] = 9.75837737E-03;
    COFD[2176] = -1.47558850E+01;
    COFD[2177] = 3.33113524E+00;
    COFD[2178] = -2.21479057E-01;
    COFD[2179] = 9.75837737E-03;
    COFD[2180] = -1.61161138E+01;
    COFD[2181] = 3.75910622E+00;
    COFD[2182] = -2.75986578E-01;
    COFD[2183] = 1.20782843E-02;
    COFD[2184] = -1.69491115E+01;
    COFD[2185] = 4.05099737E+00;
    COFD[2186] = -3.13841660E-01;
    COFD[2187] = 1.37218854E-02;
    COFD[2188] = -1.46907028E+01;
    COFD[2189] = 3.39229020E+00;
    COFD[2190] = -2.29520232E-01;
    COFD[2191] = 1.01114311E-02;
    COFD[2192] = -1.87739217E+01;
    COFD[2193] = 4.71729964E+00;
    COFD[2194] = -3.95432573E-01;
    COFD[2195] = 1.70623691E-02;
    COFD[2196] = -1.87714196E+01;
    COFD[2197] = 4.71729964E+00;
    COFD[2198] = -3.95432573E-01;
    COFD[2199] = 1.70623691E-02;
    COFD[2200] = -1.91468116E+01;
    COFD[2201] = 4.80472875E+00;
    COFD[2202] = -4.05752803E-01;
    COFD[2203] = 1.74685236E-02;
    COFD[2204] = -1.40076852E+01;
    COFD[2205] = 3.07549274E+00;
    COFD[2206] = -1.88889344E-01;
    COFD[2207] = 8.37152866E-03;
    COFD[2208] = -1.34162893E+01;
    COFD[2209] = 3.48624238E+00;
    COFD[2210] = -2.41554467E-01;
    COFD[2211] = 1.06263545E-02;
    COFD[2212] = -1.09469245E+01;
    COFD[2213] = 2.30836460E+00;
    COFD[2214] = -8.76339315E-02;
    COFD[2215] = 3.90878445E-03;
    COFD[2216] = -1.31551788E+01;
    COFD[2217] = 2.90778936E+00;
    COFD[2218] = -1.67388544E-01;
    COFD[2219] = 7.45220609E-03;
    COFD[2220] = -1.31686537E+01;
    COFD[2221] = 2.90778936E+00;
    COFD[2222] = -1.67388544E-01;
    COFD[2223] = 7.45220609E-03;
    COFD[2224] = -1.42429085E+01;
    COFD[2225] = 3.17651319E+00;
    COFD[2226] = -2.02028974E-01;
    COFD[2227] = 8.94232502E-03;
    COFD[2228] = -1.42515527E+01;
    COFD[2229] = 3.17651319E+00;
    COFD[2230] = -2.02028974E-01;
    COFD[2231] = 8.94232502E-03;
    COFD[2232] = -1.93521390E+01;
    COFD[2233] = 5.16013126E+00;
    COFD[2234] = -4.46824543E-01;
    COFD[2235] = 1.90464887E-02;
    COFD[2236] = -1.42473439E+01;
    COFD[2237] = 3.17651319E+00;
    COFD[2238] = -2.02028974E-01;
    COFD[2239] = 8.94232502E-03;
    COFD[2240] = -1.40318948E+01;
    COFD[2241] = 3.08120012E+00;
    COFD[2242] = -1.89629903E-01;
    COFD[2243] = 8.40361952E-03;
    COFD[2244] = -1.50240272E+01;
    COFD[2245] = 3.47945612E+00;
    COFD[2246] = -2.40703722E-01;
    COFD[2247] = 1.05907441E-02;
    COFD[2248] = -1.93677186E+01;
    COFD[2249] = 5.02567894E+00;
    COFD[2250] = -4.32045169E-01;
    COFD[2251] = 1.85132214E-02;
    COFD[2252] = -1.69758891E+01;
    COFD[2253] = 4.14240922E+00;
    COFD[2254] = -3.25239774E-01;
    COFD[2255] = 1.41980687E-02;
    COFD[2256] = -1.49727799E+01;
    COFD[2257] = 3.46140064E+00;
    COFD[2258] = -2.38440092E-01;
    COFD[2259] = 1.04960087E-02;
    COFD[2260] = -1.64825368E+01;
    COFD[2261] = 3.95035840E+00;
    COFD[2262] = -3.00959418E-01;
    COFD[2263] = 1.31692593E-02;
    COFD[2264] = -1.74111692E+01;
    COFD[2265] = 4.29676909E+00;
    COFD[2266] = -3.44085306E-01;
    COFD[2267] = 1.49671135E-02;
    COFD[2268] = -1.90116191E+01;
    COFD[2269] = 4.84384483E+00;
    COFD[2270] = -4.10265575E-01;
    COFD[2271] = 1.76414287E-02;
    COFD[2272] = -1.71860230E+01;
    COFD[2273] = 4.17889917E+00;
    COFD[2274] = -3.29752510E-01;
    COFD[2275] = 1.43850275E-02;
    COFD[2276] = -1.30141899E+01;
    COFD[2277] = 2.80913567E+00;
    COFD[2278] = -1.54536855E-01;
    COFD[2279] = 6.89359313E-03;
    COFD[2280] = -1.31062967E+01;
    COFD[2281] = 2.90778936E+00;
    COFD[2282] = -1.67388544E-01;
    COFD[2283] = 7.45220609E-03;
    COFD[2284] = -1.93624931E+01;
    COFD[2285] = 5.02567894E+00;
    COFD[2286] = -4.32045169E-01;
    COFD[2287] = 1.85132214E-02;
    COFD[2288] = -1.50076254E+01;
    COFD[2289] = 3.47945612E+00;
    COFD[2290] = -2.40703722E-01;
    COFD[2291] = 1.05907441E-02;
    COFD[2292] = -1.50076254E+01;
    COFD[2293] = 3.47945612E+00;
    COFD[2294] = -2.40703722E-01;
    COFD[2295] = 1.05907441E-02;
    COFD[2296] = -1.64887871E+01;
    COFD[2297] = 3.95035840E+00;
    COFD[2298] = -3.00959418E-01;
    COFD[2299] = 1.31692593E-02;
    COFD[2300] = -1.71808106E+01;
    COFD[2301] = 4.17889917E+00;
    COFD[2302] = -3.29752510E-01;
    COFD[2303] = 1.43850275E-02;
    COFD[2304] = -1.48738066E+01;
    COFD[2305] = 3.52327209E+00;
    COFD[2306] = -2.46286208E-01;
    COFD[2307] = 1.08285963E-02;
    COFD[2308] = -1.90170590E+01;
    COFD[2309] = 4.84384483E+00;
    COFD[2310] = -4.10265575E-01;
    COFD[2311] = 1.76414287E-02;
    COFD[2312] = -1.90143953E+01;
    COFD[2313] = 4.84384483E+00;
    COFD[2314] = -4.10265575E-01;
    COFD[2315] = 1.76414287E-02;
    COFD[2316] = -1.94534227E+01;
    COFD[2317] = 4.95249173E+00;
    COFD[2318] = -4.23376552E-01;
    COFD[2319] = 1.81703714E-02;
    COFD[2320] = -2.04750581E+01;
    COFD[2321] = 5.23112374E+00;
    COFD[2322] = -4.54967682E-01;
    COFD[2323] = 1.93570423E-02;
    COFD[2324] = -1.97544450E+01;
    COFD[2325] = 5.56931926E+00;
    COFD[2326] = -4.89105511E-01;
    COFD[2327] = 2.04493129E-02;
    COFD[2328] = -1.60517370E+01;
    COFD[2329] = 4.11188603E+00;
    COFD[2330] = -3.21540884E-01;
    COFD[2331] = 1.40482564E-02;
    COFD[2332] = -1.94313116E+01;
    COFD[2333] = 5.02567894E+00;
    COFD[2334] = -4.32045169E-01;
    COFD[2335] = 1.85132214E-02;
    COFD[2336] = -1.94507876E+01;
    COFD[2337] = 5.02567894E+00;
    COFD[2338] = -4.32045169E-01;
    COFD[2339] = 1.85132214E-02;
    COFD[2340] = -2.08204449E+01;
    COFD[2341] = 5.35267674E+00;
    COFD[2342] = -4.69010505E-01;
    COFD[2343] = 1.98979152E-02;
    COFD[2344] = -2.08347403E+01;
    COFD[2345] = 5.35267674E+00;
    COFD[2346] = -4.69010505E-01;
    COFD[2347] = 1.98979152E-02;
    COFD[2348] = -1.77498543E+01;
    COFD[2349] = 3.57475686E+00;
    COFD[2350] = -1.56396297E-01;
    COFD[2351] = 3.12157721E-03;
    COFD[2352] = -2.08277598E+01;
    COFD[2353] = 5.35267674E+00;
    COFD[2354] = -4.69010505E-01;
    COFD[2355] = 1.98979152E-02;
    COFD[2356] = -2.05045578E+01;
    COFD[2357] = 5.23843909E+00;
    COFD[2358] = -4.55815614E-01;
    COFD[2359] = 1.93898040E-02;
    COFD[2360] = -2.14391943E+01;
    COFD[2361] = 5.56531152E+00;
    COFD[2362] = -4.88789821E-01;
    COFD[2363] = 2.04437116E-02;
    COFD[2364] = -1.90413348E+01;
    COFD[2365] = 3.99221757E+00;
    COFD[2366] = -2.19854880E-01;
    COFD[2367] = 6.22736279E-03;
    COFD[2368] = -2.19215555E+01;
    COFD[2369] = 5.45216133E+00;
    COFD[2370] = -4.52916925E-01;
    COFD[2371] = 1.80456400E-02;
    COFD[2372] = -2.14022336E+01;
    COFD[2373] = 5.55346617E+00;
    COFD[2374] = -4.87783156E-01;
    COFD[2375] = 2.04210886E-02;
    COFD[2376] = -2.19956352E+01;
    COFD[2377] = 5.55935694E+00;
    COFD[2378] = -4.74154740E-01;
    COFD[2379] = 1.92584304E-02;
    COFD[2380] = -2.16296373E+01;
    COFD[2381] = 5.29019717E+00;
    COFD[2382] = -4.24502606E-01;
    COFD[2383] = 1.65197343E-02;
    COFD[2384] = -2.00915040E+01;
    COFD[2385] = 4.41511629E+00;
    COFD[2386] = -2.84086963E-01;
    COFD[2387] = 9.37586971E-03;
    COFD[2388] = -2.19313638E+01;
    COFD[2389] = 5.41841631E+00;
    COFD[2390] = -4.46818971E-01;
    COFD[2391] = 1.77127652E-02;
    COFD[2392] = -1.90526941E+01;
    COFD[2393] = 4.86821670E+00;
    COFD[2394] = -4.13144121E-01;
    COFD[2395] = 1.77546701E-02;
    COFD[2396] = -1.93624931E+01;
    COFD[2397] = 5.02567894E+00;
    COFD[2398] = -4.32045169E-01;
    COFD[2399] = 1.85132214E-02;
    COFD[2400] = -1.90328712E+01;
    COFD[2401] = 3.99221757E+00;
    COFD[2402] = -2.19854880E-01;
    COFD[2403] = 6.22736279E-03;
    COFD[2404] = -2.14160703E+01;
    COFD[2405] = 5.56531152E+00;
    COFD[2406] = -4.88789821E-01;
    COFD[2407] = 2.04437116E-02;
    COFD[2408] = -2.14160703E+01;
    COFD[2409] = 5.56531152E+00;
    COFD[2410] = -4.88789821E-01;
    COFD[2411] = 2.04437116E-02;
    COFD[2412] = -2.20055544E+01;
    COFD[2413] = 5.55935694E+00;
    COFD[2414] = -4.74154740E-01;
    COFD[2415] = 1.92584304E-02;
    COFD[2416] = -2.19229190E+01;
    COFD[2417] = 5.41841631E+00;
    COFD[2418] = -4.46818971E-01;
    COFD[2419] = 1.77127652E-02;
    COFD[2420] = -2.14204185E+01;
    COFD[2421] = 5.59268435E+00;
    COFD[2422] = -4.91232974E-01;
    COFD[2423] = 2.05064746E-02;
    COFD[2424] = -2.01009366E+01;
    COFD[2425] = 4.41511629E+00;
    COFD[2426] = -2.84086963E-01;
    COFD[2427] = 9.37586971E-03;
    COFD[2428] = -2.00963085E+01;
    COFD[2429] = 4.41511629E+00;
    COFD[2430] = -2.84086963E-01;
    COFD[2431] = 9.37586971E-03;
    COFD[2432] = -1.97545910E+01;
    COFD[2433] = 4.18758010E+00;
    COFD[2434] = -2.49327776E-01;
    COFD[2435] = 7.66559103E-03;
    COFD[2436] = -1.59404882E+01;
    COFD[2437] = 3.66853818E+00;
    COFD[2438] = -2.64346221E-01;
    COFD[2439] = 1.15784613E-02;
    COFD[2440] = -1.57972369E+01;
    COFD[2441] = 4.22225052E+00;
    COFD[2442] = -3.35156428E-01;
    COFD[2443] = 1.46104855E-02;
    COFD[2444] = -1.25098960E+01;
    COFD[2445] = 2.77873601E+00;
    COFD[2446] = -1.50637360E-01;
    COFD[2447] = 6.72684281E-03;
    COFD[2448] = -1.50584249E+01;
    COFD[2449] = 3.47945612E+00;
    COFD[2450] = -2.40703722E-01;
    COFD[2451] = 1.05907441E-02;
    COFD[2452] = -1.50724636E+01;
    COFD[2453] = 3.47945612E+00;
    COFD[2454] = -2.40703722E-01;
    COFD[2455] = 1.05907441E-02;
    COFD[2456] = -1.63254691E+01;
    COFD[2457] = 3.82388595E+00;
    COFD[2458] = -2.84480724E-01;
    COFD[2459] = 1.24506311E-02;
    COFD[2460] = -1.63345829E+01;
    COFD[2461] = 3.82388595E+00;
    COFD[2462] = -2.84480724E-01;
    COFD[2463] = 1.24506311E-02;
    COFD[2464] = -2.12639214E+01;
    COFD[2465] = 5.61184117E+00;
    COFD[2466] = -4.90532156E-01;
    COFD[2467] = 2.03507922E-02;
    COFD[2468] = -1.63301444E+01;
    COFD[2469] = 3.82388595E+00;
    COFD[2470] = -2.84480724E-01;
    COFD[2471] = 1.24506311E-02;
    COFD[2472] = -1.59634533E+01;
    COFD[2473] = 3.67388294E+00;
    COFD[2474] = -2.64990709E-01;
    COFD[2475] = 1.16042706E-02;
    COFD[2476] = -1.73198034E+01;
    COFD[2477] = 4.21416723E+00;
    COFD[2478] = -3.34163932E-01;
    COFD[2479] = 1.45697432E-02;
    COFD[2480] = -2.14215700E+01;
    COFD[2481] = 5.56531152E+00;
    COFD[2482] = -4.88789821E-01;
    COFD[2483] = 2.04437116E-02;
    COFD[2484] = -1.93015555E+01;
    COFD[2485] = 4.85015581E+00;
    COFD[2486] = -4.10945109E-01;
    COFD[2487] = 1.76651398E-02;
    COFD[2488] = -1.72556729E+01;
    COFD[2489] = 4.19029808E+00;
    COFD[2490] = -3.31177076E-01;
    COFD[2491] = 1.44446234E-02;
    COFD[2492] = -1.88241079E+01;
    COFD[2493] = 4.68393046E+00;
    COFD[2494] = -3.91610863E-01;
    COFD[2495] = 1.69174645E-02;
    COFD[2496] = -1.98418115E+01;
    COFD[2497] = 5.04367502E+00;
    COFD[2498] = -4.34153325E-01;
    COFD[2499] = 1.85956055E-02;
    COFD[2500] = -2.11349086E+01;
    COFD[2501] = 5.42846112E+00;
    COFD[2502] = -4.74321870E-01;
    COFD[2503] = 1.99459749E-02;
    COFD[2504] = -1.95318173E+01;
    COFD[2505] = 4.90255048E+00;
    COFD[2506] = -4.17368501E-01;
    COFD[2507] = 1.79287358E-02;
    COFD[2508] = -1.47558850E+01;
    COFD[2509] = 3.33113524E+00;
    COFD[2510] = -2.21479057E-01;
    COFD[2511] = 9.75837737E-03;
    COFD[2512] = -1.50076254E+01;
    COFD[2513] = 3.47945612E+00;
    COFD[2514] = -2.40703722E-01;
    COFD[2515] = 1.05907441E-02;
    COFD[2516] = -2.14160703E+01;
    COFD[2517] = 5.56531152E+00;
    COFD[2518] = -4.88789821E-01;
    COFD[2519] = 2.04437116E-02;
    COFD[2520] = -1.73027557E+01;
    COFD[2521] = 4.21416723E+00;
    COFD[2522] = -3.34163932E-01;
    COFD[2523] = 1.45697432E-02;
    COFD[2524] = -1.73027557E+01;
    COFD[2525] = 4.21416723E+00;
    COFD[2526] = -3.34163932E-01;
    COFD[2527] = 1.45697432E-02;
    COFD[2528] = -1.88306747E+01;
    COFD[2529] = 4.68393046E+00;
    COFD[2530] = -3.91610863E-01;
    COFD[2531] = 1.69174645E-02;
    COFD[2532] = -1.95263312E+01;
    COFD[2533] = 4.90255048E+00;
    COFD[2534] = -4.17368501E-01;
    COFD[2535] = 1.79287358E-02;
    COFD[2536] = -1.72572042E+01;
    COFD[2537] = 4.26063341E+00;
    COFD[2538] = -3.39848064E-01;
    COFD[2539] = 1.48021313E-02;
    COFD[2540] = -2.11406662E+01;
    COFD[2541] = 5.42846112E+00;
    COFD[2542] = -4.74321870E-01;
    COFD[2543] = 1.99459749E-02;
    COFD[2544] = -2.11378465E+01;
    COFD[2545] = 5.42846112E+00;
    COFD[2546] = -4.74321870E-01;
    COFD[2547] = 1.99459749E-02;
    COFD[2548] = -2.15326361E+01;
    COFD[2549] = 5.51982454E+00;
    COFD[2550] = -4.84452039E-01;
    COFD[2551] = 2.03175522E-02;
    COFD[2552] = -1.59404882E+01;
    COFD[2553] = 3.66853818E+00;
    COFD[2554] = -2.64346221E-01;
    COFD[2555] = 1.15784613E-02;
    COFD[2556] = -1.57972369E+01;
    COFD[2557] = 4.22225052E+00;
    COFD[2558] = -3.35156428E-01;
    COFD[2559] = 1.46104855E-02;
    COFD[2560] = -1.25098960E+01;
    COFD[2561] = 2.77873601E+00;
    COFD[2562] = -1.50637360E-01;
    COFD[2563] = 6.72684281E-03;
    COFD[2564] = -1.50584249E+01;
    COFD[2565] = 3.47945612E+00;
    COFD[2566] = -2.40703722E-01;
    COFD[2567] = 1.05907441E-02;
    COFD[2568] = -1.50724636E+01;
    COFD[2569] = 3.47945612E+00;
    COFD[2570] = -2.40703722E-01;
    COFD[2571] = 1.05907441E-02;
    COFD[2572] = -1.63254691E+01;
    COFD[2573] = 3.82388595E+00;
    COFD[2574] = -2.84480724E-01;
    COFD[2575] = 1.24506311E-02;
    COFD[2576] = -1.63345829E+01;
    COFD[2577] = 3.82388595E+00;
    COFD[2578] = -2.84480724E-01;
    COFD[2579] = 1.24506311E-02;
    COFD[2580] = -2.12639214E+01;
    COFD[2581] = 5.61184117E+00;
    COFD[2582] = -4.90532156E-01;
    COFD[2583] = 2.03507922E-02;
    COFD[2584] = -1.63301444E+01;
    COFD[2585] = 3.82388595E+00;
    COFD[2586] = -2.84480724E-01;
    COFD[2587] = 1.24506311E-02;
    COFD[2588] = -1.59634533E+01;
    COFD[2589] = 3.67388294E+00;
    COFD[2590] = -2.64990709E-01;
    COFD[2591] = 1.16042706E-02;
    COFD[2592] = -1.73198034E+01;
    COFD[2593] = 4.21416723E+00;
    COFD[2594] = -3.34163932E-01;
    COFD[2595] = 1.45697432E-02;
    COFD[2596] = -2.14215700E+01;
    COFD[2597] = 5.56531152E+00;
    COFD[2598] = -4.88789821E-01;
    COFD[2599] = 2.04437116E-02;
    COFD[2600] = -1.93015555E+01;
    COFD[2601] = 4.85015581E+00;
    COFD[2602] = -4.10945109E-01;
    COFD[2603] = 1.76651398E-02;
    COFD[2604] = -1.72556729E+01;
    COFD[2605] = 4.19029808E+00;
    COFD[2606] = -3.31177076E-01;
    COFD[2607] = 1.44446234E-02;
    COFD[2608] = -1.88241079E+01;
    COFD[2609] = 4.68393046E+00;
    COFD[2610] = -3.91610863E-01;
    COFD[2611] = 1.69174645E-02;
    COFD[2612] = -1.98418115E+01;
    COFD[2613] = 5.04367502E+00;
    COFD[2614] = -4.34153325E-01;
    COFD[2615] = 1.85956055E-02;
    COFD[2616] = -2.11349086E+01;
    COFD[2617] = 5.42846112E+00;
    COFD[2618] = -4.74321870E-01;
    COFD[2619] = 1.99459749E-02;
    COFD[2620] = -1.95318173E+01;
    COFD[2621] = 4.90255048E+00;
    COFD[2622] = -4.17368501E-01;
    COFD[2623] = 1.79287358E-02;
    COFD[2624] = -1.47558850E+01;
    COFD[2625] = 3.33113524E+00;
    COFD[2626] = -2.21479057E-01;
    COFD[2627] = 9.75837737E-03;
    COFD[2628] = -1.50076254E+01;
    COFD[2629] = 3.47945612E+00;
    COFD[2630] = -2.40703722E-01;
    COFD[2631] = 1.05907441E-02;
    COFD[2632] = -2.14160703E+01;
    COFD[2633] = 5.56531152E+00;
    COFD[2634] = -4.88789821E-01;
    COFD[2635] = 2.04437116E-02;
    COFD[2636] = -1.73027557E+01;
    COFD[2637] = 4.21416723E+00;
    COFD[2638] = -3.34163932E-01;
    COFD[2639] = 1.45697432E-02;
    COFD[2640] = -1.73027557E+01;
    COFD[2641] = 4.21416723E+00;
    COFD[2642] = -3.34163932E-01;
    COFD[2643] = 1.45697432E-02;
    COFD[2644] = -1.88306747E+01;
    COFD[2645] = 4.68393046E+00;
    COFD[2646] = -3.91610863E-01;
    COFD[2647] = 1.69174645E-02;
    COFD[2648] = -1.95263312E+01;
    COFD[2649] = 4.90255048E+00;
    COFD[2650] = -4.17368501E-01;
    COFD[2651] = 1.79287358E-02;
    COFD[2652] = -1.72572042E+01;
    COFD[2653] = 4.26063341E+00;
    COFD[2654] = -3.39848064E-01;
    COFD[2655] = 1.48021313E-02;
    COFD[2656] = -2.11406662E+01;
    COFD[2657] = 5.42846112E+00;
    COFD[2658] = -4.74321870E-01;
    COFD[2659] = 1.99459749E-02;
    COFD[2660] = -2.11378465E+01;
    COFD[2661] = 5.42846112E+00;
    COFD[2662] = -4.74321870E-01;
    COFD[2663] = 1.99459749E-02;
    COFD[2664] = -2.15326361E+01;
    COFD[2665] = 5.51982454E+00;
    COFD[2666] = -4.84452039E-01;
    COFD[2667] = 2.03175522E-02;
    COFD[2668] = -1.76099552E+01;
    COFD[2669] = 4.19171952E+00;
    COFD[2670] = -3.31354810E-01;
    COFD[2671] = 1.44520623E-02;
    COFD[2672] = -1.72239172E+01;
    COFD[2673] = 4.69060745E+00;
    COFD[2674] = -3.92369888E-01;
    COFD[2675] = 1.69459661E-02;
    COFD[2676] = -1.34723215E+01;
    COFD[2677] = 3.09379603E+00;
    COFD[2678] = -1.91268635E-01;
    COFD[2679] = 8.47480224E-03;
    COFD[2680] = -1.65559787E+01;
    COFD[2681] = 3.95035840E+00;
    COFD[2682] = -3.00959418E-01;
    COFD[2683] = 1.31692593E-02;
    COFD[2684] = -1.65749533E+01;
    COFD[2685] = 3.95035840E+00;
    COFD[2686] = -3.00959418E-01;
    COFD[2687] = 1.31692593E-02;
    COFD[2688] = -1.78938745E+01;
    COFD[2689] = 4.29613154E+00;
    COFD[2690] = -3.44012526E-01;
    COFD[2691] = 1.49643715E-02;
    COFD[2692] = -1.79076361E+01;
    COFD[2693] = 4.29613154E+00;
    COFD[2694] = -3.44012526E-01;
    COFD[2695] = 1.49643715E-02;
    COFD[2696] = -2.15172770E+01;
    COFD[2697] = 5.46737673E+00;
    COFD[2698] = -4.55696085E-01;
    COFD[2699] = 1.81982625E-02;
    COFD[2700] = -1.79009181E+01;
    COFD[2701] = 4.29613154E+00;
    COFD[2702] = -3.44012526E-01;
    COFD[2703] = 1.49643715E-02;
    COFD[2704] = -1.76383156E+01;
    COFD[2705] = 4.19935698E+00;
    COFD[2706] = -3.32310212E-01;
    COFD[2707] = 1.44920670E-02;
    COFD[2708] = -1.88532497E+01;
    COFD[2709] = 4.68393046E+00;
    COFD[2710] = -3.91610863E-01;
    COFD[2711] = 1.69174645E-02;
    COFD[2712] = -2.20137178E+01;
    COFD[2713] = 5.55935694E+00;
    COFD[2714] = -4.74154740E-01;
    COFD[2715] = 1.92584304E-02;
    COFD[2716] = -2.09309753E+01;
    COFD[2717] = 5.30153901E+00;
    COFD[2718] = -4.63335119E-01;
    COFD[2719] = 1.96897053E-02;
    COFD[2720] = -1.87968848E+01;
    COFD[2721] = 4.66162351E+00;
    COFD[2722] = -3.88920477E-01;
    COFD[2723] = 1.68089648E-02;
    COFD[2724] = -2.03861000E+01;
    COFD[2725] = 5.13263469E+00;
    COFD[2726] = -4.44457285E-01;
    COFD[2727] = 1.89932102E-02;
    COFD[2728] = -2.12218960E+01;
    COFD[2729] = 5.39823225E+00;
    COFD[2730] = -4.72294645E-01;
    COFD[2731] = 1.99340225E-02;
    COFD[2732] = -2.21181720E+01;
    COFD[2733] = 5.58360799E+00;
    COFD[2734] = -4.82701436E-01;
    COFD[2735] = 1.98437922E-02;
    COFD[2736] = -2.11044965E+01;
    COFD[2737] = 5.34286099E+00;
    COFD[2738] = -4.68100992E-01;
    COFD[2739] = 1.98731399E-02;
    COFD[2740] = -1.61161138E+01;
    COFD[2741] = 3.75910622E+00;
    COFD[2742] = -2.75986578E-01;
    COFD[2743] = 1.20782843E-02;
    COFD[2744] = -1.64887871E+01;
    COFD[2745] = 3.95035840E+00;
    COFD[2746] = -3.00959418E-01;
    COFD[2747] = 1.31692593E-02;
    COFD[2748] = -2.20055544E+01;
    COFD[2749] = 5.55935694E+00;
    COFD[2750] = -4.74154740E-01;
    COFD[2751] = 1.92584304E-02;
    COFD[2752] = -1.88306747E+01;
    COFD[2753] = 4.68393046E+00;
    COFD[2754] = -3.91610863E-01;
    COFD[2755] = 1.69174645E-02;
    COFD[2756] = -1.88306747E+01;
    COFD[2757] = 4.68393046E+00;
    COFD[2758] = -3.91610863E-01;
    COFD[2759] = 1.69174645E-02;
    COFD[2760] = -2.03956853E+01;
    COFD[2761] = 5.13263469E+00;
    COFD[2762] = -4.44457285E-01;
    COFD[2763] = 1.89932102E-02;
    COFD[2764] = -2.10963515E+01;
    COFD[2765] = 5.34286099E+00;
    COFD[2766] = -4.68100992E-01;
    COFD[2767] = 1.98731399E-02;
    COFD[2768] = -1.88761387E+01;
    COFD[2769] = 4.72476764E+00;
    COFD[2770] = -3.96306836E-01;
    COFD[2771] = 1.70964541E-02;
    COFD[2772] = -2.21272109E+01;
    COFD[2773] = 5.58360799E+00;
    COFD[2774] = -4.82701436E-01;
    COFD[2775] = 1.98437922E-02;
    COFD[2776] = -2.21227768E+01;
    COFD[2777] = 5.58360799E+00;
    COFD[2778] = -4.82701436E-01;
    COFD[2779] = 1.98437922E-02;
    COFD[2780] = -2.23101989E+01;
    COFD[2781] = 5.58471203E+00;
    COFD[2782] = -4.79905311E-01;
    COFD[2783] = 1.96058913E-02;
    COFD[2784] = -1.83166353E+01;
    COFD[2785] = 4.42045763E+00;
    COFD[2786] = -3.59451578E-01;
    COFD[2787] = 1.56056164E-02;
    COFD[2788] = -1.79339327E+01;
    COFD[2789] = 4.91373893E+00;
    COFD[2790] = -4.18747629E-01;
    COFD[2791] = 1.79856610E-02;
    COFD[2792] = -1.39913897E+01;
    COFD[2793] = 3.26384506E+00;
    COFD[2794] = -2.12947087E-01;
    COFD[2795] = 9.39743888E-03;
    COFD[2796] = -1.72496634E+01;
    COFD[2797] = 4.17889917E+00;
    COFD[2798] = -3.29752510E-01;
    COFD[2799] = 1.43850275E-02;
    COFD[2800] = -1.72691500E+01;
    COFD[2801] = 4.17889917E+00;
    COFD[2802] = -3.29752510E-01;
    COFD[2803] = 1.43850275E-02;
    COFD[2804] = -1.86335932E+01;
    COFD[2805] = 4.53572533E+00;
    COFD[2806] = -3.73386925E-01;
    COFD[2807] = 1.61678881E-02;
    COFD[2808] = -1.86479000E+01;
    COFD[2809] = 4.53572533E+00;
    COFD[2810] = -3.73386925E-01;
    COFD[2811] = 1.61678881E-02;
    COFD[2812] = -2.12597312E+01;
    COFD[2813] = 5.24930667E+00;
    COFD[2814] = -4.17435088E-01;
    COFD[2815] = 1.61434424E-02;
    COFD[2816] = -1.86409139E+01;
    COFD[2817] = 4.53572533E+00;
    COFD[2818] = -3.73386925E-01;
    COFD[2819] = 1.61678881E-02;
    COFD[2820] = -1.83455435E+01;
    COFD[2821] = 4.42828044E+00;
    COFD[2822] = -3.60417833E-01;
    COFD[2823] = 1.56455103E-02;
    COFD[2824] = -1.95494668E+01;
    COFD[2825] = 4.90255048E+00;
    COFD[2826] = -4.17368501E-01;
    COFD[2827] = 1.79287358E-02;
    COFD[2828] = -2.19313890E+01;
    COFD[2829] = 5.41841631E+00;
    COFD[2830] = -4.46818971E-01;
    COFD[2831] = 1.77127652E-02;
    COFD[2832] = -2.14224484E+01;
    COFD[2833] = 5.41729961E+00;
    COFD[2834] = -4.73400281E-01;
    COFD[2835] = 1.99269567E-02;
    COFD[2836] = -1.94763688E+01;
    COFD[2837] = 4.87333294E+00;
    COFD[2838] = -4.13769241E-01;
    COFD[2839] = 1.77802244E-02;
    COFD[2840] = -2.10864251E+01;
    COFD[2841] = 5.34286099E+00;
    COFD[2842] = -4.68100992E-01;
    COFD[2843] = 1.98731399E-02;
    COFD[2844] = -2.18190539E+01;
    COFD[2845] = 5.55753905E+00;
    COFD[2846] = -4.88136714E-01;
    COFD[2847] = 2.04294957E-02;
    COFD[2848] = -2.22059540E+01;
    COFD[2849] = 5.51722375E+00;
    COFD[2850] = -4.66081431E-01;
    COFD[2851] = 1.88044011E-02;
    COFD[2852] = -2.15660171E+01;
    COFD[2853] = 5.44803850E+00;
    COFD[2854] = -4.76610560E-01;
    COFD[2855] = 2.00355294E-02;
    COFD[2856] = -1.69491115E+01;
    COFD[2857] = 4.05099737E+00;
    COFD[2858] = -3.13841660E-01;
    COFD[2859] = 1.37218854E-02;
    COFD[2860] = -1.71808106E+01;
    COFD[2861] = 4.17889917E+00;
    COFD[2862] = -3.29752510E-01;
    COFD[2863] = 1.43850275E-02;
    COFD[2864] = -2.19229190E+01;
    COFD[2865] = 5.41841631E+00;
    COFD[2866] = -4.46818971E-01;
    COFD[2867] = 1.77127652E-02;
    COFD[2868] = -1.95263312E+01;
    COFD[2869] = 4.90255048E+00;
    COFD[2870] = -4.17368501E-01;
    COFD[2871] = 1.79287358E-02;
    COFD[2872] = -1.95263312E+01;
    COFD[2873] = 4.90255048E+00;
    COFD[2874] = -4.17368501E-01;
    COFD[2875] = 1.79287358E-02;
    COFD[2876] = -2.10963515E+01;
    COFD[2877] = 5.34286099E+00;
    COFD[2878] = -4.68100992E-01;
    COFD[2879] = 1.98731399E-02;
    COFD[2880] = -2.15575659E+01;
    COFD[2881] = 5.44803850E+00;
    COFD[2882] = -4.76610560E-01;
    COFD[2883] = 2.00355294E-02;
    COFD[2884] = -1.96309042E+01;
    COFD[2885] = 4.95923807E+00;
    COFD[2886] = -4.24176182E-01;
    COFD[2887] = 1.82020215E-02;
    COFD[2888] = -2.22153950E+01;
    COFD[2889] = 5.51722375E+00;
    COFD[2890] = -4.66081431E-01;
    COFD[2891] = 1.88044011E-02;
    COFD[2892] = -2.22107627E+01;
    COFD[2893] = 5.51722375E+00;
    COFD[2894] = -4.66081431E-01;
    COFD[2895] = 1.88044011E-02;
    COFD[2896] = -2.22972410E+01;
    COFD[2897] = 5.47511492E+00;
    COFD[2898] = -4.57098839E-01;
    COFD[2899] = 1.82750523E-02;
    COFD[2900] = -1.59884305E+01;
    COFD[2901] = 3.72220402E+00;
    COFD[2902] = -2.71150591E-01;
    COFD[2903] = 1.18665265E-02;
    COFD[2904] = -1.54460820E+01;
    COFD[2905] = 4.26819983E+00;
    COFD[2906] = -3.40766379E-01;
    COFD[2907] = 1.48393361E-02;
    COFD[2908] = -1.22004324E+01;
    COFD[2909] = 2.80725489E+00;
    COFD[2910] = -1.54291406E-01;
    COFD[2911] = 6.88290911E-03;
    COFD[2912] = -1.49500357E+01;
    COFD[2913] = 3.52327209E+00;
    COFD[2914] = -2.46286208E-01;
    COFD[2915] = 1.08285963E-02;
    COFD[2916] = -1.49718233E+01;
    COFD[2917] = 3.52327209E+00;
    COFD[2918] = -2.46286208E-01;
    COFD[2919] = 1.08285963E-02;
    COFD[2920] = -1.64169433E+01;
    COFD[2921] = 3.89309916E+00;
    COFD[2922] = -2.93528188E-01;
    COFD[2923] = 1.28463177E-02;
    COFD[2924] = -1.64338757E+01;
    COFD[2925] = 3.89309916E+00;
    COFD[2926] = -2.93528188E-01;
    COFD[2927] = 1.28463177E-02;
    COFD[2928] = -2.10440675E+01;
    COFD[2929] = 5.59806282E+00;
    COFD[2930] = -4.87109535E-01;
    COFD[2931] = 2.01370226E-02;
    COFD[2932] = -1.64255964E+01;
    COFD[2933] = 3.89309916E+00;
    COFD[2934] = -2.93528188E-01;
    COFD[2935] = 1.28463177E-02;
    COFD[2936] = -1.60261675E+01;
    COFD[2937] = 3.73312045E+00;
    COFD[2938] = -2.72579779E-01;
    COFD[2939] = 1.19290272E-02;
    COFD[2940] = -1.72828302E+01;
    COFD[2941] = 4.26063341E+00;
    COFD[2942] = -3.39848064E-01;
    COFD[2943] = 1.48021313E-02;
    COFD[2944] = -2.14303479E+01;
    COFD[2945] = 5.59268435E+00;
    COFD[2946] = -4.91232974E-01;
    COFD[2947] = 2.05064746E-02;
    COFD[2948] = -1.94485982E+01;
    COFD[2949] = 4.91446566E+00;
    COFD[2950] = -4.18837152E-01;
    COFD[2951] = 1.79893537E-02;
    COFD[2952] = -1.72316148E+01;
    COFD[2953] = 4.24011069E+00;
    COFD[2954] = -3.37339810E-01;
    COFD[2955] = 1.46996679E-02;
    COFD[2956] = -1.88646070E+01;
    COFD[2957] = 4.72476764E+00;
    COFD[2958] = -3.96306836E-01;
    COFD[2959] = 1.70964541E-02;
    COFD[2960] = -1.99081556E+01;
    COFD[2961] = 5.09311649E+00;
    COFD[2962] = -4.39965178E-01;
    COFD[2963] = 1.88238537E-02;
    COFD[2964] = -2.12621914E+01;
    COFD[2965] = 5.47935225E+00;
    COFD[2966] = -4.80056796E-01;
    COFD[2967] = 2.01607180E-02;
    COFD[2968] = -1.96408127E+01;
    COFD[2969] = 4.95923807E+00;
    COFD[2970] = -4.24176182E-01;
    COFD[2971] = 1.82020215E-02;
    COFD[2972] = -1.46907028E+01;
    COFD[2973] = 3.39229020E+00;
    COFD[2974] = -2.29520232E-01;
    COFD[2975] = 1.01114311E-02;
    COFD[2976] = -1.48738066E+01;
    COFD[2977] = 3.52327209E+00;
    COFD[2978] = -2.46286208E-01;
    COFD[2979] = 1.08285963E-02;
    COFD[2980] = -2.14204185E+01;
    COFD[2981] = 5.59268435E+00;
    COFD[2982] = -4.91232974E-01;
    COFD[2983] = 2.05064746E-02;
    COFD[2984] = -1.72572042E+01;
    COFD[2985] = 4.26063341E+00;
    COFD[2986] = -3.39848064E-01;
    COFD[2987] = 1.48021313E-02;
    COFD[2988] = -1.72572042E+01;
    COFD[2989] = 4.26063341E+00;
    COFD[2990] = -3.39848064E-01;
    COFD[2991] = 1.48021313E-02;
    COFD[2992] = -1.88761387E+01;
    COFD[2993] = 4.72476764E+00;
    COFD[2994] = -3.96306836E-01;
    COFD[2995] = 1.70964541E-02;
    COFD[2996] = -1.96309042E+01;
    COFD[2997] = 4.95923807E+00;
    COFD[2998] = -4.24176182E-01;
    COFD[2999] = 1.82020215E-02;
    COFD[3000] = -1.72414862E+01;
    COFD[3001] = 4.29808578E+00;
    COFD[3002] = -3.44235570E-01;
    COFD[3003] = 1.49727727E-02;
    COFD[3004] = -2.12736225E+01;
    COFD[3005] = 5.47935225E+00;
    COFD[3006] = -4.80056796E-01;
    COFD[3007] = 2.01607180E-02;
    COFD[3008] = -2.12680082E+01;
    COFD[3009] = 5.47935225E+00;
    COFD[3010] = -4.80056796E-01;
    COFD[3011] = 2.01607180E-02;
    COFD[3012] = -2.16544368E+01;
    COFD[3013] = 5.55511977E+00;
    COFD[3014] = -4.87927156E-01;
    COFD[3015] = 2.04245402E-02;
    COFD[3016] = -2.02738958E+01;
    COFD[3017] = 5.10426133E+00;
    COFD[3018] = -4.41256919E-01;
    COFD[3019] = 1.88737290E-02;
    COFD[3020] = -1.94694048E+01;
    COFD[3021] = 5.43830787E+00;
    COFD[3022] = -4.75472880E-01;
    COFD[3023] = 1.99909996E-02;
    COFD[3024] = -1.57045332E+01;
    COFD[3025] = 3.93614244E+00;
    COFD[3026] = -2.99111497E-01;
    COFD[3027] = 1.30888229E-02;
    COFD[3028] = -1.90946745E+01;
    COFD[3029] = 4.84384483E+00;
    COFD[3030] = -4.10265575E-01;
    COFD[3031] = 1.76414287E-02;
    COFD[3032] = -1.91168996E+01;
    COFD[3033] = 4.84384483E+00;
    COFD[3034] = -4.10265575E-01;
    COFD[3035] = 1.76414287E-02;
    COFD[3036] = -2.05284752E+01;
    COFD[3037] = 5.18417470E+00;
    COFD[3038] = -4.49491573E-01;
    COFD[3039] = 1.91438508E-02;
    COFD[3040] = -2.05459419E+01;
    COFD[3041] = 5.18417470E+00;
    COFD[3042] = -4.49491573E-01;
    COFD[3043] = 1.91438508E-02;
    COFD[3044] = -1.87453067E+01;
    COFD[3045] = 3.96926341E+00;
    COFD[3046] = -2.16412264E-01;
    COFD[3047] = 6.06012078E-03;
    COFD[3048] = -2.05373990E+01;
    COFD[3049] = 5.18417470E+00;
    COFD[3050] = -4.49491573E-01;
    COFD[3051] = 1.91438508E-02;
    COFD[3052] = -2.03015042E+01;
    COFD[3053] = 5.11106992E+00;
    COFD[3054] = -4.42047129E-01;
    COFD[3055] = 1.89042990E-02;
    COFD[3056] = -2.11667605E+01;
    COFD[3057] = 5.42846112E+00;
    COFD[3058] = -4.74321870E-01;
    COFD[3059] = 1.99459749E-02;
    COFD[3060] = -2.01111595E+01;
    COFD[3061] = 4.41511629E+00;
    COFD[3062] = -2.84086963E-01;
    COFD[3063] = 9.37586971E-03;
    COFD[3064] = -2.22235122E+01;
    COFD[3065] = 5.54251230E+00;
    COFD[3066] = -4.70946314E-01;
    COFD[3067] = 1.90785869E-02;
    COFD[3068] = -2.11372811E+01;
    COFD[3069] = 5.41773516E+00;
    COFD[3070] = -4.73414338E-01;
    COFD[3071] = 1.99258685E-02;
    COFD[3072] = -2.21153597E+01;
    COFD[3073] = 5.58360799E+00;
    COFD[3074] = -4.82701436E-01;
    COFD[3075] = 1.98437922E-02;
    COFD[3076] = -2.20546151E+01;
    COFD[3077] = 5.44448440E+00;
    COFD[3078] = -4.51529024E-01;
    COFD[3079] = 1.79698119E-02;
    COFD[3080] = -2.09118474E+01;
    COFD[3081] = 4.72895031E+00;
    COFD[3082] = -3.33332771E-01;
    COFD[3083] = 1.18431478E-02;
    COFD[3084] = -2.22255968E+01;
    COFD[3085] = 5.51722375E+00;
    COFD[3086] = -4.66081431E-01;
    COFD[3087] = 1.88044011E-02;
    COFD[3088] = -1.87739217E+01;
    COFD[3089] = 4.71729964E+00;
    COFD[3090] = -3.95432573E-01;
    COFD[3091] = 1.70623691E-02;
    COFD[3092] = -1.90170590E+01;
    COFD[3093] = 4.84384483E+00;
    COFD[3094] = -4.10265575E-01;
    COFD[3095] = 1.76414287E-02;
    COFD[3096] = -2.01009366E+01;
    COFD[3097] = 4.41511629E+00;
    COFD[3098] = -2.84086963E-01;
    COFD[3099] = 9.37586971E-03;
    COFD[3100] = -2.11406662E+01;
    COFD[3101] = 5.42846112E+00;
    COFD[3102] = -4.74321870E-01;
    COFD[3103] = 1.99459749E-02;
    COFD[3104] = -2.11406662E+01;
    COFD[3105] = 5.42846112E+00;
    COFD[3106] = -4.74321870E-01;
    COFD[3107] = 1.99459749E-02;
    COFD[3108] = -2.21272109E+01;
    COFD[3109] = 5.58360799E+00;
    COFD[3110] = -4.82701436E-01;
    COFD[3111] = 1.98437922E-02;
    COFD[3112] = -2.22153950E+01;
    COFD[3113] = 5.51722375E+00;
    COFD[3114] = -4.66081431E-01;
    COFD[3115] = 1.88044011E-02;
    COFD[3116] = -2.12736225E+01;
    COFD[3117] = 5.47935225E+00;
    COFD[3118] = -4.80056796E-01;
    COFD[3119] = 2.01607180E-02;
    COFD[3120] = -2.09236948E+01;
    COFD[3121] = 4.72895031E+00;
    COFD[3122] = -3.33332771E-01;
    COFD[3123] = 1.18431478E-02;
    COFD[3124] = -2.09178748E+01;
    COFD[3125] = 4.72895031E+00;
    COFD[3126] = -3.33332771E-01;
    COFD[3127] = 1.18431478E-02;
    COFD[3128] = -2.07191944E+01;
    COFD[3129] = 4.56211059E+00;
    COFD[3130] = -3.06895158E-01;
    COFD[3131] = 1.05100393E-02;
    COFD[3132] = -2.02693653E+01;
    COFD[3133] = 5.10426133E+00;
    COFD[3134] = -4.41256919E-01;
    COFD[3135] = 1.88737290E-02;
    COFD[3136] = -1.94691430E+01;
    COFD[3137] = 5.43830787E+00;
    COFD[3138] = -4.75472880E-01;
    COFD[3139] = 1.99909996E-02;
    COFD[3140] = -1.57040212E+01;
    COFD[3141] = 3.93614244E+00;
    COFD[3142] = -2.99111497E-01;
    COFD[3143] = 1.30888229E-02;
    COFD[3144] = -1.90915649E+01;
    COFD[3145] = 4.84384483E+00;
    COFD[3146] = -4.10265575E-01;
    COFD[3147] = 1.76414287E-02;
    COFD[3148] = -1.91136491E+01;
    COFD[3149] = 4.84384483E+00;
    COFD[3150] = -4.10265575E-01;
    COFD[3151] = 1.76414287E-02;
    COFD[3152] = -2.05235731E+01;
    COFD[3153] = 5.18417470E+00;
    COFD[3154] = -4.49491573E-01;
    COFD[3155] = 1.91438508E-02;
    COFD[3156] = -2.05408665E+01;
    COFD[3157] = 5.18417470E+00;
    COFD[3158] = -4.49491573E-01;
    COFD[3159] = 1.91438508E-02;
    COFD[3160] = -1.87419199E+01;
    COFD[3161] = 3.96926341E+00;
    COFD[3162] = -2.16412264E-01;
    COFD[3163] = 6.06012078E-03;
    COFD[3164] = -2.05324091E+01;
    COFD[3165] = 5.18417470E+00;
    COFD[3166] = -4.49491573E-01;
    COFD[3167] = 1.91438508E-02;
    COFD[3168] = -2.02969740E+01;
    COFD[3169] = 5.11106992E+00;
    COFD[3170] = -4.42047129E-01;
    COFD[3171] = 1.89042990E-02;
    COFD[3172] = -2.11637902E+01;
    COFD[3173] = 5.42846112E+00;
    COFD[3174] = -4.74321870E-01;
    COFD[3175] = 1.99459749E-02;
    COFD[3176] = -2.01064363E+01;
    COFD[3177] = 4.41511629E+00;
    COFD[3178] = -2.84086963E-01;
    COFD[3179] = 9.37586971E-03;
    COFD[3180] = -2.22176950E+01;
    COFD[3181] = 5.54251230E+00;
    COFD[3182] = -4.70946314E-01;
    COFD[3183] = 1.90785869E-02;
    COFD[3184] = -2.11341653E+01;
    COFD[3185] = 5.41773516E+00;
    COFD[3186] = -4.73414338E-01;
    COFD[3187] = 1.99258685E-02;
    COFD[3188] = -2.21110290E+01;
    COFD[3189] = 5.58360799E+00;
    COFD[3190] = -4.82701436E-01;
    COFD[3191] = 1.98437922E-02;
    COFD[3192] = -2.20500806E+01;
    COFD[3193] = 5.44448440E+00;
    COFD[3194] = -4.51529024E-01;
    COFD[3195] = 1.79698119E-02;
    COFD[3196] = -2.09061629E+01;
    COFD[3197] = 4.72895031E+00;
    COFD[3198] = -3.33332771E-01;
    COFD[3199] = 1.18431478E-02;
    COFD[3200] = -2.22208695E+01;
    COFD[3201] = 5.51722375E+00;
    COFD[3202] = -4.66081431E-01;
    COFD[3203] = 1.88044011E-02;
    COFD[3204] = -1.87714196E+01;
    COFD[3205] = 4.71729964E+00;
    COFD[3206] = -3.95432573E-01;
    COFD[3207] = 1.70623691E-02;
    COFD[3208] = -1.90143953E+01;
    COFD[3209] = 4.84384483E+00;
    COFD[3210] = -4.10265575E-01;
    COFD[3211] = 1.76414287E-02;
    COFD[3212] = -2.00963085E+01;
    COFD[3213] = 4.41511629E+00;
    COFD[3214] = -2.84086963E-01;
    COFD[3215] = 9.37586971E-03;
    COFD[3216] = -2.11378465E+01;
    COFD[3217] = 5.42846112E+00;
    COFD[3218] = -4.74321870E-01;
    COFD[3219] = 1.99459749E-02;
    COFD[3220] = -2.11378465E+01;
    COFD[3221] = 5.42846112E+00;
    COFD[3222] = -4.74321870E-01;
    COFD[3223] = 1.99459749E-02;
    COFD[3224] = -2.21227768E+01;
    COFD[3225] = 5.58360799E+00;
    COFD[3226] = -4.82701436E-01;
    COFD[3227] = 1.98437922E-02;
    COFD[3228] = -2.22107627E+01;
    COFD[3229] = 5.51722375E+00;
    COFD[3230] = -4.66081431E-01;
    COFD[3231] = 1.88044011E-02;
    COFD[3232] = -2.12680082E+01;
    COFD[3233] = 5.47935225E+00;
    COFD[3234] = -4.80056796E-01;
    COFD[3235] = 2.01607180E-02;
    COFD[3236] = -2.09178748E+01;
    COFD[3237] = 4.72895031E+00;
    COFD[3238] = -3.33332771E-01;
    COFD[3239] = 1.18431478E-02;
    COFD[3240] = -2.09121217E+01;
    COFD[3241] = 4.72895031E+00;
    COFD[3242] = -3.33332771E-01;
    COFD[3243] = 1.18431478E-02;
    COFD[3244] = -2.07133089E+01;
    COFD[3245] = 4.56211059E+00;
    COFD[3246] = -3.06895158E-01;
    COFD[3247] = 1.05100393E-02;
    COFD[3248] = -2.05802296E+01;
    COFD[3249] = 5.16117916E+00;
    COFD[3250] = -4.46897404E-01;
    COFD[3251] = 1.90470443E-02;
    COFD[3252] = -1.98806372E+01;
    COFD[3253] = 5.52555673E+00;
    COFD[3254] = -4.84999851E-01;
    COFD[3255] = 2.03334931E-02;
    COFD[3256] = -1.61116686E+01;
    COFD[3257] = 4.04227735E+00;
    COFD[3258] = -3.12745253E-01;
    COFD[3259] = 1.36756977E-02;
    COFD[3260] = -1.95314689E+01;
    COFD[3261] = 4.95249173E+00;
    COFD[3262] = -4.23376552E-01;
    COFD[3263] = 1.81703714E-02;
    COFD[3264] = -1.95538303E+01;
    COFD[3265] = 4.95249173E+00;
    COFD[3266] = -4.23376552E-01;
    COFD[3267] = 1.81703714E-02;
    COFD[3268] = -2.09481051E+01;
    COFD[3269] = 5.28755355E+00;
    COFD[3270] = -4.61641920E-01;
    COFD[3271] = 1.96208961E-02;
    COFD[3272] = -2.09657408E+01;
    COFD[3273] = 5.28755355E+00;
    COFD[3274] = -4.61641920E-01;
    COFD[3275] = 1.96208961E-02;
    COFD[3276] = -1.84538368E+01;
    COFD[3277] = 3.75912079E+00;
    COFD[3278] = -1.84235105E-01;
    COFD[3279] = 4.47800951E-03;
    COFD[3280] = -2.09571146E+01;
    COFD[3281] = 5.28755355E+00;
    COFD[3282] = -4.61641920E-01;
    COFD[3283] = 1.96208961E-02;
    COFD[3284] = -2.06066440E+01;
    COFD[3285] = 5.16748146E+00;
    COFD[3286] = -4.47594939E-01;
    COFD[3287] = 1.90724110E-02;
    COFD[3288] = -2.15588759E+01;
    COFD[3289] = 5.51982454E+00;
    COFD[3290] = -4.84452039E-01;
    COFD[3291] = 2.03175522E-02;
    COFD[3292] = -1.97649065E+01;
    COFD[3293] = 4.18758010E+00;
    COFD[3294] = -2.49327776E-01;
    COFD[3295] = 7.66559103E-03;
    COFD[3296] = -2.23098172E+01;
    COFD[3297] = 5.49916900E+00;
    COFD[3298] = -4.61818485E-01;
    COFD[3299] = 1.85431163E-02;
    COFD[3300] = -2.15067581E+01;
    COFD[3301] = 5.49964831E+00;
    COFD[3302] = -4.82275380E-01;
    COFD[3303] = 2.02405072E-02;
    COFD[3304] = -2.22982472E+01;
    COFD[3305] = 5.58471203E+00;
    COFD[3306] = -4.79905311E-01;
    COFD[3307] = 1.96058913E-02;
    COFD[3308] = -2.20600677E+01;
    COFD[3309] = 5.36785880E+00;
    COFD[3310] = -4.37683002E-01;
    COFD[3311] = 1.72144393E-02;
    COFD[3312] = -2.07072145E+01;
    COFD[3313] = 4.56211059E+00;
    COFD[3314] = -3.06895158E-01;
    COFD[3315] = 1.05100393E-02;
    COFD[3316] = -2.23075353E+01;
    COFD[3317] = 5.47511492E+00;
    COFD[3318] = -4.57098839E-01;
    COFD[3319] = 1.82750523E-02;
    COFD[3320] = -1.91468116E+01;
    COFD[3321] = 4.80472875E+00;
    COFD[3322] = -4.05752803E-01;
    COFD[3323] = 1.74685236E-02;
    COFD[3324] = -1.94534227E+01;
    COFD[3325] = 4.95249173E+00;
    COFD[3326] = -4.23376552E-01;
    COFD[3327] = 1.81703714E-02;
    COFD[3328] = -1.97545910E+01;
    COFD[3329] = 4.18758010E+00;
    COFD[3330] = -2.49327776E-01;
    COFD[3331] = 7.66559103E-03;
    COFD[3332] = -2.15326361E+01;
    COFD[3333] = 5.51982454E+00;
    COFD[3334] = -4.84452039E-01;
    COFD[3335] = 2.03175522E-02;
    COFD[3336] = -2.15326361E+01;
    COFD[3337] = 5.51982454E+00;
    COFD[3338] = -4.84452039E-01;
    COFD[3339] = 2.03175522E-02;
    COFD[3340] = -2.23101989E+01;
    COFD[3341] = 5.58471203E+00;
    COFD[3342] = -4.79905311E-01;
    COFD[3343] = 1.96058913E-02;
    COFD[3344] = -2.22972410E+01;
    COFD[3345] = 5.47511492E+00;
    COFD[3346] = -4.57098839E-01;
    COFD[3347] = 1.82750523E-02;
    COFD[3348] = -2.16544368E+01;
    COFD[3349] = 5.55511977E+00;
    COFD[3350] = -4.87927156E-01;
    COFD[3351] = 2.04245402E-02;
    COFD[3352] = -2.07191944E+01;
    COFD[3353] = 4.56211059E+00;
    COFD[3354] = -3.06895158E-01;
    COFD[3355] = 1.05100393E-02;
    COFD[3356] = -2.07133089E+01;
    COFD[3357] = 4.56211059E+00;
    COFD[3358] = -3.06895158E-01;
    COFD[3359] = 1.05100393E-02;
    COFD[3360] = -2.04221575E+01;
    COFD[3361] = 4.35883159E+00;
    COFD[3362] = -2.75434484E-01;
    COFD[3363] = 8.94819804E-03;
}


/*List of specs with small weight, dim NLITE */
void egtransetKTDIF(int* KTDIF) {
    KTDIF[0] = 1;
    KTDIF[1] = 2;
}


/*Poly fits for thermal diff ratios, dim NO*NLITE*KK */
void egtransetCOFTD(amrex::Real* COFTD) {
    COFTD[0] = 2.01521643E-01;
    COFTD[1] = 5.62744089E-04;
    COFTD[2] = -3.08519239E-07;
    COFTD[3] = 5.22805986E-11;
    COFTD[4] = 0.00000000E+00;
    COFTD[5] = 0.00000000E+00;
    COFTD[6] = 0.00000000E+00;
    COFTD[7] = 0.00000000E+00;
    COFTD[8] = 1.44152190E-01;
    COFTD[9] = 7.99993584E-05;
    COFTD[10] = -4.89707442E-08;
    COFTD[11] = 9.14277269E-12;
    COFTD[12] = 2.35283119E-01;
    COFTD[13] = 4.65670599E-04;
    COFTD[14] = -2.60939824E-07;
    COFTD[15] = 4.49271822E-11;
    COFTD[16] = 2.37053352E-01;
    COFTD[17] = 4.69174231E-04;
    COFTD[18] = -2.62903094E-07;
    COFTD[19] = 4.52652072E-11;
    COFTD[20] = 1.79840299E-01;
    COFTD[21] = 6.01722902E-04;
    COFTD[22] = -3.26433894E-07;
    COFTD[23] = 5.49112302E-11;
    COFTD[24] = 1.80513677E-01;
    COFTD[25] = 6.03975942E-04;
    COFTD[26] = -3.27656165E-07;
    COFTD[27] = 5.51168351E-11;
    COFTD[28] = -1.74352698E-01;
    COFTD[29] = 8.62246873E-04;
    COFTD[30] = -3.79545489E-07;
    COFTD[31] = 5.60262093E-11;
    COFTD[32] = 1.80186965E-01;
    COFTD[33] = 6.02882805E-04;
    COFTD[34] = -3.27063140E-07;
    COFTD[35] = 5.50170790E-11;
    COFTD[36] = 2.00119897E-01;
    COFTD[37] = 5.64793704E-04;
    COFTD[38] = -3.09445484E-07;
    COFTD[39] = 5.24139335E-11;
    COFTD[40] = 1.00039110E-01;
    COFTD[41] = 6.50468660E-04;
    COFTD[42] = -3.41778999E-07;
    COFTD[43] = 5.62779132E-11;
    COFTD[44] = -1.61357564E-01;
    COFTD[45] = 9.05920260E-04;
    COFTD[46] = -4.07879153E-07;
    COFTD[47] = 6.10626290E-11;
    COFTD[48] = -2.00309448E-02;
    COFTD[49] = 8.50440115E-04;
    COFTD[50] = -4.21064468E-07;
    COFTD[51] = 6.67959710E-11;
    COFTD[52] = 1.05124122E-01;
    COFTD[53] = 6.50665957E-04;
    COFTD[54] = -3.42564538E-07;
    COFTD[55] = 5.64804120E-11;
    COFTD[56] = 1.63245097E-02;
    COFTD[57] = 7.90133388E-04;
    COFTD[58] = -3.98292458E-07;
    COFTD[59] = 6.38851432E-11;
    COFTD[60] = -5.08744745E-02;
    COFTD[61] = 8.54342586E-04;
    COFTD[62] = -4.15926453E-07;
    COFTD[63] = 6.53063261E-11;
    COFTD[64] = -1.41640506E-01;
    COFTD[65] = 9.21404324E-04;
    COFTD[66] = -4.23210110E-07;
    COFTD[67] = 6.41400322E-11;
    COFTD[68] = -2.72323768E-02;
    COFTD[69] = 8.39184413E-04;
    COFTD[70] = -4.13849924E-07;
    COFTD[71] = 6.54928043E-11;
    COFTD[72] = 2.49272491E-01;
    COFTD[73] = 4.08682510E-04;
    COFTD[74] = -2.31943878E-07;
    COFTD[75] = 4.03271405E-11;
    COFTD[76] = 2.28560867E-01;
    COFTD[77] = 4.52365967E-04;
    COFTD[78] = -2.53484536E-07;
    COFTD[79] = 4.36435719E-11;
    COFTD[80] = -1.60981264E-01;
    COFTD[81] = 9.03807572E-04;
    COFTD[82] = -4.06927941E-07;
    COFTD[83] = 6.09202254E-11;
    COFTD[84] = 9.90752318E-02;
    COFTD[85] = 6.44201384E-04;
    COFTD[86] = -3.38485953E-07;
    COFTD[87] = 5.57356746E-11;
    COFTD[88] = 9.90752318E-02;
    COFTD[89] = 6.44201384E-04;
    COFTD[90] = -3.38485953E-07;
    COFTD[91] = 5.57356746E-11;
    COFTD[92] = 1.63717489E-02;
    COFTD[93] = 7.92419842E-04;
    COFTD[94] = -3.99445020E-07;
    COFTD[95] = 6.40700113E-11;
    COFTD[96] = -2.71690558E-02;
    COFTD[97] = 8.37233133E-04;
    COFTD[98] = -4.12887636E-07;
    COFTD[99] = 6.53405197E-11;
    COFTD[100] = 9.86934401E-02;
    COFTD[101] = 7.20974863E-04;
    COFTD[102] = -3.77135221E-07;
    COFTD[103] = 6.19202579E-11;
    COFTD[104] = -1.41951848E-01;
    COFTD[105] = 9.23429679E-04;
    COFTD[106] = -4.24140376E-07;
    COFTD[107] = 6.42810196E-11;
    COFTD[108] = -1.41799739E-01;
    COFTD[109] = 9.22440172E-04;
    COFTD[110] = -4.23685885E-07;
    COFTD[111] = 6.42121388E-11;
    COFTD[112] = -1.55536623E-01;
    COFTD[113] = 9.26290092E-04;
    COFTD[114] = -4.20679731E-07;
    COFTD[115] = 6.33165565E-11;
    COFTD[116] = 4.31331269E-01;
    COFTD[117] = 9.20536800E-05;
    COFTD[118] = -5.94509616E-08;
    COFTD[119] = 1.21437993E-11;
    COFTD[120] = -1.44152190E-01;
    COFTD[121] = -7.99993584E-05;
    COFTD[122] = 4.89707442E-08;
    COFTD[123] = -9.14277269E-12;
    COFTD[124] = 0.00000000E+00;
    COFTD[125] = 0.00000000E+00;
    COFTD[126] = 0.00000000E+00;
    COFTD[127] = 0.00000000E+00;
    COFTD[128] = 4.06682492E-01;
    COFTD[129] = 3.84705248E-05;
    COFTD[130] = -2.54846868E-08;
    COFTD[131] = 5.86302354E-12;
    COFTD[132] = 4.12895615E-01;
    COFTD[133] = 3.90582612E-05;
    COFTD[134] = -2.58740310E-08;
    COFTD[135] = 5.95259633E-12;
    COFTD[136] = 4.26579943E-01;
    COFTD[137] = 1.20407274E-04;
    COFTD[138] = -7.67298757E-08;
    COFTD[139] = 1.52090336E-11;
    COFTD[140] = 4.29789463E-01;
    COFTD[141] = 1.21313199E-04;
    COFTD[142] = -7.73071792E-08;
    COFTD[143] = 1.53234639E-11;
    COFTD[144] = 2.27469146E-02;
    COFTD[145] = 6.73078907E-04;
    COFTD[146] = -3.40935843E-07;
    COFTD[147] = 5.48499211E-11;
    COFTD[148] = 4.28230888E-01;
    COFTD[149] = 1.20873273E-04;
    COFTD[150] = -7.70268349E-08;
    COFTD[151] = 1.52678954E-11;
    COFTD[152] = 4.30605547E-01;
    COFTD[153] = 9.35961902E-05;
    COFTD[154] = -6.03983623E-08;
    COFTD[155] = 1.23115170E-11;
    COFTD[156] = 3.31191185E-01;
    COFTD[157] = 1.81326714E-04;
    COFTD[158] = -1.11096391E-07;
    COFTD[159] = 2.07635959E-11;
    COFTD[160] = 1.22693382E-01;
    COFTD[161] = 6.21278143E-04;
    COFTD[162] = -3.29965208E-07;
    COFTD[163] = 5.47161548E-11;
    COFTD[164] = 2.93191523E-01;
    COFTD[165] = 4.01430006E-04;
    COFTD[166] = -2.30705763E-07;
    COFTD[167] = 4.05176586E-11;
    COFTD[168] = 3.39557243E-01;
    COFTD[169] = 1.79335036E-04;
    COFTD[170] = -1.10135705E-07;
    COFTD[171] = 2.06427239E-11;
    COFTD[172] = 3.05613225E-01;
    COFTD[173] = 3.24505886E-04;
    COFTD[174] = -1.89889572E-07;
    COFTD[175] = 3.38663465E-11;
    COFTD[176] = 2.49017478E-01;
    COFTD[177] = 4.29036573E-04;
    COFTD[178] = -2.42668617E-07;
    COFTD[179] = 4.20801371E-11;
    COFTD[180] = 1.59288984E-01;
    COFTD[181] = 6.02833801E-04;
    COFTD[182] = -3.24837576E-07;
    COFTD[183] = 5.43909010E-11;
    COFTD[184] = 2.74036956E-01;
    COFTD[185] = 3.96249742E-04;
    COFTD[186] = -2.26857964E-07;
    COFTD[187] = 3.97176979E-11;
    COFTD[188] = 3.82245475E-01;
    COFTD[189] = 1.47167068E-05;
    COFTD[190] = -9.75995257E-09;
    COFTD[191] = 2.83217152E-12;
    COFTD[192] = 3.83439056E-01;
    COFTD[193] = 3.62717894E-05;
    COFTD[194] = -2.40281409E-08;
    COFTD[195] = 5.52792966E-12;
    COFTD[196] = 1.22119780E-01;
    COFTD[197] = 6.18373616E-04;
    COFTD[198] = -3.28422593E-07;
    COFTD[199] = 5.44603522E-11;
    COFTD[200] = 3.24747031E-01;
    COFTD[201] = 1.77798548E-04;
    COFTD[202] = -1.08934732E-07;
    COFTD[203] = 2.03595881E-11;
    COFTD[204] = 3.24747031E-01;
    COFTD[205] = 1.77798548E-04;
    COFTD[206] = -1.08934732E-07;
    COFTD[207] = 2.03595881E-11;
    COFTD[208] = 3.07392263E-01;
    COFTD[209] = 3.26394901E-04;
    COFTD[210] = -1.90994958E-07;
    COFTD[211] = 3.40634894E-11;
    COFTD[212] = 2.72759599E-01;
    COFTD[213] = 3.94402719E-04;
    COFTD[214] = -2.25800520E-07;
    COFTD[215] = 3.95325634E-11;
    COFTD[216] = 3.86107464E-01;
    COFTD[217] = 2.28760446E-04;
    COFTD[218] = -1.39425040E-07;
    COFTD[219] = 2.58989754E-11;
    COFTD[220] = 1.59991186E-01;
    COFTD[221] = 6.05491303E-04;
    COFTD[222] = -3.26269573E-07;
    COFTD[223] = 5.46306751E-11;
    COFTD[224] = 1.59647939E-01;
    COFTD[225] = 6.04192274E-04;
    COFTD[226] = -3.25569591E-07;
    COFTD[227] = 5.45134698E-11;
    COFTD[228] = 1.41968940E-01;
    COFTD[229] = 6.31753578E-04;
    COFTD[230] = -3.37603052E-07;
    COFTD[231] = 5.62125242E-11;
}

#ifndef AMREX_USE_GPU
/* Initializes parameter database */
void CKINIT()
{

    // (0):  2.000000 H + M <=> H2 + M
    kiv[15] = {1,2};
    nuv[15] = {-2.0,1};
    // (0):  2.000000 H + M <=> H2 + M
    fwd_A[15]     = 1.78e+18;
    fwd_beta[15]  = -1;
    fwd_Ea[15]    = 0;
    prefactor_units[15]  = 1.0000000000000002e-12;
    activation_units[15] = 0.50321666580471969;
    phase_units[15]      = pow(10,-12.000000);
    is_PD[15] = 0;
    nTB[15] = 16;
    TB[15] = (amrex::Real *) malloc(16 * sizeof(amrex::Real));
    TBid[15] = (int *) malloc(16 * sizeof(int));
    TBid[15][0] = 2; TB[15][0] = 0; // H2
    TBid[15][1] = 7; TB[15][1] = 0; // H2O
    TBid[15][2] = 18; TB[15][2] = 0; // C
    TBid[15][3] = 19; TB[15][3] = 0; // CH
    TBid[15][4] = 20; TB[15][4] = 0; // HCO
    TBid[15][5] = 21; TB[15][5] = 0; // TXCH2
    TBid[15][6] = 12; TB[15][6] = 0; // CO2
    TBid[15][7] = 22; TB[15][7] = 0; // SXCH2
    TBid[15][8] = 13; TB[15][8] = 2; // CH4
    TBid[15][9] = 23; TB[15][9] = 0; // C2H3
    TBid[15][10] = 24; TB[15][10] = 0; // C2H5
    TBid[15][11] = 25; TB[15][11] = 0; // HCCO
    TBid[15][12] = 26; TB[15][12] = 0; // CH3CHO
    TBid[15][13] = 27; TB[15][13] = 0; // CH2CHO
    TBid[15][14] = 28; TB[15][14] = 0; // C2H5O
    TBid[15][15] = 17; TB[15][15] = 3; // C2H6

    // (1):  2.000000 H + H2 <=> 2.000000 H2
    kiv[20] = {1,2,2};
    nuv[20] = {-2.0,-1,2.0};
    // (1):  2.000000 H + H2 <=> 2.000000 H2
    fwd_A[20]     = 90000000000000000;
    fwd_beta[20]  = -0.59999999999999998;
    fwd_Ea[20]    = 0;
    prefactor_units[20]  = 1.0000000000000002e-12;
    activation_units[20] = 0.50321666580471969;
    phase_units[20]      = pow(10,-18.000000);
    is_PD[20] = 0;
    nTB[20] = 0;

    // (2):  O + H2 <=> H + OH
    kiv[21] = {3,2,1,4};
    nuv[21] = {-1,-1,1,1};
    // (2):  O + H2 <=> H + OH
    fwd_A[21]     = 45900;
    fwd_beta[21]  = 2.7000000000000002;
    fwd_Ea[21]    = 6259.5600000000004;
    prefactor_units[21]  = 1.0000000000000002e-06;
    activation_units[21] = 0.50321666580471969;
    phase_units[21]      = pow(10,-12.000000);
    is_PD[21] = 0;
    nTB[21] = 0;

    // (3):  H + O + M <=> OH + M
    kiv[16] = {1,3,4};
    nuv[16] = {-1,-1,1};
    // (3):  H + O + M <=> OH + M
    fwd_A[16]     = 9.43e+18;
    fwd_beta[16]  = -1;
    fwd_Ea[16]    = 0;
    prefactor_units[16]  = 1.0000000000000002e-12;
    activation_units[16] = 0.50321666580471969;
    phase_units[16]      = pow(10,-12.000000);
    is_PD[16] = 0;
    nTB[16] = 17;
    TB[16] = (amrex::Real *) malloc(17 * sizeof(amrex::Real));
    TBid[16] = (int *) malloc(17 * sizeof(int));
    TBid[16][0] = 2; TB[16][0] = 2; // H2
    TBid[16][1] = 7; TB[16][1] = 12; // H2O
    TBid[16][2] = 18; TB[16][2] = 0; // C
    TBid[16][3] = 9; TB[16][3] = 1.75; // CO
    TBid[16][4] = 19; TB[16][4] = 0; // CH
    TBid[16][5] = 20; TB[16][5] = 0; // HCO
    TBid[16][6] = 21; TB[16][6] = 0; // TXCH2
    TBid[16][7] = 12; TB[16][7] = 3.6000000000000001; // CO2
    TBid[16][8] = 22; TB[16][8] = 0; // SXCH2
    TBid[16][9] = 13; TB[16][9] = 2; // CH4
    TBid[16][10] = 23; TB[16][10] = 0; // C2H3
    TBid[16][11] = 24; TB[16][11] = 0; // C2H5
    TBid[16][12] = 25; TB[16][12] = 0; // HCCO
    TBid[16][13] = 26; TB[16][13] = 0; // CH3CHO
    TBid[16][14] = 27; TB[16][14] = 0; // CH2CHO
    TBid[16][15] = 28; TB[16][15] = 0; // C2H5O
    TBid[16][16] = 17; TB[16][16] = 3; // C2H6

    // (4):  2.000000 O + M <=> O2 + M
    kiv[17] = {3,5};
    nuv[17] = {-2.0,1};
    // (4):  2.000000 O + M <=> O2 + M
    fwd_A[17]     = 1.2e+17;
    fwd_beta[17]  = -1;
    fwd_Ea[17]    = 0;
    prefactor_units[17]  = 1.0000000000000002e-12;
    activation_units[17] = 0.50321666580471969;
    phase_units[17]      = pow(10,-12.000000);
    is_PD[17] = 0;
    nTB[17] = 17;
    TB[17] = (amrex::Real *) malloc(17 * sizeof(amrex::Real));
    TBid[17] = (int *) malloc(17 * sizeof(int));
    TBid[17][0] = 2; TB[17][0] = 2.3999999999999999; // H2
    TBid[17][1] = 7; TB[17][1] = 15.4; // H2O
    TBid[17][2] = 18; TB[17][2] = 0; // C
    TBid[17][3] = 9; TB[17][3] = 1.75; // CO
    TBid[17][4] = 19; TB[17][4] = 0; // CH
    TBid[17][5] = 20; TB[17][5] = 0; // HCO
    TBid[17][6] = 21; TB[17][6] = 0; // TXCH2
    TBid[17][7] = 12; TB[17][7] = 3.6000000000000001; // CO2
    TBid[17][8] = 22; TB[17][8] = 0; // SXCH2
    TBid[17][9] = 13; TB[17][9] = 2; // CH4
    TBid[17][10] = 23; TB[17][10] = 0; // C2H3
    TBid[17][11] = 24; TB[17][11] = 0; // C2H5
    TBid[17][12] = 25; TB[17][12] = 0; // HCCO
    TBid[17][13] = 26; TB[17][13] = 0; // CH3CHO
    TBid[17][14] = 27; TB[17][14] = 0; // CH2CHO
    TBid[17][15] = 28; TB[17][15] = 0; // C2H5O
    TBid[17][16] = 17; TB[17][16] = 3; // C2H6

    // (5):  2.000000 OH (+M) <=> H2O2 (+M)
    kiv[0] = {4,6};
    nuv[0] = {-2.0,1};
    // (5):  2.000000 OH (+M) <=> H2O2 (+M)
    fwd_A[0]     = 111000000000000;
    fwd_beta[0]  = -0.37;
    fwd_Ea[0]    = 0;
    low_A[0]     = 2.01e+17;
    low_beta[0]  = -0.57999999999999996;
    low_Ea[0]    = -2292.0700000000002;
    troe_a[0]    = 0.73460000000000003;
    troe_Tsss[0] = 94;
    troe_Ts[0]   = 1756;
    troe_Tss[0]  = 5182;
    troe_len[0]  = 4;
    prefactor_units[0]  = 1.0000000000000002e-06;
    activation_units[0] = 0.50321666580471969;
    phase_units[0]      = pow(10,-12.000000);
    is_PD[0] = 1;
    nTB[0] = 17;
    TB[0] = (amrex::Real *) malloc(17 * sizeof(amrex::Real));
    TBid[0] = (int *) malloc(17 * sizeof(int));
    TBid[0][0] = 2; TB[0][0] = 2; // H2
    TBid[0][1] = 7; TB[0][1] = 12; // H2O
    TBid[0][2] = 18; TB[0][2] = 0; // C
    TBid[0][3] = 9; TB[0][3] = 1.75; // CO
    TBid[0][4] = 19; TB[0][4] = 0; // CH
    TBid[0][5] = 20; TB[0][5] = 0; // HCO
    TBid[0][6] = 21; TB[0][6] = 0; // TXCH2
    TBid[0][7] = 12; TB[0][7] = 3.6000000000000001; // CO2
    TBid[0][8] = 22; TB[0][8] = 0; // SXCH2
    TBid[0][9] = 13; TB[0][9] = 2; // CH4
    TBid[0][10] = 23; TB[0][10] = 0; // C2H3
    TBid[0][11] = 24; TB[0][11] = 0; // C2H5
    TBid[0][12] = 25; TB[0][12] = 0; // HCCO
    TBid[0][13] = 26; TB[0][13] = 0; // CH3CHO
    TBid[0][14] = 27; TB[0][14] = 0; // CH2CHO
    TBid[0][15] = 28; TB[0][15] = 0; // C2H5O
    TBid[0][16] = 17; TB[0][16] = 3; // C2H6

    // (6):  H + OH + M <=> H2O + M
    kiv[18] = {1,4,7};
    nuv[18] = {-1,-1,1};
    // (6):  H + OH + M <=> H2O + M
    fwd_A[18]     = 4.4e+22;
    fwd_beta[18]  = -2;
    fwd_Ea[18]    = 0;
    prefactor_units[18]  = 1.0000000000000002e-12;
    activation_units[18] = 0.50321666580471969;
    phase_units[18]      = pow(10,-12.000000);
    is_PD[18] = 0;
    nTB[18] = 17;
    TB[18] = (amrex::Real *) malloc(17 * sizeof(amrex::Real));
    TBid[18] = (int *) malloc(17 * sizeof(int));
    TBid[18][0] = 2; TB[18][0] = 2; // H2
    TBid[18][1] = 7; TB[18][1] = 6.2999999999999998; // H2O
    TBid[18][2] = 18; TB[18][2] = 0; // C
    TBid[18][3] = 9; TB[18][3] = 1.75; // CO
    TBid[18][4] = 19; TB[18][4] = 0; // CH
    TBid[18][5] = 20; TB[18][5] = 0; // HCO
    TBid[18][6] = 21; TB[18][6] = 0; // TXCH2
    TBid[18][7] = 12; TB[18][7] = 3.6000000000000001; // CO2
    TBid[18][8] = 22; TB[18][8] = 0; // SXCH2
    TBid[18][9] = 13; TB[18][9] = 2; // CH4
    TBid[18][10] = 23; TB[18][10] = 0; // C2H3
    TBid[18][11] = 24; TB[18][11] = 0; // C2H5
    TBid[18][12] = 25; TB[18][12] = 0; // HCCO
    TBid[18][13] = 26; TB[18][13] = 0; // CH3CHO
    TBid[18][14] = 27; TB[18][14] = 0; // CH2CHO
    TBid[18][15] = 28; TB[18][15] = 0; // C2H5O
    TBid[18][16] = 17; TB[18][16] = 3; // C2H6

    // (7):  2.000000 OH <=> O + H2O
    kiv[22] = {4,3,7};
    nuv[22] = {-2.0,1,1};
    // (7):  2.000000 OH <=> O + H2O
    fwd_A[22]     = 39700;
    fwd_beta[22]  = 2.3999999999999999;
    fwd_Ea[22]    = -2110.4200000000001;
    prefactor_units[22]  = 1.0000000000000002e-06;
    activation_units[22] = 0.50321666580471969;
    phase_units[22]      = pow(10,-12.000000);
    is_PD[22] = 0;
    nTB[22] = 0;

    // (8):  OH + H2 <=> H + H2O
    kiv[23] = {4,2,1,7};
    nuv[23] = {-1,-1,1,1};
    // (8):  OH + H2 <=> H + H2O
    fwd_A[23]     = 173000000;
    fwd_beta[23]  = 1.51;
    fwd_Ea[23]    = 3429.73;
    prefactor_units[23]  = 1.0000000000000002e-06;
    activation_units[23] = 0.50321666580471969;
    phase_units[23]      = pow(10,-12.000000);
    is_PD[23] = 0;
    nTB[23] = 0;

    // (9):  2.000000 H + H2O <=> H2 + H2O
    kiv[24] = {1,7,2,7};
    nuv[24] = {-2.0,-1,1,1};
    // (9):  2.000000 H + H2O <=> H2 + H2O
    fwd_A[24]     = 5.62e+19;
    fwd_beta[24]  = -1.25;
    fwd_Ea[24]    = 0;
    prefactor_units[24]  = 1.0000000000000002e-12;
    activation_units[24] = 0.50321666580471969;
    phase_units[24]      = pow(10,-18.000000);
    is_PD[24] = 0;
    nTB[24] = 0;

    // (10):  H + O2 (+M) <=> HO2 (+M)
    kiv[1] = {1,5,8};
    nuv[1] = {-1,-1,1};
    // (10):  H + O2 (+M) <=> HO2 (+M)
    fwd_A[1]     = 5120000000000;
    fwd_beta[1]  = 0.44;
    fwd_Ea[1]    = 0;
    low_A[1]     = 6.33e+19;
    low_beta[1]  = -1.3999999999999999;
    low_Ea[1]    = 0;
    troe_a[1]    = 0.5;
    troe_Tsss[1] = 0;
    troe_Ts[1]   = 10000000000;
    troe_len[1]  = 3;
    prefactor_units[1]  = 1.0000000000000002e-06;
    activation_units[1] = 0.50321666580471969;
    phase_units[1]      = pow(10,-12.000000);
    is_PD[1] = 1;
    nTB[1] = 16;
    TB[1] = (amrex::Real *) malloc(16 * sizeof(amrex::Real));
    TBid[1] = (int *) malloc(16 * sizeof(int));
    TBid[1][0] = 2; TB[1][0] = 0.75; // H2
    TBid[1][1] = 5; TB[1][1] = 0.84999999999999998; // O2
    TBid[1][2] = 7; TB[1][2] = 11.890000000000001; // H2O
    TBid[1][3] = 18; TB[1][3] = 0; // C
    TBid[1][4] = 9; TB[1][4] = 1.0900000000000001; // CO
    TBid[1][5] = 19; TB[1][5] = 0; // CH
    TBid[1][6] = 20; TB[1][6] = 0; // HCO
    TBid[1][7] = 21; TB[1][7] = 0; // TXCH2
    TBid[1][8] = 12; TB[1][8] = 2.1800000000000002; // CO2
    TBid[1][9] = 22; TB[1][9] = 0; // SXCH2
    TBid[1][10] = 23; TB[1][10] = 0; // C2H3
    TBid[1][11] = 24; TB[1][11] = 0; // C2H5
    TBid[1][12] = 25; TB[1][12] = 0; // HCCO
    TBid[1][13] = 26; TB[1][13] = 0; // CH3CHO
    TBid[1][14] = 27; TB[1][14] = 0; // CH2CHO
    TBid[1][15] = 28; TB[1][15] = 0; // C2H5O

    // (11):  H + O2 <=> O + OH
    kiv[25] = {1,5,3,4};
    nuv[25] = {-1,-1,1,1};
    // (11):  H + O2 <=> O + OH
    fwd_A[25]     = 26400000000000000;
    fwd_beta[25]  = -0.67000000000000004;
    fwd_Ea[25]    = 17041.110000000001;
    prefactor_units[25]  = 1.0000000000000002e-06;
    activation_units[25] = 0.50321666580471969;
    phase_units[25]      = pow(10,-12.000000);
    is_PD[25] = 0;
    nTB[25] = 0;

    // (12):  H2 + O2 <=> HO2 + H
    kiv[26] = {2,5,8,1};
    nuv[26] = {-1,-1,1,1};
    // (12):  H2 + O2 <=> HO2 + H
    fwd_A[26]     = 592000;
    fwd_beta[26]  = 2.4300000000000002;
    fwd_Ea[26]    = 53501.43;
    prefactor_units[26]  = 1.0000000000000002e-06;
    activation_units[26] = 0.50321666580471969;
    phase_units[26]      = pow(10,-12.000000);
    is_PD[26] = 0;
    nTB[26] = 0;

    // (13):  HO2 + OH <=> H2O + O2
    kiv[27] = {8,4,7,5};
    nuv[27] = {-1,-1,1,1};
    // (13):  HO2 + OH <=> H2O + O2
    fwd_A[27]     = 23800000000000;
    fwd_beta[27]  = 0;
    fwd_Ea[27]    = -499.51999999999998;
    prefactor_units[27]  = 1.0000000000000002e-06;
    activation_units[27] = 0.50321666580471969;
    phase_units[27]      = pow(10,-12.000000);
    is_PD[27] = 0;
    nTB[27] = 0;

    // (14):  HO2 + H <=> 2.000000 OH
    kiv[28] = {8,1,4};
    nuv[28] = {-1,-1,2.0};
    // (14):  HO2 + H <=> 2.000000 OH
    fwd_A[28]     = 74900000000000;
    fwd_beta[28]  = 0;
    fwd_Ea[28]    = 635.75999999999999;
    prefactor_units[28]  = 1.0000000000000002e-06;
    activation_units[28] = 0.50321666580471969;
    phase_units[28]      = pow(10,-12.000000);
    is_PD[28] = 0;
    nTB[28] = 0;

    // (15):  HO2 + O <=> OH + O2
    kiv[29] = {8,3,4,5};
    nuv[29] = {-1,-1,1,1};
    // (15):  HO2 + O <=> OH + O2
    fwd_A[29]     = 40000000000000;
    fwd_beta[29]  = 0;
    fwd_Ea[29]    = 0;
    prefactor_units[29]  = 1.0000000000000002e-06;
    activation_units[29] = 0.50321666580471969;
    phase_units[29]      = pow(10,-12.000000);
    is_PD[29] = 0;
    nTB[29] = 0;

    // (16):  HO2 + H <=> O + H2O
    kiv[30] = {8,1,3,7};
    nuv[30] = {-1,-1,1,1};
    // (16):  HO2 + H <=> O + H2O
    fwd_A[30]     = 3970000000000;
    fwd_beta[30]  = 0;
    fwd_Ea[30]    = 671.61000000000001;
    prefactor_units[30]  = 1.0000000000000002e-06;
    activation_units[30] = 0.50321666580471969;
    phase_units[30]      = pow(10,-12.000000);
    is_PD[30] = 0;
    nTB[30] = 0;

    // (17):  HO2 + OH <=> H2O + O2
    kiv[31] = {8,4,7,5};
    nuv[31] = {-1,-1,1,1};
    // (17):  HO2 + OH <=> H2O + O2
    fwd_A[31]     = 10000000000000000;
    fwd_beta[31]  = 0;
    fwd_Ea[31]    = 17330.310000000001;
    prefactor_units[31]  = 1.0000000000000002e-06;
    activation_units[31] = 0.50321666580471969;
    phase_units[31]      = pow(10,-12.000000);
    is_PD[31] = 0;
    nTB[31] = 0;

    // (18):  H2O2 + O <=> HO2 + OH
    kiv[32] = {6,3,8,4};
    nuv[32] = {-1,-1,1,1};
    // (18):  H2O2 + O <=> HO2 + OH
    fwd_A[32]     = 9630000;
    fwd_beta[32]  = 2;
    fwd_Ea[32]    = 3969.8899999999999;
    prefactor_units[32]  = 1.0000000000000002e-06;
    activation_units[32] = 0.50321666580471969;
    phase_units[32]      = pow(10,-12.000000);
    is_PD[32] = 0;
    nTB[32] = 0;

    // (19):  H2O2 + H <=> H2O + OH
    kiv[33] = {6,1,7,4};
    nuv[33] = {-1,-1,1,1};
    // (19):  H2O2 + H <=> H2O + OH
    fwd_A[33]     = 24100000000000;
    fwd_beta[33]  = 0;
    fwd_Ea[33]    = 3969.8899999999999;
    prefactor_units[33]  = 1.0000000000000002e-06;
    activation_units[33] = 0.50321666580471969;
    phase_units[33]      = pow(10,-12.000000);
    is_PD[33] = 0;
    nTB[33] = 0;

    // (20):  H2O2 + H <=> HO2 + H2
    kiv[34] = {6,1,8,2};
    nuv[34] = {-1,-1,1,1};
    // (20):  H2O2 + H <=> HO2 + H2
    fwd_A[34]     = 6050000;
    fwd_beta[34]  = 2;
    fwd_Ea[34]    = 5200.7600000000002;
    prefactor_units[34]  = 1.0000000000000002e-06;
    activation_units[34] = 0.50321666580471969;
    phase_units[34]      = pow(10,-12.000000);
    is_PD[34] = 0;
    nTB[34] = 0;

    // (21):  H2O2 + OH <=> HO2 + H2O
    kiv[35] = {6,4,8,7};
    nuv[35] = {-1,-1,1,1};
    // (21):  H2O2 + OH <=> HO2 + H2O
    fwd_A[35]     = 2.6700000000000001e+41;
    fwd_beta[35]  = -7;
    fwd_Ea[35]    = 37600.379999999997;
    prefactor_units[35]  = 1.0000000000000002e-06;
    activation_units[35] = 0.50321666580471969;
    phase_units[35]      = pow(10,-12.000000);
    is_PD[35] = 0;
    nTB[35] = 0;

    // (22):  H2O2 + OH <=> HO2 + H2O
    kiv[36] = {6,4,8,7};
    nuv[36] = {-1,-1,1,1};
    // (22):  H2O2 + OH <=> HO2 + H2O
    fwd_A[36]     = 2000000000000;
    fwd_beta[36]  = 0;
    fwd_Ea[36]    = 427.81999999999999;
    prefactor_units[36]  = 1.0000000000000002e-06;
    activation_units[36] = 0.50321666580471969;
    phase_units[36]      = pow(10,-12.000000);
    is_PD[36] = 0;
    nTB[36] = 0;

    // (23):  C + O2 <=> CO + O
    kiv[37] = {18,5,9,3};
    nuv[37] = {-1,-1,1,1};
    // (23):  C + O2 <=> CO + O
    fwd_A[37]     = 58000000000000;
    fwd_beta[37]  = 0;
    fwd_Ea[37]    = 576;
    prefactor_units[37]  = 1.0000000000000002e-06;
    activation_units[37] = 0.50321666580471969;
    phase_units[37]      = pow(10,-12.000000);
    is_PD[37] = 0;
    nTB[37] = 0;

    // (24):  C + OH <=> CO + H
    kiv[38] = {18,4,9,1};
    nuv[38] = {-1,-1,1,1};
    // (24):  C + OH <=> CO + H
    fwd_A[38]     = 50000000000000;
    fwd_beta[38]  = 0;
    fwd_Ea[38]    = 0;
    prefactor_units[38]  = 1.0000000000000002e-06;
    activation_units[38] = 0.50321666580471969;
    phase_units[38]      = pow(10,-12.000000);
    is_PD[38] = 0;
    nTB[38] = 0;

    // (25):  CH + OH <=> HCO + H
    kiv[39] = {19,4,20,1};
    nuv[39] = {-1,-1,1,1};
    // (25):  CH + OH <=> HCO + H
    fwd_A[39]     = 30000000000000;
    fwd_beta[39]  = 0;
    fwd_Ea[39]    = 0;
    prefactor_units[39]  = 1.0000000000000002e-06;
    activation_units[39] = 0.50321666580471969;
    phase_units[39]      = pow(10,-12.000000);
    is_PD[39] = 0;
    nTB[39] = 0;

    // (26):  CH + H2 <=> TXCH2 + H
    kiv[40] = {19,2,21,1};
    nuv[40] = {-1,-1,1,1};
    // (26):  CH + H2 <=> TXCH2 + H
    fwd_A[40]     = 108000000000000;
    fwd_beta[40]  = 0;
    fwd_Ea[40]    = 3109.46;
    prefactor_units[40]  = 1.0000000000000002e-06;
    activation_units[40] = 0.50321666580471969;
    phase_units[40]      = pow(10,-12.000000);
    is_PD[40] = 0;
    nTB[40] = 0;

    // (27):  CH + O <=> CO + H
    kiv[41] = {19,3,9,1};
    nuv[41] = {-1,-1,1,1};
    // (27):  CH + O <=> CO + H
    fwd_A[41]     = 57000000000000;
    fwd_beta[41]  = 0;
    fwd_Ea[41]    = 0;
    prefactor_units[41]  = 1.0000000000000002e-06;
    activation_units[41] = 0.50321666580471969;
    phase_units[41]      = pow(10,-12.000000);
    is_PD[41] = 0;
    nTB[41] = 0;

    // (28):  CH + O2 <=> HCO + O
    kiv[42] = {19,5,20,3};
    nuv[42] = {-1,-1,1,1};
    // (28):  CH + O2 <=> HCO + O
    fwd_A[42]     = 67100000000000;
    fwd_beta[42]  = 0;
    fwd_Ea[42]    = 0;
    prefactor_units[42]  = 1.0000000000000002e-06;
    activation_units[42] = 0.50321666580471969;
    phase_units[42]      = pow(10,-12.000000);
    is_PD[42] = 0;
    nTB[42] = 0;

    // (29):  CH + H <=> C + H2
    kiv[43] = {19,1,18,2};
    nuv[43] = {-1,-1,1,1};
    // (29):  CH + H <=> C + H2
    fwd_A[43]     = 165000000000000;
    fwd_beta[43]  = 0;
    fwd_Ea[43]    = 0;
    prefactor_units[43]  = 1.0000000000000002e-06;
    activation_units[43] = 0.50321666580471969;
    phase_units[43]      = pow(10,-12.000000);
    is_PD[43] = 0;
    nTB[43] = 0;

    // (30):  CH + H2 (+M) <=> CH3 (+M)
    kiv[2] = {19,2,10};
    nuv[2] = {-1,-1,1};
    // (30):  CH + H2 (+M) <=> CH3 (+M)
    fwd_A[2]     = 1970000000000;
    fwd_beta[2]  = 0.42999999999999999;
    fwd_Ea[2]    = -370.45999999999998;
    low_A[2]     = 4.82e+25;
    low_beta[2]  = -2.7999999999999998;
    low_Ea[2]    = 590.34000000000003;
    troe_a[2]    = 0.57799999999999996;
    troe_Tsss[2] = 122;
    troe_Ts[2]   = 2535;
    troe_Tss[2]  = 9365;
    troe_len[2]  = 4;
    prefactor_units[2]  = 1.0000000000000002e-06;
    activation_units[2] = 0.50321666580471969;
    phase_units[2]      = pow(10,-12.000000);
    is_PD[2] = 1;
    nTB[2] = 17;
    TB[2] = (amrex::Real *) malloc(17 * sizeof(amrex::Real));
    TBid[2] = (int *) malloc(17 * sizeof(int));
    TBid[2][0] = 2; TB[2][0] = 2; // H2
    TBid[2][1] = 7; TB[2][1] = 12; // H2O
    TBid[2][2] = 18; TB[2][2] = 0; // C
    TBid[2][3] = 9; TB[2][3] = 1.75; // CO
    TBid[2][4] = 19; TB[2][4] = 0; // CH
    TBid[2][5] = 20; TB[2][5] = 0; // HCO
    TBid[2][6] = 21; TB[2][6] = 0; // TXCH2
    TBid[2][7] = 12; TB[2][7] = 3.6000000000000001; // CO2
    TBid[2][8] = 22; TB[2][8] = 0; // SXCH2
    TBid[2][9] = 13; TB[2][9] = 2; // CH4
    TBid[2][10] = 23; TB[2][10] = 0; // C2H3
    TBid[2][11] = 24; TB[2][11] = 0; // C2H5
    TBid[2][12] = 25; TB[2][12] = 0; // HCCO
    TBid[2][13] = 26; TB[2][13] = 0; // CH3CHO
    TBid[2][14] = 27; TB[2][14] = 0; // CH2CHO
    TBid[2][15] = 28; TB[2][15] = 0; // C2H5O
    TBid[2][16] = 17; TB[2][16] = 3; // C2H6

    // (31):  CH + H2O <=> CH2O + H
    kiv[44] = {19,7,11,1};
    nuv[44] = {-1,-1,1,1};
    // (31):  CH + H2O <=> CH2O + H
    fwd_A[44]     = 5710000000000;
    fwd_beta[44]  = 0;
    fwd_Ea[44]    = -755.25999999999999;
    prefactor_units[44]  = 1.0000000000000002e-06;
    activation_units[44] = 0.50321666580471969;
    phase_units[44]      = pow(10,-12.000000);
    is_PD[44] = 0;
    nTB[44] = 0;

    // (32):  TXCH2 + H (+M) <=> CH3 (+M)
    kiv[3] = {21,1,10};
    nuv[3] = {-1,-1,1};
    // (32):  TXCH2 + H (+M) <=> CH3 (+M)
    fwd_A[3]     = 600000000000000;
    fwd_beta[3]  = 0;
    fwd_Ea[3]    = 0;
    low_A[3]     = 1.0399999999999999e+26;
    low_beta[3]  = -2.7599999999999998;
    low_Ea[3]    = 1598.95;
    troe_a[3]    = 0.56200000000000006;
    troe_Tsss[3] = 91;
    troe_Ts[3]   = 5836;
    troe_Tss[3]  = 8552;
    troe_len[3]  = 4;
    prefactor_units[3]  = 1.0000000000000002e-06;
    activation_units[3] = 0.50321666580471969;
    phase_units[3]      = pow(10,-12.000000);
    is_PD[3] = 1;
    nTB[3] = 17;
    TB[3] = (amrex::Real *) malloc(17 * sizeof(amrex::Real));
    TBid[3] = (int *) malloc(17 * sizeof(int));
    TBid[3][0] = 2; TB[3][0] = 2; // H2
    TBid[3][1] = 7; TB[3][1] = 12; // H2O
    TBid[3][2] = 18; TB[3][2] = 0; // C
    TBid[3][3] = 9; TB[3][3] = 1.75; // CO
    TBid[3][4] = 19; TB[3][4] = 0; // CH
    TBid[3][5] = 20; TB[3][5] = 0; // HCO
    TBid[3][6] = 21; TB[3][6] = 0; // TXCH2
    TBid[3][7] = 12; TB[3][7] = 3.6000000000000001; // CO2
    TBid[3][8] = 22; TB[3][8] = 0; // SXCH2
    TBid[3][9] = 13; TB[3][9] = 2; // CH4
    TBid[3][10] = 23; TB[3][10] = 0; // C2H3
    TBid[3][11] = 24; TB[3][11] = 0; // C2H5
    TBid[3][12] = 25; TB[3][12] = 0; // HCCO
    TBid[3][13] = 26; TB[3][13] = 0; // CH3CHO
    TBid[3][14] = 27; TB[3][14] = 0; // CH2CHO
    TBid[3][15] = 28; TB[3][15] = 0; // C2H5O
    TBid[3][16] = 17; TB[3][16] = 3; // C2H6

    // (33):  TXCH2 + O2 => OH + H + CO
    kiv[45] = {21,5,4,1,9};
    nuv[45] = {-1,-1,1,1,1};
    // (33):  TXCH2 + O2 => OH + H + CO
    fwd_A[45]     = 5000000000000;
    fwd_beta[45]  = 0;
    fwd_Ea[45]    = 1500.96;
    prefactor_units[45]  = 1.0000000000000002e-06;
    activation_units[45] = 0.50321666580471969;
    phase_units[45]      = pow(10,-12.000000);
    is_PD[45] = 0;
    nTB[45] = 0;

    // (34):  TXCH2 + O2 <=> CH2O + O
    kiv[46] = {21,5,11,3};
    nuv[46] = {-1,-1,1,1};
    // (34):  TXCH2 + O2 <=> CH2O + O
    fwd_A[46]     = 2400000000000;
    fwd_beta[46]  = 0;
    fwd_Ea[46]    = 1500.96;
    prefactor_units[46]  = 1.0000000000000002e-06;
    activation_units[46] = 0.50321666580471969;
    phase_units[46]      = pow(10,-12.000000);
    is_PD[46] = 0;
    nTB[46] = 0;

    // (35):  TXCH2 + OH <=> CH2O + H
    kiv[47] = {21,4,11,1};
    nuv[47] = {-1,-1,1,1};
    // (35):  TXCH2 + OH <=> CH2O + H
    fwd_A[47]     = 20000000000000;
    fwd_beta[47]  = 0;
    fwd_Ea[47]    = 0;
    prefactor_units[47]  = 1.0000000000000002e-06;
    activation_units[47] = 0.50321666580471969;
    phase_units[47]      = pow(10,-12.000000);
    is_PD[47] = 0;
    nTB[47] = 0;

    // (36):  TXCH2 + HO2 <=> CH2O + OH
    kiv[48] = {21,8,11,4};
    nuv[48] = {-1,-1,1,1};
    // (36):  TXCH2 + HO2 <=> CH2O + OH
    fwd_A[48]     = 20000000000000;
    fwd_beta[48]  = 0;
    fwd_Ea[48]    = 0;
    prefactor_units[48]  = 1.0000000000000002e-06;
    activation_units[48] = 0.50321666580471969;
    phase_units[48]      = pow(10,-12.000000);
    is_PD[48] = 0;
    nTB[48] = 0;

    // (37):  TXCH2 + O2 => CO2 + 2.000000 H
    kiv[49] = {21,5,12,1};
    nuv[49] = {-1,-1,1,2.0};
    // (37):  TXCH2 + O2 => CO2 + 2.000000 H
    fwd_A[49]     = 5800000000000;
    fwd_beta[49]  = 0;
    fwd_Ea[49]    = 1500.96;
    prefactor_units[49]  = 1.0000000000000002e-06;
    activation_units[49] = 0.50321666580471969;
    phase_units[49]      = pow(10,-12.000000);
    is_PD[49] = 0;
    nTB[49] = 0;

    // (38):  TXCH2 + OH <=> CH + H2O
    kiv[50] = {21,4,19,7};
    nuv[50] = {-1,-1,1,1};
    // (38):  TXCH2 + OH <=> CH + H2O
    fwd_A[50]     = 11300000;
    fwd_beta[50]  = 2;
    fwd_Ea[50]    = 2999.52;
    prefactor_units[50]  = 1.0000000000000002e-06;
    activation_units[50] = 0.50321666580471969;
    phase_units[50]      = pow(10,-12.000000);
    is_PD[50] = 0;
    nTB[50] = 0;

    // (39):  TXCH2 + O <=> HCO + H
    kiv[51] = {21,3,20,1};
    nuv[51] = {-1,-1,1,1};
    // (39):  TXCH2 + O <=> HCO + H
    fwd_A[51]     = 80000000000000;
    fwd_beta[51]  = 0;
    fwd_Ea[51]    = 0;
    prefactor_units[51]  = 1.0000000000000002e-06;
    activation_units[51] = 0.50321666580471969;
    phase_units[51]      = pow(10,-12.000000);
    is_PD[51] = 0;
    nTB[51] = 0;

    // (40):  TXCH2 + H2 <=> H + CH3
    kiv[52] = {21,2,1,10};
    nuv[52] = {-1,-1,1,1};
    // (40):  TXCH2 + H2 <=> H + CH3
    fwd_A[52]     = 500000;
    fwd_beta[52]  = 2;
    fwd_Ea[52]    = 7229.9200000000001;
    prefactor_units[52]  = 1.0000000000000002e-06;
    activation_units[52] = 0.50321666580471969;
    phase_units[52]      = pow(10,-12.000000);
    is_PD[52] = 0;
    nTB[52] = 0;

    // (41):  SXCH2 + H2O <=> TXCH2 + H2O
    kiv[53] = {22,7,21,7};
    nuv[53] = {-1,-1,1,1};
    // (41):  SXCH2 + H2O <=> TXCH2 + H2O
    fwd_A[53]     = 30000000000000;
    fwd_beta[53]  = 0;
    fwd_Ea[53]    = 0;
    prefactor_units[53]  = 1.0000000000000002e-06;
    activation_units[53] = 0.50321666580471969;
    phase_units[53]      = pow(10,-12.000000);
    is_PD[53] = 0;
    nTB[53] = 0;

    // (42):  SXCH2 + H <=> CH + H2
    kiv[54] = {22,1,19,2};
    nuv[54] = {-1,-1,1,1};
    // (42):  SXCH2 + H <=> CH + H2
    fwd_A[54]     = 30000000000000;
    fwd_beta[54]  = 0;
    fwd_Ea[54]    = 0;
    prefactor_units[54]  = 1.0000000000000002e-06;
    activation_units[54] = 0.50321666580471969;
    phase_units[54]      = pow(10,-12.000000);
    is_PD[54] = 0;
    nTB[54] = 0;

    // (43):  SXCH2 + O2 <=> H + OH + CO
    kiv[55] = {22,5,1,4,9};
    nuv[55] = {-1,-1,1,1,1};
    // (43):  SXCH2 + O2 <=> H + OH + CO
    fwd_A[55]     = 28000000000000;
    fwd_beta[55]  = 0;
    fwd_Ea[55]    = 0;
    prefactor_units[55]  = 1.0000000000000002e-06;
    activation_units[55] = 0.50321666580471969;
    phase_units[55]      = pow(10,-12.000000);
    is_PD[55] = 0;
    nTB[55] = 0;

    // (44):  SXCH2 + O <=> CO + H2
    kiv[56] = {22,3,9,2};
    nuv[56] = {-1,-1,1,1};
    // (44):  SXCH2 + O <=> CO + H2
    fwd_A[56]     = 15000000000000;
    fwd_beta[56]  = 0;
    fwd_Ea[56]    = 0;
    prefactor_units[56]  = 1.0000000000000002e-06;
    activation_units[56] = 0.50321666580471969;
    phase_units[56]      = pow(10,-12.000000);
    is_PD[56] = 0;
    nTB[56] = 0;

    // (45):  SXCH2 + O2 <=> CO + H2O
    kiv[57] = {22,5,9,7};
    nuv[57] = {-1,-1,1,1};
    // (45):  SXCH2 + O2 <=> CO + H2O
    fwd_A[57]     = 12000000000000;
    fwd_beta[57]  = 0;
    fwd_Ea[57]    = 0;
    prefactor_units[57]  = 1.0000000000000002e-06;
    activation_units[57] = 0.50321666580471969;
    phase_units[57]      = pow(10,-12.000000);
    is_PD[57] = 0;
    nTB[57] = 0;

    // (46):  SXCH2 + H2 <=> CH3 + H
    kiv[58] = {22,2,10,1};
    nuv[58] = {-1,-1,1,1};
    // (46):  SXCH2 + H2 <=> CH3 + H
    fwd_A[58]     = 70000000000000;
    fwd_beta[58]  = 0;
    fwd_Ea[58]    = 0;
    prefactor_units[58]  = 1.0000000000000002e-06;
    activation_units[58] = 0.50321666580471969;
    phase_units[58]      = pow(10,-12.000000);
    is_PD[58] = 0;
    nTB[58] = 0;

    // (47):  SXCH2 + O <=> HCO + H
    kiv[59] = {22,3,20,1};
    nuv[59] = {-1,-1,1,1};
    // (47):  SXCH2 + O <=> HCO + H
    fwd_A[59]     = 15000000000000;
    fwd_beta[59]  = 0;
    fwd_Ea[59]    = 0;
    prefactor_units[59]  = 1.0000000000000002e-06;
    activation_units[59] = 0.50321666580471969;
    phase_units[59]      = pow(10,-12.000000);
    is_PD[59] = 0;
    nTB[59] = 0;

    // (48):  SXCH2 + H2O => H2 + CH2O
    kiv[60] = {22,7,2,11};
    nuv[60] = {-1,-1,1,1};
    // (48):  SXCH2 + H2O => H2 + CH2O
    fwd_A[60]     = 68200000000;
    fwd_beta[60]  = 0.25;
    fwd_Ea[60]    = -934.50999999999999;
    prefactor_units[60]  = 1.0000000000000002e-06;
    activation_units[60] = 0.50321666580471969;
    phase_units[60]      = pow(10,-12.000000);
    is_PD[60] = 0;
    nTB[60] = 0;

    // (49):  SXCH2 + OH <=> CH2O + H
    kiv[61] = {22,4,11,1};
    nuv[61] = {-1,-1,1,1};
    // (49):  SXCH2 + OH <=> CH2O + H
    fwd_A[61]     = 30000000000000;
    fwd_beta[61]  = 0;
    fwd_Ea[61]    = 0;
    prefactor_units[61]  = 1.0000000000000002e-06;
    activation_units[61] = 0.50321666580471969;
    phase_units[61]      = pow(10,-12.000000);
    is_PD[61] = 0;
    nTB[61] = 0;

    // (50):  CH3 + OH => H2 + CH2O
    kiv[62] = {10,4,2,11};
    nuv[62] = {-1,-1,1,1};
    // (50):  CH3 + OH => H2 + CH2O
    fwd_A[62]     = 8000000000;
    fwd_beta[62]  = 0;
    fwd_Ea[62]    = -1754.3;
    prefactor_units[62]  = 1.0000000000000002e-06;
    activation_units[62] = 0.50321666580471969;
    phase_units[62]      = pow(10,-12.000000);
    is_PD[62] = 0;
    nTB[62] = 0;

    // (51):  CH3 + H2O2 <=> CH4 + HO2
    kiv[63] = {10,6,13,8};
    nuv[63] = {-1,-1,1,1};
    // (51):  CH3 + H2O2 <=> CH4 + HO2
    fwd_A[63]     = 24500;
    fwd_beta[63]  = 2.4700000000000002;
    fwd_Ea[63]    = 5179.25;
    prefactor_units[63]  = 1.0000000000000002e-06;
    activation_units[63] = 0.50321666580471969;
    phase_units[63]      = pow(10,-12.000000);
    is_PD[63] = 0;
    nTB[63] = 0;

    // (52):  CH3 + O2 <=> CH2O + OH
    kiv[64] = {10,5,11,4};
    nuv[64] = {-1,-1,1,1};
    // (52):  CH3 + O2 <=> CH2O + OH
    fwd_A[64]     = 587000000000;
    fwd_beta[64]  = 0;
    fwd_Ea[64]    = 13840.82;
    prefactor_units[64]  = 1.0000000000000002e-06;
    activation_units[64] = 0.50321666580471969;
    phase_units[64]      = pow(10,-12.000000);
    is_PD[64] = 0;
    nTB[64] = 0;

    // (53):  CH3 + CH <=> C2H3 + H
    kiv[65] = {10,19,23,1};
    nuv[65] = {-1,-1,1,1};
    // (53):  CH3 + CH <=> C2H3 + H
    fwd_A[65]     = 30000000000000;
    fwd_beta[65]  = 0;
    fwd_Ea[65]    = 0;
    prefactor_units[65]  = 1.0000000000000002e-06;
    activation_units[65] = 0.50321666580471969;
    phase_units[65]      = pow(10,-12.000000);
    is_PD[65] = 0;
    nTB[65] = 0;

    // (54):  CH3 + O <=> CH2O + H
    kiv[66] = {10,3,11,1};
    nuv[66] = {-1,-1,1,1};
    // (54):  CH3 + O <=> CH2O + H
    fwd_A[66]     = 50600000000000;
    fwd_beta[66]  = 0;
    fwd_Ea[66]    = 0;
    prefactor_units[66]  = 1.0000000000000002e-06;
    activation_units[66] = 0.50321666580471969;
    phase_units[66]      = pow(10,-12.000000);
    is_PD[66] = 0;
    nTB[66] = 0;

    // (55):  CH3 + C <=> C2H2 + H
    kiv[67] = {10,18,14,1};
    nuv[67] = {-1,-1,1,1};
    // (55):  CH3 + C <=> C2H2 + H
    fwd_A[67]     = 50000000000000;
    fwd_beta[67]  = 0;
    fwd_Ea[67]    = 0;
    prefactor_units[67]  = 1.0000000000000002e-06;
    activation_units[67] = 0.50321666580471969;
    phase_units[67]      = pow(10,-12.000000);
    is_PD[67] = 0;
    nTB[67] = 0;

    // (56):  CH3 + H (+M) <=> CH4 (+M)
    kiv[4] = {10,1,13};
    nuv[4] = {-1,-1,1};
    // (56):  CH3 + H (+M) <=> CH4 (+M)
    fwd_A[4]     = 69200000000000;
    fwd_beta[4]  = 0.17999999999999999;
    fwd_Ea[4]    = 0;
    low_A[4]     = 3.4700000000000003e+38;
    low_beta[4]  = -6.2999999999999998;
    low_Ea[4]    = 5074.0900000000001;
    troe_a[4]    = 0.78300000000000003;
    troe_Tsss[4] = 74;
    troe_Ts[4]   = 2941;
    troe_Tss[4]  = 6964;
    troe_len[4]  = 4;
    prefactor_units[4]  = 1.0000000000000002e-06;
    activation_units[4] = 0.50321666580471969;
    phase_units[4]      = pow(10,-12.000000);
    is_PD[4] = 1;
    nTB[4] = 17;
    TB[4] = (amrex::Real *) malloc(17 * sizeof(amrex::Real));
    TBid[4] = (int *) malloc(17 * sizeof(int));
    TBid[4][0] = 2; TB[4][0] = 2; // H2
    TBid[4][1] = 7; TB[4][1] = 6; // H2O
    TBid[4][2] = 18; TB[4][2] = 0; // C
    TBid[4][3] = 9; TB[4][3] = 1.5; // CO
    TBid[4][4] = 19; TB[4][4] = 0; // CH
    TBid[4][5] = 20; TB[4][5] = 0; // HCO
    TBid[4][6] = 21; TB[4][6] = 0; // TXCH2
    TBid[4][7] = 12; TB[4][7] = 2; // CO2
    TBid[4][8] = 22; TB[4][8] = 0; // SXCH2
    TBid[4][9] = 13; TB[4][9] = 3; // CH4
    TBid[4][10] = 23; TB[4][10] = 0; // C2H3
    TBid[4][11] = 24; TB[4][11] = 0; // C2H5
    TBid[4][12] = 25; TB[4][12] = 0; // HCCO
    TBid[4][13] = 26; TB[4][13] = 0; // CH3CHO
    TBid[4][14] = 27; TB[4][14] = 0; // CH2CHO
    TBid[4][15] = 28; TB[4][15] = 0; // C2H5O
    TBid[4][16] = 17; TB[4][16] = 3; // C2H6

    // (57):  CH3 + OH <=> TXCH2 + H2O
    kiv[68] = {10,4,21,7};
    nuv[68] = {-1,-1,1,1};
    // (57):  CH3 + OH <=> TXCH2 + H2O
    fwd_A[68]     = 56000000;
    fwd_beta[68]  = 1.6000000000000001;
    fwd_Ea[68]    = 5420.6499999999996;
    prefactor_units[68]  = 1.0000000000000002e-06;
    activation_units[68] = 0.50321666580471969;
    phase_units[68]      = pow(10,-12.000000);
    is_PD[68] = 0;
    nTB[68] = 0;

    // (58):  CH3 + SXCH2 <=> C2H4 + H
    kiv[69] = {10,22,15,1};
    nuv[69] = {-1,-1,1,1};
    // (58):  CH3 + SXCH2 <=> C2H4 + H
    fwd_A[69]     = 12000000000000;
    fwd_beta[69]  = 0;
    fwd_Ea[69]    = -571.22000000000003;
    prefactor_units[69]  = 1.0000000000000002e-06;
    activation_units[69] = 0.50321666580471969;
    phase_units[69]      = pow(10,-12.000000);
    is_PD[69] = 0;
    nTB[69] = 0;

    // (59):  CH3 + OH <=> SXCH2 + H2O
    kiv[70] = {10,4,22,7};
    nuv[70] = {-1,-1,1,1};
    // (59):  CH3 + OH <=> SXCH2 + H2O
    fwd_A[70]     = 6.44e+17;
    fwd_beta[70]  = -1.3400000000000001;
    fwd_Ea[70]    = 1417.3;
    prefactor_units[70]  = 1.0000000000000002e-06;
    activation_units[70] = 0.50321666580471969;
    phase_units[70]      = pow(10,-12.000000);
    is_PD[70] = 0;
    nTB[70] = 0;

    // (60):  2.000000 CH3 <=> C2H5 + H
    kiv[71] = {10,24,1};
    nuv[71] = {-2.0,1,1};
    // (60):  2.000000 CH3 <=> C2H5 + H
    fwd_A[71]     = 6840000000000;
    fwd_beta[71]  = 0.10000000000000001;
    fwd_Ea[71]    = 10599.9;
    prefactor_units[71]  = 1.0000000000000002e-06;
    activation_units[71] = 0.50321666580471969;
    phase_units[71]      = pow(10,-12.000000);
    is_PD[71] = 0;
    nTB[71] = 0;

    // (61):  CH3 + HO2 <=> CH4 + O2
    kiv[72] = {10,8,13,5};
    nuv[72] = {-1,-1,1,1};
    // (61):  CH3 + HO2 <=> CH4 + O2
    fwd_A[72]     = 3610000000000;
    fwd_beta[72]  = 0;
    fwd_Ea[72]    = 0;
    prefactor_units[72]  = 1.0000000000000002e-06;
    activation_units[72] = 0.50321666580471969;
    phase_units[72]      = pow(10,-12.000000);
    is_PD[72] = 0;
    nTB[72] = 0;

    // (62):  CH3 + TXCH2 <=> C2H4 + H
    kiv[73] = {10,21,15,1};
    nuv[73] = {-1,-1,1,1};
    // (62):  CH3 + TXCH2 <=> C2H4 + H
    fwd_A[73]     = 100000000000000;
    fwd_beta[73]  = 0;
    fwd_Ea[73]    = 0;
    prefactor_units[73]  = 1.0000000000000002e-06;
    activation_units[73] = 0.50321666580471969;
    phase_units[73]      = pow(10,-12.000000);
    is_PD[73] = 0;
    nTB[73] = 0;

    // (63):  CH3 + O => H + H2 + CO
    kiv[74] = {10,3,1,2,9};
    nuv[74] = {-1,-1,1,1,1};
    // (63):  CH3 + O => H + H2 + CO
    fwd_A[74]     = 33700000000000;
    fwd_beta[74]  = 0;
    fwd_Ea[74]    = 0;
    prefactor_units[74]  = 1.0000000000000002e-06;
    activation_units[74] = 0.50321666580471969;
    phase_units[74]      = pow(10,-12.000000);
    is_PD[74] = 0;
    nTB[74] = 0;

    // (64):  CH4 + CH <=> C2H4 + H
    kiv[75] = {13,19,15,1};
    nuv[75] = {-1,-1,1,1};
    // (64):  CH4 + CH <=> C2H4 + H
    fwd_A[75]     = 60000000000000;
    fwd_beta[75]  = 0;
    fwd_Ea[75]    = 0;
    prefactor_units[75]  = 1.0000000000000002e-06;
    activation_units[75] = 0.50321666580471969;
    phase_units[75]      = pow(10,-12.000000);
    is_PD[75] = 0;
    nTB[75] = 0;

    // (65):  CH4 + SXCH2 <=> 2.000000 CH3
    kiv[76] = {13,22,10};
    nuv[76] = {-1,-1,2.0};
    // (65):  CH4 + SXCH2 <=> 2.000000 CH3
    fwd_A[76]     = 16000000000000;
    fwd_beta[76]  = 0;
    fwd_Ea[76]    = -571.22000000000003;
    prefactor_units[76]  = 1.0000000000000002e-06;
    activation_units[76] = 0.50321666580471969;
    phase_units[76]      = pow(10,-12.000000);
    is_PD[76] = 0;
    nTB[76] = 0;

    // (66):  CH4 + O <=> CH3 + OH
    kiv[77] = {13,3,10,4};
    nuv[77] = {-1,-1,1,1};
    // (66):  CH4 + O <=> CH3 + OH
    fwd_A[77]     = 1020000000;
    fwd_beta[77]  = 1.5;
    fwd_Ea[77]    = 8599.4300000000003;
    prefactor_units[77]  = 1.0000000000000002e-06;
    activation_units[77] = 0.50321666580471969;
    phase_units[77]      = pow(10,-12.000000);
    is_PD[77] = 0;
    nTB[77] = 0;

    // (67):  CH4 + OH <=> CH3 + H2O
    kiv[78] = {13,4,10,7};
    nuv[78] = {-1,-1,1,1};
    // (67):  CH4 + OH <=> CH3 + H2O
    fwd_A[78]     = 100000000;
    fwd_beta[78]  = 1.6000000000000001;
    fwd_Ea[78]    = 3119.02;
    prefactor_units[78]  = 1.0000000000000002e-06;
    activation_units[78] = 0.50321666580471969;
    phase_units[78]      = pow(10,-12.000000);
    is_PD[78] = 0;
    nTB[78] = 0;

    // (68):  CH4 + TXCH2 <=> 2.000000 CH3
    kiv[79] = {13,21,10};
    nuv[79] = {-1,-1,2.0};
    // (68):  CH4 + TXCH2 <=> 2.000000 CH3
    fwd_A[79]     = 2460000;
    fwd_beta[79]  = 2;
    fwd_Ea[79]    = 8269.6000000000004;
    prefactor_units[79]  = 1.0000000000000002e-06;
    activation_units[79] = 0.50321666580471969;
    phase_units[79]      = pow(10,-12.000000);
    is_PD[79] = 0;
    nTB[79] = 0;

    // (69):  CH4 + H <=> CH3 + H2
    kiv[80] = {13,1,10,2};
    nuv[80] = {-1,-1,1,1};
    // (69):  CH4 + H <=> CH3 + H2
    fwd_A[80]     = 660000000;
    fwd_beta[80]  = 1.6200000000000001;
    fwd_Ea[80]    = 10841.299999999999;
    prefactor_units[80]  = 1.0000000000000002e-06;
    activation_units[80] = 0.50321666580471969;
    phase_units[80]      = pow(10,-12.000000);
    is_PD[80] = 0;
    nTB[80] = 0;

    // (70):  TXCH2 + CO (+M) <=> CH2CO (+M)
    kiv[5] = {21,9,16};
    nuv[5] = {-1,-1,1};
    // (70):  TXCH2 + CO (+M) <=> CH2CO (+M)
    fwd_A[5]     = 810000000000;
    fwd_beta[5]  = 0.5;
    fwd_Ea[5]    = 4510.04;
    low_A[5]     = 2.69e+33;
    low_beta[5]  = -5.1100000000000003;
    low_Ea[5]    = 7096.0799999999999;
    troe_a[5]    = 0.5907;
    troe_Tsss[5] = 275;
    troe_Ts[5]   = 1226;
    troe_Tss[5]  = 5185;
    troe_len[5]  = 4;
    prefactor_units[5]  = 1.0000000000000002e-06;
    activation_units[5] = 0.50321666580471969;
    phase_units[5]      = pow(10,-12.000000);
    is_PD[5] = 1;
    nTB[5] = 17;
    TB[5] = (amrex::Real *) malloc(17 * sizeof(amrex::Real));
    TBid[5] = (int *) malloc(17 * sizeof(int));
    TBid[5][0] = 2; TB[5][0] = 2; // H2
    TBid[5][1] = 7; TB[5][1] = 12; // H2O
    TBid[5][2] = 18; TB[5][2] = 0; // C
    TBid[5][3] = 9; TB[5][3] = 1.75; // CO
    TBid[5][4] = 19; TB[5][4] = 0; // CH
    TBid[5][5] = 20; TB[5][5] = 0; // HCO
    TBid[5][6] = 21; TB[5][6] = 0; // TXCH2
    TBid[5][7] = 12; TB[5][7] = 3.6000000000000001; // CO2
    TBid[5][8] = 22; TB[5][8] = 0; // SXCH2
    TBid[5][9] = 13; TB[5][9] = 2; // CH4
    TBid[5][10] = 23; TB[5][10] = 0; // C2H3
    TBid[5][11] = 24; TB[5][11] = 0; // C2H5
    TBid[5][12] = 25; TB[5][12] = 0; // HCCO
    TBid[5][13] = 26; TB[5][13] = 0; // CH3CHO
    TBid[5][14] = 27; TB[5][14] = 0; // CH2CHO
    TBid[5][15] = 28; TB[5][15] = 0; // C2H5O
    TBid[5][16] = 17; TB[5][16] = 3; // C2H6

    // (71):  SXCH2 + CO <=> TXCH2 + CO
    kiv[81] = {22,9,21,9};
    nuv[81] = {-1,-1,1,1};
    // (71):  SXCH2 + CO <=> TXCH2 + CO
    fwd_A[81]     = 9000000000000;
    fwd_beta[81]  = 0;
    fwd_Ea[81]    = 0;
    prefactor_units[81]  = 1.0000000000000002e-06;
    activation_units[81] = 0.50321666580471969;
    phase_units[81]      = pow(10,-12.000000);
    is_PD[81] = 0;
    nTB[81] = 0;

    // (72):  CO + O2 <=> CO2 + O
    kiv[82] = {9,5,12,3};
    nuv[82] = {-1,-1,1,1};
    // (72):  CO + O2 <=> CO2 + O
    fwd_A[82]     = 1120000000000;
    fwd_beta[82]  = 0;
    fwd_Ea[82]    = 47700.760000000002;
    prefactor_units[82]  = 1.0000000000000002e-06;
    activation_units[82] = 0.50321666580471969;
    phase_units[82]      = pow(10,-12.000000);
    is_PD[82] = 0;
    nTB[82] = 0;

    // (73):  CO + OH <=> CO2 + H
    kiv[83] = {9,4,12,1};
    nuv[83] = {-1,-1,1,1};
    // (73):  CO + OH <=> CO2 + H
    fwd_A[83]     = 87800000000;
    fwd_beta[83]  = 0.029999999999999999;
    fwd_Ea[83]    = -16.73;
    prefactor_units[83]  = 1.0000000000000002e-06;
    activation_units[83] = 0.50321666580471969;
    phase_units[83]      = pow(10,-12.000000);
    is_PD[83] = 0;
    nTB[83] = 0;

    // (74):  CO + H2 (+M) <=> CH2O (+M)
    kiv[6] = {9,2,11};
    nuv[6] = {-1,-1,1};
    // (74):  CO + H2 (+M) <=> CH2O (+M)
    fwd_A[6]     = 43000000;
    fwd_beta[6]  = 1.5;
    fwd_Ea[6]    = 79600.860000000001;
    low_A[6]     = 5.0699999999999998e+27;
    low_beta[6]  = -3.4199999999999999;
    low_Ea[6]    = 84349.899999999994;
    troe_a[6]    = 0.93200000000000005;
    troe_Tsss[6] = 197;
    troe_Ts[6]   = 1540;
    troe_Tss[6]  = 10300;
    troe_len[6]  = 4;
    prefactor_units[6]  = 1.0000000000000002e-06;
    activation_units[6] = 0.50321666580471969;
    phase_units[6]      = pow(10,-12.000000);
    is_PD[6] = 1;
    nTB[6] = 17;
    TB[6] = (amrex::Real *) malloc(17 * sizeof(amrex::Real));
    TBid[6] = (int *) malloc(17 * sizeof(int));
    TBid[6][0] = 2; TB[6][0] = 2; // H2
    TBid[6][1] = 7; TB[6][1] = 12; // H2O
    TBid[6][2] = 18; TB[6][2] = 0; // C
    TBid[6][3] = 9; TB[6][3] = 1.75; // CO
    TBid[6][4] = 19; TB[6][4] = 0; // CH
    TBid[6][5] = 20; TB[6][5] = 0; // HCO
    TBid[6][6] = 21; TB[6][6] = 0; // TXCH2
    TBid[6][7] = 12; TB[6][7] = 3.6000000000000001; // CO2
    TBid[6][8] = 22; TB[6][8] = 0; // SXCH2
    TBid[6][9] = 13; TB[6][9] = 2; // CH4
    TBid[6][10] = 23; TB[6][10] = 0; // C2H3
    TBid[6][11] = 24; TB[6][11] = 0; // C2H5
    TBid[6][12] = 25; TB[6][12] = 0; // HCCO
    TBid[6][13] = 26; TB[6][13] = 0; // CH3CHO
    TBid[6][14] = 27; TB[6][14] = 0; // CH2CHO
    TBid[6][15] = 28; TB[6][15] = 0; // C2H5O
    TBid[6][16] = 17; TB[6][16] = 3; // C2H6

    // (75):  CH + CO (+M) <=> HCCO (+M)
    kiv[7] = {19,9,25};
    nuv[7] = {-1,-1,1};
    // (75):  CH + CO (+M) <=> HCCO (+M)
    fwd_A[7]     = 50000000000000;
    fwd_beta[7]  = 0;
    fwd_Ea[7]    = 0;
    low_A[7]     = 2.6899999999999998e+28;
    low_beta[7]  = -3.7400000000000002;
    low_Ea[7]    = 1935.95;
    troe_a[7]    = 0.57569999999999999;
    troe_Tsss[7] = 237;
    troe_Ts[7]   = 1652;
    troe_Tss[7]  = 5069;
    troe_len[7]  = 4;
    prefactor_units[7]  = 1.0000000000000002e-06;
    activation_units[7] = 0.50321666580471969;
    phase_units[7]      = pow(10,-12.000000);
    is_PD[7] = 1;
    nTB[7] = 17;
    TB[7] = (amrex::Real *) malloc(17 * sizeof(amrex::Real));
    TBid[7] = (int *) malloc(17 * sizeof(int));
    TBid[7][0] = 2; TB[7][0] = 2; // H2
    TBid[7][1] = 7; TB[7][1] = 12; // H2O
    TBid[7][2] = 18; TB[7][2] = 0; // C
    TBid[7][3] = 9; TB[7][3] = 1.75; // CO
    TBid[7][4] = 19; TB[7][4] = 0; // CH
    TBid[7][5] = 20; TB[7][5] = 0; // HCO
    TBid[7][6] = 21; TB[7][6] = 0; // TXCH2
    TBid[7][7] = 12; TB[7][7] = 3.6000000000000001; // CO2
    TBid[7][8] = 22; TB[7][8] = 0; // SXCH2
    TBid[7][9] = 13; TB[7][9] = 2; // CH4
    TBid[7][10] = 23; TB[7][10] = 0; // C2H3
    TBid[7][11] = 24; TB[7][11] = 0; // C2H5
    TBid[7][12] = 25; TB[7][12] = 0; // HCCO
    TBid[7][13] = 26; TB[7][13] = 0; // CH3CHO
    TBid[7][14] = 27; TB[7][14] = 0; // CH2CHO
    TBid[7][15] = 28; TB[7][15] = 0; // C2H5O
    TBid[7][16] = 17; TB[7][16] = 3; // C2H6

    // (76):  CO + OH <=> CO2 + H
    kiv[84] = {9,4,12,1};
    nuv[84] = {-1,-1,1,1};
    // (76):  CO + OH <=> CO2 + H
    fwd_A[84]     = 800000000000;
    fwd_beta[84]  = 0.14000000000000001;
    fwd_Ea[84]    = 7351.8199999999997;
    prefactor_units[84]  = 1.0000000000000002e-06;
    activation_units[84] = 0.50321666580471969;
    phase_units[84]      = pow(10,-12.000000);
    is_PD[84] = 0;
    nTB[84] = 0;

    // (77):  CO + O (+M) <=> CO2 (+M)
    kiv[8] = {9,3,12};
    nuv[8] = {-1,-1,1};
    // (77):  CO + O (+M) <=> CO2 (+M)
    fwd_A[8]     = 13600000000;
    fwd_beta[8]  = 0;
    fwd_Ea[8]    = 2385.2800000000002;
    low_A[8]     = 1.17e+24;
    low_beta[8]  = -2.79;
    low_Ea[8]    = 4192.1599999999999;
    troe_a[8]    = 1;
    troe_Tsss[8] = 1;
    troe_Ts[8]   = 10000000;
    troe_Tss[8]  = 10000000;
    troe_len[8]  = 4;
    prefactor_units[8]  = 1.0000000000000002e-06;
    activation_units[8] = 0.50321666580471969;
    phase_units[8]      = pow(10,-12.000000);
    is_PD[8] = 1;
    nTB[8] = 17;
    TB[8] = (amrex::Real *) malloc(17 * sizeof(amrex::Real));
    TBid[8] = (int *) malloc(17 * sizeof(int));
    TBid[8][0] = 2; TB[8][0] = 2; // H2
    TBid[8][1] = 7; TB[8][1] = 12; // H2O
    TBid[8][2] = 18; TB[8][2] = 0; // C
    TBid[8][3] = 9; TB[8][3] = 1.75; // CO
    TBid[8][4] = 19; TB[8][4] = 0; // CH
    TBid[8][5] = 20; TB[8][5] = 0; // HCO
    TBid[8][6] = 21; TB[8][6] = 0; // TXCH2
    TBid[8][7] = 12; TB[8][7] = 3.6000000000000001; // CO2
    TBid[8][8] = 22; TB[8][8] = 0; // SXCH2
    TBid[8][9] = 13; TB[8][9] = 2; // CH4
    TBid[8][10] = 23; TB[8][10] = 0; // C2H3
    TBid[8][11] = 24; TB[8][11] = 0; // C2H5
    TBid[8][12] = 25; TB[8][12] = 0; // HCCO
    TBid[8][13] = 26; TB[8][13] = 0; // CH3CHO
    TBid[8][14] = 27; TB[8][14] = 0; // CH2CHO
    TBid[8][15] = 28; TB[8][15] = 0; // C2H5O
    TBid[8][16] = 17; TB[8][16] = 3; // C2H6

    // (78):  CO + HO2 <=> CO2 + OH
    kiv[85] = {9,8,12,4};
    nuv[85] = {-1,-1,1,1};
    // (78):  CO + HO2 <=> CO2 + OH
    fwd_A[85]     = 30100000000000;
    fwd_beta[85]  = 0;
    fwd_Ea[85]    = 22999.52;
    prefactor_units[85]  = 1.0000000000000002e-06;
    activation_units[85] = 0.50321666580471969;
    phase_units[85]      = pow(10,-12.000000);
    is_PD[85] = 0;
    nTB[85] = 0;

    // (79):  HCO + H <=> CO + H2
    kiv[86] = {20,1,9,2};
    nuv[86] = {-1,-1,1,1};
    // (79):  HCO + H <=> CO + H2
    fwd_A[86]     = 120000000000000;
    fwd_beta[86]  = 0;
    fwd_Ea[86]    = 0;
    prefactor_units[86]  = 1.0000000000000002e-06;
    activation_units[86] = 0.50321666580471969;
    phase_units[86]      = pow(10,-12.000000);
    is_PD[86] = 0;
    nTB[86] = 0;

    // (80):  HCO + H (+M) <=> CH2O (+M)
    kiv[9] = {20,1,11};
    nuv[9] = {-1,-1,1};
    // (80):  HCO + H (+M) <=> CH2O (+M)
    fwd_A[9]     = 1090000000000;
    fwd_beta[9]  = 0.47999999999999998;
    fwd_Ea[9]    = -260.51999999999998;
    low_A[9]     = 2.4700000000000001e+24;
    low_beta[9]  = -2.5699999999999998;
    low_Ea[9]    = 425.43000000000001;
    troe_a[9]    = 0.78239999999999998;
    troe_Tsss[9] = 271;
    troe_Ts[9]   = 2755;
    troe_Tss[9]  = 6570;
    troe_len[9]  = 4;
    prefactor_units[9]  = 1.0000000000000002e-06;
    activation_units[9] = 0.50321666580471969;
    phase_units[9]      = pow(10,-12.000000);
    is_PD[9] = 1;
    nTB[9] = 17;
    TB[9] = (amrex::Real *) malloc(17 * sizeof(amrex::Real));
    TBid[9] = (int *) malloc(17 * sizeof(int));
    TBid[9][0] = 2; TB[9][0] = 2; // H2
    TBid[9][1] = 7; TB[9][1] = 12; // H2O
    TBid[9][2] = 18; TB[9][2] = 0; // C
    TBid[9][3] = 9; TB[9][3] = 1.75; // CO
    TBid[9][4] = 19; TB[9][4] = 0; // CH
    TBid[9][5] = 20; TB[9][5] = 0; // HCO
    TBid[9][6] = 21; TB[9][6] = 0; // TXCH2
    TBid[9][7] = 12; TB[9][7] = 3.6000000000000001; // CO2
    TBid[9][8] = 22; TB[9][8] = 0; // SXCH2
    TBid[9][9] = 13; TB[9][9] = 2; // CH4
    TBid[9][10] = 23; TB[9][10] = 0; // C2H3
    TBid[9][11] = 24; TB[9][11] = 0; // C2H5
    TBid[9][12] = 25; TB[9][12] = 0; // HCCO
    TBid[9][13] = 26; TB[9][13] = 0; // CH3CHO
    TBid[9][14] = 27; TB[9][14] = 0; // CH2CHO
    TBid[9][15] = 28; TB[9][15] = 0; // C2H5O
    TBid[9][16] = 17; TB[9][16] = 3; // C2H6

    // (81):  CH3 + HCO <=> CH3CHO
    kiv[87] = {10,20,26};
    nuv[87] = {-1,-1,1};
    // (81):  CH3 + HCO <=> CH3CHO
    fwd_A[87]     = 50000000000000;
    fwd_beta[87]  = 0;
    fwd_Ea[87]    = 0;
    prefactor_units[87]  = 1.0000000000000002e-06;
    activation_units[87] = 0.50321666580471969;
    phase_units[87]      = pow(10,-12.000000);
    is_PD[87] = 0;
    nTB[87] = 0;

    // (82):  HCO + M <=> CO + H + M
    kiv[19] = {20,9,1};
    nuv[19] = {-1,1,1};
    // (82):  HCO + M <=> CO + H + M
    fwd_A[19]     = 1.87e+17;
    fwd_beta[19]  = -1;
    fwd_Ea[19]    = 17000.48;
    prefactor_units[19]  = 1.0000000000000002e-06;
    activation_units[19] = 0.50321666580471969;
    phase_units[19]      = pow(10,-6.000000);
    is_PD[19] = 0;
    nTB[19] = 17;
    TB[19] = (amrex::Real *) malloc(17 * sizeof(amrex::Real));
    TBid[19] = (int *) malloc(17 * sizeof(int));
    TBid[19][0] = 2; TB[19][0] = 2; // H2
    TBid[19][1] = 7; TB[19][1] = 0; // H2O
    TBid[19][2] = 18; TB[19][2] = 0; // C
    TBid[19][3] = 9; TB[19][3] = 1.75; // CO
    TBid[19][4] = 19; TB[19][4] = 0; // CH
    TBid[19][5] = 20; TB[19][5] = 0; // HCO
    TBid[19][6] = 21; TB[19][6] = 0; // TXCH2
    TBid[19][7] = 12; TB[19][7] = 3.6000000000000001; // CO2
    TBid[19][8] = 22; TB[19][8] = 0; // SXCH2
    TBid[19][9] = 13; TB[19][9] = 2; // CH4
    TBid[19][10] = 23; TB[19][10] = 0; // C2H3
    TBid[19][11] = 24; TB[19][11] = 0; // C2H5
    TBid[19][12] = 25; TB[19][12] = 0; // HCCO
    TBid[19][13] = 26; TB[19][13] = 0; // CH3CHO
    TBid[19][14] = 27; TB[19][14] = 0; // CH2CHO
    TBid[19][15] = 28; TB[19][15] = 0; // C2H5O
    TBid[19][16] = 17; TB[19][16] = 3; // C2H6

    // (83):  HCO + H2O <=> CO + H + H2O
    kiv[88] = {20,7,9,1,7};
    nuv[88] = {-1,-1,1,1,1};
    // (83):  HCO + H2O <=> CO + H + H2O
    fwd_A[88]     = 2.24e+18;
    fwd_beta[88]  = -1;
    fwd_Ea[88]    = 17000.48;
    prefactor_units[88]  = 1.0000000000000002e-06;
    activation_units[88] = 0.50321666580471969;
    phase_units[88]      = pow(10,-12.000000);
    is_PD[88] = 0;
    nTB[88] = 0;

    // (84):  HCO + O <=> CO + OH
    kiv[89] = {20,3,9,4};
    nuv[89] = {-1,-1,1,1};
    // (84):  HCO + O <=> CO + OH
    fwd_A[89]     = 30000000000000;
    fwd_beta[89]  = 0;
    fwd_Ea[89]    = 0;
    prefactor_units[89]  = 1.0000000000000002e-06;
    activation_units[89] = 0.50321666580471969;
    phase_units[89]      = pow(10,-12.000000);
    is_PD[89] = 0;
    nTB[89] = 0;

    // (85):  HCO + OH <=> CO + H2O
    kiv[90] = {20,4,9,7};
    nuv[90] = {-1,-1,1,1};
    // (85):  HCO + OH <=> CO + H2O
    fwd_A[90]     = 30200000000000;
    fwd_beta[90]  = 0;
    fwd_Ea[90]    = 0;
    prefactor_units[90]  = 1.0000000000000002e-06;
    activation_units[90] = 0.50321666580471969;
    phase_units[90]      = pow(10,-12.000000);
    is_PD[90] = 0;
    nTB[90] = 0;

    // (86):  CH3 + HCO <=> CH4 + CO
    kiv[91] = {10,20,13,9};
    nuv[91] = {-1,-1,1,1};
    // (86):  CH3 + HCO <=> CH4 + CO
    fwd_A[91]     = 26500000000000;
    fwd_beta[91]  = 0;
    fwd_Ea[91]    = 0;
    prefactor_units[91]  = 1.0000000000000002e-06;
    activation_units[91] = 0.50321666580471969;
    phase_units[91]      = pow(10,-12.000000);
    is_PD[91] = 0;
    nTB[91] = 0;

    // (87):  HCO + O <=> CO2 + H
    kiv[92] = {20,3,12,1};
    nuv[92] = {-1,-1,1,1};
    // (87):  HCO + O <=> CO2 + H
    fwd_A[92]     = 30000000000000;
    fwd_beta[92]  = 0;
    fwd_Ea[92]    = 0;
    prefactor_units[92]  = 1.0000000000000002e-06;
    activation_units[92] = 0.50321666580471969;
    phase_units[92]      = pow(10,-12.000000);
    is_PD[92] = 0;
    nTB[92] = 0;

    // (88):  HCO + O2 <=> CO + HO2
    kiv[93] = {20,5,9,8};
    nuv[93] = {-1,-1,1,1};
    // (88):  HCO + O2 <=> CO + HO2
    fwd_A[93]     = 12000000000;
    fwd_beta[93]  = 0.81000000000000005;
    fwd_Ea[93]    = -726.58000000000004;
    prefactor_units[93]  = 1.0000000000000002e-06;
    activation_units[93] = 0.50321666580471969;
    phase_units[93]      = pow(10,-12.000000);
    is_PD[93] = 0;
    nTB[93] = 0;

    // (89):  CH2O + H <=> HCO + H2
    kiv[94] = {11,1,20,2};
    nuv[94] = {-1,-1,1,1};
    // (89):  CH2O + H <=> HCO + H2
    fwd_A[94]     = 57400000;
    fwd_beta[94]  = 1.8999999999999999;
    fwd_Ea[94]    = 2741.4000000000001;
    prefactor_units[94]  = 1.0000000000000002e-06;
    activation_units[94] = 0.50321666580471969;
    phase_units[94]      = pow(10,-12.000000);
    is_PD[94] = 0;
    nTB[94] = 0;

    // (90):  CH2O + O <=> HCO + OH
    kiv[95] = {11,3,20,4};
    nuv[95] = {-1,-1,1,1};
    // (90):  CH2O + O <=> HCO + OH
    fwd_A[95]     = 39000000000000;
    fwd_beta[95]  = 0;
    fwd_Ea[95]    = 3539.6700000000001;
    prefactor_units[95]  = 1.0000000000000002e-06;
    activation_units[95] = 0.50321666580471969;
    phase_units[95]      = pow(10,-12.000000);
    is_PD[95] = 0;
    nTB[95] = 0;

    // (91):  CH3 + CH2O <=> CH4 + HCO
    kiv[96] = {10,11,13,20};
    nuv[96] = {-1,-1,1,1};
    // (91):  CH3 + CH2O <=> CH4 + HCO
    fwd_A[96]     = 3320;
    fwd_beta[96]  = 2.8100000000000001;
    fwd_Ea[96]    = 5860.4200000000001;
    prefactor_units[96]  = 1.0000000000000002e-06;
    activation_units[96] = 0.50321666580471969;
    phase_units[96]      = pow(10,-12.000000);
    is_PD[96] = 0;
    nTB[96] = 0;

    // (92):  CH2O + OH <=> HCO + H2O
    kiv[97] = {11,4,20,7};
    nuv[97] = {-1,-1,1,1};
    // (92):  CH2O + OH <=> HCO + H2O
    fwd_A[97]     = 3430000000;
    fwd_beta[97]  = 1.1799999999999999;
    fwd_Ea[97]    = -446.94;
    prefactor_units[97]  = 1.0000000000000002e-06;
    activation_units[97] = 0.50321666580471969;
    phase_units[97]      = pow(10,-12.000000);
    is_PD[97] = 0;
    nTB[97] = 0;

    // (93):  CH2O + CH <=> CH2CO + H
    kiv[98] = {11,19,16,1};
    nuv[98] = {-1,-1,1,1};
    // (93):  CH2O + CH <=> CH2CO + H
    fwd_A[98]     = 94600000000000;
    fwd_beta[98]  = 0;
    fwd_Ea[98]    = -516.25;
    prefactor_units[98]  = 1.0000000000000002e-06;
    activation_units[98] = 0.50321666580471969;
    phase_units[98]      = pow(10,-12.000000);
    is_PD[98] = 0;
    nTB[98] = 0;

    // (94):  CH2O + O2 <=> HCO + HO2
    kiv[99] = {11,5,20,8};
    nuv[99] = {-1,-1,1,1};
    // (94):  CH2O + O2 <=> HCO + HO2
    fwd_A[99]     = 100000000000000;
    fwd_beta[99]  = 0;
    fwd_Ea[99]    = 40000;
    prefactor_units[99]  = 1.0000000000000002e-06;
    activation_units[99] = 0.50321666580471969;
    phase_units[99]      = pow(10,-12.000000);
    is_PD[99] = 0;
    nTB[99] = 0;

    // (95):  CH2O + HO2 <=> HCO + H2O2
    kiv[100] = {11,8,20,6};
    nuv[100] = {-1,-1,1,1};
    // (95):  CH2O + HO2 <=> HCO + H2O2
    fwd_A[100]     = 5600000;
    fwd_beta[100]  = 2;
    fwd_Ea[100]    = 12000.48;
    prefactor_units[100]  = 1.0000000000000002e-06;
    activation_units[100] = 0.50321666580471969;
    phase_units[100]      = pow(10,-12.000000);
    is_PD[100] = 0;
    nTB[100] = 0;

    // (96):  2.000000 H + CO2 <=> H2 + CO2
    kiv[101] = {1,12,2,12};
    nuv[101] = {-2.0,-1,1,1};
    // (96):  2.000000 H + CO2 <=> H2 + CO2
    fwd_A[101]     = 5.5e+20;
    fwd_beta[101]  = -2;
    fwd_Ea[101]    = 0;
    prefactor_units[101]  = 1.0000000000000002e-12;
    activation_units[101] = 0.50321666580471969;
    phase_units[101]      = pow(10,-18.000000);
    is_PD[101] = 0;
    nTB[101] = 0;

    // (97):  SXCH2 + CO2 <=> TXCH2 + CO2
    kiv[102] = {22,12,21,12};
    nuv[102] = {-1,-1,1,1};
    // (97):  SXCH2 + CO2 <=> TXCH2 + CO2
    fwd_A[102]     = 7000000000000;
    fwd_beta[102]  = 0;
    fwd_Ea[102]    = 0;
    prefactor_units[102]  = 1.0000000000000002e-06;
    activation_units[102] = 0.50321666580471969;
    phase_units[102]      = pow(10,-12.000000);
    is_PD[102] = 0;
    nTB[102] = 0;

    // (98):  SXCH2 + CO2 <=> CH2O + CO
    kiv[103] = {22,12,11,9};
    nuv[103] = {-1,-1,1,1};
    // (98):  SXCH2 + CO2 <=> CH2O + CO
    fwd_A[103]     = 14000000000000;
    fwd_beta[103]  = 0;
    fwd_Ea[103]    = 0;
    prefactor_units[103]  = 1.0000000000000002e-06;
    activation_units[103] = 0.50321666580471969;
    phase_units[103]      = pow(10,-12.000000);
    is_PD[103] = 0;
    nTB[103] = 0;

    // (99):  CH + CO2 <=> HCO + CO
    kiv[104] = {19,12,20,9};
    nuv[104] = {-1,-1,1,1};
    // (99):  CH + CO2 <=> HCO + CO
    fwd_A[104]     = 190000000000000;
    fwd_beta[104]  = 0;
    fwd_Ea[104]    = 15791.110000000001;
    prefactor_units[104]  = 1.0000000000000002e-06;
    activation_units[104] = 0.50321666580471969;
    phase_units[104]      = pow(10,-12.000000);
    is_PD[104] = 0;
    nTB[104] = 0;

    // (100):  C2H2 + O <=> TXCH2 + CO
    kiv[105] = {14,3,21,9};
    nuv[105] = {-1,-1,1,1};
    // (100):  C2H2 + O <=> TXCH2 + CO
    fwd_A[105]     = 12500000;
    fwd_beta[105]  = 2;
    fwd_Ea[105]    = 1900.0999999999999;
    prefactor_units[105]  = 1.0000000000000002e-06;
    activation_units[105] = 0.50321666580471969;
    phase_units[105]      = pow(10,-12.000000);
    is_PD[105] = 0;
    nTB[105] = 0;

    // (101):  C2H2 + OH <=> CH3 + CO
    kiv[106] = {14,4,10,9};
    nuv[106] = {-1,-1,1,1};
    // (101):  C2H2 + OH <=> CH3 + CO
    fwd_A[106]     = 1280000000;
    fwd_beta[106]  = 0.72999999999999998;
    fwd_Ea[106]    = 2578.8699999999999;
    prefactor_units[106]  = 1.0000000000000002e-06;
    activation_units[106] = 0.50321666580471969;
    phase_units[106]      = pow(10,-12.000000);
    is_PD[106] = 0;
    nTB[106] = 0;

    // (102):  C2H2 + H (+M) <=> C2H3 (+M)
    kiv[10] = {14,1,23};
    nuv[10] = {-1,-1,1};
    // (102):  C2H2 + H (+M) <=> C2H3 (+M)
    fwd_A[10]     = 17100000000;
    fwd_beta[10]  = 1.27;
    fwd_Ea[10]    = 2707.9299999999998;
    low_A[10]     = 6.3399999999999996e+31;
    low_beta[10]  = -4.6600000000000001;
    low_Ea[10]    = 3781.0700000000002;
    troe_a[10]    = 0.2122;
    troe_Tsss[10] = 1;
    troe_Ts[10]   = -10210;
    troe_len[10]  = 3;
    prefactor_units[10]  = 1.0000000000000002e-06;
    activation_units[10] = 0.50321666580471969;
    phase_units[10]      = pow(10,-12.000000);
    is_PD[10] = 1;
    nTB[10] = 17;
    TB[10] = (amrex::Real *) malloc(17 * sizeof(amrex::Real));
    TBid[10] = (int *) malloc(17 * sizeof(int));
    TBid[10][0] = 2; TB[10][0] = 2; // H2
    TBid[10][1] = 7; TB[10][1] = 12; // H2O
    TBid[10][2] = 18; TB[10][2] = 0; // C
    TBid[10][3] = 9; TB[10][3] = 1.75; // CO
    TBid[10][4] = 19; TB[10][4] = 0; // CH
    TBid[10][5] = 20; TB[10][5] = 0; // HCO
    TBid[10][6] = 21; TB[10][6] = 0; // TXCH2
    TBid[10][7] = 12; TB[10][7] = 3.6000000000000001; // CO2
    TBid[10][8] = 22; TB[10][8] = 0; // SXCH2
    TBid[10][9] = 13; TB[10][9] = 2; // CH4
    TBid[10][10] = 23; TB[10][10] = 0; // C2H3
    TBid[10][11] = 24; TB[10][11] = 0; // C2H5
    TBid[10][12] = 25; TB[10][12] = 0; // HCCO
    TBid[10][13] = 26; TB[10][13] = 0; // CH3CHO
    TBid[10][14] = 27; TB[10][14] = 0; // CH2CHO
    TBid[10][15] = 28; TB[10][15] = 0; // C2H5O
    TBid[10][16] = 17; TB[10][16] = 3; // C2H6

    // (103):  C2H2 + OH <=> CH2CO + H
    kiv[107] = {14,4,16,1};
    nuv[107] = {-1,-1,1,1};
    // (103):  C2H2 + OH <=> CH2CO + H
    fwd_A[107]     = 7530000;
    fwd_beta[107]  = 1.55;
    fwd_Ea[107]    = 2105.6399999999999;
    prefactor_units[107]  = 1.0000000000000002e-06;
    activation_units[107] = 0.50321666580471969;
    phase_units[107]      = pow(10,-12.000000);
    is_PD[107] = 0;
    nTB[107] = 0;

    // (104):  C2H2 + O <=> HCCO + H
    kiv[108] = {14,3,25,1};
    nuv[108] = {-1,-1,1,1};
    // (104):  C2H2 + O <=> HCCO + H
    fwd_A[108]     = 8100000;
    fwd_beta[108]  = 2;
    fwd_Ea[108]    = 1900.0999999999999;
    prefactor_units[108]  = 1.0000000000000002e-06;
    activation_units[108] = 0.50321666580471969;
    phase_units[108]      = pow(10,-12.000000);
    is_PD[108] = 0;
    nTB[108] = 0;

    // (105):  C2H3 + OH <=> C2H2 + H2O
    kiv[109] = {23,4,14,7};
    nuv[109] = {-1,-1,1,1};
    // (105):  C2H3 + OH <=> C2H2 + H2O
    fwd_A[109]     = 5000000000000;
    fwd_beta[109]  = 0;
    fwd_Ea[109]    = 0;
    prefactor_units[109]  = 1.0000000000000002e-06;
    activation_units[109] = 0.50321666580471969;
    phase_units[109]      = pow(10,-12.000000);
    is_PD[109] = 0;
    nTB[109] = 0;

    // (106):  C2H3 + O2 <=> CH2CHO + O
    kiv[110] = {23,5,27,3};
    nuv[110] = {-1,-1,1,1};
    // (106):  C2H3 + O2 <=> CH2CHO + O
    fwd_A[110]     = 303000000000;
    fwd_beta[110]  = 0.28999999999999998;
    fwd_Ea[110]    = 11.949999999999999;
    prefactor_units[110]  = 1.0000000000000002e-06;
    activation_units[110] = 0.50321666580471969;
    phase_units[110]      = pow(10,-12.000000);
    is_PD[110] = 0;
    nTB[110] = 0;

    // (107):  C2H3 + O <=> CH2CHO
    kiv[111] = {23,3,27};
    nuv[111] = {-1,-1,1};
    // (107):  C2H3 + O <=> CH2CHO
    fwd_A[111]     = 10300000000000;
    fwd_beta[111]  = 0.20999999999999999;
    fwd_Ea[111]    = -427.81999999999999;
    prefactor_units[111]  = 1.0000000000000002e-06;
    activation_units[111] = 0.50321666580471969;
    phase_units[111]      = pow(10,-12.000000);
    is_PD[111] = 0;
    nTB[111] = 0;

    // (108):  C2H3 + H <=> C2H2 + H2
    kiv[112] = {23,1,14,2};
    nuv[112] = {-1,-1,1,1};
    // (108):  C2H3 + H <=> C2H2 + H2
    fwd_A[112]     = 30000000000000;
    fwd_beta[112]  = 0;
    fwd_Ea[112]    = 0;
    prefactor_units[112]  = 1.0000000000000002e-06;
    activation_units[112] = 0.50321666580471969;
    phase_units[112]      = pow(10,-12.000000);
    is_PD[112] = 0;
    nTB[112] = 0;

    // (109):  C2H3 + CH3 <=> C2H2 + CH4
    kiv[113] = {23,10,14,13};
    nuv[113] = {-1,-1,1,1};
    // (109):  C2H3 + CH3 <=> C2H2 + CH4
    fwd_A[113]     = 9030000000000;
    fwd_beta[113]  = 0;
    fwd_Ea[113]    = -764.82000000000005;
    prefactor_units[113]  = 1.0000000000000002e-06;
    activation_units[113] = 0.50321666580471969;
    phase_units[113]      = pow(10,-12.000000);
    is_PD[113] = 0;
    nTB[113] = 0;

    // (110):  C2H3 + O2 <=> HCO + CH2O
    kiv[114] = {23,5,20,11};
    nuv[114] = {-1,-1,1,1};
    // (110):  C2H3 + O2 <=> HCO + CH2O
    fwd_A[114]     = 45800000000000000;
    fwd_beta[114]  = -1.3899999999999999;
    fwd_Ea[114]    = 1015.77;
    prefactor_units[114]  = 1.0000000000000002e-06;
    activation_units[114] = 0.50321666580471969;
    phase_units[114]      = pow(10,-12.000000);
    is_PD[114] = 0;
    nTB[114] = 0;

    // (111):  C2H3 + H (+M) <=> C2H4 (+M)
    kiv[11] = {23,1,15};
    nuv[11] = {-1,-1,1};
    // (111):  C2H3 + H (+M) <=> C2H4 (+M)
    fwd_A[11]     = 6080000000000;
    fwd_beta[11]  = 0.27000000000000002;
    fwd_Ea[11]    = 279.63999999999999;
    low_A[11]     = 1.3999999999999999e+30;
    low_beta[11]  = -3.8599999999999999;
    low_Ea[11]    = 3319.79;
    troe_a[11]    = 0.78200000000000003;
    troe_Tsss[11] = 207.5;
    troe_Ts[11]   = 2663;
    troe_Tss[11]  = 6095;
    troe_len[11]  = 4;
    prefactor_units[11]  = 1.0000000000000002e-06;
    activation_units[11] = 0.50321666580471969;
    phase_units[11]      = pow(10,-12.000000);
    is_PD[11] = 1;
    nTB[11] = 17;
    TB[11] = (amrex::Real *) malloc(17 * sizeof(amrex::Real));
    TBid[11] = (int *) malloc(17 * sizeof(int));
    TBid[11][0] = 2; TB[11][0] = 2; // H2
    TBid[11][1] = 7; TB[11][1] = 12; // H2O
    TBid[11][2] = 18; TB[11][2] = 0; // C
    TBid[11][3] = 9; TB[11][3] = 1.75; // CO
    TBid[11][4] = 19; TB[11][4] = 0; // CH
    TBid[11][5] = 20; TB[11][5] = 0; // HCO
    TBid[11][6] = 21; TB[11][6] = 0; // TXCH2
    TBid[11][7] = 12; TB[11][7] = 3.6000000000000001; // CO2
    TBid[11][8] = 22; TB[11][8] = 0; // SXCH2
    TBid[11][9] = 13; TB[11][9] = 2; // CH4
    TBid[11][10] = 23; TB[11][10] = 0; // C2H3
    TBid[11][11] = 24; TB[11][11] = 0; // C2H5
    TBid[11][12] = 25; TB[11][12] = 0; // HCCO
    TBid[11][13] = 26; TB[11][13] = 0; // CH3CHO
    TBid[11][14] = 27; TB[11][14] = 0; // CH2CHO
    TBid[11][15] = 28; TB[11][15] = 0; // C2H5O
    TBid[11][16] = 17; TB[11][16] = 3; // C2H6

    // (112):  C2H3 + H2O2 <=> C2H4 + HO2
    kiv[115] = {23,6,15,8};
    nuv[115] = {-1,-1,1,1};
    // (112):  C2H3 + H2O2 <=> C2H4 + HO2
    fwd_A[115]     = 12100000000;
    fwd_beta[115]  = 0;
    fwd_Ea[115]    = -595.12;
    prefactor_units[115]  = 1.0000000000000002e-06;
    activation_units[115] = 0.50321666580471969;
    phase_units[115]      = pow(10,-12.000000);
    is_PD[115] = 0;
    nTB[115] = 0;

    // (113):  C2H3 + O2 <=> C2H2 + HO2
    kiv[116] = {23,5,14,8};
    nuv[116] = {-1,-1,1,1};
    // (113):  C2H3 + O2 <=> C2H2 + HO2
    fwd_A[116]     = 1340000;
    fwd_beta[116]  = 1.6100000000000001;
    fwd_Ea[116]    = -384.80000000000001;
    prefactor_units[116]  = 1.0000000000000002e-06;
    activation_units[116] = 0.50321666580471969;
    phase_units[116]      = pow(10,-12.000000);
    is_PD[116] = 0;
    nTB[116] = 0;

    // (114):  C2H4 + CH3 <=> C2H3 + CH4
    kiv[117] = {15,10,23,13};
    nuv[117] = {-1,-1,1,1};
    // (114):  C2H4 + CH3 <=> C2H3 + CH4
    fwd_A[117]     = 227000;
    fwd_beta[117]  = 2;
    fwd_Ea[117]    = 9199.3299999999999;
    prefactor_units[117]  = 1.0000000000000002e-06;
    activation_units[117] = 0.50321666580471969;
    phase_units[117]      = pow(10,-12.000000);
    is_PD[117] = 0;
    nTB[117] = 0;

    // (115):  C2H4 + H (+M) <=> C2H5 (+M)
    kiv[12] = {15,1,24};
    nuv[12] = {-1,-1,1};
    // (115):  C2H4 + H (+M) <=> C2H5 (+M)
    fwd_A[12]     = 1370000000;
    fwd_beta[12]  = 1.46;
    fwd_Ea[12]    = 1355.1600000000001;
    low_A[12]     = 2.0300000000000001e+39;
    low_beta[12]  = -6.6399999999999997;
    low_Ea[12]    = 5769.6000000000004;
    troe_a[12]    = -0.56899999999999995;
    troe_Tsss[12] = 299;
    troe_Ts[12]   = -9147;
    troe_Tss[12]  = 152.40000000000001;
    troe_len[12]  = 4;
    prefactor_units[12]  = 1.0000000000000002e-06;
    activation_units[12] = 0.50321666580471969;
    phase_units[12]      = pow(10,-12.000000);
    is_PD[12] = 1;
    nTB[12] = 17;
    TB[12] = (amrex::Real *) malloc(17 * sizeof(amrex::Real));
    TBid[12] = (int *) malloc(17 * sizeof(int));
    TBid[12][0] = 2; TB[12][0] = 2; // H2
    TBid[12][1] = 7; TB[12][1] = 12; // H2O
    TBid[12][2] = 18; TB[12][2] = 0; // C
    TBid[12][3] = 9; TB[12][3] = 1.75; // CO
    TBid[12][4] = 19; TB[12][4] = 0; // CH
    TBid[12][5] = 20; TB[12][5] = 0; // HCO
    TBid[12][6] = 21; TB[12][6] = 0; // TXCH2
    TBid[12][7] = 12; TB[12][7] = 3.6000000000000001; // CO2
    TBid[12][8] = 22; TB[12][8] = 0; // SXCH2
    TBid[12][9] = 13; TB[12][9] = 2; // CH4
    TBid[12][10] = 23; TB[12][10] = 0; // C2H3
    TBid[12][11] = 24; TB[12][11] = 0; // C2H5
    TBid[12][12] = 25; TB[12][12] = 0; // HCCO
    TBid[12][13] = 26; TB[12][13] = 0; // CH3CHO
    TBid[12][14] = 27; TB[12][14] = 0; // CH2CHO
    TBid[12][15] = 28; TB[12][15] = 0; // C2H5O
    TBid[12][16] = 17; TB[12][16] = 3; // C2H6

    // (116):  C2H4 + O2 => CH3 + CO2 + H
    kiv[118] = {15,5,10,12,1};
    nuv[118] = {-1,-1,1,1,1};
    // (116):  C2H4 + O2 => CH3 + CO2 + H
    fwd_A[118]     = 4900000000000;
    fwd_beta[118]  = 0.41999999999999998;
    fwd_Ea[118]    = 75800.669999999998;
    prefactor_units[118]  = 1.0000000000000002e-06;
    activation_units[118] = 0.50321666580471969;
    phase_units[118]      = pow(10,-12.000000);
    is_PD[118] = 0;
    nTB[118] = 0;

    // (117):  C2H4 + OH <=> C2H3 + H2O
    kiv[119] = {15,4,23,7};
    nuv[119] = {-1,-1,1,1};
    // (117):  C2H4 + OH <=> C2H3 + H2O
    fwd_A[119]     = 0.13100000000000001;
    fwd_beta[119]  = 4.2000000000000002;
    fwd_Ea[119]    = -860.41999999999996;
    prefactor_units[119]  = 1.0000000000000002e-06;
    activation_units[119] = 0.50321666580471969;
    phase_units[119]      = pow(10,-12.000000);
    is_PD[119] = 0;
    nTB[119] = 0;

    // (118):  C2H4 + OH <=> C2H5O
    kiv[120] = {15,4,28};
    nuv[120] = {-1,-1,1};
    // (118):  C2H4 + OH <=> C2H5O
    fwd_A[120]     = 3.7500000000000003e+36;
    fwd_beta[120]  = -7.7999999999999998;
    fwd_Ea[120]    = 7060.2299999999996;
    prefactor_units[120]  = 1.0000000000000002e-06;
    activation_units[120] = 0.50321666580471969;
    phase_units[120]      = pow(10,-12.000000);
    is_PD[120] = 0;
    nTB[120] = 0;

    // (119):  C2H4 + O <=> CH2CHO + H
    kiv[121] = {15,3,27,1};
    nuv[121] = {-1,-1,1,1};
    // (119):  C2H4 + O <=> CH2CHO + H
    fwd_A[121]     = 7660000000;
    fwd_beta[121]  = 0.88;
    fwd_Ea[121]    = 1140.0599999999999;
    prefactor_units[121]  = 1.0000000000000002e-06;
    activation_units[121] = 0.50321666580471969;
    phase_units[121]      = pow(10,-12.000000);
    is_PD[121] = 0;
    nTB[121] = 0;

    // (120):  C2H4 + O <=> CH3 + HCO
    kiv[122] = {15,3,10,20};
    nuv[122] = {-1,-1,1,1};
    // (120):  C2H4 + O <=> CH3 + HCO
    fwd_A[122]     = 389000000;
    fwd_beta[122]  = 1.3600000000000001;
    fwd_Ea[122]    = 886.71000000000004;
    prefactor_units[122]  = 1.0000000000000002e-06;
    activation_units[122] = 0.50321666580471969;
    phase_units[122]      = pow(10,-12.000000);
    is_PD[122] = 0;
    nTB[122] = 0;

    // (121):  C2H4 + O2 <=> C2H3 + HO2
    kiv[123] = {15,5,23,8};
    nuv[123] = {-1,-1,1,1};
    // (121):  C2H4 + O2 <=> C2H3 + HO2
    fwd_A[123]     = 42200000000000;
    fwd_beta[123]  = 0;
    fwd_Ea[123]    = 62100.860000000001;
    prefactor_units[123]  = 1.0000000000000002e-06;
    activation_units[123] = 0.50321666580471969;
    phase_units[123]      = pow(10,-12.000000);
    is_PD[123] = 0;
    nTB[123] = 0;

    // (122):  C2H4 + H <=> C2H3 + H2
    kiv[124] = {15,1,23,2};
    nuv[124] = {-1,-1,1,1};
    // (122):  C2H4 + H <=> C2H3 + H2
    fwd_A[124]     = 127000;
    fwd_beta[124]  = 2.75;
    fwd_Ea[124]    = 11649.139999999999;
    prefactor_units[124]  = 1.0000000000000002e-06;
    activation_units[124] = 0.50321666580471969;
    phase_units[124]      = pow(10,-12.000000);
    is_PD[124] = 0;
    nTB[124] = 0;

    // (123):  C2H4 + O <=> TXCH2 + CH2O
    kiv[125] = {15,3,21,11};
    nuv[125] = {-1,-1,1,1};
    // (123):  C2H4 + O <=> TXCH2 + CH2O
    fwd_A[125]     = 71500;
    fwd_beta[125]  = 2.4700000000000002;
    fwd_Ea[125]    = 929.73000000000002;
    prefactor_units[125]  = 1.0000000000000002e-06;
    activation_units[125] = 0.50321666580471969;
    phase_units[125]      = pow(10,-12.000000);
    is_PD[125] = 0;
    nTB[125] = 0;

    // (124):  C2H5 + HO2 <=> C2H4 + H2O2
    kiv[126] = {24,8,15,6};
    nuv[126] = {-1,-1,1,1};
    // (124):  C2H5 + HO2 <=> C2H4 + H2O2
    fwd_A[126]     = 300000000000;
    fwd_beta[126]  = 0;
    fwd_Ea[126]    = 0;
    prefactor_units[126]  = 1.0000000000000002e-06;
    activation_units[126] = 0.50321666580471969;
    phase_units[126]      = pow(10,-12.000000);
    is_PD[126] = 0;
    nTB[126] = 0;

    // (125):  C2H5 + H (+M) <=> C2H6 (+M)
    kiv[13] = {24,1,17};
    nuv[13] = {-1,-1,1};
    // (125):  C2H5 + H (+M) <=> C2H6 (+M)
    fwd_A[13]     = 5.21e+17;
    fwd_beta[13]  = -0.98999999999999999;
    fwd_Ea[13]    = 1579.8299999999999;
    low_A[13]     = 1.9900000000000001e+41;
    low_beta[13]  = -7.0800000000000001;
    low_Ea[13]    = 6684.9899999999998;
    troe_a[13]    = 0.84219999999999995;
    troe_Tsss[13] = 125;
    troe_Ts[13]   = 2219;
    troe_Tss[13]  = 6882;
    troe_len[13]  = 4;
    prefactor_units[13]  = 1.0000000000000002e-06;
    activation_units[13] = 0.50321666580471969;
    phase_units[13]      = pow(10,-12.000000);
    is_PD[13] = 1;
    nTB[13] = 17;
    TB[13] = (amrex::Real *) malloc(17 * sizeof(amrex::Real));
    TBid[13] = (int *) malloc(17 * sizeof(int));
    TBid[13][0] = 2; TB[13][0] = 2; // H2
    TBid[13][1] = 7; TB[13][1] = 12; // H2O
    TBid[13][2] = 18; TB[13][2] = 0; // C
    TBid[13][3] = 9; TB[13][3] = 1.75; // CO
    TBid[13][4] = 19; TB[13][4] = 0; // CH
    TBid[13][5] = 20; TB[13][5] = 0; // HCO
    TBid[13][6] = 21; TB[13][6] = 0; // TXCH2
    TBid[13][7] = 12; TB[13][7] = 3.6000000000000001; // CO2
    TBid[13][8] = 22; TB[13][8] = 0; // SXCH2
    TBid[13][9] = 13; TB[13][9] = 2; // CH4
    TBid[13][10] = 23; TB[13][10] = 0; // C2H3
    TBid[13][11] = 24; TB[13][11] = 0; // C2H5
    TBid[13][12] = 25; TB[13][12] = 0; // HCCO
    TBid[13][13] = 26; TB[13][13] = 0; // CH3CHO
    TBid[13][14] = 27; TB[13][14] = 0; // CH2CHO
    TBid[13][15] = 28; TB[13][15] = 0; // C2H5O
    TBid[13][16] = 17; TB[13][16] = 3; // C2H6

    // (126):  C2H5 + HO2 <=> C2H5O + OH
    kiv[127] = {24,8,28,4};
    nuv[127] = {-1,-1,1,1};
    // (126):  C2H5 + HO2 <=> C2H5O + OH
    fwd_A[127]     = 31000000000000;
    fwd_beta[127]  = 0;
    fwd_Ea[127]    = 0;
    prefactor_units[127]  = 1.0000000000000002e-06;
    activation_units[127] = 0.50321666580471969;
    phase_units[127]      = pow(10,-12.000000);
    is_PD[127] = 0;
    nTB[127] = 0;

    // (127):  C2H5 + O <=> C2H5O
    kiv[128] = {24,3,28};
    nuv[128] = {-1,-1,1};
    // (127):  C2H5 + O <=> C2H5O
    fwd_A[128]     = 31700000000000;
    fwd_beta[128]  = 0.029999999999999999;
    fwd_Ea[128]    = -394.36000000000001;
    prefactor_units[128]  = 1.0000000000000002e-06;
    activation_units[128] = 0.50321666580471969;
    phase_units[128]      = pow(10,-12.000000);
    is_PD[128] = 0;
    nTB[128] = 0;

    // (128):  C2H5 + H <=> C2H4 + H2
    kiv[129] = {24,1,15,2};
    nuv[129] = {-1,-1,1,1};
    // (128):  C2H5 + H <=> C2H4 + H2
    fwd_A[129]     = 2000000000000;
    fwd_beta[129]  = 0;
    fwd_Ea[129]    = 0;
    prefactor_units[129]  = 1.0000000000000002e-06;
    activation_units[129] = 0.50321666580471969;
    phase_units[129]      = pow(10,-12.000000);
    is_PD[129] = 0;
    nTB[129] = 0;

    // (129):  C2H5 + O2 <=> C2H4 + HO2
    kiv[130] = {24,5,15,8};
    nuv[130] = {-1,-1,1,1};
    // (129):  C2H5 + O2 <=> C2H4 + HO2
    fwd_A[130]     = 19200000;
    fwd_beta[130]  = 1.02;
    fwd_Ea[130]    = -2033.9400000000001;
    prefactor_units[130]  = 1.0000000000000002e-06;
    activation_units[130] = 0.50321666580471969;
    phase_units[130]      = pow(10,-12.000000);
    is_PD[130] = 0;
    nTB[130] = 0;

    // (130):  C2H5 + HO2 <=> C2H6 + O2
    kiv[131] = {24,8,17,5};
    nuv[131] = {-1,-1,1,1};
    // (130):  C2H5 + HO2 <=> C2H6 + O2
    fwd_A[131]     = 300000000000;
    fwd_beta[131]  = 0;
    fwd_Ea[131]    = 0;
    prefactor_units[131]  = 1.0000000000000002e-06;
    activation_units[131] = 0.50321666580471969;
    phase_units[131]      = pow(10,-12.000000);
    is_PD[131] = 0;
    nTB[131] = 0;

    // (131):  C2H5 + CH3 <=> C2H4 + CH4
    kiv[132] = {24,10,15,13};
    nuv[132] = {-1,-1,1,1};
    // (131):  C2H5 + CH3 <=> C2H4 + CH4
    fwd_A[132]     = 11800;
    fwd_beta[132]  = 2.4500000000000002;
    fwd_Ea[132]    = 2920.6500000000001;
    prefactor_units[132]  = 1.0000000000000002e-06;
    activation_units[132] = 0.50321666580471969;
    phase_units[132]      = pow(10,-12.000000);
    is_PD[132] = 0;
    nTB[132] = 0;

    // (132):  C2H6 + SXCH2 <=> C2H5 + CH3
    kiv[133] = {17,22,24,10};
    nuv[133] = {-1,-1,1,1};
    // (132):  C2H6 + SXCH2 <=> C2H5 + CH3
    fwd_A[133]     = 40000000000000;
    fwd_beta[133]  = 0;
    fwd_Ea[133]    = -549.71000000000004;
    prefactor_units[133]  = 1.0000000000000002e-06;
    activation_units[133] = 0.50321666580471969;
    phase_units[133]      = pow(10,-12.000000);
    is_PD[133] = 0;
    nTB[133] = 0;

    // (133):  C2H6 + CH3 <=> C2H5 + CH4
    kiv[134] = {17,10,24,13};
    nuv[134] = {-1,-1,1,1};
    // (133):  C2H6 + CH3 <=> C2H5 + CH4
    fwd_A[134]     = 843000000000000;
    fwd_beta[134]  = 0;
    fwd_Ea[134]    = 22256.209999999999;
    prefactor_units[134]  = 1.0000000000000002e-06;
    activation_units[134] = 0.50321666580471969;
    phase_units[134]      = pow(10,-12.000000);
    is_PD[134] = 0;
    nTB[134] = 0;

    // (134):  C2H6 + O <=> C2H5 + OH
    kiv[135] = {17,3,24,4};
    nuv[135] = {-1,-1,1,1};
    // (134):  C2H6 + O <=> C2H5 + OH
    fwd_A[135]     = 31.699999999999999;
    fwd_beta[135]  = 3.7999999999999998;
    fwd_Ea[135]    = 3130.98;
    prefactor_units[135]  = 1.0000000000000002e-06;
    activation_units[135] = 0.50321666580471969;
    phase_units[135]      = pow(10,-12.000000);
    is_PD[135] = 0;
    nTB[135] = 0;

    // (135):  C2H6 (+M) <=> 2.000000 CH3 (+M)
    kiv[14] = {17,10};
    nuv[14] = {-1,2.0};
    // (135):  C2H6 (+M) <=> 2.000000 CH3 (+M)
    fwd_A[14]     = 1.8799999999999999e+50;
    fwd_beta[14]  = -9.7200000000000006;
    fwd_Ea[14]    = 107342.25999999999;
    low_A[14]     = 3.7199999999999999e+65;
    low_beta[14]  = -13.140000000000001;
    low_Ea[14]    = 101579.83;
    troe_a[14]    = 0.39000000000000001;
    troe_Tsss[14] = 100;
    troe_Ts[14]   = 1900;
    troe_Tss[14]  = 6000;
    troe_len[14]  = 4;
    prefactor_units[14]  = 1;
    activation_units[14] = 0.50321666580471969;
    phase_units[14]      = pow(10,-6.000000);
    is_PD[14] = 1;
    nTB[14] = 17;
    TB[14] = (amrex::Real *) malloc(17 * sizeof(amrex::Real));
    TBid[14] = (int *) malloc(17 * sizeof(int));
    TBid[14][0] = 2; TB[14][0] = 2; // H2
    TBid[14][1] = 7; TB[14][1] = 12; // H2O
    TBid[14][2] = 18; TB[14][2] = 0; // C
    TBid[14][3] = 9; TB[14][3] = 1.75; // CO
    TBid[14][4] = 19; TB[14][4] = 0; // CH
    TBid[14][5] = 20; TB[14][5] = 0; // HCO
    TBid[14][6] = 21; TB[14][6] = 0; // TXCH2
    TBid[14][7] = 12; TB[14][7] = 3.6000000000000001; // CO2
    TBid[14][8] = 22; TB[14][8] = 0; // SXCH2
    TBid[14][9] = 13; TB[14][9] = 2; // CH4
    TBid[14][10] = 23; TB[14][10] = 0; // C2H3
    TBid[14][11] = 24; TB[14][11] = 0; // C2H5
    TBid[14][12] = 25; TB[14][12] = 0; // HCCO
    TBid[14][13] = 26; TB[14][13] = 0; // CH3CHO
    TBid[14][14] = 27; TB[14][14] = 0; // CH2CHO
    TBid[14][15] = 28; TB[14][15] = 0; // C2H5O
    TBid[14][16] = 17; TB[14][16] = 3; // C2H6

    // (136):  C2H6 + HO2 <=> C2H5 + H2O2
    kiv[136] = {17,8,24,6};
    nuv[136] = {-1,-1,1,1};
    // (136):  C2H6 + HO2 <=> C2H5 + H2O2
    fwd_A[136]     = 261;
    fwd_beta[136]  = 3.3700000000000001;
    fwd_Ea[136]    = 15913;
    prefactor_units[136]  = 1.0000000000000002e-06;
    activation_units[136] = 0.50321666580471969;
    phase_units[136]      = pow(10,-12.000000);
    is_PD[136] = 0;
    nTB[136] = 0;

    // (137):  C2H6 + H <=> C2H5 + H2
    kiv[137] = {17,1,24,2};
    nuv[137] = {-1,-1,1,1};
    // (137):  C2H6 + H <=> C2H5 + H2
    fwd_A[137]     = 170000;
    fwd_beta[137]  = 2.7000000000000002;
    fwd_Ea[137]    = 5740.9200000000001;
    prefactor_units[137]  = 1.0000000000000002e-06;
    activation_units[137] = 0.50321666580471969;
    phase_units[137]      = pow(10,-12.000000);
    is_PD[137] = 0;
    nTB[137] = 0;

    // (138):  C2H6 + OH <=> C2H5 + H2O
    kiv[138] = {17,4,24,7};
    nuv[138] = {-1,-1,1,1};
    // (138):  C2H6 + OH <=> C2H5 + H2O
    fwd_A[138]     = 1610000;
    fwd_beta[138]  = 2.2200000000000002;
    fwd_Ea[138]    = 740.91999999999996;
    prefactor_units[138]  = 1.0000000000000002e-06;
    activation_units[138] = 0.50321666580471969;
    phase_units[138]      = pow(10,-12.000000);
    is_PD[138] = 0;
    nTB[138] = 0;

    // (139):  HCCO + O2 <=> OH + 2.000000 CO
    kiv[139] = {25,5,4,9};
    nuv[139] = {-1,-1,1,2.0};
    // (139):  HCCO + O2 <=> OH + 2.000000 CO
    fwd_A[139]     = 42000000000;
    fwd_beta[139]  = 0;
    fwd_Ea[139]    = 853.25;
    prefactor_units[139]  = 1.0000000000000002e-06;
    activation_units[139] = 0.50321666580471969;
    phase_units[139]      = pow(10,-12.000000);
    is_PD[139] = 0;
    nTB[139] = 0;

    // (140):  HCCO + O <=> H + 2.000000 CO
    kiv[140] = {25,3,1,9};
    nuv[140] = {-1,-1,1,2.0};
    // (140):  HCCO + O <=> H + 2.000000 CO
    fwd_A[140]     = 100000000000000;
    fwd_beta[140]  = 0;
    fwd_Ea[140]    = 0;
    prefactor_units[140]  = 1.0000000000000002e-06;
    activation_units[140] = 0.50321666580471969;
    phase_units[140]      = pow(10,-12.000000);
    is_PD[140] = 0;
    nTB[140] = 0;

    // (141):  HCCO + CH3 <=> C2H4 + CO
    kiv[141] = {25,10,15,9};
    nuv[141] = {-1,-1,1,1};
    // (141):  HCCO + CH3 <=> C2H4 + CO
    fwd_A[141]     = 50000000000000;
    fwd_beta[141]  = 0;
    fwd_Ea[141]    = 0;
    prefactor_units[141]  = 1.0000000000000002e-06;
    activation_units[141] = 0.50321666580471969;
    phase_units[141]      = pow(10,-12.000000);
    is_PD[141] = 0;
    nTB[141] = 0;

    // (142):  HCCO + H <=> SXCH2 + CO
    kiv[142] = {25,1,22,9};
    nuv[142] = {-1,-1,1,1};
    // (142):  HCCO + H <=> SXCH2 + CO
    fwd_A[142]     = 100000000000000;
    fwd_beta[142]  = 0;
    fwd_Ea[142]    = 0;
    prefactor_units[142]  = 1.0000000000000002e-06;
    activation_units[142] = 0.50321666580471969;
    phase_units[142]      = pow(10,-12.000000);
    is_PD[142] = 0;
    nTB[142] = 0;

    // (143):  CH2CO + H <=> CH3 + CO
    kiv[143] = {16,1,10,9};
    nuv[143] = {-1,-1,1,1};
    // (143):  CH2CO + H <=> CH3 + CO
    fwd_A[143]     = 1500000000;
    fwd_beta[143]  = 1.3799999999999999;
    fwd_Ea[143]    = 614.24000000000001;
    prefactor_units[143]  = 1.0000000000000002e-06;
    activation_units[143] = 0.50321666580471969;
    phase_units[143]      = pow(10,-12.000000);
    is_PD[143] = 0;
    nTB[143] = 0;

    // (144):  CH2CO + TXCH2 <=> C2H4 + CO
    kiv[144] = {16,21,15,9};
    nuv[144] = {-1,-1,1,1};
    // (144):  CH2CO + TXCH2 <=> C2H4 + CO
    fwd_A[144]     = 1000000000000;
    fwd_beta[144]  = 0;
    fwd_Ea[144]    = 0;
    prefactor_units[144]  = 1.0000000000000002e-06;
    activation_units[144] = 0.50321666580471969;
    phase_units[144]      = pow(10,-12.000000);
    is_PD[144] = 0;
    nTB[144] = 0;

    // (145):  CH2CO + O <=> HCCO + OH
    kiv[145] = {16,3,25,4};
    nuv[145] = {-1,-1,1,1};
    // (145):  CH2CO + O <=> HCCO + OH
    fwd_A[145]     = 10000000000000;
    fwd_beta[145]  = 0;
    fwd_Ea[145]    = 7999.5200000000004;
    prefactor_units[145]  = 1.0000000000000002e-06;
    activation_units[145] = 0.50321666580471969;
    phase_units[145]      = pow(10,-12.000000);
    is_PD[145] = 0;
    nTB[145] = 0;

    // (146):  CH2CO + CH3 <=> HCCO + CH4
    kiv[146] = {16,10,25,13};
    nuv[146] = {-1,-1,1,1};
    // (146):  CH2CO + CH3 <=> HCCO + CH4
    fwd_A[146]     = 7500000000000;
    fwd_beta[146]  = 0;
    fwd_Ea[146]    = 12999.52;
    prefactor_units[146]  = 1.0000000000000002e-06;
    activation_units[146] = 0.50321666580471969;
    phase_units[146]      = pow(10,-12.000000);
    is_PD[146] = 0;
    nTB[146] = 0;

    // (147):  CH2CO + O <=> TXCH2 + CO2
    kiv[147] = {16,3,21,12};
    nuv[147] = {-1,-1,1,1};
    // (147):  CH2CO + O <=> TXCH2 + CO2
    fwd_A[147]     = 1750000000000;
    fwd_beta[147]  = 0;
    fwd_Ea[147]    = 1350.3800000000001;
    prefactor_units[147]  = 1.0000000000000002e-06;
    activation_units[147] = 0.50321666580471969;
    phase_units[147]      = pow(10,-12.000000);
    is_PD[147] = 0;
    nTB[147] = 0;

    // (148):  CH2CO + CH3 <=> C2H5 + CO
    kiv[148] = {16,10,24,9};
    nuv[148] = {-1,-1,1,1};
    // (148):  CH2CO + CH3 <=> C2H5 + CO
    fwd_A[148]     = 90000000000;
    fwd_beta[148]  = 0;
    fwd_Ea[148]    = 0;
    prefactor_units[148]  = 1.0000000000000002e-06;
    activation_units[148] = 0.50321666580471969;
    phase_units[148]      = pow(10,-12.000000);
    is_PD[148] = 0;
    nTB[148] = 0;

    // (149):  CH2CO + OH <=> HCCO + H2O
    kiv[149] = {16,4,25,7};
    nuv[149] = {-1,-1,1,1};
    // (149):  CH2CO + OH <=> HCCO + H2O
    fwd_A[149]     = 7500000000000;
    fwd_beta[149]  = 0;
    fwd_Ea[149]    = 2000.48;
    prefactor_units[149]  = 1.0000000000000002e-06;
    activation_units[149] = 0.50321666580471969;
    phase_units[149]      = pow(10,-12.000000);
    is_PD[149] = 0;
    nTB[149] = 0;

    // (150):  CH2CO + H <=> HCCO + H2
    kiv[150] = {16,1,25,2};
    nuv[150] = {-1,-1,1,1};
    // (150):  CH2CO + H <=> HCCO + H2
    fwd_A[150]     = 50000000000000;
    fwd_beta[150]  = 0;
    fwd_Ea[150]    = 7999.5200000000004;
    prefactor_units[150]  = 1.0000000000000002e-06;
    activation_units[150] = 0.50321666580471969;
    phase_units[150]      = pow(10,-12.000000);
    is_PD[150] = 0;
    nTB[150] = 0;

    // (151):  CH2CO + TXCH2 <=> HCCO + CH3
    kiv[151] = {16,21,25,10};
    nuv[151] = {-1,-1,1,1};
    // (151):  CH2CO + TXCH2 <=> HCCO + CH3
    fwd_A[151]     = 36000000000000;
    fwd_beta[151]  = 0;
    fwd_Ea[151]    = 10999.040000000001;
    prefactor_units[151]  = 1.0000000000000002e-06;
    activation_units[151] = 0.50321666580471969;
    phase_units[151]      = pow(10,-12.000000);
    is_PD[151] = 0;
    nTB[151] = 0;

    // (152):  CH2CHO + O <=> CH2O + HCO
    kiv[152] = {27,3,11,20};
    nuv[152] = {-1,-1,1,1};
    // (152):  CH2CHO + O <=> CH2O + HCO
    fwd_A[152]     = 31700000000000;
    fwd_beta[152]  = 0.029999999999999999;
    fwd_Ea[152]    = -394.36000000000001;
    prefactor_units[152]  = 1.0000000000000002e-06;
    activation_units[152] = 0.50321666580471969;
    phase_units[152]      = pow(10,-12.000000);
    is_PD[152] = 0;
    nTB[152] = 0;

    // (153):  CH2CHO <=> CH2CO + H
    kiv[153] = {27,16,1};
    nuv[153] = {-1,1,1};
    // (153):  CH2CHO <=> CH2CO + H
    fwd_A[153]     = 1.3199999999999999e+34;
    fwd_beta[153]  = -6.5700000000000003;
    fwd_Ea[153]    = 49457.459999999999;
    prefactor_units[153]  = 1;
    activation_units[153] = 0.50321666580471969;
    phase_units[153]      = pow(10,-6.000000);
    is_PD[153] = 0;
    nTB[153] = 0;

    // (154):  CH2CHO + OH <=> H2O + CH2CO
    kiv[154] = {27,4,7,16};
    nuv[154] = {-1,-1,1,1};
    // (154):  CH2CHO + OH <=> H2O + CH2CO
    fwd_A[154]     = 12000000000000;
    fwd_beta[154]  = 0;
    fwd_Ea[154]    = 0;
    prefactor_units[154]  = 1.0000000000000002e-06;
    activation_units[154] = 0.50321666580471969;
    phase_units[154]      = pow(10,-12.000000);
    is_PD[154] = 0;
    nTB[154] = 0;

    // (155):  CH2CHO + H <=> CH2CO + H2
    kiv[155] = {27,1,16,2};
    nuv[155] = {-1,-1,1,1};
    // (155):  CH2CHO + H <=> CH2CO + H2
    fwd_A[155]     = 11000000000000;
    fwd_beta[155]  = 0;
    fwd_Ea[155]    = 0;
    prefactor_units[155]  = 1.0000000000000002e-06;
    activation_units[155] = 0.50321666580471969;
    phase_units[155]      = pow(10,-12.000000);
    is_PD[155] = 0;
    nTB[155] = 0;

    // (156):  CH2CHO + O2 => OH + CO + CH2O
    kiv[156] = {27,5,4,9,11};
    nuv[156] = {-1,-1,1,1,1};
    // (156):  CH2CHO + O2 => OH + CO + CH2O
    fwd_A[156]     = 18100000000;
    fwd_beta[156]  = 0;
    fwd_Ea[156]    = 0;
    prefactor_units[156]  = 1.0000000000000002e-06;
    activation_units[156] = 0.50321666580471969;
    phase_units[156]      = pow(10,-12.000000);
    is_PD[156] = 0;
    nTB[156] = 0;

    // (157):  CH2CHO <=> CH3 + CO
    kiv[157] = {27,10,9};
    nuv[157] = {-1,1,1};
    // (157):  CH2CHO <=> CH3 + CO
    fwd_A[157]     = 6.5100000000000001e+34;
    fwd_beta[157]  = -6.8700000000000001;
    fwd_Ea[157]    = 47194.07;
    prefactor_units[157]  = 1;
    activation_units[157] = 0.50321666580471969;
    phase_units[157]      = pow(10,-6.000000);
    is_PD[157] = 0;
    nTB[157] = 0;

    // (158):  CH2CHO + O2 => OH + 2.000000 HCO
    kiv[158] = {27,5,4,20};
    nuv[158] = {-1,-1,1,2.0};
    // (158):  CH2CHO + O2 => OH + 2.000000 HCO
    fwd_A[158]     = 23500000000;
    fwd_beta[158]  = 0;
    fwd_Ea[158]    = 0;
    prefactor_units[158]  = 1.0000000000000002e-06;
    activation_units[158] = 0.50321666580471969;
    phase_units[158]      = pow(10,-12.000000);
    is_PD[158] = 0;
    nTB[158] = 0;

    // (159):  CH2CHO + H <=> CH3 + HCO
    kiv[159] = {27,1,10,20};
    nuv[159] = {-1,-1,1,1};
    // (159):  CH2CHO + H <=> CH3 + HCO
    fwd_A[159]     = 22000000000000;
    fwd_beta[159]  = 0;
    fwd_Ea[159]    = 0;
    prefactor_units[159]  = 1.0000000000000002e-06;
    activation_units[159] = 0.50321666580471969;
    phase_units[159]      = pow(10,-12.000000);
    is_PD[159] = 0;
    nTB[159] = 0;

    // (160):  CH3CHO + O => CH3 + CO + OH
    kiv[160] = {26,3,10,9,4};
    nuv[160] = {-1,-1,1,1,1};
    // (160):  CH3CHO + O => CH3 + CO + OH
    fwd_A[160]     = 2920000000000;
    fwd_beta[160]  = 0;
    fwd_Ea[160]    = 1809.27;
    prefactor_units[160]  = 1.0000000000000002e-06;
    activation_units[160] = 0.50321666580471969;
    phase_units[160]      = pow(10,-12.000000);
    is_PD[160] = 0;
    nTB[160] = 0;

    // (161):  CH3CHO + O2 => CH3 + CO + HO2
    kiv[161] = {26,5,10,9,8};
    nuv[161] = {-1,-1,1,1,1};
    // (161):  CH3CHO + O2 => CH3 + CO + HO2
    fwd_A[161]     = 30100000000000;
    fwd_beta[161]  = 0;
    fwd_Ea[161]    = 39149.139999999999;
    prefactor_units[161]  = 1.0000000000000002e-06;
    activation_units[161] = 0.50321666580471969;
    phase_units[161]      = pow(10,-12.000000);
    is_PD[161] = 0;
    nTB[161] = 0;

    // (162):  CH3CHO + OH => CH3 + CO + H2O
    kiv[162] = {26,4,10,9,7};
    nuv[162] = {-1,-1,1,1,1};
    // (162):  CH3CHO + OH => CH3 + CO + H2O
    fwd_A[162]     = 23400000000;
    fwd_beta[162]  = 0.72999999999999998;
    fwd_Ea[162]    = -1113.77;
    prefactor_units[162]  = 1.0000000000000002e-06;
    activation_units[162] = 0.50321666580471969;
    phase_units[162]      = pow(10,-12.000000);
    is_PD[162] = 0;
    nTB[162] = 0;

    // (163):  CH3CHO + H <=> CH2CHO + H2
    kiv[163] = {26,1,27,2};
    nuv[163] = {-1,-1,1,1};
    // (163):  CH3CHO + H <=> CH2CHO + H2
    fwd_A[163]     = 2050000000;
    fwd_beta[163]  = 1.1599999999999999;
    fwd_Ea[163]    = 2404.4000000000001;
    prefactor_units[163]  = 1.0000000000000002e-06;
    activation_units[163] = 0.50321666580471969;
    phase_units[163]      = pow(10,-12.000000);
    is_PD[163] = 0;
    nTB[163] = 0;

    // (164):  CH3CHO + H => CH3 + CO + H2
    kiv[164] = {26,1,10,9,2};
    nuv[164] = {-1,-1,1,1,1};
    // (164):  CH3CHO + H => CH3 + CO + H2
    fwd_A[164]     = 2050000000;
    fwd_beta[164]  = 1.1599999999999999;
    fwd_Ea[164]    = 2404.4000000000001;
    prefactor_units[164]  = 1.0000000000000002e-06;
    activation_units[164] = 0.50321666580471969;
    phase_units[164]      = pow(10,-12.000000);
    is_PD[164] = 0;
    nTB[164] = 0;

    // (165):  CH3CHO + O <=> CH2CHO + OH
    kiv[165] = {26,3,27,4};
    nuv[165] = {-1,-1,1,1};
    // (165):  CH3CHO + O <=> CH2CHO + OH
    fwd_A[165]     = 2920000000000;
    fwd_beta[165]  = 0;
    fwd_Ea[165]    = 1809.27;
    prefactor_units[165]  = 1.0000000000000002e-06;
    activation_units[165] = 0.50321666580471969;
    phase_units[165]      = pow(10,-12.000000);
    is_PD[165] = 0;
    nTB[165] = 0;

    // (166):  CH3CHO + CH3 => CH3 + CO + CH4
    kiv[166] = {26,10,10,9,13};
    nuv[166] = {-1,-1,1,1,1};
    // (166):  CH3CHO + CH3 => CH3 + CO + CH4
    fwd_A[166]     = 2720000;
    fwd_beta[166]  = 1.77;
    fwd_Ea[166]    = 5920.1700000000001;
    prefactor_units[166]  = 1.0000000000000002e-06;
    activation_units[166] = 0.50321666580471969;
    phase_units[166]      = pow(10,-12.000000);
    is_PD[166] = 0;
    nTB[166] = 0;

    // (167):  CH3CHO + HO2 => CH3 + CO + H2O2
    kiv[167] = {26,8,10,9,6};
    nuv[167] = {-1,-1,1,1,1};
    // (167):  CH3CHO + HO2 => CH3 + CO + H2O2
    fwd_A[167]     = 3010000000000;
    fwd_beta[167]  = 0;
    fwd_Ea[167]    = 11924;
    prefactor_units[167]  = 1.0000000000000002e-06;
    activation_units[167] = 0.50321666580471969;
    phase_units[167]      = pow(10,-12.000000);
    is_PD[167] = 0;
    nTB[167] = 0;

    // (168):  C2H5O <=> CH3 + CH2O
    kiv[168] = {28,10,11};
    nuv[168] = {-1,1,1};
    // (168):  C2H5O <=> CH3 + CH2O
    fwd_A[168]     = 1.32e+20;
    fwd_beta[168]  = -2.02;
    fwd_Ea[168]    = 20750.48;
    prefactor_units[168]  = 1;
    activation_units[168] = 0.50321666580471969;
    phase_units[168]      = pow(10,-6.000000);
    is_PD[168] = 0;
    nTB[168] = 0;

    // (169):  C2H5O <=> CH3CHO + H
    kiv[169] = {28,26,1};
    nuv[169] = {-1,1,1};
    // (169):  C2H5O <=> CH3CHO + H
    fwd_A[169]     = 5450000000000000;
    fwd_beta[169]  = -0.68999999999999995;
    fwd_Ea[169]    = 22229.919999999998;
    prefactor_units[169]  = 1;
    activation_units[169] = 0.50321666580471969;
    phase_units[169]      = pow(10,-6.000000);
    is_PD[169] = 0;
    nTB[169] = 0;

    // (170):  C2H5O + O2 <=> CH3CHO + HO2
    kiv[170] = {28,5,26,8};
    nuv[170] = {-1,-1,1,1};
    // (170):  C2H5O + O2 <=> CH3CHO + HO2
    fwd_A[170]     = 22900000000;
    fwd_beta[170]  = 0;
    fwd_Ea[170]    = 874.75999999999999;
    prefactor_units[170]  = 1.0000000000000002e-06;
    activation_units[170] = 0.50321666580471969;
    phase_units[170]      = pow(10,-12.000000);
    is_PD[170] = 0;
    nTB[170] = 0;

    // (171):  SXCH2 + N2 <=> TXCH2 + N2
    kiv[171] = {22,0,21,0};
    nuv[171] = {-1,-1,1,1};
    // (171):  SXCH2 + N2 <=> TXCH2 + N2
    fwd_A[171]     = 15000000000000;
    fwd_beta[171]  = 0;
    fwd_Ea[171]    = 599.89999999999998;
    prefactor_units[171]  = 1.0000000000000002e-06;
    activation_units[171] = 0.50321666580471969;
    phase_units[171]      = pow(10,-12.000000);
    is_PD[171] = 0;
    nTB[171] = 0;

}


/* Finalizes parameter database */
void CKFINALIZE()
{
  for (int i=0; i<172; ++i) {
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
    awt[0] = 14.006700; /*N */
    awt[1] = 1.007970; /*H */
    awt[2] = 15.999400; /*O */
    awt[3] = 12.011150; /*C */

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
    for (id = 0; id < kd * 29; ++ id) {
         ncf[id] = 0; 
    }

    /*N2 */
    ncf[ 0 * kd + 0 ] = 2; /*N */

    /*H */
    ncf[ 1 * kd + 1 ] = 1; /*H */

    /*H2 */
    ncf[ 2 * kd + 1 ] = 2; /*H */

    /*O */
    ncf[ 3 * kd + 2 ] = 1; /*O */

    /*OH */
    ncf[ 4 * kd + 2 ] = 1; /*O */
    ncf[ 4 * kd + 1 ] = 1; /*H */

    /*O2 */
    ncf[ 5 * kd + 2 ] = 2; /*O */

    /*H2O2 */
    ncf[ 6 * kd + 1 ] = 2; /*H */
    ncf[ 6 * kd + 2 ] = 2; /*O */

    /*H2O */
    ncf[ 7 * kd + 1 ] = 2; /*H */
    ncf[ 7 * kd + 2 ] = 1; /*O */

    /*HO2 */
    ncf[ 8 * kd + 1 ] = 1; /*H */
    ncf[ 8 * kd + 2 ] = 2; /*O */

    /*CO */
    ncf[ 9 * kd + 3 ] = 1; /*C */
    ncf[ 9 * kd + 2 ] = 1; /*O */

    /*CH3 */
    ncf[ 10 * kd + 3 ] = 1; /*C */
    ncf[ 10 * kd + 1 ] = 3; /*H */

    /*CH2O */
    ncf[ 11 * kd + 1 ] = 2; /*H */
    ncf[ 11 * kd + 3 ] = 1; /*C */
    ncf[ 11 * kd + 2 ] = 1; /*O */

    /*CO2 */
    ncf[ 12 * kd + 3 ] = 1; /*C */
    ncf[ 12 * kd + 2 ] = 2; /*O */

    /*CH4 */
    ncf[ 13 * kd + 3 ] = 1; /*C */
    ncf[ 13 * kd + 1 ] = 4; /*H */

    /*C2H2 */
    ncf[ 14 * kd + 3 ] = 2; /*C */
    ncf[ 14 * kd + 1 ] = 2; /*H */

    /*C2H4 */
    ncf[ 15 * kd + 3 ] = 2; /*C */
    ncf[ 15 * kd + 1 ] = 4; /*H */

    /*CH2CO */
    ncf[ 16 * kd + 3 ] = 2; /*C */
    ncf[ 16 * kd + 1 ] = 2; /*H */
    ncf[ 16 * kd + 2 ] = 1; /*O */

    /*C2H6 */
    ncf[ 17 * kd + 3 ] = 2; /*C */
    ncf[ 17 * kd + 1 ] = 6; /*H */

    /*C */
    ncf[ 18 * kd + 3 ] = 1; /*C */

    /*CH */
    ncf[ 19 * kd + 3 ] = 1; /*C */
    ncf[ 19 * kd + 1 ] = 1; /*H */

    /*HCO */
    ncf[ 20 * kd + 1 ] = 1; /*H */
    ncf[ 20 * kd + 3 ] = 1; /*C */
    ncf[ 20 * kd + 2 ] = 1; /*O */

    /*TXCH2 */
    ncf[ 21 * kd + 3 ] = 1; /*C */
    ncf[ 21 * kd + 1 ] = 2; /*H */

    /*SXCH2 */
    ncf[ 22 * kd + 3 ] = 1; /*C */
    ncf[ 22 * kd + 1 ] = 2; /*H */

    /*C2H3 */
    ncf[ 23 * kd + 3 ] = 2; /*C */
    ncf[ 23 * kd + 1 ] = 3; /*H */

    /*C2H5 */
    ncf[ 24 * kd + 3 ] = 2; /*C */
    ncf[ 24 * kd + 1 ] = 5; /*H */

    /*HCCO */
    ncf[ 25 * kd + 1 ] = 1; /*H */
    ncf[ 25 * kd + 3 ] = 2; /*C */
    ncf[ 25 * kd + 2 ] = 1; /*O */

    /*CH3CHO */
    ncf[ 26 * kd + 1 ] = 4; /*H */
    ncf[ 26 * kd + 3 ] = 2; /*C */
    ncf[ 26 * kd + 2 ] = 1; /*O */

    /*CH2CHO */
    ncf[ 27 * kd + 1 ] = 3; /*H */
    ncf[ 27 * kd + 3 ] = 2; /*C */
    ncf[ 27 * kd + 2 ] = 1; /*O */

    /*C2H5O */
    ncf[ 28 * kd + 3 ] = 2; /*C */
    ncf[ 28 * kd + 1 ] = 5; /*H */
    ncf[ 28 * kd + 2 ] = 1; /*O */

}


/* Returns the vector of strings of element names */
void CKSYME_STR(amrex::Vector<std::string>& ename)
{
    ename.resize(4);
    ename[0] = "N";
    ename[1] = "H";
    ename[2] = "O";
    ename[3] = "C";
}


/* Returns the vector of strings of species names */
void CKSYMS_STR(amrex::Vector<std::string>& kname)
{
    kname.resize(29);
    kname[0] = "N2";
    kname[1] = "H";
    kname[2] = "H2";
    kname[3] = "O";
    kname[4] = "OH";
    kname[5] = "O2";
    kname[6] = "H2O2";
    kname[7] = "H2O";
    kname[8] = "HO2";
    kname[9] = "CO";
    kname[10] = "CH3";
    kname[11] = "CH2O";
    kname[12] = "CO2";
    kname[13] = "CH4";
    kname[14] = "C2H2";
    kname[15] = "C2H4";
    kname[16] = "CH2CO";
    kname[17] = "C2H6";
    kname[18] = "C";
    kname[19] = "CH";
    kname[20] = "HCO";
    kname[21] = "TXCH2";
    kname[22] = "SXCH2";
    kname[23] = "C2H3";
    kname[24] = "C2H5";
    kname[25] = "HCCO";
    kname[26] = "CH3CHO";
    kname[27] = "CH2CHO";
    kname[28] = "C2H5O";
}

/*compute the sparsity pattern of the chemistry Jacobian */
void SPARSITY_INFO( int * nJdata, int * consP, int NCELLS)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(900);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(29);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[900];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<29; k++) {
                c_d[k] = 1.0/ 29.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<30; k++) {
        for (int l=0; l<30; l++) {
            if(J_h[ 30 * k + l] != 0.0){
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(900);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(29);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[900];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<29; k++) {
                c_d[k] = 1.0/ 29.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<30; k++) {
        for (int l=0; l<30; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J_h[ 30 * k + l] != 0.0){
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(900);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(29);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[900];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<29; k++) {
                c_d[k] = 1.0/ 29.000000 ;
            }
            aJacobian_precond(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<30; k++) {
        for (int l=0; l<30; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J_h[ 30 * k + l] != 0.0){
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

    amrex::Gpu::DeviceVector<amrex::Real> J_v(900);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(29);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[900];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<29; k++) {
                c_d[k] = 1.0/ 29.000000 ;
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
        offset_row = nc * 30;
        offset_col = nc * 30;
        for (int k=0; k<30; k++) {
            for (int l=0; l<30; l++) {
                if(J_h[30*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(900);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(29);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[900];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<29; k++) {
                c_d[k] = 1.0/ 29.000000 ;
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
            offset = nc * 30;
            for (int l=0; l<30; l++) {
                for (int k=0; k<30; k++) {
                    if(J_h[30*k + l] != 0.0) {
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
            offset = nc * 30;
            for (int l=0; l<30; l++) {
                for (int k=0; k<30; k++) {
                    if(J_h[30*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(900);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(29);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[900];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<29; k++) {
                c_d[k] = 1.0/ 29.000000 ;
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
            offset = nc * 30;
            for (int l=0; l<30; l++) {
                for (int k=0; k<30; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp-1] = l+1 + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(J_h[30*k + l] != 0.0) {
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
            offset = nc * 30;
            for (int l=0; l<30; l++) {
                for (int k=0; k<30; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp] = l + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(J_h[30*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(900);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(29);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[900];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<29; k++) {
                c_d[k] = 1.0/ 29.000000 ;
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
    for (int k=0; k<30; k++) {
        for (int l=0; l<30; l++) {
            if (k == l) {
                rowVals[nJdata_tmp] = l; 
                indx[nJdata_tmp] = 30*k + l;
                nJdata_tmp = nJdata_tmp + 1; 
            } else {
                if(J_h[30*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l; 
                    indx[nJdata_tmp] = 30*k + l;
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(900);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(29);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[900];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<29; k++) {
                c_d[k] = 1.0/ 29.000000 ;
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
        for (int l=0; l<30; l++) {
            for (int k=0; k<30; k++) {
                if (k == l) {
                    colVals[nJdata_tmp-1] = l+1; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J_h[30*k + l] != 0.0) {
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
        for (int l=0; l<30; l++) {
            for (int k=0; k<30; k++) {
                if (k == l) {
                    colVals[nJdata_tmp] = l; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J_h[30*k + l] != 0.0) {
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

