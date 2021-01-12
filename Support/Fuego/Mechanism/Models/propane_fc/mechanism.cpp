#ifndef MECHANISM_CPP
#define MECHANISM_CPP

#include "chemistry_file.H"

/* PURE CPU stuff  */
#ifndef AMREX_USE_GPU
namespace thermo
{
    amrex::Real fwd_A[167], fwd_beta[167], fwd_Ea[167];
    amrex::Real low_A[167], low_beta[167], low_Ea[167];
    amrex::Real rev_A[167], rev_beta[167], rev_Ea[167];
    amrex::Real troe_a[167],troe_Ts[167], troe_Tss[167], troe_Tsss[167];
    amrex::Real sri_a[167], sri_b[167], sri_c[167], sri_d[167], sri_e[167];
    amrex::Real activation_units[167], prefactor_units[167], phase_units[167];
    int is_PD[167], troe_len[167], sri_len[167], nTB[167], *TBid[167];
    amrex::Real *TB[167];
    std::vector<std::vector<amrex::Real>> kiv(167); 
    std::vector<std::vector<amrex::Real>> nuv(167); 
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
        if (*i > 167) {
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
    amrex::Real imw[34];

    get_imw(imw);

    for (int i=0; i<(*np); i++) {
        YOW[i] = 0.0;
    }

    for (int n=0; n<34; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] = y[n*(*np)+i] * imw[n];
            YOW[i] += x[n*(*np)+i];
        }
    }

    for (int i=0; i<(*np); i++) {
        YOW[i] = 1.0/YOW[i];
    }

    for (int n=0; n<34; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] *=  YOW[i];
        }
    }
}

/*Returns enthalpy in mass units (Eq 27.) */
void VCKHMS(int *  np, amrex::Real *  T,  amrex::Real *  hms)
{
    amrex::Real tc[5], h[34];
    amrex::Real imw[34];

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
    }

    for (int n=0; n<34; n++) {
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
    amrex::Real c[34*(*np)]; /*temporary storage */
    amrex::Real imw[34];

    get_imw(imw);

    /*See Eq 8 with an extra 1e6 so c goes to SI */
    for (int n=0; n<34; n++) {
        for (int i=0; i<(*np); i++) {
            c[n*(*np)+i] = 1.0e6 * rho[i] * y[n*(*np)+i] * imw[n];
        }
    }

    /*call productionRate */
    vproductionRate(*np, wdot, c, T);

    /*convert to chemkin units */
    for (int i=0; i<34*(*np); i++) {
        wdot[i] *= 1.0e-6;
    }
}


/*Compute P = rhoRT/W(y) */
void VCKPY(int *  np, amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  P)
{
    amrex::Real YOW[*np];
    amrex::Real imw[34];

    get_imw(imw);

    for (int i=0; i<(*np); i++) {
        YOW[i] = 0.0;
    }

    for (int n=0; n<34; n++) {
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
    amrex::Real k_f_s[167*npt], Kc_s[167*npt], mixture[npt], g_RT[34*npt];
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

    for (int n=0; n<34; n++) {
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
    vcomp_wdot_151_167(npt, wdot, mixture, sc, k_f_s, Kc_s, tc, invT, T);
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
    }
}

void vcomp_gibbs(int npt, amrex::Real *  g_RT, amrex::Real *  tc)
{
    /*compute the Gibbs free energy */
    for (int i=0; i<npt; i++) {
        amrex::Real tg[5], g[34];
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
    }
}

void vcomp_Kc(int npt, amrex::Real *  Kc_s, amrex::Real *  g_RT, amrex::Real *  invT)
{
    for (int i=0; i<npt; i++) {
        /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
        amrex::Real refC = (101325. / 8.31451) * invT[i];
        amrex::Real refCinv = 1.0 / refC;

        Kc_s[0*npt+i] = refCinv * exp((2.000000 * g_RT[3*npt+i]) - (g_RT[7*npt+i]));
        Kc_s[1*npt+i] = refCinv * exp((g_RT[2*npt+i] + g_RT[5*npt+i]) - (g_RT[8*npt+i]));
        Kc_s[2*npt+i] = refCinv * exp((g_RT[2*npt+i] + g_RT[11*npt+i]) - (g_RT[15*npt+i]));
        Kc_s[3*npt+i] = refCinv * exp((2.000000 * g_RT[11*npt+i]) - (g_RT[17*npt+i]));
        Kc_s[4*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[12*npt+i]) - (g_RT[18*npt+i]));
        Kc_s[5*npt+i] = refC * exp((g_RT[13*npt+i]) - (g_RT[2*npt+i] + g_RT[10*npt+i]));
        Kc_s[6*npt+i] = refC * exp((g_RT[24*npt+i]) - (g_RT[2*npt+i] + g_RT[22*npt+i]));
        Kc_s[7*npt+i] = refCinv * exp((g_RT[2*npt+i] + g_RT[16*npt+i]) - (g_RT[14*npt+i]));
        Kc_s[8*npt+i] = refC * exp((g_RT[16*npt+i]) - (g_RT[4*npt+i] + g_RT[22*npt+i]));
        Kc_s[9*npt+i] = refCinv * exp((g_RT[1*npt+i] + g_RT[2*npt+i]) - (g_RT[3*npt+i]));
        Kc_s[10*npt+i] = refCinv * exp((2.000000 * g_RT[1*npt+i]) - (g_RT[5*npt+i]));
        Kc_s[11*npt+i] = refCinv * exp((g_RT[2*npt+i] + g_RT[3*npt+i]) - (g_RT[6*npt+i]));
        Kc_s[12*npt+i] = refC * exp((g_RT[19*npt+i]) - (g_RT[2*npt+i] + g_RT[12*npt+i]));
        Kc_s[13*npt+i] = refCinv * exp((g_RT[2*npt+i] + g_RT[12*npt+i]) - (g_RT[19*npt+i]));
        Kc_s[14*npt+i] = refC * exp((g_RT[20*npt+i]) - (g_RT[5*npt+i] + g_RT[11*npt+i]));
        Kc_s[15*npt+i] = refCinv * exp((g_RT[5*npt+i] + g_RT[11*npt+i]) - (g_RT[20*npt+i]));
        Kc_s[16*npt+i] = exp((g_RT[1*npt+i] + g_RT[4*npt+i]) - (g_RT[2*npt+i] + g_RT[3*npt+i]));
        Kc_s[17*npt+i] = exp((g_RT[2*npt+i] + g_RT[3*npt+i]) - (g_RT[1*npt+i] + g_RT[4*npt+i]));
        Kc_s[18*npt+i] = exp((g_RT[3*npt+i] + g_RT[4*npt+i]) - (g_RT[2*npt+i] + g_RT[6*npt+i]));
        Kc_s[19*npt+i] = exp((g_RT[2*npt+i] + g_RT[6*npt+i]) - (g_RT[3*npt+i] + g_RT[4*npt+i]));
        Kc_s[20*npt+i] = exp((g_RT[1*npt+i] + g_RT[6*npt+i]) - (2.000000 * g_RT[3*npt+i]));
        Kc_s[21*npt+i] = exp((2.000000 * g_RT[3*npt+i]) - (g_RT[1*npt+i] + g_RT[6*npt+i]));
        Kc_s[22*npt+i] = exp((g_RT[2*npt+i] + g_RT[5*npt+i]) - (g_RT[1*npt+i] + g_RT[3*npt+i]));
        Kc_s[23*npt+i] = exp((g_RT[1*npt+i] + g_RT[3*npt+i]) - (g_RT[2*npt+i] + g_RT[5*npt+i]));
        Kc_s[24*npt+i] = exp((g_RT[2*npt+i] + g_RT[8*npt+i]) - (2.000000 * g_RT[3*npt+i]));
        Kc_s[25*npt+i] = exp((2.000000 * g_RT[8*npt+i]) - (g_RT[5*npt+i] + g_RT[7*npt+i]));
        Kc_s[26*npt+i] = exp((g_RT[2*npt+i] + g_RT[8*npt+i]) - (g_RT[4*npt+i] + g_RT[5*npt+i]));
        Kc_s[27*npt+i] = exp((g_RT[3*npt+i] + g_RT[8*npt+i]) - (g_RT[5*npt+i] + g_RT[6*npt+i]));
        Kc_s[28*npt+i] = exp((g_RT[5*npt+i] + g_RT[6*npt+i]) - (g_RT[3*npt+i] + g_RT[8*npt+i]));
        Kc_s[29*npt+i] = exp((2.000000 * g_RT[8*npt+i]) - (g_RT[5*npt+i] + g_RT[7*npt+i]));
        Kc_s[30*npt+i] = exp((g_RT[1*npt+i] + g_RT[8*npt+i]) - (g_RT[3*npt+i] + g_RT[5*npt+i]));
        Kc_s[31*npt+i] = exp((g_RT[3*npt+i] + g_RT[5*npt+i]) - (g_RT[1*npt+i] + g_RT[8*npt+i]));
        Kc_s[32*npt+i] = exp((g_RT[3*npt+i] + g_RT[7*npt+i]) - (g_RT[6*npt+i] + g_RT[8*npt+i]));
        Kc_s[33*npt+i] = exp((g_RT[2*npt+i] + g_RT[7*npt+i]) - (g_RT[3*npt+i] + g_RT[6*npt+i]));
        Kc_s[34*npt+i] = exp((g_RT[3*npt+i] + g_RT[7*npt+i]) - (g_RT[6*npt+i] + g_RT[8*npt+i]));
        Kc_s[35*npt+i] = exp((g_RT[6*npt+i] + g_RT[8*npt+i]) - (g_RT[3*npt+i] + g_RT[7*npt+i]));
        Kc_s[36*npt+i] = exp((g_RT[1*npt+i] + g_RT[7*npt+i]) - (g_RT[3*npt+i] + g_RT[8*npt+i]));
        Kc_s[37*npt+i] = exp((g_RT[2*npt+i] + g_RT[7*npt+i]) - (g_RT[4*npt+i] + g_RT[8*npt+i]));
        Kc_s[38*npt+i] = exp((g_RT[4*npt+i] + g_RT[8*npt+i]) - (g_RT[2*npt+i] + g_RT[7*npt+i]));
        Kc_s[39*npt+i] = exp((g_RT[3*npt+i] + g_RT[9*npt+i]) - (g_RT[2*npt+i] + g_RT[10*npt+i]));
        Kc_s[40*npt+i] = exp((g_RT[4*npt+i] + g_RT[9*npt+i]) - (g_RT[2*npt+i] + g_RT[11*npt+i]));
        Kc_s[41*npt+i] = exp((g_RT[2*npt+i] + g_RT[11*npt+i]) - (g_RT[4*npt+i] + g_RT[9*npt+i]));
        Kc_s[42*npt+i] = refC * exp((g_RT[5*npt+i] + g_RT[9*npt+i]) - (g_RT[2*npt+i] + g_RT[3*npt+i] + g_RT[12*npt+i]));
        Kc_s[43*npt+i] = refC * exp((g_RT[1*npt+i] + g_RT[9*npt+i]) - (2.000000 * g_RT[2*npt+i] + g_RT[12*npt+i]));
        Kc_s[44*npt+i] = exp((g_RT[8*npt+i] + g_RT[11*npt+i]) - (g_RT[3*npt+i] + g_RT[13*npt+i]));
        Kc_s[45*npt+i] = exp((g_RT[5*npt+i] + g_RT[11*npt+i]) - (g_RT[3*npt+i] + g_RT[10*npt+i]));
        Kc_s[46*npt+i] = exp((g_RT[5*npt+i] + g_RT[11*npt+i]) - (g_RT[1*npt+i] + g_RT[13*npt+i]));
        Kc_s[47*npt+i] = exp((g_RT[1*npt+i] + g_RT[13*npt+i]) - (g_RT[5*npt+i] + g_RT[11*npt+i]));
        Kc_s[48*npt+i] = exp((2.000000 * g_RT[11*npt+i]) - (g_RT[2*npt+i] + g_RT[14*npt+i]));
        Kc_s[49*npt+i] = exp((g_RT[8*npt+i] + g_RT[11*npt+i]) - (g_RT[5*npt+i] + g_RT[15*npt+i]));
        Kc_s[50*npt+i] = exp((g_RT[1*npt+i] + g_RT[11*npt+i]) - (g_RT[2*npt+i] + g_RT[10*npt+i]));
        Kc_s[51*npt+i] = exp((g_RT[3*npt+i] + g_RT[11*npt+i]) - (g_RT[4*npt+i] + g_RT[10*npt+i]));
        Kc_s[52*npt+i] = exp((g_RT[3*npt+i] + g_RT[11*npt+i]) - (g_RT[6*npt+i] + g_RT[9*npt+i]));
        Kc_s[53*npt+i] = exp((g_RT[6*npt+i] + g_RT[9*npt+i]) - (g_RT[3*npt+i] + g_RT[11*npt+i]));
        Kc_s[54*npt+i] = exp((g_RT[7*npt+i] + g_RT[11*npt+i]) - (g_RT[8*npt+i] + g_RT[15*npt+i]));
        Kc_s[55*npt+i] = exp((g_RT[9*npt+i] + g_RT[11*npt+i]) - (g_RT[2*npt+i] + g_RT[16*npt+i]));
        Kc_s[56*npt+i] = exp((g_RT[9*npt+i] + g_RT[15*npt+i]) - (2.000000 * g_RT[11*npt+i]));
        Kc_s[57*npt+i] = exp((2.000000 * g_RT[11*npt+i]) - (g_RT[9*npt+i] + g_RT[15*npt+i]));
        Kc_s[58*npt+i] = exp((g_RT[1*npt+i] + g_RT[15*npt+i]) - (g_RT[3*npt+i] + g_RT[11*npt+i]));
        Kc_s[59*npt+i] = exp((g_RT[3*npt+i] + g_RT[11*npt+i]) - (g_RT[1*npt+i] + g_RT[15*npt+i]));
        Kc_s[60*npt+i] = exp((g_RT[2*npt+i] + g_RT[15*npt+i]) - (g_RT[4*npt+i] + g_RT[11*npt+i]));
        Kc_s[61*npt+i] = exp((g_RT[4*npt+i] + g_RT[11*npt+i]) - (g_RT[2*npt+i] + g_RT[15*npt+i]));
        Kc_s[62*npt+i] = exp((g_RT[3*npt+i] + g_RT[15*npt+i]) - (g_RT[6*npt+i] + g_RT[11*npt+i]));
        Kc_s[63*npt+i] = exp((g_RT[6*npt+i] + g_RT[11*npt+i]) - (g_RT[3*npt+i] + g_RT[15*npt+i]));
        Kc_s[64*npt+i] = exp((g_RT[8*npt+i] + g_RT[12*npt+i]) - (g_RT[3*npt+i] + g_RT[18*npt+i]));
        Kc_s[65*npt+i] = exp((g_RT[5*npt+i] + g_RT[12*npt+i]) - (g_RT[1*npt+i] + g_RT[18*npt+i]));
        Kc_s[66*npt+i] = exp((g_RT[1*npt+i] + g_RT[18*npt+i]) - (g_RT[5*npt+i] + g_RT[12*npt+i]));
        Kc_s[67*npt+i] = exp((g_RT[3*npt+i] + g_RT[12*npt+i]) - (g_RT[2*npt+i] + g_RT[18*npt+i]));
        Kc_s[68*npt+i] = exp((g_RT[2*npt+i] + g_RT[18*npt+i]) - (g_RT[3*npt+i] + g_RT[12*npt+i]));
        Kc_s[69*npt+i] = exp((g_RT[11*npt+i] + g_RT[19*npt+i]) - (g_RT[12*npt+i] + g_RT[15*npt+i]));
        Kc_s[70*npt+i] = exp((g_RT[2*npt+i] + g_RT[19*npt+i]) - (g_RT[4*npt+i] + g_RT[12*npt+i]));
        Kc_s[71*npt+i] = exp((g_RT[5*npt+i] + g_RT[19*npt+i]) - (g_RT[8*npt+i] + g_RT[12*npt+i]));
        Kc_s[72*npt+i] = exp((g_RT[1*npt+i] + g_RT[19*npt+i]) - (g_RT[3*npt+i] + g_RT[12*npt+i]));
        Kc_s[73*npt+i] = exp((g_RT[1*npt+i] + g_RT[19*npt+i]) - (g_RT[2*npt+i] + g_RT[18*npt+i]));
        Kc_s[74*npt+i] = exp((g_RT[3*npt+i] + g_RT[19*npt+i]) - (g_RT[6*npt+i] + g_RT[12*npt+i]));
        Kc_s[75*npt+i] = exp((g_RT[2*npt+i] + g_RT[10*npt+i]) - (g_RT[4*npt+i] + g_RT[19*npt+i]));
        Kc_s[76*npt+i] = exp((g_RT[5*npt+i] + g_RT[10*npt+i]) - (g_RT[8*npt+i] + g_RT[19*npt+i]));
        Kc_s[77*npt+i] = exp((g_RT[3*npt+i] + g_RT[10*npt+i]) - (g_RT[6*npt+i] + g_RT[19*npt+i]));
        Kc_s[78*npt+i] = exp((g_RT[8*npt+i] + g_RT[10*npt+i]) - (g_RT[7*npt+i] + g_RT[19*npt+i]));
        Kc_s[79*npt+i] = exp((g_RT[1*npt+i] + g_RT[10*npt+i]) - (g_RT[3*npt+i] + g_RT[19*npt+i]));
        Kc_s[80*npt+i] = exp((g_RT[10*npt+i] + g_RT[11*npt+i]) - (g_RT[15*npt+i] + g_RT[19*npt+i]));
        Kc_s[81*npt+i] = exp((g_RT[5*npt+i] + g_RT[13*npt+i]) - (g_RT[8*npt+i] + g_RT[10*npt+i]));
        Kc_s[82*npt+i] = exp((g_RT[9*npt+i] + g_RT[18*npt+i]) - (g_RT[10*npt+i] + g_RT[12*npt+i]));
        Kc_s[83*npt+i] = exp((g_RT[11*npt+i] + g_RT[20*npt+i]) - (2.000000 * g_RT[13*npt+i]));
        Kc_s[84*npt+i] = exp((g_RT[8*npt+i] + g_RT[20*npt+i]) - (g_RT[5*npt+i] + g_RT[21*npt+i]));
        Kc_s[85*npt+i] = refC * exp((2.000000 * g_RT[20*npt+i]) - (g_RT[5*npt+i] + 2.000000 * g_RT[13*npt+i]));
        Kc_s[86*npt+i] = exp((g_RT[10*npt+i] + g_RT[20*npt+i]) - (g_RT[19*npt+i] + g_RT[21*npt+i]));
        Kc_s[87*npt+i] = refC * exp((g_RT[21*npt+i]) - (g_RT[3*npt+i] + g_RT[13*npt+i]));
        Kc_s[88*npt+i] = refCinv * exp((g_RT[3*npt+i] + g_RT[13*npt+i]) - (g_RT[21*npt+i]));
        Kc_s[89*npt+i] = exp((g_RT[5*npt+i] + g_RT[22*npt+i]) - (g_RT[3*npt+i] + g_RT[23*npt+i]));
        Kc_s[90*npt+i] = exp((g_RT[1*npt+i] + g_RT[22*npt+i]) - (g_RT[2*npt+i] + g_RT[23*npt+i]));
        Kc_s[91*npt+i] = exp((g_RT[5*npt+i] + g_RT[24*npt+i]) - (g_RT[1*npt+i] + g_RT[25*npt+i]));
        Kc_s[92*npt+i] = exp((g_RT[2*npt+i] + g_RT[24*npt+i]) - (g_RT[4*npt+i] + g_RT[22*npt+i]));
        Kc_s[93*npt+i] = exp((g_RT[5*npt+i] + g_RT[24*npt+i]) - (g_RT[10*npt+i] + g_RT[19*npt+i]));
        Kc_s[94*npt+i] = exp((g_RT[5*npt+i] + g_RT[24*npt+i]) - (g_RT[8*npt+i] + g_RT[22*npt+i]));
        Kc_s[95*npt+i] = exp((g_RT[1*npt+i] + g_RT[16*npt+i]) - (g_RT[11*npt+i] + g_RT[19*npt+i]));
        Kc_s[96*npt+i] = exp((g_RT[3*npt+i] + g_RT[16*npt+i]) - (g_RT[6*npt+i] + g_RT[24*npt+i]));
        Kc_s[97*npt+i] = exp((g_RT[6*npt+i] + g_RT[24*npt+i]) - (g_RT[3*npt+i] + g_RT[16*npt+i]));
        Kc_s[98*npt+i] = exp((g_RT[2*npt+i] + g_RT[16*npt+i]) - (g_RT[4*npt+i] + g_RT[24*npt+i]));
        Kc_s[99*npt+i] = exp((g_RT[4*npt+i] + g_RT[24*npt+i]) - (g_RT[2*npt+i] + g_RT[16*npt+i]));
        Kc_s[100*npt+i] = exp((g_RT[16*npt+i] + g_RT[20*npt+i]) - (g_RT[21*npt+i] + g_RT[24*npt+i]));
        Kc_s[101*npt+i] = exp((g_RT[11*npt+i] + g_RT[16*npt+i]) - (g_RT[15*npt+i] + g_RT[24*npt+i]));
        Kc_s[102*npt+i] = exp((g_RT[15*npt+i] + g_RT[24*npt+i]) - (g_RT[11*npt+i] + g_RT[16*npt+i]));
        Kc_s[103*npt+i] = exp((g_RT[1*npt+i] + g_RT[16*npt+i]) - (g_RT[2*npt+i] + g_RT[25*npt+i]));
        Kc_s[104*npt+i] = exp((g_RT[11*npt+i] + g_RT[14*npt+i]) - (g_RT[15*npt+i] + g_RT[16*npt+i]));
        Kc_s[105*npt+i] = exp((g_RT[5*npt+i] + g_RT[14*npt+i]) - (g_RT[8*npt+i] + g_RT[16*npt+i]));
        Kc_s[106*npt+i] = exp((g_RT[8*npt+i] + g_RT[16*npt+i]) - (g_RT[5*npt+i] + g_RT[14*npt+i]));
        Kc_s[107*npt+i] = exp((g_RT[8*npt+i] + g_RT[14*npt+i]) - (g_RT[5*npt+i] + g_RT[17*npt+i]));
        Kc_s[108*npt+i] = refCinv * exp((g_RT[2*npt+i] + g_RT[14*npt+i]) - (g_RT[17*npt+i]));
        Kc_s[109*npt+i] = exp((g_RT[2*npt+i] + g_RT[17*npt+i]) - (g_RT[4*npt+i] + g_RT[14*npt+i]));
        Kc_s[110*npt+i] = exp((g_RT[4*npt+i] + g_RT[14*npt+i]) - (g_RT[2*npt+i] + g_RT[17*npt+i]));
        Kc_s[111*npt+i] = exp((g_RT[3*npt+i] + g_RT[17*npt+i]) - (g_RT[6*npt+i] + g_RT[14*npt+i]));
        Kc_s[112*npt+i] = exp((g_RT[9*npt+i] + g_RT[17*npt+i]) - (g_RT[11*npt+i] + g_RT[14*npt+i]));
        Kc_s[113*npt+i] = exp((g_RT[1*npt+i] + g_RT[17*npt+i]) - (g_RT[3*npt+i] + g_RT[14*npt+i]));
        Kc_s[114*npt+i] = exp((g_RT[11*npt+i] + g_RT[17*npt+i]) - (g_RT[14*npt+i] + g_RT[15*npt+i]));
        Kc_s[115*npt+i] = refC * exp((g_RT[1*npt+i] + g_RT[23*npt+i]) - (g_RT[2*npt+i] + 2.000000 * g_RT[12*npt+i]));
        Kc_s[116*npt+i] = exp((g_RT[5*npt+i] + g_RT[23*npt+i]) - (g_RT[18*npt+i] + g_RT[19*npt+i]));
        Kc_s[117*npt+i] = exp((g_RT[3*npt+i] + g_RT[23*npt+i]) - (2.000000 * g_RT[19*npt+i]));
        Kc_s[118*npt+i] = exp((g_RT[2*npt+i] + g_RT[23*npt+i]) - (g_RT[9*npt+i] + g_RT[12*npt+i]));
        Kc_s[119*npt+i] = exp((g_RT[9*npt+i] + g_RT[12*npt+i]) - (g_RT[2*npt+i] + g_RT[23*npt+i]));
        Kc_s[120*npt+i] = refC * exp((g_RT[5*npt+i] + g_RT[25*npt+i]) - (g_RT[3*npt+i] + g_RT[10*npt+i] + g_RT[12*npt+i]));
        Kc_s[121*npt+i] = exp((g_RT[10*npt+i] + g_RT[26*npt+i]) - (g_RT[19*npt+i] + g_RT[27*npt+i]));
        Kc_s[122*npt+i] = exp((g_RT[19*npt+i] + g_RT[27*npt+i]) - (g_RT[10*npt+i] + g_RT[26*npt+i]));
        Kc_s[123*npt+i] = exp((g_RT[8*npt+i] + g_RT[26*npt+i]) - (g_RT[3*npt+i] + g_RT[28*npt+i]));
        Kc_s[124*npt+i] = exp((g_RT[20*npt+i] + g_RT[26*npt+i]) - (g_RT[13*npt+i] + g_RT[28*npt+i]));
        Kc_s[125*npt+i] = exp((g_RT[8*npt+i] + g_RT[26*npt+i]) - (g_RT[5*npt+i] + g_RT[27*npt+i]));
        Kc_s[126*npt+i] = exp((g_RT[5*npt+i] + g_RT[26*npt+i]) - (g_RT[10*npt+i] + g_RT[25*npt+i]));
        Kc_s[127*npt+i] = refC * exp((g_RT[5*npt+i] + g_RT[26*npt+i]) - (g_RT[3*npt+i] + g_RT[10*npt+i] + g_RT[22*npt+i]));
        Kc_s[128*npt+i] = refCinv * exp((g_RT[2*npt+i] + g_RT[26*npt+i]) - (g_RT[27*npt+i]));
        Kc_s[129*npt+i] = refC * exp((g_RT[26*npt+i]) - (g_RT[11*npt+i] + g_RT[22*npt+i]));
        Kc_s[130*npt+i] = refC * exp((g_RT[27*npt+i]) - (g_RT[11*npt+i] + g_RT[24*npt+i]));
        Kc_s[131*npt+i] = refCinv * exp((g_RT[11*npt+i] + g_RT[24*npt+i]) - (g_RT[27*npt+i]));
        Kc_s[132*npt+i] = exp((g_RT[1*npt+i] + g_RT[27*npt+i]) - (g_RT[14*npt+i] + g_RT[19*npt+i]));
        Kc_s[133*npt+i] = exp((g_RT[2*npt+i] + g_RT[27*npt+i]) - (g_RT[11*npt+i] + g_RT[16*npt+i]));
        Kc_s[134*npt+i] = exp((g_RT[11*npt+i] + g_RT[16*npt+i]) - (g_RT[2*npt+i] + g_RT[27*npt+i]));
        Kc_s[135*npt+i] = exp((g_RT[2*npt+i] + g_RT[27*npt+i]) - (g_RT[4*npt+i] + g_RT[26*npt+i]));
        Kc_s[136*npt+i] = exp((g_RT[4*npt+i] + g_RT[26*npt+i]) - (g_RT[2*npt+i] + g_RT[27*npt+i]));
        Kc_s[137*npt+i] = exp((g_RT[3*npt+i] + g_RT[27*npt+i]) - (g_RT[6*npt+i] + g_RT[26*npt+i]));
        Kc_s[138*npt+i] = exp((g_RT[1*npt+i] + g_RT[27*npt+i]) - (g_RT[3*npt+i] + g_RT[26*npt+i]));
        Kc_s[139*npt+i] = exp((g_RT[11*npt+i] + g_RT[27*npt+i]) - (g_RT[15*npt+i] + g_RT[26*npt+i]));
        Kc_s[140*npt+i] = exp((g_RT[5*npt+i] + g_RT[29*npt+i]) - (g_RT[8*npt+i] + g_RT[27*npt+i]));
        Kc_s[141*npt+i] = refC * exp((g_RT[29*npt+i]) - (g_RT[2*npt+i] + g_RT[27*npt+i]));
        Kc_s[142*npt+i] = refCinv * exp((g_RT[2*npt+i] + g_RT[27*npt+i]) - (g_RT[29*npt+i]));
        Kc_s[143*npt+i] = refC * exp((g_RT[30*npt+i]) - (g_RT[11*npt+i] + g_RT[16*npt+i]));
        Kc_s[144*npt+i] = refCinv * exp((g_RT[11*npt+i] + g_RT[16*npt+i]) - (g_RT[30*npt+i]));
        Kc_s[145*npt+i] = exp((g_RT[8*npt+i] + g_RT[30*npt+i]) - (g_RT[5*npt+i] + g_RT[31*npt+i]));
        Kc_s[146*npt+i] = exp((g_RT[5*npt+i] + g_RT[30*npt+i]) - (g_RT[8*npt+i] + g_RT[27*npt+i]));
        Kc_s[147*npt+i] = refC * exp((g_RT[30*npt+i]) - (g_RT[2*npt+i] + g_RT[27*npt+i]));
        Kc_s[148*npt+i] = refCinv * exp((g_RT[2*npt+i] + g_RT[27*npt+i]) - (g_RT[30*npt+i]));
        Kc_s[149*npt+i] = exp((g_RT[3*npt+i] + g_RT[31*npt+i]) - (g_RT[6*npt+i] + g_RT[30*npt+i]));
        Kc_s[150*npt+i] = exp((g_RT[8*npt+i] + g_RT[31*npt+i]) - (g_RT[7*npt+i] + g_RT[30*npt+i]));
        Kc_s[151*npt+i] = exp((g_RT[2*npt+i] + g_RT[31*npt+i]) - (g_RT[4*npt+i] + g_RT[30*npt+i]));
        Kc_s[152*npt+i] = exp((g_RT[3*npt+i] + g_RT[31*npt+i]) - (g_RT[6*npt+i] + g_RT[29*npt+i]));
        Kc_s[153*npt+i] = exp((g_RT[11*npt+i] + g_RT[31*npt+i]) - (g_RT[15*npt+i] + g_RT[29*npt+i]));
        Kc_s[154*npt+i] = exp((g_RT[11*npt+i] + g_RT[31*npt+i]) - (g_RT[15*npt+i] + g_RT[30*npt+i]));
        Kc_s[155*npt+i] = exp((g_RT[1*npt+i] + g_RT[31*npt+i]) - (g_RT[3*npt+i] + g_RT[29*npt+i]));
        Kc_s[156*npt+i] = exp((g_RT[8*npt+i] + g_RT[31*npt+i]) - (g_RT[7*npt+i] + g_RT[29*npt+i]));
        Kc_s[157*npt+i] = exp((g_RT[1*npt+i] + g_RT[31*npt+i]) - (g_RT[3*npt+i] + g_RT[30*npt+i]));
        Kc_s[158*npt+i] = exp((g_RT[5*npt+i] + g_RT[31*npt+i]) - (g_RT[8*npt+i] + g_RT[29*npt+i]));
        Kc_s[159*npt+i] = exp((g_RT[8*npt+i] + g_RT[29*npt+i]) - (g_RT[5*npt+i] + g_RT[31*npt+i]));
        Kc_s[160*npt+i] = exp((g_RT[2*npt+i] + g_RT[31*npt+i]) - (g_RT[4*npt+i] + g_RT[29*npt+i]));
        Kc_s[161*npt+i] = exp((g_RT[4*npt+i] + g_RT[29*npt+i]) - (g_RT[2*npt+i] + g_RT[31*npt+i]));
        Kc_s[162*npt+i] = refC * exp((g_RT[28*npt+i]) - (g_RT[10*npt+i] + g_RT[24*npt+i]));
        Kc_s[163*npt+i] = refC * exp((g_RT[32*npt+i]) - (g_RT[5*npt+i] + g_RT[29*npt+i]));
        Kc_s[164*npt+i] = refCinv * exp((g_RT[5*npt+i] + g_RT[29*npt+i]) - (g_RT[32*npt+i]));
        Kc_s[165*npt+i] = refC * exp((g_RT[33*npt+i]) - (g_RT[5*npt+i] + g_RT[30*npt+i]));
        Kc_s[166*npt+i] = refCinv * exp((g_RT[5*npt+i] + g_RT[30*npt+i]) - (g_RT[33*npt+i]));
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
        phi_f = pow(sc[3*npt+i], 2.000000);
        alpha = mixture[i] + (TB[0][0] - 1)*sc[4*npt+i] + (TB[0][1] - 1)*sc[6*npt+i] + (TB[0][2] - 1)*sc[18*npt+i] + (TB[0][3] - 1)*sc[12*npt+i];
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
        wdot[3*npt+i] -= 2.000000 * qdot;
        wdot[7*npt+i] += qdot;

        /*reaction 2: H + O2 (+M) <=> HO2 (+M) */
        phi_f = sc[2*npt+i]*sc[5*npt+i];
        alpha = mixture[i] + (TB[1][0] - 1)*sc[4*npt+i] + (TB[1][1] - 1)*sc[6*npt+i] + (TB[1][2] - 1)*sc[18*npt+i] + (TB[1][3] - 1)*sc[12*npt+i];
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
        wdot[2*npt+i] -= qdot;
        wdot[5*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;

        /*reaction 3: CH3 + H (+M) => CH4 (+M) */
        phi_f = sc[2*npt+i]*sc[11*npt+i];
        alpha = mixture[i] + (TB[2][0] - 1)*sc[4*npt+i] + (TB[2][1] - 1)*sc[6*npt+i] + (TB[2][2] - 1)*sc[18*npt+i] + (TB[2][3] - 1)*sc[12*npt+i];
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
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[11*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;

        /*reaction 4: 2.000000 CH3 (+M) <=> C2H6 (+M) */
        phi_f = pow(sc[11*npt+i], 2.000000);
        alpha = mixture[i] + (TB[3][0] - 1)*sc[4*npt+i] + (TB[3][1] - 1)*sc[6*npt+i] + (TB[3][2] - 1)*sc[18*npt+i] + (TB[3][3] - 1)*sc[12*npt+i];
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
        phi_r = sc[17*npt+i];
        Kc = Kc_s[3*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[11*npt+i] -= 2.000000 * qdot;
        wdot[17*npt+i] += qdot;

        /*reaction 5: CO + O (+M) => CO2 (+M) */
        phi_f = sc[1*npt+i]*sc[12*npt+i];
        alpha = mixture[i] + (TB[4][0] - 1)*sc[4*npt+i] + (TB[4][1] - 1)*sc[6*npt+i] + (TB[4][2] - 1)*sc[18*npt+i] + (TB[4][3] - 1)*sc[12*npt+i];
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
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[12*npt+i] -= qdot;
        wdot[18*npt+i] += qdot;

        /*reaction 6: CH3O (+M) => CH2O + H (+M) */
        phi_f = sc[13*npt+i];
        alpha = mixture[i];
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
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[10*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;

        /*reaction 7: C2H3 (+M) => H + C2H2 (+M) */
        phi_f = sc[24*npt+i];
        alpha = mixture[i] + (TB[6][0] - 1)*sc[4*npt+i] + (TB[6][1] - 1)*sc[6*npt+i] + (TB[6][2] - 1)*sc[18*npt+i] + (TB[6][3] - 1)*sc[12*npt+i];
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
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[22*npt+i] += qdot;
        wdot[24*npt+i] -= qdot;

        /*reaction 8: H + C2H4 (+M) <=> C2H5 (+M) */
        phi_f = sc[2*npt+i]*sc[16*npt+i];
        alpha = mixture[i];
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
        phi_r = sc[14*npt+i];
        Kc = Kc_s[7*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;

        /*reaction 9: C2H4 (+M) => C2H2 + H2 (+M) */
        phi_f = sc[16*npt+i];
        alpha = mixture[i];
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
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[4*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;
        wdot[22*npt+i] += qdot;

        /*reaction 10: O + H + M => OH + M */
        phi_f = sc[1*npt+i]*sc[2*npt+i];
        alpha = mixture[i] + (TB[9][0] - 1)*sc[4*npt+i] + (TB[9][1] - 1)*sc[6*npt+i] + (TB[9][2] - 1)*sc[18*npt+i] + (TB[9][3] - 1)*sc[12*npt+i];
        k_f = alpha * k_f_s[9*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;

        /*reaction 11: 2.000000 O + M => O2 + M */
        phi_f = pow(sc[1*npt+i], 2.000000);
        alpha = mixture[i] + (TB[10][0] - 1)*sc[4*npt+i] + (TB[10][1] - 1)*sc[6*npt+i] + (TB[10][2] - 1)*sc[18*npt+i] + (TB[10][3] - 1)*sc[12*npt+i];
        k_f = alpha * k_f_s[10*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= 2.000000 * qdot;
        wdot[5*npt+i] += qdot;

        /*reaction 12: H + OH + M => H2O + M */
        phi_f = sc[2*npt+i]*sc[3*npt+i];
        alpha = mixture[i] + (TB[11][0] - 1)*sc[4*npt+i] + (TB[11][1] - 1)*sc[6*npt+i] + (TB[11][2] - 1)*sc[18*npt+i] + (TB[11][3] - 1)*sc[12*npt+i];
        k_f = alpha * k_f_s[11*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;

        /*reaction 13: HCO + M => H + CO + M */
        phi_f = sc[19*npt+i];
        alpha = mixture[i] + (TB[12][0] - 1)*sc[4*npt+i] + (TB[12][1] - 1)*sc[6*npt+i] + (TB[12][2] - 1)*sc[18*npt+i] + (TB[12][3] - 1)*sc[12*npt+i];
        k_f = alpha * k_f_s[12*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[12*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 14: H + CO + M => HCO + M */
        phi_f = sc[2*npt+i]*sc[12*npt+i];
        alpha = mixture[i] + (TB[13][0] - 1)*sc[4*npt+i] + (TB[13][1] - 1)*sc[6*npt+i] + (TB[13][2] - 1)*sc[18*npt+i] + (TB[13][3] - 1)*sc[12*npt+i];
        k_f = alpha * k_f_s[13*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[12*npt+i] -= qdot;
        wdot[19*npt+i] += qdot;

        /*reaction 15: CH3O2 + M => CH3 + O2 + M */
        phi_f = sc[20*npt+i];
        alpha = mixture[i];
        k_f = alpha * k_f_s[14*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[5*npt+i] += qdot;
        wdot[11*npt+i] += qdot;
        wdot[20*npt+i] -= qdot;

        /*reaction 16: CH3 + O2 + M => CH3O2 + M */
        phi_f = sc[5*npt+i]*sc[11*npt+i];
        alpha = mixture[i];
        k_f = alpha * k_f_s[15*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[11*npt+i] -= qdot;
        wdot[20*npt+i] += qdot;

        /*reaction 17: O + H2 => H + OH */
        phi_f = sc[1*npt+i]*sc[4*npt+i];
        k_f = k_f_s[16*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;

        /*reaction 18: H + OH => O + H2 */
        phi_f = sc[2*npt+i]*sc[3*npt+i];
        k_f = k_f_s[17*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;

        /*reaction 19: OH + H2 => H + H2O */
        phi_f = sc[3*npt+i]*sc[4*npt+i];
        k_f = k_f_s[18*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;

        /*reaction 20: H + H2O => OH + H2 */
        phi_f = sc[2*npt+i]*sc[6*npt+i];
        k_f = k_f_s[19*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[4*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 21: O + H2O => 2.000000 OH */
        phi_f = sc[1*npt+i]*sc[6*npt+i];
        k_f = k_f_s[20*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] += 2.000000 * qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 22: 2.000000 OH => O + H2O */
        phi_f = pow(sc[3*npt+i], 2.000000);
        k_f = k_f_s[21*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[3*npt+i] -= 2.000000 * qdot;
        wdot[6*npt+i] += qdot;

        /*reaction 23: H + O2 => O + OH */
        phi_f = sc[2*npt+i]*sc[5*npt+i];
        k_f = k_f_s[22*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[2*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;

        /*reaction 24: O + OH => H + O2 */
        phi_f = sc[1*npt+i]*sc[3*npt+i];
        k_f = k_f_s[23*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;

        /*reaction 25: HO2 + H => 2.000000 OH */
        phi_f = sc[2*npt+i]*sc[8*npt+i];
        k_f = k_f_s[24*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[3*npt+i] += 2.000000 * qdot;
        wdot[8*npt+i] -= qdot;

        /*reaction 26: 2.000000 HO2 => H2O2 + O2 */
        phi_f = pow(sc[8*npt+i], 2.000000);
        k_f = k_f_s[25*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[5*npt+i] += qdot;
        wdot[7*npt+i] += qdot;
        wdot[8*npt+i] -= 2.000000 * qdot;

        /*reaction 27: HO2 + H => H2 + O2 */
        phi_f = sc[2*npt+i]*sc[8*npt+i];
        k_f = k_f_s[26*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[5*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;

        /*reaction 28: HO2 + OH => H2O + O2 */
        phi_f = sc[3*npt+i]*sc[8*npt+i];
        k_f = k_f_s[27*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;

        /*reaction 29: H2O + O2 => HO2 + OH */
        phi_f = sc[5*npt+i]*sc[6*npt+i];
        k_f = k_f_s[28*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[6*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;

        /*reaction 30: 2.000000 HO2 => H2O2 + O2 */
        phi_f = pow(sc[8*npt+i], 2.000000);
        k_f = k_f_s[29*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[5*npt+i] += qdot;
        wdot[7*npt+i] += qdot;
        wdot[8*npt+i] -= 2.000000 * qdot;

        /*reaction 31: HO2 + O => OH + O2 */
        phi_f = sc[1*npt+i]*sc[8*npt+i];
        k_f = k_f_s[30*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[5*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;

        /*reaction 32: OH + O2 => HO2 + O */
        phi_f = sc[3*npt+i]*sc[5*npt+i];
        k_f = k_f_s[31*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[5*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;

        /*reaction 33: H2O2 + OH => H2O + HO2 */
        phi_f = sc[3*npt+i]*sc[7*npt+i];
        k_f = k_f_s[32*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;

        /*reaction 34: H2O2 + H => H2O + OH */
        phi_f = sc[2*npt+i]*sc[7*npt+i];
        k_f = k_f_s[33*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 35: H2O2 + OH => H2O + HO2 */
        phi_f = sc[3*npt+i]*sc[7*npt+i];
        k_f = k_f_s[34*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;

        /*reaction 36: H2O + HO2 => H2O2 + OH */
        phi_f = sc[6*npt+i]*sc[8*npt+i];
        k_f = k_f_s[35*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;

        /*reaction 37: H2O2 + O => OH + HO2 */
        phi_f = sc[1*npt+i]*sc[7*npt+i];
        k_f = k_f_s[36*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;

        /*reaction 38: H2O2 + H => H2 + HO2 */
        phi_f = sc[2*npt+i]*sc[7*npt+i];
        k_f = k_f_s[37*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;

        /*reaction 39: H2 + HO2 => H2O2 + H */
        phi_f = sc[4*npt+i]*sc[8*npt+i];
        k_f = k_f_s[38*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;

        /*reaction 40: CH2GSG + OH => CH2O + H */
        phi_f = sc[3*npt+i]*sc[9*npt+i];
        k_f = k_f_s[39*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[9*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;

        /*reaction 41: CH2GSG + H2 => CH3 + H */
        phi_f = sc[4*npt+i]*sc[9*npt+i];
        k_f = k_f_s[40*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[9*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;

        /*reaction 42: CH3 + H => CH2GSG + H2 */
        phi_f = sc[2*npt+i]*sc[11*npt+i];
        k_f = k_f_s[41*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[9*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;

        /*reaction 43: CH2GSG + O2 => CO + OH + H */
        phi_f = sc[5*npt+i]*sc[9*npt+i];
        k_f = k_f_s[42*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[9*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;

        /*reaction 44: CH2GSG + O => CO + 2.000000 H */
        phi_f = sc[1*npt+i]*sc[9*npt+i];
        k_f = k_f_s[43*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += 2.000000 * qdot;
        wdot[9*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;

        /*reaction 45: CH3 + HO2 => CH3O + OH */
        phi_f = sc[8*npt+i]*sc[11*npt+i];
        k_f = k_f_s[44*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[11*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;

        /*reaction 46: CH3 + O2 => CH2O + OH */
        phi_f = sc[5*npt+i]*sc[11*npt+i];
        k_f = k_f_s[45*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;

        /*reaction 47: CH3 + O2 => CH3O + O */
        phi_f = sc[5*npt+i]*sc[11*npt+i];
        k_f = k_f_s[46*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[11*npt+i] -= qdot;
        wdot[13*npt+i] += qdot;

        /*reaction 48: CH3O + O => CH3 + O2 */
        phi_f = sc[1*npt+i]*sc[13*npt+i];
        k_f = k_f_s[47*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[11*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;

        /*reaction 49: 2.000000 CH3 <=> H + C2H5 */
        phi_f = pow(sc[11*npt+i], 2.000000);
        k_f = k_f_s[48*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[14*npt+i];
        Kc = Kc_s[48*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[11*npt+i] -= 2.000000 * qdot;
        wdot[14*npt+i] += qdot;

        /*reaction 50: CH3 + HO2 => CH4 + O2 */
        phi_f = sc[8*npt+i]*sc[11*npt+i];
        k_f = k_f_s[49*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[5*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[11*npt+i] -= qdot;
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

        /*reaction 51: CH3 + O => CH2O + H */
        phi_f = sc[1*npt+i]*sc[11*npt+i];
        k_f = k_f_s[50*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[10*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;

        /*reaction 52: CH3 + OH => CH2O + H2 */
        phi_f = sc[3*npt+i]*sc[11*npt+i];
        k_f = k_f_s[51*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[10*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;

        /*reaction 53: CH3 + OH => CH2GSG + H2O */
        phi_f = sc[3*npt+i]*sc[11*npt+i];
        k_f = k_f_s[52*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[9*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;

        /*reaction 54: CH2GSG + H2O => CH3 + OH */
        phi_f = sc[6*npt+i]*sc[9*npt+i];
        k_f = k_f_s[53*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[9*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;

        /*reaction 55: CH3 + H2O2 => CH4 + HO2 */
        phi_f = sc[7*npt+i]*sc[11*npt+i];
        k_f = k_f_s[54*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[7*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;

        /*reaction 56: CH2GSG + CH3 => C2H4 + H */
        phi_f = sc[9*npt+i]*sc[11*npt+i];
        k_f = k_f_s[55*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[9*npt+i] -= qdot;
        wdot[11*npt+i] -= qdot;
        wdot[16*npt+i] += qdot;

        /*reaction 57: CH2GSG + CH4 => 2.000000 CH3 */
        phi_f = sc[9*npt+i]*sc[15*npt+i];
        k_f = k_f_s[56*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[9*npt+i] -= qdot;
        wdot[11*npt+i] += 2.000000 * qdot;
        wdot[15*npt+i] -= qdot;

        /*reaction 58: 2.000000 CH3 => CH2GSG + CH4 */
        phi_f = pow(sc[11*npt+i], 2.000000);
        k_f = k_f_s[57*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[9*npt+i] += qdot;
        wdot[11*npt+i] -= 2.000000 * qdot;
        wdot[15*npt+i] += qdot;

        /*reaction 59: CH4 + O => CH3 + OH */
        phi_f = sc[1*npt+i]*sc[15*npt+i];
        k_f = k_f_s[58*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[11*npt+i] += qdot;
        wdot[15*npt+i] -= qdot;

        /*reaction 60: CH3 + OH => CH4 + O */
        phi_f = sc[3*npt+i]*sc[11*npt+i];
        k_f = k_f_s[59*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[11*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;

        /*reaction 61: CH4 + H => CH3 + H2 */
        phi_f = sc[2*npt+i]*sc[15*npt+i];
        k_f = k_f_s[60*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[11*npt+i] += qdot;
        wdot[15*npt+i] -= qdot;

        /*reaction 62: CH3 + H2 => CH4 + H */
        phi_f = sc[4*npt+i]*sc[11*npt+i];
        k_f = k_f_s[61*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[11*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;

        /*reaction 63: CH4 + OH => CH3 + H2O */
        phi_f = sc[3*npt+i]*sc[15*npt+i];
        k_f = k_f_s[62*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[11*npt+i] += qdot;
        wdot[15*npt+i] -= qdot;

        /*reaction 64: CH3 + H2O => CH4 + OH */
        phi_f = sc[6*npt+i]*sc[11*npt+i];
        k_f = k_f_s[63*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[11*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;

        /*reaction 65: CO + HO2 => CO2 + OH */
        phi_f = sc[8*npt+i]*sc[12*npt+i];
        k_f = k_f_s[64*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[12*npt+i] -= qdot;
        wdot[18*npt+i] += qdot;

        /*reaction 66: CO + O2 => CO2 + O */
        phi_f = sc[5*npt+i]*sc[12*npt+i];
        k_f = k_f_s[65*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[12*npt+i] -= qdot;
        wdot[18*npt+i] += qdot;

        /*reaction 67: CO2 + O => CO + O2 */
        phi_f = sc[1*npt+i]*sc[18*npt+i];
        k_f = k_f_s[66*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[12*npt+i] += qdot;
        wdot[18*npt+i] -= qdot;

        /*reaction 68: CO + OH => CO2 + H */
        phi_f = sc[3*npt+i]*sc[12*npt+i];
        k_f = k_f_s[67*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[12*npt+i] -= qdot;
        wdot[18*npt+i] += qdot;

        /*reaction 69: CO2 + H => CO + OH */
        phi_f = sc[2*npt+i]*sc[18*npt+i];
        k_f = k_f_s[68*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[12*npt+i] += qdot;
        wdot[18*npt+i] -= qdot;

        /*reaction 70: HCO + CH3 => CH4 + CO */
        phi_f = sc[11*npt+i]*sc[19*npt+i];
        k_f = k_f_s[69*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[11*npt+i] -= qdot;
        wdot[12*npt+i] += qdot;
        wdot[15*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 71: HCO + H => CO + H2 */
        phi_f = sc[2*npt+i]*sc[19*npt+i];
        k_f = k_f_s[70*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[12*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 72: HCO + O2 => CO + HO2 */
        phi_f = sc[5*npt+i]*sc[19*npt+i];
        k_f = k_f_s[71*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[12*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 73: HCO + O => CO + OH */
        phi_f = sc[1*npt+i]*sc[19*npt+i];
        k_f = k_f_s[72*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[12*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 74: HCO + O => CO2 + H */
        phi_f = sc[1*npt+i]*sc[19*npt+i];
        k_f = k_f_s[73*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[18*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 75: HCO + OH => CO + H2O */
        phi_f = sc[3*npt+i]*sc[19*npt+i];
        k_f = k_f_s[74*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[12*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;

        /*reaction 76: CH2O + H => HCO + H2 */
        phi_f = sc[2*npt+i]*sc[10*npt+i];
        k_f = k_f_s[75*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;
        wdot[19*npt+i] += qdot;

        /*reaction 77: CH2O + O2 => HCO + HO2 */
        phi_f = sc[5*npt+i]*sc[10*npt+i];
        k_f = k_f_s[76*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;
        wdot[19*npt+i] += qdot;

        /*reaction 78: CH2O + OH => HCO + H2O */
        phi_f = sc[3*npt+i]*sc[10*npt+i];
        k_f = k_f_s[77*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;
        wdot[19*npt+i] += qdot;

        /*reaction 79: CH2O + HO2 => HCO + H2O2 */
        phi_f = sc[8*npt+i]*sc[10*npt+i];
        k_f = k_f_s[78*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[7*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[10*npt+i] -= qdot;
        wdot[19*npt+i] += qdot;

        /*reaction 80: CH2O + O => HCO + OH */
        phi_f = sc[1*npt+i]*sc[10*npt+i];
        k_f = k_f_s[79*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;
        wdot[19*npt+i] += qdot;

        /*reaction 81: CH2O + CH3 => HCO + CH4 */
        phi_f = sc[10*npt+i]*sc[11*npt+i];
        k_f = k_f_s[80*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[10*npt+i] -= qdot;
        wdot[11*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;
        wdot[19*npt+i] += qdot;

        /*reaction 82: CH3O + O2 => CH2O + HO2 */
        phi_f = sc[5*npt+i]*sc[13*npt+i];
        k_f = k_f_s[81*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[10*npt+i] += qdot;
        wdot[13*npt+i] -= qdot;

        /*reaction 83: CH2GSG + CO2 => CH2O + CO */
        phi_f = sc[9*npt+i]*sc[18*npt+i];
        k_f = k_f_s[82*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[9*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;
        wdot[12*npt+i] += qdot;
        wdot[18*npt+i] -= qdot;

        /*reaction 84: CH3O2 + CH3 => 2.000000 CH3O */
        phi_f = sc[11*npt+i]*sc[20*npt+i];
        k_f = k_f_s[83*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[11*npt+i] -= qdot;
        wdot[13*npt+i] += 2.000000 * qdot;
        wdot[20*npt+i] -= qdot;

        /*reaction 85: CH3O2 + HO2 => CH3O2H + O2 */
        phi_f = sc[8*npt+i]*sc[20*npt+i];
        k_f = k_f_s[84*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[5*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[20*npt+i] -= qdot;
        wdot[21*npt+i] += qdot;

        /*reaction 86: 2.000000 CH3O2 => O2 + 2.000000 CH3O */
        phi_f = pow(sc[20*npt+i], 2.000000);
        k_f = k_f_s[85*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[5*npt+i] += qdot;
        wdot[13*npt+i] += 2.000000 * qdot;
        wdot[20*npt+i] -= 2.000000 * qdot;

        /*reaction 87: CH3O2 + CH2O => CH3O2H + HCO */
        phi_f = sc[10*npt+i]*sc[20*npt+i];
        k_f = k_f_s[86*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[10*npt+i] -= qdot;
        wdot[19*npt+i] += qdot;
        wdot[20*npt+i] -= qdot;
        wdot[21*npt+i] += qdot;

        /*reaction 88: CH3O2H => CH3O + OH */
        phi_f = sc[21*npt+i];
        k_f = k_f_s[87*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[13*npt+i] += qdot;
        wdot[21*npt+i] -= qdot;

        /*reaction 89: CH3O + OH => CH3O2H */
        phi_f = sc[3*npt+i]*sc[13*npt+i];
        k_f = k_f_s[88*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[13*npt+i] -= qdot;
        wdot[21*npt+i] += qdot;

        /*reaction 90: C2H2 + O2 => HCCO + OH */
        phi_f = sc[5*npt+i]*sc[22*npt+i];
        k_f = k_f_s[89*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[22*npt+i] -= qdot;
        wdot[23*npt+i] += qdot;

        /*reaction 91: C2H2 + O => HCCO + H */
        phi_f = sc[1*npt+i]*sc[22*npt+i];
        k_f = k_f_s[90*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[22*npt+i] -= qdot;
        wdot[23*npt+i] += qdot;

        /*reaction 92: C2H3 + O2 => CH2CHO + O */
        phi_f = sc[5*npt+i]*sc[24*npt+i];
        k_f = k_f_s[91*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[1*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[24*npt+i] -= qdot;
        wdot[25*npt+i] += qdot;

        /*reaction 93: C2H3 + H => C2H2 + H2 */
        phi_f = sc[2*npt+i]*sc[24*npt+i];
        k_f = k_f_s[92*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[22*npt+i] += qdot;
        wdot[24*npt+i] -= qdot;

        /*reaction 94: C2H3 + O2 => CH2O + HCO */
        phi_f = sc[5*npt+i]*sc[24*npt+i];
        k_f = k_f_s[93*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;
        wdot[19*npt+i] += qdot;
        wdot[24*npt+i] -= qdot;

        /*reaction 95: C2H3 + O2 => C2H2 + HO2 */
        phi_f = sc[5*npt+i]*sc[24*npt+i];
        k_f = k_f_s[94*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[22*npt+i] += qdot;
        wdot[24*npt+i] -= qdot;

        /*reaction 96: C2H4 + O => CH3 + HCO */
        phi_f = sc[1*npt+i]*sc[16*npt+i];
        k_f = k_f_s[95*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;
        wdot[19*npt+i] += qdot;

        /*reaction 97: C2H4 + OH => C2H3 + H2O */
        phi_f = sc[3*npt+i]*sc[16*npt+i];
        k_f = k_f_s[96*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;
        wdot[24*npt+i] += qdot;

        /*reaction 98: C2H3 + H2O => C2H4 + OH */
        phi_f = sc[6*npt+i]*sc[24*npt+i];
        k_f = k_f_s[97*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[16*npt+i] += qdot;
        wdot[24*npt+i] -= qdot;

        /*reaction 99: C2H4 + H => C2H3 + H2 */
        phi_f = sc[2*npt+i]*sc[16*npt+i];
        k_f = k_f_s[98*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;
        wdot[24*npt+i] += qdot;

        /*reaction 100: C2H3 + H2 => C2H4 + H */
        phi_f = sc[4*npt+i]*sc[24*npt+i];
        k_f = k_f_s[99*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[16*npt+i] += qdot;
        wdot[24*npt+i] -= qdot;
    }
}

void vcomp_wdot_101_150(int npt, amrex::Real *  wdot, amrex::Real *  mixture, amrex::Real *  sc,
		amrex::Real *  k_f_s, amrex::Real *  Kc_s,
		amrex::Real *  tc, amrex::Real *  invT, amrex::Real *  T)
{
    for (int i=0; i<npt; i++) {
        amrex::Real qdot, q_f, q_r, phi_f, phi_r, k_f, k_r, Kc;
        amrex::Real alpha;

        /*reaction 101: C2H4 + CH3O2 => C2H3 + CH3O2H */
        phi_f = sc[16*npt+i]*sc[20*npt+i];
        k_f = k_f_s[100*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[16*npt+i] -= qdot;
        wdot[20*npt+i] -= qdot;
        wdot[21*npt+i] += qdot;
        wdot[24*npt+i] += qdot;

        /*reaction 102: C2H4 + CH3 => C2H3 + CH4 */
        phi_f = sc[11*npt+i]*sc[16*npt+i];
        k_f = k_f_s[101*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[11*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;
        wdot[24*npt+i] += qdot;

        /*reaction 103: C2H3 + CH4 => C2H4 + CH3 */
        phi_f = sc[15*npt+i]*sc[24*npt+i];
        k_f = k_f_s[102*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[11*npt+i] += qdot;
        wdot[15*npt+i] -= qdot;
        wdot[16*npt+i] += qdot;
        wdot[24*npt+i] -= qdot;

        /*reaction 104: C2H4 + O => CH2CHO + H */
        phi_f = sc[1*npt+i]*sc[16*npt+i];
        k_f = k_f_s[103*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;
        wdot[25*npt+i] += qdot;

        /*reaction 105: CH3 + C2H5 => CH4 + C2H4 */
        phi_f = sc[11*npt+i]*sc[14*npt+i];
        k_f = k_f_s[104*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[11*npt+i] -= qdot;
        wdot[14*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;
        wdot[16*npt+i] += qdot;

        /*reaction 106: C2H5 + O2 => C2H4 + HO2 */
        phi_f = sc[5*npt+i]*sc[14*npt+i];
        k_f = k_f_s[105*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[14*npt+i] -= qdot;
        wdot[16*npt+i] += qdot;

        /*reaction 107: C2H4 + HO2 => C2H5 + O2 */
        phi_f = sc[8*npt+i]*sc[16*npt+i];
        k_f = k_f_s[106*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[5*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;
        wdot[16*npt+i] -= qdot;

        /*reaction 108: C2H5 + HO2 => C2H6 + O2 */
        phi_f = sc[8*npt+i]*sc[14*npt+i];
        k_f = k_f_s[107*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[5*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[14*npt+i] -= qdot;
        wdot[17*npt+i] += qdot;

        /*reaction 109: H + C2H5 => C2H6 */
        phi_f = sc[2*npt+i]*sc[14*npt+i];
        k_f = k_f_s[108*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[14*npt+i] -= qdot;
        wdot[17*npt+i] += qdot;

        /*reaction 110: C2H6 + H => C2H5 + H2 */
        phi_f = sc[2*npt+i]*sc[17*npt+i];
        k_f = k_f_s[109*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[14*npt+i] += qdot;
        wdot[17*npt+i] -= qdot;

        /*reaction 111: C2H5 + H2 => C2H6 + H */
        phi_f = sc[4*npt+i]*sc[14*npt+i];
        k_f = k_f_s[110*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[14*npt+i] -= qdot;
        wdot[17*npt+i] += qdot;

        /*reaction 112: C2H6 + OH => C2H5 + H2O */
        phi_f = sc[3*npt+i]*sc[17*npt+i];
        k_f = k_f_s[111*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[14*npt+i] += qdot;
        wdot[17*npt+i] -= qdot;

        /*reaction 113: CH2GSG + C2H6 => CH3 + C2H5 */
        phi_f = sc[9*npt+i]*sc[17*npt+i];
        k_f = k_f_s[112*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[9*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;
        wdot[14*npt+i] += qdot;
        wdot[17*npt+i] -= qdot;

        /*reaction 114: C2H6 + O => C2H5 + OH */
        phi_f = sc[1*npt+i]*sc[17*npt+i];
        k_f = k_f_s[113*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[14*npt+i] += qdot;
        wdot[17*npt+i] -= qdot;

        /*reaction 115: C2H6 + CH3 => C2H5 + CH4 */
        phi_f = sc[11*npt+i]*sc[17*npt+i];
        k_f = k_f_s[114*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[11*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;
        wdot[15*npt+i] += qdot;
        wdot[17*npt+i] -= qdot;

        /*reaction 116: HCCO + O => H + 2.000000 CO */
        phi_f = sc[1*npt+i]*sc[23*npt+i];
        k_f = k_f_s[115*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[12*npt+i] += 2.000000 * qdot;
        wdot[23*npt+i] -= qdot;

        /*reaction 117: HCCO + O2 => CO2 + HCO */
        phi_f = sc[5*npt+i]*sc[23*npt+i];
        k_f = k_f_s[116*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[18*npt+i] += qdot;
        wdot[19*npt+i] += qdot;
        wdot[23*npt+i] -= qdot;

        /*reaction 118: HCCO + OH => 2.000000 HCO */
        phi_f = sc[3*npt+i]*sc[23*npt+i];
        k_f = k_f_s[117*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[19*npt+i] += 2.000000 * qdot;
        wdot[23*npt+i] -= qdot;

        /*reaction 119: HCCO + H => CH2GSG + CO */
        phi_f = sc[2*npt+i]*sc[23*npt+i];
        k_f = k_f_s[118*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;
        wdot[12*npt+i] += qdot;
        wdot[23*npt+i] -= qdot;

        /*reaction 120: CH2GSG + CO => HCCO + H */
        phi_f = sc[9*npt+i]*sc[12*npt+i];
        k_f = k_f_s[119*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[9*npt+i] -= qdot;
        wdot[12*npt+i] -= qdot;
        wdot[23*npt+i] += qdot;

        /*reaction 121: CH2CHO + O2 => CH2O + CO + OH */
        phi_f = sc[5*npt+i]*sc[25*npt+i];
        k_f = k_f_s[120*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;
        wdot[12*npt+i] += qdot;
        wdot[25*npt+i] -= qdot;

        /*reaction 122: C3H5XA + CH2O => C3H6 + HCO */
        phi_f = sc[10*npt+i]*sc[26*npt+i];
        k_f = k_f_s[121*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[10*npt+i] -= qdot;
        wdot[19*npt+i] += qdot;
        wdot[26*npt+i] -= qdot;
        wdot[27*npt+i] += qdot;

        /*reaction 123: C3H6 + HCO => C3H5XA + CH2O */
        phi_f = sc[19*npt+i]*sc[27*npt+i];
        k_f = k_f_s[122*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[10*npt+i] += qdot;
        wdot[19*npt+i] -= qdot;
        wdot[26*npt+i] += qdot;
        wdot[27*npt+i] -= qdot;

        /*reaction 124: C3H5XA + HO2 => C3H5O + OH */
        phi_f = sc[8*npt+i]*sc[26*npt+i];
        k_f = k_f_s[123*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[26*npt+i] -= qdot;
        wdot[28*npt+i] += qdot;

        /*reaction 125: C3H5XA + CH3O2 => C3H5O + CH3O */
        phi_f = sc[20*npt+i]*sc[26*npt+i];
        k_f = k_f_s[124*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[13*npt+i] += qdot;
        wdot[20*npt+i] -= qdot;
        wdot[26*npt+i] -= qdot;
        wdot[28*npt+i] += qdot;

        /*reaction 126: C3H5XA + HO2 => C3H6 + O2 */
        phi_f = sc[8*npt+i]*sc[26*npt+i];
        k_f = k_f_s[125*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[5*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[26*npt+i] -= qdot;
        wdot[27*npt+i] += qdot;

        /*reaction 127: C3H5XA + O2 => CH2CHO + CH2O */
        phi_f = sc[5*npt+i]*sc[26*npt+i];
        k_f = k_f_s[126*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;
        wdot[25*npt+i] += qdot;
        wdot[26*npt+i] -= qdot;

        /*reaction 128: C3H5XA + O2 => C2H2 + CH2O + OH */
        phi_f = sc[5*npt+i]*sc[26*npt+i];
        k_f = k_f_s[127*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[3*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;
        wdot[22*npt+i] += qdot;
        wdot[26*npt+i] -= qdot;

        /*reaction 129: C3H5XA + H => C3H6 */
        phi_f = sc[2*npt+i]*sc[26*npt+i];
        k_f = k_f_s[128*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[26*npt+i] -= qdot;
        wdot[27*npt+i] += qdot;

        /*reaction 130: C3H5XA => C2H2 + CH3 */
        phi_f = sc[26*npt+i];
        k_f = k_f_s[129*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[11*npt+i] += qdot;
        wdot[22*npt+i] += qdot;
        wdot[26*npt+i] -= qdot;

        /*reaction 131: C3H6 => C2H3 + CH3 */
        phi_f = sc[27*npt+i];
        k_f = k_f_s[130*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[11*npt+i] += qdot;
        wdot[24*npt+i] += qdot;
        wdot[27*npt+i] -= qdot;

        /*reaction 132: C2H3 + CH3 => C3H6 */
        phi_f = sc[11*npt+i]*sc[24*npt+i];
        k_f = k_f_s[131*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[11*npt+i] -= qdot;
        wdot[24*npt+i] -= qdot;
        wdot[27*npt+i] += qdot;

        /*reaction 133: C3H6 + O => C2H5 + HCO */
        phi_f = sc[1*npt+i]*sc[27*npt+i];
        k_f = k_f_s[132*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[14*npt+i] += qdot;
        wdot[19*npt+i] += qdot;
        wdot[27*npt+i] -= qdot;

        /*reaction 134: C3H6 + H => C2H4 + CH3 */
        phi_f = sc[2*npt+i]*sc[27*npt+i];
        k_f = k_f_s[133*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[11*npt+i] += qdot;
        wdot[16*npt+i] += qdot;
        wdot[27*npt+i] -= qdot;

        /*reaction 135: C2H4 + CH3 => C3H6 + H */
        phi_f = sc[11*npt+i]*sc[16*npt+i];
        k_f = k_f_s[134*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[11*npt+i] -= qdot;
        wdot[16*npt+i] -= qdot;
        wdot[27*npt+i] += qdot;

        /*reaction 136: C3H6 + H => C3H5XA + H2 */
        phi_f = sc[2*npt+i]*sc[27*npt+i];
        k_f = k_f_s[135*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[26*npt+i] += qdot;
        wdot[27*npt+i] -= qdot;

        /*reaction 137: C3H5XA + H2 => C3H6 + H */
        phi_f = sc[4*npt+i]*sc[26*npt+i];
        k_f = k_f_s[136*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[26*npt+i] -= qdot;
        wdot[27*npt+i] += qdot;

        /*reaction 138: C3H6 + OH => C3H5XA + H2O */
        phi_f = sc[3*npt+i]*sc[27*npt+i];
        k_f = k_f_s[137*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[26*npt+i] += qdot;
        wdot[27*npt+i] -= qdot;

        /*reaction 139: C3H6 + O => C3H5XA + OH */
        phi_f = sc[1*npt+i]*sc[27*npt+i];
        k_f = k_f_s[138*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[26*npt+i] += qdot;
        wdot[27*npt+i] -= qdot;

        /*reaction 140: C3H6 + CH3 => C3H5XA + CH4 */
        phi_f = sc[11*npt+i]*sc[27*npt+i];
        k_f = k_f_s[139*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[11*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;
        wdot[26*npt+i] += qdot;
        wdot[27*npt+i] -= qdot;

        /*reaction 141: IXC3H7 + O2 => C3H6 + HO2 */
        phi_f = sc[5*npt+i]*sc[29*npt+i];
        k_f = k_f_s[140*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[27*npt+i] += qdot;
        wdot[29*npt+i] -= qdot;

        /*reaction 142: IXC3H7 => H + C3H6 */
        phi_f = sc[29*npt+i];
        k_f = k_f_s[141*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[27*npt+i] += qdot;
        wdot[29*npt+i] -= qdot;

        /*reaction 143: H + C3H6 => IXC3H7 */
        phi_f = sc[2*npt+i]*sc[27*npt+i];
        k_f = k_f_s[142*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[27*npt+i] -= qdot;
        wdot[29*npt+i] += qdot;

        /*reaction 144: NXC3H7 => CH3 + C2H4 */
        phi_f = sc[30*npt+i];
        k_f = k_f_s[143*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[11*npt+i] += qdot;
        wdot[16*npt+i] += qdot;
        wdot[30*npt+i] -= qdot;

        /*reaction 145: CH3 + C2H4 => NXC3H7 */
        phi_f = sc[11*npt+i]*sc[16*npt+i];
        k_f = k_f_s[144*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[11*npt+i] -= qdot;
        wdot[16*npt+i] -= qdot;
        wdot[30*npt+i] += qdot;

        /*reaction 146: NXC3H7 + HO2 => C3H8 + O2 */
        phi_f = sc[8*npt+i]*sc[30*npt+i];
        k_f = k_f_s[145*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[5*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[30*npt+i] -= qdot;
        wdot[31*npt+i] += qdot;

        /*reaction 147: NXC3H7 + O2 => C3H6 + HO2 */
        phi_f = sc[5*npt+i]*sc[30*npt+i];
        k_f = k_f_s[146*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[27*npt+i] += qdot;
        wdot[30*npt+i] -= qdot;

        /*reaction 148: NXC3H7 => H + C3H6 */
        phi_f = sc[30*npt+i];
        k_f = k_f_s[147*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[27*npt+i] += qdot;
        wdot[30*npt+i] -= qdot;

        /*reaction 149: H + C3H6 => NXC3H7 */
        phi_f = sc[2*npt+i]*sc[27*npt+i];
        k_f = k_f_s[148*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[27*npt+i] -= qdot;
        wdot[30*npt+i] += qdot;

        /*reaction 150: C3H8 + OH => NXC3H7 + H2O */
        phi_f = sc[3*npt+i]*sc[31*npt+i];
        k_f = k_f_s[149*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[30*npt+i] += qdot;
        wdot[31*npt+i] -= qdot;
    }
}

void vcomp_wdot_151_167(int npt, amrex::Real *  wdot, amrex::Real *  mixture, amrex::Real *  sc,
		amrex::Real *  k_f_s, amrex::Real *  Kc_s,
		amrex::Real *  tc, amrex::Real *  invT, amrex::Real *  T)
{
    for (int i=0; i<npt; i++) {
        amrex::Real qdot, q_f, q_r, phi_f, phi_r, k_f, k_r, Kc;
        amrex::Real alpha;

        /*reaction 151: C3H8 + HO2 => NXC3H7 + H2O2 */
        phi_f = sc[8*npt+i]*sc[31*npt+i];
        k_f = k_f_s[150*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[7*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[30*npt+i] += qdot;
        wdot[31*npt+i] -= qdot;

        /*reaction 152: H + C3H8 <=> H2 + NXC3H7 */
        phi_f = sc[2*npt+i]*sc[31*npt+i];
        k_f = k_f_s[151*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[30*npt+i];
        Kc = Kc_s[151*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[30*npt+i] += qdot;
        wdot[31*npt+i] -= qdot;

        /*reaction 153: C3H8 + OH => IXC3H7 + H2O */
        phi_f = sc[3*npt+i]*sc[31*npt+i];
        k_f = k_f_s[152*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;
        wdot[29*npt+i] += qdot;
        wdot[31*npt+i] -= qdot;

        /*reaction 154: CH3 + C3H8 => CH4 + IXC3H7 */
        phi_f = sc[11*npt+i]*sc[31*npt+i];
        k_f = k_f_s[153*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[11*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;
        wdot[29*npt+i] += qdot;
        wdot[31*npt+i] -= qdot;

        /*reaction 155: CH3 + C3H8 => CH4 + NXC3H7 */
        phi_f = sc[11*npt+i]*sc[31*npt+i];
        k_f = k_f_s[154*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[11*npt+i] -= qdot;
        wdot[15*npt+i] += qdot;
        wdot[30*npt+i] += qdot;
        wdot[31*npt+i] -= qdot;

        /*reaction 156: C3H8 + O => IXC3H7 + OH */
        phi_f = sc[1*npt+i]*sc[31*npt+i];
        k_f = k_f_s[155*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[29*npt+i] += qdot;
        wdot[31*npt+i] -= qdot;

        /*reaction 157: C3H8 + HO2 => IXC3H7 + H2O2 */
        phi_f = sc[8*npt+i]*sc[31*npt+i];
        k_f = k_f_s[156*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[7*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[29*npt+i] += qdot;
        wdot[31*npt+i] -= qdot;

        /*reaction 158: C3H8 + O => NXC3H7 + OH */
        phi_f = sc[1*npt+i]*sc[31*npt+i];
        k_f = k_f_s[157*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[30*npt+i] += qdot;
        wdot[31*npt+i] -= qdot;

        /*reaction 159: C3H8 + O2 => IXC3H7 + HO2 */
        phi_f = sc[5*npt+i]*sc[31*npt+i];
        k_f = k_f_s[158*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[8*npt+i] += qdot;
        wdot[29*npt+i] += qdot;
        wdot[31*npt+i] -= qdot;

        /*reaction 160: IXC3H7 + HO2 => C3H8 + O2 */
        phi_f = sc[8*npt+i]*sc[29*npt+i];
        k_f = k_f_s[159*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[5*npt+i] += qdot;
        wdot[8*npt+i] -= qdot;
        wdot[29*npt+i] -= qdot;
        wdot[31*npt+i] += qdot;

        /*reaction 161: H + C3H8 => H2 + IXC3H7 */
        phi_f = sc[2*npt+i]*sc[31*npt+i];
        k_f = k_f_s[160*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[29*npt+i] += qdot;
        wdot[31*npt+i] -= qdot;

        /*reaction 162: H2 + IXC3H7 => H + C3H8 */
        phi_f = sc[4*npt+i]*sc[29*npt+i];
        k_f = k_f_s[161*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[29*npt+i] -= qdot;
        wdot[31*npt+i] += qdot;

        /*reaction 163: C3H5O => C2H3 + CH2O */
        phi_f = sc[28*npt+i];
        k_f = k_f_s[162*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[10*npt+i] += qdot;
        wdot[24*npt+i] += qdot;
        wdot[28*npt+i] -= qdot;

        /*reaction 164: IXC3H7O2 => IXC3H7 + O2 */
        phi_f = sc[32*npt+i];
        k_f = k_f_s[163*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[5*npt+i] += qdot;
        wdot[29*npt+i] += qdot;
        wdot[32*npt+i] -= qdot;

        /*reaction 165: IXC3H7 + O2 => IXC3H7O2 */
        phi_f = sc[5*npt+i]*sc[29*npt+i];
        k_f = k_f_s[164*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[29*npt+i] -= qdot;
        wdot[32*npt+i] += qdot;

        /*reaction 166: NXC3H7O2 => NXC3H7 + O2 */
        phi_f = sc[33*npt+i];
        k_f = k_f_s[165*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[5*npt+i] += qdot;
        wdot[30*npt+i] += qdot;
        wdot[33*npt+i] -= qdot;

        /*reaction 167: NXC3H7 + O2 => NXC3H7O2 */
        phi_f = sc[5*npt+i]*sc[30*npt+i];
        k_f = k_f_s[166*npt+i];
        q_f = phi_f * k_f;
        q_r = 0.0;
        qdot = q_f - q_r;
        wdot[5*npt+i] -= qdot;
        wdot[30*npt+i] -= qdot;
        wdot[33*npt+i] += qdot;
    }
}


static amrex::Real T_save = -1;
#ifdef _OPENMP
#pragma omp threadprivate(T_save)
#endif

static amrex::Real k_f_save[167];
#ifdef _OPENMP
#pragma omp threadprivate(k_f_save)
#endif

static amrex::Real Kc_save[167];
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

    amrex::Real qdot, q_f[167], q_r[167];
    amrex::Real sc_qss[1];
    comp_qfqr(q_f, q_r, sc, sc_qss, tc, invT);

    for (int i = 0; i < 34; ++i) {
        wdot[i] = 0.0;
    }

    qdot = q_f[0]-q_r[0];
    wdot[3] -= 2.000000 * qdot;
    wdot[7] += qdot;

    qdot = q_f[1]-q_r[1];
    wdot[2] -= qdot;
    wdot[5] -= qdot;
    wdot[8] += qdot;

    qdot = q_f[2]-q_r[2];
    wdot[2] -= qdot;
    wdot[11] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[3]-q_r[3];
    wdot[11] -= 2.000000 * qdot;
    wdot[17] += qdot;

    qdot = q_f[4]-q_r[4];
    wdot[1] -= qdot;
    wdot[12] -= qdot;
    wdot[18] += qdot;

    qdot = q_f[5]-q_r[5];
    wdot[2] += qdot;
    wdot[10] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[6]-q_r[6];
    wdot[2] += qdot;
    wdot[22] += qdot;
    wdot[24] -= qdot;

    qdot = q_f[7]-q_r[7];
    wdot[2] -= qdot;
    wdot[14] += qdot;
    wdot[16] -= qdot;

    qdot = q_f[8]-q_r[8];
    wdot[4] += qdot;
    wdot[16] -= qdot;
    wdot[22] += qdot;

    qdot = q_f[9]-q_r[9];
    wdot[1] -= qdot;
    wdot[2] -= qdot;
    wdot[3] += qdot;

    qdot = q_f[10]-q_r[10];
    wdot[1] -= 2.000000 * qdot;
    wdot[5] += qdot;

    qdot = q_f[11]-q_r[11];
    wdot[2] -= qdot;
    wdot[3] -= qdot;
    wdot[6] += qdot;

    qdot = q_f[12]-q_r[12];
    wdot[2] += qdot;
    wdot[12] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[13]-q_r[13];
    wdot[2] -= qdot;
    wdot[12] -= qdot;
    wdot[19] += qdot;

    qdot = q_f[14]-q_r[14];
    wdot[5] += qdot;
    wdot[11] += qdot;
    wdot[20] -= qdot;

    qdot = q_f[15]-q_r[15];
    wdot[5] -= qdot;
    wdot[11] -= qdot;
    wdot[20] += qdot;

    qdot = q_f[16]-q_r[16];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[3] += qdot;
    wdot[4] -= qdot;

    qdot = q_f[17]-q_r[17];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[3] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[18]-q_r[18];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[4] -= qdot;
    wdot[6] += qdot;

    qdot = q_f[19]-q_r[19];
    wdot[2] -= qdot;
    wdot[3] += qdot;
    wdot[4] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[20]-q_r[20];
    wdot[1] -= qdot;
    wdot[3] += 2.000000 * qdot;
    wdot[6] -= qdot;

    qdot = q_f[21]-q_r[21];
    wdot[1] += qdot;
    wdot[3] -= 2.000000 * qdot;
    wdot[6] += qdot;

    qdot = q_f[22]-q_r[22];
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[3] += qdot;
    wdot[5] -= qdot;

    qdot = q_f[23]-q_r[23];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[24]-q_r[24];
    wdot[2] -= qdot;
    wdot[3] += 2.000000 * qdot;
    wdot[8] -= qdot;

    qdot = q_f[25]-q_r[25];
    wdot[5] += qdot;
    wdot[7] += qdot;
    wdot[8] -= 2.000000 * qdot;

    qdot = q_f[26]-q_r[26];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[5] += qdot;
    wdot[8] -= qdot;

    qdot = q_f[27]-q_r[27];
    wdot[3] -= qdot;
    wdot[5] += qdot;
    wdot[6] += qdot;
    wdot[8] -= qdot;

    qdot = q_f[28]-q_r[28];
    wdot[3] += qdot;
    wdot[5] -= qdot;
    wdot[6] -= qdot;
    wdot[8] += qdot;

    qdot = q_f[29]-q_r[29];
    wdot[5] += qdot;
    wdot[7] += qdot;
    wdot[8] -= 2.000000 * qdot;

    qdot = q_f[30]-q_r[30];
    wdot[1] -= qdot;
    wdot[3] += qdot;
    wdot[5] += qdot;
    wdot[8] -= qdot;

    qdot = q_f[31]-q_r[31];
    wdot[1] += qdot;
    wdot[3] -= qdot;
    wdot[5] -= qdot;
    wdot[8] += qdot;

    qdot = q_f[32]-q_r[32];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;
    wdot[8] += qdot;

    qdot = q_f[33]-q_r[33];
    wdot[2] -= qdot;
    wdot[3] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[34]-q_r[34];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;
    wdot[8] += qdot;

    qdot = q_f[35]-q_r[35];
    wdot[3] += qdot;
    wdot[6] -= qdot;
    wdot[7] += qdot;
    wdot[8] -= qdot;

    qdot = q_f[36]-q_r[36];
    wdot[1] -= qdot;
    wdot[3] += qdot;
    wdot[7] -= qdot;
    wdot[8] += qdot;

    qdot = q_f[37]-q_r[37];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[7] -= qdot;
    wdot[8] += qdot;

    qdot = q_f[38]-q_r[38];
    wdot[2] += qdot;
    wdot[4] -= qdot;
    wdot[7] += qdot;
    wdot[8] -= qdot;

    qdot = q_f[39]-q_r[39];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[9] -= qdot;
    wdot[10] += qdot;

    qdot = q_f[40]-q_r[40];
    wdot[2] += qdot;
    wdot[4] -= qdot;
    wdot[9] -= qdot;
    wdot[11] += qdot;

    qdot = q_f[41]-q_r[41];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[9] += qdot;
    wdot[11] -= qdot;

    qdot = q_f[42]-q_r[42];
    wdot[2] += qdot;
    wdot[3] += qdot;
    wdot[5] -= qdot;
    wdot[9] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[43]-q_r[43];
    wdot[1] -= qdot;
    wdot[2] += 2.000000 * qdot;
    wdot[9] -= qdot;
    wdot[12] += qdot;

    qdot = q_f[44]-q_r[44];
    wdot[3] += qdot;
    wdot[8] -= qdot;
    wdot[11] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[45]-q_r[45];
    wdot[3] += qdot;
    wdot[5] -= qdot;
    wdot[10] += qdot;
    wdot[11] -= qdot;

    qdot = q_f[46]-q_r[46];
    wdot[1] += qdot;
    wdot[5] -= qdot;
    wdot[11] -= qdot;
    wdot[13] += qdot;

    qdot = q_f[47]-q_r[47];
    wdot[1] -= qdot;
    wdot[5] += qdot;
    wdot[11] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[48]-q_r[48];
    wdot[2] += qdot;
    wdot[11] -= 2.000000 * qdot;
    wdot[14] += qdot;

    qdot = q_f[49]-q_r[49];
    wdot[5] += qdot;
    wdot[8] -= qdot;
    wdot[11] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[50]-q_r[50];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[10] += qdot;
    wdot[11] -= qdot;

    qdot = q_f[51]-q_r[51];
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[10] += qdot;
    wdot[11] -= qdot;

    qdot = q_f[52]-q_r[52];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[9] += qdot;
    wdot[11] -= qdot;

    qdot = q_f[53]-q_r[53];
    wdot[3] += qdot;
    wdot[6] -= qdot;
    wdot[9] -= qdot;
    wdot[11] += qdot;

    qdot = q_f[54]-q_r[54];
    wdot[7] -= qdot;
    wdot[8] += qdot;
    wdot[11] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[55]-q_r[55];
    wdot[2] += qdot;
    wdot[9] -= qdot;
    wdot[11] -= qdot;
    wdot[16] += qdot;

    qdot = q_f[56]-q_r[56];
    wdot[9] -= qdot;
    wdot[11] += 2.000000 * qdot;
    wdot[15] -= qdot;

    qdot = q_f[57]-q_r[57];
    wdot[9] += qdot;
    wdot[11] -= 2.000000 * qdot;
    wdot[15] += qdot;

    qdot = q_f[58]-q_r[58];
    wdot[1] -= qdot;
    wdot[3] += qdot;
    wdot[11] += qdot;
    wdot[15] -= qdot;

    qdot = q_f[59]-q_r[59];
    wdot[1] += qdot;
    wdot[3] -= qdot;
    wdot[11] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[60]-q_r[60];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[11] += qdot;
    wdot[15] -= qdot;

    qdot = q_f[61]-q_r[61];
    wdot[2] += qdot;
    wdot[4] -= qdot;
    wdot[11] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[62]-q_r[62];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[11] += qdot;
    wdot[15] -= qdot;

    qdot = q_f[63]-q_r[63];
    wdot[3] += qdot;
    wdot[6] -= qdot;
    wdot[11] -= qdot;
    wdot[15] += qdot;

    qdot = q_f[64]-q_r[64];
    wdot[3] += qdot;
    wdot[8] -= qdot;
    wdot[12] -= qdot;
    wdot[18] += qdot;

    qdot = q_f[65]-q_r[65];
    wdot[1] += qdot;
    wdot[5] -= qdot;
    wdot[12] -= qdot;
    wdot[18] += qdot;

    qdot = q_f[66]-q_r[66];
    wdot[1] -= qdot;
    wdot[5] += qdot;
    wdot[12] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[67]-q_r[67];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[12] -= qdot;
    wdot[18] += qdot;

    qdot = q_f[68]-q_r[68];
    wdot[2] -= qdot;
    wdot[3] += qdot;
    wdot[12] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[69]-q_r[69];
    wdot[11] -= qdot;
    wdot[12] += qdot;
    wdot[15] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[70]-q_r[70];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[12] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[71]-q_r[71];
    wdot[5] -= qdot;
    wdot[8] += qdot;
    wdot[12] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[72]-q_r[72];
    wdot[1] -= qdot;
    wdot[3] += qdot;
    wdot[12] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[73]-q_r[73];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[18] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[74]-q_r[74];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[12] += qdot;
    wdot[19] -= qdot;

    qdot = q_f[75]-q_r[75];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[10] -= qdot;
    wdot[19] += qdot;

    qdot = q_f[76]-q_r[76];
    wdot[5] -= qdot;
    wdot[8] += qdot;
    wdot[10] -= qdot;
    wdot[19] += qdot;

    qdot = q_f[77]-q_r[77];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[10] -= qdot;
    wdot[19] += qdot;

    qdot = q_f[78]-q_r[78];
    wdot[7] += qdot;
    wdot[8] -= qdot;
    wdot[10] -= qdot;
    wdot[19] += qdot;

    qdot = q_f[79]-q_r[79];
    wdot[1] -= qdot;
    wdot[3] += qdot;
    wdot[10] -= qdot;
    wdot[19] += qdot;

    qdot = q_f[80]-q_r[80];
    wdot[10] -= qdot;
    wdot[11] -= qdot;
    wdot[15] += qdot;
    wdot[19] += qdot;

    qdot = q_f[81]-q_r[81];
    wdot[5] -= qdot;
    wdot[8] += qdot;
    wdot[10] += qdot;
    wdot[13] -= qdot;

    qdot = q_f[82]-q_r[82];
    wdot[9] -= qdot;
    wdot[10] += qdot;
    wdot[12] += qdot;
    wdot[18] -= qdot;

    qdot = q_f[83]-q_r[83];
    wdot[11] -= qdot;
    wdot[13] += 2.000000 * qdot;
    wdot[20] -= qdot;

    qdot = q_f[84]-q_r[84];
    wdot[5] += qdot;
    wdot[8] -= qdot;
    wdot[20] -= qdot;
    wdot[21] += qdot;

    qdot = q_f[85]-q_r[85];
    wdot[5] += qdot;
    wdot[13] += 2.000000 * qdot;
    wdot[20] -= 2.000000 * qdot;

    qdot = q_f[86]-q_r[86];
    wdot[10] -= qdot;
    wdot[19] += qdot;
    wdot[20] -= qdot;
    wdot[21] += qdot;

    qdot = q_f[87]-q_r[87];
    wdot[3] += qdot;
    wdot[13] += qdot;
    wdot[21] -= qdot;

    qdot = q_f[88]-q_r[88];
    wdot[3] -= qdot;
    wdot[13] -= qdot;
    wdot[21] += qdot;

    qdot = q_f[89]-q_r[89];
    wdot[3] += qdot;
    wdot[5] -= qdot;
    wdot[22] -= qdot;
    wdot[23] += qdot;

    qdot = q_f[90]-q_r[90];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[22] -= qdot;
    wdot[23] += qdot;

    qdot = q_f[91]-q_r[91];
    wdot[1] += qdot;
    wdot[5] -= qdot;
    wdot[24] -= qdot;
    wdot[25] += qdot;

    qdot = q_f[92]-q_r[92];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[22] += qdot;
    wdot[24] -= qdot;

    qdot = q_f[93]-q_r[93];
    wdot[5] -= qdot;
    wdot[10] += qdot;
    wdot[19] += qdot;
    wdot[24] -= qdot;

    qdot = q_f[94]-q_r[94];
    wdot[5] -= qdot;
    wdot[8] += qdot;
    wdot[22] += qdot;
    wdot[24] -= qdot;

    qdot = q_f[95]-q_r[95];
    wdot[1] -= qdot;
    wdot[11] += qdot;
    wdot[16] -= qdot;
    wdot[19] += qdot;

    qdot = q_f[96]-q_r[96];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[16] -= qdot;
    wdot[24] += qdot;

    qdot = q_f[97]-q_r[97];
    wdot[3] += qdot;
    wdot[6] -= qdot;
    wdot[16] += qdot;
    wdot[24] -= qdot;

    qdot = q_f[98]-q_r[98];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[16] -= qdot;
    wdot[24] += qdot;

    qdot = q_f[99]-q_r[99];
    wdot[2] += qdot;
    wdot[4] -= qdot;
    wdot[16] += qdot;
    wdot[24] -= qdot;

    qdot = q_f[100]-q_r[100];
    wdot[16] -= qdot;
    wdot[20] -= qdot;
    wdot[21] += qdot;
    wdot[24] += qdot;

    qdot = q_f[101]-q_r[101];
    wdot[11] -= qdot;
    wdot[15] += qdot;
    wdot[16] -= qdot;
    wdot[24] += qdot;

    qdot = q_f[102]-q_r[102];
    wdot[11] += qdot;
    wdot[15] -= qdot;
    wdot[16] += qdot;
    wdot[24] -= qdot;

    qdot = q_f[103]-q_r[103];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[16] -= qdot;
    wdot[25] += qdot;

    qdot = q_f[104]-q_r[104];
    wdot[11] -= qdot;
    wdot[14] -= qdot;
    wdot[15] += qdot;
    wdot[16] += qdot;

    qdot = q_f[105]-q_r[105];
    wdot[5] -= qdot;
    wdot[8] += qdot;
    wdot[14] -= qdot;
    wdot[16] += qdot;

    qdot = q_f[106]-q_r[106];
    wdot[5] += qdot;
    wdot[8] -= qdot;
    wdot[14] += qdot;
    wdot[16] -= qdot;

    qdot = q_f[107]-q_r[107];
    wdot[5] += qdot;
    wdot[8] -= qdot;
    wdot[14] -= qdot;
    wdot[17] += qdot;

    qdot = q_f[108]-q_r[108];
    wdot[2] -= qdot;
    wdot[14] -= qdot;
    wdot[17] += qdot;

    qdot = q_f[109]-q_r[109];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[14] += qdot;
    wdot[17] -= qdot;

    qdot = q_f[110]-q_r[110];
    wdot[2] += qdot;
    wdot[4] -= qdot;
    wdot[14] -= qdot;
    wdot[17] += qdot;

    qdot = q_f[111]-q_r[111];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[14] += qdot;
    wdot[17] -= qdot;

    qdot = q_f[112]-q_r[112];
    wdot[9] -= qdot;
    wdot[11] += qdot;
    wdot[14] += qdot;
    wdot[17] -= qdot;

    qdot = q_f[113]-q_r[113];
    wdot[1] -= qdot;
    wdot[3] += qdot;
    wdot[14] += qdot;
    wdot[17] -= qdot;

    qdot = q_f[114]-q_r[114];
    wdot[11] -= qdot;
    wdot[14] += qdot;
    wdot[15] += qdot;
    wdot[17] -= qdot;

    qdot = q_f[115]-q_r[115];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[12] += 2.000000 * qdot;
    wdot[23] -= qdot;

    qdot = q_f[116]-q_r[116];
    wdot[5] -= qdot;
    wdot[18] += qdot;
    wdot[19] += qdot;
    wdot[23] -= qdot;

    qdot = q_f[117]-q_r[117];
    wdot[3] -= qdot;
    wdot[19] += 2.000000 * qdot;
    wdot[23] -= qdot;

    qdot = q_f[118]-q_r[118];
    wdot[2] -= qdot;
    wdot[9] += qdot;
    wdot[12] += qdot;
    wdot[23] -= qdot;

    qdot = q_f[119]-q_r[119];
    wdot[2] += qdot;
    wdot[9] -= qdot;
    wdot[12] -= qdot;
    wdot[23] += qdot;

    qdot = q_f[120]-q_r[120];
    wdot[3] += qdot;
    wdot[5] -= qdot;
    wdot[10] += qdot;
    wdot[12] += qdot;
    wdot[25] -= qdot;

    qdot = q_f[121]-q_r[121];
    wdot[10] -= qdot;
    wdot[19] += qdot;
    wdot[26] -= qdot;
    wdot[27] += qdot;

    qdot = q_f[122]-q_r[122];
    wdot[10] += qdot;
    wdot[19] -= qdot;
    wdot[26] += qdot;
    wdot[27] -= qdot;

    qdot = q_f[123]-q_r[123];
    wdot[3] += qdot;
    wdot[8] -= qdot;
    wdot[26] -= qdot;
    wdot[28] += qdot;

    qdot = q_f[124]-q_r[124];
    wdot[13] += qdot;
    wdot[20] -= qdot;
    wdot[26] -= qdot;
    wdot[28] += qdot;

    qdot = q_f[125]-q_r[125];
    wdot[5] += qdot;
    wdot[8] -= qdot;
    wdot[26] -= qdot;
    wdot[27] += qdot;

    qdot = q_f[126]-q_r[126];
    wdot[5] -= qdot;
    wdot[10] += qdot;
    wdot[25] += qdot;
    wdot[26] -= qdot;

    qdot = q_f[127]-q_r[127];
    wdot[3] += qdot;
    wdot[5] -= qdot;
    wdot[10] += qdot;
    wdot[22] += qdot;
    wdot[26] -= qdot;

    qdot = q_f[128]-q_r[128];
    wdot[2] -= qdot;
    wdot[26] -= qdot;
    wdot[27] += qdot;

    qdot = q_f[129]-q_r[129];
    wdot[11] += qdot;
    wdot[22] += qdot;
    wdot[26] -= qdot;

    qdot = q_f[130]-q_r[130];
    wdot[11] += qdot;
    wdot[24] += qdot;
    wdot[27] -= qdot;

    qdot = q_f[131]-q_r[131];
    wdot[11] -= qdot;
    wdot[24] -= qdot;
    wdot[27] += qdot;

    qdot = q_f[132]-q_r[132];
    wdot[1] -= qdot;
    wdot[14] += qdot;
    wdot[19] += qdot;
    wdot[27] -= qdot;

    qdot = q_f[133]-q_r[133];
    wdot[2] -= qdot;
    wdot[11] += qdot;
    wdot[16] += qdot;
    wdot[27] -= qdot;

    qdot = q_f[134]-q_r[134];
    wdot[2] += qdot;
    wdot[11] -= qdot;
    wdot[16] -= qdot;
    wdot[27] += qdot;

    qdot = q_f[135]-q_r[135];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[26] += qdot;
    wdot[27] -= qdot;

    qdot = q_f[136]-q_r[136];
    wdot[2] += qdot;
    wdot[4] -= qdot;
    wdot[26] -= qdot;
    wdot[27] += qdot;

    qdot = q_f[137]-q_r[137];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[26] += qdot;
    wdot[27] -= qdot;

    qdot = q_f[138]-q_r[138];
    wdot[1] -= qdot;
    wdot[3] += qdot;
    wdot[26] += qdot;
    wdot[27] -= qdot;

    qdot = q_f[139]-q_r[139];
    wdot[11] -= qdot;
    wdot[15] += qdot;
    wdot[26] += qdot;
    wdot[27] -= qdot;

    qdot = q_f[140]-q_r[140];
    wdot[5] -= qdot;
    wdot[8] += qdot;
    wdot[27] += qdot;
    wdot[29] -= qdot;

    qdot = q_f[141]-q_r[141];
    wdot[2] += qdot;
    wdot[27] += qdot;
    wdot[29] -= qdot;

    qdot = q_f[142]-q_r[142];
    wdot[2] -= qdot;
    wdot[27] -= qdot;
    wdot[29] += qdot;

    qdot = q_f[143]-q_r[143];
    wdot[11] += qdot;
    wdot[16] += qdot;
    wdot[30] -= qdot;

    qdot = q_f[144]-q_r[144];
    wdot[11] -= qdot;
    wdot[16] -= qdot;
    wdot[30] += qdot;

    qdot = q_f[145]-q_r[145];
    wdot[5] += qdot;
    wdot[8] -= qdot;
    wdot[30] -= qdot;
    wdot[31] += qdot;

    qdot = q_f[146]-q_r[146];
    wdot[5] -= qdot;
    wdot[8] += qdot;
    wdot[27] += qdot;
    wdot[30] -= qdot;

    qdot = q_f[147]-q_r[147];
    wdot[2] += qdot;
    wdot[27] += qdot;
    wdot[30] -= qdot;

    qdot = q_f[148]-q_r[148];
    wdot[2] -= qdot;
    wdot[27] -= qdot;
    wdot[30] += qdot;

    qdot = q_f[149]-q_r[149];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[30] += qdot;
    wdot[31] -= qdot;

    qdot = q_f[150]-q_r[150];
    wdot[7] += qdot;
    wdot[8] -= qdot;
    wdot[30] += qdot;
    wdot[31] -= qdot;

    qdot = q_f[151]-q_r[151];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[30] += qdot;
    wdot[31] -= qdot;

    qdot = q_f[152]-q_r[152];
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[29] += qdot;
    wdot[31] -= qdot;

    qdot = q_f[153]-q_r[153];
    wdot[11] -= qdot;
    wdot[15] += qdot;
    wdot[29] += qdot;
    wdot[31] -= qdot;

    qdot = q_f[154]-q_r[154];
    wdot[11] -= qdot;
    wdot[15] += qdot;
    wdot[30] += qdot;
    wdot[31] -= qdot;

    qdot = q_f[155]-q_r[155];
    wdot[1] -= qdot;
    wdot[3] += qdot;
    wdot[29] += qdot;
    wdot[31] -= qdot;

    qdot = q_f[156]-q_r[156];
    wdot[7] += qdot;
    wdot[8] -= qdot;
    wdot[29] += qdot;
    wdot[31] -= qdot;

    qdot = q_f[157]-q_r[157];
    wdot[1] -= qdot;
    wdot[3] += qdot;
    wdot[30] += qdot;
    wdot[31] -= qdot;

    qdot = q_f[158]-q_r[158];
    wdot[5] -= qdot;
    wdot[8] += qdot;
    wdot[29] += qdot;
    wdot[31] -= qdot;

    qdot = q_f[159]-q_r[159];
    wdot[5] += qdot;
    wdot[8] -= qdot;
    wdot[29] -= qdot;
    wdot[31] += qdot;

    qdot = q_f[160]-q_r[160];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[29] += qdot;
    wdot[31] -= qdot;

    qdot = q_f[161]-q_r[161];
    wdot[2] += qdot;
    wdot[4] -= qdot;
    wdot[29] -= qdot;
    wdot[31] += qdot;

    qdot = q_f[162]-q_r[162];
    wdot[10] += qdot;
    wdot[24] += qdot;
    wdot[28] -= qdot;

    qdot = q_f[163]-q_r[163];
    wdot[5] += qdot;
    wdot[29] += qdot;
    wdot[32] -= qdot;

    qdot = q_f[164]-q_r[164];
    wdot[5] -= qdot;
    wdot[29] -= qdot;
    wdot[32] += qdot;

    qdot = q_f[165]-q_r[165];
    wdot[5] += qdot;
    wdot[30] += qdot;
    wdot[33] -= qdot;

    qdot = q_f[166]-q_r[166];
    wdot[5] -= qdot;
    wdot[30] -= qdot;
    wdot[33] += qdot;

    return;
}

void comp_k_f(amrex::Real *  tc, amrex::Real invT, amrex::Real *  k_f)
{
    for (int i=0; i<167; ++i) {
        k_f[i] = prefactor_units[i] * fwd_A[i]
                    * exp(fwd_beta[i] * tc[0] - activation_units[i] * fwd_Ea[i] * invT);
    };
    return;
}

void comp_Kc(amrex::Real *  tc, amrex::Real invT, amrex::Real *  Kc)
{
    /*compute the Gibbs free energy */
    amrex::Real g_RT[34];
    gibbs(g_RT, tc);

    Kc[0] = 2.000000*g_RT[3] - g_RT[7];
    Kc[1] = g_RT[2] + g_RT[5] - g_RT[8];
    Kc[2] = g_RT[2] + g_RT[11] - g_RT[15];
    Kc[3] = 2.000000*g_RT[11] - g_RT[17];
    Kc[4] = g_RT[1] + g_RT[12] - g_RT[18];
    Kc[5] = -g_RT[2] - g_RT[10] + g_RT[13];
    Kc[6] = -g_RT[2] - g_RT[22] + g_RT[24];
    Kc[7] = g_RT[2] - g_RT[14] + g_RT[16];
    Kc[8] = -g_RT[4] + g_RT[16] - g_RT[22];
    Kc[9] = g_RT[1] + g_RT[2] - g_RT[3];
    Kc[10] = 2.000000*g_RT[1] - g_RT[5];
    Kc[11] = g_RT[2] + g_RT[3] - g_RT[6];
    Kc[12] = -g_RT[2] - g_RT[12] + g_RT[19];
    Kc[13] = g_RT[2] + g_RT[12] - g_RT[19];
    Kc[14] = -g_RT[5] - g_RT[11] + g_RT[20];
    Kc[15] = g_RT[5] + g_RT[11] - g_RT[20];
    Kc[16] = g_RT[1] - g_RT[2] - g_RT[3] + g_RT[4];
    Kc[17] = -g_RT[1] + g_RT[2] + g_RT[3] - g_RT[4];
    Kc[18] = -g_RT[2] + g_RT[3] + g_RT[4] - g_RT[6];
    Kc[19] = g_RT[2] - g_RT[3] - g_RT[4] + g_RT[6];
    Kc[20] = g_RT[1] - 2.000000*g_RT[3] + g_RT[6];
    Kc[21] = -g_RT[1] + 2.000000*g_RT[3] - g_RT[6];
    Kc[22] = -g_RT[1] + g_RT[2] - g_RT[3] + g_RT[5];
    Kc[23] = g_RT[1] - g_RT[2] + g_RT[3] - g_RT[5];
    Kc[24] = g_RT[2] - 2.000000*g_RT[3] + g_RT[8];
    Kc[25] = -g_RT[5] - g_RT[7] + 2.000000*g_RT[8];
    Kc[26] = g_RT[2] - g_RT[4] - g_RT[5] + g_RT[8];
    Kc[27] = g_RT[3] - g_RT[5] - g_RT[6] + g_RT[8];
    Kc[28] = -g_RT[3] + g_RT[5] + g_RT[6] - g_RT[8];
    Kc[29] = -g_RT[5] - g_RT[7] + 2.000000*g_RT[8];
    Kc[30] = g_RT[1] - g_RT[3] - g_RT[5] + g_RT[8];
    Kc[31] = -g_RT[1] + g_RT[3] + g_RT[5] - g_RT[8];
    Kc[32] = g_RT[3] - g_RT[6] + g_RT[7] - g_RT[8];
    Kc[33] = g_RT[2] - g_RT[3] - g_RT[6] + g_RT[7];
    Kc[34] = g_RT[3] - g_RT[6] + g_RT[7] - g_RT[8];
    Kc[35] = -g_RT[3] + g_RT[6] - g_RT[7] + g_RT[8];
    Kc[36] = g_RT[1] - g_RT[3] + g_RT[7] - g_RT[8];
    Kc[37] = g_RT[2] - g_RT[4] + g_RT[7] - g_RT[8];
    Kc[38] = -g_RT[2] + g_RT[4] - g_RT[7] + g_RT[8];
    Kc[39] = -g_RT[2] + g_RT[3] + g_RT[9] - g_RT[10];
    Kc[40] = -g_RT[2] + g_RT[4] + g_RT[9] - g_RT[11];
    Kc[41] = g_RT[2] - g_RT[4] - g_RT[9] + g_RT[11];
    Kc[42] = -g_RT[2] - g_RT[3] + g_RT[5] + g_RT[9] - g_RT[12];
    Kc[43] = g_RT[1] - 2.000000*g_RT[2] + g_RT[9] - g_RT[12];
    Kc[44] = -g_RT[3] + g_RT[8] + g_RT[11] - g_RT[13];
    Kc[45] = -g_RT[3] + g_RT[5] - g_RT[10] + g_RT[11];
    Kc[46] = -g_RT[1] + g_RT[5] + g_RT[11] - g_RT[13];
    Kc[47] = g_RT[1] - g_RT[5] - g_RT[11] + g_RT[13];
    Kc[48] = -g_RT[2] + 2.000000*g_RT[11] - g_RT[14];
    Kc[49] = -g_RT[5] + g_RT[8] + g_RT[11] - g_RT[15];
    Kc[50] = g_RT[1] - g_RT[2] - g_RT[10] + g_RT[11];
    Kc[51] = g_RT[3] - g_RT[4] - g_RT[10] + g_RT[11];
    Kc[52] = g_RT[3] - g_RT[6] - g_RT[9] + g_RT[11];
    Kc[53] = -g_RT[3] + g_RT[6] + g_RT[9] - g_RT[11];
    Kc[54] = g_RT[7] - g_RT[8] + g_RT[11] - g_RT[15];
    Kc[55] = -g_RT[2] + g_RT[9] + g_RT[11] - g_RT[16];
    Kc[56] = g_RT[9] - 2.000000*g_RT[11] + g_RT[15];
    Kc[57] = -g_RT[9] + 2.000000*g_RT[11] - g_RT[15];
    Kc[58] = g_RT[1] - g_RT[3] - g_RT[11] + g_RT[15];
    Kc[59] = -g_RT[1] + g_RT[3] + g_RT[11] - g_RT[15];
    Kc[60] = g_RT[2] - g_RT[4] - g_RT[11] + g_RT[15];
    Kc[61] = -g_RT[2] + g_RT[4] + g_RT[11] - g_RT[15];
    Kc[62] = g_RT[3] - g_RT[6] - g_RT[11] + g_RT[15];
    Kc[63] = -g_RT[3] + g_RT[6] + g_RT[11] - g_RT[15];
    Kc[64] = -g_RT[3] + g_RT[8] + g_RT[12] - g_RT[18];
    Kc[65] = -g_RT[1] + g_RT[5] + g_RT[12] - g_RT[18];
    Kc[66] = g_RT[1] - g_RT[5] - g_RT[12] + g_RT[18];
    Kc[67] = -g_RT[2] + g_RT[3] + g_RT[12] - g_RT[18];
    Kc[68] = g_RT[2] - g_RT[3] - g_RT[12] + g_RT[18];
    Kc[69] = g_RT[11] - g_RT[12] - g_RT[15] + g_RT[19];
    Kc[70] = g_RT[2] - g_RT[4] - g_RT[12] + g_RT[19];
    Kc[71] = g_RT[5] - g_RT[8] - g_RT[12] + g_RT[19];
    Kc[72] = g_RT[1] - g_RT[3] - g_RT[12] + g_RT[19];
    Kc[73] = g_RT[1] - g_RT[2] - g_RT[18] + g_RT[19];
    Kc[74] = g_RT[3] - g_RT[6] - g_RT[12] + g_RT[19];
    Kc[75] = g_RT[2] - g_RT[4] + g_RT[10] - g_RT[19];
    Kc[76] = g_RT[5] - g_RT[8] + g_RT[10] - g_RT[19];
    Kc[77] = g_RT[3] - g_RT[6] + g_RT[10] - g_RT[19];
    Kc[78] = -g_RT[7] + g_RT[8] + g_RT[10] - g_RT[19];
    Kc[79] = g_RT[1] - g_RT[3] + g_RT[10] - g_RT[19];
    Kc[80] = g_RT[10] + g_RT[11] - g_RT[15] - g_RT[19];
    Kc[81] = g_RT[5] - g_RT[8] - g_RT[10] + g_RT[13];
    Kc[82] = g_RT[9] - g_RT[10] - g_RT[12] + g_RT[18];
    Kc[83] = g_RT[11] - 2.000000*g_RT[13] + g_RT[20];
    Kc[84] = -g_RT[5] + g_RT[8] + g_RT[20] - g_RT[21];
    Kc[85] = -g_RT[5] - 2.000000*g_RT[13] + 2.000000*g_RT[20];
    Kc[86] = g_RT[10] - g_RT[19] + g_RT[20] - g_RT[21];
    Kc[87] = -g_RT[3] - g_RT[13] + g_RT[21];
    Kc[88] = g_RT[3] + g_RT[13] - g_RT[21];
    Kc[89] = -g_RT[3] + g_RT[5] + g_RT[22] - g_RT[23];
    Kc[90] = g_RT[1] - g_RT[2] + g_RT[22] - g_RT[23];
    Kc[91] = -g_RT[1] + g_RT[5] + g_RT[24] - g_RT[25];
    Kc[92] = g_RT[2] - g_RT[4] - g_RT[22] + g_RT[24];
    Kc[93] = g_RT[5] - g_RT[10] - g_RT[19] + g_RT[24];
    Kc[94] = g_RT[5] - g_RT[8] - g_RT[22] + g_RT[24];
    Kc[95] = g_RT[1] - g_RT[11] + g_RT[16] - g_RT[19];
    Kc[96] = g_RT[3] - g_RT[6] + g_RT[16] - g_RT[24];
    Kc[97] = -g_RT[3] + g_RT[6] - g_RT[16] + g_RT[24];
    Kc[98] = g_RT[2] - g_RT[4] + g_RT[16] - g_RT[24];
    Kc[99] = -g_RT[2] + g_RT[4] - g_RT[16] + g_RT[24];
    Kc[100] = g_RT[16] + g_RT[20] - g_RT[21] - g_RT[24];
    Kc[101] = g_RT[11] - g_RT[15] + g_RT[16] - g_RT[24];
    Kc[102] = -g_RT[11] + g_RT[15] - g_RT[16] + g_RT[24];
    Kc[103] = g_RT[1] - g_RT[2] + g_RT[16] - g_RT[25];
    Kc[104] = g_RT[11] + g_RT[14] - g_RT[15] - g_RT[16];
    Kc[105] = g_RT[5] - g_RT[8] + g_RT[14] - g_RT[16];
    Kc[106] = -g_RT[5] + g_RT[8] - g_RT[14] + g_RT[16];
    Kc[107] = -g_RT[5] + g_RT[8] + g_RT[14] - g_RT[17];
    Kc[108] = g_RT[2] + g_RT[14] - g_RT[17];
    Kc[109] = g_RT[2] - g_RT[4] - g_RT[14] + g_RT[17];
    Kc[110] = -g_RT[2] + g_RT[4] + g_RT[14] - g_RT[17];
    Kc[111] = g_RT[3] - g_RT[6] - g_RT[14] + g_RT[17];
    Kc[112] = g_RT[9] - g_RT[11] - g_RT[14] + g_RT[17];
    Kc[113] = g_RT[1] - g_RT[3] - g_RT[14] + g_RT[17];
    Kc[114] = g_RT[11] - g_RT[14] - g_RT[15] + g_RT[17];
    Kc[115] = g_RT[1] - g_RT[2] - 2.000000*g_RT[12] + g_RT[23];
    Kc[116] = g_RT[5] - g_RT[18] - g_RT[19] + g_RT[23];
    Kc[117] = g_RT[3] - 2.000000*g_RT[19] + g_RT[23];
    Kc[118] = g_RT[2] - g_RT[9] - g_RT[12] + g_RT[23];
    Kc[119] = -g_RT[2] + g_RT[9] + g_RT[12] - g_RT[23];
    Kc[120] = -g_RT[3] + g_RT[5] - g_RT[10] - g_RT[12] + g_RT[25];
    Kc[121] = g_RT[10] - g_RT[19] + g_RT[26] - g_RT[27];
    Kc[122] = -g_RT[10] + g_RT[19] - g_RT[26] + g_RT[27];
    Kc[123] = -g_RT[3] + g_RT[8] + g_RT[26] - g_RT[28];
    Kc[124] = -g_RT[13] + g_RT[20] + g_RT[26] - g_RT[28];
    Kc[125] = -g_RT[5] + g_RT[8] + g_RT[26] - g_RT[27];
    Kc[126] = g_RT[5] - g_RT[10] - g_RT[25] + g_RT[26];
    Kc[127] = -g_RT[3] + g_RT[5] - g_RT[10] - g_RT[22] + g_RT[26];
    Kc[128] = g_RT[2] + g_RT[26] - g_RT[27];
    Kc[129] = -g_RT[11] - g_RT[22] + g_RT[26];
    Kc[130] = -g_RT[11] - g_RT[24] + g_RT[27];
    Kc[131] = g_RT[11] + g_RT[24] - g_RT[27];
    Kc[132] = g_RT[1] - g_RT[14] - g_RT[19] + g_RT[27];
    Kc[133] = g_RT[2] - g_RT[11] - g_RT[16] + g_RT[27];
    Kc[134] = -g_RT[2] + g_RT[11] + g_RT[16] - g_RT[27];
    Kc[135] = g_RT[2] - g_RT[4] - g_RT[26] + g_RT[27];
    Kc[136] = -g_RT[2] + g_RT[4] + g_RT[26] - g_RT[27];
    Kc[137] = g_RT[3] - g_RT[6] - g_RT[26] + g_RT[27];
    Kc[138] = g_RT[1] - g_RT[3] - g_RT[26] + g_RT[27];
    Kc[139] = g_RT[11] - g_RT[15] - g_RT[26] + g_RT[27];
    Kc[140] = g_RT[5] - g_RT[8] - g_RT[27] + g_RT[29];
    Kc[141] = -g_RT[2] - g_RT[27] + g_RT[29];
    Kc[142] = g_RT[2] + g_RT[27] - g_RT[29];
    Kc[143] = -g_RT[11] - g_RT[16] + g_RT[30];
    Kc[144] = g_RT[11] + g_RT[16] - g_RT[30];
    Kc[145] = -g_RT[5] + g_RT[8] + g_RT[30] - g_RT[31];
    Kc[146] = g_RT[5] - g_RT[8] - g_RT[27] + g_RT[30];
    Kc[147] = -g_RT[2] - g_RT[27] + g_RT[30];
    Kc[148] = g_RT[2] + g_RT[27] - g_RT[30];
    Kc[149] = g_RT[3] - g_RT[6] - g_RT[30] + g_RT[31];
    Kc[150] = -g_RT[7] + g_RT[8] - g_RT[30] + g_RT[31];
    Kc[151] = g_RT[2] - g_RT[4] - g_RT[30] + g_RT[31];
    Kc[152] = g_RT[3] - g_RT[6] - g_RT[29] + g_RT[31];
    Kc[153] = g_RT[11] - g_RT[15] - g_RT[29] + g_RT[31];
    Kc[154] = g_RT[11] - g_RT[15] - g_RT[30] + g_RT[31];
    Kc[155] = g_RT[1] - g_RT[3] - g_RT[29] + g_RT[31];
    Kc[156] = -g_RT[7] + g_RT[8] - g_RT[29] + g_RT[31];
    Kc[157] = g_RT[1] - g_RT[3] - g_RT[30] + g_RT[31];
    Kc[158] = g_RT[5] - g_RT[8] - g_RT[29] + g_RT[31];
    Kc[159] = -g_RT[5] + g_RT[8] + g_RT[29] - g_RT[31];
    Kc[160] = g_RT[2] - g_RT[4] - g_RT[29] + g_RT[31];
    Kc[161] = -g_RT[2] + g_RT[4] + g_RT[29] - g_RT[31];
    Kc[162] = -g_RT[10] - g_RT[24] + g_RT[28];
    Kc[163] = -g_RT[5] - g_RT[29] + g_RT[32];
    Kc[164] = g_RT[5] + g_RT[29] - g_RT[32];
    Kc[165] = -g_RT[5] - g_RT[30] + g_RT[33];
    Kc[166] = g_RT[5] + g_RT[30] - g_RT[33];

    for (int i=0; i<167; ++i) {
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
    Kc[6] *= refC;
    Kc[7] *= refCinv;
    Kc[8] *= refC;
    Kc[9] *= refCinv;
    Kc[10] *= refCinv;
    Kc[11] *= refCinv;
    Kc[12] *= refC;
    Kc[13] *= refCinv;
    Kc[14] *= refC;
    Kc[15] *= refCinv;
    Kc[42] *= refC;
    Kc[43] *= refC;
    Kc[85] *= refC;
    Kc[87] *= refC;
    Kc[88] *= refCinv;
    Kc[108] *= refCinv;
    Kc[115] *= refC;
    Kc[120] *= refC;
    Kc[127] *= refC;
    Kc[128] *= refCinv;
    Kc[129] *= refC;
    Kc[130] *= refC;
    Kc[131] *= refCinv;
    Kc[141] *= refC;
    Kc[142] *= refCinv;
    Kc[143] *= refC;
    Kc[144] *= refCinv;
    Kc[147] *= refC;
    Kc[148] *= refCinv;
    Kc[162] *= refC;
    Kc[163] *= refC;
    Kc[164] *= refCinv;
    Kc[165] *= refC;
    Kc[166] *= refCinv;

    return;
}

void comp_qfqr(amrex::Real *  qf, amrex::Real *  qr, amrex::Real *  sc, amrex::Real * qss_sc, amrex::Real *  tc, amrex::Real invT)
{

    /*reaction 1: 2.000000 OH (+M) <=> H2O2 (+M) */
    qf[0] = pow(sc[3], 2.000000);
    qr[0] = sc[7];

    /*reaction 2: H + O2 (+M) <=> HO2 (+M) */
    qf[1] = sc[2]*sc[5];
    qr[1] = sc[8];

    /*reaction 3: CH3 + H (+M) => CH4 (+M) */
    qf[2] = sc[2]*sc[11];
    qr[2] = 0.0;

    /*reaction 4: 2.000000 CH3 (+M) <=> C2H6 (+M) */
    qf[3] = pow(sc[11], 2.000000);
    qr[3] = sc[17];

    /*reaction 5: CO + O (+M) => CO2 (+M) */
    qf[4] = sc[1]*sc[12];
    qr[4] = 0.0;

    /*reaction 6: CH3O (+M) => CH2O + H (+M) */
    qf[5] = sc[13];
    qr[5] = 0.0;

    /*reaction 7: C2H3 (+M) => H + C2H2 (+M) */
    qf[6] = sc[24];
    qr[6] = 0.0;

    /*reaction 8: H + C2H4 (+M) <=> C2H5 (+M) */
    qf[7] = sc[2]*sc[16];
    qr[7] = sc[14];

    /*reaction 9: C2H4 (+M) => C2H2 + H2 (+M) */
    qf[8] = sc[16];
    qr[8] = 0.0;

    /*reaction 10: O + H + M => OH + M */
    qf[9] = sc[1]*sc[2];
    qr[9] = 0.0;

    /*reaction 11: 2.000000 O + M => O2 + M */
    qf[10] = pow(sc[1], 2.000000);
    qr[10] = 0.0;

    /*reaction 12: H + OH + M => H2O + M */
    qf[11] = sc[2]*sc[3];
    qr[11] = 0.0;

    /*reaction 13: HCO + M => H + CO + M */
    qf[12] = sc[19];
    qr[12] = 0.0;

    /*reaction 14: H + CO + M => HCO + M */
    qf[13] = sc[2]*sc[12];
    qr[13] = 0.0;

    /*reaction 15: CH3O2 + M => CH3 + O2 + M */
    qf[14] = sc[20];
    qr[14] = 0.0;

    /*reaction 16: CH3 + O2 + M => CH3O2 + M */
    qf[15] = sc[5]*sc[11];
    qr[15] = 0.0;

    /*reaction 17: O + H2 => H + OH */
    qf[16] = sc[1]*sc[4];
    qr[16] = 0.0;

    /*reaction 18: H + OH => O + H2 */
    qf[17] = sc[2]*sc[3];
    qr[17] = 0.0;

    /*reaction 19: OH + H2 => H + H2O */
    qf[18] = sc[3]*sc[4];
    qr[18] = 0.0;

    /*reaction 20: H + H2O => OH + H2 */
    qf[19] = sc[2]*sc[6];
    qr[19] = 0.0;

    /*reaction 21: O + H2O => 2.000000 OH */
    qf[20] = sc[1]*sc[6];
    qr[20] = 0.0;

    /*reaction 22: 2.000000 OH => O + H2O */
    qf[21] = pow(sc[3], 2.000000);
    qr[21] = 0.0;

    /*reaction 23: H + O2 => O + OH */
    qf[22] = sc[2]*sc[5];
    qr[22] = 0.0;

    /*reaction 24: O + OH => H + O2 */
    qf[23] = sc[1]*sc[3];
    qr[23] = 0.0;

    /*reaction 25: HO2 + H => 2.000000 OH */
    qf[24] = sc[2]*sc[8];
    qr[24] = 0.0;

    /*reaction 26: 2.000000 HO2 => H2O2 + O2 */
    qf[25] = pow(sc[8], 2.000000);
    qr[25] = 0.0;

    /*reaction 27: HO2 + H => H2 + O2 */
    qf[26] = sc[2]*sc[8];
    qr[26] = 0.0;

    /*reaction 28: HO2 + OH => H2O + O2 */
    qf[27] = sc[3]*sc[8];
    qr[27] = 0.0;

    /*reaction 29: H2O + O2 => HO2 + OH */
    qf[28] = sc[5]*sc[6];
    qr[28] = 0.0;

    /*reaction 30: 2.000000 HO2 => H2O2 + O2 */
    qf[29] = pow(sc[8], 2.000000);
    qr[29] = 0.0;

    /*reaction 31: HO2 + O => OH + O2 */
    qf[30] = sc[1]*sc[8];
    qr[30] = 0.0;

    /*reaction 32: OH + O2 => HO2 + O */
    qf[31] = sc[3]*sc[5];
    qr[31] = 0.0;

    /*reaction 33: H2O2 + OH => H2O + HO2 */
    qf[32] = sc[3]*sc[7];
    qr[32] = 0.0;

    /*reaction 34: H2O2 + H => H2O + OH */
    qf[33] = sc[2]*sc[7];
    qr[33] = 0.0;

    /*reaction 35: H2O2 + OH => H2O + HO2 */
    qf[34] = sc[3]*sc[7];
    qr[34] = 0.0;

    /*reaction 36: H2O + HO2 => H2O2 + OH */
    qf[35] = sc[6]*sc[8];
    qr[35] = 0.0;

    /*reaction 37: H2O2 + O => OH + HO2 */
    qf[36] = sc[1]*sc[7];
    qr[36] = 0.0;

    /*reaction 38: H2O2 + H => H2 + HO2 */
    qf[37] = sc[2]*sc[7];
    qr[37] = 0.0;

    /*reaction 39: H2 + HO2 => H2O2 + H */
    qf[38] = sc[4]*sc[8];
    qr[38] = 0.0;

    /*reaction 40: CH2GSG + OH => CH2O + H */
    qf[39] = sc[3]*sc[9];
    qr[39] = 0.0;

    /*reaction 41: CH2GSG + H2 => CH3 + H */
    qf[40] = sc[4]*sc[9];
    qr[40] = 0.0;

    /*reaction 42: CH3 + H => CH2GSG + H2 */
    qf[41] = sc[2]*sc[11];
    qr[41] = 0.0;

    /*reaction 43: CH2GSG + O2 => CO + OH + H */
    qf[42] = sc[5]*sc[9];
    qr[42] = 0.0;

    /*reaction 44: CH2GSG + O => CO + 2.000000 H */
    qf[43] = sc[1]*sc[9];
    qr[43] = 0.0;

    /*reaction 45: CH3 + HO2 => CH3O + OH */
    qf[44] = sc[8]*sc[11];
    qr[44] = 0.0;

    /*reaction 46: CH3 + O2 => CH2O + OH */
    qf[45] = sc[5]*sc[11];
    qr[45] = 0.0;

    /*reaction 47: CH3 + O2 => CH3O + O */
    qf[46] = sc[5]*sc[11];
    qr[46] = 0.0;

    /*reaction 48: CH3O + O => CH3 + O2 */
    qf[47] = sc[1]*sc[13];
    qr[47] = 0.0;

    /*reaction 49: 2.000000 CH3 <=> H + C2H5 */
    qf[48] = pow(sc[11], 2.000000);
    qr[48] = sc[2]*sc[14];

    /*reaction 50: CH3 + HO2 => CH4 + O2 */
    qf[49] = sc[8]*sc[11];
    qr[49] = 0.0;

    /*reaction 51: CH3 + O => CH2O + H */
    qf[50] = sc[1]*sc[11];
    qr[50] = 0.0;

    /*reaction 52: CH3 + OH => CH2O + H2 */
    qf[51] = sc[3]*sc[11];
    qr[51] = 0.0;

    /*reaction 53: CH3 + OH => CH2GSG + H2O */
    qf[52] = sc[3]*sc[11];
    qr[52] = 0.0;

    /*reaction 54: CH2GSG + H2O => CH3 + OH */
    qf[53] = sc[6]*sc[9];
    qr[53] = 0.0;

    /*reaction 55: CH3 + H2O2 => CH4 + HO2 */
    qf[54] = sc[7]*sc[11];
    qr[54] = 0.0;

    /*reaction 56: CH2GSG + CH3 => C2H4 + H */
    qf[55] = sc[9]*sc[11];
    qr[55] = 0.0;

    /*reaction 57: CH2GSG + CH4 => 2.000000 CH3 */
    qf[56] = sc[9]*sc[15];
    qr[56] = 0.0;

    /*reaction 58: 2.000000 CH3 => CH2GSG + CH4 */
    qf[57] = pow(sc[11], 2.000000);
    qr[57] = 0.0;

    /*reaction 59: CH4 + O => CH3 + OH */
    qf[58] = sc[1]*sc[15];
    qr[58] = 0.0;

    /*reaction 60: CH3 + OH => CH4 + O */
    qf[59] = sc[3]*sc[11];
    qr[59] = 0.0;

    /*reaction 61: CH4 + H => CH3 + H2 */
    qf[60] = sc[2]*sc[15];
    qr[60] = 0.0;

    /*reaction 62: CH3 + H2 => CH4 + H */
    qf[61] = sc[4]*sc[11];
    qr[61] = 0.0;

    /*reaction 63: CH4 + OH => CH3 + H2O */
    qf[62] = sc[3]*sc[15];
    qr[62] = 0.0;

    /*reaction 64: CH3 + H2O => CH4 + OH */
    qf[63] = sc[6]*sc[11];
    qr[63] = 0.0;

    /*reaction 65: CO + HO2 => CO2 + OH */
    qf[64] = sc[8]*sc[12];
    qr[64] = 0.0;

    /*reaction 66: CO + O2 => CO2 + O */
    qf[65] = sc[5]*sc[12];
    qr[65] = 0.0;

    /*reaction 67: CO2 + O => CO + O2 */
    qf[66] = sc[1]*sc[18];
    qr[66] = 0.0;

    /*reaction 68: CO + OH => CO2 + H */
    qf[67] = sc[3]*sc[12];
    qr[67] = 0.0;

    /*reaction 69: CO2 + H => CO + OH */
    qf[68] = sc[2]*sc[18];
    qr[68] = 0.0;

    /*reaction 70: HCO + CH3 => CH4 + CO */
    qf[69] = sc[11]*sc[19];
    qr[69] = 0.0;

    /*reaction 71: HCO + H => CO + H2 */
    qf[70] = sc[2]*sc[19];
    qr[70] = 0.0;

    /*reaction 72: HCO + O2 => CO + HO2 */
    qf[71] = sc[5]*sc[19];
    qr[71] = 0.0;

    /*reaction 73: HCO + O => CO + OH */
    qf[72] = sc[1]*sc[19];
    qr[72] = 0.0;

    /*reaction 74: HCO + O => CO2 + H */
    qf[73] = sc[1]*sc[19];
    qr[73] = 0.0;

    /*reaction 75: HCO + OH => CO + H2O */
    qf[74] = sc[3]*sc[19];
    qr[74] = 0.0;

    /*reaction 76: CH2O + H => HCO + H2 */
    qf[75] = sc[2]*sc[10];
    qr[75] = 0.0;

    /*reaction 77: CH2O + O2 => HCO + HO2 */
    qf[76] = sc[5]*sc[10];
    qr[76] = 0.0;

    /*reaction 78: CH2O + OH => HCO + H2O */
    qf[77] = sc[3]*sc[10];
    qr[77] = 0.0;

    /*reaction 79: CH2O + HO2 => HCO + H2O2 */
    qf[78] = sc[8]*sc[10];
    qr[78] = 0.0;

    /*reaction 80: CH2O + O => HCO + OH */
    qf[79] = sc[1]*sc[10];
    qr[79] = 0.0;

    /*reaction 81: CH2O + CH3 => HCO + CH4 */
    qf[80] = sc[10]*sc[11];
    qr[80] = 0.0;

    /*reaction 82: CH3O + O2 => CH2O + HO2 */
    qf[81] = sc[5]*sc[13];
    qr[81] = 0.0;

    /*reaction 83: CH2GSG + CO2 => CH2O + CO */
    qf[82] = sc[9]*sc[18];
    qr[82] = 0.0;

    /*reaction 84: CH3O2 + CH3 => 2.000000 CH3O */
    qf[83] = sc[11]*sc[20];
    qr[83] = 0.0;

    /*reaction 85: CH3O2 + HO2 => CH3O2H + O2 */
    qf[84] = sc[8]*sc[20];
    qr[84] = 0.0;

    /*reaction 86: 2.000000 CH3O2 => O2 + 2.000000 CH3O */
    qf[85] = pow(sc[20], 2.000000);
    qr[85] = 0.0;

    /*reaction 87: CH3O2 + CH2O => CH3O2H + HCO */
    qf[86] = sc[10]*sc[20];
    qr[86] = 0.0;

    /*reaction 88: CH3O2H => CH3O + OH */
    qf[87] = sc[21];
    qr[87] = 0.0;

    /*reaction 89: CH3O + OH => CH3O2H */
    qf[88] = sc[3]*sc[13];
    qr[88] = 0.0;

    /*reaction 90: C2H2 + O2 => HCCO + OH */
    qf[89] = sc[5]*sc[22];
    qr[89] = 0.0;

    /*reaction 91: C2H2 + O => HCCO + H */
    qf[90] = sc[1]*sc[22];
    qr[90] = 0.0;

    /*reaction 92: C2H3 + O2 => CH2CHO + O */
    qf[91] = sc[5]*sc[24];
    qr[91] = 0.0;

    /*reaction 93: C2H3 + H => C2H2 + H2 */
    qf[92] = sc[2]*sc[24];
    qr[92] = 0.0;

    /*reaction 94: C2H3 + O2 => CH2O + HCO */
    qf[93] = sc[5]*sc[24];
    qr[93] = 0.0;

    /*reaction 95: C2H3 + O2 => C2H2 + HO2 */
    qf[94] = sc[5]*sc[24];
    qr[94] = 0.0;

    /*reaction 96: C2H4 + O => CH3 + HCO */
    qf[95] = sc[1]*sc[16];
    qr[95] = 0.0;

    /*reaction 97: C2H4 + OH => C2H3 + H2O */
    qf[96] = sc[3]*sc[16];
    qr[96] = 0.0;

    /*reaction 98: C2H3 + H2O => C2H4 + OH */
    qf[97] = sc[6]*sc[24];
    qr[97] = 0.0;

    /*reaction 99: C2H4 + H => C2H3 + H2 */
    qf[98] = sc[2]*sc[16];
    qr[98] = 0.0;

    /*reaction 100: C2H3 + H2 => C2H4 + H */
    qf[99] = sc[4]*sc[24];
    qr[99] = 0.0;

    /*reaction 101: C2H4 + CH3O2 => C2H3 + CH3O2H */
    qf[100] = sc[16]*sc[20];
    qr[100] = 0.0;

    /*reaction 102: C2H4 + CH3 => C2H3 + CH4 */
    qf[101] = sc[11]*sc[16];
    qr[101] = 0.0;

    /*reaction 103: C2H3 + CH4 => C2H4 + CH3 */
    qf[102] = sc[15]*sc[24];
    qr[102] = 0.0;

    /*reaction 104: C2H4 + O => CH2CHO + H */
    qf[103] = sc[1]*sc[16];
    qr[103] = 0.0;

    /*reaction 105: CH3 + C2H5 => CH4 + C2H4 */
    qf[104] = sc[11]*sc[14];
    qr[104] = 0.0;

    /*reaction 106: C2H5 + O2 => C2H4 + HO2 */
    qf[105] = sc[5]*sc[14];
    qr[105] = 0.0;

    /*reaction 107: C2H4 + HO2 => C2H5 + O2 */
    qf[106] = sc[8]*sc[16];
    qr[106] = 0.0;

    /*reaction 108: C2H5 + HO2 => C2H6 + O2 */
    qf[107] = sc[8]*sc[14];
    qr[107] = 0.0;

    /*reaction 109: H + C2H5 => C2H6 */
    qf[108] = sc[2]*sc[14];
    qr[108] = 0.0;

    /*reaction 110: C2H6 + H => C2H5 + H2 */
    qf[109] = sc[2]*sc[17];
    qr[109] = 0.0;

    /*reaction 111: C2H5 + H2 => C2H6 + H */
    qf[110] = sc[4]*sc[14];
    qr[110] = 0.0;

    /*reaction 112: C2H6 + OH => C2H5 + H2O */
    qf[111] = sc[3]*sc[17];
    qr[111] = 0.0;

    /*reaction 113: CH2GSG + C2H6 => CH3 + C2H5 */
    qf[112] = sc[9]*sc[17];
    qr[112] = 0.0;

    /*reaction 114: C2H6 + O => C2H5 + OH */
    qf[113] = sc[1]*sc[17];
    qr[113] = 0.0;

    /*reaction 115: C2H6 + CH3 => C2H5 + CH4 */
    qf[114] = sc[11]*sc[17];
    qr[114] = 0.0;

    /*reaction 116: HCCO + O => H + 2.000000 CO */
    qf[115] = sc[1]*sc[23];
    qr[115] = 0.0;

    /*reaction 117: HCCO + O2 => CO2 + HCO */
    qf[116] = sc[5]*sc[23];
    qr[116] = 0.0;

    /*reaction 118: HCCO + OH => 2.000000 HCO */
    qf[117] = sc[3]*sc[23];
    qr[117] = 0.0;

    /*reaction 119: HCCO + H => CH2GSG + CO */
    qf[118] = sc[2]*sc[23];
    qr[118] = 0.0;

    /*reaction 120: CH2GSG + CO => HCCO + H */
    qf[119] = sc[9]*sc[12];
    qr[119] = 0.0;

    /*reaction 121: CH2CHO + O2 => CH2O + CO + OH */
    qf[120] = sc[5]*sc[25];
    qr[120] = 0.0;

    /*reaction 122: C3H5XA + CH2O => C3H6 + HCO */
    qf[121] = sc[10]*sc[26];
    qr[121] = 0.0;

    /*reaction 123: C3H6 + HCO => C3H5XA + CH2O */
    qf[122] = sc[19]*sc[27];
    qr[122] = 0.0;

    /*reaction 124: C3H5XA + HO2 => C3H5O + OH */
    qf[123] = sc[8]*sc[26];
    qr[123] = 0.0;

    /*reaction 125: C3H5XA + CH3O2 => C3H5O + CH3O */
    qf[124] = sc[20]*sc[26];
    qr[124] = 0.0;

    /*reaction 126: C3H5XA + HO2 => C3H6 + O2 */
    qf[125] = sc[8]*sc[26];
    qr[125] = 0.0;

    /*reaction 127: C3H5XA + O2 => CH2CHO + CH2O */
    qf[126] = sc[5]*sc[26];
    qr[126] = 0.0;

    /*reaction 128: C3H5XA + O2 => C2H2 + CH2O + OH */
    qf[127] = sc[5]*sc[26];
    qr[127] = 0.0;

    /*reaction 129: C3H5XA + H => C3H6 */
    qf[128] = sc[2]*sc[26];
    qr[128] = 0.0;

    /*reaction 130: C3H5XA => C2H2 + CH3 */
    qf[129] = sc[26];
    qr[129] = 0.0;

    /*reaction 131: C3H6 => C2H3 + CH3 */
    qf[130] = sc[27];
    qr[130] = 0.0;

    /*reaction 132: C2H3 + CH3 => C3H6 */
    qf[131] = sc[11]*sc[24];
    qr[131] = 0.0;

    /*reaction 133: C3H6 + O => C2H5 + HCO */
    qf[132] = sc[1]*sc[27];
    qr[132] = 0.0;

    /*reaction 134: C3H6 + H => C2H4 + CH3 */
    qf[133] = sc[2]*sc[27];
    qr[133] = 0.0;

    /*reaction 135: C2H4 + CH3 => C3H6 + H */
    qf[134] = sc[11]*sc[16];
    qr[134] = 0.0;

    /*reaction 136: C3H6 + H => C3H5XA + H2 */
    qf[135] = sc[2]*sc[27];
    qr[135] = 0.0;

    /*reaction 137: C3H5XA + H2 => C3H6 + H */
    qf[136] = sc[4]*sc[26];
    qr[136] = 0.0;

    /*reaction 138: C3H6 + OH => C3H5XA + H2O */
    qf[137] = sc[3]*sc[27];
    qr[137] = 0.0;

    /*reaction 139: C3H6 + O => C3H5XA + OH */
    qf[138] = sc[1]*sc[27];
    qr[138] = 0.0;

    /*reaction 140: C3H6 + CH3 => C3H5XA + CH4 */
    qf[139] = sc[11]*sc[27];
    qr[139] = 0.0;

    /*reaction 141: IXC3H7 + O2 => C3H6 + HO2 */
    qf[140] = sc[5]*sc[29];
    qr[140] = 0.0;

    /*reaction 142: IXC3H7 => H + C3H6 */
    qf[141] = sc[29];
    qr[141] = 0.0;

    /*reaction 143: H + C3H6 => IXC3H7 */
    qf[142] = sc[2]*sc[27];
    qr[142] = 0.0;

    /*reaction 144: NXC3H7 => CH3 + C2H4 */
    qf[143] = sc[30];
    qr[143] = 0.0;

    /*reaction 145: CH3 + C2H4 => NXC3H7 */
    qf[144] = sc[11]*sc[16];
    qr[144] = 0.0;

    /*reaction 146: NXC3H7 + HO2 => C3H8 + O2 */
    qf[145] = sc[8]*sc[30];
    qr[145] = 0.0;

    /*reaction 147: NXC3H7 + O2 => C3H6 + HO2 */
    qf[146] = sc[5]*sc[30];
    qr[146] = 0.0;

    /*reaction 148: NXC3H7 => H + C3H6 */
    qf[147] = sc[30];
    qr[147] = 0.0;

    /*reaction 149: H + C3H6 => NXC3H7 */
    qf[148] = sc[2]*sc[27];
    qr[148] = 0.0;

    /*reaction 150: C3H8 + OH => NXC3H7 + H2O */
    qf[149] = sc[3]*sc[31];
    qr[149] = 0.0;

    /*reaction 151: C3H8 + HO2 => NXC3H7 + H2O2 */
    qf[150] = sc[8]*sc[31];
    qr[150] = 0.0;

    /*reaction 152: H + C3H8 <=> H2 + NXC3H7 */
    qf[151] = sc[2]*sc[31];
    qr[151] = sc[4]*sc[30];

    /*reaction 153: C3H8 + OH => IXC3H7 + H2O */
    qf[152] = sc[3]*sc[31];
    qr[152] = 0.0;

    /*reaction 154: CH3 + C3H8 => CH4 + IXC3H7 */
    qf[153] = sc[11]*sc[31];
    qr[153] = 0.0;

    /*reaction 155: CH3 + C3H8 => CH4 + NXC3H7 */
    qf[154] = sc[11]*sc[31];
    qr[154] = 0.0;

    /*reaction 156: C3H8 + O => IXC3H7 + OH */
    qf[155] = sc[1]*sc[31];
    qr[155] = 0.0;

    /*reaction 157: C3H8 + HO2 => IXC3H7 + H2O2 */
    qf[156] = sc[8]*sc[31];
    qr[156] = 0.0;

    /*reaction 158: C3H8 + O => NXC3H7 + OH */
    qf[157] = sc[1]*sc[31];
    qr[157] = 0.0;

    /*reaction 159: C3H8 + O2 => IXC3H7 + HO2 */
    qf[158] = sc[5]*sc[31];
    qr[158] = 0.0;

    /*reaction 160: IXC3H7 + HO2 => C3H8 + O2 */
    qf[159] = sc[8]*sc[29];
    qr[159] = 0.0;

    /*reaction 161: H + C3H8 => H2 + IXC3H7 */
    qf[160] = sc[2]*sc[31];
    qr[160] = 0.0;

    /*reaction 162: H2 + IXC3H7 => H + C3H8 */
    qf[161] = sc[4]*sc[29];
    qr[161] = 0.0;

    /*reaction 163: C3H5O => C2H3 + CH2O */
    qf[162] = sc[28];
    qr[162] = 0.0;

    /*reaction 164: IXC3H7O2 => IXC3H7 + O2 */
    qf[163] = sc[32];
    qr[163] = 0.0;

    /*reaction 165: IXC3H7 + O2 => IXC3H7O2 */
    qf[164] = sc[5]*sc[29];
    qr[164] = 0.0;

    /*reaction 166: NXC3H7O2 => NXC3H7 + O2 */
    qf[165] = sc[33];
    qr[165] = 0.0;

    /*reaction 167: NXC3H7 + O2 => NXC3H7O2 */
    qf[166] = sc[5]*sc[30];
    qr[166] = 0.0;

    amrex::Real T = tc[1];

    /*compute the mixture concentration */
    amrex::Real mixture = 0.0;
    for (int i = 0; i < 34; ++i) {
        mixture += sc[i];
    }

    amrex::Real Corr[167];
    for (int i = 0; i < 167; ++i) {
        Corr[i] = 1.0;
    }

    /* troe */
    {
        amrex::Real alpha[9];
        alpha[0] = mixture + (TB[0][0] - 1)*sc[4] + (TB[0][1] - 1)*sc[6] + (TB[0][2] - 1)*sc[18] + (TB[0][3] - 1)*sc[12];
        alpha[1] = mixture + (TB[1][0] - 1)*sc[4] + (TB[1][1] - 1)*sc[6] + (TB[1][2] - 1)*sc[18] + (TB[1][3] - 1)*sc[12];
        alpha[2] = mixture + (TB[2][0] - 1)*sc[4] + (TB[2][1] - 1)*sc[6] + (TB[2][2] - 1)*sc[18] + (TB[2][3] - 1)*sc[12];
        alpha[3] = mixture + (TB[3][0] - 1)*sc[4] + (TB[3][1] - 1)*sc[6] + (TB[3][2] - 1)*sc[18] + (TB[3][3] - 1)*sc[12];
        alpha[4] = mixture + (TB[4][0] - 1)*sc[4] + (TB[4][1] - 1)*sc[6] + (TB[4][2] - 1)*sc[18] + (TB[4][3] - 1)*sc[12];
        alpha[5] = mixture;
        alpha[6] = mixture + (TB[6][0] - 1)*sc[4] + (TB[6][1] - 1)*sc[6] + (TB[6][2] - 1)*sc[18] + (TB[6][3] - 1)*sc[12];
        alpha[7] = alpha[5];
        alpha[8] = alpha[5];
        for (int i=0; i<9; i++)
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
        alpha = mixture + (TB[9][0] - 1)*sc[4] + (TB[9][1] - 1)*sc[6] + (TB[9][2] - 1)*sc[18] + (TB[9][3] - 1)*sc[12];
        Corr[9] = alpha;
        alpha = mixture + (TB[10][0] - 1)*sc[4] + (TB[10][1] - 1)*sc[6] + (TB[10][2] - 1)*sc[18] + (TB[10][3] - 1)*sc[12];
        Corr[10] = alpha;
        alpha = mixture + (TB[11][0] - 1)*sc[4] + (TB[11][1] - 1)*sc[6] + (TB[11][2] - 1)*sc[18] + (TB[11][3] - 1)*sc[12];
        Corr[11] = alpha;
        alpha = mixture + (TB[12][0] - 1)*sc[4] + (TB[12][1] - 1)*sc[6] + (TB[12][2] - 1)*sc[18] + (TB[12][3] - 1)*sc[12];
        Corr[12] = alpha;
        alpha = mixture + (TB[13][0] - 1)*sc[4] + (TB[13][1] - 1)*sc[6] + (TB[13][2] - 1)*sc[18] + (TB[13][3] - 1)*sc[12];
        Corr[13] = alpha;
        alpha = mixture;
        Corr[14] = alpha;
        Corr[15] = alpha;
    }

    for (int i=0; i<167; i++)
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

    amrex::Real q_f[167], q_r[167];
    amrex::Real sc_qss[1];
    comp_qfqr(q_f, q_r, sc, sc_qss, tc, invT);

    for (int i = 0; i < 167; ++i) {
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
    amrex::Real c[34]; /*temporary storage */
    amrex::Real PORT = 1e6 * (*P)/(8.31446261815324e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 34; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRateFR(q_f, q_r, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 167; ++id) {
        q_f[id] *= 1.0e-6;
        q_r[id] *= 1.0e-6;
    }
}


/*Returns the rate of progress for each reaction */
void CKQC(amrex::Real *  T, amrex::Real *  C, amrex::Real *  qdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 34; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    progressRate(qdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 34; ++id) {
        C[id] *= 1.0e-6;
    }

    for (id = 0; id < 167; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mass fractions */
void CKQYP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y, amrex::Real *  qdot)
{
    int id; /*loop counter */
    amrex::Real c[34]; /*temporary storage */
    amrex::Real YOW = 0; 
    amrex::Real PWORT; 
    amrex::Real imw[34];

    get_imw(imw);

    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*N2 */
    YOW += y[1]*imw[1]; /*O */
    YOW += y[2]*imw[2]; /*H */
    YOW += y[3]*imw[3]; /*OH */
    YOW += y[4]*imw[4]; /*H2 */
    YOW += y[5]*imw[5]; /*O2 */
    YOW += y[6]*imw[6]; /*H2O */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*HO2 */
    YOW += y[9]*imw[9]; /*CH2GSG */
    YOW += y[10]*imw[10]; /*CH2O */
    YOW += y[11]*imw[11]; /*CH3 */
    YOW += y[12]*imw[12]; /*CO */
    YOW += y[13]*imw[13]; /*CH3O */
    YOW += y[14]*imw[14]; /*C2H5 */
    YOW += y[15]*imw[15]; /*CH4 */
    YOW += y[16]*imw[16]; /*C2H4 */
    YOW += y[17]*imw[17]; /*C2H6 */
    YOW += y[18]*imw[18]; /*CO2 */
    YOW += y[19]*imw[19]; /*HCO */
    YOW += y[20]*imw[20]; /*CH3O2 */
    YOW += y[21]*imw[21]; /*CH3O2H */
    YOW += y[22]*imw[22]; /*C2H2 */
    YOW += y[23]*imw[23]; /*HCCO */
    YOW += y[24]*imw[24]; /*C2H3 */
    YOW += y[25]*imw[25]; /*CH2CHO */
    YOW += y[26]*imw[26]; /*C3H5XA */
    YOW += y[27]*imw[27]; /*C3H6 */
    YOW += y[28]*imw[28]; /*C3H5O */
    YOW += y[29]*imw[29]; /*IXC3H7 */
    YOW += y[30]*imw[30]; /*NXC3H7 */
    YOW += y[31]*imw[31]; /*C3H8 */
    YOW += y[32]*imw[32]; /*IXC3H7O2 */
    YOW += y[33]*imw[33]; /*NXC3H7O2 */
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

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 167; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKQXP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x, amrex::Real *  qdot)
{
    int id; /*loop counter */
    amrex::Real c[34]; /*temporary storage */
    amrex::Real PORT = 1e6 * (*P)/(8.31446261815324e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 34; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 167; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mass fractions */
void CKQYR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y, amrex::Real *  qdot)
{
    int id; /*loop counter */
    amrex::Real c[34]; /*temporary storage */
    amrex::Real imw[34];

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

    /*call progressRate */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 167; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mole fractions */
void CKQXR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  x, amrex::Real *  qdot)
{
    int id; /*loop counter */
    amrex::Real c[34]; /*temporary storage */
    amrex::Real XW = 0; /*See Eq 4, 11 in CK Manual */
    amrex::Real ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*28.013400; /*N2 */
    XW += x[1]*15.999400; /*O */
    XW += x[2]*1.007970; /*H */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*2.015940; /*H2 */
    XW += x[5]*31.998800; /*O2 */
    XW += x[6]*18.015340; /*H2O */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*33.006770; /*HO2 */
    XW += x[9]*14.027090; /*CH2GSG */
    XW += x[10]*30.026490; /*CH2O */
    XW += x[11]*15.035060; /*CH3 */
    XW += x[12]*28.010550; /*CO */
    XW += x[13]*31.034460; /*CH3O */
    XW += x[14]*29.062150; /*C2H5 */
    XW += x[15]*16.043030; /*CH4 */
    XW += x[16]*28.054180; /*C2H4 */
    XW += x[17]*30.070120; /*C2H6 */
    XW += x[18]*44.009950; /*CO2 */
    XW += x[19]*29.018520; /*HCO */
    XW += x[20]*47.033860; /*CH3O2 */
    XW += x[21]*48.041830; /*CH3O2H */
    XW += x[22]*26.038240; /*C2H2 */
    XW += x[23]*41.029670; /*HCCO */
    XW += x[24]*27.046210; /*C2H3 */
    XW += x[25]*43.045610; /*CH2CHO */
    XW += x[26]*41.073300; /*C3H5XA */
    XW += x[27]*42.081270; /*C3H6 */
    XW += x[28]*57.072700; /*C3H5O */
    XW += x[29]*43.089240; /*IXC3H7 */
    XW += x[30]*43.089240; /*NXC3H7 */
    XW += x[31]*44.097210; /*C3H8 */
    XW += x[32]*75.088040; /*IXC3H7O2 */
    XW += x[33]*75.088040; /*NXC3H7O2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 34; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 167; ++id) {
        qdot[id] *= 1.0e-6;
    }
}

/*compute the reaction Jacobian on CPU */
void aJacobian(amrex::Real *  J, amrex::Real *  sc, amrex::Real T, int consP)
{
    for (int i=0; i<1225; i++) {
        J[i] = 0.0;
    }

    amrex::Real wdot[34];
    for (int k=0; k<34; k++) {
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
    for (int k = 0; k < 34; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    amrex::Real g_RT[34];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    amrex::Real h_RT[34];
    speciesEnthalpy(h_RT, tc);

    amrex::Real phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    amrex::Real dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    amrex::Real dqdci, dcdc_fac, dqdc[34];
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
    alpha = mixture + (TB[0][0] - 1)*sc[4] + (TB[0][1] - 1)*sc[6] + (TB[0][2] - 1)*sc[18] + (TB[0][3] - 1)*sc[12];
    /* forward */
    phi_f = pow(sc[3], 2.000000);
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
    Kc = refCinv * exp(2.000000*g_RT[3] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[3]) + (h_RT[7]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[3] -= 2 * q; /* OH */
    wdot[7] += q; /* H2O2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[OH] */
        dqdci =  + k_f*2.000000*sc[3];
        J[108] += -2 * dqdci;         /* dwdot[OH]/d[OH] */
        J[112] += dqdci;              /* dwdot[H2O2]/d[OH] */
        /* d()/d[H2] */
        dqdci = (TB[0][0] - 1)*dcdc_fac;
        J[143] += -2 * dqdci;         /* dwdot[OH]/d[H2] */
        J[147] += dqdci;              /* dwdot[H2O2]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[0][1] - 1)*dcdc_fac;
        J[213] += -2 * dqdci;         /* dwdot[OH]/d[H2O] */
        J[217] += dqdci;              /* dwdot[H2O2]/d[H2O] */
        /* d()/d[H2O2] */
        dqdci =  - k_r;
        J[248] += -2 * dqdci;         /* dwdot[OH]/d[H2O2] */
        J[252] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
        /* d()/d[CO] */
        dqdci = (TB[0][3] - 1)*dcdc_fac;
        J[423] += -2 * dqdci;         /* dwdot[OH]/d[CO] */
        J[427] += dqdci;              /* dwdot[H2O2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[0][2] - 1)*dcdc_fac;
        J[633] += -2 * dqdci;         /* dwdot[OH]/d[CO2] */
        J[637] += dqdci;              /* dwdot[H2O2]/d[CO2] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac + k_f*2.000000*sc[3];
        dqdc[4] = TB[0][0]*dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = TB[0][1]*dcdc_fac;
        dqdc[7] = dcdc_fac - k_r;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = TB[0][3]*dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = TB[0][2]*dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac;
        dqdc[22] = dcdc_fac;
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
        for (int k=0; k<34; k++) {
            J[35*k+3] += -2 * dqdc[k];
            J[35*k+7] += dqdc[k];
        }
    }
    J[1193] += -2 * dqdT; /* dwdot[OH]/dT */
    J[1197] += dqdT; /* dwdot[H2O2]/dT */

    /*reaction 2: H + O2 (+M) <=> HO2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[1][0] - 1)*sc[4] + (TB[1][1] - 1)*sc[6] + (TB[1][2] - 1)*sc[18] + (TB[1][3] - 1)*sc[12];
    /* forward */
    phi_f = sc[2]*sc[5];
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
    Kc = refCinv * exp(g_RT[2] + g_RT[5] - g_RT[8]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[5]) + (h_RT[8]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[5] -= q; /* O2 */
    wdot[8] += q; /* HO2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[5];
        J[72] -= dqdci;               /* dwdot[H]/d[H] */
        J[75] -= dqdci;               /* dwdot[O2]/d[H] */
        J[78] += dqdci;               /* dwdot[HO2]/d[H] */
        /* d()/d[H2] */
        dqdci = (TB[1][0] - 1)*dcdc_fac;
        J[142] -= dqdci;              /* dwdot[H]/d[H2] */
        J[145] -= dqdci;              /* dwdot[O2]/d[H2] */
        J[148] += dqdci;              /* dwdot[HO2]/d[H2] */
        /* d()/d[O2] */
        dqdci =  + k_f*sc[2];
        J[177] -= dqdci;              /* dwdot[H]/d[O2] */
        J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
        J[183] += dqdci;              /* dwdot[HO2]/d[O2] */
        /* d()/d[H2O] */
        dqdci = (TB[1][1] - 1)*dcdc_fac;
        J[212] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[215] -= dqdci;              /* dwdot[O2]/d[H2O] */
        J[218] += dqdci;              /* dwdot[HO2]/d[H2O] */
        /* d()/d[HO2] */
        dqdci =  - k_r;
        J[282] -= dqdci;              /* dwdot[H]/d[HO2] */
        J[285] -= dqdci;              /* dwdot[O2]/d[HO2] */
        J[288] += dqdci;              /* dwdot[HO2]/d[HO2] */
        /* d()/d[CO] */
        dqdci = (TB[1][3] - 1)*dcdc_fac;
        J[422] -= dqdci;              /* dwdot[H]/d[CO] */
        J[425] -= dqdci;              /* dwdot[O2]/d[CO] */
        J[428] += dqdci;              /* dwdot[HO2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[1][2] - 1)*dcdc_fac;
        J[632] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[635] -= dqdci;              /* dwdot[O2]/d[CO2] */
        J[638] += dqdci;              /* dwdot[HO2]/d[CO2] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac + k_f*sc[5];
        dqdc[3] = dcdc_fac;
        dqdc[4] = TB[1][0]*dcdc_fac;
        dqdc[5] = dcdc_fac + k_f*sc[2];
        dqdc[6] = TB[1][1]*dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac - k_r;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = TB[1][3]*dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = TB[1][2]*dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac;
        dqdc[22] = dcdc_fac;
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
        for (int k=0; k<34; k++) {
            J[35*k+2] -= dqdc[k];
            J[35*k+5] -= dqdc[k];
            J[35*k+8] += dqdc[k];
        }
    }
    J[1192] -= dqdT; /* dwdot[H]/dT */
    J[1195] -= dqdT; /* dwdot[O2]/dT */
    J[1198] += dqdT; /* dwdot[HO2]/dT */

    /*reaction 3: CH3 + H (+M) => CH4 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* FIXME: inreversible reaction in _ajac_reaction may not work*/
    /* 3-body correction factor */
    alpha = mixture + (TB[2][0] - 1)*sc[4] + (TB[2][1] - 1)*sc[6] + (TB[2][2] - 1)*sc[18] + (TB[2][3] - 1)*sc[12];
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
    /* rate of progress */
    q_nocor = k_f*phi_f;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *dlnkfdT*k_f*phi_f + dlnCorrdT*q;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[11] -= q; /* CH3 */
    wdot[15] += q; /* CH4 */
    /* for convenience */
    k_f *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[11];
        J[72] -= dqdci;               /* dwdot[H]/d[H] */
        J[81] -= dqdci;               /* dwdot[CH3]/d[H] */
        J[85] += dqdci;               /* dwdot[CH4]/d[H] */
        /* d()/d[H2] */
        dqdci = (TB[2][0] - 1)*dcdc_fac;
        J[142] -= dqdci;              /* dwdot[H]/d[H2] */
        J[151] -= dqdci;              /* dwdot[CH3]/d[H2] */
        J[155] += dqdci;              /* dwdot[CH4]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[2][1] - 1)*dcdc_fac;
        J[212] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[221] -= dqdci;              /* dwdot[CH3]/d[H2O] */
        J[225] += dqdci;              /* dwdot[CH4]/d[H2O] */
        /* d()/d[CH3] */
        dqdci =  + k_f*sc[2];
        J[387] -= dqdci;              /* dwdot[H]/d[CH3] */
        J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
        J[400] += dqdci;              /* dwdot[CH4]/d[CH3] */
        /* d()/d[CO] */
        dqdci = (TB[2][3] - 1)*dcdc_fac;
        J[422] -= dqdci;              /* dwdot[H]/d[CO] */
        J[431] -= dqdci;              /* dwdot[CH3]/d[CO] */
        J[435] += dqdci;              /* dwdot[CH4]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[2][2] - 1)*dcdc_fac;
        J[632] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[641] -= dqdci;              /* dwdot[CH3]/d[CO2] */
        J[645] += dqdci;              /* dwdot[CH4]/d[CO2] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac + k_f*sc[11];
        dqdc[3] = dcdc_fac;
        dqdc[4] = TB[2][0]*dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = TB[2][1]*dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac + k_f*sc[2];
        dqdc[12] = TB[2][3]*dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = TB[2][2]*dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac;
        dqdc[22] = dcdc_fac;
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
        for (int k=0; k<34; k++) {
            J[35*k+2] -= dqdc[k];
            J[35*k+11] -= dqdc[k];
            J[35*k+15] += dqdc[k];
        }
    }
    J[1192] -= dqdT; /* dwdot[H]/dT */
    J[1201] -= dqdT; /* dwdot[CH3]/dT */
    J[1205] += dqdT; /* dwdot[CH4]/dT */

    /*reaction 4: 2.000000 CH3 (+M) <=> C2H6 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[3][0] - 1)*sc[4] + (TB[3][1] - 1)*sc[6] + (TB[3][2] - 1)*sc[18] + (TB[3][3] - 1)*sc[12];
    /* forward */
    phi_f = pow(sc[11], 2.000000);
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
    phi_r = sc[17];
    Kc = refCinv * exp(2.000000*g_RT[11] - g_RT[17]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[11]) + (h_RT[17]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[11] -= 2 * q; /* CH3 */
    wdot[17] += q; /* C2H6 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[3][0] - 1)*dcdc_fac;
        J[151] += -2 * dqdci;         /* dwdot[CH3]/d[H2] */
        J[157] += dqdci;              /* dwdot[C2H6]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[3][1] - 1)*dcdc_fac;
        J[221] += -2 * dqdci;         /* dwdot[CH3]/d[H2O] */
        J[227] += dqdci;              /* dwdot[C2H6]/d[H2O] */
        /* d()/d[CH3] */
        dqdci =  + k_f*2.000000*sc[11];
        J[396] += -2 * dqdci;         /* dwdot[CH3]/d[CH3] */
        J[402] += dqdci;              /* dwdot[C2H6]/d[CH3] */
        /* d()/d[CO] */
        dqdci = (TB[3][3] - 1)*dcdc_fac;
        J[431] += -2 * dqdci;         /* dwdot[CH3]/d[CO] */
        J[437] += dqdci;              /* dwdot[C2H6]/d[CO] */
        /* d()/d[C2H6] */
        dqdci =  - k_r;
        J[606] += -2 * dqdci;         /* dwdot[CH3]/d[C2H6] */
        J[612] += dqdci;              /* dwdot[C2H6]/d[C2H6] */
        /* d()/d[CO2] */
        dqdci = (TB[3][2] - 1)*dcdc_fac;
        J[641] += -2 * dqdci;         /* dwdot[CH3]/d[CO2] */
        J[647] += dqdci;              /* dwdot[C2H6]/d[CO2] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = TB[3][0]*dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = TB[3][1]*dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac + k_f*2.000000*sc[11];
        dqdc[12] = TB[3][3]*dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac - k_r;
        dqdc[18] = TB[3][2]*dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac;
        dqdc[22] = dcdc_fac;
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
        for (int k=0; k<34; k++) {
            J[35*k+11] += -2 * dqdc[k];
            J[35*k+17] += dqdc[k];
        }
    }
    J[1201] += -2 * dqdT; /* dwdot[CH3]/dT */
    J[1207] += dqdT; /* dwdot[C2H6]/dT */

    /*reaction 5: CO + O (+M) => CO2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* FIXME: inreversible reaction in _ajac_reaction may not work*/
    /* 3-body correction factor */
    alpha = mixture + (TB[4][0] - 1)*sc[4] + (TB[4][1] - 1)*sc[6] + (TB[4][2] - 1)*sc[18] + (TB[4][3] - 1)*sc[12];
    /* forward */
    phi_f = sc[1]*sc[12];
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
    /* rate of progress */
    q_nocor = k_f*phi_f;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *dlnkfdT*k_f*phi_f + dlnCorrdT*q;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[12] -= q; /* CO */
    wdot[18] += q; /* CO2 */
    /* for convenience */
    k_f *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[O] */
        dqdci =  + k_f*sc[12];
        J[36] -= dqdci;               /* dwdot[O]/d[O] */
        J[47] -= dqdci;               /* dwdot[CO]/d[O] */
        J[53] += dqdci;               /* dwdot[CO2]/d[O] */
        /* d()/d[H2] */
        dqdci = (TB[4][0] - 1)*dcdc_fac;
        J[141] -= dqdci;              /* dwdot[O]/d[H2] */
        J[152] -= dqdci;              /* dwdot[CO]/d[H2] */
        J[158] += dqdci;              /* dwdot[CO2]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[4][1] - 1)*dcdc_fac;
        J[211] -= dqdci;              /* dwdot[O]/d[H2O] */
        J[222] -= dqdci;              /* dwdot[CO]/d[H2O] */
        J[228] += dqdci;              /* dwdot[CO2]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[4][3] - 1)*dcdc_fac + k_f*sc[1];
        J[421] -= dqdci;              /* dwdot[O]/d[CO] */
        J[432] -= dqdci;              /* dwdot[CO]/d[CO] */
        J[438] += dqdci;              /* dwdot[CO2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[4][2] - 1)*dcdc_fac;
        J[631] -= dqdci;              /* dwdot[O]/d[CO2] */
        J[642] -= dqdci;              /* dwdot[CO]/d[CO2] */
        J[648] += dqdci;              /* dwdot[CO2]/d[CO2] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[12];
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = TB[4][0]*dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = TB[4][1]*dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = TB[4][3]*dcdc_fac + k_f*sc[1];
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = TB[4][2]*dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac;
        dqdc[22] = dcdc_fac;
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
        for (int k=0; k<34; k++) {
            J[35*k+1] -= dqdc[k];
            J[35*k+12] -= dqdc[k];
            J[35*k+18] += dqdc[k];
        }
    }
    J[1191] -= dqdT; /* dwdot[O]/dT */
    J[1202] -= dqdT; /* dwdot[CO]/dT */
    J[1208] += dqdT; /* dwdot[CO2]/dT */

    /*reaction 6: CH3O (+M) => CH2O + H (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* FIXME: inreversible reaction in _ajac_reaction may not work*/
    /* 3-body correction factor */
    alpha = mixture;
    /* forward */
    phi_f = sc[13];
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
    /* rate of progress */
    q_nocor = k_f*phi_f;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *dlnkfdT*k_f*phi_f + dlnCorrdT*q;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[10] += q; /* CH2O */
    wdot[13] -= q; /* CH3O */
    /* for convenience */
    k_f *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[CH3O] */
        dqdci =  + k_f;
        J[457] += dqdci;              /* dwdot[H]/d[CH3O] */
        J[465] += dqdci;              /* dwdot[CH2O]/d[CH3O] */
        J[468] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = dcdc_fac + k_f;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
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
        dqdc[26] = dcdc_fac;
        dqdc[27] = dcdc_fac;
        dqdc[28] = dcdc_fac;
        dqdc[29] = dcdc_fac;
        dqdc[30] = dcdc_fac;
        dqdc[31] = dcdc_fac;
        dqdc[32] = dcdc_fac;
        dqdc[33] = dcdc_fac;
        for (int k=0; k<34; k++) {
            J[35*k+2] += dqdc[k];
            J[35*k+10] += dqdc[k];
            J[35*k+13] -= dqdc[k];
        }
    }
    J[1192] += dqdT; /* dwdot[H]/dT */
    J[1200] += dqdT; /* dwdot[CH2O]/dT */
    J[1203] -= dqdT; /* dwdot[CH3O]/dT */

    /*reaction 7: C2H3 (+M) => H + C2H2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* FIXME: inreversible reaction in _ajac_reaction may not work*/
    /* 3-body correction factor */
    alpha = mixture + (TB[6][0] - 1)*sc[4] + (TB[6][1] - 1)*sc[6] + (TB[6][2] - 1)*sc[18] + (TB[6][3] - 1)*sc[12];
    /* forward */
    phi_f = sc[24];
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
    /* rate of progress */
    q_nocor = k_f*phi_f;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *dlnkfdT*k_f*phi_f + dlnCorrdT*q;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[22] += q; /* C2H2 */
    wdot[24] -= q; /* C2H3 */
    /* for convenience */
    k_f *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[6][0] - 1)*dcdc_fac;
        J[142] += dqdci;              /* dwdot[H]/d[H2] */
        J[162] += dqdci;              /* dwdot[C2H2]/d[H2] */
        J[164] -= dqdci;              /* dwdot[C2H3]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[6][1] - 1)*dcdc_fac;
        J[212] += dqdci;              /* dwdot[H]/d[H2O] */
        J[232] += dqdci;              /* dwdot[C2H2]/d[H2O] */
        J[234] -= dqdci;              /* dwdot[C2H3]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[6][3] - 1)*dcdc_fac;
        J[422] += dqdci;              /* dwdot[H]/d[CO] */
        J[442] += dqdci;              /* dwdot[C2H2]/d[CO] */
        J[444] -= dqdci;              /* dwdot[C2H3]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[6][2] - 1)*dcdc_fac;
        J[632] += dqdci;              /* dwdot[H]/d[CO2] */
        J[652] += dqdci;              /* dwdot[C2H2]/d[CO2] */
        J[654] -= dqdci;              /* dwdot[C2H3]/d[CO2] */
        /* d()/d[C2H3] */
        dqdci =  + k_f;
        J[842] += dqdci;              /* dwdot[H]/d[C2H3] */
        J[862] += dqdci;              /* dwdot[C2H2]/d[C2H3] */
        J[864] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = TB[6][0]*dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = TB[6][1]*dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = TB[6][3]*dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac;
        dqdc[17] = dcdc_fac;
        dqdc[18] = TB[6][2]*dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac;
        dqdc[22] = dcdc_fac;
        dqdc[23] = dcdc_fac;
        dqdc[24] = dcdc_fac + k_f;
        dqdc[25] = dcdc_fac;
        dqdc[26] = dcdc_fac;
        dqdc[27] = dcdc_fac;
        dqdc[28] = dcdc_fac;
        dqdc[29] = dcdc_fac;
        dqdc[30] = dcdc_fac;
        dqdc[31] = dcdc_fac;
        dqdc[32] = dcdc_fac;
        dqdc[33] = dcdc_fac;
        for (int k=0; k<34; k++) {
            J[35*k+2] += dqdc[k];
            J[35*k+22] += dqdc[k];
            J[35*k+24] -= dqdc[k];
        }
    }
    J[1192] += dqdT; /* dwdot[H]/dT */
    J[1212] += dqdT; /* dwdot[C2H2]/dT */
    J[1214] -= dqdT; /* dwdot[C2H3]/dT */

    /*reaction 8: H + C2H4 (+M) <=> C2H5 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture;
    /* forward */
    phi_f = sc[2]*sc[16];
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
    phi_r = sc[14];
    Kc = refCinv * exp(g_RT[2] - g_RT[14] + g_RT[16]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[16]) + (h_RT[14]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[14] += q; /* C2H5 */
    wdot[16] -= q; /* C2H4 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[16];
        J[72] -= dqdci;               /* dwdot[H]/d[H] */
        J[84] += dqdci;               /* dwdot[C2H5]/d[H] */
        J[86] -= dqdci;               /* dwdot[C2H4]/d[H] */
        /* d()/d[C2H5] */
        dqdci =  - k_r;
        J[492] -= dqdci;              /* dwdot[H]/d[C2H5] */
        J[504] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
        J[506] -= dqdci;              /* dwdot[C2H4]/d[C2H5] */
        /* d()/d[C2H4] */
        dqdci =  + k_f*sc[2];
        J[562] -= dqdci;              /* dwdot[H]/d[C2H4] */
        J[574] += dqdci;              /* dwdot[C2H5]/d[C2H4] */
        J[576] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac + k_f*sc[16];
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac - k_r;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac + k_f*sc[2];
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac;
        dqdc[22] = dcdc_fac;
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
        for (int k=0; k<34; k++) {
            J[35*k+2] -= dqdc[k];
            J[35*k+14] += dqdc[k];
            J[35*k+16] -= dqdc[k];
        }
    }
    J[1192] -= dqdT; /* dwdot[H]/dT */
    J[1204] += dqdT; /* dwdot[C2H5]/dT */
    J[1206] -= dqdT; /* dwdot[C2H4]/dT */

    /*reaction 9: C2H4 (+M) => C2H2 + H2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* FIXME: inreversible reaction in _ajac_reaction may not work*/
    /* 3-body correction factor */
    alpha = mixture;
    /* forward */
    phi_f = sc[16];
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
    /* rate of progress */
    q_nocor = k_f*phi_f;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *dlnkfdT*k_f*phi_f + dlnCorrdT*q;
    /* update wdot */
    wdot[4] += q; /* H2 */
    wdot[16] -= q; /* C2H4 */
    wdot[22] += q; /* C2H2 */
    /* for convenience */
    k_f *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[C2H4] */
        dqdci =  + k_f;
        J[564] += dqdci;              /* dwdot[H2]/d[C2H4] */
        J[576] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
        J[582] += dqdci;              /* dwdot[C2H2]/d[C2H4] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = dcdc_fac;
        dqdc[14] = dcdc_fac;
        dqdc[15] = dcdc_fac;
        dqdc[16] = dcdc_fac + k_f;
        dqdc[17] = dcdc_fac;
        dqdc[18] = dcdc_fac;
        dqdc[19] = dcdc_fac;
        dqdc[20] = dcdc_fac;
        dqdc[21] = dcdc_fac;
        dqdc[22] = dcdc_fac;
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
        for (int k=0; k<34; k++) {
            J[35*k+4] += dqdc[k];
            J[35*k+16] -= dqdc[k];
            J[35*k+22] += dqdc[k];
        }
    }
    J[1194] += dqdT; /* dwdot[H2]/dT */
    J[1206] -= dqdT; /* dwdot[C2H4]/dT */
    J[1212] += dqdT; /* dwdot[C2H2]/dT */

    /*reaction 10: O + H + M => OH + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* FIXME: inreversible reaction in _ajac_reaction may not work*/
    /* 3-body correction factor */
    alpha = mixture + (TB[9][0] - 1)*sc[4] + (TB[9][1] - 1)*sc[6] + (TB[9][2] - 1)*sc[18] + (TB[9][3] - 1)*sc[12];
    /* forward */
    phi_f = sc[1]*sc[2];
    k_f = prefactor_units[9] * fwd_A[9]
                * exp(fwd_beta[9] * tc[0] - activation_units[9] * fwd_Ea[9] * invT);
    dlnkfdT = fwd_beta[9] * invT + activation_units[9] * fwd_Ea[9] * invT2;
    /* rate of progress */
    q_nocor = k_f*phi_f;
    q = alpha * q_nocor;
    dqdT = alpha * dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[2] -= q; /* H */
    wdot[3] += q; /* OH */
    /* for convenience */
    k_f *= alpha;
    k_r = 0.0;
    if (consP) {
        /* d()/d[O] */
        dqdci =  + k_f*sc[2];
        J[36] -= dqdci;               /* dwdot[O]/d[O] */
        J[37] -= dqdci;               /* dwdot[H]/d[O] */
        J[38] += dqdci;               /* dwdot[OH]/d[O] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[1];
        J[71] -= dqdci;               /* dwdot[O]/d[H] */
        J[72] -= dqdci;               /* dwdot[H]/d[H] */
        J[73] += dqdci;               /* dwdot[OH]/d[H] */
        /* d()/d[H2] */
        dqdci = (TB[9][0] - 1)*q_nocor;
        J[141] -= dqdci;              /* dwdot[O]/d[H2] */
        J[142] -= dqdci;              /* dwdot[H]/d[H2] */
        J[143] += dqdci;              /* dwdot[OH]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[9][1] - 1)*q_nocor;
        J[211] -= dqdci;              /* dwdot[O]/d[H2O] */
        J[212] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[213] += dqdci;              /* dwdot[OH]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[9][3] - 1)*q_nocor;
        J[421] -= dqdci;              /* dwdot[O]/d[CO] */
        J[422] -= dqdci;              /* dwdot[H]/d[CO] */
        J[423] += dqdci;              /* dwdot[OH]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[9][2] - 1)*q_nocor;
        J[631] -= dqdci;              /* dwdot[O]/d[CO2] */
        J[632] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[633] += dqdci;              /* dwdot[OH]/d[CO2] */
    }
    else {
        dqdc[0] = q_nocor;
        dqdc[1] = q_nocor + k_f*sc[2];
        dqdc[2] = q_nocor + k_f*sc[1];
        dqdc[3] = q_nocor;
        dqdc[4] = TB[9][0]*q_nocor;
        dqdc[5] = q_nocor;
        dqdc[6] = TB[9][1]*q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = q_nocor;
        dqdc[12] = TB[9][3]*q_nocor;
        dqdc[13] = q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = TB[9][2]*q_nocor;
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
        for (int k=0; k<34; k++) {
            J[35*k+1] -= dqdc[k];
            J[35*k+2] -= dqdc[k];
            J[35*k+3] += dqdc[k];
        }
    }
    J[1191] -= dqdT; /* dwdot[O]/dT */
    J[1192] -= dqdT; /* dwdot[H]/dT */
    J[1193] += dqdT; /* dwdot[OH]/dT */

    /*reaction 11: 2.000000 O + M => O2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* FIXME: inreversible reaction in _ajac_reaction may not work*/
    /* 3-body correction factor */
    alpha = mixture + (TB[10][0] - 1)*sc[4] + (TB[10][1] - 1)*sc[6] + (TB[10][2] - 1)*sc[18] + (TB[10][3] - 1)*sc[12];
    /* forward */
    phi_f = pow(sc[1], 2.000000);
    k_f = prefactor_units[10] * fwd_A[10]
                * exp(fwd_beta[10] * tc[0] - activation_units[10] * fwd_Ea[10] * invT);
    dlnkfdT = fwd_beta[10] * invT + activation_units[10] * fwd_Ea[10] * invT2;
    /* rate of progress */
    q_nocor = k_f*phi_f;
    q = alpha * q_nocor;
    dqdT = alpha * dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= 2 * q; /* O */
    wdot[5] += q; /* O2 */
    /* for convenience */
    k_f *= alpha;
    k_r = 0.0;
    if (consP) {
        /* d()/d[O] */
        dqdci =  + k_f*2.000000*sc[1];
        J[36] += -2 * dqdci;          /* dwdot[O]/d[O] */
        J[40] += dqdci;               /* dwdot[O2]/d[O] */
        /* d()/d[H2] */
        dqdci = (TB[10][0] - 1)*q_nocor;
        J[141] += -2 * dqdci;         /* dwdot[O]/d[H2] */
        J[145] += dqdci;              /* dwdot[O2]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[10][1] - 1)*q_nocor;
        J[211] += -2 * dqdci;         /* dwdot[O]/d[H2O] */
        J[215] += dqdci;              /* dwdot[O2]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[10][3] - 1)*q_nocor;
        J[421] += -2 * dqdci;         /* dwdot[O]/d[CO] */
        J[425] += dqdci;              /* dwdot[O2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[10][2] - 1)*q_nocor;
        J[631] += -2 * dqdci;         /* dwdot[O]/d[CO2] */
        J[635] += dqdci;              /* dwdot[O2]/d[CO2] */
    }
    else {
        dqdc[0] = q_nocor;
        dqdc[1] = q_nocor + k_f*2.000000*sc[1];
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = TB[10][0]*q_nocor;
        dqdc[5] = q_nocor;
        dqdc[6] = TB[10][1]*q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = q_nocor;
        dqdc[12] = TB[10][3]*q_nocor;
        dqdc[13] = q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = TB[10][2]*q_nocor;
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
        for (int k=0; k<34; k++) {
            J[35*k+1] += -2 * dqdc[k];
            J[35*k+5] += dqdc[k];
        }
    }
    J[1191] += -2 * dqdT; /* dwdot[O]/dT */
    J[1195] += dqdT; /* dwdot[O2]/dT */

    /*reaction 12: H + OH + M => H2O + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* FIXME: inreversible reaction in _ajac_reaction may not work*/
    /* 3-body correction factor */
    alpha = mixture + (TB[11][0] - 1)*sc[4] + (TB[11][1] - 1)*sc[6] + (TB[11][2] - 1)*sc[18] + (TB[11][3] - 1)*sc[12];
    /* forward */
    phi_f = sc[2]*sc[3];
    k_f = prefactor_units[11] * fwd_A[11]
                * exp(fwd_beta[11] * tc[0] - activation_units[11] * fwd_Ea[11] * invT);
    dlnkfdT = fwd_beta[11] * invT + activation_units[11] * fwd_Ea[11] * invT2;
    /* rate of progress */
    q_nocor = k_f*phi_f;
    q = alpha * q_nocor;
    dqdT = alpha * dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[3] -= q; /* OH */
    wdot[6] += q; /* H2O */
    /* for convenience */
    k_f *= alpha;
    k_r = 0.0;
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[3];
        J[72] -= dqdci;               /* dwdot[H]/d[H] */
        J[73] -= dqdci;               /* dwdot[OH]/d[H] */
        J[76] += dqdci;               /* dwdot[H2O]/d[H] */
        /* d()/d[OH] */
        dqdci =  + k_f*sc[2];
        J[107] -= dqdci;              /* dwdot[H]/d[OH] */
        J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
        J[111] += dqdci;              /* dwdot[H2O]/d[OH] */
        /* d()/d[H2] */
        dqdci = (TB[11][0] - 1)*q_nocor;
        J[142] -= dqdci;              /* dwdot[H]/d[H2] */
        J[143] -= dqdci;              /* dwdot[OH]/d[H2] */
        J[146] += dqdci;              /* dwdot[H2O]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[11][1] - 1)*q_nocor;
        J[212] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[213] -= dqdci;              /* dwdot[OH]/d[H2O] */
        J[216] += dqdci;              /* dwdot[H2O]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[11][3] - 1)*q_nocor;
        J[422] -= dqdci;              /* dwdot[H]/d[CO] */
        J[423] -= dqdci;              /* dwdot[OH]/d[CO] */
        J[426] += dqdci;              /* dwdot[H2O]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[11][2] - 1)*q_nocor;
        J[632] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[633] -= dqdci;              /* dwdot[OH]/d[CO2] */
        J[636] += dqdci;              /* dwdot[H2O]/d[CO2] */
    }
    else {
        dqdc[0] = q_nocor;
        dqdc[1] = q_nocor;
        dqdc[2] = q_nocor + k_f*sc[3];
        dqdc[3] = q_nocor + k_f*sc[2];
        dqdc[4] = TB[11][0]*q_nocor;
        dqdc[5] = q_nocor;
        dqdc[6] = TB[11][1]*q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = q_nocor;
        dqdc[12] = TB[11][3]*q_nocor;
        dqdc[13] = q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = TB[11][2]*q_nocor;
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
        for (int k=0; k<34; k++) {
            J[35*k+2] -= dqdc[k];
            J[35*k+3] -= dqdc[k];
            J[35*k+6] += dqdc[k];
        }
    }
    J[1192] -= dqdT; /* dwdot[H]/dT */
    J[1193] -= dqdT; /* dwdot[OH]/dT */
    J[1196] += dqdT; /* dwdot[H2O]/dT */

    /*reaction 13: HCO + M => H + CO + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* FIXME: inreversible reaction in _ajac_reaction may not work*/
    /* 3-body correction factor */
    alpha = mixture + (TB[12][0] - 1)*sc[4] + (TB[12][1] - 1)*sc[6] + (TB[12][2] - 1)*sc[18] + (TB[12][3] - 1)*sc[12];
    /* forward */
    phi_f = sc[19];
    k_f = prefactor_units[12] * fwd_A[12]
                * exp(fwd_beta[12] * tc[0] - activation_units[12] * fwd_Ea[12] * invT);
    dlnkfdT = fwd_beta[12] * invT + activation_units[12] * fwd_Ea[12] * invT2;
    /* rate of progress */
    q_nocor = k_f*phi_f;
    q = alpha * q_nocor;
    dqdT = alpha * dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[12] += q; /* CO */
    wdot[19] -= q; /* HCO */
    /* for convenience */
    k_f *= alpha;
    k_r = 0.0;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[12][0] - 1)*q_nocor;
        J[142] += dqdci;              /* dwdot[H]/d[H2] */
        J[152] += dqdci;              /* dwdot[CO]/d[H2] */
        J[159] -= dqdci;              /* dwdot[HCO]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[12][1] - 1)*q_nocor;
        J[212] += dqdci;              /* dwdot[H]/d[H2O] */
        J[222] += dqdci;              /* dwdot[CO]/d[H2O] */
        J[229] -= dqdci;              /* dwdot[HCO]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[12][3] - 1)*q_nocor;
        J[422] += dqdci;              /* dwdot[H]/d[CO] */
        J[432] += dqdci;              /* dwdot[CO]/d[CO] */
        J[439] -= dqdci;              /* dwdot[HCO]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[12][2] - 1)*q_nocor;
        J[632] += dqdci;              /* dwdot[H]/d[CO2] */
        J[642] += dqdci;              /* dwdot[CO]/d[CO2] */
        J[649] -= dqdci;              /* dwdot[HCO]/d[CO2] */
        /* d()/d[HCO] */
        dqdci =  + k_f;
        J[667] += dqdci;              /* dwdot[H]/d[HCO] */
        J[677] += dqdci;              /* dwdot[CO]/d[HCO] */
        J[684] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    }
    else {
        dqdc[0] = q_nocor;
        dqdc[1] = q_nocor;
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = TB[12][0]*q_nocor;
        dqdc[5] = q_nocor;
        dqdc[6] = TB[12][1]*q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = q_nocor;
        dqdc[12] = TB[12][3]*q_nocor;
        dqdc[13] = q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = TB[12][2]*q_nocor;
        dqdc[19] = q_nocor + k_f;
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
        for (int k=0; k<34; k++) {
            J[35*k+2] += dqdc[k];
            J[35*k+12] += dqdc[k];
            J[35*k+19] -= dqdc[k];
        }
    }
    J[1192] += dqdT; /* dwdot[H]/dT */
    J[1202] += dqdT; /* dwdot[CO]/dT */
    J[1209] -= dqdT; /* dwdot[HCO]/dT */

    /*reaction 14: H + CO + M => HCO + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* FIXME: inreversible reaction in _ajac_reaction may not work*/
    /* 3-body correction factor */
    alpha = mixture + (TB[13][0] - 1)*sc[4] + (TB[13][1] - 1)*sc[6] + (TB[13][2] - 1)*sc[18] + (TB[13][3] - 1)*sc[12];
    /* forward */
    phi_f = sc[2]*sc[12];
    k_f = prefactor_units[13] * fwd_A[13]
                * exp(fwd_beta[13] * tc[0] - activation_units[13] * fwd_Ea[13] * invT);
    dlnkfdT = fwd_beta[13] * invT + activation_units[13] * fwd_Ea[13] * invT2;
    /* rate of progress */
    q_nocor = k_f*phi_f;
    q = alpha * q_nocor;
    dqdT = alpha * dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[12] -= q; /* CO */
    wdot[19] += q; /* HCO */
    /* for convenience */
    k_f *= alpha;
    k_r = 0.0;
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[12];
        J[72] -= dqdci;               /* dwdot[H]/d[H] */
        J[82] -= dqdci;               /* dwdot[CO]/d[H] */
        J[89] += dqdci;               /* dwdot[HCO]/d[H] */
        /* d()/d[H2] */
        dqdci = (TB[13][0] - 1)*q_nocor;
        J[142] -= dqdci;              /* dwdot[H]/d[H2] */
        J[152] -= dqdci;              /* dwdot[CO]/d[H2] */
        J[159] += dqdci;              /* dwdot[HCO]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[13][1] - 1)*q_nocor;
        J[212] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[222] -= dqdci;              /* dwdot[CO]/d[H2O] */
        J[229] += dqdci;              /* dwdot[HCO]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (TB[13][3] - 1)*q_nocor + k_f*sc[2];
        J[422] -= dqdci;              /* dwdot[H]/d[CO] */
        J[432] -= dqdci;              /* dwdot[CO]/d[CO] */
        J[439] += dqdci;              /* dwdot[HCO]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[13][2] - 1)*q_nocor;
        J[632] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[642] -= dqdci;              /* dwdot[CO]/d[CO2] */
        J[649] += dqdci;              /* dwdot[HCO]/d[CO2] */
    }
    else {
        dqdc[0] = q_nocor;
        dqdc[1] = q_nocor;
        dqdc[2] = q_nocor + k_f*sc[12];
        dqdc[3] = q_nocor;
        dqdc[4] = TB[13][0]*q_nocor;
        dqdc[5] = q_nocor;
        dqdc[6] = TB[13][1]*q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = q_nocor;
        dqdc[12] = TB[13][3]*q_nocor + k_f*sc[2];
        dqdc[13] = q_nocor;
        dqdc[14] = q_nocor;
        dqdc[15] = q_nocor;
        dqdc[16] = q_nocor;
        dqdc[17] = q_nocor;
        dqdc[18] = TB[13][2]*q_nocor;
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
        for (int k=0; k<34; k++) {
            J[35*k+2] -= dqdc[k];
            J[35*k+12] -= dqdc[k];
            J[35*k+19] += dqdc[k];
        }
    }
    J[1192] -= dqdT; /* dwdot[H]/dT */
    J[1202] -= dqdT; /* dwdot[CO]/dT */
    J[1209] += dqdT; /* dwdot[HCO]/dT */

    /*reaction 15: CH3O2 + M => CH3 + O2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* FIXME: inreversible reaction in _ajac_reaction may not work*/
    /* 3-body correction factor */
    alpha = mixture;
    /* forward */
    phi_f = sc[20];
    k_f = prefactor_units[14] * fwd_A[14]
                * exp(fwd_beta[14] * tc[0] - activation_units[14] * fwd_Ea[14] * invT);
    dlnkfdT = fwd_beta[14] * invT + activation_units[14] * fwd_Ea[14] * invT2;
    /* rate of progress */
    q_nocor = k_f*phi_f;
    q = alpha * q_nocor;
    dqdT = alpha * dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[11] += q; /* CH3 */
    wdot[20] -= q; /* CH3O2 */
    /* for convenience */
    k_f *= alpha;
    k_r = 0.0;
    if (consP) {
        /* d()/d[CH3O2] */
        dqdci =  + k_f;
        J[705] += dqdci;              /* dwdot[O2]/d[CH3O2] */
        J[711] += dqdci;              /* dwdot[CH3]/d[CH3O2] */
        J[720] -= dqdci;              /* dwdot[CH3O2]/d[CH3O2] */
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
        dqdc[20] = q_nocor + k_f;
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
        for (int k=0; k<34; k++) {
            J[35*k+5] += dqdc[k];
            J[35*k+11] += dqdc[k];
            J[35*k+20] -= dqdc[k];
        }
    }
    J[1195] += dqdT; /* dwdot[O2]/dT */
    J[1201] += dqdT; /* dwdot[CH3]/dT */
    J[1210] -= dqdT; /* dwdot[CH3O2]/dT */

    /*reaction 16: CH3 + O2 + M => CH3O2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* FIXME: inreversible reaction in _ajac_reaction may not work*/
    /* 3-body correction factor */
    alpha = mixture;
    /* forward */
    phi_f = sc[5]*sc[11];
    k_f = prefactor_units[15] * fwd_A[15]
                * exp(fwd_beta[15] * tc[0] - activation_units[15] * fwd_Ea[15] * invT);
    dlnkfdT = fwd_beta[15] * invT + activation_units[15] * fwd_Ea[15] * invT2;
    /* rate of progress */
    q_nocor = k_f*phi_f;
    q = alpha * q_nocor;
    dqdT = alpha * dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[11] -= q; /* CH3 */
    wdot[20] += q; /* CH3O2 */
    /* for convenience */
    k_f *= alpha;
    k_r = 0.0;
    if (consP) {
        /* d()/d[O2] */
        dqdci =  + k_f*sc[11];
        J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
        J[186] -= dqdci;              /* dwdot[CH3]/d[O2] */
        J[195] += dqdci;              /* dwdot[CH3O2]/d[O2] */
        /* d()/d[CH3] */
        dqdci =  + k_f*sc[5];
        J[390] -= dqdci;              /* dwdot[O2]/d[CH3] */
        J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
        J[405] += dqdci;              /* dwdot[CH3O2]/d[CH3] */
    }
    else {
        dqdc[0] = q_nocor;
        dqdc[1] = q_nocor;
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor;
        dqdc[5] = q_nocor + k_f*sc[11];
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = q_nocor + k_f*sc[5];
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
        dqdc[25] = q_nocor;
        dqdc[26] = q_nocor;
        dqdc[27] = q_nocor;
        dqdc[28] = q_nocor;
        dqdc[29] = q_nocor;
        dqdc[30] = q_nocor;
        dqdc[31] = q_nocor;
        dqdc[32] = q_nocor;
        dqdc[33] = q_nocor;
        for (int k=0; k<34; k++) {
            J[35*k+5] -= dqdc[k];
            J[35*k+11] -= dqdc[k];
            J[35*k+20] += dqdc[k];
        }
    }
    J[1195] -= dqdT; /* dwdot[O2]/dT */
    J[1201] -= dqdT; /* dwdot[CH3]/dT */
    J[1210] += dqdT; /* dwdot[CH3O2]/dT */

    /*reaction 17: O + H2 => H + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[4];
    k_f = prefactor_units[16] * fwd_A[16]
                * exp(fwd_beta[16] * tc[0] - activation_units[16] * fwd_Ea[16] * invT);
    dlnkfdT = fwd_beta[16] * invT + activation_units[16] * fwd_Ea[16] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[2] += q; /* H */
    wdot[3] += q; /* OH */
    wdot[4] -= q; /* H2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[4];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[37] += dqdci;               /* dwdot[H]/d[O] */
    J[38] += dqdci;               /* dwdot[OH]/d[O] */
    J[39] -= dqdci;               /* dwdot[H2]/d[O] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[1];
    J[141] -= dqdci;              /* dwdot[O]/d[H2] */
    J[142] += dqdci;              /* dwdot[H]/d[H2] */
    J[143] += dqdci;              /* dwdot[OH]/d[H2] */
    J[144] -= dqdci;              /* dwdot[H2]/d[H2] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1194] -= dqdT;              /* dwdot[H2]/dT */

    /*reaction 18: H + OH => O + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[3];
    k_f = prefactor_units[17] * fwd_A[17]
                * exp(fwd_beta[17] * tc[0] - activation_units[17] * fwd_Ea[17] * invT);
    dlnkfdT = fwd_beta[17] * invT + activation_units[17] * fwd_Ea[17] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] += q; /* O */
    wdot[2] -= q; /* H */
    wdot[3] -= q; /* OH */
    wdot[4] += q; /* H2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[3];
    J[71] += dqdci;               /* dwdot[O]/d[H] */
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[73] -= dqdci;               /* dwdot[OH]/d[H] */
    J[74] += dqdci;               /* dwdot[H2]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[2];
    J[106] += dqdci;              /* dwdot[O]/d[OH] */
    J[107] -= dqdci;              /* dwdot[H]/d[OH] */
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[109] += dqdci;              /* dwdot[H2]/d[OH] */
    /* d()/dT */
    J[1191] += dqdT;              /* dwdot[O]/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */

    /*reaction 19: OH + H2 => H + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[4];
    k_f = prefactor_units[18] * fwd_A[18]
                * exp(fwd_beta[18] * tc[0] - activation_units[18] * fwd_Ea[18] * invT);
    dlnkfdT = fwd_beta[18] * invT + activation_units[18] * fwd_Ea[18] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[3] -= q; /* OH */
    wdot[4] -= q; /* H2 */
    wdot[6] += q; /* H2O */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[4];
    J[107] += dqdci;              /* dwdot[H]/d[OH] */
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[109] -= dqdci;              /* dwdot[H2]/d[OH] */
    J[111] += dqdci;              /* dwdot[H2O]/d[OH] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[3];
    J[142] += dqdci;              /* dwdot[H]/d[H2] */
    J[143] -= dqdci;              /* dwdot[OH]/d[H2] */
    J[144] -= dqdci;              /* dwdot[H2]/d[H2] */
    J[146] += dqdci;              /* dwdot[H2O]/d[H2] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1194] -= dqdT;              /* dwdot[H2]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */

    /*reaction 20: H + H2O => OH + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[6];
    k_f = prefactor_units[19] * fwd_A[19]
                * exp(fwd_beta[19] * tc[0] - activation_units[19] * fwd_Ea[19] * invT);
    dlnkfdT = fwd_beta[19] * invT + activation_units[19] * fwd_Ea[19] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[3] += q; /* OH */
    wdot[4] += q; /* H2 */
    wdot[6] -= q; /* H2O */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[73] += dqdci;               /* dwdot[OH]/d[H] */
    J[74] += dqdci;               /* dwdot[H2]/d[H] */
    J[76] -= dqdci;               /* dwdot[H2O]/d[H] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[2];
    J[212] -= dqdci;              /* dwdot[H]/d[H2O] */
    J[213] += dqdci;              /* dwdot[OH]/d[H2O] */
    J[214] += dqdci;              /* dwdot[H2]/d[H2O] */
    J[216] -= dqdci;              /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */
    J[1196] -= dqdT;              /* dwdot[H2O]/dT */

    /*reaction 21: O + H2O => 2.000000 OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[6];
    k_f = prefactor_units[20] * fwd_A[20]
                * exp(fwd_beta[20] * tc[0] - activation_units[20] * fwd_Ea[20] * invT);
    dlnkfdT = fwd_beta[20] * invT + activation_units[20] * fwd_Ea[20] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[3] += 2 * q; /* OH */
    wdot[6] -= q; /* H2O */
    /* d()/d[O] */
    dqdci =  + k_f*sc[6];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[38] += 2 * dqdci;           /* dwdot[OH]/d[O] */
    J[41] -= dqdci;               /* dwdot[H2O]/d[O] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[1];
    J[211] -= dqdci;              /* dwdot[O]/d[H2O] */
    J[213] += 2 * dqdci;          /* dwdot[OH]/d[H2O] */
    J[216] -= dqdci;              /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1193] += 2 * dqdT;          /* dwdot[OH]/dT */
    J[1196] -= dqdT;              /* dwdot[H2O]/dT */

    /*reaction 22: 2.000000 OH => O + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[3], 2.000000);
    k_f = prefactor_units[21] * fwd_A[21]
                * exp(fwd_beta[21] * tc[0] - activation_units[21] * fwd_Ea[21] * invT);
    dlnkfdT = fwd_beta[21] * invT + activation_units[21] * fwd_Ea[21] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] += q; /* O */
    wdot[3] -= 2 * q; /* OH */
    wdot[6] += q; /* H2O */
    /* d()/d[OH] */
    dqdci =  + k_f*2.000000*sc[3];
    J[106] += dqdci;              /* dwdot[O]/d[OH] */
    J[108] += -2 * dqdci;         /* dwdot[OH]/d[OH] */
    J[111] += dqdci;              /* dwdot[H2O]/d[OH] */
    /* d()/dT */
    J[1191] += dqdT;              /* dwdot[O]/dT */
    J[1193] += -2 * dqdT;         /* dwdot[OH]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */

    /*reaction 23: H + O2 => O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[5];
    k_f = prefactor_units[22] * fwd_A[22]
                * exp(fwd_beta[22] * tc[0] - activation_units[22] * fwd_Ea[22] * invT);
    dlnkfdT = fwd_beta[22] * invT + activation_units[22] * fwd_Ea[22] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] += q; /* O */
    wdot[2] -= q; /* H */
    wdot[3] += q; /* OH */
    wdot[5] -= q; /* O2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[5];
    J[71] += dqdci;               /* dwdot[O]/d[H] */
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[73] += dqdci;               /* dwdot[OH]/d[H] */
    J[75] -= dqdci;               /* dwdot[O2]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[2];
    J[176] += dqdci;              /* dwdot[O]/d[O2] */
    J[177] -= dqdci;              /* dwdot[H]/d[O2] */
    J[178] += dqdci;              /* dwdot[OH]/d[O2] */
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    /* d()/dT */
    J[1191] += dqdT;              /* dwdot[O]/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */

    /*reaction 24: O + OH => H + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = prefactor_units[23] * fwd_A[23]
                * exp(fwd_beta[23] * tc[0] - activation_units[23] * fwd_Ea[23] * invT);
    dlnkfdT = fwd_beta[23] * invT + activation_units[23] * fwd_Ea[23] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[2] += q; /* H */
    wdot[3] -= q; /* OH */
    wdot[5] += q; /* O2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[3];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[37] += dqdci;               /* dwdot[H]/d[O] */
    J[38] -= dqdci;               /* dwdot[OH]/d[O] */
    J[40] += dqdci;               /* dwdot[O2]/d[O] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[1];
    J[106] -= dqdci;              /* dwdot[O]/d[OH] */
    J[107] += dqdci;              /* dwdot[H]/d[OH] */
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[110] += dqdci;              /* dwdot[O2]/d[OH] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */

    /*reaction 25: HO2 + H => 2.000000 OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[8];
    k_f = prefactor_units[24] * fwd_A[24]
                * exp(fwd_beta[24] * tc[0] - activation_units[24] * fwd_Ea[24] * invT);
    dlnkfdT = fwd_beta[24] * invT + activation_units[24] * fwd_Ea[24] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[3] += 2 * q; /* OH */
    wdot[8] -= q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[8];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[73] += 2 * dqdci;           /* dwdot[OH]/d[H] */
    J[78] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[2];
    J[282] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[283] += 2 * dqdci;          /* dwdot[OH]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1193] += 2 * dqdT;          /* dwdot[OH]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */

    /*reaction 26: 2.000000 HO2 => H2O2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[8], 2.000000);
    k_f = prefactor_units[25] * fwd_A[25]
                * exp(fwd_beta[25] * tc[0] - activation_units[25] * fwd_Ea[25] * invT);
    dlnkfdT = fwd_beta[25] * invT + activation_units[25] * fwd_Ea[25] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[7] += q; /* H2O2 */
    wdot[8] -= 2 * q; /* HO2 */
    /* d()/d[HO2] */
    dqdci =  + k_f*2.000000*sc[8];
    J[285] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[287] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    J[288] += -2 * dqdci;         /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1197] += dqdT;              /* dwdot[H2O2]/dT */
    J[1198] += -2 * dqdT;         /* dwdot[HO2]/dT */

    /*reaction 27: HO2 + H => H2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[8];
    k_f = prefactor_units[26] * fwd_A[26]
                * exp(fwd_beta[26] * tc[0] - activation_units[26] * fwd_Ea[26] * invT);
    dlnkfdT = fwd_beta[26] * invT + activation_units[26] * fwd_Ea[26] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[4] += q; /* H2 */
    wdot[5] += q; /* O2 */
    wdot[8] -= q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[8];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[74] += dqdci;               /* dwdot[H2]/d[H] */
    J[75] += dqdci;               /* dwdot[O2]/d[H] */
    J[78] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[2];
    J[282] -= dqdci;              /* dwdot[H]/d[HO2] */
    J[284] += dqdci;              /* dwdot[H2]/d[HO2] */
    J[285] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */

    /*reaction 28: HO2 + OH => H2O + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[8];
    k_f = prefactor_units[27] * fwd_A[27]
                * exp(fwd_beta[27] * tc[0] - activation_units[27] * fwd_Ea[27] * invT);
    dlnkfdT = fwd_beta[27] * invT + activation_units[27] * fwd_Ea[27] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[5] += q; /* O2 */
    wdot[6] += q; /* H2O */
    wdot[8] -= q; /* HO2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[8];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[110] += dqdci;              /* dwdot[O2]/d[OH] */
    J[111] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[113] -= dqdci;              /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[3];
    J[283] -= dqdci;              /* dwdot[OH]/d[HO2] */
    J[285] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[286] += dqdci;              /* dwdot[H2O]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */

    /*reaction 29: H2O + O2 => HO2 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[6];
    k_f = prefactor_units[28] * fwd_A[28]
                * exp(fwd_beta[28] * tc[0] - activation_units[28] * fwd_Ea[28] * invT);
    dlnkfdT = fwd_beta[28] * invT + activation_units[28] * fwd_Ea[28] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] += q; /* OH */
    wdot[5] -= q; /* O2 */
    wdot[6] -= q; /* H2O */
    wdot[8] += q; /* HO2 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[6];
    J[178] += dqdci;              /* dwdot[OH]/d[O2] */
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[181] -= dqdci;              /* dwdot[H2O]/d[O2] */
    J[183] += dqdci;              /* dwdot[HO2]/d[O2] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[5];
    J[213] += dqdci;              /* dwdot[OH]/d[H2O] */
    J[215] -= dqdci;              /* dwdot[O2]/d[H2O] */
    J[216] -= dqdci;              /* dwdot[H2O]/d[H2O] */
    J[218] += dqdci;              /* dwdot[HO2]/d[H2O] */
    /* d()/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1196] -= dqdT;              /* dwdot[H2O]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */

    /*reaction 30: 2.000000 HO2 => H2O2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[8], 2.000000);
    k_f = prefactor_units[29] * fwd_A[29]
                * exp(fwd_beta[29] * tc[0] - activation_units[29] * fwd_Ea[29] * invT);
    dlnkfdT = fwd_beta[29] * invT + activation_units[29] * fwd_Ea[29] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[7] += q; /* H2O2 */
    wdot[8] -= 2 * q; /* HO2 */
    /* d()/d[HO2] */
    dqdci =  + k_f*2.000000*sc[8];
    J[285] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[287] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    J[288] += -2 * dqdci;         /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1197] += dqdT;              /* dwdot[H2O2]/dT */
    J[1198] += -2 * dqdT;         /* dwdot[HO2]/dT */

    /*reaction 31: HO2 + O => OH + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[8];
    k_f = prefactor_units[30] * fwd_A[30]
                * exp(fwd_beta[30] * tc[0] - activation_units[30] * fwd_Ea[30] * invT);
    dlnkfdT = fwd_beta[30] * invT + activation_units[30] * fwd_Ea[30] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[3] += q; /* OH */
    wdot[5] += q; /* O2 */
    wdot[8] -= q; /* HO2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[8];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[38] += dqdci;               /* dwdot[OH]/d[O] */
    J[40] += dqdci;               /* dwdot[O2]/d[O] */
    J[43] -= dqdci;               /* dwdot[HO2]/d[O] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[1];
    J[281] -= dqdci;              /* dwdot[O]/d[HO2] */
    J[283] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[285] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */

    /*reaction 32: OH + O2 => HO2 + O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[5];
    k_f = prefactor_units[31] * fwd_A[31]
                * exp(fwd_beta[31] * tc[0] - activation_units[31] * fwd_Ea[31] * invT);
    dlnkfdT = fwd_beta[31] * invT + activation_units[31] * fwd_Ea[31] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] += q; /* O */
    wdot[3] -= q; /* OH */
    wdot[5] -= q; /* O2 */
    wdot[8] += q; /* HO2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[5];
    J[106] += dqdci;              /* dwdot[O]/d[OH] */
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[110] -= dqdci;              /* dwdot[O2]/d[OH] */
    J[113] += dqdci;              /* dwdot[HO2]/d[OH] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[3];
    J[176] += dqdci;              /* dwdot[O]/d[O2] */
    J[178] -= dqdci;              /* dwdot[OH]/d[O2] */
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[183] += dqdci;              /* dwdot[HO2]/d[O2] */
    /* d()/dT */
    J[1191] += dqdT;              /* dwdot[O]/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */

    /*reaction 33: H2O2 + OH => H2O + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[7];
    k_f = prefactor_units[32] * fwd_A[32]
                * exp(fwd_beta[32] * tc[0] - activation_units[32] * fwd_Ea[32] * invT);
    dlnkfdT = fwd_beta[32] * invT + activation_units[32] * fwd_Ea[32] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[7] -= q; /* H2O2 */
    wdot[8] += q; /* HO2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[7];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[111] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[112] -= dqdci;              /* dwdot[H2O2]/d[OH] */
    J[113] += dqdci;              /* dwdot[HO2]/d[OH] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[3];
    J[248] -= dqdci;              /* dwdot[OH]/d[H2O2] */
    J[251] += dqdci;              /* dwdot[H2O]/d[H2O2] */
    J[252] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    J[253] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */
    J[1197] -= dqdT;              /* dwdot[H2O2]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */

    /*reaction 34: H2O2 + H => H2O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[7];
    k_f = prefactor_units[33] * fwd_A[33]
                * exp(fwd_beta[33] * tc[0] - activation_units[33] * fwd_Ea[33] * invT);
    dlnkfdT = fwd_beta[33] * invT + activation_units[33] * fwd_Ea[33] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[3] += q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[73] += dqdci;               /* dwdot[OH]/d[H] */
    J[76] += dqdci;               /* dwdot[H2O]/d[H] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[2];
    J[247] -= dqdci;              /* dwdot[H]/d[H2O2] */
    J[248] += dqdci;              /* dwdot[OH]/d[H2O2] */
    J[251] += dqdci;              /* dwdot[H2O]/d[H2O2] */
    J[252] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */
    J[1197] -= dqdT;              /* dwdot[H2O2]/dT */

    /*reaction 35: H2O2 + OH => H2O + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[7];
    k_f = prefactor_units[34] * fwd_A[34]
                * exp(fwd_beta[34] * tc[0] - activation_units[34] * fwd_Ea[34] * invT);
    dlnkfdT = fwd_beta[34] * invT + activation_units[34] * fwd_Ea[34] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[7] -= q; /* H2O2 */
    wdot[8] += q; /* HO2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[7];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[111] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[112] -= dqdci;              /* dwdot[H2O2]/d[OH] */
    J[113] += dqdci;              /* dwdot[HO2]/d[OH] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[3];
    J[248] -= dqdci;              /* dwdot[OH]/d[H2O2] */
    J[251] += dqdci;              /* dwdot[H2O]/d[H2O2] */
    J[252] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    J[253] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */
    J[1197] -= dqdT;              /* dwdot[H2O2]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */

    /*reaction 36: H2O + HO2 => H2O2 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[8];
    k_f = prefactor_units[35] * fwd_A[35]
                * exp(fwd_beta[35] * tc[0] - activation_units[35] * fwd_Ea[35] * invT);
    dlnkfdT = fwd_beta[35] * invT + activation_units[35] * fwd_Ea[35] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] += q; /* OH */
    wdot[6] -= q; /* H2O */
    wdot[7] += q; /* H2O2 */
    wdot[8] -= q; /* HO2 */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[8];
    J[213] += dqdci;              /* dwdot[OH]/d[H2O] */
    J[216] -= dqdci;              /* dwdot[H2O]/d[H2O] */
    J[217] += dqdci;              /* dwdot[H2O2]/d[H2O] */
    J[218] -= dqdci;              /* dwdot[HO2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[6];
    J[283] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[286] -= dqdci;              /* dwdot[H2O]/d[HO2] */
    J[287] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1196] -= dqdT;              /* dwdot[H2O]/dT */
    J[1197] += dqdT;              /* dwdot[H2O2]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */

    /*reaction 37: H2O2 + O => OH + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[7];
    k_f = prefactor_units[36] * fwd_A[36]
                * exp(fwd_beta[36] * tc[0] - activation_units[36] * fwd_Ea[36] * invT);
    dlnkfdT = fwd_beta[36] * invT + activation_units[36] * fwd_Ea[36] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[3] += q; /* OH */
    wdot[7] -= q; /* H2O2 */
    wdot[8] += q; /* HO2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[7];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[38] += dqdci;               /* dwdot[OH]/d[O] */
    J[42] -= dqdci;               /* dwdot[H2O2]/d[O] */
    J[43] += dqdci;               /* dwdot[HO2]/d[O] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[1];
    J[246] -= dqdci;              /* dwdot[O]/d[H2O2] */
    J[248] += dqdci;              /* dwdot[OH]/d[H2O2] */
    J[252] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    J[253] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1197] -= dqdT;              /* dwdot[H2O2]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */

    /*reaction 38: H2O2 + H => H2 + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[7];
    k_f = prefactor_units[37] * fwd_A[37]
                * exp(fwd_beta[37] * tc[0] - activation_units[37] * fwd_Ea[37] * invT);
    dlnkfdT = fwd_beta[37] * invT + activation_units[37] * fwd_Ea[37] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[4] += q; /* H2 */
    wdot[7] -= q; /* H2O2 */
    wdot[8] += q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[74] += dqdci;               /* dwdot[H2]/d[H] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H] */
    J[78] += dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[2];
    J[247] -= dqdci;              /* dwdot[H]/d[H2O2] */
    J[249] += dqdci;              /* dwdot[H2]/d[H2O2] */
    J[252] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    J[253] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */
    J[1197] -= dqdT;              /* dwdot[H2O2]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */

    /*reaction 39: H2 + HO2 => H2O2 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[8];
    k_f = prefactor_units[38] * fwd_A[38]
                * exp(fwd_beta[38] * tc[0] - activation_units[38] * fwd_Ea[38] * invT);
    dlnkfdT = fwd_beta[38] * invT + activation_units[38] * fwd_Ea[38] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[4] -= q; /* H2 */
    wdot[7] += q; /* H2O2 */
    wdot[8] -= q; /* HO2 */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[8];
    J[142] += dqdci;              /* dwdot[H]/d[H2] */
    J[144] -= dqdci;              /* dwdot[H2]/d[H2] */
    J[147] += dqdci;              /* dwdot[H2O2]/d[H2] */
    J[148] -= dqdci;              /* dwdot[HO2]/d[H2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[4];
    J[282] += dqdci;              /* dwdot[H]/d[HO2] */
    J[284] -= dqdci;              /* dwdot[H2]/d[HO2] */
    J[287] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1194] -= dqdT;              /* dwdot[H2]/dT */
    J[1197] += dqdT;              /* dwdot[H2O2]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */

    /*reaction 40: CH2GSG + OH => CH2O + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[9];
    k_f = prefactor_units[39] * fwd_A[39]
                * exp(fwd_beta[39] * tc[0] - activation_units[39] * fwd_Ea[39] * invT);
    dlnkfdT = fwd_beta[39] * invT + activation_units[39] * fwd_Ea[39] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[3] -= q; /* OH */
    wdot[9] -= q; /* CH2GSG */
    wdot[10] += q; /* CH2O */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[9];
    J[107] += dqdci;              /* dwdot[H]/d[OH] */
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[114] -= dqdci;              /* dwdot[CH2GSG]/d[OH] */
    J[115] += dqdci;              /* dwdot[CH2O]/d[OH] */
    /* d()/d[CH2GSG] */
    dqdci =  + k_f*sc[3];
    J[317] += dqdci;              /* dwdot[H]/d[CH2GSG] */
    J[318] -= dqdci;              /* dwdot[OH]/d[CH2GSG] */
    J[324] -= dqdci;              /* dwdot[CH2GSG]/d[CH2GSG] */
    J[325] += dqdci;              /* dwdot[CH2O]/d[CH2GSG] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1199] -= dqdT;              /* dwdot[CH2GSG]/dT */
    J[1200] += dqdT;              /* dwdot[CH2O]/dT */

    /*reaction 41: CH2GSG + H2 => CH3 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[9];
    k_f = prefactor_units[40] * fwd_A[40]
                * exp(fwd_beta[40] * tc[0] - activation_units[40] * fwd_Ea[40] * invT);
    dlnkfdT = fwd_beta[40] * invT + activation_units[40] * fwd_Ea[40] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[4] -= q; /* H2 */
    wdot[9] -= q; /* CH2GSG */
    wdot[11] += q; /* CH3 */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[9];
    J[142] += dqdci;              /* dwdot[H]/d[H2] */
    J[144] -= dqdci;              /* dwdot[H2]/d[H2] */
    J[149] -= dqdci;              /* dwdot[CH2GSG]/d[H2] */
    J[151] += dqdci;              /* dwdot[CH3]/d[H2] */
    /* d()/d[CH2GSG] */
    dqdci =  + k_f*sc[4];
    J[317] += dqdci;              /* dwdot[H]/d[CH2GSG] */
    J[319] -= dqdci;              /* dwdot[H2]/d[CH2GSG] */
    J[324] -= dqdci;              /* dwdot[CH2GSG]/d[CH2GSG] */
    J[326] += dqdci;              /* dwdot[CH3]/d[CH2GSG] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1194] -= dqdT;              /* dwdot[H2]/dT */
    J[1199] -= dqdT;              /* dwdot[CH2GSG]/dT */
    J[1201] += dqdT;              /* dwdot[CH3]/dT */

    /*reaction 42: CH3 + H => CH2GSG + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[11];
    k_f = prefactor_units[41] * fwd_A[41]
                * exp(fwd_beta[41] * tc[0] - activation_units[41] * fwd_Ea[41] * invT);
    dlnkfdT = fwd_beta[41] * invT + activation_units[41] * fwd_Ea[41] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[4] += q; /* H2 */
    wdot[9] += q; /* CH2GSG */
    wdot[11] -= q; /* CH3 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[11];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[74] += dqdci;               /* dwdot[H2]/d[H] */
    J[79] += dqdci;               /* dwdot[CH2GSG]/d[H] */
    J[81] -= dqdci;               /* dwdot[CH3]/d[H] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[2];
    J[387] -= dqdci;              /* dwdot[H]/d[CH3] */
    J[389] += dqdci;              /* dwdot[H2]/d[CH3] */
    J[394] += dqdci;              /* dwdot[CH2GSG]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */
    J[1199] += dqdT;              /* dwdot[CH2GSG]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */

    /*reaction 43: CH2GSG + O2 => CO + OH + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[9];
    k_f = prefactor_units[42] * fwd_A[42]
                * exp(fwd_beta[42] * tc[0] - activation_units[42] * fwd_Ea[42] * invT);
    dlnkfdT = fwd_beta[42] * invT + activation_units[42] * fwd_Ea[42] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[3] += q; /* OH */
    wdot[5] -= q; /* O2 */
    wdot[9] -= q; /* CH2GSG */
    wdot[12] += q; /* CO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[9];
    J[177] += dqdci;              /* dwdot[H]/d[O2] */
    J[178] += dqdci;              /* dwdot[OH]/d[O2] */
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[184] -= dqdci;              /* dwdot[CH2GSG]/d[O2] */
    J[187] += dqdci;              /* dwdot[CO]/d[O2] */
    /* d()/d[CH2GSG] */
    dqdci =  + k_f*sc[5];
    J[317] += dqdci;              /* dwdot[H]/d[CH2GSG] */
    J[318] += dqdci;              /* dwdot[OH]/d[CH2GSG] */
    J[320] -= dqdci;              /* dwdot[O2]/d[CH2GSG] */
    J[324] -= dqdci;              /* dwdot[CH2GSG]/d[CH2GSG] */
    J[327] += dqdci;              /* dwdot[CO]/d[CH2GSG] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1199] -= dqdT;              /* dwdot[CH2GSG]/dT */
    J[1202] += dqdT;              /* dwdot[CO]/dT */

    /*reaction 44: CH2GSG + O => CO + 2.000000 H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[9];
    k_f = prefactor_units[43] * fwd_A[43]
                * exp(fwd_beta[43] * tc[0] - activation_units[43] * fwd_Ea[43] * invT);
    dlnkfdT = fwd_beta[43] * invT + activation_units[43] * fwd_Ea[43] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[2] += 2 * q; /* H */
    wdot[9] -= q; /* CH2GSG */
    wdot[12] += q; /* CO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[9];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[37] += 2 * dqdci;           /* dwdot[H]/d[O] */
    J[44] -= dqdci;               /* dwdot[CH2GSG]/d[O] */
    J[47] += dqdci;               /* dwdot[CO]/d[O] */
    /* d()/d[CH2GSG] */
    dqdci =  + k_f*sc[1];
    J[316] -= dqdci;              /* dwdot[O]/d[CH2GSG] */
    J[317] += 2 * dqdci;          /* dwdot[H]/d[CH2GSG] */
    J[324] -= dqdci;              /* dwdot[CH2GSG]/d[CH2GSG] */
    J[327] += dqdci;              /* dwdot[CO]/d[CH2GSG] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1192] += 2 * dqdT;          /* dwdot[H]/dT */
    J[1199] -= dqdT;              /* dwdot[CH2GSG]/dT */
    J[1202] += dqdT;              /* dwdot[CO]/dT */

    /*reaction 45: CH3 + HO2 => CH3O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[11];
    k_f = prefactor_units[44] * fwd_A[44]
                * exp(fwd_beta[44] * tc[0] - activation_units[44] * fwd_Ea[44] * invT);
    dlnkfdT = fwd_beta[44] * invT + activation_units[44] * fwd_Ea[44] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] += q; /* OH */
    wdot[8] -= q; /* HO2 */
    wdot[11] -= q; /* CH3 */
    wdot[13] += q; /* CH3O */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[11];
    J[283] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[291] -= dqdci;              /* dwdot[CH3]/d[HO2] */
    J[293] += dqdci;              /* dwdot[CH3O]/d[HO2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[8];
    J[388] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[393] -= dqdci;              /* dwdot[HO2]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[398] += dqdci;              /* dwdot[CH3O]/d[CH3] */
    /* d()/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1203] += dqdT;              /* dwdot[CH3O]/dT */

    /*reaction 46: CH3 + O2 => CH2O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[11];
    k_f = prefactor_units[45] * fwd_A[45]
                * exp(fwd_beta[45] * tc[0] - activation_units[45] * fwd_Ea[45] * invT);
    dlnkfdT = fwd_beta[45] * invT + activation_units[45] * fwd_Ea[45] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] += q; /* OH */
    wdot[5] -= q; /* O2 */
    wdot[10] += q; /* CH2O */
    wdot[11] -= q; /* CH3 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[11];
    J[178] += dqdci;              /* dwdot[OH]/d[O2] */
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[185] += dqdci;              /* dwdot[CH2O]/d[O2] */
    J[186] -= dqdci;              /* dwdot[CH3]/d[O2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[5];
    J[388] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[390] -= dqdci;              /* dwdot[O2]/d[CH3] */
    J[395] += dqdci;              /* dwdot[CH2O]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1200] += dqdT;              /* dwdot[CH2O]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */

    /*reaction 47: CH3 + O2 => CH3O + O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[11];
    k_f = prefactor_units[46] * fwd_A[46]
                * exp(fwd_beta[46] * tc[0] - activation_units[46] * fwd_Ea[46] * invT);
    dlnkfdT = fwd_beta[46] * invT + activation_units[46] * fwd_Ea[46] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] += q; /* O */
    wdot[5] -= q; /* O2 */
    wdot[11] -= q; /* CH3 */
    wdot[13] += q; /* CH3O */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[11];
    J[176] += dqdci;              /* dwdot[O]/d[O2] */
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[186] -= dqdci;              /* dwdot[CH3]/d[O2] */
    J[188] += dqdci;              /* dwdot[CH3O]/d[O2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[5];
    J[386] += dqdci;              /* dwdot[O]/d[CH3] */
    J[390] -= dqdci;              /* dwdot[O2]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[398] += dqdci;              /* dwdot[CH3O]/d[CH3] */
    /* d()/dT */
    J[1191] += dqdT;              /* dwdot[O]/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1203] += dqdT;              /* dwdot[CH3O]/dT */

    /*reaction 48: CH3O + O => CH3 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[13];
    k_f = prefactor_units[47] * fwd_A[47]
                * exp(fwd_beta[47] * tc[0] - activation_units[47] * fwd_Ea[47] * invT);
    dlnkfdT = fwd_beta[47] * invT + activation_units[47] * fwd_Ea[47] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[5] += q; /* O2 */
    wdot[11] += q; /* CH3 */
    wdot[13] -= q; /* CH3O */
    /* d()/d[O] */
    dqdci =  + k_f*sc[13];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[40] += dqdci;               /* dwdot[O2]/d[O] */
    J[46] += dqdci;               /* dwdot[CH3]/d[O] */
    J[48] -= dqdci;               /* dwdot[CH3O]/d[O] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[1];
    J[456] -= dqdci;              /* dwdot[O]/d[CH3O] */
    J[460] += dqdci;              /* dwdot[O2]/d[CH3O] */
    J[466] += dqdci;              /* dwdot[CH3]/d[CH3O] */
    J[468] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1201] += dqdT;              /* dwdot[CH3]/dT */
    J[1203] -= dqdT;              /* dwdot[CH3O]/dT */

    /*reaction 49: 2.000000 CH3 <=> H + C2H5 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[11], 2.000000);
    k_f = prefactor_units[48] * fwd_A[48]
                * exp(fwd_beta[48] * tc[0] - activation_units[48] * fwd_Ea[48] * invT);
    dlnkfdT = fwd_beta[48] * invT + activation_units[48] * fwd_Ea[48] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[14];
    Kc = exp(-g_RT[2] + 2.000000*g_RT[11] - g_RT[14]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[11]) + (h_RT[2] + h_RT[14]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[11] -= 2 * q; /* CH3 */
    wdot[14] += q; /* C2H5 */
    /* d()/d[H] */
    dqdci =  - k_r*sc[14];
    J[72] += dqdci;               /* dwdot[H]/d[H] */
    J[81] += -2 * dqdci;          /* dwdot[CH3]/d[H] */
    J[84] += dqdci;               /* dwdot[C2H5]/d[H] */
    /* d()/d[CH3] */
    dqdci =  + k_f*2.000000*sc[11];
    J[387] += dqdci;              /* dwdot[H]/d[CH3] */
    J[396] += -2 * dqdci;         /* dwdot[CH3]/d[CH3] */
    J[399] += dqdci;              /* dwdot[C2H5]/d[CH3] */
    /* d()/d[C2H5] */
    dqdci =  - k_r*sc[2];
    J[492] += dqdci;              /* dwdot[H]/d[C2H5] */
    J[501] += -2 * dqdci;         /* dwdot[CH3]/d[C2H5] */
    J[504] += dqdci;              /* dwdot[C2H5]/d[C2H5] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1201] += -2 * dqdT;         /* dwdot[CH3]/dT */
    J[1204] += dqdT;              /* dwdot[C2H5]/dT */

    /*reaction 50: CH3 + HO2 => CH4 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[11];
    k_f = prefactor_units[49] * fwd_A[49]
                * exp(fwd_beta[49] * tc[0] - activation_units[49] * fwd_Ea[49] * invT);
    dlnkfdT = fwd_beta[49] * invT + activation_units[49] * fwd_Ea[49] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[8] -= q; /* HO2 */
    wdot[11] -= q; /* CH3 */
    wdot[15] += q; /* CH4 */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[11];
    J[285] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[291] -= dqdci;              /* dwdot[CH3]/d[HO2] */
    J[295] += dqdci;              /* dwdot[CH4]/d[HO2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[8];
    J[390] += dqdci;              /* dwdot[O2]/d[CH3] */
    J[393] -= dqdci;              /* dwdot[HO2]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[400] += dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1205] += dqdT;              /* dwdot[CH4]/dT */

    /*reaction 51: CH3 + O => CH2O + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[11];
    k_f = prefactor_units[50] * fwd_A[50]
                * exp(fwd_beta[50] * tc[0] - activation_units[50] * fwd_Ea[50] * invT);
    dlnkfdT = fwd_beta[50] * invT + activation_units[50] * fwd_Ea[50] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[2] += q; /* H */
    wdot[10] += q; /* CH2O */
    wdot[11] -= q; /* CH3 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[11];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[37] += dqdci;               /* dwdot[H]/d[O] */
    J[45] += dqdci;               /* dwdot[CH2O]/d[O] */
    J[46] -= dqdci;               /* dwdot[CH3]/d[O] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[1];
    J[386] -= dqdci;              /* dwdot[O]/d[CH3] */
    J[387] += dqdci;              /* dwdot[H]/d[CH3] */
    J[395] += dqdci;              /* dwdot[CH2O]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1200] += dqdT;              /* dwdot[CH2O]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */

    /*reaction 52: CH3 + OH => CH2O + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[11];
    k_f = prefactor_units[51] * fwd_A[51]
                * exp(fwd_beta[51] * tc[0] - activation_units[51] * fwd_Ea[51] * invT);
    dlnkfdT = fwd_beta[51] * invT + activation_units[51] * fwd_Ea[51] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[4] += q; /* H2 */
    wdot[10] += q; /* CH2O */
    wdot[11] -= q; /* CH3 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[11];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[109] += dqdci;              /* dwdot[H2]/d[OH] */
    J[115] += dqdci;              /* dwdot[CH2O]/d[OH] */
    J[116] -= dqdci;              /* dwdot[CH3]/d[OH] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[3];
    J[388] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[389] += dqdci;              /* dwdot[H2]/d[CH3] */
    J[395] += dqdci;              /* dwdot[CH2O]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */
    J[1200] += dqdT;              /* dwdot[CH2O]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */

    /*reaction 53: CH3 + OH => CH2GSG + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[11];
    k_f = prefactor_units[52] * fwd_A[52]
                * exp(fwd_beta[52] * tc[0] - activation_units[52] * fwd_Ea[52] * invT);
    dlnkfdT = fwd_beta[52] * invT + activation_units[52] * fwd_Ea[52] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[9] += q; /* CH2GSG */
    wdot[11] -= q; /* CH3 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[11];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[111] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[114] += dqdci;              /* dwdot[CH2GSG]/d[OH] */
    J[116] -= dqdci;              /* dwdot[CH3]/d[OH] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[3];
    J[388] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[391] += dqdci;              /* dwdot[H2O]/d[CH3] */
    J[394] += dqdci;              /* dwdot[CH2GSG]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */
    J[1199] += dqdT;              /* dwdot[CH2GSG]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */

    /*reaction 54: CH2GSG + H2O => CH3 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[9];
    k_f = prefactor_units[53] * fwd_A[53]
                * exp(fwd_beta[53] * tc[0] - activation_units[53] * fwd_Ea[53] * invT);
    dlnkfdT = fwd_beta[53] * invT + activation_units[53] * fwd_Ea[53] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] += q; /* OH */
    wdot[6] -= q; /* H2O */
    wdot[9] -= q; /* CH2GSG */
    wdot[11] += q; /* CH3 */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[9];
    J[213] += dqdci;              /* dwdot[OH]/d[H2O] */
    J[216] -= dqdci;              /* dwdot[H2O]/d[H2O] */
    J[219] -= dqdci;              /* dwdot[CH2GSG]/d[H2O] */
    J[221] += dqdci;              /* dwdot[CH3]/d[H2O] */
    /* d()/d[CH2GSG] */
    dqdci =  + k_f*sc[6];
    J[318] += dqdci;              /* dwdot[OH]/d[CH2GSG] */
    J[321] -= dqdci;              /* dwdot[H2O]/d[CH2GSG] */
    J[324] -= dqdci;              /* dwdot[CH2GSG]/d[CH2GSG] */
    J[326] += dqdci;              /* dwdot[CH3]/d[CH2GSG] */
    /* d()/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1196] -= dqdT;              /* dwdot[H2O]/dT */
    J[1199] -= dqdT;              /* dwdot[CH2GSG]/dT */
    J[1201] += dqdT;              /* dwdot[CH3]/dT */

    /*reaction 55: CH3 + H2O2 => CH4 + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[7]*sc[11];
    k_f = prefactor_units[54] * fwd_A[54]
                * exp(fwd_beta[54] * tc[0] - activation_units[54] * fwd_Ea[54] * invT);
    dlnkfdT = fwd_beta[54] * invT + activation_units[54] * fwd_Ea[54] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[7] -= q; /* H2O2 */
    wdot[8] += q; /* HO2 */
    wdot[11] -= q; /* CH3 */
    wdot[15] += q; /* CH4 */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[11];
    J[252] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    J[253] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[256] -= dqdci;              /* dwdot[CH3]/d[H2O2] */
    J[260] += dqdci;              /* dwdot[CH4]/d[H2O2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[7];
    J[392] -= dqdci;              /* dwdot[H2O2]/d[CH3] */
    J[393] += dqdci;              /* dwdot[HO2]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[400] += dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/dT */
    J[1197] -= dqdT;              /* dwdot[H2O2]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1205] += dqdT;              /* dwdot[CH4]/dT */

    /*reaction 56: CH2GSG + CH3 => C2H4 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[11];
    k_f = prefactor_units[55] * fwd_A[55]
                * exp(fwd_beta[55] * tc[0] - activation_units[55] * fwd_Ea[55] * invT);
    dlnkfdT = fwd_beta[55] * invT + activation_units[55] * fwd_Ea[55] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[9] -= q; /* CH2GSG */
    wdot[11] -= q; /* CH3 */
    wdot[16] += q; /* C2H4 */
    /* d()/d[CH2GSG] */
    dqdci =  + k_f*sc[11];
    J[317] += dqdci;              /* dwdot[H]/d[CH2GSG] */
    J[324] -= dqdci;              /* dwdot[CH2GSG]/d[CH2GSG] */
    J[326] -= dqdci;              /* dwdot[CH3]/d[CH2GSG] */
    J[331] += dqdci;              /* dwdot[C2H4]/d[CH2GSG] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[9];
    J[387] += dqdci;              /* dwdot[H]/d[CH3] */
    J[394] -= dqdci;              /* dwdot[CH2GSG]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[401] += dqdci;              /* dwdot[C2H4]/d[CH3] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1199] -= dqdT;              /* dwdot[CH2GSG]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1206] += dqdT;              /* dwdot[C2H4]/dT */

    /*reaction 57: CH2GSG + CH4 => 2.000000 CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[15];
    k_f = prefactor_units[56] * fwd_A[56]
                * exp(fwd_beta[56] * tc[0] - activation_units[56] * fwd_Ea[56] * invT);
    dlnkfdT = fwd_beta[56] * invT + activation_units[56] * fwd_Ea[56] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[9] -= q; /* CH2GSG */
    wdot[11] += 2 * q; /* CH3 */
    wdot[15] -= q; /* CH4 */
    /* d()/d[CH2GSG] */
    dqdci =  + k_f*sc[15];
    J[324] -= dqdci;              /* dwdot[CH2GSG]/d[CH2GSG] */
    J[326] += 2 * dqdci;          /* dwdot[CH3]/d[CH2GSG] */
    J[330] -= dqdci;              /* dwdot[CH4]/d[CH2GSG] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[9];
    J[534] -= dqdci;              /* dwdot[CH2GSG]/d[CH4] */
    J[536] += 2 * dqdci;          /* dwdot[CH3]/d[CH4] */
    J[540] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[1199] -= dqdT;              /* dwdot[CH2GSG]/dT */
    J[1201] += 2 * dqdT;          /* dwdot[CH3]/dT */
    J[1205] -= dqdT;              /* dwdot[CH4]/dT */

    /*reaction 58: 2.000000 CH3 => CH2GSG + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[11], 2.000000);
    k_f = prefactor_units[57] * fwd_A[57]
                * exp(fwd_beta[57] * tc[0] - activation_units[57] * fwd_Ea[57] * invT);
    dlnkfdT = fwd_beta[57] * invT + activation_units[57] * fwd_Ea[57] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[9] += q; /* CH2GSG */
    wdot[11] -= 2 * q; /* CH3 */
    wdot[15] += q; /* CH4 */
    /* d()/d[CH3] */
    dqdci =  + k_f*2.000000*sc[11];
    J[394] += dqdci;              /* dwdot[CH2GSG]/d[CH3] */
    J[396] += -2 * dqdci;         /* dwdot[CH3]/d[CH3] */
    J[400] += dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/dT */
    J[1199] += dqdT;              /* dwdot[CH2GSG]/dT */
    J[1201] += -2 * dqdT;         /* dwdot[CH3]/dT */
    J[1205] += dqdT;              /* dwdot[CH4]/dT */

    /*reaction 59: CH4 + O => CH3 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[15];
    k_f = prefactor_units[58] * fwd_A[58]
                * exp(fwd_beta[58] * tc[0] - activation_units[58] * fwd_Ea[58] * invT);
    dlnkfdT = fwd_beta[58] * invT + activation_units[58] * fwd_Ea[58] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[3] += q; /* OH */
    wdot[11] += q; /* CH3 */
    wdot[15] -= q; /* CH4 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[15];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[38] += dqdci;               /* dwdot[OH]/d[O] */
    J[46] += dqdci;               /* dwdot[CH3]/d[O] */
    J[50] -= dqdci;               /* dwdot[CH4]/d[O] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[1];
    J[526] -= dqdci;              /* dwdot[O]/d[CH4] */
    J[528] += dqdci;              /* dwdot[OH]/d[CH4] */
    J[536] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[540] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1201] += dqdT;              /* dwdot[CH3]/dT */
    J[1205] -= dqdT;              /* dwdot[CH4]/dT */

    /*reaction 60: CH3 + OH => CH4 + O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[11];
    k_f = prefactor_units[59] * fwd_A[59]
                * exp(fwd_beta[59] * tc[0] - activation_units[59] * fwd_Ea[59] * invT);
    dlnkfdT = fwd_beta[59] * invT + activation_units[59] * fwd_Ea[59] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] += q; /* O */
    wdot[3] -= q; /* OH */
    wdot[11] -= q; /* CH3 */
    wdot[15] += q; /* CH4 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[11];
    J[106] += dqdci;              /* dwdot[O]/d[OH] */
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[116] -= dqdci;              /* dwdot[CH3]/d[OH] */
    J[120] += dqdci;              /* dwdot[CH4]/d[OH] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[3];
    J[386] += dqdci;              /* dwdot[O]/d[CH3] */
    J[388] -= dqdci;              /* dwdot[OH]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[400] += dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/dT */
    J[1191] += dqdT;              /* dwdot[O]/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1205] += dqdT;              /* dwdot[CH4]/dT */

    /*reaction 61: CH4 + H => CH3 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[15];
    k_f = prefactor_units[60] * fwd_A[60]
                * exp(fwd_beta[60] * tc[0] - activation_units[60] * fwd_Ea[60] * invT);
    dlnkfdT = fwd_beta[60] * invT + activation_units[60] * fwd_Ea[60] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[4] += q; /* H2 */
    wdot[11] += q; /* CH3 */
    wdot[15] -= q; /* CH4 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[15];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[74] += dqdci;               /* dwdot[H2]/d[H] */
    J[81] += dqdci;               /* dwdot[CH3]/d[H] */
    J[85] -= dqdci;               /* dwdot[CH4]/d[H] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[2];
    J[527] -= dqdci;              /* dwdot[H]/d[CH4] */
    J[529] += dqdci;              /* dwdot[H2]/d[CH4] */
    J[536] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[540] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */
    J[1201] += dqdT;              /* dwdot[CH3]/dT */
    J[1205] -= dqdT;              /* dwdot[CH4]/dT */

    /*reaction 62: CH3 + H2 => CH4 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[11];
    k_f = prefactor_units[61] * fwd_A[61]
                * exp(fwd_beta[61] * tc[0] - activation_units[61] * fwd_Ea[61] * invT);
    dlnkfdT = fwd_beta[61] * invT + activation_units[61] * fwd_Ea[61] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[4] -= q; /* H2 */
    wdot[11] -= q; /* CH3 */
    wdot[15] += q; /* CH4 */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[11];
    J[142] += dqdci;              /* dwdot[H]/d[H2] */
    J[144] -= dqdci;              /* dwdot[H2]/d[H2] */
    J[151] -= dqdci;              /* dwdot[CH3]/d[H2] */
    J[155] += dqdci;              /* dwdot[CH4]/d[H2] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[4];
    J[387] += dqdci;              /* dwdot[H]/d[CH3] */
    J[389] -= dqdci;              /* dwdot[H2]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[400] += dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1194] -= dqdT;              /* dwdot[H2]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1205] += dqdT;              /* dwdot[CH4]/dT */

    /*reaction 63: CH4 + OH => CH3 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[15];
    k_f = prefactor_units[62] * fwd_A[62]
                * exp(fwd_beta[62] * tc[0] - activation_units[62] * fwd_Ea[62] * invT);
    dlnkfdT = fwd_beta[62] * invT + activation_units[62] * fwd_Ea[62] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[11] += q; /* CH3 */
    wdot[15] -= q; /* CH4 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[15];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[111] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[116] += dqdci;              /* dwdot[CH3]/d[OH] */
    J[120] -= dqdci;              /* dwdot[CH4]/d[OH] */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[3];
    J[528] -= dqdci;              /* dwdot[OH]/d[CH4] */
    J[531] += dqdci;              /* dwdot[H2O]/d[CH4] */
    J[536] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[540] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */
    J[1201] += dqdT;              /* dwdot[CH3]/dT */
    J[1205] -= dqdT;              /* dwdot[CH4]/dT */

    /*reaction 64: CH3 + H2O => CH4 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[11];
    k_f = prefactor_units[63] * fwd_A[63]
                * exp(fwd_beta[63] * tc[0] - activation_units[63] * fwd_Ea[63] * invT);
    dlnkfdT = fwd_beta[63] * invT + activation_units[63] * fwd_Ea[63] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] += q; /* OH */
    wdot[6] -= q; /* H2O */
    wdot[11] -= q; /* CH3 */
    wdot[15] += q; /* CH4 */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[11];
    J[213] += dqdci;              /* dwdot[OH]/d[H2O] */
    J[216] -= dqdci;              /* dwdot[H2O]/d[H2O] */
    J[221] -= dqdci;              /* dwdot[CH3]/d[H2O] */
    J[225] += dqdci;              /* dwdot[CH4]/d[H2O] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[6];
    J[388] += dqdci;              /* dwdot[OH]/d[CH3] */
    J[391] -= dqdci;              /* dwdot[H2O]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[400] += dqdci;              /* dwdot[CH4]/d[CH3] */
    /* d()/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1196] -= dqdT;              /* dwdot[H2O]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1205] += dqdT;              /* dwdot[CH4]/dT */

    /*reaction 65: CO + HO2 => CO2 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[12];
    k_f = prefactor_units[64] * fwd_A[64]
                * exp(fwd_beta[64] * tc[0] - activation_units[64] * fwd_Ea[64] * invT);
    dlnkfdT = fwd_beta[64] * invT + activation_units[64] * fwd_Ea[64] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] += q; /* OH */
    wdot[8] -= q; /* HO2 */
    wdot[12] -= q; /* CO */
    wdot[18] += q; /* CO2 */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[12];
    J[283] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[292] -= dqdci;              /* dwdot[CO]/d[HO2] */
    J[298] += dqdci;              /* dwdot[CO2]/d[HO2] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[8];
    J[423] += dqdci;              /* dwdot[OH]/d[CO] */
    J[428] -= dqdci;              /* dwdot[HO2]/d[CO] */
    J[432] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[438] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */
    J[1202] -= dqdT;              /* dwdot[CO]/dT */
    J[1208] += dqdT;              /* dwdot[CO2]/dT */

    /*reaction 66: CO + O2 => CO2 + O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[12];
    k_f = prefactor_units[65] * fwd_A[65]
                * exp(fwd_beta[65] * tc[0] - activation_units[65] * fwd_Ea[65] * invT);
    dlnkfdT = fwd_beta[65] * invT + activation_units[65] * fwd_Ea[65] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] += q; /* O */
    wdot[5] -= q; /* O2 */
    wdot[12] -= q; /* CO */
    wdot[18] += q; /* CO2 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[12];
    J[176] += dqdci;              /* dwdot[O]/d[O2] */
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[187] -= dqdci;              /* dwdot[CO]/d[O2] */
    J[193] += dqdci;              /* dwdot[CO2]/d[O2] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[5];
    J[421] += dqdci;              /* dwdot[O]/d[CO] */
    J[425] -= dqdci;              /* dwdot[O2]/d[CO] */
    J[432] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[438] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/dT */
    J[1191] += dqdT;              /* dwdot[O]/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1202] -= dqdT;              /* dwdot[CO]/dT */
    J[1208] += dqdT;              /* dwdot[CO2]/dT */

    /*reaction 67: CO2 + O => CO + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[18];
    k_f = prefactor_units[66] * fwd_A[66]
                * exp(fwd_beta[66] * tc[0] - activation_units[66] * fwd_Ea[66] * invT);
    dlnkfdT = fwd_beta[66] * invT + activation_units[66] * fwd_Ea[66] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[5] += q; /* O2 */
    wdot[12] += q; /* CO */
    wdot[18] -= q; /* CO2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[18];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[40] += dqdci;               /* dwdot[O2]/d[O] */
    J[47] += dqdci;               /* dwdot[CO]/d[O] */
    J[53] -= dqdci;               /* dwdot[CO2]/d[O] */
    /* d()/d[CO2] */
    dqdci =  + k_f*sc[1];
    J[631] -= dqdci;              /* dwdot[O]/d[CO2] */
    J[635] += dqdci;              /* dwdot[O2]/d[CO2] */
    J[642] += dqdci;              /* dwdot[CO]/d[CO2] */
    J[648] -= dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1202] += dqdT;              /* dwdot[CO]/dT */
    J[1208] -= dqdT;              /* dwdot[CO2]/dT */

    /*reaction 68: CO + OH => CO2 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[12];
    k_f = prefactor_units[67] * fwd_A[67]
                * exp(fwd_beta[67] * tc[0] - activation_units[67] * fwd_Ea[67] * invT);
    dlnkfdT = fwd_beta[67] * invT + activation_units[67] * fwd_Ea[67] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[3] -= q; /* OH */
    wdot[12] -= q; /* CO */
    wdot[18] += q; /* CO2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[12];
    J[107] += dqdci;              /* dwdot[H]/d[OH] */
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[117] -= dqdci;              /* dwdot[CO]/d[OH] */
    J[123] += dqdci;              /* dwdot[CO2]/d[OH] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[3];
    J[422] += dqdci;              /* dwdot[H]/d[CO] */
    J[423] -= dqdci;              /* dwdot[OH]/d[CO] */
    J[432] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[438] += dqdci;              /* dwdot[CO2]/d[CO] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1202] -= dqdT;              /* dwdot[CO]/dT */
    J[1208] += dqdT;              /* dwdot[CO2]/dT */

    /*reaction 69: CO2 + H => CO + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[18];
    k_f = prefactor_units[68] * fwd_A[68]
                * exp(fwd_beta[68] * tc[0] - activation_units[68] * fwd_Ea[68] * invT);
    dlnkfdT = fwd_beta[68] * invT + activation_units[68] * fwd_Ea[68] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[3] += q; /* OH */
    wdot[12] += q; /* CO */
    wdot[18] -= q; /* CO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[18];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[73] += dqdci;               /* dwdot[OH]/d[H] */
    J[82] += dqdci;               /* dwdot[CO]/d[H] */
    J[88] -= dqdci;               /* dwdot[CO2]/d[H] */
    /* d()/d[CO2] */
    dqdci =  + k_f*sc[2];
    J[632] -= dqdci;              /* dwdot[H]/d[CO2] */
    J[633] += dqdci;              /* dwdot[OH]/d[CO2] */
    J[642] += dqdci;              /* dwdot[CO]/d[CO2] */
    J[648] -= dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1202] += dqdT;              /* dwdot[CO]/dT */
    J[1208] -= dqdT;              /* dwdot[CO2]/dT */

    /*reaction 70: HCO + CH3 => CH4 + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[19];
    k_f = prefactor_units[69] * fwd_A[69]
                * exp(fwd_beta[69] * tc[0] - activation_units[69] * fwd_Ea[69] * invT);
    dlnkfdT = fwd_beta[69] * invT + activation_units[69] * fwd_Ea[69] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[11] -= q; /* CH3 */
    wdot[12] += q; /* CO */
    wdot[15] += q; /* CH4 */
    wdot[19] -= q; /* HCO */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[19];
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[397] += dqdci;              /* dwdot[CO]/d[CH3] */
    J[400] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[404] -= dqdci;              /* dwdot[HCO]/d[CH3] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[11];
    J[676] -= dqdci;              /* dwdot[CH3]/d[HCO] */
    J[677] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[680] += dqdci;              /* dwdot[CH4]/d[HCO] */
    J[684] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1202] += dqdT;              /* dwdot[CO]/dT */
    J[1205] += dqdT;              /* dwdot[CH4]/dT */
    J[1209] -= dqdT;              /* dwdot[HCO]/dT */

    /*reaction 71: HCO + H => CO + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[19];
    k_f = prefactor_units[70] * fwd_A[70]
                * exp(fwd_beta[70] * tc[0] - activation_units[70] * fwd_Ea[70] * invT);
    dlnkfdT = fwd_beta[70] * invT + activation_units[70] * fwd_Ea[70] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[4] += q; /* H2 */
    wdot[12] += q; /* CO */
    wdot[19] -= q; /* HCO */
    /* d()/d[H] */
    dqdci =  + k_f*sc[19];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[74] += dqdci;               /* dwdot[H2]/d[H] */
    J[82] += dqdci;               /* dwdot[CO]/d[H] */
    J[89] -= dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[2];
    J[667] -= dqdci;              /* dwdot[H]/d[HCO] */
    J[669] += dqdci;              /* dwdot[H2]/d[HCO] */
    J[677] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[684] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */
    J[1202] += dqdT;              /* dwdot[CO]/dT */
    J[1209] -= dqdT;              /* dwdot[HCO]/dT */

    /*reaction 72: HCO + O2 => CO + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[19];
    k_f = prefactor_units[71] * fwd_A[71]
                * exp(fwd_beta[71] * tc[0] - activation_units[71] * fwd_Ea[71] * invT);
    dlnkfdT = fwd_beta[71] * invT + activation_units[71] * fwd_Ea[71] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[8] += q; /* HO2 */
    wdot[12] += q; /* CO */
    wdot[19] -= q; /* HCO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[19];
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[183] += dqdci;              /* dwdot[HO2]/d[O2] */
    J[187] += dqdci;              /* dwdot[CO]/d[O2] */
    J[194] -= dqdci;              /* dwdot[HCO]/d[O2] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[5];
    J[670] -= dqdci;              /* dwdot[O2]/d[HCO] */
    J[673] += dqdci;              /* dwdot[HO2]/d[HCO] */
    J[677] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[684] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */
    J[1202] += dqdT;              /* dwdot[CO]/dT */
    J[1209] -= dqdT;              /* dwdot[HCO]/dT */

    /*reaction 73: HCO + O => CO + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[19];
    k_f = prefactor_units[72] * fwd_A[72]
                * exp(fwd_beta[72] * tc[0] - activation_units[72] * fwd_Ea[72] * invT);
    dlnkfdT = fwd_beta[72] * invT + activation_units[72] * fwd_Ea[72] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[3] += q; /* OH */
    wdot[12] += q; /* CO */
    wdot[19] -= q; /* HCO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[19];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[38] += dqdci;               /* dwdot[OH]/d[O] */
    J[47] += dqdci;               /* dwdot[CO]/d[O] */
    J[54] -= dqdci;               /* dwdot[HCO]/d[O] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[1];
    J[666] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[668] += dqdci;              /* dwdot[OH]/d[HCO] */
    J[677] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[684] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1202] += dqdT;              /* dwdot[CO]/dT */
    J[1209] -= dqdT;              /* dwdot[HCO]/dT */

    /*reaction 74: HCO + O => CO2 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[19];
    k_f = prefactor_units[73] * fwd_A[73]
                * exp(fwd_beta[73] * tc[0] - activation_units[73] * fwd_Ea[73] * invT);
    dlnkfdT = fwd_beta[73] * invT + activation_units[73] * fwd_Ea[73] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[2] += q; /* H */
    wdot[18] += q; /* CO2 */
    wdot[19] -= q; /* HCO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[19];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[37] += dqdci;               /* dwdot[H]/d[O] */
    J[53] += dqdci;               /* dwdot[CO2]/d[O] */
    J[54] -= dqdci;               /* dwdot[HCO]/d[O] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[1];
    J[666] -= dqdci;              /* dwdot[O]/d[HCO] */
    J[667] += dqdci;              /* dwdot[H]/d[HCO] */
    J[683] += dqdci;              /* dwdot[CO2]/d[HCO] */
    J[684] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1208] += dqdT;              /* dwdot[CO2]/dT */
    J[1209] -= dqdT;              /* dwdot[HCO]/dT */

    /*reaction 75: HCO + OH => CO + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[19];
    k_f = prefactor_units[74] * fwd_A[74]
                * exp(fwd_beta[74] * tc[0] - activation_units[74] * fwd_Ea[74] * invT);
    dlnkfdT = fwd_beta[74] * invT + activation_units[74] * fwd_Ea[74] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[12] += q; /* CO */
    wdot[19] -= q; /* HCO */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[19];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[111] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[117] += dqdci;              /* dwdot[CO]/d[OH] */
    J[124] -= dqdci;              /* dwdot[HCO]/d[OH] */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[3];
    J[668] -= dqdci;              /* dwdot[OH]/d[HCO] */
    J[671] += dqdci;              /* dwdot[H2O]/d[HCO] */
    J[677] += dqdci;              /* dwdot[CO]/d[HCO] */
    J[684] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */
    J[1202] += dqdT;              /* dwdot[CO]/dT */
    J[1209] -= dqdT;              /* dwdot[HCO]/dT */

    /*reaction 76: CH2O + H => HCO + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[10];
    k_f = prefactor_units[75] * fwd_A[75]
                * exp(fwd_beta[75] * tc[0] - activation_units[75] * fwd_Ea[75] * invT);
    dlnkfdT = fwd_beta[75] * invT + activation_units[75] * fwd_Ea[75] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[4] += q; /* H2 */
    wdot[10] -= q; /* CH2O */
    wdot[19] += q; /* HCO */
    /* d()/d[H] */
    dqdci =  + k_f*sc[10];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[74] += dqdci;               /* dwdot[H2]/d[H] */
    J[80] -= dqdci;               /* dwdot[CH2O]/d[H] */
    J[89] += dqdci;               /* dwdot[HCO]/d[H] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[2];
    J[352] -= dqdci;              /* dwdot[H]/d[CH2O] */
    J[354] += dqdci;              /* dwdot[H2]/d[CH2O] */
    J[360] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[369] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */
    J[1200] -= dqdT;              /* dwdot[CH2O]/dT */
    J[1209] += dqdT;              /* dwdot[HCO]/dT */

    /*reaction 77: CH2O + O2 => HCO + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[10];
    k_f = prefactor_units[76] * fwd_A[76]
                * exp(fwd_beta[76] * tc[0] - activation_units[76] * fwd_Ea[76] * invT);
    dlnkfdT = fwd_beta[76] * invT + activation_units[76] * fwd_Ea[76] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[8] += q; /* HO2 */
    wdot[10] -= q; /* CH2O */
    wdot[19] += q; /* HCO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[10];
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[183] += dqdci;              /* dwdot[HO2]/d[O2] */
    J[185] -= dqdci;              /* dwdot[CH2O]/d[O2] */
    J[194] += dqdci;              /* dwdot[HCO]/d[O2] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[5];
    J[355] -= dqdci;              /* dwdot[O2]/d[CH2O] */
    J[358] += dqdci;              /* dwdot[HO2]/d[CH2O] */
    J[360] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[369] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    /* d()/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */
    J[1200] -= dqdT;              /* dwdot[CH2O]/dT */
    J[1209] += dqdT;              /* dwdot[HCO]/dT */

    /*reaction 78: CH2O + OH => HCO + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[10];
    k_f = prefactor_units[77] * fwd_A[77]
                * exp(fwd_beta[77] * tc[0] - activation_units[77] * fwd_Ea[77] * invT);
    dlnkfdT = fwd_beta[77] * invT + activation_units[77] * fwd_Ea[77] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[10] -= q; /* CH2O */
    wdot[19] += q; /* HCO */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[10];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[111] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[115] -= dqdci;              /* dwdot[CH2O]/d[OH] */
    J[124] += dqdci;              /* dwdot[HCO]/d[OH] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[3];
    J[353] -= dqdci;              /* dwdot[OH]/d[CH2O] */
    J[356] += dqdci;              /* dwdot[H2O]/d[CH2O] */
    J[360] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[369] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */
    J[1200] -= dqdT;              /* dwdot[CH2O]/dT */
    J[1209] += dqdT;              /* dwdot[HCO]/dT */

    /*reaction 79: CH2O + HO2 => HCO + H2O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[10];
    k_f = prefactor_units[78] * fwd_A[78]
                * exp(fwd_beta[78] * tc[0] - activation_units[78] * fwd_Ea[78] * invT);
    dlnkfdT = fwd_beta[78] * invT + activation_units[78] * fwd_Ea[78] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[7] += q; /* H2O2 */
    wdot[8] -= q; /* HO2 */
    wdot[10] -= q; /* CH2O */
    wdot[19] += q; /* HCO */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[10];
    J[287] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[290] -= dqdci;              /* dwdot[CH2O]/d[HO2] */
    J[299] += dqdci;              /* dwdot[HCO]/d[HO2] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[8];
    J[357] += dqdci;              /* dwdot[H2O2]/d[CH2O] */
    J[358] -= dqdci;              /* dwdot[HO2]/d[CH2O] */
    J[360] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[369] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    /* d()/dT */
    J[1197] += dqdT;              /* dwdot[H2O2]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */
    J[1200] -= dqdT;              /* dwdot[CH2O]/dT */
    J[1209] += dqdT;              /* dwdot[HCO]/dT */

    /*reaction 80: CH2O + O => HCO + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[10];
    k_f = prefactor_units[79] * fwd_A[79]
                * exp(fwd_beta[79] * tc[0] - activation_units[79] * fwd_Ea[79] * invT);
    dlnkfdT = fwd_beta[79] * invT + activation_units[79] * fwd_Ea[79] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[3] += q; /* OH */
    wdot[10] -= q; /* CH2O */
    wdot[19] += q; /* HCO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[10];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[38] += dqdci;               /* dwdot[OH]/d[O] */
    J[45] -= dqdci;               /* dwdot[CH2O]/d[O] */
    J[54] += dqdci;               /* dwdot[HCO]/d[O] */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[1];
    J[351] -= dqdci;              /* dwdot[O]/d[CH2O] */
    J[353] += dqdci;              /* dwdot[OH]/d[CH2O] */
    J[360] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[369] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1200] -= dqdT;              /* dwdot[CH2O]/dT */
    J[1209] += dqdT;              /* dwdot[HCO]/dT */

    /*reaction 81: CH2O + CH3 => HCO + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[10]*sc[11];
    k_f = prefactor_units[80] * fwd_A[80]
                * exp(fwd_beta[80] * tc[0] - activation_units[80] * fwd_Ea[80] * invT);
    dlnkfdT = fwd_beta[80] * invT + activation_units[80] * fwd_Ea[80] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[10] -= q; /* CH2O */
    wdot[11] -= q; /* CH3 */
    wdot[15] += q; /* CH4 */
    wdot[19] += q; /* HCO */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[11];
    J[360] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[361] -= dqdci;              /* dwdot[CH3]/d[CH2O] */
    J[365] += dqdci;              /* dwdot[CH4]/d[CH2O] */
    J[369] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[10];
    J[395] -= dqdci;              /* dwdot[CH2O]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[400] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[404] += dqdci;              /* dwdot[HCO]/d[CH3] */
    /* d()/dT */
    J[1200] -= dqdT;              /* dwdot[CH2O]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1205] += dqdT;              /* dwdot[CH4]/dT */
    J[1209] += dqdT;              /* dwdot[HCO]/dT */

    /*reaction 82: CH3O + O2 => CH2O + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[13];
    k_f = prefactor_units[81] * fwd_A[81]
                * exp(fwd_beta[81] * tc[0] - activation_units[81] * fwd_Ea[81] * invT);
    dlnkfdT = fwd_beta[81] * invT + activation_units[81] * fwd_Ea[81] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[8] += q; /* HO2 */
    wdot[10] += q; /* CH2O */
    wdot[13] -= q; /* CH3O */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[13];
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[183] += dqdci;              /* dwdot[HO2]/d[O2] */
    J[185] += dqdci;              /* dwdot[CH2O]/d[O2] */
    J[188] -= dqdci;              /* dwdot[CH3O]/d[O2] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[5];
    J[460] -= dqdci;              /* dwdot[O2]/d[CH3O] */
    J[463] += dqdci;              /* dwdot[HO2]/d[CH3O] */
    J[465] += dqdci;              /* dwdot[CH2O]/d[CH3O] */
    J[468] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    /* d()/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */
    J[1200] += dqdT;              /* dwdot[CH2O]/dT */
    J[1203] -= dqdT;              /* dwdot[CH3O]/dT */

    /*reaction 83: CH2GSG + CO2 => CH2O + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[18];
    k_f = prefactor_units[82] * fwd_A[82]
                * exp(fwd_beta[82] * tc[0] - activation_units[82] * fwd_Ea[82] * invT);
    dlnkfdT = fwd_beta[82] * invT + activation_units[82] * fwd_Ea[82] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[9] -= q; /* CH2GSG */
    wdot[10] += q; /* CH2O */
    wdot[12] += q; /* CO */
    wdot[18] -= q; /* CO2 */
    /* d()/d[CH2GSG] */
    dqdci =  + k_f*sc[18];
    J[324] -= dqdci;              /* dwdot[CH2GSG]/d[CH2GSG] */
    J[325] += dqdci;              /* dwdot[CH2O]/d[CH2GSG] */
    J[327] += dqdci;              /* dwdot[CO]/d[CH2GSG] */
    J[333] -= dqdci;              /* dwdot[CO2]/d[CH2GSG] */
    /* d()/d[CO2] */
    dqdci =  + k_f*sc[9];
    J[639] -= dqdci;              /* dwdot[CH2GSG]/d[CO2] */
    J[640] += dqdci;              /* dwdot[CH2O]/d[CO2] */
    J[642] += dqdci;              /* dwdot[CO]/d[CO2] */
    J[648] -= dqdci;              /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[1199] -= dqdT;              /* dwdot[CH2GSG]/dT */
    J[1200] += dqdT;              /* dwdot[CH2O]/dT */
    J[1202] += dqdT;              /* dwdot[CO]/dT */
    J[1208] -= dqdT;              /* dwdot[CO2]/dT */

    /*reaction 84: CH3O2 + CH3 => 2.000000 CH3O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[20];
    k_f = prefactor_units[83] * fwd_A[83]
                * exp(fwd_beta[83] * tc[0] - activation_units[83] * fwd_Ea[83] * invT);
    dlnkfdT = fwd_beta[83] * invT + activation_units[83] * fwd_Ea[83] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[11] -= q; /* CH3 */
    wdot[13] += 2 * q; /* CH3O */
    wdot[20] -= q; /* CH3O2 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[20];
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[398] += 2 * dqdci;          /* dwdot[CH3O]/d[CH3] */
    J[405] -= dqdci;              /* dwdot[CH3O2]/d[CH3] */
    /* d()/d[CH3O2] */
    dqdci =  + k_f*sc[11];
    J[711] -= dqdci;              /* dwdot[CH3]/d[CH3O2] */
    J[713] += 2 * dqdci;          /* dwdot[CH3O]/d[CH3O2] */
    J[720] -= dqdci;              /* dwdot[CH3O2]/d[CH3O2] */
    /* d()/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1203] += 2 * dqdT;          /* dwdot[CH3O]/dT */
    J[1210] -= dqdT;              /* dwdot[CH3O2]/dT */

    /*reaction 85: CH3O2 + HO2 => CH3O2H + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[20];
    k_f = prefactor_units[84] * fwd_A[84]
                * exp(fwd_beta[84] * tc[0] - activation_units[84] * fwd_Ea[84] * invT);
    dlnkfdT = fwd_beta[84] * invT + activation_units[84] * fwd_Ea[84] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[8] -= q; /* HO2 */
    wdot[20] -= q; /* CH3O2 */
    wdot[21] += q; /* CH3O2H */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[20];
    J[285] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[300] -= dqdci;              /* dwdot[CH3O2]/d[HO2] */
    J[301] += dqdci;              /* dwdot[CH3O2H]/d[HO2] */
    /* d()/d[CH3O2] */
    dqdci =  + k_f*sc[8];
    J[705] += dqdci;              /* dwdot[O2]/d[CH3O2] */
    J[708] -= dqdci;              /* dwdot[HO2]/d[CH3O2] */
    J[720] -= dqdci;              /* dwdot[CH3O2]/d[CH3O2] */
    J[721] += dqdci;              /* dwdot[CH3O2H]/d[CH3O2] */
    /* d()/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */
    J[1210] -= dqdT;              /* dwdot[CH3O2]/dT */
    J[1211] += dqdT;              /* dwdot[CH3O2H]/dT */

    /*reaction 86: 2.000000 CH3O2 => O2 + 2.000000 CH3O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[20], 2.000000);
    k_f = prefactor_units[85] * fwd_A[85]
                * exp(fwd_beta[85] * tc[0] - activation_units[85] * fwd_Ea[85] * invT);
    dlnkfdT = fwd_beta[85] * invT + activation_units[85] * fwd_Ea[85] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[13] += 2 * q; /* CH3O */
    wdot[20] -= 2 * q; /* CH3O2 */
    /* d()/d[CH3O2] */
    dqdci =  + k_f*2.000000*sc[20];
    J[705] += dqdci;              /* dwdot[O2]/d[CH3O2] */
    J[713] += 2 * dqdci;          /* dwdot[CH3O]/d[CH3O2] */
    J[720] += -2 * dqdci;         /* dwdot[CH3O2]/d[CH3O2] */
    /* d()/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1203] += 2 * dqdT;          /* dwdot[CH3O]/dT */
    J[1210] += -2 * dqdT;         /* dwdot[CH3O2]/dT */

    /*reaction 87: CH3O2 + CH2O => CH3O2H + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[10]*sc[20];
    k_f = prefactor_units[86] * fwd_A[86]
                * exp(fwd_beta[86] * tc[0] - activation_units[86] * fwd_Ea[86] * invT);
    dlnkfdT = fwd_beta[86] * invT + activation_units[86] * fwd_Ea[86] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[10] -= q; /* CH2O */
    wdot[19] += q; /* HCO */
    wdot[20] -= q; /* CH3O2 */
    wdot[21] += q; /* CH3O2H */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[20];
    J[360] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[369] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[370] -= dqdci;              /* dwdot[CH3O2]/d[CH2O] */
    J[371] += dqdci;              /* dwdot[CH3O2H]/d[CH2O] */
    /* d()/d[CH3O2] */
    dqdci =  + k_f*sc[10];
    J[710] -= dqdci;              /* dwdot[CH2O]/d[CH3O2] */
    J[719] += dqdci;              /* dwdot[HCO]/d[CH3O2] */
    J[720] -= dqdci;              /* dwdot[CH3O2]/d[CH3O2] */
    J[721] += dqdci;              /* dwdot[CH3O2H]/d[CH3O2] */
    /* d()/dT */
    J[1200] -= dqdT;              /* dwdot[CH2O]/dT */
    J[1209] += dqdT;              /* dwdot[HCO]/dT */
    J[1210] -= dqdT;              /* dwdot[CH3O2]/dT */
    J[1211] += dqdT;              /* dwdot[CH3O2H]/dT */

    /*reaction 88: CH3O2H => CH3O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[21];
    k_f = prefactor_units[87] * fwd_A[87]
                * exp(fwd_beta[87] * tc[0] - activation_units[87] * fwd_Ea[87] * invT);
    dlnkfdT = fwd_beta[87] * invT + activation_units[87] * fwd_Ea[87] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] += q; /* OH */
    wdot[13] += q; /* CH3O */
    wdot[21] -= q; /* CH3O2H */
    /* d()/d[CH3O2H] */
    dqdci =  + k_f;
    J[738] += dqdci;              /* dwdot[OH]/d[CH3O2H] */
    J[748] += dqdci;              /* dwdot[CH3O]/d[CH3O2H] */
    J[756] -= dqdci;              /* dwdot[CH3O2H]/d[CH3O2H] */
    /* d()/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1203] += dqdT;              /* dwdot[CH3O]/dT */
    J[1211] -= dqdT;              /* dwdot[CH3O2H]/dT */

    /*reaction 89: CH3O + OH => CH3O2H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[13];
    k_f = prefactor_units[88] * fwd_A[88]
                * exp(fwd_beta[88] * tc[0] - activation_units[88] * fwd_Ea[88] * invT);
    dlnkfdT = fwd_beta[88] * invT + activation_units[88] * fwd_Ea[88] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[13] -= q; /* CH3O */
    wdot[21] += q; /* CH3O2H */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[13];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[118] -= dqdci;              /* dwdot[CH3O]/d[OH] */
    J[126] += dqdci;              /* dwdot[CH3O2H]/d[OH] */
    /* d()/d[CH3O] */
    dqdci =  + k_f*sc[3];
    J[458] -= dqdci;              /* dwdot[OH]/d[CH3O] */
    J[468] -= dqdci;              /* dwdot[CH3O]/d[CH3O] */
    J[476] += dqdci;              /* dwdot[CH3O2H]/d[CH3O] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1203] -= dqdT;              /* dwdot[CH3O]/dT */
    J[1211] += dqdT;              /* dwdot[CH3O2H]/dT */

    /*reaction 90: C2H2 + O2 => HCCO + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[22];
    k_f = prefactor_units[89] * fwd_A[89]
                * exp(fwd_beta[89] * tc[0] - activation_units[89] * fwd_Ea[89] * invT);
    dlnkfdT = fwd_beta[89] * invT + activation_units[89] * fwd_Ea[89] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] += q; /* OH */
    wdot[5] -= q; /* O2 */
    wdot[22] -= q; /* C2H2 */
    wdot[23] += q; /* HCCO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[22];
    J[178] += dqdci;              /* dwdot[OH]/d[O2] */
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[197] -= dqdci;              /* dwdot[C2H2]/d[O2] */
    J[198] += dqdci;              /* dwdot[HCCO]/d[O2] */
    /* d()/d[C2H2] */
    dqdci =  + k_f*sc[5];
    J[773] += dqdci;              /* dwdot[OH]/d[C2H2] */
    J[775] -= dqdci;              /* dwdot[O2]/d[C2H2] */
    J[792] -= dqdci;              /* dwdot[C2H2]/d[C2H2] */
    J[793] += dqdci;              /* dwdot[HCCO]/d[C2H2] */
    /* d()/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1212] -= dqdT;              /* dwdot[C2H2]/dT */
    J[1213] += dqdT;              /* dwdot[HCCO]/dT */

    /*reaction 91: C2H2 + O => HCCO + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[22];
    k_f = prefactor_units[90] * fwd_A[90]
                * exp(fwd_beta[90] * tc[0] - activation_units[90] * fwd_Ea[90] * invT);
    dlnkfdT = fwd_beta[90] * invT + activation_units[90] * fwd_Ea[90] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[2] += q; /* H */
    wdot[22] -= q; /* C2H2 */
    wdot[23] += q; /* HCCO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[22];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[37] += dqdci;               /* dwdot[H]/d[O] */
    J[57] -= dqdci;               /* dwdot[C2H2]/d[O] */
    J[58] += dqdci;               /* dwdot[HCCO]/d[O] */
    /* d()/d[C2H2] */
    dqdci =  + k_f*sc[1];
    J[771] -= dqdci;              /* dwdot[O]/d[C2H2] */
    J[772] += dqdci;              /* dwdot[H]/d[C2H2] */
    J[792] -= dqdci;              /* dwdot[C2H2]/d[C2H2] */
    J[793] += dqdci;              /* dwdot[HCCO]/d[C2H2] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1212] -= dqdT;              /* dwdot[C2H2]/dT */
    J[1213] += dqdT;              /* dwdot[HCCO]/dT */

    /*reaction 92: C2H3 + O2 => CH2CHO + O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[24];
    k_f = prefactor_units[91] * fwd_A[91]
                * exp(fwd_beta[91] * tc[0] - activation_units[91] * fwd_Ea[91] * invT);
    dlnkfdT = fwd_beta[91] * invT + activation_units[91] * fwd_Ea[91] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] += q; /* O */
    wdot[5] -= q; /* O2 */
    wdot[24] -= q; /* C2H3 */
    wdot[25] += q; /* CH2CHO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[24];
    J[176] += dqdci;              /* dwdot[O]/d[O2] */
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[199] -= dqdci;              /* dwdot[C2H3]/d[O2] */
    J[200] += dqdci;              /* dwdot[CH2CHO]/d[O2] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[5];
    J[841] += dqdci;              /* dwdot[O]/d[C2H3] */
    J[845] -= dqdci;              /* dwdot[O2]/d[C2H3] */
    J[864] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    J[865] += dqdci;              /* dwdot[CH2CHO]/d[C2H3] */
    /* d()/dT */
    J[1191] += dqdT;              /* dwdot[O]/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1214] -= dqdT;              /* dwdot[C2H3]/dT */
    J[1215] += dqdT;              /* dwdot[CH2CHO]/dT */

    /*reaction 93: C2H3 + H => C2H2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[24];
    k_f = prefactor_units[92] * fwd_A[92]
                * exp(fwd_beta[92] * tc[0] - activation_units[92] * fwd_Ea[92] * invT);
    dlnkfdT = fwd_beta[92] * invT + activation_units[92] * fwd_Ea[92] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[4] += q; /* H2 */
    wdot[22] += q; /* C2H2 */
    wdot[24] -= q; /* C2H3 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[24];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[74] += dqdci;               /* dwdot[H2]/d[H] */
    J[92] += dqdci;               /* dwdot[C2H2]/d[H] */
    J[94] -= dqdci;               /* dwdot[C2H3]/d[H] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[2];
    J[842] -= dqdci;              /* dwdot[H]/d[C2H3] */
    J[844] += dqdci;              /* dwdot[H2]/d[C2H3] */
    J[862] += dqdci;              /* dwdot[C2H2]/d[C2H3] */
    J[864] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */
    J[1212] += dqdT;              /* dwdot[C2H2]/dT */
    J[1214] -= dqdT;              /* dwdot[C2H3]/dT */

    /*reaction 94: C2H3 + O2 => CH2O + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[24];
    k_f = prefactor_units[93] * fwd_A[93]
                * exp(fwd_beta[93] * tc[0] - activation_units[93] * fwd_Ea[93] * invT);
    dlnkfdT = fwd_beta[93] * invT + activation_units[93] * fwd_Ea[93] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[10] += q; /* CH2O */
    wdot[19] += q; /* HCO */
    wdot[24] -= q; /* C2H3 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[24];
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[185] += dqdci;              /* dwdot[CH2O]/d[O2] */
    J[194] += dqdci;              /* dwdot[HCO]/d[O2] */
    J[199] -= dqdci;              /* dwdot[C2H3]/d[O2] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[5];
    J[845] -= dqdci;              /* dwdot[O2]/d[C2H3] */
    J[850] += dqdci;              /* dwdot[CH2O]/d[C2H3] */
    J[859] += dqdci;              /* dwdot[HCO]/d[C2H3] */
    J[864] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    /* d()/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1200] += dqdT;              /* dwdot[CH2O]/dT */
    J[1209] += dqdT;              /* dwdot[HCO]/dT */
    J[1214] -= dqdT;              /* dwdot[C2H3]/dT */

    /*reaction 95: C2H3 + O2 => C2H2 + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[24];
    k_f = prefactor_units[94] * fwd_A[94]
                * exp(fwd_beta[94] * tc[0] - activation_units[94] * fwd_Ea[94] * invT);
    dlnkfdT = fwd_beta[94] * invT + activation_units[94] * fwd_Ea[94] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[8] += q; /* HO2 */
    wdot[22] += q; /* C2H2 */
    wdot[24] -= q; /* C2H3 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[24];
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[183] += dqdci;              /* dwdot[HO2]/d[O2] */
    J[197] += dqdci;              /* dwdot[C2H2]/d[O2] */
    J[199] -= dqdci;              /* dwdot[C2H3]/d[O2] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[5];
    J[845] -= dqdci;              /* dwdot[O2]/d[C2H3] */
    J[848] += dqdci;              /* dwdot[HO2]/d[C2H3] */
    J[862] += dqdci;              /* dwdot[C2H2]/d[C2H3] */
    J[864] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    /* d()/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */
    J[1212] += dqdT;              /* dwdot[C2H2]/dT */
    J[1214] -= dqdT;              /* dwdot[C2H3]/dT */

    /*reaction 96: C2H4 + O => CH3 + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[16];
    k_f = prefactor_units[95] * fwd_A[95]
                * exp(fwd_beta[95] * tc[0] - activation_units[95] * fwd_Ea[95] * invT);
    dlnkfdT = fwd_beta[95] * invT + activation_units[95] * fwd_Ea[95] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[11] += q; /* CH3 */
    wdot[16] -= q; /* C2H4 */
    wdot[19] += q; /* HCO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[16];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[46] += dqdci;               /* dwdot[CH3]/d[O] */
    J[51] -= dqdci;               /* dwdot[C2H4]/d[O] */
    J[54] += dqdci;               /* dwdot[HCO]/d[O] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[1];
    J[561] -= dqdci;              /* dwdot[O]/d[C2H4] */
    J[571] += dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[576] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[579] += dqdci;              /* dwdot[HCO]/d[C2H4] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1201] += dqdT;              /* dwdot[CH3]/dT */
    J[1206] -= dqdT;              /* dwdot[C2H4]/dT */
    J[1209] += dqdT;              /* dwdot[HCO]/dT */

    /*reaction 97: C2H4 + OH => C2H3 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[16];
    k_f = prefactor_units[96] * fwd_A[96]
                * exp(fwd_beta[96] * tc[0] - activation_units[96] * fwd_Ea[96] * invT);
    dlnkfdT = fwd_beta[96] * invT + activation_units[96] * fwd_Ea[96] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[16] -= q; /* C2H4 */
    wdot[24] += q; /* C2H3 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[16];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[111] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[121] -= dqdci;              /* dwdot[C2H4]/d[OH] */
    J[129] += dqdci;              /* dwdot[C2H3]/d[OH] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[3];
    J[563] -= dqdci;              /* dwdot[OH]/d[C2H4] */
    J[566] += dqdci;              /* dwdot[H2O]/d[C2H4] */
    J[576] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[584] += dqdci;              /* dwdot[C2H3]/d[C2H4] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */
    J[1206] -= dqdT;              /* dwdot[C2H4]/dT */
    J[1214] += dqdT;              /* dwdot[C2H3]/dT */

    /*reaction 98: C2H3 + H2O => C2H4 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[6]*sc[24];
    k_f = prefactor_units[97] * fwd_A[97]
                * exp(fwd_beta[97] * tc[0] - activation_units[97] * fwd_Ea[97] * invT);
    dlnkfdT = fwd_beta[97] * invT + activation_units[97] * fwd_Ea[97] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] += q; /* OH */
    wdot[6] -= q; /* H2O */
    wdot[16] += q; /* C2H4 */
    wdot[24] -= q; /* C2H3 */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[24];
    J[213] += dqdci;              /* dwdot[OH]/d[H2O] */
    J[216] -= dqdci;              /* dwdot[H2O]/d[H2O] */
    J[226] += dqdci;              /* dwdot[C2H4]/d[H2O] */
    J[234] -= dqdci;              /* dwdot[C2H3]/d[H2O] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[6];
    J[843] += dqdci;              /* dwdot[OH]/d[C2H3] */
    J[846] -= dqdci;              /* dwdot[H2O]/d[C2H3] */
    J[856] += dqdci;              /* dwdot[C2H4]/d[C2H3] */
    J[864] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    /* d()/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1196] -= dqdT;              /* dwdot[H2O]/dT */
    J[1206] += dqdT;              /* dwdot[C2H4]/dT */
    J[1214] -= dqdT;              /* dwdot[C2H3]/dT */

    /*reaction 99: C2H4 + H => C2H3 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[16];
    k_f = prefactor_units[98] * fwd_A[98]
                * exp(fwd_beta[98] * tc[0] - activation_units[98] * fwd_Ea[98] * invT);
    dlnkfdT = fwd_beta[98] * invT + activation_units[98] * fwd_Ea[98] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[4] += q; /* H2 */
    wdot[16] -= q; /* C2H4 */
    wdot[24] += q; /* C2H3 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[16];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[74] += dqdci;               /* dwdot[H2]/d[H] */
    J[86] -= dqdci;               /* dwdot[C2H4]/d[H] */
    J[94] += dqdci;               /* dwdot[C2H3]/d[H] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[2];
    J[562] -= dqdci;              /* dwdot[H]/d[C2H4] */
    J[564] += dqdci;              /* dwdot[H2]/d[C2H4] */
    J[576] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[584] += dqdci;              /* dwdot[C2H3]/d[C2H4] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */
    J[1206] -= dqdT;              /* dwdot[C2H4]/dT */
    J[1214] += dqdT;              /* dwdot[C2H3]/dT */

    /*reaction 100: C2H3 + H2 => C2H4 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[24];
    k_f = prefactor_units[99] * fwd_A[99]
                * exp(fwd_beta[99] * tc[0] - activation_units[99] * fwd_Ea[99] * invT);
    dlnkfdT = fwd_beta[99] * invT + activation_units[99] * fwd_Ea[99] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[4] -= q; /* H2 */
    wdot[16] += q; /* C2H4 */
    wdot[24] -= q; /* C2H3 */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[24];
    J[142] += dqdci;              /* dwdot[H]/d[H2] */
    J[144] -= dqdci;              /* dwdot[H2]/d[H2] */
    J[156] += dqdci;              /* dwdot[C2H4]/d[H2] */
    J[164] -= dqdci;              /* dwdot[C2H3]/d[H2] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[4];
    J[842] += dqdci;              /* dwdot[H]/d[C2H3] */
    J[844] -= dqdci;              /* dwdot[H2]/d[C2H3] */
    J[856] += dqdci;              /* dwdot[C2H4]/d[C2H3] */
    J[864] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1194] -= dqdT;              /* dwdot[H2]/dT */
    J[1206] += dqdT;              /* dwdot[C2H4]/dT */
    J[1214] -= dqdT;              /* dwdot[C2H3]/dT */

    /*reaction 101: C2H4 + CH3O2 => C2H3 + CH3O2H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[16]*sc[20];
    k_f = prefactor_units[100] * fwd_A[100]
                * exp(fwd_beta[100] * tc[0] - activation_units[100] * fwd_Ea[100] * invT);
    dlnkfdT = fwd_beta[100] * invT + activation_units[100] * fwd_Ea[100] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[16] -= q; /* C2H4 */
    wdot[20] -= q; /* CH3O2 */
    wdot[21] += q; /* CH3O2H */
    wdot[24] += q; /* C2H3 */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[20];
    J[576] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[580] -= dqdci;              /* dwdot[CH3O2]/d[C2H4] */
    J[581] += dqdci;              /* dwdot[CH3O2H]/d[C2H4] */
    J[584] += dqdci;              /* dwdot[C2H3]/d[C2H4] */
    /* d()/d[CH3O2] */
    dqdci =  + k_f*sc[16];
    J[716] -= dqdci;              /* dwdot[C2H4]/d[CH3O2] */
    J[720] -= dqdci;              /* dwdot[CH3O2]/d[CH3O2] */
    J[721] += dqdci;              /* dwdot[CH3O2H]/d[CH3O2] */
    J[724] += dqdci;              /* dwdot[C2H3]/d[CH3O2] */
    /* d()/dT */
    J[1206] -= dqdT;              /* dwdot[C2H4]/dT */
    J[1210] -= dqdT;              /* dwdot[CH3O2]/dT */
    J[1211] += dqdT;              /* dwdot[CH3O2H]/dT */
    J[1214] += dqdT;              /* dwdot[C2H3]/dT */

    /*reaction 102: C2H4 + CH3 => C2H3 + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[16];
    k_f = prefactor_units[101] * fwd_A[101]
                * exp(fwd_beta[101] * tc[0] - activation_units[101] * fwd_Ea[101] * invT);
    dlnkfdT = fwd_beta[101] * invT + activation_units[101] * fwd_Ea[101] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[11] -= q; /* CH3 */
    wdot[15] += q; /* CH4 */
    wdot[16] -= q; /* C2H4 */
    wdot[24] += q; /* C2H3 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[16];
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[400] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[401] -= dqdci;              /* dwdot[C2H4]/d[CH3] */
    J[409] += dqdci;              /* dwdot[C2H3]/d[CH3] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[11];
    J[571] -= dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[575] += dqdci;              /* dwdot[CH4]/d[C2H4] */
    J[576] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[584] += dqdci;              /* dwdot[C2H3]/d[C2H4] */
    /* d()/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1205] += dqdT;              /* dwdot[CH4]/dT */
    J[1206] -= dqdT;              /* dwdot[C2H4]/dT */
    J[1214] += dqdT;              /* dwdot[C2H3]/dT */

    /*reaction 103: C2H3 + CH4 => C2H4 + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[15]*sc[24];
    k_f = prefactor_units[102] * fwd_A[102]
                * exp(fwd_beta[102] * tc[0] - activation_units[102] * fwd_Ea[102] * invT);
    dlnkfdT = fwd_beta[102] * invT + activation_units[102] * fwd_Ea[102] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[11] += q; /* CH3 */
    wdot[15] -= q; /* CH4 */
    wdot[16] += q; /* C2H4 */
    wdot[24] -= q; /* C2H3 */
    /* d()/d[CH4] */
    dqdci =  + k_f*sc[24];
    J[536] += dqdci;              /* dwdot[CH3]/d[CH4] */
    J[540] -= dqdci;              /* dwdot[CH4]/d[CH4] */
    J[541] += dqdci;              /* dwdot[C2H4]/d[CH4] */
    J[549] -= dqdci;              /* dwdot[C2H3]/d[CH4] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[15];
    J[851] += dqdci;              /* dwdot[CH3]/d[C2H3] */
    J[855] -= dqdci;              /* dwdot[CH4]/d[C2H3] */
    J[856] += dqdci;              /* dwdot[C2H4]/d[C2H3] */
    J[864] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    /* d()/dT */
    J[1201] += dqdT;              /* dwdot[CH3]/dT */
    J[1205] -= dqdT;              /* dwdot[CH4]/dT */
    J[1206] += dqdT;              /* dwdot[C2H4]/dT */
    J[1214] -= dqdT;              /* dwdot[C2H3]/dT */

    /*reaction 104: C2H4 + O => CH2CHO + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[16];
    k_f = prefactor_units[103] * fwd_A[103]
                * exp(fwd_beta[103] * tc[0] - activation_units[103] * fwd_Ea[103] * invT);
    dlnkfdT = fwd_beta[103] * invT + activation_units[103] * fwd_Ea[103] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[2] += q; /* H */
    wdot[16] -= q; /* C2H4 */
    wdot[25] += q; /* CH2CHO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[16];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[37] += dqdci;               /* dwdot[H]/d[O] */
    J[51] -= dqdci;               /* dwdot[C2H4]/d[O] */
    J[60] += dqdci;               /* dwdot[CH2CHO]/d[O] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[1];
    J[561] -= dqdci;              /* dwdot[O]/d[C2H4] */
    J[562] += dqdci;              /* dwdot[H]/d[C2H4] */
    J[576] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[585] += dqdci;              /* dwdot[CH2CHO]/d[C2H4] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1206] -= dqdT;              /* dwdot[C2H4]/dT */
    J[1215] += dqdT;              /* dwdot[CH2CHO]/dT */

    /*reaction 105: CH3 + C2H5 => CH4 + C2H4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[14];
    k_f = prefactor_units[104] * fwd_A[104]
                * exp(fwd_beta[104] * tc[0] - activation_units[104] * fwd_Ea[104] * invT);
    dlnkfdT = fwd_beta[104] * invT + activation_units[104] * fwd_Ea[104] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[11] -= q; /* CH3 */
    wdot[14] -= q; /* C2H5 */
    wdot[15] += q; /* CH4 */
    wdot[16] += q; /* C2H4 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[14];
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[399] -= dqdci;              /* dwdot[C2H5]/d[CH3] */
    J[400] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[401] += dqdci;              /* dwdot[C2H4]/d[CH3] */
    /* d()/d[C2H5] */
    dqdci =  + k_f*sc[11];
    J[501] -= dqdci;              /* dwdot[CH3]/d[C2H5] */
    J[504] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[505] += dqdci;              /* dwdot[CH4]/d[C2H5] */
    J[506] += dqdci;              /* dwdot[C2H4]/d[C2H5] */
    /* d()/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1204] -= dqdT;              /* dwdot[C2H5]/dT */
    J[1205] += dqdT;              /* dwdot[CH4]/dT */
    J[1206] += dqdT;              /* dwdot[C2H4]/dT */

    /*reaction 106: C2H5 + O2 => C2H4 + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[14];
    k_f = prefactor_units[105] * fwd_A[105]
                * exp(fwd_beta[105] * tc[0] - activation_units[105] * fwd_Ea[105] * invT);
    dlnkfdT = fwd_beta[105] * invT + activation_units[105] * fwd_Ea[105] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[8] += q; /* HO2 */
    wdot[14] -= q; /* C2H5 */
    wdot[16] += q; /* C2H4 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[14];
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[183] += dqdci;              /* dwdot[HO2]/d[O2] */
    J[189] -= dqdci;              /* dwdot[C2H5]/d[O2] */
    J[191] += dqdci;              /* dwdot[C2H4]/d[O2] */
    /* d()/d[C2H5] */
    dqdci =  + k_f*sc[5];
    J[495] -= dqdci;              /* dwdot[O2]/d[C2H5] */
    J[498] += dqdci;              /* dwdot[HO2]/d[C2H5] */
    J[504] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[506] += dqdci;              /* dwdot[C2H4]/d[C2H5] */
    /* d()/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */
    J[1204] -= dqdT;              /* dwdot[C2H5]/dT */
    J[1206] += dqdT;              /* dwdot[C2H4]/dT */

    /*reaction 107: C2H4 + HO2 => C2H5 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[16];
    k_f = prefactor_units[106] * fwd_A[106]
                * exp(fwd_beta[106] * tc[0] - activation_units[106] * fwd_Ea[106] * invT);
    dlnkfdT = fwd_beta[106] * invT + activation_units[106] * fwd_Ea[106] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[8] -= q; /* HO2 */
    wdot[14] += q; /* C2H5 */
    wdot[16] -= q; /* C2H4 */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[16];
    J[285] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[294] += dqdci;              /* dwdot[C2H5]/d[HO2] */
    J[296] -= dqdci;              /* dwdot[C2H4]/d[HO2] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[8];
    J[565] += dqdci;              /* dwdot[O2]/d[C2H4] */
    J[568] -= dqdci;              /* dwdot[HO2]/d[C2H4] */
    J[574] += dqdci;              /* dwdot[C2H5]/d[C2H4] */
    J[576] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    /* d()/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */
    J[1204] += dqdT;              /* dwdot[C2H5]/dT */
    J[1206] -= dqdT;              /* dwdot[C2H4]/dT */

    /*reaction 108: C2H5 + HO2 => C2H6 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[14];
    k_f = prefactor_units[107] * fwd_A[107]
                * exp(fwd_beta[107] * tc[0] - activation_units[107] * fwd_Ea[107] * invT);
    dlnkfdT = fwd_beta[107] * invT + activation_units[107] * fwd_Ea[107] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[8] -= q; /* HO2 */
    wdot[14] -= q; /* C2H5 */
    wdot[17] += q; /* C2H6 */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[14];
    J[285] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[294] -= dqdci;              /* dwdot[C2H5]/d[HO2] */
    J[297] += dqdci;              /* dwdot[C2H6]/d[HO2] */
    /* d()/d[C2H5] */
    dqdci =  + k_f*sc[8];
    J[495] += dqdci;              /* dwdot[O2]/d[C2H5] */
    J[498] -= dqdci;              /* dwdot[HO2]/d[C2H5] */
    J[504] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[507] += dqdci;              /* dwdot[C2H6]/d[C2H5] */
    /* d()/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */
    J[1204] -= dqdT;              /* dwdot[C2H5]/dT */
    J[1207] += dqdT;              /* dwdot[C2H6]/dT */

    /*reaction 109: H + C2H5 => C2H6 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[14];
    k_f = prefactor_units[108] * fwd_A[108]
                * exp(fwd_beta[108] * tc[0] - activation_units[108] * fwd_Ea[108] * invT);
    dlnkfdT = fwd_beta[108] * invT + activation_units[108] * fwd_Ea[108] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[14] -= q; /* C2H5 */
    wdot[17] += q; /* C2H6 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[14];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[84] -= dqdci;               /* dwdot[C2H5]/d[H] */
    J[87] += dqdci;               /* dwdot[C2H6]/d[H] */
    /* d()/d[C2H5] */
    dqdci =  + k_f*sc[2];
    J[492] -= dqdci;              /* dwdot[H]/d[C2H5] */
    J[504] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[507] += dqdci;              /* dwdot[C2H6]/d[C2H5] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1204] -= dqdT;              /* dwdot[C2H5]/dT */
    J[1207] += dqdT;              /* dwdot[C2H6]/dT */

    /*reaction 110: C2H6 + H => C2H5 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[17];
    k_f = prefactor_units[109] * fwd_A[109]
                * exp(fwd_beta[109] * tc[0] - activation_units[109] * fwd_Ea[109] * invT);
    dlnkfdT = fwd_beta[109] * invT + activation_units[109] * fwd_Ea[109] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[4] += q; /* H2 */
    wdot[14] += q; /* C2H5 */
    wdot[17] -= q; /* C2H6 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[17];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[74] += dqdci;               /* dwdot[H2]/d[H] */
    J[84] += dqdci;               /* dwdot[C2H5]/d[H] */
    J[87] -= dqdci;               /* dwdot[C2H6]/d[H] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[2];
    J[597] -= dqdci;              /* dwdot[H]/d[C2H6] */
    J[599] += dqdci;              /* dwdot[H2]/d[C2H6] */
    J[609] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[612] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */
    J[1204] += dqdT;              /* dwdot[C2H5]/dT */
    J[1207] -= dqdT;              /* dwdot[C2H6]/dT */

    /*reaction 111: C2H5 + H2 => C2H6 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[14];
    k_f = prefactor_units[110] * fwd_A[110]
                * exp(fwd_beta[110] * tc[0] - activation_units[110] * fwd_Ea[110] * invT);
    dlnkfdT = fwd_beta[110] * invT + activation_units[110] * fwd_Ea[110] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[4] -= q; /* H2 */
    wdot[14] -= q; /* C2H5 */
    wdot[17] += q; /* C2H6 */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[14];
    J[142] += dqdci;              /* dwdot[H]/d[H2] */
    J[144] -= dqdci;              /* dwdot[H2]/d[H2] */
    J[154] -= dqdci;              /* dwdot[C2H5]/d[H2] */
    J[157] += dqdci;              /* dwdot[C2H6]/d[H2] */
    /* d()/d[C2H5] */
    dqdci =  + k_f*sc[4];
    J[492] += dqdci;              /* dwdot[H]/d[C2H5] */
    J[494] -= dqdci;              /* dwdot[H2]/d[C2H5] */
    J[504] -= dqdci;              /* dwdot[C2H5]/d[C2H5] */
    J[507] += dqdci;              /* dwdot[C2H6]/d[C2H5] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1194] -= dqdT;              /* dwdot[H2]/dT */
    J[1204] -= dqdT;              /* dwdot[C2H5]/dT */
    J[1207] += dqdT;              /* dwdot[C2H6]/dT */

    /*reaction 112: C2H6 + OH => C2H5 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[17];
    k_f = prefactor_units[111] * fwd_A[111]
                * exp(fwd_beta[111] * tc[0] - activation_units[111] * fwd_Ea[111] * invT);
    dlnkfdT = fwd_beta[111] * invT + activation_units[111] * fwd_Ea[111] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[14] += q; /* C2H5 */
    wdot[17] -= q; /* C2H6 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[17];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[111] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[119] += dqdci;              /* dwdot[C2H5]/d[OH] */
    J[122] -= dqdci;              /* dwdot[C2H6]/d[OH] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[3];
    J[598] -= dqdci;              /* dwdot[OH]/d[C2H6] */
    J[601] += dqdci;              /* dwdot[H2O]/d[C2H6] */
    J[609] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[612] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */
    J[1204] += dqdT;              /* dwdot[C2H5]/dT */
    J[1207] -= dqdT;              /* dwdot[C2H6]/dT */

    /*reaction 113: CH2GSG + C2H6 => CH3 + C2H5 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[17];
    k_f = prefactor_units[112] * fwd_A[112]
                * exp(fwd_beta[112] * tc[0] - activation_units[112] * fwd_Ea[112] * invT);
    dlnkfdT = fwd_beta[112] * invT + activation_units[112] * fwd_Ea[112] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[9] -= q; /* CH2GSG */
    wdot[11] += q; /* CH3 */
    wdot[14] += q; /* C2H5 */
    wdot[17] -= q; /* C2H6 */
    /* d()/d[CH2GSG] */
    dqdci =  + k_f*sc[17];
    J[324] -= dqdci;              /* dwdot[CH2GSG]/d[CH2GSG] */
    J[326] += dqdci;              /* dwdot[CH3]/d[CH2GSG] */
    J[329] += dqdci;              /* dwdot[C2H5]/d[CH2GSG] */
    J[332] -= dqdci;              /* dwdot[C2H6]/d[CH2GSG] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[9];
    J[604] -= dqdci;              /* dwdot[CH2GSG]/d[C2H6] */
    J[606] += dqdci;              /* dwdot[CH3]/d[C2H6] */
    J[609] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[612] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[1199] -= dqdT;              /* dwdot[CH2GSG]/dT */
    J[1201] += dqdT;              /* dwdot[CH3]/dT */
    J[1204] += dqdT;              /* dwdot[C2H5]/dT */
    J[1207] -= dqdT;              /* dwdot[C2H6]/dT */

    /*reaction 114: C2H6 + O => C2H5 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[17];
    k_f = prefactor_units[113] * fwd_A[113]
                * exp(fwd_beta[113] * tc[0] - activation_units[113] * fwd_Ea[113] * invT);
    dlnkfdT = fwd_beta[113] * invT + activation_units[113] * fwd_Ea[113] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[3] += q; /* OH */
    wdot[14] += q; /* C2H5 */
    wdot[17] -= q; /* C2H6 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[17];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[38] += dqdci;               /* dwdot[OH]/d[O] */
    J[49] += dqdci;               /* dwdot[C2H5]/d[O] */
    J[52] -= dqdci;               /* dwdot[C2H6]/d[O] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[1];
    J[596] -= dqdci;              /* dwdot[O]/d[C2H6] */
    J[598] += dqdci;              /* dwdot[OH]/d[C2H6] */
    J[609] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[612] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1204] += dqdT;              /* dwdot[C2H5]/dT */
    J[1207] -= dqdT;              /* dwdot[C2H6]/dT */

    /*reaction 115: C2H6 + CH3 => C2H5 + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[17];
    k_f = prefactor_units[114] * fwd_A[114]
                * exp(fwd_beta[114] * tc[0] - activation_units[114] * fwd_Ea[114] * invT);
    dlnkfdT = fwd_beta[114] * invT + activation_units[114] * fwd_Ea[114] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[11] -= q; /* CH3 */
    wdot[14] += q; /* C2H5 */
    wdot[15] += q; /* CH4 */
    wdot[17] -= q; /* C2H6 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[17];
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[399] += dqdci;              /* dwdot[C2H5]/d[CH3] */
    J[400] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[402] -= dqdci;              /* dwdot[C2H6]/d[CH3] */
    /* d()/d[C2H6] */
    dqdci =  + k_f*sc[11];
    J[606] -= dqdci;              /* dwdot[CH3]/d[C2H6] */
    J[609] += dqdci;              /* dwdot[C2H5]/d[C2H6] */
    J[610] += dqdci;              /* dwdot[CH4]/d[C2H6] */
    J[612] -= dqdci;              /* dwdot[C2H6]/d[C2H6] */
    /* d()/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1204] += dqdT;              /* dwdot[C2H5]/dT */
    J[1205] += dqdT;              /* dwdot[CH4]/dT */
    J[1207] -= dqdT;              /* dwdot[C2H6]/dT */

    /*reaction 116: HCCO + O => H + 2.000000 CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[23];
    k_f = prefactor_units[115] * fwd_A[115]
                * exp(fwd_beta[115] * tc[0] - activation_units[115] * fwd_Ea[115] * invT);
    dlnkfdT = fwd_beta[115] * invT + activation_units[115] * fwd_Ea[115] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[2] += q; /* H */
    wdot[12] += 2 * q; /* CO */
    wdot[23] -= q; /* HCCO */
    /* d()/d[O] */
    dqdci =  + k_f*sc[23];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[37] += dqdci;               /* dwdot[H]/d[O] */
    J[47] += 2 * dqdci;           /* dwdot[CO]/d[O] */
    J[58] -= dqdci;               /* dwdot[HCCO]/d[O] */
    /* d()/d[HCCO] */
    dqdci =  + k_f*sc[1];
    J[806] -= dqdci;              /* dwdot[O]/d[HCCO] */
    J[807] += dqdci;              /* dwdot[H]/d[HCCO] */
    J[817] += 2 * dqdci;          /* dwdot[CO]/d[HCCO] */
    J[828] -= dqdci;              /* dwdot[HCCO]/d[HCCO] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1202] += 2 * dqdT;          /* dwdot[CO]/dT */
    J[1213] -= dqdT;              /* dwdot[HCCO]/dT */

    /*reaction 117: HCCO + O2 => CO2 + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[23];
    k_f = prefactor_units[116] * fwd_A[116]
                * exp(fwd_beta[116] * tc[0] - activation_units[116] * fwd_Ea[116] * invT);
    dlnkfdT = fwd_beta[116] * invT + activation_units[116] * fwd_Ea[116] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[18] += q; /* CO2 */
    wdot[19] += q; /* HCO */
    wdot[23] -= q; /* HCCO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[23];
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[193] += dqdci;              /* dwdot[CO2]/d[O2] */
    J[194] += dqdci;              /* dwdot[HCO]/d[O2] */
    J[198] -= dqdci;              /* dwdot[HCCO]/d[O2] */
    /* d()/d[HCCO] */
    dqdci =  + k_f*sc[5];
    J[810] -= dqdci;              /* dwdot[O2]/d[HCCO] */
    J[823] += dqdci;              /* dwdot[CO2]/d[HCCO] */
    J[824] += dqdci;              /* dwdot[HCO]/d[HCCO] */
    J[828] -= dqdci;              /* dwdot[HCCO]/d[HCCO] */
    /* d()/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1208] += dqdT;              /* dwdot[CO2]/dT */
    J[1209] += dqdT;              /* dwdot[HCO]/dT */
    J[1213] -= dqdT;              /* dwdot[HCCO]/dT */

    /*reaction 118: HCCO + OH => 2.000000 HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[23];
    k_f = prefactor_units[117] * fwd_A[117]
                * exp(fwd_beta[117] * tc[0] - activation_units[117] * fwd_Ea[117] * invT);
    dlnkfdT = fwd_beta[117] * invT + activation_units[117] * fwd_Ea[117] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[19] += 2 * q; /* HCO */
    wdot[23] -= q; /* HCCO */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[23];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[124] += 2 * dqdci;          /* dwdot[HCO]/d[OH] */
    J[128] -= dqdci;              /* dwdot[HCCO]/d[OH] */
    /* d()/d[HCCO] */
    dqdci =  + k_f*sc[3];
    J[808] -= dqdci;              /* dwdot[OH]/d[HCCO] */
    J[824] += 2 * dqdci;          /* dwdot[HCO]/d[HCCO] */
    J[828] -= dqdci;              /* dwdot[HCCO]/d[HCCO] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1209] += 2 * dqdT;          /* dwdot[HCO]/dT */
    J[1213] -= dqdT;              /* dwdot[HCCO]/dT */

    /*reaction 119: HCCO + H => CH2GSG + CO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[23];
    k_f = prefactor_units[118] * fwd_A[118]
                * exp(fwd_beta[118] * tc[0] - activation_units[118] * fwd_Ea[118] * invT);
    dlnkfdT = fwd_beta[118] * invT + activation_units[118] * fwd_Ea[118] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[9] += q; /* CH2GSG */
    wdot[12] += q; /* CO */
    wdot[23] -= q; /* HCCO */
    /* d()/d[H] */
    dqdci =  + k_f*sc[23];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[79] += dqdci;               /* dwdot[CH2GSG]/d[H] */
    J[82] += dqdci;               /* dwdot[CO]/d[H] */
    J[93] -= dqdci;               /* dwdot[HCCO]/d[H] */
    /* d()/d[HCCO] */
    dqdci =  + k_f*sc[2];
    J[807] -= dqdci;              /* dwdot[H]/d[HCCO] */
    J[814] += dqdci;              /* dwdot[CH2GSG]/d[HCCO] */
    J[817] += dqdci;              /* dwdot[CO]/d[HCCO] */
    J[828] -= dqdci;              /* dwdot[HCCO]/d[HCCO] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1199] += dqdT;              /* dwdot[CH2GSG]/dT */
    J[1202] += dqdT;              /* dwdot[CO]/dT */
    J[1213] -= dqdT;              /* dwdot[HCCO]/dT */

    /*reaction 120: CH2GSG + CO => HCCO + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[9]*sc[12];
    k_f = prefactor_units[119] * fwd_A[119]
                * exp(fwd_beta[119] * tc[0] - activation_units[119] * fwd_Ea[119] * invT);
    dlnkfdT = fwd_beta[119] * invT + activation_units[119] * fwd_Ea[119] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[9] -= q; /* CH2GSG */
    wdot[12] -= q; /* CO */
    wdot[23] += q; /* HCCO */
    /* d()/d[CH2GSG] */
    dqdci =  + k_f*sc[12];
    J[317] += dqdci;              /* dwdot[H]/d[CH2GSG] */
    J[324] -= dqdci;              /* dwdot[CH2GSG]/d[CH2GSG] */
    J[327] -= dqdci;              /* dwdot[CO]/d[CH2GSG] */
    J[338] += dqdci;              /* dwdot[HCCO]/d[CH2GSG] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[9];
    J[422] += dqdci;              /* dwdot[H]/d[CO] */
    J[429] -= dqdci;              /* dwdot[CH2GSG]/d[CO] */
    J[432] -= dqdci;              /* dwdot[CO]/d[CO] */
    J[443] += dqdci;              /* dwdot[HCCO]/d[CO] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1199] -= dqdT;              /* dwdot[CH2GSG]/dT */
    J[1202] -= dqdT;              /* dwdot[CO]/dT */
    J[1213] += dqdT;              /* dwdot[HCCO]/dT */

    /*reaction 121: CH2CHO + O2 => CH2O + CO + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[25];
    k_f = prefactor_units[120] * fwd_A[120]
                * exp(fwd_beta[120] * tc[0] - activation_units[120] * fwd_Ea[120] * invT);
    dlnkfdT = fwd_beta[120] * invT + activation_units[120] * fwd_Ea[120] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] += q; /* OH */
    wdot[5] -= q; /* O2 */
    wdot[10] += q; /* CH2O */
    wdot[12] += q; /* CO */
    wdot[25] -= q; /* CH2CHO */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[25];
    J[178] += dqdci;              /* dwdot[OH]/d[O2] */
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[185] += dqdci;              /* dwdot[CH2O]/d[O2] */
    J[187] += dqdci;              /* dwdot[CO]/d[O2] */
    J[200] -= dqdci;              /* dwdot[CH2CHO]/d[O2] */
    /* d()/d[CH2CHO] */
    dqdci =  + k_f*sc[5];
    J[878] += dqdci;              /* dwdot[OH]/d[CH2CHO] */
    J[880] -= dqdci;              /* dwdot[O2]/d[CH2CHO] */
    J[885] += dqdci;              /* dwdot[CH2O]/d[CH2CHO] */
    J[887] += dqdci;              /* dwdot[CO]/d[CH2CHO] */
    J[900] -= dqdci;              /* dwdot[CH2CHO]/d[CH2CHO] */
    /* d()/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1200] += dqdT;              /* dwdot[CH2O]/dT */
    J[1202] += dqdT;              /* dwdot[CO]/dT */
    J[1215] -= dqdT;              /* dwdot[CH2CHO]/dT */

    /*reaction 122: C3H5XA + CH2O => C3H6 + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[10]*sc[26];
    k_f = prefactor_units[121] * fwd_A[121]
                * exp(fwd_beta[121] * tc[0] - activation_units[121] * fwd_Ea[121] * invT);
    dlnkfdT = fwd_beta[121] * invT + activation_units[121] * fwd_Ea[121] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[10] -= q; /* CH2O */
    wdot[19] += q; /* HCO */
    wdot[26] -= q; /* C3H5XA */
    wdot[27] += q; /* C3H6 */
    /* d()/d[CH2O] */
    dqdci =  + k_f*sc[26];
    J[360] -= dqdci;              /* dwdot[CH2O]/d[CH2O] */
    J[369] += dqdci;              /* dwdot[HCO]/d[CH2O] */
    J[376] -= dqdci;              /* dwdot[C3H5XA]/d[CH2O] */
    J[377] += dqdci;              /* dwdot[C3H6]/d[CH2O] */
    /* d()/d[C3H5XA] */
    dqdci =  + k_f*sc[10];
    J[920] -= dqdci;              /* dwdot[CH2O]/d[C3H5XA] */
    J[929] += dqdci;              /* dwdot[HCO]/d[C3H5XA] */
    J[936] -= dqdci;              /* dwdot[C3H5XA]/d[C3H5XA] */
    J[937] += dqdci;              /* dwdot[C3H6]/d[C3H5XA] */
    /* d()/dT */
    J[1200] -= dqdT;              /* dwdot[CH2O]/dT */
    J[1209] += dqdT;              /* dwdot[HCO]/dT */
    J[1216] -= dqdT;              /* dwdot[C3H5XA]/dT */
    J[1217] += dqdT;              /* dwdot[C3H6]/dT */

    /*reaction 123: C3H6 + HCO => C3H5XA + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[19]*sc[27];
    k_f = prefactor_units[122] * fwd_A[122]
                * exp(fwd_beta[122] * tc[0] - activation_units[122] * fwd_Ea[122] * invT);
    dlnkfdT = fwd_beta[122] * invT + activation_units[122] * fwd_Ea[122] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[10] += q; /* CH2O */
    wdot[19] -= q; /* HCO */
    wdot[26] += q; /* C3H5XA */
    wdot[27] -= q; /* C3H6 */
    /* d()/d[HCO] */
    dqdci =  + k_f*sc[27];
    J[675] += dqdci;              /* dwdot[CH2O]/d[HCO] */
    J[684] -= dqdci;              /* dwdot[HCO]/d[HCO] */
    J[691] += dqdci;              /* dwdot[C3H5XA]/d[HCO] */
    J[692] -= dqdci;              /* dwdot[C3H6]/d[HCO] */
    /* d()/d[C3H6] */
    dqdci =  + k_f*sc[19];
    J[955] += dqdci;              /* dwdot[CH2O]/d[C3H6] */
    J[964] -= dqdci;              /* dwdot[HCO]/d[C3H6] */
    J[971] += dqdci;              /* dwdot[C3H5XA]/d[C3H6] */
    J[972] -= dqdci;              /* dwdot[C3H6]/d[C3H6] */
    /* d()/dT */
    J[1200] += dqdT;              /* dwdot[CH2O]/dT */
    J[1209] -= dqdT;              /* dwdot[HCO]/dT */
    J[1216] += dqdT;              /* dwdot[C3H5XA]/dT */
    J[1217] -= dqdT;              /* dwdot[C3H6]/dT */

    /*reaction 124: C3H5XA + HO2 => C3H5O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[26];
    k_f = prefactor_units[123] * fwd_A[123]
                * exp(fwd_beta[123] * tc[0] - activation_units[123] * fwd_Ea[123] * invT);
    dlnkfdT = fwd_beta[123] * invT + activation_units[123] * fwd_Ea[123] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] += q; /* OH */
    wdot[8] -= q; /* HO2 */
    wdot[26] -= q; /* C3H5XA */
    wdot[28] += q; /* C3H5O */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[26];
    J[283] += dqdci;              /* dwdot[OH]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[306] -= dqdci;              /* dwdot[C3H5XA]/d[HO2] */
    J[308] += dqdci;              /* dwdot[C3H5O]/d[HO2] */
    /* d()/d[C3H5XA] */
    dqdci =  + k_f*sc[8];
    J[913] += dqdci;              /* dwdot[OH]/d[C3H5XA] */
    J[918] -= dqdci;              /* dwdot[HO2]/d[C3H5XA] */
    J[936] -= dqdci;              /* dwdot[C3H5XA]/d[C3H5XA] */
    J[938] += dqdci;              /* dwdot[C3H5O]/d[C3H5XA] */
    /* d()/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */
    J[1216] -= dqdT;              /* dwdot[C3H5XA]/dT */
    J[1218] += dqdT;              /* dwdot[C3H5O]/dT */

    /*reaction 125: C3H5XA + CH3O2 => C3H5O + CH3O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[20]*sc[26];
    k_f = prefactor_units[124] * fwd_A[124]
                * exp(fwd_beta[124] * tc[0] - activation_units[124] * fwd_Ea[124] * invT);
    dlnkfdT = fwd_beta[124] * invT + activation_units[124] * fwd_Ea[124] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[13] += q; /* CH3O */
    wdot[20] -= q; /* CH3O2 */
    wdot[26] -= q; /* C3H5XA */
    wdot[28] += q; /* C3H5O */
    /* d()/d[CH3O2] */
    dqdci =  + k_f*sc[26];
    J[713] += dqdci;              /* dwdot[CH3O]/d[CH3O2] */
    J[720] -= dqdci;              /* dwdot[CH3O2]/d[CH3O2] */
    J[726] -= dqdci;              /* dwdot[C3H5XA]/d[CH3O2] */
    J[728] += dqdci;              /* dwdot[C3H5O]/d[CH3O2] */
    /* d()/d[C3H5XA] */
    dqdci =  + k_f*sc[20];
    J[923] += dqdci;              /* dwdot[CH3O]/d[C3H5XA] */
    J[930] -= dqdci;              /* dwdot[CH3O2]/d[C3H5XA] */
    J[936] -= dqdci;              /* dwdot[C3H5XA]/d[C3H5XA] */
    J[938] += dqdci;              /* dwdot[C3H5O]/d[C3H5XA] */
    /* d()/dT */
    J[1203] += dqdT;              /* dwdot[CH3O]/dT */
    J[1210] -= dqdT;              /* dwdot[CH3O2]/dT */
    J[1216] -= dqdT;              /* dwdot[C3H5XA]/dT */
    J[1218] += dqdT;              /* dwdot[C3H5O]/dT */

    /*reaction 126: C3H5XA + HO2 => C3H6 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[26];
    k_f = prefactor_units[125] * fwd_A[125]
                * exp(fwd_beta[125] * tc[0] - activation_units[125] * fwd_Ea[125] * invT);
    dlnkfdT = fwd_beta[125] * invT + activation_units[125] * fwd_Ea[125] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[8] -= q; /* HO2 */
    wdot[26] -= q; /* C3H5XA */
    wdot[27] += q; /* C3H6 */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[26];
    J[285] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[306] -= dqdci;              /* dwdot[C3H5XA]/d[HO2] */
    J[307] += dqdci;              /* dwdot[C3H6]/d[HO2] */
    /* d()/d[C3H5XA] */
    dqdci =  + k_f*sc[8];
    J[915] += dqdci;              /* dwdot[O2]/d[C3H5XA] */
    J[918] -= dqdci;              /* dwdot[HO2]/d[C3H5XA] */
    J[936] -= dqdci;              /* dwdot[C3H5XA]/d[C3H5XA] */
    J[937] += dqdci;              /* dwdot[C3H6]/d[C3H5XA] */
    /* d()/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */
    J[1216] -= dqdT;              /* dwdot[C3H5XA]/dT */
    J[1217] += dqdT;              /* dwdot[C3H6]/dT */

    /*reaction 127: C3H5XA + O2 => CH2CHO + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[26];
    k_f = prefactor_units[126] * fwd_A[126]
                * exp(fwd_beta[126] * tc[0] - activation_units[126] * fwd_Ea[126] * invT);
    dlnkfdT = fwd_beta[126] * invT + activation_units[126] * fwd_Ea[126] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[10] += q; /* CH2O */
    wdot[25] += q; /* CH2CHO */
    wdot[26] -= q; /* C3H5XA */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[26];
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[185] += dqdci;              /* dwdot[CH2O]/d[O2] */
    J[200] += dqdci;              /* dwdot[CH2CHO]/d[O2] */
    J[201] -= dqdci;              /* dwdot[C3H5XA]/d[O2] */
    /* d()/d[C3H5XA] */
    dqdci =  + k_f*sc[5];
    J[915] -= dqdci;              /* dwdot[O2]/d[C3H5XA] */
    J[920] += dqdci;              /* dwdot[CH2O]/d[C3H5XA] */
    J[935] += dqdci;              /* dwdot[CH2CHO]/d[C3H5XA] */
    J[936] -= dqdci;              /* dwdot[C3H5XA]/d[C3H5XA] */
    /* d()/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1200] += dqdT;              /* dwdot[CH2O]/dT */
    J[1215] += dqdT;              /* dwdot[CH2CHO]/dT */
    J[1216] -= dqdT;              /* dwdot[C3H5XA]/dT */

    /*reaction 128: C3H5XA + O2 => C2H2 + CH2O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[26];
    k_f = prefactor_units[127] * fwd_A[127]
                * exp(fwd_beta[127] * tc[0] - activation_units[127] * fwd_Ea[127] * invT);
    dlnkfdT = fwd_beta[127] * invT + activation_units[127] * fwd_Ea[127] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] += q; /* OH */
    wdot[5] -= q; /* O2 */
    wdot[10] += q; /* CH2O */
    wdot[22] += q; /* C2H2 */
    wdot[26] -= q; /* C3H5XA */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[26];
    J[178] += dqdci;              /* dwdot[OH]/d[O2] */
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[185] += dqdci;              /* dwdot[CH2O]/d[O2] */
    J[197] += dqdci;              /* dwdot[C2H2]/d[O2] */
    J[201] -= dqdci;              /* dwdot[C3H5XA]/d[O2] */
    /* d()/d[C3H5XA] */
    dqdci =  + k_f*sc[5];
    J[913] += dqdci;              /* dwdot[OH]/d[C3H5XA] */
    J[915] -= dqdci;              /* dwdot[O2]/d[C3H5XA] */
    J[920] += dqdci;              /* dwdot[CH2O]/d[C3H5XA] */
    J[932] += dqdci;              /* dwdot[C2H2]/d[C3H5XA] */
    J[936] -= dqdci;              /* dwdot[C3H5XA]/d[C3H5XA] */
    /* d()/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1200] += dqdT;              /* dwdot[CH2O]/dT */
    J[1212] += dqdT;              /* dwdot[C2H2]/dT */
    J[1216] -= dqdT;              /* dwdot[C3H5XA]/dT */

    /*reaction 129: C3H5XA + H => C3H6 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[26];
    k_f = prefactor_units[128] * fwd_A[128]
                * exp(fwd_beta[128] * tc[0] - activation_units[128] * fwd_Ea[128] * invT);
    dlnkfdT = fwd_beta[128] * invT + activation_units[128] * fwd_Ea[128] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[26] -= q; /* C3H5XA */
    wdot[27] += q; /* C3H6 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[26];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[96] -= dqdci;               /* dwdot[C3H5XA]/d[H] */
    J[97] += dqdci;               /* dwdot[C3H6]/d[H] */
    /* d()/d[C3H5XA] */
    dqdci =  + k_f*sc[2];
    J[912] -= dqdci;              /* dwdot[H]/d[C3H5XA] */
    J[936] -= dqdci;              /* dwdot[C3H5XA]/d[C3H5XA] */
    J[937] += dqdci;              /* dwdot[C3H6]/d[C3H5XA] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1216] -= dqdT;              /* dwdot[C3H5XA]/dT */
    J[1217] += dqdT;              /* dwdot[C3H6]/dT */

    /*reaction 130: C3H5XA => C2H2 + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[26];
    k_f = prefactor_units[129] * fwd_A[129]
                * exp(fwd_beta[129] * tc[0] - activation_units[129] * fwd_Ea[129] * invT);
    dlnkfdT = fwd_beta[129] * invT + activation_units[129] * fwd_Ea[129] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[11] += q; /* CH3 */
    wdot[22] += q; /* C2H2 */
    wdot[26] -= q; /* C3H5XA */
    /* d()/d[C3H5XA] */
    dqdci =  + k_f;
    J[921] += dqdci;              /* dwdot[CH3]/d[C3H5XA] */
    J[932] += dqdci;              /* dwdot[C2H2]/d[C3H5XA] */
    J[936] -= dqdci;              /* dwdot[C3H5XA]/d[C3H5XA] */
    /* d()/dT */
    J[1201] += dqdT;              /* dwdot[CH3]/dT */
    J[1212] += dqdT;              /* dwdot[C2H2]/dT */
    J[1216] -= dqdT;              /* dwdot[C3H5XA]/dT */

    /*reaction 131: C3H6 => C2H3 + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[27];
    k_f = prefactor_units[130] * fwd_A[130]
                * exp(fwd_beta[130] * tc[0] - activation_units[130] * fwd_Ea[130] * invT);
    dlnkfdT = fwd_beta[130] * invT + activation_units[130] * fwd_Ea[130] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[11] += q; /* CH3 */
    wdot[24] += q; /* C2H3 */
    wdot[27] -= q; /* C3H6 */
    /* d()/d[C3H6] */
    dqdci =  + k_f;
    J[956] += dqdci;              /* dwdot[CH3]/d[C3H6] */
    J[969] += dqdci;              /* dwdot[C2H3]/d[C3H6] */
    J[972] -= dqdci;              /* dwdot[C3H6]/d[C3H6] */
    /* d()/dT */
    J[1201] += dqdT;              /* dwdot[CH3]/dT */
    J[1214] += dqdT;              /* dwdot[C2H3]/dT */
    J[1217] -= dqdT;              /* dwdot[C3H6]/dT */

    /*reaction 132: C2H3 + CH3 => C3H6 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[24];
    k_f = prefactor_units[131] * fwd_A[131]
                * exp(fwd_beta[131] * tc[0] - activation_units[131] * fwd_Ea[131] * invT);
    dlnkfdT = fwd_beta[131] * invT + activation_units[131] * fwd_Ea[131] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[11] -= q; /* CH3 */
    wdot[24] -= q; /* C2H3 */
    wdot[27] += q; /* C3H6 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[24];
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[409] -= dqdci;              /* dwdot[C2H3]/d[CH3] */
    J[412] += dqdci;              /* dwdot[C3H6]/d[CH3] */
    /* d()/d[C2H3] */
    dqdci =  + k_f*sc[11];
    J[851] -= dqdci;              /* dwdot[CH3]/d[C2H3] */
    J[864] -= dqdci;              /* dwdot[C2H3]/d[C2H3] */
    J[867] += dqdci;              /* dwdot[C3H6]/d[C2H3] */
    /* d()/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1214] -= dqdT;              /* dwdot[C2H3]/dT */
    J[1217] += dqdT;              /* dwdot[C3H6]/dT */

    /*reaction 133: C3H6 + O => C2H5 + HCO */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[27];
    k_f = prefactor_units[132] * fwd_A[132]
                * exp(fwd_beta[132] * tc[0] - activation_units[132] * fwd_Ea[132] * invT);
    dlnkfdT = fwd_beta[132] * invT + activation_units[132] * fwd_Ea[132] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[14] += q; /* C2H5 */
    wdot[19] += q; /* HCO */
    wdot[27] -= q; /* C3H6 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[27];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[49] += dqdci;               /* dwdot[C2H5]/d[O] */
    J[54] += dqdci;               /* dwdot[HCO]/d[O] */
    J[62] -= dqdci;               /* dwdot[C3H6]/d[O] */
    /* d()/d[C3H6] */
    dqdci =  + k_f*sc[1];
    J[946] -= dqdci;              /* dwdot[O]/d[C3H6] */
    J[959] += dqdci;              /* dwdot[C2H5]/d[C3H6] */
    J[964] += dqdci;              /* dwdot[HCO]/d[C3H6] */
    J[972] -= dqdci;              /* dwdot[C3H6]/d[C3H6] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1204] += dqdT;              /* dwdot[C2H5]/dT */
    J[1209] += dqdT;              /* dwdot[HCO]/dT */
    J[1217] -= dqdT;              /* dwdot[C3H6]/dT */

    /*reaction 134: C3H6 + H => C2H4 + CH3 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[27];
    k_f = prefactor_units[133] * fwd_A[133]
                * exp(fwd_beta[133] * tc[0] - activation_units[133] * fwd_Ea[133] * invT);
    dlnkfdT = fwd_beta[133] * invT + activation_units[133] * fwd_Ea[133] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[11] += q; /* CH3 */
    wdot[16] += q; /* C2H4 */
    wdot[27] -= q; /* C3H6 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[27];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[81] += dqdci;               /* dwdot[CH3]/d[H] */
    J[86] += dqdci;               /* dwdot[C2H4]/d[H] */
    J[97] -= dqdci;               /* dwdot[C3H6]/d[H] */
    /* d()/d[C3H6] */
    dqdci =  + k_f*sc[2];
    J[947] -= dqdci;              /* dwdot[H]/d[C3H6] */
    J[956] += dqdci;              /* dwdot[CH3]/d[C3H6] */
    J[961] += dqdci;              /* dwdot[C2H4]/d[C3H6] */
    J[972] -= dqdci;              /* dwdot[C3H6]/d[C3H6] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1201] += dqdT;              /* dwdot[CH3]/dT */
    J[1206] += dqdT;              /* dwdot[C2H4]/dT */
    J[1217] -= dqdT;              /* dwdot[C3H6]/dT */

    /*reaction 135: C2H4 + CH3 => C3H6 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[16];
    k_f = prefactor_units[134] * fwd_A[134]
                * exp(fwd_beta[134] * tc[0] - activation_units[134] * fwd_Ea[134] * invT);
    dlnkfdT = fwd_beta[134] * invT + activation_units[134] * fwd_Ea[134] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[11] -= q; /* CH3 */
    wdot[16] -= q; /* C2H4 */
    wdot[27] += q; /* C3H6 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[16];
    J[387] += dqdci;              /* dwdot[H]/d[CH3] */
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[401] -= dqdci;              /* dwdot[C2H4]/d[CH3] */
    J[412] += dqdci;              /* dwdot[C3H6]/d[CH3] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[11];
    J[562] += dqdci;              /* dwdot[H]/d[C2H4] */
    J[571] -= dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[576] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[587] += dqdci;              /* dwdot[C3H6]/d[C2H4] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1206] -= dqdT;              /* dwdot[C2H4]/dT */
    J[1217] += dqdT;              /* dwdot[C3H6]/dT */

    /*reaction 136: C3H6 + H => C3H5XA + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[27];
    k_f = prefactor_units[135] * fwd_A[135]
                * exp(fwd_beta[135] * tc[0] - activation_units[135] * fwd_Ea[135] * invT);
    dlnkfdT = fwd_beta[135] * invT + activation_units[135] * fwd_Ea[135] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[4] += q; /* H2 */
    wdot[26] += q; /* C3H5XA */
    wdot[27] -= q; /* C3H6 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[27];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[74] += dqdci;               /* dwdot[H2]/d[H] */
    J[96] += dqdci;               /* dwdot[C3H5XA]/d[H] */
    J[97] -= dqdci;               /* dwdot[C3H6]/d[H] */
    /* d()/d[C3H6] */
    dqdci =  + k_f*sc[2];
    J[947] -= dqdci;              /* dwdot[H]/d[C3H6] */
    J[949] += dqdci;              /* dwdot[H2]/d[C3H6] */
    J[971] += dqdci;              /* dwdot[C3H5XA]/d[C3H6] */
    J[972] -= dqdci;              /* dwdot[C3H6]/d[C3H6] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */
    J[1216] += dqdT;              /* dwdot[C3H5XA]/dT */
    J[1217] -= dqdT;              /* dwdot[C3H6]/dT */

    /*reaction 137: C3H5XA + H2 => C3H6 + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[26];
    k_f = prefactor_units[136] * fwd_A[136]
                * exp(fwd_beta[136] * tc[0] - activation_units[136] * fwd_Ea[136] * invT);
    dlnkfdT = fwd_beta[136] * invT + activation_units[136] * fwd_Ea[136] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[4] -= q; /* H2 */
    wdot[26] -= q; /* C3H5XA */
    wdot[27] += q; /* C3H6 */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[26];
    J[142] += dqdci;              /* dwdot[H]/d[H2] */
    J[144] -= dqdci;              /* dwdot[H2]/d[H2] */
    J[166] -= dqdci;              /* dwdot[C3H5XA]/d[H2] */
    J[167] += dqdci;              /* dwdot[C3H6]/d[H2] */
    /* d()/d[C3H5XA] */
    dqdci =  + k_f*sc[4];
    J[912] += dqdci;              /* dwdot[H]/d[C3H5XA] */
    J[914] -= dqdci;              /* dwdot[H2]/d[C3H5XA] */
    J[936] -= dqdci;              /* dwdot[C3H5XA]/d[C3H5XA] */
    J[937] += dqdci;              /* dwdot[C3H6]/d[C3H5XA] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1194] -= dqdT;              /* dwdot[H2]/dT */
    J[1216] -= dqdT;              /* dwdot[C3H5XA]/dT */
    J[1217] += dqdT;              /* dwdot[C3H6]/dT */

    /*reaction 138: C3H6 + OH => C3H5XA + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[27];
    k_f = prefactor_units[137] * fwd_A[137]
                * exp(fwd_beta[137] * tc[0] - activation_units[137] * fwd_Ea[137] * invT);
    dlnkfdT = fwd_beta[137] * invT + activation_units[137] * fwd_Ea[137] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[26] += q; /* C3H5XA */
    wdot[27] -= q; /* C3H6 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[27];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[111] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[131] += dqdci;              /* dwdot[C3H5XA]/d[OH] */
    J[132] -= dqdci;              /* dwdot[C3H6]/d[OH] */
    /* d()/d[C3H6] */
    dqdci =  + k_f*sc[3];
    J[948] -= dqdci;              /* dwdot[OH]/d[C3H6] */
    J[951] += dqdci;              /* dwdot[H2O]/d[C3H6] */
    J[971] += dqdci;              /* dwdot[C3H5XA]/d[C3H6] */
    J[972] -= dqdci;              /* dwdot[C3H6]/d[C3H6] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */
    J[1216] += dqdT;              /* dwdot[C3H5XA]/dT */
    J[1217] -= dqdT;              /* dwdot[C3H6]/dT */

    /*reaction 139: C3H6 + O => C3H5XA + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[27];
    k_f = prefactor_units[138] * fwd_A[138]
                * exp(fwd_beta[138] * tc[0] - activation_units[138] * fwd_Ea[138] * invT);
    dlnkfdT = fwd_beta[138] * invT + activation_units[138] * fwd_Ea[138] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[3] += q; /* OH */
    wdot[26] += q; /* C3H5XA */
    wdot[27] -= q; /* C3H6 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[27];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[38] += dqdci;               /* dwdot[OH]/d[O] */
    J[61] += dqdci;               /* dwdot[C3H5XA]/d[O] */
    J[62] -= dqdci;               /* dwdot[C3H6]/d[O] */
    /* d()/d[C3H6] */
    dqdci =  + k_f*sc[1];
    J[946] -= dqdci;              /* dwdot[O]/d[C3H6] */
    J[948] += dqdci;              /* dwdot[OH]/d[C3H6] */
    J[971] += dqdci;              /* dwdot[C3H5XA]/d[C3H6] */
    J[972] -= dqdci;              /* dwdot[C3H6]/d[C3H6] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1216] += dqdT;              /* dwdot[C3H5XA]/dT */
    J[1217] -= dqdT;              /* dwdot[C3H6]/dT */

    /*reaction 140: C3H6 + CH3 => C3H5XA + CH4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[27];
    k_f = prefactor_units[139] * fwd_A[139]
                * exp(fwd_beta[139] * tc[0] - activation_units[139] * fwd_Ea[139] * invT);
    dlnkfdT = fwd_beta[139] * invT + activation_units[139] * fwd_Ea[139] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[11] -= q; /* CH3 */
    wdot[15] += q; /* CH4 */
    wdot[26] += q; /* C3H5XA */
    wdot[27] -= q; /* C3H6 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[27];
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[400] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[411] += dqdci;              /* dwdot[C3H5XA]/d[CH3] */
    J[412] -= dqdci;              /* dwdot[C3H6]/d[CH3] */
    /* d()/d[C3H6] */
    dqdci =  + k_f*sc[11];
    J[956] -= dqdci;              /* dwdot[CH3]/d[C3H6] */
    J[960] += dqdci;              /* dwdot[CH4]/d[C3H6] */
    J[971] += dqdci;              /* dwdot[C3H5XA]/d[C3H6] */
    J[972] -= dqdci;              /* dwdot[C3H6]/d[C3H6] */
    /* d()/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1205] += dqdT;              /* dwdot[CH4]/dT */
    J[1216] += dqdT;              /* dwdot[C3H5XA]/dT */
    J[1217] -= dqdT;              /* dwdot[C3H6]/dT */

    /*reaction 141: IXC3H7 + O2 => C3H6 + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[29];
    k_f = prefactor_units[140] * fwd_A[140]
                * exp(fwd_beta[140] * tc[0] - activation_units[140] * fwd_Ea[140] * invT);
    dlnkfdT = fwd_beta[140] * invT + activation_units[140] * fwd_Ea[140] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[8] += q; /* HO2 */
    wdot[27] += q; /* C3H6 */
    wdot[29] -= q; /* IXC3H7 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[29];
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[183] += dqdci;              /* dwdot[HO2]/d[O2] */
    J[202] += dqdci;              /* dwdot[C3H6]/d[O2] */
    J[204] -= dqdci;              /* dwdot[IXC3H7]/d[O2] */
    /* d()/d[IXC3H7] */
    dqdci =  + k_f*sc[5];
    J[1020] -= dqdci;             /* dwdot[O2]/d[IXC3H7] */
    J[1023] += dqdci;             /* dwdot[HO2]/d[IXC3H7] */
    J[1042] += dqdci;             /* dwdot[C3H6]/d[IXC3H7] */
    J[1044] -= dqdci;             /* dwdot[IXC3H7]/d[IXC3H7] */
    /* d()/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */
    J[1217] += dqdT;              /* dwdot[C3H6]/dT */
    J[1219] -= dqdT;              /* dwdot[IXC3H7]/dT */

    /*reaction 142: IXC3H7 => H + C3H6 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[29];
    k_f = prefactor_units[141] * fwd_A[141]
                * exp(fwd_beta[141] * tc[0] - activation_units[141] * fwd_Ea[141] * invT);
    dlnkfdT = fwd_beta[141] * invT + activation_units[141] * fwd_Ea[141] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[27] += q; /* C3H6 */
    wdot[29] -= q; /* IXC3H7 */
    /* d()/d[IXC3H7] */
    dqdci =  + k_f;
    J[1017] += dqdci;             /* dwdot[H]/d[IXC3H7] */
    J[1042] += dqdci;             /* dwdot[C3H6]/d[IXC3H7] */
    J[1044] -= dqdci;             /* dwdot[IXC3H7]/d[IXC3H7] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1217] += dqdT;              /* dwdot[C3H6]/dT */
    J[1219] -= dqdT;              /* dwdot[IXC3H7]/dT */

    /*reaction 143: H + C3H6 => IXC3H7 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[27];
    k_f = prefactor_units[142] * fwd_A[142]
                * exp(fwd_beta[142] * tc[0] - activation_units[142] * fwd_Ea[142] * invT);
    dlnkfdT = fwd_beta[142] * invT + activation_units[142] * fwd_Ea[142] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[27] -= q; /* C3H6 */
    wdot[29] += q; /* IXC3H7 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[27];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[97] -= dqdci;               /* dwdot[C3H6]/d[H] */
    J[99] += dqdci;               /* dwdot[IXC3H7]/d[H] */
    /* d()/d[C3H6] */
    dqdci =  + k_f*sc[2];
    J[947] -= dqdci;              /* dwdot[H]/d[C3H6] */
    J[972] -= dqdci;              /* dwdot[C3H6]/d[C3H6] */
    J[974] += dqdci;              /* dwdot[IXC3H7]/d[C3H6] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1217] -= dqdT;              /* dwdot[C3H6]/dT */
    J[1219] += dqdT;              /* dwdot[IXC3H7]/dT */

    /*reaction 144: NXC3H7 => CH3 + C2H4 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[30];
    k_f = prefactor_units[143] * fwd_A[143]
                * exp(fwd_beta[143] * tc[0] - activation_units[143] * fwd_Ea[143] * invT);
    dlnkfdT = fwd_beta[143] * invT + activation_units[143] * fwd_Ea[143] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[11] += q; /* CH3 */
    wdot[16] += q; /* C2H4 */
    wdot[30] -= q; /* NXC3H7 */
    /* d()/d[NXC3H7] */
    dqdci =  + k_f;
    J[1061] += dqdci;             /* dwdot[CH3]/d[NXC3H7] */
    J[1066] += dqdci;             /* dwdot[C2H4]/d[NXC3H7] */
    J[1080] -= dqdci;             /* dwdot[NXC3H7]/d[NXC3H7] */
    /* d()/dT */
    J[1201] += dqdT;              /* dwdot[CH3]/dT */
    J[1206] += dqdT;              /* dwdot[C2H4]/dT */
    J[1220] -= dqdT;              /* dwdot[NXC3H7]/dT */

    /*reaction 145: CH3 + C2H4 => NXC3H7 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[16];
    k_f = prefactor_units[144] * fwd_A[144]
                * exp(fwd_beta[144] * tc[0] - activation_units[144] * fwd_Ea[144] * invT);
    dlnkfdT = fwd_beta[144] * invT + activation_units[144] * fwd_Ea[144] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[11] -= q; /* CH3 */
    wdot[16] -= q; /* C2H4 */
    wdot[30] += q; /* NXC3H7 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[16];
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[401] -= dqdci;              /* dwdot[C2H4]/d[CH3] */
    J[415] += dqdci;              /* dwdot[NXC3H7]/d[CH3] */
    /* d()/d[C2H4] */
    dqdci =  + k_f*sc[11];
    J[571] -= dqdci;              /* dwdot[CH3]/d[C2H4] */
    J[576] -= dqdci;              /* dwdot[C2H4]/d[C2H4] */
    J[590] += dqdci;              /* dwdot[NXC3H7]/d[C2H4] */
    /* d()/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1206] -= dqdT;              /* dwdot[C2H4]/dT */
    J[1220] += dqdT;              /* dwdot[NXC3H7]/dT */

    /*reaction 146: NXC3H7 + HO2 => C3H8 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[30];
    k_f = prefactor_units[145] * fwd_A[145]
                * exp(fwd_beta[145] * tc[0] - activation_units[145] * fwd_Ea[145] * invT);
    dlnkfdT = fwd_beta[145] * invT + activation_units[145] * fwd_Ea[145] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[8] -= q; /* HO2 */
    wdot[30] -= q; /* NXC3H7 */
    wdot[31] += q; /* C3H8 */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[30];
    J[285] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[310] -= dqdci;              /* dwdot[NXC3H7]/d[HO2] */
    J[311] += dqdci;              /* dwdot[C3H8]/d[HO2] */
    /* d()/d[NXC3H7] */
    dqdci =  + k_f*sc[8];
    J[1055] += dqdci;             /* dwdot[O2]/d[NXC3H7] */
    J[1058] -= dqdci;             /* dwdot[HO2]/d[NXC3H7] */
    J[1080] -= dqdci;             /* dwdot[NXC3H7]/d[NXC3H7] */
    J[1081] += dqdci;             /* dwdot[C3H8]/d[NXC3H7] */
    /* d()/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */
    J[1220] -= dqdT;              /* dwdot[NXC3H7]/dT */
    J[1221] += dqdT;              /* dwdot[C3H8]/dT */

    /*reaction 147: NXC3H7 + O2 => C3H6 + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[30];
    k_f = prefactor_units[146] * fwd_A[146]
                * exp(fwd_beta[146] * tc[0] - activation_units[146] * fwd_Ea[146] * invT);
    dlnkfdT = fwd_beta[146] * invT + activation_units[146] * fwd_Ea[146] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[8] += q; /* HO2 */
    wdot[27] += q; /* C3H6 */
    wdot[30] -= q; /* NXC3H7 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[30];
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[183] += dqdci;              /* dwdot[HO2]/d[O2] */
    J[202] += dqdci;              /* dwdot[C3H6]/d[O2] */
    J[205] -= dqdci;              /* dwdot[NXC3H7]/d[O2] */
    /* d()/d[NXC3H7] */
    dqdci =  + k_f*sc[5];
    J[1055] -= dqdci;             /* dwdot[O2]/d[NXC3H7] */
    J[1058] += dqdci;             /* dwdot[HO2]/d[NXC3H7] */
    J[1077] += dqdci;             /* dwdot[C3H6]/d[NXC3H7] */
    J[1080] -= dqdci;             /* dwdot[NXC3H7]/d[NXC3H7] */
    /* d()/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */
    J[1217] += dqdT;              /* dwdot[C3H6]/dT */
    J[1220] -= dqdT;              /* dwdot[NXC3H7]/dT */

    /*reaction 148: NXC3H7 => H + C3H6 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[30];
    k_f = prefactor_units[147] * fwd_A[147]
                * exp(fwd_beta[147] * tc[0] - activation_units[147] * fwd_Ea[147] * invT);
    dlnkfdT = fwd_beta[147] * invT + activation_units[147] * fwd_Ea[147] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[27] += q; /* C3H6 */
    wdot[30] -= q; /* NXC3H7 */
    /* d()/d[NXC3H7] */
    dqdci =  + k_f;
    J[1052] += dqdci;             /* dwdot[H]/d[NXC3H7] */
    J[1077] += dqdci;             /* dwdot[C3H6]/d[NXC3H7] */
    J[1080] -= dqdci;             /* dwdot[NXC3H7]/d[NXC3H7] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1217] += dqdT;              /* dwdot[C3H6]/dT */
    J[1220] -= dqdT;              /* dwdot[NXC3H7]/dT */

    /*reaction 149: H + C3H6 => NXC3H7 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[27];
    k_f = prefactor_units[148] * fwd_A[148]
                * exp(fwd_beta[148] * tc[0] - activation_units[148] * fwd_Ea[148] * invT);
    dlnkfdT = fwd_beta[148] * invT + activation_units[148] * fwd_Ea[148] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[27] -= q; /* C3H6 */
    wdot[30] += q; /* NXC3H7 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[27];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[97] -= dqdci;               /* dwdot[C3H6]/d[H] */
    J[100] += dqdci;              /* dwdot[NXC3H7]/d[H] */
    /* d()/d[C3H6] */
    dqdci =  + k_f*sc[2];
    J[947] -= dqdci;              /* dwdot[H]/d[C3H6] */
    J[972] -= dqdci;              /* dwdot[C3H6]/d[C3H6] */
    J[975] += dqdci;              /* dwdot[NXC3H7]/d[C3H6] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1217] -= dqdT;              /* dwdot[C3H6]/dT */
    J[1220] += dqdT;              /* dwdot[NXC3H7]/dT */

    /*reaction 150: C3H8 + OH => NXC3H7 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[31];
    k_f = prefactor_units[149] * fwd_A[149]
                * exp(fwd_beta[149] * tc[0] - activation_units[149] * fwd_Ea[149] * invT);
    dlnkfdT = fwd_beta[149] * invT + activation_units[149] * fwd_Ea[149] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[30] += q; /* NXC3H7 */
    wdot[31] -= q; /* C3H8 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[31];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[111] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[135] += dqdci;              /* dwdot[NXC3H7]/d[OH] */
    J[136] -= dqdci;              /* dwdot[C3H8]/d[OH] */
    /* d()/d[C3H8] */
    dqdci =  + k_f*sc[3];
    J[1088] -= dqdci;             /* dwdot[OH]/d[C3H8] */
    J[1091] += dqdci;             /* dwdot[H2O]/d[C3H8] */
    J[1115] += dqdci;             /* dwdot[NXC3H7]/d[C3H8] */
    J[1116] -= dqdci;             /* dwdot[C3H8]/d[C3H8] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */
    J[1220] += dqdT;              /* dwdot[NXC3H7]/dT */
    J[1221] -= dqdT;              /* dwdot[C3H8]/dT */

    /*reaction 151: C3H8 + HO2 => NXC3H7 + H2O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[31];
    k_f = prefactor_units[150] * fwd_A[150]
                * exp(fwd_beta[150] * tc[0] - activation_units[150] * fwd_Ea[150] * invT);
    dlnkfdT = fwd_beta[150] * invT + activation_units[150] * fwd_Ea[150] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[7] += q; /* H2O2 */
    wdot[8] -= q; /* HO2 */
    wdot[30] += q; /* NXC3H7 */
    wdot[31] -= q; /* C3H8 */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[31];
    J[287] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[310] += dqdci;              /* dwdot[NXC3H7]/d[HO2] */
    J[311] -= dqdci;              /* dwdot[C3H8]/d[HO2] */
    /* d()/d[C3H8] */
    dqdci =  + k_f*sc[8];
    J[1092] += dqdci;             /* dwdot[H2O2]/d[C3H8] */
    J[1093] -= dqdci;             /* dwdot[HO2]/d[C3H8] */
    J[1115] += dqdci;             /* dwdot[NXC3H7]/d[C3H8] */
    J[1116] -= dqdci;             /* dwdot[C3H8]/d[C3H8] */
    /* d()/dT */
    J[1197] += dqdT;              /* dwdot[H2O2]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */
    J[1220] += dqdT;              /* dwdot[NXC3H7]/dT */
    J[1221] -= dqdT;              /* dwdot[C3H8]/dT */

    /*reaction 152: H + C3H8 <=> H2 + NXC3H7 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[31];
    k_f = prefactor_units[151] * fwd_A[151]
                * exp(fwd_beta[151] * tc[0] - activation_units[151] * fwd_Ea[151] * invT);
    dlnkfdT = fwd_beta[151] * invT + activation_units[151] * fwd_Ea[151] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[30];
    Kc = exp(g_RT[2] - g_RT[4] - g_RT[30] + g_RT[31]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[31]) + (h_RT[4] + h_RT[30]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[4] += q; /* H2 */
    wdot[30] += q; /* NXC3H7 */
    wdot[31] -= q; /* C3H8 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[31];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[74] += dqdci;               /* dwdot[H2]/d[H] */
    J[100] += dqdci;              /* dwdot[NXC3H7]/d[H] */
    J[101] -= dqdci;              /* dwdot[C3H8]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[30];
    J[142] -= dqdci;              /* dwdot[H]/d[H2] */
    J[144] += dqdci;              /* dwdot[H2]/d[H2] */
    J[170] += dqdci;              /* dwdot[NXC3H7]/d[H2] */
    J[171] -= dqdci;              /* dwdot[C3H8]/d[H2] */
    /* d()/d[NXC3H7] */
    dqdci =  - k_r*sc[4];
    J[1052] -= dqdci;             /* dwdot[H]/d[NXC3H7] */
    J[1054] += dqdci;             /* dwdot[H2]/d[NXC3H7] */
    J[1080] += dqdci;             /* dwdot[NXC3H7]/d[NXC3H7] */
    J[1081] -= dqdci;             /* dwdot[C3H8]/d[NXC3H7] */
    /* d()/d[C3H8] */
    dqdci =  + k_f*sc[2];
    J[1087] -= dqdci;             /* dwdot[H]/d[C3H8] */
    J[1089] += dqdci;             /* dwdot[H2]/d[C3H8] */
    J[1115] += dqdci;             /* dwdot[NXC3H7]/d[C3H8] */
    J[1116] -= dqdci;             /* dwdot[C3H8]/d[C3H8] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */
    J[1220] += dqdT;              /* dwdot[NXC3H7]/dT */
    J[1221] -= dqdT;              /* dwdot[C3H8]/dT */

    /*reaction 153: C3H8 + OH => IXC3H7 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[31];
    k_f = prefactor_units[152] * fwd_A[152]
                * exp(fwd_beta[152] * tc[0] - activation_units[152] * fwd_Ea[152] * invT);
    dlnkfdT = fwd_beta[152] * invT + activation_units[152] * fwd_Ea[152] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[6] += q; /* H2O */
    wdot[29] += q; /* IXC3H7 */
    wdot[31] -= q; /* C3H8 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[31];
    J[108] -= dqdci;              /* dwdot[OH]/d[OH] */
    J[111] += dqdci;              /* dwdot[H2O]/d[OH] */
    J[134] += dqdci;              /* dwdot[IXC3H7]/d[OH] */
    J[136] -= dqdci;              /* dwdot[C3H8]/d[OH] */
    /* d()/d[C3H8] */
    dqdci =  + k_f*sc[3];
    J[1088] -= dqdci;             /* dwdot[OH]/d[C3H8] */
    J[1091] += dqdci;             /* dwdot[H2O]/d[C3H8] */
    J[1114] += dqdci;             /* dwdot[IXC3H7]/d[C3H8] */
    J[1116] -= dqdci;             /* dwdot[C3H8]/d[C3H8] */
    /* d()/dT */
    J[1193] -= dqdT;              /* dwdot[OH]/dT */
    J[1196] += dqdT;              /* dwdot[H2O]/dT */
    J[1219] += dqdT;              /* dwdot[IXC3H7]/dT */
    J[1221] -= dqdT;              /* dwdot[C3H8]/dT */

    /*reaction 154: CH3 + C3H8 => CH4 + IXC3H7 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[31];
    k_f = prefactor_units[153] * fwd_A[153]
                * exp(fwd_beta[153] * tc[0] - activation_units[153] * fwd_Ea[153] * invT);
    dlnkfdT = fwd_beta[153] * invT + activation_units[153] * fwd_Ea[153] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[11] -= q; /* CH3 */
    wdot[15] += q; /* CH4 */
    wdot[29] += q; /* IXC3H7 */
    wdot[31] -= q; /* C3H8 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[31];
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[400] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[414] += dqdci;              /* dwdot[IXC3H7]/d[CH3] */
    J[416] -= dqdci;              /* dwdot[C3H8]/d[CH3] */
    /* d()/d[C3H8] */
    dqdci =  + k_f*sc[11];
    J[1096] -= dqdci;             /* dwdot[CH3]/d[C3H8] */
    J[1100] += dqdci;             /* dwdot[CH4]/d[C3H8] */
    J[1114] += dqdci;             /* dwdot[IXC3H7]/d[C3H8] */
    J[1116] -= dqdci;             /* dwdot[C3H8]/d[C3H8] */
    /* d()/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1205] += dqdT;              /* dwdot[CH4]/dT */
    J[1219] += dqdT;              /* dwdot[IXC3H7]/dT */
    J[1221] -= dqdT;              /* dwdot[C3H8]/dT */

    /*reaction 155: CH3 + C3H8 => CH4 + NXC3H7 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[11]*sc[31];
    k_f = prefactor_units[154] * fwd_A[154]
                * exp(fwd_beta[154] * tc[0] - activation_units[154] * fwd_Ea[154] * invT);
    dlnkfdT = fwd_beta[154] * invT + activation_units[154] * fwd_Ea[154] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[11] -= q; /* CH3 */
    wdot[15] += q; /* CH4 */
    wdot[30] += q; /* NXC3H7 */
    wdot[31] -= q; /* C3H8 */
    /* d()/d[CH3] */
    dqdci =  + k_f*sc[31];
    J[396] -= dqdci;              /* dwdot[CH3]/d[CH3] */
    J[400] += dqdci;              /* dwdot[CH4]/d[CH3] */
    J[415] += dqdci;              /* dwdot[NXC3H7]/d[CH3] */
    J[416] -= dqdci;              /* dwdot[C3H8]/d[CH3] */
    /* d()/d[C3H8] */
    dqdci =  + k_f*sc[11];
    J[1096] -= dqdci;             /* dwdot[CH3]/d[C3H8] */
    J[1100] += dqdci;             /* dwdot[CH4]/d[C3H8] */
    J[1115] += dqdci;             /* dwdot[NXC3H7]/d[C3H8] */
    J[1116] -= dqdci;             /* dwdot[C3H8]/d[C3H8] */
    /* d()/dT */
    J[1201] -= dqdT;              /* dwdot[CH3]/dT */
    J[1205] += dqdT;              /* dwdot[CH4]/dT */
    J[1220] += dqdT;              /* dwdot[NXC3H7]/dT */
    J[1221] -= dqdT;              /* dwdot[C3H8]/dT */

    /*reaction 156: C3H8 + O => IXC3H7 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[31];
    k_f = prefactor_units[155] * fwd_A[155]
                * exp(fwd_beta[155] * tc[0] - activation_units[155] * fwd_Ea[155] * invT);
    dlnkfdT = fwd_beta[155] * invT + activation_units[155] * fwd_Ea[155] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[3] += q; /* OH */
    wdot[29] += q; /* IXC3H7 */
    wdot[31] -= q; /* C3H8 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[31];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[38] += dqdci;               /* dwdot[OH]/d[O] */
    J[64] += dqdci;               /* dwdot[IXC3H7]/d[O] */
    J[66] -= dqdci;               /* dwdot[C3H8]/d[O] */
    /* d()/d[C3H8] */
    dqdci =  + k_f*sc[1];
    J[1086] -= dqdci;             /* dwdot[O]/d[C3H8] */
    J[1088] += dqdci;             /* dwdot[OH]/d[C3H8] */
    J[1114] += dqdci;             /* dwdot[IXC3H7]/d[C3H8] */
    J[1116] -= dqdci;             /* dwdot[C3H8]/d[C3H8] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1219] += dqdT;              /* dwdot[IXC3H7]/dT */
    J[1221] -= dqdT;              /* dwdot[C3H8]/dT */

    /*reaction 157: C3H8 + HO2 => IXC3H7 + H2O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[31];
    k_f = prefactor_units[156] * fwd_A[156]
                * exp(fwd_beta[156] * tc[0] - activation_units[156] * fwd_Ea[156] * invT);
    dlnkfdT = fwd_beta[156] * invT + activation_units[156] * fwd_Ea[156] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[7] += q; /* H2O2 */
    wdot[8] -= q; /* HO2 */
    wdot[29] += q; /* IXC3H7 */
    wdot[31] -= q; /* C3H8 */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[31];
    J[287] += dqdci;              /* dwdot[H2O2]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[309] += dqdci;              /* dwdot[IXC3H7]/d[HO2] */
    J[311] -= dqdci;              /* dwdot[C3H8]/d[HO2] */
    /* d()/d[C3H8] */
    dqdci =  + k_f*sc[8];
    J[1092] += dqdci;             /* dwdot[H2O2]/d[C3H8] */
    J[1093] -= dqdci;             /* dwdot[HO2]/d[C3H8] */
    J[1114] += dqdci;             /* dwdot[IXC3H7]/d[C3H8] */
    J[1116] -= dqdci;             /* dwdot[C3H8]/d[C3H8] */
    /* d()/dT */
    J[1197] += dqdT;              /* dwdot[H2O2]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */
    J[1219] += dqdT;              /* dwdot[IXC3H7]/dT */
    J[1221] -= dqdT;              /* dwdot[C3H8]/dT */

    /*reaction 158: C3H8 + O => NXC3H7 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[31];
    k_f = prefactor_units[157] * fwd_A[157]
                * exp(fwd_beta[157] * tc[0] - activation_units[157] * fwd_Ea[157] * invT);
    dlnkfdT = fwd_beta[157] * invT + activation_units[157] * fwd_Ea[157] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[1] -= q; /* O */
    wdot[3] += q; /* OH */
    wdot[30] += q; /* NXC3H7 */
    wdot[31] -= q; /* C3H8 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[31];
    J[36] -= dqdci;               /* dwdot[O]/d[O] */
    J[38] += dqdci;               /* dwdot[OH]/d[O] */
    J[65] += dqdci;               /* dwdot[NXC3H7]/d[O] */
    J[66] -= dqdci;               /* dwdot[C3H8]/d[O] */
    /* d()/d[C3H8] */
    dqdci =  + k_f*sc[1];
    J[1086] -= dqdci;             /* dwdot[O]/d[C3H8] */
    J[1088] += dqdci;             /* dwdot[OH]/d[C3H8] */
    J[1115] += dqdci;             /* dwdot[NXC3H7]/d[C3H8] */
    J[1116] -= dqdci;             /* dwdot[C3H8]/d[C3H8] */
    /* d()/dT */
    J[1191] -= dqdT;              /* dwdot[O]/dT */
    J[1193] += dqdT;              /* dwdot[OH]/dT */
    J[1220] += dqdT;              /* dwdot[NXC3H7]/dT */
    J[1221] -= dqdT;              /* dwdot[C3H8]/dT */

    /*reaction 159: C3H8 + O2 => IXC3H7 + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[31];
    k_f = prefactor_units[158] * fwd_A[158]
                * exp(fwd_beta[158] * tc[0] - activation_units[158] * fwd_Ea[158] * invT);
    dlnkfdT = fwd_beta[158] * invT + activation_units[158] * fwd_Ea[158] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[8] += q; /* HO2 */
    wdot[29] += q; /* IXC3H7 */
    wdot[31] -= q; /* C3H8 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[31];
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[183] += dqdci;              /* dwdot[HO2]/d[O2] */
    J[204] += dqdci;              /* dwdot[IXC3H7]/d[O2] */
    J[206] -= dqdci;              /* dwdot[C3H8]/d[O2] */
    /* d()/d[C3H8] */
    dqdci =  + k_f*sc[5];
    J[1090] -= dqdci;             /* dwdot[O2]/d[C3H8] */
    J[1093] += dqdci;             /* dwdot[HO2]/d[C3H8] */
    J[1114] += dqdci;             /* dwdot[IXC3H7]/d[C3H8] */
    J[1116] -= dqdci;             /* dwdot[C3H8]/d[C3H8] */
    /* d()/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1198] += dqdT;              /* dwdot[HO2]/dT */
    J[1219] += dqdT;              /* dwdot[IXC3H7]/dT */
    J[1221] -= dqdT;              /* dwdot[C3H8]/dT */

    /*reaction 160: IXC3H7 + HO2 => C3H8 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[8]*sc[29];
    k_f = prefactor_units[159] * fwd_A[159]
                * exp(fwd_beta[159] * tc[0] - activation_units[159] * fwd_Ea[159] * invT);
    dlnkfdT = fwd_beta[159] * invT + activation_units[159] * fwd_Ea[159] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[8] -= q; /* HO2 */
    wdot[29] -= q; /* IXC3H7 */
    wdot[31] += q; /* C3H8 */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[29];
    J[285] += dqdci;              /* dwdot[O2]/d[HO2] */
    J[288] -= dqdci;              /* dwdot[HO2]/d[HO2] */
    J[309] -= dqdci;              /* dwdot[IXC3H7]/d[HO2] */
    J[311] += dqdci;              /* dwdot[C3H8]/d[HO2] */
    /* d()/d[IXC3H7] */
    dqdci =  + k_f*sc[8];
    J[1020] += dqdci;             /* dwdot[O2]/d[IXC3H7] */
    J[1023] -= dqdci;             /* dwdot[HO2]/d[IXC3H7] */
    J[1044] -= dqdci;             /* dwdot[IXC3H7]/d[IXC3H7] */
    J[1046] += dqdci;             /* dwdot[C3H8]/d[IXC3H7] */
    /* d()/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1198] -= dqdT;              /* dwdot[HO2]/dT */
    J[1219] -= dqdT;              /* dwdot[IXC3H7]/dT */
    J[1221] += dqdT;              /* dwdot[C3H8]/dT */

    /*reaction 161: H + C3H8 => H2 + IXC3H7 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[31];
    k_f = prefactor_units[160] * fwd_A[160]
                * exp(fwd_beta[160] * tc[0] - activation_units[160] * fwd_Ea[160] * invT);
    dlnkfdT = fwd_beta[160] * invT + activation_units[160] * fwd_Ea[160] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] -= q; /* H */
    wdot[4] += q; /* H2 */
    wdot[29] += q; /* IXC3H7 */
    wdot[31] -= q; /* C3H8 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[31];
    J[72] -= dqdci;               /* dwdot[H]/d[H] */
    J[74] += dqdci;               /* dwdot[H2]/d[H] */
    J[99] += dqdci;               /* dwdot[IXC3H7]/d[H] */
    J[101] -= dqdci;              /* dwdot[C3H8]/d[H] */
    /* d()/d[C3H8] */
    dqdci =  + k_f*sc[2];
    J[1087] -= dqdci;             /* dwdot[H]/d[C3H8] */
    J[1089] += dqdci;             /* dwdot[H2]/d[C3H8] */
    J[1114] += dqdci;             /* dwdot[IXC3H7]/d[C3H8] */
    J[1116] -= dqdci;             /* dwdot[C3H8]/d[C3H8] */
    /* d()/dT */
    J[1192] -= dqdT;              /* dwdot[H]/dT */
    J[1194] += dqdT;              /* dwdot[H2]/dT */
    J[1219] += dqdT;              /* dwdot[IXC3H7]/dT */
    J[1221] -= dqdT;              /* dwdot[C3H8]/dT */

    /*reaction 162: H2 + IXC3H7 => H + C3H8 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[29];
    k_f = prefactor_units[161] * fwd_A[161]
                * exp(fwd_beta[161] * tc[0] - activation_units[161] * fwd_Ea[161] * invT);
    dlnkfdT = fwd_beta[161] * invT + activation_units[161] * fwd_Ea[161] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[2] += q; /* H */
    wdot[4] -= q; /* H2 */
    wdot[29] -= q; /* IXC3H7 */
    wdot[31] += q; /* C3H8 */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[29];
    J[142] += dqdci;              /* dwdot[H]/d[H2] */
    J[144] -= dqdci;              /* dwdot[H2]/d[H2] */
    J[169] -= dqdci;              /* dwdot[IXC3H7]/d[H2] */
    J[171] += dqdci;              /* dwdot[C3H8]/d[H2] */
    /* d()/d[IXC3H7] */
    dqdci =  + k_f*sc[4];
    J[1017] += dqdci;             /* dwdot[H]/d[IXC3H7] */
    J[1019] -= dqdci;             /* dwdot[H2]/d[IXC3H7] */
    J[1044] -= dqdci;             /* dwdot[IXC3H7]/d[IXC3H7] */
    J[1046] += dqdci;             /* dwdot[C3H8]/d[IXC3H7] */
    /* d()/dT */
    J[1192] += dqdT;              /* dwdot[H]/dT */
    J[1194] -= dqdT;              /* dwdot[H2]/dT */
    J[1219] -= dqdT;              /* dwdot[IXC3H7]/dT */
    J[1221] += dqdT;              /* dwdot[C3H8]/dT */

    /*reaction 163: C3H5O => C2H3 + CH2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[28];
    k_f = prefactor_units[162] * fwd_A[162]
                * exp(fwd_beta[162] * tc[0] - activation_units[162] * fwd_Ea[162] * invT);
    dlnkfdT = fwd_beta[162] * invT + activation_units[162] * fwd_Ea[162] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[10] += q; /* CH2O */
    wdot[24] += q; /* C2H3 */
    wdot[28] -= q; /* C3H5O */
    /* d()/d[C3H5O] */
    dqdci =  + k_f;
    J[990] += dqdci;              /* dwdot[CH2O]/d[C3H5O] */
    J[1004] += dqdci;             /* dwdot[C2H3]/d[C3H5O] */
    J[1008] -= dqdci;             /* dwdot[C3H5O]/d[C3H5O] */
    /* d()/dT */
    J[1200] += dqdT;              /* dwdot[CH2O]/dT */
    J[1214] += dqdT;              /* dwdot[C2H3]/dT */
    J[1218] -= dqdT;              /* dwdot[C3H5O]/dT */

    /*reaction 164: IXC3H7O2 => IXC3H7 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[32];
    k_f = prefactor_units[163] * fwd_A[163]
                * exp(fwd_beta[163] * tc[0] - activation_units[163] * fwd_Ea[163] * invT);
    dlnkfdT = fwd_beta[163] * invT + activation_units[163] * fwd_Ea[163] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[29] += q; /* IXC3H7 */
    wdot[32] -= q; /* IXC3H7O2 */
    /* d()/d[IXC3H7O2] */
    dqdci =  + k_f;
    J[1125] += dqdci;             /* dwdot[O2]/d[IXC3H7O2] */
    J[1149] += dqdci;             /* dwdot[IXC3H7]/d[IXC3H7O2] */
    J[1152] -= dqdci;             /* dwdot[IXC3H7O2]/d[IXC3H7O2] */
    /* d()/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1219] += dqdT;              /* dwdot[IXC3H7]/dT */
    J[1222] -= dqdT;              /* dwdot[IXC3H7O2]/dT */

    /*reaction 165: IXC3H7 + O2 => IXC3H7O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[29];
    k_f = prefactor_units[164] * fwd_A[164]
                * exp(fwd_beta[164] * tc[0] - activation_units[164] * fwd_Ea[164] * invT);
    dlnkfdT = fwd_beta[164] * invT + activation_units[164] * fwd_Ea[164] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[29] -= q; /* IXC3H7 */
    wdot[32] += q; /* IXC3H7O2 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[29];
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[204] -= dqdci;              /* dwdot[IXC3H7]/d[O2] */
    J[207] += dqdci;              /* dwdot[IXC3H7O2]/d[O2] */
    /* d()/d[IXC3H7] */
    dqdci =  + k_f*sc[5];
    J[1020] -= dqdci;             /* dwdot[O2]/d[IXC3H7] */
    J[1044] -= dqdci;             /* dwdot[IXC3H7]/d[IXC3H7] */
    J[1047] += dqdci;             /* dwdot[IXC3H7O2]/d[IXC3H7] */
    /* d()/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1219] -= dqdT;              /* dwdot[IXC3H7]/dT */
    J[1222] += dqdT;              /* dwdot[IXC3H7O2]/dT */

    /*reaction 166: NXC3H7O2 => NXC3H7 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[33];
    k_f = prefactor_units[165] * fwd_A[165]
                * exp(fwd_beta[165] * tc[0] - activation_units[165] * fwd_Ea[165] * invT);
    dlnkfdT = fwd_beta[165] * invT + activation_units[165] * fwd_Ea[165] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[30] += q; /* NXC3H7 */
    wdot[33] -= q; /* NXC3H7O2 */
    /* d()/d[NXC3H7O2] */
    dqdci =  + k_f;
    J[1160] += dqdci;             /* dwdot[O2]/d[NXC3H7O2] */
    J[1185] += dqdci;             /* dwdot[NXC3H7]/d[NXC3H7O2] */
    J[1188] -= dqdci;             /* dwdot[NXC3H7O2]/d[NXC3H7O2] */
    /* d()/dT */
    J[1195] += dqdT;              /* dwdot[O2]/dT */
    J[1220] += dqdT;              /* dwdot[NXC3H7]/dT */
    J[1223] -= dqdT;              /* dwdot[NXC3H7O2]/dT */

    /*reaction 167: NXC3H7 + O2 => NXC3H7O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[5]*sc[30];
    k_f = prefactor_units[166] * fwd_A[166]
                * exp(fwd_beta[166] * tc[0] - activation_units[166] * fwd_Ea[166] * invT);
    dlnkfdT = fwd_beta[166] * invT + activation_units[166] * fwd_Ea[166] * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[5] -= q; /* O2 */
    wdot[30] -= q; /* NXC3H7 */
    wdot[33] += q; /* NXC3H7O2 */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[30];
    J[180] -= dqdci;              /* dwdot[O2]/d[O2] */
    J[205] -= dqdci;              /* dwdot[NXC3H7]/d[O2] */
    J[208] += dqdci;              /* dwdot[NXC3H7O2]/d[O2] */
    /* d()/d[NXC3H7] */
    dqdci =  + k_f*sc[5];
    J[1055] -= dqdci;             /* dwdot[O2]/d[NXC3H7] */
    J[1080] -= dqdci;             /* dwdot[NXC3H7]/d[NXC3H7] */
    J[1083] += dqdci;             /* dwdot[NXC3H7O2]/d[NXC3H7] */
    /* d()/dT */
    J[1195] -= dqdT;              /* dwdot[O2]/dT */
    J[1220] -= dqdT;              /* dwdot[NXC3H7]/dT */
    J[1223] += dqdT;              /* dwdot[NXC3H7O2]/dT */

    amrex::Real c_R[34], dcRdT[34], e_RT[34];
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
    for (int k = 0; k < 34; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[1190+k];
    }

    amrex::Real cmixinv = 1.0/cmix;
    amrex::Real tmp1 = ehmix*cmixinv;
    amrex::Real tmp3 = cmixinv*T;
    amrex::Real tmp2 = tmp1*tmp3;
    amrex::Real dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 34; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 34; ++m) {
            dehmixdc += eh_RT[m]*J[k*35+m];
        }
        J[k*35+34] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[1224] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;
}

#endif
/* Transport function declarations  */


void egtransetLENIMC(int* LENIMC ) {
    *LENIMC = 138;}


void egtransetLENRMC(int* LENRMC ) {
    *LENRMC = 23154;}


void egtransetNO(int* NO ) {
    *NO = 4;}


void egtransetKK(int* KK ) {
    *KK = 34;}


void egtransetNLITE(int* NLITE ) {
    *NLITE = 2;}


/*Patm in ergs/cm3 */
void egtransetPATM(amrex::Real* PATM) {
    *PATM =   0.1013250000000000E+07;}


/*the molecular weights in g/mol */
void egtransetWT(amrex::Real* WT ) {
    WT[0] = 2.80134000E+01;
    WT[1] = 1.59994000E+01;
    WT[2] = 1.00797000E+00;
    WT[3] = 1.70073700E+01;
    WT[4] = 2.01594000E+00;
    WT[5] = 3.19988000E+01;
    WT[6] = 1.80153400E+01;
    WT[7] = 3.40147400E+01;
    WT[8] = 3.30067700E+01;
    WT[9] = 1.40270900E+01;
    WT[10] = 3.00264900E+01;
    WT[11] = 1.50350600E+01;
    WT[12] = 2.80105500E+01;
    WT[13] = 3.10344600E+01;
    WT[14] = 2.90621500E+01;
    WT[15] = 1.60430300E+01;
    WT[16] = 2.80541800E+01;
    WT[17] = 3.00701200E+01;
    WT[18] = 4.40099500E+01;
    WT[19] = 2.90185200E+01;
    WT[20] = 4.70338600E+01;
    WT[21] = 4.80418300E+01;
    WT[22] = 2.60382400E+01;
    WT[23] = 4.10296700E+01;
    WT[24] = 2.70462100E+01;
    WT[25] = 4.30456100E+01;
    WT[26] = 4.10733000E+01;
    WT[27] = 4.20812700E+01;
    WT[28] = 5.70727000E+01;
    WT[29] = 4.30892400E+01;
    WT[30] = 4.30892400E+01;
    WT[31] = 4.40972100E+01;
    WT[32] = 7.50880400E+01;
    WT[33] = 7.50880400E+01;
}


/*the lennard-jones potential well depth eps/kb in K */
void egtransetEPS(amrex::Real* EPS ) {
    EPS[0] = 9.75300000E+01;
    EPS[1] = 8.00000000E+01;
    EPS[2] = 1.45000000E+02;
    EPS[3] = 8.00000000E+01;
    EPS[4] = 3.80000000E+01;
    EPS[5] = 1.07400000E+02;
    EPS[6] = 5.72400000E+02;
    EPS[7] = 1.07400000E+02;
    EPS[8] = 1.07400000E+02;
    EPS[9] = 1.44000000E+02;
    EPS[10] = 4.98000000E+02;
    EPS[11] = 1.44000000E+02;
    EPS[12] = 9.81000000E+01;
    EPS[13] = 4.17000000E+02;
    EPS[14] = 2.47500000E+02;
    EPS[15] = 1.41400000E+02;
    EPS[16] = 2.38400000E+02;
    EPS[17] = 2.47500000E+02;
    EPS[18] = 2.44000000E+02;
    EPS[19] = 4.98000000E+02;
    EPS[20] = 4.81800000E+02;
    EPS[21] = 4.81800000E+02;
    EPS[22] = 2.65300000E+02;
    EPS[23] = 1.50000000E+02;
    EPS[24] = 2.65300000E+02;
    EPS[25] = 4.36000000E+02;
    EPS[26] = 3.16000000E+02;
    EPS[27] = 3.07800000E+02;
    EPS[28] = 4.11000000E+02;
    EPS[29] = 3.03400000E+02;
    EPS[30] = 3.03400000E+02;
    EPS[31] = 3.03400000E+02;
    EPS[32] = 4.59500000E+02;
    EPS[33] = 4.81500000E+02;
}


/*the lennard-jones collision diameter in Angstroms */
void egtransetSIG(amrex::Real* SIG ) {
    SIG[0] = 3.62100000E+00;
    SIG[1] = 2.75000000E+00;
    SIG[2] = 2.05000000E+00;
    SIG[3] = 2.75000000E+00;
    SIG[4] = 2.92000000E+00;
    SIG[5] = 3.45800000E+00;
    SIG[6] = 2.60500000E+00;
    SIG[7] = 3.45800000E+00;
    SIG[8] = 3.45800000E+00;
    SIG[9] = 3.80000000E+00;
    SIG[10] = 3.59000000E+00;
    SIG[11] = 3.80000000E+00;
    SIG[12] = 3.65000000E+00;
    SIG[13] = 3.69000000E+00;
    SIG[14] = 4.35000000E+00;
    SIG[15] = 3.74600000E+00;
    SIG[16] = 3.49600000E+00;
    SIG[17] = 4.35000000E+00;
    SIG[18] = 3.76300000E+00;
    SIG[19] = 3.59000000E+00;
    SIG[20] = 3.62600000E+00;
    SIG[21] = 3.62600000E+00;
    SIG[22] = 3.72100000E+00;
    SIG[23] = 2.50000000E+00;
    SIG[24] = 3.72100000E+00;
    SIG[25] = 3.97000000E+00;
    SIG[26] = 4.22000000E+00;
    SIG[27] = 4.14000000E+00;
    SIG[28] = 4.82000000E+00;
    SIG[29] = 4.81000000E+00;
    SIG[30] = 4.81000000E+00;
    SIG[31] = 4.81000000E+00;
    SIG[32] = 5.03600000E+00;
    SIG[33] = 4.99700000E+00;
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
    DIP[13] = 1.70000000E+00;
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
    DIP[29] = 0.00000000E+00;
    DIP[30] = 0.00000000E+00;
    DIP[31] = 0.00000000E+00;
    DIP[32] = 1.70000000E+00;
    DIP[33] = 1.70000000E+00;
}


/*the polarizability in cubic Angstroms */
void egtransetPOL(amrex::Real* POL ) {
    POL[0] = 1.76000000E+00;
    POL[1] = 0.00000000E+00;
    POL[2] = 0.00000000E+00;
    POL[3] = 0.00000000E+00;
    POL[4] = 7.90000000E-01;
    POL[5] = 1.60000000E+00;
    POL[6] = 0.00000000E+00;
    POL[7] = 0.00000000E+00;
    POL[8] = 0.00000000E+00;
    POL[9] = 0.00000000E+00;
    POL[10] = 0.00000000E+00;
    POL[11] = 0.00000000E+00;
    POL[12] = 1.95000000E+00;
    POL[13] = 0.00000000E+00;
    POL[14] = 0.00000000E+00;
    POL[15] = 2.60000000E+00;
    POL[16] = 0.00000000E+00;
    POL[17] = 0.00000000E+00;
    POL[18] = 2.65000000E+00;
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
    POL[30] = 0.00000000E+00;
    POL[31] = 0.00000000E+00;
    POL[32] = 0.00000000E+00;
    POL[33] = 0.00000000E+00;
}


/*the rotational relaxation collision number at 298 K */
void egtransetZROT(amrex::Real* ZROT ) {
    ZROT[0] = 4.00000000E+00;
    ZROT[1] = 0.00000000E+00;
    ZROT[2] = 0.00000000E+00;
    ZROT[3] = 0.00000000E+00;
    ZROT[4] = 2.80000000E+02;
    ZROT[5] = 3.80000000E+00;
    ZROT[6] = 4.00000000E+00;
    ZROT[7] = 3.80000000E+00;
    ZROT[8] = 1.00000000E+00;
    ZROT[9] = 0.00000000E+00;
    ZROT[10] = 2.00000000E+00;
    ZROT[11] = 0.00000000E+00;
    ZROT[12] = 1.80000000E+00;
    ZROT[13] = 2.00000000E+00;
    ZROT[14] = 1.50000000E+00;
    ZROT[15] = 1.30000000E+01;
    ZROT[16] = 1.50000000E+00;
    ZROT[17] = 1.50000000E+00;
    ZROT[18] = 2.10000000E+00;
    ZROT[19] = 0.00000000E+00;
    ZROT[20] = 1.00000000E+00;
    ZROT[21] = 1.00000000E+00;
    ZROT[22] = 2.50000000E+00;
    ZROT[23] = 1.00000000E+00;
    ZROT[24] = 1.00000000E+00;
    ZROT[25] = 2.00000000E+00;
    ZROT[26] = 1.00000000E+00;
    ZROT[27] = 1.00000000E+00;
    ZROT[28] = 1.00000000E+00;
    ZROT[29] = 1.00000000E+00;
    ZROT[30] = 1.00000000E+00;
    ZROT[31] = 1.00000000E+00;
    ZROT[32] = 1.00000000E+00;
    ZROT[33] = 1.00000000E+00;
}


/*0: monoatomic, 1: linear, 2: nonlinear */
void egtransetNLIN(int* NLIN) {
    NLIN[0] = 1;
    NLIN[1] = 0;
    NLIN[2] = 0;
    NLIN[3] = 1;
    NLIN[4] = 1;
    NLIN[5] = 1;
    NLIN[6] = 2;
    NLIN[7] = 2;
    NLIN[8] = 2;
    NLIN[9] = 1;
    NLIN[10] = 2;
    NLIN[11] = 1;
    NLIN[12] = 1;
    NLIN[13] = 2;
    NLIN[14] = 2;
    NLIN[15] = 2;
    NLIN[16] = 2;
    NLIN[17] = 2;
    NLIN[18] = 1;
    NLIN[19] = 2;
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
}


/*Poly fits for the viscosities, dim NO*KK */
void egtransetCOFETA(amrex::Real* COFETA) {
    COFETA[0] = -1.65695594E+01;
    COFETA[1] = 2.39056562E+00;
    COFETA[2] = -2.34558144E-01;
    COFETA[3] = 1.05024037E-02;
    COFETA[4] = -1.50926240E+01;
    COFETA[5] = 1.92606504E+00;
    COFETA[6] = -1.73487476E-01;
    COFETA[7] = 7.82572931E-03;
    COFETA[8] = -2.04078397E+01;
    COFETA[9] = 3.65436395E+00;
    COFETA[10] = -3.98339635E-01;
    COFETA[11] = 1.75883009E-02;
    COFETA[12] = -1.50620763E+01;
    COFETA[13] = 1.92606504E+00;
    COFETA[14] = -1.73487476E-01;
    COFETA[15] = 7.82572931E-03;
    COFETA[16] = -1.38347699E+01;
    COFETA[17] = 1.00106621E+00;
    COFETA[18] = -4.98105694E-02;
    COFETA[19] = 2.31450475E-03;
    COFETA[20] = -1.71618309E+01;
    COFETA[21] = 2.68036374E+00;
    COFETA[22] = -2.72570227E-01;
    COFETA[23] = 1.21650964E-02;
    COFETA[24] = -1.05420863E+01;
    COFETA[25] = -1.37777096E+00;
    COFETA[26] = 4.20502308E-01;
    COFETA[27] = -2.40627230E-02;
    COFETA[28] = -1.71312832E+01;
    COFETA[29] = 2.68036374E+00;
    COFETA[30] = -2.72570227E-01;
    COFETA[31] = 1.21650964E-02;
    COFETA[32] = -1.71463238E+01;
    COFETA[33] = 2.68036374E+00;
    COFETA[34] = -2.72570227E-01;
    COFETA[35] = 1.21650964E-02;
    COFETA[36] = -2.02663469E+01;
    COFETA[37] = 3.63241793E+00;
    COFETA[38] = -3.95581049E-01;
    COFETA[39] = 1.74725495E-02;
    COFETA[40] = -1.98330577E+01;
    COFETA[41] = 2.69480162E+00;
    COFETA[42] = -1.65880845E-01;
    COFETA[43] = 3.14504769E-03;
    COFETA[44] = -2.02316497E+01;
    COFETA[45] = 3.63241793E+00;
    COFETA[46] = -3.95581049E-01;
    COFETA[47] = 1.74725495E-02;
    COFETA[48] = -1.66188336E+01;
    COFETA[49] = 2.40307799E+00;
    COFETA[50] = -2.36167638E-01;
    COFETA[51] = 1.05714061E-02;
    COFETA[52] = -1.99945919E+01;
    COFETA[53] = 2.86923313E+00;
    COFETA[54] = -2.03325661E-01;
    COFETA[55] = 5.39056989E-03;
    COFETA[56] = -2.45602637E+01;
    COFETA[57] = 5.15878990E+00;
    COFETA[58] = -5.75274341E-01;
    COFETA[59] = 2.44975136E-02;
    COFETA[60] = -2.00094664E+01;
    COFETA[61] = 3.57220167E+00;
    COFETA[62] = -3.87936446E-01;
    COFETA[63] = 1.71483254E-02;
    COFETA[64] = -2.39690472E+01;
    COFETA[65] = 5.11436059E+00;
    COFETA[66] = -5.71999954E-01;
    COFETA[67] = 2.44581334E-02;
    COFETA[68] = -2.45432160E+01;
    COFETA[69] = 5.15878990E+00;
    COFETA[70] = -5.75274341E-01;
    COFETA[71] = 2.44975136E-02;
    COFETA[72] = -2.40014975E+01;
    COFETA[73] = 5.14359547E+00;
    COFETA[74] = -5.74269731E-01;
    COFETA[75] = 2.44937679E-02;
    COFETA[76] = -1.98501306E+01;
    COFETA[77] = 2.69480162E+00;
    COFETA[78] = -1.65880845E-01;
    COFETA[79] = 3.14504769E-03;
    COFETA[80] = -2.03725491E+01;
    COFETA[81] = 3.03946431E+00;
    COFETA[82] = -2.16994867E-01;
    COFETA[83] = 5.61394012E-03;
    COFETA[84] = -2.03619469E+01;
    COFETA[85] = 3.03946431E+00;
    COFETA[86] = -2.16994867E-01;
    COFETA[87] = 5.61394012E-03;
    COFETA[88] = -2.47697856E+01;
    COFETA[89] = 5.30039568E+00;
    COFETA[90] = -5.89273639E-01;
    COFETA[91] = 2.49261407E-02;
    COFETA[92] = -1.92183831E+01;
    COFETA[93] = 3.75164499E+00;
    COFETA[94] = -4.10390993E-01;
    COFETA[95] = 1.80861665E-02;
    COFETA[96] = -2.47507953E+01;
    COFETA[97] = 5.30039568E+00;
    COFETA[98] = -5.89273639E-01;
    COFETA[99] = 2.49261407E-02;
    COFETA[100] = -2.23277173E+01;
    COFETA[101] = 3.86433912E+00;
    COFETA[102] = -3.41553983E-01;
    COFETA[103] = 1.17083447E-02;
    COFETA[104] = -2.50402232E+01;
    COFETA[105] = 5.25451220E+00;
    COFETA[106] = -5.67228955E-01;
    COFETA[107] = 2.33156489E-02;
    COFETA[108] = -2.49727893E+01;
    COFETA[109] = 5.27067543E+00;
    COFETA[110] = -5.71909526E-01;
    COFETA[111] = 2.36230940E-02;
    COFETA[112] = -2.34114208E+01;
    COFETA[113] = 4.27461263E+00;
    COFETA[114] = -4.04617202E-01;
    COFETA[115] = 1.48351999E-02;
    COFETA[116] = -2.52462994E+01;
    COFETA[117] = 5.27749097E+00;
    COFETA[118] = -5.74219215E-01;
    COFETA[119] = 2.37811608E-02;
    COFETA[120] = -2.52462994E+01;
    COFETA[121] = 5.27749097E+00;
    COFETA[122] = -5.74219215E-01;
    COFETA[123] = 2.37811608E-02;
    COFETA[124] = -2.52347378E+01;
    COFETA[125] = 5.27749097E+00;
    COFETA[126] = -5.74219215E-01;
    COFETA[127] = 2.37811608E-02;
    COFETA[128] = -2.09096796E+01;
    COFETA[129] = 3.13792991E+00;
    COFETA[130] = -2.36145338E-01;
    COFETA[131] = 6.70688017E-03;
    COFETA[132] = -2.00330714E+01;
    COFETA[133] = 2.73318358E+00;
    COFETA[134] = -1.75653565E-01;
    COFETA[135] = 3.77126610E-03;
}


/*Poly fits for the conductivities, dim NO*KK */
void egtransetCOFLAM(amrex::Real* COFLAM) {
    COFLAM[0] = 1.29306625E+01;
    COFLAM[1] = -3.52819393E+00;
    COFLAM[2] = 6.45501923E-01;
    COFLAM[3] = -3.19376675E-02;
    COFLAM[4] = 1.69267361E+00;
    COFLAM[5] = 1.92606504E+00;
    COFLAM[6] = -1.73487476E-01;
    COFLAM[7] = 7.82572931E-03;
    COFLAM[8] = -8.57929284E-01;
    COFLAM[9] = 3.65436395E+00;
    COFLAM[10] = -3.98339635E-01;
    COFLAM[11] = 1.75883009E-02;
    COFLAM[12] = 1.00130995E+01;
    COFLAM[13] = -1.54006603E+00;
    COFLAM[14] = 3.02750813E-01;
    COFLAM[15] = -1.29859622E-02;
    COFLAM[16] = 9.13734407E+00;
    COFLAM[17] = -4.36833375E-01;
    COFLAM[18] = 1.12981765E-01;
    COFLAM[19] = -2.54610412E-03;
    COFLAM[20] = -1.93888086E+00;
    COFLAM[21] = 2.89244157E+00;
    COFLAM[22] = -2.71258557E-01;
    COFLAM[23] = 1.15340544E-02;
    COFLAM[24] = 2.35086340E+01;
    COFLAM[25] = -9.05997330E+00;
    COFLAM[26] = 1.54816025E+00;
    COFLAM[27] = -7.71639384E-02;
    COFLAM[28] = 8.99254403E-01;
    COFLAM[29] = 1.32506478E+00;
    COFLAM[30] = 1.81955930E-02;
    COFLAM[31] = -4.46691285E-03;
    COFLAM[32] = 4.76619374E+00;
    COFLAM[33] = -4.27755023E-01;
    COFLAM[34] = 2.68283274E-01;
    COFLAM[35] = -1.65411221E-02;
    COFLAM[36] = 1.83932434E+01;
    COFLAM[37] = -6.25754971E+00;
    COFLAM[38] = 1.09378595E+00;
    COFLAM[39] = -5.49435895E-02;
    COFLAM[40] = 4.84793117E+00;
    COFLAM[41] = -2.13603179E+00;
    COFLAM[42] = 6.99833906E-01;
    COFLAM[43] = -4.38222985E-02;
    COFLAM[44] = 1.36305859E+01;
    COFLAM[45] = -4.45086721E+00;
    COFLAM[46] = 8.75862692E-01;
    COFLAM[47] = -4.60275277E-02;
    COFLAM[48] = 1.15794127E+01;
    COFLAM[49] = -3.02088602E+00;
    COFLAM[50] = 5.82091962E-01;
    COFLAM[51] = -2.93406692E-02;
    COFLAM[52] = -6.14587484E+00;
    COFLAM[53] = 2.47428533E+00;
    COFLAM[54] = 6.44004931E-02;
    COFLAM[55] = -1.45368612E-02;
    COFLAM[56] = -1.32867451E+01;
    COFLAM[57] = 5.82553048E+00;
    COFLAM[58] = -4.38733487E-01;
    COFLAM[59] = 1.02209771E-02;
    COFLAM[60] = 1.29622480E+01;
    COFLAM[61] = -4.85747192E+00;
    COFLAM[62] = 1.02918185E+00;
    COFLAM[63] = -5.69931976E-02;
    COFLAM[64] = -1.34447168E+01;
    COFLAM[65] = 6.12380203E+00;
    COFLAM[66] = -4.86657425E-01;
    COFLAM[67] = 1.24614518E-02;
    COFLAM[68] = -1.44773293E+01;
    COFLAM[69] = 6.20799727E+00;
    COFLAM[70] = -4.66686188E-01;
    COFLAM[71] = 1.03037078E-02;
    COFLAM[72] = -1.15552013E+01;
    COFLAM[73] = 5.97444378E+00;
    COFLAM[74] = -5.83493959E-01;
    COFLAM[75] = 2.11390997E-02;
    COFLAM[76] = 5.62029187E+00;
    COFLAM[77] = -1.91574800E+00;
    COFLAM[78] = 5.90225152E-01;
    COFLAM[79] = -3.57616080E-02;
    COFLAM[80] = 5.28056570E+00;
    COFLAM[81] = -1.92758693E+00;
    COFLAM[82] = 6.29141169E-01;
    COFLAM[83] = -3.87203309E-02;
    COFLAM[84] = -5.85803496E+00;
    COFLAM[85] = 2.80606054E+00;
    COFLAM[86] = -2.49415114E-02;
    COFLAM[87] = -8.88253132E-03;
    COFLAM[88] = -9.20687365E+00;
    COFLAM[89] = 5.13028609E+00;
    COFLAM[90] = -4.67868863E-01;
    COFLAM[91] = 1.64674383E-02;
    COFLAM[92] = -5.31228079E-01;
    COFLAM[93] = 2.21998671E+00;
    COFLAM[94] = -1.07502048E-01;
    COFLAM[95] = 1.06754974E-03;
    COFLAM[96] = -1.36861349E+01;
    COFLAM[97] = 6.35261898E+00;
    COFLAM[98] = -5.53497228E-01;
    COFLAM[99] = 1.70958139E-02;
    COFLAM[100] = -6.27418855E+00;
    COFLAM[101] = 2.90468350E+00;
    COFLAM[102] = -4.35101878E-02;
    COFLAM[103] = -7.77942889E-03;
    COFLAM[104] = -2.08328673E+01;
    COFLAM[105] = 9.07593204E+00;
    COFLAM[106] = -8.93990863E-01;
    COFLAM[107] = 3.11142957E-02;
    COFLAM[108] = -1.54410770E+01;
    COFLAM[109] = 6.67114766E+00;
    COFLAM[110] = -5.37137624E-01;
    COFLAM[111] = 1.38051704E-02;
    COFLAM[112] = -1.68188905E+01;
    COFLAM[113] = 7.10377112E+00;
    COFLAM[114] = -6.04738536E-01;
    COFLAM[115] = 1.71995390E-02;
    COFLAM[116] = -1.69714063E+01;
    COFLAM[117] = 7.44003427E+00;
    COFLAM[118] = -6.71735783E-01;
    COFLAM[119] = 2.12052962E-02;
    COFLAM[120] = -1.59315280E+01;
    COFLAM[121] = 7.00975811E+00;
    COFLAM[122] = -6.12923034E-01;
    COFLAM[123] = 1.85439169E-02;
    COFLAM[124] = -1.96210407E+01;
    COFLAM[125] = 8.33592433E+00;
    COFLAM[126] = -7.66926846E-01;
    COFLAM[127] = 2.44576880E-02;
    COFLAM[128] = -1.19772854E+01;
    COFLAM[129] = 5.01880327E+00;
    COFLAM[130] = -3.06683992E-01;
    COFLAM[131] = 3.14889976E-03;
    COFLAM[132] = -9.04441654E+00;
    COFLAM[133] = 3.66855955E+00;
    COFLAM[134] = -1.04994889E-01;
    COFLAM[135] = -6.68408235E-03;
}


/*Poly fits for the diffusion coefficients, dim NO*KK*KK */
void egtransetCOFD(amrex::Real* COFD) {
    COFD[0] = -1.49828430E+01;
    COFD[1] = 3.25781069E+00;
    COFD[2] = -2.12199367E-01;
    COFD[3] = 9.36657283E-03;
    COFD[4] = -1.40756935E+01;
    COFD[5] = 3.07549274E+00;
    COFD[6] = -1.88889344E-01;
    COFD[7] = 8.37152866E-03;
    COFD[8] = -1.42894441E+01;
    COFD[9] = 3.67490723E+00;
    COFD[10] = -2.65114792E-01;
    COFD[11] = 1.16092671E-02;
    COFD[12] = -1.40949196E+01;
    COFD[13] = 3.07549274E+00;
    COFD[14] = -1.88889344E-01;
    COFD[15] = 8.37152866E-03;
    COFD[16] = -1.16906297E+01;
    COFD[17] = 2.47469981E+00;
    COFD[18] = -1.10436257E-01;
    COFD[19] = 4.95273813E-03;
    COFD[20] = -1.52414485E+01;
    COFD[21] = 3.35922578E+00;
    COFD[22] = -2.25181399E-01;
    COFD[23] = 9.92132878E-03;
    COFD[24] = -2.10643259E+01;
    COFD[25] = 5.53614847E+00;
    COFD[26] = -4.86046736E-01;
    COFD[27] = 2.03659188E-02;
    COFD[28] = -1.52554761E+01;
    COFD[29] = 3.35922578E+00;
    COFD[30] = -2.25181399E-01;
    COFD[31] = 9.92132878E-03;
    COFD[32] = -1.52486273E+01;
    COFD[33] = 3.35922578E+00;
    COFD[34] = -2.25181399E-01;
    COFD[35] = 9.92132878E-03;
    COFD[36] = -1.59404882E+01;
    COFD[37] = 3.66853818E+00;
    COFD[38] = -2.64346221E-01;
    COFD[39] = 1.15784613E-02;
    COFD[40] = -2.04833713E+01;
    COFD[41] = 5.23112374E+00;
    COFD[42] = -4.54967682E-01;
    COFD[43] = 1.93570423E-02;
    COFD[44] = -1.59633387E+01;
    COFD[45] = 3.66853818E+00;
    COFD[46] = -2.64346221E-01;
    COFD[47] = 1.15784613E-02;
    COFD[48] = -1.50031687E+01;
    COFD[49] = 3.26223357E+00;
    COFD[50] = -2.12746642E-01;
    COFD[51] = 9.38912883E-03;
    COFD[52] = -2.02268902E+01;
    COFD[53] = 5.13632093E+00;
    COFD[54] = -4.44839124E-01;
    COFD[55] = 1.90058354E-02;
    COFD[56] = -1.82590824E+01;
    COFD[57] = 4.39538102E+00;
    COFD[58] = -3.56367230E-01;
    COFD[59] = 1.54788461E-02;
    COFD[60] = -1.59327297E+01;
    COFD[61] = 3.65620899E+00;
    COFD[62] = -2.62933804E-01;
    COFD[63] = 1.15253223E-02;
    COFD[64] = -1.78815889E+01;
    COFD[65] = 4.34347890E+00;
    COFD[66] = -3.49890003E-01;
    COFD[67] = 1.52083459E-02;
    COFD[68] = -1.82673770E+01;
    COFD[69] = 4.39538102E+00;
    COFD[70] = -3.56367230E-01;
    COFD[71] = 1.54788461E-02;
    COFD[72] = -1.81432461E+01;
    COFD[73] = 4.37565431E+00;
    COFD[74] = -3.53906025E-01;
    COFD[75] = 1.53760786E-02;
    COFD[76] = -2.04750581E+01;
    COFD[77] = 5.23112374E+00;
    COFD[78] = -4.54967682E-01;
    COFD[79] = 1.93570423E-02;
    COFD[80] = -2.04649069E+01;
    COFD[81] = 5.18856872E+00;
    COFD[82] = -4.50001829E-01;
    COFD[83] = 1.91636142E-02;
    COFD[84] = -2.04688382E+01;
    COFD[85] = 5.18856872E+00;
    COFD[86] = -4.50001829E-01;
    COFD[87] = 1.91636142E-02;
    COFD[88] = -1.83039618E+01;
    COFD[89] = 4.47952077E+00;
    COFD[90] = -3.66569471E-01;
    COFD[91] = 1.58916129E-02;
    COFD[92] = -1.59884305E+01;
    COFD[93] = 3.72220402E+00;
    COFD[94] = -2.71150591E-01;
    COFD[95] = 1.18665265E-02;
    COFD[96] = -1.83137139E+01;
    COFD[97] = 4.47952077E+00;
    COFD[98] = -3.66569471E-01;
    COFD[99] = 1.58916129E-02;
    COFD[100] = -2.02693653E+01;
    COFD[101] = 5.10426133E+00;
    COFD[102] = -4.41256919E-01;
    COFD[103] = 1.88737290E-02;
    COFD[104] = -1.91796663E+01;
    COFD[105] = 4.70714822E+00;
    COFD[106] = -3.94261134E-01;
    COFD[107] = 1.70175169E-02;
    COFD[108] = -1.90859283E+01;
    COFD[109] = 4.68079396E+00;
    COFD[110] = -3.91231550E-01;
    COFD[111] = 1.69021170E-02;
    COFD[112] = -2.03249567E+01;
    COFD[113] = 5.03292431E+00;
    COFD[114] = -4.32887308E-01;
    COFD[115] = 1.85458102E-02;
    COFD[116] = -1.92062897E+01;
    COFD[117] = 4.66318669E+00;
    COFD[118] = -3.89108667E-01;
    COFD[119] = 1.68165377E-02;
    COFD[120] = -1.92062897E+01;
    COFD[121] = 4.66318669E+00;
    COFD[122] = -3.89108667E-01;
    COFD[123] = 1.68165377E-02;
    COFD[124] = -1.92108129E+01;
    COFD[125] = 4.66318669E+00;
    COFD[126] = -3.89108667E-01;
    COFD[127] = 1.68165377E-02;
    COFD[128] = -2.08314053E+01;
    COFD[129] = 5.16623522E+00;
    COFD[130] = -4.47455294E-01;
    COFD[131] = 1.90672494E-02;
    COFD[132] = -2.09943481E+01;
    COFD[133] = 5.22468467E+00;
    COFD[134] = -4.54220128E-01;
    COFD[135] = 1.93281042E-02;
    COFD[136] = -1.40756935E+01;
    COFD[137] = 3.07549274E+00;
    COFD[138] = -1.88889344E-01;
    COFD[139] = 8.37152866E-03;
    COFD[140] = -1.32093628E+01;
    COFD[141] = 2.90778936E+00;
    COFD[142] = -1.67388544E-01;
    COFD[143] = 7.45220609E-03;
    COFD[144] = -1.34230272E+01;
    COFD[145] = 3.48624238E+00;
    COFD[146] = -2.41554467E-01;
    COFD[147] = 1.06263545E-02;
    COFD[148] = -1.32244035E+01;
    COFD[149] = 2.90778936E+00;
    COFD[150] = -1.67388544E-01;
    COFD[151] = 7.45220609E-03;
    COFD[152] = -1.09595712E+01;
    COFD[153] = 2.30836460E+00;
    COFD[154] = -8.76339315E-02;
    COFD[155] = 3.90878445E-03;
    COFD[156] = -1.43139231E+01;
    COFD[157] = 3.17651319E+00;
    COFD[158] = -2.02028974E-01;
    COFD[159] = 8.94232502E-03;
    COFD[160] = -1.94093572E+01;
    COFD[161] = 5.16013126E+00;
    COFD[162] = -4.46824543E-01;
    COFD[163] = 1.90464887E-02;
    COFD[164] = -1.43238998E+01;
    COFD[165] = 3.17651319E+00;
    COFD[166] = -2.02028974E-01;
    COFD[167] = 8.94232502E-03;
    COFD[168] = -1.43190389E+01;
    COFD[169] = 3.17651319E+00;
    COFD[170] = -2.02028974E-01;
    COFD[171] = 8.94232502E-03;
    COFD[172] = -1.50584249E+01;
    COFD[173] = 3.47945612E+00;
    COFD[174] = -2.40703722E-01;
    COFD[175] = 1.05907441E-02;
    COFD[176] = -1.94373127E+01;
    COFD[177] = 5.02567894E+00;
    COFD[178] = -4.32045169E-01;
    COFD[179] = 1.85132214E-02;
    COFD[180] = -1.50766130E+01;
    COFD[181] = 3.47945612E+00;
    COFD[182] = -2.40703722E-01;
    COFD[183] = 1.05907441E-02;
    COFD[184] = -1.40999008E+01;
    COFD[185] = 3.08120012E+00;
    COFD[186] = -1.89629903E-01;
    COFD[187] = 8.40361952E-03;
    COFD[188] = -1.88179418E+01;
    COFD[189] = 4.79683898E+00;
    COFD[190] = -4.04829719E-01;
    COFD[191] = 1.74325475E-02;
    COFD[192] = -1.72053106E+01;
    COFD[193] = 4.15807461E+00;
    COFD[194] = -3.27178539E-01;
    COFD[195] = 1.42784349E-02;
    COFD[196] = -1.50270339E+01;
    COFD[197] = 3.46140064E+00;
    COFD[198] = -2.38440092E-01;
    COFD[199] = 1.04960087E-02;
    COFD[200] = -1.68343393E+01;
    COFD[201] = 4.11954900E+00;
    COFD[202] = -3.22470391E-01;
    COFD[203] = 1.40859564E-02;
    COFD[204] = -1.72112971E+01;
    COFD[205] = 4.15807461E+00;
    COFD[206] = -3.27178539E-01;
    COFD[207] = 1.42784349E-02;
    COFD[208] = -1.70534856E+01;
    COFD[209] = 4.14240922E+00;
    COFD[210] = -3.25239774E-01;
    COFD[211] = 1.41980687E-02;
    COFD[212] = -1.94313116E+01;
    COFD[213] = 5.02567894E+00;
    COFD[214] = -4.32045169E-01;
    COFD[215] = 1.85132214E-02;
    COFD[216] = -1.93925667E+01;
    COFD[217] = 4.98286777E+00;
    COFD[218] = -4.26970814E-01;
    COFD[219] = 1.83122917E-02;
    COFD[220] = -1.93952366E+01;
    COFD[221] = 4.98286777E+00;
    COFD[222] = -4.26970814E-01;
    COFD[223] = 1.83122917E-02;
    COFD[224] = -1.72286007E+01;
    COFD[225] = 4.24084025E+00;
    COFD[226] = -3.37428619E-01;
    COFD[227] = 1.47032793E-02;
    COFD[228] = -1.49500357E+01;
    COFD[229] = 3.52327209E+00;
    COFD[230] = -2.46286208E-01;
    COFD[231] = 1.08285963E-02;
    COFD[232] = -1.72357436E+01;
    COFD[233] = 4.24084025E+00;
    COFD[234] = -3.37428619E-01;
    COFD[235] = 1.47032793E-02;
    COFD[236] = -1.90915649E+01;
    COFD[237] = 4.84384483E+00;
    COFD[238] = -4.10265575E-01;
    COFD[239] = 1.76414287E-02;
    COFD[240] = -1.80480958E+01;
    COFD[241] = 4.45434023E+00;
    COFD[242] = -3.63584633E-01;
    COFD[243] = 1.57739270E-02;
    COFD[244] = -1.79361160E+01;
    COFD[245] = 4.42139452E+00;
    COFD[246] = -3.59567329E-01;
    COFD[247] = 1.56103969E-02;
    COFD[248] = -1.91805416E+01;
    COFD[249] = 4.78152337E+00;
    COFD[250] = -4.03052726E-01;
    COFD[251] = 1.73639773E-02;
    COFD[252] = -1.80724788E+01;
    COFD[253] = 4.40247898E+00;
    COFD[254] = -3.57238362E-01;
    COFD[255] = 1.55145651E-02;
    COFD[256] = -1.80724788E+01;
    COFD[257] = 4.40247898E+00;
    COFD[258] = -3.57238362E-01;
    COFD[259] = 1.55145651E-02;
    COFD[260] = -1.80755831E+01;
    COFD[261] = 4.40247898E+00;
    COFD[262] = -3.57238362E-01;
    COFD[263] = 1.55145651E-02;
    COFD[264] = -1.96652729E+01;
    COFD[265] = 4.92149009E+00;
    COFD[266] = -4.19702348E-01;
    COFD[267] = 1.80250370E-02;
    COFD[268] = -1.98296243E+01;
    COFD[269] = 4.98207523E+00;
    COFD[270] = -4.26877291E-01;
    COFD[271] = 1.83086094E-02;
    COFD[272] = -1.42894441E+01;
    COFD[273] = 3.67490723E+00;
    COFD[274] = -2.65114792E-01;
    COFD[275] = 1.16092671E-02;
    COFD[276] = -1.34230272E+01;
    COFD[277] = 3.48624238E+00;
    COFD[278] = -2.41554467E-01;
    COFD[279] = 1.06263545E-02;
    COFD[280] = -1.47968712E+01;
    COFD[281] = 4.23027636E+00;
    COFD[282] = -3.36139991E-01;
    COFD[283] = 1.46507621E-02;
    COFD[284] = -1.34247866E+01;
    COFD[285] = 3.48624238E+00;
    COFD[286] = -2.41554467E-01;
    COFD[287] = 1.06263545E-02;
    COFD[288] = -1.14366381E+01;
    COFD[289] = 2.78323501E+00;
    COFD[290] = -1.51214064E-01;
    COFD[291] = 6.75150012E-03;
    COFD[292] = -1.46550083E+01;
    COFD[293] = 3.83606243E+00;
    COFD[294] = -2.86076532E-01;
    COFD[295] = 1.25205829E-02;
    COFD[296] = -1.95739570E+01;
    COFD[297] = 5.61113230E+00;
    COFD[298] = -4.90190187E-01;
    COFD[299] = 2.03260675E-02;
    COFD[300] = -1.46559141E+01;
    COFD[301] = 3.83606243E+00;
    COFD[302] = -2.86076532E-01;
    COFD[303] = 1.25205829E-02;
    COFD[304] = -1.46554748E+01;
    COFD[305] = 3.83606243E+00;
    COFD[306] = -2.86076532E-01;
    COFD[307] = 1.25205829E-02;
    COFD[308] = -1.57972369E+01;
    COFD[309] = 4.22225052E+00;
    COFD[310] = -3.35156428E-01;
    COFD[311] = 1.46104855E-02;
    COFD[312] = -1.97550088E+01;
    COFD[313] = 5.56931926E+00;
    COFD[314] = -4.89105511E-01;
    COFD[315] = 2.04493129E-02;
    COFD[316] = -1.57994893E+01;
    COFD[317] = 4.22225052E+00;
    COFD[318] = -3.35156428E-01;
    COFD[319] = 1.46104855E-02;
    COFD[320] = -1.43151174E+01;
    COFD[321] = 3.68038508E+00;
    COFD[322] = -2.65779346E-01;
    COFD[323] = 1.16360771E-02;
    COFD[324] = -1.92718582E+01;
    COFD[325] = 5.41172124E+00;
    COFD[326] = -4.73213887E-01;
    COFD[327] = 1.99405473E-02;
    COFD[328] = -1.78631557E+01;
    COFD[329] = 4.88268692E+00;
    COFD[330] = -4.14917638E-01;
    COFD[331] = 1.78274298E-02;
    COFD[332] = -1.57199037E+01;
    COFD[333] = 4.19936335E+00;
    COFD[334] = -3.32311009E-01;
    COFD[335] = 1.44921003E-02;
    COFD[336] = -1.74407963E+01;
    COFD[337] = 4.83580036E+00;
    COFD[338] = -4.09383573E-01;
    COFD[339] = 1.76098175E-02;
    COFD[340] = -1.78637178E+01;
    COFD[341] = 4.88268692E+00;
    COFD[342] = -4.14917638E-01;
    COFD[343] = 1.78274298E-02;
    COFD[344] = -1.76147026E+01;
    COFD[345] = 4.86049500E+00;
    COFD[346] = -4.12200578E-01;
    COFD[347] = 1.77160971E-02;
    COFD[348] = -1.97544450E+01;
    COFD[349] = 5.56931926E+00;
    COFD[350] = -4.89105511E-01;
    COFD[351] = 2.04493129E-02;
    COFD[352] = -1.96914944E+01;
    COFD[353] = 5.54637286E+00;
    COFD[354] = -4.87070324E-01;
    COFD[355] = 2.03983467E-02;
    COFD[356] = -1.96917146E+01;
    COFD[357] = 5.54637286E+00;
    COFD[358] = -4.87070324E-01;
    COFD[359] = 2.03983467E-02;
    COFD[360] = -1.79310765E+01;
    COFD[361] = 4.98037650E+00;
    COFD[362] = -4.26676911E-01;
    COFD[363] = 1.83007231E-02;
    COFD[364] = -1.54460820E+01;
    COFD[365] = 4.26819983E+00;
    COFD[366] = -3.40766379E-01;
    COFD[367] = 1.48393361E-02;
    COFD[368] = -1.79317714E+01;
    COFD[369] = 4.98037650E+00;
    COFD[370] = -4.26676911E-01;
    COFD[371] = 1.83007231E-02;
    COFD[372] = -1.94691430E+01;
    COFD[373] = 5.43830787E+00;
    COFD[374] = -4.75472880E-01;
    COFD[375] = 1.99909996E-02;
    COFD[376] = -1.86493112E+01;
    COFD[377] = 5.16040659E+00;
    COFD[378] = -4.46843492E-01;
    COFD[379] = 1.90466181E-02;
    COFD[380] = -1.85748546E+01;
    COFD[381] = 5.14789919E+00;
    COFD[382] = -4.45930850E-01;
    COFD[383] = 1.90363341E-02;
    COFD[384] = -1.96115230E+01;
    COFD[385] = 5.40502814E+00;
    COFD[386] = -4.72746583E-01;
    COFD[387] = 1.99363615E-02;
    COFD[388] = -1.87481780E+01;
    COFD[389] = 5.13858656E+00;
    COFD[390] = -4.45075387E-01;
    COFD[391] = 1.90137309E-02;
    COFD[392] = -1.87481780E+01;
    COFD[393] = 5.13858656E+00;
    COFD[394] = -4.45075387E-01;
    COFD[395] = 1.90137309E-02;
    COFD[396] = -1.87484393E+01;
    COFD[397] = 5.13858656E+00;
    COFD[398] = -4.45075387E-01;
    COFD[399] = 1.90137309E-02;
    COFD[400] = -1.99919894E+01;
    COFD[401] = 5.50111590E+00;
    COFD[402] = -4.82434754E-01;
    COFD[403] = 2.02461869E-02;
    COFD[404] = -2.01262921E+01;
    COFD[405] = 5.54581286E+00;
    COFD[406] = -4.87014004E-01;
    COFD[407] = 2.03965482E-02;
    COFD[408] = -1.40949196E+01;
    COFD[409] = 3.07549274E+00;
    COFD[410] = -1.88889344E-01;
    COFD[411] = 8.37152866E-03;
    COFD[412] = -1.32244035E+01;
    COFD[413] = 2.90778936E+00;
    COFD[414] = -1.67388544E-01;
    COFD[415] = 7.45220609E-03;
    COFD[416] = -1.34247866E+01;
    COFD[417] = 3.48624238E+00;
    COFD[418] = -2.41554467E-01;
    COFD[419] = 1.06263545E-02;
    COFD[420] = -1.32399106E+01;
    COFD[421] = 2.90778936E+00;
    COFD[422] = -1.67388544E-01;
    COFD[423] = 7.45220609E-03;
    COFD[424] = -1.09628982E+01;
    COFD[425] = 2.30836460E+00;
    COFD[426] = -8.76339315E-02;
    COFD[427] = 3.90878445E-03;
    COFD[428] = -1.43340796E+01;
    COFD[429] = 3.17651319E+00;
    COFD[430] = -2.02028974E-01;
    COFD[431] = 8.94232502E-03;
    COFD[432] = -1.94253036E+01;
    COFD[433] = 5.16013126E+00;
    COFD[434] = -4.46824543E-01;
    COFD[435] = 1.90464887E-02;
    COFD[436] = -1.43444709E+01;
    COFD[437] = 3.17651319E+00;
    COFD[438] = -2.02028974E-01;
    COFD[439] = 8.94232502E-03;
    COFD[440] = -1.43394069E+01;
    COFD[441] = 3.17651319E+00;
    COFD[442] = -2.02028974E-01;
    COFD[443] = 8.94232502E-03;
    COFD[444] = -1.50724636E+01;
    COFD[445] = 3.47945612E+00;
    COFD[446] = -2.40703722E-01;
    COFD[447] = 1.05907441E-02;
    COFD[448] = -1.94570287E+01;
    COFD[449] = 5.02567894E+00;
    COFD[450] = -4.32045169E-01;
    COFD[451] = 1.85132214E-02;
    COFD[452] = -1.50911794E+01;
    COFD[453] = 3.47945612E+00;
    COFD[454] = -2.40703722E-01;
    COFD[455] = 1.05907441E-02;
    COFD[456] = -1.41191261E+01;
    COFD[457] = 3.08120012E+00;
    COFD[458] = -1.89629903E-01;
    COFD[459] = 8.40361952E-03;
    COFD[460] = -1.88378874E+01;
    COFD[461] = 4.79683898E+00;
    COFD[462] = -4.04829719E-01;
    COFD[463] = 1.74325475E-02;
    COFD[464] = -1.72247972E+01;
    COFD[465] = 4.15807461E+00;
    COFD[466] = -3.27178539E-01;
    COFD[467] = 1.42784349E-02;
    COFD[468] = -1.50420953E+01;
    COFD[469] = 3.46140064E+00;
    COFD[470] = -2.38440092E-01;
    COFD[471] = 1.04960087E-02;
    COFD[472] = -1.68535757E+01;
    COFD[473] = 4.11954900E+00;
    COFD[474] = -3.22470391E-01;
    COFD[475] = 1.40859564E-02;
    COFD[476] = -1.72310232E+01;
    COFD[477] = 4.15807461E+00;
    COFD[478] = -3.27178539E-01;
    COFD[479] = 1.42784349E-02;
    COFD[480] = -1.70757047E+01;
    COFD[481] = 4.14240922E+00;
    COFD[482] = -3.25239774E-01;
    COFD[483] = 1.41980687E-02;
    COFD[484] = -1.94507876E+01;
    COFD[485] = 5.02567894E+00;
    COFD[486] = -4.32045169E-01;
    COFD[487] = 1.85132214E-02;
    COFD[488] = -1.94151822E+01;
    COFD[489] = 4.98286777E+00;
    COFD[490] = -4.26970814E-01;
    COFD[491] = 1.83122917E-02;
    COFD[492] = -1.94179760E+01;
    COFD[493] = 4.98286777E+00;
    COFD[494] = -4.26970814E-01;
    COFD[495] = 1.83122917E-02;
    COFD[496] = -1.72473011E+01;
    COFD[497] = 4.24084025E+00;
    COFD[498] = -3.37428619E-01;
    COFD[499] = 1.47032793E-02;
    COFD[500] = -1.49718233E+01;
    COFD[501] = 3.52327209E+00;
    COFD[502] = -2.46286208E-01;
    COFD[503] = 1.08285963E-02;
    COFD[504] = -1.72547182E+01;
    COFD[505] = 4.24084025E+00;
    COFD[506] = -3.37428619E-01;
    COFD[507] = 1.47032793E-02;
    COFD[508] = -1.91136491E+01;
    COFD[509] = 4.84384483E+00;
    COFD[510] = -4.10265575E-01;
    COFD[511] = 1.76414287E-02;
    COFD[512] = -1.80698901E+01;
    COFD[513] = 4.45434023E+00;
    COFD[514] = -3.63584633E-01;
    COFD[515] = 1.57739270E-02;
    COFD[516] = -1.79580609E+01;
    COFD[517] = 4.42139452E+00;
    COFD[518] = -3.59567329E-01;
    COFD[519] = 1.56103969E-02;
    COFD[520] = -1.92042395E+01;
    COFD[521] = 4.78152337E+00;
    COFD[522] = -4.03052726E-01;
    COFD[523] = 1.73639773E-02;
    COFD[524] = -1.80945693E+01;
    COFD[525] = 4.40247898E+00;
    COFD[526] = -3.57238362E-01;
    COFD[527] = 1.55145651E-02;
    COFD[528] = -1.80945693E+01;
    COFD[529] = 4.40247898E+00;
    COFD[530] = -3.57238362E-01;
    COFD[531] = 1.55145651E-02;
    COFD[532] = -1.80978142E+01;
    COFD[533] = 4.40247898E+00;
    COFD[534] = -3.57238362E-01;
    COFD[535] = 1.55145651E-02;
    COFD[536] = -1.96903181E+01;
    COFD[537] = 4.92149009E+00;
    COFD[538] = -4.19702348E-01;
    COFD[539] = 1.80250370E-02;
    COFD[540] = -1.98546695E+01;
    COFD[541] = 4.98207523E+00;
    COFD[542] = -4.26877291E-01;
    COFD[543] = 1.83086094E-02;
    COFD[544] = -1.16906297E+01;
    COFD[545] = 2.47469981E+00;
    COFD[546] = -1.10436257E-01;
    COFD[547] = 4.95273813E-03;
    COFD[548] = -1.09595712E+01;
    COFD[549] = 2.30836460E+00;
    COFD[550] = -8.76339315E-02;
    COFD[551] = 3.90878445E-03;
    COFD[552] = -1.14366381E+01;
    COFD[553] = 2.78323501E+00;
    COFD[554] = -1.51214064E-01;
    COFD[555] = 6.75150012E-03;
    COFD[556] = -1.09628982E+01;
    COFD[557] = 2.30836460E+00;
    COFD[558] = -8.76339315E-02;
    COFD[559] = 3.90878445E-03;
    COFD[560] = -1.03270606E+01;
    COFD[561] = 2.19285409E+00;
    COFD[562] = -7.54492786E-02;
    COFD[563] = 3.51398213E-03;
    COFD[564] = -1.18988955E+01;
    COFD[565] = 2.57507000E+00;
    COFD[566] = -1.24033737E-01;
    COFD[567] = 5.56694959E-03;
    COFD[568] = -1.71982995E+01;
    COFD[569] = 4.63881404E+00;
    COFD[570] = -3.86139633E-01;
    COFD[571] = 1.66955081E-02;
    COFD[572] = -1.19006548E+01;
    COFD[573] = 2.57507000E+00;
    COFD[574] = -1.24033737E-01;
    COFD[575] = 5.56694959E-03;
    COFD[576] = -1.18998012E+01;
    COFD[577] = 2.57507000E+00;
    COFD[578] = -1.24033737E-01;
    COFD[579] = 5.56694959E-03;
    COFD[580] = -1.25098960E+01;
    COFD[581] = 2.77873601E+00;
    COFD[582] = -1.50637360E-01;
    COFD[583] = 6.72684281E-03;
    COFD[584] = -1.60528285E+01;
    COFD[585] = 4.11188603E+00;
    COFD[586] = -3.21540884E-01;
    COFD[587] = 1.40482564E-02;
    COFD[588] = -1.25141260E+01;
    COFD[589] = 2.77873601E+00;
    COFD[590] = -1.50637360E-01;
    COFD[591] = 6.72684281E-03;
    COFD[592] = -1.17159737E+01;
    COFD[593] = 2.48123210E+00;
    COFD[594] = -1.11322604E-01;
    COFD[595] = 4.99282389E-03;
    COFD[596] = -1.58456300E+01;
    COFD[597] = 4.02074783E+00;
    COFD[598] = -3.10018522E-01;
    COFD[599] = 1.35599552E-02;
    COFD[600] = -1.39648112E+01;
    COFD[601] = 3.24966086E+00;
    COFD[602] = -2.11199992E-01;
    COFD[603] = 9.32580661E-03;
    COFD[604] = -1.24693568E+01;
    COFD[605] = 2.76686648E+00;
    COFD[606] = -1.49120141E-01;
    COFD[607] = 6.66220432E-03;
    COFD[608] = -1.36336373E+01;
    COFD[609] = 3.22088176E+00;
    COFD[610] = -2.07623790E-01;
    COFD[611] = 9.17771542E-03;
    COFD[612] = -1.39658996E+01;
    COFD[613] = 3.24966086E+00;
    COFD[614] = -2.11199992E-01;
    COFD[615] = 9.32580661E-03;
    COFD[616] = -1.37794315E+01;
    COFD[617] = 3.23973858E+00;
    COFD[618] = -2.09989036E-01;
    COFD[619] = 9.27667906E-03;
    COFD[620] = -1.60517370E+01;
    COFD[621] = 4.11188603E+00;
    COFD[622] = -3.21540884E-01;
    COFD[623] = 1.40482564E-02;
    COFD[624] = -1.59632479E+01;
    COFD[625] = 4.07051484E+00;
    COFD[626] = -3.16303109E-01;
    COFD[627] = 1.38259377E-02;
    COFD[628] = -1.59636793E+01;
    COFD[629] = 4.07051484E+00;
    COFD[630] = -3.16303109E-01;
    COFD[631] = 1.38259377E-02;
    COFD[632] = -1.39315266E+01;
    COFD[633] = 3.30394764E+00;
    COFD[634] = -2.17920112E-01;
    COFD[635] = 9.60284243E-03;
    COFD[636] = -1.22004324E+01;
    COFD[637] = 2.80725489E+00;
    COFD[638] = -1.54291406E-01;
    COFD[639] = 6.88290911E-03;
    COFD[640] = -1.39328674E+01;
    COFD[641] = 3.30394764E+00;
    COFD[642] = -2.17920112E-01;
    COFD[643] = 9.60284243E-03;
    COFD[644] = -1.57040212E+01;
    COFD[645] = 3.93614244E+00;
    COFD[646] = -2.99111497E-01;
    COFD[647] = 1.30888229E-02;
    COFD[648] = -1.46719197E+01;
    COFD[649] = 3.52400594E+00;
    COFD[650] = -2.46379985E-01;
    COFD[651] = 1.08326032E-02;
    COFD[652] = -1.45715797E+01;
    COFD[653] = 3.49477850E+00;
    COFD[654] = -2.42635772E-01;
    COFD[655] = 1.06721490E-02;
    COFD[656] = -1.56932803E+01;
    COFD[657] = 3.84112355E+00;
    COFD[658] = -2.86739743E-01;
    COFD[659] = 1.25496520E-02;
    COFD[660] = -1.47137939E+01;
    COFD[661] = 3.48023191E+00;
    COFD[662] = -2.40800798E-01;
    COFD[663] = 1.05947990E-02;
    COFD[664] = -1.47137939E+01;
    COFD[665] = 3.48023191E+00;
    COFD[666] = -2.40800798E-01;
    COFD[667] = 1.05947990E-02;
    COFD[668] = -1.47143050E+01;
    COFD[669] = 3.48023191E+00;
    COFD[670] = -2.40800798E-01;
    COFD[671] = 1.05947990E-02;
    COFD[672] = -1.63371298E+01;
    COFD[673] = 4.06047384E+00;
    COFD[674] = -3.15033719E-01;
    COFD[675] = 1.37721355E-02;
    COFD[676] = -1.64819183E+01;
    COFD[677] = 4.11726215E+00;
    COFD[678] = -3.22193015E-01;
    COFD[679] = 1.40747074E-02;
    COFD[680] = -1.52414485E+01;
    COFD[681] = 3.35922578E+00;
    COFD[682] = -2.25181399E-01;
    COFD[683] = 9.92132878E-03;
    COFD[684] = -1.43139231E+01;
    COFD[685] = 3.17651319E+00;
    COFD[686] = -2.02028974E-01;
    COFD[687] = 8.94232502E-03;
    COFD[688] = -1.46550083E+01;
    COFD[689] = 3.83606243E+00;
    COFD[690] = -2.86076532E-01;
    COFD[691] = 1.25205829E-02;
    COFD[692] = -1.43340796E+01;
    COFD[693] = 3.17651319E+00;
    COFD[694] = -2.02028974E-01;
    COFD[695] = 8.94232502E-03;
    COFD[696] = -1.18988955E+01;
    COFD[697] = 2.57507000E+00;
    COFD[698] = -1.24033737E-01;
    COFD[699] = 5.56694959E-03;
    COFD[700] = -1.55511344E+01;
    COFD[701] = 3.48070094E+00;
    COFD[702] = -2.40859499E-01;
    COFD[703] = 1.05972514E-02;
    COFD[704] = -2.12652533E+01;
    COFD[705] = 5.59961818E+00;
    COFD[706] = -4.91624858E-01;
    COFD[707] = 2.05035550E-02;
    COFD[708] = -1.55661750E+01;
    COFD[709] = 3.48070094E+00;
    COFD[710] = -2.40859499E-01;
    COFD[711] = 1.05972514E-02;
    COFD[712] = -1.55588279E+01;
    COFD[713] = 3.48070094E+00;
    COFD[714] = -2.40859499E-01;
    COFD[715] = 1.05972514E-02;
    COFD[716] = -1.63254691E+01;
    COFD[717] = 3.82388595E+00;
    COFD[718] = -2.84480724E-01;
    COFD[719] = 1.24506311E-02;
    COFD[720] = -2.08293255E+01;
    COFD[721] = 5.35267674E+00;
    COFD[722] = -4.69010505E-01;
    COFD[723] = 1.98979152E-02;
    COFD[724] = -1.63493345E+01;
    COFD[725] = 3.82388595E+00;
    COFD[726] = -2.84480724E-01;
    COFD[727] = 1.24506311E-02;
    COFD[728] = -1.52721107E+01;
    COFD[729] = 3.36790500E+00;
    COFD[730] = -2.26321740E-01;
    COFD[731] = 9.97135055E-03;
    COFD[732] = -2.04928958E+01;
    COFD[733] = 5.22397933E+00;
    COFD[734] = -4.54138171E-01;
    COFD[735] = 1.93249285E-02;
    COFD[736] = -1.85756076E+01;
    COFD[737] = 4.51052425E+00;
    COFD[738] = -3.70301627E-01;
    COFD[739] = 1.60416153E-02;
    COFD[740] = -1.62724462E+01;
    COFD[741] = 3.79163564E+00;
    COFD[742] = -2.80257365E-01;
    COFD[743] = 1.22656902E-02;
    COFD[744] = -1.82145353E+01;
    COFD[745] = 4.46848269E+00;
    COFD[746] = -3.65269718E-01;
    COFD[747] = 1.58407652E-02;
    COFD[748] = -1.85844688E+01;
    COFD[749] = 4.51052425E+00;
    COFD[750] = -3.70301627E-01;
    COFD[751] = 1.60416153E-02;
    COFD[752] = -1.84688406E+01;
    COFD[753] = 4.49330851E+00;
    COFD[754] = -3.68208715E-01;
    COFD[755] = 1.59565402E-02;
    COFD[756] = -2.08204449E+01;
    COFD[757] = 5.35267674E+00;
    COFD[758] = -4.69010505E-01;
    COFD[759] = 1.98979152E-02;
    COFD[760] = -2.08463209E+01;
    COFD[761] = 5.32244593E+00;
    COFD[762] = -4.65829403E-01;
    COFD[763] = 1.97895274E-02;
    COFD[764] = -2.08505864E+01;
    COFD[765] = 5.32244593E+00;
    COFD[766] = -4.65829403E-01;
    COFD[767] = 1.97895274E-02;
    COFD[768] = -1.86507213E+01;
    COFD[769] = 4.60874797E+00;
    COFD[770] = -3.82368716E-01;
    COFD[771] = 1.65370164E-02;
    COFD[772] = -1.64169433E+01;
    COFD[773] = 3.89309916E+00;
    COFD[774] = -2.93528188E-01;
    COFD[775] = 1.28463177E-02;
    COFD[776] = -1.86611023E+01;
    COFD[777] = 4.60874797E+00;
    COFD[778] = -3.82368716E-01;
    COFD[779] = 1.65370164E-02;
    COFD[780] = -2.05235731E+01;
    COFD[781] = 5.18417470E+00;
    COFD[782] = -4.49491573E-01;
    COFD[783] = 1.91438508E-02;
    COFD[784] = -1.94912151E+01;
    COFD[785] = 4.81575071E+00;
    COFD[786] = -4.07042139E-01;
    COFD[787] = 1.75187504E-02;
    COFD[788] = -1.93917298E+01;
    COFD[789] = 4.78708023E+00;
    COFD[790] = -4.03693144E-01;
    COFD[791] = 1.73884817E-02;
    COFD[792] = -2.06513321E+01;
    COFD[793] = 5.14159226E+00;
    COFD[794] = -4.45395793E-01;
    COFD[795] = 1.90248052E-02;
    COFD[796] = -1.95201830E+01;
    COFD[797] = 4.77151544E+00;
    COFD[798] = -4.01882811E-01;
    COFD[799] = 1.73184814E-02;
    COFD[800] = -1.95201830E+01;
    COFD[801] = 4.77151544E+00;
    COFD[802] = -4.01882811E-01;
    COFD[803] = 1.73184814E-02;
    COFD[804] = -1.95250774E+01;
    COFD[805] = 4.77151544E+00;
    COFD[806] = -4.01882811E-01;
    COFD[807] = 1.73184814E-02;
    COFD[808] = -2.12147591E+01;
    COFD[809] = 5.29477845E+00;
    COFD[810] = -4.62516356E-01;
    COFD[811] = 1.96564008E-02;
    COFD[812] = -2.13698722E+01;
    COFD[813] = 5.34971865E+00;
    COFD[814] = -4.68771123E-01;
    COFD[815] = 1.98933811E-02;
    COFD[816] = -2.10643259E+01;
    COFD[817] = 5.53614847E+00;
    COFD[818] = -4.86046736E-01;
    COFD[819] = 2.03659188E-02;
    COFD[820] = -1.94093572E+01;
    COFD[821] = 5.16013126E+00;
    COFD[822] = -4.46824543E-01;
    COFD[823] = 1.90464887E-02;
    COFD[824] = -1.95739570E+01;
    COFD[825] = 5.61113230E+00;
    COFD[826] = -4.90190187E-01;
    COFD[827] = 2.03260675E-02;
    COFD[828] = -1.94253036E+01;
    COFD[829] = 5.16013126E+00;
    COFD[830] = -4.46824543E-01;
    COFD[831] = 1.90464887E-02;
    COFD[832] = -1.71982995E+01;
    COFD[833] = 4.63881404E+00;
    COFD[834] = -3.86139633E-01;
    COFD[835] = 1.66955081E-02;
    COFD[836] = -2.12652533E+01;
    COFD[837] = 5.59961818E+00;
    COFD[838] = -4.91624858E-01;
    COFD[839] = 2.05035550E-02;
    COFD[840] = -1.19157919E+01;
    COFD[841] = 9.28955130E-01;
    COFD[842] = 2.42107090E-01;
    COFD[843] = -1.59823963E-02;
    COFD[844] = -2.06516336E+01;
    COFD[845] = 5.41688482E+00;
    COFD[846] = -4.73387188E-01;
    COFD[847] = 1.99280175E-02;
    COFD[848] = -2.06463744E+01;
    COFD[849] = 5.41688482E+00;
    COFD[850] = -4.73387188E-01;
    COFD[851] = 1.99280175E-02;
    COFD[852] = -2.12639214E+01;
    COFD[853] = 5.61184117E+00;
    COFD[854] = -4.90532156E-01;
    COFD[855] = 2.03507922E-02;
    COFD[856] = -1.77563250E+01;
    COFD[857] = 3.57475686E+00;
    COFD[858] = -1.56396297E-01;
    COFD[859] = 3.12157721E-03;
    COFD[860] = -2.12831323E+01;
    COFD[861] = 5.61184117E+00;
    COFD[862] = -4.90532156E-01;
    COFD[863] = 2.03507922E-02;
    COFD[864] = -2.11388331E+01;
    COFD[865] = 5.55529675E+00;
    COFD[866] = -4.87942518E-01;
    COFD[867] = 2.04249054E-02;
    COFD[868] = -1.65295288E+01;
    COFD[869] = 2.97569206E+00;
    COFD[870] = -6.75652842E-02;
    COFD[871] = -1.08648422E-03;
    COFD[872] = -2.13084334E+01;
    COFD[873] = 5.27210469E+00;
    COFD[874] = -4.21419216E-01;
    COFD[875] = 1.63567178E-02;
    COFD[876] = -2.14087397E+01;
    COFD[877] = 5.57282008E+00;
    COFD[878] = -4.76690890E-01;
    COFD[879] = 1.94000719E-02;
    COFD[880] = -2.11309197E+01;
    COFD[881] = 5.32644193E+00;
    COFD[882] = -4.30581064E-01;
    COFD[883] = 1.68379725E-02;
    COFD[884] = -2.13148887E+01;
    COFD[885] = 5.27210469E+00;
    COFD[886] = -4.21419216E-01;
    COFD[887] = 1.63567178E-02;
    COFD[888] = -2.07653719E+01;
    COFD[889] = 5.01092022E+00;
    COFD[890] = -3.77985635E-01;
    COFD[891] = 1.40968645E-02;
    COFD[892] = -1.77498543E+01;
    COFD[893] = 3.57475686E+00;
    COFD[894] = -1.56396297E-01;
    COFD[895] = 3.12157721E-03;
    COFD[896] = -1.80862867E+01;
    COFD[897] = 3.69199168E+00;
    COFD[898] = -1.74005516E-01;
    COFD[899] = 3.97694372E-03;
    COFD[900] = -1.80892005E+01;
    COFD[901] = 3.69199168E+00;
    COFD[902] = -1.74005516E-01;
    COFD[903] = 3.97694372E-03;
    COFD[904] = -2.09565916E+01;
    COFD[905] = 5.18380539E+00;
    COFD[906] = -4.06234719E-01;
    COFD[907] = 1.55515345E-02;
    COFD[908] = -2.10440675E+01;
    COFD[909] = 5.59806282E+00;
    COFD[910] = -4.87109535E-01;
    COFD[911] = 2.01370226E-02;
    COFD[912] = -2.09642705E+01;
    COFD[913] = 5.18380539E+00;
    COFD[914] = -4.06234719E-01;
    COFD[915] = 1.55515345E-02;
    COFD[916] = -1.87419199E+01;
    COFD[917] = 3.96926341E+00;
    COFD[918] = -2.16412264E-01;
    COFD[919] = 6.06012078E-03;
    COFD[920] = -2.05372411E+01;
    COFD[921] = 4.83379373E+00;
    COFD[922] = -3.50008083E-01;
    COFD[923] = 1.26863426E-02;
    COFD[924] = -2.06310304E+01;
    COFD[925] = 4.89289496E+00;
    COFD[926] = -3.59346263E-01;
    COFD[927] = 1.31570901E-02;
    COFD[928] = -1.94510032E+01;
    COFD[929] = 4.17438087E+00;
    COFD[930] = -2.47331309E-01;
    COFD[931] = 7.56793055E-03;
    COFD[932] = -2.08879167E+01;
    COFD[933] = 4.92602269E+00;
    COFD[934] = -3.64572914E-01;
    COFD[935] = 1.34203681E-02;
    COFD[936] = -2.08879167E+01;
    COFD[937] = 4.92602269E+00;
    COFD[938] = -3.64572914E-01;
    COFD[939] = 1.34203681E-02;
    COFD[940] = -2.08912977E+01;
    COFD[941] = 4.92602269E+00;
    COFD[942] = -3.64572914E-01;
    COFD[943] = 1.34203681E-02;
    COFD[944] = -1.76259860E+01;
    COFD[945] = 3.29707321E+00;
    COFD[946] = -1.19066759E-01;
    COFD[947] = 1.47805843E-03;
    COFD[948] = -1.73636900E+01;
    COFD[949] = 3.17377130E+00;
    COFD[950] = -1.00394383E-01;
    COFD[951] = 5.69083899E-04;
    COFD[952] = -1.52554761E+01;
    COFD[953] = 3.35922578E+00;
    COFD[954] = -2.25181399E-01;
    COFD[955] = 9.92132878E-03;
    COFD[956] = -1.43238998E+01;
    COFD[957] = 3.17651319E+00;
    COFD[958] = -2.02028974E-01;
    COFD[959] = 8.94232502E-03;
    COFD[960] = -1.46559141E+01;
    COFD[961] = 3.83606243E+00;
    COFD[962] = -2.86076532E-01;
    COFD[963] = 1.25205829E-02;
    COFD[964] = -1.43444709E+01;
    COFD[965] = 3.17651319E+00;
    COFD[966] = -2.02028974E-01;
    COFD[967] = 8.94232502E-03;
    COFD[968] = -1.19006548E+01;
    COFD[969] = 2.57507000E+00;
    COFD[970] = -1.24033737E-01;
    COFD[971] = 5.56694959E-03;
    COFD[972] = -1.55661750E+01;
    COFD[973] = 3.48070094E+00;
    COFD[974] = -2.40859499E-01;
    COFD[975] = 1.05972514E-02;
    COFD[976] = -2.06516336E+01;
    COFD[977] = 5.41688482E+00;
    COFD[978] = -4.73387188E-01;
    COFD[979] = 1.99280175E-02;
    COFD[980] = -1.55816822E+01;
    COFD[981] = 3.48070094E+00;
    COFD[982] = -2.40859499E-01;
    COFD[983] = 1.05972514E-02;
    COFD[984] = -1.55741053E+01;
    COFD[985] = 3.48070094E+00;
    COFD[986] = -2.40859499E-01;
    COFD[987] = 1.05972514E-02;
    COFD[988] = -1.63345829E+01;
    COFD[989] = 3.82388595E+00;
    COFD[990] = -2.84480724E-01;
    COFD[991] = 1.24506311E-02;
    COFD[992] = -2.08438809E+01;
    COFD[993] = 5.35267674E+00;
    COFD[994] = -4.69010505E-01;
    COFD[995] = 1.98979152E-02;
    COFD[996] = -1.63588981E+01;
    COFD[997] = 3.82388595E+00;
    COFD[998] = -2.84480724E-01;
    COFD[999] = 1.24506311E-02;
    COFD[1000] = -1.52861376E+01;
    COFD[1001] = 3.36790500E+00;
    COFD[1002] = -2.26321740E-01;
    COFD[1003] = 9.97135055E-03;
    COFD[1004] = -2.02710316E+01;
    COFD[1005] = 5.14984081E+00;
    COFD[1006] = -4.46093018E-01;
    COFD[1007] = 1.90396647E-02;
    COFD[1008] = -1.85899144E+01;
    COFD[1009] = 4.51052425E+00;
    COFD[1010] = -3.70301627E-01;
    COFD[1011] = 1.60416153E-02;
    COFD[1012] = -1.62824412E+01;
    COFD[1013] = 3.79163564E+00;
    COFD[1014] = -2.80257365E-01;
    COFD[1015] = 1.22656902E-02;
    COFD[1016] = -1.82285740E+01;
    COFD[1017] = 4.46848269E+00;
    COFD[1018] = -3.65269718E-01;
    COFD[1019] = 1.58407652E-02;
    COFD[1020] = -1.85990352E+01;
    COFD[1021] = 4.51052425E+00;
    COFD[1022] = -3.70301627E-01;
    COFD[1023] = 1.60416153E-02;
    COFD[1024] = -1.84863000E+01;
    COFD[1025] = 4.49330851E+00;
    COFD[1026] = -3.68208715E-01;
    COFD[1027] = 1.59565402E-02;
    COFD[1028] = -2.08347403E+01;
    COFD[1029] = 5.35267674E+00;
    COFD[1030] = -4.69010505E-01;
    COFD[1031] = 1.98979152E-02;
    COFD[1032] = -2.08642748E+01;
    COFD[1033] = 5.32244593E+00;
    COFD[1034] = -4.65829403E-01;
    COFD[1035] = 1.97895274E-02;
    COFD[1036] = -2.08686970E+01;
    COFD[1037] = 5.32244593E+00;
    COFD[1038] = -4.65829403E-01;
    COFD[1039] = 1.97895274E-02;
    COFD[1040] = -1.86641962E+01;
    COFD[1041] = 4.60874797E+00;
    COFD[1042] = -3.82368716E-01;
    COFD[1043] = 1.65370164E-02;
    COFD[1044] = -1.64338757E+01;
    COFD[1045] = 3.89309916E+00;
    COFD[1046] = -2.93528188E-01;
    COFD[1047] = 1.28463177E-02;
    COFD[1048] = -1.86748638E+01;
    COFD[1049] = 4.60874797E+00;
    COFD[1050] = -3.82368716E-01;
    COFD[1051] = 1.65370164E-02;
    COFD[1052] = -2.05408665E+01;
    COFD[1053] = 5.18417470E+00;
    COFD[1054] = -4.49491573E-01;
    COFD[1055] = 1.91438508E-02;
    COFD[1056] = -1.95081555E+01;
    COFD[1057] = 4.81575071E+00;
    COFD[1058] = -4.07042139E-01;
    COFD[1059] = 1.75187504E-02;
    COFD[1060] = -1.94088529E+01;
    COFD[1061] = 4.78708023E+00;
    COFD[1062] = -4.03693144E-01;
    COFD[1063] = 1.73884817E-02;
    COFD[1064] = -2.06706896E+01;
    COFD[1065] = 5.14159226E+00;
    COFD[1066] = -4.45395793E-01;
    COFD[1067] = 1.90248052E-02;
    COFD[1068] = -1.95374840E+01;
    COFD[1069] = 4.77151544E+00;
    COFD[1070] = -4.01882811E-01;
    COFD[1071] = 1.73184814E-02;
    COFD[1072] = -1.95374840E+01;
    COFD[1073] = 4.77151544E+00;
    COFD[1074] = -4.01882811E-01;
    COFD[1075] = 1.73184814E-02;
    COFD[1076] = -1.95425515E+01;
    COFD[1077] = 4.77151544E+00;
    COFD[1078] = -4.01882811E-01;
    COFD[1079] = 1.73184814E-02;
    COFD[1080] = -2.11103966E+01;
    COFD[1081] = 5.25151362E+00;
    COFD[1082] = -4.57337706E-01;
    COFD[1083] = 1.94489055E-02;
    COFD[1084] = -2.13011157E+01;
    COFD[1085] = 5.32167660E+00;
    COFD[1086] = -4.65740624E-01;
    COFD[1087] = 1.97861081E-02;
    COFD[1088] = -1.52486273E+01;
    COFD[1089] = 3.35922578E+00;
    COFD[1090] = -2.25181399E-01;
    COFD[1091] = 9.92132878E-03;
    COFD[1092] = -1.43190389E+01;
    COFD[1093] = 3.17651319E+00;
    COFD[1094] = -2.02028974E-01;
    COFD[1095] = 8.94232502E-03;
    COFD[1096] = -1.46554748E+01;
    COFD[1097] = 3.83606243E+00;
    COFD[1098] = -2.86076532E-01;
    COFD[1099] = 1.25205829E-02;
    COFD[1100] = -1.43394069E+01;
    COFD[1101] = 3.17651319E+00;
    COFD[1102] = -2.02028974E-01;
    COFD[1103] = 8.94232502E-03;
    COFD[1104] = -1.18998012E+01;
    COFD[1105] = 2.57507000E+00;
    COFD[1106] = -1.24033737E-01;
    COFD[1107] = 5.56694959E-03;
    COFD[1108] = -1.55588279E+01;
    COFD[1109] = 3.48070094E+00;
    COFD[1110] = -2.40859499E-01;
    COFD[1111] = 1.05972514E-02;
    COFD[1112] = -2.06463744E+01;
    COFD[1113] = 5.41688482E+00;
    COFD[1114] = -4.73387188E-01;
    COFD[1115] = 1.99280175E-02;
    COFD[1116] = -1.55741053E+01;
    COFD[1117] = 3.48070094E+00;
    COFD[1118] = -2.40859499E-01;
    COFD[1119] = 1.05972514E-02;
    COFD[1120] = -1.55666415E+01;
    COFD[1121] = 3.48070094E+00;
    COFD[1122] = -2.40859499E-01;
    COFD[1123] = 1.05972514E-02;
    COFD[1124] = -1.63301444E+01;
    COFD[1125] = 3.82388595E+00;
    COFD[1126] = -2.84480724E-01;
    COFD[1127] = 1.24506311E-02;
    COFD[1128] = -2.08367725E+01;
    COFD[1129] = 5.35267674E+00;
    COFD[1130] = -4.69010505E-01;
    COFD[1131] = 1.98979152E-02;
    COFD[1132] = -1.63542394E+01;
    COFD[1133] = 3.82388595E+00;
    COFD[1134] = -2.84480724E-01;
    COFD[1135] = 1.24506311E-02;
    COFD[1136] = -1.52792891E+01;
    COFD[1137] = 3.36790500E+00;
    COFD[1138] = -2.26321740E-01;
    COFD[1139] = 9.97135055E-03;
    COFD[1140] = -2.02637994E+01;
    COFD[1141] = 5.14984081E+00;
    COFD[1142] = -4.46093018E-01;
    COFD[1143] = 1.90396647E-02;
    COFD[1144] = -1.85829283E+01;
    COFD[1145] = 4.51052425E+00;
    COFD[1146] = -3.70301627E-01;
    COFD[1147] = 1.60416153E-02;
    COFD[1148] = -1.62775714E+01;
    COFD[1149] = 3.79163564E+00;
    COFD[1150] = -2.80257365E-01;
    COFD[1151] = 1.22656902E-02;
    COFD[1152] = -1.82217198E+01;
    COFD[1153] = 4.46848269E+00;
    COFD[1154] = -3.65269718E-01;
    COFD[1155] = 1.58407652E-02;
    COFD[1156] = -1.85919214E+01;
    COFD[1157] = 4.51052425E+00;
    COFD[1158] = -3.70301627E-01;
    COFD[1159] = 1.60416153E-02;
    COFD[1160] = -1.84777607E+01;
    COFD[1161] = 4.49330851E+00;
    COFD[1162] = -3.68208715E-01;
    COFD[1163] = 1.59565402E-02;
    COFD[1164] = -2.08277598E+01;
    COFD[1165] = 5.35267674E+00;
    COFD[1166] = -4.69010505E-01;
    COFD[1167] = 1.98979152E-02;
    COFD[1168] = -2.08554914E+01;
    COFD[1169] = 5.32244593E+00;
    COFD[1170] = -4.65829403E-01;
    COFD[1171] = 1.97895274E-02;
    COFD[1172] = -2.08598363E+01;
    COFD[1173] = 5.32244593E+00;
    COFD[1174] = -4.65829403E-01;
    COFD[1175] = 1.97895274E-02;
    COFD[1176] = -1.86576191E+01;
    COFD[1177] = 4.60874797E+00;
    COFD[1178] = -3.82368716E-01;
    COFD[1179] = 1.65370164E-02;
    COFD[1180] = -1.64255964E+01;
    COFD[1181] = 3.89309916E+00;
    COFD[1182] = -2.93528188E-01;
    COFD[1183] = 1.28463177E-02;
    COFD[1184] = -1.86681459E+01;
    COFD[1185] = 4.60874797E+00;
    COFD[1186] = -3.82368716E-01;
    COFD[1187] = 1.65370164E-02;
    COFD[1188] = -2.05324091E+01;
    COFD[1189] = 5.18417470E+00;
    COFD[1190] = -4.49491573E-01;
    COFD[1191] = 1.91438508E-02;
    COFD[1192] = -1.94998722E+01;
    COFD[1193] = 4.81575071E+00;
    COFD[1194] = -4.07042139E-01;
    COFD[1195] = 1.75187504E-02;
    COFD[1196] = -1.94004795E+01;
    COFD[1197] = 4.78708023E+00;
    COFD[1198] = -4.03693144E-01;
    COFD[1199] = 1.73884817E-02;
    COFD[1200] = -2.06612128E+01;
    COFD[1201] = 5.14159226E+00;
    COFD[1202] = -4.45395793E-01;
    COFD[1203] = 1.90248052E-02;
    COFD[1204] = -1.95290229E+01;
    COFD[1205] = 4.77151544E+00;
    COFD[1206] = -4.01882811E-01;
    COFD[1207] = 1.73184814E-02;
    COFD[1208] = -1.95290229E+01;
    COFD[1209] = 4.77151544E+00;
    COFD[1210] = -4.01882811E-01;
    COFD[1211] = 1.73184814E-02;
    COFD[1212] = -1.95340050E+01;
    COFD[1213] = 4.77151544E+00;
    COFD[1214] = -4.01882811E-01;
    COFD[1215] = 1.73184814E-02;
    COFD[1216] = -2.10999968E+01;
    COFD[1217] = 5.25151362E+00;
    COFD[1218] = -4.57337706E-01;
    COFD[1219] = 1.94489055E-02;
    COFD[1220] = -2.12907159E+01;
    COFD[1221] = 5.32167660E+00;
    COFD[1222] = -4.65740624E-01;
    COFD[1223] = 1.97861081E-02;
    COFD[1224] = -1.59404882E+01;
    COFD[1225] = 3.66853818E+00;
    COFD[1226] = -2.64346221E-01;
    COFD[1227] = 1.15784613E-02;
    COFD[1228] = -1.50584249E+01;
    COFD[1229] = 3.47945612E+00;
    COFD[1230] = -2.40703722E-01;
    COFD[1231] = 1.05907441E-02;
    COFD[1232] = -1.57972369E+01;
    COFD[1233] = 4.22225052E+00;
    COFD[1234] = -3.35156428E-01;
    COFD[1235] = 1.46104855E-02;
    COFD[1236] = -1.50724636E+01;
    COFD[1237] = 3.47945612E+00;
    COFD[1238] = -2.40703722E-01;
    COFD[1239] = 1.05907441E-02;
    COFD[1240] = -1.25098960E+01;
    COFD[1241] = 2.77873601E+00;
    COFD[1242] = -1.50637360E-01;
    COFD[1243] = 6.72684281E-03;
    COFD[1244] = -1.63254691E+01;
    COFD[1245] = 3.82388595E+00;
    COFD[1246] = -2.84480724E-01;
    COFD[1247] = 1.24506311E-02;
    COFD[1248] = -2.12639214E+01;
    COFD[1249] = 5.61184117E+00;
    COFD[1250] = -4.90532156E-01;
    COFD[1251] = 2.03507922E-02;
    COFD[1252] = -1.63345829E+01;
    COFD[1253] = 3.82388595E+00;
    COFD[1254] = -2.84480724E-01;
    COFD[1255] = 1.24506311E-02;
    COFD[1256] = -1.63301444E+01;
    COFD[1257] = 3.82388595E+00;
    COFD[1258] = -2.84480724E-01;
    COFD[1259] = 1.24506311E-02;
    COFD[1260] = -1.73027557E+01;
    COFD[1261] = 4.21416723E+00;
    COFD[1262] = -3.34163932E-01;
    COFD[1263] = 1.45697432E-02;
    COFD[1264] = -2.14215700E+01;
    COFD[1265] = 5.56531152E+00;
    COFD[1266] = -4.88789821E-01;
    COFD[1267] = 2.04437116E-02;
    COFD[1268] = -1.73198034E+01;
    COFD[1269] = 4.21416723E+00;
    COFD[1270] = -3.34163932E-01;
    COFD[1271] = 1.45697432E-02;
    COFD[1272] = -1.59634533E+01;
    COFD[1273] = 3.67388294E+00;
    COFD[1274] = -2.64990709E-01;
    COFD[1275] = 1.16042706E-02;
    COFD[1276] = -2.09376196E+01;
    COFD[1277] = 5.40870099E+00;
    COFD[1278] = -4.73017610E-01;
    COFD[1279] = 1.99399066E-02;
    COFD[1280] = -1.94530250E+01;
    COFD[1281] = 4.87180830E+00;
    COFD[1282] = -4.13582958E-01;
    COFD[1283] = 1.77726094E-02;
    COFD[1284] = -1.72556729E+01;
    COFD[1285] = 4.19029808E+00;
    COFD[1286] = -3.31177076E-01;
    COFD[1287] = 1.44446234E-02;
    COFD[1288] = -1.90996795E+01;
    COFD[1289] = 4.82869066E+00;
    COFD[1290] = -4.08564514E-01;
    COFD[1291] = 1.75784675E-02;
    COFD[1292] = -1.94585111E+01;
    COFD[1293] = 4.87180830E+00;
    COFD[1294] = -4.13582958E-01;
    COFD[1295] = 1.77726094E-02;
    COFD[1296] = -1.93015555E+01;
    COFD[1297] = 4.85015581E+00;
    COFD[1298] = -4.10945109E-01;
    COFD[1299] = 1.76651398E-02;
    COFD[1300] = -2.14160703E+01;
    COFD[1301] = 5.56531152E+00;
    COFD[1302] = -4.88789821E-01;
    COFD[1303] = 2.04437116E-02;
    COFD[1304] = -2.14048982E+01;
    COFD[1305] = 5.54007827E+00;
    COFD[1306] = -4.86434511E-01;
    COFD[1307] = 2.03779006E-02;
    COFD[1308] = -2.14073140E+01;
    COFD[1309] = 5.54007827E+00;
    COFD[1310] = -4.86434511E-01;
    COFD[1311] = 2.03779006E-02;
    COFD[1312] = -1.95548230E+01;
    COFD[1313] = 4.97133070E+00;
    COFD[1314] = -4.25604177E-01;
    COFD[1315] = 1.82582594E-02;
    COFD[1316] = -1.72572042E+01;
    COFD[1317] = 4.26063341E+00;
    COFD[1318] = -3.39848064E-01;
    COFD[1319] = 1.48021313E-02;
    COFD[1320] = -1.95613899E+01;
    COFD[1321] = 4.97133070E+00;
    COFD[1322] = -4.25604177E-01;
    COFD[1323] = 1.82582594E-02;
    COFD[1324] = -2.11378465E+01;
    COFD[1325] = 5.42846112E+00;
    COFD[1326] = -4.74321870E-01;
    COFD[1327] = 1.99459749E-02;
    COFD[1328] = -2.03111230E+01;
    COFD[1329] = 5.15740122E+00;
    COFD[1330] = -4.46644818E-01;
    COFD[1331] = 1.90459001E-02;
    COFD[1332] = -2.02434438E+01;
    COFD[1333] = 5.14418672E+00;
    COFD[1334] = -4.45631004E-01;
    COFD[1335] = 1.90308403E-02;
    COFD[1336] = -2.12684314E+01;
    COFD[1337] = 5.40206122E+00;
    COFD[1338] = -4.72555229E-01;
    COFD[1339] = 1.99358199E-02;
    COFD[1340] = -2.03711787E+01;
    COFD[1341] = 5.13279789E+00;
    COFD[1342] = -4.44474174E-01;
    COFD[1343] = 1.89937678E-02;
    COFD[1344] = -2.03711787E+01;
    COFD[1345] = 5.13279789E+00;
    COFD[1346] = -4.44474174E-01;
    COFD[1347] = 1.89937678E-02;
    COFD[1348] = -2.03739934E+01;
    COFD[1349] = 5.13279789E+00;
    COFD[1350] = -4.44474174E-01;
    COFD[1351] = 1.89937678E-02;
    COFD[1352] = -2.16464729E+01;
    COFD[1353] = 5.49325617E+00;
    COFD[1354] = -4.81586493E-01;
    COFD[1355] = 2.02162204E-02;
    COFD[1356] = -2.17867314E+01;
    COFD[1357] = 5.53950393E+00;
    COFD[1358] = -4.86376204E-01;
    COFD[1359] = 2.03760106E-02;
    COFD[1360] = -2.04833713E+01;
    COFD[1361] = 5.23112374E+00;
    COFD[1362] = -4.54967682E-01;
    COFD[1363] = 1.93570423E-02;
    COFD[1364] = -1.94373127E+01;
    COFD[1365] = 5.02567894E+00;
    COFD[1366] = -4.32045169E-01;
    COFD[1367] = 1.85132214E-02;
    COFD[1368] = -1.97550088E+01;
    COFD[1369] = 5.56931926E+00;
    COFD[1370] = -4.89105511E-01;
    COFD[1371] = 2.04493129E-02;
    COFD[1372] = -1.94570287E+01;
    COFD[1373] = 5.02567894E+00;
    COFD[1374] = -4.32045169E-01;
    COFD[1375] = 1.85132214E-02;
    COFD[1376] = -1.60528285E+01;
    COFD[1377] = 4.11188603E+00;
    COFD[1378] = -3.21540884E-01;
    COFD[1379] = 1.40482564E-02;
    COFD[1380] = -2.08293255E+01;
    COFD[1381] = 5.35267674E+00;
    COFD[1382] = -4.69010505E-01;
    COFD[1383] = 1.98979152E-02;
    COFD[1384] = -1.77563250E+01;
    COFD[1385] = 3.57475686E+00;
    COFD[1386] = -1.56396297E-01;
    COFD[1387] = 3.12157721E-03;
    COFD[1388] = -2.08438809E+01;
    COFD[1389] = 5.35267674E+00;
    COFD[1390] = -4.69010505E-01;
    COFD[1391] = 1.98979152E-02;
    COFD[1392] = -2.08367725E+01;
    COFD[1393] = 5.35267674E+00;
    COFD[1394] = -4.69010505E-01;
    COFD[1395] = 1.98979152E-02;
    COFD[1396] = -2.14215700E+01;
    COFD[1397] = 5.56531152E+00;
    COFD[1398] = -4.88789821E-01;
    COFD[1399] = 2.04437116E-02;
    COFD[1400] = -1.90499441E+01;
    COFD[1401] = 3.99221757E+00;
    COFD[1402] = -2.19854880E-01;
    COFD[1403] = 6.22736279E-03;
    COFD[1404] = -2.14449559E+01;
    COFD[1405] = 5.56531152E+00;
    COFD[1406] = -4.88789821E-01;
    COFD[1407] = 2.04437116E-02;
    COFD[1408] = -2.05128705E+01;
    COFD[1409] = 5.23843909E+00;
    COFD[1410] = -4.55815614E-01;
    COFD[1411] = 1.93898040E-02;
    COFD[1412] = -2.01889168E+01;
    COFD[1413] = 4.53183330E+00;
    COFD[1414] = -3.02186760E-01;
    COFD[1415] = 1.02756490E-02;
    COFD[1416] = -2.19700270E+01;
    COFD[1417] = 5.43750833E+00;
    COFD[1418] = -4.50273329E-01;
    COFD[1419] = 1.79013718E-02;
    COFD[1420] = -2.14082453E+01;
    COFD[1421] = 5.55346617E+00;
    COFD[1422] = -4.87783156E-01;
    COFD[1423] = 2.04210886E-02;
    COFD[1424] = -2.17855148E+01;
    COFD[1425] = 5.47519298E+00;
    COFD[1426] = -4.57113040E-01;
    COFD[1427] = 1.82758312E-02;
    COFD[1428] = -2.19786173E+01;
    COFD[1429] = 5.43750833E+00;
    COFD[1430] = -4.50273329E-01;
    COFD[1431] = 1.79013718E-02;
    COFD[1432] = -2.19317743E+01;
    COFD[1433] = 5.45216133E+00;
    COFD[1434] = -4.52916925E-01;
    COFD[1435] = 1.80456400E-02;
    COFD[1436] = -1.90413348E+01;
    COFD[1437] = 3.99221757E+00;
    COFD[1438] = -2.19854880E-01;
    COFD[1439] = 6.22736279E-03;
    COFD[1440] = -1.94051843E+01;
    COFD[1441] = 4.10954793E+00;
    COFD[1442] = -2.37523329E-01;
    COFD[1443] = 7.08858141E-03;
    COFD[1444] = -1.94092888E+01;
    COFD[1445] = 4.10954793E+00;
    COFD[1446] = -2.37523329E-01;
    COFD[1447] = 7.08858141E-03;
    COFD[1448] = -2.16798265E+01;
    COFD[1449] = 5.36811769E+00;
    COFD[1450] = -4.37727086E-01;
    COFD[1451] = 1.72167686E-02;
    COFD[1452] = -2.14303479E+01;
    COFD[1453] = 5.59268435E+00;
    COFD[1454] = -4.91232974E-01;
    COFD[1455] = 2.05064746E-02;
    COFD[1456] = -2.16899073E+01;
    COFD[1457] = 5.36811769E+00;
    COFD[1458] = -4.37727086E-01;
    COFD[1459] = 1.72167686E-02;
    COFD[1460] = -2.01064363E+01;
    COFD[1461] = 4.41511629E+00;
    COFD[1462] = -2.84086963E-01;
    COFD[1463] = 9.37586971E-03;
    COFD[1464] = -2.15258568E+01;
    COFD[1465] = 5.12799307E+00;
    COFD[1466] = -3.96938732E-01;
    COFD[1467] = 1.50673195E-02;
    COFD[1468] = -2.15802788E+01;
    COFD[1469] = 5.16868516E+00;
    COFD[1470] = -4.03721581E-01;
    COFD[1471] = 1.54206640E-02;
    COFD[1472] = -2.06753031E+01;
    COFD[1473] = 4.56779427E+00;
    COFD[1474] = -3.07785839E-01;
    COFD[1475] = 1.05545767E-02;
    COFD[1476] = -2.17926864E+01;
    COFD[1477] = 5.19232842E+00;
    COFD[1478] = -4.07643284E-01;
    COFD[1479] = 1.56246434E-02;
    COFD[1480] = -2.17926864E+01;
    COFD[1481] = 5.19232842E+00;
    COFD[1482] = -4.07643284E-01;
    COFD[1483] = 1.56246434E-02;
    COFD[1484] = -2.17974021E+01;
    COFD[1485] = 5.19232842E+00;
    COFD[1486] = -4.07643284E-01;
    COFD[1487] = 1.56246434E-02;
    COFD[1488] = -2.01633664E+01;
    COFD[1489] = 4.26535185E+00;
    COFD[1490] = -2.61127236E-01;
    COFD[1491] = 8.24369619E-03;
    COFD[1492] = -1.98359760E+01;
    COFD[1493] = 4.11158627E+00;
    COFD[1494] = -2.37831519E-01;
    COFD[1495] = 7.10363413E-03;
    COFD[1496] = -1.59633387E+01;
    COFD[1497] = 3.66853818E+00;
    COFD[1498] = -2.64346221E-01;
    COFD[1499] = 1.15784613E-02;
    COFD[1500] = -1.50766130E+01;
    COFD[1501] = 3.47945612E+00;
    COFD[1502] = -2.40703722E-01;
    COFD[1503] = 1.05907441E-02;
    COFD[1504] = -1.57994893E+01;
    COFD[1505] = 4.22225052E+00;
    COFD[1506] = -3.35156428E-01;
    COFD[1507] = 1.46104855E-02;
    COFD[1508] = -1.50911794E+01;
    COFD[1509] = 3.47945612E+00;
    COFD[1510] = -2.40703722E-01;
    COFD[1511] = 1.05907441E-02;
    COFD[1512] = -1.25141260E+01;
    COFD[1513] = 2.77873601E+00;
    COFD[1514] = -1.50637360E-01;
    COFD[1515] = 6.72684281E-03;
    COFD[1516] = -1.63493345E+01;
    COFD[1517] = 3.82388595E+00;
    COFD[1518] = -2.84480724E-01;
    COFD[1519] = 1.24506311E-02;
    COFD[1520] = -2.12831323E+01;
    COFD[1521] = 5.61184117E+00;
    COFD[1522] = -4.90532156E-01;
    COFD[1523] = 2.03507922E-02;
    COFD[1524] = -1.63588981E+01;
    COFD[1525] = 3.82388595E+00;
    COFD[1526] = -2.84480724E-01;
    COFD[1527] = 1.24506311E-02;
    COFD[1528] = -1.63542394E+01;
    COFD[1529] = 3.82388595E+00;
    COFD[1530] = -2.84480724E-01;
    COFD[1531] = 1.24506311E-02;
    COFD[1532] = -1.73198034E+01;
    COFD[1533] = 4.21416723E+00;
    COFD[1534] = -3.34163932E-01;
    COFD[1535] = 1.45697432E-02;
    COFD[1536] = -2.14449559E+01;
    COFD[1537] = 5.56531152E+00;
    COFD[1538] = -4.88789821E-01;
    COFD[1539] = 2.04437116E-02;
    COFD[1540] = -1.73374529E+01;
    COFD[1541] = 4.21416723E+00;
    COFD[1542] = -3.34163932E-01;
    COFD[1543] = 1.45697432E-02;
    COFD[1544] = -1.59863030E+01;
    COFD[1545] = 3.67388294E+00;
    COFD[1546] = -2.64990709E-01;
    COFD[1547] = 1.16042706E-02;
    COFD[1548] = -2.09612557E+01;
    COFD[1549] = 5.40870099E+00;
    COFD[1550] = -4.73017610E-01;
    COFD[1551] = 1.99399066E-02;
    COFD[1552] = -1.94761606E+01;
    COFD[1553] = 4.87180830E+00;
    COFD[1554] = -4.13582958E-01;
    COFD[1555] = 1.77726094E-02;
    COFD[1556] = -1.72738845E+01;
    COFD[1557] = 4.19029808E+00;
    COFD[1558] = -3.31177076E-01;
    COFD[1559] = 1.44446234E-02;
    COFD[1560] = -1.91225414E+01;
    COFD[1561] = 4.82869066E+00;
    COFD[1562] = -4.08564514E-01;
    COFD[1563] = 1.75784675E-02;
    COFD[1564] = -1.94819080E+01;
    COFD[1565] = 4.87180830E+00;
    COFD[1566] = -4.13582958E-01;
    COFD[1567] = 1.77726094E-02;
    COFD[1568] = -1.93276434E+01;
    COFD[1569] = 4.85015581E+00;
    COFD[1570] = -4.10945109E-01;
    COFD[1571] = 1.76651398E-02;
    COFD[1572] = -2.14391943E+01;
    COFD[1573] = 5.56531152E+00;
    COFD[1574] = -4.88789821E-01;
    COFD[1575] = 2.04437116E-02;
    COFD[1576] = -2.14314090E+01;
    COFD[1577] = 5.54007827E+00;
    COFD[1578] = -4.86434511E-01;
    COFD[1579] = 2.03779006E-02;
    COFD[1580] = -2.14339566E+01;
    COFD[1581] = 5.54007827E+00;
    COFD[1582] = -4.86434511E-01;
    COFD[1583] = 2.03779006E-02;
    COFD[1584] = -1.95770968E+01;
    COFD[1585] = 4.97133070E+00;
    COFD[1586] = -4.25604177E-01;
    COFD[1587] = 1.82582594E-02;
    COFD[1588] = -1.72828302E+01;
    COFD[1589] = 4.26063341E+00;
    COFD[1590] = -3.39848064E-01;
    COFD[1591] = 1.48021313E-02;
    COFD[1592] = -1.95839648E+01;
    COFD[1593] = 4.97133070E+00;
    COFD[1594] = -4.25604177E-01;
    COFD[1595] = 1.82582594E-02;
    COFD[1596] = -2.11637902E+01;
    COFD[1597] = 5.42846112E+00;
    COFD[1598] = -4.74321870E-01;
    COFD[1599] = 1.99459749E-02;
    COFD[1600] = -2.03367561E+01;
    COFD[1601] = 5.15740122E+00;
    COFD[1602] = -4.46644818E-01;
    COFD[1603] = 1.90459001E-02;
    COFD[1604] = -2.02692384E+01;
    COFD[1605] = 5.14418672E+00;
    COFD[1606] = -4.45631004E-01;
    COFD[1607] = 1.90308403E-02;
    COFD[1608] = -2.12960899E+01;
    COFD[1609] = 5.40206122E+00;
    COFD[1610] = -4.72555229E-01;
    COFD[1611] = 1.99358199E-02;
    COFD[1612] = -2.03971290E+01;
    COFD[1613] = 5.13279789E+00;
    COFD[1614] = -4.44474174E-01;
    COFD[1615] = 1.89937678E-02;
    COFD[1616] = -2.03971290E+01;
    COFD[1617] = 5.13279789E+00;
    COFD[1618] = -4.44474174E-01;
    COFD[1619] = 1.89937678E-02;
    COFD[1620] = -2.04000941E+01;
    COFD[1621] = 5.13279789E+00;
    COFD[1622] = -4.44474174E-01;
    COFD[1623] = 1.89937678E-02;
    COFD[1624] = -2.16755464E+01;
    COFD[1625] = 5.49325617E+00;
    COFD[1626] = -4.81586493E-01;
    COFD[1627] = 2.02162204E-02;
    COFD[1628] = -2.18158049E+01;
    COFD[1629] = 5.53950393E+00;
    COFD[1630] = -4.86376204E-01;
    COFD[1631] = 2.03760106E-02;
    COFD[1632] = -1.50031687E+01;
    COFD[1633] = 3.26223357E+00;
    COFD[1634] = -2.12746642E-01;
    COFD[1635] = 9.38912883E-03;
    COFD[1636] = -1.40999008E+01;
    COFD[1637] = 3.08120012E+00;
    COFD[1638] = -1.89629903E-01;
    COFD[1639] = 8.40361952E-03;
    COFD[1640] = -1.43151174E+01;
    COFD[1641] = 3.68038508E+00;
    COFD[1642] = -2.65779346E-01;
    COFD[1643] = 1.16360771E-02;
    COFD[1644] = -1.41191261E+01;
    COFD[1645] = 3.08120012E+00;
    COFD[1646] = -1.89629903E-01;
    COFD[1647] = 8.40361952E-03;
    COFD[1648] = -1.17159737E+01;
    COFD[1649] = 2.48123210E+00;
    COFD[1650] = -1.11322604E-01;
    COFD[1651] = 4.99282389E-03;
    COFD[1652] = -1.52721107E+01;
    COFD[1653] = 3.36790500E+00;
    COFD[1654] = -2.26321740E-01;
    COFD[1655] = 9.97135055E-03;
    COFD[1656] = -2.11388331E+01;
    COFD[1657] = 5.55529675E+00;
    COFD[1658] = -4.87942518E-01;
    COFD[1659] = 2.04249054E-02;
    COFD[1660] = -1.52861376E+01;
    COFD[1661] = 3.36790500E+00;
    COFD[1662] = -2.26321740E-01;
    COFD[1663] = 9.97135055E-03;
    COFD[1664] = -1.52792891E+01;
    COFD[1665] = 3.36790500E+00;
    COFD[1666] = -2.26321740E-01;
    COFD[1667] = 9.97135055E-03;
    COFD[1668] = -1.59634533E+01;
    COFD[1669] = 3.67388294E+00;
    COFD[1670] = -2.64990709E-01;
    COFD[1671] = 1.16042706E-02;
    COFD[1672] = -2.05128705E+01;
    COFD[1673] = 5.23843909E+00;
    COFD[1674] = -4.55815614E-01;
    COFD[1675] = 1.93898040E-02;
    COFD[1676] = -1.59863030E+01;
    COFD[1677] = 3.67388294E+00;
    COFD[1678] = -2.64990709E-01;
    COFD[1679] = 1.16042706E-02;
    COFD[1680] = -1.50233475E+01;
    COFD[1681] = 3.26660767E+00;
    COFD[1682] = -2.13287177E-01;
    COFD[1683] = 9.41137857E-03;
    COFD[1684] = -2.02642227E+01;
    COFD[1685] = 5.14499740E+00;
    COFD[1686] = -4.45694430E-01;
    COFD[1687] = 1.90318646E-02;
    COFD[1688] = -1.82872310E+01;
    COFD[1689] = 4.40289649E+00;
    COFD[1690] = -3.57289765E-01;
    COFD[1691] = 1.55166804E-02;
    COFD[1692] = -1.59525102E+01;
    COFD[1693] = 3.66023858E+00;
    COFD[1694] = -2.63401043E-01;
    COFD[1695] = 1.15432000E-02;
    COFD[1696] = -1.79116531E+01;
    COFD[1697] = 4.35148286E+00;
    COFD[1698] = -3.50886647E-01;
    COFD[1699] = 1.52498573E-02;
    COFD[1700] = -1.82955252E+01;
    COFD[1701] = 4.40289649E+00;
    COFD[1702] = -3.57289765E-01;
    COFD[1703] = 1.55166804E-02;
    COFD[1704] = -1.81735763E+01;
    COFD[1705] = 4.38391495E+00;
    COFD[1706] = -3.54941287E-01;
    COFD[1707] = 1.54195107E-02;
    COFD[1708] = -2.05045578E+01;
    COFD[1709] = 5.23843909E+00;
    COFD[1710] = -4.55815614E-01;
    COFD[1711] = 1.93898040E-02;
    COFD[1712] = -2.04949373E+01;
    COFD[1713] = 5.19614628E+00;
    COFD[1714] = -4.50889164E-01;
    COFD[1715] = 1.91983328E-02;
    COFD[1716] = -2.04988684E+01;
    COFD[1717] = 5.19614628E+00;
    COFD[1718] = -4.50889164E-01;
    COFD[1719] = 1.91983328E-02;
    COFD[1720] = -1.83296965E+01;
    COFD[1721] = 4.48570999E+00;
    COFD[1722] = -3.67301524E-01;
    COFD[1723] = 1.59204254E-02;
    COFD[1724] = -1.60261675E+01;
    COFD[1725] = 3.73312045E+00;
    COFD[1726] = -2.72579779E-01;
    COFD[1727] = 1.19290272E-02;
    COFD[1728] = -1.83394481E+01;
    COFD[1729] = 4.48570999E+00;
    COFD[1730] = -3.67301524E-01;
    COFD[1731] = 1.59204254E-02;
    COFD[1732] = -2.02969740E+01;
    COFD[1733] = 5.11106992E+00;
    COFD[1734] = -4.42047129E-01;
    COFD[1735] = 1.89042990E-02;
    COFD[1736] = -1.92044492E+01;
    COFD[1737] = 4.71304783E+00;
    COFD[1738] = -3.94942083E-01;
    COFD[1739] = 1.70435959E-02;
    COFD[1740] = -1.91118445E+01;
    COFD[1741] = 4.68715685E+00;
    COFD[1742] = -3.91979493E-01;
    COFD[1743] = 1.69314004E-02;
    COFD[1744] = -2.03525027E+01;
    COFD[1745] = 5.04005588E+00;
    COFD[1746] = -4.33725091E-01;
    COFD[1747] = 1.85786663E-02;
    COFD[1748] = -1.92334028E+01;
    COFD[1749] = 4.67033934E+00;
    COFD[1750] = -3.89971551E-01;
    COFD[1751] = 1.68513441E-02;
    COFD[1752] = -1.92334028E+01;
    COFD[1753] = 4.67033934E+00;
    COFD[1754] = -3.89971551E-01;
    COFD[1755] = 1.68513441E-02;
    COFD[1756] = -1.92379258E+01;
    COFD[1757] = 4.67033934E+00;
    COFD[1758] = -3.89971551E-01;
    COFD[1759] = 1.68513441E-02;
    COFD[1760] = -2.08660787E+01;
    COFD[1761] = 5.17548471E+00;
    COFD[1762] = -4.48496897E-01;
    COFD[1763] = 1.91060173E-02;
    COFD[1764] = -2.10310742E+01;
    COFD[1765] = 5.23485505E+00;
    COFD[1766] = -4.55400362E-01;
    COFD[1767] = 1.93737680E-02;
    COFD[1768] = -2.02268902E+01;
    COFD[1769] = 5.13632093E+00;
    COFD[1770] = -4.44839124E-01;
    COFD[1771] = 1.90058354E-02;
    COFD[1772] = -1.88179418E+01;
    COFD[1773] = 4.79683898E+00;
    COFD[1774] = -4.04829719E-01;
    COFD[1775] = 1.74325475E-02;
    COFD[1776] = -1.92718582E+01;
    COFD[1777] = 5.41172124E+00;
    COFD[1778] = -4.73213887E-01;
    COFD[1779] = 1.99405473E-02;
    COFD[1780] = -1.88378874E+01;
    COFD[1781] = 4.79683898E+00;
    COFD[1782] = -4.04829719E-01;
    COFD[1783] = 1.74325475E-02;
    COFD[1784] = -1.58456300E+01;
    COFD[1785] = 4.02074783E+00;
    COFD[1786] = -3.10018522E-01;
    COFD[1787] = 1.35599552E-02;
    COFD[1788] = -2.04928958E+01;
    COFD[1789] = 5.22397933E+00;
    COFD[1790] = -4.54138171E-01;
    COFD[1791] = 1.93249285E-02;
    COFD[1792] = -1.65295288E+01;
    COFD[1793] = 2.97569206E+00;
    COFD[1794] = -6.75652842E-02;
    COFD[1795] = -1.08648422E-03;
    COFD[1796] = -2.02710316E+01;
    COFD[1797] = 5.14984081E+00;
    COFD[1798] = -4.46093018E-01;
    COFD[1799] = 1.90396647E-02;
    COFD[1800] = -2.02637994E+01;
    COFD[1801] = 5.14984081E+00;
    COFD[1802] = -4.46093018E-01;
    COFD[1803] = 1.90396647E-02;
    COFD[1804] = -2.09376196E+01;
    COFD[1805] = 5.40870099E+00;
    COFD[1806] = -4.73017610E-01;
    COFD[1807] = 1.99399066E-02;
    COFD[1808] = -2.01889168E+01;
    COFD[1809] = 4.53183330E+00;
    COFD[1810] = -3.02186760E-01;
    COFD[1811] = 1.02756490E-02;
    COFD[1812] = -2.09612557E+01;
    COFD[1813] = 5.40870099E+00;
    COFD[1814] = -4.73017610E-01;
    COFD[1815] = 1.99399066E-02;
    COFD[1816] = -2.02642227E+01;
    COFD[1817] = 5.14499740E+00;
    COFD[1818] = -4.45694430E-01;
    COFD[1819] = 1.90318646E-02;
    COFD[1820] = -1.95877017E+01;
    COFD[1821] = 4.27643051E+00;
    COFD[1822] = -2.68040901E-01;
    COFD[1823] = 8.77650113E-03;
    COFD[1824] = -2.21384805E+01;
    COFD[1825] = 5.56656297E+00;
    COFD[1826] = -4.75500048E-01;
    COFD[1827] = 1.93332291E-02;
    COFD[1828] = -2.11381508E+01;
    COFD[1829] = 5.45574440E+00;
    COFD[1830] = -4.77436155E-01;
    COFD[1831] = 2.00644596E-02;
    COFD[1832] = -2.19136842E+01;
    COFD[1833] = 5.58503445E+00;
    COFD[1834] = -4.79552117E-01;
    COFD[1835] = 1.95750393E-02;
    COFD[1836] = -2.21472114E+01;
    COFD[1837] = 5.56656297E+00;
    COFD[1838] = -4.75500048E-01;
    COFD[1839] = 1.93332291E-02;
    COFD[1840] = -2.20421041E+01;
    COFD[1841] = 5.52708332E+00;
    COFD[1842] = -4.68000808E-01;
    COFD[1843] = 1.89131908E-02;
    COFD[1844] = -2.01801667E+01;
    COFD[1845] = 4.53183330E+00;
    COFD[1846] = -3.02186760E-01;
    COFD[1847] = 1.02756490E-02;
    COFD[1848] = -2.04451935E+01;
    COFD[1849] = 4.60682543E+00;
    COFD[1850] = -3.13971634E-01;
    COFD[1851] = 1.08661011E-02;
    COFD[1852] = -2.04493813E+01;
    COFD[1853] = 4.60682543E+00;
    COFD[1854] = -3.13971634E-01;
    COFD[1855] = 1.08661011E-02;
    COFD[1856] = -2.18848136E+01;
    COFD[1857] = 5.51302074E+00;
    COFD[1858] = -4.65263979E-01;
    COFD[1859] = 1.87580679E-02;
    COFD[1860] = -2.09241647E+01;
    COFD[1861] = 5.42316225E+00;
    COFD[1862] = -4.73702801E-01;
    COFD[1863] = 1.99217718E-02;
    COFD[1864] = -2.18950505E+01;
    COFD[1865] = 5.51302074E+00;
    COFD[1866] = -4.65263979E-01;
    COFD[1867] = 1.87580679E-02;
    COFD[1868] = -2.09272429E+01;
    COFD[1869] = 4.82184721E+00;
    COFD[1870] = -3.48128875E-01;
    COFD[1871] = 1.25918978E-02;
    COFD[1872] = -2.19617977E+01;
    COFD[1873] = 5.37170913E+00;
    COFD[1874] = -4.38338667E-01;
    COFD[1875] = 1.72490835E-02;
    COFD[1876] = -2.19873532E+01;
    COFD[1877] = 5.39977369E+00;
    COFD[1878] = -4.43340854E-01;
    COFD[1879] = 1.75199613E-02;
    COFD[1880] = -2.14535736E+01;
    COFD[1881] = 4.95733546E+00;
    COFD[1882] = -3.69505821E-01;
    COFD[1883] = 1.36687047E-02;
    COFD[1884] = -2.21713935E+01;
    COFD[1885] = 5.41196486E+00;
    COFD[1886] = -4.45632422E-01;
    COFD[1887] = 1.76474237E-02;
    COFD[1888] = -2.21713935E+01;
    COFD[1889] = 5.41196486E+00;
    COFD[1890] = -4.45632422E-01;
    COFD[1891] = 1.76474237E-02;
    COFD[1892] = -2.21762017E+01;
    COFD[1893] = 5.41196486E+00;
    COFD[1894] = -4.45632422E-01;
    COFD[1895] = 1.76474237E-02;
    COFD[1896] = -2.03466999E+01;
    COFD[1897] = 4.38969526E+00;
    COFD[1898] = -2.84615363E-01;
    COFD[1899] = 9.55641838E-03;
    COFD[1900] = -2.01613414E+01;
    COFD[1901] = 4.29679630E+00;
    COFD[1902] = -2.69916064E-01;
    COFD[1903] = 8.81737046E-03;
    COFD[1904] = -1.82590824E+01;
    COFD[1905] = 4.39538102E+00;
    COFD[1906] = -3.56367230E-01;
    COFD[1907] = 1.54788461E-02;
    COFD[1908] = -1.72053106E+01;
    COFD[1909] = 4.15807461E+00;
    COFD[1910] = -3.27178539E-01;
    COFD[1911] = 1.42784349E-02;
    COFD[1912] = -1.78631557E+01;
    COFD[1913] = 4.88268692E+00;
    COFD[1914] = -4.14917638E-01;
    COFD[1915] = 1.78274298E-02;
    COFD[1916] = -1.72247972E+01;
    COFD[1917] = 4.15807461E+00;
    COFD[1918] = -3.27178539E-01;
    COFD[1919] = 1.42784349E-02;
    COFD[1920] = -1.39648112E+01;
    COFD[1921] = 3.24966086E+00;
    COFD[1922] = -2.11199992E-01;
    COFD[1923] = 9.32580661E-03;
    COFD[1924] = -1.85756076E+01;
    COFD[1925] = 4.51052425E+00;
    COFD[1926] = -3.70301627E-01;
    COFD[1927] = 1.60416153E-02;
    COFD[1928] = -2.13084334E+01;
    COFD[1929] = 5.27210469E+00;
    COFD[1930] = -4.21419216E-01;
    COFD[1931] = 1.63567178E-02;
    COFD[1932] = -1.85899144E+01;
    COFD[1933] = 4.51052425E+00;
    COFD[1934] = -3.70301627E-01;
    COFD[1935] = 1.60416153E-02;
    COFD[1936] = -1.85829283E+01;
    COFD[1937] = 4.51052425E+00;
    COFD[1938] = -3.70301627E-01;
    COFD[1939] = 1.60416153E-02;
    COFD[1940] = -1.94530250E+01;
    COFD[1941] = 4.87180830E+00;
    COFD[1942] = -4.13582958E-01;
    COFD[1943] = 1.77726094E-02;
    COFD[1944] = -2.19700270E+01;
    COFD[1945] = 5.43750833E+00;
    COFD[1946] = -4.50273329E-01;
    COFD[1947] = 1.79013718E-02;
    COFD[1948] = -1.94761606E+01;
    COFD[1949] = 4.87180830E+00;
    COFD[1950] = -4.13582958E-01;
    COFD[1951] = 1.77726094E-02;
    COFD[1952] = -1.82872310E+01;
    COFD[1953] = 4.40289649E+00;
    COFD[1954] = -3.57289765E-01;
    COFD[1955] = 1.55166804E-02;
    COFD[1956] = -2.21384805E+01;
    COFD[1957] = 5.56656297E+00;
    COFD[1958] = -4.75500048E-01;
    COFD[1959] = 1.93332291E-02;
    COFD[1960] = -2.14737305E+01;
    COFD[1961] = 5.41585806E+00;
    COFD[1962] = -4.73359323E-01;
    COFD[1963] = 1.99310239E-02;
    COFD[1964] = -1.94126575E+01;
    COFD[1965] = 4.84669430E+00;
    COFD[1966] = -4.10571455E-01;
    COFD[1967] = 1.76520543E-02;
    COFD[1968] = -2.11931178E+01;
    COFD[1969] = 5.40060531E+00;
    COFD[1970] = -4.72449699E-01;
    COFD[1971] = 1.99345817E-02;
    COFD[1972] = -2.14821817E+01;
    COFD[1973] = 5.41585806E+00;
    COFD[1974] = -4.73359323E-01;
    COFD[1975] = 1.99310239E-02;
    COFD[1976] = -2.14049543E+01;
    COFD[1977] = 5.41122754E+00;
    COFD[1978] = -4.73185889E-01;
    COFD[1979] = 1.99407905E-02;
    COFD[1980] = -2.19615570E+01;
    COFD[1981] = 5.43750833E+00;
    COFD[1982] = -4.50273329E-01;
    COFD[1983] = 1.79013718E-02;
    COFD[1984] = -2.21229141E+01;
    COFD[1985] = 5.47072190E+00;
    COFD[1986] = -4.56301261E-01;
    COFD[1987] = 1.82313566E-02;
    COFD[1988] = -2.21269367E+01;
    COFD[1989] = 5.47072190E+00;
    COFD[1990] = -4.56301261E-01;
    COFD[1991] = 1.82313566E-02;
    COFD[1992] = -2.15126310E+01;
    COFD[1993] = 5.48426911E+00;
    COFD[1994] = -4.80606512E-01;
    COFD[1995] = 2.01811046E-02;
    COFD[1996] = -1.95737308E+01;
    COFD[1997] = 4.93449043E+00;
    COFD[1998] = -4.21243802E-01;
    COFD[1999] = 1.80859966E-02;
    COFD[2000] = -2.15225573E+01;
    COFD[2001] = 5.48426911E+00;
    COFD[2002] = -4.80606512E-01;
    COFD[2003] = 2.01811046E-02;
    COFD[2004] = -2.22377812E+01;
    COFD[2005] = 5.53139819E+00;
    COFD[2006] = -4.68828555E-01;
    COFD[2007] = 1.89597887E-02;
    COFD[2008] = -2.21694197E+01;
    COFD[2009] = 5.60403905E+00;
    COFD[2010] = -4.91221691E-01;
    COFD[2011] = 2.04473483E-02;
    COFD[2012] = -2.21242889E+01;
    COFD[2013] = 5.60010742E+00;
    COFD[2014] = -4.91597429E-01;
    COFD[2015] = 2.04987718E-02;
    COFD[2016] = -2.25375432E+01;
    COFD[2017] = 5.57713269E+00;
    COFD[2018] = -4.77555529E-01;
    COFD[2019] = 1.94497781E-02;
    COFD[2020] = -2.22600844E+01;
    COFD[2021] = 5.59632316E+00;
    COFD[2022] = -4.91568011E-01;
    COFD[2023] = 2.05156966E-02;
    COFD[2024] = -2.22600844E+01;
    COFD[2025] = 5.59632316E+00;
    COFD[2026] = -4.91568011E-01;
    COFD[2027] = 2.05156966E-02;
    COFD[2028] = -2.22647093E+01;
    COFD[2029] = 5.59632316E+00;
    COFD[2030] = -4.91568011E-01;
    COFD[2031] = 2.05156966E-02;
    COFD[2032] = -2.25524874E+01;
    COFD[2033] = 5.50138987E+00;
    COFD[2034] = -4.62503139E-01;
    COFD[2035] = 1.85883715E-02;
    COFD[2036] = -2.25180193E+01;
    COFD[2037] = 5.47136127E+00;
    COFD[2038] = -4.56417141E-01;
    COFD[2039] = 1.82376994E-02;
    COFD[2040] = -1.59327297E+01;
    COFD[2041] = 3.65620899E+00;
    COFD[2042] = -2.62933804E-01;
    COFD[2043] = 1.15253223E-02;
    COFD[2044] = -1.50270339E+01;
    COFD[2045] = 3.46140064E+00;
    COFD[2046] = -2.38440092E-01;
    COFD[2047] = 1.04960087E-02;
    COFD[2048] = -1.57199037E+01;
    COFD[2049] = 4.19936335E+00;
    COFD[2050] = -3.32311009E-01;
    COFD[2051] = 1.44921003E-02;
    COFD[2052] = -1.50420953E+01;
    COFD[2053] = 3.46140064E+00;
    COFD[2054] = -2.38440092E-01;
    COFD[2055] = 1.04960087E-02;
    COFD[2056] = -1.24693568E+01;
    COFD[2057] = 2.76686648E+00;
    COFD[2058] = -1.49120141E-01;
    COFD[2059] = 6.66220432E-03;
    COFD[2060] = -1.62724462E+01;
    COFD[2061] = 3.79163564E+00;
    COFD[2062] = -2.80257365E-01;
    COFD[2063] = 1.22656902E-02;
    COFD[2064] = -2.14087397E+01;
    COFD[2065] = 5.57282008E+00;
    COFD[2066] = -4.76690890E-01;
    COFD[2067] = 1.94000719E-02;
    COFD[2068] = -1.62824412E+01;
    COFD[2069] = 3.79163564E+00;
    COFD[2070] = -2.80257365E-01;
    COFD[2071] = 1.22656902E-02;
    COFD[2072] = -1.62775714E+01;
    COFD[2073] = 3.79163564E+00;
    COFD[2074] = -2.80257365E-01;
    COFD[2075] = 1.22656902E-02;
    COFD[2076] = -1.72556729E+01;
    COFD[2077] = 4.19029808E+00;
    COFD[2078] = -3.31177076E-01;
    COFD[2079] = 1.44446234E-02;
    COFD[2080] = -2.14082453E+01;
    COFD[2081] = 5.55346617E+00;
    COFD[2082] = -4.87783156E-01;
    COFD[2083] = 2.04210886E-02;
    COFD[2084] = -1.72738845E+01;
    COFD[2085] = 4.19029808E+00;
    COFD[2086] = -3.31177076E-01;
    COFD[2087] = 1.44446234E-02;
    COFD[2088] = -1.59525102E+01;
    COFD[2089] = 3.66023858E+00;
    COFD[2090] = -2.63401043E-01;
    COFD[2091] = 1.15432000E-02;
    COFD[2092] = -2.11381508E+01;
    COFD[2093] = 5.45574440E+00;
    COFD[2094] = -4.77436155E-01;
    COFD[2095] = 2.00644596E-02;
    COFD[2096] = -1.94126575E+01;
    COFD[2097] = 4.84669430E+00;
    COFD[2098] = -4.10571455E-01;
    COFD[2099] = 1.76520543E-02;
    COFD[2100] = -1.72167708E+01;
    COFD[2101] = 4.16886779E+00;
    COFD[2102] = -3.28518156E-01;
    COFD[2103] = 1.43341626E-02;
    COFD[2104] = -1.90692595E+01;
    COFD[2105] = 4.80830699E+00;
    COFD[2106] = -4.06171933E-01;
    COFD[2107] = 1.74848791E-02;
    COFD[2108] = -1.94186547E+01;
    COFD[2109] = 4.84669430E+00;
    COFD[2110] = -4.10571455E-01;
    COFD[2111] = 1.76520543E-02;
    COFD[2112] = -1.92867554E+01;
    COFD[2113] = 4.83375900E+00;
    COFD[2114] = -4.09146560E-01;
    COFD[2115] = 1.76006599E-02;
    COFD[2116] = -2.14022336E+01;
    COFD[2117] = 5.55346617E+00;
    COFD[2118] = -4.87783156E-01;
    COFD[2119] = 2.04210886E-02;
    COFD[2120] = -2.13881945E+01;
    COFD[2121] = 5.52422470E+00;
    COFD[2122] = -4.84872944E-01;
    COFD[2123] = 2.03298213E-02;
    COFD[2124] = -2.13908698E+01;
    COFD[2125] = 5.52422470E+00;
    COFD[2126] = -4.84872944E-01;
    COFD[2127] = 2.03298213E-02;
    COFD[2128] = -1.95154079E+01;
    COFD[2129] = 4.94787350E+00;
    COFD[2130] = -4.22829292E-01;
    COFD[2131] = 1.81487163E-02;
    COFD[2132] = -1.72316148E+01;
    COFD[2133] = 4.24011069E+00;
    COFD[2134] = -3.37339810E-01;
    COFD[2135] = 1.46996679E-02;
    COFD[2136] = -1.95225629E+01;
    COFD[2137] = 4.94787350E+00;
    COFD[2138] = -4.22829292E-01;
    COFD[2139] = 1.81487163E-02;
    COFD[2140] = -2.11341653E+01;
    COFD[2141] = 5.41773516E+00;
    COFD[2142] = -4.73414338E-01;
    COFD[2143] = 1.99258685E-02;
    COFD[2144] = -2.03123540E+01;
    COFD[2145] = 5.14854169E+00;
    COFD[2146] = -4.45984343E-01;
    COFD[2147] = 1.90374217E-02;
    COFD[2148] = -2.02318658E+01;
    COFD[2149] = 5.12963391E+00;
    COFD[2150] = -4.44146826E-01;
    COFD[2151] = 1.89829640E-02;
    COFD[2152] = -2.12758591E+01;
    COFD[2153] = 5.39400772E+00;
    COFD[2154] = -4.72026046E-01;
    COFD[2155] = 1.99336599E-02;
    COFD[2156] = -2.03526104E+01;
    COFD[2157] = 5.11453301E+00;
    COFD[2158] = -4.42447016E-01;
    COFD[2159] = 1.89196698E-02;
    COFD[2160] = -2.03526104E+01;
    COFD[2161] = 5.11453301E+00;
    COFD[2162] = -4.42447016E-01;
    COFD[2163] = 1.89196698E-02;
    COFD[2164] = -2.03557208E+01;
    COFD[2165] = 5.11453301E+00;
    COFD[2166] = -4.42447016E-01;
    COFD[2167] = 1.89196698E-02;
    COFD[2168] = -2.17407707E+01;
    COFD[2169] = 5.50870827E+00;
    COFD[2170] = -4.83264910E-01;
    COFD[2171] = 2.02760941E-02;
    COFD[2172] = -2.18731920E+01;
    COFD[2173] = 5.55171660E+00;
    COFD[2174] = -4.87609504E-01;
    COFD[2175] = 2.04156590E-02;
    COFD[2176] = -1.78815889E+01;
    COFD[2177] = 4.34347890E+00;
    COFD[2178] = -3.49890003E-01;
    COFD[2179] = 1.52083459E-02;
    COFD[2180] = -1.68343393E+01;
    COFD[2181] = 4.11954900E+00;
    COFD[2182] = -3.22470391E-01;
    COFD[2183] = 1.40859564E-02;
    COFD[2184] = -1.74407963E+01;
    COFD[2185] = 4.83580036E+00;
    COFD[2186] = -4.09383573E-01;
    COFD[2187] = 1.76098175E-02;
    COFD[2188] = -1.68535757E+01;
    COFD[2189] = 4.11954900E+00;
    COFD[2190] = -3.22470391E-01;
    COFD[2191] = 1.40859564E-02;
    COFD[2192] = -1.36336373E+01;
    COFD[2193] = 3.22088176E+00;
    COFD[2194] = -2.07623790E-01;
    COFD[2195] = 9.17771542E-03;
    COFD[2196] = -1.82145353E+01;
    COFD[2197] = 4.46848269E+00;
    COFD[2198] = -3.65269718E-01;
    COFD[2199] = 1.58407652E-02;
    COFD[2200] = -2.11309197E+01;
    COFD[2201] = 5.32644193E+00;
    COFD[2202] = -4.30581064E-01;
    COFD[2203] = 1.68379725E-02;
    COFD[2204] = -1.82285740E+01;
    COFD[2205] = 4.46848269E+00;
    COFD[2206] = -3.65269718E-01;
    COFD[2207] = 1.58407652E-02;
    COFD[2208] = -1.82217198E+01;
    COFD[2209] = 4.46848269E+00;
    COFD[2210] = -3.65269718E-01;
    COFD[2211] = 1.58407652E-02;
    COFD[2212] = -1.90996795E+01;
    COFD[2213] = 4.82869066E+00;
    COFD[2214] = -4.08564514E-01;
    COFD[2215] = 1.75784675E-02;
    COFD[2216] = -2.17855148E+01;
    COFD[2217] = 5.47519298E+00;
    COFD[2218] = -4.57113040E-01;
    COFD[2219] = 1.82758312E-02;
    COFD[2220] = -1.91225414E+01;
    COFD[2221] = 4.82869066E+00;
    COFD[2222] = -4.08564514E-01;
    COFD[2223] = 1.75784675E-02;
    COFD[2224] = -1.79116531E+01;
    COFD[2225] = 4.35148286E+00;
    COFD[2226] = -3.50886647E-01;
    COFD[2227] = 1.52498573E-02;
    COFD[2228] = -2.19136842E+01;
    COFD[2229] = 5.58503445E+00;
    COFD[2230] = -4.79552117E-01;
    COFD[2231] = 1.95750393E-02;
    COFD[2232] = -2.11931178E+01;
    COFD[2233] = 5.40060531E+00;
    COFD[2234] = -4.72449699E-01;
    COFD[2235] = 1.99345817E-02;
    COFD[2236] = -1.90692595E+01;
    COFD[2237] = 4.80830699E+00;
    COFD[2238] = -4.06171933E-01;
    COFD[2239] = 1.74848791E-02;
    COFD[2240] = -2.08820897E+01;
    COFD[2241] = 5.38250415E+00;
    COFD[2242] = -4.71144140E-01;
    COFD[2243] = 1.99199779E-02;
    COFD[2244] = -2.12014186E+01;
    COFD[2245] = 5.40060531E+00;
    COFD[2246] = -4.72449699E-01;
    COFD[2247] = 1.99345817E-02;
    COFD[2248] = -2.11031143E+01;
    COFD[2249] = 5.39439999E+00;
    COFD[2250] = -4.72050184E-01;
    COFD[2251] = 1.99336257E-02;
    COFD[2252] = -2.17771954E+01;
    COFD[2253] = 5.47519298E+00;
    COFD[2254] = -4.57113040E-01;
    COFD[2255] = 1.82758312E-02;
    COFD[2256] = -2.19162360E+01;
    COFD[2257] = 5.49906960E+00;
    COFD[2258] = -4.61793001E-01;
    COFD[2259] = 1.85415189E-02;
    COFD[2260] = -2.19201709E+01;
    COFD[2261] = 5.49906960E+00;
    COFD[2262] = -4.61793001E-01;
    COFD[2263] = 1.85415189E-02;
    COFD[2264] = -2.11427744E+01;
    COFD[2265] = 5.43893233E+00;
    COFD[2266] = -4.75546039E-01;
    COFD[2267] = 1.99938690E-02;
    COFD[2268] = -1.91367023E+01;
    COFD[2269] = 4.87703209E+00;
    COFD[2270] = -4.14222202E-01;
    COFD[2271] = 1.77987878E-02;
    COFD[2272] = -2.11525334E+01;
    COFD[2273] = 5.43893233E+00;
    COFD[2274] = -4.75546039E-01;
    COFD[2275] = 1.99938690E-02;
    COFD[2276] = -2.20445411E+01;
    COFD[2277] = 5.56049839E+00;
    COFD[2278] = -4.74367872E-01;
    COFD[2279] = 1.92702787E-02;
    COFD[2280] = -2.19032561E+01;
    COFD[2281] = 5.59794138E+00;
    COFD[2282] = -4.91684532E-01;
    COFD[2283] = 2.05170953E-02;
    COFD[2284] = -2.18222696E+01;
    COFD[2285] = 5.57940140E+00;
    COFD[2286] = -4.89964112E-01;
    COFD[2287] = 2.04689539E-02;
    COFD[2288] = -2.23144051E+01;
    COFD[2289] = 5.58484440E+00;
    COFD[2290] = -4.80068319E-01;
    COFD[2291] = 1.96187346E-02;
    COFD[2292] = -2.19617258E+01;
    COFD[2293] = 5.57026255E+00;
    COFD[2294] = -4.89178491E-01;
    COFD[2295] = 2.04505218E-02;
    COFD[2296] = -2.19617258E+01;
    COFD[2297] = 5.57026255E+00;
    COFD[2298] = -4.89178491E-01;
    COFD[2299] = 2.04505218E-02;
    COFD[2300] = -2.19662531E+01;
    COFD[2301] = 5.57026255E+00;
    COFD[2302] = -4.89178491E-01;
    COFD[2303] = 2.04505218E-02;
    COFD[2304] = -2.23572782E+01;
    COFD[2305] = 5.52017185E+00;
    COFD[2306] = -4.66659932E-01;
    COFD[2307] = 1.88373100E-02;
    COFD[2308] = -2.23434237E+01;
    COFD[2309] = 5.49927389E+00;
    COFD[2310] = -4.61845436E-01;
    COFD[2311] = 1.85448066E-02;
    COFD[2312] = -1.82673770E+01;
    COFD[2313] = 4.39538102E+00;
    COFD[2314] = -3.56367230E-01;
    COFD[2315] = 1.54788461E-02;
    COFD[2316] = -1.72112971E+01;
    COFD[2317] = 4.15807461E+00;
    COFD[2318] = -3.27178539E-01;
    COFD[2319] = 1.42784349E-02;
    COFD[2320] = -1.78637178E+01;
    COFD[2321] = 4.88268692E+00;
    COFD[2322] = -4.14917638E-01;
    COFD[2323] = 1.78274298E-02;
    COFD[2324] = -1.72310232E+01;
    COFD[2325] = 4.15807461E+00;
    COFD[2326] = -3.27178539E-01;
    COFD[2327] = 1.42784349E-02;
    COFD[2328] = -1.39658996E+01;
    COFD[2329] = 3.24966086E+00;
    COFD[2330] = -2.11199992E-01;
    COFD[2331] = 9.32580661E-03;
    COFD[2332] = -1.85844688E+01;
    COFD[2333] = 4.51052425E+00;
    COFD[2334] = -3.70301627E-01;
    COFD[2335] = 1.60416153E-02;
    COFD[2336] = -2.13148887E+01;
    COFD[2337] = 5.27210469E+00;
    COFD[2338] = -4.21419216E-01;
    COFD[2339] = 1.63567178E-02;
    COFD[2340] = -1.85990352E+01;
    COFD[2341] = 4.51052425E+00;
    COFD[2342] = -3.70301627E-01;
    COFD[2343] = 1.60416153E-02;
    COFD[2344] = -1.85919214E+01;
    COFD[2345] = 4.51052425E+00;
    COFD[2346] = -3.70301627E-01;
    COFD[2347] = 1.60416153E-02;
    COFD[2348] = -1.94585111E+01;
    COFD[2349] = 4.87180830E+00;
    COFD[2350] = -4.13582958E-01;
    COFD[2351] = 1.77726094E-02;
    COFD[2352] = -2.19786173E+01;
    COFD[2353] = 5.43750833E+00;
    COFD[2354] = -4.50273329E-01;
    COFD[2355] = 1.79013718E-02;
    COFD[2356] = -1.94819080E+01;
    COFD[2357] = 4.87180830E+00;
    COFD[2358] = -4.13582958E-01;
    COFD[2359] = 1.77726094E-02;
    COFD[2360] = -1.82955252E+01;
    COFD[2361] = 4.40289649E+00;
    COFD[2362] = -3.57289765E-01;
    COFD[2363] = 1.55166804E-02;
    COFD[2364] = -2.21472114E+01;
    COFD[2365] = 5.56656297E+00;
    COFD[2366] = -4.75500048E-01;
    COFD[2367] = 1.93332291E-02;
    COFD[2368] = -2.14821817E+01;
    COFD[2369] = 5.41585806E+00;
    COFD[2370] = -4.73359323E-01;
    COFD[2371] = 1.99310239E-02;
    COFD[2372] = -1.94186547E+01;
    COFD[2373] = 4.84669430E+00;
    COFD[2374] = -4.10571455E-01;
    COFD[2375] = 1.76520543E-02;
    COFD[2376] = -2.12014186E+01;
    COFD[2377] = 5.40060531E+00;
    COFD[2378] = -4.72449699E-01;
    COFD[2379] = 1.99345817E-02;
    COFD[2380] = -2.14907782E+01;
    COFD[2381] = 5.41585806E+00;
    COFD[2382] = -4.73359323E-01;
    COFD[2383] = 1.99310239E-02;
    COFD[2384] = -2.14151520E+01;
    COFD[2385] = 5.41122754E+00;
    COFD[2386] = -4.73185889E-01;
    COFD[2387] = 1.99407905E-02;
    COFD[2388] = -2.19700018E+01;
    COFD[2389] = 5.43750833E+00;
    COFD[2390] = -4.50273329E-01;
    COFD[2391] = 1.79013718E-02;
    COFD[2392] = -2.21333822E+01;
    COFD[2393] = 5.47072190E+00;
    COFD[2394] = -4.56301261E-01;
    COFD[2395] = 1.82313566E-02;
    COFD[2396] = -2.21374903E+01;
    COFD[2397] = 5.47072190E+00;
    COFD[2398] = -4.56301261E-01;
    COFD[2399] = 1.82313566E-02;
    COFD[2400] = -2.15206146E+01;
    COFD[2401] = 5.48426911E+00;
    COFD[2402] = -4.80606512E-01;
    COFD[2403] = 2.01811046E-02;
    COFD[2404] = -1.95836394E+01;
    COFD[2405] = 4.93449043E+00;
    COFD[2406] = -4.21243802E-01;
    COFD[2407] = 1.80859966E-02;
    COFD[2408] = -2.15307023E+01;
    COFD[2409] = 5.48426911E+00;
    COFD[2410] = -4.80606512E-01;
    COFD[2411] = 2.01811046E-02;
    COFD[2412] = -2.22478879E+01;
    COFD[2413] = 5.53139819E+00;
    COFD[2414] = -4.68828555E-01;
    COFD[2415] = 1.89597887E-02;
    COFD[2416] = -2.21793326E+01;
    COFD[2417] = 5.60403905E+00;
    COFD[2418] = -4.91221691E-01;
    COFD[2419] = 2.04473483E-02;
    COFD[2420] = -2.21343023E+01;
    COFD[2421] = 5.60010742E+00;
    COFD[2422] = -4.91597429E-01;
    COFD[2423] = 2.04987718E-02;
    COFD[2424] = -2.25487737E+01;
    COFD[2425] = 5.57713269E+00;
    COFD[2426] = -4.77555529E-01;
    COFD[2427] = 1.94497781E-02;
    COFD[2428] = -2.22701953E+01;
    COFD[2429] = 5.59632316E+00;
    COFD[2430] = -4.91568011E-01;
    COFD[2431] = 2.05156966E-02;
    COFD[2432] = -2.22701953E+01;
    COFD[2433] = 5.59632316E+00;
    COFD[2434] = -4.91568011E-01;
    COFD[2435] = 2.05156966E-02;
    COFD[2436] = -2.22749151E+01;
    COFD[2437] = 5.59632316E+00;
    COFD[2438] = -4.91568011E-01;
    COFD[2439] = 2.05156966E-02;
    COFD[2440] = -2.25647193E+01;
    COFD[2441] = 5.50138987E+00;
    COFD[2442] = -4.62503139E-01;
    COFD[2443] = 1.85883715E-02;
    COFD[2444] = -2.25302512E+01;
    COFD[2445] = 5.47136127E+00;
    COFD[2446] = -4.56417141E-01;
    COFD[2447] = 1.82376994E-02;
    COFD[2448] = -1.81432461E+01;
    COFD[2449] = 4.37565431E+00;
    COFD[2450] = -3.53906025E-01;
    COFD[2451] = 1.53760786E-02;
    COFD[2452] = -1.70534856E+01;
    COFD[2453] = 4.14240922E+00;
    COFD[2454] = -3.25239774E-01;
    COFD[2455] = 1.41980687E-02;
    COFD[2456] = -1.76147026E+01;
    COFD[2457] = 4.86049500E+00;
    COFD[2458] = -4.12200578E-01;
    COFD[2459] = 1.77160971E-02;
    COFD[2460] = -1.70757047E+01;
    COFD[2461] = 4.14240922E+00;
    COFD[2462] = -3.25239774E-01;
    COFD[2463] = 1.41980687E-02;
    COFD[2464] = -1.37794315E+01;
    COFD[2465] = 3.23973858E+00;
    COFD[2466] = -2.09989036E-01;
    COFD[2467] = 9.27667906E-03;
    COFD[2468] = -1.84688406E+01;
    COFD[2469] = 4.49330851E+00;
    COFD[2470] = -3.68208715E-01;
    COFD[2471] = 1.59565402E-02;
    COFD[2472] = -2.07653719E+01;
    COFD[2473] = 5.01092022E+00;
    COFD[2474] = -3.77985635E-01;
    COFD[2475] = 1.40968645E-02;
    COFD[2476] = -1.84863000E+01;
    COFD[2477] = 4.49330851E+00;
    COFD[2478] = -3.68208715E-01;
    COFD[2479] = 1.59565402E-02;
    COFD[2480] = -1.84777607E+01;
    COFD[2481] = 4.49330851E+00;
    COFD[2482] = -3.68208715E-01;
    COFD[2483] = 1.59565402E-02;
    COFD[2484] = -1.93015555E+01;
    COFD[2485] = 4.85015581E+00;
    COFD[2486] = -4.10945109E-01;
    COFD[2487] = 1.76651398E-02;
    COFD[2488] = -2.19317743E+01;
    COFD[2489] = 5.45216133E+00;
    COFD[2490] = -4.52916925E-01;
    COFD[2491] = 1.80456400E-02;
    COFD[2492] = -1.93276434E+01;
    COFD[2493] = 4.85015581E+00;
    COFD[2494] = -4.10945109E-01;
    COFD[2495] = 1.76651398E-02;
    COFD[2496] = -1.81735763E+01;
    COFD[2497] = 4.38391495E+00;
    COFD[2498] = -3.54941287E-01;
    COFD[2499] = 1.54195107E-02;
    COFD[2500] = -2.20421041E+01;
    COFD[2501] = 5.52708332E+00;
    COFD[2502] = -4.68000808E-01;
    COFD[2503] = 1.89131908E-02;
    COFD[2504] = -2.14049543E+01;
    COFD[2505] = 5.41122754E+00;
    COFD[2506] = -4.73185889E-01;
    COFD[2507] = 1.99407905E-02;
    COFD[2508] = -1.92867554E+01;
    COFD[2509] = 4.83375900E+00;
    COFD[2510] = -4.09146560E-01;
    COFD[2511] = 1.76006599E-02;
    COFD[2512] = -2.11031143E+01;
    COFD[2513] = 5.39439999E+00;
    COFD[2514] = -4.72050184E-01;
    COFD[2515] = 1.99336257E-02;
    COFD[2516] = -2.14151520E+01;
    COFD[2517] = 5.41122754E+00;
    COFD[2518] = -4.73185889E-01;
    COFD[2519] = 1.99407905E-02;
    COFD[2520] = -2.13425698E+01;
    COFD[2521] = 5.40460130E+00;
    COFD[2522] = -4.72718910E-01;
    COFD[2523] = 1.99362717E-02;
    COFD[2524] = -2.19215555E+01;
    COFD[2525] = 5.45216133E+00;
    COFD[2526] = -4.52916925E-01;
    COFD[2527] = 1.80456400E-02;
    COFD[2528] = -2.21083035E+01;
    COFD[2529] = 5.48540187E+00;
    COFD[2530] = -4.58962148E-01;
    COFD[2531] = 1.83770355E-02;
    COFD[2532] = -2.21134005E+01;
    COFD[2533] = 5.48540187E+00;
    COFD[2534] = -4.58962148E-01;
    COFD[2535] = 1.83770355E-02;
    COFD[2536] = -2.13961414E+01;
    COFD[2537] = 5.46685775E+00;
    COFD[2538] = -4.78665416E-01;
    COFD[2539] = 2.01093915E-02;
    COFD[2540] = -1.94485982E+01;
    COFD[2541] = 4.91446566E+00;
    COFD[2542] = -4.18837152E-01;
    COFD[2543] = 1.79893537E-02;
    COFD[2544] = -2.14079882E+01;
    COFD[2545] = 5.46685775E+00;
    COFD[2546] = -4.78665416E-01;
    COFD[2547] = 2.01093915E-02;
    COFD[2548] = -2.22176950E+01;
    COFD[2549] = 5.54251230E+00;
    COFD[2550] = -4.70946314E-01;
    COFD[2551] = 1.90785869E-02;
    COFD[2552] = -2.21216828E+01;
    COFD[2553] = 5.60203389E+00;
    COFD[2554] = -4.91444416E-01;
    COFD[2555] = 2.04761886E-02;
    COFD[2556] = -2.20725883E+01;
    COFD[2557] = 5.59642965E+00;
    COFD[2558] = -4.91577716E-01;
    COFD[2559] = 2.05159582E-02;
    COFD[2560] = -2.25320751E+01;
    COFD[2561] = 5.58240011E+00;
    COFD[2562] = -4.78844918E-01;
    COFD[2563] = 1.95298191E-02;
    COFD[2564] = -2.22052004E+01;
    COFD[2565] = 5.58604166E+00;
    COFD[2566] = -4.90602184E-01;
    COFD[2567] = 2.04880352E-02;
    COFD[2568] = -2.22052004E+01;
    COFD[2569] = 5.58604166E+00;
    COFD[2570] = -4.90602184E-01;
    COFD[2571] = 2.04880352E-02;
    COFD[2572] = -2.22110089E+01;
    COFD[2573] = 5.58604166E+00;
    COFD[2574] = -4.90602184E-01;
    COFD[2575] = 2.04880352E-02;
    COFD[2576] = -2.25695646E+01;
    COFD[2577] = 5.49903771E+00;
    COFD[2578] = -4.61784834E-01;
    COFD[2579] = 1.85410072E-02;
    COFD[2580] = -2.25168081E+01;
    COFD[2581] = 5.46125558E+00;
    COFD[2582] = -4.54580949E-01;
    COFD[2583] = 1.81370928E-02;
    COFD[2584] = -2.04750581E+01;
    COFD[2585] = 5.23112374E+00;
    COFD[2586] = -4.54967682E-01;
    COFD[2587] = 1.93570423E-02;
    COFD[2588] = -1.94313116E+01;
    COFD[2589] = 5.02567894E+00;
    COFD[2590] = -4.32045169E-01;
    COFD[2591] = 1.85132214E-02;
    COFD[2592] = -1.97544450E+01;
    COFD[2593] = 5.56931926E+00;
    COFD[2594] = -4.89105511E-01;
    COFD[2595] = 2.04493129E-02;
    COFD[2596] = -1.94507876E+01;
    COFD[2597] = 5.02567894E+00;
    COFD[2598] = -4.32045169E-01;
    COFD[2599] = 1.85132214E-02;
    COFD[2600] = -1.60517370E+01;
    COFD[2601] = 4.11188603E+00;
    COFD[2602] = -3.21540884E-01;
    COFD[2603] = 1.40482564E-02;
    COFD[2604] = -2.08204449E+01;
    COFD[2605] = 5.35267674E+00;
    COFD[2606] = -4.69010505E-01;
    COFD[2607] = 1.98979152E-02;
    COFD[2608] = -1.77498543E+01;
    COFD[2609] = 3.57475686E+00;
    COFD[2610] = -1.56396297E-01;
    COFD[2611] = 3.12157721E-03;
    COFD[2612] = -2.08347403E+01;
    COFD[2613] = 5.35267674E+00;
    COFD[2614] = -4.69010505E-01;
    COFD[2615] = 1.98979152E-02;
    COFD[2616] = -2.08277598E+01;
    COFD[2617] = 5.35267674E+00;
    COFD[2618] = -4.69010505E-01;
    COFD[2619] = 1.98979152E-02;
    COFD[2620] = -2.14160703E+01;
    COFD[2621] = 5.56531152E+00;
    COFD[2622] = -4.88789821E-01;
    COFD[2623] = 2.04437116E-02;
    COFD[2624] = -1.90413348E+01;
    COFD[2625] = 3.99221757E+00;
    COFD[2626] = -2.19854880E-01;
    COFD[2627] = 6.22736279E-03;
    COFD[2628] = -2.14391943E+01;
    COFD[2629] = 5.56531152E+00;
    COFD[2630] = -4.88789821E-01;
    COFD[2631] = 2.04437116E-02;
    COFD[2632] = -2.05045578E+01;
    COFD[2633] = 5.23843909E+00;
    COFD[2634] = -4.55815614E-01;
    COFD[2635] = 1.93898040E-02;
    COFD[2636] = -2.01801667E+01;
    COFD[2637] = 4.53183330E+00;
    COFD[2638] = -3.02186760E-01;
    COFD[2639] = 1.02756490E-02;
    COFD[2640] = -2.19615570E+01;
    COFD[2641] = 5.43750833E+00;
    COFD[2642] = -4.50273329E-01;
    COFD[2643] = 1.79013718E-02;
    COFD[2644] = -2.14022336E+01;
    COFD[2645] = 5.55346617E+00;
    COFD[2646] = -4.87783156E-01;
    COFD[2647] = 2.04210886E-02;
    COFD[2648] = -2.17771954E+01;
    COFD[2649] = 5.47519298E+00;
    COFD[2650] = -4.57113040E-01;
    COFD[2651] = 1.82758312E-02;
    COFD[2652] = -2.19700018E+01;
    COFD[2653] = 5.43750833E+00;
    COFD[2654] = -4.50273329E-01;
    COFD[2655] = 1.79013718E-02;
    COFD[2656] = -2.19215555E+01;
    COFD[2657] = 5.45216133E+00;
    COFD[2658] = -4.52916925E-01;
    COFD[2659] = 1.80456400E-02;
    COFD[2660] = -1.90328712E+01;
    COFD[2661] = 3.99221757E+00;
    COFD[2662] = -2.19854880E-01;
    COFD[2663] = 6.22736279E-03;
    COFD[2664] = -1.93946947E+01;
    COFD[2665] = 4.10954793E+00;
    COFD[2666] = -2.37523329E-01;
    COFD[2667] = 7.08858141E-03;
    COFD[2668] = -1.93987136E+01;
    COFD[2669] = 4.10954793E+00;
    COFD[2670] = -2.37523329E-01;
    COFD[2671] = 7.08858141E-03;
    COFD[2672] = -2.16718247E+01;
    COFD[2673] = 5.36811769E+00;
    COFD[2674] = -4.37727086E-01;
    COFD[2675] = 1.72167686E-02;
    COFD[2676] = -2.14204185E+01;
    COFD[2677] = 5.59268435E+00;
    COFD[2678] = -4.91232974E-01;
    COFD[2679] = 2.05064746E-02;
    COFD[2680] = -2.16817439E+01;
    COFD[2681] = 5.36811769E+00;
    COFD[2682] = -4.37727086E-01;
    COFD[2683] = 1.72167686E-02;
    COFD[2684] = -2.00963085E+01;
    COFD[2685] = 4.41511629E+00;
    COFD[2686] = -2.84086963E-01;
    COFD[2687] = 9.37586971E-03;
    COFD[2688] = -2.15159231E+01;
    COFD[2689] = 5.12799307E+00;
    COFD[2690] = -3.96938732E-01;
    COFD[2691] = 1.50673195E-02;
    COFD[2692] = -2.15702446E+01;
    COFD[2693] = 5.16868516E+00;
    COFD[2694] = -4.03721581E-01;
    COFD[2695] = 1.54206640E-02;
    COFD[2696] = -2.06640503E+01;
    COFD[2697] = 4.56779427E+00;
    COFD[2698] = -3.07785839E-01;
    COFD[2699] = 1.05545767E-02;
    COFD[2700] = -2.17825544E+01;
    COFD[2701] = 5.19232842E+00;
    COFD[2702] = -4.07643284E-01;
    COFD[2703] = 1.56246434E-02;
    COFD[2704] = -2.17825544E+01;
    COFD[2705] = 5.19232842E+00;
    COFD[2706] = -4.07643284E-01;
    COFD[2707] = 1.56246434E-02;
    COFD[2708] = -2.17871751E+01;
    COFD[2709] = 5.19232842E+00;
    COFD[2710] = -4.07643284E-01;
    COFD[2711] = 1.56246434E-02;
    COFD[2712] = -2.01511112E+01;
    COFD[2713] = 4.26535185E+00;
    COFD[2714] = -2.61127236E-01;
    COFD[2715] = 8.24369619E-03;
    COFD[2716] = -1.98237209E+01;
    COFD[2717] = 4.11158627E+00;
    COFD[2718] = -2.37831519E-01;
    COFD[2719] = 7.10363413E-03;
    COFD[2720] = -2.04649069E+01;
    COFD[2721] = 5.18856872E+00;
    COFD[2722] = -4.50001829E-01;
    COFD[2723] = 1.91636142E-02;
    COFD[2724] = -1.93925667E+01;
    COFD[2725] = 4.98286777E+00;
    COFD[2726] = -4.26970814E-01;
    COFD[2727] = 1.83122917E-02;
    COFD[2728] = -1.96914944E+01;
    COFD[2729] = 5.54637286E+00;
    COFD[2730] = -4.87070324E-01;
    COFD[2731] = 2.03983467E-02;
    COFD[2732] = -1.94151822E+01;
    COFD[2733] = 4.98286777E+00;
    COFD[2734] = -4.26970814E-01;
    COFD[2735] = 1.83122917E-02;
    COFD[2736] = -1.59632479E+01;
    COFD[2737] = 4.07051484E+00;
    COFD[2738] = -3.16303109E-01;
    COFD[2739] = 1.38259377E-02;
    COFD[2740] = -2.08463209E+01;
    COFD[2741] = 5.32244593E+00;
    COFD[2742] = -4.65829403E-01;
    COFD[2743] = 1.97895274E-02;
    COFD[2744] = -1.80862867E+01;
    COFD[2745] = 3.69199168E+00;
    COFD[2746] = -1.74005516E-01;
    COFD[2747] = 3.97694372E-03;
    COFD[2748] = -2.08642748E+01;
    COFD[2749] = 5.32244593E+00;
    COFD[2750] = -4.65829403E-01;
    COFD[2751] = 1.97895274E-02;
    COFD[2752] = -2.08554914E+01;
    COFD[2753] = 5.32244593E+00;
    COFD[2754] = -4.65829403E-01;
    COFD[2755] = 1.97895274E-02;
    COFD[2756] = -2.14048982E+01;
    COFD[2757] = 5.54007827E+00;
    COFD[2758] = -4.86434511E-01;
    COFD[2759] = 2.03779006E-02;
    COFD[2760] = -1.94051843E+01;
    COFD[2761] = 4.10954793E+00;
    COFD[2762] = -2.37523329E-01;
    COFD[2763] = 7.08858141E-03;
    COFD[2764] = -2.14314090E+01;
    COFD[2765] = 5.54007827E+00;
    COFD[2766] = -4.86434511E-01;
    COFD[2767] = 2.03779006E-02;
    COFD[2768] = -2.04949373E+01;
    COFD[2769] = 5.19614628E+00;
    COFD[2770] = -4.50889164E-01;
    COFD[2771] = 1.91983328E-02;
    COFD[2772] = -2.04451935E+01;
    COFD[2773] = 4.60682543E+00;
    COFD[2774] = -3.13971634E-01;
    COFD[2775] = 1.08661011E-02;
    COFD[2776] = -2.21229141E+01;
    COFD[2777] = 5.47072190E+00;
    COFD[2778] = -4.56301261E-01;
    COFD[2779] = 1.82313566E-02;
    COFD[2780] = -2.13881945E+01;
    COFD[2781] = 5.52422470E+00;
    COFD[2782] = -4.84872944E-01;
    COFD[2783] = 2.03298213E-02;
    COFD[2784] = -2.19162360E+01;
    COFD[2785] = 5.49906960E+00;
    COFD[2786] = -4.61793001E-01;
    COFD[2787] = 1.85415189E-02;
    COFD[2788] = -2.21333822E+01;
    COFD[2789] = 5.47072190E+00;
    COFD[2790] = -4.56301261E-01;
    COFD[2791] = 1.82313566E-02;
    COFD[2792] = -2.21083035E+01;
    COFD[2793] = 5.48540187E+00;
    COFD[2794] = -4.58962148E-01;
    COFD[2795] = 1.83770355E-02;
    COFD[2796] = -1.93946947E+01;
    COFD[2797] = 4.10954793E+00;
    COFD[2798] = -2.37523329E-01;
    COFD[2799] = 7.08858141E-03;
    COFD[2800] = -1.97704178E+01;
    COFD[2801] = 4.22062499E+00;
    COFD[2802] = -2.54326872E-01;
    COFD[2803] = 7.91017784E-03;
    COFD[2804] = -1.97756908E+01;
    COFD[2805] = 4.22062499E+00;
    COFD[2806] = -2.54326872E-01;
    COFD[2807] = 7.91017784E-03;
    COFD[2808] = -2.18318278E+01;
    COFD[2809] = 5.40298848E+00;
    COFD[2810] = -4.43954594E-01;
    COFD[2811] = 1.75542998E-02;
    COFD[2812] = -2.14782277E+01;
    COFD[2813] = 5.56978987E+00;
    COFD[2814] = -4.89141980E-01;
    COFD[2815] = 2.04499210E-02;
    COFD[2816] = -2.18439681E+01;
    COFD[2817] = 5.40298848E+00;
    COFD[2818] = -4.43954594E-01;
    COFD[2819] = 1.75542998E-02;
    COFD[2820] = -2.04096182E+01;
    COFD[2821] = 4.50250781E+00;
    COFD[2822] = -2.97622106E-01;
    COFD[2823] = 1.00481473E-02;
    COFD[2824] = -2.17407419E+01;
    COFD[2825] = 5.17945041E+00;
    COFD[2826] = -4.05514689E-01;
    COFD[2827] = 1.55141412E-02;
    COFD[2828] = -2.17934580E+01;
    COFD[2829] = 5.21869603E+00;
    COFD[2830] = -4.12084772E-01;
    COFD[2831] = 1.58573035E-02;
    COFD[2832] = -2.09502015E+01;
    COFD[2833] = 4.63716925E+00;
    COFD[2834] = -3.18815070E-01;
    COFD[2835] = 1.11115446E-02;
    COFD[2836] = -2.19919464E+01;
    COFD[2837] = 5.23595129E+00;
    COFD[2838] = -4.15079064E-01;
    COFD[2839] = 1.60168286E-02;
    COFD[2840] = -2.19919464E+01;
    COFD[2841] = 5.23595129E+00;
    COFD[2842] = -4.15079064E-01;
    COFD[2843] = 1.60168286E-02;
    COFD[2844] = -2.19979468E+01;
    COFD[2845] = 5.23595129E+00;
    COFD[2846] = -4.15079064E-01;
    COFD[2847] = 1.60168286E-02;
    COFD[2848] = -2.05128005E+01;
    COFD[2849] = 4.35981364E+00;
    COFD[2850] = -2.75585363E-01;
    COFD[2851] = 8.95565009E-03;
    COFD[2852] = -2.02246117E+01;
    COFD[2853] = 4.22278378E+00;
    COFD[2854] = -2.54653500E-01;
    COFD[2855] = 7.92616085E-03;
    COFD[2856] = -2.04688382E+01;
    COFD[2857] = 5.18856872E+00;
    COFD[2858] = -4.50001829E-01;
    COFD[2859] = 1.91636142E-02;
    COFD[2860] = -1.93952366E+01;
    COFD[2861] = 4.98286777E+00;
    COFD[2862] = -4.26970814E-01;
    COFD[2863] = 1.83122917E-02;
    COFD[2864] = -1.96917146E+01;
    COFD[2865] = 5.54637286E+00;
    COFD[2866] = -4.87070324E-01;
    COFD[2867] = 2.03983467E-02;
    COFD[2868] = -1.94179760E+01;
    COFD[2869] = 4.98286777E+00;
    COFD[2870] = -4.26970814E-01;
    COFD[2871] = 1.83122917E-02;
    COFD[2872] = -1.59636793E+01;
    COFD[2873] = 4.07051484E+00;
    COFD[2874] = -3.16303109E-01;
    COFD[2875] = 1.38259377E-02;
    COFD[2876] = -2.08505864E+01;
    COFD[2877] = 5.32244593E+00;
    COFD[2878] = -4.65829403E-01;
    COFD[2879] = 1.97895274E-02;
    COFD[2880] = -1.80892005E+01;
    COFD[2881] = 3.69199168E+00;
    COFD[2882] = -1.74005516E-01;
    COFD[2883] = 3.97694372E-03;
    COFD[2884] = -2.08686970E+01;
    COFD[2885] = 5.32244593E+00;
    COFD[2886] = -4.65829403E-01;
    COFD[2887] = 1.97895274E-02;
    COFD[2888] = -2.08598363E+01;
    COFD[2889] = 5.32244593E+00;
    COFD[2890] = -4.65829403E-01;
    COFD[2891] = 1.97895274E-02;
    COFD[2892] = -2.14073140E+01;
    COFD[2893] = 5.54007827E+00;
    COFD[2894] = -4.86434511E-01;
    COFD[2895] = 2.03779006E-02;
    COFD[2896] = -1.94092888E+01;
    COFD[2897] = 4.10954793E+00;
    COFD[2898] = -2.37523329E-01;
    COFD[2899] = 7.08858141E-03;
    COFD[2900] = -2.14339566E+01;
    COFD[2901] = 5.54007827E+00;
    COFD[2902] = -4.86434511E-01;
    COFD[2903] = 2.03779006E-02;
    COFD[2904] = -2.04988684E+01;
    COFD[2905] = 5.19614628E+00;
    COFD[2906] = -4.50889164E-01;
    COFD[2907] = 1.91983328E-02;
    COFD[2908] = -2.04493813E+01;
    COFD[2909] = 4.60682543E+00;
    COFD[2910] = -3.13971634E-01;
    COFD[2911] = 1.08661011E-02;
    COFD[2912] = -2.21269367E+01;
    COFD[2913] = 5.47072190E+00;
    COFD[2914] = -4.56301261E-01;
    COFD[2915] = 1.82313566E-02;
    COFD[2916] = -2.13908698E+01;
    COFD[2917] = 5.52422470E+00;
    COFD[2918] = -4.84872944E-01;
    COFD[2919] = 2.03298213E-02;
    COFD[2920] = -2.19201709E+01;
    COFD[2921] = 5.49906960E+00;
    COFD[2922] = -4.61793001E-01;
    COFD[2923] = 1.85415189E-02;
    COFD[2924] = -2.21374903E+01;
    COFD[2925] = 5.47072190E+00;
    COFD[2926] = -4.56301261E-01;
    COFD[2927] = 1.82313566E-02;
    COFD[2928] = -2.21134005E+01;
    COFD[2929] = 5.48540187E+00;
    COFD[2930] = -4.58962148E-01;
    COFD[2931] = 1.83770355E-02;
    COFD[2932] = -1.93987136E+01;
    COFD[2933] = 4.10954793E+00;
    COFD[2934] = -2.37523329E-01;
    COFD[2935] = 7.08858141E-03;
    COFD[2936] = -1.97756908E+01;
    COFD[2937] = 4.22062499E+00;
    COFD[2938] = -2.54326872E-01;
    COFD[2939] = 7.91017784E-03;
    COFD[2940] = -1.97810200E+01;
    COFD[2941] = 4.22062499E+00;
    COFD[2942] = -2.54326872E-01;
    COFD[2943] = 7.91017784E-03;
    COFD[2944] = -2.18355800E+01;
    COFD[2945] = 5.40298848E+00;
    COFD[2946] = -4.43954594E-01;
    COFD[2947] = 1.75542998E-02;
    COFD[2948] = -2.14831394E+01;
    COFD[2949] = 5.56978987E+00;
    COFD[2950] = -4.89141980E-01;
    COFD[2951] = 2.04499210E-02;
    COFD[2952] = -2.18478129E+01;
    COFD[2953] = 5.40298848E+00;
    COFD[2954] = -4.43954594E-01;
    COFD[2955] = 1.75542998E-02;
    COFD[2956] = -2.04146565E+01;
    COFD[2957] = 4.50250781E+00;
    COFD[2958] = -2.97622106E-01;
    COFD[2959] = 1.00481473E-02;
    COFD[2960] = -2.17456564E+01;
    COFD[2961] = 5.17945041E+00;
    COFD[2962] = -4.05514689E-01;
    COFD[2963] = 1.55141412E-02;
    COFD[2964] = -2.17984365E+01;
    COFD[2965] = 5.21869603E+00;
    COFD[2966] = -4.12084772E-01;
    COFD[2967] = 1.58573035E-02;
    COFD[2968] = -2.09559859E+01;
    COFD[2969] = 4.63716925E+00;
    COFD[2970] = -3.18815070E-01;
    COFD[2971] = 1.11115446E-02;
    COFD[2972] = -2.19969874E+01;
    COFD[2973] = 5.23595129E+00;
    COFD[2974] = -4.15079064E-01;
    COFD[2975] = 1.60168286E-02;
    COFD[2976] = -2.19969874E+01;
    COFD[2977] = 5.23595129E+00;
    COFD[2978] = -4.15079064E-01;
    COFD[2979] = 1.60168286E-02;
    COFD[2980] = -2.20030490E+01;
    COFD[2981] = 5.23595129E+00;
    COFD[2982] = -4.15079064E-01;
    COFD[2983] = 1.60168286E-02;
    COFD[2984] = -2.05192927E+01;
    COFD[2985] = 4.35981364E+00;
    COFD[2986] = -2.75585363E-01;
    COFD[2987] = 8.95565009E-03;
    COFD[2988] = -2.02311039E+01;
    COFD[2989] = 4.22278378E+00;
    COFD[2990] = -2.54653500E-01;
    COFD[2991] = 7.92616085E-03;
    COFD[2992] = -1.83039618E+01;
    COFD[2993] = 4.47952077E+00;
    COFD[2994] = -3.66569471E-01;
    COFD[2995] = 1.58916129E-02;
    COFD[2996] = -1.72286007E+01;
    COFD[2997] = 4.24084025E+00;
    COFD[2998] = -3.37428619E-01;
    COFD[2999] = 1.47032793E-02;
    COFD[3000] = -1.79310765E+01;
    COFD[3001] = 4.98037650E+00;
    COFD[3002] = -4.26676911E-01;
    COFD[3003] = 1.83007231E-02;
    COFD[3004] = -1.72473011E+01;
    COFD[3005] = 4.24084025E+00;
    COFD[3006] = -3.37428619E-01;
    COFD[3007] = 1.47032793E-02;
    COFD[3008] = -1.39315266E+01;
    COFD[3009] = 3.30394764E+00;
    COFD[3010] = -2.17920112E-01;
    COFD[3011] = 9.60284243E-03;
    COFD[3012] = -1.86507213E+01;
    COFD[3013] = 4.60874797E+00;
    COFD[3014] = -3.82368716E-01;
    COFD[3015] = 1.65370164E-02;
    COFD[3016] = -2.09565916E+01;
    COFD[3017] = 5.18380539E+00;
    COFD[3018] = -4.06234719E-01;
    COFD[3019] = 1.55515345E-02;
    COFD[3020] = -1.86641962E+01;
    COFD[3021] = 4.60874797E+00;
    COFD[3022] = -3.82368716E-01;
    COFD[3023] = 1.65370164E-02;
    COFD[3024] = -1.86576191E+01;
    COFD[3025] = 4.60874797E+00;
    COFD[3026] = -3.82368716E-01;
    COFD[3027] = 1.65370164E-02;
    COFD[3028] = -1.95548230E+01;
    COFD[3029] = 4.97133070E+00;
    COFD[3030] = -4.25604177E-01;
    COFD[3031] = 1.82582594E-02;
    COFD[3032] = -2.16798265E+01;
    COFD[3033] = 5.36811769E+00;
    COFD[3034] = -4.37727086E-01;
    COFD[3035] = 1.72167686E-02;
    COFD[3036] = -1.95770968E+01;
    COFD[3037] = 4.97133070E+00;
    COFD[3038] = -4.25604177E-01;
    COFD[3039] = 1.82582594E-02;
    COFD[3040] = -1.83296965E+01;
    COFD[3041] = 4.48570999E+00;
    COFD[3042] = -3.67301524E-01;
    COFD[3043] = 1.59204254E-02;
    COFD[3044] = -2.18848136E+01;
    COFD[3045] = 5.51302074E+00;
    COFD[3046] = -4.65263979E-01;
    COFD[3047] = 1.87580679E-02;
    COFD[3048] = -2.15126310E+01;
    COFD[3049] = 5.48426911E+00;
    COFD[3050] = -4.80606512E-01;
    COFD[3051] = 2.01811046E-02;
    COFD[3052] = -1.95154079E+01;
    COFD[3053] = 4.94787350E+00;
    COFD[3054] = -4.22829292E-01;
    COFD[3055] = 1.81487163E-02;
    COFD[3056] = -2.11427744E+01;
    COFD[3057] = 5.43893233E+00;
    COFD[3058] = -4.75546039E-01;
    COFD[3059] = 1.99938690E-02;
    COFD[3060] = -2.15206146E+01;
    COFD[3061] = 5.48426911E+00;
    COFD[3062] = -4.80606512E-01;
    COFD[3063] = 2.01811046E-02;
    COFD[3064] = -2.13961414E+01;
    COFD[3065] = 5.46685775E+00;
    COFD[3066] = -4.78665416E-01;
    COFD[3067] = 2.01093915E-02;
    COFD[3068] = -2.16718247E+01;
    COFD[3069] = 5.36811769E+00;
    COFD[3070] = -4.37727086E-01;
    COFD[3071] = 1.72167686E-02;
    COFD[3072] = -2.18318278E+01;
    COFD[3073] = 5.40298848E+00;
    COFD[3074] = -4.43954594E-01;
    COFD[3075] = 1.75542998E-02;
    COFD[3076] = -2.18355800E+01;
    COFD[3077] = 5.40298848E+00;
    COFD[3078] = -4.43954594E-01;
    COFD[3079] = 1.75542998E-02;
    COFD[3080] = -2.15453676E+01;
    COFD[3081] = 5.55313619E+00;
    COFD[3082] = -4.87753729E-01;
    COFD[3083] = 2.04203421E-02;
    COFD[3084] = -1.96068586E+01;
    COFD[3085] = 5.02434088E+00;
    COFD[3086] = -4.31889635E-01;
    COFD[3087] = 1.85072024E-02;
    COFD[3088] = -2.15547727E+01;
    COFD[3089] = 5.55313619E+00;
    COFD[3090] = -4.87753729E-01;
    COFD[3091] = 2.04203421E-02;
    COFD[3092] = -2.20307777E+01;
    COFD[3093] = 5.49663315E+00;
    COFD[3094] = -4.61182837E-01;
    COFD[3095] = 1.85035558E-02;
    COFD[3096] = -2.20511271E+01;
    COFD[3097] = 5.60809037E+00;
    COFD[3098] = -4.89400803E-01;
    COFD[3099] = 2.02760802E-02;
    COFD[3100] = -2.20228343E+01;
    COFD[3101] = 5.61211028E+00;
    COFD[3102] = -4.90893171E-01;
    COFD[3103] = 2.03793118E-02;
    COFD[3104] = -2.22964873E+01;
    COFD[3105] = 5.52353335E+00;
    COFD[3106] = -4.67314589E-01;
    COFD[3107] = 1.88744228E-02;
    COFD[3108] = -2.21795362E+01;
    COFD[3109] = 5.61233637E+00;
    COFD[3110] = -4.91419253E-01;
    COFD[3111] = 2.04216738E-02;
    COFD[3112] = -2.21795362E+01;
    COFD[3113] = 5.61233637E+00;
    COFD[3114] = -4.91419253E-01;
    COFD[3115] = 2.04216738E-02;
    COFD[3116] = -2.21838598E+01;
    COFD[3117] = 5.61233637E+00;
    COFD[3118] = -4.91419253E-01;
    COFD[3119] = 2.04216738E-02;
    COFD[3120] = -2.23166846E+01;
    COFD[3121] = 5.44889462E+00;
    COFD[3122] = -4.52325002E-01;
    COFD[3123] = 1.80132629E-02;
    COFD[3124] = -2.22462130E+01;
    COFD[3125] = 5.40356304E+00;
    COFD[3126] = -4.44060256E-01;
    COFD[3127] = 1.75601121E-02;
    COFD[3128] = -1.59884305E+01;
    COFD[3129] = 3.72220402E+00;
    COFD[3130] = -2.71150591E-01;
    COFD[3131] = 1.18665265E-02;
    COFD[3132] = -1.49500357E+01;
    COFD[3133] = 3.52327209E+00;
    COFD[3134] = -2.46286208E-01;
    COFD[3135] = 1.08285963E-02;
    COFD[3136] = -1.54460820E+01;
    COFD[3137] = 4.26819983E+00;
    COFD[3138] = -3.40766379E-01;
    COFD[3139] = 1.48393361E-02;
    COFD[3140] = -1.49718233E+01;
    COFD[3141] = 3.52327209E+00;
    COFD[3142] = -2.46286208E-01;
    COFD[3143] = 1.08285963E-02;
    COFD[3144] = -1.22004324E+01;
    COFD[3145] = 2.80725489E+00;
    COFD[3146] = -1.54291406E-01;
    COFD[3147] = 6.88290911E-03;
    COFD[3148] = -1.64169433E+01;
    COFD[3149] = 3.89309916E+00;
    COFD[3150] = -2.93528188E-01;
    COFD[3151] = 1.28463177E-02;
    COFD[3152] = -2.10440675E+01;
    COFD[3153] = 5.59806282E+00;
    COFD[3154] = -4.87109535E-01;
    COFD[3155] = 2.01370226E-02;
    COFD[3156] = -1.64338757E+01;
    COFD[3157] = 3.89309916E+00;
    COFD[3158] = -2.93528188E-01;
    COFD[3159] = 1.28463177E-02;
    COFD[3160] = -1.64255964E+01;
    COFD[3161] = 3.89309916E+00;
    COFD[3162] = -2.93528188E-01;
    COFD[3163] = 1.28463177E-02;
    COFD[3164] = -1.72572042E+01;
    COFD[3165] = 4.26063341E+00;
    COFD[3166] = -3.39848064E-01;
    COFD[3167] = 1.48021313E-02;
    COFD[3168] = -2.14303479E+01;
    COFD[3169] = 5.59268435E+00;
    COFD[3170] = -4.91232974E-01;
    COFD[3171] = 2.05064746E-02;
    COFD[3172] = -1.72828302E+01;
    COFD[3173] = 4.26063341E+00;
    COFD[3174] = -3.39848064E-01;
    COFD[3175] = 1.48021313E-02;
    COFD[3176] = -1.60261675E+01;
    COFD[3177] = 3.73312045E+00;
    COFD[3178] = -2.72579779E-01;
    COFD[3179] = 1.19290272E-02;
    COFD[3180] = -2.09241647E+01;
    COFD[3181] = 5.42316225E+00;
    COFD[3182] = -4.73702801E-01;
    COFD[3183] = 1.99217718E-02;
    COFD[3184] = -1.95737308E+01;
    COFD[3185] = 4.93449043E+00;
    COFD[3186] = -4.21243802E-01;
    COFD[3187] = 1.80859966E-02;
    COFD[3188] = -1.72316148E+01;
    COFD[3189] = 4.24011069E+00;
    COFD[3190] = -3.37339810E-01;
    COFD[3191] = 1.46996679E-02;
    COFD[3192] = -1.91367023E+01;
    COFD[3193] = 4.87703209E+00;
    COFD[3194] = -4.14222202E-01;
    COFD[3195] = 1.77987878E-02;
    COFD[3196] = -1.95836394E+01;
    COFD[3197] = 4.93449043E+00;
    COFD[3198] = -4.21243802E-01;
    COFD[3199] = 1.80859966E-02;
    COFD[3200] = -1.94485982E+01;
    COFD[3201] = 4.91446566E+00;
    COFD[3202] = -4.18837152E-01;
    COFD[3203] = 1.79893537E-02;
    COFD[3204] = -2.14204185E+01;
    COFD[3205] = 5.59268435E+00;
    COFD[3206] = -4.91232974E-01;
    COFD[3207] = 2.05064746E-02;
    COFD[3208] = -2.14782277E+01;
    COFD[3209] = 5.56978987E+00;
    COFD[3210] = -4.89141980E-01;
    COFD[3211] = 2.04499210E-02;
    COFD[3212] = -2.14831394E+01;
    COFD[3213] = 5.56978987E+00;
    COFD[3214] = -4.89141980E-01;
    COFD[3215] = 2.04499210E-02;
    COFD[3216] = -1.96068586E+01;
    COFD[3217] = 5.02434088E+00;
    COFD[3218] = -4.31889635E-01;
    COFD[3219] = 1.85072024E-02;
    COFD[3220] = -1.72414862E+01;
    COFD[3221] = 4.29808578E+00;
    COFD[3222] = -3.44235570E-01;
    COFD[3223] = 1.49727727E-02;
    COFD[3224] = -1.96183903E+01;
    COFD[3225] = 5.02434088E+00;
    COFD[3226] = -4.31889635E-01;
    COFD[3227] = 1.85072024E-02;
    COFD[3228] = -2.12680082E+01;
    COFD[3229] = 5.47935225E+00;
    COFD[3230] = -4.80056796E-01;
    COFD[3231] = 2.01607180E-02;
    COFD[3232] = -2.04251023E+01;
    COFD[3233] = 5.19993608E+00;
    COFD[3234] = -4.51334924E-01;
    COFD[3235] = 1.92158646E-02;
    COFD[3236] = -2.03116266E+01;
    COFD[3237] = 5.16758304E+00;
    COFD[3238] = -4.47606321E-01;
    COFD[3239] = 1.90728318E-02;
    COFD[3240] = -2.13858410E+01;
    COFD[3241] = 5.41773320E+00;
    COFD[3242] = -4.73414274E-01;
    COFD[3243] = 1.99258733E-02;
    COFD[3244] = -2.04750337E+01;
    COFD[3245] = 5.15745622E+00;
    COFD[3246] = -4.46648283E-01;
    COFD[3247] = 1.90458987E-02;
    COFD[3248] = -2.04750337E+01;
    COFD[3249] = 5.15745622E+00;
    COFD[3250] = -4.46648283E-01;
    COFD[3251] = 1.90458987E-02;
    COFD[3252] = -2.04806396E+01;
    COFD[3253] = 5.15745622E+00;
    COFD[3254] = -4.46648283E-01;
    COFD[3255] = 1.90458987E-02;
    COFD[3256] = -2.18644512E+01;
    COFD[3257] = 5.53438558E+00;
    COFD[3258] = -4.85873757E-01;
    COFD[3259] = 2.03606278E-02;
    COFD[3260] = -2.19764159E+01;
    COFD[3261] = 5.56943713E+00;
    COFD[3262] = -4.89114655E-01;
    COFD[3263] = 2.04494661E-02;
    COFD[3264] = -1.83137139E+01;
    COFD[3265] = 4.47952077E+00;
    COFD[3266] = -3.66569471E-01;
    COFD[3267] = 1.58916129E-02;
    COFD[3268] = -1.72357436E+01;
    COFD[3269] = 4.24084025E+00;
    COFD[3270] = -3.37428619E-01;
    COFD[3271] = 1.47032793E-02;
    COFD[3272] = -1.79317714E+01;
    COFD[3273] = 4.98037650E+00;
    COFD[3274] = -4.26676911E-01;
    COFD[3275] = 1.83007231E-02;
    COFD[3276] = -1.72547182E+01;
    COFD[3277] = 4.24084025E+00;
    COFD[3278] = -3.37428619E-01;
    COFD[3279] = 1.47032793E-02;
    COFD[3280] = -1.39328674E+01;
    COFD[3281] = 3.30394764E+00;
    COFD[3282] = -2.17920112E-01;
    COFD[3283] = 9.60284243E-03;
    COFD[3284] = -1.86611023E+01;
    COFD[3285] = 4.60874797E+00;
    COFD[3286] = -3.82368716E-01;
    COFD[3287] = 1.65370164E-02;
    COFD[3288] = -2.09642705E+01;
    COFD[3289] = 5.18380539E+00;
    COFD[3290] = -4.06234719E-01;
    COFD[3291] = 1.55515345E-02;
    COFD[3292] = -1.86748638E+01;
    COFD[3293] = 4.60874797E+00;
    COFD[3294] = -3.82368716E-01;
    COFD[3295] = 1.65370164E-02;
    COFD[3296] = -1.86681459E+01;
    COFD[3297] = 4.60874797E+00;
    COFD[3298] = -3.82368716E-01;
    COFD[3299] = 1.65370164E-02;
    COFD[3300] = -1.95613899E+01;
    COFD[3301] = 4.97133070E+00;
    COFD[3302] = -4.25604177E-01;
    COFD[3303] = 1.82582594E-02;
    COFD[3304] = -2.16899073E+01;
    COFD[3305] = 5.36811769E+00;
    COFD[3306] = -4.37727086E-01;
    COFD[3307] = 1.72167686E-02;
    COFD[3308] = -1.95839648E+01;
    COFD[3309] = 4.97133070E+00;
    COFD[3310] = -4.25604177E-01;
    COFD[3311] = 1.82582594E-02;
    COFD[3312] = -1.83394481E+01;
    COFD[3313] = 4.48570999E+00;
    COFD[3314] = -3.67301524E-01;
    COFD[3315] = 1.59204254E-02;
    COFD[3316] = -2.18950505E+01;
    COFD[3317] = 5.51302074E+00;
    COFD[3318] = -4.65263979E-01;
    COFD[3319] = 1.87580679E-02;
    COFD[3320] = -2.15225573E+01;
    COFD[3321] = 5.48426911E+00;
    COFD[3322] = -4.80606512E-01;
    COFD[3323] = 2.01811046E-02;
    COFD[3324] = -1.95225629E+01;
    COFD[3325] = 4.94787350E+00;
    COFD[3326] = -4.22829292E-01;
    COFD[3327] = 1.81487163E-02;
    COFD[3328] = -2.11525334E+01;
    COFD[3329] = 5.43893233E+00;
    COFD[3330] = -4.75546039E-01;
    COFD[3331] = 1.99938690E-02;
    COFD[3332] = -2.15307023E+01;
    COFD[3333] = 5.48426911E+00;
    COFD[3334] = -4.80606512E-01;
    COFD[3335] = 2.01811046E-02;
    COFD[3336] = -2.14079882E+01;
    COFD[3337] = 5.46685775E+00;
    COFD[3338] = -4.78665416E-01;
    COFD[3339] = 2.01093915E-02;
    COFD[3340] = -2.16817439E+01;
    COFD[3341] = 5.36811769E+00;
    COFD[3342] = -4.37727086E-01;
    COFD[3343] = 1.72167686E-02;
    COFD[3344] = -2.18439681E+01;
    COFD[3345] = 5.40298848E+00;
    COFD[3346] = -4.43954594E-01;
    COFD[3347] = 1.75542998E-02;
    COFD[3348] = -2.18478129E+01;
    COFD[3349] = 5.40298848E+00;
    COFD[3350] = -4.43954594E-01;
    COFD[3351] = 1.75542998E-02;
    COFD[3352] = -2.15547727E+01;
    COFD[3353] = 5.55313619E+00;
    COFD[3354] = -4.87753729E-01;
    COFD[3355] = 2.04203421E-02;
    COFD[3356] = -1.96183903E+01;
    COFD[3357] = 5.02434088E+00;
    COFD[3358] = -4.31889635E-01;
    COFD[3359] = 1.85072024E-02;
    COFD[3360] = -2.15643580E+01;
    COFD[3361] = 5.55313619E+00;
    COFD[3362] = -4.87753729E-01;
    COFD[3363] = 2.04203421E-02;
    COFD[3364] = -2.20425255E+01;
    COFD[3365] = 5.49663315E+00;
    COFD[3366] = -4.61182837E-01;
    COFD[3367] = 1.85035558E-02;
    COFD[3368] = -2.20626636E+01;
    COFD[3369] = 5.60809037E+00;
    COFD[3370] = -4.89400803E-01;
    COFD[3371] = 2.02760802E-02;
    COFD[3372] = -2.20344803E+01;
    COFD[3373] = 5.61211028E+00;
    COFD[3374] = -4.90893171E-01;
    COFD[3375] = 2.03793118E-02;
    COFD[3376] = -2.23094501E+01;
    COFD[3377] = 5.52353335E+00;
    COFD[3378] = -4.67314589E-01;
    COFD[3379] = 1.88744228E-02;
    COFD[3380] = -2.21912885E+01;
    COFD[3381] = 5.61233637E+00;
    COFD[3382] = -4.91419253E-01;
    COFD[3383] = 2.04216738E-02;
    COFD[3384] = -2.21912885E+01;
    COFD[3385] = 5.61233637E+00;
    COFD[3386] = -4.91419253E-01;
    COFD[3387] = 2.04216738E-02;
    COFD[3388] = -2.21957154E+01;
    COFD[3389] = 5.61233637E+00;
    COFD[3390] = -4.91419253E-01;
    COFD[3391] = 2.04216738E-02;
    COFD[3392] = -2.23307159E+01;
    COFD[3393] = 5.44889462E+00;
    COFD[3394] = -4.52325002E-01;
    COFD[3395] = 1.80132629E-02;
    COFD[3396] = -2.22602443E+01;
    COFD[3397] = 5.40356304E+00;
    COFD[3398] = -4.44060256E-01;
    COFD[3399] = 1.75601121E-02;
    COFD[3400] = -2.02693653E+01;
    COFD[3401] = 5.10426133E+00;
    COFD[3402] = -4.41256919E-01;
    COFD[3403] = 1.88737290E-02;
    COFD[3404] = -1.90915649E+01;
    COFD[3405] = 4.84384483E+00;
    COFD[3406] = -4.10265575E-01;
    COFD[3407] = 1.76414287E-02;
    COFD[3408] = -1.94691430E+01;
    COFD[3409] = 5.43830787E+00;
    COFD[3410] = -4.75472880E-01;
    COFD[3411] = 1.99909996E-02;
    COFD[3412] = -1.91136491E+01;
    COFD[3413] = 4.84384483E+00;
    COFD[3414] = -4.10265575E-01;
    COFD[3415] = 1.76414287E-02;
    COFD[3416] = -1.57040212E+01;
    COFD[3417] = 3.93614244E+00;
    COFD[3418] = -2.99111497E-01;
    COFD[3419] = 1.30888229E-02;
    COFD[3420] = -2.05235731E+01;
    COFD[3421] = 5.18417470E+00;
    COFD[3422] = -4.49491573E-01;
    COFD[3423] = 1.91438508E-02;
    COFD[3424] = -1.87419199E+01;
    COFD[3425] = 3.96926341E+00;
    COFD[3426] = -2.16412264E-01;
    COFD[3427] = 6.06012078E-03;
    COFD[3428] = -2.05408665E+01;
    COFD[3429] = 5.18417470E+00;
    COFD[3430] = -4.49491573E-01;
    COFD[3431] = 1.91438508E-02;
    COFD[3432] = -2.05324091E+01;
    COFD[3433] = 5.18417470E+00;
    COFD[3434] = -4.49491573E-01;
    COFD[3435] = 1.91438508E-02;
    COFD[3436] = -2.11378465E+01;
    COFD[3437] = 5.42846112E+00;
    COFD[3438] = -4.74321870E-01;
    COFD[3439] = 1.99459749E-02;
    COFD[3440] = -2.01064363E+01;
    COFD[3441] = 4.41511629E+00;
    COFD[3442] = -2.84086963E-01;
    COFD[3443] = 9.37586971E-03;
    COFD[3444] = -2.11637902E+01;
    COFD[3445] = 5.42846112E+00;
    COFD[3446] = -4.74321870E-01;
    COFD[3447] = 1.99459749E-02;
    COFD[3448] = -2.02969740E+01;
    COFD[3449] = 5.11106992E+00;
    COFD[3450] = -4.42047129E-01;
    COFD[3451] = 1.89042990E-02;
    COFD[3452] = -2.09272429E+01;
    COFD[3453] = 4.82184721E+00;
    COFD[3454] = -3.48128875E-01;
    COFD[3455] = 1.25918978E-02;
    COFD[3456] = -2.22377812E+01;
    COFD[3457] = 5.53139819E+00;
    COFD[3458] = -4.68828555E-01;
    COFD[3459] = 1.89597887E-02;
    COFD[3460] = -2.11341653E+01;
    COFD[3461] = 5.41773516E+00;
    COFD[3462] = -4.73414338E-01;
    COFD[3463] = 1.99258685E-02;
    COFD[3464] = -2.20445411E+01;
    COFD[3465] = 5.56049839E+00;
    COFD[3466] = -4.74367872E-01;
    COFD[3467] = 1.92702787E-02;
    COFD[3468] = -2.22478879E+01;
    COFD[3469] = 5.53139819E+00;
    COFD[3470] = -4.68828555E-01;
    COFD[3471] = 1.89597887E-02;
    COFD[3472] = -2.22176950E+01;
    COFD[3473] = 5.54251230E+00;
    COFD[3474] = -4.70946314E-01;
    COFD[3475] = 1.90785869E-02;
    COFD[3476] = -2.00963085E+01;
    COFD[3477] = 4.41511629E+00;
    COFD[3478] = -2.84086963E-01;
    COFD[3479] = 9.37586971E-03;
    COFD[3480] = -2.04096182E+01;
    COFD[3481] = 4.50250781E+00;
    COFD[3482] = -2.97622106E-01;
    COFD[3483] = 1.00481473E-02;
    COFD[3484] = -2.04146565E+01;
    COFD[3485] = 4.50250781E+00;
    COFD[3486] = -2.97622106E-01;
    COFD[3487] = 1.00481473E-02;
    COFD[3488] = -2.20307777E+01;
    COFD[3489] = 5.49663315E+00;
    COFD[3490] = -4.61182837E-01;
    COFD[3491] = 1.85035558E-02;
    COFD[3492] = -2.12680082E+01;
    COFD[3493] = 5.47935225E+00;
    COFD[3494] = -4.80056796E-01;
    COFD[3495] = 2.01607180E-02;
    COFD[3496] = -2.20425255E+01;
    COFD[3497] = 5.49663315E+00;
    COFD[3498] = -4.61182837E-01;
    COFD[3499] = 1.85035558E-02;
    COFD[3500] = -2.09121217E+01;
    COFD[3501] = 4.72895031E+00;
    COFD[3502] = -3.33332771E-01;
    COFD[3503] = 1.18431478E-02;
    COFD[3504] = -2.20250551E+01;
    COFD[3505] = 5.31412694E+00;
    COFD[3506] = -4.28473898E-01;
    COFD[3507] = 1.67264841E-02;
    COFD[3508] = -2.20656222E+01;
    COFD[3509] = 5.34774760E+00;
    COFD[3510] = -4.34239753E-01;
    COFD[3511] = 1.70320676E-02;
    COFD[3512] = -2.14189878E+01;
    COFD[3513] = 4.85433278E+00;
    COFD[3514] = -3.53258974E-01;
    COFD[3515] = 1.28503760E-02;
    COFD[3516] = -2.22604974E+01;
    COFD[3517] = 5.36643605E+00;
    COFD[3518] = -4.37440735E-01;
    COFD[3519] = 1.72016388E-02;
    COFD[3520] = -2.22604974E+01;
    COFD[3521] = 5.36643605E+00;
    COFD[3522] = -4.37440735E-01;
    COFD[3523] = 1.72016388E-02;
    COFD[3524] = -2.22662419E+01;
    COFD[3525] = 5.36643605E+00;
    COFD[3526] = -4.37440735E-01;
    COFD[3527] = 1.72016388E-02;
    COFD[3528] = -2.10617067E+01;
    COFD[3529] = 4.61245983E+00;
    COFD[3530] = -3.14873896E-01;
    COFD[3531] = 1.09118770E-02;
    COFD[3532] = -2.08429322E+01;
    COFD[3533] = 4.50409026E+00;
    COFD[3534] = -2.97868419E-01;
    COFD[3535] = 1.00604224E-02;
    COFD[3536] = -1.91796663E+01;
    COFD[3537] = 4.70714822E+00;
    COFD[3538] = -3.94261134E-01;
    COFD[3539] = 1.70175169E-02;
    COFD[3540] = -1.80480958E+01;
    COFD[3541] = 4.45434023E+00;
    COFD[3542] = -3.63584633E-01;
    COFD[3543] = 1.57739270E-02;
    COFD[3544] = -1.86493112E+01;
    COFD[3545] = 5.16040659E+00;
    COFD[3546] = -4.46843492E-01;
    COFD[3547] = 1.90466181E-02;
    COFD[3548] = -1.80698901E+01;
    COFD[3549] = 4.45434023E+00;
    COFD[3550] = -3.63584633E-01;
    COFD[3551] = 1.57739270E-02;
    COFD[3552] = -1.46719197E+01;
    COFD[3553] = 3.52400594E+00;
    COFD[3554] = -2.46379985E-01;
    COFD[3555] = 1.08326032E-02;
    COFD[3556] = -1.94912151E+01;
    COFD[3557] = 4.81575071E+00;
    COFD[3558] = -4.07042139E-01;
    COFD[3559] = 1.75187504E-02;
    COFD[3560] = -2.05372411E+01;
    COFD[3561] = 4.83379373E+00;
    COFD[3562] = -3.50008083E-01;
    COFD[3563] = 1.26863426E-02;
    COFD[3564] = -1.95081555E+01;
    COFD[3565] = 4.81575071E+00;
    COFD[3566] = -4.07042139E-01;
    COFD[3567] = 1.75187504E-02;
    COFD[3568] = -1.94998722E+01;
    COFD[3569] = 4.81575071E+00;
    COFD[3570] = -4.07042139E-01;
    COFD[3571] = 1.75187504E-02;
    COFD[3572] = -2.03111230E+01;
    COFD[3573] = 5.15740122E+00;
    COFD[3574] = -4.46644818E-01;
    COFD[3575] = 1.90459001E-02;
    COFD[3576] = -2.15258568E+01;
    COFD[3577] = 5.12799307E+00;
    COFD[3578] = -3.96938732E-01;
    COFD[3579] = 1.50673195E-02;
    COFD[3580] = -2.03367561E+01;
    COFD[3581] = 5.15740122E+00;
    COFD[3582] = -4.46644818E-01;
    COFD[3583] = 1.90459001E-02;
    COFD[3584] = -1.92044492E+01;
    COFD[3585] = 4.71304783E+00;
    COFD[3586] = -3.94942083E-01;
    COFD[3587] = 1.70435959E-02;
    COFD[3588] = -2.19617977E+01;
    COFD[3589] = 5.37170913E+00;
    COFD[3590] = -4.38338667E-01;
    COFD[3591] = 1.72490835E-02;
    COFD[3592] = -2.21694197E+01;
    COFD[3593] = 5.60403905E+00;
    COFD[3594] = -4.91221691E-01;
    COFD[3595] = 2.04473483E-02;
    COFD[3596] = -2.03123540E+01;
    COFD[3597] = 5.14854169E+00;
    COFD[3598] = -4.45984343E-01;
    COFD[3599] = 1.90374217E-02;
    COFD[3600] = -2.19032561E+01;
    COFD[3601] = 5.59794138E+00;
    COFD[3602] = -4.91684532E-01;
    COFD[3603] = 2.05170953E-02;
    COFD[3604] = -2.21793326E+01;
    COFD[3605] = 5.60403905E+00;
    COFD[3606] = -4.91221691E-01;
    COFD[3607] = 2.04473483E-02;
    COFD[3608] = -2.21216828E+01;
    COFD[3609] = 5.60203389E+00;
    COFD[3610] = -4.91444416E-01;
    COFD[3611] = 2.04761886E-02;
    COFD[3612] = -2.15159231E+01;
    COFD[3613] = 5.12799307E+00;
    COFD[3614] = -3.96938732E-01;
    COFD[3615] = 1.50673195E-02;
    COFD[3616] = -2.17407419E+01;
    COFD[3617] = 5.17945041E+00;
    COFD[3618] = -4.05514689E-01;
    COFD[3619] = 1.55141412E-02;
    COFD[3620] = -2.17456564E+01;
    COFD[3621] = 5.17945041E+00;
    COFD[3622] = -4.05514689E-01;
    COFD[3623] = 1.55141412E-02;
    COFD[3624] = -2.20511271E+01;
    COFD[3625] = 5.60809037E+00;
    COFD[3626] = -4.89400803E-01;
    COFD[3627] = 2.02760802E-02;
    COFD[3628] = -2.04251023E+01;
    COFD[3629] = 5.19993608E+00;
    COFD[3630] = -4.51334924E-01;
    COFD[3631] = 1.92158646E-02;
    COFD[3632] = -2.20626636E+01;
    COFD[3633] = 5.60809037E+00;
    COFD[3634] = -4.89400803E-01;
    COFD[3635] = 2.02760802E-02;
    COFD[3636] = -2.20250551E+01;
    COFD[3637] = 5.31412694E+00;
    COFD[3638] = -4.28473898E-01;
    COFD[3639] = 1.67264841E-02;
    COFD[3640] = -2.24021886E+01;
    COFD[3641] = 5.58364149E+00;
    COFD[3642] = -4.79184111E-01;
    COFD[3643] = 1.95516164E-02;
    COFD[3644] = -2.23689627E+01;
    COFD[3645] = 5.58513878E+00;
    COFD[3646] = -4.80389524E-01;
    COFD[3647] = 1.96438689E-02;
    COFD[3648] = -2.24063517E+01;
    COFD[3649] = 5.38918452E+00;
    COFD[3650] = -4.41385864E-01;
    COFD[3651] = 1.74122637E-02;
    COFD[3652] = -2.25161211E+01;
    COFD[3653] = 5.58521783E+00;
    COFD[3654] = -4.80947522E-01;
    COFD[3655] = 1.96897222E-02;
    COFD[3656] = -2.25161211E+01;
    COFD[3657] = 5.58521783E+00;
    COFD[3658] = -4.80947522E-01;
    COFD[3659] = 1.96897222E-02;
    COFD[3660] = -2.25217301E+01;
    COFD[3661] = 5.58521783E+00;
    COFD[3662] = -4.80947522E-01;
    COFD[3663] = 1.96897222E-02;
    COFD[3664] = -2.22719261E+01;
    COFD[3665] = 5.24354480E+00;
    COFD[3666] = -4.16412182E-01;
    COFD[3667] = 1.60883090E-02;
    COFD[3668] = -2.21603646E+01;
    COFD[3669] = 5.18050127E+00;
    COFD[3670] = -4.05688517E-01;
    COFD[3671] = 1.55231713E-02;
    COFD[3672] = -1.90859283E+01;
    COFD[3673] = 4.68079396E+00;
    COFD[3674] = -3.91231550E-01;
    COFD[3675] = 1.69021170E-02;
    COFD[3676] = -1.79361160E+01;
    COFD[3677] = 4.42139452E+00;
    COFD[3678] = -3.59567329E-01;
    COFD[3679] = 1.56103969E-02;
    COFD[3680] = -1.85748546E+01;
    COFD[3681] = 5.14789919E+00;
    COFD[3682] = -4.45930850E-01;
    COFD[3683] = 1.90363341E-02;
    COFD[3684] = -1.79580609E+01;
    COFD[3685] = 4.42139452E+00;
    COFD[3686] = -3.59567329E-01;
    COFD[3687] = 1.56103969E-02;
    COFD[3688] = -1.45715797E+01;
    COFD[3689] = 3.49477850E+00;
    COFD[3690] = -2.42635772E-01;
    COFD[3691] = 1.06721490E-02;
    COFD[3692] = -1.93917298E+01;
    COFD[3693] = 4.78708023E+00;
    COFD[3694] = -4.03693144E-01;
    COFD[3695] = 1.73884817E-02;
    COFD[3696] = -2.06310304E+01;
    COFD[3697] = 4.89289496E+00;
    COFD[3698] = -3.59346263E-01;
    COFD[3699] = 1.31570901E-02;
    COFD[3700] = -1.94088529E+01;
    COFD[3701] = 4.78708023E+00;
    COFD[3702] = -4.03693144E-01;
    COFD[3703] = 1.73884817E-02;
    COFD[3704] = -1.94004795E+01;
    COFD[3705] = 4.78708023E+00;
    COFD[3706] = -4.03693144E-01;
    COFD[3707] = 1.73884817E-02;
    COFD[3708] = -2.02434438E+01;
    COFD[3709] = 5.14418672E+00;
    COFD[3710] = -4.45631004E-01;
    COFD[3711] = 1.90308403E-02;
    COFD[3712] = -2.15802788E+01;
    COFD[3713] = 5.16868516E+00;
    COFD[3714] = -4.03721581E-01;
    COFD[3715] = 1.54206640E-02;
    COFD[3716] = -2.02692384E+01;
    COFD[3717] = 5.14418672E+00;
    COFD[3718] = -4.45631004E-01;
    COFD[3719] = 1.90308403E-02;
    COFD[3720] = -1.91118445E+01;
    COFD[3721] = 4.68715685E+00;
    COFD[3722] = -3.91979493E-01;
    COFD[3723] = 1.69314004E-02;
    COFD[3724] = -2.19873532E+01;
    COFD[3725] = 5.39977369E+00;
    COFD[3726] = -4.43340854E-01;
    COFD[3727] = 1.75199613E-02;
    COFD[3728] = -2.21242889E+01;
    COFD[3729] = 5.60010742E+00;
    COFD[3730] = -4.91597429E-01;
    COFD[3731] = 2.04987718E-02;
    COFD[3732] = -2.02318658E+01;
    COFD[3733] = 5.12963391E+00;
    COFD[3734] = -4.44146826E-01;
    COFD[3735] = 1.89829640E-02;
    COFD[3736] = -2.18222696E+01;
    COFD[3737] = 5.57940140E+00;
    COFD[3738] = -4.89964112E-01;
    COFD[3739] = 2.04689539E-02;
    COFD[3740] = -2.21343023E+01;
    COFD[3741] = 5.60010742E+00;
    COFD[3742] = -4.91597429E-01;
    COFD[3743] = 2.04987718E-02;
    COFD[3744] = -2.20725883E+01;
    COFD[3745] = 5.59642965E+00;
    COFD[3746] = -4.91577716E-01;
    COFD[3747] = 2.05159582E-02;
    COFD[3748] = -2.15702446E+01;
    COFD[3749] = 5.16868516E+00;
    COFD[3750] = -4.03721581E-01;
    COFD[3751] = 1.54206640E-02;
    COFD[3752] = -2.17934580E+01;
    COFD[3753] = 5.21869603E+00;
    COFD[3754] = -4.12084772E-01;
    COFD[3755] = 1.58573035E-02;
    COFD[3756] = -2.17984365E+01;
    COFD[3757] = 5.21869603E+00;
    COFD[3758] = -4.12084772E-01;
    COFD[3759] = 1.58573035E-02;
    COFD[3760] = -2.20228343E+01;
    COFD[3761] = 5.61211028E+00;
    COFD[3762] = -4.90893171E-01;
    COFD[3763] = 2.03793118E-02;
    COFD[3764] = -2.03116266E+01;
    COFD[3765] = 5.16758304E+00;
    COFD[3766] = -4.47606321E-01;
    COFD[3767] = 1.90728318E-02;
    COFD[3768] = -2.20344803E+01;
    COFD[3769] = 5.61211028E+00;
    COFD[3770] = -4.90893171E-01;
    COFD[3771] = 2.03793118E-02;
    COFD[3772] = -2.20656222E+01;
    COFD[3773] = 5.34774760E+00;
    COFD[3774] = -4.34239753E-01;
    COFD[3775] = 1.70320676E-02;
    COFD[3776] = -2.23689627E+01;
    COFD[3777] = 5.58513878E+00;
    COFD[3778] = -4.80389524E-01;
    COFD[3779] = 1.96438689E-02;
    COFD[3780] = -2.23318349E+01;
    COFD[3781] = 5.58508387E+00;
    COFD[3782] = -4.81385216E-01;
    COFD[3783] = 1.97267369E-02;
    COFD[3784] = -2.24242090E+01;
    COFD[3785] = 5.41204570E+00;
    COFD[3786] = -4.45647651E-01;
    COFD[3787] = 1.76482725E-02;
    COFD[3788] = -2.24797372E+01;
    COFD[3789] = 5.58492389E+00;
    COFD[3790] = -4.81921515E-01;
    COFD[3791] = 1.97721229E-02;
    COFD[3792] = -2.24797372E+01;
    COFD[3793] = 5.58492389E+00;
    COFD[3794] = -4.81921515E-01;
    COFD[3795] = 1.97721229E-02;
    COFD[3796] = -2.24854161E+01;
    COFD[3797] = 5.58492389E+00;
    COFD[3798] = -4.81921515E-01;
    COFD[3799] = 1.97721229E-02;
    COFD[3800] = -2.23086362E+01;
    COFD[3801] = 5.27409919E+00;
    COFD[3802] = -4.21767988E-01;
    COFD[3803] = 1.63754073E-02;
    COFD[3804] = -2.22169882E+01;
    COFD[3805] = 5.21950983E+00;
    COFD[3806] = -4.12223195E-01;
    COFD[3807] = 1.58645894E-02;
    COFD[3808] = -2.03249567E+01;
    COFD[3809] = 5.03292431E+00;
    COFD[3810] = -4.32887308E-01;
    COFD[3811] = 1.85458102E-02;
    COFD[3812] = -1.91805416E+01;
    COFD[3813] = 4.78152337E+00;
    COFD[3814] = -4.03052726E-01;
    COFD[3815] = 1.73639773E-02;
    COFD[3816] = -1.96115230E+01;
    COFD[3817] = 5.40502814E+00;
    COFD[3818] = -4.72746583E-01;
    COFD[3819] = 1.99363615E-02;
    COFD[3820] = -1.92042395E+01;
    COFD[3821] = 4.78152337E+00;
    COFD[3822] = -4.03052726E-01;
    COFD[3823] = 1.73639773E-02;
    COFD[3824] = -1.56932803E+01;
    COFD[3825] = 3.84112355E+00;
    COFD[3826] = -2.86739743E-01;
    COFD[3827] = 1.25496520E-02;
    COFD[3828] = -2.06513321E+01;
    COFD[3829] = 5.14159226E+00;
    COFD[3830] = -4.45395793E-01;
    COFD[3831] = 1.90248052E-02;
    COFD[3832] = -1.94510032E+01;
    COFD[3833] = 4.17438087E+00;
    COFD[3834] = -2.47331309E-01;
    COFD[3835] = 7.56793055E-03;
    COFD[3836] = -2.06706896E+01;
    COFD[3837] = 5.14159226E+00;
    COFD[3838] = -4.45395793E-01;
    COFD[3839] = 1.90248052E-02;
    COFD[3840] = -2.06612128E+01;
    COFD[3841] = 5.14159226E+00;
    COFD[3842] = -4.45395793E-01;
    COFD[3843] = 1.90248052E-02;
    COFD[3844] = -2.12684314E+01;
    COFD[3845] = 5.40206122E+00;
    COFD[3846] = -4.72555229E-01;
    COFD[3847] = 1.99358199E-02;
    COFD[3848] = -2.06753031E+01;
    COFD[3849] = 4.56779427E+00;
    COFD[3850] = -3.07785839E-01;
    COFD[3851] = 1.05545767E-02;
    COFD[3852] = -2.12960899E+01;
    COFD[3853] = 5.40206122E+00;
    COFD[3854] = -4.72555229E-01;
    COFD[3855] = 1.99358199E-02;
    COFD[3856] = -2.03525027E+01;
    COFD[3857] = 5.04005588E+00;
    COFD[3858] = -4.33725091E-01;
    COFD[3859] = 1.85786663E-02;
    COFD[3860] = -2.14535736E+01;
    COFD[3861] = 4.95733546E+00;
    COFD[3862] = -3.69505821E-01;
    COFD[3863] = 1.36687047E-02;
    COFD[3864] = -2.25375432E+01;
    COFD[3865] = 5.57713269E+00;
    COFD[3866] = -4.77555529E-01;
    COFD[3867] = 1.94497781E-02;
    COFD[3868] = -2.12758591E+01;
    COFD[3869] = 5.39400772E+00;
    COFD[3870] = -4.72026046E-01;
    COFD[3871] = 1.99336599E-02;
    COFD[3872] = -2.23144051E+01;
    COFD[3873] = 5.58484440E+00;
    COFD[3874] = -4.80068319E-01;
    COFD[3875] = 1.96187346E-02;
    COFD[3876] = -2.25487737E+01;
    COFD[3877] = 5.57713269E+00;
    COFD[3878] = -4.77555529E-01;
    COFD[3879] = 1.94497781E-02;
    COFD[3880] = -2.25320751E+01;
    COFD[3881] = 5.58240011E+00;
    COFD[3882] = -4.78844918E-01;
    COFD[3883] = 1.95298191E-02;
    COFD[3884] = -2.06640503E+01;
    COFD[3885] = 4.56779427E+00;
    COFD[3886] = -3.07785839E-01;
    COFD[3887] = 1.05545767E-02;
    COFD[3888] = -2.09502015E+01;
    COFD[3889] = 4.63716925E+00;
    COFD[3890] = -3.18815070E-01;
    COFD[3891] = 1.11115446E-02;
    COFD[3892] = -2.09559859E+01;
    COFD[3893] = 4.63716925E+00;
    COFD[3894] = -3.18815070E-01;
    COFD[3895] = 1.11115446E-02;
    COFD[3896] = -2.22964873E+01;
    COFD[3897] = 5.52353335E+00;
    COFD[3898] = -4.67314589E-01;
    COFD[3899] = 1.88744228E-02;
    COFD[3900] = -2.13858410E+01;
    COFD[3901] = 5.41773320E+00;
    COFD[3902] = -4.73414274E-01;
    COFD[3903] = 1.99258733E-02;
    COFD[3904] = -2.23094501E+01;
    COFD[3905] = 5.52353335E+00;
    COFD[3906] = -4.67314589E-01;
    COFD[3907] = 1.88744228E-02;
    COFD[3908] = -2.14189878E+01;
    COFD[3909] = 4.85433278E+00;
    COFD[3910] = -3.53258974E-01;
    COFD[3911] = 1.28503760E-02;
    COFD[3912] = -2.24063517E+01;
    COFD[3913] = 5.38918452E+00;
    COFD[3914] = -4.41385864E-01;
    COFD[3915] = 1.74122637E-02;
    COFD[3916] = -2.24242090E+01;
    COFD[3917] = 5.41204570E+00;
    COFD[3918] = -4.45647651E-01;
    COFD[3919] = 1.76482725E-02;
    COFD[3920] = -2.19461461E+01;
    COFD[3921] = 4.99238651E+00;
    COFD[3922] = -3.75026436E-01;
    COFD[3923] = 1.39466941E-02;
    COFD[3924] = -2.25933800E+01;
    COFD[3925] = 5.42562556E+00;
    COFD[3926] = -4.48132861E-01;
    COFD[3927] = 1.77847329E-02;
    COFD[3928] = -2.25933800E+01;
    COFD[3929] = 5.42562556E+00;
    COFD[3930] = -4.48132861E-01;
    COFD[3931] = 1.77847329E-02;
    COFD[3932] = -2.25999351E+01;
    COFD[3933] = 5.42562556E+00;
    COFD[3934] = -4.48132861E-01;
    COFD[3935] = 1.77847329E-02;
    COFD[3936] = -2.15763465E+01;
    COFD[3937] = 4.74313639E+00;
    COFD[3938] = -3.35584534E-01;
    COFD[3939] = 1.19568362E-02;
    COFD[3940] = -2.13682139E+01;
    COFD[3941] = 4.63852747E+00;
    COFD[3942] = -3.19031243E-01;
    COFD[3943] = 1.11224841E-02;
    COFD[3944] = -1.92062897E+01;
    COFD[3945] = 4.66318669E+00;
    COFD[3946] = -3.89108667E-01;
    COFD[3947] = 1.68165377E-02;
    COFD[3948] = -1.80724788E+01;
    COFD[3949] = 4.40247898E+00;
    COFD[3950] = -3.57238362E-01;
    COFD[3951] = 1.55145651E-02;
    COFD[3952] = -1.87481780E+01;
    COFD[3953] = 5.13858656E+00;
    COFD[3954] = -4.45075387E-01;
    COFD[3955] = 1.90137309E-02;
    COFD[3956] = -1.80945693E+01;
    COFD[3957] = 4.40247898E+00;
    COFD[3958] = -3.57238362E-01;
    COFD[3959] = 1.55145651E-02;
    COFD[3960] = -1.47137939E+01;
    COFD[3961] = 3.48023191E+00;
    COFD[3962] = -2.40800798E-01;
    COFD[3963] = 1.05947990E-02;
    COFD[3964] = -1.95201830E+01;
    COFD[3965] = 4.77151544E+00;
    COFD[3966] = -4.01882811E-01;
    COFD[3967] = 1.73184814E-02;
    COFD[3968] = -2.08879167E+01;
    COFD[3969] = 4.92602269E+00;
    COFD[3970] = -3.64572914E-01;
    COFD[3971] = 1.34203681E-02;
    COFD[3972] = -1.95374840E+01;
    COFD[3973] = 4.77151544E+00;
    COFD[3974] = -4.01882811E-01;
    COFD[3975] = 1.73184814E-02;
    COFD[3976] = -1.95290229E+01;
    COFD[3977] = 4.77151544E+00;
    COFD[3978] = -4.01882811E-01;
    COFD[3979] = 1.73184814E-02;
    COFD[3980] = -2.03711787E+01;
    COFD[3981] = 5.13279789E+00;
    COFD[3982] = -4.44474174E-01;
    COFD[3983] = 1.89937678E-02;
    COFD[3984] = -2.17926864E+01;
    COFD[3985] = 5.19232842E+00;
    COFD[3986] = -4.07643284E-01;
    COFD[3987] = 1.56246434E-02;
    COFD[3988] = -2.03971290E+01;
    COFD[3989] = 5.13279789E+00;
    COFD[3990] = -4.44474174E-01;
    COFD[3991] = 1.89937678E-02;
    COFD[3992] = -1.92334028E+01;
    COFD[3993] = 4.67033934E+00;
    COFD[3994] = -3.89971551E-01;
    COFD[3995] = 1.68513441E-02;
    COFD[3996] = -2.21713935E+01;
    COFD[3997] = 5.41196486E+00;
    COFD[3998] = -4.45632422E-01;
    COFD[3999] = 1.76474237E-02;
    COFD[4000] = -2.22600844E+01;
    COFD[4001] = 5.59632316E+00;
    COFD[4002] = -4.91568011E-01;
    COFD[4003] = 2.05156966E-02;
    COFD[4004] = -2.03526104E+01;
    COFD[4005] = 5.11453301E+00;
    COFD[4006] = -4.42447016E-01;
    COFD[4007] = 1.89196698E-02;
    COFD[4008] = -2.19617258E+01;
    COFD[4009] = 5.57026255E+00;
    COFD[4010] = -4.89178491E-01;
    COFD[4011] = 2.04505218E-02;
    COFD[4012] = -2.22701953E+01;
    COFD[4013] = 5.59632316E+00;
    COFD[4014] = -4.91568011E-01;
    COFD[4015] = 2.05156966E-02;
    COFD[4016] = -2.22052004E+01;
    COFD[4017] = 5.58604166E+00;
    COFD[4018] = -4.90602184E-01;
    COFD[4019] = 2.04880352E-02;
    COFD[4020] = -2.17825544E+01;
    COFD[4021] = 5.19232842E+00;
    COFD[4022] = -4.07643284E-01;
    COFD[4023] = 1.56246434E-02;
    COFD[4024] = -2.19919464E+01;
    COFD[4025] = 5.23595129E+00;
    COFD[4026] = -4.15079064E-01;
    COFD[4027] = 1.60168286E-02;
    COFD[4028] = -2.19969874E+01;
    COFD[4029] = 5.23595129E+00;
    COFD[4030] = -4.15079064E-01;
    COFD[4031] = 1.60168286E-02;
    COFD[4032] = -2.21795362E+01;
    COFD[4033] = 5.61233637E+00;
    COFD[4034] = -4.91419253E-01;
    COFD[4035] = 2.04216738E-02;
    COFD[4036] = -2.04750337E+01;
    COFD[4037] = 5.15745622E+00;
    COFD[4038] = -4.46648283E-01;
    COFD[4039] = 1.90458987E-02;
    COFD[4040] = -2.21912885E+01;
    COFD[4041] = 5.61233637E+00;
    COFD[4042] = -4.91419253E-01;
    COFD[4043] = 2.04216738E-02;
    COFD[4044] = -2.22604974E+01;
    COFD[4045] = 5.36643605E+00;
    COFD[4046] = -4.37440735E-01;
    COFD[4047] = 1.72016388E-02;
    COFD[4048] = -2.25161211E+01;
    COFD[4049] = 5.58521783E+00;
    COFD[4050] = -4.80947522E-01;
    COFD[4051] = 1.96897222E-02;
    COFD[4052] = -2.24797372E+01;
    COFD[4053] = 5.58492389E+00;
    COFD[4054] = -4.81921515E-01;
    COFD[4055] = 1.97721229E-02;
    COFD[4056] = -2.25933800E+01;
    COFD[4057] = 5.42562556E+00;
    COFD[4058] = -4.48132861E-01;
    COFD[4059] = 1.77847329E-02;
    COFD[4060] = -2.26149345E+01;
    COFD[4061] = 5.58414475E+00;
    COFD[4062] = -4.82375215E-01;
    COFD[4063] = 1.98138565E-02;
    COFD[4064] = -2.26149345E+01;
    COFD[4065] = 5.58414475E+00;
    COFD[4066] = -4.82375215E-01;
    COFD[4067] = 1.98138565E-02;
    COFD[4068] = -2.26206818E+01;
    COFD[4069] = 5.58414475E+00;
    COFD[4070] = -4.82375215E-01;
    COFD[4071] = 1.98138565E-02;
    COFD[4072] = -2.24927410E+01;
    COFD[4073] = 5.29529316E+00;
    COFD[4074] = -4.25345414E-01;
    COFD[4075] = 1.65635428E-02;
    COFD[4076] = -2.23925982E+01;
    COFD[4077] = 5.23666690E+00;
    COFD[4078] = -4.15204403E-01;
    COFD[4079] = 1.60235416E-02;
    COFD[4080] = -1.92062897E+01;
    COFD[4081] = 4.66318669E+00;
    COFD[4082] = -3.89108667E-01;
    COFD[4083] = 1.68165377E-02;
    COFD[4084] = -1.80724788E+01;
    COFD[4085] = 4.40247898E+00;
    COFD[4086] = -3.57238362E-01;
    COFD[4087] = 1.55145651E-02;
    COFD[4088] = -1.87481780E+01;
    COFD[4089] = 5.13858656E+00;
    COFD[4090] = -4.45075387E-01;
    COFD[4091] = 1.90137309E-02;
    COFD[4092] = -1.80945693E+01;
    COFD[4093] = 4.40247898E+00;
    COFD[4094] = -3.57238362E-01;
    COFD[4095] = 1.55145651E-02;
    COFD[4096] = -1.47137939E+01;
    COFD[4097] = 3.48023191E+00;
    COFD[4098] = -2.40800798E-01;
    COFD[4099] = 1.05947990E-02;
    COFD[4100] = -1.95201830E+01;
    COFD[4101] = 4.77151544E+00;
    COFD[4102] = -4.01882811E-01;
    COFD[4103] = 1.73184814E-02;
    COFD[4104] = -2.08879167E+01;
    COFD[4105] = 4.92602269E+00;
    COFD[4106] = -3.64572914E-01;
    COFD[4107] = 1.34203681E-02;
    COFD[4108] = -1.95374840E+01;
    COFD[4109] = 4.77151544E+00;
    COFD[4110] = -4.01882811E-01;
    COFD[4111] = 1.73184814E-02;
    COFD[4112] = -1.95290229E+01;
    COFD[4113] = 4.77151544E+00;
    COFD[4114] = -4.01882811E-01;
    COFD[4115] = 1.73184814E-02;
    COFD[4116] = -2.03711787E+01;
    COFD[4117] = 5.13279789E+00;
    COFD[4118] = -4.44474174E-01;
    COFD[4119] = 1.89937678E-02;
    COFD[4120] = -2.17926864E+01;
    COFD[4121] = 5.19232842E+00;
    COFD[4122] = -4.07643284E-01;
    COFD[4123] = 1.56246434E-02;
    COFD[4124] = -2.03971290E+01;
    COFD[4125] = 5.13279789E+00;
    COFD[4126] = -4.44474174E-01;
    COFD[4127] = 1.89937678E-02;
    COFD[4128] = -1.92334028E+01;
    COFD[4129] = 4.67033934E+00;
    COFD[4130] = -3.89971551E-01;
    COFD[4131] = 1.68513441E-02;
    COFD[4132] = -2.21713935E+01;
    COFD[4133] = 5.41196486E+00;
    COFD[4134] = -4.45632422E-01;
    COFD[4135] = 1.76474237E-02;
    COFD[4136] = -2.22600844E+01;
    COFD[4137] = 5.59632316E+00;
    COFD[4138] = -4.91568011E-01;
    COFD[4139] = 2.05156966E-02;
    COFD[4140] = -2.03526104E+01;
    COFD[4141] = 5.11453301E+00;
    COFD[4142] = -4.42447016E-01;
    COFD[4143] = 1.89196698E-02;
    COFD[4144] = -2.19617258E+01;
    COFD[4145] = 5.57026255E+00;
    COFD[4146] = -4.89178491E-01;
    COFD[4147] = 2.04505218E-02;
    COFD[4148] = -2.22701953E+01;
    COFD[4149] = 5.59632316E+00;
    COFD[4150] = -4.91568011E-01;
    COFD[4151] = 2.05156966E-02;
    COFD[4152] = -2.22052004E+01;
    COFD[4153] = 5.58604166E+00;
    COFD[4154] = -4.90602184E-01;
    COFD[4155] = 2.04880352E-02;
    COFD[4156] = -2.17825544E+01;
    COFD[4157] = 5.19232842E+00;
    COFD[4158] = -4.07643284E-01;
    COFD[4159] = 1.56246434E-02;
    COFD[4160] = -2.19919464E+01;
    COFD[4161] = 5.23595129E+00;
    COFD[4162] = -4.15079064E-01;
    COFD[4163] = 1.60168286E-02;
    COFD[4164] = -2.19969874E+01;
    COFD[4165] = 5.23595129E+00;
    COFD[4166] = -4.15079064E-01;
    COFD[4167] = 1.60168286E-02;
    COFD[4168] = -2.21795362E+01;
    COFD[4169] = 5.61233637E+00;
    COFD[4170] = -4.91419253E-01;
    COFD[4171] = 2.04216738E-02;
    COFD[4172] = -2.04750337E+01;
    COFD[4173] = 5.15745622E+00;
    COFD[4174] = -4.46648283E-01;
    COFD[4175] = 1.90458987E-02;
    COFD[4176] = -2.21912885E+01;
    COFD[4177] = 5.61233637E+00;
    COFD[4178] = -4.91419253E-01;
    COFD[4179] = 2.04216738E-02;
    COFD[4180] = -2.22604974E+01;
    COFD[4181] = 5.36643605E+00;
    COFD[4182] = -4.37440735E-01;
    COFD[4183] = 1.72016388E-02;
    COFD[4184] = -2.25161211E+01;
    COFD[4185] = 5.58521783E+00;
    COFD[4186] = -4.80947522E-01;
    COFD[4187] = 1.96897222E-02;
    COFD[4188] = -2.24797372E+01;
    COFD[4189] = 5.58492389E+00;
    COFD[4190] = -4.81921515E-01;
    COFD[4191] = 1.97721229E-02;
    COFD[4192] = -2.25933800E+01;
    COFD[4193] = 5.42562556E+00;
    COFD[4194] = -4.48132861E-01;
    COFD[4195] = 1.77847329E-02;
    COFD[4196] = -2.26149345E+01;
    COFD[4197] = 5.58414475E+00;
    COFD[4198] = -4.82375215E-01;
    COFD[4199] = 1.98138565E-02;
    COFD[4200] = -2.26149345E+01;
    COFD[4201] = 5.58414475E+00;
    COFD[4202] = -4.82375215E-01;
    COFD[4203] = 1.98138565E-02;
    COFD[4204] = -2.26206818E+01;
    COFD[4205] = 5.58414475E+00;
    COFD[4206] = -4.82375215E-01;
    COFD[4207] = 1.98138565E-02;
    COFD[4208] = -2.24927410E+01;
    COFD[4209] = 5.29529316E+00;
    COFD[4210] = -4.25345414E-01;
    COFD[4211] = 1.65635428E-02;
    COFD[4212] = -2.23925982E+01;
    COFD[4213] = 5.23666690E+00;
    COFD[4214] = -4.15204403E-01;
    COFD[4215] = 1.60235416E-02;
    COFD[4216] = -1.92108129E+01;
    COFD[4217] = 4.66318669E+00;
    COFD[4218] = -3.89108667E-01;
    COFD[4219] = 1.68165377E-02;
    COFD[4220] = -1.80755831E+01;
    COFD[4221] = 4.40247898E+00;
    COFD[4222] = -3.57238362E-01;
    COFD[4223] = 1.55145651E-02;
    COFD[4224] = -1.87484393E+01;
    COFD[4225] = 5.13858656E+00;
    COFD[4226] = -4.45075387E-01;
    COFD[4227] = 1.90137309E-02;
    COFD[4228] = -1.80978142E+01;
    COFD[4229] = 4.40247898E+00;
    COFD[4230] = -3.57238362E-01;
    COFD[4231] = 1.55145651E-02;
    COFD[4232] = -1.47143050E+01;
    COFD[4233] = 3.48023191E+00;
    COFD[4234] = -2.40800798E-01;
    COFD[4235] = 1.05947990E-02;
    COFD[4236] = -1.95250774E+01;
    COFD[4237] = 4.77151544E+00;
    COFD[4238] = -4.01882811E-01;
    COFD[4239] = 1.73184814E-02;
    COFD[4240] = -2.08912977E+01;
    COFD[4241] = 4.92602269E+00;
    COFD[4242] = -3.64572914E-01;
    COFD[4243] = 1.34203681E-02;
    COFD[4244] = -1.95425515E+01;
    COFD[4245] = 4.77151544E+00;
    COFD[4246] = -4.01882811E-01;
    COFD[4247] = 1.73184814E-02;
    COFD[4248] = -1.95340050E+01;
    COFD[4249] = 4.77151544E+00;
    COFD[4250] = -4.01882811E-01;
    COFD[4251] = 1.73184814E-02;
    COFD[4252] = -2.03739934E+01;
    COFD[4253] = 5.13279789E+00;
    COFD[4254] = -4.44474174E-01;
    COFD[4255] = 1.89937678E-02;
    COFD[4256] = -2.17974021E+01;
    COFD[4257] = 5.19232842E+00;
    COFD[4258] = -4.07643284E-01;
    COFD[4259] = 1.56246434E-02;
    COFD[4260] = -2.04000941E+01;
    COFD[4261] = 5.13279789E+00;
    COFD[4262] = -4.44474174E-01;
    COFD[4263] = 1.89937678E-02;
    COFD[4264] = -1.92379258E+01;
    COFD[4265] = 4.67033934E+00;
    COFD[4266] = -3.89971551E-01;
    COFD[4267] = 1.68513441E-02;
    COFD[4268] = -2.21762017E+01;
    COFD[4269] = 5.41196486E+00;
    COFD[4270] = -4.45632422E-01;
    COFD[4271] = 1.76474237E-02;
    COFD[4272] = -2.22647093E+01;
    COFD[4273] = 5.59632316E+00;
    COFD[4274] = -4.91568011E-01;
    COFD[4275] = 2.05156966E-02;
    COFD[4276] = -2.03557208E+01;
    COFD[4277] = 5.11453301E+00;
    COFD[4278] = -4.42447016E-01;
    COFD[4279] = 1.89196698E-02;
    COFD[4280] = -2.19662531E+01;
    COFD[4281] = 5.57026255E+00;
    COFD[4282] = -4.89178491E-01;
    COFD[4283] = 2.04505218E-02;
    COFD[4284] = -2.22749151E+01;
    COFD[4285] = 5.59632316E+00;
    COFD[4286] = -4.91568011E-01;
    COFD[4287] = 2.05156966E-02;
    COFD[4288] = -2.22110089E+01;
    COFD[4289] = 5.58604166E+00;
    COFD[4290] = -4.90602184E-01;
    COFD[4291] = 2.04880352E-02;
    COFD[4292] = -2.17871751E+01;
    COFD[4293] = 5.19232842E+00;
    COFD[4294] = -4.07643284E-01;
    COFD[4295] = 1.56246434E-02;
    COFD[4296] = -2.19979468E+01;
    COFD[4297] = 5.23595129E+00;
    COFD[4298] = -4.15079064E-01;
    COFD[4299] = 1.60168286E-02;
    COFD[4300] = -2.20030490E+01;
    COFD[4301] = 5.23595129E+00;
    COFD[4302] = -4.15079064E-01;
    COFD[4303] = 1.60168286E-02;
    COFD[4304] = -2.21838598E+01;
    COFD[4305] = 5.61233637E+00;
    COFD[4306] = -4.91419253E-01;
    COFD[4307] = 2.04216738E-02;
    COFD[4308] = -2.04806396E+01;
    COFD[4309] = 5.15745622E+00;
    COFD[4310] = -4.46648283E-01;
    COFD[4311] = 1.90458987E-02;
    COFD[4312] = -2.21957154E+01;
    COFD[4313] = 5.61233637E+00;
    COFD[4314] = -4.91419253E-01;
    COFD[4315] = 2.04216738E-02;
    COFD[4316] = -2.22662419E+01;
    COFD[4317] = 5.36643605E+00;
    COFD[4318] = -4.37440735E-01;
    COFD[4319] = 1.72016388E-02;
    COFD[4320] = -2.25217301E+01;
    COFD[4321] = 5.58521783E+00;
    COFD[4322] = -4.80947522E-01;
    COFD[4323] = 1.96897222E-02;
    COFD[4324] = -2.24854161E+01;
    COFD[4325] = 5.58492389E+00;
    COFD[4326] = -4.81921515E-01;
    COFD[4327] = 1.97721229E-02;
    COFD[4328] = -2.25999351E+01;
    COFD[4329] = 5.42562556E+00;
    COFD[4330] = -4.48132861E-01;
    COFD[4331] = 1.77847329E-02;
    COFD[4332] = -2.26206818E+01;
    COFD[4333] = 5.58414475E+00;
    COFD[4334] = -4.82375215E-01;
    COFD[4335] = 1.98138565E-02;
    COFD[4336] = -2.26206818E+01;
    COFD[4337] = 5.58414475E+00;
    COFD[4338] = -4.82375215E-01;
    COFD[4339] = 1.98138565E-02;
    COFD[4340] = -2.26264961E+01;
    COFD[4341] = 5.58414475E+00;
    COFD[4342] = -4.82375215E-01;
    COFD[4343] = 1.98138565E-02;
    COFD[4344] = -2.25000560E+01;
    COFD[4345] = 5.29529316E+00;
    COFD[4346] = -4.25345414E-01;
    COFD[4347] = 1.65635428E-02;
    COFD[4348] = -2.23999132E+01;
    COFD[4349] = 5.23666690E+00;
    COFD[4350] = -4.15204403E-01;
    COFD[4351] = 1.60235416E-02;
    COFD[4352] = -2.08314053E+01;
    COFD[4353] = 5.16623522E+00;
    COFD[4354] = -4.47455294E-01;
    COFD[4355] = 1.90672494E-02;
    COFD[4356] = -1.96652729E+01;
    COFD[4357] = 4.92149009E+00;
    COFD[4358] = -4.19702348E-01;
    COFD[4359] = 1.80250370E-02;
    COFD[4360] = -1.99919894E+01;
    COFD[4361] = 5.50111590E+00;
    COFD[4362] = -4.82434754E-01;
    COFD[4363] = 2.02461869E-02;
    COFD[4364] = -1.96903181E+01;
    COFD[4365] = 4.92149009E+00;
    COFD[4366] = -4.19702348E-01;
    COFD[4367] = 1.80250370E-02;
    COFD[4368] = -1.63371298E+01;
    COFD[4369] = 4.06047384E+00;
    COFD[4370] = -3.15033719E-01;
    COFD[4371] = 1.37721355E-02;
    COFD[4372] = -2.12147591E+01;
    COFD[4373] = 5.29477845E+00;
    COFD[4374] = -4.62516356E-01;
    COFD[4375] = 1.96564008E-02;
    COFD[4376] = -1.76259860E+01;
    COFD[4377] = 3.29707321E+00;
    COFD[4378] = -1.19066759E-01;
    COFD[4379] = 1.47805843E-03;
    COFD[4380] = -2.11103966E+01;
    COFD[4381] = 5.25151362E+00;
    COFD[4382] = -4.57337706E-01;
    COFD[4383] = 1.94489055E-02;
    COFD[4384] = -2.10999968E+01;
    COFD[4385] = 5.25151362E+00;
    COFD[4386] = -4.57337706E-01;
    COFD[4387] = 1.94489055E-02;
    COFD[4388] = -2.16464729E+01;
    COFD[4389] = 5.49325617E+00;
    COFD[4390] = -4.81586493E-01;
    COFD[4391] = 2.02162204E-02;
    COFD[4392] = -2.01633664E+01;
    COFD[4393] = 4.26535185E+00;
    COFD[4394] = -2.61127236E-01;
    COFD[4395] = 8.24369619E-03;
    COFD[4396] = -2.16755464E+01;
    COFD[4397] = 5.49325617E+00;
    COFD[4398] = -4.81586493E-01;
    COFD[4399] = 2.02162204E-02;
    COFD[4400] = -2.08660787E+01;
    COFD[4401] = 5.17548471E+00;
    COFD[4402] = -4.48496897E-01;
    COFD[4403] = 1.91060173E-02;
    COFD[4404] = -2.03466999E+01;
    COFD[4405] = 4.38969526E+00;
    COFD[4406] = -2.84615363E-01;
    COFD[4407] = 9.55641838E-03;
    COFD[4408] = -2.25524874E+01;
    COFD[4409] = 5.50138987E+00;
    COFD[4410] = -4.62503139E-01;
    COFD[4411] = 1.85883715E-02;
    COFD[4412] = -2.17407707E+01;
    COFD[4413] = 5.50870827E+00;
    COFD[4414] = -4.83264910E-01;
    COFD[4415] = 2.02760941E-02;
    COFD[4416] = -2.23572782E+01;
    COFD[4417] = 5.52017185E+00;
    COFD[4418] = -4.66659932E-01;
    COFD[4419] = 1.88373100E-02;
    COFD[4420] = -2.25647193E+01;
    COFD[4421] = 5.50138987E+00;
    COFD[4422] = -4.62503139E-01;
    COFD[4423] = 1.85883715E-02;
    COFD[4424] = -2.25695646E+01;
    COFD[4425] = 5.49903771E+00;
    COFD[4426] = -4.61784834E-01;
    COFD[4427] = 1.85410072E-02;
    COFD[4428] = -2.01511112E+01;
    COFD[4429] = 4.26535185E+00;
    COFD[4430] = -2.61127236E-01;
    COFD[4431] = 8.24369619E-03;
    COFD[4432] = -2.05128005E+01;
    COFD[4433] = 4.35981364E+00;
    COFD[4434] = -2.75585363E-01;
    COFD[4435] = 8.95565009E-03;
    COFD[4436] = -2.05192927E+01;
    COFD[4437] = 4.35981364E+00;
    COFD[4438] = -2.75585363E-01;
    COFD[4439] = 8.95565009E-03;
    COFD[4440] = -2.23166846E+01;
    COFD[4441] = 5.44889462E+00;
    COFD[4442] = -4.52325002E-01;
    COFD[4443] = 1.80132629E-02;
    COFD[4444] = -2.18644512E+01;
    COFD[4445] = 5.53438558E+00;
    COFD[4446] = -4.85873757E-01;
    COFD[4447] = 2.03606278E-02;
    COFD[4448] = -2.23307159E+01;
    COFD[4449] = 5.44889462E+00;
    COFD[4450] = -4.52325002E-01;
    COFD[4451] = 1.80132629E-02;
    COFD[4452] = -2.10617067E+01;
    COFD[4453] = 4.61245983E+00;
    COFD[4454] = -3.14873896E-01;
    COFD[4455] = 1.09118770E-02;
    COFD[4456] = -2.22719261E+01;
    COFD[4457] = 5.24354480E+00;
    COFD[4458] = -4.16412182E-01;
    COFD[4459] = 1.60883090E-02;
    COFD[4460] = -2.23086362E+01;
    COFD[4461] = 5.27409919E+00;
    COFD[4462] = -4.21767988E-01;
    COFD[4463] = 1.63754073E-02;
    COFD[4464] = -2.15763465E+01;
    COFD[4465] = 4.74313639E+00;
    COFD[4466] = -3.35584534E-01;
    COFD[4467] = 1.19568362E-02;
    COFD[4468] = -2.24927410E+01;
    COFD[4469] = 5.29529316E+00;
    COFD[4470] = -4.25345414E-01;
    COFD[4471] = 1.65635428E-02;
    COFD[4472] = -2.24927410E+01;
    COFD[4473] = 5.29529316E+00;
    COFD[4474] = -4.25345414E-01;
    COFD[4475] = 1.65635428E-02;
    COFD[4476] = -2.25000560E+01;
    COFD[4477] = 5.29529316E+00;
    COFD[4478] = -4.25345414E-01;
    COFD[4479] = 1.65635428E-02;
    COFD[4480] = -2.08265545E+01;
    COFD[4481] = 4.32912651E+00;
    COFD[4482] = -2.73390486E-01;
    COFD[4483] = 8.93695796E-03;
    COFD[4484] = -2.05649896E+01;
    COFD[4485] = 4.20297208E+00;
    COFD[4486] = -2.53989253E-01;
    COFD[4487] = 7.97846764E-03;
    COFD[4488] = -2.09943481E+01;
    COFD[4489] = 5.22468467E+00;
    COFD[4490] = -4.54220128E-01;
    COFD[4491] = 1.93281042E-02;
    COFD[4492] = -1.98296243E+01;
    COFD[4493] = 4.98207523E+00;
    COFD[4494] = -4.26877291E-01;
    COFD[4495] = 1.83086094E-02;
    COFD[4496] = -2.01262921E+01;
    COFD[4497] = 5.54581286E+00;
    COFD[4498] = -4.87014004E-01;
    COFD[4499] = 2.03965482E-02;
    COFD[4500] = -1.98546695E+01;
    COFD[4501] = 4.98207523E+00;
    COFD[4502] = -4.26877291E-01;
    COFD[4503] = 1.83086094E-02;
    COFD[4504] = -1.64819183E+01;
    COFD[4505] = 4.11726215E+00;
    COFD[4506] = -3.22193015E-01;
    COFD[4507] = 1.40747074E-02;
    COFD[4508] = -2.13698722E+01;
    COFD[4509] = 5.34971865E+00;
    COFD[4510] = -4.68771123E-01;
    COFD[4511] = 1.98933811E-02;
    COFD[4512] = -1.73636900E+01;
    COFD[4513] = 3.17377130E+00;
    COFD[4514] = -1.00394383E-01;
    COFD[4515] = 5.69083899E-04;
    COFD[4516] = -2.13011157E+01;
    COFD[4517] = 5.32167660E+00;
    COFD[4518] = -4.65740624E-01;
    COFD[4519] = 1.97861081E-02;
    COFD[4520] = -2.12907159E+01;
    COFD[4521] = 5.32167660E+00;
    COFD[4522] = -4.65740624E-01;
    COFD[4523] = 1.97861081E-02;
    COFD[4524] = -2.17867314E+01;
    COFD[4525] = 5.53950393E+00;
    COFD[4526] = -4.86376204E-01;
    COFD[4527] = 2.03760106E-02;
    COFD[4528] = -1.98359760E+01;
    COFD[4529] = 4.11158627E+00;
    COFD[4530] = -2.37831519E-01;
    COFD[4531] = 7.10363413E-03;
    COFD[4532] = -2.18158049E+01;
    COFD[4533] = 5.53950393E+00;
    COFD[4534] = -4.86376204E-01;
    COFD[4535] = 2.03760106E-02;
    COFD[4536] = -2.10310742E+01;
    COFD[4537] = 5.23485505E+00;
    COFD[4538] = -4.55400362E-01;
    COFD[4539] = 1.93737680E-02;
    COFD[4540] = -2.01613414E+01;
    COFD[4541] = 4.29679630E+00;
    COFD[4542] = -2.69916064E-01;
    COFD[4543] = 8.81737046E-03;
    COFD[4544] = -2.25180193E+01;
    COFD[4545] = 5.47136127E+00;
    COFD[4546] = -4.56417141E-01;
    COFD[4547] = 1.82376994E-02;
    COFD[4548] = -2.18731920E+01;
    COFD[4549] = 5.55171660E+00;
    COFD[4550] = -4.87609504E-01;
    COFD[4551] = 2.04156590E-02;
    COFD[4552] = -2.23434237E+01;
    COFD[4553] = 5.49927389E+00;
    COFD[4554] = -4.61845436E-01;
    COFD[4555] = 1.85448066E-02;
    COFD[4556] = -2.25302512E+01;
    COFD[4557] = 5.47136127E+00;
    COFD[4558] = -4.56417141E-01;
    COFD[4559] = 1.82376994E-02;
    COFD[4560] = -2.25168081E+01;
    COFD[4561] = 5.46125558E+00;
    COFD[4562] = -4.54580949E-01;
    COFD[4563] = 1.81370928E-02;
    COFD[4564] = -1.98237209E+01;
    COFD[4565] = 4.11158627E+00;
    COFD[4566] = -2.37831519E-01;
    COFD[4567] = 7.10363413E-03;
    COFD[4568] = -2.02246117E+01;
    COFD[4569] = 4.22278378E+00;
    COFD[4570] = -2.54653500E-01;
    COFD[4571] = 7.92616085E-03;
    COFD[4572] = -2.02311039E+01;
    COFD[4573] = 4.22278378E+00;
    COFD[4574] = -2.54653500E-01;
    COFD[4575] = 7.92616085E-03;
    COFD[4576] = -2.22462130E+01;
    COFD[4577] = 5.40356304E+00;
    COFD[4578] = -4.44060256E-01;
    COFD[4579] = 1.75601121E-02;
    COFD[4580] = -2.19764159E+01;
    COFD[4581] = 5.56943713E+00;
    COFD[4582] = -4.89114655E-01;
    COFD[4583] = 2.04494661E-02;
    COFD[4584] = -2.22602443E+01;
    COFD[4585] = 5.40356304E+00;
    COFD[4586] = -4.44060256E-01;
    COFD[4587] = 1.75601121E-02;
    COFD[4588] = -2.08429322E+01;
    COFD[4589] = 4.50409026E+00;
    COFD[4590] = -2.97868419E-01;
    COFD[4591] = 1.00604224E-02;
    COFD[4592] = -2.21603646E+01;
    COFD[4593] = 5.18050127E+00;
    COFD[4594] = -4.05688517E-01;
    COFD[4595] = 1.55231713E-02;
    COFD[4596] = -2.22169882E+01;
    COFD[4597] = 5.21950983E+00;
    COFD[4598] = -4.12223195E-01;
    COFD[4599] = 1.58645894E-02;
    COFD[4600] = -2.13682139E+01;
    COFD[4601] = 4.63852747E+00;
    COFD[4602] = -3.19031243E-01;
    COFD[4603] = 1.11224841E-02;
    COFD[4604] = -2.23925982E+01;
    COFD[4605] = 5.23666690E+00;
    COFD[4606] = -4.15204403E-01;
    COFD[4607] = 1.60235416E-02;
    COFD[4608] = -2.23925982E+01;
    COFD[4609] = 5.23666690E+00;
    COFD[4610] = -4.15204403E-01;
    COFD[4611] = 1.60235416E-02;
    COFD[4612] = -2.23999132E+01;
    COFD[4613] = 5.23666690E+00;
    COFD[4614] = -4.15204403E-01;
    COFD[4615] = 1.60235416E-02;
    COFD[4616] = -2.05649896E+01;
    COFD[4617] = 4.20297208E+00;
    COFD[4618] = -2.53989253E-01;
    COFD[4619] = 7.97846764E-03;
    COFD[4620] = -2.02828056E+01;
    COFD[4621] = 4.06866060E+00;
    COFD[4622] = -2.33527101E-01;
    COFD[4623] = 6.97454219E-03;
}


/*List of specs with small weight, dim NLITE */
void egtransetKTDIF(int* KTDIF) {
    KTDIF[0] = 2;
    KTDIF[1] = 4;
}


/*Poly fits for thermal diff ratios, dim NO*NLITE*KK */
void egtransetCOFTD(amrex::Real* COFTD) {
    COFTD[0] = 2.01521643E-01;
    COFTD[1] = 5.62744089E-04;
    COFTD[2] = -3.08519239E-07;
    COFTD[3] = 5.22805986E-11;
    COFTD[4] = 2.35283119E-01;
    COFTD[5] = 4.65670599E-04;
    COFTD[6] = -2.60939824E-07;
    COFTD[7] = 4.49271822E-11;
    COFTD[8] = 0.00000000E+00;
    COFTD[9] = 0.00000000E+00;
    COFTD[10] = 0.00000000E+00;
    COFTD[11] = 0.00000000E+00;
    COFTD[12] = 2.37053352E-01;
    COFTD[13] = 4.69174231E-04;
    COFTD[14] = -2.62903094E-07;
    COFTD[15] = 4.52652072E-11;
    COFTD[16] = 1.44152190E-01;
    COFTD[17] = 7.99993584E-05;
    COFTD[18] = -4.89707442E-08;
    COFTD[19] = 9.14277269E-12;
    COFTD[20] = 1.79840299E-01;
    COFTD[21] = 6.01722902E-04;
    COFTD[22] = -3.26433894E-07;
    COFTD[23] = 5.49112302E-11;
    COFTD[24] = -1.74352698E-01;
    COFTD[25] = 8.62246873E-04;
    COFTD[26] = -3.79545489E-07;
    COFTD[27] = 5.60262093E-11;
    COFTD[28] = 1.80513677E-01;
    COFTD[29] = 6.03975942E-04;
    COFTD[30] = -3.27656165E-07;
    COFTD[31] = 5.51168351E-11;
    COFTD[32] = 1.80186965E-01;
    COFTD[33] = 6.02882805E-04;
    COFTD[34] = -3.27063140E-07;
    COFTD[35] = 5.50170790E-11;
    COFTD[36] = 9.90752318E-02;
    COFTD[37] = 6.44201384E-04;
    COFTD[38] = -3.38485953E-07;
    COFTD[39] = 5.57356746E-11;
    COFTD[40] = -1.61357564E-01;
    COFTD[41] = 9.05920260E-04;
    COFTD[42] = -4.07879153E-07;
    COFTD[43] = 6.10626290E-11;
    COFTD[44] = 1.00039110E-01;
    COFTD[45] = 6.50468660E-04;
    COFTD[46] = -3.41778999E-07;
    COFTD[47] = 5.62779132E-11;
    COFTD[48] = 2.00119897E-01;
    COFTD[49] = 5.64793704E-04;
    COFTD[50] = -3.09445484E-07;
    COFTD[51] = 5.24139335E-11;
    COFTD[52] = -1.31244519E-01;
    COFTD[53] = 9.03901384E-04;
    COFTD[54] = -4.17831507E-07;
    COFTD[55] = 6.35725667E-11;
    COFTD[56] = -2.28105944E-02;
    COFTD[57] = 8.33470403E-04;
    COFTD[58] = -4.11969112E-07;
    COFTD[59] = 6.52859371E-11;
    COFTD[60] = 1.05124122E-01;
    COFTD[61] = 6.50665957E-04;
    COFTD[62] = -3.42564538E-07;
    COFTD[63] = 5.64804120E-11;
    COFTD[64] = -1.42100396E-02;
    COFTD[65] = 8.23812102E-04;
    COFTD[66] = -4.08995515E-07;
    COFTD[67] = 6.49899310E-11;
    COFTD[68] = -2.28637575E-02;
    COFTD[69] = 8.35412914E-04;
    COFTD[70] = -4.12929260E-07;
    COFTD[71] = 6.54380945E-11;
    COFTD[72] = -2.00309448E-02;
    COFTD[73] = 8.50440115E-04;
    COFTD[74] = -4.21064468E-07;
    COFTD[75] = 6.67959710E-11;
    COFTD[76] = -1.60981264E-01;
    COFTD[77] = 9.03807572E-04;
    COFTD[78] = -4.06927941E-07;
    COFTD[79] = 6.09202254E-11;
    COFTD[80] = -1.59826932E-01;
    COFTD[81] = 9.28231324E-04;
    COFTD[82] = -4.20059750E-07;
    COFTD[83] = 6.30844146E-11;
    COFTD[84] = -1.59970790E-01;
    COFTD[85] = 9.29066816E-04;
    COFTD[86] = -4.20437842E-07;
    COFTD[87] = 6.31411962E-11;
    COFTD[88] = -3.81470765E-02;
    COFTD[89] = 8.39833490E-04;
    COFTD[90] = -4.11688915E-07;
    COFTD[91] = 6.49124952E-11;
    COFTD[92] = 9.86934401E-02;
    COFTD[93] = 7.20974863E-04;
    COFTD[94] = -3.77135221E-07;
    COFTD[95] = 6.19202579E-11;
    COFTD[96] = -3.82574649E-02;
    COFTD[97] = 8.42263764E-04;
    COFTD[98] = -4.12880242E-07;
    COFTD[99] = 6.51003362E-11;
    COFTD[100] = -1.41799739E-01;
    COFTD[101] = 9.22440172E-04;
    COFTD[102] = -4.23685885E-07;
    COFTD[103] = 6.42121388E-11;
    COFTD[104] = -7.78657454E-02;
    COFTD[105] = 8.92101015E-04;
    COFTD[106] = -4.27969255E-07;
    COFTD[107] = 6.66000503E-11;
    COFTD[108] = -7.23038994E-02;
    COFTD[109] = 8.89466098E-04;
    COFTD[110] = -4.28124818E-07;
    COFTD[111] = 6.67586244E-11;
    COFTD[112] = -1.32461416E-01;
    COFTD[113] = 9.30313107E-04;
    COFTD[114] = -4.30922433E-07;
    COFTD[115] = 6.56463431E-11;
    COFTD[116] = -6.92602151E-02;
    COFTD[117] = 8.88360172E-04;
    COFTD[118] = -4.28365765E-07;
    COFTD[119] = 6.68694606E-11;
    COFTD[120] = -6.92602151E-02;
    COFTD[121] = 8.88360172E-04;
    COFTD[122] = -4.28365765E-07;
    COFTD[123] = 6.68694606E-11;
    COFTD[124] = -6.93343621E-02;
    COFTD[125] = 8.89311212E-04;
    COFTD[126] = -4.28824355E-07;
    COFTD[127] = 6.69410482E-11;
    COFTD[128] = -1.54148352E-01;
    COFTD[129] = 9.42617019E-04;
    COFTD[130] = -4.29628100E-07;
    COFTD[131] = 6.48056821E-11;
    COFTD[132] = -1.62301128E-01;
    COFTD[133] = 9.43217155E-04;
    COFTD[134] = -4.26881994E-07;
    COFTD[135] = 6.41127358E-11;
    COFTD[136] = 4.31331269E-01;
    COFTD[137] = 9.20536800E-05;
    COFTD[138] = -5.94509616E-08;
    COFTD[139] = 1.21437993E-11;
    COFTD[140] = 4.06682492E-01;
    COFTD[141] = 3.84705248E-05;
    COFTD[142] = -2.54846868E-08;
    COFTD[143] = 5.86302354E-12;
    COFTD[144] = -1.44152190E-01;
    COFTD[145] = -7.99993584E-05;
    COFTD[146] = 4.89707442E-08;
    COFTD[147] = -9.14277269E-12;
    COFTD[148] = 4.12895615E-01;
    COFTD[149] = 3.90582612E-05;
    COFTD[150] = -2.58740310E-08;
    COFTD[151] = 5.95259633E-12;
    COFTD[152] = 0.00000000E+00;
    COFTD[153] = 0.00000000E+00;
    COFTD[154] = 0.00000000E+00;
    COFTD[155] = 0.00000000E+00;
    COFTD[156] = 4.26579943E-01;
    COFTD[157] = 1.20407274E-04;
    COFTD[158] = -7.67298757E-08;
    COFTD[159] = 1.52090336E-11;
    COFTD[160] = 2.27469146E-02;
    COFTD[161] = 6.73078907E-04;
    COFTD[162] = -3.40935843E-07;
    COFTD[163] = 5.48499211E-11;
    COFTD[164] = 4.29789463E-01;
    COFTD[165] = 1.21313199E-04;
    COFTD[166] = -7.73071792E-08;
    COFTD[167] = 1.53234639E-11;
    COFTD[168] = 4.28230888E-01;
    COFTD[169] = 1.20873273E-04;
    COFTD[170] = -7.70268349E-08;
    COFTD[171] = 1.52678954E-11;
    COFTD[172] = 3.24747031E-01;
    COFTD[173] = 1.77798548E-04;
    COFTD[174] = -1.08934732E-07;
    COFTD[175] = 2.03595881E-11;
    COFTD[176] = 1.22693382E-01;
    COFTD[177] = 6.21278143E-04;
    COFTD[178] = -3.29965208E-07;
    COFTD[179] = 5.47161548E-11;
    COFTD[180] = 3.31191185E-01;
    COFTD[181] = 1.81326714E-04;
    COFTD[182] = -1.11096391E-07;
    COFTD[183] = 2.07635959E-11;
    COFTD[184] = 4.30605547E-01;
    COFTD[185] = 9.35961902E-05;
    COFTD[186] = -6.03983623E-08;
    COFTD[187] = 1.23115170E-11;
    COFTD[188] = 1.40314191E-01;
    COFTD[189] = 6.01266129E-04;
    COFTD[190] = -3.21915137E-07;
    COFTD[191] = 5.36679068E-11;
    COFTD[192] = 2.76725963E-01;
    COFTD[193] = 3.87792818E-04;
    COFTD[194] = -2.22504581E-07;
    COFTD[195] = 3.90251143E-11;
    COFTD[196] = 3.39557243E-01;
    COFTD[197] = 1.79335036E-04;
    COFTD[198] = -1.10135705E-07;
    COFTD[199] = 2.06427239E-11;
    COFTD[200] = 2.82974392E-01;
    COFTD[201] = 3.73032949E-04;
    COFTD[202] = -2.14959161E-07;
    COFTD[203] = 3.78355155E-11;
    COFTD[204] = 2.78021896E-01;
    COFTD[205] = 3.89608886E-04;
    COFTD[206] = -2.23546590E-07;
    COFTD[207] = 3.92078724E-11;
    COFTD[208] = 2.93191523E-01;
    COFTD[209] = 4.01430006E-04;
    COFTD[210] = -2.30705763E-07;
    COFTD[211] = 4.05176586E-11;
    COFTD[212] = 1.22119780E-01;
    COFTD[213] = 6.18373616E-04;
    COFTD[214] = -3.28422593E-07;
    COFTD[215] = 5.44603522E-11;
    COFTD[216] = 1.36817715E-01;
    COFTD[217] = 6.41727473E-04;
    COFTD[218] = -3.42055963E-07;
    COFTD[219] = 5.68567648E-11;
    COFTD[220] = 1.37064455E-01;
    COFTD[221] = 6.42884781E-04;
    COFTD[222] = -3.42672835E-07;
    COFTD[223] = 5.69593018E-11;
    COFTD[224] = 2.58066832E-01;
    COFTD[225] = 4.05072593E-04;
    COFTD[226] = -2.30587443E-07;
    COFTD[227] = 4.01863841E-11;
    COFTD[228] = 3.86107464E-01;
    COFTD[229] = 2.28760446E-04;
    COFTD[230] = -1.39425040E-07;
    COFTD[231] = 2.58989754E-11;
    COFTD[232] = 2.59569092E-01;
    COFTD[233] = 4.07430603E-04;
    COFTD[234] = -2.31929740E-07;
    COFTD[235] = 4.04203173E-11;
    COFTD[236] = 1.59647939E-01;
    COFTD[237] = 6.04192274E-04;
    COFTD[238] = -3.25569591E-07;
    COFTD[239] = 5.45134698E-11;
    COFTD[240] = 2.34098762E-01;
    COFTD[241] = 4.91099181E-04;
    COFTD[242] = -2.74133967E-07;
    COFTD[243] = 4.70636702E-11;
    COFTD[244] = 2.40639006E-01;
    COFTD[245] = 4.82930111E-04;
    COFTD[246] = -2.70362190E-07;
    COFTD[247] = 4.65173265E-11;
    COFTD[248] = 1.77767674E-01;
    COFTD[249] = 5.98129745E-04;
    COFTD[250] = -3.24382563E-07;
    COFTD[251] = 5.45543332E-11;
    COFTD[252] = 2.44452926E-01;
    COFTD[253] = 4.78884724E-04;
    COFTD[254] = -2.68527379E-07;
    COFTD[255] = 4.62572763E-11;
    COFTD[256] = 2.44452926E-01;
    COFTD[257] = 4.78884724E-04;
    COFTD[258] = -2.68527379E-07;
    COFTD[259] = 4.62572763E-11;
    COFTD[260] = 2.44977451E-01;
    COFTD[261] = 4.79912272E-04;
    COFTD[262] = -2.69103560E-07;
    COFTD[263] = 4.63565309E-11;
    COFTD[264] = 1.43330409E-01;
    COFTD[265] = 6.59904424E-04;
    COFTD[266] = -3.52059953E-07;
    COFTD[267] = 5.85545554E-11;
    COFTD[268] = 1.31648645E-01;
    COFTD[269] = 6.75329826E-04;
    COFTD[270] = -3.58458833E-07;
    COFTD[271] = 5.94176903E-11;
}

#ifndef AMREX_USE_GPU
/* Initializes parameter database */
void CKINIT()
{

    // (0):  O + H + M => OH + M
    kiv[9] = {1,2,3};
    nuv[9] = {-1,-1,1};
    // (0):  O + H + M => OH + M
    fwd_A[9]     = 4.72e+18;
    fwd_beta[9]  = -1;
    fwd_Ea[9]    = 0;
    prefactor_units[9]  = 1.0000000000000002e-12;
    activation_units[9] = 0.50321666580471969;
    phase_units[9]      = pow(10,-12.000000);
    is_PD[9] = 0;
    nTB[9] = 4;
    TB[9] = (amrex::Real *) malloc(4 * sizeof(amrex::Real));
    TBid[9] = (int *) malloc(4 * sizeof(int));
    TBid[9][0] = 4; TB[9][0] = 2.5; // H2
    TBid[9][1] = 6; TB[9][1] = 12; // H2O
    TBid[9][2] = 18; TB[9][2] = 3.7999999999999998; // CO2
    TBid[9][3] = 12; TB[9][3] = 1.8999999999999999; // CO

    // (1):  O + H2 => H + OH
    kiv[16] = {1,4,2,3};
    nuv[16] = {-1,-1,1,1};
    // (1):  O + H2 => H + OH
    fwd_A[16]     = 50800;
    fwd_beta[16]  = 2.6699999999999999;
    fwd_Ea[16]    = 6292.0699999999997;
    prefactor_units[16]  = 1.0000000000000002e-06;
    activation_units[16] = 0.50321666580471969;
    phase_units[16]      = pow(10,-12.000000);
    is_PD[16] = 0;
    nTB[16] = 0;

    // (2):  H + OH => O + H2
    kiv[17] = {2,3,1,4};
    nuv[17] = {-1,-1,1,1};
    // (2):  H + OH => O + H2
    fwd_A[17]     = 22310;
    fwd_beta[17]  = 2.6699999999999999;
    fwd_Ea[17]    = 4196.9399999999996;
    prefactor_units[17]  = 1.0000000000000002e-06;
    activation_units[17] = 0.50321666580471969;
    phase_units[17]      = pow(10,-12.000000);
    is_PD[17] = 0;
    nTB[17] = 0;

    // (3):  2.000000 O + M => O2 + M
    kiv[10] = {1,5};
    nuv[10] = {-2.0,1};
    // (3):  2.000000 O + M => O2 + M
    fwd_A[10]     = 6170000000000000;
    fwd_beta[10]  = -0.5;
    fwd_Ea[10]    = 0;
    prefactor_units[10]  = 1.0000000000000002e-12;
    activation_units[10] = 0.50321666580471969;
    phase_units[10]      = pow(10,-12.000000);
    is_PD[10] = 0;
    nTB[10] = 4;
    TB[10] = (amrex::Real *) malloc(4 * sizeof(amrex::Real));
    TBid[10] = (int *) malloc(4 * sizeof(int));
    TBid[10][0] = 4; TB[10][0] = 2.5; // H2
    TBid[10][1] = 6; TB[10][1] = 12; // H2O
    TBid[10][2] = 18; TB[10][2] = 3.7999999999999998; // CO2
    TBid[10][3] = 12; TB[10][3] = 1.8999999999999999; // CO

    // (4):  OH + H2 => H + H2O
    kiv[18] = {3,4,2,6};
    nuv[18] = {-1,-1,1,1};
    // (4):  OH + H2 => H + H2O
    fwd_A[18]     = 216000000;
    fwd_beta[18]  = 1.51;
    fwd_Ea[18]    = 3429.9699999999998;
    prefactor_units[18]  = 1.0000000000000002e-06;
    activation_units[18] = 0.50321666580471969;
    phase_units[18]      = pow(10,-12.000000);
    is_PD[18] = 0;
    nTB[18] = 0;

    // (5):  H + H2O => OH + H2
    kiv[19] = {2,6,3,4};
    nuv[19] = {-1,-1,1,1};
    // (5):  H + H2O => OH + H2
    fwd_A[19]     = 935200000;
    fwd_beta[19]  = 1.51;
    fwd_Ea[19]    = 18580.07;
    prefactor_units[19]  = 1.0000000000000002e-06;
    activation_units[19] = 0.50321666580471969;
    phase_units[19]      = pow(10,-12.000000);
    is_PD[19] = 0;
    nTB[19] = 0;

    // (6):  2.000000 OH (+M) <=> H2O2 (+M)
    kiv[0] = {3,7};
    nuv[0] = {-2.0,1};
    // (6):  2.000000 OH (+M) <=> H2O2 (+M)
    fwd_A[0]     = 123600000000000;
    fwd_beta[0]  = -0.37;
    fwd_Ea[0]    = 0;
    low_A[0]     = 3.0410000000000002e+30;
    low_beta[0]  = -4.6299999999999999;
    low_Ea[0]    = 2049;
    troe_a[0]    = 0.46999999999999997;
    troe_Tsss[0] = 100;
    troe_Ts[0]   = 2000;
    troe_Tss[0]  = 1000000000000000;
    troe_len[0]  = 4;
    prefactor_units[0]  = 1.0000000000000002e-06;
    activation_units[0] = 0.50321666580471969;
    phase_units[0]      = pow(10,-12.000000);
    is_PD[0] = 1;
    nTB[0] = 4;
    TB[0] = (amrex::Real *) malloc(4 * sizeof(amrex::Real));
    TBid[0] = (int *) malloc(4 * sizeof(int));
    TBid[0][0] = 4; TB[0][0] = 2.5; // H2
    TBid[0][1] = 6; TB[0][1] = 12; // H2O
    TBid[0][2] = 18; TB[0][2] = 3.7999999999999998; // CO2
    TBid[0][3] = 12; TB[0][3] = 1.8999999999999999; // CO

    // (7):  H + OH + M => H2O + M
    kiv[11] = {2,3,6};
    nuv[11] = {-1,-1,1};
    // (7):  H + OH + M => H2O + M
    fwd_A[11]     = 2.2499999999999999e+22;
    fwd_beta[11]  = -2;
    fwd_Ea[11]    = 0;
    prefactor_units[11]  = 1.0000000000000002e-12;
    activation_units[11] = 0.50321666580471969;
    phase_units[11]      = pow(10,-12.000000);
    is_PD[11] = 0;
    nTB[11] = 4;
    TB[11] = (amrex::Real *) malloc(4 * sizeof(amrex::Real));
    TBid[11] = (int *) malloc(4 * sizeof(int));
    TBid[11][0] = 4; TB[11][0] = 2.5; // H2
    TBid[11][1] = 6; TB[11][1] = 12; // H2O
    TBid[11][2] = 18; TB[11][2] = 3.7999999999999998; // CO2
    TBid[11][3] = 12; TB[11][3] = 1.8999999999999999; // CO

    // (8):  O + H2O => 2.000000 OH
    kiv[20] = {1,6,3};
    nuv[20] = {-1,-1,2.0};
    // (8):  O + H2O => 2.000000 OH
    fwd_A[20]     = 2970000;
    fwd_beta[20]  = 2.02;
    fwd_Ea[20]    = 13400.1;
    prefactor_units[20]  = 1.0000000000000002e-06;
    activation_units[20] = 0.50321666580471969;
    phase_units[20]      = pow(10,-12.000000);
    is_PD[20] = 0;
    nTB[20] = 0;

    // (9):  2.000000 OH => O + H2O
    kiv[21] = {3,1,6};
    nuv[21] = {-2.0,1,1};
    // (9):  2.000000 OH => O + H2O
    fwd_A[21]     = 301300;
    fwd_beta[21]  = 2.02;
    fwd_Ea[21]    = -3849.9000000000001;
    prefactor_units[21]  = 1.0000000000000002e-06;
    activation_units[21] = 0.50321666580471969;
    phase_units[21]      = pow(10,-12.000000);
    is_PD[21] = 0;
    nTB[21] = 0;

    // (10):  H + O2 => O + OH
    kiv[22] = {2,5,1,3};
    nuv[22] = {-1,-1,1,1};
    // (10):  H + O2 => O + OH
    fwd_A[22]     = 197000000000000;
    fwd_beta[22]  = 0;
    fwd_Ea[22]    = 16539.91;
    prefactor_units[22]  = 1.0000000000000002e-06;
    activation_units[22] = 0.50321666580471969;
    phase_units[22]      = pow(10,-12.000000);
    is_PD[22] = 0;
    nTB[22] = 0;

    // (11):  O + OH => H + O2
    kiv[23] = {1,3,2,5};
    nuv[23] = {-1,-1,1,1};
    // (11):  O + OH => H + O2
    fwd_A[23]     = 15550000000000;
    fwd_beta[23]  = 0;
    fwd_Ea[23]    = 424.94999999999999;
    prefactor_units[23]  = 1.0000000000000002e-06;
    activation_units[23] = 0.50321666580471969;
    phase_units[23]      = pow(10,-12.000000);
    is_PD[23] = 0;
    nTB[23] = 0;

    // (12):  H + O2 (+M) <=> HO2 (+M)
    kiv[1] = {2,5,8};
    nuv[1] = {-1,-1,1};
    // (12):  H + O2 (+M) <=> HO2 (+M)
    fwd_A[1]     = 1475000000000;
    fwd_beta[1]  = 0.59999999999999998;
    fwd_Ea[1]    = 0;
    low_A[1]     = 35000000000000000;
    low_beta[1]  = -0.40999999999999998;
    low_Ea[1]    = -1115.9200000000001;
    troe_a[1]    = 0.5;
    troe_Tsss[1] = 0;
    troe_Ts[1]   = 1e+30;
    troe_Tss[1]  = 1e+100;
    troe_len[1]  = 4;
    prefactor_units[1]  = 1.0000000000000002e-06;
    activation_units[1] = 0.50321666580471969;
    phase_units[1]      = pow(10,-12.000000);
    is_PD[1] = 1;
    nTB[1] = 4;
    TB[1] = (amrex::Real *) malloc(4 * sizeof(amrex::Real));
    TBid[1] = (int *) malloc(4 * sizeof(int));
    TBid[1][0] = 4; TB[1][0] = 2.5; // H2
    TBid[1][1] = 6; TB[1][1] = 12; // H2O
    TBid[1][2] = 18; TB[1][2] = 3.7999999999999998; // CO2
    TBid[1][3] = 12; TB[1][3] = 1.8999999999999999; // CO

    // (13):  HO2 + H => 2.000000 OH
    kiv[24] = {8,2,3};
    nuv[24] = {-1,-1,2.0};
    // (13):  HO2 + H => 2.000000 OH
    fwd_A[24]     = 70800000000000;
    fwd_beta[24]  = 0;
    fwd_Ea[24]    = 299.94999999999999;
    prefactor_units[24]  = 1.0000000000000002e-06;
    activation_units[24] = 0.50321666580471969;
    phase_units[24]      = pow(10,-12.000000);
    is_PD[24] = 0;
    nTB[24] = 0;

    // (14):  2.000000 HO2 => H2O2 + O2
    kiv[25] = {8,7,5};
    nuv[25] = {-2.0,1,1};
    // (14):  2.000000 HO2 => H2O2 + O2
    fwd_A[25]     = 420000000000000;
    fwd_beta[25]  = 0;
    fwd_Ea[25]    = 11979.92;
    prefactor_units[25]  = 1.0000000000000002e-06;
    activation_units[25] = 0.50321666580471969;
    phase_units[25]      = pow(10,-12.000000);
    is_PD[25] = 0;
    nTB[25] = 0;

    // (15):  HO2 + H => H2 + O2
    kiv[26] = {8,2,4,5};
    nuv[26] = {-1,-1,1,1};
    // (15):  HO2 + H => H2 + O2
    fwd_A[26]     = 16600000000000;
    fwd_beta[26]  = 0;
    fwd_Ea[26]    = 820.02999999999997;
    prefactor_units[26]  = 1.0000000000000002e-06;
    activation_units[26] = 0.50321666580471969;
    phase_units[26]      = pow(10,-12.000000);
    is_PD[26] = 0;
    nTB[26] = 0;

    // (16):  HO2 + OH => H2O + O2
    kiv[27] = {8,3,6,5};
    nuv[27] = {-1,-1,1,1};
    // (16):  HO2 + OH => H2O + O2
    fwd_A[27]     = 28900000000000;
    fwd_beta[27]  = 0;
    fwd_Ea[27]    = -500;
    prefactor_units[27]  = 1.0000000000000002e-06;
    activation_units[27] = 0.50321666580471969;
    phase_units[27]      = pow(10,-12.000000);
    is_PD[27] = 0;
    nTB[27] = 0;

    // (17):  H2O + O2 => HO2 + OH
    kiv[28] = {6,5,8,3};
    nuv[28] = {-1,-1,1,1};
    // (17):  H2O + O2 => HO2 + OH
    fwd_A[28]     = 6888000000000000;
    fwd_beta[28]  = -0.33000000000000002;
    fwd_Ea[28]    = 72140.059999999998;
    prefactor_units[28]  = 1.0000000000000002e-06;
    activation_units[28] = 0.50321666580471969;
    phase_units[28]      = pow(10,-12.000000);
    is_PD[28] = 0;
    nTB[28] = 0;

    // (18):  2.000000 HO2 => H2O2 + O2
    kiv[29] = {8,7,5};
    nuv[29] = {-2.0,1,1};
    // (18):  2.000000 HO2 => H2O2 + O2
    fwd_A[29]     = 130000000000;
    fwd_beta[29]  = 0;
    fwd_Ea[29]    = -1629.0599999999999;
    prefactor_units[29]  = 1.0000000000000002e-06;
    activation_units[29] = 0.50321666580471969;
    phase_units[29]      = pow(10,-12.000000);
    is_PD[29] = 0;
    nTB[29] = 0;

    // (19):  HO2 + O => OH + O2
    kiv[30] = {8,1,3,5};
    nuv[30] = {-1,-1,1,1};
    // (19):  HO2 + O => OH + O2
    fwd_A[30]     = 32500000000000;
    fwd_beta[30]  = 0;
    fwd_Ea[30]    = 0;
    prefactor_units[30]  = 1.0000000000000002e-06;
    activation_units[30] = 0.50321666580471969;
    phase_units[30]      = pow(10,-12.000000);
    is_PD[30] = 0;
    nTB[30] = 0;

    // (20):  OH + O2 => HO2 + O
    kiv[31] = {3,5,8,1};
    nuv[31] = {-1,-1,1,1};
    // (20):  OH + O2 => HO2 + O
    fwd_A[31]     = 785700000000000;
    fwd_beta[31]  = -0.33000000000000002;
    fwd_Ea[31]    = 55390.059999999998;
    prefactor_units[31]  = 1.0000000000000002e-06;
    activation_units[31] = 0.50321666580471969;
    phase_units[31]      = pow(10,-12.000000);
    is_PD[31] = 0;
    nTB[31] = 0;

    // (21):  H2O2 + OH => H2O + HO2
    kiv[32] = {7,3,6,8};
    nuv[32] = {-1,-1,1,1};
    // (21):  H2O2 + OH => H2O + HO2
    fwd_A[32]     = 1000000000000;
    fwd_beta[32]  = 0;
    fwd_Ea[32]    = 0;
    prefactor_units[32]  = 1.0000000000000002e-06;
    activation_units[32] = 0.50321666580471969;
    phase_units[32]      = pow(10,-12.000000);
    is_PD[32] = 0;
    nTB[32] = 0;

    // (22):  H2O2 + H => H2O + OH
    kiv[33] = {7,2,6,3};
    nuv[33] = {-1,-1,1,1};
    // (22):  H2O2 + H => H2O + OH
    fwd_A[33]     = 24100000000000;
    fwd_beta[33]  = 0;
    fwd_Ea[33]    = 3969.8899999999999;
    prefactor_units[33]  = 1.0000000000000002e-06;
    activation_units[33] = 0.50321666580471969;
    phase_units[33]      = pow(10,-12.000000);
    is_PD[33] = 0;
    nTB[33] = 0;

    // (23):  H2O2 + OH => H2O + HO2
    kiv[34] = {7,3,6,8};
    nuv[34] = {-1,-1,1,1};
    // (23):  H2O2 + OH => H2O + HO2
    fwd_A[34]     = 580000000000000;
    fwd_beta[34]  = 0;
    fwd_Ea[34]    = 9559.9899999999998;
    prefactor_units[34]  = 1.0000000000000002e-06;
    activation_units[34] = 0.50321666580471969;
    phase_units[34]      = pow(10,-12.000000);
    is_PD[34] = 0;
    nTB[34] = 0;

    // (24):  H2O + HO2 => H2O2 + OH
    kiv[35] = {6,8,7,3};
    nuv[35] = {-1,-1,1,1};
    // (24):  H2O + HO2 => H2O2 + OH
    fwd_A[35]     = 97710000000000;
    fwd_beta[35]  = 0.33000000000000002;
    fwd_Ea[35]    = 41020.080000000002;
    prefactor_units[35]  = 1.0000000000000002e-06;
    activation_units[35] = 0.50321666580471969;
    phase_units[35]      = pow(10,-12.000000);
    is_PD[35] = 0;
    nTB[35] = 0;

    // (25):  H2O2 + O => OH + HO2
    kiv[36] = {7,1,3,8};
    nuv[36] = {-1,-1,1,1};
    // (25):  H2O2 + O => OH + HO2
    fwd_A[36]     = 9550000;
    fwd_beta[36]  = 2;
    fwd_Ea[36]    = 3969.8899999999999;
    prefactor_units[36]  = 1.0000000000000002e-06;
    activation_units[36] = 0.50321666580471969;
    phase_units[36]      = pow(10,-12.000000);
    is_PD[36] = 0;
    nTB[36] = 0;

    // (26):  H2O2 + H => H2 + HO2
    kiv[37] = {7,2,4,8};
    nuv[37] = {-1,-1,1,1};
    // (26):  H2O2 + H => H2 + HO2
    fwd_A[37]     = 48200000000000;
    fwd_beta[37]  = 0;
    fwd_Ea[37]    = 7950.0500000000002;
    prefactor_units[37]  = 1.0000000000000002e-06;
    activation_units[37] = 0.50321666580471969;
    phase_units[37]      = pow(10,-12.000000);
    is_PD[37] = 0;
    nTB[37] = 0;

    // (27):  H2 + HO2 => H2O2 + H
    kiv[38] = {4,8,7,2};
    nuv[38] = {-1,-1,1,1};
    // (27):  H2 + HO2 => H2O2 + H
    fwd_A[38]     = 1875000000000;
    fwd_beta[38]  = 0.33000000000000002;
    fwd_Ea[38]    = 24260.040000000001;
    prefactor_units[38]  = 1.0000000000000002e-06;
    activation_units[38] = 0.50321666580471969;
    phase_units[38]      = pow(10,-12.000000);
    is_PD[38] = 0;
    nTB[38] = 0;

    // (28):  CH2GSG + OH => CH2O + H
    kiv[39] = {9,3,10,2};
    nuv[39] = {-1,-1,1,1};
    // (28):  CH2GSG + OH => CH2O + H
    fwd_A[39]     = 30000000000000;
    fwd_beta[39]  = 0;
    fwd_Ea[39]    = 0;
    prefactor_units[39]  = 1.0000000000000002e-06;
    activation_units[39] = 0.50321666580471969;
    phase_units[39]      = pow(10,-12.000000);
    is_PD[39] = 0;
    nTB[39] = 0;

    // (29):  CH2GSG + H2 => CH3 + H
    kiv[40] = {9,4,11,2};
    nuv[40] = {-1,-1,1,1};
    // (29):  CH2GSG + H2 => CH3 + H
    fwd_A[40]     = 70000000000000;
    fwd_beta[40]  = 0;
    fwd_Ea[40]    = 0;
    prefactor_units[40]  = 1.0000000000000002e-06;
    activation_units[40] = 0.50321666580471969;
    phase_units[40]      = pow(10,-12.000000);
    is_PD[40] = 0;
    nTB[40] = 0;

    // (30):  CH3 + H => CH2GSG + H2
    kiv[41] = {11,2,9,4};
    nuv[41] = {-1,-1,1,1};
    // (30):  CH3 + H => CH2GSG + H2
    fwd_A[41]     = 2.482e+17;
    fwd_beta[41]  = -0.89000000000000001;
    fwd_Ea[41]    = 16130.02;
    prefactor_units[41]  = 1.0000000000000002e-06;
    activation_units[41] = 0.50321666580471969;
    phase_units[41]      = pow(10,-12.000000);
    is_PD[41] = 0;
    nTB[41] = 0;

    // (31):  CH2GSG + O2 => CO + OH + H
    kiv[42] = {9,5,12,3,2};
    nuv[42] = {-1,-1,1,1,1};
    // (31):  CH2GSG + O2 => CO + OH + H
    fwd_A[42]     = 70000000000000;
    fwd_beta[42]  = 0;
    fwd_Ea[42]    = 0;
    prefactor_units[42]  = 1.0000000000000002e-06;
    activation_units[42] = 0.50321666580471969;
    phase_units[42]      = pow(10,-12.000000);
    is_PD[42] = 0;
    nTB[42] = 0;

    // (32):  CH2GSG + O => CO + 2.000000 H
    kiv[43] = {9,1,12,2};
    nuv[43] = {-1,-1,1,2.0};
    // (32):  CH2GSG + O => CO + 2.000000 H
    fwd_A[43]     = 30000000000000;
    fwd_beta[43]  = 0;
    fwd_Ea[43]    = 0;
    prefactor_units[43]  = 1.0000000000000002e-06;
    activation_units[43] = 0.50321666580471969;
    phase_units[43]      = pow(10,-12.000000);
    is_PD[43] = 0;
    nTB[43] = 0;

    // (33):  CH3 + HO2 => CH3O + OH
    kiv[44] = {11,8,13,3};
    nuv[44] = {-1,-1,1,1};
    // (33):  CH3 + HO2 => CH3O + OH
    fwd_A[44]     = 11000000000000;
    fwd_beta[44]  = 0;
    fwd_Ea[44]    = 0;
    prefactor_units[44]  = 1.0000000000000002e-06;
    activation_units[44] = 0.50321666580471969;
    phase_units[44]      = pow(10,-12.000000);
    is_PD[44] = 0;
    nTB[44] = 0;

    // (34):  CH3 + O2 => CH2O + OH
    kiv[45] = {11,5,10,3};
    nuv[45] = {-1,-1,1,1};
    // (34):  CH3 + O2 => CH2O + OH
    fwd_A[45]     = 747000000000;
    fwd_beta[45]  = 0;
    fwd_Ea[45]    = 14250;
    prefactor_units[45]  = 1.0000000000000002e-06;
    activation_units[45] = 0.50321666580471969;
    phase_units[45]      = pow(10,-12.000000);
    is_PD[45] = 0;
    nTB[45] = 0;

    // (35):  CH3 + O2 => CH3O + O
    kiv[46] = {11,5,13,1};
    nuv[46] = {-1,-1,1,1};
    // (35):  CH3 + O2 => CH3O + O
    fwd_A[46]     = 1.995e+18;
    fwd_beta[46]  = -1.5700000000000001;
    fwd_Ea[46]    = 29210.09;
    prefactor_units[46]  = 1.0000000000000002e-06;
    activation_units[46] = 0.50321666580471969;
    phase_units[46]      = pow(10,-12.000000);
    is_PD[46] = 0;
    nTB[46] = 0;

    // (36):  CH3O + O => CH3 + O2
    kiv[47] = {13,1,11,5};
    nuv[47] = {-1,-1,1,1};
    // (36):  CH3O + O => CH3 + O2
    fwd_A[47]     = 3.585e+18;
    fwd_beta[47]  = -1.5900000000000001;
    fwd_Ea[47]    = -1630.98;
    prefactor_units[47]  = 1.0000000000000002e-06;
    activation_units[47] = 0.50321666580471969;
    phase_units[47]      = pow(10,-12.000000);
    is_PD[47] = 0;
    nTB[47] = 0;

    // (37):  2.000000 CH3 <=> H + C2H5
    kiv[48] = {11,2,14};
    nuv[48] = {-2.0,1,1};
    // (37):  2.000000 CH3 <=> H + C2H5
    fwd_A[48]     = 6840000000000;
    fwd_beta[48]  = 0.10000000000000001;
    fwd_Ea[48]    = 10599.9;
    prefactor_units[48]  = 1.0000000000000002e-06;
    activation_units[48] = 0.50321666580471969;
    phase_units[48]      = pow(10,-12.000000);
    is_PD[48] = 0;
    nTB[48] = 0;

    // (38):  CH3 + HO2 => CH4 + O2
    kiv[49] = {11,8,15,5};
    nuv[49] = {-1,-1,1,1};
    // (38):  CH3 + HO2 => CH4 + O2
    fwd_A[49]     = 3600000000000;
    fwd_beta[49]  = 0;
    fwd_Ea[49]    = 0;
    prefactor_units[49]  = 1.0000000000000002e-06;
    activation_units[49] = 0.50321666580471969;
    phase_units[49]      = pow(10,-12.000000);
    is_PD[49] = 0;
    nTB[49] = 0;

    // (39):  CH3 + O => CH2O + H
    kiv[50] = {11,1,10,2};
    nuv[50] = {-1,-1,1,1};
    // (39):  CH3 + O => CH2O + H
    fwd_A[50]     = 80000000000000;
    fwd_beta[50]  = 0;
    fwd_Ea[50]    = 0;
    prefactor_units[50]  = 1.0000000000000002e-06;
    activation_units[50] = 0.50321666580471969;
    phase_units[50]      = pow(10,-12.000000);
    is_PD[50] = 0;
    nTB[50] = 0;

    // (40):  CH3 + OH => CH2O + H2
    kiv[51] = {11,3,10,4};
    nuv[51] = {-1,-1,1,1};
    // (40):  CH3 + OH => CH2O + H2
    fwd_A[51]     = 22500000000000;
    fwd_beta[51]  = 0;
    fwd_Ea[51]    = 4299.9499999999998;
    prefactor_units[51]  = 1.0000000000000002e-06;
    activation_units[51] = 0.50321666580471969;
    phase_units[51]      = pow(10,-12.000000);
    is_PD[51] = 0;
    nTB[51] = 0;

    // (41):  CH3 + OH => CH2GSG + H2O
    kiv[52] = {11,3,9,6};
    nuv[52] = {-1,-1,1,1};
    // (41):  CH3 + OH => CH2GSG + H2O
    fwd_A[52]     = 26500000000000;
    fwd_beta[52]  = 0;
    fwd_Ea[52]    = 2185.9499999999998;
    prefactor_units[52]  = 1.0000000000000002e-06;
    activation_units[52] = 0.50321666580471969;
    phase_units[52]      = pow(10,-12.000000);
    is_PD[52] = 0;
    nTB[52] = 0;

    // (42):  CH2GSG + H2O => CH3 + OH
    kiv[53] = {9,6,11,3};
    nuv[53] = {-1,-1,1,1};
    // (42):  CH2GSG + H2O => CH3 + OH
    fwd_A[53]     = 32360000000;
    fwd_beta[53]  = 0.89000000000000001;
    fwd_Ea[53]    = 1211.04;
    prefactor_units[53]  = 1.0000000000000002e-06;
    activation_units[53] = 0.50321666580471969;
    phase_units[53]      = pow(10,-12.000000);
    is_PD[53] = 0;
    nTB[53] = 0;

    // (43):  CH3 + H2O2 => CH4 + HO2
    kiv[54] = {11,7,15,8};
    nuv[54] = {-1,-1,1,1};
    // (43):  CH3 + H2O2 => CH4 + HO2
    fwd_A[54]     = 336500000000;
    fwd_beta[54]  = -0.33000000000000002;
    fwd_Ea[54]    = 2501.9099999999999;
    prefactor_units[54]  = 1.0000000000000002e-06;
    activation_units[54] = 0.50321666580471969;
    phase_units[54]      = pow(10,-12.000000);
    is_PD[54] = 0;
    nTB[54] = 0;

    // (44):  CH2GSG + CH3 => C2H4 + H
    kiv[55] = {9,11,16,2};
    nuv[55] = {-1,-1,1,1};
    // (44):  CH2GSG + CH3 => C2H4 + H
    fwd_A[55]     = 20000000000000;
    fwd_beta[55]  = 0;
    fwd_Ea[55]    = 0;
    prefactor_units[55]  = 1.0000000000000002e-06;
    activation_units[55] = 0.50321666580471969;
    phase_units[55]      = pow(10,-12.000000);
    is_PD[55] = 0;
    nTB[55] = 0;

    // (45):  CH3 + H (+M) => CH4 (+M)
    kiv[2] = {11,2,15};
    nuv[2] = {-1,-1,1};
    // (45):  CH3 + H (+M) => CH4 (+M)
    fwd_A[2]     = 2138000000000000;
    fwd_beta[2]  = -0.40000000000000002;
    fwd_Ea[2]    = 0;
    low_A[2]     = 3.31e+30;
    low_beta[2]  = -4;
    low_Ea[2]    = 2108.0300000000002;
    troe_a[2]    = 0;
    troe_Tsss[2] = 0;
    troe_Ts[2]   = 0;
    troe_Tss[2]  = 40;
    troe_len[2]  = 4;
    prefactor_units[2]  = 1.0000000000000002e-06;
    activation_units[2] = 0.50321666580471969;
    phase_units[2]      = pow(10,-12.000000);
    is_PD[2] = 1;
    nTB[2] = 4;
    TB[2] = (amrex::Real *) malloc(4 * sizeof(amrex::Real));
    TBid[2] = (int *) malloc(4 * sizeof(int));
    TBid[2][0] = 4; TB[2][0] = 2; // H2
    TBid[2][1] = 6; TB[2][1] = 5; // H2O
    TBid[2][2] = 18; TB[2][2] = 3; // CO2
    TBid[2][3] = 12; TB[2][3] = 2; // CO

    // (46):  2.000000 CH3 (+M) <=> C2H6 (+M)
    kiv[3] = {11,17};
    nuv[3] = {-2.0,1};
    // (46):  2.000000 CH3 (+M) <=> C2H6 (+M)
    fwd_A[3]     = 92140000000000000;
    fwd_beta[3]  = -1.1699999999999999;
    fwd_Ea[3]    = 635.75999999999999;
    low_A[3]     = 1.135e+36;
    low_beta[3]  = -5.2460000000000004;
    low_Ea[3]    = 1705.0699999999999;
    troe_a[3]    = 0.40500000000000003;
    troe_Tsss[3] = 1120;
    troe_Ts[3]   = 69.599999999999994;
    troe_Tss[3]  = 1000000000000000;
    troe_len[3]  = 4;
    prefactor_units[3]  = 1.0000000000000002e-06;
    activation_units[3] = 0.50321666580471969;
    phase_units[3]      = pow(10,-12.000000);
    is_PD[3] = 1;
    nTB[3] = 4;
    TB[3] = (amrex::Real *) malloc(4 * sizeof(amrex::Real));
    TBid[3] = (int *) malloc(4 * sizeof(int));
    TBid[3][0] = 4; TB[3][0] = 2; // H2
    TBid[3][1] = 6; TB[3][1] = 5; // H2O
    TBid[3][2] = 18; TB[3][2] = 3; // CO2
    TBid[3][3] = 12; TB[3][3] = 2; // CO

    // (47):  CH2GSG + CH4 => 2.000000 CH3
    kiv[56] = {9,15,11};
    nuv[56] = {-1,-1,2.0};
    // (47):  CH2GSG + CH4 => 2.000000 CH3
    fwd_A[56]     = 40000000000000;
    fwd_beta[56]  = 0;
    fwd_Ea[56]    = 0;
    prefactor_units[56]  = 1.0000000000000002e-06;
    activation_units[56] = 0.50321666580471969;
    phase_units[56]      = pow(10,-12.000000);
    is_PD[56] = 0;
    nTB[56] = 0;

    // (48):  2.000000 CH3 => CH2GSG + CH4
    kiv[57] = {11,9,15};
    nuv[57] = {-2.0,1,1};
    // (48):  2.000000 CH3 => CH2GSG + CH4
    fwd_A[57]     = 5429000000000000;
    fwd_beta[57]  = -0.89000000000000001;
    fwd_Ea[57]    = 15650.1;
    prefactor_units[57]  = 1.0000000000000002e-06;
    activation_units[57] = 0.50321666580471969;
    phase_units[57]      = pow(10,-12.000000);
    is_PD[57] = 0;
    nTB[57] = 0;

    // (49):  CH4 + O => CH3 + OH
    kiv[58] = {15,1,11,3};
    nuv[58] = {-1,-1,1,1};
    // (49):  CH4 + O => CH3 + OH
    fwd_A[58]     = 3150000000000;
    fwd_beta[58]  = 0.5;
    fwd_Ea[58]    = 10289.91;
    prefactor_units[58]  = 1.0000000000000002e-06;
    activation_units[58] = 0.50321666580471969;
    phase_units[58]      = pow(10,-12.000000);
    is_PD[58] = 0;
    nTB[58] = 0;

    // (50):  CH3 + OH => CH4 + O
    kiv[59] = {11,3,15,1};
    nuv[59] = {-1,-1,1,1};
    // (50):  CH3 + OH => CH4 + O
    fwd_A[59]     = 52960000000;
    fwd_beta[59]  = 0.5;
    fwd_Ea[59]    = 7715.1099999999997;
    prefactor_units[59]  = 1.0000000000000002e-06;
    activation_units[59] = 0.50321666580471969;
    phase_units[59]      = pow(10,-12.000000);
    is_PD[59] = 0;
    nTB[59] = 0;

    // (51):  CH4 + H => CH3 + H2
    kiv[60] = {15,2,11,4};
    nuv[60] = {-1,-1,1,1};
    // (51):  CH4 + H => CH3 + H2
    fwd_A[60]     = 17270;
    fwd_beta[60]  = 3;
    fwd_Ea[60]    = 8223.9500000000007;
    prefactor_units[60]  = 1.0000000000000002e-06;
    activation_units[60] = 0.50321666580471969;
    phase_units[60]      = pow(10,-12.000000);
    is_PD[60] = 0;
    nTB[60] = 0;

    // (52):  CH3 + H2 => CH4 + H
    kiv[61] = {11,4,15,2};
    nuv[61] = {-1,-1,1,1};
    // (52):  CH3 + H2 => CH4 + H
    fwd_A[61]     = 661;
    fwd_beta[61]  = 3;
    fwd_Ea[61]    = 7744.0200000000004;
    prefactor_units[61]  = 1.0000000000000002e-06;
    activation_units[61] = 0.50321666580471969;
    phase_units[61]      = pow(10,-12.000000);
    is_PD[61] = 0;
    nTB[61] = 0;

    // (53):  CH4 + OH => CH3 + H2O
    kiv[62] = {15,3,11,6};
    nuv[62] = {-1,-1,1,1};
    // (53):  CH4 + OH => CH3 + H2O
    fwd_A[62]     = 193000;
    fwd_beta[62]  = 2.3999999999999999;
    fwd_Ea[62]    = 2106.1199999999999;
    prefactor_units[62]  = 1.0000000000000002e-06;
    activation_units[62] = 0.50321666580471969;
    phase_units[62]      = pow(10,-12.000000);
    is_PD[62] = 0;
    nTB[62] = 0;

    // (54):  CH3 + H2O => CH4 + OH
    kiv[63] = {11,6,15,3};
    nuv[63] = {-1,-1,1,1};
    // (54):  CH3 + H2O => CH4 + OH
    fwd_A[63]     = 31990;
    fwd_beta[63]  = 2.3999999999999999;
    fwd_Ea[63]    = 16780.110000000001;
    prefactor_units[63]  = 1.0000000000000002e-06;
    activation_units[63] = 0.50321666580471969;
    phase_units[63]      = pow(10,-12.000000);
    is_PD[63] = 0;
    nTB[63] = 0;

    // (55):  CO + HO2 => CO2 + OH
    kiv[64] = {12,8,18,3};
    nuv[64] = {-1,-1,1,1};
    // (55):  CO + HO2 => CO2 + OH
    fwd_A[64]     = 30100000000000;
    fwd_beta[64]  = 0;
    fwd_Ea[64]    = 23000;
    prefactor_units[64]  = 1.0000000000000002e-06;
    activation_units[64] = 0.50321666580471969;
    phase_units[64]      = pow(10,-12.000000);
    is_PD[64] = 0;
    nTB[64] = 0;

    // (56):  CO + O2 => CO2 + O
    kiv[65] = {12,5,18,1};
    nuv[65] = {-1,-1,1,1};
    // (56):  CO + O2 => CO2 + O
    fwd_A[65]     = 16200000000000;
    fwd_beta[65]  = 0;
    fwd_Ea[65]    = 47700.050000000003;
    prefactor_units[65]  = 1.0000000000000002e-06;
    activation_units[65] = 0.50321666580471969;
    phase_units[65]      = pow(10,-12.000000);
    is_PD[65] = 0;
    nTB[65] = 0;

    // (57):  CO2 + O => CO + O2
    kiv[66] = {18,1,12,5};
    nuv[66] = {-1,-1,1,1};
    // (57):  CO2 + O => CO + O2
    fwd_A[66]     = 143300000000000;
    fwd_beta[66]  = 0;
    fwd_Ea[66]    = 53919.93;
    prefactor_units[66]  = 1.0000000000000002e-06;
    activation_units[66] = 0.50321666580471969;
    phase_units[66]      = pow(10,-12.000000);
    is_PD[66] = 0;
    nTB[66] = 0;

    // (58):  CO + OH => CO2 + H
    kiv[67] = {12,3,18,2};
    nuv[67] = {-1,-1,1,1};
    // (58):  CO + OH => CO2 + H
    fwd_A[67]     = 140000;
    fwd_beta[67]  = 1.95;
    fwd_Ea[67]    = -1347.04;
    prefactor_units[67]  = 1.0000000000000002e-06;
    activation_units[67] = 0.50321666580471969;
    phase_units[67]      = pow(10,-12.000000);
    is_PD[67] = 0;
    nTB[67] = 0;

    // (59):  CO2 + H => CO + OH
    kiv[68] = {18,2,12,3};
    nuv[68] = {-1,-1,1,1};
    // (59):  CO2 + H => CO + OH
    fwd_A[68]     = 15680000;
    fwd_beta[68]  = 1.95;
    fwd_Ea[68]    = 20989.959999999999;
    prefactor_units[68]  = 1.0000000000000002e-06;
    activation_units[68] = 0.50321666580471969;
    phase_units[68]      = pow(10,-12.000000);
    is_PD[68] = 0;
    nTB[68] = 0;

    // (60):  CO + O (+M) => CO2 (+M)
    kiv[4] = {12,1,18};
    nuv[4] = {-1,-1,1};
    // (60):  CO + O (+M) => CO2 (+M)
    fwd_A[4]     = 18000000000;
    fwd_beta[4]  = 0;
    fwd_Ea[4]    = 2384.0799999999999;
    low_A[4]     = 1.35e+24;
    low_beta[4]  = -2.7879999999999998;
    low_Ea[4]    = 4190.9700000000003;
    troe_a[4]    = 1;
    troe_Tsss[4] = 1;
    troe_Ts[4]   = 10000000;
    troe_Tss[4]  = 10000000;
    troe_len[4]  = 4;
    prefactor_units[4]  = 1.0000000000000002e-06;
    activation_units[4] = 0.50321666580471969;
    phase_units[4]      = pow(10,-12.000000);
    is_PD[4] = 1;
    nTB[4] = 4;
    TB[4] = (amrex::Real *) malloc(4 * sizeof(amrex::Real));
    TBid[4] = (int *) malloc(4 * sizeof(int));
    TBid[4][0] = 4; TB[4][0] = 2.5; // H2
    TBid[4][1] = 6; TB[4][1] = 12; // H2O
    TBid[4][2] = 18; TB[4][2] = 3.7999999999999998; // CO2
    TBid[4][3] = 12; TB[4][3] = 1.8999999999999999; // CO

    // (61):  HCO + CH3 => CH4 + CO
    kiv[69] = {19,11,15,12};
    nuv[69] = {-1,-1,1,1};
    // (61):  HCO + CH3 => CH4 + CO
    fwd_A[69]     = 121000000000000;
    fwd_beta[69]  = 0;
    fwd_Ea[69]    = 0;
    prefactor_units[69]  = 1.0000000000000002e-06;
    activation_units[69] = 0.50321666580471969;
    phase_units[69]      = pow(10,-12.000000);
    is_PD[69] = 0;
    nTB[69] = 0;

    // (62):  HCO + H => CO + H2
    kiv[70] = {19,2,12,4};
    nuv[70] = {-1,-1,1,1};
    // (62):  HCO + H => CO + H2
    fwd_A[70]     = 73400000000000;
    fwd_beta[70]  = 0;
    fwd_Ea[70]    = 0;
    prefactor_units[70]  = 1.0000000000000002e-06;
    activation_units[70] = 0.50321666580471969;
    phase_units[70]      = pow(10,-12.000000);
    is_PD[70] = 0;
    nTB[70] = 0;

    // (63):  HCO + O2 => CO + HO2
    kiv[71] = {19,5,12,8};
    nuv[71] = {-1,-1,1,1};
    // (63):  HCO + O2 => CO + HO2
    fwd_A[71]     = 7580000000000;
    fwd_beta[71]  = 0;
    fwd_Ea[71]    = 409.88999999999999;
    prefactor_units[71]  = 1.0000000000000002e-06;
    activation_units[71] = 0.50321666580471969;
    phase_units[71]      = pow(10,-12.000000);
    is_PD[71] = 0;
    nTB[71] = 0;

    // (64):  HCO + O => CO + OH
    kiv[72] = {19,1,12,3};
    nuv[72] = {-1,-1,1,1};
    // (64):  HCO + O => CO + OH
    fwd_A[72]     = 30200000000000;
    fwd_beta[72]  = 0;
    fwd_Ea[72]    = 0;
    prefactor_units[72]  = 1.0000000000000002e-06;
    activation_units[72] = 0.50321666580471969;
    phase_units[72]      = pow(10,-12.000000);
    is_PD[72] = 0;
    nTB[72] = 0;

    // (65):  HCO + O => CO2 + H
    kiv[73] = {19,1,18,2};
    nuv[73] = {-1,-1,1,1};
    // (65):  HCO + O => CO2 + H
    fwd_A[73]     = 30000000000000;
    fwd_beta[73]  = 0;
    fwd_Ea[73]    = 0;
    prefactor_units[73]  = 1.0000000000000002e-06;
    activation_units[73] = 0.50321666580471969;
    phase_units[73]      = pow(10,-12.000000);
    is_PD[73] = 0;
    nTB[73] = 0;

    // (66):  HCO + M => H + CO + M
    kiv[12] = {19,2,12};
    nuv[12] = {-1,1,1};
    // (66):  HCO + M => H + CO + M
    fwd_A[12]     = 1.86e+17;
    fwd_beta[12]  = -1;
    fwd_Ea[12]    = 17000;
    prefactor_units[12]  = 1.0000000000000002e-06;
    activation_units[12] = 0.50321666580471969;
    phase_units[12]      = pow(10,-6.000000);
    is_PD[12] = 0;
    nTB[12] = 4;
    TB[12] = (amrex::Real *) malloc(4 * sizeof(amrex::Real));
    TBid[12] = (int *) malloc(4 * sizeof(int));
    TBid[12][0] = 4; TB[12][0] = 2.5; // H2
    TBid[12][1] = 6; TB[12][1] = 6; // H2O
    TBid[12][2] = 18; TB[12][2] = 3.7999999999999998; // CO2
    TBid[12][3] = 12; TB[12][3] = 1.8999999999999999; // CO

    // (67):  H + CO + M => HCO + M
    kiv[13] = {2,12,19};
    nuv[13] = {-1,-1,1};
    // (67):  H + CO + M => HCO + M
    fwd_A[13]     = 64670000000000;
    fwd_beta[13]  = 0;
    fwd_Ea[13]    = -441.92000000000002;
    prefactor_units[13]  = 1.0000000000000002e-12;
    activation_units[13] = 0.50321666580471969;
    phase_units[13]      = pow(10,-12.000000);
    is_PD[13] = 0;
    nTB[13] = 4;
    TB[13] = (amrex::Real *) malloc(4 * sizeof(amrex::Real));
    TBid[13] = (int *) malloc(4 * sizeof(int));
    TBid[13][0] = 4; TB[13][0] = 2.5; // H2
    TBid[13][1] = 6; TB[13][1] = 6; // H2O
    TBid[13][2] = 18; TB[13][2] = 3.7999999999999998; // CO2
    TBid[13][3] = 12; TB[13][3] = 1.8999999999999999; // CO

    // (68):  HCO + OH => CO + H2O
    kiv[74] = {19,3,12,6};
    nuv[74] = {-1,-1,1,1};
    // (68):  HCO + OH => CO + H2O
    fwd_A[74]     = 102000000000000;
    fwd_beta[74]  = 0;
    fwd_Ea[74]    = 0;
    prefactor_units[74]  = 1.0000000000000002e-06;
    activation_units[74] = 0.50321666580471969;
    phase_units[74]      = pow(10,-12.000000);
    is_PD[74] = 0;
    nTB[74] = 0;

    // (69):  CH2O + H => HCO + H2
    kiv[75] = {10,2,19,4};
    nuv[75] = {-1,-1,1,1};
    // (69):  CH2O + H => HCO + H2
    fwd_A[75]     = 933400000;
    fwd_beta[75]  = 1.5;
    fwd_Ea[75]    = 2976.0999999999999;
    prefactor_units[75]  = 1.0000000000000002e-06;
    activation_units[75] = 0.50321666580471969;
    phase_units[75]      = pow(10,-12.000000);
    is_PD[75] = 0;
    nTB[75] = 0;

    // (70):  CH2O + O2 => HCO + HO2
    kiv[76] = {10,5,19,8};
    nuv[76] = {-1,-1,1,1};
    // (70):  CH2O + O2 => HCO + HO2
    fwd_A[76]     = 20500000000000;
    fwd_beta[76]  = 0;
    fwd_Ea[76]    = 38950.050000000003;
    prefactor_units[76]  = 1.0000000000000002e-06;
    activation_units[76] = 0.50321666580471969;
    phase_units[76]      = pow(10,-12.000000);
    is_PD[76] = 0;
    nTB[76] = 0;

    // (71):  CH2O + OH => HCO + H2O
    kiv[77] = {10,3,19,6};
    nuv[77] = {-1,-1,1,1};
    // (71):  CH2O + OH => HCO + H2O
    fwd_A[77]     = 3430000000;
    fwd_beta[77]  = 1.1799999999999999;
    fwd_Ea[77]    = -446.94;
    prefactor_units[77]  = 1.0000000000000002e-06;
    activation_units[77] = 0.50321666580471969;
    phase_units[77]      = pow(10,-12.000000);
    is_PD[77] = 0;
    nTB[77] = 0;

    // (72):  CH2O + HO2 => HCO + H2O2
    kiv[78] = {10,8,19,7};
    nuv[78] = {-1,-1,1,1};
    // (72):  CH2O + HO2 => HCO + H2O2
    fwd_A[78]     = 0.0058199999999999997;
    fwd_beta[78]  = 4.5300000000000002;
    fwd_Ea[78]    = 6556.8800000000001;
    prefactor_units[78]  = 1.0000000000000002e-06;
    activation_units[78] = 0.50321666580471969;
    phase_units[78]      = pow(10,-12.000000);
    is_PD[78] = 0;
    nTB[78] = 0;

    // (73):  CH2O + O => HCO + OH
    kiv[79] = {10,1,19,3};
    nuv[79] = {-1,-1,1,1};
    // (73):  CH2O + O => HCO + OH
    fwd_A[79]     = 416000000000;
    fwd_beta[79]  = 0.56999999999999995;
    fwd_Ea[79]    = 2761.9499999999998;
    prefactor_units[79]  = 1.0000000000000002e-06;
    activation_units[79] = 0.50321666580471969;
    phase_units[79]      = pow(10,-12.000000);
    is_PD[79] = 0;
    nTB[79] = 0;

    // (74):  CH2O + CH3 => HCO + CH4
    kiv[80] = {10,11,19,15};
    nuv[80] = {-1,-1,1,1};
    // (74):  CH2O + CH3 => HCO + CH4
    fwd_A[80]     = 3.636e-06;
    fwd_beta[80]  = 5.4199999999999999;
    fwd_Ea[80]    = 998.09000000000003;
    prefactor_units[80]  = 1.0000000000000002e-06;
    activation_units[80] = 0.50321666580471969;
    phase_units[80]      = pow(10,-12.000000);
    is_PD[80] = 0;
    nTB[80] = 0;

    // (75):  CH3O (+M) => CH2O + H (+M)
    kiv[5] = {13,10,2};
    nuv[5] = {-1,1,1};
    // (75):  CH3O (+M) => CH2O + H (+M)
    fwd_A[5]     = 54500000000000;
    fwd_beta[5]  = 0;
    fwd_Ea[5]    = 13500;
    low_A[5]     = 2.344e+25;
    low_beta[5]  = -2.7000000000000002;
    low_Ea[5]    = 30599.900000000001;
    troe_a[5]    = 1;
    troe_Tsss[5] = 1;
    troe_Ts[5]   = 10000000;
    troe_Tss[5]  = 10000000;
    troe_len[5]  = 4;
    prefactor_units[5]  = 1;
    activation_units[5] = 0.50321666580471969;
    phase_units[5]      = pow(10,-6.000000);
    is_PD[5] = 1;
    nTB[5] = 0;

    // (76):  CH3O + O2 => CH2O + HO2
    kiv[81] = {13,5,10,8};
    nuv[81] = {-1,-1,1,1};
    // (76):  CH3O + O2 => CH2O + HO2
    fwd_A[81]     = 55000000000;
    fwd_beta[81]  = 0;
    fwd_Ea[81]    = 2424;
    prefactor_units[81]  = 1.0000000000000002e-06;
    activation_units[81] = 0.50321666580471969;
    phase_units[81]      = pow(10,-12.000000);
    is_PD[81] = 0;
    nTB[81] = 0;

    // (77):  CH2GSG + CO2 => CH2O + CO
    kiv[82] = {9,18,10,12};
    nuv[82] = {-1,-1,1,1};
    // (77):  CH2GSG + CO2 => CH2O + CO
    fwd_A[82]     = 3000000000000;
    fwd_beta[82]  = 0;
    fwd_Ea[82]    = 0;
    prefactor_units[82]  = 1.0000000000000002e-06;
    activation_units[82] = 0.50321666580471969;
    phase_units[82]      = pow(10,-12.000000);
    is_PD[82] = 0;
    nTB[82] = 0;

    // (78):  CH3O2 + CH3 => 2.000000 CH3O
    kiv[83] = {20,11,13};
    nuv[83] = {-1,-1,2.0};
    // (78):  CH3O2 + CH3 => 2.000000 CH3O
    fwd_A[83]     = 7000000000000;
    fwd_beta[83]  = 0;
    fwd_Ea[83]    = -1000;
    prefactor_units[83]  = 1.0000000000000002e-06;
    activation_units[83] = 0.50321666580471969;
    phase_units[83]      = pow(10,-12.000000);
    is_PD[83] = 0;
    nTB[83] = 0;

    // (79):  CH3O2 + HO2 => CH3O2H + O2
    kiv[84] = {20,8,21,5};
    nuv[84] = {-1,-1,1,1};
    // (79):  CH3O2 + HO2 => CH3O2H + O2
    fwd_A[84]     = 17500000000;
    fwd_beta[84]  = 0;
    fwd_Ea[84]    = -3275.0999999999999;
    prefactor_units[84]  = 1.0000000000000002e-06;
    activation_units[84] = 0.50321666580471969;
    phase_units[84]      = pow(10,-12.000000);
    is_PD[84] = 0;
    nTB[84] = 0;

    // (80):  2.000000 CH3O2 => O2 + 2.000000 CH3O
    kiv[85] = {20,5,13};
    nuv[85] = {-2.0,1,2.0};
    // (80):  2.000000 CH3O2 => O2 + 2.000000 CH3O
    fwd_A[85]     = 14000000000000000;
    fwd_beta[85]  = -1.6100000000000001;
    fwd_Ea[85]    = 1859.9400000000001;
    prefactor_units[85]  = 1.0000000000000002e-06;
    activation_units[85] = 0.50321666580471969;
    phase_units[85]      = pow(10,-12.000000);
    is_PD[85] = 0;
    nTB[85] = 0;

    // (81):  CH3O2 + CH2O => CH3O2H + HCO
    kiv[86] = {20,10,21,19};
    nuv[86] = {-1,-1,1,1};
    // (81):  CH3O2 + CH2O => CH3O2H + HCO
    fwd_A[86]     = 1990000000000;
    fwd_beta[86]  = 0;
    fwd_Ea[86]    = 11659.889999999999;
    prefactor_units[86]  = 1.0000000000000002e-06;
    activation_units[86] = 0.50321666580471969;
    phase_units[86]      = pow(10,-12.000000);
    is_PD[86] = 0;
    nTB[86] = 0;

    // (82):  CH3O2 + M => CH3 + O2 + M
    kiv[14] = {20,11,5};
    nuv[14] = {-1,1,1};
    // (82):  CH3O2 + M => CH3 + O2 + M
    fwd_A[14]     = 4.3430000000000002e+27;
    fwd_beta[14]  = -3.4199999999999999;
    fwd_Ea[14]    = 30469.889999999999;
    prefactor_units[14]  = 1.0000000000000002e-06;
    activation_units[14] = 0.50321666580471969;
    phase_units[14]      = pow(10,-6.000000);
    is_PD[14] = 0;
    nTB[14] = 0;

    // (83):  CH3 + O2 + M => CH3O2 + M
    kiv[15] = {11,5,20};
    nuv[15] = {-1,-1,1};
    // (83):  CH3 + O2 + M => CH3O2 + M
    fwd_A[15]     = 5.4400000000000001e+25;
    fwd_beta[15]  = -3.2999999999999998;
    fwd_Ea[15]    = 0;
    prefactor_units[15]  = 1.0000000000000002e-12;
    activation_units[15] = 0.50321666580471969;
    phase_units[15]      = pow(10,-12.000000);
    is_PD[15] = 0;
    nTB[15] = 0;

    // (84):  CH3O2H => CH3O + OH
    kiv[87] = {21,13,3};
    nuv[87] = {-1,1,1};
    // (84):  CH3O2H => CH3O + OH
    fwd_A[87]     = 631000000000000;
    fwd_beta[87]  = 0;
    fwd_Ea[87]    = 42299.949999999997;
    prefactor_units[87]  = 1;
    activation_units[87] = 0.50321666580471969;
    phase_units[87]      = pow(10,-6.000000);
    is_PD[87] = 0;
    nTB[87] = 0;

    // (85):  CH3O + OH => CH3O2H
    kiv[88] = {13,3,21};
    nuv[88] = {-1,-1,1};
    // (85):  CH3O + OH => CH3O2H
    fwd_A[88]     = 116600000000;
    fwd_beta[88]  = 0.59999999999999998;
    fwd_Ea[88]    = -1771.03;
    prefactor_units[88]  = 1.0000000000000002e-06;
    activation_units[88] = 0.50321666580471969;
    phase_units[88]      = pow(10,-12.000000);
    is_PD[88] = 0;
    nTB[88] = 0;

    // (86):  C2H2 + O2 => HCCO + OH
    kiv[89] = {22,5,23,3};
    nuv[89] = {-1,-1,1,1};
    // (86):  C2H2 + O2 => HCCO + OH
    fwd_A[89]     = 200000000;
    fwd_beta[89]  = 1.5;
    fwd_Ea[89]    = 30099.900000000001;
    prefactor_units[89]  = 1.0000000000000002e-06;
    activation_units[89] = 0.50321666580471969;
    phase_units[89]      = pow(10,-12.000000);
    is_PD[89] = 0;
    nTB[89] = 0;

    // (87):  C2H2 + O => HCCO + H
    kiv[90] = {22,1,23,2};
    nuv[90] = {-1,-1,1,1};
    // (87):  C2H2 + O => HCCO + H
    fwd_A[90]     = 14300000;
    fwd_beta[90]  = 2;
    fwd_Ea[90]    = 1900.0999999999999;
    prefactor_units[90]  = 1.0000000000000002e-06;
    activation_units[90] = 0.50321666580471969;
    phase_units[90]      = pow(10,-12.000000);
    is_PD[90] = 0;
    nTB[90] = 0;

    // (88):  C2H3 + O2 => CH2CHO + O
    kiv[91] = {24,5,25,1};
    nuv[91] = {-1,-1,1,1};
    // (88):  C2H3 + O2 => CH2CHO + O
    fwd_A[91]     = 350000000000000;
    fwd_beta[91]  = -0.60999999999999999;
    fwd_Ea[91]    = 5260.04;
    prefactor_units[91]  = 1.0000000000000002e-06;
    activation_units[91] = 0.50321666580471969;
    phase_units[91]      = pow(10,-12.000000);
    is_PD[91] = 0;
    nTB[91] = 0;

    // (89):  C2H3 + H => C2H2 + H2
    kiv[92] = {24,2,22,4};
    nuv[92] = {-1,-1,1,1};
    // (89):  C2H3 + H => C2H2 + H2
    fwd_A[92]     = 20000000000000;
    fwd_beta[92]  = 0;
    fwd_Ea[92]    = 2500;
    prefactor_units[92]  = 1.0000000000000002e-06;
    activation_units[92] = 0.50321666580471969;
    phase_units[92]      = pow(10,-12.000000);
    is_PD[92] = 0;
    nTB[92] = 0;

    // (90):  C2H3 + O2 => CH2O + HCO
    kiv[93] = {24,5,10,19};
    nuv[93] = {-1,-1,1,1};
    // (90):  C2H3 + O2 => CH2O + HCO
    fwd_A[93]     = 1.6999999999999999e+29;
    fwd_beta[93]  = -5.3099999999999996;
    fwd_Ea[93]    = 6500;
    prefactor_units[93]  = 1.0000000000000002e-06;
    activation_units[93] = 0.50321666580471969;
    phase_units[93]      = pow(10,-12.000000);
    is_PD[93] = 0;
    nTB[93] = 0;

    // (91):  C2H3 + O2 => C2H2 + HO2
    kiv[94] = {24,5,22,8};
    nuv[94] = {-1,-1,1,1};
    // (91):  C2H3 + O2 => C2H2 + HO2
    fwd_A[94]     = 2.12e-06;
    fwd_beta[94]  = 6;
    fwd_Ea[94]    = 9483.9899999999998;
    prefactor_units[94]  = 1.0000000000000002e-06;
    activation_units[94] = 0.50321666580471969;
    phase_units[94]      = pow(10,-12.000000);
    is_PD[94] = 0;
    nTB[94] = 0;

    // (92):  C2H3 (+M) => H + C2H2 (+M)
    kiv[6] = {24,2,22};
    nuv[6] = {-1,1,1};
    // (92):  C2H3 (+M) => H + C2H2 (+M)
    fwd_A[6]     = 16060000000;
    fwd_beta[6]  = 1.028;
    fwd_Ea[6]    = 40503.589999999997;
    low_A[6]     = 1.164e+39;
    low_beta[6]  = -6.8209999999999997;
    low_Ea[6]    = 44491.629999999997;
    troe_a[6]    = 1;
    troe_Tsss[6] = 0;
    troe_Ts[6]   = 675;
    troe_Tss[6]  = 1000000000000000;
    troe_len[6]  = 4;
    prefactor_units[6]  = 1;
    activation_units[6] = 0.50321666580471969;
    phase_units[6]      = pow(10,-6.000000);
    is_PD[6] = 1;
    nTB[6] = 4;
    TB[6] = (amrex::Real *) malloc(4 * sizeof(amrex::Real));
    TBid[6] = (int *) malloc(4 * sizeof(int));
    TBid[6][0] = 4; TB[6][0] = 2; // H2
    TBid[6][1] = 6; TB[6][1] = 5; // H2O
    TBid[6][2] = 18; TB[6][2] = 3; // CO2
    TBid[6][3] = 12; TB[6][3] = 2; // CO

    // (93):  C2H4 + O => CH3 + HCO
    kiv[95] = {16,1,11,19};
    nuv[95] = {-1,-1,1,1};
    // (93):  C2H4 + O => CH3 + HCO
    fwd_A[95]     = 10200000;
    fwd_beta[95]  = 1.8799999999999999;
    fwd_Ea[95]    = 179.02000000000001;
    prefactor_units[95]  = 1.0000000000000002e-06;
    activation_units[95] = 0.50321666580471969;
    phase_units[95]      = pow(10,-12.000000);
    is_PD[95] = 0;
    nTB[95] = 0;

    // (94):  H + C2H4 (+M) <=> C2H5 (+M)
    kiv[7] = {2,16,14};
    nuv[7] = {-1,-1,1};
    // (94):  H + C2H4 (+M) <=> C2H5 (+M)
    fwd_A[7]     = 1081000000000;
    fwd_beta[7]  = 0.45000000000000001;
    fwd_Ea[7]    = 1821.9400000000001;
    low_A[7]     = 1.1120000000000001e+34;
    low_beta[7]  = -5;
    low_Ea[7]    = 4447.8999999999996;
    troe_a[7]    = 1;
    troe_Tsss[7] = 0;
    troe_Ts[7]   = 95;
    troe_Tss[7]  = 200;
    troe_len[7]  = 4;
    prefactor_units[7]  = 1.0000000000000002e-06;
    activation_units[7] = 0.50321666580471969;
    phase_units[7]      = pow(10,-12.000000);
    is_PD[7] = 1;
    nTB[7] = 0;

    // (95):  C2H4 + OH => C2H3 + H2O
    kiv[96] = {16,3,24,6};
    nuv[96] = {-1,-1,1,1};
    // (95):  C2H4 + OH => C2H3 + H2O
    fwd_A[96]     = 20500000000000;
    fwd_beta[96]  = 0;
    fwd_Ea[96]    = 5950.0500000000002;
    prefactor_units[96]  = 1.0000000000000002e-06;
    activation_units[96] = 0.50321666580471969;
    phase_units[96]      = pow(10,-12.000000);
    is_PD[96] = 0;
    nTB[96] = 0;

    // (96):  C2H3 + H2O => C2H4 + OH
    kiv[97] = {24,6,16,3};
    nuv[97] = {-1,-1,1,1};
    // (96):  C2H3 + H2O => C2H4 + OH
    fwd_A[97]     = 6033000000000000;
    fwd_beta[97]  = -0.82999999999999996;
    fwd_Ea[97]    = 21760.040000000001;
    prefactor_units[97]  = 1.0000000000000002e-06;
    activation_units[97] = 0.50321666580471969;
    phase_units[97]      = pow(10,-12.000000);
    is_PD[97] = 0;
    nTB[97] = 0;

    // (97):  C2H4 + H => C2H3 + H2
    kiv[98] = {16,2,24,4};
    nuv[98] = {-1,-1,1,1};
    // (97):  C2H4 + H => C2H3 + H2
    fwd_A[98]     = 0.0084200000000000004;
    fwd_beta[98]  = 4.6200000000000001;
    fwd_Ea[98]    = 2582.9299999999998;
    prefactor_units[98]  = 1.0000000000000002e-06;
    activation_units[98] = 0.50321666580471969;
    phase_units[98]      = pow(10,-12.000000);
    is_PD[98] = 0;
    nTB[98] = 0;

    // (98):  C2H3 + H2 => C2H4 + H
    kiv[99] = {24,4,16,2};
    nuv[99] = {-1,-1,1,1};
    // (98):  C2H3 + H2 => C2H4 + H
    fwd_A[99]     = 0.57230000000000003;
    fwd_beta[99]  = 3.79;
    fwd_Ea[99]    = 3233.0300000000002;
    prefactor_units[99]  = 1.0000000000000002e-06;
    activation_units[99] = 0.50321666580471969;
    phase_units[99]      = pow(10,-12.000000);
    is_PD[99] = 0;
    nTB[99] = 0;

    // (99):  C2H4 + CH3O2 => C2H3 + CH3O2H
    kiv[100] = {16,20,24,21};
    nuv[100] = {-1,-1,1,1};
    // (99):  C2H4 + CH3O2 => C2H3 + CH3O2H
    fwd_A[100]     = 2230000000000;
    fwd_beta[100]  = 0;
    fwd_Ea[100]    = 17190.009999999998;
    prefactor_units[100]  = 1.0000000000000002e-06;
    activation_units[100] = 0.50321666580471969;
    phase_units[100]      = pow(10,-12.000000);
    is_PD[100] = 0;
    nTB[100] = 0;

    // (100):  C2H4 + CH3 => C2H3 + CH4
    kiv[101] = {16,11,24,15};
    nuv[101] = {-1,-1,1,1};
    // (100):  C2H4 + CH3 => C2H3 + CH4
    fwd_A[101]     = 6.6200000000000001;
    fwd_beta[101]  = 3.7000000000000002;
    fwd_Ea[101]    = 9500;
    prefactor_units[101]  = 1.0000000000000002e-06;
    activation_units[101] = 0.50321666580471969;
    phase_units[101]      = pow(10,-12.000000);
    is_PD[101] = 0;
    nTB[101] = 0;

    // (101):  C2H3 + CH4 => C2H4 + CH3
    kiv[102] = {24,15,16,11};
    nuv[102] = {-1,-1,1,1};
    // (101):  C2H3 + CH4 => C2H4 + CH3
    fwd_A[102]     = 1.4399999999999999;
    fwd_beta[102]  = 4.0199999999999996;
    fwd_Ea[102]    = 5472.04;
    prefactor_units[102]  = 1.0000000000000002e-06;
    activation_units[102] = 0.50321666580471969;
    phase_units[102]      = pow(10,-12.000000);
    is_PD[102] = 0;
    nTB[102] = 0;

    // (102):  C2H4 + O => CH2CHO + H
    kiv[103] = {16,1,25,2};
    nuv[103] = {-1,-1,1,1};
    // (102):  C2H4 + O => CH2CHO + H
    fwd_A[103]     = 3390000;
    fwd_beta[103]  = 1.8799999999999999;
    fwd_Ea[103]    = 179.02000000000001;
    prefactor_units[103]  = 1.0000000000000002e-06;
    activation_units[103] = 0.50321666580471969;
    phase_units[103]      = pow(10,-12.000000);
    is_PD[103] = 0;
    nTB[103] = 0;

    // (103):  C2H4 (+M) => C2H2 + H2 (+M)
    kiv[8] = {16,22,4};
    nuv[8] = {-1,1,1};
    // (103):  C2H4 (+M) => C2H2 + H2 (+M)
    fwd_A[8]     = 18000000000000;
    fwd_beta[8]  = 0;
    fwd_Ea[8]    = 76000;
    low_A[8]     = 1500000000000000;
    low_beta[8]  = 0;
    low_Ea[8]    = 55440.010000000002;
    troe_a[8]    = 1;
    troe_Tsss[8] = 1;
    troe_Ts[8]   = 10000000;
    troe_Tss[8]  = 10000000;
    troe_len[8]  = 4;
    prefactor_units[8]  = 1;
    activation_units[8] = 0.50321666580471969;
    phase_units[8]      = pow(10,-6.000000);
    is_PD[8] = 1;
    nTB[8] = 0;

    // (104):  CH3 + C2H5 => CH4 + C2H4
    kiv[104] = {11,14,15,16};
    nuv[104] = {-1,-1,1,1};
    // (104):  CH3 + C2H5 => CH4 + C2H4
    fwd_A[104]     = 19500000000000;
    fwd_beta[104]  = -0.5;
    fwd_Ea[104]    = 0;
    prefactor_units[104]  = 1.0000000000000002e-06;
    activation_units[104] = 0.50321666580471969;
    phase_units[104]      = pow(10,-12.000000);
    is_PD[104] = 0;
    nTB[104] = 0;

    // (105):  C2H5 + O2 => C2H4 + HO2
    kiv[105] = {14,5,16,8};
    nuv[105] = {-1,-1,1,1};
    // (105):  C2H5 + O2 => C2H4 + HO2
    fwd_A[105]     = 1.22e+30;
    fwd_beta[105]  = -5.7599999999999998;
    fwd_Ea[105]    = 10099.9;
    prefactor_units[105]  = 1.0000000000000002e-06;
    activation_units[105] = 0.50321666580471969;
    phase_units[105]      = pow(10,-12.000000);
    is_PD[105] = 0;
    nTB[105] = 0;

    // (106):  C2H4 + HO2 => C2H5 + O2
    kiv[106] = {16,8,14,5};
    nuv[106] = {-1,-1,1,1};
    // (106):  C2H4 + HO2 => C2H5 + O2
    fwd_A[106]     = 1.259e+30;
    fwd_beta[106]  = -5.6299999999999999;
    fwd_Ea[106]    = 22299.950000000001;
    prefactor_units[106]  = 1.0000000000000002e-06;
    activation_units[106] = 0.50321666580471969;
    phase_units[106]      = pow(10,-12.000000);
    is_PD[106] = 0;
    nTB[106] = 0;

    // (107):  C2H5 + HO2 => C2H6 + O2
    kiv[107] = {14,8,17,5};
    nuv[107] = {-1,-1,1,1};
    // (107):  C2H5 + HO2 => C2H6 + O2
    fwd_A[107]     = 267900000;
    fwd_beta[107]  = 0.89000000000000001;
    fwd_Ea[107]    = -1922.0799999999999;
    prefactor_units[107]  = 1.0000000000000002e-06;
    activation_units[107] = 0.50321666580471969;
    phase_units[107]      = pow(10,-12.000000);
    is_PD[107] = 0;
    nTB[107] = 0;

    // (108):  H + C2H5 => C2H6
    kiv[108] = {2,14,17};
    nuv[108] = {-1,-1,1};
    // (108):  H + C2H5 => C2H6
    fwd_A[108]     = 583100000000;
    fwd_beta[108]  = 0.59899999999999998;
    fwd_Ea[108]    = -2913;
    prefactor_units[108]  = 1.0000000000000002e-06;
    activation_units[108] = 0.50321666580471969;
    phase_units[108]      = pow(10,-12.000000);
    is_PD[108] = 0;
    nTB[108] = 0;

    // (109):  C2H6 + H => C2H5 + H2
    kiv[109] = {17,2,14,4};
    nuv[109] = {-1,-1,1,1};
    // (109):  C2H6 + H => C2H5 + H2
    fwd_A[109]     = 554;
    fwd_beta[109]  = 3.5;
    fwd_Ea[109]    = 5167.0699999999997;
    prefactor_units[109]  = 1.0000000000000002e-06;
    activation_units[109] = 0.50321666580471969;
    phase_units[109]      = pow(10,-12.000000);
    is_PD[109] = 0;
    nTB[109] = 0;

    // (110):  C2H5 + H2 => C2H6 + H
    kiv[110] = {14,4,17,2};
    nuv[110] = {-1,-1,1,1};
    // (110):  C2H5 + H2 => C2H6 + H
    fwd_A[110]     = 0.13550000000000001;
    fwd_beta[110]  = 4.0599999999999996;
    fwd_Ea[110]    = 8857.0699999999997;
    prefactor_units[110]  = 1.0000000000000002e-06;
    activation_units[110] = 0.50321666580471969;
    phase_units[110]      = pow(10,-12.000000);
    is_PD[110] = 0;
    nTB[110] = 0;

    // (111):  C2H6 + OH => C2H5 + H2O
    kiv[111] = {17,3,14,6};
    nuv[111] = {-1,-1,1,1};
    // (111):  C2H6 + OH => C2H5 + H2O
    fwd_A[111]     = 58000000;
    fwd_beta[111]  = 1.73;
    fwd_Ea[111]    = 1159.8900000000001;
    prefactor_units[111]  = 1.0000000000000002e-06;
    activation_units[111] = 0.50321666580471969;
    phase_units[111]      = pow(10,-12.000000);
    is_PD[111] = 0;
    nTB[111] = 0;

    // (112):  CH2GSG + C2H6 => CH3 + C2H5
    kiv[112] = {9,17,11,14};
    nuv[112] = {-1,-1,1,1};
    // (112):  CH2GSG + C2H6 => CH3 + C2H5
    fwd_A[112]     = 120000000000000;
    fwd_beta[112]  = 0;
    fwd_Ea[112]    = 0;
    prefactor_units[112]  = 1.0000000000000002e-06;
    activation_units[112] = 0.50321666580471969;
    phase_units[112]      = pow(10,-12.000000);
    is_PD[112] = 0;
    nTB[112] = 0;

    // (113):  C2H6 + O => C2H5 + OH
    kiv[113] = {17,1,14,3};
    nuv[113] = {-1,-1,1,1};
    // (113):  C2H6 + O => C2H5 + OH
    fwd_A[113]     = 13000000;
    fwd_beta[113]  = 2.1299999999999999;
    fwd_Ea[113]    = 5190.0100000000002;
    prefactor_units[113]  = 1.0000000000000002e-06;
    activation_units[113] = 0.50321666580471969;
    phase_units[113]      = pow(10,-12.000000);
    is_PD[113] = 0;
    nTB[113] = 0;

    // (114):  C2H6 + CH3 => C2H5 + CH4
    kiv[114] = {17,11,14,15};
    nuv[114] = {-1,-1,1,1};
    // (114):  C2H6 + CH3 => C2H5 + CH4
    fwd_A[114]     = 1.5099999999999999e-07;
    fwd_beta[114]  = 6;
    fwd_Ea[114]    = 6047.0799999999999;
    prefactor_units[114]  = 1.0000000000000002e-06;
    activation_units[114] = 0.50321666580471969;
    phase_units[114]      = pow(10,-12.000000);
    is_PD[114] = 0;
    nTB[114] = 0;

    // (115):  HCCO + O => H + 2.000000 CO
    kiv[115] = {23,1,2,12};
    nuv[115] = {-1,-1,1,2.0};
    // (115):  HCCO + O => H + 2.000000 CO
    fwd_A[115]     = 80000000000000;
    fwd_beta[115]  = 0;
    fwd_Ea[115]    = 0;
    prefactor_units[115]  = 1.0000000000000002e-06;
    activation_units[115] = 0.50321666580471969;
    phase_units[115]      = pow(10,-12.000000);
    is_PD[115] = 0;
    nTB[115] = 0;

    // (116):  HCCO + O2 => CO2 + HCO
    kiv[116] = {23,5,18,19};
    nuv[116] = {-1,-1,1,1};
    // (116):  HCCO + O2 => CO2 + HCO
    fwd_A[116]     = 240000000000;
    fwd_beta[116]  = 0;
    fwd_Ea[116]    = -853.97000000000003;
    prefactor_units[116]  = 1.0000000000000002e-06;
    activation_units[116] = 0.50321666580471969;
    phase_units[116]      = pow(10,-12.000000);
    is_PD[116] = 0;
    nTB[116] = 0;

    // (117):  HCCO + OH => 2.000000 HCO
    kiv[117] = {23,3,19};
    nuv[117] = {-1,-1,2.0};
    // (117):  HCCO + OH => 2.000000 HCO
    fwd_A[117]     = 10000000000000;
    fwd_beta[117]  = 0;
    fwd_Ea[117]    = 0;
    prefactor_units[117]  = 1.0000000000000002e-06;
    activation_units[117] = 0.50321666580471969;
    phase_units[117]      = pow(10,-12.000000);
    is_PD[117] = 0;
    nTB[117] = 0;

    // (118):  HCCO + H => CH2GSG + CO
    kiv[118] = {23,2,9,12};
    nuv[118] = {-1,-1,1,1};
    // (118):  HCCO + H => CH2GSG + CO
    fwd_A[118]     = 110000000000000;
    fwd_beta[118]  = 0;
    fwd_Ea[118]    = 0;
    prefactor_units[118]  = 1.0000000000000002e-06;
    activation_units[118] = 0.50321666580471969;
    phase_units[118]      = pow(10,-12.000000);
    is_PD[118] = 0;
    nTB[118] = 0;

    // (119):  CH2GSG + CO => HCCO + H
    kiv[119] = {9,12,23,2};
    nuv[119] = {-1,-1,1,1};
    // (119):  CH2GSG + CO => HCCO + H
    fwd_A[119]     = 2046000000000;
    fwd_beta[119]  = 0.89000000000000001;
    fwd_Ea[119]    = 27830.07;
    prefactor_units[119]  = 1.0000000000000002e-06;
    activation_units[119] = 0.50321666580471969;
    phase_units[119]      = pow(10,-12.000000);
    is_PD[119] = 0;
    nTB[119] = 0;

    // (120):  CH2CHO + O2 => CH2O + CO + OH
    kiv[120] = {25,5,10,12,3};
    nuv[120] = {-1,-1,1,1,1};
    // (120):  CH2CHO + O2 => CH2O + CO + OH
    fwd_A[120]     = 20000000000000;
    fwd_beta[120]  = 0;
    fwd_Ea[120]    = 4200.0500000000002;
    prefactor_units[120]  = 1.0000000000000002e-06;
    activation_units[120] = 0.50321666580471969;
    phase_units[120]      = pow(10,-12.000000);
    is_PD[120] = 0;
    nTB[120] = 0;

    // (121):  C3H5XA + CH2O => C3H6 + HCO
    kiv[121] = {26,10,27,19};
    nuv[121] = {-1,-1,1,1};
    // (121):  C3H5XA + CH2O => C3H6 + HCO
    fwd_A[121]     = 630000000;
    fwd_beta[121]  = 1.8999999999999999;
    fwd_Ea[121]    = 18190.009999999998;
    prefactor_units[121]  = 1.0000000000000002e-06;
    activation_units[121] = 0.50321666580471969;
    phase_units[121]      = pow(10,-12.000000);
    is_PD[121] = 0;
    nTB[121] = 0;

    // (122):  C3H6 + HCO => C3H5XA + CH2O
    kiv[122] = {27,19,26,10};
    nuv[122] = {-1,-1,1,1};
    // (122):  C3H6 + HCO => C3H5XA + CH2O
    fwd_A[122]     = 109700000;
    fwd_beta[122]  = 1.8899999999999999;
    fwd_Ea[122]    = 15840.110000000001;
    prefactor_units[122]  = 1.0000000000000002e-06;
    activation_units[122] = 0.50321666580471969;
    phase_units[122]      = pow(10,-12.000000);
    is_PD[122] = 0;
    nTB[122] = 0;

    // (123):  C3H5XA + HO2 => C3H5O + OH
    kiv[123] = {26,8,28,3};
    nuv[123] = {-1,-1,1,1};
    // (123):  C3H5XA + HO2 => C3H5O + OH
    fwd_A[123]     = 7000000000000;
    fwd_beta[123]  = 0;
    fwd_Ea[123]    = -1000;
    prefactor_units[123]  = 1.0000000000000002e-06;
    activation_units[123] = 0.50321666580471969;
    phase_units[123]      = pow(10,-12.000000);
    is_PD[123] = 0;
    nTB[123] = 0;

    // (124):  C3H5XA + CH3O2 => C3H5O + CH3O
    kiv[124] = {26,20,28,13};
    nuv[124] = {-1,-1,1,1};
    // (124):  C3H5XA + CH3O2 => C3H5O + CH3O
    fwd_A[124]     = 7000000000000;
    fwd_beta[124]  = 0;
    fwd_Ea[124]    = -1000;
    prefactor_units[124]  = 1.0000000000000002e-06;
    activation_units[124] = 0.50321666580471969;
    phase_units[124]      = pow(10,-12.000000);
    is_PD[124] = 0;
    nTB[124] = 0;

    // (125):  C3H5XA + HO2 => C3H6 + O2
    kiv[125] = {26,8,27,5};
    nuv[125] = {-1,-1,1,1};
    // (125):  C3H5XA + HO2 => C3H6 + O2
    fwd_A[125]     = 33320000000;
    fwd_beta[125]  = 0.34000000000000002;
    fwd_Ea[125]    = -555.92999999999995;
    prefactor_units[125]  = 1.0000000000000002e-06;
    activation_units[125] = 0.50321666580471969;
    phase_units[125]      = pow(10,-12.000000);
    is_PD[125] = 0;
    nTB[125] = 0;

    // (126):  C3H5XA + O2 => CH2CHO + CH2O
    kiv[126] = {26,5,25,10};
    nuv[126] = {-1,-1,1,1};
    // (126):  C3H5XA + O2 => CH2CHO + CH2O
    fwd_A[126]     = 7140000000000000;
    fwd_beta[126]  = -1.21;
    fwd_Ea[126]    = 21049.950000000001;
    prefactor_units[126]  = 1.0000000000000002e-06;
    activation_units[126] = 0.50321666580471969;
    phase_units[126]      = pow(10,-12.000000);
    is_PD[126] = 0;
    nTB[126] = 0;

    // (127):  C3H5XA + O2 => C2H2 + CH2O + OH
    kiv[127] = {26,5,22,10,3};
    nuv[127] = {-1,-1,1,1,1};
    // (127):  C3H5XA + O2 => C2H2 + CH2O + OH
    fwd_A[127]     = 9.7200000000000003e+29;
    fwd_beta[127]  = -5.71;
    fwd_Ea[127]    = 21450.049999999999;
    prefactor_units[127]  = 1.0000000000000002e-06;
    activation_units[127] = 0.50321666580471969;
    phase_units[127]      = pow(10,-12.000000);
    is_PD[127] = 0;
    nTB[127] = 0;

    // (128):  C3H5XA + H => C3H6
    kiv[128] = {26,2,27};
    nuv[128] = {-1,-1,1};
    // (128):  C3H5XA + H => C3H6
    fwd_A[128]     = 4.8869999999999999e+56;
    fwd_beta[128]  = -12.25;
    fwd_Ea[128]    = 28080.07;
    prefactor_units[128]  = 1.0000000000000002e-06;
    activation_units[128] = 0.50321666580471969;
    phase_units[128]      = pow(10,-12.000000);
    is_PD[128] = 0;
    nTB[128] = 0;

    // (129):  C3H5XA => C2H2 + CH3
    kiv[129] = {26,22,11};
    nuv[129] = {-1,1,1};
    // (129):  C3H5XA => C2H2 + CH3
    fwd_A[129]     = 2.397e+48;
    fwd_beta[129]  = -9.9000000000000004;
    fwd_Ea[129]    = 82080.070000000007;
    prefactor_units[129]  = 1;
    activation_units[129] = 0.50321666580471969;
    phase_units[129]      = pow(10,-6.000000);
    is_PD[129] = 0;
    nTB[129] = 0;

    // (130):  C3H6 => C2H3 + CH3
    kiv[130] = {27,24,11};
    nuv[130] = {-1,1,1};
    // (130):  C3H6 => C2H3 + CH3
    fwd_A[130]     = 2.7300000000000001e+62;
    fwd_beta[130]  = -13.279999999999999;
    fwd_Ea[130]    = 123200.05;
    prefactor_units[130]  = 1;
    activation_units[130] = 0.50321666580471969;
    phase_units[130]      = pow(10,-6.000000);
    is_PD[130] = 0;
    nTB[130] = 0;

    // (131):  C2H3 + CH3 => C3H6
    kiv[131] = {24,11,27};
    nuv[131] = {-1,-1,1};
    // (131):  C2H3 + CH3 => C3H6
    fwd_A[131]     = 4.7119999999999996e+59;
    fwd_beta[131]  = -13.19;
    fwd_Ea[131]    = 29539.91;
    prefactor_units[131]  = 1.0000000000000002e-06;
    activation_units[131] = 0.50321666580471969;
    phase_units[131]      = pow(10,-12.000000);
    is_PD[131] = 0;
    nTB[131] = 0;

    // (132):  C3H6 + O => C2H5 + HCO
    kiv[132] = {27,1,14,19};
    nuv[132] = {-1,-1,1,1};
    // (132):  C3H6 + O => C2H5 + HCO
    fwd_A[132]     = 15800000;
    fwd_beta[132]  = 1.76;
    fwd_Ea[132]    = -1216.0599999999999;
    prefactor_units[132]  = 1.0000000000000002e-06;
    activation_units[132] = 0.50321666580471969;
    phase_units[132]      = pow(10,-12.000000);
    is_PD[132] = 0;
    nTB[132] = 0;

    // (133):  C3H6 + H => C2H4 + CH3
    kiv[133] = {27,2,16,11};
    nuv[133] = {-1,-1,1,1};
    // (133):  C3H6 + H => C2H4 + CH3
    fwd_A[133]     = 4.8299999999999998e+33;
    fwd_beta[133]  = -5.8099999999999996;
    fwd_Ea[133]    = 18500;
    prefactor_units[133]  = 1.0000000000000002e-06;
    activation_units[133] = 0.50321666580471969;
    phase_units[133]      = pow(10,-12.000000);
    is_PD[133] = 0;
    nTB[133] = 0;

    // (134):  C2H4 + CH3 => C3H6 + H
    kiv[134] = {16,11,27,2};
    nuv[134] = {-1,-1,1,1};
    // (134):  C2H4 + CH3 => C3H6 + H
    fwd_A[134]     = 2.3130000000000001e+33;
    fwd_beta[134]  = -5.9000000000000004;
    fwd_Ea[134]    = 31619.98;
    prefactor_units[134]  = 1.0000000000000002e-06;
    activation_units[134] = 0.50321666580471969;
    phase_units[134]      = pow(10,-12.000000);
    is_PD[134] = 0;
    nTB[134] = 0;

    // (135):  C3H6 + H => C3H5XA + H2
    kiv[135] = {27,2,26,4};
    nuv[135] = {-1,-1,1,1};
    // (135):  C3H6 + H => C3H5XA + H2
    fwd_A[135]     = 173000;
    fwd_beta[135]  = 2.5;
    fwd_Ea[135]    = 2492.1100000000001;
    prefactor_units[135]  = 1.0000000000000002e-06;
    activation_units[135] = 0.50321666580471969;
    phase_units[135]      = pow(10,-12.000000);
    is_PD[135] = 0;
    nTB[135] = 0;

    // (136):  C3H5XA + H2 => C3H6 + H
    kiv[136] = {26,4,27,2};
    nuv[136] = {-1,-1,1,1};
    // (136):  C3H5XA + H2 => C3H6 + H
    fwd_A[136]     = 79330;
    fwd_beta[136]  = 2.5099999999999998;
    fwd_Ea[136]    = 19520.080000000002;
    prefactor_units[136]  = 1.0000000000000002e-06;
    activation_units[136] = 0.50321666580471969;
    phase_units[136]      = pow(10,-12.000000);
    is_PD[136] = 0;
    nTB[136] = 0;

    // (137):  C3H6 + OH => C3H5XA + H2O
    kiv[137] = {27,3,26,6};
    nuv[137] = {-1,-1,1,1};
    // (137):  C3H6 + OH => C3H5XA + H2O
    fwd_A[137]     = 3120000;
    fwd_beta[137]  = 2;
    fwd_Ea[137]    = -298.04000000000002;
    prefactor_units[137]  = 1.0000000000000002e-06;
    activation_units[137] = 0.50321666580471969;
    phase_units[137]      = pow(10,-12.000000);
    is_PD[137] = 0;
    nTB[137] = 0;

    // (138):  C3H6 + O => C3H5XA + OH
    kiv[138] = {27,1,26,3};
    nuv[138] = {-1,-1,1,1};
    // (138):  C3H6 + O => C3H5XA + OH
    fwd_A[138]     = 524000000000;
    fwd_beta[138]  = 0.69999999999999996;
    fwd_Ea[138]    = 5884.0799999999999;
    prefactor_units[138]  = 1.0000000000000002e-06;
    activation_units[138] = 0.50321666580471969;
    phase_units[138]      = pow(10,-12.000000);
    is_PD[138] = 0;
    nTB[138] = 0;

    // (139):  C3H6 + CH3 => C3H5XA + CH4
    kiv[139] = {27,11,26,15};
    nuv[139] = {-1,-1,1,1};
    // (139):  C3H6 + CH3 => C3H5XA + CH4
    fwd_A[139]     = 2.21;
    fwd_beta[139]  = 3.5;
    fwd_Ea[139]    = 5674.9499999999998;
    prefactor_units[139]  = 1.0000000000000002e-06;
    activation_units[139] = 0.50321666580471969;
    phase_units[139]      = pow(10,-12.000000);
    is_PD[139] = 0;
    nTB[139] = 0;

    // (140):  IXC3H7 + O2 => C3H6 + HO2
    kiv[140] = {29,5,27,8};
    nuv[140] = {-1,-1,1,1};
    // (140):  IXC3H7 + O2 => C3H6 + HO2
    fwd_A[140]     = 450000000000;
    fwd_beta[140]  = 0;
    fwd_Ea[140]    = 5020.0799999999999;
    prefactor_units[140]  = 1.0000000000000002e-06;
    activation_units[140] = 0.50321666580471969;
    phase_units[140]      = pow(10,-12.000000);
    is_PD[140] = 0;
    nTB[140] = 0;

    // (141):  IXC3H7 => H + C3H6
    kiv[141] = {29,2,27};
    nuv[141] = {-1,1,1};
    // (141):  IXC3H7 => H + C3H6
    fwd_A[141]     = 8.569e+18;
    fwd_beta[141]  = -1.5700000000000001;
    fwd_Ea[141]    = 40340.110000000001;
    prefactor_units[141]  = 1;
    activation_units[141] = 0.50321666580471969;
    phase_units[141]      = pow(10,-6.000000);
    is_PD[141] = 0;
    nTB[141] = 0;

    // (142):  H + C3H6 => IXC3H7
    kiv[142] = {2,27,29};
    nuv[142] = {-1,-1,1};
    // (142):  H + C3H6 => IXC3H7
    fwd_A[142]     = 13000000000000;
    fwd_beta[142]  = 0;
    fwd_Ea[142]    = 1559.99;
    prefactor_units[142]  = 1.0000000000000002e-06;
    activation_units[142] = 0.50321666580471969;
    phase_units[142]      = pow(10,-12.000000);
    is_PD[142] = 0;
    nTB[142] = 0;

    // (143):  NXC3H7 => CH3 + C2H4
    kiv[143] = {30,11,16};
    nuv[143] = {-1,1,1};
    // (143):  NXC3H7 => CH3 + C2H4
    fwd_A[143]     = 228400000000000;
    fwd_beta[143]  = -0.55000000000000004;
    fwd_Ea[143]    = 28400.099999999999;
    prefactor_units[143]  = 1;
    activation_units[143] = 0.50321666580471969;
    phase_units[143]      = pow(10,-6.000000);
    is_PD[143] = 0;
    nTB[143] = 0;

    // (144):  CH3 + C2H4 => NXC3H7
    kiv[144] = {11,16,30};
    nuv[144] = {-1,-1,1};
    // (144):  CH3 + C2H4 => NXC3H7
    fwd_A[144]     = 410000000000;
    fwd_beta[144]  = 0;
    fwd_Ea[144]    = 7204.1099999999997;
    prefactor_units[144]  = 1.0000000000000002e-06;
    activation_units[144] = 0.50321666580471969;
    phase_units[144]      = pow(10,-12.000000);
    is_PD[144] = 0;
    nTB[144] = 0;

    // (145):  NXC3H7 + HO2 => C3H8 + O2
    kiv[145] = {30,8,31,5};
    nuv[145] = {-1,-1,1,1};
    // (145):  NXC3H7 + HO2 => C3H8 + O2
    fwd_A[145]     = 2080000000000;
    fwd_beta[145]  = 0;
    fwd_Ea[145]    = 0;
    prefactor_units[145]  = 1.0000000000000002e-06;
    activation_units[145] = 0.50321666580471969;
    phase_units[145]      = pow(10,-12.000000);
    is_PD[145] = 0;
    nTB[145] = 0;

    // (146):  NXC3H7 + O2 => C3H6 + HO2
    kiv[146] = {30,5,27,8};
    nuv[146] = {-1,-1,1,1};
    // (146):  NXC3H7 + O2 => C3H6 + HO2
    fwd_A[146]     = 300000000000;
    fwd_beta[146]  = 0;
    fwd_Ea[146]    = 3000;
    prefactor_units[146]  = 1.0000000000000002e-06;
    activation_units[146] = 0.50321666580471969;
    phase_units[146]      = pow(10,-12.000000);
    is_PD[146] = 0;
    nTB[146] = 0;

    // (147):  NXC3H7 => H + C3H6
    kiv[147] = {30,2,27};
    nuv[147] = {-1,1,1};
    // (147):  NXC3H7 => H + C3H6
    fwd_A[147]     = 2667000000000000;
    fwd_beta[147]  = -0.64000000000000001;
    fwd_Ea[147]    = 36820.029999999999;
    prefactor_units[147]  = 1;
    activation_units[147] = 0.50321666580471969;
    phase_units[147]      = pow(10,-6.000000);
    is_PD[147] = 0;
    nTB[147] = 0;

    // (148):  H + C3H6 => NXC3H7
    kiv[148] = {2,27,30};
    nuv[148] = {-1,-1,1};
    // (148):  H + C3H6 => NXC3H7
    fwd_A[148]     = 10000000000000;
    fwd_beta[148]  = 0;
    fwd_Ea[148]    = 2500;
    prefactor_units[148]  = 1.0000000000000002e-06;
    activation_units[148] = 0.50321666580471969;
    phase_units[148]      = pow(10,-12.000000);
    is_PD[148] = 0;
    nTB[148] = 0;

    // (149):  C3H8 + OH => NXC3H7 + H2O
    kiv[149] = {31,3,30,6};
    nuv[149] = {-1,-1,1,1};
    // (149):  C3H8 + OH => NXC3H7 + H2O
    fwd_A[149]     = 10540000000;
    fwd_beta[149]  = 0.96999999999999997;
    fwd_Ea[149]    = 1586.04;
    prefactor_units[149]  = 1.0000000000000002e-06;
    activation_units[149] = 0.50321666580471969;
    phase_units[149]      = pow(10,-12.000000);
    is_PD[149] = 0;
    nTB[149] = 0;

    // (150):  C3H8 + HO2 => NXC3H7 + H2O2
    kiv[150] = {31,8,30,7};
    nuv[150] = {-1,-1,1,1};
    // (150):  C3H8 + HO2 => NXC3H7 + H2O2
    fwd_A[150]     = 16800000000000;
    fwd_beta[150]  = 0;
    fwd_Ea[150]    = 20429.970000000001;
    prefactor_units[150]  = 1.0000000000000002e-06;
    activation_units[150] = 0.50321666580471969;
    phase_units[150]      = pow(10,-12.000000);
    is_PD[150] = 0;
    nTB[150] = 0;

    // (151):  H + C3H8 <=> H2 + NXC3H7
    kiv[151] = {2,31,4,30};
    nuv[151] = {-1,-1,1,1};
    // (151):  H + C3H8 <=> H2 + NXC3H7
    fwd_A[151]     = 3972000;
    fwd_beta[151]  = 2.75;
    fwd_Ea[151]    = 6756.6899999999996;
    prefactor_units[151]  = 1.0000000000000002e-06;
    activation_units[151] = 0.50321666580471969;
    phase_units[151]      = pow(10,-12.000000);
    is_PD[151] = 0;
    nTB[151] = 0;

    // (152):  C3H8 + OH => IXC3H7 + H2O
    kiv[152] = {31,3,29,6};
    nuv[152] = {-1,-1,1,1};
    // (152):  C3H8 + OH => IXC3H7 + H2O
    fwd_A[152]     = 46700000;
    fwd_beta[152]  = 1.6100000000000001;
    fwd_Ea[152]    = -34.890000000000001;
    prefactor_units[152]  = 1.0000000000000002e-06;
    activation_units[152] = 0.50321666580471969;
    phase_units[152]      = pow(10,-12.000000);
    is_PD[152] = 0;
    nTB[152] = 0;

    // (153):  CH3 + C3H8 => CH4 + IXC3H7
    kiv[153] = {11,31,15,29};
    nuv[153] = {-1,-1,1,1};
    // (153):  CH3 + C3H8 => CH4 + IXC3H7
    fwd_A[153]     = 398000000000;
    fwd_beta[153]  = 0;
    fwd_Ea[153]    = 9500;
    prefactor_units[153]  = 1.0000000000000002e-06;
    activation_units[153] = 0.50321666580471969;
    phase_units[153]      = pow(10,-12.000000);
    is_PD[153] = 0;
    nTB[153] = 0;

    // (154):  CH3 + C3H8 => CH4 + NXC3H7
    kiv[154] = {11,31,15,30};
    nuv[154] = {-1,-1,1,1};
    // (154):  CH3 + C3H8 => CH4 + NXC3H7
    fwd_A[154]     = 1290000000000;
    fwd_beta[154]  = 0;
    fwd_Ea[154]    = 11599.9;
    prefactor_units[154]  = 1.0000000000000002e-06;
    activation_units[154] = 0.50321666580471969;
    phase_units[154]      = pow(10,-12.000000);
    is_PD[154] = 0;
    nTB[154] = 0;

    // (155):  C3H8 + O => IXC3H7 + OH
    kiv[155] = {31,1,29,3};
    nuv[155] = {-1,-1,1,1};
    // (155):  C3H8 + O => IXC3H7 + OH
    fwd_A[155]     = 28100000000000;
    fwd_beta[155]  = 0;
    fwd_Ea[155]    = 5200.0500000000002;
    prefactor_units[155]  = 1.0000000000000002e-06;
    activation_units[155] = 0.50321666580471969;
    phase_units[155]      = pow(10,-12.000000);
    is_PD[155] = 0;
    nTB[155] = 0;

    // (156):  C3H8 + HO2 => IXC3H7 + H2O2
    kiv[156] = {31,8,29,7};
    nuv[156] = {-1,-1,1,1};
    // (156):  C3H8 + HO2 => IXC3H7 + H2O2
    fwd_A[156]     = 5600000000000;
    fwd_beta[156]  = 0;
    fwd_Ea[156]    = 17700.049999999999;
    prefactor_units[156]  = 1.0000000000000002e-06;
    activation_units[156] = 0.50321666580471969;
    phase_units[156]      = pow(10,-12.000000);
    is_PD[156] = 0;
    nTB[156] = 0;

    // (157):  C3H8 + O => NXC3H7 + OH
    kiv[157] = {31,1,30,3};
    nuv[157] = {-1,-1,1,1};
    // (157):  C3H8 + O => NXC3H7 + OH
    fwd_A[157]     = 113000000000000;
    fwd_beta[157]  = 0;
    fwd_Ea[157]    = 7849.8999999999996;
    prefactor_units[157]  = 1.0000000000000002e-06;
    activation_units[157] = 0.50321666580471969;
    phase_units[157]      = pow(10,-12.000000);
    is_PD[157] = 0;
    nTB[157] = 0;

    // (158):  C3H8 + O2 => IXC3H7 + HO2
    kiv[158] = {31,5,29,8};
    nuv[158] = {-1,-1,1,1};
    // (158):  C3H8 + O2 => IXC3H7 + HO2
    fwd_A[158]     = 40000000000000;
    fwd_beta[158]  = 0;
    fwd_Ea[158]    = 47500;
    prefactor_units[158]  = 1.0000000000000002e-06;
    activation_units[158] = 0.50321666580471969;
    phase_units[158]      = pow(10,-12.000000);
    is_PD[158] = 0;
    nTB[158] = 0;

    // (159):  IXC3H7 + HO2 => C3H8 + O2
    kiv[159] = {29,8,31,5};
    nuv[159] = {-1,-1,1,1};
    // (159):  IXC3H7 + HO2 => C3H8 + O2
    fwd_A[159]     = 2080000000000;
    fwd_beta[159]  = 0;
    fwd_Ea[159]    = 0;
    prefactor_units[159]  = 1.0000000000000002e-06;
    activation_units[159] = 0.50321666580471969;
    phase_units[159]      = pow(10,-12.000000);
    is_PD[159] = 0;
    nTB[159] = 0;

    // (160):  H + C3H8 => H2 + IXC3H7
    kiv[160] = {2,31,4,29};
    nuv[160] = {-1,-1,1,1};
    // (160):  H + C3H8 => H2 + IXC3H7
    fwd_A[160]     = 1300000;
    fwd_beta[160]  = 2.3999999999999999;
    fwd_Ea[160]    = 4471.0799999999999;
    prefactor_units[160]  = 1.0000000000000002e-06;
    activation_units[160] = 0.50321666580471969;
    phase_units[160]      = pow(10,-12.000000);
    is_PD[160] = 0;
    nTB[160] = 0;

    // (161):  H2 + IXC3H7 => H + C3H8
    kiv[161] = {4,29,2,31};
    nuv[161] = {-1,-1,1,1};
    // (161):  H2 + IXC3H7 => H + C3H8
    fwd_A[161]     = 470900;
    fwd_beta[161]  = 2.1499999999999999;
    fwd_Ea[161]    = 12179.969999999999;
    prefactor_units[161]  = 1.0000000000000002e-06;
    activation_units[161] = 0.50321666580471969;
    phase_units[161]      = pow(10,-12.000000);
    is_PD[161] = 0;
    nTB[161] = 0;

    // (162):  C3H5O => C2H3 + CH2O
    kiv[162] = {28,24,10};
    nuv[162] = {-1,1,1};
    // (162):  C3H5O => C2H3 + CH2O
    fwd_A[162]     = 2028000000000;
    fwd_beta[162]  = 0.089999999999999997;
    fwd_Ea[162]    = 23559.990000000002;
    prefactor_units[162]  = 1;
    activation_units[162] = 0.50321666580471969;
    phase_units[162]      = pow(10,-6.000000);
    is_PD[162] = 0;
    nTB[162] = 0;

    // (163):  IXC3H7O2 => IXC3H7 + O2
    kiv[163] = {32,29,5};
    nuv[163] = {-1,1,1};
    // (163):  IXC3H7O2 => IXC3H7 + O2
    fwd_A[163]     = 2.803e+17;
    fwd_beta[163]  = -0.62;
    fwd_Ea[163]    = 36039.910000000003;
    prefactor_units[163]  = 1;
    activation_units[163] = 0.50321666580471969;
    phase_units[163]      = pow(10,-6.000000);
    is_PD[163] = 0;
    nTB[163] = 0;

    // (164):  IXC3H7 + O2 => IXC3H7O2
    kiv[164] = {29,5,32};
    nuv[164] = {-1,-1,1};
    // (164):  IXC3H7 + O2 => IXC3H7O2
    fwd_A[164]     = 7540000000000;
    fwd_beta[164]  = 0;
    fwd_Ea[164]    = 0;
    prefactor_units[164]  = 1.0000000000000002e-06;
    activation_units[164] = 0.50321666580471969;
    phase_units[164]      = pow(10,-12.000000);
    is_PD[164] = 0;
    nTB[164] = 0;

    // (165):  NXC3H7O2 => NXC3H7 + O2
    kiv[165] = {33,30,5};
    nuv[165] = {-1,1,1};
    // (165):  NXC3H7O2 => NXC3H7 + O2
    fwd_A[165]     = 3.364e+19;
    fwd_beta[165]  = -1.3200000000000001;
    fwd_Ea[165]    = 35760.040000000001;
    prefactor_units[165]  = 1;
    activation_units[165] = 0.50321666580471969;
    phase_units[165]      = pow(10,-6.000000);
    is_PD[165] = 0;
    nTB[165] = 0;

    // (166):  NXC3H7 + O2 => NXC3H7O2
    kiv[166] = {30,5,33};
    nuv[166] = {-1,-1,1};
    // (166):  NXC3H7 + O2 => NXC3H7O2
    fwd_A[166]     = 4520000000000;
    fwd_beta[166]  = 0;
    fwd_Ea[166]    = 0;
    prefactor_units[166]  = 1.0000000000000002e-06;
    activation_units[166] = 0.50321666580471969;
    phase_units[166]      = pow(10,-12.000000);
    is_PD[166] = 0;
    nTB[166] = 0;

}


/* Finalizes parameter database */
void CKFINALIZE()
{
  for (int i=0; i<167; ++i) {
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
    awt[1] = 15.999400; /*O */
    awt[2] = 1.007970; /*H */
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
    for (id = 0; id < kd * 34; ++ id) {
         ncf[id] = 0; 
    }

    /*N2 */
    ncf[ 0 * kd + 0 ] = 2; /*N */

    /*O */
    ncf[ 1 * kd + 1 ] = 1; /*O */

    /*H */
    ncf[ 2 * kd + 2 ] = 1; /*H */

    /*OH */
    ncf[ 3 * kd + 2 ] = 1; /*H */
    ncf[ 3 * kd + 1 ] = 1; /*O */

    /*H2 */
    ncf[ 4 * kd + 2 ] = 2; /*H */

    /*O2 */
    ncf[ 5 * kd + 1 ] = 2; /*O */

    /*H2O */
    ncf[ 6 * kd + 2 ] = 2; /*H */
    ncf[ 6 * kd + 1 ] = 1; /*O */

    /*H2O2 */
    ncf[ 7 * kd + 2 ] = 2; /*H */
    ncf[ 7 * kd + 1 ] = 2; /*O */

    /*HO2 */
    ncf[ 8 * kd + 2 ] = 1; /*H */
    ncf[ 8 * kd + 1 ] = 2; /*O */

    /*CH2GSG */
    ncf[ 9 * kd + 3 ] = 1; /*C */
    ncf[ 9 * kd + 2 ] = 2; /*H */

    /*CH2O */
    ncf[ 10 * kd + 3 ] = 1; /*C */
    ncf[ 10 * kd + 2 ] = 2; /*H */
    ncf[ 10 * kd + 1 ] = 1; /*O */

    /*CH3 */
    ncf[ 11 * kd + 3 ] = 1; /*C */
    ncf[ 11 * kd + 2 ] = 3; /*H */

    /*CO */
    ncf[ 12 * kd + 3 ] = 1; /*C */
    ncf[ 12 * kd + 1 ] = 1; /*O */

    /*CH3O */
    ncf[ 13 * kd + 3 ] = 1; /*C */
    ncf[ 13 * kd + 2 ] = 3; /*H */
    ncf[ 13 * kd + 1 ] = 1; /*O */

    /*C2H5 */
    ncf[ 14 * kd + 3 ] = 2; /*C */
    ncf[ 14 * kd + 2 ] = 5; /*H */

    /*CH4 */
    ncf[ 15 * kd + 3 ] = 1; /*C */
    ncf[ 15 * kd + 2 ] = 4; /*H */

    /*C2H4 */
    ncf[ 16 * kd + 3 ] = 2; /*C */
    ncf[ 16 * kd + 2 ] = 4; /*H */

    /*C2H6 */
    ncf[ 17 * kd + 3 ] = 2; /*C */
    ncf[ 17 * kd + 2 ] = 6; /*H */

    /*CO2 */
    ncf[ 18 * kd + 3 ] = 1; /*C */
    ncf[ 18 * kd + 1 ] = 2; /*O */

    /*HCO */
    ncf[ 19 * kd + 2 ] = 1; /*H */
    ncf[ 19 * kd + 3 ] = 1; /*C */
    ncf[ 19 * kd + 1 ] = 1; /*O */

    /*CH3O2 */
    ncf[ 20 * kd + 3 ] = 1; /*C */
    ncf[ 20 * kd + 2 ] = 3; /*H */
    ncf[ 20 * kd + 1 ] = 2; /*O */

    /*CH3O2H */
    ncf[ 21 * kd + 3 ] = 1; /*C */
    ncf[ 21 * kd + 2 ] = 4; /*H */
    ncf[ 21 * kd + 1 ] = 2; /*O */

    /*C2H2 */
    ncf[ 22 * kd + 3 ] = 2; /*C */
    ncf[ 22 * kd + 2 ] = 2; /*H */

    /*HCCO */
    ncf[ 23 * kd + 2 ] = 1; /*H */
    ncf[ 23 * kd + 3 ] = 2; /*C */
    ncf[ 23 * kd + 1 ] = 1; /*O */

    /*C2H3 */
    ncf[ 24 * kd + 3 ] = 2; /*C */
    ncf[ 24 * kd + 2 ] = 3; /*H */

    /*CH2CHO */
    ncf[ 25 * kd + 1 ] = 1; /*O */
    ncf[ 25 * kd + 2 ] = 3; /*H */
    ncf[ 25 * kd + 3 ] = 2; /*C */

    /*C3H5XA */
    ncf[ 26 * kd + 3 ] = 3; /*C */
    ncf[ 26 * kd + 2 ] = 5; /*H */

    /*C3H6 */
    ncf[ 27 * kd + 3 ] = 3; /*C */
    ncf[ 27 * kd + 2 ] = 6; /*H */

    /*C3H5O */
    ncf[ 28 * kd + 3 ] = 3; /*C */
    ncf[ 28 * kd + 2 ] = 5; /*H */
    ncf[ 28 * kd + 1 ] = 1; /*O */

    /*IXC3H7 */
    ncf[ 29 * kd + 3 ] = 3; /*C */
    ncf[ 29 * kd + 2 ] = 7; /*H */

    /*NXC3H7 */
    ncf[ 30 * kd + 3 ] = 3; /*C */
    ncf[ 30 * kd + 2 ] = 7; /*H */

    /*C3H8 */
    ncf[ 31 * kd + 3 ] = 3; /*C */
    ncf[ 31 * kd + 2 ] = 8; /*H */

    /*IXC3H7O2 */
    ncf[ 32 * kd + 3 ] = 3; /*C */
    ncf[ 32 * kd + 2 ] = 7; /*H */
    ncf[ 32 * kd + 1 ] = 2; /*O */

    /*NXC3H7O2 */
    ncf[ 33 * kd + 3 ] = 3; /*C */
    ncf[ 33 * kd + 2 ] = 7; /*H */
    ncf[ 33 * kd + 1 ] = 2; /*O */

}


/* Returns the vector of strings of element names */
void CKSYME_STR(amrex::Vector<std::string>& ename)
{
    ename.resize(4);
    ename[0] = "N";
    ename[1] = "O";
    ename[2] = "H";
    ename[3] = "C";
}


/* Returns the vector of strings of species names */
void CKSYMS_STR(amrex::Vector<std::string>& kname)
{
    kname.resize(34);
    kname[0] = "N2";
    kname[1] = "O";
    kname[2] = "H";
    kname[3] = "OH";
    kname[4] = "H2";
    kname[5] = "O2";
    kname[6] = "H2O";
    kname[7] = "H2O2";
    kname[8] = "HO2";
    kname[9] = "CH2GSG";
    kname[10] = "CH2O";
    kname[11] = "CH3";
    kname[12] = "CO";
    kname[13] = "CH3O";
    kname[14] = "C2H5";
    kname[15] = "CH4";
    kname[16] = "C2H4";
    kname[17] = "C2H6";
    kname[18] = "CO2";
    kname[19] = "HCO";
    kname[20] = "CH3O2";
    kname[21] = "CH3O2H";
    kname[22] = "C2H2";
    kname[23] = "HCCO";
    kname[24] = "C2H3";
    kname[25] = "CH2CHO";
    kname[26] = "C3H5XA";
    kname[27] = "C3H6";
    kname[28] = "C3H5O";
    kname[29] = "IXC3H7";
    kname[30] = "NXC3H7";
    kname[31] = "C3H8";
    kname[32] = "IXC3H7O2";
    kname[33] = "NXC3H7O2";
}

/*compute the sparsity pattern of the chemistry Jacobian */
void SPARSITY_INFO( int * nJdata, int * consP, int NCELLS)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(1225);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(34);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[1225];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            for (int k=0; k<34; k++) {
                c_d[k] = 1.0/ 34.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<35; k++) {
        for (int l=0; l<35; l++) {
            if(J_h[ 35 * k + l] != 0.0){
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(1225);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(34);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[1225];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            for (int k=0; k<34; k++) {
                c_d[k] = 1.0/ 34.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<35; k++) {
        for (int l=0; l<35; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J_h[ 35 * k + l] != 0.0){
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(1225);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(34);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[1225];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            for (int k=0; k<34; k++) {
                c_d[k] = 1.0/ 34.000000 ;
            }
            aJacobian_precond(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<35; k++) {
        for (int l=0; l<35; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J_h[ 35 * k + l] != 0.0){
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

    amrex::Gpu::DeviceVector<amrex::Real> J_v(1225);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(34);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[1225];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            for (int k=0; k<34; k++) {
                c_d[k] = 1.0/ 34.000000 ;
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
        offset_row = nc * 35;
        offset_col = nc * 35;
        for (int k=0; k<35; k++) {
            for (int l=0; l<35; l++) {
                if(J_h[35*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(1225);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(34);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[1225];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            for (int k=0; k<34; k++) {
                c_d[k] = 1.0/ 34.000000 ;
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
            offset = nc * 35;
            for (int l=0; l<35; l++) {
                for (int k=0; k<35; k++) {
                    if(J_h[35*k + l] != 0.0) {
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
            offset = nc * 35;
            for (int l=0; l<35; l++) {
                for (int k=0; k<35; k++) {
                    if(J_h[35*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(1225);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(34);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[1225];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            for (int k=0; k<34; k++) {
                c_d[k] = 1.0/ 34.000000 ;
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
            offset = nc * 35;
            for (int l=0; l<35; l++) {
                for (int k=0; k<35; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp-1] = l+1 + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(J_h[35*k + l] != 0.0) {
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
            offset = nc * 35;
            for (int l=0; l<35; l++) {
                for (int k=0; k<35; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp] = l + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(J_h[35*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(1225);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(34);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[1225];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            for (int k=0; k<34; k++) {
                c_d[k] = 1.0/ 34.000000 ;
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
    for (int k=0; k<35; k++) {
        for (int l=0; l<35; l++) {
            if (k == l) {
                rowVals[nJdata_tmp] = l; 
                indx[nJdata_tmp] = 35*k + l;
                nJdata_tmp = nJdata_tmp + 1; 
            } else {
                if(J_h[35*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l; 
                    indx[nJdata_tmp] = 35*k + l;
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(1225);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(34);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[1225];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            for (int k=0; k<34; k++) {
                c_d[k] = 1.0/ 34.000000 ;
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
        for (int l=0; l<35; l++) {
            for (int k=0; k<35; k++) {
                if (k == l) {
                    colVals[nJdata_tmp-1] = l+1; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J_h[35*k + l] != 0.0) {
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
        for (int l=0; l<35; l++) {
            for (int k=0; k<35; k++) {
                if (k == l) {
                    colVals[nJdata_tmp] = l; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J_h[35*k + l] != 0.0) {
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

