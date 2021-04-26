#ifndef MECHANISM_CPP
#define MECHANISM_CPP

#include "mechanism.H"

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
void productionRate_cpu(amrex::Real *  wdot, amrex::Real *  sc, amrex::Real T)
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
    comp_qfqr_cpu(q_f, q_r, sc, sc_qss, tc, invT);

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

void comp_qfqr_cpu(amrex::Real *  qf, amrex::Real *  qr, amrex::Real *  sc, amrex::Real * sc_qss, amrex::Real *  tc, amrex::Real invT)
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
    for (int i = 0; i < 0; ++i) {
        mixture += sc_qss[i];
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
    comp_qfqr_cpu(q_f, q_r, sc, sc_qss, tc, invT);

    for (int i = 0; i < 84; ++i) {
        qdot[i] = q_f[i] - q_r[i];
    }

    return;
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
void aJacobian_cpu(amrex::Real *  J, amrex::Real *  sc, amrex::Real T, int consP)
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
    } else {
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


/* Initializes parameter database */
void CKINIT()
{
#ifndef AMREX_USE_GPU
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
    fwd_A[14]     = 50000;
    fwd_beta[14]  = 2.6699999999999999;
    fwd_Ea[14]    = 6290;
    prefactor_units[14]  = 1.0000000000000002e-06;
    activation_units[14] = 0.50321666580471969;
    phase_units[14]      = pow(10,-12.000000);
    is_PD[14] = 0;
    nTB[14] = 0;
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
    fwd_A[16]     = 80000000000000;
    fwd_beta[16]  = 0;
    fwd_Ea[16]    = 0;
    prefactor_units[16]  = 1.0000000000000002e-06;
    activation_units[16] = 0.50321666580471969;
    phase_units[16]      = pow(10,-12.000000);
    is_PD[16] = 0;
    nTB[16] = 0;
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
    fwd_A[18]     = 84300000000000;
    fwd_beta[18]  = 0;
    fwd_Ea[18]    = 0;
    prefactor_units[18]  = 1.0000000000000002e-06;
    activation_units[18] = 0.50321666580471969;
    phase_units[18]      = pow(10,-12.000000);
    is_PD[18] = 0;
    nTB[18] = 0;
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
    fwd_A[20]     = 30000000000000;
    fwd_beta[20]  = 0;
    fwd_Ea[20]    = 0;
    prefactor_units[20]  = 1.0000000000000002e-06;
    activation_units[20] = 0.50321666580471969;
    phase_units[20]      = pow(10,-12.000000);
    is_PD[20] = 0;
    nTB[20] = 0;
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
    fwd_A[22]     = 39000000000000;
    fwd_beta[22]  = 0;
    fwd_Ea[22]    = 3540;
    prefactor_units[22]  = 1.0000000000000002e-06;
    activation_units[22] = 0.50321666580471969;
    phase_units[22]      = pow(10,-12.000000);
    is_PD[22] = 0;
    nTB[22] = 0;
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
    fwd_A[24]     = 132000000000000;
    fwd_beta[24]  = 0;
    fwd_Ea[24]    = 0;
    prefactor_units[24]  = 1.0000000000000002e-06;
    activation_units[24] = 0.50321666580471969;
    phase_units[24]      = pow(10,-12.000000);
    is_PD[24] = 0;
    nTB[24] = 0;
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
    fwd_A[26]     = 2500000000000;
    fwd_beta[26]  = 0;
    fwd_Ea[26]    = 47800;
    prefactor_units[26]  = 1.0000000000000002e-06;
    activation_units[26] = 0.50321666580471969;
    phase_units[26]      = pow(10,-12.000000);
    is_PD[26] = 0;
    nTB[26] = 0;
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
    fwd_A[28]     = 3e+20;
    fwd_beta[28]  = -1.72;
    fwd_Ea[28]    = 0;
    prefactor_units[28]  = 1.0000000000000002e-12;
    activation_units[28] = 0.50321666580471969;
    phase_units[28]      = pow(10,-18.000000);
    is_PD[28] = 0;
    nTB[28] = 0;
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
    fwd_A[30]     = 3.75e+20;
    fwd_beta[30]  = -1.72;
    fwd_Ea[30]    = 0;
    prefactor_units[30]  = 1.0000000000000002e-12;
    activation_units[30] = 0.50321666580471969;
    phase_units[30]      = pow(10,-18.000000);
    is_PD[30] = 0;
    nTB[30] = 0;
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
    fwd_A[32]     = 83000000000000;
    fwd_beta[32]  = 0;
    fwd_Ea[32]    = 14413;
    prefactor_units[32]  = 1.0000000000000002e-06;
    activation_units[32] = 0.50321666580471969;
    phase_units[32]      = pow(10,-12.000000);
    is_PD[32] = 0;
    nTB[32] = 0;
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
    fwd_A[33]     = 90000000000000000;
    fwd_beta[33]  = -0.59999999999999998;
    fwd_Ea[33]    = 0;
    prefactor_units[33]  = 1.0000000000000002e-12;
    activation_units[33] = 0.50321666580471969;
    phase_units[33]      = pow(10,-18.000000);
    is_PD[33] = 0;
    nTB[33] = 0;
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
    fwd_A[35]     = 5.5e+20;
    fwd_beta[35]  = -2;
    fwd_Ea[35]    = 0;
    prefactor_units[35]  = 1.0000000000000002e-12;
    activation_units[35] = 0.50321666580471969;
    phase_units[35]      = pow(10,-18.000000);
    is_PD[35] = 0;
    nTB[35] = 0;
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
    fwd_A[36]     = 28000000000000;
    fwd_beta[36]  = 0;
    fwd_Ea[36]    = 1068;
    prefactor_units[36]  = 1.0000000000000002e-06;
    activation_units[36] = 0.50321666580471969;
    phase_units[36]      = pow(10,-12.000000);
    is_PD[36] = 0;
    nTB[36] = 0;
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
    fwd_A[38]     = 660000000;
    fwd_beta[38]  = 1.6200000000000001;
    fwd_Ea[38]    = 10840;
    prefactor_units[38]  = 1.0000000000000002e-06;
    activation_units[38] = 0.50321666580471969;
    phase_units[38]      = pow(10,-12.000000);
    is_PD[38] = 0;
    nTB[38] = 0;
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
    fwd_A[39]     = 73400000000000;
    fwd_beta[39]  = 0;
    fwd_Ea[39]    = 0;
    prefactor_units[39]  = 1.0000000000000002e-06;
    activation_units[39] = 0.50321666580471969;
    phase_units[39]      = pow(10,-12.000000);
    is_PD[39] = 0;
    nTB[39] = 0;
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
    fwd_A[40]     = 23000000000;
    fwd_beta[40]  = 1.05;
    fwd_Ea[40]    = 3275;
    prefactor_units[40]  = 1.0000000000000002e-06;
    activation_units[40] = 0.50321666580471969;
    phase_units[40]      = pow(10,-12.000000);
    is_PD[40] = 0;
    nTB[40] = 0;
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
    fwd_A[42]     = 115000000;
    fwd_beta[42]  = 1.8999999999999999;
    fwd_Ea[42]    = 7530;
    prefactor_units[42]  = 1.0000000000000002e-06;
    activation_units[42] = 0.50321666580471969;
    phase_units[42]      = pow(10,-12.000000);
    is_PD[42] = 0;
    nTB[42] = 0;
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
    fwd_A[43]     = 216000000;
    fwd_beta[43]  = 1.51;
    fwd_Ea[43]    = 3430;
    prefactor_units[43]  = 1.0000000000000002e-06;
    activation_units[43] = 0.50321666580471969;
    phase_units[43]      = pow(10,-12.000000);
    is_PD[43] = 0;
    nTB[43] = 0;
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
    fwd_A[45]     = 29000000000000;
    fwd_beta[45]  = 0;
    fwd_Ea[45]    = -500;
    prefactor_units[45]  = 1.0000000000000002e-06;
    activation_units[45] = 0.50321666580471969;
    phase_units[45]      = pow(10,-12.000000);
    is_PD[45] = 0;
    nTB[45] = 0;
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
    fwd_A[47]     = 30000000000000;
    fwd_beta[47]  = 0;
    fwd_Ea[47]    = 0;
    prefactor_units[47]  = 1.0000000000000002e-06;
    activation_units[47] = 0.50321666580471969;
    phase_units[47]      = pow(10,-12.000000);
    is_PD[47] = 0;
    nTB[47] = 0;
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
    fwd_A[49]     = 25010000000000;
    fwd_beta[49]  = 0;
    fwd_Ea[49]    = 0;
    prefactor_units[49]  = 1.0000000000000002e-06;
    activation_units[49] = 0.50321666580471969;
    phase_units[49]      = pow(10,-12.000000);
    is_PD[49] = 0;
    nTB[49] = 0;
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
    fwd_A[51]     = 47600000;
    fwd_beta[51]  = 1.228;
    fwd_Ea[51]    = 70;
    prefactor_units[51]  = 1.0000000000000002e-06;
    activation_units[51] = 0.50321666580471969;
    phase_units[51]      = pow(10,-12.000000);
    is_PD[51] = 0;
    nTB[51] = 0;
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
    fwd_A[53]     = 3430000000;
    fwd_beta[53]  = 1.1799999999999999;
    fwd_Ea[53]    = -447;
    prefactor_units[53]  = 1.0000000000000002e-06;
    activation_units[53] = 0.50321666580471969;
    phase_units[53]      = pow(10,-12.000000);
    is_PD[53] = 0;
    nTB[53] = 0;
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
    fwd_A[55]     = 20000000000000;
    fwd_beta[55]  = 0;
    fwd_Ea[55]    = 0;
    prefactor_units[55]  = 1.0000000000000002e-06;
    activation_units[55] = 0.50321666580471969;
    phase_units[55]      = pow(10,-12.000000);
    is_PD[55] = 0;
    nTB[55] = 0;
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
    fwd_A[57]     = 20000000000000;
    fwd_beta[57]  = 0;
    fwd_Ea[57]    = 0;
    prefactor_units[57]  = 1.0000000000000002e-06;
    activation_units[57] = 0.50321666580471969;
    phase_units[57]      = pow(10,-12.000000);
    is_PD[57] = 0;
    nTB[57] = 0;
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
    fwd_A[59]     = 13200000000000;
    fwd_beta[59]  = 0;
    fwd_Ea[59]    = 1500;
    prefactor_units[59]  = 1.0000000000000002e-06;
    activation_units[59] = 0.50321666580471969;
    phase_units[59]      = pow(10,-12.000000);
    is_PD[59] = 0;
    nTB[59] = 0;
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
    fwd_A[61]     = 40000000000000;
    fwd_beta[61]  = 0;
    fwd_Ea[61]    = 0;
    prefactor_units[61]  = 1.0000000000000002e-06;
    activation_units[61] = 0.50321666580471969;
    phase_units[61]      = pow(10,-12.000000);
    is_PD[61] = 0;
    nTB[61] = 0;
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
    fwd_A[63]     = 15000000000000;
    fwd_beta[63]  = 0;
    fwd_Ea[63]    = 600;
    prefactor_units[63]  = 1.0000000000000002e-06;
    activation_units[63] = 0.50321666580471969;
    phase_units[63]      = pow(10,-12.000000);
    is_PD[63] = 0;
    nTB[63] = 0;
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
    fwd_A[65]     = 28000000000000;
    fwd_beta[65]  = 0;
    fwd_Ea[65]    = 0;
    prefactor_units[65]  = 1.0000000000000002e-06;
    activation_units[65] = 0.50321666580471969;
    phase_units[65]      = pow(10,-12.000000);
    is_PD[65] = 0;
    nTB[65] = 0;
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
    fwd_A[67]     = 70000000000000;
    fwd_beta[67]  = 0;
    fwd_Ea[67]    = 0;
    prefactor_units[67]  = 1.0000000000000002e-06;
    activation_units[67] = 0.50321666580471969;
    phase_units[67]      = pow(10,-12.000000);
    is_PD[67] = 0;
    nTB[67] = 0;
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
    fwd_A[69]     = 12000000000000;
    fwd_beta[69]  = 0;
    fwd_Ea[69]    = -570;
    prefactor_units[69]  = 1.0000000000000002e-06;
    activation_units[69] = 0.50321666580471969;
    phase_units[69]      = pow(10,-12.000000);
    is_PD[69] = 0;
    nTB[69] = 0;
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
    fwd_A[71]     = 9000000000000;
    fwd_beta[71]  = 0;
    fwd_Ea[71]    = 0;
    prefactor_units[71]  = 1.0000000000000002e-06;
    activation_units[71] = 0.50321666580471969;
    phase_units[71]      = pow(10,-12.000000);
    is_PD[71] = 0;
    nTB[71] = 0;
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
    fwd_A[73]     = 14000000000000;
    fwd_beta[73]  = 0;
    fwd_Ea[73]    = 0;
    prefactor_units[73]  = 1.0000000000000002e-06;
    activation_units[73] = 0.50321666580471969;
    phase_units[73]      = pow(10,-12.000000);
    is_PD[73] = 0;
    nTB[73] = 0;
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
    fwd_A[75]     = 36000000000;
    fwd_beta[75]  = 0;
    fwd_Ea[75]    = 8940;
    prefactor_units[75]  = 1.0000000000000002e-06;
    activation_units[75] = 0.50321666580471969;
    phase_units[75]      = pow(10,-12.000000);
    is_PD[75] = 0;
    nTB[75] = 0;
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
    fwd_A[76]     = 4990000000000;
    fwd_beta[76]  = 0.10000000000000001;
    fwd_Ea[76]    = 10600;
    prefactor_units[76]  = 1.0000000000000002e-06;
    activation_units[76] = 0.50321666580471969;
    phase_units[76]      = pow(10,-12.000000);
    is_PD[76] = 0;
    nTB[76] = 0;
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
    fwd_A[78]     = 3320;
    fwd_beta[78]  = 2.8100000000000001;
    fwd_Ea[78]    = 5860;
    prefactor_units[78]  = 1.0000000000000002e-06;
    activation_units[78] = 0.50321666580471969;
    phase_units[78]      = pow(10,-12.000000);
    is_PD[78] = 0;
    nTB[78] = 0;
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
    fwd_A[80]     = 2.244e+18;
    fwd_beta[80]  = -1;
    fwd_Ea[80]    = 17000;
    prefactor_units[80]  = 1.0000000000000002e-06;
    activation_units[80] = 0.50321666580471969;
    phase_units[80]      = pow(10,-12.000000);
    is_PD[80] = 0;
    nTB[80] = 0;
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
    fwd_A[81]     = 7600000000000;
    fwd_beta[81]  = 0;
    fwd_Ea[81]    = 400;
    prefactor_units[81]  = 1.0000000000000002e-06;
    activation_units[81] = 0.50321666580471969;
    phase_units[81]      = pow(10,-12.000000);
    is_PD[81] = 0;
    nTB[81] = 0;
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
    fwd_A[83]     = 840000000000;
    fwd_beta[83]  = 0;
    fwd_Ea[83]    = 3875;
    prefactor_units[83]  = 1.0000000000000002e-06;
    activation_units[83] = 0.50321666580471969;
    phase_units[83]      = pow(10,-12.000000);
    is_PD[83] = 0;
    nTB[83] = 0;
#endif
}


/* Finalizes parameter database */
void CKFINALIZE()
{
#ifndef AMREX_USE_GPU
    for (int i=0; i<84; ++i) {
        free(TB[i]); TB[i] = 0; 
        free(TBid[i]); TBid[i] = 0;
        nTB[i] = 0;
    }
#endif
}


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
