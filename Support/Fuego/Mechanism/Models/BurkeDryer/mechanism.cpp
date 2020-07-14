#include "chemistry_file.H"

#ifndef AMREX_USE_CUDA
namespace thermo
{
    double fwd_A[27], fwd_beta[27], fwd_Ea[27];
    double low_A[27], low_beta[27], low_Ea[27];
    double rev_A[27], rev_beta[27], rev_Ea[27];
    double troe_a[27],troe_Ts[27], troe_Tss[27], troe_Tsss[27];
    double sri_a[27], sri_b[27], sri_c[27], sri_d[27], sri_e[27];
    double activation_units[27], prefactor_units[27], phase_units[27];
    int is_PD[27], troe_len[27], sri_len[27], nTB[27], *TBid[27];
    double *TB[27];
    std::vector<std::vector<double>> kiv(27); 
    std::vector<std::vector<double>> nuv(27); 

    double fwd_A_DEF[27], fwd_beta_DEF[27], fwd_Ea_DEF[27];
    double low_A_DEF[27], low_beta_DEF[27], low_Ea_DEF[27];
    double rev_A_DEF[27], rev_beta_DEF[27], rev_Ea_DEF[27];
    double troe_a_DEF[27],troe_Ts_DEF[27], troe_Tss_DEF[27], troe_Tsss_DEF[27];
    double sri_a_DEF[27], sri_b_DEF[27], sri_c_DEF[27], sri_d_DEF[27], sri_e_DEF[27];
    double activation_units_DEF[27], prefactor_units_DEF[27], phase_units_DEF[27];
    int is_PD_DEF[27], troe_len_DEF[27], sri_len_DEF[27], nTB_DEF[27], *TBid_DEF[27];
    double *TB_DEF[27];
    std::vector<int> rxn_map;
};

using namespace thermo;
#endif

/* Inverse molecular weights */
/* TODO: check necessity on CPU */
static AMREX_GPU_DEVICE_MANAGED double imw[13] = {
    1.0 / 1.007970,  /*H */
    1.0 / 2.015940,  /*H2 */
    1.0 / 15.999400,  /*O */
    1.0 / 17.007370,  /*OH */
    1.0 / 18.015340,  /*H2O */
    1.0 / 31.998800,  /*O2 */
    1.0 / 33.006770,  /*HO2 */
    1.0 / 34.014740,  /*H2O2 */
    1.0 / 28.013400,  /*N2 */
    1.0 / 39.948000,  /*AR */
    1.0 / 4.002600,  /*HE */
    1.0 / 28.010550,  /*CO */
    1.0 / 44.009950};  /*CO2 */

/* Inverse molecular weights */
/* TODO: check necessity because redundant with molecularWeight */
static AMREX_GPU_DEVICE_MANAGED double molecular_weights[13] = {
    1.007970,  /*H */
    2.015940,  /*H2 */
    15.999400,  /*O */
    17.007370,  /*OH */
    18.015340,  /*H2O */
    31.998800,  /*O2 */
    33.006770,  /*HO2 */
    34.014740,  /*H2O2 */
    28.013400,  /*N2 */
    39.948000,  /*AR */
    4.002600,  /*HE */
    28.010550,  /*CO */
    44.009950};  /*CO2 */

AMREX_GPU_HOST_DEVICE
void get_imw(double imw_new[]){
    for(int i = 0; i<13; ++i) imw_new[i] = imw[i];
}

/* TODO: check necessity because redundant with CKWT */
AMREX_GPU_HOST_DEVICE
void get_mw(double mw_new[]){
    for(int i = 0; i<13; ++i) mw_new[i] = molecular_weights[i];
}


#ifndef AMREX_USE_CUDA
/* Initializes parameter database */
void CKINIT()
{

    rxn_map = {6,7,8,9,10,2,11,12,3,13,14,4,5,15,0,16,17,18,19,20,21,1,22,23,24,25,26};

    // (0):  H + O2 <=> O + OH
    kiv[6] = {0,5,2,3};
    nuv[6] = {-1,-1,1,1};
    // (0):  H + O2 <=> O + OH
    fwd_A[6]     = 104000000000000;
    fwd_beta[6]  = 0;
    fwd_Ea[6]    = 15286;
    prefactor_units[6]  = 1.0000000000000002e-06;
    activation_units[6] = 0.50321666580471969;
    phase_units[6]      = pow(10,-12.000000);
    is_PD[6] = 0;
    nTB[6] = 0;

    // (1):  O + H2 <=> H + OH
    kiv[7] = {2,1,0,3};
    nuv[7] = {-1,-1,1,1};
    // (1):  O + H2 <=> H + OH
    fwd_A[7]     = 3818000000000;
    fwd_beta[7]  = 0;
    fwd_Ea[7]    = 7948;
    prefactor_units[7]  = 1.0000000000000002e-06;
    activation_units[7] = 0.50321666580471969;
    phase_units[7]      = pow(10,-12.000000);
    is_PD[7] = 0;
    nTB[7] = 0;

    // (2):  O + H2 <=> H + OH
    kiv[8] = {2,1,0,3};
    nuv[8] = {-1,-1,1,1};
    // (2):  O + H2 <=> H + OH
    fwd_A[8]     = 879200000000000;
    fwd_beta[8]  = 0;
    fwd_Ea[8]    = 19170;
    prefactor_units[8]  = 1.0000000000000002e-06;
    activation_units[8] = 0.50321666580471969;
    phase_units[8]      = pow(10,-12.000000);
    is_PD[8] = 0;
    nTB[8] = 0;

    // (3):  H2 + OH <=> H2O + H
    kiv[9] = {1,3,4,0};
    nuv[9] = {-1,-1,1,1};
    // (3):  H2 + OH <=> H2O + H
    fwd_A[9]     = 216000000;
    fwd_beta[9]  = 1.51;
    fwd_Ea[9]    = 3430;
    prefactor_units[9]  = 1.0000000000000002e-06;
    activation_units[9] = 0.50321666580471969;
    phase_units[9]      = pow(10,-12.000000);
    is_PD[9] = 0;
    nTB[9] = 0;

    // (4):  OH + OH <=> O + H2O
    kiv[10] = {3,3,2,4};
    nuv[10] = {-1,-1,1,1};
    // (4):  OH + OH <=> O + H2O
    fwd_A[10]     = 33400;
    fwd_beta[10]  = 2.4199999999999999;
    fwd_Ea[10]    = -1930;
    prefactor_units[10]  = 1.0000000000000002e-06;
    activation_units[10] = 0.50321666580471969;
    phase_units[10]      = pow(10,-12.000000);
    is_PD[10] = 0;
    nTB[10] = 0;

    // (5):  H2 + M <=> H + H + M
    kiv[2] = {1,0,0};
    nuv[2] = {-1,1,1};
    // (5):  H2 + M <=> H + H + M
    fwd_A[2]     = 4.577e+19;
    fwd_beta[2]  = -1.3999999999999999;
    fwd_Ea[2]    = 104380;
    prefactor_units[2]  = 1.0000000000000002e-06;
    activation_units[2] = 0.50321666580471969;
    phase_units[2]      = pow(10,-6.000000);
    is_PD[2] = 0;
    nTB[2] = 6;
    TB[2] = (double *) malloc(6 * sizeof(double));
    TBid[2] = (int *) malloc(6 * sizeof(int));
    TBid[2][0] = 1; TB[2][0] = 2.5; // H2
    TBid[2][1] = 4; TB[2][1] = 12; // H2O
    TBid[2][2] = 11; TB[2][2] = 1.8999999999999999; // CO
    TBid[2][3] = 12; TB[2][3] = 3.7999999999999998; // CO2
    TBid[2][4] = 9; TB[2][4] = 0; // AR
    TBid[2][5] = 10; TB[2][5] = 0; // HE

    // (6):  H2 + AR <=> H + H + AR
    kiv[11] = {1,9,0,0,9};
    nuv[11] = {-1,-1,1,1,1};
    // (6):  H2 + AR <=> H + H + AR
    fwd_A[11]     = 5.84e+18;
    fwd_beta[11]  = -1.1000000000000001;
    fwd_Ea[11]    = 104380;
    prefactor_units[11]  = 1.0000000000000002e-06;
    activation_units[11] = 0.50321666580471969;
    phase_units[11]      = pow(10,-12.000000);
    is_PD[11] = 0;
    nTB[11] = 0;

    // (7):  H2 + HE <=> H + H + HE
    kiv[12] = {1,10,0,0,10};
    nuv[12] = {-1,-1,1,1,1};
    // (7):  H2 + HE <=> H + H + HE
    fwd_A[12]     = 5.84e+18;
    fwd_beta[12]  = -1.1000000000000001;
    fwd_Ea[12]    = 104380;
    prefactor_units[12]  = 1.0000000000000002e-06;
    activation_units[12] = 0.50321666580471969;
    phase_units[12]      = pow(10,-12.000000);
    is_PD[12] = 0;
    nTB[12] = 0;

    // (8):  O + O + M <=> O2 + M
    kiv[3] = {2,2,5};
    nuv[3] = {-1,-1,1};
    // (8):  O + O + M <=> O2 + M
    fwd_A[3]     = 6165000000000000;
    fwd_beta[3]  = -0.5;
    fwd_Ea[3]    = 0;
    prefactor_units[3]  = 1.0000000000000002e-12;
    activation_units[3] = 0.50321666580471969;
    phase_units[3]      = pow(10,-12.000000);
    is_PD[3] = 0;
    nTB[3] = 6;
    TB[3] = (double *) malloc(6 * sizeof(double));
    TBid[3] = (int *) malloc(6 * sizeof(int));
    TBid[3][0] = 1; TB[3][0] = 2.5; // H2
    TBid[3][1] = 4; TB[3][1] = 12; // H2O
    TBid[3][2] = 9; TB[3][2] = 0; // AR
    TBid[3][3] = 10; TB[3][3] = 0; // HE
    TBid[3][4] = 11; TB[3][4] = 1.8999999999999999; // CO
    TBid[3][5] = 12; TB[3][5] = 3.7999999999999998; // CO2

    // (9):  O + O + AR <=> O2 + AR
    kiv[13] = {2,2,9,5,9};
    nuv[13] = {-1,-1,-1,1,1};
    // (9):  O + O + AR <=> O2 + AR
    fwd_A[13]     = 18860000000000;
    fwd_beta[13]  = 0;
    fwd_Ea[13]    = -1788;
    prefactor_units[13]  = 1.0000000000000002e-12;
    activation_units[13] = 0.50321666580471969;
    phase_units[13]      = pow(10,-18.000000);
    is_PD[13] = 0;
    nTB[13] = 0;

    // (10):  O + O + HE <=> O2 + HE
    kiv[14] = {2,2,10,5,10};
    nuv[14] = {-1,-1,-1,1,1};
    // (10):  O + O + HE <=> O2 + HE
    fwd_A[14]     = 18860000000000;
    fwd_beta[14]  = 0;
    fwd_Ea[14]    = -1788;
    prefactor_units[14]  = 1.0000000000000002e-12;
    activation_units[14] = 0.50321666580471969;
    phase_units[14]      = pow(10,-18.000000);
    is_PD[14] = 0;
    nTB[14] = 0;

    // (11):  O + H + M <=> OH + M
    kiv[4] = {2,0,3};
    nuv[4] = {-1,-1,1};
    // (11):  O + H + M <=> OH + M
    fwd_A[4]     = 4.714e+18;
    fwd_beta[4]  = -1;
    fwd_Ea[4]    = 0;
    prefactor_units[4]  = 1.0000000000000002e-12;
    activation_units[4] = 0.50321666580471969;
    phase_units[4]      = pow(10,-12.000000);
    is_PD[4] = 0;
    nTB[4] = 6;
    TB[4] = (double *) malloc(6 * sizeof(double));
    TBid[4] = (int *) malloc(6 * sizeof(int));
    TBid[4][0] = 1; TB[4][0] = 2.5; // H2
    TBid[4][1] = 4; TB[4][1] = 12; // H2O
    TBid[4][2] = 9; TB[4][2] = 0.75; // AR
    TBid[4][3] = 10; TB[4][3] = 0.75; // HE
    TBid[4][4] = 11; TB[4][4] = 1.8999999999999999; // CO
    TBid[4][5] = 12; TB[4][5] = 3.7999999999999998; // CO2

    // (12):  H2O + M <=> H + OH + M
    kiv[5] = {4,0,3};
    nuv[5] = {-1,1,1};
    // (12):  H2O + M <=> H + OH + M
    fwd_A[5]     = 6.0640000000000002e+27;
    fwd_beta[5]  = -3.3220000000000001;
    fwd_Ea[5]    = 120790;
    prefactor_units[5]  = 1.0000000000000002e-06;
    activation_units[5] = 0.50321666580471969;
    phase_units[5]      = pow(10,-6.000000);
    is_PD[5] = 0;
    nTB[5] = 7;
    TB[5] = (double *) malloc(7 * sizeof(double));
    TBid[5] = (int *) malloc(7 * sizeof(int));
    TBid[5][0] = 1; TB[5][0] = 3; // H2
    TBid[5][1] = 4; TB[5][1] = 0; // H2O
    TBid[5][2] = 10; TB[5][2] = 1.1000000000000001; // HE
    TBid[5][3] = 8; TB[5][3] = 2; // N2
    TBid[5][4] = 5; TB[5][4] = 1.5; // O2
    TBid[5][5] = 11; TB[5][5] = 1.8999999999999999; // CO
    TBid[5][6] = 12; TB[5][6] = 3.7999999999999998; // CO2

    // (13):  H2O + H2O <=> H + OH + H2O
    kiv[15] = {4,4,0,3,4};
    nuv[15] = {-1,-1,1,1,1};
    // (13):  H2O + H2O <=> H + OH + H2O
    fwd_A[15]     = 1.006e+26;
    fwd_beta[15]  = -2.4399999999999999;
    fwd_Ea[15]    = 120180;
    prefactor_units[15]  = 1.0000000000000002e-06;
    activation_units[15] = 0.50321666580471969;
    phase_units[15]      = pow(10,-12.000000);
    is_PD[15] = 0;
    nTB[15] = 0;

    // (14):  H + O2 (+M) <=> HO2 (+M)
    kiv[0] = {0,5,6};
    nuv[0] = {-1,-1,1};
    // (14):  H + O2 (+M) <=> HO2 (+M)
    fwd_A[0]     = 4650840000000;
    fwd_beta[0]  = 0.44;
    fwd_Ea[0]    = 0;
    low_A[0]     = 6.366e+20;
    low_beta[0]  = -1.72;
    low_Ea[0]    = 524.79999999999995;
    troe_a[0]    = 0.5;
    troe_Tsss[0] = 1.0000000000000001e-30;
    troe_Ts[0]   = 1e+30;
    troe_len[0]  = 3;
    prefactor_units[0]  = 1.0000000000000002e-06;
    activation_units[0] = 0.50321666580471969;
    phase_units[0]      = pow(10,-12.000000);
    is_PD[0] = 1;
    nTB[0] = 7;
    TB[0] = (double *) malloc(7 * sizeof(double));
    TBid[0] = (int *) malloc(7 * sizeof(int));
    TBid[0][0] = 1; TB[0][0] = 2; // H2
    TBid[0][1] = 4; TB[0][1] = 14; // H2O
    TBid[0][2] = 5; TB[0][2] = 0.78000000000000003; // O2
    TBid[0][3] = 11; TB[0][3] = 1.8999999999999999; // CO
    TBid[0][4] = 12; TB[0][4] = 3.7999999999999998; // CO2
    TBid[0][5] = 9; TB[0][5] = 0.67000000000000004; // AR
    TBid[0][6] = 10; TB[0][6] = 0.80000000000000004; // HE

    // (15):  HO2 + H <=> H2 + O2
    kiv[16] = {6,0,1,5};
    nuv[16] = {-1,-1,1,1};
    // (15):  HO2 + H <=> H2 + O2
    fwd_A[16]     = 2750000;
    fwd_beta[16]  = 2.0899999999999999;
    fwd_Ea[16]    = -1451;
    prefactor_units[16]  = 1.0000000000000002e-06;
    activation_units[16] = 0.50321666580471969;
    phase_units[16]      = pow(10,-12.000000);
    is_PD[16] = 0;
    nTB[16] = 0;

    // (16):  HO2 + H <=> OH + OH
    kiv[17] = {6,0,3,3};
    nuv[17] = {-1,-1,1,1};
    // (16):  HO2 + H <=> OH + OH
    fwd_A[17]     = 70790000000000;
    fwd_beta[17]  = 0;
    fwd_Ea[17]    = 295;
    prefactor_units[17]  = 1.0000000000000002e-06;
    activation_units[17] = 0.50321666580471969;
    phase_units[17]      = pow(10,-12.000000);
    is_PD[17] = 0;
    nTB[17] = 0;

    // (17):  HO2 + O <=> O2 + OH
    kiv[18] = {6,2,5,3};
    nuv[18] = {-1,-1,1,1};
    // (17):  HO2 + O <=> O2 + OH
    fwd_A[18]     = 28500000000;
    fwd_beta[18]  = 1;
    fwd_Ea[18]    = -723.92999999999995;
    prefactor_units[18]  = 1.0000000000000002e-06;
    activation_units[18] = 0.50321666580471969;
    phase_units[18]      = pow(10,-12.000000);
    is_PD[18] = 0;
    nTB[18] = 0;

    // (18):  HO2 + OH <=> H2O + O2
    kiv[19] = {6,3,4,5};
    nuv[19] = {-1,-1,1,1};
    // (18):  HO2 + OH <=> H2O + O2
    fwd_A[19]     = 28900000000000;
    fwd_beta[19]  = 0;
    fwd_Ea[19]    = -497;
    prefactor_units[19]  = 1.0000000000000002e-06;
    activation_units[19] = 0.50321666580471969;
    phase_units[19]      = pow(10,-12.000000);
    is_PD[19] = 0;
    nTB[19] = 0;

    // (19):  HO2 + HO2 <=> H2O2 + O2
    kiv[20] = {6,6,7,5};
    nuv[20] = {-1,-1,1,1};
    // (19):  HO2 + HO2 <=> H2O2 + O2
    fwd_A[20]     = 420000000000000;
    fwd_beta[20]  = 0;
    fwd_Ea[20]    = 11982;
    prefactor_units[20]  = 1.0000000000000002e-06;
    activation_units[20] = 0.50321666580471969;
    phase_units[20]      = pow(10,-12.000000);
    is_PD[20] = 0;
    nTB[20] = 0;

    // (20):  HO2 + HO2 <=> H2O2 + O2
    kiv[21] = {6,6,7,5};
    nuv[21] = {-1,-1,1,1};
    // (20):  HO2 + HO2 <=> H2O2 + O2
    fwd_A[21]     = 130000000000;
    fwd_beta[21]  = 0;
    fwd_Ea[21]    = -1629.3;
    prefactor_units[21]  = 1.0000000000000002e-06;
    activation_units[21] = 0.50321666580471969;
    phase_units[21]      = pow(10,-12.000000);
    is_PD[21] = 0;
    nTB[21] = 0;

    // (21):  H2O2 (+M) <=> OH + OH (+M)
    kiv[1] = {7,3,3};
    nuv[1] = {-1,1,1};
    // (21):  H2O2 (+M) <=> OH + OH (+M)
    fwd_A[1]     = 2000000000000;
    fwd_beta[1]  = 0.90000000000000002;
    fwd_Ea[1]    = 48749;
    low_A[1]     = 2.49e+24;
    low_beta[1]  = -2.2999999999999998;
    low_Ea[1]    = 48749;
    troe_a[1]    = 0.42999999999999999;
    troe_Tsss[1] = 1.0000000000000001e-30;
    troe_Ts[1]   = 1e+30;
    troe_len[1]  = 3;
    prefactor_units[1]  = 1;
    activation_units[1] = 0.50321666580471969;
    phase_units[1]      = pow(10,-6.000000);
    is_PD[1] = 1;
    nTB[1] = 8;
    TB[1] = (double *) malloc(8 * sizeof(double));
    TBid[1] = (int *) malloc(8 * sizeof(int));
    TBid[1][0] = 4; TB[1][0] = 7.5; // H2O
    TBid[1][1] = 12; TB[1][1] = 1.6000000000000001; // CO2
    TBid[1][2] = 8; TB[1][2] = 1.5; // N2
    TBid[1][3] = 5; TB[1][3] = 1.2; // O2
    TBid[1][4] = 10; TB[1][4] = 0.65000000000000002; // HE
    TBid[1][5] = 7; TB[1][5] = 7.7000000000000002; // H2O2
    TBid[1][6] = 1; TB[1][6] = 3.7000000000000002; // H2
    TBid[1][7] = 11; TB[1][7] = 2.7999999999999998; // CO

    // (22):  H2O2 + H <=> H2O + OH
    kiv[22] = {7,0,4,3};
    nuv[22] = {-1,-1,1,1};
    // (22):  H2O2 + H <=> H2O + OH
    fwd_A[22]     = 24100000000000;
    fwd_beta[22]  = 0;
    fwd_Ea[22]    = 3970;
    prefactor_units[22]  = 1.0000000000000002e-06;
    activation_units[22] = 0.50321666580471969;
    phase_units[22]      = pow(10,-12.000000);
    is_PD[22] = 0;
    nTB[22] = 0;

    // (23):  H2O2 + H <=> HO2 + H2
    kiv[23] = {7,0,6,1};
    nuv[23] = {-1,-1,1,1};
    // (23):  H2O2 + H <=> HO2 + H2
    fwd_A[23]     = 48200000000000;
    fwd_beta[23]  = 0;
    fwd_Ea[23]    = 7950;
    prefactor_units[23]  = 1.0000000000000002e-06;
    activation_units[23] = 0.50321666580471969;
    phase_units[23]      = pow(10,-12.000000);
    is_PD[23] = 0;
    nTB[23] = 0;

    // (24):  H2O2 + O <=> OH + HO2
    kiv[24] = {7,2,3,6};
    nuv[24] = {-1,-1,1,1};
    // (24):  H2O2 + O <=> OH + HO2
    fwd_A[24]     = 9550000;
    fwd_beta[24]  = 2;
    fwd_Ea[24]    = 3970;
    prefactor_units[24]  = 1.0000000000000002e-06;
    activation_units[24] = 0.50321666580471969;
    phase_units[24]      = pow(10,-12.000000);
    is_PD[24] = 0;
    nTB[24] = 0;

    // (25):  H2O2 + OH <=> HO2 + H2O
    kiv[25] = {7,3,6,4};
    nuv[25] = {-1,-1,1,1};
    // (25):  H2O2 + OH <=> HO2 + H2O
    fwd_A[25]     = 1740000000000;
    fwd_beta[25]  = 0;
    fwd_Ea[25]    = 318;
    prefactor_units[25]  = 1.0000000000000002e-06;
    activation_units[25] = 0.50321666580471969;
    phase_units[25]      = pow(10,-12.000000);
    is_PD[25] = 0;
    nTB[25] = 0;

    // (26):  H2O2 + OH <=> HO2 + H2O
    kiv[26] = {7,3,6,4};
    nuv[26] = {-1,-1,1,1};
    // (26):  H2O2 + OH <=> HO2 + H2O
    fwd_A[26]     = 75900000000000;
    fwd_beta[26]  = 0;
    fwd_Ea[26]    = 7270;
    prefactor_units[26]  = 1.0000000000000002e-06;
    activation_units[26] = 0.50321666580471969;
    phase_units[26]      = pow(10,-12.000000);
    is_PD[26] = 0;
    nTB[26] = 0;

    SetAllDefaults();
}

void GET_REACTION_MAP(int *rmap)
{
    for (int i=0; i<27; ++i) {
        rmap[i] = rxn_map[i] + 1;
    }
}

#include <ReactionData.H>
double* GetParamPtr(int                reaction_id,
                    REACTION_PARAMETER param_id,
                    int                species_id,
                    int                get_default)
{
  double* ret = 0;
  if (reaction_id<0 || reaction_id>=27) {
    printf("Bad reaction id = %d",reaction_id);
    abort();
  };
  int mrid = rxn_map[reaction_id];

  if (param_id == THIRD_BODY) {
    if (species_id<0 || species_id>=13) {
      printf("GetParamPtr: Bad species id = %d",species_id);
      abort();
    }
    if (get_default) {
      for (int i=0; i<nTB_DEF[mrid]; ++i) {
        if (species_id == TBid_DEF[mrid][i]) {
          ret = &(TB_DEF[mrid][i]);
        }
      }
    }
    else {
      for (int i=0; i<nTB[mrid]; ++i) {
        if (species_id == TBid[mrid][i]) {
          ret = &(TB[mrid][i]);
        }
      }
    }
    if (ret == 0) {
      printf("GetParamPtr: No TB for reaction id = %d",reaction_id);
      abort();
    }
  }
  else {
    if (     param_id == FWD_A)     {ret = (get_default ? &(fwd_A_DEF[mrid]) : &(fwd_A[mrid]));}
      else if (param_id == FWD_BETA)  {ret = (get_default ? &(fwd_beta_DEF[mrid]) : &(fwd_beta[mrid]));}
      else if (param_id == FWD_EA)    {ret = (get_default ? &(fwd_Ea_DEF[mrid]) : &(fwd_Ea[mrid]));}
      else if (param_id == LOW_A)     {ret = (get_default ? &(low_A_DEF[mrid]) : &(low_A[mrid]));}
      else if (param_id == LOW_BETA)  {ret = (get_default ? &(low_beta_DEF[mrid]) : &(low_beta[mrid]));}
      else if (param_id == LOW_EA)    {ret = (get_default ? &(low_Ea_DEF[mrid]) : &(low_Ea[mrid]));}
      else if (param_id == REV_A)     {ret = (get_default ? &(rev_A_DEF[mrid]) : &(rev_A[mrid]));}
      else if (param_id == REV_BETA)  {ret = (get_default ? &(rev_beta_DEF[mrid]) : &(rev_beta[mrid]));}
      else if (param_id == REV_EA)    {ret = (get_default ? &(rev_Ea_DEF[mrid]) : &(rev_Ea[mrid]));}
      else if (param_id == TROE_A)    {ret = (get_default ? &(troe_a_DEF[mrid]) : &(troe_a[mrid]));}
      else if (param_id == TROE_TS)   {ret = (get_default ? &(troe_Ts_DEF[mrid]) : &(troe_Ts[mrid]));}
      else if (param_id == TROE_TSS)  {ret = (get_default ? &(troe_Tss_DEF[mrid]) : &(troe_Tss[mrid]));}
      else if (param_id == TROE_TSSS) {ret = (get_default ? &(troe_Tsss_DEF[mrid]) : &(troe_Tsss[mrid]));}
      else if (param_id == SRI_A)     {ret = (get_default ? &(sri_a_DEF[mrid]) : &(sri_a[mrid]));}
      else if (param_id == SRI_B)     {ret = (get_default ? &(sri_b_DEF[mrid]) : &(sri_b[mrid]));}
      else if (param_id == SRI_C)     {ret = (get_default ? &(sri_c_DEF[mrid]) : &(sri_c[mrid]));}
      else if (param_id == SRI_D)     {ret = (get_default ? &(sri_d_DEF[mrid]) : &(sri_d[mrid]));}
      else if (param_id == SRI_E)     {ret = (get_default ? &(sri_e_DEF[mrid]) : &(sri_e[mrid]));}
    else {
      printf("GetParamPtr: Unknown parameter id");
      abort();
    }
  }
  return ret;
}

void ResetAllParametersToDefault()
{
    for (int i=0; i<27; i++) {
        if (nTB[i] != 0) {
            nTB[i] = 0;
            free(TB[i]);
            free(TBid[i]);
        }

        fwd_A[i]    = fwd_A_DEF[i];
        fwd_beta[i] = fwd_beta_DEF[i];
        fwd_Ea[i]   = fwd_Ea_DEF[i];

        low_A[i]    = low_A_DEF[i];
        low_beta[i] = low_beta_DEF[i];
        low_Ea[i]   = low_Ea_DEF[i];

        rev_A[i]    = rev_A_DEF[i];
        rev_beta[i] = rev_beta_DEF[i];
        rev_Ea[i]   = rev_Ea_DEF[i];

        troe_a[i]    = troe_a_DEF[i];
        troe_Ts[i]   = troe_Ts_DEF[i];
        troe_Tss[i]  = troe_Tss_DEF[i];
        troe_Tsss[i] = troe_Tsss_DEF[i];

        sri_a[i] = sri_a_DEF[i];
        sri_b[i] = sri_b_DEF[i];
        sri_c[i] = sri_c_DEF[i];
        sri_d[i] = sri_d_DEF[i];
        sri_e[i] = sri_e_DEF[i];

        is_PD[i]    = is_PD_DEF[i];
        troe_len[i] = troe_len_DEF[i];
        sri_len[i]  = sri_len_DEF[i];

        activation_units[i] = activation_units_DEF[i];
        prefactor_units[i]  = prefactor_units_DEF[i];
        phase_units[i]      = phase_units_DEF[i];

        nTB[i]  = nTB_DEF[i];
        if (nTB[i] != 0) {
           TB[i] = (double *) malloc(sizeof(double) * nTB[i]);
           TBid[i] = (int *) malloc(sizeof(int) * nTB[i]);
           for (int j=0; j<nTB[i]; j++) {
             TB[i][j] = TB_DEF[i][j];
             TBid[i][j] = TBid_DEF[i][j];
           }
        }
    }
}

void SetAllDefaults()
{
    for (int i=0; i<27; i++) {
        if (nTB_DEF[i] != 0) {
            nTB_DEF[i] = 0;
            free(TB_DEF[i]);
            free(TBid_DEF[i]);
        }

        fwd_A_DEF[i]    = fwd_A[i];
        fwd_beta_DEF[i] = fwd_beta[i];
        fwd_Ea_DEF[i]   = fwd_Ea[i];

        low_A_DEF[i]    = low_A[i];
        low_beta_DEF[i] = low_beta[i];
        low_Ea_DEF[i]   = low_Ea[i];

        rev_A_DEF[i]    = rev_A[i];
        rev_beta_DEF[i] = rev_beta[i];
        rev_Ea_DEF[i]   = rev_Ea[i];

        troe_a_DEF[i]    = troe_a[i];
        troe_Ts_DEF[i]   = troe_Ts[i];
        troe_Tss_DEF[i]  = troe_Tss[i];
        troe_Tsss_DEF[i] = troe_Tsss[i];

        sri_a_DEF[i] = sri_a[i];
        sri_b_DEF[i] = sri_b[i];
        sri_c_DEF[i] = sri_c[i];
        sri_d_DEF[i] = sri_d[i];
        sri_e_DEF[i] = sri_e[i];

        is_PD_DEF[i]    = is_PD[i];
        troe_len_DEF[i] = troe_len[i];
        sri_len_DEF[i]  = sri_len[i];

        activation_units_DEF[i] = activation_units[i];
        prefactor_units_DEF[i]  = prefactor_units[i];
        phase_units_DEF[i]      = phase_units[i];

        nTB_DEF[i]  = nTB[i];
        if (nTB_DEF[i] != 0) {
           TB_DEF[i] = (double *) malloc(sizeof(double) * nTB_DEF[i]);
           TBid_DEF[i] = (int *) malloc(sizeof(int) * nTB_DEF[i]);
           for (int j=0; j<nTB_DEF[i]; j++) {
             TB_DEF[i][j] = TB[i][j];
             TBid_DEF[i][j] = TBid[i][j];
           }
        }
    }
}

/* Finalizes parameter database */
void CKFINALIZE()
{
  for (int i=0; i<27; ++i) {
    free(TB[i]); TB[i] = 0; 
    free(TBid[i]); TBid[i] = 0;
    nTB[i] = 0;

    free(TB_DEF[i]); TB_DEF[i] = 0; 
    free(TBid_DEF[i]); TBid_DEF[i] = 0;
    nTB_DEF[i] = 0;
  }
}

#else
/* TODO: Remove on GPU, right now needed by chemistry_module on FORTRAN */
AMREX_GPU_HOST_DEVICE void CKINIT()
{
}

AMREX_GPU_HOST_DEVICE void CKFINALIZE()
{
}

#endif


/*A few mechanism parameters */
void CKINDX(int * mm, int * kk, int * ii, int * nfit)
{
    *mm = 6;
    *kk = 13;
    *ii = 27;
    *nfit = -1; /*Why do you need this anyway ?  */
}



/* ckxnum... for parsing strings  */
void CKXNUM(char * line, int * nexp, int * lout, int * nval, double *  rval, int * kerr, int lenline )
{
    int n,i; /*Loop Counters */
    char cstr[1000];
    char *saveptr;
    char *p; /*String Tokens */
    /* Strip Comments  */
    for (i=0; i<lenline; ++i) {
        if (line[i]=='!') {
            break;
        }
        cstr[i] = line[i];
    }
    cstr[i] = '\0';

    p = strtok_r(cstr," ", &saveptr);
    if (!p) {
        *nval = 0;
        *kerr = 1;
        return;
    }
    for (n=0; n<*nexp; ++n) {
        rval[n] = atof(p);
        p = strtok_r(NULL, " ", &saveptr);
        if (!p) break;
    }
    *nval = n+1;
    if (*nval < *nexp) *kerr = 1;
    return;
}


/* cksnum... for parsing strings  */
void CKSNUM(char * line, int * nexp, int * lout, char * kray, int * nn, int * knum, int * nval, double *  rval, int * kerr, int lenline, int lenkray)
{
    /*Not done yet ... */
}


/* Returns the vector of strings of element names */
void CKSYME_STR(amrex::Vector<std::string>& ename)
{
    ename.resize(6);
    ename[0] = "H";
    ename[1] = "O";
    ename[2] = "N";
    ename[3] = "AR";
    ename[4] = "HE";
    ename[5] = "C";
}


/* Returns the char strings of element names */
void CKSYME(int * kname, int * plenkname )
{
    int i; /*Loop Counter */
    int lenkname = *plenkname;
    /*clear kname */
    for (i=0; i<lenkname*6; i++) {
        kname[i] = ' ';
    }

    /* H  */
    kname[ 0*lenkname + 0 ] = 'H';
    kname[ 0*lenkname + 1 ] = ' ';

    /* O  */
    kname[ 1*lenkname + 0 ] = 'O';
    kname[ 1*lenkname + 1 ] = ' ';

    /* N  */
    kname[ 2*lenkname + 0 ] = 'N';
    kname[ 2*lenkname + 1 ] = ' ';

    /* AR  */
    kname[ 3*lenkname + 0 ] = 'A';
    kname[ 3*lenkname + 1 ] = 'R';
    kname[ 3*lenkname + 2 ] = ' ';

    /* HE  */
    kname[ 4*lenkname + 0 ] = 'H';
    kname[ 4*lenkname + 1 ] = 'E';
    kname[ 4*lenkname + 2 ] = ' ';

    /* C  */
    kname[ 5*lenkname + 0 ] = 'C';
    kname[ 5*lenkname + 1 ] = ' ';

}


/* Returns the vector of strings of species names */
void CKSYMS_STR(amrex::Vector<std::string>& kname)
{
    kname.resize(13);
    kname[0] = "H";
    kname[1] = "H2";
    kname[2] = "O";
    kname[3] = "OH";
    kname[4] = "H2O";
    kname[5] = "O2";
    kname[6] = "HO2";
    kname[7] = "H2O2";
    kname[8] = "N2";
    kname[9] = "AR";
    kname[10] = "HE";
    kname[11] = "CO";
    kname[12] = "CO2";
}


/* Returns the char strings of species names */
void CKSYMS(int * kname, int * plenkname )
{
    int i; /*Loop Counter */
    int lenkname = *plenkname;
    /*clear kname */
    for (i=0; i<lenkname*13; i++) {
        kname[i] = ' ';
    }

    /* H  */
    kname[ 0*lenkname + 0 ] = 'H';
    kname[ 0*lenkname + 1 ] = ' ';

    /* H2  */
    kname[ 1*lenkname + 0 ] = 'H';
    kname[ 1*lenkname + 1 ] = '2';
    kname[ 1*lenkname + 2 ] = ' ';

    /* O  */
    kname[ 2*lenkname + 0 ] = 'O';
    kname[ 2*lenkname + 1 ] = ' ';

    /* OH  */
    kname[ 3*lenkname + 0 ] = 'O';
    kname[ 3*lenkname + 1 ] = 'H';
    kname[ 3*lenkname + 2 ] = ' ';

    /* H2O  */
    kname[ 4*lenkname + 0 ] = 'H';
    kname[ 4*lenkname + 1 ] = '2';
    kname[ 4*lenkname + 2 ] = 'O';
    kname[ 4*lenkname + 3 ] = ' ';

    /* O2  */
    kname[ 5*lenkname + 0 ] = 'O';
    kname[ 5*lenkname + 1 ] = '2';
    kname[ 5*lenkname + 2 ] = ' ';

    /* HO2  */
    kname[ 6*lenkname + 0 ] = 'H';
    kname[ 6*lenkname + 1 ] = 'O';
    kname[ 6*lenkname + 2 ] = '2';
    kname[ 6*lenkname + 3 ] = ' ';

    /* H2O2  */
    kname[ 7*lenkname + 0 ] = 'H';
    kname[ 7*lenkname + 1 ] = '2';
    kname[ 7*lenkname + 2 ] = 'O';
    kname[ 7*lenkname + 3 ] = '2';
    kname[ 7*lenkname + 4 ] = ' ';

    /* N2  */
    kname[ 8*lenkname + 0 ] = 'N';
    kname[ 8*lenkname + 1 ] = '2';
    kname[ 8*lenkname + 2 ] = ' ';

    /* AR  */
    kname[ 9*lenkname + 0 ] = 'A';
    kname[ 9*lenkname + 1 ] = 'R';
    kname[ 9*lenkname + 2 ] = ' ';

    /* HE  */
    kname[ 10*lenkname + 0 ] = 'H';
    kname[ 10*lenkname + 1 ] = 'E';
    kname[ 10*lenkname + 2 ] = ' ';

    /* CO  */
    kname[ 11*lenkname + 0 ] = 'C';
    kname[ 11*lenkname + 1 ] = 'O';
    kname[ 11*lenkname + 2 ] = ' ';

    /* CO2  */
    kname[ 12*lenkname + 0 ] = 'C';
    kname[ 12*lenkname + 1 ] = 'O';
    kname[ 12*lenkname + 2 ] = '2';
    kname[ 12*lenkname + 3 ] = ' ';

}


/* Returns R, Rc, Patm */
void CKRP(double *  ru, double *  ruc, double *  pa)
{
     *ru  = 8.31446261815324e+07; 
     *ruc = 1.98721558317399615845; 
     *pa  = 1.01325e+06; 
}


/*Compute P = rhoRT/W(x) */
void CKPX(double *  rho, double *  T, double *  x, double *  P)
{
    double XW = 0;/* To hold mean molecular wt */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*2.015940; /*H2 */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*31.998800; /*O2 */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    XW += x[9]*39.948000; /*AR */
    XW += x[10]*4.002600; /*HE */
    XW += x[11]*28.010550; /*CO */
    XW += x[12]*44.009950; /*CO2 */
    *P = *rho * 8.31446e+07 * (*T) / XW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(y) */
AMREX_GPU_HOST_DEVICE void CKPY(double *  rho, double *  T, double *  y,  double *  P)
{
    double YOW = 0;/* for computing mean MW */
    YOW += y[0]*imw[0]; /*H */
    YOW += y[1]*imw[1]; /*H2 */
    YOW += y[2]*imw[2]; /*O */
    YOW += y[3]*imw[3]; /*OH */
    YOW += y[4]*imw[4]; /*H2O */
    YOW += y[5]*imw[5]; /*O2 */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*N2 */
    YOW += y[9]*imw[9]; /*AR */
    YOW += y[10]*imw[10]; /*HE */
    YOW += y[11]*imw[11]; /*CO */
    YOW += y[12]*imw[12]; /*CO2 */
    *P = *rho * 8.31446e+07 * (*T) * YOW; /*P = rho*R*T/W */

    return;
}


#ifndef AMREX_USE_CUDA
/*Compute P = rhoRT/W(y) */
void VCKPY(int *  np, double *  rho, double *  T, double *  y,  double *  P)
{
    double YOW[*np];
    for (int i=0; i<(*np); i++) {
        YOW[i] = 0.0;
    }

    for (int n=0; n<13; n++) {
        for (int i=0; i<(*np); i++) {
            YOW[i] += y[n*(*np)+i] * imw[n];
        }
    }

    for (int i=0; i<(*np); i++) {
        P[i] = rho[i] * 8.31446e+07 * T[i] * YOW[i]; /*P = rho*R*T/W */
    }

    return;
}
#endif


/*Compute P = rhoRT/W(c) */
void CKPC(double *  rho, double *  T, double *  c,  double *  P)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    double W = 0;
    double sumC = 0;
    W += c[0]*1.007970; /*H */
    W += c[1]*2.015940; /*H2 */
    W += c[2]*15.999400; /*O */
    W += c[3]*17.007370; /*OH */
    W += c[4]*18.015340; /*H2O */
    W += c[5]*31.998800; /*O2 */
    W += c[6]*33.006770; /*HO2 */
    W += c[7]*34.014740; /*H2O2 */
    W += c[8]*28.013400; /*N2 */
    W += c[9]*39.948000; /*AR */
    W += c[10]*4.002600; /*HE */
    W += c[11]*28.010550; /*CO */
    W += c[12]*44.009950; /*CO2 */

    for (id = 0; id < 13; ++id) {
        sumC += c[id];
    }
    *P = *rho * 8.31446e+07 * (*T) * sumC / W; /*P = rho*R*T/W */

    return;
}


/*Compute rho = PW(x)/RT */
void CKRHOX(double *  P, double *  T, double *  x,  double *  rho)
{
    double XW = 0;/* To hold mean molecular wt */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*2.015940; /*H2 */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*31.998800; /*O2 */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    XW += x[9]*39.948000; /*AR */
    XW += x[10]*4.002600; /*HE */
    XW += x[11]*28.010550; /*CO */
    XW += x[12]*44.009950; /*CO2 */
    *rho = *P * XW / (8.31446e+07 * (*T)); /*rho = P*W/(R*T) */

    return;
}


/*Compute rho = P*W(y)/RT */
AMREX_GPU_HOST_DEVICE void CKRHOY(double *  P, double *  T, double *  y,  double *  rho)
{
    double YOW = 0;
    double tmp[13];

    for (int i = 0; i < 13; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 13; i++)
    {
        YOW += tmp[i];
    }

    *rho = *P / (8.31446e+07 * (*T) * YOW);/*rho = P*W/(R*T) */
    return;
}


/*Compute rho = P*W(c)/(R*T) */
void CKRHOC(double *  P, double *  T, double *  c,  double *  rho)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    double W = 0;
    double sumC = 0;
    W += c[0]*1.007970; /*H */
    W += c[1]*2.015940; /*H2 */
    W += c[2]*15.999400; /*O */
    W += c[3]*17.007370; /*OH */
    W += c[4]*18.015340; /*H2O */
    W += c[5]*31.998800; /*O2 */
    W += c[6]*33.006770; /*HO2 */
    W += c[7]*34.014740; /*H2O2 */
    W += c[8]*28.013400; /*N2 */
    W += c[9]*39.948000; /*AR */
    W += c[10]*4.002600; /*HE */
    W += c[11]*28.010550; /*CO */
    W += c[12]*44.009950; /*CO2 */

    for (id = 0; id < 13; ++id) {
        sumC += c[id];
    }
    *rho = *P * W / (sumC * (*T) * 8.31446e+07); /*rho = PW/(R*T) */

    return;
}


/*get molecular weight for all species */
void CKWT( double *  wt)
{
    get_mw(wt);
}


/*get atomic weight for all elements */
void CKAWT( double *  awt)
{
    atomicWeight(awt);
}


/*given y[species]: mass fractions */
/*returns mean molecular weight (gm/mole) */
AMREX_GPU_HOST_DEVICE void CKMMWY(double *  y,  double *  wtm)
{
    double YOW = 0;
    double tmp[13];

    for (int i = 0; i < 13; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 13; i++)
    {
        YOW += tmp[i];
    }

    *wtm = 1.0 / YOW;
    return;
}


/*given x[species]: mole fractions */
/*returns mean molecular weight (gm/mole) */
void CKMMWX(double *  x,  double *  wtm)
{
    double XW = 0;/* see Eq 4 in CK Manual */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*2.015940; /*H2 */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*31.998800; /*O2 */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    XW += x[9]*39.948000; /*AR */
    XW += x[10]*4.002600; /*HE */
    XW += x[11]*28.010550; /*CO */
    XW += x[12]*44.009950; /*CO2 */
    *wtm = XW;

    return;
}


/*given c[species]: molar concentration */
/*returns mean molecular weight (gm/mole) */
void CKMMWC(double *  c,  double *  wtm)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    double W = 0;
    double sumC = 0;
    W += c[0]*1.007970; /*H */
    W += c[1]*2.015940; /*H2 */
    W += c[2]*15.999400; /*O */
    W += c[3]*17.007370; /*OH */
    W += c[4]*18.015340; /*H2O */
    W += c[5]*31.998800; /*O2 */
    W += c[6]*33.006770; /*HO2 */
    W += c[7]*34.014740; /*H2O2 */
    W += c[8]*28.013400; /*N2 */
    W += c[9]*39.948000; /*AR */
    W += c[10]*4.002600; /*HE */
    W += c[11]*28.010550; /*CO */
    W += c[12]*44.009950; /*CO2 */

    for (id = 0; id < 13; ++id) {
        sumC += c[id];
    }
    /* CK provides no guard against divison by zero */
    *wtm = W/sumC;

    return;
}


/*convert y[species] (mass fracs) to x[species] (mole fracs) */
AMREX_GPU_HOST_DEVICE void CKYTX(double *  y,  double *  x)
{
    double YOW = 0;
    double tmp[13];

    for (int i = 0; i < 13; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 13; i++)
    {
        YOW += tmp[i];
    }

    double YOWINV = 1.0/YOW;

    for (int i = 0; i < 13; i++)
    {
        x[i] = y[i]*imw[i]*YOWINV;
    }
    return;
}


#ifndef AMREX_USE_CUDA
/*convert y[npoints*species] (mass fracs) to x[npoints*species] (mole fracs) */
void VCKYTX(int *  np, double *  y,  double *  x)
{
    double YOW[*np];
    for (int i=0; i<(*np); i++) {
        YOW[i] = 0.0;
    }

    for (int n=0; n<13; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] = y[n*(*np)+i] * imw[n];
            YOW[i] += x[n*(*np)+i];
        }
    }

    for (int i=0; i<(*np); i++) {
        YOW[i] = 1.0/YOW[i];
    }

    for (int n=0; n<13; n++) {
        for (int i=0; i<(*np); i++) {
            x[n*(*np)+i] *=  YOW[i];
        }
    }
}
#else
/*TODO: remove this on GPU */
void VCKYTX(int *  np, double *  y,  double *  x)
{
}
#endif


/*convert y[species] (mass fracs) to c[species] (molar conc) */
void CKYTCP(double *  P, double *  T, double *  y,  double *  c)
{
    double YOW = 0;
    double PWORT;

    /*Compute inverse of mean molecular wt first */
    for (int i = 0; i < 13; i++)
    {
        c[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 13; i++)
    {
        YOW += c[i];
    }

    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31446e+07 * (*T)); 
    /*Now compute conversion */

    for (int i = 0; i < 13; i++)
    {
        c[i] = PWORT * y[i] * imw[i];
    }
    return;
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
AMREX_GPU_HOST_DEVICE void CKYTCR(double *  rho, double *  T, double *  y,  double *  c)
{
    for (int i = 0; i < 13; i++)
    {
        c[i] = (*rho)  * y[i] * imw[i];
    }
}


/*convert x[species] (mole fracs) to y[species] (mass fracs) */
AMREX_GPU_HOST_DEVICE void CKXTY(double *  x,  double *  y)
{
    double XW = 0; /*See Eq 4, 9 in CK Manual */
    /*Compute mean molecular wt first */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*2.015940; /*H2 */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*31.998800; /*O2 */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    XW += x[9]*39.948000; /*AR */
    XW += x[10]*4.002600; /*HE */
    XW += x[11]*28.010550; /*CO */
    XW += x[12]*44.009950; /*CO2 */
    /*Now compute conversion */
    double XWinv = 1.0/XW;
    y[0] = x[0]*1.007970*XWinv; 
    y[1] = x[1]*2.015940*XWinv; 
    y[2] = x[2]*15.999400*XWinv; 
    y[3] = x[3]*17.007370*XWinv; 
    y[4] = x[4]*18.015340*XWinv; 
    y[5] = x[5]*31.998800*XWinv; 
    y[6] = x[6]*33.006770*XWinv; 
    y[7] = x[7]*34.014740*XWinv; 
    y[8] = x[8]*28.013400*XWinv; 
    y[9] = x[9]*39.948000*XWinv; 
    y[10] = x[10]*4.002600*XWinv; 
    y[11] = x[11]*28.010550*XWinv; 
    y[12] = x[12]*44.009950*XWinv; 

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
void CKXTCP(double *  P, double *  T, double *  x,  double *  c)
{
    int id; /*loop counter */
    double PORT = (*P)/(8.31446e+07 * (*T)); /*P/RT */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 13; ++id) {
        c[id] = x[id]*PORT;
    }

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
void CKXTCR(double *  rho, double *  T, double *  x, double *  c)
{
    int id; /*loop counter */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*2.015940; /*H2 */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*31.998800; /*O2 */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    XW += x[9]*39.948000; /*AR */
    XW += x[10]*4.002600; /*HE */
    XW += x[11]*28.010550; /*CO */
    XW += x[12]*44.009950; /*CO2 */
    ROW = (*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 13; ++id) {
        c[id] = x[id]*ROW;
    }

    return;
}


/*convert c[species] (molar conc) to x[species] (mole fracs) */
void CKCTX(double *  c, double *  x)
{
    int id; /*loop counter */
    double sumC = 0; 

    /*compute sum of c  */
    for (id = 0; id < 13; ++id) {
        sumC += c[id];
    }

    /* See Eq 13  */
    double sumCinv = 1.0/sumC;
    for (id = 0; id < 13; ++id) {
        x[id] = c[id]*sumCinv;
    }

    return;
}


/*convert c[species] (molar conc) to y[species] (mass fracs) */
void CKCTY(double *  c, double *  y)
{
    double CW = 0; /*See Eq 12 in CK Manual */
    /*compute denominator in eq 12 first */
    CW += c[0]*1.007970; /*H */
    CW += c[1]*2.015940; /*H2 */
    CW += c[2]*15.999400; /*O */
    CW += c[3]*17.007370; /*OH */
    CW += c[4]*18.015340; /*H2O */
    CW += c[5]*31.998800; /*O2 */
    CW += c[6]*33.006770; /*HO2 */
    CW += c[7]*34.014740; /*H2O2 */
    CW += c[8]*28.013400; /*N2 */
    CW += c[9]*39.948000; /*AR */
    CW += c[10]*4.002600; /*HE */
    CW += c[11]*28.010550; /*CO */
    CW += c[12]*44.009950; /*CO2 */
    /*Now compute conversion */
    double CWinv = 1.0/CW;
    y[0] = c[0]*1.007970*CWinv; 
    y[1] = c[1]*2.015940*CWinv; 
    y[2] = c[2]*15.999400*CWinv; 
    y[3] = c[3]*17.007370*CWinv; 
    y[4] = c[4]*18.015340*CWinv; 
    y[5] = c[5]*31.998800*CWinv; 
    y[6] = c[6]*33.006770*CWinv; 
    y[7] = c[7]*34.014740*CWinv; 
    y[8] = c[8]*28.013400*CWinv; 
    y[9] = c[9]*39.948000*CWinv; 
    y[10] = c[10]*4.002600*CWinv; 
    y[11] = c[11]*28.010550*CWinv; 
    y[12] = c[12]*44.009950*CWinv; 

    return;
}


/*get Cp/R as a function of T  */
/*for all species (Eq 19) */
void CKCPOR(double *  T, double *  cpor)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpor, tc);
}


/*get H/RT as a function of T  */
/*for all species (Eq 20) */
void CKHORT(double *  T, double *  hort)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEnthalpy(hort, tc);
}


/*get S/R as a function of T  */
/*for all species (Eq 21) */
void CKSOR(double *  T, double *  sor)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sor, tc);
}


/*get specific heat at constant volume as a function  */
/*of T for all species (molar units) */
void CKCVML(double *  T,  double *  cvml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cv_R(cvml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 13; ++id) {
        cvml[id] *= 8.31446e+07;
    }
}


/*get specific heat at constant pressure as a  */
/*function of T for all species (molar units) */
void CKCPML(double *  T,  double *  cpml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 13; ++id) {
        cpml[id] *= 8.31446e+07;
    }
}


/*get internal energy as a function  */
/*of T for all species (molar units) */
void CKUML(double *  T,  double *  uml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31446e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 13; ++id) {
        uml[id] *= RT;
    }
}


/*get enthalpy as a function  */
/*of T for all species (molar units) */
void CKHML(double *  T,  double *  hml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31446e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 13; ++id) {
        hml[id] *= RT;
    }
}


/*get standard-state Gibbs energy as a function  */
/*of T for all species (molar units) */
void CKGML(double *  T,  double *  gml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31446e+07*tT; /*R*T */
    gibbs(gml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 13; ++id) {
        gml[id] *= RT;
    }
}


/*get standard-state Helmholtz free energy as a  */
/*function of T for all species (molar units) */
void CKAML(double *  T,  double *  aml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31446e+07*tT; /*R*T */
    helmholtz(aml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 13; ++id) {
        aml[id] *= RT;
    }
}


/*Returns the standard-state entropies in molar units */
void CKSML(double *  T,  double *  sml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 13; ++id) {
        sml[id] *= 8.31446e+07;
    }
}


/*Returns the specific heats at constant volume */
/*in mass units (Eq. 29) */
AMREX_GPU_HOST_DEVICE void CKCVMS(double *  T,  double *  cvms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cv_R(cvms, tc);
    /*multiply by R/molecularweight */
    cvms[0] *= 8.248720317224957e+07; /*H */
    cvms[1] *= 4.124360158612479e+07; /*H2 */
    cvms[2] *= 5.196734013871295e+06; /*O */
    cvms[3] *= 4.888740950630956e+06; /*OH */
    cvms[4] *= 4.615212712140454e+06; /*H2O */
    cvms[5] *= 2.598367006935648e+06; /*O2 */
    cvms[6] *= 2.519017346487778e+06; /*HO2 */
    cvms[7] *= 2.444370475315478e+06; /*H2O2 */
    cvms[8] *= 2.968030520448514e+06; /*N2 */
    cvms[9] *= 2.081321372322329e+06; /*AR */
    cvms[10] *= 2.077265432007505e+07; /*HE */
    cvms[11] *= 2.968332509769797e+06; /*CO */
    cvms[12] *= 1.889223372931176e+06; /*CO2 */
}


/*Returns the specific heats at constant pressure */
/*in mass units (Eq. 26) */
AMREX_GPU_HOST_DEVICE void CKCPMS(double *  T,  double *  cpms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpms, tc);
    /*multiply by R/molecularweight */
    cpms[0] *= 8.248720317224957e+07; /*H */
    cpms[1] *= 4.124360158612479e+07; /*H2 */
    cpms[2] *= 5.196734013871295e+06; /*O */
    cpms[3] *= 4.888740950630956e+06; /*OH */
    cpms[4] *= 4.615212712140454e+06; /*H2O */
    cpms[5] *= 2.598367006935648e+06; /*O2 */
    cpms[6] *= 2.519017346487778e+06; /*HO2 */
    cpms[7] *= 2.444370475315478e+06; /*H2O2 */
    cpms[8] *= 2.968030520448514e+06; /*N2 */
    cpms[9] *= 2.081321372322329e+06; /*AR */
    cpms[10] *= 2.077265432007505e+07; /*HE */
    cpms[11] *= 2.968332509769797e+06; /*CO */
    cpms[12] *= 1.889223372931176e+06; /*CO2 */
}


/*Returns internal energy in mass units (Eq 30.) */
AMREX_GPU_HOST_DEVICE void CKUMS(double *  T,  double *  ums)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31446e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    for (int i = 0; i < 13; i++)
    {
        ums[i] *= RT*imw[i];
    }
}


/*Returns enthalpy in mass units (Eq 27.) */
AMREX_GPU_HOST_DEVICE void CKHMS(double *  T,  double *  hms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31446e+07*tT; /*R*T */
    speciesEnthalpy(hms, tc);
    for (int i = 0; i < 13; i++)
    {
        hms[i] *= RT*imw[i];
    }
}


#ifndef AMREX_USE_CUDA
/*Returns enthalpy in mass units (Eq 27.) */
void VCKHMS(int *  np, double *  T,  double *  hms)
{
    double tc[5], h[13];

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
    }

    for (int n=0; n<13; n++) {
        for (int i=0; i<(*np); i++) {
            hms[n*(*np)+i] *= 8.31446e+07 * T[i] * imw[n];
        }
    }
}
#else
/*TODO: remove this on GPU */
void VCKHMS(int *  np, double *  T,  double *  hms)
{
}
#endif


/*Returns gibbs in mass units (Eq 31.) */
void CKGMS(double *  T,  double *  gms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31446e+07*tT; /*R*T */
    gibbs(gms, tc);
    for (int i = 0; i < 13; i++)
    {
        gms[i] *= RT*imw[i];
    }
}


/*Returns helmholtz in mass units (Eq 32.) */
void CKAMS(double *  T,  double *  ams)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31446e+07*tT; /*R*T */
    helmholtz(ams, tc);
    for (int i = 0; i < 13; i++)
    {
        ams[i] *= RT*imw[i];
    }
}


/*Returns the entropies in mass units (Eq 28.) */
void CKSMS(double *  T,  double *  sms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sms, tc);
    /*multiply by R/molecularweight */
    sms[0] *= 8.248720317224957e+07; /*H */
    sms[1] *= 4.124360158612479e+07; /*H2 */
    sms[2] *= 5.196734013871295e+06; /*O */
    sms[3] *= 4.888740950630956e+06; /*OH */
    sms[4] *= 4.615212712140454e+06; /*H2O */
    sms[5] *= 2.598367006935648e+06; /*O2 */
    sms[6] *= 2.519017346487778e+06; /*HO2 */
    sms[7] *= 2.444370475315478e+06; /*H2O2 */
    sms[8] *= 2.968030520448514e+06; /*N2 */
    sms[9] *= 2.081321372322329e+06; /*AR */
    sms[10] *= 2.077265432007505e+07; /*HE */
    sms[11] *= 2.968332509769797e+06; /*CO */
    sms[12] *= 1.889223372931176e+06; /*CO2 */
}


/*Returns the mean specific heat at CP (Eq. 33) */
void CKCPBL(double *  T, double *  x,  double *  cpbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[13]; /* temporary storage */
    cp_R(cpor, tc);

    /*perform dot product */
    for (id = 0; id < 13; ++id) {
        result += x[id]*cpor[id];
    }

    *cpbl = result * 8.31446e+07;
}


/*Returns the mean specific heat at CP (Eq. 34) */
AMREX_GPU_HOST_DEVICE void CKCPBS(double *  T, double *  y,  double *  cpbs)
{
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[13], tresult[13]; /* temporary storage */
    cp_R(cpor, tc);
    for (int i = 0; i < 13; i++)
    {
        tresult[i] = cpor[i]*y[i]*imw[i];

    }
    for (int i = 0; i < 13; i++)
    {
        result += tresult[i];
    }

    *cpbs = result * 8.31446e+07;
}


/*Returns the mean specific heat at CV (Eq. 35) */
void CKCVBL(double *  T, double *  x,  double *  cvbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cvor[13]; /* temporary storage */
    cv_R(cvor, tc);

    /*perform dot product */
    for (id = 0; id < 13; ++id) {
        result += x[id]*cvor[id];
    }

    *cvbl = result * 8.31446e+07;
}


/*Returns the mean specific heat at CV (Eq. 36) */
AMREX_GPU_HOST_DEVICE void CKCVBS(double *  T, double *  y,  double *  cvbs)
{
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cvor[13]; /* temporary storage */
    cv_R(cvor, tc);
    /*multiply by y/molecularweight */
    result += cvor[0]*y[0]*imw[0]; /*H */
    result += cvor[1]*y[1]*imw[1]; /*H2 */
    result += cvor[2]*y[2]*imw[2]; /*O */
    result += cvor[3]*y[3]*imw[3]; /*OH */
    result += cvor[4]*y[4]*imw[4]; /*H2O */
    result += cvor[5]*y[5]*imw[5]; /*O2 */
    result += cvor[6]*y[6]*imw[6]; /*HO2 */
    result += cvor[7]*y[7]*imw[7]; /*H2O2 */
    result += cvor[8]*y[8]*imw[8]; /*N2 */
    result += cvor[9]*y[9]*imw[9]; /*AR */
    result += cvor[10]*y[10]*imw[10]; /*HE */
    result += cvor[11]*y[11]*imw[11]; /*CO */
    result += cvor[12]*y[12]*imw[12]; /*CO2 */

    *cvbs = result * 8.31446e+07;
}


/*Returns the mean enthalpy of the mixture in molar units */
void CKHBML(double *  T, double *  x,  double *  hbml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[13]; /* temporary storage */
    double RT = 8.31446e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*perform dot product */
    for (id = 0; id < 13; ++id) {
        result += x[id]*hml[id];
    }

    *hbml = result * RT;
}


/*Returns mean enthalpy of mixture in mass units */
AMREX_GPU_HOST_DEVICE void CKHBMS(double *  T, double *  y,  double *  hbms)
{
    double result = 0;
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[13], tmp[13]; /* temporary storage */
    double RT = 8.31446e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);
    int id;
    for (id = 0; id < 13; ++id) {
        tmp[id] = y[id]*hml[id]*imw[id];
    }
    for (id = 0; id < 13; ++id) {
        result += tmp[id];
    }

    *hbms = result * RT;
}


/*get mean internal energy in molar units */
void CKUBML(double *  T, double *  x,  double *  ubml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double uml[13]; /* temporary energy array */
    double RT = 8.31446e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*perform dot product */
    for (id = 0; id < 13; ++id) {
        result += x[id]*uml[id];
    }

    *ubml = result * RT;
}


/*get mean internal energy in mass units */
AMREX_GPU_HOST_DEVICE void CKUBMS(double *  T, double *  y,  double *  ubms)
{
    double result = 0;
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double ums[13]; /* temporary energy array */
    double RT = 8.31446e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    /*perform dot product + scaling by wt */
    result += y[0]*ums[0]*imw[0]; /*H */
    result += y[1]*ums[1]*imw[1]; /*H2 */
    result += y[2]*ums[2]*imw[2]; /*O */
    result += y[3]*ums[3]*imw[3]; /*OH */
    result += y[4]*ums[4]*imw[4]; /*H2O */
    result += y[5]*ums[5]*imw[5]; /*O2 */
    result += y[6]*ums[6]*imw[6]; /*HO2 */
    result += y[7]*ums[7]*imw[7]; /*H2O2 */
    result += y[8]*ums[8]*imw[8]; /*N2 */
    result += y[9]*ums[9]*imw[9]; /*AR */
    result += y[10]*ums[10]*imw[10]; /*HE */
    result += y[11]*ums[11]*imw[11]; /*CO */
    result += y[12]*ums[12]*imw[12]; /*CO2 */

    *ubms = result * RT;
}


/*get mixture entropy in molar units */
void CKSBML(double *  P, double *  T, double *  x,  double *  sbml)
{
    int id; /*loop counter */
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double sor[13]; /* temporary storage */
    speciesEntropy(sor, tc);

    /*Compute Eq 42 */
    for (id = 0; id < 13; ++id) {
        result += x[id]*(sor[id]-log((x[id]+1e-100))-logPratio);
    }

    *sbml = result * 8.31446e+07;
}


/*get mixture entropy in mass units */
void CKSBMS(double *  P, double *  T, double *  y,  double *  sbms)
{
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double sor[13]; /* temporary storage */
    double x[13]; /* need a ytx conversion */
    double YOW = 0; /*See Eq 4, 6 in CK Manual */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H */
    YOW += y[1]*imw[1]; /*H2 */
    YOW += y[2]*imw[2]; /*O */
    YOW += y[3]*imw[3]; /*OH */
    YOW += y[4]*imw[4]; /*H2O */
    YOW += y[5]*imw[5]; /*O2 */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*N2 */
    YOW += y[9]*imw[9]; /*AR */
    YOW += y[10]*imw[10]; /*HE */
    YOW += y[11]*imw[11]; /*CO */
    YOW += y[12]*imw[12]; /*CO2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(1.007970*YOW); 
    x[1] = y[1]/(2.015940*YOW); 
    x[2] = y[2]/(15.999400*YOW); 
    x[3] = y[3]/(17.007370*YOW); 
    x[4] = y[4]/(18.015340*YOW); 
    x[5] = y[5]/(31.998800*YOW); 
    x[6] = y[6]/(33.006770*YOW); 
    x[7] = y[7]/(34.014740*YOW); 
    x[8] = y[8]/(28.013400*YOW); 
    x[9] = y[9]/(39.948000*YOW); 
    x[10] = y[10]/(4.002600*YOW); 
    x[11] = y[11]/(28.010550*YOW); 
    x[12] = y[12]/(44.009950*YOW); 
    speciesEntropy(sor, tc);
    /*Perform computation in Eq 42 and 43 */
    result += x[0]*(sor[0]-log((x[0]+1e-100))-logPratio);
    result += x[1]*(sor[1]-log((x[1]+1e-100))-logPratio);
    result += x[2]*(sor[2]-log((x[2]+1e-100))-logPratio);
    result += x[3]*(sor[3]-log((x[3]+1e-100))-logPratio);
    result += x[4]*(sor[4]-log((x[4]+1e-100))-logPratio);
    result += x[5]*(sor[5]-log((x[5]+1e-100))-logPratio);
    result += x[6]*(sor[6]-log((x[6]+1e-100))-logPratio);
    result += x[7]*(sor[7]-log((x[7]+1e-100))-logPratio);
    result += x[8]*(sor[8]-log((x[8]+1e-100))-logPratio);
    result += x[9]*(sor[9]-log((x[9]+1e-100))-logPratio);
    result += x[10]*(sor[10]-log((x[10]+1e-100))-logPratio);
    result += x[11]*(sor[11]-log((x[11]+1e-100))-logPratio);
    result += x[12]*(sor[12]-log((x[12]+1e-100))-logPratio);
    /*Scale by R/W */
    *sbms = result * 8.31446e+07 * YOW;
}


/*Returns mean gibbs free energy in molar units */
void CKGBML(double *  P, double *  T, double *  x,  double *  gbml)
{
    int id; /*loop counter */
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31446e+07*tT; /*R*T */
    double gort[13]; /* temporary storage */
    /*Compute g/RT */
    gibbs(gort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 13; ++id) {
        result += x[id]*(gort[id]+log((x[id]+1e-100))+logPratio);
    }

    *gbml = result * RT;
}


/*Returns mixture gibbs free energy in mass units */
void CKGBMS(double *  P, double *  T, double *  y,  double *  gbms)
{
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31446e+07*tT; /*R*T */
    double gort[13]; /* temporary storage */
    double x[13]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H */
    YOW += y[1]*imw[1]; /*H2 */
    YOW += y[2]*imw[2]; /*O */
    YOW += y[3]*imw[3]; /*OH */
    YOW += y[4]*imw[4]; /*H2O */
    YOW += y[5]*imw[5]; /*O2 */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*N2 */
    YOW += y[9]*imw[9]; /*AR */
    YOW += y[10]*imw[10]; /*HE */
    YOW += y[11]*imw[11]; /*CO */
    YOW += y[12]*imw[12]; /*CO2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(1.007970*YOW); 
    x[1] = y[1]/(2.015940*YOW); 
    x[2] = y[2]/(15.999400*YOW); 
    x[3] = y[3]/(17.007370*YOW); 
    x[4] = y[4]/(18.015340*YOW); 
    x[5] = y[5]/(31.998800*YOW); 
    x[6] = y[6]/(33.006770*YOW); 
    x[7] = y[7]/(34.014740*YOW); 
    x[8] = y[8]/(28.013400*YOW); 
    x[9] = y[9]/(39.948000*YOW); 
    x[10] = y[10]/(4.002600*YOW); 
    x[11] = y[11]/(28.010550*YOW); 
    x[12] = y[12]/(44.009950*YOW); 
    gibbs(gort, tc);
    /*Perform computation in Eq 44 */
    result += x[0]*(gort[0]+log((x[0]+1e-100))+logPratio);
    result += x[1]*(gort[1]+log((x[1]+1e-100))+logPratio);
    result += x[2]*(gort[2]+log((x[2]+1e-100))+logPratio);
    result += x[3]*(gort[3]+log((x[3]+1e-100))+logPratio);
    result += x[4]*(gort[4]+log((x[4]+1e-100))+logPratio);
    result += x[5]*(gort[5]+log((x[5]+1e-100))+logPratio);
    result += x[6]*(gort[6]+log((x[6]+1e-100))+logPratio);
    result += x[7]*(gort[7]+log((x[7]+1e-100))+logPratio);
    result += x[8]*(gort[8]+log((x[8]+1e-100))+logPratio);
    result += x[9]*(gort[9]+log((x[9]+1e-100))+logPratio);
    result += x[10]*(gort[10]+log((x[10]+1e-100))+logPratio);
    result += x[11]*(gort[11]+log((x[11]+1e-100))+logPratio);
    result += x[12]*(gort[12]+log((x[12]+1e-100))+logPratio);
    /*Scale by RT/W */
    *gbms = result * RT * YOW;
}


/*Returns mean helmholtz free energy in molar units */
void CKABML(double *  P, double *  T, double *  x,  double *  abml)
{
    int id; /*loop counter */
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31446e+07*tT; /*R*T */
    double aort[13]; /* temporary storage */
    /*Compute g/RT */
    helmholtz(aort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 13; ++id) {
        result += x[id]*(aort[id]+log((x[id]+1e-100))+logPratio);
    }

    *abml = result * RT;
}


/*Returns mixture helmholtz free energy in mass units */
void CKABMS(double *  P, double *  T, double *  y,  double *  abms)
{
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31446e+07*tT; /*R*T */
    double aort[13]; /* temporary storage */
    double x[13]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H */
    YOW += y[1]*imw[1]; /*H2 */
    YOW += y[2]*imw[2]; /*O */
    YOW += y[3]*imw[3]; /*OH */
    YOW += y[4]*imw[4]; /*H2O */
    YOW += y[5]*imw[5]; /*O2 */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*N2 */
    YOW += y[9]*imw[9]; /*AR */
    YOW += y[10]*imw[10]; /*HE */
    YOW += y[11]*imw[11]; /*CO */
    YOW += y[12]*imw[12]; /*CO2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(1.007970*YOW); 
    x[1] = y[1]/(2.015940*YOW); 
    x[2] = y[2]/(15.999400*YOW); 
    x[3] = y[3]/(17.007370*YOW); 
    x[4] = y[4]/(18.015340*YOW); 
    x[5] = y[5]/(31.998800*YOW); 
    x[6] = y[6]/(33.006770*YOW); 
    x[7] = y[7]/(34.014740*YOW); 
    x[8] = y[8]/(28.013400*YOW); 
    x[9] = y[9]/(39.948000*YOW); 
    x[10] = y[10]/(4.002600*YOW); 
    x[11] = y[11]/(28.010550*YOW); 
    x[12] = y[12]/(44.009950*YOW); 
    helmholtz(aort, tc);
    /*Perform computation in Eq 44 */
    result += x[0]*(aort[0]+log((x[0]+1e-100))+logPratio);
    result += x[1]*(aort[1]+log((x[1]+1e-100))+logPratio);
    result += x[2]*(aort[2]+log((x[2]+1e-100))+logPratio);
    result += x[3]*(aort[3]+log((x[3]+1e-100))+logPratio);
    result += x[4]*(aort[4]+log((x[4]+1e-100))+logPratio);
    result += x[5]*(aort[5]+log((x[5]+1e-100))+logPratio);
    result += x[6]*(aort[6]+log((x[6]+1e-100))+logPratio);
    result += x[7]*(aort[7]+log((x[7]+1e-100))+logPratio);
    result += x[8]*(aort[8]+log((x[8]+1e-100))+logPratio);
    result += x[9]*(aort[9]+log((x[9]+1e-100))+logPratio);
    result += x[10]*(aort[10]+log((x[10]+1e-100))+logPratio);
    result += x[11]*(aort[11]+log((x[11]+1e-100))+logPratio);
    result += x[12]*(aort[12]+log((x[12]+1e-100))+logPratio);
    /*Scale by RT/W */
    *abms = result * RT * YOW;
}


/*compute the production rate for each species */
AMREX_GPU_HOST_DEVICE void CKWC(double *  T, double *  C,  double *  wdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 13; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    productionRate(wdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 13; ++id) {
        C[id] *= 1.0e-6;
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mass fractions */
void CKWYP(double *  P, double *  T, double *  y,  double *  wdot)
{
    int id; /*loop counter */
    double c[13]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H */
    YOW += y[1]*imw[1]; /*H2 */
    YOW += y[2]*imw[2]; /*O */
    YOW += y[3]*imw[3]; /*OH */
    YOW += y[4]*imw[4]; /*H2O */
    YOW += y[5]*imw[5]; /*O2 */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*N2 */
    YOW += y[9]*imw[9]; /*AR */
    YOW += y[10]*imw[10]; /*HE */
    YOW += y[11]*imw[11]; /*CO */
    YOW += y[12]*imw[12]; /*CO2 */
    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31446e+07 * (*T)); 
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

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 13; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mole fractions */
void CKWXP(double *  P, double *  T, double *  x,  double *  wdot)
{
    int id; /*loop counter */
    double c[13]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31446e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 13; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 13; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
AMREX_GPU_HOST_DEVICE void CKWYR(double *  rho, double *  T, double *  y,  double *  wdot)
{
    int id; /*loop counter */
    double c[13]; /*temporary storage */
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

    /*call productionRate */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 13; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
void VCKWYR(int *  np, double *  rho, double *  T,
	    double *  y,
	    double *  wdot)
{
#ifndef AMREX_USE_CUDA
    double c[13*(*np)]; /*temporary storage */
    /*See Eq 8 with an extra 1e6 so c goes to SI */
    for (int n=0; n<13; n++) {
        for (int i=0; i<(*np); i++) {
            c[n*(*np)+i] = 1.0e6 * rho[i] * y[n*(*np)+i] * imw[n];
        }
    }

    /*call productionRate */
    vproductionRate(*np, wdot, c, T);

    /*convert to chemkin units */
    for (int i=0; i<13*(*np); i++) {
        wdot[i] *= 1.0e-6;
    }
#endif
}


/*Returns the molar production rate of species */
/*Given rho, T, and mole fractions */
void CKWXR(double *  rho, double *  T, double *  x,  double *  wdot)
{
    int id; /*loop counter */
    double c[13]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*2.015940; /*H2 */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*31.998800; /*O2 */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    XW += x[9]*39.948000; /*AR */
    XW += x[10]*4.002600; /*HE */
    XW += x[11]*28.010550; /*CO */
    XW += x[12]*44.009950; /*CO2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 13; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 13; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the rate of progress for each reaction */
void CKQC(double *  T, double *  C, double *  qdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 13; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    progressRate(qdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 13; ++id) {
        C[id] *= 1.0e-6;
    }

    for (id = 0; id < 27; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKKFKR(double *  P, double *  T, double *  x, double *  q_f, double *  q_r)
{
    int id; /*loop counter */
    double c[13]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31446e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 13; ++id) {
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


/*Returns the progress rates of each reactions */
/*Given P, T, and mass fractions */
void CKQYP(double *  P, double *  T, double *  y, double *  qdot)
{
    int id; /*loop counter */
    double c[13]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H */
    YOW += y[1]*imw[1]; /*H2 */
    YOW += y[2]*imw[2]; /*O */
    YOW += y[3]*imw[3]; /*OH */
    YOW += y[4]*imw[4]; /*H2O */
    YOW += y[5]*imw[5]; /*O2 */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*N2 */
    YOW += y[9]*imw[9]; /*AR */
    YOW += y[10]*imw[10]; /*HE */
    YOW += y[11]*imw[11]; /*CO */
    YOW += y[12]*imw[12]; /*CO2 */
    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31446e+07 * (*T)); 
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

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 27; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given P, T, and mole fractions */
void CKQXP(double *  P, double *  T, double *  x, double *  qdot)
{
    int id; /*loop counter */
    double c[13]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31446e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 13; ++id) {
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
void CKQYR(double *  rho, double *  T, double *  y, double *  qdot)
{
    int id; /*loop counter */
    double c[13]; /*temporary storage */
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

    /*call progressRate */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 27; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the progress rates of each reactions */
/*Given rho, T, and mole fractions */
void CKQXR(double *  rho, double *  T, double *  x, double *  qdot)
{
    int id; /*loop counter */
    double c[13]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*1.007970; /*H */
    XW += x[1]*2.015940; /*H2 */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*17.007370; /*OH */
    XW += x[4]*18.015340; /*H2O */
    XW += x[5]*31.998800; /*O2 */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    XW += x[9]*39.948000; /*AR */
    XW += x[10]*4.002600; /*HE */
    XW += x[11]*28.010550; /*CO */
    XW += x[12]*44.009950; /*CO2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 13; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    progressRate(qdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 27; ++id) {
        qdot[id] *= 1.0e-6;
    }
}


/*Returns the stoichiometric coefficients */
/*of the reaction mechanism. (Eq 50) */
void CKNU(int * kdim,  int * nuki)
{
    int id; /*loop counter */
    int kd = (*kdim); 
    /*Zero nuki */
    for (id = 0; id < 13 * kd; ++ id) {
         nuki[id] = 0; 
    }

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    nuki[ 0 * kd + 0 ] += -1.000000 ;
    nuki[ 5 * kd + 0 ] += -1.000000 ;
    nuki[ 6 * kd + 0 ] += +1.000000 ;

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    nuki[ 7 * kd + 1 ] += -1.000000 ;
    nuki[ 3 * kd + 1 ] += +1.000000 ;
    nuki[ 3 * kd + 1 ] += +1.000000 ;

    /*reaction 3: H2 + M <=> H + H + M */
    nuki[ 1 * kd + 2 ] += -1.000000 ;
    nuki[ 0 * kd + 2 ] += +1.000000 ;
    nuki[ 0 * kd + 2 ] += +1.000000 ;

    /*reaction 4: O + O + M <=> O2 + M */
    nuki[ 2 * kd + 3 ] += -1.000000 ;
    nuki[ 2 * kd + 3 ] += -1.000000 ;
    nuki[ 5 * kd + 3 ] += +1.000000 ;

    /*reaction 5: O + H + M <=> OH + M */
    nuki[ 2 * kd + 4 ] += -1.000000 ;
    nuki[ 0 * kd + 4 ] += -1.000000 ;
    nuki[ 3 * kd + 4 ] += +1.000000 ;

    /*reaction 6: H2O + M <=> H + OH + M */
    nuki[ 4 * kd + 5 ] += -1.000000 ;
    nuki[ 0 * kd + 5 ] += +1.000000 ;
    nuki[ 3 * kd + 5 ] += +1.000000 ;

    /*reaction 7: H + O2 <=> O + OH */
    nuki[ 0 * kd + 6 ] += -1.000000 ;
    nuki[ 5 * kd + 6 ] += -1.000000 ;
    nuki[ 2 * kd + 6 ] += +1.000000 ;
    nuki[ 3 * kd + 6 ] += +1.000000 ;

    /*reaction 8: O + H2 <=> H + OH */
    nuki[ 2 * kd + 7 ] += -1.000000 ;
    nuki[ 1 * kd + 7 ] += -1.000000 ;
    nuki[ 0 * kd + 7 ] += +1.000000 ;
    nuki[ 3 * kd + 7 ] += +1.000000 ;

    /*reaction 9: O + H2 <=> H + OH */
    nuki[ 2 * kd + 8 ] += -1.000000 ;
    nuki[ 1 * kd + 8 ] += -1.000000 ;
    nuki[ 0 * kd + 8 ] += +1.000000 ;
    nuki[ 3 * kd + 8 ] += +1.000000 ;

    /*reaction 10: H2 + OH <=> H2O + H */
    nuki[ 1 * kd + 9 ] += -1.000000 ;
    nuki[ 3 * kd + 9 ] += -1.000000 ;
    nuki[ 4 * kd + 9 ] += +1.000000 ;
    nuki[ 0 * kd + 9 ] += +1.000000 ;

    /*reaction 11: OH + OH <=> O + H2O */
    nuki[ 3 * kd + 10 ] += -1.000000 ;
    nuki[ 3 * kd + 10 ] += -1.000000 ;
    nuki[ 2 * kd + 10 ] += +1.000000 ;
    nuki[ 4 * kd + 10 ] += +1.000000 ;

    /*reaction 12: H2 + AR <=> H + H + AR */
    nuki[ 1 * kd + 11 ] += -1.000000 ;
    nuki[ 9 * kd + 11 ] += -1.000000 ;
    nuki[ 0 * kd + 11 ] += +1.000000 ;
    nuki[ 0 * kd + 11 ] += +1.000000 ;
    nuki[ 9 * kd + 11 ] += +1.000000 ;

    /*reaction 13: H2 + HE <=> H + H + HE */
    nuki[ 1 * kd + 12 ] += -1.000000 ;
    nuki[ 10 * kd + 12 ] += -1.000000 ;
    nuki[ 0 * kd + 12 ] += +1.000000 ;
    nuki[ 0 * kd + 12 ] += +1.000000 ;
    nuki[ 10 * kd + 12 ] += +1.000000 ;

    /*reaction 14: O + O + AR <=> O2 + AR */
    nuki[ 2 * kd + 13 ] += -1.000000 ;
    nuki[ 2 * kd + 13 ] += -1.000000 ;
    nuki[ 9 * kd + 13 ] += -1.000000 ;
    nuki[ 5 * kd + 13 ] += +1.000000 ;
    nuki[ 9 * kd + 13 ] += +1.000000 ;

    /*reaction 15: O + O + HE <=> O2 + HE */
    nuki[ 2 * kd + 14 ] += -1.000000 ;
    nuki[ 2 * kd + 14 ] += -1.000000 ;
    nuki[ 10 * kd + 14 ] += -1.000000 ;
    nuki[ 5 * kd + 14 ] += +1.000000 ;
    nuki[ 10 * kd + 14 ] += +1.000000 ;

    /*reaction 16: H2O + H2O <=> H + OH + H2O */
    nuki[ 4 * kd + 15 ] += -1.000000 ;
    nuki[ 4 * kd + 15 ] += -1.000000 ;
    nuki[ 0 * kd + 15 ] += +1.000000 ;
    nuki[ 3 * kd + 15 ] += +1.000000 ;
    nuki[ 4 * kd + 15 ] += +1.000000 ;

    /*reaction 17: HO2 + H <=> H2 + O2 */
    nuki[ 6 * kd + 16 ] += -1.000000 ;
    nuki[ 0 * kd + 16 ] += -1.000000 ;
    nuki[ 1 * kd + 16 ] += +1.000000 ;
    nuki[ 5 * kd + 16 ] += +1.000000 ;

    /*reaction 18: HO2 + H <=> OH + OH */
    nuki[ 6 * kd + 17 ] += -1.000000 ;
    nuki[ 0 * kd + 17 ] += -1.000000 ;
    nuki[ 3 * kd + 17 ] += +1.000000 ;
    nuki[ 3 * kd + 17 ] += +1.000000 ;

    /*reaction 19: HO2 + O <=> O2 + OH */
    nuki[ 6 * kd + 18 ] += -1.000000 ;
    nuki[ 2 * kd + 18 ] += -1.000000 ;
    nuki[ 5 * kd + 18 ] += +1.000000 ;
    nuki[ 3 * kd + 18 ] += +1.000000 ;

    /*reaction 20: HO2 + OH <=> H2O + O2 */
    nuki[ 6 * kd + 19 ] += -1.000000 ;
    nuki[ 3 * kd + 19 ] += -1.000000 ;
    nuki[ 4 * kd + 19 ] += +1.000000 ;
    nuki[ 5 * kd + 19 ] += +1.000000 ;

    /*reaction 21: HO2 + HO2 <=> H2O2 + O2 */
    nuki[ 6 * kd + 20 ] += -1.000000 ;
    nuki[ 6 * kd + 20 ] += -1.000000 ;
    nuki[ 7 * kd + 20 ] += +1.000000 ;
    nuki[ 5 * kd + 20 ] += +1.000000 ;

    /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
    nuki[ 6 * kd + 21 ] += -1.000000 ;
    nuki[ 6 * kd + 21 ] += -1.000000 ;
    nuki[ 7 * kd + 21 ] += +1.000000 ;
    nuki[ 5 * kd + 21 ] += +1.000000 ;

    /*reaction 23: H2O2 + H <=> H2O + OH */
    nuki[ 7 * kd + 22 ] += -1.000000 ;
    nuki[ 0 * kd + 22 ] += -1.000000 ;
    nuki[ 4 * kd + 22 ] += +1.000000 ;
    nuki[ 3 * kd + 22 ] += +1.000000 ;

    /*reaction 24: H2O2 + H <=> HO2 + H2 */
    nuki[ 7 * kd + 23 ] += -1.000000 ;
    nuki[ 0 * kd + 23 ] += -1.000000 ;
    nuki[ 6 * kd + 23 ] += +1.000000 ;
    nuki[ 1 * kd + 23 ] += +1.000000 ;

    /*reaction 25: H2O2 + O <=> OH + HO2 */
    nuki[ 7 * kd + 24 ] += -1.000000 ;
    nuki[ 2 * kd + 24 ] += -1.000000 ;
    nuki[ 3 * kd + 24 ] += +1.000000 ;
    nuki[ 6 * kd + 24 ] += +1.000000 ;

    /*reaction 26: H2O2 + OH <=> HO2 + H2O */
    nuki[ 7 * kd + 25 ] += -1.000000 ;
    nuki[ 3 * kd + 25 ] += -1.000000 ;
    nuki[ 6 * kd + 25 ] += +1.000000 ;
    nuki[ 4 * kd + 25 ] += +1.000000 ;

    /*reaction 27: H2O2 + OH <=> HO2 + H2O */
    nuki[ 7 * kd + 26 ] += -1.000000 ;
    nuki[ 3 * kd + 26 ] += -1.000000 ;
    nuki[ 6 * kd + 26 ] += +1.000000 ;
    nuki[ 4 * kd + 26 ] += +1.000000 ;
}


#ifndef AMREX_USE_CUDA
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
#endif


/*Returns the elemental composition  */
/*of the speciesi (mdim is num of elements) */
void CKNCF(int * ncf)
{
    int id; /*loop counter */
    int kd = 6; 
    /*Zero ncf */
    for (id = 0; id < kd * 13; ++ id) {
         ncf[id] = 0; 
    }

    /*H */
    ncf[ 0 * kd + 0 ] = 1; /*H */

    /*H2 */
    ncf[ 1 * kd + 0 ] = 2; /*H */

    /*O */
    ncf[ 2 * kd + 1 ] = 1; /*O */

    /*OH */
    ncf[ 3 * kd + 1 ] = 1; /*O */
    ncf[ 3 * kd + 0 ] = 1; /*H */

    /*H2O */
    ncf[ 4 * kd + 0 ] = 2; /*H */
    ncf[ 4 * kd + 1 ] = 1; /*O */

    /*O2 */
    ncf[ 5 * kd + 1 ] = 2; /*O */

    /*HO2 */
    ncf[ 6 * kd + 0 ] = 1; /*H */
    ncf[ 6 * kd + 1 ] = 2; /*O */

    /*H2O2 */
    ncf[ 7 * kd + 0 ] = 2; /*H */
    ncf[ 7 * kd + 1 ] = 2; /*O */

    /*N2 */
    ncf[ 8 * kd + 2 ] = 2; /*N */

    /*AR */
    ncf[ 9 * kd + 3 ] = 1; /*AR */

    /*HE */
    ncf[ 10 * kd + 4 ] = 1; /*HE */

    /*CO */
    ncf[ 11 * kd + 5 ] = 1; /*C */
    ncf[ 11 * kd + 1 ] = 1; /*O */

    /*CO2 */
    ncf[ 12 * kd + 5 ] = 1; /*C */
    ncf[ 12 * kd + 1 ] = 2; /*O */

}


/*Returns the arrehenius coefficients  */
/*for all reactions */
void CKABE( double *  a, double *  b, double *  e)
{
    // (14):  H + O2 (+M) <=> HO2 (+M)
    a[0] = 4650840000000;
    b[0] = 0.44;
    e[0] = 0;

    // (21):  H2O2 (+M) <=> OH + OH (+M)
    a[1] = 2000000000000;
    b[1] = 0.90000000000000002;
    e[1] = 48749;

    // (5):  H2 + M <=> H + H + M
    a[2] = 4.577e+19;
    b[2] = -1.3999999999999999;
    e[2] = 104380;

    // (8):  O + O + M <=> O2 + M
    a[3] = 6165000000000000;
    b[3] = -0.5;
    e[3] = 0;

    // (11):  O + H + M <=> OH + M
    a[4] = 4.714e+18;
    b[4] = -1;
    e[4] = 0;

    // (12):  H2O + M <=> H + OH + M
    a[5] = 6.0640000000000002e+27;
    b[5] = -3.3220000000000001;
    e[5] = 120790;

    // (0):  H + O2 <=> O + OH
    a[6] = 104000000000000;
    b[6] = 0;
    e[6] = 15286;

    // (1):  O + H2 <=> H + OH
    a[7] = 3818000000000;
    b[7] = 0;
    e[7] = 7948;

    // (2):  O + H2 <=> H + OH
    a[8] = 879200000000000;
    b[8] = 0;
    e[8] = 19170;

    // (3):  H2 + OH <=> H2O + H
    a[9] = 216000000;
    b[9] = 1.51;
    e[9] = 3430;

    // (4):  OH + OH <=> O + H2O
    a[10] = 33400;
    b[10] = 2.4199999999999999;
    e[10] = -1930;

    // (6):  H2 + AR <=> H + H + AR
    a[11] = 5.84e+18;
    b[11] = -1.1000000000000001;
    e[11] = 104380;

    // (7):  H2 + HE <=> H + H + HE
    a[12] = 5.84e+18;
    b[12] = -1.1000000000000001;
    e[12] = 104380;

    // (9):  O + O + AR <=> O2 + AR
    a[13] = 18860000000000;
    b[13] = 0;
    e[13] = -1788;

    // (10):  O + O + HE <=> O2 + HE
    a[14] = 18860000000000;
    b[14] = 0;
    e[14] = -1788;

    // (13):  H2O + H2O <=> H + OH + H2O
    a[15] = 1.006e+26;
    b[15] = -2.4399999999999999;
    e[15] = 120180;

    // (15):  HO2 + H <=> H2 + O2
    a[16] = 2750000;
    b[16] = 2.0899999999999999;
    e[16] = -1451;

    // (16):  HO2 + H <=> OH + OH
    a[17] = 70790000000000;
    b[17] = 0;
    e[17] = 295;

    // (17):  HO2 + O <=> O2 + OH
    a[18] = 28500000000;
    b[18] = 1;
    e[18] = -723.92999999999995;

    // (18):  HO2 + OH <=> H2O + O2
    a[19] = 28900000000000;
    b[19] = 0;
    e[19] = -497;

    // (19):  HO2 + HO2 <=> H2O2 + O2
    a[20] = 420000000000000;
    b[20] = 0;
    e[20] = 11982;

    // (20):  HO2 + HO2 <=> H2O2 + O2
    a[21] = 130000000000;
    b[21] = 0;
    e[21] = -1629.3;

    // (22):  H2O2 + H <=> H2O + OH
    a[22] = 24100000000000;
    b[22] = 0;
    e[22] = 3970;

    // (23):  H2O2 + H <=> HO2 + H2
    a[23] = 48200000000000;
    b[23] = 0;
    e[23] = 7950;

    // (24):  H2O2 + O <=> OH + HO2
    a[24] = 9550000;
    b[24] = 2;
    e[24] = 3970;

    // (25):  H2O2 + OH <=> HO2 + H2O
    a[25] = 1740000000000;
    b[25] = 0;
    e[25] = 318;

    // (26):  H2O2 + OH <=> HO2 + H2O
    a[26] = 75900000000000;
    b[26] = 0;
    e[26] = 7270;


    return;
}


/*Returns the equil constants for each reaction */
void CKEQC(double *  T, double *  C, double *  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[13]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    eqcon[1] *= 1e-06; 

    /*reaction 3: H2 + M <=> H + H + M */
    eqcon[2] *= 1e-06; 

    /*reaction 4: O + O + M <=> O2 + M */
    eqcon[3] *= 1e+06; 

    /*reaction 5: O + H + M <=> OH + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: H2O + M <=> H + OH + M */
    eqcon[5] *= 1e-06; 

    /*reaction 7: H + O2 <=> O + OH */
    /*eqcon[6] *= 1;  */

    /*reaction 8: O + H2 <=> H + OH */
    /*eqcon[7] *= 1;  */

    /*reaction 9: O + H2 <=> H + OH */
    /*eqcon[8] *= 1;  */

    /*reaction 10: H2 + OH <=> H2O + H */
    /*eqcon[9] *= 1;  */

    /*reaction 11: OH + OH <=> O + H2O */
    /*eqcon[10] *= 1;  */

    /*reaction 12: H2 + AR <=> H + H + AR */
    eqcon[11] *= 1e-06; 

    /*reaction 13: H2 + HE <=> H + H + HE */
    eqcon[12] *= 1e-06; 

    /*reaction 14: O + O + AR <=> O2 + AR */
    eqcon[13] *= 1e+06; 

    /*reaction 15: O + O + HE <=> O2 + HE */
    eqcon[14] *= 1e+06; 

    /*reaction 16: H2O + H2O <=> H + OH + H2O */
    eqcon[15] *= 1e-06; 

    /*reaction 17: HO2 + H <=> H2 + O2 */
    /*eqcon[16] *= 1;  */

    /*reaction 18: HO2 + H <=> OH + OH */
    /*eqcon[17] *= 1;  */

    /*reaction 19: HO2 + O <=> O2 + OH */
    /*eqcon[18] *= 1;  */

    /*reaction 20: HO2 + OH <=> H2O + O2 */
    /*eqcon[19] *= 1;  */

    /*reaction 21: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[20] *= 1;  */

    /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[21] *= 1;  */

    /*reaction 23: H2O2 + H <=> H2O + OH */
    /*eqcon[22] *= 1;  */

    /*reaction 24: H2O2 + H <=> HO2 + H2 */
    /*eqcon[23] *= 1;  */

    /*reaction 25: H2O2 + O <=> OH + HO2 */
    /*eqcon[24] *= 1;  */

    /*reaction 26: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[25] *= 1;  */

    /*reaction 27: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[26] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mass fractions */
void CKEQYP(double *  P, double *  T, double *  y, double *  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[13]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    eqcon[1] *= 1e-06; 

    /*reaction 3: H2 + M <=> H + H + M */
    eqcon[2] *= 1e-06; 

    /*reaction 4: O + O + M <=> O2 + M */
    eqcon[3] *= 1e+06; 

    /*reaction 5: O + H + M <=> OH + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: H2O + M <=> H + OH + M */
    eqcon[5] *= 1e-06; 

    /*reaction 7: H + O2 <=> O + OH */
    /*eqcon[6] *= 1;  */

    /*reaction 8: O + H2 <=> H + OH */
    /*eqcon[7] *= 1;  */

    /*reaction 9: O + H2 <=> H + OH */
    /*eqcon[8] *= 1;  */

    /*reaction 10: H2 + OH <=> H2O + H */
    /*eqcon[9] *= 1;  */

    /*reaction 11: OH + OH <=> O + H2O */
    /*eqcon[10] *= 1;  */

    /*reaction 12: H2 + AR <=> H + H + AR */
    eqcon[11] *= 1e-06; 

    /*reaction 13: H2 + HE <=> H + H + HE */
    eqcon[12] *= 1e-06; 

    /*reaction 14: O + O + AR <=> O2 + AR */
    eqcon[13] *= 1e+06; 

    /*reaction 15: O + O + HE <=> O2 + HE */
    eqcon[14] *= 1e+06; 

    /*reaction 16: H2O + H2O <=> H + OH + H2O */
    eqcon[15] *= 1e-06; 

    /*reaction 17: HO2 + H <=> H2 + O2 */
    /*eqcon[16] *= 1;  */

    /*reaction 18: HO2 + H <=> OH + OH */
    /*eqcon[17] *= 1;  */

    /*reaction 19: HO2 + O <=> O2 + OH */
    /*eqcon[18] *= 1;  */

    /*reaction 20: HO2 + OH <=> H2O + O2 */
    /*eqcon[19] *= 1;  */

    /*reaction 21: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[20] *= 1;  */

    /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[21] *= 1;  */

    /*reaction 23: H2O2 + H <=> H2O + OH */
    /*eqcon[22] *= 1;  */

    /*reaction 24: H2O2 + H <=> HO2 + H2 */
    /*eqcon[23] *= 1;  */

    /*reaction 25: H2O2 + O <=> OH + HO2 */
    /*eqcon[24] *= 1;  */

    /*reaction 26: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[25] *= 1;  */

    /*reaction 27: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[26] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mole fractions */
void CKEQXP(double *  P, double *  T, double *  x, double *  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[13]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    eqcon[1] *= 1e-06; 

    /*reaction 3: H2 + M <=> H + H + M */
    eqcon[2] *= 1e-06; 

    /*reaction 4: O + O + M <=> O2 + M */
    eqcon[3] *= 1e+06; 

    /*reaction 5: O + H + M <=> OH + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: H2O + M <=> H + OH + M */
    eqcon[5] *= 1e-06; 

    /*reaction 7: H + O2 <=> O + OH */
    /*eqcon[6] *= 1;  */

    /*reaction 8: O + H2 <=> H + OH */
    /*eqcon[7] *= 1;  */

    /*reaction 9: O + H2 <=> H + OH */
    /*eqcon[8] *= 1;  */

    /*reaction 10: H2 + OH <=> H2O + H */
    /*eqcon[9] *= 1;  */

    /*reaction 11: OH + OH <=> O + H2O */
    /*eqcon[10] *= 1;  */

    /*reaction 12: H2 + AR <=> H + H + AR */
    eqcon[11] *= 1e-06; 

    /*reaction 13: H2 + HE <=> H + H + HE */
    eqcon[12] *= 1e-06; 

    /*reaction 14: O + O + AR <=> O2 + AR */
    eqcon[13] *= 1e+06; 

    /*reaction 15: O + O + HE <=> O2 + HE */
    eqcon[14] *= 1e+06; 

    /*reaction 16: H2O + H2O <=> H + OH + H2O */
    eqcon[15] *= 1e-06; 

    /*reaction 17: HO2 + H <=> H2 + O2 */
    /*eqcon[16] *= 1;  */

    /*reaction 18: HO2 + H <=> OH + OH */
    /*eqcon[17] *= 1;  */

    /*reaction 19: HO2 + O <=> O2 + OH */
    /*eqcon[18] *= 1;  */

    /*reaction 20: HO2 + OH <=> H2O + O2 */
    /*eqcon[19] *= 1;  */

    /*reaction 21: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[20] *= 1;  */

    /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[21] *= 1;  */

    /*reaction 23: H2O2 + H <=> H2O + OH */
    /*eqcon[22] *= 1;  */

    /*reaction 24: H2O2 + H <=> HO2 + H2 */
    /*eqcon[23] *= 1;  */

    /*reaction 25: H2O2 + O <=> OH + HO2 */
    /*eqcon[24] *= 1;  */

    /*reaction 26: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[25] *= 1;  */

    /*reaction 27: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[26] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mass fractions */
void CKEQYR(double *  rho, double *  T, double *  y, double *  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[13]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    eqcon[1] *= 1e-06; 

    /*reaction 3: H2 + M <=> H + H + M */
    eqcon[2] *= 1e-06; 

    /*reaction 4: O + O + M <=> O2 + M */
    eqcon[3] *= 1e+06; 

    /*reaction 5: O + H + M <=> OH + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: H2O + M <=> H + OH + M */
    eqcon[5] *= 1e-06; 

    /*reaction 7: H + O2 <=> O + OH */
    /*eqcon[6] *= 1;  */

    /*reaction 8: O + H2 <=> H + OH */
    /*eqcon[7] *= 1;  */

    /*reaction 9: O + H2 <=> H + OH */
    /*eqcon[8] *= 1;  */

    /*reaction 10: H2 + OH <=> H2O + H */
    /*eqcon[9] *= 1;  */

    /*reaction 11: OH + OH <=> O + H2O */
    /*eqcon[10] *= 1;  */

    /*reaction 12: H2 + AR <=> H + H + AR */
    eqcon[11] *= 1e-06; 

    /*reaction 13: H2 + HE <=> H + H + HE */
    eqcon[12] *= 1e-06; 

    /*reaction 14: O + O + AR <=> O2 + AR */
    eqcon[13] *= 1e+06; 

    /*reaction 15: O + O + HE <=> O2 + HE */
    eqcon[14] *= 1e+06; 

    /*reaction 16: H2O + H2O <=> H + OH + H2O */
    eqcon[15] *= 1e-06; 

    /*reaction 17: HO2 + H <=> H2 + O2 */
    /*eqcon[16] *= 1;  */

    /*reaction 18: HO2 + H <=> OH + OH */
    /*eqcon[17] *= 1;  */

    /*reaction 19: HO2 + O <=> O2 + OH */
    /*eqcon[18] *= 1;  */

    /*reaction 20: HO2 + OH <=> H2O + O2 */
    /*eqcon[19] *= 1;  */

    /*reaction 21: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[20] *= 1;  */

    /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[21] *= 1;  */

    /*reaction 23: H2O2 + H <=> H2O + OH */
    /*eqcon[22] *= 1;  */

    /*reaction 24: H2O2 + H <=> HO2 + H2 */
    /*eqcon[23] *= 1;  */

    /*reaction 25: H2O2 + O <=> OH + HO2 */
    /*eqcon[24] *= 1;  */

    /*reaction 26: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[25] *= 1;  */

    /*reaction 27: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[26] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mole fractions */
void CKEQXR(double *  rho, double *  T, double *  x, double *  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[13]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    eqcon[1] *= 1e-06; 

    /*reaction 3: H2 + M <=> H + H + M */
    eqcon[2] *= 1e-06; 

    /*reaction 4: O + O + M <=> O2 + M */
    eqcon[3] *= 1e+06; 

    /*reaction 5: O + H + M <=> OH + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: H2O + M <=> H + OH + M */
    eqcon[5] *= 1e-06; 

    /*reaction 7: H + O2 <=> O + OH */
    /*eqcon[6] *= 1;  */

    /*reaction 8: O + H2 <=> H + OH */
    /*eqcon[7] *= 1;  */

    /*reaction 9: O + H2 <=> H + OH */
    /*eqcon[8] *= 1;  */

    /*reaction 10: H2 + OH <=> H2O + H */
    /*eqcon[9] *= 1;  */

    /*reaction 11: OH + OH <=> O + H2O */
    /*eqcon[10] *= 1;  */

    /*reaction 12: H2 + AR <=> H + H + AR */
    eqcon[11] *= 1e-06; 

    /*reaction 13: H2 + HE <=> H + H + HE */
    eqcon[12] *= 1e-06; 

    /*reaction 14: O + O + AR <=> O2 + AR */
    eqcon[13] *= 1e+06; 

    /*reaction 15: O + O + HE <=> O2 + HE */
    eqcon[14] *= 1e+06; 

    /*reaction 16: H2O + H2O <=> H + OH + H2O */
    eqcon[15] *= 1e-06; 

    /*reaction 17: HO2 + H <=> H2 + O2 */
    /*eqcon[16] *= 1;  */

    /*reaction 18: HO2 + H <=> OH + OH */
    /*eqcon[17] *= 1;  */

    /*reaction 19: HO2 + O <=> O2 + OH */
    /*eqcon[18] *= 1;  */

    /*reaction 20: HO2 + OH <=> H2O + O2 */
    /*eqcon[19] *= 1;  */

    /*reaction 21: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[20] *= 1;  */

    /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
    /*eqcon[21] *= 1;  */

    /*reaction 23: H2O2 + H <=> H2O + OH */
    /*eqcon[22] *= 1;  */

    /*reaction 24: H2O2 + H <=> HO2 + H2 */
    /*eqcon[23] *= 1;  */

    /*reaction 25: H2O2 + O <=> OH + HO2 */
    /*eqcon[24] *= 1;  */

    /*reaction 26: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[25] *= 1;  */

    /*reaction 27: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[26] *= 1;  */
}

#ifdef AMREX_USE_CUDA
/*GPU version of productionRate: no more use of thermo namespace vectors */
/*compute the production rate for each species */
AMREX_GPU_HOST_DEVICE inline void  productionRate(double * wdot, double * sc, double T)
{
    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    double invT = 1.0 / tc[1];

    double qdot, q_f[27], q_r[27];
    comp_qfqr(q_f, q_r, sc, tc, invT);

    for (int i = 0; i < 13; ++i) {
        wdot[i] = 0.0;
    }

    qdot = q_f[0]-q_r[0];
    wdot[0] -= qdot;
    wdot[5] -= qdot;
    wdot[6] += qdot;

    qdot = q_f[1]-q_r[1];
    wdot[3] += qdot;
    wdot[3] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[2]-q_r[2];
    wdot[0] += qdot;
    wdot[0] += qdot;
    wdot[1] -= qdot;

    qdot = q_f[3]-q_r[3];
    wdot[2] -= qdot;
    wdot[2] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[4]-q_r[4];
    wdot[0] -= qdot;
    wdot[2] -= qdot;
    wdot[3] += qdot;

    qdot = q_f[5]-q_r[5];
    wdot[0] += qdot;
    wdot[3] += qdot;
    wdot[4] -= qdot;

    qdot = q_f[6]-q_r[6];
    wdot[0] -= qdot;
    wdot[2] += qdot;
    wdot[3] += qdot;
    wdot[5] -= qdot;

    qdot = q_f[7]-q_r[7];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[2] -= qdot;
    wdot[3] += qdot;

    qdot = q_f[8]-q_r[8];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[2] -= qdot;
    wdot[3] += qdot;

    qdot = q_f[9]-q_r[9];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[10]-q_r[10];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[3] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[11]-q_r[11];
    wdot[0] += qdot;
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[9] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[12]-q_r[12];
    wdot[0] += qdot;
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[10] -= qdot;
    wdot[10] += qdot;

    qdot = q_f[13]-q_r[13];
    wdot[2] -= qdot;
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[9] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[14]-q_r[14];
    wdot[2] -= qdot;
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[10] -= qdot;
    wdot[10] += qdot;

    qdot = q_f[15]-q_r[15];
    wdot[0] += qdot;
    wdot[3] += qdot;
    wdot[4] -= qdot;
    wdot[4] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[16]-q_r[16];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[5] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[17]-q_r[17];
    wdot[0] -= qdot;
    wdot[3] += qdot;
    wdot[3] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[18]-q_r[18];
    wdot[2] -= qdot;
    wdot[3] += qdot;
    wdot[5] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[19]-q_r[19];
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[5] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[20]-q_r[20];
    wdot[5] += qdot;
    wdot[6] -= qdot;
    wdot[6] -= qdot;
    wdot[7] += qdot;

    qdot = q_f[21]-q_r[21];
    wdot[5] += qdot;
    wdot[6] -= qdot;
    wdot[6] -= qdot;
    wdot[7] += qdot;

    qdot = q_f[22]-q_r[22];
    wdot[0] -= qdot;
    wdot[3] += qdot;
    wdot[4] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[23]-q_r[23];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[24]-q_r[24];
    wdot[2] -= qdot;
    wdot[3] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[25]-q_r[25];
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[26]-q_r[26];
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    return;
}

AMREX_GPU_HOST_DEVICE inline void comp_qfqr(double *  qf, double * qr, double * sc, double * tc, double invT)
{

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    qf[0] = sc[0]*sc[5];
    qr[0] = sc[6];

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    qf[1] = sc[7];
    qr[1] = sc[3]*sc[3];

    /*reaction 3: H2 + M <=> H + H + M */
    qf[2] = sc[1];
    qr[2] = sc[0]*sc[0];

    /*reaction 4: O + O + M <=> O2 + M */
    qf[3] = sc[2]*sc[2];
    qr[3] = sc[5];

    /*reaction 5: O + H + M <=> OH + M */
    qf[4] = sc[0]*sc[2];
    qr[4] = sc[3];

    /*reaction 6: H2O + M <=> H + OH + M */
    qf[5] = sc[4];
    qr[5] = sc[0]*sc[3];

    /*reaction 7: H + O2 <=> O + OH */
    qf[6] = sc[0]*sc[5];
    qr[6] = sc[2]*sc[3];

    /*reaction 8: O + H2 <=> H + OH */
    qf[7] = sc[1]*sc[2];
    qr[7] = sc[0]*sc[3];

    /*reaction 9: O + H2 <=> H + OH */
    qf[8] = sc[1]*sc[2];
    qr[8] = sc[0]*sc[3];

    /*reaction 10: H2 + OH <=> H2O + H */
    qf[9] = sc[1]*sc[3];
    qr[9] = sc[0]*sc[4];

    /*reaction 11: OH + OH <=> O + H2O */
    qf[10] = sc[3]*sc[3];
    qr[10] = sc[2]*sc[4];

    /*reaction 12: H2 + AR <=> H + H + AR */
    qf[11] = sc[1]*sc[9];
    qr[11] = sc[0]*sc[0]*sc[9];

    /*reaction 13: H2 + HE <=> H + H + HE */
    qf[12] = sc[1]*sc[10];
    qr[12] = sc[0]*sc[0]*sc[10];

    /*reaction 14: O + O + AR <=> O2 + AR */
    qf[13] = sc[2]*sc[2]*sc[9];
    qr[13] = sc[5]*sc[9];

    /*reaction 15: O + O + HE <=> O2 + HE */
    qf[14] = sc[2]*sc[2]*sc[10];
    qr[14] = sc[5]*sc[10];

    /*reaction 16: H2O + H2O <=> H + OH + H2O */
    qf[15] = sc[4]*sc[4];
    qr[15] = sc[0]*sc[3]*sc[4];

    /*reaction 17: HO2 + H <=> H2 + O2 */
    qf[16] = sc[0]*sc[6];
    qr[16] = sc[1]*sc[5];

    /*reaction 18: HO2 + H <=> OH + OH */
    qf[17] = sc[0]*sc[6];
    qr[17] = sc[3]*sc[3];

    /*reaction 19: HO2 + O <=> O2 + OH */
    qf[18] = sc[2]*sc[6];
    qr[18] = sc[3]*sc[5];

    /*reaction 20: HO2 + OH <=> H2O + O2 */
    qf[19] = sc[3]*sc[6];
    qr[19] = sc[4]*sc[5];

    /*reaction 21: HO2 + HO2 <=> H2O2 + O2 */
    qf[20] = sc[6]*sc[6];
    qr[20] = sc[5]*sc[7];

    /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
    qf[21] = sc[6]*sc[6];
    qr[21] = sc[5]*sc[7];

    /*reaction 23: H2O2 + H <=> H2O + OH */
    qf[22] = sc[0]*sc[7];
    qr[22] = sc[3]*sc[4];

    /*reaction 24: H2O2 + H <=> HO2 + H2 */
    qf[23] = sc[0]*sc[7];
    qr[23] = sc[1]*sc[6];

    /*reaction 25: H2O2 + O <=> OH + HO2 */
    qf[24] = sc[2]*sc[7];
    qr[24] = sc[3]*sc[6];

    /*reaction 26: H2O2 + OH <=> HO2 + H2O */
    qf[25] = sc[3]*sc[7];
    qr[25] = sc[4]*sc[6];

    /*reaction 27: H2O2 + OH <=> HO2 + H2O */
    qf[26] = sc[3]*sc[7];
    qr[26] = sc[4]*sc[6];

    /*compute the mixture concentration */
    double mixture = 0.0;
    for (int i = 0; i < 13; ++i) {
        mixture += sc[i];
    }

    /*compute the Gibbs free energy */
    double g_RT[13];
    gibbs(g_RT, tc);

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31446 * invT;
    double refCinv = 1 / refC;

    /* Evaluate the kfs */
    double k_f, Corr;
    double redP, F, logPred, logFcent, troe_c, troe_n, troe, F_troe;

    // (0):  H + O2 <=> O + OH
    k_f = 1.0000000000000002e-06 * 104000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (15286) * invT);
    Corr  = 1.0;
    qf[6] *= Corr * k_f;
    qr[6] *= Corr * k_f / exp(g_RT[0] - g_RT[2] - g_RT[3] + g_RT[5]);
    // (1):  O + H2 <=> H + OH
    k_f = 1.0000000000000002e-06 * 3818000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (7948) * invT);
    Corr  = 1.0;
    qf[7] *= Corr * k_f;
    qr[7] *= Corr * k_f / exp(-g_RT[0] + g_RT[1] + g_RT[2] - g_RT[3]);
    // (2):  O + H2 <=> H + OH
    k_f = 1.0000000000000002e-06 * 879200000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (19170) * invT);
    Corr  = 1.0;
    qf[8] *= Corr * k_f;
    qr[8] *= Corr * k_f / exp(-g_RT[0] + g_RT[1] + g_RT[2] - g_RT[3]);
    // (3):  H2 + OH <=> H2O + H
    k_f = 1.0000000000000002e-06 * 216000000 
               * exp(1.51 * tc[0] - 0.50321666580471969 * (3430) * invT);
    Corr  = 1.0;
    qf[9] *= Corr * k_f;
    qr[9] *= Corr * k_f / exp(-g_RT[0] + g_RT[1] + g_RT[3] - g_RT[4]);
    // (4):  OH + OH <=> O + H2O
    k_f = 1.0000000000000002e-06 * 33400 
               * exp(2.4199999999999999 * tc[0] - 0.50321666580471969 * (-1930) * invT);
    Corr  = 1.0;
    qf[10] *= Corr * k_f;
    qr[10] *= Corr * k_f / exp(-g_RT[2] + g_RT[3] + g_RT[3] - g_RT[4]);
    // (5):  H2 + M <=> H + H + M
    k_f = 1.0000000000000002e-06 * 4.577e+19 
               * exp(-1.3999999999999999 * tc[0] - 0.50321666580471969 * (104380) * invT);
    Corr  = mixture + ( 2.5 - 1)*sc[1] + ( 12 - 1)*sc[4] + ( 1.8999999999999999 - 1)*sc[11] + ( 3.7999999999999998 - 1)*sc[12] + ( 0 - 1)*sc[9] + ( 0 - 1)*sc[10];
    qf[2] *= Corr * k_f;
    qr[2] *= Corr * k_f / (exp(-g_RT[0] - g_RT[0] + g_RT[1]) * refC);
    // (6):  H2 + AR <=> H + H + AR
    k_f = 1.0000000000000002e-06 * 5.84e+18 
               * exp(-1.1000000000000001 * tc[0] - 0.50321666580471969 * (104380) * invT);
    Corr  = 1.0;
    qf[11] *= Corr * k_f;
    qr[11] *= Corr * k_f / (exp(-g_RT[0] - g_RT[0] + g_RT[1] + g_RT[9] - g_RT[9]) * refC);
    // (7):  H2 + HE <=> H + H + HE
    k_f = 1.0000000000000002e-06 * 5.84e+18 
               * exp(-1.1000000000000001 * tc[0] - 0.50321666580471969 * (104380) * invT);
    Corr  = 1.0;
    qf[12] *= Corr * k_f;
    qr[12] *= Corr * k_f / (exp(-g_RT[0] - g_RT[0] + g_RT[1] + g_RT[10] - g_RT[10]) * refC);
    // (8):  O + O + M <=> O2 + M
    k_f = 1.0000000000000002e-12 * 6165000000000000 
               * exp(-0.5 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = mixture + ( 2.5 - 1)*sc[1] + ( 12 - 1)*sc[4] + ( 0 - 1)*sc[9] + ( 0 - 1)*sc[10] + ( 1.8999999999999999 - 1)*sc[11] + ( 3.7999999999999998 - 1)*sc[12];
    qf[3] *= Corr * k_f;
    qr[3] *= Corr * k_f / (exp(g_RT[2] + g_RT[2] - g_RT[5]) * refCinv);
    // (9):  O + O + AR <=> O2 + AR
    k_f = 1.0000000000000002e-12 * 18860000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (-1788) * invT);
    Corr  = 1.0;
    qf[13] *= Corr * k_f;
    qr[13] *= Corr * k_f / (exp(g_RT[2] + g_RT[2] - g_RT[5] + g_RT[9] - g_RT[9]) * refCinv);
    // (10):  O + O + HE <=> O2 + HE
    k_f = 1.0000000000000002e-12 * 18860000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (-1788) * invT);
    Corr  = 1.0;
    qf[14] *= Corr * k_f;
    qr[14] *= Corr * k_f / (exp(g_RT[2] + g_RT[2] - g_RT[5] + g_RT[10] - g_RT[10]) * refCinv);
    // (11):  O + H + M <=> OH + M
    k_f = 1.0000000000000002e-12 * 4.714e+18 
               * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = mixture + ( 2.5 - 1)*sc[1] + ( 12 - 1)*sc[4] + ( 0.75 - 1)*sc[9] + ( 0.75 - 1)*sc[10] + ( 1.8999999999999999 - 1)*sc[11] + ( 3.7999999999999998 - 1)*sc[12];
    qf[4] *= Corr * k_f;
    qr[4] *= Corr * k_f / (exp(g_RT[0] + g_RT[2] - g_RT[3]) * refCinv);
    // (12):  H2O + M <=> H + OH + M
    k_f = 1.0000000000000002e-06 * 6.0640000000000002e+27 
               * exp(-3.3220000000000001 * tc[0] - 0.50321666580471969 * (120790) * invT);
    Corr  = mixture + ( 3 - 1)*sc[1] + ( 0 - 1)*sc[4] + ( 1.1000000000000001 - 1)*sc[10] + ( 2 - 1)*sc[8] + ( 1.5 - 1)*sc[5] + ( 1.8999999999999999 - 1)*sc[11] + ( 3.7999999999999998 - 1)*sc[12];
    qf[5] *= Corr * k_f;
    qr[5] *= Corr * k_f / (exp(-g_RT[0] - g_RT[3] + g_RT[4]) * refC);
    // (13):  H2O + H2O <=> H + OH + H2O
    k_f = 1.0000000000000002e-06 * 1.006e+26 
               * exp(-2.4399999999999999 * tc[0] - 0.50321666580471969 * (120180) * invT);
    Corr  = 1.0;
    qf[15] *= Corr * k_f;
    qr[15] *= Corr * k_f / (exp(-g_RT[0] - g_RT[3] + g_RT[4] + g_RT[4] - g_RT[4]) * refC);
    // (14):  H + O2 (+M) <=> HO2 (+M)
    k_f = 1.0000000000000002e-06 * 4650840000000 
               * exp(0.44 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = mixture + ( 2 - 1)*sc[1] + ( 14 - 1)*sc[4] + ( 0.78000000000000003 - 1)*sc[5] + ( 1.8999999999999999 - 1)*sc[11] + ( 3.7999999999999998 - 1)*sc[12] + ( 0.67000000000000004 - 1)*sc[9] + ( 0.80000000000000004 - 1)*sc[10];
    redP = Corr / k_f * 1e-12 * 6.366e+20 
               * exp(-1.72  * tc[0] - 0.50321666580471969  * (524.79999999999995) *invT);
    F = redP / (1.0 + redP);
    logPred = log10(redP);
    logFcent = log10(
        (1.-0.5)*exp(-tc[1] / 1.0000000000000001e-30) 
        + 0.5 * exp(-tc[1]/1e+30)  
        + 0.);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10., logFcent / (1.0 + troe*troe));
    Corr = F * F_troe;
    qf[0] *= Corr * k_f;
    qr[0] *= Corr * k_f / (exp(g_RT[0] + g_RT[5] - g_RT[6]) * refCinv);
    // (15):  HO2 + H <=> H2 + O2
    k_f = 1.0000000000000002e-06 * 2750000 
               * exp(2.0899999999999999 * tc[0] - 0.50321666580471969 * (-1451) * invT);
    Corr  = 1.0;
    qf[16] *= Corr * k_f;
    qr[16] *= Corr * k_f / exp(g_RT[0] - g_RT[1] - g_RT[5] + g_RT[6]);
    // (16):  HO2 + H <=> OH + OH
    k_f = 1.0000000000000002e-06 * 70790000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (295) * invT);
    Corr  = 1.0;
    qf[17] *= Corr * k_f;
    qr[17] *= Corr * k_f / exp(g_RT[0] - g_RT[3] - g_RT[3] + g_RT[6]);
    // (17):  HO2 + O <=> O2 + OH
    k_f = 1.0000000000000002e-06 * 28500000000 
               * exp(1 * tc[0] - 0.50321666580471969 * (-723.92999999999995) * invT);
    Corr  = 1.0;
    qf[18] *= Corr * k_f;
    qr[18] *= Corr * k_f / exp(g_RT[2] - g_RT[3] - g_RT[5] + g_RT[6]);
    // (18):  HO2 + OH <=> H2O + O2
    k_f = 1.0000000000000002e-06 * 28900000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (-497) * invT);
    Corr  = 1.0;
    qf[19] *= Corr * k_f;
    qr[19] *= Corr * k_f / exp(g_RT[3] - g_RT[4] - g_RT[5] + g_RT[6]);
    // (19):  HO2 + HO2 <=> H2O2 + O2
    k_f = 1.0000000000000002e-06 * 420000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (11982) * invT);
    Corr  = 1.0;
    qf[20] *= Corr * k_f;
    qr[20] *= Corr * k_f / exp(-g_RT[5] + g_RT[6] + g_RT[6] - g_RT[7]);
    // (20):  HO2 + HO2 <=> H2O2 + O2
    k_f = 1.0000000000000002e-06 * 130000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (-1629.3) * invT);
    Corr  = 1.0;
    qf[21] *= Corr * k_f;
    qr[21] *= Corr * k_f / exp(-g_RT[5] + g_RT[6] + g_RT[6] - g_RT[7]);
    // (21):  H2O2 (+M) <=> OH + OH (+M)
    k_f = 1 * 2000000000000 
               * exp(0.90000000000000002 * tc[0] - 0.50321666580471969 * (48749) * invT);
    Corr  = mixture + ( 7.5 - 1)*sc[4] + ( 1.6000000000000001 - 1)*sc[12] + ( 1.5 - 1)*sc[8] + ( 1.2 - 1)*sc[5] + ( 0.65000000000000002 - 1)*sc[10] + ( 7.7000000000000002 - 1)*sc[7] + ( 3.7000000000000002 - 1)*sc[1] + ( 2.7999999999999998 - 1)*sc[11];
    redP = Corr / k_f * 1e-6 * 2.49e+24 
               * exp(-2.2999999999999998  * tc[0] - 0.50321666580471969  * (48749) *invT);
    F = redP / (1.0 + redP);
    logPred = log10(redP);
    logFcent = log10(
        (1.-0.42999999999999999)*exp(-tc[1] / 1.0000000000000001e-30) 
        + 0.42999999999999999 * exp(-tc[1]/1e+30)  
        + 0.);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10., logFcent / (1.0 + troe*troe));
    Corr = F * F_troe;
    qf[1] *= Corr * k_f;
    qr[1] *= Corr * k_f / (exp(-g_RT[3] - g_RT[3] + g_RT[7]) * refC);
    // (22):  H2O2 + H <=> H2O + OH
    k_f = 1.0000000000000002e-06 * 24100000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (3970) * invT);
    Corr  = 1.0;
    qf[22] *= Corr * k_f;
    qr[22] *= Corr * k_f / exp(g_RT[0] - g_RT[3] - g_RT[4] + g_RT[7]);
    // (23):  H2O2 + H <=> HO2 + H2
    k_f = 1.0000000000000002e-06 * 48200000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (7950) * invT);
    Corr  = 1.0;
    qf[23] *= Corr * k_f;
    qr[23] *= Corr * k_f / exp(g_RT[0] - g_RT[1] - g_RT[6] + g_RT[7]);
    // (24):  H2O2 + O <=> OH + HO2
    k_f = 1.0000000000000002e-06 * 9550000 
               * exp(2 * tc[0] - 0.50321666580471969 * (3970) * invT);
    Corr  = 1.0;
    qf[24] *= Corr * k_f;
    qr[24] *= Corr * k_f / exp(g_RT[2] - g_RT[3] - g_RT[6] + g_RT[7]);
    // (25):  H2O2 + OH <=> HO2 + H2O
    k_f = 1.0000000000000002e-06 * 1740000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (318) * invT);
    Corr  = 1.0;
    qf[25] *= Corr * k_f;
    qr[25] *= Corr * k_f / exp(g_RT[3] - g_RT[4] - g_RT[6] + g_RT[7]);
    // (26):  H2O2 + OH <=> HO2 + H2O
    k_f = 1.0000000000000002e-06 * 75900000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (7270) * invT);
    Corr  = 1.0;
    qf[26] *= Corr * k_f;
    qr[26] *= Corr * k_f / exp(g_RT[3] - g_RT[4] - g_RT[6] + g_RT[7]);


    return;
}
#endif


#ifndef AMREX_USE_CUDA
static double T_save = -1;
#ifdef _OPENMP
#pragma omp threadprivate(T_save)
#endif

static double k_f_save[27];
#ifdef _OPENMP
#pragma omp threadprivate(k_f_save)
#endif

static double Kc_save[27];
#ifdef _OPENMP
#pragma omp threadprivate(Kc_save)
#endif


/*compute the production rate for each species pointwise on CPU */
void productionRate(double *  wdot, double *  sc, double T)
{
    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    double invT = 1.0 / tc[1];

    if (T != T_save)
    {
        T_save = T;
        comp_k_f(tc,invT,k_f_save);
        comp_Kc(tc,invT,Kc_save);
    }

    double qdot, q_f[27], q_r[27];
    comp_qfqr(q_f, q_r, sc, tc, invT);

    for (int i = 0; i < 13; ++i) {
        wdot[i] = 0.0;
    }

    qdot = q_f[0]-q_r[0];
    wdot[0] -= qdot;
    wdot[5] -= qdot;
    wdot[6] += qdot;

    qdot = q_f[1]-q_r[1];
    wdot[3] += qdot;
    wdot[3] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[2]-q_r[2];
    wdot[0] += qdot;
    wdot[0] += qdot;
    wdot[1] -= qdot;

    qdot = q_f[3]-q_r[3];
    wdot[2] -= qdot;
    wdot[2] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[4]-q_r[4];
    wdot[0] -= qdot;
    wdot[2] -= qdot;
    wdot[3] += qdot;

    qdot = q_f[5]-q_r[5];
    wdot[0] += qdot;
    wdot[3] += qdot;
    wdot[4] -= qdot;

    qdot = q_f[6]-q_r[6];
    wdot[0] -= qdot;
    wdot[2] += qdot;
    wdot[3] += qdot;
    wdot[5] -= qdot;

    qdot = q_f[7]-q_r[7];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[2] -= qdot;
    wdot[3] += qdot;

    qdot = q_f[8]-q_r[8];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[2] -= qdot;
    wdot[3] += qdot;

    qdot = q_f[9]-q_r[9];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[10]-q_r[10];
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[3] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[11]-q_r[11];
    wdot[0] += qdot;
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[9] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[12]-q_r[12];
    wdot[0] += qdot;
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[10] -= qdot;
    wdot[10] += qdot;

    qdot = q_f[13]-q_r[13];
    wdot[2] -= qdot;
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[9] -= qdot;
    wdot[9] += qdot;

    qdot = q_f[14]-q_r[14];
    wdot[2] -= qdot;
    wdot[2] -= qdot;
    wdot[5] += qdot;
    wdot[10] -= qdot;
    wdot[10] += qdot;

    qdot = q_f[15]-q_r[15];
    wdot[0] += qdot;
    wdot[3] += qdot;
    wdot[4] -= qdot;
    wdot[4] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[16]-q_r[16];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[5] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[17]-q_r[17];
    wdot[0] -= qdot;
    wdot[3] += qdot;
    wdot[3] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[18]-q_r[18];
    wdot[2] -= qdot;
    wdot[3] += qdot;
    wdot[5] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[19]-q_r[19];
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[5] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[20]-q_r[20];
    wdot[5] += qdot;
    wdot[6] -= qdot;
    wdot[6] -= qdot;
    wdot[7] += qdot;

    qdot = q_f[21]-q_r[21];
    wdot[5] += qdot;
    wdot[6] -= qdot;
    wdot[6] -= qdot;
    wdot[7] += qdot;

    qdot = q_f[22]-q_r[22];
    wdot[0] -= qdot;
    wdot[3] += qdot;
    wdot[4] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[23]-q_r[23];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[24]-q_r[24];
    wdot[2] -= qdot;
    wdot[3] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[25]-q_r[25];
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[26]-q_r[26];
    wdot[3] -= qdot;
    wdot[4] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    return;
}

void comp_k_f(double *  tc, double invT, double *  k_f)
{
#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<27; ++i) {
        k_f[i] = prefactor_units[i] * fwd_A[i]
                    * exp(fwd_beta[i] * tc[0] - activation_units[i] * fwd_Ea[i] * invT);
    };
    return;
}

void comp_Kc(double *  tc, double invT, double *  Kc)
{
    /*compute the Gibbs free energy */
    double g_RT[13];
    gibbs(g_RT, tc);

    Kc[0] = g_RT[0] + g_RT[5] - g_RT[6];
    Kc[1] = -g_RT[3] - g_RT[3] + g_RT[7];
    Kc[2] = -g_RT[0] - g_RT[0] + g_RT[1];
    Kc[3] = g_RT[2] + g_RT[2] - g_RT[5];
    Kc[4] = g_RT[0] + g_RT[2] - g_RT[3];
    Kc[5] = -g_RT[0] - g_RT[3] + g_RT[4];
    Kc[6] = g_RT[0] - g_RT[2] - g_RT[3] + g_RT[5];
    Kc[7] = -g_RT[0] + g_RT[1] + g_RT[2] - g_RT[3];
    Kc[8] = -g_RT[0] + g_RT[1] + g_RT[2] - g_RT[3];
    Kc[9] = -g_RT[0] + g_RT[1] + g_RT[3] - g_RT[4];
    Kc[10] = -g_RT[2] + g_RT[3] + g_RT[3] - g_RT[4];
    Kc[11] = -g_RT[0] - g_RT[0] + g_RT[1] + g_RT[9] - g_RT[9];
    Kc[12] = -g_RT[0] - g_RT[0] + g_RT[1] + g_RT[10] - g_RT[10];
    Kc[13] = g_RT[2] + g_RT[2] - g_RT[5] + g_RT[9] - g_RT[9];
    Kc[14] = g_RT[2] + g_RT[2] - g_RT[5] + g_RT[10] - g_RT[10];
    Kc[15] = -g_RT[0] - g_RT[3] + g_RT[4] + g_RT[4] - g_RT[4];
    Kc[16] = g_RT[0] - g_RT[1] - g_RT[5] + g_RT[6];
    Kc[17] = g_RT[0] - g_RT[3] - g_RT[3] + g_RT[6];
    Kc[18] = g_RT[2] - g_RT[3] - g_RT[5] + g_RT[6];
    Kc[19] = g_RT[3] - g_RT[4] - g_RT[5] + g_RT[6];
    Kc[20] = -g_RT[5] + g_RT[6] + g_RT[6] - g_RT[7];
    Kc[21] = -g_RT[5] + g_RT[6] + g_RT[6] - g_RT[7];
    Kc[22] = g_RT[0] - g_RT[3] - g_RT[4] + g_RT[7];
    Kc[23] = g_RT[0] - g_RT[1] - g_RT[6] + g_RT[7];
    Kc[24] = g_RT[2] - g_RT[3] - g_RT[6] + g_RT[7];
    Kc[25] = g_RT[3] - g_RT[4] - g_RT[6] + g_RT[7];
    Kc[26] = g_RT[3] - g_RT[4] - g_RT[6] + g_RT[7];

#ifdef __INTEL_COMPILER
     #pragma simd
#endif
    for (int i=0; i<27; ++i) {
        Kc[i] = exp(Kc[i]);
    };

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31446 * invT;
    double refCinv = 1 / refC;

    Kc[0] *= refCinv;
    Kc[1] *= refC;
    Kc[2] *= refC;
    Kc[3] *= refCinv;
    Kc[4] *= refCinv;
    Kc[5] *= refC;
    Kc[11] *= refC;
    Kc[12] *= refC;
    Kc[13] *= refCinv;
    Kc[14] *= refCinv;
    Kc[15] *= refC;

    return;
}

void comp_qfqr(double *  qf, double *  qr, double *  sc, double *  tc, double invT)
{

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    qf[0] = sc[0]*sc[5];
    qr[0] = sc[6];

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    qf[1] = sc[7];
    qr[1] = sc[3]*sc[3];

    /*reaction 3: H2 + M <=> H + H + M */
    qf[2] = sc[1];
    qr[2] = sc[0]*sc[0];

    /*reaction 4: O + O + M <=> O2 + M */
    qf[3] = sc[2]*sc[2];
    qr[3] = sc[5];

    /*reaction 5: O + H + M <=> OH + M */
    qf[4] = sc[0]*sc[2];
    qr[4] = sc[3];

    /*reaction 6: H2O + M <=> H + OH + M */
    qf[5] = sc[4];
    qr[5] = sc[0]*sc[3];

    /*reaction 7: H + O2 <=> O + OH */
    qf[6] = sc[0]*sc[5];
    qr[6] = sc[2]*sc[3];

    /*reaction 8: O + H2 <=> H + OH */
    qf[7] = sc[1]*sc[2];
    qr[7] = sc[0]*sc[3];

    /*reaction 9: O + H2 <=> H + OH */
    qf[8] = sc[1]*sc[2];
    qr[8] = sc[0]*sc[3];

    /*reaction 10: H2 + OH <=> H2O + H */
    qf[9] = sc[1]*sc[3];
    qr[9] = sc[0]*sc[4];

    /*reaction 11: OH + OH <=> O + H2O */
    qf[10] = sc[3]*sc[3];
    qr[10] = sc[2]*sc[4];

    /*reaction 12: H2 + AR <=> H + H + AR */
    qf[11] = sc[1]*sc[9];
    qr[11] = sc[0]*sc[0]*sc[9];

    /*reaction 13: H2 + HE <=> H + H + HE */
    qf[12] = sc[1]*sc[10];
    qr[12] = sc[0]*sc[0]*sc[10];

    /*reaction 14: O + O + AR <=> O2 + AR */
    qf[13] = sc[2]*sc[2]*sc[9];
    qr[13] = sc[5]*sc[9];

    /*reaction 15: O + O + HE <=> O2 + HE */
    qf[14] = sc[2]*sc[2]*sc[10];
    qr[14] = sc[5]*sc[10];

    /*reaction 16: H2O + H2O <=> H + OH + H2O */
    qf[15] = sc[4]*sc[4];
    qr[15] = sc[0]*sc[3]*sc[4];

    /*reaction 17: HO2 + H <=> H2 + O2 */
    qf[16] = sc[0]*sc[6];
    qr[16] = sc[1]*sc[5];

    /*reaction 18: HO2 + H <=> OH + OH */
    qf[17] = sc[0]*sc[6];
    qr[17] = sc[3]*sc[3];

    /*reaction 19: HO2 + O <=> O2 + OH */
    qf[18] = sc[2]*sc[6];
    qr[18] = sc[3]*sc[5];

    /*reaction 20: HO2 + OH <=> H2O + O2 */
    qf[19] = sc[3]*sc[6];
    qr[19] = sc[4]*sc[5];

    /*reaction 21: HO2 + HO2 <=> H2O2 + O2 */
    qf[20] = sc[6]*sc[6];
    qr[20] = sc[5]*sc[7];

    /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
    qf[21] = sc[6]*sc[6];
    qr[21] = sc[5]*sc[7];

    /*reaction 23: H2O2 + H <=> H2O + OH */
    qf[22] = sc[0]*sc[7];
    qr[22] = sc[3]*sc[4];

    /*reaction 24: H2O2 + H <=> HO2 + H2 */
    qf[23] = sc[0]*sc[7];
    qr[23] = sc[1]*sc[6];

    /*reaction 25: H2O2 + O <=> OH + HO2 */
    qf[24] = sc[2]*sc[7];
    qr[24] = sc[3]*sc[6];

    /*reaction 26: H2O2 + OH <=> HO2 + H2O */
    qf[25] = sc[3]*sc[7];
    qr[25] = sc[4]*sc[6];

    /*reaction 27: H2O2 + OH <=> HO2 + H2O */
    qf[26] = sc[3]*sc[7];
    qr[26] = sc[4]*sc[6];

    double T = tc[1];

    /*compute the mixture concentration */
    double mixture = 0.0;
    for (int i = 0; i < 13; ++i) {
        mixture += sc[i];
    }

    double Corr[27];
    for (int i = 0; i < 27; ++i) {
        Corr[i] = 1.0;
    }

    /* troe */
    {
        double alpha[2];
        alpha[0] = mixture + (TB[0][0] - 1)*sc[1] + (TB[0][1] - 1)*sc[4] + (TB[0][2] - 1)*sc[5] + (TB[0][3] - 1)*sc[11] + (TB[0][4] - 1)*sc[12] + (TB[0][5] - 1)*sc[9] + (TB[0][6] - 1)*sc[10];
        alpha[1] = mixture + (TB[1][0] - 1)*sc[4] + (TB[1][1] - 1)*sc[12] + (TB[1][2] - 1)*sc[8] + (TB[1][3] - 1)*sc[5] + (TB[1][4] - 1)*sc[10] + (TB[1][5] - 1)*sc[7] + (TB[1][6] - 1)*sc[1] + (TB[1][7] - 1)*sc[11];
        for (int i=0; i<2; i++)
        {
            double redP, F, logPred, logFcent, troe_c, troe_n, troe, F_troe;
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
        double alpha;
        alpha = mixture + (TB[2][0] - 1)*sc[1] + (TB[2][1] - 1)*sc[4] + (TB[2][2] - 1)*sc[11] + (TB[2][3] - 1)*sc[12] + (TB[2][4] - 1)*sc[9] + (TB[2][5] - 1)*sc[10];
        Corr[2] = alpha;
        alpha = mixture + (TB[3][0] - 1)*sc[1] + (TB[3][1] - 1)*sc[4] + (TB[3][2] - 1)*sc[9] + (TB[3][3] - 1)*sc[10] + (TB[3][4] - 1)*sc[11] + (TB[3][5] - 1)*sc[12];
        Corr[3] = alpha;
        alpha = mixture + (TB[4][0] - 1)*sc[1] + (TB[4][1] - 1)*sc[4] + (TB[4][2] - 1)*sc[9] + (TB[4][3] - 1)*sc[10] + (TB[4][4] - 1)*sc[11] + (TB[4][5] - 1)*sc[12];
        Corr[4] = alpha;
        alpha = mixture + (TB[5][0] - 1)*sc[1] + (TB[5][1] - 1)*sc[4] + (TB[5][2] - 1)*sc[10] + (TB[5][3] - 1)*sc[8] + (TB[5][4] - 1)*sc[5] + (TB[5][5] - 1)*sc[11] + (TB[5][6] - 1)*sc[12];
        Corr[5] = alpha;
    }

    for (int i=0; i<27; i++)
    {
        qf[i] *= Corr[i] * k_f_save[i];
        qr[i] *= Corr[i] * k_f_save[i] / Kc_save[i];
    }

    return;
}
#endif


#ifndef AMREX_USE_CUDA
/*compute the production rate for each species */
void vproductionRate(int npt, double *  wdot, double *  sc, double *  T)
{
    double k_f_s[27*npt], Kc_s[27*npt], mixture[npt], g_RT[13*npt];
    double tc[5*npt], invT[npt];

#ifdef __INTEL_COMPILER
     #pragma simd
#endif
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

    for (int n=0; n<13; n++) {
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

void vcomp_k_f(int npt, double *  k_f_s, double *  tc, double *  invT)
{
#ifdef __INTEL_COMPILER
    #pragma simd
#endif
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

void vcomp_gibbs(int npt, double *  g_RT, double *  tc)
{
    /*compute the Gibbs free energy */
    for (int i=0; i<npt; i++) {
        double tg[5], g[13];
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
    }
}

void vcomp_Kc(int npt, double *  Kc_s, double *  g_RT, double *  invT)
{
#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<npt; i++) {
        /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
        double refC = (101325. / 8.31451) * invT[i];
        double refCinv = 1.0 / refC;

        Kc_s[0*npt+i] = refCinv * exp((g_RT[0*npt+i] + g_RT[5*npt+i]) - (g_RT[6*npt+i]));
        Kc_s[1*npt+i] = refC * exp((g_RT[7*npt+i]) - (g_RT[3*npt+i] + g_RT[3*npt+i]));
        Kc_s[2*npt+i] = refC * exp((g_RT[1*npt+i]) - (g_RT[0*npt+i] + g_RT[0*npt+i]));
        Kc_s[3*npt+i] = refCinv * exp((g_RT[2*npt+i] + g_RT[2*npt+i]) - (g_RT[5*npt+i]));
        Kc_s[4*npt+i] = refCinv * exp((g_RT[0*npt+i] + g_RT[2*npt+i]) - (g_RT[3*npt+i]));
        Kc_s[5*npt+i] = refC * exp((g_RT[4*npt+i]) - (g_RT[0*npt+i] + g_RT[3*npt+i]));
        Kc_s[6*npt+i] = exp((g_RT[0*npt+i] + g_RT[5*npt+i]) - (g_RT[2*npt+i] + g_RT[3*npt+i]));
        Kc_s[7*npt+i] = exp((g_RT[1*npt+i] + g_RT[2*npt+i]) - (g_RT[0*npt+i] + g_RT[3*npt+i]));
        Kc_s[8*npt+i] = exp((g_RT[1*npt+i] + g_RT[2*npt+i]) - (g_RT[0*npt+i] + g_RT[3*npt+i]));
        Kc_s[9*npt+i] = exp((g_RT[1*npt+i] + g_RT[3*npt+i]) - (g_RT[0*npt+i] + g_RT[4*npt+i]));
        Kc_s[10*npt+i] = exp((g_RT[3*npt+i] + g_RT[3*npt+i]) - (g_RT[2*npt+i] + g_RT[4*npt+i]));
        Kc_s[11*npt+i] = refC * exp((g_RT[1*npt+i] + g_RT[9*npt+i]) - (g_RT[0*npt+i] + g_RT[0*npt+i] + g_RT[9*npt+i]));
        Kc_s[12*npt+i] = refC * exp((g_RT[1*npt+i] + g_RT[10*npt+i]) - (g_RT[0*npt+i] + g_RT[0*npt+i] + g_RT[10*npt+i]));
        Kc_s[13*npt+i] = refCinv * exp((g_RT[2*npt+i] + g_RT[2*npt+i] + g_RT[9*npt+i]) - (g_RT[5*npt+i] + g_RT[9*npt+i]));
        Kc_s[14*npt+i] = refCinv * exp((g_RT[2*npt+i] + g_RT[2*npt+i] + g_RT[10*npt+i]) - (g_RT[5*npt+i] + g_RT[10*npt+i]));
        Kc_s[15*npt+i] = refC * exp((g_RT[4*npt+i] + g_RT[4*npt+i]) - (g_RT[0*npt+i] + g_RT[3*npt+i] + g_RT[4*npt+i]));
        Kc_s[16*npt+i] = exp((g_RT[0*npt+i] + g_RT[6*npt+i]) - (g_RT[1*npt+i] + g_RT[5*npt+i]));
        Kc_s[17*npt+i] = exp((g_RT[0*npt+i] + g_RT[6*npt+i]) - (g_RT[3*npt+i] + g_RT[3*npt+i]));
        Kc_s[18*npt+i] = exp((g_RT[2*npt+i] + g_RT[6*npt+i]) - (g_RT[3*npt+i] + g_RT[5*npt+i]));
        Kc_s[19*npt+i] = exp((g_RT[3*npt+i] + g_RT[6*npt+i]) - (g_RT[4*npt+i] + g_RT[5*npt+i]));
        Kc_s[20*npt+i] = exp((g_RT[6*npt+i] + g_RT[6*npt+i]) - (g_RT[5*npt+i] + g_RT[7*npt+i]));
        Kc_s[21*npt+i] = exp((g_RT[6*npt+i] + g_RT[6*npt+i]) - (g_RT[5*npt+i] + g_RT[7*npt+i]));
        Kc_s[22*npt+i] = exp((g_RT[0*npt+i] + g_RT[7*npt+i]) - (g_RT[3*npt+i] + g_RT[4*npt+i]));
        Kc_s[23*npt+i] = exp((g_RT[0*npt+i] + g_RT[7*npt+i]) - (g_RT[1*npt+i] + g_RT[6*npt+i]));
        Kc_s[24*npt+i] = exp((g_RT[2*npt+i] + g_RT[7*npt+i]) - (g_RT[3*npt+i] + g_RT[6*npt+i]));
        Kc_s[25*npt+i] = exp((g_RT[3*npt+i] + g_RT[7*npt+i]) - (g_RT[4*npt+i] + g_RT[6*npt+i]));
        Kc_s[26*npt+i] = exp((g_RT[3*npt+i] + g_RT[7*npt+i]) - (g_RT[4*npt+i] + g_RT[6*npt+i]));
    }
}

void vcomp_wdot(int npt, double *  wdot, double *  mixture, double *  sc,
		double *  k_f_s, double *  Kc_s,
		double *  tc, double *  invT, double *  T)
{
#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<npt; i++) {
        double qdot, q_f, q_r, phi_f, phi_r, k_f, k_r, Kc;
        double alpha;
        double redP, F;
        double logPred;
        double logFcent, troe_c, troe_n, troe, F_troe;

        /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
        phi_f = sc[0*npt+i]*sc[5*npt+i];
        alpha = mixture[i] + (TB[0][0] - 1)*sc[1*npt+i] + (TB[0][1] - 1)*sc[4*npt+i] + (TB[0][2] - 1)*sc[5*npt+i] + (TB[0][3] - 1)*sc[11*npt+i] + (TB[0][4] - 1)*sc[12*npt+i] + (TB[0][5] - 1)*sc[9*npt+i] + (TB[0][6] - 1)*sc[10*npt+i];
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
        wdot[0*npt+i] -= qdot;
        wdot[5*npt+i] -= qdot;
        wdot[6*npt+i] += qdot;

        /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
        phi_f = sc[7*npt+i];
        alpha = mixture[i] + (TB[1][0] - 1)*sc[4*npt+i] + (TB[1][1] - 1)*sc[12*npt+i] + (TB[1][2] - 1)*sc[8*npt+i] + (TB[1][3] - 1)*sc[5*npt+i] + (TB[1][4] - 1)*sc[10*npt+i] + (TB[1][5] - 1)*sc[7*npt+i] + (TB[1][6] - 1)*sc[1*npt+i] + (TB[1][7] - 1)*sc[11*npt+i];
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

        /*reaction 3: H2 + M <=> H + H + M */
        phi_f = sc[1*npt+i];
        alpha = mixture[i] + (TB[2][0] - 1)*sc[1*npt+i] + (TB[2][1] - 1)*sc[4*npt+i] + (TB[2][2] - 1)*sc[11*npt+i] + (TB[2][3] - 1)*sc[12*npt+i] + (TB[2][4] - 1)*sc[9*npt+i] + (TB[2][5] - 1)*sc[10*npt+i];
        k_f = alpha * k_f_s[2*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[0*npt+i];
        Kc = Kc_s[2*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;

        /*reaction 4: O + O + M <=> O2 + M */
        phi_f = sc[2*npt+i]*sc[2*npt+i];
        alpha = mixture[i] + (TB[3][0] - 1)*sc[1*npt+i] + (TB[3][1] - 1)*sc[4*npt+i] + (TB[3][2] - 1)*sc[9*npt+i] + (TB[3][3] - 1)*sc[10*npt+i] + (TB[3][4] - 1)*sc[11*npt+i] + (TB[3][5] - 1)*sc[12*npt+i];
        k_f = alpha * k_f_s[3*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i];
        Kc = Kc_s[3*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[2*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;

        /*reaction 5: O + H + M <=> OH + M */
        phi_f = sc[0*npt+i]*sc[2*npt+i];
        alpha = mixture[i] + (TB[4][0] - 1)*sc[1*npt+i] + (TB[4][1] - 1)*sc[4*npt+i] + (TB[4][2] - 1)*sc[9*npt+i] + (TB[4][3] - 1)*sc[10*npt+i] + (TB[4][4] - 1)*sc[11*npt+i] + (TB[4][5] - 1)*sc[12*npt+i];
        k_f = alpha * k_f_s[4*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i];
        Kc = Kc_s[4*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[2*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;

        /*reaction 6: H2O + M <=> H + OH + M */
        phi_f = sc[4*npt+i];
        alpha = mixture[i] + (TB[5][0] - 1)*sc[1*npt+i] + (TB[5][1] - 1)*sc[4*npt+i] + (TB[5][2] - 1)*sc[10*npt+i] + (TB[5][3] - 1)*sc[8*npt+i] + (TB[5][4] - 1)*sc[5*npt+i] + (TB[5][5] - 1)*sc[11*npt+i] + (TB[5][6] - 1)*sc[12*npt+i];
        k_f = alpha * k_f_s[5*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[3*npt+i];
        Kc = Kc_s[5*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[3*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;

        /*reaction 7: H + O2 <=> O + OH */
        phi_f = sc[0*npt+i]*sc[5*npt+i];
        k_f = k_f_s[6*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[3*npt+i];
        Kc = Kc_s[6*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] += qdot;
        wdot[5*npt+i] -= qdot;

        /*reaction 8: O + H2 <=> H + OH */
        phi_f = sc[1*npt+i]*sc[2*npt+i];
        k_f = k_f_s[7*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[3*npt+i];
        Kc = Kc_s[7*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;

        /*reaction 9: O + H2 <=> H + OH */
        phi_f = sc[1*npt+i]*sc[2*npt+i];
        k_f = k_f_s[8*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[3*npt+i];
        Kc = Kc_s[8*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[2*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;

        /*reaction 10: H2 + OH <=> H2O + H */
        phi_f = sc[1*npt+i]*sc[3*npt+i];
        k_f = k_f_s[9*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[4*npt+i];
        Kc = Kc_s[9*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;

        /*reaction 11: OH + OH <=> O + H2O */
        phi_f = sc[3*npt+i]*sc[3*npt+i];
        k_f = k_f_s[10*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[2*npt+i]*sc[4*npt+i];
        Kc = Kc_s[10*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] += qdot;
        wdot[3*npt+i] -= qdot;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;

        /*reaction 12: H2 + AR <=> H + H + AR */
        phi_f = sc[1*npt+i]*sc[9*npt+i];
        k_f = k_f_s[11*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[0*npt+i]*sc[9*npt+i];
        Kc = Kc_s[11*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[9*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;

        /*reaction 13: H2 + HE <=> H + H + HE */
        phi_f = sc[1*npt+i]*sc[10*npt+i];
        k_f = k_f_s[12*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[0*npt+i]*sc[10*npt+i];
        Kc = Kc_s[12*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[0*npt+i] += qdot;
        wdot[1*npt+i] -= qdot;
        wdot[10*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;

        /*reaction 14: O + O + AR <=> O2 + AR */
        phi_f = sc[2*npt+i]*sc[2*npt+i]*sc[9*npt+i];
        k_f = k_f_s[13*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[9*npt+i];
        Kc = Kc_s[13*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[2*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[9*npt+i] -= qdot;
        wdot[9*npt+i] += qdot;

        /*reaction 15: O + O + HE <=> O2 + HE */
        phi_f = sc[2*npt+i]*sc[2*npt+i]*sc[10*npt+i];
        k_f = k_f_s[14*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[10*npt+i];
        Kc = Kc_s[14*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[2*npt+i] -= qdot;
        wdot[5*npt+i] += qdot;
        wdot[10*npt+i] -= qdot;
        wdot[10*npt+i] += qdot;

        /*reaction 16: H2O + H2O <=> H + OH + H2O */
        phi_f = sc[4*npt+i]*sc[4*npt+i];
        k_f = k_f_s[15*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[0*npt+i]*sc[3*npt+i]*sc[4*npt+i];
        Kc = Kc_s[15*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] += qdot;
        wdot[3*npt+i] += qdot;
        wdot[4*npt+i] -= qdot;
        wdot[4*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;

        /*reaction 17: HO2 + H <=> H2 + O2 */
        phi_f = sc[0*npt+i]*sc[6*npt+i];
        k_f = k_f_s[16*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[5*npt+i];
        Kc = Kc_s[16*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 18: HO2 + H <=> OH + OH */
        phi_f = sc[0*npt+i]*sc[6*npt+i];
        k_f = k_f_s[17*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[3*npt+i];
        Kc = Kc_s[17*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[3*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 19: HO2 + O <=> O2 + OH */
        phi_f = sc[2*npt+i]*sc[6*npt+i];
        k_f = k_f_s[18*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[5*npt+i];
        Kc = Kc_s[18*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 20: HO2 + OH <=> H2O + O2 */
        phi_f = sc[3*npt+i]*sc[6*npt+i];
        k_f = k_f_s[19*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[5*npt+i];
        Kc = Kc_s[19*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;

        /*reaction 21: HO2 + HO2 <=> H2O2 + O2 */
        phi_f = sc[6*npt+i]*sc[6*npt+i];
        k_f = k_f_s[20*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[7*npt+i];
        Kc = Kc_s[20*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[6*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;

        /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
        phi_f = sc[6*npt+i]*sc[6*npt+i];
        k_f = k_f_s[21*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[5*npt+i]*sc[7*npt+i];
        Kc = Kc_s[21*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[5*npt+i] += qdot;
        wdot[6*npt+i] -= qdot;
        wdot[6*npt+i] -= qdot;
        wdot[7*npt+i] += qdot;

        /*reaction 23: H2O2 + H <=> H2O + OH */
        phi_f = sc[0*npt+i]*sc[7*npt+i];
        k_f = k_f_s[22*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[4*npt+i];
        Kc = Kc_s[22*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[4*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 24: H2O2 + H <=> HO2 + H2 */
        phi_f = sc[0*npt+i]*sc[7*npt+i];
        k_f = k_f_s[23*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[1*npt+i]*sc[6*npt+i];
        Kc = Kc_s[23*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[0*npt+i] -= qdot;
        wdot[1*npt+i] += qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 25: H2O2 + O <=> OH + HO2 */
        phi_f = sc[2*npt+i]*sc[7*npt+i];
        k_f = k_f_s[24*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[3*npt+i]*sc[6*npt+i];
        Kc = Kc_s[24*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[2*npt+i] -= qdot;
        wdot[3*npt+i] += qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 26: H2O2 + OH <=> HO2 + H2O */
        phi_f = sc[3*npt+i]*sc[7*npt+i];
        k_f = k_f_s[25*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[6*npt+i];
        Kc = Kc_s[25*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;

        /*reaction 27: H2O2 + OH <=> HO2 + H2O */
        phi_f = sc[3*npt+i]*sc[7*npt+i];
        k_f = k_f_s[26*npt+i];
        q_f = phi_f * k_f;
        phi_r = sc[4*npt+i]*sc[6*npt+i];
        Kc = Kc_s[26*npt+i];
        k_r = k_f / Kc;
        q_r = phi_r * k_r;
        qdot = q_f - q_r;
        wdot[3*npt+i] -= qdot;
        wdot[4*npt+i] += qdot;
        wdot[6*npt+i] += qdot;
        wdot[7*npt+i] -= qdot;
    }
}
#endif

/*compute an approx to the reaction Jacobian (for preconditioning) */
AMREX_GPU_HOST_DEVICE void DWDOT_SIMPLIFIED(double *  J, double *  sc, double *  Tp, int * HP)
{
    double c[13];

    for (int k=0; k<13; k++) {
        c[k] = 1.e6 * sc[k];
    }

    aJacobian_precond(J, c, *Tp, *HP);

    /* dwdot[k]/dT */
    /* dTdot/d[X] */
    for (int k=0; k<13; k++) {
        J[182+k] *= 1.e-6;
        J[k*14+13] *= 1.e6;
    }

    return;
}

/*compute the reaction Jacobian */
AMREX_GPU_HOST_DEVICE void DWDOT(double *  J, double *  sc, double *  Tp, int * consP)
{
    double c[13];

    for (int k=0; k<13; k++) {
        c[k] = 1.e6 * sc[k];
    }

    aJacobian(J, c, *Tp, *consP);

    /* dwdot[k]/dT */
    /* dTdot/d[X] */
    for (int k=0; k<13; k++) {
        J[182+k] *= 1.e-6;
        J[k*14+13] *= 1.e6;
    }

    return;
}

/*compute the sparsity pattern of the chemistry Jacobian */
AMREX_GPU_HOST_DEVICE void SPARSITY_INFO( int * nJdata, int * consP, int NCELLS)
{
    double c[13];
    double J[196];

    for (int k=0; k<13; k++) {
        c[k] = 1.0/ 13.000000 ;
    }

    aJacobian(J, c, 1500.0, *consP);

    int nJdata_tmp = 0;
    for (int k=0; k<14; k++) {
        for (int l=0; l<14; l++) {
            if(J[ 14 * k + l] != 0.0){
                nJdata_tmp = nJdata_tmp + 1;
            }
        }
    }

    *nJdata = NCELLS * nJdata_tmp;

    return;
}



/*compute the sparsity pattern of the system Jacobian */
AMREX_GPU_HOST_DEVICE void SPARSITY_INFO_SYST( int * nJdata, int * consP, int NCELLS)
{
    double c[13];
    double J[196];

    for (int k=0; k<13; k++) {
        c[k] = 1.0/ 13.000000 ;
    }

    aJacobian(J, c, 1500.0, *consP);

    int nJdata_tmp = 0;
    for (int k=0; k<14; k++) {
        for (int l=0; l<14; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J[ 14 * k + l] != 0.0){
                    nJdata_tmp = nJdata_tmp + 1;
                }
            }
        }
    }

    *nJdata = NCELLS * nJdata_tmp;

    return;
}



/*compute the sparsity pattern of the simplified (for preconditioning) system Jacobian */
AMREX_GPU_HOST_DEVICE void SPARSITY_INFO_SYST_SIMPLIFIED( int * nJdata, int * consP)
{
    double c[13];
    double J[196];

    for (int k=0; k<13; k++) {
        c[k] = 1.0/ 13.000000 ;
    }

    aJacobian_precond(J, c, 1500.0, *consP);

    int nJdata_tmp = 0;
    for (int k=0; k<14; k++) {
        for (int l=0; l<14; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J[ 14 * k + l] != 0.0){
                    nJdata_tmp = nJdata_tmp + 1;
                }
            }
        }
    }

    nJdata[0] = nJdata_tmp;

    return;
}


/*compute the sparsity pattern of the chemistry Jacobian in CSC format -- base 0 */
AMREX_GPU_HOST_DEVICE void SPARSITY_PREPROC_CSC(int *  rowVals, int *  colPtrs, int * consP, int NCELLS)
{
    double c[13];
    double J[196];
    int offset_row;
    int offset_col;

    for (int k=0; k<13; k++) {
        c[k] = 1.0/ 13.000000 ;
    }

    aJacobian(J, c, 1500.0, *consP);

    colPtrs[0] = 0;
    int nJdata_tmp = 0;
    for (int nc=0; nc<NCELLS; nc++) {
        offset_row = nc * 14;
        offset_col = nc * 14;
        for (int k=0; k<14; k++) {
            for (int l=0; l<14; l++) {
                if(J[14*k + l] != 0.0) {
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
AMREX_GPU_HOST_DEVICE void SPARSITY_PREPROC_CSR(int * colVals, int * rowPtrs, int * consP, int NCELLS, int base)
{
    double c[13];
    double J[196];
    int offset;

    for (int k=0; k<13; k++) {
        c[k] = 1.0/ 13.000000 ;
    }

    aJacobian(J, c, 1500.0, *consP);

    if (base == 1) {
        rowPtrs[0] = 1;
        int nJdata_tmp = 1;
        for (int nc=0; nc<NCELLS; nc++) {
            offset = nc * 14;
            for (int l=0; l<14; l++) {
                for (int k=0; k<14; k++) {
                    if(J[14*k + l] != 0.0) {
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
            offset = nc * 14;
            for (int l=0; l<14; l++) {
                for (int k=0; k<14; k++) {
                    if(J[14*k + l] != 0.0) {
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
AMREX_GPU_HOST_DEVICE void SPARSITY_PREPROC_SYST_CSR(int * colVals, int * rowPtr, int * consP, int NCELLS, int base)
{
    double c[13];
    double J[196];
    int offset;

    for (int k=0; k<13; k++) {
        c[k] = 1.0/ 13.000000 ;
    }

    aJacobian(J, c, 1500.0, *consP);

    if (base == 1) {
        rowPtr[0] = 1;
        int nJdata_tmp = 1;
        for (int nc=0; nc<NCELLS; nc++) {
            offset = nc * 14;
            for (int l=0; l<14; l++) {
                for (int k=0; k<14; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp-1] = l+1 + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(J[14*k + l] != 0.0) {
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
            offset = nc * 14;
            for (int l=0; l<14; l++) {
                for (int k=0; k<14; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp] = l + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(J[14*k + l] != 0.0) {
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
AMREX_GPU_HOST_DEVICE void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(int * rowVals, int * colPtrs, int * indx, int * consP)
{
    double c[13];
    double J[196];

    for (int k=0; k<13; k++) {
        c[k] = 1.0/ 13.000000 ;
    }

    aJacobian_precond(J, c, 1500.0, *consP);

    colPtrs[0] = 0;
    int nJdata_tmp = 0;
    for (int k=0; k<14; k++) {
        for (int l=0; l<14; l++) {
            if (k == l) {
                rowVals[nJdata_tmp] = l; 
                indx[nJdata_tmp] = 14*k + l;
                nJdata_tmp = nJdata_tmp + 1; 
            } else {
                if(J[14*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l; 
                    indx[nJdata_tmp] = 14*k + l;
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
AMREX_GPU_HOST_DEVICE void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(int * colVals, int * rowPtr, int * consP, int base)
{
    double c[13];
    double J[196];

    for (int k=0; k<13; k++) {
        c[k] = 1.0/ 13.000000 ;
    }

    aJacobian_precond(J, c, 1500.0, *consP);

    if (base == 1) {
        rowPtr[0] = 1;
        int nJdata_tmp = 1;
        for (int l=0; l<14; l++) {
            for (int k=0; k<14; k++) {
                if (k == l) {
                    colVals[nJdata_tmp-1] = l+1; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J[14*k + l] != 0.0) {
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
        for (int l=0; l<14; l++) {
            for (int k=0; k<14; k++) {
                if (k == l) {
                    colVals[nJdata_tmp] = l; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J[14*k + l] != 0.0) {
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


#ifdef AMREX_USE_CUDA
/*compute the reaction Jacobian on GPU */
AMREX_GPU_HOST_DEVICE
void aJacobian(double * J, double * sc, double T, int consP)
{


    for (int i=0; i<196; i++) {
        J[i] = 0.0;
    }

    double wdot[13];
    for (int k=0; k<13; k++) {
        wdot[k] = 0.0;
    }

    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    double invT = 1.0 / tc[1];
    double invT2 = invT * invT;

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31446 / T;
    double refCinv = 1.0 / refC;

    /*compute the mixture concentration */
    double mixture = 0.0;
    for (int k = 0; k < 13; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    double g_RT[13];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    double h_RT[13];
    speciesEnthalpy(h_RT, tc);

    double phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    double dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    double dqdci, dcdc_fac, dqdc[13];
    double Pr, fPr, F, k_0, logPr;
    double logFcent, troe_c, troe_n, troePr_den, troePr, troe;
    double Fcent1, Fcent2, Fcent3, Fcent;
    double dlogFdc, dlogFdn, dlogFdcn_fac;
    double dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
    const double ln10 = log(10.0);
    const double log10e = 1.0/log(10.0);
    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[1] + ( 14 - 1)*sc[4] + ( 0.78000000000000003 - 1)*sc[5] + ( 1.8999999999999999 - 1)*sc[11] + ( 3.7999999999999998 - 1)*sc[12] + ( 0.67000000000000004 - 1)*sc[9] + ( 0.80000000000000004 - 1)*sc[10];
    /* forward */
    phi_f = sc[0]*sc[5];
    k_f = 1.0000000000000002e-06 * 4650840000000
                * exp(0.44 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0.44 * invT + 0.50321666580471969 *  0  * invT2;
    /* pressure-fall-off */
    k_0 = 6.366e+20 * exp(-1.72 * tc[0] - 0.50321666580471969 * (524.79999999999995) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -1.72 * invT + 0.50321666580471969 * (524.79999999999995) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-0.5)*exp(-T/1.0000000000000001e-30);
    Fcent2 = 0.5 * exp(-T/1e+30);
    Fcent3 = 0.;
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/1.0000000000000001e-30
        -Fcent2/1e+30
        + 0.);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[6];
    Kc = refCinv * exp(g_RT[0] + g_RT[5] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[5]) + (h_RT[6]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[5] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[5];
        J[0] -= dqdci;                /* dwdot[H]/d[H] */
        J[5] -= dqdci;                /* dwdot[O2]/d[H] */
        J[6] += dqdci;                /* dwdot[HO2]/d[H] */
        /* d()/d[H2] */
        dqdci = (2 - 1)*dcdc_fac;
        J[14] -= dqdci;               /* dwdot[H]/d[H2] */
        J[19] -= dqdci;               /* dwdot[O2]/d[H2] */
        J[20] += dqdci;               /* dwdot[HO2]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (14 - 1)*dcdc_fac;
        J[56] -= dqdci;               /* dwdot[H]/d[H2O] */
        J[61] -= dqdci;               /* dwdot[O2]/d[H2O] */
        J[62] += dqdci;               /* dwdot[HO2]/d[H2O] */
        /* d()/d[O2] */
        dqdci = (0.78000000000000003 - 1)*dcdc_fac + k_f*sc[0];
        J[70] -= dqdci;               /* dwdot[H]/d[O2] */
        J[75] -= dqdci;               /* dwdot[O2]/d[O2] */
        J[76] += dqdci;               /* dwdot[HO2]/d[O2] */
        /* d()/d[HO2] */
        dqdci =  - k_r;
        J[84] -= dqdci;               /* dwdot[H]/d[HO2] */
        J[89] -= dqdci;               /* dwdot[O2]/d[HO2] */
        J[90] += dqdci;               /* dwdot[HO2]/d[HO2] */
        /* d()/d[AR] */
        dqdci = (0.67000000000000004 - 1)*dcdc_fac;
        J[126] -= dqdci;              /* dwdot[H]/d[AR] */
        J[131] -= dqdci;              /* dwdot[O2]/d[AR] */
        J[132] += dqdci;              /* dwdot[HO2]/d[AR] */
        /* d()/d[HE] */
        dqdci = (0.80000000000000004 - 1)*dcdc_fac;
        J[140] -= dqdci;              /* dwdot[H]/d[HE] */
        J[145] -= dqdci;              /* dwdot[O2]/d[HE] */
        J[146] += dqdci;              /* dwdot[HO2]/d[HE] */
        /* d()/d[CO] */
        dqdci = (1.8999999999999999 - 1)*dcdc_fac;
        J[154] -= dqdci;              /* dwdot[H]/d[CO] */
        J[159] -= dqdci;              /* dwdot[O2]/d[CO] */
        J[160] += dqdci;              /* dwdot[HO2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (3.7999999999999998 - 1)*dcdc_fac;
        J[168] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[173] -= dqdci;              /* dwdot[O2]/d[CO2] */
        J[174] += dqdci;              /* dwdot[HO2]/d[CO2] */
    }
    else {
        dqdc[0] = dcdc_fac + k_f*sc[5];
        dqdc[1] = 2*dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = 14*dcdc_fac;
        dqdc[5] = 0.78000000000000003*dcdc_fac + k_f*sc[0];
        dqdc[6] = dcdc_fac - k_r;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = 0.67000000000000004*dcdc_fac;
        dqdc[10] = 0.80000000000000004*dcdc_fac;
        dqdc[11] = 1.8999999999999999*dcdc_fac;
        dqdc[12] = 3.7999999999999998*dcdc_fac;
        for (int k=0; k<13; k++) {
            J[14*k+0] -= dqdc[k];
            J[14*k+5] -= dqdc[k];
            J[14*k+6] += dqdc[k];
        }
    }
    J[182] -= dqdT; /* dwdot[H]/dT */
    J[187] -= dqdT; /* dwdot[O2]/dT */
    J[188] += dqdT; /* dwdot[HO2]/dT */

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 7.5 - 1)*sc[4] + ( 1.6000000000000001 - 1)*sc[12] + ( 1.5 - 1)*sc[8] + ( 1.2 - 1)*sc[5] + ( 0.65000000000000002 - 1)*sc[10] + ( 7.7000000000000002 - 1)*sc[7] + ( 3.7000000000000002 - 1)*sc[1] + ( 2.7999999999999998 - 1)*sc[11];
    /* forward */
    phi_f = sc[7];
    k_f = 1 * 2000000000000
                * exp(0.90000000000000002 * tc[0] - 0.50321666580471969 * (48749) * invT);
    dlnkfdT = 0.90000000000000002 * invT + 0.50321666580471969 *  48749  * invT2;
    /* pressure-fall-off */
    k_0 = 2.49e+24 * exp(-2.2999999999999998 * tc[0] - 0.50321666580471969 * (48749) * invT);
    Pr = 1e-6 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -2.2999999999999998 * invT + 0.50321666580471969 * (48749) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-0.42999999999999999)*exp(-T/1.0000000000000001e-30);
    Fcent2 = 0.42999999999999999 * exp(-T/1e+30);
    Fcent3 = 0.;
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/1.0000000000000001e-30
        -Fcent2/1e+30
        + 0.);
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
        dqdci = (3.7000000000000002 - 1)*dcdc_fac;
        J[17] += 2 * dqdci;           /* dwdot[OH]/d[H2] */
        J[21] -= dqdci;               /* dwdot[H2O2]/d[H2] */
        /* d()/d[OH] */
        dqdci =  - k_r*2.000000*sc[3];
        J[45] += 2 * dqdci;           /* dwdot[OH]/d[OH] */
        J[49] -= dqdci;               /* dwdot[H2O2]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (7.5 - 1)*dcdc_fac;
        J[59] += 2 * dqdci;           /* dwdot[OH]/d[H2O] */
        J[63] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
        /* d()/d[O2] */
        dqdci = (1.2 - 1)*dcdc_fac;
        J[73] += 2 * dqdci;           /* dwdot[OH]/d[O2] */
        J[77] -= dqdci;               /* dwdot[H2O2]/d[O2] */
        /* d()/d[H2O2] */
        dqdci = (7.7000000000000002 - 1)*dcdc_fac + k_f;
        J[101] += 2 * dqdci;          /* dwdot[OH]/d[H2O2] */
        J[105] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
        /* d()/d[N2] */
        dqdci = (1.5 - 1)*dcdc_fac;
        J[115] += 2 * dqdci;          /* dwdot[OH]/d[N2] */
        J[119] -= dqdci;              /* dwdot[H2O2]/d[N2] */
        /* d()/d[HE] */
        dqdci = (0.65000000000000002 - 1)*dcdc_fac;
        J[143] += 2 * dqdci;          /* dwdot[OH]/d[HE] */
        J[147] -= dqdci;              /* dwdot[H2O2]/d[HE] */
        /* d()/d[CO] */
        dqdci = (2.7999999999999998 - 1)*dcdc_fac;
        J[157] += 2 * dqdci;          /* dwdot[OH]/d[CO] */
        J[161] -= dqdci;              /* dwdot[H2O2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (1.6000000000000001 - 1)*dcdc_fac;
        J[171] += 2 * dqdci;          /* dwdot[OH]/d[CO2] */
        J[175] -= dqdci;              /* dwdot[H2O2]/d[CO2] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = 3.7000000000000002*dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac - k_r*2.000000*sc[3];
        dqdc[4] = 7.5*dcdc_fac;
        dqdc[5] = 1.2*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = 7.7000000000000002*dcdc_fac + k_f;
        dqdc[8] = 1.5*dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = 0.65000000000000002*dcdc_fac;
        dqdc[11] = 2.7999999999999998*dcdc_fac;
        dqdc[12] = 1.6000000000000001*dcdc_fac;
        for (int k=0; k<13; k++) {
            J[14*k+3] += 2 * dqdc[k];
            J[14*k+7] -= dqdc[k];
        }
    }
    J[185] += 2 * dqdT; /* dwdot[OH]/dT */
    J[189] -= dqdT; /* dwdot[H2O2]/dT */

    /*reaction 3: H2 + M <=> H + H + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 2.5 - 1)*sc[1] + ( 12 - 1)*sc[4] + ( 1.8999999999999999 - 1)*sc[11] + ( 3.7999999999999998 - 1)*sc[12] + ( 0 - 1)*sc[9] + ( 0 - 1)*sc[10];
    /* forward */
    phi_f = sc[1];
    k_f = 1.0000000000000002e-06 * 4.577e+19
                * exp(-1.3999999999999999 * tc[0] - 0.50321666580471969 * (104380) * invT);
    dlnkfdT = -1.3999999999999999 * invT + 0.50321666580471969 *  104380  * invT2;
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
        dqdci = (2.5 - 1)*q_nocor + k_f;
        J[14] += 2 * dqdci;           /* dwdot[H]/d[H2] */
        J[15] -= dqdci;               /* dwdot[H2]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (12 - 1)*q_nocor;
        J[56] += 2 * dqdci;           /* dwdot[H]/d[H2O] */
        J[57] -= dqdci;               /* dwdot[H2]/d[H2O] */
        /* d()/d[AR] */
        dqdci = (0 - 1)*q_nocor;
        J[126] += 2 * dqdci;          /* dwdot[H]/d[AR] */
        J[127] -= dqdci;              /* dwdot[H2]/d[AR] */
        /* d()/d[HE] */
        dqdci = (0 - 1)*q_nocor;
        J[140] += 2 * dqdci;          /* dwdot[H]/d[HE] */
        J[141] -= dqdci;              /* dwdot[H2]/d[HE] */
        /* d()/d[CO] */
        dqdci = (1.8999999999999999 - 1)*q_nocor;
        J[154] += 2 * dqdci;          /* dwdot[H]/d[CO] */
        J[155] -= dqdci;              /* dwdot[H2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (3.7999999999999998 - 1)*q_nocor;
        J[168] += 2 * dqdci;          /* dwdot[H]/d[CO2] */
        J[169] -= dqdci;              /* dwdot[H2]/d[CO2] */
    }
    else {
        dqdc[0] = q_nocor - k_r*2.000000*sc[0];
        dqdc[1] = 2.5*q_nocor + k_f;
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = 12*q_nocor;
        dqdc[5] = q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[11] = 1.8999999999999999*q_nocor;
        dqdc[12] = 3.7999999999999998*q_nocor;
        for (int k=0; k<13; k++) {
            J[14*k+0] += 2 * dqdc[k];
            J[14*k+1] -= dqdc[k];
        }
    }
    J[182] += 2 * dqdT; /* dwdot[H]/dT */
    J[183] -= dqdT; /* dwdot[H2]/dT */

    /*reaction 4: O + O + M <=> O2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 2.5 - 1)*sc[1] + ( 12 - 1)*sc[4] + ( 0 - 1)*sc[9] + ( 0 - 1)*sc[10] + ( 1.8999999999999999 - 1)*sc[11] + ( 3.7999999999999998 - 1)*sc[12];
    /* forward */
    phi_f = pow(sc[2], 2.000000);
    k_f = 1.0000000000000002e-12 * 6165000000000000
                * exp(-0.5 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -0.5 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[5];
    Kc = refCinv * exp(g_RT[2] + g_RT[2] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[2]) + (h_RT[5]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= 2 * q; /* O */
    wdot[5] += q; /* O2 */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (2.5 - 1)*q_nocor;
        J[16] += -2 * dqdci;          /* dwdot[O]/d[H2] */
        J[19] += dqdci;               /* dwdot[O2]/d[H2] */
        /* d()/d[O] */
        dqdci =  + k_f*2.000000*sc[2];
        J[30] += -2 * dqdci;          /* dwdot[O]/d[O] */
        J[33] += dqdci;               /* dwdot[O2]/d[O] */
        /* d()/d[H2O] */
        dqdci = (12 - 1)*q_nocor;
        J[58] += -2 * dqdci;          /* dwdot[O]/d[H2O] */
        J[61] += dqdci;               /* dwdot[O2]/d[H2O] */
        /* d()/d[O2] */
        dqdci =  - k_r;
        J[72] += -2 * dqdci;          /* dwdot[O]/d[O2] */
        J[75] += dqdci;               /* dwdot[O2]/d[O2] */
        /* d()/d[AR] */
        dqdci = (0 - 1)*q_nocor;
        J[128] += -2 * dqdci;         /* dwdot[O]/d[AR] */
        J[131] += dqdci;              /* dwdot[O2]/d[AR] */
        /* d()/d[HE] */
        dqdci = (0 - 1)*q_nocor;
        J[142] += -2 * dqdci;         /* dwdot[O]/d[HE] */
        J[145] += dqdci;              /* dwdot[O2]/d[HE] */
        /* d()/d[CO] */
        dqdci = (1.8999999999999999 - 1)*q_nocor;
        J[156] += -2 * dqdci;         /* dwdot[O]/d[CO] */
        J[159] += dqdci;              /* dwdot[O2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (3.7999999999999998 - 1)*q_nocor;
        J[170] += -2 * dqdci;         /* dwdot[O]/d[CO2] */
        J[173] += dqdci;              /* dwdot[O2]/d[CO2] */
    }
    else {
        dqdc[0] = q_nocor;
        dqdc[1] = 2.5*q_nocor;
        dqdc[2] = q_nocor + k_f*2.000000*sc[2];
        dqdc[3] = q_nocor;
        dqdc[4] = 12*q_nocor;
        dqdc[5] = q_nocor - k_r;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[11] = 1.8999999999999999*q_nocor;
        dqdc[12] = 3.7999999999999998*q_nocor;
        for (int k=0; k<13; k++) {
            J[14*k+2] += -2 * dqdc[k];
            J[14*k+5] += dqdc[k];
        }
    }
    J[184] += -2 * dqdT; /* dwdot[O]/dT */
    J[187] += dqdT; /* dwdot[O2]/dT */

    /*reaction 5: O + H + M <=> OH + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 2.5 - 1)*sc[1] + ( 12 - 1)*sc[4] + ( 0.75 - 1)*sc[9] + ( 0.75 - 1)*sc[10] + ( 1.8999999999999999 - 1)*sc[11] + ( 3.7999999999999998 - 1)*sc[12];
    /* forward */
    phi_f = sc[0]*sc[2];
    k_f = 1.0000000000000002e-12 * 4.714e+18
                * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[3];
    Kc = refCinv * exp(g_RT[0] + g_RT[2] - g_RT[3]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[2]) + (h_RT[3]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[2] -= q; /* O */
    wdot[3] += q; /* OH */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[2];
        J[0] -= dqdci;                /* dwdot[H]/d[H] */
        J[2] -= dqdci;                /* dwdot[O]/d[H] */
        J[3] += dqdci;                /* dwdot[OH]/d[H] */
        /* d()/d[H2] */
        dqdci = (2.5 - 1)*q_nocor;
        J[14] -= dqdci;               /* dwdot[H]/d[H2] */
        J[16] -= dqdci;               /* dwdot[O]/d[H2] */
        J[17] += dqdci;               /* dwdot[OH]/d[H2] */
        /* d()/d[O] */
        dqdci =  + k_f*sc[0];
        J[28] -= dqdci;               /* dwdot[H]/d[O] */
        J[30] -= dqdci;               /* dwdot[O]/d[O] */
        J[31] += dqdci;               /* dwdot[OH]/d[O] */
        /* d()/d[OH] */
        dqdci =  - k_r;
        J[42] -= dqdci;               /* dwdot[H]/d[OH] */
        J[44] -= dqdci;               /* dwdot[O]/d[OH] */
        J[45] += dqdci;               /* dwdot[OH]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (12 - 1)*q_nocor;
        J[56] -= dqdci;               /* dwdot[H]/d[H2O] */
        J[58] -= dqdci;               /* dwdot[O]/d[H2O] */
        J[59] += dqdci;               /* dwdot[OH]/d[H2O] */
        /* d()/d[AR] */
        dqdci = (0.75 - 1)*q_nocor;
        J[126] -= dqdci;              /* dwdot[H]/d[AR] */
        J[128] -= dqdci;              /* dwdot[O]/d[AR] */
        J[129] += dqdci;              /* dwdot[OH]/d[AR] */
        /* d()/d[HE] */
        dqdci = (0.75 - 1)*q_nocor;
        J[140] -= dqdci;              /* dwdot[H]/d[HE] */
        J[142] -= dqdci;              /* dwdot[O]/d[HE] */
        J[143] += dqdci;              /* dwdot[OH]/d[HE] */
        /* d()/d[CO] */
        dqdci = (1.8999999999999999 - 1)*q_nocor;
        J[154] -= dqdci;              /* dwdot[H]/d[CO] */
        J[156] -= dqdci;              /* dwdot[O]/d[CO] */
        J[157] += dqdci;              /* dwdot[OH]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (3.7999999999999998 - 1)*q_nocor;
        J[168] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[170] -= dqdci;              /* dwdot[O]/d[CO2] */
        J[171] += dqdci;              /* dwdot[OH]/d[CO2] */
    }
    else {
        dqdc[0] = q_nocor + k_f*sc[2];
        dqdc[1] = 2.5*q_nocor;
        dqdc[2] = q_nocor + k_f*sc[0];
        dqdc[3] = q_nocor - k_r;
        dqdc[4] = 12*q_nocor;
        dqdc[5] = q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = 0.75*q_nocor;
        dqdc[10] = 0.75*q_nocor;
        dqdc[11] = 1.8999999999999999*q_nocor;
        dqdc[12] = 3.7999999999999998*q_nocor;
        for (int k=0; k<13; k++) {
            J[14*k+0] -= dqdc[k];
            J[14*k+2] -= dqdc[k];
            J[14*k+3] += dqdc[k];
        }
    }
    J[182] -= dqdT; /* dwdot[H]/dT */
    J[184] -= dqdT; /* dwdot[O]/dT */
    J[185] += dqdT; /* dwdot[OH]/dT */

    /*reaction 6: H2O + M <=> H + OH + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 3 - 1)*sc[1] + ( 0 - 1)*sc[4] + ( 1.1000000000000001 - 1)*sc[10] + ( 2 - 1)*sc[8] + ( 1.5 - 1)*sc[5] + ( 1.8999999999999999 - 1)*sc[11] + ( 3.7999999999999998 - 1)*sc[12];
    /* forward */
    phi_f = sc[4];
    k_f = 1.0000000000000002e-06 * 6.0640000000000002e+27
                * exp(-3.3220000000000001 * tc[0] - 0.50321666580471969 * (120790) * invT);
    dlnkfdT = -3.3220000000000001 * invT + 0.50321666580471969 *  120790  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[3];
    Kc = refC * exp(-g_RT[0] - g_RT[3] + g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4]) + (h_RT[0] + h_RT[3]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[3] += q; /* OH */
    wdot[4] -= q; /* H2O */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H] */
        dqdci =  - k_r*sc[3];
        J[0] += dqdci;                /* dwdot[H]/d[H] */
        J[3] += dqdci;                /* dwdot[OH]/d[H] */
        J[4] -= dqdci;                /* dwdot[H2O]/d[H] */
        /* d()/d[H2] */
        dqdci = (3 - 1)*q_nocor;
        J[14] += dqdci;               /* dwdot[H]/d[H2] */
        J[17] += dqdci;               /* dwdot[OH]/d[H2] */
        J[18] -= dqdci;               /* dwdot[H2O]/d[H2] */
        /* d()/d[OH] */
        dqdci =  - k_r*sc[0];
        J[42] += dqdci;               /* dwdot[H]/d[OH] */
        J[45] += dqdci;               /* dwdot[OH]/d[OH] */
        J[46] -= dqdci;               /* dwdot[H2O]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (0 - 1)*q_nocor + k_f;
        J[56] += dqdci;               /* dwdot[H]/d[H2O] */
        J[59] += dqdci;               /* dwdot[OH]/d[H2O] */
        J[60] -= dqdci;               /* dwdot[H2O]/d[H2O] */
        /* d()/d[O2] */
        dqdci = (1.5 - 1)*q_nocor;
        J[70] += dqdci;               /* dwdot[H]/d[O2] */
        J[73] += dqdci;               /* dwdot[OH]/d[O2] */
        J[74] -= dqdci;               /* dwdot[H2O]/d[O2] */
        /* d()/d[N2] */
        dqdci = (2 - 1)*q_nocor;
        J[112] += dqdci;              /* dwdot[H]/d[N2] */
        J[115] += dqdci;              /* dwdot[OH]/d[N2] */
        J[116] -= dqdci;              /* dwdot[H2O]/d[N2] */
        /* d()/d[HE] */
        dqdci = (1.1000000000000001 - 1)*q_nocor;
        J[140] += dqdci;              /* dwdot[H]/d[HE] */
        J[143] += dqdci;              /* dwdot[OH]/d[HE] */
        J[144] -= dqdci;              /* dwdot[H2O]/d[HE] */
        /* d()/d[CO] */
        dqdci = (1.8999999999999999 - 1)*q_nocor;
        J[154] += dqdci;              /* dwdot[H]/d[CO] */
        J[157] += dqdci;              /* dwdot[OH]/d[CO] */
        J[158] -= dqdci;              /* dwdot[H2O]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (3.7999999999999998 - 1)*q_nocor;
        J[168] += dqdci;              /* dwdot[H]/d[CO2] */
        J[171] += dqdci;              /* dwdot[OH]/d[CO2] */
        J[172] -= dqdci;              /* dwdot[H2O]/d[CO2] */
    }
    else {
        dqdc[0] = q_nocor - k_r*sc[3];
        dqdc[1] = 3*q_nocor;
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor - k_r*sc[0];
        dqdc[4] =  + k_f;
        dqdc[5] = 1.5*q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = 2*q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = 1.1000000000000001*q_nocor;
        dqdc[11] = 1.8999999999999999*q_nocor;
        dqdc[12] = 3.7999999999999998*q_nocor;
        for (int k=0; k<13; k++) {
            J[14*k+0] += dqdc[k];
            J[14*k+3] += dqdc[k];
            J[14*k+4] -= dqdc[k];
        }
    }
    J[182] += dqdT; /* dwdot[H]/dT */
    J[185] += dqdT; /* dwdot[OH]/dT */
    J[186] -= dqdT; /* dwdot[H2O]/dT */

    /*reaction 7: H + O2 <=> O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[5];
    k_f = 1.0000000000000002e-06 * 104000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (15286) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  15286  * invT2;
    /* reverse */
    phi_r = sc[2]*sc[3];
    Kc = exp(g_RT[0] - g_RT[2] - g_RT[3] + g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[5]) + (h_RT[2] + h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[2] += q; /* O */
    wdot[3] += q; /* OH */
    wdot[5] -= q; /* O2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[5];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[2] += dqdci;                /* dwdot[O]/d[H] */
    J[3] += dqdci;                /* dwdot[OH]/d[H] */
    J[5] -= dqdci;                /* dwdot[O2]/d[H] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[3];
    J[28] -= dqdci;               /* dwdot[H]/d[O] */
    J[30] += dqdci;               /* dwdot[O]/d[O] */
    J[31] += dqdci;               /* dwdot[OH]/d[O] */
    J[33] -= dqdci;               /* dwdot[O2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[2];
    J[42] -= dqdci;               /* dwdot[H]/d[OH] */
    J[44] += dqdci;               /* dwdot[O]/d[OH] */
    J[45] += dqdci;               /* dwdot[OH]/d[OH] */
    J[47] -= dqdci;               /* dwdot[O2]/d[OH] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[0];
    J[70] -= dqdci;               /* dwdot[H]/d[O2] */
    J[72] += dqdci;               /* dwdot[O]/d[O2] */
    J[73] += dqdci;               /* dwdot[OH]/d[O2] */
    J[75] -= dqdci;               /* dwdot[O2]/d[O2] */
    /* d()/dT */
    J[182] -= dqdT;               /* dwdot[H]/dT */
    J[184] += dqdT;               /* dwdot[O]/dT */
    J[185] += dqdT;               /* dwdot[OH]/dT */
    J[187] -= dqdT;               /* dwdot[O2]/dT */

    /*reaction 8: O + H2 <=> H + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[2];
    k_f = 1.0000000000000002e-06 * 3818000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (7948) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  7948  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[3];
    Kc = exp(-g_RT[0] + g_RT[1] + g_RT[2] - g_RT[3]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[2]) + (h_RT[0] + h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[1] -= q; /* H2 */
    wdot[2] -= q; /* O */
    wdot[3] += q; /* OH */
    /* d()/d[H] */
    dqdci =  - k_r*sc[3];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[1] -= dqdci;                /* dwdot[H2]/d[H] */
    J[2] -= dqdci;                /* dwdot[O]/d[H] */
    J[3] += dqdci;                /* dwdot[OH]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[2];
    J[14] += dqdci;               /* dwdot[H]/d[H2] */
    J[15] -= dqdci;               /* dwdot[H2]/d[H2] */
    J[16] -= dqdci;               /* dwdot[O]/d[H2] */
    J[17] += dqdci;               /* dwdot[OH]/d[H2] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[1];
    J[28] += dqdci;               /* dwdot[H]/d[O] */
    J[29] -= dqdci;               /* dwdot[H2]/d[O] */
    J[30] -= dqdci;               /* dwdot[O]/d[O] */
    J[31] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[0];
    J[42] += dqdci;               /* dwdot[H]/d[OH] */
    J[43] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[44] -= dqdci;               /* dwdot[O]/d[OH] */
    J[45] += dqdci;               /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[182] += dqdT;               /* dwdot[H]/dT */
    J[183] -= dqdT;               /* dwdot[H2]/dT */
    J[184] -= dqdT;               /* dwdot[O]/dT */
    J[185] += dqdT;               /* dwdot[OH]/dT */

    /*reaction 9: O + H2 <=> H + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[2];
    k_f = 1.0000000000000002e-06 * 879200000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (19170) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  19170  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[3];
    Kc = exp(-g_RT[0] + g_RT[1] + g_RT[2] - g_RT[3]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[2]) + (h_RT[0] + h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[1] -= q; /* H2 */
    wdot[2] -= q; /* O */
    wdot[3] += q; /* OH */
    /* d()/d[H] */
    dqdci =  - k_r*sc[3];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[1] -= dqdci;                /* dwdot[H2]/d[H] */
    J[2] -= dqdci;                /* dwdot[O]/d[H] */
    J[3] += dqdci;                /* dwdot[OH]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[2];
    J[14] += dqdci;               /* dwdot[H]/d[H2] */
    J[15] -= dqdci;               /* dwdot[H2]/d[H2] */
    J[16] -= dqdci;               /* dwdot[O]/d[H2] */
    J[17] += dqdci;               /* dwdot[OH]/d[H2] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[1];
    J[28] += dqdci;               /* dwdot[H]/d[O] */
    J[29] -= dqdci;               /* dwdot[H2]/d[O] */
    J[30] -= dqdci;               /* dwdot[O]/d[O] */
    J[31] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[0];
    J[42] += dqdci;               /* dwdot[H]/d[OH] */
    J[43] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[44] -= dqdci;               /* dwdot[O]/d[OH] */
    J[45] += dqdci;               /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[182] += dqdT;               /* dwdot[H]/dT */
    J[183] -= dqdT;               /* dwdot[H2]/dT */
    J[184] -= dqdT;               /* dwdot[O]/dT */
    J[185] += dqdT;               /* dwdot[OH]/dT */

    /*reaction 10: H2 + OH <=> H2O + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = 1.0000000000000002e-06 * 216000000
                * exp(1.51 * tc[0] - 0.50321666580471969 * (3430) * invT);
    dlnkfdT = 1.51 * invT + 0.50321666580471969 *  3430  * invT2;
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
    wdot[3] -= q; /* OH */
    wdot[4] += q; /* H2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[4];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[1] -= dqdci;                /* dwdot[H2]/d[H] */
    J[3] -= dqdci;                /* dwdot[OH]/d[H] */
    J[4] += dqdci;                /* dwdot[H2O]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[3];
    J[14] += dqdci;               /* dwdot[H]/d[H2] */
    J[15] -= dqdci;               /* dwdot[H2]/d[H2] */
    J[17] -= dqdci;               /* dwdot[OH]/d[H2] */
    J[18] += dqdci;               /* dwdot[H2O]/d[H2] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[1];
    J[42] += dqdci;               /* dwdot[H]/d[OH] */
    J[43] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[45] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[46] += dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[0];
    J[56] += dqdci;               /* dwdot[H]/d[H2O] */
    J[57] -= dqdci;               /* dwdot[H2]/d[H2O] */
    J[59] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[60] += dqdci;               /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[182] += dqdT;               /* dwdot[H]/dT */
    J[183] -= dqdT;               /* dwdot[H2]/dT */
    J[185] -= dqdT;               /* dwdot[OH]/dT */
    J[186] += dqdT;               /* dwdot[H2O]/dT */

    /*reaction 11: OH + OH <=> O + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[3], 2.000000);
    k_f = 1.0000000000000002e-06 * 33400
                * exp(2.4199999999999999 * tc[0] - 0.50321666580471969 * (-1930) * invT);
    dlnkfdT = 2.4199999999999999 * invT + 0.50321666580471969 *  -1930  * invT2;
    /* reverse */
    phi_r = sc[2]*sc[4];
    Kc = exp(-g_RT[2] + g_RT[3] + g_RT[3] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[3]) + (h_RT[2] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* O */
    wdot[3] -= 2 * q; /* OH */
    wdot[4] += q; /* H2O */
    /* d()/d[O] */
    dqdci =  - k_r*sc[4];
    J[30] += dqdci;               /* dwdot[O]/d[O] */
    J[31] += -2 * dqdci;          /* dwdot[OH]/d[O] */
    J[32] += dqdci;               /* dwdot[H2O]/d[O] */
    /* d()/d[OH] */
    dqdci =  + k_f*2.000000*sc[3];
    J[44] += dqdci;               /* dwdot[O]/d[OH] */
    J[45] += -2 * dqdci;          /* dwdot[OH]/d[OH] */
    J[46] += dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[2];
    J[58] += dqdci;               /* dwdot[O]/d[H2O] */
    J[59] += -2 * dqdci;          /* dwdot[OH]/d[H2O] */
    J[60] += dqdci;               /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[184] += dqdT;               /* dwdot[O]/dT */
    J[185] += -2 * dqdT;          /* dwdot[OH]/dT */
    J[186] += dqdT;               /* dwdot[H2O]/dT */

    /*reaction 12: H2 + AR <=> H + H + AR */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[9];
    k_f = 1.0000000000000002e-06 * 5.84e+18
                * exp(-1.1000000000000001 * tc[0] - 0.50321666580471969 * (104380) * invT);
    dlnkfdT = -1.1000000000000001 * invT + 0.50321666580471969 *  104380  * invT2;
    /* reverse */
    phi_r = pow(sc[0], 2.000000)*sc[9];
    Kc = refC * exp(-g_RT[0] - g_RT[0] + g_RT[1] + g_RT[9] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[9]) + (2.000000*h_RT[0] + h_RT[9]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += 2 * q; /* H */
    wdot[1] -= q; /* H2 */
    /* d()/d[H] */
    dqdci =  - k_r*2.000000*sc[0]*sc[9];
    J[0] += 2 * dqdci;            /* dwdot[H]/d[H] */
    J[1] -= dqdci;                /* dwdot[H2]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[9];
    J[14] += 2 * dqdci;           /* dwdot[H]/d[H2] */
    J[15] -= dqdci;               /* dwdot[H2]/d[H2] */
    /* d()/d[AR] */
    dqdci =  + k_f*sc[1] - k_r*pow(sc[0], 2.000000);
    J[126] += 2 * dqdci;          /* dwdot[H]/d[AR] */
    J[127] -= dqdci;              /* dwdot[H2]/d[AR] */
    /* d()/dT */
    J[182] += 2 * dqdT;           /* dwdot[H]/dT */
    J[183] -= dqdT;               /* dwdot[H2]/dT */

    /*reaction 13: H2 + HE <=> H + H + HE */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[10];
    k_f = 1.0000000000000002e-06 * 5.84e+18
                * exp(-1.1000000000000001 * tc[0] - 0.50321666580471969 * (104380) * invT);
    dlnkfdT = -1.1000000000000001 * invT + 0.50321666580471969 *  104380  * invT2;
    /* reverse */
    phi_r = pow(sc[0], 2.000000)*sc[10];
    Kc = refC * exp(-g_RT[0] - g_RT[0] + g_RT[1] + g_RT[10] - g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[10]) + (2.000000*h_RT[0] + h_RT[10]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += 2 * q; /* H */
    wdot[1] -= q; /* H2 */
    /* d()/d[H] */
    dqdci =  - k_r*2.000000*sc[0]*sc[10];
    J[0] += 2 * dqdci;            /* dwdot[H]/d[H] */
    J[1] -= dqdci;                /* dwdot[H2]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[10];
    J[14] += 2 * dqdci;           /* dwdot[H]/d[H2] */
    J[15] -= dqdci;               /* dwdot[H2]/d[H2] */
    /* d()/d[HE] */
    dqdci =  + k_f*sc[1] - k_r*pow(sc[0], 2.000000);
    J[140] += 2 * dqdci;          /* dwdot[H]/d[HE] */
    J[141] -= dqdci;              /* dwdot[H2]/d[HE] */
    /* d()/dT */
    J[182] += 2 * dqdT;           /* dwdot[H]/dT */
    J[183] -= dqdT;               /* dwdot[H2]/dT */

    /*reaction 14: O + O + AR <=> O2 + AR */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[2], 2.000000)*sc[9];
    k_f = 1.0000000000000002e-12 * 18860000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-1788) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  -1788  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[9];
    Kc = refCinv * exp(g_RT[2] + g_RT[2] - g_RT[5] + g_RT[9] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[2] + h_RT[9]) + (h_RT[5] + h_RT[9]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= 2 * q; /* O */
    wdot[5] += q; /* O2 */
    /* d()/d[O] */
    dqdci =  + k_f*2.000000*sc[2]*sc[9];
    J[30] += -2 * dqdci;          /* dwdot[O]/d[O] */
    J[33] += dqdci;               /* dwdot[O2]/d[O] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[9];
    J[72] += -2 * dqdci;          /* dwdot[O]/d[O2] */
    J[75] += dqdci;               /* dwdot[O2]/d[O2] */
    /* d()/d[AR] */
    dqdci =  + k_f*pow(sc[2], 2.000000) - k_r*sc[5];
    J[128] += -2 * dqdci;         /* dwdot[O]/d[AR] */
    J[131] += dqdci;              /* dwdot[O2]/d[AR] */
    /* d()/dT */
    J[184] += -2 * dqdT;          /* dwdot[O]/dT */
    J[187] += dqdT;               /* dwdot[O2]/dT */

    /*reaction 15: O + O + HE <=> O2 + HE */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[2], 2.000000)*sc[10];
    k_f = 1.0000000000000002e-12 * 18860000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-1788) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  -1788  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[10];
    Kc = refCinv * exp(g_RT[2] + g_RT[2] - g_RT[5] + g_RT[10] - g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[2] + h_RT[10]) + (h_RT[5] + h_RT[10]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= 2 * q; /* O */
    wdot[5] += q; /* O2 */
    /* d()/d[O] */
    dqdci =  + k_f*2.000000*sc[2]*sc[10];
    J[30] += -2 * dqdci;          /* dwdot[O]/d[O] */
    J[33] += dqdci;               /* dwdot[O2]/d[O] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[10];
    J[72] += -2 * dqdci;          /* dwdot[O]/d[O2] */
    J[75] += dqdci;               /* dwdot[O2]/d[O2] */
    /* d()/d[HE] */
    dqdci =  + k_f*pow(sc[2], 2.000000) - k_r*sc[5];
    J[142] += -2 * dqdci;         /* dwdot[O]/d[HE] */
    J[145] += dqdci;              /* dwdot[O2]/d[HE] */
    /* d()/dT */
    J[184] += -2 * dqdT;          /* dwdot[O]/dT */
    J[187] += dqdT;               /* dwdot[O2]/dT */

    /*reaction 16: H2O + H2O <=> H + OH + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[4], 2.000000);
    k_f = 1.0000000000000002e-06 * 1.006e+26
                * exp(-2.4399999999999999 * tc[0] - 0.50321666580471969 * (120180) * invT);
    dlnkfdT = -2.4399999999999999 * invT + 0.50321666580471969 *  120180  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[3]*sc[4];
    Kc = refC * exp(-g_RT[0] - g_RT[3] + g_RT[4] + g_RT[4] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[4]) + (h_RT[0] + h_RT[3] + h_RT[4]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[3] += q; /* OH */
    wdot[4] -= q; /* H2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[3]*sc[4];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[3] += dqdci;                /* dwdot[OH]/d[H] */
    J[4] -= dqdci;                /* dwdot[H2O]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[0]*sc[4];
    J[42] += dqdci;               /* dwdot[H]/d[OH] */
    J[45] += dqdci;               /* dwdot[OH]/d[OH] */
    J[46] -= dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  + k_f*2.000000*sc[4] - k_r*sc[0]*sc[3];
    J[56] += dqdci;               /* dwdot[H]/d[H2O] */
    J[59] += dqdci;               /* dwdot[OH]/d[H2O] */
    J[60] -= dqdci;               /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[182] += dqdT;               /* dwdot[H]/dT */
    J[185] += dqdT;               /* dwdot[OH]/dT */
    J[186] -= dqdT;               /* dwdot[H2O]/dT */

    /*reaction 17: HO2 + H <=> H2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[6];
    k_f = 1.0000000000000002e-06 * 2750000
                * exp(2.0899999999999999 * tc[0] - 0.50321666580471969 * (-1451) * invT);
    dlnkfdT = 2.0899999999999999 * invT + 0.50321666580471969 *  -1451  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[5];
    Kc = exp(g_RT[0] - g_RT[1] - g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[6]) + (h_RT[1] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[1] += q; /* H2 */
    wdot[5] += q; /* O2 */
    wdot[6] -= q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[1] += dqdci;                /* dwdot[H2]/d[H] */
    J[5] += dqdci;                /* dwdot[O2]/d[H] */
    J[6] -= dqdci;                /* dwdot[HO2]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[5];
    J[14] -= dqdci;               /* dwdot[H]/d[H2] */
    J[15] += dqdci;               /* dwdot[H2]/d[H2] */
    J[19] += dqdci;               /* dwdot[O2]/d[H2] */
    J[20] -= dqdci;               /* dwdot[HO2]/d[H2] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[1];
    J[70] -= dqdci;               /* dwdot[H]/d[O2] */
    J[71] += dqdci;               /* dwdot[H2]/d[O2] */
    J[75] += dqdci;               /* dwdot[O2]/d[O2] */
    J[76] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[0];
    J[84] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[85] += dqdci;               /* dwdot[H2]/d[HO2] */
    J[89] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[90] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[182] -= dqdT;               /* dwdot[H]/dT */
    J[183] += dqdT;               /* dwdot[H2]/dT */
    J[187] += dqdT;               /* dwdot[O2]/dT */
    J[188] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 18: HO2 + H <=> OH + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[6];
    k_f = 1.0000000000000002e-06 * 70790000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (295) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  295  * invT2;
    /* reverse */
    phi_r = pow(sc[3], 2.000000);
    Kc = exp(g_RT[0] - g_RT[3] - g_RT[3] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[6]) + (2.000000*h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[3] += 2 * q; /* OH */
    wdot[6] -= q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[3] += 2 * dqdci;            /* dwdot[OH]/d[H] */
    J[6] -= dqdci;                /* dwdot[HO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*2.000000*sc[3];
    J[42] -= dqdci;               /* dwdot[H]/d[OH] */
    J[45] += 2 * dqdci;           /* dwdot[OH]/d[OH] */
    J[48] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[0];
    J[84] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[87] += 2 * dqdci;           /* dwdot[OH]/d[HO2] */
    J[90] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[182] -= dqdT;               /* dwdot[H]/dT */
    J[185] += 2 * dqdT;           /* dwdot[OH]/dT */
    J[188] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 19: HO2 + O <=> O2 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[6];
    k_f = 1.0000000000000002e-06 * 28500000000
                * exp(1 * tc[0] - 0.50321666580471969 * (-723.92999999999995) * invT);
    dlnkfdT = 1 * invT + 0.50321666580471969 *  -723.92999999999995  * invT2;
    /* reverse */
    phi_r = sc[3]*sc[5];
    Kc = exp(g_RT[2] - g_RT[3] - g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[6]) + (h_RT[3] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[3] += q; /* OH */
    wdot[5] += q; /* O2 */
    wdot[6] -= q; /* HO2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[6];
    J[30] -= dqdci;               /* dwdot[O]/d[O] */
    J[31] += dqdci;               /* dwdot[OH]/d[O] */
    J[33] += dqdci;               /* dwdot[O2]/d[O] */
    J[34] -= dqdci;               /* dwdot[HO2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[5];
    J[44] -= dqdci;               /* dwdot[O]/d[OH] */
    J[45] += dqdci;               /* dwdot[OH]/d[OH] */
    J[47] += dqdci;               /* dwdot[O2]/d[OH] */
    J[48] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[3];
    J[72] -= dqdci;               /* dwdot[O]/d[O2] */
    J[73] += dqdci;               /* dwdot[OH]/d[O2] */
    J[75] += dqdci;               /* dwdot[O2]/d[O2] */
    J[76] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[2];
    J[86] -= dqdci;               /* dwdot[O]/d[HO2] */
    J[87] += dqdci;               /* dwdot[OH]/d[HO2] */
    J[89] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[90] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[184] -= dqdT;               /* dwdot[O]/dT */
    J[185] += dqdT;               /* dwdot[OH]/dT */
    J[187] += dqdT;               /* dwdot[O2]/dT */
    J[188] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 20: HO2 + OH <=> H2O + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[6];
    k_f = 1.0000000000000002e-06 * 28900000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-497) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  -497  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[5];
    Kc = exp(g_RT[3] - g_RT[4] - g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[6]) + (h_RT[4] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[4] += q; /* H2O */
    wdot[5] += q; /* O2 */
    wdot[6] -= q; /* HO2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[6];
    J[45] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[46] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[47] += dqdci;               /* dwdot[O2]/d[OH] */
    J[48] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[5];
    J[59] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[60] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[61] += dqdci;               /* dwdot[O2]/d[H2O] */
    J[62] -= dqdci;               /* dwdot[HO2]/d[H2O] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[4];
    J[73] -= dqdci;               /* dwdot[OH]/d[O2] */
    J[74] += dqdci;               /* dwdot[H2O]/d[O2] */
    J[75] += dqdci;               /* dwdot[O2]/d[O2] */
    J[76] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[3];
    J[87] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[88] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[89] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[90] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[185] -= dqdT;               /* dwdot[OH]/dT */
    J[186] += dqdT;               /* dwdot[H2O]/dT */
    J[187] += dqdT;               /* dwdot[O2]/dT */
    J[188] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 21: HO2 + HO2 <=> H2O2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[6], 2.000000);
    k_f = 1.0000000000000002e-06 * 420000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (11982) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  11982  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[7];
    Kc = exp(-g_RT[5] + g_RT[6] + g_RT[6] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[6]) + (h_RT[5] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[6] -= 2 * q; /* HO2 */
    wdot[7] += q; /* H2O2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[7];
    J[75] += dqdci;               /* dwdot[O2]/d[O2] */
    J[76] += -2 * dqdci;          /* dwdot[HO2]/d[O2] */
    J[77] += dqdci;               /* dwdot[H2O2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*2.000000*sc[6];
    J[89] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[90] += -2 * dqdci;          /* dwdot[HO2]/d[HO2] */
    J[91] += dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[5];
    J[103] += dqdci;              /* dwdot[O2]/d[H2O2] */
    J[104] += -2 * dqdci;         /* dwdot[HO2]/d[H2O2] */
    J[105] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[187] += dqdT;               /* dwdot[O2]/dT */
    J[188] += -2 * dqdT;          /* dwdot[HO2]/dT */
    J[189] += dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[6], 2.000000);
    k_f = 1.0000000000000002e-06 * 130000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-1629.3) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  -1629.3  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[7];
    Kc = exp(-g_RT[5] + g_RT[6] + g_RT[6] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[6]) + (h_RT[5] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[6] -= 2 * q; /* HO2 */
    wdot[7] += q; /* H2O2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[7];
    J[75] += dqdci;               /* dwdot[O2]/d[O2] */
    J[76] += -2 * dqdci;          /* dwdot[HO2]/d[O2] */
    J[77] += dqdci;               /* dwdot[H2O2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*2.000000*sc[6];
    J[89] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[90] += -2 * dqdci;          /* dwdot[HO2]/d[HO2] */
    J[91] += dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[5];
    J[103] += dqdci;              /* dwdot[O2]/d[H2O2] */
    J[104] += -2 * dqdci;         /* dwdot[HO2]/d[H2O2] */
    J[105] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[187] += dqdT;               /* dwdot[O2]/dT */
    J[188] += -2 * dqdT;          /* dwdot[HO2]/dT */
    J[189] += dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 23: H2O2 + H <=> H2O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[7];
    k_f = 1.0000000000000002e-06 * 24100000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (3970) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  3970  * invT2;
    /* reverse */
    phi_r = sc[3]*sc[4];
    Kc = exp(g_RT[0] - g_RT[3] - g_RT[4] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[7]) + (h_RT[3] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[3] += q; /* OH */
    wdot[4] += q; /* H2O */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[3] += dqdci;                /* dwdot[OH]/d[H] */
    J[4] += dqdci;                /* dwdot[H2O]/d[H] */
    J[7] -= dqdci;                /* dwdot[H2O2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[4];
    J[42] -= dqdci;               /* dwdot[H]/d[OH] */
    J[45] += dqdci;               /* dwdot[OH]/d[OH] */
    J[46] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[49] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[3];
    J[56] -= dqdci;               /* dwdot[H]/d[H2O] */
    J[59] += dqdci;               /* dwdot[OH]/d[H2O] */
    J[60] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[63] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[0];
    J[98] -= dqdci;               /* dwdot[H]/d[H2O2] */
    J[101] += dqdci;              /* dwdot[OH]/d[H2O2] */
    J[102] += dqdci;              /* dwdot[H2O]/d[H2O2] */
    J[105] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[182] -= dqdT;               /* dwdot[H]/dT */
    J[185] += dqdT;               /* dwdot[OH]/dT */
    J[186] += dqdT;               /* dwdot[H2O]/dT */
    J[189] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 24: H2O2 + H <=> HO2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[7];
    k_f = 1.0000000000000002e-06 * 48200000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (7950) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  7950  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[6];
    Kc = exp(g_RT[0] - g_RT[1] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[7]) + (h_RT[1] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[1] += q; /* H2 */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[1] += dqdci;                /* dwdot[H2]/d[H] */
    J[6] += dqdci;                /* dwdot[HO2]/d[H] */
    J[7] -= dqdci;                /* dwdot[H2O2]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[6];
    J[14] -= dqdci;               /* dwdot[H]/d[H2] */
    J[15] += dqdci;               /* dwdot[H2]/d[H2] */
    J[20] += dqdci;               /* dwdot[HO2]/d[H2] */
    J[21] -= dqdci;               /* dwdot[H2O2]/d[H2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[1];
    J[84] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[85] += dqdci;               /* dwdot[H2]/d[HO2] */
    J[90] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[91] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[0];
    J[98] -= dqdci;               /* dwdot[H]/d[H2O2] */
    J[99] += dqdci;               /* dwdot[H2]/d[H2O2] */
    J[104] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[105] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[182] -= dqdT;               /* dwdot[H]/dT */
    J[183] += dqdT;               /* dwdot[H2]/dT */
    J[188] += dqdT;               /* dwdot[HO2]/dT */
    J[189] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 25: H2O2 + O <=> OH + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[7];
    k_f = 1.0000000000000002e-06 * 9550000
                * exp(2 * tc[0] - 0.50321666580471969 * (3970) * invT);
    dlnkfdT = 2 * invT + 0.50321666580471969 *  3970  * invT2;
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
    J[30] -= dqdci;               /* dwdot[O]/d[O] */
    J[31] += dqdci;               /* dwdot[OH]/d[O] */
    J[34] += dqdci;               /* dwdot[HO2]/d[O] */
    J[35] -= dqdci;               /* dwdot[H2O2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[6];
    J[44] -= dqdci;               /* dwdot[O]/d[OH] */
    J[45] += dqdci;               /* dwdot[OH]/d[OH] */
    J[48] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[49] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[3];
    J[86] -= dqdci;               /* dwdot[O]/d[HO2] */
    J[87] += dqdci;               /* dwdot[OH]/d[HO2] */
    J[90] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[91] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[2];
    J[100] -= dqdci;              /* dwdot[O]/d[H2O2] */
    J[101] += dqdci;              /* dwdot[OH]/d[H2O2] */
    J[104] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[105] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[184] -= dqdT;               /* dwdot[O]/dT */
    J[185] += dqdT;               /* dwdot[OH]/dT */
    J[188] += dqdT;               /* dwdot[HO2]/dT */
    J[189] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 26: H2O2 + OH <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[7];
    k_f = 1.0000000000000002e-06 * 1740000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (318) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  318  * invT2;
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
    J[45] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[46] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[48] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[49] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[6];
    J[59] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[60] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[62] += dqdci;               /* dwdot[HO2]/d[H2O] */
    J[63] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[4];
    J[87] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[88] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[90] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[91] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[3];
    J[101] -= dqdci;              /* dwdot[OH]/d[H2O2] */
    J[102] += dqdci;              /* dwdot[H2O]/d[H2O2] */
    J[104] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[105] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[185] -= dqdT;               /* dwdot[OH]/dT */
    J[186] += dqdT;               /* dwdot[H2O]/dT */
    J[188] += dqdT;               /* dwdot[HO2]/dT */
    J[189] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 27: H2O2 + OH <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[7];
    k_f = 1.0000000000000002e-06 * 75900000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (7270) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  7270  * invT2;
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
    J[45] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[46] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[48] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[49] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[6];
    J[59] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[60] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[62] += dqdci;               /* dwdot[HO2]/d[H2O] */
    J[63] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[4];
    J[87] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[88] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[90] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[91] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[3];
    J[101] -= dqdci;              /* dwdot[OH]/d[H2O2] */
    J[102] += dqdci;              /* dwdot[H2O]/d[H2O2] */
    J[104] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[105] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[185] -= dqdT;               /* dwdot[OH]/dT */
    J[186] += dqdT;               /* dwdot[H2O]/dT */
    J[188] += dqdT;               /* dwdot[HO2]/dT */
    J[189] -= dqdT;               /* dwdot[H2O2]/dT */

    double c_R[13], dcRdT[13], e_RT[13];
    double * eh_RT;
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

    double cmix = 0.0, ehmix = 0.0, dcmixdT=0.0, dehmixdT=0.0;
    for (int k = 0; k < 13; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[182+k];
    }

    double cmixinv = 1.0/cmix;
    double tmp1 = ehmix*cmixinv;
    double tmp3 = cmixinv*T;
    double tmp2 = tmp1*tmp3;
    double dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 13; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 13; ++m) {
            dehmixdc += eh_RT[m]*J[k*14+m];
        }
        J[k*14+13] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[195] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;

return;
}
#endif


#ifndef AMREX_USE_CUDA
/*compute the reaction Jacobian on CPU */
void aJacobian(double *  J, double *  sc, double T, int consP)
{
    for (int i=0; i<196; i++) {
        J[i] = 0.0;
    }

    double wdot[13];
    for (int k=0; k<13; k++) {
        wdot[k] = 0.0;
    }

    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    double invT = 1.0 / tc[1];
    double invT2 = invT * invT;

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31446 / T;
    double refCinv = 1.0 / refC;

    /*compute the mixture concentration */
    double mixture = 0.0;
    for (int k = 0; k < 13; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    double g_RT[13];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    double h_RT[13];
    speciesEnthalpy(h_RT, tc);

    double phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    double dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    double dqdci, dcdc_fac, dqdc[13];
    double Pr, fPr, F, k_0, logPr;
    double logFcent, troe_c, troe_n, troePr_den, troePr, troe;
    double Fcent1, Fcent2, Fcent3, Fcent;
    double dlogFdc, dlogFdn, dlogFdcn_fac;
    double dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
    const double ln10 = log(10.0);
    const double log10e = 1.0/log(10.0);
    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[0][0] - 1)*sc[1] + (TB[0][1] - 1)*sc[4] + (TB[0][2] - 1)*sc[5] + (TB[0][3] - 1)*sc[11] + (TB[0][4] - 1)*sc[12] + (TB[0][5] - 1)*sc[9] + (TB[0][6] - 1)*sc[10];
    /* forward */
    phi_f = sc[0]*sc[5];
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
    Kc = refCinv * exp(g_RT[0] + g_RT[5] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[5]) + (h_RT[6]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[5] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[5];
        J[0] -= dqdci;                /* dwdot[H]/d[H] */
        J[5] -= dqdci;                /* dwdot[O2]/d[H] */
        J[6] += dqdci;                /* dwdot[HO2]/d[H] */
        /* d()/d[H2] */
        dqdci = (TB[0][0] - 1)*dcdc_fac;
        J[14] -= dqdci;               /* dwdot[H]/d[H2] */
        J[19] -= dqdci;               /* dwdot[O2]/d[H2] */
        J[20] += dqdci;               /* dwdot[HO2]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[0][1] - 1)*dcdc_fac;
        J[56] -= dqdci;               /* dwdot[H]/d[H2O] */
        J[61] -= dqdci;               /* dwdot[O2]/d[H2O] */
        J[62] += dqdci;               /* dwdot[HO2]/d[H2O] */
        /* d()/d[O2] */
        dqdci = (TB[0][2] - 1)*dcdc_fac + k_f*sc[0];
        J[70] -= dqdci;               /* dwdot[H]/d[O2] */
        J[75] -= dqdci;               /* dwdot[O2]/d[O2] */
        J[76] += dqdci;               /* dwdot[HO2]/d[O2] */
        /* d()/d[HO2] */
        dqdci =  - k_r;
        J[84] -= dqdci;               /* dwdot[H]/d[HO2] */
        J[89] -= dqdci;               /* dwdot[O2]/d[HO2] */
        J[90] += dqdci;               /* dwdot[HO2]/d[HO2] */
        /* d()/d[AR] */
        dqdci = (TB[0][5] - 1)*dcdc_fac;
        J[126] -= dqdci;              /* dwdot[H]/d[AR] */
        J[131] -= dqdci;              /* dwdot[O2]/d[AR] */
        J[132] += dqdci;              /* dwdot[HO2]/d[AR] */
        /* d()/d[HE] */
        dqdci = (TB[0][6] - 1)*dcdc_fac;
        J[140] -= dqdci;              /* dwdot[H]/d[HE] */
        J[145] -= dqdci;              /* dwdot[O2]/d[HE] */
        J[146] += dqdci;              /* dwdot[HO2]/d[HE] */
        /* d()/d[CO] */
        dqdci = (TB[0][3] - 1)*dcdc_fac;
        J[154] -= dqdci;              /* dwdot[H]/d[CO] */
        J[159] -= dqdci;              /* dwdot[O2]/d[CO] */
        J[160] += dqdci;              /* dwdot[HO2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[0][4] - 1)*dcdc_fac;
        J[168] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[173] -= dqdci;              /* dwdot[O2]/d[CO2] */
        J[174] += dqdci;              /* dwdot[HO2]/d[CO2] */
    }
    else {
        dqdc[0] = dcdc_fac + k_f*sc[5];
        dqdc[1] = TB[0][0]*dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = TB[0][1]*dcdc_fac;
        dqdc[5] = TB[0][2]*dcdc_fac + k_f*sc[0];
        dqdc[6] = dcdc_fac - k_r;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = TB[0][5]*dcdc_fac;
        dqdc[10] = TB[0][6]*dcdc_fac;
        dqdc[11] = TB[0][3]*dcdc_fac;
        dqdc[12] = TB[0][4]*dcdc_fac;
        for (int k=0; k<13; k++) {
            J[14*k+0] -= dqdc[k];
            J[14*k+5] -= dqdc[k];
            J[14*k+6] += dqdc[k];
        }
    }
    J[182] -= dqdT; /* dwdot[H]/dT */
    J[187] -= dqdT; /* dwdot[O2]/dT */
    J[188] += dqdT; /* dwdot[HO2]/dT */

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + (TB[1][0] - 1)*sc[4] + (TB[1][1] - 1)*sc[12] + (TB[1][2] - 1)*sc[8] + (TB[1][3] - 1)*sc[5] + (TB[1][4] - 1)*sc[10] + (TB[1][5] - 1)*sc[7] + (TB[1][6] - 1)*sc[1] + (TB[1][7] - 1)*sc[11];
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
        dqdci = (TB[1][6] - 1)*dcdc_fac;
        J[17] += 2 * dqdci;           /* dwdot[OH]/d[H2] */
        J[21] -= dqdci;               /* dwdot[H2O2]/d[H2] */
        /* d()/d[OH] */
        dqdci =  - k_r*2.000000*sc[3];
        J[45] += 2 * dqdci;           /* dwdot[OH]/d[OH] */
        J[49] -= dqdci;               /* dwdot[H2O2]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (TB[1][0] - 1)*dcdc_fac;
        J[59] += 2 * dqdci;           /* dwdot[OH]/d[H2O] */
        J[63] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
        /* d()/d[O2] */
        dqdci = (TB[1][3] - 1)*dcdc_fac;
        J[73] += 2 * dqdci;           /* dwdot[OH]/d[O2] */
        J[77] -= dqdci;               /* dwdot[H2O2]/d[O2] */
        /* d()/d[H2O2] */
        dqdci = (TB[1][5] - 1)*dcdc_fac + k_f;
        J[101] += 2 * dqdci;          /* dwdot[OH]/d[H2O2] */
        J[105] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
        /* d()/d[N2] */
        dqdci = (TB[1][2] - 1)*dcdc_fac;
        J[115] += 2 * dqdci;          /* dwdot[OH]/d[N2] */
        J[119] -= dqdci;              /* dwdot[H2O2]/d[N2] */
        /* d()/d[HE] */
        dqdci = (TB[1][4] - 1)*dcdc_fac;
        J[143] += 2 * dqdci;          /* dwdot[OH]/d[HE] */
        J[147] -= dqdci;              /* dwdot[H2O2]/d[HE] */
        /* d()/d[CO] */
        dqdci = (TB[1][7] - 1)*dcdc_fac;
        J[157] += 2 * dqdci;          /* dwdot[OH]/d[CO] */
        J[161] -= dqdci;              /* dwdot[H2O2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[1][1] - 1)*dcdc_fac;
        J[171] += 2 * dqdci;          /* dwdot[OH]/d[CO2] */
        J[175] -= dqdci;              /* dwdot[H2O2]/d[CO2] */
    }
    else {
        dqdc[0] = dcdc_fac;
        dqdc[1] = TB[1][6]*dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac - k_r*2.000000*sc[3];
        dqdc[4] = TB[1][0]*dcdc_fac;
        dqdc[5] = TB[1][3]*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = TB[1][5]*dcdc_fac + k_f;
        dqdc[8] = TB[1][2]*dcdc_fac;
        dqdc[9] = dcdc_fac;
        dqdc[10] = TB[1][4]*dcdc_fac;
        dqdc[11] = TB[1][7]*dcdc_fac;
        dqdc[12] = TB[1][1]*dcdc_fac;
        for (int k=0; k<13; k++) {
            J[14*k+3] += 2 * dqdc[k];
            J[14*k+7] -= dqdc[k];
        }
    }
    J[185] += 2 * dqdT; /* dwdot[OH]/dT */
    J[189] -= dqdT; /* dwdot[H2O2]/dT */

    /*reaction 3: H2 + M <=> H + H + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[2][0] - 1)*sc[1] + (TB[2][1] - 1)*sc[4] + (TB[2][2] - 1)*sc[11] + (TB[2][3] - 1)*sc[12] + (TB[2][4] - 1)*sc[9] + (TB[2][5] - 1)*sc[10];
    /* forward */
    phi_f = sc[1];
    k_f = prefactor_units[2] * fwd_A[2]
                * exp(fwd_beta[2] * tc[0] - activation_units[2] * fwd_Ea[2] * invT);
    dlnkfdT = fwd_beta[2] * invT + activation_units[2] * fwd_Ea[2] * invT2;
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
        dqdci = (TB[2][0] - 1)*q_nocor + k_f;
        J[14] += 2 * dqdci;           /* dwdot[H]/d[H2] */
        J[15] -= dqdci;               /* dwdot[H2]/d[H2] */
        /* d()/d[H2O] */
        dqdci = (TB[2][1] - 1)*q_nocor;
        J[56] += 2 * dqdci;           /* dwdot[H]/d[H2O] */
        J[57] -= dqdci;               /* dwdot[H2]/d[H2O] */
        /* d()/d[AR] */
        dqdci = (TB[2][4] - 1)*q_nocor;
        J[126] += 2 * dqdci;          /* dwdot[H]/d[AR] */
        J[127] -= dqdci;              /* dwdot[H2]/d[AR] */
        /* d()/d[HE] */
        dqdci = (TB[2][5] - 1)*q_nocor;
        J[140] += 2 * dqdci;          /* dwdot[H]/d[HE] */
        J[141] -= dqdci;              /* dwdot[H2]/d[HE] */
        /* d()/d[CO] */
        dqdci = (TB[2][2] - 1)*q_nocor;
        J[154] += 2 * dqdci;          /* dwdot[H]/d[CO] */
        J[155] -= dqdci;              /* dwdot[H2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[2][3] - 1)*q_nocor;
        J[168] += 2 * dqdci;          /* dwdot[H]/d[CO2] */
        J[169] -= dqdci;              /* dwdot[H2]/d[CO2] */
    }
    else {
        dqdc[0] = q_nocor - k_r*2.000000*sc[0];
        dqdc[1] = TB[2][0]*q_nocor + k_f;
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = TB[2][1]*q_nocor;
        dqdc[5] = q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = TB[2][4]*q_nocor;
        dqdc[10] = TB[2][5]*q_nocor;
        dqdc[11] = TB[2][2]*q_nocor;
        dqdc[12] = TB[2][3]*q_nocor;
        for (int k=0; k<13; k++) {
            J[14*k+0] += 2 * dqdc[k];
            J[14*k+1] -= dqdc[k];
        }
    }
    J[182] += 2 * dqdT; /* dwdot[H]/dT */
    J[183] -= dqdT; /* dwdot[H2]/dT */

    /*reaction 4: O + O + M <=> O2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[3][0] - 1)*sc[1] + (TB[3][1] - 1)*sc[4] + (TB[3][2] - 1)*sc[9] + (TB[3][3] - 1)*sc[10] + (TB[3][4] - 1)*sc[11] + (TB[3][5] - 1)*sc[12];
    /* forward */
    phi_f = pow(sc[2], 2.000000);
    k_f = prefactor_units[3] * fwd_A[3]
                * exp(fwd_beta[3] * tc[0] - activation_units[3] * fwd_Ea[3] * invT);
    dlnkfdT = fwd_beta[3] * invT + activation_units[3] * fwd_Ea[3] * invT2;
    /* reverse */
    phi_r = sc[5];
    Kc = refCinv * exp(g_RT[2] + g_RT[2] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[2]) + (h_RT[5]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= 2 * q; /* O */
    wdot[5] += q; /* O2 */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (TB[3][0] - 1)*q_nocor;
        J[16] += -2 * dqdci;          /* dwdot[O]/d[H2] */
        J[19] += dqdci;               /* dwdot[O2]/d[H2] */
        /* d()/d[O] */
        dqdci =  + k_f*2.000000*sc[2];
        J[30] += -2 * dqdci;          /* dwdot[O]/d[O] */
        J[33] += dqdci;               /* dwdot[O2]/d[O] */
        /* d()/d[H2O] */
        dqdci = (TB[3][1] - 1)*q_nocor;
        J[58] += -2 * dqdci;          /* dwdot[O]/d[H2O] */
        J[61] += dqdci;               /* dwdot[O2]/d[H2O] */
        /* d()/d[O2] */
        dqdci =  - k_r;
        J[72] += -2 * dqdci;          /* dwdot[O]/d[O2] */
        J[75] += dqdci;               /* dwdot[O2]/d[O2] */
        /* d()/d[AR] */
        dqdci = (TB[3][2] - 1)*q_nocor;
        J[128] += -2 * dqdci;         /* dwdot[O]/d[AR] */
        J[131] += dqdci;              /* dwdot[O2]/d[AR] */
        /* d()/d[HE] */
        dqdci = (TB[3][3] - 1)*q_nocor;
        J[142] += -2 * dqdci;         /* dwdot[O]/d[HE] */
        J[145] += dqdci;              /* dwdot[O2]/d[HE] */
        /* d()/d[CO] */
        dqdci = (TB[3][4] - 1)*q_nocor;
        J[156] += -2 * dqdci;         /* dwdot[O]/d[CO] */
        J[159] += dqdci;              /* dwdot[O2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[3][5] - 1)*q_nocor;
        J[170] += -2 * dqdci;         /* dwdot[O]/d[CO2] */
        J[173] += dqdci;              /* dwdot[O2]/d[CO2] */
    }
    else {
        dqdc[0] = q_nocor;
        dqdc[1] = TB[3][0]*q_nocor;
        dqdc[2] = q_nocor + k_f*2.000000*sc[2];
        dqdc[3] = q_nocor;
        dqdc[4] = TB[3][1]*q_nocor;
        dqdc[5] = q_nocor - k_r;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = TB[3][2]*q_nocor;
        dqdc[10] = TB[3][3]*q_nocor;
        dqdc[11] = TB[3][4]*q_nocor;
        dqdc[12] = TB[3][5]*q_nocor;
        for (int k=0; k<13; k++) {
            J[14*k+2] += -2 * dqdc[k];
            J[14*k+5] += dqdc[k];
        }
    }
    J[184] += -2 * dqdT; /* dwdot[O]/dT */
    J[187] += dqdT; /* dwdot[O2]/dT */

    /*reaction 5: O + H + M <=> OH + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[4][0] - 1)*sc[1] + (TB[4][1] - 1)*sc[4] + (TB[4][2] - 1)*sc[9] + (TB[4][3] - 1)*sc[10] + (TB[4][4] - 1)*sc[11] + (TB[4][5] - 1)*sc[12];
    /* forward */
    phi_f = sc[0]*sc[2];
    k_f = prefactor_units[4] * fwd_A[4]
                * exp(fwd_beta[4] * tc[0] - activation_units[4] * fwd_Ea[4] * invT);
    dlnkfdT = fwd_beta[4] * invT + activation_units[4] * fwd_Ea[4] * invT2;
    /* reverse */
    phi_r = sc[3];
    Kc = refCinv * exp(g_RT[0] + g_RT[2] - g_RT[3]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[2]) + (h_RT[3]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[2] -= q; /* O */
    wdot[3] += q; /* OH */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[2];
        J[0] -= dqdci;                /* dwdot[H]/d[H] */
        J[2] -= dqdci;                /* dwdot[O]/d[H] */
        J[3] += dqdci;                /* dwdot[OH]/d[H] */
        /* d()/d[H2] */
        dqdci = (TB[4][0] - 1)*q_nocor;
        J[14] -= dqdci;               /* dwdot[H]/d[H2] */
        J[16] -= dqdci;               /* dwdot[O]/d[H2] */
        J[17] += dqdci;               /* dwdot[OH]/d[H2] */
        /* d()/d[O] */
        dqdci =  + k_f*sc[0];
        J[28] -= dqdci;               /* dwdot[H]/d[O] */
        J[30] -= dqdci;               /* dwdot[O]/d[O] */
        J[31] += dqdci;               /* dwdot[OH]/d[O] */
        /* d()/d[OH] */
        dqdci =  - k_r;
        J[42] -= dqdci;               /* dwdot[H]/d[OH] */
        J[44] -= dqdci;               /* dwdot[O]/d[OH] */
        J[45] += dqdci;               /* dwdot[OH]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (TB[4][1] - 1)*q_nocor;
        J[56] -= dqdci;               /* dwdot[H]/d[H2O] */
        J[58] -= dqdci;               /* dwdot[O]/d[H2O] */
        J[59] += dqdci;               /* dwdot[OH]/d[H2O] */
        /* d()/d[AR] */
        dqdci = (TB[4][2] - 1)*q_nocor;
        J[126] -= dqdci;              /* dwdot[H]/d[AR] */
        J[128] -= dqdci;              /* dwdot[O]/d[AR] */
        J[129] += dqdci;              /* dwdot[OH]/d[AR] */
        /* d()/d[HE] */
        dqdci = (TB[4][3] - 1)*q_nocor;
        J[140] -= dqdci;              /* dwdot[H]/d[HE] */
        J[142] -= dqdci;              /* dwdot[O]/d[HE] */
        J[143] += dqdci;              /* dwdot[OH]/d[HE] */
        /* d()/d[CO] */
        dqdci = (TB[4][4] - 1)*q_nocor;
        J[154] -= dqdci;              /* dwdot[H]/d[CO] */
        J[156] -= dqdci;              /* dwdot[O]/d[CO] */
        J[157] += dqdci;              /* dwdot[OH]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[4][5] - 1)*q_nocor;
        J[168] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[170] -= dqdci;              /* dwdot[O]/d[CO2] */
        J[171] += dqdci;              /* dwdot[OH]/d[CO2] */
    }
    else {
        dqdc[0] = q_nocor + k_f*sc[2];
        dqdc[1] = TB[4][0]*q_nocor;
        dqdc[2] = q_nocor + k_f*sc[0];
        dqdc[3] = q_nocor - k_r;
        dqdc[4] = TB[4][1]*q_nocor;
        dqdc[5] = q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = TB[4][2]*q_nocor;
        dqdc[10] = TB[4][3]*q_nocor;
        dqdc[11] = TB[4][4]*q_nocor;
        dqdc[12] = TB[4][5]*q_nocor;
        for (int k=0; k<13; k++) {
            J[14*k+0] -= dqdc[k];
            J[14*k+2] -= dqdc[k];
            J[14*k+3] += dqdc[k];
        }
    }
    J[182] -= dqdT; /* dwdot[H]/dT */
    J[184] -= dqdT; /* dwdot[O]/dT */
    J[185] += dqdT; /* dwdot[OH]/dT */

    /*reaction 6: H2O + M <=> H + OH + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + (TB[5][0] - 1)*sc[1] + (TB[5][1] - 1)*sc[4] + (TB[5][2] - 1)*sc[10] + (TB[5][3] - 1)*sc[8] + (TB[5][4] - 1)*sc[5] + (TB[5][5] - 1)*sc[11] + (TB[5][6] - 1)*sc[12];
    /* forward */
    phi_f = sc[4];
    k_f = prefactor_units[5] * fwd_A[5]
                * exp(fwd_beta[5] * tc[0] - activation_units[5] * fwd_Ea[5] * invT);
    dlnkfdT = fwd_beta[5] * invT + activation_units[5] * fwd_Ea[5] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[3];
    Kc = refC * exp(-g_RT[0] - g_RT[3] + g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4]) + (h_RT[0] + h_RT[3]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[3] += q; /* OH */
    wdot[4] -= q; /* H2O */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H] */
        dqdci =  - k_r*sc[3];
        J[0] += dqdci;                /* dwdot[H]/d[H] */
        J[3] += dqdci;                /* dwdot[OH]/d[H] */
        J[4] -= dqdci;                /* dwdot[H2O]/d[H] */
        /* d()/d[H2] */
        dqdci = (TB[5][0] - 1)*q_nocor;
        J[14] += dqdci;               /* dwdot[H]/d[H2] */
        J[17] += dqdci;               /* dwdot[OH]/d[H2] */
        J[18] -= dqdci;               /* dwdot[H2O]/d[H2] */
        /* d()/d[OH] */
        dqdci =  - k_r*sc[0];
        J[42] += dqdci;               /* dwdot[H]/d[OH] */
        J[45] += dqdci;               /* dwdot[OH]/d[OH] */
        J[46] -= dqdci;               /* dwdot[H2O]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (TB[5][1] - 1)*q_nocor + k_f;
        J[56] += dqdci;               /* dwdot[H]/d[H2O] */
        J[59] += dqdci;               /* dwdot[OH]/d[H2O] */
        J[60] -= dqdci;               /* dwdot[H2O]/d[H2O] */
        /* d()/d[O2] */
        dqdci = (TB[5][4] - 1)*q_nocor;
        J[70] += dqdci;               /* dwdot[H]/d[O2] */
        J[73] += dqdci;               /* dwdot[OH]/d[O2] */
        J[74] -= dqdci;               /* dwdot[H2O]/d[O2] */
        /* d()/d[N2] */
        dqdci = (TB[5][3] - 1)*q_nocor;
        J[112] += dqdci;              /* dwdot[H]/d[N2] */
        J[115] += dqdci;              /* dwdot[OH]/d[N2] */
        J[116] -= dqdci;              /* dwdot[H2O]/d[N2] */
        /* d()/d[HE] */
        dqdci = (TB[5][2] - 1)*q_nocor;
        J[140] += dqdci;              /* dwdot[H]/d[HE] */
        J[143] += dqdci;              /* dwdot[OH]/d[HE] */
        J[144] -= dqdci;              /* dwdot[H2O]/d[HE] */
        /* d()/d[CO] */
        dqdci = (TB[5][5] - 1)*q_nocor;
        J[154] += dqdci;              /* dwdot[H]/d[CO] */
        J[157] += dqdci;              /* dwdot[OH]/d[CO] */
        J[158] -= dqdci;              /* dwdot[H2O]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (TB[5][6] - 1)*q_nocor;
        J[168] += dqdci;              /* dwdot[H]/d[CO2] */
        J[171] += dqdci;              /* dwdot[OH]/d[CO2] */
        J[172] -= dqdci;              /* dwdot[H2O]/d[CO2] */
    }
    else {
        dqdc[0] = q_nocor - k_r*sc[3];
        dqdc[1] = TB[5][0]*q_nocor;
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor - k_r*sc[0];
        dqdc[4] = TB[5][1]*q_nocor + k_f;
        dqdc[5] = TB[5][4]*q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = TB[5][3]*q_nocor;
        dqdc[9] = q_nocor;
        dqdc[10] = TB[5][2]*q_nocor;
        dqdc[11] = TB[5][5]*q_nocor;
        dqdc[12] = TB[5][6]*q_nocor;
        for (int k=0; k<13; k++) {
            J[14*k+0] += dqdc[k];
            J[14*k+3] += dqdc[k];
            J[14*k+4] -= dqdc[k];
        }
    }
    J[182] += dqdT; /* dwdot[H]/dT */
    J[185] += dqdT; /* dwdot[OH]/dT */
    J[186] -= dqdT; /* dwdot[H2O]/dT */

    /*reaction 7: H + O2 <=> O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[5];
    k_f = prefactor_units[6] * fwd_A[6]
                * exp(fwd_beta[6] * tc[0] - activation_units[6] * fwd_Ea[6] * invT);
    dlnkfdT = fwd_beta[6] * invT + activation_units[6] * fwd_Ea[6] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[3];
    Kc = exp(g_RT[0] - g_RT[2] - g_RT[3] + g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[5]) + (h_RT[2] + h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[2] += q; /* O */
    wdot[3] += q; /* OH */
    wdot[5] -= q; /* O2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[5];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[2] += dqdci;                /* dwdot[O]/d[H] */
    J[3] += dqdci;                /* dwdot[OH]/d[H] */
    J[5] -= dqdci;                /* dwdot[O2]/d[H] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[3];
    J[28] -= dqdci;               /* dwdot[H]/d[O] */
    J[30] += dqdci;               /* dwdot[O]/d[O] */
    J[31] += dqdci;               /* dwdot[OH]/d[O] */
    J[33] -= dqdci;               /* dwdot[O2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[2];
    J[42] -= dqdci;               /* dwdot[H]/d[OH] */
    J[44] += dqdci;               /* dwdot[O]/d[OH] */
    J[45] += dqdci;               /* dwdot[OH]/d[OH] */
    J[47] -= dqdci;               /* dwdot[O2]/d[OH] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[0];
    J[70] -= dqdci;               /* dwdot[H]/d[O2] */
    J[72] += dqdci;               /* dwdot[O]/d[O2] */
    J[73] += dqdci;               /* dwdot[OH]/d[O2] */
    J[75] -= dqdci;               /* dwdot[O2]/d[O2] */
    /* d()/dT */
    J[182] -= dqdT;               /* dwdot[H]/dT */
    J[184] += dqdT;               /* dwdot[O]/dT */
    J[185] += dqdT;               /* dwdot[OH]/dT */
    J[187] -= dqdT;               /* dwdot[O2]/dT */

    /*reaction 8: O + H2 <=> H + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[2];
    k_f = prefactor_units[7] * fwd_A[7]
                * exp(fwd_beta[7] * tc[0] - activation_units[7] * fwd_Ea[7] * invT);
    dlnkfdT = fwd_beta[7] * invT + activation_units[7] * fwd_Ea[7] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[3];
    Kc = exp(-g_RT[0] + g_RT[1] + g_RT[2] - g_RT[3]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[2]) + (h_RT[0] + h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[1] -= q; /* H2 */
    wdot[2] -= q; /* O */
    wdot[3] += q; /* OH */
    /* d()/d[H] */
    dqdci =  - k_r*sc[3];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[1] -= dqdci;                /* dwdot[H2]/d[H] */
    J[2] -= dqdci;                /* dwdot[O]/d[H] */
    J[3] += dqdci;                /* dwdot[OH]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[2];
    J[14] += dqdci;               /* dwdot[H]/d[H2] */
    J[15] -= dqdci;               /* dwdot[H2]/d[H2] */
    J[16] -= dqdci;               /* dwdot[O]/d[H2] */
    J[17] += dqdci;               /* dwdot[OH]/d[H2] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[1];
    J[28] += dqdci;               /* dwdot[H]/d[O] */
    J[29] -= dqdci;               /* dwdot[H2]/d[O] */
    J[30] -= dqdci;               /* dwdot[O]/d[O] */
    J[31] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[0];
    J[42] += dqdci;               /* dwdot[H]/d[OH] */
    J[43] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[44] -= dqdci;               /* dwdot[O]/d[OH] */
    J[45] += dqdci;               /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[182] += dqdT;               /* dwdot[H]/dT */
    J[183] -= dqdT;               /* dwdot[H2]/dT */
    J[184] -= dqdT;               /* dwdot[O]/dT */
    J[185] += dqdT;               /* dwdot[OH]/dT */

    /*reaction 9: O + H2 <=> H + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[2];
    k_f = prefactor_units[8] * fwd_A[8]
                * exp(fwd_beta[8] * tc[0] - activation_units[8] * fwd_Ea[8] * invT);
    dlnkfdT = fwd_beta[8] * invT + activation_units[8] * fwd_Ea[8] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[3];
    Kc = exp(-g_RT[0] + g_RT[1] + g_RT[2] - g_RT[3]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[2]) + (h_RT[0] + h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[1] -= q; /* H2 */
    wdot[2] -= q; /* O */
    wdot[3] += q; /* OH */
    /* d()/d[H] */
    dqdci =  - k_r*sc[3];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[1] -= dqdci;                /* dwdot[H2]/d[H] */
    J[2] -= dqdci;                /* dwdot[O]/d[H] */
    J[3] += dqdci;                /* dwdot[OH]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[2];
    J[14] += dqdci;               /* dwdot[H]/d[H2] */
    J[15] -= dqdci;               /* dwdot[H2]/d[H2] */
    J[16] -= dqdci;               /* dwdot[O]/d[H2] */
    J[17] += dqdci;               /* dwdot[OH]/d[H2] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[1];
    J[28] += dqdci;               /* dwdot[H]/d[O] */
    J[29] -= dqdci;               /* dwdot[H2]/d[O] */
    J[30] -= dqdci;               /* dwdot[O]/d[O] */
    J[31] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[0];
    J[42] += dqdci;               /* dwdot[H]/d[OH] */
    J[43] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[44] -= dqdci;               /* dwdot[O]/d[OH] */
    J[45] += dqdci;               /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[182] += dqdT;               /* dwdot[H]/dT */
    J[183] -= dqdT;               /* dwdot[H2]/dT */
    J[184] -= dqdT;               /* dwdot[O]/dT */
    J[185] += dqdT;               /* dwdot[OH]/dT */

    /*reaction 10: H2 + OH <=> H2O + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = prefactor_units[9] * fwd_A[9]
                * exp(fwd_beta[9] * tc[0] - activation_units[9] * fwd_Ea[9] * invT);
    dlnkfdT = fwd_beta[9] * invT + activation_units[9] * fwd_Ea[9] * invT2;
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
    wdot[3] -= q; /* OH */
    wdot[4] += q; /* H2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[4];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[1] -= dqdci;                /* dwdot[H2]/d[H] */
    J[3] -= dqdci;                /* dwdot[OH]/d[H] */
    J[4] += dqdci;                /* dwdot[H2O]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[3];
    J[14] += dqdci;               /* dwdot[H]/d[H2] */
    J[15] -= dqdci;               /* dwdot[H2]/d[H2] */
    J[17] -= dqdci;               /* dwdot[OH]/d[H2] */
    J[18] += dqdci;               /* dwdot[H2O]/d[H2] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[1];
    J[42] += dqdci;               /* dwdot[H]/d[OH] */
    J[43] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[45] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[46] += dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[0];
    J[56] += dqdci;               /* dwdot[H]/d[H2O] */
    J[57] -= dqdci;               /* dwdot[H2]/d[H2O] */
    J[59] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[60] += dqdci;               /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[182] += dqdT;               /* dwdot[H]/dT */
    J[183] -= dqdT;               /* dwdot[H2]/dT */
    J[185] -= dqdT;               /* dwdot[OH]/dT */
    J[186] += dqdT;               /* dwdot[H2O]/dT */

    /*reaction 11: OH + OH <=> O + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[3], 2.000000);
    k_f = prefactor_units[10] * fwd_A[10]
                * exp(fwd_beta[10] * tc[0] - activation_units[10] * fwd_Ea[10] * invT);
    dlnkfdT = fwd_beta[10] * invT + activation_units[10] * fwd_Ea[10] * invT2;
    /* reverse */
    phi_r = sc[2]*sc[4];
    Kc = exp(-g_RT[2] + g_RT[3] + g_RT[3] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[3]) + (h_RT[2] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* O */
    wdot[3] -= 2 * q; /* OH */
    wdot[4] += q; /* H2O */
    /* d()/d[O] */
    dqdci =  - k_r*sc[4];
    J[30] += dqdci;               /* dwdot[O]/d[O] */
    J[31] += -2 * dqdci;          /* dwdot[OH]/d[O] */
    J[32] += dqdci;               /* dwdot[H2O]/d[O] */
    /* d()/d[OH] */
    dqdci =  + k_f*2.000000*sc[3];
    J[44] += dqdci;               /* dwdot[O]/d[OH] */
    J[45] += -2 * dqdci;          /* dwdot[OH]/d[OH] */
    J[46] += dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[2];
    J[58] += dqdci;               /* dwdot[O]/d[H2O] */
    J[59] += -2 * dqdci;          /* dwdot[OH]/d[H2O] */
    J[60] += dqdci;               /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[184] += dqdT;               /* dwdot[O]/dT */
    J[185] += -2 * dqdT;          /* dwdot[OH]/dT */
    J[186] += dqdT;               /* dwdot[H2O]/dT */

    /*reaction 12: H2 + AR <=> H + H + AR */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[9];
    k_f = prefactor_units[11] * fwd_A[11]
                * exp(fwd_beta[11] * tc[0] - activation_units[11] * fwd_Ea[11] * invT);
    dlnkfdT = fwd_beta[11] * invT + activation_units[11] * fwd_Ea[11] * invT2;
    /* reverse */
    phi_r = pow(sc[0], 2.000000)*sc[9];
    Kc = refC * exp(-g_RT[0] - g_RT[0] + g_RT[1] + g_RT[9] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[9]) + (2.000000*h_RT[0] + h_RT[9]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += 2 * q; /* H */
    wdot[1] -= q; /* H2 */
    /* d()/d[H] */
    dqdci =  - k_r*2.000000*sc[0]*sc[9];
    J[0] += 2 * dqdci;            /* dwdot[H]/d[H] */
    J[1] -= dqdci;                /* dwdot[H2]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[9];
    J[14] += 2 * dqdci;           /* dwdot[H]/d[H2] */
    J[15] -= dqdci;               /* dwdot[H2]/d[H2] */
    /* d()/d[AR] */
    dqdci =  + k_f*sc[1] - k_r*pow(sc[0], 2.000000);
    J[126] += 2 * dqdci;          /* dwdot[H]/d[AR] */
    J[127] -= dqdci;              /* dwdot[H2]/d[AR] */
    /* d()/dT */
    J[182] += 2 * dqdT;           /* dwdot[H]/dT */
    J[183] -= dqdT;               /* dwdot[H2]/dT */

    /*reaction 13: H2 + HE <=> H + H + HE */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[10];
    k_f = prefactor_units[12] * fwd_A[12]
                * exp(fwd_beta[12] * tc[0] - activation_units[12] * fwd_Ea[12] * invT);
    dlnkfdT = fwd_beta[12] * invT + activation_units[12] * fwd_Ea[12] * invT2;
    /* reverse */
    phi_r = pow(sc[0], 2.000000)*sc[10];
    Kc = refC * exp(-g_RT[0] - g_RT[0] + g_RT[1] + g_RT[10] - g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[10]) + (2.000000*h_RT[0] + h_RT[10]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += 2 * q; /* H */
    wdot[1] -= q; /* H2 */
    /* d()/d[H] */
    dqdci =  - k_r*2.000000*sc[0]*sc[10];
    J[0] += 2 * dqdci;            /* dwdot[H]/d[H] */
    J[1] -= dqdci;                /* dwdot[H2]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[10];
    J[14] += 2 * dqdci;           /* dwdot[H]/d[H2] */
    J[15] -= dqdci;               /* dwdot[H2]/d[H2] */
    /* d()/d[HE] */
    dqdci =  + k_f*sc[1] - k_r*pow(sc[0], 2.000000);
    J[140] += 2 * dqdci;          /* dwdot[H]/d[HE] */
    J[141] -= dqdci;              /* dwdot[H2]/d[HE] */
    /* d()/dT */
    J[182] += 2 * dqdT;           /* dwdot[H]/dT */
    J[183] -= dqdT;               /* dwdot[H2]/dT */

    /*reaction 14: O + O + AR <=> O2 + AR */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[2], 2.000000)*sc[9];
    k_f = prefactor_units[13] * fwd_A[13]
                * exp(fwd_beta[13] * tc[0] - activation_units[13] * fwd_Ea[13] * invT);
    dlnkfdT = fwd_beta[13] * invT + activation_units[13] * fwd_Ea[13] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[9];
    Kc = refCinv * exp(g_RT[2] + g_RT[2] - g_RT[5] + g_RT[9] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[2] + h_RT[9]) + (h_RT[5] + h_RT[9]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= 2 * q; /* O */
    wdot[5] += q; /* O2 */
    /* d()/d[O] */
    dqdci =  + k_f*2.000000*sc[2]*sc[9];
    J[30] += -2 * dqdci;          /* dwdot[O]/d[O] */
    J[33] += dqdci;               /* dwdot[O2]/d[O] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[9];
    J[72] += -2 * dqdci;          /* dwdot[O]/d[O2] */
    J[75] += dqdci;               /* dwdot[O2]/d[O2] */
    /* d()/d[AR] */
    dqdci =  + k_f*pow(sc[2], 2.000000) - k_r*sc[5];
    J[128] += -2 * dqdci;         /* dwdot[O]/d[AR] */
    J[131] += dqdci;              /* dwdot[O2]/d[AR] */
    /* d()/dT */
    J[184] += -2 * dqdT;          /* dwdot[O]/dT */
    J[187] += dqdT;               /* dwdot[O2]/dT */

    /*reaction 15: O + O + HE <=> O2 + HE */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[2], 2.000000)*sc[10];
    k_f = prefactor_units[14] * fwd_A[14]
                * exp(fwd_beta[14] * tc[0] - activation_units[14] * fwd_Ea[14] * invT);
    dlnkfdT = fwd_beta[14] * invT + activation_units[14] * fwd_Ea[14] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[10];
    Kc = refCinv * exp(g_RT[2] + g_RT[2] - g_RT[5] + g_RT[10] - g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[2] + h_RT[10]) + (h_RT[5] + h_RT[10]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= 2 * q; /* O */
    wdot[5] += q; /* O2 */
    /* d()/d[O] */
    dqdci =  + k_f*2.000000*sc[2]*sc[10];
    J[30] += -2 * dqdci;          /* dwdot[O]/d[O] */
    J[33] += dqdci;               /* dwdot[O2]/d[O] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[10];
    J[72] += -2 * dqdci;          /* dwdot[O]/d[O2] */
    J[75] += dqdci;               /* dwdot[O2]/d[O2] */
    /* d()/d[HE] */
    dqdci =  + k_f*pow(sc[2], 2.000000) - k_r*sc[5];
    J[142] += -2 * dqdci;         /* dwdot[O]/d[HE] */
    J[145] += dqdci;              /* dwdot[O2]/d[HE] */
    /* d()/dT */
    J[184] += -2 * dqdT;          /* dwdot[O]/dT */
    J[187] += dqdT;               /* dwdot[O2]/dT */

    /*reaction 16: H2O + H2O <=> H + OH + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[4], 2.000000);
    k_f = prefactor_units[15] * fwd_A[15]
                * exp(fwd_beta[15] * tc[0] - activation_units[15] * fwd_Ea[15] * invT);
    dlnkfdT = fwd_beta[15] * invT + activation_units[15] * fwd_Ea[15] * invT2;
    /* reverse */
    phi_r = sc[0]*sc[3]*sc[4];
    Kc = refC * exp(-g_RT[0] - g_RT[3] + g_RT[4] + g_RT[4] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[4]) + (h_RT[0] + h_RT[3] + h_RT[4]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[3] += q; /* OH */
    wdot[4] -= q; /* H2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[3]*sc[4];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[3] += dqdci;                /* dwdot[OH]/d[H] */
    J[4] -= dqdci;                /* dwdot[H2O]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[0]*sc[4];
    J[42] += dqdci;               /* dwdot[H]/d[OH] */
    J[45] += dqdci;               /* dwdot[OH]/d[OH] */
    J[46] -= dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  + k_f*2.000000*sc[4] - k_r*sc[0]*sc[3];
    J[56] += dqdci;               /* dwdot[H]/d[H2O] */
    J[59] += dqdci;               /* dwdot[OH]/d[H2O] */
    J[60] -= dqdci;               /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[182] += dqdT;               /* dwdot[H]/dT */
    J[185] += dqdT;               /* dwdot[OH]/dT */
    J[186] -= dqdT;               /* dwdot[H2O]/dT */

    /*reaction 17: HO2 + H <=> H2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[6];
    k_f = prefactor_units[16] * fwd_A[16]
                * exp(fwd_beta[16] * tc[0] - activation_units[16] * fwd_Ea[16] * invT);
    dlnkfdT = fwd_beta[16] * invT + activation_units[16] * fwd_Ea[16] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[5];
    Kc = exp(g_RT[0] - g_RT[1] - g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[6]) + (h_RT[1] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[1] += q; /* H2 */
    wdot[5] += q; /* O2 */
    wdot[6] -= q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[1] += dqdci;                /* dwdot[H2]/d[H] */
    J[5] += dqdci;                /* dwdot[O2]/d[H] */
    J[6] -= dqdci;                /* dwdot[HO2]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[5];
    J[14] -= dqdci;               /* dwdot[H]/d[H2] */
    J[15] += dqdci;               /* dwdot[H2]/d[H2] */
    J[19] += dqdci;               /* dwdot[O2]/d[H2] */
    J[20] -= dqdci;               /* dwdot[HO2]/d[H2] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[1];
    J[70] -= dqdci;               /* dwdot[H]/d[O2] */
    J[71] += dqdci;               /* dwdot[H2]/d[O2] */
    J[75] += dqdci;               /* dwdot[O2]/d[O2] */
    J[76] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[0];
    J[84] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[85] += dqdci;               /* dwdot[H2]/d[HO2] */
    J[89] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[90] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[182] -= dqdT;               /* dwdot[H]/dT */
    J[183] += dqdT;               /* dwdot[H2]/dT */
    J[187] += dqdT;               /* dwdot[O2]/dT */
    J[188] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 18: HO2 + H <=> OH + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[6];
    k_f = prefactor_units[17] * fwd_A[17]
                * exp(fwd_beta[17] * tc[0] - activation_units[17] * fwd_Ea[17] * invT);
    dlnkfdT = fwd_beta[17] * invT + activation_units[17] * fwd_Ea[17] * invT2;
    /* reverse */
    phi_r = pow(sc[3], 2.000000);
    Kc = exp(g_RT[0] - g_RT[3] - g_RT[3] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[6]) + (2.000000*h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[3] += 2 * q; /* OH */
    wdot[6] -= q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[3] += 2 * dqdci;            /* dwdot[OH]/d[H] */
    J[6] -= dqdci;                /* dwdot[HO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*2.000000*sc[3];
    J[42] -= dqdci;               /* dwdot[H]/d[OH] */
    J[45] += 2 * dqdci;           /* dwdot[OH]/d[OH] */
    J[48] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[0];
    J[84] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[87] += 2 * dqdci;           /* dwdot[OH]/d[HO2] */
    J[90] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[182] -= dqdT;               /* dwdot[H]/dT */
    J[185] += 2 * dqdT;           /* dwdot[OH]/dT */
    J[188] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 19: HO2 + O <=> O2 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[6];
    k_f = prefactor_units[18] * fwd_A[18]
                * exp(fwd_beta[18] * tc[0] - activation_units[18] * fwd_Ea[18] * invT);
    dlnkfdT = fwd_beta[18] * invT + activation_units[18] * fwd_Ea[18] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[5];
    Kc = exp(g_RT[2] - g_RT[3] - g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[6]) + (h_RT[3] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[3] += q; /* OH */
    wdot[5] += q; /* O2 */
    wdot[6] -= q; /* HO2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[6];
    J[30] -= dqdci;               /* dwdot[O]/d[O] */
    J[31] += dqdci;               /* dwdot[OH]/d[O] */
    J[33] += dqdci;               /* dwdot[O2]/d[O] */
    J[34] -= dqdci;               /* dwdot[HO2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[5];
    J[44] -= dqdci;               /* dwdot[O]/d[OH] */
    J[45] += dqdci;               /* dwdot[OH]/d[OH] */
    J[47] += dqdci;               /* dwdot[O2]/d[OH] */
    J[48] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[3];
    J[72] -= dqdci;               /* dwdot[O]/d[O2] */
    J[73] += dqdci;               /* dwdot[OH]/d[O2] */
    J[75] += dqdci;               /* dwdot[O2]/d[O2] */
    J[76] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[2];
    J[86] -= dqdci;               /* dwdot[O]/d[HO2] */
    J[87] += dqdci;               /* dwdot[OH]/d[HO2] */
    J[89] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[90] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[184] -= dqdT;               /* dwdot[O]/dT */
    J[185] += dqdT;               /* dwdot[OH]/dT */
    J[187] += dqdT;               /* dwdot[O2]/dT */
    J[188] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 20: HO2 + OH <=> H2O + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[6];
    k_f = prefactor_units[19] * fwd_A[19]
                * exp(fwd_beta[19] * tc[0] - activation_units[19] * fwd_Ea[19] * invT);
    dlnkfdT = fwd_beta[19] * invT + activation_units[19] * fwd_Ea[19] * invT2;
    /* reverse */
    phi_r = sc[4]*sc[5];
    Kc = exp(g_RT[3] - g_RT[4] - g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[6]) + (h_RT[4] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[4] += q; /* H2O */
    wdot[5] += q; /* O2 */
    wdot[6] -= q; /* HO2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[6];
    J[45] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[46] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[47] += dqdci;               /* dwdot[O2]/d[OH] */
    J[48] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[5];
    J[59] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[60] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[61] += dqdci;               /* dwdot[O2]/d[H2O] */
    J[62] -= dqdci;               /* dwdot[HO2]/d[H2O] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[4];
    J[73] -= dqdci;               /* dwdot[OH]/d[O2] */
    J[74] += dqdci;               /* dwdot[H2O]/d[O2] */
    J[75] += dqdci;               /* dwdot[O2]/d[O2] */
    J[76] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[3];
    J[87] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[88] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[89] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[90] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[185] -= dqdT;               /* dwdot[OH]/dT */
    J[186] += dqdT;               /* dwdot[H2O]/dT */
    J[187] += dqdT;               /* dwdot[O2]/dT */
    J[188] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 21: HO2 + HO2 <=> H2O2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[6], 2.000000);
    k_f = prefactor_units[20] * fwd_A[20]
                * exp(fwd_beta[20] * tc[0] - activation_units[20] * fwd_Ea[20] * invT);
    dlnkfdT = fwd_beta[20] * invT + activation_units[20] * fwd_Ea[20] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[7];
    Kc = exp(-g_RT[5] + g_RT[6] + g_RT[6] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[6]) + (h_RT[5] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[6] -= 2 * q; /* HO2 */
    wdot[7] += q; /* H2O2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[7];
    J[75] += dqdci;               /* dwdot[O2]/d[O2] */
    J[76] += -2 * dqdci;          /* dwdot[HO2]/d[O2] */
    J[77] += dqdci;               /* dwdot[H2O2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*2.000000*sc[6];
    J[89] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[90] += -2 * dqdci;          /* dwdot[HO2]/d[HO2] */
    J[91] += dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[5];
    J[103] += dqdci;              /* dwdot[O2]/d[H2O2] */
    J[104] += -2 * dqdci;         /* dwdot[HO2]/d[H2O2] */
    J[105] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[187] += dqdT;               /* dwdot[O2]/dT */
    J[188] += -2 * dqdT;          /* dwdot[HO2]/dT */
    J[189] += dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[6], 2.000000);
    k_f = prefactor_units[21] * fwd_A[21]
                * exp(fwd_beta[21] * tc[0] - activation_units[21] * fwd_Ea[21] * invT);
    dlnkfdT = fwd_beta[21] * invT + activation_units[21] * fwd_Ea[21] * invT2;
    /* reverse */
    phi_r = sc[5]*sc[7];
    Kc = exp(-g_RT[5] + g_RT[6] + g_RT[6] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[6]) + (h_RT[5] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[6] -= 2 * q; /* HO2 */
    wdot[7] += q; /* H2O2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[7];
    J[75] += dqdci;               /* dwdot[O2]/d[O2] */
    J[76] += -2 * dqdci;          /* dwdot[HO2]/d[O2] */
    J[77] += dqdci;               /* dwdot[H2O2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*2.000000*sc[6];
    J[89] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[90] += -2 * dqdci;          /* dwdot[HO2]/d[HO2] */
    J[91] += dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[5];
    J[103] += dqdci;              /* dwdot[O2]/d[H2O2] */
    J[104] += -2 * dqdci;         /* dwdot[HO2]/d[H2O2] */
    J[105] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[187] += dqdT;               /* dwdot[O2]/dT */
    J[188] += -2 * dqdT;          /* dwdot[HO2]/dT */
    J[189] += dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 23: H2O2 + H <=> H2O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[7];
    k_f = prefactor_units[22] * fwd_A[22]
                * exp(fwd_beta[22] * tc[0] - activation_units[22] * fwd_Ea[22] * invT);
    dlnkfdT = fwd_beta[22] * invT + activation_units[22] * fwd_Ea[22] * invT2;
    /* reverse */
    phi_r = sc[3]*sc[4];
    Kc = exp(g_RT[0] - g_RT[3] - g_RT[4] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[7]) + (h_RT[3] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[3] += q; /* OH */
    wdot[4] += q; /* H2O */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[3] += dqdci;                /* dwdot[OH]/d[H] */
    J[4] += dqdci;                /* dwdot[H2O]/d[H] */
    J[7] -= dqdci;                /* dwdot[H2O2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[4];
    J[42] -= dqdci;               /* dwdot[H]/d[OH] */
    J[45] += dqdci;               /* dwdot[OH]/d[OH] */
    J[46] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[49] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[3];
    J[56] -= dqdci;               /* dwdot[H]/d[H2O] */
    J[59] += dqdci;               /* dwdot[OH]/d[H2O] */
    J[60] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[63] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[0];
    J[98] -= dqdci;               /* dwdot[H]/d[H2O2] */
    J[101] += dqdci;              /* dwdot[OH]/d[H2O2] */
    J[102] += dqdci;              /* dwdot[H2O]/d[H2O2] */
    J[105] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[182] -= dqdT;               /* dwdot[H]/dT */
    J[185] += dqdT;               /* dwdot[OH]/dT */
    J[186] += dqdT;               /* dwdot[H2O]/dT */
    J[189] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 24: H2O2 + H <=> HO2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[7];
    k_f = prefactor_units[23] * fwd_A[23]
                * exp(fwd_beta[23] * tc[0] - activation_units[23] * fwd_Ea[23] * invT);
    dlnkfdT = fwd_beta[23] * invT + activation_units[23] * fwd_Ea[23] * invT2;
    /* reverse */
    phi_r = sc[1]*sc[6];
    Kc = exp(g_RT[0] - g_RT[1] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[7]) + (h_RT[1] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[1] += q; /* H2 */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[1] += dqdci;                /* dwdot[H2]/d[H] */
    J[6] += dqdci;                /* dwdot[HO2]/d[H] */
    J[7] -= dqdci;                /* dwdot[H2O2]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[6];
    J[14] -= dqdci;               /* dwdot[H]/d[H2] */
    J[15] += dqdci;               /* dwdot[H2]/d[H2] */
    J[20] += dqdci;               /* dwdot[HO2]/d[H2] */
    J[21] -= dqdci;               /* dwdot[H2O2]/d[H2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[1];
    J[84] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[85] += dqdci;               /* dwdot[H2]/d[HO2] */
    J[90] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[91] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[0];
    J[98] -= dqdci;               /* dwdot[H]/d[H2O2] */
    J[99] += dqdci;               /* dwdot[H2]/d[H2O2] */
    J[104] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[105] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[182] -= dqdT;               /* dwdot[H]/dT */
    J[183] += dqdT;               /* dwdot[H2]/dT */
    J[188] += dqdT;               /* dwdot[HO2]/dT */
    J[189] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 25: H2O2 + O <=> OH + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[7];
    k_f = prefactor_units[24] * fwd_A[24]
                * exp(fwd_beta[24] * tc[0] - activation_units[24] * fwd_Ea[24] * invT);
    dlnkfdT = fwd_beta[24] * invT + activation_units[24] * fwd_Ea[24] * invT2;
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
    J[30] -= dqdci;               /* dwdot[O]/d[O] */
    J[31] += dqdci;               /* dwdot[OH]/d[O] */
    J[34] += dqdci;               /* dwdot[HO2]/d[O] */
    J[35] -= dqdci;               /* dwdot[H2O2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[6];
    J[44] -= dqdci;               /* dwdot[O]/d[OH] */
    J[45] += dqdci;               /* dwdot[OH]/d[OH] */
    J[48] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[49] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[3];
    J[86] -= dqdci;               /* dwdot[O]/d[HO2] */
    J[87] += dqdci;               /* dwdot[OH]/d[HO2] */
    J[90] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[91] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[2];
    J[100] -= dqdci;              /* dwdot[O]/d[H2O2] */
    J[101] += dqdci;              /* dwdot[OH]/d[H2O2] */
    J[104] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[105] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[184] -= dqdT;               /* dwdot[O]/dT */
    J[185] += dqdT;               /* dwdot[OH]/dT */
    J[188] += dqdT;               /* dwdot[HO2]/dT */
    J[189] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 26: H2O2 + OH <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[7];
    k_f = prefactor_units[25] * fwd_A[25]
                * exp(fwd_beta[25] * tc[0] - activation_units[25] * fwd_Ea[25] * invT);
    dlnkfdT = fwd_beta[25] * invT + activation_units[25] * fwd_Ea[25] * invT2;
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
    J[45] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[46] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[48] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[49] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[6];
    J[59] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[60] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[62] += dqdci;               /* dwdot[HO2]/d[H2O] */
    J[63] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[4];
    J[87] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[88] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[90] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[91] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[3];
    J[101] -= dqdci;              /* dwdot[OH]/d[H2O2] */
    J[102] += dqdci;              /* dwdot[H2O]/d[H2O2] */
    J[104] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[105] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[185] -= dqdT;               /* dwdot[OH]/dT */
    J[186] += dqdT;               /* dwdot[H2O]/dT */
    J[188] += dqdT;               /* dwdot[HO2]/dT */
    J[189] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 27: H2O2 + OH <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[7];
    k_f = prefactor_units[26] * fwd_A[26]
                * exp(fwd_beta[26] * tc[0] - activation_units[26] * fwd_Ea[26] * invT);
    dlnkfdT = fwd_beta[26] * invT + activation_units[26] * fwd_Ea[26] * invT2;
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
    J[45] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[46] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[48] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[49] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[6];
    J[59] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[60] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[62] += dqdci;               /* dwdot[HO2]/d[H2O] */
    J[63] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[4];
    J[87] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[88] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[90] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[91] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[3];
    J[101] -= dqdci;              /* dwdot[OH]/d[H2O2] */
    J[102] += dqdci;              /* dwdot[H2O]/d[H2O2] */
    J[104] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[105] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[185] -= dqdT;               /* dwdot[OH]/dT */
    J[186] += dqdT;               /* dwdot[H2O]/dT */
    J[188] += dqdT;               /* dwdot[HO2]/dT */
    J[189] -= dqdT;               /* dwdot[H2O2]/dT */

    double c_R[13], dcRdT[13], e_RT[13];
    double * eh_RT;
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

    double cmix = 0.0, ehmix = 0.0, dcmixdT=0.0, dehmixdT=0.0;
    for (int k = 0; k < 13; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[182+k];
    }

    double cmixinv = 1.0/cmix;
    double tmp1 = ehmix*cmixinv;
    double tmp3 = cmixinv*T;
    double tmp2 = tmp1*tmp3;
    double dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 13; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 13; ++m) {
            dehmixdc += eh_RT[m]*J[k*14+m];
        }
        J[k*14+13] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[195] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;
}
#endif


/*compute an approx to the reaction Jacobian */
AMREX_GPU_HOST_DEVICE void aJacobian_precond(double *  J, double *  sc, double T, int HP)
{
    for (int i=0; i<196; i++) {
        J[i] = 0.0;
    }

    double wdot[13];
    for (int k=0; k<13; k++) {
        wdot[k] = 0.0;
    }

    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    double invT = 1.0 / tc[1];
    double invT2 = invT * invT;

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31446 / T;
    double refCinv = 1.0 / refC;

    /*compute the mixture concentration */
    double mixture = 0.0;
    for (int k = 0; k < 13; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    double g_RT[13];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    double h_RT[13];
    speciesEnthalpy(h_RT, tc);

    double phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    double dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    double dqdci, dcdc_fac, dqdc[13];
    double Pr, fPr, F, k_0, logPr;
    double logFcent, troe_c, troe_n, troePr_den, troePr, troe;
    double Fcent1, Fcent2, Fcent3, Fcent;
    double dlogFdc, dlogFdn, dlogFdcn_fac;
    double dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
    const double ln10 = log(10.0);
    const double log10e = 1.0/log(10.0);
    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[1] + ( 14 - 1)*sc[4] + ( 0.78000000000000003 - 1)*sc[5] + ( 1.8999999999999999 - 1)*sc[11] + ( 3.7999999999999998 - 1)*sc[12] + ( 0.67000000000000004 - 1)*sc[9] + ( 0.80000000000000004 - 1)*sc[10];
    /* forward */
    phi_f = sc[0]*sc[5];
    k_f = 1.0000000000000002e-06 * 4650840000000
                * exp(0.44 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0.44 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* pressure-fall-off */
    k_0 = 6.366e+20 * exp(-1.72 * tc[0] - 0.50321666580471969 * (524.79999999999995) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -1.72 * invT + 0.50321666580471969 * (524.79999999999995) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-0.5)*exp(-T/1.0000000000000001e-30);
    Fcent2 = 0.5 * exp(-T/1e+30);
    Fcent3 = 0.;
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/1.0000000000000001e-30
        -Fcent2/1e+30
        + 0.);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[6];
    Kc = refCinv * exp(g_RT[0] + g_RT[5] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[5]) + (h_RT[6]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[5] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = 0.0;
    dqdc[0] = dcdc_fac + k_f*sc[5];
    dqdc[1] = 2*dcdc_fac;
    dqdc[2] = dcdc_fac;
    dqdc[3] = dcdc_fac;
    dqdc[4] = 14*dcdc_fac;
    dqdc[5] = 0.78000000000000003*dcdc_fac + k_f*sc[0];
    dqdc[6] = dcdc_fac - k_r;
    dqdc[7] = dcdc_fac;
    dqdc[8] = dcdc_fac;
    dqdc[9] = 0.67000000000000004*dcdc_fac;
    dqdc[10] = 0.80000000000000004*dcdc_fac;
    dqdc[11] = 1.8999999999999999*dcdc_fac;
    dqdc[12] = 3.7999999999999998*dcdc_fac;
    for (int k=0; k<13; k++) {
        J[14*k+0] -= dqdc[k];
        J[14*k+5] -= dqdc[k];
        J[14*k+6] += dqdc[k];
    }
    J[182] -= dqdT; /* dwdot[H]/dT */
    J[187] -= dqdT; /* dwdot[O2]/dT */
    J[188] += dqdT; /* dwdot[HO2]/dT */

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 7.5 - 1)*sc[4] + ( 1.6000000000000001 - 1)*sc[12] + ( 1.5 - 1)*sc[8] + ( 1.2 - 1)*sc[5] + ( 0.65000000000000002 - 1)*sc[10] + ( 7.7000000000000002 - 1)*sc[7] + ( 3.7000000000000002 - 1)*sc[1] + ( 2.7999999999999998 - 1)*sc[11];
    /* forward */
    phi_f = sc[7];
    k_f = 1 * 2000000000000
                * exp(0.90000000000000002 * tc[0] - 0.50321666580471969 * (48749) * invT);
    dlnkfdT = 0.90000000000000002 * invT + 0.50321666580471969 *  (48749)  * invT2;
    /* pressure-fall-off */
    k_0 = 2.49e+24 * exp(-2.2999999999999998 * tc[0] - 0.50321666580471969 * (48749) * invT);
    Pr = 1e-6 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -2.2999999999999998 * invT + 0.50321666580471969 * (48749) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-0.42999999999999999)*exp(-T/1.0000000000000001e-30);
    Fcent2 = 0.42999999999999999 * exp(-T/1e+30);
    Fcent3 = 0.;
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/1.0000000000000001e-30
        -Fcent2/1e+30
        + 0.);
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
    dcdc_fac = 0.0;
    dqdc[0] = dcdc_fac;
    dqdc[1] = 3.7000000000000002*dcdc_fac;
    dqdc[2] = dcdc_fac;
    dqdc[3] = dcdc_fac - k_r*2.000000*sc[3];
    dqdc[4] = 7.5*dcdc_fac;
    dqdc[5] = 1.2*dcdc_fac;
    dqdc[6] = dcdc_fac;
    dqdc[7] = 7.7000000000000002*dcdc_fac + k_f;
    dqdc[8] = 1.5*dcdc_fac;
    dqdc[9] = dcdc_fac;
    dqdc[10] = 0.65000000000000002*dcdc_fac;
    dqdc[11] = 2.7999999999999998*dcdc_fac;
    dqdc[12] = 1.6000000000000001*dcdc_fac;
    for (int k=0; k<13; k++) {
        J[14*k+3] += 2 * dqdc[k];
        J[14*k+7] -= dqdc[k];
    }
    J[185] += 2 * dqdT; /* dwdot[OH]/dT */
    J[189] -= dqdT; /* dwdot[H2O2]/dT */

    /*reaction 3: H2 + M <=> H + H + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 2.5 - 1)*sc[1] + ( 12 - 1)*sc[4] + ( 1.8999999999999999 - 1)*sc[11] + ( 3.7999999999999998 - 1)*sc[12] + ( 0 - 1)*sc[9] + ( 0 - 1)*sc[10];
    /* forward */
    phi_f = sc[1];
    k_f = 1.0000000000000002e-06 * 4.577e+19
                * exp(-1.3999999999999999 * tc[0] - 0.50321666580471969 * (104380) * invT);
    dlnkfdT = -1.3999999999999999 * invT + 0.50321666580471969 *  (104380)  * invT2;
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
    dqdc[0] = q_nocor - k_r*2.000000*sc[0];
    dqdc[1] = 2.5*q_nocor + k_f;
    dqdc[2] = q_nocor;
    dqdc[3] = q_nocor;
    dqdc[4] = 12*q_nocor;
    dqdc[5] = q_nocor;
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    dqdc[9] = 0.0;
    dqdc[10] = 0.0;
    dqdc[11] = 1.8999999999999999*q_nocor;
    dqdc[12] = 3.7999999999999998*q_nocor;
    for (int k=0; k<13; k++) {
        J[14*k+0] += 2 * dqdc[k];
        J[14*k+1] -= dqdc[k];
    }
    J[182] += 2 * dqdT; /* dwdot[H]/dT */
    J[183] -= dqdT; /* dwdot[H2]/dT */

    /*reaction 4: O + O + M <=> O2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 2.5 - 1)*sc[1] + ( 12 - 1)*sc[4] + ( 0 - 1)*sc[9] + ( 0 - 1)*sc[10] + ( 1.8999999999999999 - 1)*sc[11] + ( 3.7999999999999998 - 1)*sc[12];
    /* forward */
    phi_f = pow(sc[2], 2.000000);
    k_f = 1.0000000000000002e-12 * 6165000000000000
                * exp(-0.5 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -0.5 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[5];
    Kc = refCinv * exp(g_RT[2] + g_RT[2] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[2]) + (h_RT[5]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= 2 * q; /* O */
    wdot[5] += q; /* O2 */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    dqdc[0] = q_nocor;
    dqdc[1] = 2.5*q_nocor;
    dqdc[2] = q_nocor + k_f*2.000000*sc[2];
    dqdc[3] = q_nocor;
    dqdc[4] = 12*q_nocor;
    dqdc[5] = q_nocor - k_r;
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    dqdc[9] = 0.0;
    dqdc[10] = 0.0;
    dqdc[11] = 1.8999999999999999*q_nocor;
    dqdc[12] = 3.7999999999999998*q_nocor;
    for (int k=0; k<13; k++) {
        J[14*k+2] += -2 * dqdc[k];
        J[14*k+5] += dqdc[k];
    }
    J[184] += -2 * dqdT; /* dwdot[O]/dT */
    J[187] += dqdT; /* dwdot[O2]/dT */

    /*reaction 5: O + H + M <=> OH + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 2.5 - 1)*sc[1] + ( 12 - 1)*sc[4] + ( 0.75 - 1)*sc[9] + ( 0.75 - 1)*sc[10] + ( 1.8999999999999999 - 1)*sc[11] + ( 3.7999999999999998 - 1)*sc[12];
    /* forward */
    phi_f = sc[0]*sc[2];
    k_f = 1.0000000000000002e-12 * 4.714e+18
                * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[3];
    Kc = refCinv * exp(g_RT[0] + g_RT[2] - g_RT[3]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[2]) + (h_RT[3]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[2] -= q; /* O */
    wdot[3] += q; /* OH */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    dqdc[0] = q_nocor + k_f*sc[2];
    dqdc[1] = 2.5*q_nocor;
    dqdc[2] = q_nocor + k_f*sc[0];
    dqdc[3] = q_nocor - k_r;
    dqdc[4] = 12*q_nocor;
    dqdc[5] = q_nocor;
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    dqdc[9] = 0.75*q_nocor;
    dqdc[10] = 0.75*q_nocor;
    dqdc[11] = 1.8999999999999999*q_nocor;
    dqdc[12] = 3.7999999999999998*q_nocor;
    for (int k=0; k<13; k++) {
        J[14*k+0] -= dqdc[k];
        J[14*k+2] -= dqdc[k];
        J[14*k+3] += dqdc[k];
    }
    J[182] -= dqdT; /* dwdot[H]/dT */
    J[184] -= dqdT; /* dwdot[O]/dT */
    J[185] += dqdT; /* dwdot[OH]/dT */

    /*reaction 6: H2O + M <=> H + OH + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 3 - 1)*sc[1] + ( 0 - 1)*sc[4] + ( 1.1000000000000001 - 1)*sc[10] + ( 2 - 1)*sc[8] + ( 1.5 - 1)*sc[5] + ( 1.8999999999999999 - 1)*sc[11] + ( 3.7999999999999998 - 1)*sc[12];
    /* forward */
    phi_f = sc[4];
    k_f = 1.0000000000000002e-06 * 6.0640000000000002e+27
                * exp(-3.3220000000000001 * tc[0] - 0.50321666580471969 * (120790) * invT);
    dlnkfdT = -3.3220000000000001 * invT + 0.50321666580471969 *  (120790)  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[3];
    Kc = refC * exp(-g_RT[0] - g_RT[3] + g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4]) + (h_RT[0] + h_RT[3]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[3] += q; /* OH */
    wdot[4] -= q; /* H2O */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    dqdc[0] = q_nocor - k_r*sc[3];
    dqdc[1] = 3*q_nocor;
    dqdc[2] = q_nocor;
    dqdc[3] = q_nocor - k_r*sc[0];
    dqdc[4] =  + k_f;
    dqdc[5] = 1.5*q_nocor;
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = 2*q_nocor;
    dqdc[9] = q_nocor;
    dqdc[10] = 1.1000000000000001*q_nocor;
    dqdc[11] = 1.8999999999999999*q_nocor;
    dqdc[12] = 3.7999999999999998*q_nocor;
    for (int k=0; k<13; k++) {
        J[14*k+0] += dqdc[k];
        J[14*k+3] += dqdc[k];
        J[14*k+4] -= dqdc[k];
    }
    J[182] += dqdT; /* dwdot[H]/dT */
    J[185] += dqdT; /* dwdot[OH]/dT */
    J[186] -= dqdT; /* dwdot[H2O]/dT */

    /*reaction 7: H + O2 <=> O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[5];
    k_f = 1.0000000000000002e-06 * 104000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (15286) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (15286)  * invT2;
    /* reverse */
    phi_r = sc[2]*sc[3];
    Kc = exp(g_RT[0] - g_RT[2] - g_RT[3] + g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[5]) + (h_RT[2] + h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[2] += q; /* O */
    wdot[3] += q; /* OH */
    wdot[5] -= q; /* O2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[5];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[2] += dqdci;                /* dwdot[O]/d[H] */
    J[3] += dqdci;                /* dwdot[OH]/d[H] */
    J[5] -= dqdci;                /* dwdot[O2]/d[H] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[3];
    J[28] -= dqdci;               /* dwdot[H]/d[O] */
    J[30] += dqdci;               /* dwdot[O]/d[O] */
    J[31] += dqdci;               /* dwdot[OH]/d[O] */
    J[33] -= dqdci;               /* dwdot[O2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[2];
    J[42] -= dqdci;               /* dwdot[H]/d[OH] */
    J[44] += dqdci;               /* dwdot[O]/d[OH] */
    J[45] += dqdci;               /* dwdot[OH]/d[OH] */
    J[47] -= dqdci;               /* dwdot[O2]/d[OH] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[0];
    J[70] -= dqdci;               /* dwdot[H]/d[O2] */
    J[72] += dqdci;               /* dwdot[O]/d[O2] */
    J[73] += dqdci;               /* dwdot[OH]/d[O2] */
    J[75] -= dqdci;               /* dwdot[O2]/d[O2] */
    /* d()/dT */
    J[182] -= dqdT;               /* dwdot[H]/dT */
    J[184] += dqdT;               /* dwdot[O]/dT */
    J[185] += dqdT;               /* dwdot[OH]/dT */
    J[187] -= dqdT;               /* dwdot[O2]/dT */

    /*reaction 8: O + H2 <=> H + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[2];
    k_f = 1.0000000000000002e-06 * 3818000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (7948) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (7948)  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[3];
    Kc = exp(-g_RT[0] + g_RT[1] + g_RT[2] - g_RT[3]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[2]) + (h_RT[0] + h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[1] -= q; /* H2 */
    wdot[2] -= q; /* O */
    wdot[3] += q; /* OH */
    /* d()/d[H] */
    dqdci =  - k_r*sc[3];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[1] -= dqdci;                /* dwdot[H2]/d[H] */
    J[2] -= dqdci;                /* dwdot[O]/d[H] */
    J[3] += dqdci;                /* dwdot[OH]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[2];
    J[14] += dqdci;               /* dwdot[H]/d[H2] */
    J[15] -= dqdci;               /* dwdot[H2]/d[H2] */
    J[16] -= dqdci;               /* dwdot[O]/d[H2] */
    J[17] += dqdci;               /* dwdot[OH]/d[H2] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[1];
    J[28] += dqdci;               /* dwdot[H]/d[O] */
    J[29] -= dqdci;               /* dwdot[H2]/d[O] */
    J[30] -= dqdci;               /* dwdot[O]/d[O] */
    J[31] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[0];
    J[42] += dqdci;               /* dwdot[H]/d[OH] */
    J[43] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[44] -= dqdci;               /* dwdot[O]/d[OH] */
    J[45] += dqdci;               /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[182] += dqdT;               /* dwdot[H]/dT */
    J[183] -= dqdT;               /* dwdot[H2]/dT */
    J[184] -= dqdT;               /* dwdot[O]/dT */
    J[185] += dqdT;               /* dwdot[OH]/dT */

    /*reaction 9: O + H2 <=> H + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[2];
    k_f = 1.0000000000000002e-06 * 879200000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (19170) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (19170)  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[3];
    Kc = exp(-g_RT[0] + g_RT[1] + g_RT[2] - g_RT[3]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[2]) + (h_RT[0] + h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[1] -= q; /* H2 */
    wdot[2] -= q; /* O */
    wdot[3] += q; /* OH */
    /* d()/d[H] */
    dqdci =  - k_r*sc[3];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[1] -= dqdci;                /* dwdot[H2]/d[H] */
    J[2] -= dqdci;                /* dwdot[O]/d[H] */
    J[3] += dqdci;                /* dwdot[OH]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[2];
    J[14] += dqdci;               /* dwdot[H]/d[H2] */
    J[15] -= dqdci;               /* dwdot[H2]/d[H2] */
    J[16] -= dqdci;               /* dwdot[O]/d[H2] */
    J[17] += dqdci;               /* dwdot[OH]/d[H2] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[1];
    J[28] += dqdci;               /* dwdot[H]/d[O] */
    J[29] -= dqdci;               /* dwdot[H2]/d[O] */
    J[30] -= dqdci;               /* dwdot[O]/d[O] */
    J[31] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[0];
    J[42] += dqdci;               /* dwdot[H]/d[OH] */
    J[43] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[44] -= dqdci;               /* dwdot[O]/d[OH] */
    J[45] += dqdci;               /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[182] += dqdT;               /* dwdot[H]/dT */
    J[183] -= dqdT;               /* dwdot[H2]/dT */
    J[184] -= dqdT;               /* dwdot[O]/dT */
    J[185] += dqdT;               /* dwdot[OH]/dT */

    /*reaction 10: H2 + OH <=> H2O + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = 1.0000000000000002e-06 * 216000000
                * exp(1.51 * tc[0] - 0.50321666580471969 * (3430) * invT);
    dlnkfdT = 1.51 * invT + 0.50321666580471969 *  (3430)  * invT2;
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
    wdot[3] -= q; /* OH */
    wdot[4] += q; /* H2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[4];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[1] -= dqdci;                /* dwdot[H2]/d[H] */
    J[3] -= dqdci;                /* dwdot[OH]/d[H] */
    J[4] += dqdci;                /* dwdot[H2O]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[3];
    J[14] += dqdci;               /* dwdot[H]/d[H2] */
    J[15] -= dqdci;               /* dwdot[H2]/d[H2] */
    J[17] -= dqdci;               /* dwdot[OH]/d[H2] */
    J[18] += dqdci;               /* dwdot[H2O]/d[H2] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[1];
    J[42] += dqdci;               /* dwdot[H]/d[OH] */
    J[43] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[45] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[46] += dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[0];
    J[56] += dqdci;               /* dwdot[H]/d[H2O] */
    J[57] -= dqdci;               /* dwdot[H2]/d[H2O] */
    J[59] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[60] += dqdci;               /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[182] += dqdT;               /* dwdot[H]/dT */
    J[183] -= dqdT;               /* dwdot[H2]/dT */
    J[185] -= dqdT;               /* dwdot[OH]/dT */
    J[186] += dqdT;               /* dwdot[H2O]/dT */

    /*reaction 11: OH + OH <=> O + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[3], 2.000000);
    k_f = 1.0000000000000002e-06 * 33400
                * exp(2.4199999999999999 * tc[0] - 0.50321666580471969 * (-1930) * invT);
    dlnkfdT = 2.4199999999999999 * invT + 0.50321666580471969 *  (-1930)  * invT2;
    /* reverse */
    phi_r = sc[2]*sc[4];
    Kc = exp(-g_RT[2] + g_RT[3] + g_RT[3] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[3]) + (h_RT[2] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* O */
    wdot[3] -= 2 * q; /* OH */
    wdot[4] += q; /* H2O */
    /* d()/d[O] */
    dqdci =  - k_r*sc[4];
    J[30] += dqdci;               /* dwdot[O]/d[O] */
    J[31] += -2 * dqdci;          /* dwdot[OH]/d[O] */
    J[32] += dqdci;               /* dwdot[H2O]/d[O] */
    /* d()/d[OH] */
    dqdci =  + k_f*2.000000*sc[3];
    J[44] += dqdci;               /* dwdot[O]/d[OH] */
    J[45] += -2 * dqdci;          /* dwdot[OH]/d[OH] */
    J[46] += dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[2];
    J[58] += dqdci;               /* dwdot[O]/d[H2O] */
    J[59] += -2 * dqdci;          /* dwdot[OH]/d[H2O] */
    J[60] += dqdci;               /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[184] += dqdT;               /* dwdot[O]/dT */
    J[185] += -2 * dqdT;          /* dwdot[OH]/dT */
    J[186] += dqdT;               /* dwdot[H2O]/dT */

    /*reaction 12: H2 + AR <=> H + H + AR */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[9];
    k_f = 1.0000000000000002e-06 * 5.84e+18
                * exp(-1.1000000000000001 * tc[0] - 0.50321666580471969 * (104380) * invT);
    dlnkfdT = -1.1000000000000001 * invT + 0.50321666580471969 *  (104380)  * invT2;
    /* reverse */
    phi_r = pow(sc[0], 2.000000)*sc[9];
    Kc = refC * exp(-g_RT[0] - g_RT[0] + g_RT[1] + g_RT[9] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[9]) + (2.000000*h_RT[0] + h_RT[9]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += 2 * q; /* H */
    wdot[1] -= q; /* H2 */
    /* d()/d[H] */
    dqdci =  - k_r*2.000000*sc[0]*sc[9];
    J[0] += 2 * dqdci;            /* dwdot[H]/d[H] */
    J[1] -= dqdci;                /* dwdot[H2]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[9];
    J[14] += 2 * dqdci;           /* dwdot[H]/d[H2] */
    J[15] -= dqdci;               /* dwdot[H2]/d[H2] */
    /* d()/d[AR] */
    dqdci =  + k_f*sc[1] - k_r*pow(sc[0], 2.000000);
    J[126] += 2 * dqdci;          /* dwdot[H]/d[AR] */
    J[127] -= dqdci;              /* dwdot[H2]/d[AR] */
    /* d()/dT */
    J[182] += 2 * dqdT;           /* dwdot[H]/dT */
    J[183] -= dqdT;               /* dwdot[H2]/dT */

    /*reaction 13: H2 + HE <=> H + H + HE */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[10];
    k_f = 1.0000000000000002e-06 * 5.84e+18
                * exp(-1.1000000000000001 * tc[0] - 0.50321666580471969 * (104380) * invT);
    dlnkfdT = -1.1000000000000001 * invT + 0.50321666580471969 *  (104380)  * invT2;
    /* reverse */
    phi_r = pow(sc[0], 2.000000)*sc[10];
    Kc = refC * exp(-g_RT[0] - g_RT[0] + g_RT[1] + g_RT[10] - g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[10]) + (2.000000*h_RT[0] + h_RT[10]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += 2 * q; /* H */
    wdot[1] -= q; /* H2 */
    /* d()/d[H] */
    dqdci =  - k_r*2.000000*sc[0]*sc[10];
    J[0] += 2 * dqdci;            /* dwdot[H]/d[H] */
    J[1] -= dqdci;                /* dwdot[H2]/d[H] */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[10];
    J[14] += 2 * dqdci;           /* dwdot[H]/d[H2] */
    J[15] -= dqdci;               /* dwdot[H2]/d[H2] */
    /* d()/d[HE] */
    dqdci =  + k_f*sc[1] - k_r*pow(sc[0], 2.000000);
    J[140] += 2 * dqdci;          /* dwdot[H]/d[HE] */
    J[141] -= dqdci;              /* dwdot[H2]/d[HE] */
    /* d()/dT */
    J[182] += 2 * dqdT;           /* dwdot[H]/dT */
    J[183] -= dqdT;               /* dwdot[H2]/dT */

    /*reaction 14: O + O + AR <=> O2 + AR */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[2], 2.000000)*sc[9];
    k_f = 1.0000000000000002e-12 * 18860000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-1788) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (-1788)  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[9];
    Kc = refCinv * exp(g_RT[2] + g_RT[2] - g_RT[5] + g_RT[9] - g_RT[9]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[2] + h_RT[9]) + (h_RT[5] + h_RT[9]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= 2 * q; /* O */
    wdot[5] += q; /* O2 */
    /* d()/d[O] */
    dqdci =  + k_f*2.000000*sc[2]*sc[9];
    J[30] += -2 * dqdci;          /* dwdot[O]/d[O] */
    J[33] += dqdci;               /* dwdot[O2]/d[O] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[9];
    J[72] += -2 * dqdci;          /* dwdot[O]/d[O2] */
    J[75] += dqdci;               /* dwdot[O2]/d[O2] */
    /* d()/d[AR] */
    dqdci =  + k_f*pow(sc[2], 2.000000) - k_r*sc[5];
    J[128] += -2 * dqdci;         /* dwdot[O]/d[AR] */
    J[131] += dqdci;              /* dwdot[O2]/d[AR] */
    /* d()/dT */
    J[184] += -2 * dqdT;          /* dwdot[O]/dT */
    J[187] += dqdT;               /* dwdot[O2]/dT */

    /*reaction 15: O + O + HE <=> O2 + HE */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[2], 2.000000)*sc[10];
    k_f = 1.0000000000000002e-12 * 18860000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-1788) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (-1788)  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[10];
    Kc = refCinv * exp(g_RT[2] + g_RT[2] - g_RT[5] + g_RT[10] - g_RT[10]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[2] + h_RT[10]) + (h_RT[5] + h_RT[10]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= 2 * q; /* O */
    wdot[5] += q; /* O2 */
    /* d()/d[O] */
    dqdci =  + k_f*2.000000*sc[2]*sc[10];
    J[30] += -2 * dqdci;          /* dwdot[O]/d[O] */
    J[33] += dqdci;               /* dwdot[O2]/d[O] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[10];
    J[72] += -2 * dqdci;          /* dwdot[O]/d[O2] */
    J[75] += dqdci;               /* dwdot[O2]/d[O2] */
    /* d()/d[HE] */
    dqdci =  + k_f*pow(sc[2], 2.000000) - k_r*sc[5];
    J[142] += -2 * dqdci;         /* dwdot[O]/d[HE] */
    J[145] += dqdci;              /* dwdot[O2]/d[HE] */
    /* d()/dT */
    J[184] += -2 * dqdT;          /* dwdot[O]/dT */
    J[187] += dqdT;               /* dwdot[O2]/dT */

    /*reaction 16: H2O + H2O <=> H + OH + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[4], 2.000000);
    k_f = 1.0000000000000002e-06 * 1.006e+26
                * exp(-2.4399999999999999 * tc[0] - 0.50321666580471969 * (120180) * invT);
    dlnkfdT = -2.4399999999999999 * invT + 0.50321666580471969 *  (120180)  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[3]*sc[4];
    Kc = refC * exp(-g_RT[0] - g_RT[3] + g_RT[4] + g_RT[4] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[4]) + (h_RT[0] + h_RT[3] + h_RT[4]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H */
    wdot[3] += q; /* OH */
    wdot[4] -= q; /* H2O */
    /* d()/d[H] */
    dqdci =  - k_r*sc[3]*sc[4];
    J[0] += dqdci;                /* dwdot[H]/d[H] */
    J[3] += dqdci;                /* dwdot[OH]/d[H] */
    J[4] -= dqdci;                /* dwdot[H2O]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[0]*sc[4];
    J[42] += dqdci;               /* dwdot[H]/d[OH] */
    J[45] += dqdci;               /* dwdot[OH]/d[OH] */
    J[46] -= dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  + k_f*2.000000*sc[4] - k_r*sc[0]*sc[3];
    J[56] += dqdci;               /* dwdot[H]/d[H2O] */
    J[59] += dqdci;               /* dwdot[OH]/d[H2O] */
    J[60] -= dqdci;               /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[182] += dqdT;               /* dwdot[H]/dT */
    J[185] += dqdT;               /* dwdot[OH]/dT */
    J[186] -= dqdT;               /* dwdot[H2O]/dT */

    /*reaction 17: HO2 + H <=> H2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[6];
    k_f = 1.0000000000000002e-06 * 2750000
                * exp(2.0899999999999999 * tc[0] - 0.50321666580471969 * (-1451) * invT);
    dlnkfdT = 2.0899999999999999 * invT + 0.50321666580471969 *  (-1451)  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[5];
    Kc = exp(g_RT[0] - g_RT[1] - g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[6]) + (h_RT[1] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[1] += q; /* H2 */
    wdot[5] += q; /* O2 */
    wdot[6] -= q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[1] += dqdci;                /* dwdot[H2]/d[H] */
    J[5] += dqdci;                /* dwdot[O2]/d[H] */
    J[6] -= dqdci;                /* dwdot[HO2]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[5];
    J[14] -= dqdci;               /* dwdot[H]/d[H2] */
    J[15] += dqdci;               /* dwdot[H2]/d[H2] */
    J[19] += dqdci;               /* dwdot[O2]/d[H2] */
    J[20] -= dqdci;               /* dwdot[HO2]/d[H2] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[1];
    J[70] -= dqdci;               /* dwdot[H]/d[O2] */
    J[71] += dqdci;               /* dwdot[H2]/d[O2] */
    J[75] += dqdci;               /* dwdot[O2]/d[O2] */
    J[76] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[0];
    J[84] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[85] += dqdci;               /* dwdot[H2]/d[HO2] */
    J[89] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[90] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[182] -= dqdT;               /* dwdot[H]/dT */
    J[183] += dqdT;               /* dwdot[H2]/dT */
    J[187] += dqdT;               /* dwdot[O2]/dT */
    J[188] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 18: HO2 + H <=> OH + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[6];
    k_f = 1.0000000000000002e-06 * 70790000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (295) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (295)  * invT2;
    /* reverse */
    phi_r = pow(sc[3], 2.000000);
    Kc = exp(g_RT[0] - g_RT[3] - g_RT[3] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[6]) + (2.000000*h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[3] += 2 * q; /* OH */
    wdot[6] -= q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[3] += 2 * dqdci;            /* dwdot[OH]/d[H] */
    J[6] -= dqdci;                /* dwdot[HO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*2.000000*sc[3];
    J[42] -= dqdci;               /* dwdot[H]/d[OH] */
    J[45] += 2 * dqdci;           /* dwdot[OH]/d[OH] */
    J[48] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[0];
    J[84] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[87] += 2 * dqdci;           /* dwdot[OH]/d[HO2] */
    J[90] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[182] -= dqdT;               /* dwdot[H]/dT */
    J[185] += 2 * dqdT;           /* dwdot[OH]/dT */
    J[188] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 19: HO2 + O <=> O2 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[6];
    k_f = 1.0000000000000002e-06 * 28500000000
                * exp(1 * tc[0] - 0.50321666580471969 * (-723.92999999999995) * invT);
    dlnkfdT = 1 * invT + 0.50321666580471969 *  (-723.92999999999995)  * invT2;
    /* reverse */
    phi_r = sc[3]*sc[5];
    Kc = exp(g_RT[2] - g_RT[3] - g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[6]) + (h_RT[3] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[3] += q; /* OH */
    wdot[5] += q; /* O2 */
    wdot[6] -= q; /* HO2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[6];
    J[30] -= dqdci;               /* dwdot[O]/d[O] */
    J[31] += dqdci;               /* dwdot[OH]/d[O] */
    J[33] += dqdci;               /* dwdot[O2]/d[O] */
    J[34] -= dqdci;               /* dwdot[HO2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[5];
    J[44] -= dqdci;               /* dwdot[O]/d[OH] */
    J[45] += dqdci;               /* dwdot[OH]/d[OH] */
    J[47] += dqdci;               /* dwdot[O2]/d[OH] */
    J[48] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[3];
    J[72] -= dqdci;               /* dwdot[O]/d[O2] */
    J[73] += dqdci;               /* dwdot[OH]/d[O2] */
    J[75] += dqdci;               /* dwdot[O2]/d[O2] */
    J[76] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[2];
    J[86] -= dqdci;               /* dwdot[O]/d[HO2] */
    J[87] += dqdci;               /* dwdot[OH]/d[HO2] */
    J[89] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[90] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[184] -= dqdT;               /* dwdot[O]/dT */
    J[185] += dqdT;               /* dwdot[OH]/dT */
    J[187] += dqdT;               /* dwdot[O2]/dT */
    J[188] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 20: HO2 + OH <=> H2O + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[6];
    k_f = 1.0000000000000002e-06 * 28900000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-497) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (-497)  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[5];
    Kc = exp(g_RT[3] - g_RT[4] - g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[6]) + (h_RT[4] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] -= q; /* OH */
    wdot[4] += q; /* H2O */
    wdot[5] += q; /* O2 */
    wdot[6] -= q; /* HO2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[6];
    J[45] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[46] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[47] += dqdci;               /* dwdot[O2]/d[OH] */
    J[48] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[5];
    J[59] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[60] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[61] += dqdci;               /* dwdot[O2]/d[H2O] */
    J[62] -= dqdci;               /* dwdot[HO2]/d[H2O] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[4];
    J[73] -= dqdci;               /* dwdot[OH]/d[O2] */
    J[74] += dqdci;               /* dwdot[H2O]/d[O2] */
    J[75] += dqdci;               /* dwdot[O2]/d[O2] */
    J[76] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[3];
    J[87] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[88] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[89] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[90] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[185] -= dqdT;               /* dwdot[OH]/dT */
    J[186] += dqdT;               /* dwdot[H2O]/dT */
    J[187] += dqdT;               /* dwdot[O2]/dT */
    J[188] -= dqdT;               /* dwdot[HO2]/dT */

    /*reaction 21: HO2 + HO2 <=> H2O2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[6], 2.000000);
    k_f = 1.0000000000000002e-06 * 420000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (11982) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (11982)  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[7];
    Kc = exp(-g_RT[5] + g_RT[6] + g_RT[6] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[6]) + (h_RT[5] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[6] -= 2 * q; /* HO2 */
    wdot[7] += q; /* H2O2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[7];
    J[75] += dqdci;               /* dwdot[O2]/d[O2] */
    J[76] += -2 * dqdci;          /* dwdot[HO2]/d[O2] */
    J[77] += dqdci;               /* dwdot[H2O2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*2.000000*sc[6];
    J[89] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[90] += -2 * dqdci;          /* dwdot[HO2]/d[HO2] */
    J[91] += dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[5];
    J[103] += dqdci;              /* dwdot[O2]/d[H2O2] */
    J[104] += -2 * dqdci;         /* dwdot[HO2]/d[H2O2] */
    J[105] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[187] += dqdT;               /* dwdot[O2]/dT */
    J[188] += -2 * dqdT;          /* dwdot[HO2]/dT */
    J[189] += dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[6], 2.000000);
    k_f = 1.0000000000000002e-06 * 130000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-1629.3) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (-1629.3)  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[7];
    Kc = exp(-g_RT[5] + g_RT[6] + g_RT[6] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[6]) + (h_RT[5] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[5] += q; /* O2 */
    wdot[6] -= 2 * q; /* HO2 */
    wdot[7] += q; /* H2O2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[7];
    J[75] += dqdci;               /* dwdot[O2]/d[O2] */
    J[76] += -2 * dqdci;          /* dwdot[HO2]/d[O2] */
    J[77] += dqdci;               /* dwdot[H2O2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*2.000000*sc[6];
    J[89] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[90] += -2 * dqdci;          /* dwdot[HO2]/d[HO2] */
    J[91] += dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[5];
    J[103] += dqdci;              /* dwdot[O2]/d[H2O2] */
    J[104] += -2 * dqdci;         /* dwdot[HO2]/d[H2O2] */
    J[105] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[187] += dqdT;               /* dwdot[O2]/dT */
    J[188] += -2 * dqdT;          /* dwdot[HO2]/dT */
    J[189] += dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 23: H2O2 + H <=> H2O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[7];
    k_f = 1.0000000000000002e-06 * 24100000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (3970) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (3970)  * invT2;
    /* reverse */
    phi_r = sc[3]*sc[4];
    Kc = exp(g_RT[0] - g_RT[3] - g_RT[4] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[7]) + (h_RT[3] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[3] += q; /* OH */
    wdot[4] += q; /* H2O */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[3] += dqdci;                /* dwdot[OH]/d[H] */
    J[4] += dqdci;                /* dwdot[H2O]/d[H] */
    J[7] -= dqdci;                /* dwdot[H2O2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[4];
    J[42] -= dqdci;               /* dwdot[H]/d[OH] */
    J[45] += dqdci;               /* dwdot[OH]/d[OH] */
    J[46] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[49] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[3];
    J[56] -= dqdci;               /* dwdot[H]/d[H2O] */
    J[59] += dqdci;               /* dwdot[OH]/d[H2O] */
    J[60] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[63] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[0];
    J[98] -= dqdci;               /* dwdot[H]/d[H2O2] */
    J[101] += dqdci;              /* dwdot[OH]/d[H2O2] */
    J[102] += dqdci;              /* dwdot[H2O]/d[H2O2] */
    J[105] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[182] -= dqdT;               /* dwdot[H]/dT */
    J[185] += dqdT;               /* dwdot[OH]/dT */
    J[186] += dqdT;               /* dwdot[H2O]/dT */
    J[189] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 24: H2O2 + H <=> HO2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[7];
    k_f = 1.0000000000000002e-06 * 48200000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (7950) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (7950)  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[6];
    Kc = exp(g_RT[0] - g_RT[1] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[7]) + (h_RT[1] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H */
    wdot[1] += q; /* H2 */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[0] -= dqdci;                /* dwdot[H]/d[H] */
    J[1] += dqdci;                /* dwdot[H2]/d[H] */
    J[6] += dqdci;                /* dwdot[HO2]/d[H] */
    J[7] -= dqdci;                /* dwdot[H2O2]/d[H] */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[6];
    J[14] -= dqdci;               /* dwdot[H]/d[H2] */
    J[15] += dqdci;               /* dwdot[H2]/d[H2] */
    J[20] += dqdci;               /* dwdot[HO2]/d[H2] */
    J[21] -= dqdci;               /* dwdot[H2O2]/d[H2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[1];
    J[84] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[85] += dqdci;               /* dwdot[H2]/d[HO2] */
    J[90] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[91] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[0];
    J[98] -= dqdci;               /* dwdot[H]/d[H2O2] */
    J[99] += dqdci;               /* dwdot[H2]/d[H2O2] */
    J[104] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[105] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[182] -= dqdT;               /* dwdot[H]/dT */
    J[183] += dqdT;               /* dwdot[H2]/dT */
    J[188] += dqdT;               /* dwdot[HO2]/dT */
    J[189] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 25: H2O2 + O <=> OH + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[7];
    k_f = 1.0000000000000002e-06 * 9550000
                * exp(2 * tc[0] - 0.50321666580471969 * (3970) * invT);
    dlnkfdT = 2 * invT + 0.50321666580471969 *  (3970)  * invT2;
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
    J[30] -= dqdci;               /* dwdot[O]/d[O] */
    J[31] += dqdci;               /* dwdot[OH]/d[O] */
    J[34] += dqdci;               /* dwdot[HO2]/d[O] */
    J[35] -= dqdci;               /* dwdot[H2O2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[6];
    J[44] -= dqdci;               /* dwdot[O]/d[OH] */
    J[45] += dqdci;               /* dwdot[OH]/d[OH] */
    J[48] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[49] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[3];
    J[86] -= dqdci;               /* dwdot[O]/d[HO2] */
    J[87] += dqdci;               /* dwdot[OH]/d[HO2] */
    J[90] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[91] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[2];
    J[100] -= dqdci;              /* dwdot[O]/d[H2O2] */
    J[101] += dqdci;              /* dwdot[OH]/d[H2O2] */
    J[104] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[105] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[184] -= dqdT;               /* dwdot[O]/dT */
    J[185] += dqdT;               /* dwdot[OH]/dT */
    J[188] += dqdT;               /* dwdot[HO2]/dT */
    J[189] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 26: H2O2 + OH <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[7];
    k_f = 1.0000000000000002e-06 * 1740000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (318) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (318)  * invT2;
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
    J[45] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[46] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[48] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[49] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[6];
    J[59] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[60] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[62] += dqdci;               /* dwdot[HO2]/d[H2O] */
    J[63] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[4];
    J[87] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[88] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[90] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[91] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[3];
    J[101] -= dqdci;              /* dwdot[OH]/d[H2O2] */
    J[102] += dqdci;              /* dwdot[H2O]/d[H2O2] */
    J[104] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[105] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[185] -= dqdT;               /* dwdot[OH]/dT */
    J[186] += dqdT;               /* dwdot[H2O]/dT */
    J[188] += dqdT;               /* dwdot[HO2]/dT */
    J[189] -= dqdT;               /* dwdot[H2O2]/dT */

    /*reaction 27: H2O2 + OH <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[3]*sc[7];
    k_f = 1.0000000000000002e-06 * 75900000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (7270) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (7270)  * invT2;
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
    J[45] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[46] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[48] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[49] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[6];
    J[59] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[60] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[62] += dqdci;               /* dwdot[HO2]/d[H2O] */
    J[63] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[4];
    J[87] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[88] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[90] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[91] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[3];
    J[101] -= dqdci;              /* dwdot[OH]/d[H2O2] */
    J[102] += dqdci;              /* dwdot[H2O]/d[H2O2] */
    J[104] += dqdci;              /* dwdot[HO2]/d[H2O2] */
    J[105] -= dqdci;              /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[185] -= dqdT;               /* dwdot[OH]/dT */
    J[186] += dqdT;               /* dwdot[H2O]/dT */
    J[188] += dqdT;               /* dwdot[HO2]/dT */
    J[189] -= dqdT;               /* dwdot[H2O2]/dT */

    double c_R[13], dcRdT[13], e_RT[13];
    double * eh_RT;
    if (HP) {
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

    double cmix = 0.0, ehmix = 0.0, dcmixdT=0.0, dehmixdT=0.0;
    for (int k = 0; k < 13; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[182+k];
    }

    double cmixinv = 1.0/cmix;
    double tmp1 = ehmix*cmixinv;
    double tmp3 = cmixinv*T;
    double tmp2 = tmp1*tmp3;
    double dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 13; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 13; ++m) {
            dehmixdc += eh_RT[m]*J[k*14+m];
        }
        J[k*14+13] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[195] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;
}


/*compute d(Cp/R)/dT and d(Cv/R)/dT at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE void dcvpRdT(double * species, double *  tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
        /*species 1: H2 */
        species[1] =
            +8.24944200e-04
            -1.62860300e-06 * tc[1]
            -2.84263020e-10 * tc[2]
            +1.65394880e-12 * tc[3];
        /*species 2: O */
        species[2] =
            -1.63816600e-03
            +4.84206400e-06 * tc[1]
            -4.80852900e-09 * tc[2]
            +1.55627840e-12 * tc[3];
        /*species 3: OH */
        species[3] =
            -3.22544939e-03
            +1.30552938e-05 * tc[1]
            -1.73956093e-08 * tc[2]
            +8.24949516e-12 * tc[3];
        /*species 4: H2O */
        species[4] =
            +3.47498200e-03
            -1.27093920e-05 * tc[1]
            +2.09057430e-08 * tc[2]
            -1.00263520e-11 * tc[3];
        /*species 5: O2 */
        species[5] =
            +1.12748600e-03
            -1.15123000e-06 * tc[1]
            +3.94163100e-09 * tc[2]
            -3.50742160e-12 * tc[3];
        /*species 6: HO2 */
        species[6] =
            -4.74912051e-03
            +4.23165782e-05 * tc[1]
            -7.28291682e-08 * tc[2]
            +3.71690050e-11 * tc[3];
        /*species 7: H2O2 */
        species[7] =
            +6.56922600e-03
            -2.97002600e-07 * tc[1]
            -1.38774180e-08 * tc[2]
            +9.88606000e-12 * tc[3];
        /*species 8: N2 */
        species[8] =
            +1.40824000e-03
            -7.92644400e-06 * tc[1]
            +1.69245450e-08 * tc[2]
            -9.77942000e-12 * tc[3];
        /*species 9: AR */
        species[9] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
        /*species 10: HE */
        species[10] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
        /*species 11: CO */
        species[11] =
            +1.51194100e-03
            -7.76351000e-06 * tc[1]
            +1.67458320e-08 * tc[2]
            -9.89980400e-12 * tc[3];
        /*species 12: CO2 */
        species[12] =
            +9.92207200e-03
            -2.08182200e-05 * tc[1]
            +2.06000610e-08 * tc[2]
            -8.46912000e-12 * tc[3];
    } else {
        /*species 0: H */
        species[0] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
        /*species 1: H2 */
        species[1] =
            +7.00064400e-04
            -1.12676580e-07 * tc[1]
            -2.76947340e-11 * tc[2]
            +6.33100800e-15 * tc[3];
        /*species 2: O */
        species[2] =
            -2.75506200e-05
            -6.20560600e-09 * tc[1]
            +1.36532010e-11 * tc[2]
            -1.74722080e-15 * tc[3];
        /*species 3: OH */
        species[3] =
            +1.05650448e-03
            -5.18165516e-07 * tc[1]
            +9.15656022e-11 * tc[2]
            -5.32783504e-15 * tc[3];
        /*species 4: H2O */
        species[4] =
            +3.05629300e-03
            -1.74605200e-06 * tc[1]
            +3.60298800e-10 * tc[2]
            -2.55664720e-14 * tc[3];
        /*species 5: O2 */
        species[5] =
            +6.13519700e-04
            -2.51768400e-07 * tc[1]
            +5.32584300e-11 * tc[2]
            -4.54574000e-15 * tc[3];
        /*species 6: HO2 */
        species[6] =
            +2.23982013e-03
            -1.26731630e-06 * tc[1]
            +3.42739110e-10 * tc[2]
            -4.31634140e-14 * tc[3];
        /*species 7: H2O2 */
        species[7] =
            +4.33613600e-03
            -2.94937800e-06 * tc[1]
            +7.04671200e-10 * tc[2]
            -5.72661600e-14 * tc[3];
        /*species 8: N2 */
        species[8] =
            +1.48797700e-03
            -1.13695220e-06 * tc[1]
            +3.02911200e-10 * tc[2]
            -2.70134040e-14 * tc[3];
        /*species 9: AR */
        species[9] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
        /*species 10: HE */
        species[10] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
        /*species 11: CO */
        species[11] =
            +1.44268900e-03
            -1.12616560e-06 * tc[1]
            +3.05574300e-10 * tc[2]
            -2.76438080e-14 * tc[3];
        /*species 12: CO2 */
        species[12] =
            +3.14016900e-03
            -2.55682200e-06 * tc[1]
            +7.18199100e-10 * tc[2]
            -6.67613200e-14 * tc[3];
    }
    return;
}


/*compute the progress rate for each reaction */
AMREX_GPU_HOST_DEVICE void progressRate(double *  qdot, double *  sc, double T)
{
    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    double invT = 1.0 / tc[1];

#ifndef AMREX_USE_CUDA
    if (T != T_save)
    {
        T_save = T;
        comp_k_f(tc,invT,k_f_save);
        comp_Kc(tc,invT,Kc_save);
    }
#endif

    double q_f[27], q_r[27];
    comp_qfqr(q_f, q_r, sc, tc, invT);

    for (int i = 0; i < 27; ++i) {
        qdot[i] = q_f[i] - q_r[i];
    }

    return;
}


/*compute the progress rate for each reaction */
AMREX_GPU_HOST_DEVICE void progressRateFR(double *  q_f, double *  q_r, double *  sc, double T)
{
    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    double invT = 1.0 / tc[1];
#ifndef AMREX_USE_CUDA

    if (T != T_save)
    {
        T_save = T;
        comp_k_f(tc,invT,k_f_save);
        comp_Kc(tc,invT,Kc_save);
    }
#endif

    comp_qfqr(q_f, q_r, sc, tc, invT);

    return;
}


/*compute the equilibrium constants for each reaction */
void equilibriumConstants(double *  kc, double *  g_RT, double T)
{
    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31446 / T;

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    kc[0] = 1.0 / (refC) * exp((g_RT[0] + g_RT[5]) - (g_RT[6]));

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    kc[1] = refC * exp((g_RT[7]) - (g_RT[3] + g_RT[3]));

    /*reaction 3: H2 + M <=> H + H + M */
    kc[2] = refC * exp((g_RT[1]) - (g_RT[0] + g_RT[0]));

    /*reaction 4: O + O + M <=> O2 + M */
    kc[3] = 1.0 / (refC) * exp((g_RT[2] + g_RT[2]) - (g_RT[5]));

    /*reaction 5: O + H + M <=> OH + M */
    kc[4] = 1.0 / (refC) * exp((g_RT[2] + g_RT[0]) - (g_RT[3]));

    /*reaction 6: H2O + M <=> H + OH + M */
    kc[5] = refC * exp((g_RT[4]) - (g_RT[0] + g_RT[3]));

    /*reaction 7: H + O2 <=> O + OH */
    kc[6] = exp((g_RT[0] + g_RT[5]) - (g_RT[2] + g_RT[3]));

    /*reaction 8: O + H2 <=> H + OH */
    kc[7] = exp((g_RT[2] + g_RT[1]) - (g_RT[0] + g_RT[3]));

    /*reaction 9: O + H2 <=> H + OH */
    kc[8] = exp((g_RT[2] + g_RT[1]) - (g_RT[0] + g_RT[3]));

    /*reaction 10: H2 + OH <=> H2O + H */
    kc[9] = exp((g_RT[1] + g_RT[3]) - (g_RT[4] + g_RT[0]));

    /*reaction 11: OH + OH <=> O + H2O */
    kc[10] = exp((g_RT[3] + g_RT[3]) - (g_RT[2] + g_RT[4]));

    /*reaction 12: H2 + AR <=> H + H + AR */
    kc[11] = refC * exp((g_RT[1] + g_RT[9]) - (g_RT[0] + g_RT[0] + g_RT[9]));

    /*reaction 13: H2 + HE <=> H + H + HE */
    kc[12] = refC * exp((g_RT[1] + g_RT[10]) - (g_RT[0] + g_RT[0] + g_RT[10]));

    /*reaction 14: O + O + AR <=> O2 + AR */
    kc[13] = 1.0 / (refC) * exp((g_RT[2] + g_RT[2] + g_RT[9]) - (g_RT[5] + g_RT[9]));

    /*reaction 15: O + O + HE <=> O2 + HE */
    kc[14] = 1.0 / (refC) * exp((g_RT[2] + g_RT[2] + g_RT[10]) - (g_RT[5] + g_RT[10]));

    /*reaction 16: H2O + H2O <=> H + OH + H2O */
    kc[15] = refC * exp((g_RT[4] + g_RT[4]) - (g_RT[0] + g_RT[3] + g_RT[4]));

    /*reaction 17: HO2 + H <=> H2 + O2 */
    kc[16] = exp((g_RT[6] + g_RT[0]) - (g_RT[1] + g_RT[5]));

    /*reaction 18: HO2 + H <=> OH + OH */
    kc[17] = exp((g_RT[6] + g_RT[0]) - (g_RT[3] + g_RT[3]));

    /*reaction 19: HO2 + O <=> O2 + OH */
    kc[18] = exp((g_RT[6] + g_RT[2]) - (g_RT[5] + g_RT[3]));

    /*reaction 20: HO2 + OH <=> H2O + O2 */
    kc[19] = exp((g_RT[6] + g_RT[3]) - (g_RT[4] + g_RT[5]));

    /*reaction 21: HO2 + HO2 <=> H2O2 + O2 */
    kc[20] = exp((g_RT[6] + g_RT[6]) - (g_RT[7] + g_RT[5]));

    /*reaction 22: HO2 + HO2 <=> H2O2 + O2 */
    kc[21] = exp((g_RT[6] + g_RT[6]) - (g_RT[7] + g_RT[5]));

    /*reaction 23: H2O2 + H <=> H2O + OH */
    kc[22] = exp((g_RT[7] + g_RT[0]) - (g_RT[4] + g_RT[3]));

    /*reaction 24: H2O2 + H <=> HO2 + H2 */
    kc[23] = exp((g_RT[7] + g_RT[0]) - (g_RT[6] + g_RT[1]));

    /*reaction 25: H2O2 + O <=> OH + HO2 */
    kc[24] = exp((g_RT[7] + g_RT[2]) - (g_RT[3] + g_RT[6]));

    /*reaction 26: H2O2 + OH <=> HO2 + H2O */
    kc[25] = exp((g_RT[7] + g_RT[3]) - (g_RT[6] + g_RT[4]));

    /*reaction 27: H2O2 + OH <=> HO2 + H2O */
    kc[26] = exp((g_RT[7] + g_RT[3]) - (g_RT[6] + g_RT[4]));

    return;
}


/*compute the g/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE void gibbs(double * species, double *  tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +2.547163000000000e+04 * invT
            +2.960117600000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 1: H2 */
        species[1] =
            -1.012521000000000e+03 * invT
            +6.592218000000000e+00
            -3.298124000000000e+00 * tc[0]
            -4.124721000000000e-04 * tc[1]
            +1.357169166666667e-07 * tc[2]
            +7.896194999999999e-12 * tc[3]
            -2.067436000000000e-14 * tc[4];
        /*species 2: O */
        species[2] =
            +2.914764000000000e+04 * invT
            -1.756599999999997e-02
            -2.946429000000000e+00 * tc[0]
            +8.190830000000000e-04 * tc[1]
            -4.035053333333333e-07 * tc[2]
            +1.335702500000000e-10 * tc[3]
            -1.945348000000000e-14 * tc[4];
        /*species 3: OH */
        species[3] =
            +3.346309130000000e+03 * invT
            +4.815738570000000e+00
            -4.125305610000000e+00 * tc[0]
            +1.612724695000000e-03 * tc[1]
            -1.087941151666667e-06 * tc[2]
            +4.832113691666666e-10 * tc[3]
            -1.031186895000000e-13 * tc[4];
        /*species 4: H2O */
        species[4] =
            -3.020811000000000e+04 * invT
            +7.966090000000001e-01
            -3.386842000000000e+00 * tc[0]
            -1.737491000000000e-03 * tc[1]
            +1.059116000000000e-06 * tc[2]
            -5.807150833333333e-10 * tc[3]
            +1.253294000000000e-13 * tc[4];
        /*species 5: O2 */
        species[5] =
            -1.005249000000000e+03 * invT
            -2.821802000000000e+00
            -3.212936000000000e+00 * tc[0]
            -5.637430000000000e-04 * tc[1]
            +9.593583333333333e-08 * tc[2]
            -1.094897500000000e-10 * tc[3]
            +4.384277000000000e-14 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +2.948080400000000e+02 * invT
            +5.851355599999999e-01
            -4.301798010000000e+00 * tc[0]
            +2.374560255000000e-03 * tc[1]
            -3.526381516666666e-06 * tc[2]
            +2.023032450000000e-09 * tc[3]
            -4.646125620000001e-13 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            -1.766315000000000e+04 * invT
            -3.396609000000000e+00
            -3.388754000000000e+00 * tc[0]
            -3.284613000000000e-03 * tc[1]
            +2.475021666666666e-08 * tc[2]
            +3.854838333333333e-10 * tc[3]
            -1.235757500000000e-13 * tc[4];
        /*species 8: N2 */
        species[8] =
            -1.020900000000000e+03 * invT
            -6.516950000000001e-01
            -3.298677000000000e+00 * tc[0]
            -7.041200000000000e-04 * tc[1]
            +6.605369999999999e-07 * tc[2]
            -4.701262500000001e-10 * tc[3]
            +1.222427500000000e-13 * tc[4];
        /*species 9: AR */
        species[9] =
            -7.453750000000000e+02 * invT
            -1.866001000000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 10: HE */
        species[10] =
            -7.453750000000000e+02 * invT
            +1.584651200000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 11: CO */
        species[11] =
            -1.431054000000000e+04 * invT
            -1.586445000000000e+00
            -3.262452000000000e+00 * tc[0]
            -7.559705000000000e-04 * tc[1]
            +6.469591666666667e-07 * tc[2]
            -4.651620000000000e-10 * tc[3]
            +1.237475500000000e-13 * tc[4];
        /*species 12: CO2 */
        species[12] =
            -4.837314000000000e+04 * invT
            -7.912765000000000e+00
            -2.275725000000000e+00 * tc[0]
            -4.961036000000000e-03 * tc[1]
            +1.734851666666667e-06 * tc[2]
            -5.722239166666667e-10 * tc[3]
            +1.058640000000000e-13 * tc[4];
    } else {
        /*species 0: H */
        species[0] =
            +2.547163000000000e+04 * invT
            +2.960117600000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 1: H2 */
        species[1] =
            -8.350340000000000e+02 * invT
            +4.346533000000000e+00
            -2.991423000000000e+00 * tc[0]
            -3.500322000000000e-04 * tc[1]
            +9.389715000000000e-09 * tc[2]
            +7.692981666666667e-13 * tc[3]
            -7.913760000000000e-17 * tc[4];
        /*species 2: O */
        species[2] =
            +2.923080000000000e+04 * invT
            -2.378248000000000e+00
            -2.542060000000000e+00 * tc[0]
            +1.377531000000000e-05 * tc[1]
            +5.171338333333333e-10 * tc[2]
            -3.792555833333334e-13 * tc[3]
            +2.184026000000000e-17 * tc[4];
        /*species 3: OH */
        species[3] =
            +3.683628750000000e+03 * invT
            -2.836911870000000e+00
            -2.864728860000000e+00 * tc[0]
            -5.282522400000000e-04 * tc[1]
            +4.318045966666667e-08 * tc[2]
            -2.543488950000000e-12 * tc[3]
            +6.659793800000000e-17 * tc[4];
        /*species 4: H2O */
        species[4] =
            -2.989921000000000e+04 * invT
            -4.190671000000000e+00
            -2.672146000000000e+00 * tc[0]
            -1.528146500000000e-03 * tc[1]
            +1.455043333333333e-07 * tc[2]
            -1.000830000000000e-11 * tc[3]
            +3.195809000000000e-16 * tc[4];
        /*species 5: O2 */
        species[5] =
            -1.233930000000000e+03 * invT
            +5.084119999999999e-01
            -3.697578000000000e+00 * tc[0]
            -3.067598500000000e-04 * tc[1]
            +2.098070000000000e-08 * tc[2]
            -1.479400833333333e-12 * tc[3]
            +5.682175000000001e-17 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +1.118567130000000e+02 * invT
            +2.321087500000001e-01
            -4.017210900000000e+00 * tc[0]
            -1.119910065000000e-03 * tc[1]
            +1.056096916666667e-07 * tc[2]
            -9.520530833333334e-12 * tc[3]
            +5.395426750000000e-16 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            -1.800696000000000e+04 * invT
            +4.072030000000000e+00
            -4.573167000000000e+00 * tc[0]
            -2.168068000000000e-03 * tc[1]
            +2.457815000000000e-07 * tc[2]
            -1.957420000000000e-11 * tc[3]
            +7.158270000000000e-16 * tc[4];
        /*species 8: N2 */
        species[8] =
            -9.227977000000000e+02 * invT
            -3.053888000000000e+00
            -2.926640000000000e+00 * tc[0]
            -7.439885000000000e-04 * tc[1]
            +9.474601666666666e-08 * tc[2]
            -8.414199999999999e-12 * tc[3]
            +3.376675500000000e-16 * tc[4];
        /*species 9: AR */
        species[9] =
            -7.453750000000000e+02 * invT
            -1.866001000000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 10: HE */
        species[10] =
            -7.453750000000000e+02 * invT
            +1.584651100000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 11: CO */
        species[11] =
            -1.426835000000000e+04 * invT
            -3.083140000000000e+00
            -3.025078000000000e+00 * tc[0]
            -7.213445000000000e-04 * tc[1]
            +9.384713333333334e-08 * tc[2]
            -8.488174999999999e-12 * tc[3]
            +3.455476000000000e-16 * tc[4];
        /*species 12: CO2 */
        species[12] =
            -4.896696000000000e+04 * invT
            +5.409018900000000e+00
            -4.453623000000000e+00 * tc[0]
            -1.570084500000000e-03 * tc[1]
            +2.130685000000000e-07 * tc[2]
            -1.994997500000000e-11 * tc[3]
            +8.345165000000000e-16 * tc[4];
    }
    return;
}


/*compute the a/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE void helmholtz(double * species, double *  tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +2.54716300e+04 * invT
            +1.96011760e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 1: H2 */
        species[1] =
            -1.01252100e+03 * invT
            +5.59221800e+00
            -3.29812400e+00 * tc[0]
            -4.12472100e-04 * tc[1]
            +1.35716917e-07 * tc[2]
            +7.89619500e-12 * tc[3]
            -2.06743600e-14 * tc[4];
        /*species 2: O */
        species[2] =
            +2.91476400e+04 * invT
            -1.01756600e+00
            -2.94642900e+00 * tc[0]
            +8.19083000e-04 * tc[1]
            -4.03505333e-07 * tc[2]
            +1.33570250e-10 * tc[3]
            -1.94534800e-14 * tc[4];
        /*species 3: OH */
        species[3] =
            +3.34630913e+03 * invT
            +3.81573857e+00
            -4.12530561e+00 * tc[0]
            +1.61272470e-03 * tc[1]
            -1.08794115e-06 * tc[2]
            +4.83211369e-10 * tc[3]
            -1.03118689e-13 * tc[4];
        /*species 4: H2O */
        species[4] =
            -3.02081100e+04 * invT
            -2.03391000e-01
            -3.38684200e+00 * tc[0]
            -1.73749100e-03 * tc[1]
            +1.05911600e-06 * tc[2]
            -5.80715083e-10 * tc[3]
            +1.25329400e-13 * tc[4];
        /*species 5: O2 */
        species[5] =
            -1.00524900e+03 * invT
            -3.82180200e+00
            -3.21293600e+00 * tc[0]
            -5.63743000e-04 * tc[1]
            +9.59358333e-08 * tc[2]
            -1.09489750e-10 * tc[3]
            +4.38427700e-14 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +2.94808040e+02 * invT
            -4.14864440e-01
            -4.30179801e+00 * tc[0]
            +2.37456025e-03 * tc[1]
            -3.52638152e-06 * tc[2]
            +2.02303245e-09 * tc[3]
            -4.64612562e-13 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            -1.76631500e+04 * invT
            -4.39660900e+00
            -3.38875400e+00 * tc[0]
            -3.28461300e-03 * tc[1]
            +2.47502167e-08 * tc[2]
            +3.85483833e-10 * tc[3]
            -1.23575750e-13 * tc[4];
        /*species 8: N2 */
        species[8] =
            -1.02090000e+03 * invT
            -1.65169500e+00
            -3.29867700e+00 * tc[0]
            -7.04120000e-04 * tc[1]
            +6.60537000e-07 * tc[2]
            -4.70126250e-10 * tc[3]
            +1.22242750e-13 * tc[4];
        /*species 9: AR */
        species[9] =
            -7.45375000e+02 * invT
            -2.86600100e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 10: HE */
        species[10] =
            -7.45375000e+02 * invT
            +5.84651200e-01
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 11: CO */
        species[11] =
            -1.43105400e+04 * invT
            -2.58644500e+00
            -3.26245200e+00 * tc[0]
            -7.55970500e-04 * tc[1]
            +6.46959167e-07 * tc[2]
            -4.65162000e-10 * tc[3]
            +1.23747550e-13 * tc[4];
        /*species 12: CO2 */
        species[12] =
            -4.83731400e+04 * invT
            -8.91276500e+00
            -2.27572500e+00 * tc[0]
            -4.96103600e-03 * tc[1]
            +1.73485167e-06 * tc[2]
            -5.72223917e-10 * tc[3]
            +1.05864000e-13 * tc[4];
    } else {
        /*species 0: H */
        species[0] =
            +2.54716300e+04 * invT
            +1.96011760e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 1: H2 */
        species[1] =
            -8.35034000e+02 * invT
            +3.34653300e+00
            -2.99142300e+00 * tc[0]
            -3.50032200e-04 * tc[1]
            +9.38971500e-09 * tc[2]
            +7.69298167e-13 * tc[3]
            -7.91376000e-17 * tc[4];
        /*species 2: O */
        species[2] =
            +2.92308000e+04 * invT
            -3.37824800e+00
            -2.54206000e+00 * tc[0]
            +1.37753100e-05 * tc[1]
            +5.17133833e-10 * tc[2]
            -3.79255583e-13 * tc[3]
            +2.18402600e-17 * tc[4];
        /*species 3: OH */
        species[3] =
            +3.68362875e+03 * invT
            -3.83691187e+00
            -2.86472886e+00 * tc[0]
            -5.28252240e-04 * tc[1]
            +4.31804597e-08 * tc[2]
            -2.54348895e-12 * tc[3]
            +6.65979380e-17 * tc[4];
        /*species 4: H2O */
        species[4] =
            -2.98992100e+04 * invT
            -5.19067100e+00
            -2.67214600e+00 * tc[0]
            -1.52814650e-03 * tc[1]
            +1.45504333e-07 * tc[2]
            -1.00083000e-11 * tc[3]
            +3.19580900e-16 * tc[4];
        /*species 5: O2 */
        species[5] =
            -1.23393000e+03 * invT
            -4.91588000e-01
            -3.69757800e+00 * tc[0]
            -3.06759850e-04 * tc[1]
            +2.09807000e-08 * tc[2]
            -1.47940083e-12 * tc[3]
            +5.68217500e-17 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +1.11856713e+02 * invT
            -7.67891250e-01
            -4.01721090e+00 * tc[0]
            -1.11991006e-03 * tc[1]
            +1.05609692e-07 * tc[2]
            -9.52053083e-12 * tc[3]
            +5.39542675e-16 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            -1.80069600e+04 * invT
            +3.07203000e+00
            -4.57316700e+00 * tc[0]
            -2.16806800e-03 * tc[1]
            +2.45781500e-07 * tc[2]
            -1.95742000e-11 * tc[3]
            +7.15827000e-16 * tc[4];
        /*species 8: N2 */
        species[8] =
            -9.22797700e+02 * invT
            -4.05388800e+00
            -2.92664000e+00 * tc[0]
            -7.43988500e-04 * tc[1]
            +9.47460167e-08 * tc[2]
            -8.41420000e-12 * tc[3]
            +3.37667550e-16 * tc[4];
        /*species 9: AR */
        species[9] =
            -7.45375000e+02 * invT
            -2.86600100e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 10: HE */
        species[10] =
            -7.45375000e+02 * invT
            +5.84651100e-01
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 11: CO */
        species[11] =
            -1.42683500e+04 * invT
            -4.08314000e+00
            -3.02507800e+00 * tc[0]
            -7.21344500e-04 * tc[1]
            +9.38471333e-08 * tc[2]
            -8.48817500e-12 * tc[3]
            +3.45547600e-16 * tc[4];
        /*species 12: CO2 */
        species[12] =
            -4.89669600e+04 * invT
            +4.40901890e+00
            -4.45362300e+00 * tc[0]
            -1.57008450e-03 * tc[1]
            +2.13068500e-07 * tc[2]
            -1.99499750e-11 * tc[3]
            +8.34516500e-16 * tc[4];
    }
    return;
}


/*compute Cv/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE void cv_R(double * species, double *  tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 1: H2 */
        species[1] =
            +2.29812400e+00
            +8.24944200e-04 * tc[1]
            -8.14301500e-07 * tc[2]
            -9.47543400e-11 * tc[3]
            +4.13487200e-13 * tc[4];
        /*species 2: O */
        species[2] =
            +1.94642900e+00
            -1.63816600e-03 * tc[1]
            +2.42103200e-06 * tc[2]
            -1.60284300e-09 * tc[3]
            +3.89069600e-13 * tc[4];
        /*species 3: OH */
        species[3] =
            +3.12530561e+00
            -3.22544939e-03 * tc[1]
            +6.52764691e-06 * tc[2]
            -5.79853643e-09 * tc[3]
            +2.06237379e-12 * tc[4];
        /*species 4: H2O */
        species[4] =
            +2.38684200e+00
            +3.47498200e-03 * tc[1]
            -6.35469600e-06 * tc[2]
            +6.96858100e-09 * tc[3]
            -2.50658800e-12 * tc[4];
        /*species 5: O2 */
        species[5] =
            +2.21293600e+00
            +1.12748600e-03 * tc[1]
            -5.75615000e-07 * tc[2]
            +1.31387700e-09 * tc[3]
            -8.76855400e-13 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +3.30179801e+00
            -4.74912051e-03 * tc[1]
            +2.11582891e-05 * tc[2]
            -2.42763894e-08 * tc[3]
            +9.29225124e-12 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +2.38875400e+00
            +6.56922600e-03 * tc[1]
            -1.48501300e-07 * tc[2]
            -4.62580600e-09 * tc[3]
            +2.47151500e-12 * tc[4];
        /*species 8: N2 */
        species[8] =
            +2.29867700e+00
            +1.40824000e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485500e-12 * tc[4];
        /*species 9: AR */
        species[9] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 10: HE */
        species[10] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 11: CO */
        species[11] =
            +2.26245200e+00
            +1.51194100e-03 * tc[1]
            -3.88175500e-06 * tc[2]
            +5.58194400e-09 * tc[3]
            -2.47495100e-12 * tc[4];
        /*species 12: CO2 */
        species[12] =
            +1.27572500e+00
            +9.92207200e-03 * tc[1]
            -1.04091100e-05 * tc[2]
            +6.86668700e-09 * tc[3]
            -2.11728000e-12 * tc[4];
    } else {
        /*species 0: H */
        species[0] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 1: H2 */
        species[1] =
            +1.99142300e+00
            +7.00064400e-04 * tc[1]
            -5.63382900e-08 * tc[2]
            -9.23157800e-12 * tc[3]
            +1.58275200e-15 * tc[4];
        /*species 2: O */
        species[2] =
            +1.54206000e+00
            -2.75506200e-05 * tc[1]
            -3.10280300e-09 * tc[2]
            +4.55106700e-12 * tc[3]
            -4.36805200e-16 * tc[4];
        /*species 3: OH */
        species[3] =
            +1.86472886e+00
            +1.05650448e-03 * tc[1]
            -2.59082758e-07 * tc[2]
            +3.05218674e-11 * tc[3]
            -1.33195876e-15 * tc[4];
        /*species 4: H2O */
        species[4] =
            +1.67214600e+00
            +3.05629300e-03 * tc[1]
            -8.73026000e-07 * tc[2]
            +1.20099600e-10 * tc[3]
            -6.39161800e-15 * tc[4];
        /*species 5: O2 */
        species[5] =
            +2.69757800e+00
            +6.13519700e-04 * tc[1]
            -1.25884200e-07 * tc[2]
            +1.77528100e-11 * tc[3]
            -1.13643500e-15 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +3.01721090e+00
            +2.23982013e-03 * tc[1]
            -6.33658150e-07 * tc[2]
            +1.14246370e-10 * tc[3]
            -1.07908535e-14 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +3.57316700e+00
            +4.33613600e-03 * tc[1]
            -1.47468900e-06 * tc[2]
            +2.34890400e-10 * tc[3]
            -1.43165400e-14 * tc[4];
        /*species 8: N2 */
        species[8] =
            +1.92664000e+00
            +1.48797700e-03 * tc[1]
            -5.68476100e-07 * tc[2]
            +1.00970400e-10 * tc[3]
            -6.75335100e-15 * tc[4];
        /*species 9: AR */
        species[9] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 10: HE */
        species[10] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 11: CO */
        species[11] =
            +2.02507800e+00
            +1.44268900e-03 * tc[1]
            -5.63082800e-07 * tc[2]
            +1.01858100e-10 * tc[3]
            -6.91095200e-15 * tc[4];
        /*species 12: CO2 */
        species[12] =
            +3.45362300e+00
            +3.14016900e-03 * tc[1]
            -1.27841100e-06 * tc[2]
            +2.39399700e-10 * tc[3]
            -1.66903300e-14 * tc[4];
    }
    return;
}


/*compute Cp/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE void cp_R(double * species, double *  tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 1: H2 */
        species[1] =
            +3.29812400e+00
            +8.24944200e-04 * tc[1]
            -8.14301500e-07 * tc[2]
            -9.47543400e-11 * tc[3]
            +4.13487200e-13 * tc[4];
        /*species 2: O */
        species[2] =
            +2.94642900e+00
            -1.63816600e-03 * tc[1]
            +2.42103200e-06 * tc[2]
            -1.60284300e-09 * tc[3]
            +3.89069600e-13 * tc[4];
        /*species 3: OH */
        species[3] =
            +4.12530561e+00
            -3.22544939e-03 * tc[1]
            +6.52764691e-06 * tc[2]
            -5.79853643e-09 * tc[3]
            +2.06237379e-12 * tc[4];
        /*species 4: H2O */
        species[4] =
            +3.38684200e+00
            +3.47498200e-03 * tc[1]
            -6.35469600e-06 * tc[2]
            +6.96858100e-09 * tc[3]
            -2.50658800e-12 * tc[4];
        /*species 5: O2 */
        species[5] =
            +3.21293600e+00
            +1.12748600e-03 * tc[1]
            -5.75615000e-07 * tc[2]
            +1.31387700e-09 * tc[3]
            -8.76855400e-13 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +4.30179801e+00
            -4.74912051e-03 * tc[1]
            +2.11582891e-05 * tc[2]
            -2.42763894e-08 * tc[3]
            +9.29225124e-12 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +3.38875400e+00
            +6.56922600e-03 * tc[1]
            -1.48501300e-07 * tc[2]
            -4.62580600e-09 * tc[3]
            +2.47151500e-12 * tc[4];
        /*species 8: N2 */
        species[8] =
            +3.29867700e+00
            +1.40824000e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485500e-12 * tc[4];
        /*species 9: AR */
        species[9] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 10: HE */
        species[10] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 11: CO */
        species[11] =
            +3.26245200e+00
            +1.51194100e-03 * tc[1]
            -3.88175500e-06 * tc[2]
            +5.58194400e-09 * tc[3]
            -2.47495100e-12 * tc[4];
        /*species 12: CO2 */
        species[12] =
            +2.27572500e+00
            +9.92207200e-03 * tc[1]
            -1.04091100e-05 * tc[2]
            +6.86668700e-09 * tc[3]
            -2.11728000e-12 * tc[4];
    } else {
        /*species 0: H */
        species[0] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 1: H2 */
        species[1] =
            +2.99142300e+00
            +7.00064400e-04 * tc[1]
            -5.63382900e-08 * tc[2]
            -9.23157800e-12 * tc[3]
            +1.58275200e-15 * tc[4];
        /*species 2: O */
        species[2] =
            +2.54206000e+00
            -2.75506200e-05 * tc[1]
            -3.10280300e-09 * tc[2]
            +4.55106700e-12 * tc[3]
            -4.36805200e-16 * tc[4];
        /*species 3: OH */
        species[3] =
            +2.86472886e+00
            +1.05650448e-03 * tc[1]
            -2.59082758e-07 * tc[2]
            +3.05218674e-11 * tc[3]
            -1.33195876e-15 * tc[4];
        /*species 4: H2O */
        species[4] =
            +2.67214600e+00
            +3.05629300e-03 * tc[1]
            -8.73026000e-07 * tc[2]
            +1.20099600e-10 * tc[3]
            -6.39161800e-15 * tc[4];
        /*species 5: O2 */
        species[5] =
            +3.69757800e+00
            +6.13519700e-04 * tc[1]
            -1.25884200e-07 * tc[2]
            +1.77528100e-11 * tc[3]
            -1.13643500e-15 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +4.01721090e+00
            +2.23982013e-03 * tc[1]
            -6.33658150e-07 * tc[2]
            +1.14246370e-10 * tc[3]
            -1.07908535e-14 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +4.57316700e+00
            +4.33613600e-03 * tc[1]
            -1.47468900e-06 * tc[2]
            +2.34890400e-10 * tc[3]
            -1.43165400e-14 * tc[4];
        /*species 8: N2 */
        species[8] =
            +2.92664000e+00
            +1.48797700e-03 * tc[1]
            -5.68476100e-07 * tc[2]
            +1.00970400e-10 * tc[3]
            -6.75335100e-15 * tc[4];
        /*species 9: AR */
        species[9] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 10: HE */
        species[10] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 11: CO */
        species[11] =
            +3.02507800e+00
            +1.44268900e-03 * tc[1]
            -5.63082800e-07 * tc[2]
            +1.01858100e-10 * tc[3]
            -6.91095200e-15 * tc[4];
        /*species 12: CO2 */
        species[12] =
            +4.45362300e+00
            +3.14016900e-03 * tc[1]
            -1.27841100e-06 * tc[2]
            +2.39399700e-10 * tc[3]
            -1.66903300e-14 * tc[4];
    }
    return;
}


/*compute the e/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE void speciesInternalEnergy(double * species, double *  tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716300e+04 * invT;
        /*species 1: H2 */
        species[1] =
            +2.29812400e+00
            +4.12472100e-04 * tc[1]
            -2.71433833e-07 * tc[2]
            -2.36885850e-11 * tc[3]
            +8.26974400e-14 * tc[4]
            -1.01252100e+03 * invT;
        /*species 2: O */
        species[2] =
            +1.94642900e+00
            -8.19083000e-04 * tc[1]
            +8.07010667e-07 * tc[2]
            -4.00710750e-10 * tc[3]
            +7.78139200e-14 * tc[4]
            +2.91476400e+04 * invT;
        /*species 3: OH */
        species[3] =
            +3.12530561e+00
            -1.61272470e-03 * tc[1]
            +2.17588230e-06 * tc[2]
            -1.44963411e-09 * tc[3]
            +4.12474758e-13 * tc[4]
            +3.34630913e+03 * invT;
        /*species 4: H2O */
        species[4] =
            +2.38684200e+00
            +1.73749100e-03 * tc[1]
            -2.11823200e-06 * tc[2]
            +1.74214525e-09 * tc[3]
            -5.01317600e-13 * tc[4]
            -3.02081100e+04 * invT;
        /*species 5: O2 */
        species[5] =
            +2.21293600e+00
            +5.63743000e-04 * tc[1]
            -1.91871667e-07 * tc[2]
            +3.28469250e-10 * tc[3]
            -1.75371080e-13 * tc[4]
            -1.00524900e+03 * invT;
        /*species 6: HO2 */
        species[6] =
            +3.30179801e+00
            -2.37456025e-03 * tc[1]
            +7.05276303e-06 * tc[2]
            -6.06909735e-09 * tc[3]
            +1.85845025e-12 * tc[4]
            +2.94808040e+02 * invT;
        /*species 7: H2O2 */
        species[7] =
            +2.38875400e+00
            +3.28461300e-03 * tc[1]
            -4.95004333e-08 * tc[2]
            -1.15645150e-09 * tc[3]
            +4.94303000e-13 * tc[4]
            -1.76631500e+04 * invT;
        /*species 8: N2 */
        species[8] =
            +2.29867700e+00
            +7.04120000e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88971000e-13 * tc[4]
            -1.02090000e+03 * invT;
        /*species 9: AR */
        species[9] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
        /*species 10: HE */
        species[10] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
        /*species 11: CO */
        species[11] =
            +2.26245200e+00
            +7.55970500e-04 * tc[1]
            -1.29391833e-06 * tc[2]
            +1.39548600e-09 * tc[3]
            -4.94990200e-13 * tc[4]
            -1.43105400e+04 * invT;
        /*species 12: CO2 */
        species[12] =
            +1.27572500e+00
            +4.96103600e-03 * tc[1]
            -3.46970333e-06 * tc[2]
            +1.71667175e-09 * tc[3]
            -4.23456000e-13 * tc[4]
            -4.83731400e+04 * invT;
    } else {
        /*species 0: H */
        species[0] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716300e+04 * invT;
        /*species 1: H2 */
        species[1] =
            +1.99142300e+00
            +3.50032200e-04 * tc[1]
            -1.87794300e-08 * tc[2]
            -2.30789450e-12 * tc[3]
            +3.16550400e-16 * tc[4]
            -8.35034000e+02 * invT;
        /*species 2: O */
        species[2] =
            +1.54206000e+00
            -1.37753100e-05 * tc[1]
            -1.03426767e-09 * tc[2]
            +1.13776675e-12 * tc[3]
            -8.73610400e-17 * tc[4]
            +2.92308000e+04 * invT;
        /*species 3: OH */
        species[3] =
            +1.86472886e+00
            +5.28252240e-04 * tc[1]
            -8.63609193e-08 * tc[2]
            +7.63046685e-12 * tc[3]
            -2.66391752e-16 * tc[4]
            +3.68362875e+03 * invT;
        /*species 4: H2O */
        species[4] =
            +1.67214600e+00
            +1.52814650e-03 * tc[1]
            -2.91008667e-07 * tc[2]
            +3.00249000e-11 * tc[3]
            -1.27832360e-15 * tc[4]
            -2.98992100e+04 * invT;
        /*species 5: O2 */
        species[5] =
            +2.69757800e+00
            +3.06759850e-04 * tc[1]
            -4.19614000e-08 * tc[2]
            +4.43820250e-12 * tc[3]
            -2.27287000e-16 * tc[4]
            -1.23393000e+03 * invT;
        /*species 6: HO2 */
        species[6] =
            +3.01721090e+00
            +1.11991006e-03 * tc[1]
            -2.11219383e-07 * tc[2]
            +2.85615925e-11 * tc[3]
            -2.15817070e-15 * tc[4]
            +1.11856713e+02 * invT;
        /*species 7: H2O2 */
        species[7] =
            +3.57316700e+00
            +2.16806800e-03 * tc[1]
            -4.91563000e-07 * tc[2]
            +5.87226000e-11 * tc[3]
            -2.86330800e-15 * tc[4]
            -1.80069600e+04 * invT;
        /*species 8: N2 */
        species[8] =
            +1.92664000e+00
            +7.43988500e-04 * tc[1]
            -1.89492033e-07 * tc[2]
            +2.52426000e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
        /*species 9: AR */
        species[9] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
        /*species 10: HE */
        species[10] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
        /*species 11: CO */
        species[11] =
            +2.02507800e+00
            +7.21344500e-04 * tc[1]
            -1.87694267e-07 * tc[2]
            +2.54645250e-11 * tc[3]
            -1.38219040e-15 * tc[4]
            -1.42683500e+04 * invT;
        /*species 12: CO2 */
        species[12] =
            +3.45362300e+00
            +1.57008450e-03 * tc[1]
            -4.26137000e-07 * tc[2]
            +5.98499250e-11 * tc[3]
            -3.33806600e-15 * tc[4]
            -4.89669600e+04 * invT;
    }
    return;
}


/*compute the h/(RT) at the given temperature (Eq 20) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE void speciesEnthalpy(double * species, double *  tc)
{

    /*temperature */
    double T = tc[1];
    double invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716300e+04 * invT;
        /*species 1: H2 */
        species[1] =
            +3.29812400e+00
            +4.12472100e-04 * tc[1]
            -2.71433833e-07 * tc[2]
            -2.36885850e-11 * tc[3]
            +8.26974400e-14 * tc[4]
            -1.01252100e+03 * invT;
        /*species 2: O */
        species[2] =
            +2.94642900e+00
            -8.19083000e-04 * tc[1]
            +8.07010667e-07 * tc[2]
            -4.00710750e-10 * tc[3]
            +7.78139200e-14 * tc[4]
            +2.91476400e+04 * invT;
        /*species 3: OH */
        species[3] =
            +4.12530561e+00
            -1.61272470e-03 * tc[1]
            +2.17588230e-06 * tc[2]
            -1.44963411e-09 * tc[3]
            +4.12474758e-13 * tc[4]
            +3.34630913e+03 * invT;
        /*species 4: H2O */
        species[4] =
            +3.38684200e+00
            +1.73749100e-03 * tc[1]
            -2.11823200e-06 * tc[2]
            +1.74214525e-09 * tc[3]
            -5.01317600e-13 * tc[4]
            -3.02081100e+04 * invT;
        /*species 5: O2 */
        species[5] =
            +3.21293600e+00
            +5.63743000e-04 * tc[1]
            -1.91871667e-07 * tc[2]
            +3.28469250e-10 * tc[3]
            -1.75371080e-13 * tc[4]
            -1.00524900e+03 * invT;
        /*species 6: HO2 */
        species[6] =
            +4.30179801e+00
            -2.37456025e-03 * tc[1]
            +7.05276303e-06 * tc[2]
            -6.06909735e-09 * tc[3]
            +1.85845025e-12 * tc[4]
            +2.94808040e+02 * invT;
        /*species 7: H2O2 */
        species[7] =
            +3.38875400e+00
            +3.28461300e-03 * tc[1]
            -4.95004333e-08 * tc[2]
            -1.15645150e-09 * tc[3]
            +4.94303000e-13 * tc[4]
            -1.76631500e+04 * invT;
        /*species 8: N2 */
        species[8] =
            +3.29867700e+00
            +7.04120000e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88971000e-13 * tc[4]
            -1.02090000e+03 * invT;
        /*species 9: AR */
        species[9] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
        /*species 10: HE */
        species[10] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
        /*species 11: CO */
        species[11] =
            +3.26245200e+00
            +7.55970500e-04 * tc[1]
            -1.29391833e-06 * tc[2]
            +1.39548600e-09 * tc[3]
            -4.94990200e-13 * tc[4]
            -1.43105400e+04 * invT;
        /*species 12: CO2 */
        species[12] =
            +2.27572500e+00
            +4.96103600e-03 * tc[1]
            -3.46970333e-06 * tc[2]
            +1.71667175e-09 * tc[3]
            -4.23456000e-13 * tc[4]
            -4.83731400e+04 * invT;
    } else {
        /*species 0: H */
        species[0] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716300e+04 * invT;
        /*species 1: H2 */
        species[1] =
            +2.99142300e+00
            +3.50032200e-04 * tc[1]
            -1.87794300e-08 * tc[2]
            -2.30789450e-12 * tc[3]
            +3.16550400e-16 * tc[4]
            -8.35034000e+02 * invT;
        /*species 2: O */
        species[2] =
            +2.54206000e+00
            -1.37753100e-05 * tc[1]
            -1.03426767e-09 * tc[2]
            +1.13776675e-12 * tc[3]
            -8.73610400e-17 * tc[4]
            +2.92308000e+04 * invT;
        /*species 3: OH */
        species[3] =
            +2.86472886e+00
            +5.28252240e-04 * tc[1]
            -8.63609193e-08 * tc[2]
            +7.63046685e-12 * tc[3]
            -2.66391752e-16 * tc[4]
            +3.68362875e+03 * invT;
        /*species 4: H2O */
        species[4] =
            +2.67214600e+00
            +1.52814650e-03 * tc[1]
            -2.91008667e-07 * tc[2]
            +3.00249000e-11 * tc[3]
            -1.27832360e-15 * tc[4]
            -2.98992100e+04 * invT;
        /*species 5: O2 */
        species[5] =
            +3.69757800e+00
            +3.06759850e-04 * tc[1]
            -4.19614000e-08 * tc[2]
            +4.43820250e-12 * tc[3]
            -2.27287000e-16 * tc[4]
            -1.23393000e+03 * invT;
        /*species 6: HO2 */
        species[6] =
            +4.01721090e+00
            +1.11991006e-03 * tc[1]
            -2.11219383e-07 * tc[2]
            +2.85615925e-11 * tc[3]
            -2.15817070e-15 * tc[4]
            +1.11856713e+02 * invT;
        /*species 7: H2O2 */
        species[7] =
            +4.57316700e+00
            +2.16806800e-03 * tc[1]
            -4.91563000e-07 * tc[2]
            +5.87226000e-11 * tc[3]
            -2.86330800e-15 * tc[4]
            -1.80069600e+04 * invT;
        /*species 8: N2 */
        species[8] =
            +2.92664000e+00
            +7.43988500e-04 * tc[1]
            -1.89492033e-07 * tc[2]
            +2.52426000e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
        /*species 9: AR */
        species[9] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
        /*species 10: HE */
        species[10] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
        /*species 11: CO */
        species[11] =
            +3.02507800e+00
            +7.21344500e-04 * tc[1]
            -1.87694267e-07 * tc[2]
            +2.54645250e-11 * tc[3]
            -1.38219040e-15 * tc[4]
            -1.42683500e+04 * invT;
        /*species 12: CO2 */
        species[12] =
            +4.45362300e+00
            +1.57008450e-03 * tc[1]
            -4.26137000e-07 * tc[2]
            +5.98499250e-11 * tc[3]
            -3.33806600e-15 * tc[4]
            -4.89669600e+04 * invT;
    }
    return;
}


/*compute the S/R at the given temperature (Eq 21) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE void speciesEntropy(double * species, double *  tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H */
        species[0] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -4.60117600e-01 ;
        /*species 1: H2 */
        species[1] =
            +3.29812400e+00 * tc[0]
            +8.24944200e-04 * tc[1]
            -4.07150750e-07 * tc[2]
            -3.15847800e-11 * tc[3]
            +1.03371800e-13 * tc[4]
            -3.29409400e+00 ;
        /*species 2: O */
        species[2] =
            +2.94642900e+00 * tc[0]
            -1.63816600e-03 * tc[1]
            +1.21051600e-06 * tc[2]
            -5.34281000e-10 * tc[3]
            +9.72674000e-14 * tc[4]
            +2.96399500e+00 ;
        /*species 3: OH */
        species[3] =
            +4.12530561e+00 * tc[0]
            -3.22544939e-03 * tc[1]
            +3.26382346e-06 * tc[2]
            -1.93284548e-09 * tc[3]
            +5.15593447e-13 * tc[4]
            -6.90432960e-01 ;
        /*species 4: H2O */
        species[4] =
            +3.38684200e+00 * tc[0]
            +3.47498200e-03 * tc[1]
            -3.17734800e-06 * tc[2]
            +2.32286033e-09 * tc[3]
            -6.26647000e-13 * tc[4]
            +2.59023300e+00 ;
        /*species 5: O2 */
        species[5] =
            +3.21293600e+00 * tc[0]
            +1.12748600e-03 * tc[1]
            -2.87807500e-07 * tc[2]
            +4.37959000e-10 * tc[3]
            -2.19213850e-13 * tc[4]
            +6.03473800e+00 ;
        /*species 6: HO2 */
        species[6] =
            +4.30179801e+00 * tc[0]
            -4.74912051e-03 * tc[1]
            +1.05791445e-05 * tc[2]
            -8.09212980e-09 * tc[3]
            +2.32306281e-12 * tc[4]
            +3.71666245e+00 ;
        /*species 7: H2O2 */
        species[7] =
            +3.38875400e+00 * tc[0]
            +6.56922600e-03 * tc[1]
            -7.42506500e-08 * tc[2]
            -1.54193533e-09 * tc[3]
            +6.17878750e-13 * tc[4]
            +6.78536300e+00 ;
        /*species 8: N2 */
        species[8] =
            +3.29867700e+00 * tc[0]
            +1.40824000e-03 * tc[1]
            -1.98161100e-06 * tc[2]
            +1.88050500e-09 * tc[3]
            -6.11213750e-13 * tc[4]
            +3.95037200e+00 ;
        /*species 9: AR */
        species[9] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +4.36600100e+00 ;
        /*species 10: HE */
        species[10] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +9.15348800e-01 ;
        /*species 11: CO */
        species[11] =
            +3.26245200e+00 * tc[0]
            +1.51194100e-03 * tc[1]
            -1.94087750e-06 * tc[2]
            +1.86064800e-09 * tc[3]
            -6.18737750e-13 * tc[4]
            +4.84889700e+00 ;
        /*species 12: CO2 */
        species[12] =
            +2.27572500e+00 * tc[0]
            +9.92207200e-03 * tc[1]
            -5.20455500e-06 * tc[2]
            +2.28889567e-09 * tc[3]
            -5.29320000e-13 * tc[4]
            +1.01884900e+01 ;
    } else {
        /*species 0: H */
        species[0] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -4.60117600e-01 ;
        /*species 1: H2 */
        species[1] =
            +2.99142300e+00 * tc[0]
            +7.00064400e-04 * tc[1]
            -2.81691450e-08 * tc[2]
            -3.07719267e-12 * tc[3]
            +3.95688000e-16 * tc[4]
            -1.35511000e+00 ;
        /*species 2: O */
        species[2] =
            +2.54206000e+00 * tc[0]
            -2.75506200e-05 * tc[1]
            -1.55140150e-09 * tc[2]
            +1.51702233e-12 * tc[3]
            -1.09201300e-16 * tc[4]
            +4.92030800e+00 ;
        /*species 3: OH */
        species[3] =
            +2.86472886e+00 * tc[0]
            +1.05650448e-03 * tc[1]
            -1.29541379e-07 * tc[2]
            +1.01739558e-11 * tc[3]
            -3.32989690e-16 * tc[4]
            +5.70164073e+00 ;
        /*species 4: H2O */
        species[4] =
            +2.67214600e+00 * tc[0]
            +3.05629300e-03 * tc[1]
            -4.36513000e-07 * tc[2]
            +4.00332000e-11 * tc[3]
            -1.59790450e-15 * tc[4]
            +6.86281700e+00 ;
        /*species 5: O2 */
        species[5] =
            +3.69757800e+00 * tc[0]
            +6.13519700e-04 * tc[1]
            -6.29421000e-08 * tc[2]
            +5.91760333e-12 * tc[3]
            -2.84108750e-16 * tc[4]
            +3.18916600e+00 ;
        /*species 6: HO2 */
        species[6] =
            +4.01721090e+00 * tc[0]
            +2.23982013e-03 * tc[1]
            -3.16829075e-07 * tc[2]
            +3.80821233e-11 * tc[3]
            -2.69771337e-15 * tc[4]
            +3.78510215e+00 ;
        /*species 7: H2O2 */
        species[7] =
            +4.57316700e+00 * tc[0]
            +4.33613600e-03 * tc[1]
            -7.37344500e-07 * tc[2]
            +7.82968000e-11 * tc[3]
            -3.57913500e-15 * tc[4]
            +5.01137000e-01 ;
        /*species 8: N2 */
        species[8] =
            +2.92664000e+00 * tc[0]
            +1.48797700e-03 * tc[1]
            -2.84238050e-07 * tc[2]
            +3.36568000e-11 * tc[3]
            -1.68833775e-15 * tc[4]
            +5.98052800e+00 ;
        /*species 9: AR */
        species[9] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +4.36600100e+00 ;
        /*species 10: HE */
        species[10] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +9.15348900e-01 ;
        /*species 11: CO */
        species[11] =
            +3.02507800e+00 * tc[0]
            +1.44268900e-03 * tc[1]
            -2.81541400e-07 * tc[2]
            +3.39527000e-11 * tc[3]
            -1.72773800e-15 * tc[4]
            +6.10821800e+00 ;
        /*species 12: CO2 */
        species[12] =
            +4.45362300e+00 * tc[0]
            +3.14016900e-03 * tc[1]
            -6.39205500e-07 * tc[2]
            +7.97999000e-11 * tc[3]
            -4.17258250e-15 * tc[4]
            -9.55395900e-01 ;
    }
    return;
}


/*save atomic weights into array */
void atomicWeight(double *  awt)
{
    awt[0] = 1.007970; /*H */
    awt[1] = 15.999400; /*O */
    awt[2] = 14.006700; /*N */
    awt[3] = 39.948000; /*AR */
    awt[4] = 4.002600; /*HE */
    awt[5] = 12.011150; /*C */

    return;
}


/* get temperature given internal energy in mass units and mass fracs */
AMREX_GPU_HOST_DEVICE void GET_T_GIVEN_EY(double *  e, double *  y, double *  t, int * ierr)
{
#ifdef CONVERGENCE
    const int maxiter = 5000;
    const double tol  = 1.e-12;
#else
    const int maxiter = 200;
    const double tol  = 1.e-6;
#endif
    double ein  = *e;
    double tmin = 90;/*max lower bound for thermo def */
    double tmax = 4000;/*min upper bound for thermo def */
    double e1,emin,emax,cv,t1,dt;
    int i;/* loop counter */
    CKUBMS(&tmin, y, &emin);
    CKUBMS(&tmax, y, &emax);
    if (ein < emin) {
        /*Linear Extrapolation below tmin */
        CKCVBS(&tmin, y, &cv);
        *t = tmin - (emin-ein)/cv;
        *ierr = 1;
        return;
    }
    if (ein > emax) {
        /*Linear Extrapolation above tmax */
        CKCVBS(&tmax, y, &cv);
        *t = tmax - (emax-ein)/cv;
        *ierr = 1;
        return;
    }
    t1 = *t;
    if (t1 < tmin || t1 > tmax) {
        t1 = tmin + (tmax-tmin)/(emax-emin)*(ein-emin);
    }
    for (i = 0; i < maxiter; ++i) {
        CKUBMS(&t1,y,&e1);
        CKCVBS(&t1,y,&cv);
        dt = (ein - e1) / cv;
        if (dt > 100.) { dt = 100.; }
        else if (dt < -100.) { dt = -100.; }
        else if (fabs(dt) < tol) break;
        else if (t1+dt == t1) break;
        t1 += dt;
    }
    *t = t1;
    *ierr = 0;
    return;
}

/* get temperature given enthalpy in mass units and mass fracs */
AMREX_GPU_HOST_DEVICE void GET_T_GIVEN_HY(double *  h, double *  y, double *  t, int * ierr)
{
#ifdef CONVERGENCE
    const int maxiter = 5000;
    const double tol  = 1.e-12;
#else
    const int maxiter = 200;
    const double tol  = 1.e-6;
#endif
    double hin  = *h;
    double tmin = 90;/*max lower bound for thermo def */
    double tmax = 4000;/*min upper bound for thermo def */
    double h1,hmin,hmax,cp,t1,dt;
    int i;/* loop counter */
    CKHBMS(&tmin, y, &hmin);
    CKHBMS(&tmax, y, &hmax);
    if (hin < hmin) {
        /*Linear Extrapolation below tmin */
        CKCPBS(&tmin, y, &cp);
        *t = tmin - (hmin-hin)/cp;
        *ierr = 1;
        return;
    }
    if (hin > hmax) {
        /*Linear Extrapolation above tmax */
        CKCPBS(&tmax, y, &cp);
        *t = tmax - (hmax-hin)/cp;
        *ierr = 1;
        return;
    }
    t1 = *t;
    if (t1 < tmin || t1 > tmax) {
        t1 = tmin + (tmax-tmin)/(hmax-hmin)*(hin-hmin);
    }
    for (i = 0; i < maxiter; ++i) {
        CKHBMS(&t1,y,&h1);
        CKCPBS(&t1,y,&cp);
        dt = (hin - h1) / cp;
        if (dt > 100.) { dt = 100.; }
        else if (dt < -100.) { dt = -100.; }
        else if (fabs(dt) < tol) break;
        else if (t1+dt == t1) break;
        t1 += dt;
    }
    *t = t1;
    *ierr = 0;
    return;
}


/*compute the critical parameters for each species */
void GET_CRITPARAMS(double *  Tci, double *  ai, double *  bi, double *  acentric_i)
{

    double   EPS[13];
    double   SIG[13];
    double    wt[13];
    double avogadro = 6.02214199e23;
    double boltzmann = 1.3806503e-16; //we work in CGS
    double Rcst = 83.144598; //in bar [CGS] !

    egtransetEPS(EPS);
    egtransetSIG(SIG);
    get_mw(wt);

    /*species 0: H */
    Tci[0] = 1.316 * EPS[0] ; 
    ai[0] = (5.55 * pow(avogadro,2.0) * EPS[0]*boltzmann * pow(1e-8*SIG[0],3.0) ) / (pow(wt[0],2.0)); 
    bi[0] = 0.855 * avogadro * pow(1e-8*SIG[0],3.0) / (wt[0]); 
    acentric_i[0] = 0.0 ;

    /*species 1: H2 */
    /*Imported from NIST */
    Tci[1] = 33.145000 ; 
    ai[1] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[1],2.0) / (pow(2.015880,2.0) * 12.964000); 
    bi[1] = 0.08664 * Rcst * Tci[1] / (2.015880 * 12.964000); 
    acentric_i[1] = -0.219000 ;

    /*species 2: O */
    Tci[2] = 1.316 * EPS[2] ; 
    ai[2] = (5.55 * pow(avogadro,2.0) * EPS[2]*boltzmann * pow(1e-8*SIG[2],3.0) ) / (pow(wt[2],2.0)); 
    bi[2] = 0.855 * avogadro * pow(1e-8*SIG[2],3.0) / (wt[2]); 
    acentric_i[2] = 0.0 ;

    /*species 3: OH */
    Tci[3] = 1.316 * EPS[3] ; 
    ai[3] = (5.55 * pow(avogadro,2.0) * EPS[3]*boltzmann * pow(1e-8*SIG[3],3.0) ) / (pow(wt[3],2.0)); 
    bi[3] = 0.855 * avogadro * pow(1e-8*SIG[3],3.0) / (wt[3]); 
    acentric_i[3] = 0.0 ;

    /*species 4: H2O */
    /*Imported from NIST */
    Tci[4] = 647.096000 ; 
    ai[4] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[4],2.0) / (pow(18.015340,2.0) * 220.640000); 
    bi[4] = 0.08664 * Rcst * Tci[4] / (18.015340 * 220.640000); 
    acentric_i[4] = 0.344300 ;

    /*species 5: O2 */
    /*Imported from NIST */
    Tci[5] = 154.581000 ; 
    ai[5] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[5],2.0) / (pow(31.998800,2.0) * 50.430466); 
    bi[5] = 0.08664 * Rcst * Tci[5] / (31.998800 * 50.430466); 
    acentric_i[5] = 0.022200 ;

    /*species 6: HO2 */
    Tci[6] = 1.316 * EPS[6] ; 
    ai[6] = (5.55 * pow(avogadro,2.0) * EPS[6]*boltzmann * pow(1e-8*SIG[6],3.0) ) / (pow(wt[6],2.0)); 
    bi[6] = 0.855 * avogadro * pow(1e-8*SIG[6],3.0) / (wt[6]); 
    acentric_i[6] = 0.0 ;

    /*species 7: H2O2 */
    Tci[7] = 1.316 * EPS[7] ; 
    ai[7] = (5.55 * pow(avogadro,2.0) * EPS[7]*boltzmann * pow(1e-8*SIG[7],3.0) ) / (pow(wt[7],2.0)); 
    bi[7] = 0.855 * avogadro * pow(1e-8*SIG[7],3.0) / (wt[7]); 
    acentric_i[7] = 0.0 ;

    /*species 8: N2 */
    /*Imported from NIST */
    Tci[8] = 126.192000 ; 
    ai[8] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[8],2.0) / (pow(28.013400,2.0) * 33.958000); 
    bi[8] = 0.08664 * Rcst * Tci[8] / (28.013400 * 33.958000); 
    acentric_i[8] = 0.037200 ;

    /*species 9: AR */
    /*Imported from NIST */
    Tci[9] = 150.860000 ; 
    ai[9] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[9],2.0) / (pow(39.948000,2.0) * 48.980000); 
    bi[9] = 0.08664 * Rcst * Tci[9] / (39.948000 * 48.980000); 
    acentric_i[9] = -0.002000 ;

    /*species 10: HE */
    Tci[10] = 1.316 * EPS[10] ; 
    ai[10] = (5.55 * pow(avogadro,2.0) * EPS[10]*boltzmann * pow(1e-8*SIG[10],3.0) ) / (pow(wt[10],2.0)); 
    bi[10] = 0.855 * avogadro * pow(1e-8*SIG[10],3.0) / (wt[10]); 
    acentric_i[10] = 0.0 ;

    /*species 11: CO */
    /*Imported from NIST */
    Tci[11] = 132.850000 ; 
    ai[11] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[11],2.0) / (pow(28.010000,2.0) * 34.940000); 
    bi[11] = 0.08664 * Rcst * Tci[11] / (28.010000 * 34.940000); 
    acentric_i[11] = 0.045000 ;

    /*species 12: CO2 */
    /*Imported from NIST */
    Tci[12] = 304.120000 ; 
    ai[12] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[12],2.0) / (pow(44.009950,2.0) * 73.740000); 
    bi[12] = 0.08664 * Rcst * Tci[12] / (44.009950 * 73.740000); 
    acentric_i[12] = 0.225000 ;

    return;
}


void egtransetLENIMC(int* LENIMC ) {
    *LENIMC = 55;}


void egtransetLENRMC(int* LENRMC ) {
    *LENRMC = 3718;}


void egtransetNO(int* NO ) {
    *NO = 4;}


void egtransetKK(int* KK ) {
    *KK = 13;}


void egtransetNLITE(int* NLITE ) {
    *NLITE = 3;}


/*Patm in ergs/cm3 */
void egtransetPATM(double* PATM) {
    *PATM =   0.1013250000000000E+07;}


/*the molecular weights in g/mol */
void egtransetWT(double* WT ) {
    WT[0] = 1.00797000E+00;
    WT[1] = 2.01594000E+00;
    WT[2] = 1.59994000E+01;
    WT[3] = 1.70073700E+01;
    WT[4] = 1.80153400E+01;
    WT[5] = 3.19988000E+01;
    WT[6] = 3.30067700E+01;
    WT[7] = 3.40147400E+01;
    WT[8] = 2.80134000E+01;
    WT[9] = 3.99480000E+01;
    WT[10] = 4.00260000E+00;
    WT[11] = 2.80105500E+01;
    WT[12] = 4.40099500E+01;
}


/*the lennard-jones potential well depth eps/kb in K */
void egtransetEPS(double* EPS ) {
    EPS[0] = 1.45000000E+02;
    EPS[1] = 3.80000000E+01;
    EPS[2] = 8.00000000E+01;
    EPS[3] = 8.00000000E+01;
    EPS[4] = 5.72400000E+02;
    EPS[5] = 1.07400000E+02;
    EPS[6] = 1.07400000E+02;
    EPS[7] = 1.07400000E+02;
    EPS[8] = 9.75300000E+01;
    EPS[9] = 1.36500000E+02;
    EPS[10] = 1.02000000E+01;
    EPS[11] = 9.81000000E+01;
    EPS[12] = 2.44000000E+02;
}


/*the lennard-jones collision diameter in Angstroms */
void egtransetSIG(double* SIG ) {
    SIG[0] = 2.05000000E+00;
    SIG[1] = 2.92000000E+00;
    SIG[2] = 2.75000000E+00;
    SIG[3] = 2.75000000E+00;
    SIG[4] = 2.60500000E+00;
    SIG[5] = 3.45800000E+00;
    SIG[6] = 3.45800000E+00;
    SIG[7] = 3.45800000E+00;
    SIG[8] = 3.62100000E+00;
    SIG[9] = 3.33000000E+00;
    SIG[10] = 2.57600000E+00;
    SIG[11] = 3.65000000E+00;
    SIG[12] = 3.76300000E+00;
}


/*the dipole moment in Debye */
void egtransetDIP(double* DIP ) {
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
    DIP[12] = 0.00000000E+00;
}


/*the polarizability in cubic Angstroms */
void egtransetPOL(double* POL ) {
    POL[0] = 0.00000000E+00;
    POL[1] = 7.90000000E-01;
    POL[2] = 0.00000000E+00;
    POL[3] = 0.00000000E+00;
    POL[4] = 0.00000000E+00;
    POL[5] = 1.60000000E+00;
    POL[6] = 0.00000000E+00;
    POL[7] = 0.00000000E+00;
    POL[8] = 1.76000000E+00;
    POL[9] = 0.00000000E+00;
    POL[10] = 0.00000000E+00;
    POL[11] = 1.95000000E+00;
    POL[12] = 2.65000000E+00;
}


/*the rotational relaxation collision number at 298 K */
void egtransetZROT(double* ZROT ) {
    ZROT[0] = 0.00000000E+00;
    ZROT[1] = 2.80000000E+02;
    ZROT[2] = 0.00000000E+00;
    ZROT[3] = 0.00000000E+00;
    ZROT[4] = 4.00000000E+00;
    ZROT[5] = 3.80000000E+00;
    ZROT[6] = 1.00000000E+00;
    ZROT[7] = 3.80000000E+00;
    ZROT[8] = 4.00000000E+00;
    ZROT[9] = 0.00000000E+00;
    ZROT[10] = 0.00000000E+00;
    ZROT[11] = 1.80000000E+00;
    ZROT[12] = 2.10000000E+00;
}


/*0: monoatomic, 1: linear, 2: nonlinear */
void egtransetNLIN(int* NLIN) {
    NLIN[0] = 0;
    NLIN[1] = 1;
    NLIN[2] = 0;
    NLIN[3] = 1;
    NLIN[4] = 2;
    NLIN[5] = 1;
    NLIN[6] = 2;
    NLIN[7] = 2;
    NLIN[8] = 1;
    NLIN[9] = 0;
    NLIN[10] = 0;
    NLIN[11] = 1;
    NLIN[12] = 1;
}


/*Poly fits for the viscosities, dim NO*KK */
void egtransetCOFETA(double* COFETA) {
    COFETA[0] = -1.98744496E+01;
    COFETA[1] = 3.41660514E+00;
    COFETA[2] = -3.63206306E-01;
    COFETA[3] = 1.58671021E-02;
    COFETA[4] = -1.37549435E+01;
    COFETA[5] = 9.65530587E-01;
    COFETA[6] = -4.45720114E-02;
    COFETA[7] = 2.05871810E-03;
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
    COFETA[20] = -1.68118868E+01;
    COFETA[21] = 2.52362554E+00;
    COFETA[22] = -2.49309128E-01;
    COFETA[23] = 1.10211025E-02;
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
    COFETA[36] = -1.86067598E+01;
    COFETA[37] = 3.27402596E+00;
    COFETA[38] = -3.45827972E-01;
    COFETA[39] = 1.51622680E-02;
    COFETA[40] = -1.11555213E+01;
    COFETA[41] = 2.18772782E-01;
    COFETA[42] = 5.60263799E-02;
    COFETA[43] = -2.36018246E-03;
    COFETA[44] = -1.63031240E+01;
    COFETA[45] = 2.26143219E+00;
    COFETA[46] = -2.15114671E-01;
    COFETA[47] = 9.53461976E-03;
    COFETA[48] = -2.36749526E+01;
    COFETA[49] = 4.99775518E+00;
    COFETA[50] = -5.52687718E-01;
    COFETA[51] = 2.34353338E-02;
}


/*Poly fits for the conductivities, dim NO*KK */
void egtransetCOFLAM(double* COFLAM) {
    COFLAM[0] = -3.24539191E-01;
    COFLAM[1] = 3.41660514E+00;
    COFLAM[2] = -3.63206306E-01;
    COFLAM[3] = 1.58671021E-02;
    COFLAM[4] = 1.11035511E+01;
    COFLAM[5] = -1.31883912E+00;
    COFLAM[6] = 2.44042473E-01;
    COFLAM[7] = -8.99836359E-03;
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
    COFLAM[20] = -2.51296685E+00;
    COFLAM[21] = 3.15165687E+00;
    COFLAM[22] = -3.10007661E-01;
    COFLAM[23] = 1.34522321E-02;
    COFLAM[24] = 5.56023763E-01;
    COFLAM[25] = 1.59073590E+00;
    COFLAM[26] = -5.28053839E-02;
    COFLAM[27] = 4.07601571E-04;
    COFLAM[28] = 1.48801095E+00;
    COFLAM[29] = 1.06176238E+00;
    COFLAM[30] = 5.72195640E-02;
    COFLAM[31] = -6.38391491E-03;
    COFLAM[32] = 1.15507419E+01;
    COFLAM[33] = -2.91453917E+00;
    COFLAM[34] = 5.55045765E-01;
    COFLAM[35] = -2.75173485E-02;
    COFLAM[36] = -2.73648952E+00;
    COFLAM[37] = 3.27402596E+00;
    COFLAM[38] = -3.45827972E-01;
    COFLAM[39] = 1.51622680E-02;
    COFLAM[40] = 7.01538340E+00;
    COFLAM[41] = 2.18772782E-01;
    COFLAM[42] = 5.60263799E-02;
    COFLAM[43] = -2.36018246E-03;
    COFLAM[44] = 9.94279402E+00;
    COFLAM[45] = -2.29161875E+00;
    COFLAM[46] = 4.74393585E-01;
    COFLAM[47] = -2.40686534E-02;
    COFLAM[48] = -1.21375633E+01;
    COFLAM[49] = 6.23624427E+00;
    COFLAM[50] = -6.22471402E-01;
    COFLAM[51] = 2.30613281E-02;
}


/*Poly fits for the diffusion coefficients, dim NO*KK*KK */
void egtransetCOFD(double* COFD) {
    COFD[0] = -1.43693056E+01;
    COFD[1] = 4.03992999E+00;
    COFD[2] = -3.08044800E-01;
    COFD[3] = 1.32757775E-02;
    COFD[4] = -1.11808682E+01;
    COFD[5] = 2.66936727E+00;
    COFD[6] = -1.34411514E-01;
    COFD[7] = 5.92957488E-03;
    COFD[8] = -1.31860117E+01;
    COFD[9] = 3.38003453E+00;
    COFD[10] = -2.25783856E-01;
    COFD[11] = 9.85028660E-03;
    COFD[12] = -1.31877711E+01;
    COFD[13] = 3.38003453E+00;
    COFD[14] = -2.25783856E-01;
    COFD[15] = 9.85028660E-03;
    COFD[16] = -1.93611051E+01;
    COFD[17] = 5.51579726E+00;
    COFD[18] = -4.76061961E-01;
    COFD[19] = 1.96329391E-02;
    COFD[20] = -1.43712864E+01;
    COFD[21] = 3.70920439E+00;
    COFD[22] = -2.67274113E-01;
    COFD[23] = 1.15967481E-02;
    COFD[24] = -1.43717529E+01;
    COFD[25] = 3.70920439E+00;
    COFD[26] = -2.67274113E-01;
    COFD[27] = 1.15967481E-02;
    COFD[28] = -1.43721922E+01;
    COFD[29] = 3.70920439E+00;
    COFD[30] = -2.67274113E-01;
    COFD[31] = 1.15967481E-02;
    COFD[32] = -1.40298830E+01;
    COFD[33] = 3.55837688E+00;
    COFD[34] = -2.47785790E-01;
    COFD[35] = 1.07555332E-02;
    COFD[36] = -1.51208119E+01;
    COFD[37] = 3.99904647E+00;
    COFD[38] = -3.03517220E-01;
    COFD[39] = 1.31117363E-02;
    COFD[40] = -9.71338331E+00;
    COFD[41] = 2.17561180E+00;
    COFD[42] = -7.28270090E-02;
    COFD[43] = 3.38302182E-03;
    COFD[44] = -1.40524065E+01;
    COFD[45] = 3.56261348E+00;
    COFD[46] = -2.48287981E-01;
    COFD[47] = 1.07752947E-02;
    COFD[48] = -1.72993972E+01;
    COFD[49] = 4.71931868E+00;
    COFD[50] = -3.91258152E-01;
    COFD[51] = 1.66866639E-02;
    COFD[52] = -1.11808682E+01;
    COFD[53] = 2.66936727E+00;
    COFD[54] = -1.34411514E-01;
    COFD[55] = 5.92957488E-03;
    COFD[56] = -1.02395222E+01;
    COFD[57] = 2.15403244E+00;
    COFD[58] = -6.97480266E-02;
    COFD[59] = 3.23666871E-03;
    COFD[60] = -1.06250182E+01;
    COFD[61] = 2.15849701E+00;
    COFD[62] = -6.53886401E-02;
    COFD[63] = 2.81453370E-03;
    COFD[64] = -1.06283453E+01;
    COFD[65] = 2.15849701E+00;
    COFD[66] = -6.53886401E-02;
    COFD[67] = 2.81453370E-03;
    COFD[68] = -1.68758926E+01;
    COFD[69] = 4.49460303E+00;
    COFD[70] = -3.64766132E-01;
    COFD[71] = 1.56457153E-02;
    COFD[72] = -1.15797750E+01;
    COFD[73] = 2.43235504E+00;
    COFD[74] = -1.02890179E-01;
    COFD[75] = 4.52903603E-03;
    COFD[76] = -1.15806808E+01;
    COFD[77] = 2.43235504E+00;
    COFD[78] = -1.02890179E-01;
    COFD[79] = 4.52903603E-03;
    COFD[80] = -1.15815344E+01;
    COFD[81] = 2.43235504E+00;
    COFD[82] = -1.02890179E-01;
    COFD[83] = 4.52903603E-03;
    COFD[84] = -1.13253458E+01;
    COFD[85] = 2.31195095E+00;
    COFD[86] = -8.63988037E-02;
    COFD[87] = 3.77573452E-03;
    COFD[88] = -1.20638601E+01;
    COFD[89] = 2.63303536E+00;
    COFD[90] = -1.29792632E-01;
    COFD[91] = 5.73363738E-03;
    COFD[92] = -9.86429034E+00;
    COFD[93] = 2.05348746E+00;
    COFD[94] = -5.90289007E-02;
    COFD[95] = 2.89596157E-03;
    COFD[96] = -1.13541075E+01;
    COFD[97] = 2.31999438E+00;
    COFD[98] = -8.75064804E-02;
    COFD[99] = 3.82656365E-03;
    COFD[100] = -1.35545239E+01;
    COFD[101] = 3.13878730E+00;
    COFD[102] = -1.94980335E-01;
    COFD[103] = 8.53744486E-03;
    COFD[104] = -1.31860117E+01;
    COFD[105] = 3.38003453E+00;
    COFD[106] = -2.25783856E-01;
    COFD[107] = 9.85028660E-03;
    COFD[108] = -1.06250182E+01;
    COFD[109] = 2.15849701E+00;
    COFD[110] = -6.53886401E-02;
    COFD[111] = 2.81453370E-03;
    COFD[112] = -1.29877365E+01;
    COFD[113] = 2.80841511E+00;
    COFD[114] = -1.52629888E-01;
    COFD[115] = 6.72604927E-03;
    COFD[116] = -1.30027772E+01;
    COFD[117] = 2.80841511E+00;
    COFD[118] = -1.52629888E-01;
    COFD[119] = 6.72604927E-03;
    COFD[120] = -1.91096797E+01;
    COFD[121] = 5.02608697E+00;
    COFD[122] = -4.26959993E-01;
    COFD[123] = 1.80709910E-02;
    COFD[124] = -1.40864894E+01;
    COFD[125] = 3.07458927E+00;
    COFD[126] = -1.86899591E-01;
    COFD[127] = 8.19829781E-03;
    COFD[128] = -1.40916052E+01;
    COFD[129] = 3.07458927E+00;
    COFD[130] = -1.86899591E-01;
    COFD[131] = 8.19829781E-03;
    COFD[132] = -1.40964661E+01;
    COFD[133] = 3.07458927E+00;
    COFD[134] = -1.86899591E-01;
    COFD[135] = 8.19829781E-03;
    COFD[136] = -1.38756407E+01;
    COFD[137] = 2.98558426E+00;
    COFD[138] = -1.75507216E-01;
    COFD[139] = 7.71173691E-03;
    COFD[140] = -1.47082523E+01;
    COFD[141] = 3.30683499E+00;
    COFD[142] = -2.16378602E-01;
    COFD[143] = 9.44670561E-03;
    COFD[144] = -9.70779324E+00;
    COFD[145] = 1.77912272E+00;
    COFD[146] = -1.67349571E-02;
    COFD[147] = 7.45446845E-04;
    COFD[148] = -1.39007410E+01;
    COFD[149] = 2.99164244E+00;
    COFD[150] = -1.76293106E-01;
    COFD[151] = 7.74575100E-03;
    COFD[152] = -1.67115577E+01;
    COFD[153] = 3.98859394E+00;
    COFD[154] = -3.02316219E-01;
    COFD[155] = 1.30661099E-02;
    COFD[156] = -1.31877711E+01;
    COFD[157] = 3.38003453E+00;
    COFD[158] = -2.25783856E-01;
    COFD[159] = 9.85028660E-03;
    COFD[160] = -1.06283453E+01;
    COFD[161] = 2.15849701E+00;
    COFD[162] = -6.53886401E-02;
    COFD[163] = 2.81453370E-03;
    COFD[164] = -1.30027772E+01;
    COFD[165] = 2.80841511E+00;
    COFD[166] = -1.52629888E-01;
    COFD[167] = 6.72604927E-03;
    COFD[168] = -1.30182843E+01;
    COFD[169] = 2.80841511E+00;
    COFD[170] = -1.52629888E-01;
    COFD[171] = 6.72604927E-03;
    COFD[172] = -1.91256261E+01;
    COFD[173] = 5.02608697E+00;
    COFD[174] = -4.26959993E-01;
    COFD[175] = 1.80709910E-02;
    COFD[176] = -1.41066459E+01;
    COFD[177] = 3.07458927E+00;
    COFD[178] = -1.86899591E-01;
    COFD[179] = 8.19829781E-03;
    COFD[180] = -1.41119732E+01;
    COFD[181] = 3.07458927E+00;
    COFD[182] = -1.86899591E-01;
    COFD[183] = 8.19829781E-03;
    COFD[184] = -1.41170372E+01;
    COFD[185] = 3.07458927E+00;
    COFD[186] = -1.86899591E-01;
    COFD[187] = 8.19829781E-03;
    COFD[188] = -1.38948667E+01;
    COFD[189] = 2.98558426E+00;
    COFD[190] = -1.75507216E-01;
    COFD[191] = 7.71173691E-03;
    COFD[192] = -1.47298720E+01;
    COFD[193] = 3.30683499E+00;
    COFD[194] = -2.16378602E-01;
    COFD[195] = 9.44670561E-03;
    COFD[196] = -9.71375861E+00;
    COFD[197] = 1.77912272E+00;
    COFD[198] = -1.67349571E-02;
    COFD[199] = 7.45446845E-04;
    COFD[200] = -1.39199663E+01;
    COFD[201] = 2.99164244E+00;
    COFD[202] = -1.76293106E-01;
    COFD[203] = 7.74575100E-03;
    COFD[204] = -1.67337768E+01;
    COFD[205] = 3.98859394E+00;
    COFD[206] = -3.02316219E-01;
    COFD[207] = 1.30661099E-02;
    COFD[208] = -1.93611051E+01;
    COFD[209] = 5.51579726E+00;
    COFD[210] = -4.76061961E-01;
    COFD[211] = 1.96329391E-02;
    COFD[212] = -1.68758926E+01;
    COFD[213] = 4.49460303E+00;
    COFD[214] = -3.64766132E-01;
    COFD[215] = 1.56457153E-02;
    COFD[216] = -1.91096797E+01;
    COFD[217] = 5.02608697E+00;
    COFD[218] = -4.26959993E-01;
    COFD[219] = 1.80709910E-02;
    COFD[220] = -1.91256261E+01;
    COFD[221] = 5.02608697E+00;
    COFD[222] = -4.26959993E-01;
    COFD[223] = 1.80709910E-02;
    COFD[224] = -1.31492641E+01;
    COFD[225] = 1.48004311E+00;
    COFD[226] = 1.60499553E-01;
    COFD[227] = -1.19765679E-02;
    COFD[228] = -2.10640014E+01;
    COFD[229] = 5.50980695E+00;
    COFD[230] = -4.78335488E-01;
    COFD[231] = 1.98515434E-02;
    COFD[232] = -2.04177482E+01;
    COFD[233] = 5.31457079E+00;
    COFD[234] = -4.58216496E-01;
    COFD[235] = 1.91825910E-02;
    COFD[236] = -2.04230073E+01;
    COFD[237] = 5.31457079E+00;
    COFD[238] = -4.58216496E-01;
    COFD[239] = 1.91825910E-02;
    COFD[240] = -2.08123325E+01;
    COFD[241] = 5.42470154E+00;
    COFD[242] = -4.69700416E-01;
    COFD[243] = 1.95706904E-02;
    COFD[244] = -2.10785324E+01;
    COFD[245] = 5.51573149E+00;
    COFD[246] = -4.78177665E-01;
    COFD[247] = 1.98082796E-02;
    COFD[248] = -1.21950642E+01;
    COFD[249] = 2.72222246E+00;
    COFD[250] = -1.41335602E-01;
    COFD[251] = 6.23222872E-03;
    COFD[252] = -2.08943798E+01;
    COFD[253] = 5.44718652E+00;
    COFD[254] = -4.72082953E-01;
    COFD[255] = 1.96531321E-02;
    COFD[256] = -2.12021508E+01;
    COFD[257] = 5.20775052E+00;
    COFD[258] = -4.07348327E-01;
    COFD[259] = 1.55473283E-02;
    COFD[260] = -1.43712864E+01;
    COFD[261] = 3.70920439E+00;
    COFD[262] = -2.67274113E-01;
    COFD[263] = 1.15967481E-02;
    COFD[264] = -1.15797750E+01;
    COFD[265] = 2.43235504E+00;
    COFD[266] = -1.02890179E-01;
    COFD[267] = 4.52903603E-03;
    COFD[268] = -1.40864894E+01;
    COFD[269] = 3.07458927E+00;
    COFD[270] = -1.86899591E-01;
    COFD[271] = 8.19829781E-03;
    COFD[272] = -1.41066459E+01;
    COFD[273] = 3.07458927E+00;
    COFD[274] = -1.86899591E-01;
    COFD[275] = 8.19829781E-03;
    COFD[276] = -2.10640014E+01;
    COFD[277] = 5.50980695E+00;
    COFD[278] = -4.78335488E-01;
    COFD[279] = 1.98515434E-02;
    COFD[280] = -1.53110708E+01;
    COFD[281] = 3.37317428E+00;
    COFD[282] = -2.24900439E-01;
    COFD[283] = 9.81228151E-03;
    COFD[284] = -1.53187643E+01;
    COFD[285] = 3.37317428E+00;
    COFD[286] = -2.24900439E-01;
    COFD[287] = 9.81228151E-03;
    COFD[288] = -1.53261114E+01;
    COFD[289] = 3.37317428E+00;
    COFD[290] = -2.24900439E-01;
    COFD[291] = 9.81228151E-03;
    COFD[292] = -1.50096240E+01;
    COFD[293] = 3.25515933E+00;
    COFD[294] = -2.09710110E-01;
    COFD[295] = 9.15941830E-03;
    COFD[296] = -1.59592184E+01;
    COFD[297] = 3.60186887E+00;
    COFD[298] = -2.53302622E-01;
    COFD[299] = 1.09893496E-02;
    COFD[300] = -1.03310318E+01;
    COFD[301] = 1.90522472E+00;
    COFD[302] = -3.44812795E-02;
    COFD[303] = 1.57640018E-03;
    COFD[304] = -1.50371784E+01;
    COFD[305] = 3.26249588E+00;
    COFD[306] = -2.10658287E-01;
    COFD[307] = 9.20032462E-03;
    COFD[308] = -1.81197354E+01;
    COFD[309] = 4.33684042E+00;
    COFD[310] = -3.44981265E-01;
    COFD[311] = 1.48142449E-02;
    COFD[312] = -1.43717529E+01;
    COFD[313] = 3.70920439E+00;
    COFD[314] = -2.67274113E-01;
    COFD[315] = 1.15967481E-02;
    COFD[316] = -1.15806808E+01;
    COFD[317] = 2.43235504E+00;
    COFD[318] = -1.02890179E-01;
    COFD[319] = 4.52903603E-03;
    COFD[320] = -1.40916052E+01;
    COFD[321] = 3.07458927E+00;
    COFD[322] = -1.86899591E-01;
    COFD[323] = 8.19829781E-03;
    COFD[324] = -1.41119732E+01;
    COFD[325] = 3.07458927E+00;
    COFD[326] = -1.86899591E-01;
    COFD[327] = 8.19829781E-03;
    COFD[328] = -2.04177482E+01;
    COFD[329] = 5.31457079E+00;
    COFD[330] = -4.58216496E-01;
    COFD[331] = 1.91825910E-02;
    COFD[332] = -1.53187643E+01;
    COFD[333] = 3.37317428E+00;
    COFD[334] = -2.24900439E-01;
    COFD[335] = 9.81228151E-03;
    COFD[336] = -1.53265780E+01;
    COFD[337] = 3.37317428E+00;
    COFD[338] = -2.24900439E-01;
    COFD[339] = 9.81228151E-03;
    COFD[340] = -1.53340417E+01;
    COFD[341] = 3.37317428E+00;
    COFD[342] = -2.24900439E-01;
    COFD[343] = 9.81228151E-03;
    COFD[344] = -1.50168028E+01;
    COFD[345] = 3.25515933E+00;
    COFD[346] = -2.09710110E-01;
    COFD[347] = 9.15941830E-03;
    COFD[348] = -1.59677692E+01;
    COFD[349] = 3.60186887E+00;
    COFD[350] = -2.53302622E-01;
    COFD[351] = 1.09893496E-02;
    COFD[352] = -1.03327323E+01;
    COFD[353] = 1.90522472E+00;
    COFD[354] = -3.44812795E-02;
    COFD[355] = 1.57640018E-03;
    COFD[356] = -1.50443569E+01;
    COFD[357] = 3.26249588E+00;
    COFD[358] = -2.10658287E-01;
    COFD[359] = 9.20032462E-03;
    COFD[360] = -1.81286555E+01;
    COFD[361] = 4.33684042E+00;
    COFD[362] = -3.44981265E-01;
    COFD[363] = 1.48142449E-02;
    COFD[364] = -1.43721922E+01;
    COFD[365] = 3.70920439E+00;
    COFD[366] = -2.67274113E-01;
    COFD[367] = 1.15967481E-02;
    COFD[368] = -1.15815344E+01;
    COFD[369] = 2.43235504E+00;
    COFD[370] = -1.02890179E-01;
    COFD[371] = 4.52903603E-03;
    COFD[372] = -1.40964661E+01;
    COFD[373] = 3.07458927E+00;
    COFD[374] = -1.86899591E-01;
    COFD[375] = 8.19829781E-03;
    COFD[376] = -1.41170372E+01;
    COFD[377] = 3.07458927E+00;
    COFD[378] = -1.86899591E-01;
    COFD[379] = 8.19829781E-03;
    COFD[380] = -2.04230073E+01;
    COFD[381] = 5.31457079E+00;
    COFD[382] = -4.58216496E-01;
    COFD[383] = 1.91825910E-02;
    COFD[384] = -1.53261114E+01;
    COFD[385] = 3.37317428E+00;
    COFD[386] = -2.24900439E-01;
    COFD[387] = 9.81228151E-03;
    COFD[388] = -1.53340417E+01;
    COFD[389] = 3.37317428E+00;
    COFD[390] = -2.24900439E-01;
    COFD[391] = 9.81228151E-03;
    COFD[392] = -1.53416186E+01;
    COFD[393] = 3.37317428E+00;
    COFD[394] = -2.24900439E-01;
    COFD[395] = 9.81228151E-03;
    COFD[396] = -1.50236516E+01;
    COFD[397] = 3.25515933E+00;
    COFD[398] = -2.09710110E-01;
    COFD[399] = 9.15941830E-03;
    COFD[400] = -1.59759490E+01;
    COFD[401] = 3.60186887E+00;
    COFD[402] = -2.53302622E-01;
    COFD[403] = 1.09893496E-02;
    COFD[404] = -1.03343373E+01;
    COFD[405] = 1.90522472E+00;
    COFD[406] = -3.44812795E-02;
    COFD[407] = 1.57640018E-03;
    COFD[408] = -1.50512053E+01;
    COFD[409] = 3.26249588E+00;
    COFD[410] = -2.10658287E-01;
    COFD[411] = 9.20032462E-03;
    COFD[412] = -1.81371948E+01;
    COFD[413] = 4.33684042E+00;
    COFD[414] = -3.44981265E-01;
    COFD[415] = 1.48142449E-02;
    COFD[416] = -1.40298830E+01;
    COFD[417] = 3.55837688E+00;
    COFD[418] = -2.47785790E-01;
    COFD[419] = 1.07555332E-02;
    COFD[420] = -1.13253458E+01;
    COFD[421] = 2.31195095E+00;
    COFD[422] = -8.63988037E-02;
    COFD[423] = 3.77573452E-03;
    COFD[424] = -1.38756407E+01;
    COFD[425] = 2.98558426E+00;
    COFD[426] = -1.75507216E-01;
    COFD[427] = 7.71173691E-03;
    COFD[428] = -1.38948667E+01;
    COFD[429] = 2.98558426E+00;
    COFD[430] = -1.75507216E-01;
    COFD[431] = 7.71173691E-03;
    COFD[432] = -2.08123325E+01;
    COFD[433] = 5.42470154E+00;
    COFD[434] = -4.69700416E-01;
    COFD[435] = 1.95706904E-02;
    COFD[436] = -1.50096240E+01;
    COFD[437] = 3.25515933E+00;
    COFD[438] = -2.09710110E-01;
    COFD[439] = 9.15941830E-03;
    COFD[440] = -1.50168028E+01;
    COFD[441] = 3.25515933E+00;
    COFD[442] = -2.09710110E-01;
    COFD[443] = 9.15941830E-03;
    COFD[444] = -1.50236516E+01;
    COFD[445] = 3.25515933E+00;
    COFD[446] = -2.09710110E-01;
    COFD[447] = 9.15941830E-03;
    COFD[448] = -1.47639290E+01;
    COFD[449] = 3.15955654E+00;
    COFD[450] = -1.97590757E-01;
    COFD[451] = 8.64692156E-03;
    COFD[452] = -1.57236706E+01;
    COFD[453] = 3.51447210E+00;
    COFD[454] = -2.42579007E-01;
    COFD[455] = 1.05506318E-02;
    COFD[456] = -1.01976409E+01;
    COFD[457] = 1.83188320E+00;
    COFD[458] = -2.40547456E-02;
    COFD[459] = 1.08399898E-03;
    COFD[460] = -1.47850486E+01;
    COFD[461] = 3.16433919E+00;
    COFD[462] = -1.98191564E-01;
    COFD[463] = 8.67209742E-03;
    COFD[464] = -1.77350592E+01;
    COFD[465] = 4.19328271E+00;
    COFD[466] = -3.26911461E-01;
    COFD[467] = 1.40520357E-02;
    COFD[468] = -1.51208119E+01;
    COFD[469] = 3.99904647E+00;
    COFD[470] = -3.03517220E-01;
    COFD[471] = 1.31117363E-02;
    COFD[472] = -1.20638601E+01;
    COFD[473] = 2.63303536E+00;
    COFD[474] = -1.29792632E-01;
    COFD[475] = 5.73363738E-03;
    COFD[476] = -1.47082523E+01;
    COFD[477] = 3.30683499E+00;
    COFD[478] = -2.16378602E-01;
    COFD[479] = 9.44670561E-03;
    COFD[480] = -1.47298720E+01;
    COFD[481] = 3.30683499E+00;
    COFD[482] = -2.16378602E-01;
    COFD[483] = 9.44670561E-03;
    COFD[484] = -2.10785324E+01;
    COFD[485] = 5.51573149E+00;
    COFD[486] = -4.78177665E-01;
    COFD[487] = 1.98082796E-02;
    COFD[488] = -1.59592184E+01;
    COFD[489] = 3.60186887E+00;
    COFD[490] = -2.53302622E-01;
    COFD[491] = 1.09893496E-02;
    COFD[492] = -1.59677692E+01;
    COFD[493] = 3.60186887E+00;
    COFD[494] = -2.53302622E-01;
    COFD[495] = 1.09893496E-02;
    COFD[496] = -1.59759490E+01;
    COFD[497] = 3.60186887E+00;
    COFD[498] = -2.53302622E-01;
    COFD[499] = 1.09893496E-02;
    COFD[500] = -1.57236706E+01;
    COFD[501] = 3.51447210E+00;
    COFD[502] = -2.42579007E-01;
    COFD[503] = 1.05506318E-02;
    COFD[504] = -1.68944722E+01;
    COFD[505] = 3.94346012E+00;
    COFD[506] = -2.96835271E-01;
    COFD[507] = 1.28438696E-02;
    COFD[508] = -1.08140177E+01;
    COFD[509] = 2.11737538E+00;
    COFD[510] = -6.46167749E-02;
    COFD[511] = 2.99827695E-03;
    COFD[512] = -1.57440433E+01;
    COFD[513] = 3.51861272E+00;
    COFD[514] = -2.43068621E-01;
    COFD[515] = 1.05698368E-02;
    COFD[516] = -1.90183510E+01;
    COFD[517] = 4.64763677E+00;
    COFD[518] = -3.82799418E-01;
    COFD[519] = 1.63539171E-02;
    COFD[520] = -9.71338331E+00;
    COFD[521] = 2.17561180E+00;
    COFD[522] = -7.28270090E-02;
    COFD[523] = 3.38302182E-03;
    COFD[524] = -9.86429034E+00;
    COFD[525] = 2.05348746E+00;
    COFD[526] = -5.90289007E-02;
    COFD[527] = 2.89596157E-03;
    COFD[528] = -9.70779324E+00;
    COFD[529] = 1.77912272E+00;
    COFD[530] = -1.67349571E-02;
    COFD[531] = 7.45446845E-04;
    COFD[532] = -9.71375861E+00;
    COFD[533] = 1.77912272E+00;
    COFD[534] = -1.67349571E-02;
    COFD[535] = 7.45446845E-04;
    COFD[536] = -1.21950642E+01;
    COFD[537] = 2.72222246E+00;
    COFD[538] = -1.41335602E-01;
    COFD[539] = 6.23222872E-03;
    COFD[540] = -1.03310318E+01;
    COFD[541] = 1.90522472E+00;
    COFD[542] = -3.44812795E-02;
    COFD[543] = 1.57640018E-03;
    COFD[544] = -1.03327323E+01;
    COFD[545] = 1.90522472E+00;
    COFD[546] = -3.44812795E-02;
    COFD[547] = 1.57640018E-03;
    COFD[548] = -1.03343373E+01;
    COFD[549] = 1.90522472E+00;
    COFD[550] = -3.44812795E-02;
    COFD[551] = 1.57640018E-03;
    COFD[552] = -1.01976409E+01;
    COFD[553] = 1.83188320E+00;
    COFD[554] = -2.40547456E-02;
    COFD[555] = 1.08399898E-03;
    COFD[556] = -1.08140177E+01;
    COFD[557] = 2.11737538E+00;
    COFD[558] = -6.46167749E-02;
    COFD[559] = 2.99827695E-03;
    COFD[560] = -7.72963289E+00;
    COFD[561] = 1.13864728E+00;
    COFD[562] = 7.22991035E-02;
    COFD[563] = -3.32491895E-03;
    COFD[564] = -1.02057322E+01;
    COFD[565] = 1.83104667E+00;
    COFD[566] = -2.39235907E-02;
    COFD[567] = 1.07741763E-03;
    COFD[568] = -1.09328506E+01;
    COFD[569] = 2.05651569E+00;
    COFD[570] = -5.19591463E-02;
    COFD[571] = 2.22384771E-03;
    COFD[572] = -1.40524065E+01;
    COFD[573] = 3.56261348E+00;
    COFD[574] = -2.48287981E-01;
    COFD[575] = 1.07752947E-02;
    COFD[576] = -1.13541075E+01;
    COFD[577] = 2.31999438E+00;
    COFD[578] = -8.75064804E-02;
    COFD[579] = 3.82656365E-03;
    COFD[580] = -1.39007410E+01;
    COFD[581] = 2.99164244E+00;
    COFD[582] = -1.76293106E-01;
    COFD[583] = 7.74575100E-03;
    COFD[584] = -1.39199663E+01;
    COFD[585] = 2.99164244E+00;
    COFD[586] = -1.76293106E-01;
    COFD[587] = 7.74575100E-03;
    COFD[588] = -2.08943798E+01;
    COFD[589] = 5.44718652E+00;
    COFD[590] = -4.72082953E-01;
    COFD[591] = 1.96531321E-02;
    COFD[592] = -1.50371784E+01;
    COFD[593] = 3.26249588E+00;
    COFD[594] = -2.10658287E-01;
    COFD[595] = 9.20032462E-03;
    COFD[596] = -1.50443569E+01;
    COFD[597] = 3.26249588E+00;
    COFD[598] = -2.10658287E-01;
    COFD[599] = 9.20032462E-03;
    COFD[600] = -1.50512053E+01;
    COFD[601] = 3.26249588E+00;
    COFD[602] = -2.10658287E-01;
    COFD[603] = 9.20032462E-03;
    COFD[604] = -1.47850486E+01;
    COFD[605] = 3.16433919E+00;
    COFD[606] = -1.98191564E-01;
    COFD[607] = 8.67209742E-03;
    COFD[608] = -1.57440433E+01;
    COFD[609] = 3.51861272E+00;
    COFD[610] = -2.43068621E-01;
    COFD[611] = 1.05698368E-02;
    COFD[612] = -1.02057322E+01;
    COFD[613] = 1.83104667E+00;
    COFD[614] = -2.39235907E-02;
    COFD[615] = 1.07741763E-03;
    COFD[616] = -1.48061490E+01;
    COFD[617] = 3.16912473E+00;
    COFD[618] = -1.98792456E-01;
    COFD[619] = 8.69726395E-03;
    COFD[620] = -1.77673000E+01;
    COFD[621] = 4.20234040E+00;
    COFD[622] = -3.28057658E-01;
    COFD[623] = 1.41006192E-02;
    COFD[624] = -1.72993972E+01;
    COFD[625] = 4.71931868E+00;
    COFD[626] = -3.91258152E-01;
    COFD[627] = 1.66866639E-02;
    COFD[628] = -1.35545239E+01;
    COFD[629] = 3.13878730E+00;
    COFD[630] = -1.94980335E-01;
    COFD[631] = 8.53744486E-03;
    COFD[632] = -1.67115577E+01;
    COFD[633] = 3.98859394E+00;
    COFD[634] = -3.02316219E-01;
    COFD[635] = 1.30661099E-02;
    COFD[636] = -1.67337768E+01;
    COFD[637] = 3.98859394E+00;
    COFD[638] = -3.02316219E-01;
    COFD[639] = 1.30661099E-02;
    COFD[640] = -2.12021508E+01;
    COFD[641] = 5.20775052E+00;
    COFD[642] = -4.07348327E-01;
    COFD[643] = 1.55473283E-02;
    COFD[644] = -1.81197354E+01;
    COFD[645] = 4.33684042E+00;
    COFD[646] = -3.44981265E-01;
    COFD[647] = 1.48142449E-02;
    COFD[648] = -1.81286555E+01;
    COFD[649] = 4.33684042E+00;
    COFD[650] = -3.44981265E-01;
    COFD[651] = 1.48142449E-02;
    COFD[652] = -1.81371948E+01;
    COFD[653] = 4.33684042E+00;
    COFD[654] = -3.44981265E-01;
    COFD[655] = 1.48142449E-02;
    COFD[656] = -1.77350592E+01;
    COFD[657] = 4.19328271E+00;
    COFD[658] = -3.26911461E-01;
    COFD[659] = 1.40520357E-02;
    COFD[660] = -1.90183510E+01;
    COFD[661] = 4.64763677E+00;
    COFD[662] = -3.82799418E-01;
    COFD[663] = 1.63539171E-02;
    COFD[664] = -1.09328506E+01;
    COFD[665] = 2.05651569E+00;
    COFD[666] = -5.19591463E-02;
    COFD[667] = 2.22384771E-03;
    COFD[668] = -1.77673000E+01;
    COFD[669] = 4.20234040E+00;
    COFD[670] = -3.28057658E-01;
    COFD[671] = 1.41006192E-02;
    COFD[672] = -2.10907727E+01;
    COFD[673] = 5.29211327E+00;
    COFD[674] = -4.56068366E-01;
    COFD[675] = 1.91195062E-02;
}


/*List of specs with small weight, dim NLITE */
void egtransetKTDIF(int* KTDIF) {
    KTDIF[0] = 1;
    KTDIF[1] = 2;
    KTDIF[2] = 11;
}


/*Poly fits for thermal diff ratios, dim NO*NLITE*KK */
void egtransetCOFTD(double* COFTD) {
    COFTD[0] = 0.00000000E+00;
    COFTD[1] = 0.00000000E+00;
    COFTD[2] = 0.00000000E+00;
    COFTD[3] = 0.00000000E+00;
    COFTD[4] = 1.52534742E-01;
    COFTD[5] = 5.46404022E-05;
    COFTD[6] = -2.93412470E-08;
    COFTD[7] = 4.87091914E-12;
    COFTD[8] = 2.70010150E-01;
    COFTD[9] = 3.61555093E-04;
    COFTD[10] = -1.80744752E-07;
    COFTD[11] = 2.75321248E-11;
    COFTD[12] = 2.72041664E-01;
    COFTD[13] = 3.64275376E-04;
    COFTD[14] = -1.82104647E-07;
    COFTD[15] = 2.77392722E-11;
    COFTD[16] = -1.41883744E-01;
    COFTD[17] = 7.66558810E-04;
    COFTD[18] = -3.06550003E-07;
    COFTD[19] = 4.02959502E-11;
    COFTD[20] = 2.20482843E-01;
    COFTD[21] = 4.80164288E-04;
    COFTD[22] = -2.32927944E-07;
    COFTD[23] = 3.46470436E-11;
    COFTD[24] = 2.20907853E-01;
    COFTD[25] = 4.81089870E-04;
    COFTD[26] = -2.33376944E-07;
    COFTD[27] = 3.47138305E-11;
    COFTD[28] = 2.21308399E-01;
    COFTD[29] = 4.81962174E-04;
    COFTD[30] = -2.33800100E-07;
    COFTD[31] = 3.47767730E-11;
    COFTD[32] = 2.40744421E-01;
    COFTD[33] = 4.45343451E-04;
    COFTD[34] = -2.18173874E-07;
    COFTD[35] = 3.26958506E-11;
    COFTD[36] = 1.65429221E-01;
    COFTD[37] = 5.61238922E-04;
    COFTD[38] = -2.65650544E-07;
    COFTD[39] = 3.88229592E-11;
    COFTD[40] = 3.40762433E-01;
    COFTD[41] = -4.04057756E-05;
    COFTD[42] = 3.27879533E-08;
    COFTD[43] = -6.27093812E-12;
    COFTD[44] = 2.39409939E-01;
    COFTD[45] = 4.47197179E-04;
    COFTD[46] = -2.18951702E-07;
    COFTD[47] = 3.27973510E-11;
    COFTD[48] = 2.44369385E-02;
    COFTD[49] = 7.18242498E-04;
    COFTD[50] = -3.19718504E-07;
    COFTD[51] = 4.48828685E-11;
    COFTD[52] = -1.52534742E-01;
    COFTD[53] = -5.46404022E-05;
    COFTD[54] = 2.93412470E-08;
    COFTD[55] = -4.87091914E-12;
    COFTD[56] = 0.00000000E+00;
    COFTD[57] = 0.00000000E+00;
    COFTD[58] = 0.00000000E+00;
    COFTD[59] = 0.00000000E+00;
    COFTD[60] = 4.15583337E-01;
    COFTD[61] = 1.09738399E-05;
    COFTD[62] = -3.96021963E-09;
    COFTD[63] = 1.14414443E-12;
    COFTD[64] = 4.21932443E-01;
    COFTD[65] = 1.11414935E-05;
    COFTD[66] = -4.02072219E-09;
    COFTD[67] = 1.16162418E-12;
    COFTD[68] = 6.02028221E-02;
    COFTD[69] = 5.61561867E-04;
    COFTD[70] = -2.55372862E-07;
    COFTD[71] = 3.63389913E-11;
    COFTD[72] = 4.42739084E-01;
    COFTD[73] = 7.11770818E-05;
    COFTD[74] = -3.84768062E-08;
    COFTD[75] = 6.86323437E-12;
    COFTD[76] = 4.44452569E-01;
    COFTD[77] = 7.14525507E-05;
    COFTD[78] = -3.86257187E-08;
    COFTD[79] = 6.88979640E-12;
    COFTD[80] = 4.46070183E-01;
    COFTD[81] = 7.17126069E-05;
    COFTD[82] = -3.87662996E-08;
    COFTD[83] = 6.91487226E-12;
    COFTD[84] = 4.45261966E-01;
    COFTD[85] = 4.94697174E-05;
    COFTD[86] = -2.63023442E-08;
    COFTD[87] = 4.90306217E-12;
    COFTD[88] = 4.22530228E-01;
    COFTD[89] = 1.32084268E-04;
    COFTD[90] = -7.12222323E-08;
    COFTD[91] = 1.19516090E-11;
    COFTD[92] = 1.61613664E-01;
    COFTD[93] = 4.74155340E-05;
    COFTD[94] = -1.67115247E-08;
    COFTD[95] = -1.88982125E-12;
    COFTD[96] = 4.44653617E-01;
    COFTD[97] = 5.06631704E-05;
    COFTD[98] = -2.69820900E-08;
    COFTD[99] = 5.01289759E-12;
    COFTD[100] = 3.25742450E-01;
    COFTD[101] = 3.03633411E-04;
    COFTD[102] = -1.55290330E-07;
    COFTD[103] = 2.41466436E-11;
    COFTD[104] = -3.40762433E-01;
    COFTD[105] = 4.04057756E-05;
    COFTD[106] = -3.27879533E-08;
    COFTD[107] = 6.27093812E-12;
    COFTD[108] = -1.61613664E-01;
    COFTD[109] = -4.74155340E-05;
    COFTD[110] = 1.67115247E-08;
    COFTD[111] = 1.88982125E-12;
    COFTD[112] = 3.31587939E-01;
    COFTD[113] = -1.96388078E-05;
    COFTD[114] = 3.02388828E-08;
    COFTD[115] = -8.44998018E-12;
    COFTD[116] = 3.42203127E-01;
    COFTD[117] = -2.02675087E-05;
    COFTD[118] = 3.12069259E-08;
    COFTD[119] = -8.72049099E-12;
    COFTD[120] = 2.84983505E-01;
    COFTD[121] = 1.15460005E-04;
    COFTD[122] = -6.17197869E-08;
    COFTD[123] = 1.01504212E-11;
    COFTD[124] = 4.40220831E-01;
    COFTD[125] = -4.83717413E-05;
    COFTD[126] = 4.66088897E-08;
    COFTD[127] = -1.02768430E-11;
    COFTD[128] = 4.43649137E-01;
    COFTD[129] = -4.87484458E-05;
    COFTD[130] = 4.69718656E-08;
    COFTD[131] = -1.03568760E-11;
    COFTD[132] = 4.46895651E-01;
    COFTD[133] = -4.91051748E-05;
    COFTD[134] = 4.73155940E-08;
    COFTD[135] = -1.04326650E-11;
    COFTD[136] = 4.22009934E-01;
    COFTD[137] = -4.14042334E-05;
    COFTD[138] = 4.38751613E-08;
    COFTD[139] = -1.02860246E-11;
    COFTD[140] = 4.66315159E-01;
    COFTD[141] = -5.60150425E-05;
    COFTD[142] = 4.65987669E-08;
    COFTD[143] = -9.13646318E-12;
    COFTD[144] = 0.00000000E+00;
    COFTD[145] = 0.00000000E+00;
    COFTD[146] = 0.00000000E+00;
    COFTD[147] = 0.00000000E+00;
    COFTD[148] = 4.22171414E-01;
    COFTD[149] = -4.17749918E-05;
    COFTD[150] = 4.39726219E-08;
    COFTD[151] = -1.02672932E-11;
    COFTD[152] = 4.59663274E-01;
    COFTD[153] = -1.74770868E-05;
    COFTD[154] = 1.42888118E-08;
    COFTD[155] = -2.03610705E-12;
}

/* Replace this routine with the one generated by the Gauss Jordan solver of DW */
AMREX_GPU_HOST_DEVICE void sgjsolve(double* A, double* x, double* b) {
    amrex::Abort("sgjsolve not implemented, choose a different solver ");
}

/* Replace this routine with the one generated by the Gauss Jordan solver of DW */
AMREX_GPU_HOST_DEVICE void sgjsolve_simplified(double* A, double* x, double* b) {
    amrex::Abort("sgjsolve_simplified not implemented, choose a different solver ");
}

