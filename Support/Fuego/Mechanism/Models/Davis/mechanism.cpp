#include "chemistry_file.H"

#ifndef AMREX_USE_CUDA
namespace thermo
{
    double fwd_A[38], fwd_beta[38], fwd_Ea[38];
    double low_A[38], low_beta[38], low_Ea[38];
    double rev_A[38], rev_beta[38], rev_Ea[38];
    double troe_a[38],troe_Ts[38], troe_Tss[38], troe_Tsss[38];
    double sri_a[38], sri_b[38], sri_c[38], sri_d[38], sri_e[38];
    double activation_units[38], prefactor_units[38], phase_units[38];
    int is_PD[38], troe_len[38], sri_len[38], nTB[38], *TBid[38];
    double *TB[38];
    std::vector<std::vector<double>> kiv(38); 
    std::vector<std::vector<double>> nuv(38); 

    double fwd_A_DEF[38], fwd_beta_DEF[38], fwd_Ea_DEF[38];
    double low_A_DEF[38], low_beta_DEF[38], low_Ea_DEF[38];
    double rev_A_DEF[38], rev_beta_DEF[38], rev_Ea_DEF[38];
    double troe_a_DEF[38],troe_Ts_DEF[38], troe_Tss_DEF[38], troe_Tsss_DEF[38];
    double sri_a_DEF[38], sri_b_DEF[38], sri_c_DEF[38], sri_d_DEF[38], sri_e_DEF[38];
    double activation_units_DEF[38], prefactor_units_DEF[38], phase_units_DEF[38];
    int is_PD_DEF[38], troe_len_DEF[38], sri_len_DEF[38], nTB_DEF[38], *TBid_DEF[38];
    double *TB_DEF[38];
    std::vector<int> rxn_map;
};

using namespace thermo;
#endif

/* Inverse molecular weights */
/* TODO: check necessity on CPU */
static AMREX_GPU_DEVICE_MANAGED double imw[14] = {
    1.0 / 2.015940,  /*H2 */
    1.0 / 1.007970,  /*H */
    1.0 / 39.948000,  /*AR */
    1.0 / 28.013400,  /*N2 */
    1.0 / 4.002600,  /*HE */
    1.0 / 15.999400,  /*O */
    1.0 / 17.007370,  /*OH */
    1.0 / 29.018520,  /*HCO */
    1.0 / 33.006770,  /*HO2 */
    1.0 / 18.015340,  /*H2O */
    1.0 / 28.010550,  /*CO */
    1.0 / 31.998800,  /*O2 */
    1.0 / 34.014740,  /*H2O2 */
    1.0 / 44.009950};  /*CO2 */

/* Inverse molecular weights */
/* TODO: check necessity because redundant with molecularWeight */
static AMREX_GPU_DEVICE_MANAGED double molecular_weights[14] = {
    2.015940,  /*H2 */
    1.007970,  /*H */
    39.948000,  /*AR */
    28.013400,  /*N2 */
    4.002600,  /*HE */
    15.999400,  /*O */
    17.007370,  /*OH */
    29.018520,  /*HCO */
    33.006770,  /*HO2 */
    18.015340,  /*H2O */
    28.010550,  /*CO */
    31.998800,  /*O2 */
    34.014740,  /*H2O2 */
    44.009950};  /*CO2 */

AMREX_GPU_HOST_DEVICE
void get_imw(double imw_new[]){
    for(int i = 0; i<14; ++i) imw_new[i] = imw[i];
}

/* TODO: check necessity because redundant with CKWT */
AMREX_GPU_HOST_DEVICE
void get_mw(double mw_new[]){
    for(int i = 0; i<14; ++i) mw_new[i] = molecular_weights[i];
}


#ifndef AMREX_USE_CUDA
/* Initializes parameter database */
void CKINIT()
{

    rxn_map = {8,9,10,11,3,12,13,14,4,5,6,0,15,1,16,17,18,19,20,21,22,23,24,25,26,27,2,28,29,30,31,32,33,34,35,7,36,37};

    // (0):  H + O2 <=> O + OH
    kiv[8] = {1,11,5,6};
    nuv[8] = {-1,-1,1,1};
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
    kiv[9] = {5,0,1,6};
    nuv[9] = {-1,-1,1,1};
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
    kiv[10] = {6,0,1,9};
    nuv[10] = {-1,-1,1,1};
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
    kiv[11] = {6,6,5,9};
    nuv[11] = {-1,-1,1,1};
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
    kiv[3] = {1,1,0};
    nuv[3] = {-1,-1,1};
    // (4):  H + H + M <=> H2 + M
    fwd_A[3]     = 1.78e+18;
    fwd_beta[3]  = -1;
    fwd_Ea[3]    = 0;
    prefactor_units[3]  = 1.0000000000000002e-12;
    activation_units[3] = 0.50321666580471969;
    phase_units[3]      = pow(10,-12.000000);
    is_PD[3] = 0;
    nTB[3] = 5;
    TB[3] = (double *) malloc(5 * sizeof(double));
    TBid[3] = (int *) malloc(5 * sizeof(int));
    TBid[3][0] = 0; TB[3][0] = 0; // H2
    TBid[3][1] = 9; TB[3][1] = 0; // H2O
    TBid[3][2] = 13; TB[3][2] = 0; // CO2
    TBid[3][3] = 2; TB[3][3] = 0.63; // AR
    TBid[3][4] = 4; TB[3][4] = 0.63; // HE

    // (5):  H + H + H2 <=> H2 + H2
    kiv[12] = {1,1,0,0,0};
    nuv[12] = {-1,-1,-1,1,1};
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
    kiv[13] = {1,1,9,0,9};
    nuv[13] = {-1,-1,-1,1,1};
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
    kiv[14] = {1,1,13,0,13};
    nuv[14] = {-1,-1,-1,1,1};
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
    kiv[4] = {1,6,9};
    nuv[4] = {-1,-1,1};
    // (8):  H + OH + M <=> H2O + M
    fwd_A[4]     = 4.4e+22;
    fwd_beta[4]  = -2;
    fwd_Ea[4]    = 0;
    prefactor_units[4]  = 1.0000000000000002e-12;
    activation_units[4] = 0.50321666580471969;
    phase_units[4]      = pow(10,-12.000000);
    is_PD[4] = 0;
    nTB[4] = 6;
    TB[4] = (double *) malloc(6 * sizeof(double));
    TBid[4] = (int *) malloc(6 * sizeof(int));
    TBid[4][0] = 0; TB[4][0] = 2; // H2
    TBid[4][1] = 9; TB[4][1] = 6.2999999999999998; // H2O
    TBid[4][2] = 10; TB[4][2] = 1.75; // CO
    TBid[4][3] = 13; TB[4][3] = 3.6000000000000001; // CO2
    TBid[4][4] = 2; TB[4][4] = 0.38; // AR
    TBid[4][5] = 4; TB[4][5] = 0.38; // HE

    // (9):  O + H + M <=> OH + M
    kiv[5] = {5,1,6};
    nuv[5] = {-1,-1,1};
    // (9):  O + H + M <=> OH + M
    fwd_A[5]     = 9.428e+18;
    fwd_beta[5]  = -1;
    fwd_Ea[5]    = 0;
    prefactor_units[5]  = 1.0000000000000002e-12;
    activation_units[5] = 0.50321666580471969;
    phase_units[5]      = pow(10,-12.000000);
    is_PD[5] = 0;
    nTB[5] = 6;
    TB[5] = (double *) malloc(6 * sizeof(double));
    TBid[5] = (int *) malloc(6 * sizeof(int));
    TBid[5][0] = 0; TB[5][0] = 2; // H2
    TBid[5][1] = 9; TB[5][1] = 12; // H2O
    TBid[5][2] = 10; TB[5][2] = 1.75; // CO
    TBid[5][3] = 13; TB[5][3] = 3.6000000000000001; // CO2
    TBid[5][4] = 2; TB[5][4] = 0.69999999999999996; // AR
    TBid[5][5] = 4; TB[5][5] = 0.69999999999999996; // HE

    // (10):  O + O + M <=> O2 + M
    kiv[6] = {5,5,11};
    nuv[6] = {-1,-1,1};
    // (10):  O + O + M <=> O2 + M
    fwd_A[6]     = 1.2e+17;
    fwd_beta[6]  = -1;
    fwd_Ea[6]    = 0;
    prefactor_units[6]  = 1.0000000000000002e-12;
    activation_units[6] = 0.50321666580471969;
    phase_units[6]      = pow(10,-12.000000);
    is_PD[6] = 0;
    nTB[6] = 6;
    TB[6] = (double *) malloc(6 * sizeof(double));
    TBid[6] = (int *) malloc(6 * sizeof(int));
    TBid[6][0] = 0; TB[6][0] = 2.3999999999999999; // H2
    TBid[6][1] = 9; TB[6][1] = 15.4; // H2O
    TBid[6][2] = 10; TB[6][2] = 1.75; // CO
    TBid[6][3] = 13; TB[6][3] = 3.6000000000000001; // CO2
    TBid[6][4] = 2; TB[6][4] = 0.82999999999999996; // AR
    TBid[6][5] = 4; TB[6][5] = 0.82999999999999996; // HE

    // (11):  H + O2 (+M) <=> HO2 (+M)
    kiv[0] = {1,11,8};
    nuv[0] = {-1,-1,1};
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
    TB[0] = (double *) malloc(7 * sizeof(double));
    TBid[0] = (int *) malloc(7 * sizeof(int));
    TBid[0][0] = 11; TB[0][0] = 0.84999999999999998; // O2
    TBid[0][1] = 9; TB[0][1] = 11.890000000000001; // H2O
    TBid[0][2] = 10; TB[0][2] = 1.0900000000000001; // CO
    TBid[0][3] = 13; TB[0][3] = 2.1800000000000002; // CO2
    TBid[0][4] = 2; TB[0][4] = 0.40000000000000002; // AR
    TBid[0][5] = 4; TB[0][5] = 0.46000000000000002; // HE
    TBid[0][6] = 0; TB[0][6] = 0.75; // H2

    // (12):  H2 + O2 <=> HO2 + H
    kiv[15] = {0,11,8,1};
    nuv[15] = {-1,-1,1,1};
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
    kiv[1] = {6,6,12};
    nuv[1] = {-1,-1,1};
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
    TB[1] = (double *) malloc(6 * sizeof(double));
    TBid[1] = (int *) malloc(6 * sizeof(int));
    TBid[1][0] = 0; TB[1][0] = 2; // H2
    TBid[1][1] = 9; TB[1][1] = 6; // H2O
    TBid[1][2] = 10; TB[1][2] = 1.75; // CO
    TBid[1][3] = 13; TB[1][3] = 3.6000000000000001; // CO2
    TBid[1][4] = 2; TB[1][4] = 0.69999999999999996; // AR
    TBid[1][5] = 4; TB[1][5] = 0.69999999999999996; // HE

    // (14):  HO2 + H <=> O + H2O
    kiv[16] = {8,1,5,9};
    nuv[16] = {-1,-1,1,1};
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
    kiv[17] = {8,1,6,6};
    nuv[17] = {-1,-1,1,1};
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
    kiv[18] = {8,5,6,11};
    nuv[18] = {-1,-1,1,1};
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
    kiv[19] = {8,6,11,9};
    nuv[19] = {-1,-1,1,1};
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
    kiv[20] = {8,6,11,9};
    nuv[20] = {-1,-1,1,1};
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
    kiv[21] = {8,8,11,12};
    nuv[21] = {-1,-1,1,1};
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
    kiv[22] = {8,8,11,12};
    nuv[22] = {-1,-1,1,1};
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
    kiv[23] = {12,1,8,0};
    nuv[23] = {-1,-1,1,1};
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
    kiv[24] = {12,1,6,9};
    nuv[24] = {-1,-1,1,1};
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
    kiv[25] = {12,5,6,8};
    nuv[25] = {-1,-1,1,1};
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
    kiv[26] = {12,6,8,9};
    nuv[26] = {-1,-1,1,1};
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
    kiv[27] = {12,6,8,9};
    nuv[27] = {-1,-1,1,1};
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
    kiv[2] = {10,5,13};
    nuv[2] = {-1,-1,1};
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
    TB[2] = (double *) malloc(6 * sizeof(double));
    TBid[2] = (int *) malloc(6 * sizeof(int));
    TBid[2][0] = 0; TB[2][0] = 2; // H2
    TBid[2][1] = 9; TB[2][1] = 12; // H2O
    TBid[2][2] = 10; TB[2][2] = 1.75; // CO
    TBid[2][3] = 13; TB[2][3] = 3.6000000000000001; // CO2
    TBid[2][4] = 2; TB[2][4] = 0.69999999999999996; // AR
    TBid[2][5] = 4; TB[2][5] = 0.69999999999999996; // HE

    // (27):  CO + OH <=> CO2 + H
    kiv[28] = {10,6,13,1};
    nuv[28] = {-1,-1,1,1};
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
    kiv[29] = {10,6,13,1};
    nuv[29] = {-1,-1,1,1};
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
    kiv[30] = {10,11,13,5};
    nuv[30] = {-1,-1,1,1};
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
    kiv[31] = {10,8,13,6};
    nuv[31] = {-1,-1,1,1};
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
    kiv[32] = {7,1,10,0};
    nuv[32] = {-1,-1,1,1};
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
    kiv[33] = {7,5,10,6};
    nuv[33] = {-1,-1,1,1};
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
    kiv[34] = {7,5,13,1};
    nuv[34] = {-1,-1,1,1};
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
    kiv[35] = {7,6,10,9};
    nuv[35] = {-1,-1,1,1};
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
    kiv[7] = {7,10,1};
    nuv[7] = {-1,1,1};
    // (35):  HCO + M <=> CO + H + M
    fwd_A[7]     = 1.87e+17;
    fwd_beta[7]  = -1;
    fwd_Ea[7]    = 17000;
    prefactor_units[7]  = 1.0000000000000002e-06;
    activation_units[7] = 0.50321666580471969;
    phase_units[7]      = pow(10,-6.000000);
    is_PD[7] = 0;
    nTB[7] = 4;
    TB[7] = (double *) malloc(4 * sizeof(double));
    TBid[7] = (int *) malloc(4 * sizeof(int));
    TBid[7][0] = 0; TB[7][0] = 2; // H2
    TBid[7][1] = 9; TB[7][1] = 0; // H2O
    TBid[7][2] = 10; TB[7][2] = 1.75; // CO
    TBid[7][3] = 13; TB[7][3] = 3.6000000000000001; // CO2

    // (36):  HCO + H2O <=> CO + H + H2O
    kiv[36] = {7,9,10,1,9};
    nuv[36] = {-1,-1,1,1,1};
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
    kiv[37] = {7,11,10,8};
    nuv[37] = {-1,-1,1,1};
    // (37):  HCO + O2 <=> CO + HO2
    fwd_A[37]     = 12040000000;
    fwd_beta[37]  = 0.80700000000000005;
    fwd_Ea[37]    = -727;
    prefactor_units[37]  = 1.0000000000000002e-06;
    activation_units[37] = 0.50321666580471969;
    phase_units[37]      = pow(10,-12.000000);
    is_PD[37] = 0;
    nTB[37] = 0;

    SetAllDefaults();
}

void GET_REACTION_MAP(int *rmap)
{
    for (int i=0; i<38; ++i) {
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
  if (reaction_id<0 || reaction_id>=38) {
    printf("Bad reaction id = %d",reaction_id);
    abort();
  };
  int mrid = rxn_map[reaction_id];

  if (param_id == THIRD_BODY) {
    if (species_id<0 || species_id>=14) {
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
    for (int i=0; i<38; i++) {
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
    for (int i=0; i<38; i++) {
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
  for (int i=0; i<38; ++i) {
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
    *kk = 14;
    *ii = 38;
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
    ename[0] = "O";
    ename[1] = "H";
    ename[2] = "C";
    ename[3] = "N";
    ename[4] = "AR";
    ename[5] = "HE";
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

    /* O  */
    kname[ 0*lenkname + 0 ] = 'O';
    kname[ 0*lenkname + 1 ] = ' ';

    /* H  */
    kname[ 1*lenkname + 0 ] = 'H';
    kname[ 1*lenkname + 1 ] = ' ';

    /* C  */
    kname[ 2*lenkname + 0 ] = 'C';
    kname[ 2*lenkname + 1 ] = ' ';

    /* N  */
    kname[ 3*lenkname + 0 ] = 'N';
    kname[ 3*lenkname + 1 ] = ' ';

    /* AR  */
    kname[ 4*lenkname + 0 ] = 'A';
    kname[ 4*lenkname + 1 ] = 'R';
    kname[ 4*lenkname + 2 ] = ' ';

    /* HE  */
    kname[ 5*lenkname + 0 ] = 'H';
    kname[ 5*lenkname + 1 ] = 'E';
    kname[ 5*lenkname + 2 ] = ' ';

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


/* Returns the char strings of species names */
void CKSYMS(int * kname, int * plenkname )
{
    int i; /*Loop Counter */
    int lenkname = *plenkname;
    /*clear kname */
    for (i=0; i<lenkname*14; i++) {
        kname[i] = ' ';
    }

    /* H2  */
    kname[ 0*lenkname + 0 ] = 'H';
    kname[ 0*lenkname + 1 ] = '2';
    kname[ 0*lenkname + 2 ] = ' ';

    /* H  */
    kname[ 1*lenkname + 0 ] = 'H';
    kname[ 1*lenkname + 1 ] = ' ';

    /* AR  */
    kname[ 2*lenkname + 0 ] = 'A';
    kname[ 2*lenkname + 1 ] = 'R';
    kname[ 2*lenkname + 2 ] = ' ';

    /* N2  */
    kname[ 3*lenkname + 0 ] = 'N';
    kname[ 3*lenkname + 1 ] = '2';
    kname[ 3*lenkname + 2 ] = ' ';

    /* HE  */
    kname[ 4*lenkname + 0 ] = 'H';
    kname[ 4*lenkname + 1 ] = 'E';
    kname[ 4*lenkname + 2 ] = ' ';

    /* O  */
    kname[ 5*lenkname + 0 ] = 'O';
    kname[ 5*lenkname + 1 ] = ' ';

    /* OH  */
    kname[ 6*lenkname + 0 ] = 'O';
    kname[ 6*lenkname + 1 ] = 'H';
    kname[ 6*lenkname + 2 ] = ' ';

    /* HCO  */
    kname[ 7*lenkname + 0 ] = 'H';
    kname[ 7*lenkname + 1 ] = 'C';
    kname[ 7*lenkname + 2 ] = 'O';
    kname[ 7*lenkname + 3 ] = ' ';

    /* HO2  */
    kname[ 8*lenkname + 0 ] = 'H';
    kname[ 8*lenkname + 1 ] = 'O';
    kname[ 8*lenkname + 2 ] = '2';
    kname[ 8*lenkname + 3 ] = ' ';

    /* H2O  */
    kname[ 9*lenkname + 0 ] = 'H';
    kname[ 9*lenkname + 1 ] = '2';
    kname[ 9*lenkname + 2 ] = 'O';
    kname[ 9*lenkname + 3 ] = ' ';

    /* CO  */
    kname[ 10*lenkname + 0 ] = 'C';
    kname[ 10*lenkname + 1 ] = 'O';
    kname[ 10*lenkname + 2 ] = ' ';

    /* O2  */
    kname[ 11*lenkname + 0 ] = 'O';
    kname[ 11*lenkname + 1 ] = '2';
    kname[ 11*lenkname + 2 ] = ' ';

    /* H2O2  */
    kname[ 12*lenkname + 0 ] = 'H';
    kname[ 12*lenkname + 1 ] = '2';
    kname[ 12*lenkname + 2 ] = 'O';
    kname[ 12*lenkname + 3 ] = '2';
    kname[ 12*lenkname + 4 ] = ' ';

    /* CO2  */
    kname[ 13*lenkname + 0 ] = 'C';
    kname[ 13*lenkname + 1 ] = 'O';
    kname[ 13*lenkname + 2 ] = '2';
    kname[ 13*lenkname + 3 ] = ' ';

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
    *P = *rho * 8.31446261815324e+07 * (*T) / XW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(y) */
AMREX_GPU_HOST_DEVICE void CKPY(double *  rho, double *  T, double *  y,  double *  P)
{
    double YOW = 0;/* for computing mean MW */
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
    *P = *rho * 8.31446261815324e+07 * (*T) * YOW; /*P = rho*R*T/W */

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
#endif


/*Compute P = rhoRT/W(c) */
void CKPC(double *  rho, double *  T, double *  c,  double *  P)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    double W = 0;
    double sumC = 0;
    W += c[0]*2.015940; /*H2 */
    W += c[1]*1.007970; /*H */
    W += c[2]*39.948000; /*AR */
    W += c[3]*28.013400; /*N2 */
    W += c[4]*4.002600; /*HE */
    W += c[5]*15.999400; /*O */
    W += c[6]*17.007370; /*OH */
    W += c[7]*29.018520; /*HCO */
    W += c[8]*33.006770; /*HO2 */
    W += c[9]*18.015340; /*H2O */
    W += c[10]*28.010550; /*CO */
    W += c[11]*31.998800; /*O2 */
    W += c[12]*34.014740; /*H2O2 */
    W += c[13]*44.009950; /*CO2 */

    for (id = 0; id < 14; ++id) {
        sumC += c[id];
    }
    *P = *rho * 8.31446261815324e+07 * (*T) * sumC / W; /*P = rho*R*T/W */

    return;
}


/*Compute rho = PW(x)/RT */
void CKRHOX(double *  P, double *  T, double *  x,  double *  rho)
{
    double XW = 0;/* To hold mean molecular wt */
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
    *rho = *P * XW / (8.31446261815324e+07 * (*T)); /*rho = P*W/(R*T) */

    return;
}


/*Compute rho = P*W(y)/RT */
AMREX_GPU_HOST_DEVICE void CKRHOY(double *  P, double *  T, double *  y,  double *  rho)
{
    double YOW = 0;
    double tmp[14];

    for (int i = 0; i < 14; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 14; i++)
    {
        YOW += tmp[i];
    }

    *rho = *P / (8.31446261815324e+07 * (*T) * YOW);/*rho = P*W/(R*T) */
    return;
}


/*Compute rho = P*W(c)/(R*T) */
void CKRHOC(double *  P, double *  T, double *  c,  double *  rho)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    double W = 0;
    double sumC = 0;
    W += c[0]*2.015940; /*H2 */
    W += c[1]*1.007970; /*H */
    W += c[2]*39.948000; /*AR */
    W += c[3]*28.013400; /*N2 */
    W += c[4]*4.002600; /*HE */
    W += c[5]*15.999400; /*O */
    W += c[6]*17.007370; /*OH */
    W += c[7]*29.018520; /*HCO */
    W += c[8]*33.006770; /*HO2 */
    W += c[9]*18.015340; /*H2O */
    W += c[10]*28.010550; /*CO */
    W += c[11]*31.998800; /*O2 */
    W += c[12]*34.014740; /*H2O2 */
    W += c[13]*44.009950; /*CO2 */

    for (id = 0; id < 14; ++id) {
        sumC += c[id];
    }
    *rho = *P * W / (sumC * (*T) * 8.31446261815324e+07); /*rho = PW/(R*T) */

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
    double tmp[14];

    for (int i = 0; i < 14; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 14; i++)
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
    W += c[0]*2.015940; /*H2 */
    W += c[1]*1.007970; /*H */
    W += c[2]*39.948000; /*AR */
    W += c[3]*28.013400; /*N2 */
    W += c[4]*4.002600; /*HE */
    W += c[5]*15.999400; /*O */
    W += c[6]*17.007370; /*OH */
    W += c[7]*29.018520; /*HCO */
    W += c[8]*33.006770; /*HO2 */
    W += c[9]*18.015340; /*H2O */
    W += c[10]*28.010550; /*CO */
    W += c[11]*31.998800; /*O2 */
    W += c[12]*34.014740; /*H2O2 */
    W += c[13]*44.009950; /*CO2 */

    for (id = 0; id < 14; ++id) {
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
    double tmp[14];

    for (int i = 0; i < 14; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 14; i++)
    {
        YOW += tmp[i];
    }

    double YOWINV = 1.0/YOW;

    for (int i = 0; i < 14; i++)
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
    for (int i = 0; i < 14; i++)
    {
        c[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 14; i++)
    {
        YOW += c[i];
    }

    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31446261815324e+07 * (*T)); 
    /*Now compute conversion */

    for (int i = 0; i < 14; i++)
    {
        c[i] = PWORT * y[i] * imw[i];
    }
    return;
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
AMREX_GPU_HOST_DEVICE void CKYTCR(double *  rho, double *  T, double *  y,  double *  c)
{
    for (int i = 0; i < 14; i++)
    {
        c[i] = (*rho)  * y[i] * imw[i];
    }
}


/*convert x[species] (mole fracs) to y[species] (mass fracs) */
AMREX_GPU_HOST_DEVICE void CKXTY(double *  x,  double *  y)
{
    double XW = 0; /*See Eq 4, 9 in CK Manual */
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
    /*Now compute conversion */
    double XWinv = 1.0/XW;
    y[0] = x[0]*2.015940*XWinv; 
    y[1] = x[1]*1.007970*XWinv; 
    y[2] = x[2]*39.948000*XWinv; 
    y[3] = x[3]*28.013400*XWinv; 
    y[4] = x[4]*4.002600*XWinv; 
    y[5] = x[5]*15.999400*XWinv; 
    y[6] = x[6]*17.007370*XWinv; 
    y[7] = x[7]*29.018520*XWinv; 
    y[8] = x[8]*33.006770*XWinv; 
    y[9] = x[9]*18.015340*XWinv; 
    y[10] = x[10]*28.010550*XWinv; 
    y[11] = x[11]*31.998800*XWinv; 
    y[12] = x[12]*34.014740*XWinv; 
    y[13] = x[13]*44.009950*XWinv; 

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
void CKXTCP(double *  P, double *  T, double *  x,  double *  c)
{
    int id; /*loop counter */
    double PORT = (*P)/(8.31446261815324e+07 * (*T)); /*P/RT */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 14; ++id) {
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
    ROW = (*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 14; ++id) {
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
    for (id = 0; id < 14; ++id) {
        sumC += c[id];
    }

    /* See Eq 13  */
    double sumCinv = 1.0/sumC;
    for (id = 0; id < 14; ++id) {
        x[id] = c[id]*sumCinv;
    }

    return;
}


/*convert c[species] (molar conc) to y[species] (mass fracs) */
void CKCTY(double *  c, double *  y)
{
    double CW = 0; /*See Eq 12 in CK Manual */
    /*compute denominator in eq 12 first */
    CW += c[0]*2.015940; /*H2 */
    CW += c[1]*1.007970; /*H */
    CW += c[2]*39.948000; /*AR */
    CW += c[3]*28.013400; /*N2 */
    CW += c[4]*4.002600; /*HE */
    CW += c[5]*15.999400; /*O */
    CW += c[6]*17.007370; /*OH */
    CW += c[7]*29.018520; /*HCO */
    CW += c[8]*33.006770; /*HO2 */
    CW += c[9]*18.015340; /*H2O */
    CW += c[10]*28.010550; /*CO */
    CW += c[11]*31.998800; /*O2 */
    CW += c[12]*34.014740; /*H2O2 */
    CW += c[13]*44.009950; /*CO2 */
    /*Now compute conversion */
    double CWinv = 1.0/CW;
    y[0] = c[0]*2.015940*CWinv; 
    y[1] = c[1]*1.007970*CWinv; 
    y[2] = c[2]*39.948000*CWinv; 
    y[3] = c[3]*28.013400*CWinv; 
    y[4] = c[4]*4.002600*CWinv; 
    y[5] = c[5]*15.999400*CWinv; 
    y[6] = c[6]*17.007370*CWinv; 
    y[7] = c[7]*29.018520*CWinv; 
    y[8] = c[8]*33.006770*CWinv; 
    y[9] = c[9]*18.015340*CWinv; 
    y[10] = c[10]*28.010550*CWinv; 
    y[11] = c[11]*31.998800*CWinv; 
    y[12] = c[12]*34.014740*CWinv; 
    y[13] = c[13]*44.009950*CWinv; 

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
    for (id = 0; id < 14; ++id) {
        cvml[id] *= 8.31446261815324e+07;
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
    for (id = 0; id < 14; ++id) {
        cpml[id] *= 8.31446261815324e+07;
    }
}


/*get internal energy as a function  */
/*of T for all species (molar units) */
void CKUML(double *  T,  double *  uml)
{
    int id; /*loop counter */
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31446261815324e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 14; ++id) {
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
    double RT = 8.31446261815324e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 14; ++id) {
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
    double RT = 8.31446261815324e+07*tT; /*R*T */
    gibbs(gml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 14; ++id) {
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
    double RT = 8.31446261815324e+07*tT; /*R*T */
    helmholtz(aml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 14; ++id) {
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
    for (id = 0; id < 14; ++id) {
        sml[id] *= 8.31446261815324e+07;
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
    cvms[0] *= 4.124360158612479e+07; /*H2 */
    cvms[1] *= 8.248720317224957e+07; /*H */
    cvms[2] *= 2.081321372322329e+06; /*AR */
    cvms[3] *= 2.968030520448514e+06; /*N2 */
    cvms[4] *= 2.077265432007505e+07; /*HE */
    cvms[5] *= 5.196734013871295e+06; /*O */
    cvms[6] *= 4.888740950630956e+06; /*OH */
    cvms[7] *= 2.865226282440744e+06; /*HCO */
    cvms[8] *= 2.519017346487778e+06; /*HO2 */
    cvms[9] *= 4.615212712140454e+06; /*H2O */
    cvms[10] *= 2.968332509769797e+06; /*CO */
    cvms[11] *= 2.598367006935648e+06; /*O2 */
    cvms[12] *= 2.444370475315478e+06; /*H2O2 */
    cvms[13] *= 1.889223372931176e+06; /*CO2 */
}


/*Returns the specific heats at constant pressure */
/*in mass units (Eq. 26) */
AMREX_GPU_HOST_DEVICE void CKCPMS(double *  T,  double *  cpms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpms, tc);
    /*multiply by R/molecularweight */
    cpms[0] *= 4.124360158612479e+07; /*H2 */
    cpms[1] *= 8.248720317224957e+07; /*H */
    cpms[2] *= 2.081321372322329e+06; /*AR */
    cpms[3] *= 2.968030520448514e+06; /*N2 */
    cpms[4] *= 2.077265432007505e+07; /*HE */
    cpms[5] *= 5.196734013871295e+06; /*O */
    cpms[6] *= 4.888740950630956e+06; /*OH */
    cpms[7] *= 2.865226282440744e+06; /*HCO */
    cpms[8] *= 2.519017346487778e+06; /*HO2 */
    cpms[9] *= 4.615212712140454e+06; /*H2O */
    cpms[10] *= 2.968332509769797e+06; /*CO */
    cpms[11] *= 2.598367006935648e+06; /*O2 */
    cpms[12] *= 2.444370475315478e+06; /*H2O2 */
    cpms[13] *= 1.889223372931176e+06; /*CO2 */
}


/*Returns internal energy in mass units (Eq 30.) */
AMREX_GPU_HOST_DEVICE void CKUMS(double *  T,  double *  ums)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31446261815324e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    for (int i = 0; i < 14; i++)
    {
        ums[i] *= RT*imw[i];
    }
}


/*Returns enthalpy in mass units (Eq 27.) */
AMREX_GPU_HOST_DEVICE void CKHMS(double *  T,  double *  hms)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31446261815324e+07*tT; /*R*T */
    speciesEnthalpy(hms, tc);
    for (int i = 0; i < 14; i++)
    {
        hms[i] *= RT*imw[i];
    }
}


#ifndef AMREX_USE_CUDA
/*Returns enthalpy in mass units (Eq 27.) */
void VCKHMS(int *  np, double *  T,  double *  hms)
{
    double tc[5], h[14];

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
    double RT = 8.31446261815324e+07*tT; /*R*T */
    gibbs(gms, tc);
    for (int i = 0; i < 14; i++)
    {
        gms[i] *= RT*imw[i];
    }
}


/*Returns helmholtz in mass units (Eq 32.) */
void CKAMS(double *  T,  double *  ams)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31446261815324e+07*tT; /*R*T */
    helmholtz(ams, tc);
    for (int i = 0; i < 14; i++)
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
    sms[0] *= 4.124360158612479e+07; /*H2 */
    sms[1] *= 8.248720317224957e+07; /*H */
    sms[2] *= 2.081321372322329e+06; /*AR */
    sms[3] *= 2.968030520448514e+06; /*N2 */
    sms[4] *= 2.077265432007505e+07; /*HE */
    sms[5] *= 5.196734013871295e+06; /*O */
    sms[6] *= 4.888740950630956e+06; /*OH */
    sms[7] *= 2.865226282440744e+06; /*HCO */
    sms[8] *= 2.519017346487778e+06; /*HO2 */
    sms[9] *= 4.615212712140454e+06; /*H2O */
    sms[10] *= 2.968332509769797e+06; /*CO */
    sms[11] *= 2.598367006935648e+06; /*O2 */
    sms[12] *= 2.444370475315478e+06; /*H2O2 */
    sms[13] *= 1.889223372931176e+06; /*CO2 */
}


/*Returns the mean specific heat at CP (Eq. 33) */
void CKCPBL(double *  T, double *  x,  double *  cpbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[14]; /* temporary storage */
    cp_R(cpor, tc);

    /*perform dot product */
    for (id = 0; id < 14; ++id) {
        result += x[id]*cpor[id];
    }

    *cpbl = result * 8.31446261815324e+07;
}


/*Returns the mean specific heat at CP (Eq. 34) */
AMREX_GPU_HOST_DEVICE void CKCPBS(double *  T, double *  y,  double *  cpbs)
{
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[14], tresult[14]; /* temporary storage */
    cp_R(cpor, tc);
    for (int i = 0; i < 14; i++)
    {
        tresult[i] = cpor[i]*y[i]*imw[i];

    }
    for (int i = 0; i < 14; i++)
    {
        result += tresult[i];
    }

    *cpbs = result * 8.31446261815324e+07;
}


/*Returns the mean specific heat at CV (Eq. 35) */
void CKCVBL(double *  T, double *  x,  double *  cvbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cvor[14]; /* temporary storage */
    cv_R(cvor, tc);

    /*perform dot product */
    for (id = 0; id < 14; ++id) {
        result += x[id]*cvor[id];
    }

    *cvbl = result * 8.31446261815324e+07;
}


/*Returns the mean specific heat at CV (Eq. 36) */
AMREX_GPU_HOST_DEVICE void CKCVBS(double *  T, double *  y,  double *  cvbs)
{
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cvor[14]; /* temporary storage */
    cv_R(cvor, tc);
    /*multiply by y/molecularweight */
    result += cvor[0]*y[0]*imw[0]; /*H2 */
    result += cvor[1]*y[1]*imw[1]; /*H */
    result += cvor[2]*y[2]*imw[2]; /*AR */
    result += cvor[3]*y[3]*imw[3]; /*N2 */
    result += cvor[4]*y[4]*imw[4]; /*HE */
    result += cvor[5]*y[5]*imw[5]; /*O */
    result += cvor[6]*y[6]*imw[6]; /*OH */
    result += cvor[7]*y[7]*imw[7]; /*HCO */
    result += cvor[8]*y[8]*imw[8]; /*HO2 */
    result += cvor[9]*y[9]*imw[9]; /*H2O */
    result += cvor[10]*y[10]*imw[10]; /*CO */
    result += cvor[11]*y[11]*imw[11]; /*O2 */
    result += cvor[12]*y[12]*imw[12]; /*H2O2 */
    result += cvor[13]*y[13]*imw[13]; /*CO2 */

    *cvbs = result * 8.31446261815324e+07;
}


/*Returns the mean enthalpy of the mixture in molar units */
void CKHBML(double *  T, double *  x,  double *  hbml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[14]; /* temporary storage */
    double RT = 8.31446261815324e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*perform dot product */
    for (id = 0; id < 14; ++id) {
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
    double hml[14], tmp[14]; /* temporary storage */
    double RT = 8.31446261815324e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);
    int id;
    for (id = 0; id < 14; ++id) {
        tmp[id] = y[id]*hml[id]*imw[id];
    }
    for (id = 0; id < 14; ++id) {
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
    double uml[14]; /* temporary energy array */
    double RT = 8.31446261815324e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*perform dot product */
    for (id = 0; id < 14; ++id) {
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
    double ums[14]; /* temporary energy array */
    double RT = 8.31446261815324e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    /*perform dot product + scaling by wt */
    result += y[0]*ums[0]*imw[0]; /*H2 */
    result += y[1]*ums[1]*imw[1]; /*H */
    result += y[2]*ums[2]*imw[2]; /*AR */
    result += y[3]*ums[3]*imw[3]; /*N2 */
    result += y[4]*ums[4]*imw[4]; /*HE */
    result += y[5]*ums[5]*imw[5]; /*O */
    result += y[6]*ums[6]*imw[6]; /*OH */
    result += y[7]*ums[7]*imw[7]; /*HCO */
    result += y[8]*ums[8]*imw[8]; /*HO2 */
    result += y[9]*ums[9]*imw[9]; /*H2O */
    result += y[10]*ums[10]*imw[10]; /*CO */
    result += y[11]*ums[11]*imw[11]; /*O2 */
    result += y[12]*ums[12]*imw[12]; /*H2O2 */
    result += y[13]*ums[13]*imw[13]; /*CO2 */

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
    double sor[14]; /* temporary storage */
    speciesEntropy(sor, tc);

    /*Compute Eq 42 */
    for (id = 0; id < 14; ++id) {
        result += x[id]*(sor[id]-log((x[id]+1e-100))-logPratio);
    }

    *sbml = result * 8.31446261815324e+07;
}


/*get mixture entropy in mass units */
void CKSBMS(double *  P, double *  T, double *  y,  double *  sbms)
{
    double result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    double logPratio = log ( *P / 1013250.0 ); 
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double sor[14]; /* temporary storage */
    double x[14]; /* need a ytx conversion */
    double YOW = 0; /*See Eq 4, 6 in CK Manual */
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
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(1.007970*YOW); 
    x[2] = y[2]/(39.948000*YOW); 
    x[3] = y[3]/(28.013400*YOW); 
    x[4] = y[4]/(4.002600*YOW); 
    x[5] = y[5]/(15.999400*YOW); 
    x[6] = y[6]/(17.007370*YOW); 
    x[7] = y[7]/(29.018520*YOW); 
    x[8] = y[8]/(33.006770*YOW); 
    x[9] = y[9]/(18.015340*YOW); 
    x[10] = y[10]/(28.010550*YOW); 
    x[11] = y[11]/(31.998800*YOW); 
    x[12] = y[12]/(34.014740*YOW); 
    x[13] = y[13]/(44.009950*YOW); 
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
    result += x[13]*(sor[13]-log((x[13]+1e-100))-logPratio);
    /*Scale by R/W */
    *sbms = result * 8.31446261815324e+07 * YOW;
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
    double RT = 8.31446261815324e+07*tT; /*R*T */
    double gort[14]; /* temporary storage */
    /*Compute g/RT */
    gibbs(gort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 14; ++id) {
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
    double RT = 8.31446261815324e+07*tT; /*R*T */
    double gort[14]; /* temporary storage */
    double x[14]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
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
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(1.007970*YOW); 
    x[2] = y[2]/(39.948000*YOW); 
    x[3] = y[3]/(28.013400*YOW); 
    x[4] = y[4]/(4.002600*YOW); 
    x[5] = y[5]/(15.999400*YOW); 
    x[6] = y[6]/(17.007370*YOW); 
    x[7] = y[7]/(29.018520*YOW); 
    x[8] = y[8]/(33.006770*YOW); 
    x[9] = y[9]/(18.015340*YOW); 
    x[10] = y[10]/(28.010550*YOW); 
    x[11] = y[11]/(31.998800*YOW); 
    x[12] = y[12]/(34.014740*YOW); 
    x[13] = y[13]/(44.009950*YOW); 
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
    result += x[13]*(gort[13]+log((x[13]+1e-100))+logPratio);
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
    double RT = 8.31446261815324e+07*tT; /*R*T */
    double aort[14]; /* temporary storage */
    /*Compute g/RT */
    helmholtz(aort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 14; ++id) {
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
    double RT = 8.31446261815324e+07*tT; /*R*T */
    double aort[14]; /* temporary storage */
    double x[14]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
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
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(1.007970*YOW); 
    x[2] = y[2]/(39.948000*YOW); 
    x[3] = y[3]/(28.013400*YOW); 
    x[4] = y[4]/(4.002600*YOW); 
    x[5] = y[5]/(15.999400*YOW); 
    x[6] = y[6]/(17.007370*YOW); 
    x[7] = y[7]/(29.018520*YOW); 
    x[8] = y[8]/(33.006770*YOW); 
    x[9] = y[9]/(18.015340*YOW); 
    x[10] = y[10]/(28.010550*YOW); 
    x[11] = y[11]/(31.998800*YOW); 
    x[12] = y[12]/(34.014740*YOW); 
    x[13] = y[13]/(44.009950*YOW); 
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
    result += x[13]*(aort[13]+log((x[13]+1e-100))+logPratio);
    /*Scale by RT/W */
    *abms = result * RT * YOW;
}


/*compute the production rate for each species */
AMREX_GPU_HOST_DEVICE void CKWC(double *  T, double *  C,  double *  wdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 14; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    productionRate(wdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 14; ++id) {
        C[id] *= 1.0e-6;
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mass fractions */
void CKWYP(double *  P, double *  T, double *  y,  double *  wdot)
{
    int id; /*loop counter */
    double c[14]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
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
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 14; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mole fractions */
void CKWXP(double *  P, double *  T, double *  x,  double *  wdot)
{
    int id; /*loop counter */
    double c[14]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31446261815324e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 14; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 14; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
AMREX_GPU_HOST_DEVICE void CKWYR(double *  rho, double *  T, double *  y,  double *  wdot)
{
    int id; /*loop counter */
    double c[14]; /*temporary storage */
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

    /*call productionRate */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 14; ++id) {
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
    double c[14*(*np)]; /*temporary storage */
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
#endif
}


/*Returns the molar production rate of species */
/*Given rho, T, and mole fractions */
void CKWXR(double *  rho, double *  T, double *  x,  double *  wdot)
{
    int id; /*loop counter */
    double c[14]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
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
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 14; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the rate of progress for each reaction */
void CKQC(double *  T, double *  C, double *  qdot)
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
/*Given P, T, and mole fractions */
void CKKFKR(double *  P, double *  T, double *  x, double *  q_f, double *  q_r)
{
    int id; /*loop counter */
    double c[14]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31446261815324e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

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


/*Returns the progress rates of each reactions */
/*Given P, T, and mass fractions */
void CKQYP(double *  P, double *  T, double *  y, double *  qdot)
{
    int id; /*loop counter */
    double c[14]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
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
void CKQXP(double *  P, double *  T, double *  x, double *  qdot)
{
    int id; /*loop counter */
    double c[14]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31446261815324e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

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
void CKQYR(double *  rho, double *  T, double *  y, double *  qdot)
{
    int id; /*loop counter */
    double c[14]; /*temporary storage */
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
void CKQXR(double *  rho, double *  T, double *  x, double *  qdot)
{
    int id; /*loop counter */
    double c[14]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
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


/*Returns the stoichiometric coefficients */
/*of the reaction mechanism. (Eq 50) */
void CKNU(int * kdim,  int * nuki)
{
    int id; /*loop counter */
    int kd = (*kdim); 
    /*Zero nuki */
    for (id = 0; id < 14 * kd; ++ id) {
         nuki[id] = 0; 
    }

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    nuki[ 1 * kd + 0 ] += -1.000000 ;
    nuki[ 11 * kd + 0 ] += -1.000000 ;
    nuki[ 8 * kd + 0 ] += +1.000000 ;

    /*reaction 2: OH + OH (+M) <=> H2O2 (+M) */
    nuki[ 6 * kd + 1 ] += -1.000000 ;
    nuki[ 6 * kd + 1 ] += -1.000000 ;
    nuki[ 12 * kd + 1 ] += +1.000000 ;

    /*reaction 3: CO + O (+M) <=> CO2 (+M) */
    nuki[ 10 * kd + 2 ] += -1.000000 ;
    nuki[ 5 * kd + 2 ] += -1.000000 ;
    nuki[ 13 * kd + 2 ] += +1.000000 ;

    /*reaction 4: H + H + M <=> H2 + M */
    nuki[ 1 * kd + 3 ] += -1.000000 ;
    nuki[ 1 * kd + 3 ] += -1.000000 ;
    nuki[ 0 * kd + 3 ] += +1.000000 ;

    /*reaction 5: H + OH + M <=> H2O + M */
    nuki[ 1 * kd + 4 ] += -1.000000 ;
    nuki[ 6 * kd + 4 ] += -1.000000 ;
    nuki[ 9 * kd + 4 ] += +1.000000 ;

    /*reaction 6: O + H + M <=> OH + M */
    nuki[ 5 * kd + 5 ] += -1.000000 ;
    nuki[ 1 * kd + 5 ] += -1.000000 ;
    nuki[ 6 * kd + 5 ] += +1.000000 ;

    /*reaction 7: O + O + M <=> O2 + M */
    nuki[ 5 * kd + 6 ] += -1.000000 ;
    nuki[ 5 * kd + 6 ] += -1.000000 ;
    nuki[ 11 * kd + 6 ] += +1.000000 ;

    /*reaction 8: HCO + M <=> CO + H + M */
    nuki[ 7 * kd + 7 ] += -1.000000 ;
    nuki[ 10 * kd + 7 ] += +1.000000 ;
    nuki[ 1 * kd + 7 ] += +1.000000 ;

    /*reaction 9: H + O2 <=> O + OH */
    nuki[ 1 * kd + 8 ] += -1.000000 ;
    nuki[ 11 * kd + 8 ] += -1.000000 ;
    nuki[ 5 * kd + 8 ] += +1.000000 ;
    nuki[ 6 * kd + 8 ] += +1.000000 ;

    /*reaction 10: O + H2 <=> H + OH */
    nuki[ 5 * kd + 9 ] += -1.000000 ;
    nuki[ 0 * kd + 9 ] += -1.000000 ;
    nuki[ 1 * kd + 9 ] += +1.000000 ;
    nuki[ 6 * kd + 9 ] += +1.000000 ;

    /*reaction 11: OH + H2 <=> H + H2O */
    nuki[ 6 * kd + 10 ] += -1.000000 ;
    nuki[ 0 * kd + 10 ] += -1.000000 ;
    nuki[ 1 * kd + 10 ] += +1.000000 ;
    nuki[ 9 * kd + 10 ] += +1.000000 ;

    /*reaction 12: OH + OH <=> O + H2O */
    nuki[ 6 * kd + 11 ] += -1.000000 ;
    nuki[ 6 * kd + 11 ] += -1.000000 ;
    nuki[ 5 * kd + 11 ] += +1.000000 ;
    nuki[ 9 * kd + 11 ] += +1.000000 ;

    /*reaction 13: H + H + H2 <=> H2 + H2 */
    nuki[ 1 * kd + 12 ] += -1.000000 ;
    nuki[ 1 * kd + 12 ] += -1.000000 ;
    nuki[ 0 * kd + 12 ] += -1.000000 ;
    nuki[ 0 * kd + 12 ] += +1.000000 ;
    nuki[ 0 * kd + 12 ] += +1.000000 ;

    /*reaction 14: H + H + H2O <=> H2 + H2O */
    nuki[ 1 * kd + 13 ] += -1.000000 ;
    nuki[ 1 * kd + 13 ] += -1.000000 ;
    nuki[ 9 * kd + 13 ] += -1.000000 ;
    nuki[ 0 * kd + 13 ] += +1.000000 ;
    nuki[ 9 * kd + 13 ] += +1.000000 ;

    /*reaction 15: H + H + CO2 <=> H2 + CO2 */
    nuki[ 1 * kd + 14 ] += -1.000000 ;
    nuki[ 1 * kd + 14 ] += -1.000000 ;
    nuki[ 13 * kd + 14 ] += -1.000000 ;
    nuki[ 0 * kd + 14 ] += +1.000000 ;
    nuki[ 13 * kd + 14 ] += +1.000000 ;

    /*reaction 16: H2 + O2 <=> HO2 + H */
    nuki[ 0 * kd + 15 ] += -1.000000 ;
    nuki[ 11 * kd + 15 ] += -1.000000 ;
    nuki[ 8 * kd + 15 ] += +1.000000 ;
    nuki[ 1 * kd + 15 ] += +1.000000 ;

    /*reaction 17: HO2 + H <=> O + H2O */
    nuki[ 8 * kd + 16 ] += -1.000000 ;
    nuki[ 1 * kd + 16 ] += -1.000000 ;
    nuki[ 5 * kd + 16 ] += +1.000000 ;
    nuki[ 9 * kd + 16 ] += +1.000000 ;

    /*reaction 18: HO2 + H <=> OH + OH */
    nuki[ 8 * kd + 17 ] += -1.000000 ;
    nuki[ 1 * kd + 17 ] += -1.000000 ;
    nuki[ 6 * kd + 17 ] += +1.000000 ;
    nuki[ 6 * kd + 17 ] += +1.000000 ;

    /*reaction 19: HO2 + O <=> OH + O2 */
    nuki[ 8 * kd + 18 ] += -1.000000 ;
    nuki[ 5 * kd + 18 ] += -1.000000 ;
    nuki[ 6 * kd + 18 ] += +1.000000 ;
    nuki[ 11 * kd + 18 ] += +1.000000 ;

    /*reaction 20: HO2 + OH <=> O2 + H2O */
    nuki[ 8 * kd + 19 ] += -1.000000 ;
    nuki[ 6 * kd + 19 ] += -1.000000 ;
    nuki[ 11 * kd + 19 ] += +1.000000 ;
    nuki[ 9 * kd + 19 ] += +1.000000 ;

    /*reaction 21: HO2 + OH <=> O2 + H2O */
    nuki[ 8 * kd + 20 ] += -1.000000 ;
    nuki[ 6 * kd + 20 ] += -1.000000 ;
    nuki[ 11 * kd + 20 ] += +1.000000 ;
    nuki[ 9 * kd + 20 ] += +1.000000 ;

    /*reaction 22: HO2 + HO2 <=> O2 + H2O2 */
    nuki[ 8 * kd + 21 ] += -1.000000 ;
    nuki[ 8 * kd + 21 ] += -1.000000 ;
    nuki[ 11 * kd + 21 ] += +1.000000 ;
    nuki[ 12 * kd + 21 ] += +1.000000 ;

    /*reaction 23: HO2 + HO2 <=> O2 + H2O2 */
    nuki[ 8 * kd + 22 ] += -1.000000 ;
    nuki[ 8 * kd + 22 ] += -1.000000 ;
    nuki[ 11 * kd + 22 ] += +1.000000 ;
    nuki[ 12 * kd + 22 ] += +1.000000 ;

    /*reaction 24: H2O2 + H <=> HO2 + H2 */
    nuki[ 12 * kd + 23 ] += -1.000000 ;
    nuki[ 1 * kd + 23 ] += -1.000000 ;
    nuki[ 8 * kd + 23 ] += +1.000000 ;
    nuki[ 0 * kd + 23 ] += +1.000000 ;

    /*reaction 25: H2O2 + H <=> OH + H2O */
    nuki[ 12 * kd + 24 ] += -1.000000 ;
    nuki[ 1 * kd + 24 ] += -1.000000 ;
    nuki[ 6 * kd + 24 ] += +1.000000 ;
    nuki[ 9 * kd + 24 ] += +1.000000 ;

    /*reaction 26: H2O2 + O <=> OH + HO2 */
    nuki[ 12 * kd + 25 ] += -1.000000 ;
    nuki[ 5 * kd + 25 ] += -1.000000 ;
    nuki[ 6 * kd + 25 ] += +1.000000 ;
    nuki[ 8 * kd + 25 ] += +1.000000 ;

    /*reaction 27: H2O2 + OH <=> HO2 + H2O */
    nuki[ 12 * kd + 26 ] += -1.000000 ;
    nuki[ 6 * kd + 26 ] += -1.000000 ;
    nuki[ 8 * kd + 26 ] += +1.000000 ;
    nuki[ 9 * kd + 26 ] += +1.000000 ;

    /*reaction 28: H2O2 + OH <=> HO2 + H2O */
    nuki[ 12 * kd + 27 ] += -1.000000 ;
    nuki[ 6 * kd + 27 ] += -1.000000 ;
    nuki[ 8 * kd + 27 ] += +1.000000 ;
    nuki[ 9 * kd + 27 ] += +1.000000 ;

    /*reaction 29: CO + OH <=> CO2 + H */
    nuki[ 10 * kd + 28 ] += -1.000000 ;
    nuki[ 6 * kd + 28 ] += -1.000000 ;
    nuki[ 13 * kd + 28 ] += +1.000000 ;
    nuki[ 1 * kd + 28 ] += +1.000000 ;

    /*reaction 30: CO + OH <=> CO2 + H */
    nuki[ 10 * kd + 29 ] += -1.000000 ;
    nuki[ 6 * kd + 29 ] += -1.000000 ;
    nuki[ 13 * kd + 29 ] += +1.000000 ;
    nuki[ 1 * kd + 29 ] += +1.000000 ;

    /*reaction 31: CO + O2 <=> CO2 + O */
    nuki[ 10 * kd + 30 ] += -1.000000 ;
    nuki[ 11 * kd + 30 ] += -1.000000 ;
    nuki[ 13 * kd + 30 ] += +1.000000 ;
    nuki[ 5 * kd + 30 ] += +1.000000 ;

    /*reaction 32: CO + HO2 <=> CO2 + OH */
    nuki[ 10 * kd + 31 ] += -1.000000 ;
    nuki[ 8 * kd + 31 ] += -1.000000 ;
    nuki[ 13 * kd + 31 ] += +1.000000 ;
    nuki[ 6 * kd + 31 ] += +1.000000 ;

    /*reaction 33: HCO + H <=> CO + H2 */
    nuki[ 7 * kd + 32 ] += -1.000000 ;
    nuki[ 1 * kd + 32 ] += -1.000000 ;
    nuki[ 10 * kd + 32 ] += +1.000000 ;
    nuki[ 0 * kd + 32 ] += +1.000000 ;

    /*reaction 34: HCO + O <=> CO + OH */
    nuki[ 7 * kd + 33 ] += -1.000000 ;
    nuki[ 5 * kd + 33 ] += -1.000000 ;
    nuki[ 10 * kd + 33 ] += +1.000000 ;
    nuki[ 6 * kd + 33 ] += +1.000000 ;

    /*reaction 35: HCO + O <=> CO2 + H */
    nuki[ 7 * kd + 34 ] += -1.000000 ;
    nuki[ 5 * kd + 34 ] += -1.000000 ;
    nuki[ 13 * kd + 34 ] += +1.000000 ;
    nuki[ 1 * kd + 34 ] += +1.000000 ;

    /*reaction 36: HCO + OH <=> CO + H2O */
    nuki[ 7 * kd + 35 ] += -1.000000 ;
    nuki[ 6 * kd + 35 ] += -1.000000 ;
    nuki[ 10 * kd + 35 ] += +1.000000 ;
    nuki[ 9 * kd + 35 ] += +1.000000 ;

    /*reaction 37: HCO + H2O <=> CO + H + H2O */
    nuki[ 7 * kd + 36 ] += -1.000000 ;
    nuki[ 9 * kd + 36 ] += -1.000000 ;
    nuki[ 10 * kd + 36 ] += +1.000000 ;
    nuki[ 1 * kd + 36 ] += +1.000000 ;
    nuki[ 9 * kd + 36 ] += +1.000000 ;

    /*reaction 38: HCO + O2 <=> CO + HO2 */
    nuki[ 7 * kd + 37 ] += -1.000000 ;
    nuki[ 11 * kd + 37 ] += -1.000000 ;
    nuki[ 10 * kd + 37 ] += +1.000000 ;
    nuki[ 8 * kd + 37 ] += +1.000000 ;
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
        if (*i > 38) {
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


/*Returns the arrehenius coefficients  */
/*for all reactions */
void CKABE( double *  a, double *  b, double *  e)
{
    // (11):  H + O2 (+M) <=> HO2 (+M)
    a[0] = 5116000000000;
    b[0] = 0.44;
    e[0] = 0;

    // (13):  OH + OH (+M) <=> H2O2 (+M)
    a[1] = 111000000000000;
    b[1] = -0.37;
    e[1] = 0;

    // (26):  CO + O (+M) <=> CO2 (+M)
    a[2] = 13620000000;
    b[2] = 0;
    e[2] = 2384;

    // (4):  H + H + M <=> H2 + M
    a[3] = 1.78e+18;
    b[3] = -1;
    e[3] = 0;

    // (8):  H + OH + M <=> H2O + M
    a[4] = 4.4e+22;
    b[4] = -2;
    e[4] = 0;

    // (9):  O + H + M <=> OH + M
    a[5] = 9.428e+18;
    b[5] = -1;
    e[5] = 0;

    // (10):  O + O + M <=> O2 + M
    a[6] = 1.2e+17;
    b[6] = -1;
    e[6] = 0;

    // (35):  HCO + M <=> CO + H + M
    a[7] = 1.87e+17;
    b[7] = -1;
    e[7] = 17000;

    // (0):  H + O2 <=> O + OH
    a[8] = 26440000000000000;
    b[8] = -0.67069999999999996;
    e[8] = 17041;

    // (1):  O + H2 <=> H + OH
    a[9] = 45890;
    b[9] = 2.7000000000000002;
    e[9] = 6260;

    // (2):  OH + H2 <=> H + H2O
    a[10] = 173400000;
    b[10] = 1.51;
    e[10] = 3430;

    // (3):  OH + OH <=> O + H2O
    a[11] = 39730;
    b[11] = 2.3999999999999999;
    e[11] = -2110;

    // (5):  H + H + H2 <=> H2 + H2
    a[12] = 90000000000000000;
    b[12] = -0.59999999999999998;
    e[12] = 0;

    // (6):  H + H + H2O <=> H2 + H2O
    a[13] = 5.624e+19;
    b[13] = -1.25;
    e[13] = 0;

    // (7):  H + H + CO2 <=> H2 + CO2
    a[14] = 5.5e+20;
    b[14] = -2;
    e[14] = 0;

    // (12):  H2 + O2 <=> HO2 + H
    a[15] = 591600;
    b[15] = 2.4329999999999998;
    e[15] = 53502;

    // (14):  HO2 + H <=> O + H2O
    a[16] = 3970000000000;
    b[16] = 0;
    e[16] = 671;

    // (15):  HO2 + H <=> OH + OH
    a[17] = 74850000000000;
    b[17] = 0;
    e[17] = 295;

    // (16):  HO2 + O <=> OH + O2
    a[18] = 40000000000000;
    b[18] = 0;
    e[18] = 0;

    // (17):  HO2 + OH <=> O2 + H2O
    a[19] = 23750000000000;
    b[19] = 0;
    e[19] = -500;

    // (18):  HO2 + OH <=> O2 + H2O
    a[20] = 10000000000000000;
    b[20] = 0;
    e[20] = 17330;

    // (19):  HO2 + HO2 <=> O2 + H2O2
    a[21] = 130000000000;
    b[21] = 0;
    e[21] = -1630;

    // (20):  HO2 + HO2 <=> O2 + H2O2
    a[22] = 365800000000000;
    b[22] = 0;
    e[22] = 12000;

    // (21):  H2O2 + H <=> HO2 + H2
    a[23] = 6050000;
    b[23] = 2;
    e[23] = 5200;

    // (22):  H2O2 + H <=> OH + H2O
    a[24] = 24100000000000;
    b[24] = 0;
    e[24] = 3970;

    // (23):  H2O2 + O <=> OH + HO2
    a[25] = 9630000;
    b[25] = 2;
    e[25] = 3970;

    // (24):  H2O2 + OH <=> HO2 + H2O
    a[26] = 2000000000000;
    b[26] = 0;
    e[26] = 427;

    // (25):  H2O2 + OH <=> HO2 + H2O
    a[27] = 2.6700000000000001e+41;
    b[27] = -7;
    e[27] = 37600;

    // (27):  CO + OH <=> CO2 + H
    a[28] = 800000000000;
    b[28] = 0.14000000000000001;
    e[28] = 7352;

    // (28):  CO + OH <=> CO2 + H
    a[29] = 87840000000;
    b[29] = 0.029999999999999999;
    e[29] = -16;

    // (29):  CO + O2 <=> CO2 + O
    a[30] = 1119000000000;
    b[30] = 0;
    e[30] = 47700;

    // (30):  CO + HO2 <=> CO2 + OH
    a[31] = 30100000000000;
    b[31] = 0;
    e[31] = 23000;

    // (31):  HCO + H <=> CO + H2
    a[32] = 120000000000000;
    b[32] = 0;
    e[32] = 0;

    // (32):  HCO + O <=> CO + OH
    a[33] = 30000000000000;
    b[33] = 0;
    e[33] = 0;

    // (33):  HCO + O <=> CO2 + H
    a[34] = 30000000000000;
    b[34] = 0;
    e[34] = 0;

    // (34):  HCO + OH <=> CO + H2O
    a[35] = 30200000000000;
    b[35] = 0;
    e[35] = 0;

    // (36):  HCO + H2O <=> CO + H + H2O
    a[36] = 2.244e+18;
    b[36] = -1;
    e[36] = 17000;

    // (37):  HCO + O2 <=> CO + HO2
    a[37] = 12040000000;
    b[37] = 0.80700000000000005;
    e[37] = -727;


    return;
}


/*Returns the equil constants for each reaction */
void CKEQC(double *  T, double *  C, double *  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[14]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: OH + OH (+M) <=> H2O2 (+M) */
    eqcon[1] *= 1e+06; 

    /*reaction 3: CO + O (+M) <=> CO2 (+M) */
    eqcon[2] *= 1e+06; 

    /*reaction 4: H + H + M <=> H2 + M */
    eqcon[3] *= 1e+06; 

    /*reaction 5: H + OH + M <=> H2O + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: O + H + M <=> OH + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: O + O + M <=> O2 + M */
    eqcon[6] *= 1e+06; 

    /*reaction 8: HCO + M <=> CO + H + M */
    eqcon[7] *= 1e-06; 

    /*reaction 9: H + O2 <=> O + OH */
    /*eqcon[8] *= 1;  */

    /*reaction 10: O + H2 <=> H + OH */
    /*eqcon[9] *= 1;  */

    /*reaction 11: OH + H2 <=> H + H2O */
    /*eqcon[10] *= 1;  */

    /*reaction 12: OH + OH <=> O + H2O */
    /*eqcon[11] *= 1;  */

    /*reaction 13: H + H + H2 <=> H2 + H2 */
    eqcon[12] *= 1e+06; 

    /*reaction 14: H + H + H2O <=> H2 + H2O */
    eqcon[13] *= 1e+06; 

    /*reaction 15: H + H + CO2 <=> H2 + CO2 */
    eqcon[14] *= 1e+06; 

    /*reaction 16: H2 + O2 <=> HO2 + H */
    /*eqcon[15] *= 1;  */

    /*reaction 17: HO2 + H <=> O + H2O */
    /*eqcon[16] *= 1;  */

    /*reaction 18: HO2 + H <=> OH + OH */
    /*eqcon[17] *= 1;  */

    /*reaction 19: HO2 + O <=> OH + O2 */
    /*eqcon[18] *= 1;  */

    /*reaction 20: HO2 + OH <=> O2 + H2O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: HO2 + OH <=> O2 + H2O */
    /*eqcon[20] *= 1;  */

    /*reaction 22: HO2 + HO2 <=> O2 + H2O2 */
    /*eqcon[21] *= 1;  */

    /*reaction 23: HO2 + HO2 <=> O2 + H2O2 */
    /*eqcon[22] *= 1;  */

    /*reaction 24: H2O2 + H <=> HO2 + H2 */
    /*eqcon[23] *= 1;  */

    /*reaction 25: H2O2 + H <=> OH + H2O */
    /*eqcon[24] *= 1;  */

    /*reaction 26: H2O2 + O <=> OH + HO2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[26] *= 1;  */

    /*reaction 28: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[27] *= 1;  */

    /*reaction 29: CO + OH <=> CO2 + H */
    /*eqcon[28] *= 1;  */

    /*reaction 30: CO + OH <=> CO2 + H */
    /*eqcon[29] *= 1;  */

    /*reaction 31: CO + O2 <=> CO2 + O */
    /*eqcon[30] *= 1;  */

    /*reaction 32: CO + HO2 <=> CO2 + OH */
    /*eqcon[31] *= 1;  */

    /*reaction 33: HCO + H <=> CO + H2 */
    /*eqcon[32] *= 1;  */

    /*reaction 34: HCO + O <=> CO + OH */
    /*eqcon[33] *= 1;  */

    /*reaction 35: HCO + O <=> CO2 + H */
    /*eqcon[34] *= 1;  */

    /*reaction 36: HCO + OH <=> CO + H2O */
    /*eqcon[35] *= 1;  */

    /*reaction 37: HCO + H2O <=> CO + H + H2O */
    eqcon[36] *= 1e-06; 

    /*reaction 38: HCO + O2 <=> CO + HO2 */
    /*eqcon[37] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mass fractions */
void CKEQYP(double *  P, double *  T, double *  y, double *  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[14]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: OH + OH (+M) <=> H2O2 (+M) */
    eqcon[1] *= 1e+06; 

    /*reaction 3: CO + O (+M) <=> CO2 (+M) */
    eqcon[2] *= 1e+06; 

    /*reaction 4: H + H + M <=> H2 + M */
    eqcon[3] *= 1e+06; 

    /*reaction 5: H + OH + M <=> H2O + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: O + H + M <=> OH + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: O + O + M <=> O2 + M */
    eqcon[6] *= 1e+06; 

    /*reaction 8: HCO + M <=> CO + H + M */
    eqcon[7] *= 1e-06; 

    /*reaction 9: H + O2 <=> O + OH */
    /*eqcon[8] *= 1;  */

    /*reaction 10: O + H2 <=> H + OH */
    /*eqcon[9] *= 1;  */

    /*reaction 11: OH + H2 <=> H + H2O */
    /*eqcon[10] *= 1;  */

    /*reaction 12: OH + OH <=> O + H2O */
    /*eqcon[11] *= 1;  */

    /*reaction 13: H + H + H2 <=> H2 + H2 */
    eqcon[12] *= 1e+06; 

    /*reaction 14: H + H + H2O <=> H2 + H2O */
    eqcon[13] *= 1e+06; 

    /*reaction 15: H + H + CO2 <=> H2 + CO2 */
    eqcon[14] *= 1e+06; 

    /*reaction 16: H2 + O2 <=> HO2 + H */
    /*eqcon[15] *= 1;  */

    /*reaction 17: HO2 + H <=> O + H2O */
    /*eqcon[16] *= 1;  */

    /*reaction 18: HO2 + H <=> OH + OH */
    /*eqcon[17] *= 1;  */

    /*reaction 19: HO2 + O <=> OH + O2 */
    /*eqcon[18] *= 1;  */

    /*reaction 20: HO2 + OH <=> O2 + H2O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: HO2 + OH <=> O2 + H2O */
    /*eqcon[20] *= 1;  */

    /*reaction 22: HO2 + HO2 <=> O2 + H2O2 */
    /*eqcon[21] *= 1;  */

    /*reaction 23: HO2 + HO2 <=> O2 + H2O2 */
    /*eqcon[22] *= 1;  */

    /*reaction 24: H2O2 + H <=> HO2 + H2 */
    /*eqcon[23] *= 1;  */

    /*reaction 25: H2O2 + H <=> OH + H2O */
    /*eqcon[24] *= 1;  */

    /*reaction 26: H2O2 + O <=> OH + HO2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[26] *= 1;  */

    /*reaction 28: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[27] *= 1;  */

    /*reaction 29: CO + OH <=> CO2 + H */
    /*eqcon[28] *= 1;  */

    /*reaction 30: CO + OH <=> CO2 + H */
    /*eqcon[29] *= 1;  */

    /*reaction 31: CO + O2 <=> CO2 + O */
    /*eqcon[30] *= 1;  */

    /*reaction 32: CO + HO2 <=> CO2 + OH */
    /*eqcon[31] *= 1;  */

    /*reaction 33: HCO + H <=> CO + H2 */
    /*eqcon[32] *= 1;  */

    /*reaction 34: HCO + O <=> CO + OH */
    /*eqcon[33] *= 1;  */

    /*reaction 35: HCO + O <=> CO2 + H */
    /*eqcon[34] *= 1;  */

    /*reaction 36: HCO + OH <=> CO + H2O */
    /*eqcon[35] *= 1;  */

    /*reaction 37: HCO + H2O <=> CO + H + H2O */
    eqcon[36] *= 1e-06; 

    /*reaction 38: HCO + O2 <=> CO + HO2 */
    /*eqcon[37] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mole fractions */
void CKEQXP(double *  P, double *  T, double *  x, double *  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[14]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: OH + OH (+M) <=> H2O2 (+M) */
    eqcon[1] *= 1e+06; 

    /*reaction 3: CO + O (+M) <=> CO2 (+M) */
    eqcon[2] *= 1e+06; 

    /*reaction 4: H + H + M <=> H2 + M */
    eqcon[3] *= 1e+06; 

    /*reaction 5: H + OH + M <=> H2O + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: O + H + M <=> OH + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: O + O + M <=> O2 + M */
    eqcon[6] *= 1e+06; 

    /*reaction 8: HCO + M <=> CO + H + M */
    eqcon[7] *= 1e-06; 

    /*reaction 9: H + O2 <=> O + OH */
    /*eqcon[8] *= 1;  */

    /*reaction 10: O + H2 <=> H + OH */
    /*eqcon[9] *= 1;  */

    /*reaction 11: OH + H2 <=> H + H2O */
    /*eqcon[10] *= 1;  */

    /*reaction 12: OH + OH <=> O + H2O */
    /*eqcon[11] *= 1;  */

    /*reaction 13: H + H + H2 <=> H2 + H2 */
    eqcon[12] *= 1e+06; 

    /*reaction 14: H + H + H2O <=> H2 + H2O */
    eqcon[13] *= 1e+06; 

    /*reaction 15: H + H + CO2 <=> H2 + CO2 */
    eqcon[14] *= 1e+06; 

    /*reaction 16: H2 + O2 <=> HO2 + H */
    /*eqcon[15] *= 1;  */

    /*reaction 17: HO2 + H <=> O + H2O */
    /*eqcon[16] *= 1;  */

    /*reaction 18: HO2 + H <=> OH + OH */
    /*eqcon[17] *= 1;  */

    /*reaction 19: HO2 + O <=> OH + O2 */
    /*eqcon[18] *= 1;  */

    /*reaction 20: HO2 + OH <=> O2 + H2O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: HO2 + OH <=> O2 + H2O */
    /*eqcon[20] *= 1;  */

    /*reaction 22: HO2 + HO2 <=> O2 + H2O2 */
    /*eqcon[21] *= 1;  */

    /*reaction 23: HO2 + HO2 <=> O2 + H2O2 */
    /*eqcon[22] *= 1;  */

    /*reaction 24: H2O2 + H <=> HO2 + H2 */
    /*eqcon[23] *= 1;  */

    /*reaction 25: H2O2 + H <=> OH + H2O */
    /*eqcon[24] *= 1;  */

    /*reaction 26: H2O2 + O <=> OH + HO2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[26] *= 1;  */

    /*reaction 28: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[27] *= 1;  */

    /*reaction 29: CO + OH <=> CO2 + H */
    /*eqcon[28] *= 1;  */

    /*reaction 30: CO + OH <=> CO2 + H */
    /*eqcon[29] *= 1;  */

    /*reaction 31: CO + O2 <=> CO2 + O */
    /*eqcon[30] *= 1;  */

    /*reaction 32: CO + HO2 <=> CO2 + OH */
    /*eqcon[31] *= 1;  */

    /*reaction 33: HCO + H <=> CO + H2 */
    /*eqcon[32] *= 1;  */

    /*reaction 34: HCO + O <=> CO + OH */
    /*eqcon[33] *= 1;  */

    /*reaction 35: HCO + O <=> CO2 + H */
    /*eqcon[34] *= 1;  */

    /*reaction 36: HCO + OH <=> CO + H2O */
    /*eqcon[35] *= 1;  */

    /*reaction 37: HCO + H2O <=> CO + H + H2O */
    eqcon[36] *= 1e-06; 

    /*reaction 38: HCO + O2 <=> CO + HO2 */
    /*eqcon[37] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mass fractions */
void CKEQYR(double *  rho, double *  T, double *  y, double *  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[14]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: OH + OH (+M) <=> H2O2 (+M) */
    eqcon[1] *= 1e+06; 

    /*reaction 3: CO + O (+M) <=> CO2 (+M) */
    eqcon[2] *= 1e+06; 

    /*reaction 4: H + H + M <=> H2 + M */
    eqcon[3] *= 1e+06; 

    /*reaction 5: H + OH + M <=> H2O + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: O + H + M <=> OH + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: O + O + M <=> O2 + M */
    eqcon[6] *= 1e+06; 

    /*reaction 8: HCO + M <=> CO + H + M */
    eqcon[7] *= 1e-06; 

    /*reaction 9: H + O2 <=> O + OH */
    /*eqcon[8] *= 1;  */

    /*reaction 10: O + H2 <=> H + OH */
    /*eqcon[9] *= 1;  */

    /*reaction 11: OH + H2 <=> H + H2O */
    /*eqcon[10] *= 1;  */

    /*reaction 12: OH + OH <=> O + H2O */
    /*eqcon[11] *= 1;  */

    /*reaction 13: H + H + H2 <=> H2 + H2 */
    eqcon[12] *= 1e+06; 

    /*reaction 14: H + H + H2O <=> H2 + H2O */
    eqcon[13] *= 1e+06; 

    /*reaction 15: H + H + CO2 <=> H2 + CO2 */
    eqcon[14] *= 1e+06; 

    /*reaction 16: H2 + O2 <=> HO2 + H */
    /*eqcon[15] *= 1;  */

    /*reaction 17: HO2 + H <=> O + H2O */
    /*eqcon[16] *= 1;  */

    /*reaction 18: HO2 + H <=> OH + OH */
    /*eqcon[17] *= 1;  */

    /*reaction 19: HO2 + O <=> OH + O2 */
    /*eqcon[18] *= 1;  */

    /*reaction 20: HO2 + OH <=> O2 + H2O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: HO2 + OH <=> O2 + H2O */
    /*eqcon[20] *= 1;  */

    /*reaction 22: HO2 + HO2 <=> O2 + H2O2 */
    /*eqcon[21] *= 1;  */

    /*reaction 23: HO2 + HO2 <=> O2 + H2O2 */
    /*eqcon[22] *= 1;  */

    /*reaction 24: H2O2 + H <=> HO2 + H2 */
    /*eqcon[23] *= 1;  */

    /*reaction 25: H2O2 + H <=> OH + H2O */
    /*eqcon[24] *= 1;  */

    /*reaction 26: H2O2 + O <=> OH + HO2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[26] *= 1;  */

    /*reaction 28: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[27] *= 1;  */

    /*reaction 29: CO + OH <=> CO2 + H */
    /*eqcon[28] *= 1;  */

    /*reaction 30: CO + OH <=> CO2 + H */
    /*eqcon[29] *= 1;  */

    /*reaction 31: CO + O2 <=> CO2 + O */
    /*eqcon[30] *= 1;  */

    /*reaction 32: CO + HO2 <=> CO2 + OH */
    /*eqcon[31] *= 1;  */

    /*reaction 33: HCO + H <=> CO + H2 */
    /*eqcon[32] *= 1;  */

    /*reaction 34: HCO + O <=> CO + OH */
    /*eqcon[33] *= 1;  */

    /*reaction 35: HCO + O <=> CO2 + H */
    /*eqcon[34] *= 1;  */

    /*reaction 36: HCO + OH <=> CO + H2O */
    /*eqcon[35] *= 1;  */

    /*reaction 37: HCO + H2O <=> CO + H + H2O */
    eqcon[36] *= 1e-06; 

    /*reaction 38: HCO + O2 <=> CO + HO2 */
    /*eqcon[37] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mole fractions */
void CKEQXR(double *  rho, double *  T, double *  x, double *  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[14]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: OH + OH (+M) <=> H2O2 (+M) */
    eqcon[1] *= 1e+06; 

    /*reaction 3: CO + O (+M) <=> CO2 (+M) */
    eqcon[2] *= 1e+06; 

    /*reaction 4: H + H + M <=> H2 + M */
    eqcon[3] *= 1e+06; 

    /*reaction 5: H + OH + M <=> H2O + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: O + H + M <=> OH + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: O + O + M <=> O2 + M */
    eqcon[6] *= 1e+06; 

    /*reaction 8: HCO + M <=> CO + H + M */
    eqcon[7] *= 1e-06; 

    /*reaction 9: H + O2 <=> O + OH */
    /*eqcon[8] *= 1;  */

    /*reaction 10: O + H2 <=> H + OH */
    /*eqcon[9] *= 1;  */

    /*reaction 11: OH + H2 <=> H + H2O */
    /*eqcon[10] *= 1;  */

    /*reaction 12: OH + OH <=> O + H2O */
    /*eqcon[11] *= 1;  */

    /*reaction 13: H + H + H2 <=> H2 + H2 */
    eqcon[12] *= 1e+06; 

    /*reaction 14: H + H + H2O <=> H2 + H2O */
    eqcon[13] *= 1e+06; 

    /*reaction 15: H + H + CO2 <=> H2 + CO2 */
    eqcon[14] *= 1e+06; 

    /*reaction 16: H2 + O2 <=> HO2 + H */
    /*eqcon[15] *= 1;  */

    /*reaction 17: HO2 + H <=> O + H2O */
    /*eqcon[16] *= 1;  */

    /*reaction 18: HO2 + H <=> OH + OH */
    /*eqcon[17] *= 1;  */

    /*reaction 19: HO2 + O <=> OH + O2 */
    /*eqcon[18] *= 1;  */

    /*reaction 20: HO2 + OH <=> O2 + H2O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: HO2 + OH <=> O2 + H2O */
    /*eqcon[20] *= 1;  */

    /*reaction 22: HO2 + HO2 <=> O2 + H2O2 */
    /*eqcon[21] *= 1;  */

    /*reaction 23: HO2 + HO2 <=> O2 + H2O2 */
    /*eqcon[22] *= 1;  */

    /*reaction 24: H2O2 + H <=> HO2 + H2 */
    /*eqcon[23] *= 1;  */

    /*reaction 25: H2O2 + H <=> OH + H2O */
    /*eqcon[24] *= 1;  */

    /*reaction 26: H2O2 + O <=> OH + HO2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[26] *= 1;  */

    /*reaction 28: H2O2 + OH <=> HO2 + H2O */
    /*eqcon[27] *= 1;  */

    /*reaction 29: CO + OH <=> CO2 + H */
    /*eqcon[28] *= 1;  */

    /*reaction 30: CO + OH <=> CO2 + H */
    /*eqcon[29] *= 1;  */

    /*reaction 31: CO + O2 <=> CO2 + O */
    /*eqcon[30] *= 1;  */

    /*reaction 32: CO + HO2 <=> CO2 + OH */
    /*eqcon[31] *= 1;  */

    /*reaction 33: HCO + H <=> CO + H2 */
    /*eqcon[32] *= 1;  */

    /*reaction 34: HCO + O <=> CO + OH */
    /*eqcon[33] *= 1;  */

    /*reaction 35: HCO + O <=> CO2 + H */
    /*eqcon[34] *= 1;  */

    /*reaction 36: HCO + OH <=> CO + H2O */
    /*eqcon[35] *= 1;  */

    /*reaction 37: HCO + H2O <=> CO + H + H2O */
    eqcon[36] *= 1e-06; 

    /*reaction 38: HCO + O2 <=> CO + HO2 */
    /*eqcon[37] *= 1;  */
}

#ifdef AMREX_USE_CUDA
/*GPU version of productionRate: no more use of thermo namespace vectors */
/*compute the production rate for each species */
AMREX_GPU_HOST_DEVICE inline void  productionRate(double * wdot, double * sc, double T)
{
    double tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    double invT = 1.0 / tc[1];

    double qdot, q_f[38], q_r[38];
    comp_qfqr(q_f, q_r, sc, tc, invT);

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

AMREX_GPU_HOST_DEVICE inline void comp_qfqr(double *  qf, double * qr, double * sc, double * tc, double invT)
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

    /*compute the mixture concentration */
    double mixture = 0.0;
    for (int i = 0; i < 14; ++i) {
        mixture += sc[i];
    }

    /*compute the Gibbs free energy */
    double g_RT[14];
    gibbs(g_RT, tc);

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31446 * invT;
    double refCinv = 1 / refC;

    /* Evaluate the kfs */
    double k_f, k_r, Corr;
    double redP, F, logPred, logFcent, troe_c, troe_n, troe, F_troe;

    // (0):  H + O2 <=> O + OH
    k_f = 1.0000000000000002e-06 * 26440000000000000 
               * exp(-0.67069999999999996 * tc[0] - 0.50321666580471969 * (17041) * invT);
    Corr  = 1.0;
    qf[8] *= Corr * k_f;
    qr[8] *= Corr * k_f / exp(g_RT[1] - g_RT[5] - g_RT[6] + g_RT[11]);
    // (1):  O + H2 <=> H + OH
    k_f = 1.0000000000000002e-06 * 45890 
               * exp(2.7000000000000002 * tc[0] - 0.50321666580471969 * (6260) * invT);
    Corr  = 1.0;
    qf[9] *= Corr * k_f;
    qr[9] *= Corr * k_f / exp(g_RT[0] - g_RT[1] + g_RT[5] - g_RT[6]);
    // (2):  OH + H2 <=> H + H2O
    k_f = 1.0000000000000002e-06 * 173400000 
               * exp(1.51 * tc[0] - 0.50321666580471969 * (3430) * invT);
    Corr  = 1.0;
    qf[10] *= Corr * k_f;
    qr[10] *= Corr * k_f / exp(g_RT[0] - g_RT[1] + g_RT[6] - g_RT[9]);
    // (3):  OH + OH <=> O + H2O
    k_f = 1.0000000000000002e-06 * 39730 
               * exp(2.3999999999999999 * tc[0] - 0.50321666580471969 * (-2110) * invT);
    Corr  = 1.0;
    qf[11] *= Corr * k_f;
    qr[11] *= Corr * k_f / exp(-g_RT[5] + g_RT[6] + g_RT[6] - g_RT[9]);
    // (4):  H + H + M <=> H2 + M
    k_f = 1.0000000000000002e-12 * 1.78e+18 
               * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = mixture + ( 0 - 1)*sc[0] + ( 0 - 1)*sc[9] + ( 0 - 1)*sc[13] + ( 0.63 - 1)*sc[2] + ( 0.63 - 1)*sc[4];
    qf[3] *= Corr * k_f;
    qr[3] *= Corr * k_f / (exp(-g_RT[0] + g_RT[1] + g_RT[1]) * refCinv);
    // (5):  H + H + H2 <=> H2 + H2
    k_f = 1.0000000000000002e-12 * 90000000000000000 
               * exp(-0.59999999999999998 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[12] *= Corr * k_f;
    qr[12] *= Corr * k_f / (exp(g_RT[0] - g_RT[0] - g_RT[0] + g_RT[1] + g_RT[1]) * refCinv);
    // (6):  H + H + H2O <=> H2 + H2O
    k_f = 1.0000000000000002e-12 * 5.624e+19 
               * exp(-1.25 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[13] *= Corr * k_f;
    qr[13] *= Corr * k_f / (exp(-g_RT[0] + g_RT[1] + g_RT[1] + g_RT[9] - g_RT[9]) * refCinv);
    // (7):  H + H + CO2 <=> H2 + CO2
    k_f = 1.0000000000000002e-12 * 5.5e+20 
               * exp(-2 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[14] *= Corr * k_f;
    qr[14] *= Corr * k_f / (exp(-g_RT[0] + g_RT[1] + g_RT[1] + g_RT[13] - g_RT[13]) * refCinv);
    // (8):  H + OH + M <=> H2O + M
    k_f = 1.0000000000000002e-12 * 4.4e+22 
               * exp(-2 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = mixture + ( 2 - 1)*sc[0] + ( 6.2999999999999998 - 1)*sc[9] + ( 1.75 - 1)*sc[10] + ( 3.6000000000000001 - 1)*sc[13] + ( 0.38 - 1)*sc[2] + ( 0.38 - 1)*sc[4];
    qf[4] *= Corr * k_f;
    qr[4] *= Corr * k_f / (exp(g_RT[1] + g_RT[6] - g_RT[9]) * refCinv);
    // (9):  O + H + M <=> OH + M
    k_f = 1.0000000000000002e-12 * 9.428e+18 
               * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = mixture + ( 2 - 1)*sc[0] + ( 12 - 1)*sc[9] + ( 1.75 - 1)*sc[10] + ( 3.6000000000000001 - 1)*sc[13] + ( 0.69999999999999996 - 1)*sc[2] + ( 0.69999999999999996 - 1)*sc[4];
    qf[5] *= Corr * k_f;
    qr[5] *= Corr * k_f / (exp(g_RT[1] + g_RT[5] - g_RT[6]) * refCinv);
    // (10):  O + O + M <=> O2 + M
    k_f = 1.0000000000000002e-12 * 1.2e+17 
               * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = mixture + ( 2.3999999999999999 - 1)*sc[0] + ( 15.4 - 1)*sc[9] + ( 1.75 - 1)*sc[10] + ( 3.6000000000000001 - 1)*sc[13] + ( 0.82999999999999996 - 1)*sc[2] + ( 0.82999999999999996 - 1)*sc[4];
    qf[6] *= Corr * k_f;
    qr[6] *= Corr * k_f / (exp(g_RT[5] + g_RT[5] - g_RT[11]) * refCinv);
    // (11):  H + O2 (+M) <=> HO2 (+M)
    k_f = 1.0000000000000002e-06 * 5116000000000 
               * exp(0.44 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = mixture + ( 0.84999999999999998 - 1)*sc[11] + ( 11.890000000000001 - 1)*sc[9] + ( 1.0900000000000001 - 1)*sc[10] + ( 2.1800000000000002 - 1)*sc[13] + ( 0.40000000000000002 - 1)*sc[2] + ( 0.46000000000000002 - 1)*sc[4] + ( 0.75 - 1)*sc[0];
    redP = Corr / k_f * 1e-12 * 6.328e+19 
               * exp(-1.3999999999999999  * tc[0] - 0.50321666580471969  * (0) *invT);
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
    qr[0] *= Corr * k_f / (exp(g_RT[1] - g_RT[8] + g_RT[11]) * refCinv);
    // (12):  H2 + O2 <=> HO2 + H
    k_f = 1.0000000000000002e-06 * 591600 
               * exp(2.4329999999999998 * tc[0] - 0.50321666580471969 * (53502) * invT);
    Corr  = 1.0;
    qf[15] *= Corr * k_f;
    qr[15] *= Corr * k_f / exp(g_RT[0] - g_RT[1] - g_RT[8] + g_RT[11]);
    // (13):  OH + OH (+M) <=> H2O2 (+M)
    k_f = 1.0000000000000002e-06 * 111000000000000 
               * exp(-0.37 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[9] + ( 1.75 - 1)*sc[10] + ( 3.6000000000000001 - 1)*sc[13] + ( 0.69999999999999996 - 1)*sc[2] + ( 0.69999999999999996 - 1)*sc[4];
    redP = Corr / k_f * 1e-12 * 2.01e+17 
               * exp(-0.58399999999999996  * tc[0] - 0.50321666580471969  * (-2293) *invT);
    F = redP / (1.0 + redP);
    logPred = log10(redP);
    logFcent = log10(
        (1.-0.73460000000000003)*exp(-tc[1] / 94) 
        + 0.73460000000000003 * exp(-tc[1]/1756)  
        + exp(-5182 * invT));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10., logFcent / (1.0 + troe*troe));
    Corr = F * F_troe;
    qf[1] *= Corr * k_f;
    qr[1] *= Corr * k_f / (exp(g_RT[6] + g_RT[6] - g_RT[12]) * refCinv);
    // (14):  HO2 + H <=> O + H2O
    k_f = 1.0000000000000002e-06 * 3970000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (671) * invT);
    Corr  = 1.0;
    qf[16] *= Corr * k_f;
    qr[16] *= Corr * k_f / exp(g_RT[1] - g_RT[5] + g_RT[8] - g_RT[9]);
    // (15):  HO2 + H <=> OH + OH
    k_f = 1.0000000000000002e-06 * 74850000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (295) * invT);
    Corr  = 1.0;
    qf[17] *= Corr * k_f;
    qr[17] *= Corr * k_f / exp(g_RT[1] - g_RT[6] - g_RT[6] + g_RT[8]);
    // (16):  HO2 + O <=> OH + O2
    k_f = 1.0000000000000002e-06 * 40000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[18] *= Corr * k_f;
    qr[18] *= Corr * k_f / exp(g_RT[5] - g_RT[6] + g_RT[8] - g_RT[11]);
    // (17):  HO2 + OH <=> O2 + H2O
    k_f = 1.0000000000000002e-06 * 23750000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (-500) * invT);
    Corr  = 1.0;
    qf[19] *= Corr * k_f;
    qr[19] *= Corr * k_f / exp(g_RT[6] + g_RT[8] - g_RT[9] - g_RT[11]);
    // (18):  HO2 + OH <=> O2 + H2O
    k_f = 1.0000000000000002e-06 * 10000000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (17330) * invT);
    Corr  = 1.0;
    qf[20] *= Corr * k_f;
    qr[20] *= Corr * k_f / exp(g_RT[6] + g_RT[8] - g_RT[9] - g_RT[11]);
    // (19):  HO2 + HO2 <=> O2 + H2O2
    k_f = 1.0000000000000002e-06 * 130000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (-1630) * invT);
    Corr  = 1.0;
    qf[21] *= Corr * k_f;
    qr[21] *= Corr * k_f / exp(g_RT[8] + g_RT[8] - g_RT[11] - g_RT[12]);
    // (20):  HO2 + HO2 <=> O2 + H2O2
    k_f = 1.0000000000000002e-06 * 365800000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (12000) * invT);
    Corr  = 1.0;
    qf[22] *= Corr * k_f;
    qr[22] *= Corr * k_f / exp(g_RT[8] + g_RT[8] - g_RT[11] - g_RT[12]);
    // (21):  H2O2 + H <=> HO2 + H2
    k_f = 1.0000000000000002e-06 * 6050000 
               * exp(2 * tc[0] - 0.50321666580471969 * (5200) * invT);
    Corr  = 1.0;
    qf[23] *= Corr * k_f;
    qr[23] *= Corr * k_f / exp(-g_RT[0] + g_RT[1] - g_RT[8] + g_RT[12]);
    // (22):  H2O2 + H <=> OH + H2O
    k_f = 1.0000000000000002e-06 * 24100000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (3970) * invT);
    Corr  = 1.0;
    qf[24] *= Corr * k_f;
    qr[24] *= Corr * k_f / exp(g_RT[1] - g_RT[6] - g_RT[9] + g_RT[12]);
    // (23):  H2O2 + O <=> OH + HO2
    k_f = 1.0000000000000002e-06 * 9630000 
               * exp(2 * tc[0] - 0.50321666580471969 * (3970) * invT);
    Corr  = 1.0;
    qf[25] *= Corr * k_f;
    qr[25] *= Corr * k_f / exp(g_RT[5] - g_RT[6] - g_RT[8] + g_RT[12]);
    // (24):  H2O2 + OH <=> HO2 + H2O
    k_f = 1.0000000000000002e-06 * 2000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (427) * invT);
    Corr  = 1.0;
    qf[26] *= Corr * k_f;
    qr[26] *= Corr * k_f / exp(g_RT[6] - g_RT[8] - g_RT[9] + g_RT[12]);
    // (25):  H2O2 + OH <=> HO2 + H2O
    k_f = 1.0000000000000002e-06 * 2.6700000000000001e+41 
               * exp(-7 * tc[0] - 0.50321666580471969 * (37600) * invT);
    Corr  = 1.0;
    qf[27] *= Corr * k_f;
    qr[27] *= Corr * k_f / exp(g_RT[6] - g_RT[8] - g_RT[9] + g_RT[12]);
    // (26):  CO + O (+M) <=> CO2 (+M)
    k_f = 1.0000000000000002e-06 * 13620000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (2384) * invT);
    Corr  = mixture + ( 2 - 1)*sc[0] + ( 12 - 1)*sc[9] + ( 1.75 - 1)*sc[10] + ( 3.6000000000000001 - 1)*sc[13] + ( 0.69999999999999996 - 1)*sc[2] + ( 0.69999999999999996 - 1)*sc[4];
    redP = Corr / k_f * 1e-12 * 1.1729999999999999e+24 
               * exp(-2.79  * tc[0] - 0.50321666580471969  * (4191) *invT);
    Corr = redP / (1. + redP);
    qf[2] *= Corr * k_f;
    qr[2] *= Corr * k_f / (exp(g_RT[5] + g_RT[10] - g_RT[13]) * refCinv);
    // (27):  CO + OH <=> CO2 + H
    k_f = 1.0000000000000002e-06 * 800000000000 
               * exp(0.14000000000000001 * tc[0] - 0.50321666580471969 * (7352) * invT);
    Corr  = 1.0;
    qf[28] *= Corr * k_f;
    qr[28] *= Corr * k_f / exp(-g_RT[1] + g_RT[6] + g_RT[10] - g_RT[13]);
    // (28):  CO + OH <=> CO2 + H
    k_f = 1.0000000000000002e-06 * 87840000000 
               * exp(0.029999999999999999 * tc[0] - 0.50321666580471969 * (-16) * invT);
    Corr  = 1.0;
    qf[29] *= Corr * k_f;
    qr[29] *= Corr * k_f / exp(-g_RT[1] + g_RT[6] + g_RT[10] - g_RT[13]);
    // (29):  CO + O2 <=> CO2 + O
    k_f = 1.0000000000000002e-06 * 1119000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (47700) * invT);
    Corr  = 1.0;
    qf[30] *= Corr * k_f;
    qr[30] *= Corr * k_f / exp(-g_RT[5] + g_RT[10] + g_RT[11] - g_RT[13]);
    // (30):  CO + HO2 <=> CO2 + OH
    k_f = 1.0000000000000002e-06 * 30100000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (23000) * invT);
    Corr  = 1.0;
    qf[31] *= Corr * k_f;
    qr[31] *= Corr * k_f / exp(-g_RT[6] + g_RT[8] + g_RT[10] - g_RT[13]);
    // (31):  HCO + H <=> CO + H2
    k_f = 1.0000000000000002e-06 * 120000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[32] *= Corr * k_f;
    qr[32] *= Corr * k_f / exp(-g_RT[0] + g_RT[1] + g_RT[7] - g_RT[10]);
    // (32):  HCO + O <=> CO + OH
    k_f = 1.0000000000000002e-06 * 30000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[33] *= Corr * k_f;
    qr[33] *= Corr * k_f / exp(g_RT[5] - g_RT[6] + g_RT[7] - g_RT[10]);
    // (33):  HCO + O <=> CO2 + H
    k_f = 1.0000000000000002e-06 * 30000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[34] *= Corr * k_f;
    qr[34] *= Corr * k_f / exp(-g_RT[1] + g_RT[5] + g_RT[7] - g_RT[13]);
    // (34):  HCO + OH <=> CO + H2O
    k_f = 1.0000000000000002e-06 * 30200000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[35] *= Corr * k_f;
    qr[35] *= Corr * k_f / exp(g_RT[6] + g_RT[7] - g_RT[9] - g_RT[10]);
    // (35):  HCO + M <=> CO + H + M
    k_f = 1.0000000000000002e-06 * 1.87e+17 
               * exp(-1 * tc[0] - 0.50321666580471969 * (17000) * invT);
    Corr  = mixture + ( 2 - 1)*sc[0] + ( 0 - 1)*sc[9] + ( 1.75 - 1)*sc[10] + ( 3.6000000000000001 - 1)*sc[13];
    qf[7] *= Corr * k_f;
    qr[7] *= Corr * k_f / (exp(-g_RT[1] + g_RT[7] - g_RT[10]) * refC);
    // (36):  HCO + H2O <=> CO + H + H2O
    k_f = 1.0000000000000002e-06 * 2.244e+18 
               * exp(-1 * tc[0] - 0.50321666580471969 * (17000) * invT);
    Corr  = 1.0;
    qf[36] *= Corr * k_f;
    qr[36] *= Corr * k_f / (exp(-g_RT[1] + g_RT[7] + g_RT[9] - g_RT[9] - g_RT[10]) * refC);
    // (37):  HCO + O2 <=> CO + HO2
    k_f = 1.0000000000000002e-06 * 12040000000 
               * exp(0.80700000000000005 * tc[0] - 0.50321666580471969 * (-727) * invT);
    Corr  = 1.0;
    qf[37] *= Corr * k_f;
    qr[37] *= Corr * k_f / exp(g_RT[7] - g_RT[8] - g_RT[10] + g_RT[11]);


    return;
}
#endif


#ifndef AMREX_USE_CUDA
static double T_save = -1;
#ifdef _OPENMP
#pragma omp threadprivate(T_save)
#endif

static double k_f_save[38];
#ifdef _OPENMP
#pragma omp threadprivate(k_f_save)
#endif

static double Kc_save[38];
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

    double qdot, q_f[38], q_r[38];
    comp_qfqr(q_f, q_r, sc, tc, invT);

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

void comp_k_f(double *  tc, double invT, double *  k_f)
{
    for (int i=0; i<38; ++i) {
        k_f[i] = prefactor_units[i] * fwd_A[i]
                    * exp(fwd_beta[i] * tc[0] - activation_units[i] * fwd_Ea[i] * invT);
    };
    return;
}

void comp_Kc(double *  tc, double invT, double *  Kc)
{
    /*compute the Gibbs free energy */
    double g_RT[14];
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
    double refC = 101325 / 8.31446 * invT;
    double refCinv = 1 / refC;

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

void comp_qfqr(double *  qf, double *  qr, double *  sc, double *  tc, double invT)
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

    double T = tc[1];

    /*compute the mixture concentration */
    double mixture = 0.0;
    for (int i = 0; i < 14; ++i) {
        mixture += sc[i];
    }

    double Corr[38];
    for (int i = 0; i < 38; ++i) {
        Corr[i] = 1.0;
    }

    /* troe */
    {
        double alpha[2];
        alpha[0] = mixture + (TB[0][0] - 1)*sc[11] + (TB[0][1] - 1)*sc[9] + (TB[0][2] - 1)*sc[10] + (TB[0][3] - 1)*sc[13] + (TB[0][4] - 1)*sc[2] + (TB[0][5] - 1)*sc[4] + (TB[0][6] - 1)*sc[0];
        alpha[1] = mixture + (TB[1][0] - 1)*sc[0] + (TB[1][1] - 1)*sc[9] + (TB[1][2] - 1)*sc[10] + (TB[1][3] - 1)*sc[13] + (TB[1][4] - 1)*sc[2] + (TB[1][5] - 1)*sc[4];
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

    /* Lindemann */
    {
        double alpha;
        alpha = mixture + (TB[2][0] - 1)*sc[0] + (TB[2][1] - 1)*sc[9] + (TB[2][2] - 1)*sc[10] + (TB[2][3] - 1)*sc[13] + (TB[2][4] - 1)*sc[2] + (TB[2][5] - 1)*sc[4];
        double redP = alpha / k_f_save[2] * phase_units[2] * low_A[2] * exp(low_beta[2] * tc[0] - activation_units[2] * low_Ea[2] * invT);
        Corr[2] = redP / (1. + redP);
    }

    /* simple three-body correction */
    {
        double alpha;
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
#endif


#ifndef AMREX_USE_CUDA
/*compute the production rate for each species */
void vproductionRate(int npt, double *  wdot, double *  sc, double *  T)
{
    double k_f_s[38*npt], Kc_s[38*npt], mixture[npt], g_RT[14*npt];
    double tc[5*npt], invT[npt];

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

void vcomp_k_f(int npt, double *  k_f_s, double *  tc, double *  invT)
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

void vcomp_gibbs(int npt, double *  g_RT, double *  tc)
{
    /*compute the Gibbs free energy */
    for (int i=0; i<npt; i++) {
        double tg[5], g[14];
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

void vcomp_Kc(int npt, double *  Kc_s, double *  g_RT, double *  invT)
{
    for (int i=0; i<npt; i++) {
        /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
        double refC = (101325. / 8.31451) * invT[i];
        double refCinv = 1.0 / refC;

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

void vcomp_wdot(int npt, double *  wdot, double *  mixture, double *  sc,
		double *  k_f_s, double *  Kc_s,
		double *  tc, double *  invT, double *  T)
{
    for (int i=0; i<npt; i++) {
        double qdot, q_f, q_r, phi_f, phi_r, k_f, k_r, Kc;
        double alpha;
        double redP, F;
        double logPred;
        double logFcent, troe_c, troe_n, troe, F_troe;

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
#endif

/*compute an approx to the reaction Jacobian (for preconditioning) */
AMREX_GPU_HOST_DEVICE void DWDOT_SIMPLIFIED(double *  J, double *  sc, double *  Tp, int * HP)
{
    double c[14];

    for (int k=0; k<14; k++) {
        c[k] = 1.e6 * sc[k];
    }

    aJacobian_precond(J, c, *Tp, *HP);

    /* dwdot[k]/dT */
    /* dTdot/d[X] */
    for (int k=0; k<14; k++) {
        J[210+k] *= 1.e-6;
        J[k*15+14] *= 1.e6;
    }

    return;
}

/*compute the reaction Jacobian */
AMREX_GPU_HOST_DEVICE void DWDOT(double *  J, double *  sc, double *  Tp, int * consP)
{
    double c[14];

    for (int k=0; k<14; k++) {
        c[k] = 1.e6 * sc[k];
    }

    aJacobian(J, c, *Tp, *consP);

    /* dwdot[k]/dT */
    /* dTdot/d[X] */
    for (int k=0; k<14; k++) {
        J[210+k] *= 1.e-6;
        J[k*15+14] *= 1.e6;
    }

    return;
}

/*compute the sparsity pattern of the chemistry Jacobian */
AMREX_GPU_HOST_DEVICE void SPARSITY_INFO( int * nJdata, int * consP, int NCELLS)
{
    double c[14];
    double J[225];

    for (int k=0; k<14; k++) {
        c[k] = 1.0/ 14.000000 ;
    }

    aJacobian(J, c, 1500.0, *consP);

    int nJdata_tmp = 0;
    for (int k=0; k<15; k++) {
        for (int l=0; l<15; l++) {
            if(J[ 15 * k + l] != 0.0){
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
    double c[14];
    double J[225];

    for (int k=0; k<14; k++) {
        c[k] = 1.0/ 14.000000 ;
    }

    aJacobian(J, c, 1500.0, *consP);

    int nJdata_tmp = 0;
    for (int k=0; k<15; k++) {
        for (int l=0; l<15; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J[ 15 * k + l] != 0.0){
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
    double c[14];
    double J[225];

    for (int k=0; k<14; k++) {
        c[k] = 1.0/ 14.000000 ;
    }

    aJacobian_precond(J, c, 1500.0, *consP);

    int nJdata_tmp = 0;
    for (int k=0; k<15; k++) {
        for (int l=0; l<15; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J[ 15 * k + l] != 0.0){
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
    double c[14];
    double J[225];
    int offset_row;
    int offset_col;

    for (int k=0; k<14; k++) {
        c[k] = 1.0/ 14.000000 ;
    }

    aJacobian(J, c, 1500.0, *consP);

    colPtrs[0] = 0;
    int nJdata_tmp = 0;
    for (int nc=0; nc<NCELLS; nc++) {
        offset_row = nc * 15;
        offset_col = nc * 15;
        for (int k=0; k<15; k++) {
            for (int l=0; l<15; l++) {
                if(J[15*k + l] != 0.0) {
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
    double c[14];
    double J[225];
    int offset;

    for (int k=0; k<14; k++) {
        c[k] = 1.0/ 14.000000 ;
    }

    aJacobian(J, c, 1500.0, *consP);

    if (base == 1) {
        rowPtrs[0] = 1;
        int nJdata_tmp = 1;
        for (int nc=0; nc<NCELLS; nc++) {
            offset = nc * 15;
            for (int l=0; l<15; l++) {
                for (int k=0; k<15; k++) {
                    if(J[15*k + l] != 0.0) {
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
                    if(J[15*k + l] != 0.0) {
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
    double c[14];
    double J[225];
    int offset;

    for (int k=0; k<14; k++) {
        c[k] = 1.0/ 14.000000 ;
    }

    aJacobian(J, c, 1500.0, *consP);

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
                        if(J[15*k + l] != 0.0) {
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
                        if(J[15*k + l] != 0.0) {
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
    double c[14];
    double J[225];

    for (int k=0; k<14; k++) {
        c[k] = 1.0/ 14.000000 ;
    }

    aJacobian_precond(J, c, 1500.0, *consP);

    colPtrs[0] = 0;
    int nJdata_tmp = 0;
    for (int k=0; k<15; k++) {
        for (int l=0; l<15; l++) {
            if (k == l) {
                rowVals[nJdata_tmp] = l; 
                indx[nJdata_tmp] = 15*k + l;
                nJdata_tmp = nJdata_tmp + 1; 
            } else {
                if(J[15*k + l] != 0.0) {
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
AMREX_GPU_HOST_DEVICE void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(int * colVals, int * rowPtr, int * consP, int base)
{
    double c[14];
    double J[225];

    for (int k=0; k<14; k++) {
        c[k] = 1.0/ 14.000000 ;
    }

    aJacobian_precond(J, c, 1500.0, *consP);

    if (base == 1) {
        rowPtr[0] = 1;
        int nJdata_tmp = 1;
        for (int l=0; l<15; l++) {
            for (int k=0; k<15; k++) {
                if (k == l) {
                    colVals[nJdata_tmp-1] = l+1; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J[15*k + l] != 0.0) {
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
                    if(J[15*k + l] != 0.0) {
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


    for (int i=0; i<225; i++) {
        J[i] = 0.0;
    }

    double wdot[14];
    for (int k=0; k<14; k++) {
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
    for (int k = 0; k < 14; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    double g_RT[14];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    double h_RT[14];
    speciesEnthalpy(h_RT, tc);

    double phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    double dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    double dqdci, dcdc_fac, dqdc[14];
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
    alpha = mixture + ( 0.84999999999999998 - 1)*sc[11] + ( 11.890000000000001 - 1)*sc[9] + ( 1.0900000000000001 - 1)*sc[10] + ( 2.1800000000000002 - 1)*sc[13] + ( 0.40000000000000002 - 1)*sc[2] + ( 0.46000000000000002 - 1)*sc[4] + ( 0.75 - 1)*sc[0];
    /* forward */
    phi_f = sc[1]*sc[11];
    k_f = 1.0000000000000002e-06 * 5116000000000
                * exp(0.44 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0.44 * invT + 0.50321666580471969 *  0  * invT2;
    /* pressure-fall-off */
    k_0 = 6.328e+19 * exp(-1.3999999999999999 * tc[0] - 0.50321666580471969 * (0) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -1.3999999999999999 * invT + 0.50321666580471969 * (0) * invT2;
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
        dqdci = (0.75 - 1)*dcdc_fac;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[8] += dqdci;                /* dwdot[HO2]/d[H2] */
        J[11] -= dqdci;               /* dwdot[O2]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[11];
        J[16] -= dqdci;               /* dwdot[H]/d[H] */
        J[23] += dqdci;               /* dwdot[HO2]/d[H] */
        J[26] -= dqdci;               /* dwdot[O2]/d[H] */
        /* d()/d[AR] */
        dqdci = (0.40000000000000002 - 1)*dcdc_fac;
        J[31] -= dqdci;               /* dwdot[H]/d[AR] */
        J[38] += dqdci;               /* dwdot[HO2]/d[AR] */
        J[41] -= dqdci;               /* dwdot[O2]/d[AR] */
        /* d()/d[HE] */
        dqdci = (0.46000000000000002 - 1)*dcdc_fac;
        J[61] -= dqdci;               /* dwdot[H]/d[HE] */
        J[68] += dqdci;               /* dwdot[HO2]/d[HE] */
        J[71] -= dqdci;               /* dwdot[O2]/d[HE] */
        /* d()/d[HO2] */
        dqdci =  - k_r;
        J[121] -= dqdci;              /* dwdot[H]/d[HO2] */
        J[128] += dqdci;              /* dwdot[HO2]/d[HO2] */
        J[131] -= dqdci;              /* dwdot[O2]/d[HO2] */
        /* d()/d[H2O] */
        dqdci = (11.890000000000001 - 1)*dcdc_fac;
        J[136] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[143] += dqdci;              /* dwdot[HO2]/d[H2O] */
        J[146] -= dqdci;              /* dwdot[O2]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (1.0900000000000001 - 1)*dcdc_fac;
        J[151] -= dqdci;              /* dwdot[H]/d[CO] */
        J[158] += dqdci;              /* dwdot[HO2]/d[CO] */
        J[161] -= dqdci;              /* dwdot[O2]/d[CO] */
        /* d()/d[O2] */
        dqdci = (0.84999999999999998 - 1)*dcdc_fac + k_f*sc[1];
        J[166] -= dqdci;              /* dwdot[H]/d[O2] */
        J[173] += dqdci;              /* dwdot[HO2]/d[O2] */
        J[176] -= dqdci;              /* dwdot[O2]/d[O2] */
        /* d()/d[CO2] */
        dqdci = (2.1800000000000002 - 1)*dcdc_fac;
        J[196] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[203] += dqdci;              /* dwdot[HO2]/d[CO2] */
        J[206] -= dqdci;              /* dwdot[O2]/d[CO2] */
    }
    else {
        dqdc[0] = 0.75*dcdc_fac;
        dqdc[1] = dcdc_fac + k_f*sc[11];
        dqdc[2] = 0.40000000000000002*dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = 0.46000000000000002*dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac - k_r;
        dqdc[9] = 11.890000000000001*dcdc_fac;
        dqdc[10] = 1.0900000000000001*dcdc_fac;
        dqdc[11] = 0.84999999999999998*dcdc_fac + k_f*sc[1];
        dqdc[12] = dcdc_fac;
        dqdc[13] = 2.1800000000000002*dcdc_fac;
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
    alpha = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[9] + ( 1.75 - 1)*sc[10] + ( 3.6000000000000001 - 1)*sc[13] + ( 0.69999999999999996 - 1)*sc[2] + ( 0.69999999999999996 - 1)*sc[4];
    /* forward */
    phi_f = pow(sc[6], 2.000000);
    k_f = 1.0000000000000002e-06 * 111000000000000
                * exp(-0.37 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -0.37 * invT + 0.50321666580471969 *  0  * invT2;
    /* pressure-fall-off */
    k_0 = 2.01e+17 * exp(-0.58399999999999996 * tc[0] - 0.50321666580471969 * (-2293) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -0.58399999999999996 * invT + 0.50321666580471969 * (-2293) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-0.73460000000000003)*exp(-T/94);
    Fcent2 = 0.73460000000000003 * exp(-T/1756);
    Fcent3 = exp(-5182 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/94
        -Fcent2/1756
        + Fcent3*5182*invT2);
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
        dqdci = (2 - 1)*dcdc_fac;
        J[6] += -2 * dqdci;           /* dwdot[OH]/d[H2] */
        J[12] += dqdci;               /* dwdot[H2O2]/d[H2] */
        /* d()/d[AR] */
        dqdci = (0.69999999999999996 - 1)*dcdc_fac;
        J[36] += -2 * dqdci;          /* dwdot[OH]/d[AR] */
        J[42] += dqdci;               /* dwdot[H2O2]/d[AR] */
        /* d()/d[HE] */
        dqdci = (0.69999999999999996 - 1)*dcdc_fac;
        J[66] += -2 * dqdci;          /* dwdot[OH]/d[HE] */
        J[72] += dqdci;               /* dwdot[H2O2]/d[HE] */
        /* d()/d[OH] */
        dqdci =  + k_f*2.000000*sc[6];
        J[96] += -2 * dqdci;          /* dwdot[OH]/d[OH] */
        J[102] += dqdci;              /* dwdot[H2O2]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (6 - 1)*dcdc_fac;
        J[141] += -2 * dqdci;         /* dwdot[OH]/d[H2O] */
        J[147] += dqdci;              /* dwdot[H2O2]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (1.75 - 1)*dcdc_fac;
        J[156] += -2 * dqdci;         /* dwdot[OH]/d[CO] */
        J[162] += dqdci;              /* dwdot[H2O2]/d[CO] */
        /* d()/d[H2O2] */
        dqdci =  - k_r;
        J[186] += -2 * dqdci;         /* dwdot[OH]/d[H2O2] */
        J[192] += dqdci;              /* dwdot[H2O2]/d[H2O2] */
        /* d()/d[CO2] */
        dqdci = (3.6000000000000001 - 1)*dcdc_fac;
        J[201] += -2 * dqdci;         /* dwdot[OH]/d[CO2] */
        J[207] += dqdci;              /* dwdot[H2O2]/d[CO2] */
    }
    else {
        dqdc[0] = 2*dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = 0.69999999999999996*dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = 0.69999999999999996*dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = dcdc_fac + k_f*2.000000*sc[6];
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = 6*dcdc_fac;
        dqdc[10] = 1.75*dcdc_fac;
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac - k_r;
        dqdc[13] = 3.6000000000000001*dcdc_fac;
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
    alpha = mixture + ( 2 - 1)*sc[0] + ( 12 - 1)*sc[9] + ( 1.75 - 1)*sc[10] + ( 3.6000000000000001 - 1)*sc[13] + ( 0.69999999999999996 - 1)*sc[2] + ( 0.69999999999999996 - 1)*sc[4];
    /* forward */
    phi_f = sc[5]*sc[10];
    k_f = 1.0000000000000002e-06 * 13620000000
                * exp(0 * tc[0] - 0.50321666580471969 * (2384) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  2384  * invT2;
    /* pressure-fall-off */
    k_0 = 1.1729999999999999e+24 * exp(-2.79 * tc[0] - 0.50321666580471969 * (4191) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -2.79 * invT + 0.50321666580471969 * (4191) * invT2;
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
        dqdci = (2 - 1)*dcdc_fac;
        J[5] -= dqdci;                /* dwdot[O]/d[H2] */
        J[10] -= dqdci;               /* dwdot[CO]/d[H2] */
        J[13] += dqdci;               /* dwdot[CO2]/d[H2] */
        /* d()/d[AR] */
        dqdci = (0.69999999999999996 - 1)*dcdc_fac;
        J[35] -= dqdci;               /* dwdot[O]/d[AR] */
        J[40] -= dqdci;               /* dwdot[CO]/d[AR] */
        J[43] += dqdci;               /* dwdot[CO2]/d[AR] */
        /* d()/d[HE] */
        dqdci = (0.69999999999999996 - 1)*dcdc_fac;
        J[65] -= dqdci;               /* dwdot[O]/d[HE] */
        J[70] -= dqdci;               /* dwdot[CO]/d[HE] */
        J[73] += dqdci;               /* dwdot[CO2]/d[HE] */
        /* d()/d[O] */
        dqdci =  + k_f*sc[10];
        J[80] -= dqdci;               /* dwdot[O]/d[O] */
        J[85] -= dqdci;               /* dwdot[CO]/d[O] */
        J[88] += dqdci;               /* dwdot[CO2]/d[O] */
        /* d()/d[H2O] */
        dqdci = (12 - 1)*dcdc_fac;
        J[140] -= dqdci;              /* dwdot[O]/d[H2O] */
        J[145] -= dqdci;              /* dwdot[CO]/d[H2O] */
        J[148] += dqdci;              /* dwdot[CO2]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (1.75 - 1)*dcdc_fac + k_f*sc[5];
        J[155] -= dqdci;              /* dwdot[O]/d[CO] */
        J[160] -= dqdci;              /* dwdot[CO]/d[CO] */
        J[163] += dqdci;              /* dwdot[CO2]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (3.6000000000000001 - 1)*dcdc_fac - k_r;
        J[200] -= dqdci;              /* dwdot[O]/d[CO2] */
        J[205] -= dqdci;              /* dwdot[CO]/d[CO2] */
        J[208] += dqdci;              /* dwdot[CO2]/d[CO2] */
    }
    else {
        dqdc[0] = 2*dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = 0.69999999999999996*dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = 0.69999999999999996*dcdc_fac;
        dqdc[5] = dcdc_fac + k_f*sc[10];
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        dqdc[9] = 12*dcdc_fac;
        dqdc[10] = 1.75*dcdc_fac + k_f*sc[5];
        dqdc[11] = dcdc_fac;
        dqdc[12] = dcdc_fac;
        dqdc[13] = 3.6000000000000001*dcdc_fac - k_r;
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
    alpha = mixture + ( 0 - 1)*sc[0] + ( 0 - 1)*sc[9] + ( 0 - 1)*sc[13] + ( 0.63 - 1)*sc[2] + ( 0.63 - 1)*sc[4];
    /* forward */
    phi_f = pow(sc[1], 2.000000);
    k_f = 1.0000000000000002e-12 * 1.78e+18
                * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1 * invT + 0.50321666580471969 *  0  * invT2;
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
        dqdci = (0 - 1)*q_nocor - k_r;
        J[0] += dqdci;                /* dwdot[H2]/d[H2] */
        J[1] += -2 * dqdci;           /* dwdot[H]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*2.000000*sc[1];
        J[15] += dqdci;               /* dwdot[H2]/d[H] */
        J[16] += -2 * dqdci;          /* dwdot[H]/d[H] */
        /* d()/d[AR] */
        dqdci = (0.63 - 1)*q_nocor;
        J[30] += dqdci;               /* dwdot[H2]/d[AR] */
        J[31] += -2 * dqdci;          /* dwdot[H]/d[AR] */
        /* d()/d[HE] */
        dqdci = (0.63 - 1)*q_nocor;
        J[60] += dqdci;               /* dwdot[H2]/d[HE] */
        J[61] += -2 * dqdci;          /* dwdot[H]/d[HE] */
        /* d()/d[H2O] */
        dqdci = (0 - 1)*q_nocor;
        J[135] += dqdci;              /* dwdot[H2]/d[H2O] */
        J[136] += -2 * dqdci;         /* dwdot[H]/d[H2O] */
        /* d()/d[CO2] */
        dqdci = (0 - 1)*q_nocor;
        J[195] += dqdci;              /* dwdot[H2]/d[CO2] */
        J[196] += -2 * dqdci;         /* dwdot[H]/d[CO2] */
    }
    else {
        dqdc[0] =  - k_r;
        dqdc[1] = q_nocor + k_f*2.000000*sc[1];
        dqdc[2] = 0.63*q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = 0.63*q_nocor;
        dqdc[5] = q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[10] = q_nocor;
        dqdc[11] = q_nocor;
        dqdc[12] = q_nocor;
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
    alpha = mixture + ( 2 - 1)*sc[0] + ( 6.2999999999999998 - 1)*sc[9] + ( 1.75 - 1)*sc[10] + ( 3.6000000000000001 - 1)*sc[13] + ( 0.38 - 1)*sc[2] + ( 0.38 - 1)*sc[4];
    /* forward */
    phi_f = sc[1]*sc[6];
    k_f = 1.0000000000000002e-12 * 4.4e+22
                * exp(-2 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -2 * invT + 0.50321666580471969 *  0  * invT2;
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
        dqdci = (2 - 1)*q_nocor;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[6] -= dqdci;                /* dwdot[OH]/d[H2] */
        J[9] += dqdci;                /* dwdot[H2O]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[6];
        J[16] -= dqdci;               /* dwdot[H]/d[H] */
        J[21] -= dqdci;               /* dwdot[OH]/d[H] */
        J[24] += dqdci;               /* dwdot[H2O]/d[H] */
        /* d()/d[AR] */
        dqdci = (0.38 - 1)*q_nocor;
        J[31] -= dqdci;               /* dwdot[H]/d[AR] */
        J[36] -= dqdci;               /* dwdot[OH]/d[AR] */
        J[39] += dqdci;               /* dwdot[H2O]/d[AR] */
        /* d()/d[HE] */
        dqdci = (0.38 - 1)*q_nocor;
        J[61] -= dqdci;               /* dwdot[H]/d[HE] */
        J[66] -= dqdci;               /* dwdot[OH]/d[HE] */
        J[69] += dqdci;               /* dwdot[H2O]/d[HE] */
        /* d()/d[OH] */
        dqdci =  + k_f*sc[1];
        J[91] -= dqdci;               /* dwdot[H]/d[OH] */
        J[96] -= dqdci;               /* dwdot[OH]/d[OH] */
        J[99] += dqdci;               /* dwdot[H2O]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (6.2999999999999998 - 1)*q_nocor - k_r;
        J[136] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[141] -= dqdci;              /* dwdot[OH]/d[H2O] */
        J[144] += dqdci;              /* dwdot[H2O]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (1.75 - 1)*q_nocor;
        J[151] -= dqdci;              /* dwdot[H]/d[CO] */
        J[156] -= dqdci;              /* dwdot[OH]/d[CO] */
        J[159] += dqdci;              /* dwdot[H2O]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (3.6000000000000001 - 1)*q_nocor;
        J[196] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[201] -= dqdci;              /* dwdot[OH]/d[CO2] */
        J[204] += dqdci;              /* dwdot[H2O]/d[CO2] */
    }
    else {
        dqdc[0] = 2*q_nocor;
        dqdc[1] = q_nocor + k_f*sc[6];
        dqdc[2] = 0.38*q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = 0.38*q_nocor;
        dqdc[5] = q_nocor;
        dqdc[6] = q_nocor + k_f*sc[1];
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = 6.2999999999999998*q_nocor - k_r;
        dqdc[10] = 1.75*q_nocor;
        dqdc[11] = q_nocor;
        dqdc[12] = q_nocor;
        dqdc[13] = 3.6000000000000001*q_nocor;
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
    alpha = mixture + ( 2 - 1)*sc[0] + ( 12 - 1)*sc[9] + ( 1.75 - 1)*sc[10] + ( 3.6000000000000001 - 1)*sc[13] + ( 0.69999999999999996 - 1)*sc[2] + ( 0.69999999999999996 - 1)*sc[4];
    /* forward */
    phi_f = sc[1]*sc[5];
    k_f = 1.0000000000000002e-12 * 9.428e+18
                * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1 * invT + 0.50321666580471969 *  0  * invT2;
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
        dqdci = (2 - 1)*q_nocor;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[5] -= dqdci;                /* dwdot[O]/d[H2] */
        J[6] += dqdci;                /* dwdot[OH]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[5];
        J[16] -= dqdci;               /* dwdot[H]/d[H] */
        J[20] -= dqdci;               /* dwdot[O]/d[H] */
        J[21] += dqdci;               /* dwdot[OH]/d[H] */
        /* d()/d[AR] */
        dqdci = (0.69999999999999996 - 1)*q_nocor;
        J[31] -= dqdci;               /* dwdot[H]/d[AR] */
        J[35] -= dqdci;               /* dwdot[O]/d[AR] */
        J[36] += dqdci;               /* dwdot[OH]/d[AR] */
        /* d()/d[HE] */
        dqdci = (0.69999999999999996 - 1)*q_nocor;
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
        dqdci = (12 - 1)*q_nocor;
        J[136] -= dqdci;              /* dwdot[H]/d[H2O] */
        J[140] -= dqdci;              /* dwdot[O]/d[H2O] */
        J[141] += dqdci;              /* dwdot[OH]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (1.75 - 1)*q_nocor;
        J[151] -= dqdci;              /* dwdot[H]/d[CO] */
        J[155] -= dqdci;              /* dwdot[O]/d[CO] */
        J[156] += dqdci;              /* dwdot[OH]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (3.6000000000000001 - 1)*q_nocor;
        J[196] -= dqdci;              /* dwdot[H]/d[CO2] */
        J[200] -= dqdci;              /* dwdot[O]/d[CO2] */
        J[201] += dqdci;              /* dwdot[OH]/d[CO2] */
    }
    else {
        dqdc[0] = 2*q_nocor;
        dqdc[1] = q_nocor + k_f*sc[5];
        dqdc[2] = 0.69999999999999996*q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = 0.69999999999999996*q_nocor;
        dqdc[5] = q_nocor + k_f*sc[1];
        dqdc[6] = q_nocor - k_r;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = 12*q_nocor;
        dqdc[10] = 1.75*q_nocor;
        dqdc[11] = q_nocor;
        dqdc[12] = q_nocor;
        dqdc[13] = 3.6000000000000001*q_nocor;
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
    alpha = mixture + ( 2.3999999999999999 - 1)*sc[0] + ( 15.4 - 1)*sc[9] + ( 1.75 - 1)*sc[10] + ( 3.6000000000000001 - 1)*sc[13] + ( 0.82999999999999996 - 1)*sc[2] + ( 0.82999999999999996 - 1)*sc[4];
    /* forward */
    phi_f = pow(sc[5], 2.000000);
    k_f = 1.0000000000000002e-12 * 1.2e+17
                * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1 * invT + 0.50321666580471969 *  0  * invT2;
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
        dqdci = (2.3999999999999999 - 1)*q_nocor;
        J[5] += -2 * dqdci;           /* dwdot[O]/d[H2] */
        J[11] += dqdci;               /* dwdot[O2]/d[H2] */
        /* d()/d[AR] */
        dqdci = (0.82999999999999996 - 1)*q_nocor;
        J[35] += -2 * dqdci;          /* dwdot[O]/d[AR] */
        J[41] += dqdci;               /* dwdot[O2]/d[AR] */
        /* d()/d[HE] */
        dqdci = (0.82999999999999996 - 1)*q_nocor;
        J[65] += -2 * dqdci;          /* dwdot[O]/d[HE] */
        J[71] += dqdci;               /* dwdot[O2]/d[HE] */
        /* d()/d[O] */
        dqdci =  + k_f*2.000000*sc[5];
        J[80] += -2 * dqdci;          /* dwdot[O]/d[O] */
        J[86] += dqdci;               /* dwdot[O2]/d[O] */
        /* d()/d[H2O] */
        dqdci = (15.4 - 1)*q_nocor;
        J[140] += -2 * dqdci;         /* dwdot[O]/d[H2O] */
        J[146] += dqdci;              /* dwdot[O2]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (1.75 - 1)*q_nocor;
        J[155] += -2 * dqdci;         /* dwdot[O]/d[CO] */
        J[161] += dqdci;              /* dwdot[O2]/d[CO] */
        /* d()/d[O2] */
        dqdci =  - k_r;
        J[170] += -2 * dqdci;         /* dwdot[O]/d[O2] */
        J[176] += dqdci;              /* dwdot[O2]/d[O2] */
        /* d()/d[CO2] */
        dqdci = (3.6000000000000001 - 1)*q_nocor;
        J[200] += -2 * dqdci;         /* dwdot[O]/d[CO2] */
        J[206] += dqdci;              /* dwdot[O2]/d[CO2] */
    }
    else {
        dqdc[0] = 2.3999999999999999*q_nocor;
        dqdc[1] = q_nocor;
        dqdc[2] = 0.82999999999999996*q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = 0.82999999999999996*q_nocor;
        dqdc[5] = q_nocor + k_f*2.000000*sc[5];
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        dqdc[9] = 15.4*q_nocor;
        dqdc[10] = 1.75*q_nocor;
        dqdc[11] = q_nocor - k_r;
        dqdc[12] = q_nocor;
        dqdc[13] = 3.6000000000000001*q_nocor;
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
    alpha = mixture + ( 2 - 1)*sc[0] + ( 0 - 1)*sc[9] + ( 1.75 - 1)*sc[10] + ( 3.6000000000000001 - 1)*sc[13];
    /* forward */
    phi_f = sc[7];
    k_f = 1.0000000000000002e-06 * 1.87e+17
                * exp(-1 * tc[0] - 0.50321666580471969 * (17000) * invT);
    dlnkfdT = -1 * invT + 0.50321666580471969 *  17000  * invT2;
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
        dqdci = (2 - 1)*q_nocor;
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
        dqdci = (0 - 1)*q_nocor;
        J[136] += dqdci;              /* dwdot[H]/d[H2O] */
        J[142] -= dqdci;              /* dwdot[HCO]/d[H2O] */
        J[145] += dqdci;              /* dwdot[CO]/d[H2O] */
        /* d()/d[CO] */
        dqdci = (1.75 - 1)*q_nocor - k_r*sc[1];
        J[151] += dqdci;              /* dwdot[H]/d[CO] */
        J[157] -= dqdci;              /* dwdot[HCO]/d[CO] */
        J[160] += dqdci;              /* dwdot[CO]/d[CO] */
        /* d()/d[CO2] */
        dqdci = (3.6000000000000001 - 1)*q_nocor;
        J[196] += dqdci;              /* dwdot[H]/d[CO2] */
        J[202] -= dqdci;              /* dwdot[HCO]/d[CO2] */
        J[205] += dqdci;              /* dwdot[CO]/d[CO2] */
    }
    else {
        dqdc[0] = 2*q_nocor;
        dqdc[1] = q_nocor - k_r*sc[10];
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor;
        dqdc[5] = q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor + k_f;
        dqdc[8] = q_nocor;
        dqdc[10] = 1.75*q_nocor - k_r*sc[1];
        dqdc[11] = q_nocor;
        dqdc[12] = q_nocor;
        dqdc[13] = 3.6000000000000001*q_nocor;
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
    k_f = 1.0000000000000002e-06 * 26440000000000000
                * exp(-0.67069999999999996 * tc[0] - 0.50321666580471969 * (17041) * invT);
    dlnkfdT = -0.67069999999999996 * invT + 0.50321666580471969 *  17041  * invT2;
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
    k_f = 1.0000000000000002e-06 * 45890
                * exp(2.7000000000000002 * tc[0] - 0.50321666580471969 * (6260) * invT);
    dlnkfdT = 2.7000000000000002 * invT + 0.50321666580471969 *  6260  * invT2;
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
    k_f = 1.0000000000000002e-06 * 173400000
                * exp(1.51 * tc[0] - 0.50321666580471969 * (3430) * invT);
    dlnkfdT = 1.51 * invT + 0.50321666580471969 *  3430  * invT2;
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
    k_f = 1.0000000000000002e-06 * 39730
                * exp(2.3999999999999999 * tc[0] - 0.50321666580471969 * (-2110) * invT);
    dlnkfdT = 2.3999999999999999 * invT + 0.50321666580471969 *  -2110  * invT2;
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
    k_f = 1.0000000000000002e-12 * 90000000000000000
                * exp(-0.59999999999999998 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -0.59999999999999998 * invT + 0.50321666580471969 *  0  * invT2;
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
    k_f = 1.0000000000000002e-12 * 5.624e+19
                * exp(-1.25 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1.25 * invT + 0.50321666580471969 *  0  * invT2;
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
    k_f = 1.0000000000000002e-12 * 5.5e+20
                * exp(-2 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -2 * invT + 0.50321666580471969 *  0  * invT2;
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
    k_f = 1.0000000000000002e-06 * 591600
                * exp(2.4329999999999998 * tc[0] - 0.50321666580471969 * (53502) * invT);
    dlnkfdT = 2.4329999999999998 * invT + 0.50321666580471969 *  53502  * invT2;
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
    k_f = 1.0000000000000002e-06 * 3970000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (671) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  671  * invT2;
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
    k_f = 1.0000000000000002e-06 * 74850000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (295) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  295  * invT2;
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
    k_f = 1.0000000000000002e-06 * 40000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
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
    k_f = 1.0000000000000002e-06 * 23750000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-500) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  -500  * invT2;
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
    k_f = 1.0000000000000002e-06 * 10000000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (17330) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  17330  * invT2;
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
    k_f = 1.0000000000000002e-06 * 130000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-1630) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  -1630  * invT2;
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
    k_f = 1.0000000000000002e-06 * 365800000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (12000) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  12000  * invT2;
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
    k_f = 1.0000000000000002e-06 * 6050000
                * exp(2 * tc[0] - 0.50321666580471969 * (5200) * invT);
    dlnkfdT = 2 * invT + 0.50321666580471969 *  5200  * invT2;
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
    k_f = 1.0000000000000002e-06 * 24100000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (3970) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  3970  * invT2;
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
    k_f = 1.0000000000000002e-06 * 9630000
                * exp(2 * tc[0] - 0.50321666580471969 * (3970) * invT);
    dlnkfdT = 2 * invT + 0.50321666580471969 *  3970  * invT2;
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
    k_f = 1.0000000000000002e-06 * 2000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (427) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  427  * invT2;
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
    k_f = 1.0000000000000002e-06 * 2.6700000000000001e+41
                * exp(-7 * tc[0] - 0.50321666580471969 * (37600) * invT);
    dlnkfdT = -7 * invT + 0.50321666580471969 *  37600  * invT2;
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
    k_f = 1.0000000000000002e-06 * 800000000000
                * exp(0.14000000000000001 * tc[0] - 0.50321666580471969 * (7352) * invT);
    dlnkfdT = 0.14000000000000001 * invT + 0.50321666580471969 *  7352  * invT2;
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
    k_f = 1.0000000000000002e-06 * 87840000000
                * exp(0.029999999999999999 * tc[0] - 0.50321666580471969 * (-16) * invT);
    dlnkfdT = 0.029999999999999999 * invT + 0.50321666580471969 *  -16  * invT2;
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
    k_f = 1.0000000000000002e-06 * 1119000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (47700) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  47700  * invT2;
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
    k_f = 1.0000000000000002e-06 * 30100000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (23000) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  23000  * invT2;
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
    k_f = 1.0000000000000002e-06 * 120000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
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
    k_f = 1.0000000000000002e-06 * 30000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
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
    k_f = 1.0000000000000002e-06 * 30000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
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
    k_f = 1.0000000000000002e-06 * 30200000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
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
    k_f = 1.0000000000000002e-06 * 2.244e+18
                * exp(-1 * tc[0] - 0.50321666580471969 * (17000) * invT);
    dlnkfdT = -1 * invT + 0.50321666580471969 *  17000  * invT2;
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
    k_f = 1.0000000000000002e-06 * 12040000000
                * exp(0.80700000000000005 * tc[0] - 0.50321666580471969 * (-727) * invT);
    dlnkfdT = 0.80700000000000005 * invT + 0.50321666580471969 *  -727  * invT2;
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

    double c_R[14], dcRdT[14], e_RT[14];
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
    for (int k = 0; k < 14; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[210+k];
    }

    double cmixinv = 1.0/cmix;
    double tmp1 = ehmix*cmixinv;
    double tmp3 = cmixinv*T;
    double tmp2 = tmp1*tmp3;
    double dehmixdc;
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

return;
}
#endif


#ifndef AMREX_USE_CUDA
/*compute the reaction Jacobian on CPU */
void aJacobian(double *  J, double *  sc, double T, int consP)
{
    for (int i=0; i<225; i++) {
        J[i] = 0.0;
    }

    double wdot[14];
    for (int k=0; k<14; k++) {
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
    for (int k = 0; k < 14; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    double g_RT[14];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    double h_RT[14];
    speciesEnthalpy(h_RT, tc);

    double phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    double dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    double dqdci, dcdc_fac, dqdc[14];
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

    double c_R[14], dcRdT[14], e_RT[14];
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
    for (int k = 0; k < 14; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[210+k];
    }

    double cmixinv = 1.0/cmix;
    double tmp1 = ehmix*cmixinv;
    double tmp3 = cmixinv*T;
    double tmp2 = tmp1*tmp3;
    double dehmixdc;
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


/*compute an approx to the reaction Jacobian */
AMREX_GPU_HOST_DEVICE void aJacobian_precond(double *  J, double *  sc, double T, int HP)
{
    for (int i=0; i<225; i++) {
        J[i] = 0.0;
    }

    double wdot[14];
    for (int k=0; k<14; k++) {
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
    for (int k = 0; k < 14; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    double g_RT[14];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    double h_RT[14];
    speciesEnthalpy(h_RT, tc);

    double phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    double dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    double dqdci, dcdc_fac, dqdc[14];
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
    alpha = mixture + ( 0.84999999999999998 - 1)*sc[11] + ( 11.890000000000001 - 1)*sc[9] + ( 1.0900000000000001 - 1)*sc[10] + ( 2.1800000000000002 - 1)*sc[13] + ( 0.40000000000000002 - 1)*sc[2] + ( 0.46000000000000002 - 1)*sc[4] + ( 0.75 - 1)*sc[0];
    /* forward */
    phi_f = sc[1]*sc[11];
    k_f = 1.0000000000000002e-06 * 5116000000000
                * exp(0.44 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0.44 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* pressure-fall-off */
    k_0 = 6.328e+19 * exp(-1.3999999999999999 * tc[0] - 0.50321666580471969 * (0) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -1.3999999999999999 * invT + 0.50321666580471969 * (0) * invT2;
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
    dcdc_fac = 0.0;
    dqdc[0] = 0.75*dcdc_fac;
    dqdc[1] = dcdc_fac + k_f*sc[11];
    dqdc[2] = 0.40000000000000002*dcdc_fac;
    dqdc[3] = dcdc_fac;
    dqdc[4] = 0.46000000000000002*dcdc_fac;
    dqdc[5] = dcdc_fac;
    dqdc[6] = dcdc_fac;
    dqdc[7] = dcdc_fac;
    dqdc[8] = dcdc_fac - k_r;
    dqdc[9] = 11.890000000000001*dcdc_fac;
    dqdc[10] = 1.0900000000000001*dcdc_fac;
    dqdc[11] = 0.84999999999999998*dcdc_fac + k_f*sc[1];
    dqdc[12] = dcdc_fac;
    dqdc[13] = 2.1800000000000002*dcdc_fac;
    for (int k=0; k<14; k++) {
        J[15*k+1] -= dqdc[k];
        J[15*k+8] += dqdc[k];
        J[15*k+11] -= dqdc[k];
    }
    J[211] -= dqdT; /* dwdot[H]/dT */
    J[218] += dqdT; /* dwdot[HO2]/dT */
    J[221] -= dqdT; /* dwdot[O2]/dT */

    /*reaction 2: OH + OH (+M) <=> H2O2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[9] + ( 1.75 - 1)*sc[10] + ( 3.6000000000000001 - 1)*sc[13] + ( 0.69999999999999996 - 1)*sc[2] + ( 0.69999999999999996 - 1)*sc[4];
    /* forward */
    phi_f = pow(sc[6], 2.000000);
    k_f = 1.0000000000000002e-06 * 111000000000000
                * exp(-0.37 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -0.37 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* pressure-fall-off */
    k_0 = 2.01e+17 * exp(-0.58399999999999996 * tc[0] - 0.50321666580471969 * (-2293) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -0.58399999999999996 * invT + 0.50321666580471969 * (-2293) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-0.73460000000000003)*exp(-T/94);
    Fcent2 = 0.73460000000000003 * exp(-T/1756);
    Fcent3 = exp(-5182 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/94
        -Fcent2/1756
        + Fcent3*5182*invT2);
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
    dcdc_fac = 0.0;
    dqdc[0] = 2*dcdc_fac;
    dqdc[1] = dcdc_fac;
    dqdc[2] = 0.69999999999999996*dcdc_fac;
    dqdc[3] = dcdc_fac;
    dqdc[4] = 0.69999999999999996*dcdc_fac;
    dqdc[5] = dcdc_fac;
    dqdc[6] = dcdc_fac + k_f*2.000000*sc[6];
    dqdc[7] = dcdc_fac;
    dqdc[8] = dcdc_fac;
    dqdc[9] = 6*dcdc_fac;
    dqdc[10] = 1.75*dcdc_fac;
    dqdc[11] = dcdc_fac;
    dqdc[12] = dcdc_fac - k_r;
    dqdc[13] = 3.6000000000000001*dcdc_fac;
    for (int k=0; k<14; k++) {
        J[15*k+6] += -2 * dqdc[k];
        J[15*k+12] += dqdc[k];
    }
    J[216] += -2 * dqdT; /* dwdot[OH]/dT */
    J[222] += dqdT; /* dwdot[H2O2]/dT */

    /*reaction 3: CO + O (+M) <=> CO2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 12 - 1)*sc[9] + ( 1.75 - 1)*sc[10] + ( 3.6000000000000001 - 1)*sc[13] + ( 0.69999999999999996 - 1)*sc[2] + ( 0.69999999999999996 - 1)*sc[4];
    /* forward */
    phi_f = sc[5]*sc[10];
    k_f = 1.0000000000000002e-06 * 13620000000
                * exp(0 * tc[0] - 0.50321666580471969 * (2384) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (2384)  * invT2;
    /* pressure-fall-off */
    k_0 = 1.1729999999999999e+24 * exp(-2.79 * tc[0] - 0.50321666580471969 * (4191) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -2.79 * invT + 0.50321666580471969 * (4191) * invT2;
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
    dcdc_fac = 0.0;
    dqdc[0] = 2*dcdc_fac;
    dqdc[1] = dcdc_fac;
    dqdc[2] = 0.69999999999999996*dcdc_fac;
    dqdc[3] = dcdc_fac;
    dqdc[4] = 0.69999999999999996*dcdc_fac;
    dqdc[5] = dcdc_fac + k_f*sc[10];
    dqdc[6] = dcdc_fac;
    dqdc[7] = dcdc_fac;
    dqdc[8] = dcdc_fac;
    dqdc[9] = 12*dcdc_fac;
    dqdc[10] = 1.75*dcdc_fac + k_f*sc[5];
    dqdc[11] = dcdc_fac;
    dqdc[12] = dcdc_fac;
    dqdc[13] = 3.6000000000000001*dcdc_fac - k_r;
    for (int k=0; k<14; k++) {
        J[15*k+5] -= dqdc[k];
        J[15*k+10] -= dqdc[k];
        J[15*k+13] += dqdc[k];
    }
    J[215] -= dqdT; /* dwdot[O]/dT */
    J[220] -= dqdT; /* dwdot[CO]/dT */
    J[223] += dqdT; /* dwdot[CO2]/dT */

    /*reaction 4: H + H + M <=> H2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 0 - 1)*sc[0] + ( 0 - 1)*sc[9] + ( 0 - 1)*sc[13] + ( 0.63 - 1)*sc[2] + ( 0.63 - 1)*sc[4];
    /* forward */
    phi_f = pow(sc[1], 2.000000);
    k_f = 1.0000000000000002e-12 * 1.78e+18
                * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1 * invT + 0.50321666580471969 *  (0)  * invT2;
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
    dqdc[0] =  - k_r;
    dqdc[1] = q_nocor + k_f*2.000000*sc[1];
    dqdc[2] = 0.63*q_nocor;
    dqdc[3] = q_nocor;
    dqdc[4] = 0.63*q_nocor;
    dqdc[5] = q_nocor;
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    dqdc[9] = 0.0;
    dqdc[10] = q_nocor;
    dqdc[11] = q_nocor;
    dqdc[12] = q_nocor;
    dqdc[13] = 0.0;
    for (int k=0; k<14; k++) {
        J[15*k+0] += dqdc[k];
        J[15*k+1] += -2 * dqdc[k];
    }
    J[210] += dqdT; /* dwdot[H2]/dT */
    J[211] += -2 * dqdT; /* dwdot[H]/dT */

    /*reaction 5: H + OH + M <=> H2O + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 6.2999999999999998 - 1)*sc[9] + ( 1.75 - 1)*sc[10] + ( 3.6000000000000001 - 1)*sc[13] + ( 0.38 - 1)*sc[2] + ( 0.38 - 1)*sc[4];
    /* forward */
    phi_f = sc[1]*sc[6];
    k_f = 1.0000000000000002e-12 * 4.4e+22
                * exp(-2 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -2 * invT + 0.50321666580471969 *  (0)  * invT2;
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
    dqdc[0] = 2*q_nocor;
    dqdc[1] = q_nocor + k_f*sc[6];
    dqdc[2] = 0.38*q_nocor;
    dqdc[3] = q_nocor;
    dqdc[4] = 0.38*q_nocor;
    dqdc[5] = q_nocor;
    dqdc[6] = q_nocor + k_f*sc[1];
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    dqdc[9] = 6.2999999999999998*q_nocor - k_r;
    dqdc[10] = 1.75*q_nocor;
    dqdc[11] = q_nocor;
    dqdc[12] = q_nocor;
    dqdc[13] = 3.6000000000000001*q_nocor;
    for (int k=0; k<14; k++) {
        J[15*k+1] -= dqdc[k];
        J[15*k+6] -= dqdc[k];
        J[15*k+9] += dqdc[k];
    }
    J[211] -= dqdT; /* dwdot[H]/dT */
    J[216] -= dqdT; /* dwdot[OH]/dT */
    J[219] += dqdT; /* dwdot[H2O]/dT */

    /*reaction 6: O + H + M <=> OH + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 12 - 1)*sc[9] + ( 1.75 - 1)*sc[10] + ( 3.6000000000000001 - 1)*sc[13] + ( 0.69999999999999996 - 1)*sc[2] + ( 0.69999999999999996 - 1)*sc[4];
    /* forward */
    phi_f = sc[1]*sc[5];
    k_f = 1.0000000000000002e-12 * 9.428e+18
                * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1 * invT + 0.50321666580471969 *  (0)  * invT2;
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
    dqdc[0] = 2*q_nocor;
    dqdc[1] = q_nocor + k_f*sc[5];
    dqdc[2] = 0.69999999999999996*q_nocor;
    dqdc[3] = q_nocor;
    dqdc[4] = 0.69999999999999996*q_nocor;
    dqdc[5] = q_nocor + k_f*sc[1];
    dqdc[6] = q_nocor - k_r;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    dqdc[9] = 12*q_nocor;
    dqdc[10] = 1.75*q_nocor;
    dqdc[11] = q_nocor;
    dqdc[12] = q_nocor;
    dqdc[13] = 3.6000000000000001*q_nocor;
    for (int k=0; k<14; k++) {
        J[15*k+1] -= dqdc[k];
        J[15*k+5] -= dqdc[k];
        J[15*k+6] += dqdc[k];
    }
    J[211] -= dqdT; /* dwdot[H]/dT */
    J[215] -= dqdT; /* dwdot[O]/dT */
    J[216] += dqdT; /* dwdot[OH]/dT */

    /*reaction 7: O + O + M <=> O2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 2.3999999999999999 - 1)*sc[0] + ( 15.4 - 1)*sc[9] + ( 1.75 - 1)*sc[10] + ( 3.6000000000000001 - 1)*sc[13] + ( 0.82999999999999996 - 1)*sc[2] + ( 0.82999999999999996 - 1)*sc[4];
    /* forward */
    phi_f = pow(sc[5], 2.000000);
    k_f = 1.0000000000000002e-12 * 1.2e+17
                * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1 * invT + 0.50321666580471969 *  (0)  * invT2;
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
    dqdc[0] = 2.3999999999999999*q_nocor;
    dqdc[1] = q_nocor;
    dqdc[2] = 0.82999999999999996*q_nocor;
    dqdc[3] = q_nocor;
    dqdc[4] = 0.82999999999999996*q_nocor;
    dqdc[5] = q_nocor + k_f*2.000000*sc[5];
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    dqdc[9] = 15.4*q_nocor;
    dqdc[10] = 1.75*q_nocor;
    dqdc[11] = q_nocor - k_r;
    dqdc[12] = q_nocor;
    dqdc[13] = 3.6000000000000001*q_nocor;
    for (int k=0; k<14; k++) {
        J[15*k+5] += -2 * dqdc[k];
        J[15*k+11] += dqdc[k];
    }
    J[215] += -2 * dqdT; /* dwdot[O]/dT */
    J[221] += dqdT; /* dwdot[O2]/dT */

    /*reaction 8: HCO + M <=> CO + H + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 0 - 1)*sc[9] + ( 1.75 - 1)*sc[10] + ( 3.6000000000000001 - 1)*sc[13];
    /* forward */
    phi_f = sc[7];
    k_f = 1.0000000000000002e-06 * 1.87e+17
                * exp(-1 * tc[0] - 0.50321666580471969 * (17000) * invT);
    dlnkfdT = -1 * invT + 0.50321666580471969 *  (17000)  * invT2;
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
    dqdc[0] = 2*q_nocor;
    dqdc[1] = q_nocor - k_r*sc[10];
    dqdc[2] = q_nocor;
    dqdc[3] = q_nocor;
    dqdc[4] = q_nocor;
    dqdc[5] = q_nocor;
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor + k_f;
    dqdc[8] = q_nocor;
    dqdc[9] = 0.0;
    dqdc[10] = 1.75*q_nocor - k_r*sc[1];
    dqdc[11] = q_nocor;
    dqdc[12] = q_nocor;
    dqdc[13] = 3.6000000000000001*q_nocor;
    for (int k=0; k<14; k++) {
        J[15*k+1] += dqdc[k];
        J[15*k+7] -= dqdc[k];
        J[15*k+10] += dqdc[k];
    }
    J[211] += dqdT; /* dwdot[H]/dT */
    J[217] -= dqdT; /* dwdot[HCO]/dT */
    J[220] += dqdT; /* dwdot[CO]/dT */

    /*reaction 9: H + O2 <=> O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[11];
    k_f = 1.0000000000000002e-06 * 26440000000000000
                * exp(-0.67069999999999996 * tc[0] - 0.50321666580471969 * (17041) * invT);
    dlnkfdT = -0.67069999999999996 * invT + 0.50321666580471969 *  (17041)  * invT2;
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
    k_f = 1.0000000000000002e-06 * 45890
                * exp(2.7000000000000002 * tc[0] - 0.50321666580471969 * (6260) * invT);
    dlnkfdT = 2.7000000000000002 * invT + 0.50321666580471969 *  (6260)  * invT2;
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
    k_f = 1.0000000000000002e-06 * 173400000
                * exp(1.51 * tc[0] - 0.50321666580471969 * (3430) * invT);
    dlnkfdT = 1.51 * invT + 0.50321666580471969 *  (3430)  * invT2;
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
    k_f = 1.0000000000000002e-06 * 39730
                * exp(2.3999999999999999 * tc[0] - 0.50321666580471969 * (-2110) * invT);
    dlnkfdT = 2.3999999999999999 * invT + 0.50321666580471969 *  (-2110)  * invT2;
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
    k_f = 1.0000000000000002e-12 * 90000000000000000
                * exp(-0.59999999999999998 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -0.59999999999999998 * invT + 0.50321666580471969 *  (0)  * invT2;
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
    k_f = 1.0000000000000002e-12 * 5.624e+19
                * exp(-1.25 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1.25 * invT + 0.50321666580471969 *  (0)  * invT2;
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
    k_f = 1.0000000000000002e-12 * 5.5e+20
                * exp(-2 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -2 * invT + 0.50321666580471969 *  (0)  * invT2;
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
    k_f = 1.0000000000000002e-06 * 591600
                * exp(2.4329999999999998 * tc[0] - 0.50321666580471969 * (53502) * invT);
    dlnkfdT = 2.4329999999999998 * invT + 0.50321666580471969 *  (53502)  * invT2;
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
    k_f = 1.0000000000000002e-06 * 3970000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (671) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (671)  * invT2;
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
    k_f = 1.0000000000000002e-06 * 74850000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (295) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (295)  * invT2;
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
    k_f = 1.0000000000000002e-06 * 40000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
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
    k_f = 1.0000000000000002e-06 * 23750000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-500) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (-500)  * invT2;
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
    k_f = 1.0000000000000002e-06 * 10000000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (17330) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (17330)  * invT2;
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
    k_f = 1.0000000000000002e-06 * 130000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-1630) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (-1630)  * invT2;
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
    k_f = 1.0000000000000002e-06 * 365800000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (12000) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (12000)  * invT2;
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
    k_f = 1.0000000000000002e-06 * 6050000
                * exp(2 * tc[0] - 0.50321666580471969 * (5200) * invT);
    dlnkfdT = 2 * invT + 0.50321666580471969 *  (5200)  * invT2;
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
    k_f = 1.0000000000000002e-06 * 24100000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (3970) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (3970)  * invT2;
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
    k_f = 1.0000000000000002e-06 * 9630000
                * exp(2 * tc[0] - 0.50321666580471969 * (3970) * invT);
    dlnkfdT = 2 * invT + 0.50321666580471969 *  (3970)  * invT2;
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
    k_f = 1.0000000000000002e-06 * 2000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (427) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (427)  * invT2;
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
    k_f = 1.0000000000000002e-06 * 2.6700000000000001e+41
                * exp(-7 * tc[0] - 0.50321666580471969 * (37600) * invT);
    dlnkfdT = -7 * invT + 0.50321666580471969 *  (37600)  * invT2;
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
    k_f = 1.0000000000000002e-06 * 800000000000
                * exp(0.14000000000000001 * tc[0] - 0.50321666580471969 * (7352) * invT);
    dlnkfdT = 0.14000000000000001 * invT + 0.50321666580471969 *  (7352)  * invT2;
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
    k_f = 1.0000000000000002e-06 * 87840000000
                * exp(0.029999999999999999 * tc[0] - 0.50321666580471969 * (-16) * invT);
    dlnkfdT = 0.029999999999999999 * invT + 0.50321666580471969 *  (-16)  * invT2;
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
    k_f = 1.0000000000000002e-06 * 1119000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (47700) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (47700)  * invT2;
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
    k_f = 1.0000000000000002e-06 * 30100000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (23000) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (23000)  * invT2;
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
    k_f = 1.0000000000000002e-06 * 120000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
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
    k_f = 1.0000000000000002e-06 * 30000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
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
    k_f = 1.0000000000000002e-06 * 30000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
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
    k_f = 1.0000000000000002e-06 * 30200000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
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
    k_f = 1.0000000000000002e-06 * 2.244e+18
                * exp(-1 * tc[0] - 0.50321666580471969 * (17000) * invT);
    dlnkfdT = -1 * invT + 0.50321666580471969 *  (17000)  * invT2;
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
    k_f = 1.0000000000000002e-06 * 12040000000
                * exp(0.80700000000000005 * tc[0] - 0.50321666580471969 * (-727) * invT);
    dlnkfdT = 0.80700000000000005 * invT + 0.50321666580471969 *  (-727)  * invT2;
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

    double c_R[14], dcRdT[14], e_RT[14];
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
    for (int k = 0; k < 14; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[210+k];
    }

    double cmixinv = 1.0/cmix;
    double tmp1 = ehmix*cmixinv;
    double tmp3 = cmixinv*T;
    double tmp2 = tmp1*tmp3;
    double dehmixdc;
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


/*compute d(Cp/R)/dT and d(Cv/R)/dT at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE void dcvpRdT(double * species, double *  tc)
{

    /*temperature */
    double T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            +7.98052075e-03
            -3.89563020e-05 * tc[1]
            +6.04716282e-08 * tc[2]
            -2.95044704e-11 * tc[3];
        /*species 1: H */
        species[1] =
            +7.05332819e-13
            -3.99183928e-15 * tc[1]
            +6.90244896e-18 * tc[2]
            -3.71092933e-21 * tc[3];
        /*species 2: AR */
        species[2] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
        /*species 3: N2 */
        species[3] =
            +1.40824040e-03
            -7.92644400e-06 * tc[1]
            +1.69245450e-08 * tc[2]
            -9.77941600e-12 * tc[3];
        /*species 4: HE */
        species[4] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
        /*species 5: O */
        species[5] =
            -3.27931884e-03
            +1.32861279e-05 * tc[1]
            -1.83841987e-08 * tc[2]
            +8.45063884e-12 * tc[3];
        /*species 6: OH */
        species[6] =
            -3.22544939e-03
            +1.30552938e-05 * tc[1]
            -1.73956093e-08 * tc[2]
            +8.24949516e-12 * tc[3];
        /*species 7: HCO */
        species[7] =
            -3.24392532e-03
            +2.75598892e-05 * tc[1]
            -3.99432279e-08 * tc[2]
            +1.73507546e-11 * tc[3];
        /*species 8: HO2 */
        species[8] =
            -4.74912051e-03
            +4.23165782e-05 * tc[1]
            -7.28291682e-08 * tc[2]
            +3.71690050e-11 * tc[3];
        /*species 9: H2O */
        species[9] =
            -2.03643410e-03
            +1.30408042e-05 * tc[1]
            -1.64639119e-08 * tc[2]
            +7.08791268e-12 * tc[3];
        /*species 10: CO */
        species[10] =
            -6.10353680e-04
            +2.03362866e-06 * tc[1]
            +2.72101765e-09 * tc[2]
            -3.61769800e-12 * tc[3];
        /*species 11: O2 */
        species[11] =
            -2.99673416e-03
            +1.96946040e-05 * tc[1]
            -2.90438853e-08 * tc[2]
            +1.29749135e-11 * tc[3];
        /*species 12: H2O2 */
        species[12] =
            -5.42822417e-04
            +3.34671402e-05 * tc[1]
            -6.47312439e-08 * tc[2]
            +3.44981745e-11 * tc[3];
        /*species 13: CO2 */
        species[13] =
            +8.98459677e-03
            -1.42471254e-05 * tc[1]
            +7.37757066e-09 * tc[2]
            -5.74798192e-13 * tc[3];
    } else {
        /*species 0: H2 */
        species[0] =
            -4.94024731e-05
            +9.98913556e-07 * tc[1]
            -5.38699182e-10 * tc[2]
            +8.01021504e-14 * tc[3];
        /*species 1: H */
        species[1] =
            -2.30842973e-11
            +3.23123896e-14 * tc[1]
            -1.42054571e-17 * tc[2]
            +1.99278943e-21 * tc[3];
        /*species 2: AR */
        species[2] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
        /*species 3: N2 */
        species[3] =
            +1.48797680e-03
            -1.13695200e-06 * tc[1]
            +3.02911140e-10 * tc[2]
            -2.70134040e-14 * tc[3];
        /*species 4: HE */
        species[4] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
        /*species 5: O */
        species[5] =
            -8.59741137e-05
            +8.38969178e-08 * tc[1]
            -3.00533397e-11 * tc[2]
            +4.91334764e-15 * tc[3];
        /*species 6: OH */
        species[6] =
            +1.05650448e-03
            -5.18165516e-07 * tc[1]
            +9.15656022e-11 * tc[2]
            -5.32783504e-15 * tc[3];
        /*species 7: HCO */
        species[7] =
            +4.95695526e-03
            -4.96891226e-06 * tc[1]
            +1.76748533e-09 * tc[2]
            -2.13403484e-13 * tc[3];
        /*species 8: HO2 */
        species[8] =
            +2.23982013e-03
            -1.26731630e-06 * tc[1]
            +3.42739110e-10 * tc[2]
            -4.31634140e-14 * tc[3];
        /*species 9: H2O */
        species[9] =
            +2.17691804e-03
            -3.28145036e-07 * tc[1]
            -2.91125961e-10 * tc[2]
            +6.72803968e-14 * tc[3];
        /*species 10: CO */
        species[10] =
            +2.06252743e-03
            -1.99765154e-06 * tc[1]
            +6.90159024e-10 * tc[2]
            -8.14590864e-14 * tc[3];
        /*species 11: O2 */
        species[11] =
            +1.48308754e-03
            -1.51593334e-06 * tc[1]
            +6.28411665e-10 * tc[2]
            -8.66871176e-14 * tc[3];
        /*species 12: H2O2 */
        species[12] =
            +4.90831694e-03
            -3.80278450e-06 * tc[1]
            +1.11355796e-09 * tc[2]
            -1.15163322e-13 * tc[3];
        /*species 13: CO2 */
        species[13] =
            +4.41437026e-03
            -4.42962808e-06 * tc[1]
            +1.57047056e-09 * tc[2]
            -1.88833666e-13 * tc[3];
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

    double q_f[38], q_r[38];
    comp_qfqr(q_f, q_r, sc, tc, invT);

    for (int i = 0; i < 38; ++i) {
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
    kc[0] = 1.0 / (refC) * exp((g_RT[1] + g_RT[11]) - (g_RT[8]));

    /*reaction 2: OH + OH (+M) <=> H2O2 (+M) */
    kc[1] = 1.0 / (refC) * exp((g_RT[6] + g_RT[6]) - (g_RT[12]));

    /*reaction 3: CO + O (+M) <=> CO2 (+M) */
    kc[2] = 1.0 / (refC) * exp((g_RT[10] + g_RT[5]) - (g_RT[13]));

    /*reaction 4: H + H + M <=> H2 + M */
    kc[3] = 1.0 / (refC) * exp((g_RT[1] + g_RT[1]) - (g_RT[0]));

    /*reaction 5: H + OH + M <=> H2O + M */
    kc[4] = 1.0 / (refC) * exp((g_RT[1] + g_RT[6]) - (g_RT[9]));

    /*reaction 6: O + H + M <=> OH + M */
    kc[5] = 1.0 / (refC) * exp((g_RT[5] + g_RT[1]) - (g_RT[6]));

    /*reaction 7: O + O + M <=> O2 + M */
    kc[6] = 1.0 / (refC) * exp((g_RT[5] + g_RT[5]) - (g_RT[11]));

    /*reaction 8: HCO + M <=> CO + H + M */
    kc[7] = refC * exp((g_RT[7]) - (g_RT[10] + g_RT[1]));

    /*reaction 9: H + O2 <=> O + OH */
    kc[8] = exp((g_RT[1] + g_RT[11]) - (g_RT[5] + g_RT[6]));

    /*reaction 10: O + H2 <=> H + OH */
    kc[9] = exp((g_RT[5] + g_RT[0]) - (g_RT[1] + g_RT[6]));

    /*reaction 11: OH + H2 <=> H + H2O */
    kc[10] = exp((g_RT[6] + g_RT[0]) - (g_RT[1] + g_RT[9]));

    /*reaction 12: OH + OH <=> O + H2O */
    kc[11] = exp((g_RT[6] + g_RT[6]) - (g_RT[5] + g_RT[9]));

    /*reaction 13: H + H + H2 <=> H2 + H2 */
    kc[12] = 1.0 / (refC) * exp((g_RT[1] + g_RT[1] + g_RT[0]) - (g_RT[0] + g_RT[0]));

    /*reaction 14: H + H + H2O <=> H2 + H2O */
    kc[13] = 1.0 / (refC) * exp((g_RT[1] + g_RT[1] + g_RT[9]) - (g_RT[0] + g_RT[9]));

    /*reaction 15: H + H + CO2 <=> H2 + CO2 */
    kc[14] = 1.0 / (refC) * exp((g_RT[1] + g_RT[1] + g_RT[13]) - (g_RT[0] + g_RT[13]));

    /*reaction 16: H2 + O2 <=> HO2 + H */
    kc[15] = exp((g_RT[0] + g_RT[11]) - (g_RT[8] + g_RT[1]));

    /*reaction 17: HO2 + H <=> O + H2O */
    kc[16] = exp((g_RT[8] + g_RT[1]) - (g_RT[5] + g_RT[9]));

    /*reaction 18: HO2 + H <=> OH + OH */
    kc[17] = exp((g_RT[8] + g_RT[1]) - (g_RT[6] + g_RT[6]));

    /*reaction 19: HO2 + O <=> OH + O2 */
    kc[18] = exp((g_RT[8] + g_RT[5]) - (g_RT[6] + g_RT[11]));

    /*reaction 20: HO2 + OH <=> O2 + H2O */
    kc[19] = exp((g_RT[8] + g_RT[6]) - (g_RT[11] + g_RT[9]));

    /*reaction 21: HO2 + OH <=> O2 + H2O */
    kc[20] = exp((g_RT[8] + g_RT[6]) - (g_RT[11] + g_RT[9]));

    /*reaction 22: HO2 + HO2 <=> O2 + H2O2 */
    kc[21] = exp((g_RT[8] + g_RT[8]) - (g_RT[11] + g_RT[12]));

    /*reaction 23: HO2 + HO2 <=> O2 + H2O2 */
    kc[22] = exp((g_RT[8] + g_RT[8]) - (g_RT[11] + g_RT[12]));

    /*reaction 24: H2O2 + H <=> HO2 + H2 */
    kc[23] = exp((g_RT[12] + g_RT[1]) - (g_RT[8] + g_RT[0]));

    /*reaction 25: H2O2 + H <=> OH + H2O */
    kc[24] = exp((g_RT[12] + g_RT[1]) - (g_RT[6] + g_RT[9]));

    /*reaction 26: H2O2 + O <=> OH + HO2 */
    kc[25] = exp((g_RT[12] + g_RT[5]) - (g_RT[6] + g_RT[8]));

    /*reaction 27: H2O2 + OH <=> HO2 + H2O */
    kc[26] = exp((g_RT[12] + g_RT[6]) - (g_RT[8] + g_RT[9]));

    /*reaction 28: H2O2 + OH <=> HO2 + H2O */
    kc[27] = exp((g_RT[12] + g_RT[6]) - (g_RT[8] + g_RT[9]));

    /*reaction 29: CO + OH <=> CO2 + H */
    kc[28] = exp((g_RT[10] + g_RT[6]) - (g_RT[13] + g_RT[1]));

    /*reaction 30: CO + OH <=> CO2 + H */
    kc[29] = exp((g_RT[10] + g_RT[6]) - (g_RT[13] + g_RT[1]));

    /*reaction 31: CO + O2 <=> CO2 + O */
    kc[30] = exp((g_RT[10] + g_RT[11]) - (g_RT[13] + g_RT[5]));

    /*reaction 32: CO + HO2 <=> CO2 + OH */
    kc[31] = exp((g_RT[10] + g_RT[8]) - (g_RT[13] + g_RT[6]));

    /*reaction 33: HCO + H <=> CO + H2 */
    kc[32] = exp((g_RT[7] + g_RT[1]) - (g_RT[10] + g_RT[0]));

    /*reaction 34: HCO + O <=> CO + OH */
    kc[33] = exp((g_RT[7] + g_RT[5]) - (g_RT[10] + g_RT[6]));

    /*reaction 35: HCO + O <=> CO2 + H */
    kc[34] = exp((g_RT[7] + g_RT[5]) - (g_RT[13] + g_RT[1]));

    /*reaction 36: HCO + OH <=> CO + H2O */
    kc[35] = exp((g_RT[7] + g_RT[6]) - (g_RT[10] + g_RT[9]));

    /*reaction 37: HCO + H2O <=> CO + H + H2O */
    kc[36] = refC * exp((g_RT[7] + g_RT[9]) - (g_RT[10] + g_RT[1] + g_RT[9]));

    /*reaction 38: HCO + O2 <=> CO + HO2 */
    kc[37] = exp((g_RT[7] + g_RT[11]) - (g_RT[10] + g_RT[8]));

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
        /*species 0: H2 */
        species[0] =
            -9.179351730000000e+02 * invT
            +1.661320882000000e+00
            -2.344331120000000e+00 * tc[0]
            -3.990260375000000e-03 * tc[1]
            +3.246358500000000e-06 * tc[2]
            -1.679767450000000e-09 * tc[3]
            +3.688058805000000e-13 * tc[4];
        /*species 1: H */
        species[1] =
            +2.547365990000000e+04 * invT
            +2.946682853000000e+00
            -2.500000000000000e+00 * tc[0]
            -3.526664095000000e-13 * tc[1]
            +3.326532733333333e-16 * tc[2]
            -1.917346933333333e-19 * tc[3]
            +4.638661660000000e-23 * tc[4];
        /*species 2: AR */
        species[2] =
            -7.453750000000000e+02 * invT
            -1.866000000000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 3: N2 */
        species[3] =
            -1.020899900000000e+03 * invT
            -6.516950000000001e-01
            -3.298677000000000e+00 * tc[0]
            -7.041202000000000e-04 * tc[1]
            +6.605369999999999e-07 * tc[2]
            -4.701262500000001e-10 * tc[3]
            +1.222427000000000e-13 * tc[4];
        /*species 4: HE */
        species[4] =
            -7.453750000000000e+02 * invT
            +1.571276026000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 5: O */
        species[5] =
            +2.912225920000000e+04 * invT
            +1.116333640000000e+00
            -3.168267100000000e+00 * tc[0]
            +1.639659420000000e-03 * tc[1]
            -1.107177326666667e-06 * tc[2]
            +5.106721866666666e-10 * tc[3]
            -1.056329855000000e-13 * tc[4];
        /*species 6: OH */
        species[6] =
            +3.381538120000000e+03 * invT
            +4.815738570000000e+00
            -4.125305610000000e+00 * tc[0]
            +1.612724695000000e-03 * tc[1]
            -1.087941151666667e-06 * tc[2]
            +4.832113691666666e-10 * tc[3]
            -1.031186895000000e-13 * tc[4];
        /*species 7: HCO */
        species[7] =
            +3.839564960000000e+03 * invT
            +8.268134100000002e-01
            -4.221185840000000e+00 * tc[0]
            +1.621962660000000e-03 * tc[1]
            -2.296657433333333e-06 * tc[2]
            +1.109534108333333e-09 * tc[3]
            -2.168844325000000e-13 * tc[4];
        /*species 8: HO2 */
        species[8] =
            +2.948080400000000e+02 * invT
            +5.851355599999999e-01
            -4.301798010000000e+00 * tc[0]
            +2.374560255000000e-03 * tc[1]
            -3.526381516666666e-06 * tc[2]
            +2.023032450000000e-09 * tc[3]
            -4.646125620000001e-13 * tc[4];
        /*species 9: H2O */
        species[9] =
            -3.029372670000000e+04 * invT
            +5.047672768000000e+00
            -4.198640560000000e+00 * tc[0]
            +1.018217050000000e-03 * tc[1]
            -1.086733685000000e-06 * tc[2]
            +4.573308850000000e-10 * tc[3]
            -8.859890850000000e-14 * tc[4];
        /*species 10: CO */
        species[10] =
            -1.434408600000000e+04 * invT
            +7.112418999999992e-02
            -3.579533470000000e+00 * tc[0]
            +3.051768400000000e-04 * tc[1]
            -1.694690550000000e-07 * tc[2]
            -7.558382366666667e-11 * tc[3]
            +4.522122495000000e-14 * tc[4];
        /*species 11: O2 */
        species[11] =
            -1.063943560000000e+03 * invT
            +1.247806300000001e-01
            -3.782456360000000e+00 * tc[0]
            +1.498367080000000e-03 * tc[1]
            -1.641217001666667e-06 * tc[2]
            +8.067745908333334e-10 * tc[3]
            -1.621864185000000e-13 * tc[4];
        /*species 12: H2O2 */
        species[12] =
            -1.770258210000000e+04 * invT
            +8.410619499999998e-01
            -4.276112690000000e+00 * tc[0]
            +2.714112085000000e-04 * tc[1]
            -2.788928350000000e-06 * tc[2]
            +1.798090108333333e-09 * tc[3]
            -4.312271815000000e-13 * tc[4];
        /*species 13: CO2 */
        species[13] =
            -4.837196970000000e+04 * invT
            -7.544278700000000e+00
            -2.356773520000000e+00 * tc[0]
            -4.492298385000000e-03 * tc[1]
            +1.187260448333333e-06 * tc[2]
            -2.049325183333333e-10 * tc[3]
            +7.184977399999999e-15 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            -9.501589220000000e+02 * invT
            +6.542302510000000e+00
            -3.337279200000000e+00 * tc[0]
            +2.470123655000000e-05 * tc[1]
            -8.324279633333333e-08 * tc[2]
            +1.496386616666667e-11 * tc[3]
            -1.001276880000000e-15 * tc[4];
        /*species 1: H */
        species[1] =
            +2.547365990000000e+04 * invT
            +2.946682924000000e+00
            -2.500000010000000e+00 * tc[0]
            +1.154214865000000e-11 * tc[1]
            -2.692699133333334e-15 * tc[2]
            +3.945960291666667e-19 * tc[3]
            -2.490986785000000e-23 * tc[4];
        /*species 2: AR */
        species[2] =
            -7.453750000000000e+02 * invT
            -1.866000000000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 3: N2 */
        species[3] =
            -9.227977000000000e+02 * invT
            -3.053888000000000e+00
            -2.926640000000000e+00 * tc[0]
            -7.439884000000000e-04 * tc[1]
            +9.474600000000001e-08 * tc[2]
            -8.414198333333333e-12 * tc[3]
            +3.376675500000000e-16 * tc[4];
        /*species 4: HE */
        species[4] =
            -7.453750000000000e+02 * invT
            +1.571276026000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 5: O */
        species[5] =
            +2.921757910000000e+04 * invT
            -2.214917859999999e+00
            -2.569420780000000e+00 * tc[0]
            +4.298705685000000e-05 * tc[1]
            -6.991409816666667e-09 * tc[2]
            +8.348149916666666e-13 * tc[3]
            -6.141684549999999e-17 * tc[4];
        /*species 6: OH */
        species[6] =
            +3.718857740000000e+03 * invT
            -2.836911870000000e+00
            -2.864728860000000e+00 * tc[0]
            -5.282522400000000e-04 * tc[1]
            +4.318045966666667e-08 * tc[2]
            -2.543488950000000e-12 * tc[3]
            +6.659793800000000e-17 * tc[4];
        /*species 7: HCO */
        species[7] =
            +4.011918150000000e+03 * invT
            -7.026170540000000e+00
            -2.772174380000000e+00 * tc[0]
            -2.478477630000000e-03 * tc[1]
            +4.140760216666667e-07 * tc[2]
            -4.909681483333334e-11 * tc[3]
            +2.667543555000000e-15 * tc[4];
        /*species 8: HO2 */
        species[8] =
            +1.118567130000000e+02 * invT
            +2.321087500000001e-01
            -4.017210900000000e+00 * tc[0]
            -1.119910065000000e-03 * tc[1]
            +1.056096916666667e-07 * tc[2]
            -9.520530833333334e-12 * tc[3]
            +5.395426750000000e-16 * tc[4];
        /*species 9: H2O */
        species[9] =
            -3.000429710000000e+04 * invT
            -1.932777610000000e+00
            -3.033992490000000e+00 * tc[0]
            -1.088459020000000e-03 * tc[1]
            +2.734541966666666e-08 * tc[2]
            +8.086832250000000e-12 * tc[3]
            -8.410049600000000e-16 * tc[4];
        /*species 10: CO */
        species[10] =
            -1.415187240000000e+04 * invT
            -5.103502110000000e+00
            -2.715185610000000e+00 * tc[0]
            -1.031263715000000e-03 * tc[1]
            +1.664709618333334e-07 * tc[2]
            -1.917108400000000e-11 * tc[3]
            +1.018238580000000e-15 * tc[4];
        /*species 11: O2 */
        species[11] =
            -1.088457720000000e+03 * invT
            -2.170693450000000e+00
            -3.282537840000000e+00 * tc[0]
            -7.415437700000000e-04 * tc[1]
            +1.263277781666667e-07 * tc[2]
            -1.745587958333333e-11 * tc[3]
            +1.083588970000000e-15 * tc[4];
        /*species 12: H2O2 */
        species[12] =
            -1.786178770000000e+04 * invT
            +1.248846229999999e+00
            -4.165002850000000e+00 * tc[0]
            -2.454158470000000e-03 * tc[1]
            +3.168987083333333e-07 * tc[2]
            -3.093216550000000e-11 * tc[3]
            +1.439541525000000e-15 * tc[4];
        /*species 13: CO2 */
        species[13] =
            -4.875916600000000e+04 * invT
            +1.585822230000000e+00
            -3.857460290000000e+00 * tc[0]
            -2.207185130000000e-03 * tc[1]
            +3.691356733333334e-07 * tc[2]
            -4.362418233333334e-11 * tc[3]
            +2.360420820000000e-15 * tc[4];
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
        /*species 0: H2 */
        species[0] =
            -9.17935173e+02 * invT
            +6.61320882e-01
            -2.34433112e+00 * tc[0]
            -3.99026037e-03 * tc[1]
            +3.24635850e-06 * tc[2]
            -1.67976745e-09 * tc[3]
            +3.68805881e-13 * tc[4];
        /*species 1: H */
        species[1] =
            +2.54736599e+04 * invT
            +1.94668285e+00
            -2.50000000e+00 * tc[0]
            -3.52666409e-13 * tc[1]
            +3.32653273e-16 * tc[2]
            -1.91734693e-19 * tc[3]
            +4.63866166e-23 * tc[4];
        /*species 2: AR */
        species[2] =
            -7.45375000e+02 * invT
            -2.86600000e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 3: N2 */
        species[3] =
            -1.02089990e+03 * invT
            -1.65169500e+00
            -3.29867700e+00 * tc[0]
            -7.04120200e-04 * tc[1]
            +6.60537000e-07 * tc[2]
            -4.70126250e-10 * tc[3]
            +1.22242700e-13 * tc[4];
        /*species 4: HE */
        species[4] =
            -7.45375000e+02 * invT
            +5.71276026e-01
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 5: O */
        species[5] =
            +2.91222592e+04 * invT
            +1.16333640e-01
            -3.16826710e+00 * tc[0]
            +1.63965942e-03 * tc[1]
            -1.10717733e-06 * tc[2]
            +5.10672187e-10 * tc[3]
            -1.05632985e-13 * tc[4];
        /*species 6: OH */
        species[6] =
            +3.38153812e+03 * invT
            +3.81573857e+00
            -4.12530561e+00 * tc[0]
            +1.61272470e-03 * tc[1]
            -1.08794115e-06 * tc[2]
            +4.83211369e-10 * tc[3]
            -1.03118689e-13 * tc[4];
        /*species 7: HCO */
        species[7] =
            +3.83956496e+03 * invT
            -1.73186590e-01
            -4.22118584e+00 * tc[0]
            +1.62196266e-03 * tc[1]
            -2.29665743e-06 * tc[2]
            +1.10953411e-09 * tc[3]
            -2.16884432e-13 * tc[4];
        /*species 8: HO2 */
        species[8] =
            +2.94808040e+02 * invT
            -4.14864440e-01
            -4.30179801e+00 * tc[0]
            +2.37456025e-03 * tc[1]
            -3.52638152e-06 * tc[2]
            +2.02303245e-09 * tc[3]
            -4.64612562e-13 * tc[4];
        /*species 9: H2O */
        species[9] =
            -3.02937267e+04 * invT
            +4.04767277e+00
            -4.19864056e+00 * tc[0]
            +1.01821705e-03 * tc[1]
            -1.08673369e-06 * tc[2]
            +4.57330885e-10 * tc[3]
            -8.85989085e-14 * tc[4];
        /*species 10: CO */
        species[10] =
            -1.43440860e+04 * invT
            -9.28875810e-01
            -3.57953347e+00 * tc[0]
            +3.05176840e-04 * tc[1]
            -1.69469055e-07 * tc[2]
            -7.55838237e-11 * tc[3]
            +4.52212249e-14 * tc[4];
        /*species 11: O2 */
        species[11] =
            -1.06394356e+03 * invT
            -8.75219370e-01
            -3.78245636e+00 * tc[0]
            +1.49836708e-03 * tc[1]
            -1.64121700e-06 * tc[2]
            +8.06774591e-10 * tc[3]
            -1.62186418e-13 * tc[4];
        /*species 12: H2O2 */
        species[12] =
            -1.77025821e+04 * invT
            -1.58938050e-01
            -4.27611269e+00 * tc[0]
            +2.71411208e-04 * tc[1]
            -2.78892835e-06 * tc[2]
            +1.79809011e-09 * tc[3]
            -4.31227182e-13 * tc[4];
        /*species 13: CO2 */
        species[13] =
            -4.83719697e+04 * invT
            -8.54427870e+00
            -2.35677352e+00 * tc[0]
            -4.49229839e-03 * tc[1]
            +1.18726045e-06 * tc[2]
            -2.04932518e-10 * tc[3]
            +7.18497740e-15 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            -9.50158922e+02 * invT
            +5.54230251e+00
            -3.33727920e+00 * tc[0]
            +2.47012365e-05 * tc[1]
            -8.32427963e-08 * tc[2]
            +1.49638662e-11 * tc[3]
            -1.00127688e-15 * tc[4];
        /*species 1: H */
        species[1] =
            +2.54736599e+04 * invT
            +1.94668292e+00
            -2.50000001e+00 * tc[0]
            +1.15421486e-11 * tc[1]
            -2.69269913e-15 * tc[2]
            +3.94596029e-19 * tc[3]
            -2.49098679e-23 * tc[4];
        /*species 2: AR */
        species[2] =
            -7.45375000e+02 * invT
            -2.86600000e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 3: N2 */
        species[3] =
            -9.22797700e+02 * invT
            -4.05388800e+00
            -2.92664000e+00 * tc[0]
            -7.43988400e-04 * tc[1]
            +9.47460000e-08 * tc[2]
            -8.41419833e-12 * tc[3]
            +3.37667550e-16 * tc[4];
        /*species 4: HE */
        species[4] =
            -7.45375000e+02 * invT
            +5.71276026e-01
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 5: O */
        species[5] =
            +2.92175791e+04 * invT
            -3.21491786e+00
            -2.56942078e+00 * tc[0]
            +4.29870569e-05 * tc[1]
            -6.99140982e-09 * tc[2]
            +8.34814992e-13 * tc[3]
            -6.14168455e-17 * tc[4];
        /*species 6: OH */
        species[6] =
            +3.71885774e+03 * invT
            -3.83691187e+00
            -2.86472886e+00 * tc[0]
            -5.28252240e-04 * tc[1]
            +4.31804597e-08 * tc[2]
            -2.54348895e-12 * tc[3]
            +6.65979380e-17 * tc[4];
        /*species 7: HCO */
        species[7] =
            +4.01191815e+03 * invT
            -8.02617054e+00
            -2.77217438e+00 * tc[0]
            -2.47847763e-03 * tc[1]
            +4.14076022e-07 * tc[2]
            -4.90968148e-11 * tc[3]
            +2.66754356e-15 * tc[4];
        /*species 8: HO2 */
        species[8] =
            +1.11856713e+02 * invT
            -7.67891250e-01
            -4.01721090e+00 * tc[0]
            -1.11991006e-03 * tc[1]
            +1.05609692e-07 * tc[2]
            -9.52053083e-12 * tc[3]
            +5.39542675e-16 * tc[4];
        /*species 9: H2O */
        species[9] =
            -3.00042971e+04 * invT
            -2.93277761e+00
            -3.03399249e+00 * tc[0]
            -1.08845902e-03 * tc[1]
            +2.73454197e-08 * tc[2]
            +8.08683225e-12 * tc[3]
            -8.41004960e-16 * tc[4];
        /*species 10: CO */
        species[10] =
            -1.41518724e+04 * invT
            -6.10350211e+00
            -2.71518561e+00 * tc[0]
            -1.03126372e-03 * tc[1]
            +1.66470962e-07 * tc[2]
            -1.91710840e-11 * tc[3]
            +1.01823858e-15 * tc[4];
        /*species 11: O2 */
        species[11] =
            -1.08845772e+03 * invT
            -3.17069345e+00
            -3.28253784e+00 * tc[0]
            -7.41543770e-04 * tc[1]
            +1.26327778e-07 * tc[2]
            -1.74558796e-11 * tc[3]
            +1.08358897e-15 * tc[4];
        /*species 12: H2O2 */
        species[12] =
            -1.78617877e+04 * invT
            +2.48846230e-01
            -4.16500285e+00 * tc[0]
            -2.45415847e-03 * tc[1]
            +3.16898708e-07 * tc[2]
            -3.09321655e-11 * tc[3]
            +1.43954153e-15 * tc[4];
        /*species 13: CO2 */
        species[13] =
            -4.87591660e+04 * invT
            +5.85822230e-01
            -3.85746029e+00 * tc[0]
            -2.20718513e-03 * tc[1]
            +3.69135673e-07 * tc[2]
            -4.36241823e-11 * tc[3]
            +2.36042082e-15 * tc[4];
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
        /*species 0: H2 */
        species[0] =
            +1.34433112e+00
            +7.98052075e-03 * tc[1]
            -1.94781510e-05 * tc[2]
            +2.01572094e-08 * tc[3]
            -7.37611761e-12 * tc[4];
        /*species 1: H */
        species[1] =
            +1.50000000e+00
            +7.05332819e-13 * tc[1]
            -1.99591964e-15 * tc[2]
            +2.30081632e-18 * tc[3]
            -9.27732332e-22 * tc[4];
        /*species 2: AR */
        species[2] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 3: N2 */
        species[3] =
            +2.29867700e+00
            +1.40824040e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485400e-12 * tc[4];
        /*species 4: HE */
        species[4] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 5: O */
        species[5] =
            +2.16826710e+00
            -3.27931884e-03 * tc[1]
            +6.64306396e-06 * tc[2]
            -6.12806624e-09 * tc[3]
            +2.11265971e-12 * tc[4];
        /*species 6: OH */
        species[6] =
            +3.12530561e+00
            -3.22544939e-03 * tc[1]
            +6.52764691e-06 * tc[2]
            -5.79853643e-09 * tc[3]
            +2.06237379e-12 * tc[4];
        /*species 7: HCO */
        species[7] =
            +3.22118584e+00
            -3.24392532e-03 * tc[1]
            +1.37799446e-05 * tc[2]
            -1.33144093e-08 * tc[3]
            +4.33768865e-12 * tc[4];
        /*species 8: HO2 */
        species[8] =
            +3.30179801e+00
            -4.74912051e-03 * tc[1]
            +2.11582891e-05 * tc[2]
            -2.42763894e-08 * tc[3]
            +9.29225124e-12 * tc[4];
        /*species 9: H2O */
        species[9] =
            +3.19864056e+00
            -2.03643410e-03 * tc[1]
            +6.52040211e-06 * tc[2]
            -5.48797062e-09 * tc[3]
            +1.77197817e-12 * tc[4];
        /*species 10: CO */
        species[10] =
            +2.57953347e+00
            -6.10353680e-04 * tc[1]
            +1.01681433e-06 * tc[2]
            +9.07005884e-10 * tc[3]
            -9.04424499e-13 * tc[4];
        /*species 11: O2 */
        species[11] =
            +2.78245636e+00
            -2.99673416e-03 * tc[1]
            +9.84730201e-06 * tc[2]
            -9.68129509e-09 * tc[3]
            +3.24372837e-12 * tc[4];
        /*species 12: H2O2 */
        species[12] =
            +3.27611269e+00
            -5.42822417e-04 * tc[1]
            +1.67335701e-05 * tc[2]
            -2.15770813e-08 * tc[3]
            +8.62454363e-12 * tc[4];
        /*species 13: CO2 */
        species[13] =
            +1.35677352e+00
            +8.98459677e-03 * tc[1]
            -7.12356269e-06 * tc[2]
            +2.45919022e-09 * tc[3]
            -1.43699548e-13 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            +2.33727920e+00
            -4.94024731e-05 * tc[1]
            +4.99456778e-07 * tc[2]
            -1.79566394e-10 * tc[3]
            +2.00255376e-14 * tc[4];
        /*species 1: H */
        species[1] =
            +1.50000001e+00
            -2.30842973e-11 * tc[1]
            +1.61561948e-14 * tc[2]
            -4.73515235e-18 * tc[3]
            +4.98197357e-22 * tc[4];
        /*species 2: AR */
        species[2] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 3: N2 */
        species[3] =
            +1.92664000e+00
            +1.48797680e-03 * tc[1]
            -5.68476000e-07 * tc[2]
            +1.00970380e-10 * tc[3]
            -6.75335100e-15 * tc[4];
        /*species 4: HE */
        species[4] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 5: O */
        species[5] =
            +1.56942078e+00
            -8.59741137e-05 * tc[1]
            +4.19484589e-08 * tc[2]
            -1.00177799e-11 * tc[3]
            +1.22833691e-15 * tc[4];
        /*species 6: OH */
        species[6] =
            +1.86472886e+00
            +1.05650448e-03 * tc[1]
            -2.59082758e-07 * tc[2]
            +3.05218674e-11 * tc[3]
            -1.33195876e-15 * tc[4];
        /*species 7: HCO */
        species[7] =
            +1.77217438e+00
            +4.95695526e-03 * tc[1]
            -2.48445613e-06 * tc[2]
            +5.89161778e-10 * tc[3]
            -5.33508711e-14 * tc[4];
        /*species 8: HO2 */
        species[8] =
            +3.01721090e+00
            +2.23982013e-03 * tc[1]
            -6.33658150e-07 * tc[2]
            +1.14246370e-10 * tc[3]
            -1.07908535e-14 * tc[4];
        /*species 9: H2O */
        species[9] =
            +2.03399249e+00
            +2.17691804e-03 * tc[1]
            -1.64072518e-07 * tc[2]
            -9.70419870e-11 * tc[3]
            +1.68200992e-14 * tc[4];
        /*species 10: CO */
        species[10] =
            +1.71518561e+00
            +2.06252743e-03 * tc[1]
            -9.98825771e-07 * tc[2]
            +2.30053008e-10 * tc[3]
            -2.03647716e-14 * tc[4];
        /*species 11: O2 */
        species[11] =
            +2.28253784e+00
            +1.48308754e-03 * tc[1]
            -7.57966669e-07 * tc[2]
            +2.09470555e-10 * tc[3]
            -2.16717794e-14 * tc[4];
        /*species 12: H2O2 */
        species[12] =
            +3.16500285e+00
            +4.90831694e-03 * tc[1]
            -1.90139225e-06 * tc[2]
            +3.71185986e-10 * tc[3]
            -2.87908305e-14 * tc[4];
        /*species 13: CO2 */
        species[13] =
            +2.85746029e+00
            +4.41437026e-03 * tc[1]
            -2.21481404e-06 * tc[2]
            +5.23490188e-10 * tc[3]
            -4.72084164e-14 * tc[4];
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
        /*species 0: H2 */
        species[0] =
            +2.34433112e+00
            +7.98052075e-03 * tc[1]
            -1.94781510e-05 * tc[2]
            +2.01572094e-08 * tc[3]
            -7.37611761e-12 * tc[4];
        /*species 1: H */
        species[1] =
            +2.50000000e+00
            +7.05332819e-13 * tc[1]
            -1.99591964e-15 * tc[2]
            +2.30081632e-18 * tc[3]
            -9.27732332e-22 * tc[4];
        /*species 2: AR */
        species[2] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 3: N2 */
        species[3] =
            +3.29867700e+00
            +1.40824040e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485400e-12 * tc[4];
        /*species 4: HE */
        species[4] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 5: O */
        species[5] =
            +3.16826710e+00
            -3.27931884e-03 * tc[1]
            +6.64306396e-06 * tc[2]
            -6.12806624e-09 * tc[3]
            +2.11265971e-12 * tc[4];
        /*species 6: OH */
        species[6] =
            +4.12530561e+00
            -3.22544939e-03 * tc[1]
            +6.52764691e-06 * tc[2]
            -5.79853643e-09 * tc[3]
            +2.06237379e-12 * tc[4];
        /*species 7: HCO */
        species[7] =
            +4.22118584e+00
            -3.24392532e-03 * tc[1]
            +1.37799446e-05 * tc[2]
            -1.33144093e-08 * tc[3]
            +4.33768865e-12 * tc[4];
        /*species 8: HO2 */
        species[8] =
            +4.30179801e+00
            -4.74912051e-03 * tc[1]
            +2.11582891e-05 * tc[2]
            -2.42763894e-08 * tc[3]
            +9.29225124e-12 * tc[4];
        /*species 9: H2O */
        species[9] =
            +4.19864056e+00
            -2.03643410e-03 * tc[1]
            +6.52040211e-06 * tc[2]
            -5.48797062e-09 * tc[3]
            +1.77197817e-12 * tc[4];
        /*species 10: CO */
        species[10] =
            +3.57953347e+00
            -6.10353680e-04 * tc[1]
            +1.01681433e-06 * tc[2]
            +9.07005884e-10 * tc[3]
            -9.04424499e-13 * tc[4];
        /*species 11: O2 */
        species[11] =
            +3.78245636e+00
            -2.99673416e-03 * tc[1]
            +9.84730201e-06 * tc[2]
            -9.68129509e-09 * tc[3]
            +3.24372837e-12 * tc[4];
        /*species 12: H2O2 */
        species[12] =
            +4.27611269e+00
            -5.42822417e-04 * tc[1]
            +1.67335701e-05 * tc[2]
            -2.15770813e-08 * tc[3]
            +8.62454363e-12 * tc[4];
        /*species 13: CO2 */
        species[13] =
            +2.35677352e+00
            +8.98459677e-03 * tc[1]
            -7.12356269e-06 * tc[2]
            +2.45919022e-09 * tc[3]
            -1.43699548e-13 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            +3.33727920e+00
            -4.94024731e-05 * tc[1]
            +4.99456778e-07 * tc[2]
            -1.79566394e-10 * tc[3]
            +2.00255376e-14 * tc[4];
        /*species 1: H */
        species[1] =
            +2.50000001e+00
            -2.30842973e-11 * tc[1]
            +1.61561948e-14 * tc[2]
            -4.73515235e-18 * tc[3]
            +4.98197357e-22 * tc[4];
        /*species 2: AR */
        species[2] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 3: N2 */
        species[3] =
            +2.92664000e+00
            +1.48797680e-03 * tc[1]
            -5.68476000e-07 * tc[2]
            +1.00970380e-10 * tc[3]
            -6.75335100e-15 * tc[4];
        /*species 4: HE */
        species[4] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 5: O */
        species[5] =
            +2.56942078e+00
            -8.59741137e-05 * tc[1]
            +4.19484589e-08 * tc[2]
            -1.00177799e-11 * tc[3]
            +1.22833691e-15 * tc[4];
        /*species 6: OH */
        species[6] =
            +2.86472886e+00
            +1.05650448e-03 * tc[1]
            -2.59082758e-07 * tc[2]
            +3.05218674e-11 * tc[3]
            -1.33195876e-15 * tc[4];
        /*species 7: HCO */
        species[7] =
            +2.77217438e+00
            +4.95695526e-03 * tc[1]
            -2.48445613e-06 * tc[2]
            +5.89161778e-10 * tc[3]
            -5.33508711e-14 * tc[4];
        /*species 8: HO2 */
        species[8] =
            +4.01721090e+00
            +2.23982013e-03 * tc[1]
            -6.33658150e-07 * tc[2]
            +1.14246370e-10 * tc[3]
            -1.07908535e-14 * tc[4];
        /*species 9: H2O */
        species[9] =
            +3.03399249e+00
            +2.17691804e-03 * tc[1]
            -1.64072518e-07 * tc[2]
            -9.70419870e-11 * tc[3]
            +1.68200992e-14 * tc[4];
        /*species 10: CO */
        species[10] =
            +2.71518561e+00
            +2.06252743e-03 * tc[1]
            -9.98825771e-07 * tc[2]
            +2.30053008e-10 * tc[3]
            -2.03647716e-14 * tc[4];
        /*species 11: O2 */
        species[11] =
            +3.28253784e+00
            +1.48308754e-03 * tc[1]
            -7.57966669e-07 * tc[2]
            +2.09470555e-10 * tc[3]
            -2.16717794e-14 * tc[4];
        /*species 12: H2O2 */
        species[12] =
            +4.16500285e+00
            +4.90831694e-03 * tc[1]
            -1.90139225e-06 * tc[2]
            +3.71185986e-10 * tc[3]
            -2.87908305e-14 * tc[4];
        /*species 13: CO2 */
        species[13] =
            +3.85746029e+00
            +4.41437026e-03 * tc[1]
            -2.21481404e-06 * tc[2]
            +5.23490188e-10 * tc[3]
            -4.72084164e-14 * tc[4];
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
        /*species 0: H2 */
        species[0] =
            +1.34433112e+00
            +3.99026037e-03 * tc[1]
            -6.49271700e-06 * tc[2]
            +5.03930235e-09 * tc[3]
            -1.47522352e-12 * tc[4]
            -9.17935173e+02 * invT;
        /*species 1: H */
        species[1] =
            +1.50000000e+00
            +3.52666409e-13 * tc[1]
            -6.65306547e-16 * tc[2]
            +5.75204080e-19 * tc[3]
            -1.85546466e-22 * tc[4]
            +2.54736599e+04 * invT;
        /*species 2: AR */
        species[2] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
        /*species 3: N2 */
        species[3] =
            +2.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 * invT;
        /*species 4: HE */
        species[4] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
        /*species 5: O */
        species[5] =
            +2.16826710e+00
            -1.63965942e-03 * tc[1]
            +2.21435465e-06 * tc[2]
            -1.53201656e-09 * tc[3]
            +4.22531942e-13 * tc[4]
            +2.91222592e+04 * invT;
        /*species 6: OH */
        species[6] =
            +3.12530561e+00
            -1.61272470e-03 * tc[1]
            +2.17588230e-06 * tc[2]
            -1.44963411e-09 * tc[3]
            +4.12474758e-13 * tc[4]
            +3.38153812e+03 * invT;
        /*species 7: HCO */
        species[7] =
            +3.22118584e+00
            -1.62196266e-03 * tc[1]
            +4.59331487e-06 * tc[2]
            -3.32860233e-09 * tc[3]
            +8.67537730e-13 * tc[4]
            +3.83956496e+03 * invT;
        /*species 8: HO2 */
        species[8] =
            +3.30179801e+00
            -2.37456025e-03 * tc[1]
            +7.05276303e-06 * tc[2]
            -6.06909735e-09 * tc[3]
            +1.85845025e-12 * tc[4]
            +2.94808040e+02 * invT;
        /*species 9: H2O */
        species[9] =
            +3.19864056e+00
            -1.01821705e-03 * tc[1]
            +2.17346737e-06 * tc[2]
            -1.37199266e-09 * tc[3]
            +3.54395634e-13 * tc[4]
            -3.02937267e+04 * invT;
        /*species 10: CO */
        species[10] =
            +2.57953347e+00
            -3.05176840e-04 * tc[1]
            +3.38938110e-07 * tc[2]
            +2.26751471e-10 * tc[3]
            -1.80884900e-13 * tc[4]
            -1.43440860e+04 * invT;
        /*species 11: O2 */
        species[11] =
            +2.78245636e+00
            -1.49836708e-03 * tc[1]
            +3.28243400e-06 * tc[2]
            -2.42032377e-09 * tc[3]
            +6.48745674e-13 * tc[4]
            -1.06394356e+03 * invT;
        /*species 12: H2O2 */
        species[12] =
            +3.27611269e+00
            -2.71411208e-04 * tc[1]
            +5.57785670e-06 * tc[2]
            -5.39427032e-09 * tc[3]
            +1.72490873e-12 * tc[4]
            -1.77025821e+04 * invT;
        /*species 13: CO2 */
        species[13] =
            +1.35677352e+00
            +4.49229839e-03 * tc[1]
            -2.37452090e-06 * tc[2]
            +6.14797555e-10 * tc[3]
            -2.87399096e-14 * tc[4]
            -4.83719697e+04 * invT;
    } else {
        /*species 0: H2 */
        species[0] =
            +2.33727920e+00
            -2.47012365e-05 * tc[1]
            +1.66485593e-07 * tc[2]
            -4.48915985e-11 * tc[3]
            +4.00510752e-15 * tc[4]
            -9.50158922e+02 * invT;
        /*species 1: H */
        species[1] =
            +1.50000001e+00
            -1.15421486e-11 * tc[1]
            +5.38539827e-15 * tc[2]
            -1.18378809e-18 * tc[3]
            +9.96394714e-23 * tc[4]
            +2.54736599e+04 * invT;
        /*species 2: AR */
        species[2] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
        /*species 3: N2 */
        species[3] =
            +1.92664000e+00
            +7.43988400e-04 * tc[1]
            -1.89492000e-07 * tc[2]
            +2.52425950e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
        /*species 4: HE */
        species[4] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
        /*species 5: O */
        species[5] =
            +1.56942078e+00
            -4.29870569e-05 * tc[1]
            +1.39828196e-08 * tc[2]
            -2.50444497e-12 * tc[3]
            +2.45667382e-16 * tc[4]
            +2.92175791e+04 * invT;
        /*species 6: OH */
        species[6] =
            +1.86472886e+00
            +5.28252240e-04 * tc[1]
            -8.63609193e-08 * tc[2]
            +7.63046685e-12 * tc[3]
            -2.66391752e-16 * tc[4]
            +3.71885774e+03 * invT;
        /*species 7: HCO */
        species[7] =
            +1.77217438e+00
            +2.47847763e-03 * tc[1]
            -8.28152043e-07 * tc[2]
            +1.47290445e-10 * tc[3]
            -1.06701742e-14 * tc[4]
            +4.01191815e+03 * invT;
        /*species 8: HO2 */
        species[8] =
            +3.01721090e+00
            +1.11991006e-03 * tc[1]
            -2.11219383e-07 * tc[2]
            +2.85615925e-11 * tc[3]
            -2.15817070e-15 * tc[4]
            +1.11856713e+02 * invT;
        /*species 9: H2O */
        species[9] =
            +2.03399249e+00
            +1.08845902e-03 * tc[1]
            -5.46908393e-08 * tc[2]
            -2.42604967e-11 * tc[3]
            +3.36401984e-15 * tc[4]
            -3.00042971e+04 * invT;
        /*species 10: CO */
        species[10] =
            +1.71518561e+00
            +1.03126372e-03 * tc[1]
            -3.32941924e-07 * tc[2]
            +5.75132520e-11 * tc[3]
            -4.07295432e-15 * tc[4]
            -1.41518724e+04 * invT;
        /*species 11: O2 */
        species[11] =
            +2.28253784e+00
            +7.41543770e-04 * tc[1]
            -2.52655556e-07 * tc[2]
            +5.23676387e-11 * tc[3]
            -4.33435588e-15 * tc[4]
            -1.08845772e+03 * invT;
        /*species 12: H2O2 */
        species[12] =
            +3.16500285e+00
            +2.45415847e-03 * tc[1]
            -6.33797417e-07 * tc[2]
            +9.27964965e-11 * tc[3]
            -5.75816610e-15 * tc[4]
            -1.78617877e+04 * invT;
        /*species 13: CO2 */
        species[13] =
            +2.85746029e+00
            +2.20718513e-03 * tc[1]
            -7.38271347e-07 * tc[2]
            +1.30872547e-10 * tc[3]
            -9.44168328e-15 * tc[4]
            -4.87591660e+04 * invT;
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
        /*species 0: H2 */
        species[0] =
            +2.34433112e+00
            +3.99026037e-03 * tc[1]
            -6.49271700e-06 * tc[2]
            +5.03930235e-09 * tc[3]
            -1.47522352e-12 * tc[4]
            -9.17935173e+02 * invT;
        /*species 1: H */
        species[1] =
            +2.50000000e+00
            +3.52666409e-13 * tc[1]
            -6.65306547e-16 * tc[2]
            +5.75204080e-19 * tc[3]
            -1.85546466e-22 * tc[4]
            +2.54736599e+04 * invT;
        /*species 2: AR */
        species[2] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
        /*species 3: N2 */
        species[3] =
            +3.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 * invT;
        /*species 4: HE */
        species[4] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
        /*species 5: O */
        species[5] =
            +3.16826710e+00
            -1.63965942e-03 * tc[1]
            +2.21435465e-06 * tc[2]
            -1.53201656e-09 * tc[3]
            +4.22531942e-13 * tc[4]
            +2.91222592e+04 * invT;
        /*species 6: OH */
        species[6] =
            +4.12530561e+00
            -1.61272470e-03 * tc[1]
            +2.17588230e-06 * tc[2]
            -1.44963411e-09 * tc[3]
            +4.12474758e-13 * tc[4]
            +3.38153812e+03 * invT;
        /*species 7: HCO */
        species[7] =
            +4.22118584e+00
            -1.62196266e-03 * tc[1]
            +4.59331487e-06 * tc[2]
            -3.32860233e-09 * tc[3]
            +8.67537730e-13 * tc[4]
            +3.83956496e+03 * invT;
        /*species 8: HO2 */
        species[8] =
            +4.30179801e+00
            -2.37456025e-03 * tc[1]
            +7.05276303e-06 * tc[2]
            -6.06909735e-09 * tc[3]
            +1.85845025e-12 * tc[4]
            +2.94808040e+02 * invT;
        /*species 9: H2O */
        species[9] =
            +4.19864056e+00
            -1.01821705e-03 * tc[1]
            +2.17346737e-06 * tc[2]
            -1.37199266e-09 * tc[3]
            +3.54395634e-13 * tc[4]
            -3.02937267e+04 * invT;
        /*species 10: CO */
        species[10] =
            +3.57953347e+00
            -3.05176840e-04 * tc[1]
            +3.38938110e-07 * tc[2]
            +2.26751471e-10 * tc[3]
            -1.80884900e-13 * tc[4]
            -1.43440860e+04 * invT;
        /*species 11: O2 */
        species[11] =
            +3.78245636e+00
            -1.49836708e-03 * tc[1]
            +3.28243400e-06 * tc[2]
            -2.42032377e-09 * tc[3]
            +6.48745674e-13 * tc[4]
            -1.06394356e+03 * invT;
        /*species 12: H2O2 */
        species[12] =
            +4.27611269e+00
            -2.71411208e-04 * tc[1]
            +5.57785670e-06 * tc[2]
            -5.39427032e-09 * tc[3]
            +1.72490873e-12 * tc[4]
            -1.77025821e+04 * invT;
        /*species 13: CO2 */
        species[13] =
            +2.35677352e+00
            +4.49229839e-03 * tc[1]
            -2.37452090e-06 * tc[2]
            +6.14797555e-10 * tc[3]
            -2.87399096e-14 * tc[4]
            -4.83719697e+04 * invT;
    } else {
        /*species 0: H2 */
        species[0] =
            +3.33727920e+00
            -2.47012365e-05 * tc[1]
            +1.66485593e-07 * tc[2]
            -4.48915985e-11 * tc[3]
            +4.00510752e-15 * tc[4]
            -9.50158922e+02 * invT;
        /*species 1: H */
        species[1] =
            +2.50000001e+00
            -1.15421486e-11 * tc[1]
            +5.38539827e-15 * tc[2]
            -1.18378809e-18 * tc[3]
            +9.96394714e-23 * tc[4]
            +2.54736599e+04 * invT;
        /*species 2: AR */
        species[2] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
        /*species 3: N2 */
        species[3] =
            +2.92664000e+00
            +7.43988400e-04 * tc[1]
            -1.89492000e-07 * tc[2]
            +2.52425950e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
        /*species 4: HE */
        species[4] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -7.45375000e+02 * invT;
        /*species 5: O */
        species[5] =
            +2.56942078e+00
            -4.29870569e-05 * tc[1]
            +1.39828196e-08 * tc[2]
            -2.50444497e-12 * tc[3]
            +2.45667382e-16 * tc[4]
            +2.92175791e+04 * invT;
        /*species 6: OH */
        species[6] =
            +2.86472886e+00
            +5.28252240e-04 * tc[1]
            -8.63609193e-08 * tc[2]
            +7.63046685e-12 * tc[3]
            -2.66391752e-16 * tc[4]
            +3.71885774e+03 * invT;
        /*species 7: HCO */
        species[7] =
            +2.77217438e+00
            +2.47847763e-03 * tc[1]
            -8.28152043e-07 * tc[2]
            +1.47290445e-10 * tc[3]
            -1.06701742e-14 * tc[4]
            +4.01191815e+03 * invT;
        /*species 8: HO2 */
        species[8] =
            +4.01721090e+00
            +1.11991006e-03 * tc[1]
            -2.11219383e-07 * tc[2]
            +2.85615925e-11 * tc[3]
            -2.15817070e-15 * tc[4]
            +1.11856713e+02 * invT;
        /*species 9: H2O */
        species[9] =
            +3.03399249e+00
            +1.08845902e-03 * tc[1]
            -5.46908393e-08 * tc[2]
            -2.42604967e-11 * tc[3]
            +3.36401984e-15 * tc[4]
            -3.00042971e+04 * invT;
        /*species 10: CO */
        species[10] =
            +2.71518561e+00
            +1.03126372e-03 * tc[1]
            -3.32941924e-07 * tc[2]
            +5.75132520e-11 * tc[3]
            -4.07295432e-15 * tc[4]
            -1.41518724e+04 * invT;
        /*species 11: O2 */
        species[11] =
            +3.28253784e+00
            +7.41543770e-04 * tc[1]
            -2.52655556e-07 * tc[2]
            +5.23676387e-11 * tc[3]
            -4.33435588e-15 * tc[4]
            -1.08845772e+03 * invT;
        /*species 12: H2O2 */
        species[12] =
            +4.16500285e+00
            +2.45415847e-03 * tc[1]
            -6.33797417e-07 * tc[2]
            +9.27964965e-11 * tc[3]
            -5.75816610e-15 * tc[4]
            -1.78617877e+04 * invT;
        /*species 13: CO2 */
        species[13] =
            +3.85746029e+00
            +2.20718513e-03 * tc[1]
            -7.38271347e-07 * tc[2]
            +1.30872547e-10 * tc[3]
            -9.44168328e-15 * tc[4]
            -4.87591660e+04 * invT;
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
        /*species 0: H2 */
        species[0] =
            +2.34433112e+00 * tc[0]
            +7.98052075e-03 * tc[1]
            -9.73907550e-06 * tc[2]
            +6.71906980e-09 * tc[3]
            -1.84402940e-12 * tc[4]
            +6.83010238e-01 ;
        /*species 1: H */
        species[1] =
            +2.50000000e+00 * tc[0]
            +7.05332819e-13 * tc[1]
            -9.97959820e-16 * tc[2]
            +7.66938773e-19 * tc[3]
            -2.31933083e-22 * tc[4]
            -4.46682853e-01 ;
        /*species 2: AR */
        species[2] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +4.36600000e+00 ;
        /*species 3: N2 */
        species[3] =
            +3.29867700e+00 * tc[0]
            +1.40824040e-03 * tc[1]
            -1.98161100e-06 * tc[2]
            +1.88050500e-09 * tc[3]
            -6.11213500e-13 * tc[4]
            +3.95037200e+00 ;
        /*species 4: HE */
        species[4] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +9.28723974e-01 ;
        /*species 5: O */
        species[5] =
            +3.16826710e+00 * tc[0]
            -3.27931884e-03 * tc[1]
            +3.32153198e-06 * tc[2]
            -2.04268875e-09 * tc[3]
            +5.28164927e-13 * tc[4]
            +2.05193346e+00 ;
        /*species 6: OH */
        species[6] =
            +4.12530561e+00 * tc[0]
            -3.22544939e-03 * tc[1]
            +3.26382346e-06 * tc[2]
            -1.93284548e-09 * tc[3]
            +5.15593447e-13 * tc[4]
            -6.90432960e-01 ;
        /*species 7: HCO */
        species[7] =
            +4.22118584e+00 * tc[0]
            -3.24392532e-03 * tc[1]
            +6.88997230e-06 * tc[2]
            -4.43813643e-09 * tc[3]
            +1.08442216e-12 * tc[4]
            +3.39437243e+00 ;
        /*species 8: HO2 */
        species[8] =
            +4.30179801e+00 * tc[0]
            -4.74912051e-03 * tc[1]
            +1.05791445e-05 * tc[2]
            -8.09212980e-09 * tc[3]
            +2.32306281e-12 * tc[4]
            +3.71666245e+00 ;
        /*species 9: H2O */
        species[9] =
            +4.19864056e+00 * tc[0]
            -2.03643410e-03 * tc[1]
            +3.26020105e-06 * tc[2]
            -1.82932354e-09 * tc[3]
            +4.42994543e-13 * tc[4]
            -8.49032208e-01 ;
        /*species 10: CO */
        species[10] =
            +3.57953347e+00 * tc[0]
            -6.10353680e-04 * tc[1]
            +5.08407165e-07 * tc[2]
            +3.02335295e-10 * tc[3]
            -2.26106125e-13 * tc[4]
            +3.50840928e+00 ;
        /*species 11: O2 */
        species[11] =
            +3.78245636e+00 * tc[0]
            -2.99673416e-03 * tc[1]
            +4.92365101e-06 * tc[2]
            -3.22709836e-09 * tc[3]
            +8.10932092e-13 * tc[4]
            +3.65767573e+00 ;
        /*species 12: H2O2 */
        species[12] =
            +4.27611269e+00 * tc[0]
            -5.42822417e-04 * tc[1]
            +8.36678505e-06 * tc[2]
            -7.19236043e-09 * tc[3]
            +2.15613591e-12 * tc[4]
            +3.43505074e+00 ;
        /*species 13: CO2 */
        species[13] =
            +2.35677352e+00 * tc[0]
            +8.98459677e-03 * tc[1]
            -3.56178134e-06 * tc[2]
            +8.19730073e-10 * tc[3]
            -3.59248870e-14 * tc[4]
            +9.90105222e+00 ;
    } else {
        /*species 0: H2 */
        species[0] =
            +3.33727920e+00 * tc[0]
            -4.94024731e-05 * tc[1]
            +2.49728389e-07 * tc[2]
            -5.98554647e-11 * tc[3]
            +5.00638440e-15 * tc[4]
            -3.20502331e+00 ;
        /*species 1: H */
        species[1] =
            +2.50000001e+00 * tc[0]
            -2.30842973e-11 * tc[1]
            +8.07809740e-15 * tc[2]
            -1.57838412e-18 * tc[3]
            +1.24549339e-22 * tc[4]
            -4.46682914e-01 ;
        /*species 2: AR */
        species[2] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +4.36600000e+00 ;
        /*species 3: N2 */
        species[3] =
            +2.92664000e+00 * tc[0]
            +1.48797680e-03 * tc[1]
            -2.84238000e-07 * tc[2]
            +3.36567933e-11 * tc[3]
            -1.68833775e-15 * tc[4]
            +5.98052800e+00 ;
        /*species 4: HE */
        species[4] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +9.28723974e-01 ;
        /*species 5: O */
        species[5] =
            +2.56942078e+00 * tc[0]
            -8.59741137e-05 * tc[1]
            +2.09742295e-08 * tc[2]
            -3.33925997e-12 * tc[3]
            +3.07084227e-16 * tc[4]
            +4.78433864e+00 ;
        /*species 6: OH */
        species[6] =
            +2.86472886e+00 * tc[0]
            +1.05650448e-03 * tc[1]
            -1.29541379e-07 * tc[2]
            +1.01739558e-11 * tc[3]
            -3.32989690e-16 * tc[4]
            +5.70164073e+00 ;
        /*species 7: HCO */
        species[7] =
            +2.77217438e+00 * tc[0]
            +4.95695526e-03 * tc[1]
            -1.24222806e-06 * tc[2]
            +1.96387259e-10 * tc[3]
            -1.33377178e-14 * tc[4]
            +9.79834492e+00 ;
        /*species 8: HO2 */
        species[8] =
            +4.01721090e+00 * tc[0]
            +2.23982013e-03 * tc[1]
            -3.16829075e-07 * tc[2]
            +3.80821233e-11 * tc[3]
            -2.69771337e-15 * tc[4]
            +3.78510215e+00 ;
        /*species 9: H2O */
        species[9] =
            +3.03399249e+00 * tc[0]
            +2.17691804e-03 * tc[1]
            -8.20362590e-08 * tc[2]
            -3.23473290e-11 * tc[3]
            +4.20502480e-15 * tc[4]
            +4.96677010e+00 ;
        /*species 10: CO */
        species[10] =
            +2.71518561e+00 * tc[0]
            +2.06252743e-03 * tc[1]
            -4.99412886e-07 * tc[2]
            +7.66843360e-11 * tc[3]
            -5.09119290e-15 * tc[4]
            +7.81868772e+00 ;
        /*species 11: O2 */
        species[11] =
            +3.28253784e+00 * tc[0]
            +1.48308754e-03 * tc[1]
            -3.78983334e-07 * tc[2]
            +6.98235183e-11 * tc[3]
            -5.41794485e-15 * tc[4]
            +5.45323129e+00 ;
        /*species 12: H2O2 */
        species[12] =
            +4.16500285e+00 * tc[0]
            +4.90831694e-03 * tc[1]
            -9.50696125e-07 * tc[2]
            +1.23728662e-10 * tc[3]
            -7.19770763e-15 * tc[4]
            +2.91615662e+00 ;
        /*species 13: CO2 */
        species[13] =
            +3.85746029e+00 * tc[0]
            +4.41437026e-03 * tc[1]
            -1.10740702e-06 * tc[2]
            +1.74496729e-10 * tc[3]
            -1.18021041e-14 * tc[4]
            +2.27163806e+00 ;
    }
    return;
}


/*save atomic weights into array */
void atomicWeight(double *  awt)
{
    awt[0] = 15.999400; /*O */
    awt[1] = 1.007970; /*H */
    awt[2] = 12.011150; /*C */
    awt[3] = 14.006700; /*N */
    awt[4] = 39.948000; /*AR */
    awt[5] = 4.002600; /*HE */

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

    double   EPS[14];
    double   SIG[14];
    double    wt[14];
    double avogadro = 6.02214199e23;
    double boltzmann = 1.3806503e-16; //we work in CGS
    double Rcst = 83.144598; //in bar [CGS] !

    egtransetEPS(EPS);
    egtransetSIG(SIG);
    get_mw(wt);

    /*species 0: H2 */
    /*Imported from NIST */
    Tci[0] = 33.145000 ; 
    ai[0] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[0],2.0) / (pow(2.015880,2.0) * 12.964000); 
    bi[0] = 0.08664 * Rcst * Tci[0] / (2.015880 * 12.964000); 
    acentric_i[0] = -0.219000 ;

    /*species 1: H */
    Tci[1] = 1.316 * EPS[1] ; 
    ai[1] = (5.55 * pow(avogadro,2.0) * EPS[1]*boltzmann * pow(1e-8*SIG[1],3.0) ) / (pow(wt[1],2.0)); 
    bi[1] = 0.855 * avogadro * pow(1e-8*SIG[1],3.0) / (wt[1]); 
    acentric_i[1] = 0.0 ;

    /*species 2: AR */
    /*Imported from NIST */
    Tci[2] = 150.860000 ; 
    ai[2] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[2],2.0) / (pow(39.948000,2.0) * 48.980000); 
    bi[2] = 0.08664 * Rcst * Tci[2] / (39.948000 * 48.980000); 
    acentric_i[2] = -0.002000 ;

    /*species 3: N2 */
    /*Imported from NIST */
    Tci[3] = 126.192000 ; 
    ai[3] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[3],2.0) / (pow(28.013400,2.0) * 33.958000); 
    bi[3] = 0.08664 * Rcst * Tci[3] / (28.013400 * 33.958000); 
    acentric_i[3] = 0.037200 ;

    /*species 4: HE */
    Tci[4] = 1.316 * EPS[4] ; 
    ai[4] = (5.55 * pow(avogadro,2.0) * EPS[4]*boltzmann * pow(1e-8*SIG[4],3.0) ) / (pow(wt[4],2.0)); 
    bi[4] = 0.855 * avogadro * pow(1e-8*SIG[4],3.0) / (wt[4]); 
    acentric_i[4] = 0.0 ;

    /*species 5: O */
    Tci[5] = 1.316 * EPS[5] ; 
    ai[5] = (5.55 * pow(avogadro,2.0) * EPS[5]*boltzmann * pow(1e-8*SIG[5],3.0) ) / (pow(wt[5],2.0)); 
    bi[5] = 0.855 * avogadro * pow(1e-8*SIG[5],3.0) / (wt[5]); 
    acentric_i[5] = 0.0 ;

    /*species 6: OH */
    Tci[6] = 1.316 * EPS[6] ; 
    ai[6] = (5.55 * pow(avogadro,2.0) * EPS[6]*boltzmann * pow(1e-8*SIG[6],3.0) ) / (pow(wt[6],2.0)); 
    bi[6] = 0.855 * avogadro * pow(1e-8*SIG[6],3.0) / (wt[6]); 
    acentric_i[6] = 0.0 ;

    /*species 7: HCO */
    Tci[7] = 1.316 * EPS[7] ; 
    ai[7] = (5.55 * pow(avogadro,2.0) * EPS[7]*boltzmann * pow(1e-8*SIG[7],3.0) ) / (pow(wt[7],2.0)); 
    bi[7] = 0.855 * avogadro * pow(1e-8*SIG[7],3.0) / (wt[7]); 
    acentric_i[7] = 0.0 ;

    /*species 8: HO2 */
    Tci[8] = 1.316 * EPS[8] ; 
    ai[8] = (5.55 * pow(avogadro,2.0) * EPS[8]*boltzmann * pow(1e-8*SIG[8],3.0) ) / (pow(wt[8],2.0)); 
    bi[8] = 0.855 * avogadro * pow(1e-8*SIG[8],3.0) / (wt[8]); 
    acentric_i[8] = 0.0 ;

    /*species 9: H2O */
    /*Imported from NIST */
    Tci[9] = 647.096000 ; 
    ai[9] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[9],2.0) / (pow(18.015340,2.0) * 220.640000); 
    bi[9] = 0.08664 * Rcst * Tci[9] / (18.015340 * 220.640000); 
    acentric_i[9] = 0.344300 ;

    /*species 10: CO */
    /*Imported from NIST */
    Tci[10] = 132.850000 ; 
    ai[10] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[10],2.0) / (pow(28.010000,2.0) * 34.940000); 
    bi[10] = 0.08664 * Rcst * Tci[10] / (28.010000 * 34.940000); 
    acentric_i[10] = 0.045000 ;

    /*species 11: O2 */
    /*Imported from NIST */
    Tci[11] = 154.581000 ; 
    ai[11] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[11],2.0) / (pow(31.998800,2.0) * 50.430466); 
    bi[11] = 0.08664 * Rcst * Tci[11] / (31.998800 * 50.430466); 
    acentric_i[11] = 0.022200 ;

    /*species 12: H2O2 */
    Tci[12] = 1.316 * EPS[12] ; 
    ai[12] = (5.55 * pow(avogadro,2.0) * EPS[12]*boltzmann * pow(1e-8*SIG[12],3.0) ) / (pow(wt[12],2.0)); 
    bi[12] = 0.855 * avogadro * pow(1e-8*SIG[12],3.0) / (wt[12]); 
    acentric_i[12] = 0.0 ;

    /*species 13: CO2 */
    /*Imported from NIST */
    Tci[13] = 304.120000 ; 
    ai[13] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[13],2.0) / (pow(44.009950,2.0) * 73.740000); 
    bi[13] = 0.08664 * Rcst * Tci[13] / (44.009950 * 73.740000); 
    acentric_i[13] = 0.225000 ;

    return;
}


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
void egtransetPATM(double* PATM) {
    *PATM =   0.1013250000000000E+07;}


/*the molecular weights in g/mol */
void egtransetWT(double* WT ) {
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
void egtransetEPS(double* EPS ) {
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
void egtransetSIG(double* SIG ) {
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
void egtransetDIP(double* DIP ) {
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
void egtransetPOL(double* POL ) {
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
void egtransetZROT(double* ZROT ) {
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
void egtransetCOFETA(double* COFETA) {
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
void egtransetCOFLAM(double* COFLAM) {
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
void egtransetCOFD(double* COFD) {
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
    KTDIF[0] = 1;
    KTDIF[1] = 2;
    KTDIF[2] = 5;
}


/*Poly fits for thermal diff ratios, dim NO*NLITE*KK */
void egtransetCOFTD(double* COFTD) {
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

/* Replace this routine with the one generated by the Gauss Jordan solver of DW */
AMREX_GPU_HOST_DEVICE void sgjsolve(double* A, double* x, double* b) {
    amrex::Abort("sgjsolve not implemented, choose a different solver ");
}

/* Replace this routine with the one generated by the Gauss Jordan solver of DW */
AMREX_GPU_HOST_DEVICE void sgjsolve_simplified(double* A, double* x, double* b) {
    amrex::Abort("sgjsolve_simplified not implemented, choose a different solver ");
}

