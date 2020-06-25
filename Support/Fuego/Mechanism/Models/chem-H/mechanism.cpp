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
static AMREX_GPU_DEVICE_MANAGED double imw[9] = {
    1.0 / 2.015940,  /*H2 */
    1.0 / 1.007970,  /*H */
    1.0 / 15.999400,  /*O */
    1.0 / 31.998800,  /*O2 */
    1.0 / 17.007370,  /*OH */
    1.0 / 18.015340,  /*H2O */
    1.0 / 33.006770,  /*HO2 */
    1.0 / 34.014740,  /*H2O2 */
    1.0 / 28.013400};  /*N2 */

/* Inverse molecular weights */
/* TODO: check necessity because redundant with molecularWeight */
static AMREX_GPU_DEVICE_MANAGED double molecular_weights[9] = {
    2.015940,  /*H2 */
    1.007970,  /*H */
    15.999400,  /*O */
    31.998800,  /*O2 */
    17.007370,  /*OH */
    18.015340,  /*H2O */
    33.006770,  /*HO2 */
    34.014740,  /*H2O2 */
    28.013400};  /*N2 */

AMREX_GPU_HOST_DEVICE
void get_imw(double imw_new[]){
    for(int i = 0; i<9; ++i) imw_new[i] = imw[i];
}

/* TODO: check necessity because redundant with CKWT */
AMREX_GPU_HOST_DEVICE
void get_mw(double mw_new[]){
    for(int i = 0; i<9; ++i) mw_new[i] = molecular_weights[i];
}


#ifndef AMREX_USE_CUDA
/* Initializes parameter database */
void CKINIT()
{

    rxn_map = {1,2,6,7,8,3,9,10,11,12,4,13,14,5,15,16,17,18,19,20,0,21,22,23,24,25,26};

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
    TB[1] = (double *) malloc(2 * sizeof(double));
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
    TB[2] = (double *) malloc(2 * sizeof(double));
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
    TB[3] = (double *) malloc(3 * sizeof(double));
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
    TB[4] = (double *) malloc(2 * sizeof(double));
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
    TB[5] = (double *) malloc(2 * sizeof(double));
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
    TB[0] = (double *) malloc(2 * sizeof(double));
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
    if (species_id<0 || species_id>=9) {
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
    *mm = 3;
    *kk = 9;
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
    ename.push_back("O");
    ename.push_back("H");
    ename.push_back("N");
}


/* Returns the char strings of element names */
void CKSYME(int * kname, int * plenkname )
{
    int i; /*Loop Counter */
    int lenkname = *plenkname;
    /*clear kname */
    for (i=0; i<lenkname*3; i++) {
        kname[i] = ' ';
    }

    /* O  */
    kname[ 0*lenkname + 0 ] = 'O';
    kname[ 0*lenkname + 1 ] = ' ';

    /* H  */
    kname[ 1*lenkname + 0 ] = 'H';
    kname[ 1*lenkname + 1 ] = ' ';

    /* N  */
    kname[ 2*lenkname + 0 ] = 'N';
    kname[ 2*lenkname + 1 ] = ' ';

}


/* Returns the vector of strings of species names */
void CKSYMS_STR(amrex::Vector<std::string>& kname)
{
    kname.push_back("H2");
    kname.push_back("H");
    kname.push_back("O");
    kname.push_back("O2");
    kname.push_back("OH");
    kname.push_back("H2O");
    kname.push_back("HO2");
    kname.push_back("H2O2");
    kname.push_back("N2");
}


/* Returns the char strings of species names */
void CKSYMS(int * kname, int * plenkname )
{
    int i; /*Loop Counter */
    int lenkname = *plenkname;
    /*clear kname */
    for (i=0; i<lenkname*9; i++) {
        kname[i] = ' ';
    }

    /* H2  */
    kname[ 0*lenkname + 0 ] = 'H';
    kname[ 0*lenkname + 1 ] = '2';
    kname[ 0*lenkname + 2 ] = ' ';

    /* H  */
    kname[ 1*lenkname + 0 ] = 'H';
    kname[ 1*lenkname + 1 ] = ' ';

    /* O  */
    kname[ 2*lenkname + 0 ] = 'O';
    kname[ 2*lenkname + 1 ] = ' ';

    /* O2  */
    kname[ 3*lenkname + 0 ] = 'O';
    kname[ 3*lenkname + 1 ] = '2';
    kname[ 3*lenkname + 2 ] = ' ';

    /* OH  */
    kname[ 4*lenkname + 0 ] = 'O';
    kname[ 4*lenkname + 1 ] = 'H';
    kname[ 4*lenkname + 2 ] = ' ';

    /* H2O  */
    kname[ 5*lenkname + 0 ] = 'H';
    kname[ 5*lenkname + 1 ] = '2';
    kname[ 5*lenkname + 2 ] = 'O';
    kname[ 5*lenkname + 3 ] = ' ';

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
    XW += x[2]*15.999400; /*O */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*18.015340; /*H2O */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    *P = *rho * 8.31446e+07 * (*T) / XW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(y) */
AMREX_GPU_HOST_DEVICE void CKPY(double *  rho, double *  T, double *  y,  double *  P)
{
    double YOW = 0;/* for computing mean MW */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*H */
    YOW += y[2]*imw[2]; /*O */
    YOW += y[3]*imw[3]; /*O2 */
    YOW += y[4]*imw[4]; /*OH */
    YOW += y[5]*imw[5]; /*H2O */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*N2 */
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

    for (int n=0; n<9; n++) {
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
    W += c[0]*2.015940; /*H2 */
    W += c[1]*1.007970; /*H */
    W += c[2]*15.999400; /*O */
    W += c[3]*31.998800; /*O2 */
    W += c[4]*17.007370; /*OH */
    W += c[5]*18.015340; /*H2O */
    W += c[6]*33.006770; /*HO2 */
    W += c[7]*34.014740; /*H2O2 */
    W += c[8]*28.013400; /*N2 */

    for (id = 0; id < 9; ++id) {
        sumC += c[id];
    }
    *P = *rho * 8.31446e+07 * (*T) * sumC / W; /*P = rho*R*T/W */

    return;
}


/*Compute rho = PW(x)/RT */
void CKRHOX(double *  P, double *  T, double *  x,  double *  rho)
{
    double XW = 0;/* To hold mean molecular wt */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*1.007970; /*H */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*18.015340; /*H2O */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    *rho = *P * XW / (8.31446e+07 * (*T)); /*rho = P*W/(R*T) */

    return;
}


/*Compute rho = P*W(y)/RT */
AMREX_GPU_HOST_DEVICE void CKRHOY(double *  P, double *  T, double *  y,  double *  rho)
{
    double YOW = 0;
    double tmp[9];

    for (int i = 0; i < 9; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 9; i++)
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
    W += c[0]*2.015940; /*H2 */
    W += c[1]*1.007970; /*H */
    W += c[2]*15.999400; /*O */
    W += c[3]*31.998800; /*O2 */
    W += c[4]*17.007370; /*OH */
    W += c[5]*18.015340; /*H2O */
    W += c[6]*33.006770; /*HO2 */
    W += c[7]*34.014740; /*H2O2 */
    W += c[8]*28.013400; /*N2 */

    for (id = 0; id < 9; ++id) {
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
    double tmp[9];

    for (int i = 0; i < 9; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 9; i++)
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
    XW += x[2]*15.999400; /*O */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*18.015340; /*H2O */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
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
    W += c[2]*15.999400; /*O */
    W += c[3]*31.998800; /*O2 */
    W += c[4]*17.007370; /*OH */
    W += c[5]*18.015340; /*H2O */
    W += c[6]*33.006770; /*HO2 */
    W += c[7]*34.014740; /*H2O2 */
    W += c[8]*28.013400; /*N2 */

    for (id = 0; id < 9; ++id) {
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
    double tmp[9];

    for (int i = 0; i < 9; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 9; i++)
    {
        YOW += tmp[i];
    }

    double YOWINV = 1.0/YOW;

    for (int i = 0; i < 9; i++)
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
    for (int i = 0; i < 9; i++)
    {
        c[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 9; i++)
    {
        YOW += c[i];
    }

    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31446e+07 * (*T)); 
    /*Now compute conversion */

    for (int i = 0; i < 9; i++)
    {
        c[i] = PWORT * y[i] * imw[i];
    }
    return;
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
AMREX_GPU_HOST_DEVICE void CKYTCR(double *  rho, double *  T, double *  y,  double *  c)
{
    for (int i = 0; i < 9; i++)
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
    XW += x[2]*15.999400; /*O */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*18.015340; /*H2O */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    /*Now compute conversion */
    double XWinv = 1.0/XW;
    y[0] = x[0]*2.015940*XWinv; 
    y[1] = x[1]*1.007970*XWinv; 
    y[2] = x[2]*15.999400*XWinv; 
    y[3] = x[3]*31.998800*XWinv; 
    y[4] = x[4]*17.007370*XWinv; 
    y[5] = x[5]*18.015340*XWinv; 
    y[6] = x[6]*33.006770*XWinv; 
    y[7] = x[7]*34.014740*XWinv; 
    y[8] = x[8]*28.013400*XWinv; 

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
void CKXTCP(double *  P, double *  T, double *  x,  double *  c)
{
    int id; /*loop counter */
    double PORT = (*P)/(8.31446e+07 * (*T)); /*P/RT */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 9; ++id) {
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
    XW += x[2]*15.999400; /*O */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*18.015340; /*H2O */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    ROW = (*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 9; ++id) {
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
    for (id = 0; id < 9; ++id) {
        sumC += c[id];
    }

    /* See Eq 13  */
    double sumCinv = 1.0/sumC;
    for (id = 0; id < 9; ++id) {
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
    CW += c[2]*15.999400; /*O */
    CW += c[3]*31.998800; /*O2 */
    CW += c[4]*17.007370; /*OH */
    CW += c[5]*18.015340; /*H2O */
    CW += c[6]*33.006770; /*HO2 */
    CW += c[7]*34.014740; /*H2O2 */
    CW += c[8]*28.013400; /*N2 */
    /*Now compute conversion */
    double CWinv = 1.0/CW;
    y[0] = c[0]*2.015940*CWinv; 
    y[1] = c[1]*1.007970*CWinv; 
    y[2] = c[2]*15.999400*CWinv; 
    y[3] = c[3]*31.998800*CWinv; 
    y[4] = c[4]*17.007370*CWinv; 
    y[5] = c[5]*18.015340*CWinv; 
    y[6] = c[6]*33.006770*CWinv; 
    y[7] = c[7]*34.014740*CWinv; 
    y[8] = c[8]*28.013400*CWinv; 

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
    for (id = 0; id < 9; ++id) {
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
    for (id = 0; id < 9; ++id) {
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
    for (id = 0; id < 9; ++id) {
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
    for (id = 0; id < 9; ++id) {
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
    for (id = 0; id < 9; ++id) {
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
    for (id = 0; id < 9; ++id) {
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
    for (id = 0; id < 9; ++id) {
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
    cvms[0] *= 4.124360158612479e+07; /*H2 */
    cvms[1] *= 8.248720317224957e+07; /*H */
    cvms[2] *= 5.196734013871295e+06; /*O */
    cvms[3] *= 2.598367006935648e+06; /*O2 */
    cvms[4] *= 4.888740950630956e+06; /*OH */
    cvms[5] *= 4.615212712140454e+06; /*H2O */
    cvms[6] *= 2.519017346487778e+06; /*HO2 */
    cvms[7] *= 2.444370475315478e+06; /*H2O2 */
    cvms[8] *= 2.968030520448514e+06; /*N2 */
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
    cpms[2] *= 5.196734013871295e+06; /*O */
    cpms[3] *= 2.598367006935648e+06; /*O2 */
    cpms[4] *= 4.888740950630956e+06; /*OH */
    cpms[5] *= 4.615212712140454e+06; /*H2O */
    cpms[6] *= 2.519017346487778e+06; /*HO2 */
    cpms[7] *= 2.444370475315478e+06; /*H2O2 */
    cpms[8] *= 2.968030520448514e+06; /*N2 */
}


/*Returns internal energy in mass units (Eq 30.) */
AMREX_GPU_HOST_DEVICE void CKUMS(double *  T,  double *  ums)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double RT = 8.31446e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    for (int i = 0; i < 9; i++)
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
    for (int i = 0; i < 9; i++)
    {
        hms[i] *= RT*imw[i];
    }
}


#ifndef AMREX_USE_CUDA
/*Returns enthalpy in mass units (Eq 27.) */
void VCKHMS(int *  np, double *  T,  double *  hms)
{
    double tc[5], h[9];

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
    for (int i = 0; i < 9; i++)
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
    for (int i = 0; i < 9; i++)
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
    sms[2] *= 5.196734013871295e+06; /*O */
    sms[3] *= 2.598367006935648e+06; /*O2 */
    sms[4] *= 4.888740950630956e+06; /*OH */
    sms[5] *= 4.615212712140454e+06; /*H2O */
    sms[6] *= 2.519017346487778e+06; /*HO2 */
    sms[7] *= 2.444370475315478e+06; /*H2O2 */
    sms[8] *= 2.968030520448514e+06; /*N2 */
}


/*Returns the mean specific heat at CP (Eq. 33) */
void CKCPBL(double *  T, double *  x,  double *  cpbl)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double cpor[9]; /* temporary storage */
    cp_R(cpor, tc);

    /*perform dot product */
    for (id = 0; id < 9; ++id) {
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
    double cpor[9], tresult[9]; /* temporary storage */
    cp_R(cpor, tc);
    for (int i = 0; i < 9; i++)
    {
        tresult[i] = cpor[i]*y[i]*imw[i];

    }
    for (int i = 0; i < 9; i++)
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
    double cvor[9]; /* temporary storage */
    cv_R(cvor, tc);

    /*perform dot product */
    for (id = 0; id < 9; ++id) {
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
    double cvor[9]; /* temporary storage */
    cv_R(cvor, tc);
    /*multiply by y/molecularweight */
    result += cvor[0]*y[0]*imw[0]; /*H2 */
    result += cvor[1]*y[1]*imw[1]; /*H */
    result += cvor[2]*y[2]*imw[2]; /*O */
    result += cvor[3]*y[3]*imw[3]; /*O2 */
    result += cvor[4]*y[4]*imw[4]; /*OH */
    result += cvor[5]*y[5]*imw[5]; /*H2O */
    result += cvor[6]*y[6]*imw[6]; /*HO2 */
    result += cvor[7]*y[7]*imw[7]; /*H2O2 */
    result += cvor[8]*y[8]*imw[8]; /*N2 */

    *cvbs = result * 8.31446e+07;
}


/*Returns the mean enthalpy of the mixture in molar units */
void CKHBML(double *  T, double *  x,  double *  hbml)
{
    int id; /*loop counter */
    double result = 0; 
    double tT = *T; /*temporary temperature */
    double tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double hml[9]; /* temporary storage */
    double RT = 8.31446e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*perform dot product */
    for (id = 0; id < 9; ++id) {
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
    double hml[9], tmp[9]; /* temporary storage */
    double RT = 8.31446e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);
    int id;
    for (id = 0; id < 9; ++id) {
        tmp[id] = y[id]*hml[id]*imw[id];
    }
    for (id = 0; id < 9; ++id) {
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
    double uml[9]; /* temporary energy array */
    double RT = 8.31446e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*perform dot product */
    for (id = 0; id < 9; ++id) {
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
    double ums[9]; /* temporary energy array */
    double RT = 8.31446e+07*tT; /*R*T */
    speciesInternalEnergy(ums, tc);
    /*perform dot product + scaling by wt */
    result += y[0]*ums[0]*imw[0]; /*H2 */
    result += y[1]*ums[1]*imw[1]; /*H */
    result += y[2]*ums[2]*imw[2]; /*O */
    result += y[3]*ums[3]*imw[3]; /*O2 */
    result += y[4]*ums[4]*imw[4]; /*OH */
    result += y[5]*ums[5]*imw[5]; /*H2O */
    result += y[6]*ums[6]*imw[6]; /*HO2 */
    result += y[7]*ums[7]*imw[7]; /*H2O2 */
    result += y[8]*ums[8]*imw[8]; /*N2 */

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
    double sor[9]; /* temporary storage */
    speciesEntropy(sor, tc);

    /*Compute Eq 42 */
    for (id = 0; id < 9; ++id) {
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
    double sor[9]; /* temporary storage */
    double x[9]; /* need a ytx conversion */
    double YOW = 0; /*See Eq 4, 6 in CK Manual */
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
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(1.007970*YOW); 
    x[2] = y[2]/(15.999400*YOW); 
    x[3] = y[3]/(31.998800*YOW); 
    x[4] = y[4]/(17.007370*YOW); 
    x[5] = y[5]/(18.015340*YOW); 
    x[6] = y[6]/(33.006770*YOW); 
    x[7] = y[7]/(34.014740*YOW); 
    x[8] = y[8]/(28.013400*YOW); 
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
    double gort[9]; /* temporary storage */
    /*Compute g/RT */
    gibbs(gort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 9; ++id) {
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
    double gort[9]; /* temporary storage */
    double x[9]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
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
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(1.007970*YOW); 
    x[2] = y[2]/(15.999400*YOW); 
    x[3] = y[3]/(31.998800*YOW); 
    x[4] = y[4]/(17.007370*YOW); 
    x[5] = y[5]/(18.015340*YOW); 
    x[6] = y[6]/(33.006770*YOW); 
    x[7] = y[7]/(34.014740*YOW); 
    x[8] = y[8]/(28.013400*YOW); 
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
    double aort[9]; /* temporary storage */
    /*Compute g/RT */
    helmholtz(aort, tc);

    /*Compute Eq 44 */
    for (id = 0; id < 9; ++id) {
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
    double aort[9]; /* temporary storage */
    double x[9]; /* need a ytx conversion */
    double YOW = 0; /*To hold 1/molecularweight */
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
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(1.007970*YOW); 
    x[2] = y[2]/(15.999400*YOW); 
    x[3] = y[3]/(31.998800*YOW); 
    x[4] = y[4]/(17.007370*YOW); 
    x[5] = y[5]/(18.015340*YOW); 
    x[6] = y[6]/(33.006770*YOW); 
    x[7] = y[7]/(34.014740*YOW); 
    x[8] = y[8]/(28.013400*YOW); 
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
    /*Scale by RT/W */
    *abms = result * RT * YOW;
}


/*compute the production rate for each species */
AMREX_GPU_HOST_DEVICE void CKWC(double *  T, double *  C,  double *  wdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 9; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    productionRate(wdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        C[id] *= 1.0e-6;
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mass fractions */
void CKWYP(double *  P, double *  T, double *  y,  double *  wdot)
{
    int id; /*loop counter */
    double c[9]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
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

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mole fractions */
void CKWXP(double *  P, double *  T, double *  x,  double *  wdot)
{
    int id; /*loop counter */
    double c[9]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31446e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 9; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
AMREX_GPU_HOST_DEVICE void CKWYR(double *  rho, double *  T, double *  y,  double *  wdot)
{
    int id; /*loop counter */
    double c[9]; /*temporary storage */
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

    /*call productionRate */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
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
    double c[9*(*np)]; /*temporary storage */
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
#endif
}


/*Returns the molar production rate of species */
/*Given rho, T, and mole fractions */
void CKWXR(double *  rho, double *  T, double *  x,  double *  wdot)
{
    int id; /*loop counter */
    double c[9]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
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
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the rate of progress for each reaction */
void CKQC(double *  T, double *  C, double *  qdot)
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
/*Given P, T, and mole fractions */
void CKKFKR(double *  P, double *  T, double *  x, double *  q_f, double *  q_r)
{
    int id; /*loop counter */
    double c[9]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31446e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

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


/*Returns the progress rates of each reactions */
/*Given P, T, and mass fractions */
void CKQYP(double *  P, double *  T, double *  y, double *  qdot)
{
    int id; /*loop counter */
    double c[9]; /*temporary storage */
    double YOW = 0; 
    double PWORT; 
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
    double c[9]; /*temporary storage */
    double PORT = 1e6 * (*P)/(8.31446e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

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
void CKQYR(double *  rho, double *  T, double *  y, double *  qdot)
{
    int id; /*loop counter */
    double c[9]; /*temporary storage */
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
void CKQXR(double *  rho, double *  T, double *  x, double *  qdot)
{
    int id; /*loop counter */
    double c[9]; /*temporary storage */
    double XW = 0; /*See Eq 4, 11 in CK Manual */
    double ROW; 
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


/*Returns the stoichiometric coefficients */
/*of the reaction mechanism. (Eq 50) */
void CKNU(int * kdim,  int * nuki)
{
    int id; /*loop counter */
    int kd = (*kdim); 
    /*Zero nuki */
    for (id = 0; id < 9 * kd; ++ id) {
         nuki[id] = 0; 
    }

    /*reaction 1: 2.000000 OH (+M) <=> H2O2 (+M) */
    nuki[ 4 * kd + 0 ] += -2.000000 ;
    nuki[ 7 * kd + 0 ] += +1.000000 ;

    /*reaction 2: 2.000000 O + M <=> O2 + M */
    nuki[ 2 * kd + 1 ] += -2.000000 ;
    nuki[ 3 * kd + 1 ] += +1.000000 ;

    /*reaction 3: O + H + M <=> OH + M */
    nuki[ 2 * kd + 2 ] += -1.000000 ;
    nuki[ 1 * kd + 2 ] += -1.000000 ;
    nuki[ 4 * kd + 2 ] += +1.000000 ;

    /*reaction 4: H + O2 + M <=> HO2 + M */
    nuki[ 1 * kd + 3 ] += -1.000000 ;
    nuki[ 3 * kd + 3 ] += -1.000000 ;
    nuki[ 6 * kd + 3 ] += +1.000000 ;

    /*reaction 5: 2.000000 H + M <=> H2 + M */
    nuki[ 1 * kd + 4 ] += -2.000000 ;
    nuki[ 0 * kd + 4 ] += +1.000000 ;

    /*reaction 6: H + OH + M <=> H2O + M */
    nuki[ 1 * kd + 5 ] += -1.000000 ;
    nuki[ 4 * kd + 5 ] += -1.000000 ;
    nuki[ 5 * kd + 5 ] += +1.000000 ;

    /*reaction 7: O + H2 <=> H + OH */
    nuki[ 2 * kd + 6 ] += -1.000000 ;
    nuki[ 0 * kd + 6 ] += -1.000000 ;
    nuki[ 1 * kd + 6 ] += +1.000000 ;
    nuki[ 4 * kd + 6 ] += +1.000000 ;

    /*reaction 8: O + HO2 <=> OH + O2 */
    nuki[ 2 * kd + 7 ] += -1.000000 ;
    nuki[ 6 * kd + 7 ] += -1.000000 ;
    nuki[ 4 * kd + 7 ] += +1.000000 ;
    nuki[ 3 * kd + 7 ] += +1.000000 ;

    /*reaction 9: O + H2O2 <=> OH + HO2 */
    nuki[ 2 * kd + 8 ] += -1.000000 ;
    nuki[ 7 * kd + 8 ] += -1.000000 ;
    nuki[ 4 * kd + 8 ] += +1.000000 ;
    nuki[ 6 * kd + 8 ] += +1.000000 ;

    /*reaction 10: H + 2.000000 O2 <=> HO2 + O2 */
    nuki[ 1 * kd + 9 ] += -1.000000 ;
    nuki[ 3 * kd + 9 ] += -2.000000 ;
    nuki[ 6 * kd + 9 ] += +1.000000 ;
    nuki[ 3 * kd + 9 ] += +1.000000 ;

    /*reaction 11: H + O2 + H2O <=> HO2 + H2O */
    nuki[ 1 * kd + 10 ] += -1.000000 ;
    nuki[ 3 * kd + 10 ] += -1.000000 ;
    nuki[ 5 * kd + 10 ] += -1.000000 ;
    nuki[ 6 * kd + 10 ] += +1.000000 ;
    nuki[ 5 * kd + 10 ] += +1.000000 ;

    /*reaction 12: H + O2 + N2 <=> HO2 + N2 */
    nuki[ 1 * kd + 11 ] += -1.000000 ;
    nuki[ 3 * kd + 11 ] += -1.000000 ;
    nuki[ 8 * kd + 11 ] += -1.000000 ;
    nuki[ 6 * kd + 11 ] += +1.000000 ;
    nuki[ 8 * kd + 11 ] += +1.000000 ;

    /*reaction 13: H + O2 <=> O + OH */
    nuki[ 1 * kd + 12 ] += -1.000000 ;
    nuki[ 3 * kd + 12 ] += -1.000000 ;
    nuki[ 2 * kd + 12 ] += +1.000000 ;
    nuki[ 4 * kd + 12 ] += +1.000000 ;

    /*reaction 14: 2.000000 H + H2 <=> 2.000000 H2 */
    nuki[ 1 * kd + 13 ] += -2.000000 ;
    nuki[ 0 * kd + 13 ] += -1.000000 ;
    nuki[ 0 * kd + 13 ] += +2.000000 ;

    /*reaction 15: 2.000000 H + H2O <=> H2 + H2O */
    nuki[ 1 * kd + 14 ] += -2.000000 ;
    nuki[ 5 * kd + 14 ] += -1.000000 ;
    nuki[ 0 * kd + 14 ] += +1.000000 ;
    nuki[ 5 * kd + 14 ] += +1.000000 ;

    /*reaction 16: H + HO2 <=> O + H2O */
    nuki[ 1 * kd + 15 ] += -1.000000 ;
    nuki[ 6 * kd + 15 ] += -1.000000 ;
    nuki[ 2 * kd + 15 ] += +1.000000 ;
    nuki[ 5 * kd + 15 ] += +1.000000 ;

    /*reaction 17: H + HO2 <=> O2 + H2 */
    nuki[ 1 * kd + 16 ] += -1.000000 ;
    nuki[ 6 * kd + 16 ] += -1.000000 ;
    nuki[ 3 * kd + 16 ] += +1.000000 ;
    nuki[ 0 * kd + 16 ] += +1.000000 ;

    /*reaction 18: H + HO2 <=> 2.000000 OH */
    nuki[ 1 * kd + 17 ] += -1.000000 ;
    nuki[ 6 * kd + 17 ] += -1.000000 ;
    nuki[ 4 * kd + 17 ] += +2.000000 ;

    /*reaction 19: H + H2O2 <=> HO2 + H2 */
    nuki[ 1 * kd + 18 ] += -1.000000 ;
    nuki[ 7 * kd + 18 ] += -1.000000 ;
    nuki[ 6 * kd + 18 ] += +1.000000 ;
    nuki[ 0 * kd + 18 ] += +1.000000 ;

    /*reaction 20: H + H2O2 <=> OH + H2O */
    nuki[ 1 * kd + 19 ] += -1.000000 ;
    nuki[ 7 * kd + 19 ] += -1.000000 ;
    nuki[ 4 * kd + 19 ] += +1.000000 ;
    nuki[ 5 * kd + 19 ] += +1.000000 ;

    /*reaction 21: OH + H2 <=> H + H2O */
    nuki[ 4 * kd + 20 ] += -1.000000 ;
    nuki[ 0 * kd + 20 ] += -1.000000 ;
    nuki[ 1 * kd + 20 ] += +1.000000 ;
    nuki[ 5 * kd + 20 ] += +1.000000 ;

    /*reaction 22: 2.000000 OH <=> O + H2O */
    nuki[ 4 * kd + 21 ] += -2.000000 ;
    nuki[ 2 * kd + 21 ] += +1.000000 ;
    nuki[ 5 * kd + 21 ] += +1.000000 ;

    /*reaction 23: OH + HO2 <=> O2 + H2O */
    nuki[ 4 * kd + 22 ] += -1.000000 ;
    nuki[ 6 * kd + 22 ] += -1.000000 ;
    nuki[ 3 * kd + 22 ] += +1.000000 ;
    nuki[ 5 * kd + 22 ] += +1.000000 ;

    /*reaction 24: OH + H2O2 <=> HO2 + H2O */
    nuki[ 4 * kd + 23 ] += -1.000000 ;
    nuki[ 7 * kd + 23 ] += -1.000000 ;
    nuki[ 6 * kd + 23 ] += +1.000000 ;
    nuki[ 5 * kd + 23 ] += +1.000000 ;

    /*reaction 25: OH + H2O2 <=> HO2 + H2O */
    nuki[ 4 * kd + 24 ] += -1.000000 ;
    nuki[ 7 * kd + 24 ] += -1.000000 ;
    nuki[ 6 * kd + 24 ] += +1.000000 ;
    nuki[ 5 * kd + 24 ] += +1.000000 ;

    /*reaction 26: 2.000000 HO2 <=> O2 + H2O2 */
    nuki[ 6 * kd + 25 ] += -2.000000 ;
    nuki[ 3 * kd + 25 ] += +1.000000 ;
    nuki[ 7 * kd + 25 ] += +1.000000 ;

    /*reaction 27: 2.000000 HO2 <=> O2 + H2O2 */
    nuki[ 6 * kd + 26 ] += -2.000000 ;
    nuki[ 3 * kd + 26 ] += +1.000000 ;
    nuki[ 7 * kd + 26 ] += +1.000000 ;
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


/*Returns the arrehenius coefficients  */
/*for all reactions */
void CKABE( double *  a, double *  b, double *  e)
{
    // (20):  2.000000 OH (+M) <=> H2O2 (+M)
    a[0] = 74000000000000;
    b[0] = -0.37;
    e[0] = 0;

    // (0):  2.000000 O + M <=> O2 + M
    a[1] = 1.2e+17;
    b[1] = -1;
    e[1] = 0;

    // (1):  O + H + M <=> OH + M
    a[2] = 5e+17;
    b[2] = -1;
    e[2] = 0;

    // (5):  H + O2 + M <=> HO2 + M
    a[3] = 2.8e+18;
    b[3] = -0.85999999999999999;
    e[3] = 0;

    // (10):  2.000000 H + M <=> H2 + M
    a[4] = 1e+18;
    b[4] = -1;
    e[4] = 0;

    // (13):  H + OH + M <=> H2O + M
    a[5] = 2.2e+22;
    b[5] = -2;
    e[5] = 0;

    // (2):  O + H2 <=> H + OH
    a[6] = 50000;
    b[6] = 2.6699999999999999;
    e[6] = 6290;

    // (3):  O + HO2 <=> OH + O2
    a[7] = 20000000000000;
    b[7] = 0;
    e[7] = 0;

    // (4):  O + H2O2 <=> OH + HO2
    a[8] = 9630000;
    b[8] = 2;
    e[8] = 4000;

    // (6):  H + 2.000000 O2 <=> HO2 + O2
    a[9] = 3e+20;
    b[9] = -1.72;
    e[9] = 0;

    // (7):  H + O2 + H2O <=> HO2 + H2O
    a[10] = 9.38e+18;
    b[10] = -0.76000000000000001;
    e[10] = 0;

    // (8):  H + O2 + N2 <=> HO2 + N2
    a[11] = 3.75e+20;
    b[11] = -1.72;
    e[11] = 0;

    // (9):  H + O2 <=> O + OH
    a[12] = 83000000000000;
    b[12] = 0;
    e[12] = 14413;

    // (11):  2.000000 H + H2 <=> 2.000000 H2
    a[13] = 90000000000000000;
    b[13] = -0.59999999999999998;
    e[13] = 0;

    // (12):  2.000000 H + H2O <=> H2 + H2O
    a[14] = 6e+19;
    b[14] = -1.25;
    e[14] = 0;

    // (14):  H + HO2 <=> O + H2O
    a[15] = 3970000000000;
    b[15] = 0;
    e[15] = 671;

    // (15):  H + HO2 <=> O2 + H2
    a[16] = 28000000000000;
    b[16] = 0;
    e[16] = 1068;

    // (16):  H + HO2 <=> 2.000000 OH
    a[17] = 134000000000000;
    b[17] = 0;
    e[17] = 635;

    // (17):  H + H2O2 <=> HO2 + H2
    a[18] = 12100000;
    b[18] = 2;
    e[18] = 5200;

    // (18):  H + H2O2 <=> OH + H2O
    a[19] = 10000000000000;
    b[19] = 0;
    e[19] = 3600;

    // (19):  OH + H2 <=> H + H2O
    a[20] = 216000000;
    b[20] = 1.51;
    e[20] = 3430;

    // (21):  2.000000 OH <=> O + H2O
    a[21] = 35700;
    b[21] = 2.3999999999999999;
    e[21] = -2110;

    // (22):  OH + HO2 <=> O2 + H2O
    a[22] = 29000000000000;
    b[22] = 0;
    e[22] = -500;

    // (23):  OH + H2O2 <=> HO2 + H2O
    a[23] = 1750000000000;
    b[23] = 0;
    e[23] = 320;

    // (24):  OH + H2O2 <=> HO2 + H2O
    a[24] = 580000000000000;
    b[24] = 0;
    e[24] = 9560;

    // (25):  2.000000 HO2 <=> O2 + H2O2
    a[25] = 130000000000;
    b[25] = 0;
    e[25] = -1630;

    // (26):  2.000000 HO2 <=> O2 + H2O2
    a[26] = 420000000000000;
    b[26] = 0;
    e[26] = 12000;


    return;
}


/*Returns the equil constants for each reaction */
void CKEQC(double *  T, double *  C, double *  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[9]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: 2.000000 OH (+M) <=> H2O2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: 2.000000 O + M <=> O2 + M */
    eqcon[1] *= 1e+06; 

    /*reaction 3: O + H + M <=> OH + M */
    eqcon[2] *= 1e+06; 

    /*reaction 4: H + O2 + M <=> HO2 + M */
    eqcon[3] *= 1e+06; 

    /*reaction 5: 2.000000 H + M <=> H2 + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: H + OH + M <=> H2O + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: O + H2 <=> H + OH */
    /*eqcon[6] *= 1;  */

    /*reaction 8: O + HO2 <=> OH + O2 */
    /*eqcon[7] *= 1;  */

    /*reaction 9: O + H2O2 <=> OH + HO2 */
    /*eqcon[8] *= 1;  */

    /*reaction 10: H + 2.000000 O2 <=> HO2 + O2 */
    eqcon[9] *= 1e+06; 

    /*reaction 11: H + O2 + H2O <=> HO2 + H2O */
    eqcon[10] *= 1e+06; 

    /*reaction 12: H + O2 + N2 <=> HO2 + N2 */
    eqcon[11] *= 1e+06; 

    /*reaction 13: H + O2 <=> O + OH */
    /*eqcon[12] *= 1;  */

    /*reaction 14: 2.000000 H + H2 <=> 2.000000 H2 */
    eqcon[13] *= 1e+06; 

    /*reaction 15: 2.000000 H + H2O <=> H2 + H2O */
    eqcon[14] *= 1e+06; 

    /*reaction 16: H + HO2 <=> O + H2O */
    /*eqcon[15] *= 1;  */

    /*reaction 17: H + HO2 <=> O2 + H2 */
    /*eqcon[16] *= 1;  */

    /*reaction 18: H + HO2 <=> 2.000000 OH */
    /*eqcon[17] *= 1;  */

    /*reaction 19: H + H2O2 <=> HO2 + H2 */
    /*eqcon[18] *= 1;  */

    /*reaction 20: H + H2O2 <=> OH + H2O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: OH + H2 <=> H + H2O */
    /*eqcon[20] *= 1;  */

    /*reaction 22: 2.000000 OH <=> O + H2O */
    /*eqcon[21] *= 1;  */

    /*reaction 23: OH + HO2 <=> O2 + H2O */
    /*eqcon[22] *= 1;  */

    /*reaction 24: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[23] *= 1;  */

    /*reaction 25: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[24] *= 1;  */

    /*reaction 26: 2.000000 HO2 <=> O2 + H2O2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: 2.000000 HO2 <=> O2 + H2O2 */
    /*eqcon[26] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mass fractions */
void CKEQYP(double *  P, double *  T, double *  y, double *  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[9]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: 2.000000 OH (+M) <=> H2O2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: 2.000000 O + M <=> O2 + M */
    eqcon[1] *= 1e+06; 

    /*reaction 3: O + H + M <=> OH + M */
    eqcon[2] *= 1e+06; 

    /*reaction 4: H + O2 + M <=> HO2 + M */
    eqcon[3] *= 1e+06; 

    /*reaction 5: 2.000000 H + M <=> H2 + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: H + OH + M <=> H2O + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: O + H2 <=> H + OH */
    /*eqcon[6] *= 1;  */

    /*reaction 8: O + HO2 <=> OH + O2 */
    /*eqcon[7] *= 1;  */

    /*reaction 9: O + H2O2 <=> OH + HO2 */
    /*eqcon[8] *= 1;  */

    /*reaction 10: H + 2.000000 O2 <=> HO2 + O2 */
    eqcon[9] *= 1e+06; 

    /*reaction 11: H + O2 + H2O <=> HO2 + H2O */
    eqcon[10] *= 1e+06; 

    /*reaction 12: H + O2 + N2 <=> HO2 + N2 */
    eqcon[11] *= 1e+06; 

    /*reaction 13: H + O2 <=> O + OH */
    /*eqcon[12] *= 1;  */

    /*reaction 14: 2.000000 H + H2 <=> 2.000000 H2 */
    eqcon[13] *= 1e+06; 

    /*reaction 15: 2.000000 H + H2O <=> H2 + H2O */
    eqcon[14] *= 1e+06; 

    /*reaction 16: H + HO2 <=> O + H2O */
    /*eqcon[15] *= 1;  */

    /*reaction 17: H + HO2 <=> O2 + H2 */
    /*eqcon[16] *= 1;  */

    /*reaction 18: H + HO2 <=> 2.000000 OH */
    /*eqcon[17] *= 1;  */

    /*reaction 19: H + H2O2 <=> HO2 + H2 */
    /*eqcon[18] *= 1;  */

    /*reaction 20: H + H2O2 <=> OH + H2O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: OH + H2 <=> H + H2O */
    /*eqcon[20] *= 1;  */

    /*reaction 22: 2.000000 OH <=> O + H2O */
    /*eqcon[21] *= 1;  */

    /*reaction 23: OH + HO2 <=> O2 + H2O */
    /*eqcon[22] *= 1;  */

    /*reaction 24: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[23] *= 1;  */

    /*reaction 25: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[24] *= 1;  */

    /*reaction 26: 2.000000 HO2 <=> O2 + H2O2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: 2.000000 HO2 <=> O2 + H2O2 */
    /*eqcon[26] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given P, T, and mole fractions */
void CKEQXP(double *  P, double *  T, double *  x, double *  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[9]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: 2.000000 OH (+M) <=> H2O2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: 2.000000 O + M <=> O2 + M */
    eqcon[1] *= 1e+06; 

    /*reaction 3: O + H + M <=> OH + M */
    eqcon[2] *= 1e+06; 

    /*reaction 4: H + O2 + M <=> HO2 + M */
    eqcon[3] *= 1e+06; 

    /*reaction 5: 2.000000 H + M <=> H2 + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: H + OH + M <=> H2O + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: O + H2 <=> H + OH */
    /*eqcon[6] *= 1;  */

    /*reaction 8: O + HO2 <=> OH + O2 */
    /*eqcon[7] *= 1;  */

    /*reaction 9: O + H2O2 <=> OH + HO2 */
    /*eqcon[8] *= 1;  */

    /*reaction 10: H + 2.000000 O2 <=> HO2 + O2 */
    eqcon[9] *= 1e+06; 

    /*reaction 11: H + O2 + H2O <=> HO2 + H2O */
    eqcon[10] *= 1e+06; 

    /*reaction 12: H + O2 + N2 <=> HO2 + N2 */
    eqcon[11] *= 1e+06; 

    /*reaction 13: H + O2 <=> O + OH */
    /*eqcon[12] *= 1;  */

    /*reaction 14: 2.000000 H + H2 <=> 2.000000 H2 */
    eqcon[13] *= 1e+06; 

    /*reaction 15: 2.000000 H + H2O <=> H2 + H2O */
    eqcon[14] *= 1e+06; 

    /*reaction 16: H + HO2 <=> O + H2O */
    /*eqcon[15] *= 1;  */

    /*reaction 17: H + HO2 <=> O2 + H2 */
    /*eqcon[16] *= 1;  */

    /*reaction 18: H + HO2 <=> 2.000000 OH */
    /*eqcon[17] *= 1;  */

    /*reaction 19: H + H2O2 <=> HO2 + H2 */
    /*eqcon[18] *= 1;  */

    /*reaction 20: H + H2O2 <=> OH + H2O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: OH + H2 <=> H + H2O */
    /*eqcon[20] *= 1;  */

    /*reaction 22: 2.000000 OH <=> O + H2O */
    /*eqcon[21] *= 1;  */

    /*reaction 23: OH + HO2 <=> O2 + H2O */
    /*eqcon[22] *= 1;  */

    /*reaction 24: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[23] *= 1;  */

    /*reaction 25: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[24] *= 1;  */

    /*reaction 26: 2.000000 HO2 <=> O2 + H2O2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: 2.000000 HO2 <=> O2 + H2O2 */
    /*eqcon[26] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mass fractions */
void CKEQYR(double *  rho, double *  T, double *  y, double *  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[9]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: 2.000000 OH (+M) <=> H2O2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: 2.000000 O + M <=> O2 + M */
    eqcon[1] *= 1e+06; 

    /*reaction 3: O + H + M <=> OH + M */
    eqcon[2] *= 1e+06; 

    /*reaction 4: H + O2 + M <=> HO2 + M */
    eqcon[3] *= 1e+06; 

    /*reaction 5: 2.000000 H + M <=> H2 + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: H + OH + M <=> H2O + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: O + H2 <=> H + OH */
    /*eqcon[6] *= 1;  */

    /*reaction 8: O + HO2 <=> OH + O2 */
    /*eqcon[7] *= 1;  */

    /*reaction 9: O + H2O2 <=> OH + HO2 */
    /*eqcon[8] *= 1;  */

    /*reaction 10: H + 2.000000 O2 <=> HO2 + O2 */
    eqcon[9] *= 1e+06; 

    /*reaction 11: H + O2 + H2O <=> HO2 + H2O */
    eqcon[10] *= 1e+06; 

    /*reaction 12: H + O2 + N2 <=> HO2 + N2 */
    eqcon[11] *= 1e+06; 

    /*reaction 13: H + O2 <=> O + OH */
    /*eqcon[12] *= 1;  */

    /*reaction 14: 2.000000 H + H2 <=> 2.000000 H2 */
    eqcon[13] *= 1e+06; 

    /*reaction 15: 2.000000 H + H2O <=> H2 + H2O */
    eqcon[14] *= 1e+06; 

    /*reaction 16: H + HO2 <=> O + H2O */
    /*eqcon[15] *= 1;  */

    /*reaction 17: H + HO2 <=> O2 + H2 */
    /*eqcon[16] *= 1;  */

    /*reaction 18: H + HO2 <=> 2.000000 OH */
    /*eqcon[17] *= 1;  */

    /*reaction 19: H + H2O2 <=> HO2 + H2 */
    /*eqcon[18] *= 1;  */

    /*reaction 20: H + H2O2 <=> OH + H2O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: OH + H2 <=> H + H2O */
    /*eqcon[20] *= 1;  */

    /*reaction 22: 2.000000 OH <=> O + H2O */
    /*eqcon[21] *= 1;  */

    /*reaction 23: OH + HO2 <=> O2 + H2O */
    /*eqcon[22] *= 1;  */

    /*reaction 24: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[23] *= 1;  */

    /*reaction 25: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[24] *= 1;  */

    /*reaction 26: 2.000000 HO2 <=> O2 + H2O2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: 2.000000 HO2 <=> O2 + H2O2 */
    /*eqcon[26] *= 1;  */
}


/*Returns the equil constants for each reaction */
/*Given rho, T, and mole fractions */
void CKEQXR(double *  rho, double *  T, double *  x, double *  eqcon)
{
    double tT = *T; /*temporary temperature */
    double tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    double gort[9]; /* temporary storage */

    /*compute the Gibbs free energy */
    gibbs(gort, tc);

    /*compute the equilibrium constants */
    equilibriumConstants(eqcon, gort, tT);

    /*reaction 1: 2.000000 OH (+M) <=> H2O2 (+M) */
    eqcon[0] *= 1e+06; 

    /*reaction 2: 2.000000 O + M <=> O2 + M */
    eqcon[1] *= 1e+06; 

    /*reaction 3: O + H + M <=> OH + M */
    eqcon[2] *= 1e+06; 

    /*reaction 4: H + O2 + M <=> HO2 + M */
    eqcon[3] *= 1e+06; 

    /*reaction 5: 2.000000 H + M <=> H2 + M */
    eqcon[4] *= 1e+06; 

    /*reaction 6: H + OH + M <=> H2O + M */
    eqcon[5] *= 1e+06; 

    /*reaction 7: O + H2 <=> H + OH */
    /*eqcon[6] *= 1;  */

    /*reaction 8: O + HO2 <=> OH + O2 */
    /*eqcon[7] *= 1;  */

    /*reaction 9: O + H2O2 <=> OH + HO2 */
    /*eqcon[8] *= 1;  */

    /*reaction 10: H + 2.000000 O2 <=> HO2 + O2 */
    eqcon[9] *= 1e+06; 

    /*reaction 11: H + O2 + H2O <=> HO2 + H2O */
    eqcon[10] *= 1e+06; 

    /*reaction 12: H + O2 + N2 <=> HO2 + N2 */
    eqcon[11] *= 1e+06; 

    /*reaction 13: H + O2 <=> O + OH */
    /*eqcon[12] *= 1;  */

    /*reaction 14: 2.000000 H + H2 <=> 2.000000 H2 */
    eqcon[13] *= 1e+06; 

    /*reaction 15: 2.000000 H + H2O <=> H2 + H2O */
    eqcon[14] *= 1e+06; 

    /*reaction 16: H + HO2 <=> O + H2O */
    /*eqcon[15] *= 1;  */

    /*reaction 17: H + HO2 <=> O2 + H2 */
    /*eqcon[16] *= 1;  */

    /*reaction 18: H + HO2 <=> 2.000000 OH */
    /*eqcon[17] *= 1;  */

    /*reaction 19: H + H2O2 <=> HO2 + H2 */
    /*eqcon[18] *= 1;  */

    /*reaction 20: H + H2O2 <=> OH + H2O */
    /*eqcon[19] *= 1;  */

    /*reaction 21: OH + H2 <=> H + H2O */
    /*eqcon[20] *= 1;  */

    /*reaction 22: 2.000000 OH <=> O + H2O */
    /*eqcon[21] *= 1;  */

    /*reaction 23: OH + HO2 <=> O2 + H2O */
    /*eqcon[22] *= 1;  */

    /*reaction 24: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[23] *= 1;  */

    /*reaction 25: OH + H2O2 <=> HO2 + H2O */
    /*eqcon[24] *= 1;  */

    /*reaction 26: 2.000000 HO2 <=> O2 + H2O2 */
    /*eqcon[25] *= 1;  */

    /*reaction 27: 2.000000 HO2 <=> O2 + H2O2 */
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

AMREX_GPU_HOST_DEVICE inline void comp_qfqr(double *  qf, double * qr, double * sc, double * tc, double invT)
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

    /*compute the mixture concentration */
    double mixture = 0.0;
    for (int i = 0; i < 9; ++i) {
        mixture += sc[i];
    }

    /*compute the Gibbs free energy */
    double g_RT[9];
    gibbs(g_RT, tc);

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    double refC = 101325 / 8.31446 * invT;
    double refCinv = 1 / refC;

    /* Evaluate the kfs */
    double k_f, Corr;
    double redP, F, logPred, logFcent, troe_c, troe_n, troe, F_troe;

    // (0):  2.000000 O + M <=> O2 + M
    k_f = 1.0000000000000002e-12 * 1.2e+17 
               * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = mixture + ( 2.3999999999999999 - 1)*sc[0] + ( 15.4 - 1)*sc[5];
    qf[1] *= Corr * k_f;
    qr[1] *= Corr * k_f / (exp(2.000000*g_RT[2] - g_RT[3]) * refCinv);
    // (1):  O + H + M <=> OH + M
    k_f = 1.0000000000000002e-12 * 5e+17 
               * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5];
    qf[2] *= Corr * k_f;
    qr[2] *= Corr * k_f / (exp(g_RT[1] + g_RT[2] - g_RT[4]) * refCinv);
    // (2):  O + H2 <=> H + OH
    k_f = 1.0000000000000002e-06 * 50000 
               * exp(2.6699999999999999 * tc[0] - 0.50321666580471969 * (6290) * invT);
    Corr  = 1.0;
    qf[6] *= Corr * k_f;
    qr[6] *= Corr * k_f / exp(g_RT[0] - g_RT[1] + g_RT[2] - g_RT[4]);
    // (3):  O + HO2 <=> OH + O2
    k_f = 1.0000000000000002e-06 * 20000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[7] *= Corr * k_f;
    qr[7] *= Corr * k_f / exp(g_RT[2] - g_RT[3] - g_RT[4] + g_RT[6]);
    // (4):  O + H2O2 <=> OH + HO2
    k_f = 1.0000000000000002e-06 * 9630000 
               * exp(2 * tc[0] - 0.50321666580471969 * (4000) * invT);
    Corr  = 1.0;
    qf[8] *= Corr * k_f;
    qr[8] *= Corr * k_f / exp(g_RT[2] - g_RT[4] - g_RT[6] + g_RT[7]);
    // (5):  H + O2 + M <=> HO2 + M
    k_f = 1.0000000000000002e-12 * 2.8e+18 
               * exp(-0.85999999999999999 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = mixture + ( 0 - 1)*sc[3] + ( 0 - 1)*sc[5] + ( 0 - 1)*sc[8];
    qf[3] *= Corr * k_f;
    qr[3] *= Corr * k_f / (exp(g_RT[1] + g_RT[3] - g_RT[6]) * refCinv);
    // (6):  H + 2.000000 O2 <=> HO2 + O2
    k_f = 1.0000000000000002e-12 * 3e+20 
               * exp(-1.72 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[9] *= Corr * k_f;
    qr[9] *= Corr * k_f / (exp(g_RT[1] + 2.000000*g_RT[3] - g_RT[3] - g_RT[6]) * refCinv);
    // (7):  H + O2 + H2O <=> HO2 + H2O
    k_f = 1.0000000000000002e-12 * 9.38e+18 
               * exp(-0.76000000000000001 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[10] *= Corr * k_f;
    qr[10] *= Corr * k_f / (exp(g_RT[1] + g_RT[3] + g_RT[5] - g_RT[5] - g_RT[6]) * refCinv);
    // (8):  H + O2 + N2 <=> HO2 + N2
    k_f = 1.0000000000000002e-12 * 3.75e+20 
               * exp(-1.72 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[11] *= Corr * k_f;
    qr[11] *= Corr * k_f / (exp(g_RT[1] + g_RT[3] - g_RT[6] + g_RT[8] - g_RT[8]) * refCinv);
    // (9):  H + O2 <=> O + OH
    k_f = 1.0000000000000002e-06 * 83000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (14413) * invT);
    Corr  = 1.0;
    qf[12] *= Corr * k_f;
    qr[12] *= Corr * k_f / exp(g_RT[1] - g_RT[2] + g_RT[3] - g_RT[4]);
    // (10):  2.000000 H + M <=> H2 + M
    k_f = 1.0000000000000002e-12 * 1e+18 
               * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = mixture + ( 0 - 1)*sc[0] + ( 0 - 1)*sc[5];
    qf[4] *= Corr * k_f;
    qr[4] *= Corr * k_f / (exp(-g_RT[0] + 2.000000*g_RT[1]) * refCinv);
    // (11):  2.000000 H + H2 <=> 2.000000 H2
    k_f = 1.0000000000000002e-12 * 90000000000000000 
               * exp(-0.59999999999999998 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[13] *= Corr * k_f;
    qr[13] *= Corr * k_f / (exp(g_RT[0] - 2.000000*g_RT[0] + 2.000000*g_RT[1]) * refCinv);
    // (12):  2.000000 H + H2O <=> H2 + H2O
    k_f = 1.0000000000000002e-12 * 6e+19 
               * exp(-1.25 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[14] *= Corr * k_f;
    qr[14] *= Corr * k_f / (exp(-g_RT[0] + 2.000000*g_RT[1] + g_RT[5] - g_RT[5]) * refCinv);
    // (13):  H + OH + M <=> H2O + M
    k_f = 1.0000000000000002e-12 * 2.2e+22 
               * exp(-2 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = mixture + ( 0.72999999999999998 - 1)*sc[0] + ( 3.6499999999999999 - 1)*sc[5];
    qf[5] *= Corr * k_f;
    qr[5] *= Corr * k_f / (exp(g_RT[1] + g_RT[4] - g_RT[5]) * refCinv);
    // (14):  H + HO2 <=> O + H2O
    k_f = 1.0000000000000002e-06 * 3970000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (671) * invT);
    Corr  = 1.0;
    qf[15] *= Corr * k_f;
    qr[15] *= Corr * k_f / exp(g_RT[1] - g_RT[2] - g_RT[5] + g_RT[6]);
    // (15):  H + HO2 <=> O2 + H2
    k_f = 1.0000000000000002e-06 * 28000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (1068) * invT);
    Corr  = 1.0;
    qf[16] *= Corr * k_f;
    qr[16] *= Corr * k_f / exp(-g_RT[0] + g_RT[1] - g_RT[3] + g_RT[6]);
    // (16):  H + HO2 <=> 2.000000 OH
    k_f = 1.0000000000000002e-06 * 134000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (635) * invT);
    Corr  = 1.0;
    qf[17] *= Corr * k_f;
    qr[17] *= Corr * k_f / exp(g_RT[1] - 2.000000*g_RT[4] + g_RT[6]);
    // (17):  H + H2O2 <=> HO2 + H2
    k_f = 1.0000000000000002e-06 * 12100000 
               * exp(2 * tc[0] - 0.50321666580471969 * (5200) * invT);
    Corr  = 1.0;
    qf[18] *= Corr * k_f;
    qr[18] *= Corr * k_f / exp(-g_RT[0] + g_RT[1] - g_RT[6] + g_RT[7]);
    // (18):  H + H2O2 <=> OH + H2O
    k_f = 1.0000000000000002e-06 * 10000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (3600) * invT);
    Corr  = 1.0;
    qf[19] *= Corr * k_f;
    qr[19] *= Corr * k_f / exp(g_RT[1] - g_RT[4] - g_RT[5] + g_RT[7]);
    // (19):  OH + H2 <=> H + H2O
    k_f = 1.0000000000000002e-06 * 216000000 
               * exp(1.51 * tc[0] - 0.50321666580471969 * (3430) * invT);
    Corr  = 1.0;
    qf[20] *= Corr * k_f;
    qr[20] *= Corr * k_f / exp(g_RT[0] - g_RT[1] + g_RT[4] - g_RT[5]);
    // (20):  2.000000 OH (+M) <=> H2O2 (+M)
    k_f = 1.0000000000000002e-06 * 74000000000000 
               * exp(-0.37 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5];
    redP = Corr / k_f * 1e-12 * 2.3e+18 
               * exp(-0.90000000000000002  * tc[0] - 0.50321666580471969  * (-1700) *invT);
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
    qf[0] *= Corr * k_f;
    qr[0] *= Corr * k_f / (exp(2.000000*g_RT[4] - g_RT[7]) * refCinv);
    // (21):  2.000000 OH <=> O + H2O
    k_f = 1.0000000000000002e-06 * 35700 
               * exp(2.3999999999999999 * tc[0] - 0.50321666580471969 * (-2110) * invT);
    Corr  = 1.0;
    qf[21] *= Corr * k_f;
    qr[21] *= Corr * k_f / exp(-g_RT[2] + 2.000000*g_RT[4] - g_RT[5]);
    // (22):  OH + HO2 <=> O2 + H2O
    k_f = 1.0000000000000002e-06 * 29000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (-500) * invT);
    Corr  = 1.0;
    qf[22] *= Corr * k_f;
    qr[22] *= Corr * k_f / exp(-g_RT[3] + g_RT[4] - g_RT[5] + g_RT[6]);
    // (23):  OH + H2O2 <=> HO2 + H2O
    k_f = 1.0000000000000002e-06 * 1750000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (320) * invT);
    Corr  = 1.0;
    qf[23] *= Corr * k_f;
    qr[23] *= Corr * k_f / exp(g_RT[4] - g_RT[5] - g_RT[6] + g_RT[7]);
    // (24):  OH + H2O2 <=> HO2 + H2O
    k_f = 1.0000000000000002e-06 * 580000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (9560) * invT);
    Corr  = 1.0;
    qf[24] *= Corr * k_f;
    qr[24] *= Corr * k_f / exp(g_RT[4] - g_RT[5] - g_RT[6] + g_RT[7]);
    // (25):  2.000000 HO2 <=> O2 + H2O2
    k_f = 1.0000000000000002e-06 * 130000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (-1630) * invT);
    Corr  = 1.0;
    qf[25] *= Corr * k_f;
    qr[25] *= Corr * k_f / exp(-g_RT[3] + 2.000000*g_RT[6] - g_RT[7]);
    // (26):  2.000000 HO2 <=> O2 + H2O2
    k_f = 1.0000000000000002e-06 * 420000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (12000) * invT);
    Corr  = 1.0;
    qf[26] *= Corr * k_f;
    qr[26] *= Corr * k_f / exp(-g_RT[3] + 2.000000*g_RT[6] - g_RT[7]);


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
    double g_RT[9];
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

void comp_qfqr(double *  qf, double *  qr, double *  sc, double *  tc, double invT)
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

    double T = tc[1];

    /*compute the mixture concentration */
    double mixture = 0.0;
    for (int i = 0; i < 9; ++i) {
        mixture += sc[i];
    }

    double Corr[27];
    for (int i = 0; i < 27; ++i) {
        Corr[i] = 1.0;
    }

    /* troe */
    {
        double alpha[1];
        alpha[0] = mixture + (TB[0][0] - 1)*sc[0] + (TB[0][1] - 1)*sc[5];
        for (int i=0; i<1; i++)
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
#endif


#ifndef AMREX_USE_CUDA
/*compute the production rate for each species */
void vproductionRate(int npt, double *  wdot, double *  sc, double *  T)
{
    double k_f_s[27*npt], Kc_s[27*npt], mixture[npt], g_RT[9*npt];
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
        double tg[5], g[9];
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

void vcomp_Kc(int npt, double *  Kc_s, double *  g_RT, double *  invT)
{
#ifdef __INTEL_COMPILER
    #pragma simd
#endif
    for (int i=0; i<npt; i++) {
        /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
        double refC = (101325. / 8.31451) * invT[i];
        double refCinv = 1.0 / refC;

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
#endif

/*compute an approx to the reaction Jacobian (for preconditioning) */
AMREX_GPU_HOST_DEVICE void DWDOT_SIMPLIFIED(double *  J, double *  sc, double *  Tp, int * HP)
{
    double c[9];

    for (int k=0; k<9; k++) {
        c[k] = 1.e6 * sc[k];
    }

    aJacobian_precond(J, c, *Tp, *HP);

    /* dwdot[k]/dT */
    /* dTdot/d[X] */
    for (int k=0; k<9; k++) {
        J[90+k] *= 1.e-6;
        J[k*10+9] *= 1.e6;
    }

    return;
}

/*compute the reaction Jacobian */
AMREX_GPU_HOST_DEVICE void DWDOT(double *  J, double *  sc, double *  Tp, int * consP)
{
    double c[9];

    for (int k=0; k<9; k++) {
        c[k] = 1.e6 * sc[k];
    }

    aJacobian(J, c, *Tp, *consP);

    /* dwdot[k]/dT */
    /* dTdot/d[X] */
    for (int k=0; k<9; k++) {
        J[90+k] *= 1.e-6;
        J[k*10+9] *= 1.e6;
    }

    return;
}

/*compute the sparsity pattern of the chemistry Jacobian */
AMREX_GPU_HOST_DEVICE void SPARSITY_INFO( int * nJdata, int * consP, int NCELLS)
{
    double c[9];
    double J[100];

    for (int k=0; k<9; k++) {
        c[k] = 1.0/ 9.000000 ;
    }

    aJacobian(J, c, 1500.0, *consP);

    int nJdata_tmp = 0;
    for (int k=0; k<10; k++) {
        for (int l=0; l<10; l++) {
            if(J[ 10 * k + l] != 0.0){
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
    double c[9];
    double J[100];

    for (int k=0; k<9; k++) {
        c[k] = 1.0/ 9.000000 ;
    }

    aJacobian(J, c, 1500.0, *consP);

    int nJdata_tmp = 0;
    for (int k=0; k<10; k++) {
        for (int l=0; l<10; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J[ 10 * k + l] != 0.0){
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
    double c[9];
    double J[100];

    for (int k=0; k<9; k++) {
        c[k] = 1.0/ 9.000000 ;
    }

    aJacobian_precond(J, c, 1500.0, *consP);

    int nJdata_tmp = 0;
    for (int k=0; k<10; k++) {
        for (int l=0; l<10; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J[ 10 * k + l] != 0.0){
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
    double c[9];
    double J[100];
    int offset_row;
    int offset_col;

    for (int k=0; k<9; k++) {
        c[k] = 1.0/ 9.000000 ;
    }

    aJacobian(J, c, 1500.0, *consP);

    colPtrs[0] = 0;
    int nJdata_tmp = 0;
    for (int nc=0; nc<NCELLS; nc++) {
        offset_row = nc * 10;
        offset_col = nc * 10;
        for (int k=0; k<10; k++) {
            for (int l=0; l<10; l++) {
                if(J[10*k + l] != 0.0) {
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
AMREX_GPU_HOST_DEVICE void SPARSITY_PREPROC_CSR(int *  colVals, int *  rowPtrs, int * consP, int NCELLS)
{
    double c[9];
    double J[100];
    int offset;

    for (int k=0; k<9; k++) {
        c[k] = 1.0/ 9.000000 ;
    }

    aJacobian(J, c, 1500.0, *consP);

    rowPtrs[0] = 0;
    int nJdata_tmp = 0;
    for (int nc=0; nc<NCELLS; nc++) {
        offset = nc * 10;
        for (int l=0; l<10; l++) {
            for (int k=0; k<10; k++) {
                if(J[10*k + l] != 0.0) {
                    colVals[nJdata_tmp] = k + offset; 
                    nJdata_tmp = nJdata_tmp + 1; 
                }
            }
            rowPtrs[offset + (l + 1)] = nJdata_tmp;
        }
    }

    return;
}

/*compute the sparsity pattern of the system Jacobian */
/*CSR format BASE is user choice */
AMREX_GPU_HOST_DEVICE void SPARSITY_PREPROC_SYST_CSR(int * colVals, int * rowPtr, int * consP, int NCELLS, int base)
{
    double c[9];
    double J[100];
    int offset;

    for (int k=0; k<9; k++) {
        c[k] = 1.0/ 9.000000 ;
    }

    aJacobian(J, c, 1500.0, *consP);

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
                        if(J[10*k + l] != 0.0) {
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
                        if(J[10*k + l] != 0.0) {
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
    double c[9];
    double J[100];

    for (int k=0; k<9; k++) {
        c[k] = 1.0/ 9.000000 ;
    }

    aJacobian_precond(J, c, 1500.0, *consP);

    colPtrs[0] = 0;
    int nJdata_tmp = 0;
    for (int k=0; k<10; k++) {
        for (int l=0; l<10; l++) {
            if (k == l) {
                rowVals[nJdata_tmp] = l; 
                indx[nJdata_tmp] = 10*k + l;
                nJdata_tmp = nJdata_tmp + 1; 
            } else {
                if(J[10*k + l] != 0.0) {
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
AMREX_GPU_HOST_DEVICE void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(int * colVals, int * rowPtr, int * consP, int base)
{
    double c[9];
    double J[100];

    for (int k=0; k<9; k++) {
        c[k] = 1.0/ 9.000000 ;
    }

    aJacobian_precond(J, c, 1500.0, *consP);

    if (base == 1) {
        rowPtr[0] = 1;
        int nJdata_tmp = 1;
        for (int l=0; l<10; l++) {
            for (int k=0; k<10; k++) {
                if (k == l) {
                    colVals[nJdata_tmp-1] = l+1; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J[10*k + l] != 0.0) {
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
                    if(J[10*k + l] != 0.0) {
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


    for (int i=0; i<100; i++) {
        J[i] = 0.0;
    }

    double wdot[9];
    for (int k=0; k<9; k++) {
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
    for (int k = 0; k < 9; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    double g_RT[9];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    double h_RT[9];
    speciesEnthalpy(h_RT, tc);

    double phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    double dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    double dqdci, dcdc_fac, dqdc[9];
    double Pr, fPr, F, k_0, logPr;
    double logFcent, troe_c, troe_n, troePr_den, troePr, troe;
    double Fcent1, Fcent2, Fcent3, Fcent;
    double dlogFdc, dlogFdn, dlogFdcn_fac;
    double dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
    const double ln10 = log(10.0);
    const double log10e = 1.0/log(10.0);
    /*reaction 1: 2.000000 OH (+M) <=> H2O2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5];
    /* forward */
    phi_f = pow(sc[4], 2.000000);
    k_f = 1.0000000000000002e-06 * 74000000000000
                * exp(-0.37 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -0.37 * invT + 0.50321666580471969 *  0  * invT2;
    /* pressure-fall-off */
    k_0 = 2.3e+18 * exp(-0.90000000000000002 * tc[0] - 0.50321666580471969 * (-1700) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -0.90000000000000002 * invT + 0.50321666580471969 * (-1700) * invT2;
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
        dqdci = (2 - 1)*dcdc_fac;
        J[4] += -2 * dqdci;           /* dwdot[OH]/d[H2] */
        J[7] += dqdci;                /* dwdot[H2O2]/d[H2] */
        /* d()/d[OH] */
        dqdci =  + k_f*2.000000*sc[4];
        J[44] += -2 * dqdci;          /* dwdot[OH]/d[OH] */
        J[47] += dqdci;               /* dwdot[H2O2]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (6 - 1)*dcdc_fac;
        J[54] += -2 * dqdci;          /* dwdot[OH]/d[H2O] */
        J[57] += dqdci;               /* dwdot[H2O2]/d[H2O] */
        /* d()/d[H2O2] */
        dqdci =  - k_r;
        J[74] += -2 * dqdci;          /* dwdot[OH]/d[H2O2] */
        J[77] += dqdci;               /* dwdot[H2O2]/d[H2O2] */
    }
    else {
        dqdc[0] = 2*dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac + k_f*2.000000*sc[4];
        dqdc[5] = 6*dcdc_fac;
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
    alpha = mixture + ( 2.3999999999999999 - 1)*sc[0] + ( 15.4 - 1)*sc[5];
    /* forward */
    phi_f = pow(sc[2], 2.000000);
    k_f = 1.0000000000000002e-12 * 1.2e+17
                * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1 * invT + 0.50321666580471969 *  0  * invT2;
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
        dqdci = (2.3999999999999999 - 1)*q_nocor;
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
        dqdci = (15.4 - 1)*q_nocor;
        J[52] += -2 * dqdci;          /* dwdot[O]/d[H2O] */
        J[53] += dqdci;               /* dwdot[O2]/d[H2O] */
    }
    else {
        dqdc[0] = 2.3999999999999999*q_nocor;
        dqdc[1] = q_nocor;
        dqdc[2] = q_nocor + k_f*2.000000*sc[2];
        dqdc[3] = q_nocor - k_r;
        dqdc[4] = q_nocor;
        dqdc[5] = 15.4*q_nocor;
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
    alpha = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5];
    /* forward */
    phi_f = sc[1]*sc[2];
    k_f = 1.0000000000000002e-12 * 5e+17
                * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1 * invT + 0.50321666580471969 *  0  * invT2;
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
        dqdci = (2 - 1)*q_nocor;
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
        dqdci = (6 - 1)*q_nocor;
        J[51] -= dqdci;               /* dwdot[H]/d[H2O] */
        J[52] -= dqdci;               /* dwdot[O]/d[H2O] */
        J[54] += dqdci;               /* dwdot[OH]/d[H2O] */
    }
    else {
        dqdc[0] = 2*q_nocor;
        dqdc[1] = q_nocor + k_f*sc[2];
        dqdc[2] = q_nocor + k_f*sc[1];
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor - k_r;
        dqdc[5] = 6*q_nocor;
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
    alpha = mixture + ( 0 - 1)*sc[3] + ( 0 - 1)*sc[5] + ( 0 - 1)*sc[8];
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = 1.0000000000000002e-12 * 2.8e+18
                * exp(-0.85999999999999999 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -0.85999999999999999 * invT + 0.50321666580471969 *  0  * invT2;
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
        dqdci = (0 - 1)*q_nocor + k_f*sc[1];
        J[31] -= dqdci;               /* dwdot[H]/d[O2] */
        J[33] -= dqdci;               /* dwdot[O2]/d[O2] */
        J[36] += dqdci;               /* dwdot[HO2]/d[O2] */
        /* d()/d[H2O] */
        dqdci = (0 - 1)*q_nocor;
        J[51] -= dqdci;               /* dwdot[H]/d[H2O] */
        J[53] -= dqdci;               /* dwdot[O2]/d[H2O] */
        J[56] += dqdci;               /* dwdot[HO2]/d[H2O] */
        /* d()/d[HO2] */
        dqdci =  - k_r;
        J[61] -= dqdci;               /* dwdot[H]/d[HO2] */
        J[63] -= dqdci;               /* dwdot[O2]/d[HO2] */
        J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
        /* d()/d[N2] */
        dqdci = (0 - 1)*q_nocor;
        J[81] -= dqdci;               /* dwdot[H]/d[N2] */
        J[83] -= dqdci;               /* dwdot[O2]/d[N2] */
        J[86] += dqdci;               /* dwdot[HO2]/d[N2] */
    }
    else {
        dqdc[0] = q_nocor;
        dqdc[1] = q_nocor + k_f*sc[3];
        dqdc[2] = q_nocor;
        dqdc[3] =  + k_f*sc[1];
        dqdc[4] = q_nocor;
        dqdc[6] = q_nocor - k_r;
        dqdc[7] = q_nocor;
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
    alpha = mixture + ( 0 - 1)*sc[0] + ( 0 - 1)*sc[5];
    /* forward */
    phi_f = pow(sc[1], 2.000000);
    k_f = 1.0000000000000002e-12 * 1e+18
                * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1 * invT + 0.50321666580471969 *  0  * invT2;
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
        dqdci = (0 - 1)*q_nocor - k_r;
        J[0] += dqdci;                /* dwdot[H2]/d[H2] */
        J[1] += -2 * dqdci;           /* dwdot[H]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*2.000000*sc[1];
        J[10] += dqdci;               /* dwdot[H2]/d[H] */
        J[11] += -2 * dqdci;          /* dwdot[H]/d[H] */
        /* d()/d[H2O] */
        dqdci = (0 - 1)*q_nocor;
        J[50] += dqdci;               /* dwdot[H2]/d[H2O] */
        J[51] += -2 * dqdci;          /* dwdot[H]/d[H2O] */
    }
    else {
        dqdc[0] =  - k_r;
        dqdc[1] = q_nocor + k_f*2.000000*sc[1];
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor;
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
    alpha = mixture + ( 0.72999999999999998 - 1)*sc[0] + ( 3.6499999999999999 - 1)*sc[5];
    /* forward */
    phi_f = sc[1]*sc[4];
    k_f = 1.0000000000000002e-12 * 2.2e+22
                * exp(-2 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -2 * invT + 0.50321666580471969 *  0  * invT2;
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
        dqdci = (0.72999999999999998 - 1)*q_nocor;
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
        dqdci = (3.6499999999999999 - 1)*q_nocor - k_r;
        J[51] -= dqdci;               /* dwdot[H]/d[H2O] */
        J[54] -= dqdci;               /* dwdot[OH]/d[H2O] */
        J[55] += dqdci;               /* dwdot[H2O]/d[H2O] */
    }
    else {
        dqdc[0] = 0.72999999999999998*q_nocor;
        dqdc[1] = q_nocor + k_f*sc[4];
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor + k_f*sc[1];
        dqdc[5] = 3.6499999999999999*q_nocor - k_r;
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
    k_f = 1.0000000000000002e-06 * 50000
                * exp(2.6699999999999999 * tc[0] - 0.50321666580471969 * (6290) * invT);
    dlnkfdT = 2.6699999999999999 * invT + 0.50321666580471969 *  6290  * invT2;
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
    k_f = 1.0000000000000002e-06 * 20000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
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
    k_f = 1.0000000000000002e-06 * 9630000
                * exp(2 * tc[0] - 0.50321666580471969 * (4000) * invT);
    dlnkfdT = 2 * invT + 0.50321666580471969 *  4000  * invT2;
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
    k_f = 1.0000000000000002e-12 * 3e+20
                * exp(-1.72 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1.72 * invT + 0.50321666580471969 *  0  * invT2;
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
    k_f = 1.0000000000000002e-12 * 9.38e+18
                * exp(-0.76000000000000001 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -0.76000000000000001 * invT + 0.50321666580471969 *  0  * invT2;
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
    k_f = 1.0000000000000002e-12 * 3.75e+20
                * exp(-1.72 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1.72 * invT + 0.50321666580471969 *  0  * invT2;
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
    k_f = 1.0000000000000002e-06 * 83000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (14413) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  14413  * invT2;
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
    k_f = 1.0000000000000002e-12 * 90000000000000000
                * exp(-0.59999999999999998 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -0.59999999999999998 * invT + 0.50321666580471969 *  0  * invT2;
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
    k_f = 1.0000000000000002e-12 * 6e+19
                * exp(-1.25 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1.25 * invT + 0.50321666580471969 *  0  * invT2;
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
    k_f = 1.0000000000000002e-06 * 3970000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (671) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  671  * invT2;
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
    k_f = 1.0000000000000002e-06 * 28000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (1068) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  1068  * invT2;
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
    k_f = 1.0000000000000002e-06 * 134000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (635) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  635  * invT2;
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
    k_f = 1.0000000000000002e-06 * 12100000
                * exp(2 * tc[0] - 0.50321666580471969 * (5200) * invT);
    dlnkfdT = 2 * invT + 0.50321666580471969 *  5200  * invT2;
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
    k_f = 1.0000000000000002e-06 * 10000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (3600) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  3600  * invT2;
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
    k_f = 1.0000000000000002e-06 * 216000000
                * exp(1.51 * tc[0] - 0.50321666580471969 * (3430) * invT);
    dlnkfdT = 1.51 * invT + 0.50321666580471969 *  3430  * invT2;
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
    k_f = 1.0000000000000002e-06 * 35700
                * exp(2.3999999999999999 * tc[0] - 0.50321666580471969 * (-2110) * invT);
    dlnkfdT = 2.3999999999999999 * invT + 0.50321666580471969 *  -2110  * invT2;
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
    k_f = 1.0000000000000002e-06 * 29000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-500) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  -500  * invT2;
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
    k_f = 1.0000000000000002e-06 * 1750000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (320) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  320  * invT2;
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
    k_f = 1.0000000000000002e-06 * 580000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (9560) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  9560  * invT2;
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
    k_f = 1.0000000000000002e-06 * 130000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-1630) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  -1630  * invT2;
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
    k_f = 1.0000000000000002e-06 * 420000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (12000) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  12000  * invT2;
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

    double c_R[9], dcRdT[9], e_RT[9];
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
    for (int k = 0; k < 9; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[90+k];
    }

    double cmixinv = 1.0/cmix;
    double tmp1 = ehmix*cmixinv;
    double tmp3 = cmixinv*T;
    double tmp2 = tmp1*tmp3;
    double dehmixdc;
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

return;
}
#endif


#ifndef AMREX_USE_CUDA
/*compute the reaction Jacobian on CPU */
void aJacobian(double *  J, double *  sc, double T, int consP)
{
    for (int i=0; i<100; i++) {
        J[i] = 0.0;
    }

    double wdot[9];
    for (int k=0; k<9; k++) {
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
    for (int k = 0; k < 9; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    double g_RT[9];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    double h_RT[9];
    speciesEnthalpy(h_RT, tc);

    double phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    double dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    double dqdci, dcdc_fac, dqdc[9];
    double Pr, fPr, F, k_0, logPr;
    double logFcent, troe_c, troe_n, troePr_den, troePr, troe;
    double Fcent1, Fcent2, Fcent3, Fcent;
    double dlogFdc, dlogFdn, dlogFdcn_fac;
    double dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
    const double ln10 = log(10.0);
    const double log10e = 1.0/log(10.0);
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

    double c_R[9], dcRdT[9], e_RT[9];
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
    for (int k = 0; k < 9; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[90+k];
    }

    double cmixinv = 1.0/cmix;
    double tmp1 = ehmix*cmixinv;
    double tmp3 = cmixinv*T;
    double tmp2 = tmp1*tmp3;
    double dehmixdc;
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


/*compute an approx to the reaction Jacobian */
AMREX_GPU_HOST_DEVICE void aJacobian_precond(double *  J, double *  sc, double T, int HP)
{
    for (int i=0; i<100; i++) {
        J[i] = 0.0;
    }

    double wdot[9];
    for (int k=0; k<9; k++) {
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
    for (int k = 0; k < 9; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    double g_RT[9];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    double h_RT[9];
    speciesEnthalpy(h_RT, tc);

    double phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    double dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    double dqdci, dcdc_fac, dqdc[9];
    double Pr, fPr, F, k_0, logPr;
    double logFcent, troe_c, troe_n, troePr_den, troePr, troe;
    double Fcent1, Fcent2, Fcent3, Fcent;
    double dlogFdc, dlogFdn, dlogFdcn_fac;
    double dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
    const double ln10 = log(10.0);
    const double log10e = 1.0/log(10.0);
    /*reaction 1: 2.000000 OH (+M) <=> H2O2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5];
    /* forward */
    phi_f = pow(sc[4], 2.000000);
    k_f = 1.0000000000000002e-06 * 74000000000000
                * exp(-0.37 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -0.37 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* pressure-fall-off */
    k_0 = 2.3e+18 * exp(-0.90000000000000002 * tc[0] - 0.50321666580471969 * (-1700) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -0.90000000000000002 * invT + 0.50321666580471969 * (-1700) * invT2;
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
    dcdc_fac = 0.0;
    dqdc[0] = 2*dcdc_fac;
    dqdc[1] = dcdc_fac;
    dqdc[2] = dcdc_fac;
    dqdc[3] = dcdc_fac;
    dqdc[4] = dcdc_fac + k_f*2.000000*sc[4];
    dqdc[5] = 6*dcdc_fac;
    dqdc[6] = dcdc_fac;
    dqdc[7] = dcdc_fac - k_r;
    dqdc[8] = dcdc_fac;
    for (int k=0; k<9; k++) {
        J[10*k+4] += -2 * dqdc[k];
        J[10*k+7] += dqdc[k];
    }
    J[94] += -2 * dqdT; /* dwdot[OH]/dT */
    J[97] += dqdT; /* dwdot[H2O2]/dT */

    /*reaction 2: 2.000000 O + M <=> O2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 2.3999999999999999 - 1)*sc[0] + ( 15.4 - 1)*sc[5];
    /* forward */
    phi_f = pow(sc[2], 2.000000);
    k_f = 1.0000000000000002e-12 * 1.2e+17
                * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1 * invT + 0.50321666580471969 *  (0)  * invT2;
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
    dqdc[0] = 2.3999999999999999*q_nocor;
    dqdc[1] = q_nocor;
    dqdc[2] = q_nocor + k_f*2.000000*sc[2];
    dqdc[3] = q_nocor - k_r;
    dqdc[4] = q_nocor;
    dqdc[5] = 15.4*q_nocor;
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    for (int k=0; k<9; k++) {
        J[10*k+2] += -2 * dqdc[k];
        J[10*k+3] += dqdc[k];
    }
    J[92] += -2 * dqdT; /* dwdot[O]/dT */
    J[93] += dqdT; /* dwdot[O2]/dT */

    /*reaction 3: O + H + M <=> OH + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5];
    /* forward */
    phi_f = sc[1]*sc[2];
    k_f = 1.0000000000000002e-12 * 5e+17
                * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1 * invT + 0.50321666580471969 *  (0)  * invT2;
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
    dqdc[0] = 2*q_nocor;
    dqdc[1] = q_nocor + k_f*sc[2];
    dqdc[2] = q_nocor + k_f*sc[1];
    dqdc[3] = q_nocor;
    dqdc[4] = q_nocor - k_r;
    dqdc[5] = 6*q_nocor;
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    for (int k=0; k<9; k++) {
        J[10*k+1] -= dqdc[k];
        J[10*k+2] -= dqdc[k];
        J[10*k+4] += dqdc[k];
    }
    J[91] -= dqdT; /* dwdot[H]/dT */
    J[92] -= dqdT; /* dwdot[O]/dT */
    J[94] += dqdT; /* dwdot[OH]/dT */

    /*reaction 4: H + O2 + M <=> HO2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 0 - 1)*sc[3] + ( 0 - 1)*sc[5] + ( 0 - 1)*sc[8];
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = 1.0000000000000002e-12 * 2.8e+18
                * exp(-0.85999999999999999 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -0.85999999999999999 * invT + 0.50321666580471969 *  (0)  * invT2;
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
    dqdc[0] = q_nocor;
    dqdc[1] = q_nocor + k_f*sc[3];
    dqdc[2] = q_nocor;
    dqdc[3] =  + k_f*sc[1];
    dqdc[4] = q_nocor;
    dqdc[5] = 0.0;
    dqdc[6] = q_nocor - k_r;
    dqdc[7] = q_nocor;
    dqdc[8] = 0.0;
    for (int k=0; k<9; k++) {
        J[10*k+1] -= dqdc[k];
        J[10*k+3] -= dqdc[k];
        J[10*k+6] += dqdc[k];
    }
    J[91] -= dqdT; /* dwdot[H]/dT */
    J[93] -= dqdT; /* dwdot[O2]/dT */
    J[96] += dqdT; /* dwdot[HO2]/dT */

    /*reaction 5: 2.000000 H + M <=> H2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 0 - 1)*sc[0] + ( 0 - 1)*sc[5];
    /* forward */
    phi_f = pow(sc[1], 2.000000);
    k_f = 1.0000000000000002e-12 * 1e+18
                * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1 * invT + 0.50321666580471969 *  (0)  * invT2;
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
    dqdc[0] =  - k_r;
    dqdc[1] = q_nocor + k_f*2.000000*sc[1];
    dqdc[2] = q_nocor;
    dqdc[3] = q_nocor;
    dqdc[4] = q_nocor;
    dqdc[5] = 0.0;
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    for (int k=0; k<9; k++) {
        J[10*k+0] += dqdc[k];
        J[10*k+1] += -2 * dqdc[k];
    }
    J[90] += dqdT; /* dwdot[H2]/dT */
    J[91] += -2 * dqdT; /* dwdot[H]/dT */

    /*reaction 6: H + OH + M <=> H2O + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 0.72999999999999998 - 1)*sc[0] + ( 3.6499999999999999 - 1)*sc[5];
    /* forward */
    phi_f = sc[1]*sc[4];
    k_f = 1.0000000000000002e-12 * 2.2e+22
                * exp(-2 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -2 * invT + 0.50321666580471969 *  (0)  * invT2;
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
    dqdc[0] = 0.72999999999999998*q_nocor;
    dqdc[1] = q_nocor + k_f*sc[4];
    dqdc[2] = q_nocor;
    dqdc[3] = q_nocor;
    dqdc[4] = q_nocor + k_f*sc[1];
    dqdc[5] = 3.6499999999999999*q_nocor - k_r;
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    for (int k=0; k<9; k++) {
        J[10*k+1] -= dqdc[k];
        J[10*k+4] -= dqdc[k];
        J[10*k+5] += dqdc[k];
    }
    J[91] -= dqdT; /* dwdot[H]/dT */
    J[94] -= dqdT; /* dwdot[OH]/dT */
    J[95] += dqdT; /* dwdot[H2O]/dT */

    /*reaction 7: O + H2 <=> H + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[2];
    k_f = 1.0000000000000002e-06 * 50000
                * exp(2.6699999999999999 * tc[0] - 0.50321666580471969 * (6290) * invT);
    dlnkfdT = 2.6699999999999999 * invT + 0.50321666580471969 *  (6290)  * invT2;
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
    k_f = 1.0000000000000002e-06 * 20000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
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
    k_f = 1.0000000000000002e-06 * 9630000
                * exp(2 * tc[0] - 0.50321666580471969 * (4000) * invT);
    dlnkfdT = 2 * invT + 0.50321666580471969 *  (4000)  * invT2;
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
    k_f = 1.0000000000000002e-12 * 3e+20
                * exp(-1.72 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1.72 * invT + 0.50321666580471969 *  (0)  * invT2;
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
    k_f = 1.0000000000000002e-12 * 9.38e+18
                * exp(-0.76000000000000001 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -0.76000000000000001 * invT + 0.50321666580471969 *  (0)  * invT2;
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
    k_f = 1.0000000000000002e-12 * 3.75e+20
                * exp(-1.72 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1.72 * invT + 0.50321666580471969 *  (0)  * invT2;
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
    k_f = 1.0000000000000002e-06 * 83000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (14413) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (14413)  * invT2;
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
    k_f = 1.0000000000000002e-12 * 90000000000000000
                * exp(-0.59999999999999998 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -0.59999999999999998 * invT + 0.50321666580471969 *  (0)  * invT2;
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
    k_f = 1.0000000000000002e-12 * 6e+19
                * exp(-1.25 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1.25 * invT + 0.50321666580471969 *  (0)  * invT2;
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
    k_f = 1.0000000000000002e-06 * 3970000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (671) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (671)  * invT2;
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
    k_f = 1.0000000000000002e-06 * 28000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (1068) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (1068)  * invT2;
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
    k_f = 1.0000000000000002e-06 * 134000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (635) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (635)  * invT2;
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
    k_f = 1.0000000000000002e-06 * 12100000
                * exp(2 * tc[0] - 0.50321666580471969 * (5200) * invT);
    dlnkfdT = 2 * invT + 0.50321666580471969 *  (5200)  * invT2;
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
    k_f = 1.0000000000000002e-06 * 10000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (3600) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (3600)  * invT2;
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
    k_f = 1.0000000000000002e-06 * 216000000
                * exp(1.51 * tc[0] - 0.50321666580471969 * (3430) * invT);
    dlnkfdT = 1.51 * invT + 0.50321666580471969 *  (3430)  * invT2;
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
    k_f = 1.0000000000000002e-06 * 35700
                * exp(2.3999999999999999 * tc[0] - 0.50321666580471969 * (-2110) * invT);
    dlnkfdT = 2.3999999999999999 * invT + 0.50321666580471969 *  (-2110)  * invT2;
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
    k_f = 1.0000000000000002e-06 * 29000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-500) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (-500)  * invT2;
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
    k_f = 1.0000000000000002e-06 * 1750000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (320) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (320)  * invT2;
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
    k_f = 1.0000000000000002e-06 * 580000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (9560) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (9560)  * invT2;
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
    k_f = 1.0000000000000002e-06 * 130000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-1630) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (-1630)  * invT2;
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
    k_f = 1.0000000000000002e-06 * 420000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (12000) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (12000)  * invT2;
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

    double c_R[9], dcRdT[9], e_RT[9];
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
    for (int k = 0; k < 9; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[90+k];
    }

    double cmixinv = 1.0/cmix;
    double tmp1 = ehmix*cmixinv;
    double tmp3 = cmixinv*T;
    double tmp2 = tmp1*tmp3;
    double dehmixdc;
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
        /*species 2: O */
        species[2] =
            -3.27931884e-03
            +1.32861279e-05 * tc[1]
            -1.83841987e-08 * tc[2]
            +8.45063884e-12 * tc[3];
        /*species 3: O2 */
        species[3] =
            -2.99673416e-03
            +1.96946040e-05 * tc[1]
            -2.90438853e-08 * tc[2]
            +1.29749135e-11 * tc[3];
        /*species 4: OH */
        species[4] =
            -2.40131752e-03
            +9.23587682e-06 * tc[1]
            -1.16434000e-08 * tc[2]
            +5.45645880e-12 * tc[3];
        /*species 5: H2O */
        species[5] =
            -2.03643410e-03
            +1.30408042e-05 * tc[1]
            -1.64639119e-08 * tc[2]
            +7.08791268e-12 * tc[3];
        /*species 6: HO2 */
        species[6] =
            -4.74912051e-03
            +4.23165782e-05 * tc[1]
            -7.28291682e-08 * tc[2]
            +3.71690050e-11 * tc[3];
        /*species 7: H2O2 */
        species[7] =
            -5.42822417e-04
            +3.34671402e-05 * tc[1]
            -6.47312439e-08 * tc[2]
            +3.44981745e-11 * tc[3];
        /*species 8: N2 */
        species[8] =
            +1.40824040e-03
            -7.92644400e-06 * tc[1]
            +1.69245450e-08 * tc[2]
            -9.77941600e-12 * tc[3];
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
        /*species 2: O */
        species[2] =
            -8.59741137e-05
            +8.38969178e-08 * tc[1]
            -3.00533397e-11 * tc[2]
            +4.91334764e-15 * tc[3];
        /*species 3: O2 */
        species[3] =
            +1.48308754e-03
            -1.51593334e-06 * tc[1]
            +6.28411665e-10 * tc[2]
            -8.66871176e-14 * tc[3];
        /*species 4: OH */
        species[4] =
            +5.48429716e-04
            +2.53010456e-07 * tc[1]
            -2.63838467e-10 * tc[2]
            +4.69649504e-14 * tc[3];
        /*species 5: H2O */
        species[5] =
            +2.17691804e-03
            -3.28145036e-07 * tc[1]
            -2.91125961e-10 * tc[2]
            +6.72803968e-14 * tc[3];
        /*species 6: HO2 */
        species[6] =
            +2.23982013e-03
            -1.26731630e-06 * tc[1]
            +3.42739110e-10 * tc[2]
            -4.31634140e-14 * tc[3];
        /*species 7: H2O2 */
        species[7] =
            +4.90831694e-03
            -3.80278450e-06 * tc[1]
            +1.11355796e-09 * tc[2]
            -1.15163322e-13 * tc[3];
        /*species 8: N2 */
        species[8] =
            +1.48797680e-03
            -1.13695200e-06 * tc[1]
            +3.02911140e-10 * tc[2]
            -2.70134040e-14 * tc[3];
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

    /*reaction 1: 2.000000 OH (+M) <=> H2O2 (+M) */
    kc[0] = 1.0 / (refC) * exp((2.000000 * g_RT[4]) - (g_RT[7]));

    /*reaction 2: 2.000000 O + M <=> O2 + M */
    kc[1] = 1.0 / (refC) * exp((2.000000 * g_RT[2]) - (g_RT[3]));

    /*reaction 3: O + H + M <=> OH + M */
    kc[2] = 1.0 / (refC) * exp((g_RT[2] + g_RT[1]) - (g_RT[4]));

    /*reaction 4: H + O2 + M <=> HO2 + M */
    kc[3] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3]) - (g_RT[6]));

    /*reaction 5: 2.000000 H + M <=> H2 + M */
    kc[4] = 1.0 / (refC) * exp((2.000000 * g_RT[1]) - (g_RT[0]));

    /*reaction 6: H + OH + M <=> H2O + M */
    kc[5] = 1.0 / (refC) * exp((g_RT[1] + g_RT[4]) - (g_RT[5]));

    /*reaction 7: O + H2 <=> H + OH */
    kc[6] = exp((g_RT[2] + g_RT[0]) - (g_RT[1] + g_RT[4]));

    /*reaction 8: O + HO2 <=> OH + O2 */
    kc[7] = exp((g_RT[2] + g_RT[6]) - (g_RT[4] + g_RT[3]));

    /*reaction 9: O + H2O2 <=> OH + HO2 */
    kc[8] = exp((g_RT[2] + g_RT[7]) - (g_RT[4] + g_RT[6]));

    /*reaction 10: H + 2.000000 O2 <=> HO2 + O2 */
    kc[9] = 1.0 / (refC) * exp((g_RT[1] + 2.000000 * g_RT[3]) - (g_RT[6] + g_RT[3]));

    /*reaction 11: H + O2 + H2O <=> HO2 + H2O */
    kc[10] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[5]) - (g_RT[6] + g_RT[5]));

    /*reaction 12: H + O2 + N2 <=> HO2 + N2 */
    kc[11] = 1.0 / (refC) * exp((g_RT[1] + g_RT[3] + g_RT[8]) - (g_RT[6] + g_RT[8]));

    /*reaction 13: H + O2 <=> O + OH */
    kc[12] = exp((g_RT[1] + g_RT[3]) - (g_RT[2] + g_RT[4]));

    /*reaction 14: 2.000000 H + H2 <=> 2.000000 H2 */
    kc[13] = 1.0 / (refC) * exp((2.000000 * g_RT[1] + g_RT[0]) - (2.000000 * g_RT[0]));

    /*reaction 15: 2.000000 H + H2O <=> H2 + H2O */
    kc[14] = 1.0 / (refC) * exp((2.000000 * g_RT[1] + g_RT[5]) - (g_RT[0] + g_RT[5]));

    /*reaction 16: H + HO2 <=> O + H2O */
    kc[15] = exp((g_RT[1] + g_RT[6]) - (g_RT[2] + g_RT[5]));

    /*reaction 17: H + HO2 <=> O2 + H2 */
    kc[16] = exp((g_RT[1] + g_RT[6]) - (g_RT[3] + g_RT[0]));

    /*reaction 18: H + HO2 <=> 2.000000 OH */
    kc[17] = exp((g_RT[1] + g_RT[6]) - (2.000000 * g_RT[4]));

    /*reaction 19: H + H2O2 <=> HO2 + H2 */
    kc[18] = exp((g_RT[1] + g_RT[7]) - (g_RT[6] + g_RT[0]));

    /*reaction 20: H + H2O2 <=> OH + H2O */
    kc[19] = exp((g_RT[1] + g_RT[7]) - (g_RT[4] + g_RT[5]));

    /*reaction 21: OH + H2 <=> H + H2O */
    kc[20] = exp((g_RT[4] + g_RT[0]) - (g_RT[1] + g_RT[5]));

    /*reaction 22: 2.000000 OH <=> O + H2O */
    kc[21] = exp((2.000000 * g_RT[4]) - (g_RT[2] + g_RT[5]));

    /*reaction 23: OH + HO2 <=> O2 + H2O */
    kc[22] = exp((g_RT[4] + g_RT[6]) - (g_RT[3] + g_RT[5]));

    /*reaction 24: OH + H2O2 <=> HO2 + H2O */
    kc[23] = exp((g_RT[4] + g_RT[7]) - (g_RT[6] + g_RT[5]));

    /*reaction 25: OH + H2O2 <=> HO2 + H2O */
    kc[24] = exp((g_RT[4] + g_RT[7]) - (g_RT[6] + g_RT[5]));

    /*reaction 26: 2.000000 HO2 <=> O2 + H2O2 */
    kc[25] = exp((2.000000 * g_RT[6]) - (g_RT[3] + g_RT[7]));

    /*reaction 27: 2.000000 HO2 <=> O2 + H2O2 */
    kc[26] = exp((2.000000 * g_RT[6]) - (g_RT[3] + g_RT[7]));

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
        /*species 2: O */
        species[2] =
            +2.912225920000000e+04 * invT
            +1.116333640000000e+00
            -3.168267100000000e+00 * tc[0]
            +1.639659420000000e-03 * tc[1]
            -1.107177326666667e-06 * tc[2]
            +5.106721866666666e-10 * tc[3]
            -1.056329855000000e-13 * tc[4];
        /*species 3: O2 */
        species[3] =
            -1.063943560000000e+03 * invT
            +1.247806300000001e-01
            -3.782456360000000e+00 * tc[0]
            +1.498367080000000e-03 * tc[1]
            -1.641217001666667e-06 * tc[2]
            +8.067745908333334e-10 * tc[3]
            -1.621864185000000e-13 * tc[4];
        /*species 4: OH */
        species[4] =
            +3.615080560000000e+03 * invT
            +4.095940888000000e+00
            -3.992015430000000e+00 * tc[0]
            +1.200658760000000e-03 * tc[1]
            -7.696564016666666e-07 * tc[2]
            +3.234277775000000e-10 * tc[3]
            -6.820573500000000e-14 * tc[4];
        /*species 5: H2O */
        species[5] =
            -3.029372670000000e+04 * invT
            +5.047672768000000e+00
            -4.198640560000000e+00 * tc[0]
            +1.018217050000000e-03 * tc[1]
            -1.086733685000000e-06 * tc[2]
            +4.573308850000000e-10 * tc[3]
            -8.859890850000000e-14 * tc[4];
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
            -1.770258210000000e+04 * invT
            +8.410619499999998e-01
            -4.276112690000000e+00 * tc[0]
            +2.714112085000000e-04 * tc[1]
            -2.788928350000000e-06 * tc[2]
            +1.798090108333333e-09 * tc[3]
            -4.312271815000000e-13 * tc[4];
        /*species 8: N2 */
        species[8] =
            -1.020899900000000e+03 * invT
            -6.516950000000001e-01
            -3.298677000000000e+00 * tc[0]
            -7.041202000000000e-04 * tc[1]
            +6.605369999999999e-07 * tc[2]
            -4.701262500000001e-10 * tc[3]
            +1.222427000000000e-13 * tc[4];
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
        /*species 2: O */
        species[2] =
            +2.921757910000000e+04 * invT
            -2.214917859999999e+00
            -2.569420780000000e+00 * tc[0]
            +4.298705685000000e-05 * tc[1]
            -6.991409816666667e-09 * tc[2]
            +8.348149916666666e-13 * tc[3]
            -6.141684549999999e-17 * tc[4];
        /*species 3: O2 */
        species[3] =
            -1.088457720000000e+03 * invT
            -2.170693450000000e+00
            -3.282537840000000e+00 * tc[0]
            -7.415437700000000e-04 * tc[1]
            +1.263277781666667e-07 * tc[2]
            -1.745587958333333e-11 * tc[3]
            +1.083588970000000e-15 * tc[4];
        /*species 4: OH */
        species[4] =
            +3.858657000000000e+03 * invT
            -1.383808430000000e+00
            -3.092887670000000e+00 * tc[0]
            -2.742148580000000e-04 * tc[1]
            -2.108420466666667e-08 * tc[2]
            +7.328846300000000e-12 * tc[3]
            -5.870618800000000e-16 * tc[4];
        /*species 5: H2O */
        species[5] =
            -3.000429710000000e+04 * invT
            -1.932777610000000e+00
            -3.033992490000000e+00 * tc[0]
            -1.088459020000000e-03 * tc[1]
            +2.734541966666666e-08 * tc[2]
            +8.086832250000000e-12 * tc[3]
            -8.410049600000000e-16 * tc[4];
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
            -1.786178770000000e+04 * invT
            +1.248846229999999e+00
            -4.165002850000000e+00 * tc[0]
            -2.454158470000000e-03 * tc[1]
            +3.168987083333333e-07 * tc[2]
            -3.093216550000000e-11 * tc[3]
            +1.439541525000000e-15 * tc[4];
        /*species 8: N2 */
        species[8] =
            -9.227977000000000e+02 * invT
            -3.053888000000000e+00
            -2.926640000000000e+00 * tc[0]
            -7.439884000000000e-04 * tc[1]
            +9.474600000000001e-08 * tc[2]
            -8.414198333333333e-12 * tc[3]
            +3.376675500000000e-16 * tc[4];
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
        /*species 2: O */
        species[2] =
            +2.91222592e+04 * invT
            +1.16333640e-01
            -3.16826710e+00 * tc[0]
            +1.63965942e-03 * tc[1]
            -1.10717733e-06 * tc[2]
            +5.10672187e-10 * tc[3]
            -1.05632985e-13 * tc[4];
        /*species 3: O2 */
        species[3] =
            -1.06394356e+03 * invT
            -8.75219370e-01
            -3.78245636e+00 * tc[0]
            +1.49836708e-03 * tc[1]
            -1.64121700e-06 * tc[2]
            +8.06774591e-10 * tc[3]
            -1.62186418e-13 * tc[4];
        /*species 4: OH */
        species[4] =
            +3.61508056e+03 * invT
            +3.09594089e+00
            -3.99201543e+00 * tc[0]
            +1.20065876e-03 * tc[1]
            -7.69656402e-07 * tc[2]
            +3.23427778e-10 * tc[3]
            -6.82057350e-14 * tc[4];
        /*species 5: H2O */
        species[5] =
            -3.02937267e+04 * invT
            +4.04767277e+00
            -4.19864056e+00 * tc[0]
            +1.01821705e-03 * tc[1]
            -1.08673369e-06 * tc[2]
            +4.57330885e-10 * tc[3]
            -8.85989085e-14 * tc[4];
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
            -1.77025821e+04 * invT
            -1.58938050e-01
            -4.27611269e+00 * tc[0]
            +2.71411208e-04 * tc[1]
            -2.78892835e-06 * tc[2]
            +1.79809011e-09 * tc[3]
            -4.31227182e-13 * tc[4];
        /*species 8: N2 */
        species[8] =
            -1.02089990e+03 * invT
            -1.65169500e+00
            -3.29867700e+00 * tc[0]
            -7.04120200e-04 * tc[1]
            +6.60537000e-07 * tc[2]
            -4.70126250e-10 * tc[3]
            +1.22242700e-13 * tc[4];
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
        /*species 2: O */
        species[2] =
            +2.92175791e+04 * invT
            -3.21491786e+00
            -2.56942078e+00 * tc[0]
            +4.29870569e-05 * tc[1]
            -6.99140982e-09 * tc[2]
            +8.34814992e-13 * tc[3]
            -6.14168455e-17 * tc[4];
        /*species 3: O2 */
        species[3] =
            -1.08845772e+03 * invT
            -3.17069345e+00
            -3.28253784e+00 * tc[0]
            -7.41543770e-04 * tc[1]
            +1.26327778e-07 * tc[2]
            -1.74558796e-11 * tc[3]
            +1.08358897e-15 * tc[4];
        /*species 4: OH */
        species[4] =
            +3.85865700e+03 * invT
            -2.38380843e+00
            -3.09288767e+00 * tc[0]
            -2.74214858e-04 * tc[1]
            -2.10842047e-08 * tc[2]
            +7.32884630e-12 * tc[3]
            -5.87061880e-16 * tc[4];
        /*species 5: H2O */
        species[5] =
            -3.00042971e+04 * invT
            -2.93277761e+00
            -3.03399249e+00 * tc[0]
            -1.08845902e-03 * tc[1]
            +2.73454197e-08 * tc[2]
            +8.08683225e-12 * tc[3]
            -8.41004960e-16 * tc[4];
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
            -1.78617877e+04 * invT
            +2.48846230e-01
            -4.16500285e+00 * tc[0]
            -2.45415847e-03 * tc[1]
            +3.16898708e-07 * tc[2]
            -3.09321655e-11 * tc[3]
            +1.43954153e-15 * tc[4];
        /*species 8: N2 */
        species[8] =
            -9.22797700e+02 * invT
            -4.05388800e+00
            -2.92664000e+00 * tc[0]
            -7.43988400e-04 * tc[1]
            +9.47460000e-08 * tc[2]
            -8.41419833e-12 * tc[3]
            +3.37667550e-16 * tc[4];
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
        /*species 2: O */
        species[2] =
            +2.16826710e+00
            -3.27931884e-03 * tc[1]
            +6.64306396e-06 * tc[2]
            -6.12806624e-09 * tc[3]
            +2.11265971e-12 * tc[4];
        /*species 3: O2 */
        species[3] =
            +2.78245636e+00
            -2.99673416e-03 * tc[1]
            +9.84730201e-06 * tc[2]
            -9.68129509e-09 * tc[3]
            +3.24372837e-12 * tc[4];
        /*species 4: OH */
        species[4] =
            +2.99201543e+00
            -2.40131752e-03 * tc[1]
            +4.61793841e-06 * tc[2]
            -3.88113333e-09 * tc[3]
            +1.36411470e-12 * tc[4];
        /*species 5: H2O */
        species[5] =
            +3.19864056e+00
            -2.03643410e-03 * tc[1]
            +6.52040211e-06 * tc[2]
            -5.48797062e-09 * tc[3]
            +1.77197817e-12 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +3.30179801e+00
            -4.74912051e-03 * tc[1]
            +2.11582891e-05 * tc[2]
            -2.42763894e-08 * tc[3]
            +9.29225124e-12 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +3.27611269e+00
            -5.42822417e-04 * tc[1]
            +1.67335701e-05 * tc[2]
            -2.15770813e-08 * tc[3]
            +8.62454363e-12 * tc[4];
        /*species 8: N2 */
        species[8] =
            +2.29867700e+00
            +1.40824040e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485400e-12 * tc[4];
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
        /*species 2: O */
        species[2] =
            +1.56942078e+00
            -8.59741137e-05 * tc[1]
            +4.19484589e-08 * tc[2]
            -1.00177799e-11 * tc[3]
            +1.22833691e-15 * tc[4];
        /*species 3: O2 */
        species[3] =
            +2.28253784e+00
            +1.48308754e-03 * tc[1]
            -7.57966669e-07 * tc[2]
            +2.09470555e-10 * tc[3]
            -2.16717794e-14 * tc[4];
        /*species 4: OH */
        species[4] =
            +2.09288767e+00
            +5.48429716e-04 * tc[1]
            +1.26505228e-07 * tc[2]
            -8.79461556e-11 * tc[3]
            +1.17412376e-14 * tc[4];
        /*species 5: H2O */
        species[5] =
            +2.03399249e+00
            +2.17691804e-03 * tc[1]
            -1.64072518e-07 * tc[2]
            -9.70419870e-11 * tc[3]
            +1.68200992e-14 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +3.01721090e+00
            +2.23982013e-03 * tc[1]
            -6.33658150e-07 * tc[2]
            +1.14246370e-10 * tc[3]
            -1.07908535e-14 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +3.16500285e+00
            +4.90831694e-03 * tc[1]
            -1.90139225e-06 * tc[2]
            +3.71185986e-10 * tc[3]
            -2.87908305e-14 * tc[4];
        /*species 8: N2 */
        species[8] =
            +1.92664000e+00
            +1.48797680e-03 * tc[1]
            -5.68476000e-07 * tc[2]
            +1.00970380e-10 * tc[3]
            -6.75335100e-15 * tc[4];
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
        /*species 2: O */
        species[2] =
            +3.16826710e+00
            -3.27931884e-03 * tc[1]
            +6.64306396e-06 * tc[2]
            -6.12806624e-09 * tc[3]
            +2.11265971e-12 * tc[4];
        /*species 3: O2 */
        species[3] =
            +3.78245636e+00
            -2.99673416e-03 * tc[1]
            +9.84730201e-06 * tc[2]
            -9.68129509e-09 * tc[3]
            +3.24372837e-12 * tc[4];
        /*species 4: OH */
        species[4] =
            +3.99201543e+00
            -2.40131752e-03 * tc[1]
            +4.61793841e-06 * tc[2]
            -3.88113333e-09 * tc[3]
            +1.36411470e-12 * tc[4];
        /*species 5: H2O */
        species[5] =
            +4.19864056e+00
            -2.03643410e-03 * tc[1]
            +6.52040211e-06 * tc[2]
            -5.48797062e-09 * tc[3]
            +1.77197817e-12 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +4.30179801e+00
            -4.74912051e-03 * tc[1]
            +2.11582891e-05 * tc[2]
            -2.42763894e-08 * tc[3]
            +9.29225124e-12 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +4.27611269e+00
            -5.42822417e-04 * tc[1]
            +1.67335701e-05 * tc[2]
            -2.15770813e-08 * tc[3]
            +8.62454363e-12 * tc[4];
        /*species 8: N2 */
        species[8] =
            +3.29867700e+00
            +1.40824040e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485400e-12 * tc[4];
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
        /*species 2: O */
        species[2] =
            +2.56942078e+00
            -8.59741137e-05 * tc[1]
            +4.19484589e-08 * tc[2]
            -1.00177799e-11 * tc[3]
            +1.22833691e-15 * tc[4];
        /*species 3: O2 */
        species[3] =
            +3.28253784e+00
            +1.48308754e-03 * tc[1]
            -7.57966669e-07 * tc[2]
            +2.09470555e-10 * tc[3]
            -2.16717794e-14 * tc[4];
        /*species 4: OH */
        species[4] =
            +3.09288767e+00
            +5.48429716e-04 * tc[1]
            +1.26505228e-07 * tc[2]
            -8.79461556e-11 * tc[3]
            +1.17412376e-14 * tc[4];
        /*species 5: H2O */
        species[5] =
            +3.03399249e+00
            +2.17691804e-03 * tc[1]
            -1.64072518e-07 * tc[2]
            -9.70419870e-11 * tc[3]
            +1.68200992e-14 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +4.01721090e+00
            +2.23982013e-03 * tc[1]
            -6.33658150e-07 * tc[2]
            +1.14246370e-10 * tc[3]
            -1.07908535e-14 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +4.16500285e+00
            +4.90831694e-03 * tc[1]
            -1.90139225e-06 * tc[2]
            +3.71185986e-10 * tc[3]
            -2.87908305e-14 * tc[4];
        /*species 8: N2 */
        species[8] =
            +2.92664000e+00
            +1.48797680e-03 * tc[1]
            -5.68476000e-07 * tc[2]
            +1.00970380e-10 * tc[3]
            -6.75335100e-15 * tc[4];
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
        /*species 2: O */
        species[2] =
            +2.16826710e+00
            -1.63965942e-03 * tc[1]
            +2.21435465e-06 * tc[2]
            -1.53201656e-09 * tc[3]
            +4.22531942e-13 * tc[4]
            +2.91222592e+04 * invT;
        /*species 3: O2 */
        species[3] =
            +2.78245636e+00
            -1.49836708e-03 * tc[1]
            +3.28243400e-06 * tc[2]
            -2.42032377e-09 * tc[3]
            +6.48745674e-13 * tc[4]
            -1.06394356e+03 * invT;
        /*species 4: OH */
        species[4] =
            +2.99201543e+00
            -1.20065876e-03 * tc[1]
            +1.53931280e-06 * tc[2]
            -9.70283332e-10 * tc[3]
            +2.72822940e-13 * tc[4]
            +3.61508056e+03 * invT;
        /*species 5: H2O */
        species[5] =
            +3.19864056e+00
            -1.01821705e-03 * tc[1]
            +2.17346737e-06 * tc[2]
            -1.37199266e-09 * tc[3]
            +3.54395634e-13 * tc[4]
            -3.02937267e+04 * invT;
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
            +3.27611269e+00
            -2.71411208e-04 * tc[1]
            +5.57785670e-06 * tc[2]
            -5.39427032e-09 * tc[3]
            +1.72490873e-12 * tc[4]
            -1.77025821e+04 * invT;
        /*species 8: N2 */
        species[8] =
            +2.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 * invT;
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
        /*species 2: O */
        species[2] =
            +1.56942078e+00
            -4.29870569e-05 * tc[1]
            +1.39828196e-08 * tc[2]
            -2.50444497e-12 * tc[3]
            +2.45667382e-16 * tc[4]
            +2.92175791e+04 * invT;
        /*species 3: O2 */
        species[3] =
            +2.28253784e+00
            +7.41543770e-04 * tc[1]
            -2.52655556e-07 * tc[2]
            +5.23676387e-11 * tc[3]
            -4.33435588e-15 * tc[4]
            -1.08845772e+03 * invT;
        /*species 4: OH */
        species[4] =
            +2.09288767e+00
            +2.74214858e-04 * tc[1]
            +4.21684093e-08 * tc[2]
            -2.19865389e-11 * tc[3]
            +2.34824752e-15 * tc[4]
            +3.85865700e+03 * invT;
        /*species 5: H2O */
        species[5] =
            +2.03399249e+00
            +1.08845902e-03 * tc[1]
            -5.46908393e-08 * tc[2]
            -2.42604967e-11 * tc[3]
            +3.36401984e-15 * tc[4]
            -3.00042971e+04 * invT;
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
            +3.16500285e+00
            +2.45415847e-03 * tc[1]
            -6.33797417e-07 * tc[2]
            +9.27964965e-11 * tc[3]
            -5.75816610e-15 * tc[4]
            -1.78617877e+04 * invT;
        /*species 8: N2 */
        species[8] =
            +1.92664000e+00
            +7.43988400e-04 * tc[1]
            -1.89492000e-07 * tc[2]
            +2.52425950e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
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
        /*species 2: O */
        species[2] =
            +3.16826710e+00
            -1.63965942e-03 * tc[1]
            +2.21435465e-06 * tc[2]
            -1.53201656e-09 * tc[3]
            +4.22531942e-13 * tc[4]
            +2.91222592e+04 * invT;
        /*species 3: O2 */
        species[3] =
            +3.78245636e+00
            -1.49836708e-03 * tc[1]
            +3.28243400e-06 * tc[2]
            -2.42032377e-09 * tc[3]
            +6.48745674e-13 * tc[4]
            -1.06394356e+03 * invT;
        /*species 4: OH */
        species[4] =
            +3.99201543e+00
            -1.20065876e-03 * tc[1]
            +1.53931280e-06 * tc[2]
            -9.70283332e-10 * tc[3]
            +2.72822940e-13 * tc[4]
            +3.61508056e+03 * invT;
        /*species 5: H2O */
        species[5] =
            +4.19864056e+00
            -1.01821705e-03 * tc[1]
            +2.17346737e-06 * tc[2]
            -1.37199266e-09 * tc[3]
            +3.54395634e-13 * tc[4]
            -3.02937267e+04 * invT;
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
            +4.27611269e+00
            -2.71411208e-04 * tc[1]
            +5.57785670e-06 * tc[2]
            -5.39427032e-09 * tc[3]
            +1.72490873e-12 * tc[4]
            -1.77025821e+04 * invT;
        /*species 8: N2 */
        species[8] =
            +3.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 * invT;
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
        /*species 2: O */
        species[2] =
            +2.56942078e+00
            -4.29870569e-05 * tc[1]
            +1.39828196e-08 * tc[2]
            -2.50444497e-12 * tc[3]
            +2.45667382e-16 * tc[4]
            +2.92175791e+04 * invT;
        /*species 3: O2 */
        species[3] =
            +3.28253784e+00
            +7.41543770e-04 * tc[1]
            -2.52655556e-07 * tc[2]
            +5.23676387e-11 * tc[3]
            -4.33435588e-15 * tc[4]
            -1.08845772e+03 * invT;
        /*species 4: OH */
        species[4] =
            +3.09288767e+00
            +2.74214858e-04 * tc[1]
            +4.21684093e-08 * tc[2]
            -2.19865389e-11 * tc[3]
            +2.34824752e-15 * tc[4]
            +3.85865700e+03 * invT;
        /*species 5: H2O */
        species[5] =
            +3.03399249e+00
            +1.08845902e-03 * tc[1]
            -5.46908393e-08 * tc[2]
            -2.42604967e-11 * tc[3]
            +3.36401984e-15 * tc[4]
            -3.00042971e+04 * invT;
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
            +4.16500285e+00
            +2.45415847e-03 * tc[1]
            -6.33797417e-07 * tc[2]
            +9.27964965e-11 * tc[3]
            -5.75816610e-15 * tc[4]
            -1.78617877e+04 * invT;
        /*species 8: N2 */
        species[8] =
            +2.92664000e+00
            +7.43988400e-04 * tc[1]
            -1.89492000e-07 * tc[2]
            +2.52425950e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
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
        /*species 2: O */
        species[2] =
            +3.16826710e+00 * tc[0]
            -3.27931884e-03 * tc[1]
            +3.32153198e-06 * tc[2]
            -2.04268875e-09 * tc[3]
            +5.28164927e-13 * tc[4]
            +2.05193346e+00 ;
        /*species 3: O2 */
        species[3] =
            +3.78245636e+00 * tc[0]
            -2.99673416e-03 * tc[1]
            +4.92365101e-06 * tc[2]
            -3.22709836e-09 * tc[3]
            +8.10932092e-13 * tc[4]
            +3.65767573e+00 ;
        /*species 4: OH */
        species[4] =
            +3.99201543e+00 * tc[0]
            -2.40131752e-03 * tc[1]
            +2.30896920e-06 * tc[2]
            -1.29371111e-09 * tc[3]
            +3.41028675e-13 * tc[4]
            -1.03925458e-01 ;
        /*species 5: H2O */
        species[5] =
            +4.19864056e+00 * tc[0]
            -2.03643410e-03 * tc[1]
            +3.26020105e-06 * tc[2]
            -1.82932354e-09 * tc[3]
            +4.42994543e-13 * tc[4]
            -8.49032208e-01 ;
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
            +4.27611269e+00 * tc[0]
            -5.42822417e-04 * tc[1]
            +8.36678505e-06 * tc[2]
            -7.19236043e-09 * tc[3]
            +2.15613591e-12 * tc[4]
            +3.43505074e+00 ;
        /*species 8: N2 */
        species[8] =
            +3.29867700e+00 * tc[0]
            +1.40824040e-03 * tc[1]
            -1.98161100e-06 * tc[2]
            +1.88050500e-09 * tc[3]
            -6.11213500e-13 * tc[4]
            +3.95037200e+00 ;
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
        /*species 2: O */
        species[2] =
            +2.56942078e+00 * tc[0]
            -8.59741137e-05 * tc[1]
            +2.09742295e-08 * tc[2]
            -3.33925997e-12 * tc[3]
            +3.07084227e-16 * tc[4]
            +4.78433864e+00 ;
        /*species 3: O2 */
        species[3] =
            +3.28253784e+00 * tc[0]
            +1.48308754e-03 * tc[1]
            -3.78983334e-07 * tc[2]
            +6.98235183e-11 * tc[3]
            -5.41794485e-15 * tc[4]
            +5.45323129e+00 ;
        /*species 4: OH */
        species[4] =
            +3.09288767e+00 * tc[0]
            +5.48429716e-04 * tc[1]
            +6.32526140e-08 * tc[2]
            -2.93153852e-11 * tc[3]
            +2.93530940e-15 * tc[4]
            +4.47669610e+00 ;
        /*species 5: H2O */
        species[5] =
            +3.03399249e+00 * tc[0]
            +2.17691804e-03 * tc[1]
            -8.20362590e-08 * tc[2]
            -3.23473290e-11 * tc[3]
            +4.20502480e-15 * tc[4]
            +4.96677010e+00 ;
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
            +4.16500285e+00 * tc[0]
            +4.90831694e-03 * tc[1]
            -9.50696125e-07 * tc[2]
            +1.23728662e-10 * tc[3]
            -7.19770763e-15 * tc[4]
            +2.91615662e+00 ;
        /*species 8: N2 */
        species[8] =
            +2.92664000e+00 * tc[0]
            +1.48797680e-03 * tc[1]
            -2.84238000e-07 * tc[2]
            +3.36567933e-11 * tc[3]
            -1.68833775e-15 * tc[4]
            +5.98052800e+00 ;
    }
    return;
}


/*save atomic weights into array */
void atomicWeight(double *  awt)
{
    awt[0] = 15.999400; /*O */
    awt[1] = 1.007970; /*H */
    awt[2] = 14.006700; /*N */

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

    double   EPS[9];
    double   SIG[9];
    double    wt[9];
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

    /*species 2: O */
    Tci[2] = 1.316 * EPS[2] ; 
    ai[2] = (5.55 * pow(avogadro,2.0) * EPS[2]*boltzmann * pow(1e-8*SIG[2],3.0) ) / (pow(wt[2],2.0)); 
    bi[2] = 0.855 * avogadro * pow(1e-8*SIG[2],3.0) / (wt[2]); 
    acentric_i[2] = 0.0 ;

    /*species 3: O2 */
    /*Imported from NIST */
    Tci[3] = 154.581000 ; 
    ai[3] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[3],2.0) / (pow(31.998800,2.0) * 50.430466); 
    bi[3] = 0.08664 * Rcst * Tci[3] / (31.998800 * 50.430466); 
    acentric_i[3] = 0.022200 ;

    /*species 4: OH */
    Tci[4] = 1.316 * EPS[4] ; 
    ai[4] = (5.55 * pow(avogadro,2.0) * EPS[4]*boltzmann * pow(1e-8*SIG[4],3.0) ) / (pow(wt[4],2.0)); 
    bi[4] = 0.855 * avogadro * pow(1e-8*SIG[4],3.0) / (wt[4]); 
    acentric_i[4] = 0.0 ;

    /*species 5: H2O */
    /*Imported from NIST */
    Tci[5] = 647.096000 ; 
    ai[5] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[5],2.0) / (pow(18.015340,2.0) * 220.640000); 
    bi[5] = 0.08664 * Rcst * Tci[5] / (18.015340 * 220.640000); 
    acentric_i[5] = 0.344300 ;

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

    return;
}


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
void egtransetPATM(double* PATM) {
    *PATM =   0.1013250000000000E+07;}


/*the molecular weights in g/mol */
void egtransetWT(double* WT ) {
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
void egtransetEPS(double* EPS ) {
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
void egtransetSIG(double* SIG ) {
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
void egtransetDIP(double* DIP ) {
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
void egtransetPOL(double* POL ) {
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
void egtransetZROT(double* ZROT ) {
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
void egtransetCOFETA(double* COFETA) {
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
void egtransetCOFLAM(double* COFLAM) {
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
void egtransetCOFD(double* COFD) {
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
    KTDIF[0] = 1;
    KTDIF[1] = 2;
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

/* Replace this routine with the one generated by the Gauss Jordan solver of DW */
AMREX_GPU_HOST_DEVICE void sgjsolve(double* A, double* x, double* b) {
    amrex::Abort("sgjsolve not implemented, choose a different solver ");
}

/* Replace this routine with the one generated by the Gauss Jordan solver of DW */
AMREX_GPU_HOST_DEVICE void sgjsolve_simplified(double* A, double* x, double* b) {
    amrex::Abort("sgjsolve_simplified not implemented, choose a different solver ");
}

