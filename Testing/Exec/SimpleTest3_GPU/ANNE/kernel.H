#ifndef _kernel_H_
#define _kernel_H_
#include <AMReX_GpuDevice.H>
#include "chemistry_file.H"

#define NUM_SPECIES 21
#define NUM_REACTIONS 84

using namespace amrex;
namespace mythermo
{
    AMREX_GPU_DEVICE_MANAGED amrex::Real fwd_A[84], fwd_beta[84], fwd_Ea[84];
    AMREX_GPU_DEVICE_MANAGED amrex::Real activation_units[84], prefactor_units[84], phase_units[84];
    AMREX_GPU_DEVICE_MANAGED amrex::GpuArray<amrex::GpuArray<int,5>, 84> kiv, nuv;
    AMREX_GPU_DEVICE_MANAGED int ki2D[420], nu2D[420];
}

using namespace mythermo;

static AMREX_GPU_DEVICE_MANAGED amrex::Real imw[21] = {
    1.0 / 2.015940,  /*H2 */
    1.0 / 1.007970,  /*H */
    1.0 / 15.999400,  /*O */
    1.0 / 31.998800,  /*O2 */
    1.0 / 17.007370,  /*OH */
    1.0 / 18.015340,  /*H2O */
    1.0 / 33.006770,  /*HO2 */
    1.0 / 14.027090,  /*CH2 */
    1.0 / 14.027090,  /*CH2(S) */
    1.0 / 15.035060,  /*CH3 */
    1.0 / 16.043030,  /*CH4 */
    1.0 / 28.010550,  /*CO */
    1.0 / 44.009950,  /*CO2 */
    1.0 / 29.018520,  /*HCO */
    1.0 / 30.026490,  /*CH2O */
    1.0 / 31.034460,  /*CH3O */
    1.0 / 28.054180,  /*C2H4 */
    1.0 / 29.062150,  /*C2H5 */
    1.0 / 30.070120,  /*C2H6 */
    1.0 / 28.013400,  /*N2 */
    1.0 / 39.948000};  /*AR */


void kinit()
{
    // (0):  O + H + M <=> OH + M
    kiv[8] = {2,1,4};
    nuv[8] = {-1,-1,1};
    fwd_A[8]     = 5e+17;
    fwd_beta[8]  = -1;
    fwd_Ea[8]    = 0;
    prefactor_units[8]  = 1.0000000000000002e-12;
    activation_units[8] = 0.50321666580471969;
    phase_units[8]      = 1e-12;

    // (1):  O + H2 <=> H + OH
    kiv[14] = {2,0,1,4};
    nuv[14] = {-1,-1,1,1};
    fwd_A[14]     = 50000;
    fwd_beta[14]  = 2.6699999999999999;
    fwd_Ea[14]    = 6290;
    prefactor_units[14]  = 1.0000000000000002e-06;
    activation_units[14] = 0.50321666580471969;
    phase_units[14]      = 1e-12;

    // (2):  O + HO2 <=> OH + O2
    kiv[15] = {2,6,4,3};
    nuv[15] = {-1,-1,1,1};
    fwd_A[15]     = 20000000000000;
    fwd_beta[15]  = 0;
    fwd_Ea[15]    = 0;
    prefactor_units[15]  = 1.0000000000000002e-06;
    activation_units[15] = 0.50321666580471969;
    phase_units[15]      = 1e-12;

    // (3):  O + CH2 <=> H + HCO
    kiv[16] = {2,7,1,13};
    nuv[16] = {-1,-1,1,1};
    fwd_A[16]     = 80000000000000;
    fwd_beta[16]  = 0;
    fwd_Ea[16]    = 0;
    prefactor_units[16]  = 1.0000000000000002e-06;
    activation_units[16] = 0.50321666580471969;
    phase_units[16]      = 1e-12;

    // (4):  O + CH2(S) <=> H + HCO
    kiv[17] = {2,8,1,13};
    nuv[17] = {-1,-1,1,1};
    fwd_A[17]     = 15000000000000;
    fwd_beta[17]  = 0;
    fwd_Ea[17]    = 0;
    prefactor_units[17]  = 1.0000000000000002e-06;
    activation_units[17] = 0.50321666580471969;
    phase_units[17]      = 1e-12;

    // (5):  O + CH3 <=> H + CH2O
    kiv[18] = {2,9,1,14};
    nuv[18] = {-1,-1,1,1};
    fwd_A[18]     = 84300000000000;
    fwd_beta[18]  = 0;
    fwd_Ea[18]    = 0;
    prefactor_units[18]  = 1.0000000000000002e-06;
    activation_units[18] = 0.50321666580471969;
    phase_units[18]      = 1e-12;

    // (6):  O + CH4 <=> OH + CH3
    kiv[19] = {2,10,4,9};
    nuv[19] = {-1,-1,1,1};
    fwd_A[19]     = 1020000000;
    fwd_beta[19]  = 1.5;
    fwd_Ea[19]    = 8600;
    prefactor_units[19]  = 1.0000000000000002e-06;
    activation_units[19] = 0.50321666580471969;
    phase_units[19]      = 1e-12;

    // (7):  O + CO + M <=> CO2 + M
    kiv[9] = {2,11,12};
    nuv[9] = {-1,-1,1};
    fwd_A[9]     = 602000000000000;
    fwd_beta[9]  = 0;
    fwd_Ea[9]    = 3000;
    prefactor_units[9]  = 1.0000000000000002e-12;
    activation_units[9] = 0.50321666580471969;
    phase_units[9]      = 1e-12;

    // (8):  O + HCO <=> OH + CO
    kiv[20] = {2,13,4,11};
    nuv[20] = {-1,-1,1,1};
    fwd_A[20]     = 30000000000000;
    fwd_beta[20]  = 0;
    fwd_Ea[20]    = 0;
    prefactor_units[20]  = 1.0000000000000002e-06;
    activation_units[20] = 0.50321666580471969;
    phase_units[20]      = 1e-12;

    // (9):  O + HCO <=> H + CO2
    kiv[21] = {2,13,1,12};
    nuv[21] = {-1,-1,1,1};
    fwd_A[21]     = 30000000000000;
    fwd_beta[21]  = 0;
    fwd_Ea[21]    = 0;
    prefactor_units[21]  = 1.0000000000000002e-06;
    activation_units[21] = 0.50321666580471969;
    phase_units[21]      = 1e-12;

    // (10):  O + CH2O <=> OH + HCO
    kiv[22] = {2,14,4,13};
    nuv[22] = {-1,-1,1,1};
    fwd_A[22]     = 39000000000000;
    fwd_beta[22]  = 0;
    fwd_Ea[22]    = 3540;
    prefactor_units[22]  = 1.0000000000000002e-06;
    activation_units[22] = 0.50321666580471969;
    phase_units[22]      = 1e-12;

    // (11):  O + C2H4 <=> CH3 + HCO
    kiv[23] = {2,16,9,13};
    nuv[23] = {-1,-1,1,1};
    fwd_A[23]     = 19200000;
    fwd_beta[23]  = 1.8300000000000001;
    fwd_Ea[23]    = 220;
    prefactor_units[23]  = 1.0000000000000002e-06;
    activation_units[23] = 0.50321666580471969;
    phase_units[23]      = 1e-12;

    // (12):  O + C2H5 <=> CH3 + CH2O
    kiv[24] = {2,17,9,14};
    nuv[24] = {-1,-1,1,1};
    fwd_A[24]     = 132000000000000;
    fwd_beta[24]  = 0;
    fwd_Ea[24]    = 0;
    prefactor_units[24]  = 1.0000000000000002e-06;
    activation_units[24] = 0.50321666580471969;
    phase_units[24]      = 1e-12;

    // (13):  O + C2H6 <=> OH + C2H5
    kiv[25] = {2,18,4,17};
    nuv[25] = {-1,-1,1,1};
    fwd_A[25]     = 89800000;
    fwd_beta[25]  = 1.9199999999999999;
    fwd_Ea[25]    = 5690;
    prefactor_units[25]  = 1.0000000000000002e-06;
    activation_units[25] = 0.50321666580471969;
    phase_units[25]      = 1e-12;

    // (14):  O2 + CO <=> O + CO2
    kiv[26] = {3,11,2,12};
    nuv[26] = {-1,-1,1,1};
    fwd_A[26]     = 2500000000000;
    fwd_beta[26]  = 0;
    fwd_Ea[26]    = 47800;
    prefactor_units[26]  = 1.0000000000000002e-06;
    activation_units[26] = 0.50321666580471969;
    phase_units[26]      = 1e-12;

    // (15):  O2 + CH2O <=> HO2 + HCO
    kiv[27] = {3,14,6,13};
    nuv[27] = {-1,-1,1,1};
    fwd_A[27]     = 100000000000000;
    fwd_beta[27]  = 0;
    fwd_Ea[27]    = 40000;
    prefactor_units[27]  = 1.0000000000000002e-06;
    activation_units[27] = 0.50321666580471969;
    phase_units[27]      = 1e-12;

    // (16):  H + O2 + M <=> HO2 + M
    kiv[10] = {1,3,6};
    nuv[10] = {-1,-1,1};
    fwd_A[10]     = 2.8e+18;
    fwd_beta[10]  = -0.85999999999999999;
    fwd_Ea[10]    = 0;
    prefactor_units[10]  = 1.0000000000000002e-12;
    activation_units[10] = 0.50321666580471969;
    phase_units[10]      = 1e-12;

    // (17):  H + 2 O2 <=> HO2 + O2
    kiv[28] = {1,3,6,3};
    nuv[28] = {-1,-2,1,1};
    fwd_A[28]     = 3e+20;
    fwd_beta[28]  = -1.72;
    fwd_Ea[28]    = 0;
    prefactor_units[28]  = 1.0000000000000002e-12;
    activation_units[28] = 0.50321666580471969;
    phase_units[28]      = 1e-18;

    // (18):  H + O2 + H2O <=> HO2 + H2O
    kiv[29] = {1,3,5,6,5};
    nuv[29] = {-1,-1,-1,1,1};
    fwd_A[29]     = 9.38e+18;
    fwd_beta[29]  = -0.76000000000000001;
    fwd_Ea[29]    = 0;
    prefactor_units[29]  = 1.0000000000000002e-12;
    activation_units[29] = 0.50321666580471969;
    phase_units[29]      = 1e-18;

    // (19):  H + O2 + N2 <=> HO2 + N2
    kiv[30] = {1,3,19,6,19};
    nuv[30] = {-1,-1,-1,1,1};
    fwd_A[30]     = 3.75e+20;
    fwd_beta[30]  = -1.72;
    fwd_Ea[30]    = 0;
    prefactor_units[30]  = 1.0000000000000002e-12;
    activation_units[30] = 0.50321666580471969;
    phase_units[30]      = 1e-18;

    // (20):  H + O2 + AR <=> HO2 + AR
    kiv[31] = {1,3,20,6,20};
    nuv[31] = {-1,-1,-1,1,1};
    fwd_A[31]     = 7e+17;
    fwd_beta[31]  = -0.80000000000000004;
    fwd_Ea[31]    = 0;
    prefactor_units[31]  = 1.0000000000000002e-12;
    activation_units[31] = 0.50321666580471969;
    phase_units[31]      = 1e-18;

    // (21):  H + O2 <=> O + OH
    kiv[32] = {1,3,2,4};
    nuv[32] = {-1,-1,1,1};
    fwd_A[32]     = 83000000000000;
    fwd_beta[32]  = 0;
    fwd_Ea[32]    = 14413;
    prefactor_units[32]  = 1.0000000000000002e-06;
    activation_units[32] = 0.50321666580471969;
    phase_units[32]      = 1e-12;

    // (22):  2 H + M <=> H2 + M
    kiv[11] = {1,0};
    nuv[11] = {-2,1};
    fwd_A[11]     = 1e+18;
    fwd_beta[11]  = -1;
    fwd_Ea[11]    = 0;
    prefactor_units[11]  = 1.0000000000000002e-12;
    activation_units[11] = 0.50321666580471969;
    phase_units[11]      = 1e-12;

    // (23):  2 H + H2 <=> 2 H2
    kiv[33] = {1,0,0};
    nuv[33] = {-2,-1,2};
    fwd_A[33]     = 90000000000000000;
    fwd_beta[33]  = -0.59999999999999998;
    fwd_Ea[33]    = 0;
    prefactor_units[33]  = 1.0000000000000002e-12;
    activation_units[33] = 0.50321666580471969;
    phase_units[33]      = 1e-18;

    // (24):  2 H + H2O <=> H2 + H2O
    kiv[34] = {1,5,0,5};
    nuv[34] = {-2,-1,1,1};
    fwd_A[34]     = 6e+19;
    fwd_beta[34]  = -1.25;
    fwd_Ea[34]    = 0;
    prefactor_units[34]  = 1.0000000000000002e-12;
    activation_units[34] = 0.50321666580471969;
    phase_units[34]      = 1e-18;

    // (25):  2 H + CO2 <=> H2 + CO2
    kiv[35] = {1,12,0,12};
    nuv[35] = {-2,-1,1,1};
    fwd_A[35]     = 5.5e+20;
    fwd_beta[35]  = -2;
    fwd_Ea[35]    = 0;
    prefactor_units[35]  = 1.0000000000000002e-12;
    activation_units[35] = 0.50321666580471969;
    phase_units[35]      = 1e-18;

    // (26):  H + OH + M <=> H2O + M
    kiv[12] = {1,4,5};
    nuv[12] = {-1,-1,1};
    fwd_A[12]     = 2.2e+22;
    fwd_beta[12]  = -2;
    fwd_Ea[12]    = 0;
    prefactor_units[12]  = 1.0000000000000002e-12;
    activation_units[12] = 0.50321666580471969;
    phase_units[12]      = 1e-12;

    // (27):  H + HO2 <=> O2 + H2
    kiv[36] = {1,6,3,0};
    nuv[36] = {-1,-1,1,1};
    fwd_A[36]     = 28000000000000;
    fwd_beta[36]  = 0;
    fwd_Ea[36]    = 1068;
    prefactor_units[36]  = 1.0000000000000002e-06;
    activation_units[36] = 0.50321666580471969;
    phase_units[36]      = 1e-12;

    // (28):  H + HO2 <=> 2 OH
    kiv[37] = {1,6,4};
    nuv[37] = {-1,-1,2};
    fwd_A[37]     = 134000000000000;
    fwd_beta[37]  = 0;
    fwd_Ea[37]    = 635;
    prefactor_units[37]  = 1.0000000000000002e-06;
    activation_units[37] = 0.50321666580471969;
    phase_units[37]      = 1e-12;

    // (29):  H + CH2 (+M) <=> CH3 (+M)
    kiv[0] = {1,7,9};
    nuv[0] = {-1,-1,1};
    fwd_A[0]     = 25000000000000000;
    fwd_beta[0]  = -0.80000000000000004;
    fwd_Ea[0]    = 0;
    prefactor_units[0]  = 1.0000000000000002e-06;
    activation_units[0] = 0.50321666580471969;
    phase_units[0]      = 1e-12;

    // (30):  H + CH3 (+M) <=> CH4 (+M)
    kiv[1] = {1,9,10};
    nuv[1] = {-1,-1,1};
    fwd_A[1]     = 12700000000000000;
    fwd_beta[1]  = -0.63;
    fwd_Ea[1]    = 383;
    prefactor_units[1]  = 1.0000000000000002e-06;
    activation_units[1] = 0.50321666580471969;
    phase_units[1]      = 1e-12;

    // (31):  H + CH4 <=> CH3 + H2
    kiv[38] = {1,10,9,0};
    nuv[38] = {-1,-1,1,1};
    fwd_A[38]     = 660000000;
    fwd_beta[38]  = 1.6200000000000001;
    fwd_Ea[38]    = 10840;
    prefactor_units[38]  = 1.0000000000000002e-06;
    activation_units[38] = 0.50321666580471969;
    phase_units[38]      = 1e-12;

    // (32):  H + HCO (+M) <=> CH2O (+M)
    kiv[2] = {1,13,14};
    nuv[2] = {-1,-1,1};
    fwd_A[2]     = 1090000000000;
    fwd_beta[2]  = 0.47999999999999998;
    fwd_Ea[2]    = -260;
    prefactor_units[2]  = 1.0000000000000002e-06;
    activation_units[2] = 0.50321666580471969;
    phase_units[2]      = 1e-12;

    // (33):  H + HCO <=> H2 + CO
    kiv[39] = {1,13,0,11};
    nuv[39] = {-1,-1,1,1};
    fwd_A[39]     = 73400000000000;
    fwd_beta[39]  = 0;
    fwd_Ea[39]    = 0;
    prefactor_units[39]  = 1.0000000000000002e-06;
    activation_units[39] = 0.50321666580471969;
    phase_units[39]      = 1e-12;

    // (34):  H + CH2O (+M) <=> CH3O (+M)
    kiv[3] = {1,14,15};
    nuv[3] = {-1,-1,1};
    fwd_A[3]     = 540000000000;
    fwd_beta[3]  = 0.45400000000000001;
    fwd_Ea[3]    = 2600;
    prefactor_units[3]  = 1.0000000000000002e-06;
    activation_units[3] = 0.50321666580471969;
    phase_units[3]      = 1e-12;

    // (35):  H + CH2O <=> HCO + H2
    kiv[40] = {1,14,13,0};
    nuv[40] = {-1,-1,1,1};
    fwd_A[40]     = 23000000000;
    fwd_beta[40]  = 1.05;
    fwd_Ea[40]    = 3275;
    prefactor_units[40]  = 1.0000000000000002e-06;
    activation_units[40] = 0.50321666580471969;
    phase_units[40]      = 1e-12;

    // (36):  H + CH3O <=> OH + CH3
    kiv[41] = {1,15,4,9};
    nuv[41] = {-1,-1,1,1};
    fwd_A[41]     = 32000000000000;
    fwd_beta[41]  = 0;
    fwd_Ea[41]    = 0;
    prefactor_units[41]  = 1.0000000000000002e-06;
    activation_units[41] = 0.50321666580471969;
    phase_units[41]      = 1e-12;

    // (37):  H + C2H4 (+M) <=> C2H5 (+M)
    kiv[4] = {1,16,17};
    nuv[4] = {-1,-1,1};
    fwd_A[4]     = 1080000000000;
    fwd_beta[4]  = 0.45400000000000001;
    fwd_Ea[4]    = 1820;
    prefactor_units[4]  = 1.0000000000000002e-06;
    activation_units[4] = 0.50321666580471969;
    phase_units[4]      = 1e-12;

    // (38):  H + C2H5 (+M) <=> C2H6 (+M)
    kiv[5] = {1,17,18};
    nuv[5] = {-1,-1,1};
    fwd_A[5]     = 5.21e+17;
    fwd_beta[5]  = -0.98999999999999999;
    fwd_Ea[5]    = 1580;
    prefactor_units[5]  = 1.0000000000000002e-06;
    activation_units[5] = 0.50321666580471969;
    phase_units[5]      = 1e-12;

    // (39):  H + C2H6 <=> C2H5 + H2
    kiv[42] = {1,18,17,0};
    nuv[42] = {-1,-1,1,1};
    fwd_A[42]     = 115000000;
    fwd_beta[42]  = 1.8999999999999999;
    fwd_Ea[42]    = 7530;
    prefactor_units[42]  = 1.0000000000000002e-06;
    activation_units[42] = 0.50321666580471969;
    phase_units[42]      = 1e-12;

    // (40):  H2 + CO (+M) <=> CH2O (+M)
    kiv[6] = {0,11,14};
    nuv[6] = {-1,-1,1};
    fwd_A[6]     = 43000000;
    fwd_beta[6]  = 1.5;
    fwd_Ea[6]    = 79600;
    prefactor_units[6]  = 1.0000000000000002e-06;
    activation_units[6] = 0.50321666580471969;
    phase_units[6]      = 1e-12;

    // (41):  OH + H2 <=> H + H2O
    kiv[43] = {4,0,1,5};
    nuv[43] = {-1,-1,1,1};
    fwd_A[43]     = 216000000;
    fwd_beta[43]  = 1.51;
    fwd_Ea[43]    = 3430;
    prefactor_units[43]  = 1.0000000000000002e-06;
    activation_units[43] = 0.50321666580471969;
    phase_units[43]      = 1e-12;

    // (42):  2 OH <=> O + H2O
    kiv[44] = {4,2,5};
    nuv[44] = {-2,1,1};
    fwd_A[44]     = 35700;
    fwd_beta[44]  = 2.3999999999999999;
    fwd_Ea[44]    = -2110;
    prefactor_units[44]  = 1.0000000000000002e-06;
    activation_units[44] = 0.50321666580471969;
    phase_units[44]      = 1e-12;

    // (43):  OH + HO2 <=> O2 + H2O
    kiv[45] = {4,6,3,5};
    nuv[45] = {-1,-1,1,1};
    fwd_A[45]     = 29000000000000;
    fwd_beta[45]  = 0;
    fwd_Ea[45]    = -500;
    prefactor_units[45]  = 1.0000000000000002e-06;
    activation_units[45] = 0.50321666580471969;
    phase_units[45]      = 1e-12;

    // (44):  OH + CH2 <=> H + CH2O
    kiv[46] = {4,7,1,14};
    nuv[46] = {-1,-1,1,1};
    fwd_A[46]     = 20000000000000;
    fwd_beta[46]  = 0;
    fwd_Ea[46]    = 0;
    prefactor_units[46]  = 1.0000000000000002e-06;
    activation_units[46] = 0.50321666580471969;
    phase_units[46]      = 1e-12;

    // (45):  OH + CH2(S) <=> H + CH2O
    kiv[47] = {4,8,1,14};
    nuv[47] = {-1,-1,1,1};
    fwd_A[47]     = 30000000000000;
    fwd_beta[47]  = 0;
    fwd_Ea[47]    = 0;
    prefactor_units[47]  = 1.0000000000000002e-06;
    activation_units[47] = 0.50321666580471969;
    phase_units[47]      = 1e-12;

    // (46):  OH + CH3 <=> CH2 + H2O
    kiv[48] = {4,9,7,5};
    nuv[48] = {-1,-1,1,1};
    fwd_A[48]     = 56000000;
    fwd_beta[48]  = 1.6000000000000001;
    fwd_Ea[48]    = 5420;
    prefactor_units[48]  = 1.0000000000000002e-06;
    activation_units[48] = 0.50321666580471969;
    phase_units[48]      = 1e-12;

    // (47):  OH + CH3 <=> CH2(S) + H2O
    kiv[49] = {4,9,8,5};
    nuv[49] = {-1,-1,1,1};
    fwd_A[49]     = 25010000000000;
    fwd_beta[49]  = 0;
    fwd_Ea[49]    = 0;
    prefactor_units[49]  = 1.0000000000000002e-06;
    activation_units[49] = 0.50321666580471969;
    phase_units[49]      = 1e-12;

    // (48):  OH + CH4 <=> CH3 + H2O
    kiv[50] = {4,10,9,5};
    nuv[50] = {-1,-1,1,1};
    fwd_A[50]     = 100000000;
    fwd_beta[50]  = 1.6000000000000001;
    fwd_Ea[50]    = 3120;
    prefactor_units[50]  = 1.0000000000000002e-06;
    activation_units[50] = 0.50321666580471969;
    phase_units[50]      = 1e-12;

    // (49):  OH + CO <=> H + CO2
    kiv[51] = {4,11,1,12};
    nuv[51] = {-1,-1,1,1};
    fwd_A[51]     = 47600000;
    fwd_beta[51]  = 1.228;
    fwd_Ea[51]    = 70;
    prefactor_units[51]  = 1.0000000000000002e-06;
    activation_units[51] = 0.50321666580471969;
    phase_units[51]      = 1e-12;

    // (50):  OH + HCO <=> H2O + CO
    kiv[52] = {4,13,5,11};
    nuv[52] = {-1,-1,1,1};
    fwd_A[52]     = 50000000000000;
    fwd_beta[52]  = 0;
    fwd_Ea[52]    = 0;
    prefactor_units[52]  = 1.0000000000000002e-06;
    activation_units[52] = 0.50321666580471969;
    phase_units[52]      = 1e-12;

    // (51):  OH + CH2O <=> HCO + H2O
    kiv[53] = {4,14,13,5};
    nuv[53] = {-1,-1,1,1};
    fwd_A[53]     = 3430000000;
    fwd_beta[53]  = 1.1799999999999999;
    fwd_Ea[53]    = -447;
    prefactor_units[53]  = 1.0000000000000002e-06;
    activation_units[53] = 0.50321666580471969;
    phase_units[53]      = 1e-12;

    // (52):  OH + C2H6 <=> C2H5 + H2O
    kiv[54] = {4,18,17,5};
    nuv[54] = {-1,-1,1,1};
    fwd_A[54]     = 3540000;
    fwd_beta[54]  = 2.1200000000000001;
    fwd_Ea[54]    = 870;
    prefactor_units[54]  = 1.0000000000000002e-06;
    activation_units[54] = 0.50321666580471969;
    phase_units[54]      = 1e-12;

    // (53):  HO2 + CH2 <=> OH + CH2O
    kiv[55] = {6,7,4,14};
    nuv[55] = {-1,-1,1,1};
    fwd_A[55]     = 20000000000000;
    fwd_beta[55]  = 0;
    fwd_Ea[55]    = 0;
    prefactor_units[55]  = 1.0000000000000002e-06;
    activation_units[55] = 0.50321666580471969;
    phase_units[55]      = 1e-12;

    // (54):  HO2 + CH3 <=> O2 + CH4
    kiv[56] = {6,9,3,10};
    nuv[56] = {-1,-1,1,1};
    fwd_A[56]     = 1000000000000;
    fwd_beta[56]  = 0;
    fwd_Ea[56]    = 0;
    prefactor_units[56]  = 1.0000000000000002e-06;
    activation_units[56] = 0.50321666580471969;
    phase_units[56]      = 1e-12;

    // (55):  HO2 + CH3 <=> OH + CH3O
    kiv[57] = {6,9,4,15};
    nuv[57] = {-1,-1,1,1};
    fwd_A[57]     = 20000000000000;
    fwd_beta[57]  = 0;
    fwd_Ea[57]    = 0;
    prefactor_units[57]  = 1.0000000000000002e-06;
    activation_units[57] = 0.50321666580471969;
    phase_units[57]      = 1e-12;

    // (56):  HO2 + CO <=> OH + CO2
    kiv[58] = {6,11,4,12};
    nuv[58] = {-1,-1,1,1};
    fwd_A[58]     = 150000000000000;
    fwd_beta[58]  = 0;
    fwd_Ea[58]    = 23600;
    prefactor_units[58]  = 1.0000000000000002e-06;
    activation_units[58] = 0.50321666580471969;
    phase_units[58]      = 1e-12;

    // (57):  CH2 + O2 <=> OH + HCO
    kiv[59] = {7,3,4,13};
    nuv[59] = {-1,-1,1,1};
    fwd_A[59]     = 13200000000000;
    fwd_beta[59]  = 0;
    fwd_Ea[59]    = 1500;
    prefactor_units[59]  = 1.0000000000000002e-06;
    activation_units[59] = 0.50321666580471969;
    phase_units[59]      = 1e-12;

    // (58):  CH2 + H2 <=> H + CH3
    kiv[60] = {7,0,1,9};
    nuv[60] = {-1,-1,1,1};
    fwd_A[60]     = 500000;
    fwd_beta[60]  = 2;
    fwd_Ea[60]    = 7230;
    prefactor_units[60]  = 1.0000000000000002e-06;
    activation_units[60] = 0.50321666580471969;
    phase_units[60]      = 1e-12;

    // (59):  CH2 + CH3 <=> H + C2H4
    kiv[61] = {7,9,1,16};
    nuv[61] = {-1,-1,1,1};
    fwd_A[61]     = 40000000000000;
    fwd_beta[61]  = 0;
    fwd_Ea[61]    = 0;
    prefactor_units[61]  = 1.0000000000000002e-06;
    activation_units[61] = 0.50321666580471969;
    phase_units[61]      = 1e-12;

    // (60):  CH2 + CH4 <=> 2 CH3
    kiv[62] = {7,10,9};
    nuv[62] = {-1,-1,2};
    fwd_A[62]     = 2460000;
    fwd_beta[62]  = 2;
    fwd_Ea[62]    = 8270;
    prefactor_units[62]  = 1.0000000000000002e-06;
    activation_units[62] = 0.50321666580471969;
    phase_units[62]      = 1e-12;

    // (61):  CH2(S) + N2 <=> CH2 + N2
    kiv[63] = {8,19,7,19};
    nuv[63] = {-1,-1,1,1};
    fwd_A[63]     = 15000000000000;
    fwd_beta[63]  = 0;
    fwd_Ea[63]    = 600;
    prefactor_units[63]  = 1.0000000000000002e-06;
    activation_units[63] = 0.50321666580471969;
    phase_units[63]      = 1e-12;

    // (62):  CH2(S) + AR <=> CH2 + AR
    kiv[64] = {8,20,7,20};
    nuv[64] = {-1,-1,1,1};
    fwd_A[64]     = 9000000000000;
    fwd_beta[64]  = 0;
    fwd_Ea[64]    = 600;
    prefactor_units[64]  = 1.0000000000000002e-06;
    activation_units[64] = 0.50321666580471969;
    phase_units[64]      = 1e-12;

    // (63):  CH2(S) + O2 <=> H + OH + CO
    kiv[65] = {8,3,1,4,11};
    nuv[65] = {-1,-1,1,1,1};
    fwd_A[65]     = 28000000000000;
    fwd_beta[65]  = 0;
    fwd_Ea[65]    = 0;
    prefactor_units[65]  = 1.0000000000000002e-06;
    activation_units[65] = 0.50321666580471969;
    phase_units[65]      = 1e-12;

    // (64):  CH2(S) + O2 <=> CO + H2O
    kiv[66] = {8,3,11,5};
    nuv[66] = {-1,-1,1,1};
    fwd_A[66]     = 12000000000000;
    fwd_beta[66]  = 0;
    fwd_Ea[66]    = 0;
    prefactor_units[66]  = 1.0000000000000002e-06;
    activation_units[66] = 0.50321666580471969;
    phase_units[66]      = 1e-12;

    // (65):  CH2(S) + H2 <=> CH3 + H
    kiv[67] = {8,0,9,1};
    nuv[67] = {-1,-1,1,1};
    fwd_A[67]     = 70000000000000;
    fwd_beta[67]  = 0;
    fwd_Ea[67]    = 0;
    prefactor_units[67]  = 1.0000000000000002e-06;
    activation_units[67] = 0.50321666580471969;
    phase_units[67]      = 1e-12;

    // (66):  CH2(S) + H2O <=> CH2 + H2O
    kiv[68] = {8,5,7,5};
    nuv[68] = {-1,-1,1,1};
    fwd_A[68]     = 30000000000000;
    fwd_beta[68]  = 0;
    fwd_Ea[68]    = 0;
    prefactor_units[68]  = 1.0000000000000002e-06;
    activation_units[68] = 0.50321666580471969;
    phase_units[68]      = 1e-12;

    // (67):  CH2(S) + CH3 <=> H + C2H4
    kiv[69] = {8,9,1,16};
    nuv[69] = {-1,-1,1,1};
    fwd_A[69]     = 12000000000000;
    fwd_beta[69]  = 0;
    fwd_Ea[69]    = -570;
    prefactor_units[69]  = 1.0000000000000002e-06;
    activation_units[69] = 0.50321666580471969;
    phase_units[69]      = 1e-12;

    // (68):  CH2(S) + CH4 <=> 2 CH3
    kiv[70] = {8,10,9};
    nuv[70] = {-1,-1,2};
    fwd_A[70]     = 16000000000000;
    fwd_beta[70]  = 0;
    fwd_Ea[70]    = -570;
    prefactor_units[70]  = 1.0000000000000002e-06;
    activation_units[70] = 0.50321666580471969;
    phase_units[70]      = 1e-12;

    // (69):  CH2(S) + CO <=> CH2 + CO
    kiv[71] = {8,11,7,11};
    nuv[71] = {-1,-1,1,1};
    fwd_A[71]     = 9000000000000;
    fwd_beta[71]  = 0;
    fwd_Ea[71]    = 0;
    prefactor_units[71]  = 1.0000000000000002e-06;
    activation_units[71] = 0.50321666580471969;
    phase_units[71]      = 1e-12;

    // (70):  CH2(S) + CO2 <=> CH2 + CO2
    kiv[72] = {8,12,7,12};
    nuv[72] = {-1,-1,1,1};
    fwd_A[72]     = 7000000000000;
    fwd_beta[72]  = 0;
    fwd_Ea[72]    = 0;
    prefactor_units[72]  = 1.0000000000000002e-06;
    activation_units[72] = 0.50321666580471969;
    phase_units[72]      = 1e-12;

    // (71):  CH2(S) + CO2 <=> CO + CH2O
    kiv[73] = {8,12,11,14};
    nuv[73] = {-1,-1,1,1};
    fwd_A[73]     = 14000000000000;
    fwd_beta[73]  = 0;
    fwd_Ea[73]    = 0;
    prefactor_units[73]  = 1.0000000000000002e-06;
    activation_units[73] = 0.50321666580471969;
    phase_units[73]      = 1e-12;

    // (72):  CH3 + O2 <=> O + CH3O
    kiv[74] = {9,3,2,15};
    nuv[74] = {-1,-1,1,1};
    fwd_A[74]     = 26750000000000;
    fwd_beta[74]  = 0;
    fwd_Ea[74]    = 28800;
    prefactor_units[74]  = 1.0000000000000002e-06;
    activation_units[74] = 0.50321666580471969;
    phase_units[74]      = 1e-12;

    // (73):  CH3 + O2 <=> OH + CH2O
    kiv[75] = {9,3,4,14};
    nuv[75] = {-1,-1,1,1};
    fwd_A[75]     = 36000000000;
    fwd_beta[75]  = 0;
    fwd_Ea[75]    = 8940;
    prefactor_units[75]  = 1.0000000000000002e-06;
    activation_units[75] = 0.50321666580471969;
    phase_units[75]      = 1e-12;

    // (74):  2 CH3 (+M) <=> C2H6 (+M)
    kiv[7] = {9,18};
    nuv[7] = {-2,1};
    fwd_A[7]     = 21200000000000000;
    fwd_beta[7]  = -0.96999999999999997;
    fwd_Ea[7]    = 620;
    prefactor_units[7]  = 1.0000000000000002e-06;
    activation_units[7] = 0.50321666580471969;
    phase_units[7]      = 1e-12;

    // (75):  2 CH3 <=> H + C2H5
    kiv[76] = {9,1,17};
    nuv[76] = {-2,1,1};
    fwd_A[76]     = 4990000000000;
    fwd_beta[76]  = 0.10000000000000001;
    fwd_Ea[76]    = 10600;
    prefactor_units[76]  = 1.0000000000000002e-06;
    activation_units[76] = 0.50321666580471969;
    phase_units[76]      = 1e-12;

    // (76):  CH3 + HCO <=> CH4 + CO
    kiv[77] = {9,13,10,11};
    nuv[77] = {-1,-1,1,1};
    fwd_A[77]     = 26480000000000;
    fwd_beta[77]  = 0;
    fwd_Ea[77]    = 0;
    prefactor_units[77]  = 1.0000000000000002e-06;
    activation_units[77] = 0.50321666580471969;
    phase_units[77]      = 1e-12;

    // (77):  CH3 + CH2O <=> HCO + CH4
    kiv[78] = {9,14,13,10};
    nuv[78] = {-1,-1,1,1};
    fwd_A[78]     = 3320;
    fwd_beta[78]  = 2.8100000000000001;
    fwd_Ea[78]    = 5860;
    prefactor_units[78]  = 1.0000000000000002e-06;
    activation_units[78] = 0.50321666580471969;
    phase_units[78]      = 1e-12;

    // (78):  CH3 + C2H6 <=> C2H5 + CH4
    kiv[79] = {9,18,17,10};
    nuv[79] = {-1,-1,1,1};
    fwd_A[79]     = 6140000;
    fwd_beta[79]  = 1.74;
    fwd_Ea[79]    = 10450;
    prefactor_units[79]  = 1.0000000000000002e-06;
    activation_units[79] = 0.50321666580471969;
    phase_units[79]      = 1e-12;

    // (79):  HCO + H2O <=> H + CO + H2O
    kiv[80] = {13,5,1,11,5};
    nuv[80] = {-1,-1,1,1,1};
    fwd_A[80]     = 2.244e+18;
    fwd_beta[80]  = -1;
    fwd_Ea[80]    = 17000;
    prefactor_units[80]  = 1.0000000000000002e-06;
    activation_units[80] = 0.50321666580471969;
    phase_units[80]      = 1e-12;

    // (80):  HCO + M <=> H + CO + M
    kiv[13] = {13,1,11};
    nuv[13] = {-1,1,1};
    fwd_A[13]     = 1.87e+17;
    fwd_beta[13]  = -1;
    fwd_Ea[13]    = 17000;
    prefactor_units[13]  = 1.0000000000000002e-06;
    activation_units[13] = 0.50321666580471969;
    phase_units[13]      = 1e-6;

    // (81):  HCO + O2 <=> HO2 + CO
    kiv[81] = {13,3,6,11};
    nuv[81] = {-1,-1,1,1};
    fwd_A[81]     = 7600000000000;
    fwd_beta[81]  = 0;
    fwd_Ea[81]    = 400;
    prefactor_units[81]  = 1.0000000000000002e-06;
    activation_units[81] = 0.50321666580471969;
    phase_units[81]      = 1e-12;

    // (82):  CH3O + O2 <=> HO2 + CH2O
    kiv[82] = {15,3,6,14};
    nuv[82] = {-1,-1,1,1};
    fwd_A[82]     = 4.2799999999999999e-13;
    fwd_beta[82]  = 7.5999999999999996;
    fwd_Ea[82]    = -3530;
    prefactor_units[82]  = 1.0000000000000002e-06;
    activation_units[82] = 0.50321666580471969;
    phase_units[82]      = 1e-12;

    // (83):  C2H5 + O2 <=> HO2 + C2H4
    kiv[83] = {17,3,6,16};
    nuv[83] = {-1,-1,1,1};
    fwd_A[83]     = 840000000000;
    fwd_beta[83]  = 0;
    fwd_Ea[83]    = 3875;
    prefactor_units[83]  = 1.0000000000000002e-06;
    activation_units[83] = 0.50321666580471969;
    phase_units[83]      = 1e-12;

    int max_nu = 5;
    for (int i=0; i<84; ++i) {
      size_t offset = max_nu * i;
      for (int j=0; j<nuv[i].size(); ++j) {
        nu2D[offset+j] = nuv[i][j];
        ki2D[offset+j] = kiv[i][j];
      }
      for (int j=nuv[i].size(); j<max_nu; ++j) {
        nu2D[offset+j] = 0;
        ki2D[offset+j] = 0;
      }
    }


#ifdef AMREX_USE_CUDA
    AMREX_CUDA_SAFE_CALL(cudaMemPrefetchAsync(fwd_A, 84, Gpu::Device::deviceId(), Gpu::gpuStream()));
    AMREX_CUDA_SAFE_CALL(cudaMemPrefetchAsync(fwd_beta, 84, Gpu::Device::deviceId(), Gpu::gpuStream()));
    AMREX_CUDA_SAFE_CALL(cudaMemPrefetchAsync(fwd_Ea, 84, Gpu::Device::deviceId(), Gpu::gpuStream()));
    AMREX_CUDA_SAFE_CALL(cudaMemPrefetchAsync(prefactor_units, 84, Gpu::Device::deviceId(), Gpu::gpuStream()));
    AMREX_CUDA_SAFE_CALL(cudaMemPrefetchAsync(activation_units, 84, Gpu::Device::deviceId(), Gpu::gpuStream()));
    AMREX_CUDA_SAFE_CALL(cudaMemPrefetchAsync(phase_units, 84, Gpu::Device::deviceId(), Gpu::gpuStream()));
    AMREX_CUDA_SAFE_CALL(cudaMemPrefetchAsync(nu2D, 420, Gpu::Device::deviceId(), Gpu::gpuStream()));
    AMREX_CUDA_SAFE_CALL(cudaMemPrefetchAsync(ki2D, 420, Gpu::Device::deviceId(), Gpu::gpuStream()));
    AMREX_CUDA_SAFE_CALL(cudaMemPrefetchAsync(imw, 21, Gpu::Device::deviceId(), Gpu::gpuStream()));
    for (int n=0; n<NUM_REACTIONS; ++n) {
      AMREX_CUDA_SAFE_CALL(cudaMemPrefetchAsync(&(kiv[n][0]), 5, Gpu::Device::deviceId(), Gpu::gpuStream()));
      AMREX_CUDA_SAFE_CALL(cudaMemPrefetchAsync(&(nuv[n][0]), 5, Gpu::Device::deviceId(), Gpu::gpuStream()));
    }

#endif

}


AMREX_GPU_DEVICE
void 
W_spec_d(int i, int j, int k,
         amrex::Array4<amrex::Real> const& rho,
         amrex::Array4<amrex::Real> const& temp,
         amrex::Array4<amrex::Real> const& Y,
         amrex::Array4<amrex::Real> const& cdots)
{
  amrex::GpuArray<amrex::Real,NUM_SPECIES> sc;
  amrex::GpuArray<amrex::Real,NUM_SPECIES> wdot;
  amrex::GpuArray<amrex::Real,NUM_REACTIONS> qf;
  amrex::Real invT, logT;

  invT = 1/temp(i,j,k);
  logT = log(temp(i,j,k));

  for (int n = 0; n < NUM_SPECIES; n++) {
    sc[n] = Y(i,j,k,n) * rho(i,j,k) * imw[n];
    wdot[n] = 0;
  }

  for (int n=0; n<NUM_REACTIONS; ++n) {
    size_t offset = n*5;
    qf[n] = 1;
    for (int m=0; m<5; ++m) {
      int nu = nu2D[offset+m];
      if (nu<0) {
        qf[n] *= pow(sc[ki2D[offset+m]],-nu);
      }
    }
    qf[n] *=  prefactor_units[n] * fwd_A[n] * exp(fwd_beta[n] * logT - activation_units[n] * fwd_Ea[n] * invT);
    for (int m=0; m<5; ++m) {
      int nu = nu2D[offset+m];
      if (nu<0) {
        wdot[ki2D[offset+m]] += nu*qf[n];
      }
    }
  }

  for (int n = 0; n < NUM_SPECIES; n++) {
    cdots(i,j,k,n) = -wdot[n];
  }
}

AMREX_GPU_DEVICE 
void 
W_spec1_d(int i, int j, int k,
		amrex::Array4<amrex::Real> const& rho,
		amrex::Array4<amrex::Real> const& temp,
		amrex::Array4<amrex::Real> const& Y,
		amrex::Array4<amrex::Real> const& cdots)
{

	amrex::GpuArray<amrex::Real,NUM_SPECIES> sc;
	amrex::GpuArray<amrex::Real,NUM_REACTIONS> qf;
	amrex::GpuArray<amrex::Real,NUM_SPECIES> wdot;
	amrex::Real imw[NUM_SPECIES];
	amrex::Real invT, logT, k_f, qdot, Corr;

	get_imw(imw);

	invT = 1/temp(i,j,k);
	logT = log(temp(i,j,k));

	for (int n = 0; n < NUM_SPECIES; n++) {
	        sc[n] = Y(i,j,k,n) * rho(i,j,k) * imw[n];
	}
#if 0
        /*reaction 1: H + CH2 (+M) <=> CH3 (+M) */
        qf[0] = sc[1]*sc[7];
        /*reaction 2: H + CH3 (+M) <=> CH4 (+M) */
        qf[1] = sc[1]*sc[9];
        /*reaction 3: H + HCO (+M) <=> CH2O (+M) */
        qf[2] = sc[1]*sc[13];
        /*reaction 4: H + CH2O (+M) <=> CH3O (+M) */
        qf[3] = sc[1]*sc[14];
        /*reaction 5: H + C2H4 (+M) <=> C2H5 (+M) */
        qf[4] = sc[1]*sc[16];
        /*reaction 6: H + C2H5 (+M) <=> C2H6 (+M) */
        qf[5] = sc[1]*sc[17];
        /*reaction 7: H2 + CO (+M) <=> CH2O (+M) */
        qf[6] = sc[0]*sc[11];
        /*reaction 8: 2 CH3 (+M) <=> C2H6 (+M) */
        qf[7] = sc[9]*sc[9];
        /*reaction 9: O + H + M <=> OH + M */
        qf[8] = sc[1]*sc[2];
        /*reaction 10: O + CO + M <=> CO2 + M */
        qf[9] = sc[2]*sc[11];
        /*reaction 11: H + O2 + M <=> HO2 + M */
        qf[10] = sc[1]*sc[3];
        /*reaction 12: 2 H + M <=> H2 + M */
        qf[11] = sc[1]*sc[1];
        /*reaction 13: H + OH + M <=> H2O + M */
        qf[12] = sc[1]*sc[4];
        /*reaction 14: HCO + M <=> H + CO + M */
        qf[13] = sc[13];
        /*reaction 15: O + H2 <=> H + OH */
        qf[14] = sc[0]*sc[2];
        /*reaction 16: O + HO2 <=> OH + O2 */
        qf[15] = sc[2]*sc[6];
        /*reaction 17: O + CH2 <=> H + HCO */
        qf[16] = sc[2]*sc[7];
        /*reaction 18: O + CH2(S) <=> H + HCO */
        qf[17] = sc[2]*sc[8];
        /*reaction 19: O + CH3 <=> H + CH2O */
        qf[18] = sc[2]*sc[9];
        /*reaction 20: O + CH4 <=> OH + CH3 */
        qf[19] = sc[2]*sc[10];
        /*reaction 21: O + HCO <=> OH + CO */
        qf[20] = sc[2]*sc[13];
        /*reaction 22: O + HCO <=> H + CO2 */
        qf[21] = sc[2]*sc[13];
        /*reaction 23: O + CH2O <=> OH + HCO */
        qf[22] = sc[2]*sc[14];
        /*reaction 24: O + C2H4 <=> CH3 + HCO */
        qf[23] = sc[2]*sc[16];
        /*reaction 25: O + C2H5 <=> CH3 + CH2O */
        qf[24] = sc[2]*sc[17];
        /*reaction 26: O + C2H6 <=> OH + C2H5 */
        qf[25] = sc[2]*sc[18];
        /*reaction 27: O2 + CO <=> O + CO2 */
        qf[26] = sc[3]*sc[11];
        /*reaction 28: O2 + CH2O <=> HO2 + HCO */
        qf[27] = sc[3]*sc[14];
        /*reaction 29: H + 2 O2 <=> HO2 + O2 */
        qf[28] = sc[1]*sc[3]*sc[3];
        /*reaction 30: H + O2 + H2O <=> HO2 + H2O */
        qf[29] = sc[1]*sc[3]*sc[5];
        /*reaction 31: H + O2 + N2 <=> HO2 + N2 */
        qf[30] = sc[1]*sc[3]*sc[19];
        /*reaction 32: H + O2 + AR <=> HO2 + AR */
        qf[31] = sc[1]*sc[3]*sc[20];
        /*reaction 33: H + O2 <=> O + OH */
        qf[32] = sc[1]*sc[3];
        /*reaction 34: 2 H + H2 <=> 2 H2 */
        qf[33] = sc[0]*sc[1]*sc[1];
        /*reaction 35: 2 H + H2O <=> H2 + H2O */
        qf[34] = sc[1]*sc[1]*sc[5];
        /*reaction 36: 2 H + CO2 <=> H2 + CO2 */
        qf[35] = sc[1]*sc[1]*sc[12];
        /*reaction 37: H + HO2 <=> O2 + H2 */
        qf[36] = sc[1]*sc[6];
        /*reaction 38: H + HO2 <=> 2 OH */
        qf[37] = sc[1]*sc[6];
        /*reaction 39: H + CH4 <=> CH3 + H2 */
        qf[38] = sc[1]*sc[10];
        /*reaction 40: H + HCO <=> H2 + CO */
        qf[39] = sc[1]*sc[13];
        /*reaction 41: H + CH2O <=> HCO + H2 */
        qf[40] = sc[1]*sc[14];
        /*reaction 42: H + CH3O <=> OH + CH3 */
        qf[41] = sc[1]*sc[15];
        /*reaction 43: H + C2H6 <=> C2H5 + H2 */
        qf[42] = sc[1]*sc[18];
        /*reaction 44: OH + H2 <=> H + H2O */
        qf[43] = sc[0]*sc[4];
        /*reaction 45: 2 OH <=> O + H2O */
        qf[44] = sc[4]*sc[4];
        /*reaction 46: OH + HO2 <=> O2 + H2O */
        qf[45] = sc[4]*sc[6];
        /*reaction 47: OH + CH2 <=> H + CH2O */
        qf[46] = sc[4]*sc[7];
        /*reaction 48: OH + CH2(S) <=> H + CH2O */
        qf[47] = sc[4]*sc[8];
        /*reaction 49: OH + CH3 <=> CH2 + H2O */
        qf[48] = sc[4]*sc[9];
        /*reaction 50: OH + CH3 <=> CH2(S) + H2O */
        qf[49] = sc[4]*sc[9];
        /*reaction 51: OH + CH4 <=> CH3 + H2O */
        qf[50] = sc[4]*sc[10];
        /*reaction 52: OH + CO <=> H + CO2 */
        qf[51] = sc[4]*sc[11];
        /*reaction 53: OH + HCO <=> H2O + CO */
        qf[52] = sc[4]*sc[13];
        /*reaction 54: OH + CH2O <=> HCO + H2O */
        qf[53] = sc[4]*sc[14];
        /*reaction 55: OH + C2H6 <=> C2H5 + H2O */
        qf[54] = sc[4]*sc[18];
        /*reaction 56: HO2 + CH2 <=> OH + CH2O */
        qf[55] = sc[6]*sc[7];
        /*reaction 57: HO2 + CH3 <=> O2 + CH4 */
        qf[56] = sc[6]*sc[9];
        /*reaction 58: HO2 + CH3 <=> OH + CH3O */
        qf[57] = sc[6]*sc[9];
        /*reaction 59: HO2 + CO <=> OH + CO2 */
        qf[58] = sc[6]*sc[11];
        /*reaction 60: CH2 + O2 <=> OH + HCO */
        qf[59] = sc[3]*sc[7];
        /*reaction 61: CH2 + H2 <=> H + CH3 */
        qf[60] = sc[0]*sc[7];
        /*reaction 62: CH2 + CH3 <=> H + C2H4 */
        qf[61] = sc[7]*sc[9];
        /*reaction 63: CH2 + CH4 <=> 2 CH3 */
        qf[62] = sc[7]*sc[10];
        /*reaction 64: CH2(S) + N2 <=> CH2 + N2 */
        qf[63] = sc[8]*sc[19];
        /*reaction 65: CH2(S) + AR <=> CH2 + AR */
        qf[64] = sc[8]*sc[20];
        /*reaction 66: CH2(S) + O2 <=> H + OH + CO */
        qf[65] = sc[3]*sc[8];
        /*reaction 67: CH2(S) + O2 <=> CO + H2O */
        qf[66] = sc[3]*sc[8];
        /*reaction 68: CH2(S) + H2 <=> CH3 + H */
        qf[67] = sc[0]*sc[8];
        /*reaction 69: CH2(S) + H2O <=> CH2 + H2O */
        qf[68] = sc[5]*sc[8];
        /*reaction 70: CH2(S) + CH3 <=> H + C2H4 */
        qf[69] = sc[8]*sc[9];
        /*reaction 71: CH2(S) + CH4 <=> 2 CH3 */
        qf[70] = sc[8]*sc[10];
        /*reaction 72: CH2(S) + CO <=> CH2 + CO */
        qf[71] = sc[8]*sc[11];
        /*reaction 73: CH2(S) + CO2 <=> CH2 + CO2 */
        qf[72] = sc[8]*sc[12];
        /*reaction 74: CH2(S) + CO2 <=> CO + CH2O */
        qf[73] = sc[8]*sc[12];
        /*reaction 75: CH3 + O2 <=> O + CH3O */
        qf[74] = sc[3]*sc[9];
        /*reaction 76: CH3 + O2 <=> OH + CH2O */
        qf[75] = sc[3]*sc[9];
        /*reaction 77: 2 CH3 <=> H + C2H5 */
        qf[76] = sc[9]*sc[9];
        /*reaction 78: CH3 + HCO <=> CH4 + CO */
        qf[77] = sc[9]*sc[13];
        /*reaction 79: CH3 + CH2O <=> HCO + CH4 */
        qf[78] = sc[9]*sc[14];
        /*reaction 80: CH3 + C2H6 <=> C2H5 + CH4 */
        qf[79] = sc[9]*sc[18];
        /*reaction 81: HCO + H2O <=> H + CO + H2O */
        qf[80] = sc[5]*sc[13];
        /*reaction 82: HCO + O2 <=> HO2 + CO */
        qf[81] = sc[3]*sc[13];
        /*reaction 83: CH3O + O2 <=> HO2 + CH2O */
        qf[82] = sc[3]*sc[15];
        /*reaction 84: C2H5 + O2 <=> HO2 + C2H4 */
        qf[83] = sc[3]*sc[17];
#else
  for (int n=0; n<NUM_REACTIONS; ++n) {
    int nu0 = nuv[n][0];
    int nu1 = nuv[n][1];
    int ki0 = kiv[n][0];
    int ki1 = kiv[n][1];
    qf[n] = (nu0<0 ? sc[ki0] : 1) * (nu0<1 ? sc[ki0] : 1)
          * (nu1<0 ? sc[ki1] : 1) * (nu1<1 ? sc[ki0] : 1);
 }
#endif

        // (0):  O + H + M <=> OH + M
        qf[8] *= 1.0000000000000002e-12 * 5e+17 
                   * exp(-1 * logT - 0.50321666580471969 * 0 * invT);
        // (1):  O + H2 <=> H + OH
        qf[14] *= 1.0000000000000002e-06 * 50000 
               * exp(2.6699999999999999 * logT - 0.50321666580471969 * 6290 * invT);
        // (2):  O + HO2 <=> OH + O2
        qf[15] *= 1.0000000000000002e-06 * 20000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        // (3):  O + CH2 <=> H + HCO
        qf[16] *= 1.0000000000000002e-06 * 80000000000000 
               * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        // (4):  O + CH2(S) <=> H + HCO
        qf[17] *= 1.0000000000000002e-06 * 15000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        // (5):  O + CH3 <=> H + CH2O
        qf[18] *= 1.0000000000000002e-06 * 84300000000000 
               * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        // (6):  O + CH4 <=> OH + CH3
        qf[19] *= 1.0000000000000002e-06 * 1020000000 
                   * exp(1.5 * logT - 0.50321666580471969 * 8600 * invT);
        // (7):  O + CO + M <=> CO2 + M
        qf[9] *= 1.0000000000000002e-12 * 602000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 3000 * invT);
        // (8):  O + HCO <=> OH + CO
        qf[20] *= 1.0000000000000002e-06 * 30000000000000 
               * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        // (9):  O + HCO <=> H + CO2
        qf[21] *= 1.0000000000000002e-06 * 30000000000000
                   * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        // (10):  O + CH2O <=> OH + HCO
        qf[22] *= 1.0000000000000002e-06 * 39000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 3540 * invT);
        // (11):  O + C2H4 <=> CH3 + HCO
        qf[23] *= 1.0000000000000002e-06 * 19200000 
                   * exp(1.8300000000000001 * logT - 0.50321666580471969 * 220 * invT);
        // (12):  O + C2H5 <=> CH3 + CH2O
        qf[24] *= 1.0000000000000002e-06 * 132000000000000 
               * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        // (13):  O + C2H6 <=> OH + C2H5
        qf[25] *= 1.0000000000000002e-06 * 89800000 
                   * exp(1.9199999999999999 * logT - 0.50321666580471969 * 5690 * invT);
        // (14):  O2 + CO <=> O + CO2
        qf[26] *= 1.0000000000000002e-06 * 2500000000000 
                   * exp(0 * logT - 0.50321666580471969 * 47800 * invT);
        // (15):  O2 + CH2O <=> HO2 + HCO
        qf[27] *= 1.0000000000000002e-06 * 100000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 40000 * invT);
        // (16):  H + O2 + M <=> HO2 + M
        qf[10] *= 1.0000000000000002e-12 * 2.8e+18 
                   * exp(-0.85999999999999999 * logT - 0.50321666580471969 * 0 * invT);
        // (17):  H + 2 O2 <=> HO2 + O2
        qf[28] *= 1.0000000000000002e-12 * 3e+20 
                   * exp(-1.72 * logT - 0.50321666580471969 * 0 * invT);
        // (18):  H + O2 + H2O <=> HO2 + H2O
        qf[29] *= 1.0000000000000002e-12 * 9.38e+18 
                   * exp(-0.76000000000000001 * logT - 0.50321666580471969 * 0 * invT);
        // (19):  H + O2 + N2 <=> HO2 + N2
	qf[30] *= 1.0000000000000002e-12 * 3.75e+20
                   * exp(-1.72 * logT - 0.50321666580471969 * 0 * invT);
        // (20):  H + O2 + AR <=> HO2 + AR
        qf[31] *= 1.0000000000000002e-12 * 7e+17 
                   * exp(-0.80000000000000004 * logT - 0.50321666580471969 * 0 * invT);
        // (21):  H + O2 <=> O + OH
        qf[32] *= 1.0000000000000002e-06 * 83000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 14413 * invT);
        // (22):  2 H + M <=> H2 + M
        qf[11] *= 1.0000000000000002e-12 * 1e+18 
                   * exp(-1 * logT - 0.50321666580471969 * 0 * invT);
        // (23):  2 H + H2 <=> 2 H2
        qf[33] *= 1.0000000000000002e-12 * 90000000000000000 
                   * exp(-0.59999999999999998 * logT - 0.50321666580471969 * 0 * invT);
        // (24):  2 H + H2O <=> H2 + H2O
        qf[34] *= 1.0000000000000002e-12 * 6e+19 
                   * exp(-1.25 * logT - 0.50321666580471969 * 0 * invT);
        // (25):  2 H + CO2 <=> H2 + CO2
        qf[35] *= 1.0000000000000002e-12 * 5.5e+20 
                   * exp(-2 * logT - 0.50321666580471969 * 0 * invT);
        // (26):  H + OH + M <=> H2O + M
        qf[12] *= 1.0000000000000002e-12 * 2.2e+22 
                   * exp(-2 * logT - 0.50321666580471969 * 0 * invT);
        // (27):  H + HO2 <=> O2 + H2
        qf[36] *= 1.0000000000000002e-06 * 28000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 1068 * invT);
        // (28):  H + HO2 <=> 2 OH
        qf[37] *= 1.0000000000000002e-06 * 134000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 635 * invT);
        // (29):  H + CH2 (+M) <=> CH3 (+M)
        qf[0] *= 1.0000000000000002e-06 * 25000000000000000 
                   * exp(-0.80000000000000004 * logT - 0.50321666580471969 * 0 * invT);
        // (30):  H + CH3 (+M) <=> CH4 (+M)
        qf[1] *= 1.0000000000000002e-06 * 12700000000000000 
                   * exp(-0.63 * logT - 0.50321666580471969 * 383 * invT);
        // (31):  H + CH4 <=> CH3 + H2
        qf[38] *= 1.0000000000000002e-06 * 660000000 
                   * exp(1.6200000000000001 * logT - 0.50321666580471969 * 10840 * invT);
        // (32):  H + HCO (+M) <=> CH2O (+M)
        qf[2] *= 1.0000000000000002e-06 * 1090000000000 
                   * exp(0.47999999999999998 * logT - 0.50321666580471969 * -260 * invT);
        // (33):  H + HCO <=> H2 + CO
        qf[39] *= 1.0000000000000002e-06 * 73400000000000 
                   * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        // (34):  H + CH2O (+M) <=> CH3O (+M)
        qf[3] *= 1.0000000000000002e-06 * 540000000000 
                   * exp(0.45400000000000001 * logT - 0.50321666580471969 * 2600 * invT);
        // (35):  H + CH2O <=> HCO + H2
        qf[40] *= 1.0000000000000002e-06 * 23000000000 
                   * exp(1.05 * logT - 0.50321666580471969 * 3275 * invT);
        // (36):  H + CH3O <=> OH + CH3
        qf[41] *= 1.0000000000000002e-06 * 32000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        // (37):  H + C2H4 (+M) <=> C2H5 (+M)
        qf[4] *= 1.0000000000000002e-06 * 1080000000000 
                   * exp(0.45400000000000001 * logT - 0.50321666580471969 * 1820 * invT);
        // (38):  H + C2H5 (+M) <=> C2H6 (+M)
        qf[5] *= 1.0000000000000002e-06 * 5.21e+17 
                   * exp(-0.98999999999999999 * logT - 0.50321666580471969 * 1580 * invT);
        // (39):  H + C2H6 <=> C2H5 + H2
        qf[42] *= 1.0000000000000002e-06 * 115000000 
                   * exp(1.8999999999999999 * logT - 0.50321666580471969 * 7530 * invT);
        // (40):  H2 + CO (+M) <=> CH2O (+M)
        qf[6] *= 1.0000000000000002e-06 * 43000000 
                   * exp(1.5 * logT - 0.50321666580471969 * 79600 * invT);
        // (41):  OH + H2 <=> H + H2O
        qf[43] *= 1.0000000000000002e-06 * 216000000 
                   * exp(1.51 * logT - 0.50321666580471969 * 3430 * invT);
        // (42):  2 OH <=> O + H2O
        qf[44] *= 1.0000000000000002e-06 * 35700 
                   * exp(2.3999999999999999 * logT - 0.50321666580471969 * -2110 * invT);
        // (43):  OH + HO2 <=> O2 + H2O
        qf[45] *= 1.0000000000000002e-06 * 29000000000000 
                   * exp(0 * logT - 0.50321666580471969 * -500 * invT);
        // (44):  OH + CH2 <=> H + CH2O
        k_f = 1.0000000000000002e-06 * 20000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        Corr  = 1.0;
        qf[46] *= Corr * k_f;
        // (45):  OH + CH2(S) <=> H + CH2O
        k_f = 1.0000000000000002e-06 * 30000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        Corr  = 1.0;
        qf[47] *= Corr * k_f;
        // (46):  OH + CH3 <=> CH2 + H2O
        k_f = 1.0000000000000002e-06 * 56000000 
                   * exp(1.6000000000000001 * logT - 0.50321666580471969 * 5420 * invT);
        Corr  = 1.0;
        qf[48] *= Corr * k_f;
        // (47):  OH + CH3 <=> CH2(S) + H2O
        k_f = 1.0000000000000002e-06 * 25010000000000 
                   * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        Corr  = 1.0;
        qf[49] *= Corr * k_f;
        // (48):  OH + CH4 <=> CH3 + H2O
        k_f = 1.0000000000000002e-06 * 100000000 
                   * exp(1.6000000000000001 * logT - 0.50321666580471969 * 3120 * invT);
        Corr  = 1.0;
        qf[50] *= Corr * k_f;
        // (49):  OH + CO <=> H + CO2
        k_f = 1.0000000000000002e-06 * 47600000 
                   * exp(1.228 * logT - 0.50321666580471969 * 70 * invT);
        Corr  = 1.0;
        qf[51] *= Corr * k_f;
        // (50):  OH + HCO <=> H2O + CO
        k_f = 1.0000000000000002e-06 * 50000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        Corr  = 1.0;
        qf[52] *= Corr * k_f;
        // (51):  OH + CH2O <=> HCO + H2O
        k_f = 1.0000000000000002e-06 * 3430000000 
                   * exp(1.1799999999999999 * logT - 0.50321666580471969 * -447 * invT);
        Corr  = 1.0;
        qf[53] *= Corr * k_f;
        // (52):  OH + C2H6 <=> C2H5 + H2O
        k_f = 1.0000000000000002e-06 * 3540000 
                   * exp(2.1200000000000001 * logT - 0.50321666580471969 * 870 * invT);
        Corr  = 1.0;
        qf[54] *= Corr * k_f;
        // (53):  HO2 + CH2 <=> OH + CH2O
        k_f = 1.0000000000000002e-06 * 20000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        Corr  = 1.0;
        qf[55] *= Corr * k_f;
        // (54):  HO2 + CH3 <=> O2 + CH4
        k_f = 1.0000000000000002e-06 * 1000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        Corr  = 1.0;
        qf[56] *= Corr * k_f;
        // (55):  HO2 + CH3 <=> OH + CH3O
        k_f = 1.0000000000000002e-06 * 20000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        Corr  = 1.0;
        qf[57] *= Corr * k_f;
        // (56):  HO2 + CO <=> OH + CO2
        k_f = 1.0000000000000002e-06 * 150000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 23600 * invT);
        Corr  = 1.0;
        qf[58] *= Corr * k_f;
        // (57):  CH2 + O2 <=> OH + HCO
        k_f = 1.0000000000000002e-06 * 13200000000000 
                   * exp(0 * logT - 0.50321666580471969 * 1500 * invT);
        Corr  = 1.0;
        qf[59] *= Corr * k_f;
        // (58):  CH2 + H2 <=> H + CH3
        k_f = 1.0000000000000002e-06 * 500000 
                   * exp(2 * logT - 0.50321666580471969 * 7230 * invT);
        Corr  = 1.0;
        qf[60] *= Corr * k_f;
        // (59):  CH2 + CH3 <=> H + C2H4
        k_f = 1.0000000000000002e-06 * 40000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        Corr  = 1.0;
        qf[61] *= Corr * k_f;
        // (60):  CH2 + CH4 <=> 2 CH3
        k_f = 1.0000000000000002e-06 * 2460000 
                   * exp(2 * logT - 0.50321666580471969 * 8270 * invT);
        Corr  = 1.0;
        qf[62] *= Corr * k_f;
        // (61):  CH2(S) + N2 <=> CH2 + N2
        k_f = 1.0000000000000002e-06 * 15000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 600 * invT);
        Corr  = 1.0;
        qf[63] *= Corr * k_f;
        // (62):  CH2(S) + AR <=> CH2 + AR
        k_f = 1.0000000000000002e-06 * 9000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 600 * invT);
        Corr  = 1.0;
        qf[64] *= Corr * k_f;
        // (63):  CH2(S) + O2 <=> H + OH + CO
        k_f = 1.0000000000000002e-06 * 28000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        Corr  = 1.0;
        qf[65] *= Corr * k_f;
        // (64):  CH2(S) + O2 <=> CO + H2O
        k_f = 1.0000000000000002e-06 * 12000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        Corr  = 1.0;
        qf[66] *= Corr * k_f;
        // (65):  CH2(S) + H2 <=> CH3 + H
        k_f = 1.0000000000000002e-06 * 70000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        Corr  = 1.0;
        qf[67] *= Corr * k_f;
        // (66):  CH2(S) + H2O <=> CH2 + H2O
        k_f = 1.0000000000000002e-06 * 30000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        Corr  = 1.0;
        qf[68] *= Corr * k_f;
        // (67):  CH2(S) + CH3 <=> H + C2H4
        k_f = 1.0000000000000002e-06 * 12000000000000 
                   * exp(0 * logT - 0.50321666580471969 * -570 * invT);
        Corr  = 1.0;
        qf[69] *= Corr * k_f;
        // (68):  CH2(S) + CH4 <=> 2 CH3
        k_f = 1.0000000000000002e-06 * 16000000000000 
                   * exp(0 * logT - 0.50321666580471969 * -570 * invT);
        Corr  = 1.0;
        qf[70] *= Corr * k_f;
        // (69):  CH2(S) + CO <=> CH2 + CO
        k_f = 1.0000000000000002e-06 * 9000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        Corr  = 1.0;
        qf[71] *= Corr * k_f;
        // (70):  CH2(S) + CO2 <=> CH2 + CO2
        k_f = 1.0000000000000002e-06 * 7000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        Corr  = 1.0;
        qf[72] *= Corr * k_f;
        // (71):  CH2(S) + CO2 <=> CO + CH2O
        k_f = 1.0000000000000002e-06 * 14000000000000 
                   * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        Corr  = 1.0;
        qf[73] *= Corr * k_f;
        // (72):  CH3 + O2 <=> O + CH3O
        k_f = 1.0000000000000002e-06 * 26750000000000 
                   * exp(0 * logT - 0.50321666580471969 * 28800 * invT);
        Corr  = 1.0;
        qf[74] *= Corr * k_f;
        // (73):  CH3 + O2 <=> OH + CH2O
        qf[75] *= 1.0000000000000002e-06 * 36000000000 
                   * exp(0 * logT - 0.50321666580471969 * 8940 * invT);
        // (74):  2 CH3 (+M) <=> C2H6 (+M)
        qf[7] *= 1.0000000000000002e-06 * 21200000000000000 
                   * exp(-0.96999999999999997 * logT - 0.50321666580471969 * 620 * invT);
        // (75):  2 CH3 <=> H + C2H5
        qf[76] *= 1.0000000000000002e-06 * 4990000000000 
                   * exp(0.10000000000000001 * logT - 0.50321666580471969 * 10600 * invT);
        // (76):  CH3 + HCO <=> CH4 + CO
        qf[77] *= 1.0000000000000002e-06 * 26480000000000 
                   * exp(0 * logT - 0.50321666580471969 * 0 * invT);
        // (77):  CH3 + CH2O <=> HCO + CH4
        qf[78] *= 1.0000000000000002e-06 * 3320 
                   * exp(2.8100000000000001 * logT - 0.50321666580471969 * 5860 * invT);
        // (78):  CH3 + C2H6 <=> C2H5 + CH4
        qf[79] *= 1.0000000000000002e-06 * 6140000 
                   * exp(1.74 * logT - 0.50321666580471969 * 10450 * invT);
        // (79):  HCO + H2O <=> H + CO + H2O
        qf[80] *= 1.0000000000000002e-06 * 2.244e+18 
                   * exp(-1 * logT - 0.50321666580471969 * 17000 * invT);
        // (80):  HCO + M <=> H + CO + M
        qf[13] *= 1.0000000000000002e-06 * 1.87e+17 
                   * exp(-1 * logT - 0.50321666580471969 * 17000 * invT);
        // (81):  HCO + O2 <=> HO2 + CO
        qf[81] *= 1.0000000000000002e-06 * 7600000000000 
                   * exp(0 * logT - 0.50321666580471969 * 400 * invT);
        // (82):  CH3O + O2 <=> HO2 + CH2O
        qf[82] *= 1.0000000000000002e-06 * 4.2799999999999999e-13 
                   * exp(7.5999999999999996 * logT - 0.50321666580471969 * -3530 * invT);
        // (83):  C2H5 + O2 <=> HO2 + C2H4
        qf[83] *= 1.0000000000000002e-06 * 840000000000 
                   * exp(0 * logT - 0.50321666580471969 * 3875 * invT);

        for (int i = 0; i < 21; ++i) {
            wdot[i] = 0.0;
        }

        qdot = qf[0];
        wdot[1] -= qdot;
        wdot[7] -= qdot;

        qdot = qf[1];
        wdot[1] -= qdot;
        wdot[9] -= qdot;

        qdot = qf[2];
        wdot[1]  -= qdot;
        wdot[13] -= qdot;

        qdot = qf[3];
        wdot[1]  -= qdot;
        wdot[14] -= qdot;

        qdot = qf[4];
        wdot[1]  -= qdot;
        wdot[16] -= qdot;

        qdot = qf[5];
        wdot[1] -= qdot;
        wdot[17] -= qdot;

        qdot = qf[6];
        wdot[0] -= qdot;
        wdot[11] -= qdot;

        qdot = qf[7];
        wdot[9] -= 2 * qdot;

        qdot = qf[8];
        wdot[1] -= qdot;
        wdot[2] -= qdot;

        qdot = qf[9];
        wdot[2] -= qdot;
        wdot[11] -= qdot;

        qdot = qf[10];
        wdot[1] -= qdot;
        wdot[3] -= qdot;

        qdot = qf[11];
        wdot[1] -= 2 * qdot;

        qdot = qf[12];
        wdot[1] -= qdot;
        wdot[4] -= qdot;

        qdot = qf[13];
        wdot[13] -= qdot;

        qdot = qf[14];
        wdot[0] -= qdot;
        wdot[2] -= qdot;

        qdot = qf[15];
        wdot[2] -= qdot;
        wdot[6] -= qdot;

        qdot = qf[16];
        wdot[2] -= qdot;
        wdot[7] -= qdot;

        qdot = qf[17];
        wdot[2] -= qdot;
        wdot[8] -= qdot;

        qdot = qf[18];
        wdot[2] -= qdot;
        wdot[9] -= qdot;

        qdot = qf[19];
        wdot[2] -= qdot;
        wdot[10] -= qdot;

        qdot = qf[20];
        wdot[2] -= qdot;
        wdot[13] -= qdot;

        qdot = qf[21];
        wdot[2] -= qdot;
        wdot[13] -= qdot;

        qdot = qf[22];
        wdot[2] -= qdot;
        wdot[14] -= qdot;

        qdot = qf[23];
        wdot[2] -= qdot;
        wdot[16] -= qdot;

        qdot = qf[24];
        wdot[2] -= qdot;
        wdot[17] -= qdot;

        qdot = qf[25];
        wdot[2] -= qdot;
        wdot[18] -= qdot;

        qdot = qf[26];
        wdot[3] -= qdot;
        wdot[11] -= qdot;

        qdot = qf[27];
        wdot[3] -= qdot;
        wdot[14] -= qdot;

        qdot = qf[28];
        wdot[1] -= qdot;
        wdot[3] -= 2 * qdot;

        qdot = qf[29];
        wdot[1] -= qdot;
        wdot[3] -= qdot;
        wdot[5] -= qdot;

        qdot = qf[30];
        wdot[1] -= qdot;
        wdot[3] -= qdot;
        wdot[19] -= qdot;

        qdot = qf[31];
        wdot[1] -= qdot;
        wdot[3] -= qdot;
        wdot[20] -= qdot;

        qdot = qf[32];
        wdot[1] -= qdot;
        wdot[3] -= qdot;

        qdot = qf[33];
        wdot[0] -= qdot;
        wdot[1] -= 2 * qdot;

        qdot = qf[34];
        wdot[1] -= 2 * qdot;
        wdot[5] -= qdot;

        qdot = qf[35];
        wdot[1] -= 2 * qdot;
        wdot[12] -= qdot;

        qdot = qf[36];
        wdot[1] -= qdot;
        wdot[6] -= qdot;

        qdot = qf[37];
        wdot[1] -= qdot;
        wdot[6] -= qdot;

        qdot = qf[38];
        wdot[1] -= qdot;
        wdot[10] -= qdot;

        qdot = qf[39];
        wdot[1] -= qdot;
        wdot[13] -= qdot;

        qdot = qf[40];
        wdot[1] -= qdot;
        wdot[14] -= qdot;

        qdot = qf[41];
        wdot[1] -= qdot;
        wdot[15] -= qdot;

        qdot = qf[42];
        wdot[1] -= qdot;
        wdot[18] -= qdot;

        qdot = qf[43];
        wdot[0] -= qdot;
        wdot[4] -= qdot;

        qdot = qf[44];
        wdot[4] -= 2 * qdot;

        qdot = qf[45];
        wdot[4] -= qdot;
        wdot[6] -= qdot;

        qdot = qf[46];
        wdot[4] -= qdot;
        wdot[7] -= qdot;

        qdot = qf[47];
        wdot[4] -= qdot;
        wdot[8] -= qdot;

        qdot = qf[48];
        wdot[4] -= qdot;
        wdot[9] -= qdot;

        qdot = qf[49];
        wdot[4] -= qdot;
        wdot[9] -= qdot;

        qdot = qf[50];
        wdot[4] -= qdot;
        wdot[10] -= qdot;

        qdot = qf[51];
        wdot[4] -= qdot;
        wdot[11] -= qdot;

        qdot = qf[52];
        wdot[4] -= qdot;
        wdot[13] -= qdot;

        qdot = qf[53];
        wdot[4] -= qdot;
        wdot[14] -= qdot;

        qdot = qf[54];
        wdot[4] -= qdot;
        wdot[18] -= qdot;

        qdot = qf[55];
        wdot[6] -= qdot;
        wdot[7] -= qdot;

        qdot = qf[56];
        wdot[6] -= qdot;
        wdot[9] -= qdot;

        qdot = qf[57];
        wdot[6] -= qdot;
        wdot[9] -= qdot;

        qdot = qf[58];
        wdot[6] -= qdot;
        wdot[11] -= qdot;

        qdot = qf[59];
        wdot[3] -= qdot;
        wdot[7] -= qdot;

        qdot = qf[60];
        wdot[0] -= qdot;
        wdot[7] -= qdot;

        qdot = qf[61];
        wdot[7] -= qdot;
        wdot[9] -= qdot;

        qdot = qf[62];
        wdot[7] -= qdot;
        wdot[10] -= qdot;

        qdot = qf[63];
        wdot[8] -= qdot;
        wdot[19] -= qdot;

        qdot = qf[64];
        wdot[8] -= qdot;
        wdot[20] -= qdot;

        qdot = qf[65];
        wdot[3] -= qdot;
        wdot[8] -= qdot;

        qdot = qf[66];
        wdot[3] -= qdot;
        wdot[8] -= qdot;

        qdot = qf[67];
        wdot[0] -= qdot;
        wdot[8] -= qdot;

        qdot = qf[68];
        wdot[5] -= qdot;
        wdot[8] -= qdot;

        qdot = qf[69];
        wdot[8] -= qdot;
        wdot[9] -= qdot;

        qdot = qf[70];
        wdot[8] -= qdot;
        wdot[10] -= qdot;

        qdot = qf[71];
        wdot[8] -= qdot;
        wdot[11] -= qdot;

        qdot = qf[72];
        wdot[8] -= qdot;
        wdot[12] -= qdot;

        qdot = qf[73];
        wdot[8] -= qdot;
        wdot[12] -= qdot;

        qdot = qf[74];
        wdot[3] -= qdot;
        wdot[9] -= qdot;

        qdot = qf[75];
        wdot[3] -= qdot;
        wdot[9] -= qdot;

        qdot = qf[76];
        wdot[9] -= 2 * qdot;

        qdot = qf[77];
        wdot[9] -= qdot;
        wdot[13] -= qdot;

        qdot = qf[78];
        wdot[9] -= qdot;
        wdot[14] -= qdot;

        qdot = qf[79];
        wdot[9] -= qdot;
        wdot[18] -= qdot;

        qdot = qf[80];
        wdot[5] -= qdot;
        wdot[13] -= qdot;

        qdot = qf[81];
        wdot[3] -= qdot;
        wdot[13] -= qdot;

        qdot = qf[82];
        wdot[3] -= qdot;
        wdot[15] -= qdot;

        qdot = qf[83];
        wdot[3] -= qdot;
        wdot[17] -= qdot;

	for (int n = 0; n < NUM_SPECIES; n++) {
		cdots(i,j,k,n) = -wdot[n];
	}
}

#endif
