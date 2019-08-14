#include <kernel.H>

namespace mythermo
{
    AMREX_GPU_DEVICE_MANAGED double fwd_A[84], fwd_beta[84], fwd_Ea[84];
    AMREX_GPU_DEVICE_MANAGED double activation_units[84], prefactor_units[84], phase_units[84];
    std::vector<std::vector<int>> kiv(84), nuv(84);
}

using namespace mythermo;

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
}

static AMREX_GPU_DEVICE_MANAGED double imw[21] = {
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

/*convert y[species] (mass fracs) to c[species] (molar conc) */
AMREX_GPU_HOST_DEVICE void molarConc(double *  rho, double *  T, double *  y,  int strideY, double *  c)
{
    for (int i = 0; i < 21; i++)
    {
        c[i] = (*rho)  * y[i*strideY] * imw[i];
    }
}

AMREX_GPU_DEVICE Real Q_reac_d(Real rho,
                               Real temp,
                               Real * Y, int strideY,
                               int idx)
{
  Real RcInv = 0.503217;
  Real k_f = fwd_A[idx] * exp(fwd_beta[idx] * temp - activation_units[idx]*RcInv * fwd_Ea[idx] / temp);
  Real C[21];
  molarConc(&rho,&temp,Y,strideY,C);
  Real Q = C[2] * C[8] * k_f;
  return Q;
}


