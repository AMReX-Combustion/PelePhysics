#include "mechanism.H"
const int rmap[134] = {
  6,   19,  24,  43,  66,  67,  106, 107, 20,  3,   5,   97,  108, 120, 131,
  132, 133, 0,   1,   2,   4,   7,   8,   9,   10,  11,  12,  13,  14,  15,
  16,  17,  18,  21,  22,  23,  25,  26,  27,  28,  29,  30,  31,  32,  33,
  34,  35,  36,  37,  38,  39,  40,  41,  42,  44,  45,  46,  47,  48,  49,
  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,  64,
  65,  68,  69,  70,  71,  72,  73,  74,  75,  76,  77,  78,  79,  80,  81,
  82,  83,  84,  85,  86,  87,  88,  89,  90,  91,  92,  93,  94,  95,  96,
  98,  99,  100, 101, 102, 103, 104, 105, 109, 110, 111, 112, 113, 114, 115,
  116, 117, 118, 119, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130};

// Returns 0-based map of reaction order
void
GET_RMAP(int* _rmap)
{
  for (int j = 0; j < 134; ++j) {
    _rmap[j] = rmap[j];
  }
}

// Returns a count of species in a reaction, and their indices
// and stoichiometric coefficients. (Eq 50)
void
CKINU(const int i, int& nspec, int ki[], int nu[])
{
  const int ns[134] = {
    4, 4, 4, 3, 3, 3, 2, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 3, 3, 3, 3, 4, 4,
    4, 3, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 5, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 3, 4,
    4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 3, 4, 5, 3, 3, 3, 4, 4, 4, 4, 4, 4,
    3, 4, 4, 5, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3};
  const int kiv[670] = {
    1,  4,  2,  6,  0, 1,  6,  2,  5,  0, 2,  3,  4,  6,  0, 2,  6,  5,  0,  0,
    5,  4,  6,  0,  0, 2,  4,  6,  0,  0, 7,  6,  0,  0,  0, 2,  7,  5,  6,  0,
    2,  7,  1,  8,  0, 7,  4,  8,  6,  0, 7,  6,  5,  8,  0, 7,  6,  5,  8,  0,
    2,  8,  6,  0,  0, 2,  8,  1,  3,  0, 8,  4,  3,  6,  0, 8,  6,  5,  3,  0,
    8,  6,  5,  3,  0, 8,  7,  3,  0,  0, 8,  7,  3,  0,  0, 2,  3,  8,  0,  0,
    9,  4,  10, 0,  0, 9,  6,  10, 2,  0, 9,  6,  10, 2,  0, 9,  8,  10, 6,  0,
    12, 2,  11, 0,  0, 11, 2,  12, 1,  0, 11, 4,  12, 6,  0, 11, 6,  12, 5,  0,
    11, 8,  12, 7,  0, 12, 8,  11, 3,  0, 13, 11, 12, 0,  0, 14, 0,  13, 0,  0,
    14, 5,  13, 5,  0, 14, 9,  13, 9,  0, 14, 10, 13, 10, 0, 14, 3,  9,  2,  6,
    14, 3,  9,  5,  0, 14, 4,  9,  1,  0, 14, 4,  2,  21, 0, 14, 1,  12, 2,  0,
    14, 2,  16, 1,  0, 14, 6,  20, 2,  0, 14, 10, 20, 9,  0, 13, 2,  12, 0,  0,
    13, 3,  21, 6,  0, 13, 3,  10, 2,  0, 13, 4,  9,  2,  0, 13, 2,  16, 1,  0,
    13, 6,  16, 5,  0, 16, 3,  21, 4,  0, 16, 4,  9,  2,  0, 16, 2,  15, 1,  0,
    16, 6,  2,  21, 0, 16, 5,  20, 2,  0, 16, 10, 9,  21, 0, 15, 6,  9,  2,  0,
    15, 3,  9,  4,  0, 12, 3,  18, 4,  0, 12, 3,  20, 6,  0, 12, 4,  20, 2,  0,
    12, 6,  14, 5,  0, 12, 6,  20, 1,  0, 12, 6,  19, 2,  0, 12, 6,  18, 2,  0,
    12, 6,  13, 5,  0, 12, 8,  18, 6,  0, 17, 12, 6,  0,  0, 17, 14, 5,  0,  0,
    17, 2,  18, 1,  0, 17, 2,  19, 1,  0, 17, 4,  18, 6,  0, 17, 4,  19, 6,  0,
    17, 6,  18, 5,  0, 17, 6,  19, 5,  0, 17, 3,  18, 8,  0, 17, 8,  18, 7,  0,
    17, 8,  19, 7,  0, 12, 17, 19, 11, 0, 12, 17, 18, 11, 0, 19, 3,  20, 8,  0,
    19, 3,  20, 8,  0, 19, 2,  20, 1,  0, 19, 8,  20, 7,  0, 19, 21, 20, 0,  0,
    19, 6,  20, 5,  0, 19, 4,  20, 6,  0, 18, 3,  20, 8,  0, 18, 2,  20, 1,  0,
    18, 8,  20, 7,  0, 12, 18, 20, 11, 0, 18, 20, 17, 0,  0, 20, 4,  21, 6,  0,
    20, 2,  1,  21, 0, 20, 6,  5,  21, 0, 20, 8,  7,  21, 0, 20, 12, 11, 21, 0,
    20, 18, 17, 21, 0, 21, 9,  2,  0,  0, 21, 3,  9,  8,  0, 21, 4,  9,  6,  0,
    2,  21, 9,  1,  0, 21, 6,  9,  5,  0, 12, 21, 11, 9,  0, 21, 20, 9,  0,  0,
    21, 4,  10, 2,  0, 21, 8,  10, 2,  6, 20, 2,  19, 0,  0, 18, 20, 2,  0,  0,
    4,  6,  8,  0,  0, 16, 6,  15, 5,  0, 13, 21, 12, 9,  0, 13, 4,  9,  1,  0,
    13, 6,  20, 2,  0, 13, 10, 20, 9,  0, 14, 4,  9,  2,  0, 14, 11, 12, 0,  0,
    17, 2,  12, 5,  0, 20, 4,  10, 2,  0, 20, 6,  10, 2,  1, 18, 9,  12, 10, 0,
    12, 4,  18, 0,  0, 21, 8,  9,  7,  0, 18, 6,  20, 5,  0, 19, 12, 20, 11, 0,
    16, 4,  22, 25, 0, 22, 5,  9,  23, 0, 25, 23, 2,  5,  0, 25, 23, 2,  6,  0,
    25, 23, 1,  6,  0, 1,  24, 25, 7,  0, 2,  24, 25, 8,  0, 25, 3,  24, 0,  0,
    25, 3,  24, 0,  0, 25, 3,  24, 0,  0};
  const int nuv[670] = {
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0,
    -1, -1, 2, 0, 0, -1, -1, 1, 0, 0, -1, 2,  0, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 2, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -2, 1,  1, 0, 0, -2, 1,  1, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 2, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 1,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 2, 0, -1, -1, 1, 2, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, 1,  1, 0, 0, -1, 1,  1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 2, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -2, 1,  1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -2, 1,  1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 1, -1, -1, 1, 0, 0, -1, 1,  1, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 2, 0, -1, -1, 2, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 2, 0, -1, -1, 1, 1, 1, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 2, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 0, 0};
  if (i < 1) {
    // Return max num species per reaction
    nspec = 5;
  } else {
    if (i > 134) {
      nspec = -1;
    } else {
      nspec = ns[i - 1];
      for (int j = 0; j < nspec; ++j) {
        ki[j] = kiv[(i - 1) * 5 + j] + 1;
        nu[j] = nuv[(i - 1) * 5 + j];
      }
    }
  }
}

// Returns the progress rates of each reactions
// Given P, T, and mole fractions
void
CKKFKR(
  const amrex::Real P,
  const amrex::Real T,
  const amrex::Real x[],
  amrex::Real q_f[],
  amrex::Real q_r[])
{
  amrex::Real c[26]; // temporary storage
  amrex::Real PORT =
    1e6 * P / (8.31446261815324e+07 * T); // 1e6 * P/RT so c goes to SI units

  // Compute conversion, see Eq 10
  for (int id = 0; id < 26; ++id) {
    c[id] = x[id] * PORT;
  }

  // convert to chemkin units
  progressRateFR(q_f, q_r, c, T);

  // convert to chemkin units
  for (int id = 0; id < 134; ++id) {
    q_f[id] *= 1.0e-6;
    q_r[id] *= 1.0e-6;
  }
}

// compute the progress rate for each reaction
// USES progressRate : todo switch to GPU
void
progressRateFR(
  amrex::Real* q_f, amrex::Real* q_r, amrex::Real* sc, amrex::Real T)
{
  const amrex::Real invT = 1.0 / T;
  const amrex::Real logT = log(T);
  // compute the Gibbs free energy
  amrex::Real g_RT[26];
  gibbs(g_RT, T);

  amrex::Real sc_qss[1];
  comp_qfqr(q_f, q_r, sc, sc_qss, T, invT, logT);
}

// save atomic weights into array
void
atomicWeight(amrex::Real* awt)
{
  awt[0] = 14.007000; // N
  awt[1] = 1.008000;  // H
  awt[2] = 15.999000; // O
  awt[3] = 12.011000; // C
  awt[4] = 0.000549;  // E
}

// get atomic weight for all elements
void
CKAWT(amrex::Real* awt)
{
  atomicWeight(awt);
}

// Returns the elemental composition
// of the speciesi (mdim is num of elements)
void
CKNCF(int* ncf)
{
  int kd = 5;
  // Zero ncf
  for (int id = 0; id < kd * 26; ++id) {
    ncf[id] = 0;
  }

  // N2
  ncf[0 * kd + 0] = 2; // N

  // H2
  ncf[1 * kd + 1] = 2; // H

  // H
  ncf[2 * kd + 1] = 1; // H

  // O2
  ncf[3 * kd + 2] = 2; // O

  // O
  ncf[4 * kd + 2] = 1; // O

  // H2O
  ncf[5 * kd + 1] = 2; // H
  ncf[5 * kd + 2] = 1; // O

  // OH
  ncf[6 * kd + 1] = 1; // H
  ncf[6 * kd + 2] = 1; // O

  // H2O2
  ncf[7 * kd + 1] = 2; // H
  ncf[7 * kd + 2] = 2; // O

  // HO2
  ncf[8 * kd + 1] = 1; // H
  ncf[8 * kd + 2] = 2; // O

  // CO
  ncf[9 * kd + 3] = 1; // C
  ncf[9 * kd + 2] = 1; // O

  // CO2
  ncf[10 * kd + 3] = 1; // C
  ncf[10 * kd + 2] = 2; // O

  // CH4
  ncf[11 * kd + 3] = 1; // C
  ncf[11 * kd + 1] = 4; // H

  // CH3
  ncf[12 * kd + 3] = 1; // C
  ncf[12 * kd + 1] = 3; // H

  // CH2
  ncf[13 * kd + 3] = 1; // C
  ncf[13 * kd + 1] = 2; // H

  // CH2(S)
  ncf[14 * kd + 3] = 1; // C
  ncf[14 * kd + 1] = 2; // H

  // C
  ncf[15 * kd + 3] = 1; // C

  // CH
  ncf[16 * kd + 3] = 1; // C
  ncf[16 * kd + 1] = 1; // H

  // CH3OH
  ncf[17 * kd + 3] = 1; // C
  ncf[17 * kd + 1] = 4; // H
  ncf[17 * kd + 2] = 1; // O

  // CH3O
  ncf[18 * kd + 3] = 1; // C
  ncf[18 * kd + 1] = 3; // H
  ncf[18 * kd + 2] = 1; // O

  // CH2OH
  ncf[19 * kd + 3] = 1; // C
  ncf[19 * kd + 1] = 3; // H
  ncf[19 * kd + 2] = 1; // O

  // CH2O
  ncf[20 * kd + 3] = 1; // C
  ncf[20 * kd + 1] = 2; // H
  ncf[20 * kd + 2] = 1; // O

  // HCO
  ncf[21 * kd + 3] = 1; // C
  ncf[21 * kd + 1] = 1; // H
  ncf[21 * kd + 2] = 1; // O

  // CHO+
  ncf[22 * kd + 3] = 1;  // C
  ncf[22 * kd + 4] = -1; // E
  ncf[22 * kd + 1] = 1;  // H
  ncf[22 * kd + 2] = 1;  // O

  // H3O+
  ncf[23 * kd + 4] = -1; // E
  ncf[23 * kd + 1] = 3;  // H
  ncf[23 * kd + 2] = 1;  // O

  // O2-
  ncf[24 * kd + 4] = 1; // E
  ncf[24 * kd + 2] = 2; // O

  // E
  ncf[25 * kd + 4] = 1; // E
}

// Returns the vector of strings of element names
void
CKSYME_STR(amrex::Vector<std::string>& ename)
{
  ename.resize(5);
  ename[0] = "N";
  ename[1] = "H";
  ename[2] = "O";
  ename[3] = "C";
  ename[4] = "E";
}

// Returns the vector of strings of species names
void
CKSYMS_STR(amrex::Vector<std::string>& kname)
{
  kname.resize(26);
  kname[0] = "N2";
  kname[1] = "H2";
  kname[2] = "H";
  kname[3] = "O2";
  kname[4] = "O";
  kname[5] = "H2O";
  kname[6] = "OH";
  kname[7] = "H2O2";
  kname[8] = "HO2";
  kname[9] = "CO";
  kname[10] = "CO2";
  kname[11] = "CH4";
  kname[12] = "CH3";
  kname[13] = "CH2";
  kname[14] = "CH2(S)";
  kname[15] = "C";
  kname[16] = "CH";
  kname[17] = "CH3OH";
  kname[18] = "CH3O";
  kname[19] = "CH2OH";
  kname[20] = "CH2O";
  kname[21] = "HCO";
  kname[22] = "CHO+";
  kname[23] = "H3O+";
  kname[24] = "O2-";
  kname[25] = "E";
}

// compute the sparsity pattern of the chemistry Jacobian
void
SPARSITY_INFO(int* nJdata, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 729> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 26> conc = {0.0};
  for (int n = 0; n < 26; n++) {
    conc[n] = 1.0 / 26.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 27; k++) {
    for (int l = 0; l < 27; l++) {
      if (Jac[27 * k + l] != 0.0) {
        nJdata_tmp = nJdata_tmp + 1;
      }
    }
  }

  *nJdata = NCELLS * nJdata_tmp;
}

// compute the sparsity pattern of the system Jacobian
void
SPARSITY_INFO_SYST(int* nJdata, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 729> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 26> conc = {0.0};
  for (int n = 0; n < 26; n++) {
    conc[n] = 1.0 / 26.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 27; k++) {
    for (int l = 0; l < 27; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[27 * k + l] != 0.0) {
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
    }
  }

  *nJdata = NCELLS * nJdata_tmp;
}

// compute the sparsity pattern of the simplified (for preconditioning) system
// Jacobian
void
SPARSITY_INFO_SYST_SIMPLIFIED(int* nJdata, const int* consP)
{
  amrex::GpuArray<amrex::Real, 729> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 26> conc = {0.0};
  for (int n = 0; n < 26; n++) {
    conc[n] = 1.0 / 26.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 27; k++) {
    for (int l = 0; l < 27; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[27 * k + l] != 0.0) {
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
    }
  }

  nJdata[0] = nJdata_tmp;
}

// compute the sparsity pattern of the chemistry Jacobian in CSC format -- base
// 0
void
SPARSITY_PREPROC_CSC(int* rowVals, int* colPtrs, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 729> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 26> conc = {0.0};
  for (int n = 0; n < 26; n++) {
    conc[n] = 1.0 / 26.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int nc = 0; nc < NCELLS; nc++) {
    int offset_row = nc * 27;
    int offset_col = nc * 27;
    for (int k = 0; k < 27; k++) {
      for (int l = 0; l < 27; l++) {
        if (Jac[27 * k + l] != 0.0) {
          rowVals[nJdata_tmp] = l + offset_row;
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
      colPtrs[offset_col + (k + 1)] = nJdata_tmp;
    }
  }
}

// compute the sparsity pattern of the chemistry Jacobian in CSR format -- base
// 0
void
SPARSITY_PREPROC_CSR(
  int* colVals, int* rowPtrs, const int* consP, int NCELLS, int base)
{
  amrex::GpuArray<amrex::Real, 729> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 26> conc = {0.0};
  for (int n = 0; n < 26; n++) {
    conc[n] = 1.0 / 26.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtrs[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 27;
      for (int l = 0; l < 27; l++) {
        for (int k = 0; k < 27; k++) {
          if (Jac[27 * k + l] != 0.0) {
            colVals[nJdata_tmp - 1] = k + 1 + offset;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
        rowPtrs[offset + (l + 1)] = nJdata_tmp;
      }
    }
  } else {
    rowPtrs[0] = 0;
    int nJdata_tmp = 0;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 27;
      for (int l = 0; l < 27; l++) {
        for (int k = 0; k < 27; k++) {
          if (Jac[27 * k + l] != 0.0) {
            colVals[nJdata_tmp] = k + offset;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
        rowPtrs[offset + (l + 1)] = nJdata_tmp;
      }
    }
  }
}

// compute the sparsity pattern of the system Jacobian
// CSR format BASE is user choice
void
SPARSITY_PREPROC_SYST_CSR(
  int* colVals, int* rowPtr, const int* consP, int NCELLS, int base)
{
  amrex::GpuArray<amrex::Real, 729> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 26> conc = {0.0};
  for (int n = 0; n < 26; n++) {
    conc[n] = 1.0 / 26.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 27;
      for (int l = 0; l < 27; l++) {
        for (int k = 0; k < 27; k++) {
          if (k == l) {
            colVals[nJdata_tmp - 1] = l + 1 + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[27 * k + l] != 0.0) {
              colVals[nJdata_tmp - 1] = k + 1 + offset;
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
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 27;
      for (int l = 0; l < 27; l++) {
        for (int k = 0; k < 27; k++) {
          if (k == l) {
            colVals[nJdata_tmp] = l + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[27 * k + l] != 0.0) {
              colVals[nJdata_tmp] = k + offset;
              nJdata_tmp = nJdata_tmp + 1;
            }
          }
        }
        rowPtr[offset + (l + 1)] = nJdata_tmp;
      }
    }
  }
}

// compute the sparsity pattern of the simplified (for precond) system Jacobian
// on CPU BASE 0
void
SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(
  int* rowVals, int* colPtrs, int* indx, const int* consP)
{
  amrex::GpuArray<amrex::Real, 729> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 26> conc = {0.0};
  for (int n = 0; n < 26; n++) {
    conc[n] = 1.0 / 26.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int k = 0; k < 27; k++) {
    for (int l = 0; l < 27; l++) {
      if (k == l) {
        rowVals[nJdata_tmp] = l;
        indx[nJdata_tmp] = 27 * k + l;
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[27 * k + l] != 0.0) {
          rowVals[nJdata_tmp] = l;
          indx[nJdata_tmp] = 27 * k + l;
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
    }
    colPtrs[k + 1] = nJdata_tmp;
  }
}

// compute the sparsity pattern of the simplified (for precond) system Jacobian
// CSR format BASE is under choice
void
SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(
  int* colVals, int* rowPtr, const int* consP, int base)
{
  amrex::GpuArray<amrex::Real, 729> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 26> conc = {0.0};
  for (int n = 0; n < 26; n++) {
    conc[n] = 1.0 / 26.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int l = 0; l < 27; l++) {
      for (int k = 0; k < 27; k++) {
        if (k == l) {
          colVals[nJdata_tmp - 1] = l + 1;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[27 * k + l] != 0.0) {
            colVals[nJdata_tmp - 1] = k + 1;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
      }
      rowPtr[l + 1] = nJdata_tmp;
    }
  } else {
    rowPtr[0] = 0;
    int nJdata_tmp = 0;
    for (int l = 0; l < 27; l++) {
      for (int k = 0; k < 27; k++) {
        if (k == l) {
          colVals[nJdata_tmp] = l;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[27 * k + l] != 0.0) {
            colVals[nJdata_tmp] = k;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
      }
      rowPtr[l + 1] = nJdata_tmp;
    }
  }
}
