#include "mechanism.H"
const int rmap[175] = {
  8,   15,  47,  48,  85,  92,  93,  94,  108, 111, 21,  4,   5,   6,   7,
  25,  35,  36,  54,  63,  118, 165, 166, 169, 170, 0,   1,   2,   3,   9,
  10,  11,  12,  13,  14,  16,  17,  18,  19,  20,  22,  23,  24,  26,  27,
  28,  29,  30,  31,  32,  33,  34,  37,  38,  39,  40,  41,  42,  43,  44,
  45,  46,  49,  50,  51,  52,  53,  55,  56,  57,  58,  59,  60,  61,  62,
  64,  65,  66,  67,  68,  69,  70,  71,  72,  73,  74,  75,  76,  77,  78,
  79,  80,  81,  82,  83,  84,  86,  87,  88,  89,  90,  91,  95,  96,  97,
  98,  99,  100, 101, 102, 103, 104, 105, 106, 107, 109, 110, 112, 113, 114,
  115, 116, 117, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130,
  131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145,
  146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160,
  161, 162, 163, 164, 167, 168, 171, 172, 173, 174};

// Returns 0-based map of reaction order
void
GET_RMAP(int* _rmap)
{
  for (int j = 0; j < 175; ++j) {
    _rmap[j] = rmap[j];
  }
}

// Returns a count of species in a reaction, and their indices
// and stoichiometric coefficients. (Eq 50)
void
CKINU(const int i, int& nspec, int ki[], int nu[])
{
  const int ns[175] = {
    4, 4, 4, 3, 2, 2, 3, 3, 3, 4, 3, 4, 4, 3, 3, 2, 4, 4, 4, 4, 4, 3, 4, 4, 4,
    3, 4, 4, 4, 4, 4, 5, 3, 4, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 2, 3, 4,
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 3, 4, 4, 3, 3, 3, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 3, 4, 4, 3, 4, 4, 3, 4, 4, 4, 4, 4, 3, 2, 4, 4, 4, 4, 4, 4,
    4, 5, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 3, 3,
    3, 3, 4, 3, 4, 2, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 5, 3, 3, 5, 5, 5, 4};
  const int kiv[875] = {
    0,  19, 5,  7,  0, 1,  5,  0,  7,  0, 1,  7,  0,  8,  0, 8,  5,  7,  0,  0,
    1,  0,  0,  0,  0, 5,  19, 0,  0,  0, 0,  5,  7,  0,  0, 0,  7,  8,  0,  0,
    0,  19, 20, 0,  0, 0,  20, 1,  19, 0, 0,  20, 7,  0,  0, 20, 5,  19, 7,  0,
    20, 7,  8,  19, 0, 20, 21, 19, 0,  0, 20, 21, 19, 0,  0, 21, 7,  0,  0,  0,
    0,  21, 8,  7,  0, 0,  21, 1,  20, 0, 21, 5,  20, 7,  0, 21, 7,  8,  20, 0,
    21, 7,  8,  20, 0, 11, 5,  22, 0,  0, 11, 19, 22, 5,  0, 11, 20, 22, 7,  0,
    11, 7,  22, 0,  0, 13, 11, 0,  0,  0, 13, 19, 11, 20, 0, 0,  13, 11, 1,  0,
    13, 5,  11, 7,  0, 13, 7,  11, 8,  0, 13, 5,  22, 0,  0, 13, 20, 22, 0,  7,
    13, 11, 1,  0,  0, 4,  13, 6,  11, 0, 13, 15, 11, 0,  0, 15, 0,  13, 0,  0,
    15, 11, 1,  0,  0, 15, 0,  1,  13, 0, 15, 5,  13, 7,  0, 15, 7,  8,  13, 0,
    15, 19, 13, 20, 0, 15, 20, 21, 13, 0, 15, 4,  6,  13, 0, 4,  5,  15, 0,  0,
    4,  19, 18, 5,  0, 4,  19, 15, 7,  0, 4,  20, 18, 7,  0, 4,  16, 0,  0,  0,
    4,  0,  6,  0,  0, 6,  0,  4,  1,  0, 6,  5,  4,  7,  0, 6,  7,  4,  8,  0,
    4,  20, 6,  19, 0, 6,  20, 4,  21, 0, 17, 15, 0,  0,  0, 17, 0,  15, 1,  0,
    17, 0,  4,  7,  0, 17, 5,  15, 7,  0, 17, 7,  15, 8,  0, 17, 19, 15, 20, 0,
    17, 19, 15, 20, 0, 17, 20, 15, 21, 0, 17, 13, 15, 0,  0, 18, 15, 0,  0,  0,
    18, 0,  4,  7,  0, 18, 5,  15, 7,  0, 18, 7,  15, 8,  0, 18, 19, 15, 20, 0,
    18, 19, 15, 20, 0, 18, 20, 15, 21, 0, 18, 11, 4,  22, 0, 4,  14, 0,  0,  0,
    2,  6,  4,  0,  0, 3,  6,  4,  0,  0, 4,  7,  2,  8,  0, 4,  7,  3,  8,  0,
    2,  4,  12, 0,  0, 3,  4,  12, 0,  0, 18, 0,  3,  8,  0, 16, 0,  14, 1,  0,
    16, 5,  14, 7,  0, 16, 7,  14, 8,  0, 16, 19, 14, 20, 0, 16, 20, 14, 21, 0,
    16, 4,  14, 6,  0, 14, 0,  16, 0,  0, 14, 0,  12, 1,  0, 14, 5,  15, 4,  0,
    14, 19, 12, 20, 0, 14, 12, 16, 0,  0, 14, 13, 16, 11, 0, 14, 5,  23, 0,  0,
    12, 9,  1,  0,  0, 12, 0,  14, 0,  0, 10, 0,  12, 0,  0, 12, 0,  10, 1,  0,
    12, 7,  10, 8,  0, 12, 4,  10, 6,  0, 12, 5,  4,  13, 0, 10, 7,  9,  8,  0,
    12, 5,  10, 7,  0, 12, 19, 10, 20, 0, 10, 0,  9,  1,  0, 10, 21, 12, 20, 0,
    10, 4,  9,  6,  0, 10, 9,  12, 0,  0, 10, 19, 15, 13, 0, 10, 19, 9,  20, 0,
    9,  0,  10, 0,  0, 9,  5,  2,  11, 0, 9,  7,  4,  11, 0, 2,  0,  4,  0,  0,
    2,  5,  0,  13, 0, 2,  7,  15, 0,  0, 2,  1,  4,  0,  0, 2,  19, 13, 7,  0,
    2,  20, 15, 7,  0, 2,  9,  1,  0,  0, 3,  2,  0,  0,  0, 3,  8,  2,  8,  0,
    3,  11, 2,  11, 0, 3,  22, 2,  22, 0, 3,  5,  11, 1,  0, 3,  5,  0,  13, 0,
    3,  7,  15, 0,  0, 3,  1,  4,  0,  0, 3,  19, 11, 0,  7, 3,  19, 11, 8,  0,
    3,  22, 15, 11, 0, 23, 4,  13, 0,  0, 26, 4,  18, 0,  0, 26, 7,  24, 8,  0,
    26, 0,  24, 1,  0, 4,  26, 24, 6,  0, 26, 5,  24, 7,  0, 26, 20, 24, 21, 0,
    26, 19, 24, 20, 0, 24, 15, 4,  0,  0, 18, 24, 15, 26, 0, 15, 24, 26, 13, 0,
    24, 20, 30, 7,  0, 30, 29, 0,  0,  0, 29, 19, 28, 20, 0, 29, 7,  28, 8,  0,
    29, 20, 28, 21, 0, 29, 5,  28, 7,  0, 29, 0,  28, 1,  0, 4,  29, 28, 6,  0,
    28, 18, 11, 0,  0, 28, 4,  22, 0,  0, 24, 19, 34, 0,  0, 34, 30, 19, 0,  0,
    34, 31, 29, 19, 0, 30, 15, 18, 0,  0, 30, 19, 29, 20, 0, 34, 35, 0,  0,  0,
    35, 15, 7,  0,  0, 35, 19, 37, 0,  0, 37, 36, 7,  0,  0, 36, 32, 7,  0,  0,
    32, 33, 0,  0,  0, 33, 11, 27, 0,  0, 33, 17, 22, 0,  0, 27, 0,  25, 0,  0,
    15, 7,  27, 0,  0, 25, 11, 8,  0,  0, 25, 22, 1,  0,  0, 25, 13, 7,  0,  0,
    25, 7,  22, 0,  8, 25, 11, 8,  0,  0, 25, 22, 1,  0,  0, 0,  25, 11, 1,  7,
    4,  25, 6,  11, 7, 25, 20, 11, 21, 7, 25, 5,  11, 7,  0};
  const int nuv[875] = {
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 2, 0, 0,
    -1, 2,  0, 0, 0, -2, 1,  0, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 2, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -2, 1,  1, 0, 0, -2, 1,  1, 0, 0, -1, 2,  0, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 1,
    -2, 2,  1, 0, 0, -1, -1, 1, 1, 0, -2, 1,  1, 0, 0, -1, 1,  1, 0, 0,
    -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -2, 1,  0, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 2, 0, 0, -1, 1,  1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -2, 1,  1, 0, 0,
    -1, -1, 2, 0, 0, -1, -1, 2, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -2, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, 1,  1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -2, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -2, 1,  1, 0, 0, -1, 1,  0, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 1, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, 1,  1, 0, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, 1,  1, 0, 0, -1, 1,  1, 0, 0, -1, -1, 1, 0, 0, -2, 2,  1, 0, 0,
    -2, 1,  1, 1, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, 1,  0, 0, 0,
    -1, 2,  1, 0, 0, -1, -1, 1, 0, 0, -1, 1,  1, 0, 0, -1, 1,  1, 0, 0,
    -1, 1,  0, 0, 0, -1, 1,  1, 0, 0, -1, 1,  1, 0, 0, -1, 1,  1, 0, 0,
    -1, -1, 1, 0, 0, -1, 1,  1, 0, 0, -1, 1,  1, 0, 0, -1, 1,  1, 0, 0,
    -1, -1, 1, 1, 1, -1, 1,  1, 0, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 1,
    -1, -1, 1, 1, 1, -1, -1, 1, 1, 1, -1, -1, 1, 2, 0};
  if (i < 1) {
    // Return max num species per reaction
    nspec = 5;
  } else {
    if (i > 175) {
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
  amrex::Real c[39]; // temporary storage
  amrex::Real PORT =
    1e6 * P / (8.31446261815324e+07 * T); // 1e6 * P/RT so c goes to SI units

  // Compute conversion, see Eq 10
  for (int id = 0; id < 39; ++id) {
    c[id] = x[id] * PORT;
  }

  // convert to chemkin units
  progressRateFR(q_f, q_r, c, T);

  // convert to chemkin units
  for (int id = 0; id < 175; ++id) {
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
  const amrex::Real tc[5] = {
    log(T), T, T * T, T * T * T, T * T * T * T}; // temperature cache
  amrex::Real invT = 1.0 / tc[1];
  // compute the Gibbs free energy
  amrex::Real g_RT[39];
  gibbs(g_RT, tc);

  amrex::Real sc_qss[1];
  comp_qfqr(q_f, q_r, sc, sc_qss, tc, invT);
}

// save atomic weights into array
void
atomicWeight(amrex::Real* awt)
{
  awt[0] = 12.011000; // C
  awt[1] = 1.008000;  // H
  awt[2] = 15.999000; // O
  awt[3] = 14.007000; // N
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
  int kd = 4;
  // Zero ncf
  for (int id = 0; id < kd * 39; ++id) {
    ncf[id] = 0;
  }

  // H
  ncf[0 * kd + 1] = 1; // H

  // H2
  ncf[1 * kd + 1] = 2; // H

  // CH2
  ncf[2 * kd + 0] = 1; // C
  ncf[2 * kd + 1] = 2; // H

  // CH2(S)
  ncf[3 * kd + 0] = 1; // C
  ncf[3 * kd + 1] = 2; // H

  // CH3
  ncf[4 * kd + 0] = 1; // C
  ncf[4 * kd + 1] = 3; // H

  // O
  ncf[5 * kd + 2] = 1; // O

  // CH4
  ncf[6 * kd + 0] = 1; // C
  ncf[6 * kd + 1] = 4; // H

  // OH
  ncf[7 * kd + 1] = 1; // H
  ncf[7 * kd + 2] = 1; // O

  // H2O
  ncf[8 * kd + 1] = 2; // H
  ncf[8 * kd + 2] = 1; // O

  // C2H2
  ncf[9 * kd + 0] = 2; // C
  ncf[9 * kd + 1] = 2; // H

  // C2H3
  ncf[10 * kd + 0] = 2; // C
  ncf[10 * kd + 1] = 3; // H

  // CO
  ncf[11 * kd + 0] = 1; // C
  ncf[11 * kd + 2] = 1; // O

  // C2H4
  ncf[12 * kd + 0] = 2; // C
  ncf[12 * kd + 1] = 4; // H

  // HCO
  ncf[13 * kd + 0] = 1; // C
  ncf[13 * kd + 1] = 1; // H
  ncf[13 * kd + 2] = 1; // O

  // C2H5
  ncf[14 * kd + 0] = 2; // C
  ncf[14 * kd + 1] = 5; // H

  // CH2O
  ncf[15 * kd + 0] = 1; // C
  ncf[15 * kd + 1] = 2; // H
  ncf[15 * kd + 2] = 1; // O

  // C2H6
  ncf[16 * kd + 0] = 2; // C
  ncf[16 * kd + 1] = 6; // H

  // CH2OH
  ncf[17 * kd + 0] = 1; // C
  ncf[17 * kd + 1] = 3; // H
  ncf[17 * kd + 2] = 1; // O

  // CH3O
  ncf[18 * kd + 0] = 1; // C
  ncf[18 * kd + 1] = 3; // H
  ncf[18 * kd + 2] = 1; // O

  // O2
  ncf[19 * kd + 2] = 2; // O

  // HO2
  ncf[20 * kd + 1] = 1; // H
  ncf[20 * kd + 2] = 2; // O

  // H2O2
  ncf[21 * kd + 1] = 2; // H
  ncf[21 * kd + 2] = 2; // O

  // CO2
  ncf[22 * kd + 0] = 1; // C
  ncf[22 * kd + 2] = 2; // O

  // CH3HCO
  ncf[23 * kd + 0] = 2; // C
  ncf[23 * kd + 1] = 4; // H
  ncf[23 * kd + 2] = 1; // O

  // CH3OCH2
  ncf[24 * kd + 0] = 2; // C
  ncf[24 * kd + 1] = 5; // H
  ncf[24 * kd + 2] = 1; // O

  // HCOOH
  ncf[25 * kd + 0] = 1; // C
  ncf[25 * kd + 1] = 2; // H
  ncf[25 * kd + 2] = 2; // O

  // CH3OCH3
  ncf[26 * kd + 0] = 2; // C
  ncf[26 * kd + 1] = 6; // H
  ncf[26 * kd + 2] = 1; // O

  // HOCH2O
  ncf[27 * kd + 0] = 1; // C
  ncf[27 * kd + 1] = 3; // H
  ncf[27 * kd + 2] = 2; // O

  // CH3OCO
  ncf[28 * kd + 0] = 2; // C
  ncf[28 * kd + 1] = 3; // H
  ncf[28 * kd + 2] = 2; // O

  // CH3OCHO
  ncf[29 * kd + 0] = 2; // C
  ncf[29 * kd + 1] = 4; // H
  ncf[29 * kd + 2] = 2; // O

  // CH3OCH2O
  ncf[30 * kd + 0] = 2; // C
  ncf[30 * kd + 1] = 5; // H
  ncf[30 * kd + 2] = 2; // O

  // CH3OCH2OH
  ncf[31 * kd + 0] = 2; // C
  ncf[31 * kd + 1] = 6; // H
  ncf[31 * kd + 2] = 2; // O

  // OCH2OCHO
  ncf[32 * kd + 0] = 2; // C
  ncf[32 * kd + 1] = 3; // H
  ncf[32 * kd + 2] = 3; // O

  // HOCH2OCO
  ncf[33 * kd + 0] = 2; // C
  ncf[33 * kd + 1] = 3; // H
  ncf[33 * kd + 2] = 3; // O

  // CH3OCH2O2
  ncf[34 * kd + 0] = 2; // C
  ncf[34 * kd + 1] = 5; // H
  ncf[34 * kd + 2] = 3; // O

  // CH2OCH2O2H
  ncf[35 * kd + 0] = 2; // C
  ncf[35 * kd + 1] = 5; // H
  ncf[35 * kd + 2] = 3; // O

  // HO2CH2OCHO
  ncf[36 * kd + 0] = 2; // C
  ncf[36 * kd + 1] = 4; // H
  ncf[36 * kd + 2] = 4; // O

  // O2CH2OCH2O2H
  ncf[37 * kd + 0] = 2; // C
  ncf[37 * kd + 1] = 5; // H
  ncf[37 * kd + 2] = 5; // O

  // N2
  ncf[38 * kd + 3] = 2; // N
}

// Returns the vector of strings of element names
void
CKSYME_STR(amrex::Vector<std::string>& ename)
{
  ename.resize(4);
  ename[0] = "C";
  ename[1] = "H";
  ename[2] = "O";
  ename[3] = "N";
}

// Returns the vector of strings of species names
void
CKSYMS_STR(amrex::Vector<std::string>& kname)
{
  kname.resize(39);
  kname[0] = "H";
  kname[1] = "H2";
  kname[2] = "CH2";
  kname[3] = "CH2(S)";
  kname[4] = "CH3";
  kname[5] = "O";
  kname[6] = "CH4";
  kname[7] = "OH";
  kname[8] = "H2O";
  kname[9] = "C2H2";
  kname[10] = "C2H3";
  kname[11] = "CO";
  kname[12] = "C2H4";
  kname[13] = "HCO";
  kname[14] = "C2H5";
  kname[15] = "CH2O";
  kname[16] = "C2H6";
  kname[17] = "CH2OH";
  kname[18] = "CH3O";
  kname[19] = "O2";
  kname[20] = "HO2";
  kname[21] = "H2O2";
  kname[22] = "CO2";
  kname[23] = "CH3HCO";
  kname[24] = "CH3OCH2";
  kname[25] = "HCOOH";
  kname[26] = "CH3OCH3";
  kname[27] = "HOCH2O";
  kname[28] = "CH3OCO";
  kname[29] = "CH3OCHO";
  kname[30] = "CH3OCH2O";
  kname[31] = "CH3OCH2OH";
  kname[32] = "OCH2OCHO";
  kname[33] = "HOCH2OCO";
  kname[34] = "CH3OCH2O2";
  kname[35] = "CH2OCH2O2H";
  kname[36] = "HO2CH2OCHO";
  kname[37] = "O2CH2OCH2O2H";
  kname[38] = "N2";
}

// compute the sparsity pattern of the chemistry Jacobian
void
SPARSITY_INFO(int* nJdata, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 1600> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 39> conc = {0.0};
  for (int n = 0; n < 39; n++) {
    conc[n] = 1.0 / 39.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 40; k++) {
    for (int l = 0; l < 40; l++) {
      if (Jac[40 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 1600> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 39> conc = {0.0};
  for (int n = 0; n < 39; n++) {
    conc[n] = 1.0 / 39.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 40; k++) {
    for (int l = 0; l < 40; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[40 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 1600> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 39> conc = {0.0};
  for (int n = 0; n < 39; n++) {
    conc[n] = 1.0 / 39.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 40; k++) {
    for (int l = 0; l < 40; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[40 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 1600> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 39> conc = {0.0};
  for (int n = 0; n < 39; n++) {
    conc[n] = 1.0 / 39.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int nc = 0; nc < NCELLS; nc++) {
    int offset_row = nc * 40;
    int offset_col = nc * 40;
    for (int k = 0; k < 40; k++) {
      for (int l = 0; l < 40; l++) {
        if (Jac[40 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 1600> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 39> conc = {0.0};
  for (int n = 0; n < 39; n++) {
    conc[n] = 1.0 / 39.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtrs[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 40;
      for (int l = 0; l < 40; l++) {
        for (int k = 0; k < 40; k++) {
          if (Jac[40 * k + l] != 0.0) {
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
      int offset = nc * 40;
      for (int l = 0; l < 40; l++) {
        for (int k = 0; k < 40; k++) {
          if (Jac[40 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 1600> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 39> conc = {0.0};
  for (int n = 0; n < 39; n++) {
    conc[n] = 1.0 / 39.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 40;
      for (int l = 0; l < 40; l++) {
        for (int k = 0; k < 40; k++) {
          if (k == l) {
            colVals[nJdata_tmp - 1] = l + 1 + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[40 * k + l] != 0.0) {
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
      int offset = nc * 40;
      for (int l = 0; l < 40; l++) {
        for (int k = 0; k < 40; k++) {
          if (k == l) {
            colVals[nJdata_tmp] = l + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[40 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 1600> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 39> conc = {0.0};
  for (int n = 0; n < 39; n++) {
    conc[n] = 1.0 / 39.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int k = 0; k < 40; k++) {
    for (int l = 0; l < 40; l++) {
      if (k == l) {
        rowVals[nJdata_tmp] = l;
        indx[nJdata_tmp] = 40 * k + l;
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[40 * k + l] != 0.0) {
          rowVals[nJdata_tmp] = l;
          indx[nJdata_tmp] = 40 * k + l;
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
  amrex::GpuArray<amrex::Real, 1600> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 39> conc = {0.0};
  for (int n = 0; n < 39; n++) {
    conc[n] = 1.0 / 39.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int l = 0; l < 40; l++) {
      for (int k = 0; k < 40; k++) {
        if (k == l) {
          colVals[nJdata_tmp - 1] = l + 1;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[40 * k + l] != 0.0) {
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
    for (int l = 0; l < 40; l++) {
      for (int k = 0; k < 40; k++) {
        if (k == l) {
          colVals[nJdata_tmp] = l;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[40 * k + l] != 0.0) {
            colVals[nJdata_tmp] = k;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
      }
      rowPtr[l + 1] = nJdata_tmp;
    }
  }
}
