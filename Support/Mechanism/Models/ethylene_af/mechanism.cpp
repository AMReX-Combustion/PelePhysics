#include "mechanism.H"
const int rmap[172] = {
  5,   10,  30,  32,  56,  70,  74,  75,  77,  80,  102, 111, 115, 125, 135,
  0,   1,   3,   4,   6,   9,   82,  83,  96,  166, 2,   7,   8,   11,  12,
  13,  14,  15,  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,
  28,  29,  31,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,
  45,  46,  47,  48,  49,  50,  51,  52,  53,  54,  55,  57,  58,  59,  60,
  61,  62,  63,  64,  65,  66,  67,  68,  69,  71,  72,  73,  76,  78,  79,
  81,  84,  85,  86,  87,  88,  89,  90,  91,  92,  93,  94,  95,  97,  98,
  99,  100, 101, 103, 104, 105, 106, 107, 108, 109, 110, 112, 113, 114, 116,
  117, 118, 119, 120, 121, 122, 123, 124, 126, 127, 128, 129, 130, 131, 132,
  133, 134, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148,
  149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163,
  164, 165, 167, 168, 169, 170, 171};

// Returns 0-based map of reaction order
void
GET_RMAP(int* _rmap)
{
  for (int j = 0; j < 172; ++j) {
    _rmap[j] = rmap[j];
  }
}

// Returns a count of species in a reaction, and their indices
// and stoichiometric coefficients. (Eq 50)
void
CKINU(const int i, int& nspec, int ki[], int nu[])
{
  const int ns[172] = {
    2, 2, 4, 3, 2, 2, 3, 3, 4, 2, 3, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 3, 4, 3, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 3, 4, 4, 5, 4, 3, 4, 4, 3, 4, 3, 4, 4, 4, 3,
    3, 4, 3, 4, 4, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 2, 4, 4, 4,
    4, 4, 3, 4, 4, 4, 4, 3, 4, 4, 4, 3, 4, 4, 4, 3, 5, 4, 3, 4, 4, 4, 4, 4, 4,
    3, 4, 3, 4, 4, 4, 4, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 3, 4, 4, 5, 3, 4, 4, 5, 5, 5, 4, 5, 4, 3, 5, 3, 3, 4, 4};
  const int kiv[860] = {
    1,  2,  0,  0,  0, 1,  2,  0,  0,  0, 2,  3,  1,  4,  0, 1,  3,  4,  0,  0,
    3,  5,  0,  0,  0, 4,  6,  0,  0,  0, 1,  4,  7,  0,  0, 4,  7,  3,  0,  0,
    2,  4,  1,  7,  0, 1,  2,  0,  0,  0, 1,  5,  8,  0,  0, 1,  5,  3,  4,  0,
    2,  5,  1,  8,  0, 8,  4,  7,  5,  0, 1,  8,  4,  0,  0, 8,  3,  5,  4,  0,
    1,  8,  7,  3,  0, 8,  4,  7,  5,  0, 6,  3,  8,  4,  0, 1,  6,  7,  4,  0,
    1,  6,  2,  8,  0, 6,  4,  7,  8,  0, 6,  4,  7,  8,  0, 18, 5,  9,  3,  0,
    18, 4,  9,  1,  0, 19, 4,  1,  20, 0, 19, 2,  1,  21, 0, 19, 3,  9,  1,  0,
    19, 5,  20, 3,  0, 19, 1,  18, 2,  0, 19, 2,  10, 0,  0, 19, 7,  11, 1,  0,
    1,  21, 10, 0,  0, 5,  21, 9,  1,  4, 5,  21, 11, 3,  0, 4,  21, 11, 1,  0,
    8,  21, 11, 4,  0, 5,  21, 12, 1,  0, 4,  21, 19, 7,  0, 3,  21, 1,  20, 0,
    2,  21, 10, 1,  0, 7,  22, 7,  21, 0, 1,  22, 19, 2,  0, 5,  22, 9,  1,  4,
    3,  22, 9,  2,  0, 5,  22, 9,  7,  0, 2,  22, 10, 1,  0, 3,  22, 1,  20, 0,
    7,  22, 11, 2,  0, 4,  22, 11, 1,  0, 10, 4,  11, 2,  0, 10, 6,  13, 8,  0,
    10, 5,  11, 4,  0, 19, 10, 23, 1,  0, 10, 3,  11, 1,  0, 18, 10, 14, 1,  0,
    10, 1,  13, 0,  0, 10, 4,  7,  21, 0, 10, 22, 15, 1,  0, 10, 4,  7,  22, 0,
    10, 24, 1,  0,  0, 10, 8,  13, 5,  0, 10, 21, 15, 1,  0, 10, 3,  9,  1,  2,
    19, 13, 15, 1,  0, 13, 22, 10, 0,  0, 13, 3,  10, 4,  0, 13, 4,  10, 7,  0,
    13, 21, 10, 0,  0, 13, 1,  10, 2,  0, 9,  21, 16, 0,  0, 9,  22, 9,  21, 0,
    9,  5,  12, 3,  0, 9,  4,  12, 1,  0, 9,  2,  11, 0,  0, 19, 9,  25, 0,  0,
    9,  4,  12, 1,  0, 9,  3,  12, 0,  0, 9,  8,  12, 4,  0, 1,  20, 9,  2,  0,
    1,  20, 11, 0,  0, 10, 20, 26, 0,  0, 20, 9,  1,  0,  0, 20, 9,  1,  0,  0,
    20, 3,  9,  4,  0, 20, 4,  9,  7,  0, 10, 20, 13, 9,  0, 20, 3,  12, 1,  0,
    20, 5,  9,  8,  0, 11, 1,  2,  20, 0, 11, 3,  20, 4,  0, 11, 10, 13, 20, 0,
    11, 4,  7,  20, 0, 19, 11, 16, 1,  0, 11, 5,  20, 8,  0, 11, 8,  6,  20, 0,
    1,  2,  0,  0,  0, 12, 22, 12, 21, 0, 12, 22, 11, 9,  0, 19, 12, 9,  20, 0,
    14, 3,  9,  21, 0, 14, 4,  10, 9,  0, 14, 1,  23, 0,  0, 14, 4,  16, 1,  0,
    14, 3,  1,  25, 0, 23, 4,  14, 7,  0, 23, 5,  27, 3,  0, 23, 3,  27, 0,  0,
    23, 1,  14, 2,  0, 23, 10, 14, 13, 0, 23, 5,  11, 20, 0, 23, 1,  15, 0,  0,
    23, 6,  15, 8,  0, 23, 5,  14, 8,  0, 15, 10, 23, 13, 0, 15, 1,  24, 0,  0,
    15, 5,  10, 12, 1, 15, 4,  23, 7,  0, 15, 4,  28, 0,  0, 15, 3,  27, 1,  0,
    15, 3,  10, 20, 0, 15, 5,  23, 8,  0, 15, 1,  23, 2,  0, 15, 3,  11, 21, 0,
    24, 8,  15, 6,  0, 24, 1,  17, 0,  0, 24, 8,  28, 4,  0, 24, 3,  28, 0,  0,
    24, 1,  15, 2,  0, 24, 5,  15, 8,  0, 24, 8,  17, 5,  0, 24, 10, 15, 13, 0,
    17, 22, 24, 10, 0, 17, 10, 24, 13, 0, 17, 3,  24, 4,  0, 17, 10, 0,  0,  0,
    17, 8,  24, 6,  0, 17, 1,  24, 2,  0, 17, 4,  24, 7,  0, 25, 5,  9,  4,  0,
    25, 3,  9,  1,  0, 10, 25, 15, 9,  0, 1,  25, 9,  22, 0, 16, 1,  10, 9,  0,
    16, 21, 15, 9,  0, 16, 3,  25, 4,  0, 16, 10, 13, 25, 0, 16, 3,  12, 21, 0,
    16, 10, 24, 9,  0, 16, 4,  7,  25, 0, 16, 1,  2,  25, 0, 16, 21, 10, 25, 0,
    27, 3,  11, 20, 0, 27, 16, 1,  0,  0, 27, 4,  16, 7,  0, 27, 1,  16, 2,  0,
    27, 5,  11, 9,  4, 27, 10, 9,  0,  0, 27, 5,  20, 4,  0, 27, 1,  10, 20, 0,
    26, 3,  10, 9,  4, 26, 5,  10, 9,  8, 26, 4,  10, 9,  7, 26, 1,  27, 2,  0,
    26, 1,  10, 9,  2, 26, 3,  27, 4,  0, 26, 13, 9,  0,  0, 26, 8,  10, 9,  6,
    28, 11, 10, 0,  0, 28, 26, 1,  0,  0, 28, 5,  26, 8,  0, 0,  22, 0,  21, 0};
  const int nuv[860] = {
    -2, 1,  0, 0, 0, -2, 1,  0, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0,
    -2, 1,  0, 0, 0, -2, 1,  0, 0, 0, -1, -1, 1, 0, 0, -2, 1,  1, 0, 0,
    -1, -1, 1, 1, 0, -2, 1,  0, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 2, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 1, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 2, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 1,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -2, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 1,
    -1, -1, 1, 1, 0, -1, -1, 2, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 2, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, 1,  1, 0, 0, -1, 1,  1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -2, 1,  0, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 1, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, 2,  0, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 2, 1, 0,
    -1, -1, 2, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 1, -1, 1,  1, 0, 0, -1, -1, 2, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 1, -1, -1, 1, 1, 1, -1, -1, 1, 1, 1, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 1, -1, -1, 1, 1, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 1,
    -1, 1,  1, 0, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0};
  if (i < 1) {
    // Return max num species per reaction
    nspec = 5;
  } else {
    if (i > 172) {
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
  amrex::Real c[29]; // temporary storage
  amrex::Real PORT =
    1e6 * P / (8.31446261815324e+07 * T); // 1e6 * P/RT so c goes to SI units

  // Compute conversion, see Eq 10
  for (int id = 0; id < 29; ++id) {
    c[id] = x[id] * PORT;
  }

  // convert to chemkin units
  progressRateFR(q_f, q_r, c, T);

  // convert to chemkin units
  for (int id = 0; id < 172; ++id) {
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
  amrex::Real invT = 1.0 / T;
  amrex::Real logT = log(T);
  // compute the Gibbs free energy
  amrex::Real g_RT[29];
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
  for (int id = 0; id < kd * 29; ++id) {
    ncf[id] = 0;
  }

  // N2
  ncf[0 * kd + 0] = 2; // N

  // H
  ncf[1 * kd + 1] = 1; // H

  // H2
  ncf[2 * kd + 1] = 2; // H

  // O
  ncf[3 * kd + 2] = 1; // O

  // OH
  ncf[4 * kd + 1] = 1; // H
  ncf[4 * kd + 2] = 1; // O

  // O2
  ncf[5 * kd + 2] = 2; // O

  // H2O2
  ncf[6 * kd + 1] = 2; // H
  ncf[6 * kd + 2] = 2; // O

  // H2O
  ncf[7 * kd + 1] = 2; // H
  ncf[7 * kd + 2] = 1; // O

  // HO2
  ncf[8 * kd + 1] = 1; // H
  ncf[8 * kd + 2] = 2; // O

  // CO
  ncf[9 * kd + 3] = 1; // C
  ncf[9 * kd + 2] = 1; // O

  // CH3
  ncf[10 * kd + 3] = 1; // C
  ncf[10 * kd + 1] = 3; // H

  // CH2O
  ncf[11 * kd + 3] = 1; // C
  ncf[11 * kd + 1] = 2; // H
  ncf[11 * kd + 2] = 1; // O

  // CO2
  ncf[12 * kd + 3] = 1; // C
  ncf[12 * kd + 2] = 2; // O

  // CH4
  ncf[13 * kd + 3] = 1; // C
  ncf[13 * kd + 1] = 4; // H

  // C2H2
  ncf[14 * kd + 3] = 2; // C
  ncf[14 * kd + 1] = 2; // H

  // C2H4
  ncf[15 * kd + 3] = 2; // C
  ncf[15 * kd + 1] = 4; // H

  // CH2CO
  ncf[16 * kd + 3] = 2; // C
  ncf[16 * kd + 1] = 2; // H
  ncf[16 * kd + 2] = 1; // O

  // C2H6
  ncf[17 * kd + 3] = 2; // C
  ncf[17 * kd + 1] = 6; // H

  // C
  ncf[18 * kd + 3] = 1; // C

  // CH
  ncf[19 * kd + 3] = 1; // C
  ncf[19 * kd + 1] = 1; // H

  // HCO
  ncf[20 * kd + 3] = 1; // C
  ncf[20 * kd + 1] = 1; // H
  ncf[20 * kd + 2] = 1; // O

  // TXCH2
  ncf[21 * kd + 3] = 1; // C
  ncf[21 * kd + 1] = 2; // H

  // SXCH2
  ncf[22 * kd + 3] = 1; // C
  ncf[22 * kd + 1] = 2; // H

  // C2H3
  ncf[23 * kd + 3] = 2; // C
  ncf[23 * kd + 1] = 3; // H

  // C2H5
  ncf[24 * kd + 3] = 2; // C
  ncf[24 * kd + 1] = 5; // H

  // HCCO
  ncf[25 * kd + 3] = 2; // C
  ncf[25 * kd + 1] = 1; // H
  ncf[25 * kd + 2] = 1; // O

  // CH3CHO
  ncf[26 * kd + 3] = 2; // C
  ncf[26 * kd + 1] = 4; // H
  ncf[26 * kd + 2] = 1; // O

  // CH2CHO
  ncf[27 * kd + 3] = 2; // C
  ncf[27 * kd + 1] = 3; // H
  ncf[27 * kd + 2] = 1; // O

  // C2H5O
  ncf[28 * kd + 3] = 2; // C
  ncf[28 * kd + 1] = 5; // H
  ncf[28 * kd + 2] = 1; // O
}

// Returns the vector of strings of element names
void
CKSYME_STR(amrex::Vector<std::string>& ename)
{
  ename.resize(4);
  ename[0] = "N";
  ename[1] = "H";
  ename[2] = "O";
  ename[3] = "C";
}

// Returns the vector of strings of species names
void
CKSYMS_STR(amrex::Vector<std::string>& kname)
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

// compute the sparsity pattern of the chemistry Jacobian
void
SPARSITY_INFO(int* nJdata, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 900> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 29> conc = {0.0};
  for (int n = 0; n < 29; n++) {
    conc[n] = 1.0 / 29.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 30; k++) {
    for (int l = 0; l < 30; l++) {
      if (Jac[30 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 900> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 29> conc = {0.0};
  for (int n = 0; n < 29; n++) {
    conc[n] = 1.0 / 29.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 30; k++) {
    for (int l = 0; l < 30; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[30 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 900> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 29> conc = {0.0};
  for (int n = 0; n < 29; n++) {
    conc[n] = 1.0 / 29.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 30; k++) {
    for (int l = 0; l < 30; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[30 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 900> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 29> conc = {0.0};
  for (int n = 0; n < 29; n++) {
    conc[n] = 1.0 / 29.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int nc = 0; nc < NCELLS; nc++) {
    int offset_row = nc * 30;
    int offset_col = nc * 30;
    for (int k = 0; k < 30; k++) {
      for (int l = 0; l < 30; l++) {
        if (Jac[30 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 900> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 29> conc = {0.0};
  for (int n = 0; n < 29; n++) {
    conc[n] = 1.0 / 29.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtrs[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 30;
      for (int l = 0; l < 30; l++) {
        for (int k = 0; k < 30; k++) {
          if (Jac[30 * k + l] != 0.0) {
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
      int offset = nc * 30;
      for (int l = 0; l < 30; l++) {
        for (int k = 0; k < 30; k++) {
          if (Jac[30 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 900> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 29> conc = {0.0};
  for (int n = 0; n < 29; n++) {
    conc[n] = 1.0 / 29.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 30;
      for (int l = 0; l < 30; l++) {
        for (int k = 0; k < 30; k++) {
          if (k == l) {
            colVals[nJdata_tmp - 1] = l + 1 + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[30 * k + l] != 0.0) {
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
      int offset = nc * 30;
      for (int l = 0; l < 30; l++) {
        for (int k = 0; k < 30; k++) {
          if (k == l) {
            colVals[nJdata_tmp] = l + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[30 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 900> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 29> conc = {0.0};
  for (int n = 0; n < 29; n++) {
    conc[n] = 1.0 / 29.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int k = 0; k < 30; k++) {
    for (int l = 0; l < 30; l++) {
      if (k == l) {
        rowVals[nJdata_tmp] = l;
        indx[nJdata_tmp] = 30 * k + l;
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[30 * k + l] != 0.0) {
          rowVals[nJdata_tmp] = l;
          indx[nJdata_tmp] = 30 * k + l;
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
  amrex::GpuArray<amrex::Real, 900> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 29> conc = {0.0};
  for (int n = 0; n < 29; n++) {
    conc[n] = 1.0 / 29.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int l = 0; l < 30; l++) {
      for (int k = 0; k < 30; k++) {
        if (k == l) {
          colVals[nJdata_tmp - 1] = l + 1;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[30 * k + l] != 0.0) {
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
    for (int l = 0; l < 30; l++) {
      for (int k = 0; k < 30; k++) {
        if (k == l) {
          colVals[nJdata_tmp] = l;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[30 * k + l] != 0.0) {
            colVals[nJdata_tmp] = k;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
      }
      rowPtr[l + 1] = nJdata_tmp;
    }
  }
}
