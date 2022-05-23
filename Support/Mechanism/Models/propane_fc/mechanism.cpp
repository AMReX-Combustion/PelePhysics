#include "mechanism.H"
const int rmap[167] = {
  6,   12,  45,  46,  60,  75,  92,  94,  103, 0,   3,   7,   66,  67,
  82,  83,  1,   2,   4,   5,   8,   9,   10,  11,  13,  14,  15,  16,
  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,  30,
  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,
  47,  48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,  61,
  62,  63,  64,  65,  68,  69,  70,  71,  72,  73,  74,  76,  77,  78,
  79,  80,  81,  84,  85,  86,  87,  88,  89,  90,  91,  93,  95,  96,
  97,  98,  99,  100, 101, 102, 104, 105, 106, 107, 108, 109, 110, 111,
  112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125,
  126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139,
  140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153,
  154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166};

// Returns 0-based map of reaction order
void
GET_RMAP(int* _rmap)
{
  for (int j = 0; j < 167; ++j) {
    _rmap[j] = rmap[j];
  }
}

// Returns a count of species in a reaction, and their indices
// and stoichiometric coefficients. (Eq 50)
void
CKINU(int* i, int* nspec, int* ki, int* nu)
{
  const int ns[167] = {
    3, 4, 4, 2, 4, 4, 2, 3, 3, 3, 4, 4, 3, 3, 3, 4, 4, 4, 3, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 5, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 3, 2, 3,
    3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4,
    4, 4, 4, 3, 4, 4, 3, 4, 3, 4, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 3, 4, 3, 4,
    4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4,
    5, 4, 4, 4, 4, 4, 4, 5, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3,
    3, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3};
  const int kiv[835] = {
    2,  1,  3,  0,  0, 4,  1,  2,  3,  0, 2,  3,  4,  1,  0, 1,  5,  0,  0,  0,
    4,  3,  2,  6,  0, 2,  6,  4,  3,  0, 3,  7,  0,  0,  0, 2,  3,  6,  0,  0,
    6,  1,  3,  0,  0, 3,  6,  1,  0,  0, 2,  5,  1,  3,  0, 1,  3,  2,  5,  0,
    2,  5,  8,  0,  0, 2,  8,  3,  0,  0, 8,  7,  5,  0,  0, 2,  8,  4,  5,  0,
    8,  3,  6,  5,  0, 6,  5,  8,  3,  0, 8,  7,  5,  0,  0, 8,  1,  5,  3,  0,
    5,  3,  8,  1,  0, 7,  3,  6,  8,  0, 2,  7,  6,  3,  0, 7,  3,  6,  8,  0,
    6,  8,  7,  3,  0, 7,  1,  8,  3,  0, 2,  7,  4,  8,  0, 4,  8,  2,  7,  0,
    9,  3,  10, 2,  0, 9,  4,  11, 2,  0, 11, 2,  9,  4,  0, 9,  5,  12, 2,  3,
    9,  1,  12, 2,  0, 11, 8,  13, 3,  0, 11, 5,  10, 3,  0, 11, 5,  13, 1,  0,
    13, 1,  11, 5,  0, 11, 14, 2,  0,  0, 11, 8,  15, 5,  0, 11, 1,  10, 2,  0,
    11, 3,  10, 4,  0, 11, 3,  9,  6,  0, 9,  6,  11, 3,  0, 11, 7,  15, 8,  0,
    9,  11, 16, 2,  0, 11, 2,  15, 0,  0, 11, 17, 0,  0,  0, 9,  15, 11, 0,  0,
    11, 9,  15, 0,  0, 15, 1,  11, 3,  0, 11, 3,  15, 1,  0, 15, 2,  11, 4,  0,
    11, 4,  15, 2,  0, 15, 3,  11, 6,  0, 11, 6,  15, 3,  0, 12, 8,  18, 3,  0,
    12, 5,  18, 1,  0, 18, 1,  12, 5,  0, 12, 3,  18, 2,  0, 18, 2,  12, 3,  0,
    12, 1,  18, 0,  0, 11, 19, 15, 12, 0, 2,  19, 12, 4,  0, 19, 5,  12, 8,  0,
    19, 1,  12, 3,  0, 19, 1,  18, 2,  0, 19, 12, 2,  0,  0, 12, 2,  19, 0,  0,
    19, 3,  12, 6,  0, 10, 2,  4,  19, 0, 10, 5,  19, 8,  0, 10, 3,  6,  19, 0,
    10, 8,  7,  19, 0, 10, 1,  19, 3,  0, 10, 11, 15, 19, 0, 13, 10, 2,  0,  0,
    13, 5,  10, 8,  0, 9,  18, 10, 12, 0, 11, 20, 13, 0,  0, 20, 8,  21, 5,  0,
    20, 13, 5,  0,  0, 10, 20, 21, 19, 0, 20, 11, 5,  0,  0, 11, 5,  20, 0,  0,
    21, 13, 3,  0,  0, 13, 3,  21, 0,  0, 22, 5,  23, 3,  0, 22, 1,  2,  23, 0,
    24, 5,  25, 1,  0, 24, 2,  22, 4,  0, 24, 5,  10, 19, 0, 24, 5,  22, 8,  0,
    24, 22, 2,  0,  0, 16, 1,  11, 19, 0, 16, 2,  14, 0,  0, 16, 3,  24, 6,  0,
    24, 6,  16, 3,  0, 16, 2,  24, 4,  0, 24, 4,  16, 2,  0, 16, 20, 24, 21, 0,
    16, 11, 24, 15, 0, 24, 15, 16, 11, 0, 16, 1,  25, 2,  0, 16, 22, 4,  0,  0,
    14, 11, 16, 15, 0, 14, 5,  16, 8,  0, 16, 8,  14, 5,  0, 14, 8,  17, 5,  0,
    14, 2,  17, 0,  0, 17, 2,  14, 4,  0, 14, 4,  17, 2,  0, 17, 3,  14, 6,  0,
    17, 9,  14, 11, 0, 17, 1,  14, 3,  0, 17, 11, 14, 15, 0, 23, 1,  12, 2,  0,
    23, 5,  18, 19, 0, 23, 3,  19, 0,  0, 2,  23, 9,  12, 0, 9,  12, 2,  23, 0,
    25, 5,  10, 12, 3, 26, 10, 27, 19, 0, 27, 19, 26, 10, 0, 26, 8,  28, 3,  0,
    26, 20, 28, 13, 0, 26, 8,  27, 5,  0, 26, 5,  25, 10, 0, 26, 5,  22, 10, 3,
    26, 2,  27, 0,  0, 26, 22, 11, 0,  0, 27, 24, 11, 0,  0, 24, 11, 27, 0,  0,
    27, 1,  14, 19, 0, 27, 2,  16, 11, 0, 16, 11, 27, 2,  0, 27, 2,  26, 4,  0,
    26, 4,  27, 2,  0, 27, 3,  26, 6,  0, 27, 1,  26, 3,  0, 27, 11, 26, 15, 0,
    29, 5,  27, 8,  0, 29, 27, 2,  0,  0, 27, 2,  29, 0,  0, 30, 16, 11, 0,  0,
    16, 11, 30, 0,  0, 8,  30, 31, 5,  0, 30, 5,  27, 8,  0, 30, 27, 2,  0,  0,
    27, 2,  30, 0,  0, 31, 3,  6,  30, 0, 31, 8,  7,  30, 0, 31, 2,  4,  30, 0,
    31, 3,  6,  29, 0, 31, 11, 15, 29, 0, 31, 11, 15, 30, 0, 31, 1,  29, 3,  0,
    31, 8,  7,  29, 0, 31, 1,  30, 3,  0, 31, 5,  8,  29, 0, 8,  29, 31, 5,  0,
    31, 2,  4,  29, 0, 4,  29, 31, 2,  0, 28, 24, 10, 0,  0, 32, 29, 5,  0,  0,
    29, 5,  32, 0,  0, 33, 30, 5,  0,  0, 30, 5,  33, 0,  0};
  const int nuv[835] = {
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -2, 1,  0, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -2, 1,  0, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 2, 0, 0, -2, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 2, 0, 0, -2, 1,  1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -2, 1,  1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 1,
    -1, -1, 1, 2, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -2, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -2, 1,  0, 0, 0, -1, -1, 2, 0, 0,
    -2, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, 1,  1, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, 1,  1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 2, 0, 0, -1, -1, 1, 1, 0,
    -2, 2,  1, 0, 0, -1, -1, 1, 1, 0, -1, 1,  1, 0, 0, -1, -1, 1, 0, 0,
    -1, 1,  1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, 1,  1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 2, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 2, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 1, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 1,
    -1, -1, 1, 0, 0, -1, 1,  1, 0, 0, -1, 1,  1, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, 1,  1, 0, 0, -1, -1, 1, 0, 0, -1, 1,  1, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, 1,  1, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, 1,  1, 0, 0, -1, 1,  1, 0, 0,
    -1, -1, 1, 0, 0, -1, 1,  1, 0, 0, -1, -1, 1, 0, 0};
  if (*i < 1) {
    // Return max num species per reaction
    *nspec = 5;
  } else {
    if (*i > 167) {
      *nspec = -1;
    } else {
      *nspec = ns[*i - 1];
      for (int j = 0; j < *nspec; ++j) {
        ki[j] = kiv[(*i - 1) * 5 + j] + 1;
        nu[j] = nuv[(*i - 1) * 5 + j];
      }
    }
  }
}

// Returns the progress rates of each reactions
// Given P, T, and mole fractions
void
CKKFKR(
  amrex::Real* P,
  amrex::Real* T,
  amrex::Real* x,
  amrex::Real* q_f,
  amrex::Real* q_r)
{
  int id;            // loop counter
  amrex::Real c[34]; // temporary storage
  amrex::Real PORT =
    1e6 * (*P) /
    (8.31446261815324e+07 * (*T)); // 1e6 * P/RT so c goes to SI units

  // Compute conversion, see Eq 10
  for (id = 0; id < 34; ++id) {
    c[id] = x[id] * PORT;
  }

  // convert to chemkin units
  progressRateFR(q_f, q_r, c, *T);

  // convert to chemkin units
  for (id = 0; id < 167; ++id) {
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
  amrex::Real g_RT[34];
  gibbs(g_RT, tc);

  amrex::Real sc_qss[1];
  amrex::Real kf_qss[0], qf_qss[0], qr_qss[0];
  comp_qfqr(q_f, q_r, sc, sc_qss, tc, invT);

  return;
}

// save atomic weights into array
void
atomicWeight(amrex::Real* awt)
{
  awt[0] = 14.007000; // N
  awt[1] = 15.999000; // O
  awt[2] = 1.008000;  // H
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
  int id; // loop counter
  int kd = 4;
  // Zero ncf
  for (id = 0; id < kd * 34; ++id) {
    ncf[id] = 0;
  }

  // N2
  ncf[0 * kd + 0] = 2; // N

  // O
  ncf[1 * kd + 1] = 1; // O

  // H
  ncf[2 * kd + 2] = 1; // H

  // OH
  ncf[3 * kd + 2] = 1; // H
  ncf[3 * kd + 1] = 1; // O

  // H2
  ncf[4 * kd + 2] = 2; // H

  // O2
  ncf[5 * kd + 1] = 2; // O

  // H2O
  ncf[6 * kd + 2] = 2; // H
  ncf[6 * kd + 1] = 1; // O

  // H2O2
  ncf[7 * kd + 2] = 2; // H
  ncf[7 * kd + 1] = 2; // O

  // HO2
  ncf[8 * kd + 2] = 1; // H
  ncf[8 * kd + 1] = 2; // O

  // CH2GSG
  ncf[9 * kd + 3] = 1; // C
  ncf[9 * kd + 2] = 2; // H

  // CH2O
  ncf[10 * kd + 3] = 1; // C
  ncf[10 * kd + 2] = 2; // H
  ncf[10 * kd + 1] = 1; // O

  // CH3
  ncf[11 * kd + 3] = 1; // C
  ncf[11 * kd + 2] = 3; // H

  // CO
  ncf[12 * kd + 3] = 1; // C
  ncf[12 * kd + 1] = 1; // O

  // CH3O
  ncf[13 * kd + 3] = 1; // C
  ncf[13 * kd + 2] = 3; // H
  ncf[13 * kd + 1] = 1; // O

  // C2H5
  ncf[14 * kd + 3] = 2; // C
  ncf[14 * kd + 2] = 5; // H

  // CH4
  ncf[15 * kd + 3] = 1; // C
  ncf[15 * kd + 2] = 4; // H

  // C2H4
  ncf[16 * kd + 3] = 2; // C
  ncf[16 * kd + 2] = 4; // H

  // C2H6
  ncf[17 * kd + 3] = 2; // C
  ncf[17 * kd + 2] = 6; // H

  // CO2
  ncf[18 * kd + 3] = 1; // C
  ncf[18 * kd + 1] = 2; // O

  // HCO
  ncf[19 * kd + 3] = 1; // C
  ncf[19 * kd + 2] = 1; // H
  ncf[19 * kd + 1] = 1; // O

  // CH3O2
  ncf[20 * kd + 3] = 1; // C
  ncf[20 * kd + 2] = 3; // H
  ncf[20 * kd + 1] = 2; // O

  // CH3O2H
  ncf[21 * kd + 3] = 1; // C
  ncf[21 * kd + 2] = 4; // H
  ncf[21 * kd + 1] = 2; // O

  // C2H2
  ncf[22 * kd + 3] = 2; // C
  ncf[22 * kd + 2] = 2; // H

  // HCCO
  ncf[23 * kd + 3] = 2; // C
  ncf[23 * kd + 2] = 1; // H
  ncf[23 * kd + 1] = 1; // O

  // C2H3
  ncf[24 * kd + 3] = 2; // C
  ncf[24 * kd + 2] = 3; // H

  // CH2CHO
  ncf[25 * kd + 3] = 2; // C
  ncf[25 * kd + 2] = 3; // H
  ncf[25 * kd + 1] = 1; // O

  // C3H5XA
  ncf[26 * kd + 3] = 3; // C
  ncf[26 * kd + 2] = 5; // H

  // C3H6
  ncf[27 * kd + 3] = 3; // C
  ncf[27 * kd + 2] = 6; // H

  // C3H5O
  ncf[28 * kd + 3] = 3; // C
  ncf[28 * kd + 2] = 5; // H
  ncf[28 * kd + 1] = 1; // O

  // IXC3H7
  ncf[29 * kd + 3] = 3; // C
  ncf[29 * kd + 2] = 7; // H

  // NXC3H7
  ncf[30 * kd + 3] = 3; // C
  ncf[30 * kd + 2] = 7; // H

  // C3H8
  ncf[31 * kd + 3] = 3; // C
  ncf[31 * kd + 2] = 8; // H

  // IXC3H7O2
  ncf[32 * kd + 3] = 3; // C
  ncf[32 * kd + 2] = 7; // H
  ncf[32 * kd + 1] = 2; // O

  // NXC3H7O2
  ncf[33 * kd + 3] = 3; // C
  ncf[33 * kd + 2] = 7; // H
  ncf[33 * kd + 1] = 2; // O
}

// Returns the vector of strings of element names
void
CKSYME_STR(amrex::Vector<std::string>& ename)
{
  ename.resize(4);
  ename[0] = "N";
  ename[1] = "O";
  ename[2] = "H";
  ename[3] = "C";
}

// Returns the vector of strings of species names
void
CKSYMS_STR(amrex::Vector<std::string>& kname)
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

// compute the sparsity pattern of the chemistry Jacobian
void
SPARSITY_INFO(int* nJdata, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 1225> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 34> conc = {0.0};
  for (int n = 0; n < 34; n++) {
    conc[n] = 1.0 / 34.000000;
  }
  aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 35; k++) {
    for (int l = 0; l < 35; l++) {
      if (Jac[35 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 1225> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 34> conc = {0.0};
  for (int n = 0; n < 34; n++) {
    conc[n] = 1.0 / 34.000000;
  }
  aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 35; k++) {
    for (int l = 0; l < 35; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[35 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 1225> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 34> conc = {0.0};
  for (int n = 0; n < 34; n++) {
    conc[n] = 1.0 / 34.000000;
  }
  aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 35; k++) {
    for (int l = 0; l < 35; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[35 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 1225> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 34> conc = {0.0};
  for (int n = 0; n < 34; n++) {
    conc[n] = 1.0 / 34.000000;
  }
  aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int nc = 0; nc < NCELLS; nc++) {
    int offset_row = nc * 35;
    int offset_col = nc * 35;
    for (int k = 0; k < 35; k++) {
      for (int l = 0; l < 35; l++) {
        if (Jac[35 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 1225> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 34> conc = {0.0};
  for (int n = 0; n < 34; n++) {
    conc[n] = 1.0 / 34.000000;
  }
  aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

  if (base == 1) {
    rowPtrs[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 35;
      for (int l = 0; l < 35; l++) {
        for (int k = 0; k < 35; k++) {
          if (Jac[35 * k + l] != 0.0) {
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
      int offset = nc * 35;
      for (int l = 0; l < 35; l++) {
        for (int k = 0; k < 35; k++) {
          if (Jac[35 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 1225> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 34> conc = {0.0};
  for (int n = 0; n < 34; n++) {
    conc[n] = 1.0 / 34.000000;
  }
  aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 35;
      for (int l = 0; l < 35; l++) {
        for (int k = 0; k < 35; k++) {
          if (k == l) {
            colVals[nJdata_tmp - 1] = l + 1 + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[35 * k + l] != 0.0) {
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
      int offset = nc * 35;
      for (int l = 0; l < 35; l++) {
        for (int k = 0; k < 35; k++) {
          if (k == l) {
            colVals[nJdata_tmp] = l + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[35 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 1225> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 34> conc = {0.0};
  for (int n = 0; n < 34; n++) {
    conc[n] = 1.0 / 34.000000;
  }
  aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int k = 0; k < 35; k++) {
    for (int l = 0; l < 35; l++) {
      if (k == l) {
        rowVals[nJdata_tmp] = l;
        indx[nJdata_tmp] = 35 * k + l;
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[35 * k + l] != 0.0) {
          rowVals[nJdata_tmp] = l;
          indx[nJdata_tmp] = 35 * k + l;
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
  amrex::GpuArray<amrex::Real, 1225> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 34> conc = {0.0};
  for (int n = 0; n < 34; n++) {
    conc[n] = 1.0 / 34.000000;
  }
  aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int l = 0; l < 35; l++) {
      for (int k = 0; k < 35; k++) {
        if (k == l) {
          colVals[nJdata_tmp - 1] = l + 1;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[35 * k + l] != 0.0) {
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
    for (int l = 0; l < 35; l++) {
      for (int k = 0; k < 35; k++) {
        if (k == l) {
          colVals[nJdata_tmp] = l;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[35 * k + l] != 0.0) {
            colVals[nJdata_tmp] = k;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
      }
      rowPtr[l + 1] = nJdata_tmp;
    }
  }
}
