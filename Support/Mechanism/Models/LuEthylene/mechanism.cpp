#include "mechanism.H"
const int rmap[206] = {
  15,  30,  38,  40,  47,  55,  70,  77,  88,  92,  114, 125, 131, 144, 147,
  154, 155, 169, 184, 189, 113, 4,   5,   6,   7,   8,   9,   10,  11,  12,
  13,  14,  28,  45,  120, 0,   1,   2,   3,   16,  17,  18,  19,  20,  21,
  22,  23,  24,  25,  26,  27,  29,  31,  32,  33,  34,  35,  36,  37,  39,
  41,  42,  43,  44,  46,  48,  49,  50,  51,  52,  53,  54,  56,  57,  58,
  59,  60,  61,  62,  63,  64,  65,  66,  67,  68,  69,  71,  72,  73,  74,
  75,  76,  78,  79,  80,  81,  82,  83,  84,  85,  86,  87,  89,  90,  91,
  93,  94,  95,  96,  97,  98,  99,  100, 101, 102, 103, 104, 105, 106, 107,
  108, 109, 110, 111, 112, 115, 116, 117, 118, 119, 121, 122, 123, 124, 126,
  127, 128, 129, 130, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142,
  143, 145, 146, 148, 149, 150, 151, 152, 153, 156, 157, 158, 159, 160, 161,
  162, 163, 164, 165, 166, 167, 168, 170, 171, 172, 173, 174, 175, 176, 177,
  178, 179, 180, 181, 182, 183, 185, 186, 187, 188, 190, 191, 192, 193, 194,
  195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205};

// Returns 0-based map of reaction order
void
GET_RMAP(int* _rmap)
{
  for (int j = 0; j < 206; ++j) {
    _rmap[j] = rmap[j];
  }
}

// Returns a count of species in a reaction, and their indices
// and stoichiometric coefficients. (Eq 50)
void
CKINU(const int i, int& nspec, int ki[], int nu[])
{
  const int ns[206] = {
    4, 4, 4, 3, 2, 2, 2, 2, 3, 3, 2, 3, 3, 3, 3, 2, 4, 4, 3, 4, 4, 3, 3,
    4, 4, 4, 4, 4, 3, 4, 3, 4, 4, 4, 4, 4, 4, 4, 3, 4, 3, 4, 4, 4, 4, 3,
    4, 3, 4, 4, 4, 4, 4, 4, 4, 3, 4, 3, 4, 4, 4, 4, 4, 4, 5, 4, 4, 4, 4,
    4, 3, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4,
    2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 3, 2, 3,
    4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 3, 4, 3, 3, 4, 4, 4, 4, 4, 5, 3, 3, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 5, 4, 4, 4, 4, 4, 4, 4,
    3, 4, 4, 5, 4, 3, 4, 4, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 4, 4};
  const int kiv[1030] = {
    1,  3,  2,  4,  0, 0,  2,  1,  4,  0, 0,  4,  1,  5,  0, 4,  5,  2,  0,  0,
    1,  0,  0,  0,  0, 1,  0,  0,  0,  0, 1,  0,  0,  0,  0, 1,  0,  0,  0,  0,
    1,  4,  5,  0,  0, 1,  2,  4,  0,  0, 2,  3,  0,  0,  0, 1,  3,  6,  0,  0,
    1,  3,  6,  0,  0, 1,  3,  6,  0,  0, 1,  3,  6,  0,  0, 4,  7,  0,  0,  0,
    1,  6,  5,  2,  0, 1,  6,  0,  3,  0, 1,  6,  4,  0,  0, 6,  2,  3,  4,  0,
    6,  4,  5,  3,  0, 6,  7,  3,  0,  0, 6,  7,  3,  0,  0, 1,  7,  0,  6,  0,
    1,  7,  5,  4,  0, 7,  2,  6,  4,  0, 7,  4,  5,  6,  0, 7,  4,  5,  6,  0,
    13, 2,  14, 0,  0, 13, 4,  14, 1,  0, 13, 0,  16, 0,  0, 13, 3,  14, 2,  0,
    13, 6,  14, 4,  0, 8,  2,  13, 1,  0, 8,  4,  1,  15, 0, 8,  0,  9,  1,  0,
    8,  5,  16, 1,  0, 8,  3,  15, 2,  0, 8,  13, 24, 0,  0, 8,  14, 13, 15, 0,
    1,  15, 16, 0,  0, 1,  15, 13, 0,  0, 15, 2,  13, 4,  0, 15, 2,  14, 1,  0,
    15, 4,  13, 5,  0, 15, 13, 1,  0,  0, 15, 3,  13, 6,  0, 9,  1,  11, 0,  0,
    9,  0,  11, 1,  0, 9,  2,  1,  15, 0, 9,  3,  15, 4,  0, 9,  3,  14, 1,  0,
    9,  4,  16, 1,  0, 9,  4,  8,  5,  0, 9,  6,  16, 4,  0, 9,  13, 25, 0,  0,
    8,  9,  18, 1,  0, 9,  18, 0,  0,  0, 10, 31, 9,  31, 0, 10, 1,  8,  0,  0,
    10, 2,  13, 0,  0, 10, 2,  1,  15, 0, 10, 4,  16, 1,  0, 10, 0,  11, 1,  0,
    10, 3,  13, 1,  4, 10, 3,  13, 5,  0, 10, 5,  9,  5,  0, 10, 13, 9,  13, 0,
    10, 14, 9,  14, 0, 10, 14, 16, 13, 0, 16, 1,  17, 0,  0, 16, 1,  0,  15, 0,
    16, 2,  15, 4,  0, 16, 4,  5,  15, 0, 16, 3,  15, 6,  0, 16, 6,  7,  15, 0,
    8,  16, 25, 1,  0, 11, 1,  12, 0,  0, 11, 2,  16, 1,  0, 11, 4,  9,  5,  0,
    11, 4,  10, 5,  0, 11, 3,  17, 2,  0, 11, 3,  16, 4,  0, 11, 6,  12, 3,  0,
    11, 6,  17, 4,  0, 11, 7,  12, 6,  0, 8,  11, 20, 1,  0, 11, 15, 12, 13, 0,
    11, 15, 27, 0,  0, 16, 11, 12, 15, 0, 9,  11, 21, 1,  0, 10, 11, 21, 1,  0,
    11, 23, 0,  0,  0, 11, 22, 1,  0,  0, 11, 24, 21, 13, 0, 17, 1,  16, 0,  0,
    17, 1,  11, 4,  0, 17, 1,  10, 5,  0, 17, 2,  16, 4,  0, 17, 4,  16, 5,  0,
    17, 3,  16, 6,  0, 12, 1,  11, 0,  0, 12, 2,  11, 4,  0, 12, 4,  11, 5,  0,
    8,  12, 21, 1,  0, 9,  12, 11, 0,  0, 10, 12, 11, 0,  0, 1,  24, 10, 13, 0,
    24, 2,  13, 1,  0, 24, 3,  13, 4,  0, 8,  24, 18, 13, 0, 9,  24, 20, 13, 0,
    24, 18, 13, 0,  0, 18, 19, 0,  0,  0, 20, 18, 1,  0,  0, 18, 2,  1,  24, 0,
    18, 2,  9,  13, 0, 18, 4,  25, 1,  0, 18, 4,  11, 13, 0, 18, 15, 20, 13, 0,
    18, 11, 28, 0,  0, 1,  19, 18, 1,  0, 19, 2,  9,  13, 0, 19, 4,  25, 1,  0,
    19, 3,  9,  14, 0, 25, 1,  26, 0,  0, 25, 1,  0,  24, 0, 25, 1,  11, 13, 0,
    25, 2,  24, 4,  0, 25, 2,  9,  14, 0, 25, 4,  5,  24, 0, 20, 1,  21, 0,  0,
    20, 1,  18, 0,  0, 20, 1,  0,  19, 0, 20, 2,  25, 1,  0, 20, 2,  11, 13, 0,
    20, 4,  18, 5,  0, 20, 3,  18, 6,  0, 20, 3,  26, 2,  0, 20, 3,  16, 15, 0,
    20, 6,  26, 4,  0, 20, 7,  21, 6,  0, 20, 15, 21, 13, 0, 20, 11, 18, 12, 0,
    20, 11, 29, 0,  0, 20, 11, 1,  28, 0, 26, 11, 13, 0,  0, 26, 1,  27, 0,  0,
    26, 1,  11, 15, 0, 26, 1,  25, 0,  0, 26, 2,  25, 4,  0, 26, 4,  25, 5,  0,
    26, 3,  25, 6,  0, 26, 3,  16, 13, 4, 21, 0,  19, 0,  0, 21, 1,  22, 0,  0,
    21, 1,  20, 0,  0, 21, 2,  20, 4,  0, 21, 2,  11, 15, 0, 21, 2,  9,  16, 0,
    21, 4,  20, 5,  0, 21, 3,  20, 6,  0, 21, 6,  27, 4,  0, 21, 15, 22, 13, 0,
    21, 9,  1,  28, 0, 21, 10, 12, 19, 0, 21, 10, 1,  28, 0, 21, 11, 20, 12, 0,
    21, 11, 30, 0,  0, 22, 1,  23, 0,  0, 22, 1,  21, 0,  0, 22, 2,  16, 11, 0,
    22, 2,  27, 1,  0, 22, 3,  21, 6,  0, 22, 6,  23, 3,  0, 22, 6,  21, 7,  0,
    22, 6,  16, 11, 4, 22, 7,  23, 6,  0, 22, 15, 23, 13, 0, 23, 1,  22, 0,  0,
    23, 2,  22, 4,  0, 23, 4,  22, 5,  0, 23, 10, 22, 11, 0, 23, 11, 22, 12, 0,
    1,  28, 29, 0,  0, 1,  28, 12, 19, 0, 6,  28, 29, 3,  0, 6,  28, 20, 16, 4,
    15, 28, 29, 13, 0, 29, 1,  30, 0,  0, 29, 1,  21, 11, 0, 29, 1,  0,  28, 0,
    29, 2,  25, 11, 1, 29, 2,  22, 15, 0, 29, 2,  4,  28, 0, 29, 4,  5,  28, 0,
    29, 6,  7,  28, 0, 29, 11, 12, 28, 0, 1,  30, 22, 11, 0, 1,  30, 29, 0,  0,
    2,  30, 22, 16, 0, 4,  30, 29, 5,  0, 3,  30, 29, 6,  0, 6,  30, 22, 16, 4,
    11, 30, 29, 12, 0, 20, 22, 11, 28, 0};
  const int nuv[1030] = {
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -2, 1,  1, 0, 0,
    -2, 1,  0, 0, 0, -2, 1,  0, 0, 0, -2, 1,  0, 0, 0, -2, 1,  0, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -2, 1,  0, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -2, 1,  0, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 2, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -2, 1,  1, 0, 0, -2, 1,  1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 2, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -2, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 1, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -2, 1,  0, 0, 0, -2, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 2, 0, 0, -1, -1, 2, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 2, 1, 0, -1, -1, 2, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -2, 1,  2, 0, 0, -1, 1,  0, 0, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, 1,  1, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 1, -1, 1,  1, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 1, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 1,
    -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 1, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 1,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0};
  if (i < 1) {
    // Return max num species per reaction
    nspec = 5;
  } else {
    if (i > 206) {
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
  int id;            // loop counter
  amrex::Real c[32]; // temporary storage
  amrex::Real PORT =
    1e6 * P / (8.31446261815324e+07 * T); // 1e6 * P/RT so c goes to SI units

  // Compute conversion, see Eq 10
  for (id = 0; id < 32; ++id) {
    c[id] = x[id] * PORT;
  }

  // convert to chemkin units
  progressRateFR(q_f, q_r, c, T);

  // convert to chemkin units
  for (id = 0; id < 206; ++id) {
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
  amrex::Real g_RT[32];
  gibbs(g_RT, tc);

  amrex::Real sc_qss[1];
  comp_qfqr(q_f, q_r, sc, sc_qss, tc, invT);
}

// save atomic weights into array
void
atomicWeight(amrex::Real* awt)
{
  awt[0] = 15.999000; // O
  awt[1] = 1.008000;  // H
  awt[2] = 12.011000; // C
  awt[3] = 14.007000; // N
  awt[4] = 39.950000; // Ar
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
  int kd = 5;
  // Zero ncf
  for (id = 0; id < kd * 32; ++id) {
    ncf[id] = 0;
  }

  // H2
  ncf[0 * kd + 1] = 2; // H

  // H
  ncf[1 * kd + 1] = 1; // H

  // O
  ncf[2 * kd + 0] = 1; // O

  // O2
  ncf[3 * kd + 0] = 2; // O

  // OH
  ncf[4 * kd + 1] = 1; // H
  ncf[4 * kd + 0] = 1; // O

  // H2O
  ncf[5 * kd + 1] = 2; // H
  ncf[5 * kd + 0] = 1; // O

  // HO2
  ncf[6 * kd + 1] = 1; // H
  ncf[6 * kd + 0] = 2; // O

  // H2O2
  ncf[7 * kd + 1] = 2; // H
  ncf[7 * kd + 0] = 2; // O

  // CH
  ncf[8 * kd + 2] = 1; // C
  ncf[8 * kd + 1] = 1; // H

  // CH2
  ncf[9 * kd + 2] = 1; // C
  ncf[9 * kd + 1] = 2; // H

  // CH2*
  ncf[10 * kd + 2] = 1; // C
  ncf[10 * kd + 1] = 2; // H

  // CH3
  ncf[11 * kd + 2] = 1; // C
  ncf[11 * kd + 1] = 3; // H

  // CH4
  ncf[12 * kd + 2] = 1; // C
  ncf[12 * kd + 1] = 4; // H

  // CO
  ncf[13 * kd + 2] = 1; // C
  ncf[13 * kd + 0] = 1; // O

  // CO2
  ncf[14 * kd + 2] = 1; // C
  ncf[14 * kd + 0] = 2; // O

  // HCO
  ncf[15 * kd + 2] = 1; // C
  ncf[15 * kd + 1] = 1; // H
  ncf[15 * kd + 0] = 1; // O

  // CH2O
  ncf[16 * kd + 2] = 1; // C
  ncf[16 * kd + 1] = 2; // H
  ncf[16 * kd + 0] = 1; // O

  // CH3O
  ncf[17 * kd + 2] = 1; // C
  ncf[17 * kd + 1] = 3; // H
  ncf[17 * kd + 0] = 1; // O

  // C2H2
  ncf[18 * kd + 2] = 2; // C
  ncf[18 * kd + 1] = 2; // H

  // H2CC
  ncf[19 * kd + 2] = 2; // C
  ncf[19 * kd + 1] = 2; // H

  // C2H3
  ncf[20 * kd + 2] = 2; // C
  ncf[20 * kd + 1] = 3; // H

  // C2H4
  ncf[21 * kd + 2] = 2; // C
  ncf[21 * kd + 1] = 4; // H

  // C2H5
  ncf[22 * kd + 2] = 2; // C
  ncf[22 * kd + 1] = 5; // H

  // C2H6
  ncf[23 * kd + 2] = 2; // C
  ncf[23 * kd + 1] = 6; // H

  // HCCO
  ncf[24 * kd + 2] = 2; // C
  ncf[24 * kd + 1] = 1; // H
  ncf[24 * kd + 0] = 1; // O

  // CH2CO
  ncf[25 * kd + 2] = 2; // C
  ncf[25 * kd + 1] = 2; // H
  ncf[25 * kd + 0] = 1; // O

  // CH2CHO
  ncf[26 * kd + 2] = 2; // C
  ncf[26 * kd + 1] = 3; // H
  ncf[26 * kd + 0] = 1; // O

  // CH3CHO
  ncf[27 * kd + 2] = 2; // C
  ncf[27 * kd + 1] = 4; // H
  ncf[27 * kd + 0] = 1; // O

  // aC3H5
  ncf[28 * kd + 2] = 3; // C
  ncf[28 * kd + 1] = 5; // H

  // C3H6
  ncf[29 * kd + 2] = 3; // C
  ncf[29 * kd + 1] = 6; // H

  // nC3H7
  ncf[30 * kd + 2] = 3; // C
  ncf[30 * kd + 1] = 7; // H

  // N2
  ncf[31 * kd + 3] = 2; // N
}

// Returns the vector of strings of element names
void
CKSYME_STR(amrex::Vector<std::string>& ename)
{
  ename.resize(5);
  ename[0] = "O";
  ename[1] = "H";
  ename[2] = "C";
  ename[3] = "N";
  ename[4] = "Ar";
}

// Returns the vector of strings of species names
void
CKSYMS_STR(amrex::Vector<std::string>& kname)
{
  kname.resize(32);
  kname[0] = "H2";
  kname[1] = "H";
  kname[2] = "O";
  kname[3] = "O2";
  kname[4] = "OH";
  kname[5] = "H2O";
  kname[6] = "HO2";
  kname[7] = "H2O2";
  kname[8] = "CH";
  kname[9] = "CH2";
  kname[10] = "CH2*";
  kname[11] = "CH3";
  kname[12] = "CH4";
  kname[13] = "CO";
  kname[14] = "CO2";
  kname[15] = "HCO";
  kname[16] = "CH2O";
  kname[17] = "CH3O";
  kname[18] = "C2H2";
  kname[19] = "H2CC";
  kname[20] = "C2H3";
  kname[21] = "C2H4";
  kname[22] = "C2H5";
  kname[23] = "C2H6";
  kname[24] = "HCCO";
  kname[25] = "CH2CO";
  kname[26] = "CH2CHO";
  kname[27] = "CH3CHO";
  kname[28] = "aC3H5";
  kname[29] = "C3H6";
  kname[30] = "nC3H7";
  kname[31] = "N2";
}

// compute the sparsity pattern of the chemistry Jacobian
void
SPARSITY_INFO(int* nJdata, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 1089> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 32> conc = {0.0};
  for (int n = 0; n < 32; n++) {
    conc[n] = 1.0 / 32.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 33; k++) {
    for (int l = 0; l < 33; l++) {
      if (Jac[33 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 1089> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 32> conc = {0.0};
  for (int n = 0; n < 32; n++) {
    conc[n] = 1.0 / 32.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 33; k++) {
    for (int l = 0; l < 33; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[33 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 1089> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 32> conc = {0.0};
  for (int n = 0; n < 32; n++) {
    conc[n] = 1.0 / 32.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 33; k++) {
    for (int l = 0; l < 33; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[33 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 1089> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 32> conc = {0.0};
  for (int n = 0; n < 32; n++) {
    conc[n] = 1.0 / 32.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int nc = 0; nc < NCELLS; nc++) {
    int offset_row = nc * 33;
    int offset_col = nc * 33;
    for (int k = 0; k < 33; k++) {
      for (int l = 0; l < 33; l++) {
        if (Jac[33 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 1089> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 32> conc = {0.0};
  for (int n = 0; n < 32; n++) {
    conc[n] = 1.0 / 32.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtrs[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 33;
      for (int l = 0; l < 33; l++) {
        for (int k = 0; k < 33; k++) {
          if (Jac[33 * k + l] != 0.0) {
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
      int offset = nc * 33;
      for (int l = 0; l < 33; l++) {
        for (int k = 0; k < 33; k++) {
          if (Jac[33 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 1089> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 32> conc = {0.0};
  for (int n = 0; n < 32; n++) {
    conc[n] = 1.0 / 32.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 33;
      for (int l = 0; l < 33; l++) {
        for (int k = 0; k < 33; k++) {
          if (k == l) {
            colVals[nJdata_tmp - 1] = l + 1 + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[33 * k + l] != 0.0) {
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
      int offset = nc * 33;
      for (int l = 0; l < 33; l++) {
        for (int k = 0; k < 33; k++) {
          if (k == l) {
            colVals[nJdata_tmp] = l + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[33 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 1089> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 32> conc = {0.0};
  for (int n = 0; n < 32; n++) {
    conc[n] = 1.0 / 32.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int k = 0; k < 33; k++) {
    for (int l = 0; l < 33; l++) {
      if (k == l) {
        rowVals[nJdata_tmp] = l;
        indx[nJdata_tmp] = 33 * k + l;
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[33 * k + l] != 0.0) {
          rowVals[nJdata_tmp] = l;
          indx[nJdata_tmp] = 33 * k + l;
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
  amrex::GpuArray<amrex::Real, 1089> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 32> conc = {0.0};
  for (int n = 0; n < 32; n++) {
    conc[n] = 1.0 / 32.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int l = 0; l < 33; l++) {
      for (int k = 0; k < 33; k++) {
        if (k == l) {
          colVals[nJdata_tmp - 1] = l + 1;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[33 * k + l] != 0.0) {
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
    for (int l = 0; l < 33; l++) {
      for (int k = 0; k < 33; k++) {
        if (k == l) {
          colVals[nJdata_tmp] = l;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[33 * k + l] != 0.0) {
            colVals[nJdata_tmp] = k;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
      }
      rowPtr[l + 1] = nJdata_tmp;
    }
  }
}
