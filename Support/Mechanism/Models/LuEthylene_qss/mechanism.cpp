#include "mechanism.H"
const int rmap[206] = {
  0,   1,   2,   3,   4,   5,   6,   7,   8,   9,   10,  11,  12,  13,  14,
  15,  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  29,
  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,  42,  43,  44,
  45,  46,  47,  48,  49,  50,  51,  52,  53,  54,  55,  56,  57,  58,  59,
  60,  61,  62,  63,  64,  65,  66,  67,  68,  69,  70,  71,  72,  73,  74,
  75,  76,  77,  78,  79,  80,  81,  82,  83,  84,  85,  86,  87,  88,  89,
  90,  91,  92,  93,  94,  95,  96,  97,  98,  99,  100, 101, 102, 103, 104,
  105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119,
  120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134,
  135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149,
  150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164,
  165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179,
  180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194,
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
    2, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 2,
    2, 2, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 4, 4, 4, 3, 4, 4, 3, 4, 4, 3, 3,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 3, 3, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 5, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 5, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 5, 4, 4, 4, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 4, 4};
  const int kiv[1030] = {
    4,  7,  0,  0,  0, 10, 0,  12, 0,  0, 22, 10, 16, 0,  0, 1,  25, 12, 0,  0,
    23, 1,  8,  0,  0, 23, 10, 17, 0,  0, 12, 1,  26, 0,  0, 8,  1,  9,  0,  0,
    8,  25, 18, 0,  0, 8,  15, 0,  0,  0, 28, 13, 1,  0,  0, 17, 1,  30, 0,  0,
    28, 1,  14, 0,  0, 28, 8,  20, 0,  0, 30, 1,  18, 0,  0, 14, 0,  27, 0,  0,
    14, 1,  29, 0,  0, 29, 1,  15, 0,  0, 1,  19, 20, 0,  0, 20, 1,  31, 0,  0,
    13, 27, 0,  0,  0, 1,  0,  0,  0,  0, 1,  0,  0,  0,  0, 1,  0,  0,  0,  0,
    1,  0,  0,  0,  0, 1,  4,  5,  0,  0, 1,  2,  4,  0,  0, 2,  3,  0,  0,  0,
    1,  3,  6,  0,  0, 1,  3,  6,  0,  0, 1,  3,  6,  0,  0, 1,  3,  6,  0,  0,
    10, 2,  11, 0,  0, 25, 10, 1,  0,  0, 13, 8,  19, 0,  0, 1,  3,  2,  4,  0,
    0,  2,  1,  4,  0, 0,  4,  1,  5,  0, 4,  5,  2,  0,  0, 1,  6,  5,  2,  0,
    1,  6,  0,  3,  0, 1,  6,  4,  0,  0, 6,  2,  3,  4,  0, 6,  4,  5,  3,  0,
    6,  7,  3,  0,  0, 6,  7,  3,  0,  0, 1,  7,  0,  6,  0, 1,  7,  5,  4,  0,
    7,  2,  6,  4,  0, 7,  4,  5,  6,  0, 7,  4,  5,  6,  0, 10, 4,  11, 1,  0,
    10, 3,  11, 2,  0, 10, 6,  11, 4,  0, 22, 2,  10, 1,  0, 22, 4,  1,  25, 0,
    22, 0,  23, 1,  0, 22, 5,  12, 1,  0, 22, 3,  25, 2,  0, 22, 11, 10, 25, 0,
    1,  25, 10, 0,  0, 25, 2,  10, 4,  0, 25, 2,  11, 1,  0, 25, 4,  10, 5,  0,
    25, 3,  10, 6,  0, 23, 0,  8,  1,  0, 23, 2,  1,  25, 0, 23, 3,  25, 4,  0,
    23, 3,  11, 1,  0, 23, 4,  12, 1,  0, 23, 4,  22, 5,  0, 23, 6,  12, 4,  0,
    22, 23, 13, 1,  0, 23, 13, 0,  0,  0, 24, 21, 23, 21, 0, 24, 1,  22, 0,  0,
    24, 2,  10, 0,  0, 24, 2,  1,  25, 0, 24, 4,  12, 1,  0, 24, 0,  8,  1,  0,
    24, 3,  10, 1,  4, 24, 3,  10, 5,  0, 24, 5,  23, 5,  0, 24, 10, 23, 10, 0,
    24, 11, 23, 11, 0, 24, 11, 12, 10, 0, 12, 1,  0,  25, 0, 12, 2,  25, 4,  0,
    12, 4,  5,  25, 0, 12, 3,  25, 6,  0, 12, 6,  7,  25, 0, 22, 12, 17, 1,  0,
    8,  2,  12, 1,  0, 8,  4,  23, 5,  0, 8,  4,  24, 5,  0, 8,  3,  26, 2,  0,
    8,  3,  12, 4,  0, 8,  6,  9,  3,  0, 8,  6,  26, 4,  0, 8,  7,  9,  6,  0,
    22, 8,  28, 1,  0, 8,  25, 9,  10, 0, 12, 8,  9,  25, 0, 23, 8,  14, 1,  0,
    24, 8,  14, 1,  0, 8,  29, 1,  0,  0, 8,  16, 14, 10, 0, 26, 1,  12, 0,  0,
    26, 1,  8,  4,  0, 26, 1,  24, 5,  0, 26, 2,  12, 4,  0, 26, 4,  12, 5,  0,
    26, 3,  12, 6,  0, 9,  1,  8,  0,  0, 9,  2,  8,  4,  0, 9,  4,  8,  5,  0,
    22, 9,  14, 1,  0, 23, 9,  8,  0,  0, 24, 9,  8,  0,  0, 1,  16, 24, 10, 0,
    16, 2,  10, 1,  0, 16, 3,  10, 4,  0, 22, 16, 13, 10, 0, 23, 16, 28, 10, 0,
    16, 13, 10, 0,  0, 13, 2,  1,  16, 0, 13, 2,  23, 10, 0, 13, 4,  17, 1,  0,
    13, 4,  8,  10, 0, 13, 25, 28, 10, 0, 1,  27, 13, 1,  0, 27, 2,  23, 10, 0,
    27, 4,  17, 1,  0, 27, 3,  23, 11, 0, 17, 1,  0,  16, 0, 17, 1,  8,  10, 0,
    17, 2,  16, 4,  0, 17, 2,  23, 11, 0, 17, 4,  5,  16, 0, 28, 1,  13, 0,  0,
    28, 1,  0,  27, 0, 28, 2,  17, 1,  0, 28, 2,  8,  10, 0, 28, 4,  13, 5,  0,
    28, 3,  13, 6,  0, 28, 3,  30, 2,  0, 28, 3,  12, 25, 0, 28, 6,  30, 4,  0,
    28, 7,  14, 6,  0, 28, 25, 14, 10, 0, 28, 8,  13, 9,  0, 28, 8,  1,  19, 0,
    30, 8,  10, 0,  0, 30, 1,  8,  25, 0, 30, 1,  17, 0,  0, 30, 2,  17, 4,  0,
    30, 4,  17, 5,  0, 30, 3,  17, 6,  0, 30, 3,  12, 10, 4, 14, 1,  28, 0,  0,
    14, 2,  28, 4,  0, 14, 2,  8,  25, 0, 14, 2,  23, 12, 0, 14, 4,  28, 5,  0,
    14, 3,  28, 6,  0, 14, 6,  18, 4,  0, 14, 25, 29, 10, 0, 14, 23, 1,  19, 0,
    14, 24, 9,  27, 0, 14, 24, 1,  19, 0, 14, 8,  28, 9,  0, 14, 8,  31, 0,  0,
    29, 1,  14, 0,  0, 29, 2,  12, 8,  0, 29, 2,  18, 1,  0, 29, 3,  14, 6,  0,
    29, 6,  15, 3,  0, 29, 6,  14, 7,  0, 29, 6,  12, 8,  4, 29, 7,  15, 6,  0,
    29, 25, 15, 10, 0, 15, 1,  29, 0,  0, 15, 2,  29, 4,  0, 15, 4,  29, 5,  0,
    15, 24, 29, 8,  0, 15, 8,  29, 9,  0, 1,  19, 9,  27, 0, 6,  19, 20, 3,  0,
    6,  19, 28, 12, 4, 25, 19, 20, 10, 0, 20, 1,  14, 8,  0, 20, 1,  0,  19, 0,
    20, 2,  17, 8,  1, 20, 2,  29, 25, 0, 20, 2,  4,  19, 0, 20, 4,  5,  19, 0,
    20, 6,  7,  19, 0, 20, 8,  9,  19, 0, 1,  31, 29, 8,  0, 1,  31, 20, 0,  0,
    2,  31, 29, 12, 0, 4,  31, 20, 5,  0, 3,  31, 20, 6,  0, 6,  31, 29, 12, 4,
    8,  31, 20, 9,  0, 28, 29, 8,  19, 0};
  const int nuv[1030] = {
    -2, 1,  0, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 0, 0, -2, 1,  0, 0, 0, -1, 1,  1, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, 1,  1, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0,
    -1, 1,  0, 0, 0, -2, 1,  0, 0, 0, -2, 1,  0, 0, 0, -2, 1,  0, 0, 0,
    -2, 1,  0, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -2, 1,  0, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 0, 0, -1, 1,  1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -2, 1,  1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 2, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -2, 1,  1, 0, 0, -2, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 2, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -2, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 1, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -2, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 2, 0, 0, -1, -1, 2, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 2, 1, 0, -1, -1, 2, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -2, 1,  2, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 1, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 1, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 1, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
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
  amrex::Real c[22]; // temporary storage
  amrex::Real PORT =
    1e6 * P / (8.31446261815324e+07 * T); // 1e6 * P/RT so c goes to SI units

  // Compute conversion, see Eq 10
  for (int id = 0; id < 22; ++id) {
    c[id] = x[id] * PORT;
  }

  // convert to chemkin units
  progressRateFR(q_f, q_r, c, T);

  // convert to chemkin units
  for (int id = 0; id < 206; ++id) {
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
  amrex::Real g_RT[22];
  gibbs(g_RT, T);
  amrex::Real g_RT_qss[10];
  gibbs_qss(g_RT_qss, T);

  amrex::Real sc_qss[10];
  // Fill sc_qss here
  amrex::Real kf_qss[142], qf_qss[142], qr_qss[142];
  comp_k_f_qss(T, invT, logT, kf_qss);
  comp_qss_coeff(kf_qss, qf_qss, qr_qss, sc, T, g_RT, g_RT_qss);
  comp_sc_qss(sc_qss, qf_qss, qr_qss);
  comp_qfqr(q_f, q_r, sc, sc_qss, T, invT, logT);
}

// save atomic weights into array
void
atomicWeight(amrex::Real* awt)
{
  awt[0] = 1.008000;  // H
  awt[1] = 15.999000; // O
  awt[2] = 12.011000; // C
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
  for (int id = 0; id < kd * 22; ++id) {
    ncf[id] = 0;
  }

  // H2
  ncf[0 * kd + 0] = 2; // H

  // H
  ncf[1 * kd + 0] = 1; // H

  // O
  ncf[2 * kd + 1] = 1; // O

  // O2
  ncf[3 * kd + 1] = 2; // O

  // OH
  ncf[4 * kd + 0] = 1; // H
  ncf[4 * kd + 1] = 1; // O

  // H2O
  ncf[5 * kd + 0] = 2; // H
  ncf[5 * kd + 1] = 1; // O

  // HO2
  ncf[6 * kd + 0] = 1; // H
  ncf[6 * kd + 1] = 2; // O

  // H2O2
  ncf[7 * kd + 0] = 2; // H
  ncf[7 * kd + 1] = 2; // O

  // CH3
  ncf[8 * kd + 2] = 1; // C
  ncf[8 * kd + 0] = 3; // H

  // CH4
  ncf[9 * kd + 2] = 1; // C
  ncf[9 * kd + 0] = 4; // H

  // CO
  ncf[10 * kd + 2] = 1; // C
  ncf[10 * kd + 1] = 1; // O

  // CO2
  ncf[11 * kd + 2] = 1; // C
  ncf[11 * kd + 1] = 2; // O

  // CH2O
  ncf[12 * kd + 2] = 1; // C
  ncf[12 * kd + 0] = 2; // H
  ncf[12 * kd + 1] = 1; // O

  // C2H2
  ncf[13 * kd + 2] = 2; // C
  ncf[13 * kd + 0] = 2; // H

  // C2H4
  ncf[14 * kd + 2] = 2; // C
  ncf[14 * kd + 0] = 4; // H

  // C2H6
  ncf[15 * kd + 2] = 2; // C
  ncf[15 * kd + 0] = 6; // H

  // HCCO
  ncf[16 * kd + 2] = 2; // C
  ncf[16 * kd + 0] = 1; // H
  ncf[16 * kd + 1] = 1; // O

  // CH2CO
  ncf[17 * kd + 2] = 2; // C
  ncf[17 * kd + 0] = 2; // H
  ncf[17 * kd + 1] = 1; // O

  // CH3CHO
  ncf[18 * kd + 2] = 2; // C
  ncf[18 * kd + 0] = 4; // H
  ncf[18 * kd + 1] = 1; // O

  // aC3H5
  ncf[19 * kd + 2] = 3; // C
  ncf[19 * kd + 0] = 5; // H

  // C3H6
  ncf[20 * kd + 2] = 3; // C
  ncf[20 * kd + 0] = 6; // H

  // N2
  ncf[21 * kd + 3] = 2; // N
}

// Returns the vector of strings of element names
void
CKSYME_STR(amrex::Vector<std::string>& ename)
{
  ename.resize(4);
  ename[0] = "H";
  ename[1] = "O";
  ename[2] = "C";
  ename[3] = "N";
}

// Returns the vector of strings of species names
void
CKSYMS_STR(amrex::Vector<std::string>& kname)
{
  kname.resize(22);
  kname[0] = "H2";
  kname[1] = "H";
  kname[2] = "O";
  kname[3] = "O2";
  kname[4] = "OH";
  kname[5] = "H2O";
  kname[6] = "HO2";
  kname[7] = "H2O2";
  kname[8] = "CH3";
  kname[9] = "CH4";
  kname[10] = "CO";
  kname[11] = "CO2";
  kname[12] = "CH2O";
  kname[13] = "C2H2";
  kname[14] = "C2H4";
  kname[15] = "C2H6";
  kname[16] = "HCCO";
  kname[17] = "CH2CO";
  kname[18] = "CH3CHO";
  kname[19] = "aC3H5";
  kname[20] = "C3H6";
  kname[21] = "N2";
}

// compute the sparsity pattern of the chemistry Jacobian
void
SPARSITY_INFO(int* nJdata, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 529> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 22> conc = {0.0};
  for (int n = 0; n < 22; n++) {
    conc[n] = 1.0 / 22.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 23; k++) {
    for (int l = 0; l < 23; l++) {
      if (Jac[23 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 529> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 22> conc = {0.0};
  for (int n = 0; n < 22; n++) {
    conc[n] = 1.0 / 22.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 23; k++) {
    for (int l = 0; l < 23; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[23 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 529> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 22> conc = {0.0};
  for (int n = 0; n < 22; n++) {
    conc[n] = 1.0 / 22.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 23; k++) {
    for (int l = 0; l < 23; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[23 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 529> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 22> conc = {0.0};
  for (int n = 0; n < 22; n++) {
    conc[n] = 1.0 / 22.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int nc = 0; nc < NCELLS; nc++) {
    int offset_row = nc * 23;
    int offset_col = nc * 23;
    for (int k = 0; k < 23; k++) {
      for (int l = 0; l < 23; l++) {
        if (Jac[23 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 529> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 22> conc = {0.0};
  for (int n = 0; n < 22; n++) {
    conc[n] = 1.0 / 22.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtrs[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 23;
      for (int l = 0; l < 23; l++) {
        for (int k = 0; k < 23; k++) {
          if (Jac[23 * k + l] != 0.0) {
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
      int offset = nc * 23;
      for (int l = 0; l < 23; l++) {
        for (int k = 0; k < 23; k++) {
          if (Jac[23 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 529> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 22> conc = {0.0};
  for (int n = 0; n < 22; n++) {
    conc[n] = 1.0 / 22.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 23;
      for (int l = 0; l < 23; l++) {
        for (int k = 0; k < 23; k++) {
          if (k == l) {
            colVals[nJdata_tmp - 1] = l + 1 + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[23 * k + l] != 0.0) {
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
      int offset = nc * 23;
      for (int l = 0; l < 23; l++) {
        for (int k = 0; k < 23; k++) {
          if (k == l) {
            colVals[nJdata_tmp] = l + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[23 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 529> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 22> conc = {0.0};
  for (int n = 0; n < 22; n++) {
    conc[n] = 1.0 / 22.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int k = 0; k < 23; k++) {
    for (int l = 0; l < 23; l++) {
      if (k == l) {
        rowVals[nJdata_tmp] = l;
        indx[nJdata_tmp] = 23 * k + l;
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[23 * k + l] != 0.0) {
          rowVals[nJdata_tmp] = l;
          indx[nJdata_tmp] = 23 * k + l;
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
  amrex::GpuArray<amrex::Real, 529> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 22> conc = {0.0};
  for (int n = 0; n < 22; n++) {
    conc[n] = 1.0 / 22.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int l = 0; l < 23; l++) {
      for (int k = 0; k < 23; k++) {
        if (k == l) {
          colVals[nJdata_tmp - 1] = l + 1;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[23 * k + l] != 0.0) {
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
    for (int l = 0; l < 23; l++) {
      for (int k = 0; k < 23; k++) {
        if (k == l) {
          colVals[nJdata_tmp] = l;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[23 * k + l] != 0.0) {
            colVals[nJdata_tmp] = k;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
      }
      rowPtr[l + 1] = nJdata_tmp;
    }
  }
}
