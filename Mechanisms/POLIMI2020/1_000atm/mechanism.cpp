#include "mechanism.H"
const int rmap[NUM_REACTIONS] = {
  8,   21,  47,  140, 182, 201, 187, 0,   3,   5,   7,   22,  60,  163, 177,
  1,   2,   4,   6,   9,   10,  11,  12,  13,  14,  15,  16,  17,  18,  19,
  20,  23,  24,  25,  26,  27,  28,  29,  30,  31,  32,  33,  34,  35,  36,
  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  48,  49,  50,  51,  52,
  53,  54,  55,  56,  57,  58,  59,  61,  62,  63,  64,  65,  66,  67,  68,
  69,  70,  71,  72,  73,  74,  75,  76,  77,  78,  79,  80,  81,  82,  83,
  84,  85,  86,  87,  88,  89,  90,  91,  92,  93,  94,  95,  96,  97,  98,
  99,  100, 101, 102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113,
  114, 115, 116, 117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128,
  129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 141, 142, 143, 144,
  145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155, 156, 157, 158, 159,
  160, 161, 162, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175,
  176, 178, 179, 180, 181, 183, 184, 185, 186, 188, 189, 190, 191, 192, 193,
  194, 195, 196, 197, 198, 199, 200, 202};

// Returns 0-based map of reaction order
void
GET_RMAP(int* _rmap)
{
  for (int j = 0; j < NUM_REACTIONS; ++j) {
    _rmap[j] = rmap[j];
  }
}

// Returns a count of gas species in a gas reaction, and their indices
// and stoichiometric coefficients. (Eq 50)
void
CKINU(const int i, int& nspec, int ki[], int nu[])
{
  const int ns[NUM_GAS_REACTIONS] = {
    2, 4, 4, 2, 4, 3, 3, 3, 2, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 3, 3, 3, 3,
    3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 2, 3,
    3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3,
    4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4,
    4, 4, 4, 2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 2, 3, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 5, 3, 3, 3, 4, 4, 3, 4,
    4, 3, 3, 3, 4, 4, 3, 4, 4, 4, 4, 4, 2, 4, 4, 4, 3, 3, 4};
  const int kiv[NUM_GAS_REACTIONS * 5] = {
    3,  4,  0,  0,  0, 3,  6,  4,  8,  0, 3,  8,  4,  7,  0, 6,  5,  0,  0,  0,
    4,  5,  6,  8,  0, 4,  8,  7,  0,  0, 7,  6,  8,  0,  0, 4,  6,  8,  0,  0,
    9,  8,  0,  0,  0, 4,  9,  7,  8,  0, 4,  9,  3,  10, 0, 9,  6,  10, 8,  0,
    9,  8,  7,  10, 0, 9,  8,  7,  10, 0, 4,  10, 8,  0,  0, 4,  10, 3,  5,  0,
    10, 6,  5,  8,  0, 10, 8,  7,  5,  0, 10, 8,  7,  5,  0, 10, 9,  5,  0,  0,
    10, 9,  5,  0,  0, 4,  5,  10, 0,  0, 6,  8,  10, 0,  0, 22, 4,  28, 0,  0,
    4,  28, 3,  26, 0, 4,  22, 3,  28, 0, 22, 8,  7,  28, 0, 22, 6,  28, 8,  0,
    10, 22, 9,  28, 0, 22, 5,  10, 28, 0, 28, 6,  4,  14, 0, 28, 6,  4,  14, 0,
    28, 6,  26, 8,  0, 28, 6,  26, 8,  0, 28, 8,  7,  26, 0, 28, 5,  14, 8,  0,
    28, 5,  29, 6,  0, 10, 28, 29, 8,  0, 26, 28, 4,  18, 0, 26, 28, 24, 22, 0,
    26, 24, 28, 0,  0, 26, 3,  1,  0,  0, 26, 4,  1,  0,  0, 28, 26, 22, 0,  0,
    28, 23, 0,  0,  0, 28, 4,  30, 0,  0, 28, 3,  19, 0,  0, 20, 28, 8,  0,  0,
    4,  20, 3,  21, 0, 4,  20, 3,  29, 0, 20, 6,  21, 8,  0, 20, 6,  29, 8,  0,
    20, 8,  7,  21, 0, 20, 8,  29, 7,  0, 28, 20, 21, 22, 0, 28, 20, 29, 22, 0,
    26, 20, 21, 28, 0, 26, 20, 29, 28, 0, 10, 20, 9,  21, 0, 10, 20, 29, 9,  0,
    21, 4,  14, 0,  0, 21, 5,  14, 10, 0, 4,  21, 28, 8,  0, 4,  21, 3,  14, 0,
    21, 6,  14, 8,  0, 21, 6,  14, 8,  0, 21, 8,  7,  14, 0, 21, 10, 9,  14, 0,
    21, 10, 20, 5,  0, 21, 28, 14, 22, 0, 21, 28, 30, 8,  0, 21, 28, 19, 7,  0,
    21, 13, 14, 16, 0, 28, 13, 29, 11, 0, 28, 13, 7,  12, 0, 28, 11, 7,  1,  0,
    28, 11, 27, 8,  0, 4,  26, 3,  24, 0, 26, 6,  4,  11, 0, 26, 8,  4,  14, 0,
    26, 8,  7,  24, 0, 26, 5,  14, 6,  0, 26, 5,  11, 8,  0, 24, 26, 4,  1,  0,
    26, 11, 4,  12, 0, 26, 11, 1,  8,  0, 26, 13, 12, 8,  0, 26, 13, 14, 11, 0,
    24, 8,  4,  11, 0, 24, 5,  11, 6,  0, 24, 11, 1,  6,  0, 23, 3,  19, 0,  0,
    4,  23, 3,  30, 0, 23, 6,  30, 8,  0, 23, 6,  7,  18, 0, 23, 8,  7,  30, 0,
    23, 28, 30, 22, 0, 23, 11, 14, 30, 0, 23, 13, 16, 30, 0, 23, 13, 15, 30, 0,
    30, 4,  18, 0,  0, 4,  30, 3,  18, 0, 30, 6,  14, 28, 0, 30, 6,  18, 8,  0,
    30, 8,  7,  18, 0, 30, 8,  19, 7,  0, 30, 28, 18, 22, 0, 30, 28, 19, 22, 0,
    10, 30, 9,  18, 0, 10, 30, 23, 5,  0, 18, 4,  27, 0,  0, 18, 4,  27, 0,  0,
    4,  18, 3,  27, 0, 18, 6,  27, 8,  0, 18, 8,  7,  27, 0, 18, 11, 12, 28, 0,
    18, 26, 28, 27, 0, 18, 28, 22, 27, 0, 18, 19, 0,  0,  0, 19, 4,  27, 0,  0,
    19, 4,  27, 0,  0, 19, 5,  28, 13, 0, 4,  19, 4,  18, 0, 4,  19, 3,  27, 0,
    19, 6,  28, 11, 0, 19, 6,  27, 8,  0, 19, 8,  7,  27, 0, 19, 28, 22, 27, 0,
    19, 10, 9,  27, 0, 27, 4,  1,  0,  0, 4,  27, 3,  1,  0, 27, 6,  4,  12, 0,
    27, 6,  26, 11, 0, 27, 6,  1,  8,  0, 27, 8,  7,  1,  0, 27, 5,  10, 1,  0,
    28, 27, 1,  22, 0, 10, 27, 9,  1,  0, 27, 11, 14, 1,  0, 10, 11, 13, 8,  0,
    11, 6,  13, 0,  0, 11, 8,  16, 0,  0, 14, 4,  11, 0,  0, 4,  14, 3,  11, 0,
    14, 6,  11, 8,  0, 14, 8,  4,  16, 0, 14, 8,  7,  11, 0, 14, 5,  10, 11, 0,
    14, 28, 22, 11, 0, 14, 11, 12, 8,  0, 14, 13, 16, 11, 0, 14, 13, 15, 11, 0,
    14, 7,  12, 0,  0, 14, 10, 15, 8,  0, 4,  16, 3,  13, 0, 4,  15, 3,  13, 0,
    4,  16, 7,  11, 0, 16, 6,  13, 8,  0, 16, 8,  7,  13, 0, 16, 26, 28, 13, 0,
    16, 28, 22, 13, 0, 16, 7,  11, 13, 0, 29, 21, 0,  0,  0, 29, 4,  14, 0,  0,
    4,  29, 3,  14, 0, 4,  29, 28, 8,  0, 29, 6,  14, 8,  0, 29, 8,  7,  14, 0,
    29, 11, 14, 0,  0, 29, 13, 14, 16, 0, 29, 28, 14, 22, 0, 29, 5,  14, 10, 0,
    29, 10, 9,  14, 0, 4,  25, 13, 8,  0, 25, 6,  13, 5,  0, 25, 8,  10, 13, 0,
    10, 25, 13, 5,  8, 25, 11, 5,  0,  0, 25, 11, 5,  0,  0, 25, 13, 5,  0,  0,
    4,  13, 11, 8,  0, 13, 6,  11, 5,  0, 13, 6,  25, 0,  0, 10, 13, 16, 5,  0,
    10, 13, 15, 5,  0, 13, 11, 5,  0,  0, 13, 11, 25, 0,  0, 12, 1,  6,  0,  0,
    4,  12, 1,  8,  0, 4,  12, 1,  8,  0, 12, 6,  11, 0,  0, 12, 6,  1,  5,  0,
    12, 8,  10, 1,  0, 12, 11, 1,  13, 0, 15, 6,  13, 8,  0, 15, 8,  7,  13, 0,
    15, 16, 0,  0,  0, 15, 28, 22, 13, 0, 4,  15, 7,  11, 0, 4,  15, 14, 8,  0,
    11, 8,  15, 0,  0, 13, 8,  17, 0,  0, 17, 8,  7,  25, 0};
  const int nuv[NUM_GAS_REACTIONS * 5] = {
    -1, 2,  0, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -2, 1,  0, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 2, 0, 0, -1, -1, 1, 0, 0,
    -1, 2,  0, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 2, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -2, 1,  1, 0, 0,
    -2, 1,  1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, 1,  1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -2, 1,  1, 0, 0, -2, 1,  1, 0, 0, -2, 2,  1, 0, 0, -2, 1,  1, 0, 0,
    -2, 1,  0, 0, 0, -2, 1,  1, 0, 0, -2, 1,  1, 0, 0, -1, 1,  1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, 1,  1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, 1,  1, 0, 0, -1, 1,  1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, 1,  0, 0, 0, -1, 1,  1, 0, 0,
    -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -2, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -2, 1,  1, 1, 0, -1, 1,  0, 0, 0, -1, 1,  1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 2, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 1, -1, 1,  1, 0, 0, -1, 1,  1, 0, 0, -2, 2,  1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -2, 2,  1, 0, 0, -2, 1,  1, 0, 0, -1, 1,  1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 2, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, 1,  0, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0};
  if (i < 1) {
    // Return max num species per reaction
    nspec = 5;
  } else {
    if (i > NUM_GAS_REACTIONS) {
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
  amrex::Real c[31]; // temporary storage
  amrex::Real PORT =
    1e6 * P / (8.31446261815324e+07 * T); // 1e6 * P/RT so c goes to SI units

  // Compute conversion, see Eq 10
  for (int id = 0; id < 31; ++id) {
    c[id] = x[id] * PORT;
  }

  // convert to chemkin units
  progressRateFR(q_f, q_r, c, T);

  // convert to chemkin units
  for (int id = 0; id < 203; ++id) {
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
  amrex::Real g_RT[31];
  gibbs(g_RT, T);

  amrex::Real sc_qss[1];
  comp_qfqr(q_f, q_r, sc, sc_qss, T, invT, logT);
}

// save atomic weights into array
void
atomicWeight(amrex::Real* awt)
{
  awt[0] = 12.011000; // C
  awt[1] = 1.008000;  // H
  awt[2] = 14.007000; // N
  awt[3] = 15.999000; // O
  awt[4] = 39.950000; // Ar
  awt[5] = 4.002602;  // He
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
  int kd = 6;
  // Zero ncf
  for (int id = 0; id < kd * 31; ++id) {
    ncf[id] = 0;
  }

  // AR
  ncf[0 * kd + 4] = 1; // Ar

  // N2
  ncf[1 * kd + 2] = 2; // N

  // HE
  ncf[2 * kd + 5] = 1; // He

  // H2
  ncf[3 * kd + 1] = 2; // H

  // H
  ncf[4 * kd + 1] = 1; // H

  // O2
  ncf[5 * kd + 3] = 2; // O

  // O
  ncf[6 * kd + 3] = 1; // O

  // H2O
  ncf[7 * kd + 1] = 2; // H
  ncf[7 * kd + 3] = 1; // O

  // OH
  ncf[8 * kd + 1] = 1; // H
  ncf[8 * kd + 3] = 1; // O

  // H2O2
  ncf[9 * kd + 1] = 2; // H
  ncf[9 * kd + 3] = 2; // O

  // HO2
  ncf[10 * kd + 1] = 1; // H
  ncf[10 * kd + 3] = 2; // O

  // NO
  ncf[11 * kd + 2] = 1; // N
  ncf[11 * kd + 3] = 1; // O

  // N2O
  ncf[12 * kd + 2] = 2; // N
  ncf[12 * kd + 3] = 1; // O

  // NO2
  ncf[13 * kd + 2] = 1; // N
  ncf[13 * kd + 3] = 2; // O

  // HNO
  ncf[14 * kd + 1] = 1; // H
  ncf[14 * kd + 2] = 1; // N
  ncf[14 * kd + 3] = 1; // O

  // HNO2
  ncf[15 * kd + 1] = 1; // H
  ncf[15 * kd + 2] = 1; // N
  ncf[15 * kd + 3] = 2; // O

  // HONO
  ncf[16 * kd + 1] = 1; // H
  ncf[16 * kd + 2] = 1; // N
  ncf[16 * kd + 3] = 2; // O

  // HONO2
  ncf[17 * kd + 1] = 1; // H
  ncf[17 * kd + 2] = 1; // N
  ncf[17 * kd + 3] = 3; // O

  // N2H2
  ncf[18 * kd + 1] = 2; // H
  ncf[18 * kd + 2] = 2; // N

  // H2NN
  ncf[19 * kd + 1] = 2; // H
  ncf[19 * kd + 2] = 2; // N

  // NH2OH
  ncf[20 * kd + 1] = 3; // H
  ncf[20 * kd + 2] = 1; // N
  ncf[20 * kd + 3] = 1; // O

  // HNOH
  ncf[21 * kd + 1] = 2; // H
  ncf[21 * kd + 2] = 1; // N
  ncf[21 * kd + 3] = 1; // O

  // NH3
  ncf[22 * kd + 1] = 3; // H
  ncf[22 * kd + 2] = 1; // N

  // N2H4
  ncf[23 * kd + 1] = 4; // H
  ncf[23 * kd + 2] = 2; // N

  // N
  ncf[24 * kd + 2] = 1; // N

  // NO3
  ncf[25 * kd + 2] = 1; // N
  ncf[25 * kd + 3] = 3; // O

  // NH
  ncf[26 * kd + 1] = 1; // H
  ncf[26 * kd + 2] = 1; // N

  // NNH
  ncf[27 * kd + 1] = 1; // H
  ncf[27 * kd + 2] = 2; // N

  // NH2
  ncf[28 * kd + 1] = 2; // H
  ncf[28 * kd + 2] = 1; // N

  // H2NO
  ncf[29 * kd + 1] = 2; // H
  ncf[29 * kd + 2] = 1; // N
  ncf[29 * kd + 3] = 1; // O

  // N2H3
  ncf[30 * kd + 1] = 3; // H
  ncf[30 * kd + 2] = 2; // N
}

// Returns the vector of strings of element names
void
CKSYME_STR(amrex::Vector<std::string>& ename)
{
  ename.resize(6);
  ename[0] = "C";
  ename[1] = "H";
  ename[2] = "N";
  ename[3] = "O";
  ename[4] = "Ar";
  ename[5] = "He";
}

// Returns the vector of strings of species names
void
CKSYMS_STR(amrex::Vector<std::string>& kname)
{
  kname.resize(31);
  kname[0] = "AR";
  kname[1] = "N2";
  kname[2] = "HE";
  kname[3] = "H2";
  kname[4] = "H";
  kname[5] = "O2";
  kname[6] = "O";
  kname[7] = "H2O";
  kname[8] = "OH";
  kname[9] = "H2O2";
  kname[10] = "HO2";
  kname[11] = "NO";
  kname[12] = "N2O";
  kname[13] = "NO2";
  kname[14] = "HNO";
  kname[15] = "HNO2";
  kname[16] = "HONO";
  kname[17] = "HONO2";
  kname[18] = "N2H2";
  kname[19] = "H2NN";
  kname[20] = "NH2OH";
  kname[21] = "HNOH";
  kname[22] = "NH3";
  kname[23] = "N2H4";
  kname[24] = "N";
  kname[25] = "NO3";
  kname[26] = "NH";
  kname[27] = "NNH";
  kname[28] = "NH2";
  kname[29] = "H2NO";
  kname[30] = "N2H3";
}

// compute the sparsity pattern of the chemistry Jacobian
void
SPARSITY_INFO(int* nJdata, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 1024> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 31> conc = {0.0};
  for (int n = 0; n < 31; n++) {
    conc[n] = 1.0 / 31.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 32; k++) {
    for (int l = 0; l < 32; l++) {
      if (Jac[32 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 1024> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 31> conc = {0.0};
  for (int n = 0; n < 31; n++) {
    conc[n] = 1.0 / 31.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 32; k++) {
    for (int l = 0; l < 32; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[32 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 1024> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 31> conc = {0.0};
  for (int n = 0; n < 31; n++) {
    conc[n] = 1.0 / 31.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 32; k++) {
    for (int l = 0; l < 32; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[32 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 1024> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 31> conc = {0.0};
  for (int n = 0; n < 31; n++) {
    conc[n] = 1.0 / 31.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int nc = 0; nc < NCELLS; nc++) {
    int offset_row = nc * 32;
    int offset_col = nc * 32;
    for (int k = 0; k < 32; k++) {
      for (int l = 0; l < 32; l++) {
        if (Jac[32 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 1024> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 31> conc = {0.0};
  for (int n = 0; n < 31; n++) {
    conc[n] = 1.0 / 31.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtrs[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 32;
      for (int l = 0; l < 32; l++) {
        for (int k = 0; k < 32; k++) {
          if (Jac[32 * k + l] != 0.0) {
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
      int offset = nc * 32;
      for (int l = 0; l < 32; l++) {
        for (int k = 0; k < 32; k++) {
          if (Jac[32 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 1024> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 31> conc = {0.0};
  for (int n = 0; n < 31; n++) {
    conc[n] = 1.0 / 31.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 32;
      for (int l = 0; l < 32; l++) {
        for (int k = 0; k < 32; k++) {
          if (k == l) {
            colVals[nJdata_tmp - 1] = l + 1 + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[32 * k + l] != 0.0) {
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
      int offset = nc * 32;
      for (int l = 0; l < 32; l++) {
        for (int k = 0; k < 32; k++) {
          if (k == l) {
            colVals[nJdata_tmp] = l + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[32 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 1024> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 31> conc = {0.0};
  for (int n = 0; n < 31; n++) {
    conc[n] = 1.0 / 31.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int k = 0; k < 32; k++) {
    for (int l = 0; l < 32; l++) {
      if (k == l) {
        rowVals[nJdata_tmp] = l;
        indx[nJdata_tmp] = 32 * k + l;
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[32 * k + l] != 0.0) {
          rowVals[nJdata_tmp] = l;
          indx[nJdata_tmp] = 32 * k + l;
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
  amrex::GpuArray<amrex::Real, 1024> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 31> conc = {0.0};
  for (int n = 0; n < 31; n++) {
    conc[n] = 1.0 / 31.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int l = 0; l < 32; l++) {
      for (int k = 0; k < 32; k++) {
        if (k == l) {
          colVals[nJdata_tmp - 1] = l + 1;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[32 * k + l] != 0.0) {
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
    for (int l = 0; l < 32; l++) {
      for (int k = 0; k < 32; k++) {
        if (k == l) {
          colVals[nJdata_tmp] = l;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[32 * k + l] != 0.0) {
            colVals[nJdata_tmp] = k;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
      }
      rowPtr[l + 1] = nJdata_tmp;
    }
  }
}
