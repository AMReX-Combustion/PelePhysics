#include "mechanism.H"
const int rmap[217] = {
  48,  50,  52,  54,  55,  57,  61,  68,  69,  70,  72,  74,  81,  83,  93,
  129, 138, 144, 155, 171, 180, 195, 203, 209, 211, 11,  0,   1,   32,  33,
  34,  35,  37,  38,  39,  40,  41,  163, 164, 194, 2,   3,   4,   5,   6,
  7,   8,   9,   10,  12,  13,  14,  15,  16,  17,  18,  19,  20,  21,  22,
  23,  24,  25,  26,  27,  28,  29,  30,  31,  36,  42,  43,  44,  45,  46,
  47,  49,  51,  53,  56,  58,  59,  60,  62,  63,  64,  65,  66,  67,  71,
  73,  75,  76,  77,  78,  79,  80,  82,  84,  85,  86,  87,  88,  89,  90,
  91,  92,  94,  95,  96,  97,  98,  99,  100, 101, 102, 103, 104, 105, 106,
  107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121,
  122, 123, 124, 125, 126, 127, 128, 130, 131, 132, 133, 134, 135, 136, 137,
  139, 140, 141, 142, 143, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154,
  156, 157, 158, 159, 160, 161, 162, 165, 166, 167, 168, 169, 170, 172, 173,
  174, 175, 176, 177, 178, 179, 181, 182, 183, 184, 185, 186, 187, 188, 189,
  190, 191, 192, 193, 196, 197, 198, 199, 200, 201, 202, 204, 205, 206, 207,
  208, 210, 212, 213, 214, 215, 216};

// Returns 0-based map of reaction order
void
GET_RMAP(int* _rmap)
{
  for (int j = 0; j < 217; ++j) {
    _rmap[j] = rmap[j];
  }
}

// Returns a count of species in a reaction, and their indices
// and stoichiometric coefficients. (Eq 50)
void
CKINU(const int i, int& nspec, int ki[], int nu[])
{
  const int ns[217] = {
    2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 4, 2, 2, 2, 2, 3, 4, 4, 3, 4, 4, 4, 3, 4,
    3, 4, 3, 4, 3, 3, 4, 3, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 3, 3, 3, 4, 3, 4, 3,
    4, 4, 4, 4, 4, 4, 3, 4, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 3, 4, 4, 4, 5, 4, 3, 4, 3, 3, 4, 4, 5, 4, 4, 3, 4, 4, 3, 4, 4,
    4, 4, 4, 4, 4, 2, 3, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 3, 4, 4, 3,
    5, 4, 4, 4, 4, 3, 4, 4, 3, 4, 4, 4, 4, 5, 5, 4, 5, 5, 5, 3, 3, 5, 5, 4, 4,
    4, 4, 4, 3, 4, 4, 4, 4, 4, 3, 4, 3, 4, 4, 4, 5, 3};
  const int kiv[1085] = {
    2,  3,  0,  0,  0, 1,  2,  4,  0,  0, 0,  2,  1,  4,  0, 6,  2,  3,  4,  0,
    7,  2,  6,  4,  0, 9,  2,  14, 1,  0, 10, 2,  1,  16, 0, 11, 2,  14, 0,  0,
    11, 2,  1,  16, 0, 12, 2,  17, 1,  0, 13, 2,  12, 4,  0, 14, 2,  15, 0,  0,
    16, 2,  14, 4,  0, 16, 2,  15, 1,  0, 17, 2,  16, 4,  0, 18, 2,  17, 4,  0,
    19, 2,  17, 4,  0, 20, 2,  18, 4,  0, 20, 2,  19, 4,  0, 21, 2,  9,  14, 0,
    22, 2,  1,  27, 0, 22, 2,  21, 4,  0, 22, 2,  10, 14, 0, 23, 2,  28, 1,  0,
    24, 2,  12, 16, 0, 25, 2,  17, 12, 0, 26, 2,  25, 4,  0, 27, 2,  14, 1,  0,
    28, 2,  27, 4,  0, 28, 2,  10, 15, 0, 14, 3,  15, 2,  0, 17, 3,  16, 6,  0,
    1,  3,  6,  0,  0, 1,  3,  6,  0,  0, 1,  3,  6,  0,  0, 1,  3,  6,  0,  0,
    1,  3,  2,  4,  0, 1,  0,  0,  0,  0, 1,  0,  0,  0,  0, 1,  0,  0,  0,  0,
    1,  0,  0,  0,  0, 1,  4,  5,  0,  0, 1,  6,  5,  2,  0, 1,  6,  0,  3,  0,
    1,  6,  4,  0,  0, 1,  7,  0,  6,  0, 1,  7,  5,  4,  0, 9,  1,  8,  0,  0,
    10, 1,  12, 0,  0, 11, 1,  9,  0,  0, 12, 1,  13, 0,  0, 13, 1,  12, 0,  0,
    1,  16, 17, 0,  0, 1,  16, 14, 0,  0, 17, 1,  18, 0,  0, 17, 1,  19, 0,  0,
    17, 1,  0,  16, 0, 18, 1,  20, 0,  0, 18, 1,  17, 0,  0, 18, 1,  12, 4,  0,
    18, 1,  11, 5,  0, 19, 1,  20, 0,  0, 19, 1,  18, 1,  0, 19, 1,  17, 0,  0,
    19, 1,  12, 4,  0, 19, 1,  11, 5,  0, 20, 1,  18, 0,  0, 20, 1,  19, 0,  0,
    21, 1,  22, 0,  0, 22, 1,  23, 0,  0, 23, 1,  24, 0,  0, 23, 1,  22, 0,  0,
    24, 1,  25, 0,  0, 24, 1,  23, 0,  0, 25, 1,  26, 0,  0, 25, 1,  24, 0,  0,
    26, 1,  25, 0,  0, 1,  27, 11, 14, 0, 28, 1,  0,  27, 0, 28, 1,  12, 14, 0,
    1,  29, 28, 1,  0, 14, 0,  17, 0,  0, 0,  4,  1,  5,  0, 4,  7,  0,  0,  0,
    4,  5,  2,  0,  0, 6,  4,  5,  3,  0, 7,  4,  5,  6,  0, 7,  4,  5,  6,  0,
    8,  4,  14, 1,  0, 9,  4,  1,  16, 0, 10, 4,  17, 1,  0, 10, 4,  9,  5,  0,
    11, 4,  17, 1,  0, 12, 4,  20, 0,  0, 12, 4,  10, 5,  0, 12, 4,  11, 5,  0,
    13, 4,  12, 5,  0, 14, 4,  15, 1,  0, 16, 4,  14, 5,  0, 17, 4,  5,  16, 0,
    18, 4,  17, 5,  0, 19, 4,  17, 5,  0, 20, 4,  18, 5,  0, 20, 4,  19, 5,  0,
    21, 4,  1,  27, 0, 22, 4,  28, 1,  0, 22, 4,  1,  29, 0, 22, 4,  21, 5,  0,
    22, 4,  12, 14, 0, 23, 4,  22, 5,  0, 24, 4,  23, 5,  0, 26, 4,  25, 5,  0,
    28, 4,  5,  27, 0, 6,  7,  3,  0,  0, 6,  7,  3,  0,  0, 10, 6,  17, 4,  0,
    12, 6,  13, 3,  0, 12, 6,  19, 4,  0, 14, 6,  15, 4,  0, 17, 6,  7,  16, 0,
    8,  3,  14, 2,  0, 8,  10, 21, 1,  0, 8,  12, 22, 1,  0, 9,  3,  16, 2,  0,
    9,  0,  10, 1,  0, 9,  5,  17, 1,  0, 9,  10, 22, 1,  0, 9,  12, 23, 1,  0,
    9,  13, 24, 1,  0, 9,  14, 27, 0,  0, 9,  15, 14, 16, 0, 9,  17, 28, 1,  0,
    9,  27, 22, 14, 0, 10, 3,  14, 1,  4, 10, 0,  12, 1,  0, 10, 22, 0,  0,  0,
    10, 12, 24, 1,  0, 10, 13, 12, 0,  0, 10, 14, 28, 0,  0, 10, 27, 23, 14, 0,
    11, 30, 10, 30, 0, 11, 3,  14, 1,  4, 11, 3,  14, 5,  0, 11, 0,  12, 1,  0,
    11, 5,  20, 0,  0, 11, 5,  10, 5,  0, 11, 12, 24, 1,  0, 11, 13, 12, 0,  0,
    11, 14, 10, 14, 0, 11, 15, 10, 15, 0, 11, 15, 17, 14, 0, 26, 11, 25, 12, 0,
    12, 3,  19, 2,  0, 12, 3,  17, 4,  0, 12, 7,  13, 6,  0, 12, 26, 0,  0,  0,
    12, 25, 1,  0,  0, 12, 16, 13, 14, 0, 17, 12, 13, 16, 0, 12, 20, 18, 13, 0,
    12, 20, 19, 13, 0, 24, 12, 23, 13, 0, 26, 12, 25, 13, 0, 16, 14, 1,  0,  0,
    16, 14, 1,  0,  0, 16, 3,  14, 6,  0, 18, 3,  17, 6,  0, 19, 3,  17, 6,  0,
    21, 3,  14, 16, 0, 21, 0,  22, 1,  0, 23, 3,  17, 16, 0, 24, 22, 0,  0,  0,
    25, 3,  24, 6,  0, 27, 3,  14, 4,  0, 27, 22, 14, 0,  0, 12, 2,  14, 1,  0,
    24, 2,  33, 1,  0, 25, 2,  34, 1,  0, 6,  4,  5,  3,  0, 12, 4,  17, 0,  0,
    9,  0,  12, 0,  0, 10, 3,  15, 1,  0, 10, 3,  17, 2,  0, 10, 22, 1,  0,  0,
    11, 5,  17, 0,  0, 23, 3,  33, 2,  0, 23, 3,  22, 6,  0, 34, 2,  33, 4,  0,
    34, 2,  12, 14, 4, 34, 3,  12, 14, 6, 34, 1,  33, 0,  0, 34, 1,  12, 14, 0,
    34, 4,  12, 14, 5, 34, 6,  12, 14, 7, 34, 13, 14, 0,  0, 28, 1,  33, 0,  0,
    33, 2,  10, 15, 1, 33, 3,  17, 14, 4, 33, 3,  16, 4,  0, 33, 1,  12, 16, 0,
    33, 1,  28, 0,  0, 33, 4,  28, 5,  0, 33, 4,  18, 16, 0, 25, 12, 32, 0,  0,
    32, 2,  31, 4,  0, 32, 1,  31, 0,  0, 32, 4,  31, 5,  0, 31, 7,  32, 6,  0,
    32, 12, 31, 13, 0, 24, 12, 31, 0,  0, 31, 2,  25, 17, 0, 31, 1,  32, 0,  0,
    31, 1,  25, 12, 0, 31, 4,  25, 18, 0, 31, 6,  32, 3,  0, 31, 6,  25, 17, 4,
    31, 12, 25, 0,  0};
  const int nuv[1085] = {
    -2, 1,  0, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 2, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -2, 1,  0, 0, 0, -2, 1,  0, 0, 0, -2, 1,  0, 0, 0,
    -2, 1,  0, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 2, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -2, 1,  0, 0, 0,
    -2, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -2, 1,  1, 0, 0, -2, 1,  1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 1, -1, -1, 1, 1, 0, -2, 1,  1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 2, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 1, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 2, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -2, 1,  0, 0, 0,
    -2, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, 1,  1, 0, 0,
    -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, 1,  1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 2, 1, 0, -2, 1,  2, 0, 0, -1, -1, 1, 1, 1,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 2, 0, -1, -1, 1, 1, 0, -2, 1,  2, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 1, -1, -1, 1, 1, 1, -1, -1, 1, 1, 0, -1, -1, 1, 1, 1,
    -1, -1, 1, 1, 1, -1, -1, 1, 1, 1, -1, 1,  1, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 1, -1, -1, 1, 1, 1, -1, -1, 2, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 1,
    -1, -1, 2, 0, 0};
  if (i < 1) {
    // Return max num species per reaction
    nspec = 5;
  } else {
    if (i > 217) {
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
  amrex::Real c[35]; // temporary storage
  amrex::Real PORT =
    1e6 * P / (8.31446261815324e+07 * T); // 1e6 * P/RT so c goes to SI units

  // Compute conversion, see Eq 10
  for (int id = 0; id < 35; ++id) {
    c[id] = x[id] * PORT;
  }

  // convert to chemkin units
  progressRateFR(q_f, q_r, c, T);

  // convert to chemkin units
  for (int id = 0; id < 217; ++id) {
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
  amrex::Real g_RT[35];
  gibbs(g_RT, T);

  amrex::Real sc_qss[1];
  comp_qfqr(q_f, q_r, sc, sc_qss, T, invT, logT);
}

// save atomic weights into array
void
atomicWeight(amrex::Real* awt)
{
  awt[0] = 15.999000; // O
  awt[1] = 1.008000;  // H
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
  for (int id = 0; id < kd * 35; ++id) {
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

  // C
  ncf[8 * kd + 2] = 1; // C

  // CH
  ncf[9 * kd + 2] = 1; // C
  ncf[9 * kd + 1] = 1; // H

  // CH2
  ncf[10 * kd + 2] = 1; // C
  ncf[10 * kd + 1] = 2; // H

  // CH2(S)
  ncf[11 * kd + 2] = 1; // C
  ncf[11 * kd + 1] = 2; // H

  // CH3
  ncf[12 * kd + 2] = 1; // C
  ncf[12 * kd + 1] = 3; // H

  // CH4
  ncf[13 * kd + 2] = 1; // C
  ncf[13 * kd + 1] = 4; // H

  // CO
  ncf[14 * kd + 2] = 1; // C
  ncf[14 * kd + 0] = 1; // O

  // CO2
  ncf[15 * kd + 2] = 1; // C
  ncf[15 * kd + 0] = 2; // O

  // HCO
  ncf[16 * kd + 2] = 1; // C
  ncf[16 * kd + 1] = 1; // H
  ncf[16 * kd + 0] = 1; // O

  // CH2O
  ncf[17 * kd + 2] = 1; // C
  ncf[17 * kd + 1] = 2; // H
  ncf[17 * kd + 0] = 1; // O

  // CH2OH
  ncf[18 * kd + 2] = 1; // C
  ncf[18 * kd + 1] = 3; // H
  ncf[18 * kd + 0] = 1; // O

  // CH3O
  ncf[19 * kd + 2] = 1; // C
  ncf[19 * kd + 1] = 3; // H
  ncf[19 * kd + 0] = 1; // O

  // CH3OH
  ncf[20 * kd + 2] = 1; // C
  ncf[20 * kd + 1] = 4; // H
  ncf[20 * kd + 0] = 1; // O

  // C2H
  ncf[21 * kd + 2] = 2; // C
  ncf[21 * kd + 1] = 1; // H

  // C2H2
  ncf[22 * kd + 2] = 2; // C
  ncf[22 * kd + 1] = 2; // H

  // C2H3
  ncf[23 * kd + 2] = 2; // C
  ncf[23 * kd + 1] = 3; // H

  // C2H4
  ncf[24 * kd + 2] = 2; // C
  ncf[24 * kd + 1] = 4; // H

  // C2H5
  ncf[25 * kd + 2] = 2; // C
  ncf[25 * kd + 1] = 5; // H

  // C2H6
  ncf[26 * kd + 2] = 2; // C
  ncf[26 * kd + 1] = 6; // H

  // HCCO
  ncf[27 * kd + 2] = 2; // C
  ncf[27 * kd + 1] = 1; // H
  ncf[27 * kd + 0] = 1; // O

  // CH2CO
  ncf[28 * kd + 2] = 2; // C
  ncf[28 * kd + 1] = 2; // H
  ncf[28 * kd + 0] = 1; // O

  // HCCOH
  ncf[29 * kd + 2] = 2; // C
  ncf[29 * kd + 1] = 2; // H
  ncf[29 * kd + 0] = 1; // O

  // N2
  ncf[30 * kd + 3] = 2; // N

  // C3H7
  ncf[31 * kd + 2] = 3; // C
  ncf[31 * kd + 1] = 7; // H

  // C3H8
  ncf[32 * kd + 2] = 3; // C
  ncf[32 * kd + 1] = 8; // H

  // CH2CHO
  ncf[33 * kd + 2] = 2; // C
  ncf[33 * kd + 1] = 3; // H
  ncf[33 * kd + 0] = 1; // O

  // CH3CHO
  ncf[34 * kd + 2] = 2; // C
  ncf[34 * kd + 1] = 4; // H
  ncf[34 * kd + 0] = 1; // O
}

// Returns the vector of strings of element names
void
CKSYME_STR(amrex::Vector<std::string>& ename)
{
  ename.resize(4);
  ename[0] = "O";
  ename[1] = "H";
  ename[2] = "C";
  ename[3] = "N";
}

// Returns the vector of strings of species names
void
CKSYMS_STR(amrex::Vector<std::string>& kname)
{
  kname.resize(35);
  kname[0] = "H2";
  kname[1] = "H";
  kname[2] = "O";
  kname[3] = "O2";
  kname[4] = "OH";
  kname[5] = "H2O";
  kname[6] = "HO2";
  kname[7] = "H2O2";
  kname[8] = "C";
  kname[9] = "CH";
  kname[10] = "CH2";
  kname[11] = "CH2(S)";
  kname[12] = "CH3";
  kname[13] = "CH4";
  kname[14] = "CO";
  kname[15] = "CO2";
  kname[16] = "HCO";
  kname[17] = "CH2O";
  kname[18] = "CH2OH";
  kname[19] = "CH3O";
  kname[20] = "CH3OH";
  kname[21] = "C2H";
  kname[22] = "C2H2";
  kname[23] = "C2H3";
  kname[24] = "C2H4";
  kname[25] = "C2H5";
  kname[26] = "C2H6";
  kname[27] = "HCCO";
  kname[28] = "CH2CO";
  kname[29] = "HCCOH";
  kname[30] = "N2";
  kname[31] = "C3H7";
  kname[32] = "C3H8";
  kname[33] = "CH2CHO";
  kname[34] = "CH3CHO";
}

// compute the sparsity pattern of the chemistry Jacobian
void
SPARSITY_INFO(int* nJdata, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 1296> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 35> conc = {0.0};
  for (int n = 0; n < 35; n++) {
    conc[n] = 1.0 / 35.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 36; k++) {
    for (int l = 0; l < 36; l++) {
      if (Jac[36 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 1296> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 35> conc = {0.0};
  for (int n = 0; n < 35; n++) {
    conc[n] = 1.0 / 35.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 36; k++) {
    for (int l = 0; l < 36; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[36 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 1296> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 35> conc = {0.0};
  for (int n = 0; n < 35; n++) {
    conc[n] = 1.0 / 35.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 36; k++) {
    for (int l = 0; l < 36; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[36 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 1296> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 35> conc = {0.0};
  for (int n = 0; n < 35; n++) {
    conc[n] = 1.0 / 35.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int nc = 0; nc < NCELLS; nc++) {
    int offset_row = nc * 36;
    int offset_col = nc * 36;
    for (int k = 0; k < 36; k++) {
      for (int l = 0; l < 36; l++) {
        if (Jac[36 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 1296> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 35> conc = {0.0};
  for (int n = 0; n < 35; n++) {
    conc[n] = 1.0 / 35.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtrs[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 36;
      for (int l = 0; l < 36; l++) {
        for (int k = 0; k < 36; k++) {
          if (Jac[36 * k + l] != 0.0) {
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
      int offset = nc * 36;
      for (int l = 0; l < 36; l++) {
        for (int k = 0; k < 36; k++) {
          if (Jac[36 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 1296> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 35> conc = {0.0};
  for (int n = 0; n < 35; n++) {
    conc[n] = 1.0 / 35.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 36;
      for (int l = 0; l < 36; l++) {
        for (int k = 0; k < 36; k++) {
          if (k == l) {
            colVals[nJdata_tmp - 1] = l + 1 + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[36 * k + l] != 0.0) {
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
      int offset = nc * 36;
      for (int l = 0; l < 36; l++) {
        for (int k = 0; k < 36; k++) {
          if (k == l) {
            colVals[nJdata_tmp] = l + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[36 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 1296> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 35> conc = {0.0};
  for (int n = 0; n < 35; n++) {
    conc[n] = 1.0 / 35.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int k = 0; k < 36; k++) {
    for (int l = 0; l < 36; l++) {
      if (k == l) {
        rowVals[nJdata_tmp] = l;
        indx[nJdata_tmp] = 36 * k + l;
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[36 * k + l] != 0.0) {
          rowVals[nJdata_tmp] = l;
          indx[nJdata_tmp] = 36 * k + l;
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
  amrex::GpuArray<amrex::Real, 1296> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 35> conc = {0.0};
  for (int n = 0; n < 35; n++) {
    conc[n] = 1.0 / 35.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int l = 0; l < 36; l++) {
      for (int k = 0; k < 36; k++) {
        if (k == l) {
          colVals[nJdata_tmp - 1] = l + 1;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[36 * k + l] != 0.0) {
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
    for (int l = 0; l < 36; l++) {
      for (int k = 0; k < 36; k++) {
        if (k == l) {
          colVals[nJdata_tmp] = l;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[36 * k + l] != 0.0) {
            colVals[nJdata_tmp] = k;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
      }
      rowPtr[l + 1] = nJdata_tmp;
    }
  }
}
