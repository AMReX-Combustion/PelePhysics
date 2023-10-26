#include "mechanism.H"
const int rmap[224] = {
  6,   29,  35,  51,  75,  83,  93,  101, 102, 108, 113, 149, 32,  2,   3,
  60,  61,  0,   1,   4,   5,   7,   8,   9,   10,  11,  12,  13,  14,  15,
  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  30,  31,
  33,  34,  36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,  48,
  49,  50,  52,  53,  54,  55,  56,  57,  58,  59,  62,  63,  64,  65,  66,
  67,  68,  69,  70,  71,  72,  73,  74,  76,  77,  78,  79,  80,  81,  82,
  84,  85,  86,  87,  88,  89,  90,  91,  92,  94,  95,  96,  97,  98,  99,
  100, 103, 104, 105, 106, 107, 109, 110, 111, 112, 114, 115, 116, 117, 118,
  119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133,
  134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148,
  150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164,
  165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179,
  180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194,
  195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209,
  210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223};

// Returns 0-based map of reaction order
void
GET_RMAP(int* _rmap)
{
  for (int j = 0; j < 224; ++j) {
    _rmap[j] = rmap[j];
  }
}

// Returns a count of species in a reaction, and their indices
// and stoichiometric coefficients. (Eq 50)
void
CKINU(const int i, int& nspec, int ki[], int nu[])
{
  const int ns[224] = {
    4, 4, 3, 3, 4, 3, 3, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    5, 4, 5, 4, 4, 4, 4, 3, 4, 5, 4, 4, 4, 5, 4, 5, 4, 4, 4, 3, 3, 4, 4, 4, 4,
    3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    3, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 5, 4, 4, 4, 4, 3, 4, 5, 4, 4, 4, 4,
    4, 3, 3, 4, 4, 4, 4, 5, 2, 4, 4, 4, 4, 3, 4, 4, 4, 3, 4, 4, 3, 3, 3, 4, 5,
    4, 5, 2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 2, 4, 4, 4, 4, 3, 4, 4, 4, 5, 4, 3,
    3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 3, 5, 6, 5, 6, 3, 4, 4, 4, 4, 4, 3, 3,
    3, 3, 3, 5, 5, 4, 4, 5, 3, 3, 4, 4, 4, 5, 5, 3, 4, 4, 5, 4, 4, 3, 3, 4, 4,
    4, 3, 4, 3, 3, 3, 4, 4, 4, 4, 4, 3, 5, 5, 3, 4, 3, 4, 4, 3, 4, 3, 4, 4};
  const int kiv[1344] = {
    0,  1,  0,  2,  0,  0,  4,  3,  5,  6,  0,  0,  5,  3,  6,  0,  0,  0,  5,
    6,  7,  0,  0,  0,  4,  6,  5,  7,  0,  0,  6,  7,  3,  0,  0,  0,  5,  8,
    9,  0,  0,  0,  5,  8,  3,  6,  0,  0,  5,  9,  4,  8,  0,  0,  5,  9,  7,
    3,  0,  0,  5,  9,  6,  0,  0,  0,  9,  3,  8,  6,  0,  0,  9,  6,  7,  8,
    0,  0,  9,  6,  7,  8,  0,  0,  10, 3,  11, 5,  0,  0,  10, 6,  5,  12, 0,
    0,  10, 4,  5,  2,  0,  0,  10, 7,  13, 5,  0,  0,  10, 8,  12, 3,  0,  0,
    3,  2,  5,  12, 0,  0,  6,  2,  13, 5,  0,  0,  6,  2,  10, 7,  0,  0,  4,
    2,  14, 5,  0,  0,  8,  2,  15, 5,  0,  0,  8,  2,  13, 3,  0,  0,  8,  2,
    11, 5,  6,  0,  4,  1,  14, 5,  0,  0,  8,  1,  11, 5,  6,  0,  8,  1,  11,
    7,  0,  0,  7,  1,  13, 4,  0,  0,  7,  1,  7,  2,  0,  0,  7,  1,  13, 4,
    0,  0,  14, 5,  16, 0,  0,  0,  14, 3,  13, 5,  0,  0,  14, 3,  11, 5,  4,
    0,  14, 6,  13, 4,  0,  0,  14, 6,  7,  2,  0,  0,  14, 6,  7,  1,  0,  0,
    14, 8,  13, 5,  3,  0,  14, 8,  13, 6,  0,  0,  14, 9,  13, 5,  6,  0,  14,
    9,  16, 8,  0,  0,  10, 14, 17, 5,  0,  0,  14, 2,  18, 5,  0,  0,  14, 19,
    5,  0,  0,  0,  14, 16, 1,  0,  0,  0,  16, 5,  14, 4,  0,  0,  16, 3,  14,
    6,  0,  0,  16, 6,  14, 7,  0,  0,  10, 16, 18, 5,  0,  0,  16, 2,  14, 0,
    0,  0,  11, 3,  15, 0,  0,  0,  11, 6,  15, 5,  0,  0,  11, 6,  15, 5,  0,
    0,  11, 9,  15, 6,  0,  0,  11, 1,  11, 2,  0,  0,  5,  12, 11, 4,  0,  0,
    12, 3,  11, 6,  0,  0,  12, 3,  15, 5,  0,  0,  12, 6,  11, 7,  0,  0,  12,
    11, 5,  0,  0,  0,  12, 11, 5,  0,  0,  0,  12, 8,  11, 9,  0,  0,  14, 12,
    16, 11, 0,  0,  13, 5,  4,  12, 0,  0,  13, 3,  12, 6,  0,  0,  13, 6,  7,
    12, 0,  0,  13, 14, 16, 12, 0,  0,  10, 15, 11, 12, 0,  0,  15, 1,  15, 2,
    0,  0,  15, 1,  13, 11, 0,  0,  20, 3,  10, 11, 0,  0,  20, 6,  5,  21, 0,
    0,  20, 8,  11, 12, 0,  0,  20, 4,  22, 5,  0,  0,  22, 5,  17, 0,  0,  0,
    22, 3,  5,  21, 0,  0,  22, 3,  11, 2,  0,  0,  22, 3,  20, 6,  0,  0,  22,
    6,  20, 7,  0,  0,  22, 6,  14, 11, 0,  0,  22, 6,  14, 11, 0,  0,  22, 1,
    23, 5,  0,  0,  17, 5,  18, 0,  0,  0,  17, 5,  22, 4,  0,  0,  17, 3,  14,
    11, 0,  0,  17, 6,  22, 7,  0,  0,  17, 8,  22, 9,  0,  0,  17, 8,  14, 11,
    3,  0,  17, 8,  13, 12, 0,  0,  17, 12, 18, 11, 0,  0,  17, 14, 22, 16, 0,
    0,  17, 14, 24, 5,  0,  0,  18, 5,  19, 0,  0,  0,  18, 5,  17, 4,  0,  0,
    18, 3,  14, 11, 5,  0,  18, 3,  13, 2,  0,  0,  18, 3,  14, 12, 0,  0,  18,
    6,  17, 7,  0,  0,  18, 6,  13, 14, 0,  0,  18, 14, 17, 16, 0,  0,  18, 14,
    25, 0,  0,  0,  19, 5,  26, 0,  0,  0,  19, 3,  13, 14, 0,  0,  19, 8,  18,
    9,  0,  0,  19, 12, 26, 11, 0,  0,  19, 9,  26, 8,  0,  0,  19, 9,  13, 14,
    6,  0,  26, 14, 0,  0,  0,  0,  26, 5,  19, 4,  0,  0,  26, 3,  19, 6,  0,
    0,  26, 6,  19, 7,  0,  0,  26, 14, 19, 16, 0,  0,  21, 10, 11, 0,  0,  0,
    5,  21, 11, 1,  0,  0,  21, 3,  11, 5,  0,  0,  21, 8,  11, 6,  0,  0,  21,
    22, 11, 0,  0,  0,  22, 21, 23, 11, 0,  0,  14, 21, 18, 11, 0,  0,  21, 6,
    12, 0,  0,  0,  23, 5,  27, 0,  0,  0,  23, 5,  28, 0,  0,  0,  23, 6,  18,
    11, 0,  0,  23, 3,  22, 11, 5,  0,  23, 8,  14, 11, 0,  0,  23, 9,  17, 11,
    6,  0,  23, 29, 0,  0,  0,  0,  23, 30, 5,  0,  0,  0,  22, 23, 31, 0,  0,
    0,  5,  27, 22, 14, 0,  0,  5,  27, 23, 4,  0,  0,  6,  27, 23, 7,  0,  0,
    14, 27, 23, 16, 0,  0,  3,  27, 14, 21, 0,  0,  3,  27, 18, 11, 0,  0,  6,
    27, 19, 11, 0,  0,  28, 5,  22, 14, 0,  0,  28, 27, 0,  0,  0,  0,  28, 5,
    5,  27, 0,  0,  28, 5,  23, 4,  0,  0,  28, 6,  23, 7,  0,  0,  28, 14, 23,
    16, 0,  0,  28, 5,  24, 0,  0,  0,  24, 5,  28, 4,  0,  0,  24, 6,  28, 7,
    0,  0,  24, 14, 28, 16, 0,  0,  24, 8,  22, 13, 6,  0,  24, 23, 29, 5,  0,
    0,  32, 5,  25, 0,  0,  0,  32, 17, 14, 0,  0,  0,  24, 5,  32, 0,  0,  0,
    32, 5,  18, 14, 0,  0,  32, 3,  14, 11, 0,  0,  32, 3,  19, 12, 0,  0,  32,
    5,  24, 4,  0,  0,  32, 3,  24, 6,  0,  0,  32, 6,  24, 7,  0,  0,  32, 14,
    24, 16, 0,  0,  33, 5,  18, 19, 0,  0,  33, 24, 14, 0,  0,  0,  31, 9,  34,
    8,  0,  0,  31, 35, 5,  0,  0,  0,  31, 3,  22, 17, 11, 0,  31, 9,  22, 17,
    11, 6,  31, 9,  22, 11, 7,  0,  31, 6,  22, 17, 11, 5,  34, 31, 5,  0,  0,
    0,  34, 5,  31, 4,  0,  0,  34, 5,  24, 22, 0,  0,  34, 3,  31, 6,  0,  0,
    34, 6,  31, 7,  0,  0,  34, 14, 31, 16, 0,  0,  36, 5,  37, 0,  0,  0,  36,
    24, 19, 0,  0,  0,  37, 18, 25, 0,  0,  0,  37, 36, 5,  0,  0,  0,  37, 19,
    32, 0,  0,  0,  30, 8,  31, 11, 3,  0,  30, 8,  22, 11, 5,  0,  30, 3,  31,
    11, 0,  0,  30, 6,  34, 11, 0,  0,  30, 9,  31, 11, 6,  0,  30, 22, 38, 0,
    0,  0,  29, 30, 5,  0,  0,  0,  29, 5,  30, 4,  0,  0,  29, 6,  30, 7,  0,
    0,  29, 14, 30, 16, 0,  0,  29, 3,  31, 11, 5,  0,  29, 6,  34, 11, 5,  0,
    39, 22, 31, 0,  0,  0,  39, 5,  30, 14, 0,  0,  39, 3,  40, 5,  0,  0,  39,
    9,  40, 5,  6,  0,  39, 23, 35, 5,  0,  0,  41, 5,  29, 14, 0,  0,  41, 39,
    5,  0,  0,  0,  41, 30, 14, 0,  0,  0,  41, 5,  39, 4,  0,  0,  41, 6,  39,
    7,  0,  0,  41, 14, 39, 16, 0,  0,  42, 18, 37, 0,  0,  0,  42, 18, 19, 32,
    0,  0,  42, 33, 25, 0,  0,  0,  42, 19, 36, 0,  0,  0,  43, 19, 37, 0,  0,
    0,  43, 18, 19, 25, 0,  0,  5,  43, 42, 4,  0,  0,  43, 3,  42, 6,  0,  0,
    43, 6,  42, 7,  0,  0,  14, 43, 42, 16, 0,  0,  40, 30, 12, 0,  0,  0,  40,
    5,  30, 11, 4,  0,  40, 14, 30, 16, 11, 0,  44, 22, 45, 0,  0,  0,  44, 18,
    35, 5,  0,  0,  46, 44, 5,  0,  0,  0,  46, 5,  44, 4,  0,  0,  46, 6,  44,
    7,  0,  0,  38, 46, 5,  0,  0,  0,  38, 22, 35, 5,  0,  0,  35, 45, 5,  0,
    0,  0,  35, 5,  45, 4,  0,  0,  35, 6,  45, 7,  0,  0};
  const int nuv[1344] = {
    -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  0,  0,  0,  -1,
    -1, 1,  0,  0,  0,  -1, -1, 1,  1,  0,  0,  -2, 1,  1,  0,  0,  0,  -1, -1,
    1,  0,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,
    1,  0,  0,  -1, -1, 2,  0,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,
    0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,
    0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,
    -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1,
    -1, 1,  1,  0,  0,  -1, -1, 1,  2,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1,
    1,  1,  1,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  1,  0,  -1, -1, 1,
    1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,
    0,  0,  -1, -1, 1,  0,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  1,
    0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,
    -1, -1, 1,  1,  1,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  1,  0,  -1,
    -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -2, 1,
    1,  0,  0,  0,  -2, 1,  1,  0,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,
    1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 2,  0,
    0,  0,  -1, -1, 1,  0,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,
    0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,
    -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1,
    1,  1,  0,  0,  0,  -1, 1,  1,  0,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1,
    1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,
    1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,
    0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,
    0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  0,  0,  0,
    -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1,
    -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1,
    1,  1,  0,  0,  -1, -1, 1,  0,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,
    1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,
    1,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,
    0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  0,  0,  0,  -1, -1, 1,  1,  0,  0,
    -1, -1, 1,  1,  1,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1,
    -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1,
    1,  0,  0,  0,  -1, -1, 1,  0,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,
    1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,
    1,  0,  -1, 2,  0,  0,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,
    0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, 1,  1,  0,  0,  0,
    -1, -1, 1,  1,  0,  0,  -1, -1, 2,  1,  0,  0,  -1, -1, 2,  1,  0,  0,  -2,
    1,  2,  0,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1,
    2,  0,  0,  0,  -1, -1, 1,  0,  0,  0,  -1, -1, 1,  0,  0,  0,  -1, -1, 1,
    1,  0,  0,  -1, -1, 1,  1,  1,  0,  -1, -1, 1,  2,  0,  0,  -1, -1, 1,  1,
    1,  0,  -2, 1,  0,  0,  0,  0,  -2, 1,  1,  0,  0,  0,  -1, -1, 1,  0,  0,
    0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,
    -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1,
    -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, 1,  0,  0,  0,  0,  -1, -1,
    1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,
    1,  0,  0,  -1, -1, 1,  0,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,
    0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  1,  0,  -1, -1, 1,  2,  0,
    0,  -1, -1, 1,  0,  0,  0,  -1, 1,  1,  0,  0,  0,  -1, -1, 1,  0,  0,  0,
    -1, -1, 1,  1,  0,  0,  -1, -1, 2,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1,
    -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1,
    1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, 1,  1,  0,  0,  0,  -1, -1, 1,
    1,  0,  0,  -2, 1,  2,  0,  0,  0,  -1, -1, 1,  1,  1,  0,  -1, -1, 1,  1,
    1,  1,  -1, -1, 2,  1,  1,  0,  -1, -1, 1,  1,  1,  1,  -1, 1,  1,  0,  0,
    0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,
    -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  0,  0,  0,  -1,
    1,  1,  0,  0,  0,  -1, 1,  1,  0,  0,  0,  -1, 1,  1,  0,  0,  0,  -1, 1,
    1,  0,  0,  0,  -1, -1, 1,  1,  1,  0,  -1, -1, 2,  2,  1,  0,  -1, -1, 1,
    1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  1,  0,  -1, -1, 1,  0,
    0,  0,  -1, 1,  1,  0,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,
    0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  1,  0,  -1, -1, 1,  1,  1,  0,
    -1, 1,  1,  0,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1,
    -1, 1,  1,  1,  0,  -1, -1, 1,  2,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, 1,
    1,  0,  0,  0,  -1, 1,  1,  0,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,
    1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, 1,  1,  0,  0,  0,  -1, 1,  1,  1,
    0,  0,  -1, 1,  1,  0,  0,  0,  -1, 1,  1,  0,  0,  0,  -1, 1,  1,  0,  0,
    0,  -1, 1,  1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,
    -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, 1,  1,  0,  0,  0,  -1,
    -1, 1,  1,  1,  0,  -1, -1, 1,  1,  1,  0,  -1, -1, 1,  0,  0,  0,  -1, -1,
    1,  1,  0,  0,  -1, 1,  1,  0,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,
    1,  0,  0,  -1, 1,  1,  0,  0,  0,  -1, -1, 1,  1,  0,  0,  -1, 1,  1,  0,
    0,  0,  -1, -1, 1,  1,  0,  0,  -1, -1, 1,  1,  0,  0};
  if (i < 1) {
    // Return max num species per reaction
    nspec = 6;
  } else {
    if (i > 224) {
      nspec = -1;
    } else {
      nspec = ns[i - 1];
      for (int j = 0; j < nspec; ++j) {
        ki[j] = kiv[(i - 1) * 6 + j] + 1;
        nu[j] = nuv[(i - 1) * 6 + j];
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
  amrex::Real c[47]; // temporary storage
  amrex::Real PORT =
    1e6 * P / (8.31446261815324e+07 * T); // 1e6 * P/RT so c goes to SI units

  // Compute conversion, see Eq 10
  for (int id = 0; id < 47; ++id) {
    c[id] = x[id] * PORT;
  }

  // convert to chemkin units
  progressRateFR(q_f, q_r, c, T);

  // convert to chemkin units
  for (int id = 0; id < 224; ++id) {
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
  amrex::Real g_RT[47];
  gibbs(g_RT, T);

  amrex::Real sc_qss[1];
  comp_qfqr(q_f, q_r, sc, sc_qss, T, invT, logT);
}

// save atomic weights into array
void
atomicWeight(amrex::Real* awt)
{
  awt[0] = 14.007000; // N
  awt[1] = 12.011000; // C
  awt[2] = 1.008000;  // H
  awt[3] = 15.999000; // O
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
  for (int id = 0; id < kd * 47; ++id) {
    ncf[id] = 0;
  }

  // N2
  ncf[0 * kd + 0] = 2; // N

  // S-CH2
  ncf[1 * kd + 1] = 1; // C
  ncf[1 * kd + 2] = 2; // H

  // T-CH2
  ncf[2 * kd + 1] = 1; // C
  ncf[2 * kd + 2] = 2; // H

  // O
  ncf[3 * kd + 3] = 1; // O

  // H2
  ncf[4 * kd + 2] = 2; // H

  // H
  ncf[5 * kd + 2] = 1; // H

  // OH
  ncf[6 * kd + 2] = 1; // H
  ncf[6 * kd + 3] = 1; // O

  // H2O
  ncf[7 * kd + 2] = 2; // H
  ncf[7 * kd + 3] = 1; // O

  // O2
  ncf[8 * kd + 3] = 2; // O

  // HO2
  ncf[9 * kd + 2] = 1; // H
  ncf[9 * kd + 3] = 2; // O

  // CH
  ncf[10 * kd + 1] = 1; // C
  ncf[10 * kd + 2] = 1; // H

  // CO
  ncf[11 * kd + 1] = 1; // C
  ncf[11 * kd + 3] = 1; // O

  // HCO
  ncf[12 * kd + 1] = 1; // C
  ncf[12 * kd + 2] = 1; // H
  ncf[12 * kd + 3] = 1; // O

  // CH2O
  ncf[13 * kd + 1] = 1; // C
  ncf[13 * kd + 2] = 2; // H
  ncf[13 * kd + 3] = 1; // O

  // CH3
  ncf[14 * kd + 1] = 1; // C
  ncf[14 * kd + 2] = 3; // H

  // CO2
  ncf[15 * kd + 1] = 1; // C
  ncf[15 * kd + 3] = 2; // O

  // CH4
  ncf[16 * kd + 1] = 1; // C
  ncf[16 * kd + 2] = 4; // H

  // C2H3
  ncf[17 * kd + 1] = 2; // C
  ncf[17 * kd + 2] = 3; // H

  // C2H4
  ncf[18 * kd + 1] = 2; // C
  ncf[18 * kd + 2] = 4; // H

  // C2H5
  ncf[19 * kd + 1] = 2; // C
  ncf[19 * kd + 2] = 5; // H

  // C2H
  ncf[20 * kd + 1] = 2; // C
  ncf[20 * kd + 2] = 1; // H

  // HCCO
  ncf[21 * kd + 1] = 2; // C
  ncf[21 * kd + 2] = 1; // H
  ncf[21 * kd + 3] = 1; // O

  // C2H2
  ncf[22 * kd + 1] = 2; // C
  ncf[22 * kd + 2] = 2; // H

  // C3H3
  ncf[23 * kd + 1] = 3; // C
  ncf[23 * kd + 2] = 3; // H

  // A-C3H5
  ncf[24 * kd + 1] = 3; // C
  ncf[24 * kd + 2] = 5; // H

  // N-C3H7
  ncf[25 * kd + 1] = 3; // C
  ncf[25 * kd + 2] = 7; // H

  // C2H6
  ncf[26 * kd + 1] = 2; // C
  ncf[26 * kd + 2] = 6; // H

  // P-C3H4
  ncf[27 * kd + 1] = 3; // C
  ncf[27 * kd + 2] = 4; // H

  // A-C3H4
  ncf[28 * kd + 1] = 3; // C
  ncf[28 * kd + 2] = 4; // H

  // A1
  ncf[29 * kd + 1] = 6; // C
  ncf[29 * kd + 2] = 6; // H

  // A1-
  ncf[30 * kd + 1] = 6; // C
  ncf[30 * kd + 2] = 5; // H

  // C5H5
  ncf[31 * kd + 1] = 5; // C
  ncf[31 * kd + 2] = 5; // H

  // C3H6
  ncf[32 * kd + 1] = 3; // C
  ncf[32 * kd + 2] = 6; // H

  // C4H8
  ncf[33 * kd + 1] = 4; // C
  ncf[33 * kd + 2] = 8; // H

  // C5H6
  ncf[34 * kd + 1] = 5; // C
  ncf[34 * kd + 2] = 6; // H

  // A2
  ncf[35 * kd + 1] = 10; // C
  ncf[35 * kd + 2] = 8;  // H

  // C5H10
  ncf[36 * kd + 1] = 5;  // C
  ncf[36 * kd + 2] = 10; // H

  // C5H11
  ncf[37 * kd + 1] = 5;  // C
  ncf[37 * kd + 2] = 11; // H

  // A1C2H2
  ncf[38 * kd + 1] = 8; // C
  ncf[38 * kd + 2] = 7; // H

  // A1CH2
  ncf[39 * kd + 1] = 7; // C
  ncf[39 * kd + 2] = 7; // H

  // A1CHO
  ncf[40 * kd + 1] = 7; // C
  ncf[40 * kd + 2] = 6; // H
  ncf[40 * kd + 3] = 1; // O

  // A1CH3
  ncf[41 * kd + 1] = 7; // C
  ncf[41 * kd + 2] = 8; // H

  // C7H15
  ncf[42 * kd + 1] = 7;  // C
  ncf[42 * kd + 2] = 15; // H

  // N-C7H16
  ncf[43 * kd + 1] = 7;  // C
  ncf[43 * kd + 2] = 16; // H

  // A1C2H*
  ncf[44 * kd + 1] = 8; // C
  ncf[44 * kd + 2] = 5; // H

  // A2-
  ncf[45 * kd + 1] = 10; // C
  ncf[45 * kd + 2] = 7;  // H

  // A1C2H
  ncf[46 * kd + 1] = 8; // C
  ncf[46 * kd + 2] = 6; // H
}

// Returns the vector of strings of element names
void
CKSYME_STR(amrex::Vector<std::string>& ename)
{
  ename.resize(4);
  ename[0] = "N";
  ename[1] = "C";
  ename[2] = "H";
  ename[3] = "O";
}

// Returns the vector of strings of species names
void
CKSYMS_STR(amrex::Vector<std::string>& kname)
{
  kname.resize(47);
  kname[0] = "N2";
  kname[1] = "S-CH2";
  kname[2] = "T-CH2";
  kname[3] = "O";
  kname[4] = "H2";
  kname[5] = "H";
  kname[6] = "OH";
  kname[7] = "H2O";
  kname[8] = "O2";
  kname[9] = "HO2";
  kname[10] = "CH";
  kname[11] = "CO";
  kname[12] = "HCO";
  kname[13] = "CH2O";
  kname[14] = "CH3";
  kname[15] = "CO2";
  kname[16] = "CH4";
  kname[17] = "C2H3";
  kname[18] = "C2H4";
  kname[19] = "C2H5";
  kname[20] = "C2H";
  kname[21] = "HCCO";
  kname[22] = "C2H2";
  kname[23] = "C3H3";
  kname[24] = "A-C3H5";
  kname[25] = "N-C3H7";
  kname[26] = "C2H6";
  kname[27] = "P-C3H4";
  kname[28] = "A-C3H4";
  kname[29] = "A1";
  kname[30] = "A1-";
  kname[31] = "C5H5";
  kname[32] = "C3H6";
  kname[33] = "C4H8";
  kname[34] = "C5H6";
  kname[35] = "A2";
  kname[36] = "C5H10";
  kname[37] = "C5H11";
  kname[38] = "A1C2H2";
  kname[39] = "A1CH2";
  kname[40] = "A1CHO";
  kname[41] = "A1CH3";
  kname[42] = "C7H15";
  kname[43] = "N-C7H16";
  kname[44] = "A1C2H*";
  kname[45] = "A2-";
  kname[46] = "A1C2H";
}

// compute the sparsity pattern of the chemistry Jacobian
void
SPARSITY_INFO(int* nJdata, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 2304> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 47> conc = {0.0};
  for (int n = 0; n < 47; n++) {
    conc[n] = 1.0 / 47.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 48; k++) {
    for (int l = 0; l < 48; l++) {
      if (Jac[48 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 2304> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 47> conc = {0.0};
  for (int n = 0; n < 47; n++) {
    conc[n] = 1.0 / 47.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 48; k++) {
    for (int l = 0; l < 48; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[48 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 2304> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 47> conc = {0.0};
  for (int n = 0; n < 47; n++) {
    conc[n] = 1.0 / 47.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 48; k++) {
    for (int l = 0; l < 48; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[48 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 2304> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 47> conc = {0.0};
  for (int n = 0; n < 47; n++) {
    conc[n] = 1.0 / 47.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int nc = 0; nc < NCELLS; nc++) {
    int offset_row = nc * 48;
    int offset_col = nc * 48;
    for (int k = 0; k < 48; k++) {
      for (int l = 0; l < 48; l++) {
        if (Jac[48 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 2304> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 47> conc = {0.0};
  for (int n = 0; n < 47; n++) {
    conc[n] = 1.0 / 47.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtrs[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 48;
      for (int l = 0; l < 48; l++) {
        for (int k = 0; k < 48; k++) {
          if (Jac[48 * k + l] != 0.0) {
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
      int offset = nc * 48;
      for (int l = 0; l < 48; l++) {
        for (int k = 0; k < 48; k++) {
          if (Jac[48 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 2304> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 47> conc = {0.0};
  for (int n = 0; n < 47; n++) {
    conc[n] = 1.0 / 47.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 48;
      for (int l = 0; l < 48; l++) {
        for (int k = 0; k < 48; k++) {
          if (k == l) {
            colVals[nJdata_tmp - 1] = l + 1 + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[48 * k + l] != 0.0) {
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
      int offset = nc * 48;
      for (int l = 0; l < 48; l++) {
        for (int k = 0; k < 48; k++) {
          if (k == l) {
            colVals[nJdata_tmp] = l + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[48 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 2304> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 47> conc = {0.0};
  for (int n = 0; n < 47; n++) {
    conc[n] = 1.0 / 47.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int k = 0; k < 48; k++) {
    for (int l = 0; l < 48; l++) {
      if (k == l) {
        rowVals[nJdata_tmp] = l;
        indx[nJdata_tmp] = 48 * k + l;
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[48 * k + l] != 0.0) {
          rowVals[nJdata_tmp] = l;
          indx[nJdata_tmp] = 48 * k + l;
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
  amrex::GpuArray<amrex::Real, 2304> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 47> conc = {0.0};
  for (int n = 0; n < 47; n++) {
    conc[n] = 1.0 / 47.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int l = 0; l < 48; l++) {
      for (int k = 0; k < 48; k++) {
        if (k == l) {
          colVals[nJdata_tmp - 1] = l + 1;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[48 * k + l] != 0.0) {
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
    for (int l = 0; l < 48; l++) {
      for (int k = 0; k < 48; k++) {
        if (k == l) {
          colVals[nJdata_tmp] = l;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[48 * k + l] != 0.0) {
            colVals[nJdata_tmp] = k;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
      }
      rowPtr[l + 1] = nJdata_tmp;
    }
  }
}
