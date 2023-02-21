#include "mechanism.H"
const int rmap[269] = {
  11,  13,  41,  50,  52,  53,  61,  73,  78,  85,  87,  101, 104, 128, 139,
  145, 158, 172, 173, 188, 196, 203, 213, 220, 222, 29,  127, 4,   5,   6,
  7,   8,   9,   10,  38,  39,  245, 252, 0,   1,   2,   3,   12,  14,  15,
  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,  28,  30,  31,
  32,  33,  34,  35,  36,  37,  40,  42,  43,  44,  45,  46,  47,  48,  49,
  51,  54,  55,  56,  57,  58,  59,  60,  62,  63,  64,  65,  66,  67,  68,
  69,  70,  71,  72,  74,  75,  76,  77,  79,  80,  81,  82,  83,  84,  86,
  88,  89,  90,  91,  92,  93,  94,  95,  96,  97,  98,  99,  100, 102, 103,
  105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119,
  120, 121, 122, 123, 124, 125, 126, 129, 130, 131, 132, 133, 134, 135, 136,
  137, 138, 140, 141, 142, 143, 144, 146, 147, 148, 149, 150, 151, 152, 153,
  154, 155, 156, 157, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169,
  170, 171, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186,
  187, 189, 190, 191, 192, 193, 194, 195, 197, 198, 199, 200, 201, 202, 204,
  205, 206, 207, 208, 209, 210, 211, 212, 214, 215, 216, 217, 218, 219, 221,
  223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237,
  238, 239, 240, 241, 242, 243, 244, 246, 247, 248, 249, 250, 251, 253, 254,
  255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267, 268};

// Returns 0-based map of reaction order
void
GET_RMAP(int* _rmap)
{
  for (int j = 0; j < 269; ++j) {
    _rmap[j] = rmap[j];
  }
}

// Returns a count of species in a reaction, and their indices
// and stoichiometric coefficients. (Eq 50)
void
CKINU(const int i, int& nspec, int ki[], int nu[])
{
  const int ns[269] = {
    4, 4, 4, 3, 2, 2, 2, 2, 3, 3, 2, 3, 4, 2, 4, 3, 4, 3, 3, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4,
    3, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 3, 4, 3, 4, 4, 4, 4, 4, 4, 4, 5, 4, 3, 4,
    4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 3, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 2, 3, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 3, 2, 3, 4, 4, 4, 4, 4, 4, 3, 4, 4, 3, 3, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 3, 3, 4, 4, 4, 4, 4, 5, 3, 3, 4, 4, 3, 3, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 3, 4, 4, 4, 4, 4, 5, 4, 3, 4, 4, 4,
    4, 4, 4, 3, 5, 4, 3, 3, 4, 4, 4, 3, 4, 3, 4, 4, 4, 4, 5, 4, 3, 4, 3, 4, 4,
    5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 5,
    4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};
  const int kiv[1345] = {
    2,  9,  3,  4,  0, 6,  3,  2,  4,  0,  6,  4,  2,  7,  0, 4,  7,  3,  0,  0,
    2,  6,  0,  0,  0, 2,  6,  0,  0,  0,  2,  6,  0,  0,  0, 2,  6,  0,  0,  0,
    2,  4,  7,  0,  0, 2,  3,  4,  0,  0,  3,  9,  0,  0,  0, 2,  9,  5,  0,  0,
    6,  9,  2,  5,  0, 4,  8,  0,  0,  0,  2,  5,  7,  3,  0, 2,  5,  4,  0,  0,
    5,  3,  9,  4,  0, 5,  8,  9,  0,  0,  5,  8,  9,  0,  0, 5,  4,  7,  9,  0,
    5,  4,  7,  9,  0, 5,  4,  7,  9,  0,  5,  4,  7,  9,  0, 5,  4,  7,  9,  0,
    2,  8,  6,  5,  0, 2,  8,  7,  4,  0,  8,  3,  5,  4,  0, 8,  4,  7,  5,  0,
    8,  4,  7,  5,  0, 20, 3,  21, 0,  0,  20, 4,  21, 2,  0, 20, 4,  21, 2,  0,
    20, 9,  21, 3,  0, 20, 5,  21, 4,  0,  2,  16, 20, 6,  0, 16, 3,  20, 4,  0,
    16, 3,  21, 2,  0, 16, 4,  20, 7,  0,  16, 20, 2,  0,  0, 16, 20, 2,  0,  0,
    16, 9,  20, 5,  0, 20, 6,  17, 0,  0,  10, 4,  20, 2,  0, 10, 9,  20, 3,  0,
    11, 2,  10, 6,  0, 11, 3,  20, 2,  0,  11, 4,  2,  16, 0, 11, 6,  12, 2,  0,
    11, 7,  17, 2,  0, 11, 9,  16, 3,  0,  11, 20, 28, 0,  0, 11, 21, 20, 16, 0,
    2,  16, 17, 0,  0, 12, 2,  14, 0,  0,  12, 3,  2,  16, 0, 12, 4,  17, 2,  0,
    12, 4,  11, 7,  0, 12, 6,  14, 2,  0,  12, 9,  16, 4,  0, 12, 9,  21, 2,  0,
    12, 5,  17, 4,  0, 12, 20, 29, 0,  0,  11, 12, 22, 2,  0, 12, 22, 6,  0,  0,
    13, 1,  12, 1,  0, 0,  13, 0,  12, 0,  13, 2,  11, 6,  0, 13, 3,  20, 6,  0,
    13, 3,  2,  16, 0, 13, 4,  17, 2,  0,  13, 6,  14, 2,  0, 13, 9,  20, 2,  4,
    13, 9,  20, 7,  0, 13, 7,  19, 0,  0,  13, 7,  12, 7,  0, 13, 20, 12, 20, 0,
    13, 21, 12, 21, 0, 13, 21, 17, 20, 0,  17, 2,  18, 0,  0, 17, 2,  6,  16, 0,
    17, 3,  16, 4,  0, 17, 4,  7,  16, 0,  17, 9,  16, 5,  0, 17, 5,  8,  16, 0,
    11, 17, 29, 2,  0, 14, 2,  15, 0,  0,  14, 3,  17, 2,  0, 14, 4,  19, 0,  0,
    14, 4,  12, 7,  0, 14, 4,  13, 7,  0,  14, 9,  18, 3,  0, 14, 9,  17, 4,  0,
    14, 5,  15, 9,  0, 14, 5,  18, 4,  0,  14, 8,  15, 5,  0, 10, 14, 22, 2,  0,
    11, 14, 24, 2,  0, 14, 16, 15, 20, 0,  17, 14, 15, 16, 0, 12, 14, 25, 2,  0,
    13, 14, 25, 2,  0, 14, 27, 0,  0,  0,  14, 26, 2,  0,  0, 14, 28, 25, 20, 0,
    18, 2,  19, 0,  0, 18, 2,  17, 6,  0,  18, 2,  14, 4,  0, 18, 2,  13, 7,  0,
    18, 3,  17, 4,  0, 18, 4,  17, 7,  0,  18, 9,  17, 5,  0, 15, 2,  14, 6,  0,
    15, 3,  14, 4,  0, 15, 4,  14, 7,  0,  11, 15, 25, 2,  0, 12, 15, 14, 0,  0,
    13, 15, 14, 0,  0, 19, 2,  18, 6,  0,  19, 3,  18, 4,  0, 19, 4,  18, 7,  0,
    14, 19, 18, 15, 0, 2,  28, 13, 20, 0,  28, 3,  20, 2,  0, 28, 9,  20, 4,  0,
    11, 28, 22, 20, 0, 12, 28, 24, 20, 0,  28, 22, 20, 0,  0, 22, 23, 0,  0,  0,
    24, 22, 2,  0,  0, 22, 3,  12, 20, 0,  22, 3,  2,  28, 0, 22, 4,  29, 2,  0,
    22, 4,  14, 20, 0, 22, 16, 24, 20, 0,  22, 14, 2,  32, 0, 22, 14, 33, 0,  0,
    2,  23, 22, 2,  0, 23, 4,  29, 2,  0,  23, 9,  16, 0,  0, 29, 2,  30, 0,  0,
    29, 2,  6,  28, 0, 29, 2,  14, 20, 0,  29, 3,  28, 4,  0, 29, 3,  12, 21, 0,
    29, 4,  7,  28, 0, 24, 2,  25, 0,  0,  24, 2,  22, 6,  0, 24, 2,  6,  23, 0,
    24, 3,  29, 2,  0, 24, 3,  14, 20, 0,  24, 4,  22, 7,  0, 24, 9,  22, 5,  0,
    24, 9,  30, 3,  0, 24, 9,  17, 16, 0,  24, 5,  30, 4,  0, 24, 8,  25, 5,  0,
    24, 16, 25, 20, 0, 24, 14, 22, 15, 0,  24, 14, 34, 0,  0, 24, 14, 2,  33, 0,
    24, 22, 25, 0,  0, 30, 14, 20, 0,  0,  30, 2,  14, 16, 0, 30, 2,  29, 6,  0,
    30, 3,  29, 4,  0, 30, 4,  29, 7,  0,  30, 9,  29, 5,  0, 30, 9,  17, 20, 4,
    31, 14, 16, 0,  0, 31, 15, 20, 0,  0,  31, 2,  24, 7,  0, 31, 2,  25, 4,  0,
    25, 6,  23, 0,  0, 25, 2,  26, 0,  0,  25, 2,  24, 6,  0, 25, 3,  24, 4,  0,
    25, 3,  14, 16, 0, 25, 3,  12, 17, 0,  25, 4,  24, 7,  0, 25, 16, 26, 20, 0,
    25, 11, 2,  32, 0, 25, 12, 2,  33, 0,  25, 13, 15, 23, 0, 25, 13, 2,  33, 0,
    25, 14, 24, 15, 0, 25, 14, 35, 0,  0,  25, 9,  24, 5,  0, 25, 5,  31, 4,  0,
    26, 2,  27, 0,  0, 26, 2,  25, 6,  0,  26, 3,  17, 14, 0, 26, 9,  25, 5,  0,
    26, 5,  27, 9,  0, 26, 5,  25, 8,  0,  26, 5,  17, 14, 4, 26, 8,  27, 5,  0,
    24, 26, 36, 0,  0, 24, 26, 14, 33, 0,  27, 2,  26, 6,  0, 27, 3,  26, 4,  0,
    27, 4,  26, 7,  0, 27, 13, 26, 14, 0,  27, 14, 26, 15, 0, 14, 9,  43, 0,  0,
    43, 5,  17, 7,  9, 43, 17, 19, 9,  0,  43, 18, 9,  0,  0, 14, 43, 18, 0,  0,
    43, 4,  19, 9,  0, 43, 2,  18, 4,  0,  43, 3,  18, 9,  0, 2,  32, 33, 0,  0,
    3,  32, 25, 20, 0, 2,  33, 34, 0,  0,  2,  33, 6,  32, 0, 4,  33, 7,  32, 0,
    9,  33, 5,  32, 0, 5,  33, 34, 9,  0,  5,  33, 24, 17, 4, 16, 33, 34, 20, 0,
    14, 33, 36, 0,  0, 14, 33, 15, 32, 0,  34, 2,  35, 0,  0, 34, 2,  25, 14, 0,
    34, 2,  6,  33, 0, 34, 3,  29, 14, 2,  34, 3,  26, 16, 0, 34, 3,  4,  33, 0,
    34, 4,  7,  33, 0, 34, 5,  8,  33, 0,  34, 14, 15, 33, 0, 2,  35, 26, 14, 0,
    2,  35, 34, 6,  0, 3,  35, 26, 17, 0,  4,  35, 34, 7,  0, 9,  35, 34, 5,  0,
    5,  35, 26, 17, 4, 14, 35, 34, 15, 0,  36, 2,  25, 26, 0, 36, 2,  34, 14, 0,
    36, 3,  16, 35, 0, 37, 38, 1,  3,  0,  37, 9,  38, 3,  0, 37, 4,  2,  38, 0,
    5,  38, 39, 4,  0, 38, 3,  39, 0,  0,  39, 3,  38, 9,  0, 2,  39, 38, 4,  0,
    42, 3,  20, 38, 0, 42, 4,  20, 2,  38, 37, 42, 20, 1,  0, 42, 9,  21, 38, 0,
    42, 20, 37, 0,  0, 42, 38, 21, 1,  0,  10, 38, 20, 37, 0, 11, 38, 2,  42, 0,
    11, 38, 16, 37, 0, 12, 38, 2,  41, 0,  12, 38, 2,  40, 0, 13, 38, 2,  41, 0,
    13, 38, 2,  40, 0, 41, 3,  42, 4,  0,  2,  41, 6,  42, 0, 41, 4,  7,  42, 0,
    2,  40, 2,  41, 0, 28, 38, 20, 40, 0,  21, 37, 20, 38, 0, 14, 39, 18, 38, 0,
    43, 38, 18, 39, 0};
  const int nuv[1345] = {
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -2, 1,  1, 0, 0,
    -2, 1,  0, 0, 0, -2, 1,  0, 0, 0, -2, 1,  0, 0, 0, -2, 1,  0, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -2, 1,  0, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -2, 1,  0, 0, 0, -1, -1, 1, 1, 0, -1, -1, 2, 0, 0,
    -1, -1, 1, 1, 0, -2, 1,  1, 0, 0, -2, 1,  1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, 1,  1, 0, 0, -1, 1,  1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 2, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -2, 1,  1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 1,
    -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -2, 1,  0, 0, 0, -2, 1,  1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 2, 0, 0,
    -1, -1, 2, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 2, 1, 0, -1, -1, 2, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -2, 1,  2, 0, 0, -1, 1,  0, 0, 0,
    -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 2, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0,
    -2, 1,  1, 0, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 1,
    -1, 1,  1, 0, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, 1,  1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 1, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 1, -2, 1,  1, 1, 0, -2, 2,  1, 0, 0, -1, -1, 2, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 1, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 1, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 1, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 1, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0};
  if (i < 1) {
    // Return max num species per reaction
    nspec = 5;
  } else {
    if (i > 269) {
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
  amrex::Real c[44]; // temporary storage
  amrex::Real PORT =
    1e6 * P / (8.31446261815324e+07 * T); // 1e6 * P/RT so c goes to SI units

  // Compute conversion, see Eq 10
  for (int id = 0; id < 44; ++id) {
    c[id] = x[id] * PORT;
  }

  // convert to chemkin units
  progressRateFR(q_f, q_r, c, T);

  // convert to chemkin units
  for (int id = 0; id < 269; ++id) {
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
  amrex::Real g_RT[44];
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
  int kd = 5;
  // Zero ncf
  for (int id = 0; id < kd * 44; ++id) {
    ncf[id] = 0;
  }

  // AR
  ncf[0 * kd + 4] = 1; // Ar

  // N2
  ncf[1 * kd + 3] = 2; // N

  // H
  ncf[2 * kd + 1] = 1; // H

  // O
  ncf[3 * kd + 0] = 1; // O

  // OH
  ncf[4 * kd + 1] = 1; // H
  ncf[4 * kd + 0] = 1; // O

  // HO2
  ncf[5 * kd + 1] = 1; // H
  ncf[5 * kd + 0] = 2; // O

  // H2
  ncf[6 * kd + 1] = 2; // H

  // H2O
  ncf[7 * kd + 1] = 2; // H
  ncf[7 * kd + 0] = 1; // O

  // H2O2
  ncf[8 * kd + 1] = 2; // H
  ncf[8 * kd + 0] = 2; // O

  // O2
  ncf[9 * kd + 0] = 2; // O

  // C
  ncf[10 * kd + 2] = 1; // C

  // CH
  ncf[11 * kd + 2] = 1; // C
  ncf[11 * kd + 1] = 1; // H

  // CH2
  ncf[12 * kd + 2] = 1; // C
  ncf[12 * kd + 1] = 2; // H

  // CH2*
  ncf[13 * kd + 2] = 1; // C
  ncf[13 * kd + 1] = 2; // H

  // CH3
  ncf[14 * kd + 2] = 1; // C
  ncf[14 * kd + 1] = 3; // H

  // CH4
  ncf[15 * kd + 2] = 1; // C
  ncf[15 * kd + 1] = 4; // H

  // HCO
  ncf[16 * kd + 2] = 1; // C
  ncf[16 * kd + 1] = 1; // H
  ncf[16 * kd + 0] = 1; // O

  // CH2O
  ncf[17 * kd + 2] = 1; // C
  ncf[17 * kd + 1] = 2; // H
  ncf[17 * kd + 0] = 1; // O

  // CH3O
  ncf[18 * kd + 2] = 1; // C
  ncf[18 * kd + 1] = 3; // H
  ncf[18 * kd + 0] = 1; // O

  // CH3OH
  ncf[19 * kd + 2] = 1; // C
  ncf[19 * kd + 1] = 4; // H
  ncf[19 * kd + 0] = 1; // O

  // CO
  ncf[20 * kd + 2] = 1; // C
  ncf[20 * kd + 0] = 1; // O

  // CO2
  ncf[21 * kd + 2] = 1; // C
  ncf[21 * kd + 0] = 2; // O

  // C2H2
  ncf[22 * kd + 2] = 2; // C
  ncf[22 * kd + 1] = 2; // H

  // H2CC
  ncf[23 * kd + 2] = 2; // C
  ncf[23 * kd + 1] = 2; // H

  // C2H3
  ncf[24 * kd + 2] = 2; // C
  ncf[24 * kd + 1] = 3; // H

  // C2H4
  ncf[25 * kd + 2] = 2; // C
  ncf[25 * kd + 1] = 4; // H

  // C2H5
  ncf[26 * kd + 2] = 2; // C
  ncf[26 * kd + 1] = 5; // H

  // C2H6
  ncf[27 * kd + 2] = 2; // C
  ncf[27 * kd + 1] = 6; // H

  // HCCO
  ncf[28 * kd + 2] = 2; // C
  ncf[28 * kd + 1] = 1; // H
  ncf[28 * kd + 0] = 1; // O

  // CH2CO
  ncf[29 * kd + 2] = 2; // C
  ncf[29 * kd + 1] = 2; // H
  ncf[29 * kd + 0] = 1; // O

  // CH2CHO
  ncf[30 * kd + 2] = 2; // C
  ncf[30 * kd + 1] = 3; // H
  ncf[30 * kd + 0] = 1; // O

  // CH2OCH2
  ncf[31 * kd + 2] = 2; // C
  ncf[31 * kd + 1] = 4; // H
  ncf[31 * kd + 0] = 1; // O

  // aC3H4
  ncf[32 * kd + 2] = 3; // C
  ncf[32 * kd + 1] = 4; // H

  // aC3H5
  ncf[33 * kd + 2] = 3; // C
  ncf[33 * kd + 1] = 5; // H

  // C3H6
  ncf[34 * kd + 2] = 3; // C
  ncf[34 * kd + 1] = 6; // H

  // nC3H7
  ncf[35 * kd + 2] = 3; // C
  ncf[35 * kd + 1] = 7; // H

  // C4H81
  ncf[36 * kd + 2] = 4; // C
  ncf[36 * kd + 1] = 8; // H

  // N
  ncf[37 * kd + 3] = 1; // N

  // NO
  ncf[38 * kd + 3] = 1; // N
  ncf[38 * kd + 0] = 1; // O

  // NO2
  ncf[39 * kd + 3] = 1; // N
  ncf[39 * kd + 0] = 2; // O

  // HCNO
  ncf[40 * kd + 2] = 1; // C
  ncf[40 * kd + 1] = 1; // H
  ncf[40 * kd + 3] = 1; // N
  ncf[40 * kd + 0] = 1; // O

  // HNCO
  ncf[41 * kd + 2] = 1; // C
  ncf[41 * kd + 1] = 1; // H
  ncf[41 * kd + 3] = 1; // N
  ncf[41 * kd + 0] = 1; // O

  // NCO
  ncf[42 * kd + 2] = 1; // C
  ncf[42 * kd + 3] = 1; // N
  ncf[42 * kd + 0] = 1; // O

  // CH3O2
  ncf[43 * kd + 2] = 1; // C
  ncf[43 * kd + 1] = 3; // H
  ncf[43 * kd + 0] = 2; // O
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
  kname.resize(44);
  kname[0] = "AR";
  kname[1] = "N2";
  kname[2] = "H";
  kname[3] = "O";
  kname[4] = "OH";
  kname[5] = "HO2";
  kname[6] = "H2";
  kname[7] = "H2O";
  kname[8] = "H2O2";
  kname[9] = "O2";
  kname[10] = "C";
  kname[11] = "CH";
  kname[12] = "CH2";
  kname[13] = "CH2*";
  kname[14] = "CH3";
  kname[15] = "CH4";
  kname[16] = "HCO";
  kname[17] = "CH2O";
  kname[18] = "CH3O";
  kname[19] = "CH3OH";
  kname[20] = "CO";
  kname[21] = "CO2";
  kname[22] = "C2H2";
  kname[23] = "H2CC";
  kname[24] = "C2H3";
  kname[25] = "C2H4";
  kname[26] = "C2H5";
  kname[27] = "C2H6";
  kname[28] = "HCCO";
  kname[29] = "CH2CO";
  kname[30] = "CH2CHO";
  kname[31] = "CH2OCH2";
  kname[32] = "aC3H4";
  kname[33] = "aC3H5";
  kname[34] = "C3H6";
  kname[35] = "nC3H7";
  kname[36] = "C4H81";
  kname[37] = "N";
  kname[38] = "NO";
  kname[39] = "NO2";
  kname[40] = "HCNO";
  kname[41] = "HNCO";
  kname[42] = "NCO";
  kname[43] = "CH3O2";
}

// compute the sparsity pattern of the chemistry Jacobian
void
SPARSITY_INFO(int* nJdata, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 2025> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 44> conc = {0.0};
  for (int n = 0; n < 44; n++) {
    conc[n] = 1.0 / 44.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 45; k++) {
    for (int l = 0; l < 45; l++) {
      if (Jac[45 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 2025> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 44> conc = {0.0};
  for (int n = 0; n < 44; n++) {
    conc[n] = 1.0 / 44.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 45; k++) {
    for (int l = 0; l < 45; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[45 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 2025> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 44> conc = {0.0};
  for (int n = 0; n < 44; n++) {
    conc[n] = 1.0 / 44.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 45; k++) {
    for (int l = 0; l < 45; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[45 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 2025> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 44> conc = {0.0};
  for (int n = 0; n < 44; n++) {
    conc[n] = 1.0 / 44.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int nc = 0; nc < NCELLS; nc++) {
    int offset_row = nc * 45;
    int offset_col = nc * 45;
    for (int k = 0; k < 45; k++) {
      for (int l = 0; l < 45; l++) {
        if (Jac[45 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 2025> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 44> conc = {0.0};
  for (int n = 0; n < 44; n++) {
    conc[n] = 1.0 / 44.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtrs[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 45;
      for (int l = 0; l < 45; l++) {
        for (int k = 0; k < 45; k++) {
          if (Jac[45 * k + l] != 0.0) {
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
      int offset = nc * 45;
      for (int l = 0; l < 45; l++) {
        for (int k = 0; k < 45; k++) {
          if (Jac[45 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 2025> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 44> conc = {0.0};
  for (int n = 0; n < 44; n++) {
    conc[n] = 1.0 / 44.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 45;
      for (int l = 0; l < 45; l++) {
        for (int k = 0; k < 45; k++) {
          if (k == l) {
            colVals[nJdata_tmp - 1] = l + 1 + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[45 * k + l] != 0.0) {
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
      int offset = nc * 45;
      for (int l = 0; l < 45; l++) {
        for (int k = 0; k < 45; k++) {
          if (k == l) {
            colVals[nJdata_tmp] = l + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[45 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 2025> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 44> conc = {0.0};
  for (int n = 0; n < 44; n++) {
    conc[n] = 1.0 / 44.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int k = 0; k < 45; k++) {
    for (int l = 0; l < 45; l++) {
      if (k == l) {
        rowVals[nJdata_tmp] = l;
        indx[nJdata_tmp] = 45 * k + l;
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[45 * k + l] != 0.0) {
          rowVals[nJdata_tmp] = l;
          indx[nJdata_tmp] = 45 * k + l;
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
  amrex::GpuArray<amrex::Real, 2025> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 44> conc = {0.0};
  for (int n = 0; n < 44; n++) {
    conc[n] = 1.0 / 44.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int l = 0; l < 45; l++) {
      for (int k = 0; k < 45; k++) {
        if (k == l) {
          colVals[nJdata_tmp - 1] = l + 1;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[45 * k + l] != 0.0) {
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
    for (int l = 0; l < 45; l++) {
      for (int k = 0; k < 45; k++) {
        if (k == l) {
          colVals[nJdata_tmp] = l;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[45 * k + l] != 0.0) {
            colVals[nJdata_tmp] = k;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
      }
      rowPtr[l + 1] = nJdata_tmp;
    }
  }
}
