#include "mechanism.H"
const int rmap[269] = {
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
  195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209,
  210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224,
  225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
  240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254,
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
    3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 2, 2, 2, 2, 2, 3, 3, 2, 3, 3, 3, 3, 4, 4, 4, 3, 4, 4, 3, 4, 3, 3, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 5, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 3, 4,
    4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4,
    4, 4, 5, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4,
    4, 5, 4, 4, 4, 4, 4, 4, 4, 5, 4, 3, 3, 4, 4, 4, 3, 4, 4, 4, 4, 4, 5, 4, 4,
    4, 4, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4};
  const int kiv[1345] = {
    2,  9,  5,  0,  0, 4,  8,  0,  0,  0, 16, 6,  13, 0,  0,  39, 16, 23, 0,  0,
    2,  12, 13, 0,  0, 40, 2,  10, 0,  0, 40, 16, 24, 0,  0,  41, 7,  15, 0,  0,
    13, 2,  14, 0,  0, 10, 2,  11, 0,  0, 10, 4,  15, 0,  0,  10, 22, 0,  0,  0,
    14, 2,  15, 0,  0, 19, 18, 2,  0,  0, 24, 2,  25, 0,  0,  19, 2,  20, 0,  0,
    19, 10, 29, 0,  0, 20, 6,  42, 0,  0, 20, 2,  21, 0,  0,  21, 2,  22, 0,  0,
    19, 21, 31, 0,  0, 10, 9,  37, 0,  0, 2,  28, 29, 0,  0,  10, 28, 31, 0,  0,
    29, 2,  30, 0,  0, 16, 3,  17, 0,  0, 18, 42, 0,  0,  0,  2,  6,  0,  0,  0,
    2,  6,  0,  0,  0, 2,  6,  0,  0,  0, 2,  6,  0,  0,  0,  2,  4,  7,  0,  0,
    2,  3,  4,  0,  0, 3,  9,  0,  0,  0, 12, 16, 2,  0,  0,  12, 16, 2,  0,  0,
    43, 3,  33, 0,  0, 36, 16, 32, 0,  0, 2,  9,  3,  4,  0,  6,  3,  2,  4,  0,
    6,  4,  2,  7,  0, 4,  7,  3,  0,  0, 6,  9,  2,  5,  0,  2,  5,  7,  3,  0,
    2,  5,  4,  0,  0, 5,  3,  9,  4,  0, 5,  8,  9,  0,  0,  5,  8,  9,  0,  0,
    5,  4,  7,  9,  0, 5,  4,  7,  9,  0, 5,  4,  7,  9,  0,  5,  4,  7,  9,  0,
    5,  4,  7,  9,  0, 2,  8,  6,  5,  0, 2,  8,  7,  4,  0,  8,  3,  5,  4,  0,
    8,  4,  7,  5,  0, 8,  4,  7,  5,  0, 16, 4,  17, 2,  0,  16, 4,  17, 2,  0,
    16, 9,  17, 3,  0, 16, 5,  17, 4,  0, 2,  12, 16, 6,  0,  12, 3,  16, 4,  0,
    12, 3,  17, 2,  0, 12, 4,  16, 7,  0, 12, 9,  16, 5,  0,  38, 4,  16, 2,  0,
    38, 9,  16, 3,  0, 39, 2,  38, 6,  0, 39, 3,  16, 2,  0,  39, 4,  2,  12, 0,
    39, 6,  40, 2,  0, 39, 7,  13, 2,  0, 39, 9,  12, 3,  0,  39, 17, 16, 12, 0,
    40, 3,  2,  12, 0, 40, 4,  13, 2,  0, 40, 4,  39, 7,  0,  40, 6,  10, 2,  0,
    40, 9,  12, 4,  0, 40, 9,  17, 2,  0, 40, 5,  13, 4,  0,  39, 40, 18, 2,  0,
    40, 18, 6,  0,  0, 41, 1,  40, 1,  0, 0,  41, 0,  40, 0,  41, 2,  39, 6,  0,
    41, 3,  16, 6,  0, 41, 3,  2,  12, 0, 41, 4,  13, 2,  0,  41, 6,  10, 2,  0,
    41, 9,  16, 2,  4, 41, 9,  16, 7,  0, 41, 7,  40, 7,  0,  41, 16, 40, 16, 0,
    41, 17, 40, 17, 0, 41, 17, 13, 16, 0, 13, 2,  6,  12, 0,  13, 3,  12, 4,  0,
    13, 4,  7,  12, 0, 13, 9,  12, 5,  0, 13, 5,  8,  12, 0,  39, 13, 24, 2,  0,
    10, 3,  13, 2,  0, 10, 4,  40, 7,  0, 10, 4,  41, 7,  0,  10, 9,  14, 3,  0,
    10, 9,  13, 4,  0, 10, 5,  11, 9,  0, 10, 5,  14, 4,  0,  10, 8,  11, 5,  0,
    38, 10, 18, 2,  0, 39, 10, 19, 2,  0, 10, 12, 11, 16, 0,  13, 10, 11, 12, 0,
    40, 10, 20, 2,  0, 41, 10, 20, 2,  0, 10, 21, 2,  0,  0,  10, 23, 20, 16, 0,
    14, 2,  13, 6,  0, 14, 2,  10, 4,  0, 14, 2,  41, 7,  0,  14, 3,  13, 4,  0,
    14, 4,  13, 7,  0, 14, 9,  13, 5,  0, 11, 2,  10, 6,  0,  11, 3,  10, 4,  0,
    11, 4,  10, 7,  0, 39, 11, 20, 2,  0, 40, 11, 10, 0,  0,  41, 11, 10, 0,  0,
    15, 2,  14, 6,  0, 15, 3,  14, 4,  0, 15, 4,  14, 7,  0,  10, 15, 14, 11, 0,
    2,  23, 41, 16, 0, 23, 3,  16, 2,  0, 23, 9,  16, 4,  0,  39, 23, 18, 16, 0,
    40, 23, 19, 16, 0, 23, 18, 16, 0,  0, 18, 3,  40, 16, 0,  18, 3,  2,  23, 0,
    18, 4,  24, 2,  0, 18, 4,  10, 16, 0, 18, 12, 19, 16, 0,  18, 10, 2,  27, 0,
    18, 10, 28, 0,  0, 2,  42, 18, 2,  0, 42, 4,  24, 2,  0,  42, 9,  12, 0,  0,
    24, 2,  6,  23, 0, 24, 2,  10, 16, 0, 24, 3,  23, 4,  0,  24, 3,  40, 17, 0,
    24, 4,  7,  23, 0, 19, 2,  18, 6,  0, 19, 2,  6,  42, 0,  19, 3,  24, 2,  0,
    19, 3,  10, 16, 0, 19, 4,  18, 7,  0, 19, 9,  18, 5,  0,  19, 9,  25, 3,  0,
    19, 9,  13, 12, 0, 19, 5,  25, 4,  0, 19, 8,  20, 5,  0,  19, 12, 20, 16, 0,
    19, 10, 18, 11, 0, 19, 10, 2,  28, 0, 19, 18, 20, 0,  0,  25, 10, 16, 0,  0,
    25, 2,  10, 12, 0, 25, 2,  24, 6,  0, 25, 3,  24, 4,  0,  25, 4,  24, 7,  0,
    25, 9,  24, 5,  0, 25, 9,  13, 16, 4, 26, 10, 12, 0,  0,  26, 11, 16, 0,  0,
    26, 2,  19, 7,  0, 26, 2,  20, 4,  0, 20, 2,  19, 6,  0,  20, 3,  19, 4,  0,
    20, 3,  10, 12, 0, 20, 3,  40, 13, 0, 20, 4,  19, 7,  0,  20, 12, 21, 16, 0,
    20, 39, 2,  27, 0, 20, 40, 2,  28, 0, 20, 41, 11, 42, 0,  20, 41, 2,  28, 0,
    20, 10, 19, 11, 0, 20, 10, 30, 0,  0, 20, 9,  19, 5,  0,  20, 5,  26, 4,  0,
    21, 2,  20, 6,  0, 21, 3,  13, 10, 0, 21, 9,  20, 5,  0,  21, 5,  22, 9,  0,
    21, 5,  20, 8,  0, 21, 5,  13, 10, 4, 21, 8,  22, 5,  0,  19, 21, 10, 28, 0,
    22, 2,  21, 6,  0, 22, 3,  21, 4,  0, 22, 4,  21, 7,  0,  22, 41, 21, 10, 0,
    22, 10, 21, 11, 0, 37, 5,  13, 7,  9, 37, 13, 15, 9,  0,  37, 14, 9,  0,  0,
    10, 37, 14, 0,  0, 37, 4,  15, 9,  0, 37, 2,  14, 4,  0,  37, 3,  14, 9,  0,
    2,  27, 28, 0,  0, 3,  27, 20, 16, 0, 2,  28, 6,  27, 0,  4,  28, 7,  27, 0,
    9,  28, 5,  27, 0, 5,  28, 29, 9,  0, 5,  28, 19, 13, 4,  12, 28, 29, 16, 0,
    10, 28, 11, 27, 0, 29, 2,  20, 10, 0, 29, 2,  6,  28, 0,  29, 3,  24, 10, 2,
    29, 3,  21, 12, 0, 29, 3,  4,  28, 0, 29, 4,  7,  28, 0,  29, 5,  8,  28, 0,
    29, 10, 11, 28, 0, 2,  30, 21, 10, 0, 2,  30, 29, 6,  0,  3,  30, 21, 13, 0,
    4,  30, 29, 7,  0, 9,  30, 29, 5,  0, 5,  30, 21, 13, 4,  10, 30, 29, 11, 0,
    31, 2,  20, 21, 0, 31, 2,  29, 10, 0, 31, 3,  12, 30, 0,  32, 43, 1,  3,  0,
    32, 9,  43, 3,  0, 32, 4,  2,  43, 0, 5,  43, 33, 4,  0,  33, 3,  43, 9,  0,
    2,  33, 43, 4,  0, 36, 3,  16, 43, 0, 36, 4,  16, 2,  43, 32, 36, 16, 1,  0,
    36, 9,  17, 43, 0, 36, 43, 17, 1,  0, 38, 43, 16, 32, 0,  39, 43, 2,  36, 0,
    39, 43, 12, 32, 0, 40, 43, 2,  35, 0, 40, 43, 2,  34, 0,  41, 43, 2,  35, 0,
    41, 43, 2,  34, 0, 35, 3,  36, 4,  0, 2,  35, 6,  36, 0,  35, 4,  7,  36, 0,
    2,  34, 2,  35, 0, 23, 43, 16, 34, 0, 17, 32, 16, 43, 0,  10, 33, 14, 43, 0,
    37, 43, 14, 33, 0};
  const int nuv[1345] = {
    -1, -1, 1, 0, 0, -2, 1,  0, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -2, 1,  0, 0, 0,
    -1, -1, 1, 0, 0, -1, 1,  1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 0, 0, -1, 1,  1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, 1,  0, 0, 0, -2, 1,  0, 0, 0,
    -2, 1,  0, 0, 0, -2, 1,  0, 0, 0, -2, 1,  0, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 0, 0, -2, 1,  0, 0, 0, -1, 1,  1, 0, 0, -1, 1,  1, 0, 0,
    -1, -1, 1, 0, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -2, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 2, 0, 0, -1, -1, 1, 1, 0, -2, 1,  1, 0, 0, -2, 1,  1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 2, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -2, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 1, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -2, 1,  1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 2, 0, 0, -1, -1, 2, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 2, 1, 0, -1, -1, 2, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -2, 1,  2, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 2, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -2, 1,  1, 0, 0, -1, 1,  1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 1, -1, 1,  1, 0, 0, -1, 1,  1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 1, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 1, -2, 1,  1, 1, 0, -2, 2,  1, 0, 0,
    -1, -1, 2, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 1, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 1,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 1, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 1, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
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
  int id;            // loop counter
  amrex::Real c[38]; // temporary storage
  amrex::Real PORT =
    1e6 * P / (8.31446261815324e+07 * T); // 1e6 * P/RT so c goes to SI units

  // Compute conversion, see Eq 10
  for (id = 0; id < 38; ++id) {
    c[id] = x[id] * PORT;
  }

  // convert to chemkin units
  progressRateFR(q_f, q_r, c, T);

  // convert to chemkin units
  for (id = 0; id < 269; ++id) {
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
  amrex::Real g_RT[38];
  gibbs(g_RT, tc);
  amrex::Real g_RT_qss[6];
  gibbs_qss(g_RT_qss, tc);

  amrex::Real sc_qss[6];
  // Fill sc_qss here
  amrex::Real kf_qss[85], qf_qss[85], qr_qss[85];
  comp_k_f_qss(tc, invT, kf_qss);
  comp_qss_coeff(kf_qss, qf_qss, qr_qss, sc, tc, g_RT, g_RT_qss);
  comp_sc_qss(sc_qss, qf_qss, qr_qss);
  comp_qfqr(q_f, q_r, sc, sc_qss, tc, invT);
}

// save atomic weights into array
void
atomicWeight(amrex::Real* awt)
{
  awt[0] = 39.950000; // Ar
  awt[1] = 14.007000; // N
  awt[2] = 1.008000;  // H
  awt[3] = 15.999000; // O
  awt[4] = 12.011000; // C
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
  for (id = 0; id < kd * 38; ++id) {
    ncf[id] = 0;
  }

  // AR
  ncf[0 * kd + 0] = 1; // Ar

  // N2
  ncf[1 * kd + 1] = 2; // N

  // H
  ncf[2 * kd + 2] = 1; // H

  // O
  ncf[3 * kd + 3] = 1; // O

  // OH
  ncf[4 * kd + 2] = 1; // H
  ncf[4 * kd + 3] = 1; // O

  // HO2
  ncf[5 * kd + 2] = 1; // H
  ncf[5 * kd + 3] = 2; // O

  // H2
  ncf[6 * kd + 2] = 2; // H

  // H2O
  ncf[7 * kd + 2] = 2; // H
  ncf[7 * kd + 3] = 1; // O

  // H2O2
  ncf[8 * kd + 2] = 2; // H
  ncf[8 * kd + 3] = 2; // O

  // O2
  ncf[9 * kd + 3] = 2; // O

  // CH3
  ncf[10 * kd + 4] = 1; // C
  ncf[10 * kd + 2] = 3; // H

  // CH4
  ncf[11 * kd + 4] = 1; // C
  ncf[11 * kd + 2] = 4; // H

  // HCO
  ncf[12 * kd + 4] = 1; // C
  ncf[12 * kd + 2] = 1; // H
  ncf[12 * kd + 3] = 1; // O

  // CH2O
  ncf[13 * kd + 4] = 1; // C
  ncf[13 * kd + 2] = 2; // H
  ncf[13 * kd + 3] = 1; // O

  // CH3O
  ncf[14 * kd + 4] = 1; // C
  ncf[14 * kd + 2] = 3; // H
  ncf[14 * kd + 3] = 1; // O

  // CH3OH
  ncf[15 * kd + 4] = 1; // C
  ncf[15 * kd + 2] = 4; // H
  ncf[15 * kd + 3] = 1; // O

  // CO
  ncf[16 * kd + 4] = 1; // C
  ncf[16 * kd + 3] = 1; // O

  // CO2
  ncf[17 * kd + 4] = 1; // C
  ncf[17 * kd + 3] = 2; // O

  // C2H2
  ncf[18 * kd + 4] = 2; // C
  ncf[18 * kd + 2] = 2; // H

  // C2H3
  ncf[19 * kd + 4] = 2; // C
  ncf[19 * kd + 2] = 3; // H

  // C2H4
  ncf[20 * kd + 4] = 2; // C
  ncf[20 * kd + 2] = 4; // H

  // C2H5
  ncf[21 * kd + 4] = 2; // C
  ncf[21 * kd + 2] = 5; // H

  // C2H6
  ncf[22 * kd + 4] = 2; // C
  ncf[22 * kd + 2] = 6; // H

  // HCCO
  ncf[23 * kd + 4] = 2; // C
  ncf[23 * kd + 2] = 1; // H
  ncf[23 * kd + 3] = 1; // O

  // CH2CO
  ncf[24 * kd + 4] = 2; // C
  ncf[24 * kd + 2] = 2; // H
  ncf[24 * kd + 3] = 1; // O

  // CH2CHO
  ncf[25 * kd + 4] = 2; // C
  ncf[25 * kd + 2] = 3; // H
  ncf[25 * kd + 3] = 1; // O

  // CH2OCH2
  ncf[26 * kd + 4] = 2; // C
  ncf[26 * kd + 2] = 4; // H
  ncf[26 * kd + 3] = 1; // O

  // aC3H4
  ncf[27 * kd + 4] = 3; // C
  ncf[27 * kd + 2] = 4; // H

  // aC3H5
  ncf[28 * kd + 4] = 3; // C
  ncf[28 * kd + 2] = 5; // H

  // C3H6
  ncf[29 * kd + 4] = 3; // C
  ncf[29 * kd + 2] = 6; // H

  // nC3H7
  ncf[30 * kd + 4] = 3; // C
  ncf[30 * kd + 2] = 7; // H

  // C4H81
  ncf[31 * kd + 4] = 4; // C
  ncf[31 * kd + 2] = 8; // H

  // N
  ncf[32 * kd + 1] = 1; // N

  // NO2
  ncf[33 * kd + 1] = 1; // N
  ncf[33 * kd + 3] = 2; // O

  // HCNO
  ncf[34 * kd + 4] = 1; // C
  ncf[34 * kd + 2] = 1; // H
  ncf[34 * kd + 1] = 1; // N
  ncf[34 * kd + 3] = 1; // O

  // HNCO
  ncf[35 * kd + 4] = 1; // C
  ncf[35 * kd + 2] = 1; // H
  ncf[35 * kd + 1] = 1; // N
  ncf[35 * kd + 3] = 1; // O

  // NCO
  ncf[36 * kd + 4] = 1; // C
  ncf[36 * kd + 1] = 1; // N
  ncf[36 * kd + 3] = 1; // O

  // CH3O2
  ncf[37 * kd + 4] = 1; // C
  ncf[37 * kd + 2] = 3; // H
  ncf[37 * kd + 3] = 2; // O
}

// Returns the vector of strings of element names
void
CKSYME_STR(amrex::Vector<std::string>& ename)
{
  ename.resize(5);
  ename[0] = "Ar";
  ename[1] = "N";
  ename[2] = "H";
  ename[3] = "O";
  ename[4] = "C";
}

// Returns the vector of strings of species names
void
CKSYMS_STR(amrex::Vector<std::string>& kname)
{
  kname.resize(38);
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
  kname[10] = "CH3";
  kname[11] = "CH4";
  kname[12] = "HCO";
  kname[13] = "CH2O";
  kname[14] = "CH3O";
  kname[15] = "CH3OH";
  kname[16] = "CO";
  kname[17] = "CO2";
  kname[18] = "C2H2";
  kname[19] = "C2H3";
  kname[20] = "C2H4";
  kname[21] = "C2H5";
  kname[22] = "C2H6";
  kname[23] = "HCCO";
  kname[24] = "CH2CO";
  kname[25] = "CH2CHO";
  kname[26] = "CH2OCH2";
  kname[27] = "aC3H4";
  kname[28] = "aC3H5";
  kname[29] = "C3H6";
  kname[30] = "nC3H7";
  kname[31] = "C4H81";
  kname[32] = "N";
  kname[33] = "NO2";
  kname[34] = "HCNO";
  kname[35] = "HNCO";
  kname[36] = "NCO";
  kname[37] = "CH3O2";
}

// compute the sparsity pattern of the chemistry Jacobian
void
SPARSITY_INFO(int* nJdata, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 1521> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 38> conc = {0.0};
  for (int n = 0; n < 38; n++) {
    conc[n] = 1.0 / 38.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 39; k++) {
    for (int l = 0; l < 39; l++) {
      if (Jac[39 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 1521> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 38> conc = {0.0};
  for (int n = 0; n < 38; n++) {
    conc[n] = 1.0 / 38.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 39; k++) {
    for (int l = 0; l < 39; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[39 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 1521> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 38> conc = {0.0};
  for (int n = 0; n < 38; n++) {
    conc[n] = 1.0 / 38.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 39; k++) {
    for (int l = 0; l < 39; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[39 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 1521> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 38> conc = {0.0};
  for (int n = 0; n < 38; n++) {
    conc[n] = 1.0 / 38.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int nc = 0; nc < NCELLS; nc++) {
    int offset_row = nc * 39;
    int offset_col = nc * 39;
    for (int k = 0; k < 39; k++) {
      for (int l = 0; l < 39; l++) {
        if (Jac[39 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 1521> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 38> conc = {0.0};
  for (int n = 0; n < 38; n++) {
    conc[n] = 1.0 / 38.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtrs[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 39;
      for (int l = 0; l < 39; l++) {
        for (int k = 0; k < 39; k++) {
          if (Jac[39 * k + l] != 0.0) {
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
      int offset = nc * 39;
      for (int l = 0; l < 39; l++) {
        for (int k = 0; k < 39; k++) {
          if (Jac[39 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 1521> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 38> conc = {0.0};
  for (int n = 0; n < 38; n++) {
    conc[n] = 1.0 / 38.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 39;
      for (int l = 0; l < 39; l++) {
        for (int k = 0; k < 39; k++) {
          if (k == l) {
            colVals[nJdata_tmp - 1] = l + 1 + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[39 * k + l] != 0.0) {
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
      int offset = nc * 39;
      for (int l = 0; l < 39; l++) {
        for (int k = 0; k < 39; k++) {
          if (k == l) {
            colVals[nJdata_tmp] = l + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[39 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 1521> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 38> conc = {0.0};
  for (int n = 0; n < 38; n++) {
    conc[n] = 1.0 / 38.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int k = 0; k < 39; k++) {
    for (int l = 0; l < 39; l++) {
      if (k == l) {
        rowVals[nJdata_tmp] = l;
        indx[nJdata_tmp] = 39 * k + l;
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[39 * k + l] != 0.0) {
          rowVals[nJdata_tmp] = l;
          indx[nJdata_tmp] = 39 * k + l;
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
  amrex::GpuArray<amrex::Real, 1521> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 38> conc = {0.0};
  for (int n = 0; n < 38; n++) {
    conc[n] = 1.0 / 38.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int l = 0; l < 39; l++) {
      for (int k = 0; k < 39; k++) {
        if (k == l) {
          colVals[nJdata_tmp - 1] = l + 1;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[39 * k + l] != 0.0) {
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
    for (int l = 0; l < 39; l++) {
      for (int k = 0; k < 39; k++) {
        if (k == l) {
          colVals[nJdata_tmp] = l;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[39 * k + l] != 0.0) {
            colVals[nJdata_tmp] = k;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
      }
      rowPtr[l + 1] = nJdata_tmp;
    }
  }
}
