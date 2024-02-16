#include "mechanism.H"
const int rmap[NUM_REACTIONS] = {
  7,   17,  41,  43,  46,  55,  56,  72,  96,  100, 125, 4,   29,  30,  64,
  77,  80,  86,  87,  126, 0,   1,   2,   3,   5,   6,   8,   9,   10,  11,
  12,  13,  14,  15,  16,  18,  19,  20,  21,  22,  23,  24,  25,  26,  27,
  28,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  42,  44,  45,  47,
  48,  49,  50,  51,  52,  53,  54,  57,  58,  59,  60,  61,  62,  63,  65,
  66,  67,  68,  69,  70,  71,  73,  74,  75,  76,  78,  79,  81,  82,  83,
  84,  85,  88,  89,  90,  91,  92,  93,  94,  95,  97,  98,  99,  101, 102,
  103, 104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116, 117,
  118, 119, 120, 121, 122, 123, 124, 127, 128, 129, 130, 131, 132, 133, 134,
  135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149,
  150, 151, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164,
  165, 166, 167, 168, 169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179,
  180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194,
  195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208, 209,
  210, 211, 212, 213, 214, 215, 216, 217};

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
    4, 4, 4, 4, 3, 3, 3, 3, 4, 4, 4, 4, 3, 4, 3, 3, 4, 2, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 2, 2, 4, 4, 5, 4, 4, 4, 4, 4, 4, 4, 3, 4, 3, 4, 4, 2, 3, 4, 4,
    4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 3, 4, 3, 4, 4,
    4, 5, 3, 3, 3, 3, 5, 3, 3, 4, 4, 3, 3, 3, 4, 4, 4, 4, 3, 4, 4, 3, 4, 4, 4,
    3, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 3, 3,
    3, 3, 3, 3, 3, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 4, 3, 4, 4,
    4, 4, 5, 4, 4, 3, 3, 3, 3, 3, 3, 2, 2, 4, 4, 4, 4, 4, 3, 4, 4, 3, 3, 4, 4,
    4, 3, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 3, 3, 3, 4, 4, 3, 3, 4, 3, 3,
    3, 3, 4, 3, 3, 3, 3, 3, 4, 4, 4, 3, 4, 3, 4, 4, 4, 4};
  const int kiv[NUM_GAS_REACTIONS * 5] = {
    2,  1,  3,  4,  0, 3,  4,  2,  1,  0, 2,  4,  3,  5,  0,  3,  5,  2,  4,  0,
    3,  4,  5,  0,  0, 5,  1,  4,  0,  0, 4,  5,  1,  0,  0,  3,  6,  7,  0,  0,
    3,  6,  1,  4,  0, 1,  4,  3,  6,  0, 7,  4,  5,  6,  0,  7,  1,  6,  4,  0,
    7,  8,  6,  0,  0, 3,  7,  2,  6,  0, 7,  8,  6,  0,  0,  3,  7,  4,  0,  0,
    8,  4,  5,  7,  0, 8,  4,  0,  0,  0, 9,  6,  10, 1,  0,  11, 6,  12, 3,  0,
    11, 6,  13, 5,  0, 11, 1,  13, 3,  0, 11, 6,  14, 1,  0,  11, 3,  9,  2,  0,
    9,  2,  11, 3,  0, 11, 4,  9,  5,  0, 9,  5,  11, 4,  0,  11, 6,  12, 2,  0,
    15, 3,  9,  2,  0, 15, 11, 0,  0,  0, 11, 15, 0,  0,  0,  15, 2,  16, 3,  0,
    16, 3,  15, 2,  0, 15, 6,  13, 3,  4, 15, 4,  14, 3,  0,  16, 4,  14, 2,  0,
    16, 4,  15, 5,  0, 15, 5,  16, 4,  0, 16, 1,  14, 3,  0,  16, 7,  17, 4,  0,
    16, 7,  18, 6,  0, 16, 4,  19, 0,  0, 16, 6,  14, 4,  0,  16, 3,  18, 0,  0,
    16, 3,  11, 2,  0, 11, 2,  16, 3,  0, 16, 20, 0,  0,  0,  16, 21, 3,  0,  0,
    16, 4,  11, 5,  0, 11, 5,  16, 4,  0, 18, 1,  16, 4,  0,  18, 3,  16, 2,  0,
    16, 2,  18, 3,  0, 18, 4,  16, 5,  0, 16, 5,  18, 4,  0,  11, 13, 22, 0,  0,
    13, 1,  12, 0,  0, 13, 4,  12, 3,  0, 12, 3,  13, 4,  0,  10, 6,  13, 7,  0,
    10, 1,  12, 3,  0, 10, 4,  13, 5,  0, 3,  10, 13, 2,  0,  10, 1,  13, 4,  0,
    10, 13, 3,  0,  0, 16, 10, 18, 13, 0, 14, 4,  5,  10, 0,  14, 1,  10, 4,  0,
    14, 3,  2,  10, 0, 14, 16, 18, 10, 0, 17, 14, 19, 0,  0,  17, 6,  14, 7,  0,
    17, 14, 3,  0,  0, 17, 2,  19, 3,  0, 19, 4,  17, 5,  0,  15, 12, 14, 13, 0,
    3,  23, 13, 2,  4, 23, 13, 5,  0,  0, 23, 10, 4,  0,  0,  10, 4,  23, 0,  0,
    23, 12, 2,  0,  0, 23, 4,  12, 3,  5, 24, 17, 6,  0,  0,  16, 24, 17, 0,  0,
    24, 14, 19, 6,  0, 24, 7,  25, 6,  0, 24, 16, 6,  0,  0,  16, 6,  24, 0,  0,
    25, 17, 4,  0,  0, 26, 1,  11, 13, 0, 26, 1,  3,  27, 0,  28, 3,  26, 2,  0,
    28, 6,  29, 1,  0, 28, 16, 30, 0,  0, 28, 6,  26, 7,  0,  28, 6,  14, 10, 0,
    28, 26, 3,  0,  0, 31, 16, 28, 18, 0, 31, 1,  16, 10, 0,  31, 4,  28, 5,  0,
    31, 3,  21, 0,  0, 31, 1,  29, 3,  0, 31, 3,  28, 2,  0,  28, 2,  31, 3,  0,
    21, 3,  20, 0,  0, 21, 24, 32, 17, 0, 21, 7,  32, 4,  0,  21, 6,  31, 7,  0,
    20, 1,  21, 4,  0, 20, 4,  21, 5,  0, 20, 3,  21, 2,  0,  27, 1,  13, 3,  0,
    27, 4,  10, 0,  0, 27, 6,  12, 10, 0, 3,  27, 15, 13, 0,  15, 13, 3,  27, 0,
    22, 1,  27, 4,  0, 22, 3,  2,  27, 0, 2,  27, 22, 3,  0,  22, 3,  16, 13, 0,
    22, 1,  11, 12, 0, 22, 4,  5,  27, 0, 29, 6,  14, 13, 4,  29, 22, 3,  0,  0,
    22, 3,  29, 0,  0, 33, 16, 13, 0,  0, 32, 14, 16, 0,  0,  34, 21, 6,  0,  0,
    21, 6,  34, 0,  0, 34, 31, 7,  0,  0, 35, 6,  13, 3,  27, 35, 4,  26, 10, 0,
    36, 6,  22, 10, 0, 36, 7,  37, 6,  0, 36, 3,  35, 2,  0,  36, 4,  35, 5,  0,
    35, 5,  36, 4,  0, 37, 3,  36, 2,  0, 37, 4,  36, 5,  0,  37, 1,  31, 13, 0,
    38, 3,  37, 2,  0, 38, 7,  30, 6,  0, 38, 3,  30, 0,  0,  38, 26, 16, 0,  0,
    38, 37, 3,  0,  0, 37, 3,  38, 0,  0, 38, 14, 30, 10, 0,  38, 37, 30, 0,  0,
    30, 3,  31, 16, 0, 30, 3,  38, 2,  0, 30, 1,  21, 10, 0,  30, 1,  38, 4,  0,
    30, 1,  22, 16, 3, 30, 4,  38, 5,  0, 39, 6,  30, 7,  0,  39, 31, 16, 0,  0,
    31, 16, 39, 0,  0, 39, 30, 3,  0,  0, 30, 3,  39, 0,  0,  40, 39, 6,  0,  0,
    39, 6,  40, 0,  0, 41, 28, 0,  0,  0, 28, 41, 0,  0,  0,  41, 4,  38, 14, 0,
    41, 4,  21, 22, 0, 41, 1,  31, 22, 0, 41, 3,  28, 31, 0,  41, 1,  37, 14, 0,
    42, 3,  43, 0,  0, 38, 42, 30, 41, 0, 21, 42, 20, 41, 0,  42, 41, 3,  0,  0,
    41, 3,  42, 0,  0, 42, 16, 41, 18, 0, 42, 7,  43, 6,  0,  42, 6,  41, 7,  0,
    42, 28, 31, 0,  0, 42, 3,  41, 2,  0, 43, 3,  42, 2,  0,  43, 4,  14, 39, 0,
    43, 4,  20, 33, 0, 43, 1,  21, 33, 0, 43, 1,  30, 14, 0,  43, 4,  42, 5,  0,
    43, 38, 16, 0,  0, 38, 16, 43, 0,  0, 44, 43, 3,  0,  0,  43, 3,  44, 0,  0,
    44, 31, 21, 0,  0, 45, 6,  44, 0,  0, 6,  44, 45, 0,  0,  46, 41, 16, 0,  0,
    46, 31, 38, 0,  0, 47, 4,  46, 5,  0, 47, 3,  46, 2,  0,  47, 21, 38, 0,  0,
    21, 38, 47, 0,  0, 47, 1,  46, 4,  0, 48, 21, 30, 0,  0,  48, 31, 39, 0,  0,
    48, 47, 3,  0,  0, 49, 38, 39, 0,  0, 49, 4,  48, 14, 0,  50, 49, 16, 0,  0,
    50, 30, 44, 0,  0, 50, 43, 39, 0,  0, 50, 31, 48, 0,  0,  50, 21, 47, 0,  0,
    50, 7,  51, 6,  0, 24, 51, 50, 25, 0, 3,  51, 50, 2,  0,  51, 39, 44, 0,  0,
    7,  51, 50, 8,  0, 51, 21, 48, 0,  0, 17, 51, 50, 19, 0,  51, 1,  50, 4,  0,
    51, 4,  50, 5,  0, 16, 51, 50, 18, 0};
  const int nuv[NUM_GAS_REACTIONS * 5] = {
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 2, 0, 0, -2, 1,  1, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -2, 1,  1, 0, 0, -1, -1, 1, 1, 0, -2, 1,  1, 0, 0, -1, -1, 2, 0, 0,
    -1, -1, 1, 1, 0, -1, 2,  0, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 2, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 2, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, 1,  0, 0, 0, -1, 1,  0, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 1, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -2, 1,  0, 0, 0, -2, 1,  1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -2, 1,  1, 0, 0, -1, -1, 1, 1, 0,
    -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 1, -1, 1,  1, 0, 0, -1, 1,  1, 0, 0, -1, -1, 1, 0, 0,
    -1, 1,  1, 0, 0, -1, -1, 1, 1, 1, -2, 2,  1, 0, 0, -1, -1, 2, 0, 0,
    -2, 1,  1, 1, 0, -1, -1, 1, 1, 0, -1, 1,  1, 0, 0, -1, -1, 1, 0, 0,
    -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 2, 1, 0,
    -1, -1, 2, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 1, -1, 1,  1, 0, 0,
    -1, -1, 1, 0, 0, -1, 1,  1, 0, 0, -1, 1,  1, 0, 0, -1, 1,  1, 0, 0,
    -1, -1, 1, 0, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 1, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, 1,  1, 0, 0,
    -1, 1,  1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -2, 1,  1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 1, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, 1,  1, 0, 0,
    -1, -1, 1, 0, 0, -1, 1,  1, 0, 0, -1, -1, 1, 0, 0, -1, 1,  1, 0, 0,
    -1, -1, 1, 0, 0, -1, 2,  0, 0, 0, -2, 1,  0, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, 1,  1, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, 1,  1, 0, 0, -1, -1, 1, 0, 0, -1, 1,  1, 0, 0, -1, -1, 1, 0, 0,
    -1, 1,  1, 0, 0, -1, 1,  1, 0, 0, -1, -1, 1, 0, 0, -1, 1,  1, 0, 0,
    -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, 1,  1, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, 1,  1, 0, 0, -1, 1,  1, 0, 0,
    -1, 1,  1, 0, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, 1,  1, 0, 0,
    -1, 1,  1, 0, 0, -1, 1,  1, 0, 0, -1, 1,  1, 0, 0, -1, 1,  1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, 1,  1, 0, 0,
    -1, -1, 1, 1, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0};
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
  amrex::Real c[52]; // temporary storage
  amrex::Real PORT =
    1e6 * P / (8.31446261815324e+07 * T); // 1e6 * P/RT so c goes to SI units

  // Compute conversion, see Eq 10
  for (int id = 0; id < 52; ++id) {
    c[id] = x[id] * PORT;
  }

  // convert to chemkin units
  progressRateFR(q_f, q_r, c, T);

  // convert to chemkin units
  for (int id = 0; id < 218; ++id) {
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
  amrex::Real g_RT[52];
  gibbs(g_RT, T);

  amrex::Real sc_qss[1];
  comp_qfqr(q_f, q_r, sc, sc_qss, T, invT, logT);
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
  int kd = 4;
  // Zero ncf
  for (int id = 0; id < kd * 52; ++id) {
    ncf[id] = 0;
  }

  // N2
  ncf[0 * kd + 0] = 2; // N

  // O
  ncf[1 * kd + 1] = 1; // O

  // H2
  ncf[2 * kd + 2] = 2; // H

  // H
  ncf[3 * kd + 2] = 1; // H

  // OH
  ncf[4 * kd + 2] = 1; // H
  ncf[4 * kd + 1] = 1; // O

  // H2O
  ncf[5 * kd + 2] = 2; // H
  ncf[5 * kd + 1] = 1; // O

  // O2
  ncf[6 * kd + 1] = 2; // O

  // HO2
  ncf[7 * kd + 2] = 1; // H
  ncf[7 * kd + 1] = 2; // O

  // H2O2
  ncf[8 * kd + 2] = 2; // H
  ncf[8 * kd + 1] = 2; // O

  // CH
  ncf[9 * kd + 3] = 1; // C
  ncf[9 * kd + 2] = 1; // H

  // HCO
  ncf[10 * kd + 3] = 1; // C
  ncf[10 * kd + 2] = 1; // H
  ncf[10 * kd + 1] = 1; // O

  // CH2
  ncf[11 * kd + 3] = 1; // C
  ncf[11 * kd + 2] = 2; // H

  // CO2
  ncf[12 * kd + 3] = 1; // C
  ncf[12 * kd + 1] = 2; // O

  // CO
  ncf[13 * kd + 3] = 1; // C
  ncf[13 * kd + 1] = 1; // O

  // CH2O
  ncf[14 * kd + 3] = 1; // C
  ncf[14 * kd + 2] = 2; // H
  ncf[14 * kd + 1] = 1; // O

  // CH2GSG
  ncf[15 * kd + 3] = 1; // C
  ncf[15 * kd + 2] = 2; // H

  // CH3
  ncf[16 * kd + 3] = 1; // C
  ncf[16 * kd + 2] = 3; // H

  // CH3O
  ncf[17 * kd + 3] = 1; // C
  ncf[17 * kd + 2] = 3; // H
  ncf[17 * kd + 1] = 1; // O

  // CH4
  ncf[18 * kd + 3] = 1; // C
  ncf[18 * kd + 2] = 4; // H

  // CH3OH
  ncf[19 * kd + 3] = 1; // C
  ncf[19 * kd + 2] = 4; // H
  ncf[19 * kd + 1] = 1; // O

  // C2H6
  ncf[20 * kd + 3] = 2; // C
  ncf[20 * kd + 2] = 6; // H

  // C2H5
  ncf[21 * kd + 3] = 2; // C
  ncf[21 * kd + 2] = 5; // H

  // CH2CO
  ncf[22 * kd + 3] = 2; // C
  ncf[22 * kd + 2] = 2; // H
  ncf[22 * kd + 1] = 1; // O

  // HOCHO
  ncf[23 * kd + 3] = 1; // C
  ncf[23 * kd + 2] = 2; // H
  ncf[23 * kd + 1] = 2; // O

  // CH3O2
  ncf[24 * kd + 3] = 1; // C
  ncf[24 * kd + 2] = 3; // H
  ncf[24 * kd + 1] = 2; // O

  // CH3O2H
  ncf[25 * kd + 3] = 1; // C
  ncf[25 * kd + 2] = 4; // H
  ncf[25 * kd + 1] = 2; // O

  // C2H2
  ncf[26 * kd + 3] = 2; // C
  ncf[26 * kd + 2] = 2; // H

  // HCCO
  ncf[27 * kd + 3] = 2; // C
  ncf[27 * kd + 2] = 1; // H
  ncf[27 * kd + 1] = 1; // O

  // C2H3
  ncf[28 * kd + 3] = 2; // C
  ncf[28 * kd + 2] = 3; // H

  // CH2CHO
  ncf[29 * kd + 3] = 2; // C
  ncf[29 * kd + 2] = 3; // H
  ncf[29 * kd + 1] = 1; // O

  // C3H6
  ncf[30 * kd + 3] = 3; // C
  ncf[30 * kd + 2] = 6; // H

  // C2H4
  ncf[31 * kd + 3] = 2; // C
  ncf[31 * kd + 2] = 4; // H

  // C2H5O
  ncf[32 * kd + 3] = 2; // C
  ncf[32 * kd + 2] = 5; // H
  ncf[32 * kd + 1] = 1; // O

  // CH3CO
  ncf[33 * kd + 3] = 2; // C
  ncf[33 * kd + 2] = 3; // H
  ncf[33 * kd + 1] = 1; // O

  // C2H5O2
  ncf[34 * kd + 3] = 2; // C
  ncf[34 * kd + 2] = 5; // H
  ncf[34 * kd + 1] = 2; // O

  // C3H2
  ncf[35 * kd + 3] = 3; // C
  ncf[35 * kd + 2] = 2; // H

  // C3H3
  ncf[36 * kd + 3] = 3; // C
  ncf[36 * kd + 2] = 3; // H

  // C3H4XA
  ncf[37 * kd + 3] = 3; // C
  ncf[37 * kd + 2] = 4; // H

  // C3H5XA
  ncf[38 * kd + 3] = 3; // C
  ncf[38 * kd + 2] = 5; // H

  // NXC3H7
  ncf[39 * kd + 3] = 3; // C
  ncf[39 * kd + 2] = 7; // H

  // NXC3H7O2
  ncf[40 * kd + 3] = 3; // C
  ncf[40 * kd + 2] = 7; // H
  ncf[40 * kd + 1] = 2; // O

  // C4H6
  ncf[41 * kd + 3] = 4; // C
  ncf[41 * kd + 2] = 6; // H

  // C4H7
  ncf[42 * kd + 3] = 4; // C
  ncf[42 * kd + 2] = 7; // H

  // C4H8X1
  ncf[43 * kd + 3] = 4; // C
  ncf[43 * kd + 2] = 8; // H

  // PXC4H9
  ncf[44 * kd + 3] = 4; // C
  ncf[44 * kd + 2] = 9; // H

  // PXC4H9O2
  ncf[45 * kd + 3] = 4; // C
  ncf[45 * kd + 2] = 9; // H
  ncf[45 * kd + 1] = 2; // O

  // C5H9
  ncf[46 * kd + 3] = 5; // C
  ncf[46 * kd + 2] = 9; // H

  // C5H10X1
  ncf[47 * kd + 3] = 5;  // C
  ncf[47 * kd + 2] = 10; // H

  // C5H11X1
  ncf[48 * kd + 3] = 5;  // C
  ncf[48 * kd + 2] = 11; // H

  // C6H12X1
  ncf[49 * kd + 3] = 6;  // C
  ncf[49 * kd + 2] = 12; // H

  // C7H15X2
  ncf[50 * kd + 3] = 7;  // C
  ncf[50 * kd + 2] = 15; // H

  // NXC7H16
  ncf[51 * kd + 3] = 7;  // C
  ncf[51 * kd + 2] = 16; // H
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
  kname.resize(52);
  kname[0] = "N2";
  kname[1] = "O";
  kname[2] = "H2";
  kname[3] = "H";
  kname[4] = "OH";
  kname[5] = "H2O";
  kname[6] = "O2";
  kname[7] = "HO2";
  kname[8] = "H2O2";
  kname[9] = "CH";
  kname[10] = "HCO";
  kname[11] = "CH2";
  kname[12] = "CO2";
  kname[13] = "CO";
  kname[14] = "CH2O";
  kname[15] = "CH2GSG";
  kname[16] = "CH3";
  kname[17] = "CH3O";
  kname[18] = "CH4";
  kname[19] = "CH3OH";
  kname[20] = "C2H6";
  kname[21] = "C2H5";
  kname[22] = "CH2CO";
  kname[23] = "HOCHO";
  kname[24] = "CH3O2";
  kname[25] = "CH3O2H";
  kname[26] = "C2H2";
  kname[27] = "HCCO";
  kname[28] = "C2H3";
  kname[29] = "CH2CHO";
  kname[30] = "C3H6";
  kname[31] = "C2H4";
  kname[32] = "C2H5O";
  kname[33] = "CH3CO";
  kname[34] = "C2H5O2";
  kname[35] = "C3H2";
  kname[36] = "C3H3";
  kname[37] = "C3H4XA";
  kname[38] = "C3H5XA";
  kname[39] = "NXC3H7";
  kname[40] = "NXC3H7O2";
  kname[41] = "C4H6";
  kname[42] = "C4H7";
  kname[43] = "C4H8X1";
  kname[44] = "PXC4H9";
  kname[45] = "PXC4H9O2";
  kname[46] = "C5H9";
  kname[47] = "C5H10X1";
  kname[48] = "C5H11X1";
  kname[49] = "C6H12X1";
  kname[50] = "C7H15X2";
  kname[51] = "NXC7H16";
}

// compute the sparsity pattern of the chemistry Jacobian
void
SPARSITY_INFO(int* nJdata, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 2809> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 52> conc = {0.0};
  for (int n = 0; n < 52; n++) {
    conc[n] = 1.0 / 52.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 53; k++) {
    for (int l = 0; l < 53; l++) {
      if (Jac[53 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 2809> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 52> conc = {0.0};
  for (int n = 0; n < 52; n++) {
    conc[n] = 1.0 / 52.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 53; k++) {
    for (int l = 0; l < 53; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[53 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 2809> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 52> conc = {0.0};
  for (int n = 0; n < 52; n++) {
    conc[n] = 1.0 / 52.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 53; k++) {
    for (int l = 0; l < 53; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[53 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 2809> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 52> conc = {0.0};
  for (int n = 0; n < 52; n++) {
    conc[n] = 1.0 / 52.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int nc = 0; nc < NCELLS; nc++) {
    int offset_row = nc * 53;
    int offset_col = nc * 53;
    for (int k = 0; k < 53; k++) {
      for (int l = 0; l < 53; l++) {
        if (Jac[53 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 2809> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 52> conc = {0.0};
  for (int n = 0; n < 52; n++) {
    conc[n] = 1.0 / 52.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtrs[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 53;
      for (int l = 0; l < 53; l++) {
        for (int k = 0; k < 53; k++) {
          if (Jac[53 * k + l] != 0.0) {
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
      int offset = nc * 53;
      for (int l = 0; l < 53; l++) {
        for (int k = 0; k < 53; k++) {
          if (Jac[53 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 2809> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 52> conc = {0.0};
  for (int n = 0; n < 52; n++) {
    conc[n] = 1.0 / 52.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 53;
      for (int l = 0; l < 53; l++) {
        for (int k = 0; k < 53; k++) {
          if (k == l) {
            colVals[nJdata_tmp - 1] = l + 1 + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[53 * k + l] != 0.0) {
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
      int offset = nc * 53;
      for (int l = 0; l < 53; l++) {
        for (int k = 0; k < 53; k++) {
          if (k == l) {
            colVals[nJdata_tmp] = l + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[53 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 2809> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 52> conc = {0.0};
  for (int n = 0; n < 52; n++) {
    conc[n] = 1.0 / 52.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int k = 0; k < 53; k++) {
    for (int l = 0; l < 53; l++) {
      if (k == l) {
        rowVals[nJdata_tmp] = l;
        indx[nJdata_tmp] = 53 * k + l;
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[53 * k + l] != 0.0) {
          rowVals[nJdata_tmp] = l;
          indx[nJdata_tmp] = 53 * k + l;
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
  amrex::GpuArray<amrex::Real, 2809> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 52> conc = {0.0};
  for (int n = 0; n < 52; n++) {
    conc[n] = 1.0 / 52.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int l = 0; l < 53; l++) {
      for (int k = 0; k < 53; k++) {
        if (k == l) {
          colVals[nJdata_tmp - 1] = l + 1;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[53 * k + l] != 0.0) {
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
    for (int l = 0; l < 53; l++) {
      for (int k = 0; k < 53; k++) {
        if (k == l) {
          colVals[nJdata_tmp] = l;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[53 * k + l] != 0.0) {
            colVals[nJdata_tmp] = k;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
      }
      rowPtr[l + 1] = nJdata_tmp;
    }
  }
}
