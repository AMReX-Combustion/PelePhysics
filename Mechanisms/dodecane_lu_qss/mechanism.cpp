#include "mechanism.H"
const int rmap[NUM_REACTIONS] = {
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
  255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267};

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
    3, 2, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
    3, 2, 3, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4, 3, 3, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 5, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 3, 3,
    4, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 3, 4, 4, 4, 4, 4, 5, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 4, 4, 4, 4, 4, 5,
    4, 4, 5, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 5, 4, 3, 3, 4, 4,
    3, 4, 4, 3, 4, 4, 3, 4, 4, 3, 4, 4, 3, 4, 4, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    3, 3, 3, 3, 3, 3, 2, 2, 4, 4, 4, 4, 4, 4, 3, 3, 3, 5};
  const int kiv[NUM_GAS_REACTIONS * 5] = {
    1,  8,  4,  0,  0,  3,  7,  0,  0,  0,  12, 5,  11, 0,  0,  1,  37, 11, 0,
    0,  35, 1,  9,  0,  0,  11, 1,  38, 0,  0,  9,  1,  10, 0,  0,  9,  16, 0,
    0,  0,  39, 14, 1,  0,  0,  39, 1,  15, 0,  0,  39, 9,  19, 0,  0,  15, 1,
    40, 0,  0,  40, 1,  16, 0,  0,  39, 40, 22, 0,  0,  1,  18, 19, 0,  0,  9,
    18, 22, 0,  0,  19, 1,  41, 0,  0,  21, 1,  22, 0,  0,  22, 1,  42, 0,  0,
    24, 1,  43, 0,  0,  25, 1,  44, 0,  0,  26, 1,  45, 0,  0,  27, 1,  46, 0,
    0,  28, 1,  29, 0,  0,  30, 1,  47, 0,  0,  12, 2,  13, 0,  0,  1,  5,  0,
    0,  0,  1,  3,  6,  0,  0,  2,  8,  0,  0,  0,  1,  5,  0,  0,  0,  1,  5,
    0,  0,  0,  1,  5,  0,  0,  0,  1,  2,  3,  0,  0,  37, 12, 1,  0,  0,  37,
    12, 1,  0,  0,  1,  8,  2,  3,  0,  5,  2,  1,  3,  0,  5,  3,  1,  6,  0,
    3,  6,  2,  0,  0,  1,  4,  3,  0,  0,  5,  8,  1,  4,  0,  4,  3,  6,  8,
    0,  1,  4,  6,  2,  0,  4,  2,  8,  3,  0,  4,  7,  8,  0,  0,  4,  7,  8,
    0,  0,  1,  7,  6,  3,  0,  1,  7,  5,  4,  0,  7,  2,  4,  3,  0,  7,  3,
    6,  4,  0,  7,  3,  6,  4,  0,  12, 3,  13, 1,  0,  12, 3,  13, 1,  0,  12,
    4,  13, 3,  0,  12, 8,  13, 2,  0,  1,  37, 12, 5,  0,  37, 2,  12, 3,  0,
    37, 2,  13, 1,  0,  37, 3,  12, 6,  0,  37, 8,  12, 4,  0,  35, 2,  1,  37,
    0,  35, 3,  11, 1,  0,  35, 5,  9,  1,  0,  35, 8,  37, 3,  0,  35, 8,  13,
    1,  0,  35, 4,  11, 3,  0,  35, 14, 5,  0,  0,  36, 34, 35, 34, 0,  36, 2,
    12, 5,  0,  36, 2,  1,  37, 0,  36, 3,  11, 1,  0,  36, 5,  9,  1,  0,  36,
    8,  12, 1,  3,  36, 8,  12, 6,  0,  36, 6,  35, 6,  0,  36, 12, 35, 12, 0,
    36, 13, 35, 13, 0,  36, 13, 11, 12, 0,  11, 1,  5,  37, 0,  11, 2,  37, 3,
    0,  11, 3,  6,  37, 0,  11, 8,  37, 4,  0,  11, 4,  7,  37, 0,  9,  2,  11,
    1,  0,  9,  3,  35, 6,  0,  9,  3,  36, 6,  0,  9,  8,  38, 2,  0,  9,  8,
    11, 3,  0,  9,  4,  10, 8,  0,  9,  4,  38, 3,  0,  9,  7,  10, 4,  0,  9,
    37, 10, 12, 0,  11, 9,  10, 37, 0,  35, 9,  15, 1,  0,  36, 9,  15, 1,  0,
    9,  40, 1,  0,  0,  38, 1,  11, 5,  0,  38, 1,  9,  3,  0,  38, 1,  36, 6,
    0,  38, 2,  11, 3,  0,  38, 3,  11, 6,  0,  38, 8,  11, 4,  0,  10, 1,  9,
    5,  0,  10, 2,  9,  3,  0,  10, 3,  9,  6,  0,  35, 10, 9,  0,  0,  36, 10,
    9,  0,  0,  14, 2,  35, 12, 0,  14, 3,  9,  12, 0,  14, 37, 39, 12, 0,  14,
    9,  18, 0,  0,  39, 1,  14, 5,  0,  39, 2,  9,  12, 0,  39, 3,  14, 6,  0,
    39, 8,  14, 4,  0,  39, 8,  17, 2,  0,  39, 8,  11, 37, 0,  39, 4,  17, 3,
    0,  39, 7,  15, 4,  0,  39, 37, 15, 12, 0,  39, 37, 20, 0,  0,  39, 9,  14,
    10, 0,  39, 9,  1,  18, 0,  39, 14, 15, 0,  0,  17, 9,  12, 0,  0,  17, 1,
    9,  37, 0,  17, 8,  11, 12, 3,  15, 1,  39, 5,  0,  15, 2,  39, 3,  0,  15,
    2,  9,  37, 0,  15, 2,  35, 11, 0,  15, 3,  39, 6,  0,  15, 37, 40, 12, 0,
    15, 35, 1,  18, 0,  15, 36, 1,  18, 0,  15, 9,  39, 10, 0,  41, 15, 9,  0,
    0,  15, 8,  39, 4,  0,  39, 15, 21, 0,  0,  40, 1,  15, 5,  0,  40, 2,  11,
    9,  0,  40, 8,  15, 4,  0,  40, 4,  16, 8,  0,  40, 4,  15, 7,  0,  40, 4,
    11, 9,  3,  40, 7,  16, 4,  0,  39, 40, 9,  18, 0,  16, 1,  40, 5,  0,  16,
    2,  40, 3,  0,  16, 3,  40, 6,  0,  16, 36, 40, 9,  0,  16, 9,  40, 10, 0,
    2,  18, 20, 1,  0,  3,  18, 20, 1,  0,  8,  18, 20, 3,  0,  4,  18, 19, 8,
    0,  4,  18, 39, 11, 3,  37, 18, 19, 12, 0,  19, 1,  15, 9,  0,  19, 1,  5,
    18, 0,  19, 2,  20, 1,  0,  19, 2,  40, 37, 0,  19, 2,  3,  18, 0,  19, 3,
    6,  18, 0,  19, 4,  7,  18, 0,  19, 9,  10, 18, 0,  20, 1,  15, 37, 0,  20,
    2,  39, 12, 3,  20, 3,  39, 12, 6,  1,  41, 40, 9,  0,  1,  41, 19, 5,  0,
    2,  41, 40, 11, 0,  3,  41, 19, 6,  0,  8,  41, 19, 4,  0,  4,  41, 40, 11,
    3,  9,  41, 19, 10, 0,  21, 1,  9,  18, 0,  21, 4,  11, 3,  18, 21, 37, 22,
    12, 0,  22, 1,  15, 40, 0,  22, 1,  19, 9,  0,  22, 1,  21, 5,  0,  22, 2,
    37, 41, 0,  22, 2,  21, 3,  0,  22, 2,  21, 3,  0,  22, 3,  21, 6,  0,  22,
    8,  21, 4,  0,  22, 4,  21, 7,  0,  22, 9,  21, 10, 0,  1,  42, 40, 0,  0,
    1,  42, 22, 5,  0,  2,  42, 11, 41, 0,  3,  42, 22, 6,  0,  8,  42, 22, 4,
    0,  4,  42, 11, 3,  41, 9,  42, 22, 10, 0,  23, 15, 18, 0,  0,  23, 39, 19,
    0,  0,  24, 1,  15, 41, 0,  24, 1,  40, 19, 0,  15, 41, 43, 0,  0,  25, 1,
    15, 42, 0,  25, 1,  19, 41, 0,  15, 42, 44, 0,  0,  26, 1,  15, 43, 0,  26,
    1,  19, 42, 0,  15, 43, 45, 0,  0,  27, 1,  15, 44, 0,  27, 1,  19, 43, 0,
    15, 44, 46, 0,  0,  28, 1,  15, 45, 0,  28, 1,  19, 44, 0,  15, 45, 29, 0,
    0,  30, 1,  15, 46, 0,  30, 1,  19, 45, 0,  15, 46, 47, 0,  0,  31, 23, 45,
    0,  0,  15, 47, 48, 0,  0,  48, 50, 0,  0,  0,  19, 29, 49, 0,  0,  22, 46,
    49, 0,  0,  24, 45, 50, 0,  0,  30, 40, 50, 0,  0,  25, 44, 50, 0,  0,  28,
    41, 50, 0,  0,  26, 43, 50, 0,  0,  27, 42, 50, 0,  0,  40, 47, 0,  0,  0,
    29, 41, 0,  0,  0,  46, 42, 0,  0,  0,  43, 45, 0,  0,  0,  44, 0,  0,  0,
    0,  1,  0,  5,  48, 0,  1,  0,  5,  49, 0,  1,  0,  5,  50, 0,  0,  2,  3,
    48, 0,  0,  2,  3,  49, 0,  0,  2,  3,  50, 0,  0,  3,  6,  48, 0,  0,  3,
    6,  49, 0,  0,  3,  6,  50, 0,  0,  8,  4,  48, 0,  0,  8,  4,  49, 0,  0,
    8,  4,  50, 0,  4,  0,  7,  48, 0,  4,  0,  7,  49, 0,  4,  0,  7,  50, 0,
    9,  0,  10, 48, 0,  9,  0,  10, 49, 0,  9,  0,  10, 50, 0,  8,  48, 32, 0,
    0,  32, 8,  48, 0,  0,  8,  49, 32, 0,  0,  32, 8,  49, 0,  0,  8,  50, 32,
    0,  0,  32, 8,  50, 0,  0,  32, 51, 0,  0,  0,  51, 32, 0,  0,  0,  8,  48,
    31, 4,  0,  31, 4,  8,  48, 0,  8,  49, 31, 4,  0,  31, 4,  8,  49, 0,  8,
    50, 31, 4,  0,  31, 4,  8,  50, 0,  51, 8,  52, 0,  0,  52, 51, 8,  0,  0,
    52, 33, 3,  0,  0,  33, 15, 40, 17, 3};
  const int nuv[NUM_GAS_REACTIONS * 5] = {
    -1, -1, 1, 0, 0, -2, 1,  0, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -2, 1,  0, 0, 0,
    -1, 1,  1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -2, 1,  0, 0, 0, -1, -1, 1, 0, 0,
    -2, 1,  0, 0, 0, -2, 1,  0, 0, 0, -2, 1,  0, 0, 0, -2, 1,  0, 0, 0,
    -1, -1, 1, 0, 0, -1, 1,  1, 0, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -2, 1,  1, 0, 0, -1, -1, 2, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -2, 1,  1, 0, 0, -2, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 2, 0, -1, -1, 1, 1, 0, -2, 1,  1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 1, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -2, 1,  1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 2, 0, 0, -1, -1, 2, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -2, 1,  1, 0, 0,
    -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 1, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 1, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 2, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 1, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 2, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 1,
    -1, -1, 1, 1, 1, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 1, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 1, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 2, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 1, -1, -1, 1, 1, 0,
    -1, 1,  1, 0, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0,
    -1, 1,  1, 0, 0, -1, -1, 1, 0, 0, -1, 1,  0, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -2, 1,  0, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, 1,  1, 0, 0,
    -1, -1, 1, 0, 0, -1, 1,  1, 0, 0, -1, -1, 1, 0, 0, -1, 1,  1, 0, 0,
    -1, 1,  0, 0, 0, -1, 1,  0, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, 1,  1, 0, 0, -1, 1,  1, 0, 0, -1, 3,  1, 2, 1};
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

// Returns the progress rates of each reaction
// Given P, T, and mole fractions
void
CKKFKR(
  const amrex::Real P,
  const amrex::Real T,
  const amrex::Real x[],
  amrex::Real q_f[],
  amrex::Real q_r[])
{
  amrex::Real c[NUM_SPECIES]; // temporary storage
  amrex::Real PORT =
    1e6 * P / (8.31446261815324e+07 * T); // convert to SI (mol/cm^3 to mol/m^3)

  // Compute conversion, see Eq 10
  for (int id = 0; id < NUM_GAS_SPECIES; ++id) {
    c[id] = x[id] * PORT;
  }

  // convert to chemkin units
  progressRateFR(q_f, q_r, c, T);

  // convert to chemkin units
  for (int id = 0; id < NUM_GAS_REACTIONS; ++id) {
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
  amrex::Real g_RT[NUM_SPECIES];
  gibbs(g_RT, T);
  amrex::Real g_RT_qss[NUM_GAS_SPECIES];
  gibbs_qss(g_RT_qss, T);

  amrex::Real sc_qss[18];
  // Fill sc_qss here
  amrex::Real kf_qss[198], qf_qss[198], qr_qss[198];
  comp_k_f_qss(T, invT, logT, kf_qss);
  comp_qss_coeff(kf_qss, qf_qss, qr_qss, sc, T, g_RT, g_RT_qss);
  comp_sc_qss(sc_qss, qf_qss, qr_qss);
  comp_qfqr(q_f, q_r, sc, sc_qss, T, invT, logT);
}

// save atomic weights into array
void
atomicWeight(amrex::Real* awt)
{
  awt[0] = 12.011000; // C
  awt[1] = 1.008000;  // H
  awt[2] = 15.999000; // O
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
  for (int id = 0; id < kd * NUM_GAS_SPECIES; ++id) {
    ncf[id] = 0;
  }

  // NC12H26
  ncf[0 * kd + 0] = 12; // C
  ncf[0 * kd + 1] = 26; // H

  // H
  ncf[1 * kd + 1] = 1; // H

  // O
  ncf[2 * kd + 2] = 1; // O

  // OH
  ncf[3 * kd + 1] = 1; // H
  ncf[3 * kd + 2] = 1; // O

  // HO2
  ncf[4 * kd + 1] = 1; // H
  ncf[4 * kd + 2] = 2; // O

  // H2
  ncf[5 * kd + 1] = 2; // H

  // H2O
  ncf[6 * kd + 1] = 2; // H
  ncf[6 * kd + 2] = 1; // O

  // H2O2
  ncf[7 * kd + 1] = 2; // H
  ncf[7 * kd + 2] = 2; // O

  // O2
  ncf[8 * kd + 2] = 2; // O

  // CH3
  ncf[9 * kd + 0] = 1; // C
  ncf[9 * kd + 1] = 3; // H

  // CH4
  ncf[10 * kd + 0] = 1; // C
  ncf[10 * kd + 1] = 4; // H

  // CH2O
  ncf[11 * kd + 0] = 1; // C
  ncf[11 * kd + 1] = 2; // H
  ncf[11 * kd + 2] = 1; // O

  // CO
  ncf[12 * kd + 0] = 1; // C
  ncf[12 * kd + 2] = 1; // O

  // CO2
  ncf[13 * kd + 0] = 1; // C
  ncf[13 * kd + 2] = 2; // O

  // C2H2
  ncf[14 * kd + 0] = 2; // C
  ncf[14 * kd + 1] = 2; // H

  // C2H4
  ncf[15 * kd + 0] = 2; // C
  ncf[15 * kd + 1] = 4; // H

  // C2H6
  ncf[16 * kd + 0] = 2; // C
  ncf[16 * kd + 1] = 6; // H

  // CH2CHO
  ncf[17 * kd + 0] = 2; // C
  ncf[17 * kd + 1] = 3; // H
  ncf[17 * kd + 2] = 1; // O

  // aC3H5
  ncf[18 * kd + 0] = 3; // C
  ncf[18 * kd + 1] = 5; // H

  // C3H6
  ncf[19 * kd + 0] = 3; // C
  ncf[19 * kd + 1] = 6; // H

  // C2H3CHO
  ncf[20 * kd + 0] = 3; // C
  ncf[20 * kd + 1] = 4; // H
  ncf[20 * kd + 2] = 1; // O

  // C4H7
  ncf[21 * kd + 0] = 4; // C
  ncf[21 * kd + 1] = 7; // H

  // C4H81
  ncf[22 * kd + 0] = 4; // C
  ncf[22 * kd + 1] = 8; // H

  // C5H9
  ncf[23 * kd + 0] = 5; // C
  ncf[23 * kd + 1] = 9; // H

  // C5H10
  ncf[24 * kd + 0] = 5;  // C
  ncf[24 * kd + 1] = 10; // H

  // C6H12
  ncf[25 * kd + 0] = 6;  // C
  ncf[25 * kd + 1] = 12; // H

  // C7H14
  ncf[26 * kd + 0] = 7;  // C
  ncf[26 * kd + 1] = 14; // H

  // C8H16
  ncf[27 * kd + 0] = 8;  // C
  ncf[27 * kd + 1] = 16; // H

  // C9H18
  ncf[28 * kd + 0] = 9;  // C
  ncf[28 * kd + 1] = 18; // H

  // PXC9H19
  ncf[29 * kd + 0] = 9;  // C
  ncf[29 * kd + 1] = 19; // H

  // C10H20
  ncf[30 * kd + 0] = 10; // C
  ncf[30 * kd + 1] = 20; // H

  // C12H24
  ncf[31 * kd + 0] = 12; // C
  ncf[31 * kd + 1] = 24; // H

  // C12H25O2
  ncf[32 * kd + 0] = 12; // C
  ncf[32 * kd + 1] = 25; // H
  ncf[32 * kd + 2] = 2;  // O

  // OC12H23OOH
  ncf[33 * kd + 0] = 12; // C
  ncf[33 * kd + 1] = 24; // H
  ncf[33 * kd + 2] = 3;  // O

  // N2
  ncf[34 * kd + 3] = 2; // N
}

// Returns the vector of strings of element names
void
CKSYME_STR(amrex::Vector<std::string>& ename)
{
  ename.resize(NUM_ELEMENTS);
  ename[0] = "C";
  ename[1] = "H";
  ename[2] = "O";
  ename[3] = "N";
}

// Returns the vector of strings of species names
void
CKSYMS_STR(amrex::Vector<std::string>& kname)
{
  kname.resize(NUM_SPECIES);
  kname[0] = "NC12H26";
  kname[1] = "H";
  kname[2] = "O";
  kname[3] = "OH";
  kname[4] = "HO2";
  kname[5] = "H2";
  kname[6] = "H2O";
  kname[7] = "H2O2";
  kname[8] = "O2";
  kname[9] = "CH3";
  kname[10] = "CH4";
  kname[11] = "CH2O";
  kname[12] = "CO";
  kname[13] = "CO2";
  kname[14] = "C2H2";
  kname[15] = "C2H4";
  kname[16] = "C2H6";
  kname[17] = "CH2CHO";
  kname[18] = "aC3H5";
  kname[19] = "C3H6";
  kname[20] = "C2H3CHO";
  kname[21] = "C4H7";
  kname[22] = "C4H81";
  kname[23] = "C5H9";
  kname[24] = "C5H10";
  kname[25] = "C6H12";
  kname[26] = "C7H14";
  kname[27] = "C8H16";
  kname[28] = "C9H18";
  kname[29] = "PXC9H19";
  kname[30] = "C10H20";
  kname[31] = "C12H24";
  kname[32] = "C12H25O2";
  kname[33] = "OC12H23OOH";
  kname[34] = "N2";
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
