#include "mechanism.H"
const int rmap[NUM_REACTIONS] = {
  9,   11,  48,  55,  60,  61,  68,  70,  81,  111, 124, 128, 139, 148, 151,
  156, 157, 171, 179, 199, 207, 211, 212, 279, 282, 24,  4,   5,   6,   7,
  8,   32,  33,  0,   1,   2,   3,   10,  12,  13,  14,  15,  16,  17,  18,
  19,  20,  21,  22,  23,  25,  26,  27,  28,  29,  30,  31,  34,  35,  36,
  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,  49,  50,  51,  52,
  53,  54,  56,  57,  58,  59,  62,  63,  64,  65,  66,  67,  69,  71,  72,
  73,  74,  75,  76,  77,  78,  79,  80,  82,  83,  84,  85,  86,  87,  88,
  89,  90,  91,  92,  93,  94,  95,  96,  97,  98,  99,  100, 101, 102, 103,
  104, 105, 106, 107, 108, 109, 110, 112, 113, 114, 115, 116, 117, 118, 119,
  120, 121, 122, 123, 125, 126, 127, 129, 130, 131, 132, 133, 134, 135, 136,
  137, 138, 140, 141, 142, 143, 144, 145, 146, 147, 149, 150, 152, 153, 154,
  155, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 172,
  173, 174, 175, 176, 177, 178, 180, 181, 182, 183, 184, 185, 186, 187, 188,
  189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 200, 201, 202, 203, 204,
  205, 206, 208, 209, 210, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222,
  223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237,
  238, 239, 240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252,
  253, 254, 255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267,
  268, 269, 270, 271, 272, 273, 274, 275, 276, 277, 278, 280, 281, 283, 284,
  285, 286, 287, 288};

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
    4, 4, 4, 3, 2, 2, 3, 3, 2, 3, 4, 2, 4, 3, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 3,
    4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4,
    4, 4, 4, 5, 4, 3, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 3, 4, 3, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3,
    4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 3, 4, 4, 4, 4, 4, 5, 3, 4,
    5, 3, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 4, 4, 4,
    4, 4, 4, 5, 3, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 3, 3, 4, 4, 2, 4, 3, 4, 4, 3,
    4, 4, 4, 4, 4, 5, 4, 3, 4, 4, 5, 3, 3, 4, 4, 4, 5, 4, 4, 4, 4, 4, 4, 4, 5,
    4, 4, 4, 4, 4, 5, 4, 4, 4, 4, 5, 4, 3, 4, 4, 2, 4, 4, 3, 3, 3, 4, 5, 4, 4,
    4, 4, 3, 4, 4, 5, 4, 3, 3, 5, 6, 6, 6, 6, 5, 6, 6, 6, 6, 5, 6, 6, 6, 6, 5,
    6, 6, 6, 6, 4, 4, 4, 4, 4, 5, 5, 6, 5, 4};
  const int kiv[NUM_GAS_REACTIONS * 6] = {
    0,  7,  1,  2,  0,  0,  4,  1,  0,  2,  0,  0,  4,  2,  0,  5,  0,  0,  2,
    5,  1,  0,  0,  0,  0,  4,  0,  0,  0,  0,  0,  4,  0,  0,  0,  0,  0,  2,
    5,  0,  0,  0,  0,  1,  2,  0,  0,  0,  1,  7,  0,  0,  0,  0,  0,  7,  3,
    0,  0,  0,  4,  7,  0,  3,  0,  0,  2,  6,  0,  0,  0,  0,  0,  3,  5,  1,
    0,  0,  0,  3,  2,  0,  0,  0,  3,  1,  7,  2,  0,  0,  3,  2,  5,  7,  0,
    0,  3,  2,  5,  7,  0,  0,  3,  6,  7,  0,  0,  0,  3,  6,  7,  0,  0,  0,
    0,  6,  4,  3,  0,  0,  0,  6,  5,  2,  0,  0,  6,  1,  3,  2,  0,  0,  6,
    2,  5,  3,  0,  0,  6,  2,  5,  3,  0,  0,  18, 1,  19, 0,  0,  0,  18, 2,
    19, 0,  0,  0,  18, 2,  19, 0,  0,  0,  18, 3,  19, 2,  0,  0,  0,  13, 18,
    4,  0,  0,  13, 1,  18, 2,  0,  0,  13, 1,  19, 0,  0,  0,  13, 2,  18, 5,
    0,  0,  13, 18, 0,  0,  0,  0,  13, 18, 0,  0,  0,  0,  13, 7,  18, 3,  0,
    0,  8,  1,  18, 0,  0,  0,  8,  2,  0,  13, 0,  0,  8,  4,  9,  0,  0,  0,
    8,  5,  14, 0,  0,  0,  8,  7,  13, 1,  0,  0,  8,  19, 18, 13, 0,  0,  9,
    1,  0,  13, 0,  0,  9,  2,  14, 0,  0,  0,  9,  2,  8,  5,  0,  0,  9,  4,
    11, 0,  0,  0,  9,  7,  13, 2,  0,  0,  9,  7,  19, 0,  0,  0,  9,  3,  14,
    2,  0,  0,  9,  18, 27, 0,  0,  0,  10, 55, 9,  55, 0,  0,  10, 0,  8,  4,
    0,  0,  10, 2,  14, 0,  0,  0,  10, 4,  11, 0,  0,  0,  10, 7,  18, 0,  2,
    0,  10, 7,  18, 5,  0,  0,  10, 5,  17, 0,  0,  0,  10, 5,  9,  5,  0,  0,
    10, 18, 9,  18, 0,  0,  10, 19, 9,  19, 0,  0,  10, 19, 14, 18, 0,  0,  14,
    0,  16, 0,  0,  0,  14, 0,  15, 0,  0,  0,  14, 0,  4,  13, 0,  0,  14, 1,
    13, 2,  0,  0,  14, 2,  5,  13, 0,  0,  14, 7,  13, 3,  0,  0,  14, 3,  6,
    13, 0,  0,  8,  14, 27, 0,  0,  0,  11, 0,  12, 0,  0,  0,  11, 1,  14, 0,
    0,  0,  11, 2,  17, 0,  0,  0,  11, 2,  9,  5,  0,  0,  11, 2,  10, 5,  0,
    0,  11, 7,  15, 1,  0,  0,  11, 7,  14, 2,  0,  0,  11, 3,  12, 7,  0,  0,
    11, 3,  15, 2,  0,  0,  8,  11, 22, 0,  0,  0,  11, 13, 12, 18, 0,  0,  14,
    11, 12, 13, 0,  0,  9,  11, 23, 0,  0,  0,  11, 25, 0,  0,  0,  0,  11, 24,
    0,  0,  0,  0,  11, 26, 23, 18, 0,  0,  15, 0,  14, 4,  0,  0,  15, 0,  11,
    2,  0,  0,  15, 0,  10, 5,  0,  0,  15, 2,  14, 5,  0,  0,  15, 7,  14, 3,
    0,  0,  16, 0,  14, 4,  0,  0,  16, 0,  11, 2,  0,  0,  16, 0,  10, 5,  0,
    0,  16, 7,  14, 3,  0,  0,  12, 0,  11, 4,  0,  0,  12, 1,  11, 2,  0,  0,
    12, 2,  11, 5,  0,  0,  8,  12, 23, 0,  0,  0,  9,  12, 11, 0,  0,  0,  10,
    12, 11, 0,  0,  0,  17, 0,  16, 4,  0,  0,  17, 0,  15, 4,  0,  0,  17, 1,
    16, 2,  0,  0,  17, 2,  16, 5,  0,  0,  17, 2,  15, 5,  0,  0,  20, 1,  8,
    18, 0,  0,  20, 2,  0,  26, 0,  0,  20, 7,  18, 13, 0,  0,  20, 4,  21, 0,
    0,  0,  0,  26, 10, 18, 0,  0,  26, 1,  18, 0,  0,  0,  26, 7,  18, 2,  0,
    0,  22, 21, 0,  0,  0,  0,  21, 1,  20, 2,  0,  0,  21, 1,  9,  18, 0,  0,
    21, 1,  0,  26, 0,  0,  21, 2,  27, 0,  0,  0,  21, 2,  27, 0,  0,  0,  21,
    2,  20, 5,  0,  0,  21, 13, 22, 18, 0,  0,  21, 9,  31, 0,  0,  0,  21, 10,
    31, 0,  0,  0,  20, 21, 40, 0,  0,  0,  21, 11, 0,  32, 0,  0,  21, 11, 0,
    33, 0,  0,  27, 0,  29, 0,  0,  0,  27, 0,  4,  26, 0,  0,  27, 0,  11, 18,
    0,  0,  27, 2,  5,  26, 0,  0,  22, 0,  23, 0,  0,  0,  22, 0,  21, 4,  0,
    0,  22, 1,  27, 0,  0,  0,  22, 1,  11, 18, 0,  0,  22, 2,  21, 5,  0,  0,
    22, 7,  21, 3,  0,  0,  22, 7,  29, 1,  0,  0,  22, 7,  14, 13, 0,  0,  22,
    3,  29, 2,  0,  0,  22, 13, 23, 18, 0,  0,  22, 13, 39, 0,  0,  0,  22, 11,
    36, 0,  0,  0,  22, 11, 0,  34, 0,  0,  29, 11, 18, 0,  0,  0,  29, 0,  28,
    0,  0,  0,  29, 0,  11, 13, 0,  0,  29, 0,  27, 4,  0,  0,  29, 2,  27, 5,
    0,  0,  29, 7,  27, 3,  0,  0,  29, 7,  14, 18, 2,  0,  11, 18, 28, 0,  0,
    0,  28, 0,  11, 13, 0,  0,  28, 3,  11, 19, 2,  0,  11, 13, 30, 0,  0,  0,
    30, 0,  28, 4,  0,  0,  30, 2,  28, 5,  0,  0,  11, 30, 28, 12, 0,  0,  30,
    7,  28, 3,  0,  0,  23, 21, 4,  0,  0,  0,  23, 0,  24, 0,  0,  0,  23, 0,
    22, 4,  0,  0,  23, 1,  22, 2,  0,  0,  23, 1,  11, 13, 0,  0,  23, 1,  9,
    14, 0,  0,  23, 2,  22, 5,  0,  0,  23, 13, 24, 18, 0,  0,  23, 8,  0,  33,
    0,  0,  23, 8,  0,  32, 0,  0,  23, 9,  0,  34, 0,  0,  23, 10, 0,  34, 0,
    0,  23, 11, 22, 12, 0,  0,  23, 11, 37, 0,  0,  0,  22, 23, 48, 0,  0,  0,
    24, 0,  25, 0,  0,  0,  24, 0,  23, 4,  0,  0,  24, 1,  14, 11, 0,  0,  24,
    1,  30, 0,  0,  0,  24, 7,  23, 3,  0,  0,  24, 3,  25, 7,  0,  0,  24, 3,
    23, 6,  0,  0,  24, 3,  14, 11, 2,  0,  22, 24, 49, 0,  0,  0,  22, 24, 11,
    34, 0,  0,  25, 0,  24, 4,  0,  0,  25, 1,  24, 2,  0,  0,  25, 2,  24, 5,
    0,  0,  25, 10, 24, 11, 0,  0,  25, 11, 24, 12, 0,  0,  31, 0,  32, 0,  0,
    0,  31, 1,  20, 14, 0,  0,  31, 7,  27, 13, 0,  0,  31, 3,  7,  32, 0,  0,
    0,  33, 35, 0,  0,  0,  0,  33, 34, 0,  0,  0,  1,  33, 23, 18, 0,  0,  2,
    33, 31, 5,  0,  0,  32, 33, 0,  0,  0,  0,  0,  32, 0,  33, 0,  0,  0,  32,
    35, 0,  0,  0,  1,  32, 23, 18, 0,  0,  2,  32, 31, 5,  0,  0,  0,  34, 36,
    0,  0,  0,  0,  34, 4,  33, 0,  0,  1,  34, 39, 0,  0,  0,  2,  34, 39, 0,
    0,  0,  2,  34, 5,  33, 0,  0,  3,  34, 36, 7,  0,  0,  3,  34, 22, 14, 2,
    0,  13, 34, 36, 18, 0,  0,  11, 34, 49, 0,  0,  0,  11, 34, 12, 33, 0,  0,
    35, 7,  14, 28, 0,  0,  35, 3,  27, 11, 2,  0,  36, 0,  37, 0,  0,  0,  36,
    0,  38, 0,  0,  0,  36, 0,  23, 11, 0,  0,  36, 0,  4,  34, 0,  0,  36, 0,
    35, 4,  0,  0,  36, 1,  27, 11, 0,  0,  36, 1,  39, 0,  0,  0,  36, 1,  24,
    13, 0,  0,  36, 1,  2,  34, 0,  0,  36, 1,  35, 2,  0,  0,  36, 2,  5,  34,
    0,  0,  36, 2,  35, 5,  0,  0,  36, 11, 12, 34, 0,  0,  39, 1,  22, 18, 2,
    0,  39, 1,  27, 14, 0,  0,  0,  38, 24, 11, 0,  0,  1,  38, 11, 30, 0,  0,
    2,  38, 36, 5,  0,  0,  7,  38, 36, 3,  0,  0,  3,  38, 11, 30, 2,  0,  11,
    38, 36, 12, 0,  0,  0,  37, 24, 11, 0,  0,  2,  37, 36, 5,  0,  0,  7,  37,
    36, 3,  0,  0,  3,  37, 24, 14, 2,  0,  11, 37, 36, 12, 0,  0,  40, 0,  41,
    0,  0,  0,  0,  41, 40, 4,  0,  0,  42, 2,  5,  41, 0,  0,  44, 43, 0,  0,
    0,  0,  45, 0,  22, 23, 0,  0,  45, 2,  5,  43, 0,  0,  46, 0,  43, 0,  0,
    0,  47, 44, 0,  0,  0,  0,  48, 45, 0,  0,  0,  0,  48, 7,  45, 3,  0,  0,
    48, 3,  14, 2,  34, 0,  49, 0,  23, 24, 0,  0,  49, 0,  36, 11, 0,  0,  49,
    0,  48, 4,  0,  0,  49, 1,  13, 37, 0,  0,  23, 24, 50, 0,  0,  0,  2,  50,
    49, 5,  0,  0,  7,  50, 49, 3,  0,  0,  3,  50, 14, 2,  37, 0,  11, 50, 49,
    12, 0,  0,  51, 23, 37, 0,  0,  0,  51, 23, 50, 0,  0,  0,  0,  51, 23, 4,
    50, 0,  0,  51, 23, 49, 4,  50, 0,  51, 36, 52, 4,  37, 0,  51, 23, 54, 4,
    37, 0,  51, 23, 52, 4,  50, 11, 51, 23, 12, 50, 0,  11, 51, 23, 49, 12, 50,
    11, 51, 36, 52, 12, 37, 11, 51, 23, 54, 12, 37, 11, 51, 23, 52, 12, 50, 51,
    1,  23, 2,  50, 0,  51, 1,  23, 49, 2,  50, 51, 1,  36, 52, 2,  37, 51, 1,
    23, 54, 2,  37, 51, 1,  23, 52, 2,  50, 51, 2,  23, 5,  50, 0,  51, 2,  23,
    49, 5,  50, 51, 2,  36, 52, 5,  37, 51, 2,  23, 54, 5,  37, 51, 2,  23, 52,
    5,  50, 52, 0,  36, 37, 0,  0,  52, 0,  23, 50, 0,  0,  52, 0,  53, 4,  0,
    0,  54, 0,  24, 36, 0,  0,  54, 0,  23, 37, 0,  0,  54, 0,  23, 4,  34, 0,
    53, 0,  23, 11, 34, 0,  53, 3,  23, 14, 2,  34, 52, 1,  23, 13, 37, 0,  54,
    1,  13, 50, 0,  0};
  const int nuv[NUM_GAS_REACTIONS * 6] = {
    -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -2, 1,  1,
    0,  0,  0, -2, 1,  0, 0,  0,  0, -2, 1,  0, 0,  0,  0, -1, -1, 1, 0,  0,  0,
    -1, -1, 1, 0,  0,  0, -2, 1,  0, 0,  0,  0, -1, -1, 1, 0,  0,  0, -1, -1, 1,
    1,  0,  0, -2, 1,  0, 0,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 2, 0,  0,  0,
    -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -2, 1,  1,
    0,  0,  0, -2, 1,  1, 0,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0,
    -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1,
    0,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0,
    -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1,
    1,  0,  0, -1, 1,  1, 0,  0,  0, -1, 1,  1, 0,  0,  0, -1, -1, 1, 1,  0,  0,
    -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1,
    1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0,
    -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1,
    1,  0,  0, -1, -1, 1, 2,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 0,  0,  0,
    -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1,
    1,  0,  0, -1, -1, 1, 1,  1,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 0,  0,  0,
    -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1,
    1,  0,  0, -1, -1, 1, 0,  0,  0, -1, -1, 1, 0,  0,  0, -1, -1, 1, 1,  0,  0,
    -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1,
    1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 0,  0,  0, -1, -1, 1, 1,  0,  0,
    -1, -1, 1, 0,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1,
    1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0,
    -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1,
    1,  0,  0, -2, 1,  0, 0,  0,  0, -2, 1,  1, 0,  0,  0, -1, -1, 1, 1,  0,  0,
    -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1,
    1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0,
    -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1,
    1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 2, 0,  0,  0,
    -1, -1, 2, 0,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1,
    1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0,
    -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1,
    1,  0,  0, -1, -1, 2, 1,  0,  0, -1, -1, 2, 1,  0,  0, -1, 1,  1, 0,  0,  0,
    -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1,
    1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0,
    -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1,
    1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 0,  0,  0, -1, -1, 1, 1,  0,  0,
    -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 0,  0,  0, -1, -1, 1,
    1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0,
    -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1,
    1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 0,  0,  0, -1, -1, 1, 0,  0,  0,
    -1, -1, 1, 1,  0,  0, -1, 1,  1, 0,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1,
    1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0,
    -1, -1, 1, 1,  1,  0, -1, -1, 1, 0,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1,
    1,  1,  0, -1, -1, 1, 0,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0,
    -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, 1,  1, 0,  0,  0, -1, -1, 1,
    0,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0,
    -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1,
    1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0,
    -1, -1, 1, 1,  0,  0, -1, -1, 1, 0,  0,  0, -1, -1, 1, 0,  0,  0, -1, -1, 1,
    0,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0,
    -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1,
    1,  1,  0, -1, -1, 1, 0,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0,
    -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1,
    1,  0,  0, -1, -1, 1, 0,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0,
    -1, -1, 1, 1,  0,  0, -1, -1, 1, 0,  0,  0, -1, -1, 1, 0,  0,  0, -1, -1, 1,
    1,  0,  0, -1, -1, 1, 1,  0,  0, -1, 1,  0, 0,  0,  0, -1, -1, 1, 1,  0,  0,
    -1, -1, 1, 0,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1,
    0,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 2,  0,  0,
    -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  1,  0, -1, -1, 1,
    1,  0,  0, -1, -1, 1, 0,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0,
    -1, -1, 1, 1,  1,  0, -1, -1, 1, 0,  0,  0, -1, -1, 1, 0,  0,  0, -1, -1, 1,
    1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  1,  0,
    -1, -1, 1, 2,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1,
    1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0,
    -1, -1, 1, 1,  1,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1,
    1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  1,  0,
    -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1,
    1,  0,  0, -1, -1, 1, 1,  1,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 0,  0,  0,
    -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, 1,  0, 0,  0,  0, -1, -1, 1,
    1,  0,  0, -1, -1, 1, 1,  0,  0, -1, 1,  1, 0,  0,  0, -1, 1,  1, 0,  0,  0,
    -1, 1,  1, 0,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  1,  0, -1, -1, 1,
    1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0,
    -1, -1, 1, 0,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1,
    1,  1,  0, -1, -1, 1, 1,  0,  0, -1, 3,  2, 0,  0,  0, -1, 2,  2, 0,  0,  0,
    -1, -1, 4, 1,  1,  0, -1, -1, 2, 1,  1,  1, -1, -1, 1, 1,  1,  1, -1, -1, 2,
    1,  1,  1, -1, -1, 1, 1,  1,  1, -1, -1, 4, 1,  1,  0, -1, -1, 2, 1,  1,  1,
    -1, -1, 1, 1,  1,  1, -1, -1, 2, 1,  1,  1, -1, -1, 1, 1,  1,  1, -1, -1, 4,
    1,  1,  0, -1, -1, 2, 1,  1,  1, -1, -1, 1, 1,  1,  1, -1, -1, 2, 1,  1,  1,
    -1, -1, 1, 1,  1,  1, -1, -1, 4, 1,  1,  0, -1, -1, 2, 1,  1,  1, -1, -1, 1,
    1,  1,  1, -1, -1, 2, 1,  1,  1, -1, -1, 1, 1,  1,  1, -1, -1, 1, 1,  0,  0,
    -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1, 1,  0,  0, -1, -1, 1,
    1,  0,  0, -1, -1, 1, 1,  1,  0, -1, -1, 1, 1,  1,  0, -1, -1, 1, 1,  1,  1,
    -1, -1, 1, 1,  1,  0, -1, -1, 1, 1,  0,  0};
  if (i < 1) {
    // Return max num species per reaction
    nspec = 6;
  } else {
    if (i > NUM_GAS_REACTIONS) {
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
  for (int id = 0; id < kd * NUM_GAS_SPECIES; ++id) {
    ncf[id] = 0;
  }

  // H
  ncf[0 * kd + 1] = 1; // H

  // O
  ncf[1 * kd + 0] = 1; // O

  // OH
  ncf[2 * kd + 1] = 1; // H
  ncf[2 * kd + 0] = 1; // O

  // HO2
  ncf[3 * kd + 1] = 1; // H
  ncf[3 * kd + 0] = 2; // O

  // H2
  ncf[4 * kd + 1] = 2; // H

  // H2O
  ncf[5 * kd + 1] = 2; // H
  ncf[5 * kd + 0] = 1; // O

  // H2O2
  ncf[6 * kd + 1] = 2; // H
  ncf[6 * kd + 0] = 2; // O

  // O2
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

  // HCO
  ncf[13 * kd + 2] = 1; // C
  ncf[13 * kd + 1] = 1; // H
  ncf[13 * kd + 0] = 1; // O

  // CH2O
  ncf[14 * kd + 2] = 1; // C
  ncf[14 * kd + 1] = 2; // H
  ncf[14 * kd + 0] = 1; // O

  // CH3O
  ncf[15 * kd + 2] = 1; // C
  ncf[15 * kd + 1] = 3; // H
  ncf[15 * kd + 0] = 1; // O

  // CH2OH
  ncf[16 * kd + 2] = 1; // C
  ncf[16 * kd + 1] = 3; // H
  ncf[16 * kd + 0] = 1; // O

  // CH3OH
  ncf[17 * kd + 2] = 1; // C
  ncf[17 * kd + 1] = 4; // H
  ncf[17 * kd + 0] = 1; // O

  // CO
  ncf[18 * kd + 2] = 1; // C
  ncf[18 * kd + 0] = 1; // O

  // CO2
  ncf[19 * kd + 2] = 1; // C
  ncf[19 * kd + 0] = 2; // O

  // C2H
  ncf[20 * kd + 2] = 2; // C
  ncf[20 * kd + 1] = 1; // H

  // C2H2
  ncf[21 * kd + 2] = 2; // C
  ncf[21 * kd + 1] = 2; // H

  // C2H3
  ncf[22 * kd + 2] = 2; // C
  ncf[22 * kd + 1] = 3; // H

  // C2H4
  ncf[23 * kd + 2] = 2; // C
  ncf[23 * kd + 1] = 4; // H

  // C2H5
  ncf[24 * kd + 2] = 2; // C
  ncf[24 * kd + 1] = 5; // H

  // C2H6
  ncf[25 * kd + 2] = 2; // C
  ncf[25 * kd + 1] = 6; // H

  // HCCO
  ncf[26 * kd + 2] = 2; // C
  ncf[26 * kd + 1] = 1; // H
  ncf[26 * kd + 0] = 1; // O

  // CH2CO
  ncf[27 * kd + 2] = 2; // C
  ncf[27 * kd + 1] = 2; // H
  ncf[27 * kd + 0] = 1; // O

  // CH3CO
  ncf[28 * kd + 2] = 2; // C
  ncf[28 * kd + 1] = 3; // H
  ncf[28 * kd + 0] = 1; // O

  // CH2CHO
  ncf[29 * kd + 2] = 2; // C
  ncf[29 * kd + 1] = 3; // H
  ncf[29 * kd + 0] = 1; // O

  // CH3CHO
  ncf[30 * kd + 2] = 2; // C
  ncf[30 * kd + 1] = 4; // H
  ncf[30 * kd + 0] = 1; // O

  // C3H3
  ncf[31 * kd + 2] = 3; // C
  ncf[31 * kd + 1] = 3; // H

  // pC3H4
  ncf[32 * kd + 2] = 3; // C
  ncf[32 * kd + 1] = 4; // H

  // aC3H4
  ncf[33 * kd + 2] = 3; // C
  ncf[33 * kd + 1] = 4; // H

  // aC3H5
  ncf[34 * kd + 2] = 3; // C
  ncf[34 * kd + 1] = 5; // H

  // CH3CCH2
  ncf[35 * kd + 2] = 3; // C
  ncf[35 * kd + 1] = 5; // H

  // C3H6
  ncf[36 * kd + 2] = 3; // C
  ncf[36 * kd + 1] = 6; // H

  // nC3H7
  ncf[37 * kd + 2] = 3; // C
  ncf[37 * kd + 1] = 7; // H

  // iC3H7
  ncf[38 * kd + 2] = 3; // C
  ncf[38 * kd + 1] = 7; // H

  // C2H3CHO
  ncf[39 * kd + 2] = 3; // C
  ncf[39 * kd + 1] = 4; // H
  ncf[39 * kd + 0] = 1; // O

  // C4H2
  ncf[40 * kd + 2] = 4; // C
  ncf[40 * kd + 1] = 2; // H

  // iC4H3
  ncf[41 * kd + 2] = 4; // C
  ncf[41 * kd + 1] = 3; // H

  // C4H4
  ncf[42 * kd + 2] = 4; // C
  ncf[42 * kd + 1] = 4; // H

  // iC4H5
  ncf[43 * kd + 2] = 4; // C
  ncf[43 * kd + 1] = 5; // H

  // C4H5-2
  ncf[44 * kd + 2] = 4; // C
  ncf[44 * kd + 1] = 5; // H

  // C4H6
  ncf[45 * kd + 2] = 4; // C
  ncf[45 * kd + 1] = 6; // H

  // C4H612
  ncf[46 * kd + 2] = 4; // C
  ncf[46 * kd + 1] = 6; // H

  // C4H6-2
  ncf[47 * kd + 2] = 4; // C
  ncf[47 * kd + 1] = 6; // H

  // C4H7
  ncf[48 * kd + 2] = 4; // C
  ncf[48 * kd + 1] = 7; // H

  // C4H81
  ncf[49 * kd + 2] = 4; // C
  ncf[49 * kd + 1] = 8; // H

  // pC4H9
  ncf[50 * kd + 2] = 4; // C
  ncf[50 * kd + 1] = 9; // H

  // NC12H26
  ncf[51 * kd + 2] = 12; // C
  ncf[51 * kd + 1] = 26; // H

  // C6H12
  ncf[52 * kd + 2] = 6;  // C
  ncf[52 * kd + 1] = 12; // H

  // C6H11
  ncf[53 * kd + 2] = 6;  // C
  ncf[53 * kd + 1] = 11; // H

  // C5H10
  ncf[54 * kd + 2] = 5;  // C
  ncf[54 * kd + 1] = 10; // H

  // N2
  ncf[55 * kd + 3] = 2; // N
}

// Returns the vector of strings of element names
void
CKSYME_STR(amrex::Vector<std::string>& ename)
{
  ename.resize(NUM_ELEMENTS);
  ename[0] = "O";
  ename[1] = "H";
  ename[2] = "C";
  ename[3] = "N";
}

// Returns the vector of strings of species names
void
CKSYMS_STR(amrex::Vector<std::string>& kname)
{
  kname.resize(NUM_SPECIES);
  kname[0] = "H";
  kname[1] = "O";
  kname[2] = "OH";
  kname[3] = "HO2";
  kname[4] = "H2";
  kname[5] = "H2O";
  kname[6] = "H2O2";
  kname[7] = "O2";
  kname[8] = "CH";
  kname[9] = "CH2";
  kname[10] = "CH2*";
  kname[11] = "CH3";
  kname[12] = "CH4";
  kname[13] = "HCO";
  kname[14] = "CH2O";
  kname[15] = "CH3O";
  kname[16] = "CH2OH";
  kname[17] = "CH3OH";
  kname[18] = "CO";
  kname[19] = "CO2";
  kname[20] = "C2H";
  kname[21] = "C2H2";
  kname[22] = "C2H3";
  kname[23] = "C2H4";
  kname[24] = "C2H5";
  kname[25] = "C2H6";
  kname[26] = "HCCO";
  kname[27] = "CH2CO";
  kname[28] = "CH3CO";
  kname[29] = "CH2CHO";
  kname[30] = "CH3CHO";
  kname[31] = "C3H3";
  kname[32] = "pC3H4";
  kname[33] = "aC3H4";
  kname[34] = "aC3H5";
  kname[35] = "CH3CCH2";
  kname[36] = "C3H6";
  kname[37] = "nC3H7";
  kname[38] = "iC3H7";
  kname[39] = "C2H3CHO";
  kname[40] = "C4H2";
  kname[41] = "iC4H3";
  kname[42] = "C4H4";
  kname[43] = "iC4H5";
  kname[44] = "C4H5-2";
  kname[45] = "C4H6";
  kname[46] = "C4H612";
  kname[47] = "C4H6-2";
  kname[48] = "C4H7";
  kname[49] = "C4H81";
  kname[50] = "pC4H9";
  kname[51] = "NC12H26";
  kname[52] = "C6H12";
  kname[53] = "C6H11";
  kname[54] = "C5H10";
  kname[55] = "N2";
}

// compute the sparsity pattern of the chemistry Jacobian
void
SPARSITY_INFO(int* nJdata, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 3249> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 56> conc = {0.0};
  for (int n = 0; n < 56; n++) {
    conc[n] = 1.0 / 56.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 57; k++) {
    for (int l = 0; l < 57; l++) {
      if (Jac[57 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 3249> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 56> conc = {0.0};
  for (int n = 0; n < 56; n++) {
    conc[n] = 1.0 / 56.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 57; k++) {
    for (int l = 0; l < 57; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[57 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 3249> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 56> conc = {0.0};
  for (int n = 0; n < 56; n++) {
    conc[n] = 1.0 / 56.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 57; k++) {
    for (int l = 0; l < 57; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[57 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 3249> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 56> conc = {0.0};
  for (int n = 0; n < 56; n++) {
    conc[n] = 1.0 / 56.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int nc = 0; nc < NCELLS; nc++) {
    int offset_row = nc * 57;
    int offset_col = nc * 57;
    for (int k = 0; k < 57; k++) {
      for (int l = 0; l < 57; l++) {
        if (Jac[57 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 3249> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 56> conc = {0.0};
  for (int n = 0; n < 56; n++) {
    conc[n] = 1.0 / 56.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtrs[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 57;
      for (int l = 0; l < 57; l++) {
        for (int k = 0; k < 57; k++) {
          if (Jac[57 * k + l] != 0.0) {
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
      int offset = nc * 57;
      for (int l = 0; l < 57; l++) {
        for (int k = 0; k < 57; k++) {
          if (Jac[57 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 3249> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 56> conc = {0.0};
  for (int n = 0; n < 56; n++) {
    conc[n] = 1.0 / 56.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 57;
      for (int l = 0; l < 57; l++) {
        for (int k = 0; k < 57; k++) {
          if (k == l) {
            colVals[nJdata_tmp - 1] = l + 1 + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[57 * k + l] != 0.0) {
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
      int offset = nc * 57;
      for (int l = 0; l < 57; l++) {
        for (int k = 0; k < 57; k++) {
          if (k == l) {
            colVals[nJdata_tmp] = l + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[57 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 3249> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 56> conc = {0.0};
  for (int n = 0; n < 56; n++) {
    conc[n] = 1.0 / 56.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int k = 0; k < 57; k++) {
    for (int l = 0; l < 57; l++) {
      if (k == l) {
        rowVals[nJdata_tmp] = l;
        indx[nJdata_tmp] = 57 * k + l;
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[57 * k + l] != 0.0) {
          rowVals[nJdata_tmp] = l;
          indx[nJdata_tmp] = 57 * k + l;
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
  amrex::GpuArray<amrex::Real, 3249> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 56> conc = {0.0};
  for (int n = 0; n < 56; n++) {
    conc[n] = 1.0 / 56.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int l = 0; l < 57; l++) {
      for (int k = 0; k < 57; k++) {
        if (k == l) {
          colVals[nJdata_tmp - 1] = l + 1;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[57 * k + l] != 0.0) {
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
    for (int l = 0; l < 57; l++) {
      for (int k = 0; k < 57; k++) {
        if (k == l) {
          colVals[nJdata_tmp] = l;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[57 * k + l] != 0.0) {
            colVals[nJdata_tmp] = k;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
      }
      rowPtr[l + 1] = nJdata_tmp;
    }
  }
}
