#include "mechanism.H"
const int rmap[268] = {
  4,   17,  37,  38,  39,  58,  64,  77,  90,  95,  107, 113, 126, 134, 141,
  148, 149, 168, 172, 192, 196, 200, 204, 208, 212, 28,  18,  19,  20,  21,
  22,  23,  24,  30,  36,  0,   1,   2,   3,   5,   6,   7,   8,   9,   10,
  11,  12,  13,  14,  15,  16,  25,  26,  27,  29,  31,  32,  33,  34,  35,
  40,  41,  42,  43,  44,  45,  46,  47,  48,  49,  50,  51,  52,  53,  54,
  55,  56,  57,  59,  60,  61,  62,  63,  65,  66,  67,  68,  69,  70,  71,
  72,  73,  74,  75,  76,  78,  79,  80,  81,  82,  83,  84,  85,  86,  87,
  88,  89,  91,  92,  93,  94,  96,  97,  98,  99,  100, 101, 102, 103, 104,
  105, 106, 108, 109, 110, 111, 112, 114, 115, 116, 117, 118, 119, 120, 121,
  122, 123, 124, 125, 127, 128, 129, 130, 131, 132, 133, 135, 136, 137, 138,
  139, 140, 142, 143, 144, 145, 146, 147, 150, 151, 152, 153, 154, 155, 156,
  157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 169, 170, 171, 173,
  174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188,
  189, 190, 191, 193, 194, 195, 197, 198, 199, 201, 202, 203, 205, 206, 207,
  209, 210, 211, 213, 214, 215, 216, 217, 218, 219, 220, 221, 222, 223, 224,
  225, 226, 227, 228, 229, 230, 231, 232, 233, 234, 235, 236, 237, 238, 239,
  240, 241, 242, 243, 244, 245, 246, 247, 248, 249, 250, 251, 252, 253, 254,
  255, 256, 257, 258, 259, 260, 261, 262, 263, 264, 265, 266, 267};

// Returns 0-based map of reaction order
void
GET_RMAP(int* _rmap)
{
  for (int j = 0; j < 268; ++j) {
    _rmap[j] = rmap[j];
  }
}

// Returns a count of species in a reaction, and their indices
// and stoichiometric coefficients. (Eq 50)
void
CKINU(int* i, int* nspec, int* ki, int* nu)
{
  const int ns[268] = {
    4, 4, 4, 3, 3, 3, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 2, 2, 3, 2, 2, 2, 2, 3,
    4, 4, 4, 3, 4, 3, 4, 4, 4, 4, 4, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4,
    4, 4, 5, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 2, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 4, 4, 4, 3, 3, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 3, 4, 3, 4, 3, 3, 4, 5, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4,
    3, 3, 4, 4, 4, 4, 4, 5, 4, 3, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 5, 4, 3, 3,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 4, 4, 4, 4, 4, 5, 4, 3, 4, 5, 4, 3, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 5, 4, 3, 3, 3, 4, 4, 3, 3, 4, 4, 3,
    3, 4, 4, 3, 3, 4, 4, 3, 3, 4, 4, 3, 3, 4, 4, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3,
    3, 3, 3, 3, 3, 3, 2, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    3, 3, 3, 3, 3, 3, 2, 2, 4, 4, 4, 4, 4, 4, 3, 3, 3, 5};
  const int kiv[1340] = {
    1,  8,  2,  3,  0,  5,  2,  1,  3,  0,  5,  3,  1,  6,  0,  3,  6,  2,  0,
    0,  1,  8,  4,  0,  0,  1,  4,  3,  0,  0,  5,  8,  1,  4,  0,  4,  3,  6,
    8,  0,  1,  4,  6,  2,  0,  4,  2,  8,  3,  0,  4,  7,  8,  0,  0,  4,  7,
    8,  0,  0,  1,  7,  6,  3,  0,  1,  7,  5,  4,  0,  7,  2,  4,  3,  0,  7,
    3,  6,  4,  0,  7,  3,  6,  4,  0,  3,  7,  0,  0,  0,  1,  5,  0,  0,  0,
    1,  3,  6,  0,  0,  2,  8,  0,  0,  0,  1,  5,  0,  0,  0,  1,  5,  0,  0,
    0,  1,  5,  0,  0,  0,  1,  2,  3,  0,  0,  16, 3,  17, 1,  0,  16, 3,  17,
    1,  0,  16, 4,  17, 3,  0,  16, 2,  17, 0,  0,  16, 8,  17, 2,  0,  13, 16,
    1,  0,  0,  1,  13, 16, 5,  0,  13, 2,  16, 3,  0,  13, 2,  17, 1,  0,  13,
    3,  16, 6,  0,  13, 8,  16, 4,  0,  13, 16, 1,  0,  0,  16, 5,  14, 0,  0,
    1,  13, 14, 0,  0,  9,  1,  11, 0,  0,  9,  2,  1,  13, 0,  9,  3,  14, 1,
    0,  9,  5,  11, 1,  0,  9,  8,  13, 3,  0,  9,  8,  17, 1,  0,  9,  4,  14,
    3,  0,  9,  18, 5,  0,  0,  10, 52, 9,  52, 0,  10, 2,  16, 5,  0,  10, 2,
    1,  13, 0,  10, 3,  14, 1,  0,  10, 5,  11, 1,  0,  10, 8,  16, 1,  3,  10,
    8,  16, 6,  0,  10, 6,  9,  6,  0,  10, 16, 9,  16, 0,  10, 17, 9,  17, 0,
    10, 17, 14, 16, 0,  14, 1,  15, 0,  0,  14, 1,  5,  13, 0,  14, 2,  13, 3,
    0,  14, 3,  6,  13, 0,  14, 8,  13, 4,  0,  14, 4,  7,  13, 0,  11, 1,  12,
    0,  0,  11, 2,  14, 1,  0,  11, 3,  9,  6,  0,  11, 3,  10, 6,  0,  11, 8,
    15, 2,  0,  11, 8,  14, 3,  0,  11, 4,  12, 8,  0,  11, 4,  15, 3,  0,  11,
    7,  12, 4,  0,  11, 13, 12, 16, 0,  14, 11, 12, 13, 0,  9,  11, 20, 1,  0,
    10, 11, 20, 1,  0,  11, 22, 0,  0,  0,  11, 21, 1,  0,  0,  15, 1,  14, 5,
    0,  15, 1,  11, 3,  0,  15, 1,  10, 6,  0,  15, 2,  14, 3,  0,  15, 3,  14,
    6,  0,  15, 8,  14, 4,  0,  12, 1,  11, 5,  0,  12, 2,  11, 3,  0,  12, 3,
    11, 6,  0,  9,  12, 11, 0,  0,  10, 12, 11, 0,  0,  19, 18, 1,  0,  0,  18,
    2,  9,  16, 0,  18, 3,  11, 16, 0,  18, 13, 19, 16, 0,  18, 11, 24, 0,  0,
    19, 1,  20, 0,  0,  19, 1,  18, 5,  0,  19, 2,  11, 16, 0,  19, 3,  18, 6,
    0,  19, 8,  18, 4,  0,  19, 8,  23, 2,  0,  19, 8,  14, 13, 0,  19, 4,  23,
    3,  0,  19, 7,  20, 4,  0,  19, 13, 20, 16, 0,  19, 13, 27, 0,  0,  19, 11,
    18, 12, 0,  19, 11, 25, 0,  0,  19, 11, 1,  24, 0,  19, 18, 20, 0,  0,  23,
    11, 16, 0,  0,  23, 1,  11, 13, 0,  23, 8,  14, 16, 3,  20, 1,  21, 0,  0,
    20, 1,  19, 5,  0,  20, 2,  19, 3,  0,  20, 2,  11, 13, 0,  20, 2,  9,  14,
    0,  20, 3,  19, 6,  0,  20, 13, 21, 16, 0,  20, 9,  1,  24, 0,  20, 10, 1,
    24, 0,  20, 11, 19, 12, 0,  26, 20, 11, 0,  0,  20, 8,  19, 4,  0,  19, 20,
    28, 0,  0,  21, 1,  22, 0,  0,  21, 1,  20, 5,  0,  21, 2,  14, 11, 0,  21,
    8,  20, 4,  0,  21, 4,  22, 8,  0,  21, 4,  20, 7,  0,  21, 4,  14, 11, 3,
    21, 7,  22, 4,  0,  19, 21, 29, 0,  0,  19, 21, 11, 24, 0,  22, 1,  21, 5,
    0,  22, 2,  21, 3,  0,  22, 3,  21, 6,  0,  22, 10, 21, 11, 0,  22, 11, 21,
    12, 0,  1,  24, 25, 0,  0,  2,  24, 27, 1,  0,  3,  24, 27, 1,  0,  8,  24,
    27, 3,  0,  4,  24, 25, 8,  0,  4,  24, 19, 14, 3,  13, 24, 25, 16, 0,  11,
    24, 29, 0,  0,  25, 1,  26, 0,  0,  25, 1,  20, 11, 0,  25, 1,  5,  24, 0,
    25, 2,  27, 1,  0,  25, 2,  21, 13, 0,  25, 2,  3,  24, 0,  25, 3,  6,  24,
    0,  25, 4,  7,  24, 0,  25, 11, 12, 24, 0,  27, 1,  20, 13, 0,  27, 2,  19,
    16, 3,  27, 3,  19, 16, 6,  1,  26, 21, 11, 0,  1,  26, 25, 5,  0,  2,  26,
    21, 14, 0,  3,  26, 25, 6,  0,  8,  26, 25, 4,  0,  4,  26, 21, 14, 3,  11,
    26, 25, 12, 0,  28, 1,  29, 0,  0,  28, 1,  11, 24, 0,  28, 4,  14, 3,  24,
    28, 13, 29, 16, 0,  29, 1,  30, 0,  0,  29, 1,  20, 21, 0,  29, 1,  25, 11,
    0,  29, 1,  28, 5,  0,  29, 2,  13, 26, 0,  29, 2,  28, 3,  0,  29, 2,  28,
    3,  0,  29, 3,  28, 6,  0,  29, 8,  28, 4,  0,  29, 4,  28, 7,  0,  29, 11,
    28, 12, 0,  1,  30, 21, 0,  0,  1,  30, 29, 5,  0,  2,  30, 14, 26, 0,  3,
    30, 29, 6,  0,  8,  30, 29, 4,  0,  4,  30, 14, 3,  26, 11, 30, 29, 12, 0,
    31, 20, 24, 0,  0,  31, 19, 25, 0,  0,  32, 1,  33, 0,  0,  32, 1,  20, 26,
    0,  32, 1,  21, 25, 0,  20, 26, 33, 0,  0,  34, 1,  35, 0,  0,  34, 1,  20,
    30, 0,  34, 1,  25, 26, 0,  20, 30, 35, 0,  0,  36, 1,  37, 0,  0,  36, 1,
    20, 33, 0,  36, 1,  25, 30, 0,  20, 33, 37, 0,  0,  38, 1,  39, 0,  0,  38,
    1,  20, 35, 0,  38, 1,  25, 33, 0,  20, 35, 39, 0,  0,  40, 1,  41, 0,  0,
    40, 1,  20, 37, 0,  40, 1,  25, 35, 0,  20, 37, 41, 0,  0,  42, 1,  43, 0,
    0,  42, 1,  20, 39, 0,  42, 1,  25, 37, 0,  20, 39, 43, 0,  0,  47, 31, 37,
    0,  0,  20, 43, 44, 0,  0,  44, 46, 0,  0,  0,  25, 41, 45, 0,  0,  29, 39,
    45, 0,  0,  32, 37, 46, 0,  0,  42, 21, 46, 0,  0,  34, 35, 46, 0,  0,  40,
    26, 46, 0,  0,  36, 33, 46, 0,  0,  38, 30, 46, 0,  0,  21, 43, 0,  0,  0,
    41, 26, 0,  0,  0,  39, 30, 0,  0,  0,  33, 37, 0,  0,  0,  35, 0,  0,  0,
    0,  1,  0,  5,  44, 0,  1,  0,  5,  45, 0,  1,  0,  5,  46, 0,  0,  2,  3,
    44, 0,  0,  2,  3,  45, 0,  0,  2,  3,  46, 0,  0,  3,  6,  44, 0,  0,  3,
    6,  45, 0,  0,  3,  6,  46, 0,  0,  8,  4,  44, 0,  0,  8,  4,  45, 0,  0,
    8,  4,  46, 0,  4,  0,  7,  44, 0,  4,  0,  7,  45, 0,  4,  0,  7,  46, 0,
    11, 0,  12, 44, 0,  11, 0,  12, 45, 0,  11, 0,  12, 46, 0,  8,  44, 48, 0,
    0,  48, 8,  44, 0,  0,  8,  45, 48, 0,  0,  48, 8,  45, 0,  0,  8,  46, 48,
    0,  0,  48, 8,  46, 0,  0,  48, 49, 0,  0,  0,  49, 48, 0,  0,  0,  8,  44,
    47, 4,  0,  47, 4,  8,  44, 0,  8,  45, 47, 4,  0,  47, 4,  8,  45, 0,  8,
    46, 47, 4,  0,  47, 4,  8,  46, 0,  49, 8,  50, 0,  0,  50, 49, 8,  0,  0,
    50, 51, 3,  0,  0,  51, 20, 21, 23, 3};
  const int nuv[1340] = {
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -2, 1,  1, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 2, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -2, 1,  1, 0, 0, -2, 1,  1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -2, 1,  0, 0, 0, -2, 1,  0, 0, 0, -1, -1, 1, 0, 0,
    -2, 1,  0, 0, 0, -2, 1,  0, 0, 0, -2, 1,  0, 0, 0, -2, 1,  0, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, 1,  1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 2, 0, -1, -1, 1, 1, 0, -2, 1,  1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 1, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -2, 1,  0, 0, 0, -2, 1,  1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 2, 0, 0, -1, -1, 2, 0, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 1, 0, -2, 1,  1, 0, 0, -1, 1,  1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 1, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, 1,  1, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 1, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 2, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 1, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 2, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 1,
    -1, -1, 1, 1, 1, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 1, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 1, -1, -1, 1, 1, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 2, 0, 0,
    -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0,
    -1, -1, 1, 1, 1, -1, -1, 1, 1, 0, -1, 1,  1, 0, 0, -1, 1,  1, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0,
    -1, -1, 1, 0, 0, -1, -1, 1, 1, 0, -1, -1, 1, 1, 0, -1, -1, 1, 0, 0,
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
  if (*i < 1) {
    // Return max num species per reaction
    *nspec = 5;
  } else {
    if (*i > 268) {
      *nspec = -1;
    } else {
      *nspec = ns[*i - 1];
      for (int j = 0; j < *nspec; ++j) {
        ki[j] = kiv[(*i - 1) * 5 + j] + 1;
        nu[j] = nuv[(*i - 1) * 5 + j];
      }
    }
  }
}

// Returns the progress rates of each reactions
// Given P, T, and mole fractions
void
CKKFKR(
  amrex::Real* P,
  amrex::Real* T,
  amrex::Real* x,
  amrex::Real* q_f,
  amrex::Real* q_r)
{
  int id;            // loop counter
  amrex::Real c[53]; // temporary storage
  amrex::Real PORT =
    1e6 * (*P) /
    (8.31446261815324e+07 * (*T)); // 1e6 * P/RT so c goes to SI units

  // Compute conversion, see Eq 10
  for (id = 0; id < 53; ++id) {
    c[id] = x[id] * PORT;
  }

  // convert to chemkin units
  progressRateFR(q_f, q_r, c, *T);

  // convert to chemkin units
  for (id = 0; id < 268; ++id) {
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
  amrex::Real g_RT[53];
  gibbs(g_RT, tc);

  amrex::Real sc_qss[1];
  amrex::Real kf_qss[0], qf_qss[0], qr_qss[0];
  comp_qfqr(q_f, q_r, sc, sc_qss, tc, invT);

  return;
}

// save atomic weights into array
void
atomicWeight(amrex::Real* awt)
{
  awt[0] = 15.999400; // O
  awt[1] = 1.007970;  // H
  awt[2] = 12.011150; // C
  awt[3] = 14.006700; // N
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
  int kd = 4;
  // Zero ncf
  for (id = 0; id < kd * 53; ++id) {
    ncf[id] = 0;
  }

  // NC12H26
  ncf[0 * kd + 2] = 12; // C
  ncf[0 * kd + 1] = 26; // H

  // H
  ncf[1 * kd + 1] = 1; // H

  // O
  ncf[2 * kd + 0] = 1; // O

  // OH
  ncf[3 * kd + 1] = 1; // H
  ncf[3 * kd + 0] = 1; // O

  // HO2
  ncf[4 * kd + 1] = 1; // H
  ncf[4 * kd + 0] = 2; // O

  // H2
  ncf[5 * kd + 1] = 2; // H

  // H2O
  ncf[6 * kd + 1] = 2; // H
  ncf[6 * kd + 0] = 1; // O

  // H2O2
  ncf[7 * kd + 1] = 2; // H
  ncf[7 * kd + 0] = 2; // O

  // O2
  ncf[8 * kd + 0] = 2; // O

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

  // CO
  ncf[16 * kd + 2] = 1; // C
  ncf[16 * kd + 0] = 1; // O

  // CO2
  ncf[17 * kd + 2] = 1; // C
  ncf[17 * kd + 0] = 2; // O

  // C2H2
  ncf[18 * kd + 2] = 2; // C
  ncf[18 * kd + 1] = 2; // H

  // C2H3
  ncf[19 * kd + 2] = 2; // C
  ncf[19 * kd + 1] = 3; // H

  // C2H4
  ncf[20 * kd + 2] = 2; // C
  ncf[20 * kd + 1] = 4; // H

  // C2H5
  ncf[21 * kd + 2] = 2; // C
  ncf[21 * kd + 1] = 5; // H

  // C2H6
  ncf[22 * kd + 2] = 2; // C
  ncf[22 * kd + 1] = 6; // H

  // CH2CHO
  ncf[23 * kd + 2] = 2; // C
  ncf[23 * kd + 1] = 3; // H
  ncf[23 * kd + 0] = 1; // O

  // aC3H5
  ncf[24 * kd + 2] = 3; // C
  ncf[24 * kd + 1] = 5; // H

  // C3H6
  ncf[25 * kd + 2] = 3; // C
  ncf[25 * kd + 1] = 6; // H

  // nC3H7
  ncf[26 * kd + 2] = 3; // C
  ncf[26 * kd + 1] = 7; // H

  // C2H3CHO
  ncf[27 * kd + 2] = 3; // C
  ncf[27 * kd + 1] = 4; // H
  ncf[27 * kd + 0] = 1; // O

  // C4H7
  ncf[28 * kd + 2] = 4; // C
  ncf[28 * kd + 1] = 7; // H

  // C4H81
  ncf[29 * kd + 2] = 4; // C
  ncf[29 * kd + 1] = 8; // H

  // pC4H9
  ncf[30 * kd + 2] = 4; // C
  ncf[30 * kd + 1] = 9; // H

  // C5H9
  ncf[31 * kd + 2] = 5; // C
  ncf[31 * kd + 1] = 9; // H

  // C5H10
  ncf[32 * kd + 2] = 5;  // C
  ncf[32 * kd + 1] = 10; // H

  // PXC5H11
  ncf[33 * kd + 2] = 5;  // C
  ncf[33 * kd + 1] = 11; // H

  // C6H12
  ncf[34 * kd + 2] = 6;  // C
  ncf[34 * kd + 1] = 12; // H

  // PXC6H13
  ncf[35 * kd + 2] = 6;  // C
  ncf[35 * kd + 1] = 13; // H

  // C7H14
  ncf[36 * kd + 2] = 7;  // C
  ncf[36 * kd + 1] = 14; // H

  // PXC7H15
  ncf[37 * kd + 2] = 7;  // C
  ncf[37 * kd + 1] = 15; // H

  // C8H16
  ncf[38 * kd + 2] = 8;  // C
  ncf[38 * kd + 1] = 16; // H

  // PXC8H17
  ncf[39 * kd + 2] = 8;  // C
  ncf[39 * kd + 1] = 17; // H

  // C9H18
  ncf[40 * kd + 2] = 9;  // C
  ncf[40 * kd + 1] = 18; // H

  // PXC9H19
  ncf[41 * kd + 2] = 9;  // C
  ncf[41 * kd + 1] = 19; // H

  // C10H20
  ncf[42 * kd + 2] = 10; // C
  ncf[42 * kd + 1] = 20; // H

  // PXC10H21
  ncf[43 * kd + 2] = 10; // C
  ncf[43 * kd + 1] = 21; // H

  // PXC12H25
  ncf[44 * kd + 2] = 12; // C
  ncf[44 * kd + 1] = 25; // H

  // SXC12H25
  ncf[45 * kd + 2] = 12; // C
  ncf[45 * kd + 1] = 25; // H

  // S3XC12H25
  ncf[46 * kd + 2] = 12; // C
  ncf[46 * kd + 1] = 25; // H

  // C12H24
  ncf[47 * kd + 2] = 12; // C
  ncf[47 * kd + 1] = 24; // H

  // C12H25O2
  ncf[48 * kd + 2] = 12; // C
  ncf[48 * kd + 1] = 25; // H
  ncf[48 * kd + 0] = 2;  // O

  // C12OOH
  ncf[49 * kd + 2] = 12; // C
  ncf[49 * kd + 1] = 25; // H
  ncf[49 * kd + 0] = 2;  // O

  // O2C12H24OOH
  ncf[50 * kd + 2] = 12; // C
  ncf[50 * kd + 1] = 25; // H
  ncf[50 * kd + 0] = 4;  // O

  // OC12H23OOH
  ncf[51 * kd + 2] = 12; // C
  ncf[51 * kd + 1] = 24; // H
  ncf[51 * kd + 0] = 3;  // O

  // N2
  ncf[52 * kd + 3] = 2; // N
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
  kname.resize(53);
  kname[0] = "NC12H26";
  kname[1] = "H";
  kname[2] = "O";
  kname[3] = "OH";
  kname[4] = "HO2";
  kname[5] = "H2";
  kname[6] = "H2O";
  kname[7] = "H2O2";
  kname[8] = "O2";
  kname[9] = "CH2";
  kname[10] = "CH2*";
  kname[11] = "CH3";
  kname[12] = "CH4";
  kname[13] = "HCO";
  kname[14] = "CH2O";
  kname[15] = "CH3O";
  kname[16] = "CO";
  kname[17] = "CO2";
  kname[18] = "C2H2";
  kname[19] = "C2H3";
  kname[20] = "C2H4";
  kname[21] = "C2H5";
  kname[22] = "C2H6";
  kname[23] = "CH2CHO";
  kname[24] = "aC3H5";
  kname[25] = "C3H6";
  kname[26] = "nC3H7";
  kname[27] = "C2H3CHO";
  kname[28] = "C4H7";
  kname[29] = "C4H81";
  kname[30] = "pC4H9";
  kname[31] = "C5H9";
  kname[32] = "C5H10";
  kname[33] = "PXC5H11";
  kname[34] = "C6H12";
  kname[35] = "PXC6H13";
  kname[36] = "C7H14";
  kname[37] = "PXC7H15";
  kname[38] = "C8H16";
  kname[39] = "PXC8H17";
  kname[40] = "C9H18";
  kname[41] = "PXC9H19";
  kname[42] = "C10H20";
  kname[43] = "PXC10H21";
  kname[44] = "PXC12H25";
  kname[45] = "SXC12H25";
  kname[46] = "S3XC12H25";
  kname[47] = "C12H24";
  kname[48] = "C12H25O2";
  kname[49] = "C12OOH";
  kname[50] = "O2C12H24OOH";
  kname[51] = "OC12H23OOH";
  kname[52] = "N2";
}

// compute the sparsity pattern of the chemistry Jacobian
void
SPARSITY_INFO(int* nJdata, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 2916> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 53> conc = {0.0};
  for (int n = 0; n < 53; n++) {
    conc[n] = 1.0 / 53.000000;
  }
  aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 54; k++) {
    for (int l = 0; l < 54; l++) {
      if (Jac[54 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 2916> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 53> conc = {0.0};
  for (int n = 0; n < 53; n++) {
    conc[n] = 1.0 / 53.000000;
  }
  aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 54; k++) {
    for (int l = 0; l < 54; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[54 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 2916> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 53> conc = {0.0};
  for (int n = 0; n < 53; n++) {
    conc[n] = 1.0 / 53.000000;
  }
  aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 54; k++) {
    for (int l = 0; l < 54; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[54 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 2916> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 53> conc = {0.0};
  for (int n = 0; n < 53; n++) {
    conc[n] = 1.0 / 53.000000;
  }
  aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int nc = 0; nc < NCELLS; nc++) {
    int offset_row = nc * 54;
    int offset_col = nc * 54;
    for (int k = 0; k < 54; k++) {
      for (int l = 0; l < 54; l++) {
        if (Jac[54 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 2916> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 53> conc = {0.0};
  for (int n = 0; n < 53; n++) {
    conc[n] = 1.0 / 53.000000;
  }
  aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

  if (base == 1) {
    rowPtrs[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 54;
      for (int l = 0; l < 54; l++) {
        for (int k = 0; k < 54; k++) {
          if (Jac[54 * k + l] != 0.0) {
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
      int offset = nc * 54;
      for (int l = 0; l < 54; l++) {
        for (int k = 0; k < 54; k++) {
          if (Jac[54 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 2916> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 53> conc = {0.0};
  for (int n = 0; n < 53; n++) {
    conc[n] = 1.0 / 53.000000;
  }
  aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 54;
      for (int l = 0; l < 54; l++) {
        for (int k = 0; k < 54; k++) {
          if (k == l) {
            colVals[nJdata_tmp - 1] = l + 1 + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[54 * k + l] != 0.0) {
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
      int offset = nc * 54;
      for (int l = 0; l < 54; l++) {
        for (int k = 0; k < 54; k++) {
          if (k == l) {
            colVals[nJdata_tmp] = l + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[54 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 2916> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 53> conc = {0.0};
  for (int n = 0; n < 53; n++) {
    conc[n] = 1.0 / 53.000000;
  }
  aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int k = 0; k < 54; k++) {
    for (int l = 0; l < 54; l++) {
      if (k == l) {
        rowVals[nJdata_tmp] = l;
        indx[nJdata_tmp] = 54 * k + l;
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[54 * k + l] != 0.0) {
          rowVals[nJdata_tmp] = l;
          indx[nJdata_tmp] = 54 * k + l;
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
  amrex::GpuArray<amrex::Real, 2916> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 53> conc = {0.0};
  for (int n = 0; n < 53; n++) {
    conc[n] = 1.0 / 53.000000;
  }
  aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int l = 0; l < 54; l++) {
      for (int k = 0; k < 54; k++) {
        if (k == l) {
          colVals[nJdata_tmp - 1] = l + 1;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[54 * k + l] != 0.0) {
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
    for (int l = 0; l < 54; l++) {
      for (int k = 0; k < 54; k++) {
        if (k == l) {
          colVals[nJdata_tmp] = l;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[54 * k + l] != 0.0) {
            colVals[nJdata_tmp] = k;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
      }
      rowPtr[l + 1] = nJdata_tmp;
    }
  }
}
