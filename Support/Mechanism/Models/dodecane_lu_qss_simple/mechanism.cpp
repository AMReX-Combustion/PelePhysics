#include "mechanism.H"
const int rmap[268] = {
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
  const int kiv[1340] = {
    1,  8,  4,  0,  0,  3,  7,  0,  0,  0,  15, 5,  14, 0,  0,  1,  13, 14, 0,
    0,  9,  1,  11, 0,  0,  14, 1,  50, 0,  0,  11, 1,  12, 0,  0,  11, 21, 0,
    0,  0,  18, 17, 1,  0,  0,  18, 1,  19, 0,  0,  18, 11, 24, 0,  0,  19, 1,
    20, 0,  0,  20, 1,  21, 0,  0,  18, 20, 28, 0,  0,  1,  23, 24, 0,  0,  11,
    23, 28, 0,  0,  24, 1,  25, 0,  0,  27, 1,  28, 0,  0,  28, 1,  29, 0,  0,
    31, 1,  32, 0,  0,  33, 1,  34, 0,  0,  35, 1,  36, 0,  0,  37, 1,  38, 0,
    0,  39, 1,  40, 0,  0,  41, 1,  42, 0,  0,  15, 2,  16, 0,  0,  1,  5,  0,
    0,  0,  1,  3,  6,  0,  0,  2,  8,  0,  0,  0,  1,  5,  0,  0,  0,  1,  5,
    0,  0,  0,  1,  5,  0,  0,  0,  1,  2,  3,  0,  0,  13, 15, 1,  0,  0,  13,
    15, 1,  0,  0,  1,  8,  2,  3,  0,  5,  2,  1,  3,  0,  5,  3,  1,  6,  0,
    3,  6,  2,  0,  0,  1,  4,  3,  0,  0,  5,  8,  1,  4,  0,  4,  3,  6,  8,
    0,  1,  4,  6,  2,  0,  4,  2,  8,  3,  0,  4,  7,  8,  0,  0,  4,  7,  8,
    0,  0,  1,  7,  6,  3,  0,  1,  7,  5,  4,  0,  7,  2,  4,  3,  0,  7,  3,
    6,  4,  0,  7,  3,  6,  4,  0,  15, 3,  16, 1,  0,  15, 3,  16, 1,  0,  15,
    4,  16, 3,  0,  15, 8,  16, 2,  0,  1,  13, 15, 5,  0,  13, 2,  15, 3,  0,
    13, 2,  16, 1,  0,  13, 3,  15, 6,  0,  13, 8,  15, 4,  0,  9,  2,  1,  13,
    0,  9,  3,  14, 1,  0,  9,  5,  11, 1,  0,  9,  8,  13, 3,  0,  9,  8,  16,
    1,  0,  9,  4,  14, 3,  0,  9,  17, 5,  0,  0,  10, 49, 9,  49, 0,  10, 2,
    15, 5,  0,  10, 2,  1,  13, 0,  10, 3,  14, 1,  0,  10, 5,  11, 1,  0,  10,
    8,  15, 1,  3,  10, 8,  15, 6,  0,  10, 6,  9,  6,  0,  10, 15, 9,  15, 0,
    10, 16, 9,  16, 0,  10, 16, 14, 15, 0,  14, 1,  5,  13, 0,  14, 2,  13, 3,
    0,  14, 3,  6,  13, 0,  14, 8,  13, 4,  0,  14, 4,  7,  13, 0,  11, 2,  14,
    1,  0,  11, 3,  9,  6,  0,  11, 3,  10, 6,  0,  11, 8,  50, 2,  0,  11, 8,
    14, 3,  0,  11, 4,  12, 8,  0,  11, 4,  50, 3,  0,  11, 7,  12, 4,  0,  11,
    13, 12, 15, 0,  14, 11, 12, 13, 0,  9,  11, 19, 1,  0,  10, 11, 19, 1,  0,
    11, 20, 1,  0,  0,  50, 1,  14, 5,  0,  50, 1,  11, 3,  0,  50, 1,  10, 6,
    0,  50, 2,  14, 3,  0,  50, 3,  14, 6,  0,  50, 8,  14, 4,  0,  12, 1,  11,
    5,  0,  12, 2,  11, 3,  0,  12, 3,  11, 6,  0,  9,  12, 11, 0,  0,  10, 12,
    11, 0,  0,  17, 2,  9,  15, 0,  17, 3,  11, 15, 0,  17, 13, 18, 15, 0,  17,
    11, 23, 0,  0,  18, 1,  17, 5,  0,  18, 2,  11, 15, 0,  18, 3,  17, 6,  0,
    18, 8,  17, 4,  0,  18, 8,  22, 2,  0,  18, 8,  14, 13, 0,  18, 4,  22, 3,
    0,  18, 7,  19, 4,  0,  18, 13, 19, 15, 0,  18, 13, 26, 0,  0,  18, 11, 17,
    12, 0,  18, 11, 1,  23, 0,  18, 17, 19, 0,  0,  22, 11, 15, 0,  0,  22, 1,
    11, 13, 0,  22, 8,  14, 15, 3,  19, 1,  18, 5,  0,  19, 2,  18, 3,  0,  19,
    2,  11, 13, 0,  19, 2,  9,  14, 0,  19, 3,  18, 6,  0,  19, 13, 20, 15, 0,
    19, 9,  1,  23, 0,  19, 10, 1,  23, 0,  19, 11, 18, 12, 0,  25, 19, 11, 0,
    0,  19, 8,  18, 4,  0,  18, 19, 27, 0,  0,  20, 1,  19, 5,  0,  20, 2,  14,
    11, 0,  20, 8,  19, 4,  0,  20, 4,  21, 8,  0,  20, 4,  19, 7,  0,  20, 4,
    14, 11, 3,  20, 7,  21, 4,  0,  18, 20, 11, 23, 0,  21, 1,  20, 5,  0,  21,
    2,  20, 3,  0,  21, 3,  20, 6,  0,  21, 10, 20, 11, 0,  21, 11, 20, 12, 0,
    2,  23, 26, 1,  0,  3,  23, 26, 1,  0,  8,  23, 26, 3,  0,  4,  23, 24, 8,
    0,  4,  23, 18, 14, 3,  13, 23, 24, 15, 0,  24, 1,  19, 11, 0,  24, 1,  5,
    23, 0,  24, 2,  26, 1,  0,  24, 2,  20, 13, 0,  24, 2,  3,  23, 0,  24, 3,
    6,  23, 0,  24, 4,  7,  23, 0,  24, 11, 12, 23, 0,  26, 1,  19, 13, 0,  26,
    2,  18, 15, 3,  26, 3,  18, 15, 6,  1,  25, 20, 11, 0,  1,  25, 24, 5,  0,
    2,  25, 20, 14, 0,  3,  25, 24, 6,  0,  8,  25, 24, 4,  0,  4,  25, 20, 14,
    3,  11, 25, 24, 12, 0,  27, 1,  11, 23, 0,  27, 4,  14, 3,  23, 27, 13, 28,
    15, 0,  28, 1,  19, 20, 0,  28, 1,  24, 11, 0,  28, 1,  27, 5,  0,  28, 2,
    13, 25, 0,  28, 2,  27, 3,  0,  28, 2,  27, 3,  0,  28, 3,  27, 6,  0,  28,
    8,  27, 4,  0,  28, 4,  27, 7,  0,  28, 11, 27, 12, 0,  1,  29, 20, 0,  0,
    1,  29, 28, 5,  0,  2,  29, 14, 25, 0,  3,  29, 28, 6,  0,  8,  29, 28, 4,
    0,  4,  29, 14, 3,  25, 11, 29, 28, 12, 0,  30, 19, 23, 0,  0,  30, 18, 24,
    0,  0,  31, 1,  19, 25, 0,  31, 1,  20, 24, 0,  19, 25, 32, 0,  0,  33, 1,
    19, 29, 0,  33, 1,  24, 25, 0,  19, 29, 34, 0,  0,  35, 1,  19, 32, 0,  35,
    1,  24, 29, 0,  19, 32, 36, 0,  0,  37, 1,  19, 34, 0,  37, 1,  24, 32, 0,
    19, 34, 38, 0,  0,  39, 1,  19, 36, 0,  39, 1,  24, 34, 0,  19, 36, 40, 0,
    0,  41, 1,  19, 38, 0,  41, 1,  24, 36, 0,  19, 38, 42, 0,  0,  46, 30, 36,
    0,  0,  19, 42, 43, 0,  0,  43, 45, 0,  0,  0,  24, 40, 44, 0,  0,  28, 38,
    44, 0,  0,  31, 36, 45, 0,  0,  41, 20, 45, 0,  0,  33, 34, 45, 0,  0,  39,
    25, 45, 0,  0,  35, 32, 45, 0,  0,  37, 29, 45, 0,  0,  20, 42, 0,  0,  0,
    40, 25, 0,  0,  0,  38, 29, 0,  0,  0,  32, 36, 0,  0,  0,  34, 0,  0,  0,
    0,  1,  0,  5,  43, 0,  1,  0,  5,  44, 0,  1,  0,  5,  45, 0,  0,  2,  3,
    43, 0,  0,  2,  3,  44, 0,  0,  2,  3,  45, 0,  0,  3,  6,  43, 0,  0,  3,
    6,  44, 0,  0,  3,  6,  45, 0,  0,  8,  4,  43, 0,  0,  8,  4,  44, 0,  0,
    8,  4,  45, 0,  4,  0,  7,  43, 0,  4,  0,  7,  44, 0,  4,  0,  7,  45, 0,
    11, 0,  12, 43, 0,  11, 0,  12, 44, 0,  11, 0,  12, 45, 0,  8,  43, 47, 0,
    0,  47, 8,  43, 0,  0,  8,  44, 47, 0,  0,  47, 8,  44, 0,  0,  8,  45, 47,
    0,  0,  47, 8,  45, 0,  0,  47, 51, 0,  0,  0,  51, 47, 0,  0,  0,  8,  43,
    46, 4,  0,  46, 4,  8,  43, 0,  8,  44, 46, 4,  0,  46, 4,  8,  44, 0,  8,
    45, 46, 4,  0,  46, 4,  8,  45, 0,  51, 8,  52, 0,  0,  52, 51, 8,  0,  0,
    52, 48, 3,  0,  0,  48, 19, 20, 22, 3};
  const int nuv[1340] = {
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
  amrex::Real c[50]; // temporary storage
  amrex::Real PORT =
    1e6 * (*P) /
    (8.31446261815324e+07 * (*T)); // 1e6 * P/RT so c goes to SI units

  // Compute conversion, see Eq 10
  for (id = 0; id < 50; ++id) {
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
  amrex::Real g_RT[50];
  gibbs(g_RT, tc);
  amrex::Real g_RT_qss[3];
  gibbs_qss(g_RT_qss, tc);

  amrex::Real sc_qss[3];
  amrex::Real kf_qss[14], qf_qss[14], qr_qss[14];
  // Fill sc_qss here
  comp_k_f_qss(tc, invT, kf_qss);
  comp_qss_coeff(kf_qss, qf_qss, qr_qss, sc, tc, g_RT, g_RT_qss);
  comp_sc_qss(sc_qss, qf_qss, qr_qss);
  comp_qfqr(q_f, q_r, sc, sc_qss, tc, invT);

  return;
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
  int id; // loop counter
  int kd = 4;
  // Zero ncf
  for (id = 0; id < kd * 50; ++id) {
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

  // CH2
  ncf[9 * kd + 0] = 1; // C
  ncf[9 * kd + 1] = 2; // H

  // CH2*
  ncf[10 * kd + 0] = 1; // C
  ncf[10 * kd + 1] = 2; // H

  // CH3
  ncf[11 * kd + 0] = 1; // C
  ncf[11 * kd + 1] = 3; // H

  // CH4
  ncf[12 * kd + 0] = 1; // C
  ncf[12 * kd + 1] = 4; // H

  // HCO
  ncf[13 * kd + 0] = 1; // C
  ncf[13 * kd + 1] = 1; // H
  ncf[13 * kd + 2] = 1; // O

  // CH2O
  ncf[14 * kd + 0] = 1; // C
  ncf[14 * kd + 1] = 2; // H
  ncf[14 * kd + 2] = 1; // O

  // CO
  ncf[15 * kd + 0] = 1; // C
  ncf[15 * kd + 2] = 1; // O

  // CO2
  ncf[16 * kd + 0] = 1; // C
  ncf[16 * kd + 2] = 2; // O

  // C2H2
  ncf[17 * kd + 0] = 2; // C
  ncf[17 * kd + 1] = 2; // H

  // C2H3
  ncf[18 * kd + 0] = 2; // C
  ncf[18 * kd + 1] = 3; // H

  // C2H4
  ncf[19 * kd + 0] = 2; // C
  ncf[19 * kd + 1] = 4; // H

  // C2H5
  ncf[20 * kd + 0] = 2; // C
  ncf[20 * kd + 1] = 5; // H

  // C2H6
  ncf[21 * kd + 0] = 2; // C
  ncf[21 * kd + 1] = 6; // H

  // CH2CHO
  ncf[22 * kd + 0] = 2; // C
  ncf[22 * kd + 1] = 3; // H
  ncf[22 * kd + 2] = 1; // O

  // aC3H5
  ncf[23 * kd + 0] = 3; // C
  ncf[23 * kd + 1] = 5; // H

  // C3H6
  ncf[24 * kd + 0] = 3; // C
  ncf[24 * kd + 1] = 6; // H

  // nC3H7
  ncf[25 * kd + 0] = 3; // C
  ncf[25 * kd + 1] = 7; // H

  // C2H3CHO
  ncf[26 * kd + 0] = 3; // C
  ncf[26 * kd + 1] = 4; // H
  ncf[26 * kd + 2] = 1; // O

  // C4H7
  ncf[27 * kd + 0] = 4; // C
  ncf[27 * kd + 1] = 7; // H

  // C4H81
  ncf[28 * kd + 0] = 4; // C
  ncf[28 * kd + 1] = 8; // H

  // pC4H9
  ncf[29 * kd + 0] = 4; // C
  ncf[29 * kd + 1] = 9; // H

  // C5H9
  ncf[30 * kd + 0] = 5; // C
  ncf[30 * kd + 1] = 9; // H

  // C5H10
  ncf[31 * kd + 0] = 5;  // C
  ncf[31 * kd + 1] = 10; // H

  // PXC5H11
  ncf[32 * kd + 0] = 5;  // C
  ncf[32 * kd + 1] = 11; // H

  // C6H12
  ncf[33 * kd + 0] = 6;  // C
  ncf[33 * kd + 1] = 12; // H

  // PXC6H13
  ncf[34 * kd + 0] = 6;  // C
  ncf[34 * kd + 1] = 13; // H

  // C7H14
  ncf[35 * kd + 0] = 7;  // C
  ncf[35 * kd + 1] = 14; // H

  // PXC7H15
  ncf[36 * kd + 0] = 7;  // C
  ncf[36 * kd + 1] = 15; // H

  // C8H16
  ncf[37 * kd + 0] = 8;  // C
  ncf[37 * kd + 1] = 16; // H

  // PXC8H17
  ncf[38 * kd + 0] = 8;  // C
  ncf[38 * kd + 1] = 17; // H

  // C9H18
  ncf[39 * kd + 0] = 9;  // C
  ncf[39 * kd + 1] = 18; // H

  // PXC9H19
  ncf[40 * kd + 0] = 9;  // C
  ncf[40 * kd + 1] = 19; // H

  // C10H20
  ncf[41 * kd + 0] = 10; // C
  ncf[41 * kd + 1] = 20; // H

  // PXC10H21
  ncf[42 * kd + 0] = 10; // C
  ncf[42 * kd + 1] = 21; // H

  // PXC12H25
  ncf[43 * kd + 0] = 12; // C
  ncf[43 * kd + 1] = 25; // H

  // SXC12H25
  ncf[44 * kd + 0] = 12; // C
  ncf[44 * kd + 1] = 25; // H

  // S3XC12H25
  ncf[45 * kd + 0] = 12; // C
  ncf[45 * kd + 1] = 25; // H

  // C12H24
  ncf[46 * kd + 0] = 12; // C
  ncf[46 * kd + 1] = 24; // H

  // C12H25O2
  ncf[47 * kd + 0] = 12; // C
  ncf[47 * kd + 1] = 25; // H
  ncf[47 * kd + 2] = 2;  // O

  // OC12H23OOH
  ncf[48 * kd + 0] = 12; // C
  ncf[48 * kd + 1] = 24; // H
  ncf[48 * kd + 2] = 3;  // O

  // N2
  ncf[49 * kd + 3] = 2; // N
}

// Returns the vector of strings of element names
void
CKSYME_STR(amrex::Vector<std::string>& ename)
{
  ename.resize(4);
  ename[0] = "C";
  ename[1] = "H";
  ename[2] = "O";
  ename[3] = "N";
}

// Returns the vector of strings of species names
void
CKSYMS_STR(amrex::Vector<std::string>& kname)
{
  kname.resize(50);
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
  kname[15] = "CO";
  kname[16] = "CO2";
  kname[17] = "C2H2";
  kname[18] = "C2H3";
  kname[19] = "C2H4";
  kname[20] = "C2H5";
  kname[21] = "C2H6";
  kname[22] = "CH2CHO";
  kname[23] = "aC3H5";
  kname[24] = "C3H6";
  kname[25] = "nC3H7";
  kname[26] = "C2H3CHO";
  kname[27] = "C4H7";
  kname[28] = "C4H81";
  kname[29] = "pC4H9";
  kname[30] = "C5H9";
  kname[31] = "C5H10";
  kname[32] = "PXC5H11";
  kname[33] = "C6H12";
  kname[34] = "PXC6H13";
  kname[35] = "C7H14";
  kname[36] = "PXC7H15";
  kname[37] = "C8H16";
  kname[38] = "PXC8H17";
  kname[39] = "C9H18";
  kname[40] = "PXC9H19";
  kname[41] = "C10H20";
  kname[42] = "PXC10H21";
  kname[43] = "PXC12H25";
  kname[44] = "SXC12H25";
  kname[45] = "S3XC12H25";
  kname[46] = "C12H24";
  kname[47] = "C12H25O2";
  kname[48] = "OC12H23OOH";
  kname[49] = "N2";
}

// compute the sparsity pattern of the chemistry Jacobian
void
SPARSITY_INFO(int* nJdata, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 2601> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 50> conc = {0.0};
  for (int n = 0; n < 50; n++) {
    conc[n] = 1.0 / 50.000000;
  }
  aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 51; k++) {
    for (int l = 0; l < 51; l++) {
      if (Jac[51 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 2601> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 50> conc = {0.0};
  for (int n = 0; n < 50; n++) {
    conc[n] = 1.0 / 50.000000;
  }
  aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 51; k++) {
    for (int l = 0; l < 51; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[51 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 2601> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 50> conc = {0.0};
  for (int n = 0; n < 50; n++) {
    conc[n] = 1.0 / 50.000000;
  }
  aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 51; k++) {
    for (int l = 0; l < 51; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[51 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 2601> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 50> conc = {0.0};
  for (int n = 0; n < 50; n++) {
    conc[n] = 1.0 / 50.000000;
  }
  aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int nc = 0; nc < NCELLS; nc++) {
    int offset_row = nc * 51;
    int offset_col = nc * 51;
    for (int k = 0; k < 51; k++) {
      for (int l = 0; l < 51; l++) {
        if (Jac[51 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 2601> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 50> conc = {0.0};
  for (int n = 0; n < 50; n++) {
    conc[n] = 1.0 / 50.000000;
  }
  aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

  if (base == 1) {
    rowPtrs[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 51;
      for (int l = 0; l < 51; l++) {
        for (int k = 0; k < 51; k++) {
          if (Jac[51 * k + l] != 0.0) {
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
      int offset = nc * 51;
      for (int l = 0; l < 51; l++) {
        for (int k = 0; k < 51; k++) {
          if (Jac[51 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 2601> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 50> conc = {0.0};
  for (int n = 0; n < 50; n++) {
    conc[n] = 1.0 / 50.000000;
  }
  aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 51;
      for (int l = 0; l < 51; l++) {
        for (int k = 0; k < 51; k++) {
          if (k == l) {
            colVals[nJdata_tmp - 1] = l + 1 + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[51 * k + l] != 0.0) {
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
      int offset = nc * 51;
      for (int l = 0; l < 51; l++) {
        for (int k = 0; k < 51; k++) {
          if (k == l) {
            colVals[nJdata_tmp] = l + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[51 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 2601> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 50> conc = {0.0};
  for (int n = 0; n < 50; n++) {
    conc[n] = 1.0 / 50.000000;
  }
  aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int k = 0; k < 51; k++) {
    for (int l = 0; l < 51; l++) {
      if (k == l) {
        rowVals[nJdata_tmp] = l;
        indx[nJdata_tmp] = 51 * k + l;
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[51 * k + l] != 0.0) {
          rowVals[nJdata_tmp] = l;
          indx[nJdata_tmp] = 51 * k + l;
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
  amrex::GpuArray<amrex::Real, 2601> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 50> conc = {0.0};
  for (int n = 0; n < 50; n++) {
    conc[n] = 1.0 / 50.000000;
  }
  aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int l = 0; l < 51; l++) {
      for (int k = 0; k < 51; k++) {
        if (k == l) {
          colVals[nJdata_tmp - 1] = l + 1;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[51 * k + l] != 0.0) {
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
    for (int l = 0; l < 51; l++) {
      for (int k = 0; k < 51; k++) {
        if (k == l) {
          colVals[nJdata_tmp] = l;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[51 * k + l] != 0.0) {
            colVals[nJdata_tmp] = k;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
      }
      rowPtr[l + 1] = nJdata_tmp;
    }
  }
}
