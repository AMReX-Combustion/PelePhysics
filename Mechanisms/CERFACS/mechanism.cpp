#include "mechanism.H"
const int rmap[NUM_REACTIONS] = {
  8,   21,  112, 143, 0,   3,   5,   7,   22,  130, 1,   2,   4,   6,   9,
  10,  11,  12,  13,  14,  15,  16,  17,  18,  19,  20,  23,  24,  25,  26,
  27,  28,  29,  30,  31,  32,  33,  34,  35,  36,  37,  38,  39,  40,  41,
  42,  43,  44,  45,  46,  47,  48,  49,  50,  51,  52,  53,  54,  55,  56,
  57,  58,  59,  60,  61,  62,  63,  64,  65,  66,  67,  68,  69,  70,  71,
  72,  73,  74,  75,  76,  77,  78,  79,  80,  81,  82,  83,  84,  85,  86,
  87,  88,  89,  90,  91,  92,  93,  94,  95,  96,  97,  98,  99,  100, 101,
  102, 103, 104, 105, 106, 107, 108, 109, 110, 111, 113, 114, 115, 116, 117,
  118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129, 131, 132, 133,
  134, 135, 136, 137, 138, 139, 140, 141, 142, 144, 145, 146, 147, 148, 149};

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
    2, 4, 4, 2, 4, 3, 3, 3, 2, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 3, 3, 3, 3, 3, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 2, 3, 3, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 3, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 4, 4, 4, 4, 2, 3, 3, 4, 4, 4, 4, 4, 4, 4,
    4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 3, 3, 4, 4, 3, 4, 4, 4};
  const int kiv[NUM_GAS_REACTIONS * 4] = {
    1,  2,  0,  0,  1,  4,  2,  5,  1,  5,  2,  8,  4,  3,  0,  0,  2,  3,  4,
    5,  2,  5,  8,  0,  8,  4,  5,  0,  2,  4,  5,  0,  6,  5,  0,  0,  2,  6,
    8,  5,  2,  6,  1,  7,  6,  4,  7,  5,  6,  5,  8,  7,  6,  5,  8,  7,  2,
    7,  5,  0,  2,  7,  1,  3,  7,  4,  3,  5,  7,  5,  8,  3,  7,  5,  8,  3,
    7,  6,  3,  0,  7,  6,  3,  0,  2,  3,  7,  0,  4,  5,  7,  0,  13, 2,  16,
    0,  2,  16, 1,  15, 2,  13, 1,  16, 13, 5,  8,  16, 13, 4,  16, 5,  7,  13,
    6,  16, 13, 3,  7,  16, 16, 4,  2,  18, 16, 4,  2,  18, 16, 4,  15, 5,  16,
    4,  15, 5,  16, 5,  8,  15, 16, 3,  18, 5,  16, 3,  22, 4,  7,  16, 22, 5,
    15, 16, 2,  12, 15, 16, 20, 13, 15, 20, 16, 0,  15, 1,  0,  0,  15, 2,  0,
    0,  16, 15, 13, 0,  16, 14, 0,  0,  16, 2,  23, 0,  16, 1,  19, 0,  16, 10,
    22, 9,  16, 10, 8,  11, 16, 9,  8,  0,  16, 9,  21, 5,  2,  15, 1,  20, 15,
    4,  2,  9,  15, 5,  2,  18, 15, 5,  8,  20, 15, 3,  18, 4,  15, 3,  9,  5,
    20, 15, 2,  0,  15, 9,  2,  11, 15, 9,  0,  5,  15, 10, 11, 5,  15, 10, 18,
    9,  20, 5,  2,  9,  20, 3,  9,  4,  20, 9,  0,  4,  14, 1,  19, 0,  2,  14,
    1,  23, 14, 4,  23, 5,  14, 4,  8,  12, 14, 5,  8,  23, 14, 16, 23, 13, 14,
    10, 17, 23, 23, 2,  12, 0,  2,  23, 1,  12, 23, 4,  18, 16, 23, 4,  12, 5,
    23, 5,  8,  12, 23, 5,  19, 8,  23, 16, 12, 13, 23, 16, 19, 13, 7,  23, 6,
    12, 7,  23, 14, 3,  12, 2,  21, 0,  12, 2,  21, 0,  2,  12, 1,  21, 12, 4,
    21, 5,  12, 5,  8,  21, 12, 9,  11, 16, 12, 15, 16, 21, 12, 16, 13, 21, 12,
    19, 0,  0,  19, 2,  21, 0,  19, 2,  21, 0,  19, 3,  16, 10, 2,  19, 2,  12,
    2,  19, 1,  21, 19, 4,  16, 9,  19, 4,  21, 5,  19, 5,  8,  21, 19, 16, 13,
    21, 19, 7,  6,  21, 21, 2,  0,  0,  2,  21, 1,  0,  21, 4,  2,  11, 21, 4,
    15, 9,  21, 4,  0,  5,  21, 5,  8,  0,  21, 3,  7,  0,  16, 21, 0,  13, 7,
    21, 6,  0,  21, 9,  18, 0,  7,  9,  10, 5,  9,  4,  10, 0,  9,  5,  17, 0,
    18, 2,  9,  0,  2,  18, 1,  9,  18, 4,  9,  5,  18, 5,  2,  17, 18, 5,  8,
    9,  18, 3,  7,  9,  18, 16, 13, 9,  18, 9,  11, 5,  18, 10, 17, 9,  2,  17,
    1,  10, 2,  17, 8,  9,  17, 4,  10, 5,  17, 5,  8,  10, 17, 15, 16, 10, 17,
    16, 13, 10, 17, 8,  9,  10, 22, 2,  18, 0,  2,  22, 1,  18, 2,  22, 16, 5,
    22, 4,  18, 5,  22, 5,  8,  18, 22, 10, 18, 17, 22, 16, 18, 13, 22, 3,  18,
    7,  22, 7,  6,  18, 2,  10, 9,  5,  10, 4,  9,  3,  7,  10, 17, 3,  10, 9,
    3,  0,  11, 0,  4,  0,  2,  11, 0,  5,  2,  11, 0,  5,  11, 4,  9,  0,  11,
    4,  0,  3,  11, 5,  7,  0,  11, 9,  0,  10};
  const int nuv[NUM_GAS_REACTIONS * 4] = {
    -1, 2,  0, 0, -1, -1, 1, 1, -1, -1, 1, 1, -2, 1,  0, 0, -1, -1, 1, 1,
    -1, -1, 1, 0, -1, -1, 2, 0, -1, -1, 1, 0, -1, 2,  0, 0, -1, -1, 1, 1,
    -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 2, 0,
    -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -2, 1,  1, 0,
    -2, 1,  1, 0, -1, -1, 1, 0, -1, -1, 1, 0, -1, 1,  1, 0, -1, -1, 1, 1,
    -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1,
    -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1,
    -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1,
    -2, 1,  1, 0, -2, 1,  1, 0, -2, 2,  1, 0, -2, 1,  1, 0, -2, 1,  0, 0,
    -2, 1,  1, 0, -2, 1,  1, 0, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1,
    -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1,
    -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1,
    -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1,
    -1, 1,  1, 0, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1,
    -1, -1, 1, 1, -1, -1, 1, 1, -1, 1,  1, 0, -1, -1, 1, 1, -1, -1, 1, 1,
    -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1,
    -1, -1, 1, 1, -1, -1, 1, 1, -1, 1,  1, 0, -1, 1,  1, 0, -1, -1, 1, 1,
    -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1,
    -1, 1,  0, 0, -1, 1,  1, 0, -1, 1,  1, 0, -1, -1, 1, 1, -1, -1, 1, 1,
    -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1,
    -1, -1, 1, 1, -1, 1,  1, 0, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1,
    -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1,
    -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 0, -1, -1, 1, 0, -1, 1,  1, 0,
    -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1,
    -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1,
    -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -2, 1,  1, 1,
    -1, 1,  1, 0, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1,
    -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1,
    -1, -1, 1, 1, -1, -1, 1, 1, -2, 2,  1, 0, -1, 1,  1, 0, -1, -1, 1, 1,
    -1, -1, 1, 1, -1, -1, 2, 0, -1, -1, 1, 1, -1, -1, 1, 1, -1, -1, 1, 1};
  if (i < 1) {
    // Return max num species per reaction
    nspec = 4;
  } else {
    if (i > NUM_GAS_REACTIONS) {
      nspec = -1;
    } else {
      nspec = ns[i - 1];
      for (int j = 0; j < nspec; ++j) {
        ki[j] = kiv[(i - 1) * 4 + j] + 1;
        nu[j] = nuv[(i - 1) * 4 + j];
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
  amrex::Real c[24]; // temporary storage
  amrex::Real PORT =
    1e6 * P / (8.31446261815324e+07 * T); // 1e6 * P/RT so c goes to SI units

  // Compute conversion, see Eq 10
  for (int id = 0; id < 24; ++id) {
    c[id] = x[id] * PORT;
  }

  // convert to chemkin units
  progressRateFR(q_f, q_r, c, T);

  // convert to chemkin units
  for (int id = 0; id < 150; ++id) {
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
  amrex::Real g_RT[24];
  gibbs(g_RT, T);

  amrex::Real sc_qss[1];
  comp_qfqr(q_f, q_r, sc, sc_qss, T, invT, logT);
}

// save atomic weights into array
void
atomicWeight(amrex::Real* awt)
{
  awt[0] = 4.002602;  // He
  awt[1] = 39.950000; // Ar
  awt[2] = 14.007000; // N
  awt[3] = 15.999000; // O
  awt[4] = 1.008000;  // H
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
  for (int id = 0; id < kd * 24; ++id) {
    ncf[id] = 0;
  }

  // N2
  ncf[0 * kd + 2] = 2; // N

  // H2
  ncf[1 * kd + 4] = 2; // H

  // H
  ncf[2 * kd + 4] = 1; // H

  // O2
  ncf[3 * kd + 3] = 2; // O

  // O
  ncf[4 * kd + 3] = 1; // O

  // OH
  ncf[5 * kd + 4] = 1; // H
  ncf[5 * kd + 3] = 1; // O

  // H2O2
  ncf[6 * kd + 4] = 2; // H
  ncf[6 * kd + 3] = 2; // O

  // HO2
  ncf[7 * kd + 4] = 1; // H
  ncf[7 * kd + 3] = 2; // O

  // H2O
  ncf[8 * kd + 4] = 2; // H
  ncf[8 * kd + 3] = 1; // O

  // NO
  ncf[9 * kd + 2] = 1; // N
  ncf[9 * kd + 3] = 1; // O

  // NO2
  ncf[10 * kd + 2] = 1; // N
  ncf[10 * kd + 3] = 2; // O

  // N2O
  ncf[11 * kd + 2] = 2; // N
  ncf[11 * kd + 3] = 1; // O

  // N2H2
  ncf[12 * kd + 4] = 2; // H
  ncf[12 * kd + 2] = 2; // N

  // NH3
  ncf[13 * kd + 4] = 3; // H
  ncf[13 * kd + 2] = 1; // N

  // N2H4
  ncf[14 * kd + 4] = 4; // H
  ncf[14 * kd + 2] = 2; // N

  // NH
  ncf[15 * kd + 4] = 1; // H
  ncf[15 * kd + 2] = 1; // N

  // NH2
  ncf[16 * kd + 4] = 2; // H
  ncf[16 * kd + 2] = 1; // N

  // HONO
  ncf[17 * kd + 4] = 1; // H
  ncf[17 * kd + 2] = 1; // N
  ncf[17 * kd + 3] = 2; // O

  // HNO
  ncf[18 * kd + 4] = 1; // H
  ncf[18 * kd + 2] = 1; // N
  ncf[18 * kd + 3] = 1; // O

  // H2NN
  ncf[19 * kd + 4] = 2; // H
  ncf[19 * kd + 2] = 2; // N

  // N
  ncf[20 * kd + 2] = 1; // N

  // NNH
  ncf[21 * kd + 4] = 1; // H
  ncf[21 * kd + 2] = 2; // N

  // H2NO
  ncf[22 * kd + 4] = 2; // H
  ncf[22 * kd + 2] = 1; // N
  ncf[22 * kd + 3] = 1; // O

  // N2H3
  ncf[23 * kd + 4] = 3; // H
  ncf[23 * kd + 2] = 2; // N
}

// Returns the vector of strings of element names
void
CKSYME_STR(amrex::Vector<std::string>& ename)
{
  ename.resize(5);
  ename[0] = "He";
  ename[1] = "Ar";
  ename[2] = "N";
  ename[3] = "O";
  ename[4] = "H";
}

// Returns the vector of strings of species names
void
CKSYMS_STR(amrex::Vector<std::string>& kname)
{
  kname.resize(24);
  kname[0] = "N2";
  kname[1] = "H2";
  kname[2] = "H";
  kname[3] = "O2";
  kname[4] = "O";
  kname[5] = "OH";
  kname[6] = "H2O2";
  kname[7] = "HO2";
  kname[8] = "H2O";
  kname[9] = "NO";
  kname[10] = "NO2";
  kname[11] = "N2O";
  kname[12] = "N2H2";
  kname[13] = "NH3";
  kname[14] = "N2H4";
  kname[15] = "NH";
  kname[16] = "NH2";
  kname[17] = "HONO";
  kname[18] = "HNO";
  kname[19] = "H2NN";
  kname[20] = "N";
  kname[21] = "NNH";
  kname[22] = "H2NO";
  kname[23] = "N2H3";
}

// compute the sparsity pattern of the chemistry Jacobian
void
SPARSITY_INFO(int* nJdata, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 625> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 24> conc = {0.0};
  for (int n = 0; n < 24; n++) {
    conc[n] = 1.0 / 24.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 25; k++) {
    for (int l = 0; l < 25; l++) {
      if (Jac[25 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 625> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 24> conc = {0.0};
  for (int n = 0; n < 24; n++) {
    conc[n] = 1.0 / 24.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 25; k++) {
    for (int l = 0; l < 25; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[25 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 625> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 24> conc = {0.0};
  for (int n = 0; n < 24; n++) {
    conc[n] = 1.0 / 24.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 25; k++) {
    for (int l = 0; l < 25; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[25 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 625> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 24> conc = {0.0};
  for (int n = 0; n < 24; n++) {
    conc[n] = 1.0 / 24.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int nc = 0; nc < NCELLS; nc++) {
    int offset_row = nc * 25;
    int offset_col = nc * 25;
    for (int k = 0; k < 25; k++) {
      for (int l = 0; l < 25; l++) {
        if (Jac[25 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 625> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 24> conc = {0.0};
  for (int n = 0; n < 24; n++) {
    conc[n] = 1.0 / 24.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtrs[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 25;
      for (int l = 0; l < 25; l++) {
        for (int k = 0; k < 25; k++) {
          if (Jac[25 * k + l] != 0.0) {
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
      int offset = nc * 25;
      for (int l = 0; l < 25; l++) {
        for (int k = 0; k < 25; k++) {
          if (Jac[25 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 625> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 24> conc = {0.0};
  for (int n = 0; n < 24; n++) {
    conc[n] = 1.0 / 24.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 25;
      for (int l = 0; l < 25; l++) {
        for (int k = 0; k < 25; k++) {
          if (k == l) {
            colVals[nJdata_tmp - 1] = l + 1 + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[25 * k + l] != 0.0) {
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
      int offset = nc * 25;
      for (int l = 0; l < 25; l++) {
        for (int k = 0; k < 25; k++) {
          if (k == l) {
            colVals[nJdata_tmp] = l + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[25 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 625> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 24> conc = {0.0};
  for (int n = 0; n < 24; n++) {
    conc[n] = 1.0 / 24.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int k = 0; k < 25; k++) {
    for (int l = 0; l < 25; l++) {
      if (k == l) {
        rowVals[nJdata_tmp] = l;
        indx[nJdata_tmp] = 25 * k + l;
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[25 * k + l] != 0.0) {
          rowVals[nJdata_tmp] = l;
          indx[nJdata_tmp] = 25 * k + l;
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
  amrex::GpuArray<amrex::Real, 625> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 24> conc = {0.0};
  for (int n = 0; n < 24; n++) {
    conc[n] = 1.0 / 24.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int l = 0; l < 25; l++) {
      for (int k = 0; k < 25; k++) {
        if (k == l) {
          colVals[nJdata_tmp - 1] = l + 1;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[25 * k + l] != 0.0) {
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
    for (int l = 0; l < 25; l++) {
      for (int k = 0; k < 25; k++) {
        if (k == l) {
          colVals[nJdata_tmp] = l;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[25 * k + l] != 0.0) {
            colVals[nJdata_tmp] = k;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
      }
      rowPtr[l + 1] = nJdata_tmp;
    }
  }
}
