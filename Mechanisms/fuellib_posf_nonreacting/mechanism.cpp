#include "mechanism.H"

// Returns 0-based map of reaction order
void
GET_RMAP(int* /*_rmap*/)
{
}

// Returns a count of gas species in a gas reaction, and their indices
// and stoichiometric coefficients. (Eq 50)
void
CKINU(const int i, int& nspec, int* /*ki*/, int* /*nu*/)
{
  if (i < 1) {
    // Return max num species per reaction
    nspec = 0;
  } else {
    nspec = -1;
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
  amrex::Real c[93]; // temporary storage
  amrex::Real PORT =
    1e6 * P / (8.31446261815324e+07 * T); // 1e6 * P/RT so c goes to SI units

  // Compute conversion, see Eq 10
  for (int id = 0; id < 93; ++id) {
    c[id] = x[id] * PORT;
  }

  // convert to chemkin units
  progressRateFR(q_f, q_r, c, T);
}

// compute the progress rate for each reaction
// USES progressRate : todo switch to GPU
void
progressRateFR(
  amrex::Real* /*q_f*/,
  amrex::Real* /*q_r*/,
  amrex::Real* /*sc*/,
  amrex::Real /*T*/)
{
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
  awt[5] = 4.002602;  // He
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
  int kd = 6;
  // Zero ncf
  for (int id = 0; id < kd * 93; ++id) {
    ncf[id] = 0;
  }

  // O2
  ncf[0 * kd + 0] = 2; // O

  // N2
  ncf[1 * kd + 3] = 2; // N

  // H2O
  ncf[2 * kd + 1] = 2; // H
  ncf[2 * kd + 0] = 1; // O

  // CO2
  ncf[3 * kd + 2] = 1; // C
  ncf[3 * kd + 0] = 2; // O

  // C6H5CH3
  ncf[4 * kd + 2] = 7; // C
  ncf[4 * kd + 1] = 8; // H

  // C6H5C2H5
  ncf[5 * kd + 2] = 8;  // C
  ncf[5 * kd + 1] = 10; // H

  // C6H5C3H7
  ncf[6 * kd + 2] = 9;  // C
  ncf[6 * kd + 1] = 12; // H

  // C6H5C4H9
  ncf[7 * kd + 2] = 10; // C
  ncf[7 * kd + 1] = 14; // H

  // C6H5C5H11
  ncf[8 * kd + 2] = 11; // C
  ncf[8 * kd + 1] = 16; // H

  // C6H5C6H13
  ncf[9 * kd + 2] = 12; // C
  ncf[9 * kd + 1] = 18; // H

  // C6H5C7H15
  ncf[10 * kd + 2] = 13; // C
  ncf[10 * kd + 1] = 20; // H

  // C6H5C8H17
  ncf[11 * kd + 2] = 14; // C
  ncf[11 * kd + 1] = 22; // H

  // C6H5C9H19
  ncf[12 * kd + 2] = 15; // C
  ncf[12 * kd + 1] = 24; // H

  // C6H5C10H21
  ncf[13 * kd + 2] = 16; // C
  ncf[13 * kd + 1] = 26; // H

  // NAPH
  ncf[14 * kd + 2] = 10; // C
  ncf[14 * kd + 1] = 8;  // H

  // METHYLNAPH-1
  ncf[15 * kd + 2] = 11; // C
  ncf[15 * kd + 1] = 10; // H

  // ETHYLNAPH-1
  ncf[16 * kd + 2] = 12; // C
  ncf[16 * kd + 1] = 12; // H

  // PROPYLNAPH-1
  ncf[17 * kd + 2] = 13; // C
  ncf[17 * kd + 1] = 14; // H

  // BUTYLNAPH-1
  ncf[18 * kd + 2] = 14; // C
  ncf[18 * kd + 1] = 16; // H

  // INDANE
  ncf[19 * kd + 2] = 9;  // C
  ncf[19 * kd + 1] = 10; // H

  // TETRA
  ncf[20 * kd + 2] = 10; // C
  ncf[20 * kd + 1] = 12; // H

  // METHYLTETRA-2
  ncf[21 * kd + 2] = 11; // C
  ncf[21 * kd + 1] = 14; // H

  // ETHYLTETRA-2
  ncf[22 * kd + 2] = 12; // C
  ncf[22 * kd + 1] = 16; // H

  // PROPYLTETRA-2
  ncf[23 * kd + 2] = 13; // C
  ncf[23 * kd + 1] = 18; // H

  // BUTYLTETRA-2
  ncf[24 * kd + 2] = 14; // C
  ncf[24 * kd + 1] = 20; // H

  // PENTYLTETRA-2
  ncf[25 * kd + 2] = 15; // C
  ncf[25 * kd + 1] = 22; // H

  // C7H16-2
  ncf[26 * kd + 2] = 7;  // C
  ncf[26 * kd + 1] = 16; // H

  // C8H18-2
  ncf[27 * kd + 2] = 8;  // C
  ncf[27 * kd + 1] = 18; // H

  // C9H20-2
  ncf[28 * kd + 2] = 9;  // C
  ncf[28 * kd + 1] = 20; // H

  // C10H22-2
  ncf[29 * kd + 2] = 10; // C
  ncf[29 * kd + 1] = 22; // H

  // C11H24-2
  ncf[30 * kd + 2] = 11; // C
  ncf[30 * kd + 1] = 24; // H

  // C12H26-2
  ncf[31 * kd + 2] = 12; // C
  ncf[31 * kd + 1] = 26; // H

  // C13H28-2
  ncf[32 * kd + 2] = 13; // C
  ncf[32 * kd + 1] = 28; // H

  // C14H30-2
  ncf[33 * kd + 2] = 14; // C
  ncf[33 * kd + 1] = 30; // H

  // C15H32-2
  ncf[34 * kd + 2] = 15; // C
  ncf[34 * kd + 1] = 32; // H

  // C16H34-2
  ncf[35 * kd + 2] = 16; // C
  ncf[35 * kd + 1] = 34; // H

  // C17H36-2
  ncf[36 * kd + 2] = 17; // C
  ncf[36 * kd + 1] = 36; // H

  // C18H38-2
  ncf[37 * kd + 2] = 18; // C
  ncf[37 * kd + 1] = 38; // H

  // C19H40-2
  ncf[38 * kd + 2] = 19; // C
  ncf[38 * kd + 1] = 40; // H

  // C20H42-2
  ncf[39 * kd + 2] = 20; // C
  ncf[39 * kd + 1] = 42; // H

  // C21H44-2
  ncf[40 * kd + 2] = 21; // C
  ncf[40 * kd + 1] = 44; // H

  // C22H46-2
  ncf[41 * kd + 2] = 22; // C
  ncf[41 * kd + 1] = 46; // H

  // C23H48-2
  ncf[42 * kd + 2] = 23; // C
  ncf[42 * kd + 1] = 48; // H

  // C24H50-2
  ncf[43 * kd + 2] = 24; // C
  ncf[43 * kd + 1] = 50; // H

  // NC7H16
  ncf[44 * kd + 2] = 7;  // C
  ncf[44 * kd + 1] = 16; // H

  // NC8H18
  ncf[45 * kd + 2] = 8;  // C
  ncf[45 * kd + 1] = 18; // H

  // NC9H20
  ncf[46 * kd + 2] = 9;  // C
  ncf[46 * kd + 1] = 20; // H

  // NC10H22
  ncf[47 * kd + 2] = 10; // C
  ncf[47 * kd + 1] = 22; // H

  // NC11H24
  ncf[48 * kd + 2] = 11; // C
  ncf[48 * kd + 1] = 24; // H

  // NC12H26
  ncf[49 * kd + 2] = 12; // C
  ncf[49 * kd + 1] = 26; // H

  // NC13H28
  ncf[50 * kd + 2] = 13; // C
  ncf[50 * kd + 1] = 28; // H

  // NC14H30
  ncf[51 * kd + 2] = 14; // C
  ncf[51 * kd + 1] = 30; // H

  // NC15H32
  ncf[52 * kd + 2] = 15; // C
  ncf[52 * kd + 1] = 32; // H

  // NC16H34
  ncf[53 * kd + 2] = 16; // C
  ncf[53 * kd + 1] = 34; // H

  // NC17H36
  ncf[54 * kd + 2] = 17; // C
  ncf[54 * kd + 1] = 36; // H

  // NC18H38
  ncf[55 * kd + 2] = 18; // C
  ncf[55 * kd + 1] = 38; // H

  // NC19H40
  ncf[56 * kd + 2] = 19; // C
  ncf[56 * kd + 1] = 40; // H

  // NC20H42
  ncf[57 * kd + 2] = 20; // C
  ncf[57 * kd + 1] = 42; // H

  // NC21H44
  ncf[58 * kd + 2] = 21; // C
  ncf[58 * kd + 1] = 44; // H

  // NC22H46
  ncf[59 * kd + 2] = 22; // C
  ncf[59 * kd + 1] = 46; // H

  // NC23H48
  ncf[60 * kd + 2] = 23; // C
  ncf[60 * kd + 1] = 48; // H

  // C6H11CH3
  ncf[61 * kd + 2] = 7;  // C
  ncf[61 * kd + 1] = 14; // H

  // C6H11C2H5
  ncf[62 * kd + 2] = 8;  // C
  ncf[62 * kd + 1] = 16; // H

  // C6H11C3H7
  ncf[63 * kd + 2] = 9;  // C
  ncf[63 * kd + 1] = 18; // H

  // C6H11C4H9
  ncf[64 * kd + 2] = 10; // C
  ncf[64 * kd + 1] = 20; // H

  // C6H11C5H11
  ncf[65 * kd + 2] = 11; // C
  ncf[65 * kd + 1] = 22; // H

  // C6H11C6H13
  ncf[66 * kd + 2] = 12; // C
  ncf[66 * kd + 1] = 24; // H

  // C6H11C7H15
  ncf[67 * kd + 2] = 13; // C
  ncf[67 * kd + 1] = 26; // H

  // C6H11C8H17
  ncf[68 * kd + 2] = 14; // C
  ncf[68 * kd + 1] = 28; // H

  // C6H11C9H19
  ncf[69 * kd + 2] = 15; // C
  ncf[69 * kd + 1] = 30; // H

  // C6H11C10H21
  ncf[70 * kd + 2] = 16; // C
  ncf[70 * kd + 1] = 32; // H

  // C6H11C11H23
  ncf[71 * kd + 2] = 17; // C
  ncf[71 * kd + 1] = 34; // H

  // C6H11C12H25
  ncf[72 * kd + 2] = 18; // C
  ncf[72 * kd + 1] = 36; // H

  // C6H11C13H27
  ncf[73 * kd + 2] = 19; // C
  ncf[73 * kd + 1] = 38; // H

  // OHPEN
  ncf[74 * kd + 2] = 8;  // C
  ncf[74 * kd + 1] = 14; // H

  // HYDRIND
  ncf[75 * kd + 2] = 9;  // C
  ncf[75 * kd + 1] = 16; // H

  // DECALIN
  ncf[76 * kd + 2] = 10; // C
  ncf[76 * kd + 1] = 18; // H

  // METHYLDECALIN-2
  ncf[77 * kd + 2] = 11; // C
  ncf[77 * kd + 1] = 20; // H

  // ETHYLDECALIN-2
  ncf[78 * kd + 2] = 12; // C
  ncf[78 * kd + 1] = 22; // H

  // PROPYLDECALIN-2
  ncf[79 * kd + 2] = 13; // C
  ncf[79 * kd + 1] = 24; // H

  // BUTYLDECALIN-2
  ncf[80 * kd + 2] = 14; // C
  ncf[80 * kd + 1] = 26; // H

  // PENTYLDECALIN-2
  ncf[81 * kd + 2] = 15; // C
  ncf[81 * kd + 1] = 28; // H

  // HEXYLDECALIN-2
  ncf[82 * kd + 2] = 16; // C
  ncf[82 * kd + 1] = 30; // H

  // HEPTYLDECALIN-2
  ncf[83 * kd + 2] = 17; // C
  ncf[83 * kd + 1] = 32; // H

  // TRICYCLO-C10
  ncf[84 * kd + 2] = 10; // C
  ncf[84 * kd + 1] = 16; // H

  // TRICYCLO-C11
  ncf[85 * kd + 2] = 11; // C
  ncf[85 * kd + 1] = 18; // H

  // TRICYCLO-C12
  ncf[86 * kd + 2] = 12; // C
  ncf[86 * kd + 1] = 20; // H

  // TRICYCLO-C13
  ncf[87 * kd + 2] = 13; // C
  ncf[87 * kd + 1] = 22; // H

  // TRICYCLO-C14
  ncf[88 * kd + 2] = 14; // C
  ncf[88 * kd + 1] = 24; // H

  // DECENE-1
  ncf[89 * kd + 2] = 10; // C
  ncf[89 * kd + 1] = 20; // H

  // DODECENE-1
  ncf[90 * kd + 2] = 12; // C
  ncf[90 * kd + 1] = 24; // H

  // TETRADECENE-1
  ncf[91 * kd + 2] = 14; // C
  ncf[91 * kd + 1] = 28; // H

  // HEXADECENE-1
  ncf[92 * kd + 2] = 16; // C
  ncf[92 * kd + 1] = 32; // H
}

// Returns the vector of strings of element names
void
CKSYME_STR(amrex::Vector<std::string>& ename)
{
  ename.resize(6);
  ename[0] = "O";
  ename[1] = "H";
  ename[2] = "C";
  ename[3] = "N";
  ename[4] = "Ar";
  ename[5] = "He";
}

// Returns the vector of strings of species names
void
CKSYMS_STR(amrex::Vector<std::string>& kname)
{
  kname.resize(93);
  kname[0] = "O2";
  kname[1] = "N2";
  kname[2] = "H2O";
  kname[3] = "CO2";
  kname[4] = "C6H5CH3";
  kname[5] = "C6H5C2H5";
  kname[6] = "C6H5C3H7";
  kname[7] = "C6H5C4H9";
  kname[8] = "C6H5C5H11";
  kname[9] = "C6H5C6H13";
  kname[10] = "C6H5C7H15";
  kname[11] = "C6H5C8H17";
  kname[12] = "C6H5C9H19";
  kname[13] = "C6H5C10H21";
  kname[14] = "NAPH";
  kname[15] = "METHYLNAPH-1";
  kname[16] = "ETHYLNAPH-1";
  kname[17] = "PROPYLNAPH-1";
  kname[18] = "BUTYLNAPH-1";
  kname[19] = "INDANE";
  kname[20] = "TETRA";
  kname[21] = "METHYLTETRA-2";
  kname[22] = "ETHYLTETRA-2";
  kname[23] = "PROPYLTETRA-2";
  kname[24] = "BUTYLTETRA-2";
  kname[25] = "PENTYLTETRA-2";
  kname[26] = "C7H16-2";
  kname[27] = "C8H18-2";
  kname[28] = "C9H20-2";
  kname[29] = "C10H22-2";
  kname[30] = "C11H24-2";
  kname[31] = "C12H26-2";
  kname[32] = "C13H28-2";
  kname[33] = "C14H30-2";
  kname[34] = "C15H32-2";
  kname[35] = "C16H34-2";
  kname[36] = "C17H36-2";
  kname[37] = "C18H38-2";
  kname[38] = "C19H40-2";
  kname[39] = "C20H42-2";
  kname[40] = "C21H44-2";
  kname[41] = "C22H46-2";
  kname[42] = "C23H48-2";
  kname[43] = "C24H50-2";
  kname[44] = "NC7H16";
  kname[45] = "NC8H18";
  kname[46] = "NC9H20";
  kname[47] = "NC10H22";
  kname[48] = "NC11H24";
  kname[49] = "NC12H26";
  kname[50] = "NC13H28";
  kname[51] = "NC14H30";
  kname[52] = "NC15H32";
  kname[53] = "NC16H34";
  kname[54] = "NC17H36";
  kname[55] = "NC18H38";
  kname[56] = "NC19H40";
  kname[57] = "NC20H42";
  kname[58] = "NC21H44";
  kname[59] = "NC22H46";
  kname[60] = "NC23H48";
  kname[61] = "C6H11CH3";
  kname[62] = "C6H11C2H5";
  kname[63] = "C6H11C3H7";
  kname[64] = "C6H11C4H9";
  kname[65] = "C6H11C5H11";
  kname[66] = "C6H11C6H13";
  kname[67] = "C6H11C7H15";
  kname[68] = "C6H11C8H17";
  kname[69] = "C6H11C9H19";
  kname[70] = "C6H11C10H21";
  kname[71] = "C6H11C11H23";
  kname[72] = "C6H11C12H25";
  kname[73] = "C6H11C13H27";
  kname[74] = "OHPEN";
  kname[75] = "HYDRIND";
  kname[76] = "DECALIN";
  kname[77] = "METHYLDECALIN-2";
  kname[78] = "ETHYLDECALIN-2";
  kname[79] = "PROPYLDECALIN-2";
  kname[80] = "BUTYLDECALIN-2";
  kname[81] = "PENTYLDECALIN-2";
  kname[82] = "HEXYLDECALIN-2";
  kname[83] = "HEPTYLDECALIN-2";
  kname[84] = "TRICYCLO-C10";
  kname[85] = "TRICYCLO-C11";
  kname[86] = "TRICYCLO-C12";
  kname[87] = "TRICYCLO-C13";
  kname[88] = "TRICYCLO-C14";
  kname[89] = "DECENE-1";
  kname[90] = "DODECENE-1";
  kname[91] = "TETRADECENE-1";
  kname[92] = "HEXADECENE-1";
}

// compute the sparsity pattern of the chemistry Jacobian
void
SPARSITY_INFO(int* nJdata, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 8836> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 93> conc = {0.0};
  for (int n = 0; n < 93; n++) {
    conc[n] = 1.0 / 93.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 94; k++) {
    for (int l = 0; l < 94; l++) {
      if (Jac[94 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 8836> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 93> conc = {0.0};
  for (int n = 0; n < 93; n++) {
    conc[n] = 1.0 / 93.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 94; k++) {
    for (int l = 0; l < 94; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[94 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 8836> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 93> conc = {0.0};
  for (int n = 0; n < 93; n++) {
    conc[n] = 1.0 / 93.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 94; k++) {
    for (int l = 0; l < 94; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[94 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 8836> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 93> conc = {0.0};
  for (int n = 0; n < 93; n++) {
    conc[n] = 1.0 / 93.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int nc = 0; nc < NCELLS; nc++) {
    int offset_row = nc * 94;
    int offset_col = nc * 94;
    for (int k = 0; k < 94; k++) {
      for (int l = 0; l < 94; l++) {
        if (Jac[94 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 8836> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 93> conc = {0.0};
  for (int n = 0; n < 93; n++) {
    conc[n] = 1.0 / 93.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtrs[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 94;
      for (int l = 0; l < 94; l++) {
        for (int k = 0; k < 94; k++) {
          if (Jac[94 * k + l] != 0.0) {
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
      int offset = nc * 94;
      for (int l = 0; l < 94; l++) {
        for (int k = 0; k < 94; k++) {
          if (Jac[94 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 8836> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 93> conc = {0.0};
  for (int n = 0; n < 93; n++) {
    conc[n] = 1.0 / 93.000000;
  }
  aJacobian(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 94;
      for (int l = 0; l < 94; l++) {
        for (int k = 0; k < 94; k++) {
          if (k == l) {
            colVals[nJdata_tmp - 1] = l + 1 + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[94 * k + l] != 0.0) {
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
      int offset = nc * 94;
      for (int l = 0; l < 94; l++) {
        for (int k = 0; k < 94; k++) {
          if (k == l) {
            colVals[nJdata_tmp] = l + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[94 * k + l] != 0.0) {
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
  amrex::GpuArray<amrex::Real, 8836> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 93> conc = {0.0};
  for (int n = 0; n < 93; n++) {
    conc[n] = 1.0 / 93.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int k = 0; k < 94; k++) {
    for (int l = 0; l < 94; l++) {
      if (k == l) {
        rowVals[nJdata_tmp] = l;
        indx[nJdata_tmp] = 94 * k + l;
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[94 * k + l] != 0.0) {
          rowVals[nJdata_tmp] = l;
          indx[nJdata_tmp] = 94 * k + l;
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
  amrex::GpuArray<amrex::Real, 8836> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 93> conc = {0.0};
  for (int n = 0; n < 93; n++) {
    conc[n] = 1.0 / 93.000000;
  }
  aJacobian_precond(Jac.data(), conc.data(), 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int l = 0; l < 94; l++) {
      for (int k = 0; k < 94; k++) {
        if (k == l) {
          colVals[nJdata_tmp - 1] = l + 1;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[94 * k + l] != 0.0) {
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
    for (int l = 0; l < 94; l++) {
      for (int k = 0; k < 94; k++) {
        if (k == l) {
          colVals[nJdata_tmp] = l;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[94 * k + l] != 0.0) {
            colVals[nJdata_tmp] = k;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
      }
      rowPtr[l + 1] = nJdata_tmp;
    }
  }
}
