#ifndef MECHANISM_CPP
#define MECHANISM_CPP

#include "mechanism.H"

/*save atomic weights into array */
void
atomicWeight(amrex::Real* awt)
{
  awt[0] = 12.011150; /*C */
  awt[1] = 1.007970;  /*H */
  awt[2] = 15.999400; /*O */
  awt[3] = 14.006700; /*N */
  awt[4] = 39.948000; /*AR */
  awt[5] = 4.002600;  /*HE */
}

/*get atomic weight for all elements */
void
CKAWT(amrex::Real* awt)
{
  atomicWeight(awt);
}

/*Returns the elemental composition  */
/*of the speciesi (mdim is num of elements) */
void
CKNCF(int* ncf)
{
  int id; /*loop counter */
  int kd = 6;
  /*Zero ncf */
  for (id = 0; id < kd * 130; ++id) {
    ncf[id] = 0;
  }

  /*H */
  ncf[0 * kd + 1] = 1; /*H */

  /*H2 */
  ncf[1 * kd + 1] = 2; /*H */

  /*CH2 */
  ncf[2 * kd + 0] = 1; /*C */
  ncf[2 * kd + 1] = 2; /*H */

  /*CH2(S) */
  ncf[3 * kd + 0] = 1; /*C */
  ncf[3 * kd + 1] = 2; /*H */

  /*CH3 */
  ncf[4 * kd + 0] = 1; /*C */
  ncf[4 * kd + 1] = 3; /*H */

  /*O */
  ncf[5 * kd + 2] = 1; /*O */

  /*CH4 */
  ncf[6 * kd + 0] = 1; /*C */
  ncf[6 * kd + 1] = 4; /*H */

  /*OH */
  ncf[7 * kd + 2] = 1; /*O */
  ncf[7 * kd + 1] = 1; /*H */

  /*H2O */
  ncf[8 * kd + 1] = 2; /*H */
  ncf[8 * kd + 2] = 1; /*O */

  /*C2H */
  ncf[9 * kd + 0] = 2; /*C */
  ncf[9 * kd + 1] = 1; /*H */

  /*C2H2 */
  ncf[10 * kd + 0] = 2; /*C */
  ncf[10 * kd + 1] = 2; /*H */

  /*C2H3 */
  ncf[11 * kd + 0] = 2; /*C */
  ncf[11 * kd + 1] = 3; /*H */

  /*CO */
  ncf[12 * kd + 0] = 1; /*C */
  ncf[12 * kd + 2] = 1; /*O */

  /*N2 */
  ncf[13 * kd + 3] = 2; /*N */

  /*C2H4 */
  ncf[14 * kd + 0] = 2; /*C */
  ncf[14 * kd + 1] = 4; /*H */

  /*HCO */
  ncf[15 * kd + 0] = 1; /*C */
  ncf[15 * kd + 1] = 1; /*H */
  ncf[15 * kd + 2] = 1; /*O */

  /*C2H5 */
  ncf[16 * kd + 0] = 2; /*C */
  ncf[16 * kd + 1] = 5; /*H */

  /*CH2O */
  ncf[17 * kd + 1] = 2; /*H */
  ncf[17 * kd + 0] = 1; /*C */
  ncf[17 * kd + 2] = 1; /*O */

  /*C2H6 */
  ncf[18 * kd + 0] = 2; /*C */
  ncf[18 * kd + 1] = 6; /*H */

  /*CH2OH */
  ncf[19 * kd + 0] = 1; /*C */
  ncf[19 * kd + 1] = 3; /*H */
  ncf[19 * kd + 2] = 1; /*O */

  /*CH3O */
  ncf[20 * kd + 0] = 1; /*C */
  ncf[20 * kd + 1] = 3; /*H */
  ncf[20 * kd + 2] = 1; /*O */

  /*O2 */
  ncf[21 * kd + 2] = 2; /*O */

  /*CH3OH */
  ncf[22 * kd + 0] = 1; /*C */
  ncf[22 * kd + 1] = 4; /*H */
  ncf[22 * kd + 2] = 1; /*O */

  /*HO2 */
  ncf[23 * kd + 1] = 1; /*H */
  ncf[23 * kd + 2] = 2; /*O */

  /*H2O2 */
  ncf[24 * kd + 1] = 2; /*H */
  ncf[24 * kd + 2] = 2; /*O */

  /*HCCO */
  ncf[25 * kd + 1] = 1; /*H */
  ncf[25 * kd + 0] = 2; /*C */
  ncf[25 * kd + 2] = 1; /*O */

  /*CH2CO */
  ncf[26 * kd + 0] = 2; /*C */
  ncf[26 * kd + 1] = 2; /*H */
  ncf[26 * kd + 2] = 1; /*O */

  /*HCCOH */
  ncf[27 * kd + 1] = 2; /*H */
  ncf[27 * kd + 0] = 2; /*C */
  ncf[27 * kd + 2] = 1; /*O */

  /*CH2HCO */
  ncf[28 * kd + 2] = 1; /*O */
  ncf[28 * kd + 1] = 3; /*H */
  ncf[28 * kd + 0] = 2; /*C */

  /*CH3CO */
  ncf[29 * kd + 0] = 2; /*C */
  ncf[29 * kd + 1] = 3; /*H */
  ncf[29 * kd + 2] = 1; /*O */

  /*CO2 */
  ncf[30 * kd + 0] = 1; /*C */
  ncf[30 * kd + 2] = 2; /*O */

  /*CH3HCO */
  ncf[31 * kd + 0] = 2; /*C */
  ncf[31 * kd + 1] = 4; /*H */
  ncf[31 * kd + 2] = 1; /*O */

  /*OCHO */
  ncf[32 * kd + 0] = 1; /*C */
  ncf[32 * kd + 1] = 1; /*H */
  ncf[32 * kd + 2] = 2; /*O */

  /*CH3CHOH */
  ncf[33 * kd + 0] = 2; /*C */
  ncf[33 * kd + 2] = 1; /*O */
  ncf[33 * kd + 1] = 5; /*H */

  /*C2H4OH */
  ncf[34 * kd + 1] = 5; /*H */
  ncf[34 * kd + 0] = 2; /*C */
  ncf[34 * kd + 2] = 1; /*O */

  /*CH3CH2O */
  ncf[35 * kd + 2] = 1; /*O */
  ncf[35 * kd + 0] = 2; /*C */
  ncf[35 * kd + 1] = 5; /*H */

  /*CH3OCH2 */
  ncf[36 * kd + 0] = 2; /*C */
  ncf[36 * kd + 1] = 5; /*H */
  ncf[36 * kd + 2] = 1; /*O */

  /*HOCHO */
  ncf[37 * kd + 1] = 2; /*H */
  ncf[37 * kd + 0] = 1; /*C */
  ncf[37 * kd + 2] = 2; /*O */

  /*CH3OCH3 */
  ncf[38 * kd + 0] = 2; /*C */
  ncf[38 * kd + 1] = 6; /*H */
  ncf[38 * kd + 2] = 1; /*O */

  /*C2H5OH */
  ncf[39 * kd + 0] = 2; /*C */
  ncf[39 * kd + 1] = 6; /*H */
  ncf[39 * kd + 2] = 1; /*O */

  /*HOCH2O */
  ncf[40 * kd + 0] = 1; /*C */
  ncf[40 * kd + 1] = 3; /*H */
  ncf[40 * kd + 2] = 2; /*O */

  /*CH3OCO */
  ncf[41 * kd + 0] = 2; /*C */
  ncf[41 * kd + 1] = 3; /*H */
  ncf[41 * kd + 2] = 2; /*O */

  /*CH3OCHO */
  ncf[42 * kd + 0] = 2; /*C */
  ncf[42 * kd + 1] = 4; /*H */
  ncf[42 * kd + 2] = 2; /*O */

  /*CH3OCH2O */
  ncf[43 * kd + 0] = 2; /*C */
  ncf[43 * kd + 1] = 5; /*H */
  ncf[43 * kd + 2] = 2; /*O */

  /*CH3OCH2OH */
  ncf[44 * kd + 0] = 2; /*C */
  ncf[44 * kd + 1] = 6; /*H */
  ncf[44 * kd + 2] = 2; /*O */

  /*OCH2OCHO */
  ncf[45 * kd + 0] = 2; /*C */
  ncf[45 * kd + 1] = 3; /*H */
  ncf[45 * kd + 2] = 3; /*O */

  /*HOCH2OCO */
  ncf[46 * kd + 0] = 2; /*C */
  ncf[46 * kd + 1] = 3; /*H */
  ncf[46 * kd + 2] = 3; /*O */

  /*HOC2H4O2 */
  ncf[47 * kd + 0] = 2; /*C */
  ncf[47 * kd + 1] = 5; /*H */
  ncf[47 * kd + 2] = 3; /*O */

  /*CH3OCH2O2 */
  ncf[48 * kd + 0] = 2; /*C */
  ncf[48 * kd + 1] = 5; /*H */
  ncf[48 * kd + 2] = 3; /*O */

  /*CH2OCH2O2H */
  ncf[49 * kd + 0] = 2; /*C */
  ncf[49 * kd + 1] = 5; /*H */
  ncf[49 * kd + 2] = 3; /*O */

  /*CH3OCH2O2H */
  ncf[50 * kd + 0] = 2; /*C */
  ncf[50 * kd + 1] = 6; /*H */
  ncf[50 * kd + 2] = 3; /*O */

  /*HO2CH2OCHO */
  ncf[51 * kd + 0] = 2; /*C */
  ncf[51 * kd + 1] = 4; /*H */
  ncf[51 * kd + 2] = 4; /*O */

  /*O2CH2OCH2O2H */
  ncf[52 * kd + 0] = 2; /*C */
  ncf[52 * kd + 1] = 5; /*H */
  ncf[52 * kd + 2] = 5; /*O */

  /*AR */
  ncf[53 * kd + 4] = 1; /*AR */

  /*HE */
  ncf[54 * kd + 5] = 1; /*HE */

  /*C */
  ncf[55 * kd + 0] = 1; /*C */

  /*CH */
  ncf[56 * kd + 0] = 1; /*C */
  ncf[56 * kd + 1] = 1; /*H */

  /*HCOH */
  ncf[57 * kd + 0] = 1; /*C */
  ncf[57 * kd + 1] = 2; /*H */
  ncf[57 * kd + 2] = 1; /*O */

  /*C2O */
  ncf[58 * kd + 0] = 2; /*C */
  ncf[58 * kd + 2] = 1; /*O */

  /*CH3O2 */
  ncf[59 * kd + 0] = 1; /*C */
  ncf[59 * kd + 1] = 3; /*H */
  ncf[59 * kd + 2] = 2; /*O */

  /*CH3O2H */
  ncf[60 * kd + 0] = 1; /*C */
  ncf[60 * kd + 1] = 4; /*H */
  ncf[60 * kd + 2] = 2; /*O */

  /*C3H2 */
  ncf[61 * kd + 1] = 2; /*H */
  ncf[61 * kd + 0] = 3; /*C */

  /*C3H2SING */
  ncf[62 * kd + 0] = 3; /*C */
  ncf[62 * kd + 1] = 2; /*H */

  /*c-C3H2 */
  ncf[63 * kd + 0] = 3; /*C */
  ncf[63 * kd + 1] = 2; /*H */

  /*C3H3 */
  ncf[64 * kd + 0] = 3; /*C */
  ncf[64 * kd + 1] = 3; /*H */

  /*C3H6 */
  ncf[65 * kd + 0] = 3; /*C */
  ncf[65 * kd + 1] = 6; /*H */

  /*C4H2 */
  ncf[66 * kd + 0] = 4; /*C */
  ncf[66 * kd + 1] = 2; /*H */

  /*pC3H4 */
  ncf[67 * kd + 1] = 4; /*H */
  ncf[67 * kd + 0] = 3; /*C */

  /*aC3H4 */
  ncf[68 * kd + 1] = 4; /*H */
  ncf[68 * kd + 0] = 3; /*C */

  /*c-C3H4 */
  ncf[69 * kd + 0] = 3; /*C */
  ncf[69 * kd + 1] = 4; /*H */

  /*aC3H5 */
  ncf[70 * kd + 0] = 3; /*C */
  ncf[70 * kd + 1] = 5; /*H */

  /*sC3H5 */
  ncf[71 * kd + 0] = 3; /*C */
  ncf[71 * kd + 1] = 5; /*H */

  /*tC3H5 */
  ncf[72 * kd + 0] = 3; /*C */
  ncf[72 * kd + 1] = 5; /*H */

  /*nC4H3 */
  ncf[73 * kd + 0] = 4; /*C */
  ncf[73 * kd + 1] = 3; /*H */

  /*iC4H3 */
  ncf[74 * kd + 0] = 4; /*C */
  ncf[74 * kd + 1] = 3; /*H */

  /*C4H4 */
  ncf[75 * kd + 0] = 4; /*C */
  ncf[75 * kd + 1] = 4; /*H */

  /*nC4H5 */
  ncf[76 * kd + 0] = 4; /*C */
  ncf[76 * kd + 1] = 5; /*H */

  /*iC4H5 */
  ncf[77 * kd + 0] = 4; /*C */
  ncf[77 * kd + 1] = 5; /*H */

  /*C4H6-13 */
  ncf[78 * kd + 0] = 4; /*C */
  ncf[78 * kd + 1] = 6; /*H */

  /*C4H6-12 */
  ncf[79 * kd + 0] = 4; /*C */
  ncf[79 * kd + 1] = 6; /*H */

  /*C4H6-2 */
  ncf[80 * kd + 0] = 4; /*C */
  ncf[80 * kd + 1] = 6; /*H */

  /*H2CC */
  ncf[81 * kd + 1] = 2; /*H */
  ncf[81 * kd + 0] = 2; /*C */

  /*CH2CHOH */
  ncf[82 * kd + 0] = 2; /*C */
  ncf[82 * kd + 1] = 4; /*H */
  ncf[82 * kd + 2] = 1; /*O */

  /*C5H2 */
  ncf[83 * kd + 0] = 5; /*C */
  ncf[83 * kd + 1] = 2; /*H */

  /*HCCCHCCH */
  ncf[84 * kd + 0] = 5; /*C */
  ncf[84 * kd + 1] = 3; /*H */

  /*H2CCCCCH */
  ncf[85 * kd + 0] = 5; /*C */
  ncf[85 * kd + 1] = 3; /*H */

  /*C4H5C2H */
  ncf[86 * kd + 0] = 6; /*C */
  ncf[86 * kd + 1] = 6; /*H */

  /*FC6H6 */
  ncf[87 * kd + 0] = 6; /*C */
  ncf[87 * kd + 1] = 6; /*H */

  /*C6H6 */
  ncf[88 * kd + 0] = 6; /*C */
  ncf[88 * kd + 1] = 6; /*H */

  /*C6H5 */
  ncf[89 * kd + 0] = 6; /*C */
  ncf[89 * kd + 1] = 5; /*H */

  /*C2H3CHO */
  ncf[90 * kd + 0] = 3; /*C */
  ncf[90 * kd + 1] = 4; /*H */
  ncf[90 * kd + 2] = 1; /*O */

  /*C5H6 */
  ncf[91 * kd + 0] = 5; /*C */
  ncf[91 * kd + 1] = 6; /*H */

  /*C2H3CO */
  ncf[92 * kd + 0] = 3; /*C */
  ncf[92 * kd + 1] = 3; /*H */
  ncf[92 * kd + 2] = 1; /*O */

  /*C4H */
  ncf[93 * kd + 0] = 4; /*C */
  ncf[93 * kd + 1] = 1; /*H */

  /*H2C4O */
  ncf[94 * kd + 1] = 2; /*H */
  ncf[94 * kd + 0] = 4; /*C */
  ncf[94 * kd + 2] = 1; /*O */

  /*C6H2 */
  ncf[95 * kd + 0] = 6; /*C */
  ncf[95 * kd + 1] = 2; /*H */

  /*C6H3 */
  ncf[96 * kd + 0] = 6; /*C */
  ncf[96 * kd + 1] = 3; /*H */

  /*l-C6H4 */
  ncf[97 * kd + 0] = 6; /*C */
  ncf[97 * kd + 1] = 4; /*H */

  /*o-C6H4 */
  ncf[98 * kd + 0] = 6; /*C */
  ncf[98 * kd + 1] = 4; /*H */

  /*CH3CHCHCO */
  ncf[99 * kd + 0] = 4; /*C */
  ncf[99 * kd + 1] = 5; /*H */
  ncf[99 * kd + 2] = 1; /*O */

  /*CH2CHCHCHO */
  ncf[100 * kd + 0] = 4; /*C */
  ncf[100 * kd + 1] = 5; /*H */
  ncf[100 * kd + 2] = 1; /*O */

  /*C5H5 */
  ncf[101 * kd + 0] = 5; /*C */
  ncf[101 * kd + 1] = 5; /*H */

  /*C6H5O */
  ncf[102 * kd + 0] = 6; /*C */
  ncf[102 * kd + 1] = 5; /*H */
  ncf[102 * kd + 2] = 1; /*O */

  /*C6H4O2 */
  ncf[103 * kd + 0] = 6; /*C */
  ncf[103 * kd + 1] = 4; /*H */
  ncf[103 * kd + 2] = 2; /*O */

  /*C6H5OH */
  ncf[104 * kd + 0] = 6; /*C */
  ncf[104 * kd + 1] = 6; /*H */
  ncf[104 * kd + 2] = 1; /*O */

  /*C5H4O */
  ncf[105 * kd + 0] = 5; /*C */
  ncf[105 * kd + 1] = 4; /*H */
  ncf[105 * kd + 2] = 1; /*O */

  /*C5H5O */
  ncf[106 * kd + 0] = 5; /*C */
  ncf[106 * kd + 1] = 5; /*H */
  ncf[106 * kd + 2] = 1; /*O */

  /*OSING */
  ncf[107 * kd + 2] = 1; /*O */

  /*O2SING */
  ncf[108 * kd + 2] = 2; /*O */

  /*O3 */
  ncf[109 * kd + 2] = 3; /*O */

  /*C2H4O2 */
  ncf[110 * kd + 0] = 2; /*C */
  ncf[110 * kd + 1] = 4; /*H */
  ncf[110 * kd + 2] = 2; /*O */

  /*CH2OCHO */
  ncf[111 * kd + 0] = 2; /*C */
  ncf[111 * kd + 1] = 3; /*H */
  ncf[111 * kd + 2] = 2; /*O */

  /*HOCO */
  ncf[112 * kd + 0] = 1; /*C */
  ncf[112 * kd + 2] = 2; /*O */
  ncf[112 * kd + 1] = 1; /*H */

  /*C2H5O */
  ncf[113 * kd + 0] = 2; /*C */
  ncf[113 * kd + 1] = 5; /*H */
  ncf[113 * kd + 2] = 1; /*O */

  /*C2H5O2 */
  ncf[114 * kd + 0] = 2; /*C */
  ncf[114 * kd + 1] = 5; /*H */
  ncf[114 * kd + 2] = 2; /*O */

  /*C2H5O2H */
  ncf[115 * kd + 0] = 2; /*C */
  ncf[115 * kd + 1] = 6; /*H */
  ncf[115 * kd + 2] = 2; /*O */

  /*C2H4O2H */
  ncf[116 * kd + 0] = 2; /*C */
  ncf[116 * kd + 1] = 5; /*H */
  ncf[116 * kd + 2] = 2; /*O */

  /*C2H5CO */
  ncf[117 * kd + 0] = 3; /*C */
  ncf[117 * kd + 1] = 5; /*H */
  ncf[117 * kd + 2] = 1; /*O */

  /*C2H5CHO */
  ncf[118 * kd + 0] = 3; /*C */
  ncf[118 * kd + 1] = 6; /*H */
  ncf[118 * kd + 2] = 1; /*O */

  /*CH3CO2 */
  ncf[119 * kd + 0] = 2; /*C */
  ncf[119 * kd + 1] = 3; /*H */
  ncf[119 * kd + 2] = 2; /*O */

  /*CH3CO3 */
  ncf[120 * kd + 0] = 2; /*C */
  ncf[120 * kd + 1] = 3; /*H */
  ncf[120 * kd + 2] = 3; /*O */

  /*CH3CO3H */
  ncf[121 * kd + 0] = 2; /*C */
  ncf[121 * kd + 1] = 4; /*H */
  ncf[121 * kd + 2] = 3; /*O */

  /*OCH2O2H */
  ncf[122 * kd + 0] = 1; /*C */
  ncf[122 * kd + 1] = 3; /*H */
  ncf[122 * kd + 2] = 3; /*O */

  /*HOCH2O2 */
  ncf[123 * kd + 0] = 1; /*C */
  ncf[123 * kd + 1] = 3; /*H */
  ncf[123 * kd + 2] = 3; /*O */

  /*HOCH2O2H */
  ncf[124 * kd + 0] = 1; /*C */
  ncf[124 * kd + 1] = 4; /*H */
  ncf[124 * kd + 2] = 3; /*O */

  /*O2CHO */
  ncf[125 * kd + 0] = 1; /*C */
  ncf[125 * kd + 1] = 1; /*H */
  ncf[125 * kd + 2] = 3; /*O */

  /*HO2CHO */
  ncf[126 * kd + 0] = 1; /*C */
  ncf[126 * kd + 1] = 2; /*H */
  ncf[126 * kd + 2] = 3; /*O */

  /*nC3H7 */
  ncf[127 * kd + 0] = 3; /*C */
  ncf[127 * kd + 1] = 7; /*H */

  /*iC3H7 */
  ncf[128 * kd + 0] = 3; /*C */
  ncf[128 * kd + 1] = 7; /*H */

  /*C3H8 */
  ncf[129 * kd + 0] = 3; /*C */
  ncf[129 * kd + 1] = 8; /*H */
}

/* Returns the vector of strings of element names */
void
CKSYME_STR(amrex::Vector<std::string>& ename)
{
  ename.resize(6);
  ename[0] = "C";
  ename[1] = "H";
  ename[2] = "O";
  ename[3] = "N";
  ename[4] = "AR";
  ename[5] = "HE";
}

/* Returns the vector of strings of species names */
void
CKSYMS_STR(amrex::Vector<std::string>& kname)
{
  kname.resize(130);
  kname[0] = "H";
  kname[1] = "H2";
  kname[2] = "CH2";
  kname[3] = "CH2(S)";
  kname[4] = "CH3";
  kname[5] = "O";
  kname[6] = "CH4";
  kname[7] = "OH";
  kname[8] = "H2O";
  kname[9] = "C2H";
  kname[10] = "C2H2";
  kname[11] = "C2H3";
  kname[12] = "CO";
  kname[13] = "N2";
  kname[14] = "C2H4";
  kname[15] = "HCO";
  kname[16] = "C2H5";
  kname[17] = "CH2O";
  kname[18] = "C2H6";
  kname[19] = "CH2OH";
  kname[20] = "CH3O";
  kname[21] = "O2";
  kname[22] = "CH3OH";
  kname[23] = "HO2";
  kname[24] = "H2O2";
  kname[25] = "HCCO";
  kname[26] = "CH2CO";
  kname[27] = "HCCOH";
  kname[28] = "CH2HCO";
  kname[29] = "CH3CO";
  kname[30] = "CO2";
  kname[31] = "CH3HCO";
  kname[32] = "OCHO";
  kname[33] = "CH3CHOH";
  kname[34] = "C2H4OH";
  kname[35] = "CH3CH2O";
  kname[36] = "CH3OCH2";
  kname[37] = "HOCHO";
  kname[38] = "CH3OCH3";
  kname[39] = "C2H5OH";
  kname[40] = "HOCH2O";
  kname[41] = "CH3OCO";
  kname[42] = "CH3OCHO";
  kname[43] = "CH3OCH2O";
  kname[44] = "CH3OCH2OH";
  kname[45] = "OCH2OCHO";
  kname[46] = "HOCH2OCO";
  kname[47] = "HOC2H4O2";
  kname[48] = "CH3OCH2O2";
  kname[49] = "CH2OCH2O2H";
  kname[50] = "CH3OCH2O2H";
  kname[51] = "HO2CH2OCHO";
  kname[52] = "O2CH2OCH2O2H";
  kname[53] = "AR";
  kname[54] = "HE";
  kname[55] = "C";
  kname[56] = "CH";
  kname[57] = "HCOH";
  kname[58] = "C2O";
  kname[59] = "CH3O2";
  kname[60] = "CH3O2H";
  kname[61] = "C3H2";
  kname[62] = "C3H2SING";
  kname[63] = "c-C3H2";
  kname[64] = "C3H3";
  kname[65] = "C3H6";
  kname[66] = "C4H2";
  kname[67] = "pC3H4";
  kname[68] = "aC3H4";
  kname[69] = "c-C3H4";
  kname[70] = "aC3H5";
  kname[71] = "sC3H5";
  kname[72] = "tC3H5";
  kname[73] = "nC4H3";
  kname[74] = "iC4H3";
  kname[75] = "C4H4";
  kname[76] = "nC4H5";
  kname[77] = "iC4H5";
  kname[78] = "C4H6-13";
  kname[79] = "C4H6-12";
  kname[80] = "C4H6-2";
  kname[81] = "H2CC";
  kname[82] = "CH2CHOH";
  kname[83] = "C5H2";
  kname[84] = "HCCCHCCH";
  kname[85] = "H2CCCCCH";
  kname[86] = "C4H5C2H";
  kname[87] = "FC6H6";
  kname[88] = "C6H6";
  kname[89] = "C6H5";
  kname[90] = "C2H3CHO";
  kname[91] = "C5H6";
  kname[92] = "C2H3CO";
  kname[93] = "C4H";
  kname[94] = "H2C4O";
  kname[95] = "C6H2";
  kname[96] = "C6H3";
  kname[97] = "l-C6H4";
  kname[98] = "o-C6H4";
  kname[99] = "CH3CHCHCO";
  kname[100] = "CH2CHCHCHO";
  kname[101] = "C5H5";
  kname[102] = "C6H5O";
  kname[103] = "C6H4O2";
  kname[104] = "C6H5OH";
  kname[105] = "C5H4O";
  kname[106] = "C5H5O";
  kname[107] = "OSING";
  kname[108] = "O2SING";
  kname[109] = "O3";
  kname[110] = "C2H4O2";
  kname[111] = "CH2OCHO";
  kname[112] = "HOCO";
  kname[113] = "C2H5O";
  kname[114] = "C2H5O2";
  kname[115] = "C2H5O2H";
  kname[116] = "C2H4O2H";
  kname[117] = "C2H5CO";
  kname[118] = "C2H5CHO";
  kname[119] = "CH3CO2";
  kname[120] = "CH3CO3";
  kname[121] = "CH3CO3H";
  kname[122] = "OCH2O2H";
  kname[123] = "HOCH2O2";
  kname[124] = "HOCH2O2H";
  kname[125] = "O2CHO";
  kname[126] = "HO2CHO";
  kname[127] = "nC3H7";
  kname[128] = "iC3H7";
  kname[129] = "C3H8";
}

/*compute the sparsity pattern of the chemistry Jacobian */
void
SPARSITY_INFO(int* nJdata, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 17161> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 130> conc = {0.0};
  for (int n = 0; n < 130; n++) {
    conc[n] = 1.0 / 130.000000;
  }
  aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 131; k++) {
    for (int l = 0; l < 131; l++) {
      if (Jac[131 * k + l] != 0.0) {
        nJdata_tmp = nJdata_tmp + 1;
      }
    }
  }

  *nJdata = NCELLS * nJdata_tmp;
}

/*compute the sparsity pattern of the system Jacobian */
void
SPARSITY_INFO_SYST(int* nJdata, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 17161> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 130> conc = {0.0};
  for (int n = 0; n < 130; n++) {
    conc[n] = 1.0 / 130.000000;
  }
  aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 131; k++) {
    for (int l = 0; l < 131; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[131 * k + l] != 0.0) {
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
    }
  }

  *nJdata = NCELLS * nJdata_tmp;
}

/*compute the sparsity pattern of the simplified (for preconditioning) system
 * Jacobian */
void
SPARSITY_INFO_SYST_SIMPLIFIED(int* nJdata, const int* consP)
{
  amrex::GpuArray<amrex::Real, 17161> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 130> conc = {0.0};
  for (int n = 0; n < 130; n++) {
    conc[n] = 1.0 / 130.000000;
  }
  aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

  int nJdata_tmp = 0;
  for (int k = 0; k < 131; k++) {
    for (int l = 0; l < 131; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[131 * k + l] != 0.0) {
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
    }
  }

  nJdata[0] = nJdata_tmp;
}

/*compute the sparsity pattern of the chemistry Jacobian in CSC format -- base 0
 */
void
SPARSITY_PREPROC_CSC(int* rowVals, int* colPtrs, const int* consP, int NCELLS)
{
  amrex::GpuArray<amrex::Real, 17161> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 130> conc = {0.0};
  for (int n = 0; n < 130; n++) {
    conc[n] = 1.0 / 130.000000;
  }
  aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int nc = 0; nc < NCELLS; nc++) {
    int offset_row = nc * 131;
    int offset_col = nc * 131;
    for (int k = 0; k < 131; k++) {
      for (int l = 0; l < 131; l++) {
        if (Jac[131 * k + l] != 0.0) {
          rowVals[nJdata_tmp] = l + offset_row;
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
      colPtrs[offset_col + (k + 1)] = nJdata_tmp;
    }
  }
}

/*compute the sparsity pattern of the chemistry Jacobian in CSR format -- base 0
 */
void
SPARSITY_PREPROC_CSR(
  int* colVals, int* rowPtrs, const int* consP, int NCELLS, int base)
{
  amrex::GpuArray<amrex::Real, 17161> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 130> conc = {0.0};
  for (int n = 0; n < 130; n++) {
    conc[n] = 1.0 / 130.000000;
  }
  aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

  if (base == 1) {
    rowPtrs[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 131;
      for (int l = 0; l < 131; l++) {
        for (int k = 0; k < 131; k++) {
          if (Jac[131 * k + l] != 0.0) {
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
      int offset = nc * 131;
      for (int l = 0; l < 131; l++) {
        for (int k = 0; k < 131; k++) {
          if (Jac[131 * k + l] != 0.0) {
            colVals[nJdata_tmp] = k + offset;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
        rowPtrs[offset + (l + 1)] = nJdata_tmp;
      }
    }
  }
}

/*compute the sparsity pattern of the system Jacobian */
/*CSR format BASE is user choice */
void
SPARSITY_PREPROC_SYST_CSR(
  int* colVals, int* rowPtr, const int* consP, int NCELLS, int base)
{
  amrex::GpuArray<amrex::Real, 17161> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 130> conc = {0.0};
  for (int n = 0; n < 130; n++) {
    conc[n] = 1.0 / 130.000000;
  }
  aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      int offset = nc * 131;
      for (int l = 0; l < 131; l++) {
        for (int k = 0; k < 131; k++) {
          if (k == l) {
            colVals[nJdata_tmp - 1] = l + 1 + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[131 * k + l] != 0.0) {
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
      int offset = nc * 131;
      for (int l = 0; l < 131; l++) {
        for (int k = 0; k < 131; k++) {
          if (k == l) {
            colVals[nJdata_tmp] = l + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (Jac[131 * k + l] != 0.0) {
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

/*compute the sparsity pattern of the simplified (for precond) system Jacobian
 * on CPU */
/*BASE 0 */
void
SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(
  int* rowVals, int* colPtrs, int* indx, const int* consP)
{
  amrex::GpuArray<amrex::Real, 17161> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 130> conc = {0.0};
  for (int n = 0; n < 130; n++) {
    conc[n] = 1.0 / 130.000000;
  }
  aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int k = 0; k < 131; k++) {
    for (int l = 0; l < 131; l++) {
      if (k == l) {
        rowVals[nJdata_tmp] = l;
        indx[nJdata_tmp] = 131 * k + l;
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (Jac[131 * k + l] != 0.0) {
          rowVals[nJdata_tmp] = l;
          indx[nJdata_tmp] = 131 * k + l;
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
    }
    colPtrs[k + 1] = nJdata_tmp;
  }
}

/*compute the sparsity pattern of the simplified (for precond) system Jacobian
 */
/*CSR format BASE is under choice */
void
SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(
  int* colVals, int* rowPtr, const int* consP, int base)
{
  amrex::GpuArray<amrex::Real, 17161> Jac = {0.0};
  amrex::GpuArray<amrex::Real, 130> conc = {0.0};
  for (int n = 0; n < 130; n++) {
    conc[n] = 1.0 / 130.000000;
  }
  aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int l = 0; l < 131; l++) {
      for (int k = 0; k < 131; k++) {
        if (k == l) {
          colVals[nJdata_tmp - 1] = l + 1;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[131 * k + l] != 0.0) {
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
    for (int l = 0; l < 131; l++) {
      for (int k = 0; k < 131; k++) {
        if (k == l) {
          colVals[nJdata_tmp] = l;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (Jac[131 * k + l] != 0.0) {
            colVals[nJdata_tmp] = k;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
      }
      rowPtr[l + 1] = nJdata_tmp;
    }
  }
}

#endif
