#ifndef MECHANISM_CPP
#define MECHANISM_CPP

#include "chemistry_file.H"

/*save atomic weights into array */
void
atomicWeight(amrex::Real* awt)
{
  awt[0] = 12.011150; /*C */
  awt[1] = 1.007970;  /*H */
  awt[2] = 15.999400; /*O */
  awt[3] = 14.006700; /*N */

  return;
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
  int kd = 4;
  /*Zero ncf */
  for (id = 0; id < kd * 143; ++id) {
    ncf[id] = 0;
  }

  /*H */
  ncf[0 * kd + 1] = 1; /*H */

  /*H2 */
  ncf[1 * kd + 1] = 2; /*H */

  /*O */
  ncf[2 * kd + 2] = 1; /*O */

  /*O2 */
  ncf[3 * kd + 2] = 2; /*O */

  /*OH */
  ncf[4 * kd + 1] = 1; /*H */
  ncf[4 * kd + 2] = 1; /*O */

  /*H2O */
  ncf[5 * kd + 1] = 2; /*H */
  ncf[5 * kd + 2] = 1; /*O */

  /*HO2 */
  ncf[6 * kd + 1] = 1; /*H */
  ncf[6 * kd + 2] = 2; /*O */

  /*H2O2 */
  ncf[7 * kd + 1] = 2; /*H */
  ncf[7 * kd + 2] = 2; /*O */

  /*CO */
  ncf[8 * kd + 0] = 1; /*C */
  ncf[8 * kd + 2] = 1; /*O */

  /*CO2 */
  ncf[9 * kd + 0] = 1; /*C */
  ncf[9 * kd + 2] = 2; /*O */

  /*CH2O */
  ncf[10 * kd + 0] = 1; /*C */
  ncf[10 * kd + 1] = 2; /*H */
  ncf[10 * kd + 2] = 1; /*O */

  /*HCO */
  ncf[11 * kd + 1] = 1; /*H */
  ncf[11 * kd + 0] = 1; /*C */
  ncf[11 * kd + 2] = 1; /*O */

  /*HO2CHO */
  ncf[12 * kd + 0] = 1; /*C */
  ncf[12 * kd + 1] = 2; /*H */
  ncf[12 * kd + 2] = 3; /*O */

  /*O2CHO */
  ncf[13 * kd + 0] = 1; /*C */
  ncf[13 * kd + 1] = 1; /*H */
  ncf[13 * kd + 2] = 3; /*O */

  /*OCHO */
  ncf[14 * kd + 0] = 1; /*C */
  ncf[14 * kd + 1] = 1; /*H */
  ncf[14 * kd + 2] = 2; /*O */

  /*CH2OH */
  ncf[15 * kd + 0] = 1; /*C */
  ncf[15 * kd + 1] = 3; /*H */
  ncf[15 * kd + 2] = 1; /*O */

  /*CH3O */
  ncf[16 * kd + 0] = 1; /*C */
  ncf[16 * kd + 1] = 3; /*H */
  ncf[16 * kd + 2] = 1; /*O */

  /*CH3O2H */
  ncf[17 * kd + 0] = 1; /*C */
  ncf[17 * kd + 1] = 4; /*H */
  ncf[17 * kd + 2] = 2; /*O */

  /*CH3O2 */
  ncf[18 * kd + 0] = 1; /*C */
  ncf[18 * kd + 1] = 3; /*H */
  ncf[18 * kd + 2] = 2; /*O */

  /*CH4 */
  ncf[19 * kd + 1] = 4; /*H */
  ncf[19 * kd + 0] = 1; /*C */

  /*CH3 */
  ncf[20 * kd + 0] = 1; /*C */
  ncf[20 * kd + 1] = 3; /*H */

  /*CH2 */
  ncf[21 * kd + 0] = 1; /*C */
  ncf[21 * kd + 1] = 2; /*H */

  /*CH2(S) */
  ncf[22 * kd + 0] = 1; /*C */
  ncf[22 * kd + 1] = 2; /*H */

  /*C2H6 */
  ncf[23 * kd + 0] = 2; /*C */
  ncf[23 * kd + 1] = 6; /*H */

  /*C2H5 */
  ncf[24 * kd + 0] = 2; /*C */
  ncf[24 * kd + 1] = 5; /*H */

  /*C2H4 */
  ncf[25 * kd + 0] = 2; /*C */
  ncf[25 * kd + 1] = 4; /*H */

  /*C2H3 */
  ncf[26 * kd + 0] = 2; /*C */
  ncf[26 * kd + 1] = 3; /*H */

  /*C2H2 */
  ncf[27 * kd + 0] = 2; /*C */
  ncf[27 * kd + 1] = 2; /*H */

  /*CH3CHO */
  ncf[28 * kd + 0] = 2; /*C */
  ncf[28 * kd + 1] = 4; /*H */
  ncf[28 * kd + 2] = 1; /*O */

  /*CH3CO */
  ncf[29 * kd + 0] = 2; /*C */
  ncf[29 * kd + 1] = 3; /*H */
  ncf[29 * kd + 2] = 1; /*O */

  /*CH2CO */
  ncf[30 * kd + 0] = 2; /*C */
  ncf[30 * kd + 1] = 2; /*H */
  ncf[30 * kd + 2] = 1; /*O */

  /*HCCO */
  ncf[31 * kd + 1] = 1; /*H */
  ncf[31 * kd + 0] = 2; /*C */
  ncf[31 * kd + 2] = 1; /*O */

  /*C2H5O */
  ncf[32 * kd + 0] = 2; /*C */
  ncf[32 * kd + 1] = 5; /*H */
  ncf[32 * kd + 2] = 1; /*O */

  /*C2H5O2 */
  ncf[33 * kd + 0] = 2; /*C */
  ncf[33 * kd + 1] = 5; /*H */
  ncf[33 * kd + 2] = 2; /*O */

  /*C2H3O1-2 */
  ncf[34 * kd + 0] = 2; /*C */
  ncf[34 * kd + 1] = 3; /*H */
  ncf[34 * kd + 2] = 1; /*O */

  /*CH3COCH3 */
  ncf[35 * kd + 0] = 3; /*C */
  ncf[35 * kd + 1] = 6; /*H */
  ncf[35 * kd + 2] = 1; /*O */

  /*CH3COCH2 */
  ncf[36 * kd + 0] = 3; /*C */
  ncf[36 * kd + 1] = 5; /*H */
  ncf[36 * kd + 2] = 1; /*O */

  /*C2H3CHO */
  ncf[37 * kd + 0] = 3; /*C */
  ncf[37 * kd + 1] = 4; /*H */
  ncf[37 * kd + 2] = 1; /*O */

  /*C2H3CO */
  ncf[38 * kd + 0] = 3; /*C */
  ncf[38 * kd + 1] = 3; /*H */
  ncf[38 * kd + 2] = 1; /*O */

  /*C2H5CO */
  ncf[39 * kd + 0] = 3; /*C */
  ncf[39 * kd + 1] = 5; /*H */
  ncf[39 * kd + 2] = 1; /*O */

  /*IC3H7 */
  ncf[40 * kd + 0] = 3; /*C */
  ncf[40 * kd + 1] = 7; /*H */

  /*C3H6 */
  ncf[41 * kd + 0] = 3; /*C */
  ncf[41 * kd + 1] = 6; /*H */

  /*C3H5-A */
  ncf[42 * kd + 0] = 3; /*C */
  ncf[42 * kd + 1] = 5; /*H */

  /*C3H5-S */
  ncf[43 * kd + 0] = 3; /*C */
  ncf[43 * kd + 1] = 5; /*H */

  /*C3H5-T */
  ncf[44 * kd + 0] = 3; /*C */
  ncf[44 * kd + 1] = 5; /*H */

  /*C3H4-P */
  ncf[45 * kd + 1] = 4; /*H */
  ncf[45 * kd + 0] = 3; /*C */

  /*C3H4-A */
  ncf[46 * kd + 1] = 4; /*H */
  ncf[46 * kd + 0] = 3; /*C */

  /*C3H3 */
  ncf[47 * kd + 0] = 3; /*C */
  ncf[47 * kd + 1] = 3; /*H */

  /*C3H2 */
  ncf[48 * kd + 1] = 2; /*H */
  ncf[48 * kd + 0] = 3; /*C */

  /*C3H5O */
  ncf[49 * kd + 0] = 3; /*C */
  ncf[49 * kd + 1] = 5; /*H */
  ncf[49 * kd + 2] = 1; /*O */

  /*C3H6OOH2-1 */
  ncf[50 * kd + 0] = 3; /*C */
  ncf[50 * kd + 1] = 7; /*H */
  ncf[50 * kd + 2] = 2; /*O */

  /*C3H6OOH2-1O2 */
  ncf[51 * kd + 0] = 3; /*C */
  ncf[51 * kd + 1] = 7; /*H */
  ncf[51 * kd + 2] = 4; /*O */

  /*IC3H7O2 */
  ncf[52 * kd + 0] = 3; /*C */
  ncf[52 * kd + 1] = 7; /*H */
  ncf[52 * kd + 2] = 2; /*O */

  /*C3KET21 */
  ncf[53 * kd + 0] = 3; /*C */
  ncf[53 * kd + 1] = 6; /*H */
  ncf[53 * kd + 2] = 3; /*O */

  /*CH3CHCO */
  ncf[54 * kd + 0] = 3; /*C */
  ncf[54 * kd + 1] = 4; /*H */
  ncf[54 * kd + 2] = 1; /*O */

  /*SC4H9 */
  ncf[55 * kd + 0] = 4; /*C */
  ncf[55 * kd + 1] = 9; /*H */

  /*IC4H9 */
  ncf[56 * kd + 0] = 4; /*C */
  ncf[56 * kd + 1] = 9; /*H */

  /*TC4H9 */
  ncf[57 * kd + 0] = 4; /*C */
  ncf[57 * kd + 1] = 9; /*H */

  /*IC4H8 */
  ncf[58 * kd + 0] = 4; /*C */
  ncf[58 * kd + 1] = 8; /*H */

  /*IC4H7 */
  ncf[59 * kd + 0] = 4; /*C */
  ncf[59 * kd + 1] = 7; /*H */

  /*TC4H9O2 */
  ncf[60 * kd + 0] = 4; /*C */
  ncf[60 * kd + 1] = 9; /*H */
  ncf[60 * kd + 2] = 2; /*O */

  /*IC4H9O2 */
  ncf[61 * kd + 0] = 4; /*C */
  ncf[61 * kd + 1] = 9; /*H */
  ncf[61 * kd + 2] = 2; /*O */

  /*IC4H8O2H-I */
  ncf[62 * kd + 0] = 4; /*C */
  ncf[62 * kd + 1] = 9; /*H */
  ncf[62 * kd + 2] = 2; /*O */

  /*TC4H9O */
  ncf[63 * kd + 0] = 4; /*C */
  ncf[63 * kd + 1] = 9; /*H */
  ncf[63 * kd + 2] = 1; /*O */

  /*TC4H9O2H */
  ncf[64 * kd + 0] = 4;  /*C */
  ncf[64 * kd + 1] = 10; /*H */
  ncf[64 * kd + 2] = 2;  /*O */

  /*IC4H7O */
  ncf[65 * kd + 0] = 4; /*C */
  ncf[65 * kd + 1] = 7; /*H */
  ncf[65 * kd + 2] = 1; /*O */

  /*IC3H7CHO */
  ncf[66 * kd + 0] = 4; /*C */
  ncf[66 * kd + 1] = 8; /*H */
  ncf[66 * kd + 2] = 1; /*O */

  /*TC3H6CHO */
  ncf[67 * kd + 0] = 4; /*C */
  ncf[67 * kd + 1] = 7; /*H */
  ncf[67 * kd + 2] = 1; /*O */

  /*IC4H8OOH-IO2 */
  ncf[68 * kd + 0] = 4; /*C */
  ncf[68 * kd + 1] = 9; /*H */
  ncf[68 * kd + 2] = 4; /*O */

  /*IC4KETII */
  ncf[69 * kd + 0] = 4; /*C */
  ncf[69 * kd + 1] = 8; /*H */
  ncf[69 * kd + 2] = 3; /*O */

  /*IC4H7OH */
  ncf[70 * kd + 0] = 4; /*C */
  ncf[70 * kd + 1] = 8; /*H */
  ncf[70 * kd + 2] = 1; /*O */

  /*IC4H6OH */
  ncf[71 * kd + 0] = 4; /*C */
  ncf[71 * kd + 1] = 7; /*H */
  ncf[71 * kd + 2] = 1; /*O */

  /*IC3H5CHO */
  ncf[72 * kd + 0] = 4; /*C */
  ncf[72 * kd + 1] = 6; /*H */
  ncf[72 * kd + 2] = 1; /*O */

  /*IC3H5CO */
  ncf[73 * kd + 0] = 4; /*C */
  ncf[73 * kd + 1] = 5; /*H */
  ncf[73 * kd + 2] = 1; /*O */

  /*TC3H6OCHO */
  ncf[74 * kd + 0] = 4; /*C */
  ncf[74 * kd + 1] = 7; /*H */
  ncf[74 * kd + 2] = 2; /*O */

  /*IC3H6CO */
  ncf[75 * kd + 0] = 4; /*C */
  ncf[75 * kd + 1] = 6; /*H */
  ncf[75 * kd + 2] = 1; /*O */

  /*IC4H7OOH */
  ncf[76 * kd + 0] = 4; /*C */
  ncf[76 * kd + 1] = 8; /*H */
  ncf[76 * kd + 2] = 2; /*O */

  /*TC3H6O2CHO */
  ncf[77 * kd + 0] = 4; /*C */
  ncf[77 * kd + 1] = 7; /*H */
  ncf[77 * kd + 2] = 3; /*O */

  /*CH2CCH2OH */
  ncf[78 * kd + 0] = 3; /*C */
  ncf[78 * kd + 1] = 5; /*H */
  ncf[78 * kd + 2] = 1; /*O */

  /*BC5H11 */
  ncf[79 * kd + 0] = 5;  /*C */
  ncf[79 * kd + 1] = 11; /*H */

  /*AC5H10 */
  ncf[80 * kd + 0] = 5;  /*C */
  ncf[80 * kd + 1] = 10; /*H */

  /*BC5H10 */
  ncf[81 * kd + 0] = 5;  /*C */
  ncf[81 * kd + 1] = 10; /*H */

  /*CC5H10 */
  ncf[82 * kd + 0] = 5;  /*C */
  ncf[82 * kd + 1] = 10; /*H */

  /*AC5H9-C */
  ncf[83 * kd + 0] = 5; /*C */
  ncf[83 * kd + 1] = 9; /*H */

  /*CC5H9-B */
  ncf[84 * kd + 0] = 5; /*C */
  ncf[84 * kd + 1] = 9; /*H */

  /*AC5H9O-C */
  ncf[85 * kd + 0] = 5; /*C */
  ncf[85 * kd + 1] = 9; /*H */
  ncf[85 * kd + 2] = 1; /*O */

  /*CC5H9O-B */
  ncf[86 * kd + 0] = 5; /*C */
  ncf[86 * kd + 1] = 9; /*H */
  ncf[86 * kd + 2] = 1; /*O */

  /*CH3CHCHO */
  ncf[87 * kd + 0] = 3; /*C */
  ncf[87 * kd + 1] = 5; /*H */
  ncf[87 * kd + 2] = 1; /*O */

  /*BC6H12 */
  ncf[88 * kd + 0] = 6;  /*C */
  ncf[88 * kd + 1] = 12; /*H */

  /*CC6H12 */
  ncf[89 * kd + 0] = 6;  /*C */
  ncf[89 * kd + 1] = 12; /*H */

  /*C5H10-2 */
  ncf[90 * kd + 0] = 5;  /*C */
  ncf[90 * kd + 1] = 10; /*H */

  /*IC4H7-I1 */
  ncf[91 * kd + 0] = 4; /*C */
  ncf[91 * kd + 1] = 7; /*H */

  /*YC7H15 */
  ncf[92 * kd + 0] = 7;  /*C */
  ncf[92 * kd + 1] = 15; /*H */

  /*XC7H14 */
  ncf[93 * kd + 0] = 7;  /*C */
  ncf[93 * kd + 1] = 14; /*H */

  /*YC7H14 */
  ncf[94 * kd + 0] = 7;  /*C */
  ncf[94 * kd + 1] = 14; /*H */

  /*XC7H13-Z */
  ncf[95 * kd + 0] = 7;  /*C */
  ncf[95 * kd + 1] = 13; /*H */

  /*YC7H13-Y2 */
  ncf[96 * kd + 0] = 7;  /*C */
  ncf[96 * kd + 1] = 13; /*H */

  /*YC7H13O-Y2 */
  ncf[97 * kd + 0] = 7;  /*C */
  ncf[97 * kd + 1] = 13; /*H */
  ncf[97 * kd + 2] = 1;  /*O */

  /*YC7H15O2 */
  ncf[98 * kd + 0] = 7;  /*C */
  ncf[98 * kd + 1] = 15; /*H */
  ncf[98 * kd + 2] = 2;  /*O */

  /*ACC6H10 */
  ncf[99 * kd + 0] = 6;  /*C */
  ncf[99 * kd + 1] = 10; /*H */

  /*ACC6H9-A */
  ncf[100 * kd + 0] = 6; /*C */
  ncf[100 * kd + 1] = 9; /*H */

  /*ACC6H9-D */
  ncf[101 * kd + 0] = 6; /*C */
  ncf[101 * kd + 1] = 9; /*H */

  /*NEOC5H11 */
  ncf[102 * kd + 0] = 5;  /*C */
  ncf[102 * kd + 1] = 11; /*H */

  /*NEOC5H11O2 */
  ncf[103 * kd + 0] = 5;  /*C */
  ncf[103 * kd + 1] = 11; /*H */
  ncf[103 * kd + 2] = 2;  /*O */

  /*NEOC5H10OOH */
  ncf[104 * kd + 0] = 5;  /*C */
  ncf[104 * kd + 1] = 11; /*H */
  ncf[104 * kd + 2] = 2;  /*O */

  /*TC4H9CHO */
  ncf[105 * kd + 0] = 5;  /*C */
  ncf[105 * kd + 1] = 10; /*H */
  ncf[105 * kd + 2] = 1;  /*O */

  /*TC4H9CO */
  ncf[106 * kd + 0] = 5; /*C */
  ncf[106 * kd + 1] = 9; /*H */
  ncf[106 * kd + 2] = 1; /*O */

  /*IC8H18 */
  ncf[107 * kd + 0] = 8;  /*C */
  ncf[107 * kd + 1] = 18; /*H */

  /*AC8H17 */
  ncf[108 * kd + 0] = 8;  /*C */
  ncf[108 * kd + 1] = 17; /*H */

  /*BC8H17 */
  ncf[109 * kd + 0] = 8;  /*C */
  ncf[109 * kd + 1] = 17; /*H */

  /*CC8H17 */
  ncf[110 * kd + 0] = 8;  /*C */
  ncf[110 * kd + 1] = 17; /*H */

  /*DC8H17 */
  ncf[111 * kd + 0] = 8;  /*C */
  ncf[111 * kd + 1] = 17; /*H */

  /*IC8H16 */
  ncf[112 * kd + 0] = 8;  /*C */
  ncf[112 * kd + 1] = 16; /*H */

  /*JC8H16 */
  ncf[113 * kd + 0] = 8;  /*C */
  ncf[113 * kd + 1] = 16; /*H */

  /*AC8H17O2 */
  ncf[114 * kd + 0] = 8;  /*C */
  ncf[114 * kd + 1] = 17; /*H */
  ncf[114 * kd + 2] = 2;  /*O */

  /*BC8H17O2 */
  ncf[115 * kd + 0] = 8;  /*C */
  ncf[115 * kd + 1] = 17; /*H */
  ncf[115 * kd + 2] = 2;  /*O */

  /*CC8H17O2 */
  ncf[116 * kd + 0] = 8;  /*C */
  ncf[116 * kd + 1] = 17; /*H */
  ncf[116 * kd + 2] = 2;  /*O */

  /*DC8H17O2 */
  ncf[117 * kd + 0] = 8;  /*C */
  ncf[117 * kd + 1] = 17; /*H */
  ncf[117 * kd + 2] = 2;  /*O */

  /*CC8H17O2H */
  ncf[118 * kd + 0] = 8;  /*C */
  ncf[118 * kd + 1] = 18; /*H */
  ncf[118 * kd + 2] = 2;  /*O */

  /*CC8H17O */
  ncf[119 * kd + 0] = 8;  /*C */
  ncf[119 * kd + 1] = 17; /*H */
  ncf[119 * kd + 2] = 1;  /*O */

  /*AC8H16OOH-A */
  ncf[120 * kd + 0] = 8;  /*C */
  ncf[120 * kd + 1] = 17; /*H */
  ncf[120 * kd + 2] = 2;  /*O */

  /*AC8H16OOH-B */
  ncf[121 * kd + 0] = 8;  /*C */
  ncf[121 * kd + 1] = 17; /*H */
  ncf[121 * kd + 2] = 2;  /*O */

  /*AC8H16OOH-C */
  ncf[122 * kd + 0] = 8;  /*C */
  ncf[122 * kd + 1] = 17; /*H */
  ncf[122 * kd + 2] = 2;  /*O */

  /*BC8H16OOH-A */
  ncf[123 * kd + 0] = 8;  /*C */
  ncf[123 * kd + 1] = 17; /*H */
  ncf[123 * kd + 2] = 2;  /*O */

  /*BC8H16OOH-D */
  ncf[124 * kd + 0] = 8;  /*C */
  ncf[124 * kd + 1] = 17; /*H */
  ncf[124 * kd + 2] = 2;  /*O */

  /*CC8H16OOH-A */
  ncf[125 * kd + 0] = 8;  /*C */
  ncf[125 * kd + 1] = 17; /*H */
  ncf[125 * kd + 2] = 2;  /*O */

  /*DC8H16OOH-C */
  ncf[126 * kd + 0] = 8;  /*C */
  ncf[126 * kd + 1] = 17; /*H */
  ncf[126 * kd + 2] = 2;  /*O */

  /*DC8H16OOH-B */
  ncf[127 * kd + 0] = 8;  /*C */
  ncf[127 * kd + 1] = 17; /*H */
  ncf[127 * kd + 2] = 2;  /*O */

  /*IC8ETERAB */
  ncf[128 * kd + 0] = 8;  /*C */
  ncf[128 * kd + 1] = 16; /*H */
  ncf[128 * kd + 2] = 1;  /*O */

  /*IC8ETERAC */
  ncf[129 * kd + 0] = 8;  /*C */
  ncf[129 * kd + 1] = 16; /*H */
  ncf[129 * kd + 2] = 1;  /*O */

  /*IC8ETERBD */
  ncf[130 * kd + 0] = 8;  /*C */
  ncf[130 * kd + 1] = 16; /*H */
  ncf[130 * kd + 2] = 1;  /*O */

  /*AC8H16OOH-BO2 */
  ncf[131 * kd + 0] = 8;  /*C */
  ncf[131 * kd + 1] = 17; /*H */
  ncf[131 * kd + 2] = 4;  /*O */

  /*BC8H16OOH-AO2 */
  ncf[132 * kd + 0] = 8;  /*C */
  ncf[132 * kd + 1] = 17; /*H */
  ncf[132 * kd + 2] = 4;  /*O */

  /*BC8H16OOH-DO2 */
  ncf[133 * kd + 0] = 8;  /*C */
  ncf[133 * kd + 1] = 17; /*H */
  ncf[133 * kd + 2] = 4;  /*O */

  /*CC8H16OOH-AO2 */
  ncf[134 * kd + 0] = 8;  /*C */
  ncf[134 * kd + 1] = 17; /*H */
  ncf[134 * kd + 2] = 4;  /*O */

  /*DC8H16OOH-BO2 */
  ncf[135 * kd + 0] = 8;  /*C */
  ncf[135 * kd + 1] = 17; /*H */
  ncf[135 * kd + 2] = 4;  /*O */

  /*IC8KETAB */
  ncf[136 * kd + 0] = 8;  /*C */
  ncf[136 * kd + 1] = 16; /*H */
  ncf[136 * kd + 2] = 3;  /*O */

  /*IC8KETBA */
  ncf[137 * kd + 0] = 8;  /*C */
  ncf[137 * kd + 1] = 16; /*H */
  ncf[137 * kd + 2] = 3;  /*O */

  /*IC8KETBD */
  ncf[138 * kd + 0] = 8;  /*C */
  ncf[138 * kd + 1] = 16; /*H */
  ncf[138 * kd + 2] = 3;  /*O */

  /*IC8KETDB */
  ncf[139 * kd + 0] = 8;  /*C */
  ncf[139 * kd + 1] = 16; /*H */
  ncf[139 * kd + 2] = 3;  /*O */

  /*IC3H7COC3H6-T */
  ncf[140 * kd + 0] = 7;  /*C */
  ncf[140 * kd + 1] = 13; /*H */
  ncf[140 * kd + 2] = 1;  /*O */

  /*TC4H9COC2H4S */
  ncf[141 * kd + 0] = 7;  /*C */
  ncf[141 * kd + 1] = 13; /*H */
  ncf[141 * kd + 2] = 1;  /*O */

  /*N2 */
  ncf[142 * kd + 3] = 2; /*N */
}

/* Returns the vector of strings of element names */
void
CKSYME_STR(amrex::Vector<std::string>& ename)
{
  ename.resize(4);
  ename[0] = "C";
  ename[1] = "H";
  ename[2] = "O";
  ename[3] = "N";
}

/* Returns the vector of strings of species names */
void
CKSYMS_STR(amrex::Vector<std::string>& kname)
{
  kname.resize(143);
  kname[0] = "H";
  kname[1] = "H2";
  kname[2] = "O";
  kname[3] = "O2";
  kname[4] = "OH";
  kname[5] = "H2O";
  kname[6] = "HO2";
  kname[7] = "H2O2";
  kname[8] = "CO";
  kname[9] = "CO2";
  kname[10] = "CH2O";
  kname[11] = "HCO";
  kname[12] = "HO2CHO";
  kname[13] = "O2CHO";
  kname[14] = "OCHO";
  kname[15] = "CH2OH";
  kname[16] = "CH3O";
  kname[17] = "CH3O2H";
  kname[18] = "CH3O2";
  kname[19] = "CH4";
  kname[20] = "CH3";
  kname[21] = "CH2";
  kname[22] = "CH2(S)";
  kname[23] = "C2H6";
  kname[24] = "C2H5";
  kname[25] = "C2H4";
  kname[26] = "C2H3";
  kname[27] = "C2H2";
  kname[28] = "CH3CHO";
  kname[29] = "CH3CO";
  kname[30] = "CH2CO";
  kname[31] = "HCCO";
  kname[32] = "C2H5O";
  kname[33] = "C2H5O2";
  kname[34] = "C2H3O1-2";
  kname[35] = "CH3COCH3";
  kname[36] = "CH3COCH2";
  kname[37] = "C2H3CHO";
  kname[38] = "C2H3CO";
  kname[39] = "C2H5CO";
  kname[40] = "IC3H7";
  kname[41] = "C3H6";
  kname[42] = "C3H5-A";
  kname[43] = "C3H5-S";
  kname[44] = "C3H5-T";
  kname[45] = "C3H4-P";
  kname[46] = "C3H4-A";
  kname[47] = "C3H3";
  kname[48] = "C3H2";
  kname[49] = "C3H5O";
  kname[50] = "C3H6OOH2-1";
  kname[51] = "C3H6OOH2-1O2";
  kname[52] = "IC3H7O2";
  kname[53] = "C3KET21";
  kname[54] = "CH3CHCO";
  kname[55] = "SC4H9";
  kname[56] = "IC4H9";
  kname[57] = "TC4H9";
  kname[58] = "IC4H8";
  kname[59] = "IC4H7";
  kname[60] = "TC4H9O2";
  kname[61] = "IC4H9O2";
  kname[62] = "IC4H8O2H-I";
  kname[63] = "TC4H9O";
  kname[64] = "TC4H9O2H";
  kname[65] = "IC4H7O";
  kname[66] = "IC3H7CHO";
  kname[67] = "TC3H6CHO";
  kname[68] = "IC4H8OOH-IO2";
  kname[69] = "IC4KETII";
  kname[70] = "IC4H7OH";
  kname[71] = "IC4H6OH";
  kname[72] = "IC3H5CHO";
  kname[73] = "IC3H5CO";
  kname[74] = "TC3H6OCHO";
  kname[75] = "IC3H6CO";
  kname[76] = "IC4H7OOH";
  kname[77] = "TC3H6O2CHO";
  kname[78] = "CH2CCH2OH";
  kname[79] = "BC5H11";
  kname[80] = "AC5H10";
  kname[81] = "BC5H10";
  kname[82] = "CC5H10";
  kname[83] = "AC5H9-C";
  kname[84] = "CC5H9-B";
  kname[85] = "AC5H9O-C";
  kname[86] = "CC5H9O-B";
  kname[87] = "CH3CHCHO";
  kname[88] = "BC6H12";
  kname[89] = "CC6H12";
  kname[90] = "C5H10-2";
  kname[91] = "IC4H7-I1";
  kname[92] = "YC7H15";
  kname[93] = "XC7H14";
  kname[94] = "YC7H14";
  kname[95] = "XC7H13-Z";
  kname[96] = "YC7H13-Y2";
  kname[97] = "YC7H13O-Y2";
  kname[98] = "YC7H15O2";
  kname[99] = "ACC6H10";
  kname[100] = "ACC6H9-A";
  kname[101] = "ACC6H9-D";
  kname[102] = "NEOC5H11";
  kname[103] = "NEOC5H11O2";
  kname[104] = "NEOC5H10OOH";
  kname[105] = "TC4H9CHO";
  kname[106] = "TC4H9CO";
  kname[107] = "IC8H18";
  kname[108] = "AC8H17";
  kname[109] = "BC8H17";
  kname[110] = "CC8H17";
  kname[111] = "DC8H17";
  kname[112] = "IC8H16";
  kname[113] = "JC8H16";
  kname[114] = "AC8H17O2";
  kname[115] = "BC8H17O2";
  kname[116] = "CC8H17O2";
  kname[117] = "DC8H17O2";
  kname[118] = "CC8H17O2H";
  kname[119] = "CC8H17O";
  kname[120] = "AC8H16OOH-A";
  kname[121] = "AC8H16OOH-B";
  kname[122] = "AC8H16OOH-C";
  kname[123] = "BC8H16OOH-A";
  kname[124] = "BC8H16OOH-D";
  kname[125] = "CC8H16OOH-A";
  kname[126] = "DC8H16OOH-C";
  kname[127] = "DC8H16OOH-B";
  kname[128] = "IC8ETERAB";
  kname[129] = "IC8ETERAC";
  kname[130] = "IC8ETERBD";
  kname[131] = "AC8H16OOH-BO2";
  kname[132] = "BC8H16OOH-AO2";
  kname[133] = "BC8H16OOH-DO2";
  kname[134] = "CC8H16OOH-AO2";
  kname[135] = "DC8H16OOH-BO2";
  kname[136] = "IC8KETAB";
  kname[137] = "IC8KETBA";
  kname[138] = "IC8KETBD";
  kname[139] = "IC8KETDB";
  kname[140] = "IC3H7COC3H6-T";
  kname[141] = "TC4H9COC2H4S";
  kname[142] = "N2";
}

/*compute the sparsity pattern of the chemistry Jacobian */
void
SPARSITY_INFO(int* nJdata, int* consP, int NCELLS)
{
  amrex::Gpu::DeviceVector<amrex::Real> J_v(20736);
  amrex::Gpu::DeviceVector<amrex::Real> c_v(143);
  amrex::Real* J_d = J_v.data();
  amrex::Real* c_d = c_v.data();

  amrex::Real J_h[20736];

  amrex::IntVect iv(AMREX_D_DECL(0, 0, 0));
  amrex::ParallelFor(
    amrex::Box(iv, iv),
    [=] AMREX_GPU_HOST_DEVICE(int /*i*/, int /*j*/, int /*k*/) noexcept {
      for (int l = 0; l < 143; l++) {
        c_d[l] = 1.0 / 143.000000;
      }
      aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
  amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
  std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

  int nJdata_tmp = 0;
  for (int k = 0; k < 144; k++) {
    for (int l = 0; l < 144; l++) {
      if (J_h[144 * k + l] != 0.0) {
        nJdata_tmp = nJdata_tmp + 1;
      }
    }
  }

  *nJdata = NCELLS * nJdata_tmp;

  return;
}

/*compute the sparsity pattern of the system Jacobian */
void
SPARSITY_INFO_SYST(int* nJdata, int* consP, int NCELLS)
{
  amrex::Gpu::DeviceVector<amrex::Real> J_v(20736);
  amrex::Gpu::DeviceVector<amrex::Real> c_v(143);
  amrex::Real* J_d = J_v.data();
  amrex::Real* c_d = c_v.data();

  amrex::Real J_h[20736];

  amrex::IntVect iv(AMREX_D_DECL(0, 0, 0));
  amrex::ParallelFor(
    amrex::Box(iv, iv),
    [=] AMREX_GPU_HOST_DEVICE(int /*i*/, int /*j*/, int /*k*/) noexcept {
      for (int k = 0; k < 143; k++) {
        c_d[k] = 1.0 / 143.000000;
      }
      aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
  amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
  std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

  int nJdata_tmp = 0;
  for (int k = 0; k < 144; k++) {
    for (int l = 0; l < 144; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (J_h[144 * k + l] != 0.0) {
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
    }
  }

  *nJdata = NCELLS * nJdata_tmp;

  return;
}

/*compute the sparsity pattern of the simplified (for preconditioning) system
 * Jacobian */
void
SPARSITY_INFO_SYST_SIMPLIFIED(int* nJdata, int* consP)
{
  amrex::Gpu::DeviceVector<amrex::Real> J_v(20736);
  amrex::Gpu::DeviceVector<amrex::Real> c_v(143);
  amrex::Real* J_d = J_v.data();
  amrex::Real* c_d = c_v.data();

  amrex::Real J_h[20736];

  amrex::IntVect iv(AMREX_D_DECL(0, 0, 0));
  amrex::ParallelFor(
    amrex::Box(iv, iv),
    [=] AMREX_GPU_HOST_DEVICE(int /*i*/, int /*j*/, int /*k*/) noexcept {
      for (int k = 0; k < 143; k++) {
        c_d[k] = 1.0 / 143.000000;
      }
      aJacobian_precond(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
  amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
  std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

  int nJdata_tmp = 0;
  for (int k = 0; k < 144; k++) {
    for (int l = 0; l < 144; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (J_h[144 * k + l] != 0.0) {
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
    }
  }

  nJdata[0] = nJdata_tmp;

  return;
}

/*compute the sparsity pattern of the chemistry Jacobian in CSC format -- base 0
 */
void
SPARSITY_PREPROC_CSC(int* rowVals, int* colPtrs, int* consP, int NCELLS)
{
  int offset_row;
  int offset_col;

  amrex::Gpu::DeviceVector<amrex::Real> J_v(20736);
  amrex::Gpu::DeviceVector<amrex::Real> c_v(143);
  amrex::Real* J_d = J_v.data();
  amrex::Real* c_d = c_v.data();

  amrex::Real J_h[20736];

  amrex::IntVect iv(AMREX_D_DECL(0, 0, 0));
  amrex::ParallelFor(
    amrex::Box(iv, iv),
    [=] AMREX_GPU_HOST_DEVICE(int /*i*/, int /*j*/, int /*k*/) noexcept {
      for (int k = 0; k < 143; k++) {
        c_d[k] = 1.0 / 143.000000;
      }
      aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
  amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
  std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int nc = 0; nc < NCELLS; nc++) {
    offset_row = nc * 144;
    offset_col = nc * 144;
    for (int k = 0; k < 144; k++) {
      for (int l = 0; l < 144; l++) {
        if (J_h[144 * k + l] != 0.0) {
          rowVals[nJdata_tmp] = l + offset_row;
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
      colPtrs[offset_col + (k + 1)] = nJdata_tmp;
    }
  }

  return;
}

/*compute the sparsity pattern of the chemistry Jacobian in CSR format -- base 0
 */
void
SPARSITY_PREPROC_CSR(
  int* colVals, int* rowPtrs, int* consP, int NCELLS, int base)
{
  int offset;
  amrex::Gpu::DeviceVector<amrex::Real> J_v(20736);
  amrex::Gpu::DeviceVector<amrex::Real> c_v(143);
  amrex::Real* J_d = J_v.data();
  amrex::Real* c_d = c_v.data();

  amrex::Real J_h[20736];

  amrex::IntVect iv(AMREX_D_DECL(0, 0, 0));
  amrex::ParallelFor(
    amrex::Box(iv, iv),
    [=] AMREX_GPU_HOST_DEVICE(int /*i*/, int /*j*/, int /*k*/) noexcept {
      for (int k = 0; k < 143; k++) {
        c_d[k] = 1.0 / 143.000000;
      }
      aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
  amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
  std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

  if (base == 1) {
    rowPtrs[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      offset = nc * 144;
      for (int l = 0; l < 144; l++) {
        for (int k = 0; k < 144; k++) {
          if (J_h[144 * k + l] != 0.0) {
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
      offset = nc * 144;
      for (int l = 0; l < 144; l++) {
        for (int k = 0; k < 144; k++) {
          if (J_h[144 * k + l] != 0.0) {
            colVals[nJdata_tmp] = k + offset;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
        rowPtrs[offset + (l + 1)] = nJdata_tmp;
      }
    }
  }

  return;
}

/*compute the sparsity pattern of the system Jacobian */
/*CSR format BASE is user choice */
void
SPARSITY_PREPROC_SYST_CSR(
  int* colVals, int* rowPtr, int* consP, int NCELLS, int base)
{
  int offset;
  amrex::Gpu::DeviceVector<amrex::Real> J_v(20736);
  amrex::Gpu::DeviceVector<amrex::Real> c_v(143);
  amrex::Real* J_d = J_v.data();
  amrex::Real* c_d = c_v.data();

  amrex::Real J_h[20736];

  amrex::IntVect iv(AMREX_D_DECL(0, 0, 0));
  amrex::ParallelFor(
    amrex::Box(iv, iv),
    [=] AMREX_GPU_HOST_DEVICE(int /*i*/, int /*j*/, int /*k*/) noexcept {
      for (int k = 0; k < 143; k++) {
        c_d[k] = 1.0 / 143.000000;
      }
      aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
  amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
  std::memcpy(&J_h, J_d, sizeof(J_h));
#endif
  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int nc = 0; nc < NCELLS; nc++) {
      offset = nc * 144;
      for (int l = 0; l < 144; l++) {
        for (int k = 0; k < 144; k++) {
          if (k == l) {
            colVals[nJdata_tmp - 1] = l + 1 + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (J_h[144 * k + l] != 0.0) {
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
      offset = nc * 144;
      for (int l = 0; l < 144; l++) {
        for (int k = 0; k < 144; k++) {
          if (k == l) {
            colVals[nJdata_tmp] = l + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (J_h[144 * k + l] != 0.0) {
              colVals[nJdata_tmp] = k + offset;
              nJdata_tmp = nJdata_tmp + 1;
            }
          }
        }
        rowPtr[offset + (l + 1)] = nJdata_tmp;
      }
    }
  }

  return;
}

/*compute the sparsity pattern of the simplified (for precond) system Jacobian
 * on CPU */
/*BASE 0 */
void
SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(
  int* rowVals, int* colPtrs, int* indx, int* consP)
{
  amrex::Gpu::DeviceVector<amrex::Real> J_v(20736);
  amrex::Gpu::DeviceVector<amrex::Real> c_v(143);
  amrex::Real* J_d = J_v.data();
  amrex::Real* c_d = c_v.data();

  amrex::Real J_h[20736];

  amrex::IntVect iv(AMREX_D_DECL(0, 0, 0));
  amrex::ParallelFor(
    amrex::Box(iv, iv),
    [=] AMREX_GPU_HOST_DEVICE(int /*i*/, int /*j*/, int /*k*/) noexcept {
      for (int k = 0; k < 143; k++) {
        c_d[k] = 1.0 / 143.000000;
      }
      aJacobian_precond(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
  amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
  std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

  colPtrs[0] = 0;
  int nJdata_tmp = 0;
  for (int k = 0; k < 144; k++) {
    for (int l = 0; l < 144; l++) {
      if (k == l) {
        rowVals[nJdata_tmp] = l;
        indx[nJdata_tmp] = 144 * k + l;
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (J_h[144 * k + l] != 0.0) {
          rowVals[nJdata_tmp] = l;
          indx[nJdata_tmp] = 144 * k + l;
          nJdata_tmp = nJdata_tmp + 1;
        }
      }
    }
    colPtrs[k + 1] = nJdata_tmp;
  }

  return;
}

/*compute the sparsity pattern of the simplified (for precond) system Jacobian
 */
/*CSR format BASE is under choice */
void
SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(
  int* colVals, int* rowPtr, int* consP, int base)
{
  amrex::Gpu::DeviceVector<amrex::Real> J_v(20736);
  amrex::Gpu::DeviceVector<amrex::Real> c_v(143);
  amrex::Real* J_d = J_v.data();
  amrex::Real* c_d = c_v.data();

  amrex::Real J_h[20736];

  amrex::IntVect iv(AMREX_D_DECL(0, 0, 0));
  amrex::ParallelFor(
    amrex::Box(iv, iv),
    [=] AMREX_GPU_HOST_DEVICE(int /*i*/, int /*j*/, int /*k*/) noexcept {
      for (int k = 0; k < 143; k++) {
        c_d[k] = 1.0 / 143.000000;
      }
      aJacobian_precond(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
  amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
  std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

  if (base == 1) {
    rowPtr[0] = 1;
    int nJdata_tmp = 1;
    for (int l = 0; l < 144; l++) {
      for (int k = 0; k < 144; k++) {
        if (k == l) {
          colVals[nJdata_tmp - 1] = l + 1;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (J_h[144 * k + l] != 0.0) {
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
    for (int l = 0; l < 144; l++) {
      for (int k = 0; k < 144; k++) {
        if (k == l) {
          colVals[nJdata_tmp] = l;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (J_h[144 * k + l] != 0.0) {
            colVals[nJdata_tmp] = k;
            nJdata_tmp = nJdata_tmp + 1;
          }
        }
      }
      rowPtr[l + 1] = nJdata_tmp;
    }
  }

  return;
}

#endif
