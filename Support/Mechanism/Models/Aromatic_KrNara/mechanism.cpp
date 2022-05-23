#ifndef MECHANISM_CPP
#define MECHANISM_CPP

#include "chemistry_file.H"

/*save atomic weights into array */
void
atomicWeight(amrex::Real* awt)
{
  awt[0] = 14.006700; /*N */
  awt[1] = 1.007970;  /*H */
  awt[2] = 15.999400; /*O */
  awt[3] = 12.011150; /*C */
  awt[4] = 39.948000; /*AR */

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
  int kd = 5;
  /*Zero ncf */
  for (id = 0; id < kd * 158; ++id) {
    ncf[id] = 0;
  }

  /*H */
  ncf[0 * kd + 1] = 1; /*H */

  /*O2 */
  ncf[1 * kd + 2] = 2; /*O */

  /*O */
  ncf[2 * kd + 2] = 1; /*O */

  /*OH */
  ncf[3 * kd + 2] = 1; /*O */
  ncf[3 * kd + 1] = 1; /*H */

  /*H2 */
  ncf[4 * kd + 1] = 2; /*H */

  /*H2O */
  ncf[5 * kd + 1] = 2; /*H */
  ncf[5 * kd + 2] = 1; /*O */

  /*CO2 */
  ncf[6 * kd + 3] = 1; /*C */
  ncf[6 * kd + 2] = 2; /*O */

  /*HO2 */
  ncf[7 * kd + 1] = 1; /*H */
  ncf[7 * kd + 2] = 2; /*O */

  /*H2O2 */
  ncf[8 * kd + 1] = 2; /*H */
  ncf[8 * kd + 2] = 2; /*O */

  /*CO */
  ncf[9 * kd + 3] = 1; /*C */
  ncf[9 * kd + 2] = 1; /*O */

  /*HCO */
  ncf[10 * kd + 1] = 1; /*H */
  ncf[10 * kd + 3] = 1; /*C */
  ncf[10 * kd + 2] = 1; /*O */

  /*C */
  ncf[11 * kd + 3] = 1; /*C */

  /*CH */
  ncf[12 * kd + 3] = 1; /*C */
  ncf[12 * kd + 1] = 1; /*H */

  /*T-CH2 */
  ncf[13 * kd + 3] = 1; /*C */
  ncf[13 * kd + 1] = 2; /*H */

  /*CH3 */
  ncf[14 * kd + 3] = 1; /*C */
  ncf[14 * kd + 1] = 3; /*H */

  /*CH2O */
  ncf[15 * kd + 1] = 2; /*H */
  ncf[15 * kd + 3] = 1; /*C */
  ncf[15 * kd + 2] = 1; /*O */

  /*HCCO */
  ncf[16 * kd + 1] = 1; /*H */
  ncf[16 * kd + 3] = 2; /*C */
  ncf[16 * kd + 2] = 1; /*O */

  /*C2H */
  ncf[17 * kd + 3] = 2; /*C */
  ncf[17 * kd + 1] = 1; /*H */

  /*CH2CO */
  ncf[18 * kd + 3] = 2; /*C */
  ncf[18 * kd + 1] = 2; /*H */
  ncf[18 * kd + 2] = 1; /*O */

  /*C2H2 */
  ncf[19 * kd + 3] = 2; /*C */
  ncf[19 * kd + 1] = 2; /*H */

  /*S-CH2 */
  ncf[20 * kd + 3] = 1; /*C */
  ncf[20 * kd + 1] = 2; /*H */

  /*AR */
  ncf[21 * kd + 4] = 1; /*AR */

  /*CH3OH */
  ncf[22 * kd + 3] = 1; /*C */
  ncf[22 * kd + 1] = 4; /*H */
  ncf[22 * kd + 2] = 1; /*O */

  /*CH2OH */
  ncf[23 * kd + 3] = 1; /*C */
  ncf[23 * kd + 1] = 3; /*H */
  ncf[23 * kd + 2] = 1; /*O */

  /*CH3O */
  ncf[24 * kd + 3] = 1; /*C */
  ncf[24 * kd + 1] = 3; /*H */
  ncf[24 * kd + 2] = 1; /*O */

  /*CH4 */
  ncf[25 * kd + 3] = 1; /*C */
  ncf[25 * kd + 1] = 4; /*H */

  /*CH3O2 */
  ncf[26 * kd + 3] = 1; /*C */
  ncf[26 * kd + 1] = 3; /*H */
  ncf[26 * kd + 2] = 2; /*O */

  /*C2H3 */
  ncf[27 * kd + 3] = 2; /*C */
  ncf[27 * kd + 1] = 3; /*H */

  /*C2H4 */
  ncf[28 * kd + 3] = 2; /*C */
  ncf[28 * kd + 1] = 4; /*H */

  /*C2H5 */
  ncf[29 * kd + 3] = 2; /*C */
  ncf[29 * kd + 1] = 5; /*H */

  /*HCCOH */
  ncf[30 * kd + 3] = 2; /*C */
  ncf[30 * kd + 2] = 1; /*O */
  ncf[30 * kd + 1] = 2; /*H */

  /*CH2CHO */
  ncf[31 * kd + 1] = 3; /*H */
  ncf[31 * kd + 3] = 2; /*C */
  ncf[31 * kd + 2] = 1; /*O */

  /*CH3CHO */
  ncf[32 * kd + 1] = 4; /*H */
  ncf[32 * kd + 3] = 2; /*C */
  ncf[32 * kd + 2] = 1; /*O */

  /*H2C2 */
  ncf[33 * kd + 1] = 2; /*H */
  ncf[33 * kd + 3] = 2; /*C */

  /*C2H5O */
  ncf[34 * kd + 3] = 2; /*C */
  ncf[34 * kd + 1] = 5; /*H */
  ncf[34 * kd + 2] = 1; /*O */

  /*N-C3H7 */
  ncf[35 * kd + 3] = 3; /*C */
  ncf[35 * kd + 1] = 7; /*H */

  /*C2H6 */
  ncf[36 * kd + 3] = 2; /*C */
  ncf[36 * kd + 1] = 6; /*H */

  /*C3H8 */
  ncf[37 * kd + 3] = 3; /*C */
  ncf[37 * kd + 1] = 8; /*H */

  /*C3H6 */
  ncf[38 * kd + 1] = 6; /*H */
  ncf[38 * kd + 3] = 3; /*C */

  /*C3H3 */
  ncf[39 * kd + 1] = 3; /*H */
  ncf[39 * kd + 3] = 3; /*C */

  /*P-C3H4 */
  ncf[40 * kd + 1] = 4; /*H */
  ncf[40 * kd + 3] = 3; /*C */

  /*A-C3H4 */
  ncf[41 * kd + 1] = 4; /*H */
  ncf[41 * kd + 3] = 3; /*C */

  /*S-C3H5 */
  ncf[42 * kd + 1] = 5; /*H */
  ncf[42 * kd + 3] = 3; /*C */

  /*N-C4H3 */
  ncf[43 * kd + 1] = 3; /*H */
  ncf[43 * kd + 3] = 4; /*C */

  /*C2H3CHO */
  ncf[44 * kd + 3] = 3; /*C */
  ncf[44 * kd + 1] = 4; /*H */
  ncf[44 * kd + 2] = 1; /*O */

  /*A-C3H5 */
  ncf[45 * kd + 1] = 5; /*H */
  ncf[45 * kd + 3] = 3; /*C */

  /*C2O */
  ncf[46 * kd + 3] = 2; /*C */
  ncf[46 * kd + 2] = 1; /*O */

  /*C4H4 */
  ncf[47 * kd + 1] = 4; /*H */
  ncf[47 * kd + 3] = 4; /*C */

  /*C3H2 */
  ncf[48 * kd + 1] = 2; /*H */
  ncf[48 * kd + 3] = 3; /*C */

  /*C3H2O */
  ncf[49 * kd + 1] = 2; /*H */
  ncf[49 * kd + 3] = 3; /*C */
  ncf[49 * kd + 2] = 1; /*O */

  /*C4H2 */
  ncf[50 * kd + 1] = 2; /*H */
  ncf[50 * kd + 3] = 4; /*C */

  /*I-C4H3 */
  ncf[51 * kd + 1] = 3; /*H */
  ncf[51 * kd + 3] = 4; /*C */

  /*T-C3H5 */
  ncf[52 * kd + 1] = 5; /*H */
  ncf[52 * kd + 3] = 3; /*C */

  /*C3H5O */
  ncf[53 * kd + 3] = 3; /*C */
  ncf[53 * kd + 1] = 5; /*H */
  ncf[53 * kd + 2] = 1; /*O */

  /*C4H */
  ncf[54 * kd + 1] = 1; /*H */
  ncf[54 * kd + 3] = 4; /*C */

  /*C8H2 */
  ncf[55 * kd + 3] = 8; /*C */
  ncf[55 * kd + 1] = 2; /*H */

  /*C6H2 */
  ncf[56 * kd + 3] = 6; /*C */
  ncf[56 * kd + 1] = 2; /*H */

  /*C4H6 */
  ncf[57 * kd + 1] = 6; /*H */
  ncf[57 * kd + 3] = 4; /*C */

  /*N-C4H5 */
  ncf[58 * kd + 1] = 5; /*H */
  ncf[58 * kd + 3] = 4; /*C */

  /*I-C4H5 */
  ncf[59 * kd + 1] = 5; /*H */
  ncf[59 * kd + 3] = 4; /*C */

  /*A1 */
  ncf[60 * kd + 1] = 6; /*H */
  ncf[60 * kd + 3] = 6; /*C */

  /*N-C7H16 */
  ncf[61 * kd + 3] = 7;  /*C */
  ncf[61 * kd + 1] = 16; /*H */

  /*C5H11 */
  ncf[62 * kd + 3] = 5;  /*C */
  ncf[62 * kd + 1] = 11; /*H */

  /*P-C4H9 */
  ncf[63 * kd + 3] = 4; /*C */
  ncf[63 * kd + 1] = 9; /*H */

  /*C7H15 */
  ncf[64 * kd + 3] = 7;  /*C */
  ncf[64 * kd + 1] = 15; /*H */

  /*P-C4H8 */
  ncf[65 * kd + 3] = 4; /*C */
  ncf[65 * kd + 1] = 8; /*H */

  /*C5H10 */
  ncf[66 * kd + 3] = 5;  /*C */
  ncf[66 * kd + 1] = 10; /*H */

  /*C7H14 */
  ncf[67 * kd + 3] = 7;  /*C */
  ncf[67 * kd + 1] = 14; /*H */

  /*C7H15O */
  ncf[68 * kd + 3] = 7;  /*C */
  ncf[68 * kd + 1] = 15; /*H */
  ncf[68 * kd + 2] = 1;  /*O */

  /*C3H7CHO */
  ncf[69 * kd + 3] = 4; /*C */
  ncf[69 * kd + 1] = 8; /*H */
  ncf[69 * kd + 2] = 1; /*O */

  /*C4H7 */
  ncf[70 * kd + 3] = 4; /*C */
  ncf[70 * kd + 1] = 7; /*H */

  /*C7H13 */
  ncf[71 * kd + 3] = 7;  /*C */
  ncf[71 * kd + 1] = 13; /*H */

  /*C5H9 */
  ncf[72 * kd + 3] = 5; /*C */
  ncf[72 * kd + 1] = 9; /*H */

  /*C4H7O */
  ncf[73 * kd + 3] = 4; /*C */
  ncf[73 * kd + 1] = 7; /*H */
  ncf[73 * kd + 2] = 1; /*O */

  /*N-C3H7O */
  ncf[74 * kd + 3] = 3; /*C */
  ncf[74 * kd + 1] = 7; /*H */
  ncf[74 * kd + 2] = 1; /*O */

  /*I-C8H18 */
  ncf[75 * kd + 3] = 8;  /*C */
  ncf[75 * kd + 1] = 18; /*H */

  /*Y-C7H15 */
  ncf[76 * kd + 3] = 7;  /*C */
  ncf[76 * kd + 1] = 15; /*H */

  /*I-C4H8 */
  ncf[77 * kd + 3] = 4; /*C */
  ncf[77 * kd + 1] = 8; /*H */

  /*I-C3H7 */
  ncf[78 * kd + 3] = 3; /*C */
  ncf[78 * kd + 1] = 7; /*H */

  /*T-C4H9 */
  ncf[79 * kd + 3] = 4; /*C */
  ncf[79 * kd + 1] = 9; /*H */

  /*C-C8H17 */
  ncf[80 * kd + 3] = 8;  /*C */
  ncf[80 * kd + 1] = 17; /*H */

  /*Y-C7H14 */
  ncf[81 * kd + 3] = 7;  /*C */
  ncf[81 * kd + 1] = 14; /*H */

  /*D-C8H17O */
  ncf[82 * kd + 3] = 8;  /*C */
  ncf[82 * kd + 1] = 17; /*H */
  ncf[82 * kd + 2] = 1;  /*O */

  /*CH3COCH3 */
  ncf[83 * kd + 3] = 3; /*C */
  ncf[83 * kd + 1] = 6; /*H */
  ncf[83 * kd + 2] = 1; /*O */

  /*I-C4H7 */
  ncf[84 * kd + 3] = 4; /*C */
  ncf[84 * kd + 1] = 7; /*H */

  /*X-C7H13 */
  ncf[85 * kd + 3] = 7;  /*C */
  ncf[85 * kd + 1] = 13; /*H */

  /*I-C3H5CHO */
  ncf[86 * kd + 3] = 4; /*C */
  ncf[86 * kd + 1] = 6; /*H */
  ncf[86 * kd + 2] = 1; /*O */

  /*T-C4H9O */
  ncf[87 * kd + 3] = 4; /*C */
  ncf[87 * kd + 1] = 9; /*H */
  ncf[87 * kd + 2] = 1; /*O */

  /*I-C4H7O */
  ncf[88 * kd + 3] = 4; /*C */
  ncf[88 * kd + 1] = 7; /*H */
  ncf[88 * kd + 2] = 1; /*O */

  /*C5H4CH2 */
  ncf[89 * kd + 1] = 6; /*H */
  ncf[89 * kd + 3] = 6; /*C */

  /*A1- */
  ncf[90 * kd + 1] = 5; /*H */
  ncf[90 * kd + 3] = 6; /*C */

  /*A1C2H2 */
  ncf[91 * kd + 1] = 7; /*H */
  ncf[91 * kd + 3] = 8; /*C */

  /*A1C2H3 */
  ncf[92 * kd + 1] = 8; /*H */
  ncf[92 * kd + 3] = 8; /*C */

  /*A1C2H */
  ncf[93 * kd + 1] = 6; /*H */
  ncf[93 * kd + 3] = 8; /*C */

  /*A1C2H* */
  ncf[94 * kd + 1] = 5; /*H */
  ncf[94 * kd + 3] = 8; /*C */

  /*A1C2H3* */
  ncf[95 * kd + 1] = 7; /*H */
  ncf[95 * kd + 3] = 8; /*C */

  /*A2- */
  ncf[96 * kd + 1] = 7;  /*H */
  ncf[96 * kd + 3] = 10; /*C */

  /*A2 */
  ncf[97 * kd + 1] = 8;  /*H */
  ncf[97 * kd + 3] = 10; /*C */

  /*A2* */
  ncf[98 * kd + 1] = 7;  /*H */
  ncf[98 * kd + 3] = 10; /*C */

  /*A2C2H2A */
  ncf[99 * kd + 1] = 9;  /*H */
  ncf[99 * kd + 3] = 12; /*C */

  /*A2C2H2B */
  ncf[100 * kd + 1] = 9;  /*H */
  ncf[100 * kd + 3] = 12; /*C */

  /*A2C2HA */
  ncf[101 * kd + 1] = 8;  /*H */
  ncf[101 * kd + 3] = 12; /*C */

  /*A2C2HB */
  ncf[102 * kd + 1] = 8;  /*H */
  ncf[102 * kd + 3] = 12; /*C */

  /*A2C2HA* */
  ncf[103 * kd + 1] = 7;  /*H */
  ncf[103 * kd + 3] = 12; /*C */

  /*A2C2HB* */
  ncf[104 * kd + 1] = 7;  /*H */
  ncf[104 * kd + 3] = 12; /*C */

  /*A2R5 */
  ncf[105 * kd + 1] = 8;  /*H */
  ncf[105 * kd + 3] = 12; /*C */

  /*A2R5- */
  ncf[106 * kd + 1] = 7;  /*H */
  ncf[106 * kd + 3] = 12; /*C */

  /*A2R5C2H2 */
  ncf[107 * kd + 1] = 9;  /*H */
  ncf[107 * kd + 3] = 14; /*C */

  /*A2R5C2H */
  ncf[108 * kd + 1] = 8;  /*H */
  ncf[108 * kd + 3] = 14; /*C */

  /*A2R5C2H* */
  ncf[109 * kd + 1] = 7;  /*H */
  ncf[109 * kd + 3] = 14; /*C */

  /*P2 */
  ncf[110 * kd + 1] = 10; /*H */
  ncf[110 * kd + 3] = 12; /*C */

  /*P2- */
  ncf[111 * kd + 1] = 9;  /*H */
  ncf[111 * kd + 3] = 12; /*C */

  /*A3- */
  ncf[112 * kd + 1] = 9;  /*H */
  ncf[112 * kd + 3] = 14; /*C */

  /*A3 */
  ncf[113 * kd + 1] = 10; /*H */
  ncf[113 * kd + 3] = 14; /*C */

  /*A3* */
  ncf[114 * kd + 1] = 9;  /*H */
  ncf[114 * kd + 3] = 14; /*C */

  /*A3R5- */
  ncf[115 * kd + 1] = 9;  /*H */
  ncf[115 * kd + 3] = 16; /*C */

  /*A3R5 */
  ncf[116 * kd + 1] = 10; /*H */
  ncf[116 * kd + 3] = 16; /*C */

  /*A4 */
  ncf[117 * kd + 1] = 10; /*H */
  ncf[117 * kd + 3] = 16; /*C */

  /*A4- */
  ncf[118 * kd + 1] = 9;  /*H */
  ncf[118 * kd + 3] = 16; /*C */

  /*A4R5 */
  ncf[119 * kd + 1] = 10; /*H */
  ncf[119 * kd + 3] = 18; /*C */

  /*FLTN */
  ncf[120 * kd + 1] = 10; /*H */
  ncf[120 * kd + 3] = 16; /*C */

  /*C5H6 */
  ncf[121 * kd + 1] = 6; /*H */
  ncf[121 * kd + 3] = 5; /*C */

  /*C5H5 */
  ncf[122 * kd + 1] = 5; /*H */
  ncf[122 * kd + 3] = 5; /*C */

  /*T-C5H5O */
  ncf[123 * kd + 3] = 5; /*C */
  ncf[123 * kd + 1] = 5; /*H */
  ncf[123 * kd + 2] = 1; /*O */

  /*C5H4O */
  ncf[124 * kd + 1] = 4; /*H */
  ncf[124 * kd + 3] = 5; /*C */
  ncf[124 * kd + 2] = 1; /*O */

  /*S-C5H5O */
  ncf[125 * kd + 3] = 5; /*C */
  ncf[125 * kd + 1] = 5; /*H */
  ncf[125 * kd + 2] = 1; /*O */

  /*C9H8 */
  ncf[126 * kd + 1] = 8; /*H */
  ncf[126 * kd + 3] = 9; /*C */

  /*C9H7 */
  ncf[127 * kd + 1] = 7; /*H */
  ncf[127 * kd + 3] = 9; /*C */

  /*A1CH2 */
  ncf[128 * kd + 1] = 7; /*H */
  ncf[128 * kd + 3] = 7; /*C */

  /*C9H6O */
  ncf[129 * kd + 1] = 6; /*H */
  ncf[129 * kd + 3] = 9; /*C */
  ncf[129 * kd + 2] = 1; /*O */

  /*O-C6H4 */
  ncf[130 * kd + 1] = 4; /*H */
  ncf[130 * kd + 3] = 6; /*C */

  /*A1CH3 */
  ncf[131 * kd + 1] = 8; /*H */
  ncf[131 * kd + 3] = 7; /*C */

  /*A1OH */
  ncf[132 * kd + 1] = 6; /*H */
  ncf[132 * kd + 3] = 6; /*C */
  ncf[132 * kd + 2] = 1; /*O */

  /*HOA1CH3 */
  ncf[133 * kd + 1] = 8; /*H */
  ncf[133 * kd + 3] = 7; /*C */
  ncf[133 * kd + 2] = 1; /*O */

  /*OA1CH3 */
  ncf[134 * kd + 1] = 7; /*H */
  ncf[134 * kd + 3] = 7; /*C */
  ncf[134 * kd + 2] = 1; /*O */

  /*A1CH2O */
  ncf[135 * kd + 1] = 7; /*H */
  ncf[135 * kd + 3] = 7; /*C */
  ncf[135 * kd + 2] = 1; /*O */

  /*A1CH2OH */
  ncf[136 * kd + 3] = 7; /*C */
  ncf[136 * kd + 1] = 8; /*H */
  ncf[136 * kd + 2] = 1; /*O */

  /*A1CHO */
  ncf[137 * kd + 1] = 6; /*H */
  ncf[137 * kd + 3] = 7; /*C */
  ncf[137 * kd + 2] = 1; /*O */

  /*A1O */
  ncf[138 * kd + 1] = 5; /*H */
  ncf[138 * kd + 3] = 6; /*C */
  ncf[138 * kd + 2] = 1; /*O */

  /*A1CH3* */
  ncf[139 * kd + 1] = 7; /*H */
  ncf[139 * kd + 3] = 7; /*C */

  /*A1C2H4 */
  ncf[140 * kd + 3] = 8; /*C */
  ncf[140 * kd + 1] = 9; /*H */

  /*A1C2H5 */
  ncf[141 * kd + 3] = 8;  /*C */
  ncf[141 * kd + 1] = 10; /*H */

  /*C8H9O2 */
  ncf[142 * kd + 1] = 9; /*H */
  ncf[142 * kd + 3] = 8; /*C */
  ncf[142 * kd + 2] = 2; /*O */

  /*C8H8OOH */
  ncf[143 * kd + 1] = 9; /*H */
  ncf[143 * kd + 3] = 8; /*C */
  ncf[143 * kd + 2] = 2; /*O */

  /*OC8H7OOH */
  ncf[144 * kd + 1] = 8; /*H */
  ncf[144 * kd + 3] = 8; /*C */
  ncf[144 * kd + 2] = 3; /*O */

  /*A1CH3CH3 */
  ncf[145 * kd + 1] = 10; /*H */
  ncf[145 * kd + 3] = 8;  /*C */

  /*A1CH3CH2 */
  ncf[146 * kd + 1] = 9; /*H */
  ncf[146 * kd + 3] = 8; /*C */

  /*A1CH3CHO */
  ncf[147 * kd + 1] = 8; /*H */
  ncf[147 * kd + 3] = 8; /*C */
  ncf[147 * kd + 2] = 1; /*O */

  /*A2CH3 */
  ncf[148 * kd + 1] = 10; /*H */
  ncf[148 * kd + 3] = 11; /*C */

  /*A1CHOCH2 */
  ncf[149 * kd + 1] = 7; /*H */
  ncf[149 * kd + 3] = 8; /*C */
  ncf[149 * kd + 2] = 1; /*O */

  /*A1CHOCHO */
  ncf[150 * kd + 1] = 6; /*H */
  ncf[150 * kd + 3] = 8; /*C */
  ncf[150 * kd + 2] = 2; /*O */

  /*A2OH */
  ncf[151 * kd + 3] = 10; /*C */
  ncf[151 * kd + 1] = 8;  /*H */
  ncf[151 * kd + 2] = 1;  /*O */

  /*A2CH2 */
  ncf[152 * kd + 1] = 9;  /*H */
  ncf[152 * kd + 3] = 11; /*C */

  /*A2CH2O */
  ncf[153 * kd + 1] = 9;  /*H */
  ncf[153 * kd + 3] = 11; /*C */
  ncf[153 * kd + 2] = 1;  /*O */

  /*A2CHO */
  ncf[154 * kd + 1] = 8;  /*H */
  ncf[154 * kd + 3] = 11; /*C */
  ncf[154 * kd + 2] = 1;  /*O */

  /*A2O */
  ncf[155 * kd + 3] = 10; /*C */
  ncf[155 * kd + 1] = 7;  /*H */
  ncf[155 * kd + 2] = 1;  /*O */

  /*OC6H4O */
  ncf[156 * kd + 1] = 4; /*H */
  ncf[156 * kd + 3] = 6; /*C */
  ncf[156 * kd + 2] = 2; /*O */

  /*N2 */
  ncf[157 * kd + 0] = 2; /*N */
}

/* Returns the vector of strings of element names */
void
CKSYME_STR(amrex::Vector<std::string>& ename)
{
  ename.resize(5);
  ename[0] = "N";
  ename[1] = "H";
  ename[2] = "O";
  ename[3] = "C";
  ename[4] = "AR";
}

/* Returns the vector of strings of species names */
void
CKSYMS_STR(amrex::Vector<std::string>& kname)
{
  kname.resize(158);
  kname[0] = "H";
  kname[1] = "O2";
  kname[2] = "O";
  kname[3] = "OH";
  kname[4] = "H2";
  kname[5] = "H2O";
  kname[6] = "CO2";
  kname[7] = "HO2";
  kname[8] = "H2O2";
  kname[9] = "CO";
  kname[10] = "HCO";
  kname[11] = "C";
  kname[12] = "CH";
  kname[13] = "T-CH2";
  kname[14] = "CH3";
  kname[15] = "CH2O";
  kname[16] = "HCCO";
  kname[17] = "C2H";
  kname[18] = "CH2CO";
  kname[19] = "C2H2";
  kname[20] = "S-CH2";
  kname[21] = "AR";
  kname[22] = "CH3OH";
  kname[23] = "CH2OH";
  kname[24] = "CH3O";
  kname[25] = "CH4";
  kname[26] = "CH3O2";
  kname[27] = "C2H3";
  kname[28] = "C2H4";
  kname[29] = "C2H5";
  kname[30] = "HCCOH";
  kname[31] = "CH2CHO";
  kname[32] = "CH3CHO";
  kname[33] = "H2C2";
  kname[34] = "C2H5O";
  kname[35] = "N-C3H7";
  kname[36] = "C2H6";
  kname[37] = "C3H8";
  kname[38] = "C3H6";
  kname[39] = "C3H3";
  kname[40] = "P-C3H4";
  kname[41] = "A-C3H4";
  kname[42] = "S-C3H5";
  kname[43] = "N-C4H3";
  kname[44] = "C2H3CHO";
  kname[45] = "A-C3H5";
  kname[46] = "C2O";
  kname[47] = "C4H4";
  kname[48] = "C3H2";
  kname[49] = "C3H2O";
  kname[50] = "C4H2";
  kname[51] = "I-C4H3";
  kname[52] = "T-C3H5";
  kname[53] = "C3H5O";
  kname[54] = "C4H";
  kname[55] = "C8H2";
  kname[56] = "C6H2";
  kname[57] = "C4H6";
  kname[58] = "N-C4H5";
  kname[59] = "I-C4H5";
  kname[60] = "A1";
  kname[61] = "N-C7H16";
  kname[62] = "C5H11";
  kname[63] = "P-C4H9";
  kname[64] = "C7H15";
  kname[65] = "P-C4H8";
  kname[66] = "C5H10";
  kname[67] = "C7H14";
  kname[68] = "C7H15O";
  kname[69] = "C3H7CHO";
  kname[70] = "C4H7";
  kname[71] = "C7H13";
  kname[72] = "C5H9";
  kname[73] = "C4H7O";
  kname[74] = "N-C3H7O";
  kname[75] = "I-C8H18";
  kname[76] = "Y-C7H15";
  kname[77] = "I-C4H8";
  kname[78] = "I-C3H7";
  kname[79] = "T-C4H9";
  kname[80] = "C-C8H17";
  kname[81] = "Y-C7H14";
  kname[82] = "D-C8H17O";
  kname[83] = "CH3COCH3";
  kname[84] = "I-C4H7";
  kname[85] = "X-C7H13";
  kname[86] = "I-C3H5CHO";
  kname[87] = "T-C4H9O";
  kname[88] = "I-C4H7O";
  kname[89] = "C5H4CH2";
  kname[90] = "A1-";
  kname[91] = "A1C2H2";
  kname[92] = "A1C2H3";
  kname[93] = "A1C2H";
  kname[94] = "A1C2H*";
  kname[95] = "A1C2H3*";
  kname[96] = "A2-";
  kname[97] = "A2";
  kname[98] = "A2*";
  kname[99] = "A2C2H2A";
  kname[100] = "A2C2H2B";
  kname[101] = "A2C2HA";
  kname[102] = "A2C2HB";
  kname[103] = "A2C2HA*";
  kname[104] = "A2C2HB*";
  kname[105] = "A2R5";
  kname[106] = "A2R5-";
  kname[107] = "A2R5C2H2";
  kname[108] = "A2R5C2H";
  kname[109] = "A2R5C2H*";
  kname[110] = "P2";
  kname[111] = "P2-";
  kname[112] = "A3-";
  kname[113] = "A3";
  kname[114] = "A3*";
  kname[115] = "A3R5-";
  kname[116] = "A3R5";
  kname[117] = "A4";
  kname[118] = "A4-";
  kname[119] = "A4R5";
  kname[120] = "FLTN";
  kname[121] = "C5H6";
  kname[122] = "C5H5";
  kname[123] = "T-C5H5O";
  kname[124] = "C5H4O";
  kname[125] = "S-C5H5O";
  kname[126] = "C9H8";
  kname[127] = "C9H7";
  kname[128] = "A1CH2";
  kname[129] = "C9H6O";
  kname[130] = "O-C6H4";
  kname[131] = "A1CH3";
  kname[132] = "A1OH";
  kname[133] = "HOA1CH3";
  kname[134] = "OA1CH3";
  kname[135] = "A1CH2O";
  kname[136] = "A1CH2OH";
  kname[137] = "A1CHO";
  kname[138] = "A1O";
  kname[139] = "A1CH3*";
  kname[140] = "A1C2H4";
  kname[141] = "A1C2H5";
  kname[142] = "C8H9O2";
  kname[143] = "C8H8OOH";
  kname[144] = "OC8H7OOH";
  kname[145] = "A1CH3CH3";
  kname[146] = "A1CH3CH2";
  kname[147] = "A1CH3CHO";
  kname[148] = "A2CH3";
  kname[149] = "A1CHOCH2";
  kname[150] = "A1CHOCHO";
  kname[151] = "A2OH";
  kname[152] = "A2CH2";
  kname[153] = "A2CH2O";
  kname[154] = "A2CHO";
  kname[155] = "A2O";
  kname[156] = "OC6H4O";
  kname[157] = "N2";
}

/*compute the sparsity pattern of the chemistry Jacobian */
void
SPARSITY_INFO(int* nJdata, int* consP, int NCELLS)
{
  amrex::Gpu::DeviceVector<amrex::Real> J_v(25281);
  amrex::Gpu::DeviceVector<amrex::Real> c_v(158);
  amrex::Real* J_d = J_v.data();
  amrex::Real* c_d = c_v.data();

  amrex::Real J_h[25281];

  amrex::IntVect iv(AMREX_D_DECL(0, 0, 0));
  amrex::ParallelFor(
    amrex::Box(iv, iv),
    [=] AMREX_GPU_HOST_DEVICE(int /*i*/, int /*j*/, int /*k*/) noexcept {
      for (int l = 0; l < 158; l++) {
        c_d[l] = 1.0 / 158.000000;
      }
      aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
  amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
  std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

  int nJdata_tmp = 0;
  for (int k = 0; k < 159; k++) {
    for (int l = 0; l < 159; l++) {
      if (J_h[159 * k + l] != 0.0) {
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
  amrex::Gpu::DeviceVector<amrex::Real> J_v(25281);
  amrex::Gpu::DeviceVector<amrex::Real> c_v(158);
  amrex::Real* J_d = J_v.data();
  amrex::Real* c_d = c_v.data();

  amrex::Real J_h[25281];

  amrex::IntVect iv(AMREX_D_DECL(0, 0, 0));
  amrex::ParallelFor(
    amrex::Box(iv, iv),
    [=] AMREX_GPU_HOST_DEVICE(int /*i*/, int /*j*/, int /*k*/) noexcept {
      for (int k = 0; k < 158; k++) {
        c_d[k] = 1.0 / 158.000000;
      }
      aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
  amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
  std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

  int nJdata_tmp = 0;
  for (int k = 0; k < 159; k++) {
    for (int l = 0; l < 159; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (J_h[159 * k + l] != 0.0) {
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
  amrex::Gpu::DeviceVector<amrex::Real> J_v(25281);
  amrex::Gpu::DeviceVector<amrex::Real> c_v(158);
  amrex::Real* J_d = J_v.data();
  amrex::Real* c_d = c_v.data();

  amrex::Real J_h[25281];

  amrex::IntVect iv(AMREX_D_DECL(0, 0, 0));
  amrex::ParallelFor(
    amrex::Box(iv, iv),
    [=] AMREX_GPU_HOST_DEVICE(int /*i*/, int /*j*/, int /*k*/) noexcept {
      for (int k = 0; k < 158; k++) {
        c_d[k] = 1.0 / 158.000000;
      }
      aJacobian_precond(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
  amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
  std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

  int nJdata_tmp = 0;
  for (int k = 0; k < 159; k++) {
    for (int l = 0; l < 159; l++) {
      if (k == l) {
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (J_h[159 * k + l] != 0.0) {
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

  amrex::Gpu::DeviceVector<amrex::Real> J_v(25281);
  amrex::Gpu::DeviceVector<amrex::Real> c_v(158);
  amrex::Real* J_d = J_v.data();
  amrex::Real* c_d = c_v.data();

  amrex::Real J_h[25281];

  amrex::IntVect iv(AMREX_D_DECL(0, 0, 0));
  amrex::ParallelFor(
    amrex::Box(iv, iv),
    [=] AMREX_GPU_HOST_DEVICE(int /*i*/, int /*j*/, int /*k*/) noexcept {
      for (int k = 0; k < 158; k++) {
        c_d[k] = 1.0 / 158.000000;
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
    offset_row = nc * 159;
    offset_col = nc * 159;
    for (int k = 0; k < 159; k++) {
      for (int l = 0; l < 159; l++) {
        if (J_h[159 * k + l] != 0.0) {
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
  amrex::Gpu::DeviceVector<amrex::Real> J_v(25281);
  amrex::Gpu::DeviceVector<amrex::Real> c_v(158);
  amrex::Real* J_d = J_v.data();
  amrex::Real* c_d = c_v.data();

  amrex::Real J_h[25281];

  amrex::IntVect iv(AMREX_D_DECL(0, 0, 0));
  amrex::ParallelFor(
    amrex::Box(iv, iv),
    [=] AMREX_GPU_HOST_DEVICE(int /*i*/, int /*j*/, int /*k*/) noexcept {
      for (int k = 0; k < 158; k++) {
        c_d[k] = 1.0 / 158.000000;
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
      offset = nc * 159;
      for (int l = 0; l < 159; l++) {
        for (int k = 0; k < 159; k++) {
          if (J_h[159 * k + l] != 0.0) {
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
      offset = nc * 159;
      for (int l = 0; l < 159; l++) {
        for (int k = 0; k < 159; k++) {
          if (J_h[159 * k + l] != 0.0) {
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
  amrex::Gpu::DeviceVector<amrex::Real> J_v(25281);
  amrex::Gpu::DeviceVector<amrex::Real> c_v(158);
  amrex::Real* J_d = J_v.data();
  amrex::Real* c_d = c_v.data();

  amrex::Real J_h[25281];

  amrex::IntVect iv(AMREX_D_DECL(0, 0, 0));
  amrex::ParallelFor(
    amrex::Box(iv, iv),
    [=] AMREX_GPU_HOST_DEVICE(int /*i*/, int /*j*/, int /*k*/) noexcept {
      for (int k = 0; k < 158; k++) {
        c_d[k] = 1.0 / 158.000000;
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
      offset = nc * 159;
      for (int l = 0; l < 159; l++) {
        for (int k = 0; k < 159; k++) {
          if (k == l) {
            colVals[nJdata_tmp - 1] = l + 1 + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (J_h[159 * k + l] != 0.0) {
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
      offset = nc * 159;
      for (int l = 0; l < 159; l++) {
        for (int k = 0; k < 159; k++) {
          if (k == l) {
            colVals[nJdata_tmp] = l + offset;
            nJdata_tmp = nJdata_tmp + 1;
          } else {
            if (J_h[159 * k + l] != 0.0) {
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
  amrex::Gpu::DeviceVector<amrex::Real> J_v(25281);
  amrex::Gpu::DeviceVector<amrex::Real> c_v(158);
  amrex::Real* J_d = J_v.data();
  amrex::Real* c_d = c_v.data();

  amrex::Real J_h[25281];

  amrex::IntVect iv(AMREX_D_DECL(0, 0, 0));
  amrex::ParallelFor(
    amrex::Box(iv, iv),
    [=] AMREX_GPU_HOST_DEVICE(int /*i*/, int /*j*/, int /*k*/) noexcept {
      for (int k = 0; k < 158; k++) {
        c_d[k] = 1.0 / 158.000000;
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
  for (int k = 0; k < 159; k++) {
    for (int l = 0; l < 159; l++) {
      if (k == l) {
        rowVals[nJdata_tmp] = l;
        indx[nJdata_tmp] = 159 * k + l;
        nJdata_tmp = nJdata_tmp + 1;
      } else {
        if (J_h[159 * k + l] != 0.0) {
          rowVals[nJdata_tmp] = l;
          indx[nJdata_tmp] = 159 * k + l;
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
  amrex::Gpu::DeviceVector<amrex::Real> J_v(25281);
  amrex::Gpu::DeviceVector<amrex::Real> c_v(158);
  amrex::Real* J_d = J_v.data();
  amrex::Real* c_d = c_v.data();

  amrex::Real J_h[25281];

  amrex::IntVect iv(AMREX_D_DECL(0, 0, 0));
  amrex::ParallelFor(
    amrex::Box(iv, iv),
    [=] AMREX_GPU_HOST_DEVICE(int /*i*/, int /*j*/, int /*k*/) noexcept {
      for (int k = 0; k < 158; k++) {
        c_d[k] = 1.0 / 158.000000;
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
    for (int l = 0; l < 159; l++) {
      for (int k = 0; k < 159; k++) {
        if (k == l) {
          colVals[nJdata_tmp - 1] = l + 1;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (J_h[159 * k + l] != 0.0) {
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
    for (int l = 0; l < 159; l++) {
      for (int k = 0; k < 159; k++) {
        if (k == l) {
          colVals[nJdata_tmp] = l;
          nJdata_tmp = nJdata_tmp + 1;
        } else {
          if (J_h[159 * k + l] != 0.0) {
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
