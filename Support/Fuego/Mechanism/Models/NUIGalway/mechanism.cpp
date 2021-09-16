#ifndef MECHANISM_CPP
#define MECHANISM_CPP

#include "mechanism.H"



/*save atomic weights into array */
void atomicWeight(amrex::Real *  awt)
{
    awt[0] = 12.011150; /*C */
    awt[1] = 1.007970; /*H */
    awt[2] = 14.006700; /*N */
    awt[3] = 15.999400; /*O */
    awt[4] = 39.948000; /*AR */
    awt[5] = 4.002600; /*HE */
}



/*get atomic weight for all elements */
void CKAWT( amrex::Real *  awt)
{
    atomicWeight(awt);
}



/*Returns the elemental composition  */
/*of the speciesi (mdim is num of elements) */
void CKNCF(int * ncf)
{
    int id; /*loop counter */
    int kd = 6; 
    /*Zero ncf */
    for (id = 0; id < kd * 118; ++ id) {
         ncf[id] = 0; 
    }

    /*H */
    ncf[ 0 * kd + 1 ] = 1; /*H */

    /*H2 */
    ncf[ 1 * kd + 1 ] = 2; /*H */

    /*C */
    ncf[ 2 * kd + 0 ] = 1; /*C */

    /*CH3COCH2O2H */
    ncf[ 3 * kd + 0 ] = 3; /*C */
    ncf[ 3 * kd + 1 ] = 6; /*H */
    ncf[ 3 * kd + 3 ] = 3; /*O */

    /*C3H2 */
    ncf[ 4 * kd + 1 ] = 2; /*H */
    ncf[ 4 * kd + 0 ] = 3; /*C */

    /*O */
    ncf[ 5 * kd + 3 ] = 1; /*O */

    /*CH */
    ncf[ 6 * kd + 0 ] = 1; /*C */
    ncf[ 6 * kd + 1 ] = 1; /*H */

    /*CH3COCH2O */
    ncf[ 7 * kd + 0 ] = 3; /*C */
    ncf[ 7 * kd + 1 ] = 5; /*H */
    ncf[ 7 * kd + 3 ] = 2; /*O */

    /*C3H5O */
    ncf[ 8 * kd + 0 ] = 3; /*C */
    ncf[ 8 * kd + 1 ] = 5; /*H */
    ncf[ 8 * kd + 3 ] = 1; /*O */

    /*O2 */
    ncf[ 9 * kd + 3 ] = 2; /*O */

    /*C2H6 */
    ncf[ 10 * kd + 0 ] = 2; /*C */
    ncf[ 10 * kd + 1 ] = 6; /*H */

    /*C2H3CHO */
    ncf[ 11 * kd + 0 ] = 3; /*C */
    ncf[ 11 * kd + 1 ] = 4; /*H */
    ncf[ 11 * kd + 3 ] = 1; /*O */

    /*C3H6OOH1-2 */
    ncf[ 12 * kd + 0 ] = 3; /*C */
    ncf[ 12 * kd + 1 ] = 7; /*H */
    ncf[ 12 * kd + 3 ] = 2; /*O */

    /*OH */
    ncf[ 13 * kd + 1 ] = 1; /*H */
    ncf[ 13 * kd + 3 ] = 1; /*O */

    /*C2H5 */
    ncf[ 14 * kd + 0 ] = 2; /*C */
    ncf[ 14 * kd + 1 ] = 5; /*H */

    /*C2H3CO */
    ncf[ 15 * kd + 0 ] = 3; /*C */
    ncf[ 15 * kd + 1 ] = 3; /*H */
    ncf[ 15 * kd + 3 ] = 1; /*O */

    /*C3H6OOH1-3 */
    ncf[ 16 * kd + 0 ] = 3; /*C */
    ncf[ 16 * kd + 1 ] = 7; /*H */
    ncf[ 16 * kd + 3 ] = 2; /*O */

    /*H2O */
    ncf[ 17 * kd + 1 ] = 2; /*H */
    ncf[ 17 * kd + 3 ] = 1; /*O */

    /*C2H4 */
    ncf[ 18 * kd + 0 ] = 2; /*C */
    ncf[ 18 * kd + 1 ] = 4; /*H */

    /*C2H5CHO */
    ncf[ 19 * kd + 0 ] = 3; /*C */
    ncf[ 19 * kd + 1 ] = 6; /*H */
    ncf[ 19 * kd + 3 ] = 1; /*O */

    /*C3H6OOH2-1 */
    ncf[ 20 * kd + 0 ] = 3; /*C */
    ncf[ 20 * kd + 1 ] = 7; /*H */
    ncf[ 20 * kd + 3 ] = 2; /*O */

    /*C2H3 */
    ncf[ 21 * kd + 0 ] = 2; /*C */
    ncf[ 21 * kd + 1 ] = 3; /*H */

    /*C2H5CO */
    ncf[ 22 * kd + 0 ] = 3; /*C */
    ncf[ 22 * kd + 1 ] = 5; /*H */
    ncf[ 22 * kd + 3 ] = 1; /*O */

    /*C3H6OOH1-2O2 */
    ncf[ 23 * kd + 0 ] = 3; /*C */
    ncf[ 23 * kd + 1 ] = 7; /*H */
    ncf[ 23 * kd + 3 ] = 4; /*O */

    /*HO2 */
    ncf[ 24 * kd + 1 ] = 1; /*H */
    ncf[ 24 * kd + 3 ] = 2; /*O */

    /*C2H2 */
    ncf[ 25 * kd + 0 ] = 2; /*C */
    ncf[ 25 * kd + 1 ] = 2; /*H */

    /*CH3OCH3 */
    ncf[ 26 * kd + 0 ] = 2; /*C */
    ncf[ 26 * kd + 1 ] = 6; /*H */
    ncf[ 26 * kd + 3 ] = 1; /*O */

    /*C3H6OOH1-3O2 */
    ncf[ 27 * kd + 0 ] = 3; /*C */
    ncf[ 27 * kd + 1 ] = 7; /*H */
    ncf[ 27 * kd + 3 ] = 4; /*O */

    /*H2O2 */
    ncf[ 28 * kd + 1 ] = 2; /*H */
    ncf[ 28 * kd + 3 ] = 2; /*O */

    /*C2H */
    ncf[ 29 * kd + 0 ] = 2; /*C */
    ncf[ 29 * kd + 1 ] = 1; /*H */

    /*CH3OCH2 */
    ncf[ 30 * kd + 0 ] = 2; /*C */
    ncf[ 30 * kd + 1 ] = 5; /*H */
    ncf[ 30 * kd + 3 ] = 1; /*O */

    /*C3H6OOH2-1O2 */
    ncf[ 31 * kd + 0 ] = 3; /*C */
    ncf[ 31 * kd + 1 ] = 7; /*H */
    ncf[ 31 * kd + 3 ] = 4; /*O */

    /*CH3CHO */
    ncf[ 32 * kd + 0 ] = 2; /*C */
    ncf[ 32 * kd + 3 ] = 1; /*O */
    ncf[ 32 * kd + 1 ] = 4; /*H */

    /*CH3OCH2O2 */
    ncf[ 33 * kd + 0 ] = 2; /*C */
    ncf[ 33 * kd + 1 ] = 5; /*H */
    ncf[ 33 * kd + 3 ] = 3; /*O */

    /*C3H6OOH2-2 */
    ncf[ 34 * kd + 0 ] = 3; /*C */
    ncf[ 34 * kd + 1 ] = 7; /*H */
    ncf[ 34 * kd + 3 ] = 2; /*O */

    /*CO */
    ncf[ 35 * kd + 0 ] = 1; /*C */
    ncf[ 35 * kd + 3 ] = 1; /*O */

    /*CH3CO */
    ncf[ 36 * kd + 0 ] = 2; /*C */
    ncf[ 36 * kd + 1 ] = 3; /*H */
    ncf[ 36 * kd + 3 ] = 1; /*O */

    /*CH2OCH2O2H */
    ncf[ 37 * kd + 0 ] = 2; /*C */
    ncf[ 37 * kd + 1 ] = 5; /*H */
    ncf[ 37 * kd + 3 ] = 3; /*O */

    /*NC3H7O2H */
    ncf[ 38 * kd + 0 ] = 3; /*C */
    ncf[ 38 * kd + 1 ] = 8; /*H */
    ncf[ 38 * kd + 3 ] = 2; /*O */

    /*CO2 */
    ncf[ 39 * kd + 0 ] = 1; /*C */
    ncf[ 39 * kd + 3 ] = 2; /*O */

    /*CH2CHO */
    ncf[ 40 * kd + 3 ] = 1; /*O */
    ncf[ 40 * kd + 1 ] = 3; /*H */
    ncf[ 40 * kd + 0 ] = 2; /*C */

    /*CH3OCH2O2H */
    ncf[ 41 * kd + 0 ] = 2; /*C */
    ncf[ 41 * kd + 1 ] = 6; /*H */
    ncf[ 41 * kd + 3 ] = 3; /*O */

    /*IC3H7O2H */
    ncf[ 42 * kd + 0 ] = 3; /*C */
    ncf[ 42 * kd + 1 ] = 8; /*H */
    ncf[ 42 * kd + 3 ] = 2; /*O */

    /*CH2O */
    ncf[ 43 * kd + 0 ] = 1; /*C */
    ncf[ 43 * kd + 1 ] = 2; /*H */
    ncf[ 43 * kd + 3 ] = 1; /*O */

    /*CH2CO */
    ncf[ 44 * kd + 0 ] = 2; /*C */
    ncf[ 44 * kd + 1 ] = 2; /*H */
    ncf[ 44 * kd + 3 ] = 1; /*O */

    /*CH3OCH2O */
    ncf[ 45 * kd + 0 ] = 2; /*C */
    ncf[ 45 * kd + 1 ] = 5; /*H */
    ncf[ 45 * kd + 3 ] = 2; /*O */

    /*NC3H7O2 */
    ncf[ 46 * kd + 0 ] = 3; /*C */
    ncf[ 46 * kd + 1 ] = 7; /*H */
    ncf[ 46 * kd + 3 ] = 2; /*O */

    /*HCO */
    ncf[ 47 * kd + 1 ] = 1; /*H */
    ncf[ 47 * kd + 0 ] = 1; /*C */
    ncf[ 47 * kd + 3 ] = 1; /*O */

    /*HCCO */
    ncf[ 48 * kd + 1 ] = 1; /*H */
    ncf[ 48 * kd + 0 ] = 2; /*C */
    ncf[ 48 * kd + 3 ] = 1; /*O */

    /*O2CH2OCH2O2H */
    ncf[ 49 * kd + 0 ] = 2; /*C */
    ncf[ 49 * kd + 1 ] = 5; /*H */
    ncf[ 49 * kd + 3 ] = 5; /*O */

    /*IC3H7O2 */
    ncf[ 50 * kd + 0 ] = 3; /*C */
    ncf[ 50 * kd + 1 ] = 7; /*H */
    ncf[ 50 * kd + 3 ] = 2; /*O */

    /*HO2CHO */
    ncf[ 51 * kd + 0 ] = 1; /*C */
    ncf[ 51 * kd + 1 ] = 2; /*H */
    ncf[ 51 * kd + 3 ] = 3; /*O */

    /*HCCOH */
    ncf[ 52 * kd + 1 ] = 2; /*H */
    ncf[ 52 * kd + 0 ] = 2; /*C */
    ncf[ 52 * kd + 3 ] = 1; /*O */

    /*HO2CH2OCHO */
    ncf[ 53 * kd + 0 ] = 2; /*C */
    ncf[ 53 * kd + 1 ] = 4; /*H */
    ncf[ 53 * kd + 3 ] = 4; /*O */

    /*NC3H7O */
    ncf[ 54 * kd + 0 ] = 3; /*C */
    ncf[ 54 * kd + 1 ] = 7; /*H */
    ncf[ 54 * kd + 3 ] = 1; /*O */

    /*O2CHO */
    ncf[ 55 * kd + 0 ] = 1; /*C */
    ncf[ 55 * kd + 1 ] = 1; /*H */
    ncf[ 55 * kd + 3 ] = 3; /*O */

    /*CH3CO3H */
    ncf[ 56 * kd + 0 ] = 2; /*C */
    ncf[ 56 * kd + 1 ] = 4; /*H */
    ncf[ 56 * kd + 3 ] = 3; /*O */

    /*OCH2OCHO */
    ncf[ 57 * kd + 0 ] = 2; /*C */
    ncf[ 57 * kd + 1 ] = 3; /*H */
    ncf[ 57 * kd + 3 ] = 3; /*O */

    /*IC3H7O */
    ncf[ 58 * kd + 0 ] = 3; /*C */
    ncf[ 58 * kd + 1 ] = 7; /*H */
    ncf[ 58 * kd + 3 ] = 1; /*O */

    /*HOCHO */
    ncf[ 59 * kd + 0 ] = 1; /*C */
    ncf[ 59 * kd + 1 ] = 2; /*H */
    ncf[ 59 * kd + 3 ] = 2; /*O */

    /*CH3CO3 */
    ncf[ 60 * kd + 0 ] = 2; /*C */
    ncf[ 60 * kd + 1 ] = 3; /*H */
    ncf[ 60 * kd + 3 ] = 3; /*O */

    /*HOCH2OCO */
    ncf[ 61 * kd + 0 ] = 2; /*C */
    ncf[ 61 * kd + 1 ] = 3; /*H */
    ncf[ 61 * kd + 3 ] = 3; /*O */

    /*C3H6O1-2 */
    ncf[ 62 * kd + 0 ] = 3; /*C */
    ncf[ 62 * kd + 1 ] = 6; /*H */
    ncf[ 62 * kd + 3 ] = 1; /*O */

    /*OCHO */
    ncf[ 63 * kd + 0 ] = 1; /*C */
    ncf[ 63 * kd + 1 ] = 1; /*H */
    ncf[ 63 * kd + 3 ] = 2; /*O */

    /*CH3CO2 */
    ncf[ 64 * kd + 0 ] = 2; /*C */
    ncf[ 64 * kd + 1 ] = 3; /*H */
    ncf[ 64 * kd + 3 ] = 2; /*O */

    /*CH3OCHO */
    ncf[ 65 * kd + 0 ] = 2; /*C */
    ncf[ 65 * kd + 1 ] = 4; /*H */
    ncf[ 65 * kd + 3 ] = 2; /*O */

    /*C3H6O1-3 */
    ncf[ 66 * kd + 0 ] = 3; /*C */
    ncf[ 66 * kd + 1 ] = 6; /*H */
    ncf[ 66 * kd + 3 ] = 1; /*O */

    /*HOCH2O2H */
    ncf[ 67 * kd + 0 ] = 1; /*C */
    ncf[ 67 * kd + 1 ] = 4; /*H */
    ncf[ 67 * kd + 3 ] = 3; /*O */

    /*C2H5OH */
    ncf[ 68 * kd + 0 ] = 2; /*C */
    ncf[ 68 * kd + 1 ] = 6; /*H */
    ncf[ 68 * kd + 3 ] = 1; /*O */

    /*CH3OCO */
    ncf[ 69 * kd + 0 ] = 2; /*C */
    ncf[ 69 * kd + 1 ] = 3; /*H */
    ncf[ 69 * kd + 3 ] = 2; /*O */

    /*C3KET12 */
    ncf[ 70 * kd + 0 ] = 3; /*C */
    ncf[ 70 * kd + 1 ] = 6; /*H */
    ncf[ 70 * kd + 3 ] = 3; /*O */

    /*HOCH2O2 */
    ncf[ 71 * kd + 0 ] = 1; /*C */
    ncf[ 71 * kd + 1 ] = 3; /*H */
    ncf[ 71 * kd + 3 ] = 3; /*O */

    /*C2H5O */
    ncf[ 72 * kd + 0 ] = 2; /*C */
    ncf[ 72 * kd + 1 ] = 5; /*H */
    ncf[ 72 * kd + 3 ] = 1; /*O */

    /*CH2OCHO */
    ncf[ 73 * kd + 0 ] = 2; /*C */
    ncf[ 73 * kd + 1 ] = 3; /*H */
    ncf[ 73 * kd + 3 ] = 2; /*O */

    /*C3KET13 */
    ncf[ 74 * kd + 0 ] = 3; /*C */
    ncf[ 74 * kd + 1 ] = 6; /*H */
    ncf[ 74 * kd + 3 ] = 3; /*O */

    /*OCH2O2H */
    ncf[ 75 * kd + 0 ] = 1; /*C */
    ncf[ 75 * kd + 1 ] = 3; /*H */
    ncf[ 75 * kd + 3 ] = 3; /*O */

    /*PC2H4OH */
    ncf[ 76 * kd + 0 ] = 2; /*C */
    ncf[ 76 * kd + 1 ] = 5; /*H */
    ncf[ 76 * kd + 3 ] = 1; /*O */

    /*C3KET21 */
    ncf[ 77 * kd + 0 ] = 3; /*C */
    ncf[ 77 * kd + 1 ] = 6; /*H */
    ncf[ 77 * kd + 3 ] = 3; /*O */

    /*HOCH2O */
    ncf[ 78 * kd + 0 ] = 1; /*C */
    ncf[ 78 * kd + 1 ] = 3; /*H */
    ncf[ 78 * kd + 3 ] = 2; /*O */

    /*SC2H4OH */
    ncf[ 79 * kd + 0 ] = 2; /*C */
    ncf[ 79 * kd + 1 ] = 5; /*H */
    ncf[ 79 * kd + 3 ] = 1; /*O */

    /*C3H8 */
    ncf[ 80 * kd + 0 ] = 3; /*C */
    ncf[ 80 * kd + 1 ] = 8; /*H */

    /*C3H51-23OOH */
    ncf[ 81 * kd + 0 ] = 3; /*C */
    ncf[ 81 * kd + 1 ] = 7; /*H */
    ncf[ 81 * kd + 3 ] = 4; /*O */

    /*CH3OH */
    ncf[ 82 * kd + 0 ] = 1; /*C */
    ncf[ 82 * kd + 1 ] = 4; /*H */
    ncf[ 82 * kd + 3 ] = 1; /*O */

    /*O2C2H4OH */
    ncf[ 83 * kd + 0 ] = 2; /*C */
    ncf[ 83 * kd + 1 ] = 5; /*H */
    ncf[ 83 * kd + 3 ] = 3; /*O */

    /*IC3H7 */
    ncf[ 84 * kd + 0 ] = 3; /*C */
    ncf[ 84 * kd + 1 ] = 7; /*H */

    /*C3H52-13OOH */
    ncf[ 85 * kd + 0 ] = 3; /*C */
    ncf[ 85 * kd + 1 ] = 7; /*H */
    ncf[ 85 * kd + 3 ] = 4; /*O */

    /*CH2OH */
    ncf[ 86 * kd + 1 ] = 3; /*H */
    ncf[ 86 * kd + 0 ] = 1; /*C */
    ncf[ 86 * kd + 3 ] = 1; /*O */

    /*C2H5O2H */
    ncf[ 87 * kd + 0 ] = 2; /*C */
    ncf[ 87 * kd + 1 ] = 6; /*H */
    ncf[ 87 * kd + 3 ] = 2; /*O */

    /*NC3H7 */
    ncf[ 88 * kd + 0 ] = 3; /*C */
    ncf[ 88 * kd + 1 ] = 7; /*H */

    /*C3H6OH */
    ncf[ 89 * kd + 0 ] = 3; /*C */
    ncf[ 89 * kd + 1 ] = 7; /*H */
    ncf[ 89 * kd + 3 ] = 1; /*O */

    /*CH3O */
    ncf[ 90 * kd + 0 ] = 1; /*C */
    ncf[ 90 * kd + 1 ] = 3; /*H */
    ncf[ 90 * kd + 3 ] = 1; /*O */

    /*C2H5O2 */
    ncf[ 91 * kd + 0 ] = 2; /*C */
    ncf[ 91 * kd + 1 ] = 5; /*H */
    ncf[ 91 * kd + 3 ] = 2; /*O */

    /*C3H6 */
    ncf[ 92 * kd + 0 ] = 3; /*C */
    ncf[ 92 * kd + 1 ] = 6; /*H */

    /*HOC3H6O2 */
    ncf[ 93 * kd + 0 ] = 3; /*C */
    ncf[ 93 * kd + 1 ] = 7; /*H */
    ncf[ 93 * kd + 3 ] = 3; /*O */

    /*CH3O2H */
    ncf[ 94 * kd + 0 ] = 1; /*C */
    ncf[ 94 * kd + 1 ] = 4; /*H */
    ncf[ 94 * kd + 3 ] = 2; /*O */

    /*C2H4O2H */
    ncf[ 95 * kd + 0 ] = 2; /*C */
    ncf[ 95 * kd + 1 ] = 5; /*H */
    ncf[ 95 * kd + 3 ] = 2; /*O */

    /*C3H5-A */
    ncf[ 96 * kd + 0 ] = 3; /*C */
    ncf[ 96 * kd + 1 ] = 5; /*H */

    /*CH3CHCO */
    ncf[ 97 * kd + 0 ] = 3; /*C */
    ncf[ 97 * kd + 1 ] = 4; /*H */
    ncf[ 97 * kd + 3 ] = 1; /*O */

    /*CH3O2 */
    ncf[ 98 * kd + 0 ] = 1; /*C */
    ncf[ 98 * kd + 1 ] = 3; /*H */
    ncf[ 98 * kd + 3 ] = 2; /*O */

    /*C2H4O1-2 */
    ncf[ 99 * kd + 0 ] = 2; /*C */
    ncf[ 99 * kd + 1 ] = 4; /*H */
    ncf[ 99 * kd + 3 ] = 1; /*O */

    /*C3H5-S */
    ncf[ 100 * kd + 0 ] = 3; /*C */
    ncf[ 100 * kd + 1 ] = 5; /*H */

    /*AC3H5OOH */
    ncf[ 101 * kd + 0 ] = 3; /*C */
    ncf[ 101 * kd + 1 ] = 6; /*H */
    ncf[ 101 * kd + 3 ] = 2; /*O */

    /*CH4 */
    ncf[ 102 * kd + 0 ] = 1; /*C */
    ncf[ 102 * kd + 1 ] = 4; /*H */

    /*C2H3O1-2 */
    ncf[ 103 * kd + 0 ] = 2; /*C */
    ncf[ 103 * kd + 1 ] = 3; /*H */
    ncf[ 103 * kd + 3 ] = 1; /*O */

    /*C3H5-T */
    ncf[ 104 * kd + 0 ] = 3; /*C */
    ncf[ 104 * kd + 1 ] = 5; /*H */

    /*C2H3OOH */
    ncf[ 105 * kd + 0 ] = 2; /*C */
    ncf[ 105 * kd + 1 ] = 4; /*H */
    ncf[ 105 * kd + 3 ] = 2; /*O */

    /*CH3 */
    ncf[ 106 * kd + 0 ] = 1; /*C */
    ncf[ 106 * kd + 1 ] = 3; /*H */

    /*CH3COCH3 */
    ncf[ 107 * kd + 0 ] = 3; /*C */
    ncf[ 107 * kd + 1 ] = 6; /*H */
    ncf[ 107 * kd + 3 ] = 1; /*O */

    /*C3H4-P */
    ncf[ 108 * kd + 1 ] = 4; /*H */
    ncf[ 108 * kd + 0 ] = 3; /*C */

    /*CH2 */
    ncf[ 109 * kd + 0 ] = 1; /*C */
    ncf[ 109 * kd + 1 ] = 2; /*H */

    /*CH3COCH2 */
    ncf[ 110 * kd + 0 ] = 3; /*C */
    ncf[ 110 * kd + 1 ] = 5; /*H */
    ncf[ 110 * kd + 3 ] = 1; /*O */

    /*C3H4-A */
    ncf[ 111 * kd + 1 ] = 4; /*H */
    ncf[ 111 * kd + 0 ] = 3; /*C */

    /*CH2(S) */
    ncf[ 112 * kd + 0 ] = 1; /*C */
    ncf[ 112 * kd + 1 ] = 2; /*H */

    /*CH3COCH2O2 */
    ncf[ 113 * kd + 0 ] = 3; /*C */
    ncf[ 113 * kd + 1 ] = 5; /*H */
    ncf[ 113 * kd + 3 ] = 3; /*O */

    /*C3H3 */
    ncf[ 114 * kd + 0 ] = 3; /*C */
    ncf[ 114 * kd + 1 ] = 3; /*H */

    /*HE */
    ncf[ 115 * kd + 5 ] = 1; /*HE */

    /*AR */
    ncf[ 116 * kd + 4 ] = 1; /*AR */

    /*N2 */
    ncf[ 117 * kd + 2 ] = 2; /*N */

}


/* Returns the vector of strings of element names */
void CKSYME_STR(amrex::Vector<std::string>& ename)
{
    ename.resize(6);
    ename[0] = "C";
    ename[1] = "H";
    ename[2] = "N";
    ename[3] = "O";
    ename[4] = "AR";
    ename[5] = "HE";
}


/* Returns the vector of strings of species names */
void CKSYMS_STR(amrex::Vector<std::string>& kname)
{
    kname.resize(118);
    kname[0] = "H";
    kname[1] = "H2";
    kname[2] = "C";
    kname[3] = "CH3COCH2O2H";
    kname[4] = "C3H2";
    kname[5] = "O";
    kname[6] = "CH";
    kname[7] = "CH3COCH2O";
    kname[8] = "C3H5O";
    kname[9] = "O2";
    kname[10] = "C2H6";
    kname[11] = "C2H3CHO";
    kname[12] = "C3H6OOH1-2";
    kname[13] = "OH";
    kname[14] = "C2H5";
    kname[15] = "C2H3CO";
    kname[16] = "C3H6OOH1-3";
    kname[17] = "H2O";
    kname[18] = "C2H4";
    kname[19] = "C2H5CHO";
    kname[20] = "C3H6OOH2-1";
    kname[21] = "C2H3";
    kname[22] = "C2H5CO";
    kname[23] = "C3H6OOH1-2O2";
    kname[24] = "HO2";
    kname[25] = "C2H2";
    kname[26] = "CH3OCH3";
    kname[27] = "C3H6OOH1-3O2";
    kname[28] = "H2O2";
    kname[29] = "C2H";
    kname[30] = "CH3OCH2";
    kname[31] = "C3H6OOH2-1O2";
    kname[32] = "CH3CHO";
    kname[33] = "CH3OCH2O2";
    kname[34] = "C3H6OOH2-2";
    kname[35] = "CO";
    kname[36] = "CH3CO";
    kname[37] = "CH2OCH2O2H";
    kname[38] = "NC3H7O2H";
    kname[39] = "CO2";
    kname[40] = "CH2CHO";
    kname[41] = "CH3OCH2O2H";
    kname[42] = "IC3H7O2H";
    kname[43] = "CH2O";
    kname[44] = "CH2CO";
    kname[45] = "CH3OCH2O";
    kname[46] = "NC3H7O2";
    kname[47] = "HCO";
    kname[48] = "HCCO";
    kname[49] = "O2CH2OCH2O2H";
    kname[50] = "IC3H7O2";
    kname[51] = "HO2CHO";
    kname[52] = "HCCOH";
    kname[53] = "HO2CH2OCHO";
    kname[54] = "NC3H7O";
    kname[55] = "O2CHO";
    kname[56] = "CH3CO3H";
    kname[57] = "OCH2OCHO";
    kname[58] = "IC3H7O";
    kname[59] = "HOCHO";
    kname[60] = "CH3CO3";
    kname[61] = "HOCH2OCO";
    kname[62] = "C3H6O1-2";
    kname[63] = "OCHO";
    kname[64] = "CH3CO2";
    kname[65] = "CH3OCHO";
    kname[66] = "C3H6O1-3";
    kname[67] = "HOCH2O2H";
    kname[68] = "C2H5OH";
    kname[69] = "CH3OCO";
    kname[70] = "C3KET12";
    kname[71] = "HOCH2O2";
    kname[72] = "C2H5O";
    kname[73] = "CH2OCHO";
    kname[74] = "C3KET13";
    kname[75] = "OCH2O2H";
    kname[76] = "PC2H4OH";
    kname[77] = "C3KET21";
    kname[78] = "HOCH2O";
    kname[79] = "SC2H4OH";
    kname[80] = "C3H8";
    kname[81] = "C3H51-23OOH";
    kname[82] = "CH3OH";
    kname[83] = "O2C2H4OH";
    kname[84] = "IC3H7";
    kname[85] = "C3H52-13OOH";
    kname[86] = "CH2OH";
    kname[87] = "C2H5O2H";
    kname[88] = "NC3H7";
    kname[89] = "C3H6OH";
    kname[90] = "CH3O";
    kname[91] = "C2H5O2";
    kname[92] = "C3H6";
    kname[93] = "HOC3H6O2";
    kname[94] = "CH3O2H";
    kname[95] = "C2H4O2H";
    kname[96] = "C3H5-A";
    kname[97] = "CH3CHCO";
    kname[98] = "CH3O2";
    kname[99] = "C2H4O1-2";
    kname[100] = "C3H5-S";
    kname[101] = "AC3H5OOH";
    kname[102] = "CH4";
    kname[103] = "C2H3O1-2";
    kname[104] = "C3H5-T";
    kname[105] = "C2H3OOH";
    kname[106] = "CH3";
    kname[107] = "CH3COCH3";
    kname[108] = "C3H4-P";
    kname[109] = "CH2";
    kname[110] = "CH3COCH2";
    kname[111] = "C3H4-A";
    kname[112] = "CH2(S)";
    kname[113] = "CH3COCH2O2";
    kname[114] = "C3H3";
    kname[115] = "HE";
    kname[116] = "AR";
    kname[117] = "N2";
}

/*compute the sparsity pattern of the chemistry Jacobian */
void SPARSITY_INFO( int * nJdata, const int * consP, int NCELLS)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(14161);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(118);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[14161];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int l=0; l<118; l++) {
                c_d[l] = 1.0/ 118.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<119; k++) {
        for (int l=0; l<119; l++) {
            if(J_h[ 119 * k + l] != 0.0){
                nJdata_tmp = nJdata_tmp + 1;
            }
        }
    }

    *nJdata = NCELLS * nJdata_tmp;
}



/*compute the sparsity pattern of the system Jacobian */
void SPARSITY_INFO_SYST( int * nJdata, const int * consP, int NCELLS)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(14161);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(118);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[14161];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<118; k++) {
                c_d[k] = 1.0/ 118.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<119; k++) {
        for (int l=0; l<119; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J_h[ 119 * k + l] != 0.0){
                    nJdata_tmp = nJdata_tmp + 1;
                }
            }
        }
    }

    *nJdata = NCELLS * nJdata_tmp;
}



/*compute the sparsity pattern of the simplified (for preconditioning) system Jacobian */
void SPARSITY_INFO_SYST_SIMPLIFIED( int * nJdata, const int * consP)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(14161);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(118);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[14161];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<118; k++) {
                c_d[k] = 1.0/ 118.000000 ;
            }
            aJacobian_precond(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<119; k++) {
        for (int l=0; l<119; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J_h[ 119 * k + l] != 0.0){
                    nJdata_tmp = nJdata_tmp + 1;
                }
            }
        }
    }

    nJdata[0] = nJdata_tmp;
}


/*compute the sparsity pattern of the chemistry Jacobian in CSC format -- base 0 */
void SPARSITY_PREPROC_CSC(int *  rowVals, int *  colPtrs, const int * consP, int NCELLS)
{
    int offset_row;
    int offset_col;

    amrex::Gpu::DeviceVector<amrex::Real> J_v(14161);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(118);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[14161];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<118; k++) {
                c_d[k] = 1.0/ 118.000000 ;
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
    for (int nc=0; nc<NCELLS; nc++) {
        offset_row = nc * 119;
        offset_col = nc * 119;
        for (int k=0; k<119; k++) {
            for (int l=0; l<119; l++) {
                if(J_h[119*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l + offset_row; 
                    nJdata_tmp = nJdata_tmp + 1; 
                }
            }
            colPtrs[offset_col + (k + 1)] = nJdata_tmp;
        }
    }
}

/*compute the sparsity pattern of the chemistry Jacobian in CSR format -- base 0 */
void SPARSITY_PREPROC_CSR(int * colVals, int * rowPtrs, const int * consP, int NCELLS, int base)
{
    int offset;
    amrex::Gpu::DeviceVector<amrex::Real> J_v(14161);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(118);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[14161];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<118; k++) {
                c_d[k] = 1.0/ 118.000000 ;
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
        for (int nc=0; nc<NCELLS; nc++) {
            offset = nc * 119;
            for (int l=0; l<119; l++) {
                for (int k=0; k<119; k++) {
                    if(J_h[119*k + l] != 0.0) {
                        colVals[nJdata_tmp-1] = k+1 + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    }
                }
                rowPtrs[offset + (l + 1)] = nJdata_tmp;
            }
        }
    } else {
        rowPtrs[0] = 0;
        int nJdata_tmp = 0;
        for (int nc=0; nc<NCELLS; nc++) {
            offset = nc * 119;
            for (int l=0; l<119; l++) {
                for (int k=0; k<119; k++) {
                    if(J_h[119*k + l] != 0.0) {
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
void SPARSITY_PREPROC_SYST_CSR(int * colVals, int * rowPtr, const int * consP, int NCELLS, int base)
{
    int offset;
    amrex::Gpu::DeviceVector<amrex::Real> J_v(14161);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(118);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[14161];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<118; k++) {
                c_d[k] = 1.0/ 118.000000 ;
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
        for (int nc=0; nc<NCELLS; nc++) {
            offset = nc * 119;
            for (int l=0; l<119; l++) {
                for (int k=0; k<119; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp-1] = l+1 + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(J_h[119*k + l] != 0.0) {
                            colVals[nJdata_tmp-1] = k+1 + offset; 
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
        for (int nc=0; nc<NCELLS; nc++) {
            offset = nc * 119;
            for (int l=0; l<119; l++) {
                for (int k=0; k<119; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp] = l + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(J_h[119*k + l] != 0.0) {
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

/*compute the sparsity pattern of the simplified (for precond) system Jacobian on CPU */
/*BASE 0 */
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(int * rowVals, int * colPtrs, int * indx, const int * consP)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(14161);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(118);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[14161];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<118; k++) {
                c_d[k] = 1.0/ 118.000000 ;
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
    for (int k=0; k<119; k++) {
        for (int l=0; l<119; l++) {
            if (k == l) {
                rowVals[nJdata_tmp] = l; 
                indx[nJdata_tmp] = 119*k + l;
                nJdata_tmp = nJdata_tmp + 1; 
            } else {
                if(J_h[119*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l; 
                    indx[nJdata_tmp] = 119*k + l;
                    nJdata_tmp = nJdata_tmp + 1; 
                }
            }
        }
        colPtrs[k+1] = nJdata_tmp;
    }
}

/*compute the sparsity pattern of the simplified (for precond) system Jacobian */
/*CSR format BASE is under choice */
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(int * colVals, int * rowPtr, const int * consP, int base)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(14161);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(118);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[14161];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<118; k++) {
                c_d[k] = 1.0/ 118.000000 ;
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
        for (int l=0; l<119; l++) {
            for (int k=0; k<119; k++) {
                if (k == l) {
                    colVals[nJdata_tmp-1] = l+1; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J_h[119*k + l] != 0.0) {
                        colVals[nJdata_tmp-1] = k+1; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    }
                }
            }
            rowPtr[l+1] = nJdata_tmp;
        }
    } else {
        rowPtr[0] = 0;
        int nJdata_tmp = 0;
        for (int l=0; l<119; l++) {
            for (int k=0; k<119; k++) {
                if (k == l) {
                    colVals[nJdata_tmp] = l; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J_h[119*k + l] != 0.0) {
                        colVals[nJdata_tmp] = k; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    }
                }
            }
            rowPtr[l+1] = nJdata_tmp;
        }
    }
}

#endif
