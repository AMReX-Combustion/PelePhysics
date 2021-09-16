#ifndef MECHANISM_CPP
#define MECHANISM_CPP

#include "mechanism.H"



/*save atomic weights into array */
void atomicWeight(amrex::Real *  awt)
{
    awt[0] = 1.007970; /*H */
    awt[1] = 12.011150; /*C */
    awt[2] = 15.999400; /*O */
    awt[3] = 14.006700; /*N */
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
    int kd = 4; 
    /*Zero ncf */
    for (id = 0; id < kd * 88; ++ id) {
         ncf[id] = 0; 
    }

    /*H */
    ncf[ 0 * kd + 0 ] = 1; /*H */

    /*H2 */
    ncf[ 1 * kd + 0 ] = 2; /*H */

    /*O */
    ncf[ 2 * kd + 2 ] = 1; /*O */

    /*O2 */
    ncf[ 3 * kd + 2 ] = 2; /*O */

    /*OH */
    ncf[ 4 * kd + 0 ] = 1; /*H */
    ncf[ 4 * kd + 2 ] = 1; /*O */

    /*H2O */
    ncf[ 5 * kd + 0 ] = 2; /*H */
    ncf[ 5 * kd + 2 ] = 1; /*O */

    /*CO */
    ncf[ 6 * kd + 1 ] = 1; /*C */
    ncf[ 6 * kd + 2 ] = 1; /*O */

    /*HCO */
    ncf[ 7 * kd + 0 ] = 1; /*H */
    ncf[ 7 * kd + 1 ] = 1; /*C */
    ncf[ 7 * kd + 2 ] = 1; /*O */

    /*CO2 */
    ncf[ 8 * kd + 1 ] = 1; /*C */
    ncf[ 8 * kd + 2 ] = 2; /*O */

    /*CH3 */
    ncf[ 9 * kd + 1 ] = 1; /*C */
    ncf[ 9 * kd + 0 ] = 3; /*H */

    /*CH4 */
    ncf[ 10 * kd + 1 ] = 1; /*C */
    ncf[ 10 * kd + 0 ] = 4; /*H */

    /*HO2 */
    ncf[ 11 * kd + 0 ] = 1; /*H */
    ncf[ 11 * kd + 2 ] = 2; /*O */

    /*H2O2 */
    ncf[ 12 * kd + 0 ] = 2; /*H */
    ncf[ 12 * kd + 2 ] = 2; /*O */

    /*CH2O */
    ncf[ 13 * kd + 1 ] = 1; /*C */
    ncf[ 13 * kd + 0 ] = 2; /*H */
    ncf[ 13 * kd + 2 ] = 1; /*O */

    /*CH3O */
    ncf[ 14 * kd + 1 ] = 1; /*C */
    ncf[ 14 * kd + 0 ] = 3; /*H */
    ncf[ 14 * kd + 2 ] = 1; /*O */

    /*C2H6 */
    ncf[ 15 * kd + 1 ] = 2; /*C */
    ncf[ 15 * kd + 0 ] = 6; /*H */

    /*C2H4 */
    ncf[ 16 * kd + 1 ] = 2; /*C */
    ncf[ 16 * kd + 0 ] = 4; /*H */

    /*C2H5 */
    ncf[ 17 * kd + 1 ] = 2; /*C */
    ncf[ 17 * kd + 0 ] = 5; /*H */

    /*C2H2 */
    ncf[ 18 * kd + 1 ] = 2; /*C */
    ncf[ 18 * kd + 0 ] = 2; /*H */

    /*C2H3 */
    ncf[ 19 * kd + 1 ] = 2; /*C */
    ncf[ 19 * kd + 0 ] = 3; /*H */

    /*CH2CO */
    ncf[ 20 * kd + 1 ] = 2; /*C */
    ncf[ 20 * kd + 0 ] = 2; /*H */
    ncf[ 20 * kd + 2 ] = 1; /*O */

    /*HCCO */
    ncf[ 21 * kd + 0 ] = 1; /*H */
    ncf[ 21 * kd + 1 ] = 2; /*C */
    ncf[ 21 * kd + 2 ] = 1; /*O */

    /*CH3CO */
    ncf[ 22 * kd + 1 ] = 2; /*C */
    ncf[ 22 * kd + 0 ] = 3; /*H */
    ncf[ 22 * kd + 2 ] = 1; /*O */

    /*CH2CHO */
    ncf[ 23 * kd + 2 ] = 1; /*O */
    ncf[ 23 * kd + 0 ] = 3; /*H */
    ncf[ 23 * kd + 1 ] = 2; /*C */

    /*CH3CHO */
    ncf[ 24 * kd + 1 ] = 2; /*C */
    ncf[ 24 * kd + 2 ] = 1; /*O */
    ncf[ 24 * kd + 0 ] = 4; /*H */

    /*C3H4-A */
    ncf[ 25 * kd + 0 ] = 4; /*H */
    ncf[ 25 * kd + 1 ] = 3; /*C */

    /*C3H6 */
    ncf[ 26 * kd + 1 ] = 3; /*C */
    ncf[ 26 * kd + 0 ] = 6; /*H */

    /*C4H6 */
    ncf[ 27 * kd + 1 ] = 4; /*C */
    ncf[ 27 * kd + 0 ] = 6; /*H */

    /*NC3H7 */
    ncf[ 28 * kd + 1 ] = 3; /*C */
    ncf[ 28 * kd + 0 ] = 7; /*H */

    /*C4H7 */
    ncf[ 29 * kd + 1 ] = 4; /*C */
    ncf[ 29 * kd + 0 ] = 7; /*H */

    /*C4H8-1 */
    ncf[ 30 * kd + 1 ] = 4; /*C */
    ncf[ 30 * kd + 0 ] = 8; /*H */

    /*PC4H9 */
    ncf[ 31 * kd + 1 ] = 4; /*C */
    ncf[ 31 * kd + 0 ] = 9; /*H */

    /*CH3COCH2 */
    ncf[ 32 * kd + 1 ] = 3; /*C */
    ncf[ 32 * kd + 0 ] = 5; /*H */
    ncf[ 32 * kd + 2 ] = 1; /*O */

    /*C2H5CHO */
    ncf[ 33 * kd + 1 ] = 3; /*C */
    ncf[ 33 * kd + 0 ] = 6; /*H */
    ncf[ 33 * kd + 2 ] = 1; /*O */

    /*C2H5CO */
    ncf[ 34 * kd + 1 ] = 3; /*C */
    ncf[ 34 * kd + 0 ] = 5; /*H */
    ncf[ 34 * kd + 2 ] = 1; /*O */

    /*C5H9 */
    ncf[ 35 * kd + 1 ] = 5; /*C */
    ncf[ 35 * kd + 0 ] = 9; /*H */

    /*C5H10-1 */
    ncf[ 36 * kd + 1 ] = 5; /*C */
    ncf[ 36 * kd + 0 ] = 10; /*H */

    /*C2H5O */
    ncf[ 37 * kd + 1 ] = 2; /*C */
    ncf[ 37 * kd + 0 ] = 5; /*H */
    ncf[ 37 * kd + 2 ] = 1; /*O */

    /*CH3O2 */
    ncf[ 38 * kd + 1 ] = 1; /*C */
    ncf[ 38 * kd + 0 ] = 3; /*H */
    ncf[ 38 * kd + 2 ] = 2; /*O */

    /*CH3O2H */
    ncf[ 39 * kd + 1 ] = 1; /*C */
    ncf[ 39 * kd + 0 ] = 4; /*H */
    ncf[ 39 * kd + 2 ] = 2; /*O */

    /*C2H3CO */
    ncf[ 40 * kd + 1 ] = 3; /*C */
    ncf[ 40 * kd + 0 ] = 3; /*H */
    ncf[ 40 * kd + 2 ] = 1; /*O */

    /*C2H3CHO */
    ncf[ 41 * kd + 1 ] = 3; /*C */
    ncf[ 41 * kd + 0 ] = 4; /*H */
    ncf[ 41 * kd + 2 ] = 1; /*O */

    /*C3H5O */
    ncf[ 42 * kd + 1 ] = 3; /*C */
    ncf[ 42 * kd + 0 ] = 5; /*H */
    ncf[ 42 * kd + 2 ] = 1; /*O */

    /*C4H7O */
    ncf[ 43 * kd + 1 ] = 4; /*C */
    ncf[ 43 * kd + 0 ] = 7; /*H */
    ncf[ 43 * kd + 2 ] = 1; /*O */

    /*C4H8OOH1-3O2 */
    ncf[ 44 * kd + 1 ] = 4; /*C */
    ncf[ 44 * kd + 0 ] = 9; /*H */
    ncf[ 44 * kd + 2 ] = 4; /*O */

    /*C4H8OOH1-3 */
    ncf[ 45 * kd + 1 ] = 4; /*C */
    ncf[ 45 * kd + 0 ] = 9; /*H */
    ncf[ 45 * kd + 2 ] = 2; /*O */

    /*PC4H9O2 */
    ncf[ 46 * kd + 1 ] = 4; /*C */
    ncf[ 46 * kd + 0 ] = 9; /*H */
    ncf[ 46 * kd + 2 ] = 2; /*O */

    /*C3H5-A */
    ncf[ 47 * kd + 1 ] = 3; /*C */
    ncf[ 47 * kd + 0 ] = 5; /*H */

    /*C3H3 */
    ncf[ 48 * kd + 1 ] = 3; /*C */
    ncf[ 48 * kd + 0 ] = 3; /*H */

    /*C3H2 */
    ncf[ 49 * kd + 0 ] = 2; /*H */
    ncf[ 49 * kd + 1 ] = 3; /*C */

    /*CH2(S) */
    ncf[ 50 * kd + 1 ] = 1; /*C */
    ncf[ 50 * kd + 0 ] = 2; /*H */

    /*NC4KET13 */
    ncf[ 51 * kd + 1 ] = 4; /*C */
    ncf[ 51 * kd + 0 ] = 8; /*H */
    ncf[ 51 * kd + 2 ] = 3; /*O */

    /*NC3H7CHO */
    ncf[ 52 * kd + 1 ] = 4; /*C */
    ncf[ 52 * kd + 0 ] = 8; /*H */
    ncf[ 52 * kd + 2 ] = 1; /*O */

    /*NC3H7CO */
    ncf[ 53 * kd + 1 ] = 4; /*C */
    ncf[ 53 * kd + 0 ] = 7; /*H */
    ncf[ 53 * kd + 2 ] = 1; /*O */

    /*C2H5COCH2 */
    ncf[ 54 * kd + 1 ] = 4; /*C */
    ncf[ 54 * kd + 0 ] = 7; /*H */
    ncf[ 54 * kd + 2 ] = 1; /*O */

    /*NC3H7COCH2 */
    ncf[ 55 * kd + 1 ] = 5; /*C */
    ncf[ 55 * kd + 0 ] = 9; /*H */
    ncf[ 55 * kd + 2 ] = 1; /*O */

    /*NC4H9CHO */
    ncf[ 56 * kd + 1 ] = 5; /*C */
    ncf[ 56 * kd + 0 ] = 10; /*H */
    ncf[ 56 * kd + 2 ] = 1; /*O */

    /*NC4H9CO */
    ncf[ 57 * kd + 1 ] = 5; /*C */
    ncf[ 57 * kd + 0 ] = 9; /*H */
    ncf[ 57 * kd + 2 ] = 1; /*O */

    /*NC7H16 */
    ncf[ 58 * kd + 1 ] = 7; /*C */
    ncf[ 58 * kd + 0 ] = 16; /*H */

    /*C7H15-1 */
    ncf[ 59 * kd + 1 ] = 7; /*C */
    ncf[ 59 * kd + 0 ] = 15; /*H */

    /*C7H15-2 */
    ncf[ 60 * kd + 1 ] = 7; /*C */
    ncf[ 60 * kd + 0 ] = 15; /*H */

    /*C7H15-3 */
    ncf[ 61 * kd + 1 ] = 7; /*C */
    ncf[ 61 * kd + 0 ] = 15; /*H */

    /*C7H15-4 */
    ncf[ 62 * kd + 1 ] = 7; /*C */
    ncf[ 62 * kd + 0 ] = 15; /*H */

    /*C7H14-2 */
    ncf[ 63 * kd + 1 ] = 7; /*C */
    ncf[ 63 * kd + 0 ] = 14; /*H */

    /*C7H14-3 */
    ncf[ 64 * kd + 1 ] = 7; /*C */
    ncf[ 64 * kd + 0 ] = 14; /*H */

    /*C7H15O2-1 */
    ncf[ 65 * kd + 1 ] = 7; /*C */
    ncf[ 65 * kd + 0 ] = 15; /*H */
    ncf[ 65 * kd + 2 ] = 2; /*O */

    /*C7H15O2-2 */
    ncf[ 66 * kd + 1 ] = 7; /*C */
    ncf[ 66 * kd + 0 ] = 15; /*H */
    ncf[ 66 * kd + 2 ] = 2; /*O */

    /*C7H15O2-3 */
    ncf[ 67 * kd + 1 ] = 7; /*C */
    ncf[ 67 * kd + 0 ] = 15; /*H */
    ncf[ 67 * kd + 2 ] = 2; /*O */

    /*C7H15O2-4 */
    ncf[ 68 * kd + 1 ] = 7; /*C */
    ncf[ 68 * kd + 0 ] = 15; /*H */
    ncf[ 68 * kd + 2 ] = 2; /*O */

    /*C7H14OOH1-3 */
    ncf[ 69 * kd + 1 ] = 7; /*C */
    ncf[ 69 * kd + 0 ] = 15; /*H */
    ncf[ 69 * kd + 2 ] = 2; /*O */

    /*C7H14OOH2-3 */
    ncf[ 70 * kd + 1 ] = 7; /*C */
    ncf[ 70 * kd + 0 ] = 15; /*H */
    ncf[ 70 * kd + 2 ] = 2; /*O */

    /*C7H14OOH2-4 */
    ncf[ 71 * kd + 1 ] = 7; /*C */
    ncf[ 71 * kd + 0 ] = 15; /*H */
    ncf[ 71 * kd + 2 ] = 2; /*O */

    /*C7H14OOH3-2 */
    ncf[ 72 * kd + 1 ] = 7; /*C */
    ncf[ 72 * kd + 0 ] = 15; /*H */
    ncf[ 72 * kd + 2 ] = 2; /*O */

    /*C7H14OOH3-4 */
    ncf[ 73 * kd + 1 ] = 7; /*C */
    ncf[ 73 * kd + 0 ] = 15; /*H */
    ncf[ 73 * kd + 2 ] = 2; /*O */

    /*C7H14OOH3-5 */
    ncf[ 74 * kd + 1 ] = 7; /*C */
    ncf[ 74 * kd + 0 ] = 15; /*H */
    ncf[ 74 * kd + 2 ] = 2; /*O */

    /*C7H14OOH4-2 */
    ncf[ 75 * kd + 1 ] = 7; /*C */
    ncf[ 75 * kd + 0 ] = 15; /*H */
    ncf[ 75 * kd + 2 ] = 2; /*O */

    /*C7H14OOH4-3 */
    ncf[ 76 * kd + 1 ] = 7; /*C */
    ncf[ 76 * kd + 0 ] = 15; /*H */
    ncf[ 76 * kd + 2 ] = 2; /*O */

    /*C7H14OOH1-3O2 */
    ncf[ 77 * kd + 1 ] = 7; /*C */
    ncf[ 77 * kd + 0 ] = 15; /*H */
    ncf[ 77 * kd + 2 ] = 4; /*O */

    /*C7H14OOH2-4O2 */
    ncf[ 78 * kd + 1 ] = 7; /*C */
    ncf[ 78 * kd + 0 ] = 15; /*H */
    ncf[ 78 * kd + 2 ] = 4; /*O */

    /*C7H14OOH3-5O2 */
    ncf[ 79 * kd + 1 ] = 7; /*C */
    ncf[ 79 * kd + 0 ] = 15; /*H */
    ncf[ 79 * kd + 2 ] = 4; /*O */

    /*C7H14OOH4-2O2 */
    ncf[ 80 * kd + 1 ] = 7; /*C */
    ncf[ 80 * kd + 0 ] = 15; /*H */
    ncf[ 80 * kd + 2 ] = 4; /*O */

    /*C7H14O1-3 */
    ncf[ 81 * kd + 1 ] = 7; /*C */
    ncf[ 81 * kd + 0 ] = 14; /*H */
    ncf[ 81 * kd + 2 ] = 1; /*O */

    /*C7H14O2-4 */
    ncf[ 82 * kd + 1 ] = 7; /*C */
    ncf[ 82 * kd + 0 ] = 14; /*H */
    ncf[ 82 * kd + 2 ] = 1; /*O */

    /*NC7KET13 */
    ncf[ 83 * kd + 1 ] = 7; /*C */
    ncf[ 83 * kd + 0 ] = 14; /*H */
    ncf[ 83 * kd + 2 ] = 3; /*O */

    /*NC7KET24 */
    ncf[ 84 * kd + 1 ] = 7; /*C */
    ncf[ 84 * kd + 0 ] = 14; /*H */
    ncf[ 84 * kd + 2 ] = 3; /*O */

    /*NC7KET35 */
    ncf[ 85 * kd + 1 ] = 7; /*C */
    ncf[ 85 * kd + 0 ] = 14; /*H */
    ncf[ 85 * kd + 2 ] = 3; /*O */

    /*NC7KET42 */
    ncf[ 86 * kd + 1 ] = 7; /*C */
    ncf[ 86 * kd + 0 ] = 14; /*H */
    ncf[ 86 * kd + 2 ] = 3; /*O */

    /*N2 */
    ncf[ 87 * kd + 3 ] = 2; /*N */

}


/* Returns the vector of strings of element names */
void CKSYME_STR(amrex::Vector<std::string>& ename)
{
    ename.resize(4);
    ename[0] = "H";
    ename[1] = "C";
    ename[2] = "O";
    ename[3] = "N";
}


/* Returns the vector of strings of species names */
void CKSYMS_STR(amrex::Vector<std::string>& kname)
{
    kname.resize(88);
    kname[0] = "H";
    kname[1] = "H2";
    kname[2] = "O";
    kname[3] = "O2";
    kname[4] = "OH";
    kname[5] = "H2O";
    kname[6] = "CO";
    kname[7] = "HCO";
    kname[8] = "CO2";
    kname[9] = "CH3";
    kname[10] = "CH4";
    kname[11] = "HO2";
    kname[12] = "H2O2";
    kname[13] = "CH2O";
    kname[14] = "CH3O";
    kname[15] = "C2H6";
    kname[16] = "C2H4";
    kname[17] = "C2H5";
    kname[18] = "C2H2";
    kname[19] = "C2H3";
    kname[20] = "CH2CO";
    kname[21] = "HCCO";
    kname[22] = "CH3CO";
    kname[23] = "CH2CHO";
    kname[24] = "CH3CHO";
    kname[25] = "C3H4-A";
    kname[26] = "C3H6";
    kname[27] = "C4H6";
    kname[28] = "NC3H7";
    kname[29] = "C4H7";
    kname[30] = "C4H8-1";
    kname[31] = "PC4H9";
    kname[32] = "CH3COCH2";
    kname[33] = "C2H5CHO";
    kname[34] = "C2H5CO";
    kname[35] = "C5H9";
    kname[36] = "C5H10-1";
    kname[37] = "C2H5O";
    kname[38] = "CH3O2";
    kname[39] = "CH3O2H";
    kname[40] = "C2H3CO";
    kname[41] = "C2H3CHO";
    kname[42] = "C3H5O";
    kname[43] = "C4H7O";
    kname[44] = "C4H8OOH1-3O2";
    kname[45] = "C4H8OOH1-3";
    kname[46] = "PC4H9O2";
    kname[47] = "C3H5-A";
    kname[48] = "C3H3";
    kname[49] = "C3H2";
    kname[50] = "CH2(S)";
    kname[51] = "NC4KET13";
    kname[52] = "NC3H7CHO";
    kname[53] = "NC3H7CO";
    kname[54] = "C2H5COCH2";
    kname[55] = "NC3H7COCH2";
    kname[56] = "NC4H9CHO";
    kname[57] = "NC4H9CO";
    kname[58] = "NC7H16";
    kname[59] = "C7H15-1";
    kname[60] = "C7H15-2";
    kname[61] = "C7H15-3";
    kname[62] = "C7H15-4";
    kname[63] = "C7H14-2";
    kname[64] = "C7H14-3";
    kname[65] = "C7H15O2-1";
    kname[66] = "C7H15O2-2";
    kname[67] = "C7H15O2-3";
    kname[68] = "C7H15O2-4";
    kname[69] = "C7H14OOH1-3";
    kname[70] = "C7H14OOH2-3";
    kname[71] = "C7H14OOH2-4";
    kname[72] = "C7H14OOH3-2";
    kname[73] = "C7H14OOH3-4";
    kname[74] = "C7H14OOH3-5";
    kname[75] = "C7H14OOH4-2";
    kname[76] = "C7H14OOH4-3";
    kname[77] = "C7H14OOH1-3O2";
    kname[78] = "C7H14OOH2-4O2";
    kname[79] = "C7H14OOH3-5O2";
    kname[80] = "C7H14OOH4-2O2";
    kname[81] = "C7H14O1-3";
    kname[82] = "C7H14O2-4";
    kname[83] = "NC7KET13";
    kname[84] = "NC7KET24";
    kname[85] = "NC7KET35";
    kname[86] = "NC7KET42";
    kname[87] = "N2";
}

/*compute the sparsity pattern of the chemistry Jacobian */
void SPARSITY_INFO( int * nJdata, const int * consP, int NCELLS)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(7921);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(88);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[7921];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int l=0; l<88; l++) {
                c_d[l] = 1.0/ 88.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<89; k++) {
        for (int l=0; l<89; l++) {
            if(J_h[ 89 * k + l] != 0.0){
                nJdata_tmp = nJdata_tmp + 1;
            }
        }
    }

    *nJdata = NCELLS * nJdata_tmp;
}



/*compute the sparsity pattern of the system Jacobian */
void SPARSITY_INFO_SYST( int * nJdata, const int * consP, int NCELLS)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(7921);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(88);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[7921];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<88; k++) {
                c_d[k] = 1.0/ 88.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<89; k++) {
        for (int l=0; l<89; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J_h[ 89 * k + l] != 0.0){
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(7921);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(88);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[7921];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<88; k++) {
                c_d[k] = 1.0/ 88.000000 ;
            }
            aJacobian_precond(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<89; k++) {
        for (int l=0; l<89; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J_h[ 89 * k + l] != 0.0){
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

    amrex::Gpu::DeviceVector<amrex::Real> J_v(7921);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(88);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[7921];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<88; k++) {
                c_d[k] = 1.0/ 88.000000 ;
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
        offset_row = nc * 89;
        offset_col = nc * 89;
        for (int k=0; k<89; k++) {
            for (int l=0; l<89; l++) {
                if(J_h[89*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(7921);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(88);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[7921];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<88; k++) {
                c_d[k] = 1.0/ 88.000000 ;
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
            offset = nc * 89;
            for (int l=0; l<89; l++) {
                for (int k=0; k<89; k++) {
                    if(J_h[89*k + l] != 0.0) {
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
            offset = nc * 89;
            for (int l=0; l<89; l++) {
                for (int k=0; k<89; k++) {
                    if(J_h[89*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(7921);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(88);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[7921];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<88; k++) {
                c_d[k] = 1.0/ 88.000000 ;
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
            offset = nc * 89;
            for (int l=0; l<89; l++) {
                for (int k=0; k<89; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp-1] = l+1 + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(J_h[89*k + l] != 0.0) {
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
            offset = nc * 89;
            for (int l=0; l<89; l++) {
                for (int k=0; k<89; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp] = l + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(J_h[89*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(7921);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(88);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[7921];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<88; k++) {
                c_d[k] = 1.0/ 88.000000 ;
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
    for (int k=0; k<89; k++) {
        for (int l=0; l<89; l++) {
            if (k == l) {
                rowVals[nJdata_tmp] = l; 
                indx[nJdata_tmp] = 89*k + l;
                nJdata_tmp = nJdata_tmp + 1; 
            } else {
                if(J_h[89*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l; 
                    indx[nJdata_tmp] = 89*k + l;
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(7921);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(88);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[7921];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<88; k++) {
                c_d[k] = 1.0/ 88.000000 ;
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
        for (int l=0; l<89; l++) {
            for (int k=0; k<89; k++) {
                if (k == l) {
                    colVals[nJdata_tmp-1] = l+1; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J_h[89*k + l] != 0.0) {
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
        for (int l=0; l<89; l++) {
            for (int k=0; k<89; k++) {
                if (k == l) {
                    colVals[nJdata_tmp] = l; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J_h[89*k + l] != 0.0) {
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
