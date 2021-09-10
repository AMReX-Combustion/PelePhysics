#ifndef MECHANISM_CPP
#define MECHANISM_CPP

#include "mechanism.H"



/*save atomic weights into array */
void atomicWeight(amrex::Real *  awt)
{
    awt[0] = 15.999400; /*O */
    awt[1] = 1.007970; /*H */
    awt[2] = 12.011150; /*C */
    awt[3] = 14.006700; /*N */
    awt[4] = 39.948000; /*AR */
    awt[5] = 18.998400; /*F */
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
    for (id = 0; id < kd * 72; ++ id) {
         ncf[id] = 0; 
    }

    /*H2 */
    ncf[ 0 * kd + 1 ] = 2; /*H */

    /*H */
    ncf[ 1 * kd + 1 ] = 1; /*H */

    /*O */
    ncf[ 2 * kd + 0 ] = 1; /*O */

    /*O2 */
    ncf[ 3 * kd + 0 ] = 2; /*O */

    /*OH */
    ncf[ 4 * kd + 0 ] = 1; /*O */
    ncf[ 4 * kd + 1 ] = 1; /*H */

    /*H2O */
    ncf[ 5 * kd + 1 ] = 2; /*H */
    ncf[ 5 * kd + 0 ] = 1; /*O */

    /*HO2 */
    ncf[ 6 * kd + 1 ] = 1; /*H */
    ncf[ 6 * kd + 0 ] = 2; /*O */

    /*H2O2 */
    ncf[ 7 * kd + 1 ] = 2; /*H */
    ncf[ 7 * kd + 0 ] = 2; /*O */

    /*C */
    ncf[ 8 * kd + 2 ] = 1; /*C */

    /*CH */
    ncf[ 9 * kd + 2 ] = 1; /*C */
    ncf[ 9 * kd + 1 ] = 1; /*H */

    /*CH2 */
    ncf[ 10 * kd + 2 ] = 1; /*C */
    ncf[ 10 * kd + 1 ] = 2; /*H */

    /*CH2(S) */
    ncf[ 11 * kd + 2 ] = 1; /*C */
    ncf[ 11 * kd + 1 ] = 2; /*H */

    /*CH3 */
    ncf[ 12 * kd + 2 ] = 1; /*C */
    ncf[ 12 * kd + 1 ] = 3; /*H */

    /*CH4 */
    ncf[ 13 * kd + 2 ] = 1; /*C */
    ncf[ 13 * kd + 1 ] = 4; /*H */

    /*CO */
    ncf[ 14 * kd + 2 ] = 1; /*C */
    ncf[ 14 * kd + 0 ] = 1; /*O */

    /*CO2 */
    ncf[ 15 * kd + 2 ] = 1; /*C */
    ncf[ 15 * kd + 0 ] = 2; /*O */

    /*HCO */
    ncf[ 16 * kd + 1 ] = 1; /*H */
    ncf[ 16 * kd + 2 ] = 1; /*C */
    ncf[ 16 * kd + 0 ] = 1; /*O */

    /*CH2O */
    ncf[ 17 * kd + 1 ] = 2; /*H */
    ncf[ 17 * kd + 2 ] = 1; /*C */
    ncf[ 17 * kd + 0 ] = 1; /*O */

    /*CH2OH */
    ncf[ 18 * kd + 2 ] = 1; /*C */
    ncf[ 18 * kd + 1 ] = 3; /*H */
    ncf[ 18 * kd + 0 ] = 1; /*O */

    /*CH3O */
    ncf[ 19 * kd + 2 ] = 1; /*C */
    ncf[ 19 * kd + 1 ] = 3; /*H */
    ncf[ 19 * kd + 0 ] = 1; /*O */

    /*CH3OH */
    ncf[ 20 * kd + 2 ] = 1; /*C */
    ncf[ 20 * kd + 1 ] = 4; /*H */
    ncf[ 20 * kd + 0 ] = 1; /*O */

    /*C2H */
    ncf[ 21 * kd + 2 ] = 2; /*C */
    ncf[ 21 * kd + 1 ] = 1; /*H */

    /*C2H2 */
    ncf[ 22 * kd + 2 ] = 2; /*C */
    ncf[ 22 * kd + 1 ] = 2; /*H */

    /*C2H3 */
    ncf[ 23 * kd + 2 ] = 2; /*C */
    ncf[ 23 * kd + 1 ] = 3; /*H */

    /*C2H4 */
    ncf[ 24 * kd + 2 ] = 2; /*C */
    ncf[ 24 * kd + 1 ] = 4; /*H */

    /*C2H5 */
    ncf[ 25 * kd + 2 ] = 2; /*C */
    ncf[ 25 * kd + 1 ] = 5; /*H */

    /*C2H6 */
    ncf[ 26 * kd + 2 ] = 2; /*C */
    ncf[ 26 * kd + 1 ] = 6; /*H */

    /*HCCO */
    ncf[ 27 * kd + 1 ] = 1; /*H */
    ncf[ 27 * kd + 2 ] = 2; /*C */
    ncf[ 27 * kd + 0 ] = 1; /*O */

    /*CH2CO */
    ncf[ 28 * kd + 2 ] = 2; /*C */
    ncf[ 28 * kd + 1 ] = 2; /*H */
    ncf[ 28 * kd + 0 ] = 1; /*O */

    /*HCCOH */
    ncf[ 29 * kd + 2 ] = 2; /*C */
    ncf[ 29 * kd + 0 ] = 1; /*O */
    ncf[ 29 * kd + 1 ] = 2; /*H */

    /*N */
    ncf[ 30 * kd + 3 ] = 1; /*N */

    /*NH */
    ncf[ 31 * kd + 3 ] = 1; /*N */
    ncf[ 31 * kd + 1 ] = 1; /*H */

    /*NH2 */
    ncf[ 32 * kd + 3 ] = 1; /*N */
    ncf[ 32 * kd + 1 ] = 2; /*H */

    /*NH3 */
    ncf[ 33 * kd + 3 ] = 1; /*N */
    ncf[ 33 * kd + 1 ] = 3; /*H */

    /*NNH */
    ncf[ 34 * kd + 3 ] = 2; /*N */
    ncf[ 34 * kd + 1 ] = 1; /*H */

    /*NO */
    ncf[ 35 * kd + 3 ] = 1; /*N */
    ncf[ 35 * kd + 0 ] = 1; /*O */

    /*NO2 */
    ncf[ 36 * kd + 3 ] = 1; /*N */
    ncf[ 36 * kd + 0 ] = 2; /*O */

    /*N2O */
    ncf[ 37 * kd + 3 ] = 2; /*N */
    ncf[ 37 * kd + 0 ] = 1; /*O */

    /*HNO */
    ncf[ 38 * kd + 1 ] = 1; /*H */
    ncf[ 38 * kd + 3 ] = 1; /*N */
    ncf[ 38 * kd + 0 ] = 1; /*O */

    /*CN */
    ncf[ 39 * kd + 2 ] = 1; /*C */
    ncf[ 39 * kd + 3 ] = 1; /*N */

    /*HCN */
    ncf[ 40 * kd + 1 ] = 1; /*H */
    ncf[ 40 * kd + 2 ] = 1; /*C */
    ncf[ 40 * kd + 3 ] = 1; /*N */

    /*H2CN */
    ncf[ 41 * kd + 1 ] = 2; /*H */
    ncf[ 41 * kd + 2 ] = 1; /*C */
    ncf[ 41 * kd + 3 ] = 1; /*N */

    /*HCNN */
    ncf[ 42 * kd + 2 ] = 1; /*C */
    ncf[ 42 * kd + 3 ] = 2; /*N */
    ncf[ 42 * kd + 1 ] = 1; /*H */

    /*HCNO */
    ncf[ 43 * kd + 1 ] = 1; /*H */
    ncf[ 43 * kd + 3 ] = 1; /*N */
    ncf[ 43 * kd + 2 ] = 1; /*C */
    ncf[ 43 * kd + 0 ] = 1; /*O */

    /*HOCN */
    ncf[ 44 * kd + 1 ] = 1; /*H */
    ncf[ 44 * kd + 3 ] = 1; /*N */
    ncf[ 44 * kd + 2 ] = 1; /*C */
    ncf[ 44 * kd + 0 ] = 1; /*O */

    /*HNCO */
    ncf[ 45 * kd + 1 ] = 1; /*H */
    ncf[ 45 * kd + 3 ] = 1; /*N */
    ncf[ 45 * kd + 2 ] = 1; /*C */
    ncf[ 45 * kd + 0 ] = 1; /*O */

    /*NCO */
    ncf[ 46 * kd + 3 ] = 1; /*N */
    ncf[ 46 * kd + 2 ] = 1; /*C */
    ncf[ 46 * kd + 0 ] = 1; /*O */

    /*N2 */
    ncf[ 47 * kd + 3 ] = 2; /*N */

    /*AR */
    ncf[ 48 * kd + 4 ] = 1; /*AR */

    /*C3H7 */
    ncf[ 49 * kd + 2 ] = 3; /*C */
    ncf[ 49 * kd + 1 ] = 7; /*H */

    /*C3H8 */
    ncf[ 50 * kd + 2 ] = 3; /*C */
    ncf[ 50 * kd + 1 ] = 8; /*H */

    /*CH2CHO */
    ncf[ 51 * kd + 0 ] = 1; /*O */
    ncf[ 51 * kd + 1 ] = 3; /*H */
    ncf[ 51 * kd + 2 ] = 2; /*C */

    /*CH3CHO */
    ncf[ 52 * kd + 2 ] = 2; /*C */
    ncf[ 52 * kd + 1 ] = 4; /*H */
    ncf[ 52 * kd + 0 ] = 1; /*O */

    /*HF */
    ncf[ 53 * kd + 1 ] = 1; /*H */
    ncf[ 53 * kd + 5 ] = 1; /*F */

    /*F */
    ncf[ 54 * kd + 5 ] = 1; /*F */

    /*F2 */
    ncf[ 55 * kd + 5 ] = 2; /*F */

    /*CH3F */
    ncf[ 56 * kd + 2 ] = 1; /*C */
    ncf[ 56 * kd + 1 ] = 3; /*H */
    ncf[ 56 * kd + 5 ] = 1; /*F */

    /*CH2F2 */
    ncf[ 57 * kd + 2 ] = 1; /*C */
    ncf[ 57 * kd + 1 ] = 2; /*H */
    ncf[ 57 * kd + 5 ] = 2; /*F */

    /*CHF3 */
    ncf[ 58 * kd + 2 ] = 1; /*C */
    ncf[ 58 * kd + 1 ] = 1; /*H */
    ncf[ 58 * kd + 5 ] = 3; /*F */

    /*CF4 */
    ncf[ 59 * kd + 2 ] = 1; /*C */
    ncf[ 59 * kd + 5 ] = 4; /*F */

    /*CH2F */
    ncf[ 60 * kd + 2 ] = 1; /*C */
    ncf[ 60 * kd + 1 ] = 2; /*H */
    ncf[ 60 * kd + 5 ] = 1; /*F */

    /*CHF2 */
    ncf[ 61 * kd + 2 ] = 1; /*C */
    ncf[ 61 * kd + 1 ] = 1; /*H */
    ncf[ 61 * kd + 5 ] = 2; /*F */

    /*CF3 */
    ncf[ 62 * kd + 2 ] = 1; /*C */
    ncf[ 62 * kd + 5 ] = 3; /*F */

    /*CHF */
    ncf[ 63 * kd + 2 ] = 1; /*C */
    ncf[ 63 * kd + 1 ] = 1; /*H */
    ncf[ 63 * kd + 5 ] = 1; /*F */

    /*CF2 */
    ncf[ 64 * kd + 2 ] = 1; /*C */
    ncf[ 64 * kd + 5 ] = 2; /*F */

    /*CF */
    ncf[ 65 * kd + 2 ] = 1; /*C */
    ncf[ 65 * kd + 5 ] = 1; /*F */

    /*CF3O */
    ncf[ 66 * kd + 2 ] = 1; /*C */
    ncf[ 66 * kd + 0 ] = 1; /*O */
    ncf[ 66 * kd + 5 ] = 3; /*F */

    /*CF2O */
    ncf[ 67 * kd + 2 ] = 1; /*C */
    ncf[ 67 * kd + 0 ] = 1; /*O */
    ncf[ 67 * kd + 5 ] = 2; /*F */

    /*CHFO */
    ncf[ 68 * kd + 2 ] = 1; /*C */
    ncf[ 68 * kd + 1 ] = 1; /*H */
    ncf[ 68 * kd + 0 ] = 1; /*O */
    ncf[ 68 * kd + 5 ] = 1; /*F */

    /*CFO */
    ncf[ 69 * kd + 2 ] = 1; /*C */
    ncf[ 69 * kd + 0 ] = 1; /*O */
    ncf[ 69 * kd + 5 ] = 1; /*F */

    /*iC3H7 */
    ncf[ 70 * kd + 2 ] = 3; /*C */
    ncf[ 70 * kd + 1 ] = 7; /*H */

    /*nC3H7 */
    ncf[ 71 * kd + 2 ] = 3; /*C */
    ncf[ 71 * kd + 1 ] = 7; /*H */

}


/* Returns the vector of strings of element names */
void CKSYME_STR(amrex::Vector<std::string>& ename)
{
    ename.resize(6);
    ename[0] = "O";
    ename[1] = "H";
    ename[2] = "C";
    ename[3] = "N";
    ename[4] = "AR";
    ename[5] = "F";
}


/* Returns the vector of strings of species names */
void CKSYMS_STR(amrex::Vector<std::string>& kname)
{
    kname.resize(72);
    kname[0] = "H2";
    kname[1] = "H";
    kname[2] = "O";
    kname[3] = "O2";
    kname[4] = "OH";
    kname[5] = "H2O";
    kname[6] = "HO2";
    kname[7] = "H2O2";
    kname[8] = "C";
    kname[9] = "CH";
    kname[10] = "CH2";
    kname[11] = "CH2(S)";
    kname[12] = "CH3";
    kname[13] = "CH4";
    kname[14] = "CO";
    kname[15] = "CO2";
    kname[16] = "HCO";
    kname[17] = "CH2O";
    kname[18] = "CH2OH";
    kname[19] = "CH3O";
    kname[20] = "CH3OH";
    kname[21] = "C2H";
    kname[22] = "C2H2";
    kname[23] = "C2H3";
    kname[24] = "C2H4";
    kname[25] = "C2H5";
    kname[26] = "C2H6";
    kname[27] = "HCCO";
    kname[28] = "CH2CO";
    kname[29] = "HCCOH";
    kname[30] = "N";
    kname[31] = "NH";
    kname[32] = "NH2";
    kname[33] = "NH3";
    kname[34] = "NNH";
    kname[35] = "NO";
    kname[36] = "NO2";
    kname[37] = "N2O";
    kname[38] = "HNO";
    kname[39] = "CN";
    kname[40] = "HCN";
    kname[41] = "H2CN";
    kname[42] = "HCNN";
    kname[43] = "HCNO";
    kname[44] = "HOCN";
    kname[45] = "HNCO";
    kname[46] = "NCO";
    kname[47] = "N2";
    kname[48] = "AR";
    kname[49] = "C3H7";
    kname[50] = "C3H8";
    kname[51] = "CH2CHO";
    kname[52] = "CH3CHO";
    kname[53] = "HF";
    kname[54] = "F";
    kname[55] = "F2";
    kname[56] = "CH3F";
    kname[57] = "CH2F2";
    kname[58] = "CHF3";
    kname[59] = "CF4";
    kname[60] = "CH2F";
    kname[61] = "CHF2";
    kname[62] = "CF3";
    kname[63] = "CHF";
    kname[64] = "CF2";
    kname[65] = "CF";
    kname[66] = "CF3O";
    kname[67] = "CF2O";
    kname[68] = "CHFO";
    kname[69] = "CFO";
    kname[70] = "iC3H7";
    kname[71] = "nC3H7";
}

/*compute the sparsity pattern of the chemistry Jacobian */
void SPARSITY_INFO( int * nJdata, const int * consP, int NCELLS)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(5329);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(72);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[5329];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int l=0; l<72; l++) {
                c_d[l] = 1.0/ 72.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<73; k++) {
        for (int l=0; l<73; l++) {
            if(J_h[ 73 * k + l] != 0.0){
                nJdata_tmp = nJdata_tmp + 1;
            }
        }
    }

    *nJdata = NCELLS * nJdata_tmp;
}
#endif



/*compute the sparsity pattern of the system Jacobian */
void SPARSITY_INFO_SYST( int * nJdata, const int * consP, int NCELLS)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(5329);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(72);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[5329];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<72; k++) {
                c_d[k] = 1.0/ 72.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<73; k++) {
        for (int l=0; l<73; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J_h[ 73 * k + l] != 0.0){
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(5329);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(72);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[5329];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<72; k++) {
                c_d[k] = 1.0/ 72.000000 ;
            }
            aJacobian_precond(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<73; k++) {
        for (int l=0; l<73; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J_h[ 73 * k + l] != 0.0){
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

    amrex::Gpu::DeviceVector<amrex::Real> J_v(5329);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(72);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[5329];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<72; k++) {
                c_d[k] = 1.0/ 72.000000 ;
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
        offset_row = nc * 73;
        offset_col = nc * 73;
        for (int k=0; k<73; k++) {
            for (int l=0; l<73; l++) {
                if(J_h[73*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(5329);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(72);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[5329];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<72; k++) {
                c_d[k] = 1.0/ 72.000000 ;
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
            offset = nc * 73;
            for (int l=0; l<73; l++) {
                for (int k=0; k<73; k++) {
                    if(J_h[73*k + l] != 0.0) {
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
            offset = nc * 73;
            for (int l=0; l<73; l++) {
                for (int k=0; k<73; k++) {
                    if(J_h[73*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(5329);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(72);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[5329];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<72; k++) {
                c_d[k] = 1.0/ 72.000000 ;
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
            offset = nc * 73;
            for (int l=0; l<73; l++) {
                for (int k=0; k<73; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp-1] = l+1 + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(J_h[73*k + l] != 0.0) {
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
            offset = nc * 73;
            for (int l=0; l<73; l++) {
                for (int k=0; k<73; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp] = l + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(J_h[73*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(5329);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(72);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[5329];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<72; k++) {
                c_d[k] = 1.0/ 72.000000 ;
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
    for (int k=0; k<73; k++) {
        for (int l=0; l<73; l++) {
            if (k == l) {
                rowVals[nJdata_tmp] = l; 
                indx[nJdata_tmp] = 73*k + l;
                nJdata_tmp = nJdata_tmp + 1; 
            } else {
                if(J_h[73*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l; 
                    indx[nJdata_tmp] = 73*k + l;
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(5329);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(72);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[5329];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<72; k++) {
                c_d[k] = 1.0/ 72.000000 ;
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
        for (int l=0; l<73; l++) {
            for (int k=0; k<73; k++) {
                if (k == l) {
                    colVals[nJdata_tmp-1] = l+1; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J_h[73*k + l] != 0.0) {
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
        for (int l=0; l<73; l++) {
            for (int k=0; k<73; k++) {
                if (k == l) {
                    colVals[nJdata_tmp] = l; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J_h[73*k + l] != 0.0) {
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
