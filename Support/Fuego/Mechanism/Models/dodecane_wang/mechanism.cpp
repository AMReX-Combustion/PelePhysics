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
    for (id = 0; id < kd * 56; ++ id) {
         ncf[id] = 0; 
    }

    /*H */
    ncf[ 0 * kd + 1 ] = 1; /*H */

    /*O */
    ncf[ 1 * kd + 0 ] = 1; /*O */

    /*OH */
    ncf[ 2 * kd + 0 ] = 1; /*O */
    ncf[ 2 * kd + 1 ] = 1; /*H */

    /*HO2 */
    ncf[ 3 * kd + 1 ] = 1; /*H */
    ncf[ 3 * kd + 0 ] = 2; /*O */

    /*H2 */
    ncf[ 4 * kd + 1 ] = 2; /*H */

    /*H2O */
    ncf[ 5 * kd + 1 ] = 2; /*H */
    ncf[ 5 * kd + 0 ] = 1; /*O */

    /*H2O2 */
    ncf[ 6 * kd + 1 ] = 2; /*H */
    ncf[ 6 * kd + 0 ] = 2; /*O */

    /*O2 */
    ncf[ 7 * kd + 0 ] = 2; /*O */

    /*CH */
    ncf[ 8 * kd + 2 ] = 1; /*C */
    ncf[ 8 * kd + 1 ] = 1; /*H */

    /*CH2 */
    ncf[ 9 * kd + 2 ] = 1; /*C */
    ncf[ 9 * kd + 1 ] = 2; /*H */

    /*CH2* */
    ncf[ 10 * kd + 2 ] = 1; /*C */
    ncf[ 10 * kd + 1 ] = 2; /*H */

    /*CH3 */
    ncf[ 11 * kd + 2 ] = 1; /*C */
    ncf[ 11 * kd + 1 ] = 3; /*H */

    /*CH4 */
    ncf[ 12 * kd + 2 ] = 1; /*C */
    ncf[ 12 * kd + 1 ] = 4; /*H */

    /*HCO */
    ncf[ 13 * kd + 1 ] = 1; /*H */
    ncf[ 13 * kd + 2 ] = 1; /*C */
    ncf[ 13 * kd + 0 ] = 1; /*O */

    /*CH2O */
    ncf[ 14 * kd + 1 ] = 2; /*H */
    ncf[ 14 * kd + 2 ] = 1; /*C */
    ncf[ 14 * kd + 0 ] = 1; /*O */

    /*CH3O */
    ncf[ 15 * kd + 2 ] = 1; /*C */
    ncf[ 15 * kd + 1 ] = 3; /*H */
    ncf[ 15 * kd + 0 ] = 1; /*O */

    /*CH2OH */
    ncf[ 16 * kd + 2 ] = 1; /*C */
    ncf[ 16 * kd + 1 ] = 3; /*H */
    ncf[ 16 * kd + 0 ] = 1; /*O */

    /*CH3OH */
    ncf[ 17 * kd + 2 ] = 1; /*C */
    ncf[ 17 * kd + 1 ] = 4; /*H */
    ncf[ 17 * kd + 0 ] = 1; /*O */

    /*CO */
    ncf[ 18 * kd + 2 ] = 1; /*C */
    ncf[ 18 * kd + 0 ] = 1; /*O */

    /*CO2 */
    ncf[ 19 * kd + 2 ] = 1; /*C */
    ncf[ 19 * kd + 0 ] = 2; /*O */

    /*C2H */
    ncf[ 20 * kd + 2 ] = 2; /*C */
    ncf[ 20 * kd + 1 ] = 1; /*H */

    /*C2H2 */
    ncf[ 21 * kd + 2 ] = 2; /*C */
    ncf[ 21 * kd + 1 ] = 2; /*H */

    /*C2H3 */
    ncf[ 22 * kd + 2 ] = 2; /*C */
    ncf[ 22 * kd + 1 ] = 3; /*H */

    /*C2H4 */
    ncf[ 23 * kd + 2 ] = 2; /*C */
    ncf[ 23 * kd + 1 ] = 4; /*H */

    /*C2H5 */
    ncf[ 24 * kd + 2 ] = 2; /*C */
    ncf[ 24 * kd + 1 ] = 5; /*H */

    /*C2H6 */
    ncf[ 25 * kd + 2 ] = 2; /*C */
    ncf[ 25 * kd + 1 ] = 6; /*H */

    /*HCCO */
    ncf[ 26 * kd + 1 ] = 1; /*H */
    ncf[ 26 * kd + 2 ] = 2; /*C */
    ncf[ 26 * kd + 0 ] = 1; /*O */

    /*CH2CO */
    ncf[ 27 * kd + 2 ] = 2; /*C */
    ncf[ 27 * kd + 1 ] = 2; /*H */
    ncf[ 27 * kd + 0 ] = 1; /*O */

    /*CH3CO */
    ncf[ 28 * kd + 2 ] = 2; /*C */
    ncf[ 28 * kd + 1 ] = 3; /*H */
    ncf[ 28 * kd + 0 ] = 1; /*O */

    /*CH2CHO */
    ncf[ 29 * kd + 0 ] = 1; /*O */
    ncf[ 29 * kd + 1 ] = 3; /*H */
    ncf[ 29 * kd + 2 ] = 2; /*C */

    /*CH3CHO */
    ncf[ 30 * kd + 2 ] = 2; /*C */
    ncf[ 30 * kd + 1 ] = 4; /*H */
    ncf[ 30 * kd + 0 ] = 1; /*O */

    /*C3H3 */
    ncf[ 31 * kd + 2 ] = 3; /*C */
    ncf[ 31 * kd + 1 ] = 3; /*H */

    /*pC3H4 */
    ncf[ 32 * kd + 1 ] = 4; /*H */
    ncf[ 32 * kd + 2 ] = 3; /*C */

    /*aC3H4 */
    ncf[ 33 * kd + 2 ] = 3; /*C */
    ncf[ 33 * kd + 1 ] = 4; /*H */

    /*aC3H5 */
    ncf[ 34 * kd + 2 ] = 3; /*C */
    ncf[ 34 * kd + 1 ] = 5; /*H */

    /*CH3CCH2 */
    ncf[ 35 * kd + 2 ] = 3; /*C */
    ncf[ 35 * kd + 1 ] = 5; /*H */

    /*C3H6 */
    ncf[ 36 * kd + 2 ] = 3; /*C */
    ncf[ 36 * kd + 1 ] = 6; /*H */

    /*nC3H7 */
    ncf[ 37 * kd + 2 ] = 3; /*C */
    ncf[ 37 * kd + 1 ] = 7; /*H */

    /*iC3H7 */
    ncf[ 38 * kd + 2 ] = 3; /*C */
    ncf[ 38 * kd + 1 ] = 7; /*H */

    /*C2H3CHO */
    ncf[ 39 * kd + 2 ] = 3; /*C */
    ncf[ 39 * kd + 1 ] = 4; /*H */
    ncf[ 39 * kd + 0 ] = 1; /*O */

    /*C4H2 */
    ncf[ 40 * kd + 2 ] = 4; /*C */
    ncf[ 40 * kd + 1 ] = 2; /*H */

    /*iC4H3 */
    ncf[ 41 * kd + 2 ] = 4; /*C */
    ncf[ 41 * kd + 1 ] = 3; /*H */

    /*C4H4 */
    ncf[ 42 * kd + 2 ] = 4; /*C */
    ncf[ 42 * kd + 1 ] = 4; /*H */

    /*iC4H5 */
    ncf[ 43 * kd + 2 ] = 4; /*C */
    ncf[ 43 * kd + 1 ] = 5; /*H */

    /*C4H5-2 */
    ncf[ 44 * kd + 2 ] = 4; /*C */
    ncf[ 44 * kd + 1 ] = 5; /*H */

    /*C4H6 */
    ncf[ 45 * kd + 2 ] = 4; /*C */
    ncf[ 45 * kd + 1 ] = 6; /*H */

    /*C4H612 */
    ncf[ 46 * kd + 2 ] = 4; /*C */
    ncf[ 46 * kd + 1 ] = 6; /*H */

    /*C4H6-2 */
    ncf[ 47 * kd + 2 ] = 4; /*C */
    ncf[ 47 * kd + 1 ] = 6; /*H */

    /*C4H7 */
    ncf[ 48 * kd + 2 ] = 4; /*C */
    ncf[ 48 * kd + 1 ] = 7; /*H */

    /*C4H81 */
    ncf[ 49 * kd + 2 ] = 4; /*C */
    ncf[ 49 * kd + 1 ] = 8; /*H */

    /*pC4H9 */
    ncf[ 50 * kd + 2 ] = 4; /*C */
    ncf[ 50 * kd + 1 ] = 9; /*H */

    /*NC12H26 */
    ncf[ 51 * kd + 2 ] = 12; /*C */
    ncf[ 51 * kd + 1 ] = 26; /*H */

    /*C6H12 */
    ncf[ 52 * kd + 2 ] = 6; /*C */
    ncf[ 52 * kd + 1 ] = 12; /*H */

    /*C6H11 */
    ncf[ 53 * kd + 2 ] = 6; /*C */
    ncf[ 53 * kd + 1 ] = 11; /*H */

    /*C5H10 */
    ncf[ 54 * kd + 2 ] = 5; /*C */
    ncf[ 54 * kd + 1 ] = 10; /*H */

    /*N2 */
    ncf[ 55 * kd + 3 ] = 2; /*N */

}


/* Returns the vector of strings of element names */
void CKSYME_STR(amrex::Vector<std::string>& ename)
{
    ename.resize(4);
    ename[0] = "O";
    ename[1] = "H";
    ename[2] = "C";
    ename[3] = "N";
}


/* Returns the vector of strings of species names */
void CKSYMS_STR(amrex::Vector<std::string>& kname)
{
    kname.resize(56);
    kname[0] = "H";
    kname[1] = "O";
    kname[2] = "OH";
    kname[3] = "HO2";
    kname[4] = "H2";
    kname[5] = "H2O";
    kname[6] = "H2O2";
    kname[7] = "O2";
    kname[8] = "CH";
    kname[9] = "CH2";
    kname[10] = "CH2*";
    kname[11] = "CH3";
    kname[12] = "CH4";
    kname[13] = "HCO";
    kname[14] = "CH2O";
    kname[15] = "CH3O";
    kname[16] = "CH2OH";
    kname[17] = "CH3OH";
    kname[18] = "CO";
    kname[19] = "CO2";
    kname[20] = "C2H";
    kname[21] = "C2H2";
    kname[22] = "C2H3";
    kname[23] = "C2H4";
    kname[24] = "C2H5";
    kname[25] = "C2H6";
    kname[26] = "HCCO";
    kname[27] = "CH2CO";
    kname[28] = "CH3CO";
    kname[29] = "CH2CHO";
    kname[30] = "CH3CHO";
    kname[31] = "C3H3";
    kname[32] = "pC3H4";
    kname[33] = "aC3H4";
    kname[34] = "aC3H5";
    kname[35] = "CH3CCH2";
    kname[36] = "C3H6";
    kname[37] = "nC3H7";
    kname[38] = "iC3H7";
    kname[39] = "C2H3CHO";
    kname[40] = "C4H2";
    kname[41] = "iC4H3";
    kname[42] = "C4H4";
    kname[43] = "iC4H5";
    kname[44] = "C4H5-2";
    kname[45] = "C4H6";
    kname[46] = "C4H612";
    kname[47] = "C4H6-2";
    kname[48] = "C4H7";
    kname[49] = "C4H81";
    kname[50] = "pC4H9";
    kname[51] = "NC12H26";
    kname[52] = "C6H12";
    kname[53] = "C6H11";
    kname[54] = "C5H10";
    kname[55] = "N2";
}

/*compute the sparsity pattern of the chemistry Jacobian */
void SPARSITY_INFO( int * nJdata, const int * consP, int NCELLS)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(3249);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(56);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[3249];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int l=0; l<56; l++) {
                c_d[l] = 1.0/ 56.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<57; k++) {
        for (int l=0; l<57; l++) {
            if(J_h[ 57 * k + l] != 0.0){
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(3249);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(56);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[3249];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<56; k++) {
                c_d[k] = 1.0/ 56.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<57; k++) {
        for (int l=0; l<57; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J_h[ 57 * k + l] != 0.0){
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(3249);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(56);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[3249];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<56; k++) {
                c_d[k] = 1.0/ 56.000000 ;
            }
            aJacobian_precond(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<57; k++) {
        for (int l=0; l<57; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J_h[ 57 * k + l] != 0.0){
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

    amrex::Gpu::DeviceVector<amrex::Real> J_v(3249);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(56);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[3249];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<56; k++) {
                c_d[k] = 1.0/ 56.000000 ;
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
        offset_row = nc * 57;
        offset_col = nc * 57;
        for (int k=0; k<57; k++) {
            for (int l=0; l<57; l++) {
                if(J_h[57*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(3249);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(56);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[3249];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<56; k++) {
                c_d[k] = 1.0/ 56.000000 ;
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
            offset = nc * 57;
            for (int l=0; l<57; l++) {
                for (int k=0; k<57; k++) {
                    if(J_h[57*k + l] != 0.0) {
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
            offset = nc * 57;
            for (int l=0; l<57; l++) {
                for (int k=0; k<57; k++) {
                    if(J_h[57*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(3249);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(56);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[3249];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<56; k++) {
                c_d[k] = 1.0/ 56.000000 ;
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
            offset = nc * 57;
            for (int l=0; l<57; l++) {
                for (int k=0; k<57; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp-1] = l+1 + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(J_h[57*k + l] != 0.0) {
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
            offset = nc * 57;
            for (int l=0; l<57; l++) {
                for (int k=0; k<57; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp] = l + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(J_h[57*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(3249);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(56);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[3249];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<56; k++) {
                c_d[k] = 1.0/ 56.000000 ;
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
    for (int k=0; k<57; k++) {
        for (int l=0; l<57; l++) {
            if (k == l) {
                rowVals[nJdata_tmp] = l; 
                indx[nJdata_tmp] = 57*k + l;
                nJdata_tmp = nJdata_tmp + 1; 
            } else {
                if(J_h[57*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l; 
                    indx[nJdata_tmp] = 57*k + l;
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(3249);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(56);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[3249];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<56; k++) {
                c_d[k] = 1.0/ 56.000000 ;
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
        for (int l=0; l<57; l++) {
            for (int k=0; k<57; k++) {
                if (k == l) {
                    colVals[nJdata_tmp-1] = l+1; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J_h[57*k + l] != 0.0) {
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
        for (int l=0; l<57; l++) {
            for (int k=0; k<57; k++) {
                if (k == l) {
                    colVals[nJdata_tmp] = l; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J_h[57*k + l] != 0.0) {
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
