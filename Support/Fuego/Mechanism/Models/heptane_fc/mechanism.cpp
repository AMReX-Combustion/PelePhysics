#ifndef MECHANISM_CPP
#define MECHANISM_CPP

#include "mechanism.H"



/*save atomic weights into array */
void atomicWeight(amrex::Real *  awt)
{
    awt[0] = 14.006700; /*N */
    awt[1] = 15.999400; /*O */
    awt[2] = 1.007970; /*H */
    awt[3] = 12.011150; /*C */
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
    for (id = 0; id < kd * 52; ++ id) {
         ncf[id] = 0; 
    }

    /*N2 */
    ncf[ 0 * kd + 0 ] = 2; /*N */

    /*O */
    ncf[ 1 * kd + 1 ] = 1; /*O */

    /*H2 */
    ncf[ 2 * kd + 2 ] = 2; /*H */

    /*H */
    ncf[ 3 * kd + 2 ] = 1; /*H */

    /*OH */
    ncf[ 4 * kd + 2 ] = 1; /*H */
    ncf[ 4 * kd + 1 ] = 1; /*O */

    /*H2O */
    ncf[ 5 * kd + 2 ] = 2; /*H */
    ncf[ 5 * kd + 1 ] = 1; /*O */

    /*O2 */
    ncf[ 6 * kd + 1 ] = 2; /*O */

    /*HO2 */
    ncf[ 7 * kd + 2 ] = 1; /*H */
    ncf[ 7 * kd + 1 ] = 2; /*O */

    /*H2O2 */
    ncf[ 8 * kd + 2 ] = 2; /*H */
    ncf[ 8 * kd + 1 ] = 2; /*O */

    /*CH */
    ncf[ 9 * kd + 3 ] = 1; /*C */
    ncf[ 9 * kd + 2 ] = 1; /*H */

    /*HCO */
    ncf[ 10 * kd + 2 ] = 1; /*H */
    ncf[ 10 * kd + 3 ] = 1; /*C */
    ncf[ 10 * kd + 1 ] = 1; /*O */

    /*CH2 */
    ncf[ 11 * kd + 3 ] = 1; /*C */
    ncf[ 11 * kd + 2 ] = 2; /*H */

    /*CO2 */
    ncf[ 12 * kd + 3 ] = 1; /*C */
    ncf[ 12 * kd + 1 ] = 2; /*O */

    /*CO */
    ncf[ 13 * kd + 3 ] = 1; /*C */
    ncf[ 13 * kd + 1 ] = 1; /*O */

    /*CH2O */
    ncf[ 14 * kd + 3 ] = 1; /*C */
    ncf[ 14 * kd + 2 ] = 2; /*H */
    ncf[ 14 * kd + 1 ] = 1; /*O */

    /*CH2GSG */
    ncf[ 15 * kd + 3 ] = 1; /*C */
    ncf[ 15 * kd + 2 ] = 2; /*H */

    /*CH3 */
    ncf[ 16 * kd + 3 ] = 1; /*C */
    ncf[ 16 * kd + 2 ] = 3; /*H */

    /*CH3O */
    ncf[ 17 * kd + 3 ] = 1; /*C */
    ncf[ 17 * kd + 2 ] = 3; /*H */
    ncf[ 17 * kd + 1 ] = 1; /*O */

    /*CH4 */
    ncf[ 18 * kd + 3 ] = 1; /*C */
    ncf[ 18 * kd + 2 ] = 4; /*H */

    /*CH3OH */
    ncf[ 19 * kd + 3 ] = 1; /*C */
    ncf[ 19 * kd + 2 ] = 4; /*H */
    ncf[ 19 * kd + 1 ] = 1; /*O */

    /*C2H6 */
    ncf[ 20 * kd + 3 ] = 2; /*C */
    ncf[ 20 * kd + 2 ] = 6; /*H */

    /*C2H5 */
    ncf[ 21 * kd + 3 ] = 2; /*C */
    ncf[ 21 * kd + 2 ] = 5; /*H */

    /*CH2CO */
    ncf[ 22 * kd + 3 ] = 2; /*C */
    ncf[ 22 * kd + 2 ] = 2; /*H */
    ncf[ 22 * kd + 1 ] = 1; /*O */

    /*HOCHO */
    ncf[ 23 * kd + 3 ] = 1; /*C */
    ncf[ 23 * kd + 2 ] = 2; /*H */
    ncf[ 23 * kd + 1 ] = 2; /*O */

    /*CH3O2 */
    ncf[ 24 * kd + 3 ] = 1; /*C */
    ncf[ 24 * kd + 2 ] = 3; /*H */
    ncf[ 24 * kd + 1 ] = 2; /*O */

    /*CH3O2H */
    ncf[ 25 * kd + 3 ] = 1; /*C */
    ncf[ 25 * kd + 2 ] = 4; /*H */
    ncf[ 25 * kd + 1 ] = 2; /*O */

    /*C2H2 */
    ncf[ 26 * kd + 3 ] = 2; /*C */
    ncf[ 26 * kd + 2 ] = 2; /*H */

    /*HCCO */
    ncf[ 27 * kd + 2 ] = 1; /*H */
    ncf[ 27 * kd + 3 ] = 2; /*C */
    ncf[ 27 * kd + 1 ] = 1; /*O */

    /*C2H3 */
    ncf[ 28 * kd + 3 ] = 2; /*C */
    ncf[ 28 * kd + 2 ] = 3; /*H */

    /*CH2CHO */
    ncf[ 29 * kd + 1 ] = 1; /*O */
    ncf[ 29 * kd + 2 ] = 3; /*H */
    ncf[ 29 * kd + 3 ] = 2; /*C */

    /*C3H6 */
    ncf[ 30 * kd + 3 ] = 3; /*C */
    ncf[ 30 * kd + 2 ] = 6; /*H */

    /*C2H4 */
    ncf[ 31 * kd + 3 ] = 2; /*C */
    ncf[ 31 * kd + 2 ] = 4; /*H */

    /*C2H5O */
    ncf[ 32 * kd + 3 ] = 2; /*C */
    ncf[ 32 * kd + 2 ] = 5; /*H */
    ncf[ 32 * kd + 1 ] = 1; /*O */

    /*CH3CO */
    ncf[ 33 * kd + 3 ] = 2; /*C */
    ncf[ 33 * kd + 2 ] = 3; /*H */
    ncf[ 33 * kd + 1 ] = 1; /*O */

    /*C2H5O2 */
    ncf[ 34 * kd + 3 ] = 2; /*C */
    ncf[ 34 * kd + 2 ] = 5; /*H */
    ncf[ 34 * kd + 1 ] = 2; /*O */

    /*C3H2 */
    ncf[ 35 * kd + 2 ] = 2; /*H */
    ncf[ 35 * kd + 3 ] = 3; /*C */

    /*C3H3 */
    ncf[ 36 * kd + 3 ] = 3; /*C */
    ncf[ 36 * kd + 2 ] = 3; /*H */

    /*C3H4XA */
    ncf[ 37 * kd + 2 ] = 4; /*H */
    ncf[ 37 * kd + 3 ] = 3; /*C */

    /*C3H5XA */
    ncf[ 38 * kd + 3 ] = 3; /*C */
    ncf[ 38 * kd + 2 ] = 5; /*H */

    /*NXC3H7 */
    ncf[ 39 * kd + 3 ] = 3; /*C */
    ncf[ 39 * kd + 2 ] = 7; /*H */

    /*NXC3H7O2 */
    ncf[ 40 * kd + 3 ] = 3; /*C */
    ncf[ 40 * kd + 2 ] = 7; /*H */
    ncf[ 40 * kd + 1 ] = 2; /*O */

    /*C4H6 */
    ncf[ 41 * kd + 3 ] = 4; /*C */
    ncf[ 41 * kd + 2 ] = 6; /*H */

    /*C4H7 */
    ncf[ 42 * kd + 3 ] = 4; /*C */
    ncf[ 42 * kd + 2 ] = 7; /*H */

    /*C4H8X1 */
    ncf[ 43 * kd + 3 ] = 4; /*C */
    ncf[ 43 * kd + 2 ] = 8; /*H */

    /*PXC4H9 */
    ncf[ 44 * kd + 3 ] = 4; /*C */
    ncf[ 44 * kd + 2 ] = 9; /*H */

    /*PXC4H9O2 */
    ncf[ 45 * kd + 3 ] = 4; /*C */
    ncf[ 45 * kd + 2 ] = 9; /*H */
    ncf[ 45 * kd + 1 ] = 2; /*O */

    /*C5H9 */
    ncf[ 46 * kd + 3 ] = 5; /*C */
    ncf[ 46 * kd + 2 ] = 9; /*H */

    /*C5H10X1 */
    ncf[ 47 * kd + 3 ] = 5; /*C */
    ncf[ 47 * kd + 2 ] = 10; /*H */

    /*C5H11X1 */
    ncf[ 48 * kd + 3 ] = 5; /*C */
    ncf[ 48 * kd + 2 ] = 11; /*H */

    /*C6H12X1 */
    ncf[ 49 * kd + 3 ] = 6; /*C */
    ncf[ 49 * kd + 2 ] = 12; /*H */

    /*C7H15X2 */
    ncf[ 50 * kd + 3 ] = 7; /*C */
    ncf[ 50 * kd + 2 ] = 15; /*H */

    /*NXC7H16 */
    ncf[ 51 * kd + 3 ] = 7; /*C */
    ncf[ 51 * kd + 2 ] = 16; /*H */

}


/* Returns the vector of strings of element names */
void CKSYME_STR(amrex::Vector<std::string>& ename)
{
    ename.resize(4);
    ename[0] = "N";
    ename[1] = "O";
    ename[2] = "H";
    ename[3] = "C";
}


/* Returns the vector of strings of species names */
void CKSYMS_STR(amrex::Vector<std::string>& kname)
{
    kname.resize(52);
    kname[0] = "N2";
    kname[1] = "O";
    kname[2] = "H2";
    kname[3] = "H";
    kname[4] = "OH";
    kname[5] = "H2O";
    kname[6] = "O2";
    kname[7] = "HO2";
    kname[8] = "H2O2";
    kname[9] = "CH";
    kname[10] = "HCO";
    kname[11] = "CH2";
    kname[12] = "CO2";
    kname[13] = "CO";
    kname[14] = "CH2O";
    kname[15] = "CH2GSG";
    kname[16] = "CH3";
    kname[17] = "CH3O";
    kname[18] = "CH4";
    kname[19] = "CH3OH";
    kname[20] = "C2H6";
    kname[21] = "C2H5";
    kname[22] = "CH2CO";
    kname[23] = "HOCHO";
    kname[24] = "CH3O2";
    kname[25] = "CH3O2H";
    kname[26] = "C2H2";
    kname[27] = "HCCO";
    kname[28] = "C2H3";
    kname[29] = "CH2CHO";
    kname[30] = "C3H6";
    kname[31] = "C2H4";
    kname[32] = "C2H5O";
    kname[33] = "CH3CO";
    kname[34] = "C2H5O2";
    kname[35] = "C3H2";
    kname[36] = "C3H3";
    kname[37] = "C3H4XA";
    kname[38] = "C3H5XA";
    kname[39] = "NXC3H7";
    kname[40] = "NXC3H7O2";
    kname[41] = "C4H6";
    kname[42] = "C4H7";
    kname[43] = "C4H8X1";
    kname[44] = "PXC4H9";
    kname[45] = "PXC4H9O2";
    kname[46] = "C5H9";
    kname[47] = "C5H10X1";
    kname[48] = "C5H11X1";
    kname[49] = "C6H12X1";
    kname[50] = "C7H15X2";
    kname[51] = "NXC7H16";
}

/*compute the sparsity pattern of the chemistry Jacobian */
void SPARSITY_INFO( int * nJdata, const int * consP, int NCELLS)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(2809);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(52);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[2809];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int l=0; l<52; l++) {
                c_d[l] = 1.0/ 52.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<53; k++) {
        for (int l=0; l<53; l++) {
            if(J_h[ 53 * k + l] != 0.0){
                nJdata_tmp = nJdata_tmp + 1;
            }
        }
    }

    *nJdata = NCELLS * nJdata_tmp;
}



/*compute the sparsity pattern of the system Jacobian */
void SPARSITY_INFO_SYST( int * nJdata, const int * consP, int NCELLS)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(2809);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(52);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[2809];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<52; k++) {
                c_d[k] = 1.0/ 52.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<53; k++) {
        for (int l=0; l<53; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J_h[ 53 * k + l] != 0.0){
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(2809);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(52);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[2809];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<52; k++) {
                c_d[k] = 1.0/ 52.000000 ;
            }
            aJacobian_precond(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<53; k++) {
        for (int l=0; l<53; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J_h[ 53 * k + l] != 0.0){
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

    amrex::Gpu::DeviceVector<amrex::Real> J_v(2809);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(52);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[2809];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<52; k++) {
                c_d[k] = 1.0/ 52.000000 ;
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
        offset_row = nc * 53;
        offset_col = nc * 53;
        for (int k=0; k<53; k++) {
            for (int l=0; l<53; l++) {
                if(J_h[53*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(2809);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(52);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[2809];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<52; k++) {
                c_d[k] = 1.0/ 52.000000 ;
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
            offset = nc * 53;
            for (int l=0; l<53; l++) {
                for (int k=0; k<53; k++) {
                    if(J_h[53*k + l] != 0.0) {
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
            offset = nc * 53;
            for (int l=0; l<53; l++) {
                for (int k=0; k<53; k++) {
                    if(J_h[53*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(2809);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(52);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[2809];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<52; k++) {
                c_d[k] = 1.0/ 52.000000 ;
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
            offset = nc * 53;
            for (int l=0; l<53; l++) {
                for (int k=0; k<53; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp-1] = l+1 + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(J_h[53*k + l] != 0.0) {
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
            offset = nc * 53;
            for (int l=0; l<53; l++) {
                for (int k=0; k<53; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp] = l + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(J_h[53*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(2809);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(52);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[2809];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<52; k++) {
                c_d[k] = 1.0/ 52.000000 ;
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
    for (int k=0; k<53; k++) {
        for (int l=0; l<53; l++) {
            if (k == l) {
                rowVals[nJdata_tmp] = l; 
                indx[nJdata_tmp] = 53*k + l;
                nJdata_tmp = nJdata_tmp + 1; 
            } else {
                if(J_h[53*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l; 
                    indx[nJdata_tmp] = 53*k + l;
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(2809);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(52);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[2809];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<52; k++) {
                c_d[k] = 1.0/ 52.000000 ;
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
        for (int l=0; l<53; l++) {
            for (int k=0; k<53; k++) {
                if (k == l) {
                    colVals[nJdata_tmp-1] = l+1; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J_h[53*k + l] != 0.0) {
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
        for (int l=0; l<53; l++) {
            for (int k=0; k<53; k++) {
                if (k == l) {
                    colVals[nJdata_tmp] = l; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J_h[53*k + l] != 0.0) {
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
