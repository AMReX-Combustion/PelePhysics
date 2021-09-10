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
    for (id = 0; id < kd * 34; ++ id) {
         ncf[id] = 0; 
    }

    /*N2 */
    ncf[ 0 * kd + 0 ] = 2; /*N */

    /*O */
    ncf[ 1 * kd + 1 ] = 1; /*O */

    /*H */
    ncf[ 2 * kd + 2 ] = 1; /*H */

    /*OH */
    ncf[ 3 * kd + 2 ] = 1; /*H */
    ncf[ 3 * kd + 1 ] = 1; /*O */

    /*H2 */
    ncf[ 4 * kd + 2 ] = 2; /*H */

    /*O2 */
    ncf[ 5 * kd + 1 ] = 2; /*O */

    /*H2O */
    ncf[ 6 * kd + 2 ] = 2; /*H */
    ncf[ 6 * kd + 1 ] = 1; /*O */

    /*H2O2 */
    ncf[ 7 * kd + 2 ] = 2; /*H */
    ncf[ 7 * kd + 1 ] = 2; /*O */

    /*HO2 */
    ncf[ 8 * kd + 2 ] = 1; /*H */
    ncf[ 8 * kd + 1 ] = 2; /*O */

    /*CH2GSG */
    ncf[ 9 * kd + 3 ] = 1; /*C */
    ncf[ 9 * kd + 2 ] = 2; /*H */

    /*CH2O */
    ncf[ 10 * kd + 3 ] = 1; /*C */
    ncf[ 10 * kd + 2 ] = 2; /*H */
    ncf[ 10 * kd + 1 ] = 1; /*O */

    /*CH3 */
    ncf[ 11 * kd + 3 ] = 1; /*C */
    ncf[ 11 * kd + 2 ] = 3; /*H */

    /*CO */
    ncf[ 12 * kd + 3 ] = 1; /*C */
    ncf[ 12 * kd + 1 ] = 1; /*O */

    /*CH3O */
    ncf[ 13 * kd + 3 ] = 1; /*C */
    ncf[ 13 * kd + 2 ] = 3; /*H */
    ncf[ 13 * kd + 1 ] = 1; /*O */

    /*C2H5 */
    ncf[ 14 * kd + 3 ] = 2; /*C */
    ncf[ 14 * kd + 2 ] = 5; /*H */

    /*CH4 */
    ncf[ 15 * kd + 3 ] = 1; /*C */
    ncf[ 15 * kd + 2 ] = 4; /*H */

    /*C2H4 */
    ncf[ 16 * kd + 3 ] = 2; /*C */
    ncf[ 16 * kd + 2 ] = 4; /*H */

    /*C2H6 */
    ncf[ 17 * kd + 3 ] = 2; /*C */
    ncf[ 17 * kd + 2 ] = 6; /*H */

    /*CO2 */
    ncf[ 18 * kd + 3 ] = 1; /*C */
    ncf[ 18 * kd + 1 ] = 2; /*O */

    /*HCO */
    ncf[ 19 * kd + 2 ] = 1; /*H */
    ncf[ 19 * kd + 3 ] = 1; /*C */
    ncf[ 19 * kd + 1 ] = 1; /*O */

    /*CH3O2 */
    ncf[ 20 * kd + 3 ] = 1; /*C */
    ncf[ 20 * kd + 2 ] = 3; /*H */
    ncf[ 20 * kd + 1 ] = 2; /*O */

    /*CH3O2H */
    ncf[ 21 * kd + 3 ] = 1; /*C */
    ncf[ 21 * kd + 2 ] = 4; /*H */
    ncf[ 21 * kd + 1 ] = 2; /*O */

    /*C2H2 */
    ncf[ 22 * kd + 3 ] = 2; /*C */
    ncf[ 22 * kd + 2 ] = 2; /*H */

    /*HCCO */
    ncf[ 23 * kd + 2 ] = 1; /*H */
    ncf[ 23 * kd + 3 ] = 2; /*C */
    ncf[ 23 * kd + 1 ] = 1; /*O */

    /*C2H3 */
    ncf[ 24 * kd + 3 ] = 2; /*C */
    ncf[ 24 * kd + 2 ] = 3; /*H */

    /*CH2CHO */
    ncf[ 25 * kd + 1 ] = 1; /*O */
    ncf[ 25 * kd + 2 ] = 3; /*H */
    ncf[ 25 * kd + 3 ] = 2; /*C */

    /*C3H5XA */
    ncf[ 26 * kd + 3 ] = 3; /*C */
    ncf[ 26 * kd + 2 ] = 5; /*H */

    /*C3H6 */
    ncf[ 27 * kd + 3 ] = 3; /*C */
    ncf[ 27 * kd + 2 ] = 6; /*H */

    /*C3H5O */
    ncf[ 28 * kd + 3 ] = 3; /*C */
    ncf[ 28 * kd + 2 ] = 5; /*H */
    ncf[ 28 * kd + 1 ] = 1; /*O */

    /*IXC3H7 */
    ncf[ 29 * kd + 3 ] = 3; /*C */
    ncf[ 29 * kd + 2 ] = 7; /*H */

    /*NXC3H7 */
    ncf[ 30 * kd + 3 ] = 3; /*C */
    ncf[ 30 * kd + 2 ] = 7; /*H */

    /*C3H8 */
    ncf[ 31 * kd + 3 ] = 3; /*C */
    ncf[ 31 * kd + 2 ] = 8; /*H */

    /*IXC3H7O2 */
    ncf[ 32 * kd + 3 ] = 3; /*C */
    ncf[ 32 * kd + 2 ] = 7; /*H */
    ncf[ 32 * kd + 1 ] = 2; /*O */

    /*NXC3H7O2 */
    ncf[ 33 * kd + 3 ] = 3; /*C */
    ncf[ 33 * kd + 2 ] = 7; /*H */
    ncf[ 33 * kd + 1 ] = 2; /*O */

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
    kname.resize(34);
    kname[0] = "N2";
    kname[1] = "O";
    kname[2] = "H";
    kname[3] = "OH";
    kname[4] = "H2";
    kname[5] = "O2";
    kname[6] = "H2O";
    kname[7] = "H2O2";
    kname[8] = "HO2";
    kname[9] = "CH2GSG";
    kname[10] = "CH2O";
    kname[11] = "CH3";
    kname[12] = "CO";
    kname[13] = "CH3O";
    kname[14] = "C2H5";
    kname[15] = "CH4";
    kname[16] = "C2H4";
    kname[17] = "C2H6";
    kname[18] = "CO2";
    kname[19] = "HCO";
    kname[20] = "CH3O2";
    kname[21] = "CH3O2H";
    kname[22] = "C2H2";
    kname[23] = "HCCO";
    kname[24] = "C2H3";
    kname[25] = "CH2CHO";
    kname[26] = "C3H5XA";
    kname[27] = "C3H6";
    kname[28] = "C3H5O";
    kname[29] = "IXC3H7";
    kname[30] = "NXC3H7";
    kname[31] = "C3H8";
    kname[32] = "IXC3H7O2";
    kname[33] = "NXC3H7O2";
}

/*compute the sparsity pattern of the chemistry Jacobian */
void SPARSITY_INFO( int * nJdata, const int * consP, int NCELLS)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(1225);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(34);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[1225];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int l=0; l<34; l++) {
                c_d[l] = 1.0/ 34.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<35; k++) {
        for (int l=0; l<35; l++) {
            if(J_h[ 35 * k + l] != 0.0){
                nJdata_tmp = nJdata_tmp + 1;
            }
        }
    }

    *nJdata = NCELLS * nJdata_tmp;
}



/*compute the sparsity pattern of the system Jacobian */
void SPARSITY_INFO_SYST( int * nJdata, const int * consP, int NCELLS)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(1225);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(34);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[1225];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<34; k++) {
                c_d[k] = 1.0/ 34.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<35; k++) {
        for (int l=0; l<35; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J_h[ 35 * k + l] != 0.0){
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(1225);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(34);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[1225];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<34; k++) {
                c_d[k] = 1.0/ 34.000000 ;
            }
            aJacobian_precond(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<35; k++) {
        for (int l=0; l<35; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J_h[ 35 * k + l] != 0.0){
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

    amrex::Gpu::DeviceVector<amrex::Real> J_v(1225);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(34);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[1225];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<34; k++) {
                c_d[k] = 1.0/ 34.000000 ;
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
        offset_row = nc * 35;
        offset_col = nc * 35;
        for (int k=0; k<35; k++) {
            for (int l=0; l<35; l++) {
                if(J_h[35*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(1225);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(34);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[1225];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<34; k++) {
                c_d[k] = 1.0/ 34.000000 ;
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
            offset = nc * 35;
            for (int l=0; l<35; l++) {
                for (int k=0; k<35; k++) {
                    if(J_h[35*k + l] != 0.0) {
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
            offset = nc * 35;
            for (int l=0; l<35; l++) {
                for (int k=0; k<35; k++) {
                    if(J_h[35*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(1225);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(34);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[1225];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<34; k++) {
                c_d[k] = 1.0/ 34.000000 ;
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
            offset = nc * 35;
            for (int l=0; l<35; l++) {
                for (int k=0; k<35; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp-1] = l+1 + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(J_h[35*k + l] != 0.0) {
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
            offset = nc * 35;
            for (int l=0; l<35; l++) {
                for (int k=0; k<35; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp] = l + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(J_h[35*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(1225);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(34);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[1225];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<34; k++) {
                c_d[k] = 1.0/ 34.000000 ;
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
    for (int k=0; k<35; k++) {
        for (int l=0; l<35; l++) {
            if (k == l) {
                rowVals[nJdata_tmp] = l; 
                indx[nJdata_tmp] = 35*k + l;
                nJdata_tmp = nJdata_tmp + 1; 
            } else {
                if(J_h[35*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l; 
                    indx[nJdata_tmp] = 35*k + l;
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(1225);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(34);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[1225];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<34; k++) {
                c_d[k] = 1.0/ 34.000000 ;
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
        for (int l=0; l<35; l++) {
            for (int k=0; k<35; k++) {
                if (k == l) {
                    colVals[nJdata_tmp-1] = l+1; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J_h[35*k + l] != 0.0) {
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
        for (int l=0; l<35; l++) {
            for (int k=0; k<35; k++) {
                if (k == l) {
                    colVals[nJdata_tmp] = l; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J_h[35*k + l] != 0.0) {
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
