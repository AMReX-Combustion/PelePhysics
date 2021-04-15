#ifndef MECHANISM_CPP
#define MECHANISM_CPP

#include "chemistry_file.H"



/*save atomic weights into array */
void atomicWeight(amrex::Real *  awt)
{
    awt[0] = 1.007970; /*H */
    awt[1] = 15.999400; /*O */
    awt[2] = 14.006700; /*N */
    awt[3] = 12.011150; /*C */

    return;
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
    for (id = 0; id < kd * 16; ++ id) {
         ncf[id] = 0; 
    }

    /*H */
    ncf[ 0 * kd + 0 ] = 1; /*H */

    /*O2 */
    ncf[ 1 * kd + 1 ] = 2; /*O */

    /*O */
    ncf[ 2 * kd + 1 ] = 1; /*O */

    /*H2O */
    ncf[ 3 * kd + 0 ] = 2; /*H */
    ncf[ 3 * kd + 1 ] = 1; /*O */

    /*OH */
    ncf[ 4 * kd + 1 ] = 1; /*O */
    ncf[ 4 * kd + 0 ] = 1; /*H */

    /*H2O2 */
    ncf[ 5 * kd + 0 ] = 2; /*H */
    ncf[ 5 * kd + 1 ] = 2; /*O */

    /*HO2 */
    ncf[ 6 * kd + 0 ] = 1; /*H */
    ncf[ 6 * kd + 1 ] = 2; /*O */

    /*CO */
    ncf[ 7 * kd + 3 ] = 1; /*C */
    ncf[ 7 * kd + 1 ] = 1; /*O */

    /*CO2 */
    ncf[ 8 * kd + 3 ] = 1; /*C */
    ncf[ 8 * kd + 1 ] = 2; /*O */

    /*CH4 */
    ncf[ 9 * kd + 3 ] = 1; /*C */
    ncf[ 9 * kd + 0 ] = 4; /*H */

    /*CH3 */
    ncf[ 10 * kd + 3 ] = 1; /*C */
    ncf[ 10 * kd + 0 ] = 3; /*H */

    /*CH3O2H */
    ncf[ 11 * kd + 3 ] = 1; /*C */
    ncf[ 11 * kd + 0 ] = 4; /*H */
    ncf[ 11 * kd + 1 ] = 2; /*O */

    /*CH3O2 */
    ncf[ 12 * kd + 0 ] = 3; /*H */
    ncf[ 12 * kd + 3 ] = 1; /*C */
    ncf[ 12 * kd + 1 ] = 2; /*O */

    /*CH3O */
    ncf[ 13 * kd + 3 ] = 1; /*C */
    ncf[ 13 * kd + 0 ] = 3; /*H */
    ncf[ 13 * kd + 1 ] = 1; /*O */

    /*CH2O */
    ncf[ 14 * kd + 0 ] = 2; /*H */
    ncf[ 14 * kd + 3 ] = 1; /*C */
    ncf[ 14 * kd + 1 ] = 1; /*O */

    /*HCO */
    ncf[ 15 * kd + 3 ] = 1; /*C */
    ncf[ 15 * kd + 0 ] = 1; /*H */
    ncf[ 15 * kd + 1 ] = 1; /*O */

}


/* Returns the vector of strings of element names */
void CKSYME_STR(amrex::Vector<std::string>& ename)
{
    ename.resize(4);
    ename[0] = "H";
    ename[1] = "O";
    ename[2] = "N";
    ename[3] = "C";
}


/* Returns the vector of strings of species names */
void CKSYMS_STR(amrex::Vector<std::string>& kname)
{
    kname.resize(16);
    kname[0] = "H";
    kname[1] = "O2";
    kname[2] = "O";
    kname[3] = "H2O";
    kname[4] = "OH";
    kname[5] = "H2O2";
    kname[6] = "HO2";
    kname[7] = "CO";
    kname[8] = "CO2";
    kname[9] = "CH4";
    kname[10] = "CH3";
    kname[11] = "CH3O2H";
    kname[12] = "CH3O2";
    kname[13] = "CH3O";
    kname[14] = "CH2O";
    kname[15] = "HCO";
}

/*compute the sparsity pattern of the chemistry Jacobian */
void SPARSITY_INFO( int * nJdata, int * consP, int NCELLS)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(289);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(16);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[289];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int l=0; l<16; l++) {
                c_d[l] = 1.0/ 16.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<17; k++) {
        for (int l=0; l<17; l++) {
            if(J_h[ 17 * k + l] != 0.0){
                nJdata_tmp = nJdata_tmp + 1;
            }
        }
    }

    *nJdata = NCELLS * nJdata_tmp;

    return;
}



/*compute the sparsity pattern of the system Jacobian */
void SPARSITY_INFO_SYST( int * nJdata, int * consP, int NCELLS)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(289);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(16);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[289];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<16; k++) {
                c_d[k] = 1.0/ 16.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<17; k++) {
        for (int l=0; l<17; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J_h[ 17 * k + l] != 0.0){
                    nJdata_tmp = nJdata_tmp + 1;
                }
            }
        }
    }

    *nJdata = NCELLS * nJdata_tmp;

    return;
}



/*compute the sparsity pattern of the simplified (for preconditioning) system Jacobian */
void SPARSITY_INFO_SYST_SIMPLIFIED( int * nJdata, int * consP)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(289);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(16);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[289];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<16; k++) {
                c_d[k] = 1.0/ 16.000000 ;
            }
            aJacobian_precond(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<17; k++) {
        for (int l=0; l<17; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J_h[ 17 * k + l] != 0.0){
                    nJdata_tmp = nJdata_tmp + 1;
                }
            }
        }
    }

    nJdata[0] = nJdata_tmp;

    return;
}


/*compute the sparsity pattern of the chemistry Jacobian in CSC format -- base 0 */
void SPARSITY_PREPROC_CSC(int *  rowVals, int *  colPtrs, int * consP, int NCELLS)
{
    int offset_row;
    int offset_col;

    amrex::Gpu::DeviceVector<amrex::Real> J_v(289);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(16);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[289];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<16; k++) {
                c_d[k] = 1.0/ 16.000000 ;
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
        offset_row = nc * 17;
        offset_col = nc * 17;
        for (int k=0; k<17; k++) {
            for (int l=0; l<17; l++) {
                if(J_h[17*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l + offset_row; 
                    nJdata_tmp = nJdata_tmp + 1; 
                }
            }
            colPtrs[offset_col + (k + 1)] = nJdata_tmp;
        }
    }

    return;
}

/*compute the sparsity pattern of the chemistry Jacobian in CSR format -- base 0 */
void SPARSITY_PREPROC_CSR(int * colVals, int * rowPtrs, int * consP, int NCELLS, int base)
{
    int offset;
    amrex::Gpu::DeviceVector<amrex::Real> J_v(289);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(16);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[289];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<16; k++) {
                c_d[k] = 1.0/ 16.000000 ;
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
            offset = nc * 17;
            for (int l=0; l<17; l++) {
                for (int k=0; k<17; k++) {
                    if(J_h[17*k + l] != 0.0) {
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
            offset = nc * 17;
            for (int l=0; l<17; l++) {
                for (int k=0; k<17; k++) {
                    if(J_h[17*k + l] != 0.0) {
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
void SPARSITY_PREPROC_SYST_CSR(int * colVals, int * rowPtr, int * consP, int NCELLS, int base)
{
    int offset;
    amrex::Gpu::DeviceVector<amrex::Real> J_v(289);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(16);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[289];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<16; k++) {
                c_d[k] = 1.0/ 16.000000 ;
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
            offset = nc * 17;
            for (int l=0; l<17; l++) {
                for (int k=0; k<17; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp-1] = l+1 + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(J_h[17*k + l] != 0.0) {
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
            offset = nc * 17;
            for (int l=0; l<17; l++) {
                for (int k=0; k<17; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp] = l + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(J_h[17*k + l] != 0.0) {
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

/*compute the sparsity pattern of the simplified (for precond) system Jacobian on CPU */
/*BASE 0 */
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(int * rowVals, int * colPtrs, int * indx, int * consP)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(289);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(16);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[289];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<16; k++) {
                c_d[k] = 1.0/ 16.000000 ;
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
    for (int k=0; k<17; k++) {
        for (int l=0; l<17; l++) {
            if (k == l) {
                rowVals[nJdata_tmp] = l; 
                indx[nJdata_tmp] = 17*k + l;
                nJdata_tmp = nJdata_tmp + 1; 
            } else {
                if(J_h[17*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l; 
                    indx[nJdata_tmp] = 17*k + l;
                    nJdata_tmp = nJdata_tmp + 1; 
                }
            }
        }
        colPtrs[k+1] = nJdata_tmp;
    }

    return;
}

/*compute the sparsity pattern of the simplified (for precond) system Jacobian */
/*CSR format BASE is under choice */
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(int * colVals, int * rowPtr, int * consP, int base)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(289);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(16);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[289];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<16; k++) {
                c_d[k] = 1.0/ 16.000000 ;
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
        for (int l=0; l<17; l++) {
            for (int k=0; k<17; k++) {
                if (k == l) {
                    colVals[nJdata_tmp-1] = l+1; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J_h[17*k + l] != 0.0) {
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
        for (int l=0; l<17; l++) {
            for (int k=0; k<17; k++) {
                if (k == l) {
                    colVals[nJdata_tmp] = l; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J_h[17*k + l] != 0.0) {
                        colVals[nJdata_tmp] = k; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    }
                }
            }
            rowPtr[l+1] = nJdata_tmp;
        }
    }

    return;
}

#endif

