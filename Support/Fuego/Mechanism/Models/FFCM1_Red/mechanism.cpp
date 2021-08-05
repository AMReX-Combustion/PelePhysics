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
    awt[4] = 4.002600; /*HE */

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
    int kd = 5; 
    /*Zero ncf */
    for (id = 0; id < kd * 21; ++ id) {
         ncf[id] = 0; 
    }

    /*N2 */
    ncf[ 0 * kd + 3 ] = 2; /*N */

    /*H2 */
    ncf[ 1 * kd + 1 ] = 2; /*H */

    /*H */
    ncf[ 2 * kd + 1 ] = 1; /*H */

    /*O */
    ncf[ 3 * kd + 0 ] = 1; /*O */

    /*O2 */
    ncf[ 4 * kd + 0 ] = 2; /*O */

    /*OH */
    ncf[ 5 * kd + 0 ] = 1; /*O */
    ncf[ 5 * kd + 1 ] = 1; /*H */

    /*H2O */
    ncf[ 6 * kd + 1 ] = 2; /*H */
    ncf[ 6 * kd + 0 ] = 1; /*O */

    /*HO2 */
    ncf[ 7 * kd + 1 ] = 1; /*H */
    ncf[ 7 * kd + 0 ] = 2; /*O */

    /*CO */
    ncf[ 8 * kd + 2 ] = 1; /*C */
    ncf[ 8 * kd + 0 ] = 1; /*O */

    /*CO2 */
    ncf[ 9 * kd + 2 ] = 1; /*C */
    ncf[ 9 * kd + 0 ] = 2; /*O */

    /*CH */
    ncf[ 10 * kd + 2 ] = 1; /*C */
    ncf[ 10 * kd + 1 ] = 1; /*H */

    /*CH2 */
    ncf[ 11 * kd + 2 ] = 1; /*C */
    ncf[ 11 * kd + 1 ] = 2; /*H */

    /*CH2(S) */
    ncf[ 12 * kd + 2 ] = 1; /*C */
    ncf[ 12 * kd + 1 ] = 2; /*H */

    /*CH3 */
    ncf[ 13 * kd + 2 ] = 1; /*C */
    ncf[ 13 * kd + 1 ] = 3; /*H */

    /*CH4 */
    ncf[ 14 * kd + 2 ] = 1; /*C */
    ncf[ 14 * kd + 1 ] = 4; /*H */

    /*HCO */
    ncf[ 15 * kd + 2 ] = 1; /*C */
    ncf[ 15 * kd + 1 ] = 1; /*H */
    ncf[ 15 * kd + 0 ] = 1; /*O */

    /*CH2O */
    ncf[ 16 * kd + 1 ] = 2; /*H */
    ncf[ 16 * kd + 2 ] = 1; /*C */
    ncf[ 16 * kd + 0 ] = 1; /*O */

    /*CH2OH */
    ncf[ 17 * kd + 2 ] = 1; /*C */
    ncf[ 17 * kd + 1 ] = 3; /*H */
    ncf[ 17 * kd + 0 ] = 1; /*O */

    /*CH3O */
    ncf[ 18 * kd + 2 ] = 1; /*C */
    ncf[ 18 * kd + 1 ] = 3; /*H */
    ncf[ 18 * kd + 0 ] = 1; /*O */

    /*CH3OH */
    ncf[ 19 * kd + 2 ] = 1; /*C */
    ncf[ 19 * kd + 1 ] = 4; /*H */
    ncf[ 19 * kd + 0 ] = 1; /*O */

    /*HE */
    ncf[ 20 * kd + 4 ] = 1; /*HE */

}


/* Returns the vector of strings of element names */
void CKSYME_STR(amrex::Vector<std::string>& ename)
{
    ename.resize(5);
    ename[0] = "O";
    ename[1] = "H";
    ename[2] = "C";
    ename[3] = "N";
    ename[4] = "HE";
}


/* Returns the vector of strings of species names */
void CKSYMS_STR(amrex::Vector<std::string>& kname)
{
    kname.resize(21);
    kname[0] = "N2";
    kname[1] = "H2";
    kname[2] = "H";
    kname[3] = "O";
    kname[4] = "O2";
    kname[5] = "OH";
    kname[6] = "H2O";
    kname[7] = "HO2";
    kname[8] = "CO";
    kname[9] = "CO2";
    kname[10] = "CH";
    kname[11] = "CH2";
    kname[12] = "CH2(S)";
    kname[13] = "CH3";
    kname[14] = "CH4";
    kname[15] = "HCO";
    kname[16] = "CH2O";
    kname[17] = "CH2OH";
    kname[18] = "CH3O";
    kname[19] = "CH3OH";
    kname[20] = "HE";
}

/*compute the sparsity pattern of the chemistry Jacobian */
#ifdef COMPILE_JACOBIAN
void SPARSITY_INFO( int * nJdata, int * consP, int NCELLS)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(484);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(21);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[484];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int l=0; l<21; l++) {
                c_d[l] = 1.0/ 21.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<22; k++) {
        for (int l=0; l<22; l++) {
            if(J_h[ 22 * k + l] != 0.0){
                nJdata_tmp = nJdata_tmp + 1;
            }
        }
    }

    *nJdata = NCELLS * nJdata_tmp;

    return;
}
#endif



/*compute the sparsity pattern of the system Jacobian */
#ifdef COMPILE_JACOBIAN
void SPARSITY_INFO_SYST( int * nJdata, int * consP, int NCELLS)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(484);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(21);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[484];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<21; k++) {
                c_d[k] = 1.0/ 21.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<22; k++) {
        for (int l=0; l<22; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J_h[ 22 * k + l] != 0.0){
                    nJdata_tmp = nJdata_tmp + 1;
                }
            }
        }
    }

    *nJdata = NCELLS * nJdata_tmp;

    return;
}
#endif



/*compute the sparsity pattern of the simplified (for preconditioning) system Jacobian */
#ifdef COMPILE_JACOBIAN
void SPARSITY_INFO_SYST_SIMPLIFIED( int * nJdata, int * consP)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(484);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(21);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[484];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<21; k++) {
                c_d[k] = 1.0/ 21.000000 ;
            }
            aJacobian_precond(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<22; k++) {
        for (int l=0; l<22; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J_h[ 22 * k + l] != 0.0){
                    nJdata_tmp = nJdata_tmp + 1;
                }
            }
        }
    }

    nJdata[0] = nJdata_tmp;

    return;
}
#endif


/*compute the sparsity pattern of the chemistry Jacobian in CSC format -- base 0 */
#ifdef COMPILE_JACOBIAN
void SPARSITY_PREPROC_CSC(int *  rowVals, int *  colPtrs, int * consP, int NCELLS)
{
    int offset_row;
    int offset_col;

    amrex::Gpu::DeviceVector<amrex::Real> J_v(484);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(21);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[484];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<21; k++) {
                c_d[k] = 1.0/ 21.000000 ;
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
        offset_row = nc * 22;
        offset_col = nc * 22;
        for (int k=0; k<22; k++) {
            for (int l=0; l<22; l++) {
                if(J_h[22*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l + offset_row; 
                    nJdata_tmp = nJdata_tmp + 1; 
                }
            }
            colPtrs[offset_col + (k + 1)] = nJdata_tmp;
        }
    }

    return;
}
#endif

/*compute the sparsity pattern of the chemistry Jacobian in CSR format -- base 0 */
#ifdef COMPILE_JACOBIAN
void SPARSITY_PREPROC_CSR(int * colVals, int * rowPtrs, int * consP, int NCELLS, int base)
{
    int offset;
    amrex::Gpu::DeviceVector<amrex::Real> J_v(484);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(21);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[484];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<21; k++) {
                c_d[k] = 1.0/ 21.000000 ;
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
            offset = nc * 22;
            for (int l=0; l<22; l++) {
                for (int k=0; k<22; k++) {
                    if(J_h[22*k + l] != 0.0) {
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
            offset = nc * 22;
            for (int l=0; l<22; l++) {
                for (int k=0; k<22; k++) {
                    if(J_h[22*k + l] != 0.0) {
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
#endif

/*compute the sparsity pattern of the system Jacobian */
/*CSR format BASE is user choice */
#ifdef COMPILE_JACOBIAN
void SPARSITY_PREPROC_SYST_CSR(int * colVals, int * rowPtr, int * consP, int NCELLS, int base)
{
    int offset;
    amrex::Gpu::DeviceVector<amrex::Real> J_v(484);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(21);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[484];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<21; k++) {
                c_d[k] = 1.0/ 21.000000 ;
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
            offset = nc * 22;
            for (int l=0; l<22; l++) {
                for (int k=0; k<22; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp-1] = l+1 + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(J_h[22*k + l] != 0.0) {
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
            offset = nc * 22;
            for (int l=0; l<22; l++) {
                for (int k=0; k<22; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp] = l + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(J_h[22*k + l] != 0.0) {
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
#endif

/*compute the sparsity pattern of the simplified (for precond) system Jacobian on CPU */
/*BASE 0 */
#ifdef COMPILE_JACOBIAN
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(int * rowVals, int * colPtrs, int * indx, int * consP)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(484);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(21);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[484];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<21; k++) {
                c_d[k] = 1.0/ 21.000000 ;
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
    for (int k=0; k<22; k++) {
        for (int l=0; l<22; l++) {
            if (k == l) {
                rowVals[nJdata_tmp] = l; 
                indx[nJdata_tmp] = 22*k + l;
                nJdata_tmp = nJdata_tmp + 1; 
            } else {
                if(J_h[22*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l; 
                    indx[nJdata_tmp] = 22*k + l;
                    nJdata_tmp = nJdata_tmp + 1; 
                }
            }
        }
        colPtrs[k+1] = nJdata_tmp;
    }

    return;
}
#endif

/*compute the sparsity pattern of the simplified (for precond) system Jacobian */
/*CSR format BASE is under choice */
#ifdef COMPILE_JACOBIAN
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(int * colVals, int * rowPtr, int * consP, int base)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(484);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(21);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[484];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<21; k++) {
                c_d[k] = 1.0/ 21.000000 ;
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
        for (int l=0; l<22; l++) {
            for (int k=0; k<22; k++) {
                if (k == l) {
                    colVals[nJdata_tmp-1] = l+1; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J_h[22*k + l] != 0.0) {
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
        for (int l=0; l<22; l++) {
            for (int k=0; k<22; k++) {
                if (k == l) {
                    colVals[nJdata_tmp] = l; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J_h[22*k + l] != 0.0) {
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

#endif
