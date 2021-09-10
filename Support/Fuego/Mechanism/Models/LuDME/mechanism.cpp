#ifndef MECHANISM_CPP
#define MECHANISM_CPP

#include "mechanism.H"



/*save atomic weights into array */
void atomicWeight(amrex::Real *  awt)
{
    awt[0] = 12.011150; /*C */
    awt[1] = 1.007970; /*H */
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
    for (id = 0; id < kd * 39; ++ id) {
         ncf[id] = 0; 
    }

    /*H */
    ncf[ 0 * kd + 1 ] = 1; /*H */

    /*H2 */
    ncf[ 1 * kd + 1 ] = 2; /*H */

    /*CH2 */
    ncf[ 2 * kd + 0 ] = 1; /*C */
    ncf[ 2 * kd + 1 ] = 2; /*H */

    /*CH2(S) */
    ncf[ 3 * kd + 0 ] = 1; /*C */
    ncf[ 3 * kd + 1 ] = 2; /*H */

    /*CH3 */
    ncf[ 4 * kd + 0 ] = 1; /*C */
    ncf[ 4 * kd + 1 ] = 3; /*H */

    /*O */
    ncf[ 5 * kd + 2 ] = 1; /*O */

    /*CH4 */
    ncf[ 6 * kd + 0 ] = 1; /*C */
    ncf[ 6 * kd + 1 ] = 4; /*H */

    /*OH */
    ncf[ 7 * kd + 2 ] = 1; /*O */
    ncf[ 7 * kd + 1 ] = 1; /*H */

    /*H2O */
    ncf[ 8 * kd + 1 ] = 2; /*H */
    ncf[ 8 * kd + 2 ] = 1; /*O */

    /*C2H2 */
    ncf[ 9 * kd + 0 ] = 2; /*C */
    ncf[ 9 * kd + 1 ] = 2; /*H */

    /*C2H3 */
    ncf[ 10 * kd + 0 ] = 2; /*C */
    ncf[ 10 * kd + 1 ] = 3; /*H */

    /*CO */
    ncf[ 11 * kd + 0 ] = 1; /*C */
    ncf[ 11 * kd + 2 ] = 1; /*O */

    /*C2H4 */
    ncf[ 12 * kd + 0 ] = 2; /*C */
    ncf[ 12 * kd + 1 ] = 4; /*H */

    /*HCO */
    ncf[ 13 * kd + 1 ] = 1; /*H */
    ncf[ 13 * kd + 0 ] = 1; /*C */
    ncf[ 13 * kd + 2 ] = 1; /*O */

    /*C2H5 */
    ncf[ 14 * kd + 0 ] = 2; /*C */
    ncf[ 14 * kd + 1 ] = 5; /*H */

    /*CH2O */
    ncf[ 15 * kd + 0 ] = 1; /*C */
    ncf[ 15 * kd + 1 ] = 2; /*H */
    ncf[ 15 * kd + 2 ] = 1; /*O */

    /*C2H6 */
    ncf[ 16 * kd + 0 ] = 2; /*C */
    ncf[ 16 * kd + 1 ] = 6; /*H */

    /*CH2OH */
    ncf[ 17 * kd + 0 ] = 1; /*C */
    ncf[ 17 * kd + 1 ] = 3; /*H */
    ncf[ 17 * kd + 2 ] = 1; /*O */

    /*CH3O */
    ncf[ 18 * kd + 0 ] = 1; /*C */
    ncf[ 18 * kd + 1 ] = 3; /*H */
    ncf[ 18 * kd + 2 ] = 1; /*O */

    /*O2 */
    ncf[ 19 * kd + 2 ] = 2; /*O */

    /*HO2 */
    ncf[ 20 * kd + 1 ] = 1; /*H */
    ncf[ 20 * kd + 2 ] = 2; /*O */

    /*H2O2 */
    ncf[ 21 * kd + 1 ] = 2; /*H */
    ncf[ 21 * kd + 2 ] = 2; /*O */

    /*CO2 */
    ncf[ 22 * kd + 0 ] = 1; /*C */
    ncf[ 22 * kd + 2 ] = 2; /*O */

    /*CH3HCO */
    ncf[ 23 * kd + 0 ] = 2; /*C */
    ncf[ 23 * kd + 1 ] = 4; /*H */
    ncf[ 23 * kd + 2 ] = 1; /*O */

    /*CH3OCH2 */
    ncf[ 24 * kd + 0 ] = 2; /*C */
    ncf[ 24 * kd + 1 ] = 5; /*H */
    ncf[ 24 * kd + 2 ] = 1; /*O */

    /*HCOOH */
    ncf[ 25 * kd + 0 ] = 1; /*C */
    ncf[ 25 * kd + 1 ] = 2; /*H */
    ncf[ 25 * kd + 2 ] = 2; /*O */

    /*CH3OCH3 */
    ncf[ 26 * kd + 0 ] = 2; /*C */
    ncf[ 26 * kd + 1 ] = 6; /*H */
    ncf[ 26 * kd + 2 ] = 1; /*O */

    /*HOCH2O */
    ncf[ 27 * kd + 0 ] = 1; /*C */
    ncf[ 27 * kd + 1 ] = 3; /*H */
    ncf[ 27 * kd + 2 ] = 2; /*O */

    /*CH3OCO */
    ncf[ 28 * kd + 0 ] = 2; /*C */
    ncf[ 28 * kd + 1 ] = 3; /*H */
    ncf[ 28 * kd + 2 ] = 2; /*O */

    /*CH3OCHO */
    ncf[ 29 * kd + 0 ] = 2; /*C */
    ncf[ 29 * kd + 1 ] = 4; /*H */
    ncf[ 29 * kd + 2 ] = 2; /*O */

    /*CH3OCH2O */
    ncf[ 30 * kd + 0 ] = 2; /*C */
    ncf[ 30 * kd + 1 ] = 5; /*H */
    ncf[ 30 * kd + 2 ] = 2; /*O */

    /*CH3OCH2OH */
    ncf[ 31 * kd + 0 ] = 2; /*C */
    ncf[ 31 * kd + 1 ] = 6; /*H */
    ncf[ 31 * kd + 2 ] = 2; /*O */

    /*OCH2OCHO */
    ncf[ 32 * kd + 0 ] = 2; /*C */
    ncf[ 32 * kd + 1 ] = 3; /*H */
    ncf[ 32 * kd + 2 ] = 3; /*O */

    /*HOCH2OCO */
    ncf[ 33 * kd + 0 ] = 2; /*C */
    ncf[ 33 * kd + 1 ] = 3; /*H */
    ncf[ 33 * kd + 2 ] = 3; /*O */

    /*CH3OCH2O2 */
    ncf[ 34 * kd + 0 ] = 2; /*C */
    ncf[ 34 * kd + 1 ] = 5; /*H */
    ncf[ 34 * kd + 2 ] = 3; /*O */

    /*CH2OCH2O2H */
    ncf[ 35 * kd + 0 ] = 2; /*C */
    ncf[ 35 * kd + 1 ] = 5; /*H */
    ncf[ 35 * kd + 2 ] = 3; /*O */

    /*HO2CH2OCHO */
    ncf[ 36 * kd + 0 ] = 2; /*C */
    ncf[ 36 * kd + 1 ] = 4; /*H */
    ncf[ 36 * kd + 2 ] = 4; /*O */

    /*O2CH2OCH2O2H */
    ncf[ 37 * kd + 0 ] = 2; /*C */
    ncf[ 37 * kd + 1 ] = 5; /*H */
    ncf[ 37 * kd + 2 ] = 5; /*O */

    /*N2 */
    ncf[ 38 * kd + 3 ] = 2; /*N */

}


/* Returns the vector of strings of element names */
void CKSYME_STR(amrex::Vector<std::string>& ename)
{
    ename.resize(4);
    ename[0] = "C";
    ename[1] = "H";
    ename[2] = "O";
    ename[3] = "N";
}


/* Returns the vector of strings of species names */
void CKSYMS_STR(amrex::Vector<std::string>& kname)
{
    kname.resize(39);
    kname[0] = "H";
    kname[1] = "H2";
    kname[2] = "CH2";
    kname[3] = "CH2(S)";
    kname[4] = "CH3";
    kname[5] = "O";
    kname[6] = "CH4";
    kname[7] = "OH";
    kname[8] = "H2O";
    kname[9] = "C2H2";
    kname[10] = "C2H3";
    kname[11] = "CO";
    kname[12] = "C2H4";
    kname[13] = "HCO";
    kname[14] = "C2H5";
    kname[15] = "CH2O";
    kname[16] = "C2H6";
    kname[17] = "CH2OH";
    kname[18] = "CH3O";
    kname[19] = "O2";
    kname[20] = "HO2";
    kname[21] = "H2O2";
    kname[22] = "CO2";
    kname[23] = "CH3HCO";
    kname[24] = "CH3OCH2";
    kname[25] = "HCOOH";
    kname[26] = "CH3OCH3";
    kname[27] = "HOCH2O";
    kname[28] = "CH3OCO";
    kname[29] = "CH3OCHO";
    kname[30] = "CH3OCH2O";
    kname[31] = "CH3OCH2OH";
    kname[32] = "OCH2OCHO";
    kname[33] = "HOCH2OCO";
    kname[34] = "CH3OCH2O2";
    kname[35] = "CH2OCH2O2H";
    kname[36] = "HO2CH2OCHO";
    kname[37] = "O2CH2OCH2O2H";
    kname[38] = "N2";
}

/*compute the sparsity pattern of the chemistry Jacobian */
void SPARSITY_INFO( int * nJdata, const int * consP, int NCELLS)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(1600);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(39);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[1600];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int l=0; l<39; l++) {
                c_d[l] = 1.0/ 39.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<40; k++) {
        for (int l=0; l<40; l++) {
            if(J_h[ 40 * k + l] != 0.0){
                nJdata_tmp = nJdata_tmp + 1;
            }
        }
    }

    *nJdata = NCELLS * nJdata_tmp;
}



/*compute the sparsity pattern of the system Jacobian */
void SPARSITY_INFO_SYST( int * nJdata, const int * consP, int NCELLS)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(1600);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(39);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[1600];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<39; k++) {
                c_d[k] = 1.0/ 39.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<40; k++) {
        for (int l=0; l<40; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J_h[ 40 * k + l] != 0.0){
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(1600);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(39);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[1600];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<39; k++) {
                c_d[k] = 1.0/ 39.000000 ;
            }
            aJacobian_precond(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<40; k++) {
        for (int l=0; l<40; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J_h[ 40 * k + l] != 0.0){
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

    amrex::Gpu::DeviceVector<amrex::Real> J_v(1600);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(39);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[1600];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<39; k++) {
                c_d[k] = 1.0/ 39.000000 ;
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
        offset_row = nc * 40;
        offset_col = nc * 40;
        for (int k=0; k<40; k++) {
            for (int l=0; l<40; l++) {
                if(J_h[40*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(1600);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(39);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[1600];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<39; k++) {
                c_d[k] = 1.0/ 39.000000 ;
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
            offset = nc * 40;
            for (int l=0; l<40; l++) {
                for (int k=0; k<40; k++) {
                    if(J_h[40*k + l] != 0.0) {
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
            offset = nc * 40;
            for (int l=0; l<40; l++) {
                for (int k=0; k<40; k++) {
                    if(J_h[40*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(1600);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(39);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[1600];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<39; k++) {
                c_d[k] = 1.0/ 39.000000 ;
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
            offset = nc * 40;
            for (int l=0; l<40; l++) {
                for (int k=0; k<40; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp-1] = l+1 + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(J_h[40*k + l] != 0.0) {
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
            offset = nc * 40;
            for (int l=0; l<40; l++) {
                for (int k=0; k<40; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp] = l + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(J_h[40*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(1600);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(39);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[1600];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<39; k++) {
                c_d[k] = 1.0/ 39.000000 ;
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
    for (int k=0; k<40; k++) {
        for (int l=0; l<40; l++) {
            if (k == l) {
                rowVals[nJdata_tmp] = l; 
                indx[nJdata_tmp] = 40*k + l;
                nJdata_tmp = nJdata_tmp + 1; 
            } else {
                if(J_h[40*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l; 
                    indx[nJdata_tmp] = 40*k + l;
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(1600);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(39);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[1600];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<39; k++) {
                c_d[k] = 1.0/ 39.000000 ;
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
        for (int l=0; l<40; l++) {
            for (int k=0; k<40; k++) {
                if (k == l) {
                    colVals[nJdata_tmp-1] = l+1; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J_h[40*k + l] != 0.0) {
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
        for (int l=0; l<40; l++) {
            for (int k=0; k<40; k++) {
                if (k == l) {
                    colVals[nJdata_tmp] = l; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J_h[40*k + l] != 0.0) {
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
