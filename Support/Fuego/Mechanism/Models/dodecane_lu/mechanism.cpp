#ifndef MECHANISM_CPP
#define MECHANISM_CPP

#include "chemistry_file.H"



/*save atomic weights into array */
void atomicWeight(amrex::Real *  awt)
{
    awt[0] = 15.999400; /*O */
    awt[1] = 1.007970; /*H */
    awt[2] = 12.011150; /*C */
    awt[3] = 14.006700; /*N */

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
    for (id = 0; id < kd * 53; ++ id) {
         ncf[id] = 0; 
    }

    /*NC12H26 */
    ncf[ 0 * kd + 2 ] = 12; /*C */
    ncf[ 0 * kd + 1 ] = 26; /*H */

    /*H */
    ncf[ 1 * kd + 1 ] = 1; /*H */

    /*O */
    ncf[ 2 * kd + 0 ] = 1; /*O */

    /*OH */
    ncf[ 3 * kd + 0 ] = 1; /*O */
    ncf[ 3 * kd + 1 ] = 1; /*H */

    /*HO2 */
    ncf[ 4 * kd + 1 ] = 1; /*H */
    ncf[ 4 * kd + 0 ] = 2; /*O */

    /*H2 */
    ncf[ 5 * kd + 1 ] = 2; /*H */

    /*H2O */
    ncf[ 6 * kd + 1 ] = 2; /*H */
    ncf[ 6 * kd + 0 ] = 1; /*O */

    /*H2O2 */
    ncf[ 7 * kd + 1 ] = 2; /*H */
    ncf[ 7 * kd + 0 ] = 2; /*O */

    /*O2 */
    ncf[ 8 * kd + 0 ] = 2; /*O */

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

    /*CO */
    ncf[ 16 * kd + 2 ] = 1; /*C */
    ncf[ 16 * kd + 0 ] = 1; /*O */

    /*CO2 */
    ncf[ 17 * kd + 2 ] = 1; /*C */
    ncf[ 17 * kd + 0 ] = 2; /*O */

    /*C2H2 */
    ncf[ 18 * kd + 2 ] = 2; /*C */
    ncf[ 18 * kd + 1 ] = 2; /*H */

    /*C2H3 */
    ncf[ 19 * kd + 2 ] = 2; /*C */
    ncf[ 19 * kd + 1 ] = 3; /*H */

    /*C2H4 */
    ncf[ 20 * kd + 2 ] = 2; /*C */
    ncf[ 20 * kd + 1 ] = 4; /*H */

    /*C2H5 */
    ncf[ 21 * kd + 2 ] = 2; /*C */
    ncf[ 21 * kd + 1 ] = 5; /*H */

    /*C2H6 */
    ncf[ 22 * kd + 2 ] = 2; /*C */
    ncf[ 22 * kd + 1 ] = 6; /*H */

    /*CH2CHO */
    ncf[ 23 * kd + 0 ] = 1; /*O */
    ncf[ 23 * kd + 1 ] = 3; /*H */
    ncf[ 23 * kd + 2 ] = 2; /*C */

    /*aC3H5 */
    ncf[ 24 * kd + 2 ] = 3; /*C */
    ncf[ 24 * kd + 1 ] = 5; /*H */

    /*C3H6 */
    ncf[ 25 * kd + 2 ] = 3; /*C */
    ncf[ 25 * kd + 1 ] = 6; /*H */

    /*nC3H7 */
    ncf[ 26 * kd + 2 ] = 3; /*C */
    ncf[ 26 * kd + 1 ] = 7; /*H */

    /*C2H3CHO */
    ncf[ 27 * kd + 2 ] = 3; /*C */
    ncf[ 27 * kd + 1 ] = 4; /*H */
    ncf[ 27 * kd + 0 ] = 1; /*O */

    /*C4H7 */
    ncf[ 28 * kd + 2 ] = 4; /*C */
    ncf[ 28 * kd + 1 ] = 7; /*H */

    /*C4H81 */
    ncf[ 29 * kd + 2 ] = 4; /*C */
    ncf[ 29 * kd + 1 ] = 8; /*H */

    /*pC4H9 */
    ncf[ 30 * kd + 2 ] = 4; /*C */
    ncf[ 30 * kd + 1 ] = 9; /*H */

    /*C5H9 */
    ncf[ 31 * kd + 1 ] = 9; /*H */
    ncf[ 31 * kd + 2 ] = 5; /*C */

    /*C5H10 */
    ncf[ 32 * kd + 2 ] = 5; /*C */
    ncf[ 32 * kd + 1 ] = 10; /*H */

    /*PXC5H11 */
    ncf[ 33 * kd + 2 ] = 5; /*C */
    ncf[ 33 * kd + 1 ] = 11; /*H */

    /*C6H12 */
    ncf[ 34 * kd + 2 ] = 6; /*C */
    ncf[ 34 * kd + 1 ] = 12; /*H */

    /*PXC6H13 */
    ncf[ 35 * kd + 2 ] = 6; /*C */
    ncf[ 35 * kd + 1 ] = 13; /*H */

    /*C7H14 */
    ncf[ 36 * kd + 2 ] = 7; /*C */
    ncf[ 36 * kd + 1 ] = 14; /*H */

    /*PXC7H15 */
    ncf[ 37 * kd + 2 ] = 7; /*C */
    ncf[ 37 * kd + 1 ] = 15; /*H */

    /*C8H16 */
    ncf[ 38 * kd + 2 ] = 8; /*C */
    ncf[ 38 * kd + 1 ] = 16; /*H */

    /*PXC8H17 */
    ncf[ 39 * kd + 2 ] = 8; /*C */
    ncf[ 39 * kd + 1 ] = 17; /*H */

    /*C9H18 */
    ncf[ 40 * kd + 2 ] = 9; /*C */
    ncf[ 40 * kd + 1 ] = 18; /*H */

    /*PXC9H19 */
    ncf[ 41 * kd + 2 ] = 9; /*C */
    ncf[ 41 * kd + 1 ] = 19; /*H */

    /*C10H20 */
    ncf[ 42 * kd + 2 ] = 10; /*C */
    ncf[ 42 * kd + 1 ] = 20; /*H */

    /*PXC10H21 */
    ncf[ 43 * kd + 2 ] = 10; /*C */
    ncf[ 43 * kd + 1 ] = 21; /*H */

    /*PXC12H25 */
    ncf[ 44 * kd + 2 ] = 12; /*C */
    ncf[ 44 * kd + 1 ] = 25; /*H */

    /*SXC12H25 */
    ncf[ 45 * kd + 2 ] = 12; /*C */
    ncf[ 45 * kd + 1 ] = 25; /*H */

    /*S3XC12H25 */
    ncf[ 46 * kd + 2 ] = 12; /*C */
    ncf[ 46 * kd + 1 ] = 25; /*H */

    /*C12H24 */
    ncf[ 47 * kd + 2 ] = 12; /*C */
    ncf[ 47 * kd + 1 ] = 24; /*H */

    /*C12H25O2 */
    ncf[ 48 * kd + 1 ] = 25; /*H */
    ncf[ 48 * kd + 0 ] = 2; /*O */
    ncf[ 48 * kd + 2 ] = 12; /*C */

    /*C12OOH */
    ncf[ 49 * kd + 1 ] = 25; /*H */
    ncf[ 49 * kd + 0 ] = 2; /*O */
    ncf[ 49 * kd + 2 ] = 12; /*C */

    /*O2C12H24OOH */
    ncf[ 50 * kd + 1 ] = 25; /*H */
    ncf[ 50 * kd + 0 ] = 4; /*O */
    ncf[ 50 * kd + 2 ] = 12; /*C */

    /*OC12H23OOH */
    ncf[ 51 * kd + 1 ] = 24; /*H */
    ncf[ 51 * kd + 0 ] = 3; /*O */
    ncf[ 51 * kd + 2 ] = 12; /*C */

    /*N2 */
    ncf[ 52 * kd + 3 ] = 2; /*N */

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
    kname.resize(53);
    kname[0] = "NC12H26";
    kname[1] = "H";
    kname[2] = "O";
    kname[3] = "OH";
    kname[4] = "HO2";
    kname[5] = "H2";
    kname[6] = "H2O";
    kname[7] = "H2O2";
    kname[8] = "O2";
    kname[9] = "CH2";
    kname[10] = "CH2*";
    kname[11] = "CH3";
    kname[12] = "CH4";
    kname[13] = "HCO";
    kname[14] = "CH2O";
    kname[15] = "CH3O";
    kname[16] = "CO";
    kname[17] = "CO2";
    kname[18] = "C2H2";
    kname[19] = "C2H3";
    kname[20] = "C2H4";
    kname[21] = "C2H5";
    kname[22] = "C2H6";
    kname[23] = "CH2CHO";
    kname[24] = "aC3H5";
    kname[25] = "C3H6";
    kname[26] = "nC3H7";
    kname[27] = "C2H3CHO";
    kname[28] = "C4H7";
    kname[29] = "C4H81";
    kname[30] = "pC4H9";
    kname[31] = "C5H9";
    kname[32] = "C5H10";
    kname[33] = "PXC5H11";
    kname[34] = "C6H12";
    kname[35] = "PXC6H13";
    kname[36] = "C7H14";
    kname[37] = "PXC7H15";
    kname[38] = "C8H16";
    kname[39] = "PXC8H17";
    kname[40] = "C9H18";
    kname[41] = "PXC9H19";
    kname[42] = "C10H20";
    kname[43] = "PXC10H21";
    kname[44] = "PXC12H25";
    kname[45] = "SXC12H25";
    kname[46] = "S3XC12H25";
    kname[47] = "C12H24";
    kname[48] = "C12H25O2";
    kname[49] = "C12OOH";
    kname[50] = "O2C12H24OOH";
    kname[51] = "OC12H23OOH";
    kname[52] = "N2";
}

/*compute the sparsity pattern of the chemistry Jacobian */
void SPARSITY_INFO( int * nJdata, int * consP, int NCELLS)
{
    amrex::Gpu::DeviceVector<amrex::Real> J_v(2916);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(53);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[2916];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int l=0; l<53; l++) {
                c_d[l] = 1.0/ 53.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<54; k++) {
        for (int l=0; l<54; l++) {
            if(J_h[ 54 * k + l] != 0.0){
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(2916);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(53);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[2916];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<53; k++) {
                c_d[k] = 1.0/ 53.000000 ;
            }
            aJacobian(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<54; k++) {
        for (int l=0; l<54; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J_h[ 54 * k + l] != 0.0){
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(2916);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(53);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[2916];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<53; k++) {
                c_d[k] = 1.0/ 53.000000 ;
            }
            aJacobian_precond(J_d, c_d, 1500.0, *consP);
    });

#ifdef AMREX_USE_GPU
    amrex::Gpu::dtoh_memcpy(J_h, J_d, sizeof(J_d));
#else
    std::memcpy(&J_h, J_d, sizeof(J_h));
#endif

    int nJdata_tmp = 0;
    for (int k=0; k<54; k++) {
        for (int l=0; l<54; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(J_h[ 54 * k + l] != 0.0){
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

    amrex::Gpu::DeviceVector<amrex::Real> J_v(2916);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(53);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[2916];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<53; k++) {
                c_d[k] = 1.0/ 53.000000 ;
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
        offset_row = nc * 54;
        offset_col = nc * 54;
        for (int k=0; k<54; k++) {
            for (int l=0; l<54; l++) {
                if(J_h[54*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(2916);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(53);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[2916];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<53; k++) {
                c_d[k] = 1.0/ 53.000000 ;
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
            offset = nc * 54;
            for (int l=0; l<54; l++) {
                for (int k=0; k<54; k++) {
                    if(J_h[54*k + l] != 0.0) {
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
            offset = nc * 54;
            for (int l=0; l<54; l++) {
                for (int k=0; k<54; k++) {
                    if(J_h[54*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(2916);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(53);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[2916];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<53; k++) {
                c_d[k] = 1.0/ 53.000000 ;
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
            offset = nc * 54;
            for (int l=0; l<54; l++) {
                for (int k=0; k<54; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp-1] = l+1 + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(J_h[54*k + l] != 0.0) {
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
            offset = nc * 54;
            for (int l=0; l<54; l++) {
                for (int k=0; k<54; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp] = l + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(J_h[54*k + l] != 0.0) {
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(2916);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(53);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[2916];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<53; k++) {
                c_d[k] = 1.0/ 53.000000 ;
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
    for (int k=0; k<54; k++) {
        for (int l=0; l<54; l++) {
            if (k == l) {
                rowVals[nJdata_tmp] = l; 
                indx[nJdata_tmp] = 54*k + l;
                nJdata_tmp = nJdata_tmp + 1; 
            } else {
                if(J_h[54*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l; 
                    indx[nJdata_tmp] = 54*k + l;
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
    amrex::Gpu::DeviceVector<amrex::Real> J_v(2916);
    amrex::Gpu::DeviceVector<amrex::Real> c_v(53);
    amrex::Real * J_d = J_v.data();
    amrex::Real * c_d = c_v.data();

    amrex::Real J_h[2916];

    amrex::IntVect iv(AMREX_D_DECL(0,0,0));
    amrex::ParallelFor(amrex::Box(iv,iv),
        [=] AMREX_GPU_HOST_DEVICE (int /*i*/, int /*j*/, int /*k*/) noexcept {
            for (int k=0; k<53; k++) {
                c_d[k] = 1.0/ 53.000000 ;
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
        for (int l=0; l<54; l++) {
            for (int k=0; k<54; k++) {
                if (k == l) {
                    colVals[nJdata_tmp-1] = l+1; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J_h[54*k + l] != 0.0) {
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
        for (int l=0; l<54; l++) {
            for (int k=0; k<54; k++) {
                if (k == l) {
                    colVals[nJdata_tmp] = l; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(J_h[54*k + l] != 0.0) {
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

