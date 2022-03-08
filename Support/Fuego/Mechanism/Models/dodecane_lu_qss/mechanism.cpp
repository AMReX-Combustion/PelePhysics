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
    for (id = 0; id < kd * 35; ++ id) {
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

    /*CH3 */
    ncf[ 9 * kd + 2 ] = 1; /*C */
    ncf[ 9 * kd + 1 ] = 3; /*H */

    /*CH4 */
    ncf[ 10 * kd + 2 ] = 1; /*C */
    ncf[ 10 * kd + 1 ] = 4; /*H */

    /*CH2O */
    ncf[ 11 * kd + 1 ] = 2; /*H */
    ncf[ 11 * kd + 2 ] = 1; /*C */
    ncf[ 11 * kd + 0 ] = 1; /*O */

    /*CO */
    ncf[ 12 * kd + 2 ] = 1; /*C */
    ncf[ 12 * kd + 0 ] = 1; /*O */

    /*CO2 */
    ncf[ 13 * kd + 2 ] = 1; /*C */
    ncf[ 13 * kd + 0 ] = 2; /*O */

    /*C2H2 */
    ncf[ 14 * kd + 2 ] = 2; /*C */
    ncf[ 14 * kd + 1 ] = 2; /*H */

    /*C2H4 */
    ncf[ 15 * kd + 2 ] = 2; /*C */
    ncf[ 15 * kd + 1 ] = 4; /*H */

    /*C2H6 */
    ncf[ 16 * kd + 2 ] = 2; /*C */
    ncf[ 16 * kd + 1 ] = 6; /*H */

    /*CH2CHO */
    ncf[ 17 * kd + 0 ] = 1; /*O */
    ncf[ 17 * kd + 1 ] = 3; /*H */
    ncf[ 17 * kd + 2 ] = 2; /*C */

    /*aC3H5 */
    ncf[ 18 * kd + 2 ] = 3; /*C */
    ncf[ 18 * kd + 1 ] = 5; /*H */

    /*C3H6 */
    ncf[ 19 * kd + 2 ] = 3; /*C */
    ncf[ 19 * kd + 1 ] = 6; /*H */

    /*C2H3CHO */
    ncf[ 20 * kd + 2 ] = 3; /*C */
    ncf[ 20 * kd + 1 ] = 4; /*H */
    ncf[ 20 * kd + 0 ] = 1; /*O */

    /*C4H7 */
    ncf[ 21 * kd + 2 ] = 4; /*C */
    ncf[ 21 * kd + 1 ] = 7; /*H */

    /*C4H81 */
    ncf[ 22 * kd + 2 ] = 4; /*C */
    ncf[ 22 * kd + 1 ] = 8; /*H */

    /*C5H9 */
    ncf[ 23 * kd + 1 ] = 9; /*H */
    ncf[ 23 * kd + 2 ] = 5; /*C */

    /*C5H10 */
    ncf[ 24 * kd + 2 ] = 5; /*C */
    ncf[ 24 * kd + 1 ] = 10; /*H */

    /*C6H12 */
    ncf[ 25 * kd + 2 ] = 6; /*C */
    ncf[ 25 * kd + 1 ] = 12; /*H */

    /*C7H14 */
    ncf[ 26 * kd + 2 ] = 7; /*C */
    ncf[ 26 * kd + 1 ] = 14; /*H */

    /*C8H16 */
    ncf[ 27 * kd + 2 ] = 8; /*C */
    ncf[ 27 * kd + 1 ] = 16; /*H */

    /*C9H18 */
    ncf[ 28 * kd + 2 ] = 9; /*C */
    ncf[ 28 * kd + 1 ] = 18; /*H */

    /*PXC9H19 */
    ncf[ 29 * kd + 2 ] = 9; /*C */
    ncf[ 29 * kd + 1 ] = 19; /*H */

    /*C10H20 */
    ncf[ 30 * kd + 2 ] = 10; /*C */
    ncf[ 30 * kd + 1 ] = 20; /*H */

    /*C12H24 */
    ncf[ 31 * kd + 2 ] = 12; /*C */
    ncf[ 31 * kd + 1 ] = 24; /*H */

    /*C12H25O2 */
    ncf[ 32 * kd + 1 ] = 25; /*H */
    ncf[ 32 * kd + 0 ] = 2; /*O */
    ncf[ 32 * kd + 2 ] = 12; /*C */

    /*OC12H23OOH */
    ncf[ 33 * kd + 1 ] = 24; /*H */
    ncf[ 33 * kd + 0 ] = 3; /*O */
    ncf[ 33 * kd + 2 ] = 12; /*C */

    /*N2 */
    ncf[ 34 * kd + 3 ] = 2; /*N */

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
    kname.resize(35);
    kname[0] = "NC12H26";
    kname[1] = "H";
    kname[2] = "O";
    kname[3] = "OH";
    kname[4] = "HO2";
    kname[5] = "H2";
    kname[6] = "H2O";
    kname[7] = "H2O2";
    kname[8] = "O2";
    kname[9] = "CH3";
    kname[10] = "CH4";
    kname[11] = "CH2O";
    kname[12] = "CO";
    kname[13] = "CO2";
    kname[14] = "C2H2";
    kname[15] = "C2H4";
    kname[16] = "C2H6";
    kname[17] = "CH2CHO";
    kname[18] = "aC3H5";
    kname[19] = "C3H6";
    kname[20] = "C2H3CHO";
    kname[21] = "C4H7";
    kname[22] = "C4H81";
    kname[23] = "C5H9";
    kname[24] = "C5H10";
    kname[25] = "C6H12";
    kname[26] = "C7H14";
    kname[27] = "C8H16";
    kname[28] = "C9H18";
    kname[29] = "PXC9H19";
    kname[30] = "C10H20";
    kname[31] = "C12H24";
    kname[32] = "C12H25O2";
    kname[33] = "OC12H23OOH";
    kname[34] = "N2";
}

/*compute the sparsity pattern of the chemistry Jacobian */
void SPARSITY_INFO( int * nJdata, const int * consP, int NCELLS)
{
    amrex::GpuArray<amrex::Real,1296> Jac = {0.0};
    amrex::GpuArray<amrex::Real,35> conc = {0.0};
    for (int n=0; n<35; n++) {
        conc[n] = 1.0/ 35.000000 ;
    }
    aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

    int nJdata_tmp = 0;
    for (int k=0; k<36; k++) {
        for (int l=0; l<36; l++) {
            if(Jac[ 36 * k + l] != 0.0){
                nJdata_tmp = nJdata_tmp + 1;
            }
        }
    }

    *nJdata = NCELLS * nJdata_tmp;
}



/*compute the sparsity pattern of the system Jacobian */
void SPARSITY_INFO_SYST( int * nJdata, const int * consP, int NCELLS)
{
    amrex::GpuArray<amrex::Real,1296> Jac = {0.0};
    amrex::GpuArray<amrex::Real,35> conc = {0.0};
    for (int n=0; n<35; n++) {
        conc[n] = 1.0/ 35.000000 ;
    }
    aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

    int nJdata_tmp = 0;
    for (int k=0; k<36; k++) {
        for (int l=0; l<36; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(Jac[ 36 * k + l] != 0.0){
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
    amrex::GpuArray<amrex::Real,1296> Jac = {0.0};
    amrex::GpuArray<amrex::Real,35> conc = {0.0};
    for (int n=0; n<35; n++) {
        conc[n] = 1.0/ 35.000000 ;
    }
    aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

    int nJdata_tmp = 0;
    for (int k=0; k<36; k++) {
        for (int l=0; l<36; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(Jac[ 36 * k + l] != 0.0){
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
    amrex::GpuArray<amrex::Real,1296> Jac = {0.0};
    amrex::GpuArray<amrex::Real,35> conc = {0.0};
    for (int n=0; n<35; n++) {
        conc[n] = 1.0/ 35.000000 ;
    }
    aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

    colPtrs[0] = 0;
    int nJdata_tmp = 0;
    for (int nc=0; nc<NCELLS; nc++) {
        int offset_row = nc * 36;
        int offset_col = nc * 36;
        for (int k=0; k<36; k++) {
            for (int l=0; l<36; l++) {
                if(Jac[36*k + l] != 0.0) {
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
    amrex::GpuArray<amrex::Real,1296> Jac = {0.0};
    amrex::GpuArray<amrex::Real,35> conc = {0.0};
    for (int n=0; n<35; n++) {
        conc[n] = 1.0/ 35.000000 ;
    }
    aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

    if (base == 1) {
        rowPtrs[0] = 1;
        int nJdata_tmp = 1;
        for (int nc=0; nc<NCELLS; nc++) {
            int offset = nc * 36;
            for (int l=0; l<36; l++) {
                for (int k=0; k<36; k++) {
                    if(Jac[36*k + l] != 0.0) {
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
            int offset = nc * 36;
            for (int l=0; l<36; l++) {
                for (int k=0; k<36; k++) {
                    if(Jac[36*k + l] != 0.0) {
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
    amrex::GpuArray<amrex::Real,1296> Jac = {0.0};
    amrex::GpuArray<amrex::Real,35> conc = {0.0};
    for (int n=0; n<35; n++) {
        conc[n] = 1.0/ 35.000000 ;
    }
    aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

    if (base == 1) {
        rowPtr[0] = 1;
        int nJdata_tmp = 1;
        for (int nc=0; nc<NCELLS; nc++) {
            int offset = nc * 36;
            for (int l=0; l<36; l++) {
                for (int k=0; k<36; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp-1] = l+1 + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(Jac[36*k + l] != 0.0) {
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
            int offset = nc * 36;
            for (int l=0; l<36; l++) {
                for (int k=0; k<36; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp] = l + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(Jac[36*k + l] != 0.0) {
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
    amrex::GpuArray<amrex::Real,1296> Jac = {0.0};
    amrex::GpuArray<amrex::Real,35> conc = {0.0};
    for (int n=0; n<35; n++) {
        conc[n] = 1.0/ 35.000000 ;
    }
    aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

    colPtrs[0] = 0;
    int nJdata_tmp = 0;
    for (int k=0; k<36; k++) {
        for (int l=0; l<36; l++) {
            if (k == l) {
                rowVals[nJdata_tmp] = l; 
                indx[nJdata_tmp] = 36*k + l;
                nJdata_tmp = nJdata_tmp + 1; 
            } else {
                if(Jac[36*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l; 
                    indx[nJdata_tmp] = 36*k + l;
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
    amrex::GpuArray<amrex::Real,1296> Jac = {0.0};
    amrex::GpuArray<amrex::Real,35> conc = {0.0};
    for (int n=0; n<35; n++) {
        conc[n] = 1.0/ 35.000000 ;
    }
    aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

    if (base == 1) {
        rowPtr[0] = 1;
        int nJdata_tmp = 1;
        for (int l=0; l<36; l++) {
            for (int k=0; k<36; k++) {
                if (k == l) {
                    colVals[nJdata_tmp-1] = l+1; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(Jac[36*k + l] != 0.0) {
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
        for (int l=0; l<36; l++) {
            for (int k=0; k<36; k++) {
                if (k == l) {
                    colVals[nJdata_tmp] = l; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(Jac[36*k + l] != 0.0) {
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

/* End of file  */
