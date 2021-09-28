#ifndef MECHANISM_CPP
#define MECHANISM_CPP

#include "mechanism.H"



/*save atomic weights into array */
void atomicWeight(amrex::Real *  awt)
{
    awt[0] = 14.006700; /*N */
    awt[1] = 12.011150; /*C */
    awt[2] = 1.007970; /*H */
    awt[3] = 15.999400; /*O */
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
    for (id = 0; id < kd * 47; ++ id) {
         ncf[id] = 0; 
    }

    /*N2 */
    ncf[ 0 * kd + 0 ] = 2; /*N */

    /*S-CH2 */
    ncf[ 1 * kd + 1 ] = 1; /*C */
    ncf[ 1 * kd + 2 ] = 2; /*H */

    /*T-CH2 */
    ncf[ 2 * kd + 1 ] = 1; /*C */
    ncf[ 2 * kd + 2 ] = 2; /*H */

    /*O */
    ncf[ 3 * kd + 3 ] = 1; /*O */

    /*H2 */
    ncf[ 4 * kd + 2 ] = 2; /*H */

    /*H */
    ncf[ 5 * kd + 2 ] = 1; /*H */

    /*OH */
    ncf[ 6 * kd + 3 ] = 1; /*O */
    ncf[ 6 * kd + 2 ] = 1; /*H */

    /*H2O */
    ncf[ 7 * kd + 2 ] = 2; /*H */
    ncf[ 7 * kd + 3 ] = 1; /*O */

    /*O2 */
    ncf[ 8 * kd + 3 ] = 2; /*O */

    /*HO2 */
    ncf[ 9 * kd + 2 ] = 1; /*H */
    ncf[ 9 * kd + 3 ] = 2; /*O */

    /*CH */
    ncf[ 10 * kd + 1 ] = 1; /*C */
    ncf[ 10 * kd + 2 ] = 1; /*H */

    /*CO */
    ncf[ 11 * kd + 1 ] = 1; /*C */
    ncf[ 11 * kd + 3 ] = 1; /*O */

    /*HCO */
    ncf[ 12 * kd + 2 ] = 1; /*H */
    ncf[ 12 * kd + 1 ] = 1; /*C */
    ncf[ 12 * kd + 3 ] = 1; /*O */

    /*CH2O */
    ncf[ 13 * kd + 2 ] = 2; /*H */
    ncf[ 13 * kd + 1 ] = 1; /*C */
    ncf[ 13 * kd + 3 ] = 1; /*O */

    /*CH3 */
    ncf[ 14 * kd + 1 ] = 1; /*C */
    ncf[ 14 * kd + 2 ] = 3; /*H */

    /*CO2 */
    ncf[ 15 * kd + 1 ] = 1; /*C */
    ncf[ 15 * kd + 3 ] = 2; /*O */

    /*CH4 */
    ncf[ 16 * kd + 1 ] = 1; /*C */
    ncf[ 16 * kd + 2 ] = 4; /*H */

    /*C2H3 */
    ncf[ 17 * kd + 1 ] = 2; /*C */
    ncf[ 17 * kd + 2 ] = 3; /*H */

    /*C2H4 */
    ncf[ 18 * kd + 1 ] = 2; /*C */
    ncf[ 18 * kd + 2 ] = 4; /*H */

    /*C2H5 */
    ncf[ 19 * kd + 1 ] = 2; /*C */
    ncf[ 19 * kd + 2 ] = 5; /*H */

    /*C2H */
    ncf[ 20 * kd + 1 ] = 2; /*C */
    ncf[ 20 * kd + 2 ] = 1; /*H */

    /*HCCO */
    ncf[ 21 * kd + 2 ] = 1; /*H */
    ncf[ 21 * kd + 1 ] = 2; /*C */
    ncf[ 21 * kd + 3 ] = 1; /*O */

    /*C2H2 */
    ncf[ 22 * kd + 1 ] = 2; /*C */
    ncf[ 22 * kd + 2 ] = 2; /*H */

    /*C3H3 */
    ncf[ 23 * kd + 2 ] = 3; /*H */
    ncf[ 23 * kd + 1 ] = 3; /*C */

    /*A-C3H5 */
    ncf[ 24 * kd + 2 ] = 5; /*H */
    ncf[ 24 * kd + 1 ] = 3; /*C */

    /*N-C3H7 */
    ncf[ 25 * kd + 1 ] = 3; /*C */
    ncf[ 25 * kd + 2 ] = 7; /*H */

    /*C2H6 */
    ncf[ 26 * kd + 1 ] = 2; /*C */
    ncf[ 26 * kd + 2 ] = 6; /*H */

    /*P-C3H4 */
    ncf[ 27 * kd + 2 ] = 4; /*H */
    ncf[ 27 * kd + 1 ] = 3; /*C */

    /*A-C3H4 */
    ncf[ 28 * kd + 2 ] = 4; /*H */
    ncf[ 28 * kd + 1 ] = 3; /*C */

    /*A1 */
    ncf[ 29 * kd + 2 ] = 6; /*H */
    ncf[ 29 * kd + 1 ] = 6; /*C */

    /*A1- */
    ncf[ 30 * kd + 2 ] = 5; /*H */
    ncf[ 30 * kd + 1 ] = 6; /*C */

    /*C5H5 */
    ncf[ 31 * kd + 2 ] = 5; /*H */
    ncf[ 31 * kd + 1 ] = 5; /*C */

    /*C3H6 */
    ncf[ 32 * kd + 2 ] = 6; /*H */
    ncf[ 32 * kd + 1 ] = 3; /*C */

    /*C4H8 */
    ncf[ 33 * kd + 1 ] = 4; /*C */
    ncf[ 33 * kd + 2 ] = 8; /*H */

    /*C5H6 */
    ncf[ 34 * kd + 2 ] = 6; /*H */
    ncf[ 34 * kd + 1 ] = 5; /*C */

    /*A2 */
    ncf[ 35 * kd + 2 ] = 8; /*H */
    ncf[ 35 * kd + 1 ] = 10; /*C */

    /*C5H10 */
    ncf[ 36 * kd + 1 ] = 5; /*C */
    ncf[ 36 * kd + 2 ] = 10; /*H */

    /*C5H11 */
    ncf[ 37 * kd + 1 ] = 5; /*C */
    ncf[ 37 * kd + 2 ] = 11; /*H */

    /*A1C2H2 */
    ncf[ 38 * kd + 2 ] = 7; /*H */
    ncf[ 38 * kd + 1 ] = 8; /*C */

    /*A1CH2 */
    ncf[ 39 * kd + 2 ] = 7; /*H */
    ncf[ 39 * kd + 1 ] = 7; /*C */

    /*A1CHO */
    ncf[ 40 * kd + 2 ] = 6; /*H */
    ncf[ 40 * kd + 1 ] = 7; /*C */
    ncf[ 40 * kd + 3 ] = 1; /*O */

    /*A1CH3 */
    ncf[ 41 * kd + 2 ] = 8; /*H */
    ncf[ 41 * kd + 1 ] = 7; /*C */

    /*C7H15 */
    ncf[ 42 * kd + 1 ] = 7; /*C */
    ncf[ 42 * kd + 2 ] = 15; /*H */

    /*N-C7H16 */
    ncf[ 43 * kd + 1 ] = 7; /*C */
    ncf[ 43 * kd + 2 ] = 16; /*H */

    /*A1C2H* */
    ncf[ 44 * kd + 2 ] = 5; /*H */
    ncf[ 44 * kd + 1 ] = 8; /*C */

    /*A2- */
    ncf[ 45 * kd + 2 ] = 7; /*H */
    ncf[ 45 * kd + 1 ] = 10; /*C */

    /*A1C2H */
    ncf[ 46 * kd + 2 ] = 6; /*H */
    ncf[ 46 * kd + 1 ] = 8; /*C */

}


/* Returns the vector of strings of element names */
void CKSYME_STR(amrex::Vector<std::string>& ename)
{
    ename.resize(4);
    ename[0] = "N";
    ename[1] = "C";
    ename[2] = "H";
    ename[3] = "O";
}


/* Returns the vector of strings of species names */
void CKSYMS_STR(amrex::Vector<std::string>& kname)
{
    kname.resize(47);
    kname[0] = "N2";
    kname[1] = "S-CH2";
    kname[2] = "T-CH2";
    kname[3] = "O";
    kname[4] = "H2";
    kname[5] = "H";
    kname[6] = "OH";
    kname[7] = "H2O";
    kname[8] = "O2";
    kname[9] = "HO2";
    kname[10] = "CH";
    kname[11] = "CO";
    kname[12] = "HCO";
    kname[13] = "CH2O";
    kname[14] = "CH3";
    kname[15] = "CO2";
    kname[16] = "CH4";
    kname[17] = "C2H3";
    kname[18] = "C2H4";
    kname[19] = "C2H5";
    kname[20] = "C2H";
    kname[21] = "HCCO";
    kname[22] = "C2H2";
    kname[23] = "C3H3";
    kname[24] = "A-C3H5";
    kname[25] = "N-C3H7";
    kname[26] = "C2H6";
    kname[27] = "P-C3H4";
    kname[28] = "A-C3H4";
    kname[29] = "A1";
    kname[30] = "A1-";
    kname[31] = "C5H5";
    kname[32] = "C3H6";
    kname[33] = "C4H8";
    kname[34] = "C5H6";
    kname[35] = "A2";
    kname[36] = "C5H10";
    kname[37] = "C5H11";
    kname[38] = "A1C2H2";
    kname[39] = "A1CH2";
    kname[40] = "A1CHO";
    kname[41] = "A1CH3";
    kname[42] = "C7H15";
    kname[43] = "N-C7H16";
    kname[44] = "A1C2H*";
    kname[45] = "A2-";
    kname[46] = "A1C2H";
}

/*compute the sparsity pattern of the chemistry Jacobian */
void SPARSITY_INFO( int * nJdata, const int * consP, int NCELLS)
{
    amrex::GpuArray<amrex::Real,2304> Jac = {0.0};
    amrex::GpuArray<amrex::Real,47> conc = {0.0};
    for (int n=0; n<47; n++) {
        conc[n] = 1.0/ 47.000000 ;
    }
    aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

    int nJdata_tmp = 0;
    for (int k=0; k<48; k++) {
        for (int l=0; l<48; l++) {
            if(Jac[ 48 * k + l] != 0.0){
                nJdata_tmp = nJdata_tmp + 1;
            }
        }
    }

    *nJdata = NCELLS * nJdata_tmp;
}



/*compute the sparsity pattern of the system Jacobian */
void SPARSITY_INFO_SYST( int * nJdata, const int * consP, int NCELLS)
{
    amrex::GpuArray<amrex::Real,2304> Jac = {0.0};
    amrex::GpuArray<amrex::Real,47> conc = {0.0};
    for (int n=0; n<47; n++) {
        conc[n] = 1.0/ 47.000000 ;
    }
    aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

    int nJdata_tmp = 0;
    for (int k=0; k<48; k++) {
        for (int l=0; l<48; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(Jac[ 48 * k + l] != 0.0){
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
    amrex::GpuArray<amrex::Real,2304> Jac = {0.0};
    amrex::GpuArray<amrex::Real,47> conc = {0.0};
    for (int n=0; n<47; n++) {
        conc[n] = 1.0/ 47.000000 ;
    }
    aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

    int nJdata_tmp = 0;
    for (int k=0; k<48; k++) {
        for (int l=0; l<48; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(Jac[ 48 * k + l] != 0.0){
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
    amrex::GpuArray<amrex::Real,2304> Jac = {0.0};
    amrex::GpuArray<amrex::Real,47> conc = {0.0};
    for (int n=0; n<47; n++) {
        conc[n] = 1.0/ 47.000000 ;
    }
    aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

    colPtrs[0] = 0;
    int nJdata_tmp = 0;
    for (int nc=0; nc<NCELLS; nc++) {
        int offset_row = nc * 48;
        int offset_col = nc * 48;
        for (int k=0; k<48; k++) {
            for (int l=0; l<48; l++) {
                if(Jac[48*k + l] != 0.0) {
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
    amrex::GpuArray<amrex::Real,2304> Jac = {0.0};
    amrex::GpuArray<amrex::Real,47> conc = {0.0};
    for (int n=0; n<47; n++) {
        conc[n] = 1.0/ 47.000000 ;
    }
    aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

    if (base == 1) {
        rowPtrs[0] = 1;
        int nJdata_tmp = 1;
        for (int nc=0; nc<NCELLS; nc++) {
            int offset = nc * 48;
            for (int l=0; l<48; l++) {
                for (int k=0; k<48; k++) {
                    if(Jac[48*k + l] != 0.0) {
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
            int offset = nc * 48;
            for (int l=0; l<48; l++) {
                for (int k=0; k<48; k++) {
                    if(Jac[48*k + l] != 0.0) {
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
    amrex::GpuArray<amrex::Real,2304> Jac = {0.0};
    amrex::GpuArray<amrex::Real,47> conc = {0.0};
    for (int n=0; n<47; n++) {
        conc[n] = 1.0/ 47.000000 ;
    }
    aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

    if (base == 1) {
        rowPtr[0] = 1;
        int nJdata_tmp = 1;
        for (int nc=0; nc<NCELLS; nc++) {
            int offset = nc * 48;
            for (int l=0; l<48; l++) {
                for (int k=0; k<48; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp-1] = l+1 + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(Jac[48*k + l] != 0.0) {
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
            int offset = nc * 48;
            for (int l=0; l<48; l++) {
                for (int k=0; k<48; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp] = l + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(Jac[48*k + l] != 0.0) {
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
    amrex::GpuArray<amrex::Real,2304> Jac = {0.0};
    amrex::GpuArray<amrex::Real,47> conc = {0.0};
    for (int n=0; n<47; n++) {
        conc[n] = 1.0/ 47.000000 ;
    }
    aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

    colPtrs[0] = 0;
    int nJdata_tmp = 0;
    for (int k=0; k<48; k++) {
        for (int l=0; l<48; l++) {
            if (k == l) {
                rowVals[nJdata_tmp] = l; 
                indx[nJdata_tmp] = 48*k + l;
                nJdata_tmp = nJdata_tmp + 1; 
            } else {
                if(Jac[48*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l; 
                    indx[nJdata_tmp] = 48*k + l;
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
    amrex::GpuArray<amrex::Real,2304> Jac = {0.0};
    amrex::GpuArray<amrex::Real,47> conc = {0.0};
    for (int n=0; n<47; n++) {
        conc[n] = 1.0/ 47.000000 ;
    }
    aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

    if (base == 1) {
        rowPtr[0] = 1;
        int nJdata_tmp = 1;
        for (int l=0; l<48; l++) {
            for (int k=0; k<48; k++) {
                if (k == l) {
                    colVals[nJdata_tmp-1] = l+1; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(Jac[48*k + l] != 0.0) {
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
        for (int l=0; l<48; l++) {
            for (int k=0; k<48; k++) {
                if (k == l) {
                    colVals[nJdata_tmp] = l; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(Jac[48*k + l] != 0.0) {
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
