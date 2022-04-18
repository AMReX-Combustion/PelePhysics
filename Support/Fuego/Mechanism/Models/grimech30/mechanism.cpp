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
    for (id = 0; id < kd * 53; ++ id) {
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

}


/* Returns the vector of strings of element names */
void CKSYME_STR(amrex::Vector<std::string>& ename)
{
    ename.resize(5);
    ename[0] = "O";
    ename[1] = "H";
    ename[2] = "C";
    ename[3] = "N";
    ename[4] = "AR";
}


/* Returns the vector of strings of species names */
void CKSYMS_STR(amrex::Vector<std::string>& kname)
{
    kname.resize(53);
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
}

/*compute the sparsity pattern of the chemistry Jacobian */
void SPARSITY_INFO( int * nJdata, const int * consP, int NCELLS)
{
    amrex::GpuArray<amrex::Real,2916> Jac = {0.0};
    amrex::GpuArray<amrex::Real,53> conc = {0.0};
    for (int n=0; n<53; n++) {
        conc[n] = 1.0/ 53.000000 ;
    }
    aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

    int nJdata_tmp = 0;
    for (int k=0; k<54; k++) {
        for (int l=0; l<54; l++) {
            if(Jac[ 54 * k + l] != 0.0){
                nJdata_tmp = nJdata_tmp + 1;
            }
        }
    }

    *nJdata = NCELLS * nJdata_tmp;
}



/*compute the sparsity pattern of the system Jacobian */
void SPARSITY_INFO_SYST( int * nJdata, const int * consP, int NCELLS)
{
    amrex::GpuArray<amrex::Real,2916> Jac = {0.0};
    amrex::GpuArray<amrex::Real,53> conc = {0.0};
    for (int n=0; n<53; n++) {
        conc[n] = 1.0/ 53.000000 ;
    }
    aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

    int nJdata_tmp = 0;
    for (int k=0; k<54; k++) {
        for (int l=0; l<54; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(Jac[ 54 * k + l] != 0.0){
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
    amrex::GpuArray<amrex::Real,2916> Jac = {0.0};
    amrex::GpuArray<amrex::Real,53> conc = {0.0};
    for (int n=0; n<53; n++) {
        conc[n] = 1.0/ 53.000000 ;
    }
    aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

    int nJdata_tmp = 0;
    for (int k=0; k<54; k++) {
        for (int l=0; l<54; l++) {
            if(k == l){
                nJdata_tmp = nJdata_tmp + 1;
            } else {
                if(Jac[ 54 * k + l] != 0.0){
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
    amrex::GpuArray<amrex::Real,2916> Jac = {0.0};
    amrex::GpuArray<amrex::Real,53> conc = {0.0};
    for (int n=0; n<53; n++) {
        conc[n] = 1.0/ 53.000000 ;
    }
    aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

    colPtrs[0] = 0;
    int nJdata_tmp = 0;
    for (int nc=0; nc<NCELLS; nc++) {
        int offset_row = nc * 54;
        int offset_col = nc * 54;
        for (int k=0; k<54; k++) {
            for (int l=0; l<54; l++) {
                if(Jac[54*k + l] != 0.0) {
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
    amrex::GpuArray<amrex::Real,2916> Jac = {0.0};
    amrex::GpuArray<amrex::Real,53> conc = {0.0};
    for (int n=0; n<53; n++) {
        conc[n] = 1.0/ 53.000000 ;
    }
    aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

    if (base == 1) {
        rowPtrs[0] = 1;
        int nJdata_tmp = 1;
        for (int nc=0; nc<NCELLS; nc++) {
            int offset = nc * 54;
            for (int l=0; l<54; l++) {
                for (int k=0; k<54; k++) {
                    if(Jac[54*k + l] != 0.0) {
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
            int offset = nc * 54;
            for (int l=0; l<54; l++) {
                for (int k=0; k<54; k++) {
                    if(Jac[54*k + l] != 0.0) {
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
    amrex::GpuArray<amrex::Real,2916> Jac = {0.0};
    amrex::GpuArray<amrex::Real,53> conc = {0.0};
    for (int n=0; n<53; n++) {
        conc[n] = 1.0/ 53.000000 ;
    }
    aJacobian(&Jac[0], &conc[0], 1500.0, *consP);

    if (base == 1) {
        rowPtr[0] = 1;
        int nJdata_tmp = 1;
        for (int nc=0; nc<NCELLS; nc++) {
            int offset = nc * 54;
            for (int l=0; l<54; l++) {
                for (int k=0; k<54; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp-1] = l+1 + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(Jac[54*k + l] != 0.0) {
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
            int offset = nc * 54;
            for (int l=0; l<54; l++) {
                for (int k=0; k<54; k++) {
                    if (k == l) {
                        colVals[nJdata_tmp] = l + offset; 
                        nJdata_tmp = nJdata_tmp + 1; 
                    } else {
                        if(Jac[54*k + l] != 0.0) {
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
    amrex::GpuArray<amrex::Real,2916> Jac = {0.0};
    amrex::GpuArray<amrex::Real,53> conc = {0.0};
    for (int n=0; n<53; n++) {
        conc[n] = 1.0/ 53.000000 ;
    }
    aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

    colPtrs[0] = 0;
    int nJdata_tmp = 0;
    for (int k=0; k<54; k++) {
        for (int l=0; l<54; l++) {
            if (k == l) {
                rowVals[nJdata_tmp] = l; 
                indx[nJdata_tmp] = 54*k + l;
                nJdata_tmp = nJdata_tmp + 1; 
            } else {
                if(Jac[54*k + l] != 0.0) {
                    rowVals[nJdata_tmp] = l; 
                    indx[nJdata_tmp] = 54*k + l;
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
    amrex::GpuArray<amrex::Real,2916> Jac = {0.0};
    amrex::GpuArray<amrex::Real,53> conc = {0.0};
    for (int n=0; n<53; n++) {
        conc[n] = 1.0/ 53.000000 ;
    }
    aJacobian_precond(&Jac[0], &conc[0], 1500.0, *consP);

    if (base == 1) {
        rowPtr[0] = 1;
        int nJdata_tmp = 1;
        for (int l=0; l<54; l++) {
            for (int k=0; k<54; k++) {
                if (k == l) {
                    colVals[nJdata_tmp-1] = l+1; 
                    nJdata_tmp = nJdata_tmp + 1; 
                } else {
                    if(Jac[54*k + l] != 0.0) {
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
                    if(Jac[54*k + l] != 0.0) {
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
