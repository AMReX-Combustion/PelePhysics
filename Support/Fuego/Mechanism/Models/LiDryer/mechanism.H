#ifndef MECHANISM_H
#define MECHANISM_H

#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>

/* Elements
0  H
1  O
2  N
*/

// Species
#define H2_ID 0
#define O2_ID 1
#define H2O_ID 2
#define H_ID 3
#define O_ID 4
#define OH_ID 5
#define HO2_ID 6
#define H2O2_ID 7
#define N2_ID 8

#define NUM_ELEMENTS 3
#define NUM_SPECIES 9
#define NUM_REACTIONS 21

#define NUM_FIT 4


/* ALWAYS on CPU stuff -- can have different def depending on if we are CPU or GPU based. Defined in mechanism.cpp  */
void atomicWeight(amrex::Real *  awt);
/* MISC  */
void CKAWT(amrex::Real *  awt);
void CKNCF(int * ncf);
void CKSYME_STR(amrex::Vector<std::string>& ename);
void CKSYMS_STR(amrex::Vector<std::string>& kname);
/* SPARSE INFORMATION  */
void SPARSITY_INFO(int * nJdata, const int * consP, int NCELLS);
void SPARSITY_INFO_SYST(int * nJdata, const int * consP, int NCELLS);
void SPARSITY_INFO_SYST_SIMPLIFIED(int * nJdata, const int * consP);
void SPARSITY_PREPROC_CSC(int * rowVals, int * colPtrs, const int * consP, int NCELLS);
void SPARSITY_PREPROC_CSR(int * colVals, int * rowPtrs, const int * consP, int NCELLS, int base);
void SPARSITY_PREPROC_SYST_CSR(int * colVals, int * rowPtrs, const int * consP, int NCELLS, int base);
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(int * rowVals, int * colPtrs, int * indx, const int * consP);
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(int * colVals, int * rowPtr, const int * consP, int base);



/*A few mechanism parameters */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKINDX(int * mm, int * kk, int * ii, int * nfit)
{
    *mm = 3;
    *kk = 9;
    *ii = 21;
    *nfit = -1; /*Why do you need this anyway ?  */
}


/* inverse molecular weights  */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void get_imw(amrex::Real *imw_new){
    imw_new[0] = 1.0/2.015940;/*H2 */
    imw_new[1] = 1.0/31.998800;/*O2 */
    imw_new[2] = 1.0/18.015340;/*H2O */
    imw_new[3] = 1.0/1.007970;/*H */
    imw_new[4] = 1.0/15.999400;/*O */
    imw_new[5] = 1.0/17.007370;/*OH */
    imw_new[6] = 1.0/33.006770;/*HO2 */
    imw_new[7] = 1.0/34.014740;/*H2O2 */
    imw_new[8] = 1.0/28.013400;/*N2 */
}

/* molecular weights  */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void get_mw(amrex::Real *mw_new){
    mw_new[0] = 2.015940;/*H2 */
    mw_new[1] = 31.998800;/*O2 */
    mw_new[2] = 18.015340;/*H2O */
    mw_new[3] = 1.007970;/*H */
    mw_new[4] = 15.999400;/*O */
    mw_new[5] = 17.007370;/*OH */
    mw_new[6] = 33.006770;/*HO2 */
    mw_new[7] = 34.014740;/*H2O2 */
    mw_new[8] = 28.013400;/*N2 */
}



/* Returns R, Rc, Patm */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKRP(amrex::Real *  ru, amrex::Real *  ruc, amrex::Real *  pa)
{
     *ru  = 8.31446261815324e+07; 
     *ruc = 1.98721558317399615845; 
     *pa  = 1.01325e+06; 
}


/*compute Cv/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void cv_R(amrex::Real * species, const amrex::Real *  tc)
{

    /*temperature */
    const amrex::Real T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            +2.29812431e+00
            +8.24944174e-04 * tc[1]
            -8.14301529e-07 * tc[2]
            -9.47543433e-11 * tc[3]
            +4.13487224e-13 * tc[4];
        /*species 1: O2 */
        species[1] =
            +2.21293640e+00
            +1.12748635e-03 * tc[1]
            -5.75615047e-07 * tc[2]
            +1.31387723e-09 * tc[3]
            -8.76855392e-13 * tc[4];
        /*species 2: H2O */
        species[2] =
            +2.38684249e+00
            +3.47498246e-03 * tc[1]
            -6.35469633e-06 * tc[2]
            +6.96858127e-09 * tc[3]
            -2.50658847e-12 * tc[4];
        /*species 3: H */
        species[3] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 4: O */
        species[4] =
            +1.94642878e+00
            -1.63816649e-03 * tc[1]
            +2.42103170e-06 * tc[2]
            -1.60284319e-09 * tc[3]
            +3.89069636e-13 * tc[4];
        /*species 5: OH */
        species[5] =
            +3.12530561e+00
            -3.22544939e-03 * tc[1]
            +6.52764691e-06 * tc[2]
            -5.79853643e-09 * tc[3]
            +2.06237379e-12 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +3.30179801e+00
            -4.74912051e-03 * tc[1]
            +2.11582891e-05 * tc[2]
            -2.42763894e-08 * tc[3]
            +9.29225124e-12 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +2.38875365e+00
            +6.56922581e-03 * tc[1]
            -1.48501258e-07 * tc[2]
            -4.62580552e-09 * tc[3]
            +2.47151475e-12 * tc[4];
        /*species 8: N2 */
        species[8] =
            +2.29867700e+00
            +1.40824000e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485500e-12 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            +1.99142337e+00
            +7.00064411e-04 * tc[1]
            -5.63382869e-08 * tc[2]
            -9.23157818e-12 * tc[3]
            +1.58275179e-15 * tc[4];
        /*species 1: O2 */
        species[1] =
            +2.69757819e+00
            +6.13519689e-04 * tc[1]
            -1.25884199e-07 * tc[2]
            +1.77528148e-11 * tc[3]
            -1.13643531e-15 * tc[4];
        /*species 2: H2O */
        species[2] =
            +1.67214561e+00
            +3.05629289e-03 * tc[1]
            -8.73026011e-07 * tc[2]
            +1.20099639e-10 * tc[3]
            -6.39161787e-15 * tc[4];
        /*species 3: H */
        species[3] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 4: O */
        species[4] =
            +1.54205966e+00
            -2.75506191e-05 * tc[1]
            -3.10280335e-09 * tc[2]
            +4.55106742e-12 * tc[3]
            -4.36805150e-16 * tc[4];
        /*species 5: OH */
        species[5] =
            +1.86472886e+00
            +1.05650448e-03 * tc[1]
            -2.59082758e-07 * tc[2]
            +3.05218674e-11 * tc[3]
            -1.33195876e-15 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +3.01721090e+00
            +2.23982013e-03 * tc[1]
            -6.33658150e-07 * tc[2]
            +1.14246370e-10 * tc[3]
            -1.07908535e-14 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +3.57316685e+00
            +4.33613639e-03 * tc[1]
            -1.47468882e-06 * tc[2]
            +2.34890357e-10 * tc[3]
            -1.43165356e-14 * tc[4];
        /*species 8: N2 */
        species[8] =
            +1.92664000e+00
            +1.48797700e-03 * tc[1]
            -5.68476100e-07 * tc[2]
            +1.00970400e-10 * tc[3]
            -6.75335100e-15 * tc[4];
    }
    return;
}


/*compute Cp/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void cp_R(amrex::Real * species, const amrex::Real *  tc)
{

    /*temperature */
    const amrex::Real T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            +3.29812431e+00
            +8.24944174e-04 * tc[1]
            -8.14301529e-07 * tc[2]
            -9.47543433e-11 * tc[3]
            +4.13487224e-13 * tc[4];
        /*species 1: O2 */
        species[1] =
            +3.21293640e+00
            +1.12748635e-03 * tc[1]
            -5.75615047e-07 * tc[2]
            +1.31387723e-09 * tc[3]
            -8.76855392e-13 * tc[4];
        /*species 2: H2O */
        species[2] =
            +3.38684249e+00
            +3.47498246e-03 * tc[1]
            -6.35469633e-06 * tc[2]
            +6.96858127e-09 * tc[3]
            -2.50658847e-12 * tc[4];
        /*species 3: H */
        species[3] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 4: O */
        species[4] =
            +2.94642878e+00
            -1.63816649e-03 * tc[1]
            +2.42103170e-06 * tc[2]
            -1.60284319e-09 * tc[3]
            +3.89069636e-13 * tc[4];
        /*species 5: OH */
        species[5] =
            +4.12530561e+00
            -3.22544939e-03 * tc[1]
            +6.52764691e-06 * tc[2]
            -5.79853643e-09 * tc[3]
            +2.06237379e-12 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +4.30179801e+00
            -4.74912051e-03 * tc[1]
            +2.11582891e-05 * tc[2]
            -2.42763894e-08 * tc[3]
            +9.29225124e-12 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +3.38875365e+00
            +6.56922581e-03 * tc[1]
            -1.48501258e-07 * tc[2]
            -4.62580552e-09 * tc[3]
            +2.47151475e-12 * tc[4];
        /*species 8: N2 */
        species[8] =
            +3.29867700e+00
            +1.40824000e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485500e-12 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            +2.99142337e+00
            +7.00064411e-04 * tc[1]
            -5.63382869e-08 * tc[2]
            -9.23157818e-12 * tc[3]
            +1.58275179e-15 * tc[4];
        /*species 1: O2 */
        species[1] =
            +3.69757819e+00
            +6.13519689e-04 * tc[1]
            -1.25884199e-07 * tc[2]
            +1.77528148e-11 * tc[3]
            -1.13643531e-15 * tc[4];
        /*species 2: H2O */
        species[2] =
            +2.67214561e+00
            +3.05629289e-03 * tc[1]
            -8.73026011e-07 * tc[2]
            +1.20099639e-10 * tc[3]
            -6.39161787e-15 * tc[4];
        /*species 3: H */
        species[3] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4];
        /*species 4: O */
        species[4] =
            +2.54205966e+00
            -2.75506191e-05 * tc[1]
            -3.10280335e-09 * tc[2]
            +4.55106742e-12 * tc[3]
            -4.36805150e-16 * tc[4];
        /*species 5: OH */
        species[5] =
            +2.86472886e+00
            +1.05650448e-03 * tc[1]
            -2.59082758e-07 * tc[2]
            +3.05218674e-11 * tc[3]
            -1.33195876e-15 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +4.01721090e+00
            +2.23982013e-03 * tc[1]
            -6.33658150e-07 * tc[2]
            +1.14246370e-10 * tc[3]
            -1.07908535e-14 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +4.57316685e+00
            +4.33613639e-03 * tc[1]
            -1.47468882e-06 * tc[2]
            +2.34890357e-10 * tc[3]
            -1.43165356e-14 * tc[4];
        /*species 8: N2 */
        species[8] =
            +2.92664000e+00
            +1.48797700e-03 * tc[1]
            -5.68476100e-07 * tc[2]
            +1.00970400e-10 * tc[3]
            -6.75335100e-15 * tc[4];
    }
    return;
}


/*compute the g/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void gibbs(amrex::Real * species, const amrex::Real *  tc)
{

    /*temperature */
    const amrex::Real T = tc[1];
    const amrex::Real invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            -1.012520870000000e+03 * invT
            +6.592218400000000e+00
            -3.298124310000000e+00 * tc[0]
            -4.124720870000000e-04 * tc[1]
            +1.357169215000000e-07 * tc[2]
            +7.896195275000000e-12 * tc[3]
            -2.067436120000000e-14 * tc[4];
        /*species 1: O2 */
        species[1] =
            -1.005249020000000e+03 * invT
            -2.821801190000000e+00
            -3.212936400000000e+00 * tc[0]
            -5.637431750000000e-04 * tc[1]
            +9.593584116666666e-08 * tc[2]
            -1.094897691666667e-10 * tc[3]
            +4.384276960000000e-14 * tc[4];
        /*species 2: H2O */
        species[2] =
            -3.020811330000000e+04 * invT
            +7.966096399999998e-01
            -3.386842490000000e+00 * tc[0]
            -1.737491230000000e-03 * tc[1]
            +1.059116055000000e-06 * tc[2]
            -5.807151058333333e-10 * tc[3]
            +1.253294235000000e-13 * tc[4];
        /*species 3: H */
        species[3] =
            +2.547162700000000e+04 * invT
            +2.960117608000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 4: O */
        species[4] =
            +2.914764450000000e+04 * invT
            -1.756619999999964e-02
            -2.946428780000000e+00 * tc[0]
            +8.190832450000000e-04 * tc[1]
            -4.035052833333333e-07 * tc[2]
            +1.335702658333333e-10 * tc[3]
            -1.945348180000000e-14 * tc[4];
        /*species 5: OH */
        species[5] =
            +3.346309130000000e+03 * invT
            +4.815738570000000e+00
            -4.125305610000000e+00 * tc[0]
            +1.612724695000000e-03 * tc[1]
            -1.087941151666667e-06 * tc[2]
            +4.832113691666666e-10 * tc[3]
            -1.031186895000000e-13 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +2.948080400000000e+02 * invT
            +5.851355599999999e-01
            -4.301798010000000e+00 * tc[0]
            +2.374560255000000e-03 * tc[1]
            -3.526381516666666e-06 * tc[2]
            +2.023032450000000e-09 * tc[3]
            -4.646125620000001e-13 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            -1.766314650000000e+04 * invT
            -3.396609550000000e+00
            -3.388753650000000e+00 * tc[0]
            -3.284612905000000e-03 * tc[1]
            +2.475020966666667e-08 * tc[2]
            +3.854837933333333e-10 * tc[3]
            -1.235757375000000e-13 * tc[4];
        /*species 8: N2 */
        species[8] =
            -1.020900000000000e+03 * invT
            -6.516950000000001e-01
            -3.298677000000000e+00 * tc[0]
            -7.041200000000000e-04 * tc[1]
            +6.605369999999999e-07 * tc[2]
            -4.701262500000001e-10 * tc[3]
            +1.222427500000000e-13 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            -8.350339970000000e+02 * invT
            +4.346533540000000e+00
            -2.991423370000000e+00 * tc[0]
            -3.500322055000000e-04 * tc[1]
            +9.389714483333333e-09 * tc[2]
            +7.692981816666667e-13 * tc[3]
            -7.913758950000000e-17 * tc[4];
        /*species 1: O2 */
        species[1] =
            -1.233930180000000e+03 * invT
            +5.084126000000002e-01
            -3.697578190000000e+00 * tc[0]
            -3.067598445000000e-04 * tc[1]
            +2.098069983333333e-08 * tc[2]
            -1.479401233333333e-12 * tc[3]
            +5.682176550000000e-17 * tc[4];
        /*species 2: H2O */
        species[2] =
            -2.989920900000000e+04 * invT
            -4.190671200000001e+00
            -2.672145610000000e+00 * tc[0]
            -1.528146445000000e-03 * tc[1]
            +1.455043351666667e-07 * tc[2]
            -1.000830325000000e-11 * tc[3]
            +3.195808935000000e-16 * tc[4];
        /*species 3: H */
        species[3] =
            +2.547162700000000e+04 * invT
            +2.960117638000000e+00
            -2.500000000000000e+00 * tc[0]
            -0.000000000000000e+00 * tc[1]
            -0.000000000000000e+00 * tc[2]
            -0.000000000000000e+00 * tc[3]
            -0.000000000000000e+00 * tc[4];
        /*species 4: O */
        species[4] =
            +2.923080270000000e+04 * invT
            -2.378248450000000e+00
            -2.542059660000000e+00 * tc[0]
            +1.377530955000000e-05 * tc[1]
            +5.171338916666667e-10 * tc[2]
            -3.792556183333333e-13 * tc[3]
            +2.184025750000000e-17 * tc[4];
        /*species 5: OH */
        species[5] =
            +3.683628750000000e+03 * invT
            -2.836911870000000e+00
            -2.864728860000000e+00 * tc[0]
            -5.282522400000000e-04 * tc[1]
            +4.318045966666667e-08 * tc[2]
            -2.543488950000000e-12 * tc[3]
            +6.659793800000000e-17 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +1.118567130000000e+02 * invT
            +2.321087500000001e-01
            -4.017210900000000e+00 * tc[0]
            -1.119910065000000e-03 * tc[1]
            +1.056096916666667e-07 * tc[2]
            -9.520530833333334e-12 * tc[3]
            +5.395426750000000e-16 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            -1.800696090000000e+04 * invT
            +4.072029891000000e+00
            -4.573166850000000e+00 * tc[0]
            -2.168068195000000e-03 * tc[1]
            +2.457814700000000e-07 * tc[2]
            -1.957419641666667e-11 * tc[3]
            +7.158267800000000e-16 * tc[4];
        /*species 8: N2 */
        species[8] =
            -9.227977000000000e+02 * invT
            -3.053888000000000e+00
            -2.926640000000000e+00 * tc[0]
            -7.439885000000000e-04 * tc[1]
            +9.474601666666666e-08 * tc[2]
            -8.414199999999999e-12 * tc[3]
            +3.376675500000000e-16 * tc[4];
    }
    return;
}


/*compute the a/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void helmholtz(amrex::Real * species, const amrex::Real *  tc)
{

    /*temperature */
    const amrex::Real T = tc[1];
    const amrex::Real invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            -1.01252087e+03 * invT
            +5.59221840e+00
            -3.29812431e+00 * tc[0]
            -4.12472087e-04 * tc[1]
            +1.35716922e-07 * tc[2]
            +7.89619527e-12 * tc[3]
            -2.06743612e-14 * tc[4];
        /*species 1: O2 */
        species[1] =
            -1.00524902e+03 * invT
            -3.82180119e+00
            -3.21293640e+00 * tc[0]
            -5.63743175e-04 * tc[1]
            +9.59358412e-08 * tc[2]
            -1.09489769e-10 * tc[3]
            +4.38427696e-14 * tc[4];
        /*species 2: H2O */
        species[2] =
            -3.02081133e+04 * invT
            -2.03390360e-01
            -3.38684249e+00 * tc[0]
            -1.73749123e-03 * tc[1]
            +1.05911606e-06 * tc[2]
            -5.80715106e-10 * tc[3]
            +1.25329424e-13 * tc[4];
        /*species 3: H */
        species[3] =
            +2.54716270e+04 * invT
            +1.96011761e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 4: O */
        species[4] =
            +2.91476445e+04 * invT
            -1.01756620e+00
            -2.94642878e+00 * tc[0]
            +8.19083245e-04 * tc[1]
            -4.03505283e-07 * tc[2]
            +1.33570266e-10 * tc[3]
            -1.94534818e-14 * tc[4];
        /*species 5: OH */
        species[5] =
            +3.34630913e+03 * invT
            +3.81573857e+00
            -4.12530561e+00 * tc[0]
            +1.61272470e-03 * tc[1]
            -1.08794115e-06 * tc[2]
            +4.83211369e-10 * tc[3]
            -1.03118689e-13 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +2.94808040e+02 * invT
            -4.14864440e-01
            -4.30179801e+00 * tc[0]
            +2.37456025e-03 * tc[1]
            -3.52638152e-06 * tc[2]
            +2.02303245e-09 * tc[3]
            -4.64612562e-13 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            -1.76631465e+04 * invT
            -4.39660955e+00
            -3.38875365e+00 * tc[0]
            -3.28461290e-03 * tc[1]
            +2.47502097e-08 * tc[2]
            +3.85483793e-10 * tc[3]
            -1.23575738e-13 * tc[4];
        /*species 8: N2 */
        species[8] =
            -1.02090000e+03 * invT
            -1.65169500e+00
            -3.29867700e+00 * tc[0]
            -7.04120000e-04 * tc[1]
            +6.60537000e-07 * tc[2]
            -4.70126250e-10 * tc[3]
            +1.22242750e-13 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            -8.35033997e+02 * invT
            +3.34653354e+00
            -2.99142337e+00 * tc[0]
            -3.50032206e-04 * tc[1]
            +9.38971448e-09 * tc[2]
            +7.69298182e-13 * tc[3]
            -7.91375895e-17 * tc[4];
        /*species 1: O2 */
        species[1] =
            -1.23393018e+03 * invT
            -4.91587400e-01
            -3.69757819e+00 * tc[0]
            -3.06759845e-04 * tc[1]
            +2.09806998e-08 * tc[2]
            -1.47940123e-12 * tc[3]
            +5.68217655e-17 * tc[4];
        /*species 2: H2O */
        species[2] =
            -2.98992090e+04 * invT
            -5.19067120e+00
            -2.67214561e+00 * tc[0]
            -1.52814644e-03 * tc[1]
            +1.45504335e-07 * tc[2]
            -1.00083033e-11 * tc[3]
            +3.19580894e-16 * tc[4];
        /*species 3: H */
        species[3] =
            +2.54716270e+04 * invT
            +1.96011764e+00
            -2.50000000e+00 * tc[0]
            -0.00000000e+00 * tc[1]
            -0.00000000e+00 * tc[2]
            -0.00000000e+00 * tc[3]
            -0.00000000e+00 * tc[4];
        /*species 4: O */
        species[4] =
            +2.92308027e+04 * invT
            -3.37824845e+00
            -2.54205966e+00 * tc[0]
            +1.37753096e-05 * tc[1]
            +5.17133892e-10 * tc[2]
            -3.79255618e-13 * tc[3]
            +2.18402575e-17 * tc[4];
        /*species 5: OH */
        species[5] =
            +3.68362875e+03 * invT
            -3.83691187e+00
            -2.86472886e+00 * tc[0]
            -5.28252240e-04 * tc[1]
            +4.31804597e-08 * tc[2]
            -2.54348895e-12 * tc[3]
            +6.65979380e-17 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +1.11856713e+02 * invT
            -7.67891250e-01
            -4.01721090e+00 * tc[0]
            -1.11991006e-03 * tc[1]
            +1.05609692e-07 * tc[2]
            -9.52053083e-12 * tc[3]
            +5.39542675e-16 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            -1.80069609e+04 * invT
            +3.07202989e+00
            -4.57316685e+00 * tc[0]
            -2.16806820e-03 * tc[1]
            +2.45781470e-07 * tc[2]
            -1.95741964e-11 * tc[3]
            +7.15826780e-16 * tc[4];
        /*species 8: N2 */
        species[8] =
            -9.22797700e+02 * invT
            -4.05388800e+00
            -2.92664000e+00 * tc[0]
            -7.43988500e-04 * tc[1]
            +9.47460167e-08 * tc[2]
            -8.41420000e-12 * tc[3]
            +3.37667550e-16 * tc[4];
    }
    return;
}


/*compute the e/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void speciesInternalEnergy(amrex::Real * species, const amrex::Real *  tc)
{

    /*temperature */
    const amrex::Real T = tc[1];
    const amrex::Real invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            +2.29812431e+00
            +4.12472087e-04 * tc[1]
            -2.71433843e-07 * tc[2]
            -2.36885858e-11 * tc[3]
            +8.26974448e-14 * tc[4]
            -1.01252087e+03 * invT;
        /*species 1: O2 */
        species[1] =
            +2.21293640e+00
            +5.63743175e-04 * tc[1]
            -1.91871682e-07 * tc[2]
            +3.28469308e-10 * tc[3]
            -1.75371078e-13 * tc[4]
            -1.00524902e+03 * invT;
        /*species 2: H2O */
        species[2] =
            +2.38684249e+00
            +1.73749123e-03 * tc[1]
            -2.11823211e-06 * tc[2]
            +1.74214532e-09 * tc[3]
            -5.01317694e-13 * tc[4]
            -3.02081133e+04 * invT;
        /*species 3: H */
        species[3] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716270e+04 * invT;
        /*species 4: O */
        species[4] =
            +1.94642878e+00
            -8.19083245e-04 * tc[1]
            +8.07010567e-07 * tc[2]
            -4.00710797e-10 * tc[3]
            +7.78139272e-14 * tc[4]
            +2.91476445e+04 * invT;
        /*species 5: OH */
        species[5] =
            +3.12530561e+00
            -1.61272470e-03 * tc[1]
            +2.17588230e-06 * tc[2]
            -1.44963411e-09 * tc[3]
            +4.12474758e-13 * tc[4]
            +3.34630913e+03 * invT;
        /*species 6: HO2 */
        species[6] =
            +3.30179801e+00
            -2.37456025e-03 * tc[1]
            +7.05276303e-06 * tc[2]
            -6.06909735e-09 * tc[3]
            +1.85845025e-12 * tc[4]
            +2.94808040e+02 * invT;
        /*species 7: H2O2 */
        species[7] =
            +2.38875365e+00
            +3.28461290e-03 * tc[1]
            -4.95004193e-08 * tc[2]
            -1.15645138e-09 * tc[3]
            +4.94302950e-13 * tc[4]
            -1.76631465e+04 * invT;
        /*species 8: N2 */
        species[8] =
            +2.29867700e+00
            +7.04120000e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88971000e-13 * tc[4]
            -1.02090000e+03 * invT;
    } else {
        /*species 0: H2 */
        species[0] =
            +1.99142337e+00
            +3.50032206e-04 * tc[1]
            -1.87794290e-08 * tc[2]
            -2.30789455e-12 * tc[3]
            +3.16550358e-16 * tc[4]
            -8.35033997e+02 * invT;
        /*species 1: O2 */
        species[1] =
            +2.69757819e+00
            +3.06759845e-04 * tc[1]
            -4.19613997e-08 * tc[2]
            +4.43820370e-12 * tc[3]
            -2.27287062e-16 * tc[4]
            -1.23393018e+03 * invT;
        /*species 2: H2O */
        species[2] =
            +1.67214561e+00
            +1.52814644e-03 * tc[1]
            -2.91008670e-07 * tc[2]
            +3.00249098e-11 * tc[3]
            -1.27832357e-15 * tc[4]
            -2.98992090e+04 * invT;
        /*species 3: H */
        species[3] =
            +1.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716270e+04 * invT;
        /*species 4: O */
        species[4] =
            +1.54205966e+00
            -1.37753096e-05 * tc[1]
            -1.03426778e-09 * tc[2]
            +1.13776685e-12 * tc[3]
            -8.73610300e-17 * tc[4]
            +2.92308027e+04 * invT;
        /*species 5: OH */
        species[5] =
            +1.86472886e+00
            +5.28252240e-04 * tc[1]
            -8.63609193e-08 * tc[2]
            +7.63046685e-12 * tc[3]
            -2.66391752e-16 * tc[4]
            +3.68362875e+03 * invT;
        /*species 6: HO2 */
        species[6] =
            +3.01721090e+00
            +1.11991006e-03 * tc[1]
            -2.11219383e-07 * tc[2]
            +2.85615925e-11 * tc[3]
            -2.15817070e-15 * tc[4]
            +1.11856713e+02 * invT;
        /*species 7: H2O2 */
        species[7] =
            +3.57316685e+00
            +2.16806820e-03 * tc[1]
            -4.91562940e-07 * tc[2]
            +5.87225893e-11 * tc[3]
            -2.86330712e-15 * tc[4]
            -1.80069609e+04 * invT;
        /*species 8: N2 */
        species[8] =
            +1.92664000e+00
            +7.43988500e-04 * tc[1]
            -1.89492033e-07 * tc[2]
            +2.52426000e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
    }
    return;
}


/*compute the h/(RT) at the given temperature (Eq 20) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void speciesEnthalpy(amrex::Real * species, const amrex::Real *  tc)
{

    /*temperature */
    const amrex::Real T = tc[1];
    const amrex::Real invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            +3.29812431e+00
            +4.12472087e-04 * tc[1]
            -2.71433843e-07 * tc[2]
            -2.36885858e-11 * tc[3]
            +8.26974448e-14 * tc[4]
            -1.01252087e+03 * invT;
        /*species 1: O2 */
        species[1] =
            +3.21293640e+00
            +5.63743175e-04 * tc[1]
            -1.91871682e-07 * tc[2]
            +3.28469308e-10 * tc[3]
            -1.75371078e-13 * tc[4]
            -1.00524902e+03 * invT;
        /*species 2: H2O */
        species[2] =
            +3.38684249e+00
            +1.73749123e-03 * tc[1]
            -2.11823211e-06 * tc[2]
            +1.74214532e-09 * tc[3]
            -5.01317694e-13 * tc[4]
            -3.02081133e+04 * invT;
        /*species 3: H */
        species[3] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716270e+04 * invT;
        /*species 4: O */
        species[4] =
            +2.94642878e+00
            -8.19083245e-04 * tc[1]
            +8.07010567e-07 * tc[2]
            -4.00710797e-10 * tc[3]
            +7.78139272e-14 * tc[4]
            +2.91476445e+04 * invT;
        /*species 5: OH */
        species[5] =
            +4.12530561e+00
            -1.61272470e-03 * tc[1]
            +2.17588230e-06 * tc[2]
            -1.44963411e-09 * tc[3]
            +4.12474758e-13 * tc[4]
            +3.34630913e+03 * invT;
        /*species 6: HO2 */
        species[6] =
            +4.30179801e+00
            -2.37456025e-03 * tc[1]
            +7.05276303e-06 * tc[2]
            -6.06909735e-09 * tc[3]
            +1.85845025e-12 * tc[4]
            +2.94808040e+02 * invT;
        /*species 7: H2O2 */
        species[7] =
            +3.38875365e+00
            +3.28461290e-03 * tc[1]
            -4.95004193e-08 * tc[2]
            -1.15645138e-09 * tc[3]
            +4.94302950e-13 * tc[4]
            -1.76631465e+04 * invT;
        /*species 8: N2 */
        species[8] =
            +3.29867700e+00
            +7.04120000e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88971000e-13 * tc[4]
            -1.02090000e+03 * invT;
    } else {
        /*species 0: H2 */
        species[0] =
            +2.99142337e+00
            +3.50032206e-04 * tc[1]
            -1.87794290e-08 * tc[2]
            -2.30789455e-12 * tc[3]
            +3.16550358e-16 * tc[4]
            -8.35033997e+02 * invT;
        /*species 1: O2 */
        species[1] =
            +3.69757819e+00
            +3.06759845e-04 * tc[1]
            -4.19613997e-08 * tc[2]
            +4.43820370e-12 * tc[3]
            -2.27287062e-16 * tc[4]
            -1.23393018e+03 * invT;
        /*species 2: H2O */
        species[2] =
            +2.67214561e+00
            +1.52814644e-03 * tc[1]
            -2.91008670e-07 * tc[2]
            +3.00249098e-11 * tc[3]
            -1.27832357e-15 * tc[4]
            -2.98992090e+04 * invT;
        /*species 3: H */
        species[3] =
            +2.50000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            +2.54716270e+04 * invT;
        /*species 4: O */
        species[4] =
            +2.54205966e+00
            -1.37753096e-05 * tc[1]
            -1.03426778e-09 * tc[2]
            +1.13776685e-12 * tc[3]
            -8.73610300e-17 * tc[4]
            +2.92308027e+04 * invT;
        /*species 5: OH */
        species[5] =
            +2.86472886e+00
            +5.28252240e-04 * tc[1]
            -8.63609193e-08 * tc[2]
            +7.63046685e-12 * tc[3]
            -2.66391752e-16 * tc[4]
            +3.68362875e+03 * invT;
        /*species 6: HO2 */
        species[6] =
            +4.01721090e+00
            +1.11991006e-03 * tc[1]
            -2.11219383e-07 * tc[2]
            +2.85615925e-11 * tc[3]
            -2.15817070e-15 * tc[4]
            +1.11856713e+02 * invT;
        /*species 7: H2O2 */
        species[7] =
            +4.57316685e+00
            +2.16806820e-03 * tc[1]
            -4.91562940e-07 * tc[2]
            +5.87225893e-11 * tc[3]
            -2.86330712e-15 * tc[4]
            -1.80069609e+04 * invT;
        /*species 8: N2 */
        species[8] =
            +2.92664000e+00
            +7.43988500e-04 * tc[1]
            -1.89492033e-07 * tc[2]
            +2.52426000e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
    }
    return;
}


/*compute the S/R at the given temperature (Eq 21) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void speciesEntropy(amrex::Real * species, const amrex::Real *  tc)
{

    /*temperature */
    const amrex::Real T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            +3.29812431e+00 * tc[0]
            +8.24944174e-04 * tc[1]
            -4.07150765e-07 * tc[2]
            -3.15847811e-11 * tc[3]
            +1.03371806e-13 * tc[4]
            -3.29409409e+00 ;
        /*species 1: O2 */
        species[1] =
            +3.21293640e+00 * tc[0]
            +1.12748635e-03 * tc[1]
            -2.87807523e-07 * tc[2]
            +4.37959077e-10 * tc[3]
            -2.19213848e-13 * tc[4]
            +6.03473759e+00 ;
        /*species 2: H2O */
        species[2] =
            +3.38684249e+00 * tc[0]
            +3.47498246e-03 * tc[1]
            -3.17734817e-06 * tc[2]
            +2.32286042e-09 * tc[3]
            -6.26647117e-13 * tc[4]
            +2.59023285e+00 ;
        /*species 3: H */
        species[3] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -4.60117608e-01 ;
        /*species 4: O */
        species[4] =
            +2.94642878e+00 * tc[0]
            -1.63816649e-03 * tc[1]
            +1.21051585e-06 * tc[2]
            -5.34281063e-10 * tc[3]
            +9.72674090e-14 * tc[4]
            +2.96399498e+00 ;
        /*species 5: OH */
        species[5] =
            +4.12530561e+00 * tc[0]
            -3.22544939e-03 * tc[1]
            +3.26382346e-06 * tc[2]
            -1.93284548e-09 * tc[3]
            +5.15593447e-13 * tc[4]
            -6.90432960e-01 ;
        /*species 6: HO2 */
        species[6] =
            +4.30179801e+00 * tc[0]
            -4.74912051e-03 * tc[1]
            +1.05791445e-05 * tc[2]
            -8.09212980e-09 * tc[3]
            +2.32306281e-12 * tc[4]
            +3.71666245e+00 ;
        /*species 7: H2O2 */
        species[7] =
            +3.38875365e+00 * tc[0]
            +6.56922581e-03 * tc[1]
            -7.42506290e-08 * tc[2]
            -1.54193517e-09 * tc[3]
            +6.17878688e-13 * tc[4]
            +6.78536320e+00 ;
        /*species 8: N2 */
        species[8] =
            +3.29867700e+00 * tc[0]
            +1.40824000e-03 * tc[1]
            -1.98161100e-06 * tc[2]
            +1.88050500e-09 * tc[3]
            -6.11213750e-13 * tc[4]
            +3.95037200e+00 ;
    } else {
        /*species 0: H2 */
        species[0] =
            +2.99142337e+00 * tc[0]
            +7.00064411e-04 * tc[1]
            -2.81691434e-08 * tc[2]
            -3.07719273e-12 * tc[3]
            +3.95687948e-16 * tc[4]
            -1.35511017e+00 ;
        /*species 1: O2 */
        species[1] =
            +3.69757819e+00 * tc[0]
            +6.13519689e-04 * tc[1]
            -6.29420995e-08 * tc[2]
            +5.91760493e-12 * tc[3]
            -2.84108828e-16 * tc[4]
            +3.18916559e+00 ;
        /*species 2: H2O */
        species[2] =
            +2.67214561e+00 * tc[0]
            +3.05629289e-03 * tc[1]
            -4.36513005e-07 * tc[2]
            +4.00332130e-11 * tc[3]
            -1.59790447e-15 * tc[4]
            +6.86281681e+00 ;
        /*species 3: H */
        species[3] =
            +2.50000000e+00 * tc[0]
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3]
            +0.00000000e+00 * tc[4]
            -4.60117638e-01 ;
        /*species 4: O */
        species[4] =
            +2.54205966e+00 * tc[0]
            -2.75506191e-05 * tc[1]
            -1.55140167e-09 * tc[2]
            +1.51702247e-12 * tc[3]
            -1.09201287e-16 * tc[4]
            +4.92030811e+00 ;
        /*species 5: OH */
        species[5] =
            +2.86472886e+00 * tc[0]
            +1.05650448e-03 * tc[1]
            -1.29541379e-07 * tc[2]
            +1.01739558e-11 * tc[3]
            -3.32989690e-16 * tc[4]
            +5.70164073e+00 ;
        /*species 6: HO2 */
        species[6] =
            +4.01721090e+00 * tc[0]
            +2.23982013e-03 * tc[1]
            -3.16829075e-07 * tc[2]
            +3.80821233e-11 * tc[3]
            -2.69771337e-15 * tc[4]
            +3.78510215e+00 ;
        /*species 7: H2O2 */
        species[7] =
            +4.57316685e+00 * tc[0]
            +4.33613639e-03 * tc[1]
            -7.37344410e-07 * tc[2]
            +7.82967857e-11 * tc[3]
            -3.57913390e-15 * tc[4]
            +5.01136959e-01 ;
        /*species 8: N2 */
        species[8] =
            +2.92664000e+00 * tc[0]
            +1.48797700e-03 * tc[1]
            -2.84238050e-07 * tc[2]
            +3.36568000e-11 * tc[3]
            -1.68833775e-15 * tc[4]
            +5.98052800e+00 ;
    }
    return;
}


/*Returns the mean specific heat at CP (Eq. 33) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCPBL(amrex::Real *  T, amrex::Real *  x,  amrex::Real *  cpbl)
{
    int id; /*loop counter */
    amrex::Real result = 0; 
    amrex::Real tT = *T; /*temporary temperature */
    const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real cpor[9]; /* temporary storage */
    cp_R(cpor, tc);

    /*perform dot product */
    for (id = 0; id < 9; ++id) {
        result += x[id]*cpor[id];
    }

    *cpbl = result * 8.31446261815324e+07;
}


/*Returns the mean specific heat at CP (Eq. 34) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCPBS(amrex::Real *  T, amrex::Real *  y,  amrex::Real *  cpbs)
{
    amrex::Real result = 0; 
    amrex::Real tT = *T; /*temporary temperature */
    const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real cpor[9], tresult[9]; /* temporary storage */
    amrex::Real imw[9];

    get_imw(imw);

    cp_R(cpor, tc);

    for (int i = 0; i < 9; i++)
    {
        tresult[i] = cpor[i]*y[i]*imw[i];

    }
    for (int i = 0; i < 9; i++)
    {
        result += tresult[i];
    }

    *cpbs = result * 8.31446261815324e+07;
}


/*Returns the mean specific heat at CV (Eq. 35) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCVBL(amrex::Real *  T, amrex::Real *  x,  amrex::Real *  cvbl)
{
    int id; /*loop counter */
    amrex::Real result = 0; 
    amrex::Real tT = *T; /*temporary temperature */
    const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real cvor[9]; /* temporary storage */
    cv_R(cvor, tc);

    /*perform dot product */
    for (id = 0; id < 9; ++id) {
        result += x[id]*cvor[id];
    }

    *cvbl = result * 8.31446261815324e+07;
}


/*Returns the mean specific heat at CV (Eq. 36) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCVBS(amrex::Real *  T, amrex::Real *  y,  amrex::Real *  cvbs)
{
    amrex::Real result = 0; 
    amrex::Real tT = *T; /*temporary temperature */
    const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real cvor[9]; /* temporary storage */
    amrex::Real imw[9];

    get_imw(imw);

    cv_R(cvor, tc);

    /*multiply by y/molecularweight */
    result += cvor[0]*y[0]*imw[0]; /*H2 */
    result += cvor[1]*y[1]*imw[1]; /*O2 */
    result += cvor[2]*y[2]*imw[2]; /*H2O */
    result += cvor[3]*y[3]*imw[3]; /*H */
    result += cvor[4]*y[4]*imw[4]; /*O */
    result += cvor[5]*y[5]*imw[5]; /*OH */
    result += cvor[6]*y[6]*imw[6]; /*HO2 */
    result += cvor[7]*y[7]*imw[7]; /*H2O2 */
    result += cvor[8]*y[8]*imw[8]; /*N2 */

    *cvbs = result * 8.31446261815324e+07;
}


/*Returns the mean enthalpy of the mixture in molar units */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKHBML(amrex::Real *  T, amrex::Real *  x,  amrex::Real *  hbml)
{
    int id; /*loop counter */
    amrex::Real result = 0; 
    amrex::Real tT = *T; /*temporary temperature */
    const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real hml[9]; /* temporary storage */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*perform dot product */
    for (id = 0; id < 9; ++id) {
        result += x[id]*hml[id];
    }

    *hbml = result * RT;
}


/*Returns mean enthalpy of mixture in mass units */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKHBMS(amrex::Real *  T, amrex::Real *  y,  amrex::Real *  hbms)
{
    amrex::Real result = 0;
    amrex::Real tT = *T; /*temporary temperature */
    const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real hml[9], tmp[9]; /* temporary storage */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    amrex::Real imw[9];

    get_imw(imw);

    speciesEnthalpy(hml, tc);

    int id;
    for (id = 0; id < 9; ++id) {
        tmp[id] = y[id]*hml[id]*imw[id];
    }
    for (id = 0; id < 9; ++id) {
        result += tmp[id];
    }

    *hbms = result * RT;
}


/*get mean internal energy in molar units */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKUBML(amrex::Real *  T, amrex::Real *  x,  amrex::Real *  ubml)
{
    int id; /*loop counter */
    amrex::Real result = 0; 
    amrex::Real tT = *T; /*temporary temperature */
    const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real uml[9]; /* temporary energy array */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*perform dot product */
    for (id = 0; id < 9; ++id) {
        result += x[id]*uml[id];
    }

    *ubml = result * RT;
}


/*get mean internal energy in mass units */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKUBMS(amrex::Real *  T, amrex::Real *  y,  amrex::Real *  ubms)
{
    amrex::Real result = 0;
    amrex::Real tT = *T; /*temporary temperature */
    const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real ums[9]; /* temporary energy array */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    amrex::Real imw[9];

    get_imw(imw);

    speciesInternalEnergy(ums, tc);

    /*perform dot product + scaling by wt */
    result += y[0]*ums[0]*imw[0]; /*H2 */
    result += y[1]*ums[1]*imw[1]; /*O2 */
    result += y[2]*ums[2]*imw[2]; /*H2O */
    result += y[3]*ums[3]*imw[3]; /*H */
    result += y[4]*ums[4]*imw[4]; /*O */
    result += y[5]*ums[5]*imw[5]; /*OH */
    result += y[6]*ums[6]*imw[6]; /*HO2 */
    result += y[7]*ums[7]*imw[7]; /*H2O2 */
    result += y[8]*ums[8]*imw[8]; /*N2 */

    *ubms = result * RT;
}


/*get mixture entropy in molar units */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKSBML(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x,  amrex::Real *  sbml)
{
    int id; /*loop counter */
    amrex::Real result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    amrex::Real logPratio = log ( *P / 1013250.0 ); 
    amrex::Real tT = *T; /*temporary temperature */
    const amrex::Real tc[5] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real sor[9]; /* temporary storage */
    speciesEntropy(sor, tc);

    /*Compute Eq 42 */
    for (id = 0; id < 9; ++id) {
        result += x[id]*(sor[id]-log((x[id]+1e-100))-logPratio);
    }

    *sbml = result * 8.31446261815324e+07;
}


/*get mixture entropy in mass units */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKSBMS(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  sbms)
{
    amrex::Real result = 0; 
    /*Log of normalized pressure in cgs units dynes/cm^2 by Patm */
    amrex::Real logPratio = log ( *P / 1013250.0 ); 
    amrex::Real tT = *T; /*temporary temperature */
    const amrex::Real tc[5] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real sor[9]; /* temporary storage */
    amrex::Real x[9]; /* need a ytx conversion */
    amrex::Real YOW = 0; /*See Eq 4, 6 in CK Manual */
    amrex::Real imw[9];

    get_imw(imw);

    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*H2O */
    YOW += y[3]*imw[3]; /*H */
    YOW += y[4]*imw[4]; /*O */
    YOW += y[5]*imw[5]; /*OH */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(31.998800*YOW); 
    x[2] = y[2]/(18.015340*YOW); 
    x[3] = y[3]/(1.007970*YOW); 
    x[4] = y[4]/(15.999400*YOW); 
    x[5] = y[5]/(17.007370*YOW); 
    x[6] = y[6]/(33.006770*YOW); 
    x[7] = y[7]/(34.014740*YOW); 
    x[8] = y[8]/(28.013400*YOW); 
    speciesEntropy(sor, tc);
    /*Perform computation in Eq 42 and 43 */
    result += x[0]*(sor[0]-log((x[0]+1e-100))-logPratio);
    result += x[1]*(sor[1]-log((x[1]+1e-100))-logPratio);
    result += x[2]*(sor[2]-log((x[2]+1e-100))-logPratio);
    result += x[3]*(sor[3]-log((x[3]+1e-100))-logPratio);
    result += x[4]*(sor[4]-log((x[4]+1e-100))-logPratio);
    result += x[5]*(sor[5]-log((x[5]+1e-100))-logPratio);
    result += x[6]*(sor[6]-log((x[6]+1e-100))-logPratio);
    result += x[7]*(sor[7]-log((x[7]+1e-100))-logPratio);
    result += x[8]*(sor[8]-log((x[8]+1e-100))-logPratio);
    /*Scale by R/W */
    *sbms = result * 8.31446261815324e+07 * YOW;
}

/* get temperature given internal energy in mass units and mass fracs */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void GET_T_GIVEN_EY(amrex::Real *  e, amrex::Real *  y, amrex::Real *  t, int * ierr)
{
#ifdef CONVERGENCE
    const int maxiter = 5000;
    const amrex::Real tol  = 1.e-12;
#else
    const int maxiter = 200;
    const amrex::Real tol  = 1.e-6;
#endif
    amrex::Real ein  = *e;
    amrex::Real tmin = 90;/*max lower bound for thermo def */
    amrex::Real tmax = 4000;/*min upper bound for thermo def */
    amrex::Real e1,emin,emax,cv,t1,dt;
    int i;/* loop counter */
    CKUBMS(&tmin, y, &emin);
    CKUBMS(&tmax, y, &emax);
    if (ein < emin) {
        /*Linear Extrapolation below tmin */
        CKCVBS(&tmin, y, &cv);
        *t = tmin - (emin-ein)/cv;
        *ierr = 1;
        return;
    }
    if (ein > emax) {
        /*Linear Extrapolation above tmax */
        CKCVBS(&tmax, y, &cv);
        *t = tmax - (emax-ein)/cv;
        *ierr = 1;
        return;
    }
    t1 = *t;
    if (t1 < tmin || t1 > tmax) {
        t1 = tmin + (tmax-tmin)/(emax-emin)*(ein-emin);
    }
    for (i = 0; i < maxiter; ++i) {
        CKUBMS(&t1,y,&e1);
        CKCVBS(&t1,y,&cv);
        dt = (ein - e1) / cv;
        if (dt > 100.) { dt = 100.; }
        else if (dt < -100.) { dt = -100.; }
        else if (fabs(dt) < tol) break;
        else if (t1+dt == t1) break;
        t1 += dt;
    }
    *t = t1;
    *ierr = 0;
    return;
}

/* get temperature given enthalpy in mass units and mass fracs */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void GET_T_GIVEN_HY(amrex::Real *  h, amrex::Real *  y, amrex::Real *  t, int * ierr)
{
#ifdef CONVERGENCE
    const int maxiter = 5000;
    const amrex::Real tol  = 1.e-12;
#else
    const int maxiter = 200;
    const amrex::Real tol  = 1.e-6;
#endif
    amrex::Real hin  = *h;
    amrex::Real tmin = 90;/*max lower bound for thermo def */
    amrex::Real tmax = 4000;/*min upper bound for thermo def */
    amrex::Real h1,hmin,hmax,cp,t1,dt;
    int i;/* loop counter */
    CKHBMS(&tmin, y, &hmin);
    CKHBMS(&tmax, y, &hmax);
    if (hin < hmin) {
        /*Linear Extrapolation below tmin */
        CKCPBS(&tmin, y, &cp);
        *t = tmin - (hmin-hin)/cp;
        *ierr = 1;
        return;
    }
    if (hin > hmax) {
        /*Linear Extrapolation above tmax */
        CKCPBS(&tmax, y, &cp);
        *t = tmax - (hmax-hin)/cp;
        *ierr = 1;
        return;
    }
    t1 = *t;
    if (t1 < tmin || t1 > tmax) {
        t1 = tmin + (tmax-tmin)/(hmax-hmin)*(hin-hmin);
    }
    for (i = 0; i < maxiter; ++i) {
        CKHBMS(&t1,y,&h1);
        CKCPBS(&t1,y,&cp);
        dt = (hin - h1) / cp;
        if (dt > 100.) { dt = 100.; }
        else if (dt < -100.) { dt = -100.; }
        else if (fabs(dt) < tol) break;
        else if (t1+dt == t1) break;
        t1 += dt;
    }
    *t = t1;
    *ierr = 0;
    return;
}


/*Compute P = rhoRT/W(x) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKPX(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  x, amrex::Real *  P)
{
    amrex::Real XW = 0;/* To hold mean molecular wt */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*18.015340; /*H2O */
    XW += x[3]*1.007970; /*H */
    XW += x[4]*15.999400; /*O */
    XW += x[5]*17.007370; /*OH */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    *P = *rho * 8.31446261815324e+07 * (*T) / XW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(y) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKPY(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  P)
{
    amrex::Real YOW = 0;/* for computing mean MW */
    amrex::Real imw[9];

    get_imw(imw);

    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*H2O */
    YOW += y[3]*imw[3]; /*H */
    YOW += y[4]*imw[4]; /*O */
    YOW += y[5]*imw[5]; /*OH */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*N2 */
    *P = *rho * 8.31446261815324e+07 * (*T) * YOW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(c) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKPC(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  c,  amrex::Real *  P)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    amrex::Real W = 0;
    amrex::Real sumC = 0;
    W += c[0]*2.015940; /*H2 */
    W += c[1]*31.998800; /*O2 */
    W += c[2]*18.015340; /*H2O */
    W += c[3]*1.007970; /*H */
    W += c[4]*15.999400; /*O */
    W += c[5]*17.007370; /*OH */
    W += c[6]*33.006770; /*HO2 */
    W += c[7]*34.014740; /*H2O2 */
    W += c[8]*28.013400; /*N2 */

    for (id = 0; id < 9; ++id) {
        sumC += c[id];
    }
    *P = *rho * 8.31446261815324e+07 * (*T) * sumC / W; /*P = rho*R*T/W */

    return;
}


/*Compute rho = PW(x)/RT */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKRHOX(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x,  amrex::Real *  rho)
{
    amrex::Real XW = 0;/* To hold mean molecular wt */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*18.015340; /*H2O */
    XW += x[3]*1.007970; /*H */
    XW += x[4]*15.999400; /*O */
    XW += x[5]*17.007370; /*OH */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    *rho = *P * XW / (8.31446261815324e+07 * (*T)); /*rho = P*W/(R*T) */

    return;
}


/*Compute rho = P*W(y)/RT */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKRHOY(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  rho)
{
    amrex::Real YOW = 0;
    amrex::Real tmp[9];
    amrex::Real imw[9];

    get_imw(imw);

    for (int i = 0; i < 9; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 9; i++)
    {
        YOW += tmp[i];
    }

    *rho = *P / (8.31446261815324e+07 * (*T) * YOW);/*rho = P*W/(R*T) */
    return;
}


/*Compute rho = P*W(c)/(R*T) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKRHOC(amrex::Real *  P, amrex::Real *  T, amrex::Real *  c,  amrex::Real *  rho)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    amrex::Real W = 0;
    amrex::Real sumC = 0;
    W += c[0]*2.015940; /*H2 */
    W += c[1]*31.998800; /*O2 */
    W += c[2]*18.015340; /*H2O */
    W += c[3]*1.007970; /*H */
    W += c[4]*15.999400; /*O */
    W += c[5]*17.007370; /*OH */
    W += c[6]*33.006770; /*HO2 */
    W += c[7]*34.014740; /*H2O2 */
    W += c[8]*28.013400; /*N2 */

    for (id = 0; id < 9; ++id) {
        sumC += c[id];
    }
    *rho = *P * W / (sumC * (*T) * 8.31446261815324e+07); /*rho = PW/(R*T) */

    return;
}


/*get molecular weight for all species */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWT( amrex::Real *  wt)
{
    get_mw(wt);
}

/*given y[species]: mass fractions */
/*returns mean molecular weight (gm/mole) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKMMWY(amrex::Real *  y,  amrex::Real *  wtm)
{
    amrex::Real YOW = 0;
    amrex::Real tmp[9];
    amrex::Real imw[9];

    get_imw(imw);

    for (int i = 0; i < 9; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 9; i++)
    {
        YOW += tmp[i];
    }

    *wtm = 1.0 / YOW;
    return;
}


/*given x[species]: mole fractions */
/*returns mean molecular weight (gm/mole) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKMMWX(amrex::Real *  x,  amrex::Real *  wtm)
{
    amrex::Real XW = 0;/* see Eq 4 in CK Manual */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*18.015340; /*H2O */
    XW += x[3]*1.007970; /*H */
    XW += x[4]*15.999400; /*O */
    XW += x[5]*17.007370; /*OH */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    *wtm = XW;

    return;
}

/*given c[species]: molar concentration */
/*returns mean molecular weight (gm/mole) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKMMWC(amrex::Real *  c,  amrex::Real *  wtm)
{
    int id; /*loop counter */
    /*See Eq 5 in CK Manual */
    amrex::Real W = 0;
    amrex::Real sumC = 0;
    W += c[0]*2.015940; /*H2 */
    W += c[1]*31.998800; /*O2 */
    W += c[2]*18.015340; /*H2O */
    W += c[3]*1.007970; /*H */
    W += c[4]*15.999400; /*O */
    W += c[5]*17.007370; /*OH */
    W += c[6]*33.006770; /*HO2 */
    W += c[7]*34.014740; /*H2O2 */
    W += c[8]*28.013400; /*N2 */

    for (id = 0; id < 9; ++id) {
        sumC += c[id];
    }
    /* CK provides no guard against divison by zero */
    *wtm = W/sumC;

    return;
}


/*get Cp/R as a function of T  */
/*for all species (Eq 19) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCPOR(amrex::Real *  T, amrex::Real *  cpor)
{
    amrex::Real tT = *T; /*temporary temperature */
    const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpor, tc);
}


/*get H/RT as a function of T  */
/*for all species (Eq 20) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKHORT(amrex::Real *  T, amrex::Real *  hort)
{
    amrex::Real tT = *T; /*temporary temperature */
    const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEnthalpy(hort, tc);
}


/*get S/R as a function of T  */
/*for all species (Eq 21) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKSOR(amrex::Real *  T, amrex::Real *  sor)
{
    amrex::Real tT = *T; /*temporary temperature */
    const amrex::Real tc[5] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sor, tc);
}

/*convert y[species] (mass fracs) to x[species] (mole fracs) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKYTX(amrex::Real *  y,  amrex::Real *  x)
{
    amrex::Real YOW = 0;
    amrex::Real tmp[9];
    amrex::Real imw[9];

    get_imw(imw);

    for (int i = 0; i < 9; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 9; i++)
    {
        YOW += tmp[i];
    }

    amrex::Real YOWINV = 1.0/YOW;

    for (int i = 0; i < 9; i++)
    {
        x[i] = y[i]*imw[i]*YOWINV;
    }
    return;
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKYTCP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  c)
{
    amrex::Real YOW = 0;
    amrex::Real PWORT;
    amrex::Real imw[9];

    get_imw(imw);

    /*Compute inverse of mean molecular wt first */
    for (int i = 0; i < 9; i++)
    {
        c[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 9; i++)
    {
        YOW += c[i];
    }

    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31446261815324e+07 * (*T)); 
    /*Now compute conversion */

    for (int i = 0; i < 9; i++)
    {
        c[i] = PWORT * y[i] * imw[i];
    }
    return;
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKYTCR(amrex::Real *  rho, amrex::Real * /*T*/, amrex::Real * y,  amrex::Real * c)
{
    amrex::Real imw[9];

    get_imw(imw);

    for (int i = 0; i < 9; i++)
    {
        c[i] = (*rho)  * y[i] * imw[i];
    }
}


/*convert x[species] (mole fracs) to y[species] (mass fracs) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKXTY(amrex::Real *  x,  amrex::Real *  y)
{
    amrex::Real XW = 0; /*See Eq 4, 9 in CK Manual */
    /*Compute mean molecular wt first */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*18.015340; /*H2O */
    XW += x[3]*1.007970; /*H */
    XW += x[4]*15.999400; /*O */
    XW += x[5]*17.007370; /*OH */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    /*Now compute conversion */
    amrex::Real XWinv = 1.0/XW;
    y[0] = x[0]*2.015940*XWinv; 
    y[1] = x[1]*31.998800*XWinv; 
    y[2] = x[2]*18.015340*XWinv; 
    y[3] = x[3]*1.007970*XWinv; 
    y[4] = x[4]*15.999400*XWinv; 
    y[5] = x[5]*17.007370*XWinv; 
    y[6] = x[6]*33.006770*XWinv; 
    y[7] = x[7]*34.014740*XWinv; 
    y[8] = x[8]*28.013400*XWinv; 

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKXTCP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x,  amrex::Real *  c)
{
    int id; /*loop counter */
    amrex::Real PORT = (*P)/(8.31446261815324e+07 * (*T)); /*P/RT */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 9; ++id) {
        c[id] = x[id]*PORT;
    }

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKXTCR(amrex::Real *  rho, amrex::Real * /*T*/, amrex::Real *  x, amrex::Real *  c)
{
    int id; /*loop counter */
    amrex::Real XW = 0; /*See Eq 4, 11 in CK Manual */
    amrex::Real ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*18.015340; /*H2O */
    XW += x[3]*1.007970; /*H */
    XW += x[4]*15.999400; /*O */
    XW += x[5]*17.007370; /*OH */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    ROW = (*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 9; ++id) {
        c[id] = x[id]*ROW;
    }

    return;
}


/*convert c[species] (molar conc) to x[species] (mole fracs) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCTX(amrex::Real *  c, amrex::Real *  x)
{
    int id; /*loop counter */
    amrex::Real sumC = 0; 

    /*compute sum of c  */
    for (id = 0; id < 9; ++id) {
        sumC += c[id];
    }

    /* See Eq 13  */
    amrex::Real sumCinv = 1.0/sumC;
    for (id = 0; id < 9; ++id) {
        x[id] = c[id]*sumCinv;
    }

    return;
}


/*convert c[species] (molar conc) to y[species] (mass fracs) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCTY(amrex::Real *  c, amrex::Real *  y)
{
    amrex::Real CW = 0; /*See Eq 12 in CK Manual */
    /*compute denominator in eq 12 first */
    CW += c[0]*2.015940; /*H2 */
    CW += c[1]*31.998800; /*O2 */
    CW += c[2]*18.015340; /*H2O */
    CW += c[3]*1.007970; /*H */
    CW += c[4]*15.999400; /*O */
    CW += c[5]*17.007370; /*OH */
    CW += c[6]*33.006770; /*HO2 */
    CW += c[7]*34.014740; /*H2O2 */
    CW += c[8]*28.013400; /*N2 */
    /*Now compute conversion */
    amrex::Real CWinv = 1.0/CW;
    y[0] = c[0]*2.015940*CWinv; 
    y[1] = c[1]*31.998800*CWinv; 
    y[2] = c[2]*18.015340*CWinv; 
    y[3] = c[3]*1.007970*CWinv; 
    y[4] = c[4]*15.999400*CWinv; 
    y[5] = c[5]*17.007370*CWinv; 
    y[6] = c[6]*33.006770*CWinv; 
    y[7] = c[7]*34.014740*CWinv; 
    y[8] = c[8]*28.013400*CWinv; 

    return;
}


/*get specific heat at constant volume as a function  */
/*of T for all species (molar units) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCVML(amrex::Real *  T,  amrex::Real *  cvml)
{
    int id; /*loop counter */
    amrex::Real tT = *T; /*temporary temperature */
    const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cv_R(cvml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        cvml[id] *= 8.31446261815324e+07;
    }
}


/*get specific heat at constant pressure as a  */
/*function of T for all species (molar units) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCPML(amrex::Real *  T,  amrex::Real *  cpml)
{
    int id; /*loop counter */
    amrex::Real tT = *T; /*temporary temperature */
    const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        cpml[id] *= 8.31446261815324e+07;
    }
}


/*get internal energy as a function  */
/*of T for all species (molar units) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKUML(amrex::Real *  T,  amrex::Real *  uml)
{
    int id; /*loop counter */
    amrex::Real tT = *T; /*temporary temperature */
    const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        uml[id] *= RT;
    }
}


/*get enthalpy as a function  */
/*of T for all species (molar units) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKHML(amrex::Real *  T,  amrex::Real *  hml)
{
    int id; /*loop counter */
    amrex::Real tT = *T; /*temporary temperature */
    const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        hml[id] *= RT;
    }
}


/*Returns the standard-state entropies in molar units */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKSML(amrex::Real *  T,  amrex::Real *  sml)
{
    int id; /*loop counter */
    amrex::Real tT = *T; /*temporary temperature */
    const amrex::Real tc[5] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sml, tc);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        sml[id] *= 8.31446261815324e+07;
    }
}


/*Returns the specific heats at constant volume */
/*in mass units (Eq. 29) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCVMS(amrex::Real *  T,  amrex::Real *  cvms)
{
    amrex::Real tT = *T; /*temporary temperature */
    const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cv_R(cvms, tc);
    /*multiply by R/molecularweight */
    cvms[0] *= 4.124360158612479e+07; /*H2 */
    cvms[1] *= 2.598367006935648e+06; /*O2 */
    cvms[2] *= 4.615212712140454e+06; /*H2O */
    cvms[3] *= 8.248720317224957e+07; /*H */
    cvms[4] *= 5.196734013871295e+06; /*O */
    cvms[5] *= 4.888740950630956e+06; /*OH */
    cvms[6] *= 2.519017346487778e+06; /*HO2 */
    cvms[7] *= 2.444370475315478e+06; /*H2O2 */
    cvms[8] *= 2.968030520448514e+06; /*N2 */
}


/*Returns the specific heats at constant pressure */
/*in mass units (Eq. 26) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCPMS(amrex::Real *  T,  amrex::Real *  cpms)
{
    amrex::Real tT = *T; /*temporary temperature */
    const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpms, tc);
    /*multiply by R/molecularweight */
    cpms[0] *= 4.124360158612479e+07; /*H2 */
    cpms[1] *= 2.598367006935648e+06; /*O2 */
    cpms[2] *= 4.615212712140454e+06; /*H2O */
    cpms[3] *= 8.248720317224957e+07; /*H */
    cpms[4] *= 5.196734013871295e+06; /*O */
    cpms[5] *= 4.888740950630956e+06; /*OH */
    cpms[6] *= 2.519017346487778e+06; /*HO2 */
    cpms[7] *= 2.444370475315478e+06; /*H2O2 */
    cpms[8] *= 2.968030520448514e+06; /*N2 */
}


/*Returns internal energy in mass units (Eq 30.) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKUMS(amrex::Real *  T,  amrex::Real *  ums)
{
    amrex::Real tT = *T; /*temporary temperature */
    const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    amrex::Real imw[9];

    get_imw(imw);

    speciesInternalEnergy(ums, tc);

    for (int i = 0; i < 9; i++)
    {
        ums[i] *= RT*imw[i];
    }
}


/*Returns enthalpy in mass units (Eq 27.) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKHMS(amrex::Real *  T,  amrex::Real *  hms)
{
    amrex::Real tT = *T; /*temporary temperature */
    const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    amrex::Real imw[9];

    get_imw(imw);

    speciesEnthalpy(hms, tc);

    for (int i = 0; i < 9; i++)
    {
        hms[i] *= RT*imw[i];
    }
}


/*Returns the entropies in mass units (Eq 28.) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKSMS(amrex::Real *  T,  amrex::Real *  sms)
{
    amrex::Real tT = *T; /*temporary temperature */
    const amrex::Real tc[5] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sms, tc);
    /*multiply by R/molecularweight */
    sms[0] *= 4.124360158612479e+07; /*H2 */
    sms[1] *= 2.598367006935648e+06; /*O2 */
    sms[2] *= 4.615212712140454e+06; /*H2O */
    sms[3] *= 8.248720317224957e+07; /*H */
    sms[4] *= 5.196734013871295e+06; /*O */
    sms[5] *= 4.888740950630956e+06; /*OH */
    sms[6] *= 2.519017346487778e+06; /*HO2 */
    sms[7] *= 2.444370475315478e+06; /*H2O2 */
    sms[8] *= 2.968030520448514e+06; /*N2 */
}


/* GPU version of productionRate: no more use of thermo namespace vectors  */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void comp_qfqr(amrex::Real *  qf, amrex::Real * qr, amrex::Real * sc, amrex::Real * sc_qss, amrex::Real * tc, amrex::Real invT)
{

    /*reaction 9: H + O2 (+M) <=> HO2 (+M) */
    qf[0] = sc[1]*sc[3];
    qr[0] = sc[6];

    /*reaction 16: H2O2 (+M) <=> OH + OH (+M) */
    qf[1] = sc[7];
    qr[1] = sc[5]*sc[5];

    /*reaction 5: H2 + M <=> H + H + M */
    qf[2] = sc[0];
    qr[2] = sc[3]*sc[3];

    /*reaction 6: O + O + M <=> O2 + M */
    qf[3] = sc[4]*sc[4];
    qr[3] = sc[1];

    /*reaction 7: O + H + M <=> OH + M */
    qf[4] = sc[3]*sc[4];
    qr[4] = sc[5];

    /*reaction 8: H + OH + M <=> H2O + M */
    qf[5] = sc[3]*sc[5];
    qr[5] = sc[2];

    /*reaction 1: H + O2 <=> O + OH */
    qf[6] = sc[1]*sc[3];
    qr[6] = sc[4]*sc[5];

    /*reaction 2: O + H2 <=> H + OH */
    qf[7] = sc[0]*sc[4];
    qr[7] = sc[3]*sc[5];

    /*reaction 3: H2 + OH <=> H2O + H */
    qf[8] = sc[0]*sc[5];
    qr[8] = sc[2]*sc[3];

    /*reaction 4: O + H2O <=> OH + OH */
    qf[9] = sc[2]*sc[4];
    qr[9] = sc[5]*sc[5];

    /*reaction 10: HO2 + H <=> H2 + O2 */
    qf[10] = sc[3]*sc[6];
    qr[10] = sc[0]*sc[1];

    /*reaction 11: HO2 + H <=> OH + OH */
    qf[11] = sc[3]*sc[6];
    qr[11] = sc[5]*sc[5];

    /*reaction 12: HO2 + O <=> O2 + OH */
    qf[12] = sc[4]*sc[6];
    qr[12] = sc[1]*sc[5];

    /*reaction 13: HO2 + OH <=> H2O + O2 */
    qf[13] = sc[5]*sc[6];
    qr[13] = sc[1]*sc[2];

    /*reaction 14: HO2 + HO2 <=> H2O2 + O2 */
    qf[14] = sc[6]*sc[6];
    qr[14] = sc[1]*sc[7];

    /*reaction 15: HO2 + HO2 <=> H2O2 + O2 */
    qf[15] = sc[6]*sc[6];
    qr[15] = sc[1]*sc[7];

    /*reaction 17: H2O2 + H <=> H2O + OH */
    qf[16] = sc[3]*sc[7];
    qr[16] = sc[2]*sc[5];

    /*reaction 18: H2O2 + H <=> HO2 + H2 */
    qf[17] = sc[3]*sc[7];
    qr[17] = sc[0]*sc[6];

    /*reaction 19: H2O2 + O <=> OH + HO2 */
    qf[18] = sc[4]*sc[7];
    qr[18] = sc[5]*sc[6];

    /*reaction 20: H2O2 + OH <=> HO2 + H2O */
    qf[19] = sc[5]*sc[7];
    qr[19] = sc[2]*sc[6];

    /*reaction 21: H2O2 + OH <=> HO2 + H2O */
    qf[20] = sc[5]*sc[7];
    qr[20] = sc[2]*sc[6];

    /*compute the mixture concentration */
    amrex::Real mixture = 0.0;
    for (int i = 0; i < 9; ++i) {
        mixture += sc[i];
    }

    /*compute the Gibbs free energy */
    amrex::Real g_RT[9];
    gibbs(g_RT, tc);

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    amrex::Real refC = 101325 / 8.31446 * invT;
    amrex::Real refCinv = 1 / refC;

    // Evaluate the kfs
    amrex::Real k_f, k_r, Corr;
    amrex::Real redP, F, logPred, logFcent, troe_c, troe_n, troe, F_troe;

    // (0):  H + O2 <=> O + OH
    k_f = 3547000000
               * exp((-0.406) * tc[0] - (8352.89343569254) * invT);
    qf[6] *= k_f;
    qr[6] *= k_f * exp(-(g_RT[1] + g_RT[3] - g_RT[4] - g_RT[5]));
    // (1):  O + H2 <=> H + OH
    k_f = 0.0508
               * exp((2.67) * tc[0] - (3165.23282791169) * invT);
    qf[7] *= k_f;
    qr[7] *= k_f * exp(-(g_RT[0] - g_RT[3] + g_RT[4] - g_RT[5]));
    // (2):  H2 + OH <=> H2O + H
    k_f = 216
               * exp((1.51) * tc[0] - (1726.03316371019) * invT);
    qf[8] *= k_f;
    qr[8] *= k_f * exp(-(g_RT[0] - g_RT[2] - g_RT[3] + g_RT[5]));
    // (3):  O + H2O <=> OH + OH
    k_f = 2.97
               * exp((2.02) * tc[0] - (6743.10332178324) * invT);
    qf[9] *= k_f;
    qr[9] *= k_f * exp(-(g_RT[2] + g_RT[4] - g_RT[5] - g_RT[5]));
    // (4):  H2 + M <=> H + H + M
    k_f = 45770000000000
               * exp((-1.4) * tc[0] - (52525.7555766966) * invT);
    Corr  = mixture + (1.5)*sc[0] + (11)*sc[2];
    qf[2] *= Corr * k_f;
    qr[2] *= Corr * k_f * exp(-(g_RT[0] - g_RT[3] - g_RT[3])) * (refCinv);
    // (5):  O + O + M <=> O2 + M
    k_f = 6165
               * exp((-0.5) * tc[0]);
    Corr  = mixture + (1.5)*sc[0] + (11)*sc[2];
    qf[3] *= Corr * k_f;
    qr[3] *= Corr * k_f * exp(-(-g_RT[1] + g_RT[4] + g_RT[4])) * (refC);
    // (6):  O + H + M <=> OH + M
    k_f = 4714000
               * exp((-1) * tc[0]);
    Corr  = mixture + (1.5)*sc[0] + (11)*sc[2];
    qf[4] *= Corr * k_f;
    qr[4] *= Corr * k_f * exp(-(g_RT[3] + g_RT[4] - g_RT[5])) * (refC);
    // (7):  H + OH + M <=> H2O + M
    k_f = 38000000000
               * exp((-2) * tc[0]);
    Corr  = mixture + (1.5)*sc[0] + (11)*sc[2];
    qf[5] *= Corr * k_f;
    qr[5] *= Corr * k_f * exp(-(-g_RT[2] + g_RT[3] + g_RT[5])) * (refC);
    // (8):  H + O2 (+M) <=> HO2 (+M)
    k_f = 1475000
               * exp((0.6) * tc[0]);
    Corr  = mixture + sc[0] + (10)*sc[2] + (-0.22)*sc[1];
    redP = Corr / k_f * 1e-12 * 6.366e+20 
               * exp(-1.72  * tc[0] - 0.50321666580472  * (524.8) *invT);
    F = redP / (1.0 + redP);
    logPred = log10(redP);
    logFcent = log10(
        (0.2)*exp(-tc[1] * 1e+30)
        + 0.8 * exp(-tc[1] * 1e-30)
        + 0.0);
    troe_c = -0.4 - 0.67 * logFcent;
    troe_n = 0.75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - 0.14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    Corr = F * F_troe;
    qf[0] *= Corr * k_f;
    qr[0] *= Corr * k_f * exp(-(g_RT[1] + g_RT[3] - g_RT[6])) * (refC);
    // (9):  HO2 + H <=> H2 + O2
    k_f = 16600000
               * exp(-(414.147315957284) * invT);
    qf[10] *= k_f;
    qr[10] *= k_f * exp(-(-g_RT[0] - g_RT[1] + g_RT[3] + g_RT[6]));
    // (10):  HO2 + H <=> OH + OH
    k_f = 70790000
               * exp(-(148.448916412392) * invT);
    qf[11] *= k_f;
    qr[11] *= k_f * exp(-(g_RT[3] - g_RT[5] - g_RT[5] + g_RT[6]));
    // (11):  HO2 + O <=> O2 + OH
    k_f = 32500000
               ;
    qf[12] *= k_f;
    qr[12] *= k_f * exp(-(-g_RT[1] + g_RT[4] - g_RT[5] + g_RT[6]));
    // (12):  HO2 + OH <=> H2O + O2
    k_f = 28900000
               * exp(-(-250.098682904946) * invT);
    qf[13] *= k_f;
    qr[13] *= k_f * exp(-(-g_RT[1] - g_RT[2] + g_RT[5] + g_RT[6]));
    // (13):  HO2 + HO2 <=> H2O2 + O2
    k_f = 420000000
               * exp(-(6029.54208967215) * invT);
    qf[14] *= k_f;
    qr[14] *= k_f * exp(-(-g_RT[1] + g_RT[6] + g_RT[6] - g_RT[7]));
    // (14):  HO2 + HO2 <=> H2O2 + O2
    k_f = 130000
               * exp(-(-819.89091359563) * invT);
    qf[15] *= k_f;
    qr[15] *= k_f * exp(-(-g_RT[1] + g_RT[6] + g_RT[6] - g_RT[7]));
    // (15):  H2O2 (+M) <=> OH + OH (+M)
    k_f = 295100000000000
               * exp(-(24370.7831249226) * invT);
    Corr  = mixture + (1.5)*sc[0] + (11)*sc[2];
    redP = Corr / k_f * 1e-6 * 1.202e+17 
               * exp(0  * tc[0] - 0.50321666580472  * (45500) *invT);
    F = redP / (1.0 + redP);
    logPred = log10(redP);
    logFcent = log10(
        (0.5)*exp(-tc[1] * 1e+30)
        + 0.5 * exp(-tc[1] * 1e-30)
        + 0.0);
    troe_c = -0.4 - 0.67 * logFcent;
    troe_n = 0.75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - 0.14*(troe_c + logPred));
    F_troe = pow(10, logFcent / (1.0 + troe*troe));
    Corr = F * F_troe;
    qf[1] *= Corr * k_f;
    qr[1] *= Corr * k_f * exp(-(-g_RT[5] - g_RT[5] + g_RT[7])) * (refCinv);
    // (16):  H2O2 + H <=> H2O + OH
    k_f = 24100000
               * exp(-(1997.77016324474) * invT);
    qf[16] *= k_f;
    qr[16] *= k_f * exp(-(-g_RT[2] + g_RT[3] - g_RT[5] + g_RT[7]));
    // (17):  H2O2 + H <=> HO2 + H2
    k_f = 48200000
               * exp(-(4000.57249314752) * invT);
    qf[17] *= k_f;
    qr[17] *= k_f * exp(-(-g_RT[0] + g_RT[3] - g_RT[6] + g_RT[7]));
    // (18):  H2O2 + O <=> OH + HO2
    k_f = 9.55
               * exp((2) * tc[0] - (1997.77016324474) * invT);
    qf[18] *= k_f;
    qr[18] *= k_f * exp(-(g_RT[4] - g_RT[5] - g_RT[6] + g_RT[7]));
    // (19):  H2O2 + OH <=> HO2 + H2O
    k_f = 1000000
               ;
    qf[19] *= k_f;
    qr[19] *= k_f * exp(-(-g_RT[2] + g_RT[5] - g_RT[6] + g_RT[7]));
    // (20):  H2O2 + OH <=> HO2 + H2O
    k_f = 580000000
               * exp(-(4809.24167509571) * invT);
    qf[20] *= k_f;
    qr[20] *= k_f * exp(-(-g_RT[2] + g_RT[5] - g_RT[6] + g_RT[7]));


    return;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void  productionRate(amrex::Real * wdot, amrex::Real * sc, amrex::Real T)
{
    const amrex::Real tc[5] = { log(T), T, T*T, T*T*T, T*T*T*T }; // temperature cache
    const amrex::Real invT = 1.0 / tc[1];

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    const amrex::Real refC = 101325 / 8.31446 * invT;
    const amrex::Real refCinv = 1 / refC;

    for (int i = 0; i < 9; ++i) {
        wdot[i] = 0.0;
    }

    /*compute the mixture concentration */
    amrex::Real mixture = 0.0;
    for (int i = 0; i < 9; ++i) {
        mixture += sc[i];
    }

    /*compute the Gibbs free energy */
    amrex::Real g_RT[9];
    gibbs(g_RT, tc);

    {
        // (8):  H + O2 (+M) <=> HO2 (+M)
        const amrex::Real k_f = 1475000
                   * exp((0.6) * tc[0]);
        amrex::Real Corr = mixture + sc[0] + (10)*sc[2] + (-0.22)*sc[1];
        const amrex::Real redP = Corr / k_f * 636600000 
                   * exp(-1.72 * tc[0] - 264.088106214317 * invT);
        const amrex::Real F = redP / (1.0 + redP);
        const amrex::Real logPred = log10(redP);
        const amrex::Real logFcent = log10(
            0.2 * exp(-tc[1] * 1e+30)
            + 0.8 * exp(-tc[1] * 1e-30)
            + 0.0);
        const amrex::Real troe_c = -0.4 - 0.67 * logFcent;
        const amrex::Real troe_n = 0.75 - 1.27 * logFcent;
        const amrex::Real troe = (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
        const amrex::Real F_troe = pow(10, logFcent / (1.0 + troe * troe));
        Corr = F * F_troe;
        const amrex::Real qf = Corr * k_f * (sc[1]*sc[3]);
        const amrex::Real qr = Corr * k_f * exp(-(g_RT[1] + g_RT[3] - g_RT[6])) * (refC) * (sc[6]);
        const amrex::Real qdot = qf - qr;
        wdot[1] -= qdot;
        wdot[3] -= qdot;
        wdot[6] += qdot;
    }

    {
        // (15):  H2O2 (+M) <=> OH + OH (+M)
        const amrex::Real k_f = 295100000000000
                   * exp(-(24370.7831249226) * invT);
        amrex::Real Corr = mixture + (1.5)*sc[0] + (11)*sc[2];
        const amrex::Real redP = Corr / k_f * 120200000000 
                   * exp(0 * tc[0] - 22896.3582941147 * invT);
        const amrex::Real F = redP / (1.0 + redP);
        const amrex::Real logPred = log10(redP);
        const amrex::Real logFcent = log10(
            0.5 * exp(-tc[1] * 1e+30)
            + 0.5 * exp(-tc[1] * 1e-30)
            + 0.0);
        const amrex::Real troe_c = -0.4 - 0.67 * logFcent;
        const amrex::Real troe_n = 0.75 - 1.27 * logFcent;
        const amrex::Real troe = (troe_c + logPred) / (troe_n - 0.14 * (troe_c + logPred));
        const amrex::Real F_troe = pow(10, logFcent / (1.0 + troe * troe));
        Corr = F * F_troe;
        const amrex::Real qf = Corr * k_f * (sc[7]);
        const amrex::Real qr = Corr * k_f * exp(-(-g_RT[5] - g_RT[5] + g_RT[7])) * (refCinv) * (sc[5]*sc[5]);
        const amrex::Real qdot = qf - qr;
        wdot[5] += qdot;
        wdot[5] += qdot;
        wdot[7] -= qdot;
    }

    {
        // (4):  H2 + M <=> H + H + M
        const amrex::Real k_f = 45770000000000
                   * exp((-1.4) * tc[0] - (52525.7555766966) * invT);
        const amrex::Real Corr = mixture + (1.5)*sc[0] + (11)*sc[2];
        const amrex::Real qf = Corr * k_f * (sc[0]);
        const amrex::Real qr = Corr * k_f * exp(-(g_RT[0] - g_RT[3] - g_RT[3])) * (refCinv) * (sc[3]*sc[3]);
        const amrex::Real qdot = qf - qr;
        wdot[0] -= qdot;
        wdot[3] += qdot;
        wdot[3] += qdot;
    }

    {
        // (5):  O + O + M <=> O2 + M
        const amrex::Real k_f = 6165
                   * exp((-0.5) * tc[0]);
        const amrex::Real Corr = mixture + (1.5)*sc[0] + (11)*sc[2];
        const amrex::Real qf = Corr * k_f * (sc[4]*sc[4]);
        const amrex::Real qr = Corr * k_f * exp(-(-g_RT[1] + g_RT[4] + g_RT[4])) * (refC) * (sc[1]);
        const amrex::Real qdot = qf - qr;
        wdot[1] += qdot;
        wdot[4] -= qdot;
        wdot[4] -= qdot;
    }

    {
        // (6):  O + H + M <=> OH + M
        const amrex::Real k_f = 4714000
                   * exp((-1) * tc[0]);
        const amrex::Real Corr = mixture + (1.5)*sc[0] + (11)*sc[2];
        const amrex::Real qf = Corr * k_f * (sc[3]*sc[4]);
        const amrex::Real qr = Corr * k_f * exp(-(g_RT[3] + g_RT[4] - g_RT[5])) * (refC) * (sc[5]);
        const amrex::Real qdot = qf - qr;
        wdot[3] -= qdot;
        wdot[4] -= qdot;
        wdot[5] += qdot;
    }

    {
        // (7):  H + OH + M <=> H2O + M
        const amrex::Real k_f = 38000000000
                   * exp((-2) * tc[0]);
        const amrex::Real Corr = mixture + (1.5)*sc[0] + (11)*sc[2];
        const amrex::Real qf = Corr * k_f * (sc[3]*sc[5]);
        const amrex::Real qr = Corr * k_f * exp(-(-g_RT[2] + g_RT[3] + g_RT[5])) * (refC) * (sc[2]);
        const amrex::Real qdot = qf - qr;
        wdot[2] += qdot;
        wdot[3] -= qdot;
        wdot[5] -= qdot;
    }

    {
        // (0):  H + O2 <=> O + OH
        const amrex::Real k_f = 3547000000
                   * exp((-0.406) * tc[0] - (8352.89343569254) * invT);
        const amrex::Real qf = k_f * (sc[1]*sc[3]);
        const amrex::Real qr = k_f * exp(-(g_RT[1] + g_RT[3] - g_RT[4] - g_RT[5])) * (sc[4]*sc[5]);
        const amrex::Real qdot = qf - qr;
        wdot[1] -= qdot;
        wdot[3] -= qdot;
        wdot[4] += qdot;
        wdot[5] += qdot;
    }

    {
        // (1):  O + H2 <=> H + OH
        const amrex::Real k_f = 0.0508
                   * exp((2.67) * tc[0] - (3165.23282791169) * invT);
        const amrex::Real qf = k_f * (sc[0]*sc[4]);
        const amrex::Real qr = k_f * exp(-(g_RT[0] - g_RT[3] + g_RT[4] - g_RT[5])) * (sc[3]*sc[5]);
        const amrex::Real qdot = qf - qr;
        wdot[0] -= qdot;
        wdot[3] += qdot;
        wdot[4] -= qdot;
        wdot[5] += qdot;
    }

    {
        // (2):  H2 + OH <=> H2O + H
        const amrex::Real k_f = 216
                   * exp((1.51) * tc[0] - (1726.03316371019) * invT);
        const amrex::Real qf = k_f * (sc[0]*sc[5]);
        const amrex::Real qr = k_f * exp(-(g_RT[0] - g_RT[2] - g_RT[3] + g_RT[5])) * (sc[2]*sc[3]);
        const amrex::Real qdot = qf - qr;
        wdot[0] -= qdot;
        wdot[2] += qdot;
        wdot[3] += qdot;
        wdot[5] -= qdot;
    }

    {
        // (3):  O + H2O <=> OH + OH
        const amrex::Real k_f = 2.97
                   * exp((2.02) * tc[0] - (6743.10332178324) * invT);
        const amrex::Real qf = k_f * (sc[2]*sc[4]);
        const amrex::Real qr = k_f * exp(-(g_RT[2] + g_RT[4] - g_RT[5] - g_RT[5])) * (sc[5]*sc[5]);
        const amrex::Real qdot = qf - qr;
        wdot[2] -= qdot;
        wdot[4] -= qdot;
        wdot[5] += qdot;
        wdot[5] += qdot;
    }

    {
        // (9):  HO2 + H <=> H2 + O2
        const amrex::Real k_f = 16600000
                   * exp(-(414.147315957284) * invT);
        const amrex::Real qf = k_f * (sc[3]*sc[6]);
        const amrex::Real qr = k_f * exp(-(-g_RT[0] - g_RT[1] + g_RT[3] + g_RT[6])) * (sc[0]*sc[1]);
        const amrex::Real qdot = qf - qr;
        wdot[0] += qdot;
        wdot[1] += qdot;
        wdot[3] -= qdot;
        wdot[6] -= qdot;
    }

    {
        // (10):  HO2 + H <=> OH + OH
        const amrex::Real k_f = 70790000
                   * exp(-(148.448916412392) * invT);
        const amrex::Real qf = k_f * (sc[3]*sc[6]);
        const amrex::Real qr = k_f * exp(-(g_RT[3] - g_RT[5] - g_RT[5] + g_RT[6])) * (sc[5]*sc[5]);
        const amrex::Real qdot = qf - qr;
        wdot[3] -= qdot;
        wdot[5] += qdot;
        wdot[5] += qdot;
        wdot[6] -= qdot;
    }

    {
        // (11):  HO2 + O <=> O2 + OH
        const amrex::Real k_f = 32500000
                   ;
        const amrex::Real qf = k_f * (sc[4]*sc[6]);
        const amrex::Real qr = k_f * exp(-(-g_RT[1] + g_RT[4] - g_RT[5] + g_RT[6])) * (sc[1]*sc[5]);
        const amrex::Real qdot = qf - qr;
        wdot[1] += qdot;
        wdot[4] -= qdot;
        wdot[5] += qdot;
        wdot[6] -= qdot;
    }

    {
        // (12):  HO2 + OH <=> H2O + O2
        const amrex::Real k_f = 28900000
                   * exp(-(-250.098682904946) * invT);
        const amrex::Real qf = k_f * (sc[5]*sc[6]);
        const amrex::Real qr = k_f * exp(-(-g_RT[1] - g_RT[2] + g_RT[5] + g_RT[6])) * (sc[1]*sc[2]);
        const amrex::Real qdot = qf - qr;
        wdot[1] += qdot;
        wdot[2] += qdot;
        wdot[5] -= qdot;
        wdot[6] -= qdot;
    }

    {
        // (13):  HO2 + HO2 <=> H2O2 + O2
        const amrex::Real k_f = 420000000
                   * exp(-(6029.54208967215) * invT);
        const amrex::Real qf = k_f * (sc[6]*sc[6]);
        const amrex::Real qr = k_f * exp(-(-g_RT[1] + g_RT[6] + g_RT[6] - g_RT[7])) * (sc[1]*sc[7]);
        const amrex::Real qdot = qf - qr;
        wdot[1] += qdot;
        wdot[6] -= qdot;
        wdot[6] -= qdot;
        wdot[7] += qdot;
    }

    {
        // (14):  HO2 + HO2 <=> H2O2 + O2
        const amrex::Real k_f = 130000
                   * exp(-(-819.89091359563) * invT);
        const amrex::Real qf = k_f * (sc[6]*sc[6]);
        const amrex::Real qr = k_f * exp(-(-g_RT[1] + g_RT[6] + g_RT[6] - g_RT[7])) * (sc[1]*sc[7]);
        const amrex::Real qdot = qf - qr;
        wdot[1] += qdot;
        wdot[6] -= qdot;
        wdot[6] -= qdot;
        wdot[7] += qdot;
    }

    {
        // (16):  H2O2 + H <=> H2O + OH
        const amrex::Real k_f = 24100000
                   * exp(-(1997.77016324474) * invT);
        const amrex::Real qf = k_f * (sc[3]*sc[7]);
        const amrex::Real qr = k_f * exp(-(-g_RT[2] + g_RT[3] - g_RT[5] + g_RT[7])) * (sc[2]*sc[5]);
        const amrex::Real qdot = qf - qr;
        wdot[2] += qdot;
        wdot[3] -= qdot;
        wdot[5] += qdot;
        wdot[7] -= qdot;
    }

    {
        // (17):  H2O2 + H <=> HO2 + H2
        const amrex::Real k_f = 48200000
                   * exp(-(4000.57249314752) * invT);
        const amrex::Real qf = k_f * (sc[3]*sc[7]);
        const amrex::Real qr = k_f * exp(-(-g_RT[0] + g_RT[3] - g_RT[6] + g_RT[7])) * (sc[0]*sc[6]);
        const amrex::Real qdot = qf - qr;
        wdot[0] += qdot;
        wdot[3] -= qdot;
        wdot[6] += qdot;
        wdot[7] -= qdot;
    }

    {
        // (18):  H2O2 + O <=> OH + HO2
        const amrex::Real k_f = 9.55
                   * exp((2) * tc[0] - (1997.77016324474) * invT);
        const amrex::Real qf = k_f * (sc[4]*sc[7]);
        const amrex::Real qr = k_f * exp(-(g_RT[4] - g_RT[5] - g_RT[6] + g_RT[7])) * (sc[5]*sc[6]);
        const amrex::Real qdot = qf - qr;
        wdot[4] -= qdot;
        wdot[5] += qdot;
        wdot[6] += qdot;
        wdot[7] -= qdot;
    }

    {
        // (19):  H2O2 + OH <=> HO2 + H2O
        const amrex::Real k_f = 1000000
                   ;
        const amrex::Real qf = k_f * (sc[5]*sc[7]);
        const amrex::Real qr = k_f * exp(-(-g_RT[2] + g_RT[5] - g_RT[6] + g_RT[7])) * (sc[2]*sc[6]);
        const amrex::Real qdot = qf - qr;
        wdot[2] += qdot;
        wdot[5] -= qdot;
        wdot[6] += qdot;
        wdot[7] -= qdot;
    }

    {
        // (20):  H2O2 + OH <=> HO2 + H2O
        const amrex::Real k_f = 580000000
                   * exp(-(4809.24167509571) * invT);
        const amrex::Real qf = k_f * (sc[5]*sc[7]);
        const amrex::Real qr = k_f * exp(-(-g_RT[2] + g_RT[5] - g_RT[6] + g_RT[7])) * (sc[2]*sc[6]);
        const amrex::Real qdot = qf - qr;
        wdot[2] += qdot;
        wdot[5] -= qdot;
        wdot[6] += qdot;
        wdot[7] -= qdot;
    }


    return;
}



/*compute the production rate for each species */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWC(amrex::Real *  T, amrex::Real *  C,  amrex::Real *  wdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 9; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    productionRate(wdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        C[id] *= 1.0e-6;
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mass fractions */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWYP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  wdot)
{
    int id; /*loop counter */
    amrex::Real c[9]; /*temporary storage */
    amrex::Real YOW = 0; 
    amrex::Real PWORT; 
    amrex::Real imw[9];

    get_imw(imw);

    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*O2 */
    YOW += y[2]*imw[2]; /*H2O */
    YOW += y[3]*imw[3]; /*H */
    YOW += y[4]*imw[4]; /*O */
    YOW += y[5]*imw[5]; /*OH */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*N2 */
    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31446261815324e+07 * (*T)); 
    /*multiply by 1e6 so c goes to SI */
    PWORT *= 1e6; 
    /*Now compute conversion (and go to SI) */
    c[0] = PWORT * y[0]*imw[0]; 
    c[1] = PWORT * y[1]*imw[1]; 
    c[2] = PWORT * y[2]*imw[2]; 
    c[3] = PWORT * y[3]*imw[3]; 
    c[4] = PWORT * y[4]*imw[4]; 
    c[5] = PWORT * y[5]*imw[5]; 
    c[6] = PWORT * y[6]*imw[6]; 
    c[7] = PWORT * y[7]*imw[7]; 
    c[8] = PWORT * y[8]*imw[8]; 

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mole fractions */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWXP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x,  amrex::Real *  wdot)
{
    int id; /*loop counter */
    amrex::Real c[9]; /*temporary storage */
    amrex::Real PORT = 1e6 * (*P)/(8.31446261815324e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 9; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWYR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  wdot)
{
    int id; /*loop counter */
    amrex::Real c[9]; /*temporary storage */
    amrex::Real imw[9];

    get_imw(imw);

    /*See Eq 8 with an extra 1e6 so c goes to SI */
    c[0] = 1e6 * (*rho) * y[0]*imw[0]; 
    c[1] = 1e6 * (*rho) * y[1]*imw[1]; 
    c[2] = 1e6 * (*rho) * y[2]*imw[2]; 
    c[3] = 1e6 * (*rho) * y[3]*imw[3]; 
    c[4] = 1e6 * (*rho) * y[4]*imw[4]; 
    c[5] = 1e6 * (*rho) * y[5]*imw[5]; 
    c[6] = 1e6 * (*rho) * y[6]*imw[6]; 
    c[7] = 1e6 * (*rho) * y[7]*imw[7]; 
    c[8] = 1e6 * (*rho) * y[8]*imw[8]; 

    /*call productionRate */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mole fractions */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWXR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  x,  amrex::Real *  wdot)
{
    int id; /*loop counter */
    amrex::Real c[9]; /*temporary storage */
    amrex::Real XW = 0; /*See Eq 4, 11 in CK Manual */
    amrex::Real ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*31.998800; /*O2 */
    XW += x[2]*18.015340; /*H2O */
    XW += x[3]*1.007970; /*H */
    XW += x[4]*15.999400; /*O */
    XW += x[5]*17.007370; /*OH */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 9; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 9; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*compute d(Cp/R)/dT and d(Cv/R)/dT at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void dcvpRdT(amrex::Real * species, const amrex::Real *  tc)
{

    /*temperature */
    const amrex::Real T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            +8.24944174e-04
            -1.62860306e-06 * tc[1]
            -2.84263030e-10 * tc[2]
            +1.65394890e-12 * tc[3];
        /*species 1: O2 */
        species[1] =
            +1.12748635e-03
            -1.15123009e-06 * tc[1]
            +3.94163169e-09 * tc[2]
            -3.50742157e-12 * tc[3];
        /*species 2: H2O */
        species[2] =
            +3.47498246e-03
            -1.27093927e-05 * tc[1]
            +2.09057438e-08 * tc[2]
            -1.00263539e-11 * tc[3];
        /*species 3: H */
        species[3] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
        /*species 4: O */
        species[4] =
            -1.63816649e-03
            +4.84206340e-06 * tc[1]
            -4.80852957e-09 * tc[2]
            +1.55627854e-12 * tc[3];
        /*species 5: OH */
        species[5] =
            -3.22544939e-03
            +1.30552938e-05 * tc[1]
            -1.73956093e-08 * tc[2]
            +8.24949516e-12 * tc[3];
        /*species 6: HO2 */
        species[6] =
            -4.74912051e-03
            +4.23165782e-05 * tc[1]
            -7.28291682e-08 * tc[2]
            +3.71690050e-11 * tc[3];
        /*species 7: H2O2 */
        species[7] =
            +6.56922581e-03
            -2.97002516e-07 * tc[1]
            -1.38774166e-08 * tc[2]
            +9.88605900e-12 * tc[3];
        /*species 8: N2 */
        species[8] =
            +1.40824000e-03
            -7.92644400e-06 * tc[1]
            +1.69245450e-08 * tc[2]
            -9.77942000e-12 * tc[3];
    } else {
        /*species 0: H2 */
        species[0] =
            +7.00064411e-04
            -1.12676574e-07 * tc[1]
            -2.76947345e-11 * tc[2]
            +6.33100716e-15 * tc[3];
        /*species 1: O2 */
        species[1] =
            +6.13519689e-04
            -2.51768398e-07 * tc[1]
            +5.32584444e-11 * tc[2]
            -4.54574124e-15 * tc[3];
        /*species 2: H2O */
        species[2] =
            +3.05629289e-03
            -1.74605202e-06 * tc[1]
            +3.60298917e-10 * tc[2]
            -2.55664715e-14 * tc[3];
        /*species 3: H */
        species[3] =
            +0.00000000e+00
            +0.00000000e+00 * tc[1]
            +0.00000000e+00 * tc[2]
            +0.00000000e+00 * tc[3];
        /*species 4: O */
        species[4] =
            -2.75506191e-05
            -6.20560670e-09 * tc[1]
            +1.36532023e-11 * tc[2]
            -1.74722060e-15 * tc[3];
        /*species 5: OH */
        species[5] =
            +1.05650448e-03
            -5.18165516e-07 * tc[1]
            +9.15656022e-11 * tc[2]
            -5.32783504e-15 * tc[3];
        /*species 6: HO2 */
        species[6] =
            +2.23982013e-03
            -1.26731630e-06 * tc[1]
            +3.42739110e-10 * tc[2]
            -4.31634140e-14 * tc[3];
        /*species 7: H2O2 */
        species[7] =
            +4.33613639e-03
            -2.94937764e-06 * tc[1]
            +7.04671071e-10 * tc[2]
            -5.72661424e-14 * tc[3];
        /*species 8: N2 */
        species[8] =
            +1.48797700e-03
            -1.13695220e-06 * tc[1]
            +3.02911200e-10 * tc[2]
            -2.70134040e-14 * tc[3];
    }
    return;
}

/*compute an approx to the reaction Jacobian */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void aJacobian_precond(amrex::Real *  J, amrex::Real *  sc, amrex::Real T, int HP)
{
    for (int i=0; i<100; i++) {
        J[i] = 0.0;
    }

    amrex::Real wdot[9];
    for (int k=0; k<9; k++) {
        wdot[k] = 0.0;
    }

    const amrex::Real tc[5] = { log(T), T, T*T, T*T*T, T*T*T*T }; // temperature cache
    amrex::Real invT = 1.0 / tc[1];
    amrex::Real invT2 = invT * invT;

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    amrex::Real refC = 101325 / 8.31446 / T;
    amrex::Real refCinv = 1.0 / refC;

    /*compute the mixture concentration */
    amrex::Real mixture = 0.0;
    for (int k = 0; k < 9; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    amrex::Real g_RT[9];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    amrex::Real h_RT[9];
    speciesEnthalpy(h_RT, tc);

    amrex::Real phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    amrex::Real dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    amrex::Real dqdci, dcdc_fac, dqdc[9];
    amrex::Real Pr, fPr, F, k_0, logPr;
    amrex::Real logFcent, troe_c, troe_n, troePr_den, troePr, troe;
    amrex::Real Fcent1, Fcent2, Fcent3, Fcent;
    amrex::Real dlogFdc, dlogFdn, dlogFdcn_fac;
    amrex::Real dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
    const amrex::Real ln10 = log(10.0);
    const amrex::Real log10e = 1.0/log(10.0);
    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    /*a pressure-fall-off reaction */
    // also 3-body
    // 3-body correction factor
    alpha = mixture + ( 2 - 1)*sc[0] + ( 11 - 1)*sc[2] + ( 0.78 - 1)*sc[1];
    // forward
    phi_f = sc[1]*sc[3];
    k_f = 1e-06 * 1475000000000
                * exp(0.6 * tc[0] - 0.50321666580472 * (0) * invT);
    dlnkfdT = 0.6 * invT + 0.50321666580472 *  (0)  * invT2;
    // pressure-fall-off
    k_0 = 6.366e+20 * exp(-1.72 * tc[0] - 0.50321666580472 * (524.8) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -1.72 * invT + 0.50321666580472 * (524.8) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    // Troe form
    logPr = log10(Pr);
    Fcent1 = (1.-0.8)*exp(-T/1e-30);
    Fcent2 = 0.8 * exp(-T/1e+30);
    Fcent3 = 0.;
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/1e-30
        -Fcent2/1e+30
    );
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    // reverse
    phi_r = sc[6];
    Kc = refCinv * exp(g_RT[1] + g_RT[3] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3]) + (h_RT[6]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    // rate of progress
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    // update wdot
    wdot[1] -= q; // O2
    wdot[3] -= q; // H
    wdot[6] += q; // HO2
    // for convenience
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = 0.0;
    dqdc[0] = 2*dcdc_fac;
    dqdc[1] = 0.78*dcdc_fac + k_f*sc[3];
    dqdc[2] = 11*dcdc_fac;
    dqdc[3] = dcdc_fac + k_f*sc[1];
    dqdc[4] = dcdc_fac;
    dqdc[5] = dcdc_fac;
    dqdc[6] = dcdc_fac - k_r;
    dqdc[7] = dcdc_fac;
    dqdc[8] = dcdc_fac;
    for (int k=0; k<9; k++) {
        J[10*k+1] -= dqdc[k];
        J[10*k+3] -= dqdc[k];
        J[10*k+6] += dqdc[k];
    }
    J[91] -= dqdT; // dwdot[O2]/dT
    J[93] -= dqdT; // dwdot[H]/dT
    J[96] += dqdT; // dwdot[HO2]/dT

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    /*a pressure-fall-off reaction */
    // also 3-body
    // 3-body correction factor
    alpha = mixture + ( 2.5 - 1)*sc[0] + ( 12 - 1)*sc[2];
    // forward
    phi_f = sc[7];
    k_f = 1 * 295100000000000
                * exp(0 * tc[0] - 0.50321666580472 * (48430) * invT);
    dlnkfdT = 0 * invT + 0.50321666580472 *  (48430)  * invT2;
    // pressure-fall-off
    k_0 = 1.202e+17 * exp(0 * tc[0] - 0.50321666580472 * (45500) * invT);
    Pr = 1e-6 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = 0 * invT + 0.50321666580472 * (45500) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    // Troe form
    logPr = log10(Pr);
    Fcent1 = (1.-0.5)*exp(-T/1e-30);
    Fcent2 = 0.5 * exp(-T/1e+30);
    Fcent3 = 0.;
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/1e-30
        -Fcent2/1e+30
    );
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    // reverse
    phi_r = (sc[5] * sc[5]);
    Kc = refC * exp(-g_RT[5] - g_RT[5] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7]) + (2.000000*h_RT[5]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    // rate of progress
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    // update wdot
    wdot[5] += 2 * q; // OH
    wdot[7] -= q; // H2O2
    // for convenience
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = 0.0;
    dqdc[0] = 2.5*dcdc_fac;
    dqdc[1] = dcdc_fac;
    dqdc[2] = 12*dcdc_fac;
    dqdc[3] = dcdc_fac;
    dqdc[4] = dcdc_fac;
    dqdc[5] = dcdc_fac - k_r*2.000000*sc[5];
    dqdc[6] = dcdc_fac;
    dqdc[7] = dcdc_fac + k_f;
    dqdc[8] = dcdc_fac;
    for (int k=0; k<9; k++) {
        J[10*k+5] += 2 * dqdc[k];
        J[10*k+7] -= dqdc[k];
    }
    J[95] += 2 * dqdT; // dwdot[OH]/dT
    J[97] -= dqdT; // dwdot[H2O2]/dT

    /*reaction 3: H2 + M <=> H + H + M */
    /*a third-body and non-pressure-fall-off reaction */
    // 3-body correction factor
    alpha = mixture + ( 2.5 - 1)*sc[0] + ( 12 - 1)*sc[2];
    // forward
    phi_f = sc[0];
    k_f = 1e-06 * 4.577e+19
                * exp(-1.4 * tc[0] - 0.50321666580472 * (104380) * invT);
    dlnkfdT = -1.4 * invT + 0.50321666580472 *  (104380)  * invT2;
    // reverse
    phi_r = (sc[3] * sc[3]);
    Kc = refC * exp(g_RT[0] - g_RT[3] - g_RT[3]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0]) + (2.000000*h_RT[3]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    // rate of progress
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    // update wdot
    wdot[0] -= q; // H2
    wdot[3] += 2 * q; // H
    // for convenience
    k_f *= alpha;
    k_r *= alpha;
    dqdc[0] = 2.5*q_nocor + k_f;
    dqdc[1] = q_nocor;
    dqdc[2] = 12*q_nocor;
    dqdc[3] = q_nocor - k_r*2.000000*sc[3];
    dqdc[4] = q_nocor;
    dqdc[5] = q_nocor;
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    for (int k=0; k<9; k++) {
        J[10*k+0] -= dqdc[k];
        J[10*k+3] += 2 * dqdc[k];
    }
    J[90] -= dqdT; // dwdot[H2]/dT
    J[93] += 2 * dqdT; // dwdot[H]/dT

    /*reaction 4: O + O + M <=> O2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    // 3-body correction factor
    alpha = mixture + ( 2.5 - 1)*sc[0] + ( 12 - 1)*sc[2];
    // forward
    phi_f = (sc[4] * sc[4]);
    k_f = 1e-12 * 6.165e+15
                * exp(-0.5 * tc[0] - 0.50321666580472 * (0) * invT);
    dlnkfdT = -0.5 * invT + 0.50321666580472 *  (0)  * invT2;
    // reverse
    phi_r = sc[1];
    Kc = refCinv * exp(-g_RT[1] + g_RT[4] + g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[4]) + (h_RT[1]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    // rate of progress
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    // update wdot
    wdot[1] += q; // O2
    wdot[4] -= 2 * q; // O
    // for convenience
    k_f *= alpha;
    k_r *= alpha;
    dqdc[0] = 2.5*q_nocor;
    dqdc[1] = q_nocor - k_r;
    dqdc[2] = 12*q_nocor;
    dqdc[3] = q_nocor;
    dqdc[4] = q_nocor + k_f*2.000000*sc[4];
    dqdc[5] = q_nocor;
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    for (int k=0; k<9; k++) {
        J[10*k+1] += dqdc[k];
        J[10*k+4] += -2 * dqdc[k];
    }
    J[91] += dqdT; // dwdot[O2]/dT
    J[94] += -2 * dqdT; // dwdot[O]/dT

    /*reaction 5: O + H + M <=> OH + M */
    /*a third-body and non-pressure-fall-off reaction */
    // 3-body correction factor
    alpha = mixture + ( 2.5 - 1)*sc[0] + ( 12 - 1)*sc[2];
    // forward
    phi_f = sc[3]*sc[4];
    k_f = 1e-12 * 4.714e+18
                * exp(-1 * tc[0] - 0.50321666580472 * (0) * invT);
    dlnkfdT = -1 * invT + 0.50321666580472 *  (0)  * invT2;
    // reverse
    phi_r = sc[5];
    Kc = refCinv * exp(g_RT[3] + g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[4]) + (h_RT[5]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    // rate of progress
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    // update wdot
    wdot[3] -= q; // H
    wdot[4] -= q; // O
    wdot[5] += q; // OH
    // for convenience
    k_f *= alpha;
    k_r *= alpha;
    dqdc[0] = 2.5*q_nocor;
    dqdc[1] = q_nocor;
    dqdc[2] = 12*q_nocor;
    dqdc[3] = q_nocor + k_f*sc[4];
    dqdc[4] = q_nocor + k_f*sc[3];
    dqdc[5] = q_nocor - k_r;
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    for (int k=0; k<9; k++) {
        J[10*k+3] -= dqdc[k];
        J[10*k+4] -= dqdc[k];
        J[10*k+5] += dqdc[k];
    }
    J[93] -= dqdT; // dwdot[H]/dT
    J[94] -= dqdT; // dwdot[O]/dT
    J[95] += dqdT; // dwdot[OH]/dT

    /*reaction 6: H + OH + M <=> H2O + M */
    /*a third-body and non-pressure-fall-off reaction */
    // 3-body correction factor
    alpha = mixture + ( 2.5 - 1)*sc[0] + ( 12 - 1)*sc[2];
    // forward
    phi_f = sc[3]*sc[5];
    k_f = 1e-12 * 3.8e+22
                * exp(-2 * tc[0] - 0.50321666580472 * (0) * invT);
    dlnkfdT = -2 * invT + 0.50321666580472 *  (0)  * invT2;
    // reverse
    phi_r = sc[2];
    Kc = refCinv * exp(-g_RT[2] + g_RT[3] + g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[5]) + (h_RT[2]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    // rate of progress
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    // update wdot
    wdot[2] += q; // H2O
    wdot[3] -= q; // H
    wdot[5] -= q; // OH
    // for convenience
    k_f *= alpha;
    k_r *= alpha;
    dqdc[0] = 2.5*q_nocor;
    dqdc[1] = q_nocor;
    dqdc[2] = 12*q_nocor - k_r;
    dqdc[3] = q_nocor + k_f*sc[5];
    dqdc[4] = q_nocor;
    dqdc[5] = q_nocor + k_f*sc[3];
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    for (int k=0; k<9; k++) {
        J[10*k+2] += dqdc[k];
        J[10*k+3] -= dqdc[k];
        J[10*k+5] -= dqdc[k];
    }
    J[92] += dqdT; // dwdot[H2O]/dT
    J[93] -= dqdT; // dwdot[H]/dT
    J[95] -= dqdT; // dwdot[OH]/dT

    /*reaction 7: H + O2 <=> O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    // forward
    phi_f = sc[1]*sc[3];
    k_f = 1e-06 * 3.547e+15
                * exp(-0.406 * tc[0] - 0.50321666580472 * (16599) * invT);
    dlnkfdT = -0.406 * invT + 0.50321666580472 *  (16599)  * invT2;
    // reverse
    phi_r = sc[4]*sc[5];
    Kc = exp(g_RT[1] + g_RT[3] - g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3]) + (h_RT[4] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    // rate of progress
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    // update wdot
    wdot[1] -= q; // O2
    wdot[3] -= q; // H
    wdot[4] += q; // O
    wdot[5] += q; // OH
    // d()/d[O2]
    dqdci =  + k_f*sc[3];
    J[11] -= dqdci;               // dwdot[O2]/d[O2]
    J[13] -= dqdci;               // dwdot[H]/d[O2]
    J[14] += dqdci;               // dwdot[O]/d[O2]
    J[15] += dqdci;               // dwdot[OH]/d[O2]
    // d()/d[H]
    dqdci =  + k_f*sc[1];
    J[31] -= dqdci;               // dwdot[O2]/d[H]
    J[33] -= dqdci;               // dwdot[H]/d[H]
    J[34] += dqdci;               // dwdot[O]/d[H]
    J[35] += dqdci;               // dwdot[OH]/d[H]
    // d()/d[O]
    dqdci =  - k_r*sc[5];
    J[41] -= dqdci;               // dwdot[O2]/d[O]
    J[43] -= dqdci;               // dwdot[H]/d[O]
    J[44] += dqdci;               // dwdot[O]/d[O]
    J[45] += dqdci;               // dwdot[OH]/d[O]
    // d()/d[OH]
    dqdci =  - k_r*sc[4];
    J[51] -= dqdci;               // dwdot[O2]/d[OH]
    J[53] -= dqdci;               // dwdot[H]/d[OH]
    J[54] += dqdci;               // dwdot[O]/d[OH]
    J[55] += dqdci;               // dwdot[OH]/d[OH]
    // d()/dT
    J[91] -= dqdT;                // dwdot[O2]/dT
    J[93] -= dqdT;                // dwdot[H]/dT
    J[94] += dqdT;                // dwdot[O]/dT
    J[95] += dqdT;                // dwdot[OH]/dT

    /*reaction 8: O + H2 <=> H + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    // forward
    phi_f = sc[0]*sc[4];
    k_f = 1e-06 * 50800
                * exp(2.67 * tc[0] - 0.50321666580472 * (6290) * invT);
    dlnkfdT = 2.67 * invT + 0.50321666580472 *  (6290)  * invT2;
    // reverse
    phi_r = sc[3]*sc[5];
    Kc = exp(g_RT[0] - g_RT[3] + g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[4]) + (h_RT[3] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    // rate of progress
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    // update wdot
    wdot[0] -= q; // H2
    wdot[3] += q; // H
    wdot[4] -= q; // O
    wdot[5] += q; // OH
    // d()/d[H2]
    dqdci =  + k_f*sc[4];
    J[0] -= dqdci;                // dwdot[H2]/d[H2]
    J[3] += dqdci;                // dwdot[H]/d[H2]
    J[4] -= dqdci;                // dwdot[O]/d[H2]
    J[5] += dqdci;                // dwdot[OH]/d[H2]
    // d()/d[H]
    dqdci =  - k_r*sc[5];
    J[30] -= dqdci;               // dwdot[H2]/d[H]
    J[33] += dqdci;               // dwdot[H]/d[H]
    J[34] -= dqdci;               // dwdot[O]/d[H]
    J[35] += dqdci;               // dwdot[OH]/d[H]
    // d()/d[O]
    dqdci =  + k_f*sc[0];
    J[40] -= dqdci;               // dwdot[H2]/d[O]
    J[43] += dqdci;               // dwdot[H]/d[O]
    J[44] -= dqdci;               // dwdot[O]/d[O]
    J[45] += dqdci;               // dwdot[OH]/d[O]
    // d()/d[OH]
    dqdci =  - k_r*sc[3];
    J[50] -= dqdci;               // dwdot[H2]/d[OH]
    J[53] += dqdci;               // dwdot[H]/d[OH]
    J[54] -= dqdci;               // dwdot[O]/d[OH]
    J[55] += dqdci;               // dwdot[OH]/d[OH]
    // d()/dT
    J[90] -= dqdT;                // dwdot[H2]/dT
    J[93] += dqdT;                // dwdot[H]/dT
    J[94] -= dqdT;                // dwdot[O]/dT
    J[95] += dqdT;                // dwdot[OH]/dT

    /*reaction 9: H2 + OH <=> H2O + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    // forward
    phi_f = sc[0]*sc[5];
    k_f = 1e-06 * 216000000
                * exp(1.51 * tc[0] - 0.50321666580472 * (3430) * invT);
    dlnkfdT = 1.51 * invT + 0.50321666580472 *  (3430)  * invT2;
    // reverse
    phi_r = sc[2]*sc[3];
    Kc = exp(g_RT[0] - g_RT[2] - g_RT[3] + g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[5]) + (h_RT[2] + h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    // rate of progress
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    // update wdot
    wdot[0] -= q; // H2
    wdot[2] += q; // H2O
    wdot[3] += q; // H
    wdot[5] -= q; // OH
    // d()/d[H2]
    dqdci =  + k_f*sc[5];
    J[0] -= dqdci;                // dwdot[H2]/d[H2]
    J[2] += dqdci;                // dwdot[H2O]/d[H2]
    J[3] += dqdci;                // dwdot[H]/d[H2]
    J[5] -= dqdci;                // dwdot[OH]/d[H2]
    // d()/d[H2O]
    dqdci =  - k_r*sc[3];
    J[20] -= dqdci;               // dwdot[H2]/d[H2O]
    J[22] += dqdci;               // dwdot[H2O]/d[H2O]
    J[23] += dqdci;               // dwdot[H]/d[H2O]
    J[25] -= dqdci;               // dwdot[OH]/d[H2O]
    // d()/d[H]
    dqdci =  - k_r*sc[2];
    J[30] -= dqdci;               // dwdot[H2]/d[H]
    J[32] += dqdci;               // dwdot[H2O]/d[H]
    J[33] += dqdci;               // dwdot[H]/d[H]
    J[35] -= dqdci;               // dwdot[OH]/d[H]
    // d()/d[OH]
    dqdci =  + k_f*sc[0];
    J[50] -= dqdci;               // dwdot[H2]/d[OH]
    J[52] += dqdci;               // dwdot[H2O]/d[OH]
    J[53] += dqdci;               // dwdot[H]/d[OH]
    J[55] -= dqdci;               // dwdot[OH]/d[OH]
    // d()/dT
    J[90] -= dqdT;                // dwdot[H2]/dT
    J[92] += dqdT;                // dwdot[H2O]/dT
    J[93] += dqdT;                // dwdot[H]/dT
    J[95] -= dqdT;                // dwdot[OH]/dT

    /*reaction 10: O + H2O <=> OH + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    // forward
    phi_f = sc[2]*sc[4];
    k_f = 1e-06 * 2970000
                * exp(2.02 * tc[0] - 0.50321666580472 * (13400) * invT);
    dlnkfdT = 2.02 * invT + 0.50321666580472 *  (13400)  * invT2;
    // reverse
    phi_r = (sc[5] * sc[5]);
    Kc = exp(g_RT[2] + g_RT[4] - g_RT[5] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[4]) + (2.000000*h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    // rate of progress
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    // update wdot
    wdot[2] -= q; // H2O
    wdot[4] -= q; // O
    wdot[5] += 2 * q; // OH
    // d()/d[H2O]
    dqdci =  + k_f*sc[4];
    J[22] -= dqdci;               // dwdot[H2O]/d[H2O]
    J[24] -= dqdci;               // dwdot[O]/d[H2O]
    J[25] += 2 * dqdci;           // dwdot[OH]/d[H2O]
    // d()/d[O]
    dqdci =  + k_f*sc[2];
    J[42] -= dqdci;               // dwdot[H2O]/d[O]
    J[44] -= dqdci;               // dwdot[O]/d[O]
    J[45] += 2 * dqdci;           // dwdot[OH]/d[O]
    // d()/d[OH]
    dqdci =  - k_r*2.000000*sc[5];
    J[52] -= dqdci;               // dwdot[H2O]/d[OH]
    J[54] -= dqdci;               // dwdot[O]/d[OH]
    J[55] += 2 * dqdci;           // dwdot[OH]/d[OH]
    // d()/dT
    J[92] -= dqdT;                // dwdot[H2O]/dT
    J[94] -= dqdT;                // dwdot[O]/dT
    J[95] += 2 * dqdT;            // dwdot[OH]/dT

    /*reaction 11: HO2 + H <=> H2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    // forward
    phi_f = sc[3]*sc[6];
    k_f = 1e-06 * 16600000000000
                * exp(0 * tc[0] - 0.50321666580472 * (823) * invT);
    dlnkfdT = 0 * invT + 0.50321666580472 *  (823)  * invT2;
    // reverse
    phi_r = sc[0]*sc[1];
    Kc = exp(-g_RT[0] - g_RT[1] + g_RT[3] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[6]) + (h_RT[0] + h_RT[1]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    // rate of progress
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    // update wdot
    wdot[0] += q; // H2
    wdot[1] += q; // O2
    wdot[3] -= q; // H
    wdot[6] -= q; // HO2
    // d()/d[H2]
    dqdci =  - k_r*sc[1];
    J[0] += dqdci;                // dwdot[H2]/d[H2]
    J[1] += dqdci;                // dwdot[O2]/d[H2]
    J[3] -= dqdci;                // dwdot[H]/d[H2]
    J[6] -= dqdci;                // dwdot[HO2]/d[H2]
    // d()/d[O2]
    dqdci =  - k_r*sc[0];
    J[10] += dqdci;               // dwdot[H2]/d[O2]
    J[11] += dqdci;               // dwdot[O2]/d[O2]
    J[13] -= dqdci;               // dwdot[H]/d[O2]
    J[16] -= dqdci;               // dwdot[HO2]/d[O2]
    // d()/d[H]
    dqdci =  + k_f*sc[6];
    J[30] += dqdci;               // dwdot[H2]/d[H]
    J[31] += dqdci;               // dwdot[O2]/d[H]
    J[33] -= dqdci;               // dwdot[H]/d[H]
    J[36] -= dqdci;               // dwdot[HO2]/d[H]
    // d()/d[HO2]
    dqdci =  + k_f*sc[3];
    J[60] += dqdci;               // dwdot[H2]/d[HO2]
    J[61] += dqdci;               // dwdot[O2]/d[HO2]
    J[63] -= dqdci;               // dwdot[H]/d[HO2]
    J[66] -= dqdci;               // dwdot[HO2]/d[HO2]
    // d()/dT
    J[90] += dqdT;                // dwdot[H2]/dT
    J[91] += dqdT;                // dwdot[O2]/dT
    J[93] -= dqdT;                // dwdot[H]/dT
    J[96] -= dqdT;                // dwdot[HO2]/dT

    /*reaction 12: HO2 + H <=> OH + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    // forward
    phi_f = sc[3]*sc[6];
    k_f = 1e-06 * 70790000000000
                * exp(0 * tc[0] - 0.50321666580472 * (295) * invT);
    dlnkfdT = 0 * invT + 0.50321666580472 *  (295)  * invT2;
    // reverse
    phi_r = (sc[5] * sc[5]);
    Kc = exp(g_RT[3] - g_RT[5] - g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[6]) + (2.000000*h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    // rate of progress
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    // update wdot
    wdot[3] -= q; // H
    wdot[5] += 2 * q; // OH
    wdot[6] -= q; // HO2
    // d()/d[H]
    dqdci =  + k_f*sc[6];
    J[33] -= dqdci;               // dwdot[H]/d[H]
    J[35] += 2 * dqdci;           // dwdot[OH]/d[H]
    J[36] -= dqdci;               // dwdot[HO2]/d[H]
    // d()/d[OH]
    dqdci =  - k_r*2.000000*sc[5];
    J[53] -= dqdci;               // dwdot[H]/d[OH]
    J[55] += 2 * dqdci;           // dwdot[OH]/d[OH]
    J[56] -= dqdci;               // dwdot[HO2]/d[OH]
    // d()/d[HO2]
    dqdci =  + k_f*sc[3];
    J[63] -= dqdci;               // dwdot[H]/d[HO2]
    J[65] += 2 * dqdci;           // dwdot[OH]/d[HO2]
    J[66] -= dqdci;               // dwdot[HO2]/d[HO2]
    // d()/dT
    J[93] -= dqdT;                // dwdot[H]/dT
    J[95] += 2 * dqdT;            // dwdot[OH]/dT
    J[96] -= dqdT;                // dwdot[HO2]/dT

    /*reaction 13: HO2 + O <=> O2 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    // forward
    phi_f = sc[4]*sc[6];
    k_f = 1e-06 * 32500000000000
                * exp(0 * tc[0] - 0.50321666580472 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580472 *  (0)  * invT2;
    // reverse
    phi_r = sc[1]*sc[5];
    Kc = exp(-g_RT[1] + g_RT[4] - g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[6]) + (h_RT[1] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    // rate of progress
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    // update wdot
    wdot[1] += q; // O2
    wdot[4] -= q; // O
    wdot[5] += q; // OH
    wdot[6] -= q; // HO2
    // d()/d[O2]
    dqdci =  - k_r*sc[5];
    J[11] += dqdci;               // dwdot[O2]/d[O2]
    J[14] -= dqdci;               // dwdot[O]/d[O2]
    J[15] += dqdci;               // dwdot[OH]/d[O2]
    J[16] -= dqdci;               // dwdot[HO2]/d[O2]
    // d()/d[O]
    dqdci =  + k_f*sc[6];
    J[41] += dqdci;               // dwdot[O2]/d[O]
    J[44] -= dqdci;               // dwdot[O]/d[O]
    J[45] += dqdci;               // dwdot[OH]/d[O]
    J[46] -= dqdci;               // dwdot[HO2]/d[O]
    // d()/d[OH]
    dqdci =  - k_r*sc[1];
    J[51] += dqdci;               // dwdot[O2]/d[OH]
    J[54] -= dqdci;               // dwdot[O]/d[OH]
    J[55] += dqdci;               // dwdot[OH]/d[OH]
    J[56] -= dqdci;               // dwdot[HO2]/d[OH]
    // d()/d[HO2]
    dqdci =  + k_f*sc[4];
    J[61] += dqdci;               // dwdot[O2]/d[HO2]
    J[64] -= dqdci;               // dwdot[O]/d[HO2]
    J[65] += dqdci;               // dwdot[OH]/d[HO2]
    J[66] -= dqdci;               // dwdot[HO2]/d[HO2]
    // d()/dT
    J[91] += dqdT;                // dwdot[O2]/dT
    J[94] -= dqdT;                // dwdot[O]/dT
    J[95] += dqdT;                // dwdot[OH]/dT
    J[96] -= dqdT;                // dwdot[HO2]/dT

    /*reaction 14: HO2 + OH <=> H2O + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    // forward
    phi_f = sc[5]*sc[6];
    k_f = 1e-06 * 28900000000000
                * exp(0 * tc[0] - 0.50321666580472 * (-497) * invT);
    dlnkfdT = 0 * invT + 0.50321666580472 *  (-497)  * invT2;
    // reverse
    phi_r = sc[1]*sc[2];
    Kc = exp(-g_RT[1] - g_RT[2] + g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[6]) + (h_RT[1] + h_RT[2]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    // rate of progress
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    // update wdot
    wdot[1] += q; // O2
    wdot[2] += q; // H2O
    wdot[5] -= q; // OH
    wdot[6] -= q; // HO2
    // d()/d[O2]
    dqdci =  - k_r*sc[2];
    J[11] += dqdci;               // dwdot[O2]/d[O2]
    J[12] += dqdci;               // dwdot[H2O]/d[O2]
    J[15] -= dqdci;               // dwdot[OH]/d[O2]
    J[16] -= dqdci;               // dwdot[HO2]/d[O2]
    // d()/d[H2O]
    dqdci =  - k_r*sc[1];
    J[21] += dqdci;               // dwdot[O2]/d[H2O]
    J[22] += dqdci;               // dwdot[H2O]/d[H2O]
    J[25] -= dqdci;               // dwdot[OH]/d[H2O]
    J[26] -= dqdci;               // dwdot[HO2]/d[H2O]
    // d()/d[OH]
    dqdci =  + k_f*sc[6];
    J[51] += dqdci;               // dwdot[O2]/d[OH]
    J[52] += dqdci;               // dwdot[H2O]/d[OH]
    J[55] -= dqdci;               // dwdot[OH]/d[OH]
    J[56] -= dqdci;               // dwdot[HO2]/d[OH]
    // d()/d[HO2]
    dqdci =  + k_f*sc[5];
    J[61] += dqdci;               // dwdot[O2]/d[HO2]
    J[62] += dqdci;               // dwdot[H2O]/d[HO2]
    J[65] -= dqdci;               // dwdot[OH]/d[HO2]
    J[66] -= dqdci;               // dwdot[HO2]/d[HO2]
    // d()/dT
    J[91] += dqdT;                // dwdot[O2]/dT
    J[92] += dqdT;                // dwdot[H2O]/dT
    J[95] -= dqdT;                // dwdot[OH]/dT
    J[96] -= dqdT;                // dwdot[HO2]/dT

    /*reaction 15: HO2 + HO2 <=> H2O2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    // forward
    phi_f = (sc[6] * sc[6]);
    k_f = 1e-06 * 420000000000000
                * exp(0 * tc[0] - 0.50321666580472 * (11982) * invT);
    dlnkfdT = 0 * invT + 0.50321666580472 *  (11982)  * invT2;
    // reverse
    phi_r = sc[1]*sc[7];
    Kc = exp(-g_RT[1] + g_RT[6] + g_RT[6] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[6]) + (h_RT[1] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    // rate of progress
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    // update wdot
    wdot[1] += q; // O2
    wdot[6] -= 2 * q; // HO2
    wdot[7] += q; // H2O2
    // d()/d[O2]
    dqdci =  - k_r*sc[7];
    J[11] += dqdci;               // dwdot[O2]/d[O2]
    J[16] += -2 * dqdci;          // dwdot[HO2]/d[O2]
    J[17] += dqdci;               // dwdot[H2O2]/d[O2]
    // d()/d[HO2]
    dqdci =  + k_f*2.000000*sc[6];
    J[61] += dqdci;               // dwdot[O2]/d[HO2]
    J[66] += -2 * dqdci;          // dwdot[HO2]/d[HO2]
    J[67] += dqdci;               // dwdot[H2O2]/d[HO2]
    // d()/d[H2O2]
    dqdci =  - k_r*sc[1];
    J[71] += dqdci;               // dwdot[O2]/d[H2O2]
    J[76] += -2 * dqdci;          // dwdot[HO2]/d[H2O2]
    J[77] += dqdci;               // dwdot[H2O2]/d[H2O2]
    // d()/dT
    J[91] += dqdT;                // dwdot[O2]/dT
    J[96] += -2 * dqdT;           // dwdot[HO2]/dT
    J[97] += dqdT;                // dwdot[H2O2]/dT

    /*reaction 16: HO2 + HO2 <=> H2O2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    // forward
    phi_f = (sc[6] * sc[6]);
    k_f = 1e-06 * 130000000000
                * exp(0 * tc[0] - 0.50321666580472 * (-1629.3) * invT);
    dlnkfdT = 0 * invT + 0.50321666580472 *  (-1629.3)  * invT2;
    // reverse
    phi_r = sc[1]*sc[7];
    Kc = exp(-g_RT[1] + g_RT[6] + g_RT[6] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[6]) + (h_RT[1] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    // rate of progress
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    // update wdot
    wdot[1] += q; // O2
    wdot[6] -= 2 * q; // HO2
    wdot[7] += q; // H2O2
    // d()/d[O2]
    dqdci =  - k_r*sc[7];
    J[11] += dqdci;               // dwdot[O2]/d[O2]
    J[16] += -2 * dqdci;          // dwdot[HO2]/d[O2]
    J[17] += dqdci;               // dwdot[H2O2]/d[O2]
    // d()/d[HO2]
    dqdci =  + k_f*2.000000*sc[6];
    J[61] += dqdci;               // dwdot[O2]/d[HO2]
    J[66] += -2 * dqdci;          // dwdot[HO2]/d[HO2]
    J[67] += dqdci;               // dwdot[H2O2]/d[HO2]
    // d()/d[H2O2]
    dqdci =  - k_r*sc[1];
    J[71] += dqdci;               // dwdot[O2]/d[H2O2]
    J[76] += -2 * dqdci;          // dwdot[HO2]/d[H2O2]
    J[77] += dqdci;               // dwdot[H2O2]/d[H2O2]
    // d()/dT
    J[91] += dqdT;                // dwdot[O2]/dT
    J[96] += -2 * dqdT;           // dwdot[HO2]/dT
    J[97] += dqdT;                // dwdot[H2O2]/dT

    /*reaction 17: H2O2 + H <=> H2O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    // forward
    phi_f = sc[3]*sc[7];
    k_f = 1e-06 * 24100000000000
                * exp(0 * tc[0] - 0.50321666580472 * (3970) * invT);
    dlnkfdT = 0 * invT + 0.50321666580472 *  (3970)  * invT2;
    // reverse
    phi_r = sc[2]*sc[5];
    Kc = exp(-g_RT[2] + g_RT[3] - g_RT[5] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[7]) + (h_RT[2] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    // rate of progress
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    // update wdot
    wdot[2] += q; // H2O
    wdot[3] -= q; // H
    wdot[5] += q; // OH
    wdot[7] -= q; // H2O2
    // d()/d[H2O]
    dqdci =  - k_r*sc[5];
    J[22] += dqdci;               // dwdot[H2O]/d[H2O]
    J[23] -= dqdci;               // dwdot[H]/d[H2O]
    J[25] += dqdci;               // dwdot[OH]/d[H2O]
    J[27] -= dqdci;               // dwdot[H2O2]/d[H2O]
    // d()/d[H]
    dqdci =  + k_f*sc[7];
    J[32] += dqdci;               // dwdot[H2O]/d[H]
    J[33] -= dqdci;               // dwdot[H]/d[H]
    J[35] += dqdci;               // dwdot[OH]/d[H]
    J[37] -= dqdci;               // dwdot[H2O2]/d[H]
    // d()/d[OH]
    dqdci =  - k_r*sc[2];
    J[52] += dqdci;               // dwdot[H2O]/d[OH]
    J[53] -= dqdci;               // dwdot[H]/d[OH]
    J[55] += dqdci;               // dwdot[OH]/d[OH]
    J[57] -= dqdci;               // dwdot[H2O2]/d[OH]
    // d()/d[H2O2]
    dqdci =  + k_f*sc[3];
    J[72] += dqdci;               // dwdot[H2O]/d[H2O2]
    J[73] -= dqdci;               // dwdot[H]/d[H2O2]
    J[75] += dqdci;               // dwdot[OH]/d[H2O2]
    J[77] -= dqdci;               // dwdot[H2O2]/d[H2O2]
    // d()/dT
    J[92] += dqdT;                // dwdot[H2O]/dT
    J[93] -= dqdT;                // dwdot[H]/dT
    J[95] += dqdT;                // dwdot[OH]/dT
    J[97] -= dqdT;                // dwdot[H2O2]/dT

    /*reaction 18: H2O2 + H <=> HO2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    // forward
    phi_f = sc[3]*sc[7];
    k_f = 1e-06 * 48200000000000
                * exp(0 * tc[0] - 0.50321666580472 * (7950) * invT);
    dlnkfdT = 0 * invT + 0.50321666580472 *  (7950)  * invT2;
    // reverse
    phi_r = sc[0]*sc[6];
    Kc = exp(-g_RT[0] + g_RT[3] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[7]) + (h_RT[0] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    // rate of progress
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    // update wdot
    wdot[0] += q; // H2
    wdot[3] -= q; // H
    wdot[6] += q; // HO2
    wdot[7] -= q; // H2O2
    // d()/d[H2]
    dqdci =  - k_r*sc[6];
    J[0] += dqdci;                // dwdot[H2]/d[H2]
    J[3] -= dqdci;                // dwdot[H]/d[H2]
    J[6] += dqdci;                // dwdot[HO2]/d[H2]
    J[7] -= dqdci;                // dwdot[H2O2]/d[H2]
    // d()/d[H]
    dqdci =  + k_f*sc[7];
    J[30] += dqdci;               // dwdot[H2]/d[H]
    J[33] -= dqdci;               // dwdot[H]/d[H]
    J[36] += dqdci;               // dwdot[HO2]/d[H]
    J[37] -= dqdci;               // dwdot[H2O2]/d[H]
    // d()/d[HO2]
    dqdci =  - k_r*sc[0];
    J[60] += dqdci;               // dwdot[H2]/d[HO2]
    J[63] -= dqdci;               // dwdot[H]/d[HO2]
    J[66] += dqdci;               // dwdot[HO2]/d[HO2]
    J[67] -= dqdci;               // dwdot[H2O2]/d[HO2]
    // d()/d[H2O2]
    dqdci =  + k_f*sc[3];
    J[70] += dqdci;               // dwdot[H2]/d[H2O2]
    J[73] -= dqdci;               // dwdot[H]/d[H2O2]
    J[76] += dqdci;               // dwdot[HO2]/d[H2O2]
    J[77] -= dqdci;               // dwdot[H2O2]/d[H2O2]
    // d()/dT
    J[90] += dqdT;                // dwdot[H2]/dT
    J[93] -= dqdT;                // dwdot[H]/dT
    J[96] += dqdT;                // dwdot[HO2]/dT
    J[97] -= dqdT;                // dwdot[H2O2]/dT

    /*reaction 19: H2O2 + O <=> OH + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    // forward
    phi_f = sc[4]*sc[7];
    k_f = 1e-06 * 9550000
                * exp(2 * tc[0] - 0.50321666580472 * (3970) * invT);
    dlnkfdT = 2 * invT + 0.50321666580472 *  (3970)  * invT2;
    // reverse
    phi_r = sc[5]*sc[6];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[7]) + (h_RT[5] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    // rate of progress
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    // update wdot
    wdot[4] -= q; // O
    wdot[5] += q; // OH
    wdot[6] += q; // HO2
    wdot[7] -= q; // H2O2
    // d()/d[O]
    dqdci =  + k_f*sc[7];
    J[44] -= dqdci;               // dwdot[O]/d[O]
    J[45] += dqdci;               // dwdot[OH]/d[O]
    J[46] += dqdci;               // dwdot[HO2]/d[O]
    J[47] -= dqdci;               // dwdot[H2O2]/d[O]
    // d()/d[OH]
    dqdci =  - k_r*sc[6];
    J[54] -= dqdci;               // dwdot[O]/d[OH]
    J[55] += dqdci;               // dwdot[OH]/d[OH]
    J[56] += dqdci;               // dwdot[HO2]/d[OH]
    J[57] -= dqdci;               // dwdot[H2O2]/d[OH]
    // d()/d[HO2]
    dqdci =  - k_r*sc[5];
    J[64] -= dqdci;               // dwdot[O]/d[HO2]
    J[65] += dqdci;               // dwdot[OH]/d[HO2]
    J[66] += dqdci;               // dwdot[HO2]/d[HO2]
    J[67] -= dqdci;               // dwdot[H2O2]/d[HO2]
    // d()/d[H2O2]
    dqdci =  + k_f*sc[4];
    J[74] -= dqdci;               // dwdot[O]/d[H2O2]
    J[75] += dqdci;               // dwdot[OH]/d[H2O2]
    J[76] += dqdci;               // dwdot[HO2]/d[H2O2]
    J[77] -= dqdci;               // dwdot[H2O2]/d[H2O2]
    // d()/dT
    J[94] -= dqdT;                // dwdot[O]/dT
    J[95] += dqdT;                // dwdot[OH]/dT
    J[96] += dqdT;                // dwdot[HO2]/dT
    J[97] -= dqdT;                // dwdot[H2O2]/dT

    /*reaction 20: H2O2 + OH <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    // forward
    phi_f = sc[5]*sc[7];
    k_f = 1e-06 * 1000000000000
                * exp(0 * tc[0] - 0.50321666580472 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580472 *  (0)  * invT2;
    // reverse
    phi_r = sc[2]*sc[6];
    Kc = exp(-g_RT[2] + g_RT[5] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[7]) + (h_RT[2] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    // rate of progress
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    // update wdot
    wdot[2] += q; // H2O
    wdot[5] -= q; // OH
    wdot[6] += q; // HO2
    wdot[7] -= q; // H2O2
    // d()/d[H2O]
    dqdci =  - k_r*sc[6];
    J[22] += dqdci;               // dwdot[H2O]/d[H2O]
    J[25] -= dqdci;               // dwdot[OH]/d[H2O]
    J[26] += dqdci;               // dwdot[HO2]/d[H2O]
    J[27] -= dqdci;               // dwdot[H2O2]/d[H2O]
    // d()/d[OH]
    dqdci =  + k_f*sc[7];
    J[52] += dqdci;               // dwdot[H2O]/d[OH]
    J[55] -= dqdci;               // dwdot[OH]/d[OH]
    J[56] += dqdci;               // dwdot[HO2]/d[OH]
    J[57] -= dqdci;               // dwdot[H2O2]/d[OH]
    // d()/d[HO2]
    dqdci =  - k_r*sc[2];
    J[62] += dqdci;               // dwdot[H2O]/d[HO2]
    J[65] -= dqdci;               // dwdot[OH]/d[HO2]
    J[66] += dqdci;               // dwdot[HO2]/d[HO2]
    J[67] -= dqdci;               // dwdot[H2O2]/d[HO2]
    // d()/d[H2O2]
    dqdci =  + k_f*sc[5];
    J[72] += dqdci;               // dwdot[H2O]/d[H2O2]
    J[75] -= dqdci;               // dwdot[OH]/d[H2O2]
    J[76] += dqdci;               // dwdot[HO2]/d[H2O2]
    J[77] -= dqdci;               // dwdot[H2O2]/d[H2O2]
    // d()/dT
    J[92] += dqdT;                // dwdot[H2O]/dT
    J[95] -= dqdT;                // dwdot[OH]/dT
    J[96] += dqdT;                // dwdot[HO2]/dT
    J[97] -= dqdT;                // dwdot[H2O2]/dT

    /*reaction 21: H2O2 + OH <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    // forward
    phi_f = sc[5]*sc[7];
    k_f = 1e-06 * 580000000000000
                * exp(0 * tc[0] - 0.50321666580472 * (9557) * invT);
    dlnkfdT = 0 * invT + 0.50321666580472 *  (9557)  * invT2;
    // reverse
    phi_r = sc[2]*sc[6];
    Kc = exp(-g_RT[2] + g_RT[5] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[7]) + (h_RT[2] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    // rate of progress
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    // update wdot
    wdot[2] += q; // H2O
    wdot[5] -= q; // OH
    wdot[6] += q; // HO2
    wdot[7] -= q; // H2O2
    // d()/d[H2O]
    dqdci =  - k_r*sc[6];
    J[22] += dqdci;               // dwdot[H2O]/d[H2O]
    J[25] -= dqdci;               // dwdot[OH]/d[H2O]
    J[26] += dqdci;               // dwdot[HO2]/d[H2O]
    J[27] -= dqdci;               // dwdot[H2O2]/d[H2O]
    // d()/d[OH]
    dqdci =  + k_f*sc[7];
    J[52] += dqdci;               // dwdot[H2O]/d[OH]
    J[55] -= dqdci;               // dwdot[OH]/d[OH]
    J[56] += dqdci;               // dwdot[HO2]/d[OH]
    J[57] -= dqdci;               // dwdot[H2O2]/d[OH]
    // d()/d[HO2]
    dqdci =  - k_r*sc[2];
    J[62] += dqdci;               // dwdot[H2O]/d[HO2]
    J[65] -= dqdci;               // dwdot[OH]/d[HO2]
    J[66] += dqdci;               // dwdot[HO2]/d[HO2]
    J[67] -= dqdci;               // dwdot[H2O2]/d[HO2]
    // d()/d[H2O2]
    dqdci =  + k_f*sc[5];
    J[72] += dqdci;               // dwdot[H2O]/d[H2O2]
    J[75] -= dqdci;               // dwdot[OH]/d[H2O2]
    J[76] += dqdci;               // dwdot[HO2]/d[H2O2]
    J[77] -= dqdci;               // dwdot[H2O2]/d[H2O2]
    // d()/dT
    J[92] += dqdT;                // dwdot[H2O]/dT
    J[95] -= dqdT;                // dwdot[OH]/dT
    J[96] += dqdT;                // dwdot[HO2]/dT
    J[97] -= dqdT;                // dwdot[H2O2]/dT

    amrex::Real c_R[9], dcRdT[9], e_RT[9];
    amrex::Real * eh_RT;
    if (HP) {
        cp_R(c_R, tc);
        dcvpRdT(dcRdT, tc);
        eh_RT = &h_RT[0];
    }
    else {
        cv_R(c_R, tc);
        dcvpRdT(dcRdT, tc);
        speciesInternalEnergy(e_RT, tc);
        eh_RT = &e_RT[0];
    }

    amrex::Real cmix = 0.0, ehmix = 0.0, dcmixdT=0.0, dehmixdT=0.0;
    for (int k = 0; k < 9; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[90+k];
    }

    amrex::Real cmixinv = 1.0/cmix;
    amrex::Real tmp1 = ehmix*cmixinv;
    amrex::Real tmp3 = cmixinv*T;
    amrex::Real tmp2 = tmp1*tmp3;
    amrex::Real dehmixdc;
    // dTdot/d[X]
    for (int k = 0; k < 9; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 9; ++m) {
            dehmixdc += eh_RT[m]*J[k*10+m];
        }
        J[k*10+9] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    // dTdot/dT
    J[99] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;
}

/*compute an approx to the reaction Jacobian (for preconditioning) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void DWDOT_SIMPLIFIED(amrex::Real *  J, amrex::Real *  sc, amrex::Real *  Tp, const int * HP)
{
    amrex::Real c[9];

    for (int k=0; k<9; k++) {
        c[k] = 1.e6 * sc[k];
    }

    aJacobian_precond(J, c, *Tp, *HP);

    // dwdot[k]/dT
    // dTdot/d[X]
    for (int k=0; k<9; k++) {
        J[90+k] *= 1.e-6;
        J[k*10+9] *= 1.e6;
    }

    return;
}

/*compute the reaction Jacobian on GPU */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void aJacobian(amrex::Real * J, amrex::Real * sc, amrex::Real T, const int consP)
{

    for (int i=0; i<100; i++) {
        J[i] = 0.0;
    }

    amrex::Real wdot[9];
    for (int k=0; k<9; k++) {
        wdot[k] = 0.0;
    }

    const amrex::Real tc[5] = { log(T), T, T*T, T*T*T, T*T*T*T }; // temperature cache
    amrex::Real invT = 1.0 / tc[1];
    amrex::Real invT2 = invT * invT;


    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    amrex::Real refC = 101325 / 8.31446 / T;
    amrex::Real refCinv = 1.0 / refC;

    /*compute the mixture concentration */
    amrex::Real mixture = 0.0;
    for (int k = 0; k < 9; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    amrex::Real g_RT[9];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    amrex::Real h_RT[9];
    speciesEnthalpy(h_RT, tc);

    amrex::Real phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    amrex::Real dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    amrex::Real dqdci, dcdc_fac, dqdc[9];
    amrex::Real Pr, fPr, F, k_0, logPr;
    amrex::Real logFcent, troe_c, troe_n, troePr_den, troePr, troe;
    amrex::Real Fcent1, Fcent2, Fcent3, Fcent;
    amrex::Real dlogFdc, dlogFdn, dlogFdcn_fac;
    amrex::Real dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
    const amrex::Real ln10 = log(10.0);
    const amrex::Real log10e = 1.0/log(10.0);
    /*reaction 1: H + O2 (+M) <=> HO2 (+M) */
    /*a pressure-fall-off reaction */
    // also 3-body
    // 3-body correction factor
    alpha = mixture + ( 2 - 1)*sc[0] + ( 11 - 1)*sc[2] + ( 0.78 - 1)*sc[1];
    // forward
    phi_f = sc[1]*sc[3];
    k_f = 1e-06 * 1475000000000
                * exp(0.6 * tc[0] - 0.50321666580472 * (0) * invT);
    dlnkfdT = 0.6 * invT + 0.50321666580472 *  0  * invT2;
    // pressure-fall-off
    k_0 = 6.366e+20 * exp(-1.72 * tc[0] - 0.50321666580472 * (524.8) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -1.72 * invT + 0.50321666580472 * (524.8) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    // Troe form
    logPr = log10(Pr);
    Fcent1 = (1.-0.8)*exp(-T/1e-30);
    Fcent2 = 0.8 * exp(-T/1e+30);
    Fcent3 = 0.;
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/1e-30
        -Fcent2/1e+30
    );
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    // reverse
    phi_r = sc[6];
    Kc = refCinv * exp(g_RT[1] + g_RT[3] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3]) + (h_RT[6]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    // rate of progress
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    // update wdot
    wdot[1] -= q; // O2
    wdot[3] -= q; // H
    wdot[6] += q; // HO2
    // for convenience
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        // d()/d[H2]
        dqdci = (2 - 1)*dcdc_fac;
        J[1] -= dqdci;                // dwdot[O2]/d[H2]
        J[3] -= dqdci;                // dwdot[H]/d[H2]
        J[6] += dqdci;                // dwdot[HO2]/d[H2]
        // d()/d[O2]
        dqdci = (0.78 - 1)*dcdc_fac + k_f*sc[3];
        J[11] -= dqdci;               // dwdot[O2]/d[O2]
        J[13] -= dqdci;               // dwdot[H]/d[O2]
        J[16] += dqdci;               // dwdot[HO2]/d[O2]
        // d()/d[H2O]
        dqdci = (11 - 1)*dcdc_fac;
        J[21] -= dqdci;               // dwdot[O2]/d[H2O]
        J[23] -= dqdci;               // dwdot[H]/d[H2O]
        J[26] += dqdci;               // dwdot[HO2]/d[H2O]
        // d()/d[H]
        dqdci =  + k_f*sc[1];
        J[31] -= dqdci;               // dwdot[O2]/d[H]
        J[33] -= dqdci;               // dwdot[H]/d[H]
        J[36] += dqdci;               // dwdot[HO2]/d[H]
        // d()/d[HO2]
        dqdci =  - k_r;
        J[61] -= dqdci;               // dwdot[O2]/d[HO2]
        J[63] -= dqdci;               // dwdot[H]/d[HO2]
        J[66] += dqdci;               // dwdot[HO2]/d[HO2]
    }
    else {
        dqdc[0] = 2*dcdc_fac;
        dqdc[1] = 0.78*dcdc_fac + k_f*sc[3];
        dqdc[2] = 11*dcdc_fac;
        dqdc[3] = dcdc_fac + k_f*sc[1];
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac;
        dqdc[6] = dcdc_fac - k_r;
        dqdc[7] = dcdc_fac;
        dqdc[8] = dcdc_fac;
        for (int k=0; k<9; k++) {
            J[10*k+1] -= dqdc[k];
            J[10*k+3] -= dqdc[k];
            J[10*k+6] += dqdc[k];
        }
    }
    J[91] -= dqdT; // dwdot[O2]/dT
    J[93] -= dqdT; // dwdot[H]/dT
    J[96] += dqdT; // dwdot[HO2]/dT

    /*reaction 2: H2O2 (+M) <=> OH + OH (+M) */
    /*a pressure-fall-off reaction */
    // also 3-body
    // 3-body correction factor
    alpha = mixture + ( 2.5 - 1)*sc[0] + ( 12 - 1)*sc[2];
    // forward
    phi_f = sc[7];
    k_f = 1 * 295100000000000
                * exp(0 * tc[0] - 0.50321666580472 * (48430) * invT);
    dlnkfdT = 0 * invT + 0.50321666580472 *  48430  * invT2;
    // pressure-fall-off
    k_0 = 1.202e+17 * exp(0 * tc[0] - 0.50321666580472 * (45500) * invT);
    Pr = 1e-6 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = 0 * invT + 0.50321666580472 * (45500) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    // Troe form
    logPr = log10(Pr);
    Fcent1 = (1.-0.5)*exp(-T/1e-30);
    Fcent2 = 0.5 * exp(-T/1e+30);
    Fcent3 = 0.;
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/1e-30
        -Fcent2/1e+30
    );
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    // reverse
    phi_r = (sc[5] * sc[5]);
    Kc = refC * exp(-g_RT[5] - g_RT[5] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[7]) + (2.000000*h_RT[5]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    // rate of progress
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    // update wdot
    wdot[5] += 2 * q; // OH
    wdot[7] -= q; // H2O2
    // for convenience
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        // d()/d[H2]
        dqdci = (2.5 - 1)*dcdc_fac;
        J[5] += 2 * dqdci;            // dwdot[OH]/d[H2]
        J[7] -= dqdci;                // dwdot[H2O2]/d[H2]
        // d()/d[H2O]
        dqdci = (12 - 1)*dcdc_fac;
        J[25] += 2 * dqdci;           // dwdot[OH]/d[H2O]
        J[27] -= dqdci;               // dwdot[H2O2]/d[H2O]
        // d()/d[OH]
        dqdci =  - k_r*2.000000*sc[5];
        J[55] += 2 * dqdci;           // dwdot[OH]/d[OH]
        J[57] -= dqdci;               // dwdot[H2O2]/d[OH]
        // d()/d[H2O2]
        dqdci =  + k_f;
        J[75] += 2 * dqdci;           // dwdot[OH]/d[H2O2]
        J[77] -= dqdci;               // dwdot[H2O2]/d[H2O2]
    }
    else {
        dqdc[0] = 2.5*dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = 12*dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac;
        dqdc[5] = dcdc_fac - k_r*2.000000*sc[5];
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac + k_f;
        dqdc[8] = dcdc_fac;
        for (int k=0; k<9; k++) {
            J[10*k+5] += 2 * dqdc[k];
            J[10*k+7] -= dqdc[k];
        }
    }
    J[95] += 2 * dqdT; // dwdot[OH]/dT
    J[97] -= dqdT; // dwdot[H2O2]/dT

    /*reaction 3: H2 + M <=> H + H + M */
    /*a third-body and non-pressure-fall-off reaction */
    // 3-body correction factor
    alpha = mixture + ( 2.5 - 1)*sc[0] + ( 12 - 1)*sc[2];
    // forward
    phi_f = sc[0];
    k_f = 1e-06 * 4.577e+19
                * exp(-1.4 * tc[0] - 0.50321666580472 * (104380) * invT);
    dlnkfdT = -1.4 * invT + 0.50321666580472 *  104380  * invT2;
    // reverse
    phi_r = (sc[3] * sc[3]);
    Kc = refC * exp(g_RT[0] - g_RT[3] - g_RT[3]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0]) + (2.000000*h_RT[3]) - 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    // rate of progress
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    // update wdot
    wdot[0] -= q; // H2
    wdot[3] += 2 * q; // H
    // for convenience
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        // d()/d[H2]
        dqdci = (2.5 - 1)*q_nocor + k_f;
        J[0] -= dqdci;                // dwdot[H2]/d[H2]
        J[3] += 2 * dqdci;            // dwdot[H]/d[H2]
        // d()/d[H2O]
        dqdci = (12 - 1)*q_nocor;
        J[20] -= dqdci;               // dwdot[H2]/d[H2O]
        J[23] += 2 * dqdci;           // dwdot[H]/d[H2O]
        // d()/d[H]
        dqdci =  - k_r*2.000000*sc[3];
        J[30] -= dqdci;               // dwdot[H2]/d[H]
        J[33] += 2 * dqdci;           // dwdot[H]/d[H]
    }
    else {
        dqdc[0] = 2.5*q_nocor + k_f;
        dqdc[1] = q_nocor;
        dqdc[2] = 12*q_nocor;
        dqdc[3] = q_nocor - k_r*2.000000*sc[3];
        dqdc[4] = q_nocor;
        dqdc[5] = q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        for (int k=0; k<9; k++) {
            J[10*k+0] -= dqdc[k];
            J[10*k+3] += 2 * dqdc[k];
        }
    }
    J[90] -= dqdT; // dwdot[H2]/dT
    J[93] += 2 * dqdT; // dwdot[H]/dT

    /*reaction 4: O + O + M <=> O2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    // 3-body correction factor
    alpha = mixture + ( 2.5 - 1)*sc[0] + ( 12 - 1)*sc[2];
    // forward
    phi_f = (sc[4] * sc[4]);
    k_f = 1e-12 * 6.165e+15
                * exp(-0.5 * tc[0] - 0.50321666580472 * (0) * invT);
    dlnkfdT = -0.5 * invT + 0.50321666580472 *  0  * invT2;
    // reverse
    phi_r = sc[1];
    Kc = refCinv * exp(-g_RT[1] + g_RT[4] + g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[4]) + (h_RT[1]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    // rate of progress
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    // update wdot
    wdot[1] += q; // O2
    wdot[4] -= 2 * q; // O
    // for convenience
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        // d()/d[H2]
        dqdci = (2.5 - 1)*q_nocor;
        J[1] += dqdci;                // dwdot[O2]/d[H2]
        J[4] += -2 * dqdci;           // dwdot[O]/d[H2]
        // d()/d[O2]
        dqdci =  - k_r;
        J[11] += dqdci;               // dwdot[O2]/d[O2]
        J[14] += -2 * dqdci;          // dwdot[O]/d[O2]
        // d()/d[H2O]
        dqdci = (12 - 1)*q_nocor;
        J[21] += dqdci;               // dwdot[O2]/d[H2O]
        J[24] += -2 * dqdci;          // dwdot[O]/d[H2O]
        // d()/d[O]
        dqdci =  + k_f*2.000000*sc[4];
        J[41] += dqdci;               // dwdot[O2]/d[O]
        J[44] += -2 * dqdci;          // dwdot[O]/d[O]
    }
    else {
        dqdc[0] = 2.5*q_nocor;
        dqdc[1] = q_nocor - k_r;
        dqdc[2] = 12*q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor + k_f*2.000000*sc[4];
        dqdc[5] = q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        for (int k=0; k<9; k++) {
            J[10*k+1] += dqdc[k];
            J[10*k+4] += -2 * dqdc[k];
        }
    }
    J[91] += dqdT; // dwdot[O2]/dT
    J[94] += -2 * dqdT; // dwdot[O]/dT

    /*reaction 5: O + H + M <=> OH + M */
    /*a third-body and non-pressure-fall-off reaction */
    // 3-body correction factor
    alpha = mixture + ( 2.5 - 1)*sc[0] + ( 12 - 1)*sc[2];
    // forward
    phi_f = sc[3]*sc[4];
    k_f = 1e-12 * 4.714e+18
                * exp(-1 * tc[0] - 0.50321666580472 * (0) * invT);
    dlnkfdT = -1 * invT + 0.50321666580472 *  0  * invT2;
    // reverse
    phi_r = sc[5];
    Kc = refCinv * exp(g_RT[3] + g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[4]) + (h_RT[5]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    // rate of progress
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    // update wdot
    wdot[3] -= q; // H
    wdot[4] -= q; // O
    wdot[5] += q; // OH
    // for convenience
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        // d()/d[H2]
        dqdci = (2.5 - 1)*q_nocor;
        J[3] -= dqdci;                // dwdot[H]/d[H2]
        J[4] -= dqdci;                // dwdot[O]/d[H2]
        J[5] += dqdci;                // dwdot[OH]/d[H2]
        // d()/d[H2O]
        dqdci = (12 - 1)*q_nocor;
        J[23] -= dqdci;               // dwdot[H]/d[H2O]
        J[24] -= dqdci;               // dwdot[O]/d[H2O]
        J[25] += dqdci;               // dwdot[OH]/d[H2O]
        // d()/d[H]
        dqdci =  + k_f*sc[4];
        J[33] -= dqdci;               // dwdot[H]/d[H]
        J[34] -= dqdci;               // dwdot[O]/d[H]
        J[35] += dqdci;               // dwdot[OH]/d[H]
        // d()/d[O]
        dqdci =  + k_f*sc[3];
        J[43] -= dqdci;               // dwdot[H]/d[O]
        J[44] -= dqdci;               // dwdot[O]/d[O]
        J[45] += dqdci;               // dwdot[OH]/d[O]
        // d()/d[OH]
        dqdci =  - k_r;
        J[53] -= dqdci;               // dwdot[H]/d[OH]
        J[54] -= dqdci;               // dwdot[O]/d[OH]
        J[55] += dqdci;               // dwdot[OH]/d[OH]
    }
    else {
        dqdc[0] = 2.5*q_nocor;
        dqdc[1] = q_nocor;
        dqdc[2] = 12*q_nocor;
        dqdc[3] = q_nocor + k_f*sc[4];
        dqdc[4] = q_nocor + k_f*sc[3];
        dqdc[5] = q_nocor - k_r;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        for (int k=0; k<9; k++) {
            J[10*k+3] -= dqdc[k];
            J[10*k+4] -= dqdc[k];
            J[10*k+5] += dqdc[k];
        }
    }
    J[93] -= dqdT; // dwdot[H]/dT
    J[94] -= dqdT; // dwdot[O]/dT
    J[95] += dqdT; // dwdot[OH]/dT

    /*reaction 6: H + OH + M <=> H2O + M */
    /*a third-body and non-pressure-fall-off reaction */
    // 3-body correction factor
    alpha = mixture + ( 2.5 - 1)*sc[0] + ( 12 - 1)*sc[2];
    // forward
    phi_f = sc[3]*sc[5];
    k_f = 1e-12 * 3.8e+22
                * exp(-2 * tc[0] - 0.50321666580472 * (0) * invT);
    dlnkfdT = -2 * invT + 0.50321666580472 *  0  * invT2;
    // reverse
    phi_r = sc[2];
    Kc = refCinv * exp(-g_RT[2] + g_RT[3] + g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[5]) + (h_RT[2]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    // rate of progress
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    // update wdot
    wdot[2] += q; // H2O
    wdot[3] -= q; // H
    wdot[5] -= q; // OH
    // for convenience
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        // d()/d[H2]
        dqdci = (2.5 - 1)*q_nocor;
        J[2] += dqdci;                // dwdot[H2O]/d[H2]
        J[3] -= dqdci;                // dwdot[H]/d[H2]
        J[5] -= dqdci;                // dwdot[OH]/d[H2]
        // d()/d[H2O]
        dqdci = (12 - 1)*q_nocor - k_r;
        J[22] += dqdci;               // dwdot[H2O]/d[H2O]
        J[23] -= dqdci;               // dwdot[H]/d[H2O]
        J[25] -= dqdci;               // dwdot[OH]/d[H2O]
        // d()/d[H]
        dqdci =  + k_f*sc[5];
        J[32] += dqdci;               // dwdot[H2O]/d[H]
        J[33] -= dqdci;               // dwdot[H]/d[H]
        J[35] -= dqdci;               // dwdot[OH]/d[H]
        // d()/d[OH]
        dqdci =  + k_f*sc[3];
        J[52] += dqdci;               // dwdot[H2O]/d[OH]
        J[53] -= dqdci;               // dwdot[H]/d[OH]
        J[55] -= dqdci;               // dwdot[OH]/d[OH]
    }
    else {
        dqdc[0] = 2.5*q_nocor;
        dqdc[1] = q_nocor;
        dqdc[2] = 12*q_nocor - k_r;
        dqdc[3] = q_nocor + k_f*sc[5];
        dqdc[4] = q_nocor;
        dqdc[5] = q_nocor + k_f*sc[3];
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        for (int k=0; k<9; k++) {
            J[10*k+2] += dqdc[k];
            J[10*k+3] -= dqdc[k];
            J[10*k+5] -= dqdc[k];
        }
    }
    J[92] += dqdT; // dwdot[H2O]/dT
    J[93] -= dqdT; // dwdot[H]/dT
    J[95] -= dqdT; // dwdot[OH]/dT

    /*reaction 7: H + O2 <=> O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    // forward
    phi_f = sc[1]*sc[3];
    k_f = 1e-06 * 3.547e+15
                * exp(-0.406 * tc[0] - 0.50321666580472 * (16599) * invT);
    dlnkfdT = -0.406 * invT + 0.50321666580472 *  16599  * invT2;
    // reverse
    phi_r = sc[4]*sc[5];
    Kc = exp(g_RT[1] + g_RT[3] - g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3]) + (h_RT[4] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    // rate of progress
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    // update wdot
    wdot[1] -= q; // O2
    wdot[3] -= q; // H
    wdot[4] += q; // O
    wdot[5] += q; // OH
    // d()/d[O2]
    dqdci =  + k_f*sc[3];
    J[11] -= dqdci;               // dwdot[O2]/d[O2]
    J[13] -= dqdci;               // dwdot[H]/d[O2]
    J[14] += dqdci;               // dwdot[O]/d[O2]
    J[15] += dqdci;               // dwdot[OH]/d[O2]
    // d()/d[H]
    dqdci =  + k_f*sc[1];
    J[31] -= dqdci;               // dwdot[O2]/d[H]
    J[33] -= dqdci;               // dwdot[H]/d[H]
    J[34] += dqdci;               // dwdot[O]/d[H]
    J[35] += dqdci;               // dwdot[OH]/d[H]
    // d()/d[O]
    dqdci =  - k_r*sc[5];
    J[41] -= dqdci;               // dwdot[O2]/d[O]
    J[43] -= dqdci;               // dwdot[H]/d[O]
    J[44] += dqdci;               // dwdot[O]/d[O]
    J[45] += dqdci;               // dwdot[OH]/d[O]
    // d()/d[OH]
    dqdci =  - k_r*sc[4];
    J[51] -= dqdci;               // dwdot[O2]/d[OH]
    J[53] -= dqdci;               // dwdot[H]/d[OH]
    J[54] += dqdci;               // dwdot[O]/d[OH]
    J[55] += dqdci;               // dwdot[OH]/d[OH]
    // d()/dT
    J[91] -= dqdT;                // dwdot[O2]/dT
    J[93] -= dqdT;                // dwdot[H]/dT
    J[94] += dqdT;                // dwdot[O]/dT
    J[95] += dqdT;                // dwdot[OH]/dT

    /*reaction 8: O + H2 <=> H + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    // forward
    phi_f = sc[0]*sc[4];
    k_f = 1e-06 * 50800
                * exp(2.67 * tc[0] - 0.50321666580472 * (6290) * invT);
    dlnkfdT = 2.67 * invT + 0.50321666580472 *  6290  * invT2;
    // reverse
    phi_r = sc[3]*sc[5];
    Kc = exp(g_RT[0] - g_RT[3] + g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[4]) + (h_RT[3] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    // rate of progress
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    // update wdot
    wdot[0] -= q; // H2
    wdot[3] += q; // H
    wdot[4] -= q; // O
    wdot[5] += q; // OH
    // d()/d[H2]
    dqdci =  + k_f*sc[4];
    J[0] -= dqdci;                // dwdot[H2]/d[H2]
    J[3] += dqdci;                // dwdot[H]/d[H2]
    J[4] -= dqdci;                // dwdot[O]/d[H2]
    J[5] += dqdci;                // dwdot[OH]/d[H2]
    // d()/d[H]
    dqdci =  - k_r*sc[5];
    J[30] -= dqdci;               // dwdot[H2]/d[H]
    J[33] += dqdci;               // dwdot[H]/d[H]
    J[34] -= dqdci;               // dwdot[O]/d[H]
    J[35] += dqdci;               // dwdot[OH]/d[H]
    // d()/d[O]
    dqdci =  + k_f*sc[0];
    J[40] -= dqdci;               // dwdot[H2]/d[O]
    J[43] += dqdci;               // dwdot[H]/d[O]
    J[44] -= dqdci;               // dwdot[O]/d[O]
    J[45] += dqdci;               // dwdot[OH]/d[O]
    // d()/d[OH]
    dqdci =  - k_r*sc[3];
    J[50] -= dqdci;               // dwdot[H2]/d[OH]
    J[53] += dqdci;               // dwdot[H]/d[OH]
    J[54] -= dqdci;               // dwdot[O]/d[OH]
    J[55] += dqdci;               // dwdot[OH]/d[OH]
    // d()/dT
    J[90] -= dqdT;                // dwdot[H2]/dT
    J[93] += dqdT;                // dwdot[H]/dT
    J[94] -= dqdT;                // dwdot[O]/dT
    J[95] += dqdT;                // dwdot[OH]/dT

    /*reaction 9: H2 + OH <=> H2O + H */
    /*a non-third-body and non-pressure-fall-off reaction */
    // forward
    phi_f = sc[0]*sc[5];
    k_f = 1e-06 * 216000000
                * exp(1.51 * tc[0] - 0.50321666580472 * (3430) * invT);
    dlnkfdT = 1.51 * invT + 0.50321666580472 *  3430  * invT2;
    // reverse
    phi_r = sc[2]*sc[3];
    Kc = exp(g_RT[0] - g_RT[2] - g_RT[3] + g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[5]) + (h_RT[2] + h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    // rate of progress
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    // update wdot
    wdot[0] -= q; // H2
    wdot[2] += q; // H2O
    wdot[3] += q; // H
    wdot[5] -= q; // OH
    // d()/d[H2]
    dqdci =  + k_f*sc[5];
    J[0] -= dqdci;                // dwdot[H2]/d[H2]
    J[2] += dqdci;                // dwdot[H2O]/d[H2]
    J[3] += dqdci;                // dwdot[H]/d[H2]
    J[5] -= dqdci;                // dwdot[OH]/d[H2]
    // d()/d[H2O]
    dqdci =  - k_r*sc[3];
    J[20] -= dqdci;               // dwdot[H2]/d[H2O]
    J[22] += dqdci;               // dwdot[H2O]/d[H2O]
    J[23] += dqdci;               // dwdot[H]/d[H2O]
    J[25] -= dqdci;               // dwdot[OH]/d[H2O]
    // d()/d[H]
    dqdci =  - k_r*sc[2];
    J[30] -= dqdci;               // dwdot[H2]/d[H]
    J[32] += dqdci;               // dwdot[H2O]/d[H]
    J[33] += dqdci;               // dwdot[H]/d[H]
    J[35] -= dqdci;               // dwdot[OH]/d[H]
    // d()/d[OH]
    dqdci =  + k_f*sc[0];
    J[50] -= dqdci;               // dwdot[H2]/d[OH]
    J[52] += dqdci;               // dwdot[H2O]/d[OH]
    J[53] += dqdci;               // dwdot[H]/d[OH]
    J[55] -= dqdci;               // dwdot[OH]/d[OH]
    // d()/dT
    J[90] -= dqdT;                // dwdot[H2]/dT
    J[92] += dqdT;                // dwdot[H2O]/dT
    J[93] += dqdT;                // dwdot[H]/dT
    J[95] -= dqdT;                // dwdot[OH]/dT

    /*reaction 10: O + H2O <=> OH + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    // forward
    phi_f = sc[2]*sc[4];
    k_f = 1e-06 * 2970000
                * exp(2.02 * tc[0] - 0.50321666580472 * (13400) * invT);
    dlnkfdT = 2.02 * invT + 0.50321666580472 *  13400  * invT2;
    // reverse
    phi_r = (sc[5] * sc[5]);
    Kc = exp(g_RT[2] + g_RT[4] - g_RT[5] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[4]) + (2.000000*h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    // rate of progress
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    // update wdot
    wdot[2] -= q; // H2O
    wdot[4] -= q; // O
    wdot[5] += 2 * q; // OH
    // d()/d[H2O]
    dqdci =  + k_f*sc[4];
    J[22] -= dqdci;               // dwdot[H2O]/d[H2O]
    J[24] -= dqdci;               // dwdot[O]/d[H2O]
    J[25] += 2 * dqdci;           // dwdot[OH]/d[H2O]
    // d()/d[O]
    dqdci =  + k_f*sc[2];
    J[42] -= dqdci;               // dwdot[H2O]/d[O]
    J[44] -= dqdci;               // dwdot[O]/d[O]
    J[45] += 2 * dqdci;           // dwdot[OH]/d[O]
    // d()/d[OH]
    dqdci =  - k_r*2.000000*sc[5];
    J[52] -= dqdci;               // dwdot[H2O]/d[OH]
    J[54] -= dqdci;               // dwdot[O]/d[OH]
    J[55] += 2 * dqdci;           // dwdot[OH]/d[OH]
    // d()/dT
    J[92] -= dqdT;                // dwdot[H2O]/dT
    J[94] -= dqdT;                // dwdot[O]/dT
    J[95] += 2 * dqdT;            // dwdot[OH]/dT

    /*reaction 11: HO2 + H <=> H2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    // forward
    phi_f = sc[3]*sc[6];
    k_f = 1e-06 * 16600000000000
                * exp(0 * tc[0] - 0.50321666580472 * (823) * invT);
    dlnkfdT = 0 * invT + 0.50321666580472 *  823  * invT2;
    // reverse
    phi_r = sc[0]*sc[1];
    Kc = exp(-g_RT[0] - g_RT[1] + g_RT[3] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[6]) + (h_RT[0] + h_RT[1]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    // rate of progress
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    // update wdot
    wdot[0] += q; // H2
    wdot[1] += q; // O2
    wdot[3] -= q; // H
    wdot[6] -= q; // HO2
    // d()/d[H2]
    dqdci =  - k_r*sc[1];
    J[0] += dqdci;                // dwdot[H2]/d[H2]
    J[1] += dqdci;                // dwdot[O2]/d[H2]
    J[3] -= dqdci;                // dwdot[H]/d[H2]
    J[6] -= dqdci;                // dwdot[HO2]/d[H2]
    // d()/d[O2]
    dqdci =  - k_r*sc[0];
    J[10] += dqdci;               // dwdot[H2]/d[O2]
    J[11] += dqdci;               // dwdot[O2]/d[O2]
    J[13] -= dqdci;               // dwdot[H]/d[O2]
    J[16] -= dqdci;               // dwdot[HO2]/d[O2]
    // d()/d[H]
    dqdci =  + k_f*sc[6];
    J[30] += dqdci;               // dwdot[H2]/d[H]
    J[31] += dqdci;               // dwdot[O2]/d[H]
    J[33] -= dqdci;               // dwdot[H]/d[H]
    J[36] -= dqdci;               // dwdot[HO2]/d[H]
    // d()/d[HO2]
    dqdci =  + k_f*sc[3];
    J[60] += dqdci;               // dwdot[H2]/d[HO2]
    J[61] += dqdci;               // dwdot[O2]/d[HO2]
    J[63] -= dqdci;               // dwdot[H]/d[HO2]
    J[66] -= dqdci;               // dwdot[HO2]/d[HO2]
    // d()/dT
    J[90] += dqdT;                // dwdot[H2]/dT
    J[91] += dqdT;                // dwdot[O2]/dT
    J[93] -= dqdT;                // dwdot[H]/dT
    J[96] -= dqdT;                // dwdot[HO2]/dT

    /*reaction 12: HO2 + H <=> OH + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    // forward
    phi_f = sc[3]*sc[6];
    k_f = 1e-06 * 70790000000000
                * exp(0 * tc[0] - 0.50321666580472 * (295) * invT);
    dlnkfdT = 0 * invT + 0.50321666580472 *  295  * invT2;
    // reverse
    phi_r = (sc[5] * sc[5]);
    Kc = exp(g_RT[3] - g_RT[5] - g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[6]) + (2.000000*h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    // rate of progress
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    // update wdot
    wdot[3] -= q; // H
    wdot[5] += 2 * q; // OH
    wdot[6] -= q; // HO2
    // d()/d[H]
    dqdci =  + k_f*sc[6];
    J[33] -= dqdci;               // dwdot[H]/d[H]
    J[35] += 2 * dqdci;           // dwdot[OH]/d[H]
    J[36] -= dqdci;               // dwdot[HO2]/d[H]
    // d()/d[OH]
    dqdci =  - k_r*2.000000*sc[5];
    J[53] -= dqdci;               // dwdot[H]/d[OH]
    J[55] += 2 * dqdci;           // dwdot[OH]/d[OH]
    J[56] -= dqdci;               // dwdot[HO2]/d[OH]
    // d()/d[HO2]
    dqdci =  + k_f*sc[3];
    J[63] -= dqdci;               // dwdot[H]/d[HO2]
    J[65] += 2 * dqdci;           // dwdot[OH]/d[HO2]
    J[66] -= dqdci;               // dwdot[HO2]/d[HO2]
    // d()/dT
    J[93] -= dqdT;                // dwdot[H]/dT
    J[95] += 2 * dqdT;            // dwdot[OH]/dT
    J[96] -= dqdT;                // dwdot[HO2]/dT

    /*reaction 13: HO2 + O <=> O2 + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    // forward
    phi_f = sc[4]*sc[6];
    k_f = 1e-06 * 32500000000000
                * exp(0 * tc[0] - 0.50321666580472 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580472 *  0  * invT2;
    // reverse
    phi_r = sc[1]*sc[5];
    Kc = exp(-g_RT[1] + g_RT[4] - g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[6]) + (h_RT[1] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    // rate of progress
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    // update wdot
    wdot[1] += q; // O2
    wdot[4] -= q; // O
    wdot[5] += q; // OH
    wdot[6] -= q; // HO2
    // d()/d[O2]
    dqdci =  - k_r*sc[5];
    J[11] += dqdci;               // dwdot[O2]/d[O2]
    J[14] -= dqdci;               // dwdot[O]/d[O2]
    J[15] += dqdci;               // dwdot[OH]/d[O2]
    J[16] -= dqdci;               // dwdot[HO2]/d[O2]
    // d()/d[O]
    dqdci =  + k_f*sc[6];
    J[41] += dqdci;               // dwdot[O2]/d[O]
    J[44] -= dqdci;               // dwdot[O]/d[O]
    J[45] += dqdci;               // dwdot[OH]/d[O]
    J[46] -= dqdci;               // dwdot[HO2]/d[O]
    // d()/d[OH]
    dqdci =  - k_r*sc[1];
    J[51] += dqdci;               // dwdot[O2]/d[OH]
    J[54] -= dqdci;               // dwdot[O]/d[OH]
    J[55] += dqdci;               // dwdot[OH]/d[OH]
    J[56] -= dqdci;               // dwdot[HO2]/d[OH]
    // d()/d[HO2]
    dqdci =  + k_f*sc[4];
    J[61] += dqdci;               // dwdot[O2]/d[HO2]
    J[64] -= dqdci;               // dwdot[O]/d[HO2]
    J[65] += dqdci;               // dwdot[OH]/d[HO2]
    J[66] -= dqdci;               // dwdot[HO2]/d[HO2]
    // d()/dT
    J[91] += dqdT;                // dwdot[O2]/dT
    J[94] -= dqdT;                // dwdot[O]/dT
    J[95] += dqdT;                // dwdot[OH]/dT
    J[96] -= dqdT;                // dwdot[HO2]/dT

    /*reaction 14: HO2 + OH <=> H2O + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    // forward
    phi_f = sc[5]*sc[6];
    k_f = 1e-06 * 28900000000000
                * exp(0 * tc[0] - 0.50321666580472 * (-497) * invT);
    dlnkfdT = 0 * invT + 0.50321666580472 *  -497  * invT2;
    // reverse
    phi_r = sc[1]*sc[2];
    Kc = exp(-g_RT[1] - g_RT[2] + g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[6]) + (h_RT[1] + h_RT[2]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    // rate of progress
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    // update wdot
    wdot[1] += q; // O2
    wdot[2] += q; // H2O
    wdot[5] -= q; // OH
    wdot[6] -= q; // HO2
    // d()/d[O2]
    dqdci =  - k_r*sc[2];
    J[11] += dqdci;               // dwdot[O2]/d[O2]
    J[12] += dqdci;               // dwdot[H2O]/d[O2]
    J[15] -= dqdci;               // dwdot[OH]/d[O2]
    J[16] -= dqdci;               // dwdot[HO2]/d[O2]
    // d()/d[H2O]
    dqdci =  - k_r*sc[1];
    J[21] += dqdci;               // dwdot[O2]/d[H2O]
    J[22] += dqdci;               // dwdot[H2O]/d[H2O]
    J[25] -= dqdci;               // dwdot[OH]/d[H2O]
    J[26] -= dqdci;               // dwdot[HO2]/d[H2O]
    // d()/d[OH]
    dqdci =  + k_f*sc[6];
    J[51] += dqdci;               // dwdot[O2]/d[OH]
    J[52] += dqdci;               // dwdot[H2O]/d[OH]
    J[55] -= dqdci;               // dwdot[OH]/d[OH]
    J[56] -= dqdci;               // dwdot[HO2]/d[OH]
    // d()/d[HO2]
    dqdci =  + k_f*sc[5];
    J[61] += dqdci;               // dwdot[O2]/d[HO2]
    J[62] += dqdci;               // dwdot[H2O]/d[HO2]
    J[65] -= dqdci;               // dwdot[OH]/d[HO2]
    J[66] -= dqdci;               // dwdot[HO2]/d[HO2]
    // d()/dT
    J[91] += dqdT;                // dwdot[O2]/dT
    J[92] += dqdT;                // dwdot[H2O]/dT
    J[95] -= dqdT;                // dwdot[OH]/dT
    J[96] -= dqdT;                // dwdot[HO2]/dT

    /*reaction 15: HO2 + HO2 <=> H2O2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    // forward
    phi_f = (sc[6] * sc[6]);
    k_f = 1e-06 * 420000000000000
                * exp(0 * tc[0] - 0.50321666580472 * (11982) * invT);
    dlnkfdT = 0 * invT + 0.50321666580472 *  11982  * invT2;
    // reverse
    phi_r = sc[1]*sc[7];
    Kc = exp(-g_RT[1] + g_RT[6] + g_RT[6] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[6]) + (h_RT[1] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    // rate of progress
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    // update wdot
    wdot[1] += q; // O2
    wdot[6] -= 2 * q; // HO2
    wdot[7] += q; // H2O2
    // d()/d[O2]
    dqdci =  - k_r*sc[7];
    J[11] += dqdci;               // dwdot[O2]/d[O2]
    J[16] += -2 * dqdci;          // dwdot[HO2]/d[O2]
    J[17] += dqdci;               // dwdot[H2O2]/d[O2]
    // d()/d[HO2]
    dqdci =  + k_f*2.000000*sc[6];
    J[61] += dqdci;               // dwdot[O2]/d[HO2]
    J[66] += -2 * dqdci;          // dwdot[HO2]/d[HO2]
    J[67] += dqdci;               // dwdot[H2O2]/d[HO2]
    // d()/d[H2O2]
    dqdci =  - k_r*sc[1];
    J[71] += dqdci;               // dwdot[O2]/d[H2O2]
    J[76] += -2 * dqdci;          // dwdot[HO2]/d[H2O2]
    J[77] += dqdci;               // dwdot[H2O2]/d[H2O2]
    // d()/dT
    J[91] += dqdT;                // dwdot[O2]/dT
    J[96] += -2 * dqdT;           // dwdot[HO2]/dT
    J[97] += dqdT;                // dwdot[H2O2]/dT

    /*reaction 16: HO2 + HO2 <=> H2O2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    // forward
    phi_f = (sc[6] * sc[6]);
    k_f = 1e-06 * 130000000000
                * exp(0 * tc[0] - 0.50321666580472 * (-1629.3) * invT);
    dlnkfdT = 0 * invT + 0.50321666580472 *  -1629.3  * invT2;
    // reverse
    phi_r = sc[1]*sc[7];
    Kc = exp(-g_RT[1] + g_RT[6] + g_RT[6] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[6]) + (h_RT[1] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    // rate of progress
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    // update wdot
    wdot[1] += q; // O2
    wdot[6] -= 2 * q; // HO2
    wdot[7] += q; // H2O2
    // d()/d[O2]
    dqdci =  - k_r*sc[7];
    J[11] += dqdci;               // dwdot[O2]/d[O2]
    J[16] += -2 * dqdci;          // dwdot[HO2]/d[O2]
    J[17] += dqdci;               // dwdot[H2O2]/d[O2]
    // d()/d[HO2]
    dqdci =  + k_f*2.000000*sc[6];
    J[61] += dqdci;               // dwdot[O2]/d[HO2]
    J[66] += -2 * dqdci;          // dwdot[HO2]/d[HO2]
    J[67] += dqdci;               // dwdot[H2O2]/d[HO2]
    // d()/d[H2O2]
    dqdci =  - k_r*sc[1];
    J[71] += dqdci;               // dwdot[O2]/d[H2O2]
    J[76] += -2 * dqdci;          // dwdot[HO2]/d[H2O2]
    J[77] += dqdci;               // dwdot[H2O2]/d[H2O2]
    // d()/dT
    J[91] += dqdT;                // dwdot[O2]/dT
    J[96] += -2 * dqdT;           // dwdot[HO2]/dT
    J[97] += dqdT;                // dwdot[H2O2]/dT

    /*reaction 17: H2O2 + H <=> H2O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    // forward
    phi_f = sc[3]*sc[7];
    k_f = 1e-06 * 24100000000000
                * exp(0 * tc[0] - 0.50321666580472 * (3970) * invT);
    dlnkfdT = 0 * invT + 0.50321666580472 *  3970  * invT2;
    // reverse
    phi_r = sc[2]*sc[5];
    Kc = exp(-g_RT[2] + g_RT[3] - g_RT[5] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[7]) + (h_RT[2] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    // rate of progress
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    // update wdot
    wdot[2] += q; // H2O
    wdot[3] -= q; // H
    wdot[5] += q; // OH
    wdot[7] -= q; // H2O2
    // d()/d[H2O]
    dqdci =  - k_r*sc[5];
    J[22] += dqdci;               // dwdot[H2O]/d[H2O]
    J[23] -= dqdci;               // dwdot[H]/d[H2O]
    J[25] += dqdci;               // dwdot[OH]/d[H2O]
    J[27] -= dqdci;               // dwdot[H2O2]/d[H2O]
    // d()/d[H]
    dqdci =  + k_f*sc[7];
    J[32] += dqdci;               // dwdot[H2O]/d[H]
    J[33] -= dqdci;               // dwdot[H]/d[H]
    J[35] += dqdci;               // dwdot[OH]/d[H]
    J[37] -= dqdci;               // dwdot[H2O2]/d[H]
    // d()/d[OH]
    dqdci =  - k_r*sc[2];
    J[52] += dqdci;               // dwdot[H2O]/d[OH]
    J[53] -= dqdci;               // dwdot[H]/d[OH]
    J[55] += dqdci;               // dwdot[OH]/d[OH]
    J[57] -= dqdci;               // dwdot[H2O2]/d[OH]
    // d()/d[H2O2]
    dqdci =  + k_f*sc[3];
    J[72] += dqdci;               // dwdot[H2O]/d[H2O2]
    J[73] -= dqdci;               // dwdot[H]/d[H2O2]
    J[75] += dqdci;               // dwdot[OH]/d[H2O2]
    J[77] -= dqdci;               // dwdot[H2O2]/d[H2O2]
    // d()/dT
    J[92] += dqdT;                // dwdot[H2O]/dT
    J[93] -= dqdT;                // dwdot[H]/dT
    J[95] += dqdT;                // dwdot[OH]/dT
    J[97] -= dqdT;                // dwdot[H2O2]/dT

    /*reaction 18: H2O2 + H <=> HO2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    // forward
    phi_f = sc[3]*sc[7];
    k_f = 1e-06 * 48200000000000
                * exp(0 * tc[0] - 0.50321666580472 * (7950) * invT);
    dlnkfdT = 0 * invT + 0.50321666580472 *  7950  * invT2;
    // reverse
    phi_r = sc[0]*sc[6];
    Kc = exp(-g_RT[0] + g_RT[3] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[3] + h_RT[7]) + (h_RT[0] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    // rate of progress
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    // update wdot
    wdot[0] += q; // H2
    wdot[3] -= q; // H
    wdot[6] += q; // HO2
    wdot[7] -= q; // H2O2
    // d()/d[H2]
    dqdci =  - k_r*sc[6];
    J[0] += dqdci;                // dwdot[H2]/d[H2]
    J[3] -= dqdci;                // dwdot[H]/d[H2]
    J[6] += dqdci;                // dwdot[HO2]/d[H2]
    J[7] -= dqdci;                // dwdot[H2O2]/d[H2]
    // d()/d[H]
    dqdci =  + k_f*sc[7];
    J[30] += dqdci;               // dwdot[H2]/d[H]
    J[33] -= dqdci;               // dwdot[H]/d[H]
    J[36] += dqdci;               // dwdot[HO2]/d[H]
    J[37] -= dqdci;               // dwdot[H2O2]/d[H]
    // d()/d[HO2]
    dqdci =  - k_r*sc[0];
    J[60] += dqdci;               // dwdot[H2]/d[HO2]
    J[63] -= dqdci;               // dwdot[H]/d[HO2]
    J[66] += dqdci;               // dwdot[HO2]/d[HO2]
    J[67] -= dqdci;               // dwdot[H2O2]/d[HO2]
    // d()/d[H2O2]
    dqdci =  + k_f*sc[3];
    J[70] += dqdci;               // dwdot[H2]/d[H2O2]
    J[73] -= dqdci;               // dwdot[H]/d[H2O2]
    J[76] += dqdci;               // dwdot[HO2]/d[H2O2]
    J[77] -= dqdci;               // dwdot[H2O2]/d[H2O2]
    // d()/dT
    J[90] += dqdT;                // dwdot[H2]/dT
    J[93] -= dqdT;                // dwdot[H]/dT
    J[96] += dqdT;                // dwdot[HO2]/dT
    J[97] -= dqdT;                // dwdot[H2O2]/dT

    /*reaction 19: H2O2 + O <=> OH + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    // forward
    phi_f = sc[4]*sc[7];
    k_f = 1e-06 * 9550000
                * exp(2 * tc[0] - 0.50321666580472 * (3970) * invT);
    dlnkfdT = 2 * invT + 0.50321666580472 *  3970  * invT2;
    // reverse
    phi_r = sc[5]*sc[6];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[7]) + (h_RT[5] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    // rate of progress
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    // update wdot
    wdot[4] -= q; // O
    wdot[5] += q; // OH
    wdot[6] += q; // HO2
    wdot[7] -= q; // H2O2
    // d()/d[O]
    dqdci =  + k_f*sc[7];
    J[44] -= dqdci;               // dwdot[O]/d[O]
    J[45] += dqdci;               // dwdot[OH]/d[O]
    J[46] += dqdci;               // dwdot[HO2]/d[O]
    J[47] -= dqdci;               // dwdot[H2O2]/d[O]
    // d()/d[OH]
    dqdci =  - k_r*sc[6];
    J[54] -= dqdci;               // dwdot[O]/d[OH]
    J[55] += dqdci;               // dwdot[OH]/d[OH]
    J[56] += dqdci;               // dwdot[HO2]/d[OH]
    J[57] -= dqdci;               // dwdot[H2O2]/d[OH]
    // d()/d[HO2]
    dqdci =  - k_r*sc[5];
    J[64] -= dqdci;               // dwdot[O]/d[HO2]
    J[65] += dqdci;               // dwdot[OH]/d[HO2]
    J[66] += dqdci;               // dwdot[HO2]/d[HO2]
    J[67] -= dqdci;               // dwdot[H2O2]/d[HO2]
    // d()/d[H2O2]
    dqdci =  + k_f*sc[4];
    J[74] -= dqdci;               // dwdot[O]/d[H2O2]
    J[75] += dqdci;               // dwdot[OH]/d[H2O2]
    J[76] += dqdci;               // dwdot[HO2]/d[H2O2]
    J[77] -= dqdci;               // dwdot[H2O2]/d[H2O2]
    // d()/dT
    J[94] -= dqdT;                // dwdot[O]/dT
    J[95] += dqdT;                // dwdot[OH]/dT
    J[96] += dqdT;                // dwdot[HO2]/dT
    J[97] -= dqdT;                // dwdot[H2O2]/dT

    /*reaction 20: H2O2 + OH <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    // forward
    phi_f = sc[5]*sc[7];
    k_f = 1e-06 * 1000000000000
                * exp(0 * tc[0] - 0.50321666580472 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580472 *  0  * invT2;
    // reverse
    phi_r = sc[2]*sc[6];
    Kc = exp(-g_RT[2] + g_RT[5] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[7]) + (h_RT[2] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    // rate of progress
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    // update wdot
    wdot[2] += q; // H2O
    wdot[5] -= q; // OH
    wdot[6] += q; // HO2
    wdot[7] -= q; // H2O2
    // d()/d[H2O]
    dqdci =  - k_r*sc[6];
    J[22] += dqdci;               // dwdot[H2O]/d[H2O]
    J[25] -= dqdci;               // dwdot[OH]/d[H2O]
    J[26] += dqdci;               // dwdot[HO2]/d[H2O]
    J[27] -= dqdci;               // dwdot[H2O2]/d[H2O]
    // d()/d[OH]
    dqdci =  + k_f*sc[7];
    J[52] += dqdci;               // dwdot[H2O]/d[OH]
    J[55] -= dqdci;               // dwdot[OH]/d[OH]
    J[56] += dqdci;               // dwdot[HO2]/d[OH]
    J[57] -= dqdci;               // dwdot[H2O2]/d[OH]
    // d()/d[HO2]
    dqdci =  - k_r*sc[2];
    J[62] += dqdci;               // dwdot[H2O]/d[HO2]
    J[65] -= dqdci;               // dwdot[OH]/d[HO2]
    J[66] += dqdci;               // dwdot[HO2]/d[HO2]
    J[67] -= dqdci;               // dwdot[H2O2]/d[HO2]
    // d()/d[H2O2]
    dqdci =  + k_f*sc[5];
    J[72] += dqdci;               // dwdot[H2O]/d[H2O2]
    J[75] -= dqdci;               // dwdot[OH]/d[H2O2]
    J[76] += dqdci;               // dwdot[HO2]/d[H2O2]
    J[77] -= dqdci;               // dwdot[H2O2]/d[H2O2]
    // d()/dT
    J[92] += dqdT;                // dwdot[H2O]/dT
    J[95] -= dqdT;                // dwdot[OH]/dT
    J[96] += dqdT;                // dwdot[HO2]/dT
    J[97] -= dqdT;                // dwdot[H2O2]/dT

    /*reaction 21: H2O2 + OH <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    // forward
    phi_f = sc[5]*sc[7];
    k_f = 1e-06 * 580000000000000
                * exp(0 * tc[0] - 0.50321666580472 * (9557) * invT);
    dlnkfdT = 0 * invT + 0.50321666580472 *  9557  * invT2;
    // reverse
    phi_r = sc[2]*sc[6];
    Kc = exp(-g_RT[2] + g_RT[5] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[5] + h_RT[7]) + (h_RT[2] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    // rate of progress
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    // update wdot
    wdot[2] += q; // H2O
    wdot[5] -= q; // OH
    wdot[6] += q; // HO2
    wdot[7] -= q; // H2O2
    // d()/d[H2O]
    dqdci =  - k_r*sc[6];
    J[22] += dqdci;               // dwdot[H2O]/d[H2O]
    J[25] -= dqdci;               // dwdot[OH]/d[H2O]
    J[26] += dqdci;               // dwdot[HO2]/d[H2O]
    J[27] -= dqdci;               // dwdot[H2O2]/d[H2O]
    // d()/d[OH]
    dqdci =  + k_f*sc[7];
    J[52] += dqdci;               // dwdot[H2O]/d[OH]
    J[55] -= dqdci;               // dwdot[OH]/d[OH]
    J[56] += dqdci;               // dwdot[HO2]/d[OH]
    J[57] -= dqdci;               // dwdot[H2O2]/d[OH]
    // d()/d[HO2]
    dqdci =  - k_r*sc[2];
    J[62] += dqdci;               // dwdot[H2O]/d[HO2]
    J[65] -= dqdci;               // dwdot[OH]/d[HO2]
    J[66] += dqdci;               // dwdot[HO2]/d[HO2]
    J[67] -= dqdci;               // dwdot[H2O2]/d[HO2]
    // d()/d[H2O2]
    dqdci =  + k_f*sc[5];
    J[72] += dqdci;               // dwdot[H2O]/d[H2O2]
    J[75] -= dqdci;               // dwdot[OH]/d[H2O2]
    J[76] += dqdci;               // dwdot[HO2]/d[H2O2]
    J[77] -= dqdci;               // dwdot[H2O2]/d[H2O2]
    // d()/dT
    J[92] += dqdT;                // dwdot[H2O]/dT
    J[95] -= dqdT;                // dwdot[OH]/dT
    J[96] += dqdT;                // dwdot[HO2]/dT
    J[97] -= dqdT;                // dwdot[H2O2]/dT

    amrex::Real c_R[9], dcRdT[9], e_RT[9];
    amrex::Real * eh_RT;
    if (consP) {
        cp_R(c_R, tc);
        dcvpRdT(dcRdT, tc);
        eh_RT = &h_RT[0];
    }
    else {
        cv_R(c_R, tc);
        dcvpRdT(dcRdT, tc);
        speciesInternalEnergy(e_RT, tc);
        eh_RT = &e_RT[0];
    }

    amrex::Real cmix = 0.0, ehmix = 0.0, dcmixdT=0.0, dehmixdT=0.0;
    for (int k = 0; k < 9; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[90+k];
    }

    amrex::Real cmixinv = 1.0/cmix;
    amrex::Real tmp1 = ehmix*cmixinv;
    amrex::Real tmp3 = cmixinv*T;
    amrex::Real tmp2 = tmp1*tmp3;
    amrex::Real dehmixdc;
    // dTdot/d[X]
    for (int k = 0; k < 9; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 9; ++m) {
            dehmixdc += eh_RT[m]*J[k*10+m];
        }
        J[k*10+9] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    // dTdot/dT
    J[99] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;

return;
}


/*compute the reaction Jacobian */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void DWDOT(amrex::Real *  J, amrex::Real *  sc, amrex::Real *  Tp, const int * consP)
{
    amrex::Real c[9];

    for (int k=0; k<9; k++) {
        c[k] = 1.e6 * sc[k];
    }

    aJacobian(J, c, *Tp, *consP);

    // dwdot[k]/dT
    // dTdot/d[X]
    for (int k=0; k<9; k++) {
        J[90+k] *= 1.e-6;
        J[k*10+9] *= 1.e6;
    }

    return;
}


/* Transport function declarations  */


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetLENIMC(int* LENIMC ) {
    *LENIMC = 38;}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetLENRMC(int* LENRMC ) {
    *LENRMC = 1854;}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetNO(int* NO ) {
    *NO = 4;}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetKK(int* KK ) {
    *KK = 9;}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetNLITE(int* NLITE ) {
    *NLITE = 2;}


/*Patm in ergs/cm3 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetPATM(amrex::Real* PATM) {
    *PATM =   0.1013250000000000E+07;}


/*the molecular weights in g/mol */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetWT(amrex::Real* WT ) {
    WT[0] = 2.01594000E+00;
    WT[1] = 3.19988000E+01;
    WT[2] = 1.80153400E+01;
    WT[3] = 1.00797000E+00;
    WT[4] = 1.59994000E+01;
    WT[5] = 1.70073700E+01;
    WT[6] = 3.30067700E+01;
    WT[7] = 3.40147400E+01;
    WT[8] = 2.80134000E+01;
}


/*the lennard-jones potential well depth eps/kb in K */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetEPS(amrex::Real* EPS ) {
    EPS[0] = 3.80000000E+01;
    EPS[1] = 1.07400000E+02;
    EPS[2] = 5.72400000E+02;
    EPS[3] = 1.45000000E+02;
    EPS[4] = 8.00000000E+01;
    EPS[5] = 8.00000000E+01;
    EPS[6] = 1.07400000E+02;
    EPS[7] = 1.07400000E+02;
    EPS[8] = 9.75300000E+01;
}


/*the lennard-jones collision diameter in Angstroms */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetSIG(amrex::Real* SIG ) {
    SIG[0] = 2.92000000E+00;
    SIG[1] = 3.45800000E+00;
    SIG[2] = 2.60500000E+00;
    SIG[3] = 2.05000000E+00;
    SIG[4] = 2.75000000E+00;
    SIG[5] = 2.75000000E+00;
    SIG[6] = 3.45800000E+00;
    SIG[7] = 3.45800000E+00;
    SIG[8] = 3.62100000E+00;
}


/*the dipole moment in Debye */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetDIP(amrex::Real* DIP ) {
    DIP[0] = 0.00000000E+00;
    DIP[1] = 0.00000000E+00;
    DIP[2] = 1.84400000E+00;
    DIP[3] = 0.00000000E+00;
    DIP[4] = 0.00000000E+00;
    DIP[5] = 0.00000000E+00;
    DIP[6] = 0.00000000E+00;
    DIP[7] = 0.00000000E+00;
    DIP[8] = 0.00000000E+00;
}


/*the polarizability in cubic Angstroms */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetPOL(amrex::Real* POL ) {
    POL[0] = 7.90000000E-01;
    POL[1] = 1.60000000E+00;
    POL[2] = 0.00000000E+00;
    POL[3] = 0.00000000E+00;
    POL[4] = 0.00000000E+00;
    POL[5] = 0.00000000E+00;
    POL[6] = 0.00000000E+00;
    POL[7] = 0.00000000E+00;
    POL[8] = 1.76000000E+00;
}


/*the rotational relaxation collision number at 298 K */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetZROT(amrex::Real* ZROT ) {
    ZROT[0] = 2.80000000E+02;
    ZROT[1] = 3.80000000E+00;
    ZROT[2] = 4.00000000E+00;
    ZROT[3] = 0.00000000E+00;
    ZROT[4] = 0.00000000E+00;
    ZROT[5] = 0.00000000E+00;
    ZROT[6] = 1.00000000E+00;
    ZROT[7] = 3.80000000E+00;
    ZROT[8] = 4.00000000E+00;
}


/*0: monoatomic, 1: linear, 2: nonlinear */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetNLIN(int* NLIN) {
    NLIN[0] = 1;
    NLIN[1] = 1;
    NLIN[2] = 2;
    NLIN[3] = 0;
    NLIN[4] = 0;
    NLIN[5] = 1;
    NLIN[6] = 2;
    NLIN[7] = 2;
    NLIN[8] = 1;
}


/*Poly fits for the viscosities, dim NO*KK */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetCOFETA(amrex::Real* COFETA) {
    COFETA[0] = -1.37549435E+01;
    COFETA[1] = 9.65530587E-01;
    COFETA[2] = -4.45720114E-02;
    COFETA[3] = 2.05871810E-03;
    COFETA[4] = -1.68118868E+01;
    COFETA[5] = 2.52362554E+00;
    COFETA[6] = -2.49309128E-01;
    COFETA[7] = 1.10211025E-02;
    COFETA[8] = -1.17770937E+01;
    COFETA[9] = -8.26742721E-01;
    COFETA[10] = 3.39009079E-01;
    COFETA[11] = -2.00674327E-02;
    COFETA[12] = -1.98744496E+01;
    COFETA[13] = 3.41660514E+00;
    COFETA[14] = -3.63206306E-01;
    COFETA[15] = 1.58671021E-02;
    COFETA[16] = -1.48001581E+01;
    COFETA[17] = 1.79491990E+00;
    COFETA[18] = -1.54008440E-01;
    COFETA[19] = 6.86719439E-03;
    COFETA[20] = -1.47696103E+01;
    COFETA[21] = 1.79491990E+00;
    COFETA[22] = -1.54008440E-01;
    COFETA[23] = 6.86719439E-03;
    COFETA[24] = -1.67963797E+01;
    COFETA[25] = 2.52362554E+00;
    COFETA[26] = -2.49309128E-01;
    COFETA[27] = 1.10211025E-02;
    COFETA[28] = -1.67813391E+01;
    COFETA[29] = 2.52362554E+00;
    COFETA[30] = -2.49309128E-01;
    COFETA[31] = 1.10211025E-02;
    COFETA[32] = -1.62526779E+01;
    COFETA[33] = 2.24839597E+00;
    COFETA[34] = -2.13428438E-01;
    COFETA[35] = 9.46192413E-03;
}


/*Poly fits for the conductivities, dim NO*KK */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetCOFLAM(amrex::Real* COFLAM) {
    COFLAM[0] = 1.11035551E+01;
    COFLAM[1] = -1.31884089E+00;
    COFLAM[2] = 2.44042735E-01;
    COFLAM[3] = -8.99837633E-03;
    COFLAM[4] = -2.51299170E+00;
    COFLAM[5] = 3.15166792E+00;
    COFLAM[6] = -3.10009273E-01;
    COFLAM[7] = 1.34523096E-02;
    COFLAM[8] = 2.21730076E+01;
    COFLAM[9] = -8.46933644E+00;
    COFLAM[10] = 1.46153515E+00;
    COFLAM[11] = -7.29500936E-02;
    COFLAM[12] = -3.24539191E-01;
    COFLAM[13] = 3.41660514E+00;
    COFLAM[14] = -3.63206306E-01;
    COFLAM[15] = 1.58671021E-02;
    COFLAM[16] = 1.98513952E+00;
    COFLAM[17] = 1.79491990E+00;
    COFLAM[18] = -1.54008440E-01;
    COFLAM[19] = 6.86719439E-03;
    COFLAM[20] = 1.60618734E+01;
    COFLAM[21] = -4.10626869E+00;
    COFLAM[22] = 6.63571339E-01;
    COFLAM[23] = -2.97906324E-02;
    COFLAM[24] = 5.56023763E-01;
    COFLAM[25] = 1.59073590E+00;
    COFLAM[26] = -5.28053839E-02;
    COFLAM[27] = 4.07601571E-04;
    COFLAM[28] = 1.48801857E+00;
    COFLAM[29] = 1.06175905E+00;
    COFLAM[30] = 5.72200266E-02;
    COFLAM[31] = -6.38393531E-03;
    COFLAM[32] = 1.15507419E+01;
    COFLAM[33] = -2.91453917E+00;
    COFLAM[34] = 5.55045765E-01;
    COFLAM[35] = -2.75173485E-02;
}


/*Poly fits for the diffusion coefficients, dim NO*KK*KK */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetCOFD(amrex::Real* COFD) {
    COFD[0] = -1.02395222E+01;
    COFD[1] = 2.15403244E+00;
    COFD[2] = -6.97480266E-02;
    COFD[3] = 3.23666871E-03;
    COFD[4] = -1.15797750E+01;
    COFD[5] = 2.43235504E+00;
    COFD[6] = -1.02890179E-01;
    COFD[7] = 4.52903603E-03;
    COFD[8] = -1.68758926E+01;
    COFD[9] = 4.49460303E+00;
    COFD[10] = -3.64766132E-01;
    COFD[11] = 1.56457153E-02;
    COFD[12] = -1.11808682E+01;
    COFD[13] = 2.66936727E+00;
    COFD[14] = -1.34411514E-01;
    COFD[15] = 5.92957488E-03;
    COFD[16] = -1.06250182E+01;
    COFD[17] = 2.15849701E+00;
    COFD[18] = -6.53886401E-02;
    COFD[19] = 2.81453370E-03;
    COFD[20] = -1.06283453E+01;
    COFD[21] = 2.15849701E+00;
    COFD[22] = -6.53886401E-02;
    COFD[23] = 2.81453370E-03;
    COFD[24] = -1.15806808E+01;
    COFD[25] = 2.43235504E+00;
    COFD[26] = -1.02890179E-01;
    COFD[27] = 4.52903603E-03;
    COFD[28] = -1.15815344E+01;
    COFD[29] = 2.43235504E+00;
    COFD[30] = -1.02890179E-01;
    COFD[31] = 4.52903603E-03;
    COFD[32] = -1.13253458E+01;
    COFD[33] = 2.31195095E+00;
    COFD[34] = -8.63988037E-02;
    COFD[35] = 3.77573452E-03;
    COFD[36] = -1.15797750E+01;
    COFD[37] = 2.43235504E+00;
    COFD[38] = -1.02890179E-01;
    COFD[39] = 4.52903603E-03;
    COFD[40] = -1.53110708E+01;
    COFD[41] = 3.37317428E+00;
    COFD[42] = -2.24900439E-01;
    COFD[43] = 9.81228151E-03;
    COFD[44] = -2.10640014E+01;
    COFD[45] = 5.50980695E+00;
    COFD[46] = -4.78335488E-01;
    COFD[47] = 1.98515434E-02;
    COFD[48] = -1.43712864E+01;
    COFD[49] = 3.70920439E+00;
    COFD[50] = -2.67274113E-01;
    COFD[51] = 1.15967481E-02;
    COFD[52] = -1.40864894E+01;
    COFD[53] = 3.07458927E+00;
    COFD[54] = -1.86899591E-01;
    COFD[55] = 8.19829781E-03;
    COFD[56] = -1.41066459E+01;
    COFD[57] = 3.07458927E+00;
    COFD[58] = -1.86899591E-01;
    COFD[59] = 8.19829781E-03;
    COFD[60] = -1.53187643E+01;
    COFD[61] = 3.37317428E+00;
    COFD[62] = -2.24900439E-01;
    COFD[63] = 9.81228151E-03;
    COFD[64] = -1.53261114E+01;
    COFD[65] = 3.37317428E+00;
    COFD[66] = -2.24900439E-01;
    COFD[67] = 9.81228151E-03;
    COFD[68] = -1.50096240E+01;
    COFD[69] = 3.25515933E+00;
    COFD[70] = -2.09710110E-01;
    COFD[71] = 9.15941830E-03;
    COFD[72] = -1.68758926E+01;
    COFD[73] = 4.49460303E+00;
    COFD[74] = -3.64766132E-01;
    COFD[75] = 1.56457153E-02;
    COFD[76] = -2.10640014E+01;
    COFD[77] = 5.50980695E+00;
    COFD[78] = -4.78335488E-01;
    COFD[79] = 1.98515434E-02;
    COFD[80] = -1.31492641E+01;
    COFD[81] = 1.48004311E+00;
    COFD[82] = 1.60499553E-01;
    COFD[83] = -1.19765679E-02;
    COFD[84] = -1.93611051E+01;
    COFD[85] = 5.51579726E+00;
    COFD[86] = -4.76061961E-01;
    COFD[87] = 1.96329391E-02;
    COFD[88] = -1.91096797E+01;
    COFD[89] = 5.02608697E+00;
    COFD[90] = -4.26959993E-01;
    COFD[91] = 1.80709910E-02;
    COFD[92] = -1.91256261E+01;
    COFD[93] = 5.02608697E+00;
    COFD[94] = -4.26959993E-01;
    COFD[95] = 1.80709910E-02;
    COFD[96] = -2.04177482E+01;
    COFD[97] = 5.31457079E+00;
    COFD[98] = -4.58216496E-01;
    COFD[99] = 1.91825910E-02;
    COFD[100] = -2.04230073E+01;
    COFD[101] = 5.31457079E+00;
    COFD[102] = -4.58216496E-01;
    COFD[103] = 1.91825910E-02;
    COFD[104] = -2.08123325E+01;
    COFD[105] = 5.42470154E+00;
    COFD[106] = -4.69700416E-01;
    COFD[107] = 1.95706904E-02;
    COFD[108] = -1.11808682E+01;
    COFD[109] = 2.66936727E+00;
    COFD[110] = -1.34411514E-01;
    COFD[111] = 5.92957488E-03;
    COFD[112] = -1.43712864E+01;
    COFD[113] = 3.70920439E+00;
    COFD[114] = -2.67274113E-01;
    COFD[115] = 1.15967481E-02;
    COFD[116] = -1.93611051E+01;
    COFD[117] = 5.51579726E+00;
    COFD[118] = -4.76061961E-01;
    COFD[119] = 1.96329391E-02;
    COFD[120] = -1.43693056E+01;
    COFD[121] = 4.03992999E+00;
    COFD[122] = -3.08044800E-01;
    COFD[123] = 1.32757775E-02;
    COFD[124] = -1.31860117E+01;
    COFD[125] = 3.38003453E+00;
    COFD[126] = -2.25783856E-01;
    COFD[127] = 9.85028660E-03;
    COFD[128] = -1.31877711E+01;
    COFD[129] = 3.38003453E+00;
    COFD[130] = -2.25783856E-01;
    COFD[131] = 9.85028660E-03;
    COFD[132] = -1.43717529E+01;
    COFD[133] = 3.70920439E+00;
    COFD[134] = -2.67274113E-01;
    COFD[135] = 1.15967481E-02;
    COFD[136] = -1.43721922E+01;
    COFD[137] = 3.70920439E+00;
    COFD[138] = -2.67274113E-01;
    COFD[139] = 1.15967481E-02;
    COFD[140] = -1.40298830E+01;
    COFD[141] = 3.55837688E+00;
    COFD[142] = -2.47785790E-01;
    COFD[143] = 1.07555332E-02;
    COFD[144] = -1.06250182E+01;
    COFD[145] = 2.15849701E+00;
    COFD[146] = -6.53886401E-02;
    COFD[147] = 2.81453370E-03;
    COFD[148] = -1.40864894E+01;
    COFD[149] = 3.07458927E+00;
    COFD[150] = -1.86899591E-01;
    COFD[151] = 8.19829781E-03;
    COFD[152] = -1.91096797E+01;
    COFD[153] = 5.02608697E+00;
    COFD[154] = -4.26959993E-01;
    COFD[155] = 1.80709910E-02;
    COFD[156] = -1.31860117E+01;
    COFD[157] = 3.38003453E+00;
    COFD[158] = -2.25783856E-01;
    COFD[159] = 9.85028660E-03;
    COFD[160] = -1.29877365E+01;
    COFD[161] = 2.80841511E+00;
    COFD[162] = -1.52629888E-01;
    COFD[163] = 6.72604927E-03;
    COFD[164] = -1.30027772E+01;
    COFD[165] = 2.80841511E+00;
    COFD[166] = -1.52629888E-01;
    COFD[167] = 6.72604927E-03;
    COFD[168] = -1.40916052E+01;
    COFD[169] = 3.07458927E+00;
    COFD[170] = -1.86899591E-01;
    COFD[171] = 8.19829781E-03;
    COFD[172] = -1.40964661E+01;
    COFD[173] = 3.07458927E+00;
    COFD[174] = -1.86899591E-01;
    COFD[175] = 8.19829781E-03;
    COFD[176] = -1.38756407E+01;
    COFD[177] = 2.98558426E+00;
    COFD[178] = -1.75507216E-01;
    COFD[179] = 7.71173691E-03;
    COFD[180] = -1.06283453E+01;
    COFD[181] = 2.15849701E+00;
    COFD[182] = -6.53886401E-02;
    COFD[183] = 2.81453370E-03;
    COFD[184] = -1.41066459E+01;
    COFD[185] = 3.07458927E+00;
    COFD[186] = -1.86899591E-01;
    COFD[187] = 8.19829781E-03;
    COFD[188] = -1.91256261E+01;
    COFD[189] = 5.02608697E+00;
    COFD[190] = -4.26959993E-01;
    COFD[191] = 1.80709910E-02;
    COFD[192] = -1.31877711E+01;
    COFD[193] = 3.38003453E+00;
    COFD[194] = -2.25783856E-01;
    COFD[195] = 9.85028660E-03;
    COFD[196] = -1.30027772E+01;
    COFD[197] = 2.80841511E+00;
    COFD[198] = -1.52629888E-01;
    COFD[199] = 6.72604927E-03;
    COFD[200] = -1.30182843E+01;
    COFD[201] = 2.80841511E+00;
    COFD[202] = -1.52629888E-01;
    COFD[203] = 6.72604927E-03;
    COFD[204] = -1.41119732E+01;
    COFD[205] = 3.07458927E+00;
    COFD[206] = -1.86899591E-01;
    COFD[207] = 8.19829781E-03;
    COFD[208] = -1.41170372E+01;
    COFD[209] = 3.07458927E+00;
    COFD[210] = -1.86899591E-01;
    COFD[211] = 8.19829781E-03;
    COFD[212] = -1.38948667E+01;
    COFD[213] = 2.98558426E+00;
    COFD[214] = -1.75507216E-01;
    COFD[215] = 7.71173691E-03;
    COFD[216] = -1.15806808E+01;
    COFD[217] = 2.43235504E+00;
    COFD[218] = -1.02890179E-01;
    COFD[219] = 4.52903603E-03;
    COFD[220] = -1.53187643E+01;
    COFD[221] = 3.37317428E+00;
    COFD[222] = -2.24900439E-01;
    COFD[223] = 9.81228151E-03;
    COFD[224] = -2.04177482E+01;
    COFD[225] = 5.31457079E+00;
    COFD[226] = -4.58216496E-01;
    COFD[227] = 1.91825910E-02;
    COFD[228] = -1.43717529E+01;
    COFD[229] = 3.70920439E+00;
    COFD[230] = -2.67274113E-01;
    COFD[231] = 1.15967481E-02;
    COFD[232] = -1.40916052E+01;
    COFD[233] = 3.07458927E+00;
    COFD[234] = -1.86899591E-01;
    COFD[235] = 8.19829781E-03;
    COFD[236] = -1.41119732E+01;
    COFD[237] = 3.07458927E+00;
    COFD[238] = -1.86899591E-01;
    COFD[239] = 8.19829781E-03;
    COFD[240] = -1.53265780E+01;
    COFD[241] = 3.37317428E+00;
    COFD[242] = -2.24900439E-01;
    COFD[243] = 9.81228151E-03;
    COFD[244] = -1.53340417E+01;
    COFD[245] = 3.37317428E+00;
    COFD[246] = -2.24900439E-01;
    COFD[247] = 9.81228151E-03;
    COFD[248] = -1.50168028E+01;
    COFD[249] = 3.25515933E+00;
    COFD[250] = -2.09710110E-01;
    COFD[251] = 9.15941830E-03;
    COFD[252] = -1.15815344E+01;
    COFD[253] = 2.43235504E+00;
    COFD[254] = -1.02890179E-01;
    COFD[255] = 4.52903603E-03;
    COFD[256] = -1.53261114E+01;
    COFD[257] = 3.37317428E+00;
    COFD[258] = -2.24900439E-01;
    COFD[259] = 9.81228151E-03;
    COFD[260] = -2.04230073E+01;
    COFD[261] = 5.31457079E+00;
    COFD[262] = -4.58216496E-01;
    COFD[263] = 1.91825910E-02;
    COFD[264] = -1.43721922E+01;
    COFD[265] = 3.70920439E+00;
    COFD[266] = -2.67274113E-01;
    COFD[267] = 1.15967481E-02;
    COFD[268] = -1.40964661E+01;
    COFD[269] = 3.07458927E+00;
    COFD[270] = -1.86899591E-01;
    COFD[271] = 8.19829781E-03;
    COFD[272] = -1.41170372E+01;
    COFD[273] = 3.07458927E+00;
    COFD[274] = -1.86899591E-01;
    COFD[275] = 8.19829781E-03;
    COFD[276] = -1.53340417E+01;
    COFD[277] = 3.37317428E+00;
    COFD[278] = -2.24900439E-01;
    COFD[279] = 9.81228151E-03;
    COFD[280] = -1.53416186E+01;
    COFD[281] = 3.37317428E+00;
    COFD[282] = -2.24900439E-01;
    COFD[283] = 9.81228151E-03;
    COFD[284] = -1.50236516E+01;
    COFD[285] = 3.25515933E+00;
    COFD[286] = -2.09710110E-01;
    COFD[287] = 9.15941830E-03;
    COFD[288] = -1.13253458E+01;
    COFD[289] = 2.31195095E+00;
    COFD[290] = -8.63988037E-02;
    COFD[291] = 3.77573452E-03;
    COFD[292] = -1.50096240E+01;
    COFD[293] = 3.25515933E+00;
    COFD[294] = -2.09710110E-01;
    COFD[295] = 9.15941830E-03;
    COFD[296] = -2.08123325E+01;
    COFD[297] = 5.42470154E+00;
    COFD[298] = -4.69700416E-01;
    COFD[299] = 1.95706904E-02;
    COFD[300] = -1.40298830E+01;
    COFD[301] = 3.55837688E+00;
    COFD[302] = -2.47785790E-01;
    COFD[303] = 1.07555332E-02;
    COFD[304] = -1.38756407E+01;
    COFD[305] = 2.98558426E+00;
    COFD[306] = -1.75507216E-01;
    COFD[307] = 7.71173691E-03;
    COFD[308] = -1.38948667E+01;
    COFD[309] = 2.98558426E+00;
    COFD[310] = -1.75507216E-01;
    COFD[311] = 7.71173691E-03;
    COFD[312] = -1.50168028E+01;
    COFD[313] = 3.25515933E+00;
    COFD[314] = -2.09710110E-01;
    COFD[315] = 9.15941830E-03;
    COFD[316] = -1.50236516E+01;
    COFD[317] = 3.25515933E+00;
    COFD[318] = -2.09710110E-01;
    COFD[319] = 9.15941830E-03;
    COFD[320] = -1.47639290E+01;
    COFD[321] = 3.15955654E+00;
    COFD[322] = -1.97590757E-01;
    COFD[323] = 8.64692156E-03;
}


/*List of specs with small weight, dim NLITE */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetKTDIF(int* KTDIF) {
    KTDIF[0] = 0;
    KTDIF[1] = 3;
}


/*Poly fits for thermal diff ratios, dim NO*NLITE*KK */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetCOFTD(amrex::Real* COFTD) {
    COFTD[0] = 0.00000000E+00;
    COFTD[1] = 0.00000000E+00;
    COFTD[2] = 0.00000000E+00;
    COFTD[3] = 0.00000000E+00;
    COFTD[4] = 4.42739084E-01;
    COFTD[5] = 7.11770818E-05;
    COFTD[6] = -3.84768062E-08;
    COFTD[7] = 6.86323437E-12;
    COFTD[8] = 6.02028221E-02;
    COFTD[9] = 5.61561867E-04;
    COFTD[10] = -2.55372862E-07;
    COFTD[11] = 3.63389913E-11;
    COFTD[12] = -1.52534742E-01;
    COFTD[13] = -5.46404022E-05;
    COFTD[14] = 2.93412470E-08;
    COFTD[15] = -4.87091914E-12;
    COFTD[16] = 4.15583337E-01;
    COFTD[17] = 1.09738399E-05;
    COFTD[18] = -3.96021963E-09;
    COFTD[19] = 1.14414443E-12;
    COFTD[20] = 4.21932443E-01;
    COFTD[21] = 1.11414935E-05;
    COFTD[22] = -4.02072219E-09;
    COFTD[23] = 1.16162418E-12;
    COFTD[24] = 4.44452569E-01;
    COFTD[25] = 7.14525507E-05;
    COFTD[26] = -3.86257187E-08;
    COFTD[27] = 6.88979640E-12;
    COFTD[28] = 4.46070183E-01;
    COFTD[29] = 7.17126069E-05;
    COFTD[30] = -3.87662996E-08;
    COFTD[31] = 6.91487226E-12;
    COFTD[32] = 4.45261966E-01;
    COFTD[33] = 4.94697174E-05;
    COFTD[34] = -2.63023442E-08;
    COFTD[35] = 4.90306217E-12;
    COFTD[36] = 1.52534742E-01;
    COFTD[37] = 5.46404022E-05;
    COFTD[38] = -2.93412470E-08;
    COFTD[39] = 4.87091914E-12;
    COFTD[40] = 2.20482843E-01;
    COFTD[41] = 4.80164288E-04;
    COFTD[42] = -2.32927944E-07;
    COFTD[43] = 3.46470436E-11;
    COFTD[44] = -1.41883744E-01;
    COFTD[45] = 7.66558810E-04;
    COFTD[46] = -3.06550003E-07;
    COFTD[47] = 4.02959502E-11;
    COFTD[48] = 0.00000000E+00;
    COFTD[49] = 0.00000000E+00;
    COFTD[50] = 0.00000000E+00;
    COFTD[51] = 0.00000000E+00;
    COFTD[52] = 2.70010150E-01;
    COFTD[53] = 3.61555093E-04;
    COFTD[54] = -1.80744752E-07;
    COFTD[55] = 2.75321248E-11;
    COFTD[56] = 2.72041664E-01;
    COFTD[57] = 3.64275376E-04;
    COFTD[58] = -1.82104647E-07;
    COFTD[59] = 2.77392722E-11;
    COFTD[60] = 2.20907853E-01;
    COFTD[61] = 4.81089870E-04;
    COFTD[62] = -2.33376944E-07;
    COFTD[63] = 3.47138305E-11;
    COFTD[64] = 2.21308399E-01;
    COFTD[65] = 4.81962174E-04;
    COFTD[66] = -2.33800100E-07;
    COFTD[67] = 3.47767730E-11;
    COFTD[68] = 2.40744421E-01;
    COFTD[69] = 4.45343451E-04;
    COFTD[70] = -2.18173874E-07;
    COFTD[71] = 3.26958506E-11;
}


/*compute the critical parameters for each species */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void GET_CRITPARAMS(amrex::Real *  Tci, amrex::Real *  ai, amrex::Real *  bi, amrex::Real *  acentric_i)
{

    amrex::Real   EPS[9];
    amrex::Real   SIG[9];
    amrex::Real    wt[9];
    amrex::Real avogadro = 6.02214199e23;
    amrex::Real boltzmann = 1.3806503e-16; //we work in CGS
    amrex::Real Rcst = 83.144598; //in bar [CGS] !

    egtransetEPS(EPS);
    egtransetSIG(SIG);
    get_mw(wt);

    /*species 0: H2 */
    /*Imported from NIST */
    Tci[0] = 33.145000 ; 
    ai[0] = 1e6 * 0.42748 * Rcst * Rcst * Tci[0] * Tci[0] / (2.015880 * 2.015880 * 12.964000); 
    bi[0] = 0.08664 * Rcst * Tci[0] / (2.015880 * 12.964000); 
    acentric_i[0] = -0.219000 ;

    /*species 1: O2 */
    /*Imported from NIST */
    Tci[1] = 154.581000 ; 
    ai[1] = 1e6 * 0.42748 * Rcst * Rcst * Tci[1] * Tci[1] / (31.998800 * 31.998800 * 50.430466); 
    bi[1] = 0.08664 * Rcst * Tci[1] / (31.998800 * 50.430466); 
    acentric_i[1] = 0.022200 ;

    /*species 2: H2O */
    /*Imported from NIST */
    Tci[2] = 647.096000 ; 
    ai[2] = 1e6 * 0.42748 * Rcst * Rcst * Tci[2] * Tci[2] / (18.015340 * 18.015340 * 220.640000); 
    bi[2] = 0.08664 * Rcst * Tci[2] / (18.015340 * 220.640000); 
    acentric_i[2] = 0.344300 ;

    /*species 3: H */
    Tci[3] = 1.316 * EPS[3] ; 
    ai[3] = (5.55 * avogadro * avogadro * EPS[3]*boltzmann * pow(1e-8*SIG[3],3.0) ) / (wt[3] * wt[3]); 
    bi[3] = 0.855 * avogadro * pow(1e-8*SIG[3],3.0) / (wt[3]); 
    acentric_i[3] = 0.0 ;

    /*species 4: O */
    Tci[4] = 1.316 * EPS[4] ; 
    ai[4] = (5.55 * avogadro * avogadro * EPS[4]*boltzmann * pow(1e-8*SIG[4],3.0) ) / (wt[4] * wt[4]); 
    bi[4] = 0.855 * avogadro * pow(1e-8*SIG[4],3.0) / (wt[4]); 
    acentric_i[4] = 0.0 ;

    /*species 5: OH */
    Tci[5] = 1.316 * EPS[5] ; 
    ai[5] = (5.55 * avogadro * avogadro * EPS[5]*boltzmann * pow(1e-8*SIG[5],3.0) ) / (wt[5] * wt[5]); 
    bi[5] = 0.855 * avogadro * pow(1e-8*SIG[5],3.0) / (wt[5]); 
    acentric_i[5] = 0.0 ;

    /*species 6: HO2 */
    Tci[6] = 1.316 * EPS[6] ; 
    ai[6] = (5.55 * avogadro * avogadro * EPS[6]*boltzmann * pow(1e-8*SIG[6],3.0) ) / (wt[6] * wt[6]); 
    bi[6] = 0.855 * avogadro * pow(1e-8*SIG[6],3.0) / (wt[6]); 
    acentric_i[6] = 0.0 ;

    /*species 7: H2O2 */
    Tci[7] = 1.316 * EPS[7] ; 
    ai[7] = (5.55 * avogadro * avogadro * EPS[7]*boltzmann * pow(1e-8*SIG[7],3.0) ) / (wt[7] * wt[7]); 
    bi[7] = 0.855 * avogadro * pow(1e-8*SIG[7],3.0) / (wt[7]); 
    acentric_i[7] = 0.0 ;

    /*species 8: N2 */
    /*Imported from NIST */
    Tci[8] = 126.192000 ; 
    ai[8] = 1e6 * 0.42748 * Rcst * Rcst * Tci[8] * Tci[8] / (28.013400 * 28.013400 * 33.958000); 
    bi[8] = 0.08664 * Rcst * Tci[8] / (28.013400 * 33.958000); 
    acentric_i[8] = 0.037200 ;

    return;
}

/* gauss-jordan solver external routine */
/* Replace this routine with the one generated by the Gauss Jordan solver of DW */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void sgjsolve(amrex::Real* /*A*/, amrex::Real* /*x*/, amrex::Real* /*b*/) {
    amrex::Abort("sgjsolve not implemented, choose a different solver ");
}

/* Replace this routine with the one generated by the Gauss Jordan solver of DW */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void sgjsolve_simplified(amrex::Real* /*A*/, amrex::Real* /*x*/, amrex::Real* /*b*/) {
    amrex::Abort("sgjsolve_simplified not implemented, choose a different solver ");
}

#endif

/* End of file  */
