#ifndef CHEMISTRY_FILE_H
#define CHEMISTRY_FILE_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <vector>
#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>
#include <mechanism.h>



/* ALWAYS on CPU stuff -- can have different def depending on if we are CPU or GPU based. Defined in mechanism.cpp  */
/* INIT and FINALIZE stuff  */
void CKINIT();
void CKFINALIZE();
void atomicWeight(amrex::Real *  awt);
/* MISC  */
void CKNCF(int * ncf);
void CKSYME_STR(amrex::Vector<std::string>& ename);
void CKSYMS_STR(amrex::Vector<std::string>& kname);
void CKAWT(amrex::Real *  awt);
/* SPARSE INFORMATION  */
void SPARSITY_INFO(int * nJdata, int * consP, int NCELLS);
void SPARSITY_INFO_SYST(int * nJdata, int * consP, int NCELLS);
void SPARSITY_INFO_SYST_SIMPLIFIED(int * nJdata, int * consP);
void SPARSITY_PREPROC_CSC(int * rowVals, int * colPtrs, int * consP, int NCELLS);
void SPARSITY_PREPROC_CSR(int * colVals, int * rowPtrs, int * consP, int NCELLS, int base);
void SPARSITY_PREPROC_SYST_CSR(int * colVals, int * rowPtrs, int * consP, int NCELLS, int base);
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(int * rowVals, int * colPtrs, int * indx, int * consP);
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(int * colVals, int * rowPtr, int * consP, int base);

/* GPU stuff that default to CPU if not compiled with GPU  */


/*A few mechanism parameters */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKINDX(int * mm, int * kk, int * ii, int * nfit)
{
    *mm = 3;
    *kk = 9;
    *ii = 27;
    *nfit = -1; /*Why do you need this anyway ?  */
}


/* inverse molecular weights  */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void get_imw(amrex::Real *imw_new){
    imw_new[0] = 1.0/2.015940;/*H2 */
    imw_new[1] = 1.0/1.007970;/*H */
    imw_new[2] = 1.0/15.999400;/*O */
    imw_new[3] = 1.0/31.998800;/*O2 */
    imw_new[4] = 1.0/17.007370;/*OH */
    imw_new[5] = 1.0/18.015340;/*H2O */
    imw_new[6] = 1.0/33.006770;/*HO2 */
    imw_new[7] = 1.0/34.014740;/*H2O2 */
    imw_new[8] = 1.0/28.013400;/*N2 */
}

/* molecular weights  */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void get_mw(amrex::Real *mw_new){
    mw_new[0] = 2.015940;/*H2 */
    mw_new[1] = 1.007970;/*H */
    mw_new[2] = 15.999400;/*O */
    mw_new[3] = 31.998800;/*O2 */
    mw_new[4] = 17.007370;/*OH */
    mw_new[5] = 18.015340;/*H2O */
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
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void cv_R(amrex::Real * species, amrex::Real *  tc)
{

    /*temperature */
    amrex::Real T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            +1.34433112e+00
            +7.98052075e-03 * tc[1]
            -1.94781510e-05 * tc[2]
            +2.01572094e-08 * tc[3]
            -7.37611761e-12 * tc[4];
        /*species 1: H */
        species[1] =
            +1.50000000e+00
            +7.05332819e-13 * tc[1]
            -1.99591964e-15 * tc[2]
            +2.30081632e-18 * tc[3]
            -9.27732332e-22 * tc[4];
        /*species 2: O */
        species[2] =
            +2.16826710e+00
            -3.27931884e-03 * tc[1]
            +6.64306396e-06 * tc[2]
            -6.12806624e-09 * tc[3]
            +2.11265971e-12 * tc[4];
        /*species 3: O2 */
        species[3] =
            +2.78245636e+00
            -2.99673416e-03 * tc[1]
            +9.84730201e-06 * tc[2]
            -9.68129509e-09 * tc[3]
            +3.24372837e-12 * tc[4];
        /*species 4: OH */
        species[4] =
            +2.99201543e+00
            -2.40131752e-03 * tc[1]
            +4.61793841e-06 * tc[2]
            -3.88113333e-09 * tc[3]
            +1.36411470e-12 * tc[4];
        /*species 5: H2O */
        species[5] =
            +3.19864056e+00
            -2.03643410e-03 * tc[1]
            +6.52040211e-06 * tc[2]
            -5.48797062e-09 * tc[3]
            +1.77197817e-12 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +3.30179801e+00
            -4.74912051e-03 * tc[1]
            +2.11582891e-05 * tc[2]
            -2.42763894e-08 * tc[3]
            +9.29225124e-12 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +3.27611269e+00
            -5.42822417e-04 * tc[1]
            +1.67335701e-05 * tc[2]
            -2.15770813e-08 * tc[3]
            +8.62454363e-12 * tc[4];
        /*species 8: N2 */
        species[8] =
            +2.29867700e+00
            +1.40824040e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485400e-12 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            +2.33727920e+00
            -4.94024731e-05 * tc[1]
            +4.99456778e-07 * tc[2]
            -1.79566394e-10 * tc[3]
            +2.00255376e-14 * tc[4];
        /*species 1: H */
        species[1] =
            +1.50000001e+00
            -2.30842973e-11 * tc[1]
            +1.61561948e-14 * tc[2]
            -4.73515235e-18 * tc[3]
            +4.98197357e-22 * tc[4];
        /*species 2: O */
        species[2] =
            +1.56942078e+00
            -8.59741137e-05 * tc[1]
            +4.19484589e-08 * tc[2]
            -1.00177799e-11 * tc[3]
            +1.22833691e-15 * tc[4];
        /*species 3: O2 */
        species[3] =
            +2.28253784e+00
            +1.48308754e-03 * tc[1]
            -7.57966669e-07 * tc[2]
            +2.09470555e-10 * tc[3]
            -2.16717794e-14 * tc[4];
        /*species 4: OH */
        species[4] =
            +2.09288767e+00
            +5.48429716e-04 * tc[1]
            +1.26505228e-07 * tc[2]
            -8.79461556e-11 * tc[3]
            +1.17412376e-14 * tc[4];
        /*species 5: H2O */
        species[5] =
            +2.03399249e+00
            +2.17691804e-03 * tc[1]
            -1.64072518e-07 * tc[2]
            -9.70419870e-11 * tc[3]
            +1.68200992e-14 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +3.01721090e+00
            +2.23982013e-03 * tc[1]
            -6.33658150e-07 * tc[2]
            +1.14246370e-10 * tc[3]
            -1.07908535e-14 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +3.16500285e+00
            +4.90831694e-03 * tc[1]
            -1.90139225e-06 * tc[2]
            +3.71185986e-10 * tc[3]
            -2.87908305e-14 * tc[4];
        /*species 8: N2 */
        species[8] =
            +1.92664000e+00
            +1.48797680e-03 * tc[1]
            -5.68476000e-07 * tc[2]
            +1.00970380e-10 * tc[3]
            -6.75335100e-15 * tc[4];
    }
    return;
}


/*compute Cp/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void cp_R(amrex::Real * species, amrex::Real *  tc)
{

    /*temperature */
    amrex::Real T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            +2.34433112e+00
            +7.98052075e-03 * tc[1]
            -1.94781510e-05 * tc[2]
            +2.01572094e-08 * tc[3]
            -7.37611761e-12 * tc[4];
        /*species 1: H */
        species[1] =
            +2.50000000e+00
            +7.05332819e-13 * tc[1]
            -1.99591964e-15 * tc[2]
            +2.30081632e-18 * tc[3]
            -9.27732332e-22 * tc[4];
        /*species 2: O */
        species[2] =
            +3.16826710e+00
            -3.27931884e-03 * tc[1]
            +6.64306396e-06 * tc[2]
            -6.12806624e-09 * tc[3]
            +2.11265971e-12 * tc[4];
        /*species 3: O2 */
        species[3] =
            +3.78245636e+00
            -2.99673416e-03 * tc[1]
            +9.84730201e-06 * tc[2]
            -9.68129509e-09 * tc[3]
            +3.24372837e-12 * tc[4];
        /*species 4: OH */
        species[4] =
            +3.99201543e+00
            -2.40131752e-03 * tc[1]
            +4.61793841e-06 * tc[2]
            -3.88113333e-09 * tc[3]
            +1.36411470e-12 * tc[4];
        /*species 5: H2O */
        species[5] =
            +4.19864056e+00
            -2.03643410e-03 * tc[1]
            +6.52040211e-06 * tc[2]
            -5.48797062e-09 * tc[3]
            +1.77197817e-12 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +4.30179801e+00
            -4.74912051e-03 * tc[1]
            +2.11582891e-05 * tc[2]
            -2.42763894e-08 * tc[3]
            +9.29225124e-12 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +4.27611269e+00
            -5.42822417e-04 * tc[1]
            +1.67335701e-05 * tc[2]
            -2.15770813e-08 * tc[3]
            +8.62454363e-12 * tc[4];
        /*species 8: N2 */
        species[8] =
            +3.29867700e+00
            +1.40824040e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485400e-12 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            +3.33727920e+00
            -4.94024731e-05 * tc[1]
            +4.99456778e-07 * tc[2]
            -1.79566394e-10 * tc[3]
            +2.00255376e-14 * tc[4];
        /*species 1: H */
        species[1] =
            +2.50000001e+00
            -2.30842973e-11 * tc[1]
            +1.61561948e-14 * tc[2]
            -4.73515235e-18 * tc[3]
            +4.98197357e-22 * tc[4];
        /*species 2: O */
        species[2] =
            +2.56942078e+00
            -8.59741137e-05 * tc[1]
            +4.19484589e-08 * tc[2]
            -1.00177799e-11 * tc[3]
            +1.22833691e-15 * tc[4];
        /*species 3: O2 */
        species[3] =
            +3.28253784e+00
            +1.48308754e-03 * tc[1]
            -7.57966669e-07 * tc[2]
            +2.09470555e-10 * tc[3]
            -2.16717794e-14 * tc[4];
        /*species 4: OH */
        species[4] =
            +3.09288767e+00
            +5.48429716e-04 * tc[1]
            +1.26505228e-07 * tc[2]
            -8.79461556e-11 * tc[3]
            +1.17412376e-14 * tc[4];
        /*species 5: H2O */
        species[5] =
            +3.03399249e+00
            +2.17691804e-03 * tc[1]
            -1.64072518e-07 * tc[2]
            -9.70419870e-11 * tc[3]
            +1.68200992e-14 * tc[4];
        /*species 6: HO2 */
        species[6] =
            +4.01721090e+00
            +2.23982013e-03 * tc[1]
            -6.33658150e-07 * tc[2]
            +1.14246370e-10 * tc[3]
            -1.07908535e-14 * tc[4];
        /*species 7: H2O2 */
        species[7] =
            +4.16500285e+00
            +4.90831694e-03 * tc[1]
            -1.90139225e-06 * tc[2]
            +3.71185986e-10 * tc[3]
            -2.87908305e-14 * tc[4];
        /*species 8: N2 */
        species[8] =
            +2.92664000e+00
            +1.48797680e-03 * tc[1]
            -5.68476000e-07 * tc[2]
            +1.00970380e-10 * tc[3]
            -6.75335100e-15 * tc[4];
    }
    return;
}


/*compute the g/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void gibbs(amrex::Real * species, amrex::Real *  tc)
{

    /*temperature */
    amrex::Real T = tc[1];
    amrex::Real invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            -9.179351730000000e+02 * invT
            +1.661320882000000e+00
            -2.344331120000000e+00 * tc[0]
            -3.990260375000000e-03 * tc[1]
            +3.246358500000000e-06 * tc[2]
            -1.679767450000000e-09 * tc[3]
            +3.688058805000000e-13 * tc[4];
        /*species 1: H */
        species[1] =
            +2.547365990000000e+04 * invT
            +2.946682853000000e+00
            -2.500000000000000e+00 * tc[0]
            -3.526664095000000e-13 * tc[1]
            +3.326532733333333e-16 * tc[2]
            -1.917346933333333e-19 * tc[3]
            +4.638661660000000e-23 * tc[4];
        /*species 2: O */
        species[2] =
            +2.912225920000000e+04 * invT
            +1.116333640000000e+00
            -3.168267100000000e+00 * tc[0]
            +1.639659420000000e-03 * tc[1]
            -1.107177326666667e-06 * tc[2]
            +5.106721866666666e-10 * tc[3]
            -1.056329855000000e-13 * tc[4];
        /*species 3: O2 */
        species[3] =
            -1.063943560000000e+03 * invT
            +1.247806300000001e-01
            -3.782456360000000e+00 * tc[0]
            +1.498367080000000e-03 * tc[1]
            -1.641217001666667e-06 * tc[2]
            +8.067745908333334e-10 * tc[3]
            -1.621864185000000e-13 * tc[4];
        /*species 4: OH */
        species[4] =
            +3.615080560000000e+03 * invT
            +4.095940888000000e+00
            -3.992015430000000e+00 * tc[0]
            +1.200658760000000e-03 * tc[1]
            -7.696564016666666e-07 * tc[2]
            +3.234277775000000e-10 * tc[3]
            -6.820573500000000e-14 * tc[4];
        /*species 5: H2O */
        species[5] =
            -3.029372670000000e+04 * invT
            +5.047672768000000e+00
            -4.198640560000000e+00 * tc[0]
            +1.018217050000000e-03 * tc[1]
            -1.086733685000000e-06 * tc[2]
            +4.573308850000000e-10 * tc[3]
            -8.859890850000000e-14 * tc[4];
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
            -1.770258210000000e+04 * invT
            +8.410619499999998e-01
            -4.276112690000000e+00 * tc[0]
            +2.714112085000000e-04 * tc[1]
            -2.788928350000000e-06 * tc[2]
            +1.798090108333333e-09 * tc[3]
            -4.312271815000000e-13 * tc[4];
        /*species 8: N2 */
        species[8] =
            -1.020899900000000e+03 * invT
            -6.516950000000001e-01
            -3.298677000000000e+00 * tc[0]
            -7.041202000000000e-04 * tc[1]
            +6.605369999999999e-07 * tc[2]
            -4.701262500000001e-10 * tc[3]
            +1.222427000000000e-13 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            -9.501589220000000e+02 * invT
            +6.542302510000000e+00
            -3.337279200000000e+00 * tc[0]
            +2.470123655000000e-05 * tc[1]
            -8.324279633333333e-08 * tc[2]
            +1.496386616666667e-11 * tc[3]
            -1.001276880000000e-15 * tc[4];
        /*species 1: H */
        species[1] =
            +2.547365990000000e+04 * invT
            +2.946682924000000e+00
            -2.500000010000000e+00 * tc[0]
            +1.154214865000000e-11 * tc[1]
            -2.692699133333334e-15 * tc[2]
            +3.945960291666667e-19 * tc[3]
            -2.490986785000000e-23 * tc[4];
        /*species 2: O */
        species[2] =
            +2.921757910000000e+04 * invT
            -2.214917859999999e+00
            -2.569420780000000e+00 * tc[0]
            +4.298705685000000e-05 * tc[1]
            -6.991409816666667e-09 * tc[2]
            +8.348149916666666e-13 * tc[3]
            -6.141684549999999e-17 * tc[4];
        /*species 3: O2 */
        species[3] =
            -1.088457720000000e+03 * invT
            -2.170693450000000e+00
            -3.282537840000000e+00 * tc[0]
            -7.415437700000000e-04 * tc[1]
            +1.263277781666667e-07 * tc[2]
            -1.745587958333333e-11 * tc[3]
            +1.083588970000000e-15 * tc[4];
        /*species 4: OH */
        species[4] =
            +3.858657000000000e+03 * invT
            -1.383808430000000e+00
            -3.092887670000000e+00 * tc[0]
            -2.742148580000000e-04 * tc[1]
            -2.108420466666667e-08 * tc[2]
            +7.328846300000000e-12 * tc[3]
            -5.870618800000000e-16 * tc[4];
        /*species 5: H2O */
        species[5] =
            -3.000429710000000e+04 * invT
            -1.932777610000000e+00
            -3.033992490000000e+00 * tc[0]
            -1.088459020000000e-03 * tc[1]
            +2.734541966666666e-08 * tc[2]
            +8.086832250000000e-12 * tc[3]
            -8.410049600000000e-16 * tc[4];
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
            -1.786178770000000e+04 * invT
            +1.248846229999999e+00
            -4.165002850000000e+00 * tc[0]
            -2.454158470000000e-03 * tc[1]
            +3.168987083333333e-07 * tc[2]
            -3.093216550000000e-11 * tc[3]
            +1.439541525000000e-15 * tc[4];
        /*species 8: N2 */
        species[8] =
            -9.227977000000000e+02 * invT
            -3.053888000000000e+00
            -2.926640000000000e+00 * tc[0]
            -7.439884000000000e-04 * tc[1]
            +9.474600000000001e-08 * tc[2]
            -8.414198333333333e-12 * tc[3]
            +3.376675500000000e-16 * tc[4];
    }
    return;
}


/*compute the a/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void helmholtz(amrex::Real * species, amrex::Real *  tc)
{

    /*temperature */
    amrex::Real T = tc[1];
    amrex::Real invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            -9.17935173e+02 * invT
            +6.61320882e-01
            -2.34433112e+00 * tc[0]
            -3.99026037e-03 * tc[1]
            +3.24635850e-06 * tc[2]
            -1.67976745e-09 * tc[3]
            +3.68805881e-13 * tc[4];
        /*species 1: H */
        species[1] =
            +2.54736599e+04 * invT
            +1.94668285e+00
            -2.50000000e+00 * tc[0]
            -3.52666409e-13 * tc[1]
            +3.32653273e-16 * tc[2]
            -1.91734693e-19 * tc[3]
            +4.63866166e-23 * tc[4];
        /*species 2: O */
        species[2] =
            +2.91222592e+04 * invT
            +1.16333640e-01
            -3.16826710e+00 * tc[0]
            +1.63965942e-03 * tc[1]
            -1.10717733e-06 * tc[2]
            +5.10672187e-10 * tc[3]
            -1.05632985e-13 * tc[4];
        /*species 3: O2 */
        species[3] =
            -1.06394356e+03 * invT
            -8.75219370e-01
            -3.78245636e+00 * tc[0]
            +1.49836708e-03 * tc[1]
            -1.64121700e-06 * tc[2]
            +8.06774591e-10 * tc[3]
            -1.62186418e-13 * tc[4];
        /*species 4: OH */
        species[4] =
            +3.61508056e+03 * invT
            +3.09594089e+00
            -3.99201543e+00 * tc[0]
            +1.20065876e-03 * tc[1]
            -7.69656402e-07 * tc[2]
            +3.23427778e-10 * tc[3]
            -6.82057350e-14 * tc[4];
        /*species 5: H2O */
        species[5] =
            -3.02937267e+04 * invT
            +4.04767277e+00
            -4.19864056e+00 * tc[0]
            +1.01821705e-03 * tc[1]
            -1.08673369e-06 * tc[2]
            +4.57330885e-10 * tc[3]
            -8.85989085e-14 * tc[4];
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
            -1.77025821e+04 * invT
            -1.58938050e-01
            -4.27611269e+00 * tc[0]
            +2.71411208e-04 * tc[1]
            -2.78892835e-06 * tc[2]
            +1.79809011e-09 * tc[3]
            -4.31227182e-13 * tc[4];
        /*species 8: N2 */
        species[8] =
            -1.02089990e+03 * invT
            -1.65169500e+00
            -3.29867700e+00 * tc[0]
            -7.04120200e-04 * tc[1]
            +6.60537000e-07 * tc[2]
            -4.70126250e-10 * tc[3]
            +1.22242700e-13 * tc[4];
    } else {
        /*species 0: H2 */
        species[0] =
            -9.50158922e+02 * invT
            +5.54230251e+00
            -3.33727920e+00 * tc[0]
            +2.47012365e-05 * tc[1]
            -8.32427963e-08 * tc[2]
            +1.49638662e-11 * tc[3]
            -1.00127688e-15 * tc[4];
        /*species 1: H */
        species[1] =
            +2.54736599e+04 * invT
            +1.94668292e+00
            -2.50000001e+00 * tc[0]
            +1.15421486e-11 * tc[1]
            -2.69269913e-15 * tc[2]
            +3.94596029e-19 * tc[3]
            -2.49098679e-23 * tc[4];
        /*species 2: O */
        species[2] =
            +2.92175791e+04 * invT
            -3.21491786e+00
            -2.56942078e+00 * tc[0]
            +4.29870569e-05 * tc[1]
            -6.99140982e-09 * tc[2]
            +8.34814992e-13 * tc[3]
            -6.14168455e-17 * tc[4];
        /*species 3: O2 */
        species[3] =
            -1.08845772e+03 * invT
            -3.17069345e+00
            -3.28253784e+00 * tc[0]
            -7.41543770e-04 * tc[1]
            +1.26327778e-07 * tc[2]
            -1.74558796e-11 * tc[3]
            +1.08358897e-15 * tc[4];
        /*species 4: OH */
        species[4] =
            +3.85865700e+03 * invT
            -2.38380843e+00
            -3.09288767e+00 * tc[0]
            -2.74214858e-04 * tc[1]
            -2.10842047e-08 * tc[2]
            +7.32884630e-12 * tc[3]
            -5.87061880e-16 * tc[4];
        /*species 5: H2O */
        species[5] =
            -3.00042971e+04 * invT
            -2.93277761e+00
            -3.03399249e+00 * tc[0]
            -1.08845902e-03 * tc[1]
            +2.73454197e-08 * tc[2]
            +8.08683225e-12 * tc[3]
            -8.41004960e-16 * tc[4];
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
            -1.78617877e+04 * invT
            +2.48846230e-01
            -4.16500285e+00 * tc[0]
            -2.45415847e-03 * tc[1]
            +3.16898708e-07 * tc[2]
            -3.09321655e-11 * tc[3]
            +1.43954153e-15 * tc[4];
        /*species 8: N2 */
        species[8] =
            -9.22797700e+02 * invT
            -4.05388800e+00
            -2.92664000e+00 * tc[0]
            -7.43988400e-04 * tc[1]
            +9.47460000e-08 * tc[2]
            -8.41419833e-12 * tc[3]
            +3.37667550e-16 * tc[4];
    }
    return;
}


/*compute the e/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void speciesInternalEnergy(amrex::Real * species, amrex::Real *  tc)
{

    /*temperature */
    amrex::Real T = tc[1];
    amrex::Real invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            +1.34433112e+00
            +3.99026037e-03 * tc[1]
            -6.49271700e-06 * tc[2]
            +5.03930235e-09 * tc[3]
            -1.47522352e-12 * tc[4]
            -9.17935173e+02 * invT;
        /*species 1: H */
        species[1] =
            +1.50000000e+00
            +3.52666409e-13 * tc[1]
            -6.65306547e-16 * tc[2]
            +5.75204080e-19 * tc[3]
            -1.85546466e-22 * tc[4]
            +2.54736599e+04 * invT;
        /*species 2: O */
        species[2] =
            +2.16826710e+00
            -1.63965942e-03 * tc[1]
            +2.21435465e-06 * tc[2]
            -1.53201656e-09 * tc[3]
            +4.22531942e-13 * tc[4]
            +2.91222592e+04 * invT;
        /*species 3: O2 */
        species[3] =
            +2.78245636e+00
            -1.49836708e-03 * tc[1]
            +3.28243400e-06 * tc[2]
            -2.42032377e-09 * tc[3]
            +6.48745674e-13 * tc[4]
            -1.06394356e+03 * invT;
        /*species 4: OH */
        species[4] =
            +2.99201543e+00
            -1.20065876e-03 * tc[1]
            +1.53931280e-06 * tc[2]
            -9.70283332e-10 * tc[3]
            +2.72822940e-13 * tc[4]
            +3.61508056e+03 * invT;
        /*species 5: H2O */
        species[5] =
            +3.19864056e+00
            -1.01821705e-03 * tc[1]
            +2.17346737e-06 * tc[2]
            -1.37199266e-09 * tc[3]
            +3.54395634e-13 * tc[4]
            -3.02937267e+04 * invT;
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
            +3.27611269e+00
            -2.71411208e-04 * tc[1]
            +5.57785670e-06 * tc[2]
            -5.39427032e-09 * tc[3]
            +1.72490873e-12 * tc[4]
            -1.77025821e+04 * invT;
        /*species 8: N2 */
        species[8] =
            +2.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 * invT;
    } else {
        /*species 0: H2 */
        species[0] =
            +2.33727920e+00
            -2.47012365e-05 * tc[1]
            +1.66485593e-07 * tc[2]
            -4.48915985e-11 * tc[3]
            +4.00510752e-15 * tc[4]
            -9.50158922e+02 * invT;
        /*species 1: H */
        species[1] =
            +1.50000001e+00
            -1.15421486e-11 * tc[1]
            +5.38539827e-15 * tc[2]
            -1.18378809e-18 * tc[3]
            +9.96394714e-23 * tc[4]
            +2.54736599e+04 * invT;
        /*species 2: O */
        species[2] =
            +1.56942078e+00
            -4.29870569e-05 * tc[1]
            +1.39828196e-08 * tc[2]
            -2.50444497e-12 * tc[3]
            +2.45667382e-16 * tc[4]
            +2.92175791e+04 * invT;
        /*species 3: O2 */
        species[3] =
            +2.28253784e+00
            +7.41543770e-04 * tc[1]
            -2.52655556e-07 * tc[2]
            +5.23676387e-11 * tc[3]
            -4.33435588e-15 * tc[4]
            -1.08845772e+03 * invT;
        /*species 4: OH */
        species[4] =
            +2.09288767e+00
            +2.74214858e-04 * tc[1]
            +4.21684093e-08 * tc[2]
            -2.19865389e-11 * tc[3]
            +2.34824752e-15 * tc[4]
            +3.85865700e+03 * invT;
        /*species 5: H2O */
        species[5] =
            +2.03399249e+00
            +1.08845902e-03 * tc[1]
            -5.46908393e-08 * tc[2]
            -2.42604967e-11 * tc[3]
            +3.36401984e-15 * tc[4]
            -3.00042971e+04 * invT;
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
            +3.16500285e+00
            +2.45415847e-03 * tc[1]
            -6.33797417e-07 * tc[2]
            +9.27964965e-11 * tc[3]
            -5.75816610e-15 * tc[4]
            -1.78617877e+04 * invT;
        /*species 8: N2 */
        species[8] =
            +1.92664000e+00
            +7.43988400e-04 * tc[1]
            -1.89492000e-07 * tc[2]
            +2.52425950e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
    }
    return;
}


/*compute the h/(RT) at the given temperature (Eq 20) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void speciesEnthalpy(amrex::Real * species, amrex::Real *  tc)
{

    /*temperature */
    amrex::Real T = tc[1];
    amrex::Real invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            +2.34433112e+00
            +3.99026037e-03 * tc[1]
            -6.49271700e-06 * tc[2]
            +5.03930235e-09 * tc[3]
            -1.47522352e-12 * tc[4]
            -9.17935173e+02 * invT;
        /*species 1: H */
        species[1] =
            +2.50000000e+00
            +3.52666409e-13 * tc[1]
            -6.65306547e-16 * tc[2]
            +5.75204080e-19 * tc[3]
            -1.85546466e-22 * tc[4]
            +2.54736599e+04 * invT;
        /*species 2: O */
        species[2] =
            +3.16826710e+00
            -1.63965942e-03 * tc[1]
            +2.21435465e-06 * tc[2]
            -1.53201656e-09 * tc[3]
            +4.22531942e-13 * tc[4]
            +2.91222592e+04 * invT;
        /*species 3: O2 */
        species[3] =
            +3.78245636e+00
            -1.49836708e-03 * tc[1]
            +3.28243400e-06 * tc[2]
            -2.42032377e-09 * tc[3]
            +6.48745674e-13 * tc[4]
            -1.06394356e+03 * invT;
        /*species 4: OH */
        species[4] =
            +3.99201543e+00
            -1.20065876e-03 * tc[1]
            +1.53931280e-06 * tc[2]
            -9.70283332e-10 * tc[3]
            +2.72822940e-13 * tc[4]
            +3.61508056e+03 * invT;
        /*species 5: H2O */
        species[5] =
            +4.19864056e+00
            -1.01821705e-03 * tc[1]
            +2.17346737e-06 * tc[2]
            -1.37199266e-09 * tc[3]
            +3.54395634e-13 * tc[4]
            -3.02937267e+04 * invT;
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
            +4.27611269e+00
            -2.71411208e-04 * tc[1]
            +5.57785670e-06 * tc[2]
            -5.39427032e-09 * tc[3]
            +1.72490873e-12 * tc[4]
            -1.77025821e+04 * invT;
        /*species 8: N2 */
        species[8] =
            +3.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 * invT;
    } else {
        /*species 0: H2 */
        species[0] =
            +3.33727920e+00
            -2.47012365e-05 * tc[1]
            +1.66485593e-07 * tc[2]
            -4.48915985e-11 * tc[3]
            +4.00510752e-15 * tc[4]
            -9.50158922e+02 * invT;
        /*species 1: H */
        species[1] =
            +2.50000001e+00
            -1.15421486e-11 * tc[1]
            +5.38539827e-15 * tc[2]
            -1.18378809e-18 * tc[3]
            +9.96394714e-23 * tc[4]
            +2.54736599e+04 * invT;
        /*species 2: O */
        species[2] =
            +2.56942078e+00
            -4.29870569e-05 * tc[1]
            +1.39828196e-08 * tc[2]
            -2.50444497e-12 * tc[3]
            +2.45667382e-16 * tc[4]
            +2.92175791e+04 * invT;
        /*species 3: O2 */
        species[3] =
            +3.28253784e+00
            +7.41543770e-04 * tc[1]
            -2.52655556e-07 * tc[2]
            +5.23676387e-11 * tc[3]
            -4.33435588e-15 * tc[4]
            -1.08845772e+03 * invT;
        /*species 4: OH */
        species[4] =
            +3.09288767e+00
            +2.74214858e-04 * tc[1]
            +4.21684093e-08 * tc[2]
            -2.19865389e-11 * tc[3]
            +2.34824752e-15 * tc[4]
            +3.85865700e+03 * invT;
        /*species 5: H2O */
        species[5] =
            +3.03399249e+00
            +1.08845902e-03 * tc[1]
            -5.46908393e-08 * tc[2]
            -2.42604967e-11 * tc[3]
            +3.36401984e-15 * tc[4]
            -3.00042971e+04 * invT;
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
            +4.16500285e+00
            +2.45415847e-03 * tc[1]
            -6.33797417e-07 * tc[2]
            +9.27964965e-11 * tc[3]
            -5.75816610e-15 * tc[4]
            -1.78617877e+04 * invT;
        /*species 8: N2 */
        species[8] =
            +2.92664000e+00
            +7.43988400e-04 * tc[1]
            -1.89492000e-07 * tc[2]
            +2.52425950e-11 * tc[3]
            -1.35067020e-15 * tc[4]
            -9.22797700e+02 * invT;
    }
    return;
}


/*compute the S/R at the given temperature (Eq 21) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void speciesEntropy(amrex::Real * species, amrex::Real *  tc)
{

    /*temperature */
    amrex::Real T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            +2.34433112e+00 * tc[0]
            +7.98052075e-03 * tc[1]
            -9.73907550e-06 * tc[2]
            +6.71906980e-09 * tc[3]
            -1.84402940e-12 * tc[4]
            +6.83010238e-01 ;
        /*species 1: H */
        species[1] =
            +2.50000000e+00 * tc[0]
            +7.05332819e-13 * tc[1]
            -9.97959820e-16 * tc[2]
            +7.66938773e-19 * tc[3]
            -2.31933083e-22 * tc[4]
            -4.46682853e-01 ;
        /*species 2: O */
        species[2] =
            +3.16826710e+00 * tc[0]
            -3.27931884e-03 * tc[1]
            +3.32153198e-06 * tc[2]
            -2.04268875e-09 * tc[3]
            +5.28164927e-13 * tc[4]
            +2.05193346e+00 ;
        /*species 3: O2 */
        species[3] =
            +3.78245636e+00 * tc[0]
            -2.99673416e-03 * tc[1]
            +4.92365101e-06 * tc[2]
            -3.22709836e-09 * tc[3]
            +8.10932092e-13 * tc[4]
            +3.65767573e+00 ;
        /*species 4: OH */
        species[4] =
            +3.99201543e+00 * tc[0]
            -2.40131752e-03 * tc[1]
            +2.30896920e-06 * tc[2]
            -1.29371111e-09 * tc[3]
            +3.41028675e-13 * tc[4]
            -1.03925458e-01 ;
        /*species 5: H2O */
        species[5] =
            +4.19864056e+00 * tc[0]
            -2.03643410e-03 * tc[1]
            +3.26020105e-06 * tc[2]
            -1.82932354e-09 * tc[3]
            +4.42994543e-13 * tc[4]
            -8.49032208e-01 ;
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
            +4.27611269e+00 * tc[0]
            -5.42822417e-04 * tc[1]
            +8.36678505e-06 * tc[2]
            -7.19236043e-09 * tc[3]
            +2.15613591e-12 * tc[4]
            +3.43505074e+00 ;
        /*species 8: N2 */
        species[8] =
            +3.29867700e+00 * tc[0]
            +1.40824040e-03 * tc[1]
            -1.98161100e-06 * tc[2]
            +1.88050500e-09 * tc[3]
            -6.11213500e-13 * tc[4]
            +3.95037200e+00 ;
    } else {
        /*species 0: H2 */
        species[0] =
            +3.33727920e+00 * tc[0]
            -4.94024731e-05 * tc[1]
            +2.49728389e-07 * tc[2]
            -5.98554647e-11 * tc[3]
            +5.00638440e-15 * tc[4]
            -3.20502331e+00 ;
        /*species 1: H */
        species[1] =
            +2.50000001e+00 * tc[0]
            -2.30842973e-11 * tc[1]
            +8.07809740e-15 * tc[2]
            -1.57838412e-18 * tc[3]
            +1.24549339e-22 * tc[4]
            -4.46682914e-01 ;
        /*species 2: O */
        species[2] =
            +2.56942078e+00 * tc[0]
            -8.59741137e-05 * tc[1]
            +2.09742295e-08 * tc[2]
            -3.33925997e-12 * tc[3]
            +3.07084227e-16 * tc[4]
            +4.78433864e+00 ;
        /*species 3: O2 */
        species[3] =
            +3.28253784e+00 * tc[0]
            +1.48308754e-03 * tc[1]
            -3.78983334e-07 * tc[2]
            +6.98235183e-11 * tc[3]
            -5.41794485e-15 * tc[4]
            +5.45323129e+00 ;
        /*species 4: OH */
        species[4] =
            +3.09288767e+00 * tc[0]
            +5.48429716e-04 * tc[1]
            +6.32526140e-08 * tc[2]
            -2.93153852e-11 * tc[3]
            +2.93530940e-15 * tc[4]
            +4.47669610e+00 ;
        /*species 5: H2O */
        species[5] =
            +3.03399249e+00 * tc[0]
            +2.17691804e-03 * tc[1]
            -8.20362590e-08 * tc[2]
            -3.23473290e-11 * tc[3]
            +4.20502480e-15 * tc[4]
            +4.96677010e+00 ;
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
            +4.16500285e+00 * tc[0]
            +4.90831694e-03 * tc[1]
            -9.50696125e-07 * tc[2]
            +1.23728662e-10 * tc[3]
            -7.19770763e-15 * tc[4]
            +2.91615662e+00 ;
        /*species 8: N2 */
        species[8] =
            +2.92664000e+00 * tc[0]
            +1.48797680e-03 * tc[1]
            -2.84238000e-07 * tc[2]
            +3.36567933e-11 * tc[3]
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
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
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
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
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
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
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
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real cvor[9]; /* temporary storage */
    amrex::Real imw[9];

    get_imw(imw);

    cv_R(cvor, tc);

    /*multiply by y/molecularweight */
    result += cvor[0]*y[0]*imw[0]; /*H2 */
    result += cvor[1]*y[1]*imw[1]; /*H */
    result += cvor[2]*y[2]*imw[2]; /*O */
    result += cvor[3]*y[3]*imw[3]; /*O2 */
    result += cvor[4]*y[4]*imw[4]; /*OH */
    result += cvor[5]*y[5]*imw[5]; /*H2O */
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
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
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
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
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
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
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
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real ums[9]; /* temporary energy array */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    amrex::Real imw[9];

    get_imw(imw);

    speciesInternalEnergy(ums, tc);

    /*perform dot product + scaling by wt */
    result += y[0]*ums[0]*imw[0]; /*H2 */
    result += y[1]*ums[1]*imw[1]; /*H */
    result += y[2]*ums[2]*imw[2]; /*O */
    result += y[3]*ums[3]*imw[3]; /*O2 */
    result += y[4]*ums[4]*imw[4]; /*OH */
    result += y[5]*ums[5]*imw[5]; /*H2O */
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
    amrex::Real tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
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
    amrex::Real tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real sor[9]; /* temporary storage */
    amrex::Real x[9]; /* need a ytx conversion */
    amrex::Real YOW = 0; /*See Eq 4, 6 in CK Manual */
    amrex::Real imw[9];

    get_imw(imw);

    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*H2 */
    YOW += y[1]*imw[1]; /*H */
    YOW += y[2]*imw[2]; /*O */
    YOW += y[3]*imw[3]; /*O2 */
    YOW += y[4]*imw[4]; /*OH */
    YOW += y[5]*imw[5]; /*H2O */
    YOW += y[6]*imw[6]; /*HO2 */
    YOW += y[7]*imw[7]; /*H2O2 */
    YOW += y[8]*imw[8]; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(2.015940*YOW); 
    x[1] = y[1]/(1.007970*YOW); 
    x[2] = y[2]/(15.999400*YOW); 
    x[3] = y[3]/(31.998800*YOW); 
    x[4] = y[4]/(17.007370*YOW); 
    x[5] = y[5]/(18.015340*YOW); 
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
    ai[0] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[0],2.0) / (pow(2.015880,2.0) * 12.964000); 
    bi[0] = 0.08664 * Rcst * Tci[0] / (2.015880 * 12.964000); 
    acentric_i[0] = -0.219000 ;

    /*species 1: H */
    Tci[1] = 1.316 * EPS[1] ; 
    ai[1] = (5.55 * pow(avogadro,2.0) * EPS[1]*boltzmann * pow(1e-8*SIG[1],3.0) ) / (pow(wt[1],2.0)); 
    bi[1] = 0.855 * avogadro * pow(1e-8*SIG[1],3.0) / (wt[1]); 
    acentric_i[1] = 0.0 ;

    /*species 2: O */
    Tci[2] = 1.316 * EPS[2] ; 
    ai[2] = (5.55 * pow(avogadro,2.0) * EPS[2]*boltzmann * pow(1e-8*SIG[2],3.0) ) / (pow(wt[2],2.0)); 
    bi[2] = 0.855 * avogadro * pow(1e-8*SIG[2],3.0) / (wt[2]); 
    acentric_i[2] = 0.0 ;

    /*species 3: O2 */
    /*Imported from NIST */
    Tci[3] = 154.581000 ; 
    ai[3] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[3],2.0) / (pow(31.998800,2.0) * 50.430466); 
    bi[3] = 0.08664 * Rcst * Tci[3] / (31.998800 * 50.430466); 
    acentric_i[3] = 0.022200 ;

    /*species 4: OH */
    Tci[4] = 1.316 * EPS[4] ; 
    ai[4] = (5.55 * pow(avogadro,2.0) * EPS[4]*boltzmann * pow(1e-8*SIG[4],3.0) ) / (pow(wt[4],2.0)); 
    bi[4] = 0.855 * avogadro * pow(1e-8*SIG[4],3.0) / (wt[4]); 
    acentric_i[4] = 0.0 ;

    /*species 5: H2O */
    /*Imported from NIST */
    Tci[5] = 647.096000 ; 
    ai[5] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[5],2.0) / (pow(18.015340,2.0) * 220.640000); 
    bi[5] = 0.08664 * Rcst * Tci[5] / (18.015340 * 220.640000); 
    acentric_i[5] = 0.344300 ;

    /*species 6: HO2 */
    Tci[6] = 1.316 * EPS[6] ; 
    ai[6] = (5.55 * pow(avogadro,2.0) * EPS[6]*boltzmann * pow(1e-8*SIG[6],3.0) ) / (pow(wt[6],2.0)); 
    bi[6] = 0.855 * avogadro * pow(1e-8*SIG[6],3.0) / (wt[6]); 
    acentric_i[6] = 0.0 ;

    /*species 7: H2O2 */
    Tci[7] = 1.316 * EPS[7] ; 
    ai[7] = (5.55 * pow(avogadro,2.0) * EPS[7]*boltzmann * pow(1e-8*SIG[7],3.0) ) / (pow(wt[7],2.0)); 
    bi[7] = 0.855 * avogadro * pow(1e-8*SIG[7],3.0) / (wt[7]); 
    acentric_i[7] = 0.0 ;

    /*species 8: N2 */
    /*Imported from NIST */
    Tci[8] = 126.192000 ; 
    ai[8] = 1e6 * 0.42748 * pow(Rcst,2.0) * pow(Tci[8],2.0) / (pow(28.013400,2.0) * 33.958000); 
    bi[8] = 0.08664 * Rcst * Tci[8] / (28.013400 * 33.958000); 
    acentric_i[8] = 0.037200 ;

    return;
}


/*Compute P = rhoRT/W(x) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKPX(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  x, amrex::Real *  P)
{
    amrex::Real XW = 0;/* To hold mean molecular wt */
    XW += x[0]*2.015940; /*H2 */
    XW += x[1]*1.007970; /*H */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*18.015340; /*H2O */
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
    YOW += y[1]*imw[1]; /*H */
    YOW += y[2]*imw[2]; /*O */
    YOW += y[3]*imw[3]; /*O2 */
    YOW += y[4]*imw[4]; /*OH */
    YOW += y[5]*imw[5]; /*H2O */
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
    W += c[1]*1.007970; /*H */
    W += c[2]*15.999400; /*O */
    W += c[3]*31.998800; /*O2 */
    W += c[4]*17.007370; /*OH */
    W += c[5]*18.015340; /*H2O */
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
    XW += x[1]*1.007970; /*H */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*18.015340; /*H2O */
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
    W += c[1]*1.007970; /*H */
    W += c[2]*15.999400; /*O */
    W += c[3]*31.998800; /*O2 */
    W += c[4]*17.007370; /*OH */
    W += c[5]*18.015340; /*H2O */
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
    XW += x[1]*1.007970; /*H */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*18.015340; /*H2O */
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
    W += c[1]*1.007970; /*H */
    W += c[2]*15.999400; /*O */
    W += c[3]*31.998800; /*O2 */
    W += c[4]*17.007370; /*OH */
    W += c[5]*18.015340; /*H2O */
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
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpor, tc);
}


/*get H/RT as a function of T  */
/*for all species (Eq 20) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKHORT(amrex::Real *  T, amrex::Real *  hort)
{
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEnthalpy(hort, tc);
}


/*get S/R as a function of T  */
/*for all species (Eq 21) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKSOR(amrex::Real *  T, amrex::Real *  sor)
{
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
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
    XW += x[1]*1.007970; /*H */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*18.015340; /*H2O */
    XW += x[6]*33.006770; /*HO2 */
    XW += x[7]*34.014740; /*H2O2 */
    XW += x[8]*28.013400; /*N2 */
    /*Now compute conversion */
    amrex::Real XWinv = 1.0/XW;
    y[0] = x[0]*2.015940*XWinv; 
    y[1] = x[1]*1.007970*XWinv; 
    y[2] = x[2]*15.999400*XWinv; 
    y[3] = x[3]*31.998800*XWinv; 
    y[4] = x[4]*17.007370*XWinv; 
    y[5] = x[5]*18.015340*XWinv; 
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
    XW += x[1]*1.007970; /*H */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*18.015340; /*H2O */
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
    CW += c[1]*1.007970; /*H */
    CW += c[2]*15.999400; /*O */
    CW += c[3]*31.998800; /*O2 */
    CW += c[4]*17.007370; /*OH */
    CW += c[5]*18.015340; /*H2O */
    CW += c[6]*33.006770; /*HO2 */
    CW += c[7]*34.014740; /*H2O2 */
    CW += c[8]*28.013400; /*N2 */
    /*Now compute conversion */
    amrex::Real CWinv = 1.0/CW;
    y[0] = c[0]*2.015940*CWinv; 
    y[1] = c[1]*1.007970*CWinv; 
    y[2] = c[2]*15.999400*CWinv; 
    y[3] = c[3]*31.998800*CWinv; 
    y[4] = c[4]*17.007370*CWinv; 
    y[5] = c[5]*18.015340*CWinv; 
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
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
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
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
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
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
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
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
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
    amrex::Real tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
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
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cv_R(cvms, tc);
    /*multiply by R/molecularweight */
    cvms[0] *= 4.124360158612479e+07; /*H2 */
    cvms[1] *= 8.248720317224957e+07; /*H */
    cvms[2] *= 5.196734013871295e+06; /*O */
    cvms[3] *= 2.598367006935648e+06; /*O2 */
    cvms[4] *= 4.888740950630956e+06; /*OH */
    cvms[5] *= 4.615212712140454e+06; /*H2O */
    cvms[6] *= 2.519017346487778e+06; /*HO2 */
    cvms[7] *= 2.444370475315478e+06; /*H2O2 */
    cvms[8] *= 2.968030520448514e+06; /*N2 */
}


/*Returns the specific heats at constant pressure */
/*in mass units (Eq. 26) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCPMS(amrex::Real *  T,  amrex::Real *  cpms)
{
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpms, tc);
    /*multiply by R/molecularweight */
    cpms[0] *= 4.124360158612479e+07; /*H2 */
    cpms[1] *= 8.248720317224957e+07; /*H */
    cpms[2] *= 5.196734013871295e+06; /*O */
    cpms[3] *= 2.598367006935648e+06; /*O2 */
    cpms[4] *= 4.888740950630956e+06; /*OH */
    cpms[5] *= 4.615212712140454e+06; /*H2O */
    cpms[6] *= 2.519017346487778e+06; /*HO2 */
    cpms[7] *= 2.444370475315478e+06; /*H2O2 */
    cpms[8] *= 2.968030520448514e+06; /*N2 */
}


/*Returns internal energy in mass units (Eq 30.) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKUMS(amrex::Real *  T,  amrex::Real *  ums)
{
    amrex::Real tT = *T; /*temporary temperature */
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
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
    amrex::Real tc[] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
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
    amrex::Real tc[] = { log(tT), tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    speciesEntropy(sms, tc);
    /*multiply by R/molecularweight */
    sms[0] *= 4.124360158612479e+07; /*H2 */
    sms[1] *= 8.248720317224957e+07; /*H */
    sms[2] *= 5.196734013871295e+06; /*O */
    sms[3] *= 2.598367006935648e+06; /*O2 */
    sms[4] *= 4.888740950630956e+06; /*OH */
    sms[5] *= 4.615212712140454e+06; /*H2O */
    sms[6] *= 2.519017346487778e+06; /*HO2 */
    sms[7] *= 2.444370475315478e+06; /*H2O2 */
    sms[8] *= 2.968030520448514e+06; /*N2 */
}


/* GPU version of productionRate: no more use of thermo namespace vectors  */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void comp_qfqr(amrex::Real *  qf, amrex::Real * qr, amrex::Real * sc, amrex::Real * qss_sc, amrex::Real * tc, amrex::Real invT)
{

    /*reaction 21: 2.000000 OH (+M) <=> H2O2 (+M) */
    qf[0] = pow(sc[4], 2.000000);
    qr[0] = sc[7];

    /*reaction 1: 2.000000 O + M <=> O2 + M */
    qf[1] = pow(sc[2], 2.000000);
    qr[1] = sc[3];

    /*reaction 2: O + H + M <=> OH + M */
    qf[2] = sc[1]*sc[2];
    qr[2] = sc[4];

    /*reaction 6: H + O2 + M <=> HO2 + M */
    qf[3] = sc[1]*sc[3];
    qr[3] = sc[6];

    /*reaction 11: 2.000000 H + M <=> H2 + M */
    qf[4] = pow(sc[1], 2.000000);
    qr[4] = sc[0];

    /*reaction 14: H + OH + M <=> H2O + M */
    qf[5] = sc[1]*sc[4];
    qr[5] = sc[5];

    /*reaction 3: O + H2 <=> H + OH */
    qf[6] = sc[0]*sc[2];
    qr[6] = sc[1]*sc[4];

    /*reaction 4: O + HO2 <=> OH + O2 */
    qf[7] = sc[2]*sc[6];
    qr[7] = sc[3]*sc[4];

    /*reaction 5: O + H2O2 <=> OH + HO2 */
    qf[8] = sc[2]*sc[7];
    qr[8] = sc[4]*sc[6];

    /*reaction 7: H + 2.000000 O2 <=> HO2 + O2 */
    qf[9] = sc[1]*pow(sc[3], 2.000000);
    qr[9] = sc[3]*sc[6];

    /*reaction 8: H + O2 + H2O <=> HO2 + H2O */
    qf[10] = sc[1]*sc[3]*sc[5];
    qr[10] = sc[5]*sc[6];

    /*reaction 9: H + O2 + N2 <=> HO2 + N2 */
    qf[11] = sc[1]*sc[3]*sc[8];
    qr[11] = sc[6]*sc[8];

    /*reaction 10: H + O2 <=> O + OH */
    qf[12] = sc[1]*sc[3];
    qr[12] = sc[2]*sc[4];

    /*reaction 12: 2.000000 H + H2 <=> 2.000000 H2 */
    qf[13] = sc[0]*pow(sc[1], 2.000000);
    qr[13] = pow(sc[0], 2.000000);

    /*reaction 13: 2.000000 H + H2O <=> H2 + H2O */
    qf[14] = pow(sc[1], 2.000000)*sc[5];
    qr[14] = sc[0]*sc[5];

    /*reaction 15: H + HO2 <=> O + H2O */
    qf[15] = sc[1]*sc[6];
    qr[15] = sc[2]*sc[5];

    /*reaction 16: H + HO2 <=> O2 + H2 */
    qf[16] = sc[1]*sc[6];
    qr[16] = sc[0]*sc[3];

    /*reaction 17: H + HO2 <=> 2.000000 OH */
    qf[17] = sc[1]*sc[6];
    qr[17] = pow(sc[4], 2.000000);

    /*reaction 18: H + H2O2 <=> HO2 + H2 */
    qf[18] = sc[1]*sc[7];
    qr[18] = sc[0]*sc[6];

    /*reaction 19: H + H2O2 <=> OH + H2O */
    qf[19] = sc[1]*sc[7];
    qr[19] = sc[4]*sc[5];

    /*reaction 20: OH + H2 <=> H + H2O */
    qf[20] = sc[0]*sc[4];
    qr[20] = sc[1]*sc[5];

    /*reaction 22: 2.000000 OH <=> O + H2O */
    qf[21] = pow(sc[4], 2.000000);
    qr[21] = sc[2]*sc[5];

    /*reaction 23: OH + HO2 <=> O2 + H2O */
    qf[22] = sc[4]*sc[6];
    qr[22] = sc[3]*sc[5];

    /*reaction 24: OH + H2O2 <=> HO2 + H2O */
    qf[23] = sc[4]*sc[7];
    qr[23] = sc[5]*sc[6];

    /*reaction 25: OH + H2O2 <=> HO2 + H2O */
    qf[24] = sc[4]*sc[7];
    qr[24] = sc[5]*sc[6];

    /*reaction 26: 2.000000 HO2 <=> O2 + H2O2 */
    qf[25] = pow(sc[6], 2.000000);
    qr[25] = sc[3]*sc[7];

    /*reaction 27: 2.000000 HO2 <=> O2 + H2O2 */
    qf[26] = pow(sc[6], 2.000000);
    qr[26] = sc[3]*sc[7];

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

    /* Evaluate the kfs */
    amrex::Real k_f, k_r, Corr;
    amrex::Real redP, F, logPred, logFcent, troe_c, troe_n, troe, F_troe;

    // (0):  2.000000 O + M <=> O2 + M
    k_f = 1.0000000000000002e-12 * 1.2e+17 
               * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = mixture + ( 2.3999999999999999 - 1)*sc[0] + ( 15.4 - 1)*sc[5];
    qf[1] *= Corr * k_f;
    qr[1] *= Corr * k_f / (exp(2.000000*g_RT[2] - g_RT[3]) * refCinv);
    // (1):  O + H + M <=> OH + M
    k_f = 1.0000000000000002e-12 * 5e+17 
               * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5];
    qf[2] *= Corr * k_f;
    qr[2] *= Corr * k_f / (exp(g_RT[1] + g_RT[2] - g_RT[4]) * refCinv);
    // (2):  O + H2 <=> H + OH
    k_f = 1.0000000000000002e-06 * 50000 
               * exp(2.6699999999999999 * tc[0] - 0.50321666580471969 * (6290) * invT);
    Corr  = 1.0;
    qf[6] *= Corr * k_f;
    qr[6] *= Corr * k_f / exp(g_RT[0] - g_RT[1] + g_RT[2] - g_RT[4]);
    // (3):  O + HO2 <=> OH + O2
    k_f = 1.0000000000000002e-06 * 20000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[7] *= Corr * k_f;
    qr[7] *= Corr * k_f / exp(g_RT[2] - g_RT[3] - g_RT[4] + g_RT[6]);
    // (4):  O + H2O2 <=> OH + HO2
    k_f = 1.0000000000000002e-06 * 9630000 
               * exp(2 * tc[0] - 0.50321666580471969 * (4000) * invT);
    Corr  = 1.0;
    qf[8] *= Corr * k_f;
    qr[8] *= Corr * k_f / exp(g_RT[2] - g_RT[4] - g_RT[6] + g_RT[7]);
    // (5):  H + O2 + M <=> HO2 + M
    k_f = 1.0000000000000002e-12 * 2.8e+18 
               * exp(-0.85999999999999999 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = mixture + ( 0 - 1)*sc[3] + ( 0 - 1)*sc[5] + ( 0 - 1)*sc[8];
    qf[3] *= Corr * k_f;
    qr[3] *= Corr * k_f / (exp(g_RT[1] + g_RT[3] - g_RT[6]) * refCinv);
    // (6):  H + 2.000000 O2 <=> HO2 + O2
    k_f = 1.0000000000000002e-12 * 3e+20 
               * exp(-1.72 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[9] *= Corr * k_f;
    qr[9] *= Corr * k_f / (exp(g_RT[1] + 2.000000*g_RT[3] - g_RT[3] - g_RT[6]) * refCinv);
    // (7):  H + O2 + H2O <=> HO2 + H2O
    k_f = 1.0000000000000002e-12 * 9.38e+18 
               * exp(-0.76000000000000001 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[10] *= Corr * k_f;
    qr[10] *= Corr * k_f / (exp(g_RT[1] + g_RT[3] + g_RT[5] - g_RT[5] - g_RT[6]) * refCinv);
    // (8):  H + O2 + N2 <=> HO2 + N2
    k_f = 1.0000000000000002e-12 * 3.75e+20 
               * exp(-1.72 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[11] *= Corr * k_f;
    qr[11] *= Corr * k_f / (exp(g_RT[1] + g_RT[3] - g_RT[6] + g_RT[8] - g_RT[8]) * refCinv);
    // (9):  H + O2 <=> O + OH
    k_f = 1.0000000000000002e-06 * 83000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (14413) * invT);
    Corr  = 1.0;
    qf[12] *= Corr * k_f;
    qr[12] *= Corr * k_f / exp(g_RT[1] - g_RT[2] + g_RT[3] - g_RT[4]);
    // (10):  2.000000 H + M <=> H2 + M
    k_f = 1.0000000000000002e-12 * 1e+18 
               * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = mixture + ( 0 - 1)*sc[0] + ( 0 - 1)*sc[5];
    qf[4] *= Corr * k_f;
    qr[4] *= Corr * k_f / (exp(-g_RT[0] + 2.000000*g_RT[1]) * refCinv);
    // (11):  2.000000 H + H2 <=> 2.000000 H2
    k_f = 1.0000000000000002e-12 * 90000000000000000 
               * exp(-0.59999999999999998 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[13] *= Corr * k_f;
    qr[13] *= Corr * k_f / (exp(g_RT[0] - 2.000000*g_RT[0] + 2.000000*g_RT[1]) * refCinv);
    // (12):  2.000000 H + H2O <=> H2 + H2O
    k_f = 1.0000000000000002e-12 * 6e+19 
               * exp(-1.25 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = 1.0;
    qf[14] *= Corr * k_f;
    qr[14] *= Corr * k_f / (exp(-g_RT[0] + 2.000000*g_RT[1] + g_RT[5] - g_RT[5]) * refCinv);
    // (13):  H + OH + M <=> H2O + M
    k_f = 1.0000000000000002e-12 * 2.2e+22 
               * exp(-2 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = mixture + ( 0.72999999999999998 - 1)*sc[0] + ( 3.6499999999999999 - 1)*sc[5];
    qf[5] *= Corr * k_f;
    qr[5] *= Corr * k_f / (exp(g_RT[1] + g_RT[4] - g_RT[5]) * refCinv);
    // (14):  H + HO2 <=> O + H2O
    k_f = 1.0000000000000002e-06 * 3970000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (671) * invT);
    Corr  = 1.0;
    qf[15] *= Corr * k_f;
    qr[15] *= Corr * k_f / exp(g_RT[1] - g_RT[2] - g_RT[5] + g_RT[6]);
    // (15):  H + HO2 <=> O2 + H2
    k_f = 1.0000000000000002e-06 * 28000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (1068) * invT);
    Corr  = 1.0;
    qf[16] *= Corr * k_f;
    qr[16] *= Corr * k_f / exp(-g_RT[0] + g_RT[1] - g_RT[3] + g_RT[6]);
    // (16):  H + HO2 <=> 2.000000 OH
    k_f = 1.0000000000000002e-06 * 134000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (635) * invT);
    Corr  = 1.0;
    qf[17] *= Corr * k_f;
    qr[17] *= Corr * k_f / exp(g_RT[1] - 2.000000*g_RT[4] + g_RT[6]);
    // (17):  H + H2O2 <=> HO2 + H2
    k_f = 1.0000000000000002e-06 * 12100000 
               * exp(2 * tc[0] - 0.50321666580471969 * (5200) * invT);
    Corr  = 1.0;
    qf[18] *= Corr * k_f;
    qr[18] *= Corr * k_f / exp(-g_RT[0] + g_RT[1] - g_RT[6] + g_RT[7]);
    // (18):  H + H2O2 <=> OH + H2O
    k_f = 1.0000000000000002e-06 * 10000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (3600) * invT);
    Corr  = 1.0;
    qf[19] *= Corr * k_f;
    qr[19] *= Corr * k_f / exp(g_RT[1] - g_RT[4] - g_RT[5] + g_RT[7]);
    // (19):  OH + H2 <=> H + H2O
    k_f = 1.0000000000000002e-06 * 216000000 
               * exp(1.51 * tc[0] - 0.50321666580471969 * (3430) * invT);
    Corr  = 1.0;
    qf[20] *= Corr * k_f;
    qr[20] *= Corr * k_f / exp(g_RT[0] - g_RT[1] + g_RT[4] - g_RT[5]);
    // (20):  2.000000 OH (+M) <=> H2O2 (+M)
    k_f = 1.0000000000000002e-06 * 74000000000000 
               * exp(-0.37 * tc[0] - 0.50321666580471969 * (0) * invT);
    Corr  = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5];
    redP = Corr / k_f * 1e-12 * 2.3e+18 
               * exp(-0.90000000000000002  * tc[0] - 0.50321666580471969  * (-1700) *invT);
    F = redP / (1.0 + redP);
    logPred = log10(redP);
    logFcent = log10(
        (1.-0.73460000000000003)*exp(-tc[1] / 94) 
        + 0.73460000000000003 * exp(-tc[1]/1756)  
        + exp(-5182 * invT));
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troe = (troe_c + logPred) / (troe_n - .14*(troe_c + logPred));
    F_troe = pow(10., logFcent / (1.0 + troe*troe));
    Corr = F * F_troe;
    qf[0] *= Corr * k_f;
    qr[0] *= Corr * k_f / (exp(2.000000*g_RT[4] - g_RT[7]) * refCinv);
    // (21):  2.000000 OH <=> O + H2O
    k_f = 1.0000000000000002e-06 * 35700 
               * exp(2.3999999999999999 * tc[0] - 0.50321666580471969 * (-2110) * invT);
    Corr  = 1.0;
    qf[21] *= Corr * k_f;
    qr[21] *= Corr * k_f / exp(-g_RT[2] + 2.000000*g_RT[4] - g_RT[5]);
    // (22):  OH + HO2 <=> O2 + H2O
    k_f = 1.0000000000000002e-06 * 29000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (-500) * invT);
    Corr  = 1.0;
    qf[22] *= Corr * k_f;
    qr[22] *= Corr * k_f / exp(-g_RT[3] + g_RT[4] - g_RT[5] + g_RT[6]);
    // (23):  OH + H2O2 <=> HO2 + H2O
    k_f = 1.0000000000000002e-06 * 1750000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (320) * invT);
    Corr  = 1.0;
    qf[23] *= Corr * k_f;
    qr[23] *= Corr * k_f / exp(g_RT[4] - g_RT[5] - g_RT[6] + g_RT[7]);
    // (24):  OH + H2O2 <=> HO2 + H2O
    k_f = 1.0000000000000002e-06 * 580000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (9560) * invT);
    Corr  = 1.0;
    qf[24] *= Corr * k_f;
    qr[24] *= Corr * k_f / exp(g_RT[4] - g_RT[5] - g_RT[6] + g_RT[7]);
    // (25):  2.000000 HO2 <=> O2 + H2O2
    k_f = 1.0000000000000002e-06 * 130000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (-1630) * invT);
    Corr  = 1.0;
    qf[25] *= Corr * k_f;
    qr[25] *= Corr * k_f / exp(-g_RT[3] + 2.000000*g_RT[6] - g_RT[7]);
    // (26):  2.000000 HO2 <=> O2 + H2O2
    k_f = 1.0000000000000002e-06 * 420000000000000 
               * exp(0 * tc[0] - 0.50321666580471969 * (12000) * invT);
    Corr  = 1.0;
    qf[26] *= Corr * k_f;
    qr[26] *= Corr * k_f / exp(-g_RT[3] + 2.000000*g_RT[6] - g_RT[7]);


    return;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void  productionRate(amrex::Real * wdot, amrex::Real * sc, amrex::Real T)
{
    amrex::Real tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    amrex::Real invT = 1.0 / tc[1];

    amrex::Real qdot, q_f[27], q_r[27];
    amrex::Real sc_qss[1];
    comp_qfqr(q_f, q_r, sc, sc_qss, tc, invT);

    for (int i = 0; i < 9; ++i) {
        wdot[i] = 0.0;
    }

    qdot = q_f[0]-q_r[0];
    wdot[4] -= 2.000000 * qdot;
    wdot[7] += qdot;

    qdot = q_f[1]-q_r[1];
    wdot[2] -= 2.000000 * qdot;
    wdot[3] += qdot;

    qdot = q_f[2]-q_r[2];
    wdot[1] -= qdot;
    wdot[2] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[3]-q_r[3];
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[6] += qdot;

    qdot = q_f[4]-q_r[4];
    wdot[0] += qdot;
    wdot[1] -= 2.000000 * qdot;

    qdot = q_f[5]-q_r[5];
    wdot[1] -= qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[6]-q_r[6];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[2] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[7]-q_r[7];
    wdot[2] -= qdot;
    wdot[3] += qdot;
    wdot[4] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[8]-q_r[8];
    wdot[2] -= qdot;
    wdot[4] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[9]-q_r[9];
    wdot[1] -= qdot;
    wdot[3] -= 2.000000 * qdot;
    wdot[3] += qdot;
    wdot[6] += qdot;

    qdot = q_f[10]-q_r[10];
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[5] -= qdot;
    wdot[5] += qdot;
    wdot[6] += qdot;

    qdot = q_f[11]-q_r[11];
    wdot[1] -= qdot;
    wdot[3] -= qdot;
    wdot[6] += qdot;
    wdot[8] -= qdot;
    wdot[8] += qdot;

    qdot = q_f[12]-q_r[12];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[3] -= qdot;
    wdot[4] += qdot;

    qdot = q_f[13]-q_r[13];
    wdot[0] -= qdot;
    wdot[0] += 2.000000 * qdot;
    wdot[1] -= 2.000000 * qdot;

    qdot = q_f[14]-q_r[14];
    wdot[0] += qdot;
    wdot[1] -= 2.000000 * qdot;
    wdot[5] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[15]-q_r[15];
    wdot[1] -= qdot;
    wdot[2] += qdot;
    wdot[5] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[16]-q_r[16];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[3] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[17]-q_r[17];
    wdot[1] -= qdot;
    wdot[4] += 2.000000 * qdot;
    wdot[6] -= qdot;

    qdot = q_f[18]-q_r[18];
    wdot[0] += qdot;
    wdot[1] -= qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[19]-q_r[19];
    wdot[1] -= qdot;
    wdot[4] += qdot;
    wdot[5] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[20]-q_r[20];
    wdot[0] -= qdot;
    wdot[1] += qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;

    qdot = q_f[21]-q_r[21];
    wdot[2] += qdot;
    wdot[4] -= 2.000000 * qdot;
    wdot[5] += qdot;

    qdot = q_f[22]-q_r[22];
    wdot[3] += qdot;
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[6] -= qdot;

    qdot = q_f[23]-q_r[23];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[24]-q_r[24];
    wdot[4] -= qdot;
    wdot[5] += qdot;
    wdot[6] += qdot;
    wdot[7] -= qdot;

    qdot = q_f[25]-q_r[25];
    wdot[3] += qdot;
    wdot[6] -= 2.000000 * qdot;
    wdot[7] += qdot;

    qdot = q_f[26]-q_r[26];
    wdot[3] += qdot;
    wdot[6] -= 2.000000 * qdot;
    wdot[7] += qdot;

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
    YOW += y[1]*imw[1]; /*H */
    YOW += y[2]*imw[2]; /*O */
    YOW += y[3]*imw[3]; /*O2 */
    YOW += y[4]*imw[4]; /*OH */
    YOW += y[5]*imw[5]; /*H2O */
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
    XW += x[1]*1.007970; /*H */
    XW += x[2]*15.999400; /*O */
    XW += x[3]*31.998800; /*O2 */
    XW += x[4]*17.007370; /*OH */
    XW += x[5]*18.015340; /*H2O */
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
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void dcvpRdT(amrex::Real * species, amrex::Real *  tc)
{

    /*temperature */
    amrex::Real T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: H2 */
        species[0] =
            +7.98052075e-03
            -3.89563020e-05 * tc[1]
            +6.04716282e-08 * tc[2]
            -2.95044704e-11 * tc[3];
        /*species 1: H */
        species[1] =
            +7.05332819e-13
            -3.99183928e-15 * tc[1]
            +6.90244896e-18 * tc[2]
            -3.71092933e-21 * tc[3];
        /*species 2: O */
        species[2] =
            -3.27931884e-03
            +1.32861279e-05 * tc[1]
            -1.83841987e-08 * tc[2]
            +8.45063884e-12 * tc[3];
        /*species 3: O2 */
        species[3] =
            -2.99673416e-03
            +1.96946040e-05 * tc[1]
            -2.90438853e-08 * tc[2]
            +1.29749135e-11 * tc[3];
        /*species 4: OH */
        species[4] =
            -2.40131752e-03
            +9.23587682e-06 * tc[1]
            -1.16434000e-08 * tc[2]
            +5.45645880e-12 * tc[3];
        /*species 5: H2O */
        species[5] =
            -2.03643410e-03
            +1.30408042e-05 * tc[1]
            -1.64639119e-08 * tc[2]
            +7.08791268e-12 * tc[3];
        /*species 6: HO2 */
        species[6] =
            -4.74912051e-03
            +4.23165782e-05 * tc[1]
            -7.28291682e-08 * tc[2]
            +3.71690050e-11 * tc[3];
        /*species 7: H2O2 */
        species[7] =
            -5.42822417e-04
            +3.34671402e-05 * tc[1]
            -6.47312439e-08 * tc[2]
            +3.44981745e-11 * tc[3];
        /*species 8: N2 */
        species[8] =
            +1.40824040e-03
            -7.92644400e-06 * tc[1]
            +1.69245450e-08 * tc[2]
            -9.77941600e-12 * tc[3];
    } else {
        /*species 0: H2 */
        species[0] =
            -4.94024731e-05
            +9.98913556e-07 * tc[1]
            -5.38699182e-10 * tc[2]
            +8.01021504e-14 * tc[3];
        /*species 1: H */
        species[1] =
            -2.30842973e-11
            +3.23123896e-14 * tc[1]
            -1.42054571e-17 * tc[2]
            +1.99278943e-21 * tc[3];
        /*species 2: O */
        species[2] =
            -8.59741137e-05
            +8.38969178e-08 * tc[1]
            -3.00533397e-11 * tc[2]
            +4.91334764e-15 * tc[3];
        /*species 3: O2 */
        species[3] =
            +1.48308754e-03
            -1.51593334e-06 * tc[1]
            +6.28411665e-10 * tc[2]
            -8.66871176e-14 * tc[3];
        /*species 4: OH */
        species[4] =
            +5.48429716e-04
            +2.53010456e-07 * tc[1]
            -2.63838467e-10 * tc[2]
            +4.69649504e-14 * tc[3];
        /*species 5: H2O */
        species[5] =
            +2.17691804e-03
            -3.28145036e-07 * tc[1]
            -2.91125961e-10 * tc[2]
            +6.72803968e-14 * tc[3];
        /*species 6: HO2 */
        species[6] =
            +2.23982013e-03
            -1.26731630e-06 * tc[1]
            +3.42739110e-10 * tc[2]
            -4.31634140e-14 * tc[3];
        /*species 7: H2O2 */
        species[7] =
            +4.90831694e-03
            -3.80278450e-06 * tc[1]
            +1.11355796e-09 * tc[2]
            -1.15163322e-13 * tc[3];
        /*species 8: N2 */
        species[8] =
            +1.48797680e-03
            -1.13695200e-06 * tc[1]
            +3.02911140e-10 * tc[2]
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

    amrex::Real tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
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
    /*reaction 1: 2.000000 OH (+M) <=> H2O2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5];
    /* forward */
    phi_f = pow(sc[4], 2.000000);
    k_f = 1.0000000000000002e-06 * 74000000000000
                * exp(-0.37 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -0.37 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* pressure-fall-off */
    k_0 = 2.3e+18 * exp(-0.90000000000000002 * tc[0] - 0.50321666580471969 * (-1700) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -0.90000000000000002 * invT + 0.50321666580471969 * (-1700) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-0.73460000000000003)*exp(-T/94);
    Fcent2 = 0.73460000000000003 * exp(-T/1756);
    Fcent3 = exp(-5182 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/94
        -Fcent2/1756
        + Fcent3*5182*invT2);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[7];
    Kc = refCinv * exp(2.000000*g_RT[4] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[4]) + (h_RT[7]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[4] -= 2 * q; /* OH */
    wdot[7] += q; /* H2O2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = 0.0;
    dqdc[0] = 2*dcdc_fac;
    dqdc[1] = dcdc_fac;
    dqdc[2] = dcdc_fac;
    dqdc[3] = dcdc_fac;
    dqdc[4] = dcdc_fac + k_f*2.000000*sc[4];
    dqdc[5] = 6*dcdc_fac;
    dqdc[6] = dcdc_fac;
    dqdc[7] = dcdc_fac - k_r;
    dqdc[8] = dcdc_fac;
    for (int k=0; k<9; k++) {
        J[10*k+4] += -2 * dqdc[k];
        J[10*k+7] += dqdc[k];
    }
    J[94] += -2 * dqdT; /* dwdot[OH]/dT */
    J[97] += dqdT; /* dwdot[H2O2]/dT */

    /*reaction 2: 2.000000 O + M <=> O2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 2.3999999999999999 - 1)*sc[0] + ( 15.4 - 1)*sc[5];
    /* forward */
    phi_f = pow(sc[2], 2.000000);
    k_f = 1.0000000000000002e-12 * 1.2e+17
                * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[3];
    Kc = refCinv * exp(2.000000*g_RT[2] - g_RT[3]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[2]) + (h_RT[3]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= 2 * q; /* O */
    wdot[3] += q; /* O2 */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    dqdc[0] = 2.3999999999999999*q_nocor;
    dqdc[1] = q_nocor;
    dqdc[2] = q_nocor + k_f*2.000000*sc[2];
    dqdc[3] = q_nocor - k_r;
    dqdc[4] = q_nocor;
    dqdc[5] = 15.4*q_nocor;
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    for (int k=0; k<9; k++) {
        J[10*k+2] += -2 * dqdc[k];
        J[10*k+3] += dqdc[k];
    }
    J[92] += -2 * dqdT; /* dwdot[O]/dT */
    J[93] += dqdT; /* dwdot[O2]/dT */

    /*reaction 3: O + H + M <=> OH + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5];
    /* forward */
    phi_f = sc[1]*sc[2];
    k_f = 1.0000000000000002e-12 * 5e+17
                * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[4];
    Kc = refCinv * exp(g_RT[1] + g_RT[2] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[2]) + (h_RT[4]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    dqdc[0] = 2*q_nocor;
    dqdc[1] = q_nocor + k_f*sc[2];
    dqdc[2] = q_nocor + k_f*sc[1];
    dqdc[3] = q_nocor;
    dqdc[4] = q_nocor - k_r;
    dqdc[5] = 6*q_nocor;
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    for (int k=0; k<9; k++) {
        J[10*k+1] -= dqdc[k];
        J[10*k+2] -= dqdc[k];
        J[10*k+4] += dqdc[k];
    }
    J[91] -= dqdT; /* dwdot[H]/dT */
    J[92] -= dqdT; /* dwdot[O]/dT */
    J[94] += dqdT; /* dwdot[OH]/dT */

    /*reaction 4: H + O2 + M <=> HO2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 0 - 1)*sc[3] + ( 0 - 1)*sc[5] + ( 0 - 1)*sc[8];
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = 1.0000000000000002e-12 * 2.8e+18
                * exp(-0.85999999999999999 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -0.85999999999999999 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[6];
    Kc = refCinv * exp(g_RT[1] + g_RT[3] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3]) + (h_RT[6]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    dqdc[0] = q_nocor;
    dqdc[1] = q_nocor + k_f*sc[3];
    dqdc[2] = q_nocor;
    dqdc[3] =  + k_f*sc[1];
    dqdc[4] = q_nocor;
    dqdc[5] = 0.0;
    dqdc[6] = q_nocor - k_r;
    dqdc[7] = q_nocor;
    dqdc[8] = 0.0;
    for (int k=0; k<9; k++) {
        J[10*k+1] -= dqdc[k];
        J[10*k+3] -= dqdc[k];
        J[10*k+6] += dqdc[k];
    }
    J[91] -= dqdT; /* dwdot[H]/dT */
    J[93] -= dqdT; /* dwdot[O2]/dT */
    J[96] += dqdT; /* dwdot[HO2]/dT */

    /*reaction 5: 2.000000 H + M <=> H2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 0 - 1)*sc[0] + ( 0 - 1)*sc[5];
    /* forward */
    phi_f = pow(sc[1], 2.000000);
    k_f = 1.0000000000000002e-12 * 1e+18
                * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[0];
    Kc = refCinv * exp(-g_RT[0] + 2.000000*g_RT[1]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[1]) + (h_RT[0]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= 2 * q; /* H */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    dqdc[0] =  - k_r;
    dqdc[1] = q_nocor + k_f*2.000000*sc[1];
    dqdc[2] = q_nocor;
    dqdc[3] = q_nocor;
    dqdc[4] = q_nocor;
    dqdc[5] = 0.0;
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    for (int k=0; k<9; k++) {
        J[10*k+0] += dqdc[k];
        J[10*k+1] += -2 * dqdc[k];
    }
    J[90] += dqdT; /* dwdot[H2]/dT */
    J[91] += -2 * dqdT; /* dwdot[H]/dT */

    /*reaction 6: H + OH + M <=> H2O + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 0.72999999999999998 - 1)*sc[0] + ( 3.6499999999999999 - 1)*sc[5];
    /* forward */
    phi_f = sc[1]*sc[4];
    k_f = 1.0000000000000002e-12 * 2.2e+22
                * exp(-2 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -2 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[5];
    Kc = refCinv * exp(g_RT[1] + g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[4]) + (h_RT[5]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    dqdc[0] = 0.72999999999999998*q_nocor;
    dqdc[1] = q_nocor + k_f*sc[4];
    dqdc[2] = q_nocor;
    dqdc[3] = q_nocor;
    dqdc[4] = q_nocor + k_f*sc[1];
    dqdc[5] = 3.6499999999999999*q_nocor - k_r;
    dqdc[6] = q_nocor;
    dqdc[7] = q_nocor;
    dqdc[8] = q_nocor;
    for (int k=0; k<9; k++) {
        J[10*k+1] -= dqdc[k];
        J[10*k+4] -= dqdc[k];
        J[10*k+5] += dqdc[k];
    }
    J[91] -= dqdT; /* dwdot[H]/dT */
    J[94] -= dqdT; /* dwdot[OH]/dT */
    J[95] += dqdT; /* dwdot[H2O]/dT */

    /*reaction 7: O + H2 <=> H + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[2];
    k_f = 1.0000000000000002e-06 * 50000
                * exp(2.6699999999999999 * tc[0] - 0.50321666580471969 * (6290) * invT);
    dlnkfdT = 2.6699999999999999 * invT + 0.50321666580471969 *  (6290)  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[4];
    Kc = exp(g_RT[0] - g_RT[1] + g_RT[2] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[2]) + (h_RT[1] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[1] += q; /* H */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[2];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += dqdci;                /* dwdot[H]/d[H2] */
    J[2] -= dqdci;                /* dwdot[O]/d[H2] */
    J[4] += dqdci;                /* dwdot[OH]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[4];
    J[10] -= dqdci;               /* dwdot[H2]/d[H] */
    J[11] += dqdci;               /* dwdot[H]/d[H] */
    J[12] -= dqdci;               /* dwdot[O]/d[H] */
    J[14] += dqdci;               /* dwdot[OH]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[0];
    J[20] -= dqdci;               /* dwdot[H2]/d[O] */
    J[21] += dqdci;               /* dwdot[H]/d[O] */
    J[22] -= dqdci;               /* dwdot[O]/d[O] */
    J[24] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[1];
    J[40] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[41] += dqdci;               /* dwdot[H]/d[OH] */
    J[42] -= dqdci;               /* dwdot[O]/d[OH] */
    J[44] += dqdci;               /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[90] -= dqdT;                /* dwdot[H2]/dT */
    J[91] += dqdT;                /* dwdot[H]/dT */
    J[92] -= dqdT;                /* dwdot[O]/dT */
    J[94] += dqdT;                /* dwdot[OH]/dT */

    /*reaction 8: O + HO2 <=> OH + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[6];
    k_f = 1.0000000000000002e-06 * 20000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[3]*sc[4];
    Kc = exp(g_RT[2] - g_RT[3] - g_RT[4] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[6]) + (h_RT[3] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[3] += q; /* O2 */
    wdot[4] += q; /* OH */
    wdot[6] -= q; /* HO2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[6];
    J[22] -= dqdci;               /* dwdot[O]/d[O] */
    J[23] += dqdci;               /* dwdot[O2]/d[O] */
    J[24] += dqdci;               /* dwdot[OH]/d[O] */
    J[26] -= dqdci;               /* dwdot[HO2]/d[O] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[4];
    J[32] -= dqdci;               /* dwdot[O]/d[O2] */
    J[33] += dqdci;               /* dwdot[O2]/d[O2] */
    J[34] += dqdci;               /* dwdot[OH]/d[O2] */
    J[36] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[3];
    J[42] -= dqdci;               /* dwdot[O]/d[OH] */
    J[43] += dqdci;               /* dwdot[O2]/d[OH] */
    J[44] += dqdci;               /* dwdot[OH]/d[OH] */
    J[46] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[2];
    J[62] -= dqdci;               /* dwdot[O]/d[HO2] */
    J[63] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[64] += dqdci;               /* dwdot[OH]/d[HO2] */
    J[66] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[92] -= dqdT;                /* dwdot[O]/dT */
    J[93] += dqdT;                /* dwdot[O2]/dT */
    J[94] += dqdT;                /* dwdot[OH]/dT */
    J[96] -= dqdT;                /* dwdot[HO2]/dT */

    /*reaction 9: O + H2O2 <=> OH + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[7];
    k_f = 1.0000000000000002e-06 * 9630000
                * exp(2 * tc[0] - 0.50321666580471969 * (4000) * invT);
    dlnkfdT = 2 * invT + 0.50321666580471969 *  (4000)  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[6];
    Kc = exp(g_RT[2] - g_RT[4] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[7]) + (h_RT[4] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[7];
    J[22] -= dqdci;               /* dwdot[O]/d[O] */
    J[24] += dqdci;               /* dwdot[OH]/d[O] */
    J[26] += dqdci;               /* dwdot[HO2]/d[O] */
    J[27] -= dqdci;               /* dwdot[H2O2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[6];
    J[42] -= dqdci;               /* dwdot[O]/d[OH] */
    J[44] += dqdci;               /* dwdot[OH]/d[OH] */
    J[46] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[47] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[4];
    J[62] -= dqdci;               /* dwdot[O]/d[HO2] */
    J[64] += dqdci;               /* dwdot[OH]/d[HO2] */
    J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[67] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[2];
    J[72] -= dqdci;               /* dwdot[O]/d[H2O2] */
    J[74] += dqdci;               /* dwdot[OH]/d[H2O2] */
    J[76] += dqdci;               /* dwdot[HO2]/d[H2O2] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[92] -= dqdT;                /* dwdot[O]/dT */
    J[94] += dqdT;                /* dwdot[OH]/dT */
    J[96] += dqdT;                /* dwdot[HO2]/dT */
    J[97] -= dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 10: H + 2.000000 O2 <=> HO2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*pow(sc[3], 2.000000);
    k_f = 1.0000000000000002e-12 * 3e+20
                * exp(-1.72 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1.72 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[3]*sc[6];
    Kc = refCinv * exp(g_RT[1] + 2.000000*g_RT[3] - g_RT[3] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + 2.000000*h_RT[3]) + (h_RT[3] + h_RT[6]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*pow(sc[3], 2.000000);
    J[11] -= dqdci;               /* dwdot[H]/d[H] */
    J[13] -= dqdci;               /* dwdot[O2]/d[H] */
    J[16] += dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1]*2.000000*sc[3] - k_r*sc[6];
    J[31] -= dqdci;               /* dwdot[H]/d[O2] */
    J[33] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[36] += dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[3];
    J[61] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[63] -= dqdci;               /* dwdot[O2]/d[HO2] */
    J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[91] -= dqdT;                /* dwdot[H]/dT */
    J[93] -= dqdT;                /* dwdot[O2]/dT */
    J[96] += dqdT;                /* dwdot[HO2]/dT */

    /*reaction 11: H + O2 + H2O <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3]*sc[5];
    k_f = 1.0000000000000002e-12 * 9.38e+18
                * exp(-0.76000000000000001 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -0.76000000000000001 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[6];
    Kc = refCinv * exp(g_RT[1] + g_RT[3] + g_RT[5] - g_RT[5] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3] + h_RT[5]) + (h_RT[5] + h_RT[6]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[3]*sc[5];
    J[11] -= dqdci;               /* dwdot[H]/d[H] */
    J[13] -= dqdci;               /* dwdot[O2]/d[H] */
    J[16] += dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1]*sc[5];
    J[31] -= dqdci;               /* dwdot[H]/d[O2] */
    J[33] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[36] += dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[1]*sc[3] - k_r*sc[6];
    J[51] -= dqdci;               /* dwdot[H]/d[H2O] */
    J[53] -= dqdci;               /* dwdot[O2]/d[H2O] */
    J[56] += dqdci;               /* dwdot[HO2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[5];
    J[61] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[63] -= dqdci;               /* dwdot[O2]/d[HO2] */
    J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[91] -= dqdT;                /* dwdot[H]/dT */
    J[93] -= dqdT;                /* dwdot[O2]/dT */
    J[96] += dqdT;                /* dwdot[HO2]/dT */

    /*reaction 12: H + O2 + N2 <=> HO2 + N2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3]*sc[8];
    k_f = 1.0000000000000002e-12 * 3.75e+20
                * exp(-1.72 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1.72 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[6]*sc[8];
    Kc = refCinv * exp(g_RT[1] + g_RT[3] - g_RT[6] + g_RT[8] - g_RT[8]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3] + h_RT[8]) + (h_RT[6] + h_RT[8]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[3]*sc[8];
    J[11] -= dqdci;               /* dwdot[H]/d[H] */
    J[13] -= dqdci;               /* dwdot[O2]/d[H] */
    J[16] += dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1]*sc[8];
    J[31] -= dqdci;               /* dwdot[H]/d[O2] */
    J[33] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[36] += dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[8];
    J[61] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[63] -= dqdci;               /* dwdot[O2]/d[HO2] */
    J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/d[N2] */
    dqdci =  + k_f*sc[1]*sc[3] - k_r*sc[6];
    J[81] -= dqdci;               /* dwdot[H]/d[N2] */
    J[83] -= dqdci;               /* dwdot[O2]/d[N2] */
    J[86] += dqdci;               /* dwdot[HO2]/d[N2] */
    /* d()/dT */
    J[91] -= dqdT;                /* dwdot[H]/dT */
    J[93] -= dqdT;                /* dwdot[O2]/dT */
    J[96] += dqdT;                /* dwdot[HO2]/dT */

    /*reaction 13: H + O2 <=> O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = 1.0000000000000002e-06 * 83000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (14413) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (14413)  * invT2;
    /* reverse */
    phi_r = sc[2]*sc[4];
    Kc = exp(g_RT[1] - g_RT[2] + g_RT[3] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3]) + (h_RT[2] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[2] += q; /* O */
    wdot[3] -= q; /* O2 */
    wdot[4] += q; /* OH */
    /* d()/d[H] */
    dqdci =  + k_f*sc[3];
    J[11] -= dqdci;               /* dwdot[H]/d[H] */
    J[12] += dqdci;               /* dwdot[O]/d[H] */
    J[13] -= dqdci;               /* dwdot[O2]/d[H] */
    J[14] += dqdci;               /* dwdot[OH]/d[H] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[4];
    J[21] -= dqdci;               /* dwdot[H]/d[O] */
    J[22] += dqdci;               /* dwdot[O]/d[O] */
    J[23] -= dqdci;               /* dwdot[O2]/d[O] */
    J[24] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1];
    J[31] -= dqdci;               /* dwdot[H]/d[O2] */
    J[32] += dqdci;               /* dwdot[O]/d[O2] */
    J[33] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[34] += dqdci;               /* dwdot[OH]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[2];
    J[41] -= dqdci;               /* dwdot[H]/d[OH] */
    J[42] += dqdci;               /* dwdot[O]/d[OH] */
    J[43] -= dqdci;               /* dwdot[O2]/d[OH] */
    J[44] += dqdci;               /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[91] -= dqdT;                /* dwdot[H]/dT */
    J[92] += dqdT;                /* dwdot[O]/dT */
    J[93] -= dqdT;                /* dwdot[O2]/dT */
    J[94] += dqdT;                /* dwdot[OH]/dT */

    /*reaction 14: 2.000000 H + H2 <=> 2.000000 H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*pow(sc[1], 2.000000);
    k_f = 1.0000000000000002e-12 * 90000000000000000
                * exp(-0.59999999999999998 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -0.59999999999999998 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = pow(sc[0], 2.000000);
    Kc = refCinv * exp(g_RT[0] - 2.000000*g_RT[0] + 2.000000*g_RT[1]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + 2.000000*h_RT[1]) + (2.000000*h_RT[0]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= 2 * q; /* H */
    /* d()/d[H2] */
    dqdci =  + k_f*pow(sc[1], 2.000000) - k_r*2.000000*sc[0];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += -2 * dqdci;           /* dwdot[H]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[0]*2.000000*sc[1];
    J[10] += dqdci;               /* dwdot[H2]/d[H] */
    J[11] += -2 * dqdci;          /* dwdot[H]/d[H] */
    /* d()/dT */
    J[90] += dqdT;                /* dwdot[H2]/dT */
    J[91] += -2 * dqdT;           /* dwdot[H]/dT */

    /*reaction 15: 2.000000 H + H2O <=> H2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[1], 2.000000)*sc[5];
    k_f = 1.0000000000000002e-12 * 6e+19
                * exp(-1.25 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1.25 * invT + 0.50321666580471969 *  (0)  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[5];
    Kc = refCinv * exp(-g_RT[0] + 2.000000*g_RT[1] + g_RT[5] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[1] + h_RT[5]) + (h_RT[0] + h_RT[5]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= 2 * q; /* H */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[5];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += -2 * dqdci;           /* dwdot[H]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*2.000000*sc[1]*sc[5];
    J[10] += dqdci;               /* dwdot[H2]/d[H] */
    J[11] += -2 * dqdci;          /* dwdot[H]/d[H] */
    /* d()/d[H2O] */
    dqdci =  + k_f*pow(sc[1], 2.000000) - k_r*sc[0];
    J[50] += dqdci;               /* dwdot[H2]/d[H2O] */
    J[51] += -2 * dqdci;          /* dwdot[H]/d[H2O] */
    /* d()/dT */
    J[90] += dqdT;                /* dwdot[H2]/dT */
    J[91] += -2 * dqdT;           /* dwdot[H]/dT */

    /*reaction 16: H + HO2 <=> O + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[6];
    k_f = 1.0000000000000002e-06 * 3970000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (671) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (671)  * invT2;
    /* reverse */
    phi_r = sc[2]*sc[5];
    Kc = exp(g_RT[1] - g_RT[2] - g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[6]) + (h_RT[2] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[2] += q; /* O */
    wdot[5] += q; /* H2O */
    wdot[6] -= q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[11] -= dqdci;               /* dwdot[H]/d[H] */
    J[12] += dqdci;               /* dwdot[O]/d[H] */
    J[15] += dqdci;               /* dwdot[H2O]/d[H] */
    J[16] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[5];
    J[21] -= dqdci;               /* dwdot[H]/d[O] */
    J[22] += dqdci;               /* dwdot[O]/d[O] */
    J[25] += dqdci;               /* dwdot[H2O]/d[O] */
    J[26] -= dqdci;               /* dwdot[HO2]/d[O] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[2];
    J[51] -= dqdci;               /* dwdot[H]/d[H2O] */
    J[52] += dqdci;               /* dwdot[O]/d[H2O] */
    J[55] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[56] -= dqdci;               /* dwdot[HO2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[1];
    J[61] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[62] += dqdci;               /* dwdot[O]/d[HO2] */
    J[65] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[66] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[91] -= dqdT;                /* dwdot[H]/dT */
    J[92] += dqdT;                /* dwdot[O]/dT */
    J[95] += dqdT;                /* dwdot[H2O]/dT */
    J[96] -= dqdT;                /* dwdot[HO2]/dT */

    /*reaction 17: H + HO2 <=> O2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[6];
    k_f = 1.0000000000000002e-06 * 28000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (1068) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (1068)  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[3];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[3] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[6]) + (h_RT[0] + h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[3] += q; /* O2 */
    wdot[6] -= q; /* HO2 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[3];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[3] += dqdci;                /* dwdot[O2]/d[H2] */
    J[6] -= dqdci;                /* dwdot[HO2]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[10] += dqdci;               /* dwdot[H2]/d[H] */
    J[11] -= dqdci;               /* dwdot[H]/d[H] */
    J[13] += dqdci;               /* dwdot[O2]/d[H] */
    J[16] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[0];
    J[30] += dqdci;               /* dwdot[H2]/d[O2] */
    J[31] -= dqdci;               /* dwdot[H]/d[O2] */
    J[33] += dqdci;               /* dwdot[O2]/d[O2] */
    J[36] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[1];
    J[60] += dqdci;               /* dwdot[H2]/d[HO2] */
    J[61] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[63] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[66] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[90] += dqdT;                /* dwdot[H2]/dT */
    J[91] -= dqdT;                /* dwdot[H]/dT */
    J[93] += dqdT;                /* dwdot[O2]/dT */
    J[96] -= dqdT;                /* dwdot[HO2]/dT */

    /*reaction 18: H + HO2 <=> 2.000000 OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[6];
    k_f = 1.0000000000000002e-06 * 134000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (635) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (635)  * invT2;
    /* reverse */
    phi_r = pow(sc[4], 2.000000);
    Kc = exp(g_RT[1] - 2.000000*g_RT[4] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[6]) + (2.000000*h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[4] += 2 * q; /* OH */
    wdot[6] -= q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[11] -= dqdci;               /* dwdot[H]/d[H] */
    J[14] += 2 * dqdci;           /* dwdot[OH]/d[H] */
    J[16] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*2.000000*sc[4];
    J[41] -= dqdci;               /* dwdot[H]/d[OH] */
    J[44] += 2 * dqdci;           /* dwdot[OH]/d[OH] */
    J[46] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[1];
    J[61] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[64] += 2 * dqdci;           /* dwdot[OH]/d[HO2] */
    J[66] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[91] -= dqdT;                /* dwdot[H]/dT */
    J[94] += 2 * dqdT;            /* dwdot[OH]/dT */
    J[96] -= dqdT;                /* dwdot[HO2]/dT */

    /*reaction 19: H + H2O2 <=> HO2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[7];
    k_f = 1.0000000000000002e-06 * 12100000
                * exp(2 * tc[0] - 0.50321666580471969 * (5200) * invT);
    dlnkfdT = 2 * invT + 0.50321666580471969 *  (5200)  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[6];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[7]) + (h_RT[0] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[6];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[6] += dqdci;                /* dwdot[HO2]/d[H2] */
    J[7] -= dqdci;                /* dwdot[H2O2]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[10] += dqdci;               /* dwdot[H2]/d[H] */
    J[11] -= dqdci;               /* dwdot[H]/d[H] */
    J[16] += dqdci;               /* dwdot[HO2]/d[H] */
    J[17] -= dqdci;               /* dwdot[H2O2]/d[H] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[0];
    J[60] += dqdci;               /* dwdot[H2]/d[HO2] */
    J[61] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[67] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[1];
    J[70] += dqdci;               /* dwdot[H2]/d[H2O2] */
    J[71] -= dqdci;               /* dwdot[H]/d[H2O2] */
    J[76] += dqdci;               /* dwdot[HO2]/d[H2O2] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[90] += dqdT;                /* dwdot[H2]/dT */
    J[91] -= dqdT;                /* dwdot[H]/dT */
    J[96] += dqdT;                /* dwdot[HO2]/dT */
    J[97] -= dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 20: H + H2O2 <=> OH + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[7];
    k_f = 1.0000000000000002e-06 * 10000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (3600) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (3600)  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[5];
    Kc = exp(g_RT[1] - g_RT[4] - g_RT[5] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[7]) + (h_RT[4] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[4] += q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[11] -= dqdci;               /* dwdot[H]/d[H] */
    J[14] += dqdci;               /* dwdot[OH]/d[H] */
    J[15] += dqdci;               /* dwdot[H2O]/d[H] */
    J[17] -= dqdci;               /* dwdot[H2O2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[5];
    J[41] -= dqdci;               /* dwdot[H]/d[OH] */
    J[44] += dqdci;               /* dwdot[OH]/d[OH] */
    J[45] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[47] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[4];
    J[51] -= dqdci;               /* dwdot[H]/d[H2O] */
    J[54] += dqdci;               /* dwdot[OH]/d[H2O] */
    J[55] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[57] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[1];
    J[71] -= dqdci;               /* dwdot[H]/d[H2O2] */
    J[74] += dqdci;               /* dwdot[OH]/d[H2O2] */
    J[75] += dqdci;               /* dwdot[H2O]/d[H2O2] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[91] -= dqdT;                /* dwdot[H]/dT */
    J[94] += dqdT;                /* dwdot[OH]/dT */
    J[95] += dqdT;                /* dwdot[H2O]/dT */
    J[97] -= dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 21: OH + H2 <=> H + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[4];
    k_f = 1.0000000000000002e-06 * 216000000
                * exp(1.51 * tc[0] - 0.50321666580471969 * (3430) * invT);
    dlnkfdT = 1.51 * invT + 0.50321666580471969 *  (3430)  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[5];
    Kc = exp(g_RT[0] - g_RT[1] + g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[4]) + (h_RT[1] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[1] += q; /* H */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[4];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += dqdci;                /* dwdot[H]/d[H2] */
    J[4] -= dqdci;                /* dwdot[OH]/d[H2] */
    J[5] += dqdci;                /* dwdot[H2O]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[5];
    J[10] -= dqdci;               /* dwdot[H2]/d[H] */
    J[11] += dqdci;               /* dwdot[H]/d[H] */
    J[14] -= dqdci;               /* dwdot[OH]/d[H] */
    J[15] += dqdci;               /* dwdot[H2O]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[0];
    J[40] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[41] += dqdci;               /* dwdot[H]/d[OH] */
    J[44] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[45] += dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[1];
    J[50] -= dqdci;               /* dwdot[H2]/d[H2O] */
    J[51] += dqdci;               /* dwdot[H]/d[H2O] */
    J[54] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[55] += dqdci;               /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[90] -= dqdT;                /* dwdot[H2]/dT */
    J[91] += dqdT;                /* dwdot[H]/dT */
    J[94] -= dqdT;                /* dwdot[OH]/dT */
    J[95] += dqdT;                /* dwdot[H2O]/dT */

    /*reaction 22: 2.000000 OH <=> O + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[4], 2.000000);
    k_f = 1.0000000000000002e-06 * 35700
                * exp(2.3999999999999999 * tc[0] - 0.50321666580471969 * (-2110) * invT);
    dlnkfdT = 2.3999999999999999 * invT + 0.50321666580471969 *  (-2110)  * invT2;
    /* reverse */
    phi_r = sc[2]*sc[5];
    Kc = exp(-g_RT[2] + 2.000000*g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[4]) + (h_RT[2] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* O */
    wdot[4] -= 2 * q; /* OH */
    wdot[5] += q; /* H2O */
    /* d()/d[O] */
    dqdci =  - k_r*sc[5];
    J[22] += dqdci;               /* dwdot[O]/d[O] */
    J[24] += -2 * dqdci;          /* dwdot[OH]/d[O] */
    J[25] += dqdci;               /* dwdot[H2O]/d[O] */
    /* d()/d[OH] */
    dqdci =  + k_f*2.000000*sc[4];
    J[42] += dqdci;               /* dwdot[O]/d[OH] */
    J[44] += -2 * dqdci;          /* dwdot[OH]/d[OH] */
    J[45] += dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[2];
    J[52] += dqdci;               /* dwdot[O]/d[H2O] */
    J[54] += -2 * dqdci;          /* dwdot[OH]/d[H2O] */
    J[55] += dqdci;               /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[92] += dqdT;                /* dwdot[O]/dT */
    J[94] += -2 * dqdT;           /* dwdot[OH]/dT */
    J[95] += dqdT;                /* dwdot[H2O]/dT */

    /*reaction 23: OH + HO2 <=> O2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[6];
    k_f = 1.0000000000000002e-06 * 29000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-500) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (-500)  * invT2;
    /* reverse */
    phi_r = sc[3]*sc[5];
    Kc = exp(-g_RT[3] + g_RT[4] - g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[6]) + (h_RT[3] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += q; /* O2 */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[6] -= q; /* HO2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[5];
    J[33] += dqdci;               /* dwdot[O2]/d[O2] */
    J[34] -= dqdci;               /* dwdot[OH]/d[O2] */
    J[35] += dqdci;               /* dwdot[H2O]/d[O2] */
    J[36] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[6];
    J[43] += dqdci;               /* dwdot[O2]/d[OH] */
    J[44] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[45] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[46] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[3];
    J[53] += dqdci;               /* dwdot[O2]/d[H2O] */
    J[54] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[55] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[56] -= dqdci;               /* dwdot[HO2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[4];
    J[63] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[64] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[65] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[66] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[93] += dqdT;                /* dwdot[O2]/dT */
    J[94] -= dqdT;                /* dwdot[OH]/dT */
    J[95] += dqdT;                /* dwdot[H2O]/dT */
    J[96] -= dqdT;                /* dwdot[HO2]/dT */

    /*reaction 24: OH + H2O2 <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[7];
    k_f = 1.0000000000000002e-06 * 1750000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (320) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (320)  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[6];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[7]) + (h_RT[5] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[7];
    J[44] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[45] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[46] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[47] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[6];
    J[54] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[55] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[56] += dqdci;               /* dwdot[HO2]/d[H2O] */
    J[57] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[5];
    J[64] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[65] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[67] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[4];
    J[74] -= dqdci;               /* dwdot[OH]/d[H2O2] */
    J[75] += dqdci;               /* dwdot[H2O]/d[H2O2] */
    J[76] += dqdci;               /* dwdot[HO2]/d[H2O2] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[94] -= dqdT;                /* dwdot[OH]/dT */
    J[95] += dqdT;                /* dwdot[H2O]/dT */
    J[96] += dqdT;                /* dwdot[HO2]/dT */
    J[97] -= dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 25: OH + H2O2 <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[7];
    k_f = 1.0000000000000002e-06 * 580000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (9560) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (9560)  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[6];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[7]) + (h_RT[5] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[7];
    J[44] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[45] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[46] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[47] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[6];
    J[54] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[55] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[56] += dqdci;               /* dwdot[HO2]/d[H2O] */
    J[57] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[5];
    J[64] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[65] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[67] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[4];
    J[74] -= dqdci;               /* dwdot[OH]/d[H2O2] */
    J[75] += dqdci;               /* dwdot[H2O]/d[H2O2] */
    J[76] += dqdci;               /* dwdot[HO2]/d[H2O2] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[94] -= dqdT;                /* dwdot[OH]/dT */
    J[95] += dqdT;                /* dwdot[H2O]/dT */
    J[96] += dqdT;                /* dwdot[HO2]/dT */
    J[97] -= dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 26: 2.000000 HO2 <=> O2 + H2O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[6], 2.000000);
    k_f = 1.0000000000000002e-06 * 130000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-1630) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (-1630)  * invT2;
    /* reverse */
    phi_r = sc[3]*sc[7];
    Kc = exp(-g_RT[3] + 2.000000*g_RT[6] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[6]) + (h_RT[3] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += q; /* O2 */
    wdot[6] -= 2 * q; /* HO2 */
    wdot[7] += q; /* H2O2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[7];
    J[33] += dqdci;               /* dwdot[O2]/d[O2] */
    J[36] += -2 * dqdci;          /* dwdot[HO2]/d[O2] */
    J[37] += dqdci;               /* dwdot[H2O2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*2.000000*sc[6];
    J[63] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[66] += -2 * dqdci;          /* dwdot[HO2]/d[HO2] */
    J[67] += dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[3];
    J[73] += dqdci;               /* dwdot[O2]/d[H2O2] */
    J[76] += -2 * dqdci;          /* dwdot[HO2]/d[H2O2] */
    J[77] += dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[93] += dqdT;                /* dwdot[O2]/dT */
    J[96] += -2 * dqdT;           /* dwdot[HO2]/dT */
    J[97] += dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 27: 2.000000 HO2 <=> O2 + H2O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[6], 2.000000);
    k_f = 1.0000000000000002e-06 * 420000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (12000) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  (12000)  * invT2;
    /* reverse */
    phi_r = sc[3]*sc[7];
    Kc = exp(-g_RT[3] + 2.000000*g_RT[6] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[6]) + (h_RT[3] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += q; /* O2 */
    wdot[6] -= 2 * q; /* HO2 */
    wdot[7] += q; /* H2O2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[7];
    J[33] += dqdci;               /* dwdot[O2]/d[O2] */
    J[36] += -2 * dqdci;          /* dwdot[HO2]/d[O2] */
    J[37] += dqdci;               /* dwdot[H2O2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*2.000000*sc[6];
    J[63] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[66] += -2 * dqdci;          /* dwdot[HO2]/d[HO2] */
    J[67] += dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[3];
    J[73] += dqdci;               /* dwdot[O2]/d[H2O2] */
    J[76] += -2 * dqdci;          /* dwdot[HO2]/d[H2O2] */
    J[77] += dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[93] += dqdT;                /* dwdot[O2]/dT */
    J[96] += -2 * dqdT;           /* dwdot[HO2]/dT */
    J[97] += dqdT;                /* dwdot[H2O2]/dT */

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
    /* dTdot/d[X] */
    for (int k = 0; k < 9; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 9; ++m) {
            dehmixdc += eh_RT[m]*J[k*10+m];
        }
        J[k*10+9] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[99] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;
}

/*compute an approx to the reaction Jacobian (for preconditioning) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void DWDOT_SIMPLIFIED(amrex::Real *  J, amrex::Real *  sc, amrex::Real *  Tp, int * HP)
{
    amrex::Real c[9];

    for (int k=0; k<9; k++) {
        c[k] = 1.e6 * sc[k];
    }

    aJacobian_precond(J, c, *Tp, *HP);

    /* dwdot[k]/dT */
    /* dTdot/d[X] */
    for (int k=0; k<9; k++) {
        J[90+k] *= 1.e-6;
        J[k*10+9] *= 1.e6;
    }

    return;
}

/*compute the reaction Jacobian on GPU */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void aJacobian(amrex::Real * J, amrex::Real * sc, amrex::Real T, int consP)
{

    for (int i=0; i<100; i++) {
        J[i] = 0.0;
    }

    amrex::Real wdot[9];
    for (int k=0; k<9; k++) {
        wdot[k] = 0.0;
    }

    amrex::Real tc[] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
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
    /*reaction 1: 2.000000 OH (+M) <=> H2O2 (+M) */
    /*a pressure-fall-off reaction */
    /* also 3-body */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5];
    /* forward */
    phi_f = pow(sc[4], 2.000000);
    k_f = 1.0000000000000002e-06 * 74000000000000
                * exp(-0.37 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -0.37 * invT + 0.50321666580471969 *  0  * invT2;
    /* pressure-fall-off */
    k_0 = 2.3e+18 * exp(-0.90000000000000002 * tc[0] - 0.50321666580471969 * (-1700) * invT);
    Pr = 1e-12 * alpha / k_f * k_0;
    fPr = Pr / (1.0+Pr);
    dlnk0dT = -0.90000000000000002 * invT + 0.50321666580471969 * (-1700) * invT2;
    dlogPrdT = log10e*(dlnk0dT - dlnkfdT);
    dlogfPrdT = dlogPrdT / (1.0+Pr);
    /* Troe form */
    logPr = log10(Pr);
    Fcent1 = (1.-0.73460000000000003)*exp(-T/94);
    Fcent2 = 0.73460000000000003 * exp(-T/1756);
    Fcent3 = exp(-5182 * invT);
    Fcent = Fcent1 + Fcent2 + Fcent3;
    logFcent = log10(Fcent);
    troe_c = -.4 - .67 * logFcent;
    troe_n = .75 - 1.27 * logFcent;
    troePr_den = 1.0 / (troe_n - .14*(troe_c + logPr));
    troePr = (troe_c + logPr) * troePr_den;
    troe = 1.0 / (1.0 + troePr*troePr);
    F = pow(10.0, logFcent * troe);
    dlogFcentdT = log10e/Fcent*( 
        -Fcent1/94
        -Fcent2/1756
        + Fcent3*5182*invT2);
    dlogFdcn_fac = 2.0 * logFcent * troe*troe * troePr * troePr_den;
    dlogFdc = -troe_n * dlogFdcn_fac * troePr_den;
    dlogFdn = dlogFdcn_fac * troePr;
    dlogFdlogPr = dlogFdc;
    dlogFdT = dlogFcentdT*(troe - 0.67*dlogFdc - 1.27*dlogFdn) + dlogFdlogPr * dlogPrdT;
    /* reverse */
    phi_r = sc[7];
    Kc = refCinv * exp(2.000000*g_RT[4] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[4]) + (h_RT[7]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    Corr = fPr * F;
    q = Corr * q_nocor;
    dlnCorrdT = ln10*(dlogfPrdT + dlogFdT);
    dqdT = Corr *(dlnkfdT*k_f*phi_f - dkrdT*phi_r) + dlnCorrdT*q;
    /* update wdot */
    wdot[4] -= 2 * q; /* OH */
    wdot[7] += q; /* H2O2 */
    /* for convenience */
    k_f *= Corr;
    k_r *= Corr;
    dcdc_fac = q/alpha*(1.0/(Pr+1.0) + dlogFdlogPr);
    if (consP) {
        /* d()/d[H2] */
        dqdci = (2 - 1)*dcdc_fac;
        J[4] += -2 * dqdci;           /* dwdot[OH]/d[H2] */
        J[7] += dqdci;                /* dwdot[H2O2]/d[H2] */
        /* d()/d[OH] */
        dqdci =  + k_f*2.000000*sc[4];
        J[44] += -2 * dqdci;          /* dwdot[OH]/d[OH] */
        J[47] += dqdci;               /* dwdot[H2O2]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (6 - 1)*dcdc_fac;
        J[54] += -2 * dqdci;          /* dwdot[OH]/d[H2O] */
        J[57] += dqdci;               /* dwdot[H2O2]/d[H2O] */
        /* d()/d[H2O2] */
        dqdci =  - k_r;
        J[74] += -2 * dqdci;          /* dwdot[OH]/d[H2O2] */
        J[77] += dqdci;               /* dwdot[H2O2]/d[H2O2] */
    }
    else {
        dqdc[0] = 2*dcdc_fac;
        dqdc[1] = dcdc_fac;
        dqdc[2] = dcdc_fac;
        dqdc[3] = dcdc_fac;
        dqdc[4] = dcdc_fac + k_f*2.000000*sc[4];
        dqdc[5] = 6*dcdc_fac;
        dqdc[6] = dcdc_fac;
        dqdc[7] = dcdc_fac - k_r;
        dqdc[8] = dcdc_fac;
        for (int k=0; k<9; k++) {
            J[10*k+4] += -2 * dqdc[k];
            J[10*k+7] += dqdc[k];
        }
    }
    J[94] += -2 * dqdT; /* dwdot[OH]/dT */
    J[97] += dqdT; /* dwdot[H2O2]/dT */

    /*reaction 2: 2.000000 O + M <=> O2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 2.3999999999999999 - 1)*sc[0] + ( 15.4 - 1)*sc[5];
    /* forward */
    phi_f = pow(sc[2], 2.000000);
    k_f = 1.0000000000000002e-12 * 1.2e+17
                * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[3];
    Kc = refCinv * exp(2.000000*g_RT[2] - g_RT[3]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[2]) + (h_RT[3]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= 2 * q; /* O */
    wdot[3] += q; /* O2 */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (2.3999999999999999 - 1)*q_nocor;
        J[2] += -2 * dqdci;           /* dwdot[O]/d[H2] */
        J[3] += dqdci;                /* dwdot[O2]/d[H2] */
        /* d()/d[O] */
        dqdci =  + k_f*2.000000*sc[2];
        J[22] += -2 * dqdci;          /* dwdot[O]/d[O] */
        J[23] += dqdci;               /* dwdot[O2]/d[O] */
        /* d()/d[O2] */
        dqdci =  - k_r;
        J[32] += -2 * dqdci;          /* dwdot[O]/d[O2] */
        J[33] += dqdci;               /* dwdot[O2]/d[O2] */
        /* d()/d[H2O] */
        dqdci = (15.4 - 1)*q_nocor;
        J[52] += -2 * dqdci;          /* dwdot[O]/d[H2O] */
        J[53] += dqdci;               /* dwdot[O2]/d[H2O] */
    }
    else {
        dqdc[0] = 2.3999999999999999*q_nocor;
        dqdc[1] = q_nocor;
        dqdc[2] = q_nocor + k_f*2.000000*sc[2];
        dqdc[3] = q_nocor - k_r;
        dqdc[4] = q_nocor;
        dqdc[5] = 15.4*q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        for (int k=0; k<9; k++) {
            J[10*k+2] += -2 * dqdc[k];
            J[10*k+3] += dqdc[k];
        }
    }
    J[92] += -2 * dqdT; /* dwdot[O]/dT */
    J[93] += dqdT; /* dwdot[O2]/dT */

    /*reaction 3: O + H + M <=> OH + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 2 - 1)*sc[0] + ( 6 - 1)*sc[5];
    /* forward */
    phi_f = sc[1]*sc[2];
    k_f = 1.0000000000000002e-12 * 5e+17
                * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[4];
    Kc = refCinv * exp(g_RT[1] + g_RT[2] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[2]) + (h_RT[4]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (2 - 1)*q_nocor;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[2] -= dqdci;                /* dwdot[O]/d[H2] */
        J[4] += dqdci;                /* dwdot[OH]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[2];
        J[11] -= dqdci;               /* dwdot[H]/d[H] */
        J[12] -= dqdci;               /* dwdot[O]/d[H] */
        J[14] += dqdci;               /* dwdot[OH]/d[H] */
        /* d()/d[O] */
        dqdci =  + k_f*sc[1];
        J[21] -= dqdci;               /* dwdot[H]/d[O] */
        J[22] -= dqdci;               /* dwdot[O]/d[O] */
        J[24] += dqdci;               /* dwdot[OH]/d[O] */
        /* d()/d[OH] */
        dqdci =  - k_r;
        J[41] -= dqdci;               /* dwdot[H]/d[OH] */
        J[42] -= dqdci;               /* dwdot[O]/d[OH] */
        J[44] += dqdci;               /* dwdot[OH]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (6 - 1)*q_nocor;
        J[51] -= dqdci;               /* dwdot[H]/d[H2O] */
        J[52] -= dqdci;               /* dwdot[O]/d[H2O] */
        J[54] += dqdci;               /* dwdot[OH]/d[H2O] */
    }
    else {
        dqdc[0] = 2*q_nocor;
        dqdc[1] = q_nocor + k_f*sc[2];
        dqdc[2] = q_nocor + k_f*sc[1];
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor - k_r;
        dqdc[5] = 6*q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        for (int k=0; k<9; k++) {
            J[10*k+1] -= dqdc[k];
            J[10*k+2] -= dqdc[k];
            J[10*k+4] += dqdc[k];
        }
    }
    J[91] -= dqdT; /* dwdot[H]/dT */
    J[92] -= dqdT; /* dwdot[O]/dT */
    J[94] += dqdT; /* dwdot[OH]/dT */

    /*reaction 4: H + O2 + M <=> HO2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 0 - 1)*sc[3] + ( 0 - 1)*sc[5] + ( 0 - 1)*sc[8];
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = 1.0000000000000002e-12 * 2.8e+18
                * exp(-0.85999999999999999 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -0.85999999999999999 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[6];
    Kc = refCinv * exp(g_RT[1] + g_RT[3] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3]) + (h_RT[6]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H] */
        dqdci =  + k_f*sc[3];
        J[11] -= dqdci;               /* dwdot[H]/d[H] */
        J[13] -= dqdci;               /* dwdot[O2]/d[H] */
        J[16] += dqdci;               /* dwdot[HO2]/d[H] */
        /* d()/d[O2] */
        dqdci = (0 - 1)*q_nocor + k_f*sc[1];
        J[31] -= dqdci;               /* dwdot[H]/d[O2] */
        J[33] -= dqdci;               /* dwdot[O2]/d[O2] */
        J[36] += dqdci;               /* dwdot[HO2]/d[O2] */
        /* d()/d[H2O] */
        dqdci = (0 - 1)*q_nocor;
        J[51] -= dqdci;               /* dwdot[H]/d[H2O] */
        J[53] -= dqdci;               /* dwdot[O2]/d[H2O] */
        J[56] += dqdci;               /* dwdot[HO2]/d[H2O] */
        /* d()/d[HO2] */
        dqdci =  - k_r;
        J[61] -= dqdci;               /* dwdot[H]/d[HO2] */
        J[63] -= dqdci;               /* dwdot[O2]/d[HO2] */
        J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
        /* d()/d[N2] */
        dqdci = (0 - 1)*q_nocor;
        J[81] -= dqdci;               /* dwdot[H]/d[N2] */
        J[83] -= dqdci;               /* dwdot[O2]/d[N2] */
        J[86] += dqdci;               /* dwdot[HO2]/d[N2] */
    }
    else {
        dqdc[0] = q_nocor;
        dqdc[1] = q_nocor + k_f*sc[3];
        dqdc[2] = q_nocor;
        dqdc[3] =  + k_f*sc[1];
        dqdc[4] = q_nocor;
        dqdc[6] = q_nocor - k_r;
        dqdc[7] = q_nocor;
        for (int k=0; k<9; k++) {
            J[10*k+1] -= dqdc[k];
            J[10*k+3] -= dqdc[k];
            J[10*k+6] += dqdc[k];
        }
    }
    J[91] -= dqdT; /* dwdot[H]/dT */
    J[93] -= dqdT; /* dwdot[O2]/dT */
    J[96] += dqdT; /* dwdot[HO2]/dT */

    /*reaction 5: 2.000000 H + M <=> H2 + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 0 - 1)*sc[0] + ( 0 - 1)*sc[5];
    /* forward */
    phi_f = pow(sc[1], 2.000000);
    k_f = 1.0000000000000002e-12 * 1e+18
                * exp(-1 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[0];
    Kc = refCinv * exp(-g_RT[0] + 2.000000*g_RT[1]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[1]) + (h_RT[0]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= 2 * q; /* H */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (0 - 1)*q_nocor - k_r;
        J[0] += dqdci;                /* dwdot[H2]/d[H2] */
        J[1] += -2 * dqdci;           /* dwdot[H]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*2.000000*sc[1];
        J[10] += dqdci;               /* dwdot[H2]/d[H] */
        J[11] += -2 * dqdci;          /* dwdot[H]/d[H] */
        /* d()/d[H2O] */
        dqdci = (0 - 1)*q_nocor;
        J[50] += dqdci;               /* dwdot[H2]/d[H2O] */
        J[51] += -2 * dqdci;          /* dwdot[H]/d[H2O] */
    }
    else {
        dqdc[0] =  - k_r;
        dqdc[1] = q_nocor + k_f*2.000000*sc[1];
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        for (int k=0; k<9; k++) {
            J[10*k+0] += dqdc[k];
            J[10*k+1] += -2 * dqdc[k];
        }
    }
    J[90] += dqdT; /* dwdot[H2]/dT */
    J[91] += -2 * dqdT; /* dwdot[H]/dT */

    /*reaction 6: H + OH + M <=> H2O + M */
    /*a third-body and non-pressure-fall-off reaction */
    /* 3-body correction factor */
    alpha = mixture + ( 0.72999999999999998 - 1)*sc[0] + ( 3.6499999999999999 - 1)*sc[5];
    /* forward */
    phi_f = sc[1]*sc[4];
    k_f = 1.0000000000000002e-12 * 2.2e+22
                * exp(-2 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -2 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[5];
    Kc = refCinv * exp(g_RT[1] + g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[4]) + (h_RT[5]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q_nocor = k_f*phi_f - k_r*phi_r;
    q = alpha * q_nocor;
    dqdT = alpha * (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    /* for convenience */
    k_f *= alpha;
    k_r *= alpha;
    if (consP) {
        /* d()/d[H2] */
        dqdci = (0.72999999999999998 - 1)*q_nocor;
        J[1] -= dqdci;                /* dwdot[H]/d[H2] */
        J[4] -= dqdci;                /* dwdot[OH]/d[H2] */
        J[5] += dqdci;                /* dwdot[H2O]/d[H2] */
        /* d()/d[H] */
        dqdci =  + k_f*sc[4];
        J[11] -= dqdci;               /* dwdot[H]/d[H] */
        J[14] -= dqdci;               /* dwdot[OH]/d[H] */
        J[15] += dqdci;               /* dwdot[H2O]/d[H] */
        /* d()/d[OH] */
        dqdci =  + k_f*sc[1];
        J[41] -= dqdci;               /* dwdot[H]/d[OH] */
        J[44] -= dqdci;               /* dwdot[OH]/d[OH] */
        J[45] += dqdci;               /* dwdot[H2O]/d[OH] */
        /* d()/d[H2O] */
        dqdci = (3.6499999999999999 - 1)*q_nocor - k_r;
        J[51] -= dqdci;               /* dwdot[H]/d[H2O] */
        J[54] -= dqdci;               /* dwdot[OH]/d[H2O] */
        J[55] += dqdci;               /* dwdot[H2O]/d[H2O] */
    }
    else {
        dqdc[0] = 0.72999999999999998*q_nocor;
        dqdc[1] = q_nocor + k_f*sc[4];
        dqdc[2] = q_nocor;
        dqdc[3] = q_nocor;
        dqdc[4] = q_nocor + k_f*sc[1];
        dqdc[5] = 3.6499999999999999*q_nocor - k_r;
        dqdc[6] = q_nocor;
        dqdc[7] = q_nocor;
        dqdc[8] = q_nocor;
        for (int k=0; k<9; k++) {
            J[10*k+1] -= dqdc[k];
            J[10*k+4] -= dqdc[k];
            J[10*k+5] += dqdc[k];
        }
    }
    J[91] -= dqdT; /* dwdot[H]/dT */
    J[94] -= dqdT; /* dwdot[OH]/dT */
    J[95] += dqdT; /* dwdot[H2O]/dT */

    /*reaction 7: O + H2 <=> H + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[2];
    k_f = 1.0000000000000002e-06 * 50000
                * exp(2.6699999999999999 * tc[0] - 0.50321666580471969 * (6290) * invT);
    dlnkfdT = 2.6699999999999999 * invT + 0.50321666580471969 *  6290  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[4];
    Kc = exp(g_RT[0] - g_RT[1] + g_RT[2] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[2]) + (h_RT[1] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[1] += q; /* H */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[2];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += dqdci;                /* dwdot[H]/d[H2] */
    J[2] -= dqdci;                /* dwdot[O]/d[H2] */
    J[4] += dqdci;                /* dwdot[OH]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[4];
    J[10] -= dqdci;               /* dwdot[H2]/d[H] */
    J[11] += dqdci;               /* dwdot[H]/d[H] */
    J[12] -= dqdci;               /* dwdot[O]/d[H] */
    J[14] += dqdci;               /* dwdot[OH]/d[H] */
    /* d()/d[O] */
    dqdci =  + k_f*sc[0];
    J[20] -= dqdci;               /* dwdot[H2]/d[O] */
    J[21] += dqdci;               /* dwdot[H]/d[O] */
    J[22] -= dqdci;               /* dwdot[O]/d[O] */
    J[24] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[1];
    J[40] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[41] += dqdci;               /* dwdot[H]/d[OH] */
    J[42] -= dqdci;               /* dwdot[O]/d[OH] */
    J[44] += dqdci;               /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[90] -= dqdT;                /* dwdot[H2]/dT */
    J[91] += dqdT;                /* dwdot[H]/dT */
    J[92] -= dqdT;                /* dwdot[O]/dT */
    J[94] += dqdT;                /* dwdot[OH]/dT */

    /*reaction 8: O + HO2 <=> OH + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[6];
    k_f = 1.0000000000000002e-06 * 20000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[3]*sc[4];
    Kc = exp(g_RT[2] - g_RT[3] - g_RT[4] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[6]) + (h_RT[3] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[3] += q; /* O2 */
    wdot[4] += q; /* OH */
    wdot[6] -= q; /* HO2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[6];
    J[22] -= dqdci;               /* dwdot[O]/d[O] */
    J[23] += dqdci;               /* dwdot[O2]/d[O] */
    J[24] += dqdci;               /* dwdot[OH]/d[O] */
    J[26] -= dqdci;               /* dwdot[HO2]/d[O] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[4];
    J[32] -= dqdci;               /* dwdot[O]/d[O2] */
    J[33] += dqdci;               /* dwdot[O2]/d[O2] */
    J[34] += dqdci;               /* dwdot[OH]/d[O2] */
    J[36] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[3];
    J[42] -= dqdci;               /* dwdot[O]/d[OH] */
    J[43] += dqdci;               /* dwdot[O2]/d[OH] */
    J[44] += dqdci;               /* dwdot[OH]/d[OH] */
    J[46] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[2];
    J[62] -= dqdci;               /* dwdot[O]/d[HO2] */
    J[63] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[64] += dqdci;               /* dwdot[OH]/d[HO2] */
    J[66] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[92] -= dqdT;                /* dwdot[O]/dT */
    J[93] += dqdT;                /* dwdot[O2]/dT */
    J[94] += dqdT;                /* dwdot[OH]/dT */
    J[96] -= dqdT;                /* dwdot[HO2]/dT */

    /*reaction 9: O + H2O2 <=> OH + HO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[2]*sc[7];
    k_f = 1.0000000000000002e-06 * 9630000
                * exp(2 * tc[0] - 0.50321666580471969 * (4000) * invT);
    dlnkfdT = 2 * invT + 0.50321666580471969 *  4000  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[6];
    Kc = exp(g_RT[2] - g_RT[4] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[2] + h_RT[7]) + (h_RT[4] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] -= q; /* O */
    wdot[4] += q; /* OH */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[O] */
    dqdci =  + k_f*sc[7];
    J[22] -= dqdci;               /* dwdot[O]/d[O] */
    J[24] += dqdci;               /* dwdot[OH]/d[O] */
    J[26] += dqdci;               /* dwdot[HO2]/d[O] */
    J[27] -= dqdci;               /* dwdot[H2O2]/d[O] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[6];
    J[42] -= dqdci;               /* dwdot[O]/d[OH] */
    J[44] += dqdci;               /* dwdot[OH]/d[OH] */
    J[46] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[47] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[4];
    J[62] -= dqdci;               /* dwdot[O]/d[HO2] */
    J[64] += dqdci;               /* dwdot[OH]/d[HO2] */
    J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[67] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[2];
    J[72] -= dqdci;               /* dwdot[O]/d[H2O2] */
    J[74] += dqdci;               /* dwdot[OH]/d[H2O2] */
    J[76] += dqdci;               /* dwdot[HO2]/d[H2O2] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[92] -= dqdT;                /* dwdot[O]/dT */
    J[94] += dqdT;                /* dwdot[OH]/dT */
    J[96] += dqdT;                /* dwdot[HO2]/dT */
    J[97] -= dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 10: H + 2.000000 O2 <=> HO2 + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*pow(sc[3], 2.000000);
    k_f = 1.0000000000000002e-12 * 3e+20
                * exp(-1.72 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1.72 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[3]*sc[6];
    Kc = refCinv * exp(g_RT[1] + 2.000000*g_RT[3] - g_RT[3] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + 2.000000*h_RT[3]) + (h_RT[3] + h_RT[6]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*pow(sc[3], 2.000000);
    J[11] -= dqdci;               /* dwdot[H]/d[H] */
    J[13] -= dqdci;               /* dwdot[O2]/d[H] */
    J[16] += dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1]*2.000000*sc[3] - k_r*sc[6];
    J[31] -= dqdci;               /* dwdot[H]/d[O2] */
    J[33] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[36] += dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[3];
    J[61] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[63] -= dqdci;               /* dwdot[O2]/d[HO2] */
    J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[91] -= dqdT;                /* dwdot[H]/dT */
    J[93] -= dqdT;                /* dwdot[O2]/dT */
    J[96] += dqdT;                /* dwdot[HO2]/dT */

    /*reaction 11: H + O2 + H2O <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3]*sc[5];
    k_f = 1.0000000000000002e-12 * 9.38e+18
                * exp(-0.76000000000000001 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -0.76000000000000001 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[6];
    Kc = refCinv * exp(g_RT[1] + g_RT[3] + g_RT[5] - g_RT[5] - g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3] + h_RT[5]) + (h_RT[5] + h_RT[6]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[3]*sc[5];
    J[11] -= dqdci;               /* dwdot[H]/d[H] */
    J[13] -= dqdci;               /* dwdot[O2]/d[H] */
    J[16] += dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1]*sc[5];
    J[31] -= dqdci;               /* dwdot[H]/d[O2] */
    J[33] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[36] += dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[H2O] */
    dqdci =  + k_f*sc[1]*sc[3] - k_r*sc[6];
    J[51] -= dqdci;               /* dwdot[H]/d[H2O] */
    J[53] -= dqdci;               /* dwdot[O2]/d[H2O] */
    J[56] += dqdci;               /* dwdot[HO2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[5];
    J[61] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[63] -= dqdci;               /* dwdot[O2]/d[HO2] */
    J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[91] -= dqdT;                /* dwdot[H]/dT */
    J[93] -= dqdT;                /* dwdot[O2]/dT */
    J[96] += dqdT;                /* dwdot[HO2]/dT */

    /*reaction 12: H + O2 + N2 <=> HO2 + N2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3]*sc[8];
    k_f = 1.0000000000000002e-12 * 3.75e+20
                * exp(-1.72 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1.72 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[6]*sc[8];
    Kc = refCinv * exp(g_RT[1] + g_RT[3] - g_RT[6] + g_RT[8] - g_RT[8]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3] + h_RT[8]) + (h_RT[6] + h_RT[8]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[3] -= q; /* O2 */
    wdot[6] += q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[3]*sc[8];
    J[11] -= dqdci;               /* dwdot[H]/d[H] */
    J[13] -= dqdci;               /* dwdot[O2]/d[H] */
    J[16] += dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1]*sc[8];
    J[31] -= dqdci;               /* dwdot[H]/d[O2] */
    J[33] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[36] += dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[8];
    J[61] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[63] -= dqdci;               /* dwdot[O2]/d[HO2] */
    J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/d[N2] */
    dqdci =  + k_f*sc[1]*sc[3] - k_r*sc[6];
    J[81] -= dqdci;               /* dwdot[H]/d[N2] */
    J[83] -= dqdci;               /* dwdot[O2]/d[N2] */
    J[86] += dqdci;               /* dwdot[HO2]/d[N2] */
    /* d()/dT */
    J[91] -= dqdT;                /* dwdot[H]/dT */
    J[93] -= dqdT;                /* dwdot[O2]/dT */
    J[96] += dqdT;                /* dwdot[HO2]/dT */

    /*reaction 13: H + O2 <=> O + OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[3];
    k_f = 1.0000000000000002e-06 * 83000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (14413) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  14413  * invT2;
    /* reverse */
    phi_r = sc[2]*sc[4];
    Kc = exp(g_RT[1] - g_RT[2] + g_RT[3] - g_RT[4]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[3]) + (h_RT[2] + h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[2] += q; /* O */
    wdot[3] -= q; /* O2 */
    wdot[4] += q; /* OH */
    /* d()/d[H] */
    dqdci =  + k_f*sc[3];
    J[11] -= dqdci;               /* dwdot[H]/d[H] */
    J[12] += dqdci;               /* dwdot[O]/d[H] */
    J[13] -= dqdci;               /* dwdot[O2]/d[H] */
    J[14] += dqdci;               /* dwdot[OH]/d[H] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[4];
    J[21] -= dqdci;               /* dwdot[H]/d[O] */
    J[22] += dqdci;               /* dwdot[O]/d[O] */
    J[23] -= dqdci;               /* dwdot[O2]/d[O] */
    J[24] += dqdci;               /* dwdot[OH]/d[O] */
    /* d()/d[O2] */
    dqdci =  + k_f*sc[1];
    J[31] -= dqdci;               /* dwdot[H]/d[O2] */
    J[32] += dqdci;               /* dwdot[O]/d[O2] */
    J[33] -= dqdci;               /* dwdot[O2]/d[O2] */
    J[34] += dqdci;               /* dwdot[OH]/d[O2] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[2];
    J[41] -= dqdci;               /* dwdot[H]/d[OH] */
    J[42] += dqdci;               /* dwdot[O]/d[OH] */
    J[43] -= dqdci;               /* dwdot[O2]/d[OH] */
    J[44] += dqdci;               /* dwdot[OH]/d[OH] */
    /* d()/dT */
    J[91] -= dqdT;                /* dwdot[H]/dT */
    J[92] += dqdT;                /* dwdot[O]/dT */
    J[93] -= dqdT;                /* dwdot[O2]/dT */
    J[94] += dqdT;                /* dwdot[OH]/dT */

    /*reaction 14: 2.000000 H + H2 <=> 2.000000 H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*pow(sc[1], 2.000000);
    k_f = 1.0000000000000002e-12 * 90000000000000000
                * exp(-0.59999999999999998 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -0.59999999999999998 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = pow(sc[0], 2.000000);
    Kc = refCinv * exp(g_RT[0] - 2.000000*g_RT[0] + 2.000000*g_RT[1]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + 2.000000*h_RT[1]) + (2.000000*h_RT[0]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= 2 * q; /* H */
    /* d()/d[H2] */
    dqdci =  + k_f*pow(sc[1], 2.000000) - k_r*2.000000*sc[0];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += -2 * dqdci;           /* dwdot[H]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[0]*2.000000*sc[1];
    J[10] += dqdci;               /* dwdot[H2]/d[H] */
    J[11] += -2 * dqdci;          /* dwdot[H]/d[H] */
    /* d()/dT */
    J[90] += dqdT;                /* dwdot[H2]/dT */
    J[91] += -2 * dqdT;           /* dwdot[H]/dT */

    /*reaction 15: 2.000000 H + H2O <=> H2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[1], 2.000000)*sc[5];
    k_f = 1.0000000000000002e-12 * 6e+19
                * exp(-1.25 * tc[0] - 0.50321666580471969 * (0) * invT);
    dlnkfdT = -1.25 * invT + 0.50321666580471969 *  0  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[5];
    Kc = refCinv * exp(-g_RT[0] + 2.000000*g_RT[1] + g_RT[5] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[1] + h_RT[5]) + (h_RT[0] + h_RT[5]) + 1.000000);
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= 2 * q; /* H */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[5];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += -2 * dqdci;           /* dwdot[H]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*2.000000*sc[1]*sc[5];
    J[10] += dqdci;               /* dwdot[H2]/d[H] */
    J[11] += -2 * dqdci;          /* dwdot[H]/d[H] */
    /* d()/d[H2O] */
    dqdci =  + k_f*pow(sc[1], 2.000000) - k_r*sc[0];
    J[50] += dqdci;               /* dwdot[H2]/d[H2O] */
    J[51] += -2 * dqdci;          /* dwdot[H]/d[H2O] */
    /* d()/dT */
    J[90] += dqdT;                /* dwdot[H2]/dT */
    J[91] += -2 * dqdT;           /* dwdot[H]/dT */

    /*reaction 16: H + HO2 <=> O + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[6];
    k_f = 1.0000000000000002e-06 * 3970000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (671) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  671  * invT2;
    /* reverse */
    phi_r = sc[2]*sc[5];
    Kc = exp(g_RT[1] - g_RT[2] - g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[6]) + (h_RT[2] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[2] += q; /* O */
    wdot[5] += q; /* H2O */
    wdot[6] -= q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[11] -= dqdci;               /* dwdot[H]/d[H] */
    J[12] += dqdci;               /* dwdot[O]/d[H] */
    J[15] += dqdci;               /* dwdot[H2O]/d[H] */
    J[16] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O] */
    dqdci =  - k_r*sc[5];
    J[21] -= dqdci;               /* dwdot[H]/d[O] */
    J[22] += dqdci;               /* dwdot[O]/d[O] */
    J[25] += dqdci;               /* dwdot[H2O]/d[O] */
    J[26] -= dqdci;               /* dwdot[HO2]/d[O] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[2];
    J[51] -= dqdci;               /* dwdot[H]/d[H2O] */
    J[52] += dqdci;               /* dwdot[O]/d[H2O] */
    J[55] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[56] -= dqdci;               /* dwdot[HO2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[1];
    J[61] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[62] += dqdci;               /* dwdot[O]/d[HO2] */
    J[65] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[66] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[91] -= dqdT;                /* dwdot[H]/dT */
    J[92] += dqdT;                /* dwdot[O]/dT */
    J[95] += dqdT;                /* dwdot[H2O]/dT */
    J[96] -= dqdT;                /* dwdot[HO2]/dT */

    /*reaction 17: H + HO2 <=> O2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[6];
    k_f = 1.0000000000000002e-06 * 28000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (1068) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  1068  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[3];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[3] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[6]) + (h_RT[0] + h_RT[3]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[3] += q; /* O2 */
    wdot[6] -= q; /* HO2 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[3];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[3] += dqdci;                /* dwdot[O2]/d[H2] */
    J[6] -= dqdci;                /* dwdot[HO2]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[10] += dqdci;               /* dwdot[H2]/d[H] */
    J[11] -= dqdci;               /* dwdot[H]/d[H] */
    J[13] += dqdci;               /* dwdot[O2]/d[H] */
    J[16] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[0];
    J[30] += dqdci;               /* dwdot[H2]/d[O2] */
    J[31] -= dqdci;               /* dwdot[H]/d[O2] */
    J[33] += dqdci;               /* dwdot[O2]/d[O2] */
    J[36] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[1];
    J[60] += dqdci;               /* dwdot[H2]/d[HO2] */
    J[61] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[63] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[66] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[90] += dqdT;                /* dwdot[H2]/dT */
    J[91] -= dqdT;                /* dwdot[H]/dT */
    J[93] += dqdT;                /* dwdot[O2]/dT */
    J[96] -= dqdT;                /* dwdot[HO2]/dT */

    /*reaction 18: H + HO2 <=> 2.000000 OH */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[6];
    k_f = 1.0000000000000002e-06 * 134000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (635) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  635  * invT2;
    /* reverse */
    phi_r = pow(sc[4], 2.000000);
    Kc = exp(g_RT[1] - 2.000000*g_RT[4] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[6]) + (2.000000*h_RT[4]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[4] += 2 * q; /* OH */
    wdot[6] -= q; /* HO2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[6];
    J[11] -= dqdci;               /* dwdot[H]/d[H] */
    J[14] += 2 * dqdci;           /* dwdot[OH]/d[H] */
    J[16] -= dqdci;               /* dwdot[HO2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*2.000000*sc[4];
    J[41] -= dqdci;               /* dwdot[H]/d[OH] */
    J[44] += 2 * dqdci;           /* dwdot[OH]/d[OH] */
    J[46] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[1];
    J[61] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[64] += 2 * dqdci;           /* dwdot[OH]/d[HO2] */
    J[66] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[91] -= dqdT;                /* dwdot[H]/dT */
    J[94] += 2 * dqdT;            /* dwdot[OH]/dT */
    J[96] -= dqdT;                /* dwdot[HO2]/dT */

    /*reaction 19: H + H2O2 <=> HO2 + H2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[7];
    k_f = 1.0000000000000002e-06 * 12100000
                * exp(2 * tc[0] - 0.50321666580471969 * (5200) * invT);
    dlnkfdT = 2 * invT + 0.50321666580471969 *  5200  * invT2;
    /* reverse */
    phi_r = sc[0]*sc[6];
    Kc = exp(-g_RT[0] + g_RT[1] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[7]) + (h_RT[0] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] += q; /* H2 */
    wdot[1] -= q; /* H */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H2] */
    dqdci =  - k_r*sc[6];
    J[0] += dqdci;                /* dwdot[H2]/d[H2] */
    J[1] -= dqdci;                /* dwdot[H]/d[H2] */
    J[6] += dqdci;                /* dwdot[HO2]/d[H2] */
    J[7] -= dqdci;                /* dwdot[H2O2]/d[H2] */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[10] += dqdci;               /* dwdot[H2]/d[H] */
    J[11] -= dqdci;               /* dwdot[H]/d[H] */
    J[16] += dqdci;               /* dwdot[HO2]/d[H] */
    J[17] -= dqdci;               /* dwdot[H2O2]/d[H] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[0];
    J[60] += dqdci;               /* dwdot[H2]/d[HO2] */
    J[61] -= dqdci;               /* dwdot[H]/d[HO2] */
    J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[67] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[1];
    J[70] += dqdci;               /* dwdot[H2]/d[H2O2] */
    J[71] -= dqdci;               /* dwdot[H]/d[H2O2] */
    J[76] += dqdci;               /* dwdot[HO2]/d[H2O2] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[90] += dqdT;                /* dwdot[H2]/dT */
    J[91] -= dqdT;                /* dwdot[H]/dT */
    J[96] += dqdT;                /* dwdot[HO2]/dT */
    J[97] -= dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 20: H + H2O2 <=> OH + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[1]*sc[7];
    k_f = 1.0000000000000002e-06 * 10000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (3600) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  3600  * invT2;
    /* reverse */
    phi_r = sc[4]*sc[5];
    Kc = exp(g_RT[1] - g_RT[4] - g_RT[5] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[1] + h_RT[7]) + (h_RT[4] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[1] -= q; /* H */
    wdot[4] += q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[H] */
    dqdci =  + k_f*sc[7];
    J[11] -= dqdci;               /* dwdot[H]/d[H] */
    J[14] += dqdci;               /* dwdot[OH]/d[H] */
    J[15] += dqdci;               /* dwdot[H2O]/d[H] */
    J[17] -= dqdci;               /* dwdot[H2O2]/d[H] */
    /* d()/d[OH] */
    dqdci =  - k_r*sc[5];
    J[41] -= dqdci;               /* dwdot[H]/d[OH] */
    J[44] += dqdci;               /* dwdot[OH]/d[OH] */
    J[45] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[47] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[4];
    J[51] -= dqdci;               /* dwdot[H]/d[H2O] */
    J[54] += dqdci;               /* dwdot[OH]/d[H2O] */
    J[55] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[57] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[1];
    J[71] -= dqdci;               /* dwdot[H]/d[H2O2] */
    J[74] += dqdci;               /* dwdot[OH]/d[H2O2] */
    J[75] += dqdci;               /* dwdot[H2O]/d[H2O2] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[91] -= dqdT;                /* dwdot[H]/dT */
    J[94] += dqdT;                /* dwdot[OH]/dT */
    J[95] += dqdT;                /* dwdot[H2O]/dT */
    J[97] -= dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 21: OH + H2 <=> H + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*sc[4];
    k_f = 1.0000000000000002e-06 * 216000000
                * exp(1.51 * tc[0] - 0.50321666580471969 * (3430) * invT);
    dlnkfdT = 1.51 * invT + 0.50321666580471969 *  3430  * invT2;
    /* reverse */
    phi_r = sc[1]*sc[5];
    Kc = exp(g_RT[0] - g_RT[1] + g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[0] + h_RT[4]) + (h_RT[1] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[0] -= q; /* H2 */
    wdot[1] += q; /* H */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    /* d()/d[H2] */
    dqdci =  + k_f*sc[4];
    J[0] -= dqdci;                /* dwdot[H2]/d[H2] */
    J[1] += dqdci;                /* dwdot[H]/d[H2] */
    J[4] -= dqdci;                /* dwdot[OH]/d[H2] */
    J[5] += dqdci;                /* dwdot[H2O]/d[H2] */
    /* d()/d[H] */
    dqdci =  - k_r*sc[5];
    J[10] -= dqdci;               /* dwdot[H2]/d[H] */
    J[11] += dqdci;               /* dwdot[H]/d[H] */
    J[14] -= dqdci;               /* dwdot[OH]/d[H] */
    J[15] += dqdci;               /* dwdot[H2O]/d[H] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[0];
    J[40] -= dqdci;               /* dwdot[H2]/d[OH] */
    J[41] += dqdci;               /* dwdot[H]/d[OH] */
    J[44] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[45] += dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[1];
    J[50] -= dqdci;               /* dwdot[H2]/d[H2O] */
    J[51] += dqdci;               /* dwdot[H]/d[H2O] */
    J[54] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[55] += dqdci;               /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[90] -= dqdT;                /* dwdot[H2]/dT */
    J[91] += dqdT;                /* dwdot[H]/dT */
    J[94] -= dqdT;                /* dwdot[OH]/dT */
    J[95] += dqdT;                /* dwdot[H2O]/dT */

    /*reaction 22: 2.000000 OH <=> O + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[4], 2.000000);
    k_f = 1.0000000000000002e-06 * 35700
                * exp(2.3999999999999999 * tc[0] - 0.50321666580471969 * (-2110) * invT);
    dlnkfdT = 2.3999999999999999 * invT + 0.50321666580471969 *  -2110  * invT2;
    /* reverse */
    phi_r = sc[2]*sc[5];
    Kc = exp(-g_RT[2] + 2.000000*g_RT[4] - g_RT[5]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[4]) + (h_RT[2] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[2] += q; /* O */
    wdot[4] -= 2 * q; /* OH */
    wdot[5] += q; /* H2O */
    /* d()/d[O] */
    dqdci =  - k_r*sc[5];
    J[22] += dqdci;               /* dwdot[O]/d[O] */
    J[24] += -2 * dqdci;          /* dwdot[OH]/d[O] */
    J[25] += dqdci;               /* dwdot[H2O]/d[O] */
    /* d()/d[OH] */
    dqdci =  + k_f*2.000000*sc[4];
    J[42] += dqdci;               /* dwdot[O]/d[OH] */
    J[44] += -2 * dqdci;          /* dwdot[OH]/d[OH] */
    J[45] += dqdci;               /* dwdot[H2O]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[2];
    J[52] += dqdci;               /* dwdot[O]/d[H2O] */
    J[54] += -2 * dqdci;          /* dwdot[OH]/d[H2O] */
    J[55] += dqdci;               /* dwdot[H2O]/d[H2O] */
    /* d()/dT */
    J[92] += dqdT;                /* dwdot[O]/dT */
    J[94] += -2 * dqdT;           /* dwdot[OH]/dT */
    J[95] += dqdT;                /* dwdot[H2O]/dT */

    /*reaction 23: OH + HO2 <=> O2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[6];
    k_f = 1.0000000000000002e-06 * 29000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-500) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  -500  * invT2;
    /* reverse */
    phi_r = sc[3]*sc[5];
    Kc = exp(-g_RT[3] + g_RT[4] - g_RT[5] + g_RT[6]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[6]) + (h_RT[3] + h_RT[5]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += q; /* O2 */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[6] -= q; /* HO2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[5];
    J[33] += dqdci;               /* dwdot[O2]/d[O2] */
    J[34] -= dqdci;               /* dwdot[OH]/d[O2] */
    J[35] += dqdci;               /* dwdot[H2O]/d[O2] */
    J[36] -= dqdci;               /* dwdot[HO2]/d[O2] */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[6];
    J[43] += dqdci;               /* dwdot[O2]/d[OH] */
    J[44] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[45] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[46] -= dqdci;               /* dwdot[HO2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[3];
    J[53] += dqdci;               /* dwdot[O2]/d[H2O] */
    J[54] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[55] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[56] -= dqdci;               /* dwdot[HO2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  + k_f*sc[4];
    J[63] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[64] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[65] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[66] -= dqdci;               /* dwdot[HO2]/d[HO2] */
    /* d()/dT */
    J[93] += dqdT;                /* dwdot[O2]/dT */
    J[94] -= dqdT;                /* dwdot[OH]/dT */
    J[95] += dqdT;                /* dwdot[H2O]/dT */
    J[96] -= dqdT;                /* dwdot[HO2]/dT */

    /*reaction 24: OH + H2O2 <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[7];
    k_f = 1.0000000000000002e-06 * 1750000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (320) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  320  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[6];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[7]) + (h_RT[5] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[7];
    J[44] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[45] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[46] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[47] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[6];
    J[54] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[55] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[56] += dqdci;               /* dwdot[HO2]/d[H2O] */
    J[57] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[5];
    J[64] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[65] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[67] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[4];
    J[74] -= dqdci;               /* dwdot[OH]/d[H2O2] */
    J[75] += dqdci;               /* dwdot[H2O]/d[H2O2] */
    J[76] += dqdci;               /* dwdot[HO2]/d[H2O2] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[94] -= dqdT;                /* dwdot[OH]/dT */
    J[95] += dqdT;                /* dwdot[H2O]/dT */
    J[96] += dqdT;                /* dwdot[HO2]/dT */
    J[97] -= dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 25: OH + H2O2 <=> HO2 + H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[4]*sc[7];
    k_f = 1.0000000000000002e-06 * 580000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (9560) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  9560  * invT2;
    /* reverse */
    phi_r = sc[5]*sc[6];
    Kc = exp(g_RT[4] - g_RT[5] - g_RT[6] + g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(h_RT[4] + h_RT[7]) + (h_RT[5] + h_RT[6]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[4] -= q; /* OH */
    wdot[5] += q; /* H2O */
    wdot[6] += q; /* HO2 */
    wdot[7] -= q; /* H2O2 */
    /* d()/d[OH] */
    dqdci =  + k_f*sc[7];
    J[44] -= dqdci;               /* dwdot[OH]/d[OH] */
    J[45] += dqdci;               /* dwdot[H2O]/d[OH] */
    J[46] += dqdci;               /* dwdot[HO2]/d[OH] */
    J[47] -= dqdci;               /* dwdot[H2O2]/d[OH] */
    /* d()/d[H2O] */
    dqdci =  - k_r*sc[6];
    J[54] -= dqdci;               /* dwdot[OH]/d[H2O] */
    J[55] += dqdci;               /* dwdot[H2O]/d[H2O] */
    J[56] += dqdci;               /* dwdot[HO2]/d[H2O] */
    J[57] -= dqdci;               /* dwdot[H2O2]/d[H2O] */
    /* d()/d[HO2] */
    dqdci =  - k_r*sc[5];
    J[64] -= dqdci;               /* dwdot[OH]/d[HO2] */
    J[65] += dqdci;               /* dwdot[H2O]/d[HO2] */
    J[66] += dqdci;               /* dwdot[HO2]/d[HO2] */
    J[67] -= dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  + k_f*sc[4];
    J[74] -= dqdci;               /* dwdot[OH]/d[H2O2] */
    J[75] += dqdci;               /* dwdot[H2O]/d[H2O2] */
    J[76] += dqdci;               /* dwdot[HO2]/d[H2O2] */
    J[77] -= dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[94] -= dqdT;                /* dwdot[OH]/dT */
    J[95] += dqdT;                /* dwdot[H2O]/dT */
    J[96] += dqdT;                /* dwdot[HO2]/dT */
    J[97] -= dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 26: 2.000000 HO2 <=> O2 + H2O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[6], 2.000000);
    k_f = 1.0000000000000002e-06 * 130000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (-1630) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  -1630  * invT2;
    /* reverse */
    phi_r = sc[3]*sc[7];
    Kc = exp(-g_RT[3] + 2.000000*g_RT[6] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[6]) + (h_RT[3] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += q; /* O2 */
    wdot[6] -= 2 * q; /* HO2 */
    wdot[7] += q; /* H2O2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[7];
    J[33] += dqdci;               /* dwdot[O2]/d[O2] */
    J[36] += -2 * dqdci;          /* dwdot[HO2]/d[O2] */
    J[37] += dqdci;               /* dwdot[H2O2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*2.000000*sc[6];
    J[63] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[66] += -2 * dqdci;          /* dwdot[HO2]/d[HO2] */
    J[67] += dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[3];
    J[73] += dqdci;               /* dwdot[O2]/d[H2O2] */
    J[76] += -2 * dqdci;          /* dwdot[HO2]/d[H2O2] */
    J[77] += dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[93] += dqdT;                /* dwdot[O2]/dT */
    J[96] += -2 * dqdT;           /* dwdot[HO2]/dT */
    J[97] += dqdT;                /* dwdot[H2O2]/dT */

    /*reaction 27: 2.000000 HO2 <=> O2 + H2O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[6], 2.000000);
    k_f = 1.0000000000000002e-06 * 420000000000000
                * exp(0 * tc[0] - 0.50321666580471969 * (12000) * invT);
    dlnkfdT = 0 * invT + 0.50321666580471969 *  12000  * invT2;
    /* reverse */
    phi_r = sc[3]*sc[7];
    Kc = exp(-g_RT[3] + 2.000000*g_RT[6] - g_RT[7]);
    k_r = k_f / Kc;
    dlnKcdT = invT * (-(2.000000*h_RT[6]) + (h_RT[3] + h_RT[7]));
    dkrdT = (dlnkfdT - dlnKcdT)*k_r;
    /* rate of progress */
    q = k_f*phi_f - k_r*phi_r;
    dqdT = (dlnkfdT*k_f*phi_f - dkrdT*phi_r);
    /* update wdot */
    wdot[3] += q; /* O2 */
    wdot[6] -= 2 * q; /* HO2 */
    wdot[7] += q; /* H2O2 */
    /* d()/d[O2] */
    dqdci =  - k_r*sc[7];
    J[33] += dqdci;               /* dwdot[O2]/d[O2] */
    J[36] += -2 * dqdci;          /* dwdot[HO2]/d[O2] */
    J[37] += dqdci;               /* dwdot[H2O2]/d[O2] */
    /* d()/d[HO2] */
    dqdci =  + k_f*2.000000*sc[6];
    J[63] += dqdci;               /* dwdot[O2]/d[HO2] */
    J[66] += -2 * dqdci;          /* dwdot[HO2]/d[HO2] */
    J[67] += dqdci;               /* dwdot[H2O2]/d[HO2] */
    /* d()/d[H2O2] */
    dqdci =  - k_r*sc[3];
    J[73] += dqdci;               /* dwdot[O2]/d[H2O2] */
    J[76] += -2 * dqdci;          /* dwdot[HO2]/d[H2O2] */
    J[77] += dqdci;               /* dwdot[H2O2]/d[H2O2] */
    /* d()/dT */
    J[93] += dqdT;                /* dwdot[O2]/dT */
    J[96] += -2 * dqdT;           /* dwdot[HO2]/dT */
    J[97] += dqdT;                /* dwdot[H2O2]/dT */

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
    /* dTdot/d[X] */
    for (int k = 0; k < 9; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 9; ++m) {
            dehmixdc += eh_RT[m]*J[k*10+m];
        }
        J[k*10+9] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[99] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;

return;
}


/*compute the reaction Jacobian */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void DWDOT(amrex::Real *  J, amrex::Real *  sc, amrex::Real *  Tp, int * consP)
{
    amrex::Real c[9];

    for (int k=0; k<9; k++) {
        c[k] = 1.e6 * sc[k];
    }

    aJacobian(J, c, *Tp, *consP);

    /* dwdot[k]/dT */
    /* dTdot/d[X] */
    for (int k=0; k<9; k++) {
        J[90+k] *= 1.e-6;
        J[k*10+9] *= 1.e6;
    }

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
