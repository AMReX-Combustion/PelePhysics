#ifndef MECHANISM_H
#define MECHANISM_H

#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>

#if 0
/* Elements
0  O
1  H
2  C
3  N
*/
#endif

/* Species */
#define O2_ID 0
#define H2O_ID 1
#define CH4_ID 2
#define CO_ID 3
#define CO2_ID 4
#define N2_ID 5

#define NUM_ELEMENTS 4
#define NUM_SPECIES 6
#define NUM_REACTIONS 3

#define NUM_FIT 4


/* ALWAYS on CPU stuff -- can have different def depending on if we are CPU or GPU based. Defined in mechanism.cpp  */
void atomicWeight(amrex::Real *  awt);
/* MISC  */
void CKAWT(amrex::Real *  awt);
void CKNCF(int * ncf);
void CKSYME_STR(amrex::Vector<std::string>& ename);
void CKSYMS_STR(amrex::Vector<std::string>& kname);
/* SPARSE INFORMATION  */
void SPARSITY_INFO(int * nJdata, int * consP, int NCELLS);
void SPARSITY_INFO_SYST(int * nJdata, int * consP, int NCELLS);
void SPARSITY_INFO_SYST_SIMPLIFIED(int * nJdata, int * consP);
void SPARSITY_PREPROC_CSC(int * rowVals, int * colPtrs, int * consP, int NCELLS);
void SPARSITY_PREPROC_CSR(int * colVals, int * rowPtrs, int * consP, int NCELLS, int base);
void SPARSITY_PREPROC_SYST_CSR(int * colVals, int * rowPtrs, int * consP, int NCELLS, int base);
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(int * rowVals, int * colPtrs, int * indx, int * consP);
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(int * colVals, int * rowPtr, int * consP, int base);



/*A few mechanism parameters */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKINDX(int * mm, int * kk, int * ii, int * nfit)
{
    *mm = 4;
    *kk = 6;
    *ii = 3;
    *nfit = -1; /*Why do you need this anyway ?  */
}


/* inverse molecular weights  */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void get_imw(amrex::Real *imw_new){
    imw_new[0] = 1.0/31.998800;/*O2 */
    imw_new[1] = 1.0/18.015340;/*H2O */
    imw_new[2] = 1.0/16.043030;/*CH4 */
    imw_new[3] = 1.0/28.010550;/*CO */
    imw_new[4] = 1.0/44.009950;/*CO2 */
    imw_new[5] = 1.0/28.013400;/*N2 */
}

/* molecular weights  */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void get_mw(amrex::Real *mw_new){
    mw_new[0] = 31.998800;/*O2 */
    mw_new[1] = 18.015340;/*H2O */
    mw_new[2] = 16.043030;/*CH4 */
    mw_new[3] = 28.010550;/*CO */
    mw_new[4] = 44.009950;/*CO2 */
    mw_new[5] = 28.013400;/*N2 */
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
        /*species 0: O2 */
        species[0] =
            +2.78245636e+00
            -2.99673416e-03 * tc[1]
            +9.84730201e-06 * tc[2]
            -9.68129509e-09 * tc[3]
            +3.24372837e-12 * tc[4];
        /*species 1: H2O */
        species[1] =
            +3.19864056e+00
            -2.03643410e-03 * tc[1]
            +6.52040211e-06 * tc[2]
            -5.48797062e-09 * tc[3]
            +1.77197817e-12 * tc[4];
        /*species 2: CH4 */
        species[2] =
            +4.14987613e+00
            -1.36709788e-02 * tc[1]
            +4.91800599e-05 * tc[2]
            -4.84743026e-08 * tc[3]
            +1.66693956e-11 * tc[4];
        /*species 3: CO */
        species[3] =
            +2.57953347e+00
            -6.10353680e-04 * tc[1]
            +1.01681433e-06 * tc[2]
            +9.07005884e-10 * tc[3]
            -9.04424499e-13 * tc[4];
        /*species 4: CO2 */
        species[4] =
            +1.35677352e+00
            +8.98459677e-03 * tc[1]
            -7.12356269e-06 * tc[2]
            +2.45919022e-09 * tc[3]
            -1.43699548e-13 * tc[4];
        /*species 5: N2 */
        species[5] =
            +2.29867700e+00
            +1.40824040e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485400e-12 * tc[4];
    } else {
        /*species 0: O2 */
        species[0] =
            +2.28253784e+00
            +1.48308754e-03 * tc[1]
            -7.57966669e-07 * tc[2]
            +2.09470555e-10 * tc[3]
            -2.16717794e-14 * tc[4];
        /*species 1: H2O */
        species[1] =
            +2.03399249e+00
            +2.17691804e-03 * tc[1]
            -1.64072518e-07 * tc[2]
            -9.70419870e-11 * tc[3]
            +1.68200992e-14 * tc[4];
        /*species 2: CH4 */
        species[2] =
            -9.25148505e-01
            +1.33909467e-02 * tc[1]
            -5.73285809e-06 * tc[2]
            +1.22292535e-09 * tc[3]
            -1.01815230e-13 * tc[4];
        /*species 3: CO */
        species[3] =
            +1.71518561e+00
            +2.06252743e-03 * tc[1]
            -9.98825771e-07 * tc[2]
            +2.30053008e-10 * tc[3]
            -2.03647716e-14 * tc[4];
        /*species 4: CO2 */
        species[4] =
            +2.85746029e+00
            +4.41437026e-03 * tc[1]
            -2.21481404e-06 * tc[2]
            +5.23490188e-10 * tc[3]
            -4.72084164e-14 * tc[4];
        /*species 5: N2 */
        species[5] =
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
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void cp_R(amrex::Real * species, const amrex::Real *  tc)
{

    /*temperature */
    const amrex::Real T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: O2 */
        species[0] =
            +3.78245636e+00
            -2.99673416e-03 * tc[1]
            +9.84730201e-06 * tc[2]
            -9.68129509e-09 * tc[3]
            +3.24372837e-12 * tc[4];
        /*species 1: H2O */
        species[1] =
            +4.19864056e+00
            -2.03643410e-03 * tc[1]
            +6.52040211e-06 * tc[2]
            -5.48797062e-09 * tc[3]
            +1.77197817e-12 * tc[4];
        /*species 2: CH4 */
        species[2] =
            +5.14987613e+00
            -1.36709788e-02 * tc[1]
            +4.91800599e-05 * tc[2]
            -4.84743026e-08 * tc[3]
            +1.66693956e-11 * tc[4];
        /*species 3: CO */
        species[3] =
            +3.57953347e+00
            -6.10353680e-04 * tc[1]
            +1.01681433e-06 * tc[2]
            +9.07005884e-10 * tc[3]
            -9.04424499e-13 * tc[4];
        /*species 4: CO2 */
        species[4] =
            +2.35677352e+00
            +8.98459677e-03 * tc[1]
            -7.12356269e-06 * tc[2]
            +2.45919022e-09 * tc[3]
            -1.43699548e-13 * tc[4];
        /*species 5: N2 */
        species[5] =
            +3.29867700e+00
            +1.40824040e-03 * tc[1]
            -3.96322200e-06 * tc[2]
            +5.64151500e-09 * tc[3]
            -2.44485400e-12 * tc[4];
    } else {
        /*species 0: O2 */
        species[0] =
            +3.28253784e+00
            +1.48308754e-03 * tc[1]
            -7.57966669e-07 * tc[2]
            +2.09470555e-10 * tc[3]
            -2.16717794e-14 * tc[4];
        /*species 1: H2O */
        species[1] =
            +3.03399249e+00
            +2.17691804e-03 * tc[1]
            -1.64072518e-07 * tc[2]
            -9.70419870e-11 * tc[3]
            +1.68200992e-14 * tc[4];
        /*species 2: CH4 */
        species[2] =
            +7.48514950e-02
            +1.33909467e-02 * tc[1]
            -5.73285809e-06 * tc[2]
            +1.22292535e-09 * tc[3]
            -1.01815230e-13 * tc[4];
        /*species 3: CO */
        species[3] =
            +2.71518561e+00
            +2.06252743e-03 * tc[1]
            -9.98825771e-07 * tc[2]
            +2.30053008e-10 * tc[3]
            -2.03647716e-14 * tc[4];
        /*species 4: CO2 */
        species[4] =
            +3.85746029e+00
            +4.41437026e-03 * tc[1]
            -2.21481404e-06 * tc[2]
            +5.23490188e-10 * tc[3]
            -4.72084164e-14 * tc[4];
        /*species 5: N2 */
        species[5] =
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
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void gibbs(amrex::Real * species, const amrex::Real *  tc)
{

    /*temperature */
    const amrex::Real T = tc[1];
    const amrex::Real invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: O2 */
        species[0] =
            -1.063943560000000e+03 * invT
            +1.247806300000001e-01
            -3.782456360000000e+00 * tc[0]
            +1.498367080000000e-03 * tc[1]
            -1.641217001666667e-06 * tc[2]
            +8.067745908333334e-10 * tc[3]
            -1.621864185000000e-13 * tc[4];
        /*species 1: H2O */
        species[1] =
            -3.029372670000000e+04 * invT
            +5.047672768000000e+00
            -4.198640560000000e+00 * tc[0]
            +1.018217050000000e-03 * tc[1]
            -1.086733685000000e-06 * tc[2]
            +4.573308850000000e-10 * tc[3]
            -8.859890850000000e-14 * tc[4];
        /*species 2: CH4 */
        species[2] =
            -1.024664760000000e+04 * invT
            +9.791179889999999e+00
            -5.149876130000000e+00 * tc[0]
            +6.835489400000000e-03 * tc[1]
            -8.196676650000000e-06 * tc[2]
            +4.039525216666667e-09 * tc[3]
            -8.334697800000000e-13 * tc[4];
        /*species 3: CO */
        species[3] =
            -1.434408600000000e+04 * invT
            +7.112418999999992e-02
            -3.579533470000000e+00 * tc[0]
            +3.051768400000000e-04 * tc[1]
            -1.694690550000000e-07 * tc[2]
            -7.558382366666667e-11 * tc[3]
            +4.522122495000000e-14 * tc[4];
        /*species 4: CO2 */
        species[4] =
            -4.837196970000000e+04 * invT
            -7.544278700000000e+00
            -2.356773520000000e+00 * tc[0]
            -4.492298385000000e-03 * tc[1]
            +1.187260448333333e-06 * tc[2]
            -2.049325183333333e-10 * tc[3]
            +7.184977399999999e-15 * tc[4];
        /*species 5: N2 */
        species[5] =
            -1.020899900000000e+03 * invT
            -6.516950000000001e-01
            -3.298677000000000e+00 * tc[0]
            -7.041202000000000e-04 * tc[1]
            +6.605369999999999e-07 * tc[2]
            -4.701262500000001e-10 * tc[3]
            +1.222427000000000e-13 * tc[4];
    } else {
        /*species 0: O2 */
        species[0] =
            -1.088457720000000e+03 * invT
            -2.170693450000000e+00
            -3.282537840000000e+00 * tc[0]
            -7.415437700000000e-04 * tc[1]
            +1.263277781666667e-07 * tc[2]
            -1.745587958333333e-11 * tc[3]
            +1.083588970000000e-15 * tc[4];
        /*species 1: H2O */
        species[1] =
            -3.000429710000000e+04 * invT
            -1.932777610000000e+00
            -3.033992490000000e+00 * tc[0]
            -1.088459020000000e-03 * tc[1]
            +2.734541966666666e-08 * tc[2]
            +8.086832250000000e-12 * tc[3]
            -8.410049600000000e-16 * tc[4];
        /*species 2: CH4 */
        species[2] =
            -9.468344590000001e+03 * invT
            -1.836246650500000e+01
            -7.485149500000000e-02 * tc[0]
            -6.695473350000000e-03 * tc[1]
            +9.554763483333333e-07 * tc[2]
            -1.019104458333333e-10 * tc[3]
            +5.090761500000000e-15 * tc[4];
        /*species 3: CO */
        species[3] =
            -1.415187240000000e+04 * invT
            -5.103502110000000e+00
            -2.715185610000000e+00 * tc[0]
            -1.031263715000000e-03 * tc[1]
            +1.664709618333334e-07 * tc[2]
            -1.917108400000000e-11 * tc[3]
            +1.018238580000000e-15 * tc[4];
        /*species 4: CO2 */
        species[4] =
            -4.875916600000000e+04 * invT
            +1.585822230000000e+00
            -3.857460290000000e+00 * tc[0]
            -2.207185130000000e-03 * tc[1]
            +3.691356733333334e-07 * tc[2]
            -4.362418233333334e-11 * tc[3]
            +2.360420820000000e-15 * tc[4];
        /*species 5: N2 */
        species[5] =
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
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void helmholtz(amrex::Real * species, const amrex::Real *  tc)
{

    /*temperature */
    const amrex::Real T = tc[1];
    const amrex::Real invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: O2 */
        species[0] =
            -1.06394356e+03 * invT
            -8.75219370e-01
            -3.78245636e+00 * tc[0]
            +1.49836708e-03 * tc[1]
            -1.64121700e-06 * tc[2]
            +8.06774591e-10 * tc[3]
            -1.62186418e-13 * tc[4];
        /*species 1: H2O */
        species[1] =
            -3.02937267e+04 * invT
            +4.04767277e+00
            -4.19864056e+00 * tc[0]
            +1.01821705e-03 * tc[1]
            -1.08673369e-06 * tc[2]
            +4.57330885e-10 * tc[3]
            -8.85989085e-14 * tc[4];
        /*species 2: CH4 */
        species[2] =
            -1.02466476e+04 * invT
            +8.79117989e+00
            -5.14987613e+00 * tc[0]
            +6.83548940e-03 * tc[1]
            -8.19667665e-06 * tc[2]
            +4.03952522e-09 * tc[3]
            -8.33469780e-13 * tc[4];
        /*species 3: CO */
        species[3] =
            -1.43440860e+04 * invT
            -9.28875810e-01
            -3.57953347e+00 * tc[0]
            +3.05176840e-04 * tc[1]
            -1.69469055e-07 * tc[2]
            -7.55838237e-11 * tc[3]
            +4.52212249e-14 * tc[4];
        /*species 4: CO2 */
        species[4] =
            -4.83719697e+04 * invT
            -8.54427870e+00
            -2.35677352e+00 * tc[0]
            -4.49229839e-03 * tc[1]
            +1.18726045e-06 * tc[2]
            -2.04932518e-10 * tc[3]
            +7.18497740e-15 * tc[4];
        /*species 5: N2 */
        species[5] =
            -1.02089990e+03 * invT
            -1.65169500e+00
            -3.29867700e+00 * tc[0]
            -7.04120200e-04 * tc[1]
            +6.60537000e-07 * tc[2]
            -4.70126250e-10 * tc[3]
            +1.22242700e-13 * tc[4];
    } else {
        /*species 0: O2 */
        species[0] =
            -1.08845772e+03 * invT
            -3.17069345e+00
            -3.28253784e+00 * tc[0]
            -7.41543770e-04 * tc[1]
            +1.26327778e-07 * tc[2]
            -1.74558796e-11 * tc[3]
            +1.08358897e-15 * tc[4];
        /*species 1: H2O */
        species[1] =
            -3.00042971e+04 * invT
            -2.93277761e+00
            -3.03399249e+00 * tc[0]
            -1.08845902e-03 * tc[1]
            +2.73454197e-08 * tc[2]
            +8.08683225e-12 * tc[3]
            -8.41004960e-16 * tc[4];
        /*species 2: CH4 */
        species[2] =
            -9.46834459e+03 * invT
            -1.93624665e+01
            -7.48514950e-02 * tc[0]
            -6.69547335e-03 * tc[1]
            +9.55476348e-07 * tc[2]
            -1.01910446e-10 * tc[3]
            +5.09076150e-15 * tc[4];
        /*species 3: CO */
        species[3] =
            -1.41518724e+04 * invT
            -6.10350211e+00
            -2.71518561e+00 * tc[0]
            -1.03126372e-03 * tc[1]
            +1.66470962e-07 * tc[2]
            -1.91710840e-11 * tc[3]
            +1.01823858e-15 * tc[4];
        /*species 4: CO2 */
        species[4] =
            -4.87591660e+04 * invT
            +5.85822230e-01
            -3.85746029e+00 * tc[0]
            -2.20718513e-03 * tc[1]
            +3.69135673e-07 * tc[2]
            -4.36241823e-11 * tc[3]
            +2.36042082e-15 * tc[4];
        /*species 5: N2 */
        species[5] =
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
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void speciesInternalEnergy(amrex::Real * species, const amrex::Real *  tc)
{

    /*temperature */
    const amrex::Real T = tc[1];
    const amrex::Real invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: O2 */
        species[0] =
            +2.78245636e+00
            -1.49836708e-03 * tc[1]
            +3.28243400e-06 * tc[2]
            -2.42032377e-09 * tc[3]
            +6.48745674e-13 * tc[4]
            -1.06394356e+03 * invT;
        /*species 1: H2O */
        species[1] =
            +3.19864056e+00
            -1.01821705e-03 * tc[1]
            +2.17346737e-06 * tc[2]
            -1.37199266e-09 * tc[3]
            +3.54395634e-13 * tc[4]
            -3.02937267e+04 * invT;
        /*species 2: CH4 */
        species[2] =
            +4.14987613e+00
            -6.83548940e-03 * tc[1]
            +1.63933533e-05 * tc[2]
            -1.21185757e-08 * tc[3]
            +3.33387912e-12 * tc[4]
            -1.02466476e+04 * invT;
        /*species 3: CO */
        species[3] =
            +2.57953347e+00
            -3.05176840e-04 * tc[1]
            +3.38938110e-07 * tc[2]
            +2.26751471e-10 * tc[3]
            -1.80884900e-13 * tc[4]
            -1.43440860e+04 * invT;
        /*species 4: CO2 */
        species[4] =
            +1.35677352e+00
            +4.49229839e-03 * tc[1]
            -2.37452090e-06 * tc[2]
            +6.14797555e-10 * tc[3]
            -2.87399096e-14 * tc[4]
            -4.83719697e+04 * invT;
        /*species 5: N2 */
        species[5] =
            +2.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 * invT;
    } else {
        /*species 0: O2 */
        species[0] =
            +2.28253784e+00
            +7.41543770e-04 * tc[1]
            -2.52655556e-07 * tc[2]
            +5.23676387e-11 * tc[3]
            -4.33435588e-15 * tc[4]
            -1.08845772e+03 * invT;
        /*species 1: H2O */
        species[1] =
            +2.03399249e+00
            +1.08845902e-03 * tc[1]
            -5.46908393e-08 * tc[2]
            -2.42604967e-11 * tc[3]
            +3.36401984e-15 * tc[4]
            -3.00042971e+04 * invT;
        /*species 2: CH4 */
        species[2] =
            -9.25148505e-01
            +6.69547335e-03 * tc[1]
            -1.91095270e-06 * tc[2]
            +3.05731338e-10 * tc[3]
            -2.03630460e-14 * tc[4]
            -9.46834459e+03 * invT;
        /*species 3: CO */
        species[3] =
            +1.71518561e+00
            +1.03126372e-03 * tc[1]
            -3.32941924e-07 * tc[2]
            +5.75132520e-11 * tc[3]
            -4.07295432e-15 * tc[4]
            -1.41518724e+04 * invT;
        /*species 4: CO2 */
        species[4] =
            +2.85746029e+00
            +2.20718513e-03 * tc[1]
            -7.38271347e-07 * tc[2]
            +1.30872547e-10 * tc[3]
            -9.44168328e-15 * tc[4]
            -4.87591660e+04 * invT;
        /*species 5: N2 */
        species[5] =
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
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void speciesEnthalpy(amrex::Real * species, const amrex::Real *  tc)
{

    /*temperature */
    const amrex::Real T = tc[1];
    const amrex::Real invT = 1 / T;

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: O2 */
        species[0] =
            +3.78245636e+00
            -1.49836708e-03 * tc[1]
            +3.28243400e-06 * tc[2]
            -2.42032377e-09 * tc[3]
            +6.48745674e-13 * tc[4]
            -1.06394356e+03 * invT;
        /*species 1: H2O */
        species[1] =
            +4.19864056e+00
            -1.01821705e-03 * tc[1]
            +2.17346737e-06 * tc[2]
            -1.37199266e-09 * tc[3]
            +3.54395634e-13 * tc[4]
            -3.02937267e+04 * invT;
        /*species 2: CH4 */
        species[2] =
            +5.14987613e+00
            -6.83548940e-03 * tc[1]
            +1.63933533e-05 * tc[2]
            -1.21185757e-08 * tc[3]
            +3.33387912e-12 * tc[4]
            -1.02466476e+04 * invT;
        /*species 3: CO */
        species[3] =
            +3.57953347e+00
            -3.05176840e-04 * tc[1]
            +3.38938110e-07 * tc[2]
            +2.26751471e-10 * tc[3]
            -1.80884900e-13 * tc[4]
            -1.43440860e+04 * invT;
        /*species 4: CO2 */
        species[4] =
            +2.35677352e+00
            +4.49229839e-03 * tc[1]
            -2.37452090e-06 * tc[2]
            +6.14797555e-10 * tc[3]
            -2.87399096e-14 * tc[4]
            -4.83719697e+04 * invT;
        /*species 5: N2 */
        species[5] =
            +3.29867700e+00
            +7.04120200e-04 * tc[1]
            -1.32107400e-06 * tc[2]
            +1.41037875e-09 * tc[3]
            -4.88970800e-13 * tc[4]
            -1.02089990e+03 * invT;
    } else {
        /*species 0: O2 */
        species[0] =
            +3.28253784e+00
            +7.41543770e-04 * tc[1]
            -2.52655556e-07 * tc[2]
            +5.23676387e-11 * tc[3]
            -4.33435588e-15 * tc[4]
            -1.08845772e+03 * invT;
        /*species 1: H2O */
        species[1] =
            +3.03399249e+00
            +1.08845902e-03 * tc[1]
            -5.46908393e-08 * tc[2]
            -2.42604967e-11 * tc[3]
            +3.36401984e-15 * tc[4]
            -3.00042971e+04 * invT;
        /*species 2: CH4 */
        species[2] =
            +7.48514950e-02
            +6.69547335e-03 * tc[1]
            -1.91095270e-06 * tc[2]
            +3.05731338e-10 * tc[3]
            -2.03630460e-14 * tc[4]
            -9.46834459e+03 * invT;
        /*species 3: CO */
        species[3] =
            +2.71518561e+00
            +1.03126372e-03 * tc[1]
            -3.32941924e-07 * tc[2]
            +5.75132520e-11 * tc[3]
            -4.07295432e-15 * tc[4]
            -1.41518724e+04 * invT;
        /*species 4: CO2 */
        species[4] =
            +3.85746029e+00
            +2.20718513e-03 * tc[1]
            -7.38271347e-07 * tc[2]
            +1.30872547e-10 * tc[3]
            -9.44168328e-15 * tc[4]
            -4.87591660e+04 * invT;
        /*species 5: N2 */
        species[5] =
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
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void speciesEntropy(amrex::Real * species, const amrex::Real *  tc)
{

    /*temperature */
    const amrex::Real T = tc[1];

    /*species with midpoint at T=1000 kelvin */
    if (T < 1000) {
        /*species 0: O2 */
        species[0] =
            +3.78245636e+00 * tc[0]
            -2.99673416e-03 * tc[1]
            +4.92365101e-06 * tc[2]
            -3.22709836e-09 * tc[3]
            +8.10932092e-13 * tc[4]
            +3.65767573e+00 ;
        /*species 1: H2O */
        species[1] =
            +4.19864056e+00 * tc[0]
            -2.03643410e-03 * tc[1]
            +3.26020105e-06 * tc[2]
            -1.82932354e-09 * tc[3]
            +4.42994543e-13 * tc[4]
            -8.49032208e-01 ;
        /*species 2: CH4 */
        species[2] =
            +5.14987613e+00 * tc[0]
            -1.36709788e-02 * tc[1]
            +2.45900299e-05 * tc[2]
            -1.61581009e-08 * tc[3]
            +4.16734890e-12 * tc[4]
            -4.64130376e+00 ;
        /*species 3: CO */
        species[3] =
            +3.57953347e+00 * tc[0]
            -6.10353680e-04 * tc[1]
            +5.08407165e-07 * tc[2]
            +3.02335295e-10 * tc[3]
            -2.26106125e-13 * tc[4]
            +3.50840928e+00 ;
        /*species 4: CO2 */
        species[4] =
            +2.35677352e+00 * tc[0]
            +8.98459677e-03 * tc[1]
            -3.56178134e-06 * tc[2]
            +8.19730073e-10 * tc[3]
            -3.59248870e-14 * tc[4]
            +9.90105222e+00 ;
        /*species 5: N2 */
        species[5] =
            +3.29867700e+00 * tc[0]
            +1.40824040e-03 * tc[1]
            -1.98161100e-06 * tc[2]
            +1.88050500e-09 * tc[3]
            -6.11213500e-13 * tc[4]
            +3.95037200e+00 ;
    } else {
        /*species 0: O2 */
        species[0] =
            +3.28253784e+00 * tc[0]
            +1.48308754e-03 * tc[1]
            -3.78983334e-07 * tc[2]
            +6.98235183e-11 * tc[3]
            -5.41794485e-15 * tc[4]
            +5.45323129e+00 ;
        /*species 1: H2O */
        species[1] =
            +3.03399249e+00 * tc[0]
            +2.17691804e-03 * tc[1]
            -8.20362590e-08 * tc[2]
            -3.23473290e-11 * tc[3]
            +4.20502480e-15 * tc[4]
            +4.96677010e+00 ;
        /*species 2: CH4 */
        species[2] =
            +7.48514950e-02 * tc[0]
            +1.33909467e-02 * tc[1]
            -2.86642905e-06 * tc[2]
            +4.07641783e-10 * tc[3]
            -2.54538075e-14 * tc[4]
            +1.84373180e+01 ;
        /*species 3: CO */
        species[3] =
            +2.71518561e+00 * tc[0]
            +2.06252743e-03 * tc[1]
            -4.99412886e-07 * tc[2]
            +7.66843360e-11 * tc[3]
            -5.09119290e-15 * tc[4]
            +7.81868772e+00 ;
        /*species 4: CO2 */
        species[4] =
            +3.85746029e+00 * tc[0]
            +4.41437026e-03 * tc[1]
            -1.10740702e-06 * tc[2]
            +1.74496729e-10 * tc[3]
            -1.18021041e-14 * tc[4]
            +2.27163806e+00 ;
        /*species 5: N2 */
        species[5] =
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
    const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real cpor[6]; /* temporary storage */
    cp_R(cpor, tc);

    /*perform dot product */
    for (id = 0; id < 6; ++id) {
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
    amrex::Real cpor[6], tresult[6]; /* temporary storage */
    amrex::Real imw[6];

    get_imw(imw);

    cp_R(cpor, tc);

    for (int i = 0; i < 6; i++)
    {
        tresult[i] = cpor[i]*y[i]*imw[i];

    }
    for (int i = 0; i < 6; i++)
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
    amrex::Real cvor[6]; /* temporary storage */
    cv_R(cvor, tc);

    /*perform dot product */
    for (id = 0; id < 6; ++id) {
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
    amrex::Real cvor[6]; /* temporary storage */
    amrex::Real imw[6];

    get_imw(imw);

    cv_R(cvor, tc);

    /*multiply by y/molecularweight */
    result += cvor[0]*y[0]*imw[0]; /*O2 */
    result += cvor[1]*y[1]*imw[1]; /*H2O */
    result += cvor[2]*y[2]*imw[2]; /*CH4 */
    result += cvor[3]*y[3]*imw[3]; /*CO */
    result += cvor[4]*y[4]*imw[4]; /*CO2 */
    result += cvor[5]*y[5]*imw[5]; /*N2 */

    *cvbs = result * 8.31446261815324e+07;
}


/*Returns the mean enthalpy of the mixture in molar units */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKHBML(amrex::Real *  T, amrex::Real *  x,  amrex::Real *  hbml)
{
    int id; /*loop counter */
    amrex::Real result = 0; 
    amrex::Real tT = *T; /*temporary temperature */
    const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real hml[6]; /* temporary storage */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    speciesEnthalpy(hml, tc);

    /*perform dot product */
    for (id = 0; id < 6; ++id) {
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
    amrex::Real hml[6], tmp[6]; /* temporary storage */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    amrex::Real imw[6];

    get_imw(imw);

    speciesEnthalpy(hml, tc);

    int id;
    for (id = 0; id < 6; ++id) {
        tmp[id] = y[id]*hml[id]*imw[id];
    }
    for (id = 0; id < 6; ++id) {
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
    amrex::Real uml[6]; /* temporary energy array */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    speciesInternalEnergy(uml, tc);

    /*perform dot product */
    for (id = 0; id < 6; ++id) {
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
    amrex::Real ums[6]; /* temporary energy array */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    amrex::Real imw[6];

    get_imw(imw);

    speciesInternalEnergy(ums, tc);

    /*perform dot product + scaling by wt */
    result += y[0]*ums[0]*imw[0]; /*O2 */
    result += y[1]*ums[1]*imw[1]; /*H2O */
    result += y[2]*ums[2]*imw[2]; /*CH4 */
    result += y[3]*ums[3]*imw[3]; /*CO */
    result += y[4]*ums[4]*imw[4]; /*CO2 */
    result += y[5]*ums[5]*imw[5]; /*N2 */

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
    amrex::Real sor[6]; /* temporary storage */
    speciesEntropy(sor, tc);

    /*Compute Eq 42 */
    for (id = 0; id < 6; ++id) {
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
    amrex::Real sor[6]; /* temporary storage */
    amrex::Real x[6]; /* need a ytx conversion */
    amrex::Real YOW = 0; /*See Eq 4, 6 in CK Manual */
    amrex::Real imw[6];

    get_imw(imw);

    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*O2 */
    YOW += y[1]*imw[1]; /*H2O */
    YOW += y[2]*imw[2]; /*CH4 */
    YOW += y[3]*imw[3]; /*CO */
    YOW += y[4]*imw[4]; /*CO2 */
    YOW += y[5]*imw[5]; /*N2 */
    /*Now compute y to x conversion */
    x[0] = y[0]/(31.998800*YOW); 
    x[1] = y[1]/(18.015340*YOW); 
    x[2] = y[2]/(16.043030*YOW); 
    x[3] = y[3]/(28.010550*YOW); 
    x[4] = y[4]/(44.009950*YOW); 
    x[5] = y[5]/(28.013400*YOW); 
    speciesEntropy(sor, tc);
    /*Perform computation in Eq 42 and 43 */
    result += x[0]*(sor[0]-log((x[0]+1e-100))-logPratio);
    result += x[1]*(sor[1]-log((x[1]+1e-100))-logPratio);
    result += x[2]*(sor[2]-log((x[2]+1e-100))-logPratio);
    result += x[3]*(sor[3]-log((x[3]+1e-100))-logPratio);
    result += x[4]*(sor[4]-log((x[4]+1e-100))-logPratio);
    result += x[5]*(sor[5]-log((x[5]+1e-100))-logPratio);
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
    XW += x[0]*31.998800; /*O2 */
    XW += x[1]*18.015340; /*H2O */
    XW += x[2]*16.043030; /*CH4 */
    XW += x[3]*28.010550; /*CO */
    XW += x[4]*44.009950; /*CO2 */
    XW += x[5]*28.013400; /*N2 */
    *P = *rho * 8.31446261815324e+07 * (*T) / XW; /*P = rho*R*T/W */

    return;
}


/*Compute P = rhoRT/W(y) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKPY(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  P)
{
    amrex::Real YOW = 0;/* for computing mean MW */
    amrex::Real imw[6];

    get_imw(imw);

    YOW += y[0]*imw[0]; /*O2 */
    YOW += y[1]*imw[1]; /*H2O */
    YOW += y[2]*imw[2]; /*CH4 */
    YOW += y[3]*imw[3]; /*CO */
    YOW += y[4]*imw[4]; /*CO2 */
    YOW += y[5]*imw[5]; /*N2 */
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
    W += c[0]*31.998800; /*O2 */
    W += c[1]*18.015340; /*H2O */
    W += c[2]*16.043030; /*CH4 */
    W += c[3]*28.010550; /*CO */
    W += c[4]*44.009950; /*CO2 */
    W += c[5]*28.013400; /*N2 */

    for (id = 0; id < 6; ++id) {
        sumC += c[id];
    }
    *P = *rho * 8.31446261815324e+07 * (*T) * sumC / W; /*P = rho*R*T/W */

    return;
}


/*Compute rho = PW(x)/RT */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKRHOX(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x,  amrex::Real *  rho)
{
    amrex::Real XW = 0;/* To hold mean molecular wt */
    XW += x[0]*31.998800; /*O2 */
    XW += x[1]*18.015340; /*H2O */
    XW += x[2]*16.043030; /*CH4 */
    XW += x[3]*28.010550; /*CO */
    XW += x[4]*44.009950; /*CO2 */
    XW += x[5]*28.013400; /*N2 */
    *rho = *P * XW / (8.31446261815324e+07 * (*T)); /*rho = P*W/(R*T) */

    return;
}


/*Compute rho = P*W(y)/RT */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKRHOY(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  rho)
{
    amrex::Real YOW = 0;
    amrex::Real tmp[6];
    amrex::Real imw[6];

    get_imw(imw);

    for (int i = 0; i < 6; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 6; i++)
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
    W += c[0]*31.998800; /*O2 */
    W += c[1]*18.015340; /*H2O */
    W += c[2]*16.043030; /*CH4 */
    W += c[3]*28.010550; /*CO */
    W += c[4]*44.009950; /*CO2 */
    W += c[5]*28.013400; /*N2 */

    for (id = 0; id < 6; ++id) {
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
    amrex::Real tmp[6];
    amrex::Real imw[6];

    get_imw(imw);

    for (int i = 0; i < 6; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 6; i++)
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
    XW += x[0]*31.998800; /*O2 */
    XW += x[1]*18.015340; /*H2O */
    XW += x[2]*16.043030; /*CH4 */
    XW += x[3]*28.010550; /*CO */
    XW += x[4]*44.009950; /*CO2 */
    XW += x[5]*28.013400; /*N2 */
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
    W += c[0]*31.998800; /*O2 */
    W += c[1]*18.015340; /*H2O */
    W += c[2]*16.043030; /*CH4 */
    W += c[3]*28.010550; /*CO */
    W += c[4]*44.009950; /*CO2 */
    W += c[5]*28.013400; /*N2 */

    for (id = 0; id < 6; ++id) {
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
    amrex::Real tmp[6];
    amrex::Real imw[6];

    get_imw(imw);

    for (int i = 0; i < 6; i++)
    {
        tmp[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 6; i++)
    {
        YOW += tmp[i];
    }

    amrex::Real YOWINV = 1.0/YOW;

    for (int i = 0; i < 6; i++)
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
    amrex::Real imw[6];

    get_imw(imw);

    /*Compute inverse of mean molecular wt first */
    for (int i = 0; i < 6; i++)
    {
        c[i] = y[i]*imw[i];
    }
    for (int i = 0; i < 6; i++)
    {
        YOW += c[i];
    }

    /*PW/RT (see Eq. 7) */
    PWORT = (*P)/(YOW * 8.31446261815324e+07 * (*T)); 
    /*Now compute conversion */

    for (int i = 0; i < 6; i++)
    {
        c[i] = PWORT * y[i] * imw[i];
    }
    return;
}


/*convert y[species] (mass fracs) to c[species] (molar conc) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKYTCR(amrex::Real *  rho, amrex::Real * /*T*/, amrex::Real * y,  amrex::Real * c)
{
    amrex::Real imw[6];

    get_imw(imw);

    for (int i = 0; i < 6; i++)
    {
        c[i] = (*rho)  * y[i] * imw[i];
    }
}


/*convert x[species] (mole fracs) to y[species] (mass fracs) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKXTY(amrex::Real *  x,  amrex::Real *  y)
{
    amrex::Real XW = 0; /*See Eq 4, 9 in CK Manual */
    /*Compute mean molecular wt first */
    XW += x[0]*31.998800; /*O2 */
    XW += x[1]*18.015340; /*H2O */
    XW += x[2]*16.043030; /*CH4 */
    XW += x[3]*28.010550; /*CO */
    XW += x[4]*44.009950; /*CO2 */
    XW += x[5]*28.013400; /*N2 */
    /*Now compute conversion */
    amrex::Real XWinv = 1.0/XW;
    y[0] = x[0]*31.998800*XWinv; 
    y[1] = x[1]*18.015340*XWinv; 
    y[2] = x[2]*16.043030*XWinv; 
    y[3] = x[3]*28.010550*XWinv; 
    y[4] = x[4]*44.009950*XWinv; 
    y[5] = x[5]*28.013400*XWinv; 

    return;
}


/*convert x[species] (mole fracs) to c[species] (molar conc) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKXTCP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x,  amrex::Real *  c)
{
    int id; /*loop counter */
    amrex::Real PORT = (*P)/(8.31446261815324e+07 * (*T)); /*P/RT */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 6; ++id) {
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
    XW += x[0]*31.998800; /*O2 */
    XW += x[1]*18.015340; /*H2O */
    XW += x[2]*16.043030; /*CH4 */
    XW += x[3]*28.010550; /*CO */
    XW += x[4]*44.009950; /*CO2 */
    XW += x[5]*28.013400; /*N2 */
    ROW = (*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 6; ++id) {
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
    for (id = 0; id < 6; ++id) {
        sumC += c[id];
    }

    /* See Eq 13  */
    amrex::Real sumCinv = 1.0/sumC;
    for (id = 0; id < 6; ++id) {
        x[id] = c[id]*sumCinv;
    }

    return;
}


/*convert c[species] (molar conc) to y[species] (mass fracs) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCTY(amrex::Real *  c, amrex::Real *  y)
{
    amrex::Real CW = 0; /*See Eq 12 in CK Manual */
    /*compute denominator in eq 12 first */
    CW += c[0]*31.998800; /*O2 */
    CW += c[1]*18.015340; /*H2O */
    CW += c[2]*16.043030; /*CH4 */
    CW += c[3]*28.010550; /*CO */
    CW += c[4]*44.009950; /*CO2 */
    CW += c[5]*28.013400; /*N2 */
    /*Now compute conversion */
    amrex::Real CWinv = 1.0/CW;
    y[0] = c[0]*31.998800*CWinv; 
    y[1] = c[1]*18.015340*CWinv; 
    y[2] = c[2]*16.043030*CWinv; 
    y[3] = c[3]*28.010550*CWinv; 
    y[4] = c[4]*44.009950*CWinv; 
    y[5] = c[5]*28.013400*CWinv; 

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
    for (id = 0; id < 6; ++id) {
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
    for (id = 0; id < 6; ++id) {
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
    for (id = 0; id < 6; ++id) {
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
    for (id = 0; id < 6; ++id) {
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
    for (id = 0; id < 6; ++id) {
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
    cvms[0] *= 2.598367006935648e+06; /*O2 */
    cvms[1] *= 4.615212712140454e+06; /*H2O */
    cvms[2] *= 5.182601178301878e+06; /*CH4 */
    cvms[3] *= 2.968332509769797e+06; /*CO */
    cvms[4] *= 1.889223372931176e+06; /*CO2 */
    cvms[5] *= 2.968030520448514e+06; /*N2 */
}


/*Returns the specific heats at constant pressure */
/*in mass units (Eq. 26) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKCPMS(amrex::Real *  T,  amrex::Real *  cpms)
{
    amrex::Real tT = *T; /*temporary temperature */
    const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    cp_R(cpms, tc);
    /*multiply by R/molecularweight */
    cpms[0] *= 2.598367006935648e+06; /*O2 */
    cpms[1] *= 4.615212712140454e+06; /*H2O */
    cpms[2] *= 5.182601178301878e+06; /*CH4 */
    cpms[3] *= 2.968332509769797e+06; /*CO */
    cpms[4] *= 1.889223372931176e+06; /*CO2 */
    cpms[5] *= 2.968030520448514e+06; /*N2 */
}


/*Returns internal energy in mass units (Eq 30.) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKUMS(amrex::Real *  T,  amrex::Real *  ums)
{
    amrex::Real tT = *T; /*temporary temperature */
    const amrex::Real tc[5] = { 0, tT, tT*tT, tT*tT*tT, tT*tT*tT*tT }; /*temperature cache */
    amrex::Real RT = 8.31446261815324e+07*tT; /*R*T */
    amrex::Real imw[6];

    get_imw(imw);

    speciesInternalEnergy(ums, tc);

    for (int i = 0; i < 6; i++)
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
    amrex::Real imw[6];

    get_imw(imw);

    speciesEnthalpy(hms, tc);

    for (int i = 0; i < 6; i++)
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
    sms[0] *= 2.598367006935648e+06; /*O2 */
    sms[1] *= 4.615212712140454e+06; /*H2O */
    sms[2] *= 5.182601178301878e+06; /*CH4 */
    sms[3] *= 2.968332509769797e+06; /*CO */
    sms[4] *= 1.889223372931176e+06; /*CO2 */
    sms[5] *= 2.968030520448514e+06; /*N2 */
}


/* GPU version of productionRate: no more use of thermo namespace vectors  */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void comp_qfqr(amrex::Real *  qf, amrex::Real * qr, amrex::Real * sc, amrex::Real * sc_qss, amrex::Real * tc, amrex::Real invT)
{

    /*reaction 1: 2.000000 CH4 + 3.000000 O2 => 2.000000 CO + 4.000000 H2O */
    qf[0] = sc[0]*sc[2];
    qr[0] = 0.0;

    /*reaction 2: 2.000000 CO + O2 => 2.000000 CO2 */
    qf[1] = pow(sc[0], 0.250000)*pow(sc[1], 0.500000)*sc[3];
    qr[1] = 0.0;

    /*reaction 3: 2.000000 CO2 => 2.000000 CO + O2 */
    qf[2] = sc[4];
    qr[2] = 0.0;

    /*compute the mixture concentration */
    amrex::Real mixture = 0.0;
    for (int i = 0; i < 6; ++i) {
        mixture += sc[i];
    }

    /*compute the Gibbs free energy */
    amrex::Real g_RT[6];
    gibbs(g_RT, tc);

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    amrex::Real refC = 101325 / 8.31446 * invT;
    amrex::Real refCinv = 1 / refC;

    /* Evaluate the kfs */
    amrex::Real k_f, k_r, Corr;

    // (0):  2.000000 CH4 + 3.000000 O2 => 2.000000 CO + 4.000000 H2O
    k_f = 1.545e-10
               * exp((0.5) * tc[0] - (20075.8288822793) * invT);
    qf[0] *= k_f;
    qr[0] *= k_f * exp(-(3.000000*g_RT[0] - 4.000000*g_RT[1] + 2.000000*g_RT[2] - 2.000000*g_RT[3])) * (refCinv);
    // (1):  2.000000 CO + O2 => 2.000000 CO2
    k_f = 199.1
               * exp(-(20128.6666321888) * invT);
    qf[1] *= k_f;
    qr[1] *= k_f * exp(-(g_RT[0] + 2.000000*g_RT[3] - 2.000000*g_RT[4])) * (refC);
    // (2):  2.000000 CO2 => 2.000000 CO + O2
    k_f = 250
               * exp(-(20128.6666321888) * invT);
    qf[2] *= k_f;
    qr[2] *= k_f * exp(-(-g_RT[0] - 2.000000*g_RT[3] + 2.000000*g_RT[4])) * (refCinv);


    return;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void  productionRate(amrex::Real * wdot, amrex::Real * sc, amrex::Real T)
{
    const amrex::Real tc[5] = { log(T), T, T*T, T*T*T, T*T*T*T }; // temperature cache
    const amrex::Real invT = 1.0 / tc[1];

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    const amrex::Real refC = 101325 / 8.31446 * invT;
    const amrex::Real refCinv = 1 / refC;

    for (int i = 0; i < 6; ++i) {
        wdot[i] = 0.0;
    }

    /*compute the mixture concentration */
    amrex::Real mixture = 0.0;
    for (int i = 0; i < 6; ++i) {
        mixture += sc[i];
    }

    /*compute the Gibbs free energy */
    amrex::Real g_RT[6];
    gibbs(g_RT, tc);

    {
        // (0):  2.000000 CH4 + 3.000000 O2 => 2.000000 CO + 4.000000 H2O
        const amrex::Real k_f = 1.545e-10
                   * exp((0.5) * tc[0] - (20075.8288822793) * invT);
        const amrex::Real qf = k_f * (sc[0]*sc[2]);
        const amrex::Real qr = k_f * exp(-(3.000000*g_RT[0] - 4.000000*g_RT[1] + 2.000000*g_RT[2] - 2.000000*g_RT[3])) * (refCinv) * (0.0);
        const amrex::Real qdot = qf - qr;
        wdot[0] -= 3.000000 * qdot;
        wdot[1] += 4.000000 * qdot;
        wdot[2] -= 2.000000 * qdot;
        wdot[3] += 2.000000 * qdot;
    }

    {
        // (1):  2.000000 CO + O2 => 2.000000 CO2
        const amrex::Real k_f = 199.1
                   * exp(-(20128.6666321888) * invT);
        const amrex::Real qf = k_f * (pow(sc[0], 0.250000)*pow(sc[1], 0.500000)*sc[3]);
        const amrex::Real qr = k_f * exp(-(g_RT[0] + 2.000000*g_RT[3] - 2.000000*g_RT[4])) * (refC) * (0.0);
        const amrex::Real qdot = qf - qr;
        wdot[0] -= qdot;
        wdot[3] -= 2.000000 * qdot;
        wdot[4] += 2.000000 * qdot;
    }

    {
        // (2):  2.000000 CO2 => 2.000000 CO + O2
        const amrex::Real k_f = 250
                   * exp(-(20128.6666321888) * invT);
        const amrex::Real qf = k_f * (sc[4]);
        const amrex::Real qr = k_f * exp(-(-g_RT[0] - 2.000000*g_RT[3] + 2.000000*g_RT[4])) * (refCinv) * (0.0);
        const amrex::Real qdot = qf - qr;
        wdot[0] += qdot;
        wdot[3] += 2.000000 * qdot;
        wdot[4] -= 2.000000 * qdot;
    }


    return;
}



/*compute the production rate for each species */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWC(amrex::Real *  T, amrex::Real *  C,  amrex::Real *  wdot)
{
    int id; /*loop counter */

    /*convert to SI */
    for (id = 0; id < 6; ++id) {
        C[id] *= 1.0e6;
    }

    /*convert to chemkin units */
    productionRate(wdot, C, *T);

    /*convert to chemkin units */
    for (id = 0; id < 6; ++id) {
        C[id] *= 1.0e-6;
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mass fractions */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWYP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  wdot)
{
    int id; /*loop counter */
    amrex::Real c[6]; /*temporary storage */
    amrex::Real YOW = 0; 
    amrex::Real PWORT; 
    amrex::Real imw[6];

    get_imw(imw);

    /*Compute inverse of mean molecular wt first */
    YOW += y[0]*imw[0]; /*O2 */
    YOW += y[1]*imw[1]; /*H2O */
    YOW += y[2]*imw[2]; /*CH4 */
    YOW += y[3]*imw[3]; /*CO */
    YOW += y[4]*imw[4]; /*CO2 */
    YOW += y[5]*imw[5]; /*N2 */
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

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 6; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given P, T, and mole fractions */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWXP(amrex::Real *  P, amrex::Real *  T, amrex::Real *  x,  amrex::Real *  wdot)
{
    int id; /*loop counter */
    amrex::Real c[6]; /*temporary storage */
    amrex::Real PORT = 1e6 * (*P)/(8.31446261815324e+07 * (*T)); /*1e6 * P/RT so c goes to SI units */

    /*Compute conversion, see Eq 10 */
    for (id = 0; id < 6; ++id) {
        c[id] = x[id]*PORT;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 6; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mass fractions */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWYR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  y,  amrex::Real *  wdot)
{
    int id; /*loop counter */
    amrex::Real c[6]; /*temporary storage */
    amrex::Real imw[6];

    get_imw(imw);

    /*See Eq 8 with an extra 1e6 so c goes to SI */
    c[0] = 1e6 * (*rho) * y[0]*imw[0]; 
    c[1] = 1e6 * (*rho) * y[1]*imw[1]; 
    c[2] = 1e6 * (*rho) * y[2]*imw[2]; 
    c[3] = 1e6 * (*rho) * y[3]*imw[3]; 
    c[4] = 1e6 * (*rho) * y[4]*imw[4]; 
    c[5] = 1e6 * (*rho) * y[5]*imw[5]; 

    /*call productionRate */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 6; ++id) {
        wdot[id] *= 1.0e-6;
    }
}


/*Returns the molar production rate of species */
/*Given rho, T, and mole fractions */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void CKWXR(amrex::Real *  rho, amrex::Real *  T, amrex::Real *  x,  amrex::Real *  wdot)
{
    int id; /*loop counter */
    amrex::Real c[6]; /*temporary storage */
    amrex::Real XW = 0; /*See Eq 4, 11 in CK Manual */
    amrex::Real ROW; 
    /*Compute mean molecular wt first */
    XW += x[0]*31.998800; /*O2 */
    XW += x[1]*18.015340; /*H2O */
    XW += x[2]*16.043030; /*CH4 */
    XW += x[3]*28.010550; /*CO */
    XW += x[4]*44.009950; /*CO2 */
    XW += x[5]*28.013400; /*N2 */
    /*Extra 1e6 factor to take c to SI */
    ROW = 1e6*(*rho) / XW;

    /*Compute conversion, see Eq 11 */
    for (id = 0; id < 6; ++id) {
        c[id] = x[id]*ROW;
    }

    /*convert to chemkin units */
    productionRate(wdot, c, *T);

    /*convert to chemkin units */
    for (id = 0; id < 6; ++id) {
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
        /*species 0: O2 */
        species[0] =
            -2.99673416e-03
            +1.96946040e-05 * tc[1]
            -2.90438853e-08 * tc[2]
            +1.29749135e-11 * tc[3];
        /*species 1: H2O */
        species[1] =
            -2.03643410e-03
            +1.30408042e-05 * tc[1]
            -1.64639119e-08 * tc[2]
            +7.08791268e-12 * tc[3];
        /*species 2: CH4 */
        species[2] =
            -1.36709788e-02
            +9.83601198e-05 * tc[1]
            -1.45422908e-07 * tc[2]
            +6.66775824e-11 * tc[3];
        /*species 3: CO */
        species[3] =
            -6.10353680e-04
            +2.03362866e-06 * tc[1]
            +2.72101765e-09 * tc[2]
            -3.61769800e-12 * tc[3];
        /*species 4: CO2 */
        species[4] =
            +8.98459677e-03
            -1.42471254e-05 * tc[1]
            +7.37757066e-09 * tc[2]
            -5.74798192e-13 * tc[3];
        /*species 5: N2 */
        species[5] =
            +1.40824040e-03
            -7.92644400e-06 * tc[1]
            +1.69245450e-08 * tc[2]
            -9.77941600e-12 * tc[3];
    } else {
        /*species 0: O2 */
        species[0] =
            +1.48308754e-03
            -1.51593334e-06 * tc[1]
            +6.28411665e-10 * tc[2]
            -8.66871176e-14 * tc[3];
        /*species 1: H2O */
        species[1] =
            +2.17691804e-03
            -3.28145036e-07 * tc[1]
            -2.91125961e-10 * tc[2]
            +6.72803968e-14 * tc[3];
        /*species 2: CH4 */
        species[2] =
            +1.33909467e-02
            -1.14657162e-05 * tc[1]
            +3.66877605e-09 * tc[2]
            -4.07260920e-13 * tc[3];
        /*species 3: CO */
        species[3] =
            +2.06252743e-03
            -1.99765154e-06 * tc[1]
            +6.90159024e-10 * tc[2]
            -8.14590864e-14 * tc[3];
        /*species 4: CO2 */
        species[4] =
            +4.41437026e-03
            -4.42962808e-06 * tc[1]
            +1.57047056e-09 * tc[2]
            -1.88833666e-13 * tc[3];
        /*species 5: N2 */
        species[5] =
            +1.48797680e-03
            -1.13695200e-06 * tc[1]
            +3.02911140e-10 * tc[2]
            -2.70134040e-14 * tc[3];
    }
    return;
}

/*compute an approx to the reaction Jacobian */
#ifdef COMPILE_JACOBIAN
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void aJacobian_precond(amrex::Real *  J, amrex::Real *  sc, amrex::Real T, int HP)
{
    for (int i=0; i<49; i++) {
        J[i] = 0.0;
    }

    amrex::Real wdot[6];
    for (int k=0; k<6; k++) {
        wdot[k] = 0.0;
    }

    const amrex::Real tc[5] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    amrex::Real invT = 1.0 / tc[1];
    amrex::Real invT2 = invT * invT;

    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    amrex::Real refC = 101325 / 8.31446 / T;
    amrex::Real refCinv = 1.0 / refC;

    /*compute the mixture concentration */
    amrex::Real mixture = 0.0;
    for (int k = 0; k < 6; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    amrex::Real g_RT[6];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    amrex::Real h_RT[6];
    speciesEnthalpy(h_RT, tc);

    amrex::Real phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    amrex::Real dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    amrex::Real dqdci, dcdc_fac, dqdc[6];
    amrex::Real Pr, fPr, F, k_0, logPr;
    amrex::Real logFcent, troe_c, troe_n, troePr_den, troePr, troe;
    amrex::Real Fcent1, Fcent2, Fcent3, Fcent;
    amrex::Real dlogFdc, dlogFdn, dlogFdcn_fac;
    amrex::Real dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
    const amrex::Real ln10 = log(10.0);
    const amrex::Real log10e = 1.0/log(10.0);
    /*reaction 1: 2.000000 CH4 + 3.000000 O2 => 2.000000 CO + 4.000000 H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[0], 3.000000)*(sc[2] * sc[2]);
    k_f = 1e-24 * 154500000000000
                * exp(0.5 * tc[0] - 0.50321666580472 * (39895) * invT);
    dlnkfdT = 0.5 * invT + 0.50321666580472 *  (39895)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[0] -= 3 * q; /* O2 */
    wdot[1] += 4 * q; /* H2O */
    wdot[2] -= 2 * q; /* CH4 */
    wdot[3] += 2 * q; /* CO */
    /* d()/d[O2] */
    dqdci =  + k_f*3.000000*pow(sc[0],2.000000)*pow(sc[2], 2.000000);
    J[0] += -3 * dqdci;           /* dwdot[O2]/d[O2] */
    J[1] += 4 * dqdci;            /* dwdot[H2O]/d[O2] */
    J[2] += -2 * dqdci;           /* dwdot[CH4]/d[O2] */
    J[3] += 2 * dqdci;            /* dwdot[CO]/d[O2] */
    /* d()/d[CH4] */
    dqdci =  + k_f*pow(sc[0], 3.000000)*2.000000*sc[2];
    J[14] += -3 * dqdci;          /* dwdot[O2]/d[CH4] */
    J[15] += 4 * dqdci;           /* dwdot[H2O]/d[CH4] */
    J[16] += -2 * dqdci;          /* dwdot[CH4]/d[CH4] */
    J[17] += 2 * dqdci;           /* dwdot[CO]/d[CH4] */
    /* d()/dT */
    J[42] += -3 * dqdT;           /* dwdot[O2]/dT */
    J[43] += 4 * dqdT;            /* dwdot[H2O]/dT */
    J[44] += -2 * dqdT;           /* dwdot[CH4]/dT */
    J[45] += 2 * dqdT;            /* dwdot[CO]/dT */

    /*reaction 2: 2.000000 CO + O2 => 2.000000 CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*(sc[3] * sc[3]);
    k_f = 1e-12 * 199100000000000
                * exp(0 * tc[0] - 0.50321666580472 * (40000) * invT);
    dlnkfdT = 0 * invT + 0.50321666580472 *  (40000)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[0] -= q; /* O2 */
    wdot[3] -= 2 * q; /* CO */
    wdot[4] += 2 * q; /* CO2 */
    /* d()/d[O2] */
    dqdci =  + k_f*pow(sc[3], 2.000000);
    J[0] -= dqdci;                /* dwdot[O2]/d[O2] */
    J[3] += -2 * dqdci;           /* dwdot[CO]/d[O2] */
    J[4] += 2 * dqdci;            /* dwdot[CO2]/d[O2] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[0]*2.000000*sc[3];
    J[21] -= dqdci;               /* dwdot[O2]/d[CO] */
    J[24] += -2 * dqdci;          /* dwdot[CO]/d[CO] */
    J[25] += 2 * dqdci;           /* dwdot[CO2]/d[CO] */
    /* d()/dT */
    J[42] -= dqdT;                /* dwdot[O2]/dT */
    J[45] += -2 * dqdT;           /* dwdot[CO]/dT */
    J[46] += 2 * dqdT;            /* dwdot[CO2]/dT */

    /*reaction 3: 2.000000 CO2 => 2.000000 CO + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = (sc[4] * sc[4]);
    k_f = 1e-06 * 250000000
                * exp(0 * tc[0] - 0.50321666580472 * (40000) * invT);
    dlnkfdT = 0 * invT + 0.50321666580472 *  (40000)  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[0] += q; /* O2 */
    wdot[3] += 2 * q; /* CO */
    wdot[4] -= 2 * q; /* CO2 */
    /* d()/d[CO2] */
    dqdci =  + k_f*2.000000*sc[4];
    J[28] += dqdci;               /* dwdot[O2]/d[CO2] */
    J[31] += 2 * dqdci;           /* dwdot[CO]/d[CO2] */
    J[32] += -2 * dqdci;          /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[42] += dqdT;                /* dwdot[O2]/dT */
    J[45] += 2 * dqdT;            /* dwdot[CO]/dT */
    J[46] += -2 * dqdT;           /* dwdot[CO2]/dT */

    amrex::Real c_R[6], dcRdT[6], e_RT[6];
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
    for (int k = 0; k < 6; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[42+k];
    }

    amrex::Real cmixinv = 1.0/cmix;
    amrex::Real tmp1 = ehmix*cmixinv;
    amrex::Real tmp3 = cmixinv*T;
    amrex::Real tmp2 = tmp1*tmp3;
    amrex::Real dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 6; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 6; ++m) {
            dehmixdc += eh_RT[m]*J[k*7+m];
        }
        J[k*7+6] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[48] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;
}
#endif

/*compute an approx to the reaction Jacobian (for preconditioning) */
#ifdef COMPILE_JACOBIAN
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void DWDOT_SIMPLIFIED(amrex::Real *  J, amrex::Real *  sc, amrex::Real *  Tp, int * HP)
{
    amrex::Real c[6];

    for (int k=0; k<6; k++) {
        c[k] = 1.e6 * sc[k];
    }

    aJacobian_precond(J, c, *Tp, *HP);

    /* dwdot[k]/dT */
    /* dTdot/d[X] */
    for (int k=0; k<6; k++) {
        J[42+k] *= 1.e-6;
        J[k*7+6] *= 1.e6;
    }

    return;
}
#endif

/*compute the reaction Jacobian on GPU */
#ifdef COMPILE_JACOBIAN
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void aJacobian(amrex::Real * J, amrex::Real * sc, amrex::Real T, int consP)
{

    for (int i=0; i<49; i++) {
        J[i] = 0.0;
    }

    amrex::Real wdot[6];
    for (int k=0; k<6; k++) {
        wdot[k] = 0.0;
    }

    const amrex::Real tc[5] = { log(T), T, T*T, T*T*T, T*T*T*T }; /*temperature cache */
    amrex::Real invT = 1.0 / tc[1];
    amrex::Real invT2 = invT * invT;


    /*reference concentration: P_atm / (RT) in inverse mol/m^3 */
    amrex::Real refC = 101325 / 8.31446 / T;
    amrex::Real refCinv = 1.0 / refC;

    /*compute the mixture concentration */
    amrex::Real mixture = 0.0;
    for (int k = 0; k < 6; ++k) {
        mixture += sc[k];
    }

    /*compute the Gibbs free energy */
    amrex::Real g_RT[6];
    gibbs(g_RT, tc);

    /*compute the species enthalpy */
    amrex::Real h_RT[6];
    speciesEnthalpy(h_RT, tc);

    amrex::Real phi_f, k_f, k_r, phi_r, Kc, q, q_nocor, Corr, alpha;
    amrex::Real dlnkfdT, dlnk0dT, dlnKcdT, dkrdT, dqdT;
    amrex::Real dqdci, dcdc_fac, dqdc[6];
    amrex::Real Pr, fPr, F, k_0, logPr;
    amrex::Real logFcent, troe_c, troe_n, troePr_den, troePr, troe;
    amrex::Real Fcent1, Fcent2, Fcent3, Fcent;
    amrex::Real dlogFdc, dlogFdn, dlogFdcn_fac;
    amrex::Real dlogPrdT, dlogfPrdT, dlogFdT, dlogFcentdT, dlogFdlogPr, dlnCorrdT;
    const amrex::Real ln10 = log(10.0);
    const amrex::Real log10e = 1.0/log(10.0);
    /*reaction 1: 2.000000 CH4 + 3.000000 O2 => 2.000000 CO + 4.000000 H2O */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = pow(sc[0], 3.000000)*(sc[2] * sc[2]);
    k_f = 1e-24 * 154500000000000
                * exp(0.5 * tc[0] - 0.50321666580472 * (39895) * invT);
    dlnkfdT = 0.5 * invT + 0.50321666580472 *  39895  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[0] -= 3 * q; /* O2 */
    wdot[1] += 4 * q; /* H2O */
    wdot[2] -= 2 * q; /* CH4 */
    wdot[3] += 2 * q; /* CO */
    /* d()/d[O2] */
    dqdci =  + k_f*3.000000*pow(sc[0],2.000000)*pow(sc[2], 2.000000);
    J[0] += -3 * dqdci;           /* dwdot[O2]/d[O2] */
    J[1] += 4 * dqdci;            /* dwdot[H2O]/d[O2] */
    J[2] += -2 * dqdci;           /* dwdot[CH4]/d[O2] */
    J[3] += 2 * dqdci;            /* dwdot[CO]/d[O2] */
    /* d()/d[CH4] */
    dqdci =  + k_f*pow(sc[0], 3.000000)*2.000000*sc[2];
    J[14] += -3 * dqdci;          /* dwdot[O2]/d[CH4] */
    J[15] += 4 * dqdci;           /* dwdot[H2O]/d[CH4] */
    J[16] += -2 * dqdci;          /* dwdot[CH4]/d[CH4] */
    J[17] += 2 * dqdci;           /* dwdot[CO]/d[CH4] */
    /* d()/dT */
    J[42] += -3 * dqdT;           /* dwdot[O2]/dT */
    J[43] += 4 * dqdT;            /* dwdot[H2O]/dT */
    J[44] += -2 * dqdT;           /* dwdot[CH4]/dT */
    J[45] += 2 * dqdT;            /* dwdot[CO]/dT */

    /*reaction 2: 2.000000 CO + O2 => 2.000000 CO2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = sc[0]*(sc[3] * sc[3]);
    k_f = 1e-12 * 199100000000000
                * exp(0 * tc[0] - 0.50321666580472 * (40000) * invT);
    dlnkfdT = 0 * invT + 0.50321666580472 *  40000  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[0] -= q; /* O2 */
    wdot[3] -= 2 * q; /* CO */
    wdot[4] += 2 * q; /* CO2 */
    /* d()/d[O2] */
    dqdci =  + k_f*pow(sc[3], 2.000000);
    J[0] -= dqdci;                /* dwdot[O2]/d[O2] */
    J[3] += -2 * dqdci;           /* dwdot[CO]/d[O2] */
    J[4] += 2 * dqdci;            /* dwdot[CO2]/d[O2] */
    /* d()/d[CO] */
    dqdci =  + k_f*sc[0]*2.000000*sc[3];
    J[21] -= dqdci;               /* dwdot[O2]/d[CO] */
    J[24] += -2 * dqdci;          /* dwdot[CO]/d[CO] */
    J[25] += 2 * dqdci;           /* dwdot[CO2]/d[CO] */
    /* d()/dT */
    J[42] -= dqdT;                /* dwdot[O2]/dT */
    J[45] += -2 * dqdT;           /* dwdot[CO]/dT */
    J[46] += 2 * dqdT;            /* dwdot[CO2]/dT */

    /*reaction 3: 2.000000 CO2 => 2.000000 CO + O2 */
    /*a non-third-body and non-pressure-fall-off reaction */
    /* forward */
    phi_f = (sc[4] * sc[4]);
    k_f = 1e-06 * 250000000
                * exp(0 * tc[0] - 0.50321666580472 * (40000) * invT);
    dlnkfdT = 0 * invT + 0.50321666580472 *  40000  * invT2;
    /* rate of progress */
    q = k_f*phi_f;
    dqdT = dlnkfdT*k_f*phi_f;
    /* update wdot */
    wdot[0] += q; /* O2 */
    wdot[3] += 2 * q; /* CO */
    wdot[4] -= 2 * q; /* CO2 */
    /* d()/d[CO2] */
    dqdci =  + k_f*2.000000*sc[4];
    J[28] += dqdci;               /* dwdot[O2]/d[CO2] */
    J[31] += 2 * dqdci;           /* dwdot[CO]/d[CO2] */
    J[32] += -2 * dqdci;          /* dwdot[CO2]/d[CO2] */
    /* d()/dT */
    J[42] += dqdT;                /* dwdot[O2]/dT */
    J[45] += 2 * dqdT;            /* dwdot[CO]/dT */
    J[46] += -2 * dqdT;           /* dwdot[CO2]/dT */

    amrex::Real c_R[6], dcRdT[6], e_RT[6];
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
    for (int k = 0; k < 6; ++k) {
        cmix += c_R[k]*sc[k];
        dcmixdT += dcRdT[k]*sc[k];
        ehmix += eh_RT[k]*wdot[k];
        dehmixdT += invT*(c_R[k]-eh_RT[k])*wdot[k] + eh_RT[k]*J[42+k];
    }

    amrex::Real cmixinv = 1.0/cmix;
    amrex::Real tmp1 = ehmix*cmixinv;
    amrex::Real tmp3 = cmixinv*T;
    amrex::Real tmp2 = tmp1*tmp3;
    amrex::Real dehmixdc;
    /* dTdot/d[X] */
    for (int k = 0; k < 6; ++k) {
        dehmixdc = 0.0;
        for (int m = 0; m < 6; ++m) {
            dehmixdc += eh_RT[m]*J[k*7+m];
        }
        J[k*7+6] = tmp2*c_R[k] - tmp3*dehmixdc;
    }
    /* dTdot/dT */
    J[48] = -tmp1 + tmp2*dcmixdT - tmp3*dehmixdT;

return;
}
#endif


/*compute the reaction Jacobian */
#ifdef COMPILE_JACOBIAN
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void DWDOT(amrex::Real *  J, amrex::Real *  sc, amrex::Real *  Tp, int * consP)
{
    amrex::Real c[6];

    for (int k=0; k<6; k++) {
        c[k] = 1.e6 * sc[k];
    }

    aJacobian(J, c, *Tp, *consP);

    /* dwdot[k]/dT */
    /* dTdot/d[X] */
    for (int k=0; k<6; k++) {
        J[42+k] *= 1.e-6;
        J[k*7+6] *= 1.e6;
    }

    return;
}
#endif


/* Transport function declarations  */


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetLENIMC(int* LENIMC ) {
    *LENIMC = 24;}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetLENRMC(int* LENRMC ) {
    *LENRMC = 846;}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetNO(int* NO ) {
    *NO = 4;}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetKK(int* KK ) {
    *KK = 6;}


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetNLITE(int* NLITE ) {
    *NLITE = 0;}


/*Patm in ergs/cm3 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetPATM(amrex::Real* PATM) {
    *PATM =   0.1013250000000000E+07;}


/*the molecular weights in g/mol */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetWT(amrex::Real* WT ) {
    WT[0] = 3.19988000E+01;
    WT[1] = 1.80153400E+01;
    WT[2] = 1.60430300E+01;
    WT[3] = 2.80105500E+01;
    WT[4] = 4.40099500E+01;
    WT[5] = 2.80134000E+01;
}


/*the lennard-jones potential well depth eps/kb in K */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetEPS(amrex::Real* EPS ) {
    EPS[0] = 1.07400000E+02;
    EPS[1] = 5.72400000E+02;
    EPS[2] = 1.41400000E+02;
    EPS[3] = 9.81000000E+01;
    EPS[4] = 2.44000000E+02;
    EPS[5] = 9.75300000E+01;
}


/*the lennard-jones collision diameter in Angstroms */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetSIG(amrex::Real* SIG ) {
    SIG[0] = 3.45800000E+00;
    SIG[1] = 2.60500000E+00;
    SIG[2] = 3.74600000E+00;
    SIG[3] = 3.65000000E+00;
    SIG[4] = 3.76300000E+00;
    SIG[5] = 3.62100000E+00;
}


/*the dipole moment in Debye */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetDIP(amrex::Real* DIP ) {
    DIP[0] = 0.00000000E+00;
    DIP[1] = 1.84400000E+00;
    DIP[2] = 0.00000000E+00;
    DIP[3] = 0.00000000E+00;
    DIP[4] = 0.00000000E+00;
    DIP[5] = 0.00000000E+00;
}


/*the polarizability in cubic Angstroms */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetPOL(amrex::Real* POL ) {
    POL[0] = 1.60000000E+00;
    POL[1] = 0.00000000E+00;
    POL[2] = 2.60000000E+00;
    POL[3] = 1.95000000E+00;
    POL[4] = 2.65000000E+00;
    POL[5] = 1.76000000E+00;
}


/*the rotational relaxation collision number at 298 K */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetZROT(amrex::Real* ZROT ) {
    ZROT[0] = 3.80000000E+00;
    ZROT[1] = 4.00000000E+00;
    ZROT[2] = 1.30000000E+01;
    ZROT[3] = 1.80000000E+00;
    ZROT[4] = 2.10000000E+00;
    ZROT[5] = 4.00000000E+00;
}


/*0: monoatomic, 1: linear, 2: nonlinear */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetNLIN(int* NLIN) {
    NLIN[0] = 1;
    NLIN[1] = 2;
    NLIN[2] = 2;
    NLIN[3] = 1;
    NLIN[4] = 1;
    NLIN[5] = 1;
}


/*Poly fits for the viscosities, dim NO*KK */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetCOFETA(amrex::Real* COFETA) {
    COFETA[0] = -1.68118868E+01;
    COFETA[1] = 2.52362554E+00;
    COFETA[2] = -2.49309128E-01;
    COFETA[3] = 1.10211025E-02;
    COFETA[4] = -1.17770937E+01;
    COFETA[5] = -8.26742721E-01;
    COFETA[6] = 3.39009079E-01;
    COFETA[7] = -2.00674327E-02;
    COFETA[8] = -1.95453421E+01;
    COFETA[9] = 3.36385478E+00;
    COFETA[10] = -3.56948469E-01;
    COFETA[11] = 1.56210922E-02;
    COFETA[12] = -1.63031240E+01;
    COFETA[13] = 2.26143219E+00;
    COFETA[14] = -2.15114671E-01;
    COFETA[15] = 9.53461976E-03;
    COFETA[16] = -2.36749526E+01;
    COFETA[17] = 4.99775518E+00;
    COFETA[18] = -5.52687718E-01;
    COFETA[19] = 2.34353338E-02;
    COFETA[20] = -1.62526779E+01;
    COFETA[21] = 2.24839597E+00;
    COFETA[22] = -2.13428438E-01;
    COFETA[23] = 9.46192413E-03;
}


/*Poly fits for the conductivities, dim NO*KK */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetCOFLAM(amrex::Real* COFLAM) {
    COFLAM[0] = -3.01284291E+00;
    COFLAM[1] = 3.37554994E+00;
    COFLAM[2] = -3.43353119E-01;
    COFLAM[3] = 1.51043444E-02;
    COFLAM[4] = 2.28195645E+01;
    COFLAM[5] = -8.72278946E+00;
    COFLAM[6] = 1.49300487E+00;
    COFLAM[7] = -7.41524047E-02;
    COFLAM[8] = 8.43903811E+00;
    COFLAM[9] = -2.78020298E+00;
    COFLAM[10] = 7.09313623E-01;
    COFLAM[11] = -4.04300978E-02;
    COFLAM[12] = 9.92459822E+00;
    COFLAM[13] = -2.28318157E+00;
    COFLAM[14] = 4.73113746E-01;
    COFLAM[15] = -2.40056652E-02;
    COFLAM[16] = -1.24047589E+01;
    COFLAM[17] = 6.34783131E+00;
    COFLAM[18] = -6.37857884E-01;
    COFLAM[19] = 2.37613820E-02;
    COFLAM[20] = 1.15507063E+01;
    COFLAM[21] = -2.91452379E+00;
    COFLAM[22] = 5.55043580E-01;
    COFLAM[23] = -2.75172461E-02;
}


/*Poly fits for the diffusion coefficients, dim NO*KK*KK */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetCOFD(amrex::Real* COFD) {
    COFD[0] = -1.53110708E+01;
    COFD[1] = 3.37317428E+00;
    COFD[2] = -2.24900439E-01;
    COFD[3] = 9.81228151E-03;
    COFD[4] = -2.10640014E+01;
    COFD[5] = 5.50980695E+00;
    COFD[6] = -4.78335488E-01;
    COFD[7] = 1.98515434E-02;
    COFD[8] = -1.59878769E+01;
    COFD[9] = 3.66478157E+00;
    COFD[10] = -2.61506432E-01;
    COFD[11] = 1.13466231E-02;
    COFD[12] = -1.50371784E+01;
    COFD[13] = 3.26249588E+00;
    COFD[14] = -2.10658287E-01;
    COFD[15] = 9.20032462E-03;
    COFD[16] = -1.81197354E+01;
    COFD[17] = 4.33684042E+00;
    COFD[18] = -3.44981265E-01;
    COFD[19] = 1.48142449E-02;
    COFD[20] = -1.50096240E+01;
    COFD[21] = 3.25515933E+00;
    COFD[22] = -2.09710110E-01;
    COFD[23] = 9.15941830E-03;
    COFD[24] = -2.10640014E+01;
    COFD[25] = 5.50980695E+00;
    COFD[26] = -4.78335488E-01;
    COFD[27] = 1.98515434E-02;
    COFD[28] = -1.31492641E+01;
    COFD[29] = 1.48004311E+00;
    COFD[30] = 1.60499553E-01;
    COFD[31] = -1.19765679E-02;
    COFD[32] = -2.12953660E+01;
    COFD[33] = 5.52385268E+00;
    COFD[34] = -4.69683833E-01;
    COFD[35] = 1.90677265E-02;
    COFD[36] = -2.08943798E+01;
    COFD[37] = 5.44718652E+00;
    COFD[38] = -4.72082953E-01;
    COFD[39] = 1.96531321E-02;
    COFD[40] = -2.12021508E+01;
    COFD[41] = 5.20775052E+00;
    COFD[42] = -4.07348327E-01;
    COFD[43] = 1.55473283E-02;
    COFD[44] = -2.08123325E+01;
    COFD[45] = 5.42470154E+00;
    COFD[46] = -4.69700416E-01;
    COFD[47] = 1.95706904E-02;
    COFD[48] = -1.59878769E+01;
    COFD[49] = 3.66478157E+00;
    COFD[50] = -2.61506432E-01;
    COFD[51] = 1.13466231E-02;
    COFD[52] = -2.12953660E+01;
    COFD[53] = 5.52385268E+00;
    COFD[54] = -4.69683833E-01;
    COFD[55] = 1.90677265E-02;
    COFD[56] = -1.68532868E+01;
    COFD[57] = 4.00572564E+00;
    COFD[58] = -3.04255586E-01;
    COFD[59] = 1.31384083E-02;
    COFD[60] = -1.56953813E+01;
    COFD[61] = 3.54443207E+00;
    COFD[62] = -2.46133003E-01;
    COFD[63] = 1.06905018E-02;
    COFD[64] = -1.89594784E+01;
    COFD[65] = 4.68735220E+00;
    COFD[66] = -3.87449514E-01;
    COFD[67] = 1.65352261E-02;
    COFD[68] = -1.56756076E+01;
    COFD[69] = 3.54035770E+00;
    COFD[70] = -2.45653442E-01;
    COFD[71] = 1.06717969E-02;
    COFD[72] = -1.50371784E+01;
    COFD[73] = 3.26249588E+00;
    COFD[74] = -2.10658287E-01;
    COFD[75] = 9.20032462E-03;
    COFD[76] = -2.08943798E+01;
    COFD[77] = 5.44718652E+00;
    COFD[78] = -4.72082953E-01;
    COFD[79] = 1.96531321E-02;
    COFD[80] = -1.56953813E+01;
    COFD[81] = 3.54443207E+00;
    COFD[82] = -2.46133003E-01;
    COFD[83] = 1.06905018E-02;
    COFD[84] = -1.48061490E+01;
    COFD[85] = 3.16912473E+00;
    COFD[86] = -1.98792456E-01;
    COFD[87] = 8.69726395E-03;
    COFD[88] = -1.77673000E+01;
    COFD[89] = 4.20234040E+00;
    COFD[90] = -3.28057658E-01;
    COFD[91] = 1.41006192E-02;
    COFD[92] = -1.47850486E+01;
    COFD[93] = 3.16433919E+00;
    COFD[94] = -1.98191564E-01;
    COFD[95] = 8.67209742E-03;
    COFD[96] = -1.81197354E+01;
    COFD[97] = 4.33684042E+00;
    COFD[98] = -3.44981265E-01;
    COFD[99] = 1.48142449E-02;
    COFD[100] = -2.12021508E+01;
    COFD[101] = 5.20775052E+00;
    COFD[102] = -4.07348327E-01;
    COFD[103] = 1.55473283E-02;
    COFD[104] = -1.89594784E+01;
    COFD[105] = 4.68735220E+00;
    COFD[106] = -3.87449514E-01;
    COFD[107] = 1.65352261E-02;
    COFD[108] = -1.77673000E+01;
    COFD[109] = 4.20234040E+00;
    COFD[110] = -3.28057658E-01;
    COFD[111] = 1.41006192E-02;
    COFD[112] = -2.10907727E+01;
    COFD[113] = 5.29211327E+00;
    COFD[114] = -4.56068366E-01;
    COFD[115] = 1.91195062E-02;
    COFD[116] = -1.77350592E+01;
    COFD[117] = 4.19328271E+00;
    COFD[118] = -3.26911461E-01;
    COFD[119] = 1.40520357E-02;
    COFD[120] = -1.50096240E+01;
    COFD[121] = 3.25515933E+00;
    COFD[122] = -2.09710110E-01;
    COFD[123] = 9.15941830E-03;
    COFD[124] = -2.08123325E+01;
    COFD[125] = 5.42470154E+00;
    COFD[126] = -4.69700416E-01;
    COFD[127] = 1.95706904E-02;
    COFD[128] = -1.56756076E+01;
    COFD[129] = 3.54035770E+00;
    COFD[130] = -2.45653442E-01;
    COFD[131] = 1.06717969E-02;
    COFD[132] = -1.47850486E+01;
    COFD[133] = 3.16433919E+00;
    COFD[134] = -1.98191564E-01;
    COFD[135] = 8.67209742E-03;
    COFD[136] = -1.77350592E+01;
    COFD[137] = 4.19328271E+00;
    COFD[138] = -3.26911461E-01;
    COFD[139] = 1.40520357E-02;
    COFD[140] = -1.47639290E+01;
    COFD[141] = 3.15955654E+00;
    COFD[142] = -1.97590757E-01;
    COFD[143] = 8.64692156E-03;
}


/*List of specs with small weight, dim NLITE */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetKTDIF(int* KTDIF) {
}


/*Poly fits for thermal diff ratios, dim NO*NLITE*KK */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetCOFTD(amrex::Real* COFTD) {
}


/*compute the critical parameters for each species */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void GET_CRITPARAMS(amrex::Real *  Tci, amrex::Real *  ai, amrex::Real *  bi, amrex::Real *  acentric_i)
{

    amrex::Real   EPS[6];
    amrex::Real   SIG[6];
    amrex::Real    wt[6];
    amrex::Real avogadro = 6.02214199e23;
    amrex::Real boltzmann = 1.3806503e-16; //we work in CGS
    amrex::Real Rcst = 83.144598; //in bar [CGS] !

    egtransetEPS(EPS);
    egtransetSIG(SIG);
    get_mw(wt);

    /*species 0: O2 */
    /*Imported from NIST */
    Tci[0] = 154.581000 ; 
    ai[0] = 1e6 * 0.42748 * Rcst * Rcst * Tci[0] * Tci[0] / (31.998800 * 31.998800 * 50.430466); 
    bi[0] = 0.08664 * Rcst * Tci[0] / (31.998800 * 50.430466); 
    acentric_i[0] = 0.022200 ;

    /*species 1: H2O */
    /*Imported from NIST */
    Tci[1] = 647.096000 ; 
    ai[1] = 1e6 * 0.42748 * Rcst * Rcst * Tci[1] * Tci[1] / (18.015340 * 18.015340 * 220.640000); 
    bi[1] = 0.08664 * Rcst * Tci[1] / (18.015340 * 220.640000); 
    acentric_i[1] = 0.344300 ;

    /*species 2: CH4 */
    /*Imported from NIST */
    Tci[2] = 190.560000 ; 
    ai[2] = 1e6 * 0.42748 * Rcst * Rcst * Tci[2] * Tci[2] / (16.043030 * 16.043030 * 45.990000); 
    bi[2] = 0.08664 * Rcst * Tci[2] / (16.043030 * 45.990000); 
    acentric_i[2] = 0.011000 ;

    /*species 3: CO */
    /*Imported from NIST */
    Tci[3] = 132.850000 ; 
    ai[3] = 1e6 * 0.42748 * Rcst * Rcst * Tci[3] * Tci[3] / (28.010000 * 28.010000 * 34.940000); 
    bi[3] = 0.08664 * Rcst * Tci[3] / (28.010000 * 34.940000); 
    acentric_i[3] = 0.045000 ;

    /*species 4: CO2 */
    /*Imported from NIST */
    Tci[4] = 304.120000 ; 
    ai[4] = 1e6 * 0.42748 * Rcst * Rcst * Tci[4] * Tci[4] / (44.009950 * 44.009950 * 73.740000); 
    bi[4] = 0.08664 * Rcst * Tci[4] / (44.009950 * 73.740000); 
    acentric_i[4] = 0.225000 ;

    /*species 5: N2 */
    /*Imported from NIST */
    Tci[5] = 126.192000 ; 
    ai[5] = 1e6 * 0.42748 * Rcst * Rcst * Tci[5] * Tci[5] / (28.013400 * 28.013400 * 33.958000); 
    bi[5] = 0.08664 * Rcst * Tci[5] / (28.013400 * 33.958000); 
    acentric_i[5] = 0.037200 ;

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
