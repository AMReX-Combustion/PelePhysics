#ifndef MECHANISM_h
#define MECHANISM_h
#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>

#if 0
/* Elements
0  C
1  H
2  O
3  N
*/
#endif

/* Species */
#define NC10H22_ID 0
#define O2_ID 1
#define N2_ID 2

#define NUM_ELEMENTS 4
#define NUM_SPECIES 3
#define NUM_REACTIONS 0

#define NUM_FIT 4
/* Transport function declarations  */


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetLENIMC(int* LENIMC ) {
    *LENIMC = 12;};


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetLENRMC(int* LENRMC ) {
    *LENRMC = 252;};


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetNO(int* NO ) {
    *NO = 4;};


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetKK(int* KK ) {
    *KK = 3;};


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetNLITE(int* NLITE ) {
    *NLITE = 0;};


/*Patm in ergs/cm3 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetPATM(amrex::Real* PATM) {
    *PATM =   0.1013250000000000E+07;};


/*the molecular weights in g/mol */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetWT(amrex::Real* WT ) {
    WT[0] = 1.42286840E+02;
    WT[1] = 3.19988000E+01;
    WT[2] = 2.80134000E+01;
};


/*the lennard-jones potential well depth eps/kb in K */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetEPS(amrex::Real* EPS ) {
    EPS[0] = 7.04917000E+02;
    EPS[1] = 1.07400000E+02;
    EPS[2] = 9.75300000E+01;
};


/*the lennard-jones collision diameter in Angstroms */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetSIG(amrex::Real* SIG ) {
    SIG[0] = 6.67500000E+00;
    SIG[1] = 3.45800000E+00;
    SIG[2] = 3.62100000E+00;
};


/*the dipole moment in Debye */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetDIP(amrex::Real* DIP ) {
    DIP[0] = 0.00000000E+00;
    DIP[1] = 0.00000000E+00;
    DIP[2] = 0.00000000E+00;
};


/*the polarizability in cubic Angstroms */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetPOL(amrex::Real* POL ) {
    POL[0] = 0.00000000E+00;
    POL[1] = 1.60000000E+00;
    POL[2] = 1.76000000E+00;
};


/*the rotational relaxation collision number at 298 K */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetZROT(amrex::Real* ZROT ) {
    ZROT[0] = 1.00000000E+00;
    ZROT[1] = 3.80000000E+00;
    ZROT[2] = 4.00000000E+00;
};


/*0: monoatomic, 1: linear, 2: nonlinear */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetNLIN(int* NLIN) {
    NLIN[0] = 2;
    NLIN[1] = 1;
    NLIN[2] = 1;
};


/*Poly fits for the viscosities, dim NO*KK */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetCOFETA(amrex::Real* COFETA) {
    COFETA[0] = -1.29489842E+01;
    COFETA[1] = -6.04647944E-01;
    COFETA[2] = 3.07320970E-01;
    COFETA[3] = -1.89227139E-02;
    COFETA[4] = -1.68118868E+01;
    COFETA[5] = 2.52362554E+00;
    COFETA[6] = -2.49309128E-01;
    COFETA[7] = 1.10211025E-02;
    COFETA[8] = -1.62526779E+01;
    COFETA[9] = 2.24839597E+00;
    COFETA[10] = -2.13428438E-01;
    COFETA[11] = 9.46192413E-03;
};


/*Poly fits for the conductivities, dim NO*KK */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetCOFLAM(amrex::Real* COFLAM) {
    COFLAM[0] = -1.38317203E+01;
    COFLAM[1] = 5.32397257E+00;
    COFLAM[2] = -3.02248916E-01;
    COFLAM[3] = 1.51838934E-03;
    COFLAM[4] = -3.01284291E+00;
    COFLAM[5] = 3.37554994E+00;
    COFLAM[6] = -3.43353119E-01;
    COFLAM[7] = 1.51043444E-02;
    COFLAM[8] = 1.15507063E+01;
    COFLAM[9] = -2.91452379E+00;
    COFLAM[10] = 5.55043580E-01;
    COFLAM[11] = -2.75172461E-02;
};


/*Poly fits for the diffusion coefficients, dim NO*KK*KK */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetCOFD(amrex::Real* COFD) {
    COFD[0] = -1.65247199E+01;
    COFD[1] = 1.86376261E+00;
    COFD[2] = 9.20701178E-02;
    COFD[3] = -8.51470516E-03;
    COFD[4] = -2.24760016E+01;
    COFD[5] = 5.50822096E+00;
    COFD[6] = -4.78234659E-01;
    COFD[7] = 1.98512785E-02;
    COFD[8] = -2.21643696E+01;
    COFD[9] = 5.41977231E+00;
    COFD[10] = -4.69168499E-01;
    COFD[11] = 1.95518107E-02;
    COFD[12] = -2.24760016E+01;
    COFD[13] = 5.50822096E+00;
    COFD[14] = -4.78234659E-01;
    COFD[15] = 1.98512785E-02;
    COFD[16] = -1.53110708E+01;
    COFD[17] = 3.37317428E+00;
    COFD[18] = -2.24900439E-01;
    COFD[19] = 9.81228151E-03;
    COFD[20] = -1.50096240E+01;
    COFD[21] = 3.25515933E+00;
    COFD[22] = -2.09710110E-01;
    COFD[23] = 9.15941830E-03;
    COFD[24] = -2.21643696E+01;
    COFD[25] = 5.41977231E+00;
    COFD[26] = -4.69168499E-01;
    COFD[27] = 1.95518107E-02;
    COFD[28] = -1.50096240E+01;
    COFD[29] = 3.25515933E+00;
    COFD[30] = -2.09710110E-01;
    COFD[31] = 9.15941830E-03;
    COFD[32] = -1.47639290E+01;
    COFD[33] = 3.15955654E+00;
    COFD[34] = -1.97590757E-01;
    COFD[35] = 8.64692156E-03;
};


/*List of specs with small weight, dim NLITE */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetKTDIF(int* KTDIF) {
};


/*Poly fits for thermal diff ratios, dim NO*NLITE*KK */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetCOFTD(amrex::Real* COFTD) {
};
#endif

/* End of file  */
