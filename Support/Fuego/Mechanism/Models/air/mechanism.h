#ifndef MECHANISM_h
#define MECHANISM_h
#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>

#if 0
/* Elements
0  O
1  N
*/
#endif

/* Species */
#define O2_ID 0
#define N2_ID 1

#define NUM_ELEMENTS 2
#define NUM_SPECIES 2
#define NUM_REACTIONS 0

#define NUM_FIT 4
/* Transport function declarations  */


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetLENIMC(int* LENIMC ) {
    *LENIMC = 8;};


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetLENRMC(int* LENRMC ) {
    *LENRMC = 130;};


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetNO(int* NO ) {
    *NO = 4;};


AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetKK(int* KK ) {
    *KK = 2;};


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
    WT[0] = 3.19988000E+01;
    WT[1] = 2.80134000E+01;
};


/*the lennard-jones potential well depth eps/kb in K */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetEPS(amrex::Real* EPS ) {
    EPS[0] = 1.07400000E+02;
    EPS[1] = 9.75300000E+01;
};


/*the lennard-jones collision diameter in Angstroms */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetSIG(amrex::Real* SIG ) {
    SIG[0] = 3.45800000E+00;
    SIG[1] = 3.62100000E+00;
};


/*the dipole moment in Debye */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetDIP(amrex::Real* DIP ) {
    DIP[0] = 0.00000000E+00;
    DIP[1] = 0.00000000E+00;
};


/*the polarizability in cubic Angstroms */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetPOL(amrex::Real* POL ) {
    POL[0] = 1.60000000E+00;
    POL[1] = 1.76000000E+00;
};


/*the rotational relaxation collision number at 298 K */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetZROT(amrex::Real* ZROT ) {
    ZROT[0] = 3.80000000E+00;
    ZROT[1] = 4.00000000E+00;
};


/*0: monoatomic, 1: linear, 2: nonlinear */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetNLIN(int* NLIN) {
    NLIN[0] = 1;
    NLIN[1] = 1;
};


/*Poly fits for the viscosities, dim NO*KK */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetCOFETA(amrex::Real* COFETA) {
    COFETA[0] = -1.60066324E+01;
    COFETA[1] = 2.16753735E+00;
    COFETA[2] = -1.97226850E-01;
    COFETA[3] = 8.50065468E-03;
    COFETA[4] = -1.55270326E+01;
    COFETA[5] = 1.92766908E+00;
    COFETA[6] = -1.66518287E-01;
    COFETA[7] = 7.19100649E-03;
};


/*Poly fits for the conductivities, dim NO*KK */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetCOFLAM(amrex::Real* COFLAM) {
    COFLAM[0] = -2.11869892E+00;
    COFLAM[1] = 2.98568651E+00;
    COFLAM[2] = -2.86879123E-01;
    COFLAM[3] = 1.23850873E-02;
    COFLAM[4] = 7.60997504E+00;
    COFLAM[5] = -1.18418698E+00;
    COFLAM[6] = 3.03558703E-01;
    COFLAM[7] = -1.54159597E-02;
};


/*Poly fits for the diffusion coefficients, dim NO*KK*KK */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void egtransetCOFD(amrex::Real* COFD) {
    COFD[0] = -1.47079646E+01;
    COFD[1] = 3.10657376E+00;
    COFD[2] = -1.85922460E-01;
    COFD[3] = 7.92680827E-03;
    COFD[4] = -1.44285949E+01;
    COFD[5] = 2.99858376E+00;
    COFD[6] = -1.72232643E-01;
    COFD[7] = 7.34804765E-03;
    COFD[8] = -1.44285949E+01;
    COFD[9] = 2.99858376E+00;
    COFD[10] = -1.72232643E-01;
    COFD[11] = 7.34804765E-03;
    COFD[12] = -1.42056656E+01;
    COFD[13] = 2.91297621E+00;
    COFD[14] = -1.61544771E-01;
    COFD[15] = 6.90271324E-03;
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
