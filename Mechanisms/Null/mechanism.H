#ifndef MECHANISM_H
#define MECHANISM_H

#include <AMReX_Gpu.H>
#include <AMReX_REAL.H>

#if 0
/* Elements
0  X
*/
#endif

/* Species */
#define XO_ID 0

#define NUM_ELEMENTS 1
#define NUM_SPECIES 1
#define NUM_REACTIONS 0
#define NUM_FIT 4

/* ALWAYS on CPU stuff -- can have different def depending on if we are CPU or
 * GPU based. Defined in mechanism.cpp  */
/* MISC  */
void CKAWT(amrex::Real* awt);
void CKNCF(int* ncf);
void CKSYME_STR(amrex::Vector<std::string>& ename);
void CKSYMS_STR(amrex::Vector<std::string>& kname);
/* SPARSE INFORMATION  */
void SPARSITY_INFO(int* nJdata, const int* consP, int NCELLS);
void SPARSITY_INFO_SYST(int* nJdata, const int* consP, int NCELLS);
void SPARSITY_INFO_SYST_SIMPLIFIED(int* nJdata, const int* consP);
void
SPARSITY_PREPROC_CSC(int* rowVals, int* colPtrs, const int* consP, int NCELLS);
void SPARSITY_PREPROC_CSR(
  int* colVals, int* rowPtrs, const int* consP, int NCELLS, int base);
void SPARSITY_PREPROC_SYST_CSR(
  int* colVals, int* rowPtrs, const int* consP, int NCELLS, int base);
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(
  int* rowVals, int* colPtrs, int* indx, const int* consP);
void SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(
  int* colVals, int* rowPtr, const int* consP, int base);

//  inverse molecular weights
#ifdef AMREX_USE_GPU
AMREX_GPU_CONSTANT const amrex::Real global_imw[1] = {
  1.0,
};
#endif
const amrex::Real h_global_imw[1] = {
  1.0,
};

//  molecular weights
#ifdef AMREX_USE_GPU
AMREX_GPU_CONSTANT const amrex::Real global_mw[1] = {
  1.0,
};
#endif
const amrex::Real h_global_mw[1] = {
  1.0,
};

/* inverse molecular weights  */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
get_imw(amrex::Real imw_new[])
{
  imw_new[0] = 1.0;
}
//  inverse molecular weight
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
imw(const int n)
{
#if AMREX_DEVICE_COMPILE
  return global_imw[n];
#else
  return h_global_imw[n];
#endif
}

/* molecular weights  */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
get_mw(amrex::Real mw_new[])
{
  mw_new[0] = 1.0;
}

//  molecular weight
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
mw(const int n)
{
#if AMREX_DEVICE_COMPILE
  return global_mw[n];
#else
  return h_global_mw[n];
#endif
}

/* Returns R, Rc, Patm */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKRP(amrex::Real& ru, amrex::Real& ruc, amrex::Real& pa)
{
  ru = 8.31446261815324e+07;
  ruc = 1.98721558317399615845;
  pa = 1.01325e+06;
}

/*compute the production rate for each species */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKWC(const amrex::Real /*T*/, const amrex::Real* /*C*/, amrex::Real* wdot)
{
  wdot[0] = 0.0e0;
}

/*Returns the mean specific heat at CP (Eq. 34) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCPBS(const amrex::Real /*T*/, const amrex::Real* /*y*/, amrex::Real& /*cpbs*/)
{
  amrex::Abort("CKCPBS not implemented for Null chemistry.");
}

/*Returns the specific heats at constant pressure */
/*in mass units (Eq. 26) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCPMS(const amrex::Real /*T*/, amrex::Real* /*cpms*/)
{
  amrex::Abort("CKCPMS not implemented for Null chemistry.");
}

/*given y[species]: mass fractions */
/*returns mean molecular weight (gm/mole) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKMMWY(const amrex::Real* /*y*/, amrex::Real& /*wtm*/)
{
  amrex::Abort("CKMMWY not implemented for Null chemistry.");
}

/* get temperature given internal energy in mass units and mass fracs */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
GET_T_GIVEN_EY(
  const amrex::Real /*e*/,
  const amrex::Real* /*y[]*/,
  amrex::Real& /*t*/,
  int& /*ierr*/)
{
  amrex::Abort("GET_T_GIVEN_EY not implemented for Null chemistry.");
}
/* get temperature given enthalpy in mass units and mass fracs */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
GET_T_GIVEN_HY(
  const amrex::Real /*h*/,
  const amrex::Real* /*y*/,
  amrex::Real& /*t*/,
  int& /*ierr*/)
{
  amrex::Abort("GET_T_GIVEN_HY not implemented for Null chemistry.");
}

/*Compute P = rhoRT/W(y) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKPY(
  const amrex::Real /*rho*/,
  const amrex::Real /*T*/,
  const amrex::Real* /*y*/,
  amrex::Real& /*P*/)
{
  amrex::Abort("CKPY not implemented for Null chemistry.");
}

/*convert y[species] (mass fracs) to c[species] (molar conc) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKYTCR(
  const amrex::Real /*rho*/,
  const amrex::Real /*T*/,
  const amrex::Real* /*y*/,
  amrex::Real* /*c*/)
{
  amrex::Abort("CKYTCR not implemented for Null chemistry.");
}

/*compute the reaction Jacobian */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
DWDOT(
  amrex::Real* /*J*/,
  const amrex::Real* /*sc*/,
  const amrex::Real* /*Tp*/,
  const int* /*consP*/)
{
  amrex::Abort("DWDOT not implemented for Null chemistry.");
}

/*compute an approx to the reaction Jacobian (for preconditioning) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
DWDOT_SIMPLIFIED(
  amrex::Real* /*J*/,
  const amrex::Real* /*sc*/,
  const amrex::Real* /*Tp*/,
  const int* /*HP*/)
{
  amrex::Abort("DWDOT_SIMPLIFIED not implemented for Null chemistry.");
}

/*Returns internal energy in mass units (Eq 30.) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKUMS(const amrex::Real /*T*/, amrex::Real* /*ums*/)
{
  amrex::Abort("CKUMS not implemented for Null chemistry.");
}

// get mean internal energy in mass units
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKUBMS(
  const amrex::Real /*T*/, const amrex::Real* /*y[]*/, amrex::Real& /*ubms*/)
{
  amrex::Abort("CKUMBS is not implemented for Null chemsitry");
}

/*convert y[species] (mass fracs) to x[species] (mole fracs) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKYTX(const amrex::Real* /*y*/, amrex::Real* /*x*/)
{
  amrex::Abort("CKYTX not implemented for Null chemistry.");
}

/*Compute rho = P*W(y)/RT */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKRHOY(
  const amrex::Real /*P*/,
  const amrex::Real /*T*/,
  const amrex::Real* /*y*/,
  amrex::Real& /*rho*/)
{
  amrex::Abort("CKRHOY not implemented for Null chemistry.");
}

/*Returns the mean specific heat at CV (Eq. 36) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCVBS(const amrex::Real /*T*/, const amrex::Real* /*y*/, amrex::Real& /*cvbs*/)
{
  amrex::Abort("CKCVBS not implemented for Null chemistry.");
}

/*Returns enthalpy in mass units (Eq 27.) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKHMS(const amrex::Real /*T*/, amrex::Real* /*hms*/)
{
  amrex::Abort("CKHMS not implemented for Null chemistry.");
}

/*Returns the specific heats at constant volume */
/*in mass units (Eq. 29) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKCVMS(const amrex::Real /*T*/, amrex::Real* /*cvms*/)
{
  amrex::Abort("CKCVMS not implemented for Null chemistry.");
}

/*convert x[species] (mole fracs) to y[species] (mass fracs) */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKXTY(const amrex::Real* /*x*/, amrex::Real* /*y*/)
{
  amrex::Abort("CKXTY not implemented for Null chemistry.");
}

/*compute the critical parameters for each species */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
GET_CRITPARAMS(
  amrex::Real* /*Tci*/,
  amrex::Real* /*ai*/,
  amrex::Real* /*bi*/,
  amrex::Real* /*acentric_i*/)
{
  amrex::Abort("GET_CRITPARAMS not implemented for Null chemistry.");
}

// compute the critical parameter quantities for each species for SRK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
GET_CRITPARAMS_SRK(
  amrex::Real* /*sqrtOneOverTc*/,
  amrex::Real* /*sqrtAsti*/,
  amrex::Real* /*Bi*/,
  amrex::Real* /*Fomega*/)
{
  amrex::Abort("GET_CRITPARAMS_SRK not implemented for Null chemistry.");
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetLENIMC(int* /*LENIMC*/)
{
  amrex::Abort("egtransetLENIMC not implemented for Null chemistry.");
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetLENRMC(int* /*LENRMC*/)
{
  amrex::Abort("egtransetLENRMC not implemented for Null chemistry.");
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetNO(int* /*NO*/)
{
  amrex::Abort("egtransetNO not implemented for Null chemistry.");
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetKK(int* /*KK*/)
{
  amrex::Abort("egtransetKK not implemented for Null chemistry.");
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetNLITE(int* /*NLITE*/)
{
  amrex::Abort("egtransetNLITE not implemented for Null chemistry.");
}

// Patm in ergs/cm3
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetPATM(amrex::Real* PATM)
{
  *PATM = 0.1013250000000000E+07;
}

/*the molecular weights in g/mol */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetWT(amrex::Real* /*WT*/)
{
  amrex::Abort("egtransetWT not implemented for Null chemistry.");
}

/*the lennard-jones potential well depth eps/kb in K */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetEPS(amrex::Real* /*EPS*/)
{
  amrex::Abort("egtransetEPS not implemented for Null chemistry.");
}

/*the lennard-jones collision diameter in Angstroms */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetSIG(amrex::Real* /*SIG*/)
{
  amrex::Abort("egtransetSIG not implemented for Null chemistry.");
}

/*the dipole moment in Debye */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetDIP(amrex::Real* /*DIP*/)
{
  amrex::Abort("egtransetDIP not implemented for Null chemistry.");
}

/*the polarizability in cubic Angstroms */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetPOL(amrex::Real* /*POL*/)
{
  amrex::Abort("egtransetPOL not implemented for Null chemistry.");
}

/*the rotational relaxation collision number at 298 K */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetZROT(amrex::Real* /*ZROT*/)
{
  amrex::Abort("egtransetZROT not implemented for Null chemistry.");
}

/*0: monoatomic, 1: linear, 2: nonlinear */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetNLIN(int* /*NLIN*/)
{
  amrex::Abort("egtransetNLIN not implemented for Null chemistry.");
}

/*Poly fits for the viscosities, dim NO*KK */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetCOFETA(amrex::Real* /*COFETA*/)
{
  amrex::Abort("egtransetCOFETA not implemented for Null chemistry.");
}

/*Poly fits for the conductivities, dim NO*KK */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetCOFLAM(amrex::Real* /*COFLAM*/)
{
  amrex::Abort("egtransetCOFLAM not implemented for Null chemistry.");
}

/*Poly fits for the diffusion coefficients, dim NO*KK*KK */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetCOFD(amrex::Real* /*COFD*/)
{
  amrex::Abort("egtransetCOFD not implemented for Null chemistry.");
}

// List of specs with small weight, dim NLITE
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetKTDIF(int* /*KTDIF*/)
{
  amrex::Abort("egtransetKTDIF not implemented for Null chemistry.");
}

// Poly fits for thermal diff ratios, dim NO*NLITE*KK
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
egtransetCOFTD(amrex::Real* /*COFTD*/)
{
  amrex::Abort("egtransetCOFTD not implemented for Null chemistry.");
}

/* gauss-jordan solver external routine */
/* Replace this routine with the one generated by the Gauss Jordan solver of DW
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
sgjsolve(amrex::Real* /*A*/, amrex::Real* /*x*/, amrex::Real* /*b*/)
{
  amrex::Abort("sgjsolve not implemented, choose a different solver ");
}

/* Replace this routine with the one generated by the Gauss Jordan solver of DW
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
sgjsolve_simplified(amrex::Real* /*A*/, amrex::Real* /*x*/, amrex::Real* /*b*/)
{
  amrex::Abort(
    "sgjsolve_simplified not implemented, choose a different solver ");
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
CKINDX(int& mm, int& kk, int& ii, int& nfit)
{
  mm = 1;
  kk = 1;
  ii = 0;
  nfit = -1; /*Why do you need this anyway ?  */
}

#endif

/* End of file  */
