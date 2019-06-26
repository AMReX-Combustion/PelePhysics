#ifndef _eos_type_cpp_H
#define _eos_type_cpp_H

#include <AMReX_REAL.H>   
#include <AMReX_IntVect.H> 

#ifndef BL_LANG_FORT

struct  EOSTYPE
{
	amrex::Real  rho
	amrex::Real  T
	amrex::Real  p
	amrex::Real  e
	amrex::Real  h
	amrex::Real  s
        amrex::Real  f

        amrex::Real  *massfrac
        amrex::Real  *molefrac
        amrex::Real  *aux

        amrex::Real  dpdT
        amrex::Real  dpdr
        amrex::Real  dedT
        amrex::Real  dedr
        amrex::Real  dhdT
        amrex::Real  dhdr
        amrex::Real  dsdT
        amrex::Real  dsdr
        amrex::Real  dpde
        amrex::Real  dpdr_e

        amrex::Real  cv
        amrex::Real  cp
        amrex::Real  *cpi
        amrex::Real  *cvi
        amrex::Real  *hi
        amrex::Real  *ei
        amrex::Real  *si
        amrex::Real  wbar
        amrex::Real  *mui
        amrex::Real  *Acti
        amrex::Real  *dedY
        amrex::Real  *dpdY
        amrex::Real  *dhdY
        amrex::Real  gam1
        amrex::Real  cs

        /* Quantities used for non-Ideal EOS */
        amrex::Real  am
        amrex::Real  bm
        amrex::Real  *damdYk
        amrex::Real  *d2amdYkdT
        amrex::Real  *dPdYk
        amrex::Real  damdT
        amrex::Real  d2amdT2
        amrex::Real  dpdtau
        amrex::Real  Z
        amrex::Real  *taui
        amrex::Real  *diP
        amrex::Real  *dijY
};

//extern "C" 
//{
//  void eos_build_cpp(EOSTYPE* EOST) 
//  {
//	 = (double *)  malloc((NCELLS*NEQ)*sizeof(double));
//  };

#endif

#endif
