#ifndef _EOS_H_
#define _EOS_H_

#include <AMReX_ArrayLim.H>
#include <AMReX_BLFort.H>

extern "C" {

    void eos_re_ext(amrex::Real *rho, amrex::Real *Y, amrex::Real *T, amrex::Real *e, amrex::Real *ei, amrex::Real *cv);

}
#endif 
