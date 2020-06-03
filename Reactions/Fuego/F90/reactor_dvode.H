#ifndef _REACTOR_DVODE_H_
#define _REACTOR_DVODE_H_

#include <AMReX_ArrayLim.H>
#include <AMReX_BLFort.H>

extern "C" {

int reactor_init(const int* cvode_iE, const int* Ncells);

int react(amrex::Real *rY_in, amrex::Real *rY_src_in, 
		amrex::Real *rX_in, amrex::Real *rX_src_in, 
		amrex::Real *P_in, 
		amrex::Real *dt_react, amrex::Real *time);

}
#endif 
