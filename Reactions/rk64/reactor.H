#ifndef _REACTOR_H_
#define _REACTOR_H_

#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>

#include "mechanism.H"
#include "PelePhysics.H"

typedef struct
{
  int ncells;
  int iverbose;
  int ireactor_type;
  amrex::Real errtol;
  int nsubsteps_guess;
  int nsubsteps_min;
  int nsubsteps_max;
} * UserData;

int reactor_init(
  int reactor_type,
  int Ncells,
  amrex::Real rk64_errtol = 1e-16,
  int rk64_nsubsteps_guess = 10,
  int rk64_nsusbteps_min = 5,
  int rk64_nsubsteps_max = 500);

int react(
  amrex::Real* rY_in,
  amrex::Real* rY_src_in,
  amrex::Real* rX_in,
  amrex::Real* rX_src_in,
  amrex::Real& dt_react,
  amrex::Real& time);

void reactor_close();

void FreeUserData(UserData data);

void fKernelSpec(
  amrex::Real* dt, amrex::Real* yvec_d, amrex::Real* ydot_d, void* user_data);

UserData AllocUserData(
  int reactor_type,
  int num_cells,
  amrex::Real rk64_errtol,
  int rk64_nsubsteps_guess,
  int rk64_nsubsteps_min,
  int rk64_nsubsteps_max);

void SetTypValsODE(const std::vector<amrex::Real>& ExtTypVals);

void SetTolFactODE(amrex::Real relative_tol, amrex::Real absolute_tol);

#endif
