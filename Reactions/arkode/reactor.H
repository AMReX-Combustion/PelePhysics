#ifndef _REACTOR_H_
#define _REACTOR_H_

#include <AMReX_Gpu.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_erkstep.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

#include "mechanism.H"
#include "PelePhysics.H"
#include "reactor_utilities.H"

#ifdef AMREX_USE_GPU
#include "AMReX_SUNMemory.H"
#endif
#if defined(AMREX_USE_CUDA)
#include <nvector/nvector_cuda.h>
#elif defined(AMREX_USE_HIP)
#include <nvector/nvector_hip.h>
#endif

typedef struct ARKODEUserData
{
  bool reactor_arkode_initialized;
  int ncells_d;
  int neqs_per_cell;
  int iverbose;
  int ireactor_type;
  amrex::Real dt_save;
  amrex::Real* rhoe_init_d = NULL;
  amrex::Real* rhoesrc_ext_d = NULL;
  amrex::Real* rYsrc_d = NULL;
#ifdef AMREX_USE_GPU
  amrex::gpuStream_t stream;
#endif
  int nbBlocks;
  int nbThreads;
} * UserData;

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void fKernelSpec(
  int ncells,
  amrex::Real dt_save,
  int reactor_type,
  realtype* yvec_d,
  realtype* ydot_d,
  amrex::Real* rhoX_init,
  amrex::Real* rhoXsrc_ext,
  amrex::Real* rYs);

int reactor_init(int cvode_iE, int Ncells);

int cF_RHS(realtype t, N_Vector y_in, N_Vector ydot, void* user_data);

int check_flag(void* flagvalue, const char* funcname, int opt);

// React for array4
int react(
  const amrex::Box& box,
  amrex::Array4<amrex::Real> const& rY_in,
  amrex::Array4<amrex::Real> const& rY_src_in,
  amrex::Array4<amrex::Real> const& T_in,
  amrex::Array4<amrex::Real> const& rEner_in,
  amrex::Array4<amrex::Real> const& rEner_src_in,
  amrex::Array4<amrex::Real> const& FC_in,
  amrex::Array4<int> const& mask,
  amrex::Real& dt_react,
  amrex::Real& time,
  const int& reactor_type
#ifdef AMREX_USE_GPU
  ,
  amrex::gpuStream_t stream
#endif
);

// React for 1d array
int react(
  realtype* rY_in,
  realtype* rY_src_in,
  realtype* rX_in,
  realtype* rX_src_in,
  realtype& dt_react,
  realtype& time,
  int reactor_type,
  int Ncells
#ifdef AMREX_USE_GPU
  ,
  amrex::gpuStream_t stream
#endif
);

#endif
