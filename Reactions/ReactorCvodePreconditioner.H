#ifndef _REACTORCVODEPRECONDITIONER_H_
#define _REACTORCVODEPRECONDITIONER_H_

#include <AMReX_REAL.H>
#include <nvector/nvector_serial.h>
#include <sunlinsol/sunlinsol_dense.h>
#include "PelePhysics.H"
#include "mechanism.H"
#include "ReactorCvodeUtils.H"

#ifdef AMREX_USE_HIP
#include <nvector/nvector_hip.h>
#endif

#ifdef AMREX_USE_CUDA
#include <nvector/nvector_cuda.h>
#endif

namespace pele {
namespace physics {
namespace reactions {
namespace cvode {
// Basic preconditioner functions
int Precond(
  amrex::Real tn,
  N_Vector u,
  N_Vector fu,
  booleantype jok,
  booleantype* jcurPtr,
  amrex::Real gamma,
  void* user_data);

int PSolve(
  amrex::Real tn,
  N_Vector u,
  N_Vector fu,
  N_Vector r,
  N_Vector z,
  amrex::Real gamma,
  amrex::Real delta,
  int lr,
  void* user_data);

#ifndef AMREX_USE_GPU
// Custom SGS preconditioner functions
int PSolve_custom(
  amrex::Real tn,
  N_Vector u,
  N_Vector fu,
  N_Vector r,
  N_Vector z,
  amrex::Real gamma,
  amrex::Real delta,
  int lr,
  void* user_data);

int Precond_custom(
  amrex::Real tn,
  N_Vector u,
  N_Vector fu,
  booleantype jok,
  booleantype* jcurPtr,
  amrex::Real gamma,
  void* user_data);

#ifdef PELE_USE_KLU
// KLU preconditioner functions
int PSolve_sparse(
  amrex::Real tn,
  N_Vector u,
  N_Vector fu,
  N_Vector r,
  N_Vector z,
  amrex::Real gamma,
  amrex::Real delta,
  int lr,
  void* user_data);

int Precond_sparse(
  amrex::Real tn,
  N_Vector u,
  N_Vector fu,
  booleantype jok,
  booleantype* jcurPtr,
  amrex::Real gamma,
  void* user_data);
#endif
#endif
} // namespace cvode
} // namespace reactions
} // namespace physics
} // namespace pele
#endif
