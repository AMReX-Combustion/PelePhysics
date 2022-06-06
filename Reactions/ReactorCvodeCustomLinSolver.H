#ifndef REACTORCVODECUSTOMLINSOLVER_H
#define REACTORCVODECUSTOMLINSOLVER_H

#include <AMReX_BLProfiler.H>
#include <AMReX_REAL.H>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_sparse.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#ifdef AMREX_USE_CUDA
#include <nvector/nvector_cuda.h>
#include <sunmatrix/sunmatrix_cusparse.h>
#include <sunlinsol/sunlinsol_cusolversp_batchqr.h>
#endif
#include "mechanism.H"

namespace pele {
namespace physics {
namespace reactions {
namespace cvode {
#ifdef AMREX_USE_GPU
struct _SUNLinearSolverContent_Dense_custom
{
  sunindextype last_flag;
  int nsubsys;
  int subsys_size;
  int subsys_nnz;
  int nbBlocks;
  int nbThreads;
  amrex::gpuStream_t stream;
};

using SUNLinearSolverContent_Dense_custom =
  _SUNLinearSolverContent_Dense_custom*;

SUNLinearSolver SUNLinSol_dense_custom(
  N_Vector y, SUNMatrix A, amrex::gpuStream_t stream, SUNContext sunctx);

SUNLinearSolver_Type SUNLinSolGetType_Dense_custom(SUNLinearSolver S);

int SUNLinSolSolve_Dense_custom(
  SUNLinearSolver S, SUNMatrix A, N_Vector x, N_Vector b, amrex::Real tol);

int SUNLinSolSetup_Dense_custom(SUNLinearSolver S, SUNMatrix A);

int SUNLinSolFree_Dense_custom(SUNLinearSolver S);

#else
struct SUNLinearSolverContent_Sparse_custom_struct
{
  sunindextype last_flag;
  int reactor_type;
  int nsubsys;     // number of subsystems
  int subsys_size; // size of each subsystem
  int subsys_nnz;
};

using SUNLinearSolverContent_Sparse_custom =
  SUNLinearSolverContent_Sparse_custom_struct*;

SUNLinearSolver SUNLinSol_sparse_custom(
  N_Vector y,
  SUNMatrix A,
  int reactor_type,
  int nsubsys,
  int subsys_size,
  int subsys_nnz,
  SUNContext sunctx);

SUNLinearSolver_Type SUNLinSolGetType_Sparse_custom(SUNLinearSolver S);

int SUNLinSolSolve_Sparse_custom(
  SUNLinearSolver S, SUNMatrix A, N_Vector x, N_Vector b, amrex::Real tol);
#endif
} // namespace cvode
} // namespace reactions
} // namespace physics
} // namespace pele
#endif
