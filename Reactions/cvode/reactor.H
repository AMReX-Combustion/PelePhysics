#ifndef _REACTOR_H_
#define _REACTOR_H_

#include <AMReX_Gpu.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>

#include <cvode/cvode.h>
#include <cvode/cvode_direct.h>
#include <cvode/cvode_spils.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>
#include <sunmatrix/sunmatrix_sparse.h>

#ifdef USE_KLU_PP
#include "klu.h"
#include <sunlinsol/sunlinsol_klu.h>
#endif

#include "reactor_utilities.H"
#include "PelePhysics.H"
#include "chemistry_file.H"
#include "mechanism.h"

#ifdef AMREX_USE_CUDA
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cusolverSp.h>
#include <cusparse.h>
#include <cuda_runtime_api.h>
#include <nvector/nvector_cuda.h>
#include <sunmatrix/sunmatrix_cusparse.h>
#include <sunlinsol/sunlinsol_cusolversp_batchqr.h>
#endif

#ifdef AMREX_USE_HIP
#include <nvector/nvector_hip.h>
#endif

#include "AMReX_SUNMemory.H"

typedef struct CVodeUserData
{
  amrex::Real gamma;
  amrex::Real dt_save;
  int ncells;
  int neqs_per_cell;
  int ireactor_type;
  int ianalytical_jacobian;
  int isolve_type;
  amrex::Real* rhoe_init_d = NULL;
  amrex::Real* rhoesrc_ext_d = NULL;
  amrex::Real* rYsrc_d = NULL;
  int iverbose;
  int NNZ;
  int* csr_row_count_h;
  int* csr_col_index_h;
  int* csr_row_count_d;
  int* csr_col_index_d;
  amrex::Real* csr_val_h;
  amrex::Real* csr_jac_h;
  amrex::Real* csr_val_d;
  amrex::Real* csr_jac_d;
  SUNMatrix R = NULL;
#ifdef AMREX_USE_CUDA
  void* buffer_qr = NULL;
  csrqrInfo_t info;
  cusparseMatDescr_t descrA;
  cusolverSpHandle_t cusolverHandle;
  cusparseHandle_t cuSPHandle;
#endif
#ifdef AMREX_USE_GPU
  amrex::gpuStream_t stream;
#endif
  int nbBlocks;
  int nbThreads;
} * UserData;

static int cF_RHS(realtype t, N_Vector y_in, N_Vector ydot, void* user_data);

int reactor_init(int cvode_iE, int Ncells);

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

#ifdef AMREX_USE_CUDA
static int Precond(
  realtype tn,
  N_Vector u,
  N_Vector fu,
  booleantype jok,
  booleantype* jcurPtr,
  realtype gamma,
  void* user_data);

static int PSolve(
  realtype tn,
  N_Vector u,
  N_Vector fu,
  N_Vector r,
  N_Vector z,
  realtype gamma,
  realtype delta,
  int lr,
  void* user_data);

static int cJac(
  realtype t,
  N_Vector y,
  N_Vector fy,
  SUNMatrix J,
  void* user_data,
  N_Vector tmp1,
  N_Vector tmp2,
  N_Vector tmp3);
#endif

void reactor_close();

static int check_flag(void* flagvalue, const char* funcname, int opt);

static void PrintFinalStats(void* cvode_mem);

// void SetTypValsODE(const std::vector<amrex::Real>& ExtTypVals);

// void SetTolFactODE(amrex::Real relative_tol, amrex::Real absolute_tol);

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

// All the following functions are for CUDA only implementation for now

#ifdef AMREX_USE_CUDA
// JACOBIANS
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void fKernelComputeallAJ(
  int ncells, void* user_data, realtype* u_d, realtype* csr_val);

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void fKernelComputeAJsys(
  int ncells, void* user_data, realtype* u_d, realtype* csr_val);

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void fKernelComputeAJchem(
  int ncells, void* user_data, realtype* u_d, realtype* Jdata);

// CUSTOM
__global__ void fKernelDenseSolve(
  int ncells,
  realtype* x_d,
  realtype* b_d,
  int subsys_size,
  int subsys_nnz,
  realtype* csr_val);

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

typedef struct _SUNLinearSolverContent_Dense_custom*
  SUNLinearSolverContent_Dense_custom;

SUNLinearSolver
SUNLinSol_dense_custom(N_Vector y, SUNMatrix A, amrex::gpuStream_t stream);

SUNLinearSolver_Type SUNLinSolGetType_Dense_custom(SUNLinearSolver S);

int SUNLinSolSolve_Dense_custom(
  SUNLinearSolver S, SUNMatrix A, N_Vector x, N_Vector b, realtype tol);

int SUNLinSolSetup_Dense_custom(SUNLinearSolver S, SUNMatrix A);

int SUNLinSolFree_Dense_custom(SUNLinearSolver S);
#endif

#endif
