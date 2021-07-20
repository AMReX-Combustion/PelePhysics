#ifndef _REACTOR_H_
#define _REACTOR_H_

//#include <AMREX_misc.H>

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

#include "PelePhysics.H"
#include "mechanism.H"

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

#ifdef AMREX_USE_GPU

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
box_flatten(
  int icell,
  int i,
  int j,
  int k,
  int ireactor_type,
  amrex::Array4<amrex::Real> const& rhoY,
  amrex::Array4<amrex::Real> const& frcExt,
  amrex::Array4<amrex::Real> const& temperature,
  amrex::Array4<amrex::Real> const& rhoE,
  amrex::Array4<amrex::Real> const& frcEExt,
  amrex::Real* tmp_vect,
  amrex::Real* tmp_src_vect,
  amrex::Real* tmp_vect_energy,
  amrex::Real* tmp_src_vect_energy) noexcept
{
  amrex::Real mass_frac[NUM_SPECIES];
  amrex::Real rho = 0.0;
  amrex::Real Enrg_loc;
  amrex::Real temp;

  int offset_yvec = icell * (NUM_SPECIES + 1);
  int offset_spec = icell * (NUM_SPECIES);

  for (int n = 0; n < NUM_SPECIES; n++) {
    tmp_vect[offset_yvec + n] = rhoY(i, j, k, n);
    tmp_src_vect[offset_spec + n] = frcExt(i, j, k, n);
    rho += tmp_vect[offset_yvec + n];
  }
  amrex::Real rho_inv = 1.0 / rho;
  temp = temperature(i, j, k, 0);
  tmp_vect_energy[icell] = rhoE(i, j, k, 0);
  tmp_src_vect_energy[icell] = frcEExt(i, j, k, 0);

  for (int n = 0; n < NUM_SPECIES; n++) {
    mass_frac[n] = tmp_vect[offset_yvec + n] * rho_inv;
  }
  Enrg_loc = tmp_vect_energy[icell] / rho;
  auto eos = pele::physics::PhysicsType::eos();
  if (ireactor_type == 1) {
    eos.EY2T(Enrg_loc, mass_frac, temp);
  } else {
    eos.HY2T(Enrg_loc, mass_frac, temp);
  }
  tmp_vect[offset_yvec + NUM_SPECIES] = temp;
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
box_unflatten(
  int icell,
  int i,
  int j,
  int k,
  int ireactor_type,
  amrex::Array4<amrex::Real> const& rhoY,
  amrex::Array4<amrex::Real> const& T_in,
  amrex::Array4<amrex::Real> const& rhoE,
  amrex::Array4<amrex::Real> const& frcEExt,
  amrex::Array4<amrex::Real> const& FC_in,
  amrex::Real* tmp_vect,
  amrex::Real* tmp_vect_energy,
  long int const& tmp_FCunt,
  amrex::Real dt) noexcept
{
  amrex::Real mass_frac[NUM_SPECIES];
  amrex::Real rho = 0.0;
  amrex::Real Enrg_loc;
  amrex::Real temp;

  int offset_yvec = icell * (NUM_SPECIES + 1);

  for (int n = 0; n < NUM_SPECIES; n++) {
    rhoY(i, j, k, n) = tmp_vect[offset_yvec + n];
    rho += tmp_vect[offset_yvec + n];
  }
  amrex::Real rho_inv = 1.0 / rho;
  temp = tmp_vect[offset_yvec + NUM_SPECIES];
  rhoE(i, j, k, 0) = tmp_vect_energy[icell] + dt * frcEExt(i, j, k, 0);

  for (int n = 0; n < NUM_SPECIES; n++) {
    mass_frac[n] = rhoY(i, j, k, n) * rho_inv;
  }
  Enrg_loc = rhoE(i, j, k, 0) / rho;
  auto eos = pele::physics::PhysicsType::eos();
  if (ireactor_type == 1) {
    eos.EY2T(Enrg_loc, mass_frac, temp);
  } else {
    eos.HY2T(Enrg_loc, mass_frac, temp);
  }
  T_in(i, j, k, 0) = temp;
  FC_in(i, j, k, 0) = tmp_FCunt;
}

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

void SetTypValsODE(const std::vector<amrex::Real>& ExtTypVals);

void SetTolFactODE(amrex::Real relative_tol, amrex::Real absolute_tol);

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

#else

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
box_flatten(
  int icell,
  int i,
  int j,
  int k,
  int ireactor_type,
  amrex::Array4<amrex::Real> const& rhoY,
  amrex::Array4<amrex::Real> const& frcExt,
  amrex::Array4<amrex::Real> const& temperature,
  amrex::Array4<amrex::Real> const& rhoE,
  amrex::Array4<amrex::Real> const& frcEExt,
  amrex::Real* tmp_vect,
  amrex::Real* tmp_src_vect,
  amrex::Real* tmp_vect_energy,
  amrex::Real* tmp_src_vect_energy) noexcept
{
  amrex::Real mass_frac[NUM_SPECIES];
  amrex::Real rho = 0.0;
  amrex::Real Enrg_loc;
  amrex::Real temp;

  int offset_yvec = icell * (NUM_SPECIES + 1);
  int offset_spec = icell * (NUM_SPECIES);

  for (int n = 0; n < NUM_SPECIES; n++) {
    tmp_vect[offset_yvec + n] = rhoY(i, j, k, n);
    tmp_src_vect[offset_spec + n] = frcExt(i, j, k, n);
    rho += tmp_vect[offset_yvec + n];
  }
  amrex::Real rho_inv = 1.0 / rho;
  temp = temperature(i, j, k);
  tmp_vect_energy[icell] = rhoE(i, j, k);
  tmp_src_vect_energy[icell] = frcEExt(i, j, k);

  for (int n = 0; n < NUM_SPECIES; n++) {
    mass_frac[n] = tmp_vect[offset_yvec + n] * rho_inv;
  }
  Enrg_loc = tmp_vect_energy[icell] / rho;
  auto eos = pele::physics::PhysicsType::eos();
  if (ireactor_type == 1) {
    eos.EY2T(Enrg_loc, mass_frac, temp);
  } else {
    eos.HY2T(Enrg_loc, mass_frac, temp);
  }
  tmp_vect[offset_yvec + NUM_SPECIES] = temp;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
box_unflatten(
  int icell,
  int i,
  int j,
  int k,
  int ireactor_type,
  amrex::Array4<amrex::Real> const& rhoY,
  amrex::Array4<amrex::Real> const& T_in,
  amrex::Array4<amrex::Real> const& rhoE,
  amrex::Array4<amrex::Real> const& frcEExt,
  amrex::Array4<amrex::Real> const& FC_in,
  amrex::Real* tmp_vect,
  amrex::Real* tmp_vect_energy,
  long int const& tmp_FCunt,
  amrex::Real dt) noexcept
{
  amrex::Real mass_frac[NUM_SPECIES];
  amrex::Real rho = 0.0;
  amrex::Real Enrg_loc;
  amrex::Real temp;

  int offset_yvec = icell * (NUM_SPECIES + 1);

  for (int n = 0; n < NUM_SPECIES; n++) {
    rhoY(i, j, k, n) = tmp_vect[offset_yvec + n];
    rho += tmp_vect[offset_yvec + n];
  }
  amrex::Real rho_inv = 1.0 / rho;
  temp = tmp_vect[offset_yvec + NUM_SPECIES];
  rhoE(i, j, k) = tmp_vect_energy[icell] + dt * frcEExt(i, j, k, 0);

  for (int n = 0; n < NUM_SPECIES; n++) {
    mass_frac[n] = rhoY(i, j, k, n) * rho_inv;
  }
  Enrg_loc = rhoE(i, j, k) / rho;
  auto eos = pele::physics::PhysicsType::eos();
  if (ireactor_type == 1) {
    eos.EY2T(Enrg_loc, mass_frac, temp);
  } else {
    eos.HY2T(Enrg_loc, mass_frac, temp);
  }
  T_in(i, j, k) = temp;

  FC_in(i, j, k) = tmp_FCunt;
}

typedef struct
{
  bool FirstTimePrecond;
  bool reactor_cvode_initialized;
  bool actual_ok_to_react;
  int ncells;
  int iverbose;
  int isolve_type;
  int ianalytical_jacobian;
  int ireactor_type;
  int boxcell;
  amrex::Real* Yvect_full = NULL;
  amrex::Real* rhoX_init = NULL;
  amrex::Real* rhoXsrc_ext = NULL;
  amrex::Real* rYsrc = NULL;
  int* FCunt = NULL;
  int* mask = NULL;

  int NNZ;
  SUNMatrix* PS;
  realtype** Jdata = NULL;
  int** rowVals = NULL;
  int** rowPtrs = NULL;
  int** colPtrs = NULL;
  int** colVals = NULL;
  int* indx = NULL;
  realtype** JSPSmat = NULL;
#ifdef USE_KLU_PP
  klu_common* Common;
  klu_symbolic** Symbolic;
  klu_numeric** Numeric;
#else
  realtype**(**Jbd);
  realtype**(**P);
  sunindextype*(**pivot);
#endif
  SUNMatrix PSc;
  int* colVals_c;
  int* rowPtrs_c;
} * UserData;

int cF_RHS(realtype t, N_Vector y_in, N_Vector ydot, void* user_data);

int cJac(
  realtype tn,
  N_Vector y,
  N_Vector fy,
  SUNMatrix J,
  void* user_data,
  N_Vector tmp1,
  N_Vector tmp2,
  N_Vector tmp3);

int cJac_sps(
  realtype tn,
  N_Vector y,
  N_Vector fy,
  SUNMatrix J,
  void* user_data,
  N_Vector tmp1,
  N_Vector tmp2,
  N_Vector tmp3);

int PSolve_custom(
  realtype tn,
  N_Vector u,
  N_Vector fu,
  N_Vector r,
  N_Vector z,
  realtype gamma,
  realtype delta,
  int lr,
  void* user_data);

int Precond_custom(
  realtype tn,
  N_Vector u,
  N_Vector fu,
  booleantype jok,
  booleantype* jcurPtr,
  realtype gamma,
  void* user_data);

#ifdef USE_KLU_PP
int cJac_KLU(
  realtype tn,
  N_Vector y,
  N_Vector fy,
  SUNMatrix J,
  void* user_data,
  N_Vector tmp1,
  N_Vector tmp2,
  N_Vector tmp3);

int PSolve_sparse(
  realtype tn,
  N_Vector u,
  N_Vector fu,
  N_Vector r,
  N_Vector z,
  realtype gamma,
  realtype delta,
  int lr,
  void* user_data);
int Precond_sparse(
  realtype tn,
  N_Vector u,
  N_Vector fu,
  booleantype jok,
  booleantype* jcurPtr,
  realtype gamma,
  void* user_data);
#else
int PSolve(
  realtype tn,
  N_Vector u,
  N_Vector fu,
  N_Vector r,
  N_Vector z,
  realtype gamma,
  realtype delta,
  int lr,
  void* user_data);

int Precond(
  realtype tn,
  N_Vector u,
  N_Vector fu,
  booleantype jok,
  booleantype* jcurPtr,
  realtype gamma,
  void* user_data);
#endif

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
  const int& reactor_type);

int react_2(
  const amrex::Box& box,
  amrex::Array4<amrex::Real> const& rY_in,
  amrex::Array4<amrex::Real> const& rY_src_in,
  amrex::Array4<amrex::Real> const& T_in,
  amrex::Array4<amrex::Real> const& rEner_in,
  amrex::Array4<amrex::Real> const& rEner_src_in,
  amrex::Array4<amrex::Real> const& FC_in,
  amrex::Array4<int> const& mask,
  amrex::Real& dt_react,
  amrex::Real& time);

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

void reactor_close();

int check_flag(void* flagvalue, const char* funcname, int opt);

void PrintFinalStats(void* cvodeMem, amrex::Real Temp);

UserData AllocUserData(int iE, int num_cells);

void FreeUserData(UserData data);

void check_state(N_Vector yvec);

void SetTypValsODE(const std::vector<double>& ExtTypVals);

void SetTolFactODE(double relative_tol, double absolute_tol);

void ReSetTolODE();

void
fKernelSpec(realtype* dt, realtype* yvec_d, realtype* ydot_d, void* user_data);

struct _SUNLinearSolverContent_Sparse_custom
{
  sunindextype last_flag;
  int reactor_type;
  int nsubsys;     // number of subsystems
  int subsys_size; // size of each subsystem
  int subsys_nnz;
};

typedef struct _SUNLinearSolverContent_Sparse_custom*
  SUNLinearSolverContent_Sparse_custom;

SUNLinearSolver SUNLinSol_sparse_custom(
  N_Vector y,
  SUNMatrix A,
  int reactor_type,
  int nsubsys,
  int subsys_size,
  int subsys_nnz);

SUNLinearSolver_Type SUNLinSolGetType_Sparse_custom(SUNLinearSolver S);

int SUNLinSolSolve_Sparse_custom(
  SUNLinearSolver S, SUNMatrix A, N_Vector x, N_Vector b, realtype tol);

#endif

#endif
