#ifndef _REACTOR_H_
#define _REACTOR_H_

//#include <AMREX_misc.H>

#include <AMReX_Gpu.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>

#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_erkstep.h>
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
#include "reactor_utils.H"

extern amrex::Array<double, NUM_SPECIES + 1> typVals;
extern int eint_rho; // in/out = rhoE/rhoY
extern int enth_rho; // in/out = rhoH/rhoY
extern amrex::Real relTol;
extern amrex::Real absTol;
extern amrex::Real time_init;
extern int dense_solve;
extern int sparse_solve;
extern int iterative_gmres_solve;
extern int sparse_solve_custom;
extern int iterative_gmres_solve_custom;
extern int hack_dump_sparsity_pattern;
extern int sparse_cusolver_solve;
extern int sparse_custom_solve;
extern int use_erkstep;
extern int rk_method;
extern int rk_controller;

#ifdef AMREX_USE_OMP
#pragma omp threadprivate(time_init)
#pragma omp threadprivate(typVals)
#pragma omp threadprivate(relTol, absTol)
#endif

//===================================================================
//initialize and close reactor functions
//===================================================================
int reactor_init(int reactor_type, int Ncells);
void reactor_close();
//===================================================================


//===================================================================
//User data to feed into CPU/GPU sundials handle
//===================================================================
#ifdef AMREX_USE_GPU
 struct CVODEUserData
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
  amrex::gpuStream_t stream;
  int nbBlocks;
  int nbThreads;
};
#else
 struct CVODEUserData
{
  amrex::Real dt_save;
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
};
#endif

 struct ARKODEUserData
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
};

//===================================================================
//Box flattening functions
//===================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void box_flatten(
  int icell,
  int ncells,
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
  amrex::Real* tmp_src_vect_energy)
{
  amrex::Real mass_frac[NUM_SPECIES];
  amrex::Real rho = 0.0;
  for (int n = 0; n < NUM_SPECIES; n++) {
    const int idx = sunvec_index(n, icell, ncells);
    const int idx2 = sunspec_index(n, icell, ncells);
    const amrex::Real rhoYn = rhoY(i, j, k, n);
    tmp_vect[idx] = rhoYn;
    tmp_src_vect[idx2] = frcExt(i, j, k, n);
    rho += rhoYn;
    mass_frac[n] = rhoYn;
  }
  const amrex::Real rho_inv = 1.0 / rho;
  for (int n = 0; n < NUM_SPECIES; n++) {
    mass_frac[n] *= rho_inv;
  }

  amrex::Real temp = temperature(i, j, k, 0);
  tmp_vect_energy[icell] = rhoE(i, j, k, 0);
  tmp_src_vect_energy[icell] = frcEExt(i, j, k, 0);

  amrex::Real Enrg_loc = tmp_vect_energy[icell] * rho_inv;
  auto eos = pele::physics::PhysicsType::eos();
  if (ireactor_type == 1) {
    eos.EY2T(Enrg_loc, mass_frac, temp);
  } else {
    eos.HY2T(Enrg_loc, mass_frac, temp);
  }
  tmp_vect[sunvec_index(NUM_SPECIES, icell, ncells)] = temp;
}
//===================================================================
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
box_unflatten(
  int icell,
  int ncells,
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
  amrex::Real dt) 
{
  amrex::Real mass_frac[NUM_SPECIES];
  amrex::Real rho = 0.0;
  for (int n = 0; n < NUM_SPECIES; n++) {
    const amrex::Real rhoYn = tmp_vect[sunvec_index(n, icell, ncells)];
    rhoY(i, j, k, n) = rhoYn;
    rho += rhoYn;
    mass_frac[n] = rhoYn;
  }
  const amrex::Real rho_inv = 1.0 / rho;
  for (int n = 0; n < NUM_SPECIES; n++) {
    mass_frac[n] *= rho_inv;
  }

  amrex::Real temp = tmp_vect[sunvec_index(NUM_SPECIES, icell, ncells)];
  rhoE(i, j, k, 0) = tmp_vect_energy[icell] + dt * frcEExt(i, j, k, 0);

  amrex::Real Enrg_loc = rhoE(i, j, k, 0) * rho_inv;
  auto eos = pele::physics::PhysicsType::eos();
  if (ireactor_type == 1) {
    eos.EY2T(Enrg_loc, mass_frac, temp);
  } else {
    eos.HY2T(Enrg_loc, mass_frac, temp);
  }
  T_in(i, j, k, 0) = temp;
  FC_in(i, j, k, 0) = tmp_FCunt;
}
//===================================================================

//===================================================================
//Tolerance setting with typical values
//===================================================================
void SetTypValsODE(const std::vector<amrex::Real>& ExtTypVals);
void SetTolFactODE(amrex::Real relative_tol, amrex::Real absolute_tol);


//===================================================================
//RHS computation functions
//===================================================================
static int cF_RHS(realtype t, N_Vector y_in, N_Vector ydot, void* user_data);

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void fKernelSpec(int icell, int ncells,
  double dt_save,
  int reactor_type,
  realtype* yvec_d,
  realtype* ydot_d,
  amrex::Real* rhoe_init,
  amrex::Real* rhoesrc_ext,
  amrex::Real* rYs)
{
  amrex::Real rho_pt = 0.0;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {0.0};
  for (int n = 0; n < NUM_SPECIES; n++) {
    massfrac[n] = yvec_d[sunvec_index(n, icell, ncells)];
    rho_pt += massfrac[n];
  }
  const amrex::Real rho_pt_inv = 1.0 / rho_pt;

  for (int n = 0; n < NUM_SPECIES; n++) {
    massfrac[n] *= rho_pt_inv;
  }

  const amrex::Real nrg_pt =
    (rhoe_init[icell] + rhoesrc_ext[icell] * dt_save) * rho_pt_inv;

  amrex::Real temp_pt = yvec_d[sunvec_index(NUM_SPECIES, icell, ncells)];

  amrex::Real Cv_pt = 0.0;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> ei_pt = {0.0};
  auto eos = pele::physics::PhysicsType::eos();
  if (reactor_type == 1) {
    eos.EY2T(nrg_pt, massfrac.arr, temp_pt);
    eos.T2Ei(temp_pt, ei_pt.arr);
    eos.TY2Cv(temp_pt, massfrac.arr, Cv_pt);
  } else {
    eos.HY2T(nrg_pt, massfrac.arr, temp_pt);
    eos.TY2Cp(temp_pt, massfrac.arr, Cv_pt);
    eos.T2Hi(temp_pt, ei_pt.arr);
  }

  amrex::GpuArray<amrex::Real, NUM_SPECIES> cdots_pt = {0.0};
  eos.RTY2WDOT(rho_pt, temp_pt, massfrac.arr, cdots_pt.arr);

  amrex::Real rhoesrc = rhoesrc_ext[icell];
  for (int n = 0; n < NUM_SPECIES; n++) {
    const amrex::Real cdot_rYs = cdots_pt[n] + rYs[sunspec_index(n, icell, ncells)];
    ydot_d[sunvec_index(n, icell, ncells)] = cdot_rYs;
    rhoesrc -= cdot_rYs * ei_pt[n];
  }
  ydot_d[sunvec_index(NUM_SPECIES, icell, ncells)] = rhoesrc * (rho_pt_inv / Cv_pt);
}
//===================================================================

//===================================================================
//Different variants of react calls
//===================================================================
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
//===================================================================

//===================================================================
//CUDA specific things for cvode
//===================================================================
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
//===================================================================

//===================================================================
//Stats and checks
//===================================================================
int check_flag(void* flagvalue, const char* funcname, int opt);

//FIXME: make them the same
#ifdef AMREX_USE_GPU
static void PrintFinalStats(void* cvode_mem);
#else
static void PrintFinalStats(void* cvode_mem,amrex::Real Temp);
#endif
//===================================================================


//===================================================================
//CPU specific cvode functions
//===================================================================
#ifndef AMREX_USE_GPU
CVODEUserData* AllocUserData(int reactor_type, int num_cells);

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

void FreeUserData(CVODEUserData *data);

void check_state(N_Vector yvec);

void SetTypValsODE(const std::vector<double>& ExtTypVals);

void ReSetTolODE();

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
//===================================================================

#endif
