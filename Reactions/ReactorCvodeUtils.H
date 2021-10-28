#ifndef _REACTORUTILSCVODE_H_
#define _REACTORUTILSCVODE_H_
#include "ReactorTypes.H"

#include <AMReX.H>
#include <AMReX_Print.H>
#include <AMReX_REAL.H>
#include <cvode/cvode.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_sparse.h>
#include "ReactorUtils.H"

#ifdef AMREX_USE_CUDA
#include <cusolverSp.h>
#endif

#ifdef PELE_USE_KLU
#include "klu.h"
#include <sunlinsol/sunlinsol_klu.h>
#endif

namespace pele {
namespace physics {
namespace reactions {

namespace cvode {

#ifdef AMREX_USE_GPU
enum linSolveType { sparseDirect, customDirect, magmaDirect, GMRES, precGMRES };
enum precondType { sparseSimpleAJac };
#else
enum linSolveType {
  denseFDDirect,
  denseDirect,
  sparseDirect,
  customDirect,
  GMRES,
  precGMRES,
  hackDumpSparsePattern
};
enum precondType { denseSimpleAJac, sparseSimpleAJac, customSimpleAJac };
#endif

// Error function for CVODE
void cvodeErrHandler(
  int error_code,
  const char* /*module*/,
  const char* /*function*/,
  char* msg,
  void* /*eh_data*/);

void printFinalStats(void* cvodeMem);

#ifdef AMREX_USE_GPU
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
fKernelComputeallAJ(
  int ncell, void* user_data, amrex::Real* u_d, amrex::Real* csr_val_arg)
{
  CVODEUserData* udata = static_cast<CVODEUserData*>(user_data);
  int neqs = NUM_SPECIES + 1;
  auto NNZ = udata->NNZ;
  auto reactor_type = udata->reactor_type;
  auto csr_jac_d = udata->csr_jac_d;
  auto csr_row_count_d = udata->csr_row_count_d;
  auto csr_col_index_d = udata->csr_col_index_d;
  auto gamma = udata->gamma;

  amrex::GpuArray<amrex::Real, NUM_SPECIES> activity{0.0};
  amrex::GpuArray<amrex::Real, (NUM_SPECIES + 1) * (NUM_SPECIES + 1)> Jmat_pt{
    0.0};

  int u_offset = ncell * neqs;
  amrex::Real* u_curr = u_d + u_offset;

  amrex::Real mw[NUM_SPECIES] = {0.0};
  get_mw(mw);

  amrex::Real rho_pt = 0.0;
  for (int n = 0; n < NUM_SPECIES; n++) {
    rho_pt = rho_pt + u_curr[n];
  }

  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac{0.0};
  for (int i = 0; i < NUM_SPECIES; i++) {
    massfrac[i] = u_curr[i] / rho_pt;
  }

  amrex::Real temp_pt = u_curr[NUM_SPECIES];

  auto eos = pele::physics::PhysicsType::eos();
  eos.RTY2C(rho_pt, temp_pt, massfrac.arr, activity.arr);

  int consP = reactor_type == ReactorTypes::h_reactor_type;
  DWDOT_SIMPLIFIED(Jmat_pt.arr, activity.arr, &temp_pt, &consP);

  int jac_offset = ncell * NNZ;
  amrex::Real* csr_jac_cell = csr_jac_d + jac_offset;
  amrex::Real* csr_val_cell = csr_val_arg + jac_offset;
  for (int i = 0; i < NUM_SPECIES; i++) {
    for (int k = 0; k < NUM_SPECIES; k++) {
      Jmat_pt[k * neqs + i] *= mw[i] / mw[k];
    }
    Jmat_pt[i * neqs + NUM_SPECIES] /= mw[i];
    Jmat_pt[NUM_SPECIES * neqs + i] *= mw[i];
  }
  for (int i = 1; i < NUM_SPECIES + 2; i++) {
    int nbVals = csr_row_count_d[i] - csr_row_count_d[i - 1];
    for (int j = 0; j < nbVals; j++) {
      int idx = csr_col_index_d[csr_row_count_d[i - 1] + j - 1] - 1;
      if (idx == (i - 1)) {
        csr_val_cell[csr_row_count_d[i - 1] + j - 1] =
          1.0 - gamma * Jmat_pt[idx * neqs + idx];
        csr_jac_cell[csr_row_count_d[i - 1] + j - 1] =
          Jmat_pt[idx * neqs + idx];
      } else {
        csr_val_cell[csr_row_count_d[i - 1] + j - 1] =
          -gamma * Jmat_pt[idx * neqs + i - 1];
        csr_jac_cell[csr_row_count_d[i - 1] + j - 1] =
          Jmat_pt[idx * neqs + i - 1];
      }
    }
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
fKernelComputeAJsys(
  int ncell, void* user_data, amrex::Real* u_d, amrex::Real* csr_val_arg)
{
  CVODEUserData* udata = static_cast<CVODEUserData*>(user_data);
  auto NNZ = udata->NNZ;
  auto csr_jac_d = udata->csr_jac_d;
  auto csr_row_count_d = udata->csr_row_count_d;
  auto csr_col_index_d = udata->csr_col_index_d;
  auto gamma = udata->gamma;

  int jac_offset = ncell * NNZ;
  amrex::Real* csr_jac_cell = csr_jac_d + jac_offset;
  amrex::Real* csr_val_cell = csr_val_arg + jac_offset;

  for (int i = 1; i < NUM_SPECIES + 2; i++) {
    int nbVals = csr_row_count_d[i] - csr_row_count_d[i - 1];
    for (int j = 0; j < nbVals; j++) {
      int idx = csr_col_index_d[csr_row_count_d[i - 1] + j - 1] - 1;
      if (idx == (i - 1)) {
        csr_val_cell[csr_row_count_d[i - 1] + j - 1] =
          1.0 - gamma * csr_jac_cell[csr_row_count_d[i - 1] + j - 1];
      } else {
        csr_val_cell[csr_row_count_d[i - 1] + j - 1] =
          -gamma * csr_jac_cell[csr_row_count_d[i - 1] + j - 1];
      }
    }
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
fKernelComputeAJchem(
  int ncell, void* user_data, amrex::Real* u_d, amrex::Real* Jdata)
{
  CVODEUserData* udata = static_cast<CVODEUserData*>(user_data);
  const int neqs = NUM_SPECIES + 1;
  auto NNZ = udata->NNZ;
  auto reactor_type = udata->reactor_type;
  auto csr_row_count_d = udata->csr_row_count_d;
  auto csr_col_index_d = udata->csr_col_index_d;

  int u_offset = ncell * neqs;
  amrex::Real* u_curr = u_d + u_offset;

  amrex::Real rho_pt = 0.0;
  for (int n = 0; n < NUM_SPECIES; n++) {
    rho_pt = rho_pt + u_curr[n];
  }

  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {0.0};
  for (int n = 0; n < NUM_SPECIES; n++) {
    massfrac[n] = u_curr[n] / rho_pt;
  }
  amrex::Real temp_pt = u_curr[NUM_SPECIES];

  const int consP = reactor_type == ReactorTypes::h_reactor_type;
  amrex::GpuArray<amrex::Real, neqs* neqs> Jmat_pt = {0.0};
  auto eos = pele::physics::PhysicsType::eos();
  eos.RTY2JAC(rho_pt, temp_pt, massfrac.arr, Jmat_pt.arr, consP);

  // Scale Jacobian and pass into outgoing data ptr.
  amrex::Real mw[NUM_SPECIES] = {0.0};
  get_mw(mw);

  int jac_offset = ncell * NNZ;
  amrex::Real* csr_jac_cell = Jdata + jac_offset;

  for (int i = 0; i < NUM_SPECIES; i++) {
    for (int k = 0; k < NUM_SPECIES; k++) {
      Jmat_pt[k * neqs + i] *= mw[i] / mw[k];
    }
    Jmat_pt[i * neqs + NUM_SPECIES] /= mw[i];
    Jmat_pt[NUM_SPECIES * neqs + i] *= mw[i];
  }
  for (int i = 1; i < NUM_SPECIES + 2; i++) {
    int nbVals = csr_row_count_d[i] - csr_row_count_d[i - 1];
    for (int j = 0; j < nbVals; j++) {
      int idx_cell = csr_col_index_d[csr_row_count_d[i - 1] + j];
      csr_jac_cell[csr_row_count_d[i - 1] + j] =
        Jmat_pt[idx_cell * neqs + i - 1];
    }
  }
}

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
fKernelDenseAJchem(
  int icell, void* user_data, amrex::Real* u_d, amrex::Real* Jdata)
{
  CVODEUserData* udata = static_cast<CVODEUserData*>(user_data);
  const int neqs = NUM_SPECIES + 1;
  auto reactor_type = udata->reactor_type;

  int u_offset = icell * neqs;
  amrex::Real* u_curr = u_d + u_offset;

  amrex::Real rho_pt = 0.0;
  for (int n = 0; n < NUM_SPECIES; n++) {
    rho_pt = rho_pt + u_curr[n];
  }

  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {0.0};
  for (int n = 0; n < NUM_SPECIES; n++) {
    massfrac[n] = u_curr[n] / rho_pt;
  }
  amrex::Real temp_pt = u_curr[NUM_SPECIES];

  const int consP = reactor_type == ReactorTypes::h_reactor_type;
  amrex::GpuArray<amrex::Real, neqs* neqs> Jmat_pt = {0.0};
  auto eos = pele::physics::PhysicsType::eos();
  eos.RTY2JAC(rho_pt, temp_pt, massfrac.arr, Jmat_pt.arr, consP);

  // Scale Jacobian and pass into outgoing data ptr.
  amrex::Real mw[NUM_SPECIES] = {0.0};
  get_mw(mw);

  int jac_offset = icell * (neqs * neqs);
  amrex::Real* Jcurr = Jdata + jac_offset;

  for (int i = 0; i < NUM_SPECIES; i++) {
    for (int k = 0; k < NUM_SPECIES; k++) {
      Jcurr[k * neqs + i] = Jmat_pt[k * neqs + i] * mw[i] / mw[k];
    }
    Jcurr[i * neqs + NUM_SPECIES] = Jmat_pt[i * neqs + NUM_SPECIES] / mw[i];
    Jcurr[NUM_SPECIES * neqs + i] = Jmat_pt[NUM_SPECIES * neqs + i] * mw[i];
  }
}

#ifdef AMREX_USE_CUDA
/*
__global__ void
fKernelDenseSolve( int ncell,
                   amrex::Real* x_d,
                   amrex::Real* b_d,
                   int subsys_size,
                   int subsys_nnz,
                   amrex::Real* csr_val)
{
  int stride = blockDim.x * gridDim.x;

  for (int icell = blockDim.x * blockIdx.x + threadIdx.x; icell < ncell;
       icell += stride) {
    int offset = icell * subsys_size;
    int offset_A = icell * subsys_nnz;

    amrex::Real* csr_val_cell = csr_val + offset_A;
    amrex::Real* x_cell = x_d + offset;
    amrex::Real* b_cell = b_d + offset;

    sgjsolve(csr_val_cell, x_cell, b_cell);
  }
}
*/
#endif

#endif
} // namespace cvode
} // namespace reactions
} // namespace physics
} // namespace pele
#endif
