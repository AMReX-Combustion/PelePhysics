#ifndef _REACTORCVODEJACOBIAN_H_
#define _REACTORCVODEJACOBIAN_H_

#include <AMReX_REAL.H>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_sparse.h>
#include "PelePhysics.H"
#include "mechanism.H"
#include "ReactorTypes.H"
#include "ReactorCvodeUtils.H"

#ifdef AMREX_USE_HIP
#include <nvector/nvector_hip.h>
#endif

#ifdef AMREX_USE_DPCPP
#include <nvector/nvector_sycl.h>
#endif

#ifdef AMREX_USE_CUDA
#include <nvector/nvector_cuda.h>
#include <sunmatrix/sunmatrix_cusparse.h>
#include <sunlinsol/sunlinsol_cusolversp_batchqr.h>
#endif

#ifdef AMREX_USE_GPU
#ifdef PELE_USE_MAGMA
#include <sunmatrix/sunmatrix_magmadense.h>
#include <sunlinsol/sunlinsol_magmadense.h>
#endif
#endif

namespace pele {
namespace physics {
namespace reactions {
namespace cvode {
int cJac(
  amrex::Real tn,
  N_Vector u,
  N_Vector fy,
  SUNMatrix J,
  void* user_data,
  N_Vector tmp1,
  N_Vector tmp2,
  N_Vector tmp3);

#ifndef AMREX_USE_GPU
int cJac_sps(
  amrex::Real tn,
  N_Vector u,
  N_Vector fy,
  SUNMatrix J,
  void* user_data,
  N_Vector tmp1,
  N_Vector tmp2,
  N_Vector tmp3);

#ifdef PELE_USE_KLU
int cJac_KLU(
  amrex::Real tn,
  N_Vector y,
  N_Vector fy,
  SUNMatrix J,
  void* user_data,
  N_Vector tmp1,
  N_Vector tmp2,
  N_Vector tmp3);
#endif
#endif
} // namespace cvode
} // namespace reactions
} // namespace physics
} // namespace pele
#endif
