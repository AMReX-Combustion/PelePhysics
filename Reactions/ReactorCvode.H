#ifndef _REACTORCVODE_H_
#define _REACTORCVODE_H_

#include "AMReX_Sundials.H"
#include "ReactorBase.H"
#include "ReactorCvodeUtils.H"
#include "ReactorCvodePreconditioner.H"
#include "ReactorCvodeJacobian.H"
#include "ReactorCvodeCustomLinSolver.H"
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

#ifdef AMREX_USE_HIP
#include <nvector/nvector_hip.h>
#endif

#ifdef AMREX_USE_DPCPP
#include <nvector/nvector_sycl.h>
#endif

#ifdef PELE_USE_KLU
#include "klu.h"
#include <sunlinsol/sunlinsol_klu.h>
#endif

#ifdef PELE_USE_MAGMA
#include <sunmatrix/sunmatrix_magmadense.h>
#include <sunlinsol/sunlinsol_magmadense.h>
#endif

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

namespace pele {
namespace physics {
namespace reactions {

class ReactorCvode : public ReactorBase::Register<ReactorCvode>
{
public:
  static std::string identifier() { return "ReactorCvode"; }

  using Ordering = utils::YCOrder;

  int init(int reactor_type, int ncells) override;

  void close() override;

  void print_final_stats(void* cvodemem, bool print_ls_stats);
  void print_final_stats(void* cvodemem) override
  {
    print_final_stats(cvodemem, true);
  }

  int react(
    const amrex::Box& box,
    amrex::Array4<amrex::Real> const& rY_in,
    amrex::Array4<amrex::Real> const& rYsrc_in,
    amrex::Array4<amrex::Real> const& T_in,
    amrex::Array4<amrex::Real> const& rEner_in,
    amrex::Array4<amrex::Real> const& rEner_src_in,
    amrex::Array4<amrex::Real> const& FC_in,
    amrex::Array4<int> const& mask,
    amrex::Real& dt_react,
    amrex::Real& time
#ifdef AMREX_USE_GPU
    ,
    amrex::gpuStream_t stream
#endif
    ) override;

  int react(
    amrex::Real* rY_in,
    amrex::Real* rYsrc_in,
    amrex::Real* rX_in,
    amrex::Real* rX_src_in,
    amrex::Real& dt_react,
    amrex::Real& time,
    int ncells
#ifdef AMREX_USE_GPU
    ,
    amrex::gpuStream_t stream
#endif
    ) override;

  static int
  cF_RHS(amrex::Real t, N_Vector y_in, N_Vector ydot, void* user_data);

  void flatten(
    const amrex::Box box,
    const int ncells,
    amrex::Array4<const amrex::Real> const& rhoY,
    amrex::Array4<const amrex::Real> const& frcExt,
    amrex::Array4<const amrex::Real> const& temperature,
    amrex::Array4<const amrex::Real> const& rhoE,
    amrex::Array4<const amrex::Real> const& frcEExt,
    amrex::Real* y_vect,
    amrex::Real* src_vect,
    amrex::Real* vect_energy,
    amrex::Real* src_vect_energy) override
  {
    flatten_ops.flatten(
      box, ncells, m_reactor_type, rhoY, frcExt, temperature, rhoE, frcEExt,
      y_vect, src_vect, vect_energy, src_vect_energy);
  }

  void unflatten(
    const amrex::Box box,
    const int ncells,
    amrex::Array4<amrex::Real> const& rhoY,
    amrex::Array4<amrex::Real> const& temperature,
    amrex::Array4<amrex::Real> const& rhoE,
    amrex::Array4<amrex::Real> const& frcEExt,
    amrex::Array4<amrex::Real> const& FC_in,
    amrex::Real* y_vect,
    amrex::Real* vect_energy,
    long int* FCunt,
    amrex::Real dt) override
  {
    flatten_ops.unflatten(
      box, ncells, m_reactor_type, rhoY, temperature, rhoE, frcEExt, FC_in,
      y_vect, vect_energy, FCunt, dt);
  }

private:
  void checkCvodeOptions() const;
  void allocUserData(
    CVODEUserData* udata,
    int a_ncells
#ifdef AMREX_USE_GPU
    ,
    SUNMatrix& a_A,
    amrex::gpuStream_t stream
#endif
  ) const;
  void freeUserData(CVODEUserData* data_wk);

  amrex::Real relTol{1e-6};
  amrex::Real absTol{1e-10};
  utils::FlattenOps<Ordering> flatten_ops;
  int m_reactor_type{0};
#if defined(AMREX_USE_HIP)
  int atomic_reductions{0};
#else
  int atomic_reductions{1};
#endif

#ifndef AMREX_USE_GPU
  N_Vector y = nullptr;
  SUNLinearSolver LS = nullptr;
  SUNLinearSolver NLS = nullptr;
  SUNMatrix A = nullptr;
  void* cvode_mem = nullptr;
  CVODEUserData* udata_g = nullptr;
#endif
};
} // namespace reactions
} // namespace physics
} // namespace pele
#endif
