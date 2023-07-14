#ifndef REACTORCVODE_H
#define REACTORCVODE_H

#include "AMReX_Sundials.H"
#include "ReactorBase.H"
#include "ReactorCvodeUtils.H"
#include "ReactorCvodePreconditioner.H"
#include "ReactorCvodeJacobian.H"
#include "ReactorCvodeCustomLinSolver.H"
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_sparse.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunlinsol/sunlinsol_spgmr.h>
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>

#ifdef AMREX_USE_HIP
#include <nvector/nvector_hip.h>
#endif

#ifdef AMREX_USE_SYCL
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

namespace pele::physics::reactions {

class ReactorCvode : public ReactorBase::Register<ReactorCvode>
{
public:
  static std::string identifier() { return "ReactorCvode"; }

#if defined(AMREX_USE_GPU) && !defined(PELE_USE_MAGMA)
  using Ordering = utils::CYOrder;
#else
  using Ordering = utils::YCOrder;
#endif

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
    const amrex::Box& box,
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
      box, ncells, m_reactor_type, m_clean_init_massfrac, rhoY, frcExt,
      temperature, rhoE, frcEExt, y_vect, src_vect, vect_energy,
      src_vect_energy);
  }

  void unflatten(
    const amrex::Box& box,
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
      box, ncells, m_reactor_type, m_clean_init_massfrac, rhoY, temperature,
      rhoE, frcEExt, FC_in, y_vect, vect_energy, FCunt, dt);
  }

private:
  void checkCvodeOptions(
    const std::string& a_solve_type_str,
    const std::string& a_precond_type_str,
    int& a_solve_type,
    int& a_ajac,
    int& a_precond_type) const;
  int initCvode(
    N_Vector& a_y,
    SUNMatrix& a_A,
    CVODEUserData* a_udata,
    SUNNonlinearSolver& a_NLS,
    SUNLinearSolver& a_LS,
    void* a_cvode_mem,
#ifdef AMREX_USE_GPU
    amrex::gpuStream_t stream,
#endif
    const amrex::Real& a_time,
    int ncells);
  void allocUserData(
    CVODEUserData* udata,
    int a_ncells
#ifdef AMREX_USE_GPU
    ,
    SUNMatrix& a_A,
    amrex::gpuStream_t stream
#endif
  ) const;
  static void freeUserData(CVODEUserData* data_wk);

  amrex::Real relTol{1e-6};
  amrex::Real absTol{1e-10};
  utils::FlattenOps<Ordering> flatten_ops;
  int m_reactor_type{0};
#if defined(AMREX_USE_HIP)
  int atomic_reductions{0};
#else
  int atomic_reductions{1};
#endif
  int max_nls_iters{3};
  int max_fp_accel{2};
  int m_print_profiling{0};
  int m_cvode_maxorder{2};
  int m_solve_type{-1};
  int m_analytical_jacobian{-1};
  int m_precond_type{-1};
  bool m_clean_init_massfrac{false};
};
} // namespace pele::physics::reactions
#endif
