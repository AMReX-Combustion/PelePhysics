#ifndef REACTORARKODE_H
#define REACTORARKODE_H

#include "ReactorBase.H"
#include "AMReX_Sundials.H"
#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_erkstep.h>
#include <nvector/nvector_serial.h>

#ifdef AMREX_USE_DPCPP
#include <nvector/nvector_sycl.h>
#endif

#ifdef AMREX_USE_HIP
#include <nvector/nvector_hip.h>
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

struct ARKODEUserData
{
  int ncells;
  int neq;
  int verbose;
  int reactor_type;
  amrex::Real dt_save;
  amrex::Real* rhoe_init = nullptr;
  amrex::Real* rhoesrc_ext = nullptr;
  amrex::Real* rYsrc_ext = nullptr;
};

class ReactorArkode : public ReactorBase::Register<ReactorArkode>
{
public:
  static std::string identifier() { return "ReactorArkode"; }

#ifdef AMREX_USE_GPU
  using Ordering = utils::CYOrder;
#else
  using Ordering = utils::YCOrder;
#endif

  int init(int reactor_type, int ncells) override;

  void close() override {}

  void print_final_stats(void* arkode_mem) override;

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
    const& amrex::Box box,
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
    const& amrex::Box box,
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
  amrex::Real relTol{1e-6};
  amrex::Real absTol{1e-10};
  int use_erkstep{1};
  int rk_method{40};
  int rk_controller{0};
  utils::FlattenOps<Ordering> flatten_ops;
  int m_reactor_type{0};
#if defined(AMREX_USE_HIP)
  int atomic_reductions{0};
#else
  int atomic_reductions{1};
#endif
  bool m_clean_init_massfrac{false};
};

} // namespace reactions
} // namespace physics
} // namespace pele

#endif
