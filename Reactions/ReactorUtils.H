#ifndef _REACTORUTILS_H_
#define _REACTORUTILS_H_

#include <utility>
#include <AMReX_Gpu.H>
#include "mechanism.H"
#include "PelePhysics.H"
#include "ReactorTypes.H"
#include <cvode/cvode.h>
#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_erkstep.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunmatrix/sunmatrix_sparse.h>
#include <nvector/nvector_serial.h>

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

struct CVODEUserData
{
  amrex::Real dt_save;     // Internal cvode dt holder
  int ncells;              // Number of cells in solve
  int verbose;             // Verbose
  int maxOrder;            // CVODE max order
  int reactor_type;        // Either HP (LM) or UV (C)
  int analytical_jacobian; // Analytical Jacobian 'On' flag
  int solve_type;          // Type of linear solve for Newton direction
  int precond_type;        // Type of preconditioner (if iterative solve)
  int NNZ; // Number of non-zero entry for sparse representations
  amrex::Real* rhoe_init = nullptr; // Initial energy (rhoE for C, rhoH for LM)
  amrex::Real* rhoesrc_ext = nullptr; // External energy forcing
  amrex::Real* rYsrc_ext = nullptr;   // External species forcing
  int* mask =
    nullptr; // Masking tagging cells where integration should not be perfomed
  int* FCunt = nullptr; // Number of RHS evaluations (not used on GPU)
  amrex::Real gamma;    // System Chem. jacobian coefficient
  int nbBlocks;         // GPU kernel launch parameter
  int nbThreads;        // GPU kernel launch parameter

#ifdef AMREX_USE_GPU
  // Matrix data
  int* csr_row_count_h;   // Host CSR sparse mat
  int* csr_col_index_h;   // Host CSR sparse mat
  int* csr_row_count_d;   // Device CSR sparse mat
  int* csr_col_index_d;   // Device CSR sparse mat
  amrex::Real* csr_val_d; // Syst. Jac. CSR sparse mat values
  amrex::Real* csr_jac_d; // Chem. Jac. CSR sparse mat values

#ifdef AMREX_USE_CUDA
  // Batched QR preconditioner data
  void* buffer_qr = NULL;    // cuSolver internal work arrays
  csrqrInfo_t info;          // QR solve info
  cusparseMatDescr_t descrA; // Sparse mat. descriptor
  cusolverSpHandle_t cusolverHandle;
  cusparseHandle_t cuSPHandle;
#endif
  amrex::gpuStream_t stream;

#else
  bool FirstTimePrecond;
  amrex::Real** Jdata = NULL;
  // Sparse direct solve data
  // Ptrs to CSR/CSC matrix for each submatrices (cells)
  SUNMatrix* PS;
  int** rowVals = NULL;
  int** rowPtrs = NULL;
  int** colPtrs = NULL;
  int** colVals = NULL;
  int* indx = NULL;
  amrex::Real** JSPSmat = NULL;
#ifdef PELE_USE_KLU
  klu_common* Common;
  klu_symbolic** Symbolic;
  klu_numeric** Numeric;
#endif
  // Matrix data for denseSimpleAjac preconditioner
  amrex::Real**** Jbd;
  amrex::Real**** P;
  sunindextype*** pivot;
  // Custom direct linear solve: a CSR SUNMatrix and ptrs
  SUNMatrix PSc;
  int* colVals_c;
  int* rowPtrs_c;
#endif
};

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
namespace utils {
int check_flag(void* flagvalue, const char* funcname, int opt);

struct YCOrder
{
  // default: [ [Y_0,...Y_n]_{cell 0}, ..., [Y_0,...Y_n]_{cell m}]
  using ordering_type = YCOrder;
};

struct CYOrder
{
  // [ [Y_{cell 0},...Y_m]_{species 0}, ..., [Y_{cell 0},...Y_m]_{species n}]
  using ordering_type = CYOrder;
};

template <typename T>
struct is_ordering_type
{
  static constexpr bool value =
    std::is_base_of<typename T::ordering_type, T>::value;
};

template <typename OrderType>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE int
vec_index(const int specie_idx, const int icell, const int /*ncells*/)
{
  return icell * (NUM_SPECIES + 1) + specie_idx;
}

template <typename OrderType>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE int
spec_index(const int specie_idx, const int icell, const int /*ncells*/)
{
  return icell * NUM_SPECIES + specie_idx;
}

template <>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE int
vec_index<CYOrder>(const int specie_idx, const int icell, const int ncells)
{
  return specie_idx * ncells + icell;
}

template <>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE int
spec_index<CYOrder>(const int specie_idx, const int icell, const int ncells)
{
  return vec_index<CYOrder>(specie_idx, icell, ncells);
}

template <typename OrderType>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
box_flatten(
  const int icell,
  const int i,
  const int j,
  const int k,
  const int ncells,
  const int reactor_type,
  amrex::Array4<const amrex::Real> const& rhoY,
  amrex::Array4<const amrex::Real> const& frcExt,
  amrex::Array4<const amrex::Real> const& temperature,
  amrex::Array4<const amrex::Real> const& rhoE,
  amrex::Array4<const amrex::Real> const& frcEExt,
  amrex::Real* y_vect,
  amrex::Real* src_vect,
  amrex::Real* vect_energy,
  amrex::Real* src_vect_energy)
{
  amrex::Real mass_frac[NUM_SPECIES];
  amrex::Real rho = 0.0;
  for (int n = 0; n < NUM_SPECIES; n++) {
    const int idx = vec_index<OrderType>(n, icell, ncells);
    const int idx2 = spec_index<OrderType>(n, icell, ncells);
    const amrex::Real rhoYn = rhoY(i, j, k, n);
    y_vect[idx] = rhoYn;
    src_vect[idx2] = frcExt(i, j, k, n);
    rho += rhoYn;
    mass_frac[n] = rhoYn;
  }
  const amrex::Real rho_inv = 1.0 / rho;
  for (int n = 0; n < NUM_SPECIES; n++) {
    mass_frac[n] *= rho_inv;
  }

  amrex::Real temp = temperature(i, j, k, 0);
  vect_energy[icell] = rhoE(i, j, k, 0);
  src_vect_energy[icell] = frcEExt(i, j, k, 0);

  amrex::Real Enrg_loc = vect_energy[icell] * rho_inv;
  auto eos = pele::physics::PhysicsType::eos();
  if (reactor_type == ReactorTypes::e_reactor_type) {
    eos.REY2T(rho, Enrg_loc, mass_frac, temp);
  } else if (reactor_type == ReactorTypes::h_reactor_type) {
    eos.RHY2T(rho, Enrg_loc, mass_frac, temp);
  } else {
    amrex::Abort("Wrong reactor type. Choose between 1 (e) or 2 (h).");
  }
  y_vect[vec_index<OrderType>(NUM_SPECIES, icell, ncells)] = temp;
}

template <typename OrderType>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
box_unflatten(
  const int icell,
  const int i,
  const int j,
  const int k,
  const int ncells,
  const int reactor_type,
  amrex::Array4<amrex::Real> const& rhoY,
  amrex::Array4<amrex::Real> const& temperature,
  amrex::Array4<amrex::Real> const& rhoE,
  amrex::Array4<amrex::Real> const& frcEExt,
  amrex::Array4<amrex::Real> const& FC_in,
  amrex::Real* y_vect,
  amrex::Real* vect_energy,
  long int const& FCunt,
  amrex::Real dt)
{
  amrex::Real mass_frac[NUM_SPECIES];
  amrex::Real rho = 0.0;
  for (int n = 0; n < NUM_SPECIES; n++) {
    const amrex::Real rhoYn = y_vect[vec_index<OrderType>(n, icell, ncells)];
    rhoY(i, j, k, n) = rhoYn;
    rho += rhoYn;
    mass_frac[n] = rhoYn;
  }
  const amrex::Real rho_inv = 1.0 / rho;
  for (int n = 0; n < NUM_SPECIES; n++) {
    mass_frac[n] *= rho_inv;
  }

  amrex::Real temp = y_vect[vec_index<OrderType>(NUM_SPECIES, icell, ncells)];
  rhoE(i, j, k, 0) = vect_energy[icell] + dt * frcEExt(i, j, k, 0);

  amrex::Real Enrg_loc = rhoE(i, j, k, 0) * rho_inv;
  auto eos = pele::physics::PhysicsType::eos();
  if (reactor_type == ReactorTypes::e_reactor_type) {
    eos.REY2T(rho, Enrg_loc, mass_frac, temp);
  } else if (reactor_type == ReactorTypes::h_reactor_type) {
    eos.RHY2T(rho, Enrg_loc, mass_frac, temp);
  } else {
    amrex::Abort("Wrong reactor type. Choose between 1 (e) or 2 (h).");
  }
  temperature(i, j, k, 0) = temp;
  FC_in(i, j, k, 0) = FCunt;
}

template <typename OrderType>
class FlattenOps
{
public:
  static_assert(
    is_ordering_type<OrderType>::value,
    "FlattenOps must have ordering type as its template argument");

  void flatten(
    const amrex::Box box,
    const int ncells,
    const int reactor_type,
    amrex::Array4<const amrex::Real> const& rhoY,
    amrex::Array4<const amrex::Real> const& frcExt,
    amrex::Array4<const amrex::Real> const& temperature,
    amrex::Array4<const amrex::Real> const& rhoE,
    amrex::Array4<const amrex::Real> const& frcEExt,
    amrex::Real* y_vect,
    amrex::Real* src_vect,
    amrex::Real* vect_energy,
    amrex::Real* src_vect_energy)
  {
    BL_PROFILE("Pele::flatten()");
    const auto len = amrex::length(box);
    const auto lo = amrex::lbound(box);
    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const int icell =
        (k - lo.z) * len.x * len.y + (j - lo.y) * len.x + (i - lo.x);
      box_flatten<OrderType>(
        icell, i, j, k, ncells, reactor_type, rhoY, frcExt, temperature, rhoE,
        frcEExt, y_vect, src_vect, vect_energy, src_vect_energy);
    });
  }

  void unflatten(
    const amrex::Box box,
    const int ncells,
    const int reactor_type,
    amrex::Array4<amrex::Real> const& rhoY,
    amrex::Array4<amrex::Real> const& temperature,
    amrex::Array4<amrex::Real> const& rhoE,
    amrex::Array4<amrex::Real> const& frcEExt,
    amrex::Array4<amrex::Real> const& FC_in,
    amrex::Real* y_vect,
    amrex::Real* vect_energy,
    long int* FCunt,
    amrex::Real dt)
  {
    BL_PROFILE("Pele::unflatten()");
    const auto len = amrex::length(box);
    const auto lo = amrex::lbound(box);
    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      const int icell =
        (k - lo.z) * len.x * len.y + (j - lo.y) * len.x + (i - lo.x);
      box_unflatten<OrderType>(
        icell, i, j, k, ncells, reactor_type, rhoY, temperature, rhoE, frcEExt,
        FC_in, y_vect, vect_energy, FCunt[icell], dt);
    });
  }
};

template <typename OrderType>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
fKernelSpec(
  const int icell,
  const int ncells,
  const double dt_save,
  const int reactor_type,
  amrex::Real* yvec_d,
  amrex::Real* ydot_d,
  amrex::Real* rhoe_init,
  amrex::Real* rhoesrc_ext,
  amrex::Real* rYs)
{
  amrex::Real rho_pt = 0.0;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {0.0};
  for (int n = 0; n < NUM_SPECIES; n++) {
    massfrac[n] = yvec_d[vec_index<OrderType>(n, icell, ncells)];
    rho_pt += massfrac[n];
  }
  const amrex::Real rho_pt_inv = 1.0 / rho_pt;

  for (int n = 0; n < NUM_SPECIES; n++) {
    massfrac[n] *= rho_pt_inv;
  }

  const amrex::Real nrg_pt =
    (rhoe_init[icell] + rhoesrc_ext[icell] * dt_save) * rho_pt_inv;

  amrex::Real temp_pt =
    yvec_d[vec_index<OrderType>(NUM_SPECIES, icell, ncells)];

  amrex::Real Cv_pt = 0.0;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> ei_pt = {0.0};
  auto eos = pele::physics::PhysicsType::eos();
  if (reactor_type == ReactorTypes::e_reactor_type) {
    eos.REY2T(rho_pt, nrg_pt, massfrac.arr, temp_pt);
    eos.RTY2Ei(rho_pt, temp_pt, massfrac.arr, ei_pt.arr);
    eos.RTY2Cv(rho_pt, temp_pt, massfrac.arr, Cv_pt);
  } else if (reactor_type == ReactorTypes::h_reactor_type) {
    eos.RHY2T(rho_pt, nrg_pt, massfrac.arr, temp_pt);
    eos.RTY2Hi(rho_pt, temp_pt, massfrac.arr, ei_pt.arr);
    eos.RTY2Cp(rho_pt, temp_pt, massfrac.arr, Cv_pt);
  } else {
    amrex::Abort("Wrong reactor type. Choose between 1 (e) or 2 (h).");
  }

  amrex::GpuArray<amrex::Real, NUM_SPECIES> cdots_pt = {0.0};
  eos.RTY2WDOT(rho_pt, temp_pt, massfrac.arr, cdots_pt.arr);

  amrex::Real rhoesrc = rhoesrc_ext[icell];
  for (int n = 0; n < NUM_SPECIES; n++) {
    const amrex::Real cdot_rYs =
      cdots_pt[n] + rYs[spec_index<OrderType>(n, icell, ncells)];
    ydot_d[vec_index<OrderType>(n, icell, ncells)] = cdot_rYs;
    rhoesrc -= cdot_rYs * ei_pt[n];
  }
  ydot_d[vec_index<OrderType>(NUM_SPECIES, icell, ncells)] =
    rhoesrc * (rho_pt_inv / Cv_pt);
}
} // namespace utils
} // namespace reactions
} // namespace physics
} // namespace pele
#endif
