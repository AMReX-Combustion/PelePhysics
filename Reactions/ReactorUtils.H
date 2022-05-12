#ifndef _REACTORUTILS_H_
#define _REACTORUTILS_H_

#include <utility>
#include <AMReX_Gpu.H>
#include "mechanism.H"
#include "PelePhysics.H"
#include "ReactorTypes.H"
#include "AMReX_Sundials.H"

#ifdef AMREX_USE_HIP
#include <nvector/nvector_hip.h>
#endif

#ifdef AMREX_USE_DPCPP
#include <nvector/nvector_sycl.h>
#endif

#ifdef AMREX_USE_CUDA
#include <nvector/nvector_cuda.h>
#endif

namespace pele {
namespace physics {
namespace reactions {
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
  const bool JBBhack,
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
    const amrex::Real rhoYn = (JBBhack) ? amrex::max(0.0,rhoY(i, j, k, n)) : rhoY(i, j, k, n);
    y_vect[idx] = rhoYn;
    src_vect[idx2] = frcExt(i, j, k, n);
    rho += rhoYn;
    mass_frac[n] = rhoYn;
  }
  const amrex::Real rho_inv = 1.0 / rho;
  for (amrex::Real& Yn : mass_frac) {
    Yn *= rho_inv;
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
  const bool JBBhack,
  amrex::Array4<amrex::Real> const& rhoY,
  amrex::Array4<amrex::Real> const& temperature,
  amrex::Array4<amrex::Real> const& rhoE,
  amrex::Array4<amrex::Real> const& frcEExt,
  amrex::Array4<amrex::Real> const& FC_in,
  const amrex::Real* y_vect,
  const amrex::Real* vect_energy,
  long int const& FCunt,
  amrex::Real dt)
{
  // Get the old Ys & rho as well as the clipped version
  // that went into CVODE
  amrex::Real rhoY_oldclipped[NUM_SPECIES];
  amrex::Real rho_old = 0.0;
  amrex::Real rho_oldclipped = 0.0;
  if (JBBhack) {
    for (int n = 0; n < NUM_SPECIES; n++) {
       rhoY_oldclipped[n] = amrex::max(0.0,rhoY(i, j, k, n));
       rho_old += rhoY(i, j, k, n);
       rho_oldclipped += rhoY_oldclipped[n];
    }
  }

  amrex::Real mass_frac[NUM_SPECIES];
  amrex::Real rho = 0.0;
  for (int n = 0; n < NUM_SPECIES; n++) {
    // With JBBhack, new rhoY equal old one + CVODE increment
    const amrex::Real rhoYn = (JBBhack)? rhoY(i, j, k, n)
                                + (  y_vect[vec_index<OrderType>(n, icell, ncells)]
                                     - rhoY_oldclipped[n]) * rho_old / rho_oldclipped : y_vect[vec_index<OrderType>(n, icell, ncells)];
    rhoY(i, j, k, n) = rhoYn;
    rho += rhoYn;
    mass_frac[n] = rhoYn;
  }
  const amrex::Real rho_inv = 1.0 / rho;
  for (amrex::Real& Yn : mass_frac) {
    Yn *= rho_inv;
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
    const bool JBBhack,
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
        icell, i, j, k, ncells, reactor_type, JBBhack, rhoY, frcExt, temperature, rhoE,
        frcEExt, y_vect, src_vect, vect_energy, src_vect_energy);
    });
  }

  void unflatten(
    const amrex::Box box,
    const int ncells,
    const int reactor_type,
    const bool JBBhack,
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
        icell, i, j, k, ncells, reactor_type, JBBhack, rhoY, temperature, rhoE, frcEExt,
        FC_in, y_vect, vect_energy, FCunt[icell], dt);
    });
  }
};

template <typename OrderType>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
fKernelSpec(
  const int icell,
  const int ncells,
  const amrex::Real dt_save,
  const int reactor_type,
  const amrex::Real* yvec_d,
  amrex::Real* ydot_d, // NOLINT(readability-non-const-parameter)
  const amrex::Real* rhoe_init,
  const amrex::Real* rhoesrc_ext,
  const amrex::Real* rYs)
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

#ifdef AMREX_USE_GPU
N_Vector
setNVectorGPU(int nvsize, int atomic_reduction, amrex::gpuStream_t stream);
#endif
} // namespace utils
} // namespace reactions
} // namespace physics
} // namespace pele
#endif
