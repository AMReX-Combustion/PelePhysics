#ifndef _REACTORUTILS_H_
#define _REACTORUTILS_H_

#include <AMReX_Gpu.H>
#include "mechanism.H"
#include "PelePhysics.H"

namespace pele {
namespace physics {
namespace reactions {

int check_flag(void* flagvalue, const char* funcname, int opt);

template <typename ReactorType>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE int
sunvec_index(const int specie_idx, const int icell, const int /*ncells*/)
{
  return icell * (NUM_SPECIES + 1) + specie_idx;
}

template <typename ReactorType>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE int
sunspec_index(const int specie_idx, const int icell, const int /*ncells*/)
{
  return icell * NUM_SPECIES + specie_idx;
}

class ReactorArkode;

template <>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE int
sunvec_index<ReactorArkode>(
  const int specie_idx, const int icell, const int ncells)
{
  return specie_idx * ncells + icell;
}

template <>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE int
sunspec_index<ReactorArkode>(
  const int specie_idx, const int icell, const int ncells)
{
  return sunvec_index<ReactorArkode>(specie_idx, icell, ncells);
}

template <typename ReactorType>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
box_flatten(
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
    const int idx = sunvec_index<ReactorType>(n, icell, ncells);
    const int idx2 = sunspec_index<ReactorType>(n, icell, ncells);
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
    eos.REY2T(rho, Enrg_loc, mass_frac, temp);
  } else {
    eos.RHY2T(rho, Enrg_loc, mass_frac, temp);
  }
  tmp_vect[sunvec_index<ReactorType>(NUM_SPECIES, icell, ncells)] = temp;
}

template <typename ReactorType>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
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
    const amrex::Real rhoYn =
      tmp_vect[sunvec_index<ReactorType>(n, icell, ncells)];
    rhoY(i, j, k, n) = rhoYn;
    rho += rhoYn;
    mass_frac[n] = rhoYn;
  }
  const amrex::Real rho_inv = 1.0 / rho;
  for (int n = 0; n < NUM_SPECIES; n++) {
    mass_frac[n] *= rho_inv;
  }

  amrex::Real temp =
    tmp_vect[sunvec_index<ReactorType>(NUM_SPECIES, icell, ncells)];
  rhoE(i, j, k, 0) = tmp_vect_energy[icell] + dt * frcEExt(i, j, k, 0);

  amrex::Real Enrg_loc = rhoE(i, j, k, 0) * rho_inv;
  auto eos = pele::physics::PhysicsType::eos();
  if (ireactor_type == 1) {
    eos.REY2T(rho, Enrg_loc, mass_frac, temp);
  } else {
    eos.RHY2T(rho, Enrg_loc, mass_frac, temp);
  }
  T_in(i, j, k, 0) = temp;
  FC_in(i, j, k, 0) = tmp_FCunt;
}

template <typename ReactorType>
AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
fKernelSpec(
  int icell,
  int ncells,
  double dt_save,
  int reactor_type,
  amrex::Real* yvec_d,
  amrex::Real* ydot_d,
  amrex::Real* rhoe_init,
  amrex::Real* rhoesrc_ext,
  amrex::Real* rYs)
{
  amrex::Real rho_pt = 0.0;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> massfrac = {0.0};
  for (int n = 0; n < NUM_SPECIES; n++) {
    massfrac[n] = yvec_d[sunvec_index<ReactorType>(n, icell, ncells)];
    rho_pt += massfrac[n];
  }
  const amrex::Real rho_pt_inv = 1.0 / rho_pt;

  for (int n = 0; n < NUM_SPECIES; n++) {
    massfrac[n] *= rho_pt_inv;
  }

  const amrex::Real nrg_pt =
    (rhoe_init[icell] + rhoesrc_ext[icell] * dt_save) * rho_pt_inv;

  amrex::Real temp_pt =
    yvec_d[sunvec_index<ReactorType>(NUM_SPECIES, icell, ncells)];

  amrex::Real Cv_pt = 0.0;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> ei_pt = {0.0};
  auto eos = pele::physics::PhysicsType::eos();
  if (reactor_type == 1) {
    eos.REY2T(rho_pt, nrg_pt, massfrac.arr, temp_pt);
    eos.RTY2Ei(rho_pt, temp_pt, massfrac.arr, ei_pt.arr);
    eos.RTY2Cv(rho_pt, temp_pt, massfrac.arr, Cv_pt);
  } else {
    eos.RHY2T(rho_pt, nrg_pt, massfrac.arr, temp_pt);
    eos.RTY2Hi(rho_pt, temp_pt, massfrac.arr, ei_pt.arr);
    eos.RTY2Cp(rho_pt, temp_pt, massfrac.arr, Cv_pt);
  }

  amrex::GpuArray<amrex::Real, NUM_SPECIES> cdots_pt = {0.0};
  eos.RTY2WDOT(rho_pt, temp_pt, massfrac.arr, cdots_pt.arr);

  amrex::Real rhoesrc = rhoesrc_ext[icell];
  for (int n = 0; n < NUM_SPECIES; n++) {
    const amrex::Real cdot_rYs =
      cdots_pt[n] + rYs[sunspec_index<ReactorType>(n, icell, ncells)];
    ydot_d[sunvec_index<ReactorType>(n, icell, ncells)] = cdot_rYs;
    rhoesrc -= cdot_rYs * ei_pt[n];
  }
  ydot_d[sunvec_index<ReactorType>(NUM_SPECIES, icell, ncells)] =
    rhoesrc * (rho_pt_inv / Cv_pt);
}
} // namespace reactions
} // namespace physics
} // namespace pele
#endif
