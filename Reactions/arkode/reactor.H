#ifndef _REACTOR_H_
#define _REACTOR_H_

#include <AMReX_Gpu.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_erkstep.h>
#include <nvector/nvector_serial.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>
#include <sunnonlinsol/sunnonlinsol_newton.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

#include "mechanism.H"
#include "PelePhysics.H"

#ifdef AMREX_USE_GPU
#include "AMReX_SUNMemory.H"
#endif
#if defined(AMREX_USE_CUDA)
#include <nvector/nvector_cuda.h>
#elif defined(AMREX_USE_HIP)
#include <nvector/nvector_hip.h>
#endif

typedef struct ARKODEUserData
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
} * UserData;

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

int reactor_init(int cvode_iE, int Ncells);

int cF_RHS(realtype t, N_Vector y_in, N_Vector ydot, void* user_data);

int check_flag(void* flagvalue, const char* funcname, int opt);

// React for array4
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

// React for 1d array
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

#endif
