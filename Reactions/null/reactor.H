#ifndef _REACTOR_H_
#define _REACTOR_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

#include "PelePhysics.H"

#ifdef AMREX_USE_GPU

int reactor_info(int cvode_iE, int Ncells);

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
  const int& reactor_type,
  amrex::gpuStream_t stream);

int react(
  amrex::Real* rY_in,
  amrex::Real* rY_src_in,
  amrex::Real* rX_in,
  amrex::Real* rX_src_in,
  amrex::Real& dt_react,
  amrex::Real& time,
  int cvode_iE,
  int Ncells,
  amrex::gpuStream_t stream);

void reactor_close();

void SetTypValsODE(const std::vector<amrex::Real>& ExtTypVals);

void SetTolFactODE(amrex::Real relative_tol, amrex::Real absolute_tol);

#else

int reactor_init(int cvode_iE, int Ncells);

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
  amrex::Real& time);

int react(
  amrex::Real* rY_in,
  amrex::Real* rY_src_in,
  amrex::Real* rX_in,
  amrex::Real* rX_src_in,
  amrex::Real& dt_react,
  amrex::Real& time);

void reactor_close();

void SetTypValsODE(const std::vector<amrex::Real>& ExtTypVals);

void SetTolFactODE(amrex::Real relative_tol, amrex::Real absolute_tol);

void ReSetTolODE();

#endif
