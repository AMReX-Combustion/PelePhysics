
#ifndef _SOOTMODEL_DERIVE_H_
#define _SOOTMODEL_DERIVE_H_
#include <AMReX_FArrayBox.H>
#include <SootModel.H>

void soot_largeparticledata(
  const amrex::Box& bx,
  amrex::FArrayBox& slfab,
  const int dcomp,
  const int ncomp,
  const amrex::FArrayBox& datafab,
  const amrex::Geometry& geomdata,
  const amrex::Real time,
  const int* bcrec,
  const int level);

void soot_genvars(
  const amrex::Box& bx,
  amrex::FArrayBox& slfab,
  const int dcomp,
  const int ncomp,
  const amrex::FArrayBox& datafab,
  const amrex::Geometry& geomdata,
  const amrex::Real time,
  const int* bcrec,
  const int level);

#endif
