#ifndef _REACTORBASE_H_
#define _REACTORBASE_H_

#include <AMReX_ParmParse.H>

#include "Factory.H"
#include "PelePhysics.H"
#include "mechanism.H"
#include "ReactorUtils.H"

namespace pele {
namespace physics {
namespace reactions {

class ReactorBase : public Factory<ReactorBase>
{
public:
  amrex::Array<double, NUM_SPECIES + 1> typVals;

  static const std::string base_identifier() { return "ReactorBase"; }

  virtual int init(int reactor_type, int ncells) = 0;

  virtual void close() = 0;

  virtual int react(
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
    ) = 0;

  virtual int react(
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
    ) = 0;

  virtual void flatten(
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
    amrex::Real* src_vect_energy) = 0;

  virtual void unflatten(
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
    amrex::Real dt) = 0;

  void SetTypValsODE(const std::vector<amrex::Real>& ExtTypVals);

  void setSundialsSolverTols(void* sundials_mem, int ncells, int verbose, amrex::Real relTol,
          amrex::Real absTol, std::string solvername);

  virtual ~ReactorBase() = default;
};
} // namespace reactions
} // namespace physics
} // namespace pele
#endif
