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
  static const std::string base_identifier() { return "Reactor"; }

  virtual int init(int reactor_type, int Ncells) = 0;

  virtual void close() = 0;

  virtual int react(
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
    ) = 0;

  virtual int react(
    amrex::Real* rY_in,
    amrex::Real* rY_src_in,
    amrex::Real* rX_in,
    amrex::Real* rX_src_in,
    amrex::Real& dt_react,
    amrex::Real& time,
    int reactor_type,
    int Ncells
#ifdef AMREX_USE_GPU
    ,
    amrex::gpuStream_t stream
#endif
    ) = 0;

  virtual ~ReactorBase() = default;
};

} // namespace reactions
} // namespace physics
} // namespace pele
#endif
