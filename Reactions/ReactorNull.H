#ifndef _REACTORNULL_H_
#define _REACTORNULL_H_

#include "ReactorBase.H"

namespace pele {
namespace physics {
namespace reactions {

class ReactorNull : public ReactorBase::Register<ReactorNull>
{
public:
  static std::string identifier() { return "ReactorNull"; }

  int init(int reactor_type, int Ncells) override { return 0; }

  void close() override {}

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
    ) override;

  int react(
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
    ) override;

  void SetTypValsODE(const std::vector<amrex::Real>& ExtTypVals) override;
};
} // namespace reactions
} // namespace physics
} // namespace pele
#endif
