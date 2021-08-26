#ifndef _REACTORCVODE_H_
#define _REACTORCVODE_H_

#include "ReactorBase.H"
#include <nvector/nvector_serial.h>

namespace pele {
namespace physics {
namespace reactions {

class ReactorCvode : public ReactorBase::Register<ReactorCvode>
{
public:
  static std::string identifier() { return "ReactorCvode"; }

  int init(int reactor_type, int Ncells) override;

  void close() override;

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

  static int
  cF_RHS(amrex::Real t, N_Vector y_in, N_Vector ydot, void* user_data);

  void SetTypValsODE(const std::vector<amrex::Real>& ExtTypVals);
  void SetTolFactODE(amrex::Real relative_tol, amrex::Real absolute_tol);

private:
  amrex::Real relTol = 1e-6;
  amrex::Real absTol = 1e-10;
  amrex::Array<double, NUM_SPECIES + 1> typVals;
};
} // namespace reactions
} // namespace physics
} // namespace pele
#endif
