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

  using Ordering = YCOrder;

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
    amrex::Real* src_vect_energy) override
  {
    flatten_ops.flatten(
      box, ncells, reactor_type, rhoY, frcExt, temperature, rhoE, frcEExt,
      y_vect, src_vect, vect_energy, src_vect_energy);
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
    amrex::Real dt) override
  {
    flatten_ops.unflatten(
      box, ncells, reactor_type, rhoY, temperature, rhoE, frcEExt, FC_in,
      y_vect, vect_energy, FCunt, dt);
  }

  void SetTypValsODE(const std::vector<amrex::Real>& ExtTypVals) override;

private:
  amrex::Real relTol = 1e-6;
  amrex::Real absTol = 1e-10;
  amrex::Array<double, NUM_SPECIES + 1> typVals;
  FlattenOps<Ordering> flatten_ops;
};
} // namespace reactions
} // namespace physics
} // namespace pele
#endif
