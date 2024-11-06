#ifndef REACTORBASE_H
#define REACTORBASE_H

#include <AMReX_ParmParse.H>
#include <AMReX_Gpu.H>
#include <sundials/sundials_context.h>
#include "Factory.H"
#include "PelePhysics.H"
#include "mechanism.H"
#include "ReactorUtils.H"

namespace pele::physics::reactions {

class ReactorBase : public Factory<ReactorBase>
{
public:
  static std::string base_identifier() { return "ReactorBase"; }

  virtual int init(int reactor_type, int ncells) = 0;

  virtual void close() = 0;

  virtual void print_final_stats(void* sundials_mem) = 0;

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
    const amrex::Box& box,
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
    const amrex::Box& box,
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

  void set_typ_vals_ode(const std::vector<amrex::Real>& ExtTypVals);

  // Manifold EOS needs an eosparm - right now that is only propagated through
  // the RK64 reactor
  virtual void set_eos_parm(
    const pele::physics::eos::EosParm<pele::physics::PhysicsType::eos_type>*
      h_eosparm,
    const pele::physics::eos::EosParm<pele::physics::PhysicsType::eos_type>*
      d_eosparm)
  {
    auto eos = pele::physics::PhysicsType::eos(h_eosparm);
    if (eos.identifier() == "Manifold") {
      amrex::Abort(
        "Manifold EOS only supported with ReactorRK64 and ReactorNull for now");
    } else {
      m_h_eosparm = h_eosparm;
      m_d_eosparm = d_eosparm;
    }
  }

  ~ReactorBase() override = default;

protected:
  int verbose{0};
  amrex::GpuArray<amrex::Real, NUM_SPECIES + 1> m_typ_vals = {0.0};
  const pele::physics::eos::EosParm<pele::physics::PhysicsType::eos_type>*
    m_h_eosparm;
  const pele::physics::eos::EosParm<pele::physics::PhysicsType::eos_type>*
    m_d_eosparm;
};
} // namespace pele::physics::reactions
#endif
