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

  using Ordering = utils::YCOrder;

  int init(int reactor_type, int ncells) override;

  void close() override {}

  void print_final_stats(void* /*mem*/) override {}

  int react(
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
    ) override;

  int react(
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
    ) override;

  void flatten(
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
    amrex::Real* src_vect_energy) override
  {
    flatten_ops.flatten(
      box, ncells, m_reactor_type, m_JBBhack, rhoY, frcExt, temperature, rhoE, frcEExt,
      y_vect, src_vect, vect_energy, src_vect_energy);
  }

  void unflatten(
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
    amrex::Real dt) override
  {
    flatten_ops.unflatten(
      box, ncells, m_reactor_type, m_JBBhack, rhoY, temperature, rhoE, frcEExt, FC_in,
      y_vect, vect_energy, FCunt, dt);
  }

private:
  utils::FlattenOps<Ordering> flatten_ops;
  int m_reactor_type{0};
  bool m_JBBhack{false};
};
} // namespace reactions
} // namespace physics
} // namespace pele
#endif
