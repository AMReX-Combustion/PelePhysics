#ifndef REACTORRK64_H
#define REACTORRK64_H

#include "ReactorBase.H"

namespace pele::physics::reactions {

// stages and coefficients for rk64
struct RK64Params
{
  const int nstages_rk64 = 6;
  const amrex::Real alpha_rk64[6] = {
    0.218150805229859,  // 3296351145737.0/15110423921029.0,
    0.256702469801519,  // 1879360555526.0/ 7321162733569.0,
    0.527402592007520,  // 10797097731880.0/20472212111779.0,
    0.0484864267224467, // 754636544611.0/15563872110659.0,
    1.24517071533530,   // 3260218886217.0/ 2618290685819.0,
    0.412366034843237,  // 5069185909380.0/12292927838509.0
  };

  const amrex::Real beta_rk64[6] = {
    -0.113554138044166,  //-1204558336989.0/10607789004752.0,
    -0.215118587818400,  //-3028468927040.0/14078136890693.0,
    -0.0510152146250577, //-455570672869.0/ 8930094212428.0,
    -1.07992686223881,   //-17275898420483.0/15997285579755.0,
    -0.248664241213447,  //-2453906524165.0/ 9868353053862.0,
    0.0};

  const amrex::Real err_rk64[6] = {
    -0.0554699315064507, //-530312978447.0/ 9560368366154.0,
    0.158481845574980,   // 473021958881.0/ 2984707536468.0,
    -0.0905918835751907, //-947229622805.0/10456009803779.0,
    -0.219084567203338,  //-2921473878215.0/13334914072261.0,
    0.164022338959433,   // 1519535112975.0/ 9264196100452.0,
    0.0426421977505659   // 167623581683.0/ 3930932046784.0
  };
  const amrex::Real exp1_rk64 = 0.25;
  const amrex::Real exp2_rk64 = 0.2;
  const amrex::Real betaerr_rk64 = 1.0;
};

class ReactorRK64 : public ReactorBase::Register<ReactorRK64>
{
public:
  static std::string identifier() { return "ReactorRK64"; }

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
    amrex::Real* src_vect_energy) override
  {
    flatten_ops.flatten(
      box, ncells, m_reactor_type, m_clean_init_massfrac, rhoY, frcExt,
      temperature, rhoE, frcEExt, y_vect, src_vect, vect_energy,
      src_vect_energy);
  }

  void unflatten(
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
    amrex::Real dt) override
  {
    flatten_ops.unflatten(
      box, ncells, m_reactor_type, m_clean_init_massfrac, rhoY, temperature,
      rhoE, frcEExt, FC_in, y_vect, vect_energy, FCunt, dt);
  }

  void set_eos_parm(
    const pele::physics::eos::EosParm<pele::physics::PhysicsType::eos_type>*
      h_eosparm,
    const pele::physics::eos::EosParm<pele::physics::PhysicsType::eos_type>*
      d_eosparm) override
  {
    m_h_eosparm = h_eosparm;
    m_d_eosparm = d_eosparm;
  }

private:
  amrex::Real absTol{1e-10};
  int rk64_nsubsteps_guess{10};
  int rk64_nsubsteps_min{5};
  int rk64_nsubsteps_max{500};
  utils::FlattenOps<Ordering> flatten_ops;
  int m_reactor_type{0};
  bool m_clean_init_massfrac{false};
};
} // namespace pele::physics::reactions
#endif
