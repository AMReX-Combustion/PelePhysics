#ifndef REACTORBDF_H
#define REACTORBDF_H

#define BDF1SCHEME 0
#define TRPZSCHEME 1
#define BDF2SCHEME 2
#define BDF3SCHEME 3

#define NBDFSCHEMES 4
#define NBDFSTENCIL 4

#include "ReactorBase.H"
// stages and coefficients for rk64
struct BDFParams
{
   //3 schemes are implemented
   //BDF1, TRPZ, BDF2
   
   //mass matrix multipliers
   const amrex::Real TCOEFFMAT[NBDFSCHEMES][NBDFSTENCIL]
    = { {1.0000,-1.0000, 0.0000,  0.0000},
        {1.0000,-1.0000, 0.0000,  0.0000},
        {1.0000,-1.3333, 0.3333,  0.0000},
        {1.0000,-1.6363, 0.8181, -0.1818}
      };

    //forcing multipliers
    const amrex::Real FCOEFFMAT[NBDFSCHEMES][NBDFSTENCIL]
     = { {1.0000, 0.0000, 0.0000, 0.0000},
         {0.5000, 0.5000, 0.0000, 0.0000},
         {0.6666, 0.0000, 0.0000, 0.0000},
         {0.5454, 0.0000, 0.0000, 0.0000}
       };
};

namespace pele::physics::reactions {

class ReactorBDF : public ReactorBase::Register<ReactorBDF>
{
public:
  static std::string identifier() { return "ReactorBDF"; }

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

private:
  amrex::Real m_nonlin_tol{1e-6};
  int m_nsubsteps{1};
  int m_nonlinear_iters{1};
  int m_gmres_restarts{10};
  amrex::Real m_gmres_tol{1e-10};
  int m_reactor_type{0};
  bool m_clean_init_massfrac{false};
  utils::FlattenOps<Ordering> flatten_ops;
  int m_gmres_precond{2};
  int m_tstepscheme{BDF3SCHEME};
};
} // namespace pele::physics::reactions
#endif
