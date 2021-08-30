#ifndef _REACTORARKODE_H_
#define _REACTORARKODE_H_

#include "ReactorBase.H"
#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_erkstep.h>
#include <nvector/nvector_serial.h>

namespace pele {
namespace physics {
namespace reactions {

struct ARKODEUserData
{
  bool reactor_arkode_initialized;
  int ncells_d;
  int neqs_per_cell;
  int iverbose;
  int ireactor_type;
  amrex::Real dt_save;
  amrex::Real* rhoe_init_d = NULL;
  amrex::Real* rhoesrc_ext_d = NULL;
  amrex::Real* rYsrc_d = NULL;
#ifdef AMREX_USE_GPU
  amrex::gpuStream_t stream;
#endif
  int nbBlocks;
  int nbThreads;
};

class ReactorArkode : public ReactorBase::Register<ReactorArkode>
{
public:
  static std::string identifier() { return "ReactorArkode"; }

  int init(int reactor_type, int Ncells) override;

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

  static int
  cF_RHS(amrex::Real t, N_Vector y_in, N_Vector ydot, void* user_data);

  void SetTypValsODE(const std::vector<amrex::Real>& ExtTypVals) override;

private:
  amrex::Real relTol = 1e-6;
  amrex::Real absTol = 1e-10;
  int use_erkstep = 0;
  int rk_method = 40;
  int rk_controller = 0;
};

} // namespace reactions
} // namespace physics
} // namespace pele

#endif
