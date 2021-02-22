#ifndef TRANSPORT_PARAMS_H
#define TRANSPORT_PARAMS_H

#include <AMReX_REAL.H>
#include "mechanism.h"

namespace pele {
namespace physics {
namespace transport {
struct TransParm
{
  // FIXME potentially move these into the transport struct constructor
  amrex::Real* trans_wt = nullptr;
  amrex::Real* trans_iwt = nullptr;
  amrex::Real* trans_eps = nullptr;
  amrex::Real* trans_sig = nullptr;
  amrex::Real* trans_dip = nullptr;
  amrex::Real* trans_pol = nullptr;
  amrex::Real* trans_zrot = nullptr;
  amrex::Real* trans_fitmu = nullptr;
  amrex::Real* trans_fitlam = nullptr;
  amrex::Real* trans_fitdbin = nullptr;
  int* trans_nlin = nullptr;
  amrex::Real* Afac = nullptr;
  amrex::Real* Bfac = nullptr;
  amrex::Real* Sigmaij = nullptr;
  amrex::Real* Kappai = nullptr;
  amrex::Real* omega = nullptr;
};

extern TransParm trans_parm;
extern TransParm* trans_parm_g;
} // namespace transport
} // namespace physics
} // namespace pele
#endif
