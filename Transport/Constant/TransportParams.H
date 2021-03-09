#ifndef TRANSPORT_PARAMS_H
#define TRANSPORT_PARAMS_H

#include <AMReX_REAL.H>
namespace pele {
namespace physics {
namespace transport {

struct TransParm
{
  amrex::Real const_viscosity;
  amrex::Real const_bulk_viscosity;
  amrex::Real const_diffusivity;
  amrex::Real const_conductivity;
};

extern TransParm* trans_parm_g;
} // namespace transport
} // namespace physics
} // namespace pele
#endif
