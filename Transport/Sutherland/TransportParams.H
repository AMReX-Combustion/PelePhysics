#ifndef TRANSPORT_PARAMS_H
#define TRANSPORT_PARAMS_H

#include <AMReX_REAL.H>

struct TransParm {
    amrex::Real Prandtl_number;
    amrex::Real viscosity_mu_ref;
    amrex::Real viscosity_T_ref;
    amrex::Real viscosity_S;
    amrex::Real const_bulk_viscosity;
    amrex::Real const_diffusivity;
};

extern TransParm *trans_parm_g;
#endif
