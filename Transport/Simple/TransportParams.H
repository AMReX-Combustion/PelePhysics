#ifndef TRANSPORT_PARAMS_H
#define TRANSPORT_PARAMS_H

#include <AMReX_REAL.H>
#include "mechanism.h"

struct TransParm {
    amrex::Real* trans_wt;
    amrex::Real* trans_iwt;
    amrex::Real* trans_eps;
    amrex::Real* trans_sig;
    amrex::Real* trans_dip;
    amrex::Real* trans_pol;
    amrex::Real* trans_zrot;
    amrex::Real* trans_fitmu;
    amrex::Real* trans_fitlam;
    amrex::Real* trans_fitdbin;
    int* trans_nlin;
    int array_size = NUM_SPECIES;
    int fit_length = NUM_FIT;
};

extern TransParm *trans_parm_g;
#endif
