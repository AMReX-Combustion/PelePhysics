#include <math.h>
#include <iostream>
#include <cstring>
#include <chrono>
#include <cassert>
#include <assert.h>

#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_erkstep.h>

#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */

#ifdef AMREX_USE_CUDA
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cusolverSp.h>
#include <cusparse.h>
#include <cuda_runtime_api.h>
#include <nvector/nvector_cuda.h>
#endif

#ifdef AMREX_USE_HIP
#include <nvector/nvector_hip.h>
#endif

#include <AMReX_Print.H>
/**********************************/

typedef struct ARKODEUserData 
{
    /* Checks */
    bool reactor_arkode_initialized;
    /* Base items */
    int ncells_d;
    int neqs_per_cell;
    int iverbose;
    int ireactor_type;
    double dt_save;

    double *rhoe_init_d = NULL;
    double *rhoesrc_ext_d = NULL;
    double *rYsrc_d = NULL;
    amrex::gpuStream_t stream;
    int nbBlocks;
    int nbThreads;
} *UserData;

int reactor_info(int cvode_iE, int Ncells);

static int cF_RHS(realtype t, N_Vector y_in, N_Vector ydot, void *user_data);

void reactor_close();

int react(realtype *rY_in, realtype *rY_src_in,
          realtype *rX_in, realtype *rX_src_in,
          realtype &dt_react, realtype &time,
          int cvode_iE, int Ncells, amrex::gpuStream_t stream);

int react(const amrex::Box& box,
          amrex::Array4<amrex::Real> const& rY_in,
          amrex::Array4<amrex::Real> const& rY_src_in,
          amrex::Array4<amrex::Real> const& T_in,
          amrex::Array4<amrex::Real> const& rEner_in,
          amrex::Array4<amrex::Real> const& rEner_src_in,
          amrex::Array4<amrex::Real> const& FC_in,
          amrex::Array4<int> const& mask,
          amrex::Real &dt_react,
          amrex::Real &time,
          const int &reactor_type,
          amrex::gpuStream_t stream);

AMREX_GPU_DEVICE inline void
fKernelSpec(int ncells, double dt_save, int reactor_type,
            realtype *yvec_d, realtype *ydot_d,
            double *rhoX_init, double *rhoXsrc_ext, double *rYs);

void SetTypValsODE(const std::vector<double>& ExtTypVals);

void SetTolFactODE(double relative_tol,double absolute_tol);

static int check_flag(void *flagvalue, const char *funcname, int opt);
