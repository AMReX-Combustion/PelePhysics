#include <math.h>
#include <iostream>
#include <cstring>
#include <chrono>

#include <arkode/arkode_arkstep.h>
#include <arkode/arkode_erkstep.h>

#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sunnonlinsol/sunnonlinsol_newton.h>
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#include <sundials/sundials_math.h>

#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

#include <AMReX_Print.H>
#include <EOS.H>
/**********************************/

typedef struct {
      /* Checks */
      bool reactor_arkode_initialized;
      /* Base items */
      int ncells;
      int iverbose;
      int ianalytical_jacobian;
      int ireactor_type;
      int iimplicit_solve = 0;
      int iuse_erkstep = 0;
      int boxcell = 0;

      amrex::Real *Yvect_full = NULL;
      amrex::Real *rhoX_init = NULL;
      amrex::Real *rhoXsrc_ext = NULL;
      amrex::Real *rYsrc = NULL;
      int *FCunt = NULL;
      int *mask = NULL;
} *UserData;


/**********************************/
/* Functions Called by the Solver */
int cF_RHS(realtype t, N_Vector y_in, N_Vector ydot, void *user_data);

int cJac(realtype tn, N_Vector y, N_Vector fy, SUNMatrix J,
		void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

/**********************************/
/* Functions Called by the Program */
int reactor_init(int cvode_iE, int Ncells);

int react(realtype *rY_in, realtype *rY_src_in,
        realtype *rX_in, realtype *rX_src_in,
        realtype &dt_react, realtype &time);

int react(const amrex::Box& box,
        amrex::Array4<amrex::Real> const& rY_in,
        amrex::Array4<amrex::Real> const& rY_src_in,
        amrex::Array4<amrex::Real> const& T_in,
        amrex::Array4<amrex::Real> const& rEner_in,
        amrex::Array4<amrex::Real> const& rEner_src_in,
        amrex::Array4<amrex::Real> const& FC_in,
        amrex::Array4<int> const& mask,
        amrex::Real &dt_react,
        amrex::Real &time); 

void reactor_close();


/**********************************/
/* Helper functions */
int check_flag(void *flagvalue, const char *funcname, int opt);

void PrintFinalStats(void *arkodeMem, realtype Temp);

UserData AllocUserData(int iE, int num_cells);

void FreeUserData(UserData data);

void SetTypValsODE(const std::vector<double>& ExtTypVals);

void SetTolFactODE(double relative_tol,double absolute_tol);

void ReSetTolODE();

/**********************************/
/* Main Kernel fct called in solver RHS */
void fKernelSpec(realtype *dt, realtype *yvec_d, realtype *ydot_d,
        void *user_data);


