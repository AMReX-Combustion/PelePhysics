#include <math.h>
#include <iostream>
#include <cstring>
#include <chrono>

#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sunlinsol/sunlinsol_spgmr.h> /* access to SPGMR SUNLinearSolver     */
#include <cvode/cvode_direct.h>        /* access to CVDls interface            */
#include <cvode/cvode_spils.h>         /* access to CVSpils interface */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#include <sundials/sundials_math.h>
#include <sunmatrix/sunmatrix_sparse.h>

#ifdef USE_KLU_PP 
#include "klu.h"
#include <sunlinsol/sunlinsol_klu.h>
#endif

#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>

#include <AMReX_Print.H>
#include <EOS.H>
/**********************************/

typedef struct {
      /* hacks */
      bool FirstTimePrecond;
      /* Checks */
      bool reactor_cvode_initialized;
      bool actual_ok_to_react;
      /* Base items */
      int ncells;
      int iverbose;
      int isolve_type;
      int ianalytical_jacobian;
      int ireactor_type;
      int boxcell;
      /* external forcing */
      amrex::Gpu::ManagedVector<amrex::Real> Yvect_full; 
      amrex::Gpu::ManagedVector<amrex::Real> rhoX_init;
      amrex::Gpu::ManagedVector<amrex::Real> rhoXsrc_ext;
      amrex::Gpu::ManagedVector<amrex::Real> rYsrc;
      amrex::Gpu::ManagedVector<int> FCunt;
      /* Options */
      int NNZ; 
      /* Sparse Matrices for KLU-related solve */
      SUNMatrix *PS;
      /* SUNSparseMatrix_Data */
      realtype **Jdata  = NULL;
      /* SUNSparseMatrix_IndexValues */
      int **rowVals     = NULL;
      int **rowPtrs     = NULL;
      /* SUNSparseMatrix_IndexPointers */
      int **colPtrs     = NULL;
      int **colVals     = NULL;
      /* Holder for sparse matrix in Fuego fetches */
      int *indx = NULL;
      realtype **JSPSmat = NULL;
#ifdef USE_KLU_PP 
      /* KLU objects */
      klu_common *Common;
      klu_symbolic **Symbolic;
      klu_numeric **Numeric;
#else
      realtype **(**Jbd);
      realtype **(**P);
      sunindextype *(**pivot);
#endif
      /* Sparse custom */
      SUNMatrix PSc;
      int *colVals_c;
      int *rowPtrs_c;
} *UserData;


/**********************************/
/* Functions Called by the Solver */
int cF_RHS(realtype t, N_Vector y_in, N_Vector ydot, void *user_data);

int cJac(realtype tn, N_Vector y, N_Vector fy, SUNMatrix J,
		void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

int cJac_sps(realtype tn, N_Vector y, N_Vector fy, SUNMatrix J,
		void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

int PSolve_custom(realtype tn, N_Vector u, N_Vector fu, N_Vector r, 
		N_Vector z, realtype gamma, realtype delta, int lr, void *user_data);

int Precond_custom(realtype tn, N_Vector u, N_Vector fu, booleantype jok,
		booleantype *jcurPtr, realtype gamma, void *user_data);

#ifdef USE_KLU_PP 
int cJac_KLU(realtype tn, N_Vector y, N_Vector fy, SUNMatrix J,
		void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

int PSolve_sparse(realtype tn, N_Vector u, N_Vector fu, N_Vector r, 
		N_Vector z, realtype gamma, realtype delta, int lr, void *user_data);
int Precond_sparse(realtype tn, N_Vector u, N_Vector fu, booleantype jok,
		booleantype *jcurPtr, realtype gamma, void *user_data);
#else
int PSolve(realtype tn, N_Vector u, N_Vector fu, N_Vector r, 
		N_Vector z, realtype gamma, realtype delta, int lr, void *user_data);

int Precond(realtype tn, N_Vector u, N_Vector fu, booleantype jok,
		booleantype *jcurPtr, realtype gamma, void *user_data);
#endif



/**********************************/
/* Functions Called by the Program */
int reactor_init(const int cvode_iE, const int Ncells);

int react(const amrex::Box& box,
          amrex::Array4<amrex::Real> const& rY_in,
          amrex::Array4<amrex::Real> const& rY_src_in, 
          amrex::Array4<amrex::Real> const& T_in,
          amrex::Array4<amrex::Real> const& rEner_in,  
          amrex::Array4<amrex::Real> const& rEner_src_in,
          amrex::Array4<amrex::Real> const& FC_in,
          amrex::Array4<int> const& mask,
          int box_ncells,
          amrex::Real &dt_react,
          amrex::Real &time);

int react(realtype *rY_in, realtype *rY_src_in, 
	      realtype *rX_in, realtype *rX_src_in, 
	      realtype &dt_react, realtype &time);

void reactor_close();


/**********************************/
/* Helper functions */
int check_flag(void *flagvalue, const char *funcname, int opt);

void PrintFinalStats(void *cvodeMem, realtype Temp);

UserData AllocUserData(int iE, int num_cells);

void FreeUserData(UserData data);

void check_state(N_Vector yvec);

void SetTypValsODE(std::vector<double> ExtTypVals);

void SetTolFactODE(double relative_tol,double absolute_tol);

void ReSetTolODE();

/**********************************/
/* Main Kernel fct called in solver RHS */
void fKernelSpec(realtype *dt, realtype *yvec_d, realtype *ydot_d,
		void *user_data);

/**********************************/
/* custom solver */
struct _SUNLinearSolverContent_Sparse_custom {
	sunindextype       last_flag;
	int                reactor_type;
	int                nsubsys;       /* number of subsystems */
	int                subsys_size;   /* size of each subsystem */
	int                subsys_nnz;

};

typedef struct _SUNLinearSolverContent_Sparse_custom *SUNLinearSolverContent_Sparse_custom; 

SUNLinearSolver SUNLinSol_sparse_custom(N_Vector y, SUNMatrix A, int reactor_type,
		int nsubsys, int subsys_size, int subsys_nnz);

SUNLinearSolver_Type SUNLinSolGetType_Sparse_custom(SUNLinearSolver S); 

int SUNLinSolSolve_Sparse_custom(SUNLinearSolver S, SUNMatrix A, N_Vector x,
		N_Vector b, realtype tol);

