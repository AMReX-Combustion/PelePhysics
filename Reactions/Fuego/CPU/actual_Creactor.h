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

#ifdef USE_KLU_PP 
#include "klu.h"
#include <sunlinsol/sunlinsol_klu.h>
#include <sunmatrix/sunmatrix_sparse.h>
#endif

#include <AMReX_Print.H>
#include <Fuego_EOS.H>
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
      int iDense_Creact;
      int iJac_Creact;
      int iE_Creact;
#ifdef USE_KLU_PP 
      int NNZ; 
      /* Sparse Matrices for KLU-related solve */
      SUNMatrix *PS;
      /* SUNSparseMatrix_Data */
      realtype **Jdata  = NULL;
      /* SUNSparseMatrix_IndexValues */
      int **rowVals     = NULL;
      /* SUNSparseMatrix_IndexPointers */
      int **colPtrs     = NULL;
      /* Holder for sparse matrix in Fuego fetches */
      int *indx = NULL;
      realtype **JSPSmat = NULL;
      /* KLU objects */
      klu_common *Common;
      klu_symbolic **Symbolic;
      klu_numeric **Numeric;
#else
      realtype **(**Jbd);
      realtype **(**P);
      sunindextype *(**pivot);
#endif
} *UserData;


/**********************************/
/* Functions Called by the Solver */
int cF_RHS(realtype t, N_Vector y_in, N_Vector ydot, void *user_data);

int cJac(realtype tn, N_Vector y, N_Vector fy, SUNMatrix J,
		void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

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
extern "C"
{
    int reactor_init(const int* cvode_iE, const int* Ncells);

    int react(realtype *rY_in, realtype *rY_src_in, 
		realtype *rX_in, realtype *rX_src_in, 
		realtype *dt_react, realtype *time);

    void reactor_close();
}



/**********************************/
/* Helper functions */
int check_flag(void *flagvalue, const char *funcname, int opt);

void PrintFinalStats(void *cvodeMem, realtype Temp);

UserData AllocUserData(int iE, int num_cells);

void FreeUserData(UserData data);

void check_state(N_Vector yvec);


/**********************************/
/* Main Kernel fct called in solver RHS */
void fKernelSpec(realtype *dt, realtype *yvec_d, realtype *ydot_d,
		void *user_data);


