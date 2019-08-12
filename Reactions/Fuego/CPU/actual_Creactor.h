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
#include <cvode/cvode_impl.h>

#ifdef USE_KLU 
#include "klu.h"
#include <sunlinsol/sunlinsol_klu.h>
#include <sundials/sundials_sparse.h>
#include <sunmatrix/sunmatrix_sparse.h>
#endif


#include <AMReX_Print.H>
#include <Fuego_EOS.H>

/**********************************/

typedef struct {
  realtype **(**P), **(**Jbd);
  sunindextype *(**pivot);
#ifdef USE_KLU 
  int NNZ; 
  SUNMatrix *PS;
  realtype **Jdata = NULL;
  int **rowVals = NULL;
  int **colPtrs = NULL;
  klu_common *Common;
  klu_symbolic **Symbolic;
  klu_numeric **Numeric;
#endif
} *UserData;


/**********************************/
/* Functions Called by the Solver */
static int cF_RHS(realtype t, N_Vector y_in, N_Vector ydot, void *user_data);

static int cJac(realtype tn, N_Vector y, N_Vector fy, SUNMatrix J,
		void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

#ifdef USE_KLU 
static int cJac_KLU(realtype tn, N_Vector y, N_Vector fy, SUNMatrix J,
		void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
#endif

//static int jtv(N_Vector v, N_Vector Jv, realtype t, N_Vector u,
//		N_Vector fu, void *user_data, N_Vector tmp);

static int PSolve(realtype tn, N_Vector u, N_Vector fu, N_Vector r, 
		N_Vector z, realtype gamma, realtype delta, int lr, void *user_data);

#ifdef USE_KLU 
static int PSolve_sparse(realtype tn, N_Vector u, N_Vector fu, N_Vector r, 
		N_Vector z, realtype gamma, realtype delta, int lr, void *user_data);
static int Precond_sparse(realtype tn, N_Vector u, N_Vector fu, booleantype jok,
		booleantype *jcurPtr, realtype gamma, void *user_data);
#endif

static int Precond(realtype tn, N_Vector u, N_Vector fu, booleantype jok,
		booleantype *jcurPtr, realtype gamma, void *user_data);


/**********************************/
/* Functions Called by the Program */
int reactor_init(const int* cvode_iE, const int* Ncells);

int react(realtype *rY_in, realtype *rY_src_in, 
		realtype *rX_in, realtype *rX_src_in, 
		realtype *P_in, realtype *dt_react, realtype *time, int *Init);

void reactor_close();


/**********************************/
/* Helper functions */
static int check_flag(void *flagvalue, const char *funcname, int opt);

static void PrintFinalStats(void *cvodeMem, realtype Temp);

static UserData AllocUserData(void);

static void FreeUserData(UserData data);

static void check_state(N_Vector yvec);


/**********************************/
/* Main Kernel fct called in solver RHS */
void fKernelSpec(realtype *dt, realtype *yvec_d, realtype *ydot_d,
		double *rhoX_init, double *rhoXsrc_ext, double *rYs);


