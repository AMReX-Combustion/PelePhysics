#include <math.h>
#include <iostream>
#include <cassert>
#include <chrono>
#include <ctime>
#include <assert.h>

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

#include <nvector/nvector_cuda.h>

//#include <cusolver/cvode_cusolver_spqr.h>

#include <sundials/sundials_sparse.h>
#include <sunmatrix/sunmatrix_sparse.h>
#include <AMReX_Print.H>

#include <cuda_runtime.h>
#include <cublas_v2.h>
#include <cusolverSp.h>
#include <cusparse.h>
#include <cuda_runtime_api.h>

/**********************************/
typedef struct CVodeUserData {
    /* LS gamma */
    double gamma_d;
    /* dt on device */
    double dt_save;
    /* nb of cells to integrate */
    int ncells_d[1]; 
    /* nb of eq per cell */
    int neqs_per_cell[1];
    /* HP/UV react */
    int flagP;
    /* Are we using a AJ */
    int iJac;
    /* energy related variables */
    double *rhoe_init = NULL;
    double *rhoesrc_ext = NULL;
    double *rYsrc = NULL;
    /* verbose level */
    int iverbose;
    // Sparse
    /* Precond sparse stuff */
    int NNZ; 
    int* csr_row_count_d;
    int* csr_col_index_d;
    double* csr_val_d;
    double* csr_jac_d;
    /* CUDA cusolver */
    void *buffer_qr = NULL;
    csrqrInfo_t info;
    cusparseMatDescr_t descrA;
    cusolverSpHandle_t cusolverHandle;
    cudaStream_t stream;
}* UserData;

/* Functions Called by the Solver */
static int cF_RHS(realtype t, N_Vector y_in, N_Vector ydot, void *user_data);

/**********************************/
/* Functions Called by the main program */
int reactor_info(const int* cvode_iE, const int* Ncells); 

int react(realtype *rY_in, realtype *rY_src_in, 
		realtype *rX_in, realtype *rX_src_in, 
		realtype *P_in, 
		realtype *dt_react, realtype *time, int *Init,
                const int* cvode_iE, const int* Ncells, cudaStream_t stream);

//static int Precond(realtype tn, N_Vector u, N_Vector fu, booleantype jok,
//		booleantype *jcurPtr, realtype gamma, void *user_data);

//static int PSolve(realtype tn, N_Vector u, N_Vector fu, N_Vector r, 
//		N_Vector z, realtype gamma, realtype delta, int lr, void *user_data);

void reactor_close();

/**********************************/
/* Additional useful functions */
//static void initialize_chemistry_device(UserData user_data);

static int check_flag(void *flagvalue, const char *funcname, int opt);

static void PrintFinalStats(void *cvode_mem);

/**********************************/
/* Device crap               */

/* Main Kernel fct called in solver RHS */
AMREX_GPU_DEVICE
inline
void
fKernelSpec(int ncells, void *user_data, 
		            realtype *yvec_d, realtype *ydot_d,  
		            double *rhoX_init, double *rhoXsrc_ext, double *rYs);


//__global__ void fKernelJacCSR(realtype t, void *user_data,
//                                          realtype *yvec_d, realtype *ydot_d,
//                                          realtype* csr_jac, 
//                                          const int size, const int nnz, 
//                                          const int nbatched);

//__global__ void fKernelComputeAJ(void *user_data, realtype *u_d, realtype *udot_d, realtype *csr_val);

//__global__ void fKernelFillJB(void *user_data, realtype *gamma);


