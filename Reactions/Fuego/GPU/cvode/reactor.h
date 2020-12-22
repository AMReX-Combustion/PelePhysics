#include <math.h>
#include <iostream>
#include <cassert>
#include <chrono>
#include <ctime>
#include <assert.h>

#include <cvode/cvode.h>               /* prototypes for CVODE fcts., consts.  */
#include <nvector/nvector_serial.h>    /* access to serial N_Vector            */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <nvector/nvector_cuda.h>
#include <sunmatrix/sunmatrix_cusparse.h>
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <sunlinsol/sunlinsol_spgmr.h> /* access to SPGMR SUNLinearSolver     */
#include <sunlinsol/sunlinsol_cusolversp_batchqr.h>
#include <cvode/cvode_direct.h>        /* access to CVDls interface            */
#include <sundials/sundials_types.h>   /* defs. of realtype, sunindextype      */
#include <sundials/sundials_math.h>

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
    int ireactor_type;
    /* Are we using a AJ */
    int ianalytical_jacobian;
    /* Are we using a IS or DS */
    int isolve_type;
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
    SUNMatrix R = NULL;
    /* CUDA cusolver */
    void *buffer_qr = NULL;
    csrqrInfo_t info;
    cusparseMatDescr_t descrA;
    cusolverSpHandle_t cusolverHandle;
    cusparseHandle_t cuSPHandle;
    cudaStream_t stream;
    int nbBlocks;
    int nbThreads;
}* UserData;

/* Functions Called by the Solver */
static int cF_RHS(realtype t, N_Vector y_in, N_Vector ydot, void *user_data);

/**********************************/
/* Functions Called by the main program */
int reactor_info(int cvode_iE, int Ncells);

int react(realtype *rY_in, realtype *rY_src_in,
          realtype *rX_in, realtype *rX_src_in,
          realtype &dt_react, realtype &time,
          int cvode_iE, int Ncells,
          cudaStream_t stream);

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
          const int &cvode_iE,
          cudaStream_t stream);

static int Precond(realtype tn, N_Vector u, N_Vector fu, booleantype jok,
                   booleantype *jcurPtr, realtype gamma, void *user_data);

static int PSolve(realtype tn, N_Vector u, N_Vector fu, N_Vector r, 
                  N_Vector z, realtype gamma, realtype delta, int lr, void *user_data);

static int cJac(realtype t, N_Vector y, N_Vector fy, SUNMatrix J,
                void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

void reactor_close();

/**********************************/
/* Additional useful functions */

static int check_flag(void *flagvalue, const char *funcname, int opt);

static void PrintFinalStats(void *cvode_mem);

void SetTypValsODE(const std::vector<double>& ExtTypVals);

void SetTolFactODE(double relative_tol,double absolute_tol);

/**********************************/
/* Device crap               */

// RHS kernel
AMREX_GPU_DEVICE
inline
void
fKernelSpec(int ncells, void *user_data, 
            realtype *yvec_d, realtype *ydot_d,  
            double *rhoX_init, double *rhoXsrc_ext, double *rYs);

// JACOBIANS
AMREX_GPU_DEVICE
inline
void 
fKernelComputeallAJ(int ncells, void *user_data, realtype *u_d, realtype *csr_val);

AMREX_GPU_DEVICE
inline
void 
fKernelComputeAJsys(int ncells, void *user_data, realtype *u_d, realtype *csr_val);

AMREX_GPU_DEVICE
inline
void 
fKernelComputeAJchem(int ncells, void *user_data, realtype *u_d, realtype *Jdata);

// CUSTOM
__global__
void 
fKernelDenseSolve(int ncells, realtype *x_d, realtype *b_d,
                  int subsys_size, int subsys_nnz, realtype *csr_val);

/**********************************/
/* Custom solver stuff */
struct _SUNLinearSolverContent_Dense_custom {
    sunindextype       last_flag;
    int                nsubsys;       /* number of subsystems */
    int                subsys_size;   /* size of each subsystem */
    int                subsys_nnz;
    int                nbBlocks;
    int                nbThreads;
    cudaStream_t       stream;
};

typedef struct _SUNLinearSolverContent_Dense_custom *SUNLinearSolverContent_Dense_custom; 

SUNLinearSolver SUNLinSol_dense_custom(N_Vector y, SUNMatrix A, 
                                       cudaStream_t stream);

SUNLinearSolver_Type SUNLinSolGetType_Dense_custom(SUNLinearSolver S); 

int SUNLinSolSolve_Dense_custom(SUNLinearSolver S, SUNMatrix A, N_Vector x,
                                N_Vector b, realtype tol);

int SUNLinSolSetup_Dense_custom(SUNLinearSolver S, SUNMatrix A);

int SUNLinSolFree_Dense_custom(SUNLinearSolver S);




