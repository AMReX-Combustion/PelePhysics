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

static int Precond(realtype tn, N_Vector u, N_Vector fu, booleantype jok,
		booleantype *jcurPtr, realtype gamma, void *user_data);

static int PSolve(realtype tn, N_Vector u, N_Vector fu, N_Vector r, 
		N_Vector z, realtype gamma, realtype delta, int lr, void *user_data);

void reactor_close();

/**********************************/
/* Additional useful functions */
//static void initialize_chemistry_device(UserData user_data);

static int check_flag(void *flagvalue, const char *funcname, int opt);

static void PrintFinalStats(void *cvode_mem);

/**********************************/
/* Device crap               */

/* Main Kernel fct called in solver RHS */
__global__ void fKernelSpec(void *user_data, 
		            realtype *yvec_d, realtype *ydot_d,  
		            double *rhoX_init, double *rhoXsrc_ext, double *rYs);


__global__ void fKernelJacCSR(realtype t, void *user_data,
                                          realtype *yvec_d, realtype *ydot_d,
                                          realtype* csr_jac, 
                                          const int size, const int nnz, 
                                          const int nbatched);

__global__ void fKernelComputeAJ(void *user_data, realtype *u_d, realtype *udot_d, realtype *csr_val);

__global__ void fKernelFillJB(void *user_data, realtype *gamma);


/* FROM FUEGO */
/*save inv molecular weights into array */
__device__ void imolecularWeight_d(double * iwt);
/* Returns R, Rc, Patm */
__device__ void ckrp_d( double * ru, double * ruc, double * pa);
/*Compute P = rhoRT/W(x) */
__device__ void ckpx_d(double * rho, double * T, double * x, double * P);
/*Compute P = rhoRT/W(y) */
__device__ void ckpy_d(double * rho, double * T, double * y_wk, double * P);
/*Compute rho = P*W(y)/RT */
__device__ void ckrhoy_d(double * P, double * T, double * y_wk, double * rho);
/*convert y[species] (mass fracs) to c[species] (molar conc) */
__device__ void ckytcr_d(double * rho, double * T, double * y_wk, double * c);
/*Returns the specific heats at constant volume */
/*in mass units (Eq. 29) */
__device__ void ckcvms_d(double * T, double * cvms);
/*Returns the specific heats at constant pressure */
/*in mass units (Eq. 26) */
__device__ void ckcpms_d(double * T, double * cpms);
/*Returns internal energy in mass units (Eq 30.) */
__device__ void ckums_d(double * T, double * ums);
/*Returns enthalpy in mass units (Eq 27.) */
__device__ void ckhms_d(double * T, double * hms);
/*Returns the mean specific heat at CP (Eq. 34) */
__device__ void ckcpbs_d(double * T, double * y_wk, double * cpbs);
/*Returns the mean specific heat at CV (Eq. 36) */
__device__ void ckcvbs_d(double * T, double * y_wk, double * cvbs);
/*Returns mean enthalpy of mixture in mass units */
__device__ void ckhbms_d(double * T, double * y_wk, double * hbms);
/*get mean internal energy in mass units */
__device__ void ckubms_d(double * T, double * y_wk, double * ubms);
/*compute the production rate for each species */
__device__ void ckwc_d(double * T, double * C, double * wdot);
/*compute the production rate for each species */
__device__ void productionRate_d(double * wdot, double * sc, double T);
__device__ void comp_qfqr_d(double *  qf, double * qr, double * sc, double * tc, double invT);
/*compute the production rate for each species */
__device__ void dwdot_d(double * J, double * sc, double * Tp, int * consP);
/*compute the reaction Jacobian */
__device__ void ajacobian_d(double * J, double * sc, double T, int consP);
/*compute the g/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
__device__ void gibbs_d(double * species, double *  tc);
/*compute Cv/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
__device__ void cv_R_d(double * species, double *  tc);
/*compute Cp/R at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
__device__ void cp_R_d(double * species, double *  tc);
/*compute d(Cp/R)/dT and d(Cv/R)/dT at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
__device__ void dcvpRdT_d(double * species, double *  tc);
/*compute the e/(RT) at the given temperature */
/*tc contains precomputed powers of T, tc[0] = log(T) */
__device__ void speciesInternalEnergy_d(double * species, double *  tc);
/*compute the h/(RT) at the given temperature (Eq 20) */
/*tc contains precomputed powers of T, tc[0] = log(T) */
__device__ void speciesEnthalpy_d(double * species, double *  tc);
/*save molecular weights into array */
__device__ void molecularWeight_d(double * wt);
/* get temperature given internal energy in mass units and mass fracs */
__device__ void get_t_given_ey_d_(double * e, double * y_wk, double * t, int * ierr);
/* get temperature given enthalpy in mass units and mass fracs */
__device__ void get_t_given_hy_d_(double * h, double * y_wk, double * t, int * ierr);
