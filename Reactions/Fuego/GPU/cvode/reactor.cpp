#include <reactor.h> 
#include <AMReX_ParmParse.H>
#include <chemistry_file.H>
#include <EOS.H>
#include "mechanism.h"
#include <AMReX_Gpu.H>

using namespace amrex;

#define SUN_CUSP_CONTENT(S)        ( (SUNLinearSolverContent_Dense_custom)(S->content) )
#define SUN_CUSP_SUBSYS_SIZE(S)    ( SUN_CUSP_CONTENT(S)->subsys_size )
#define SUN_CUSP_NUM_SUBSYS(S)     ( SUN_CUSP_CONTENT(S)->nsubsys )
#define SUN_CUSP_SUBSYS_NNZ(S)     ( SUN_CUSP_CONTENT(S)->subsys_nnz ) 

#define SUN_CUSP_DCOLIND(S)        ( SUN_CUSP_CONTENT(S)->d_colind ) 
#define SUN_CUSP_DROWPTR(S)        ( SUN_CUSP_CONTENT(S)->d_rowptr )
#define SUN_CUSP_DVALUES(S)        ( SUN_CUSP_CONTENT(S)->d_values )

#define SUN_CUSP_LASTFLAG(S)       ( SUN_CUSP_CONTENT(S)->last_flag )
#define SUN_CUSP_STREAM(S)         ( SUN_CUSP_CONTENT(S)->stream )
#define SUN_CUSP_NBLOCK(S)         ( SUN_CUSP_CONTENT(S)->nbBlocks )
#define SUN_CUSP_NTHREAD(S)        ( SUN_CUSP_CONTENT(S)->nbThreads )

/**********************************/
/* Global Variables */

AMREX_GPU_DEVICE_MANAGED int sparse_solve          = 1;
AMREX_GPU_DEVICE_MANAGED int sparse_cusolver_solve = 5;
AMREX_GPU_DEVICE_MANAGED int iterative_gmres_solve = 99;
AMREX_GPU_DEVICE_MANAGED int eint_rho = 1; // in/out = rhoE/rhoY
AMREX_GPU_DEVICE_MANAGED int enth_rho = 2; // in/out = rhoH/rhoY
/**********************************/

/**********************************/
/* Infos to print once */
int reactor_info(const int* reactor_type,const int* Ncells){ 

    /* ParmParse from the inputs file */ 
    amrex::ParmParse pp("ode");
    int ianalytical_jacobian = 0;
    pp.query("analytical_jacobian",ianalytical_jacobian);
    int iverbose = 1;
    pp.query("verbose",iverbose);

    if (iverbose > 0) {
        amrex::Print() << "Nb of spec in mech is " << NUM_SPECIES << "\n";    

        amrex::Print() << "Ncells in one solve is " << *Ncells << "\n";
    }

    std::string  solve_type_str = "none"; 
    amrex::ParmParse ppcv("cvode");
    ppcv.query("solve_type", solve_type_str);
    /* options are:
    sparse_solve          = 1;
    sparse_cusolver_solve = 5;
    iterative_gmres_solve = 99;
    */
    int isolve_type = iterative_gmres_solve;
    if (solve_type_str == "sparse_custom") {
        isolve_type = sparse_solve; 
    } else if (solve_type_str == "sparse") {
        isolve_type = sparse_cusolver_solve;
    } else if (solve_type_str == "GMRES") {
        isolve_type = iterative_gmres_solve; 
    }

    /* Checks */
    if (isolve_type == iterative_gmres_solve) {
        if (ianalytical_jacobian == 1) { 
            if (iverbose > 0) {
                amrex::Print() <<"Using an Iterative GMRES Solver with sparse simplified preconditioning \n";
            }
            int nJdata;
            int HP;
            if (*reactor_type == eint_rho) {
                HP = 0;
            } else {
                HP = 1;
            }
            /* Precond data */ 
            SPARSITY_INFO_SYST_SIMPLIFIED(&nJdata,&HP);
            if (iverbose > 0) {
                amrex::Print() << "--> SPARSE Preconditioner -- non zero entries: " << nJdata << ", which represents "<< nJdata/float((NUM_SPECIES+1) * (NUM_SPECIES+1)) *100.0 <<" % fill-in pattern\n";
            }         
        } else {
            if (iverbose > 0) {
                amrex::Print() <<"Using an Iterative GMRES Solver without preconditionning \n";
            }
        }

    } else if (isolve_type == sparse_solve) {
        if (ianalytical_jacobian == 1) {
            if (iverbose > 0) {
                amrex::Print() <<"Using a custom Sparse Direct solver (with Analytical Jacobian) \n";
            }
            int nJdata;
            int HP;
            if (*reactor_type == eint_rho) {
                HP = 0;
            } else {
                HP = 1;
            }
            /* Jac data */ 
            SPARSITY_INFO_SYST(&nJdata,&HP,*Ncells);
            if (iverbose > 0) {
                amrex::Print() << "--> SPARSE Solver -- non zero entries: " << nJdata << ", which represents "<< nJdata/float(*Ncells * (NUM_SPECIES+1) * (NUM_SPECIES+1)) * 100.0 <<" % fill-in pattern\n";
            }         
        } else {
            amrex::Abort("\n--> When using a custom direct sparse solver, specify an AJ \n");
        }

    } else if (isolve_type == sparse_cusolver_solve) {
        if (ianalytical_jacobian == 1) {
            amrex::Print() <<"Using a Sparse Direct Solver based on cuSolver \n";
            int nJdata;
            int HP;
            if (*reactor_type == eint_rho) {
                HP = 0;
            } else {
                HP = 1;
            }
            /* Jac data */ 
            SPARSITY_INFO_SYST(&nJdata,&HP,*Ncells);
            if (iverbose > 0) {
                amrex::Print() << "--> SPARSE Solver -- non zero entries: " << nJdata << ", which represents "<< nJdata/float(*Ncells * (NUM_SPECIES+1) * (NUM_SPECIES+1)) * 100.0 <<" % fill-in pattern\n";
            }         
        } else {
                amrex::Abort("\n--> When using a SDS, specify an AJ \n");
        }

    } else {
        amrex::Abort("\n--> Wrong choice of input parameters (cvode) \n");
    }

    amrex::Print() << "\n--> DONE WITH INITIALIZATION (GPU)" << *reactor_type << "\n";

    return(0);
}


/* Main routine for external looping */
int react(realtype *rY_in, realtype *rY_src_in, 
          realtype *rX_in, realtype *rX_src_in,
          realtype *dt_react, realtype *time,
          const int* reactor_type,const int* Ncells, cudaStream_t stream) {

    /* CVODE */
    N_Vector y         = NULL;
    SUNLinearSolver LS = NULL;
    SUNMatrix A        = NULL;
    void *cvode_mem    = NULL;
    /* Misc */
    int flag;
    int NCELLS, NEQ, neq_tot;
    realtype reltol;
    N_Vector atol;
    realtype *ratol;

    /* ParmParse from the inputs file */ 
    amrex::ParmParse pp("ode");
    int ianalytical_jacobian = 0;
    pp.query("analytical_jacobian",ianalytical_jacobian);
    int iverbose = 1;
    pp.query("verbose",iverbose);

    std::string  solve_type_str = "none"; 
    amrex::ParmParse ppcv("cvode");
    ppcv.query("solve_type", solve_type_str);
    /* options are:
    sparse_solve          = 1;
    sparse_cusolver_solve = 5;
    iterative_gmres_solve = 99;
    */
    int isolve_type = iterative_gmres_solve;
    if (solve_type_str == "sparse_custom") {
        isolve_type = sparse_solve; 
    } else if (solve_type_str == "sparse") {
        isolve_type = sparse_cusolver_solve;
    } else if (solve_type_str == "GMRES") {
        isolve_type = iterative_gmres_solve; 
    }

    NEQ = NUM_SPECIES;

    /* Args */
    NCELLS         = *Ncells;
    neq_tot        = (NEQ + 1) * NCELLS;

    /* User data */
    UserData user_data;
    BL_PROFILE_VAR("AllocsInCVODE", AllocsCVODE);
    cudaMallocManaged(&user_data, sizeof(struct CVodeUserData));
    BL_PROFILE_VAR_STOP(AllocsCVODE);
    user_data->ncells_d[0]          = NCELLS;
    user_data->neqs_per_cell[0]     = NEQ;
    user_data->ireactor_type        = *reactor_type; 
    user_data->ianalytical_jacobian = ianalytical_jacobian;
    user_data->isolve_type          = isolve_type; 
    user_data->iverbose             = iverbose;
    user_data->stream               = stream;
    user_data->nbBlocks             = std::max(1,NCELLS/32);
    user_data->nbThreads            = 32;

    if (user_data->ianalytical_jacobian == 1) { 
        int HP;
        if (user_data->ireactor_type == 1) {
            HP = 0;
        } else {
            HP = 1;
        }

        /* Find sparsity pattern to fill structure of sparse matrix */
        BL_PROFILE_VAR("SparsityFuegoStuff", SparsityStuff);
        BL_PROFILE_VAR_STOP(SparsityStuff);
        if (user_data->isolve_type == iterative_gmres_solve) {
            BL_PROFILE_VAR_START(SparsityStuff);
            SPARSITY_INFO_SYST_SIMPLIFIED(&(user_data->NNZ),&HP);
	    BL_PROFILE_VAR_STOP(SparsityStuff);

	    BL_PROFILE_VAR_START(AllocsCVODE);
            cudaMallocManaged(&(user_data->csr_row_count_d), (NEQ+2) * sizeof(int));
            cudaMallocManaged(&(user_data->csr_col_index_d), user_data->NNZ * sizeof(int));
            cudaMallocManaged(&(user_data->csr_jac_d), user_data->NNZ * NCELLS * sizeof(double));
            cudaMallocManaged(&(user_data->csr_val_d), user_data->NNZ * NCELLS * sizeof(double));
	    BL_PROFILE_VAR_STOP(AllocsCVODE);

            BL_PROFILE_VAR_START(SparsityStuff);
            SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(user_data->csr_col_index_d, user_data->csr_row_count_d, &HP,1);
	    BL_PROFILE_VAR_STOP(SparsityStuff);

        } else if (isolve_type == sparse_cusolver_solve) {
            BL_PROFILE_VAR_START(SparsityStuff);
            SPARSITY_INFO_SYST(&(user_data->NNZ),&HP,1);
	    BL_PROFILE_VAR_STOP(SparsityStuff);

	    BL_PROFILE_VAR_START(AllocsCVODE);
            cudaMallocManaged(&(user_data->csr_row_count_d), (NEQ+2) * sizeof(int));
            cudaMallocManaged(&(user_data->csr_col_index_d), user_data->NNZ * sizeof(int));
            cudaMallocManaged(&(user_data->csr_jac_d), user_data->NNZ * NCELLS * sizeof(double));
	    BL_PROFILE_VAR_STOP(AllocsCVODE);

            cusolverStatus_t cusolver_status = CUSOLVER_STATUS_SUCCESS;
            cusparseStatus_t cusparse_status = CUSPARSE_STATUS_SUCCESS;
            cudaError_t      cuda_status     = cudaSuccess;
	    int retval;

            cusolver_status = cusolverSpCreate(&(user_data->cusolverHandle));
            assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);

	    cusparse_status = cusparseCreate(&(user_data->cuSPHandle));
	    assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);

	    A = SUNMatrix_cuSparse_NewBlockCSR(NCELLS, (NEQ + 1), (NEQ + 1), user_data->NNZ, user_data->cuSPHandle);
	    if (check_flag((void *)A, "SUNMatrix_cuSparse_NewBlockCSR", 0)) return(1);

	    retval = SUNMatrix_cuSparse_SetFixedPattern(A, 1); 
	    if(check_flag(&retval, "SUNMatrix_cuSparse_SetFixedPattern", 1)) return(1);

            BL_PROFILE_VAR_START(SparsityStuff);
	    SPARSITY_PREPROC_SYST_CSR(user_data->csr_col_index_d, user_data->csr_row_count_d, &HP, 1, 0); 
            SUNMatrix_cuSparse_CopyToDevice(A, NULL, user_data->csr_row_count_d, user_data->csr_col_index_d);
            cuda_status = cudaDeviceSynchronize();
            assert(cuda_status == cudaSuccess);
	    //SPARSITY_PREPROC_CSR(SUNMatrix_cuSparse_IndexValues(A), SUNMatrix_cuSparse_IndexPointers(A), &HP, 1, 0); 
	    //SPARSITY_PREPROC_CSR(user_data->csr_col_index_d, user_data->csr_row_count_d, &HP, 1, 1); 
	    BL_PROFILE_VAR_STOP(SparsityStuff);
        } else {
            BL_PROFILE_VAR_START(SparsityStuff);
            SPARSITY_INFO_SYST(&(user_data->NNZ),&HP,1);
	    BL_PROFILE_VAR_STOP(SparsityStuff);

	    BL_PROFILE_VAR_START(AllocsCVODE);
            cudaMallocManaged(&(user_data->csr_row_count_d), (NEQ+2) * sizeof(int));
            cudaMallocManaged(&(user_data->csr_col_index_d), user_data->NNZ * sizeof(int));
            cudaMallocManaged(&(user_data->csr_jac_d), user_data->NNZ * NCELLS * sizeof(double));
	    BL_PROFILE_VAR_STOP(AllocsCVODE);

            A = SUNSparseMatrix(neq_tot, neq_tot, user_data->NNZ * NCELLS, CSR_MAT);
            if (check_flag((void *)A, "SUNSparseMatrix", 0)) return(1);

            BL_PROFILE_VAR_START(SparsityStuff);
	    SPARSITY_PREPROC_SYST_CSR(user_data->csr_col_index_d, user_data->csr_row_count_d, &HP, 1, 1); 
	    SPARSITY_PREPROC_SYST_CSR( (int*)SUNSparseMatrix_IndexValues(A), (int*)SUNSparseMatrix_IndexPointers(A), &HP, NCELLS, 0); 
	    BL_PROFILE_VAR_STOP(SparsityStuff);
        }

    }


    // Create Sparse batch QR solver
    // qr info and matrix descriptor
    BL_PROFILE_VAR("CuSolverInit", CuSolverInit);
    if (user_data->isolve_type == iterative_gmres_solve) {
        size_t workspaceInBytes, internalDataInBytes;
        cusolverStatus_t cusolver_status = CUSOLVER_STATUS_SUCCESS;
        cusparseStatus_t cusparse_status = CUSPARSE_STATUS_SUCCESS;

        workspaceInBytes = 0;
        internalDataInBytes = 0;

        cusolver_status = cusolverSpCreate(&(user_data->cusolverHandle));
        assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);

        cusparse_status = cusparseCreateMatDescr(&(user_data->descrA)); 
        assert(cusparse_status == CUSPARSE_STATUS_SUCCESS);

        cusparse_status = cusparseSetMatType(user_data->descrA, CUSPARSE_MATRIX_TYPE_GENERAL);
        assert(cusparse_status == CUSPARSE_STATUS_SUCCESS);

        cusparse_status = cusparseSetMatIndexBase(user_data->descrA, CUSPARSE_INDEX_BASE_ONE);
        assert(cusparse_status == CUSPARSE_STATUS_SUCCESS);

        cusolver_status = cusolverSpCreateCsrqrInfo(&(user_data->info));
        assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);

        // symbolic analysis
        cusolver_status = cusolverSpXcsrqrAnalysisBatched(user_data->cusolverHandle,
                                   NEQ+1, // size per subsystem
                                   NEQ+1, // size per subsystem
                                   user_data->NNZ,
                                   user_data->descrA,
                                   user_data->csr_row_count_d,
                                   user_data->csr_col_index_d,
                                   user_data->info);
        assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);
           
        /*
        size_t free_mem = 0;
        size_t total_mem = 0;
        cudaStat1 = cudaMemGetInfo( &free_mem, &total_mem );
        assert( cudaSuccess == cudaStat1 );
        std::cout<<"(AFTER SA) Free: "<< free_mem<< " Tot: "<<total_mem<<std::endl;
        */

        // allocate working space 
        cusolver_status = cusolverSpDcsrqrBufferInfoBatched(user_data->cusolverHandle,
                                                  NEQ+1, // size per subsystem
                                                  NEQ+1, // size per subsystem
                                                  user_data->NNZ,
                                                  user_data->descrA,
                                                  user_data->csr_val_d,
                                                  user_data->csr_row_count_d,
                                                  user_data->csr_col_index_d,
                                                  NCELLS,
                                                  user_data->info,
                                                  &internalDataInBytes,
                                                  &workspaceInBytes);
        assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);

	cudaError_t cudaStat1            = cudaSuccess;
	cudaStat1 = cudaMalloc((void**)&(user_data->buffer_qr), workspaceInBytes);
	assert(cudaStat1 == cudaSuccess);
    }
    BL_PROFILE_VAR_STOP(CuSolverInit);

    /* Definition of main vector */
    y = N_VNewManaged_Cuda(neq_tot);
    if(check_flag((void*)y, "N_VNewManaged_Cuda", 0)) return(1);

    /* Use a non-default cuda stream for kernel execution */
    N_VSetCudaStream_Cuda(y, &stream);

    /* Call CVodeCreate to create the solver memory and specify the
     * Backward Differentiation Formula and the use of a Newton iteration */
    cvode_mem = CVodeCreate(CV_BDF);
    if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);

    flag = CVodeSetUserData(cvode_mem, static_cast<void*>(user_data));

    BL_PROFILE_VAR_START(AllocsCVODE);
    /* Define vectors to be used later in creact */
    cudaMalloc(&(user_data->rhoe_init), NCELLS*sizeof(double));
    cudaMalloc(&(user_data->rhoesrc_ext), NCELLS*sizeof(double));
    cudaMalloc(&(user_data->rYsrc), (NCELLS*NEQ)*sizeof(double));
    BL_PROFILE_VAR_STOP(AllocsCVODE);

    /* Get Device MemCpy of in arrays */
    /* Get Device pointer of solution vector */
    realtype *yvec_d      = N_VGetDeviceArrayPointer_Cuda(y);
    BL_PROFILE_VAR("AsyncCpy", AsyncCpy);
    // rhoY,T
    cudaMemcpy(yvec_d, rY_in, sizeof(realtype) * ((NEQ+1)*NCELLS), cudaMemcpyHostToDevice);
    // rhoY_src_ext
    cudaMemcpy(user_data->rYsrc, rY_src_in, (NEQ*NCELLS)*sizeof(double), cudaMemcpyHostToDevice);
    // rhoE/rhoH
    cudaMemcpy(user_data->rhoe_init, rX_in, sizeof(realtype) * NCELLS, cudaMemcpyHostToDevice);
    cudaMemcpy(user_data->rhoesrc_ext, rX_src_in, sizeof(realtype) * NCELLS, cudaMemcpyHostToDevice);
    BL_PROFILE_VAR_STOP(AsyncCpy);

    realtype time_init, time_out ;
    time_init = *time;
    time_out  = *time + *dt_react;

    /* Call CVodeInit to initialize the integrator memory and specify the
     * user's right hand side function, the inital time, and 
     * initial dependent variable vector y. */
    flag = CVodeInit(cvode_mem, cF_RHS, time_init, y);
    if (check_flag(&flag, "CVodeInit", 1)) return(1);
    
    /* Definition of tolerances: one for each species */
    reltol = 1.0e-10;
    atol  = N_VNew_Cuda(neq_tot);
    ratol = N_VGetHostArrayPointer_Cuda(atol);
    for (int i=0; i<neq_tot; i++) {
        ratol[i] = 1.0e-10;
    }
    N_VCopyToDevice_Cuda(atol);
    /* Call CVodeSVtolerances to specify the scalar relative tolerance
     * and vector absolute tolerances */
    flag = CVodeSVtolerances(cvode_mem, reltol, atol);
    if (check_flag(&flag, "CVodeSVtolerances", 1)) return(1);

    /* Create the linear solver object */
    if (user_data->isolve_type == iterative_gmres_solve) {
        if (user_data->ianalytical_jacobian == 0) { 
            LS = SUNSPGMR(y, PREC_NONE, 0);
            if(check_flag((void *)LS, "SUNDenseLinearSolver", 0)) return(1);
        } else { 
            LS = SUNSPGMR(y, PREC_LEFT, 0);
            if(check_flag((void *)LS, "SUNDenseLinearSolver", 0)) return(1);
        }

        /* Set matrix and linear solver to Cvode */
        flag = CVodeSetLinearSolver(cvode_mem, LS, NULL);
        if(check_flag(&flag, "CVodeSetLinearSolver", 1)) return(1);

	/* Set the JAcobian-times-vector function */
	flag = CVodeSetJacTimes(cvode_mem, NULL, NULL);
	if(check_flag(&flag, "CVodeSetJacTimes", 1)) return(1);

	if (user_data->ianalytical_jacobian == 1) {
	    /* Set the preconditioner solve and setup functions */
	    flag = CVodeSetPreconditioner(cvode_mem, Precond, PSolve);
	    if(check_flag(&flag, "CVodeSetPreconditioner", 1)) return(1);
	}

    } else if (user_data->isolve_type == sparse_cusolver_solve) {

        LS = SUNLinSol_cuSolverSp_batchQR(y, A, user_data->cusolverHandle);
	if(check_flag((void *)LS, "SUNLinSol_cuSolverSp_batchQR", 0)) return(1);

        /* Set matrix and linear solver to Cvode */
        flag = CVodeSetLinearSolver(cvode_mem, LS, A);
        if(check_flag(&flag, "CVodeSetLinearSolver", 1)) return(1);

	/* Set the user-supplied Jacobian routine Jac */
        flag = CVodeSetJacFn(cvode_mem, cJac);
	if(check_flag(&flag, "CVodeSetJacFn", 1)) return(1); 

    } else {

        /* Create dense SUNLinearSolver object for use by CVode */ 
	LS = SUNLinSol_dense_custom(y, A, NCELLS, (NEQ+1), user_data->NNZ, stream);
	if(check_flag((void *)LS, "SUNDenseLinearSolver", 0)) return(1);

	/* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */
	flag = CVodeSetLinearSolver(cvode_mem, LS, A); 
	if(check_flag(&flag, "CVodeSetLinearSolver", 1)) return(1); 

	/* Set the user-supplied Jacobian routine Jac */
        flag = CVodeSetJacFn(cvode_mem, cJac);
	if(check_flag(&flag, "CVodeSetJacFn", 1)) return(1);
    }

    /* Set the max number of time steps */ 
    flag = CVodeSetMaxNumSteps(cvode_mem, 100000);
    if(check_flag(&flag, "CVodeSetMaxNumSteps", 1)) return(1);

    /* Set the max order */
    flag = CVodeSetMaxOrd(cvode_mem, 2);
    if(check_flag(&flag, "CVodeSetMaxOrd", 1)) return(1);

    BL_PROFILE_VAR("AroundCVODE", AroundCVODE);
    flag = CVode(cvode_mem, time_out, y, &time_init, CV_NORMAL);
    if (check_flag(&flag, "CVode", 1)) return(1);
    BL_PROFILE_VAR_STOP(AroundCVODE);

#ifdef MOD_REACTOR
    /* ONLY FOR PP */
    /*If reactor mode is activated, update time */
    *dt_react = time_init - *time;
    *time = time_init;
#endif

    /* Pack data to return in main routine external */
    BL_PROFILE_VAR_START(AsyncCpy);
    cudaMemcpyAsync(rY_in, yvec_d, ((NEQ+1)*NCELLS)*sizeof(realtype), cudaMemcpyDeviceToHost,stream);

    for  (int i = 0; i < NCELLS; i++) {
        rX_in[i] = rX_in[i] + (*dt_react) * rX_src_in[i];
    }
    BL_PROFILE_VAR_STOP(AsyncCpy);

    long int nfe;
    flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);

    if (user_data->iverbose > 1) {
        PrintFinalStats(cvode_mem);
    }
	
    N_VDestroy(y);          /* Free the y vector */
    CVodeFree(&cvode_mem);

    SUNLinSolFree(LS);

    if (user_data->isolve_type != iterative_gmres_solve)
    {
        SUNMatDestroy(A);
    }

    cudaFree(user_data->rhoe_init);
    cudaFree(user_data->rhoesrc_ext);
    cudaFree(user_data->rYsrc);

    if (user_data->ianalytical_jacobian == 1) {
        cudaFree(user_data->csr_row_count_d);
        cudaFree(user_data->csr_col_index_d);
        cudaFree(user_data->csr_jac_d);
        if (user_data->isolve_type == iterative_gmres_solve) {
            cudaFree(user_data->csr_val_d);

            cusolverStatus_t cusolver_status = cusolverSpDestroy(user_data->cusolverHandle);
            assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);

            cusolver_status = cusolverSpDestroyCsrqrInfo(user_data->info);
            assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);
 
	    cudaFree(user_data->buffer_qr);
	} else if (user_data->isolve_type == sparse_cusolver_solve) {
            cusolverStatus_t cusolver_status = cusolverSpDestroy(user_data->cusolverHandle);
            assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);

	    cusparseStatus_t cusparse_status = cusparseDestroy(user_data->cuSPHandle);
	    assert(cusparse_status == CUSPARSE_STATUS_SUCCESS);
	}
    }

    cudaFree(user_data);

    N_VDestroy(atol);          /* Free the atol vector */

    return nfe;
}
/**********************************/


/*
 * CPU routines
 */
/* RHS routine used in CVODE */
static int cF_RHS(realtype t, N_Vector y_in, N_Vector ydot_in, 
                  void *user_data){

    BL_PROFILE_VAR("fKernelSpec()", fKernelSpec);

    cudaError_t cuda_status = cudaSuccess;

    /* Get Device pointers for Kernel call */
    realtype *yvec_d      = N_VGetDeviceArrayPointer_Cuda(y_in);
    realtype *ydot_d      = N_VGetDeviceArrayPointer_Cuda(ydot_in);
    
    // allocate working space 
    UserData udata = static_cast<CVodeUserData*>(user_data);
    udata->dt_save = t;

    const auto ec = Gpu::ExecutionConfig(udata->ncells_d[0]);   
    //amrex::launch_global<<<ec.numBlocks, ec.numThreads, ec.sharedMem, udata->stream>>>(
    amrex::launch_global<<<udata->nbBlocks, udata->nbThreads, ec.sharedMem, udata->stream>>>(
        [=] AMREX_GPU_DEVICE () noexcept {
            for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
                     icell < udata->ncells_d[0]; icell += stride) {
                         fKernelSpec(icell, user_data, yvec_d, ydot_d, udata->rhoe_init, 
                                     udata->rhoesrc_ext, udata->rYsrc);    
        }
    }); 

    cuda_status = cudaStreamSynchronize(udata->stream);  
    assert(cuda_status == cudaSuccess);

    BL_PROFILE_VAR_STOP(fKernelSpec);
    
    return(0);
}


static int Precond(realtype tn, N_Vector u, N_Vector fu, booleantype jok,
               booleantype *jcurPtr, realtype gamma, void *user_data) {

    BL_PROFILE_VAR("Precond()", Precond);

    cudaError_t cuda_status = cudaSuccess;
    size_t workspaceInBytes, internalDataInBytes;
    cusolverStatus_t cusolver_status = CUSOLVER_STATUS_SUCCESS;

    workspaceInBytes = 0;
    internalDataInBytes = 0;

    /* Get Device pointers for Kernel call */
    realtype *u_d      = N_VGetDeviceArrayPointer_Cuda(u);
    realtype *udot_d   = N_VGetDeviceArrayPointer_Cuda(fu);

    // allocate working space 
    UserData udata = static_cast<CVodeUserData*>(user_data);
    udata->gamma_d = gamma;

    BL_PROFILE_VAR("fKernelComputeAJ()", fKernelComputeAJ);
    if (jok) {
        /* GPU tests */
        const auto ec = Gpu::ExecutionConfig(udata->ncells_d[0]);   
	amrex::launch_global<<<udata->nbBlocks, udata->nbThreads, ec.sharedMem, udata->stream>>>(
	    [=] AMREX_GPU_DEVICE () noexcept {
	            for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
	    	    icell < udata->ncells_d[0]; icell += stride) {
	    	    fKernelComputeAJsys(icell, user_data, u_d, udata->csr_val_d);
	    	}
        }); 
        *jcurPtr = SUNFALSE;
    } else {
        const auto ec = Gpu::ExecutionConfig(udata->ncells_d[0]);   
        amrex::launch_global<<<udata->nbBlocks, udata->nbThreads, ec.sharedMem, udata->stream>>>(
	    [=] AMREX_GPU_DEVICE () noexcept {
	            for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
	    	    icell < udata->ncells_d[0]; icell += stride) {
	    	    fKernelComputeallAJ(icell, user_data, u_d, udata->csr_val_d);
	    	}
        }); 
        *jcurPtr = SUNTRUE;
    }

    cuda_status = cudaStreamSynchronize(udata->stream);  
    assert(cuda_status == cudaSuccess);
    BL_PROFILE_VAR_STOP(fKernelComputeAJ);

    BL_PROFILE_VAR("InfoBatched(inPrecond)", InfoBatched);
    cusolver_status = cusolverSpDcsrqrBufferInfoBatched(udata->cusolverHandle,udata->neqs_per_cell[0]+1,udata->neqs_per_cell[0]+1, 
                                (udata->NNZ),
                                udata->descrA,
                                udata->csr_val_d,
                                udata->csr_row_count_d,
                                udata->csr_col_index_d,
                                udata->ncells_d[0],
                                udata->info,
                                &internalDataInBytes,
                                &workspaceInBytes);

    assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);

    cuda_status = cudaDeviceSynchronize();  
    assert(cuda_status == cudaSuccess);

    BL_PROFILE_VAR_STOP(InfoBatched);

    BL_PROFILE_VAR_STOP(Precond);

    return(0);
}



static int PSolve(realtype tn, N_Vector u, N_Vector fu, N_Vector r, N_Vector z,
                  realtype gamma, realtype delta, int lr, void *user_data)
{
    BL_PROFILE_VAR("Psolve()", cusolverPsolve);

    cudaError_t cuda_status = cudaSuccess;
    cusolverStatus_t cusolver_status = CUSOLVER_STATUS_SUCCESS;

    UserData udata = static_cast<CVodeUserData*>(user_data);

    realtype *z_d      = N_VGetDeviceArrayPointer_Cuda(z);
    realtype *r_d      = N_VGetDeviceArrayPointer_Cuda(r);

    cusolver_status = cusolverSpDcsrqrsvBatched(udata->cusolverHandle,udata->neqs_per_cell[0]+1,udata->neqs_per_cell[0]+1,
                               (udata->NNZ),
                               udata->descrA,
                               udata->csr_val_d,
                               udata->csr_row_count_d,
                               udata->csr_col_index_d,
                               r_d, 
                               z_d,
                               udata->ncells_d[0],
                               udata->info,
                               udata->buffer_qr);

        assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);

        cuda_status = cudaDeviceSynchronize();  
        assert(cuda_status == cudaSuccess);

        N_VCopyFromDevice_Cuda(z);
        N_VCopyFromDevice_Cuda(r);

	BL_PROFILE_VAR_STOP(cusolverPsolve);

        /* Checks */
        //if (udata->iverbose > 4) {
        //    for(int batchId = 0 ; batchId < udata->ncells_d[0]; batchId++){
        //        // measure |bj - Aj*xj|
        //        realtype *csrValAj = (udata->csr_val_d) + batchId * (udata->NNZ);
        //        double *xj       = N_VGetHostArrayPointer_Cuda(z) + batchId * (udata->neqs_per_cell[0]+1);
        //        double *bj       = N_VGetHostArrayPointer_Cuda(r) + batchId * (udata->neqs_per_cell[0]+1);
        //        // sup| bj - Aj*xj|
        //        double sup_res = 0;
        //        for(int row = 0 ; row < (udata->neqs_per_cell[0]+1) ; row++){
        //            printf("\n     row %d: ", row);
        //            const int start = udata->csr_row_count_d[row] - 1;
        //            const int end = udata->csr_row_count_d[row +1] - 1;
        //            double Ax = 0.0; // Aj(row,:)*xj
        //            for(int colidx = start ; colidx < end ; colidx++){
        //                const int col = udata->csr_col_index_d[colidx] - 1;
        //                const double Areg = csrValAj[colidx];
        //                const double xreg = xj[col];
        //                printf("  (%d, %14.8e, %14.8e, %14.8e) ", col,Areg,xreg,bj[row] );
        //                Ax = Ax + Areg * xreg;
        //            }
        //            double rresidi = bj[row] - Ax;
        //            sup_res = (sup_res > fabs(rresidi))? sup_res : fabs(rresidi);
        //        }
        //        printf("batchId %d: sup|bj - Aj*xj| = %E \n", batchId, sup_res);
        //    }
        //}

        return(0);
}

/* Will not work for cuSolver_sparse_solve right now */
static int cJac(realtype t, N_Vector y_in, N_Vector fy, SUNMatrix J,
                void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
        cudaError_t cuda_status = cudaSuccess;

        /* allocate working space */
        UserData udata = static_cast<CVodeUserData*>(user_data);

	/* Get Device pointers for Kernel call */
	realtype *yvec_d       = N_VGetDeviceArrayPointer_Cuda(y_in);
        realtype *ydot_d       = N_VGetDeviceArrayPointer_Cuda(fy);

        /* Jdata */
        realtype *Jdata;

        /* Fixed Indices and Pointers for Jacobian Matrix */
	BL_PROFILE_VAR("cJac::SparsityStuff",cJacSparsityStuff);
        int consP;
        if (udata->ireactor_type == 1){
            consP = 0 ;
        } else {
            consP = 1;
        }
	
        if (udata->isolve_type == sparse_cusolver_solve) {
	    Jdata   = SUNMatrix_cuSparse_Data(J);
            if ((SUNMatrix_cuSparse_Rows(J) != (udata->neqs_per_cell[0]+1)*(udata->ncells_d[0])) || 
	    	(SUNMatrix_cuSparse_Columns(J) != (udata->neqs_per_cell[0]+1)*(udata->ncells_d[0])) ||
                (SUNMatrix_cuSparse_NNZ(J) != udata->ncells_d[0] * udata->NNZ )) {
                        printf("Jac error: matrix is wrong size!\n");
                        return 1;
            }
	} else {
	    SPARSITY_PREPROC_SYST_CSR( (int*) SUNSparseMatrix_IndexValues(J), (int*) SUNSparseMatrix_IndexPointers(J), 
                                     &consP, udata->ncells_d[0], 0); 
	    /* Create empty chem Jacobian matrix (if not done already) */
	    //if (udata->R == NULL) {
	    //	udata->R = SUNSparseMatrix(SUNSparseMatrix_Rows(J),
	    //			SUNSparseMatrix_Columns(J),
	    //			SUNSparseMatrix_NNZ(J), CSR_MAT);
	    //}
	    /* check that vector/matrix dimensions match up */
            if ((SUNSparseMatrix_Rows(J) != (udata->neqs_per_cell[0]+1)*(udata->ncells_d[0])) || 
	    	(SUNSparseMatrix_Columns(J) != (udata->neqs_per_cell[0]+1)*(udata->ncells_d[0])) ||
                (SUNSparseMatrix_NNZ(J) != udata->ncells_d[0] * udata->NNZ )) {
                        printf("Jac error: matrix is wrong size!\n");
                        return 1;
            }
	}
	BL_PROFILE_VAR_STOP(cJacSparsityStuff);

	BL_PROFILE_VAR("Jacobian()", fKernelJac );
        if (udata->isolve_type == sparse_cusolver_solve) {
	    /* GPU tests */
            const auto ec = Gpu::ExecutionConfig(udata->ncells_d[0]);   
	    amrex::launch_global<<<udata->nbBlocks, udata->nbThreads, ec.sharedMem, udata->stream>>>(
	    [=] AMREX_GPU_DEVICE () noexcept {
	        for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
	    	icell < udata->ncells_d[0]; icell += stride) {
	            fKernelComputeAJchemCuSolver(icell, user_data, yvec_d, Jdata);
	    	}
            }); 
	} else {
	    /* GPU tests */
            const auto ec = Gpu::ExecutionConfig(udata->ncells_d[0]);   
	    //amrex::launch_global<<<ec.numBlocks, ec.numThreads, ec.sharedMem, udata->stream>>>(
	    amrex::launch_global<<<udata->nbBlocks, udata->nbThreads, ec.sharedMem, udata->stream>>>(
	    [=] AMREX_GPU_DEVICE () noexcept {
	        for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
	    	icell < udata->ncells_d[0]; icell += stride) {
	            fKernelComputeAJchem(icell, user_data, yvec_d);
	    	}
            }); 
	}
        cuda_status = cudaStreamSynchronize(udata->stream);  
        assert(cuda_status == cudaSuccess);
	BL_PROFILE_VAR_STOP(fKernelJac);

        if (udata->isolve_type == sparse_solve) {
            BL_PROFILE_VAR("cJac::memcpy",cJacmemcpy);
	    realtype *Jdata = SUNSparseMatrix_Data(J); 
            std::memcpy(Jdata, udata->csr_jac_d, sizeof(realtype)*udata->NNZ*(udata->ncells_d[0]));
	    BL_PROFILE_VAR_STOP(cJacmemcpy);
	}

        return(0);

}
/**********************************/


/*
 * CUDA kernels
 */
AMREX_GPU_DEVICE
inline
void 
fKernelSpec(int icell, void *user_data, 
            realtype *yvec_d, realtype *ydot_d,  
            double *rhoe_init, double *rhoesrc_ext, double *rYs)
{
  UserData udata = static_cast<CVodeUserData*>(user_data);

  amrex::Real mw[NUM_SPECIES];
  amrex::GpuArray<amrex::Real,NUM_SPECIES> massfrac;
  amrex::GpuArray<amrex::Real,NUM_SPECIES> ei_pt;
  amrex::GpuArray<amrex::Real,NUM_SPECIES> cdots_pt;
  amrex::Real Cv_pt, rho_pt, temp_pt, nrg_pt;

  int offset = icell * (NUM_SPECIES + 1); 

  /* MW CGS */
  get_mw(mw);

  /* rho */ 
  rho_pt = 0.0;
  for (int n = 0; n < NUM_SPECIES; n++) {
      rho_pt = rho_pt + yvec_d[offset + n];
  }

  /* Yks, C CGS*/
  for (int i = 0; i < NUM_SPECIES; i++){
      massfrac[i] = yvec_d[offset + i] / rho_pt;
  }

  /* NRG CGS */
  nrg_pt = (rhoe_init[icell] + rhoesrc_ext[icell]*(udata->dt_save)) /rho_pt;

  /* temp */
  temp_pt = yvec_d[offset + NUM_SPECIES];

  /* Additional var needed */
  if (udata->ireactor_type == eint_rho){
      /* UV REACTOR */
      EOS::EY2T(nrg_pt, massfrac.arr, temp_pt);
      EOS::T2Ei(temp_pt, ei_pt.arr);
      EOS::TY2Cv(temp_pt, massfrac.arr, Cv_pt);
  } else {
      /* HP REACTOR */
      EOS::HY2T(nrg_pt, massfrac.arr, temp_pt);
      EOS::TY2Cp(temp_pt, massfrac.arr, Cv_pt);
      EOS::T2Hi(temp_pt, ei_pt.arr);
  }

  EOS::RTY2WDOT(rho_pt, temp_pt, massfrac.arr, cdots_pt.arr);

  /* Fill ydot vect */
  ydot_d[offset + NUM_SPECIES] = rhoesrc_ext[icell];
  for (int i = 0; i < NUM_SPECIES; i++){
      ydot_d[offset + i]           = cdots_pt[i] + rYs[icell * NUM_SPECIES + i];
      ydot_d[offset + NUM_SPECIES] = ydot_d[offset + NUM_SPECIES]  - ydot_d[offset + i] * ei_pt[i];
  }
  ydot_d[offset + NUM_SPECIES] = ydot_d[offset + NUM_SPECIES] /(rho_pt * Cv_pt);
}


AMREX_GPU_DEVICE
inline
void 
fKernelComputeAJchemCuSolver(int ncell, void *user_data, realtype *u_d, realtype *Jdata)
{
  UserData udata = static_cast<CVodeUserData*>(user_data);

  amrex::Real mw[NUM_SPECIES];
  amrex::GpuArray<amrex::Real,NUM_SPECIES> massfrac;
  amrex::GpuArray<amrex::Real,(NUM_SPECIES+1)*(NUM_SPECIES+1)> Jmat_pt;
  amrex::Real rho_pt, temp_pt;

  int u_offset      = ncell * (NUM_SPECIES + 1); 
  int jac_offset = ncell * (udata->NNZ); 

  realtype* u_curr  = u_d + u_offset;
  realtype* csr_jac_cell         = Jdata + jac_offset;
  
  /* MW CGS */
  get_mw(mw);

  /* rho */ 
  rho_pt = 0.0;
  for (int n = 0; n < NUM_SPECIES; n++) {
      rho_pt = rho_pt + u_curr[n];
  }

  /* Yks, C CGS*/
  for (int i = 0; i < NUM_SPECIES; i++){
      massfrac[i] = u_curr[i] / rho_pt;
  }

  /* temp */
  temp_pt = u_curr[NUM_SPECIES];

  /* Additional var needed */
  int consP;
  if (udata->ireactor_type == 1){
      consP = 0 ;
  } else {
      consP = 1;
  }
  EOS::RTY2JAC(rho_pt, temp_pt, massfrac.arr, Jmat_pt.arr, consP); 

  /* renorm the DenseMat */
  for (int i = 0; i < udata->neqs_per_cell[0]; i++){
      for (int k = 0; k < udata->neqs_per_cell[0]; k++){
          Jmat_pt[k*(udata->neqs_per_cell[0]+1)+i] = Jmat_pt[k*(udata->neqs_per_cell[0]+1)+i] * mw[i] / mw[k];
      }
      Jmat_pt[i*(udata->neqs_per_cell[0]+1)+udata->neqs_per_cell[0]] = Jmat_pt[i*(udata->neqs_per_cell[0]+1)+udata->neqs_per_cell[0]] / mw[i]; 
      Jmat_pt[udata->neqs_per_cell[0]*(udata->neqs_per_cell[0]+1)+i] = Jmat_pt[udata->neqs_per_cell[0]*(udata->neqs_per_cell[0]+1)+i] * mw[i]; 
  }
  /* Fill the Sps Mat */
  int nbVals;
  for (int i = 1; i < udata->neqs_per_cell[0]+2; i++) {
      nbVals = udata->csr_row_count_d[i]-udata->csr_row_count_d[i-1];
      for (int j = 0; j < nbVals; j++) {
	      int idx_cell = udata->csr_col_index_d[ udata->csr_row_count_d[i-1] + j ];
              csr_jac_cell[ udata->csr_row_count_d[i-1] + j ] = Jmat_pt[ idx_cell * (udata->neqs_per_cell[0]+1) + i - 1 ]; 
      }
  }

}



AMREX_GPU_DEVICE
inline
void 
fKernelComputeAJchem(int ncell, void *user_data, realtype *u_d)
{
  UserData udata = static_cast<CVodeUserData*>(user_data);

  amrex::Real mw[NUM_SPECIES];
  amrex::GpuArray<amrex::Real,NUM_SPECIES> massfrac;
  amrex::GpuArray<amrex::Real,(NUM_SPECIES+1)*(NUM_SPECIES+1)> Jmat_pt;
  amrex::Real rho_pt, temp_pt;

  int u_offset   = ncell * (NUM_SPECIES + 1); 
  int jac_offset = ncell * (udata->NNZ); 

  realtype* u_curr               = u_d + u_offset;
  realtype* csr_jac_cell         = udata->csr_jac_d + jac_offset;
  
  /* MW CGS */
  get_mw(mw);

  /* rho */ 
  rho_pt = 0.0;
  for (int n = 0; n < NUM_SPECIES; n++) {
      rho_pt = rho_pt + u_curr[n];
  }

  /* Yks, C CGS*/
  for (int i = 0; i < NUM_SPECIES; i++){
      massfrac[i] = u_curr[i] / rho_pt;
  }

  /* temp */
  temp_pt = u_curr[NUM_SPECIES];

  /* Additional var needed */
  int consP;
  if (udata->ireactor_type == 1){
      consP = 0 ;
  } else {
      consP = 1;
  }
  EOS::RTY2JAC(rho_pt, temp_pt, massfrac.arr, Jmat_pt.arr, consP); 

  /* renorm the DenseMat */
  for (int i = 0; i < udata->neqs_per_cell[0]; i++){
      for (int k = 0; k < udata->neqs_per_cell[0]; k++){
          Jmat_pt[k*(udata->neqs_per_cell[0]+1)+i] = Jmat_pt[k*(udata->neqs_per_cell[0]+1)+i] * mw[i] / mw[k];
      }
      Jmat_pt[i*(udata->neqs_per_cell[0]+1)+udata->neqs_per_cell[0]] = Jmat_pt[i*(udata->neqs_per_cell[0]+1)+udata->neqs_per_cell[0]] / mw[i]; 
      Jmat_pt[udata->neqs_per_cell[0]*(udata->neqs_per_cell[0]+1)+i] = Jmat_pt[udata->neqs_per_cell[0]*(udata->neqs_per_cell[0]+1)+i] * mw[i]; 
  }
  /* Fill the Sps Mat */
  int nbVals;
  for (int i = 1; i < udata->neqs_per_cell[0]+2; i++) {
      nbVals = udata->csr_row_count_d[i]-udata->csr_row_count_d[i-1];
      for (int j = 0; j < nbVals; j++) {
          int idx_cell = udata->csr_col_index_d[ udata->csr_row_count_d[i-1] + j - 1] - 1 ;
              csr_jac_cell[ udata->csr_row_count_d[i-1] + j - 1 ] = Jmat_pt[ idx_cell * (udata->neqs_per_cell[0]+1) + i - 1 ]; 
      }
  }

}


AMREX_GPU_DEVICE
inline
void 
fKernelComputeallAJ(int ncell, void *user_data, realtype *u_d, double * csr_val_arg)
{
  UserData udata = static_cast<CVodeUserData*>(user_data);

  amrex::Real mw[NUM_SPECIES];
  amrex::GpuArray<amrex::Real,NUM_SPECIES> massfrac, activity;
  amrex::GpuArray<amrex::Real,(NUM_SPECIES+1)*(NUM_SPECIES+1)> Jmat_pt;
  amrex::Real rho_pt, temp_pt;

  int u_offset   = ncell * (NUM_SPECIES + 1); 
  int jac_offset = ncell * (udata->NNZ); 

  realtype* u_curr = u_d + u_offset;
  realtype* csr_jac_cell = udata->csr_jac_d + jac_offset;
  realtype* csr_val_cell = csr_val_arg + jac_offset;
  
  /* MW CGS */
  get_mw(mw);

  /* rho */ 
  rho_pt = 0.0;
  for (int n = 0; n < NUM_SPECIES; n++) {
      rho_pt = rho_pt + u_curr[n];
  }

  /* Yks, C CGS*/
  for (int i = 0; i < NUM_SPECIES; i++){
      massfrac[i] = u_curr[i] / rho_pt;
  }

  /* temp */
  temp_pt = u_curr[NUM_SPECIES];

  /* Activities */
  EOS::RTY2C(rho_pt, temp_pt, massfrac.arr, activity.arr);

  /* Additional var needed */
  int consP;
  if (udata->ireactor_type == eint_rho){
      consP = 0 ;
  } else {
      consP = 1;
  }
  DWDOT_SIMPLIFIED(Jmat_pt.arr, activity.arr, &temp_pt, &consP);

  /* renorm the DenseMat */
  for (int i = 0; i < udata->neqs_per_cell[0]; i++){
      for (int k = 0; k < udata->neqs_per_cell[0]; k++){
          Jmat_pt[k*(udata->neqs_per_cell[0]+1)+i] = Jmat_pt[k*(udata->neqs_per_cell[0]+1)+i] * mw[i] / mw[k];
      }
      Jmat_pt[i*(udata->neqs_per_cell[0]+1)+udata->neqs_per_cell[0]] = Jmat_pt[i*(udata->neqs_per_cell[0]+1)+udata->neqs_per_cell[0]] / mw[i]; 
      Jmat_pt[udata->neqs_per_cell[0]*(udata->neqs_per_cell[0]+1)+i] = Jmat_pt[udata->neqs_per_cell[0]*(udata->neqs_per_cell[0]+1)+i] * mw[i]; 
  }
  /* Fill the Sps Mat */
  int nbVals;
  for (int i = 1; i < udata->neqs_per_cell[0]+2; i++) {
      nbVals = udata->csr_row_count_d[i]-udata->csr_row_count_d[i-1];
      for (int j = 0; j < nbVals; j++) {
          int idx = udata->csr_col_index_d[ udata->csr_row_count_d[i-1] + j - 1 ] - 1;
          /* Scale by -gamma */
          /* Add identity matrix */
          if (idx == (i-1)) {
              csr_val_cell[ udata->csr_row_count_d[i-1] + j - 1 ] = 1.0 - (udata->gamma_d) * Jmat_pt[ idx * (udata->neqs_per_cell[0]+1) + idx ]; 
              csr_jac_cell[ udata->csr_row_count_d[i-1] + j - 1 ] = Jmat_pt[ idx * (udata->neqs_per_cell[0]+1) + idx ]; 
          } else {
              csr_val_cell[ udata->csr_row_count_d[i-1] + j - 1 ] = - (udata->gamma_d) * Jmat_pt[ idx * (udata->neqs_per_cell[0]+1) + i-1 ]; 
              csr_jac_cell[ udata->csr_row_count_d[i-1] + j - 1 ] = Jmat_pt[ idx * (udata->neqs_per_cell[0]+1) + i-1 ]; 
          }
      }
  }

}

AMREX_GPU_DEVICE
inline
void 
fKernelComputeAJsys(int ncell, void *user_data, realtype *u_d, double * csr_val_arg)
{
  UserData udata = static_cast<CVodeUserData*>(user_data);

  int jac_offset = ncell * (udata->NNZ); 

  realtype* csr_jac_cell = udata->csr_jac_d + jac_offset;
  realtype* csr_val_cell = csr_val_arg + jac_offset;
  
  int nbVals;
  for (int i = 1; i < udata->neqs_per_cell[0]+2; i++) {
      nbVals = udata->csr_row_count_d[i]-udata->csr_row_count_d[i-1];
      for (int j = 0; j < nbVals; j++) {
          int idx = udata->csr_col_index_d[ udata->csr_row_count_d[i-1] + j - 1 ] - 1;
          /* Scale by -gamma */
          /* Add identity matrix */
          if (idx == (i-1)) {
              csr_val_cell[ udata->csr_row_count_d[i-1] + j - 1 ] = 1.0 - (udata->gamma_d) * csr_jac_cell[ udata->csr_row_count_d[i-1] + j - 1 ];
          } else {
              csr_val_cell[ udata->csr_row_count_d[i-1] + j - 1 ] = - (udata->gamma_d) * csr_jac_cell[ udata->csr_row_count_d[i-1] + j - 1 ];
          }
      }
  }

}


__global__
void 
fKernelDenseSolve(int ncell, realtype *x_d, realtype *b_d,
          int subsys_size, int subsys_nnz, realtype *csr_val)
{

  int stride = blockDim.x*gridDim.x;
  
  for (int icell = blockDim.x*blockIdx.x+threadIdx.x;
           icell < ncell; icell += stride) {   
               int offset   = icell * subsys_size; 
               int offset_A = icell * subsys_nnz;

               realtype* csr_val_cell = csr_val + offset_A;
               realtype* x_cell       = x_d + offset;
               realtype* b_cell       = b_d + offset;

               /* Solve the subsystem of the cell */
               sgjsolve(csr_val_cell, x_cell, b_cell);
  }  
}

/**********************************/


/* 
 * CUSTOM SOLVER STUFF
 */
SUNLinearSolver SUNLinSol_dense_custom(N_Vector y, SUNMatrix A, 
        int nsubsys, int subsys_size, int subsys_nnz, cudaStream_t stream) 
{
  SUNLinearSolver S;
  SUNLinearSolverContent_Dense_custom content;

  /* Check that required arguments are not NULL */
  if (y == NULL || A == NULL) return(NULL);

  /* Check compatibility with supplied SUNMatrix and N_Vector */
  if (N_VGetVectorID(y) != SUNDIALS_NVEC_CUDA) return(NULL);
  if (SUNMatGetID(A) != SUNMATRIX_SPARSE) return(NULL); 

  /* Matrix should be square */
  if (SUNSparseMatrix_Columns(A) != SUNSparseMatrix_Rows(A)) return(NULL);

  /* Check that it is a CSR matrix */
  if (SUNSparseMatrix_SparseType(A) != CSR_MAT) return(NULL);

  /* Check that the vector is using managed memory */
  if (!N_VIsManagedMemory_Cuda(y)) return(NULL);

  /* Matrix and vector dimensions must agree */
  if (N_VGetLength(y) != SUNSparseMatrix_Columns(A)) return(NULL);

  /* All subsystems must be the same size */ 
  if (SUNSparseMatrix_Columns(A) != (subsys_size * nsubsys)) return(NULL);

  /* Number of nonzeros per subsys must be the same */
  if (SUNSparseMatrix_NNZ(A) != (subsys_nnz * nsubsys)) return(NULL);

  /* Allocate device memory for the matrix */
  int *d_colind, *d_rowptr;
  N_Vector d_values;

  d_colind = NULL;
  d_rowptr = NULL;
  d_values = NULL;

  //cudaError_t cuerr;
  //cuerr = cudaMalloc(&d_colind, sizeof(int) * (subsys_nnz * nsubsys));
  //if (cuerr != cudaSuccess) return(NULL); 

  //cuerr = cudaMalloc(&d_rowptr, sizeof(int) * (subsys_size * nsubsys + 1));
  //if (cuerr != cudaSuccess) { cudaFree(d_colind); return(NULL); }

  d_values = N_VNewManaged_Cuda(subsys_nnz * nsubsys);

  /* copy matrix stuff to the device */
  //cuerr = cudaMemcpy(d_colind, SUNSparseMatrix_IndexValues(A),
  //          sizeof(int) * subsys_nnz * nsubsys, cudaMemcpyHostToDevice);
  //if (cuerr != cudaSuccess) { 
  //        cudaFree(d_rowptr); 
  //        cudaFree(d_colind); 
  //        N_VDestroy(d_values);
  //        return(NULL); 
  //};

  //cuerr = cudaMemcpy(d_rowptr, SUNSparseMatrix_IndexPointers(A),
  //          sizeof(int) * (subsys_size * nsubsys + 1), cudaMemcpyHostToDevice);
  //if (cuerr != cudaSuccess) { 
  //        cudaFree(d_rowptr); 
  //        cudaFree(d_colind); 
  //        N_VDestroy(d_values);
  //        return(NULL); 
  //};

  /* Create an empty linear solver */
  S = NULL;
  S = SUNLinSolNewEmpty(); 
  if (S == NULL) { 
     //cudaFree(d_rowptr); 
     //cudaFree(d_colind); 
     N_VDestroy(d_values);
     return(NULL);
  }

  /* Attach operations */ 
  S->ops->gettype    = SUNLinSolGetType_Dense_custom;
  S->ops->setup      = SUNLinSolSetup_Dense_custom;  
  S->ops->solve      = SUNLinSolSolve_Dense_custom;
  S->ops->free       = SUNLinSolFree_Dense_custom;

  /* Create content */
  content = NULL; 
  content = (SUNLinearSolverContent_Dense_custom) malloc(sizeof *content);
  if (content == NULL) { 
      //cudaFree(d_rowptr); 
      //cudaFree(d_colind); 
      N_VDestroy(d_values);
      SUNLinSolFree(S); 
      return(NULL); 
  }

  /* Attach content */
  S->content = content; 

  /* Fill content */ 
  content->last_flag   = 0;
  content->nsubsys     = nsubsys;
  content->subsys_size = subsys_size;
  content->subsys_nnz  = subsys_nnz;
  content->d_values    = d_values;
  content->d_colind    = d_colind;
  content->d_rowptr    = d_rowptr; 
  content->nbBlocks    = std::max(1,nsubsys/32); 
  content->nbThreads   = 32;
  content->stream      = stream; 

  return(S);
}


SUNLinearSolver_Type SUNLinSolGetType_Dense_custom(SUNLinearSolver S) 
{
  return(SUNLINEARSOLVER_DIRECT);
}

int SUNLinSolSetup_Dense_custom(SUNLinearSolver S, SUNMatrix A) 
{
  cudaError_t cuerr = cudaSuccess;

  realtype *data_d   = N_VGetDeviceArrayPointer_Cuda(SUN_CUSP_DVALUES(S));

  BL_PROFILE_VAR("SETUP::MemCpyMatrix",MemCpyMatrix);
  cuerr = cudaMemcpy(data_d, SUNSparseMatrix_Data(A),
                     sizeof(double) * (SUN_CUSP_SUBSYS_NNZ(S) * SUN_CUSP_NUM_SUBSYS(S)), cudaMemcpyHostToDevice);
  if (cuerr != cudaSuccess) {
      SUN_CUSP_LASTFLAG(S) = SUNLS_MEM_FAIL;
      amrex::Abort("\nPB MEMCPY\n");
  }
  BL_PROFILE_VAR_STOP(MemCpyMatrix);

  return(SUNLS_SUCCESS);
}

int SUNLinSolSolve_Dense_custom(SUNLinearSolver S, SUNMatrix A, N_Vector x,
        N_Vector b, realtype tol)
{
  cudaError_t cuda_status = cudaSuccess;

  /* Get Device pointers for Kernel call */
  realtype *x_d      = N_VGetDeviceArrayPointer_Cuda(x);
  realtype *b_d      = N_VGetDeviceArrayPointer_Cuda(b);

  realtype *data_d   = N_VGetDeviceArrayPointer_Cuda(SUN_CUSP_DVALUES(S));

  BL_PROFILE_VAR("fKernelDenseSolve()", fKernelDenseSolve);
  const auto ec = Gpu::ExecutionConfig(SUN_CUSP_NUM_SUBSYS(S));  
  // TODO: why is this AMREX version NOT working ?
  //amrex::launch_global<<<ec.numBlocks, ec.numThreads, ec.sharedMem, SUN_CUSP_STREAM(S)>>>(
  //    [=] AMREX_GPU_DEVICE () noexcept {
  //        for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
  //            icell < SUN_CUSP_NUM_SUBSYS(S); icell += stride) {
  //            fKernelDenseSolve(icell, x_d, b_d,
  //                  SUN_CUSP_SUBSYS_SIZE(S), SUN_CUSP_SUBSYS_NNZ(S), data_d);
  //        }
  //    }); 
  //fKernelDenseSolve<<<ec.numBlocks, ec.numThreads, ec.sharedMem, SUN_CUSP_STREAM(S)>>>(SUN_CUSP_NUM_SUBSYS(S), x_d, b_d, 
  fKernelDenseSolve<<<SUN_CUSP_NBLOCK(S), SUN_CUSP_NTHREAD(S), ec.sharedMem, SUN_CUSP_STREAM(S)>>>
                   (SUN_CUSP_NUM_SUBSYS(S), x_d, b_d, SUN_CUSP_SUBSYS_SIZE(S), SUN_CUSP_SUBSYS_NNZ(S), data_d);

  cuda_status = cudaStreamSynchronize(SUN_CUSP_STREAM(S));  
  assert(cuda_status == cudaSuccess);

  BL_PROFILE_VAR_STOP(fKernelDenseSolve);

  return(SUNLS_SUCCESS);
}

int SUNLinSolFree_Dense_custom(SUNLinearSolver S) 
{
  /* return with success if already freed */
  if (S == NULL) return(SUNLS_SUCCESS);

  /* free stuff in the content structure */
  //cudaFree(SUN_CUSP_DCOLIND(S));
  //cudaFree(SUN_CUSP_DROWPTR(S));
  N_VDestroy(SUN_CUSP_DVALUES(S)); 

  /* free content structure */
  if (S->content) {
      free(S->content);
      S->content = NULL;
  }

  /* free ops structure */
  if (S->ops) {
      free(S->ops);
      S->ops = NULL;
  }

  /* free the actual SUNLinSol */
  free(S);
  S = NULL;

  return(SUNLS_SUCCESS);
}

/**********************************/


/* 
 * OTHERS
*/

/* Get and print some final statistics */
static void PrintFinalStats(void *cvodeMem)
{
  long lenrw, leniw ;
  long lenrwLS, leniwLS;
  long int nst, nfe, nsetups, nni, ncfn, netf;
  long int nli, npe, nps, ncfl, nfeLS;
  int flag;

  flag = CVodeGetWorkSpace(cvodeMem, &lenrw, &leniw);
  check_flag(&flag, "CVodeGetWorkSpace", 1);
  flag = CVodeGetNumSteps(cvodeMem, &nst);
  check_flag(&flag, "CVodeGetNumSteps", 1);
  flag = CVodeGetNumRhsEvals(cvodeMem, &nfe);
  check_flag(&flag, "CVodeGetNumRhsEvals", 1);
  flag = CVodeGetNumLinSolvSetups(cvodeMem, &nsetups);
  check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
  flag = CVodeGetNumErrTestFails(cvodeMem, &netf);
  check_flag(&flag, "CVodeGetNumErrTestFails", 1);
  flag = CVodeGetNumNonlinSolvIters(cvodeMem, &nni);
  check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvodeMem, &ncfn);
  check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);

  flag = CVodeGetLinWorkSpace(cvodeMem, &lenrwLS, &leniwLS);
  check_flag(&flag, "CVodeGetLinWorkSpace", 1);
  flag = CVodeGetNumLinIters(cvodeMem, &nli);
  check_flag(&flag, "CVodeGetNumLinIters", 1);
  //flag = CVodeGetNumJacEvals(cvodeMem, &nje);
  //check_flag(&flag, "CVodeGetNumJacEvals", 1);
  flag = CVodeGetNumLinRhsEvals(cvodeMem, &nfeLS);
  check_flag(&flag, "CVodeGetNumLinRhsEvals", 1);

  flag = CVodeGetNumPrecEvals(cvodeMem, &npe);
  check_flag(&flag, "CVodeGetNumPrecEvals", 1);
  flag = CVodeGetNumPrecSolves(cvodeMem, &nps);
  check_flag(&flag, "CVodeGetNumPrecSolves", 1);
  
  flag = CVodeGetNumLinConvFails(cvodeMem, &ncfl);
  check_flag(&flag, "CVodeGetNumLinConvFails", 1);

  amrex::Print() <<"\nFinal Statistics: \n";
  amrex::Print() <<"lenrw      = " << lenrw   <<"    leniw          = " << leniw   << "\n";
  amrex::Print() <<"lenrwLS    = " << lenrwLS <<"    leniwLS        = " << leniwLS << "\n";
  amrex::Print() <<"nSteps     = " << nst     <<"\n";
  amrex::Print() <<"nRHSeval   = " << nfe     <<"     nLinRHSeval   = " << nfeLS   << "\n";
  amrex::Print() <<"nnLinIt    = " << nni     <<"     nLinIt        = " << nli     << "\n";
  amrex::Print() <<"nLinsetups = " << nsetups <<"     nErrtf        = " << netf    << "\n";
  amrex::Print() <<"nPreceval  = " << npe     <<"     nPrecsolve    = " << nps     << "\n";
  amrex::Print() <<"nConvfail  = " << ncfn    <<"     nLinConvfail  = " << ncfl    << "\n\n";

}

/* Check function return value...
     opt == 0 means SUNDIALS function allocates memory so check if
              returned NULL pointer
     opt == 1 means SUNDIALS function returns a flag so check if
              flag >= 0
     opt == 2 means function allocates memory so check if returned
              NULL pointer */
static int check_flag(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
      funcname);
      return(1); 
  }
  /* Check if flag < 0 */
  else if (opt == 1) {
      errflag = (int *) flagvalue;
      if (*errflag < 0) {
          fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
          funcname, *errflag);
          return(1); 
      }
  }
  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
      fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
      funcname);
      return(1); 
  }

  return(0);
}

/* End of file  */
