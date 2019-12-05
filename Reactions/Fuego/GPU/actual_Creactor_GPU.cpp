#include <actual_Creactor_GPU.h> 
#include <AMReX_ParmParse.H>

#include <chemistry_file.H>
#include <Fuego_EOS.H>
#include "mechanism.h"

#include <AMReX_Gpu.H>

using namespace amrex;

/**********************************/
/* Infos to print once */
int reactor_info(const int* cvode_iE,const int* Ncells){ 

	int mm, ii;

        printf("Nb of spec is %d \n", NUM_SPECIES);

	/* Args */
	printf("Ncells in one solve is %d\n",*Ncells);
	printf("Reactor type is %d\n",*cvode_iE);

	/* ParmParse from the inputs file */ 
	amrex::ParmParse pp("ns");
	pp.query("cvode_iJac",mm);
	pp.query("cvode_iDense", ii);

        /* User data */
        if ((ii == 99) && (mm == 1)) { 
            printf("Using an Iterative GMRES Solver with sparse simplified preconditioning \n");
	    int nJdata;
            int HP;
            if (*cvode_iE == 1) {
                HP = 0;
            } else {
                HP = 1;
            }
            /* Precond data */ 
            SPARSITY_INFO_SYST_SIMPLIFIED(&nJdata,&HP);
            printf("--> SPARSE Preconditioner -- non zero entries %d represents %f %% fill pattern.\n", nJdata, nJdata/float((NUM_SPECIES+1) * (NUM_SPECIES+1)) *100.0);
        } else if ((ii == 99) && (mm == 0 )) {
            printf("\n--> Using an Iterative Solver without preconditionning \n");
	} else {
	    amrex::Abort("\n--> Only solver implemented is an iterative GMRES ...\n");
	}

        printf(" --> DONE WITH INITIALIZATION (GPU) %d \n", *cvode_iE);

	return(0);
}


/* Main routine for external looping */
int react(realtype *rY_in, realtype *rY_src_in, 
		realtype *rX_in, realtype *rX_src_in,
                realtype *dt_react, realtype *time,
                const int* cvode_iE,const int* Ncells, cudaStream_t stream) {

        /* CVODE */
        N_Vector y         = NULL;
        SUNLinearSolver LS = NULL;
        void *cvode_mem    = NULL;
        /* Misc */
	int flag;
        int NCELLS, NEQ, neq_tot;
	realtype reltol;
	N_Vector atol;
	realtype *ratol;
        /* tmp vects */
        int iE_Creact, iJac_Creact;

        /* cuSolver */
        size_t workspaceInBytes, internalDataInBytes;
        cusolverStatus_t cusolver_status = CUSOLVER_STATUS_SUCCESS;
        cusparseStatus_t cusparse_status = CUSPARSE_STATUS_SUCCESS;
        cudaError_t cudaStat1            = cudaSuccess;

        workspaceInBytes = 0;
        internalDataInBytes = 0;

	NEQ = NUM_SPECIES;

	/* ParmParse from the inputs file */ 
	amrex::ParmParse pp("ns");
	pp.query("cvode_iJac",iJac_Creact);

	/* Args */
	iE_Creact      = *cvode_iE;
	NCELLS         = *Ncells;
        neq_tot        = (NEQ + 1) * NCELLS;

        /* User data */
        UserData user_data;
	BL_PROFILE_VAR("AllocsInCVODE", AllocsCVODE);
        cudaMallocManaged(&user_data, sizeof(struct CVodeUserData));
	BL_PROFILE_VAR_STOP(AllocsCVODE);
        user_data->ncells_d[0]      = NCELLS;
        user_data->neqs_per_cell[0] = NEQ;
        user_data->flagP            = iE_Creact; 
        user_data->iJac             = iJac_Creact;
        user_data->iverbose         = 1;
        user_data->stream           = stream;

        if (iJac_Creact == 1) { 
            int HP;
            if (iE_Creact == 1) {
                HP = 0;
            } else {
                HP = 1;
            }
            // Find sparsity pattern to fill structure of sparse matrix
	    BL_PROFILE_VAR("SparsityFuegoStuff", SparsityStuff);
            SPARSITY_INFO_SYST_SIMPLIFIED(&(user_data->NNZ),&HP);
	    BL_PROFILE_VAR_STOP(SparsityStuff);

	    BL_PROFILE_VAR_START(AllocsCVODE);
            cudaMallocManaged(&(user_data->csr_row_count_d), (NEQ+2) * sizeof(int));
            cudaMallocManaged(&(user_data->csr_col_index_d), user_data->NNZ * sizeof(int));
            cudaMallocManaged(&(user_data->csr_jac_d), user_data->NNZ * NCELLS * sizeof(double));
            cudaMallocManaged(&(user_data->csr_val_d), user_data->NNZ * NCELLS * sizeof(double));
	    BL_PROFILE_VAR_STOP(AllocsCVODE);

	    BL_PROFILE_VAR_START(SparsityStuff);
            SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(user_data->csr_col_index_d, user_data->csr_row_count_d, &HP);
	    BL_PROFILE_VAR_STOP(SparsityStuff);

            // Create Sparse batch QR solver
            // qr info and matrix descriptor
	    BL_PROFILE_VAR("CuSolverInit", CuSolverInit);
            cusolver_status = cusolverSpCreate(&(user_data->cusolverHandle));
            assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);

            cusparse_status = cusparseCreateMatDescr(&(user_data->descrA)); 
            assert(cusparse_status == CUSPARSE_STATUS_SUCCESS);

            cusparse_status = cusparseSetMatType(user_data->descrA, CUSPARSE_MATRIX_TYPE_GENERAL);
            cusparseSetMatIndexBase(user_data->descrA, CUSPARSE_INDEX_BASE_ONE);
 
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
	    BL_PROFILE_VAR_STOP(CuSolverInit);
            
	    BL_PROFILE_VAR_START(AllocsCVODE);
            cudaStat1 = cudaMalloc((void**)&(user_data->buffer_qr), workspaceInBytes);
            assert(cudaStat1 == cudaSuccess);
	    BL_PROFILE_VAR_STOP(AllocsCVODE);

	    /*
            free_mem = 0;
            total_mem = 0;
            cudaStat1 = cudaMemGetInfo( &free_mem, &total_mem );
            assert( cudaSuccess == cudaStat1 );
            std::cout<<"(AFTER AWS) Free: "<< free_mem<< " Tot: "<<total_mem<<std::endl;
            std::cout<<" WORKSPACE "<< workspaceInBytes<<" INTERNAL DATA "<< internalDataInBytes<<"\n"<<std::endl;
	    */
        }

	/* Definition of main vector */
	y = N_VNew_Cuda(neq_tot);
	if(check_flag((void*)y, "N_VNew_Cuda", 0)) return(1);

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
	cudaMemcpyAsync(yvec_d, rY_in, sizeof(realtype) * ((NEQ+1)*NCELLS), cudaMemcpyHostToDevice,stream);
	// rhoY_src_ext
	cudaMemcpyAsync(user_data->rYsrc, rY_src_in, (NEQ*NCELLS)*sizeof(double), cudaMemcpyHostToDevice,stream);
	// rhoE/rhoH
	cudaMemcpyAsync(user_data->rhoe_init, rX_in, sizeof(realtype) * NCELLS, cudaMemcpyHostToDevice, stream);
	cudaMemcpyAsync(user_data->rhoesrc_ext, rX_src_in, sizeof(realtype) * NCELLS, cudaMemcpyHostToDevice,stream);
	BL_PROFILE_VAR_STOP(AsyncCpy)

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
        if (iJac_Creact == 0) { 
	        LS = SUNSPGMR(y, PREC_NONE, 0);
	        if(check_flag((void *)LS, "SUNDenseLinearSolver", 0)) return(1);
        } else { 
                LS = SUNSPGMR(y, PREC_LEFT, 0);
                if(check_flag((void *)LS, "SUNDenseLinearSolver", 0)) return(1);
        }

        /* Set CVSpils linear solver to LS */
        flag = CVSpilsSetLinearSolver(cvode_mem, LS);
        if(check_flag(&flag, "CVSpilsSetLinearSolver", 1)) return(1);

	if (user_data->iJac == 1) {
	    /* Set the JAcobian-times-vector function */
	    flag = CVSpilsSetJacTimes(cvode_mem, NULL, NULL);
	    if(check_flag(&flag, "CVSpilsSetJacTimes", 1)) return(1);

	    /* Set the preconditioner solve and setup functions */
	    flag = CVodeSetPreconditioner(cvode_mem, Precond, PSolve);
	    if(check_flag(&flag, "CVSpilsSetPreconditioner", 1)) return(1);
	}

        /* Set the max number of time steps */ 
	flag = CVodeSetMaxNumSteps(cvode_mem, 100000);
	if(check_flag(&flag, "CVodeSetMaxNumSteps", 1)) return(1);

        /* Set the max order */
        flag = CVodeSetMaxOrd(cvode_mem, 2);
        if(check_flag(&flag, "CVodeSetMaxOrd", 1)) return(1);

	/* Call CVODE: ReInit for convergence */
        //CVodeReInit(cvode_mem, time_init, y);

	BL_PROFILE_VAR("AroundCVODE", AroundCVODE);
	flag = CVode(cvode_mem, time_out, y, &time_init, CV_NORMAL);
	if (check_flag(&flag, "CVode", 1)) return(1);
	BL_PROFILE_VAR_STOP(AroundCVODE);

        /* ONLY FOR PP */
#ifdef MOD_REACTOR
	/* If reactor mode is activated, update time */
        *dt_react = time_init - *time;
        *time = time_init;
#endif

	/* Pack data to return in main routine external */
	BL_PROFILE_VAR_START(AsyncCpy)
	cudaMemcpyAsync(rY_in, yvec_d, ((NEQ+1)*NCELLS)*sizeof(realtype), cudaMemcpyDeviceToHost,stream);

	for  (int i = 0; i < NCELLS; i++) {
            rX_in[i] = rX_in[i] + (*dt_react) * rX_src_in[i];
	}
	BL_PROFILE_VAR_STOP(AsyncCpy)

        long int nfe;
	flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);

	SUNLinSolFree(LS);
	N_VDestroy(y);          /* Free the y vector */
	CVodeFree(&cvode_mem);

	cudaFree(user_data->rhoe_init);
        cudaFree(user_data->rhoesrc_ext);
	cudaFree(user_data->rYsrc);
	cudaFree(user_data->csr_row_count_d);
	cudaFree(user_data->csr_col_index_d);
	cudaFree(user_data->csr_jac_d);
	cudaFree(user_data->csr_val_d);

        cusolver_status = cusolverSpDestroy(user_data->cusolverHandle);
        assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);

        cusolver_status = cusolverSpDestroyCsrqrInfo(user_data->info);
        assert(cusolver_status == CUSOLVER_STATUS_SUCCESS);
 
	cudaFree(user_data->buffer_qr);

	cudaFree(user_data);

	N_VDestroy(atol);          /* Free the atol vector */

        //cuerr = cudaStreamDestroy(stream); /* Free and cleanup the CUDA stream */    
        //if(cuerr != cudaSuccess) { printf("Error: cudaStreamDestroy() failed\n"); return(1); }

	return nfe;
}



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

	// UV !!
	/* GPU tests */
        const auto ec = Gpu::ExecutionConfig(udata->ncells_d[0]);   
	amrex::launch_global<<<ec.numBlocks, ec.numThreads, ec.sharedMem, udata->stream>>>(
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

AMREX_GPU_DEVICE
inline
void 
fKernelSpec(int icell, void *user_data, 
            realtype *yvec_d, realtype *ydot_d,  
            double *rhoe_init, double *rhoesrc_ext, double *rYs)
{
  UserData udata = static_cast<CVodeUserData*>(user_data);

  EOS eos;

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
  if (udata->flagP == 1){
      /* UV REACTOR */
      eos.eos_EY2T(massfrac.arr, nrg_pt, temp_pt);
      eos.eos_T2EI(temp_pt, ei_pt.arr);
      eos.eos_TY2Cv(temp_pt, massfrac.arr, &Cv_pt);
  }else {
      /* HP REACTOR */
      eos.eos_HY2T(massfrac.arr, nrg_pt, temp_pt);
      eos.eos_TY2Cp(temp_pt, massfrac.arr, &Cv_pt);
      eos.eos_T2HI(temp_pt, ei_pt.arr);
  }

  eos.eos_RTY2W(rho_pt, temp_pt, massfrac.arr, cdots_pt.arr);

  /* Fill ydot vect */
  ydot_d[offset + NUM_SPECIES] = rhoesrc_ext[icell];
  for (int i = 0; i < NUM_SPECIES; i++){
      ydot_d[offset + i]           = cdots_pt[i] * mw[i] + rYs[icell * NUM_SPECIES + i];
      ydot_d[offset + NUM_SPECIES] = ydot_d[offset + NUM_SPECIES]  - ydot_d[offset + i] * ei_pt[i];
  }
  ydot_d[offset + NUM_SPECIES] = ydot_d[offset + NUM_SPECIES] /(rho_pt * Cv_pt);
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
	    amrex::launch_global<<<ec.numBlocks, ec.numThreads, ec.sharedMem, udata->stream>>>(
	    [=] AMREX_GPU_DEVICE () noexcept {
	            for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
	    	    icell < udata->ncells_d[0]; icell += stride) {
	    	    fKernelComputeAJ(icell, user_data, u_d, udot_d, udata->csr_val_d);
	    	}
            }); 
            cuda_status = cudaStreamSynchronize(udata->stream);  
            assert(cuda_status == cudaSuccess);
            *jcurPtr = SUNFALSE;
        } else {
            const auto ec = Gpu::ExecutionConfig(udata->ncells_d[0]);   
	    amrex::launch_global<<<ec.numBlocks, ec.numThreads, ec.sharedMem, udata->stream>>>(
	    [=] AMREX_GPU_DEVICE () noexcept {
	            for (int icell = blockDim.x*blockIdx.x+threadIdx.x, stride = blockDim.x*gridDim.x;
	    	    icell < udata->ncells_d[0]; icell += stride) {
	    	    fKernelComputeAJ(icell, user_data, u_d, udot_d, udata->csr_val_d);
	    	}
            }); 
            //cuda_status = cudaDeviceSynchronize();  
            cuda_status = cudaStreamSynchronize(udata->stream);  
            assert(cuda_status == cudaSuccess);
            *jcurPtr = SUNTRUE;
        }
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

        cuda_status = cudaDeviceSynchronize();  
        assert(cuda_status == cudaSuccess);

	BL_PROFILE_VAR_STOP(cusolverPsolve);


        /* Checks */
        N_VCopyFromDevice_Cuda(z);
        N_VCopyFromDevice_Cuda(r);

        //if (udata->iverbose > 4) {
        //    for(int batchId = 0 ; batchId < udata->ncells_d[0]; batchId++){
        //        // measure |bj - Aj*xj|
        //        double *csrValAj = (udata->csr_val_d) + batchId * (udata->NNZ);
        //        double *xj       = z_h + batchId * (udata->neqs_per_cell[0]+1);
        //        double *bj       = r_h + batchId * (udata->neqs_per_cell[0]+1);
        //        // sup| bj - Aj*xj|
        //        double sup_res = 0;
        //        for(int row = 0 ; row < (udata->neqs_per_cell[0]+1) ; row++){
        //            const int start = udata->csr_row_count_d[row] - 1;
        //            const int end = udata->csr_row_count_d[row +1] - 1;
        //            //printf("     row %d =  %d values \n", row, end - start);
        //            double Ax = 0.0; // Aj(row,:)*xj
        //            for(int colidx = start ; colidx < end ; colidx++){
        //                const int col = udata->csr_col_index_d[colidx] - 1;
        //                const double Areg = csrValAj[colidx];
        //                const double xreg = xj[col];
        //                //printf("        Value %d = col is %d and A is %E \n", colidx, col, Areg);
        //                Ax = Ax + Areg * xreg;
        //            }
        //            double r = bj[row] - Ax;
        //            sup_res = (sup_res > fabs(r))? sup_res : fabs(r);
        //        }
        //        printf("batchId %d: sup|bj - Aj*xj| = %E \n", batchId, sup_res);
        //    }
        //}

        return(0);
}


/*
 * CUDA kernels
 */
AMREX_GPU_DEVICE
inline
void 
fKernelComputeAJ(int ncell, void *user_data, realtype *u_d, realtype *udot_d, double * csr_val_arg)
{
  UserData udata = static_cast<CVodeUserData*>(user_data);

  EOS eos;

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
  eos.eos_RTY2C(rho_pt, temp_pt, massfrac.arr, activity.arr);

  /* Additional var needed */
  int consP;
  if (udata->flagP == 1){
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

//__global__ void fKernelJacCSR(realtype t, void *user_data,
//                                          realtype *yvec_d, realtype *ydot_d,
//                                          realtype* csr_jac,
//                                          const int size, const int nnz, 
//                                          const int nbatched)
//{
//
//    UserData udata = static_cast<CVodeUserData*>(user_data);
//    int tid = blockIdx.x * blockDim.x + threadIdx.x;
//
//    if (tid < nbatched) {
//        realtype activity[21];
//        realtype mw[21];
//        realtype temp;
//        realtype Jmat_pt[484];
//
//        int jac_offset = tid * nnz;
//        int y_offset = tid * size;
//
//        realtype* csr_jac_cell = csr_jac + jac_offset;
//        realtype* actual_y = yvec_d + y_offset;
//
//        /* MW CGS */
//        molecularWeight_d(mw);
//        /* rho */ 
//        //realtype rho = 0.0;
//        //for (int i = 0; i < udata->neqs_per_cell[0]; i++){
//        //    rho = rho + actual_y[i];
//        //}
//        /* temp */
//        temp = actual_y[udata->neqs_per_cell[0]];
//        /* Yks, C CGS*/
//        for (int i = 0; i < udata->neqs_per_cell[0]; i++){
//	    activity[i] = actual_y[i]/(mw[i]);
//        }
//        /* Fuego calls on device 
//         * NB to be more accurate should use energy to
//         * recompute temp ...      */
//        if (udata->flagP == 1){
//            int consP = 0 ;
//            dwdot_d(Jmat_pt, activity, &temp, &consP);
//        } else {
//            int consP = 1 ;
//            dwdot_d(Jmat_pt, activity, &temp, &consP);
//        }
//        /* fill the sunMat */
//        for (int k = 0; k < udata->neqs_per_cell[0]; k++){
//	    for (int i = 0; i < udata->neqs_per_cell[0]; i++){
//                csr_jac_cell[k*(udata->neqs_per_cell[0]+1)+i] = Jmat_pt[i*(udata->neqs_per_cell[0]+1)+k] * mw[k] / mw[i];
//	    }
//	    csr_jac_cell[k*(udata->neqs_per_cell[0]+1)+udata->neqs_per_cell[0]] = Jmat_pt[udata->neqs_per_cell[0]*(udata->neqs_per_cell[0]+1)+k] * mw[k]; 
//        }
//        for (int i = 0; i < udata->neqs_per_cell[0]; i++){
//            csr_jac_cell[udata->neqs_per_cell[0]*(udata->neqs_per_cell[0]+1)+i] = Jmat_pt[i*(udata->neqs_per_cell[0]+1)+udata->neqs_per_cell[0]] / mw[i]; 
//        }
//    }
//}


/* Get and print some final statistics */
static void PrintFinalStats(void *cvodeMem)
{
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
  int flag;

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

  flag = CVDlsGetNumJacEvals(cvodeMem, &nje);
  check_flag(&flag, "CVDlsGetNumJacEvals", 1);
  flag = CVDlsGetNumRhsEvals(cvodeMem, &nfeLS);
  check_flag(&flag, "CVDlsGetNumRhsEvals", 1);

  flag = CVodeGetNumGEvals(cvodeMem, &nge);
  check_flag(&flag, "CVodeGetNumGEvals", 1);

  printf("\nFinal Statistics:\n");
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld\n",
	 nst, nfe, nsetups, nfeLS, nje);
  printf("nni = %-6ld ncfn = %-6ld netf = %-6ld nge = %ld\n \n",
	 nni, ncfn, netf, nge);
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
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
              funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
            funcname);
    return(1); }

  return(0);
}

/* End of file  */
