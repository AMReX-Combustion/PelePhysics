#include <CPU/actual_Creactor.h> 
#include <AMReX_ParmParse.H>
#include <chemistry_file.H>
#include "mechanism.h"
#include <eos.H>

/**********************************/
/* Global Variables */
  N_Vector y         = NULL;
  SUNLinearSolver LS = NULL;
  SUNMatrix A        = NULL;
  void *cvode_mem    = NULL;
  /* User data */
  UserData data      = NULL;
/* OPTIONS */
  /* energy */
  double *rhoX_init   = NULL;
  double *rhoXsrc_ext = NULL;
  double *rYsrc       = NULL;

#ifdef _OPENMP
#pragma omp threadprivate(y,LS,A)
#pragma omp threadprivate(cvode_mem,data)
#pragma omp threadprivate(rhoX_init,rhoXsrc_ext,rYsrc)
#endif
/**********************************/

/**********************************/
/* Initialization routine, called once at the begining of the problem */
int reactor_init(const int* cvode_iE, const int* Ncells) {
        /* CVODE return Flag  */
	int flag;
	/* CVODE initial time - 0 */
	realtype time;
	/* CVODE tolerances */
	realtype reltol;
	N_Vector atol;
	realtype *ratol;
	/* Tot numb of eq to integrate */
	int neq_tot;
#ifdef _OPENMP
	int omp_thread;

	/* omp thread if applicable */
	omp_thread = omp_get_thread_num(); 
#endif
	/* Total number of eq to integrate */
        neq_tot        = (NUM_SPECIES + 1) * (*Ncells);

	/* Definition of main vector */
	y = N_VNew_Serial(neq_tot);
	if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);

	/* Call CVodeCreate to create the solver memory and specify the
	 * Backward Differentiation Formula and the use of a Newton iteration */
	cvode_mem = CVodeCreate(CV_BDF);
	if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);

        /* Does not work for more than 1 cell right now */
	data = AllocUserData(*cvode_iE, *Ncells);
	if(check_flag((void *)data, "AllocUserData", 2)) return(1);

	/* Nb of species and cells in mechanism */
#ifdef _OPENMP
        if ((data->iverbose > 0) && (omp_thread == 0)) {
#else
        if (data->iverbose > 0) {
#endif
		amrex::Print() << "Nb of spec in mech is " << NUM_SPECIES << "\n";    
		amrex::Print() << "Ncells in one solve is " << data->ncells << "\n";
	}

	/* Set the pointer to user-defined data */
	flag = CVodeSetUserData(cvode_mem, data);
	if(check_flag(&flag, "CVodeSetUserData", 1)) return(1);   

        time = 0.0e+0;
	/* Call CVodeInit to initialize the integrator memory and specify the
	 * user's right hand side function, the inital time, and 
	 * initial dependent variable vector y. */
	flag = CVodeInit(cvode_mem, cF_RHS, time, y);
	if (check_flag(&flag, "CVodeInit", 1)) return(1);
	
	/* Definition of tolerances: one for each species */
	/* TODO in fct of variable !! */
	reltol = 1.0e-10;
        atol  = N_VNew_Serial(neq_tot);
	ratol = N_VGetArrayPointer(atol);
        for (int i=0; i<neq_tot; i++) {
            ratol[i] = 1.0e-10;
        }
	/* Call CVodeSVtolerances to specify the scalar relative tolerance
	 * and vector absolute tolerances */
	flag = CVodeSVtolerances(cvode_mem, reltol, atol);
	if (check_flag(&flag, "CVodeSVtolerances", 1)) return(1);

	//flag = CVodeSetNonlinConvCoef(cvode_mem, 1.0e-1);
	//if (check_flag(&flag, "CVodeSetNonlinConvCoef", 1)) return(1);

	flag = CVodeSetMaxNonlinIters(cvode_mem, 50);
	if (check_flag(&flag, "CVodeSetMaxNonlinIters", 1)) return(1);

	flag = CVodeSetMaxErrTestFails(cvode_mem, 100);
	if (check_flag(&flag, "CVodeSetMaxErrTestFails", 1)) return(1);

	if (data->iDense_Creact == 1) {
#ifdef _OPENMP
            if ((data->iverbose > 0) && (omp_thread == 0)) {
#else
            if (data->iverbose > 0) {
#endif
	    	amrex::Print() << "\n--> Using a Direct Dense Solver\n";    
	    }

            /* Create dense SUNMatrix for use in linear solves */
	    A = SUNDenseMatrix(neq_tot, neq_tot);
	    if(check_flag((void *)A, "SUNDenseMatrix", 0)) return(1);

	    /* Create dense SUNLinearSolver object for use by CVode */
	    LS = SUNDenseLinearSolver(y, A);
	    if(check_flag((void *)LS, "SUNDenseLinearSolver", 0)) return(1);

	    /* Call CVDlsSetLinearSolver to attach the matrix and linear solver to CVode */
	    flag = CVDlsSetLinearSolver(cvode_mem, LS, A);
	    if(check_flag(&flag, "CVDlsSetLinearSolver", 1)) return(1);

	} else if (data->iDense_Creact == 5) {
#ifdef USE_KLU_PP 
#ifdef _OPENMP
            if ((data->iverbose > 0) && (omp_thread == 0)) {
#else
            if (data->iverbose > 0) {
#endif
	    	amrex::Print() << "\n--> Using a Direct Sparse Solver\n";    
	    }
	    /* Create sparse SUNMatrix for use in linear solves */
	    A = SUNSparseMatrix(neq_tot, neq_tot, (data->NNZ)*data->ncells, CSC_MAT);
            if(check_flag((void *)A, "SUNSparseMatrix", 0)) return(1);

	    /* Create KLU solver object for use by CVode */
	    LS = SUNLinSol_KLU(y, A);
	    if(check_flag((void *)LS, "SUNLinSol_KLU", 0)) return(1);

	    /* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */
	    flag = CVodeSetLinearSolver(cvode_mem, LS, A);
	    if(check_flag(&flag, "CVodeSetLinearSolver", 1)) return(1);
#else
            if (data->iverbose > 0) {
                amrex::Abort("Sparse solver not valid without KLU solver.");
	    }
#endif

	} else if (data->iDense_Creact == 99) {
#ifdef _OPENMP
            if ((data->iverbose > 0) && (omp_thread == 0)) {
#else
            if (data->iverbose > 0) {
#endif
	    	amrex::Print() << "\n--> Using an Iterative Solver\n";    
	    }

            /* Create the linear solver object */
	    if (data->iJac_Creact == 0) {
	        LS = SUNSPGMR(y, PREC_NONE, 0);
	        if(check_flag((void *)LS, "SUNDenseLinearSolver", 0)) return(1);
	    } else {
	        LS = SUNSPGMR(y, PREC_LEFT, 0);
	        if(check_flag((void *)LS, "SUNDenseLinearSolver", 0)) return(1);
	    }

	    /* Set CVSpils linear solver to LS */
	    flag = CVSpilsSetLinearSolver(cvode_mem, LS);
	    if(check_flag(&flag, "CVSpilsSetLinearSolver", 1)) return(1);
	} else {
            if (data->iverbose > 0) {
                amrex::Abort("Linear solvers availables are: Direct Dense (1), Direct Sparse (5) or Iterative GMRES (99)");
	    }
	}

	if (data->iJac_Creact == 0) {
#ifdef _OPENMP
            if ((data->iverbose > 0) && (omp_thread == 0)) {
#else
            if (data->iverbose > 0) {
#endif
		    amrex::Print() << "    Without Analytical J/Preconditioner\n";
	    }
#ifdef USE_KLU_PP 
	    if (data->iDense_Creact == 5) {
		if (data->iverbose > 0) {
	            amrex::Abort("A Sparse Solver should have an Analytical J");
		}
	    }
#endif
	} else {
	    if (data->iDense_Creact == 99) {
	        /* Set the JAcobian-times-vector function */
	        flag = CVSpilsSetJacTimes(cvode_mem, NULL, NULL);
	        if(check_flag(&flag, "CVSpilsSetJacTimes", 1)) return(1);
#ifdef USE_KLU_PP 
#ifdef _OPENMP
                if ((data->iverbose > 0) && (omp_thread == 0)) {
#else
                if (data->iverbose > 0) {
#endif
			amrex::Print() << "    With a Sparse Preconditioner\n";
		}
	        /* Set the preconditioner solve and setup functions */
	        flag = CVSpilsSetPreconditioner(cvode_mem, Precond_sparse, PSolve_sparse);
	        if(check_flag(&flag, "CVSpilsSetPreconditioner", 1)) return(1);
#else
#ifdef _OPENMP
                if ((data->iverbose > 0) && (omp_thread == 0)) {
#else
                if (data->iverbose > 0) {
#endif
			amrex::Print() << "    With a Preconditioner\n";
		}
	        /* Set the preconditioner solve and setup functions */
	        flag = CVSpilsSetPreconditioner(cvode_mem, Precond, PSolve);
	        if(check_flag(&flag, "CVSpilsSetPreconditioner", 1)) return(1);
#endif
#ifdef USE_KLU_PP 
	    } else if (data->iDense_Creact == 5){
#ifdef _OPENMP
                if ((data->iverbose > 0) && (omp_thread == 0)) {
#else
                if (data->iverbose > 0) {
#endif
                    amrex::Print() << "    With a Sparse Analytical J\n";
	        }
		/* Set the user-supplied Jacobian routine Jac */
		flag = CVodeSetJacFn(cvode_mem, cJac_KLU);
		if(check_flag(&flag, "CVodeSetJacFn", 1)) return(1); 
#endif
	    } else {
#ifdef _OPENMP
                if ((data->iverbose > 0) && (omp_thread == 0)) {
#else
                if (data->iverbose > 0) {
#endif
                    amrex::Print() << "    With Analytical J\n";
	        }
	        /* Set the user-supplied Jacobian routine Jac */
                flag = CVodeSetJacFn(cvode_mem, cJac);
		if(check_flag(&flag, "CVodeSetJacFn", 1)) return(1);
	    }
	}

        /* Set the max number of time steps */ 
	flag = CVodeSetMaxNumSteps(cvode_mem, 10000);
	if(check_flag(&flag, "CVodeSetMaxNumSteps", 1)) return(1);

        /* Set the max order */ 
        flag = CVodeSetMaxOrd(cvode_mem, 2);
	if(check_flag(&flag, "CVodeSetMaxOrd", 1)) return(1);

        /* Set the num of steps to wait inbetween Jac evals */ 
	flag = CVodeSetMaxStepsBetweenJac(cvode_mem, 100);
	if(check_flag(&flag, "CVodeSetMaxStepsBetweenJac", 1)) return(1);

	/* Define vectors to be used later in creact */
	rhoX_init = (double *) malloc(data->ncells*sizeof(double));
	rhoXsrc_ext = (double *) malloc( data->ncells*sizeof(double));
	rYsrc       = (double *)  malloc((data->ncells*NUM_SPECIES)*sizeof(double));

	/* Free the atol vector */
	N_VDestroy(atol); 

	/* Ok we're done ...*/
#ifdef _OPENMP
        if ((data->iverbose > 1) && (omp_thread == 0)) {
#else
        if (data->iverbose > 1) {
#endif
		amrex::Print() << "\n--> DONE WITH INITIALIZATION (CPU)" << data->iE_Creact << "\n";
	}

	/* Reactor is now initialized */
	data->reactor_cvode_initialized = true;

	return(0);
}


/* Main CVODE call routine */
int react(realtype *rY_in, realtype *rY_src_in, 
		realtype *rX_in, realtype *rX_src_in,
                realtype *dt_react, realtype *time){

	realtype time_init, time_out, dummy_time;
	int flag;
#ifdef _OPENMP
	int omp_thread;

	/* omp thread if applicable */
	omp_thread = omp_get_thread_num(); 
#endif

#ifdef _OPENMP
	if ((data->iverbose > 1) && (omp_thread == 0)) {
#else
        if (data->iverbose > 1) {
#endif
	    amrex::Print() <<"\n -------------------------------------\n";
	}

	/* Initial time and time to reach after integration */
        time_init = *time;
	time_out  = *time + (*dt_react);

#ifdef _OPENMP
	if ((data->iverbose > 3) && (omp_thread == 0)) {
#else
        if (data->iverbose > 3) {
#endif
	    amrex::Print() <<"BEG : time curr is "<< time_init << " and dt_react is " << *dt_react << " and final time should be " << time_out << "\n";
	}

	/* Get Device MemCpy of in arrays */
	/* Get Device pointer of solution vector */
	realtype *yvec_d      = N_VGetArrayPointer(y);
	/* rhoY,T */
	std::memcpy(yvec_d, rY_in, sizeof(realtype) * ((NUM_SPECIES+1)*data->ncells));
	/* rhoY_src_ext */
	std::memcpy(rYsrc, rY_src_in, (NUM_SPECIES * data->ncells)*sizeof(double));
	/* rhoE/rhoH */
	std::memcpy(rhoX_init, rX_in, sizeof(realtype) * data->ncells);
	std::memcpy(rhoXsrc_ext, rX_src_in, sizeof(realtype) * data->ncells);

	/* Check if y is within physical bounds */
	check_state(y);
	if (!(data->actual_ok_to_react))  { 
	    amrex::Abort("\n Check_state failed: state is out of react bounds \n");
	}

	/* ReInit CVODE is faster */
	CVodeReInit(cvode_mem, time_init, y);

	flag = CVode(cvode_mem, time_out, y, &dummy_time, CV_NORMAL);
	/* ONE STEP MODE FOR DEBUGGING */
	//flag = CVode(cvode_mem, time_out, y, &dummy_time, CV_ONE_STEP);
	if (check_flag(&flag, "CVode", 1)) return(1);

	/* Update dt_react with real time step taken ... */
	*dt_react = dummy_time - time_init;
#ifdef MOD_REACTOR
	/* If reactor mode is activated, update time */
	*time  = time_init + (*dt_react);
#endif

#ifdef _OPENMP
	if ((data->iverbose > 3) && (omp_thread == 0)) {
#else
        if (data->iverbose > 3) {
#endif
	    amrex::Print() <<"END : time curr is "<< dummy_time << " and actual dt_react is " << *dt_react << "\n";
	}

	/* Pack data to return in main routine external */
	std::memcpy(rY_in, yvec_d, ((NUM_SPECIES+1)*data->ncells)*sizeof(realtype));
	for  (int i = 0; i < data->ncells; i++) {
            rX_in[i] = rX_in[i] + (*dt_react) * rX_src_in[i];
	}

#ifdef _OPENMP
	if ((data->iverbose > 1) && (omp_thread == 0)) {
#else
	if (data->iverbose > 1) {
#endif
	    amrex::Print() <<"Additional verbose info --\n";
	    PrintFinalStats(cvode_mem, rY_in[NUM_SPECIES]);
	    amrex::Print() <<"\n -------------------------------------\n";
	}

	/* Get estimate of how hard the integration process was */
        long int nfe;
	flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
	return nfe;
}


/* RHS routine */
int cF_RHS(realtype t, N_Vector y_in, N_Vector ydot_in, 
		void *user_data){

	realtype *y_d      = N_VGetArrayPointer(y_in);
	realtype *ydot_d   = N_VGetArrayPointer(ydot_in);

        fKernelSpec(&t, y_d, ydot_d, 
		    user_data);

	return(0);
}
/**********************************/


/*
 * kernels
 */

/* RHS source terms evaluation */
void fKernelSpec(realtype *dt, realtype *yvec_d, realtype *ydot_d,  
			    void *user_data)
{
  /* Make local copies of pointers in user_data (cell M)*/
  UserData data_wk;
  data_wk = (UserData) user_data;   

  /* Tmp vars */
  int tid;

  /* Loop on packed cells */
  for (tid = 0; tid < data_wk->ncells; tid ++) {
      /* Tmp vars */
      realtype massfrac[NUM_SPECIES];
      realtype Xi[NUM_SPECIES];
      realtype cdot[NUM_SPECIES], molecular_weight[NUM_SPECIES];
      realtype cX;
      realtype temp, energy;
      /* EOS object in cpp */
      EOS eos;

      /* Offset in case several cells */
      int offset = tid * (NUM_SPECIES + 1); 
      
      /* MW CGS */
      CKWT(molecular_weight);

      /* rho MKS */ 
      realtype rho = 0.0;
      for (int i = 0; i < NUM_SPECIES; i++){
          rho = rho + yvec_d[offset + i];
      }

      /* temp */
      temp = yvec_d[offset + NUM_SPECIES];

      /* Yks */
      for (int i = 0; i < NUM_SPECIES; i++){
          massfrac[i] = yvec_d[offset + i] / rho;
      }

      /* NRG CGS */
      energy = (rhoX_init[tid] + rhoXsrc_ext[tid]*(*dt)) /rho;

      if (data_wk->iE_Creact == 1){
          /* UV REACTOR */
          eos.eos_EY2T(massfrac, energy, temp);
          eos.eos_TY2Cv(temp, massfrac, &cX);
          eos.eos_T2EI(temp, Xi);
      } else {
          /* HP REACTOR */
          eos.eos_HY2T(massfrac, energy, temp);
          eos.eos_TY2Cp(temp, massfrac, &cX);
          eos.eos_T2HI(temp, Xi);
      }
      eos.eos_RTY2W(rho, temp, massfrac, cdot);

      /* Fill ydot vect */
      ydot_d[offset + NUM_SPECIES] = rhoXsrc_ext[tid];
      for (int i = 0; i < NUM_SPECIES; i++){
          ydot_d[offset + i] = cdot[i] * molecular_weight[i] + rYsrc[tid * (NUM_SPECIES) + i];
          ydot_d[offset + NUM_SPECIES] = ydot_d[offset + NUM_SPECIES]  - ydot_d[offset + i] * Xi[i];
      }
      ydot_d[offset + NUM_SPECIES] = ydot_d[offset + NUM_SPECIES] /(rho * cX);
  }
}


/* Analytical Jacobian evaluation */
int cJac(realtype tn, N_Vector u, N_Vector fu, SUNMatrix J,
		void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){

  /* Make local copies of pointers to input data (big M) */
  realtype *ydata  = N_VGetArrayPointer(u);

  /* Make local copies of pointers in user_data (cell M)*/
  UserData data_wk;
  data_wk = (UserData) user_data;   

  int tid;
  for (tid = 0; tid < data_wk->ncells; tid ++) {
      /* Tmp vars */
      realtype *J_col_k;
      realtype massfrac[NUM_SPECIES], molecular_weight[NUM_SPECIES];
      realtype temp; 
      realtype Jmat_tmp[(NUM_SPECIES+1)*(NUM_SPECIES+1)];
      /* EOS object in cpp */
      EOS eos;

      /* Offset in case several cells */
      int offset = tid * (NUM_SPECIES + 1); 

      /* MW CGS */
      CKWT(molecular_weight);

      /* rho MKS */ 
      realtype rho = 0.0;
      for (int i = 0; i < NUM_SPECIES; i++){
          rho = rho + ydata[offset + i];
      }

      /* temp */
      temp = ydata[offset + NUM_SPECIES];

      /* Yks */
      for (int i = 0; i < NUM_SPECIES; i++){
          massfrac[i] = ydata[offset + i] / rho;
      }

      /* Jac */
      int consP;
      if (data_wk->iE_Creact == 1) {
	  consP = 0;
      } else {
          consP = 1;
      }
      eos.eos_RTY2JAC(rho, temp, massfrac, Jmat_tmp, consP);
      /* fill the sunMat */
      for (int k = 0; k < NUM_SPECIES; k++){
	  J_col_k = SM_COLUMN_D(J,offset + k);
	  for (int i = 0; i < NUM_SPECIES; i++){
	        J_col_k[offset + i] = Jmat_tmp[k*(NUM_SPECIES+1)+i] * molecular_weight[i] / molecular_weight[k]; 
          }
	  J_col_k[offset + NUM_SPECIES] = Jmat_tmp[k*(NUM_SPECIES+1)+NUM_SPECIES] / molecular_weight[k]; 
      }
      J_col_k = SM_COLUMN_D(J,offset + NUM_SPECIES);
      for (int i = 0; i < NUM_SPECIES; i++){
          J_col_k[offset + i] = Jmat_tmp[NUM_SPECIES*(NUM_SPECIES+1)+i] * molecular_weight[i]; 
      }
  }

  return(0);

}


#ifdef USE_KLU_PP 
/* Analytical SPARSE Jacobian evaluation */
int cJac_KLU(realtype tn, N_Vector u, N_Vector fu, SUNMatrix J,
		void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){

  /* Make local copies of pointers to input data (big M) */
  realtype *ydata           = N_VGetArrayPointer(u);
  sunindextype *colptrs_tmp = SUNSparseMatrix_IndexPointers(J);
  sunindextype *rowvals_tmp = SUNSparseMatrix_IndexValues(J);
  realtype *Jdata           = SUNSparseMatrix_Data(J);

  /* Make local copies of pointers in user_data (cell M)*/
  UserData data_wk;
  data_wk = (UserData) user_data;   

  /* MW CGS */
  realtype molecular_weight[NUM_SPECIES];
  CKWT(molecular_weight);

  /* Fixed RowVals */
  for (int i=0;i<data_wk->NNZ;i++) {
      rowvals_tmp[i] = data_wk->rowVals[0][i];
  }
  /* Fixed colPtrs */
  colptrs_tmp[0] = data_wk->colPtrs[0][0];
  for (int i=0;i<data_wk->ncells*(NUM_SPECIES + 1);i++) {
      colptrs_tmp[i+1] = data_wk->colPtrs[0][i+1];
  }

  /* Temp vectors */
  realtype temp_save_lcl, temp;
  realtype massfrac[NUM_SPECIES];
  realtype Jmat_tmp[(NUM_SPECIES+1)*(NUM_SPECIES+1)];
  /* EOS object in cpp */
  EOS eos;
  /* Idx for sparsity */
  int tid, offset, nbVals, idx;
  /* Save Jac from cell to cell if more than one */
  temp_save_lcl = 0.0;
  for (tid = 0; tid < data_wk->ncells; tid ++) {
      /* Offset in case several cells */
      offset = tid * (NUM_SPECIES + 1); 
      /* rho MKS */ 
      realtype rho = 0.0;
      for (int i = 0; i < NUM_SPECIES; i++){
          rho = rho + ydata[offset + i];
      }
      /* Yks */
      for (int i = 0; i < NUM_SPECIES; i++){
          massfrac[i] = ydata[offset + i] / rho;
      }
      /* temp */
      temp = ydata[offset + NUM_SPECIES];
      /* Do we recompute Jac ? */
      if (fabs(temp - temp_save_lcl) > 1.0) {
          /* NRG CGS */
          int consP;
          if (data_wk->iE_Creact == 1) {
              consP = 0;
          } else {
              consP = 1;
          }
          eos.eos_RTY2JAC(rho, temp, massfrac, Jmat_tmp, consP);
	  temp_save_lcl = temp;
	  /* rescale */
          for (int i = 0; i < NUM_SPECIES; i++) {
              for (int k = 0; k < NUM_SPECIES; k++) {
                  Jmat_tmp[k*(NUM_SPECIES+1) + i] = Jmat_tmp[k*(NUM_SPECIES+1) + i] * molecular_weight[i] / molecular_weight[k];
              }
              Jmat_tmp[i*(NUM_SPECIES+1) + NUM_SPECIES] = Jmat_tmp[i*(NUM_SPECIES+1) + NUM_SPECIES] / molecular_weight[i];
          }
          for (int i = 0; i < NUM_SPECIES; i++) {
              Jmat_tmp[NUM_SPECIES*(NUM_SPECIES+1) + i] = Jmat_tmp[NUM_SPECIES*(NUM_SPECIES+1) + i] * molecular_weight[i];
          }
      }
      /* Go from Dense to Sparse */
      for (int i = 1; i < NUM_SPECIES+2; i++) {
	  nbVals = data_wk->colPtrs[0][i]-data_wk->colPtrs[0][i - 1];
	  for (int j = 0; j < nbVals; j++) {
	          idx = data_wk->rowVals[0][ data_wk->colPtrs[0][i - 1] + j ];
	              Jdata[ data_wk->colPtrs[0][offset + i - 1] + j ] = Jmat_tmp[(i - 1) * (NUM_SPECIES + 1) + idx];
	  }
      }
  }

  return(0);

}
#endif


/* Jacobian-times-vector routine.
 * Currently not used !!
int jtv(N_Vector v, N_Vector Jv, realtype t, N_Vector u, N_Vector fu, 
			  void *user_data, N_Vector tmp)
{
  realtype *udata, *vdata, *Jvdata;
  vdata  = N_VGetArrayPointer(v); 
  udata  = N_VGetArrayPointer(u);
  Jvdata = N_VGetArrayPointer(Jv);

  int tid;
  for (tid = 0; tid < data->ncells; tid ++) {
	realtype temp;
	realtype activity[NUM_SPECIES], molecular_weight[NUM_SPECIES];
	realtype J[(NUM_SPECIES+1)*(NUM_SPECIES+1)];

        int offset = tid * (NUM_SPECIES + 1); 

        // MW CGS 
	CKWT(molecular_weight);
	for (int i = 0; i < NUM_SPECIES; i++){
            activity[i] = udata[offset + i]/(molecular_weight[i]);
	}
        // temp 
	temp = udata[offset + NUM_SPECIES];
        // NRG CGS 
        if (iE_Creact == 1) {
            int consP = 0;
            DWDOT(J, activity, &temp, &consP);
        } else {
	    int consP = 1 ;
            DWDOT(J, activity, &temp, &consP);
        }

	// PRINT JAC INFO: debug mode
	//for (int i = 0; i < NUM_SPECIES+1; i++){
	//    for (int j = 0; j < NUM_SPECIES+1; j++){
        //        printf(" %4.8e ", (J[j*(NUM_SPECIES+1)+i]));
	//    }
        //    printf(" \n");
	//}
	//amrex::Abort("\n--> ABORT\n");

	// reinit 
	for (int i = 0; i < NUM_SPECIES+1; i++){
	    Jvdata[offset + i] = 0.0;
	}
	// COmpute 
	for (int i = 0; i < NUM_SPECIES; i++){
            for (int j = 0; j < NUM_SPECIES; j++){
                Jvdata[offset + i] = Jvdata[offset + i] + J[j*(NUM_SPECIES+1)+i] * vdata[offset + j] * molecular_weight[i] / molecular_weight[j];
	    }
            Jvdata[offset + i] = Jvdata[offset + i] + J[NUM_SPECIES*(NUM_SPECIES+1)+i] * vdata[offset + NUM_SPECIES] * molecular_weight[i];
        }
	for (int j = 0; j < NUM_SPECIES; j++){
	    Jvdata[offset + NUM_SPECIES] = Jvdata[offset + NUM_SPECIES] + J[j*(NUM_SPECIES+1)+NUM_SPECIES] * vdata[offset + j] / molecular_weight[j];
	}
	Jvdata[offset + NUM_SPECIES] = Jvdata[offset + NUM_SPECIES]  + J[(NUM_SPECIES+1)*(NUM_SPECIES+1)-1] * vdata[offset + NUM_SPECIES] ;  
  }

  return(0);
}
 */


#ifdef USE_KLU_PP 
/* Preconditioner setup routine for GMRES solver when KLU sparse mode is activated 
 * Generate and preprocess P
*/
int Precond_sparse(realtype tn, N_Vector u, N_Vector fu, booleantype jok, 
		booleantype *jcurPtr, realtype gamma, void *user_data)
{
  /* Make local copies of pointers to input data (big M) */
  realtype *udata   = N_VGetArrayPointer(u);
  /* Make local copies of pointers in user_data (cell M)*/
  UserData data_wk;
  data_wk = (UserData) user_data;   
  /* Tmp array */
  int ok,tid;

  /* MW CGS */
  realtype molecular_weight[NUM_SPECIES];
  CKWT(molecular_weight);

  /* Check if Jac is stale */
  if (jok) {
        /* jok = SUNTRUE: Copy Jbd to P */
        *jcurPtr = SUNFALSE;
  } else {
        /* Temp vectors */
        realtype temp, temp_save_lcl;
        realtype activity[NUM_SPECIES], massfrac[NUM_SPECIES];
        /* EOS object in cpp */
        EOS eos;
        /* Save Jac from cell to cell if more than one */
        temp_save_lcl = 0.0;
        for (tid = 0; tid < data_wk->ncells; tid ++) {
            /* Offset in case several cells */
            int offset = tid * (NUM_SPECIES + 1); 
            /* rho MKS */ 
            realtype rho = 0.0;
            for (int i = 0; i < NUM_SPECIES; i++){
                rho = rho + udata[offset + i];
            }
            /* Yks */
            for (int i = 0; i < NUM_SPECIES; i++){
                massfrac[i] = udata[offset + i] / rho;
            }
            /* temp */
            temp = udata[offset + NUM_SPECIES];
            /* Activities */
	    eos.eos_RTY2C(rho, temp, massfrac, activity);
            /* Do we recompute Jac ? */
            if (fabs(temp - temp_save_lcl) > 1.0) {
                /* Formalism */
                int consP;
                if (data_wk->iE_Creact == 1) {
                    consP = 0;
                } else {
                    consP = 1;
                }
                DWDOT_PRECOND(data_wk->JSPSmat[tid], activity, &temp, &consP);

                for (int i = 0; i < NUM_SPECIES; i++) {
                    for (int k = 0; k < NUM_SPECIES; k++) {
                        (data_wk->JSPSmat[tid])[k*(NUM_SPECIES+1) + i] = (data_wk->JSPSmat[tid])[k*(NUM_SPECIES+1) + i] * molecular_weight[i] / molecular_weight[k];
                    }
                    (data_wk->JSPSmat[tid])[i*(NUM_SPECIES+1) + NUM_SPECIES] = (data_wk->JSPSmat[tid])[i*(NUM_SPECIES+1) + NUM_SPECIES] / molecular_weight[i];
                }
                for (int i = 0; i < NUM_SPECIES; i++) {
                    (data_wk->JSPSmat[tid])[NUM_SPECIES*(NUM_SPECIES+1) + i] = (data_wk->JSPSmat[tid])[NUM_SPECIES*(NUM_SPECIES+1) + i] * molecular_weight[i];
                }
	        temp_save_lcl = temp;
	    } else {
		/* if not: copy the one from prev cell */
		for (int i = 0; i < NUM_SPECIES+1; i++) {
		    for (int k = 0; k < NUM_SPECIES+1; k++) {
		        (data_wk->JSPSmat[tid])[k*(NUM_SPECIES+1) + i] = (data_wk->JSPSmat[tid-1])[k*(NUM_SPECIES+1) + i];
		    }
		}
	    }
	}

        *jcurPtr = SUNTRUE;
  }

  int nbVals;
  for (int i = 1; i < NUM_SPECIES+2; i++) {
      /* nb non zeros elem should be the same for all cells */
      nbVals = data_wk->colPtrs[0][i]-data_wk->colPtrs[0][i-1];
      //printf("Rows %d : we have %d nonzero values \n", i-1, nbVals);
      for (int j = 0; j < nbVals; j++) {
    	  /* row of non zero elem should be the same for all cells */
    	  int idx = data_wk->rowVals[0][ data_wk->colPtrs[0][i-1] + j ];
          /* Scale by -gamma */
          /* Add identity matrix */
          for (tid = 0; tid < data_wk->ncells; tid ++) {
    	      if (idx == (i-1)) {
                  data_wk->Jdata[tid][ data_wk->colPtrs[tid][i-1] + j ] = 1.0 - gamma * (data_wk->JSPSmat[tid])[ idx * (NUM_SPECIES+1) + idx]; 
    	      } else {
                  data_wk->Jdata[tid][ data_wk->colPtrs[tid][i-1] + j ] = - gamma * (data_wk->JSPSmat[tid])[ (i-1) * (NUM_SPECIES+1) + idx ]; 
    	      }
          }
      }
  }
  
  if (!(data_wk->FirstTimePrecond)) {
      for (tid = 0; tid < data_wk->ncells; tid ++) {
          ok = klu_refactor(data_wk->colPtrs[tid], data_wk->rowVals[tid], data_wk->Jdata[tid], data_wk->Symbolic[tid], data_wk->Numeric[tid], &(data_wk->Common[tid]));
      }
  } else {
      for (tid = 0; tid < data_wk->ncells; tid ++) {
          data_wk->Numeric[tid] = klu_factor(data_wk->colPtrs[tid], data_wk->rowVals[tid], data_wk->Jdata[tid], data_wk->Symbolic[tid], &(data_wk->Common[tid])) ; 
      }
      data_wk->FirstTimePrecond = false;
  }

  return(0);
}

#else 
/* Preconditioner setup routine for GMRES solver when no sparse mode is activated 
 * Generate and preprocess P
*/
int Precond(realtype tn, N_Vector u, N_Vector fu, booleantype jok, 
		booleantype *jcurPtr, realtype gamma, void *user_data)
{
  /* Make local copies of pointers to input data (big M) */
  realtype *udata = N_VGetArrayPointer(u);

  /* Make local copies of pointers in user_data */
  UserData data_wk;
  data_wk = (UserData) user_data;   
  realtype **(**P), **(**Jbd);
  sunindextype *(**pivot);
  P     = (data_wk->P);
  Jbd   = (data_wk->Jbd);
  pivot = (data_wk->pivot);

  /* Tmp arrays */
  realtype molecular_weight[NUM_SPECIES];
  realtype Jmat[(NUM_SPECIES+1)*(NUM_SPECIES+1)];
  realtype massfrac[NUM_SPECIES], activity[NUM_SPECIES];
  realtype temp;
  sunindextype ierr;

  /* MW CGS */
  CKWT(molecular_weight);

  if (jok) {
      /* jok = SUNTRUE: Copy Jbd to P */
      denseCopy(Jbd[0][0], P[0][0], NUM_SPECIES+1, NUM_SPECIES+1);
      *jcurPtr = SUNFALSE;
  } else {
      /* EOS object in cpp */   
      EOS eos;
      /* rho MKS */ 
      realtype rho = 0.0;
      for (int i = 0; i < NUM_SPECIES; i++){
          rho = rho + udata[i];
      }
      /* Yks */
      for (int i = 0; i < NUM_SPECIES; i++){
          massfrac[i] = udata[i] / rho;
      }
      /* temp */
      temp = udata[NUM_SPECIES];
      /* Activities */
      eos.eos_RTY2C(rho, temp, massfrac, activity);
      /* jok = SUNFALSE: Generate Jbd from scratch and copy to P */
      /* Make local copies of problem variables, for efficiency. */
      int consP;
      if (data_wk->iE_Creact == 1) { 
          consP = 0;
      } else {
          consP = 1;
      }
      DWDOT_PRECOND(Jmat, activity, &temp, &consP);

      /* Compute Jacobian.  Load into P. */
      denseScale(0.0, Jbd[0][0], NUM_SPECIES+1, NUM_SPECIES+1);
      for (int i = 0; i < NUM_SPECIES; i++) {
          for (int k = 0; k < NUM_SPECIES; k++) {
              (Jbd[0][0])[k][i] = Jmat[k*(NUM_SPECIES+1) + i] * molecular_weight[i] / molecular_weight[k];
          }
          (Jbd[0][0])[i][NUM_SPECIES] = Jmat[i*(NUM_SPECIES+1) + NUM_SPECIES] / molecular_weight[i];
      }
      for (int i = 0; i < NUM_SPECIES; i++) {
          (Jbd[0][0])[NUM_SPECIES][i] = Jmat[NUM_SPECIES*(NUM_SPECIES+1) + i] * molecular_weight[i];
      }
      (Jbd[0][0])[NUM_SPECIES][NUM_SPECIES] = Jmat[(NUM_SPECIES+1)*(NUM_SPECIES+1)-1];

      denseCopy(Jbd[0][0], P[0][0], NUM_SPECIES+1, NUM_SPECIES+1);

      *jcurPtr = SUNTRUE;
  }
  
  /* Scale by -gamma */
  denseScale(-gamma, P[0][0], NUM_SPECIES+1, NUM_SPECIES+1);
  //denseScale(0.0, P[0][0], NUM_SPECIES+1, NUM_SPECIES+1);
  
  /* Add identity matrix and do LU decompositions on blocks in place. */
  denseAddIdentity(P[0][0], NUM_SPECIES+1);
  ierr = denseGETRF(P[0][0], NUM_SPECIES+1, NUM_SPECIES+1, pivot[0][0]);
  if (ierr != 0) return(1);
  
  return(0);
}
#endif


#ifdef USE_KLU_PP 
/* PSolve for GMRES solver when KLU sparse mode is activated */
int PSolve_sparse(realtype tn, N_Vector u, N_Vector fu, N_Vector r, N_Vector z,
                  realtype gamma, realtype delta, int lr, void *user_data)
{
  /* Make local copies of pointers in user_data */
  UserData data_wk;
  data_wk = (UserData) user_data;

  /* Make local copies of pointers to input data (big M) */
  realtype *zdata = N_VGetArrayPointer(z);

  N_VScale(1.0, r, z);

  /* Solve the block-diagonal system Pz = r using LU factors stored
     in P and pivot data in pivot, and return the solution in z. */
  int tid, offset_beg, offset_end;
  realtype zdata_cell[NUM_SPECIES+1];
  for (tid = 0; tid < data_wk->ncells; tid ++) {
      offset_beg = tid * (NUM_SPECIES + 1); 
      offset_end = (tid + 1) * (NUM_SPECIES + 1);
      std::memcpy(zdata_cell, zdata+offset_beg, (NUM_SPECIES+1)*sizeof(realtype));
      klu_solve(data_wk->Symbolic[tid], data_wk->Numeric[tid], NUM_SPECIES+1, 1, zdata_cell, &(data_wk->Common[tid])) ; 
      std::memcpy(zdata+offset_beg, zdata_cell, (NUM_SPECIES+1)*sizeof(realtype));
  }

  return(0);
}

#else
/* PSolve for GMRES solver when no sparse mode is activated */
int PSolve(realtype tn, N_Vector u, N_Vector fu, N_Vector r, N_Vector z,
                  realtype gamma, realtype delta, int lr, void *user_data)
{
  /* Make local copies of pointers to input data (big M) */
  realtype *zdata = N_VGetArrayPointer(z);

  /* Extract the P and pivot arrays from user_data. */
  UserData data_wk;
  data_wk = (UserData) user_data;
  realtype **(**P);
  sunindextype *(**pivot);
  P     = data_wk->P;
  pivot = data_wk->pivot;
  
  N_VScale(1.0, r, z);
  
  /* Solve the block-diagonal system Pz = r using LU factors stored
     in P and pivot data in pivot, and return the solution in z. */
  realtype *v = zdata;
  denseGETRS(P[0][0], NUM_SPECIES+1, pivot[0][0], v);

  return(0);
}
#endif


/* 
 * OTHERS
*/

void check_state(N_Vector yvec) 
{
  realtype *ydata = N_VGetArrayPointer(yvec);

  data->actual_ok_to_react = true;

  realtype Temp;
  int offset;
  for (int tid = 0; tid < data->ncells; tid ++) {
      /* Offset in case several cells */
      offset = tid * (NUM_SPECIES + 1); 
      /* rho MKS */ 
      realtype rho = 0.0;
      for (int k = 0; k < NUM_SPECIES; k ++) {
          rho =  rho + ydata[offset + k];
      }
      /* temp */
      Temp = ydata[offset + NUM_SPECIES];
      if ((rho < 1.0e-10) || (rho > 1.e10)) {
          data->actual_ok_to_react = false;
      }
      if ((Temp < 200.0) || (Temp > 5000.0)) {
          data->actual_ok_to_react = false; 
      }
  }

}

/* Get and print some final statistics */
void PrintFinalStats(void *cvodeMem, realtype Temp)
{
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn;
  long int nli, npe, nps, ncfl, netfails;
  int flag;
  realtype hlast, hinused, hcur;

  flag = CVodeGetNumSteps(cvodeMem, &nst);
  check_flag(&flag, "CVodeGetNumSteps", 1);
  flag = CVodeGetNumRhsEvals(cvodeMem, &nfe);
  check_flag(&flag, "CVodeGetNumRhsEvals", 1);
  flag = CVodeGetNumNonlinSolvIters(cvodeMem, &nni);
  check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumLinSolvSetups(cvodeMem, &nsetups);
  check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvodeMem, &ncfn);
  check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);
  flag = CVodeGetNumErrTestFails(cvodeMem, &netfails);
  check_flag(&flag, "CVodeGetNumErrTestFails", 1);
  flag = CVodeGetLastStep(cvodeMem, &hlast);
  check_flag(&flag, "CVodeGetLastStep", 1);
  flag = CVodeGetActualInitStep(cvodeMem, &hinused);
  check_flag(&flag, "CVodeGetActualInitStep", 1);
  flag = CVodeGetCurrentTime(cvodeMem, &hcur);
  check_flag(&flag, "CVodeGetCurrentTime", 1);

  if (data->iDense_Creact == 1){
      flag = CVDlsGetNumRhsEvals(cvodeMem, &nfeLS);
      check_flag(&flag, "CVDlsGetNumRhsEvals", 1);
      flag = CVDlsGetNumJacEvals(cvodeMem, &nje);
      check_flag(&flag, "CVDlsGetNumJacEvals", 1);
  } else if (data->iDense_Creact == 99){
      flag = CVSpilsGetNumRhsEvals(cvode_mem, &nfeLS);
      check_flag(&flag, "CVSpilsGetNumRhsEvals", 1);
      flag = CVSpilsGetNumJtimesEvals(cvodeMem, &nje);
      //flag = CVSpilsGetNumJTSetupEvals(cvodeMem, &nje);
      check_flag(&flag, "CVSpilsGetNumJTSetupEvals", 1);
      flag = CVSpilsGetNumPrecEvals(cvodeMem, &npe);
      check_flag(&flag, "CVSpilsGetNumPrecEvals", 1);
      flag = CVSpilsGetNumPrecSolves(cvodeMem, &nps);
      check_flag(&flag, "CVSpilsGetNumPrecSolves", 1);
      flag = CVSpilsGetNumLinIters(cvodeMem, &nli);
      check_flag(&flag, "CVSpilsGetNumLinIters", 1);
      flag = CVSpilsGetNumConvFails(cvodeMem, &ncfl); 
      check_flag(&flag, "CVSpilsGetNumConvFails", 1);
  }

  amrex::Print() << "-- Final Statistics --\n";
  amrex::Print() << "NonLinear (Newton) related --\n";
  amrex::Print() << Temp << " |DT(dt, dtcur) = " << nst << "(" << hlast << "," << hcur << "), RHS = " << nfe << ", Iterations = " << nni << ", ErrTestFails = " << netfails << ", LinSolvSetups = " << nsetups << "\n";
  if (data->iDense_Creact == 1){
	  amrex::Print() <<"Linear (Dense Direct Solve) related --\n";
	  amrex::Print()<<Temp << " |FD RHS = "<< nfeLS<<", NumJacEvals = "<< nje <<" \n";
  } else if (data->iDense_Creact == 99){
	  // LinSolvSetups actually reflects the number of time the LinSolver has been called. 
	  // NonLinIterations can be taken without the need for LinItes
      amrex::Print() << "Linear (Krylov GMRES Solve) related --\n";
      amrex::Print() << Temp << " |RHSeval = "<< nfeLS << ", jtvEval = "<<nje << ", NumPrecEvals = "<< npe << ", NumPrecSolves = "<< nps <<"\n"; 
      amrex::Print() <<Temp << " |Iterations = "<< nli <<", ConvFails = "<< ncfl<<"\n"; 
  }
}


/* Check function return value...
     opt == 0 means SUNDIALS function allocates memory so check if
              returned NULL pointer
     opt == 1 means SUNDIALS function returns a flag so check if
              flag >= 0
     opt == 2 means function allocates memory so check if returned
              NULL pointer */

int check_flag(void *flagvalue, const char *funcname, int opt)
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


/* Alloc Data for CVODE */
UserData AllocUserData(int iE, int num_cells)
{
  /* Make local copies of pointers in user_data */
  UserData data_wk;
  data_wk = (UserData) malloc(sizeof *data_wk);
#ifdef _OPENMP
  int omp_thread;

  /* omp thread if applicable */
  omp_thread = omp_get_thread_num(); 
#endif

  /* ParmParse from the inputs file */
  /* TODO change that in the future */ 
  amrex::ParmParse pp("ns");
  pp.query("cvode_iJac",data_wk->iJac_Creact);
  pp.query("cvode_iDense", data_wk->iDense_Creact);
  (data_wk->iE_Creact)      = iE;

  (data_wk->ncells)                    = num_cells;

  (data_wk->iverbose)                  = 1;

  (data_wk->FirstTimePrecond)          = true;
  (data_wk->reactor_cvode_initialized) = false;
  (data_wk->actual_ok_to_react)        = true; 

#ifndef USE_KLU_PP
  if (data_wk->iDense_Creact == 99) {
      /* Precond data */
      (data_wk->P)     = new realtype***[data_wk->ncells];
      (data_wk->Jbd)   = new realtype***[data_wk->ncells];
      (data_wk->pivot) = new sunindextype**[data_wk->ncells];
      for(int i = 0; i < data_wk->ncells; ++i) {
              (data_wk->P)[i]     = new realtype**[data_wk->ncells];
              (data_wk->Jbd)[i]   = new realtype**[data_wk->ncells];
              (data_wk->pivot)[i] = new sunindextype*[data_wk->ncells];
      }

      for(int i = 0; i < data_wk->ncells; ++i) {
          (data_wk->P)[i][i]     = newDenseMat(NUM_SPECIES+1, NUM_SPECIES+1);
          (data_wk->Jbd)[i][i]   = newDenseMat(NUM_SPECIES+1, NUM_SPECIES+1);
          (data_wk->pivot)[i][i] = newIndexArray(NUM_SPECIES+1);
      }
  } 

#else
  /* Sparse Direct and Sparse (It) Precond data */
  data_wk->colPtrs = new int*[data_wk->ncells];
  data_wk->rowVals = new int*[data_wk->ncells];
  data_wk->Jdata   = new realtype*[data_wk->ncells];

  int HP;
  if (data_wk->iE_Creact == 1) {
      HP = 0;
  } else {
      HP = 1;
  }
  if (data_wk->iDense_Creact == 5) {
      /* Sparse Matrix for Direct Sparse KLU solver */
      (data_wk->PS) = new SUNMatrix[1];
      SPARSITY_INFO(&(data_wk->NNZ),&HP,data_wk->ncells);
#ifdef _OPENMP
      if ((data_wk->iverbose > 0) && (omp_thread == 0)) {
#else
      if (data_wk->iverbose > 0) {
#endif
          amrex::Print() << "--> SPARSE solver -- non zero entries: " << data_wk->NNZ << ", which represents "<< data_wk->NNZ/float((NUM_SPECIES+1) * (NUM_SPECIES+1) * (data_wk->ncells) * (data_wk->ncells)) *100.0 <<" % fill-in pattern\n";
      }
      (data_wk->PS)[0] = SUNSparseMatrix((NUM_SPECIES+1)*data_wk->ncells, (NUM_SPECIES+1)*data_wk->ncells, data_wk->NNZ, CSC_MAT);
      data_wk->colPtrs[0] = (int*) SUNSparseMatrix_IndexPointers((data_wk->PS)[0]); 
      data_wk->rowVals[0] = (int*) SUNSparseMatrix_IndexValues((data_wk->PS)[0]);
      data_wk->Jdata[0] = SUNSparseMatrix_Data((data_wk->PS)[0]);
      SPARSITY_PREPROC(data_wk->rowVals[0],data_wk->colPtrs[0],&HP,data_wk->ncells);

  } else if (data_wk->iDense_Creact == 99) {
      /* KLU internal storage */
      data_wk->Common   = new klu_common[data_wk->ncells];
      data_wk->Symbolic = new klu_symbolic*[data_wk->ncells];
      data_wk->Numeric  = new klu_numeric*[data_wk->ncells];
      /* Sparse Matrices for It Sparse KLU block-solve */
      data_wk->PS = new SUNMatrix[data_wk->ncells];
      /* Nb of non zero elements*/
      SPARSITY_INFO_PRECOND(&(data_wk->NNZ),&HP);
#ifdef _OPENMP
      if ((data_wk->iverbose > 0) && (omp_thread == 0) && (data_wk->iJac_Creact != 0)) {
#else
      if ((data_wk->iverbose > 0) && (data_wk->iJac_Creact != 0)) {
#endif
          amrex::Print() << "--> SPARSE Preconditioner -- non zero entries: " << data_wk->NNZ << ", which represents "<< data_wk->NNZ/float((NUM_SPECIES+1) * (NUM_SPECIES+1)) *100.0 <<" % fill-in pattern\n";
      }
      /* Not used yet. TODO use to fetch sparse Mat */
      data_wk->indx      = new int[data_wk->NNZ];
      data_wk->JSPSmat   = new realtype*[data_wk->ncells];
      for(int i = 0; i < data_wk->ncells; ++i) {
          (data_wk->PS)[i]    = SUNSparseMatrix(NUM_SPECIES+1, NUM_SPECIES+1, data_wk->NNZ, CSC_MAT);
          data_wk->colPtrs[i] = (int*) SUNSparseMatrix_IndexPointers((data_wk->PS)[i]); 
          data_wk->rowVals[i] = (int*) SUNSparseMatrix_IndexValues((data_wk->PS)[i]);
          data_wk->Jdata[i]   = SUNSparseMatrix_Data((data_wk->PS)[i]);
	  /* indx not used YET */
          SPARSITY_PREPROC_PRECOND(data_wk->rowVals[i],data_wk->colPtrs[i],data_wk->indx,&HP);
          data_wk->JSPSmat[i] = new realtype[(NUM_SPECIES+1)*(NUM_SPECIES+1)];
          klu_defaults (&(data_wk->Common[i]));
          //data_wk->Common.btf = 0;
          //(data_wk->Common[i]).maxwork = 15;
          //data_wk->Common.ordering = 1;
          data_wk->Symbolic[i] = klu_analyze (NUM_SPECIES+1, data_wk->colPtrs[i], data_wk->rowVals[i], &(data_wk->Common[i])) ; 
      }
  }
#endif

  return(data_wk);
}


/* Free memory */
void reactor_close(){

  CVodeFree(&cvode_mem);
  SUNLinSolFree(LS);

  if (data->iDense_Creact == 1) {
    SUNMatDestroy(A);
  }

  N_VDestroy(y); 
  FreeUserData(data);

  free(rhoX_init);
  free(rhoXsrc_ext);
  free(rYsrc);
}


/* Free data memory 
 * Probably not complete, how about the stuff allocated in KLU mode ? */
void FreeUserData(UserData data_wk)
{
#ifndef USE_KLU_PP
  if (data_wk->iDense_Creact == 99) {
      for(int i = 0; i < data_wk->ncells; ++i) {
          destroyMat((data_wk->P)[i][i]);
          destroyMat((data_wk->Jbd)[i][i]);
          destroyArray((data_wk->pivot)[i][i]);
      }
  }
#else
  free(data_wk->colPtrs);
  free(data_wk->rowVals);
  free(data_wk->Jdata);
  // Destroy data_wk->PS ?
#endif
  free(data_wk);
} 









