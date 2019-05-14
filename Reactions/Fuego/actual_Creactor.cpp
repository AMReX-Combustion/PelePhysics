#include <actual_Creactor.h> 
#include <AMReX_ParmParse.H>
#include <chemistry_file.H>

/**********************************/

/* Global Variables */
/* FOR CVODE */
  N_Vector y         = NULL;
  SUNLinearSolver LS = NULL;
  SUNMatrix A        = NULL;
  void *cvode_mem    = NULL;
  /* counters */
  long int nsetups_old = 0;
  long int nst_old     = 0;
  long int nfe_old     = 0;
  long int nni_old     = 0;
  long int ncfn_old    = 0;
  long int nfeLS_old   = 0;
  long int nje_old = 0;
  long int npe_old = 0;
  long int nps_old = 0;
  /* User data */
  UserData data = NULL;
  int NEQ       = 0;
  int NCELLS    = 0;
/* OPTIONS */
  /* base */
  int iDense_Creact = 1;
  int iJac_Creact   = 0;
  int iE_Creact     = 1;
  int iverbose      = 1;
  /* energy */
  double *rhoe_init = NULL;
  double *rhoh_init = NULL;
  double *rhoesrc_ext = NULL;
  double *rhohsrc_ext = NULL;
  double *rYsrc       = NULL;
  /* hacks */
  double temp_old  = 0;
  double temp_save = 0;
  bool InitPartial      = false;
  bool FirstTimePrecond = true;
  /* Checks */
  bool reactor_cvode_initialized = false;
  bool actual_ok_to_react = true;
/* TIMERS */
  std::chrono::duration<double> elapsed_seconds;
  std::chrono::duration<double> elapsed_seconds_RHS;
  std::chrono::duration<double> elapsed_seconds_Pcond;
  std::chrono::duration<double> elapsed_seconds_JacFuego;

/**********************************/
/* Definitions */
/* Initialization routine, called once at the begining of the problem */
int reactor_init(const int* cvode_iE, const int* Ncells) {

	int flag;
	realtype reltol, time;
	N_Vector atol;
	realtype *ratol;
	int mm, ii, nfit;
	int neq_tot;

	CKINDX(&mm,&NEQ,&ii,&nfit);
        if (iverbose > 0) {
	    printf("Nb of spec is %d \n", NEQ);
	}

	/* ParmParse from the inputs file */ 
	amrex::ParmParse pp("ns");
	pp.query("cvode_iJac",iJac_Creact);
	pp.query("cvode_iDense", iDense_Creact);

	/* Args */
	iE_Creact      = *cvode_iE;
	NCELLS         = *Ncells;
        neq_tot        = (NEQ + 1) * NCELLS;

        if (iverbose > 0) {
	    printf("Ncells in one solve ? %d\n",NCELLS);
	}

	/* Definition of main vector */
	y = N_VNew_Serial(neq_tot);
	if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);

	/* Call CVodeCreate to create the solver memory and specify the
	 * Backward Differentiation Formula and the use of a Newton iteration */
	//if ((*cvode_meth == 2) && (*cvode_itmeth == 2))
	//{
	cvode_mem = CVodeCreate(CV_BDF);
	if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);
	//} else {
	//    amrex::Abort("\n--> Weird inputs to CVodeCreate. Viable options are CV_BDF (=2), CV_NEWTON (=2)\n");
	//}

        /* Does not work for more than 1 cell right now */
	data = AllocUserData();
	if(check_flag((void *)data, "AllocUserData", 2)) return(1);

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

	flag = CVodeSetInitStep(cvode_mem, 1.0e-09);
	if (check_flag(&flag, "CVodeSetInitStep", 1)) return(1);

	//flag = CVodeSetNonlinConvCoef(cvode_mem, 1.0e-03);
	//if (check_flag(&flag, "CVodeSetNonlinConvCoef", 1)) return(1);

	flag = CVodeSetMaxNonlinIters(cvode_mem, 50);
	if (check_flag(&flag, "CVodeSetMaxNonlinIters", 1)) return(1);

	flag = CVodeSetMaxErrTestFails(cvode_mem, 100);
	if (check_flag(&flag, "CVodeSetMaxErrTestFails", 1)) return(1);

	if (iDense_Creact == 1) {
            printf("\n--> Using a Direct Dense Solver \n");

            /* Create dense SUNMatrix for use in linear solves */
	    A = SUNDenseMatrix(neq_tot, neq_tot);
	    if(check_flag((void *)A, "SUNDenseMatrix", 0)) return(1);

	    /* Create dense SUNLinearSolver object for use by CVode */
	    LS = SUNDenseLinearSolver(y, A);
	    if(check_flag((void *)LS, "SUNDenseLinearSolver", 0)) return(1);

	    /* Call CVDlsSetLinearSolver to attach the matrix and linear solver to CVode */
	    flag = CVDlsSetLinearSolver(cvode_mem, LS, A);
	    if(check_flag(&flag, "CVDlsSetLinearSolver", 1)) return(1);

#ifdef USE_KLU 
	} else if (iDense_Creact == 5) {

	    /* Create sparse SUNMatrix for use in linear solves */
	    A = SUNSparseMatrix(neq_tot, neq_tot, (data->NNZ)*NCELLS, CSC_MAT);
            if(check_flag((void *)A, "SUNSparseMatrix", 0)) return(1);

	    /* Create KLU solver object for use by CVode */
	    LS = SUNLinSol_KLU(y, A);
	    if(check_flag((void *)LS, "SUNLinSol_KLU", 0)) return(1);

	    /* Call CVodeSetLinearSolver to attach the matrix and linear solver to CVode */
	    flag = CVodeSetLinearSolver(cvode_mem, LS, A);
	    if(check_flag(&flag, "CVodeSetLinearSolver", 1)) return(1);
#endif

	} else if (iDense_Creact == 99) {
            printf("\n--> Using an Iterative Solver \n");

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
	} else {
	    int nJdata;
	    int HP = 0;
	    SPARSITY_INFO(&nJdata, &HP, 1);
            printf("--> SPARSE solver -- non zero entries %d represents %f %% sparsity pattern.", nJdata, nJdata/float((NEQ+1) * (NEQ+1)) *100.0);
	    amrex::Abort("\n \n");
	}

	if (iJac_Creact == 0) {
            printf("\n--> Without Analytical J\n");
#ifdef USE_KLU 
	    if (iDense_Creact == 5) {
	        amrex::Abort("\n--> SPARSE SOLVER SHOULD HAVE AN AJ...\n");
	    }
#endif
	} else {
            printf("\n--> With Analytical J\n");
	    if (iDense_Creact == 99) {
                if (iverbose > 0) {
                    printf("\n    (99)\n");
		}
	        /* Set the JAcobian-times-vector function */
	        flag = CVSpilsSetJacTimes(cvode_mem, NULL, NULL);
	        if(check_flag(&flag, "CVSpilsSetJacTimes", 1)) return(1);
	        /* Set the preconditioner solve and setup functions */
#ifdef USE_KLU 
	        flag = CVSpilsSetPreconditioner(cvode_mem, Precond_sparse, PSolve_sparse);
	        if(check_flag(&flag, "CVSpilsSetPreconditioner", 1)) return(1);
#else
	        flag = CVSpilsSetPreconditioner(cvode_mem, Precond, PSolve);
	        if(check_flag(&flag, "CVSpilsSetPreconditioner", 1)) return(1);
#endif
#ifdef USE_KLU 
	    } else if (iDense_Creact == 5){
		/* Set the user-supplied Jacobian routine Jac */
		flag = CVodeSetJacFn(cvode_mem, cJac_KLU);
		if(check_flag(&flag, "CVodeSetJacFn", 1)) return(1); 
#endif
	    } else {
                if (iverbose > 0) {
                    printf("\n    (1)\n");
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

	/* Define vectors to be used later in creact */
	if (iE_Creact == 1) { 
	    rhoe_init = (double *) malloc(NCELLS*sizeof(double));
	    rhoesrc_ext = (double *) malloc( NCELLS*sizeof(double));
	} else {
	    rhoh_init = (double *) malloc(NCELLS*sizeof(double));
	    rhohsrc_ext = (double *) malloc( NCELLS*sizeof(double));
	}
	rYsrc = (double *)  malloc((NCELLS*NEQ)*sizeof(double));

	N_VDestroy(atol);          /* Free the atol vector */

	/* Ok we're done ...*/
        if (iverbose > 0) {
	    printf(" --> DONE WITH INITIALIZATION (CPU) %d \n", iE_Creact);
	}
	reactor_cvode_initialized = true;

	return(0);
}

/* Main CVODE call routine */
int react(realtype *rY_in, realtype *rY_src_in, 
		realtype *rX_in, realtype *rX_src_in,
		realtype *P_in, 
                realtype *dt_react, realtype *time, int *Init){

        std::chrono::time_point<std::chrono::system_clock> start, end;		
        std::chrono::duration<double> total_elapsed;
	realtype time_init, time_out, dummy_time, temperature_save ;
	int flag;

        //FirstTimePrecond = true;

        time_init = *time;
	time_out  = *time + (*dt_react);
	start = std::chrono::system_clock::now();
        elapsed_seconds = start - start;
        elapsed_seconds_RHS = start - start;
        elapsed_seconds_Pcond = start - start;
	elapsed_seconds_JacFuego = start - start;

        if (iverbose > 3) {
	    printf("BEG : time curr is %14.6e and dt_react is %14.6e and final time should be %14.6e \n", time_init, *dt_react, time_out);
	}

	/* Get Device MemCpy of in arrays */
	/* Get Device pointer of solution vector */
	realtype *yvec_d      = N_VGetArrayPointer(y);
	// rhoY,T
	std::memcpy(yvec_d, rY_in, sizeof(realtype) * ((NEQ+1)*NCELLS));
	temperature_save = rY_in[NEQ];
	// rhoY_src_ext
	std::memcpy(rYsrc, rY_src_in, (NEQ*NCELLS)*sizeof(double));
	// rhoE/rhoH
	if (iE_Creact == 1) { 
	    std::memcpy(rhoe_init, rX_in, sizeof(realtype) * NCELLS);
	    std::memcpy(rhoesrc_ext, rX_src_in, sizeof(realtype) * NCELLS);
	} else {
	    std::memcpy(rhoh_init, rX_in, sizeof(realtype) * NCELLS);
	    std::memcpy(rhohsrc_ext, rX_src_in, sizeof(realtype) * NCELLS);
	}

	check_state(y);
	if (!actual_ok_to_react)  { 
	    amrex::Abort("\n Check_state failed: state is out of react bounds \n");
	}

	/* Call CVODE: ReInitPartial for faster start */
        if (iverbose > 1) {
            printf("\n -------------------------------------\n");
	}
	//if (*Init == 1) {
            if (iverbose > 1) {
                printf("ReInit always \n");
	    }
	    CVodeReInit(cvode_mem, time_init, y);
	    InitPartial = false;
	//} else {
	//    temp_old = fabs(rY_in[NEQ] - temp_old);
	//    // Sloppy but I can't think of anything better now
        //    if (temp_old > 50.0) {
        //        if (iverbose > 1) {
	//            printf("ReInit delta_T = %f \n", temp_old);
	//	}
	//        CVodeReInit(cvode_mem, time_init, y);
	//	InitPartial = false;
	//    } else {
        //        if (iverbose > 1) {
	//            printf("ReInit Partial delta_T = %f \n", temp_old);
	//	}
	//        CVodeReInitPartial(cvode_mem, time_init, y);
	//	InitPartial = true;
	//    }
	//}

	flag = CVode(cvode_mem, time_out, y, &dummy_time, CV_NORMAL);
	/* ONE STEP MODE FOR DEBUGGING */
	//flag = CVode(cvode_mem, time_out, y, &dummy_time, CV_ONE_STEP);
	if (check_flag(&flag, "CVode", 1)) return(1);

	*dt_react = dummy_time - time_init;
        if (iverbose > 3) {
	    printf("END : time curr is %14.6e and actual dt_react is %14.6e \n", dummy_time, *dt_react);
	}

	/* Pack data to return in main routine external */
	std::memcpy(rY_in, yvec_d, ((NEQ+1)*NCELLS)*sizeof(realtype));
	for  (int i = 0; i < NCELLS; i++) {
            rX_in[i] = rX_in[i] + (*dt_react) * rX_src_in[i];
	}

        end = std::chrono::system_clock::now();
        total_elapsed = end - start;

	if (*Init != 1) {
	    temp_old = rY_in[NEQ];
	}

	/* VERBOSE MODE */
        if (iverbose > 5) {
            for (int tid = 0; tid < NCELLS; tid ++) {
	        double rhov, energy, temp, energy2;
	        double MF[NEQ];
                //realtype activity[NEQ], cdot[NEQ], molecular_weight[NEQ];
	        int  lierr;
	        rhov = 0.0;
                int offset = tid * (NEQ + 1); 
                for (int k = 0; k < NEQ; k ++) {
	    	rhov =  rhov + rY_in[offset + k];
	        }
                //CKWT(molecular_weight);
                for (int k = 0; k < NEQ; k ++) {
	    	    MF[k] = rY_in[offset + k]/rhov;
	            //activity[k] = rY_in[offset + k]/(molecular_weight[k]);
	        }
	        energy = rX_in[tid]/rhov ;
	        if (iE_Creact == 1) { 
	            GET_T_GIVEN_EY(&energy, MF, &temp, &lierr);
	            CKHBMS(&temp, MF, &energy2);
	            CKUBMS(&temp, MF, &energy);
	    	    CKPY(&rhov, &temp, MF, P_in);
	            printf("e,h,p,rho ? %4.16e %4.16e %4.16e %4.16e \n",energy, energy2, *P_in, rhov);
	        } else {
	            GET_T_GIVEN_HY(&energy, MF, &temp, &lierr);
	            CKHBMS(&temp, MF, &energy);
	            CKUBMS(&temp, MF, &energy2);
	    	    CKPY(&rhov, &temp, MF, P_in);
	            printf("e,h,p,rho ? %4.16e %4.16e %4.16e %4.16e\n",energy2, energy, *P_in, rhov);
	        }
	        //rY_in[offset + NEQ] =  temp;
	        // DEBUG CHEKS
                //CKWC(&temp, activity, cdot);
	        // *P_in = cdot[2];
	    }

	    PrintFinalStats(cvode_mem, rY_in[NEQ]);

	} else if (iverbose > 2) {
	    printf("\nAdditional verbose info --\n");
	    PrintFinalStats(cvode_mem, temperature_save);
	    //std::cout << "Time spent computing Jac ? " << elapsed_seconds_JacFuego.count() << " " << elapsed_seconds_JacFuego.count() / total_elapsed.count() * 100.0 <<  std::endl; 
	    //std::cout << "Temp, chemistry solve, RHSeval, PSolve, Precond = " << rY_in[NEQ] << " " << total_elapsed.count() << " "<< elapsed_seconds_RHS.count() << " " << elapsed_seconds.count() << " " << elapsed_seconds_Pcond.count() << std::endl; 
	    //std::cout << "Temp, RHSeval represnts, PSolve represents, Precond represents = " << rY_in[NEQ] << " " << elapsed_seconds_RHS.count() / total_elapsed.count() * 100.0 <<  " " << elapsed_seconds.count()/total_elapsed.count() * 100.0 << " " << elapsed_seconds_Pcond.count() /total_elapsed.count() * 100.0 << std::endl;
            printf(" -------------------------------------\n");
	}

        long int nfe;
	flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
	return nfe;
}


/* RHS routine */
static int cF_RHS(realtype t, N_Vector y_in, N_Vector ydot_in, 
		void *user_data){

        std::chrono::time_point<std::chrono::system_clock> start, end;		
	start = std::chrono::system_clock::now();

	realtype *y_d      = N_VGetArrayPointer(y_in);
	realtype *ydot_d   = N_VGetArrayPointer(ydot_in);

        if (iE_Creact == 1) {
	    fKernelSpec(&t, y_d, ydot_d, 
			    rhoe_init, rhoesrc_ext, rYsrc);
	} else {
	    fKernelSpec(&t, y_d, ydot_d, 
			    rhoh_init, rhohsrc_ext, rYsrc);
	}
	end = std::chrono::system_clock::now();
        elapsed_seconds_RHS = elapsed_seconds_RHS + end - start;
	return(0);
}


/*
 * kernels
 */

/* RHS source terms evaluation */
void fKernelSpec(realtype *dt, realtype *yvec_d, realtype *ydot_d,  
		            double *rhoX_init, double *rhoXsrc_ext, double *rYs)
{
  int tid;
  for (tid = 0; tid < NCELLS; tid ++) {
      realtype massfrac[NEQ],activity[NEQ];
      realtype Xi[NEQ], cXi[NEQ];
      realtype cdot[NEQ], molecular_weight[NEQ];
      realtype temp, energy;
      int lierr;

      int offset = tid * (NEQ + 1); 

      /* MW CGS */
      CKWT(molecular_weight);
      /* rho MKS */ 
      realtype rho = 0.0;
      for (int i = 0; i < NEQ; i++){
          rho = rho + yvec_d[offset + i];
      }
      /* temp */
      temp = yvec_d[offset + NEQ];
      /* Yks, C CGS*/
      for (int i = 0; i < NEQ; i++){
          massfrac[i] = yvec_d[offset + i] / rho;
	  activity[i] = yvec_d[offset + i]/(molecular_weight[i]);
      }
      /* NRG CGS */
      energy = (rhoX_init[tid] + rhoXsrc_ext[tid]*(*dt)) /rho;

      /* Fuego calls on device */
      if (iE_Creact == 1){
          GET_T_GIVEN_EY(&energy, massfrac, &temp, &lierr);
          CKUMS(&temp, Xi);
          CKCVMS(&temp, cXi);
      } else {
          GET_T_GIVEN_HY(&energy, massfrac, &temp, &lierr);
          CKHMS(&temp, Xi);
          CKCPMS(&temp, cXi);
      }
      CKWC(&temp, activity, cdot);
      int cX = 0.0;
      for (int i = 0; i < NEQ; i++){
          cX = cX + massfrac[i] * cXi[i];
      }

      /* Fill ydot vect */
      ydot_d[offset + NEQ] = rhoXsrc_ext[tid];
      for (int i = 0; i < NEQ; i++){
          ydot_d[offset + i] = cdot[i] * molecular_weight[i] + rYs[tid * (NEQ) + i];
          ydot_d[offset + NEQ] = ydot_d[offset + NEQ]  - ydot_d[offset + i] * Xi[i];
      }
      ydot_d[offset + NEQ] = ydot_d[offset + NEQ] /(rho * cX);
  }
}


/* Analytical Jacobian evaluation */
static int cJac(realtype tn, N_Vector u, N_Vector fu, SUNMatrix J,
		void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){

  realtype *ydata  = N_VGetArrayPointer(u);

  int tid;
  for (tid = 0; tid < NCELLS; tid ++) {
      realtype *J_col_k;
      realtype  temp; 
      realtype activity[NEQ], molecular_weight[NEQ];
      realtype Jmat_tmp[(NEQ+1)*(NEQ+1)];

      int offset = tid * (NEQ + 1); 

      /* MW CGS */
      CKWT(molecular_weight);
      /* temp */
      temp = ydata[offset + NEQ];
      for (int i = 0; i < NEQ; i++){
          activity[i] = ydata[offset + i]/(molecular_weight[i]);
      }
      /* NRG CGS */
      int consP;
      if (iE_Creact == 1) {
	  consP = 0;
          DWDOT(Jmat_tmp, activity, &temp, &consP);
      } else {
          consP = 1;
          DWDOT(Jmat_tmp, activity, &temp, &consP);
      }
      /* fill the sunMat */
      for (int k = 0; k < NEQ; k++){
	  J_col_k = SM_COLUMN_D(J,offset + k);
	  for (int i = 0; i < NEQ; i++){
	        J_col_k[offset + i] = Jmat_tmp[k*(NEQ+1)+i] * molecular_weight[i] / molecular_weight[k]; 
          }
	  J_col_k[offset + NEQ] = Jmat_tmp[k*(NEQ+1)+NEQ] / molecular_weight[k]; 
      }
      J_col_k = SM_COLUMN_D(J,offset + NEQ);
      for (int i = 0; i < NEQ; i++){
          J_col_k[offset + i] = Jmat_tmp[NEQ*(NEQ+1)+i] * molecular_weight[i]; 
      }
  }

  return(0);

}


#ifdef USE_KLU 
/* Analytical SPARSE Jacobian evaluation */
static int cJac_KLU(realtype tn, N_Vector u, N_Vector fu, SUNMatrix J,
		void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){

  /* Make local copies of pointers to input data (big M) */
  realtype *ydata  = N_VGetArrayPointer(u);
  sunindextype *colptrs_tmp = SUNSparseMatrix_IndexPointers(J);
  sunindextype *rowvals_tmp = SUNSparseMatrix_IndexValues(J);
  realtype *data = SUNSparseMatrix_Data(J);

  /* Make local copies of pointers in user_data (cell M)*/
  UserData data_wk;
  data_wk = (UserData) user_data;   

  /* MW CGS */
  realtype molecular_weight[NEQ];
  CKWT(molecular_weight);

  /* Fixed RowVals */
  for (int i=0;i<data_wk->NNZ;i++) {
      rowvals_tmp[i] = data_wk->rowVals[0][i];
  }
  /* Fixed colPtrs */
  colptrs_tmp[0] = data_wk->colPtrs[0][0];
  for (int i=0;i<NCELLS*(NEQ + 1);i++) {
      colptrs_tmp[i+1] = data_wk->colPtrs[0][i+1];
  }

  /* Temp vectors */
  realtype temp_save_lcl, temp;
  realtype activity[NEQ];
  realtype Jmat_tmp[(NEQ+1)*(NEQ+1)];
  int tid, offset, nbVals, idx;
  temp_save_lcl = 0.0;
  for (tid = 0; tid < NCELLS; tid ++) {
      offset = tid * (NEQ + 1); 
      /* temp */
      temp = ydata[offset + NEQ];
      /* Do we recompute Jac ? */
      if (fabs(temp - temp_save_lcl) > 1.0) {
          for (int i = 0; i < NEQ; i++){
              activity[i] = ydata[offset + i]/(molecular_weight[i]);
          }
          /* NRG CGS */
          int consP;
          if (iE_Creact == 1) {
              consP = 0;
              DWDOT(Jmat_tmp, activity, &temp, &consP);
          } else {
              consP = 1;
              DWDOT(Jmat_tmp, activity, &temp, &consP);
          }
	  temp_save_lcl = temp;
	  /* rescale */
          for (int i = 0; i < NEQ; i++) {
              for (int k = 0; k < NEQ; k++) {
                  Jmat_tmp[k*(NEQ+1) + i] = Jmat_tmp[k*(NEQ+1) + i] * molecular_weight[i] / molecular_weight[k];
              }
              Jmat_tmp[i*(NEQ+1) + NEQ] = Jmat_tmp[i*(NEQ+1) + NEQ] / molecular_weight[i];
          }
          for (int i = 0; i < NEQ; i++) {
              Jmat_tmp[NEQ*(NEQ+1) + i] = Jmat_tmp[NEQ*(NEQ+1) + i] * molecular_weight[i];
          }
      }
      /* Go from Dense to Sparse */
      //printf("\n (new cell) \n");
      for (int i = 1; i < NEQ+2; i++) {
	  nbVals = data_wk->colPtrs[0][i]-data_wk->colPtrs[0][i - 1];
	  //printf(" Col %d has %d vals \n", i-1, nbVals);
	  for (int j = 0; j < nbVals; j++) {
	          idx = data_wk->rowVals[0][ data_wk->colPtrs[0][i - 1] + j ];
	              data[ data_wk->colPtrs[0][offset + i - 1] + j ] = Jmat_tmp[(i - 1) * (NEQ + 1) + idx];
	  }
      }
  }

  return(0);

}
#endif


/* Jacobian-times-vector routine.
 * Currently not used !!
 */
static int jtv(N_Vector v, N_Vector Jv, realtype t, N_Vector u, N_Vector fu, 
			  void *user_data, N_Vector tmp)
{
  realtype *udata, *vdata, *Jvdata;
  vdata  = N_VGetArrayPointer(v); 
  udata  = N_VGetArrayPointer(u);
  Jvdata = N_VGetArrayPointer(Jv);

  int tid;
  for (tid = 0; tid < NCELLS; tid ++) {
	realtype temp;
	realtype activity[NEQ], molecular_weight[NEQ];
	realtype J[(NEQ+1)*(NEQ+1)];

        int offset = tid * (NEQ + 1); 

        /* MW CGS */
	CKWT(molecular_weight);
	for (int i = 0; i < NEQ; i++){
            activity[i] = udata[offset + i]/(molecular_weight[i]);
	}
        /* temp */
	temp = udata[offset + NEQ];
        /* NRG CGS */
        if (iE_Creact == 1) {
            int consP = 0;
            DWDOT(J, activity, &temp, &consP);
        } else {
	    int consP = 1 ;
            DWDOT(J, activity, &temp, &consP);
        }

	/* PRINT JAC INFO: debug mode
	for (int i = 0; i < NEQ+1; i++){
	    for (int j = 0; j < NEQ+1; j++){
                printf(" %4.8e ", (J[j*(NEQ+1)+i]));
	    }
            printf(" \n");
	}
	amrex::Abort("\n--> ABORT\n");
	*/

	/* reinit */
	for (int i = 0; i < NEQ+1; i++){
	    Jvdata[offset + i] = 0.0;
	}
	/* COmpute */
	for (int i = 0; i < NEQ; i++){
            for (int j = 0; j < NEQ; j++){
                Jvdata[offset + i] = Jvdata[offset + i] + J[j*(NEQ+1)+i] * vdata[offset + j] * molecular_weight[i] / molecular_weight[j];
	    }
            Jvdata[offset + i] = Jvdata[offset + i] + J[NEQ*(NEQ+1)+i] * vdata[offset + NEQ] * molecular_weight[i];
        }
	for (int j = 0; j < NEQ; j++){
	    Jvdata[offset + NEQ] = Jvdata[offset + NEQ] + J[j*(NEQ+1)+NEQ] * vdata[offset + j] / molecular_weight[j];
	}
	Jvdata[offset + NEQ] = Jvdata[offset + NEQ]  + J[(NEQ+1)*(NEQ+1)-1] * vdata[offset + NEQ] ;  
  }

  return(0);
}


#ifdef USE_KLU 
/* Preconditioner setup routine for GMRES solver when KLU sparse mode is activated 
 * Generate and preprocess P
*/
static int Precond_sparse(realtype tn, N_Vector u, N_Vector fu, booleantype jok, 
		booleantype *jcurPtr, realtype gamma, void *user_data)
{
  //std::chrono::time_point<std::chrono::system_clock> start, end;		
  //std::chrono::time_point<std::chrono::system_clock> start_Jcomp, end_Jcomp;		
  //std::chrono::time_point<std::chrono::system_clock> start_JNcomp, end_JNcomp;		
  //std::chrono::time_point<std::chrono::system_clock> start_LUfac, end_LUfac;		
  
  //std::chrono::duration<double> elapsed_seconds_Jcomp;
  //std::chrono::duration<double> elapsed_seconds_JNcomp;
  //std::chrono::duration<double> elapsed_seconds_LUfac;
  //std::chrono::duration<double> elapsed_seconds_Pcond_prov;

  //start = std::chrono::system_clock::now();

  int ok,tid;

  /* MW CGS */
  realtype molecular_weight[NEQ];
  CKWT(molecular_weight);

  /* Formalism */
  int consP;
  if (iE_Creact == 1) { 
       consP = 0;
  } else {
       consP = 1;
  }

  /* Make local copies of pointers in user_data, and of pointer to u's data */
  UserData data_wk;
  data_wk = (UserData) user_data;   
  realtype **(**Jbd);
  Jbd = (data_wk->Jbd);
  realtype *udata;
  udata = N_VGetArrayPointer(u);

  //start_Jcomp = std::chrono::system_clock::now();
  //printf("temp_save %f \n",temp_save);
  if (jok) {
        /* jok = SUNTRUE: Copy Jbd to P */
        *jcurPtr = SUNFALSE;
  } else {
        /* Temp vectors */
        realtype temp, temp_save_lcl;
        realtype Jmat[(NEQ+1)*(NEQ+1)];
        realtype activity[NEQ];
	int offset,nbVals,idx;
        temp_save_lcl = 0.0;
        for (tid = 0; tid < NCELLS; tid ++) {
            offset = tid * (NEQ + 1); 
            /* temp */
            temp = udata[offset + NEQ];
            /* Do we recompute Jac ? */
            if (fabs(temp - temp_save_lcl) > 1.0) {
                for (int i = 0; i < NEQ; i++){
		    activity[i] = udata[offset + i]/(molecular_weight[i]);
                }
                DWDOT_PRECOND(Jmat, activity, &temp, &consP);

                /* Compute Jacobian.  Load into P. */
                denseScale(0.0, Jbd[tid][tid], NEQ+1, NEQ+1);

                for (int i = 0; i < NEQ; i++) {
                    for (int k = 0; k < NEQ; k++) {
                        (Jbd[tid][tid])[k][i] = Jmat[k*(NEQ+1) + i] * molecular_weight[i] / molecular_weight[k];
                    }
                    (Jbd[tid][tid])[i][NEQ] = Jmat[i*(NEQ+1) + NEQ] / molecular_weight[i];
                }
                for (int i = 0; i < NEQ; i++) {
                    (Jbd[tid][tid])[NEQ][i] = Jmat[NEQ*(NEQ+1) + i] * molecular_weight[i];
                }
	        (Jbd[tid][tid])[NEQ][NEQ] = Jmat[(NEQ+1)*(NEQ+1)-1];
	        temp_save_lcl = temp;
	    } else {
		/* if not: copy the one from prev cell */
		for (int i = 0; i < NEQ+1; i++) {
		    for (int k = 0; k < NEQ+1; k++) {
		        (Jbd[tid][tid])[i][k] = (Jbd[tid-1][tid-1])[i][k];
		    }
		}
	    }
	}

        *jcurPtr = SUNTRUE;
  }
  //end_Jcomp = std::chrono::system_clock::now();
  //elapsed_seconds_Jcomp = end_Jcomp - start_Jcomp;

  //start_JNcomp = std::chrono::system_clock::now();
  int nbVals;
  for (int i = 1; i < NEQ+2; i++) {
      /* nb non zeros elem should be the same for all cells */
      nbVals = data_wk->colPtrs[0][i]-data_wk->colPtrs[0][i-1];
      //printf("Rows %d : we have %d nonzero values \n", i-1, nbVals);
      for (int j = 0; j < nbVals; j++) {
    	  /* row of non zero elem should be the same for all cells */
    	  int idx = data_wk->rowVals[0][ data_wk->colPtrs[0][i-1] + j ];
          /* Scale by -gamma */
          /* Add identity matrix */
          for (tid = 0; tid < NCELLS; tid ++) {
    	      if (idx == (i-1)) {
                  data_wk->Jdata[tid][ data_wk->colPtrs[tid][i-1] + j ] = 1.0 - gamma * (Jbd[tid][tid])[ i-1 ][ idx ]; 
    	      } else {
                  data_wk->Jdata[tid][ data_wk->colPtrs[tid][i-1] + j ] = - gamma * (Jbd[tid][tid])[ i-1 ][ idx ]; 
    	      }
          }
      }
  }
  
  //end_JNcomp = std::chrono::system_clock::now();
  //elapsed_seconds_JNcomp = end_JNcomp - start_JNcomp;

  //start_LUfac = std::chrono::system_clock::now();
  if (!FirstTimePrecond) {
      //printf("and reuse pivots ...");
      for (tid = 0; tid < NCELLS; tid ++) {
          ok = klu_refactor(data_wk->colPtrs[tid], data_wk->rowVals[tid], data_wk->Jdata[tid], data_wk->Symbolic[tid], data_wk->Numeric[tid], &(data_wk->Common[tid]));
      }
  } else {
      //printf("and compute pivots ...");
      for (tid = 0; tid < NCELLS; tid ++) {
          data_wk->Numeric[tid] = klu_factor(data_wk->colPtrs[tid], data_wk->rowVals[tid], data_wk->Jdata[tid], data_wk->Symbolic[tid], &(data_wk->Common[tid])) ; 
      }
      FirstTimePrecond = false;
  }
  //end_LUfac = std::chrono::system_clock::now();
  //elapsed_seconds_LUfac = end_LUfac - start_LUfac;

  //end = std::chrono::system_clock::now();
  //elapsed_seconds_Pcond = elapsed_seconds_Pcond + end - start;

  //elapsed_seconds_Pcond_prov = end - start;
  //std::cout << " stats (Jcomp,JNcomp,pivots,TOTAL) = " << elapsed_seconds_Jcomp.count() <<  " " << elapsed_seconds_JNcomp.count() << " " << elapsed_seconds_LUfac.count() << " " << elapsed_seconds_Pcond_prov.count() << std::endl;

  return(0);
}
#endif


/* Preconditioner setup routine for GMRES solver when no sparse mode is activated 
 * Generate and preprocess P
*/
static int Precond(realtype tn, N_Vector u, N_Vector fu, booleantype jok, 
		booleantype *jcurPtr, realtype gamma, void *user_data)
{
  std::chrono::time_point<std::chrono::system_clock> start, end;		
  start = std::chrono::system_clock::now();

  UserData data_wk;
  realtype **(**P), **(**Jbd);
  sunindextype *(**pivot), ierr;
  realtype *udata; //, **a, **j;
  realtype temp;
  realtype Jmat[(NEQ+1)*(NEQ+1)];
  realtype activity[NEQ];

  /* MW CGS */
  realtype molecular_weight[NEQ];
  CKWT(molecular_weight);

  /* Make local copies of pointers in user_data, and of pointer to u's data */
  data_wk = (UserData) user_data;   
  P = (data_wk->P);
  Jbd = (data_wk->Jbd);
  pivot = (data_wk->pivot);
  udata = N_VGetArrayPointer(u);

  if (jok) {
      /* jok = SUNTRUE: Copy Jbd to P */
      denseCopy(Jbd[0][0], P[0][0], NEQ+1, NEQ+1);
      *jcurPtr = SUNFALSE;
  } else {
      /* jok = SUNFALSE: Generate Jbd from scratch and copy to P */
      /* Make local copies of problem variables, for efficiency. */
      for (int i = 0; i < NEQ; i++){
          activity[i] = udata[i]/(molecular_weight[i]);
      }
      temp = udata[NEQ];

      // C in mol/cm3
      int consP;
      if (iE_Creact == 1) { 
          consP = 0;
      } else {
          consP = 1;
      }
      DWDOT_PRECOND(Jmat, activity, &temp, &consP);
      //DWDOT(Jmat, activity, &temp, &consP);
      /* Compute Jacobian.  Load into P. */
      denseScale(0.0, Jbd[0][0], NEQ+1, NEQ+1);
      for (int i = 0; i < NEQ; i++) {
          for (int k = 0; k < NEQ; k++) {
              (Jbd[0][0])[k][i] = Jmat[k*(NEQ+1) + i] * molecular_weight[i] / molecular_weight[k];
          }
          (Jbd[0][0])[i][NEQ] = Jmat[i*(NEQ+1) + NEQ] / molecular_weight[i];
      }
      for (int i = 0; i < NEQ; i++) {
          (Jbd[0][0])[NEQ][i] = Jmat[NEQ*(NEQ+1) + i] * molecular_weight[i];
      }
      (Jbd[0][0])[NEQ][NEQ] = Jmat[(NEQ+1)*(NEQ+1)-1];

      denseCopy(Jbd[0][0], P[0][0], NEQ+1, NEQ+1);

      *jcurPtr = SUNTRUE;
  }
  
  /* Scale by -gamma */
  denseScale(-gamma, P[0][0], NEQ+1, NEQ+1);
  //denseScale(0.0, P[0][0], NEQ+1, NEQ+1);
  
  /* Add identity matrix and do LU decompositions on blocks in place. */
  denseAddIdentity(P[0][0], NEQ+1);
  ierr = denseGETRF(P[0][0], NEQ+1, NEQ+1, pivot[0][0]);
  if (ierr != 0) return(1);
  
  end = std::chrono::system_clock::now();
  elapsed_seconds_Pcond = elapsed_seconds_Pcond + end - start;

  return(0);
}


#ifdef USE_KLU 
/* PSolve for GMRES solver when KLU sparse mode is activated */
static int PSolve_sparse(realtype tn, N_Vector u, N_Vector fu, N_Vector r, N_Vector z,
                  realtype gamma, realtype delta, int lr, void *user_data)
{
  std::chrono::time_point<std::chrono::system_clock> start, end;		
  start = std::chrono::system_clock::now();

  UserData data_wk;
  data_wk = (UserData) user_data;

  realtype *zdata;
  zdata = N_VGetArrayPointer(z);

  N_VScale(1.0, r, z);

  /* Solve the block-diagonal system Pz = r using LU factors stored
     in P and pivot data in pivot, and return the solution in z. */
  int tid, offset_beg, offset_end;
  realtype zdata_cell[NEQ+1];
  for (tid = 0; tid < NCELLS; tid ++) {
      offset_beg = tid * (NEQ + 1); 
      offset_end = (tid + 1) * (NEQ + 1);
      std::memcpy(zdata_cell, zdata+offset_beg, (NEQ+1)*sizeof(realtype));
      klu_solve(data_wk->Symbolic[tid], data_wk->Numeric[tid], NEQ+1, 1, zdata_cell, &(data_wk->Common[tid])) ; 
      std::memcpy(zdata+offset_beg, zdata_cell, (NEQ+1)*sizeof(realtype));
  }

  end = std::chrono::system_clock::now();
  elapsed_seconds = elapsed_seconds + end - start;

  return(0);
}
#endif


/* PSolve for GMRES solver when no sparse mode is activated */
static int PSolve(realtype tn, N_Vector u, N_Vector fu, N_Vector r, N_Vector z,
                  realtype gamma, realtype delta, int lr, void *user_data)
{
  std::chrono::time_point<std::chrono::system_clock> start, end;		
  start = std::chrono::system_clock::now();

  realtype **(**P);
  sunindextype *(**pivot);
  realtype *zdata, *v;
  UserData data_wk;

  /* Extract the P and pivot arrays from user_data. */

  data_wk = (UserData) user_data;
  P = data_wk->P;
  pivot = data_wk->pivot;
  zdata = N_VGetArrayPointer(z);
  
  N_VScale(1.0, r, z);
  
  /* Solve the block-diagonal system Pz = r using LU factors stored
     in P and pivot data in pivot, and return the solution in z. */
  v = zdata;
  denseGETRS(P[0][0], NEQ+1, pivot[0][0], v);

  end = std::chrono::system_clock::now();
  elapsed_seconds = elapsed_seconds + end - start;
  //std::cout << " RHS duration " << elapsed_seconds.count() << std::endl;

  return(0);
}


/* 
 * OTHERS
*/


static void check_state(N_Vector yvec) 
{
  realtype *ydata;
  ydata = N_VGetArrayPointer(yvec);

  actual_ok_to_react = true;

  double rho, Temp;
  int offset;
  for (int tid = 0; tid < NCELLS; tid ++) {
      rho = 0.0;
      offset = tid * (NEQ + 1); 
      //amrex::Print() << " cell, offset, T " << tid << " " << offset << " " << ydata[offset + NEQ] << std::endl; 
      for (int k = 0; k < NEQ; k ++) {
          rho =  rho + ydata[offset + k];
      }
      Temp = ydata[offset + NEQ];
      if ((rho < 1.0e-10) || (rho > 1.e10)) {
          actual_ok_to_react = false;
      }
      if ((Temp < 200.0) || (Temp > 5000.0)) {
          actual_ok_to_react = false; 
      }
  }

}

/* Get and print some final statistics */
static void PrintFinalStats(void *cvodeMem, realtype Temp)
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

  if (iDense_Creact == 1){
      flag = CVDlsGetNumRhsEvals(cvodeMem, &nfeLS);
      check_flag(&flag, "CVDlsGetNumRhsEvals", 1);
      flag = CVDlsGetNumJacEvals(cvodeMem, &nje);
      check_flag(&flag, "CVDlsGetNumJacEvals", 1);
  } else if (iDense_Creact == 99){
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

  printf("-- Final Statistics --\n");
  printf("NonLinear (Newton) related --\n");
  if (InitPartial) {
         printf("    DT(dtL, dtI), RHS, Iterations, ConvFails, LinSolvSetups = %f %-6ld(%14.6e, %14.6e) %-6ld %-6ld %-6ld %-6ld \n",
	 Temp, nst-nst_old, hlast, hinused, nfe-nfe_old, nni-nni_old, ncfn-ncfn_old, nsetups-nsetups_old);
  }else{
         printf("    DT(dt, dtcur), RHS, Iterations, ErrTestFails, LinSolvSetups = %f %-6ld(%14.6e %14.6e) %-6ld %-6ld %-6ld %-6ld \n",
	 Temp, nst, hlast, hcur, nfe, nni, netfails, nsetups);
  }
  /* RESET */
  nsetups_old = nsetups;
  nst_old = nst;
  nfe_old = nfe;
  nni_old = nni;
  ncfn_old = ncfn;

  if (iDense_Creact == 1){
      printf("Linear (Dense Direct Solve) related --\n");
      if (InitPartial) {
          printf("    FD RHS, NumJacEvals                           = %f %-6ld %-6ld \n", Temp, nfeLS-nfeLS_old, nje-nje_old);
      } else {
          printf("    FD RHS, NumJacEvals                           = %f %-6ld %-6ld \n", Temp, nfeLS, nje);
      }
      /* RESET */
      nfeLS_old = nfeLS;
      nje_old = nje;
  } else if (iDense_Creact == 99){
	  // LinSolvSetups actually reflects the number of time the LinSolver has been called. 
	  // NonLinIterations can be taken without the need for LinItes
      printf("Linear (Krylov GMRES Solve) related --\n");
      if (InitPartial) {
          printf("    RHSeval, jtvEval, NumPrecEvals, NumPrecSolves = %f %-6ld %-6ld %-6ld %-6ld \n", 
            	      Temp, nfeLS-nfeLS_old, nje-nje_old, npe-npe_old, nps-nps_old);
      } else {
          printf("    RHSeval, jtvEval, NumPrecEvals, NumPrecSolves = %f %-6ld %-6ld %-6ld %-6ld \n", 
            	      Temp, nfeLS, nje, npe, nps);
          printf("    Iterations, ConvFails = %f %-6ld %-6ld \n", 
            	      Temp, nli, ncfl );
      }
      /* RESET */
      nfeLS_old = nfeLS;
      nje_old = nje;
      npe_old = npe;
      nps_old = nps;
  }
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


/* Alloc Data for CVODE */
static UserData AllocUserData(void)
{
  UserData data_wk;

  data_wk = (UserData) malloc(sizeof *data_wk);

  if (iDense_Creact == 99) {
      /* Precond data */
      printf("Alloc stuff for Precond \n");
      (data_wk->P) = new realtype***[NCELLS];
      (data_wk->Jbd) = new realtype***[NCELLS];
      (data_wk->pivot) = new sunindextype**[NCELLS];
      for(int i = 0; i < NCELLS; ++i) {
              (data_wk->P)[i] = new realtype**[NCELLS];
              (data_wk->Jbd)[i] = new realtype**[NCELLS];
              (data_wk->pivot)[i] = new sunindextype*[NCELLS];
      }

      for(int i = 0; i < NCELLS; ++i) {
          (data_wk->P)[i][i] = newDenseMat(NEQ+1, NEQ+1);
          (data_wk->Jbd)[i][i] = newDenseMat(NEQ+1, NEQ+1);
          (data_wk->pivot)[i][i] = newIndexArray(NEQ+1);
      }
  } 

#ifdef USE_KLU 
  /* Sparse Direct and Sparse (It) Precond data */
  data_wk->colPtrs = new int*[NCELLS];
  data_wk->rowVals = new int*[NCELLS];
  data_wk->Jdata = new realtype*[NCELLS];

  int HP;
  if (iE_Creact == 1) {
      HP = 0;
  } else {
      HP = 1;
  }
  if (iDense_Creact == 5) {
      /* Sparse Matrix for Direct Sparse KLU solver */
      (data_wk->PS) = new SUNMatrix[1];
      SPARSITY_INFO(&(data_wk->NNZ),&HP,NCELLS);
      printf("--> SPARSE solver -- non zero entries %d represents %f %% fill pattern.\n", data_wk->NNZ, data_wk->NNZ/float((NEQ+1) * (NEQ+1) * NCELLS * NCELLS) *100.0);
      //for(int i = 0; i < NCELLS; ++i) {
          (data_wk->PS)[0] = SUNSparseMatrix((NEQ+1)*NCELLS, (NEQ+1)*NCELLS, data_wk->NNZ, CSC_MAT);
          data_wk->colPtrs[0] = (int*) SUNSparseMatrix_IndexPointers((data_wk->PS)[0]); 
          data_wk->rowVals[0] = (int*) SUNSparseMatrix_IndexValues((data_wk->PS)[0]);
          data_wk->Jdata[0] = SUNSparseMatrix_Data((data_wk->PS)[0]);
          SPARSITY_PREPROC(data_wk->rowVals[0],data_wk->colPtrs[0],&HP,NCELLS);
      //}

  } else if (iDense_Creact == 99) {
      /* KLU internal storage */
      data_wk->Common = new klu_common[NCELLS];
      data_wk->Symbolic = new klu_symbolic*[NCELLS];
      data_wk->Numeric = new klu_numeric*[NCELLS];
      /* Sparse Matrices for It Sparse KLU block-solve */
      data_wk->PS = new SUNMatrix[NCELLS];
      /* Nb of non zero elements*/
      SPARSITY_INFO_PRECOND(&(data_wk->NNZ),&HP);
      printf("--> SPARSE Preconditioner -- non zero entries %d represents %f %% fill pattern.\n", data_wk->NNZ, data_wk->NNZ/float((NEQ+1) * (NEQ+1)) *100.0);
      for(int i = 0; i < NCELLS; ++i) {
          (data_wk->PS)[i] = SUNSparseMatrix(NEQ+1, NEQ+1, data_wk->NNZ, CSC_MAT);
          data_wk->colPtrs[i] = (int*) SUNSparseMatrix_IndexPointers((data_wk->PS)[i]); 
          data_wk->rowVals[i] = (int*) SUNSparseMatrix_IndexValues((data_wk->PS)[i]);
          data_wk->Jdata[i] = SUNSparseMatrix_Data((data_wk->PS)[i]);
          SPARSITY_PREPROC_PRECOND(data_wk->rowVals[i],data_wk->colPtrs[i],&HP);
          klu_defaults (&(data_wk->Common[i]));
          //data_wk->Common.btf = 0;
          //(data_wk->Common[i]).maxwork = 15;
          //data_wk->Common.ordering = 1;
          data_wk->Symbolic[i] = klu_analyze (NEQ+1, data_wk->colPtrs[i], data_wk->rowVals[i], &(data_wk->Common[i])) ; 
      }
  }
#endif

  return(data_wk);
}


/* Free memory */
void reactor_close(){

  CVodeFree(&cvode_mem);
  SUNLinSolFree(LS);
  if (iDense_Creact == 1) {
    SUNMatDestroy(A);
  }
  N_VDestroy(y); 
  FreeUserData(data);

  if (iE_Creact == 1) { 
    free(rhoe_init);
    free(rhoesrc_ext);
  } else {
    free(rhoh_init);
    free(rhohsrc_ext);
  }
  free(rYsrc);

  reactor_cvode_initialized = false;
}


/* Free data memory 
 * Probably not complete, how about the stuff allocated in KLU mode ? */
static void FreeUserData(UserData data_wk)
{
  if (iDense_Creact == 99) {
      for(int i = 0; i < NCELLS; ++i) {
          destroyMat((data_wk->P)[i][i]);
          destroyMat((data_wk->Jbd)[i][i]);
          destroyArray((data_wk->pivot)[i][i]);
      }
  }
  free(data_wk);
} 









