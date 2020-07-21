#include <reactor.h> 
#include <AMReX_ParmParse.H>
#include <chemistry_file.H>
#include "mechanism.h"
#include <AMREX_misc.H>

#define SUN_CUSP_CONTENT(S)        ( (SUNLinearSolverContent_Sparse_custom)(S->content) )
#define SUN_CUSP_REACTYPE(S)       ( SUN_CUSP_CONTENT(S)->reactor_type )
#define SUN_CUSP_NUM_SUBSYS(S)     ( SUN_CUSP_CONTENT(S)->nsubsys )
#define SUN_CUSP_SUBSYS_NNZ(S)     ( SUN_CUSP_CONTENT(S)->subsys_nnz )
#define SUN_CUSP_SUBSYS_SIZE(S)     ( SUN_CUSP_CONTENT(S)->subsys_size )

using namespace amrex;

/**********************************/
/* Global Variables */
  N_Vector         y = NULL; 
  SUNLinearSolver LS = NULL;
  SUNMatrix A        = NULL;
  void *cvode_mem    = NULL;
  /* User data */
  UserData data      = NULL;
/* OPTIONS */
  Real time_init    = 0.0;
  Array<double,NUM_SPECIES+1> typVals = {-1};
  double relTol       = 1.0e-10;
  double absTol       = 1.0e-10;
/* REMOVE MAYBE LATER */
  int dense_solve           = 1;
  int sparse_solve          = 5;
  int iterative_gmres_solve = 99;
  int sparse_solve_custom   = 101;
  int iterative_gmres_solve_custom = 199;
  int hack_dump_sparsity_pattern = -5;
  int eint_rho = 1; // in/out = rhoE/rhoY
  int enth_rho = 2; // in/out = rhoH/rhoY 

#ifdef _OPENMP
#pragma omp threadprivate(y,LS,A)
#pragma omp threadprivate(cvode_mem,data)
#pragma omp threadprivate(time_init)
#pragma omp threadprivate(typVals)
#pragma omp threadprivate(relTol,absTol)
#endif
/**********************************/

/**********************************/
/* Set or update typVals */
void SetTypValsODE(const std::vector<double>& ExtTypVals) {
    int size_ETV = (NUM_SPECIES + 1);
    Vector<std::string> kname;
    EOS::speciesNames(kname);
    int omp_thread = 0;
#ifdef _OPENMP
    omp_thread = omp_get_thread_num();
#endif

    for (int i=0; i<size_ETV-1; i++) {
      typVals[i] = ExtTypVals[i];
    }
    typVals[size_ETV-1] = ExtTypVals[size_ETV-1];
    if (omp_thread == 0){
        Print() << "Set the typVals in PelePhysics: \n  ";
        for (int i=0; i<size_ETV-1; i++) {
            Print() << kname[i] << ":" << typVals[i] << "  ";
        }
        Print() << "Temp:"<< typVals[size_ETV-1] <<  " \n";
    }
}


/* Set or update the rel/abs tolerances  */
void SetTolFactODE(double relative_tol,double absolute_tol) {
    relTol = relative_tol;
    absTol = absolute_tol;
    int omp_thread = 0;
#ifdef _OPENMP
    omp_thread = omp_get_thread_num();
#endif

    if (omp_thread == 0){
        Print() << "Set RTOL, ATOL = "<<relTol<< " "<<absTol<<  " in PelePhysics\n";
    }
}


/* Function to ReSet the tol of the cvode object directly */
void ReSetTolODE() {
    int omp_thread = 0;
#ifdef _OPENMP
    omp_thread = omp_get_thread_num();
#endif

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(data != NULL, "Reactor object is not initialized !!");

    int neq_tot     = (NUM_SPECIES + 1) * data->ncells;
    N_Vector atol   = N_VNew_Serial(neq_tot);
    realtype *ratol = N_VGetArrayPointer(atol);

    if (typVals[0] > 0) {
        if ((data->iverbose > 0) && (omp_thread == 0)) {
            Print() << "Setting CVODE tolerances rtol = " << relTol << " atolfact = " << absTol << " in PelePhysics \n";
        }
        for  (int i = 0; i < data->ncells; i++) {
            int offset = i * (NUM_SPECIES + 1);
            for  (int k = 0; k < NUM_SPECIES + 1; k++) {
                ratol[offset + k] = typVals[k]*absTol;
            }
        }
    } else {
        if ((data->iverbose > 0) && (omp_thread == 0)) {
            Print() << "Setting CVODE tolerances rtol = " << relTol << " atol = " << absTol << " in PelePhysics \n";
        }
        for (int i=0; i<neq_tot; i++) {
            ratol[i] = absTol;
        }
    }
    /* Call CVodeSVtolerances to specify the scalar relative tolerance
     * and vector absolute tolerances */
    int flag = CVodeSVtolerances(cvode_mem, relTol, atol);
    if (check_flag(&flag, "CVodeSVtolerances", 1)) { 
        Abort("Problem in ReSetTolODE");
    }

    N_VDestroy(atol);
}


/* Initialization routine, called once at the begining of the problem */
int reactor_init(int reactor_type, int ode_ncells) {

    BL_PROFILE_VAR("reactInit", reactInit);

    int omp_thread = 0;
#ifdef _OPENMP
    omp_thread = omp_get_thread_num();
#endif
    /* Total number of eq to integrate */
    int neq_tot = (NUM_SPECIES + 1) * ode_ncells;

    /* Definition of main vector */
    y = N_VNew_Serial(neq_tot);
    if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);

    /* Call CVodeCreate to create the solver memory and specify the
     * Backward Differentiation Formula and the use of a Newton iteration */
    cvode_mem = CVodeCreate(CV_BDF);
    if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);

    /* Does not work for more than 1 cell right now */
    data = AllocUserData(reactor_type, ode_ncells);
    if(check_flag((void *)data, "AllocUserData", 2)) return(1);

    /* Number of species and cells in mechanism */
    if ((data->iverbose > 0) && (omp_thread == 0)) {
        Print() << "Number of species in mech is " << NUM_SPECIES << "\n";
        Print() << "Number of cells in one solve is " << data->ncells << "\n";
    }

    /* Set the pointer to user-defined data */
    int flag = CVodeSetUserData(cvode_mem, data);
    if(check_flag(&flag, "CVodeSetUserData", 1)) return(1);   

    realtype time = 0.0e+0;
    /* Call CVodeInit to initialize the integrator memory and specify the
     * user's right hand side function, the inital time, and 
     * initial dependent variable vector y. */
    flag = CVodeInit(cvode_mem, cF_RHS, time, y);
    if (check_flag(&flag, "CVodeInit", 1)) return(1);
    
    /* Definition of tolerances: one for each species */
    N_Vector atol = N_VNew_Serial(neq_tot);
    realtype *ratol = N_VGetArrayPointer(atol);
    if (typVals[0] > 0) {
        if ((data->iverbose > 0) && (omp_thread == 0)) {
            Print() << "Setting CVODE tolerances rtol = " << relTol << " atolfact = " << absTol << " in PelePhysics \n";
        }
        for  (int i = 0; i < data->ncells; i++) {
            int offset = i * (NUM_SPECIES + 1);
            for  (int k = 0; k < NUM_SPECIES + 1; k++) {
                //ratol[offset + k] = std::max(typVals[k]*absTol,relTol);
                ratol[offset + k] = typVals[k]*absTol;
            }
        }
    } else {
        if ((data->iverbose > 0) && (omp_thread == 0)) {
            Print() << "Setting CVODE tolerances rtol = " << relTol << " atol = " << absTol << " in PelePhysics \n";
        }
        for (int i=0; i<neq_tot; i++) {
            ratol[i] = absTol;
        }
    }
    /* Call CVodeSVtolerances to specify the scalar relative tolerance
     * and vector absolute tolerances */
    flag = CVodeSVtolerances(cvode_mem, relTol, atol);
    if (check_flag(&flag, "CVodeSVtolerances", 1)) return(1);

    //flag = CVodeSetNonlinConvCoef(cvode_mem, 1.0e-1);
    //if (check_flag(&flag, "CVodeSetNonlinConvCoef", 1)) return(1);

    flag = CVodeSetMaxNonlinIters(cvode_mem, 50);
    if (check_flag(&flag, "CVodeSetMaxNonlinIters", 1)) return(1);

    flag = CVodeSetMaxErrTestFails(cvode_mem, 100);
    if (check_flag(&flag, "CVodeSetMaxErrTestFails", 1)) return(1);

    if (data->isolve_type == dense_solve) {
        if ((data->iverbose > 0) && (omp_thread == 0)) {
            Print() << "\n--> Using a Direct Dense Solver\n";
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

    } else if (data->isolve_type == sparse_solve_custom) {
        if ((data->iverbose > 0) && (omp_thread == 0)) {
            Print() << "\n--> Using a custom Direct Sparse Solver\n";
        }
        /* Create dense SUNMatrix for use in linear solves */
        A = SUNSparseMatrix(neq_tot, neq_tot, (data->NNZ)*data->ncells, CSR_MAT);
        if(check_flag((void *)A, "SUNDenseMatrix", 0)) return(1);

        /* Create dense SUNLinearSolver object for use by CVode */
        LS = SUNLinSol_sparse_custom(y, A, reactor_type, data->ncells, (NUM_SPECIES+1), data->NNZ);
        if(check_flag((void *)LS, "SUNDenseLinearSolver", 0)) return(1);

        /* Call CVDlsSetLinearSolver to attach the matrix and linear solver to CVode */
        flag = CVDlsSetLinearSolver(cvode_mem, LS, A);
        if(check_flag(&flag, "CVDlsSetLinearSolver", 1)) return(1);

    } else if (data->isolve_type == sparse_solve) {
#ifdef USE_KLU_PP 
        if ((data->iverbose > 0) && (omp_thread == 0)) {
            Print() << "\n--> Using a Direct Sparse Solver\n";
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
        Abort("Sparse solver not valid without KLU solver.");
#endif

    } else if ((data->isolve_type == iterative_gmres_solve) 
            || (data->isolve_type == iterative_gmres_solve_custom)) {
            if ((data->iverbose > 0) && (omp_thread == 0)) {
            Print() << "\n--> Using an Iterative Solver ("<<data->isolve_type<<")\n";
        }

            /* Create the linear solver object */
        if (data->ianalytical_jacobian == 0) {
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
        Abort("Wrong choice of linear solver...");
    }

    if (data->ianalytical_jacobian == 0) {
        if ((data->iverbose > 0) && (omp_thread == 0)) {
            Print() << "    Without Analytical J/Preconditioner\n";
        }
#ifdef USE_KLU_PP 
        if (data->isolve_type == sparse_solve) {
            Abort("Sparse Solver requires an Analytical J");
        }
#endif
        if (data->isolve_type == sparse_solve_custom) {
            Abort("Custom sparse solver requires an Analytical J");
        }
    } else {
        if (data->isolve_type == iterative_gmres_solve_custom) {
            /* Set the JAcobian-times-vector function */
            flag = CVSpilsSetJacTimes(cvode_mem, NULL, NULL);
            if(check_flag(&flag, "CVSpilsSetJacTimes", 1)) return(1);

            if ((data->iverbose > 0) && (omp_thread == 0)) {
                Print() << "    With a custom Sparse Preconditioner\n";
            }
            /* Set the preconditioner solve and setup functions */
            flag = CVSpilsSetPreconditioner(cvode_mem, Precond_custom, PSolve_custom);
            if(check_flag(&flag, "CVSpilsSetPreconditioner", 1)) return(1);

        } else if (data->isolve_type == iterative_gmres_solve) {
            /* Set the JAcobian-times-vector function */
            flag = CVSpilsSetJacTimes(cvode_mem, NULL, NULL);
            if(check_flag(&flag, "CVSpilsSetJacTimes", 1)) return(1);
#ifdef USE_KLU_PP 
            if ((data->iverbose > 0) && (omp_thread == 0)) {
                Print() << "    With a Sparse Preconditioner\n";
            }
            /* Set the preconditioner solve and setup functions */
            flag = CVSpilsSetPreconditioner(cvode_mem, Precond_sparse, PSolve_sparse);
            if(check_flag(&flag, "CVSpilsSetPreconditioner", 1)) return(1);
#else
            if ((data->iverbose > 0) && (omp_thread == 0)) {
                Print() << "    With a Preconditioner\n";
            }
            /* Set the preconditioner solve and setup functions */
            flag = CVSpilsSetPreconditioner(cvode_mem, Precond, PSolve);
            if(check_flag(&flag, "CVSpilsSetPreconditioner", 1)) return(1);
#endif
#ifdef USE_KLU_PP 
        } else if (data->isolve_type == sparse_solve){
            if ((data->iverbose > 0) && (omp_thread == 0)) {
                Print() << "    With a Sparse Analytical J\n";
            }
            /* Set the user-supplied Jacobian routine Jac */
            flag = CVodeSetJacFn(cvode_mem, cJac_KLU);
            if(check_flag(&flag, "CVodeSetJacFn", 1)) return(1); 
#endif
        } else if (data->isolve_type == dense_solve){
            if ((data->iverbose > 0) && (omp_thread == 0)) {
                Print() << "    With Analytical J\n";
            }
            /* Set the user-supplied Jacobian routine Jac */
            flag = CVodeSetJacFn(cvode_mem, cJac);
            if(check_flag(&flag, "CVodeSetJacFn", 1)) return(1);

        }  else if (data->isolve_type == sparse_solve_custom) {
            if ((data->iverbose > 0) && (omp_thread == 0)) {
                Print() << "    With a Sparse Analytical J\n";
            }
            /* Set the user-supplied Jacobian routine Jac */
            flag = CVodeSetJacFn(cvode_mem, cJac_sps);
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

    /* Free the atol vector */
    N_VDestroy(atol);

    /* Ok we're done ...*/
    if ((data->iverbose > 1) && (omp_thread == 0)) {
        Print() << "\n--> DONE WITH INITIALIZATION (CPU)" << data->ireactor_type << "\n";
    }

    /* Reactor is now initialized */
    data->reactor_cvode_initialized = true;

    BL_PROFILE_VAR_STOP(reactInit);

    return(0);
}

/* Main routine for CVode integration: integrate a Box version 1*/
int react(const Box& box,
          Array4<Real> const& rY_in,
          Array4<Real> const& rY_src_in,
          Array4<Real> const& T_in,
          Array4<Real> const& rEner_in,
          Array4<Real> const& rEner_src_in,
          Array4<Real> const& FC_in,
          Array4<int> const& mask,
          Real &dt_react,
          Real &time) {

    int omp_thread = 0;
#ifdef _OPENMP
    omp_thread = omp_get_thread_num(); 
#endif

    if ((data->iverbose > 1) && (omp_thread == 0)) {
        Print() <<"\n -------------------------------------\n";
    }

    /* Initial time and time to reach after integration */
    time_init = time;

    if ((data->iverbose > 3) && (omp_thread == 0)) {
        Print() <<"BEG : time curr is "<< time_init << " and dt_react is " << dt_react << " and final time should be " << time_init + dt_react << "\n";
    }

    if (data->ncells != 1) {
        Abort("CVODE react can only integrate one cell at a time");
    }
    int box_ncells  = box.numPts();
    data->boxcell   = 0; 

    if ((data->iverbose > 2) && (omp_thread == 0)) {
        Print() <<"Ncells in the box = "<<  box_ncells  << "\n";
    }

    BL_PROFILE_VAR("reactor::ExtForcingAlloc", ExtForcingAlloc);
    /* External forcing crap */
    if ((data->rhoX_init).size() != data->ncells) {
        (data->rhoX_init).resize(data->ncells);
        (data->rhoXsrc_ext).resize(data->ncells);
        (data->rYsrc).resize(data->ncells*NUM_SPECIES);
    }
    BL_PROFILE_VAR_STOP(ExtForcingAlloc);

    /* Perform integration one cell at a time */
    ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {

        Real mass_frac[NUM_SPECIES];
        Real rho = 0.0;
        Real rho_inv;
        Real Enrg_loc;
        Real temp;

        realtype *yvec_d      = N_VGetArrayPointer(y);
        
        BL_PROFILE_VAR("reactor::FlatStuff", FlatStuff);
        for (int n = 0; n < NUM_SPECIES; n++) {
            yvec_d[n]        = rY_in(i,j,k,n);
            (data->rYsrc)[n] = rY_src_in(i,j,k,n);
            rho += yvec_d[n]; 
        }
        rho_inv                 = 1.0 / rho;
        temp                    = T_in(i,j,k,0);
        (data->rhoX_init)[0]    = rEner_in(i,j,k,0); 
        (data->rhoXsrc_ext)[0]  = rEner_src_in(i,j,k,0);

        /* T update with energy and Y */
        for (int n = 0; n < NUM_SPECIES; n++) {
            mass_frac[n] = yvec_d[n] * rho_inv;
        }
        Enrg_loc = (data->rhoX_init)[0] / rho;
        if (data->ireactor_type == 1){
            EOS::EY2T(Enrg_loc,mass_frac,temp);
        } else {
            EOS::HY2T(Enrg_loc,mass_frac,temp);
        }
        yvec_d[NUM_SPECIES] = temp;
        BL_PROFILE_VAR_STOP(FlatStuff);

        /* ReInit CVODE is faster */
        CVodeReInit(cvode_mem, time_init, y);

        /* Time to reach after integration */
        Real time_out_lcl  = time_init + dt_react;

        /* Integration */
        Real dummy_time;
        BL_PROFILE_VAR("reactor::AroundCVODE", AroundCVODE);
        int flag = CVode(cvode_mem, time_out_lcl, y, &dummy_time, CV_NORMAL);
        //if (check_flag(&flag, "CVode", 1)) return(1);
        BL_PROFILE_VAR_STOP(AroundCVODE);

        if ((data->iverbose > 1) && (omp_thread == 0)) {
            Print() <<"Additional verbose info --\n";
            PrintFinalStats(cvode_mem, yvec_d[NUM_SPECIES]);
            Print() <<"\n -------------------------------------\n";
        }

        /* Get estimate of how hard the integration process was */
        long int nfe,nfeLS;
        flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
        flag = CVodeGetNumLinRhsEvals(cvode_mem, &nfeLS);
        FC_in(i,j,k,0) = nfe+nfeLS;

        BL_PROFILE_VAR_START(FlatStuff);
        rho = 0.0;
        for (int n = 0; n < NUM_SPECIES; n++) {
            rY_in(i,j,k,n) = yvec_d[n];
            rho += yvec_d[n]; 
        }
        rho_inv    = 1.0 / rho; 
        temp       = yvec_d[NUM_SPECIES];

        /* T update with energy and Y */
        for (int n = 0; n < NUM_SPECIES; n++) {
            mass_frac[n] = yvec_d[n] * rho_inv;
        }
        //Enrg_loc = ((data->rhoX_init)[0] + (dummy_time - time_init) * rEner_src_in(i,j,k,0)) / rho;
        Enrg_loc = ((data->rhoX_init)[0] + (dummy_time - time_init) * (data->rhoXsrc_ext)[0]) /rho;
        if (data->ireactor_type == 1){
            EOS::EY2T(Enrg_loc,mass_frac,temp);
        } else {
            EOS::HY2T(Enrg_loc,mass_frac,temp);
        }
        T_in(i,j,k,0) = temp;
        BL_PROFILE_VAR_STOP(FlatStuff);

        if ((data->iverbose > 3) && (omp_thread == 0)) {
            Print() <<"END : time curr is "<< dummy_time << " and actual dt_react is " << (dummy_time - time_init) << "\n";
        }
    });

    /* Update dt_react with real time step taken ... 
       should be very similar to input dt_react */
    //dt_react = dummy_time - time_init;
#ifdef MOD_REACTOR
    /* If reactor mode is activated, update time to perform subcycling */
    time  = time_init + dt_react;
#endif


    /* Get estimate of how hard the integration process was */
    return 20;
}


/* Main routine for CVode integration: integrate a Box version 2*/
int react_2(const Box& box,
          Array4<Real> const& rY_in,
          Array4<Real> const& rY_src_in,
          Array4<Real> const& T_in,
          Array4<Real> const& rEner_in,
          Array4<Real> const& rEner_src_in,
          Array4<Real> const& FC_in,
          Array4<int> const& mask,
          Real &dt_react,
          Real &time) {

    realtype dummy_time;
    int flag;
    int omp_thread = 0;
#ifdef _OPENMP
    omp_thread = omp_get_thread_num(); 
#endif

    if ((data->iverbose > 1) && (omp_thread == 0)) {
        Print() <<"\n -------------------------------------\n";
    }

    /* Initial time and time to reach after integration */
    time_init = time;
    realtype time_out  = time + dt_react;

    if ((data->iverbose > 3) && (omp_thread == 0)) {
        Print() <<"BEG : time curr is "<< time_init << " and dt_react is " << dt_react << " and final time should be " << time_out << "\n";
    }

    /* Define full box_ncells length vectors to be integrated piece by piece
       by CVode */
    int box_ncells  = box.numPts();
    if ((data->iverbose > 2) && (omp_thread == 0)) {
        Print() <<"Ncells in the box = "<<  box_ncells  << "\n";
    }
    BL_PROFILE_VAR("reactor::ExtForcingAlloc", ExtForcingAlloc);
    if ((data->rhoX_init).size() != box_ncells) {
        (data->Yvect_full).resize(box_ncells*(NUM_SPECIES+1));
        (data->rhoX_init).resize(box_ncells);
        (data->rhoXsrc_ext).resize(box_ncells);
        (data->rYsrc).resize(box_ncells*NUM_SPECIES);
        (data->FCunt).resize(box_ncells);
    }
    BL_PROFILE_VAR_STOP(ExtForcingAlloc);

    BL_PROFILE_VAR("reactor::FlatStuff", FlatStuff);
    /* Fill the full box_ncells length vectors from input Array4*/
    const auto len        = length(box);
    const auto lo         = lbound(box);
    ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        int icell = (k-lo.z)*len.x*len.y + (j-lo.y)*len.x + (i-lo.x);
        box_flatten(icell, i, j, k, data->ireactor_type,
                    rY_in, rY_src_in, T_in, 
                    rEner_in, rEner_src_in,
                    data->Yvect_full, data->rYsrc, data->rhoX_init, data->rhoXsrc_ext);
    });
    BL_PROFILE_VAR_STOP(FlatStuff);

    BL_PROFILE_VAR("reactor::AroundCVODE", AroundCVODE);
    BL_PROFILE_VAR_STOP(AroundCVODE);

    /* We may need extra cells to fill the fixed data->ncells in this case 
       since we do not Init each time */
    int extra_cells = box_ncells - box_ncells / (data->ncells) * (data->ncells);
    if ((data->iverbose > 2) && (omp_thread == 0)) {
        Print() <<" Extra cells = "<< extra_cells  << "\n";
    }
    
    /* Integrate data->ncells at a time with CVode 
       The extra cell machinery is not ope yet and most likely produce
       out of bound errors */
    realtype *yvec_d      = N_VGetArrayPointer(y);
    for  (int i = 0; i < box_ncells+extra_cells; i+=data->ncells) {
        //Print() <<" dealing with cell " << i <<  "\n";
        int offset = i * (NUM_SPECIES + 1);
        data->boxcell = i; 
        for  (int k = 0; k < data->ncells*(NUM_SPECIES+1); k++) {
            yvec_d[k] = data->Yvect_full[offset + k];
        }
        
        /* ReInit CVODE is faster */
        CVodeReInit(cvode_mem, time_init, y);

        BL_PROFILE_VAR_START(AroundCVODE);
        /* Integration */
        flag = CVode(cvode_mem, time_out, y, &dummy_time, CV_NORMAL);
        if (check_flag(&flag, "CVode", 1)) return(1);
        BL_PROFILE_VAR_STOP(AroundCVODE);

        /* Update full box length vector */
        for  (int k = 0; k < data->ncells*(NUM_SPECIES+1); k++) {
            data->Yvect_full[offset + k] = yvec_d[k];
        }

        /* Get estimate of how hard the integration process was */
        long int nfe,nfeLS;
        flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
        flag = CVodeGetNumLinRhsEvals(cvode_mem, &nfeLS);
        for  (int k = 0; k < data->ncells; k++) {
            data->FCunt[i + k] = nfe+nfeLS;
        }

        if ((data->iverbose > 3) && (omp_thread == 0)) {
            Print() <<"END : time curr is "<< dummy_time << " and actual dt_react is " << (dummy_time - time_init) << "\n";
        }

    }

#ifdef MOD_REACTOR
    /* If reactor mode is activated, update time to perform subcycling */
    time  = time_init + dt_react;
#endif

    BL_PROFILE_VAR_START(FlatStuff);
    /* Update the input/output Array4 rY_in and rEner_in*/
    ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        int icell = (k-lo.z)*len.x*len.y + (j-lo.y)*len.x + (i-lo.x);
        box_unflatten(icell, i, j, k, data->ireactor_type,
                    rY_in, T_in, rEner_in, rEner_src_in, FC_in,
                    data->Yvect_full, data->rhoX_init, data->FCunt, dt_react);
    });
    BL_PROFILE_VAR_STOP(FlatStuff);

    if ((data->iverbose > 1) && (omp_thread == 0)) {
        Print() <<"Additional verbose info --\n";
        PrintFinalStats(cvode_mem, yvec_d[NUM_SPECIES]);
        Print() <<"\n -------------------------------------\n";
    }

    /* Get estimate of how hard the integration process was */
    long int nfe,nfeLS;
    flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
    flag = CVodeGetNumLinRhsEvals(cvode_mem, &nfeLS);
    return nfe+nfeLS;
}


/* Main routine for CVode integration: classic version */
int react(realtype *rY_in, realtype *rY_src_in, 
          realtype *rX_in, realtype *rX_src_in,
          realtype &dt_react, realtype &time){

    realtype dummy_time;
    int flag;
    int omp_thread = 0;
#ifdef _OPENMP
    omp_thread = omp_get_thread_num();
#endif

    if ((data->iverbose > 1) && (omp_thread == 0)) {
        Print() <<"\n -------------------------------------\n";
    }

    /* Initial time and time to reach after integration */
    time_init = time;
    realtype time_out  = time + dt_react;

    if ((data->iverbose > 3) && (omp_thread == 0)) {
        Print() <<"BEG : time curr is "<< time_init << " and dt_react is " << dt_react << " and final time should be " << time_out << "\n";
    }

    /* Define full box_ncells length vectors to be integrated piece by piece
       by CVode */
    if ((data->iverbose > 2) && (omp_thread == 0)) {
        Print() <<"Ncells in the box = "<<  data->ncells  << "\n";
    }
    BL_PROFILE_VAR("reactor::ExtForcingAlloc", ExtForcingAlloc);
    if ((data->rhoX_init).size() != data->ncells) {
        (data->Yvect_full).resize(data->ncells*(NUM_SPECIES+1));
        (data->rYsrc).resize(data->ncells*NUM_SPECIES);
        (data->rhoX_init).resize(data->ncells);
        (data->rhoXsrc_ext).resize(data->ncells);
    }
    BL_PROFILE_VAR_STOP(ExtForcingAlloc);

    BL_PROFILE_VAR("reactor::FlatStuff", FlatStuff);
    /* Get Device MemCpy of in arrays */
    /* Get Device pointer of solution vector */
    realtype *yvec_d      = N_VGetArrayPointer(y);
    /* rhoY,T */
    std::memcpy(yvec_d,                     rY_in,     sizeof(Real) * ((NUM_SPECIES+1)*data->ncells));
    /* rhoY_src_ext */
    std::memcpy((data->rYsrc).data(),       rY_src_in, sizeof(Real) * (NUM_SPECIES * data->ncells));
    /* rhoE/rhoH */
    std::memcpy((data->rhoX_init).data(),   rX_in,     sizeof(Real) * data->ncells);
    std::memcpy((data->rhoXsrc_ext).data(), rX_src_in, sizeof(Real) * data->ncells);
    BL_PROFILE_VAR_STOP(FlatStuff);

    /* Check if y is within physical bounds
       we may remove that eventually */
    check_state(y);
    if (!(data->actual_ok_to_react))  { 
#ifdef MOD_REACTOR
        /* If reactor mode is activated, update time */
        time  = time_out;
#endif
        return 0;
    }

    BL_PROFILE_VAR_START(FlatStuff);
    /* T update with energy and Y */
    int offset;
    realtype rho, rho_inv, nrg_loc, temp;
    for  (int i = 0; i < data->ncells; i++) {
        offset = i * (NUM_SPECIES + 1);
        realtype* mass_frac = rY_in + offset;
        // get rho
        rho = 0;
        for  (int kk = 0; kk < NUM_SPECIES; kk++) {
            rho += mass_frac[kk];
        }
        rho_inv = 1 / rho; 
        // get Yks
        for  (int kk = 0; kk < NUM_SPECIES; kk++) {
            mass_frac[kk] = mass_frac[kk] * rho_inv;
        }
        // get energy
        nrg_loc = rX_in[i] * rho_inv;
        // recompute T
        if (data->ireactor_type == eint_rho){
            EOS::EY2T(nrg_loc,mass_frac,temp);
        } else {
            EOS::HY2T(nrg_loc,mass_frac,temp);
        }
        // store T in y
        yvec_d[offset + NUM_SPECIES] = temp;
    }
    BL_PROFILE_VAR_STOP(FlatStuff);

    /* ReInit CVODE is faster */
    CVodeReInit(cvode_mem, time_init, y);
    
    /* There should be no internal looping of CVOde */
    data->boxcell = 0;

    BL_PROFILE_VAR("reactor::AroundCVODE", AroundCVODE);
    flag = CVode(cvode_mem, time_out, y, &dummy_time, CV_NORMAL);
    /* ONE STEP MODE FOR DEBUGGING */
    //flag = CVode(cvode_mem, time_out, y, &dummy_time, CV_ONE_STEP);
    if (check_flag(&flag, "CVode", 1)) return(1);
    BL_PROFILE_VAR_STOP(AroundCVODE);

    /* Update dt_react with real time step taken ... 
       should be very similar to input dt_react */
    dt_react = dummy_time - time_init;
#ifdef MOD_REACTOR
    /* If reactor mode is activated, update time */
    time  = time_init + dt_react;
#endif

    if ((data->iverbose > 3) && (omp_thread == 0)) {
        Print() <<"END : time curr is "<< dummy_time << " and actual dt_react is " << dt_react << "\n";
    }

    BL_PROFILE_VAR_START(FlatStuff);
    /* Pack data to return in main routine external */
    std::memcpy(rY_in, yvec_d, ((NUM_SPECIES+1)*data->ncells)*sizeof(realtype));
    for  (int i = 0; i < data->ncells; i++) {
        rX_in[i] = rX_in[i] + dt_react * rX_src_in[i];
    }

    /* T update with energy and Y */
    for  (int i = 0; i < data->ncells; i++) {
        offset = i * (NUM_SPECIES + 1);
        realtype* mass_frac = yvec_d + offset;
        // get rho
        rho = 0;
        for  (int kk = 0; kk < NUM_SPECIES; kk++) {
            rho += mass_frac[kk];
        }
        rho_inv = 1 / rho;
        // get Yks
        for  (int kk = 0; kk < NUM_SPECIES; kk++) {
            mass_frac[kk] = mass_frac[kk] * rho_inv;
        }
        // get energy
        nrg_loc = rX_in[i] * rho_inv;
        // recompute T
        if (data->ireactor_type == eint_rho){
            EOS::EY2T(nrg_loc,mass_frac,temp);
        } else {
            EOS::HY2T(nrg_loc,mass_frac,temp);
        }
        // store T in rY_in
        rY_in[offset + NUM_SPECIES] = temp;
    }
    BL_PROFILE_VAR_STOP(FlatStuff);

    if ((data->iverbose > 1) && (omp_thread == 0)) {
        Print() <<"Additional verbose info --\n";
        PrintFinalStats(cvode_mem, rY_in[NUM_SPECIES]);
        Print() <<"\n -------------------------------------\n";
    }

    /* Get estimate of how hard the integration process was */
    long int nfe,nfeLS;
    flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
    flag = CVodeGetNumLinRhsEvals(cvode_mem, &nfeLS);
    return nfe+nfeLS;
}


/* RHS routine */
int cF_RHS(realtype t, N_Vector y_in, N_Vector ydot_in, 
           void *user_data){

    realtype *y_d      = N_VGetArrayPointer(y_in);
    realtype *ydot_d   = N_VGetArrayPointer(ydot_in);

    BL_PROFILE_VAR("fKernelSpec()", fKernelSpec);
    fKernelSpec(&t, y_d, ydot_d, user_data);
    BL_PROFILE_VAR_STOP(fKernelSpec);

    return(0);
}
/**********************************/



/*
 * kernels
 */

/* RHS source terms evaluation */
void fKernelSpec(realtype *t, realtype *yvec_d, realtype *ydot_d,  
                 void *user_data)
{
  /* Make local copies of pointers in user_data (cell M)*/
  UserData data_wk = (UserData) user_data;

  /* Loop on packed cells */
  for (int tid = 0; tid < data_wk->ncells; tid ++) {
      /* Tmp vars */
      realtype massfrac[NUM_SPECIES];
      realtype Xi[NUM_SPECIES];
      realtype cdot[NUM_SPECIES], molecular_weight[NUM_SPECIES];
      realtype cX;
      realtype temp, energy;
      realtype dt;

      /* dt is curr time - time init */
      dt = *t - time_init;

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
      energy = (data_wk->rhoX_init[data->boxcell + tid] + data_wk->rhoXsrc_ext[data_wk->boxcell + tid] * dt) /rho;

      if (data_wk->ireactor_type == eint_rho){
          /* UV REACTOR */
          EOS::EY2T(energy, massfrac, temp);
          EOS::TY2Cv(temp, massfrac, cX);
          EOS::T2Ei(temp, Xi);
      } else if (data_wk->ireactor_type == enth_rho) {
          /* HP REACTOR */
          EOS::HY2T(energy, massfrac, temp);
          EOS::TY2Cp(temp, massfrac, cX);
          EOS::T2Hi(temp, Xi);
      }
      EOS::RTY2WDOT(rho, temp, massfrac, cdot);

      /* Fill ydot vect */
      ydot_d[offset + NUM_SPECIES] = data_wk->rhoXsrc_ext[data_wk->boxcell + tid];
      for (int i = 0; i < NUM_SPECIES; i++){
          ydot_d[offset + i] = cdot[i] + data_wk->rYsrc[(data_wk->boxcell + tid) * (NUM_SPECIES) + i];
          ydot_d[offset + NUM_SPECIES] = ydot_d[offset + NUM_SPECIES]  - ydot_d[offset + i] * Xi[i];
      }
      ydot_d[offset + NUM_SPECIES] = ydot_d[offset + NUM_SPECIES] /(rho * cX);
  }
}


/*
 * Auxiliary routines
 */

/* Analytical Jacobian evaluation */
int cJac(realtype tn, N_Vector u, N_Vector fu, SUNMatrix J,
         void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){

  /* Make local copies of pointers to input data (big M) */
  realtype *ydata  = N_VGetArrayPointer(u);

  /* Make local copies of pointers in user_data (cell M)*/
  UserData data_wk = (UserData) user_data;

  BL_PROFILE_VAR("DenseJac", DenseJac);
  for (int tid = 0; tid < data_wk->ncells; tid ++) {
      /* Tmp vars */
      realtype *J_col_k;
      realtype massfrac[NUM_SPECIES], molecular_weight[NUM_SPECIES];
      realtype temp; 
      realtype Jmat_tmp[(NUM_SPECIES+1)*(NUM_SPECIES+1)];

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
      if (data_wk->ireactor_type == eint_rho) {
          consP = 0;
      } else {
          consP = 1;
      }
      EOS::RTY2JAC(rho, temp, massfrac, Jmat_tmp, consP);
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
      J_col_k = SM_COLUMN_D(J,offset);
  }
  BL_PROFILE_VAR_STOP(DenseJac);

  return(0);

}

/* Analytical SPARSE Jacobian evaluation */
int cJac_sps(realtype tn, N_Vector u, N_Vector fu, SUNMatrix J,
             void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){
  /* Make local copies of pointers to input data (big M) */
  realtype *ydata           = N_VGetArrayPointer(u);
  sunindextype *rowPtrs_tmp = SUNSparseMatrix_IndexPointers(J);
  sunindextype *colIndx_tmp = SUNSparseMatrix_IndexValues(J);
  realtype *Jdata           = SUNSparseMatrix_Data(J);

  /* Make local copies of pointers in user_data (cell M)*/
  UserData data_wk = (UserData) user_data;

  /* MW CGS */
  realtype molecular_weight[NUM_SPECIES];
  CKWT(molecular_weight);

  BL_PROFILE_VAR("FillSparseData", FillSpsData);
  /* Fixed colVals*/
  for (int i=0;i<data_wk->NNZ*data_wk->ncells;i++) {
      colIndx_tmp[i] = (sunindextype) data_wk->colVals_c[i];
  }
  rowPtrs_tmp[0] = (sunindextype) data_wk->rowPtrs_c[0];
  /* Fixed rowPtrs */
  for (int i=0;i<data_wk->ncells*(NUM_SPECIES + 1);i++) {
      rowPtrs_tmp[i+1] = (sunindextype) data_wk->rowPtrs_c[i+1];
  }
  BL_PROFILE_VAR_STOP(FillSpsData);


  BL_PROFILE_VAR("SparseJac", SpsJac);
  /* Temp vectors */
  realtype temp_save_lcl, temp;
  realtype massfrac[NUM_SPECIES];
  realtype Jmat_tmp[(NUM_SPECIES+1)*(NUM_SPECIES+1)];
  /* Save Jac from cell to cell if more than one */
  temp_save_lcl  = 0.0;
  for (int tid = 0; tid < data_wk->ncells; tid ++) {
      /* Offset in case several cells */
      int offset   = tid * (NUM_SPECIES + 1);
      int offset_J = tid * data_wk->NNZ;
      /* rho MKS */ 
      realtype rho = 0.0;
      for (int i = 0; i < NUM_SPECIES; i++){
          rho = rho + ydata[offset + i];
      }
      /* Yks */
      realtype rhoinv = 1.0 / rho;
      for (int i = 0; i < NUM_SPECIES; i++){
          massfrac[i] = ydata[offset + i] * rhoinv;
      }
      /* temp */
      temp = ydata[offset + NUM_SPECIES];
      /* Do we recompute Jac ? */
      if (fabs(temp - temp_save_lcl) > 1.0) {
          int consP = data_wk->ireactor_type == eint_rho ? 0 : 1;
          EOS::RTY2JAC(rho, temp, massfrac, Jmat_tmp, consP);
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
          int nbVals = data_wk->rowPtrs_c[i]-data_wk->rowPtrs_c[i - 1];
          for (int j = 0; j < nbVals; j++) {
              int idx = data_wk->colVals_c[ data_wk->rowPtrs_c[i - 1] + j ];
              Jdata[ offset_J + data_wk->rowPtrs_c[i - 1] + j ] = Jmat_tmp[(i - 1) + (NUM_SPECIES + 1)*idx];
          }
      }
  }
  BL_PROFILE_VAR_STOP(SpsJac);

  return(0);
}


#ifdef USE_KLU_PP 
/* Analytical SPARSE Jacobian evaluation */
int cJac_KLU(realtype tn, N_Vector u, N_Vector fu, SUNMatrix J,
             void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3){

  BL_PROFILE_VAR("SparseKLUJac", SpsKLUJac);
  /* Make local copies of pointers to input data (big M) */
  realtype *ydata           = N_VGetArrayPointer(u);
  sunindextype *colptrs_tmp = SUNSparseMatrix_IndexPointers(J);
  sunindextype *rowvals_tmp = SUNSparseMatrix_IndexValues(J);
  realtype *Jdata           = SUNSparseMatrix_Data(J);

  /* Make local copies of pointers in user_data (cell M)*/
  UserData data_wk = (UserData) user_data;

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
  /* Save Jac from cell to cell if more than one */
  temp_save_lcl = 0.0;
  for (int tid = 0; tid < data_wk->ncells; tid ++) {
      /* Offset in case several cells */
      int offset = tid * (NUM_SPECIES + 1);
      /* rho MKS */ 
      realtype rho = 0.0;
      for (int i = 0; i < NUM_SPECIES; i++){
          rho = rho + ydata[offset + i];
      }
      /* Yks */
      realtype rhoinv = 1.0 / rho;
      for (int i = 0; i < NUM_SPECIES; i++){
          massfrac[i] = ydata[offset + i] * rhoinv;
      }
      /* temp */
      temp = ydata[offset + NUM_SPECIES];
      /* Do we recompute Jac ? */
      if (fabs(temp - temp_save_lcl) > 1.0) {
          /* NRG CGS */
          int consP = data_wk->ireactor_type == eint_rho ? 0 : 1;
          EOS::RTY2JAC(rho, temp, massfrac, Jmat_tmp, consP);
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
      BL_PROFILE_VAR("DensetoSps", DtoS);
      for (int i = 1; i < NUM_SPECIES+2; i++) {
          int nbVals = data_wk->colPtrs[0][i]-data_wk->colPtrs[0][i - 1];
          for (int j = 0; j < nbVals; j++) {
              int idx = data_wk->rowVals[0][ data_wk->colPtrs[0][i - 1] + j ];
              Jdata[ data_wk->colPtrs[0][offset + i - 1] + j ] = Jmat_tmp[(i - 1) * (NUM_SPECIES + 1) + idx];
          }
      }
      BL_PROFILE_VAR_STOP(DtoS);
  }

  BL_PROFILE_VAR_STOP(SpsKLUJac);

  return(0);

}
#endif

/* Preconditioner setup routine for GMRES solver when custom sparse mode is activated 
 * Generate and preprocess P 
 */
int Precond_custom(realtype tn, N_Vector u, N_Vector fu, booleantype jok, 
                   booleantype *jcurPtr, realtype gamma, void *user_data)
{
  /* Make local copies of pointers to input data (big M) */
  realtype *udata   = N_VGetArrayPointer(u);
  /* Make local copies of pointers in user_data (cell M)*/
  UserData data_wk = (UserData) user_data;

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
      /* Save Jac from cell to cell if more than one */
      temp_save_lcl = 0.0;
      for (int tid = 0; tid < data_wk->ncells; tid ++) {
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
          EOS::RTY2C(rho, temp, massfrac, activity);
          /* Do we recompute Jac ? */
          if (fabs(temp - temp_save_lcl) > 1.0) {
              /* Formalism */
              int consP;
              if (data_wk->ireactor_type == eint_rho) {
                  consP = 0;
              } else {
                  consP = 1;
              }
              DWDOT_SIMPLIFIED(data_wk->JSPSmat[tid], activity, &temp, &consP);

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
      nbVals = data_wk->rowPtrs[0][i]-data_wk->rowPtrs[0][i-1];
      for (int j = 0; j < nbVals; j++) {
          /* row of non zero elem should be the same for all cells */
          int idx = data_wk->colVals[0][ data_wk->rowPtrs[0][i-1] + j ];
          /* Scale by -gamma */
          /* Add identity matrix */
          for (int tid = 0; tid < data_wk->ncells; tid ++) {
              if (idx == (i-1)) {
                  data_wk->Jdata[tid][ data_wk->rowPtrs[tid][i-1] + j ] = 1.0 - gamma * (data_wk->JSPSmat[tid])[ idx * (NUM_SPECIES+1) + idx]; 
              } else {
                  data_wk->Jdata[tid][ data_wk->rowPtrs[tid][i-1] + j ] = - gamma * (data_wk->JSPSmat[tid])[ (i - 1) + (NUM_SPECIES+1)*idx ]; 
              }
          }
      }
  }

  return(0);
}


#ifdef USE_KLU_PP 
/* Preconditioner setup routine for GMRES solver when KLU sparse mode is activated 
 * Generate and preprocess P
*/
int Precond_sparse(realtype tn, N_Vector u, N_Vector fu, booleantype jok, 
                   booleantype *jcurPtr, realtype gamma, void *user_data)
{
  /* Make local copies of pointers to input data (big M) */
  realtype *udata = N_VGetArrayPointer(u);
  /* Make local copies of pointers in user_data (cell M)*/
  UserData data_wk = (UserData) user_data;

  /* MW CGS */
  realtype molecular_weight[NUM_SPECIES];
  CKWT(molecular_weight);

  /* Check if Jac is stale */
  if (jok) {
      /* jok = SUNTRUE: Copy Jbd to P */
      *jcurPtr = SUNFALSE;
  } else {
      /* Temp vectors */
      realtype activity[NUM_SPECIES], massfrac[NUM_SPECIES];
      /* Save Jac from cell to cell if more than one */
      realtype temp_save_lcl = 0.0;
      for (int tid = 0; tid < data_wk->ncells; tid ++) {
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
          realtype temp = udata[offset + NUM_SPECIES];
          /* Activities */
          EOS::RTY2C(rho, temp, massfrac, activity);
          /* Do we recompute Jac ? */
          if (fabs(temp - temp_save_lcl) > 1.0) {
              /* Formalism */
              int consP = data_wk->ireactor_type == eint_rho ? 0 : 1;
              DWDOT_SIMPLIFIED(data_wk->JSPSmat[tid], activity, &temp, &consP);

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

  for (int i = 1; i < NUM_SPECIES+2; i++) {
      /* nb non zeros elem should be the same for all cells */
      int nbVals = data_wk->colPtrs[0][i]-data_wk->colPtrs[0][i-1];
      for (int j = 0; j < nbVals; j++) {
          /* row of non zero elem should be the same for all cells */
          int idx = data_wk->rowVals[0][ data_wk->colPtrs[0][i-1] + j ];
          /* Scale by -gamma */
          /* Add identity matrix */
          for (int tid = 0; tid < data_wk->ncells; tid ++) {
              if (idx == (i-1)) {
                  data_wk->Jdata[tid][ data_wk->colPtrs[tid][i-1] + j ] = 1.0 - gamma * (data_wk->JSPSmat[tid])[ idx * (NUM_SPECIES+1) + idx]; 
              } else {
                  data_wk->Jdata[tid][ data_wk->colPtrs[tid][i-1] + j ] = - gamma * (data_wk->JSPSmat[tid])[ (i-1) * (NUM_SPECIES+1) + idx ]; 
              }
          }
      }
  }
  
  BL_PROFILE_VAR("KLU_factorization", KLU_factor);
  if (!(data_wk->FirstTimePrecond)) {
      for (int tid = 0; tid < data_wk->ncells; tid ++) {
          int ok = klu_refactor(data_wk->colPtrs[tid], data_wk->rowVals[tid], data_wk->Jdata[tid], data_wk->Symbolic[tid], data_wk->Numeric[tid], &(data_wk->Common[tid]));
      }
  } else {
      for (int tid = 0; tid < data_wk->ncells; tid ++) {
          data_wk->Numeric[tid] = klu_factor(data_wk->colPtrs[tid], data_wk->rowVals[tid], data_wk->Jdata[tid], data_wk->Symbolic[tid], &(data_wk->Common[tid])) ; 
      }
      data_wk->FirstTimePrecond = false;
  }
  BL_PROFILE_VAR_STOP(KLU_factor);

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
  UserData data_wk = (UserData) user_data;
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
      EOS::RTY2C(rho, temp, massfrac, activity);
      /* jok = SUNFALSE: Generate Jbd from scratch and copy to P */
      /* Make local copies of problem variables, for efficiency. */
      int consP;
      if (data_wk->ireactor_type == eint_rho) { 
          consP = 0;
      } else {
          consP = 1;
      }
      DWDOT_SIMPLIFIED(Jmat, activity, &temp, &consP);

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

/* PSolve for GMRES solver when custom sparse mode is activated */
int PSolve_custom(realtype tn, N_Vector u, N_Vector fu, N_Vector r, N_Vector z,
                  realtype gamma, realtype delta, int lr, void *user_data)
{
  /* Make local copies of pointers in user_data */
  UserData data_wk = (UserData) user_data;

  /* Make local copies of pointers to input data (big M) */
  realtype *zdata = N_VGetArrayPointer(z);
  realtype *rdata = N_VGetArrayPointer(r);
  
  N_VScale(1.0, r, z);

  /* Solve the block-diagonal system Pz = r using LU factors stored
     in P and pivot data in pivot, and return the solution in z. */
  BL_PROFILE_VAR("GaussSolver", GaussSolver);
  for (int tid = 0; tid < data_wk->ncells; tid ++) {
      int offset          = tid * (NUM_SPECIES + 1);
      double *z_d_offset  = zdata  + offset;
      double *r_d_offset  = rdata  + offset;
      sgjsolve_simplified(data_wk->Jdata[tid], z_d_offset, r_d_offset);
  }
  BL_PROFILE_VAR_STOP(GaussSolver);

  return(0);
}

#ifdef USE_KLU_PP 
/* PSolve for GMRES solver when KLU sparse mode is activated */
int PSolve_sparse(realtype tn, N_Vector u, N_Vector fu, N_Vector r, N_Vector z,
                  realtype gamma, realtype delta, int lr, void *user_data)
{
  /* Make local copies of pointers in user_data */
  UserData data_wk = (UserData) user_data;

  /* Make local copies of pointers to input data (big M) */
  realtype *zdata = N_VGetArrayPointer(z);
  realtype *rdata = N_VGetArrayPointer(r);

  BL_PROFILE_VAR("KLU_inversion", PSolve_sparse);
  N_VScale(1.0, r, z);

  /* Solve the block-diagonal system Pz = r using LU factors stored
     in P and pivot data in pivot, and return the solution in z. */
  realtype zdata_cell[NUM_SPECIES+1];
  for (int tid = 0; tid < data_wk->ncells; tid ++) {
      int offset_beg = tid * (NUM_SPECIES + 1);
      int offset_end = (tid + 1) * (NUM_SPECIES + 1);
      std::memcpy(zdata_cell, zdata+offset_beg, (NUM_SPECIES+1)*sizeof(realtype));
      klu_solve(data_wk->Symbolic[tid], data_wk->Numeric[tid], NUM_SPECIES+1, 1, zdata_cell, &(data_wk->Common[tid])) ; 
      std::memcpy(zdata+offset_beg, zdata_cell, (NUM_SPECIES+1)*sizeof(realtype));
  }
  BL_PROFILE_VAR_STOP(PSolve_sparse);

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
  UserData data_wk = (UserData) user_data;
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
/**********************************/

/* 
 * CUSTOM SOLVER STUFF
 */
SUNLinearSolver SUNLinSol_sparse_custom(N_Vector y, SUNMatrix A, int reactor_type,
        int nsubsys, int subsys_size, int subsys_nnz) 
{
  SUNLinearSolver S;
  SUNLinearSolverContent_Sparse_custom content;

  /* Check that required arguments are not NULL */
  if (y == NULL || A == NULL) return(NULL);
  if (SUNMatGetID(A) != SUNMATRIX_SPARSE) return(NULL); 

  /* Matrix should be square */
  if (SUNSparseMatrix_Columns(A) != SUNSparseMatrix_Rows(A)) return(NULL);

  /* Check that it is a CSR matrix */
  if (SUNSparseMatrix_SparseType(A) != CSR_MAT) return(NULL);

  /* Matrix and vector dimensions must agree */
  if (N_VGetLength(y) != SUNSparseMatrix_Columns(A)) return(NULL);

  /* All subsystems must be the same size */ 
  if (SUNSparseMatrix_Columns(A) != (subsys_size * nsubsys)) return(NULL);

  /* Number of nonzeros per subsys must be the same */
  if (SUNSparseMatrix_NNZ(A) != (subsys_nnz * nsubsys)) return(NULL);

  /* Create an empty linear solver */
  S = SUNLinSolNewEmpty(); 
  if (S == NULL) return(NULL);

  /* Attach operations */ 
  S->ops->gettype    = SUNLinSolGetType_Sparse_custom;
  S->ops->solve      = SUNLinSolSolve_Sparse_custom;

  /* Create content */
  content = (SUNLinearSolverContent_Sparse_custom) malloc(sizeof *content);
  if (content == NULL) { SUNLinSolFree(S); return(NULL); }

  /* Attach content */
  S->content = content; 

  /* Fill content */ 
  content->last_flag    = 0;
  content->reactor_type = reactor_type;
  content->nsubsys      = nsubsys;
  content->subsys_size  = subsys_size;
  content->subsys_nnz   = subsys_nnz;

  return(S);
}


SUNLinearSolver_Type SUNLinSolGetType_Sparse_custom(SUNLinearSolver S) 
{
  return(SUNLINEARSOLVER_DIRECT);
}

int SUNLinSolSolve_Sparse_custom(SUNLinearSolver S, SUNMatrix A, N_Vector x,
        N_Vector b, realtype tol)
{
  realtype *x_d      = N_VGetArrayPointer(x);
  realtype *b_d      = N_VGetArrayPointer(b);

  double *Data = (double*) SUNSparseMatrix_Data(A);

  BL_PROFILE_VAR("GaussSolver", GaussSolver);
  for (int tid = 0; tid < SUN_CUSP_NUM_SUBSYS(S); tid ++) {
      int offset          = tid * SUN_CUSP_SUBSYS_NNZ(S);
      int offset_RHS      = tid * SUN_CUSP_SUBSYS_SIZE(S);
      double *Data_offset = Data + offset;
      double *x_d_offset  = x_d  + offset_RHS;
      double *b_d_offset  = b_d  + offset_RHS;
      sgjsolve(Data_offset, x_d_offset, b_d_offset);
  }
  BL_PROFILE_VAR_STOP(GaussSolver);

  return(SUNLS_SUCCESS);
}


/**********************************/
/* 
 * OTHERS
*/

void check_state(N_Vector yvec) 
{
  realtype *ydata = N_VGetArrayPointer(yvec);

  data->actual_ok_to_react = true;

  for (int tid = 0; tid < data->ncells; tid ++) {
      /* Offset in case several cells */
      int offset = tid * (NUM_SPECIES + 1);
      /* rho MKS */ 
      realtype rho = 0.0;
      for (int k = 0; k < NUM_SPECIES; k ++) {
          rho =  rho + ydata[offset + k];
      }
      /* temp */
      realtype Temp = ydata[offset + NUM_SPECIES];
      if ((rho < 1.0e-10) || (rho > 1.e10)) {
          data->actual_ok_to_react = false;
          Print() <<"rho "<< rho << "\n";
      }
      if ((Temp < 200.0) || (Temp > 5000.0)) {
          data->actual_ok_to_react = false; 
          Print() <<"Temp "<< Temp << "\n";
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

  if (data->isolve_type == dense_solve){
      flag = CVDlsGetNumRhsEvals(cvodeMem, &nfeLS);
      check_flag(&flag, "CVDlsGetNumRhsEvals", 1);
      flag = CVDlsGetNumJacEvals(cvodeMem, &nje);
      check_flag(&flag, "CVDlsGetNumJacEvals", 1);
  } else if (data->isolve_type == iterative_gmres_solve){
      flag = CVSpilsGetNumRhsEvals(cvode_mem, &nfeLS);
      check_flag(&flag, "CVSpilsGetNumRhsEvals", 1);
      flag = CVSpilsGetNumJtimesEvals(cvodeMem, &nje);
      //flag = CVSpilsGetNumJTSetupEvals(cvodeMem, &nje);
      check_flag(&flag, "CVSpilsGetNumJtimesEvals", 1);
      flag = CVSpilsGetNumPrecEvals(cvodeMem, &npe);
      check_flag(&flag, "CVSpilsGetNumPrecEvals", 1);
      flag = CVSpilsGetNumPrecSolves(cvodeMem, &nps);
      check_flag(&flag, "CVSpilsGetNumPrecSolves", 1);
      flag = CVSpilsGetNumLinIters(cvodeMem, &nli);
      check_flag(&flag, "CVSpilsGetNumLinIters", 1);
      flag = CVSpilsGetNumConvFails(cvodeMem, &ncfl); 
      check_flag(&flag, "CVSpilsGetNumConvFails", 1);
  }

  Print() << "-- Final Statistics --\n";
  Print() << "NonLinear (Newton) related --\n";
  Print() << Temp << " |DT(dt, dtcur) = " << nst << "(" << hlast << "," << hcur << "), RHS = " << nfe << ", Iterations = " << nni << ", ErrTestFails = " << netfails << ", LinSolvSetups = " << nsetups << "\n";
  if (data->isolve_type == dense_solve){
      Print() <<"Linear (Dense Direct Solve) related --\n";
      Print()<<Temp << " |FD RHS = "<< nfeLS<<", NumJacEvals = "<< nje <<" \n";
  } else if (data->isolve_type == iterative_gmres_solve){
      // LinSolvSetups actually reflects the number of time the LinSolver has been called. 
      // NonLinIterations can be taken without the need for LinItes
      Print() << "Linear (Krylov GMRES Solve) related --\n";
      Print() << Temp << " |RHSeval = "<< nfeLS << ", jtvEval = "<<nje << ", NumPrecEvals = "<< npe << ", NumPrecSolves = "<< nps <<"\n";
      Print() <<Temp << " |Iterations = "<< nli <<", ConvFails = "<< ncfl<<"\n";
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
      if (ParallelDescriptor::IOProcessor()) {
          fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
                  funcname);
      }
      return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
      errflag = (int *) flagvalue;
      if (*errflag < 0) {
          if (ParallelDescriptor::IOProcessor()) {
              fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
                      funcname, *errflag);
          }
          return(1); 
      }
  }

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
      if (ParallelDescriptor::IOProcessor()) {
          fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
                  funcname);
      }
      return(1); 
  }

  return(0);
}


/* Alloc Data for CVODE */
UserData AllocUserData(int reactor_type, int num_cells)
{
  /* Make local copies of pointers in user_data */
  UserData data_wk = (UserData) malloc(sizeof *data_wk);
  int omp_thread = 0;
#ifdef _OPENMP
  omp_thread = omp_get_thread_num();
#endif

  /* ParmParse from the inputs file: only done once */
  ParmParse pp("ode");
  pp.query("analytical_jacobian",data_wk->ianalytical_jacobian);
  data_wk->iverbose = 1;
  pp.query("verbose",data_wk->iverbose);

  std::string  solve_type_str = "none";
  ParmParse ppcv("cvode");
  ppcv.query("solve_type", solve_type_str);
  /* options are: 
  dense_solve           = 1;
  sparse_solve          = 5;
  iterative_gmres_solve = 99;
  sparse_solve_custom   = 101;
  iterative_gmres_solve_custom = 199;
  hack_dump_sparsity_pattern = -5;
  */
  if (solve_type_str == "dense") {
      data_wk->isolve_type = dense_solve; 
  } else if (solve_type_str == "sparse") {
      data_wk->isolve_type = sparse_solve;
  } else if (solve_type_str == "GMRES") {
      data_wk->isolve_type = iterative_gmres_solve;
  } else if (solve_type_str == "sparse_custom") {
      data_wk->isolve_type = sparse_solve_custom;
  } else if (solve_type_str == "GMRES_custom") { 
      data_wk->isolve_type = iterative_gmres_solve_custom;
  } else if (solve_type_str == "diag") {
      data_wk->isolve_type = hack_dump_sparsity_pattern;
  } else {
      Abort("Wrong solve_type. Options are: dense, sparse, GMRES, sparse_custom, GMRES_custom");
  }

  (data_wk->ireactor_type)             = reactor_type;

  (data_wk->ncells)                    = num_cells;

  (data_wk->FirstTimePrecond)          = true;
  (data_wk->reactor_cvode_initialized) = false;
  (data_wk->actual_ok_to_react)        = true; 

  int HP;
  if (data_wk->ireactor_type == eint_rho) {
      HP = 0;
  } else {
      HP = 1;
  }

#ifndef USE_KLU_PP
  if (data_wk->isolve_type == iterative_gmres_solve) {
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
  //} 

#else
  /* Sparse Direct and Sparse (It) Precond data */
  data_wk->colPtrs = new int*[data_wk->ncells];
  data_wk->rowVals = new int*[data_wk->ncells];
  data_wk->Jdata   = new realtype*[data_wk->ncells];

  if (data_wk->isolve_type == sparse_solve) {
      /* Sparse Matrix for Direct Sparse KLU solver */
      (data_wk->PS) = new SUNMatrix[1];
      SPARSITY_INFO(&(data_wk->NNZ),&HP,data_wk->ncells);
      if ((data_wk->iverbose > 0) && (omp_thread == 0)) {
          Print() << "--> SPARSE solver -- non zero entries: " << data_wk->NNZ << ", which represents "<< data_wk->NNZ/float((NUM_SPECIES+1) * (NUM_SPECIES+1) * (data_wk->ncells) * (data_wk->ncells)) *100.0 <<" % fill-in pattern\n";
      }
      (data_wk->PS)[0] = SUNSparseMatrix((NUM_SPECIES+1)*data_wk->ncells, (NUM_SPECIES+1)*data_wk->ncells, data_wk->NNZ, CSC_MAT);
      data_wk->colPtrs[0] = (int*) SUNSparseMatrix_IndexPointers((data_wk->PS)[0]); 
      data_wk->rowVals[0] = (int*) SUNSparseMatrix_IndexValues((data_wk->PS)[0]);
      data_wk->Jdata[0] = SUNSparseMatrix_Data((data_wk->PS)[0]);
      SPARSITY_PREPROC_CSC(data_wk->rowVals[0],data_wk->colPtrs[0],&HP,data_wk->ncells);

  } else if (data_wk->isolve_type == iterative_gmres_solve) {
      /* KLU internal storage */
      data_wk->Common   = new klu_common[data_wk->ncells];
      data_wk->Symbolic = new klu_symbolic*[data_wk->ncells];
      data_wk->Numeric  = new klu_numeric*[data_wk->ncells];
      /* Sparse Matrices for It Sparse KLU block-solve */
      data_wk->PS = new SUNMatrix[data_wk->ncells];
      /* Number of non zero elements*/
      SPARSITY_INFO_SYST_SIMPLIFIED(&(data_wk->NNZ),&HP);
      if ((data_wk->iverbose > 0) && (omp_thread == 0) && (data_wk->ianalytical_jacobian != 0)) {
          Print() << "--> SPARSE Preconditioner -- non zero entries: " << data_wk->NNZ << ", which represents "<< data_wk->NNZ/float((NUM_SPECIES+1) * (NUM_SPECIES+1)) *100.0 <<" % fill-in pattern\n";
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
          SPARSITY_PREPROC_SYST_SIMPLIFIED_CSC(data_wk->rowVals[i],data_wk->colPtrs[i],data_wk->indx,&HP);
          data_wk->JSPSmat[i] = new realtype[(NUM_SPECIES+1)*(NUM_SPECIES+1)];
          klu_defaults (&(data_wk->Common[i]));
          //data_wk->Common.btf = 0;
          //(data_wk->Common[i]).maxwork = 15;
          //data_wk->Common.ordering = 1;
          data_wk->Symbolic[i] = klu_analyze (NUM_SPECIES+1, data_wk->colPtrs[i], data_wk->rowVals[i], &(data_wk->Common[i])) ; 
      }
  //}
#endif

  } else if (data_wk->isolve_type == iterative_gmres_solve_custom) {
      /* Sparse Direct and Sparse (It) Precond data */
      data_wk->colVals = new int*[data_wk->ncells];
      data_wk->rowPtrs = new int*[data_wk->ncells];
      data_wk->Jdata   = new realtype*[data_wk->ncells];
      /* Matrices for It Sparse custom block-solve */
      data_wk->PS         = new SUNMatrix[data_wk->ncells];
      data_wk->JSPSmat    = new realtype*[data_wk->ncells];
      /* Number of non zero elements*/
      SPARSITY_INFO_SYST_SIMPLIFIED(&(data_wk->NNZ),&HP);
      for(int i = 0; i < data_wk->ncells; ++i) {
          (data_wk->PS)[i]       = SUNSparseMatrix(NUM_SPECIES+1, NUM_SPECIES+1, data_wk->NNZ, CSR_MAT);
          data_wk->rowPtrs[i]    = (int*) SUNSparseMatrix_IndexPointers((data_wk->PS)[i]);
          data_wk->colVals[i]    = (int*) SUNSparseMatrix_IndexValues((data_wk->PS)[i]);
          data_wk->Jdata[i]      = SUNSparseMatrix_Data((data_wk->PS)[i]);
          SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(data_wk->colVals[i],data_wk->rowPtrs[i],&HP,0);
          data_wk->JSPSmat[i]    = new realtype[(NUM_SPECIES+1)*(NUM_SPECIES+1)];
      }
      if ((data_wk->iverbose > 0) && (omp_thread == 0)) {
          Print() << "--> SPARSE Preconditioner -- non zero entries: " << data_wk->NNZ*data_wk->ncells << ", which represents "<< data_wk->NNZ/float((NUM_SPECIES+1) * (NUM_SPECIES+1) * data_wk->ncells) *100.0 <<" % fill-in pattern\n";
      }
  } else if (data_wk->isolve_type == sparse_solve_custom) {
      /* Number of non zero elements*/
      SPARSITY_INFO_SYST(&(data_wk->NNZ),&HP,1);
      data_wk->PSc          = SUNSparseMatrix((NUM_SPECIES+1)*data_wk->ncells, (NUM_SPECIES+1)*data_wk->ncells, data_wk->NNZ*data_wk->ncells, CSR_MAT);
      data_wk->rowPtrs_c    = (int*) SUNSparseMatrix_IndexPointers(data_wk->PSc); 
      data_wk->colVals_c    = (int*) SUNSparseMatrix_IndexValues(data_wk->PSc);
      SPARSITY_PREPROC_SYST_CSR(data_wk->colVals_c,data_wk->rowPtrs_c,&HP,data_wk->ncells,0);
      if ((data_wk->iverbose > 0) && (omp_thread == 0)) {
          Print() << "--> SPARSE solver -- non zero entries: " << data_wk->NNZ*data_wk->ncells << ", which represents "<< data_wk->NNZ/float((NUM_SPECIES+1) * (NUM_SPECIES+1) * data_wk->ncells) *100.0 <<" % fill-in pattern\n";
      }
  }  else if (data_wk->isolve_type == hack_dump_sparsity_pattern) {
      /* Debug mode, makes no sense to call with OMP/MPI activated */
      int counter;

      /* CHEMISTRY JAC */
      SPARSITY_INFO(&(data_wk->NNZ),&HP,1);
      Print() << "--> Chem Jac -- non zero entries: " << data_wk->NNZ << ", which represents "<< data_wk->NNZ/float((NUM_SPECIES+1) * (NUM_SPECIES+1)) *100.0 <<" % fill-in pattern\n";
      SUNMatrix PS;
      PS = SUNSparseMatrix((NUM_SPECIES+1), (NUM_SPECIES+1), data_wk->NNZ, CSR_MAT);
      int *colIdx, *rowCount;
      rowCount = (int*) SUNSparseMatrix_IndexPointers(PS); 
      colIdx   = (int*) SUNSparseMatrix_IndexValues(PS);
      SPARSITY_PREPROC_CSR(colIdx,rowCount,&HP,1, 0);
      std::cout <<" " << std::endl;
      std::cout << "*** Treating CHEM Jac (CSR symbolic analysis)***" << std::endl;
      std::cout <<" " << std::endl;
      int nbVals;
      counter = 0;
      for (int i = 0; i < NUM_SPECIES+1; i++) {
          nbVals         = rowCount[i+1] - rowCount[i];
          int idx_arr[nbVals];
          std::fill_n(idx_arr, nbVals, -1);
          std::memcpy(idx_arr, colIdx + rowCount[i], nbVals*sizeof(int));
          int idx        = 0;
          for (int j = 0; j < NUM_SPECIES+1; j++) {
              if ((j == idx_arr[idx]) && (nbVals > 0)) {
                  std::cout << 1 << " ";
                  idx = idx + 1;
                  counter = counter + 1;
              } else {
                  std::cout << 0 << " ";
              }
          }
          std::cout << std::endl;
      }
      std::cout << " There was " << counter << " non zero elems (compare to the "<<data_wk->NNZ<< " we need)" << std::endl;

      /* SYST JAC */
      SPARSITY_INFO_SYST(&(data_wk->NNZ),&HP,1);
      Print() << "--> Syst Jac -- non zero entries: " << data_wk->NNZ << ", which represents "<< data_wk->NNZ/float((NUM_SPECIES+1) * (NUM_SPECIES+1)) *100.0 <<" % fill-in pattern\n";
      PS = SUNSparseMatrix((NUM_SPECIES+1), (NUM_SPECIES+1), data_wk->NNZ, CSR_MAT);
      rowCount = (int*) SUNSparseMatrix_IndexPointers(PS); 
      colIdx   = (int*) SUNSparseMatrix_IndexValues(PS);
      SPARSITY_PREPROC_SYST_CSR(colIdx,rowCount,&HP,1,1);
      /* CHEMISTRY JAC */
      std::cout <<" " << std::endl;
      std::cout << "*** Treating SYST Jac (CSR symbolic analysis)***" << std::endl;
      std::cout <<" " << std::endl;
      counter = 0;
      for (int i = 0; i < NUM_SPECIES+1; i++) {
          nbVals         = rowCount[i+1] - rowCount[i];
          int idx_arr[nbVals];
          std::fill_n(idx_arr, nbVals, -1);
          std::memcpy(idx_arr, colIdx + (rowCount[i] - 1), nbVals*sizeof(int));
          int idx        = 0;
          for (int j = 0; j < NUM_SPECIES+1; j++) {
              if ((j == idx_arr[idx] - 1) && ((nbVals-idx) > 0)) {
                  std::cout << 1 << " ";
                  idx = idx + 1;
                  counter = counter + 1;
              } else {
                  std::cout << 0 << " ";
              }
          }
          std::cout << std::endl;
      }
      std::cout << " There was " << counter << " non zero elems (compare to the "<<data_wk->NNZ<< " we need)" << std::endl;

      /* SYST JAC SIMPLIFIED*/
      SPARSITY_INFO_SYST_SIMPLIFIED(&(data_wk->NNZ),&HP);
      Print() << "--> Simplified Syst Jac (for Precond) -- non zero entries: " << data_wk->NNZ << ", which represents "<< data_wk->NNZ/float((NUM_SPECIES+1) * (NUM_SPECIES+1)) *100.0 <<" % fill-in pattern\n";
      PS = SUNSparseMatrix((NUM_SPECIES+1), (NUM_SPECIES+1), data_wk->NNZ, CSR_MAT);
      rowCount = (int*) SUNSparseMatrix_IndexPointers(PS); 
      colIdx   = (int*) SUNSparseMatrix_IndexValues(PS);
      SPARSITY_PREPROC_SYST_SIMPLIFIED_CSR(colIdx,rowCount,&HP,1);
      /* CHEMISTRY JAC */
      std::cout <<" " << std::endl;
      std::cout << "*** Treating simplified SYST Jac (CSR symbolic analysis)***" << std::endl;
      std::cout <<" " << std::endl;
      counter = 0;
      for (int i = 0; i < NUM_SPECIES+1; i++) {
          nbVals         = rowCount[i+1] - rowCount[i];
          int idx_arr[nbVals];
          std::fill_n(idx_arr, nbVals, -1);
          std::memcpy(idx_arr, colIdx + (rowCount[i] - 1), nbVals*sizeof(int));
          int idx        = 0;
          for (int j = 0; j < NUM_SPECIES+1; j++) {
              if ((j == idx_arr[idx] - 1) && ((nbVals-idx) > 0)) {
                  std::cout << 1 << " ";
                  idx = idx + 1;
                  counter = counter + 1;
              } else {
                  std::cout << 0 << " ";
              }
          }
          std::cout << std::endl;
      }
      std::cout << " There was " << counter << " non zero elems (compare to the "<<data_wk->NNZ<< " we need)" << std::endl;

      Abort("Dump Sparsity Patern of different Jacobians in CSR format.");
  }

  return(data_wk);
}


/* Free memory */
void reactor_close(){

  CVodeFree(&cvode_mem);
  SUNLinSolFree(LS);

  if (data->isolve_type == dense_solve) {
      SUNMatDestroy(A);
  }

  N_VDestroy(y); 
  FreeUserData(data);
}


/* Free data memory 
 * Probably not complete, how about the stuff allocated in KLU mode ? */
void FreeUserData(UserData data_wk)
{
#ifndef USE_KLU_PP
  if (data_wk->isolve_type == iterative_gmres_solve) {
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


/* End of file  */
