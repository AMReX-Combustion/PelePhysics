#include <reactor.h> 
#include <AMReX_ParmParse.H>
#include <chemistry_file.H>
#include "mechanism.h"
#include <EOS.H>
#include <AMREX_misc.H>

using namespace amrex;

/**********************************/
/* Global Variables */
N_Vector y         = NULL;
SUNLinearSolver LS = NULL;
SUNNonlinearSolver NLS = NULL;
SUNMatrix A        = NULL;
void *arkode_mem    = NULL;
/* User data */
UserData data      = NULL;
Real time_init    = 0.0;
Array<Real,NUM_SPECIES+1> typVals = {-1};
Real relTol       = 1.0e-6;
Real absTol       = 1.0e-10;
/* REMOVE MAYBE LATER */
int dense_solve           = 1;
int eint_rho = 1; // in/out = rhoE/rhoY
int enth_rho = 2; // in/out = rhoH/rhoY 

#ifdef _OPENMP
#pragma omp threadprivate(y,LS,NLS,A)
#pragma omp threadprivate(arkode_mem,data)
#pragma omp threadprivate(typVals)
#pragma omp threadprivate(time_init)
#pragma omp threadprivate(relTol,absTol)
#pragma omp threadprivate(dense_solve,eint_rho,enth_rho)
#endif
/**********************************/

/**********************************/
/* Initialization of typVals */
void SetTypValsODE(const std::vector<Real>& ExtTypVals) 
{
    int size_ETV = NUM_SPECIES + 1;
    Vector<std::string> kname;
    EOS::speciesNames(kname);
    int omp_thread = 0;

#ifdef _OPENMP
    omp_thread = omp_get_thread_num();
#endif

    for (int i=0; i<size_ETV-1; i++) 
    {
        typVals[i] = ExtTypVals[i];
    }
    typVals[size_ETV-1] = ExtTypVals[size_ETV-1];
    if (omp_thread == 0)
    {
        Print() << "Set the typVals in PelePhysics: \n  ";
        for (int i=0; i<size_ETV-1; i++) 
        {
            Print() << kname[i] << ":" << typVals[i] << "  ";
        }
        Print() << "Temp:"<< typVals[size_ETV-1] <<  " \n";
    }
}
//===========================================================
/* Set or update the rel/abs tolerances  */
void SetTolFactODE(Real relative_tol,Real absolute_tol) 
{
    relTol = relative_tol;
    absTol = absolute_tol;

    int omp_thread = 0;
#ifdef _OPENMP
    omp_thread = omp_get_thread_num();
#endif
/*
    if (omp_thread == 0)
    {
        Print() << "Set RTOL, ATOL = "<<relTol<< " "<<absTol<<  " in PelePhysics\n";
    }*/
}
//===========================================================
/* Function to ReSet the Tolerances */
void ReSetTolODE() 
{

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(data != NULL, "Reactor object is not initialized !!");

    int neq_tot;
    N_Vector atol;
    realtype *ratol;
    neq_tot = (NUM_SPECIES + 1) * data->ncells;
    atol    = N_VNew_Serial(neq_tot);
    ratol   = N_VGetArrayPointer(atol);

    int omp_thread = 0;
#ifdef _OPENMP
    omp_thread = omp_get_thread_num(); 
#endif

    int offset;
    if (typVals[0] > 0) 
    {
        if ((data->iverbose > 0) && (omp_thread == 0)) 
        {
            Print() << "Setting ARK/ERKODE tolerances rtol = " << relTol 
                << " atol = " << absTol << " in PelePhysics \n";
        }
        for  (int i = 0; i < data->ncells; i++) 
        {
            offset = i * (NUM_SPECIES + 1);
            for  (int k = 0; k < NUM_SPECIES + 1; k++) 
            {
                ratol[offset + k] = typVals[k]*absTol;
            }
        }
    } 
    else 
    {
        if ((data->iverbose > 0) && (omp_thread == 0)) 
        {
            Print() << "Setting ARK/ERKODE tolerances rtol = " << 
                relTol << " atol = " << absTol << " in PelePhysics \n";
        }
        for (int i=0; i<neq_tot; i++) 
        {
            ratol[i] = absTol;
        }
    }
    if (data->iuse_erkstep == 1) 
    {
        int flag = ERKStepSVtolerances(arkode_mem, relTol, atol); 
        if (check_flag(&flag, "ERKStepSVtolerances", 1)) Abort("Problem in ReSetTolODE");
    } else 
    {
        int flag = ARKStepSVtolerances(arkode_mem, relTol, atol); 
        if (check_flag(&flag, "ARKStepSVtolerances", 1)) Abort("Problem in ReSetTolODE");
    }

    N_VDestroy(atol);
}
//===========================================================
/* Initialization routine, called once at the begining of the problem */
int reactor_init(int reactor_type, int Ncells) 
{
    /* return Flag  */
    int flag;
    /* ARKODE initial time - 0 */
    realtype time;
    /* ARKODE tolerances */
    N_Vector atol; 
    realtype *ratol;
    /* Tot numb of eq to integrate */
    int neq_tot;
    int omp_thread = 0;
#ifdef _OPENMP
    omp_thread = omp_get_thread_num(); 
#endif
    /* Total number of eq to integrate */
    neq_tot        = (NUM_SPECIES + 1) * Ncells;

    /* Definition of main vector */
    y = N_VNew_Serial(neq_tot);
    if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);

    /* Does not work for more than 1 cell right now */
    data = AllocUserData(reactor_type, Ncells);
    if(check_flag((void *)data, "AllocUserData", 2)) return(1);

    /* Just a sanity check */
    if ( (data->iimplicit_solve == 1) && (data->iuse_erkstep == 1)) 
    {
        Abort("ERK ODE is for explicit updates, cannot do implicit");
    }

    /* Number of species and cells in mechanism */
    if ((data->iverbose > 0) && (omp_thread == 0)) 
    {
        Print() << "Number of species in mech is " << NUM_SPECIES << "\n";
        Print() << "Number of cells in one solve is " << data->ncells << "\n";
    }

    /* Create the solver memory, specify RHS function 
     * and set the pointer to user-defined data*/
    time = 0.0e+0;
    if (data->iimplicit_solve == 1) 
    {
        arkode_mem = ARKStepCreate(NULL, cF_RHS, time, y);
        if (check_flag((void *)arkode_mem, "ARKStepCreate", 0)) return 1;
        flag = ARKStepSetUserData(arkode_mem, data);
        if (check_flag(&flag, "ARKStepSetUserData", 1)) return 1;
    } 
    else 
    {
        if (data->iuse_erkstep == 1) 
        {
            arkode_mem = ERKStepCreate(cF_RHS, time, y);
            if (check_flag((void *)arkode_mem, "ERKStepCreate", 0)) return 1;
            flag = ERKStepSetUserData(arkode_mem, data);
            if (check_flag(&flag, "ERKStepSetUserData", 1)) return 1;
        } 
        else 
        {
            arkode_mem = ARKStepCreate(cF_RHS, NULL, time, y);
            if (check_flag((void *)arkode_mem, "ARKStepCreate", 0)) return 1;
            flag = ARKStepSetUserData(arkode_mem, data);
            if (check_flag(&flag, "ARKStepSetUserData", 1)) return 1;
        }
    }

    /* Definition of tolerances */
    atol  = N_VNew_Serial(neq_tot);
    ratol = N_VGetArrayPointer(atol);
    int offset;
    if (typVals[0] > 0) 
    {
        if ((data->iverbose > 0) && (omp_thread == 0)) 
        {
            Print() << "Setting ARK/ERKODE tolerances rtol = " << relTol 
                << " atolfact = " << absTol << " in PelePhysics \n";
        }
        for  (int i = 0; i < data->ncells; i++) 
        {
            offset = i * (NUM_SPECIES + 1);
            for  (int k = 0; k < NUM_SPECIES + 1; k++) 
            {
                ratol[offset + k] = typVals[k]*absTol;
            }
        }
    } 
    else 
    {
        if ((data->iverbose > 0) && (omp_thread == 0)) 
        {
            Print() << "Setting ARK/ERKODE tolerances rtol = " 
                << relTol << " atol = " << absTol << " in PelePhysics \n";
        }
        for (int i=0; i<neq_tot; i++) 
        {
            ratol[i] = absTol;
        }
    }
    if (data->iuse_erkstep == 1) 
    {
        flag = ERKStepSVtolerances(arkode_mem, relTol, atol); 
        if (check_flag(&flag, "ERKStepSVtolerances", 1)) return 1;
    } 
    else 
    {
        flag = ARKStepSVtolerances(arkode_mem, relTol, atol); 
        if (check_flag(&flag, "ARKStepSVtolerances", 1)) return 1;
    }

    if(data->iimplicit_solve == 1)
    {
        if ((data->iverbose > 0) && (omp_thread == 0)) 
        {
            Print() << "\n--> Using the ARKStep implicit solver: "<<
                "RK solver with Newton solve and Dense Direct Linear solve \n";
        }

        /* Create dense SUNMatrix for use in linear solves */
        A = SUNDenseMatrix(neq_tot, neq_tot);
        if(check_flag((void *)A, "SUNDenseMatrix", 0)) return(1);

        /* Create dense SUNLinearSolver object */
        LS = SUNLinSol_Dense(y, A);
        if(check_flag((void *)LS, "SUNLinSol_Dense", 0)) return(1);

        /* Linear solver interface */
        flag = ARKStepSetLinearSolver(arkode_mem, LS, A); 
        if (check_flag(&flag, "ARKStepSetLinearSolver", 1)) return(1);

        /* Setup Newton solve */
        NLS = SUNNonlinSol_Newton(y);
        if (check_flag((void *)NLS, "SUNNonlinSol_FixedPoint", 0)) return(1);

        flag = ARKStepSetNonlinearSolver(arkode_mem, NLS);
        if (check_flag(&flag, "ARKStepSetNonlinearSolver", 1)) return(1);

        flag = ARKStepSetMaxNonlinIters(arkode_mem, 100); 
        if (check_flag(&flag, "ARKStepSetMaxNonlinIters", 1)) return(1);

        if (data->ianalytical_jacobian == 0) 
        {
            if ((data->iverbose > 0) && (omp_thread == 0)) 
            {
                Print() << "    Without Analytical J/Preconditioner\n";
            }
        } 
        else 
        {
            if ((data->iverbose > 0) && (omp_thread == 0)) 
            {
                Print() << "    With Analytical J\n";
            }
            /* Set the user-supplied Jacobian routine Jac */
            flag = ARKStepSetJacFn(arkode_mem, cJac);                 /* Set Jacobian routine */
            if (check_flag(&flag, "ARKStepSetJacFn", 1)) return 1;
        }
    } 
    else 
    {
        if ((data->iverbose > 0) && (omp_thread == 0)) 
        {
            if (data->iuse_erkstep == 1) 
            {
                Print() << "\n--> Using the ERKStep explicit solver (RK solver) \n";
            } else 
            {
                Print() << "\n--> Using the ARKStep explicit solver (RK solver) \n";
            }
        }
    }

    /* Free the atol vector */
    N_VDestroy(atol); 

    // FIXME: Only output errors to IO rank, but turn off all error output for now.
    // Otherwise large runs can output GBs of text.
    ARKStepSetErrFile(arkode_mem, NULL);

    /* Ok we're done ...*/
    if ((data->iverbose > 1) && (omp_thread == 0)) 
    {
        Print() << "\n--> DONE WITH INITIALIZATION (CPU)" << data->ireactor_type << "\n";
    }

    /* Reactor is now initialized */
    data->reactor_arkode_initialized = true;

    return(0);
}
//===========================================================
/* react with array4s */
int react(const Box& box,
        Array4<Real> const& rY_in,
        Array4<Real> const& rY_src_in,
        Array4<Real> const& T_in,
        Array4<Real> const& rEner_in,
        Array4<Real> const& rEner_src_in,
        Array4<Real> const& FC_in,
        Array4<int> const& mask,
        Real &dt_react,
        Real &time) 
{

    int omp_thread = 0;
#ifdef _OPENMP
    omp_thread = omp_get_thread_num(); 
#endif

    if ((data->iverbose > 1) && (omp_thread == 0)) 
    {
        Print() <<"\n -------------------------------------\n";
    }

    /* Initial time and time to reach after integration */
    time_init = time;

    if ((data->iverbose > 3) && (omp_thread == 0)) 
    {
        Print() <<"BEG : time curr is "<< time_init << " and dt_react is " 
            << dt_react << " and final time should be " << time_init+dt_react << "\n";
    }

    int box_ncells  = box.numPts();
    data->boxcell   = 0; 

    if ((data->iverbose > 2) && (omp_thread == 0)) 
    {
        Print() <<"Ncells in the box = "<<  box_ncells  << "\n";
    }

    /* Perform integration one cell at a time */
    amrex::ParallelFor(box,
    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        if (mask(i,j,k) != -1)
        {
            Real mass_frac[NUM_SPECIES];
            Real rho = 0.0;
            Real rho_inv;
            Real Enrg_loc; //local energy
            Real temp;

            realtype *yvec_d      = N_VGetArrayPointer(y);

            BL_PROFILE_VAR("reactor::FlatStuff", FlatStuff);
            for (int n = 0; n < NUM_SPECIES; n++) 
            {
                yvec_d[n]        = rY_in(i,j,k,n);
                data->rYsrc[n]   = rY_src_in(i,j,k,n);
                rho += yvec_d[n]; 
            }
            rho_inv                 = 1.0 / rho;
            temp                    = T_in(i,j,k,0);
            data->rhoX_init[0]      = rEner_in(i,j,k,0); 
            data->rhoXsrc_ext[0]    = rEner_src_in(i,j,k,0);

            // T update with energy and Y 
            for (int n = 0; n < NUM_SPECIES; n++) 
            {
                mass_frac[n] = yvec_d[n] * rho_inv;
            }
            Enrg_loc = data->rhoX_init[0] / rho;
            if (data->ireactor_type == 1)
            {
                EOS::EY2T(Enrg_loc,mass_frac,temp);
            } else 
            {
                EOS::HY2T(Enrg_loc,mass_frac,temp);
            }
            yvec_d[NUM_SPECIES] = temp;
            BL_PROFILE_VAR_STOP(FlatStuff);


            if (data->iimplicit_solve == 1) 
            {
                //set explicit rhs to null
                ARKStepReInit(arkode_mem, NULL, cF_RHS, time_init, y);
            } 
            else 
            {
                if (data->iuse_erkstep == 1) 
                {
                    ERKStepReInit(arkode_mem, cF_RHS, time_init, y);
                } 
                else 
                {
                    //set implicit rhs to null
                    ARKStepReInit(arkode_mem, cF_RHS, NULL, time_init, y);
                }
            }

            // Time to reach after integration local to cell
            Real time_out_lcl  = time_init + dt_react;

            // Integration 
            Real dummy_time;
            BL_PROFILE_VAR("reactor::AroundARKODE", AroundARKODE);
            if (data->iuse_erkstep == 1) 
            { 
                // call integrator 
                int flag = ERKStepEvolve(arkode_mem, time_out_lcl, y, &dummy_time, ARK_NORMAL);   
            } 
            else 
            {
                // call integrator 
                int flag = ARKStepEvolve(arkode_mem, time_out_lcl, y, &dummy_time, ARK_NORMAL); 
            }
            BL_PROFILE_VAR_STOP(AroundARKODE);

            if ((data->iverbose > 1) && (omp_thread == 0)) 
            {
                Print() <<"Additional verbose info --\n";
                PrintFinalStats(arkode_mem, T_in(i,j,k,0));
                Print() <<"\n -------------------------------------\n";
            }

            // Get estimate of how hard the integration process was 
            long int nfe, nfi, nf;
            if (data->iuse_erkstep == 1) 
            {
                int flag = ERKStepGetNumRhsEvals(arkode_mem, &nfe);
                nf=nfe;
            } 
            else 
            {
                int flag = ARKStepGetNumRhsEvals(arkode_mem, &nfe, &nfi);
                nf=nfi;
            }

            FC_in(i,j,k,0) = nf;

            BL_PROFILE_VAR_START(FlatStuff);
            rho = 0.0;
            for (int n = 0; n < NUM_SPECIES; n++) 
            {
                rY_in(i,j,k,n) = yvec_d[n];
                rho += yvec_d[n]; 
            }
            rho_inv    = 1.0 / rho; 
            temp       = yvec_d[NUM_SPECIES];

            // T update with energy and Y 
            for (int n = 0; n < NUM_SPECIES; n++) 
            {
                mass_frac[n] = yvec_d[n] * rho_inv;
            }
            //find energy
            rEner_in(i,j,k,0) = data->rhoX_init[0] 
                + (dummy_time - time_init) * data->rhoXsrc_ext[0]; 

            Enrg_loc          = rEner_in(i,j,k,0) * rho_inv;
            if (data->ireactor_type == 1)
            {
                EOS::EY2T(Enrg_loc,mass_frac,temp);
            } else 
            {
                EOS::HY2T(Enrg_loc,mass_frac,temp);
            }
            T_in(i,j,k,0) = temp;
            BL_PROFILE_VAR_STOP(FlatStuff);

            if ((data->iverbose > 3) && (omp_thread == 0)) 
            {
                Print() <<"END : time curr is "<< dummy_time << " and actual "<< 
                    "dt_react is " << (dummy_time - time_init) << "\n";
            }
        } 
        else 
        {
            FC_in(i,j,k,0) = 0.0;
        }
    });


     /* Update dt_react with real time step taken ... */
#ifdef MOD_REACTOR
     /* If reactor mode is activated, update time */
     time  = time_init + dt_react;
#endif
    return(1);
}
//===========================================================
/* RHS routine */
int cF_RHS(realtype t, N_Vector y_in, N_Vector ydot_in, 
        void *user_data)
{

    realtype *y_d      = N_VGetArrayPointer(y_in);
    realtype *ydot_d   = N_VGetArrayPointer(ydot_in);

    fKernelSpec(&t, y_d, ydot_d,user_data);

    return(0);
}
//===========================================================
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
      energy = (data_wk->rhoX_init[data_wk->boxcell + tid] + data_wk->rhoXsrc_ext[data_wk->boxcell + tid] * dt) /rho;

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
//===========================================================
/* Free memory */
void reactor_close()
{

    if (data->iuse_erkstep == 1) 
    {
        ERKStepFree(&arkode_mem);    /* Free integrator memory */
    } 
    else 
    {
        ARKStepFree(&arkode_mem);    /* Free integrator memory */
        if (data->iimplicit_solve == 1) 
        {
            SUNLinSolFree(LS);
            SUNMatDestroy(A);
            SUNNonlinSolFree(NLS); 
        }
    }

    N_VDestroy(y); 
    FreeUserData(data);
}
//===========================================================
// react call without array4s
int react(realtype *rY_in, realtype *rY_src_in, 
        realtype *rX_in, realtype *rX_src_in,
        realtype &dt_react, realtype &time)
{

    int omp_thread = 0;
#ifdef _OPENMP
    omp_thread = omp_get_thread_num(); 
#endif

    if ((data->iverbose > 1) && (omp_thread == 0)) 
    {
        Print() <<"\n -------------------------------------\n";
    }

    /* Initial time and time to reach after integration */
    time_init = time;
    realtype time_out  = time + dt_react;

    if ((data->iverbose > 3) && (omp_thread == 0)) 
    {
        Print() <<"BEG : time curr is "<< time_init << " and dt_react is " 
            << dt_react << " and final time should be " << time_out << "\n";
    }

    /* Get Device MemCpy of in arrays */
    /* Get Device pointer of solution vector */
    realtype *yvec_d      = N_VGetArrayPointer(y);
    /* rhoY,T */
    std::memcpy(yvec_d, rY_in, sizeof(realtype) * ((NUM_SPECIES+1)*data->ncells));
    /* rhoY_src_ext */
    std::memcpy(data->rYsrc, rY_src_in, (NUM_SPECIES * data->ncells)*sizeof(Real));
    /* rhoE/rhoH */
    std::memcpy(data->rhoX_init, rX_in, sizeof(realtype) * data->ncells);
    std::memcpy(data->rhoXsrc_ext, rX_src_in, sizeof(realtype) * data->ncells);

    data->boxcell   = 0; 
    if (data->iimplicit_solve == 1) 
    {
        //set explicit rhs to null
        ARKStepReInit(arkode_mem, NULL, cF_RHS, time_init, y);
    } 
    else 
    {
        if (data->iuse_erkstep == 1) 
        {
            ERKStepReInit(arkode_mem, cF_RHS, time_init, y);
        } 
        else 
        {
            //set implicit rhs to null
            ARKStepReInit(arkode_mem, cF_RHS, NULL, time_init, y);
        }
    }

    Real dummy_time = 0;
    if (data->iuse_erkstep == 1) 
    { 
        //call integrator
        int flag = ERKStepEvolve(arkode_mem, time_out, y, &dummy_time, ARK_NORMAL);     
        if (check_flag(&flag, "ERKStepEvolve", 1)) return 1;
    } 
    else 
    {
        //call integrator
        int flag = ARKStepEvolve(arkode_mem, time_out, y, &dummy_time, ARK_NORMAL);
        if (check_flag(&flag, "ARKStepEvolve", 1)) return 1;
    }

    /* Update dt_react with real time step taken ... */
    dt_react = dummy_time - time_init;
#ifdef MOD_REACTOR
    /* If reactor mode is activated, update time */
    time  = time_init + dt_react;
#endif

    if ((data->iverbose > 3) && (omp_thread == 0)) 
    {
        Print() <<"END : time curr is "<< 
            dummy_time << " and actual dt_react is " << dt_react << "\n";
    }

    /* Pack data to return in main routine external */
    std::memcpy(rY_in, yvec_d, ((NUM_SPECIES+1)*data->ncells)*sizeof(realtype));
    for  (int i = 0; i < data->ncells; i++) 
    {
        rX_in[i] = rX_in[i] + dt_react * rX_src_in[i];
    }

    if ((data->iverbose > 1) && (omp_thread == 0)) 
    {
        Print() <<"Additional verbose info --\n";
        PrintFinalStats(arkode_mem, rY_in[NUM_SPECIES]);
        Print() <<"\n -------------------------------------\n";
    }

    /* Get estimate of how hard the integration process was */
    long int nfe, nfi, nf;
    if (data->iuse_erkstep == 1) 
    {
        int flag = ERKStepGetNumRhsEvals(arkode_mem, &nfe);
        nf=nfe;
    } 
    else 
    {
        int flag = ARKStepGetNumRhsEvals(arkode_mem, &nfe, &nfi);
        nf=nfi;
    }

    return nf;
}
//===========================================================
/* Analytical Jacobian evaluation */
int cJac(realtype tn, N_Vector u, N_Vector fu, SUNMatrix J,
        void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{

    /* Make local copies of pointers to input data (big M) */
    realtype *ydata  = N_VGetArrayPointer(u);

    /* Make local copies of pointers in user_data (cell M)*/
    UserData data_wk;
    data_wk = (UserData) user_data;   

    int tid;
    for (tid = 0; tid < data_wk->ncells; tid ++) 
    {
        /* Tmp vars */
        realtype *J_col_k;
        realtype massfrac[NUM_SPECIES], molecular_weight[NUM_SPECIES];
        realtype temp; 
        realtype Jmat_tmp[(NUM_SPECIES+1)*(NUM_SPECIES+1)];
        /* EOS object in cpp */

        /* Offset in case several cells */
        int offset = tid * (NUM_SPECIES + 1); 

        /* MW CGS */
        CKWT(molecular_weight);

        /* rho MKS */ 
        realtype rho = 0.0;
        for (int i = 0; i < NUM_SPECIES; i++)
        {
            rho = rho + ydata[offset + i];
        }

        /* temp */
        temp = ydata[offset + NUM_SPECIES];

        /* Yks */
        for (int i = 0; i < NUM_SPECIES; i++)
        {
            massfrac[i] = ydata[offset + i] / rho;
        }

        /* Jac */
        int consP;
        if (data_wk->ireactor_type == eint_rho) 
        {
            consP = 0;
        } else 
        {
            consP = 1;
        }
        EOS::RTY2JAC(rho, temp, massfrac, Jmat_tmp, consP);
        /* fill the sunMat */
        for (int k = 0; k < NUM_SPECIES; k++)
        {
            J_col_k = SM_COLUMN_D(J,offset + k);
            for (int i = 0; i < NUM_SPECIES; i++)
            {
                J_col_k[offset + i] = Jmat_tmp[k*(NUM_SPECIES+1)+i] 
                    * molecular_weight[i] / molecular_weight[k]; 
            }
            J_col_k[offset + NUM_SPECIES] = Jmat_tmp[k*(NUM_SPECIES+1)+NUM_SPECIES] 
                / molecular_weight[k]; 
        }
        J_col_k = SM_COLUMN_D(J,offset + NUM_SPECIES);
        for (int i = 0; i < NUM_SPECIES; i++)
        {
            J_col_k[offset + i] = Jmat_tmp[NUM_SPECIES*(NUM_SPECIES+1)+i] * molecular_weight[i]; 
        }
    }

    return(0);

}
//===========================================================
/* 
 * OTHERS
 */
/* Get and print some final statistics */
void PrintFinalStats(void *arkode_mem, realtype Temp)
{
    long int nst, nst_a, nfe, nfi, nsetups, nje, nfeLS, nni, ncfn, netf;
    int flag;

    if (data->iuse_erkstep == 1) 
    {
        flag = ERKStepGetNumSteps(arkode_mem, &nst);
        check_flag(&flag, "ERKStepGetNumSteps", 1);
        flag = ERKStepGetNumStepAttempts(arkode_mem, &nst_a);
        check_flag(&flag, "ERKStepGetNumStepAttempts", 1);
        flag = ERKStepGetNumRhsEvals(arkode_mem, &nfe);
        check_flag(&flag, "ERKStepGetNumRhsEvals", 1);

        Print() << "\nFinal Solver Statistics:\n";
        Print() << "   Internal solver steps = " << nst << " (attempted = " << nst_a << ")\n";
        Print() << "   Total RHS evals:  Fe = " << nfe << "\n";
    } 
    else 
    {
        flag = ARKStepGetNumSteps(arkode_mem, &nst);
        check_flag(&flag, "ARKStepGetNumSteps", 1);
        flag = ARKStepGetNumStepAttempts(arkode_mem, &nst_a);
        check_flag(&flag, "ARKStepGetNumStepAttempts", 1);
        flag = ARKStepGetNumRhsEvals(arkode_mem, &nfe, &nfi);
        check_flag(&flag, "ARKStepGetNumRhsEvals", 1);
        flag = ARKStepGetNumLinSolvSetups(arkode_mem, &nsetups);
        check_flag(&flag, "ARKStepGetNumLinSolvSetups", 1);
        flag = ARKStepGetNumErrTestFails(arkode_mem, &netf);
        check_flag(&flag, "ARKStepGetNumErrTestFails", 1);
        flag = ARKStepGetNumNonlinSolvIters(arkode_mem, &nni);
        check_flag(&flag, "ARKStepGetNumNonlinSolvIters", 1);
        flag = ARKStepGetNumNonlinSolvConvFails(arkode_mem, &ncfn);
        check_flag(&flag, "ARKStepGetNumNonlinSolvConvFails", 1);
        flag = ARKStepGetNumJacEvals(arkode_mem, &nje);
        check_flag(&flag, "ARKStepGetNumJacEvals", 1);
        flag = ARKStepGetNumLinRhsEvals(arkode_mem, &nfeLS);
        check_flag(&flag, "ARKStepGetNumLinRhsEvals", 1);

        Print() << "\nFinal Solver Statistics:\n";
        Print() << "   Internal solver steps = " << nst << " (attempted = " << nst_a << ")\n";
        Print() << "   Total RHS evals:  Fe = " << nfe << ",  Fi = " << nfi << "\n";
        Print() << "   Total linear solver setups = " << nsetups << "\n";
        Print() << "   Total RHS evals for setting up the linear system = " << nfeLS << "\n";
        Print() << "   Total number of Jacobian evaluations = " << nje << "\n";
        Print() << "   Total number of Newton iterations = " << nni << "\n";
        Print() << "   Total number of linear solver convergence failures = " << ncfn << "\n";
        Print() << "   Total number of error test failures = " << netf << "\n\n";
    }

}
//===========================================================
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
    if (opt == 0 && flagvalue == NULL) 
    {
        if (ParallelDescriptor::IOProcessor()) 
        {
            fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
                    funcname);
        }
        return(1); 
    }

    /* Check if flag < 0 */
    else if (opt == 1) 
    {
        errflag = (int *) flagvalue;
        if (*errflag < 0) 
        {
            if (ParallelDescriptor::IOProcessor()) 
            {
                fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
                        funcname, *errflag);
            }
            return(1); 
        }
    }

    /* Check if function returned NULL pointer - no memory allocated */
    else if (opt == 2 && flagvalue == NULL) 
    {
        if (ParallelDescriptor::IOProcessor()) 
        {
            fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
                    funcname);
        }
        return(1); 
    }

    return(0);
}
//===========================================================
UserData AllocUserData(int reactor_type, int num_cells)
{
    /* Make local copies of pointers in user_data */
    UserData data_wk = (UserData) malloc(sizeof *data_wk);
    int omp_thread = 0;

#ifdef _OPENMP
    omp_thread = omp_get_thread_num();
    if(omp_thread==0)
    {
       Print() << "Allocating data for ARKODE\n";
    }
#else
    Print() << "Allocating data for ARKODE\n";
#endif

    /* ParmParse from the inputs file */
    /* TODO change that in the future */ 
    ParmParse pp("ode");
    pp.query("analytical_jacobian",data_wk->ianalytical_jacobian);
    pp.query("implicit_solve", data_wk->iimplicit_solve);
    pp.query("use_erkstep", data_wk->iuse_erkstep);
    pp.query("rtol",relTol);
    pp.query("atol",absTol);

    (data_wk->ireactor_type)      = reactor_type;

    (data_wk->ncells)                    = num_cells;

    (data_wk->iverbose)                  = 1;

    (data_wk->reactor_arkode_initialized) = false;

    (data_wk->Yvect_full)  = new  amrex::Real[data_wk->ncells*(NUM_SPECIES+1)];
    (data_wk->rYsrc)       = new  amrex::Real[data_wk->ncells*(NUM_SPECIES)];
    (data_wk->rhoX_init)   = new  amrex::Real[data_wk->ncells];
    (data_wk->rhoXsrc_ext) = new  amrex::Real[data_wk->ncells];
    (data_wk->FCunt)       = new  int[data_wk->ncells];
    (data_wk->mask)        = new  int[data_wk->ncells];
    
    SetTolFactODE(relative_tol,absolute_tol);

    return(data_wk);
}
//===========================================================
/* Free data memory */
void FreeUserData(UserData data_wk)
{
    free(data_wk);
} 
//===========================================================
