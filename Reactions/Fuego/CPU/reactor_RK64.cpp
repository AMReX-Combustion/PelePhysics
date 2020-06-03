#include <CPU/reactor_RK64.h> 
#include <AMReX_ParmParse.H>
#include <chemistry_file.H>
#include "mechanism.h"
#include <EOS.H>

/**********************************/
UserData data      = NULL;
/* Global Variables */
/* energy */
/* REMOVE MAYBE LATER */
int eint_rho = 1; // in/out = rhoE/rhoY
int enth_rho = 2; // in/out = rhoH/rhoY 
double *rhoX_init   = NULL;
double *rhoXsrc_ext = NULL;
double *rYsrc       = NULL;
/**********************************/

#ifdef _OPENMP
#pragma omp threadprivate(data)
#pragma omp threadprivate(rhoX_init,rhoXsrc_ext,rYsrc)
#endif

/******************************************************************************************/
/* Initialization routine, called once at the begining of the problem */
int reactor_init(const int* reactor_type, const int* Ncells,double rk64_errtol,
        int rk64_nsubsteps_guess,int rk64_nsubsteps_min,int rk64_nsubsteps_max) {

#ifdef _OPENMP
    int omp_thread;
    /* omp thread if applicable */
    omp_thread = omp_get_thread_num(); 
#endif
    data = AllocUserData(*reactor_type, *Ncells,rk64_errtol,rk64_nsubsteps_guess,
            rk64_nsubsteps_min,rk64_nsubsteps_max);

    /* Define vectors to be used later in creact */
    rhoX_init   = (double *) malloc(data->ncells*sizeof(double));
    rhoXsrc_ext = (double *) malloc( data->ncells*sizeof(double));
    rYsrc       = (double *)  malloc((data->ncells*NUM_SPECIES)*sizeof(double));

    return(0);
}
/******************************************************************************************/
UserData AllocUserData(int reactor_type, int num_cells,double rk64_errtol,
        int rk64_nsubsteps_guess,int rk64_nsubsteps_min,int rk64_nsubsteps_max)
{
    /* Make local copies of pointers in user_data */
    UserData data_wk;
    data_wk = (UserData) malloc(sizeof *data_wk);
#ifdef _OPENMP
    int omp_thread;

    /* omp thread if applicable */
    omp_thread = omp_get_thread_num(); 
#endif

    (data_wk->ncells)        = num_cells;
    (data_wk->iverbose)      = 1;
    (data_wk->ireactor_type) = reactor_type;

    (data_wk->errtol) = rk64_errtol;
    (data_wk->nsubsteps_guess) = rk64_nsubsteps_guess;
    (data_wk->nsubsteps_min)   = rk64_nsubsteps_min;
    (data_wk->nsubsteps_max)   = rk64_nsubsteps_max;

    return(data_wk);
}
/******************************************************************************************/
/* Main call routine */
int react(double *rY_in, double *rY_src_in, 
        double *rX_in, double *rX_src_in,
        double *dt_react, double *time)
{

    double time_init = *time;
    double time_out  = *time + (*dt_react);
    double current_time = time_init;
    double *soln_reg,*carryover_reg,*error_reg,*zero_reg,*rhs;
    double dt_rk,dt_rk_min,dt_rk_max,change_factor;
    const double exp1=0.25;
    const double exp2=0.2;
    const double beta=1.0;
    const double tinyval=1e-50;

    int neq_tot        = (NUM_SPECIES + 1) * (data->ncells);

    const int nstages_rk64=6;
    const amrex::Real alpha_rk64[6] = {
        0.218150805229859,  //            3296351145737.0/15110423921029.0,
        0.256702469801519,  //            1879360555526.0/ 7321162733569.0,
        0.527402592007520,  //            10797097731880.0/20472212111779.0,
        0.0484864267224467, //            754636544611.0/15563872110659.0,
        1.24517071533530,   //            3260218886217.0/ 2618290685819.0,
        0.412366034843237,  //            5069185909380.0/12292927838509.0
    };

    const amrex::Real beta_rk64[6] = {
        -0.113554138044166,  //-1204558336989.0/10607789004752.0,
        -0.215118587818400,  //-3028468927040.0/14078136890693.0,
        -0.0510152146250577, //-455570672869.0/ 8930094212428.0,
        -1.07992686223881,   //-17275898420483.0/15997285579755.0,
        -0.248664241213447,  //-2453906524165.0/ 9868353053862.0,
        0.0};

    const amrex::Real err_rk64[6] = {
        -0.0554699315064507, //-530312978447.0/ 9560368366154.0,
        0.158481845574980,   // 473021958881.0/ 2984707536468.0,
        -0.0905918835751907, //-947229622805.0/10456009803779.0,
        -0.219084567203338,  //-2921473878215.0/13334914072261.0,
        0.164022338959433,   // 1519535112975.0/ 9264196100452.0,
        0.0426421977505659   // 167623581683.0/ 3930932046784.0
    };

#ifdef _OPENMP
    int omp_thread;

    /* omp thread if applicable */
    omp_thread = omp_get_thread_num(); 
#endif

#ifdef _OPENMP
    if ((data->iverbose > 1) && (omp_thread == 0)) 
    {
#else
    if (data->iverbose > 1) 
    {
#endif
	    amrex::Print() <<"\n -------------------------------------\n";
    }

#ifdef _OPENMP
    if ((data->iverbose > 3) && (omp_thread == 0)) 
    {
#else
    if (data->iverbose > 3) 
    {
#endif
	amrex::Print() <<"BEG : time curr is "<< time_init << 
            " and dt_react is " << *dt_react << 
            " and final time should be " << time_out << "\n";
    }


    soln_reg      = (double *) calloc(neq_tot,sizeof(double));
    carryover_reg = (double *) calloc(neq_tot,sizeof(double));
    error_reg     = (double *) calloc(neq_tot,sizeof(double));
    zero_reg      = (double *) calloc(neq_tot,sizeof(double));
    rhs           = (double *) calloc(neq_tot,sizeof(double));


    /* Get Device MemCpy of in arrays */
    /* Get Device pointer of solution vector */
    /* rhoY,T */
    std::memcpy(soln_reg,      rY_in, sizeof(double) * neq_tot);
    std::memcpy(carryover_reg, rY_in, sizeof(double) * neq_tot);
    /* rhoY_src_ext */
    std::memcpy(rYsrc, rY_src_in, (NUM_SPECIES * data->ncells)*sizeof(double));
    /* rhoE/rhoH */
    std::memcpy(rhoX_init, rX_in, sizeof(double) * data->ncells);
    std::memcpy(rhoXsrc_ext, rX_src_in, sizeof(double) * data->ncells);

    dt_rk     = *dt_react/double(data->nsubsteps_guess);
    dt_rk_min = *dt_react/double(data->nsubsteps_max);
    dt_rk_max = *dt_react/double(data->nsubsteps_min);

    int nsteps=0;

    while(current_time < time_out)
    {
        current_time += dt_rk;
        nsteps++;
        std::memcpy(carryover_reg, soln_reg, sizeof(double) * neq_tot);
        std::memcpy(error_reg,     zero_reg, sizeof(double) * neq_tot);

        for(int stage=0;stage<nstages_rk64;stage++)
        {
            std::memcpy(rhs,zero_reg, sizeof(double) * neq_tot);
            fKernelSpec(&current_time, soln_reg, rhs, data);

            for(int i=0;i<neq_tot;i++)
            {
                error_reg[i]    += err_rk64[stage]*dt_rk*rhs[i];
                soln_reg[i]      = carryover_reg[i] + alpha_rk64[stage]*dt_rk*rhs[i];
                carryover_reg[i] = soln_reg[i]      + beta_rk64[stage] *dt_rk*rhs[i];
            }
        } 

        //adapt time-step
        double max_err=tinyval;
        for(int i=0;i<neq_tot;i++)
        {
            if(fabs(error_reg[i]) > max_err)
            {
                max_err=fabs(error_reg[i]);
            }
        }

        //chance to increase time step
        if(max_err < data->errtol)
        {
            change_factor = beta*pow((data->errtol/max_err),exp1);
            dt_rk = std::min(dt_rk_max,dt_rk*change_factor);
        }
        //reduce time step (error is high!)
        else
        {
            change_factor=beta*pow((data->errtol/max_err),exp2);
            dt_rk = std::max(dt_rk_min,dt_rk*change_factor);
        }
    }

    //update guess from the current update
    (data->nsubsteps_guess) = nsteps;

    #ifdef MOD_REACTOR
    /* If reactor mode is activated, update time */
    *time  = time_init + (*dt_react);
    #endif

    #ifdef _OPENMP
    if ((data->iverbose > 3) && (omp_thread == 0)) 
    {
    #else
    if (data->iverbose > 3) 
    {
    #endif
          amrex::Print() <<"END : time curr is "<< current_time << 
              " and actual dt_react is " << *dt_react << "\n";
    }

    /* Pack data to return in main routine external */
    std::memcpy(rY_in, soln_reg, neq_tot * sizeof(double));
    for  (int i = 0; i < data->ncells; i++) 
    {
        rX_in[i] = rX_in[i] + (*dt_react) * rX_src_in[i];
    }

    return nsteps;
}
/******************************************************************************************/
/*
 * kernels
 */
/* RHS source terms evaluation */
void fKernelSpec(double *dt, double *yvec_d, double *ydot_d,  
        void *user_data)
{
    /* Make local copies of pointers in user_data (cell M)*/
    UserData data_wk;
    data_wk = (UserData) user_data;   

    /* Tmp vars */
    int tid;

    /* Loop on packed cells */
    for (tid = 0; tid < data_wk->ncells; tid ++) 
    {
        /* Tmp vars */
        double massfrac[NUM_SPECIES];
        double Xi[NUM_SPECIES];
        double cdot[NUM_SPECIES], molecular_weight[NUM_SPECIES];
        double cX;
        double temp, energy;
        /* EOS object in cpp */

        /* Offset in case several cells */
        int offset = tid * (NUM_SPECIES + 1); 

        /* MW CGS */
        CKWT(molecular_weight);

        /* rho MKS */ 
        double rho = 0.0;
        for (int i = 0; i < NUM_SPECIES; i++)
        {
            rho = rho + yvec_d[offset + i];
        }

        /* temp */
        temp = yvec_d[offset + NUM_SPECIES];

        /* Yks */
        for (int i = 0; i < NUM_SPECIES; i++)
        {
            massfrac[i] = yvec_d[offset + i] / rho;
        }
        

        /* NRG CGS */
        energy = (rhoX_init[tid] + rhoXsrc_ext[tid]*(*dt)) /rho;

        if (data_wk->ireactor_type == eint_rho)
        {
            /* UV REACTOR */
            EOS::EY2T(energy, massfrac, temp);
            EOS::TY2Cv(temp, massfrac, cX);
            EOS::T2Ei(temp, Xi);
        } 
        else 
        {
            /* HP REACTOR */
            EOS::HY2T(energy, massfrac, temp);
            EOS::TY2Cp(temp, massfrac, cX);
            EOS::T2Hi(temp, Xi);
        }
        EOS::RTY2WDOT(rho, temp, massfrac, cdot);

        /* Fill ydot vect */
        ydot_d[offset + NUM_SPECIES] = rhoXsrc_ext[tid];
        for (int i = 0; i < NUM_SPECIES; i++)
        {
            ydot_d[offset + i] = cdot[i] + rYsrc[tid * (NUM_SPECIES) + i];
            ydot_d[offset + NUM_SPECIES] = ydot_d[offset + NUM_SPECIES]  - ydot_d[offset + i] * Xi[i];
        }
        ydot_d[offset + NUM_SPECIES] = ydot_d[offset + NUM_SPECIES] /(rho * cX);
    }
}
/******************************************************************************************/
/* Free memory */
void reactor_close()
{

    FreeUserData(data);
    free(rhoX_init);
    free(rhoXsrc_ext);
    free(rYsrc);
}
/******************************************************************************************/
/* Free data memory */
void FreeUserData(UserData data_wk)
{
    free(data_wk);
} 
/******************************************************************************************/
