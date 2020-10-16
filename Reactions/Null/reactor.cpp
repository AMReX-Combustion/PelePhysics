#include <reactor.h> 

using namespace amrex;

/**********************************/
/* Global Variables */
  int ireactor_type;
  int eint_rho = 1; // in/out = rhoE/rhoY
  int enth_rho = 2; // in/out = rhoH/rhoY 
/**********************************/

/**********************************/
/* Set or update typVals */
void SetTypValsODE(const std::vector<double>& ExtTypVals) {}

/* Set or update the rel/abs tolerances  */
void SetTolFactODE(double relative_tol,double absolute_tol) {}

/* Function to ReSet the tol of the cvode object directly */
void ReSetTolODE() {}

/* Initialization routine, called once at the begining of the problem */
int reactor_init(int reactor_type, int ode_ncells) { 
    ireactor_type = reactor_type;
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

    /* Initial time and time to reach after integration */
    Real time_init = time;

    /* Perform integration one cell at a time */
    ParallelFor(box,
    [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept
    {
        // Integration of external source term
        Real renergy_loc  = rEner_in(i,j,k,0) + rEner_src_in(i,j,k,0) * dt_react; 
        rEner_in(i,j,k,0) = renergy_loc;
        Real rY_loc[NUM_SPECIES];
        for (int n = 0; n < NUM_SPECIES; n++){
            rY_loc[n] = rY_in(i,j,k,n) + rY_src_in(i,j,k,n) * dt_react;
            rY_in(i,j,k,n) = rY_loc[n];
        }
        /* T update with energy and Y */
        // Get updated rho
        Real rho_loc = 0.0;
        for (int n = 0; n < NUM_SPECIES; n++) {
            rho_loc += rY_loc[n]; 
        }
        // Get updated Ys
        Real Y_loc[NUM_SPECIES];
        for (int n = 0; n < NUM_SPECIES; n++) {
            Y_loc[n] = rY_loc[n] / rho_loc;
        }
        // Get updated energy
        Real energy_loc = renergy_loc / rho_loc;
        // Get curr estimate of T
        Real T_loc    = T_in(i,j,k,0);
        if (ireactor_type == eint_rho){
            EOS::EY2T(energy_loc,Y_loc,T_loc);
        } else {
            EOS::HY2T(energy_loc,Y_loc,T_loc);
        }
        T_in(i,j,k,0) = T_loc;
        /* Dummy */
        FC_in(i,j,k,0) = 0.0;
    });  

#ifdef MOD_REACTOR
    /* If reactor mode is activated, update time to perform subcycling */
    time  = time_init + dt_react;
#endif
    /* Get estimate of how hard the integration process was */
    return 0;
}

/* Main routine for CVode integration: classic version */
int react(amrex::Real *rY_in, amrex::Real *rY_src_in, 
          amrex::Real *rX_in, amrex::Real *rX_src_in,
          amrex::Real &dt_react, amrex::Real &time) {

    /* Initial time and time to reach after integration */
    Real time_init = time;

#ifdef MOD_REACTOR
    /* If reactor mode is activated, update time */
    time  = time_init + dt_react;
#endif
    /* Get estimate of how hard the integration process was */
    return 0;
}

/* Free memory */
void reactor_close() {}


/* End of file  */
