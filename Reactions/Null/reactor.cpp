#include <reactor.h> 

using namespace amrex;

/**********************************/
/* Set or update typVals */
void SetTypValsODE(const std::vector<double>& ExtTypVals) {}

/* Set or update the rel/abs tolerances  */
void SetTolFactODE(double relative_tol,double absolute_tol) {}

/* Function to ReSet the tol of the cvode object directly */
void ReSetTolODE() {}

/* Initialization routine, called once at the begining of the problem */
int reactor_init(int reactor_type, int ode_ncells) { return(0); }

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

    /* Initial time and time to reach after integration */
    Real time_init = time;

#ifdef MOD_REACTOR
    /* If reactor mode is activated, update time to perform subcycling */
    time  = time_init + dt_react;
#endif
    /* Get estimate of how hard the integration process was */
    return 0;
}

/* Main routine for CVode integration: integrate a Box version 2*/
int react_2(const Box& box,
          Array4<Real> const& rY_in,
          Array4<Real> const& rY_src_in,
          Array4<Real> const& T_in,
          Array4<Real> const& rEner_in,
          Array4<Real> const& rEner_src_in,
          Array4<Real> const& FC_in,
          Array4<int> const& mask_in,
          Real &dt_react,
          Real &time) {

    /* Initial time and time to reach after integration */
    Real time_init = time;

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
