#include <math.h>

#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <EOS.H>

/**********************************/
/* Functions Called by the Program */
int reactor_init(int cvode_iE, int Ncells);
#ifdef AMREX_USE_GPU
int reactor_info(int cvode_iE, int Ncells);
#endif

int react(const amrex::Box& box,
          amrex::Array4<amrex::Real> const& rY_in,
          amrex::Array4<amrex::Real> const& rY_src_in, 
          amrex::Array4<amrex::Real> const& T_in, 
          amrex::Array4<amrex::Real> const& rEner_in,  
          amrex::Array4<amrex::Real> const& rEner_src_in,
          amrex::Array4<amrex::Real> const& FC_in,
          amrex::Array4<int> const& mask, 
          amrex::Real &dt_react,
          amrex::Real &time
#ifdef AMREX_USE_GPU
          , const int reactor_type
          , cudaStream_t stream
#endif
          );


int react(amrex::Real *rY_in, amrex::Real *rY_src_in, 
	      amrex::Real *rX_in, amrex::Real *rX_src_in, 
	      amrex::Real &dt_react, amrex::Real &time);

void reactor_close();


/**********************************/
/* Helper functions */

void SetTypValsODE(const std::vector<double>& ExtTypVals);

void SetTolFactODE(double relative_tol,double absolute_tol);

void ReSetTolODE();
