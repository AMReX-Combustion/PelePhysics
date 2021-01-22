#include <math.h>
#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <EOS.H>

/**********************************/
/* Functions Called by the Program */
int reactor_info(int cvode_iE, int Ncells);

int react(const amrex::Box& box,
          amrex::Array4<amrex::Real> const& rY_in,
          amrex::Array4<amrex::Real> const& rY_src_in, 
          amrex::Array4<amrex::Real> const& T_in, 
          amrex::Array4<amrex::Real> const& rEner_in,  
          amrex::Array4<amrex::Real> const& rEner_src_in,
          amrex::Array4<amrex::Real> const& FC_in,
          amrex::Array4<int> const& mask, 
          amrex::Real &dt_react,
          amrex::Real &time, 
          const int &reactor_type, 
          amrex::gpuStream_t stream);


int react(amrex::Real *rY_in, amrex::Real *rY_src_in, 
	      amrex::Real *rX_in, amrex::Real *rX_src_in, 
	      amrex::Real &dt_react, amrex::Real &time,
         int cvode_iE, int Ncells,
         cudaStream_t stream);

void reactor_close();


/**********************************/
/* Helper functions */

void SetTypValsODE(const std::vector<double>& ExtTypVals);

void SetTolFactODE(double relative_tol,double absolute_tol);
