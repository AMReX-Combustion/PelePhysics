#include "EOS.H"

namespace EOS {

  bool initialized=false;
  
  amrex::Real Tc[NUM_SPECIES];
  amrex::Real Bi[NUM_SPECIES];
  amrex::Real oneOverTc[NUM_SPECIES];
  amrex::Real sqrtAsti[NUM_SPECIES];
  amrex::Real Fomega[NUM_SPECIES];
  
void init() 
{
  amrex::Real Asti[NUM_SPECIES];
  amrex::Real omega[NUM_SPECIES];
  EOS::initialized=true;
  GET_CRITPARAMS(EOS::Tc, Asti, EOS::Bi, omega);
  for (int ii = 0; ii<NUM_SPECIES; ii++){
    EOS::oneOverTc[ii] = 1.0/Tc[ii];
    EOS::sqrtAsti[ii] = std::sqrt(Asti[ii]);
    EOS::Fomega[ii] = EOS::f0 + omega[ii]*(EOS::f1 + EOS::f2*omega[ii]);
  }
    
  CKINIT();
}

void close()
{
  CKFINALIZE();
}

} // namespace EOS
