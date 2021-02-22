#include "EOS.H"

namespace EOS {

  AMREX_GPU_DEVICE_MANAGED bool initialized=false;
  
  AMREX_GPU_DEVICE_MANAGED amrex::Real Tc[NUM_SPECIES];
  AMREX_GPU_DEVICE_MANAGED amrex::Real Bi[NUM_SPECIES];
  AMREX_GPU_DEVICE_MANAGED amrex::Real oneOverTc[NUM_SPECIES];
  AMREX_GPU_DEVICE_MANAGED amrex::Real sqrtOneOverTc[NUM_SPECIES];
  AMREX_GPU_DEVICE_MANAGED amrex::Real sqrtAsti[NUM_SPECIES];
  AMREX_GPU_DEVICE_MANAGED amrex::Real Fomega[NUM_SPECIES];
  
void init() 
{
  amrex::Real Asti[NUM_SPECIES];
  amrex::Real omega[NUM_SPECIES];
  EOS::initialized=true;
  GET_CRITPARAMS(EOS::Tc, Asti, EOS::Bi, omega);
  for (int ii = 0; ii<NUM_SPECIES; ii++){
    EOS::oneOverTc[ii] = 1.0/Tc[ii];
    EOS::sqrtOneOverTc[ii] = std::sqrt(oneOverTc[ii]);
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
