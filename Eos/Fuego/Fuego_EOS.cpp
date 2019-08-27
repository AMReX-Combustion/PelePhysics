#include "Fuego_EOS.H"
#include "mechanism.h"
#include <chemistry_file.H>
#include <cmath>



AMREX_GPU_HOST_DEVICE EOS::EOS()
{}

AMREX_GPU_HOST_DEVICE EOS::~EOS()
{}


AMREX_GPU_HOST_DEVICE
void EOS::eos_EY2T(amrex::Real *Y, amrex::Real E, amrex::Real T)
{
    int lierr=0; 
    GET_T_GIVEN_EY(&E, Y, &T, &lierr);

}

AMREX_GPU_HOST_DEVICE
void EOS::eos_T2EI(amrex::Real T, amrex::Real *ei)
{
    CKUMS(&T,  ei); 
}

AMREX_GPU_HOST_DEVICE
void EOS::eos_TY2Cv(amrex::Real T, amrex::Real *Y, amrex::Real *Cv)
{
    amrex::Real cvi[NUM_SPECIES]; 

    CKCVMS(&T,  cvi); 

    *Cv = 0.0; 
    for(int i = 0; i < NUM_SPECIES; ++i){
        *Cv += Y[i] * cvi[i];
    }
}

AMREX_GPU_HOST_DEVICE
void EOS::eos_HY2T(amrex::Real *Y, amrex::Real H, amrex::Real T)
{
    int lierr=0; 
    GET_T_GIVEN_HY(&H, Y, &T, &lierr);

}

AMREX_GPU_HOST_DEVICE
void EOS::eos_T2HI(amrex::Real T, amrex::Real *hi)
{
    CKHMS(&T,  hi); 
}

AMREX_GPU_HOST_DEVICE
void EOS::eos_TY2Cp(amrex::Real T, amrex::Real *Y, amrex::Real *Cp)
{
    amrex::Real cpi[NUM_SPECIES]; 

    CKCPMS(&T,  cpi); 

    *Cp = 0.0; 
    for(int i = 0; i < NUM_SPECIES; ++i){
        *Cp += Y[i] * cpi[i];
    }
}

AMREX_GPU_HOST_DEVICE
void EOS::eos_RTY2C(amrex::Real rho, amrex::Real T, amrex::Real *Y, amrex::Real *acti)
{
    CKYTCR(&rho, &T, Y, acti);
}

AMREX_GPU_HOST_DEVICE      
void EOS::eos_RTY2W(amrex::Real rho, amrex::Real T, amrex::Real *Y, amrex::Real *wdot)
{
    amrex::Real C[NUM_SPECIES]; 

    CKYTCR(&rho, &T, Y, C); 
    CKWC(&T, C, wdot);
}

/* Should not be here but I am not sure wether the CKYTCR needs to be wraped. */
AMREX_GPU_HOST_DEVICE      
void EOS::eos_RTY2JAC(amrex::Real rho, amrex::Real T, amrex::Real *Y, amrex::Real *Jac, int HP)
{
    amrex::Real C[NUM_SPECIES]; 

    CKYTCR(&rho, &T, Y, C); 
    DWDOT(Jac,C,&T,&HP);
}
