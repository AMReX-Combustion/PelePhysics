#include "Fuego_EOS.H"
#include "mechanism.h"
#include <chemistry_file.H>
#include <cmath>



AMREX_GPU_HOST_DEVICE EOS::EOS()
{}

AMREX_GPU_HOST_DEVICE EOS::~EOS()
{}

//AMREX_GPU_HOST_DEVICE 
//void EOS::eos_bottom()
//{
//    CKCVMS(&T,  cvi);
//    CKCPMS(&T,  cpi); 
//    CKHMS(&T,   hi);
//    cv = 0.e0, cp = 0.e0, h = 0.e0; 
//    for(int i = 0; i < NUM_SPECIES; ++i){
//         cv+=massfrac[i]*cvi[i];
//         cp+=massfrac[i]*cpi[i]; 
//         h +=massfrac[i]* hi[i]; 
//    }
//    amrex::Real Cvx = wbar*cv; 
//    gam1 = (Cvx + Ru)/Cvx; 
//    cs = std::sqrt(gam1*p/rho); 
//    dpdr_e = p/rho;
//    dpde = (gam1 - 1.0)*rho; 
//    s = 1.e0; 
//    dpdr = 0.e0;
//}
//
//
//AMREX_GPU_HOST_DEVICE
//void EOS::eos_wb()
//{
//    amrex::Real imw[NUM_SPECIES]; 
//    get_imw(imw);
//    amrex::Real summ =0.0; 
//#pragma unroll 
//    for(int i = 0; i < NUM_SPECIES; ++i) summ+= massfrac[i]*imw[i]; 
//    wbar = 1.0/summ; 
//}


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
void EOS::eos_RTY2W(amrex::Real rho, amrex::Real T, amrex::Real *Y, amrex::Real *wdot)
{
    amrex::Real C[NUM_SPECIES]; 

    CKYTCR(&rho, &T, Y, C); 
    CKWC(&T, C, wdot);
}


//AMREX_GPU_HOST_DEVICE
//void EOS::eos_EY2T()
//{
//    int lierr=0; 
//    //eos_wb();
//    GET_T_GIVEN_EY(&e, massfrac, &T, &lierr);
//    //CKUMS(&T,  ei); 
//    //CKPY(&rho, &T, massfrac, &p);
//
//    //eos_bottom(); 
//}
//
//AMREX_GPU_HOST_DEVICE
//void EOS::eos_rt()
//{
//    eos_wb(); 
//    CKPY(&rho, &T, massfrac, &p); 
//    CKUMS(&T, ei);
//    for(int i = 0; i < NUM_SPECIES; ++i) e+= massfrac[i]*ei[i]; 
//    
//    eos_bottom(); 
//}
//
//
//AMREX_GPU_HOST_DEVICE
//void EOS::eos_mpr2wdot(amrex::Real wdot[])
//{
//    CKWYR(&rho, &T, massfrac, wdot);
//    eos_rt(); 
//    amrex::Real mw[NUM_SPECIES]; 
//    get_mw(mw); 
//    for(int n = 0; n < NUM_SPECIES; n++){
//        wdot[n] *= mw[n];
//    }  
//}
//
//
//AMREX_GPU_HOST_DEVICE
//void EOS::eos_rp()
//{
//    eos_wb(); 
//    T = p*wbar/(rho*Ru);
//    CKUMS(&T,  ei);
//    e = 0.0;  
//#pragma unroll
//    for(int i = 0; i < NUM_SPECIES; ++i) e += massfrac[i]*ei[i]; 
//    eos_bottom();     
//}
//
//AMREX_GPU_HOST_DEVICE
//void EOS::eos_ytx()
//{
//    CKYTX(massfrac, molefrac); 
//}
//
//AMREX_GPU_HOST_DEVICE
//void EOS::eos_hi()
//{
//   CKHMS( &T,   hi);
//}
//AMREX_GPU_HOST_DEVICE
//void EOS::eos_cv()
//{
//    CKCVBS(&T, massfrac, &cv); 
//}
//
