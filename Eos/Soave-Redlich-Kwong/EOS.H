#ifndef _EOS_H_
#define _EOS_H_

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>

#include "mechanism.h"
#include "chemistry_file.H"

/* List of Speedups
 Save oneOverSqrtTc instead of oneOverTc
 Put diagonal of double loops into single loops, only do upper triang of double loop
 Precalculate Amloc and Amlocder
 Parmparse MaxIter and ConvCrit
 Raise error if fail MaxIter or ConvCrit
*/

namespace EOS {

  // Are these needed for SRK?
  constexpr amrex::Real RU   = 8.31446261815324e7;
  constexpr amrex::Real RUC  = 1.98721558317399615845;
  constexpr amrex::Real PATM = 1.01325e+06;

  // Constants for SRK
  constexpr amrex::Real f0 = 0.48508e+0;
  constexpr amrex::Real f1 = 1.5517e+0;
  constexpr amrex::Real f2 = -0.151613e+0;
  constexpr amrex::Real convCrit = 1e-4;
  constexpr int maxIter  = 2000;
  
  extern bool initialized;
  extern amrex::Real Tc[NUM_SPECIES];
  extern amrex::Real Bi[NUM_SPECIES];
  extern amrex::Real oneOverTc[NUM_SPECIES];
  extern amrex::Real sqrtAsti[NUM_SPECIES];
  extern amrex::Real Fomega[NUM_SPECIES];
  
void init();

void close();

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
MixingRuleAmBm(amrex::Real T, amrex::Real Y[], amrex::Real &am, amrex::Real &bm)
{
  am = 0.0;
  bm = 0.0;
  amrex::Real Tr ;
  amrex::Real amloc[NUM_SPECIES];

  // Combine as follows: add diagonal to am in first loop, loop upper triang in second loop
  
  for (int ii = 0; ii<NUM_SPECIES; ii++) {
    Tr = T * oneOverTc[ii];
    amloc[ii] = (1.0 + Fomega[ii] *(1.0 - std::sqrt(Tr))) * sqrtAsti[ii];
    bm += Y[ii] * Bi[ii];
  }

  for (int ii = 0; ii<NUM_SPECIES; ii++) {
    for (int jj = 0; jj<NUM_SPECIES; jj++) {
      am += Y[ii] * Y[jj] * amloc[ii] * amloc[jj];
    }
  }
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
MixingRuleAm(amrex::Real T, amrex::Real Y[], amrex::Real &am)
{
  am = 0.0;
  amrex::Real Tr ;
  amrex::Real amloc[NUM_SPECIES];

  // Combine as follows: add diagonal to am in first loop, loop upper triang in second loop
  
  for (int ii = 0; ii<NUM_SPECIES; ii++) {
    Tr = T * oneOverTc[ii];
    amloc[ii] = (1.0 + Fomega[ii] *(1.0 - std::sqrt(Tr))) * sqrtAsti[ii];
  }

  for (int ii = 0; ii<NUM_SPECIES; ii++) {
    for (int jj = 0; jj<NUM_SPECIES; jj++) {
      am += Y[ii] * Y[jj] * amloc[ii] * amloc[jj];
    }
  }
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
MixingRuleBm(amrex::Real Y[], amrex::Real &bm)
{
  bm = 0.0;
  for (int ii = 0; ii<NUM_SPECIES; ii++) {
    bm += Y[ii] * Bi[ii];
  }  
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
Calc_dAmdT(amrex::Real T, amrex::Real Y[], amrex::Real &dAmdT)
{
  dAmdT = 0.0;
  amrex::Real oneOverT = 1.0/T;
  amrex::Real Tr ;
  amrex::Real amloc[NUM_SPECIES];
  amrex::Real amlocder[NUM_SPECIES];

  for (int ii = 0; ii<NUM_SPECIES; ii++) {
    Tr = T * oneOverTc[ii];
    amloc[ii] = (1.0 + Fomega[ii]*(1.0 - std::sqrt(Tr)))*sqrtAsti[ii];
    amlocder[ii] = -0.5 * Fomega[ii] * sqrtAsti[ii]
      * oneOverT * oneOverTc[ii] * std::sqrt(T*Tc[ii])  ;
  }

  for (int ii = 0; ii<NUM_SPECIES; ii++) {
    for (int jj = 0; jj<NUM_SPECIES; jj++) {
      dAmdT += Y[ii] * Y[jj] * (amloc[ii] * amlocder[jj] + amloc[jj]*amlocder[ii]);
      //      std::cout << "test" << dAmdT << " " << Y[ii] << " " << Y[jj] << " " << amloc[ii]<< " " << amloc[ii]<< " " << amlocder[ii]<< " " << amlocder[ii] << std::endl;
    }
  } 
}
  

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
Calc_dAmdY(amrex::Real T, amrex::Real Y[], amrex::Real dAmdY[])
{
  amrex::Real Tr ;
  amrex::Real amloc[NUM_SPECIES];

  for (int ii = 0; ii<NUM_SPECIES; ii++) {
    dAmdY[ii] = 0.0;
    Tr = T * oneOverTc[ii];
    amloc[ii] = (1.0 + Fomega[ii]*(1.0 - std::sqrt(Tr)))*sqrtAsti[ii];
  }

  for (int ii = 0; ii<NUM_SPECIES; ii++) {
    for (int jj = 0; jj<NUM_SPECIES; jj++) {
      dAmdY[ii] += 2.0 * Y[jj] * (amloc[ii] * amloc[jj]);
    }
  }
}


AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
Calc_d2AmdY2(amrex::Real T, amrex::Real Y[], amrex::Real d2AmdY2[][NUM_SPECIES])
{
  amrex::Real Tr ;
  amrex::Real amloc[NUM_SPECIES];

  for (int ii = 0; ii<NUM_SPECIES; ii++) {
    Tr = T * oneOverTc[ii];
    amloc[ii] = (1.0 + Fomega[ii]*(1.0 - std::sqrt(Tr)))*sqrtAsti[ii];
  }

  for (int ii = 0; ii<NUM_SPECIES; ii++) {
    for (int jj = 0; jj<NUM_SPECIES; jj++) {
      d2AmdY2[ii][jj] = 2.0 * amloc[ii] * amloc[jj] ;
    }
  }
}


AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
Calc_d2AmdTY(amrex::Real T, amrex::Real Y[], amrex::Real d2AmdTY[])
{
  amrex::Real oneOverT = 1.0/T;
  amrex::Real Tr ;
  amrex::Real amloc[NUM_SPECIES];
  amrex::Real amlocder[NUM_SPECIES];

  for (int ii = 0; ii<NUM_SPECIES; ii++) {
    d2AmdTY[ii] = 0.0 ;
    Tr = T * oneOverTc[ii];
    amloc[ii] = (1.0 + Fomega[ii]*(1.0 - std::sqrt(Tr)))*sqrtAsti[ii];
    amlocder[ii] = -0.5 * Fomega[ii] * sqrtAsti[ii]
      * oneOverT * oneOverTc[ii] * std::sqrt(T*Tc[ii])  ;
  }

  for (int ii = 0; ii<NUM_SPECIES; ii++) {
    for (int jj = 0; jj<NUM_SPECIES; jj++) {
      d2AmdTY[ii] += 2.0 * Y[jj] * (amloc[ii] * amlocder[jj] + amloc[jj]*amlocder[ii]);
    }
  }
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
Calc_CompressFactor_Z(amrex::Real &Z,
		      amrex::Real am,
		      amrex::Real bm,
		      amrex::Real P,
		      amrex::Real T,
		      amrex::Real Wbar)
{
  amrex::Real theta, Z1, Z2, Z3, sqrtQ,third;
  amrex::Real RmT = EOS::RU/Wbar*T;
  amrex::Real B1 = bm*P/RmT;
  amrex::Real R1 = RmT;
  amrex::Real R2 = R1*RmT;
  amrex::Real R3 = R2*RmT;
  amrex::Real alpha = -1.0;
  amrex::Real beta = (am*P - bm*P*bm*P)/R2 - B1;
  amrex::Real gamma = -(am*bm*P*P)/R3;
  amrex::Real Q = (alpha * alpha - 3.0*beta)/9.0;
  amrex::Real R = (2.0*alpha*alpha*alpha - 9.0*alpha*beta + 27.0*gamma)/54.0;
  amrex::Real PI = 3.1415926535897932;
  // Multiple roots of cubic
  third = 1.0/3.0;
  if ((Q*Q*Q - R*R) > 0) {
    sqrtQ = std::sqrt(Q);
    theta = std::acos(R/(Q*sqrtQ));
    Z1 = -2.0 * sqrtQ *std::cos(theta*third) - alpha*third;
    Z2 = -2.0 * sqrtQ *std::cos((theta+ 2.0*PI)*third) - alpha*third;
    Z3 = -2.0 * sqrtQ *std::cos((theta+ 4.0*PI)*third) - alpha*third;
    Z = std::max(Z1,Z2);
    Z = std::max(Z ,Z3);
  } else {
    Z = -amrex::Math::copysign(1.0,R) *
      ( std::pow((std::sqrt(R*R - Q*Q*Q) + std::abs(R)),third) +
	Q/(pow(std::sqrt(R*R - Q*Q*Q) + std::abs(R),third )) ) - alpha*third ;
  }
}


AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
Calc_d2AmdT2(amrex::Real T, amrex::Real Y[], amrex::Real &d2AmdT2)
{
  amrex::Real oneOverT = 1.0/T;
  amrex::Real tmp1 = -0.5 *oneOverT;
  amrex::Real Tr ;
  amrex::Real amloc[NUM_SPECIES];
  amrex::Real amlocder[NUM_SPECIES];

  d2AmdT2 = 0.0;
  
  for (int ii = 0; ii<NUM_SPECIES; ii++) {
    Tr = T * oneOverTc[ii];
    amloc[ii] = (1.0 + Fomega[ii]*(1.0 - std::sqrt(Tr)))*sqrtAsti[ii];
    amlocder[ii] = -0.5 * Fomega[ii] * sqrtAsti[ii]
      * oneOverT * oneOverTc[ii] * std::sqrt(T*Tc[ii])  ;
  }

  for (int ii = 0; ii<NUM_SPECIES; ii++) {
    for (int jj = 0; jj<NUM_SPECIES; jj++) {
      d2AmdT2 += tmp1 * Y[ii] * Y[jj] * (-4.0 * T * amlocder[ii] * amlocder[jj] + amloc[ii] * amlocder[jj] + amloc[jj] * amlocder[ii]);
    }
  }
}
  
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
molecular_weight(amrex::Real mw[])
{
  get_mw(mw);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
inv_molecular_weight(amrex::Real imw[])
{
  get_imw(imw);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
TY2Cp(amrex::Real T, amrex::Real Y[], amrex::Real& Cp)
{
  amrex::Error("TY2Cp not physically possible for this EoS");
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
RTY2Cp(amrex::Real R, amrex::Real T, amrex::Real Y[], amrex::Real& Cp)
{
  CKCPBS(&T, Y, &Cp);
}
  
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
TY2Cv(amrex::Real T, amrex::Real Y[], amrex::Real& Cv)
{
  amrex::Error("TY2Cv not physically possible for this EoS");
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
RTY2Cv(amrex::Real R, amrex::Real T, amrex::Real Y[], amrex::Real& Cv)
{
  CKCVBS(&T, Y, &Cv);
}
  
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
T2Cpi(amrex::Real T, amrex::Real Cpi[])
{
  CKCPMS(&T, Cpi);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
RPY2Cs(amrex::Real R, amrex::Real P, amrex::Real Y[], amrex::Real& Cs)
{
  amrex::Real tmp[NUM_SPECIES];
  amrex::Real wbar;
  CKMMWY(&Y[0], &wbar);
  amrex::Real T = P * wbar / (R * EOS::RU);
  CKCVMS(&T, tmp);
  amrex::Real Cv = 0.0;
  for (int i = 0; i < NUM_SPECIES; i++)
    Cv += Y[i] * tmp[i];
  amrex::Real G = (wbar * Cv + EOS::RU) / (wbar * Cv);
  Cs = std::sqrt(G * P / R);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
RTY2Cs(amrex::Real R, amrex::Real T, amrex::Real Y[], amrex::Real& Cs)
{
  amrex::Real tmp[NUM_SPECIES];
  amrex::Real P;
  CKPY(&R, &T, Y, &P);
  CKCVMS(&T, tmp);
  amrex::Real Cv = 0.0;
  for (int i = 0; i < NUM_SPECIES; ++i)
    Cv += Y[i] * tmp[i];
  amrex::Real wbar;
  CKMMWY(&Y[0], &wbar);
  amrex::Real G = (wbar * Cv + EOS::RU) / (wbar * Cv);
  Cs = std::sqrt(G * P / R);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
EY2T(amrex::Real E, amrex::Real Y[], amrex::Real& T)
{
  // For Fuego this function is really just a wrapper for GET_T_GIVEN_EY
  // In SRK this will be different probably
  amrex::Error("EY2T not physically possible for this EoS");
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
REY2T(amrex::Real R, amrex::Real E, amrex::Real Y[], amrex::Real& T)
{
  amrex::Real tau = 1.0/R;
  amrex::Real Tn = T;
  amrex::Real bm, am, dAmdT, d2AmdT2;
  amrex::Real K1, ei_ig[NUM_SPECIES], Eig, Cv, dT;
  int nIter = 0;
  int lierr = 0;
  amrex::Real fzero = 1.0;
  
  // Run some tests
  //  First print inputs
  Tn = 501.0;
  //std::cout << "Testing" <<std::endl;
  //std::cout << "R " << R << std::endl;
  //std::cout << "E " << E << std::endl;
  //std::cout << "Y ";
  for (int ii =0; ii<NUM_SPECIES; ii++) { std::cout << Y[ii]; };
  //std::cout << std::endl;
  //std::cout << "T " << Tn << std::endl;

  // Now run tests
  //std::cout << " Running tests:" << std::endl;
  // EOS::MixingRuleBm(Y, bm);  
  // std::cout << "bm " << bm << std::endl;
  // EOS::MixingRuleAm(Tn, Y, am);
  // std::cout << "am " << am << std::endl;
  // EOS::MixingRuleAmBm(Tn, Y, am, bm);
  // std::cout << "ambm " << am << " " << bm << std::endl;
  // EOS::Calc_dAmdT(Tn, Y, dAmdT);
  // std::cout << "dAmdT " << dAmdT << std::endl;
  // EOS::Calc_d2AmdT2(Tn, Y, d2AmdT2);
  // std::cout << "d2AmdT2" << d2AmdT2 << std::endl;
  
  // amrex::Real temp_vec[NUM_SPECIES];
  // EOS::Calc_dAmdY(Tn, Y, temp_vec);
  // std::cout << "dAmdY ";
  // for (int ii =0; ii<NUM_SPECIES; ii++) { std::cout << temp_vec[ii] << " "; };
  // std::cout << std::endl;
  
  // EOS::Calc_d2AmdTY(Tn, Y, temp_vec);
  // std::cout << "d2AmdTY ";
  // for (int ii =0; ii<NUM_SPECIES; ii++) { std::cout << temp_vec[ii] << " "; };
  // std::cout << std::endl;
  
  // amrex::Real temp_vec2[NUM_SPECIES][NUM_SPECIES];
  // EOS::Calc_d2AmdY2(Tn, Y, temp_vec2);
  // std::cout << "d2AmdTY ";
  // for (int ii =0; ii<NUM_SPECIES; ii++)
  //   {for (int jj =0; jj<NUM_SPECIES; jj++){
  // 	std::cout << temp_vec2[ii][jj] << " ";
  //     }}
  // std::cout << std::endl;
  // std::cout << std::endl;
  
  EOS::MixingRuleBm(Y, bm);
  K1 = (1.0 / bm) * log(1.0 + bm/tau); // log or log10?	\

  // Use input T as initial guess, even though it may be garbage
  
  while ( amrex::Math::abs(fzero) > convCrit && nIter < maxIter) {
    nIter++;
    EOS::MixingRuleAm(Tn, Y, am);
    EOS::Calc_dAmdT(Tn, Y, dAmdT);
    EOS::Calc_d2AmdT2(Tn, Y, d2AmdT2);
    // ideal gas internal energy
    CKUMS(&Tn, ei_ig);
    Eig = 0.0; for (int ii = 0; ii< NUM_SPECIES; ii++) {Eig += Y[ii]*ei_ig[ii];};
    // ideal gas heat capacity
    CKCVBS(&Tn, Y, &Cv);
    //std::cout << "Cv1 " <<Cv << std::endl;
    //std::cout << "Eig " << Eig << std::endl;
    // real gas heat capacity
    Cv -= Tn * d2AmdT2 * K1;

    // Take difference between E(Tn) and E and iterate
    fzero = -E + Eig + (Tn * dAmdT - am)*K1;
    //std::cout << fzero << std::endl;
    //std::cout << d2AmdT2 << std::endl;
    //std::cout << Cv << std::endl;
    dT = fzero / Cv;
    Tn -= dT;
    //std::cout << dT << std::endl;
    //std::cout << "T " << nIter << " " << Tn << std::endl;
  }
  T = Tn;
  std::cout << "T " << nIter << " " << Tn << std::endl;
  //std::cout << "E " << E << " Eig " << Eig << std::endl << std::endl;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
TPY2R(amrex::Real T, amrex::Real P, amrex::Real Y[],amrex::Real& R)
{
  amrex::Real am, bm, Z, wbar, dAmdT, K1;
  amrex::Real ei_ig[NUM_SPECIES];
  EOS::MixingRuleAmBm(T, Y, am, bm);
  CKMMWY(&Y[0], &wbar);
  EOS::Calc_CompressFactor_Z(Z, am, bm, P, T, wbar);
  R = P * wbar/(Z*EOS::RU*T); 
  //std::cout << "Z " << Z << std::endl;
  //std::cout << "R " << R << std::endl; 
  //std::cout << "T " << T << std::endl; 
  //std::cout << "P " << P << std::endl;
  //std::cout << "RU " << EOS::RU << std::endl;
  //std::cout << "wbar " << wbar << std::endl;   
}
  
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
HY2T(amrex::Real H, amrex::Real Y[], amrex::Real& T)
{
  // For Fuego this function is really just a wrapper for GET_T_GIVEN_HY
  // In SRK this will be different probably
  int lierr = 0;
  GET_T_GIVEN_HY(&H, Y, &T, &lierr);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
RYET2P(
  amrex::Real R,
  amrex::Real Y[],
  amrex::Real& E,
  amrex::Real& T,
  amrex::Real& P)
{
  int lierr = 0;
  GET_T_GIVEN_EY(&E, Y, &T, &lierr);
  CKPY(&R, &T, Y, &P);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
RTY2P(amrex::Real R, amrex::Real T, amrex::Real Y[], amrex::Real& P)
{
  CKPY(&R, &T, Y, &P);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
RYP2T(amrex::Real R, amrex::Real Y[], amrex::Real P, amrex::Real& T)
{
  amrex::Real wbar;
  CKMMWY(&Y[0], &wbar);
  T = P * wbar / (R * EOS::RU);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
RTY2WDOT(amrex::Real R, amrex::Real T, amrex::Real Y[], amrex::Real WDOT[])
{
  //CKWYR(&R, &T, Y, WDOT);
  amrex::Real C[NUM_SPECIES];
  CKYTCR(&R, &T, Y, C);
  CKWC(&T, C, WDOT);

  amrex::Real mw[NUM_SPECIES];
  get_mw(mw);
  for (int n = 0; n < NUM_SPECIES; n++)
    WDOT[n] *= mw[n];
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
RTY2JAC(amrex::Real R, amrex::Real T, amrex::Real Y[], amrex::Real JAC[], int HP)
{
  amrex::Real C[NUM_SPECIES];
  CKYTCR(&R, &T, Y, C);
  DWDOT(JAC, C, &T,&HP);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
RTY2C(amrex::Real R, amrex::Real T, amrex::Real Y[], amrex::Real acti[])
{
  CKYTCR(&R, &T,  Y, acti); 
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
T2Ei(amrex::Real T, amrex::Real Ei[])
{
  CKUMS(&T, Ei);
}

// Function added for non-ideal
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
RTY2E(amrex::Real R, amrex::Real T, amrex::Real Y[], amrex::Real& E)
{
  amrex::Real Ei[NUM_SPECIES];
  amrex::Real am, bm, dAmdT, K1;
  
  // Calculate ideal gas portion 
  E = 0.0;
  CKUMS(&T, Ei);
  for (int ii = 0; ii < NUM_SPECIES; ii++){
    E += Ei[ii] * Y[ii];
  }

  // Add in non-ideal portion
  MixingRuleAmBm(T, Y, am, bm);
  Calc_dAmdT(T, Y, dAmdT);
  // below was log(1 + bm/tau), tau = 1/R so this is simpler)
  K1 = (1.0 / bm) * log(1.0 + bm*R);
  E += (T * dAmdT - am)*K1;

  {
    amrex::Real wbar, tau, P;
  CKMMWY(&Y[0], &wbar);
  tau = 1/R;
  P = (EOS::RU/wbar)* T/(tau-bm) - am /(tau*(tau+bm));
  // std::cout << "Recalculated P " << P << std::endl;
  }
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
Y2X(amrex::Real Y[], amrex::Real X[])
{
  CKYTX(Y, X);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
X2Y(amrex::Real X[], amrex::Real Y[])
{
  CKXTY(X, Y);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
PYT2RE(
  amrex::Real P, amrex::Real Y[], amrex::Real T, amrex::Real& R, amrex::Real& E)
{
  CKRHOY(&P, &T, Y, &R);
  amrex::Real ei[NUM_SPECIES];
  CKUMS(&T, ei);
  E = 0.0;
  for (int n = 0; n < NUM_SPECIES; n++)
    E += Y[n] * ei[n];
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
PYT2R(amrex::Real P, amrex::Real Y[], amrex::Real T, amrex::Real &R) 
{
  CKRHOY(&P, &T, Y, &R);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
RYP2E(amrex::Real R, amrex::Real Y[], amrex::Real P, amrex::Real& E)
{
  amrex::Real wbar;
  CKMMWY(&Y[0], &wbar);
  amrex::Real T = P * wbar / (R * EOS::RU);
  amrex::Real ei[NUM_SPECIES];
  CKUMS(&T, ei);
  E = 0.0;
  for (int n = 0; n < NUM_SPECIES; n++)
    E += Y[n] * ei[n];
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
T2Hi(amrex::Real T, amrex::Real Hi[])
{
  CKHMS(&T, Hi);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
S(amrex::Real& s)
{
  s = 1.0;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
TY2G(amrex::Real T, amrex::Real Y[NUM_SPECIES], amrex::Real& G)
{
  amrex::Real wbar, Cv, Cvx;
  EOS::TY2Cv(T, Y, Cv);
  CKMMWY(&Y[0], &wbar);
  Cvx = wbar * Cv;
  G = (Cvx + EOS::RU) / Cvx;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
TY2H(amrex::Real T, amrex::Real Y[NUM_SPECIES], amrex::Real& H)
{
  amrex::Real hi[NUM_SPECIES];
  CKHMS(&T, hi);
  H = 0.0;
  for (int n = 0; n < NUM_SPECIES; n++)
    H += Y[n] * hi[n];
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
Y2WBAR(amrex::Real Y[], amrex::Real& WBAR)
{
  amrex::Real tmp[NUM_SPECIES];
  get_imw(tmp);
  amrex::Real summ = 0.0;
  for (int i = 0; i < NUM_SPECIES; i++)
    summ += Y[i] * tmp[i];
  WBAR = 1.0 / summ;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
RPE2dpdr_e(amrex::Real R, amrex::Real P, amrex::Real E, amrex::Real& dpdr_e)
{
  dpdr_e = P / R;
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
RG2dpde(amrex::Real R, amrex::Real G, amrex::Real& dpde)
{
  dpde = (G - 1.0) * R;
}




  
}; // namespace EOS

#endif
