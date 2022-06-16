#include <iostream>
#include <vector>

#include <AMReX_MultiFab.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>

#ifdef AMREX_USE_GPU
#include <AMReX_SUNMemory.H>
#endif

#include "mechanism.H"
#include <GPU_misc.H>

#include <PelePhysics.H>
#include <ReactorBase.H>
#include <random>

// Leanify main script
#include "utils/printFunctions.H"
#include "utils/typeFunction.H"
#include "utils/helperFunctions.H"
#include "utils/initFunctions.H"
#include "utils/plotFunctions.H"
#include "utils/reactFunctions.H"

namespace {
const std::string level_prefix{"Level_"};
}

void
GotoNextLine(std::istream& is)
{
  constexpr std::streamsize bl_ignore_max{100000};
  is.ignore(bl_ignore_max, '\n');
}

void 
perturb_sc(amrex::Real* sc1, amrex::Real* sc2, amrex::Real* sc, amrex::Real perturb, int index)
{
    for (int i = 0; i < NUM_SPECIES; ++i) {
        sc1[i] = sc[i];
        sc2[i] = sc[i];
    }
    sc1[index] = sc1[index] + perturb;
    sc2[index] = sc2[index] - perturb;
}
void 
perturb_T(amrex::Real* T1, amrex::Real* T2, amrex::Real T, amrex::Real perturbT)
{
    *T1 = T + perturbT;
    *T2 = T - perturbT;
}
void 
getTdot(amrex::Real T, amrex::Real P, amrex::Real* sc, amrex::Real* Tdot)
{
    const amrex::Real tc[5] = {
    log(T), T, T * T, T * T * T, T * T * T * T}; // temperature cache

    amrex::Real rho; 
    amrex::Real cp;
    amrex::Real hdot;
    amrex::Real h;
    amrex::Real R;
    R = 8.31446261815324e+07;

    amrex::Real cpml[NUM_SPECIES], dcRdT[NUM_SPECIES], e_RT[NUM_SPECIES], h_RT[NUM_SPECIES], wt[NUM_SPECIES];
    amrex::Real wdot[NUM_SPECIES], x[NUM_SPECIES], y[NUM_SPECIES], cpms[NUM_SPECIES];
    amrex::Real hml[NUM_SPECIES];
   
    // Get mole fractions
    CKCTX(sc, x);
    // Get maxx fractions
    CKXTY(x, y);
    // Get wdot
    productionRate(wdot, sc, T);

    //Compute W
    CKWT(wt); 
    //Compute rho [kg/m3]
    CKRHOC(&P, &T, sc, &rho);
    rho *= 1000;
    //Compute cp [J/(K kg)]
    CKCPMS(&T, cpms);
    for (int k = 0; k < NUM_SPECIES; ++k) {
        cpms[k] *= 1e-4;
    } 
    cp = 0;
    for (int k = 0; k < NUM_SPECIES; ++k) {
        cp += cpms[k] * y[k];
    }
    //Compute h [J/mol] hml [J/mol]
    CKHML(&T, hml);
    for (int k = 0; k < NUM_SPECIES; ++k) {
        hml[k] *= 1e-7;
    } 

    //h = 0;
    //for (int k = 0; k < NUM_SPECIES; ++k) {
    //    h += hml[k] * x[k];
    //}
    //std::cout << "h = " << h << "\n";

    // hdot
    hdot = 0;
    for (int k = 0; k < NUM_SPECIES; ++k) {
        hdot += hml[k] * wdot[k];
    }

    // Tdot
    //std::cout << "cp = " << cp << "\n";
    //std::cout << "rho = " << rho << "\n";
    //std::cout << "hdot = " << hdot << "\n";
    *Tdot = -(1/(rho*cp)) * hdot;
}

int
main(int argc, char* argv[])
{
    
    int NUM_QSS_SPECIES = 18;
    int NUM_QSS_REACTIONS = 198;
    int NUM_PERT = 5;
    int NUM_PERT_T = 5;
    int NUM_JAC_ENTRIES = (NUM_SPECIES+1)*(NUM_SPECIES+1);
    amrex::Real T, P, T_pert1, T_pert2;
    T = 1000.0;
    P = 5.0*1013250.0; 
    const amrex::Real tc[5] = {
    log(T), T, T * T, T * T * T, T * T * T * T}; // temperature cache
    const amrex::Real invT = 1.0 / tc[1];
    amrex::Real * sc = new amrex::Real [NUM_SPECIES];
    amrex::Real * x = new amrex::Real [NUM_SPECIES];
    amrex::Real * sc_pert1 = new amrex::Real [NUM_SPECIES];
    amrex::Real * sc_pert2 = new amrex::Real [NUM_SPECIES];
    amrex::Real * wdot_pert1 = new amrex::Real [NUM_SPECIES];
    amrex::Real * wdot_pert2 = new amrex::Real [NUM_SPECIES];
    amrex::Real Tdot_pert1, Tdot_pert2;
    amrex::Real * wdot = new amrex::Real [NUM_SPECIES];
    amrex::Real * wdot_sym = new amrex::Real [NUM_SPECIES];
    amrex::Real * J_fuego_consP = new amrex::Real [NUM_JAC_ENTRIES];
    amrex::Real * J_fuego_noconsP = new amrex::Real [NUM_JAC_ENTRIES];
    amrex::Real * J_sympy = new amrex::Real [1];
    const amrex::Real pertMag[5] = {
    1e-2, 1e-3, 1e-4, 1e-5, 1e-6};//pert
    const amrex::Real pertMagT[5] = {
    10, 1, 1e-1, 1e-2, 1e-3};//pert
    amrex::Real * error_consP = new amrex::Real [NUM_PERT];
    amrex::Real * error_consP_Tdot = new amrex::Real [NUM_PERT];
    amrex::Real * error_consP_Tdot_T = new amrex::Real [NUM_PERT_T];
    amrex::Real * error_noconsP = new amrex::Real [NUM_PERT];
    amrex::Real * error_noconsP_Tdot = new amrex::Real [NUM_PERT];
    amrex::Real * error_noconsP_Tdot_T = new amrex::Real [NUM_PERT_T];
    amrex::Real * J_approx = new amrex::Real [NUM_PERT];
    amrex::Real * J_approx_Tdot = new amrex::Real [NUM_PERT];
    amrex::Real * J_approx_Tdot_T = new amrex::Real [NUM_PERT_T];
    int index_sc;
    int index_wdot;
    int index_J;
    int index_J_Tdot;
    int index_J_Tdot_T;
    const int consP = 1;
    const int noconsP = 0;

    // SET SC 
    srand (time(NULL));
    for (int i = 0; i < NUM_SPECIES; ++i) {
        // dummy numbers
        //sc[i] = ((double) rand() / (RAND_MAX)) + 10;//(i+1)*0.02;
        //sc[i] = (i+1)*0.02;
        sc[i] = 0;
    }
    // fuel
    sc[0] = 0.1; 
    // o2
    sc[8] = 0.1;

    // Compute analytically
    //ajac_term_debug(J_sympy, sc, T);
    aJacobian(J_fuego_consP, sc, T, consP);
    aJacobian(J_fuego_noconsP, sc, T, noconsP);

    // Perturb
    index_sc = 0;
    index_wdot = 0;
    index_J = index_sc * (NUM_SPECIES+1) + index_wdot;

    std::cout << "Computing dwdot" << index_wdot << "/dsc"<< index_sc << "\n";    
    for (int i=0; i < NUM_PERT; ++i){
        perturb_sc(sc_pert1, sc_pert2, sc, pertMag[i], index_sc);
        productionRate(wdot_pert1, sc_pert1, T);
        productionRate(wdot_pert2, sc_pert2, T);
        J_approx[i] = (wdot_pert1[index_wdot] - wdot_pert2[index_wdot])/(2.0*pertMag[i]);
        error_consP[i] = std::abs(J_approx[i]-J_fuego_consP[index_J]);
        error_noconsP[i] = std::abs(J_approx[i]-J_fuego_noconsP[index_J]);
    } 
   
    std::cout << "Computing dTdot/dsc"<< index_sc << "\n";    
    index_J_Tdot = (NUM_SPECIES) * (NUM_SPECIES+1) + index_sc;
    for (int i=0; i < NUM_PERT; ++i){
        perturb_sc(sc_pert1, sc_pert2, sc, pertMag[i], index_sc);
        getTdot(T, P, sc_pert1, &Tdot_pert1);
        getTdot(T, P, sc_pert2, &Tdot_pert2);
        //std::cout << "Tdot_pert1 = " << Tdot_pert1 << "\n";
        //std::cout << "Tdot_pert2 = " << Tdot_pert2 << "\n";
        J_approx_Tdot[i] = (Tdot_pert1 - Tdot_pert2)/(2.0*pertMag[i]);
        error_consP_Tdot[i] = std::abs(J_approx_Tdot[i]-J_fuego_consP[index_J_Tdot]);
        error_noconsP_Tdot[i] = std::abs(J_approx_Tdot[i]-J_fuego_noconsP[index_J_Tdot]);
    } 
    std::cout << "Computing dTdot/dT"<< "\n";    
    index_J_Tdot_T = (NUM_SPECIES+1) * (NUM_SPECIES+1) - 1;
    for (int i=0; i < NUM_PERT_T; ++i){
        perturb_T(&T_pert1, &T_pert2, T, pertMagT[i]);
        getTdot(T_pert1, P, sc, &Tdot_pert1);
        getTdot(T_pert2, P, sc, &Tdot_pert2);
        //std::cout << "Tdot_pert1 = " << Tdot_pert1 << "\n";
        //std::cout << "Tdot_pert2 = " << Tdot_pert2 << "\n";
        J_approx_Tdot_T[i] = (Tdot_pert1 - Tdot_pert2)/(2.0*pertMagT[i]);
        error_consP_Tdot_T[i] = std::abs(J_approx_Tdot_T[i]-J_fuego_consP[index_J_Tdot_T]);
        error_noconsP_Tdot_T[i] = std::abs(J_approx_Tdot_T[i]-J_fuego_noconsP[index_J_Tdot_T]);
    } 
   

    //print<double>(error_consP,NUM_PERT,"error_consP");
    //print<double>(error_noconsP,NUM_PERT,"error_noconsP");
    std::cout << "Printing dwdot" << index_wdot << "/dsc"<< index_sc << "\n";    
    print<double>(J_approx,NUM_PERT,"J_approx");
    std::cout << "J_fuego_consP = " << J_fuego_consP[index_J] << "\n";
    std::cout << "J_fuego_noconsP = " << J_fuego_noconsP[index_J] << "\n";
    
    //print<double>(error_consP_Tdot,NUM_PERT,"error_consP_Tdot");
    //print<double>(error_noconsP_Tdot,NUM_PERT,"error_noconsP_Tdot");
    std::cout << "Printing dTdot/dsc"<< index_sc << "\n";    
    print<double>(J_approx_Tdot,NUM_PERT,"J_approx_Tdot");
    std::cout << "J_fuego_consP_Tdot = " << J_fuego_consP[index_J_Tdot] << "\n";
    std::cout << "J_fuego_noconsP_Tdot = " << J_fuego_noconsP[index_J_Tdot] << "\n";
     
    std::cout << "Printing dTdot/dT" << "\n";    
    print<double>(J_approx_Tdot_T,NUM_PERT_T,"J_approx_Tdot");
    std::cout << "J_fuego_consP_Tdot_T = " << J_fuego_consP[index_J_Tdot_T] << "\n";
    std::cout << "J_fuego_noconsP_Tdot_T = " << J_fuego_noconsP[index_J_Tdot_T] << "\n";
 
    
    //amrex::Real * allJfuegoTdot = new amrex::Real [NUM_SPECIES+1];
    //for (int k=0; k<NUM_SPECIES+1; ++k) {
    //   allJfuegoTdot[k] = J_fuego_noconsP[(NUM_SPECIES) * (NUM_SPECIES+1) + k];
    //}
    //print<double>(allJfuegoTdot,NUM_SPECIES+1,"all fuego Tdot");

   

  return 0;
}

 
