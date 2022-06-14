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

int
main(int argc, char* argv[])
{
    
    int NUM_QSS_SPECIES = 18;
    int NUM_QSS_REACTIONS = 198;
    int NUM_PERT = 5;
    int NUM_JAC_ENTRIES = (NUM_SPECIES+1)*(NUM_SPECIES+1);
    amrex::Real T;
    T = 1000.0;
    const amrex::Real tc[5] = {
    log(T), T, T * T, T * T * T, T * T * T * T}; // temperature cache
    const amrex::Real invT = 1.0 / tc[1];
    amrex::Real * sc = new amrex::Real [NUM_SPECIES];
    amrex::Real * sc_pert1 = new amrex::Real [NUM_SPECIES];
    amrex::Real * sc_pert2 = new amrex::Real [NUM_SPECIES];
    amrex::Real * wdot_pert1 = new amrex::Real [NUM_SPECIES];
    amrex::Real * wdot_pert2 = new amrex::Real [NUM_SPECIES];
    amrex::Real * wdot = new amrex::Real [NUM_SPECIES];
    amrex::Real * wdot_sym = new amrex::Real [NUM_SPECIES];
    amrex::Real * J_fuego_consP = new amrex::Real [NUM_JAC_ENTRIES];
    amrex::Real * J_fuego_noconsP = new amrex::Real [NUM_JAC_ENTRIES];
    amrex::Real * J_sympy = new amrex::Real [1];
    const amrex::Real pertMag[5] = {
    1e-3, 1e-4, 1e-5, 1e-6, 1e-7};//pert
    amrex::Real * error = new amrex::Real [NUM_PERT];
    amrex::Real * J_approx = new amrex::Real [NUM_PERT];
    int index_sc;
    int index_wdot;
    int index_J;
    const int consP = 1;
    const int noconsP = 1;

    // SET SC 
    srand (time(NULL));
    for (int i = 0; i < NUM_SPECIES; ++i) {
        // dummy numbers
        sc[i] = ((double) rand() / (RAND_MAX)) + 1e-3;//(i+1)*0.02;
    }

    // Compute analytically
    ajac_term_debug(J_sympy, sc, T);
    aJacobian(J_fuego_consP, sc, T, consP);
    aJacobian(J_fuego_noconsP, sc, T, noconsP);

    // Perturb
    index_sc = 0;
    index_wdot = 0;
    index_J = index_sc * (NUM_SPECIES+1) + index_wdot;

    for (int i=0; i < NUM_PERT; ++i){
        perturb_sc(sc_pert1, sc_pert2, sc, pertMag[i], index_sc);
        productionRate(wdot_pert1, sc_pert1, T);
        productionRate(wdot_pert2, sc_pert2, T);
        J_approx[i] = (wdot_pert1[index_wdot] - wdot_pert2[index_wdot])/(2.0*pertMag[i]);
        error[i] = std::abs(J_approx[i]-J_sympy[index_J]);
    } 

    print<double>(error,NUM_PERT,"error");
    print<double>(J_approx,NUM_PERT,"J_approx");
    std::cout << "J_sympy = " << J_sympy[index_J] << "\n";
    std::cout << "J_fuego_consP = " << J_fuego_consP[index_J] << "\n";
    std::cout << "J_fuego_noconsP = " << J_fuego_noconsP[index_J] << "\n";
    
     

  return 0;
}

 
