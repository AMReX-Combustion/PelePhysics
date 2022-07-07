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
perturb_T(amrex::Real * T1, amrex::Real * T2, amrex::Real T, amrex::Real perturb)
{
    *T1 = T + perturb;
    *T2 = T - perturb;
}

int
main(int argc, char* argv[])
{
    
    int NUM_QSS_SPECIES = 18;
    int NUM_QSS_REACTIONS = 198;
    int NUM_PERT = 5;
    int NUM_JAC_ENTRIES = (NUM_SPECIES+1)*(NUM_SPECIES+1);
    amrex::Real T, T_pert1, T_pert2;
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
    amrex::Real * J_sympy = new amrex::Real [NUM_JAC_ENTRIES];
    amrex::Real * J_FD = new amrex::Real [NUM_JAC_ENTRIES];
    amrex::Real pertMag = 1e-5;
    amrex::Real pertMagT = 1e-2;
    amrex::Real * error = new amrex::Real [NUM_JAC_ENTRIES];
    int index_sc;
    int index_T;
    int index_wdot;
    int index_J, index_J_T;
    const int consP = 1;
    const int noconsP = 0;

    // SET SC 
    //srand (time(NULL));
    srand (42);
    for (int i = 0; i < NUM_SPECIES; ++i) {
        // dummy numbers
        sc[i] = ((double) rand() / (RAND_MAX)) + 1e-3;//(i+1)*0.02;
    }

    // Compute analytically
    aJacobian(J_sympy, sc, T, consP);

    // Compute dwdot/dsc by finite difference
    for (int i=0; i < NUM_SPECIES; ++i){
        for (int j=0; j < NUM_SPECIES; ++j){

            index_sc = i;
            index_wdot = j;
            index_J = index_sc * (NUM_SPECIES+1) + index_wdot;


            perturb_sc(sc_pert1, sc_pert2, sc, pertMag, index_sc);
            productionRate(wdot_pert1, sc_pert1, T);
            productionRate(wdot_pert2, sc_pert2, T);
            J_FD[index_J] = (wdot_pert1[index_wdot] - wdot_pert2[index_wdot])/(2.0*pertMag);
            error[index_J] = std::abs(J_FD[index_J]-J_sympy[index_J])/std::max(std::abs(J_FD[index_J]),1e-16);
        }
    } 
    // Compute dwdot/dT by finite difference
    for (int j=0; j < NUM_SPECIES; ++j){

        index_T = NUM_SPECIES;
        index_wdot = j;
        index_J = index_T * (NUM_SPECIES+1) + index_wdot;
        perturb_T(&T_pert1, &T_pert2, T, pertMagT);
        productionRate(wdot_pert1, sc, T_pert1);
        productionRate(wdot_pert2, sc, T_pert2);
        J_FD[index_J] = (wdot_pert1[index_wdot] - wdot_pert2[index_wdot])/(2.0*pertMagT);
        error[index_J] = std::abs(J_FD[index_J]-J_sympy[index_J])/std::max(std::abs(J_FD[index_J]),1e-16);
    }

    std::cout << "error" << " : "  << "\n";
    for (int i=0; i < NUM_SPECIES; ++i){
        for (int j=0; j < NUM_SPECIES; ++j){
            index_sc = i;
            index_wdot = j;
            index_J = index_sc * (NUM_SPECIES+1) + index_wdot;
            //if (error[index_J] > 1e-5 and (std::abs(J_sympy[index_J])>1e-15 or std::abs(J_sympy[index_J])>1e-15) and error[index_J]< 0.99){
            if (error[index_J] > 1e-3 and (std::abs(J_sympy[index_J])>1e-15 or std::abs(J_sympy[index_J])>1e-15) ){
                std::cout << "\t wdot " << index_wdot << "\t sc " << index_sc  << " : " <<  error[index_J] << "\n";
                std::cout << "\t \t FD : " << J_FD[index_J] << "\n";
                std::cout << "\t \t SYMPY : " <<  J_sympy[index_J]  << "\n";
            }
        }
    }
    std::cout << "errorT" << " : "  << "\n";
    for (int j=0; j < NUM_SPECIES; ++j){
        index_T = NUM_SPECIES;
        index_wdot = j;
        index_J = index_T * (NUM_SPECIES+1) + index_wdot;
        if (error[index_J] > 1e-2 and (std::abs(J_sympy[index_J])>1e-15 or std::abs(J_sympy[index_J])>1e-15) ){
            std::cout << "\t wdot " << index_wdot  << " : " <<  error[index_J] << "\n";
            std::cout << "\t \t FD : " << J_FD[index_J] << "\n";
            std::cout << "\t \t SYMPY : " <<  J_sympy[index_J]  << "\n";
        }
    }

  return 0;
}

 
