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
    amrex::Real * sc_qss_pert1 = new amrex::Real [NUM_QSS_SPECIES];
    amrex::Real * sc_qss_pert2 = new amrex::Real [NUM_QSS_SPECIES];

    amrex::Real * h_RT_qss = new amrex::Real [NUM_QSS_SPECIES];
    amrex::Real * h_RT = new amrex::Real [NUM_SPECIES];
    amrex::Real * g_RT_qss = new amrex::Real [NUM_QSS_SPECIES];
    amrex::Real * g_RT = new amrex::Real [NUM_SPECIES];
    amrex::Real * kf_qss = new amrex::Real [NUM_QSS_REACTIONS];
    amrex::Real * qf_qss = new amrex::Real [NUM_QSS_REACTIONS];
    amrex::Real * qr_qss = new amrex::Real [NUM_QSS_REACTIONS];

    amrex::Real * qf_qss_pert1 = new amrex::Real [NUM_QSS_REACTIONS];
    amrex::Real * qr_qss_pert1 = new amrex::Real [NUM_QSS_REACTIONS];

    amrex::Real * qf_qss_pert2 = new amrex::Real [NUM_QSS_REACTIONS];
    amrex::Real * qr_qss_pert2 = new amrex::Real [NUM_QSS_REACTIONS];

    amrex::Real pertMag = 1e-5;//pert
    amrex::Real * error = new amrex::Real [NUM_QSS_SPECIES*NUM_SPECIES];
    amrex::Real * gradient_sympy = new amrex::Real [NUM_QSS_SPECIES*NUM_SPECIES];
    amrex::Real * gradient_FD = new amrex::Real [NUM_QSS_SPECIES*NUM_SPECIES];
    int index_sc;
    int index_sc_qss;
    int index_gradient;
    int consP, noconsP;
    consP = 1; 
    noconsP = 0;

    // SET SC 
    srand (time(NULL));
    srand (42);
    for (int i = 0; i < NUM_SPECIES; ++i) {
        // dummy numbers
        sc[i] = ((double) rand() / (RAND_MAX)) + 1e-3;//(i+1)*0.02;
    }

    // Compute dscqss/dsc with sympy
    dscqss_dsc_fast_debug(gradient_sympy, sc, T);

    // Set up thermo and coeffs
    gibbs(g_RT, tc);
    gibbs_qss(g_RT_qss, tc);
    speciesEnthalpy(h_RT, tc);
    speciesEnthalpy_qss(h_RT_qss, tc);
    comp_k_f_qss(tc, invT, kf_qss);

    // Compute dscqss/dsc by finite difference
    for (int i=0; i < NUM_QSS_SPECIES; ++i){
        for (int j=0; j < NUM_SPECIES; ++j){
            index_sc_qss = i;
            index_sc = j;
            index_gradient = index_sc + index_sc_qss * NUM_SPECIES;
            perturb_sc(sc_pert1, sc_pert2, sc, pertMag, index_sc);
            comp_qss_coeff(kf_qss, qf_qss_pert1, qr_qss_pert1, sc_pert1, tc, g_RT, g_RT_qss);
            comp_sc_qss(sc_qss_pert1, qf_qss_pert1, qr_qss_pert1);
            comp_qss_coeff(kf_qss, qf_qss_pert2, qr_qss_pert2, sc_pert2, tc, g_RT, g_RT_qss);
            comp_sc_qss(sc_qss_pert2, qf_qss_pert2, qr_qss_pert2);
            gradient_FD[index_gradient] = (sc_qss_pert1[index_sc_qss] - sc_qss_pert2[index_sc_qss])/(sc_pert1[index_sc] - sc_pert2[index_sc]);
            error[index_gradient] = std::abs(gradient_FD[index_gradient] - gradient_sympy[index_gradient]) / std::max(std::abs(gradient_FD[index_gradient]), 1e-16);
      }
    } 

    
    std::cout << "error" << " : "  << "\n";
    for (int i=0; i < NUM_QSS_SPECIES; ++i){
        for (int j=0; j < NUM_SPECIES; ++j){
                index_sc_qss = i;
                index_sc = j;
                index_gradient = index_sc + index_sc_qss * NUM_SPECIES;
                if (error[index_gradient] > 1e-5 and (std::abs(gradient_sympy[index_gradient])>1e-15 or std::abs(gradient_sympy[index_gradient])>1e-15) and error[index_gradient]< 0.99){
                    std::cout << "\t scqss " << index_sc_qss << "\t sc " << index_sc  << " : " <<  error[index_gradient] << "\n";
                    std::cout << "\t \t FD : " << gradient_FD[index_gradient] << "\n";
                    std::cout << "\t \t SYMPY : " <<  gradient_sympy[index_gradient]  << "\n";
                    std::cout << "\t \t index : " <<  index_gradient << "\n";
                } 
        } 
    }


  return 0;
}

 
