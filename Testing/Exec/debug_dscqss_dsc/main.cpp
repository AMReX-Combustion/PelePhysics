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

    amrex::Real * g_RT_qss_pert1 = new amrex::Real [NUM_QSS_SPECIES];
    amrex::Real * g_RT_pert1 = new amrex::Real [NUM_SPECIES];
    amrex::Real * h_RT_qss_pert1 = new amrex::Real [NUM_QSS_SPECIES];
    amrex::Real * h_RT_pert1 = new amrex::Real [NUM_SPECIES];
    amrex::Real * kf_qss_pert1 = new amrex::Real [NUM_QSS_REACTIONS];
    amrex::Real * qf_qss_pert1 = new amrex::Real [NUM_QSS_REACTIONS];
    amrex::Real * qr_qss_pert1 = new amrex::Real [NUM_QSS_REACTIONS];

    amrex::Real * h_RT_qss_pert2 = new amrex::Real [NUM_QSS_SPECIES];
    amrex::Real * h_RT_pert2 = new amrex::Real [NUM_SPECIES];
    amrex::Real * g_RT_qss_pert2 = new amrex::Real [NUM_QSS_SPECIES];
    amrex::Real * g_RT_pert2 = new amrex::Real [NUM_SPECIES];
    amrex::Real * kf_qss_pert2 = new amrex::Real [NUM_QSS_REACTIONS];
    amrex::Real * qf_qss_pert2 = new amrex::Real [NUM_QSS_REACTIONS];
    amrex::Real * qr_qss_pert2 = new amrex::Real [NUM_QSS_REACTIONS];

    const amrex::Real pertMag[5] = {
    1e-3, 1e-4, 1e-5, 1e-6, 1e-7};//pert
    amrex::Real * error = new amrex::Real [NUM_PERT];
    amrex::Real * gradient_approx = new amrex::Real [NUM_PERT];
    amrex::Real * gradient_sympy = new amrex::Real [NUM_JAC_ENTRIES];
    int index_sc;
    int index_sc_qss;
    int consP, noconsP;
    consP = 1; 
    noconsP = 0;

    // SET SC 
    srand (time(NULL));
    for (int i = 0; i < NUM_SPECIES; ++i) {
        // dummy numbers
        sc[i] = ((double) rand() / (RAND_MAX)) + 1e-3;//(i+1)*0.02;
    }

    // Compute dscqss5dsc0
    ajac_term_debug(gradient_sympy, sc, T, consP);
    // Compute scqss by finite diff
    gibbs(g_RT, tc);
    gibbs_qss(g_RT_qss, tc);
    speciesEnthalpy(h_RT, tc);
    speciesEnthalpy_qss(h_RT_qss, tc);
    comp_k_f_qss(tc, invT, kf_qss);

    // Perturb
    index_sc_qss = 5;
    index_sc = 0;

    for (int i=0; i < NUM_PERT; ++i){
        perturb_sc(sc_pert1, sc_pert2, sc, pertMag[i], index_sc);
        comp_qss_coeff(kf_qss, qf_qss_pert1, qr_qss_pert1, sc_pert1, tc, g_RT, g_RT_qss);
        comp_sc_qss(sc_qss_pert1, qf_qss_pert1, qr_qss_pert1);
        comp_qss_coeff(kf_qss, qf_qss_pert2, qr_qss_pert2, sc_pert2, tc, g_RT, g_RT_qss);
        comp_sc_qss(sc_qss_pert2, qf_qss_pert2, qr_qss_pert2);
        gradient_approx[i] = (sc_qss_pert1[index_sc_qss] - sc_qss_pert2[index_sc_qss])/(2.0*pertMag[i]);
        error[i] = std::abs(gradient_approx[i]-gradient_sympy[0])/std::abs(gradient_approx[i]);
    } 

    std::cout << "Computing dscqss" << index_sc_qss << "/dsc"<< index_sc << "\n";
    print<double>(error,NUM_PERT,"Relative error");
    print<double>(gradient_approx,NUM_PERT,"gradient_approx");
    std::cout << "gradient_sympy = " << gradient_sympy[0] << "\n";
     

  return 0;
}

 
