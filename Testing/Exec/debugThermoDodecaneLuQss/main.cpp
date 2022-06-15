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

int
main(int argc, char* argv[])
{
    
    int NUM_QSS_SPECIES = 18;
    int NUM_QSS_REACTIONS = 198;
    amrex::Real T;
    T = 1000.0;
    const amrex::Real tc[5] = {
    log(T), T, T * T, T * T * T, T * T * T * T}; // temperature cache
    const amrex::Real invT = 1.0 / tc[1];
    amrex::Real * sc = new amrex::Real [NUM_SPECIES];

    amrex::Real * sc_qss_sym = new amrex::Real [NUM_QSS_SPECIES];
    amrex::Real * sc_qss = new amrex::Real [NUM_QSS_SPECIES];
    amrex::Real * kf_qss = new amrex::Real [NUM_QSS_REACTIONS];
    amrex::Real * kf_qss_sym = new amrex::Real [NUM_QSS_REACTIONS];
    amrex::Real * qf_qss = new amrex::Real [NUM_QSS_REACTIONS];
    amrex::Real * qr_qss = new amrex::Real [NUM_QSS_REACTIONS];
    amrex::Real * qf_qss_sym = new amrex::Real [NUM_QSS_REACTIONS];
    amrex::Real * qr_qss_sym = new amrex::Real [NUM_QSS_REACTIONS];
    amrex::Real * g_RT = new amrex::Real [NUM_SPECIES];
    amrex::Real * g_RT_sym = new amrex::Real [NUM_SPECIES];
    amrex::Real * g_RT_qss = new amrex::Real [NUM_QSS_SPECIES];
    amrex::Real * g_RT_qss_sym = new amrex::Real [NUM_QSS_SPECIES];
    amrex::Real * h_RT = new amrex::Real [NUM_SPECIES];
    amrex::Real * h_RT_sym = new amrex::Real [NUM_SPECIES];
    amrex::Real * h_RT_qss = new amrex::Real [NUM_QSS_SPECIES];
    amrex::Real * h_RT_qss_sym = new amrex::Real [NUM_QSS_SPECIES];
    for (int i = 0; i < NUM_SPECIES; ++i) {
        // dummy numbers
        sc[i] = (i+1)*0.02;
    }

    gibbs(g_RT, tc);
    gibbs_debug(g_RT_sym, tc);
    gibbs_qss(g_RT_qss, tc);
    gibbs_qss_debug(g_RT_qss_sym, tc);
    speciesEnthalpy(h_RT, tc);
    speciesEnthalpy_debug(h_RT_sym, tc);
    speciesEnthalpy_qss(h_RT_qss, tc);
    speciesEnthalpy_qss_debug(h_RT_qss_sym, tc);
    comp_k_f_qss(tc, invT, kf_qss);
    comp_k_f_qss_debug(tc, invT, kf_qss_sym);
 
    std::string filenameBase="logBase.txt";
    std::string filenameSym="logSym.txt";

    print<double>(filenameBase,g_RT,NUM_SPECIES,"g_RT");
    print<double>(filenameSym,g_RT_sym,NUM_SPECIES,"g_RT");
    print<double>(filenameBase,h_RT,NUM_SPECIES,"h_RT");
    print<double>(filenameSym,h_RT_sym,NUM_SPECIES,"h_RT");

    print<double>(filenameBase,g_RT_qss,NUM_QSS_SPECIES,"g_RT_qss");
    print<double>(filenameSym,g_RT_qss_sym,NUM_QSS_SPECIES,"g_RT_qss");
    print<double>(filenameBase,h_RT_qss,NUM_QSS_SPECIES,"h_RT_qss");
    print<double>(filenameSym,h_RT_qss_sym,NUM_QSS_SPECIES,"h_RT_qss");
    print<double>(filenameBase,kf_qss,NUM_QSS_REACTIONS,"kf_qss");
    print<double>(filenameSym,kf_qss_sym,NUM_QSS_REACTIONS,"kf_qss");


  return 0;
}

 
