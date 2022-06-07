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
    amrex::Real * qf_qss = new amrex::Real [NUM_QSS_REACTIONS];
    amrex::Real * qr_qss = new amrex::Real [NUM_QSS_REACTIONS];
    amrex::Real * g_RT = new amrex::Real [NUM_SPECIES];
    amrex::Real * g_RT_qss = new amrex::Real [NUM_QSS_SPECIES];
    for (int i = 0; i < NUM_SPECIES; ++i) {
        // dummy numbers
        sc[i] = (i+1)*0.001;
    }
    comp_sc_qss_debug(sc, sc_qss_sym, T);

    gibbs(g_RT, tc);
    gibbs_qss(g_RT_qss, tc);
    comp_k_f_qss(tc, invT, kf_qss);
    comp_qss_coeff(kf_qss, qf_qss, qr_qss, sc, tc, g_RT, g_RT_qss);
    comp_sc_qss(sc_qss, qf_qss, qr_qss);

    print<double>(sc_qss,NUM_QSS_SPECIES,"sc_qss");
    print<double>(sc_qss_sym,NUM_QSS_SPECIES,"sc_qss_sym");


  return 0;
}

 
