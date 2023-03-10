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

    amrex::Real T;
    T = 900.0;
    amrex::Real P;
    P = 1.4*1013250;
    const amrex::Real tc[5] = {
    log(T), T, T * T, T * T * T, T * T * T * T}; // temperature cache
    const amrex::Real invT = 1.0 / tc[1];

    amrex::Real * x = new amrex::Real [NUM_SPECIES];
    amrex::Real * y = new amrex::Real [NUM_SPECIES];
    amrex::Real * sc = new amrex::Real [NUM_SPECIES];
    amrex::Real * sc_ct = new amrex::Real [NUM_SPECIES];
    amrex::Real * wdot = new amrex::Real [NUM_SPECIES];
    amrex::Real * wdot_ct = new amrex::Real [NUM_SPECIES];
    amrex::Real * sc_qss = new amrex::Real [1];
    amrex::Real * qf = new amrex::Real [NUM_REACTIONS];
    amrex::Real * qf_ct = new amrex::Real [NUM_REACTIONS];
    amrex::Real * qr = new amrex::Real [NUM_REACTIONS];
    amrex::Real * qr_ct = new amrex::Real [NUM_REACTIONS];

    // INIT X
    x[O2_ID] = 0.20153;
    x[POSF11498_ID] = 0.040307;
    x[N2_ID] = 0.75816;
    // GET SC
    CKXTCP(P, T, x, sc);

    //for (int i = 0; i < NUM_SPECIES; ++i) {
    //    // dummy numbers
    //    sc[i] = (i+1)*0.001;
    //}
    
    comp_qfqr(qf, qr, sc, sc_qss, tc, invT);
    productionRate(wdot, sc, T);
 

    // Convert to cantera units
    print<double>(sc,NUM_SPECIES,"sc");
    print<double>(wdot,NUM_SPECIES,"wdot");
    print<double>(qf,NUM_REACTIONS,"qf");
    print<double>(qr,NUM_REACTIONS,"qr");

  return 0;
}
