#include <iostream>
#include <vector>

#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Print.H>
#include <AMReX_VisMF.H>
#include <AMReX_iMultiFab.H>

#ifdef AMREX_USE_GPU
#include <AMReX_SUNMemory.H>
#endif

#include "mechanism.H"
#include <PMFData.H>
#include <data_K.H>

#include <PelePhysics.H>
#include <ReactorBase.H>

using namespace amrex;

int main(int argc, char *argv[]) {
  Initialize(argc, argv);

  {
    Real strt_time = amrex::ParallelDescriptor::second();

    BL_PROFILE_VAR("main::main()", pmain);

    // Init PMF
    pele::physics::PMF::PmfData pmf_data;
    pmf_data.initialize();

    // Get standoff
    ParmParse ppp("pele");
    Real standoff = 0.0;
    ppp.query("standoff", standoff);

    // Define the geometry
    ParmParse pp("amr");

    Vector<int> n_cell(AMREX_SPACEDIM);
    pp.getarr("n_cell", n_cell, 0, AMREX_SPACEDIM);

    IntVect lo(IntVect::TheZeroVector()), hi(n_cell);
    hi -= IntVect::TheUnitVector();
    Box domain(lo, hi);

    amrex::RealBox *real_box = nullptr;
    int coord = -1;
    int *is_periodic = nullptr;

    Geometry geom(domain, real_box, coord, &is_periodic[0]);

    // Define BoxArray / Dmap
    int max_grid_size = 16;
    pp.query("max_grid_size", max_grid_size);
    auto grids = BoxArray(domain);
    grids.maxSize(max_grid_size);
    auto dmaps = DistributionMapping(grids, ParallelDescriptor::NProcs());

    // Create MF
    int num_grow = 0;
    //               vel              rho rhoYs        e/h  T
    const int NVAR = AMREX_SPACEDIM + 1 + NUM_SPECIES + 1 + 1;
    MultiFab stateMF(grids, dmaps, NVAR, num_grow);

    FabArrayBase::mfiter_tile_size =
      IntVect(AMREX_D_DECL(1024, 1024, 1024));

    // Initialize data from PMF
    const auto geomdata = geom.data();
    pele::physics::PMF::PmfData::DataContainer const* lpmfdata = pmf_data.getDeviceData();
    auto const& sma = stateMF.arrays();
    amrex::ParallelFor(stateMF,
    [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
    {    
       initdata(i, j, k, sma[box_no], standoff, geomdata, lpmfdata);
    });

    // Print data 
    std::string outfile = "pltInitData";
    Vector<int> isteps(1, 0);
    Vector<IntVect> refRatios(1, {AMREX_D_DECL(2, 2, 2)});
    Vector<std::string> plt_VarsName;
    AMREX_D_TERM(plt_VarsName.push_back("velocity_x");,
                 plt_VarsName.push_back("velocity_y");,
                 plt_VarsName.push_back("velocity_z"));
    plt_VarsName.push_back("density");
    Vector<std::string> specNames;
    pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(specNames);
    for (int n = 0; n < NUM_SPECIES; ++n) {
      plt_VarsName.push_back("rhoY(" + specNames[n] + ")");
    }
    plt_VarsName.push_back("RhoH");
    plt_VarsName.push_back("Temp");

    const MultiFab *mf = &stateMF;
    WriteMultiLevelPlotfile(
      outfile, 1, {mf}, plt_VarsName, {geom},
      0.0, isteps, refRatios);

    BL_PROFILE_VAR_STOP(pmain);

    Real run_time = ParallelDescriptor::second() - strt_time;
    ParallelDescriptor::ReduceRealMax(
      run_time, ParallelDescriptor::IOProcessorNumber());
    Print() << " \n >> PMFReader::main() " << run_time << "\n\n";
  }
  Finalize();

  return 0;
}
