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
  amrex::Initialize(argc, argv);

#ifdef AMREX_USE_GPU
  amrex::sundials::Initialize();
#endif
  {

    amrex::Real strt_time = amrex::ParallelDescriptor::second();
    BL_PROFILE_VAR("main::main()", pmain);

    // ~~~~ Init: Read input, initialize geom, data
    // Parse the relevant inputs
    std::string fuel_name;
    std::string chem_integrator;
    bool do_plt;
    std::string pltfile;
    std::string outputFolderHR;
    int initFromChk, reactFunc, ode_ncells, ndt, ode_iE, use_typ_vals,
      max_grid_size;
    std::string chkfile, reactFormat;
    amrex::Real dt, rtol, atol;
    amrex::Real t0, equiv_ratio, press;
    std::array<int, 3> ncells;
    amrex::ParmParse pp;
    amrex::ParmParse ppode("ode");
    amrex::ParmParse pphr("hr");
    parseInput(
      pp, ppode, pphr, fuel_name, chem_integrator, do_plt, pltfile, initFromChk,
      chkfile, reactFormat, reactFunc, ode_ncells, dt, ndt, ode_iE, rtol, atol,
      use_typ_vals, ncells, max_grid_size, t0, equiv_ratio, press,
      outputFolderHR);

    // initialize eosparm (only really needed for manifold EOS)
    pele::physics::PeleParams<
      pele::physics::eos::EosParm<pele::physics::PhysicsType::eos_type>>
      eos_parms;
    eos_parms.initialize();
    auto const* leosparm = eos_parms.device_parm();

    // Assign Fuel ID - don't need to do this for manifold
    int fuel_idx = -1;
    if (pele::physics::PhysicsType::eos_type::identifier() != "Manifold") {
      getFuelID(fuel_name, fuel_idx);
    }

    // Initialize reactor object inside OMP region, including tolerances
    BL_PROFILE_VAR("main::reactor_info()", reactInfo);
    std::unique_ptr<pele::physics::reactions::ReactorBase> reactor =
      pele::physics::reactions::ReactorBase::create(chem_integrator);
    reactor->init(ode_iE, ode_ncells);
    reactor->set_eos_parm(
      &(eos_parms.host_parm()), leosparm); // only needed for manifold
    BL_PROFILE_VAR_STOP(reactInfo);

    // Initialize Geometry
    int finest_level = 0;
    amrex::Vector<amrex::Geometry> geoms;
    amrex::Vector<amrex::BoxArray> grids;
    amrex::Vector<amrex::DistributionMapping> dmaps;
    initializeGeom(
      geoms, grids, dmaps, finest_level, ncells, ndt, dt, max_grid_size);

    // Initialize Data
    int num_grow = 0;
    amrex::Vector<amrex::MultiFab> mf(finest_level + 1);
    amrex::Vector<amrex::MultiFab> rY_source_ext(finest_level + 1);
    amrex::Vector<amrex::MultiFab> mfE(finest_level + 1);
    amrex::Vector<amrex::MultiFab> rY_source_energy_ext(finest_level + 1);
    amrex::Vector<amrex::MultiFab> fctCount(finest_level + 1);
    amrex::Vector<amrex::iMultiFab> dummyMask(finest_level + 1);
    initializeData(
      num_grow, mf, rY_source_ext, mfE, rY_source_energy_ext, t0, equiv_ratio,
      press, fctCount, dummyMask, finest_level, geoms, grids, dmaps, fuel_idx,
      ode_iE, leosparm);

    // ~~~~ Reac
    amrex::Print() << " \n STARTING THE ADVANCE \n";

    for (int lev = 0; lev <= finest_level; ++lev) {
      amrex::Real lvl_strt = amrex::ParallelDescriptor::second();
      BL_PROFILE_VAR("Advance_Level" + std::to_string(lev), Advance);
#ifdef AMREX_USE_OMP
      const auto tiling = amrex::MFItInfo().SetDynamic(true);
#pragma omp parallel
#else
      const bool tiling = amrex::TilingIfNotGPU();
#endif
      for (amrex::MFIter mfi(mf[lev], tiling); mfi.isValid(); ++mfi) {

        int omp_thread = 0;
#ifdef AMREX_USE_OMP
        omp_thread = omp_get_thread_num();
#endif
        // Reaction with Array4
        if (reactFunc == 1) {
          doReact_Array4(
            lev, dt, ndt, omp_thread, mfi, mf, rY_source_ext, mfE,
            rY_source_energy_ext, fctCount, dummyMask, reactor, outputFolderHR);
          // Reaction with 1dArray
        } else if (reactFunc == 2) {
          doReact_1dArray(
            lev, dt, ndt, omp_thread, ode_ncells, mfi, mf, rY_source_ext, mfE,
            rY_source_energy_ext, fctCount, dummyMask, reactor, outputFolderHR);
        }
      }
      BL_PROFILE_VAR_STOP(Advance);
      amrex::Real lvl_run_time = amrex::ParallelDescriptor::second() - lvl_strt;
      amrex::ParallelDescriptor::ReduceRealMax(
        lvl_run_time, amrex::ParallelDescriptor::IOProcessorNumber());
      amrex::Print() << "   >> Level " << lev << " advance: " << lvl_run_time
                     << "\n";
    }

    // TODO multilevel max.
    {
      amrex::Vector<double> typ_vals(NUM_SPECIES + 1);
      amrex::Print() << "ode.typ_vals= ";
      for (int i = 0; i < NUM_SPECIES + 1; ++i) {
        amrex::Print() << std::max(1.e-10, mf[0].max(i)) << " ";
      }
      amrex::Print() << std::endl;
    }

    // ~~~~ Plot
    plotResult(do_plt, pltfile, finest_level, mf, geoms);

    // ~~~~ Finalize
    BL_PROFILE_VAR_STOP(pmain);
    amrex::Real run_time = amrex::ParallelDescriptor::second() - strt_time;
    amrex::ParallelDescriptor::ReduceRealMax(
      run_time, amrex::ParallelDescriptor::IOProcessorNumber());
    amrex::Print() << " \n >> React::main() " << run_time << "\n\n";

    eos_parms.deallocate();
  }
#ifdef AMREX_USE_GPU
  amrex::sundials::Finalize();
#endif
  amrex::Finalize();

  return 0;
}
