#include <iostream>
#include <vector>

#include <AMReX_MultiFab.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>

#include "mechanism.H"
#include <GPU_misc.H>

#include <PelePhysics.H>

int
main(int argc, char* argv[])
{
  amrex::Initialize(argc, argv);
  {
    pele::physics::eos::EosParams<pele::physics::PhysicsType::eos_type>
      eos_parms;
    amrex::Print() << " Initialization of EOS (CPP)... \n";
#ifdef USE_MANIFOLD_EOS
    static std::unique_ptr<pele::physics::ManFuncParams> manfunc_par;

    amrex::ParmParse ppm("manifold");
    std::string manifold_model;
    ppm.get("model", manifold_model);
    if (manifold_model == "Table") {
      manfunc_par.reset(new pele::physics::TabFuncParams());
      amrex::Print() << " Initialization of Table (CPP)... \n";
      manfunc_par->initialize();
      eos_parms.allocate(manfunc_par->device_manfunc_data());
    } else if (manifold_model == "NeuralNet") {
      manfunc_par.reset(new pele::physics::NNFuncParams());
      amrex::Print() << " Initialization of Neural Net Func. (CPP)... \n";
      manfunc_par->initialize();
      eos_parms.allocate(manfunc_par->device_manfunc_data());
    } else {
      amrex::Error("Invalid manifold model!");
    }
#else
    eos_parms.allocate();
#endif
    auto const* leosparm = eos_parms.device_eos_parm();

    amrex::ParmParse pp;

    int gridsize = 128;
    pp.query("gridsize", gridsize);

    // Define geometry
    amrex::Array<int, AMREX_SPACEDIM> npts{AMREX_D_DECL(1, 1, 1)};

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      npts[i] = gridsize;
    }

    amrex::Box domain(
      amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
      amrex::IntVect(AMREX_D_DECL(npts[0] - 1, npts[1] - 1, npts[2] - 1)));

    amrex::RealBox real_box(
      {AMREX_D_DECL(-1.0, -1.0, -1.0)}, {AMREX_D_DECL(1.0, 1.0, 1.0)});

    int coord = 0;

    amrex::Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1, 1, 1)};

    amrex::Geometry geom(domain, real_box, coord, is_periodic);

    // Define BoxArray
    int max_size = 32;
    pp.query("max_size", max_size);
    amrex::BoxArray ba(domain);
    ba.maxSize(max_size);

    amrex::ParmParse ppa("amr");
    bool do_plot = true;
    ppa.query("do_plot", do_plot);
    std::string pltfile("plt");
    ppa.query("plot_file", pltfile);

    amrex::DistributionMapping dm{ba};
    int num_grow = 0;

    // Data MFs
    amrex::MultiFab mass_frac(ba, dm, NUM_SPECIES, num_grow);
    amrex::MultiFab temperature(ba, dm, 1, num_grow);
    amrex::MultiFab density(ba, dm, 1, num_grow);
    amrex::MultiFab energy(ba, dm, 1, num_grow);

    const auto geomdata = geom.data();
    {
      BL_PROFILE("Pele::init()");
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (amrex::MFIter mfi(mass_frac, amrex::TilingIfNotGPU()); mfi.isValid();
           ++mfi) {

        const amrex::Box& bx = mfi.tilebox();
        auto const& Y_a = mass_frac.array(mfi);
        auto const& T_a = temperature.array(mfi);
        auto const& rho_a = density.array(mfi);
        auto const& e_a = energy.array(mfi);
        amrex::ParallelFor(
          bx, [Y_a, T_a, rho_a, e_a, geomdata,
               leosparm] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            initialize_data(i, j, k, Y_a, T_a, rho_a, e_a, geomdata, leosparm);
          });
      }
    }
    amrex::MultiFab cp(ba, dm, 1, num_grow);
    amrex::MultiFab wdot(ba, dm, NUM_SPECIES, num_grow);
    {
      BL_PROFILE("Pele::cp()");
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (amrex::MFIter mfi(mass_frac, amrex::TilingIfNotGPU()); mfi.isValid();
           ++mfi) {

        const amrex::Box& box = mfi.tilebox();

        auto const& Y_a = mass_frac.const_array(mfi);
        auto const& T_a = temperature.const_array(mfi);
        auto const& cp_a = cp.array(mfi);
        auto const& rho_a = density.array(mfi);
        amrex::ParallelFor(
          box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            get_cp(i, j, k, Y_a, T_a, rho_a, cp_a, leosparm);
          });
      }
    }
    {
      BL_PROFILE("Pele::cv()");
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (amrex::MFIter mfi(mass_frac, amrex::TilingIfNotGPU()); mfi.isValid();
           ++mfi) {

        const amrex::Box& box = mfi.tilebox();

        auto const& Y_a = mass_frac.const_array(mfi);
        auto const& T_a = temperature.const_array(mfi);
        auto const& wdot_a = wdot.array(mfi);
        auto const& rho_a = density.array(mfi);
        amrex::ParallelFor(
          box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            get_wdot(i, j, k, Y_a, T_a, rho_a, wdot_a, leosparm);
          });
      }
    }

    {
      BL_PROFILE("Pele::getE()");
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
      for (amrex::MFIter mfi(mass_frac, amrex::TilingIfNotGPU()); mfi.isValid();
           ++mfi) {

        const amrex::Box& box = mfi.tilebox();

        auto const& Y_a = mass_frac.const_array(mfi);
        auto const& e_a = energy.const_array(mfi);
        auto const& T_a = temperature.array(mfi);
        auto const& rho_a = density.array(mfi);
        amrex::ParallelFor(
          box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
            get_T_from_EY(i, j, k, Y_a, T_a, rho_a, e_a, leosparm);
          });
      }
    }
    if (do_plot) {
      amrex::MultiFab VarPlt(ba, dm, 4 + 2 * NUM_SPECIES, num_grow);
      amrex::MultiFab::Copy(VarPlt, density, 0, 0, 1, num_grow);
      amrex::MultiFab::Copy(VarPlt, mass_frac, 0, 1, NUM_SPECIES, num_grow);
      amrex::MultiFab::Copy(VarPlt, cp, 0, 1 + NUM_SPECIES, 1, num_grow);
      amrex::MultiFab::Copy(
        VarPlt, wdot, 0, 2 + NUM_SPECIES, NUM_SPECIES, num_grow);
      amrex::MultiFab::Copy(
        VarPlt, temperature, 0, 2 + 2 * NUM_SPECIES, 1, num_grow);
      amrex::MultiFab::Copy(
        VarPlt, energy, 0, 3 + 2 * NUM_SPECIES, 1, num_grow);

      std::string outfile = amrex::Concatenate(pltfile, 1);
      amrex::Vector<std::string> plt_VarsName, spec_VarsName;
      pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(
        spec_VarsName);

      plt_VarsName.push_back("rho");
      for (int n = 0; n < NUM_SPECIES; ++n) {
        plt_VarsName.push_back("Y(" + spec_VarsName[n] + ")");
      }
      plt_VarsName.push_back("cp");
      for (int n = 0; n < NUM_SPECIES; ++n) {
        plt_VarsName.push_back("wdot(" + spec_VarsName[n] + ")");
      }
      plt_VarsName.push_back("temperature");
      plt_VarsName.push_back("energy");

      amrex::WriteSingleLevelPlotfile(
        outfile, VarPlt, plt_VarsName, geom, 0.0, 0);
    }
  }

  amrex::Finalize();

  return 0;
}
