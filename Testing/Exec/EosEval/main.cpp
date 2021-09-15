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

    amrex::ParmParse pp;

    // Define geometry
    amrex::Array<int, AMREX_SPACEDIM> npts{AMREX_D_DECL(1, 1, 1)};

    for (int i = 0; i < AMREX_SPACEDIM; ++i) {
      npts[i] = 128;
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

#ifdef _OPENMP
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
        bx, [Y_a, T_a, rho_a, e_a,
             geomdata] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          initialize_data(i, j, k, Y_a, T_a, rho_a, e_a, geomdata);
        });
    }

    amrex::MultiFab VarPlt(ba, dm, 4, num_grow);
    amrex::MultiFab cp(ba, dm, 1, num_grow);
    amrex::MultiFab cv(ba, dm, 1, num_grow);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(mass_frac, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {

      const amrex::Box& box = mfi.tilebox();

      auto const& Y_a = mass_frac.const_array(mfi);
      auto const& T_a = temperature.const_array(mfi);
      auto const& cp_a = cp.array(mfi);
      amrex::ParallelFor(
        box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          get_cp(i, j, k, Y_a, T_a, cp_a);
        });
    }
    amrex::MultiFab::Copy(VarPlt, cp, 0, 0, 1, num_grow);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(mass_frac, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {

      const amrex::Box& box = mfi.tilebox();

      auto const& Y_a = mass_frac.const_array(mfi);
      auto const& T_a = temperature.const_array(mfi);
      auto const& cv_a = cv.array(mfi);
      amrex::ParallelFor(
        box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          get_cv(i, j, k, Y_a, T_a, cv_a);
        });
    }
    amrex::MultiFab::Copy(VarPlt, cv, 0, 1, 1, num_grow);

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(mass_frac, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {

      const amrex::Box& box = mfi.tilebox();

      auto const& Y_a = mass_frac.const_array(mfi);
      auto const& e_a = energy.const_array(mfi);
      auto const& T_a = temperature.array(mfi);
      amrex::ParallelFor(
        box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          get_T_from_EY(i, j, k, Y_a, T_a, e_a);
        });
    }
    amrex::MultiFab::Copy(VarPlt, temperature, 0, 2, 1, num_grow);
    amrex::MultiFab::Copy(VarPlt, energy, 0, 3, 1, num_grow);

    std::string outfile = amrex::Concatenate(pltfile, 1);
    // TODO: add fct count to this output
    amrex::Vector<std::string> plt_VarsName;
    plt_VarsName.push_back("cp");
    plt_VarsName.push_back("cv");
    plt_VarsName.push_back("temperature");
    plt_VarsName.push_back("energy");
    amrex::WriteSingleLevelPlotfile(
      outfile, VarPlt, plt_VarsName, geom, 0.0, 0);
  }

  amrex::Finalize();

  return 0;
}
