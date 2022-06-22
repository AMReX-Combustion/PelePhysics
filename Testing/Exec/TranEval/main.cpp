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

    pele::physics::transport::TransportParams<
      pele::physics::PhysicsType::transport_type>
      trans_parms;
    trans_parms.allocate();

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

    const auto geomdata = geom.data();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(mass_frac, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {

      const amrex::Box& bx = mfi.tilebox();

      auto const& Y_a = mass_frac.array(mfi);
      auto const& T_a = temperature.array(mfi);
      auto const& rho_a = density.array(mfi);

      amrex::ParallelFor(
        bx, [Y_a, T_a, rho_a,
             geomdata] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          initialize_data(i, j, k, Y_a, T_a, rho_a, geomdata);
        });
    }

    amrex::MultiFab D(ba, dm, NUM_SPECIES, num_grow);
    amrex::MultiFab mu(ba, dm, 1, num_grow);
    amrex::MultiFab xi(ba, dm, 1, num_grow);
    amrex::MultiFab lam(ba, dm, 1, num_grow);

    // Get the transport data pointer
    auto const* ltransparm = trans_parms.device_trans_parm();

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(mass_frac, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {

      const amrex::Box& gbox = mfi.tilebox();

      amrex::Array4<amrex::Real> const& Y_a = mass_frac.array(mfi);
      amrex::Array4<amrex::Real> const& T_a = temperature.array(mfi);
      amrex::Array4<amrex::Real> const& rho_a = density.array(mfi);
      amrex::Array4<amrex::Real> const& D_a = D.array(mfi);
      amrex::Array4<amrex::Real> const& mu_a = mu.array(mfi);
      amrex::Array4<amrex::Real> const& xi_a = xi.array(mfi);
      amrex::Array4<amrex::Real> const& lam_a = lam.array(mfi);

      amrex::launch(gbox, [=] AMREX_GPU_DEVICE(amrex::Box const& tbx) {
        auto trans = pele::physics::PhysicsType::transport();
        trans.get_transport_coeffs(
          tbx, Y_a, T_a, rho_a, D_a, mu_a, xi_a, lam_a, ltransparm);
      });
    }

    trans_parms.deallocate();

    amrex::MultiFab VarPlt(ba, dm, NUM_SPECIES + 3, num_grow);
    amrex::MultiFab::Copy(VarPlt, D, 0, 0, NUM_SPECIES, num_grow);
    amrex::MultiFab::Copy(VarPlt, mu, 0, NUM_SPECIES, 1, num_grow);
    amrex::MultiFab::Copy(VarPlt, xi, 0, NUM_SPECIES + 1, 1, num_grow);
    amrex::MultiFab::Copy(VarPlt, lam, 0, NUM_SPECIES + 2, 1, num_grow);

    std::string outfile = amrex::Concatenate(pltfile, 1);
    amrex::Vector<std::string> plt_VarsName;
    for (int k = 0; k < NUM_SPECIES; ++k) {
      plt_VarsName.push_back("D" + std::to_string(k));
    }
    plt_VarsName.push_back("mu");
    plt_VarsName.push_back("xi");
    plt_VarsName.push_back("lam");
    amrex::WriteSingleLevelPlotfile(
      outfile, VarPlt, plt_VarsName, geom, 0.0, 0);
  }

  amrex::Finalize();

  return 0;
}
