#include <iostream>
#include <vector>

#include <AMReX_MultiFab.H>
#include <AMReX_Print.H>
#include <AMReX_VisMF.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>

#include <mechanism.H>
#include <PelePhysics.H>
#include <AMReX_GpuDevice.H>

#include <ReactorBase.H>

AMREX_GPU_DEVICE
inline void
initialize_data(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& mf,
  amrex::Array4<amrex::Real> const& temp,
  amrex::Array4<amrex::Real> const& rho,
  const amrex::GpuArray<amrex::Real, 3>& dx,
  const amrex::GpuArray<amrex::Real, 3>& plo,
  const amrex::GpuArray<amrex::Real, 3>& phi) noexcept
{
  amrex::Real dTemp = 5.0;
  amrex::Real dRho = 0.005;
  amrex::Real y = plo[1] + (j + 0.5) * dx[1];
  amrex::Real x = plo[0] + (i + 0.5) * dx[0];
  amrex::Real pi = 3.1415926535897932;
  amrex::Real L[3];
  amrex::Real P[3];
  amrex::Real Y_lo[NUM_SPECIES];
  amrex::Real Y_hi[NUM_SPECIES];

  for (int n = 0; n < 3; n++) {
    L[n] = phi[n] - plo[n];
    P[n] = L[n] / 4.0;
  }
  for (int n = 0; n < NUM_SPECIES; n++) {
    Y_lo[n] = 0.0;
    Y_hi[n] = 1.0 / NUM_SPECIES;
  }
  Y_lo[0] = 1.0;

  // T, Yk, rho
  amrex::Real Tavg = 500;
  amrex::Real Ravg = 0.01;
#if (AMREX_SPACEDIM == 1)
  temp(i, j, k) = Tavg;
  rho(i, j, k) = Ravg;
#else
  temp(i, j, k) = Tavg + dTemp * std::sin(2.0 * pi * y / P[1]);
  rho(i, j, k) = Ravg + dRho * std::sin(2.0 * pi * y / P[1]);
#endif
  for (int n = 0; n < NUM_SPECIES; n++) {
    mf(i, j, k, n) = Y_lo[n] + (Y_hi[n] - Y_lo[n]) * x / L[0];
  }
  // corr Yk
  amrex::Real dummy = 0.0;
  for (int n = 0; n < NUM_SPECIES - 1; n++) {
    dummy = dummy + mf(i, j, k, n);
  }
  mf(i, j, k, NUM_SPECIES - 1) = 1.0 - dummy;
}

int
main(int argc, char* argv[])
{
  amrex::Initialize(argc, argv);
  {
    amrex::Print() << " Initialization of EOS (CPP)... \n";

    amrex::ParmParse pp;
    int nc = 512;
    pp.query("nc", nc);
    std::vector<int> npts(3, nc);

    amrex::Box domain(
      amrex::IntVect(D_DECL(0, 0, 0)),
      amrex::IntVect(D_DECL(npts[0] - 1, npts[1] - 1, npts[2] - 1)));

    amrex::GpuArray<amrex::Real, 3> plo, phi, dx;
    for (int i = 0; i < BL_SPACEDIM; ++i) {
      phi[i] = domain.length(i);
      dx[i] = (phi[i] - plo[i]) / domain.length(i);
    }

    amrex::RealBox real_box(
      {AMREX_D_DECL(plo[0], plo[1], plo[2])},
      {AMREX_D_DECL(phi[0], phi[1], phi[2])});

    int coord = 0;

    amrex::Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1, 1, 1)};

    amrex::Geometry geom(domain, real_box, coord, is_periodic);

    int max_size = 128;
    pp.query("max_size", max_size);
    amrex::BoxArray ba(domain);
    ba.maxSize(max_size);

    const int num_spec = NUM_SPECIES;

    amrex::DistributionMapping dm{ba};

    int num_grow = 0;
    amrex::MultiFab mass_frac(ba, dm, num_spec, num_grow);
    amrex::MultiFab temperature(ba, dm, 1, num_grow);
    amrex::MultiFab density(ba, dm, 1, num_grow);

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
    for (amrex::MFIter mfi(mass_frac, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
      const amrex::Box& gbox = mfi.tilebox();

      amrex::Array4<amrex::Real> const& Y_a = mass_frac.array(mfi);
      amrex::Array4<amrex::Real> const& T_a = temperature.array(mfi);
      amrex::Array4<amrex::Real> const& rho_a = density.array(mfi);

      amrex::ParallelFor(
        gbox, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          initialize_data(i, j, k, Y_a, T_a, rho_a, dx, plo, phi);
        });
    }

    bool write_output = pp.countval("plot_file") > 0;
    amrex::MultiFab wdots(ba, dm, num_spec, num_grow);

#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(mass_frac, amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {
      const amrex::Box& box = mfi.tilebox();

      const auto mf = mass_frac.array(mfi);
      const auto temp = temperature.array(mfi);
      const auto rho = density.array(mfi);
      const auto wdot = wdots.array(mfi);

      ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
        amrex::Real Yl[NUM_SPECIES];
        amrex::Real Wl[NUM_SPECIES];
        for (int n = 0; n < NUM_SPECIES; ++n)
          Yl[n] = mf(i, j, k, n);
        auto eos = pele::physics::PhysicsType::eos();
        eos.RTY2WDOT(rho(i, j, k), temp(i, j, k), Yl, Wl);
        for (int n = 0; n < NUM_SPECIES; ++n)
          wdot(i, j, k, n) = Wl[n];
      });
    }

    if (write_output) {
      std::string pltfile;
      pp.get("plot_file", pltfile);

      std::string outfile = amrex::Concatenate(pltfile, 1);
      amrex::Vector<std::string> plt_VarsName;
      for (int k = 0; k < NUM_SPECIES; ++k) {
        plt_VarsName.push_back("wdot_" + std::to_string(k));
      }
      amrex::WriteSingleLevelPlotfile(
        outfile, wdots, plt_VarsName, geom, 0.0, 0);
    }
  }

  amrex::Finalize();

  return 0;
}
