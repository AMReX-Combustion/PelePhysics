#include <AMReX.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <POneSingleEB.H>

#include <AMReX_EB2.H>
#include <AMReX_EB2_IF.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_EBFabFactory.H>

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
actual_init_coefs_eb(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& phi,
  amrex::Array4<amrex::Real> const& phi_exact,
  amrex::Array4<amrex::Real> const& rhs,
  amrex::Array4<amrex::Real> const& acoef,
  amrex::Array4<amrex::Real> const& bx,
  amrex::Array4<amrex::EBCellFlag const> const& flag,
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
  amrex::Box const& vbx)
{
  amrex::Real const L = std::sqrt(2.0);
  amrex::Real const n = 3.0;
  amrex::Real npioverL = n * M_PI / L;

  amrex::Real pioverfour = M_PI / 4.0;
  amrex::Real cospioverfour = std::cos(pioverfour);
  amrex::Real sinpioverfour = std::sin(pioverfour);

  if (vbx.contains(i, j, k)) {
    bx(i, j, k) = 1.0;
    phi(i, j, k) = 0.0;

    if (flag(i, j, k).isCovered()) {
      rhs(i, j, k) = 0.0;
      phi_exact(i, j, k) = 0.0;
    } else {
      amrex::Real x = dx[0] * (i + 0.5) - 0.5;
      amrex::Real y = dx[1] * (j + 0.5) - 0.5;
      amrex::Real z = dx[2] * k;

      // rotation
      amrex::Real xp = x * cospioverfour + y * sinpioverfour;
      amrex::Real yp = -x * sinpioverfour + y * cospioverfour;

      amrex::Real sincossin = std::sin(npioverL * xp) *
                              std::cos(npioverL * yp) * std::sin(npioverL * z);

      rhs(i, j, k) = (1.0 + npioverL * npioverL) * bx(i, j, k) * sincossin;
      acoef(i, j, k) = 1.0;

      phi_exact(i, j, k) = sincossin;
    }
  }
}

void
initProbABecLaplacian(
  amrex::Geometry& geom,
  std::unique_ptr<amrex::EBFArrayBoxFactory>& factory,
  amrex::MultiFab& solution,
  amrex::MultiFab& rhs,
  amrex::MultiFab& exact_solution,
  amrex::MultiFab& acoef,
  amrex::MultiFab& bcoef)
{

  amrex::FabArray<amrex::EBCellFlagFab> const& flags =
    factory->getMultiEBCellFlagFab();
  //    amrex::MultiCutFab const& bcent = factory->getBndryCent();
  //    amrex::MultiCutFab const& cent  = factory->getCentroid();

  auto const dx = geom.CellSizeArray();
  amrex::Box const& domainbox = geom.Domain();

  for (amrex::MFIter mfi(rhs); mfi.isValid(); ++mfi) {
    amrex::Box const& bx = mfi.validbox();
    amrex::Box const& nbx = amrex::surroundingNodes(bx);
    //        amrex::Box const& gbx = amrex::grow(bx, 1);
    amrex::Array4<amrex::Real> const& phi_arr = solution.array(mfi);
    amrex::Array4<amrex::Real> const& phi_ex_arr = exact_solution.array(mfi);
    amrex::Array4<amrex::Real> const& rhs_arr = rhs.array(mfi);

    amrex::Array4<amrex::Real> const& bx_arr = bcoef.array(mfi);

    auto fabtyp = flags[mfi].getType(bx);
    if (fabtyp == amrex::FabType::covered) {
      std::cout << " amrex::FabType::covered == fabtyp \n";
    } else if (fabtyp == amrex::FabType::regular) {
      std::cout << " amrex::FabType::regular == fabtyp \n";
    } else {
      //            std::cout << " amrex::FabType  else \n";
      amrex::Array4<amrex::Real> const& acoef_arr = acoef.array(mfi);
      amrex::Array4<amrex::EBCellFlag const> const& flag_arr =
        flags.const_array(mfi);
      // amrex::Array4<amrex::Real const> const& cent_arr
      //     = cent.const_array(mfi);
      // amrex::Array4<amrex::Real const> const& bcent_arr
      //     = bcent.const_array(mfi);

      amrex::ParallelFor(
        nbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          actual_init_coefs_eb(
            i, j, k, phi_arr, phi_ex_arr, rhs_arr, acoef_arr, bx_arr, flag_arr,
            dx, bx);
        });
    }
  }
}

amrex::Real
check_norm(
  amrex::MultiFab const& phi,
  amrex::MultiFab const& exact,
  amrex::MultiFab const& vfrc)
{
  amrex::MultiFab mf(phi.boxArray(), phi.DistributionMap(), 1, 0);
  amrex::MultiFab::Copy(mf, phi, 0, 0, 1, 0);

  amrex::MultiFab::Subtract(mf, exact, 0, 0, 1, 0);
  amrex::MultiFab::Multiply(mf, vfrc, 0, 0, 1, 0);

  amrex::Real L0norm = mf.norm0();
  amrex::Real L1norm = mf.norm1();
  std::cout << " L0 norm:" << L0norm << ", L1 norm:" << L1norm << std::endl;

  return L1norm;
}

int
main(int argc, char* argv[])
{
  amrex::Initialize(argc, argv);
  amrex::ParmParse pp;
  PeleRad::AMRParam amrpp(pp);
  PeleRad::MLMGParam mlmgpp(pp);

  bool const write = false;
  int const n_cell = amrpp.n_cell_;

  amrex::Geometry geom;
  amrex::BoxArray grids;
  amrex::DistributionMapping dmap;

  amrex::MultiFab solution;
  amrex::MultiFab exact_solution;
  amrex::MultiFab rhs;
  amrex::MultiFab acoef;
  amrex::MultiFab bcoef;

  amrex::MultiFab robin_a;
  amrex::MultiFab robin_b;
  amrex::MultiFab robin_f;

  std::unique_ptr<amrex::EBFArrayBoxFactory> factory;

  int const max_grid_size = amrpp.max_grid_size_;

  amrex::RealBox rb(
    {AMREX_D_DECL(-1.0, -1.0, -1.0)}, {AMREX_D_DECL(1.0, 1.0, 1.0)});
  amrex::Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1, 1, 1)};
  amrex::Geometry::Setup(&rb, 0, is_periodic.data());
  amrex::Box domain0(
    amrex::IntVect{AMREX_D_DECL(0, 0, 0)},
    amrex::IntVect{AMREX_D_DECL(n_cell - 1, n_cell - 1, n_cell - 1)});
  //    geom.define(domain0);
  geom.define(domain0, rb, amrex::CoordSys::cartesian, is_periodic);

  grids.define(domain0);
  grids.maxSize(max_grid_size);

  // rotated box
  int const max_coarsening_level = mlmgpp.max_coarsening_level_;
  amrex::Real const la = std::sqrt(2.0) / 2.0;
  amrex::Real const shift = std::sqrt(2.0) / 3.0;
  amrex::EB2::BoxIF box(
    {AMREX_D_DECL(-la, -la, -1.0)}, {AMREX_D_DECL(la, la, -1 + 4.0 * shift)},
    true);
  auto gshop = amrex::EB2::makeShop(amrex::EB2::translate(
    amrex::EB2::rotate(
      amrex::EB2::translate(box, {AMREX_D_DECL(-0.0, -0.0, -0.0)}),
      std::atan(1.0) * 1.0, 2),
    {AMREX_D_DECL(0.0, 0.0, 0.0)}));
  amrex::EB2::Build(gshop, geom, 0, max_coarsening_level);

  amrex::IntVect ng = amrex::IntVect{1};

  dmap.define(grids);

  amrex::EB2::IndexSpace const& eb_is = amrex::EB2::IndexSpace::top();
  amrex::EB2::Level const& eb_level = eb_is.getLevel(geom);
  factory = std::make_unique<amrex::EBFArrayBoxFactory>(
    eb_level, geom, grids, dmap, amrex::Vector<int>{2, 2, 2},
    amrex::EBSupport::full);

  amrex::MultiFab const& vfrc = factory->getVolFrac();

  solution.define(grids, dmap, 1, 1, amrex::MFInfo(), *factory);
  exact_solution.define(grids, dmap, 1, 0, amrex::MFInfo(), *factory);
  rhs.define(grids, dmap, 1, 0, amrex::MFInfo(), *factory);
  acoef.define(grids, dmap, 1, 0, amrex::MFInfo(), *factory);
  bcoef.define(grids, dmap, 1, 1, amrex::MFInfo(), *factory);
  robin_a.define(grids, dmap, 1, 0, amrex::MFInfo(), *factory);
  robin_b.define(grids, dmap, 1, 0, amrex::MFInfo(), *factory);
  robin_f.define(grids, dmap, 1, 0, amrex::MFInfo(), *factory);

  solution.setVal(0.0);
  rhs.setVal(0.0);
  acoef.setVal(1.0);
  bcoef.setVal(1.0);
  robin_a.setVal(0.0);
  robin_b.setVal(0.0);
  robin_f.setVal(0.0);

  initProbABecLaplacian(
    geom, factory, solution, rhs, exact_solution, acoef, bcoef);

  //     std::cout << "initialize data ... \n";
  //    initMeshandData(amrpp, mlmgpp, geom, grids, dmap, factory, solution,
  //    rhs,
  //        exact_solution, acoef, bcoef);

  std::cout << "construct the PDE ... \n";
  PeleRad::POneSingleEB rte(
    mlmgpp, geom, grids, dmap, factory, solution, rhs, acoef, bcoef, robin_a,
    robin_b, robin_f);

  // std::cout << "solve the PDE ... \n";
  rte.solve();

  auto eps = check_norm(solution, exact_solution, vfrc);
  eps /= static_cast<amrex::Real>(n_cell * n_cell * n_cell);
  std::cout << "n_cell=" << n_cell << ", normalized L1 norm:" << eps
            << std::endl;

  // plot results
  if (write) {
    std::cout << "write the results ... \n";
    amrex::MultiFab plotmf(grids, dmap, 6, 0);
    amrex::MultiFab::Copy(plotmf, solution, 0, 0, 1, 0);
    amrex::MultiFab::Copy(plotmf, exact_solution, 0, 1, 1, 0);
    amrex::MultiFab::Copy(plotmf, rhs, 0, 2, 1, 0);
    amrex::MultiFab::Copy(plotmf, acoef, 0, 3, 1, 0);
    amrex::MultiFab::Copy(plotmf, bcoef, 0, 4, 1, 0);
    amrex::MultiFab::Copy(plotmf, vfrc, 0, 5, 1, 0);

    auto const plot_file_name = amrpp.plot_file_name_;
    amrex::WriteSingleLevelPlotfile(
      plot_file_name, plotmf,
      {"phi", "exact", "rhs", "acoef", "bcoefx", "vfrac"}, geom, 0.0, 0);

    // for amrvis
    /*
    amrex::writeFabs(solution, "solution");
    amrex::writeFabs(bcoef, "bcoef");
    amrex::writeFabs(robin_a, "robin_a");
    amrex::writeFabs(robin_b, "robin_b");
    amrex::writeFabs(robin_f, "robin_f");
    */
  }

  amrex::Finalize();

  if (eps < 5e-3) {
    return (0);
  } else {
    return (-1);
  }
}
