#include <AMReX.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_PlotFileUtil.H>
#include <POneSingle.H>

AMREX_GPU_DEVICE AMREX_FORCE_INLINE void
actual_init_coefs(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& rhs,
  amrex::Array4<amrex::Real> const& sol,
  amrex::Array4<amrex::Real> const& alpha,
  amrex::Array4<amrex::Real> const& beta,
  amrex::Array4<amrex::Real> const& robin_a,
  amrex::Array4<amrex::Real> const& robin_b,
  amrex::Array4<amrex::Real> const& robin_f,
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& prob_lo,
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& prob_hi,
  amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx,
  amrex::Dim3 const& dlo,
  amrex::Dim3 const& dhi,
  amrex::Box const& vbx)
{
  amrex::Real const L = 2.0;
  amrex::Real const n = 3.0;
  amrex::Real const npioverL = n * M_PI / L;

  amrex::Real x = prob_lo[0] + dx[0] * (i + 0.5);
  amrex::Real y = prob_lo[1] + dx[1] * (j + 0.5);
  amrex::Real z = prob_lo[2] + dx[2] * (k + 0.5);

  beta(i, j, k) = 1.0;

  x = amrex::min(amrex::max(x, prob_lo[0]), prob_hi[0]);
  y = amrex::min(amrex::max(y, prob_lo[1]), prob_hi[1]);
  z = amrex::min(amrex::max(z, prob_lo[2]), prob_hi[2]);

  amrex::Real sincossin =
    std::sin(npioverL * x) * std::cos(npioverL * y) * std::sin(npioverL * z);

  amrex::Real coscossin =
    std::cos(npioverL * x) * std::cos(npioverL * y) * std::sin(npioverL * z);

  sol(i, j, k) = sincossin;

  if (vbx.contains(i, j, k)) {
    rhs(i, j, k) = (1.0 + npioverL * npioverL) * beta(i, j, k) * sincossin;
    alpha(i, j, k) = 1.0;
  }

  // Robin BC
  if (j >= dlo.y && j <= dhi.y && k >= dlo.z && k <= dhi.z) {
    if (i > dhi.x || i < dlo.x) {
      robin_a(i, j, k) = -1.0;
      robin_b(i, j, k) = -2.0 / 3.0;
      robin_f(i, j, k) = robin_a(i, j, k) * sol(i, j, k) +
                         robin_b(i, j, k) * npioverL * coscossin;
    }
  } else if (i >= dlo.x && i <= dhi.x && k >= dlo.z && k <= dhi.z) {
    if (j > dhi.y || j < dlo.y) {
      robin_a(i, j, k) = -1.0;
      robin_b(i, j, k) = -2.0 / 3.0;

      amrex::Real msinsinsin = -std::sin(npioverL * x) *
                               std::sin(npioverL * y) * std::sin(npioverL * z);
      robin_f(i, j, k) = robin_a(i, j, k) * sol(i, j, k) +
                         robin_b(i, j, k) * npioverL * msinsinsin;
    }
  } else if (i >= dlo.x && i <= dhi.x && j >= dlo.y && j <= dhi.y) {
    if (k > dhi.z || k < dlo.z) {
      robin_a(i, j, k) = -1.0;
      robin_b(i, j, k) = -2.0 / 3.0;

      amrex::Real sincoscos = -std::sin(npioverL * x) * std::cos(npioverL * y) *
                              std::cos(npioverL * z);

      robin_f(i, j, k) = robin_a(i, j, k) * sol(i, j, k) +
                         robin_b(i, j, k) * npioverL * sincoscos;
    }
  }
}

void
initProbABecLaplacian(
  amrex::Geometry& geom,
  amrex::MultiFab& solution,
  amrex::MultiFab& rhs,
  amrex::MultiFab& exact_solution,
  amrex::MultiFab& acoef,
  amrex::MultiFab& bcoef,
  amrex::MultiFab& robin_a,
  amrex::MultiFab& robin_b,
  amrex::MultiFab& robin_f)
{
  auto const prob_lo = geom.ProbLoArray();
  auto const prob_hi = geom.ProbHiArray();
  auto const dx = geom.CellSizeArray();
  auto const dlo = amrex::lbound(geom.Domain());
  auto const dhi = amrex::ubound(geom.Domain());
  for (amrex::MFIter mfi(rhs); mfi.isValid(); ++mfi) {
    amrex::Box const& bx = mfi.validbox();
    amrex::Box const& gbx = amrex::grow(bx, 1);
    auto const& rhsfab = rhs.array(mfi);
    auto const& solfab = solution.array(mfi);
    auto const& acfab = acoef.array(mfi);
    auto const& bcfab = bcoef.array(mfi);
    auto const& rafab = robin_a.array(mfi);
    auto const& rbfab = robin_b.array(mfi);
    auto const& rffab = robin_f.array(mfi);
    amrex::ParallelFor(gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      actual_init_coefs(
        i, j, k, rhsfab, solfab, acfab, bcfab, rafab, rbfab, rffab, prob_lo,
        prob_hi, dx, dlo, dhi, bx);
    });
  }

  amrex::MultiFab::Copy(exact_solution, solution, 0, 0, 1, 0);
  solution.setVal(0.0, 0, 1, amrex::IntVect(0));
}

void
initMeshandData(
  PeleRad::AMRParam const& amrpp,
  amrex::Geometry& geom,
  amrex::BoxArray& grids,
  amrex::DistributionMapping& dmap,
  amrex::MultiFab& solution,
  amrex::MultiFab& rhs,
  amrex::MultiFab& exact_solution,
  amrex::MultiFab& acoef,
  amrex::MultiFab& bcoef,
  amrex::MultiFab& robin_a,
  amrex::MultiFab& robin_b,
  amrex::MultiFab& robin_f)
{
  int const n_cell = amrpp.n_cell_;
  int const max_grid_size = amrpp.max_grid_size_;

  amrex::RealBox rb(
    {AMREX_D_DECL(-1.0, -1.0, -1.0)}, {AMREX_D_DECL(1.0, 1.0, 1.0)});
  amrex::Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0, 0, 0)};
  amrex::Geometry::Setup(&rb, 0, is_periodic.data());
  amrex::Box domain0(
    amrex::IntVect{AMREX_D_DECL(0, 0, 0)},
    amrex::IntVect{AMREX_D_DECL(n_cell - 1, n_cell - 1, n_cell - 1)});
  geom.define(domain0);

  grids.define(domain0);
  grids.maxSize(max_grid_size);

  amrex::IntVect ng = amrex::IntVect{1};

  dmap.define(grids);
  solution.define(grids, dmap, 1, ng);
  rhs.define(grids, dmap, 1, 0);
  exact_solution.define(grids, dmap, 1, ng);
  acoef.define(grids, dmap, 1, 0);
  bcoef.define(grids, dmap, 1, ng);
  robin_a.define(grids, dmap, 1, ng);
  robin_b.define(grids, dmap, 1, ng);
  robin_f.define(grids, dmap, 1, ng);

  initProbABecLaplacian(
    geom, solution, rhs, exact_solution, acoef, bcoef, robin_a, robin_b,
    robin_f);
}

amrex::Real
check_norm(amrex::MultiFab const& phi, amrex::MultiFab const& exact)
{
  amrex::MultiFab mf(phi.boxArray(), phi.DistributionMap(), 1, 0);
  amrex::MultiFab::Copy(mf, phi, 0, 0, 1, 0);

  amrex::MultiFab::Subtract(mf, exact, 0, 0, 1, 0);

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
  amrex::MultiFab rhs;
  amrex::MultiFab exact_solution;

  amrex::MultiFab acoef;
  amrex::MultiFab bcoef;
  amrex::MultiFab robin_a;
  amrex::MultiFab robin_b;
  amrex::MultiFab robin_f;

  // std::cout << "initialize data ... \n";
  initMeshandData(
    amrpp, geom, grids, dmap, solution, rhs, exact_solution, acoef, bcoef,
    robin_a, robin_b, robin_f);

  // std::cout << "construct the PDE ... \n";
  PeleRad::POneSingle rte(
    mlmgpp, geom, grids, dmap, solution, rhs, acoef, bcoef, robin_a, robin_b,
    robin_f);

  // std::cout << "solve the PDE ... \n";
  rte.solve();

  auto eps = check_norm(solution, exact_solution);
  eps /= static_cast<amrex::Real>(n_cell * n_cell * n_cell);
  std::cout << "n_cell=" << n_cell << ", normalized L1 norm:" << eps
            << std::endl;

  // plot results
  if (write) {
    std::cout << "write the results ... \n";
    amrex::MultiFab plotmf(grids, dmap, 4, 0);
    amrex::MultiFab::Copy(plotmf, solution, 0, 0, 1, 0);
    amrex::MultiFab::Copy(plotmf, rhs, 0, 1, 1, 0);
    amrex::MultiFab::Copy(plotmf, exact_solution, 0, 2, 1, 0);
    amrex::MultiFab::Copy(plotmf, solution, 0, 3, 1, 0);
    amrex::MultiFab::Subtract(plotmf, plotmf, 2, 3, 1, 0);

    auto const plot_file_name = amrpp.plot_file_name_;
    amrex::WriteSingleLevelPlotfile(
      plot_file_name, plotmf, {"phi", "rhs", "exact", "error"}, geom, 0.0, 0);

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
  if (eps < 1e-3) {
    return (0);
  } else {
    return (-1);
  }
}
