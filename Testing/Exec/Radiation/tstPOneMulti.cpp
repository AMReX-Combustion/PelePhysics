#include <AMReX.H>
#include <AMReX_PlotFileUtil.H>
#include <POneMulti.H>
#include <POneMultiLevbyLev.H>

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
  amrex::Vector<amrex::Geometry>& geom,
  amrex::Vector<amrex::MultiFab>& solution,
  amrex::Vector<amrex::MultiFab>& rhs,
  amrex::Vector<amrex::MultiFab>& exact_solution,
  amrex::Vector<amrex::MultiFab>& acoef,
  amrex::Vector<amrex::MultiFab>& bcoef,
  amrex::Vector<amrex::MultiFab>& robin_a,
  amrex::Vector<amrex::MultiFab>& robin_b,
  amrex::Vector<amrex::MultiFab>& robin_f)
{
  auto nlevels = geom.size();
  for (int ilev = 0; ilev < nlevels; ++ilev) {
    auto const prob_lo = geom[ilev].ProbLoArray();
    auto const prob_hi = geom[ilev].ProbHiArray();
    auto const dx = geom[ilev].CellSizeArray();
    auto const dlo = amrex::lbound(geom[ilev].Domain());
    auto const dhi = amrex::ubound(geom[ilev].Domain());
    for (amrex::MFIter mfi(rhs[ilev]); mfi.isValid(); ++mfi) {
      amrex::Box const& bx = mfi.validbox();
      amrex::Box const& gbx = amrex::grow(bx, 1);
      auto const& rhsfab = rhs[ilev].array(mfi);
      auto const& solfab = solution[ilev].array(mfi);
      auto const& acfab = acoef[ilev].array(mfi);
      auto const& bcfab = bcoef[ilev].array(mfi);
      auto const& rafab = robin_a[ilev].array(mfi);
      auto const& rbfab = robin_b[ilev].array(mfi);
      auto const& rffab = robin_f[ilev].array(mfi);
      amrex::ParallelFor(
        gbx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          actual_init_coefs(
            i, j, k, rhsfab, solfab, acfab, bcfab, rafab, rbfab, rffab, prob_lo,
            prob_hi, dx, dlo, dhi, bx);
        });
    }

    amrex::MultiFab::Copy(exact_solution[ilev], solution[ilev], 0, 0, 1, 0);
    solution[ilev].setVal(0.0, 0, 1, amrex::IntVect(0));
  }
}

void
initMeshandData(
  PeleRad::AMRParam const& amrpp,
  amrex::Vector<amrex::Geometry>& geom,
  amrex::Vector<amrex::BoxArray>& grids,
  amrex::Vector<amrex::DistributionMapping>& dmap,
  amrex::Vector<amrex::MultiFab>& solution,
  amrex::Vector<amrex::MultiFab>& rhs,
  amrex::Vector<amrex::MultiFab>& exact_solution,
  amrex::Vector<amrex::MultiFab>& acoef,
  amrex::Vector<amrex::MultiFab>& bcoef,
  amrex::Vector<amrex::MultiFab>& robin_a,
  amrex::Vector<amrex::MultiFab>& robin_b,
  amrex::Vector<amrex::MultiFab>& robin_f)
{
  int const nlevels = amrpp.max_level_ + 1;
  int const ref_ratio = amrpp.ref_ratio_;
  int const n_cell = amrpp.n_cell_;
  int const max_grid_size = amrpp.max_grid_size_;

  // initialize mesh
  // std::cout << "initialize the mesh" << std::endl;
  geom.resize(nlevels);
  grids.resize(nlevels);
  dmap.resize(nlevels);

  amrex::RealBox rb(
    {AMREX_D_DECL(-1.0, -1.0, -1.0)}, {AMREX_D_DECL(1.0, 1.0, 1.0)});
  amrex::Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(0, 0, 0)};
  amrex::Geometry::Setup(&rb, 0, is_periodic.data());
  amrex::Box domain0(
    amrex::IntVect{AMREX_D_DECL(0, 0, 0)},
    amrex::IntVect{AMREX_D_DECL(n_cell - 1, n_cell - 1, n_cell - 1)});

  amrex::Box domain = domain0;

  for (int ilev = 0; ilev < nlevels; ++ilev) {
    geom[ilev].define(domain);
    domain.refine(ref_ratio);
  }

  domain = domain0;
  for (int ilev = 0; ilev < nlevels; ++ilev) {
    grids[ilev].define(domain);
    grids[ilev].maxSize(max_grid_size);
    domain.refine(ref_ratio);
  }

  // initialize variables
  // std::cout << "initialize the data" << std::endl;

  solution.resize(nlevels);
  rhs.resize(nlevels);
  exact_solution.resize(nlevels);
  acoef.resize(nlevels);
  bcoef.resize(nlevels);
  robin_a.resize(nlevels);
  robin_b.resize(nlevels);
  robin_f.resize(nlevels);

  amrex::IntVect ng = amrex::IntVect{1};

  for (int ilev = 0; ilev < nlevels; ++ilev) {
    dmap[ilev].define(grids[ilev]);
    solution[ilev].define(grids[ilev], dmap[ilev], 1, ng);
    rhs[ilev].define(grids[ilev], dmap[ilev], 1, 0);
    exact_solution[ilev].define(grids[ilev], dmap[ilev], 1, ng);
    acoef[ilev].define(grids[ilev], dmap[ilev], 1, 0);
    bcoef[ilev].define(grids[ilev], dmap[ilev], 1, ng);
    robin_a[ilev].define(grids[ilev], dmap[ilev], 1, ng);
    robin_b[ilev].define(grids[ilev], dmap[ilev], 1, ng);
    robin_f[ilev].define(grids[ilev], dmap[ilev], 1, ng);
  }

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
  int const nlevels = amrpp.max_level_ + 1;
  int const ref_ratio = amrpp.ref_ratio_;
  int const composite_solve = mlmgpp.composite_solve_;

  amrex::Vector<amrex::Geometry> geom;
  amrex::Vector<amrex::BoxArray> grids;
  amrex::Vector<amrex::DistributionMapping> dmap;

  amrex::Vector<amrex::MultiFab> solution;
  amrex::Vector<amrex::MultiFab> rhs;
  amrex::Vector<amrex::MultiFab> exact_solution;

  amrex::Vector<amrex::MultiFab> acoef;
  amrex::Vector<amrex::MultiFab> bcoef;
  amrex::Vector<amrex::MultiFab> robin_a;
  amrex::Vector<amrex::MultiFab> robin_b;
  amrex::Vector<amrex::MultiFab> robin_f;

  //    std::cout << "initialize data ... \n";
  initMeshandData(
    amrpp, geom, grids, dmap, solution, rhs, exact_solution, acoef, bcoef,
    robin_a, robin_b, robin_f);
  //    std::cout << "construct the PDE ... \n";
  if (composite_solve) {

    std::cout << "composite solve ... \n";
    PeleRad::POneMulti rte(
      mlmgpp, geom, grids, dmap, solution, rhs, acoef, bcoef, robin_a, robin_b,
      robin_f);
    rte.solve();
  } else {
    std::cout << "level by level solve ... \n";
    PeleRad::POneMultiLevbyLev rte(
      mlmgpp, ref_ratio, geom, grids, dmap, solution, rhs, acoef, bcoef,
      robin_a, robin_b, robin_f);
    rte.solve();
  }

  amrex::Real eps = 0.0;
  amrex::Real eps_max = 0.0;
  for (int ilev = 0; ilev < nlevels; ++ilev) {
    auto dx = geom[ilev].CellSize();
    amrex::Real dvol = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);
    eps = check_norm(solution[ilev], exact_solution[ilev]);
    eps *= dvol;
    std::cout << "Level=" << ilev << ", normalized L1 norm:" << eps
              << std::endl;
    if (eps > eps_max)
      eps_max = eps;
  }

  // plot results
  if (write) {
    amrex::Vector<amrex::MultiFab> plotmf(nlevels);

    for (int ilev = 0; ilev < nlevels; ++ilev) {
      // std::cout << "write the results ... \n";
      plotmf[ilev].define(grids[ilev], dmap[ilev], 4, 0);
      amrex::MultiFab::Copy(plotmf[ilev], solution[ilev], 0, 0, 1, 0);
      amrex::MultiFab::Copy(plotmf[ilev], rhs[ilev], 0, 1, 1, 0);
      amrex::MultiFab::Copy(plotmf[ilev], exact_solution[ilev], 0, 2, 1, 0);
      amrex::MultiFab::Copy(plotmf[ilev], solution[ilev], 0, 3, 1, 0);
      amrex::MultiFab::Subtract(plotmf[ilev], plotmf[ilev], 2, 3, 1, 0);

      // For amrvis
      /*
      amrex::writeFabs(solution[ilev],
      "solution_lev"+std::to_string(ilev)); amrex::writeFabs(bcoef[ilev],
      "bcoef_lev"+std::to_string(ilev)); amrex::writeFabs(robin_a[ilev],
      "robin_a_lev"+std::to_string(ilev)); amrex::writeFabs(robin_b[ilev],
      "robin_b_lev"+std::to_string(ilev)); amrex::writeFabs(robin_f[ilev],
      "robin_f_lev"+std::to_string(ilev));
      */
    }

    auto const plot_file_name = amrpp.plot_file_name_;
    amrex::WriteMultiLevelPlotfile(
      plot_file_name, nlevels, amrex::GetVecOfConstPtrs(plotmf),
      {"phi", "rhs", "exact", "error"}, geom, 0.0,
      amrex::Vector<int>(nlevels, 0),
      amrex::Vector<amrex::IntVect>(nlevels, amrex::IntVect{ref_ratio}));
  }

  amrex::Finalize();

  if (eps_max < 1e-2) {
    return (0);
  } else {
    return (-1);
  }
}
