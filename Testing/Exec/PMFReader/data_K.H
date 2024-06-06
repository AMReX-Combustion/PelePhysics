#ifndef PROB_H
#define PROB_H

#include <AMReX_Geometry.H>
#include <AMReX_SPACE.H>
#include <PMF.H>
#include <PMFData.H>
#include <PelePhysics.H>

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
initdata(
  int i,
  int j,
  int k,
  amrex::Array4<amrex::Real> const& state,
  const amrex::Real& standoff,
  amrex::GeometryData const& geomdata,
  pele::physics::PMF::DataContainer const* pmf_data)
{
  const amrex::Real* prob_lo = geomdata.ProbLo();
  const amrex::Real* prob_hi = geomdata.ProbHi();
  const amrex::Real* dx = geomdata.CellSize();

  AMREX_D_TERM(const amrex::Real x = prob_lo[0] + (i + 0.5) * dx[0];
               , const amrex::Real y = prob_lo[1] + (j + 0.5) * dx[1];
               , const amrex::Real z = prob_lo[2] + (k + 0.5) * dx[2];);

  AMREX_D_TERM(const amrex::Real Lx = prob_hi[0] - prob_lo[0];
               , const amrex::Real Ly = prob_hi[1] - prob_lo[1];
               , const amrex::Real Lz = prob_hi[2] - prob_lo[2]);

  AMREX_D_TERM(const amrex::Real xc = prob_lo[0] + Lx / 2.0;
               , const amrex::Real yc = prob_lo[1] + Ly / 2.0;
               , const amrex::Real zc = prob_lo[2] + Lz / 2.0);

  constexpr amrex::Real Pi = 3.14159265358979323846264338327950288;
  amrex::Real molefrac[NUM_SPECIES] = {0.0};
  amrex::Real massfrac[NUM_SPECIES] = {0.0};

  // Data layout
  AMREX_D_TERM(const int Xvel = 0;, const int Yvel = 1;, const int Zvel = 2);
  const int Density = AMREX_SPACEDIM;
  const int FirstSpec = AMREX_SPACEDIM + 1;
  const int RhoH = FirstSpec + NUM_SPECIES;
  const int Temp = RhoH + 1;

  amrex::GpuArray<amrex::Real, NUM_SPECIES + 4> pmf_vals = {0.0};
  auto eos = pele::physics::PhysicsType::eos();

#if (AMREX_SPACEDIM == 2)
  amrex::Real y1 = (y - 0.5 * dx[1] - standoff) * 100;
  amrex::Real y2 = (y + 0.5 * dx[1] - standoff) * 100;
#elif (AMREX_SPACEDIM == 3)
  amrex::Real y1 = (z - 0.5 * dx[2] - standoff) * 100;
  amrex::Real y2 = (z + 0.5 * dx[2] - standoff) * 100;
#endif

  pele::physics::PMF::pmf(pmf_data, y1, y2, pmf_vals);

  state(i, j, k, Xvel) = 0;
#if (AMREX_SPACEDIM == 2)
  state(i, j, k, Yvel) = pmf_vals[1] * 1e-2;
#elif (AMREX_SPACEDIM == 3)
  state(i, j, k, Yvel) = 0;
  state(i, j, k, Zvel) = pmf_vals[1] * 1e-2;
#endif

  state(i, j, k, Temp) = pmf_vals[0];

  for (int n = 0; n < NUM_SPECIES; n++) {
    molefrac[n] = pmf_vals[3 + n];
  }
  eos.X2Y(molefrac, massfrac);

  amrex::Real P_cgs = 101325.0 * 10.0;

  amrex::Real rho_cgs;
  eos.PYT2R(P_cgs, massfrac, state(i, j, k, Temp), rho_cgs);
  state(i, j, k, Density) = rho_cgs * 1.0e3; // CGS -> MKS conversion

  eos.TY2H(state(i, j, k, Temp), massfrac, state(i, j, k, RhoH));
  state(i, j, k, RhoH) = state(i, j, k, RhoH) * 1.0e-4 *
                         state(i, j, k, Density); // CGS -> MKS conversion

  for (int n = 0; n < NUM_SPECIES; n++) {
    state(i, j, k, FirstSpec + n) = massfrac[n] * state(i, j, k, Density);
  }
}

#endif
