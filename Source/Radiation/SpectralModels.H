#ifndef SPECTRAL_MODELS_H
#define SPECTRAL_MODELS_H

#include <AMReX.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Gpu.H>

namespace PeleRad {

namespace RadProp {
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
interpT(amrex::Real const& T, int& TindexL, amrex::Real& weight)
{
  if (T < 300) {
    TindexL = 0;
    weight = 0.0;
    return;
  }
  if (T > 2800) {
    TindexL = 125;
    weight = 1.0;
    return;
  }

  amrex::Real TindexReal = (T - 300.0) / 20.0;
  amrex::Real TindexInte = floor((T - 300.0) / 20.0);
  TindexL = static_cast<int>(TindexInte);
  weight = TindexReal - TindexInte;

  assert(weight <= 1 && weight >= 0);
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
interpk(
  int const& TindexL,
  amrex::Real const& weight,
  amrex::GpuArray<amrex::Real, 126ul> const& k)
{
  return (1.0 - weight) * k[TindexL] + weight * k[TindexL + 1];
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
getRadPropGas(
  int i,
  int j,
  int k,
  amrex::Array4<const amrex::Real> const& yco2,
  amrex::Array4<const amrex::Real> const& yh2o,
  amrex::Array4<const amrex::Real> const& yco,
  amrex::Array4<const amrex::Real> const& temp,
  amrex::Array4<const amrex::Real> const& pressure,
  amrex::Array4<amrex::Real> const& absc,
  amrex::GpuArray<amrex::Real, 126ul> const& kdataco2,
  amrex::GpuArray<amrex::Real, 126ul> const& kdatah2o,
  amrex::GpuArray<amrex::Real, 126ul> const& kdataco)
{
  int TindexL = 0;
  amrex::Real weight = 1.0;
  interpT(temp(i, j, k), TindexL, weight);

  amrex::Real kp_co2 = interpk(TindexL, weight, kdataco2);
  amrex::Real kp_h2o = interpk(TindexL, weight, kdatah2o);
  amrex::Real kp_co = interpk(TindexL, weight, kdataco);

  absc(i, j, k) =
    yco2(i, j, k) * kp_co2 + yh2o(i, j, k) * kp_h2o + yco(i, j, k) * kp_co;

  // absc(i, j, k) *= pressure(i, j, k) / 1.0e5 * 100.0; //si, in m-1
  absc(i, j, k) *=
    pressure(i, j, k) /
    1.0e5; // cgs, in cm-1; if P is in bar(PeleLM), no correction
           // is needed, otherwise the absc needs to be corrected.
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
getRadPropSoot(
  int i,
  int j,
  int k,
  amrex::Array4<const amrex::Real> const& fv,
  amrex::Array4<const amrex::Real> const& temp,
  amrex::Array4<amrex::Real> const& absc,
  amrex::GpuArray<amrex::Real, 126ul> const& kdatasoot)
{
  int TindexL = 0;
  amrex::Real weight = 1.0;

  interpT(temp(i, j, k), TindexL, weight);

  amrex::Real kp_soot = interpk(TindexL, weight, kdatasoot);

  // absc(i, j, k) += fv(i, j, k) * kp_soot * 100.0; //si, in m-1
  absc(i, j, k) += fv(i, j, k) * kp_soot; // cgs, in cm-1
}

} // namespace RadProp

} // namespace PeleRad
#endif
