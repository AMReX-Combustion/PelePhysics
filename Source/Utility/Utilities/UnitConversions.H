#ifndef UNITCONVERSIONS_H
#define UNITCONVERSIONS_H

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>

namespace pele::physics::utilities {
namespace cgs2mks {
// CGS to MKS conversions
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Length(const amrex::Real L_cgs)
{
  return L_cgs * 1.0e-2; // cm to m
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Length(const amrex::Real L_cgs, const amrex::Real n)
{
  return L_cgs * std::pow(1.0e-2, n); // cm^n to m^n
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Mass(const amrex::Real M_cgs)
{
  return M_cgs * 1.0e-3; // g to kg
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Mass(const amrex::Real M_cgs, const amrex::Real n)
{
  return M_cgs * std::pow(1.0e-3, n); // g^n to kg^n
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Energy(const amrex::Real E_cgs)
{
  return E_cgs * 1.0e-7; // erg to J
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Rho(const amrex::Real rho_cgs)
{
  return rho_cgs * 1.0e3; // g/cm^3 to kg/m^3
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
U(const amrex::Real u_cgs)
{
  return u_cgs * 1.0e-2; // cm/s to m/s
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
RhoU(const amrex::Real rhou_cgs)
{
  return rhou_cgs * 1.0e1; // g/cm^2/s to kg/m^2/s
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
P(const amrex::Real p_cgs)
{
  return p_cgs * 1.0e-1; // dyne/cm^2 to Pa
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
H(const amrex::Real h_cgs)
{
  return h_cgs * 1.0e-4; // erg/g to J/kg
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
RhoH(const amrex::Real rhoh_cgs)
{
  return rhoh_cgs * 1.0e1; // erg/cm^3 to J/m^3
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Nu(const amrex::Real nu_cgs)
{
  return nu_cgs * 1.0e-4; // cm^2/s to m^2/s
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Mu(const amrex::Real mu_cgs)
{
  return mu_cgs * 1.0e-1; // g/cm/s to kg/m/s
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Cp(const amrex::Real cp_cgs)
{
  return cp_cgs * 1.0e-4; // erg/g/K to J/kg/K
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Alpha(const amrex::Real a_cgs)
{
  return a_cgs * 1.0e-4; // cm^2/s to m^2/s
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Lambda(const amrex::Real l_cgs)
{
  return l_cgs * 1.0e-5; // g/cm/s^3/K to W/m/K
}
} // namespace cgs2mks

namespace mks2cgs {
// MKS to CGS conversions
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Length(const amrex::Real L_mks)
{
  return L_mks * 1.0e2; // m to cm
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Length(const amrex::Real L_mks, const amrex::Real n)
{
  return L_mks * std::pow(1.0e2, n); // m^n to cm^n
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Mass(const amrex::Real M_mks)
{
  return M_mks * 1.0e3; // kg to g
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Mass(const amrex::Real M_mks, const amrex::Real n)
{
  return M_mks * std::pow(1.0e3, n); // kg^n to g^n
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Energy(const amrex::Real E_mks)
{
  return E_mks * 1.0e7; // J to erg
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Rho(const amrex::Real Rho_mks)
{
  return Rho_mks * 1.0e-3; // kg/m^3 to g/cm^3
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
U(const amrex::Real u_mks)
{
  return u_mks * 1.0e2; // m/s to cm/s
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
RhoU(const amrex::Real Rhou_mks)
{
  return Rhou_mks * 1.0e-1; // kg/m^2/s to g/cm^2/s
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
P(const amrex::Real p_mks)
{
  return p_mks * 1.0e1; // Pa to dyne/cm^2
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
H(const amrex::Real h_mks)
{
  return h_mks * 1.0e4; // J/kg to erg/g
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
RhoH(const amrex::Real Rhoh_mks)
{
  return Rhoh_mks * 1.0e1; // J/m^3 to erg/cm^3
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Nu(const amrex::Real Nu_mks)
{
  return Nu_mks * 1.0e4; // m^2/s to cm^2/s
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Mu(const amrex::Real Mu_mks)
{
  return Mu_mks * 1.0e1; // kg/m/s to g/cm/s
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Cp(const amrex::Real cp_mks)
{
  return cp_mks * 1.0e4; // J/kg/K to erg/g/K
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Alpha(const amrex::Real a_mks)
{
  return a_mks * 1.0e4; // m^2/s to cm^2/s
}

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
Lambda(const amrex::Real l_mks)
{
  return l_mks * 1.0e5; // W/m/K to g/cm/s^3/K
}
} // namespace mks2cgs
} // namespace pele::physics::utilities
#endif