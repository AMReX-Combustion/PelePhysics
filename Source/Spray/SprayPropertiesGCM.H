#ifndef SPRAYPROPERTIESGCM_H
#define SPRAYPROPERTIESGCM_H

#include "SprayProperties.H"

namespace pele::physics::SprayProps {

// Struct of properties for GCM
namespace utils = pele::physics::utilities;
struct GCMLiqProps : public BaseLiqProps
{
  // Only implemented for MKS units currently
  amrex::Real ref_T = 298.15; // Reference temperature for liquid properties
  amrex::GpuArray<int, SPRAY_FUEL_NUM> family;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> critP;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> critV;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> omega;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> molar_vol;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> cp_a;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> cp_b;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> cp_c;

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real molar_vol_i(const amrex::Real T, const int spf) const
  {
    const amrex::Real Tc_i = critT_i(spf);
    amrex::Real phi_i;
    if (T > Tc_i) {
      phi_i = -std::pow(1 - (ref_T / Tc_i), 2. / 7.);
    } else {
      phi_i = std::pow(1 - (T / Tc_i), 2. / 7.) -
              std::pow(1 - (ref_T / Tc_i), 2. / 7.);
    }
    const amrex::Real z_i = 0.29056 - 0.08775 * omega[spf];
    return molar_vol[spf] * std::pow(z_i, phi_i);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real cp_i(const amrex::Real T, const int spf) const
  {
    const amrex::Real theta = (T - 298.) / 700.;
    const amrex::Real cp =
      cp_a[spf] + cp_b[spf] * theta + cp_c[spf] * theta * theta;
    return cp;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real rho_i(const amrex::Real T, const int spf) const
  {
    return mw_i(spf) / molar_vol_i(T, spf);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real lambda_i(const amrex::Real T, const int spf) const
  {
    // Uses Latini et al. method (Poling equation 10-9.1)
    const amrex::Real alpha = 1.2;
    const amrex::Real gamma = 0.167;

    // Default values for saturated hydrocarbons
    amrex::Real Astar = 0.00350;
    amrex::Real beta = 0.5;
    if (family[spf] == 1) {
      // Aromatics
      Astar = 0.0346;
      beta = 1.;
    } else if (family[spf] == 2) {
      // Cycloparaffins
      Astar = 0.0310;
      beta = 1.;
    } else if (family[spf] == 3) {
      // Olefins
      Astar = 0.0361;
      beta = 1.;
    }

    const amrex::Real mw_beta = std::pow(mw_i(spf) * 1e3, beta); // mw in g/mol
    const amrex::Real Tr = T / critT_i(spf);
    const amrex::Real A = Astar * std::pow(boilT_i(spf), alpha) /
                          (mw_beta * std::pow(critT_i(spf), gamma));
    const amrex::Real lambda =
      A * std::pow(1 - Tr, 0.38) / std::pow(Tr, 1. / 6.);

    return lambda;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real sigma_i(const amrex::Real T, const int spf) const
  {
    // Brock-Bird method
    // Predictcs surface tension in dyn/cm and is converted to SI units
    const amrex::Real Tc_i = critT_i(spf);
    const amrex::Real Pc_i = critP[spf] * 1e-5; // Pa to bar
    const amrex::Real Tb_i = boilT_i(spf);
    const amrex::Real Tr = T / Tc_i;
    const amrex::Real Tbr = Tb_i / Tc_i;
    const amrex::Real Q =
      0.1196 * (1. + (Tbr * std::log(Pc_i / 1.01325)) / (1.0 - Tbr)) - 0.279;

    const amrex::Real sigma = std::pow(Pc_i, 2. / 3.) *
                              std::pow(Tc_i, 1. / 3.) * Q *
                              std::pow(1 - Tr, 11. / 9.);
    return sigma * 1e-3; // Convert from dyn/cm to N/m;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real nu_i(const amrex::Real T, const int spf) const
  {
    // Dutt's equation (4.23) from "Viscosity of Liquids"
    // Predicts viscosity in mm^2/s and is converted to SI units
    const amrex::Real T_cels = utils::K2C(T);
    const amrex::Real Tb_cels = utils::K2C(boilT_i(spf));
    const amrex::Real rhs =
      -3.0171 + (442.78 + 1.6452 * Tb_cels) / (T_cels + 239. - 0.19 * Tb_cels);
    const amrex::Real nu = std::exp(rhs) * 1.0e-6; // Convert to m^2/s
    return nu;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real mu_i(const amrex::Real T, const int spf) const
  {
    const amrex::Real nu = nu_i(T, spf);
    const amrex::Real rho = rho_i(T, spf);
    return nu * rho;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real
  psat_i(const amrex::Real T, int spf, const amrex::Real /*part_latent*/) const
  {
    // Lee-Kesler correlation
    const amrex::Real Tr = T / critT_i(spf);
    const amrex::Real f0 = 5.92714 - (6.09648 / Tr) - 1.28862 * std::log(Tr) +
                           0.169347 * std::pow(Tr, 6);
    const amrex::Real f1 = 15.2518 - (15.6875 / Tr) - 13.4721 * std::log(Tr) +
                           0.43577 * std::pow(Tr, 6);
    const amrex::Real rhs = std::exp(f0 + omega[spf] * f1);
    return critP[spf] * rhs;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real latent_i(
    const amrex::Real T, const amrex::Real h_gas_spf, const int spf) const
  {
    amrex::Real latent_heat_i =
      h_gas_spf + latentRef_minus_gasRefH_i[spf] - cp_i(T, spf) * (T - ref_T);
    return latent_heat_i;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real latent_correlation_i(const amrex::Real T, const int spf) const
  {
    // Correlation to be used with manifold model
    const amrex::Real Tc_i = critT_i(spf);
    const amrex::Real Tb_i = boilT_i(spf);
    const amrex::Real Tr = T / Tc_i;
    const amrex::Real Trb = Tb_i / Tc_i;
    amrex::Real latent_heat_i = 0.;
    if (T > Tc_i) {
      latent_heat_i = 0.0;
    } else {
      latent_heat_i = latent[spf] * std::pow((1.0 - Tr) / (1.0 - Trb), 0.38);
    }
    return latent_heat_i;
  }

  // Mixture functions use RealArrayLike to allow for different types of arrays
  // with an [] operator that returns amrex::Real

  template <typename RealArrayLike>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
  mw_mix(const RealArrayLike& Y) const
  {
    amrex::Real mw_part = 0.;
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      mw_part += Y[spf] / mw_i(spf);
    }
    mw_part = 1. / mw_part;
    return mw_part;
  }

  template <typename RealArrayLike>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
  critT_mix(const RealArrayLike& Y) const
  {
    amrex::Real Tcrit = 0.;
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      Tcrit += Y[spf] * critT_i(spf);
    }
    return Tcrit;
  }

  template <typename RealArrayLike>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
  boilT_mix(const RealArrayLike& Y, const amrex::Real* cBoilT) const
  {
    amrex::Real Tboil = 0.;
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      Tboil += Y[spf] * cBoilT[spf];
    }
    return Tboil;
  }

  template <typename RealArrayLike>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
  cp_mix(const RealArrayLike& Y, const amrex::Real T, const amrex::Real* cBoilT)
    const
  {
    amrex::Real Cp = 0.;
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      Cp += Y[spf] * cp_i(amrex::min(T, cBoilT[spf]), spf);
    }
    return Cp;
  }

  // rho_mix can be called as rho_mix(Y,T) or rho_mix(Y,T,cBoilT)
  template <typename RealArrayLike>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
  rho_mix(const RealArrayLike& Y, const amrex::Real T) const
  {
    amrex::Real rho = 0.;
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      rho += Y[spf] / rho_i(T, spf);
    }
    rho = 1. / rho;
    return rho;
  }

  template <typename RealArrayLike>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real rho_mix(
    const RealArrayLike& Y,
    const amrex::Real T,
    const amrex::Real* cBoilT) const
  {
    amrex::Real rho = 0.;
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      rho += Y[spf] / rho_i(amrex::min(T, cBoilT[spf]), spf);
    }
    rho = 1. / rho;
    return rho;
  }

  template <typename RealArrayLike>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real lambda_mix(
    const RealArrayLike& Y,
    const amrex::Real T,
    const amrex::Real* cBoilT) const
  {

    amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> X = {{0.}};
    Y2X(Y, X);

    amrex::Real sum = 0.;
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      sum += X[spf] * std::pow(lambda_i(amrex::min(T, cBoilT[spf]), spf), -2.);
    }
    return std::pow(sum, -0.5);
  }

  template <typename RealArrayLike>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
  mu_mix(const RealArrayLike& Y, const amrex::Real T, const amrex::Real* cBoilT)
    const
  {
    // Kendall-Monroe mixing rule for kinematic viscosity
    amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> X = {{0.}};
    Y2X(Y, X);
    amrex::Real nu = 0.;
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      nu += X[spf] * std::pow(nu_i(amrex::min(T, cBoilT[spf]), spf), 1. / 3.);
    }
    nu = std::pow(nu, 3.);
    return nu * rho_mix(Y, T, cBoilT);
  }

  template <typename RealArrayLike>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real sigma_mix(
    const RealArrayLike& Y,
    const amrex::Real T,
    const amrex::Real* cBoilT) const
  {
    // Method recommended by Hugill and van Welsenes (1986).
    amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> X = {{0.}};
    Y2X(Y, X);

    amrex::Real sigma = 0.;
    for (int i = 0; i < SPRAY_FUEL_NUM; ++i) {
      for (int j = 0; j < SPRAY_FUEL_NUM; ++j) {
        sigma += 0.5 * X[i] * X[j] *
                 (sigma_i(amrex::min(T, cBoilT[i]), i) +
                  sigma_i(amrex::min(T, cBoilT[j]), j));
      }
    }

    return sigma;
  }
};

// Define the type of liquid properties based on the compilation flag
#ifdef PELE_SPRAY_GCM
#ifndef PELELM_USE_SPRAY
static_assert(
  false, "GCM liquid property model can currently only be used with PeleLMeX");
#endif
using LiqPropType = GCMLiqProps;
#endif

// Specialization for GCMLiqProps
template <>
struct InitLiqProps<GCMLiqProps>
{
  void
  operator()(GCMLiqProps* ldata, const std::vector<std::string>& fuel_names)
  {
    amrex::ParmParse pp("particles");
    getInpVal(ldata->family.data(), pp, fuel_names.data(), "family");
    getInpVal(ldata->critT.data(), pp, fuel_names.data(), "crit_temp");
    getInpVal(ldata->critP.data(), pp, fuel_names.data(), "crit_press");
    getInpVal(ldata->critV.data(), pp, fuel_names.data(), "crit_vol");
    getInpVal(ldata->boilT.data(), pp, fuel_names.data(), "boil_temp");
    getInpVal(ldata->omega.data(), pp, fuel_names.data(), "acentric_factor");
    getInpVal(ldata->molar_vol.data(), pp, fuel_names.data(), "molar_vol");
    getInpVal(ldata->cp_a.data(), pp, fuel_names.data(), "cp_a");
    getInpVal(ldata->cp_b.data(), pp, fuel_names.data(), "cp_b");
    getInpVal(ldata->cp_c.data(), pp, fuel_names.data(), "cp_c");
    getInpVal(ldata->latent.data(), pp, fuel_names.data(), "latent");
  }

  static void init_breakupsplash(
    GCMLiqProps* /*ldata*/, std::vector<std::string>& /*fuel_names*/)
  {
    // GCM calculates surface tension and viscosity from other properties
  }
};

} // namespace pele::physics::SprayProps
#endif