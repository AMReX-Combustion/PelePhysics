#ifndef SPRAYPROPERTIES_H
#define SPRAYPROPERTIES_H
#include "UnitConversions.H"

// Structs for liquid fuel properties
namespace pele::physics::SprayProps {

// Original spray properties ported from PeleMP to PelePhysics
struct MPLiqProps
{
  amrex::Real ref_T = 300.; // Reference temperature for liquid properties
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> critT;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> boilT;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> cp;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> latent;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> latentRef_minus_gasRefH_i;
  amrex::Real sigma = -1.; // Surface tension
  // 3 coefficients for Antoine equation and conversion to appropriate units
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM * 4> psat_coef;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM * 4> rho_coef;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM * 4> lambda_coef;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM * 4> mu_coef;

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real critT_i(const int spf) const { return critT[spf]; }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real boilT_i(const int spf) const { return boilT[spf]; }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real cp_i(const amrex::Real /*T*/, const int spf) const
  {
    return cp[spf];
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real rho_i(const amrex::Real T, const int spf) const
  {
    amrex::Real a = rho_coef[4 * spf];
    amrex::Real b = rho_coef[4 * spf + 1];
    amrex::Real c = rho_coef[4 * spf + 2];
    amrex::Real d = rho_coef[4 * spf + 3];
    return a + T * (b + T * (c + T * d));
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real lambda_i(const amrex::Real T, const int spf) const
  {
    amrex::Real a = lambda_coef[4 * spf];
    amrex::Real b = lambda_coef[4 * spf + 1];
    amrex::Real c = lambda_coef[4 * spf + 2];
    amrex::Real d = lambda_coef[4 * spf + 3];
    return a + T * (b + T * (c + T * d));
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real mu_i(const amrex::Real T, const int spf) const
  {
    amrex::Real a = mu_coef[4 * spf];
    amrex::Real b = mu_coef[4 * spf + 1];
    amrex::Real c = mu_coef[4 * spf + 2];
    amrex::Real d = mu_coef[4 * spf + 3];
    return a + ((d / T + c) / T + b) / T;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real psat_i(const amrex::Real T, const int spf) const
  {
    amrex::Real a = psat_coef[4 * spf];
    amrex::Real b = psat_coef[4 * spf + 1];
    amrex::Real c = psat_coef[4 * spf + 2];
    amrex::Real d = psat_coef[4 * spf + 3];
    return d * std::pow(10., a - b / (T + c));
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real latent_i(const amrex::Real T, const int spf) const
  {
    // Since we only know the latent heat at the reference temperature,
    // modify Watsons power law to find latent heat at any temperature
    return latent[spf] *
           std::pow((critT_i(spf) - ref_T) / (critT_i(spf) - T), -0.38);
  }

  // Mixture functions use RealArrayLike to allow for different types of arrays
  // with an [] operator that returns amrex::Real

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
  template <
    typename RealArrayLike,
    typename RealArrayLike2 = amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM>>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real rho_mix(
    const RealArrayLike& Y,
    const amrex::Real T,
    const RealArrayLike2& cBoilT = [] {
      amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> arr;
      arr.fill(AMREX_REAL_MAX);
      return arr;
    }()) const
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
    amrex::Real lambda = 0.;
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      lambda += Y[spf] * lambda_i(amrex::min(T, cBoilT[spf]), spf);
    }
    return lambda;
  }

  template <typename RealArrayLike>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
  mu_mix(const RealArrayLike& Y, const amrex::Real T, const amrex::Real* cBoilT)
    const
  {
    amrex::Real mu = 0.;
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      mu += Y[spf] * mu_i(amrex::min(T, cBoilT[spf]), spf);
    }
    return mu;
  }

  template <typename RealArrayLike>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real sigma_mix(
    const RealArrayLike& /*Y*/,
    const amrex::Real /*T*/,
    const amrex::Real* /*cBoilT*/) const
  {
    return sigma;
  }
};

// Struct of properties for GCM
namespace utils = pele::physics::utilities;
struct GCMLiqProps
{
  // Only implemented for MKS units currently
  amrex::Real ref_T = 298.15; // Reference temperature for liquid properties
  amrex::GpuArray<int, SPRAY_FUEL_NUM> family;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> molar_wt;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> critT;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> critP;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> critV;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> boilT;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> omega;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> molar_vol;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> cp_a;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> cp_b;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> cp_c;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> latent;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> latentRef_minus_gasRefH_i;

  template <typename RealArrayLike>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
  Y2X(const RealArrayLike& Y, RealArrayLike& X)
  {
    // EOS::Y2X uses mw of gas species, which may not include all liquid species
    amrex::Real Mbar = 0.;
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      Mbar += Y[spf] / molar_wt_i(spf);
    }
    Mbar = 1. / Mbar;
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      X[spf] = Mbar * Y[spf] / molar_wt_i(spf);
    }
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real molar_wt_i(const int spf) const { return molar_wt[spf]; }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real critT_i(const int spf) const { return critT[spf]; }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real boilT_i(const int spf) const { return boilT[spf]; }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real molar_vol_i(const amrex::Real T, const int spf) const
  {
    amrex::Real Tc_i = critT_i(spf);
    amrex::Real phi_i;
    if (T > Tc_i) {
      phi_i = -std::pow(1 - (ref_T / Tc_i), 2. / 7.);
    } else {
      phi_i = std::pow(1 - (T / Tc_i), 2. / 7.) -
              std::pow(1 - (ref_T / Tc_i), 2. / 7.);
    }
    amrex::Real z_i = 0.29056 - 0.08775 * omega[spf];
    return molar_vol[spf] * std::pow(z_i, phi_i);
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real cp_i(const amrex::Real T, const int spf) const
  {
    amrex::Real theta = (T - 298.) / 700.;
    amrex::Real cp = cp_a[spf] + cp_b[spf] * theta + cp_c[spf] * theta * theta;
    return cp;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real rho_i(const amrex::Real T, const int spf) const
  {
    return molar_wt_i(spf) / molar_vol_i(T, spf);
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

    amrex::Real mw_beta = std::pow(molar_wt_i(spf) * 1e3, beta); // mw in g/mol
    amrex::Real Tr = T / critT_i(spf);
    amrex::Real A = Astar * std::pow(boilT_i(spf), alpha) /
                    (mw_beta * std::pow(critT_i(spf), gamma));
    amrex::Real lambda = A * std::pow(1 - Tr, 0.38) / std::pow(Tr, 1. / 6.);

    return lambda;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real sigma_i(const amrex::Real T, const int spf) const
  {
    // Brock-Bird method
    // Predictcs surface tension in dyn/cm and is converted to SI units
    amrex::Real Tc_i = critT_i(spf);
    amrex::Real Pc_i = critP[spf] * 1e-5; // Pa to bar
    amrex::Real Tb_i = boilT_i(spf);
    amrex::Real omega_i = omega[spf];
    amrex::Real Tr = T / Tc_i;
    amrex::Real Tbr = Tb_i / Tc_i;
    amrex::Real Q =
      0.1196 * (1. + (Tbr * std::log(Pc_i / 1.01325)) / (1.0 - Tbr)) - 0.279;

    amrex::Real sigma = std::pow(Pc_i, 2. / 3.) * std::pow(Tc_i, 1. / 3.) * Q *
                        std::pow(1 - Tr, 11. / 9.);
    return sigma *= 1e-3; // Convert from dyn/cm to N/m;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real nu_i(const amrex::Real T, const int spf) const
  {
    // Dutt's equation (4.23) from "Viscosity of Liquids"
    // Predicts viscosity in mm^2/s and is converted to SI units
    amrex::Real T_cels = utils::K2C(T);
    amrex::Real Tb_cels = utils::K2C(
      boilT_i(spf)); // Flag: Should we use cBoilT here to account for high P?
    amrex::Real rhs =
      -3.0171 + (442.78 + 1.6452 * Tb_cels) / (T_cels + 239. - 0.19 * Tb_cels);
    amrex::Real nu = std::exp(rhs) * 1.0e-6; // Convert from mm^2/s to m^2/s
    return nu;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real mu_i(const amrex::Real T, const int spf) const
  {
    amrex::Real nu =
      nu_i(T, spf); // Flag: Should we use cBoilT here to account for high P?
    amrex::Real rho = rho_i(T, spf);
    return nu * rho;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real psat_i(const amrex::Real T, const int spf) const
  {
    // Lee-Kesler correlation
    amrex::Real Tr = T / critT_i(spf);
    amrex::Real f0 = 5.92714 - (6.09648 / Tr) - 1.28862 * std::log(Tr) +
                     0.169347 * std::pow(Tr, 6);
    amrex::Real f1 = 15.2518 - (15.6875 / Tr) - 13.4721 * std::log(Tr) +
                     0.43577 * std::pow(Tr, 6);
    amrex::Real rhs = std::exp(f0 + omega[spf] * f1);
    return critP[spf] * rhs;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real latent_i(const amrex::Real T, const int spf) const
  {
    amrex::Real Tc_i = critT_i(spf);
    amrex::Real Tb_i = boilT_i(spf);
    amrex::Real Tr = T / Tc_i;
    amrex::Real Trb = Tb_i / Tc_i;

    if (T > Tc_i) {
      return 0.0;
    } else {
      return latent[spf] * std::pow((1.0 - Tr) / (1.0 - Trb), 0.38);
    }
  }

  // Mixture functions use RealArrayLike to allow for different types of arrays
  // with an [] operator that returns amrex::Real

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
  template <
    typename RealArrayLike,
    typename RealArrayLike2 = amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM>>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real rho_mix(
    const RealArrayLike& Y,
    const amrex::Real T,
    const RealArrayLike2& cBoilT = [] {
      amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> arr;
      arr.fill(AMREX_REAL_MAX);
      return arr;
    }()) const
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
#ifndef SPRAY_GCM
using LiqPropType = MPLiqProps;
#else
using LiqPropType = GCMLiqProps;
#endif

template <typename LiqPropType>
struct InitLiqProps
{
  void
  operator()(LiqPropType* /*ldata*/, std::vector<std::string>& /*fuel_names*/)
  {
    amrex::Abort("InitLiqProps not implemented for this LiqPropType.");
  }
};

// Specialization for MPLiqProps
template <>
struct InitLiqProps<MPLiqProps>
{
  void operator()(MPLiqProps* ldata, std::vector<std::string>& fuel_names)
  {
    amrex::ParmParse pp("particles");
    getInpVal(ldata->critT.data(), pp, fuel_names.data(), "crit_temp");
    getInpVal(ldata->boilT.data(), pp, fuel_names.data(), "boil_temp");
    getInpVal(ldata->cp.data(), pp, fuel_names.data(), "cp");
    getInpVal(ldata->latent.data(), pp, fuel_names.data(), "latent");
    pp.get("fuel_ref_temp", ldata->ref_T);

    getInpCoef(ldata->lambda_coef.data(), pp, fuel_names.data(), "lambda");
    getInpCoef(ldata->psat_coef.data(), pp, fuel_names.data(), "psat");
    getInpCoef(ldata->rho_coef.data(), pp, fuel_names.data(), "rho", true);
    getInpCoef(ldata->mu_coef.data(), pp, fuel_names.data(), "mu");
  }

  static void getInpCoef(
    amrex::Real* coef,
    const amrex::ParmParse& ppp,
    const std::string* fuel_names,
    const std::string& varname,
    bool is_required = false)
  {
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      std::string var_read = fuel_names[spf] + "_" + varname;
      int numvals = ppp.countval(var_read.c_str());
      if (numvals == 4) {
        std::vector<amrex::Real> inp_coef(4, 0.);
        if (is_required) {
          ppp.getarr(var_read.c_str(), inp_coef);
        } else {
          ppp.queryarr(var_read.c_str(), inp_coef);
        }
        for (int i = 0; i < 4; ++i) {
          coef[4 * spf + i] = inp_coef[i];
        }
      } else if (numvals == 1) {
        amrex::Real inp_coef = 0.;
        for (int i = 0; i < 4; ++i) {
          coef[4 * spf + i] = 0.;
        }
        if (is_required) {
          ppp.get(var_read.c_str(), inp_coef);
        } else {
          ppp.query(var_read.c_str(), inp_coef);
        }
        coef[4 * spf] = inp_coef;
      }
    }
  }

  static void getInpVal(
    amrex::Real* coef,
    const amrex::ParmParse& ppp,
    const std::string* fuel_names,
    const std::string& varname)
  {
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      std::string var_read = fuel_names[spf] + "_" + varname;
      ppp.get(var_read.c_str(), coef[spf]);
    }
  }
};

// Specialization for GCMLiqProps
template <>
struct InitLiqProps<GCMLiqProps>
{
  void
  operator()(GCMLiqProps* /*ldata*/, std::vector<std::string>& /*fuel_names*/)
  {
  }
};

} // namespace pele::physics::SprayProps

#endif