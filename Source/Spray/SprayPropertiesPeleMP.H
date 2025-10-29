#ifndef SPRAYPROPERTIESPELEMP_H
#define SPRAYPROPERTIESPELEMP_H
#include "SprayProperties.H"

namespace pele::physics::SprayProps {

// Original spray properties ported from PeleMP to PelePhysics
struct MPLiqProps : public BaseLiqProps
{
  amrex::Real ref_T = 300.; // Reference temperature for liquid properties
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> cp;
  amrex::Real sigma = -1.; // Surface tension
  // 3 coefficients for Antoine equation and conversion to appropriate units
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM * 4> psat_coef;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM * 4> rho_coef;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM * 4> lambda_coef;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM * 4> mu_coef;

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
    const amrex::Real a = rho_coef[4 * spf];
    const amrex::Real b = rho_coef[4 * spf + 1];
    const amrex::Real c = rho_coef[4 * spf + 2];
    const amrex::Real d = rho_coef[4 * spf + 3];
    return a + T * (b + T * (c + T * d));
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real lambda_i(const amrex::Real T, const int spf) const
  {
    const amrex::Real a = lambda_coef[4 * spf];
    const amrex::Real b = lambda_coef[4 * spf + 1];
    const amrex::Real c = lambda_coef[4 * spf + 2];
    const amrex::Real d = lambda_coef[4 * spf + 3];
    return a + T * (b + T * (c + T * d));
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real mu_i(const amrex::Real T, const int spf) const
  {
    const amrex::Real a = mu_coef[4 * spf];
    const amrex::Real b = mu_coef[4 * spf + 1];
    const amrex::Real c = mu_coef[4 * spf + 2];
    const amrex::Real d = mu_coef[4 * spf + 3];
    return a + ((d / T + c) / T + b) / T;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real
  psat_i(const amrex::Real T, int spf, const amrex::Real part_latent) const
  {
    amrex::Real pres_sat_i = 0.;
    if (psat_coef[4 * spf + 3] == 0.) {
      // Using the Clasius-Clapeyron relation
      const amrex::Real RU = pele::physics::Constants::RU * SprayUnits::ru_conv;
      const amrex::Real PATM =
        pele::physics::Constants::PATM * SprayUnits::pres_conv;
      const amrex::Real boilT_ref = boilT_i(spf);
      pres_sat_i =
        PATM *
        std::exp(part_latent * mw_i(spf) / RU * (1. / boilT_ref - 1. / T));
    } else {
      // Using the Antoine equation
      const amrex::Real a = psat_coef[4 * spf];
      const amrex::Real b = psat_coef[4 * spf + 1];
      const amrex::Real c = psat_coef[4 * spf + 2];
      const amrex::Real d = psat_coef[4 * spf + 3];
      pres_sat_i = d * std::pow(10., a - b / (T + c));
    }
    return pres_sat_i;
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

// Define the type of liquid properties based on the compilation flag
#ifndef PELE_SPRAY_GCM
using LiqPropType = MPLiqProps;
#endif

// Specialization for MPLiqProps
template <>
struct InitLiqProps<MPLiqProps>
{
  void operator()(MPLiqProps* ldata, const std::vector<std::string>& fuel_names)
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
      const std::string var_read = fuel_names[spf] + "_" + varname;
      const int numvals = ppp.countval(var_read);
      if (numvals == 4) {
        std::vector<amrex::Real> inp_coef(4, 0.);
        if (is_required) {
          ppp.getarr(var_read, inp_coef);
        } else {
          ppp.queryarr(var_read, inp_coef);
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
          ppp.get(var_read, inp_coef);
        } else {
          ppp.query(var_read, inp_coef);
        }
        coef[4 * spf] = inp_coef;
      }
    }
  }

  static void
  init_breakupsplash(MPLiqProps* ldata, std::vector<std::string>& fuel_names)
  {
    // Check for sigma and mu coeffs (copied from PeleMP implementation)
    amrex::ParmParse pp("particles");
    bool wrong_data = false;
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      const std::string var_read = fuel_names[spf] + "_mu";
      if (!pp.contains(var_read)) {
        wrong_data = true;
      }
    }
    if (wrong_data || !pp.contains("fuel_sigma")) {
      amrex::Abort(
        "fuel_sigma and mu coeffs must be set for splash or breakup model.");
    }
    pp.get("fuel_sigma", ldata->sigma);
  }
};

} // namespace pele::physics::SprayProps
#endif
