#ifndef SPRAYPROPERTIES_H
#define SPRAYPROPERTIES_H

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
  amrex::Real cp_i(const amrex::Real& /*T*/, const int spf) const
  {
    return cp[spf];
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real rho_i(const amrex::Real& T, const int spf) const
  {
    amrex::Real a = rho_coef[4 * spf];
    amrex::Real b = rho_coef[4 * spf + 1];
    amrex::Real c = rho_coef[4 * spf + 2];
    amrex::Real d = rho_coef[4 * spf + 3];
    return a + T * (b + T * (c + T * d));
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real lambda_i(const amrex::Real& T, const int spf) const
  {
    amrex::Real a = lambda_coef[4 * spf];
    amrex::Real b = lambda_coef[4 * spf + 1];
    amrex::Real c = lambda_coef[4 * spf + 2];
    amrex::Real d = lambda_coef[4 * spf + 3];
    return a + T * (b + T * (c + T * d));
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real mu_i(const amrex::Real& T, const int spf) const
  {
    amrex::Real a = mu_coef[4 * spf];
    amrex::Real b = mu_coef[4 * spf + 1];
    amrex::Real c = mu_coef[4 * spf + 2];
    amrex::Real d = mu_coef[4 * spf + 3];
    return a + ((d / T + c) / T + b) / T;
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real psat_i(const amrex::Real& T, const int spf) const
  {
    amrex::Real a = psat_coef[4 * spf];
    amrex::Real b = psat_coef[4 * spf + 1];
    amrex::Real c = psat_coef[4 * spf + 2];
    amrex::Real d = psat_coef[4 * spf + 3];
    return d * std::pow(10., a - b / (T + c));
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real latent_i(const amrex::Real& T, const int spf) const
  {
    // Since we only know the latent heat at the reference temperature,
    // modify Watsons power law to find latent heat at any temperature
    return latent[spf] *
           std::pow((critT_i(spf) - ref_T) / (critT_i(spf) - T), -0.38);
  }

  // RealArrayLike can be anything with an [] operator that returns
  // amrex::Real&
  template <typename RealArrayLike>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE amrex::Real
  rho_mix(const RealArrayLike& Y, const amrex::Real& T) const
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
    const amrex::Real& T,
    const amrex::Real* cBoilT) const
  {
    amrex::Real rho = 0.;
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      rho += Y[spf] / rho_i(amrex::min(T, cBoilT[spf]), spf);
    }
    rho = 1. / rho;
    return rho;
  }
};

// Struct of properties for GCM
struct GCMLiqProps
{
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
