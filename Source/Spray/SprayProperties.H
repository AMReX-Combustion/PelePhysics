#ifndef SPRAYPROPERTIES_H
#define SPRAYPROPERTIES_H
#include "UnitConversions.H"

namespace c2m = pele::physics::utilities::cgs2mks;

// Units for sprays
struct SprayUnits
{
#ifdef PELELM_USE_SPRAY
  // For converting CGS to MKS
  static constexpr amrex::Real ru_conv = c2m::Energy(1.0);
  static constexpr amrex::Real eng_conv = c2m::H(1.0);
  static constexpr amrex::Real rho_conv = c2m::Rho(1.0);
  static constexpr amrex::Real mass_conv = c2m::Mass(1.0);
  static constexpr amrex::Real rhod_conv =
    c2m::Rho(1.0) * c2m::Alpha(1.0); // rho D conversion
  static constexpr amrex::Real mu_conv = c2m::Mu(1.0);
  static constexpr amrex::Real lambda_conv = c2m::Lambda(1.0);
  static constexpr amrex::Real pres_conv = c2m::P(1.0);
  static constexpr amrex::Real len_conv = c2m::Length(1.0);
#else
  static constexpr amrex::Real ru_conv = 1.;
  static constexpr amrex::Real eng_conv = 1.;
  static constexpr amrex::Real rho_conv = 1.;
  static constexpr amrex::Real mass_conv = 1.;
  static constexpr amrex::Real rhod_conv = 1.;
  static constexpr amrex::Real mu_conv = 1.;
  static constexpr amrex::Real lambda_conv = 1.;
  static constexpr amrex::Real pres_conv = 1.;
  static constexpr amrex::Real len_conv = 1.;
#endif
  static constexpr amrex::Real min_mass = 1.E-13 * mass_conv;
};

// Structs for liquid fuel properties
namespace pele::physics::SprayProps {

template <typename parsedqty>
// Function for reading particle input values
static void
getInpVal(
  parsedqty* val,
  const amrex::ParmParse& ppp,
  const std::string* fuel_names,
  const std::string& varname)
{
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    const std::string var_read = fuel_names[spf] + "_" + varname;
    ppp.get(var_read, val[spf]);
  }
}

// Base struct for liquid properties
struct BaseLiqProps
{
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> mw;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> critT;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> boilT;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> latent;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> latentRef_minus_gasRefH_i;

  template <typename RealArrayLike>
  AMREX_FORCE_INLINE void init_mw(const RealArrayLike& mw_gas, const int* indx)
  {
#if SPRAY_FUEL_NUM <= NUM_SPECIES
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      const int fspec = indx[spf];
      mw[spf] = mw_gas[fspec] * SprayUnits::mass_conv;
    }
#ifdef PELE_SPRAY_GCM
    {
      amrex::ParmParse pp("particles");
      std::vector<std::string> fuel_names;
      pp.getarr("fuel_species", fuel_names);
      amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> mw_fl;
      getInpVal(mw_fl.data(), pp, fuel_names.data(), "molar_weight");
      for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
        const amrex::Real diff = std::abs(mw_fl[spf] - mw[spf]);
        if (diff > 1.E-3) {
          amrex::Warning(
            "Warning: GCM species " + fuel_names[spf] + " mw " +
            std::to_string(mw_fl[spf]) +
            " does not match gas-phase species mw " + std::to_string(mw[spf]));
        }
      }
    }
#endif // PELE_SPRAY_GCM
#else
    amrex::Abort("SPRAY_FUEL_NUM > NUM_SPECIES");
    // Will require mw to be read from input file in many-to-one case
#endif
  }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real mw_i(const int spf) const { return mw[spf]; }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real critT_i(const int spf) const { return critT[spf]; }

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  amrex::Real boilT_i(const int spf) const { return boilT[spf]; }

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

  template <typename RealArrayLike, typename RealArrayLike2>
  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void
  Y2X(const RealArrayLike& Y, RealArrayLike2& X) const
  {
    // EOS::Y2X uses mw of gas species, which may not include all liquid species
    amrex::Real Mbar = 0.;
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      Mbar += Y[spf] / mw_i(spf);
    }
    Mbar = 1. / Mbar;
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      X[spf] = Mbar * Y[spf] / mw_i(spf);
    }
  }
};

template <typename LiqPropType>
struct InitLiqProps
{
  void
  operator()(LiqPropType* /*ldata*/, std::vector<std::string>& /*fuel_names*/)
  {
    amrex::Abort("InitLiqProps not implemented for this LiqPropType.");
  }
};

} // namespace pele::physics::SprayProps

#endif
