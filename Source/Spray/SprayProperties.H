#ifndef SPRAYPROPERTIES_H
#define SPRAYPROPERTIES_H

using namespace amrex;

// Structs for liquid physical data
namespace pele::physics{

  namespace SprayProps{

    // Original spray properties ported from PeleMP to PelePhysics
    struct MPLiqProps
    {
      amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> critT;
      amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> boilT;
      amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> cp;
      amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> latent;
      amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> ref_latent;
      // 3 coefficients for Antoine equation and conversion to appropriate units
      amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM * 4> psat_coef;
      amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM * 4> rho_coef;
      amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM * 4> lambda_coef;
      amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM * 4> mu_coef;
      amrex::GpuArray<int, SPRAY_FUEL_NUM> indx = {{-1}};
      amrex::GpuArray<int, SPRAY_FUEL_NUM> dep_indx = {{-1}};

      
      AMREX_GPU_HOST_DEVICE
      AMREX_FORCE_INLINE
      amrex::Real rhoL(const amrex::Real& T, const int spf) const
      {
        amrex::Real a = rho_coef[4 * spf];
        amrex::Real b = rho_coef[4 * spf + 1];
        amrex::Real c = rho_coef[4 * spf + 2];
        amrex::Real d = rho_coef[4 * spf + 3];
        return a + T * (b + T * (c + T * d));
      }

      
      AMREX_GPU_HOST_DEVICE
      AMREX_FORCE_INLINE
      amrex::Real lambdaL(const amrex::Real& T, const int spf) const
      {
        amrex::Real a = lambda_coef[4 * spf];
        amrex::Real b = lambda_coef[4 * spf + 1];
        amrex::Real c = lambda_coef[4 * spf + 2];
        amrex::Real d = lambda_coef[4 * spf + 3];
        return a + T * (b + T * (c + T * d));
      }

      AMREX_GPU_HOST_DEVICE
      AMREX_FORCE_INLINE
      amrex::Real muL(const amrex::Real& T, const int spf) const
      {
        amrex::Real a = mu_coef[4 * spf];
        amrex::Real b = mu_coef[4 * spf + 1];
        amrex::Real c = mu_coef[4 * spf + 2];
        amrex::Real d = mu_coef[4 * spf + 3];
        return a + ((d / T + c) / T + b) / T;
      }

      AMREX_GPU_HOST_DEVICE
      AMREX_FORCE_INLINE
      amrex::Real psat(const amrex::Real& T, const int spf) const
      {
        amrex::Real a = psat_coef[4 * spf];
        amrex::Real b = psat_coef[4 * spf + 1];
        amrex::Real c = psat_coef[4 * spf + 2];
        amrex::Real d = psat_coef[4 * spf + 3];
        return d * std::pow(10., a - b / (T + c));
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
    struct InitLiqProps {
        void operator()(LiqPropType* ldata, std::vector<std::string>& fuel_names) {
            amrex::Abort("InitLiqProps not implemented for this LiqPropType.");
        }
    };

    // Specialization for MPLiqProps
    template <>
    struct InitLiqProps<MPLiqProps> {
        void operator()(MPLiqProps* ldata, std::vector<std::string>& fuel_names) {
            ParmParse pp("particles");
            getInpVal(ldata->critT.data(), pp, fuel_names.data(), "crit_temp");
            getInpVal(ldata->boilT.data(), pp, fuel_names.data(), "boil_temp");
            getInpVal(ldata->cp.data(), pp, fuel_names.data(), "cp");
            getInpVal(ldata->ref_latent.data(), pp, fuel_names.data(), "latent");

            getInpCoef(
              ldata->lambda_coef.data(), pp, fuel_names.data(), "lambda");
            getInpCoef(ldata->psat_coef.data(), pp, fuel_names.data(), "psat");
            getInpCoef(
              ldata->rho_coef.data(), pp, fuel_names.data(), "rho", true);
            getInpCoef(ldata->mu_coef.data(), pp, fuel_names.data(), "mu");
        }

        void getInpCoef(
            Real* coef,
            const amrex::ParmParse& ppp,
            const std::string* fuel_names,
            const std::string& varname,
            bool is_required = false) {
            for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
                std::string var_read = fuel_names[spf] + "_" + varname;
                int numvals = ppp.countval(var_read.c_str());
                if (numvals == 4) {
                    std::vector<Real> inp_coef(4, 0.);
                    if (is_required) {
                        ppp.getarr(var_read.c_str(), inp_coef);
                    } else {
                        ppp.queryarr(var_read.c_str(), inp_coef);
                    }
                    for (int i = 0; i < 4; ++i) {
                        coef[4 * spf + i] = inp_coef[i];
                    }
                } else if (numvals == 1) {
                    Real inp_coef = 0.;
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

        void getInpVal(
            Real* coef,
            const amrex::ParmParse& ppp,
            const std::string* fuel_names,
            const std::string& varname) {
            for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
                std::string var_read = fuel_names[spf] + "_" + varname;
                ppp.get(var_read.c_str(), coef[spf]);
            }
        }
    };

    /*
    // Specialization for GCMLiqProps
    template <>
    struct InitLiqProps<GCMLiqProps> {
        void operator()(GCMLiqProps* ldata, std::vector<std::string>& fuel_names) {
            amrex::ParmParse pp("particles");

            
        }
    };
    */
  } // namespace SprayProps
} // namespace pele::physics

#endif