#ifndef EOS_PARAMS_H
#define EOS_PARAMS_H

#include <AMReX_REAL.H>
#include <AMReX_ParmParse.H>

#include "mechanism.H"
#include "BlackBoxFunction.H"
#include "PhysicsConstants.H"
#include "PeleParamsGeneric.H"

namespace pele::physics {
namespace eos {

enum class density_lookup_type { linear, log, inverse };

struct GammaLaw;
struct Fuego;
struct SRK;
struct Manifold;

template <typename EOSType>
struct EosParm
{
};

template <>
struct EosParm<GammaLaw>
{
  amrex::Real gamma{Constants::gamma};
};

#ifndef MANIFOLD_DIM
#define MANIFOLD_DIM 1
#endif

#ifndef AMREX_USE_SYCL
template <>
struct EosParm<Manifold>
{
  const BlackBoxFunctionData* manf_data;
  amrex::Real Pnom_cgs;
  density_lookup_type dens_lookup;
  bool compute_temperature{false};
  int idx_density{0};
  int idx_T{0};
  int idx_Wdot[MANIFOLD_DIM];
  bool has_mani_src{true};
};
#endif

} // namespace eos

template <>
struct InitParm<eos::EosParm<eos::GammaLaw>>
{
  static void host_initialize(PeleParams<eos::EosParm<eos::GammaLaw>>* parm_in)
  {
    amrex::ParmParse pp("eos");
    pp.query("gamma", parm_in->m_h_parm.gamma);
    if (parm_in->m_h_parm.gamma != Constants::gamma) {
      amrex::Warning(
        "Runtime gamma does not match compile time gamma, use caution because "
        "runtime specification of gamma may not yet be supported everywhere");
    }
  }

  static void
  host_deallocate(PeleParams<eos::EosParm<eos::GammaLaw>>* /*parm_in*/)
  {
  }
};

#ifndef AMREX_USE_SYCL
template <>
struct HostOnlyParm<eos::EosParm<eos::Manifold>>
{
  std::shared_ptr<PeleParamsGeneric<BlackBoxFunctionData>> manfunc_par;
};

template <>
struct InitParm<eos::EosParm<eos::Manifold>>
{
  static void host_initialize(PeleParams<eos::EosParm<eos::Manifold>>* parm_in)
  {
    amrex::ParmParse pp("manifold");
    std::string manifold_model;
    pp.get("model", manifold_model);
    parm_in->m_host_only_parm.manfunc_par =
      PeleParamsGeneric<BlackBoxFunctionData>::create(manifold_model);
    parm_in->m_host_only_parm.manfunc_par->host_only_parm().parm_parse_prefix =
      "manifold";
    parm_in->m_host_only_parm.manfunc_par->initialize();
    const BlackBoxFunctionData* d_manf_data_in =
      parm_in->m_host_only_parm.manfunc_par->device_parm();
    BlackBoxFunctionData* h_manf_data_in =
      &(parm_in->m_host_only_parm.manfunc_par->host_parm());

    parm_in->m_h_parm.manf_data = d_manf_data_in;
    // First (N-1) species are table dimensions (last species corresponds to
    // density)
    AMREX_ALWAYS_ASSERT(h_manf_data_in->Ndim == NUM_SPECIES - 1);

    pp.get("nominal_pressure_cgs", parm_in->m_h_parm.Pnom_cgs);
    pp.query("has_mani_src", parm_in->m_h_parm.has_mani_src);
    pp.get("compute_temperature", parm_in->m_h_parm.compute_temperature);

    // Setup density lookups
    std::string density_lookup_type_string{"linear"};
    pp.query("density_lookup_type", density_lookup_type_string);
    if (density_lookup_type_string == "linear") {
      parm_in->m_h_parm.dens_lookup = eos::density_lookup_type::linear;
      parm_in->m_h_parm.idx_density = get_var_index("RHO", h_manf_data_in);
      amrex::Print() << "Manifold EOS: Using linear density lookups : index = "
                     << parm_in->m_h_parm.idx_density << std::endl;
    } else if (density_lookup_type_string == "log") {
      parm_in->m_h_parm.dens_lookup = eos::density_lookup_type::log;
      parm_in->m_h_parm.idx_density = get_var_index("lnRHO", h_manf_data_in);
      amrex::Print()
        << "Manifold EOS: Using logarithmic density lookups : index = "
        << parm_in->m_h_parm.idx_density << std::endl;
    } else if (density_lookup_type_string == "inverse") {
      parm_in->m_h_parm.dens_lookup = eos::density_lookup_type::inverse;
      parm_in->m_h_parm.idx_density = get_var_index("invRHO", h_manf_data_in);
      amrex::Print() << "Manifold EOS: Using inverse density lookups : index = "
                     << parm_in->m_h_parm.idx_density << std::endl;
    } else {
      amrex::Abort("Invalid density lookup type supplied");
    }

    // Get important indices
    parm_in->m_h_parm.idx_T = get_var_index("T", h_manf_data_in);

    // For manifold table parameter source terms, assume if index not found,
    // source term is 0 For neural net, require a definition to be supplied for
    // each manifold parameter
    if (parm_in->m_h_parm.has_mani_src) {
      for (int idim = 0; idim < h_manf_data_in->Ndim; idim++) {
        const std::string dimname = std::string(
          &h_manf_data_in->dimnames[idim * h_manf_data_in->len_str],
          h_manf_data_in->len_str);
        const std::string dim_src =
          (dimname.rfind("Y-") == 0)
            ? "SRC_" + amrex::trim(dimname).substr(2U, std::string::npos)
            : "SRC_" + amrex::trim(dimname);
        parm_in->m_h_parm.idx_Wdot[idim] =
          get_var_index(dim_src.c_str(), h_manf_data_in, false);
        if (parm_in->m_h_parm.idx_Wdot[idim] < 0) {
          amrex::Print()
            << "Warning: No source term found for manifold parameter "
            << amrex::trim(dimname) << ", assuming SRC_" << amrex::trim(dimname)
            << " = 0" << std::endl;
        }
      }
    } else {
      for (int idim = 0; idim < h_manf_data_in->Ndim; idim++) {
        // The info file / metadata loader will ensure we have a definition for
        // each
        parm_in->m_h_parm.idx_Wdot[idim] = -1;
        amrex::Abort("Computation of manifold src terms from species src terms "
                     "needs to be reimplemented");
      }
    }
  }

  static void host_deallocate(PeleParams<eos::EosParm<eos::Manifold>>* parm_in)
  {
    parm_in->m_host_only_parm.manfunc_par->deallocate();
    parm_in->m_host_only_parm.manfunc_par = nullptr;
  }
};
#endif

} // namespace pele::physics
#endif
