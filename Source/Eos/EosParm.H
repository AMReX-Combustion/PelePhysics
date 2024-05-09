#ifndef EOS_PARAMS_H
#define EOS_PARAMS_H

#include <AMReX_REAL.H>
#include <AMReX_ParmParse.H>

#include "mechanism.H"
#include "ManifoldFunc.H"
#include "Table.H"
#include "NeuralNetHomerolled.H"
#include "PhysicsConstants.H"
#include "PeleParamsGeneric.H"

// TODO: Single general interface for eosparm, pmfdata, transparm, etc.


namespace pele {
namespace physics {
namespace eos{

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
  
template <>
struct EosParm<Manifold>
{
  ManFuncParams::ManFuncData* manf_data;
  amrex::Real Pnom_cgs;
  density_lookup_type dens_lookup;
  bool compute_temperature{false};
  int idx_density{0};
  int idx_T{0};
  int idx_Wdot[NUM_SPECIES - 1];
  bool has_mani_src = true;
};
  
} // namespace eos

  template<>
  struct HostOnlyParm<eos::EosParm<eos::Manifold>>
  {
    std::unique_ptr<ManFuncParams> manfunc_par;
  };

  
  template<>
void PeleParams<eos::EosParm<eos::GammaLaw>>::host_initialize()
{
  amrex::ParmParse pp("eos");
  pp.query("gamma", m_h_parm.gamma);
  if (m_h_parm.gamma != Constants::gamma) {
    amrex::Warning(
                   "Runtime gamma does not match compile time gamma, use caution because "
                   "runtime specification of gamma may not yet be supported everywhere");
  }
};
  
template<>
void PeleParams<eos::EosParm<eos::Manifold>>::host_initialize()
{ 
  amrex::ParmParse ppm("manifold");
  std::string manifold_model;
  ppm.get("model", manifold_model);
  if (manifold_model == "Table") {
    m_host_only_parm.manfunc_par.reset(new TabFuncParams());
    amrex::Print() << " Initialization of Table (CPP)... \n";
  } else if (manifold_model == "NeuralNet") {
    m_host_only_parm.manfunc_par.reset(new NNFuncParams());
    amrex::Print() << " Initialization of Neural Net Func. (CPP)... \n";
  } else {
    amrex::Error("Invalid manifold model!");
  }
  m_host_only_parm.manfunc_par->initialize();
  
  ManFuncParams::ManFuncData* d_manf_data_in = m_host_only_parm.manfunc_par->device_manfunc_data();
  ManFuncParams::ManFuncData* h_manf_data_in = &(m_host_only_parm.manfunc_par->host_manfunc_data());
    
    m_h_parm.manf_data = d_manf_data_in;
    // First (N-1) species are table dimensions (last species corresponds to
    // density)
    AMREX_ALWAYS_ASSERT(h_manf_data_in->Ndim == NUM_SPECIES - 1);

    amrex::ParmParse pp("eos");
    pp.get("nominal_pressure_cgs", m_h_parm.Pnom_cgs);
    pp.query("has_mani_src", m_h_parm.has_mani_src);
    pp.get("compute_temperature", m_h_parm.compute_temperature);

    // Setup density lookups
    std::string density_lookup_type_string{"linear"};
    pp.query("density_lookup_type", density_lookup_type_string);
    if (density_lookup_type_string == "linear") {
      m_h_parm.dens_lookup = eos::density_lookup_type::linear;
      m_h_parm.idx_density = get_var_index("RHO", h_manf_data_in);
      amrex::Print() << "Manifold EOS: Using linear density lookups : index = "
                     << m_h_parm.idx_density << std::endl;
    } else if (density_lookup_type_string == "log") {
      m_h_parm.dens_lookup = eos::density_lookup_type::log;
      m_h_parm.idx_density = get_var_index("lnRHO", h_manf_data_in);
      amrex::Print()
        << "Manifold EOS: Using logarithmic density lookups : index = "
        << m_h_parm.idx_density << std::endl;
    } else if (density_lookup_type_string == "inverse") {
      m_h_parm.dens_lookup = eos::density_lookup_type::inverse;
      m_h_parm.idx_density = get_var_index("invRHO", h_manf_data_in);
      amrex::Print() << "Manifold EOS: Using inverse density lookups : index = "
                     << m_h_parm.idx_density << std::endl;
    } else {
      amrex::Abort("Invalid density lookup type supplied");
    }

    // Get important indices
    m_h_parm.idx_T = get_var_index("T", h_manf_data_in);

    // For manifold table parameter source terms, assume if index not found,
    // source term is 0 For neural net, require a definition to be supplied for
    // each manifold parameter
    if (m_h_parm.has_mani_src) {
      for (int idim = 0; idim < h_manf_data_in->Ndim; idim++) {
        std::string dimname = std::string(
          &h_manf_data_in->dimnames[idim * h_manf_data_in->len_str],
          h_manf_data_in->len_str);
        std::string dim_src;
        if (dimname.rfind("Y-") == 0) {
          dim_src = "SRC_" + amrex::trim(dimname).substr(2u, std::string::npos);
        } else {
          dim_src = "SRC_" + amrex::trim(dimname);
        }
        m_h_parm.idx_Wdot[idim] =
          get_var_index(dim_src.c_str(), h_manf_data_in, false);
        if (m_h_parm.idx_Wdot[idim] < 0) {
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
        m_h_parm.idx_Wdot[idim] = -1;
        amrex::Abort("Computation of manifold src terms from species src terms "
                     "needs to be reimplemented");
      }
    }
};

template<>
void PeleParams<eos::EosParm<eos::Manifold>>::host_deallocate () {
  m_host_only_parm.manfunc_par->deallocate();
};
  
} // namespace physics
} // namespace pele
#endif
