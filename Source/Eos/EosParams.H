#ifndef EOS_PARAMS_H
#define EOS_PARAMS_H

#include <AMReX_REAL.H>
#include <AMReX_ParmParse.H>

#include "Table.H"

// TODO: Single general interface for eosparm, pmfdata, transparm, etc.

namespace pele {
namespace physics {
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

template <typename EOSType>
struct InitEosParm
{
  void operator()(EosParm<EOSType>* /*eosparm*/) {}
};

template <>
struct EosParm<GammaLaw>
{
  amrex::Real gamma{Constants::gamma};
};

template <>
struct InitEosParm<GammaLaw>
{
  void operator()(EosParm<GammaLaw>* eosparm)
  {
    amrex::ParmParse pp("eos");
    pp.query("gamma", eosparm->gamma);
    if (eosparm->gamma != Constants::gamma) {
      amrex::Warning(
        "Runtime gamma does not match compile time gamma, use caution because "
        "runtime specification of gamma may not yet be supported everywhere");
    }
  }
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
  bool has_mani_src = false;
};

template <>
struct InitEosParm<Manifold>
{
  void operator()(EosParm<Manifold>* eosparm)
  {
    amrex::Error("Manifold EOS cannot be initialized without tabulated or "
                 "neural network data");
  }

  void operator()(
    EosParm<Manifold>* eosparm, ManFuncParams::ManFuncData* manf_data_in)
  {
    eosparm->manf_data = manf_data_in;
    // First (N-1) species are table dimensions (last species corresponds to
    // density)
    AMREX_ALWAYS_ASSERT(manf_data_in->Ndim == NUM_SPECIES - 1);

    amrex::ParmParse pp("eos");
    pp.get("nominal_pressure_cgs", eosparm->Pnom_cgs);
    pp.query("has_mani_src", eosparm->has_mani_src);
    pp.get("compute_temperature", eosparm->compute_temperature);

    // Setup density lookups
    std::string density_lookup_type_string{"linear"};
    pp.query("density_lookup_type", density_lookup_type_string);
    if (density_lookup_type_string == "linear") {
      eosparm->dens_lookup = density_lookup_type::linear;
      eosparm->idx_density = get_var_index("RHO", manf_data_in);
      amrex::Print() << "Manifold EOS: Using linear density lookups : index = "
                     << eosparm->idx_density << std::endl;
    } else if (density_lookup_type_string == "log") {
      eosparm->dens_lookup = density_lookup_type::log;
      eosparm->idx_density = get_var_index("lnRHO", manf_data_in);
      amrex::Print()
        << "Manifold EOS: Using logarithmic density lookups : index = "
        << eosparm->idx_density << std::endl;
    } else if (density_lookup_type_string == "inverse") {
      eosparm->dens_lookup = density_lookup_type::inverse;
      eosparm->idx_density = get_var_index("invRHO", manf_data_in);
      amrex::Print() << "Manifold EOS: Using inverse density lookups : index = "
                     << eosparm->idx_density << std::endl;
    } else {
      amrex::Abort("Invalid density lookup type supplied");
    }

    // Get important indices
    eosparm->idx_T = get_var_index("T", manf_data_in);

    // For manifold table parameter source terms, assume if index not found,
    // source term is 0 For neural net, require a definition to be supplied for
    // each manifold parameter
    if (eosparm->has_mani_src) {
      for (int idim = 0; idim < manf_data_in->Ndim; idim++) {
        std::string dimname = std::string(
          &manf_data_in->dimnames[idim * manf_data_in->len_str],
          manf_data_in->len_str);
        std::string dim_src;
        if (dimname.rfind("Y-") == 0) {
          dim_src = "SRC_" + amrex::trim(dimname).substr(2u, std::string::npos);
        } else {
          dim_src = "SRC_" + amrex::trim(dimname);
        }
        eosparm->idx_Wdot[idim] =
          get_var_index(dim_src.c_str(), manf_data_in, false);
        if (eosparm->idx_Wdot[idim] < 0) {
          amrex::Print()
            << "Warning: No source term found for manifold parameter "
            << amrex::trim(dimname) << ", assuming SRC_" << amrex::trim(dimname)
            << " = 0" << std::endl;
        }
      }
    } else {
      for (int idim = 0; idim < manf_data_in->Ndim; idim++) {
        // The info file / metadata loader will ensure we have a definition for
        // each
        eosparm->idx_Wdot[idim] = -1;
      }
    }
  }
};

template <typename EOSType>
class EosParams
{
public:
  EosParams() {}

  ~EosParams() {}

  void allocate()
  {
    InitEosParm<EOSType>()(&m_h_eos_parm);
    if (!m_device_allocated) {
      m_d_eos_parm = (EosParm<EOSType>*)amrex::The_Device_Arena()->alloc(
        sizeof(m_h_eos_parm));
      m_device_allocated = true;
      sync_to_device();
    }
  }

  void allocate(ManFuncParams::ManFuncData* manf_data)
  {
    InitEosParm<EOSType>()(&m_h_eos_parm, manf_data);
    if (!m_device_allocated) {
      m_d_eos_parm = (EosParm<EOSType>*)amrex::The_Device_Arena()->alloc(
        sizeof(m_h_eos_parm));
      m_device_allocated = true;
      sync_to_device();
    }
  }

  void deallocate()
  {
    if (m_device_allocated)
      amrex::The_Device_Arena()->free(m_d_eos_parm);
  }

  void sync_to_device()
  {
    if (!m_device_allocated) {
      amrex::Abort("Device params not allocated yet");
    } else {
      amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, &m_h_eos_parm, &m_h_eos_parm + 1,
        m_d_eos_parm);
    }
  }

  EosParm<EOSType>& host_eos_parm() { return m_h_eos_parm; }

  const EosParm<EOSType>* device_eos_parm() { return m_d_eos_parm; }

private:
  EosParm<EOSType> m_h_eos_parm;
  EosParm<EOSType>* m_d_eos_parm;
  bool m_device_allocated{false};
};

} // namespace eos
} // namespace physics
} // namespace pele
#endif
