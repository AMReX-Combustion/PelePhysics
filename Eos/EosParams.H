#ifndef EOS_PARAMS_H
#define EOS_PARAMS_H

#include <AMReX_REAL.H>
#include <AMReX_ParmParse.H>

#include "Table.H"

//TODO: Single general interface for eosparm, pmfdata, transparm, etc.

namespace pele {
namespace physics {
namespace eos {

AMREX_FORCE_INLINE
int get_var_index(const char* varname, const char* varlist, const int nvar, const int len_str, const bool require = true)
{
  std::string var(varname);
  int ii;
  for (ii = 0; ii < nvar; ii++) {
    std::string vartest = std::string(&varlist[ii*len_str], len_str);
    if (amrex::trim(var) == amrex::trim(vartest)) break;
  }
  if (ii == nvar && require) {
    amrex::Abort("Requested variable not found in table: " + var);
  }
  else if (ii == nvar) {
    ii = -1; // less than 0 = not found
  }
  return ii;
}
  
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
  amrex::Real gamma{1.4};
};

template <>
struct InitEosParm<GammaLaw>
{
  void operator()(EosParm<GammaLaw>* eosparm)
  {
    amrex::ParmParse pp("eos");
    pp.query("gamma", eosparm->gamma);
    //TODO: Actually integrate this into gamma law EOS code
    amrex::Error("Gamma cannot yet be specified at runtime");
  }
};

template <>
struct EosParm<Manifold>
{
  TabFuncParams::TabFuncData* tf_data;
  amrex::Real Pnom_cgs;
  int idx_density{0};
  int idx_T{0};
  int idx_Wdot[NUM_SPECIES - 1] = {0};
};

template <>
struct InitEosParm<Manifold>
{
  void operator()(EosParm<Manifold>* eosparm)
  {
    amrex::Error("Manifold EOS cannot be initialized without tabulated data");
  }

  void operator()(EosParm<Manifold>* eosparm, TabFuncParams::TabFuncData* tf_data_in)
  {
    eosparm->tf_data = tf_data_in;
    // First (N-1) species are table dimensions (last species corresponds to density)
    AMREX_ALWAYS_ASSERT(tf_data_in->Ndim == NUM_SPECIES - 1);
      
    amrex::ParmParse pp("eos");
    pp.get("nominal_pressure_cgs", eosparm->Pnom_cgs);
      
    // Get important indices
    eosparm->idx_density = get_var_index("RHO", tf_data_in->varnames, tf_data_in->Nvar, tf_data_in->len_str);
    eosparm->idx_T = get_var_index("T", tf_data_in->varnames, tf_data_in->Nvar, tf_data_in->len_str);

    // For manifold parameter source terms, assume if index not found, source term is 0
    for (int idim = 0; idim < tf_data_in->Ndim; idim++) {
      std::string dimname = std::string(&tf_data_in->dimnames[idim*tf_data_in->len_str], tf_data_in->len_str);
      eosparm->idx_Wdot[idim] = get_var_index(("SRC_"+amrex::trim(dimname)).c_str(),
					       tf_data_in->varnames, tf_data_in->Nvar, tf_data_in->len_str, false);
      if( eosparm->idx_Wdot[idim] < 0) {
	amrex::Print() << "Warning: No source term found for manifold parameter " << tf_data_in->dimnames[idim]
		       << ", assuming SRC_" << tf_data_in->dimnames[idim] << " = 0";
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
      m_d_eos_parm =
        (EosParm<EOSType>*)amrex::The_Device_Arena()->alloc(
          sizeof(m_h_eos_parm));
      m_device_allocated = true;
      sync_to_device();
    }
  }

  void allocate(TabFuncParams::TabFuncData* tf_data)
  {
    InitEosParm<EOSType>()(&m_h_eos_parm, tf_data);
    if (!m_device_allocated) {
      m_d_eos_parm =
        (EosParm<EOSType>*)amrex::The_Device_Arena()->alloc(
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

  EosParm<EOSType>& host_eos_parm()
  {
    return m_h_eos_parm;
  }
  
  const EosParm<EOSType>* device_eos_parm()
  {
    return m_d_eos_parm;
  }

private:
  EosParm<EOSType> m_h_eos_parm;
  EosParm<EOSType>* m_d_eos_parm;
  bool m_device_allocated{false};
};

} // namespace eos
} // namespace physics
} // namespace pele
#endif
