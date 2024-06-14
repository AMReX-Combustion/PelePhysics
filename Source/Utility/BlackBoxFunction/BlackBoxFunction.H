#ifndef BLACK_BOX_FUNC_H
#define BLACK_BOX_FUNC_H

#include <AMReX_Utility.H>
#include "PeleParamsGeneric.H"

namespace pele::physics {

enum class BlackBoxModel { TABLE = 0, NEURAL_NET = 1 };

struct BlackBoxFunctionData
{
  BlackBoxModel bbmodel;
  static constexpr int len_str = 64; // Length of strings in the char arrays
  char model_name[len_str + 1];      // Name of the neural network model
  int Ndim;                          // Input dimension (# of manifold params)
  int Nvar;       // Output dimension (# of quantities of interest)
  int Nmanpar;    // # of manifold parameters (discounting variables passed to
                  // input layer)
  char* dimnames; // Names of the manifold parameters (including the passed
                  // variables)
  char* varnames; // Names of the output variables
};

template <>
struct HostOnlyParm<BlackBoxFunctionData>
{
  std::string parm_parse_prefix;
};

AMREX_FORCE_INLINE
int
get_var_index(
  const char* varname,
  const BlackBoxFunctionData* manf_data,
  const bool require = true)
{
  std::string var(varname);
  int ii;
  for (ii = 0; ii < manf_data->Nvar; ii++) {
    std::string vartest = std::string(
      &(manf_data->varnames)[ii * manf_data->len_str], manf_data->len_str);
    if (amrex::trim(var) == amrex::trim(vartest)) {
      break;
    }
  }
  if (ii == manf_data->Nvar && require) {
    amrex::Abort("Requested variable not found in table: " + var);
  } else if (ii == manf_data->Nvar) {
    ii = -1; // less than 0 = not found
  }
  return ii;
}

class BlackBoxFunction
{
public:
  AMREX_GPU_HOST_DEVICE
  virtual void
  get_value(const int ivar, const amrex::Real indata[], amrex::Real& out) = 0;

  AMREX_GPU_HOST_DEVICE
  virtual void get_values(
    const int nvar,
    const int ivar[],
    const amrex::Real indata[],
    amrex::Real out[]) = 0;

  AMREX_GPU_HOST_DEVICE
  virtual void
  get_all_values(const amrex::Real indata[], amrex::Real out[]) = 0;

  AMREX_GPU_HOST_DEVICE
  virtual void get_derivs(
    const int ivar, const amrex::Real indata[], amrex::Real derivs[]) = 0;

  AMREX_GPU_HOST_DEVICE
  virtual BlackBoxModel model() = 0;

  AMREX_GPU_HOST_DEVICE
  virtual ~BlackBoxFunction() = default;

}; // class BlackBoxFunction

} // namespace pele::physics
#endif
