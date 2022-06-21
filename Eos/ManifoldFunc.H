#ifndef _MANIFOLD_FUNC_H_
#define _MANIFOLD_FUNC_H_

namespace pele {
namespace physics {
  
enum class ManifoldModel
{
  TABLE = 0,
  NEURAL_NET = 1
};

class ManFuncParams
{
public:

  virtual void initialize() = 0;
  virtual void allocate() = 0;
  virtual void deallocate() = 0;
  virtual void sync_to_device() = 0;
  
  struct ManFuncData
  {
    ManifoldModel manmodel;
    static constexpr int len_str = 64; // Length of strings in the char arrays
    char model_name[len_str+1]; // Name of the neural network model
    int Ndim; // Input dimension (# of manifold params)
    int Nvar; // Output dimension (# of quantities of interest)
    int Nmanpar; // # of manifold parameters (discounting variables passed to input layer)
    char* dimnames; // Names of the manifold parameters (including the passed variables)
    char* varnames; // Names of the output variables
  };
  
  virtual ManFuncData& host_manfunc_data() = 0;
  virtual ManFuncData* device_manfunc_data() = 0;
  
  virtual ~ManFuncParams() {}

}; // class ManFuncParams

AMREX_FORCE_INLINE
int get_var_index(const char* varname, const ManFuncParams::ManFuncData* manf_data,
    const bool require = true)
{      
  std::string var(varname);
  int ii;
  for (ii = 0; ii < manf_data->Nvar; ii++) {
    std::string vartest = std::string(&(manf_data->varnames)[ii*manf_data->len_str], manf_data->len_str);
    if (amrex::trim(var) == amrex::trim(vartest)) break;
  }
  if (ii == manf_data->Nvar && require) {
    amrex::Abort("Requested variable not found in table: " + var);
  }
  else if (ii == manf_data->Nvar) {
    ii = -1; // less than 0 = not found
  }
  return ii;
}

class ManifoldFunc
{
public:

  AMREX_GPU_HOST_DEVICE
  virtual void get_value(const int ivar, const amrex::Real indata[], amrex::Real& out) = 0;

  AMREX_GPU_HOST_DEVICE
  virtual void get_derivs(const int ivar, const amrex::Real indata[], amrex::Real derivs[]) = 0;
  
  AMREX_GPU_HOST_DEVICE
  virtual void calculate_Wdot(const int paramidx, const amrex::Real indata[], amrex::Real& out) = 0;
  
  AMREX_GPU_HOST_DEVICE
  virtual ManifoldModel model() = 0;
  
  virtual ~ManifoldFunc() {}
  
}; // class ManifoldFunc
}  // namespace physics
}  // namespace pele

#endif
