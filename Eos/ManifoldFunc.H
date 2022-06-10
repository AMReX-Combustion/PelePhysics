#ifndef _MANIFOLD_FUNC_H_
#define _MANIFOLD_FUNC_H_

namespace pele {
namespace physics {

class ManFuncParams
{
public:

  virtual void initialize();
  virtual void allocate();
  virtual void deallocate();
  virtual void sync_to_device();
  
  struct ManFuncData
  {
    static constexpr int len_str = 64;
    int Ndim; // Input dimension (# of manifold params)
    int Nvar; // Output dimension (# of quantities of interest)
    char model_name[len_str+1];
    char* varnames;
    char* dimnames;
  };

}; // class ManFuncParams

AMREX_FORCE_INLINE
int get_var_index(const char* varname, const ManFuncParams::ManFuncData* manf_data,
    const bool require = true)
{      
  std::string var(varname);
  int ii;
  for (ii = 0; ii < manf_data->Nvar; ii++) {
    std::string vartest = std::string(&(manf_data->varnames)[ii*manf_data->len_str],
        manf_data->len_str);
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
  AMREX_FORCE_INLINE
  virtual void get_value(const int ivar, const amrex::Real indata[], amrex::Real& out);

  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  virtual void get_derivs(const int ivar, const amrex::Real indata[], amrex::Real derivs[]);

}; // class ManifoldFunc
}  // namespace physics
}  // namespace pele

#endif
