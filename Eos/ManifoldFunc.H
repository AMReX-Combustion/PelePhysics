#ifndef _MANIFOLD_FUNC_H_
#define _MANIFOLD_FUNC_H_

namespace pele {
namespace physics {

struct ManFuncData
{
  static constexpr int len_str = 64;
  int Ndim;
  int Nvar;
  char model_name[len_str+1];
  char* varnames;
  char* dimnames;
};

class ManFuncParams
{
public:

  virtual void initialize();
  virtual void allocate();
  virtual void deallocate();
  virtual void sync_to_device();

}; // class ManFuncParams

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
