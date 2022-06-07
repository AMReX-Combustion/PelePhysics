#ifndef _MANIFOLD_FUNC_H_
#define _MANIFOLD_FUNC_H_

namespace pele {
namespace physics {

class ManFuncParams
{
public:

  ManifoldFunc() = delete;

  virtual void initialize();
  virtual void allocate();
  virtual void deallocate();
  virtual void sync_to_device();

}; // class ManFuncParams

class ManifoldFunc
{
public:

  ManifoldFunc() = delete;
  virtual void get_value(const int ivar, amrex::Real indata[], amrex::Real& out);
  virtual void get_derivs(const int ivar, amrex::Real indata[], amrex::Real derivs[]);

}; // class ManifoldFunc
}  // namespace physics
}  // namespace pele

#endif
