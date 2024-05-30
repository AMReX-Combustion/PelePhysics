#ifndef PELEPARAMSGENERIC_H
#define PELEPARAMSGENERIC_H

#include <AMReX_Gpu.H>
#include "Factory.H"

namespace pele::physics {

// ParmType: Generic parameters to be copied from host to device
// Template for HostOnlyParm: Data kept only by Host
template <typename ParmType>
struct HostOnlyParm
{
};

// Forward Declaration of PeleParams
template <typename ParmType, typename BaseParmType = ParmType>
class PeleParams;

// Template for initializer/destructor functions
template <typename ParmType, typename BaseParmType = ParmType>
struct InitParm
{
  static void host_initialize(PeleParams<ParmType, BaseParmType>*) {}

  static void host_deallocate(PeleParams<ParmType, BaseParmType>*) {}
};

// Define generic interface
// This gets a bit complicated to handle the case of base and derived ParmType
// classes
template <typename BaseParmType>
class PeleParamsGeneric : public Factory<PeleParamsGeneric<BaseParmType>>
{
public:
  static std::string base_identifier() { return "pele_params_base_generic"; }

  virtual void initialize() = 0;
  virtual void device_allocate() = 0;
  virtual void sync_to_device() = 0;
  virtual void deallocate() = 0;
  virtual BaseParmType& host_parm() = 0;
  virtual const BaseParmType* device_parm() = 0;
  virtual HostOnlyParm<BaseParmType>& host_only_parm() = 0;
};

// Template for class that allocates/deallocates a generic ParmType
// For most use cases ParmType doesn't have inheritance to worry about
// and uses the default BaseParmType = ParmType defined in the forward
// declaration above
template <typename ParmType, typename BaseParmType>
class PeleParams : public PeleParamsGeneric<BaseParmType>
{

  friend struct InitParm<ParmType, BaseParmType>;

public:
  static_assert(std::is_base_of_v<BaseParmType, ParmType> == true);

  PeleParams() {}

  virtual ~PeleParams() {}

  void initialize() override
  {
    InitParm<ParmType, BaseParmType>::host_initialize(this);
    device_allocate();
  }

  void device_allocate() override
  {
    if (!m_device_allocated) {
      m_d_parm = (ParmType*)amrex::The_Device_Arena()->alloc(sizeof(m_h_parm));
      m_device_allocated = true;
      sync_to_device();
    }
  }

  void sync_to_device() override
  {
    if (!m_device_allocated) {
      amrex::Abort("Device params not allocated yet");
    } else {
      amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, &m_h_parm, &m_h_parm + 1, m_d_parm);
    }
  }

  void deallocate() override
  {
    InitParm<ParmType, BaseParmType>::host_deallocate(this);
    if (m_device_allocated)
      amrex::The_Device_Arena()->free(m_d_parm);
  }

  BaseParmType& host_parm() override { return m_h_parm; }

  const BaseParmType* device_parm() override { return m_d_parm; }

  HostOnlyParm<BaseParmType>& host_only_parm() override
  {
    return m_host_only_parm;
  }

private:
  HostOnlyParm<BaseParmType> m_host_only_parm;
  ParmType m_h_parm;
  ParmType* m_d_parm;
  bool m_device_allocated{false};
};

} // namespace pele::physics
#endif
