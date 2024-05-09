#ifndef PELEPARAMSGENERIC_H
#define PELEPARAMSGENERIC_H

#include <AMReX_Gpu.H>
#include "Factory.H"

namespace pele::physics {

template <typename ParmType>
struct HostOnlyParm {};

template <typename ParmType>
class PeleParams;
  
template <typename ParmType>
struct InitParm {
  static void host_initialize(PeleParams<ParmType>*) {}

  static void host_deallocate(PeleParams<ParmType>*) {}
};
  
template <typename ParmType>
class PeleParams : public Factory<PeleParams<ParmType>>
{

  friend class InitParm<ParmType>;
  
public: 
  static std::string base_identifier() { return "pele_params_base_generic";}

  PeleParams() {}

  virtual ~PeleParams() {}

  virtual void initialize()
  {
    std::cout << "starting initialization" << std::endl;
    InitParm<ParmType>::host_initialize(this);
    device_allocate();
  }
  
  virtual void device_allocate()
  {
    if (!m_device_allocated) {
      m_d_parm = (ParmType*)amrex::The_Device_Arena()->alloc(
        sizeof(m_h_parm));
      m_device_allocated = true;
      sync_to_device();
    }
  }
  
  virtual void sync_to_device()
  {
    if (!m_device_allocated) {
      amrex::Abort("Device params not allocated yet");
    } else {
      amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, &m_h_parm, &m_h_parm + 1,
        m_d_parm);
    }
  }
  
  virtual void deallocate()
  {
    InitParm<ParmType>::host_deallocate(this);
    if (m_device_allocated)
      amrex::The_Device_Arena()->free(m_d_parm);
  }
  
  virtual ParmType& host_parm() { return m_h_parm; }

  virtual const ParmType* device_parm() { return m_d_parm; }

  private:
  HostOnlyParm<ParmType> m_host_only_parm;
  ParmType m_h_parm;
  ParmType* m_d_parm;
  bool m_device_allocated{false};
};

}
#endif
