#ifndef SPRAYFUELDATA_H_
#define SPRAYFUELDATA_H_

#include <AMReX_Gpu.H>
#include "mechanism.h"

using namespace amrex;

struct SprayData
{
  SprayData (const Gpu::DeviceVector<Real>& critT,
             const Gpu::DeviceVector<Real>& boilT,
             const Gpu::DeviceVector<Real>& fuelCp,
             const Gpu::DeviceVector<Real>& fuelLatent,
             const Gpu::DeviceVector<int>&  fuelIndx)
    : m_critT_ptr(critT.dataPtr()),
      m_boilT_ptr(boilT.dataPtr()),
      m_fuelCp_ptr(fuelCp.dataPtr()),
      m_fuelLatent_ptr(fuelLatent.dataPtr()),
      m_fuelIndx_ptr(fuelIndx.dataPtr())
  {}

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  Real critT (const int fuelSpecIndx) const
  {
    return m_critT_ptr[fuelSpecIndx];
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  Real boilT (const int fuelSpecIndx) const
  {
    return m_boilT_ptr[fuelSpecIndx];
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  Real cp (const int fuelSpecIndx) const
  {
    return m_fuelCp_ptr[fuelSpecIndx];
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  Real latent (const int fuelSpecIndx) const
  {
    return m_fuelLatent_ptr[fuelSpecIndx];
  }

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  int indx (const int fuelSpecIndx) const
  {
    return m_fuelIndx_ptr[fuelSpecIndx];
  }

private:

  const Real * m_critT_ptr;
  const Real * m_boilT_ptr;
  const Real * m_fuelCp_ptr;
  const Real * m_fuelLatent_ptr;
  const int * m_fuelIndx_ptr;
};

class SprayDataContainer
{
public:

  void build (const Gpu::HostVector<Real>& critT,
              const Gpu::HostVector<Real>& boilT,
              const Gpu::HostVector<Real>& fuelCp,
              const Gpu::HostVector<Real>& refFuelH,
              const Gpu::HostVector<Real>& fuelLatent,
              const Gpu::HostVector<int>&  fuelIndx)
  {
    const int nfspec = SPRAY_FUEL_NUM;
    m_critT.resize(nfspec);
    m_boilT.resize(nfspec);
    m_fuelCp.resize(nfspec);
    m_fuelLatent.resize(nfspec);
    m_fuelIndx.resize(nfspec);
    Gpu::HostVector<Real> tmpFuelLatent(nfspec);
    for (int i = 0; i != nfspec; ++i) {
      const int specIndx = fuelIndx[i];
      tmpFuelLatent[i] = fuelLatent[i];
#ifndef LEGACY_SPRAY
      tmpFuelLatent[i] -= refFuelH[specIndx];
#endif
    }
    Gpu::copy(Gpu::hostToDevice, critT.begin(), critT.end(),
              m_critT.begin());
    Gpu::copy(Gpu::hostToDevice, boilT.begin(), boilT.end(),
              m_boilT.begin());
    Gpu::copy(Gpu::hostToDevice, fuelCp.begin(), fuelCp.end(),
              m_fuelCp.begin());
    Gpu::copy(Gpu::hostToDevice, tmpFuelLatent.begin(), tmpFuelLatent.end(),
              m_fuelLatent.begin());
    Gpu::copy(Gpu::hostToDevice, fuelIndx.begin(), fuelIndx.end(),
              m_fuelIndx.begin());
  }

  SprayData getSprayData()
  {
    return SprayData(m_critT, m_boilT, m_fuelCp, m_fuelLatent, m_fuelIndx);
  }

private:

  // Vectors with constant reference values
  Gpu::DeviceVector<Real> m_critT;
  Gpu::DeviceVector<Real> m_boilT;
  Gpu::DeviceVector<Real> m_fuelCp;
  Gpu::DeviceVector<Real> m_fuelLatent;
  Gpu::DeviceVector<int>  m_fuelIndx;
};

#endif
