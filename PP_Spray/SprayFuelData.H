#ifndef SPRAYFUELDATA_H_
#define SPRAYFUELDATA_H_

#include <AMReX_Gpu.H>
#include "mechanism.h"

using namespace amrex;

class SprayFuelData
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
    Gpu::HostVector<Real> invBoilT(nfspec);
    Gpu::HostVector<Real> tmpFuelLatent(nfspec);
    m_invBoilT.resize(nfspec);
    for (int i = 0; i != nfspec; ++i) {
      const int specIndx = fuelIndx[i];
      tmpFuelLatent[i] = fuelLatent[i];
#ifndef LEGACY_SPRAY
      tmpFuelLatent[i] -= refFuelH[specIndx];
#endif
      invBoilT[i] = 1./boilT[i];
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
    Gpu::copy(Gpu::hostToDevice, invBoilT.begin(), invBoilT.end(),
	      m_invBoilT.begin());
  }

  //! \brief returns the const pointer to the critical temperature for the liquid fuels
  const Real* critT() const noexcept { return m_critT.dataPtr(); }

  //! \brief returns the const pointer to the boiling temeperature for the liquid fuels
  const Real* boilT() const noexcept { return m_boilT.dataPtr(); }

  //! \brief returns the const pointer to the specific heat for the liquid fuels
  const Real* fuelCp() const noexcept { return m_fuelCp.dataPtr(); }

  //! \brief returns the const pointer to the latent heat for the spray fuel
  const Real* fuelLatent() const noexcept { return m_fuelLatent.dataPtr(); }

  //! \brief returns the const pointer to the index of the fuel species relative to the gas species
  const int* fuelIndx() const noexcept { return m_fuelIndx.dataPtr(); }

  //! \brief returns the const pointer to the inverse of the boiling temperature
  const Real* invBoilT() const noexcept { return m_invBoilT.dataPtr(); }

private:

  // Vectors with constant reference values
  Gpu::DeviceVector<Real> m_critT;
  Gpu::DeviceVector<Real> m_boilT;
  Gpu::DeviceVector<Real> m_fuelCp;
  Gpu::DeviceVector<Real> m_fuelLatent;
  Gpu::DeviceVector<int>  m_fuelIndx;
  Gpu::DeviceVector<Real> m_invBoilT;
};

#endif
