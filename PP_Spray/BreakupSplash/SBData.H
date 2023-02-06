#ifndef SBDATA_H
#define SBDATA_H

// This contains data SB (splashing or breakup) used for creating new droplets
// using data that is generated on device. Variables phi1 and phi2 will differ
// between if the droplet is splashing or breaking up.
struct SBPtrs
{
  amrex::Real* norm = nullptr;
  amrex::Real* vel = nullptr;
  amrex::Real* loc = nullptr;
  amrex::Real* dtpp = nullptr;
  amrex::Real* T0 = nullptr;
  amrex::Real* d0 = nullptr;
  amrex::Real* Y0 = nullptr;
  amrex::Real* numDens = nullptr;
  amrex::Real* phi1 = nullptr;
  amrex::Real* phi2 = nullptr;
};

struct SBVects
{
  // Normal vector of wall (for splashing)
  amrex::Gpu::HostVector<amrex::Real> norm_h;
  amrex::Gpu::DeviceVector<amrex::Real> norm_d;
  // Velocity of droplet
  amrex::Gpu::HostVector<amrex::Real> vel_h;
  amrex::Gpu::DeviceVector<amrex::Real> vel_d;
  // Location of droplet (placed at wall for splashing)
  amrex::Gpu::HostVector<amrex::Real> loc_h;
  amrex::Gpu::DeviceVector<amrex::Real> loc_d;
  // Dt fraction remaining after splashing/breakup
  amrex::Gpu::HostVector<amrex::Real> dtpp_h;
  amrex::Gpu::DeviceVector<amrex::Real> dtpp_d;
  // Droplet temperature
  amrex::Gpu::HostVector<amrex::Real> T0_h;
  amrex::Gpu::DeviceVector<amrex::Real> T0_d;
  // Droplet diameter
  // Splashing: Original droplet diameter
  // Breakup: Final droplet diameter after breakup
  amrex::Gpu::HostVector<amrex::Real> d0_h;
  amrex::Gpu::DeviceVector<amrex::Real> d0_d;
  // Droplet mass fractions
  amrex::Gpu::HostVector<amrex::Real> Y0_h;
  amrex::Gpu::DeviceVector<amrex::Real> Y0_d;
  // Variable
  // Splashing: Kv
  // Breakup: Utan, tangential velocity magnitude from breakup
  amrex::Gpu::HostVector<amrex::Real> phi1_h;
  amrex::Gpu::DeviceVector<amrex::Real> phi1_d;
  // Variable
  // Splashing: ms, amount of mass that splashes
  // TAB Breakup: TABY value
  // RT Breakup: distance from jet
  // KH Breakup: unused
  amrex::Gpu::HostVector<amrex::Real> phi2_h;
  amrex::Gpu::DeviceVector<amrex::Real> phi2_d;
  // Droplet number density
  amrex::Gpu::HostVector<amrex::Real> numDens_h;
  amrex::Gpu::DeviceVector<amrex::Real> numDens_d;
  SBVects(const int Np)
  {
    norm_h.assign(AMREX_SPACEDIM * Np, 0.);
    vel_h.assign(AMREX_SPACEDIM * Np, 0.);
    loc_h.assign(AMREX_SPACEDIM * Np, 0.);
    dtpp_h.assign(Np, 0.);
    T0_h.assign(Np, 0.);
    d0_h.assign(Np, 0.);
    numDens_h.assign(Np, 0.);
    phi1_h.assign(Np, 0.);
    phi2_h.assign(Np, 0.);

    norm_d.resize(AMREX_SPACEDIM * Np);
    vel_d.resize(AMREX_SPACEDIM * Np);
    loc_d.resize(AMREX_SPACEDIM * Np);
    dtpp_d.resize(Np);
    T0_d.resize(Np);
    numDens_d.resize(Np);
    d0_d.resize(Np);
    phi1_d.resize(Np);
    phi2_d.resize(Np);
#if SPRAY_FUEL_NUM > 1
    Y0_h.assign(SPRAY_FUEL_NUM * Np, 0.);
    Y0_d.resize(SPRAY_FUEL_NUM * Np);
    amrex::Gpu::copyAsync(amrex::Gpu::hostToDevice, Y0_h.begin(), Y0_h.end(), Y0_d.begin());
#endif
    amrex::Gpu::copyAsync(
      amrex::Gpu::hostToDevice, norm_h.begin(), norm_h.end(), norm_d.begin());
    amrex::Gpu::copyAsync(
      amrex::Gpu::hostToDevice, vel_h.begin(), vel_h.end(), vel_d.begin());
    amrex::Gpu::copyAsync(
      amrex::Gpu::hostToDevice, loc_h.begin(), loc_h.end(), loc_d.begin());
    amrex::Gpu::copyAsync(
      amrex::Gpu::hostToDevice, dtpp_h.begin(), dtpp_h.end(), dtpp_d.begin());
    amrex::Gpu::copyAsync(
      amrex::Gpu::hostToDevice, T0_h.begin(), T0_h.end(), T0_d.begin());
    amrex::Gpu::copyAsync(
      amrex::Gpu::hostToDevice, numDens_h.begin(), numDens_h.end(),
      numDens_d.begin());
    amrex::Gpu::copyAsync(
      amrex::Gpu::hostToDevice, d0_h.begin(), d0_h.end(), d0_d.begin());
    amrex::Gpu::copyAsync(
      amrex::Gpu::hostToDevice, phi1_h.begin(), phi1_h.end(), phi1_d.begin());
    amrex::Gpu::copyAsync(
      amrex::Gpu::hostToDevice, phi2_h.begin(), phi2_h.end(), phi2_d.begin());
  }
  SBVects() = delete;
  SBVects(const SBVects&) = delete;

  void retrieve_data()
  {
#if SPRAY_FUEL_NUM > 1
    amrex::Gpu::copy(
      amrex::Gpu::deviceToHost, Y0_d.begin(), Y0_d.end(), Y0_h.begin());
#endif
    amrex::Gpu::copy(
      amrex::Gpu::deviceToHost, norm_d.begin(), norm_d.end(), norm_h.begin());
    amrex::Gpu::copy(
      amrex::Gpu::deviceToHost, vel_d.begin(), vel_d.end(), vel_h.begin());
    amrex::Gpu::copy(
      amrex::Gpu::deviceToHost, loc_d.begin(), loc_d.end(), loc_h.begin());
    amrex::Gpu::copy(
      amrex::Gpu::deviceToHost, dtpp_d.begin(), dtpp_d.end(), dtpp_h.begin());
    amrex::Gpu::copy(
      amrex::Gpu::deviceToHost, T0_d.begin(), T0_d.end(), T0_h.begin());
    amrex::Gpu::copy(
      amrex::Gpu::deviceToHost, numDens_d.begin(), numDens_d.end(),
      numDens_h.begin());
    amrex::Gpu::copy(
      amrex::Gpu::deviceToHost, d0_d.begin(), d0_d.end(), d0_h.begin());
    amrex::Gpu::copy(
      amrex::Gpu::deviceToHost, phi1_d.begin(), phi1_d.end(), phi1_h.begin());
    amrex::Gpu::copy(
      amrex::Gpu::deviceToHost, phi2_d.begin(), phi2_d.end(), phi2_h.begin());
  }

  void fillPtrs_d(SBPtrs& rf)
  {
    rf.norm = norm_d.data();
    rf.vel = vel_d.data();
    rf.loc = loc_d.data();
    rf.dtpp = dtpp_d.data();
    rf.T0 = T0_d.data();
    rf.numDens = numDens_d.data();
    rf.d0 = d0_d.data();
    rf.Y0 = Y0_d.data();
    rf.phi1 = phi1_d.data();
    rf.phi2 = phi2_d.data();
  }
  void fillPtrs_h(SBPtrs& rf)
  {
    rf.norm = norm_h.data();
    rf.vel = vel_h.data();
    rf.loc = loc_h.data();
    rf.dtpp = dtpp_h.data();
    rf.T0 = T0_h.data();
    rf.numDens = numDens_h.data();
    rf.d0 = d0_h.data();
    rf.Y0 = Y0_h.data();
    rf.phi1 = phi1_h.data();
    rf.phi2 = phi2_h.data();
  }
};

#endif
