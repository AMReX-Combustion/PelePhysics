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
  amrex::Real* T0 = nullptr;
  amrex::Real* ref_dia = nullptr;
  amrex::Real* Y0 = nullptr;
  amrex::Real* num_dens = nullptr;
  amrex::Real* phi1 = nullptr;
  amrex::Real* phi2 = nullptr;
  amrex::Real* phi3 = nullptr;
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
  // Droplet temperature
  amrex::Gpu::HostVector<amrex::Real> T0_h;
  amrex::Gpu::DeviceVector<amrex::Real> T0_d;
  // Droplet diameter
  // Splashing: Original droplet diameter
  // Breakup: Final droplet diameter after breakup
  amrex::Gpu::HostVector<amrex::Real> ref_dia_h;
  amrex::Gpu::DeviceVector<amrex::Real> ref_dia_d;
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
  // KH-RT Breakup: parent droplet original number density
  amrex::Gpu::HostVector<amrex::Real> phi2_h;
  amrex::Gpu::DeviceVector<amrex::Real> phi2_d;
  // Variable
  // Splashing: film thickness / drop diameter
  // TAB Breakup: TABY_dot value
  // KH-RT Breakup: Unused
  amrex::Gpu::HostVector<amrex::Real> phi3_h;
  amrex::Gpu::DeviceVector<amrex::Real> phi3_d;
  // Variable
  // Splashing: Original parcel number density
  // Breakup: Total number of created droplets
  amrex::Gpu::HostVector<amrex::Real> num_dens_h;
  amrex::Gpu::DeviceVector<amrex::Real> num_dens_d;

  SBVects() = default;

  void build(const int Np)
  {
    norm_h.assign(AMREX_SPACEDIM * Np, 0.);
    vel_h.assign(AMREX_SPACEDIM * Np, 0.);
    loc_h.assign(AMREX_SPACEDIM * Np, 0.);
    T0_h.assign(Np, 0.);
    ref_dia_h.assign(Np, 0.);
    num_dens_h.assign(Np, 0.);
    phi1_h.assign(Np, 0.);
    phi2_h.assign(Np, 0.);
    phi3_h.assign(Np, 0.);

    norm_d.resize(AMREX_SPACEDIM * Np);
    vel_d.resize(AMREX_SPACEDIM * Np);
    loc_d.resize(AMREX_SPACEDIM * Np);
    T0_d.resize(Np);
    num_dens_d.resize(Np);
    ref_dia_d.resize(Np);
    phi1_d.resize(Np);
    phi2_d.resize(Np);
    phi3_d.resize(Np);
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
      amrex::Gpu::hostToDevice, T0_h.begin(), T0_h.end(), T0_d.begin());
    amrex::Gpu::copyAsync(
      amrex::Gpu::hostToDevice, num_dens_h.begin(), num_dens_h.end(),
      num_dens_d.begin());
    amrex::Gpu::copyAsync(
      amrex::Gpu::hostToDevice, ref_dia_h.begin(), ref_dia_h.end(),
      ref_dia_d.begin());
    amrex::Gpu::copyAsync(
      amrex::Gpu::hostToDevice, phi1_h.begin(), phi1_h.end(), phi1_d.begin());
    amrex::Gpu::copyAsync(
      amrex::Gpu::hostToDevice, phi2_h.begin(), phi2_h.end(), phi2_d.begin());
    amrex::Gpu::copyAsync(
      amrex::Gpu::hostToDevice, phi3_h.begin(), phi3_h.end(), phi3_d.begin());
  }

  SBVects(const SBVects&) = delete;

  void retrieve_data()
  {
#if SPRAY_FUEL_NUM > 1
    amrex::Gpu::copyAsync(
      amrex::Gpu::deviceToHost, Y0_d.begin(), Y0_d.end(), Y0_h.begin());
#endif
    amrex::Gpu::copyAsync(
      amrex::Gpu::deviceToHost, norm_d.begin(), norm_d.end(), norm_h.begin());
    amrex::Gpu::copyAsync(
      amrex::Gpu::deviceToHost, vel_d.begin(), vel_d.end(), vel_h.begin());
    amrex::Gpu::copyAsync(
      amrex::Gpu::deviceToHost, loc_d.begin(), loc_d.end(), loc_h.begin());
    amrex::Gpu::copyAsync(
      amrex::Gpu::deviceToHost, T0_d.begin(), T0_d.end(), T0_h.begin());
    amrex::Gpu::copyAsync(
      amrex::Gpu::deviceToHost, num_dens_d.begin(), num_dens_d.end(),
      num_dens_h.begin());
    amrex::Gpu::copyAsync(
      amrex::Gpu::deviceToHost, ref_dia_d.begin(), ref_dia_d.end(),
      ref_dia_h.begin());
    amrex::Gpu::copyAsync(
      amrex::Gpu::deviceToHost, phi1_d.begin(), phi1_d.end(), phi1_h.begin());
    amrex::Gpu::copyAsync(
      amrex::Gpu::deviceToHost, phi2_d.begin(), phi2_d.end(), phi2_h.begin());
    amrex::Gpu::copyAsync(
      amrex::Gpu::deviceToHost, phi3_d.begin(), phi3_d.end(), phi3_h.begin());
    amrex::Gpu::Device::streamSynchronize();
  }

  void fillPtrs_d(SBPtrs& rf)
  {
    rf.norm = norm_d.data();
    rf.vel = vel_d.data();
    rf.loc = loc_d.data();
    rf.T0 = T0_d.data();
    rf.num_dens = num_dens_d.data();
    rf.ref_dia = ref_dia_d.data();
    rf.Y0 = Y0_d.data();
    rf.phi1 = phi1_d.data();
    rf.phi2 = phi2_d.data();
    rf.phi3 = phi3_d.data();
  }
  void fillPtrs_h(SBPtrs& rf)
  {
    rf.norm = norm_h.data();
    rf.vel = vel_h.data();
    rf.loc = loc_h.data();
    rf.T0 = T0_h.data();
    rf.num_dens = num_dens_h.data();
    rf.ref_dia = ref_dia_h.data();
    rf.Y0 = Y0_h.data();
    rf.phi1 = phi1_h.data();
    rf.phi2 = phi2_h.data();
    rf.phi3 = phi3_h.data();
  }
};

#endif
