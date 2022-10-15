#ifndef SPRAYREFLDATA_H_
#define SPRAYREFLDATA_H_

struct ReflPtrs
{
  amrex::Real* norm = nullptr;
  amrex::Real* vel = nullptr;
  amrex::Real* loc = nullptr;
  amrex::Real* alpha = nullptr;
  amrex::Real* Kv = nullptr;
  amrex::Real* dtpp = nullptr;
  amrex::Real* ms = nullptr;
  amrex::Real* T0 = nullptr;
  amrex::Real* d0 = nullptr;
  amrex::Real* Y0 = nullptr;
};

struct ReflVects
{
  // Normal vector
  amrex::Gpu::HostVector<amrex::Real> norm_h;
  amrex::Gpu::DeviceVector<amrex::Real> norm_d;
  // Velocity of splashing droplet
  amrex::Gpu::HostVector<amrex::Real> vel_h;
  amrex::Gpu::DeviceVector<amrex::Real> vel_d;
  // Location of splashing droplet (at wall)
  amrex::Gpu::HostVector<amrex::Real> loc_h;
  amrex::Gpu::DeviceVector<amrex::Real> loc_d;
  // Alpha value from Kuhnke model
  amrex::Gpu::HostVector<amrex::Real> alpha_h;
  amrex::Gpu::DeviceVector<amrex::Real> alpha_d;
  // Weber number
  amrex::Gpu::HostVector<amrex::Real> kv_h;
  amrex::Gpu::DeviceVector<amrex::Real> kv_d;
  // Dt fraction remaining after splashing
  amrex::Gpu::HostVector<amrex::Real> dtpp_h;
  amrex::Gpu::DeviceVector<amrex::Real> dtpp_d;
  // T* from Kuhnke model
  amrex::Gpu::HostVector<amrex::Real> ms_h;
  amrex::Gpu::DeviceVector<amrex::Real> ms_d;
  // Splashing droplet temperature
  amrex::Gpu::HostVector<amrex::Real> T0_h;
  amrex::Gpu::DeviceVector<amrex::Real> T0_d;
  // Splashing droplet diameter
  amrex::Gpu::HostVector<amrex::Real> d0_h;
  amrex::Gpu::DeviceVector<amrex::Real> d0_d;
  // Splashing droplet mass fractions
  amrex::Gpu::HostVector<amrex::Real> Y0_h;
  amrex::Gpu::DeviceVector<amrex::Real> Y0_d;
  ReflVects(const int Np)
  {
    norm_h.assign(AMREX_SPACEDIM * Np, 0.);
    vel_h.assign(AMREX_SPACEDIM * Np, 0.);
    loc_h.assign(AMREX_SPACEDIM * Np, 0.);
    alpha_h.assign(Np, 0.);
    kv_h.assign(Np, 0.);
    dtpp_h.assign(Np, 0.);
    ms_h.assign(Np, 0.);
    T0_h.assign(Np, 0.);
    d0_h.assign(Np, 0.);

    norm_d.resize(AMREX_SPACEDIM * Np);
    vel_d.resize(AMREX_SPACEDIM * Np);
    loc_d.resize(AMREX_SPACEDIM * Np);
    alpha_d.resize(Np);
    kv_d.resize(Np);
    dtpp_d.resize(Np);
    ms_d.resize(Np);
    T0_d.resize(Np);
    d0_d.resize(Np);
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
      amrex::Gpu::hostToDevice, alpha_h.begin(), alpha_h.end(),
      alpha_d.begin());
    amrex::Gpu::copyAsync(
      amrex::Gpu::hostToDevice, kv_h.begin(), kv_h.end(), kv_d.begin());
    amrex::Gpu::copyAsync(
      amrex::Gpu::hostToDevice, dtpp_h.begin(), dtpp_h.end(), dtpp_d.begin());
    amrex::Gpu::copyAsync(
      amrex::Gpu::hostToDevice, ms_h.begin(), ms_h.end(), ms_d.begin());
    amrex::Gpu::copyAsync(
      amrex::Gpu::hostToDevice, T0_h.begin(), T0_h.end(), T0_d.begin());
    amrex::Gpu::copyAsync(
      amrex::Gpu::hostToDevice, d0_h.begin(), d0_h.end(), d0_d.begin());
  }
  ReflVects() = delete;
  ReflVects (const ReflVects&) = delete;

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
      amrex::Gpu::deviceToHost, alpha_d.begin(), alpha_d.end(),
      alpha_h.begin());
    amrex::Gpu::copy(
      amrex::Gpu::deviceToHost, kv_d.begin(), kv_d.end(), kv_h.begin());
    amrex::Gpu::copy(
      amrex::Gpu::deviceToHost, dtpp_d.begin(), dtpp_d.end(), dtpp_h.begin());
    amrex::Gpu::copy(
      amrex::Gpu::deviceToHost, ms_d.begin(), ms_d.end(), ms_h.begin());
    amrex::Gpu::copy(
      amrex::Gpu::deviceToHost, T0_d.begin(), T0_d.end(), T0_h.begin());
    amrex::Gpu::copy(
      amrex::Gpu::deviceToHost, d0_d.begin(), d0_d.end(), d0_h.begin());
  }

  void fillPtrs_d(ReflPtrs& rf)
  {
    rf.norm = norm_d.data();
    rf.vel = vel_d.data();
    rf.loc = loc_d.data();
    rf.alpha = alpha_d.data();
    rf.Kv = kv_d.data();
    rf.dtpp = dtpp_d.data();
    rf.ms = ms_d.data();
    rf.T0 = T0_d.data();
    rf.d0 = d0_d.data();
    rf.Y0 = Y0_d.data();
  }
  void fillPtrs_h(ReflPtrs& rf)
  {
    rf.norm = norm_h.data();
    rf.vel = vel_h.data();
    rf.loc = loc_h.data();
    rf.alpha = alpha_h.data();
    rf.Kv = kv_h.data();
    rf.dtpp = dtpp_h.data();
    rf.ms = ms_h.data();
    rf.T0 = T0_h.data();
    rf.d0 = d0_h.data();
    rf.Y0 = Y0_h.data();
  }
};

#endif
