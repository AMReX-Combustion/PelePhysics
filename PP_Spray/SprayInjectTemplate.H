#ifndef SPRAYINJECTTEMPLATE_H
#define SPRAYINJECTTEMPLATE_H
#include "SprayParticles.H"

// This provides a function to be called in the SprayParticlesInitInsert.cpp
// file included is a standard jet injection calls
#if AMREX_SPACEDIM == 3
amrex::RealVect
jet_vel_comp(
  amrex::Real theta_1,
  amrex::Real phi_1,
  amrex::Real theta_2,
  amrex::Real phi_2)
{
  amrex::Real st1 = std::sin(theta_1);
  amrex::Real ct1 = std::cos(theta_1);
  amrex::Real st2 = std::sin(theta_2);
  amrex::Real ct2 = std::cos(theta_2);
  amrex::Real sp1 = std::sin(phi_1);
  amrex::Real cp1 = std::cos(phi_1);
  amrex::Real sp2 = std::sin(phi_2);
  amrex::Real cp2 = std::cos(phi_2);
  amrex::Real v1 = st1 * ct2 + st2 * cp2 * ct1;
  amrex::RealVect vel(
    cp1 * v1 - sp1 * sp2 * st2, sp1 * v1 + sp2 * st2 * cp1,
    ct1 * ct2 - st1 * st2 * cp2);
  return vel;
}
#endif

amrex::RealVect
part_loc_trans(
  amrex::Real radius, amrex::Real theta_1, amrex::Real phi_1, amrex::Real phi_2)
{
  amrex::Real ct1 = std::cos(theta_1);
  amrex::Real sp1 = std::sin(phi_1);
  amrex::Real cp1 = std::cos(phi_1);
  amrex::Real sp2 = std::sin(phi_2);
  amrex::Real cp2 = std::cos(phi_2);
  amrex::RealVect loc(AMREX_D_DECL(
    radius * (cp1 * cp2 * ct1 - sp1 * sp2),
    radius * (sp1 * cp2 * ct1 + cp1 * sp2), -radius * std::sin(theta_1) * cp2));
  return loc;
}

/*
 This is generalized injection routine.
 We have 2 spherical coordinate systems
 1: Jet spherical coordinates for the center of the jet relative to
 Cartesian. This is the jet_angle and azimuthal angles
 2: Spread spherical coordinates for spread of particle relative
 center of jet. This uses the spread_angle
 Within these coordinates, theta is the inclination angle and phi is
 the azimuthal
*/

void
SprayParticleContainer::sprayInjection(
  SprayJet* spray_jet,
  const amrex::Real sim_dt,
  const int level,
  int startInjProc,
  int numInjProcs,
  const int min_parcel)
{
  // Don't try to use more procs than we have
  numInjProcs = amrex::min(amrex::ParallelDescriptor::NProcs(), numInjProcs);
  startInjProc =
    amrex::min(amrex::ParallelDescriptor::NProcs() - numInjProcs, startInjProc);
#if AMREX_SPACEDIM == 2
  // Do not parallelize injection for 2D
  numInjProcs = 1;
#endif
  int curProc = amrex::ParallelDescriptor::MyProc();
  if (curProc < startInjProc || curProc >= startInjProc + numInjProcs) {
    return;
  }
  const int pstateVel = m_sprayIndx.pstateVel;
  const int pstateT = m_sprayIndx.pstateT;
  const int pstateDia = m_sprayIndx.pstateDia;
  const int pstateY = m_sprayIndx.pstateY;
  const SprayUnits SPU;
  const SprayData* fdat = m_sprayData;
  amrex::Real rho_part = 0.;
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    rho_part += spray_jet->Y(spf) / fdat->rhoL(spray_jet->T(), spf);
  }
  rho_part = 1. / rho_part;
  const amrex::Real Pi_six = M_PI / 6.;
  // Number of particles per parcel
  const amrex::Real num_ppp = fdat->num_ppp;
  // Check if mass must be injected across multiple timesteps
  const amrex::Real avg_dia = spray_jet->get_avg_dia();
  const amrex::Real avg_mass = Pi_six * rho_part * std::pow(avg_dia, 3);
  const amrex::Real min_dia = std::cbrt(6. * SPU.min_mass / (M_PI * rho_part));
  if (avg_dia < min_dia || avg_mass < SPU.min_mass) {
    amrex::Abort(
      "Average droplet size too small, floating point issues expected");
  }
  // Minimum number of parcels we want injected at a time
  const amrex::Real min_inj_parcel = static_cast<amrex::Real>(min_parcel);
  amrex::Real dt = sim_dt;
  amrex::Real inject_mass = spray_jet->mass_flow_rate() * sim_dt;
  // See if there is uninjected mass from previous time steps
  if (spray_jet->m_sumInjMass > 0.) {
    dt += spray_jet->m_sumInjTime;
    inject_mass += spray_jet->m_sumInjMass;
  }
  if (inject_mass / (num_ppp * avg_mass) < min_inj_parcel) {
    spray_jet->m_sumInjMass = inject_mass;
    spray_jet->m_sumInjTime = dt;
    return;
  }

  amrex::ParticleLocData pld;
  std::map<std::pair<int, int>, amrex::Gpu::HostVector<ParticleType>>
    host_particles;
  amrex::Real cur_mass = 0.;
  amrex::RealVect jet_norm = spray_jet->jet_norm();
  amrex::RealVect jet_cent = spray_jet->jet_cent();
  amrex::Real jet_vel = spray_jet->jet_vel();
  amrex::Real injProcs = static_cast<amrex::Real>(numInjProcs);
#if AMREX_SPACEDIM == 3
  amrex::Real norm = jet_norm.vectorLength();
  amrex::Real theta_1 = std::acos(jet_norm[2] / norm);
  amrex::Real phi_1 = std::atan2(jet_norm[1] / norm, jet_norm[0] / norm);
  amrex::Real dphi = 2. * M_PI / injProcs;
  amrex::Real phi_2_lo = amrex::ParallelDescriptor::MyProc() * dphi;
#else
  amrex::Real theta_1 = std::atan2(jet_norm[1], jet_norm[0]) + M_PI / 2.;
#endif
  amrex::Real mass_perc = inject_mass / injProcs;
  while (cur_mass < mass_perc) {
    // Pick random percentage from 0 to 1
    amrex::Real radp = amrex::Random();
#if AMREX_SPACEDIM == 3
    if (spray_jet->hollow_spray()) {
      radp = 1.;
    }
    // This determines the radial location of the particle within the jet inlet
    amrex::Real cur_rad = radp * spray_jet->jet_dia() / 2.;
    // The closer the particle is to the edge of the injection circle, the
    // greater the spread angle
    amrex::Real theta_2 = radp * spray_jet->spread_angle() / 2.;
    // Pick random azimuthal angle relative to jet inlet CS
    // Each processor is only injecting a portion of the circle
    amrex::Real phi_2 = phi_2_lo + amrex::Random() * dphi;
    amrex::RealVect part_loc =
      jet_cent + part_loc_trans(cur_rad, theta_1, phi_1, phi_2);
    amrex::RealVect part_vel = jet_vel_comp(theta_1, phi_1, theta_2, phi_2);
#else
    if (spray_jet->hollow_spray()) {
      if (radp <= 0.5) {
        radp = 0.;
      } else {
        radp = 1.;
      }
    }
    amrex::Real cur_rad = (radp - 0.5) * spray_jet->jet_dia();
    amrex::Real theta_2 = -(radp - 0.5) * spray_jet->spread_angle();
    amrex::Real st1 = std::sin(theta_1);
    amrex::Real ct1 = std::cos(theta_1);
    amrex::Real st2 = std::sin(theta_2);
    amrex::Real ct2 = std::cos(theta_2);
    amrex::RealVect part_loc(
      jet_cent[0] + cur_rad * ct1, jet_cent[1] + cur_rad * st1);
    amrex::RealVect part_vel(st1 * ct2 - st2 * ct1, -ct1 * ct2 - st1 * st2);
#endif
    ParticleType p;
    p.id() = ParticleType::NextID();
    p.cpu() = amrex::ParallelDescriptor::MyProc();
    AMREX_D_TERM(p.rdata(pstateVel) = jet_vel * part_vel[0];
                 , p.rdata(pstateVel + 1) = jet_vel * part_vel[1];
                 , p.rdata(pstateVel + 2) = jet_vel * part_vel[2];);
    // Never add particle with less than minimum mass
    amrex::Real cur_dia = amrex::max(min_dia, spray_jet->get_dia());
    // Add particles as if they have advanced some random portion of
    // dt
    amrex::Real pmov = amrex::Random();
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      p.pos(dir) = part_loc[dir] + pmov * dt * p.rdata(pstateVel + dir);
    }
    p.rdata(pstateT) = spray_jet->T();
    p.rdata(pstateDia) = cur_dia;
    for (int sp = 0; sp < SPRAY_FUEL_NUM; ++sp) {
      p.rdata(pstateY + sp) = spray_jet->Y(sp);
    }
    amrex::Real pmass = Pi_six * rho_part * std::pow(cur_dia, 3);
    bool where = Where(p, pld);
    if (!where) {
      amrex::Abort("Bad injection particle");
    }
    std::pair<int, int> ind(pld.m_grid, pld.m_tile);
    host_particles[ind].push_back(p);
    cur_mass += num_ppp * pmass;
  }
  for (auto& kv : host_particles) {
    auto grid = kv.first.first;
    auto tile = kv.first.second;
    const auto& src_tile = kv.second;
    auto& dst_tile = GetParticles(level)[std::make_pair(grid, tile)];
    auto old_size = dst_tile.GetArrayOfStructs().size();
    auto new_size = old_size + src_tile.size();
    dst_tile.resize(new_size);
    // Copy the AoS part of the host particles to the GPU
    amrex::Gpu::copy(
      amrex::Gpu::hostToDevice, src_tile.begin(), src_tile.end(),
      dst_tile.GetArrayOfStructs().begin() + old_size);
  }
  spray_jet->reset_sum();
}
#endif
