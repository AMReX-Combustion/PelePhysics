#ifndef SPRAYINJECTTEMPLATE_H
#define SPRAYINJECTTEMPLATE_H
#include "SprayParticles.H"

/*
 This is generalized injection routine.
 We have 2 spherical coordinate systems
 1: Jet spherical coordinates for the center of the jet relative to
 Cartesian. This is the jet_angle (theta_1) and azimuthal angles (phi_1)
 2: Spread spherical coordinates for spread of particle relative
 center of jet. This uses the spread_angle
 Within these coordinates, theta is the inclination angle and phi is
 the azimuthal
*/

void
SprayParticleContainer::sprayInjection(
  const amrex::Real time,
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
  SprayUnits SPU;
  const SprayData* fdat = m_sprayData;
  amrex::Real rho_part = 0.;
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    rho_part +=
      spray_jet->get_avg_Y(spf) / fdat->rhoL(spray_jet->get_avg_T(), spf);
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
  amrex::Real mass_perc = inject_mass;
#if AMREX_SPACEDIM == 3
  amrex::Real injProcs = static_cast<amrex::Real>(numInjProcs);
  amrex::Real dphi = 2. * M_PI / injProcs;
  amrex::Real phi_rad_lo = amrex::ParallelDescriptor::MyProc() * dphi;
  mass_perc = inject_mass / injProcs;
#endif
  while (cur_mass < mass_perc) {
    // Pick random percentage from 0 to 1
    amrex::Real radp = amrex::Random();
#if AMREX_SPACEDIM == 3
    if (spray_jet->hollow_spray()) {
      radp = 1.;
    }
    amrex::Real phi_radial = phi_rad_lo + amrex::Random() * dphi;
    // This determines the radial location of the particle within the jet inlet
    amrex::Real cur_rad = radp * spray_jet->jet_dia() / 2.;
#else
    if (spray_jet->hollow_spray()) {
      if (radp <= 0.5) {
        radp = 0.;
      } else {
        radp = 1.;
      }
    }
    amrex::Real phi_radial = 0.;
    amrex::Real cur_rad = (radp - 0.5) * spray_jet->jet_dia();
#endif
    amrex::Real umag, theta_spread, phi_swirl, dia_part, T_part;
    amrex::Real Y_part[SPRAY_FUEL_NUM];
    bool make_new_part = spray_jet->get_new_particle(
      time, phi_radial, cur_rad, umag, theta_spread, phi_swirl, dia_part,
      T_part, Y_part);
    if (make_new_part && dia_part > min_dia) {
      amrex::RealVect part_loc, part_vel;
      spray_jet->transform_loc_vel(
        theta_spread, phi_radial, cur_rad, umag, phi_swirl, part_vel, part_loc);
      ParticleType p;
      p.id() = ParticleType::NextID();
      p.cpu() = amrex::ParallelDescriptor::MyProc();
      AMREX_D_TERM(p.rdata(pstateVel) = part_vel[0];
                   , p.rdata(pstateVel + 1) = part_vel[1];
                   , p.rdata(pstateVel + 2) = part_vel[2];);
      p.rdata(pstateT) = T_part;
      // Never add particle with less than minimum mass
      p.rdata(pstateDia) = dia_part;
      amrex::Real rho_part = 0.;
      if (SPRAY_FUEL_NUM > 1) {
        for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
          p.rdata(pstateY + spf) = Y_part[spf];
          rho_part += Y_part[spf] / fdat->rhoL(T_part, spf);
        }
        rho_part = 1. / rho_part;
      } else {
        rho_part = fdat->rhoL(T_part, 0);
        p.rdata(pstateY) = 1.;
      }
      // Add particles as if they have advanced some random portion of
      // dt
      amrex::Real pmov = amrex::Random();
      for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        p.pos(dir) = part_loc[dir] + pmov * dt * p.rdata(pstateVel + dir);
      }
      amrex::Real pmass = Pi_six * rho_part * std::pow(dia_part, 3);
      bool where = Where(p, pld);
      if (!where) {
        amrex::Abort("Bad injection particle");
      }
      std::pair<int, int> ind(pld.m_grid, pld.m_tile);
      host_particles[ind].push_back(p);
      cur_mass += num_ppp * pmass;
    }
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
