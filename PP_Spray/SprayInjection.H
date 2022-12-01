#ifndef SPRAYINJECTION_H
#define SPRAYINJECTION_H
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
  const int pstateTABY = m_sprayIndx.pstateTABY;
  const int pstateTABYdot = m_sprayIndx.pstateTABYdot;
  SprayUnits SPU;
  const SprayData* fdat = m_sprayData;
  amrex::Real rho_avg = 0.;
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    rho_avg +=
      spray_jet->get_avg_Y(spf) / fdat->rhoL(spray_jet->get_avg_T(), spf);
  }
  rho_avg = 1. / rho_avg;
  const amrex::Real rho_gas = spray_jet->get_rho_gas();
  const amrex::Real jet_C = spray_jet->get_jet_C();
  const amrex::Real Wec = 12.;
  if (fdat->use_ETAB && rho_gas < 0.) {
    amrex::Abort("Must call 'set_rho_gas()' if ETAB model is to be used");
  }
  const amrex::Real Pi_six = M_PI / 6.;
  // Number of particles per parcel
  const amrex::Real num_ppp = fdat->num_ppp;
  // Check if mass must be injected across multiple timesteps
  const amrex::Real avg_dia = spray_jet->get_avg_dia();
  const amrex::Real avg_mass = Pi_six * rho_avg * std::pow(avg_dia, 3);
  const amrex::Real min_dia = std::cbrt(6. * SPU.min_mass / (M_PI * rho_avg));
  if (avg_dia < min_dia || avg_mass < SPU.min_mass) {
    amrex::Abort(
      "Average droplet size too small, floating point issues expected");
  }
  // Minimum number of parcels we want injected at a time
  const auto min_inj_parcel = static_cast<amrex::Real>(min_parcel);
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
#if AMREX_SPACEDIM == 3
  auto injProcs = static_cast<amrex::Real>(numInjProcs);
  amrex::Real dphi = 2. * M_PI / injProcs;
  amrex::Real phi_rad_lo = amrex::ParallelDescriptor::MyProc() * dphi;
  amrex::Real mass_perc = inject_mass / injProcs;
#else
  amrex::Real mass_perc = inject_mass;
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
      p.rdata(pstateTABY) = 0.;
      // If ETAB model is applied
      if (fdat->use_ETAB) {
        amrex::Real velMag = vel_part.radSquared();
        amrex::Real omega =
          std::sqrt(8. * fdat->sigma / (rho_part * std::pow(dia_part / 2., 3)));
        amrex::Real We = rho_gas * dia_part / 2. * velMag / fdat->sigma;
        amrex::Real otbu = 4. * jet_C * std::sqrt(2. / We);
        p.rdata(pstateTABYdot) =
          (1. - We / Wec * (1. - std::cos(otbu))) * omega / std::sin(otbu);
      } else {
        p.rdata(pstateTABYdot) = 0.;
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

amrex::IntVect
unflatten_particles(const amrex::ULong idx, const amrex::IntVect& max_parts)
{
  amrex::IntVect indx;
  amrex::ULong cidx = idx;
  amrex::ULong d1 = max_parts[0];
#if AMREX_SPACEDIM > 2
  amrex::ULong d2 = max_parts[1];
  indx[2] = int(cidx / (d1 * d2));
  cidx -= amrex::ULong(indx[2]) * d1 * d2;
#endif
  indx[1] = int(cidx / d1);
  indx[0] = int(cidx % d1);
  return indx;
}

void
SprayParticleContainer::uniformSprayInit(
  const amrex::IntVect num_part,
  const amrex::RealVect vel_part,
  const amrex::Real dia_part,
  const amrex::Real T_part,
  const amrex::Real* Y_part,
  const int level,
  const int numRedist)
{
  const int MyProc = amrex::ParallelDescriptor::MyProc();
  const int NProcs = amrex::ParallelDescriptor::NProcs();
  int NRedist = numRedist;
  // TODO: This might be overkill but issues persisted at high Summit node
  // counts
  if (NRedist < 0) {
    NRedist = 1;
    if (NProcs <= 1024) {
      NRedist = 2;
    } else if (NProcs <= 2048) {
      NRedist = 4;
    } else if (NProcs <= 4096) {
      NRedist = 8;
    }
  }
  const int pstateVel = m_sprayIndx.pstateVel;
  const int pstateDia = m_sprayIndx.pstateDia;
  const int pstateT = m_sprayIndx.pstateT;
  const int pstateY = m_sprayIndx.pstateY;
  const int pstatePb = m_sprayIndx.pstatePb;
  const int pstatePbdot = m_sprayIndx.pstatePbdot;
  // Reference values for the particles
  amrex::Real part_vals[NAR_SPR + NSR_SPR];
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    part_vals[pstateVel + dir] = vel_part[dir];
  }
  part_vals[pstateT] = T_part;
  part_vals[pstateDia] = dia_part;
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    part_vals[pstateY + spf] = Y_part[spf];
  }
  part_vals[pstatePb] = 0.;
  part_vals[pstatePbdot] = 0.;
  const amrex::RealVect dx_part(AMREX_D_DECL(
    Geom(level).ProbLength(0) / amrex::Real(num_part[0]),
    Geom(level).ProbLength(1) / amrex::Real(num_part[1]),
    Geom(level).ProbLength(2) / amrex::Real(num_part[2])));
  AMREX_D_TERM(amrex::ULong np0 = num_part[0];, amrex::ULong np1 = num_part[1];
               , amrex::ULong np2 = num_part[2];);
  const amrex::ULong total_part_num = AMREX_D_TERM(np0, *np1, *np2);
  amrex::ULong parts_pp = total_part_num / amrex::ULong(NProcs);
  // Number of particles per processor to be initialized
  amrex::ULong cur_parts_pp = parts_pp;
  // Give any remaining particles to the last processor
  if (MyProc == NProcs - 1) {
    cur_parts_pp += (total_part_num % amrex::ULong(NProcs));
  }
  // Starting particle for this processor
  const amrex::ULong first_part = amrex::ULong(MyProc) * parts_pp;
  amrex::Gpu::HostVector<ParticleType> nparticles;
  amrex::Vector<amrex::Gpu::HostVector<amrex::Real>> nreals;
  if (NAR_SPR > 0) {
    nreals.resize(NAR_SPR);
  }
  for (amrex::ULong prc = 0; prc < cur_parts_pp; ++prc) {
    amrex::ULong cur_part = first_part + prc;
    amrex::IntVect indx = unflatten_particles(cur_part, num_part);
    ParticleType p;
    p.id() = ParticleType::NextID();
    p.cpu() = amrex::ParallelDescriptor::MyProc();
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      p.pos(dir) = (amrex::Real(indx[dir]) + 0.5) * dx_part[dir];
    }
    for (int n = 0; n < NSR_SPR; ++n) {
      p.rdata(n) = part_vals[n];
    }
    for (int n = 0; n < NAR_SPR; ++n) {
      nreals[n].push_back(part_vals[n]);
    }
    nparticles.push_back(p);
  }
  amrex::ParticleLocData pld;
  // Only copy particle data for certain processors at a time
  int NRchunk = NProcs / NRedist;
  for (int nr = 0; nr < NRedist; ++nr) {
    std::map<std::pair<int, int>, amrex::Gpu::HostVector<ParticleType>>
      host_particles;
    std::map<
      std::pair<int, int>,
      std::array<amrex::Gpu::HostVector<amrex::Real>, NAR_SPR>>
      host_real_attribs;
    if (m_verbose > 0) {
      amrex::Print() << "Redistributing from processor " << nr * NRchunk
                     << " to " << (nr + 1) * NRchunk - 1 << '\n';
    }
    for (int which = nr * NRchunk; which < (nr + 1) * NRchunk; ++which) {
      if (which == MyProc) {
        while (!nparticles.empty()) {
          // Retrieve the last particle entry and add it to host_particles
          ParticleType& p = nparticles.back();
          bool where = Where(p, pld);
          if (!where) {
            amrex::Abort("Bad particle");
          }
          std::pair<int, int> ind(pld.m_grid, pld.m_tile);
          host_particles[ind].push_back(p);
          for (int n = 0; n < NAR_SPR; ++n) {
            host_real_attribs[ind][n].push_back(nreals[n].back());
          }
          // Remove the particle just read
          nparticles.pop_back();
          for (int n = 0; n < NAR_SPR; ++n) {
            nreals[n].pop_back();
          }
        }
      } // if (which == MyProc)
    }   // for (int which ...
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
      for (int i = 0; i < NAR_SPR; ++i) {
        amrex::Gpu::copy(
          amrex::Gpu::hostToDevice,
          host_real_attribs[std::make_pair(grid, tile)][i].begin(),
          host_real_attribs[std::make_pair(grid, tile)][i].end(),
          dst_tile.GetStructOfArrays().GetRealData(i).begin() + old_size);
      }
    }
    Redistribute();
  } // for (int nr ...
  // Now copy over any remaining processors
  for (int which = NRedist * NRchunk; which < NProcs; ++which) {
    std::map<std::pair<int, int>, amrex::Gpu::HostVector<ParticleType>>
      host_particles;
    std::map<
      std::pair<int, int>,
      std::array<amrex::Gpu::HostVector<amrex::Real>, NAR_SPR>>
      host_real_attribs;
    if (m_verbose > 0) {
      amrex::Print() << "Redistributing from processor " << NRedist * NRchunk
                     << " to " << NProcs << '\n';
    }
    if (which == MyProc) {
      while (!nparticles.empty()) {
        // Retrieve the last particle entry and add it to host_particles
        ParticleType& p = nparticles.back();
        Where(p, pld);
        std::pair<int, int> ind(pld.m_grid, pld.m_tile);
        host_particles[ind].push_back(p);
        for (int n = 0; n < NAR_SPR; ++n) {
          host_real_attribs[ind][n].push_back(nreals[n].back());
        }
        // Remove the particle just read
        nparticles.pop_back();
        for (int n = 0; n < NAR_SPR; ++n) {
          nreals[n].pop_back();
        }
      }
    } // if (which == MyProc)
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
      for (int i = 0; i < NAR_SPR; ++i) {
        amrex::Gpu::copy(
          amrex::Gpu::hostToDevice,
          host_real_attribs[std::make_pair(grid, tile)][i].begin(),
          host_real_attribs[std::make_pair(grid, tile)][i].end(),
          dst_tile.GetStructOfArrays().GetRealData(i).begin() + old_size);
      }
    }
    Redistribute();
  } // for (int which ...
  amrex::Gpu::streamSynchronize();
}
#endif
