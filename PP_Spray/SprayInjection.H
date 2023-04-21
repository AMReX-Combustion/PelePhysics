#ifndef SPRAYINJECTION_H
#define SPRAYINJECTION_H
#include "SprayParticles.H"

/*
 These are generalized initialization and injection routines to be called in a
 problem-specific SprayParticlesInitInsert.cpp file, when necessary.
*/

void
SprayParticleContainer::sprayInjection(
  const amrex::Real time,
  SprayJet* spray_jet,
  const amrex::Real sim_dt,
  const int level)
{
  if (!spray_jet->jet_active(time) || sim_dt <= 0.) {
    return;
  }
  int curProc = amrex::ParallelDescriptor::MyProc();
  int injProc = spray_jet->Proc();
  if (curProc != injProc) {
    return;
  }
  SprayUnits SPU;
  const SprayData* fdat = m_sprayData;
  amrex::Real rho_avg = 0.;
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    rho_avg +=
      spray_jet->get_avg_Y(spf) / fdat->rhoL(spray_jet->get_avg_T(), spf);
  }
  rho_avg = 1. / rho_avg;
  const amrex::Real Pi_six = M_PI / 6.;
  // Number of particles per parcel
  amrex::Real num_ppp = spray_jet->num_ppp();
  const amrex::Real min_ppp = 0.01;
  // Check if mass must be injected across multiple timesteps
  const amrex::Real avg_dia = spray_jet->get_avg_dia();
  const amrex::Real avg_mass = Pi_six * rho_avg * std::pow(avg_dia, 3);
  const amrex::Real min_dia = std::cbrt(6. * SPU.min_mass / (M_PI * rho_avg));
  if (avg_dia < min_dia || avg_mass < SPU.min_mass) {
    amrex::Abort(
      "Average droplet size too small, floating point issues expected");
  }
  // Minimum number of parcels we want injected at a time
  amrex::Real min_inj_parcel = spray_jet->m_minParcel;
  amrex::Real dt = sim_dt;
  amrex::Real mdot = spray_jet->mass_flow_rate();
  amrex::Real inject_mass = mdot * sim_dt;
  // See if there is uninjected mass from previous time steps
  if (spray_jet->m_sumInjMass > 0.) {
    dt += spray_jet->m_sumInjTime;
    inject_mass += spray_jet->m_sumInjMass;
  }
  // If we can't get enough mass to create the desired minimum parcel injection,
  // accumulate mass injection over more time steps
  if (inject_mass / (num_ppp * avg_mass) < min_inj_parcel) {
    spray_jet->m_sumInjMass = inject_mass;
    spray_jet->m_sumInjTime = dt;
    return;
  }

  amrex::ParticleLocData pld;
  std::map<std::pair<int, int>, amrex::Gpu::HostVector<ParticleType>>
    host_particles;
  amrex::Real cur_mass = 0.;
  while (cur_mass < inject_mass) {
    // Pick random percentage from 0 to 1
    amrex::Real radp = amrex::Random();
#if AMREX_SPACEDIM == 3
    if (spray_jet->hollow_spray()) {
      radp = 1.;
    }
    amrex::Real phi_radial = amrex::Random() * 2. * M_PI;
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
    // phi_radial - azimuthal angle, called phi_j in docs
    // theta_spread - inclination angle, called theta_j in docs
    // phi_swirl - azimuthal swirling component, called phi_s in docs
    amrex::Real umag, theta_spread, phi_swirl, dia_part, T_part;
    amrex::Real Y_part[SPRAY_FUEL_NUM];
    bool make_new_part = spray_jet->get_new_particle(
      time, phi_radial, cur_rad, umag, theta_spread, phi_swirl, dia_part,
      T_part, Y_part);
    if (make_new_part && dia_part > min_dia) {
      amrex::RealVect part_loc, vel_part;
      spray_jet->transform_loc_vel(
        theta_spread, phi_radial, cur_rad, umag, phi_swirl, vel_part, part_loc);
      ParticleType p;
      p.id() = ParticleType::NextID();
      p.cpu() = amrex::ParallelDescriptor::MyProc();
      AMREX_D_TERM(p.rdata(SprayComps::pstateVel) = part_vel[0];
                   , p.rdata(SprayComps::pstateVel + 1) = part_vel[1];
                   , p.rdata(SprayComps::pstateVel + 2) = part_vel[2];);
      p.rdata(SprayComps::pstateT) = T_part;
      // Never add particle with less than minimum mass
      p.rdata(SprayComps::pstateDia) = dia_part;
      amrex::Real rho_part = 0.;
      if (SPRAY_FUEL_NUM > 1) {
        for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
          p.rdata(SprayComps::pstateY + spf) = Y_part[spf];
          rho_part += Y_part[spf] / fdat->rhoL(T_part, spf);
        }
        rho_part = 1. / rho_part;
      } else {
        rho_part = fdat->rhoL(T_part, 0);
        p.rdata(SprayComps::pstateY) = 1.;
      }
      // Add particles as if they have advanced some random portion of
      // dt
      amrex::Real pmov = amrex::Random();
      for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        p.pos(dir) =
          part_loc[dir] + pmov * dt * p.rdata(SprayComps::pstateVel + dir);
      }
      // If KHRT breakup is used, BM1 is starting diameter
      if (m_sprayData->do_breakup == 2) {
        p.rdata(SprayComps::pstateBM1) = dia_part;
      } else {
        // If TAB is used, BM1 is y
        p.rdata(SprayComps::pstateBM1) = 0.;
      }
      // If KHRT is used, BM2 is RT time
      // If TAB is used, BM2 is y dot
      p.rdata(SprayComps::pstateBM2) = 0.;
      p.rdata(SprayComps::pstateFilmHght) = 0.;
      p.rdata(SprayComps::pstateNumDens) = num_ppp;
      amrex::Real pmass = Pi_six * rho_part * std::pow(dia_part, 3);
      amrex::Real new_mass = cur_mass + num_ppp * pmass;
      bool where = Where(p, pld);
      if (!where) {
        amrex::Abort("Bad injection particle");
      }
      std::pair<int, int> ind(pld.m_grid, pld.m_tile);
      host_particles[ind].push_back(p);
      cur_mass = new_mass;
    }
  }
  amrex::Real est_mdot = cur_mass / dt;
  // If we are over-injecting mass, we can either decrease the number density or increase the minimum parcel
  if (est_mdot - mdot > 0.1 * mdot) {
    // Increase the minimum injection parcels
    spray_jet->m_minParcel += 1.;
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
  const int num_redist)
{
  const int MyProc = amrex::ParallelDescriptor::MyProc();
  const int NProcs = amrex::ParallelDescriptor::NProcs();
  int NRedist = num_redist;
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
  // Reference values for the particles
  amrex::Real part_vals[NSR_SPR];
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    part_vals[SprayComps::pstateVel + dir] = vel_part[dir];
  }
  part_vals[SprayComps::pstateT] = T_part;
  part_vals[SprayComps::pstateDia] = dia_part;
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    part_vals[SprayComps::pstateY + spf] = Y_part[spf];
  }
  part_vals[SprayComps::pstateNumDens] = num_ppp;
  if (m_sprayData->do_breakup == 2) {
    part_vals[SprayComps::pstateBM1] = dia_part;
  } else {
    part_vals[SprayComps::pstateBM1] = 0.;
  }
  part_vals[SprayComps::pstateBM2] = 0.;
  part_vals[SprayComps::pstateFilmHght] = 0.;
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
    nparticles.push_back(p);
  }
  amrex::ParticleLocData pld;
  // Only copy particle data for certain processors at a time
  int NRchunk = NProcs / NRedist;
  for (int nr = 0; nr < NRedist; ++nr) {
    std::map<std::pair<int, int>, amrex::Gpu::HostVector<ParticleType>>
      host_particles;
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
          // Remove the particle just read
          nparticles.pop_back();
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
    }
    Redistribute();
  } // for (int nr ...
  // Now copy over any remaining processors
  for (int which = NRedist * NRchunk; which < NProcs; ++which) {
    std::map<std::pair<int, int>, amrex::Gpu::HostVector<ParticleType>>
      host_particles;
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
        // Remove the particle just read
        nparticles.pop_back();
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
    }
    Redistribute();
  } // for (int which ...
  amrex::Gpu::streamSynchronize();
}
#endif
