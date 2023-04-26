
#include "SprayParticles.H"
#include "Drag.H"
#include "SprayInterpolation.H"
#include "Transport.H"
#include "WallFunctions.H"
#include "TABBreakup.H"
#include "ReitzKHRT.H"
#include "WallFilm.H"
#ifdef AMREX_USE_EB
#include <AMReX_EBFArrayBox.H>
#endif

using namespace amrex;

void
SprayParticleContainer::init_bcs()
{
  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    if (
      phys_bc->lo(dir) == Symmetry || phys_bc->lo(dir) == SlipWall ||
      phys_bc->lo(dir) == NoSlipWall) {
      reflect_lo[dir] = true;
    } else {
      reflect_lo[dir] = false;
    }
    if (
      phys_bc->hi(dir) == Symmetry || phys_bc->hi(dir) == SlipWall ||
      phys_bc->hi(dir) == NoSlipWall) {
      reflect_hi[dir] = true;
    } else {
      reflect_hi[dir] = false;
    }
  }
}

void
SprayParticleContainer::moveKick(
  MultiFab& state,
  MultiFab& source,
  const int level,
  const Real& dt,
  const Real time,
  const bool isVirtualPart,
  const bool isGhostPart,
  const int state_ghosts,
  const int source_ghosts,
  pele::physics::transport::TransParm<
    pele::physics::EosType,
    pele::physics::TransportType> const* ltransparm,
  const Real spray_cfl_lev)
{
  bool do_move = false;
  moveKickDrift(
    state, source, level, dt, time, isVirtualPart, isGhostPart, state_ghosts,
    source_ghosts, do_move, ltransparm, spray_cfl_lev);
}

void
SprayParticleContainer::moveKickDrift(
  MultiFab& state,
  MultiFab& source,
  const int level,
  const Real& dt,
  const Real time,
  const bool isVirtualPart,
  const bool isGhostPart,
  const int state_ghosts,
  const int source_ghosts,
  const bool do_move,
  pele::physics::transport::TransParm<
    pele::physics::EosType,
    pele::physics::TransportType> const* ltransparm,
  const Real spray_cfl_lev)
{
  BL_PROFILE("ParticleContainer::moveKickDrift()");
  AMREX_ASSERT(level >= 0);

  // If there are no particles at this level
  if (level >= this->GetParticles().size()) {
    return;
  }

  updateParticles(
    level, state, source, dt, time, state_ghosts, source_ghosts, isVirtualPart,
    isGhostPart, do_move, ltransparm, spray_cfl_lev);
}

Real
SprayParticleContainer::estTimestep(int level, Real cfl) const
{
  BL_PROFILE("ParticleContainer::estTimestep()");
  Real dt = std::numeric_limits<Real>::max();
  if (level >= this->GetParticles().size() || m_sprayData->fixed_parts) {
    return -1.;
  }

  const auto dx = Geom(level).CellSizeArray();
  const auto dxi = Geom(level).InvCellSizeArray();
  {
    ReduceOps<ReduceOpMin> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    for (MyParConstIter pti(*this, level); pti.isValid(); ++pti) {
      const AoS& pbox = pti.GetArrayOfStructs();
      const ParticleType* pstruct = pbox().data();
      const Long n = pbox.numParticles();
      reduce_op.eval(
        n, reduce_data, [=] AMREX_GPU_DEVICE(const Long i) -> ReduceTuple {
          const ParticleType& p = pstruct[i];
          if (p.id() > 0) {
            const Real max_mag_vdx = amrex::max(AMREX_D_DECL(
              std::abs(p.rdata(SprayComps::pstateVel)) * dxi[0],
              std::abs(p.rdata(SprayComps::pstateVel + 1)) * dxi[1],
              std::abs(p.rdata(SprayComps::pstateVel + 2)) * dxi[2]));
            Real dt_part = (max_mag_vdx > 0.) ? (cfl / max_mag_vdx) : 1.E50;
            return dt_part;
          }
          return -1.;
        });
    }
    ReduceTuple hv = reduce_data.value();
    Real ldt_cpu = amrex::get<0>(hv);
    dt = amrex::min(dt, ldt_cpu);
  }
  ParallelDescriptor::ReduceRealMin(dt);
  // Check if the velocity of particles being injected is greater than existing
  // particle velocities
  if (m_injectVel > 0.) {
    dt = amrex::min(dt, cfl * dx[0] / m_injectVel);
  }
  return dt;
}

void
SprayParticleContainer::updateParticles(
  const int& level,
  MultiFab& state,
  MultiFab& source,
  const Real& flow_dt,
  const Real& /*time*/,
  const int state_ghosts,
  const int source_ghosts,
  const bool isVirt,
  const bool isGhost,
  const bool do_move,
  pele::physics::transport::TransParm<
    pele::physics::EosType,
    pele::physics::TransportType> const* ltransparm,
  const Real spray_cfl_lev)
{
  BL_PROFILE("SprayParticleContainer::updateParticles()");
  AMREX_ASSERT(OnSameGrids(level, state));
  AMREX_ASSERT(OnSameGrids(level, source));
  bool isActive = !(isVirt || isGhost);
  bool do_splash = (m_sprayData->do_splash && isActive && do_move);
  bool do_breakup = (m_sprayData->do_breakup > 0 && isActive && do_move);
  Real B0 = B0_KHRT;
  Real B1 = B1_KHRT;
  Real C3 = C3_KHRT;
  Real max_ppp = max_num_ppp;
  const auto dxiarr = this->Geom(level).InvCellSizeArray();
  const auto dxarr = this->Geom(level).CellSizeArray();
  const auto ploarr = this->Geom(level).ProbLoArray();
  const auto phiarr = this->Geom(level).ProbHiArray();
  const RealVect dxi(AMREX_D_DECL(dxiarr[0], dxiarr[1], dxiarr[2]));
  const RealVect dx(AMREX_D_DECL(dxarr[0], dxarr[1], dxarr[2]));
  const RealVect plo(AMREX_D_DECL(ploarr[0], ploarr[1], ploarr[2]));
  const RealVect phi(AMREX_D_DECL(phiarr[0], phiarr[1], phiarr[2]));
  const auto domain = this->Geom(level).Domain();
#ifdef AMREX_USE_EB
  const auto& factory =
    dynamic_cast<EBFArrayBoxFactory const&>(state.Factory());
  const auto* cellcent = &(factory.getCentroid());
  const auto* bndrycent = &(factory.getBndryCent());
  const auto* bndrynorm = &(factory.getBndryNormal());
  const auto* volfrac = &(factory.getVolFrac());
#endif
  IntVect bndry_lo; // Designation for boundary types
  IntVect bndry_hi; // 0 - Periodic, 1 - Reflective, -1 - Non-reflective
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    if (!this->Geom(level).isPeriodic(dir)) {
      if (reflect_lo[dir]) {
        bndry_lo[dir] = 1;
      } else {
        bndry_lo[dir] = -1;
      }
      if (reflect_hi[dir]) {
        bndry_hi[dir] = 1;
      } else {
        bndry_hi[dir] = -1;
      }
    } else {
      bndry_lo[dir] = 0;
      bndry_hi[dir] = 0;
    }
  }
  const Real vol = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);
  const Real inv_vol = 1. / vol;
  // If particle subcycling is being done, determine the number of subcycles
  // Note: this is different than the AMR subcycling
  Real sub_cfl = 0.5; // CFL for each subcycle
  Real sub_dt = flow_dt;
  int num_iter = 1;
  if (do_move && spray_cfl_lev > sub_cfl) {
    num_iter = static_cast<int>(std::ceil(spray_cfl_lev / sub_cfl));
    sub_dt = flow_dt / static_cast<Real>(num_iter);
  }
  Real avg_inject_d3 = 0.;
  if (isActive && m_sprayData->do_breakup == 2) {
    int numJets = static_cast<int>(m_sprayJets.size());
    for (int jindx = 0; jindx < numJets; ++jindx) {
      Real injDia = m_sprayJets[jindx].get()->get_avg_dia();
      avg_inject_d3 += std::pow(injDia, 3) / static_cast<Real>(numJets);
    }
  }
  // Particle components indices
  SprayComps SPI = m_sprayIndx;
  // Start the ParIter, which loops over separate sets of particles in different
  // boxes
  for (MyParIter pti(*this, level); pti.isValid(); ++pti) {
    const Box tile_box = pti.tilebox();
    const Box src_box = pti.growntilebox(source_ghosts);
    const Box state_box = pti.growntilebox(state_ghosts);
    bool at_bounds = tile_at_bndry(tile_box, bndry_lo, bndry_hi, domain);
    const Long Np = pti.numParticles();
    AoS& pbox = pti.GetArrayOfStructs();
    ParticleType* pstruct = pbox().data();
    const SprayData* fdat = d_sprayData;
    Array4<const Real> const& Tarr = state.array(pti, SPI.utempIndx);
    Array4<const Real> const& rhoYarr = state.array(pti, SPI.specIndx);
    Array4<const Real> const& rhoarr = state.array(pti, SPI.rhoIndx);
    Array4<const Real> const& momarr = state.array(pti, SPI.momIndx);
    Array4<const Real> const& engarr = state.array(pti, SPI.engIndx);
    Array4<Real> const& rhoYSrcarr = source.array(pti, SPI.specSrcIndx);
    Array4<Real> const& rhoSrcarr = source.array(pti, SPI.rhoSrcIndx);
    Array4<Real> const& momSrcarr = source.array(pti, SPI.momSrcIndx);
    Array4<Real> const& engSrcarr = source.array(pti, SPI.engSrcIndx);
    bool eb_in_box = false;

#ifdef AMREX_USE_EB
    eb_in_box = true;
    const auto& interp_fab = static_cast<EBFArrayBox const&>(state[pti]);
    const EBCellFlagFab& flags = interp_fab.getEBCellFlagFab();
    Array4<const Real> ccent_fab;
    Array4<const Real> bcent_fab;
    Array4<const Real> bnorm_fab;
    Array4<const Real> volfrac_fab;
    const auto& flags_array = flags.array();
    if (flags.getType(state_box) == FabType::regular) {
      eb_in_box = false;
    } else {
      // Cell centroids
      ccent_fab = cellcent->array(pti);
      // Centroid of EB
      bcent_fab = bndrycent->array(pti);
      // Normal of EB
      bnorm_fab = bndrynorm->array(pti);
      volfrac_fab = volfrac->array(pti);
    }
#endif
    amrex::ParallelFor(Np, [=] AMREX_GPU_DEVICE(int pid) noexcept {
      ParticleType& p = pstruct[pid];
      if (p.id() > 0) {
        auto eos = pele::physics::PhysicsType::eos();
        SprayUnits SPU;
        GasPhaseVals gpv;
        eos.molecular_weight(gpv.mw.data());
        for (int n = 0; n < NUM_SPECIES; ++n) {
          gpv.mw[n] *= SPU.mass_conv;
        }
        GpuArray<IntVect, AMREX_D_PICK(2, 4, 8)>
          indx_array; // array of adjacent cells
        GpuArray<Real, AMREX_D_PICK(2, 4, 8)>
          weights; // array of corresponding weights
        RealVect lx = (p.pos() - plo) * dxi + 0.5;
        IntVect ijk = lx.floor(); // Upper cell center
        RealVect lxc = (p.pos() - plo) * dxi;
        IntVect ijkc = lxc.floor(); // Cell with particle
        IntVect ijkc_prev = ijkc;
        IntVect bflags(IntVect::TheZeroVector());
        if (at_bounds) {
          // Check if particle has left the domain or is boundary adjacent
          bool left_dom =
            check_bounds(p.pos(), plo, phi, dx, bndry_lo, bndry_hi, bflags);
          if (left_dom) {
            Abort("Particle has incorrectly left the domain");
          }
        }
        // Used for ETAB breakup model
        Real Utan_total = 0.;
        Real Reyn_d = 0.;
        bool is_film = false;
        RealVect film_normal;
        // Gather wall film values
        if (p.rdata(SprayComps::pstateFilmHght) > 0.) {
          is_film = true;
        }
        // Subcycle loop
        int cur_iter = 0;
        while (p.id() > 0 && cur_iter < num_iter) {
          // Flag for whether we are near EB boundaries
          bool do_fe_interp = false;
          if (is_film) {
            do_fe_interp = interpolateFilm(
              ijkc, bflags, film_normal,
#ifdef AMREX_USE_EB
              eb_in_box, flags_array, bnorm_fab,
#endif
              indx_array.data(), weights.data());
          }
#ifdef AMREX_USE_EB
          else if (eb_in_box) {
            do_fe_interp = eb_interp(
              p, isVirt, ijkc, ijk, dx, dxi, lx, plo, bflags, flags_array,
              ccent_fab, bcent_fab, bnorm_fab, volfrac_fab, fdat->min_eb_vfrac,
              indx_array.data(), weights.data());
          }
#endif
          else {
            trilinear_interp(
              ijk, lx, indx_array.data(), weights.data(), bflags);
          }
          // Interpolate fluid state
          gpv.reset();
          InterpolateGasPhase(
            gpv, state_box, rhoarr, rhoYarr, Tarr, momarr, engarr,
            indx_array.data(), weights.data());
          // Solve for avg mw and pressure at droplet location
          gpv.define();
          fdat->calcBoilT(gpv, cBoilT.data());
          if (is_film) {
            calculateFilmSource(
              sub_dt, gpv, *fdat, p, cBoilT.data(), film_normal, ltransparm);
          } else {
            Reyn_d = calculateSpraySource(
              sub_dt, gpv, *fdat, p, cBoilT.data(), ltransparm);
          }
          IntVect cur_indx = ijkc;
          Real cvol = inv_vol;
          if (p.id() > 0 && do_breakup) {
            // Update breakup variables and determine if breakup occurs
            if (fdat->do_breakup == 1) {
              Utan_total +=
                updateBreakupTAB(Reyn_d, sub_dt, cBoilT.data(), gpv, *fdat, p);
            }
            if (cur_iter == num_iter - 1) {
              if (fdat->do_breakup == 1) {
                // Determine if parcel must be split into multiple parcels
                splitDropletTAB(pid, p, max_ppp, N_SB, rf_d, Utan_total);
              } else {
                // Update breakup for KH-RT model
                updateBreakupKHRT(
                  pid, p, Reyn_d, flow_dt, cBoilT.data(), avg_inject_d3, B0, B1,
                  C3, gpv, *fdat, N_SB, rf_d);
              }
            }
          }
#ifdef AMREX_USE_EB
          if (flags_array(cur_indx).isSingleValued()) {
            if (volfrac_fab(cur_indx) < fdat->min_eb_vfrac) {
              Real min_dis = 1.E12;
              for (int aindx = 0; aindx < AMREX_D_PICK(2, 4, 8); ++aindx) {
                IntVect cindx = indx_array[aindx];
                if (volfrac_fab(cindx) > fdat->min_eb_vfrac) {
                  Real dis = 0.;
                  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                    dis += std::pow(
                      p.pos(dir) - plo[dir] -
                        static_cast<Real>(cindx[dir]) * dx[dir],
                      2);
                  }
                  if (dis < min_dis) {
                    min_dis = dis;
                    cur_indx = cindx;
                  }
                }
              }
            }
            cvol *= 1. / (volfrac_fab(cur_indx));
          }
#endif
          Real cur_coef = -cvol * sub_dt / flow_dt;
          if (!src_box.contains(cur_indx)) {
            if (!isGhost) {
              Abort("SprayParticleContainer::updateParticles() -- source box "
                    "too small");
            }
          }
          if (fdat->mom_trans) {
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
              Gpu::Atomic::Add(
                &momSrcarr(cur_indx, dir), cur_coef * gpv.fluid_mom_src[dir]);
            }
          }
          if (fdat->mass_trans) {
            Gpu::Atomic::Add(
              &rhoSrcarr(cur_indx), cur_coef * gpv.fluid_mass_src);
            for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
              Gpu::Atomic::Add(
                &rhoYSrcarr(cur_indx, spf), cur_coef * gpv.fluid_Y_dot[spf]);
            }
          }
          Gpu::Atomic::Add(&engSrcarr(cur_indx), cur_coef * gpv.fluid_eng_src);
          // Modify particle position by whole time step
          if (do_move && !fdat->fixed_parts && p.id() > 0 && !is_film) {
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
              const Real cvel = p.rdata(SprayComps::pstateVel + dir);
              p.pos(dir) += sub_dt * cvel;
            }
            if (at_bounds || do_fe_interp) {
              // First check if particle has exited the domain through a
              // Cartesian boundary
              bool left_dom =
                check_bounds(p.pos(), plo, phi, dx, bndry_lo, bndry_hi, bflags);
              if (left_dom) {
                p.id() = -1;
              } else {
                // TODO: Add methods for determining this
                Real film_h = 0.;
                // Next reflect particles off BC or EB walls if necessary
                impose_wall(
                  p, dx, plo, phi, bndry_lo, bndry_hi, bflags, eb_in_box,
#ifdef AMREX_USE_EB
                  flags_array, bcent_fab, bnorm_fab, volfrac_fab,
                  fdat->min_eb_vfrac,
#endif
                  ijkc, N_SB, rf_d, film_h);
              }
            } // if (at_bounds || fe_interp)
            // Update indices
            lx = (p.pos() - plo) * dxi + 0.5;
            ijk = lx.floor();
            lxc = (p.pos() - plo) * dxi;
            ijkc = lxc.floor(); // New cell center
          }
          if (isGhost && !src_box.contains(ijkc)) {
            p.id() = -1;
          }
        } // End of subcycle loop
      }   // End of p.id() > 0 check
    });   // End of loop over particles
    if (do_splash_box || do_breakup) {
      Gpu::copy(
        Gpu::deviceToHost, N_SB_d.begin(), N_SB_d.end(), N_SB_h.begin());
      bool get_new_parts = false;
      for (Long n = 0; n < Np; n++) {
        if (N_SB_h[n] != splash_breakup::no_change) {
          get_new_parts = true;
        }
      }
      if (get_new_parts) {
        refv.retrieve_data();
        SBPtrs rfh;
        refv.fillPtrs_h(rfh);
        CreateSBDroplets(Np, N_SB_h.data(), rfh, level);
      }
    }
  } // for (int MyParIter pti..
}
