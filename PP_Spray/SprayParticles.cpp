
#include "SprayParticles.H"
#include <AMReX_ParticleReduce.H>
#include <AMReX_Particles.H>
#ifdef SPRAY_PELE_LM
#include "PeleLM.H"
#endif
#include "Drag.H"
#include "SprayInterpolation.H"
#include "Transport.H"
#include "WallFunctions.H"
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
  MultiFab* u_mac)
{
  bool do_move = false;
  int width = 0;
  moveKickDrift(
    state, source, level, dt, time, isVirtualPart, isGhostPart, state_ghosts,
    source_ghosts, do_move, width, u_mac);
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
  const int /*where_width*/,
  MultiFab* u_mac)
{
  BL_PROFILE("ParticleContainer::moveKickDrift()");
  AMREX_ASSERT(u_mac == nullptr || u_mac[0].nGrow() >= 1);
  AMREX_ASSERT(level >= 0);
  AMREX_ASSERT(state.nGrow() >= 2);

  // If there are no particles at this level
  if (level >= this->GetParticles().size())
    return;

  bool isActive = (isVirtualPart || isGhostPart) ? false : true;

  BL_PROFILE_VAR("SprayParticles::updateParticles()", UPD_PART);
  updateParticles(
    level, state, source, dt, time, state_ghosts, source_ghosts, isActive,
    do_move, u_mac);
  BL_PROFILE_VAR_STOP(UPD_PART);

  // Fill ghost cells after we've synced up ..
  // TODO: Check to see if this is needed at all
  // if (level > 0)
  //   source.FillBoundary(Geom(level).periodicity());

  // ********************************************************************************

  // ********************************************************************************
}

Real
SprayParticleContainer::estTimestep(int level, Real cfl) const
{
  BL_PROFILE("ParticleContainer::estTimestep()");
  // TODO: Clean up this mess and bring the num particle functionality back
  Real dt = std::numeric_limits<Real>::max();
  if (level >= this->GetParticles().size() || m_sprayIndx.mom_tran == 0)
    return -1.;

  const Geometry& geom = this->m_gdb->Geom(level);
  const auto dx = Geom(level).CellSizeArray();
  const auto dxi = Geom(level).InvCellSizeArray();
  {
    amrex::ReduceOps<amrex::ReduceOpMin> reduce_op;
    amrex::ReduceData<amrex::Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    for (MyParConstIter pti(*this, level); pti.isValid(); ++pti) {
      const AoS& pbox = pti.GetArrayOfStructs();
      const ParticleType* pstruct = pbox().data();
      const Long n = pbox.numParticles();
#ifdef USE_SPRAY_SOA
      auto& attribs = pti.GetAttribs();
      AMREX_D_TERM(const Real* up = attribs[0].data();
                   , const Real* vp = attribs[1].data();
                   , const Real* wp = attribs[2].data(););
#endif
      reduce_op.eval(
        n, reduce_data, [=] AMREX_GPU_DEVICE(const Long i) -> ReduceTuple {
          const ParticleType& p = pstruct[i];
          // TODO: This assumes that pstateVel = 0 and dxi[0] = dxi[1] =
          // dxi[2]
          if (p.id() > 0) {
#ifdef USE_SPRAY_SOA
            const Real max_mag_vdx =
              amrex::max(AMREX_D_DECL(
                std::abs(up[i]), std::abs(vp[i]), std::abs(wp[i]))) *
              dxi[0];
#else
          const Real max_mag_vdx = amrex::max(AMREX_D_DECL(std::abs(p.rdata(0)),
                                                           std::abs(p.rdata(1)),
                                                           std::abs(p.rdata(2))))*dxi[0];
#endif
            Real dt_part = (max_mag_vdx > 0.) ? (cfl / max_mag_vdx) : 1.E50;
#ifdef SPRAY_PELE_LM
            // Conversion since particle velocities are in cm and dx is in m
            // for PeleLM
            dt_part *= 100.;
#endif
            return dt_part;
          }
          return 1.E50;
        });
    }
    ReduceTuple hv = reduce_data.value();
    Real ldt_cpu = amrex::get<0>(hv);
    dt = amrex::min(dt, ldt_cpu);
  }
  ParallelDescriptor::ReduceRealMin(dt);
  // Check if the velocity of particles being injected
  // is greater existing particle velocities
  if (m_injectVel > 0.)
    dt = amrex::min(dt, cfl * dx[0] / m_injectVel);

  return dt;
}

void
SprayParticleContainer::updateParticles(
  const int& level,
  MultiFab& state,
  MultiFab& source,
  const Real& flow_dt,
  const Real& time,
  const int state_ghosts,
  const int source_ghosts,
  const bool isActive,
  const bool do_move,
  MultiFab* u_mac)
{
  AMREX_ASSERT(OnSameGrids(level, state));
  AMREX_ASSERT(OnSameGrids(level, source));
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
  const auto& flagmf = factory.getMultiEBCellFlagFab();
  const auto cellcent = &(factory.getCentroid());
  const auto bndrycent = &(factory.getBndryCent());
  const auto bndryarea = &(factory.getBndryArea());
  const auto bndrynorm = &(factory.getBndryNormal());
  const auto volfrac = &(factory.getVolFrac());
#endif
  IntVect bndry_lo; // Designation for boundary types
  IntVect bndry_hi; // 0 - Periodic, 1 - Reflective, -1 - Non-reflective
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    if (!this->Geom(level).isPeriodic(dir)) {
      if (reflect_lo[dir])
        bndry_lo[dir] = 1;
      else
        bndry_lo[dir] = -1;
      if (reflect_hi[dir])
        bndry_hi[dir] = 1;
      else
        bndry_hi[dir] = -1;
    } else {
      bndry_lo[dir] = 0;
      bndry_hi[dir] = 0;
    }
  }
  const Real wallT = m_wallT;
  const Real vol = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);
  const Real inv_vol = 1. / vol;
  // Particle components indices
  SprayComps SPI = m_sprayIndx;
  SprayUnits SPU;
  pele::physics::transport::TransParm const* ltransparm =
    pele::physics::transport::trans_parm_g;
  // Start the ParIter, which loops over separate sets of particles in different
  // boxes
  for (MyParIter pti(*this, level); pti.isValid(); ++pti) {
    const Box tile_box = pti.tilebox();
    const Box state_box = pti.growntilebox(state_ghosts);
    const Box src_box = pti.growntilebox(source_ghosts);
    bool at_bounds = tile_at_bndry(tile_box, bndry_lo, bndry_hi, domain);
    const Long Np = pti.numParticles();
    ParticleType* pstruct = &(pti.GetArrayOfStructs()[0]);
#ifdef USE_SPRAY_SOA
    // Get particle attributes if StructOfArrays are used
    auto& attribs = pti.GetAttribs();
#endif
    const SprayData* fdat = d_sprayData;
    Array4<const Real> const& statearr = state.array(pti);
    Array4<Real> const& sourcearr = source.array(pti);
#ifdef AMREX_USE_EB
    bool eb_in_box = true;
    const EBFArrayBox& interp_fab = static_cast<EBFArrayBox const&>(state[pti]);
    const EBCellFlagFab& flags = interp_fab.getEBCellFlagFab();
    Array4<const Real> ccent_fab;
    Array4<const Real> bcent_fab;
    Array4<const Real> bnorm_fab;
    Array4<const Real> barea_fab;
    Array4<const Real> volfrac_fab;
    const auto& flags_array = flags.array();
    if (flags.getType(state_box) == FabType::regular) {
      eb_in_box = false;
    } else {
      // Cell centroids
      ccent_fab = cellcent->array(pti);
      // Centroid of EB
      bcent_fab = bndrycent->array(pti);
      // Area of EB face
      barea_fab = bndryarea->array(pti);
      // Normal of EB
      bnorm_fab = bndrynorm->array(pti);
      volfrac_fab = volfrac->array(pti);
    }
#endif
    // #ifdef SPRAY_PELE_LM
    //     GpuArray<
    //       Array4<const Real>, AMREX_SPACEDIM> const
    //       umac{AMREX_D_DECL(u_mac[0].array(pti), u_mac[1].array(pti),
    //       u_mac[2].array(pti))};
    // #endif
    amrex::ParallelFor(
      Np, [pstruct, statearr, sourcearr, plo, phi, dx, dxi, do_move, SPI, SPU,
           fdat, src_box, state_box, bndry_hi, bndry_lo, flow_dt, inv_vol,
           ltransparm, at_bounds, wallT, isActive
#ifdef USE_SPRAY_SOA
           ,
           attribs
#endif
#ifdef AMREX_USE_EB
           ,
           flags_array, ccent_fab, bcent_fab, bnorm_fab, barea_fab, volfrac_fab,
           eb_in_box
#endif
    ] AMREX_GPU_DEVICE(int pid) noexcept {
        auto eos = pele::physics::PhysicsType::eos();
        GpuArray<Real, NUM_SPECIES> mw_fluid;
        GpuArray<Real, NUM_SPECIES> invmw;
        eos.molecular_weight(mw_fluid.data());
        eos.inv_molecular_weight(invmw.data());
        ParticleType& p = pstruct[pid];
        if (p.id() > 0) {
          GpuArray<IntVect, AMREX_D_PICK(2, 4, 8)>
            indx_array; // Array of adjacent cells
          GpuArray<Real, AMREX_D_PICK(2, 4, 8)>
            weights; // Array of corresponding weights
          bool remove_particle = false;
          RealVect lx = (p.pos() - plo) * dxi + 0.5;
          IntVect ijk = lx.floor(); // Upper cell center
          bool is_wall_film = false;
          Real face_area = 0.;
          // If the temperature is a negative value, the particle is a wall film
          if (p.rdata(SPI.pstateT) < 0.)
            is_wall_film = true;
          if (at_bounds) {
            IntVect bflags(IntVect::TheZeroVector());
            // Check if particle has left the domain, is wall film,
            // or is boundary adjacent and must be shifted
            bool left_dom = check_bounds(
              p.pos(), plo, phi, dx, bndry_lo, bndry_hi, ijk, bflags);
            if (left_dom)
              Abort("Particle has incorrectly left domain");
            if (is_wall_film) {
              // TODO: Assumes grid spacing is uniform in all directions
              face_area = AMREX_D_TERM(1., *dx[0], *dx[0]);
            }
          }
          // Length from cell center to boundary face center
          Real diff_cent = 0.5 * dx[0];
          bool do_fe_interp = true;
#ifdef AMREX_USE_EB
          // Cell containing particle centroid
          AMREX_D_TERM(const int ip = static_cast<int>(
                         amrex::Math::floor((p.pos(0) - plo[0]) * dxi[0]));
                       , const int jp = static_cast<int>(
                           amrex::Math::floor((p.pos(1) - plo[1]) * dxi[1]));
                       , const int kp = static_cast<int>(
                           amrex::Math::floor((p.pos(2) - plo[2]) * dxi[2])););
          AMREX_D_TERM(const int i = static_cast<int>(amrex::Math::floor(
                         (p.pos(0) - plo[0]) * dxi[0] + 0.5));
                       , const int j = static_cast<int>(amrex::Math::floor(
                           (p.pos(1) - plo[1]) * dxi[1] + 0.5));
                       , const int k = static_cast<int>(amrex::Math::floor(
                           (p.pos(2) - plo[2]) * dxi[2] + 0.5)););
          if (!eb_in_box) {
            do_fe_interp = false;
          } else {
            // All cells in the stencil are regular. Use
            // traditional trilinear interpolation
            if (
              flags_array(i - 1, j - 1, k - 1).isRegular() and
              flags_array(i, j - 1, k - 1).isRegular() and
              flags_array(i - 1, j, k - 1).isRegular() and
              flags_array(i, j, k - 1).isRegular() and
              flags_array(i - 1, j - 1, k).isRegular() and
              flags_array(i, j - 1, k).isRegular() and
              flags_array(i - 1, j, k).isRegular() and
              flags_array(i, j, k).isRegular())
              do_fe_interp = false;
          }
          if (do_fe_interp) {
            fe_interp(
              p.pos(), ip, jp, kp, dx, dxi, plo, flags_array, ccent_fab,
              bcent_fab, bnorm_fab, volfrac_fab, indx_array.data(),
              weights.data());
          } else {
            trilinear_interp(ijk, lx, indx_array.data(), weights.data());
          }
          bool eb_wall_film = false;
          if (is_wall_film && flags_array(ip, jp, kp).isSingleValued()) {
            eb_wall_film = true;
            face_area = barea_fab(ip, jp, kp);
            diff_cent = 0.;
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
              diff_cent += std::pow(
                bcent_fab(ip, jp, kp, dir) - ccent_fab(ip, jp, kp, dir), 2);
            diff_cent = std::sqrt(diff_cent);
          }
#else
          trilinear_interp(ijk, lx, indx_array.data(), weights.data());
#endif // AMREX_USE_EB
       // Interpolate fluid state
          Real T_fluid = 0.;
          Real rho_fluid = 0.;
          GpuArray<Real, NUM_SPECIES> Y_fluid;
          for (int n = 0; n < NUM_SPECIES; ++n)
            Y_fluid[n] = 0.;
          RealVect vel_fluid(RealVect::TheZeroVector());
          {
            GpuArray<Real, NUM_SPECIES> mass_frac;
            for (int aindx = 0.; aindx < AMREX_D_PICK(2, 4, 8); ++aindx) {
              IntVect cur_indx = indx_array[aindx];
              Real cw = weights[aindx];
#ifdef AMREX_DEBUG
              if (!state_box.contains(cur_indx))
                Abort("SprayParticleContainer::updateParticles() -- state box "
                      "too small");
#endif
              Real cur_rho = statearr(cur_indx, SPI.rhoIndx);
              rho_fluid += cw * cur_rho;
              Real inv_rho = 1. / cur_rho;
              for (int n = 0; n < NUM_SPECIES; ++n) {
                int sp = n + SPI.specIndx;
                Real cur_mf = statearr(cur_indx, sp) * inv_rho;
                Y_fluid[n] += cw * cur_mf;
                mass_frac[n] = cur_mf;
              }
#ifdef SPRAY_PELE_LM
              inv_rho = 1.;
#endif
              Real ke = 0.;
              for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                Real vel = statearr(cur_indx, SPI.momIndx + dir) * inv_rho *
                           SPU.vel_conv;
                vel_fluid[dir] += cw * vel;
                ke += vel * vel / 2.;
              }
              Real T_i = statearr(cur_indx, SPI.utempIndx);
#ifndef SPRAY_PELE_LM
              Real intEng = statearr(cur_indx, SPI.engIndx) * inv_rho - ke;
              eos.EY2T(intEng, mass_frac.data(), T_i);
#endif
              T_fluid += cw * T_i;
            }
            rho_fluid *= SPU.rho_conv;
          }
          GasPhaseVals gpv(
            vel_fluid, T_fluid, rho_fluid, Y_fluid.data(), mw_fluid.data(),
            invmw.data());
          if (!is_wall_film) {
            remove_particle = calculateSpraySource(
              flow_dt, do_move, gpv, SPI, *fdat, p,
#ifdef USE_SPRAY_SOA
              attribs, pid,
#endif
              ltransparm);
            // Modify particle position by whole time step
            if (do_move && SPI.mom_tran) {
              for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
#ifdef USE_SPRAY_SOA
                const Real cvel = attribs[SPI.pstateVel + dir].data()[pid];
#else
                const Real cvel = p.rdata(SPI.pstateVel+dir);
#endif
                Gpu::Atomic::Add(&p.pos(dir), flow_dt * cvel * SPU.pos_conv);
              }
            }
          } else {
            remove_particle = calculateWallFilmSource(
              flow_dt, gpv, SPI, *fdat, p,
#ifdef USE_SPRAY_SOA
              attribs, pid,
#endif
              wallT, face_area, diff_cent, ltransparm);
          }
          if (remove_particle)
            p.id() = -1;
          for (int aindx = 0; aindx < AMREX_D_PICK(2, 4, 8); ++aindx) {
            IntVect cur_indx = indx_array[aindx];
            Real cvol = inv_vol;
#ifdef AMREX_USE_EB
            if (!flags_array(cur_indx).isRegular())
              cvol *= 1. / (volfrac_fab(cur_indx));
#endif
            Real cur_coef = -weights[aindx] * fdat->num_ppp * cvol;
#ifdef AMREX_DEBUG
            if (!src_box.contains(cur_indx))
              Abort(
                "SprayParticleContainer::updateParticles() -- source box too "
                "small");
#endif
            if (SPI.mom_tran) {
              for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                const int nf = SPI.momIndx + dir;
                Gpu::Atomic::Add(
                  &sourcearr(cur_indx, nf),
                  cur_coef * gpv.fluid_mom_src[dir] * SPU.mom_src_conv);
              }
            }
            if (SPI.mass_tran) {
              Gpu::Atomic::Add(
                &sourcearr(cur_indx, SPI.rhoIndx),
                cur_coef * gpv.fluid_mass_src * SPU.mass_src_conv);
              for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
                const int nf = SPI.specIndx + fdat->indx[spf];
                Gpu::Atomic::Add(
                  &sourcearr(cur_indx, nf),
                  cur_coef * gpv.fluid_Y_dot[spf] * SPU.mass_src_conv);
              }
            }
            Gpu::Atomic::Add(
              &sourcearr(cur_indx, SPI.engIndx),
              cur_coef * gpv.fluid_eng_src * SPU.eng_src_conv);
          }
          // Solve for splash model/wall film formation
          if (at_bounds && do_move) {
            IntVect ijk_prev = ijk;
            lx = (p.pos() - plo) * dxi + 0.5;
            ijk = lx.floor();
            const Real T_part = p.rdata(SPI.pstateT);
            const Real dia_part = p.rdata(SPI.pstateDia);
            IntVect bloc(ijk);
            RealVect normal;
            RealVect bcentv;
            bool dry_wall = false; // TODO: Implement this check
            bool wall_check = check_wall(
              p, ijk, ijk_prev, dx, plo, phi,
#ifdef AMREX_USE_EB
              eb_in_box, flags_array, bcent_fab, bnorm_fab,
#endif
              bndry_lo, bndry_hi, bloc, normal, bcentv);
            if (wall_check) {
              splash_type splash_flag = splash_type::no_impact;
              if (T_part > 0.) {
                SprayRefl SPRF;
                SPRF.pos_refl = p.pos();
                for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf)
                  SPRF.Y_refl[spf] = p.rdata(SPI.pstateY + spf);
                splash_flag = impose_wall(
                  p, SPI, SPU, *fdat, dx, plo, phi, wallT, bloc, normal, bcentv,
                  SPRF, isActive, dry_wall);
              } // TODO: Add check for if it is wall film
            }
          } // if (at_bounds...
        }   // End of p.id() > 0 check
      });   // End of loop over particles
  }         // for (int MyParIter pti..
}
