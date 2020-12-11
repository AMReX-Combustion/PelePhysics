
#include "SprayParticles.H"
#include <AMReX_Particles.H>
#include <AMReX_ParticleReduce.H>
#ifdef SPRAY_PELE_LM
#include "PeleLM.H"
#endif
#include "Transport.H"
#include "Drag.H"

using namespace amrex;

void
SprayParticleContainer::init_bcs()
{
  for (int dir = 0; dir < AMREX_SPACEDIM; dir++) {
    if (phys_bc->lo(dir) == Symmetry   ||
        phys_bc->lo(dir) == SlipWall   ||
        phys_bc->lo(dir) == NoSlipWall) {
      reflect_lo[dir] = true;
    } else {
      reflect_lo[dir] = false;
    }
    if (phys_bc->hi(dir) == Symmetry   ||
        phys_bc->hi(dir) == SlipWall   ||
        phys_bc->hi(dir) == NoSlipWall) {
      reflect_hi[dir] = true;
    } else {
      reflect_hi[dir] = false;
    }
  }
}

void
SprayParticleContainer::moveKick (MultiFab&   state,
                                  MultiFab&   source,
                                  const int   level,
                                  const Real& dt,
                                  const Real  time,
                                  const bool  isVirtualPart,
                                  const bool  isGhostPart,
                                  const int   state_ghosts,
                                  const int   source_ghosts,
                                  MultiFab*   u_mac)
{
  bool do_move = false;
  int width = 0;
  moveKickDrift(state, source, level, dt, time, isVirtualPart, isGhostPart,
                state_ghosts, source_ghosts, do_move, width, u_mac);
}

void
SprayParticleContainer::moveKickDrift (MultiFab&   state,
                                       MultiFab&   source,
                                       const int   level,
                                       const Real& dt,
                                       const Real  time,
                                       const bool  isVirtualPart,
                                       const bool  isGhostPart,
                                       const int   state_ghosts,
                                       const int   source_ghosts,
                                       const bool  do_move,
                                       const int   where_width,
                                       MultiFab*   u_mac)
{
  BL_PROFILE("ParticleContainer::moveKickDrift()");
  AMREX_ASSERT(u_mac == nullptr || u_mac[0].nGrow() >= 1);
  AMREX_ASSERT(level >= 0);
  AMREX_ASSERT(state.nGrow() >= 2);

  //If there are no particles at this level
  if (level >= this->GetParticles().size())
    return;

  const Real strttime = ParallelDescriptor::second();

  BL_PROFILE_VAR("SprayParticles::updateParticles()", UPD_PART);
  updateParticles(level, state, source, dt, time, state_ghosts, source_ghosts, do_move, u_mac);
  BL_PROFILE_VAR_STOP(UPD_PART);

  // Fill ghost cells after we've synced up ..
  // TODO: Check to see if this is needed at all
  // if (level > 0)
  //   source.FillBoundary(Geom(level).periodicity());

  // ********************************************************************************

  // ********************************************************************************

  if (this->m_verbose > 1) {
    Real stoptime = ParallelDescriptor::second() - strttime;
    ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());
    if (ParallelDescriptor::IOProcessor()) {
      if (do_move) {
        Print() << "SprayParticleContainer::moveKickDrift() time: "
                << stoptime << '\n';
      } else {
        Print() << "SprayParticleContainer::moveKick() time: "
                << stoptime << '\n';
      }
    }
  }
}

Real
SprayParticleContainer::estTimestep (int level, Real cfl) const
{
  BL_PROFILE("ParticleContainer::estTimestep()");
  AMREX_ASSERT(m_setFuelData);
  // TODO: Clean up this mess and bring the num particle functionality back
  Real dt = std::numeric_limits<Real>::max();
  if (level >= this->GetParticles().size() ||
      m_sprayIndx.mom_tran == 0)
    return -1.;

  const Real strttime = ParallelDescriptor::second();
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
      AMREX_D_TERM(const Real * up = attribs[0].data();,
                   const Real * vp = attribs[1].data();,
                   const Real * wp = attribs[2].data(););
#endif
      reduce_op.eval(n, reduce_data,
                     [=] AMREX_GPU_DEVICE (const Long i) -> ReduceTuple
      {
        const ParticleType& p = pstruct[i];
        // TODO: This assumes that pstateVel = 0 and dxi[0] = dxi[1] = dxi[2]
        if (p.id() > 0) {
#ifdef USE_SPRAY_SOA
          const Real max_mag_vdx = amrex::max(AMREX_D_DECL(std::abs(up[i]),
                                                           std::abs(vp[i]),
                                                           std::abs(wp[i])))*dxi[0];
#else
          const Real max_mag_vdx = amrex::max(AMREX_D_DECL(std::abs(p.rdata(0)),
                                                           std::abs(p.rdata(1)),
                                                           std::abs(p.rdata(2))))*dxi[0];
#endif
          Real dt_part = (max_mag_vdx > 0.) ? (cfl/max_mag_vdx) : 1.E50;
#ifdef SPRAY_PELE_LM
          // Conversion since particle velocities are in cm and dx is in m for PeleLM
          dt_part *= 100.;
#endif
          return dt_part;
        }
        return 1.E50;
      });
    }
    ReduceTuple hv = reduce_data.value();
    Real ldt_cpu = amrex::get<0>(hv);
    dt = amrex::min(dt,ldt_cpu);
  }
  ParallelDescriptor::ReduceRealMin(dt);
  // Check if the velocity of particles being injected
  // is greater existing particle velocities
  if (m_injectVel > 0.)
    dt = amrex::min(dt, cfl*dx[0]/m_injectVel);

  if (this->m_verbose > 1) {
    Real stoptime = ParallelDescriptor::second() - strttime;
    ParallelDescriptor::ReduceRealMax(stoptime,
                                      ParallelDescriptor::IOProcessorNumber());
    if (ParallelDescriptor::IOProcessor())
      std::cout << "SprayParticleContainer::estTimestep() time: "
                << stoptime << '\n';
  }
  return dt;
}

void
SprayParticleContainer::updateParticles(const int&  level,
                                        MultiFab&   state,
                                        MultiFab&   source,
                                        const Real& flow_dt,
                                        const Real& time,
                                        const int   state_ghosts,
                                        const int   source_ghosts,
                                        const bool  do_move,
                                        MultiFab*   u_mac)
{
  AMREX_ASSERT(m_setFuelData);
  AMREX_ASSERT(OnSameGrids(level, state));
  AMREX_ASSERT(OnSameGrids(level, source));
  const auto dxi = this->Geom(level).InvCellSizeArray();
  const auto dx = this->Geom(level).CellSizeArray();
  const auto plo = this->Geom(level).ProbLoArray();
  const auto phi = this->Geom(level).ProbHiArray();
  const auto domain = this->Geom(level).Domain();
  IntVect dom_lo = domain.smallEnd();
  IntVect dom_hi = domain.bigEnd();
  IntVect lo_bound;
  IntVect hi_bound;
  for (int dir = 0; dir != AMREX_SPACEDIM; ++dir) {
    if (!this->Geom(level).isPeriodic(dir)) {
      if (reflect_lo[dir]) lo_bound[dir] = 1;
      else lo_bound[dir] = -1;
      if (reflect_hi[dir]) hi_bound[dir] = 1;
      else hi_bound[dir] = -1;
    } else {
      lo_bound[dir] = 0;
      hi_bound[dir] = 0;
    }
    // Only concerned with reflective boundaries
    if (!reflect_lo[dir]) dom_lo[dir] -= 100;
    if (!reflect_hi[dir]) dom_hi[dir] += 100;
  }
  const Real vol = AMREX_D_TERM(dx[0],*dx[1],*dx[2]);
  const Real inv_vol = 1./vol;
  // Set all constants
  const Real num_ppp = m_parcelSize;
  Real mw_fluid[NUM_SPECIES];
  Real invmw[NUM_SPECIES];
  EOS::molecular_weight(mw_fluid);
  EOS::inv_molecular_weight(invmw);
  const Real ref_T = m_sprayRefT;
  // Particle components indices
  SprayComps SPI = m_sprayIndx;
  Real vel_conv = 1.;
  Real pos_conv = 1.;
  Real rho_conv = 1.;
  Real eng_conv = 1.;
  Real mom_src_conv = 1.;
  Real mass_src_conv = 1.;
  Real eng_src_conv = 1.;
#ifdef SPRAY_PELE_LM
  vel_conv = 100.;  // Turn m/s to cm/s
  rho_conv = 0.001; // Turn kg/m^3 to g/cm^3
  pos_conv = 0.01;  // Turn cm to m for updating position
  eng_conv = 1.E4; // For converting enthalpy to CGS
  // This makes no sense, conversions should be independent
  // of dimensions but numerical tests show this isn't the case
#if AMREX_SPACEDIM == 2
  mom_src_conv = 1.E-3;
  mass_src_conv = 1.E-1;
  eng_src_conv = 1.E-5;
#elif AMREX_SPACEDIM == 3
  mom_src_conv = 1.E-5;
  mass_src_conv = 1.E-3;
  eng_src_conv = 1.E-7;
#endif
#endif

  // Start the ParIter, which loops over separate sets of particles in different boxes
  for (MyParIter pti(*this, level); pti.isValid(); ++pti) {
    const Box& tile_box = pti.tilebox();
    const Box& state_box = pti.growntilebox(state_ghosts);
    const Box& src_box = pti.growntilebox(source_ghosts);
    const Long Np = pti.numParticles();
    ParticleType* pstruct = &(pti.GetArrayOfStructs()[0]);
    // Get particle attributes if StructOfArrays are used
#ifdef USE_SPRAY_SOA
    auto& attribs = pti.GetAttribs();
    Real * velp[AMREX_SPACEDIM];
    AMREX_D_TERM(velp[0] = attribs[SPI.pstateVel].data();,
                 velp[1] = attribs[SPI.pstateVel+1].data();,
                 velp[2] = attribs[SPI.pstateVel+2].data(););
    Real * Tp = attribs[SPI.pstateT].data();
    Real * diap = attribs[SPI.pstateDia].data();
    Real * rhop = attribs[SPI.pstateRho].data();
    std::array<Real *, SPRAY_FUEL_NUM> Yp;
    for (int spf = 0; spf != SPRAY_FUEL_NUM; ++spf) {
      Yp[spf] = attribs[SPI.pstateY+spf].data();
    }
#endif
    SprayData fdat = m_fuelData.getSprayData();
    Array4<const Real> const& statearr = state.array(pti);
    Array4<Real> const& sourcearr = source.array(pti);
// #ifdef SPRAY_PELE_LM
//     GpuArray<
//       Array4<const Real>, AMREX_SPACEDIM> const
//       umac{AMREX_D_DECL(u_mac[0].array(pti), u_mac[1].array(pti), u_mac[2].array(pti))};
// #endif
    AMREX_FOR_1D ( Np, i,
    {
      ParticleType& p = pstruct[i];
      if (p.id() > 0) {
        // TODO: I was hoping to not have to instantiate this everytime
        Real Y_fluid[NUM_SPECIES];
        Real mass_frac[NUM_SPECIES];
        // Weights for interpolation
        Real coef[AMREX_D_PICK(2, 4, 8)];
        // Indices of adjacent cells
        IntVect indx_array[AMREX_D_PICK(2, 4, 8)];
        // Cell-based length of particle in domain
        RealVect len(AMREX_D_DECL((p.pos(0) - plo[0])*dxi[0] + 0.5,
                                  (p.pos(1) - plo[1])*dxi[1] + 0.5,
                                  (p.pos(2) - plo[2])*dxi[2] + 0.5));
        RealVect vel_fluid(RealVect::TheZeroVector());
// #ifdef SPRAY_PELE_LM
//         InterpolateFaceVelocity(len, dom_lo, dom_hi, umac, vel_fluid);
// #endif
        Real T_fluid = 0.;
        Real rho_fluid = 0.;
        for (int sp = 0; sp != NUM_SPECIES; ++sp) Y_fluid[sp] = 0.;
        // Do initial interpolation
        AdjIndexWeights(len, indx_array, coef, dom_lo, dom_hi);
        // Extract adjacent values and interpolate fluid at particle location
        for (int aindx = 0; aindx != AMREX_D_PICK(2, 4, 8); ++aindx) {
          IntVect cur_indx = indx_array[aindx];
#ifdef AMREX_DEBUG
          if (!state_box.contains(cur_indx))
            Abort("SprayParticleContainer::updateParticles() -- state box too small");
#endif
          Real cur_coef = coef[aindx];
          Real cur_rho = statearr(cur_indx, SPI.rhoIndx);
          rho_fluid += cur_coef*cur_rho*rho_conv;
          Real inv_rho = 1./cur_rho;
          for (int sp = 0; sp != NUM_SPECIES; ++sp) {
            int mf_indx = sp + SPI.specIndx;
            Real cur_mf = statearr(cur_indx, mf_indx)*inv_rho;
            Y_fluid[sp] += cur_coef*cur_mf;
            mass_frac[sp] = cur_mf;
          }
#ifdef SPRAY_PELE_LM
          inv_rho = 1.; // Since velocity is provided instead of momentum
#endif
          AMREX_D_TERM(Real velx = statearr(cur_indx, SPI.momIndx)*inv_rho*vel_conv;
                       Real ke = 0.5*velx*velx;
                       vel_fluid[0] += cur_coef*velx;,
                       Real vely = statearr(cur_indx, SPI.momIndx+1)*inv_rho*vel_conv;
                       ke += 0.5*vely*vely;
                       vel_fluid[1] += cur_coef*vely;,
                       Real velz = statearr(cur_indx, SPI.momIndx+2)*inv_rho*vel_conv;
                       ke += 0.5*velz*velz;
                       vel_fluid[2] += cur_coef*velz;);
          Real T_val = statearr(cur_indx, SPI.utempIndx);
#ifndef SPRAY_PELE_LM
          Real intEng = statearr(cur_indx, SPI.engIndx)*inv_rho - ke;
          EOS::EY2T(intEng, mass_frac, T_val);
#endif
          T_fluid += cur_coef*T_val;
        }
        GasPhaseVals gpv(vel_fluid, T_fluid, rho_fluid, Y_fluid,
                         mw_fluid, invmw, ref_T);
        bool remove_particle = calculateSpraySource(flow_dt, do_move, gpv, SPI, fdat, p
#ifdef USE_SPRAY_SOA
                                                    , attribs, i
#endif
                                                    );
        // Modify particle position by whole time step
        if (do_move) {
          for (int dir = 0; dir != AMREX_SPACEDIM; ++dir) {
#ifdef USE_SPRAY_SOA
            Real& cvel = velp[dir][i];
#else
            Real& cvel = p.rdata(SPI.pstateVel+dir);
#endif
            Gpu::Atomic::Add(&p.pos(dir), flow_dt*cvel*pos_conv);
            remove_particle = checkWall(p.pos(dir), cvel,
                                        phi[dir], plo[dir], hi_bound[dir], lo_bound[dir]);
          }
        }
        if (remove_particle) p.id() = -1;
        for (int aindx = 0; aindx != AMREX_D_PICK(2, 4, 8); ++aindx) {
          Real cur_coef = -coef[aindx]*inv_vol*num_ppp;
          IntVect cur_indx = indx_array[aindx];
#ifdef AMREX_DEBUG
          if (!src_box.contains(cur_indx))
            Abort("SprayParticleContainer::updateParticles() -- source box too small");
#endif
          if (SPI.mom_tran) {
            for (int dir = 0; dir != AMREX_SPACEDIM; ++dir) {
              const int nf = SPI.momIndx + dir;
              Gpu::Atomic::Add(&sourcearr(cur_indx, nf), cur_coef*gpv.fluid_mom_src[dir]*mom_src_conv);
            }
          }
          if (SPI.mass_tran) {
            Gpu::Atomic::Add(&sourcearr(cur_indx, SPI.rhoIndx), cur_coef*gpv.fluid_mass_src*mass_src_conv);
            for (int spf = 0; spf != SPRAY_FUEL_NUM; ++spf) {
              const int nf = SPI.specIndx + fdat.indx(spf);
              Gpu::Atomic::Add(&sourcearr(cur_indx, nf), cur_coef*gpv.fluid_Y_dot[spf]*mass_src_conv);
            }
          }
          Gpu::Atomic::Add(&sourcearr(cur_indx, SPI.engIndx), cur_coef*gpv.fluid_eng_src*eng_src_conv);
        }
      } // End of p.id() > 0 check
    }); // End of loop over particles
  }
}
