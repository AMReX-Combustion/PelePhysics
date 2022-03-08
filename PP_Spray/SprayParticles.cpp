
#include "SprayParticles.H"
#include <AMReX_ParticleReduce.H>
#include <AMReX_Particles.H>
#include "Drag.H"
#include "SprayInterpolation.H"
#include "Transport.H"
#include "WallFunctions.H"
#ifdef AMREX_USE_EB
#include <AMReX_EBFArrayBox.H>
#endif

using namespace amrex;

void
getPSatCoef(
  Real* psat_coef, ParmParse& ppp, std::string fuel_name, const int spf)
{
  std::string psat_read = fuel_name + "_psat";
  std::vector<Real> inp_coef(4, 0.);
  ppp.queryarr(psat_read.c_str(), inp_coef);
  for (int i = 0; i < 4; ++i) {
    psat_coef[4 * spf + i] = inp_coef[i];
  }
}

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
SprayParticleContainer::readSprayParams(
  int& particle_verbose,
  Real& particle_cfl,
  Real& wall_temp,
  int& mass_trans,
  int& mom_trans,
  int& write_spray_ascii_files,
  int& plot_spray_src,
  int& init_function,
  std::string& init_file,
  SprayData& sprayData,
  std::string* sprayFuelNames)
//  Vector<std::string>& sprayFuelNames)
{
  amrex::ParmParse pp("particles");
  //
  // Control the verbosity of the Particle class
  pp.query("v", particle_verbose);

  pp.get("mass_transfer", mass_trans);
  pp.get("mom_transfer", mom_trans);
  pp.query("cfl", particle_cfl);
  if (particle_cfl > 0.5) {
    amrex::Abort("particles.cfl must be <= 0.5");
  }
  // Number of fuel species in spray droplets
  // Must match the number specified at compile time
  const int nfuel = pp.countval("fuel_species");
  if (nfuel != SPRAY_FUEL_NUM) {
    amrex::Abort(
      "Number of fuel species in input file must match SPRAY_FUEL_NUM");
  }

  std::vector<std::string> fuel_names;
  std::vector<Real> crit_T;
  std::vector<Real> boil_T;
  std::vector<Real> spraycp;
  std::vector<Real> latent;
  std::vector<Real> sprayrho;
  std::vector<Real> mu(nfuel, 0.);
  std::vector<Real> lambda(nfuel, 0.);
  {
    //sprayFuelNames.assign(nfuel, "");
    pp.getarr("fuel_species", fuel_names);
    pp.getarr("fuel_crit_temp", crit_T);
    pp.getarr("fuel_boil_temp", boil_T);
    pp.getarr("fuel_cp", spraycp);
    pp.getarr("fuel_latent", latent);
    pp.getarr("fuel_rho", sprayrho);
    pp.queryarr("fuel_mu", mu);
    pp.queryarr("fuel_lambda", lambda);
    for (int i = 0; i < nfuel; ++i) {
      sprayFuelNames[i] = fuel_names[i];
      sprayData.critT[i] = crit_T[i];
      sprayData.boilT[i] = boil_T[i];
      sprayData.cp[i] = spraycp[i];
      sprayData.latent[i] = latent[i];
      sprayData.ref_latent[i] = latent[i];
      sprayData.rho[i] = sprayrho[i];
      sprayData.mu[i] = mu[i];
      sprayData.lambda[i] = lambda[i];
      getPSatCoef(sprayData.psat_coef.data(), pp, fuel_names[i], i);
    }
  }

  Real parcel_size = 1;
  Real spray_ref_T = 300.;
  Real spray_sigma = -1.;
  bool splash_model = false;
  //
  // Set the number of particles per parcel
  //
  pp.query("parcel_size", parcel_size);
  pp.query("use_splash_model", splash_model);
  if (splash_model) {
    if (!pp.contains("fuel_sigma") || !pp.contains("wall_temp")) {
      Print() << "fuel_sigma and wall_temp must be set for splash model. "
              << "Set use_splash_model = false to turn off splash model"
              << std::endl;
      Abort();
    }
    // Set the fuel surface tension and contact angle
    pp.get("fuel_sigma", spray_sigma);
    // TODO: Have this retrieved from proper boundary data
    pp.get("wall_temp", wall_temp);
  }

  // Must use same reference temperature for all fuels
  // TODO: This means the reference temperature must be the same for all fuel
  // species
  pp.get("fuel_ref_temp", spray_ref_T);
  //
  // Set if spray ascii files should be written
  //
  pp.query("write_spray_ascii_files", write_spray_ascii_files);
  //
  // Set if gas phase spray source term should be written
  //
  pp.query("plot_src", plot_spray_src);
  //
  // Used in initData() on startup to read in a file of particles.
  //
  pp.query("init_file", init_file);
  //
  // Used in initData() on startup to set the particle field using the
  // SprayParticlesInitInsert.cpp problem specific function
  //
  pp.query("init_function", init_function);

  sprayData.num_ppp = parcel_size;
  sprayData.ref_T = spray_ref_T;
  sprayData.sigma = spray_sigma;

  if (particle_verbose && ParallelDescriptor::IOProcessor()) {
    amrex::Print() << "Spray fuel species " << sprayFuelNames[0];
    for (int i = 1; i < SPRAY_FUEL_NUM; ++i) {
      amrex::Print() << ", " << sprayFuelNames[i];
    }
    amrex::Print() << std::endl;
    amrex::Print() << "Number of particles per parcel " << parcel_size
                   << std::endl;
  }
  //
  // Force other processors to wait till directory is built.
  //
  ParallelDescriptor::Barrier();
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
  MultiFab* u_mac)
{
  bool do_move = false;
  moveKickDrift(
    state, source, level, dt, time, isVirtualPart, isGhostPart, state_ghosts,
    source_ghosts, do_move, ltransparm, u_mac);
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
  MultiFab* u_mac)
{
  BL_PROFILE("ParticleContainer::moveKickDrift()");
  AMREX_ASSERT(u_mac == nullptr || u_mac[0].nGrow() >= 1);
  AMREX_ASSERT(level >= 0);
  AMREX_ASSERT(state.nGrow() >= 2);

  // If there are no particles at this level
  if (level >= this->GetParticles().size()) {
    return;
  }

  bool isActive = !(isVirtualPart || isGhostPart);

  updateParticles(
    level, state, source, dt, time, state_ghosts, source_ghosts, isActive,
    do_move, ltransparm, u_mac);

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
  if (level >= this->GetParticles().size() || m_sprayIndx.mom_tran == 0) {
    return -1.;
  }

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
  const bool isActive,
  const bool do_move,
  pele::physics::transport::TransParm<
  pele::physics::EosType,
  pele::physics::TransportType> const* ltransparm,
  MultiFab* /*u_mac*/)
{
  BL_PROFILE("SprayParticleContainer::updateParticles()");
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
  const auto* cellcent = &(factory.getCentroid());
  const auto* bndrycent = &(factory.getBndryCent());
  const auto* bndryarea = &(factory.getBndryArea());
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
  const Real wallT = m_wallT;
  const Real vol = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);
  const Real inv_vol = 1. / vol;
  // Particle components indices
  SprayComps SPI = m_sprayIndx;
  // Start the ParIter, which loops over separate sets of particles in different
  // boxes
  for (MyParIter pti(*this, level); pti.isValid(); ++pti) {
    const Box tile_box = pti.tilebox();
#ifdef AMREX_DEBUG
    const Box src_box = pti.growntilebox(source_ghosts);
#endif
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
    const Box state_box = pti.growntilebox(state_ghosts);
    bool eb_in_box = true;
    const auto& interp_fab = static_cast<EBFArrayBox const&>(state[pti]);
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
      Np, [pstruct, statearr, sourcearr, plo, phi, dx, dxi, do_move, SPI, fdat,
           bndry_hi, bndry_lo, flow_dt, inv_vol, ltransparm, at_bounds, wallT,
           isActive
#ifdef AMREX_DEBUG
           ,
           src_box
#endif
#ifdef USE_SPRAY_SOA
           ,
           attribs
#endif
#ifdef AMREX_USE_EB
           ,
           flags_array, ccent_fab, bcent_fab, bnorm_fab, barea_fab,
           volfrac_fab, eb_in_box
#endif
    ] AMREX_GPU_DEVICE(int pid) noexcept {
        ParticleType& p = pstruct[pid];
        if (p.id() > 0) {
          auto eos = pele::physics::PhysicsType::eos();
          SprayUnits SPU;
          GasPhaseVals gpv;
          eos.molecular_weight(gpv.mw_fluid.data());
          eos.inv_molecular_weight(gpv.invmw.data());
          for (int n = 0; n < NUM_SPECIES; ++n) {
            gpv.mw_fluid[n] *= SPU.mass_conv;
            gpv.invmw[n] /= SPU.mass_conv;
          }
          GpuArray<IntVect, AMREX_D_PICK(2, 4, 8)>
            indx_array; // Array of adjacent cells
          GpuArray<Real, AMREX_D_PICK(2, 4, 8)>
            weights; // Array of corresponding weights
          RealVect lx = (p.pos() - plo) * dxi + 0.5;
          IntVect ijk = lx.floor(); // Upper cell center
          RealVect lxc = (p.pos() - plo) * dxi;
          IntVect ijkc = lxc.floor(); // Cell with particle
          bool is_wall_film = false;
          Real face_area = 0.;
          IntVect bflags(IntVect::TheZeroVector());
          // If the temperature is a negative value, the particle is a wall film
          if (p.rdata(SPI.pstateT) < 0.) {
            is_wall_film = true;
          }
          if (at_bounds) {
            // Check if particle has left the domain, is wall film,
            // or is boundary adjacent and must be shifted
            bool left_dom = check_bounds(
              p.pos(), plo, phi, dx, bndry_lo, bndry_hi, ijk, bflags);
            if (left_dom) {
              Abort("Particle has incorrectly left domain");
            }
            if (is_wall_film) {
              // TODO: Assumes grid spacing is uniform in all directions
              face_area = AMREX_D_TERM(1., *dx[0], *dx[0]);
            }
          }
          // Length from cell center to boundary face center
          Real diff_cent = 0.5 * dx[0];
          bool do_fe_interp = false;
#ifdef AMREX_USE_EB
          do_fe_interp = true;
          int i = 0, j = 0, k = 0, ip = 0, jp = 0, kp = 0;
          // Cell containing particle centroid
          AMREX_D_TERM(ip = ijkc[0];, jp = ijkc[1];
                       , kp = ijkc[2];);
          AMREX_D_TERM(i = ijk[0];, j = ijk[1];
                       , k = ijk[2];);
          if (!eb_in_box) {
            do_fe_interp = false;
          } else {
            // All cells in the stencil are regular. Use
            // traditional trilinear interpolation
            if (
#if AMREX_SPACEDIM == 3
              flags_array(i - 1, j - 1, k - 1).isRegular() &&
              flags_array(i, j - 1, k - 1).isRegular() &&
              flags_array(i - 1, j, k - 1).isRegular() &&
              flags_array(i, j, k - 1).isRegular() &&
              flags_array(i - 1, j - 1, k).isRegular() &&
#endif
              flags_array(i, j - 1, k).isRegular() &&
              flags_array(i - 1, j, k).isRegular() &&
              flags_array(i, j, k).isRegular()) {
              do_fe_interp = false;
            }
          }
          if (do_fe_interp) {
            fe_interp(
              p.pos(), ip, jp, kp, dx, dxi, plo, flags_array, ccent_fab,
              bcent_fab, bnorm_fab, volfrac_fab, indx_array.data(),
              weights.data());
          } else {
            trilinear_interp(
              ijk, lx, indx_array.data(), weights.data(), bflags);
          }
          bool eb_wall_film = false;
          if (is_wall_film && flags_array(ip, jp, kp).isSingleValued()) {
            eb_wall_film = true;
            face_area = barea_fab(ip, jp, kp);
            diff_cent = 0.;
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
              diff_cent += std::pow(
                bcent_fab(ip, jp, kp, dir) - ccent_fab(ip, jp, kp, dir), 2);
            }
            diff_cent = std::sqrt(diff_cent);
          }
#else
          trilinear_interp(ijk, lx, indx_array.data(), weights.data(), bflags);
#endif // AMREX_USE_EB
       // Interpolate fluid state
          {
            GpuArray<Real, NUM_SPECIES> mass_frac;
            for (int aindx = 0.; aindx < AMREX_D_PICK(2, 4, 8); ++aindx) {
              IntVect cur_indx = indx_array[aindx];
              Real cw = weights[aindx];
              Real cur_rho = statearr(cur_indx, SPI.rhoIndx);
              gpv.rho_fluid += cw * cur_rho;
              Real inv_rho = 1. / cur_rho;
              for (int n = 0; n < NUM_SPECIES; ++n) {
                int sp = n + SPI.specIndx;
                Real cur_mf = statearr(cur_indx, sp) * inv_rho;
                gpv.Y_fluid[n] += cw * cur_mf;
                mass_frac[n] = cur_mf;
              }
#ifdef SPRAY_PELE_LM
              inv_rho = 1.;
#endif
              Real ke = 0.;
              for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                Real vel = statearr(cur_indx, SPI.momIndx + dir) * inv_rho;
                gpv.vel_fluid[dir] += cw * vel;
                ke += vel * vel / 2.;
              }
              Real T_i = statearr(cur_indx, SPI.utempIndx);
#ifndef SPRAY_PELE_LM
              Real intEng = statearr(cur_indx, SPI.engIndx) * inv_rho - ke;
              eos.EY2T(intEng, mass_frac.data(), T_i);
#endif
              gpv.T_fluid += cw * T_i;
            }
          }
          // Solve for avg mw and pressure at droplet location
          gpv.define();
          if (!is_wall_film) {
            calculateSpraySource(
              flow_dt, gpv, SPI, *fdat, p,
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
                const Real cvel = p.rdata(SPI.pstateVel + dir);
#endif
                Gpu::Atomic::Add(&p.pos(dir), flow_dt * cvel);
              }
            }
          } else {
            calculateWallFilmSource(
              flow_dt, gpv, SPI, *fdat, p,
#ifdef USE_SPRAY_SOA
              attribs, pid,
#endif
              wallT, face_area, diff_cent, ltransparm);
          }
          for (int aindx = 0; aindx < AMREX_D_PICK(2, 4, 8); ++aindx) {
            IntVect cur_indx = indx_array[aindx];
            Real cvol = inv_vol;
#ifdef AMREX_USE_EB
            if (!flags_array(cur_indx).isRegular()) {
              cvol *= 1. / (volfrac_fab(cur_indx));
            }
#endif
            Real cur_coef = -weights[aindx] * fdat->num_ppp * cvol;
#ifdef AMREX_DEBUG
            if (!src_box.contains(cur_indx)) {
              Abort(
                "SprayParticleContainer::updateParticles() -- source box too "
                "small");
            }
#endif
            if (SPI.mom_tran) {
              for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                const int nf = SPI.momSrcIndx + dir;
                Gpu::Atomic::Add(
                  &sourcearr(cur_indx, nf), cur_coef * gpv.fluid_mom_src[dir]);
              }
            }
            if (SPI.mass_tran) {
              Gpu::Atomic::Add(
                &sourcearr(cur_indx, SPI.rhoSrcIndx),
                cur_coef * gpv.fluid_mass_src);
              for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
                const int nf = SPI.specSrcIndx + spf;
                Gpu::Atomic::Add(
                  &sourcearr(cur_indx, nf), cur_coef * gpv.fluid_Y_dot[spf]);
              }
            }
            Gpu::Atomic::Add(
              &sourcearr(cur_indx, SPI.engSrcIndx),
              cur_coef * gpv.fluid_eng_src);
          }
          // Solve for splash model/wall film formation
          if ((at_bounds || do_fe_interp) && do_move) {
#ifdef AMREX_USE_EB
            IntVect ijkc_prev = ijkc;
#endif
            lx = (p.pos() - plo) * dxi + 0.5;
            ijk = lx.floor();
            lxc = (p.pos() - plo) * dxi;
            ijkc = lxc.floor(); // New cell center
            const Real T_part = p.rdata(SPI.pstateT);
            IntVect bloc(ijkc);
            RealVect normal;
            RealVect bcentv;
            bool dry_wall = false; // TODO: Implement this check
            // First check if particle has exited the domain through a Cartesian
            // boundary
            bool left_dom = check_bounds(
              p.pos(), plo, phi, dx, bndry_lo, bndry_hi, ijk, bflags);
            if (left_dom) {
              p.id() = -1;
            } else {
              bool wall_check = check_wall(
                bflags,
#ifdef AMREX_USE_EB
                ijkc, ijkc_prev, eb_in_box, flags_array, bcent_fab, bnorm_fab,
                bloc,
#endif
                normal, bcentv);
              if (wall_check) {
                splash_type splash_flag = splash_type::no_impact;
                if (T_part > 0.) {
                  SprayRefl SPRF;
                  SPRF.pos_refl = p.pos();
                  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
                    SPRF.Y_refl[spf] = p.rdata(SPI.pstateY + spf);
                  }
                  splash_flag = impose_wall(
                    p, SPI, *fdat, dx, plo, wallT, bloc, normal, bcentv, SPRF,
                    isActive, dry_wall);
                }
              } // if (wall_check)
            }   // if (left_dom)
          }     // if (at_bounds...
        }       // End of p.id() > 0 check
      });       // End of loop over particles
  }             // for (int MyParIter pti..
}
