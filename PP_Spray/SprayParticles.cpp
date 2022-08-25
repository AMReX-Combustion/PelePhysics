
#include "SprayParticles.H"
#include "Drag.H"
#include "SprayInterpolation.H"
#include "Transport.H"
#include "WallFunctions.H"
#ifdef AMREX_USE_EB
#include <AMReX_EBFArrayBox.H>
#endif

using namespace amrex;

void
getInpCoef(
  Real* coef,
  const ParmParse& ppp,
  const std::string& fuel_name,
  const std::string& varname,
  const int spf)
{
  std::string psat_read = fuel_name + "_" + varname;
  std::vector<Real> inp_coef(4, 0.);
  ppp.queryarr(psat_read.c_str(), inp_coef);
  for (int i = 0; i < 4; ++i) {
    coef[4 * spf + i] = inp_coef[i];
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
  int& write_spray_ascii_files,
  int& plot_spray_src,
  int& init_function,
  std::string& init_file,
  SprayData& sprayData,
  std::string* sprayFuelNames,
  Vector<std::string>& derivePlotVars,
  const Real& max_cfl)
{
  amrex::ParmParse pp("particles");
  //
  // Control the verbosity of the Particle class
  pp.query("v", particle_verbose);

  pp.query("mass_transfer", sprayData.mass_trans);
  pp.query("mom_transfer", sprayData.mom_trans);
  pp.query("fixed_parts", sprayData.fixed_parts);
  pp.query("cfl", particle_cfl);
  if (particle_cfl > max_cfl) {
    std::string errorstr =
      "particles.cfl must be <= " + std::to_string(max_cfl);
    Abort(errorstr);
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
      getInpCoef(sprayData.psat_coef.data(), pp, fuel_names[i], "psat", i);
      getInpCoef(sprayData.rho_coef.data(), pp, fuel_names[i], "rho", i);
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
    Abort("Splash model is not fully implemented");
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
#ifdef AMREX_USE_EB
  //
  // Spray source terms are only added to cells with a volume fraction higher
  // than this value
  //
  pp.query("min_eb_vfrac", sprayData.min_eb_vfrac);
#endif

  sprayData.num_ppp = parcel_size;
  sprayData.ref_T = spray_ref_T;
  sprayData.sigma = spray_sigma;

  // List of known derived spray quantities
  std::vector<std::string> derive_names = {
    "spray_mass",      // Total liquid mass in a cell
    "spray_num",       // Number of spray droplets in a cell
    "spray_vol",       // Total liquid volume in a cell
    "spray_surf_area", // Total liquid surface area in a cell
    "spray_vol_frac",  // Volume fraction of liquid in cell
    "d10",             // Average diameter
    "d32",             // SMD
    "spray_temp",      // Mass-weighted average temperature
    AMREX_D_DECL("spray_x_vel", "spray_y_vel", "spray_z_vel")};
  int derive_plot_vars = 1;
  pp.query("derive_plot_vars", derive_plot_vars);
  int derive_plot_species = 1;
  pp.query("derive_plot_species", derive_plot_species);
  // If derive_spray_vars if present, add above spray quantities in the same
  // order
  for (const auto& derive_name : derive_names) {
    derivePlotVars.push_back(derive_name);
  }
  if (derive_plot_species == 1 && SPRAY_FUEL_NUM > 1) {
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      derivePlotVars.push_back("spray_mass_" + sprayFuelNames[spf]);
    }
  }

  if (particle_verbose >= 1 && ParallelDescriptor::IOProcessor()) {
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
  const Real spray_cfl_lev,
  MultiFab* u_mac)
{
  bool do_move = false;
  moveKickDrift(
    state, source, level, dt, time, isVirtualPart, isGhostPart, state_ghosts,
    source_ghosts, do_move, ltransparm, spray_cfl_lev, u_mac);
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
  const Real spray_cfl_lev,
  MultiFab* u_mac)
{
  BL_PROFILE("ParticleContainer::moveKickDrift()");
  AMREX_ASSERT(u_mac == nullptr || u_mac[0].nGrow() >= 1);
  AMREX_ASSERT(level >= 0);

  // If there are no particles at this level
  if (level >= this->GetParticles().size()) {
    return;
  }

  updateParticles(
    level, state, source, dt, time, state_ghosts, source_ghosts, isVirtualPart,
    isGhostPart, do_move, ltransparm, spray_cfl_lev, u_mac);

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
  if (level >= this->GetParticles().size() || m_sprayData->fixed_parts) {
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
      reduce_op.eval(
        n, reduce_data, [=] AMREX_GPU_DEVICE(const Long i) -> ReduceTuple {
          const ParticleType& p = pstruct[i];
          // TODO: This assumes that pstateVel = 0 and dxi[0] = dxi[1] =
          // dxi[2]
          if (p.id() > 0) {
            const Real max_mag_vdx =
              amrex::max(AMREX_D_DECL(
                amrex::Math::abs(p.rdata(0)), amrex::Math::abs(p.rdata(1)),
                amrex::Math::abs(p.rdata(2)))) *
              dxi[0];
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
  const bool isVirt,
  const bool isGhost,
  const bool do_move,
  pele::physics::transport::TransParm<
    pele::physics::EosType,
    pele::physics::TransportType> const* ltransparm,
  const Real spray_cfl_lev,
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
  // If particle subcycling is being done, determine the number of subcycles
  // Note: this is different than the AMR subcycling
  Real sub_cfl = 0.5; // CFL for each subcycle
  Real sub_dt = flow_dt;
  int num_iter = 1;
  if (do_move && spray_cfl_lev > sub_cfl) {
    num_iter = int(std::ceil(spray_cfl_lev / sub_cfl));
    sub_dt = flow_dt / Real(num_iter);
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
    ParticleType* pstruct = &(pti.GetArrayOfStructs()[0]);
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

#ifdef AMREX_USE_EB
    bool eb_in_box = true;
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
    // #ifdef PELELM_USE_SPRAY
    //     GpuArray<
    //       Array4<const Real>, AMREX_SPACEDIM> const
    //       umac{AMREX_D_DECL(u_mac[0].array(pti), u_mac[1].array(pti),
    //       u_mac[2].array(pti))};
    // #endif
    amrex::ParallelFor(
      Np,
      [pstruct, Tarr, rhoYarr, rhoarr, momarr, engarr, rhoYSrcarr, rhoSrcarr,
       momSrcarr, engSrcarr, plo, phi, dx, dxi, do_move, SPI, fdat, bndry_hi,
       bndry_lo, flow_dt, inv_vol, ltransparm, at_bounds, wallT, isGhost,
       isVirt, src_box, state_box, sub_cfl, num_iter, sub_dt, spray_cfl_lev
#ifdef AMREX_USE_EB
       ,
       flags_array, ccent_fab, bcent_fab, bnorm_fab, volfrac_fab, eb_in_box
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
            // and must be shifted
            bool left_dom = check_bounds(
              p.pos(), plo, phi, dx, bndry_lo, bndry_hi, ijk, bflags);
            if (left_dom) {
              Abort("Particle has incorrectly left the domain");
            }
          }
          // Subcycle loop
          Real ctime = 0.;
          Real cur_dt = sub_dt;
          Real cur_cfl = spray_cfl_lev;
          int cur_iter = 0;
          while (p.id() > 0 && cur_iter < num_iter) {
            cur_cfl -= sub_cfl;
            // Flag for whether we are near EB boundaries
            bool do_fe_interp = false;
#ifdef AMREX_USE_EB
            if (eb_in_box) {
              do_fe_interp = eb_interp(
                p, SPI, isVirt, ijkc, ijk, dx, dxi, lx, plo, bflags,
                flags_array, ccent_fab, bcent_fab, bnorm_fab, volfrac_fab,
                fdat->min_eb_vfrac, indx_array.data(), weights.data());
            } else
#endif
            {
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
            calculateSpraySource(cur_dt, gpv, SPI, *fdat, p, ltransparm);
            for (int aindx = 0; aindx < AMREX_D_PICK(2, 4, 8); ++aindx) {
              IntVect cur_indx = indx_array[aindx];
              Real cvol = inv_vol;
#ifdef AMREX_USE_EB
              if (flags_array(cur_indx).isSingleValued()) {
                cvol *= 1. / (volfrac_fab(cur_indx));
              }
#endif
              Real cur_coef =
                -weights[aindx] * fdat->num_ppp * cvol * cur_dt / flow_dt;
              if (!src_box.contains(cur_indx)) {
                if (!isGhost) {
                  Abort("SprayParticleContainer::updateParticles() -- source "
                        "box too small");
                } else {
                  continue;
                }
              }
              if (fdat->mom_trans) {
                for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                  Gpu::Atomic::Add(
                    &momSrcarr(cur_indx, dir),
                    cur_coef * gpv.fluid_mom_src[dir]);
                }
              }
              if (fdat->mass_trans) {
                Gpu::Atomic::Add(
                  &rhoSrcarr(cur_indx), cur_coef * gpv.fluid_mass_src);
                for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
                  Gpu::Atomic::Add(
                    &rhoYSrcarr(cur_indx, spf),
                    cur_coef * gpv.fluid_Y_dot[spf]);
                }
              }
              Gpu::Atomic::Add(
                &engSrcarr(cur_indx), cur_coef * gpv.fluid_eng_src);
            }
            // Modify particle position by whole time step
            if (do_move && !fdat->fixed_parts) {
              for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                const Real cvel = p.rdata(SPI.pstateVel + dir);
                p.pos(dir) += cur_dt * cvel;
              }
              // Update indices
              ijkc_prev = ijkc;
              lx = (p.pos() - plo) * dxi + 0.5;
              ijk = lx.floor();
              lxc = (p.pos() - plo) * dxi;
              ijkc = lxc.floor(); // New cell center
              if (at_bounds || do_fe_interp) {
                // First check if particle has exited the domain through a
                // Cartesian boundary
                bool left_dom = check_bounds(
                  p.pos(), plo, phi, dx, bndry_lo, bndry_hi, ijk, bflags);
                if (left_dom) {
                  p.id() = -1;
                } else {
                  // Next reflect particles off BC or EB walls if necessary
                  impose_wall(
                    p, SPI, dx, plo, bflags,
#ifdef AMREX_USE_EB
                    eb_in_box, flags_array, bcent_fab, bnorm_fab, vfrac_fab,
                    fdat->min_eb_vfrac,
#endif
                    ijkc, ijkc_prev);
                  lx = (p.pos() - plo) * dxi + 0.5;
                  ijk = lx.floor();
                  lxc = (p.pos() - plo) * dxi;
                  ijkc = lxc.floor();
                }
              } // if (at_bounds || fe_interp)
            }   // if (do_move)
            cur_iter++;
            ctime += cur_dt;
            if (isGhost && !src_box.contains(ijkc)) {
              p.id() = -1;
            }
          } // End of subcycle loop
        }   // End of p.id() > 0 check
      });   // End of loop over particles
  }         // for (int MyParIter pti..
}
