
#include <SprayParticles.H>
#include <AMReX_Particles.H>
#include <AMReX_ParticleReduce.H>
#include "Constants.H"
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

// void
// SprayParticleContainer::SetAll (Real val, int pstate_idx, int lev)
// {
//   BL_ASSERT(lev >= 0 && lev < GetParticles().size());
//   ParticleLevel& plev = GetParticles(lev);
//   for (auto& kv : plev) {
//     AoS& particles = kv.second.GetArrayOfStructs();
//     for (auto& p : particles) {
//       if (p.id() > 0)
// 	p.rdata(pstate_idx) = val;
//     }
//   }
//   return;
// }

void
SprayParticleContainer::moveKick (MultiFab&   state,
				  MultiFab&   source,
				  const int   lev,
				  const Real& dt,
				  const Real  time,
				  const bool  isVirtual,
				  const bool  isGhost,
				  const int   tmp_src_width)
{
  bool do_move = false;
  int width = 0;
  moveKickDrift(state, source, lev, dt, time, isVirtual, isGhost,
		tmp_src_width, do_move, width);
}

void
SprayParticleContainer::moveKickDrift (MultiFab&   state,
				       MultiFab&   source,
				       const int   lev,
				       const Real& dt,
				       const Real  time,
				       const bool  isVirtual,
				       const bool  isGhost,
				       const int   tmp_src_width,
				       const bool  do_move,
				       const int   where_width)
{
  BL_PROFILE("ParticleContainer::moveKickDrift()");
  AMREX_ASSERT(lev >= 0);
  AMREX_ASSERT(state.nGrow() >= 2);

  //If there are no particles at this level
  if (lev >= this->GetParticles().size())
    return;

  const Real strttime = ParallelDescriptor::second();

  MultiFab* state_ptr;

  // ********************************************************************************
  // We only make a new state_ptr if the boxArray of the state differs from the
  // boxArray onto which the particles are decomposed
  // ********************************************************************************
  bool tempState = false;

  if (this->OnSameGrids(lev, state)) {
    state_ptr = &state;
  } else {
    state_ptr = new MultiFab(this->m_gdb->ParticleBoxArray(lev),
                             this->m_gdb->ParticleDistributionMap(lev),
                             state.nComp(), state.nGrow());
    state_ptr->setVal(0.);
    state_ptr->copy(state,0,0,state.nComp());
    state_ptr->FillBoundary(Geom(lev).periodicity());
    tempState = true;
  }

  // ********************************************************************************
  // We make a temporary MultiFab for the source here and initialize it to zero
  // because if we use one that already has values in it, the SumBoundary call
  // will end up duplicating those values with each call
  // ********************************************************************************
  bool tempSource = false;

  MultiFab* tmp_src_ptr;
  if (this->OnSameGrids(lev, source) && !isVirtual && !isGhost
      && source.nGrow() >= tmp_src_width) {
    tmp_src_ptr = &source;
  } else {
    tmp_src_ptr = new MultiFab(this->m_gdb->ParticleBoxArray(lev),
                               this->m_gdb->ParticleDistributionMap(lev),
                               source.nComp(), tmp_src_width);
    tmp_src_ptr->setVal(0.);
    tempSource = true;
  }
  BL_PROFILE_VAR("SprayParticles::updateParticles()", UPD_PART);
  updateParticles(lev, (*state_ptr), (*tmp_src_ptr), dt, time, tmp_src_width, do_move);
  BL_PROFILE_VAR_STOP(UPD_PART);

  // ********************************************************************************
  // Make sure the momentum put into ghost cells of each grid is added to both
  //      valid regions AND the ghost cells of other grids.  If at level = 0
  //      we can accomplish this with SumBoundary; however if this level has ghost
  //      cells not covered at this level then we need to do more.
  // ********************************************************************************
  if (lev > 0) {
    IntVect ghostVect(tmp_src_width*IntVect::TheUnitVector());
    tmp_src_ptr->
      SumBoundary(0, tmp_src_ptr->nComp(), ghostVect, Geom(lev).periodicity());
  } else {
      tmp_src_ptr->SumBoundary(Geom(lev).periodicity());
  }

  // Add new sources into source *after* we have called SumBoundary
  if (tempSource) {
    MultiFab::Add(source, *tmp_src_ptr, 0, 0, source.nComp(),
                  amrex::min(source.nGrow(), tmp_src_width));
    delete tmp_src_ptr;
  }

  // Fill ghost cells after we've synced up ..
  // TODO: Check to see if this is needed at all
  if (lev > 0)
    source.FillBoundary(Geom(lev).periodicity());

  // Only delete this if in fact we created it.  Note we didn't change state_ptr so
  //  we don't need to copy anything back into state
  if (tempState) delete state_ptr;

  // ********************************************************************************

  if (lev > 0 && sub_cycle && do_move && !isVirtual) {
    ParticleLocData pld;
    auto& pmap = GetParticles(lev);
    for (auto& kv : pmap) {
      auto& pbox = kv.second.GetArrayOfStructs();
      const long Np = pbox.numParticles();
      for (int k = 0; k < Np; ++k) {
        ParticleType& p = pbox[k];
        //  TODO: Double check this for correctness and figure out what it is doing
        if (p.id() > 0) {
          if (!this->Where(p, pld, lev, lev, where_width)) {
            if (p.id() == GhostParticleID) {
              p.id() = -1;
            } else {
              Abort("Trying to remove non-ghost particle");
            }
          }
        }
      }
    }
  }

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
SprayParticleContainer::estTimestep (int lev, Real cfl) const
{
  BL_PROFILE("ParticleContainer::estTimestep()");
  // TODO: Clean up this mess and bring the num particle functionality back
  Real dt = std::numeric_limits<Real>::max();
  if (lev >= this->GetParticles().size() || PeleC::particle_mom_tran == 0)
    return dt;

  const Real strttime = ParallelDescriptor::second();
  const Geometry& geom = this->m_gdb->Geom(lev);
  const auto dxi = geom.InvCellSizeArray();
  {
    amrex::ReduceOps<amrex::ReduceOpMin> reduce_op;
    amrex::ReduceData<amrex::Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    for (MyParConstIter pti(*this, lev); pti.isValid(); ++pti) {
      const AoS& pbox = pti.GetArrayOfStructs();
      const ParticleType* pstruct = pbox().data();
      const int n = pbox.numParticles();
#ifdef USE_SPRAY_SOA
      auto& attribs = pti.GetAttribs();
      AMREX_D_TERM(const Real * up = attribs[0].data();,
                   const Real * vp = attribs[1].data();,
                   const Real * wp = attribs[2].data(););
#endif
      reduce_op.eval(n, reduce_data,
                     [=] AMREX_GPU_DEVICE (const int i) -> ReduceTuple
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
SprayParticleContainer::updateParticles(const int&  lev,
                                        MultiFab&   state,
                                        MultiFab&   source,
                                        const Real& flow_dt,
                                        const Real& time,
                                        const int   numGhost,
                                        const bool  do_move)
{
  AMREX_ASSERT(OnSameGrids(lev, state));
  AMREX_ASSERT(OnSameGrids(lev, source));
  const auto dxi = this->Geom(lev).InvCellSizeArray();
  const auto plo = this->Geom(lev).ProbLoArray();
  const auto phi = this->Geom(lev).ProbHiArray();
  const auto domain = this->Geom(lev).Domain();
  IntVect dom_lo = domain.smallEnd();
  IntVect dom_hi = domain.bigEnd();
  // Vector to determine boundary type: 1 - reflective,
  // -1 - not reflective, 0 - periodic
  int lo_bound[AMREX_SPACEDIM];
  int hi_bound[AMREX_SPACEDIM];
  // We are only concerned with domain boundaries that are reflective
  for (int dir = 0; dir != AMREX_SPACEDIM; ++dir) {
    if (!this->Geom(lev).isPeriodic(dir)) {
      if (reflect_lo[dir]) {
        lo_bound[dir] = 1;
      } else {
        lo_bound[dir] = -1;
      }
      if (reflect_hi[dir]) {
        hi_bound[dir] = 1;
      } else {
        hi_bound[dir] = -1;
      }
    } else {
      lo_bound[dir] = 0;
      hi_bound[dir] = 0;
    }
    // Only concered with reflective boundaries
    if (!reflect_lo[dir]) dom_lo[dir] -= 100;
    if (!reflect_hi[dir]) dom_hi[dir] += 100;
  }
  const auto dx = this->Geom(lev).CellSizeArray();
  const Real vol = AMREX_D_TERM(dx[0],*dx[1],*dx[2]);
  const Real inv_vol = 1./vol;
  // Set all constants
  const Real C_eps = 1.E-15;
  const Real B_eps = 1.E-7;
  const Real dia_eps = 2.E-6;
  const int nSubMax = 100;
  const Real third = 1./3.;
  const Real rule = third; // This determines how average values for the vapor are approximated
  const Real Pi_six = M_PI/6.;
  Real mw_fluid[NUM_SPECIES];
  Real invmw[NUM_SPECIES];
  EOS::molecular_weight(mw_fluid);
  EOS::inv_molecular_weight(invmw);
  // Extract control parameters for mass, heat, and momentum transfer
  const int heat_trans = PeleC::particle_heat_tran;
  const int mass_trans = PeleC::particle_mass_tran;
  const int mom_trans = PeleC::particle_mom_tran;
  const Real inv_Ru = 1./EOS::RU;
  const Real ref_T = PeleC::sprayRefT;
  // Particle components indices
  const int pstateVel = PeleC::pstateVel;
  const int pstateT   = PeleC::pstateT;
  const int pstateRho = PeleC::pstateRho;
  const int pstateDia = PeleC::pstateDia;
  const int pstateY   = PeleC::pstateY;
  bool get_xi = false;
  bool get_Ddiag = true;
  bool get_lambda = true;
  bool get_mu = true;
  if (!mass_trans && !heat_trans) {
    get_Ddiag = false;
    get_lambda = false;
  }

  // Component indices for conservative state
  const int rhoIndx = PeleC::Density;
  const int momIndx = PeleC::Xmom;
  const int engIndx = PeleC::Eden;
  const int specIndx = PeleC::FirstSpec;
  // Start the ParIter, which loops over separate sets of particles in different boxes
  for (MyParIter pti(*this, lev); pti.isValid(); ++pti) {
    const Box& tile_box = pti.growntilebox(numGhost);
    const long Np = pti.numParticles();
    ParticleType* pstruct = &(pti.GetArrayOfStructs()[0]);
    // Get particle attributes if StructOfArrays are used
#ifdef USE_SPRAY_SOA
    auto& attribs = pti.GetAttribs();
    Real * velp[AMREX_SPACEDIM];
    AMREX_D_TERM(velp[0] = attribs[pstateVel].data();,
                 velp[1] = attribs[pstateVel+1].data();,
                 velp[2] = attribs[pstateVel+2].data(););
    Real * Tp = attribs[pstateT].data();
    Real * diap = attribs[pstateDia].data();
    Real * rhop = attribs[pstateRho].data();
    std::array<Real *, SPRAY_FUEL_NUM> Yp;
    for (int spf = 0; spf != SPRAY_FUEL_NUM; ++spf) {
      Yp[spf] = attribs[pstateY+spf].data();
    }
#endif
    auto const crit_T = m_fuelData.critT();
    auto const boil_T = m_fuelData.boilT();
    auto const invBoilT = m_fuelData.invBoilT();
    auto const fuel_cp = m_fuelData.fuelCp();
    auto const fuel_latent = m_fuelData.fuelLatent();
    auto const fuel_indx = m_fuelData.fuelIndx();
    Array4<const Real> const& statearr = state.array(pti);
    Array4<Real> const& sourcearr = source.array(pti);
    AMREX_FOR_1D ( Np, i,
    {
      ParticleType& p = pstruct[i];
      if (p.id() > 0) {
        Real dt = flow_dt;
        Real sub_source = inv_vol;
        // TODO: I was hoping to not have to instantiate this everytime
        Real Y_fluid[NUM_SPECIES];
        Real Y_skin[NUM_SPECIES];
        Real h_skin[NUM_SPECIES];
        Real cp_n[NUM_SPECIES];
        Real mass_frac[NUM_SPECIES];
        Real Ddiag[NUM_SPECIES];
        Real B_M_num[SPRAY_FUEL_NUM];
        Real Sh_num[SPRAY_FUEL_NUM];
        Real Y_dot[SPRAY_FUEL_NUM];
        Real L_fuel[SPRAY_FUEL_NUM];
        // Weights for interpolation
        Real coef[AMREX_D_PICK(2, 4, 8)];
        // Indices of adjacent cells
        IntVect indx_array[AMREX_D_PICK(2, 4, 8)];
        // Storage for fluid info in cells adjacent to the particle
        RealVect len(AMREX_D_DECL((p.pos(0) - plo[0])*dxi[0] + 0.5,
                                  (p.pos(1) - plo[1])*dxi[1] + 0.5,
                                  (p.pos(2) - plo[2])*dxi[2] + 0.5));
        // Do initial interpolation and save corresponding adjacent fluid information
        AdjIndexWeights(len, indx_array, coef, dom_lo, dom_hi);
        RealVect vel_fluid(RealVect::TheZeroVector());
        Real T_fluid = 0.;
        Real rho_fluid = 0.;
        Real T_val = 300.;
        for (int sp = 0; sp != NUM_SPECIES; ++sp) Y_fluid[sp] = 0.;
        // Extract adjacent values and interpolate fluid at particle location
        for (int aindx = 0; aindx != AMREX_D_PICK(2, 4, 8); ++aindx) {
          IntVect cur_indx = indx_array[aindx];
          AMREX_ASSERT(tile_box.contains(cur_indx));
          Real cur_coef = coef[aindx];
          Real cur_rho = statearr(cur_indx, rhoIndx);
          rho_fluid += cur_coef*cur_rho;
          Real inv_rho = 1./cur_rho;
          AMREX_D_TERM(Real velx = statearr(cur_indx, momIndx)*inv_rho;
                       Real ke = 0.5*velx*velx;
                       vel_fluid[0] += cur_coef*velx;,
                       Real vely = statearr(cur_indx, momIndx+1)*inv_rho;
                       ke += 0.5*vely*vely;
                       vel_fluid[1] += cur_coef*vely;,
                       Real velz = statearr(cur_indx, momIndx+2)*inv_rho;
                       ke += 0.5*velz*velz;
                       vel_fluid[2] += cur_coef*velz;);
          for (int sp = 0; sp != NUM_SPECIES; ++sp) {
            int mf_indx = sp + specIndx;
            Real cur_mf = statearr(cur_indx, mf_indx)*inv_rho;
            Y_fluid[sp] += cur_coef*cur_mf;
            mass_frac[sp] = cur_mf;
          }
          Real intEng = statearr(cur_indx, engIndx)*inv_rho - ke;
          EOS::EY2T(intEng, mass_frac, T_val);
          T_fluid += cur_coef*T_val;
        }
        int isub = 1; // Initialize the number of sub-cycles
        int nsub = 1; // This is set in the first run through the loop
        while (isub <= nsub) {
#ifdef USE_SPRAY_SOA
          RealVect vel_part(AMREX_D_DECL(velp[0][i], velp[1][i], velp[2][i]));
          Real T_part = Tp[i];
          Real dia_part = diap[i];
          Real rho_part = rhop[i];
#else
          RealVect vel_part(AMREX_D_DECL(p.rdata(pstateVel),
                                         p.rdata(pstateVel+1),
                                         p.rdata(pstateVel+2)));
          Real T_part = p.rdata(pstateT);
          Real dia_part = p.rdata(pstateDia);
          Real rho_part = p.rdata(pstateRho);
#endif
          Real dia2_part = dia_part*dia_part;
          Real pmass = Pi_six*rho_part*dia_part*dia2_part;
          Real part_ke = 0.5*vel_part.radSquared();
          // If multiple sub-cycle iterations are needed, we might need to
          // re-interpolate values at the particle
          // However, since we don't allow the particle to move more than
          // 5% of a cell (particle cfl = 0.05), the adjacent fluid values
          // should remain the same
          // Model the fuel vapor using the one-third rule
          Real delT = amrex::max(T_fluid - T_part, 0.);
          Real T_skin = T_part + rule*delT;
          // Calculate the C_p at the skin temperature for each species
          EOS::T2Cpi(T_skin, cp_n);
          EOS::T2Hi(T_part, h_skin);
          Real mw_mix = 0.;  // Average molar mass of gas mixture
          if (heat_trans || mass_trans) {
            for (int sp = 0; sp != NUM_SPECIES; ++sp) {
              mw_mix += Y_fluid[sp]*invmw[sp];
              Y_skin[sp] = 0.;
            }
          } else {
            for (int sp = 0; sp != NUM_SPECIES; ++sp) {
              mw_mix += Y_fluid[sp]*invmw[sp];
              Y_skin[sp] = Y_fluid[sp];
            }
          }
          Real p_fluid = rho_fluid*EOS::RU*mw_mix*T_fluid;
          mw_mix = 1./mw_mix;

          // Solve for state of the vapor and mass transfer coefficient B_M
          Real sumYSkin = 0.; // Mass fraction of the fuel in skin film, uses one-thirds rule
          Real sumYFuel = 0.; // Mass fraction of the fuel in the gas phase
          Real cp_skin = 0.; // Averaged C_p at particle surface
          Real cp_L_av = 0.;  // Cp of the liquid state
          if (heat_trans || mass_trans) {
            for (int spf = 0; spf != SPRAY_FUEL_NUM; ++spf) {
              const int fspec = fuel_indx[spf];
              const Real mw_fuel = mw_fluid[fspec];
              // Compute latent heat
#ifdef LEGACY_SPRAY
              Real part_latent = fuel_latent[spf]*
                std::pow(amrex::max((crit_T[spf] - T_part)/
                                    (crit_T[spf] - boil_T[spf]), 0.), 0.38);
#else
              Real part_latent = h_skin[fspec] + fuel_latent[spf]
                - fuel_cp[spf]*(T_part - ref_T);
#endif
              L_fuel[spf] = part_latent;
              // Compute the mass fraction of the fuel vapor at droplet surface
              Real pres_sat = EOS::PATM*std::exp(part_latent*inv_Ru*mw_fuel*
                                                 (invBoilT[spf] - 1./T_part)) + C_eps;
              Real Yfv = mw_fuel*pres_sat/(mw_mix*p_fluid + (mw_fuel - mw_mix)*pres_sat);
              Yfv = amrex::max(0., amrex::min(1. - C_eps, Yfv));
              B_M_num[spf] = amrex::max(C_eps, (Yfv - Y_fluid[fspec])/(1. - Yfv));
#ifndef LEGACY_SPRAY
              Y_skin[fspec] = Yfv + rule*(Y_fluid[fspec] - Yfv);
              sumYSkin += Y_skin[fspec];
#endif
#ifdef USE_SPRAY_SOA
              cp_L_av += Yp[spf][i]*fuel_cp[spf];
#else
              cp_L_av += p.rdata(pstateY+spf)*fuel_cp[spf];
#endif
              sumYFuel += Y_fluid[fspec];
            }
            const Real restYSkin = 1. - sumYSkin;
            for (int sp = 0; sp != NUM_SPECIES; ++sp) {
#ifdef LEGACY_SPRAY
              Y_skin[sp] = Y_fluid[sp];
#else
              Y_skin[sp] += restYSkin*Y_fluid[sp];
#endif
              cp_skin += Y_skin[sp]*cp_n[sp];
            }
          }
          Real lambda_skin = 0.;
          Real mu_skin = 0.;
          Real xi_skin = 0.;
          transport(get_xi, get_mu, get_lambda, get_Ddiag,
                    T_fluid, rho_fluid, Y_skin, Ddiag,
                    mu_skin, xi_skin, lambda_skin);
          // Ensure gas is not all fuel to allow evaporation
          bool evap_fuel = (sumYFuel >= 1.) ? false : true;
          RealVect diff_vel = vel_fluid - vel_part;
          Real diff_vel_mag = diff_vel.vectorLength();
          // Local Reynolds number
          Real Reyn = rho_fluid*diff_vel_mag*dia_part/mu_skin;
          Real Nu_0 = 1.;

          // Solve mass transfer source terms
          Real m_dot = 0.;
          Real d_dot = 0.;
          if ((mass_trans || heat_trans) && evap_fuel) {
            Real Pr_skin = mu_skin*cp_skin/lambda_skin;
            Real powR = amrex::max(std::pow(Reyn, 0.077), 1.);
            Nu_0 = 1. + powR*std::cbrt(1. + Reyn*Pr_skin);
            for (int spf = 0; spf != SPRAY_FUEL_NUM; ++spf) {
              const int fspec = fuel_indx[spf];
              const Real rhoD = Ddiag[fspec];
              const Real Sc_skin = mu_skin/rhoD;
              const Real B_M = B_M_num[spf];
              Real logB = std::log(1. + B_M);
              // Calculate Sherwood number and evaporation rate
              Real invFM = B_M/(logB*std::pow(1. + B_M, 0.7));
              Real Sh_0 = 1. + powR*std::cbrt(1. + Reyn*Sc_skin);
              Sh_num[spf] = 2. + (Sh_0 - 2.)*invFM;
              if (mass_trans) {
                Y_dot[spf] = -amrex::max(M_PI*rhoD*dia_part*Sh_num[spf]*logB, 0.);
                m_dot += Y_dot[spf];
              }
            }
            d_dot = m_dot/(0.5*M_PI*rho_part*dia2_part);
          }

          // Solve for momentum source terms
          const Real inv_pmass = 1./pmass;
          RealVect fluid_mom_src(RealVect::TheZeroVector());
          RealVect part_mom_src(RealVect::TheZeroVector());
          Real fluid_eng_src = 0.;
          if (mom_trans) {
#ifdef LEGACY_SPRAY
            Real drag_force = 3.*M_PI*mu_skin*dia_part*(1. + 0.15*std::pow(Reyn, 0.687));
#else
            Real drag_coef =
              (Reyn > 1.) ? 24./Reyn*(1. + std::cbrt(Reyn*Reyn)/6.) : 24./Reyn;
            Real drag_force = 0.125*rho_fluid*drag_coef*M_PI*dia2_part*diff_vel_mag;
#endif
            part_mom_src = drag_force*diff_vel;
            fluid_mom_src = part_mom_src + vel_part*m_dot;
            // s_d,mu dot u_d
            Real S_dmu_dot_u = part_mom_src.dotProduct(vel_part);
            fluid_eng_src += S_dmu_dot_u + m_dot*part_ke;
            Real inv_tau_var = drag_force*inv_pmass;
            if (isub == 1)
              nsub = amrex::min(int(flow_dt*inv_tau_var) + 1, nSubMax);
          }

          // Solve for energy source terms
          Real part_temp_src = 0.;
          if (heat_trans && evap_fuel) {
            const Real inv_pm_cp = inv_pmass/cp_L_av;
            Real coeff_heat = 0.;
            for (int spf = 0; spf != SPRAY_FUEL_NUM; ++spf) {
              const int fspec = fuel_indx[spf];
              Real ratio = cp_n[fspec]*Sh_num[spf]*Ddiag[fspec]/lambda_skin;
              Real heatC = calcHeatCoeff(ratio, B_M_num[spf], B_eps, C_eps, Nu_0);
              // Convection term
              coeff_heat += heatC;
#ifdef LEGACY_SPRAY
              Real htmp = Y_dot[spf]*(cp_n[fspec]*(T_skin - T_part) + L_fuel[spf]);
              fluid_eng_src += htmp;
              part_temp_src += htmp;
#else
              fluid_eng_src += Y_dot[spf]*h_skin[fspec];
              part_temp_src += Y_dot[spf]*L_fuel[spf];
#endif
            }
            Real conv_src = M_PI*lambda_skin*dia_part*delT*coeff_heat;
            fluid_eng_src += conv_src;
            part_temp_src += conv_src;
            part_temp_src *= inv_pm_cp;
            if (isub == 1 && delT > C_eps) {
              Real inv_tau_T = conv_src*inv_pm_cp/delT;
              nsub = amrex::min(amrex::max(nsub, int(flow_dt*inv_tau_T) + 1), nSubMax);
            }
          }
          if (isub == 1) {
            sub_source /= Real(nsub);
            dt = flow_dt/Real(nsub);
          }
          const Real part_dt = 0.5*dt;
          if (mom_trans || mass_trans || heat_trans) {
            bool remove_particle = false;
            // Modify particle velocity by half time step
#ifdef USE_SPRAY_SOA
            if (mom_trans) {
              AMREX_D_TERM(Gpu::Atomic::Add(&velp[0][i], part_dt*part_mom_src[0]*inv_pmass);,
                           Gpu::Atomic::Add(&velp[1][i], part_dt*part_mom_src[1]*inv_pmass);,
                           Gpu::Atomic::Add(&velp[2][i], part_dt*part_mom_src[2]*inv_pmass););
              if (do_move) {
                for (int dir = 0; dir != AMREX_SPACEDIM; ++dir) {
                  Gpu::Atomic::Add(&p.pos(dir), dt*velp[dir][i]);
                  remove_particle = checkWall(p.pos(dir), velp[dir][i],
                                              phi[dir], plo[dir], hi_bound[dir], lo_bound[dir]);
                }
              }
            }
            if (heat_trans) {
              // Modify particle temperature
              Gpu::Atomic::Add(&Tp[i], part_dt*part_temp_src);
            }
            if (mass_trans) {
              // Compute new particle diameter
              Real new_dia = dia_part + part_dt*d_dot;
              if (new_dia < dia_eps) {
                remove_particle = true;
              } else {
                diap[i] = new_dia;
              }
            }
#else
            if (mom_trans) {
              AMREX_D_TERM(Gpu::Atomic::Add(&p.rdata(pstateVel), part_dt*part_mom_src[0]*inv_pmass);,
                           Gpu::Atomic::Add(&p.rdata(pstateVel+1), part_dt*part_mom_src[1]*inv_pmass);,
                           Gpu::Atomic::Add(&p.rdata(pstateVel+2), part_dt*part_mom_src[2]*inv_pmass););
              // Modify particle position by whole time step
              if (do_move) {
                for (int dir = 0; dir != AMREX_SPACEDIM; ++dir) {
                  Gpu::Atomic::Add(&p.pos(dir), dt*p.rdata(pstateVel+dir));
                  remove_particle = checkWall(p.pos(dir), p.rdata(pstateVel+dir),
                                              phi[dir], plo[dir], hi_bound[dir], lo_bound[dir]);
                }
              }
            }
            if (heat_trans) {
              // Modify particle temperature
              Gpu::Atomic::Add(&p.rdata(pstateT), part_dt*part_temp_src);
            }
            if (mass_trans) {
              // Compute new particle diameter
              Real new_dia = dia_part + part_dt*d_dot;
              if (new_dia < dia_eps) {
                remove_particle = true;
              } else {
                p.rdata(pstateDia) = new_dia;
              }
            }
#endif // End of check for SOA
            if (remove_particle) {
              p.id() = -1;
              isub = nsub + 1;
              continue;
            }
            // Add the spray source terms to the Eulerian fluid
            for (int aindx = 0; aindx != AMREX_D_PICK(2, 4, 8); ++aindx) {
              Real cur_coef = -coef[aindx]*sub_source;
              IntVect cur_indx = indx_array[aindx];
              AMREX_ASSERT(tile_box.contains(cur_indx));
              if (mom_trans) {
                for (int dir = 0; dir != AMREX_SPACEDIM; ++dir) {
                  const int nf = momIndx + dir;
                  Gpu::Atomic::Add(&sourcearr(cur_indx, nf), cur_coef*fluid_mom_src[dir]);
                }
              }
              if (mass_trans) {
                Gpu::Atomic::Add(&sourcearr(cur_indx, rhoIndx), cur_coef*m_dot);
                for (int spf = 0; spf != SPRAY_FUEL_NUM; ++spf) {
                  const int nf = specIndx + fuel_indx[spf];
                  Gpu::Atomic::Add(&sourcearr(cur_indx, nf), cur_coef*Y_dot[spf]);
                }
              }
              if (mass_trans || heat_trans)
                Gpu::Atomic::Add(&sourcearr(cur_indx, engIndx), cur_coef*fluid_eng_src);
            }
          }
          // Increment sub-iteration
          ++isub;
        } // End of isub while loop
      } // End of p.id() > 0 check
    }); // End of loop over particles
  }
}
