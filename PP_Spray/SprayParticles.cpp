
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
  for (int dir = 0; dir < BL_SPACEDIM; dir++) {
    if (phys_bc->lo(dir) == Symmetry   ||
	phys_bc->lo(dir) == SlipWall   ||
	phys_bc->lo(dir) == NoSlipWall) {
      reflect_lo[dir] = 1;
    } else {
      reflect_lo[dir] = 0;
    }
    if (phys_bc->hi(dir) == Symmetry   ||
	phys_bc->hi(dir) == SlipWall   ||
	phys_bc->hi(dir) == NoSlipWall) {
      reflect_hi[dir] = 1;
    } else {
      reflect_hi[dir] = 0;
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
				  const int   tmp_src_width)
{
  bool do_move = false;
  int width = 0;
  moveKickDrift(state, source, lev, dt, time, tmp_src_width,
		do_move, width);
}

void
SprayParticleContainer::moveKickDrift (MultiFab&   state,
				       MultiFab&   source,
				       const int   lev,
				       const Real& dt,
				       const Real  time,
				       const int   tmp_src_width,
				       const bool  do_move,
				       const int   where_width)
{  
  BL_PROFILE("ParticleContainer::moveKickDrift()");
  BL_ASSERT(lev >= 0);
  BL_ASSERT(state.nGrow() >= 2);

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

  MultiFab* tmp_src_ptr = new MultiFab(this->m_gdb->ParticleBoxArray(lev),
				       this->m_gdb->ParticleDistributionMap(lev),
				       source.nComp(), tmp_src_width);
  tmp_src_ptr->setVal(0.);
  updateParticles(lev, (*state_ptr), (*tmp_src_ptr), dt, time, do_move);

  // ********************************************************************************
  // Make sure the momentum put into ghost cells of each grid is added to both 
  //      valid regions AND the ghost cells of other grids.  If at level = 0 
  //      we can accomplish this with SumBoundary; however if this level has ghost
  //      cells not covered at this level then we need to do more.
  // ********************************************************************************
  if (lev > 0) {
    int ncomp = tmp_src_ptr->nComp();
    MultiFab tmp(this->m_gdb->ParticleBoxArray(lev),
		 this->m_gdb->ParticleDistributionMap(lev),
		 ncomp, tmp_src_ptr->nGrow());
    tmp.setVal(0.);
    tmp.copy(*tmp_src_ptr, 0, 0, ncomp, tmp_src_ptr->nGrow(), tmp_src_ptr->nGrow(),
	     Geom(lev).periodicity(), FabArrayBase::ADD);
    tmp_src_ptr->copy(tmp, 0, 0, ncomp, tmp_src_width, tmp_src_width,
		      Geom(lev).periodicity(), FabArrayBase::COPY);
  } else {
    tmp_src_ptr->SumBoundary(Geom(lev).periodicity());
  }

  // Add new sources into source *after* we have called SumBoundary
  MultiFab::Add(source,*tmp_src_ptr,0,0, source.nComp(),
		amrex::min(source.nGrow(), tmp_src_width));
  delete tmp_src_ptr;

  // Fill ghost cells after we've synced up ..
  source.FillBoundary(Geom(lev).periodicity());

  // Only delete this if in fact we created it.  Note we didn't change state_ptr so 
  //  we don't need to copy anything back into state
  if (tempState) delete state_ptr;

  // ********************************************************************************

  if (lev > 0 && sub_cycle && do_move) {
    ParticleLocData pld;
    auto& pmap = GetParticles(lev);
    for (auto& kv : pmap) {
      auto& pbox = kv.second.GetArrayOfStructs();
      const int Np = pbox.size();
      for (int k = 0; k < Np; ++k) {
  	ParticleType& p = pbox[k];
  	//  TODO: Double check this for correctness and figure out what it is doing
  	if (p.id() <= 0) continue;
  	if (!this->Where(p, pld, lev, lev, where_width)) {
  	  if (p.id() == GhostParticleID) {
  	    p.id() = -1;
  	  } else {
  	    p.id() = -1;
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
  // TODO: Clean up this mess and bring the num particle functionality back
  Real dt = std::numeric_limits<Real>::max();

  if (this->GetParticles().size() == 0 || PeleC::particle_mom_tran == 0)
    return dt;

  const Real strttime = ParallelDescriptor::second();
  const Geometry& geom = this->m_gdb->Geom(lev);
  const auto dxi = geom.InvCellSizeArray();
  // TODO: This assumes that pstateVel = 0 and dxi[0] = dxi[1] = dxi[2]
  using PType = typename SprayParticleContainer::SuperParticleType;
  Real max_mag_vdx = ReduceMax(*this, lev, [=] AMREX_GPU_HOST_DEVICE (const PType& p) ->
  			       Real { return amrex::max(AMREX_D_DECL(std::abs(p.rdata(0)),
  								     std::abs(p.rdata(1)),
  								     std::abs(p.rdata(2))));})*dxi[0];
  Real global_dt = (max_mag_vdx > 0.) ? (cfl/max_mag_vdx) : 1.E50;
  dt = global_dt;

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
					const bool  do_move)
{
  AMREX_ASSERT(OnSameGrids(lev, state));
  AMREX_ASSERT(OnSameGrids(lev, source));
  const auto dxi = this->Geom(lev).InvCellSizeArray();
  const auto plo = this->Geom(lev).ProbLoArray();
  const auto domain = this->Geom(lev).Domain();
  const auto dx = this->Geom(lev).CellSizeArray();
  const Real vol = AMREX_D_TERM(dx[0],*dx[1],*dx[2]);
  const Real inv_vol = 1./vol;
  // Set all constants
  const Real eps = 1.E-15;
  const Real B_eps = 1.E-7;
  const Real dia_eps = 2.E-6;
  const int nSubMax = 100;
  const Real third = 1./3.;
  const Real rule = third; // This determines how average values for the vapor are approximated
  const Real Pi_six = M_PI/6.;
  Real mw_fluid[NUM_SPECIES];
  Real invmw[NUM_SPECIES];
  EOS::get_mw_eos(mw_fluid);
  EOS::get_imw_eos(invmw);
  // Extract control parameters for mass, heat, and momentum transfer
  const int heat_trans = PeleC::particle_heat_tran;
  const int mass_trans = PeleC::particle_mass_tran;
  const int mom_trans = PeleC::particle_mom_tran;
  const Real p0 = 1.01325E6;
  const Real inv_Ru = 1./RU;
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
    int Np = pti.numParticles();
    ParticleType* pstruct = &(pti.GetArrayOfStructs()[0]);
    auto const ref_h = m_fuelData.refFuelH();
    auto const crit_T = m_fuelData.critT();
    auto const boil_T = m_fuelData.boilT();
    auto const invBoilT = m_fuelData.invBoilT();
    auto const fuel_cp = m_fuelData.fuelCp();
    auto const fuel_latent = m_fuelData.fuelLatent();
    auto const fuel_indx = m_fuelData.fuelIndx();
    auto const statearr = state.array(pti);
    auto sourcearr = source.array(pti);
    AMREX_FOR_1D ( Np, i,
    {
      ParticleType& p = pstruct[i];
      if (p.id() == -1) continue;
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
      Real rho_fluid_interp[AMREX_D_PICK(2, 4, 8)];
      Real T_fluid_interp[AMREX_D_PICK(2, 4, 8)];
      Real vel_fluid_interp[AMREX_D_PICK(2, 4, 8)][AMREX_SPACEDIM];
      Real Y_fluid_interp[AMREX_D_PICK(2, 4, 8)][NUM_SPECIES];
      RealVect len(AMREX_D_DECL((p.pos(0) - plo[0])*dxi[0] + 0.5,
				(p.pos(1) - plo[1])*dxi[1] + 0.5,
				(p.pos(2) - plo[2])*dxi[2] + 0.5));
      // Do initial interpolation and save corresponding adjacent fluid information
      AdjIndexWeights(len, indx_array, coef);
      RealVect vel_fluid = RealVect::TheZeroVector();
      Real T_fluid = 0.;
      Real rho_fluid = 0.;
      for (int sp = 0; sp != NUM_SPECIES; ++sp) Y_fluid[sp] = 0.;
      // Extract adjacent values and interpolate fluid at particle location
      for (int aindx = 0; aindx != AMREX_D_PICK(2, 4, 8); ++aindx) {
	IntVect cur_indx = indx_array[aindx];
	Real cur_coef = coef[aindx];
	Real cur_rho = statearr(cur_indx, rhoIndx);
	rho_fluid_interp[aindx] = cur_rho;
	rho_fluid += cur_coef*cur_rho;
	Real inv_rho = 1./cur_rho;
	Real ke = 0.;
	for (int dir = 0; dir != AMREX_SPACEDIM; ++dir) {
	  int nf = momIndx + dir;
	  Real vel = statearr(cur_indx, nf)*inv_rho;
	  vel_fluid_interp[aindx][dir] = vel;
	  vel_fluid[dir] += cur_coef*vel;
	  ke += 0.5*vel*vel;
	}
	for (int sp = 0; sp != NUM_SPECIES; ++sp) {
	  int mf_indx = sp + specIndx;
	  Real cur_mf = statearr(cur_indx, mf_indx)*inv_rho;
	  Y_fluid_interp[aindx][sp] = cur_mf;
	  Y_fluid[sp] += cur_coef*cur_mf;
	  mass_frac[sp] = cur_mf;
	}
	Real intEng = statearr(cur_indx, engIndx)*inv_rho - ke;
	Real T_val;
	EOS::cmpT(intEng, mass_frac, T_val);
	T_fluid_interp[aindx] = T_val;
	T_fluid += cur_coef*T_val;
      }
      int isub = 1; // Initialize the number of sub-cycles
      int nsub = 1; // This is set in the first run through the loop
      while (isub <= nsub) {
	RealVect vel_part(AMREX_D_DECL(p.rdata(pstateVel),
				       p.rdata(pstateVel+1),
				       p.rdata(pstateVel+2)));
	Real T_part = p.rdata(pstateT);
	Real rho_part = p.rdata(pstateRho);
	Real dia_part = p.rdata(pstateDia);
	Real dia2_part = dia_part*dia_part;
	Real pmass = Pi_six*rho_part*dia_part*dia2_part;
	Real part_ke = 0.5*vel_part.radSquared();
	// If multiple sub-cycle iterations are needed, we might need to
	// re-interpolate values at the particle
	// However, since we don't allow the particle to move more than
	// 5% of a cell (particle cfl = 0.05), the adjacent fluid values
	// should remain the same
	if (isub > 1 && do_move) {
	  // Reset values to zero for interpolation
	  for (int sp = 0; sp != NUM_SPECIES; ++sp) Y_fluid[sp] = 0.;
	  rho_fluid = 0.;
	  T_fluid = 0.;
	  vel_fluid = RealVect::TheZeroVector();
	  AMREX_D_TERM(len[0] = (p.pos(0) - plo[0])*dxi[0] + 0.5;,
		       len[1] =	(p.pos(1) - plo[1])*dxi[1] + 0.5;,
		       len[2] = (p.pos(2) - plo[2])*dxi[2] + 0.5;);
	  // Recompute the weights
	  AdjIndexWeights(len, indx_array, coef);
	  // Re-interpolate the fluid state at the new particle location
	  for (int aindx = 0; aindx != AMREX_D_PICK(2, 4, 8); ++aindx) {
	    Real cur_coef = coef[aindx];
	    rho_fluid += cur_coef*rho_fluid_interp[aindx];
	    T_fluid += cur_coef*T_fluid_interp[aindx];
	    AMREX_D_TERM(vel_fluid[0] += cur_coef*vel_fluid_interp[aindx][0];,
	  		 vel_fluid[1] += cur_coef*vel_fluid_interp[aindx][1];,
	  		 vel_fluid[2] += cur_coef*vel_fluid_interp[aindx][2];);
	    for (int sp = 0; sp != NUM_SPECIES; ++sp)
	      Y_fluid[sp] += cur_coef*Y_fluid_interp[aindx][sp];
	  }
	}
	// Model the fuel vapor using the one-third rule
	Real delT = amrex::max(T_fluid - T_part, 0.);
	Real T_skin = T_part + rule*delT;
	// Calculate the C_p at the skin temperature for each species
	EOS::get_cpi(T_skin, cp_n);
#ifdef PELEC_EOS_FUEGO
	EOS::get_hi(nullptr, T_part, h_skin);
#else
	amrex::Real masstmp[NUM_SPECIES] = {28.97};
	EOS::get_hi(masstmp, T_part, h_skin);
#endif
	Real cp_skin = 0.; // Averaged C_p at particle surface
	Real mw_mix = 0.;  // Average molar mass of gas mixture
	for (int sp = 0; sp != NUM_SPECIES; ++sp) {
	  Real Y_n = Y_fluid[sp];
	  mw_mix += Y_n*invmw[sp];
	  Y_skin[sp] = 0.;
	}
	Real R_fluid = RU*mw_mix;
	Real p_fluid = rho_fluid*R_fluid*T_fluid;
	mw_mix = 1./mw_mix;

	// Solve for state of the vapor and mass transfer coefficient B_M
	Real sumYSkin = 0.; // Mass fraction of the fuel in skin film, uses one-thirds rule
	Real sumYFuel = 0.; // Mass fraction of the fuel in the gas phase
	Real cp_L_av = 0.;  // Cp of the liquid state
	if (heat_trans || mass_trans) {
	  for (int spf = 0; spf != SPRAY_FUEL_NUM; ++spf) {
	    const int fspec = fuel_indx[spf];
	    // Compute latent heat
#ifdef LEGACY_SPRAY
	    Real part_latent = fuel_latent[spf]*
	      std::pow(amrex::max((crit_T[spf] - T_part)/
				  (crit_T[spf] - boil_T[spf]), 0.), 0.38);
#else
	    Real part_latent = h_skin[fspec] - ref_h[fspec]
	      + fuel_latent[spf] - fuel_cp[spf]*(T_part - ref_T);
#endif
	    L_fuel[spf] = part_latent;
	    // Compute the mass fraction of the fuel vapor at droplet surface
	    Real Yfv = calcVaporMF(part_latent, T_part, p_fluid,
				   mw_mix, mw_fluid[fspec],
				   invBoilT[spf], inv_Ru, p0, eps);
	    B_M_num[spf] = (Yfv - Y_fluid[fspec])/(1. - Yfv);
#ifndef LEGACY_SPRAY
	    Y_skin[fspec] += Yfv + rule*(Y_fluid[fspec] - Yfv);
#endif
	    sumYFuel += Y_fluid[fspec];
	    sumYSkin += Y_skin[fspec];
	    cp_L_av += p.rdata(pstateY+spf)*fuel_cp[spf];
	    Y_dot[spf] = 0.;
	  }
	  const Real restYSkin = 1. - sumYSkin;
	  for (int sp = 0; sp != NUM_SPECIES; ++sp) {
	    Y_skin[sp] += restYSkin*Y_fluid[sp];
	    cp_skin += Y_skin[sp]*cp_n[sp];
	  }
	} else {
	  for (int spf = 0; spf != SPRAY_FUEL_NUM; ++spf) Y_dot[spf] = 0.;
	  for (int sp = 0; sp != NUM_SPECIES; ++sp) Y_skin[sp] = Y_fluid[sp];
	}
	Real lambda_skin = 0.;
	Real mu_skin = 0.;
	Real xi_skin = 0.;
	pc_actual_transport(get_xi, get_mu, get_lambda, get_Ddiag,
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
	    Real logB = std::log(1. + B_M_num[spf]);
	    // Calculate Sherwood number and evaporation rate
	    Sh_num[spf] =
	      calcSpecEvapRate(dia_part, B_M_num[spf], logB, Reyn, powR, Sc_skin, rhoD);
	    Y_dot[spf] = -amrex::max(M_PI*rhoD*dia_part*Sh_num[spf]*logB, 0.);
	    m_dot += Y_dot[spf];
	  }
	  d_dot = m_dot/(0.5*M_PI*rho_part*dia2_part);
	  if (!mass_trans) {
	    d_dot = 0.;
	    m_dot = 0.;
	    for (int spf = 0; spf != SPRAY_FUEL_NUM; ++spf) Y_dot[spf] = 0.;
	  }
	}

	// Solve for momentum source terms
	RealVect fluid_mom_src = RealVect::TheZeroVector();
	RealVect part_mom_src = RealVect::TheZeroVector();
	Real fluid_eng_src = 0.;
	if (mom_trans) {
#ifdef LEGACY_SPRAY
          Real drag_force = 3.*M_PI*mu_skin*dia_part*(1. + 0.15*std::pow(Reyn, 0.687));
#else
          Real drag_coef = 24./Reyn;
          if (Reyn > 1.) drag_coef *= (1. + std::cbrt(Reyn*Reyn)/6.);
          Real drag_force = 0.125*rho_fluid*drag_coef*M_PI*dia2_part*diff_vel_mag;
#endif
	  part_mom_src = drag_force*diff_vel;
	  fluid_mom_src = part_mom_src + vel_part*m_dot;
	  // s_d,mu dot u_d
	  Real S_dmu_dot_u = part_mom_src.dotProduct(vel_part);
	  fluid_eng_src += S_dmu_dot_u + m_dot*part_ke;
	  Real inv_tau_var = drag_force/pmass;
	  if (isub == 1)
	    nsub = amrex::min(int(flow_dt*inv_tau_var) + 1, nSubMax);
	}

	// Solve for energy source terms
	Real part_temp_src = 0.;
	if (heat_trans && evap_fuel) {
	  Real coeff_heat = 0.;
	  for (int spf = 0; spf != SPRAY_FUEL_NUM; ++spf) {
	    const int fspec = fuel_indx[spf];
	    Real ratio = cp_n[fspec]*Sh_num[spf]*Ddiag[fspec]/lambda_skin;
	    Real heatC = calcHeatCoeff(ratio, B_M_num[spf], B_eps, Nu_0);
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
	  part_temp_src = (part_temp_src + conv_src)/(pmass*cp_L_av);
	  if (isub == 1) {
	    Real inv_tau_T = conv_src/(pmass*cp_L_av*delT);
	    nsub = amrex::min(amrex::max(nsub, int(flow_dt*inv_tau_T) + 1), nSubMax);
	  }
	}
	if (isub == 1) {
	  sub_source /= Real(nsub);
	  dt = flow_dt/Real(nsub);
	}
	if (mom_trans || mass_trans || heat_trans) {
	  for (int aindx = 0; aindx != AMREX_D_PICK(2, 4, 8); ++aindx) {
	    Real cur_coef = -coef[aindx]*sub_source;
	    IntVect cur_indx = indx_array[aindx];
	    for (int dir = 0; dir != AMREX_SPACEDIM; ++dir) {
	      const int nd = momIndx + dir;
	      amrex::Gpu::Atomic::Add(&sourcearr(cur_indx, nd),
	      			      cur_coef*fluid_mom_src[dir]);
	    }
	    Gpu::Atomic::Add(&sourcearr(cur_indx, rhoIndx), cur_coef*m_dot);
	    for (int spf = 0; spf != SPRAY_FUEL_NUM; ++spf) {
	      const int nf = specIndx + fuel_indx[spf];
	      Gpu::Atomic::Add(&sourcearr(cur_indx, nf), cur_coef*Y_dot[spf]);
	    }
	    Gpu::Atomic::Add(&sourcearr(cur_indx, engIndx), cur_coef*fluid_eng_src);
	  }
	  if (mass_trans) {
	    // Compute new particle diameter
	    Real new_dia = dia_part + 0.5*dt*d_dot;
	    if (new_dia < dia_eps) {
	      p.id() = -1; // Particle is considered completely evaporated
	      isub = nsub + 1; // Make sure to break out of while loop
	    } else {
	      p.rdata(pstateDia) = new_dia;
	    }
	  }
	  if (mom_trans) {
	    // Modify particle velocity by half time step
	    for (int dir = 0; dir != AMREX_SPACEDIM; ++dir) {
	      Gpu::Atomic::Add(&p.rdata(pstateVel+dir), 0.5*dt*part_mom_src[dir]/pmass);
	    }
	    // Modify particle position by whole time step
	    if (do_move) {
	      for (int dir = 0; dir != AMREX_SPACEDIM; ++dir) {
		Gpu::Atomic::Add(&p.pos(dir), dt*p.rdata(pstateVel+dir));
	      }
	    }
	  }
	  if (heat_trans) {
	    // Modify particle temperature
	    Gpu::Atomic::Add(&p.rdata(pstateT), 0.5*dt*part_temp_src);
	  }
	}
	// Increment sub-iteration
	++isub;
      } // End of isub while loop
    }); // End of loop over particles
  }
}
