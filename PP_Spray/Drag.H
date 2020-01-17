
#include <AMReX_REAL.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_AmrParticles.H>
#include <AMReX_Particles.H>
#include <PeleC_F.H>
#include <PeleC.H>
#include <Spray_F.H>
#include <chemistry_file.H>
//#include <eos.H>
#include <Fuego_EOS.H>
#include <SprayParticles.H>

using namespace amrex;


// Determine adjacent index locations and interpolation weightings
void
AdjIndexWeights(const IntVect&  indx,
		const RealVect& len,
		IntVect         indx_array[],
		Real            coef[],
		const int       indx_fact,
		const Real      wt_fact)
{
  // Compute distances for trilinear interpolation
  RealVect w_hi(AMREX_D_DECL(len[0] - indx[0], len[1] - indx[1], len[2] - indx[2]));
  RealVect w_lo = RealVect::Unit - w_hi;
  AMREX_D_TERM(indx_array[0] = indx;
	       indx_array[1] = indx + indx_fact*IntVect(AMREX_D_DECL(1,0,0));,
	       indx_array[2] = indx + indx_fact*IntVect(AMREX_D_DECL(0,1,0));
	       indx_array[3] = indx + indx_fact*IntVect(AMREX_D_DECL(1,1,0));,
	       indx_array[4] = indx + indx_fact*IntVect(0,0,1);
	       indx_array[5] = indx + indx_fact*IntVect(1,0,1);
	       indx_array[6] = indx + indx_fact*IntVect(0,1,1);
	       indx_array[7] = indx + indx_fact*IntVect(1,1,1););
  AMREX_D_TERM(coef[0] = AMREX_D_TERM(w_hi[0],*w_hi[1],*w_hi[2])*wt_fact;
	       coef[1] = AMREX_D_TERM(w_lo[0],*w_hi[1],*w_hi[2])*wt_fact;,
	       coef[2] = AMREX_D_TERM(w_hi[0],*w_lo[1],*w_hi[2])*wt_fact;
	       coef[3] = AMREX_D_TERM(w_lo[0],*w_lo[1],*w_hi[2])*wt_fact;,
	       coef[4] = w_hi[0]*w_hi[1]*w_lo[2]*wt_fact;
	       coef[5] = w_lo[0]*w_hi[1]*w_lo[2]*wt_fact;
	       coef[6] = w_hi[0]*w_lo[1]*w_lo[2]*wt_fact;
	       coef[7] = w_lo[0]*w_lo[1]*w_lo[2]*wt_fact;);
}

// Interpolate fluid state to particle locations
void
InterpStateAtParticle(Real&            rho_fluid,
		      Real&            T_fluid,
		      RealVect&        vel_fluid,
		      Vector<Real>&    Y_fluid,
		      IntVect          indx_array[],
		      Real             coef[],
		      const FArrayBox& state,
		      const IntVect&   indx,
		      const RealVect&  len,
		      const int&       nspecies,
		      EOS&             eos,
		      Vector<Real>&    mass_frac)
{
  const int num_adj = AMREX_D_PICK(2, 4, 8);
  AdjIndexWeights(indx, len, indx_array, coef, -1, 1.);
  
  // Loop over adjacent cell-centers
  for (int aindx = 0; aindx != num_adj; ++aindx)
    {
      IntVect cur_indx = indx_array[aindx];
      Real cur_coef = coef[aindx];
      Real cur_rho = state(cur_indx, PeleC::Density);
      Real inv_rho = 1./cur_rho;
      rho_fluid += cur_coef*cur_rho;
      Real ke = 0.;
      for (int dir = 0; dir != AMREX_SPACEDIM; ++dir)
	{
	  int nf = PeleC::Xmom + dir;
	  Real vel = state(cur_indx, nf)*inv_rho;
	  vel_fluid[dir] += cur_coef*vel;
	  ke += vel*vel;
	}
      ke *= 0.5;
      for (int ns = 0; ns != nspecies; ++ns)
	{
	  int mf_indx = ns + PeleC::FirstSpec;
	  Real cur_mf = state(cur_indx, mf_indx)*inv_rho;
	  Y_fluid[ns] += cur_coef*cur_mf;
	  mass_frac[ns] = cur_mf;
	}
      Real intEng = state(cur_indx, PeleC::Eden)*inv_rho - ke;
      Real T_val;
      eos.eos_EY2T(mass_frac.data(), intEng, T_val);
      T_fluid += cur_coef*T_val;
    }
}  

// Compute the inverse of the correction factors
Real
calcInvCorrectionFactor(const Real& B,
			const Real& logB)
{
  return B/(logB*std::pow(1. + B, 0.7));
}

// Calculate the fuel vapor mass fraction
Real
calcVaporMF(const Real& h_L,
	    const Real& T_part,
	    const Real& p_fluid,
	    const Real& mw_mix,
	    const Real& mw_fuel,
	    const Real& inv_boil_temp,
	    const Real& inv_Ru,
	    const Real& p_ref,
	    const Real& eps)
{
  Real pres_sat = p_ref*std::exp(h_L*inv_Ru*mw_fuel*
				 (inv_boil_temp - 1./T_part)) + eps;
  Real Yfv = mw_fuel*pres_sat/(mw_mix*p_fluid +
			       (mw_fuel - mw_mix)*pres_sat);
  return std::min(1. - eps, Yfv);
}

// Calculate evaporation rate and Sherwood number
void
calcSpecEvapRate(const Real& dia_part,
		 const Real& B_M_num,
		 const Real& Reyn,
		 const Real& powR,
		 const Real& Sc_skin,
		 const Real& D_skin,
		 Real&       Sh_num,
		 Real&       Y_dot)
{
  Real logB = std::log(1. + B_M_num);
  Real invFM = calcInvCorrectionFactor(B_M_num, logB);
  Real Sh_0 = 1. + powR*std::pow(1. + Reyn*Sc_skin, 1./3.);
  Sh_num = 2. + (Sh_0 - 2.)*invFM;
  Y_dot = -std::max(M_PI*D_skin*dia_part*Sh_num*logB, 0.);
}

// Compute the corrected Nusselt number and B_T value
Real
calcThermalB(const Real& ratio,
	     const Real& B_M,
	     const Real& B_eps,
	     Real&       Nu_num) // Nu_num should be initialized as uncorrected
{
  Real Nu_0 = Nu_num;
  Real phi = ratio/Nu_0;
  Real BM1 = 1. + B_M;
  Real B_T = std::pow(BM1, phi) - 1.;
  const int maxIter = 100;
  int k = 0;
  Real B_T_old, logB, invFT, Nu_tmp;
  Real error = B_eps*1000.;
  while (k < maxIter && error > B_eps)
    {
      B_T_old = B_T;
      logB = std::log(1. + B_T);
      invFT = calcInvCorrectionFactor(B_T, logB);
      Nu_tmp = 2. + (Nu_0 - 2.)*invFT;
      phi = ratio/Nu_tmp;
      B_T = std::pow(BM1, phi) - 1.;
      error = std::abs(B_T - B_T_old);
      k++;
    }
  B_T = std::min(std::max(B_T, 1.E-10), 20.);
  logB = std::log(1. + B_T);
  invFT = calcInvCorrectionFactor(B_T, logB);
  Nu_num = 2. + (Nu_0 - 2.)*invFT;
  return B_T;
}
      
template<class AoS>
void
update_particlesCPP(const int&      lev,
		    AoS&            particles,
		    FArrayBox&      state,
		    FArrayBox&      source,
		    const Box&      domain,
		    const RealBox&  prob_dom,
		    const IntVect&  reflect_lo,
		    const IntVect&  reflect_hi,
		    const RealVect& dx,
		    const Real&     flow_dt,
		    const int&      do_move)
{
  // Tolerance for Spalding number calculation
  const Real eps = 1.E-15;
  // Threshold for removing particle
  const Real dia_eps = 2.E-6;
  // B_T threshold
  const Real B_eps = 1.E-6;
  // Extract control parameters for mass, heat, and momentum transfer
  int heat_trans = PeleC::particle_heat_tran;
  int mass_trans = PeleC::particle_mass_tran;
  int mom_trans = PeleC::particle_mom_tran;
  // Number of particles
  const int Np = particles.size();
  const RealVect inv_dx = 1./dx;
  const Real third = 1./3.; // FIXME: should exist elsewhere
  const Real Pi_six = M_PI/6.;
  // Maximum number of sub-iterations
  const int nSubMax = 100;
  const Real vol = AMREX_D_TERM(dx[0],*dx[1],*dx[2]);
  Real inv_vol = 1./vol;
  // Get universal gas constant and reference pressure (1 atm)
  Real Ru, p0, ruc;
  CKRP(&Ru, &ruc, &p0);
  const Real inv_Ru = 1./Ru;

  const Box& state_box = state.box();
  const Box& src_box = source.box();

  // Low/hi positions of the problem domain
  const RealVect plo(AMREX_D_DECL(prob_dom.lo(0), prob_dom.lo(1), prob_dom.lo(2)));
  const RealVect phi(AMREX_D_DECL(prob_dom.hi(0), prob_dom.hi(1), prob_dom.hi(2)));

  EOS eos;
  // Retrieve the number of species for the gas phase
  int nspecies;
  get_num_spec(&nspecies);
  Vector<Real> mw_fluid(nspecies);
  get_mw(mw_fluid.data());

  // Get fuel property values
  int nspec_f = PeleC::particle_nfuel_species;
  
  // Compute the inverses of fuel properties
  Real ref_T = PeleC::partf_ref_T;
  Vector<Real> inv_boil_temp(nspec_f);
  Vector<Real> inv_diff_temp(nspec_f);
  Vector<Real> inv_fuel_mw(nspec_f);
  // Set initial C_p the same for all droplets
  // TODO: These should both be variables for each particle
  // TODO: This actually only needs to be computed once and stored in memory
  Vector<Real> ref_h(nspecies, 0.);
  eos.eos_T2HI(ref_T, ref_h.dataPtr());
  Real cp_L_av = 0.;
  Real ref_latent_heat = 0.;
  for (int L = 0; L != nspec_f; ++L)
    {
      const int fspec = PeleC::partf_indx[L];
      inv_fuel_mw[L] = 1./PeleC::partf_molwt[L];
      inv_boil_temp[L] = 1./PeleC::partf_boil_T[L];
      inv_diff_temp[L] = 1./(PeleC::partf_crit_T[L] - PeleC::partf_boil_T[L]);
      cp_L_av += PeleC::partf_mf[L]*PeleC::partf_cp[L];
      ref_latent_heat += PeleC::partf_mf[L]*PeleC::partf_latent[L];
    }

  // Other fuel properties to be filled in loops below
  Vector<Real> D_skin(nspec_f);
  Vector<Real> Sc_skin(nspec_f);
  Vector<Real> cp_f(nspec_f);
  Vector<Real> h_f(nspec_f);
  Vector<Real> Y_dot(nspec_f);
  Vector<Real> B_M_num(nspec_f); // Spalding number
  Vector<Real> B_T_num(nspec_f); // Spalding heat transfer number
  Vector<Real> Sh_num(nspec_f); // Sherwood number
  Vector<Real> L_fuel(nspec_f);

  RealVect vel_fluid(RealVect::Zero);
  Vector<Real> Y_fluid(nspecies, 0.);
  Vector<Real> vapor_h(nspecies, 0.);

  // Array of indices of adjacent cells
  IntVect indx_array[AMREX_D_PICK(2, 4, 8)];
  // Array of weightings for interpolation
  Real coef[AMREX_D_PICK(2, 4, 8)];

  // Make a dummy box that is sized as the number of particles
  Box dummyBox(IntVect::Zero, IntVect::Zero);
  // TODO: These should be using QVAR, cQRHO, cQU, etc but they are protected
  const int numQ = PeleC::NUM_STATE; // Should be PeleC::QVAR
  const int rhoIndx = PeleC::Density; // Should be PeleC::cQRHO
  const int velIndx = PeleC::Xmom;    // Should be PeleC::cQU
  const int tempIndx = PeleC::Temp;  // Should be PeleC::cQTEMP
  const int specIndx = PeleC::FirstSpec; // Should be PeleC::cQFS
  // Make FABs that contains the interpolated fluid data at
  // each particle location. The spatial index of the FAB
  // designates the corresponding particle index
  FArrayBox interp_state(dummyBox, numQ);

  // Component indices for transport properties
  int dComp_rhoD = 0; // Mass diffusion, rho*D_n
  int dComp_mu = dComp_rhoD + nspecies; // Dynamic viscosity
  int dComp_xi = dComp_mu + 1; // Volume viscosity
  int dComp_lambda = dComp_xi + 1; // Thermal conductivity
  int nCompTr = dComp_lambda + 1;
  // Interpolated transport properties at particle locations
  FArrayBox interp_trans(dummyBox, nCompTr);
  // Number of adjacent cells
  const int num_adj = AMREX_D_PICK(2, 4, 8);

  // Define one time for use later
  Vector<Real> mass_frac(nspecies, 0.);
  Vector<Real> cpi(nspecies, 0.);

  // Loop over number of particles
  for (int n = 0; n != Np; ++n)
    {
      // For referencing particle location in interp FABs
      IntVect part_IV(AMREX_D_DECL(0,0,0));
      // Compute the forcing at each particle location
      int isub = 1; // Initialize the number of sub-cycles
      int nsub = 1; // This is assigned in the first run through the loop
      Real dt = flow_dt;
      if (particles[n].id() == -1)
      	{
      	  // Particle does not exist, skip the remainder of the logic
      	  continue;
      	}
      // Since we are doing sub-iterations for the particles, we only want to
      // add a portion of the fluid source term back for each sub-iteration
      Real sub_source = inv_vol;
      Real rho_fluid = 0.;
      Real T_fluid = 0.;
      Real T_prev = 0.;
      while (isub <= nsub)
	{
	  // Particle data, temperature, diameter, diameter^2, density, and velocity
	  Real T_part = particles[n].rdata(PeleC::pstate_T);
	  Real dia_part = particles[n].rdata(PeleC::pstate_dia);
	  Real dia2_part = dia_part*dia_part;
	  Real rho_part = particles[n].rdata(PeleC::pstate_rho);
	  RealVect vel_part(AMREX_D_DECL(particles[n].rdata(PeleC::pstate_vel),
					 particles[n].rdata(PeleC::pstate_vel+1),
					 particles[n].rdata(PeleC::pstate_vel+2)));
	  Real part_ke = 0.5*std::sqrt(vel_part.radSquared());
	  // rho*D^3*Pi/6
	  Real pmass = Pi_six*rho_part*dia_part*dia2_part;
	  // Find the length in index space
	  RealVect len(AMREX_D_DECL((particles[n].pos(0) - plo[0])*inv_dx[0] + 0.5,
				    (particles[n].pos(1) - plo[1])*inv_dx[1] + 0.5,
				    (particles[n].pos(2) - plo[2])*inv_dx[2] + 0.5));

	  // Find the nearest cell center in the positive directions
	  IntVect indx(AMREX_D_DECL(std::floor(len[0]),
				    std::floor(len[1]),
				    std::floor(len[2])));
	  // Only interpolate fluid values at particle location once unless the particle has moved
	  if (isub == 1 || (do_move && mom_trans))
	    {
	      // Interpolated primitive state at the particle location
	      rho_fluid = 0.;
	      T_fluid = 0.;
	      vel_fluid = RealVect::Zero;
	      Y_fluid.assign(nspecies, 0.);
	      // Compute the interpolated fluid state of rho, T, u, and Y into
	      // rho_fluid, T_fluid, vel_fluid, and Y_fluid, respectively
	      InterpStateAtParticle(rho_fluid, T_fluid, vel_fluid, Y_fluid,
				    indx_array, coef, state, indx, len, nspecies,
				    eos, mass_frac);
	      // Plug interpolated data into FAB for computing transport coefs
	      interp_state(part_IV, rhoIndx) = rho_fluid;
	      AMREX_D_TERM(interp_state(part_IV, velIndx) = vel_fluid[0];,
			   interp_state(part_IV, velIndx+1) = vel_fluid[1];,
			   interp_state(part_IV, velIndx+2) = vel_fluid[2];);
	    }
	  Real delT = std::max(T_fluid - T_part, 0.);
	  // Compute the skin temperature of the droplet using the 1/3 rule
	  Real T_skin = T_part + third*delT;
	  interp_state(part_IV, tempIndx) = T_skin;
	  // Compute the enthalpy of the gas at the particle temperature
	  // TODO: This is inefficient to solve for the enthalpy of every
	  // species when we only need the enthalpy for the fuel species
	  // TODO: Make an eos function to do this
	  CKCPMS(&T_skin, cpi.dataPtr());
	  eos.eos_T2HI(T_part, vapor_h.dataPtr());
	  Real cp_vapor = 0.;
	  Real mw_mix = 0.;
	  for (int ns = 0; ns != nspecies; ++ns)
	    {
	      int mf_indx = ns + specIndx;
	      interp_state(part_IV, mf_indx) = Y_fluid[ns];
	      Real ymw = Y_fluid[ns]/mw_fluid[ns];
	      mw_mix += ymw;
	      cp_vapor += Y_fluid[ns]*cpi[ns];
	    }
	  Real R_fluid = Ru*mw_mix;
	  Real p_fluid = rho_fluid*R_fluid*T_fluid;
	  mw_mix = 1./mw_mix;
	  Real Y_fuel = 0.;
	  for (int L = 0; L != nspec_f; ++L)
	    {
	      const int fspec = PeleC::partf_indx[L];
	      const int mf_indx = fspec + specIndx;
	      Real part_latent_heat = vapor_h[fspec] - ref_h[fspec]
		+ ref_latent_heat - cp_L_av*(T_part - ref_T);
// 	      Real T_crit = PeleC::partf_crit_T[L];
// 	      Real T_boil = PeleC::partf_boil_T[L];
// 	      Real part_latent_heat = ref_latent_heat*
// 		std::pow((T_crit - T_part)*inv_diff_temp[L], 0.38);
	      L_fuel[L] = part_latent_heat;
	      Real Yfv = calcVaporMF(part_latent_heat, T_part, p_fluid,
				     mw_mix, PeleC::partf_molwt[L],
				     inv_boil_temp[L], inv_Ru, p0, eps);
	      B_M_num[L] = std::min(20., std::max((Yfv - Y_fluid[fspec])/(1. - Yfv), 0.));
	      cp_f[L] = cpi[fspec];
	      h_f[L] = vapor_h[fspec];
	      Y_fuel += Y_fluid[fspec];
	    }
	  bool all_fuel_vapor = false;
	  if (Y_fuel >= 1.)
	    {
	      all_fuel_vapor = true;
	    }
	  Real cp_skin = cp_vapor;
	  // Compute the interpolated transport properties
	  get_transport_coeffs(ARLIM_3D(dummyBox.loVect()),
			       ARLIM_3D(dummyBox.hiVect()),
			       BL_TO_FORTRAN_N_3D(interp_state, specIndx),
			       BL_TO_FORTRAN_N_3D(interp_state, tempIndx),
			       BL_TO_FORTRAN_N_3D(interp_state, rhoIndx),
			       BL_TO_FORTRAN_N_3D(interp_trans, dComp_rhoD),
			       BL_TO_FORTRAN_N_3D(interp_trans, dComp_mu),
			       BL_TO_FORTRAN_N_3D(interp_trans, dComp_xi),
			       BL_TO_FORTRAN_N_3D(interp_trans, dComp_lambda));
	  // Compute mass transport coefficient at skin temperature and reset vectors
	  for (int L = 0; L != nspec_f; ++L)
	    {
	      const int fspec = dComp_rhoD + PeleC::partf_indx[L];
	      D_skin[L] = interp_trans(part_IV, fspec);
	      Y_dot[L] = 0.;
	    }
	  Real mu_skin = interp_trans(part_IV, dComp_mu);
	  Real xi_skin = interp_trans(part_IV, dComp_xi);
	  Real lambda_skin = interp_trans(part_IV, dComp_lambda);
	  RealVect diff_vel = vel_fluid - vel_part;
	  Real diff_vel_mag = diff_vel.vectorLength();
	  // Local Reynolds number
	  Real Reyn = rho_fluid*diff_vel_mag*dia_part/mu_skin;
	  Real tau_var = mu_skin/(rho_part*dia2_part);
	  Real inv_tau_v = 18.*tau_var;
	  if (isub == 1)
	    {
	      nsub = std::min(int(flow_dt*inv_tau_v) + 1, nSubMax);
	      dt = flow_dt/Real(nsub);
	    }
	  Real m_dot = 0.; // Total mass transfer (g/s)
	  Real d_dot = 0.; // Diameter rate of change
	  // Calculate skin Prandtl number
	  Real Pr_skin = std::min(mu_skin*cp_skin/lambda_skin, 3.);
	  // Calculate the uncorrected Nusselt number
	  Real powR = std::max(std::pow(Reyn, 0.077), 1.);
	  // Uncorrected Nusselt Number
	  Real Nu_0 = 1. + powR*std::pow(1. + Reyn*Pr_skin, third);
	  // If mass transfer is applied 
	  if (mass_trans && !all_fuel_vapor)
	    {
	      for (int L = 0; L != nspec_f; ++L)
		{
		  const int fspec = PeleC::partf_indx[L];
		  // Compute the skin Schmidt number
		  Sc_skin[L] = std::min(mu_skin/D_skin[L], 3.);
		  // Calculate the Sherwood number and the evaporation rate
		  calcSpecEvapRate(dia_part, B_M_num[L], Reyn, powR, Sc_skin[L],
				   D_skin[L], Sh_num[L], Y_dot[L]);
		  m_dot += Y_dot[L];
		}
	      // Diameter rate time scale from Borghesi et al. 2012
	      if (isub == 1)
		{
		  Real inv_tau_p = std::abs(m_dot)/pmass;
		  nsub = std::min(std::max(nsub,int(flow_dt*inv_tau_p) + 1), nSubMax);
		  dt = flow_dt/Real(nsub);
		}
	      d_dot = m_dot/(0.5*M_PI*rho_part*dia2_part);
	      if (d_dot != d_dot)
		{
		  Print() << "Error in d_dot for particle ID " << n << std::endl;
		  Abort();
		}
	    }
	  // Fluid momentum source term
	  RealVect fluid_mom_src = RealVect::Zero;
	  // Particle momentum source term
	  RealVect part_mom_src = RealVect::Zero;
	  // Energy source term for the fluid
	  Real fluid_eng_src = 0.;
	  // Temperature source term for the particle
	  Real part_temp_src = 0.;
	  if (mom_trans)
	    {
	      Real drag_coef = 1. + 0.15*std::pow(Reyn, 0.687);
	      Real drag_force = 0.125*rho_fluid*drag_coef*M_PI*dia2_part*diff_vel_mag;
	      part_mom_src = diff_vel*drag_force;
	      fluid_mom_src = part_mom_src + vel_part*m_dot;
	      // s_d,mu dot u_d 
	      Real S_dmu_dot_u = part_mom_src.dotProduct(vel_part);
	      fluid_eng_src += S_dmu_dot_u + m_dot*part_ke;
	    }
	  if ((heat_trans || mass_trans) && !all_fuel_vapor)
	    {
	      // Calculate energy transfer due to heat convection and evaporation
	      Real coeff_heat = 0.;
	      for (int L = 0; L != nspec_f; ++L)
		{
		  Real Nu_num = Nu_0;
		  // Calculate Spalding heat transfer number and the
		  // corrected Nusselt number
		  Real ratio = cp_f[L]*Pr_skin*Sh_num[L]/(cp_skin*Sc_skin[L]);
		  B_T_num[L] = calcThermalB(ratio, B_M_num[L], B_eps, Nu_num);
		  // Convection term
		  coeff_heat += Nu_num*std::log(1. + B_T_num[L])/B_T_num[L];
		  Real S_dmh = M_PI*lambda_skin*dia_part*Nu_num*delT*
		    std::log(1. + B_T_num[L])/B_T_num[L];
		  fluid_eng_src += Y_dot[L]*h_f[L];
		  part_temp_src += Y_dot[L]*L_fuel[L];
		}
	      Real conv_src = M_PI*lambda_skin*dia_part*delT*coeff_heat;
	      fluid_eng_src += conv_src;
	      part_temp_src = (part_temp_src + conv_src)/(pmass*cp_L_av);
	      if (isub == 1)
		{
		  Real inv_tau_T = conv_src/(pmass*cp_L_av*delT);
		  nsub = std::min(std::max(nsub, int(flow_dt*inv_tau_T) + 1), nSubMax);
		  dt = flow_dt/Real(nsub);
		}
	    }
	  // Multiply the fluid source terms by a percentage of the sub-iteration
	  if (isub == 1)
	    {
	      sub_source /= Real(nsub);
	    }
	  // Put the same forcing term on the grid
	  if (mom_trans || mass_trans || heat_trans)
	    {
	      RealVect len2 = len - RealVect::Unit;
	      IntVect indx2(AMREX_D_DECL(std::floor(len2[0]),
					 std::floor(len2[1]),
					 std::floor(len2[2])));
	      AdjIndexWeights(indx2, len2, indx_array, coef, 1, sub_source);
	      for (int aindx = 0; aindx != num_adj; ++aindx)
		{
		  Real cur_coef = coef[aindx];
		  IntVect cur_indx = indx_array[aindx];
		  BL_ASSERT(src_box.contains(cur_indx));
		  for (int dir = 0; dir != AMREX_SPACEDIM; ++dir)
		    {
		      const int nf = PeleC::Xmom + dir;
		      source(cur_indx, nf) -= cur_coef*fluid_mom_src[dir];
		    }
		  source(cur_indx, PeleC::Density) -= cur_coef*m_dot;
		  for (int L = 0; L != nspec_f; ++L)
		    {
		      const int nf = PeleC::FirstSpec + PeleC::partf_indx[L];
		      source(cur_indx, nf) -= cur_coef*Y_dot[L];
		    }
		  source(cur_indx, PeleC::Eden) -= cur_coef*fluid_eng_src;
		}
	      // Compute new particle diameter
	      Real new_dia = particles[n].rdata(PeleC::pstate_dia) + 0.5*dt*d_dot;
	      if (new_dia < dia_eps)
		{
		  // Particle is considered completely evaporated
		  particles[n].id() = -1;
		}
	      else
		{
		  particles[n].rdata(PeleC::pstate_dia) = new_dia;
		}
	      // Modify particle velocity and postion
	      for (int dir = 0; dir != AMREX_SPACEDIM; ++dir)
		{
		  particles[n].rdata(PeleC::pstate_vel+dir) += 
		    0.5*dt*part_mom_src[dir]/pmass;
		}
	      if (do_move)
		{
		  for (int dir = 0; dir != AMREX_SPACEDIM; ++dir)
		    {
		      particles[n].pos(dir) += dt*
			particles[n].rdata(PeleC::pstate_vel+dir);
		    }
		}
	      // Modify particle temperature
	      particles[n].rdata(PeleC::pstate_T) += 0.5*dt*part_temp_src;
	    }
	  // TODO: Confirm this is how particles are eliminated
	  if (particles[n].id() == -1)
	    {
	      // Particle does not exist, break out of while loop
	      isub = nsub + 1;
	    }
	  isub += 1;
	} // End of sub-cycles while loop
    } // End of loop over particles
  if (do_move)
    {
      for (int dir = 0; dir != AMREX_SPACEDIM; ++dir)
	{
	  if (reflect_lo[dir])
	    {
	      // Loop over all particles
	      for (int n = 0; n != Np; ++n)
		{
		  if (particles[n].id() != -1 && 
		      particles[n].pos(dir) < plo[dir])
		    {
		      particles[n].pos(dir) *= -1.;
		      particles[n].pos(dir) += 2.*plo[dir];
		      particles[n].rdata(PeleC::pstate_vel + dir) *= -1.;
		    }
		}
	    }
	  if (reflect_hi[dir])
	    {
	      for (int n = 0; n != Np; ++n)
		{
		  if (particles[n].id() != -1 && 
		      particles[n].pos(dir) > phi[dir])
		    {
		      particles[n].pos(dir) *= -1.;
		      particles[n].pos(dir) += 2.*phi[dir];
		      particles[n].rdata(PeleC::pstate_vel + dir) *= -1.;
		    }
		}
	    }
	}
    }
  for (int dir = 0; dir != AMREX_SPACEDIM; ++dir)
    {
      if (reflect_lo[dir])
	{
	  Box extDomBox = adjCellLo(domain, dir);
	  // Extract box that is outside the domain but part of src_box
	  Box extCells = src_box & extDomBox;
	  if (!extCells.isEmpty())
	    {
	      // Extract cells adjacent to the domain edges
	      Box intCells = extCells;
	      intCells.shift(dir, 1);
	      source.atomicAdd(source, extCells, intCells,
			       0, 0, PeleC::NUM_STATE);
	    }
	}
      if (reflect_hi[dir])
	{
	  Box extDomBox = adjCellHi(domain, dir);
	  // Extract box that is outside the domain but part of src_box
	  Box extCells = src_box & extDomBox;
	  if (!extCells.isEmpty())
	    {
	      // Extract cells adjacent to the domain edges
	      Box intCells = extCells;
	      intCells.shift(dir, -1);
	      source.atomicAdd(source, extCells, intCells,
			       0, 0, PeleC::NUM_STATE);
	    }
	}
    }
}
