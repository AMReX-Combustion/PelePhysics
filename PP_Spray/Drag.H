#ifndef DRAG_H_
#define DRAG_H_

#include <AMReX_REAL.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_AmrParticles.H>
#include <AMReX_Particles.H>
#include <AMReX_ParIter.H>
#include "SprayParticles.H"

using namespace amrex;

// Compute the heat transfer coefficient using the
// corrected Nusselt number and B_T value
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
Real calcHeatCoeff (const Real& ratio,
                    const Real& B_M,
                    const Real& B_eps,
                    const Real& C_eps,
                    const Real& Nu_0)
{
  if (B_M <= C_eps) return 0.;
  const int maxIter = 100;
  const Real NU2 = Nu_0 - 2.;
  const Real BM1 = 1. + B_M;
  Real phi = ratio/Nu_0;
  Real B_T_old = std::pow(BM1, phi) - 1.;
  Real BT1 = 1. + B_T_old;
  Real logB = std::log(BT1);
  Real invFT = B_T_old/(logB*std::pow(BT1, 0.7));
  Real Nu_num = 2. + NU2*invFT;
  phi = ratio/Nu_num;
  Real B_T = std::pow(BM1, phi) - 1.;
  Real error = std::abs(B_T - B_T_old);
  int k = 0;
  while (k < maxIter && error > B_eps) {
    B_T_old = B_T;
    BT1 = 1. + B_T;
    logB = std::log(BT1);
    invFT = B_T/(logB*std::pow(BT1, 0.7));
    Nu_num = 2. + NU2*invFT;
    phi = ratio/Nu_num;
    B_T = std::pow(BM1, phi) - 1.;
    error = std::abs(B_T - B_T_old);
    k++;
  }
  logB = std::log(1. + B_T);
  invFT = B_T/(logB*std::pow(1. + B_T, 0.7));
  Nu_num = 2. + NU2*invFT;
  return Nu_num*logB/B_T;
}

// Compute source terms and update particles
// Return true if particle should be removed
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
bool calculateSpraySource (const Real    flow_dt,
                           const bool    do_move,
                           GasPhaseVals& gpv,
                           SprayComps    SPI,
                           SprayData     fdat,
                           SprayParticleContainer::ParticleType& p,
#ifdef USE_SPRAY_SOA
                           const std::array<SprayParticleContainer::SoA::RealVector, NAR_SPR>& attribs,
#endif
                           const int pid,
                           TransParm const* trans_parm)
{
  const Real rule = 1./3.;
  const Real inv_Ru = 1./EOS::RU;
  const Real C_eps = 1.E-15;
  const Real B_eps = 1.E-7;
  const Real dia_eps = 2.E-6;
  const int nSubMax = 100;
  bool get_xi = false;
  bool get_Ddiag = true;
  bool get_lambda = true;
  bool get_mu = true;
  if (!SPI.mass_tran && !SPI.heat_tran) {
    get_Ddiag = false;
    get_lambda = false;
  }
  GpuArray<Real, NUM_SPECIES> Y_skin;
  GpuArray<Real, NUM_SPECIES> h_skin;
  GpuArray<Real, NUM_SPECIES> cp_n;
  GpuArray<Real, NUM_SPECIES> Ddiag;
  GpuArray<Real, SPRAY_FUEL_NUM> B_M_num;
  GpuArray<Real, SPRAY_FUEL_NUM> Sh_num;
  GpuArray<Real, SPRAY_FUEL_NUM> L_fuel;
  GpuArray<Real, SPRAY_FUEL_NUM> mi_dot;
  GpuArray<Real, SPRAY_FUEL_NUM> Y_part;
#ifdef USE_SPRAY_SOA
  RealVect vel_part(AMREX_D_DECL(attribs[SPI.pstateVel].data()[pid],
                                 attribs[SPI.pstateVel+1].data()[pid],
                                 attribs[SPI.pstateVel+2].data()[pid]));
  Real T_part = attribs[SPI.pstateT].data()[pid];
  Real dia_part = attribs[SPI.pstateDia].data()[pid];
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf)
    Y_part[spf] = attribs[SPI.pstateY+spf].data()[pid];
#else
  RealVect vel_part(AMREX_D_DECL(p.rdata(SPI.pstateVel),
                                 p.rdata(SPI.pstateVel+1),
                                 p.rdata(SPI.pstateVel+2)));
  Real T_part = p.rdata(SPI.pstateT);
  Real dia_part = p.rdata(SPI.pstateDia);
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf)
    Y_part[spf] = p.rdata(SPI.pstateY+spf);
#endif
  Real dt = flow_dt;
  int isub = 1;
  int nsub = 1;
  while (isub <= nsub) {
    Real rho_part = 0.;
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf)
      rho_part += Y_part[spf]*fdat.rho[spf];
    Real dia2_part = dia_part*dia_part;
    Real pmass = M_PI/6.*rho_part*dia_part*dia2_part;
    Real part_ke = 0.5*vel_part.radSquared();
    // Model the fuel vapor using the one-third rule
    Real delT = amrex::max(gpv.T_fluid - T_part, 0.);
    Real T_skin = T_part + rule*delT;
    // Calculate the C_p at the skin temperature for each species
    EOS::T2Cpi(T_skin, cp_n.data());
    EOS::T2Hi(T_part, h_skin.data());
    if (SPI.heat_tran || SPI.mass_tran) {
      for (int n = 0; n < NUM_SPECIES; ++n)
        Y_skin[n] = 0.;
    } else {
      for (int n = 0; n < NUM_SPECIES; ++n)
        Y_skin[n] = gpv.Y_fluid[n];
    }
    // Solve for state of the vapor and mass transfer coefficient B_M
    Real sumYSkin = 0.; // Mass fraction of the fuel in skin film, uses one-thirds rule
    Real sumYFuel = 0.; // Mass fraction of the fuel in the gas phase
    Real cp_skin = 0.; // Averaged C_p at particle surface
    Real cp_L_av = 0.; // Cp of the liquid state
    Real mw_vap = 0.; // Average molar mass of vapor mixture
    if (SPI.heat_tran || SPI.mass_tran) {
      GpuArray<Real, SPRAY_FUEL_NUM> x_v;
      Real sum_xv = 0.; // Sum x_v
      Real sum_mw_xv = 0.; // Sum mw_v*x_v
      Real nt = 0.; // Total moles of liquid
      // This computes the mass fractions in the vapor
      // The equation is a normalized version of y_v = x_vi*mw_i/((1 - sum_k x_vk) + sum_j x_vj mw_j)
      // x_vi = x_li*p_sat/(p_gas*mw_gas)
      // x_li = (y_li/mw_i)/(sum_j y_lj/mw_j)
      for (int spf = 0; spf != SPRAY_FUEL_NUM; ++spf) {
        const int fspec = fdat.indx[spf];
        const Real mw_fuel = gpv.mw_fluid[fspec];
        // Compute latent heat
        Real part_latent = h_skin[fspec] + fdat.latent[spf]
          - fdat.cp[spf]*(T_part - fdat.ref_T);
        L_fuel[spf] = part_latent;
        // Compute the mass fraction of the fuel vapor at droplet surface
        Real pres_sat = EOS::PATM*std::exp(part_latent*inv_Ru*mw_fuel*
                                           (1./fdat.boilT[spf] - 1./T_part)) + C_eps;
        Real x_l = Y_part[spf]/mw_fuel;
        Real x_vc = x_l*pres_sat;
        nt += x_l;
        sum_xv += x_vc;
        x_v[spf] = mw_fuel*x_vc;
        sum_mw_xv += x_v[spf];
      }
      Real totalmwx = gpv.mw_mix*(nt*gpv.p_fluid - sum_xv) + sum_mw_xv;
      for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
        const int fspec = fdat.indx[spf];
        Real Yfv = x_v[spf]/totalmwx;
        Yfv = amrex::max(0., amrex::min(1. - C_eps, Yfv));
        B_M_num[spf] = amrex::max(C_eps, (Yfv - gpv.Y_fluid[fspec])/(1. - Yfv));
        Y_skin[fspec] = Yfv + rule*(gpv.Y_fluid[fspec] - Yfv);
        sumYSkin += Y_skin[fspec];
        cp_L_av += Y_part[spf]*fdat.cp[spf];
        sumYFuel += gpv.Y_fluid[fspec];
      }
      const Real restYSkin = 1. - sumYSkin;
      for (int sp = 0; sp != NUM_SPECIES; ++sp) {
        Y_skin[sp] += restYSkin*gpv.Y_fluid[sp];
        cp_skin += Y_skin[sp]*cp_n[sp];
        mw_vap += Y_skin[sp]*gpv.invmw[sp];
      }
      mw_vap = 1./mw_vap;
    } else {
      mw_vap = gpv.mw_mix;
    }
    Real lambda_skin = 0.;
    Real mu_skin = 0.;
    Real xi_skin = 0.;
    transport(get_xi, get_mu, get_lambda, get_Ddiag, T_skin,
              gpv.rho_fluid, Y_skin.data(), Ddiag.data(), mu_skin,
              xi_skin, lambda_skin, trans_parm);
    // Ensure gas is not all fuel to allow evaporation
    bool evap_fuel = (sumYFuel >= 1.) ? false : true;
    RealVect diff_vel = gpv.vel_fluid - vel_part;
    Real diff_vel_mag = diff_vel.vectorLength();
    // Local Reynolds number
    Real Reyn = gpv.rho_fluid*diff_vel_mag*dia_part/mu_skin;
    Real Nu_0 = 1.;
    // Solve mass transfer source terms
    Real m_dot = 0.;
    Real d_dot = 0.;
    if ((SPI.mass_tran || SPI.heat_tran) && evap_fuel) {
      Real Pr_skin = mu_skin*cp_skin/lambda_skin;
      Real powR = amrex::max(std::pow(Reyn, 0.077), 1.);
      Nu_0 = 1. + powR*std::cbrt(1. + Reyn*Pr_skin);
      for (int spf = 0; spf != SPRAY_FUEL_NUM; ++spf) {
        const int fspec = fdat.indx[spf];
        // Convert mass diffusion coefficient from mixture average
        // to binary for fuel only, not concerned with other species
        Ddiag[fspec] *= mw_vap*gpv.invmw[fspec];
        const Real rhoD = Ddiag[fspec];
        const Real Sc_skin = mu_skin/rhoD;
        const Real B_M = B_M_num[spf];
        Real logB = std::log(1. + B_M);
        // Calculate Sherwood number and evaporation rate
        Real invFM = B_M/(logB*std::pow(1. + B_M, 0.7));
        Real Sh_0 = 1. + powR*std::cbrt(1. + Reyn*Sc_skin);
        Sh_num[spf] = 2. + (Sh_0 - 2.)*invFM;
        if (SPI.mass_tran) {
          mi_dot[spf] = -amrex::max(M_PI*rhoD*dia_part*Sh_num[spf]*logB, 0.);
          gpv.fluid_Y_dot[spf] += mi_dot[spf];
          m_dot += mi_dot[spf];
        } else {
          mi_dot[spf] = 0.;
        }
      }
      d_dot = m_dot/(0.5*M_PI*rho_part*dia2_part);
      Real inv_tau_d = -m_dot/(3.*pmass);
      if (isub == 1)
        nsub = amrex::min(int(flow_dt*inv_tau_d) + 1, nSubMax);
    }
    // Solve momentum source terms
    const Real inv_pmass = 1./pmass;
    RealVect part_mom_src(RealVect::TheZeroVector());
    if (SPI.mom_tran) {
      Real drag_coef = 0.;
      if (Reyn > 0.)
        drag_coef =
          (Reyn > 1.) ? 24./Reyn*(1. + std::cbrt(Reyn*Reyn)/6.) : 24./Reyn;
      Real drag_force = 0.125*gpv.rho_fluid*drag_coef*M_PI*dia2_part*diff_vel_mag;
      part_mom_src = drag_force*diff_vel;
      gpv.fluid_mom_src += part_mom_src + vel_part*m_dot;
      // s_d,mu dot u_d
      Real S_dmu_dot_u = part_mom_src.dotProduct(vel_part);
#ifndef SPRAY_PELE_LM
      gpv.fluid_eng_src += S_dmu_dot_u + m_dot*part_ke;
#endif
      Real inv_tau_var = drag_force*inv_pmass;
      if (isub == 1)
        nsub = amrex::min(amrex::max(nsub, int(flow_dt*inv_tau_var) + 1), nSubMax);
    }
    // Solve for energy source terms
    Real part_temp_src = 0.;
    if (evap_fuel && (SPI.heat_tran || SPI.mass_tran)) {
      const Real inv_pm_cp = inv_pmass/cp_L_av;
      Real coeff_heat = 0.;
      for (int spf = 0; spf != SPRAY_FUEL_NUM; ++spf) {
        const int fspec = fdat.indx[spf];
        Real ratio = cp_n[fspec]*Sh_num[spf]*Ddiag[fspec]/lambda_skin;
        Real heatC = calcHeatCoeff(ratio, B_M_num[spf], B_eps, C_eps, Nu_0);
        // Convection term
        coeff_heat += heatC;
        gpv.fluid_eng_src += mi_dot[spf]*h_skin[fspec];
        part_temp_src += mi_dot[spf]*L_fuel[spf];
      }
      Real conv_src = M_PI*lambda_skin*dia_part*delT*coeff_heat;
      if (SPI.heat_tran) {
        gpv.fluid_eng_src += conv_src;
        part_temp_src += conv_src;
      }
      part_temp_src *= inv_pm_cp;
      if (isub == 1 && delT > C_eps) {
        Real inv_tau_T = conv_src*inv_pm_cp/delT;
        nsub = amrex::min(amrex::max(nsub, int(flow_dt*inv_tau_T) + 1), nSubMax);
      }
    }
    if (isub == 1)
      dt = flow_dt/Real(nsub);
    const Real part_dt = 0.5*dt;
    AMREX_D_TERM(vel_part[0] += part_dt*part_mom_src[0]*inv_pmass;,
                 vel_part[1] += part_dt*part_mom_src[1]*inv_pmass;,
                 vel_part[2] += part_dt*part_mom_src[2]*inv_pmass;);
    T_part += part_dt*part_temp_src;
    Real new_dia = dia_part + part_dt*d_dot;
    if (new_dia > dia_eps) {
      dia_part = new_dia;
      Real new_mass = pmass + m_dot*part_dt;
      for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf)
        Y_part[spf] = (Y_part[spf]*pmass + mi_dot[spf]*part_dt)/new_mass;
    } else {
      for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf)
        mi_dot[spf] += Y_part[spf]*pmass;
      m_dot += pmass;
      p.id() = -1;
      isub = nsub + 1;
    }
    gpv.fluid_mass_src += m_dot;
    ++isub;
  }
#ifdef USE_SPRAY_SOA
  AMREX_D_TERM(attribs[SPI.pstateVel].data()[pid] = vel_part[0];,
               attribs[SPI.pstateVel+1].data()[pid] = vel_part[1];,
               attribs[SPI.pstateVel+2].data()[pid] = vel_part[2];);
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf)
    attribs[SPI.pstateY+spf].data()[pid] = Y_part[spf];
  attribs[SPI.pstateT].data()[pid] = T_part;
  attribs[SPI.pstateDia].data()[pid] = dia_part;
#else
  AMREX_D_TERM(p.rdata(SPI.pstateVel) = vel_part[0];,
               p.rdata(SPI.pstateVel+1) = vel_part[1];,
               p.rdata(SPI.pstateVel+2) = vel_part[2];);
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf)
    p.rdata(SPI.pstateY+spf) = Y_part[spf];
  p.rdata(SPI.pstateT) = T_part;
  p.rdata(SPI.pstateDia) = dia_part;
#endif
  return false;
}

#endif
