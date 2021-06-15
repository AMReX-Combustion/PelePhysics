#ifndef DRAG_H_
#define DRAG_H_

#include "SprayParticles.H"
#include <AMReX_AmrParticles.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParIter.H>
#include <AMReX_Particles.H>
#include <AMReX_REAL.H>
#include <Transport.H>

using namespace amrex;

// Compute the heat transfer coefficient using the
// corrected Nusselt number and B_T value
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
Real
calcHeatCoeff(
  const Real& ratio,
  const Real& B_M,
  const Real& B_eps,
  const Real& C_eps,
  const Real& Nu_0)
{
  if (B_M <= C_eps)
    return 0.;
  const int maxIter = 100;
  const Real NU2 = Nu_0 - 2.;
  const Real BM1 = 1. + B_M;
  Real phi = ratio / Nu_0;
  Real B_T_old = std::pow(BM1, phi) - 1.;
  Real BT1 = 1. + B_T_old;
  Real logB = std::log(BT1);
  Real invFT = B_T_old / (logB * std::pow(BT1, 0.7));
  Real Nu_num = 2. + NU2 * invFT;
  phi = ratio / Nu_num;
  Real B_T = std::pow(BM1, phi) - 1.;
  Real error = std::abs(B_T - B_T_old);
  int k = 0;
  while (k < maxIter && error > B_eps) {
    B_T_old = B_T;
    BT1 = 1. + B_T;
    logB = std::log(BT1);
    invFT = B_T / (logB * std::pow(BT1, 0.7));
    Nu_num = 2. + NU2 * invFT;
    phi = ratio / Nu_num;
    B_T = std::pow(BM1, phi) - 1.;
    error = std::abs(B_T - B_T_old);
    k++;
  }
  logB = std::log(1. + B_T);
  invFT = B_T / (logB * std::pow(1. + B_T, 0.7));
  Nu_num = 2. + NU2 * invFT;
  return Nu_num * logB / B_T;
}

// Compute the flash boiling coefficient
// From Adachi 1997
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
Real
calcAlpha(const Real& delTb)
{
  SprayUnits SPU;
  Real alpha = 0.;
  if (delTb > 25.) {
    alpha = 13800.0 * std::pow(delTb, 0.39);
  } else if (delTb > 5.) {
    alpha = 27.0 * std::pow(delTb, 2.33);
  } else {
    alpha = 760.0 * std::pow(delTb, 0.26);
  }
  // Normally in J/(s m^2 K), convert when necessary
  return alpha * 1.E3 * SPU.mass_conv;
}

// Compute the flash boiling vaporization rate
// From Zuo, Gomes, and Rutland 2000
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
Real
calcFlashVaporRate(const Real& dh, const Real& coeff, const Real& Gf)
{
  if (dh < 0.)
    return 0.;
  Real tol = 1.E-4;
  Real Grat = 1.E-5;
  Real G = 0.;
  Real error = 100.;
  while (error > tol) {
    Real Gold = Grat;
    G = coeff / (1. + Grat) * std::log(1. + (1. + Grat) * dh);
    Grat = Gf / G;
    error = std::abs(Grat - Gold) / Gold;
  }
  return G;
}

// Estimate the boil temperature
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
calcBoilT(const SprayData& fdat, const GasPhaseVals& gpv, Real* cBoilT)
{
  SprayUnits SPU;
  Real RU = pele::physics::Constants::RU * SPU.ru_conv;
  Real PATM = pele::physics::Constants::PATM * SPU.pres_conv;
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    const int fspec = fdat.indx[spf];
    const Real critT = fdat.critT[spf];
    const Real boilT_ref = fdat.boilT[spf];
    const Real mw_fuel = gpv.mw_fluid[fspec];
    // Since we only know the latent heat at the reference temperature,
    // modify Watsons power law to find latent heat at boiling conditions
    Real Hboil_ref =
      fdat.ref_latent[spf] *
      std::pow((critT - fdat.ref_T) / (critT - boilT_ref), -0.38);
    // Estimate the boiling temperature at the gas phase pressure using
    // Clasius-Clapeyron relation
    cBoilT[spf] =
      1. / (std::log(PATM / gpv.p_fluid) * RU / (Hboil_ref * mw_fuel) +
            1. / boilT_ref);
    cBoilT[spf] = amrex::min(critT, cBoilT[spf]);
  }
}

// Compute the vapor mass fractions and
// latent heat
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
calcVaporY(
  const SprayData& fdat,
  const GasPhaseVals& gpv,
  const Real& T_in,
  const Real& C_eps,
  const Real* Y_l,
  const Real* h_part,
  const Real* cBoilT,
  Real* Y_vapor,
  Real* L_fuel)
{
  SprayUnits SPU;
  Real RU = pele::physics::Constants::RU * SPU.ru_conv;
  Real PATM = pele::physics::Constants::PATM * SPU.pres_conv;
  Real sum_xv = 0.;    // Sum x_v
  Real sum_mw_xv = 0.; // Sum mw_v*x_v
  Real nt = 0.;        // Total moles of liquid
  // This computes the mass fractions in the vapor
  // The equation is a normalized version of y_v = x_vi*mw_i/((1 - sum_k x_vk) +
  // sum_j x_vj mw_j)
  // x_vi = x_li*p_sat/(p_gas*mw_gas)
  // x_li = (y_li/mw_i)/(sum_j y_lj/mw_j)
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    const int fspec = fdat.indx[spf];
    const Real boilT_ref = fdat.boilT[spf];
    const Real mw_fuel = gpv.mw_fluid[fspec];
    Real T_part = amrex::min(T_in, cBoilT[spf]);
    // Compute latent heat
    Real part_latent =
      h_part[fspec] + fdat.latent[spf] - fdat.cp[spf] * (T_part - fdat.ref_T);
    L_fuel[spf] = part_latent;
    // Compute the saturation pressure at the current particle temperature
    Real a = fdat.psat_coef[4 * spf];
    Real b = fdat.psat_coef[4 * spf + 1];
    Real c = fdat.psat_coef[4 * spf + 2];
    Real d = fdat.psat_coef[4 * spf + 3]; // For converting to Pa or dyne/cm^2
    Real pres_sat = 0.;
    // Using the Clasius-Clapeyron relation
    if (d == 0.) {
      pres_sat =
        PATM *
        std::exp(part_latent * mw_fuel / RU * (1. / boilT_ref - 1. / T_part));
      // Using the Antoine equation
    } else {
      pres_sat = d * std::pow(10., a - b / (T_part + c));
    }
    Real x_l = Y_l[spf] / mw_fuel;
    Real x_vc = x_l * pres_sat;
    nt += x_l;
    sum_xv += x_vc;
    Y_vapor[spf] = mw_fuel * x_vc;
    sum_mw_xv += Y_vapor[spf];
  }
  Real totalmwx = gpv.mw_mix * (nt * gpv.p_fluid - sum_xv) + sum_mw_xv;
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    Real Yfv = Y_vapor[spf] / totalmwx;
    Y_vapor[spf] = amrex::max(0., amrex::min(1. - C_eps, Yfv));
  }
}

// Compute source terms and update particles
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
calculateSpraySource(
  const Real flow_dt,
  const bool do_move,
  GasPhaseVals& gpv,
  SprayComps SPI,
  SprayData fdat,
  SprayParticleContainer::ParticleType& p,
#ifdef USE_SPRAY_SOA
  const std::array<SprayParticleContainer::SoA::RealVector, NAR_SPR>& attribs,
  const int pid,
#endif
  pele::physics::transport::TransParm const* trans_parm)
{
  auto eos = pele::physics::PhysicsType::eos();
  SprayUnits SPU;
  const Real rule = 1. / 3.;
  const Real C_eps = 1.E-15;
  const Real B_eps = 1.E-7;
  const Real Y_eps = 1.E-4;
  const Real mass_eps = 8.E-18 * SPU.mass_conv;
  const int nSubMax = 100;
  bool get_xi = false;
  bool get_Ddiag = true;
  bool get_lambda = true;
  bool get_mu = true;
  if (!SPI.mass_tran) {
    get_Ddiag = false;
    get_lambda = false;
  }
  GpuArray<Real, NUM_SPECIES> Y_skin;
  GpuArray<Real, NUM_SPECIES> h_part;
  GpuArray<Real, NUM_SPECIES> h_fluid;
  GpuArray<Real, NUM_SPECIES> cp_n;
  GpuArray<Real, NUM_SPECIES> Ddiag;
  GpuArray<Real, SPRAY_FUEL_NUM> B_M_num;
  GpuArray<Real, SPRAY_FUEL_NUM> Sh_num;
  GpuArray<Real, SPRAY_FUEL_NUM> L_fuel;
  GpuArray<Real, SPRAY_FUEL_NUM> mi_dot;
  GpuArray<Real, SPRAY_FUEL_NUM> Y_part;
  GpuArray<Real, SPRAY_FUEL_NUM> Y_vapor;
  GpuArray<Real, SPRAY_FUEL_NUM>
    cBoilT; // Boiling temperature at current pressure
#ifdef USE_SPRAY_SOA
  RealVect vel_part(AMREX_D_DECL(
    attribs[SPI.pstateVel].data()[pid], attribs[SPI.pstateVel + 1].data()[pid],
    attribs[SPI.pstateVel + 2].data()[pid]));
  Real T_part = attribs[SPI.pstateT].data()[pid];
  Real dia_part = attribs[SPI.pstateDia].data()[pid];
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf)
    Y_part[spf] = attribs[SPI.pstateY + spf].data()[pid];
#else
  RealVect vel_part(AMREX_D_DECL(
    p.rdata(SPI.pstateVel), p.rdata(SPI.pstateVel + 1),
    p.rdata(SPI.pstateVel + 2)));
  Real T_part = p.rdata(SPI.pstateT);
  Real dia_part = p.rdata(SPI.pstateDia);
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf)
    Y_part[spf] = p.rdata(SPI.pstateY + spf);
#endif
  Real dt = flow_dt;
  int isub = 1;
  int nsub = 1;
  bool dosub = true; // Flag to turn subcycling off
  Real rho_part = 0.;
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf)
    rho_part += Y_part[spf] / fdat.rho[spf];
  rho_part = 1. / rho_part;
  Real pmass = M_PI / 6. * rho_part * std::pow(dia_part, 3);
  Real startmass = pmass;
  eos.T2Hi(gpv.T_fluid, h_fluid.data());
  calcBoilT(fdat, gpv, cBoilT.data());
  while (isub <= nsub) {
    // Model the fuel vapor using the one-third rule
    Real delT = amrex::max(gpv.T_fluid - T_part, 0.);
    Real T_skin = T_part + rule * delT;
    // Calculate the C_p at the skin temperature for each species
    eos.T2Cpi(T_skin, cp_n.data());
    eos.T2Hi(T_part, h_part.data());
    Real hg_s = 0.; // Enthalpy at particle surface temperature
    Real hg_g = 0.; // Enthalpy at gas temperature
    if (SPI.mass_tran) {
      for (int n = 0; n < NUM_SPECIES; ++n) {
        h_part[n] *= SPU.eng_conv;
        h_fluid[n] *= SPU.eng_conv;
        cp_n[n] *= SPU.eng_conv;
        Y_skin[n] = 0.;
        hg_s += gpv.Y_fluid[n] * h_part[n];
        hg_g += gpv.Y_fluid[n] * h_fluid[n];
      }
    } else {
      for (int n = 0; n < NUM_SPECIES; ++n) {
        Y_skin[n] = gpv.Y_fluid[n];
        h_part[n] *= SPU.eng_conv;
        cp_n[n] *= SPU.eng_conv;
      }
    }
    // Solve for state of the vapor and mass transfer coefficient B_M
    Real sumYSkin = 0.; // Mass fraction of the fuel in the vapor phase
    Real sumYFuel = 0.; // Mass fraction of the fuel in the gas phase
    Real cp_skin = 0.;  // Averaged C_p at particle surface
    Real cp_part = 0.;  // Cp of the liquid state
    Real mw_vap = 0.;   // Average molar mass of vapor mixture
    if (SPI.mass_tran) {
      calcVaporY(
        fdat, gpv, T_part, C_eps, Y_part.data(), h_part.data(), cBoilT.data(),
        Y_vapor.data(), L_fuel.data());
      for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
        const int fspec = fdat.indx[spf];
        Real Yfv = Y_vapor[spf];
        mi_dot[spf] = 0.;
        B_M_num[spf] =
          amrex::max(C_eps, (Yfv - gpv.Y_fluid[fspec]) / (1. - Yfv));
        Y_skin[fspec] = Yfv + rule * (gpv.Y_fluid[fspec] - Yfv);
        sumYSkin += Y_skin[fspec];
        cp_part += Y_part[spf] * fdat.cp[spf];
        sumYFuel += gpv.Y_fluid[fspec];
      }
      Real restYfluid = 1. - sumYFuel;
      Real restYSkin = 1. - sumYSkin;
      Real renorm = restYSkin / restYfluid;
      for (int sp = 0; sp < NUM_SPECIES; ++sp) {
        if (Y_skin[sp] == 0.) {
          Y_skin[sp] = gpv.Y_fluid[sp] * renorm;
        }
        cp_skin += Y_skin[sp] * cp_n[sp];
        mw_vap += Y_skin[sp] * gpv.invmw[sp];
      }
      mw_vap = 1. / mw_vap;
    } else {
      mw_vap = gpv.mw_mix;
    }
    Real lambda_skin = 0.;
    Real mu_skin = 0.;
    Real xi_skin = 0.;
    Real rho_skin = gpv.rho_fluid; // TODO: Check if this should be modeled
    Real rho_cgs = rho_skin / SPU.rho_conv;
    auto trans = pele::physics::PhysicsType::transport();
    trans.transport(
      get_xi, get_mu, get_lambda, get_Ddiag, T_skin, rho_cgs, Y_skin.data(),
      Ddiag.data(), mu_skin, xi_skin, lambda_skin, trans_parm);
    mu_skin *= SPU.mu_conv;
    lambda_skin *= SPU.lambda_conv;
    // Ensure gas is not all fuel to allow evaporation
    bool evap_fuel = (sumYFuel >= 1.) ? false : true;
    RealVect diff_vel = gpv.vel_fluid - vel_part;
    Real diff_vel_mag = diff_vel.vectorLength();
    // Local Reynolds number
    Real Reyn = rho_skin * diff_vel_mag * dia_part / mu_skin;
    Real Nu_0 = 1.;
    // Solve mass transfer source terms
    Real m_dot = 0.;
    if (SPI.mass_tran && evap_fuel) {
      Real Pr_skin = mu_skin * cp_skin / lambda_skin;
      Real powR = amrex::max(std::pow(Reyn, 0.077), 1.);
      Nu_0 = 1. + powR * std::cbrt(1. + Reyn * Pr_skin);
      for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
        if (Y_part[spf] > 0.) {
          const int fspec = fdat.indx[spf];
          // Convert mass diffusion coefficient from mixture average
          // to binary for fuel only, not concerned with other species
          Ddiag[fspec] *= mw_vap * gpv.invmw[fspec] * SPU.rhod_conv;
          const Real rhoD = Ddiag[fspec];
          const Real Sc_skin = mu_skin / rhoD;
          const Real B_M = B_M_num[spf];
          Real logB = std::log(1. + B_M);
          // Calculate Sherwood number and evaporation rate
          Real invFM = B_M / (logB * std::pow(1. + B_M, 0.7));
          Real Sh_0 = 1. + powR * std::cbrt(1. + Reyn * Sc_skin);
          Sh_num[spf] = 2. + (Sh_0 - 2.) * invFM;
          Real Tboil = cBoilT[spf];
          Real Tcrit = fdat.critT[spf];
          // If temperature exceeds critical temperature, mass is instantly
          // evaporated
          if (T_part > Tcrit) {
            mi_dot[spf] = -Y_part[spf] * pmass / (0.5 * flow_dt);
            dosub = false;
          } else if (T_part > Tboil) {
            Real delTb = T_part - Tboil;
            Real alpha = calcAlpha(delTb);
            Real Gf = M_PI * dia_part * dia_part * alpha * delTb / L_fuel[spf];
            Real dh = (hg_g - hg_s) / L_fuel[spf];
            Real coeff =
              M_PI * lambda_skin / cp_skin * dia_part * Sh_num[spf] * logB;
            Real G = calcFlashVaporRate(dh, coeff, Gf);
            mi_dot[spf] = -amrex::max(G + Gf, 0.);
          } else {
            mi_dot[spf] =
              -amrex::max(M_PI * rhoD * dia_part * Sh_num[spf] * logB, 0.);
          }
          m_dot += mi_dot[spf];
        }
        if (isub == 1 && dosub) {
          Real inv_tau_d = -m_dot / (3. * pmass);
          nsub =
            amrex::max(nsub, amrex::min(int(flow_dt * inv_tau_d) + 1, nSubMax));
        }
      }
    }
    // Solve momentum source terms
    const Real inv_pmass = 1. / pmass;
    RealVect part_mom_src(RealVect::TheZeroVector());
    if (SPI.mom_tran) {
      Real drag_coef = 0.;
      if (Reyn > 0.)
        drag_coef = (Reyn > 1.)
                      ? 24. / Reyn * (1. + std::cbrt(Reyn * Reyn) / 6.)
                      : 24. / Reyn;
      Real drag_force = 0.125 * rho_skin * drag_coef * M_PI * dia_part *
                        dia_part * diff_vel_mag;
      part_mom_src = drag_force * diff_vel;
      gpv.fluid_mom_src += part_mom_src;
      // s_d,mu dot u_d
      Real S_dmu_dot_u = part_mom_src.dotProduct(vel_part);
#ifndef SPRAY_PELE_LM
      gpv.fluid_eng_src += S_dmu_dot_u;
#endif
      if (isub == 1 && dosub) {
        Real inv_tau_var = drag_force * inv_pmass;
        nsub =
          amrex::min(amrex::max(nsub, int(flow_dt * inv_tau_var) + 1), nSubMax);
      }
    }
    // Solve for energy source terms
    Real part_temp_src = 0.;
    if (evap_fuel && SPI.mass_tran) {
      const Real inv_pm_cp = inv_pmass / cp_part;
      Real coeff_heat = 0.;
      for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
        if (Y_part[spf] > 0.) {
          const int fspec = fdat.indx[spf];
          Real ratio = cp_n[fspec] * Sh_num[spf] * Ddiag[fspec] / lambda_skin;
          Real heatC = calcHeatCoeff(ratio, B_M_num[spf], B_eps, C_eps, Nu_0);
          // Convection term
          coeff_heat += heatC;
          part_temp_src += mi_dot[spf] * L_fuel[spf];
        }
      }
      Real conv_src = M_PI * lambda_skin * dia_part * delT * coeff_heat;
      gpv.fluid_eng_src += conv_src;
      part_temp_src = (part_temp_src + conv_src) * inv_pm_cp;
      if (isub == 1 && delT > C_eps && dosub) {
        Real inv_tau_T = conv_src * inv_pm_cp / delT;
        nsub =
          amrex::min(amrex::max(nsub, int(flow_dt * inv_tau_T) + 1), nSubMax);
      }
    }
    if (isub == 1)
      dt = flow_dt / Real(nsub);
    const Real part_dt = 0.5 * dt;
    // Update particle components
    AMREX_D_TERM(vel_part[0] += part_dt * part_mom_src[0] * inv_pmass;
                 , vel_part[1] += part_dt * part_mom_src[1] * inv_pmass;
                 , vel_part[2] += part_dt * part_mom_src[2] * inv_pmass;);
    T_part += part_dt * part_temp_src;
    Real new_mass = pmass + m_dot * part_dt;
    if (new_mass > mass_eps) {
      rho_part = 0.;
      for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
        Y_part[spf] = amrex::min(
          1., amrex::max(
                0., (Y_part[spf] * pmass + mi_dot[spf] * part_dt) / new_mass));
        rho_part += Y_part[spf] / fdat.rho[spf];
      }
      rho_part = 1. / rho_part;
      pmass = new_mass;
      dia_part = std::cbrt(6. * pmass / (M_PI * rho_part));
    } else {
      m_dot -= new_mass / part_dt;
      pmass = 0.;
      for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
        Y_part[spf] = 0.;
      }
      p.id() = -1;
      isub = nsub + 1;
    }
    ++isub;
  }
  if (nsub > 1) {
    gpv.fluid_eng_src /= Real(nsub);
    gpv.fluid_mom_src /= Real(nsub);
  }
  // Must add any mass related sources at the end in case
  // some species disappear completely
  Real mdot_total = (pmass - startmass) / (0.5 * flow_dt);
  gpv.fluid_mass_src = mdot_total;
  Real part_ke = 0.5 * vel_part.radSquared();
  gpv.fluid_eng_src += part_ke * mdot_total;
  gpv.fluid_mom_src += vel_part * mdot_total;
  Real RT = pele::physics::Constants::RU * SPU.ru_conv * gpv.T_fluid;
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    const int fspec = fdat.indx[spf];
    Real oldY = p.rdata(SPI.pstateY + spf);
    Real newY = Y_part[spf];
    Real midot = (newY * pmass - oldY * startmass) / (0.5 * flow_dt);
    gpv.fluid_Y_dot[spf] = midot;
    gpv.fluid_eng_src += midot * (h_part[fspec] - RT * gpv.invmw[fspec]);
  }
#ifdef USE_SPRAY_SOA
  AMREX_D_TERM(attribs[SPI.pstateVel].data()[pid] = vel_part[0];
               , attribs[SPI.pstateVel + 1].data()[pid] = vel_part[1];
               , attribs[SPI.pstateVel + 2].data()[pid] = vel_part[2];);
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf)
    attribs[SPI.pstateY + spf].data()[pid] = Y_part[spf];
  attribs[SPI.pstateT].data()[pid] = T_part;
  attribs[SPI.pstateDia].data()[pid] = dia_part;
#else
  AMREX_D_TERM(p.rdata(SPI.pstateVel) = vel_part[0];
               , p.rdata(SPI.pstateVel + 1) = vel_part[1];
               , p.rdata(SPI.pstateVel + 2) = vel_part[2];);
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf)
    p.rdata(SPI.pstateY + spf) = Y_part[spf];
  p.rdata(SPI.pstateT) = T_part;
  p.rdata(SPI.pstateDia) = dia_part;
#endif
}

#endif
