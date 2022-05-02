#ifndef DRAG_H_
#define DRAG_H_

#include "SprayParticles.H"
#include <AMReX_AmrParticles.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_ParIter.H>
#include <AMReX_Particles.H>
#include <AMReX_REAL.H>
#include <Transport.H>

// Compute the heat transfer coefficient using the
// corrected Nusselt number and B_T value
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
calcHeatCoeff(
  const amrex::Real& ratio,
  const amrex::Real& B_M,
  const amrex::Real& B_eps,
  const amrex::Real& C_eps,
  const amrex::Real& Nu_0)
{
  if (B_M <= C_eps) {
    return 0.;
  }
  const int maxIter = 100;
  const amrex::Real NU2 = Nu_0 - 2.;
  const amrex::Real BM1 = 1. + B_M;
  amrex::Real phi = ratio / Nu_0;
  amrex::Real B_T_old = std::pow(BM1, phi) - 1.;
  amrex::Real BT1 = 1. + B_T_old;
  amrex::Real logB = std::log(BT1);
  amrex::Real invFT = B_T_old / (logB * std::pow(BT1, 0.7));
  amrex::Real Nu_num = 2. + NU2 * invFT;
  phi = ratio / Nu_num;
  amrex::Real B_T = std::pow(BM1, phi) - 1.;
  amrex::Real error = std::abs(B_T - B_T_old);
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

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
rhoL(const SprayData& fdat, const amrex::Real& T, const int spf)
{
  if (fdat.rho_coef[4 * spf] == 0.) {
    return fdat.rho[spf];
  } else {
    amrex::Real a = fdat.rho_coef[4 * spf];
    amrex::Real b = fdat.rho_coef[4 * spf + 1];
    amrex::Real c = fdat.rho_coef[4 * spf + 2];
    amrex::Real d = fdat.rho_coef[4 * spf + 3];
    return a + T * (b + T * (c + T * d));
  }
}

// Compute the flash boiling coefficient
// From Adachi 1997
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
calcAlpha(const amrex::Real& delTb)
{
  SprayUnits SPU;
  amrex::Real alpha = 0.;
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
amrex::Real
calcFlashVaporRate(
  const amrex::Real& dh, const amrex::Real& coeff, const amrex::Real& Gf)
{
  if (dh <= 1.E-5) {
    return 0.;
  }
  amrex::Real tol = 1.E-4;
  amrex::Real Grat = 1.E-5;
  amrex::Real G = 0.;
  amrex::Real error = 100.;
  while (error > tol) {
    amrex::Real Gold = Grat;
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
calcBoilT(const SprayData& fdat, const GasPhaseVals& gpv, amrex::Real* cBoilT)
{
  SprayUnits SPU;
  amrex::Real RU = pele::physics::Constants::RU * SPU.ru_conv;
  amrex::Real PATM = pele::physics::Constants::PATM * SPU.pres_conv;
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    const int fspec = fdat.indx[spf];
    const amrex::Real critT = fdat.critT[spf];
    const amrex::Real boilT_ref = fdat.boilT[spf];
    const amrex::Real mw_fuel = gpv.mw_fluid[fspec];
    // Since we only know the latent heat at the reference temperature,
    // modify Watsons power law to find latent heat at boiling conditions
    amrex::Real Hboil_ref =
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
  const amrex::Real& T_in,
  const amrex::Real& C_eps,
  const amrex::Real* Y_l,
  const amrex::Real* h_part,
  const amrex::Real* cBoilT,
  amrex::Real* Y_vapor,
  amrex::Real* Psat,
  amrex::Real* L_fuel)
{
  SprayUnits SPU;
  amrex::Real RU = pele::physics::Constants::RU * SPU.ru_conv;
  amrex::Real PATM = pele::physics::Constants::PATM * SPU.pres_conv;
  amrex::Real sum_xv = 0.;    // Sum x_v
  amrex::Real sum_mw_xv = 0.; // Sum mw_v*x_v
  amrex::Real nt = 0.;        // Total moles of liquid
  // This computes the mass fractions in the vapor
  // The vapor mass fraction is y_vi = x_vi * mw_i / mw_vg
  // where mw_vg = mw_gas * (1 - sum_k x_vk) + sum_k x_vk * mw_k
  // and x_vk = x_lk * psat_k / p_gas
  // Note: this function seems more convuluted because we start with
  // the liquid mass fractions, not molar fractions
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    const int fspec = fdat.indx[spf];
    const amrex::Real boilT_ref = fdat.boilT[spf];
    const amrex::Real mw_fuel = gpv.mw_fluid[fspec];
    amrex::Real T_part = amrex::min(T_in, 0.999 * cBoilT[spf]);
    // Compute latent heat
    amrex::Real part_latent =
      h_part[fspec] + fdat.latent[spf] - fdat.cp[spf] * (T_part - fdat.ref_T);
    L_fuel[spf] = part_latent;
    // Compute the saturation pressure at the current particle temperature
    amrex::Real a = fdat.psat_coef[4 * spf];
    amrex::Real b = fdat.psat_coef[4 * spf + 1];
    amrex::Real c = fdat.psat_coef[4 * spf + 2];
    amrex::Real d =
      fdat.psat_coef[4 * spf + 3]; // For converting to Pa or dyne/cm^2
    amrex::Real pres_sat = 0.;
    // Using the Clasius-Clapeyron relation
    if (d == 0.) {
      pres_sat =
        PATM *
        std::exp(part_latent * mw_fuel / RU * (1. / boilT_ref - 1. / T_part));
      // Using the Antoine equation
    } else {
      pres_sat = d * std::pow(10., a - b / (T_part + c));
    }
    Psat[spf] = pres_sat;
    amrex::Real x_l = Y_l[spf] / mw_fuel;
    amrex::Real x_vc = x_l * pres_sat;
    nt += x_l;
    sum_xv += x_vc;
    Y_vapor[spf] = mw_fuel * x_vc;
    sum_mw_xv += Y_vapor[spf];
  }
  amrex::Real totalmwx = gpv.mw_mix * (nt * gpv.p_fluid - sum_xv) + sum_mw_xv;
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    amrex::Real Yfv = Y_vapor[spf] / totalmwx;
    Y_vapor[spf] = amrex::max(0., amrex::min(1. - C_eps, Yfv));
  }
}

// Compute source terms and update particles
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
calculateSpraySource(
  const amrex::Real flow_dt,
  GasPhaseVals& gpv,
  SprayComps SPI,
  SprayData fdat,
  SprayParticleContainer::ParticleType& p,
  pele::physics::transport::TransParm<
    pele::physics::EosType,
    pele::physics::TransportType> const* trans_parm)
{
  auto eos = pele::physics::PhysicsType::eos();
  SprayUnits SPU;
  const amrex::Real rule = 1. / 3.;
  // Advance half dt per function call
  const amrex::Real dtmod = 0.5;
  const amrex::Real C_eps = 1.E-15;
  const amrex::Real B_eps = 1.E-7;
  const amrex::Real mass_eps = 8.E-18 * SPU.mass_conv;
  const int nSubMax = 100;
  bool get_xi = false;
  bool get_Ddiag = true;
  bool get_lambda = true;
  bool get_mu = true;
  if (!SPI.mass_tran) {
    get_Ddiag = false;
    get_lambda = false;
  }
  amrex::GpuArray<amrex::Real, NUM_SPECIES> Y_skin;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> h_part;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> cp_n;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> Ddiag;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> L_fuel = {{0.0}};
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> mi_dot = {{0.0}};
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> Y_part; // Liquid mass fractions
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> Y_vapor;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> Psat = {{0.0}};
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM>
    cBoilT; // Boiling temperature at current pressure
  amrex::RealVect vel_part(AMREX_D_DECL(
    p.rdata(SPI.pstateVel), p.rdata(SPI.pstateVel + 1),
    p.rdata(SPI.pstateVel + 2)));
  amrex::Real T_part = p.rdata(SPI.pstateT);
  calcBoilT(fdat, gpv, cBoilT.data());
  amrex::Real dia_part = p.rdata(SPI.pstateDia);
  amrex::Real rho_part = 0.;
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    Y_part[spf] = p.rdata(SPI.pstateY + spf);
    rho_part += Y_part[spf] / rhoL(fdat, amrex::min(T_part, cBoilT[spf]), spf);
  }
  rho_part = 1. / rho_part;
  amrex::Real dt = flow_dt;
  int isub = 1;
  int nsub = 1;
  amrex::Real pmass = M_PI / 6. * rho_part * std::pow(dia_part, 3);
  amrex::Real startmass = pmass;
  while (isub <= nsub) {
    amrex::Real cp_part = 0.; // Cp of the liquid state
    amrex::Real Tcrit = 0.;   // Liquid mixture critical temperature
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      Tcrit += Y_part[spf] * fdat.critT[spf];
      cp_part += Y_part[spf] * fdat.cp[spf];
    }
    T_part = amrex::min(0.999 * Tcrit, T_part);
    // Model the fuel vapor using the one-third rule
    amrex::Real delT = amrex::max(gpv.T_fluid - T_part, 0.);
    amrex::Real T_skin = T_part + rule * delT;
    // Calculate the C_p at the skin temperature for each species
    eos.T2Cpi(T_skin, cp_n.data());
    eos.T2Hi(T_part, h_part.data());
    if (SPI.mass_tran) {
      for (int n = 0; n < NUM_SPECIES; ++n) {
        h_part[n] *= SPU.eng_conv;
        cp_n[n] *= SPU.eng_conv;
        Y_skin[n] = 0.;
      }
    } else {
      for (int n = 0; n < NUM_SPECIES; ++n) {
        Y_skin[n] = gpv.Y_fluid[n];
        h_part[n] *= SPU.eng_conv;
        cp_n[n] *= SPU.eng_conv;
      }
    }
    // Solve for state of the vapor and mass transfer coefficient B_M
    amrex::Real sumYSkin =
      0.; // Mass fraction of fuel in the modeled skin phase
    amrex::Real sumYfFluid = 0.; // Mass fraction of fuel in the gas phase
    amrex::Real sumYVap = 0.;    // Mass fraction of fuel in the vapor phase
    amrex::Real cp_skin = 0.;    // Average C_p in modeled skin phase
    amrex::Real mw_skin = 0.;    // Average molar mass of vapor mixture
    amrex::Real B_M = 0.;        // Mass Spalding number
    amrex::Real Sh_num = 0.;     // Sherwood number
    amrex::Real psattotal = 0.;  // Sum of Y_L Psat_f / mw_f
    if (SPI.mass_tran) {
      calcVaporY(
        fdat, gpv, T_part, C_eps, Y_part.data(), h_part.data(), cBoilT.data(),
        Y_vapor.data(), Psat.data(), L_fuel.data());
      for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
        const int fspec = fdat.indx[spf];
        amrex::Real Yfv = Y_vapor[spf];
        sumYVap += Yfv;
        sumYfFluid += gpv.Y_fluid[fspec];
        Y_skin[fspec] = Yfv + rule * (gpv.Y_fluid[fspec] - Yfv);
        sumYSkin += Y_skin[fspec];
        psattotal += Psat[spf] * Y_part[spf] * gpv.invmw[fspec];
      }
      // Normalize skin mass fractions to ensure they sum to 1
      amrex::Real renorm = (1. - sumYSkin) / (1. - sumYfFluid);
      for (int sp = 0; sp < NUM_SPECIES; ++sp) {
        if (Y_skin[sp] == 0.) {
          Y_skin[sp] = gpv.Y_fluid[sp] * renorm;
        }
        cp_skin += Y_skin[sp] * cp_n[sp];
        mw_skin += Y_skin[sp] * gpv.invmw[sp];
      }
      mw_skin = 1. / mw_skin;
      B_M = (sumYVap - sumYfFluid) / (1. - sumYVap);
      B_M = amrex::min(20., B_M);
    } else {
      mw_skin = gpv.mw_mix;
    }
    amrex::Real lambda_skin = 0.;
    amrex::Real mu_skin = 0.;
    amrex::Real xi_skin = 0.;
    amrex::Real rho_skin =
      mw_skin * gpv.p_fluid /
      (pele::physics::Constants::RU * T_skin * SPU.ru_conv);
    amrex::Real rho_cgs = rho_skin / SPU.rho_conv;
    auto trans = pele::physics::PhysicsType::transport();
    trans.transport(
      get_xi, get_mu, get_lambda, get_Ddiag, T_skin, rho_cgs, Y_skin.data(),
      Ddiag.data(), mu_skin, xi_skin, lambda_skin, trans_parm);
    mu_skin *= SPU.mu_conv;
    lambda_skin *= SPU.lambda_conv;
    // Ensure gas is not all fuel to allow evaporation
    bool evap_fuel = (sumYfFluid >= 1.) ? false : true;
    amrex::RealVect diff_vel = gpv.vel_fluid - vel_part;
    amrex::Real diff_vel_mag = diff_vel.vectorLength();
    // Local Reynolds number
    amrex::Real Reyn = rho_skin * diff_vel_mag * dia_part / mu_skin;
    Reyn = amrex::min(Reyn, 400.);
    amrex::Real Nu_0 = 1.;
    // Solve mass transfer source terms
    amrex::Real m_dot = 0.;
    amrex::Real rhoDtotal = 0.;
    if (SPI.mass_tran && evap_fuel) {
      amrex::Real Pr_skin = mu_skin * cp_skin / lambda_skin;
      amrex::Real powR = amrex::max(std::pow(Reyn, 0.077), 1.);
      Nu_0 = 1. + powR * std::cbrt(1. + Reyn * Pr_skin);
      for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
        const int fspec = fdat.indx[spf];
        // Convert mass diffusion coefficient from mixture average
        // to binary for fuel only, not concerned with other species
        Ddiag[fspec] *= mw_skin * gpv.invmw[fspec] * SPU.rhod_conv;
        // Get normalized vapor molar fraction
        amrex::Real X_norm =
          Y_part[spf] * gpv.invmw[fspec] * Psat[spf] / psattotal;
        Ddiag[fspec] *= X_norm;
        rhoDtotal += Ddiag[fspec];
      }
      amrex::Real Sc_skin = mu_skin / rhoDtotal;
      amrex::Real logB = std::log(1. + B_M);
      amrex::Real invFM = B_M / (logB * std::pow(1. + B_M, 0.7));
      amrex::Real Sh_0 = 1. + powR * std::cbrt(1. + Reyn * Sc_skin);
      Sh_num = 2. + (Sh_0 - 2.) * invFM;
      amrex::Real mdotcoeff = M_PI * dia_part * Sh_num * logB;
      // TODO: Determine if this is usefule
      /* if (T_part > Tboil) { */
      /*   amrex::Real hg_s = 0.; */
      /*   eos.TY2H(T_part, Y_skin.data(), hg_s); */
      /*   hg_s *= SPU.eng_conv; */
      /*   amrex::Real delTb = amrex::max(0., T_part - Tboil); */
      /*   amrex::Real alpha = calcAlpha(delTb); */
      /*   amrex::Real Gf = */
      /*     M_PI * dia_part * dia_part * alpha * delTb / L_fuel; */
      /*   amrex::Real dh = (h_fluid - hg_s) / Lfuel; */
      /*   amrex::Real coeff = */
      /*     M_PI * lambda_skin / cp_skin * dia_part * Sh_num; */
      /*   amrex::Real G = calcFlashVaporRate(dh, coeff, Gf); */
      /*   mdotboil = amrex::max(G + Gf, 0.); */
      /* } */
      for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
        const int fspec = fdat.indx[spf];
        if (Y_part[spf] > 0. && Psat[spf] > gpv.X_fluid[fspec] * gpv.p_fluid) {
          mi_dot[spf] = -amrex::max(Ddiag[fspec] * mdotcoeff, 0.);
          m_dot += mi_dot[spf];
        }
        if (isub == 1) {
          amrex::Real inv_tau_d = -m_dot / (3. * pmass);
          nsub =
            amrex::max(nsub, amrex::min(int(flow_dt * inv_tau_d) + 1, nSubMax));
        }
      }
    }
    // Solve momentum source terms
    const amrex::Real inv_pmass = 1. / pmass;
    amrex::RealVect part_mom_src(amrex::RealVect::TheZeroVector());
    if (SPI.mom_tran) {
      amrex::Real drag_coef = 0.;
      if (Reyn > 0.) {
        drag_coef = (Reyn > 1.)
                      ? 24. / Reyn * (1. + std::cbrt(Reyn * Reyn) / 6.)
                      : 24. / Reyn;
      }
      amrex::Real drag_force = 0.125 * rho_skin * drag_coef * M_PI * dia_part *
                               dia_part * diff_vel_mag;
      part_mom_src = drag_force * diff_vel;
      gpv.fluid_mom_src += part_mom_src;
#ifndef PELELM_USE_SPRAY
      // s_d,mu dot u_d
      amrex::Real S_dmu_dot_u = part_mom_src.dotProduct(vel_part);
      gpv.fluid_eng_src += S_dmu_dot_u;
#endif
      if (isub == 1) {
        amrex::Real inv_tau_var = drag_force * inv_pmass;
        nsub =
          amrex::min(amrex::max(nsub, int(flow_dt * inv_tau_var) + 1), nSubMax);
      }
    }
    // Solve for energy source terms
    amrex::Real part_temp_src = 0.;
    if (evap_fuel && SPI.mass_tran) {
      const amrex::Real inv_pm_cp = inv_pmass / cp_part;
      amrex::Real ratio = cp_skin * Sh_num * rhoDtotal / lambda_skin;
      amrex::Real coeff_heat = calcHeatCoeff(ratio, B_M, B_eps, C_eps, Nu_0);
      for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
        part_temp_src += mi_dot[spf] * L_fuel[spf];
      }
      amrex::Real conv_src = M_PI * lambda_skin * dia_part * delT * coeff_heat;
      gpv.fluid_eng_src += conv_src;
      part_temp_src = (part_temp_src + conv_src) * inv_pm_cp;
      if (isub == 1 && delT > C_eps) {
        amrex::Real inv_tau_T = conv_src * inv_pm_cp / delT;
        nsub =
          amrex::min(amrex::max(nsub, int(flow_dt * inv_tau_T) + 1), nSubMax);
      }
    }
    if (isub == 1) {
      dt = flow_dt / amrex::Real(nsub);
    }
    const amrex::Real part_dt = dtmod * dt;
    // Update particle components
    AMREX_D_TERM(vel_part[0] += part_dt * part_mom_src[0] * inv_pmass;
                 , vel_part[1] += part_dt * part_mom_src[1] * inv_pmass;
                 , vel_part[2] += part_dt * part_mom_src[2] * inv_pmass;);
    if (SPI.mass_tran) {
      T_part += part_dt * part_temp_src;
      T_part = amrex::min(Tcrit, T_part);
      amrex::Real new_mass = pmass + m_dot * part_dt;
      if (new_mass > mass_eps) {
        rho_part = 0.;
        for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
          Y_part[spf] = amrex::min(
            1.,
            amrex::max(
              0., (Y_part[spf] * pmass + mi_dot[spf] * part_dt) / new_mass));
          rho_part +=
            Y_part[spf] / rhoL(fdat, amrex::min(T_part, cBoilT[spf]), spf);
        }
        rho_part = 1. / rho_part;
        pmass = new_mass;
        dia_part = std::cbrt(6. * pmass / (M_PI * rho_part));
      } else {
        pmass = 0.;
        p.id() = -1;
        nsub = isub;
        isub = nsub + 1;
      }
    }
    ++isub;
  }
  if (nsub > 1) {
    gpv.fluid_eng_src /= amrex::Real(nsub);
    gpv.fluid_mom_src /= amrex::Real(nsub);
  }
  // Must add any mass related sources at the end in case
  // some species disappear completely
  amrex::Real mdot_total = (pmass - startmass) / (dtmod * flow_dt);
  gpv.fluid_mass_src = mdot_total;
#ifndef PELELM_USE_SPRAY
  amrex::Real part_ke = 0.5 * vel_part.radSquared();
  gpv.fluid_eng_src += part_ke * mdot_total;
#endif
  if (SPI.mom_tran) {
    gpv.fluid_mom_src += vel_part * mdot_total;
  }
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    const int fspec = fdat.indx[spf];
    amrex::Real oldY = p.rdata(SPI.pstateY + spf);
    amrex::Real newY = Y_part[spf];
    amrex::Real midot = (newY * pmass - oldY * startmass) / (dtmod * flow_dt);
    gpv.fluid_Y_dot[spf] = midot;
    gpv.fluid_eng_src += midot * h_part[fspec];
  }
  AMREX_D_TERM(p.rdata(SPI.pstateVel) = vel_part[0];
               , p.rdata(SPI.pstateVel + 1) = vel_part[1];
               , p.rdata(SPI.pstateVel + 2) = vel_part[2];);
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    p.rdata(SPI.pstateY + spf) = Y_part[spf];
  }
  p.rdata(SPI.pstateT) = T_part;
  p.rdata(SPI.pstateDia) = dia_part;
}

#endif
