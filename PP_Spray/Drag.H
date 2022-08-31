#ifndef DRAG_H
#define DRAG_H

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
  amrex::Real logB = std::log1p(B_T_old);
  amrex::Real invFT = B_T_old / (logB * std::pow(1. + B_T_old, 0.7));
  amrex::Real Nu_num = 2. + NU2 * invFT;
  phi = ratio / Nu_num;
  amrex::Real B_T = std::pow(BM1, phi) - 1.;
  amrex::Real error = amrex::Math::abs(B_T - B_T_old);
  int k = 0;
  while (k < maxIter && error > B_eps) {
    B_T_old = B_T;
    logB = std::log1p(B_T);
    invFT = B_T / (logB * std::pow(1. + B_T, 0.7));
    Nu_num = 2. + NU2 * invFT;
    phi = ratio / Nu_num;
    B_T = std::pow(BM1, phi) - 1.;
    error = amrex::Math::abs(B_T - B_T_old);
    k++;
  }
  logB = std::log1p(B_T);
  invFT = B_T / (logB * std::pow(1. + B_T, 0.7));
  Nu_num = 2. + NU2 * invFT;
  return Nu_num * logB / B_T;
}

// Compute the state in the vapor and skin phase
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void
calcVaporState(
  const SprayData& fdat,
  const GasPhaseVals& gpv,
  const amrex::Real& rule,
  const amrex::Real& T_in,
  const amrex::Real& C_eps,
  const amrex::Real* Y_l,
  const amrex::Real* h_part,
  const amrex::Real* cp_n,
  const amrex::Real* cBoilT,
  amrex::Real* Y_skin,
  amrex::Real* Psat,
  amrex::Real* L_fuel,
  amrex::Real& B_M,
  amrex::Real& sumXVap,
  amrex::Real& cp_skin,
  amrex::Real& mw_skin)
{
  SprayUnits SPU;
  amrex::Real RU = pele::physics::Constants::RU * SPU.ru_conv;
  amrex::Real PATM = pele::physics::Constants::PATM * SPU.pres_conv;

  amrex::Real sum1 = 0.;
  amrex::Real sum2 = 0.;
  amrex::Real sum3 = 0.;
  /*
    This computes the mass fractions in the vapor. The vapor mass fraction is
    y_vi = x_vi * mw_i / mw_vg where mw_vg = mw_gas * (1 - sum_k x_vk) + sum_k
    x_vk * mw_k and x_vk = x_lk * psat_k / p_gas Note: this function seems more
    convuluted because we start with the liquid mass fractions, not molar
    fractions
  */
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    if (Y_l[spf] > 1.E-12) {
      const int fspec = fdat.indx[spf];
      const amrex::Real boilT_ref = fdat.boilT[spf];
      const amrex::Real mw_fuel = gpv.mw_fluid[fspec];
      amrex::Real T_part = amrex::min(T_in, 0.999 * cBoilT[spf]);
      // Compute latent heat
      amrex::Real part_latent =
        h_part[fspec] + fdat.latent[spf] - fdat.cp[spf] * (T_part - fdat.ref_T);
      L_fuel[spf] = part_latent;
      amrex::Real pres_sat = 0.;
      // Using the Clasius-Clapeyron relation
      if (fdat.psat_coef[4 * spf + 3] == 0.) {
        pres_sat =
          PATM *
          std::exp(part_latent * mw_fuel / RU * (1. / boilT_ref - 1. / T_part));
        // Using the Antoine equation
      } else {
        pres_sat = fdat.psat(T_part, spf);
      }
      Psat[spf] = pres_sat;
      sum1 += Y_l[spf] * pres_sat / mw_fuel;
      sum2 += Y_l[spf] * pres_sat;
      sum3 += Y_l[spf] / mw_fuel;
    } else {
      Psat[spf] = 0.;
      L_fuel[spf] = 0.;
    }
  }
  amrex::Real totalmwx = gpv.mw_mix * (sum3 * gpv.p_fluid - sum1) + sum2;
  amrex::Real sumYSkin = 0.; // Mass fraction of fuel in the modeled skin phase
  amrex::Real sumYfFluid = 0.; // Mass fraction of fuel in the gas phase
  amrex::Real sumYVap = 0.;    // Mass fraction of fuel in the vapor phase
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    const int fspec = fdat.indx[spf];
    const amrex::Real mw_fuel = gpv.mw_fluid[fspec];
    amrex::Real Yfv = Y_l[spf] * Psat[spf] / totalmwx;
    if (Psat[spf] > gpv.X_fluid[fspec] * gpv.p_fluid) {
      sumYVap += Yfv;
      sumYfFluid += gpv.Y_fluid[fspec];
      amrex::Real Ysk = Yfv + rule * (gpv.Y_fluid[fspec] - Yfv);
      Y_skin[fspec] = Ysk;
      sumYSkin += Ysk;
      sumXVap += Psat[spf] * Y_l[spf] / mw_fuel;
    }
  }
  // Normalize skin mass fractions to ensure they sum to 1
  amrex::Real renorm = (1. - sumYSkin) / (1. - sumYfFluid);
  for (int n = 0; n < NUM_SPECIES; ++n) {
    if (Y_skin[n] == 0.) {
      Y_skin[n] = gpv.Y_fluid[n] * renorm;
    }
    cp_skin += Y_skin[n] * cp_n[n];
    mw_skin += Y_skin[n] * gpv.invmw[n];
  }
  mw_skin = 1. / mw_skin;
  B_M = (sumYVap - sumYfFluid) / amrex::max(C_eps, (1. - sumYVap));
  B_M = amrex::min(20., B_M);
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
  const amrex::Real C_eps = 1.E-15;
  const amrex::Real B_eps = 1.E-7;
  const amrex::Real min_mass = SPU.min_mass;
  const int nSubMax = 100;
  bool get_xi = false;
  bool get_Ddiag = true;
  bool get_lambda = true;
  bool get_mu = true;
  if (!fdat.mass_trans) {
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
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> Psat = {{0.0}};
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM>
    cBoilT; // Boiling temperature at current pressure
  amrex::RealVect vel_part(AMREX_D_DECL(
    p.rdata(SPI.pstateVel), p.rdata(SPI.pstateVel + 1),
    p.rdata(SPI.pstateVel + 2)));
  // If particle is fixed in place, make velocity zero
  if (fdat.fixed_parts) {
    vel_part = amrex::RealVect::TheZeroVector();
  }
  amrex::Real T_part = p.rdata(SPI.pstateT);
  fdat.calcBoilT(gpv, cBoilT.data());
  amrex::Real dia_part = p.rdata(SPI.pstateDia);
  amrex::Real rho_part = 0.;
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    Y_part[spf] = p.rdata(SPI.pstateY + spf);
    rho_part += Y_part[spf] / fdat.rhoL(amrex::min(T_part, cBoilT[spf]), spf);
  }
  rho_part = 1. / rho_part;
  amrex::Real dt = flow_dt;
  int isub = 1;
  int nsub = 1;
  amrex::Real pmass = M_PI / 6. * rho_part * std::pow(dia_part, 3);
  amrex::Real startmass = pmass;
  amrex::RealVect part_mom_src;
  while (isub <= nsub) {
    amrex::Real cp_part = 0.; // Cp of the liquid state
    amrex::Real Tcrit = 0.;   // Liquid mixture critical temperature
    amrex::Real Tboil = 0.;
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      Tcrit += Y_part[spf] * fdat.critT[spf];
      Tboil += Y_part[spf] * cBoilT[spf];
      cp_part += Y_part[spf] * fdat.cp[spf];
    }
    T_part = amrex::min(0.999 * Tcrit, T_part);
    // Model the fuel vapor using the one-third rule
    amrex::Real delT = amrex::max(gpv.T_fluid - T_part, 0.);
    amrex::Real T_skin = gpv.T_fluid;
    if (fdat.mass_trans) {
      T_skin = T_part + rule * delT;
    }
    // Calculate the C_p at the skin temperature for each species
    eos.T2Cpi(T_skin, cp_n.data());
    eos.T2Hi(T_part, h_part.data());
    for (int n = 0; n < NUM_SPECIES; ++n) {
      Y_skin[n] = 0.;
      h_part[n] *= SPU.eng_conv;
      cp_n[n] *= SPU.eng_conv;
    }
    amrex::Real cp_skin = 0.; // Average C_p in modeled skin phase
    amrex::Real mw_skin = 0.; // Average molar mass of skin phase
    amrex::Real B_M = 0.;     // Mass Spalding number
    amrex::Real sumXVap = 0.; // Sum of Y_L Psat_f / mw_f
    if (fdat.mass_trans) {
      calcVaporState(
        fdat, gpv, rule, T_part, C_eps, Y_part.data(), h_part.data(),
        cp_n.data(), cBoilT.data(), Y_skin.data(), Psat.data(), L_fuel.data(),
        B_M, sumXVap, cp_skin, mw_skin);
    } else {
      for (int n = 0; n < NUM_SPECIES; ++n) {
        Y_skin[n] = gpv.Y_fluid[n];
        cp_skin += gpv.Y_fluid[n] * cp_n[n];
      }
      mw_skin = gpv.mw_mix;
    }
    amrex::Real lambda_skin = 0.;
    amrex::Real mu_skin = 0.;
    amrex::Real xi_skin = 0.;
    amrex::Real rho_skin = gpv.rho_fluid;
    if (fdat.mass_trans) {
      rho_skin = mw_skin * gpv.p_fluid /
                 (pele::physics::Constants::RU * T_skin * SPU.ru_conv);
    }
    amrex::Real rho_cgs = rho_skin / SPU.rho_conv;
    auto trans = pele::physics::PhysicsType::transport();
    trans.transport(
      get_xi, get_mu, get_lambda, get_Ddiag, T_skin, rho_cgs, Y_skin.data(),
      Ddiag.data(), mu_skin, xi_skin, lambda_skin, trans_parm);
    mu_skin *= SPU.mu_conv;
    lambda_skin *= SPU.lambda_conv;
    amrex::RealVect diff_vel = gpv.vel_fluid - vel_part;
    amrex::Real diff_vel_mag = diff_vel.vectorLength();
    // Local Reynolds number
    amrex::Real Reyn = rho_skin * diff_vel_mag * dia_part / mu_skin;
    Reyn = amrex::min(Reyn, 1000.);

    // Solve momentum source terms
    const amrex::Real inv_pmass = 1. / pmass;
    amrex::Real drag_coef = 0.;
    if (Reyn > 0.) {
      drag_coef = (Reyn > 1.) ? 24. / Reyn * (1. + std::cbrt(Reyn * Reyn) / 6.)
                              : 24. / Reyn;
    }
    amrex::Real drag_force =
      0.125 * rho_skin * drag_coef * M_PI * dia_part * dia_part * diff_vel_mag;
    amrex::Real drag_force_p =
      0.75 * rho_skin * drag_coef * dia_part * diff_vel_mag / rho_part;
    part_mom_src = drag_force_p * diff_vel;
    if (fdat.mom_trans) {
      gpv.fluid_mom_src += drag_force * diff_vel;
#ifndef PELELM_USE_SPRAY
      // s_d,mu dot u_d
      amrex::Real S_dmu_dot_u = diff_vel.dotProduct(vel_part);
      gpv.fluid_eng_src += drag_force * S_dmu_dot_u;
#endif
      if (isub == 1) {
        amrex::Real inv_tau_var = drag_force_p;
        nsub =
          amrex::min(amrex::max(nsub, int(flow_dt * inv_tau_var) + 1), nSubMax);
      }
    }

    // Solve mass and energy transfer source terms
    amrex::Real part_temp_src = 0.;
    amrex::Real m_dot = 0.;
    // Only transfer heat and mass if gas phase is not saturated
    if (fdat.mass_trans && sumXVap > 0.) {
      amrex::Real Pr_skin = mu_skin * cp_skin / lambda_skin;
      amrex::Real powR =
        amrex::max(std::pow(amrex::min(400., Reyn), 0.077), 1.);
      amrex::Real rhoDtotal = 0.;
      amrex::Real maxPsat = 0.;
      amrex::Real Nu_num = 0.;
      int maxY = -1;
      for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
        const int fspec = fdat.indx[spf];
        // Convert mass diffusion coefficient from mixture average
        // to binary for fuel only, not concerned with other species
        Ddiag[fspec] *= mw_skin * gpv.invmw[fspec] * SPU.rhod_conv;
        // Get normalized vapor molar fraction
        amrex::Real X_norm =
          Y_part[spf] * gpv.invmw[fspec] * Psat[spf] / sumXVap;
        Ddiag[fspec] *= X_norm;
        rhoDtotal += Ddiag[fspec];
        if (Psat[spf] > maxPsat) {
          maxY = spf;
          maxPsat = Psat[spf];
        }
      }
      amrex::Real sumL = 0.;
      amrex::Real Nu_0 = 1. + powR * std::cbrt(1. + Reyn * Pr_skin);
      if (T_part > Tboil) {
        // Normal boiling of only the species with the highest saturation
        // pressure
        amrex::Real B_T = cp_skin * (gpv.T_fluid - T_part) / L_fuel[maxY];
        amrex::Real logBT = std::log1p(B_T);
        amrex::Real invFT = B_T / (logBT * std::pow(1. + B_T, 0.7));
        amrex::Real Nu_star = 2. + (Nu_0 - 2.) * invFT;
        amrex::Real coeff = M_PI * lambda_skin / cp_skin * dia_part * Nu_star;
        mi_dot[maxY] = -amrex::max(coeff * logBT, 0.);
        Nu_num = 0.; // Assume Q_L = 0
        m_dot += mi_dot[maxY];
        sumL += mi_dot[maxY] * L_fuel[maxY];
      } else {
        // Normal evaporation
        amrex::Real Sc_skin = mu_skin / rhoDtotal;
        amrex::Real logB = std::log1p(B_M);
        amrex::Real invFM = B_M / (logB * std::pow(1. + B_M, 0.7));
        amrex::Real Sh_0 = 1. + powR * std::cbrt(1. + Reyn * Sc_skin);
        amrex::Real Sh_num = 2. + (Sh_0 - 2.) * invFM;
        amrex::Real mdotcoeff = M_PI * dia_part * Sh_num * logB;
        amrex::Real ratio = cp_skin * Sh_num * rhoDtotal / lambda_skin;
        Nu_num = calcHeatCoeff(ratio, B_M, B_eps, C_eps, Nu_0);
        for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
          const int fspec = fdat.indx[spf];
          if (
            Y_part[spf] > 1.E-12 &&
            Psat[spf] > gpv.X_fluid[fspec] * gpv.p_fluid) {
            mi_dot[spf] = -amrex::max(Ddiag[fspec] * mdotcoeff, 0.);
            m_dot += mi_dot[spf];
            sumL += mi_dot[spf] * L_fuel[spf];
          }
        }
      }
      amrex::Real conv_src = M_PI * lambda_skin * dia_part * delT * Nu_num;
      gpv.fluid_eng_src += conv_src;
      part_temp_src = (sumL + conv_src) * inv_pmass / cp_part;
      if (isub == 1 && delT > C_eps) {
        // Limit dt so change in mass does not exceed 10%
        amrex::Real inv_tau_d = -m_dot / (0.2 * pmass);
        amrex::Real inv_tau_T = conv_src * inv_pmass / (cp_part * delT);
        nsub = amrex::min(
          amrex::max(
            nsub,
            amrex::max(int(flow_dt * inv_tau_d), int(flow_dt * inv_tau_T)) + 1),
          nSubMax);
      }
    }
    if (isub == 1) {
      dt = flow_dt / amrex::Real(nsub);
    }
    const amrex::Real part_dt = fdat.dtmod * dt;
    if (!fdat.fixed_parts) {
      // Update particle velocity
      AMREX_D_TERM(vel_part[0] += part_dt * part_mom_src[0];
                   , vel_part[1] += part_dt * part_mom_src[1];
                   , vel_part[2] += part_dt * part_mom_src[2];);
    }
    if (fdat.mass_trans) {
      T_part += part_dt * part_temp_src;
      T_part = amrex::min(0.999 * Tcrit, T_part);
      amrex::Real new_mass = pmass + m_dot * part_dt;
      if (new_mass > min_mass) {
        if (SPRAY_FUEL_NUM > 1) {
          rho_part = 0.;
          amrex::Real sumY = 0.;
          for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
            Y_part[spf] = amrex::min(
              1.,
              amrex::max(
                0., (Y_part[spf] * pmass + mi_dot[spf] * part_dt) / new_mass));
            if (Y_part[spf] <= 1.E-12) {
              Y_part[spf] = 0.;
            }
            sumY += Y_part[spf];
          }
          for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
            Y_part[spf] /= sumY;
            rho_part +=
              Y_part[spf] / fdat.rhoL(amrex::min(T_part, cBoilT[spf]), spf);
          }
          rho_part = 1. / rho_part;
        } else {
          rho_part = fdat.rhoL(amrex::min(T_part, cBoilT[0]), 0);
        }
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
  amrex::Real mdot_total = (pmass - startmass) / (fdat.dtmod * flow_dt);
  gpv.fluid_mass_src = mdot_total;
#ifndef PELELM_USE_SPRAY
  amrex::Real part_ke = 0.5 * vel_part.radSquared();
  gpv.fluid_eng_src += part_ke * mdot_total;
#endif
  if (fdat.mom_trans) {
    gpv.fluid_mom_src += vel_part * mdot_total;
  }
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    const int fspec = fdat.indx[spf];
    amrex::Real oldY = p.rdata(SPI.pstateY + spf);
    amrex::Real newY = Y_part[spf];
    amrex::Real midot =
      (newY * pmass - oldY * startmass) / (fdat.dtmod * flow_dt);
    gpv.fluid_Y_dot[spf] = midot;
    gpv.fluid_eng_src += midot * h_part[fspec];
    p.rdata(SPI.pstateY + spf) = Y_part[spf];
  }
  AMREX_D_TERM(p.rdata(SPI.pstateVel) = vel_part[0];
               , p.rdata(SPI.pstateVel + 1) = vel_part[1];
               , p.rdata(SPI.pstateVel + 2) = vel_part[2];);
  p.rdata(SPI.pstateT) = T_part;
  p.rdata(SPI.pstateDia) = dia_part;
}

#endif
