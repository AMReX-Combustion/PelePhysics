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
  amrex::Real error = std::abs(B_T - B_T_old);
  int k = 0;
  while (k < maxIter && error > B_eps) {
    B_T_old = B_T;
    logB = std::log1p(B_T);
    invFT = B_T / (logB * std::pow(1. + B_T, 0.7));
    Nu_num = 2. + NU2 * invFT;
    phi = ratio / Nu_num;
    B_T = std::pow(BM1, phi) - 1.;
    error = std::abs(B_T - B_T_old);
    k++;
  }
  logB = std::log1p(B_T);
  invFT = B_T / (logB * std::pow(1. + B_T, 0.7));
  Nu_num = 2. + NU2 * invFT;
  return Nu_num * logB / B_T;
}

// Compute the state in the vapor and skin phase
//
template <typename EOSType>
struct calcVaporState
{
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void operator()(
    const SprayData& fdat,
    const GasPhaseVals& gpv,
    const amrex::Real& rule,
    const amrex::Real& T_in,
    const amrex::Real& C_eps,
    const amrex::Real& mw_part,
    const amrex::Real* Y_l,
    const amrex::Real T_skin,
    const amrex::Real* cBoilT,
    amrex::Real* Y_skin,
    amrex::Real* X_vapor,
    amrex::Real* L_fuel,
    amrex::Real& B_M,
    amrex::Real& sumXVap,
    amrex::Real& cp_skin,
    amrex::Real& mw_skin)
  {
    auto eos = pele::physics::PhysicsType::eos(fdat.eosparm);
    amrex::Real h_gas[NUM_SPECIES];
    eos.T2Hi(T_in, h_gas);

    amrex::Real mbar_vap = 0.; // Average molar mass of evaporated fuel in vapor
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      const int fspec = fdat.indx[spf];
      const int fdspec = fdat.dep_indx[spf];
      const amrex::Real mw_li = fdat.liqprops.mw_i(spf);
      amrex::Real T_part = amrex::min<amrex::Real>(T_in, 0.99 * cBoilT[spf]);

      // Compute latent heat
      amrex::Real part_latent = fdat.liqprops.latent_i(
        T_part, h_gas[fspec] * SprayUnits::eng_conv, spf);
      L_fuel[spf] = part_latent;
      amrex::Real pres_sat = fdat.liqprops.psat_i(T_part, spf, part_latent);

      amrex::Real X_fluid = gpv.Y_fluid[fdspec] * gpv.mw_mix / mw_li;
      if (pres_sat > X_fluid * gpv.p_fluid) {
        amrex::Real Xl = Y_l[spf] * mw_part / mw_li;
        amrex::Real Xv = Xl * pres_sat / gpv.p_fluid;
        mbar_vap += Xv * mw_li;
        sumXVap += Xv;
        X_vapor[spf] = Xv;
      } else {
        X_vapor[spf] = 0.;
      }
    }
    amrex::Real totalmwx = mbar_vap + (1. - sumXVap) * gpv.mw_mix;
    amrex::Real sumYSkin =
      0.; // Mass fraction of fuel in the modeled skin phase
    amrex::Real sumYfFluid = 0.; // Mass fraction of fuel in the gas phase
    amrex::Real sumYVap = 0.;    // Mass fraction of fuel in the vapor phase
    amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> Yfv{0.0};
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      if (X_vapor[spf] > 0.) {
        const amrex::Real mw_li = fdat.liqprops.mw_i(spf);
        Yfv[spf] = X_vapor[spf] * mw_li / totalmwx;
        sumYVap += Yfv[spf];
      }
    }
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      if (X_vapor[spf] > 0.) {
        const int fspec = fdat.indx[spf];
        const int fdspec = fdat.dep_indx[spf];
        amrex::Real gasY = gpv.Y_fluid[fdspec];
        // If multiple species share a dep spec, split gas phase across them
        // in proportion to their contribution at equilibrium
        if (fdat.liquid_spec_share_gas_dep) {
          amrex::Real Yvap_sum = 0.0;
          for (int spf2 = 0; spf2 < SPRAY_FUEL_NUM; ++spf2) {
            if (fdat.dep_indx[spf2] == fdspec) {
              Yvap_sum += Yfv[spf2];
            }
          }
          gasY *= (Yfv[spf] / Yvap_sum);
        }
        sumYfFluid += gasY;
        amrex::Real Ysk = Yfv[spf] + rule * (gasY - Yfv[spf]);
        Y_skin[fspec] = Ysk;
        sumYSkin += Ysk;
      }
    }
    // Normalize skin mass fractions to ensure they sum to 1
    amrex::Real renorm = (1. - sumYSkin) / (1. - sumYfFluid);
    for (int n = 0; n < NUM_SPECIES; ++n) {
      if (Y_skin[n] == 0.) {
        Y_skin[n] = gpv.Y_fluid[n] * renorm;
      }
    }
    // eos and unit conversions
    eos.Y2WBAR(Y_skin, mw_skin);
    eos.TY2Cp(T_skin, Y_skin, cp_skin);
    mw_skin *= SprayUnits::mass_conv;
    cp_skin *= SprayUnits::eng_conv;
    B_M =
      (sumYVap - sumYfFluid) / amrex::max<amrex::Real>(C_eps, (1. - sumYVap));
    B_M = amrex::min<amrex::Real>(20., B_M);
  }
};

template <>
struct calcVaporState<pele::physics::eos::Manifold>
{
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void operator()(
    const SprayData& fdat,
    const GasPhaseVals& gpv,
    const amrex::Real& rule,
    const amrex::Real& T_in,
    const amrex::Real& C_eps,
    const amrex::Real& mw_part,
    const amrex::Real* Y_l,
    const amrex::Real T_skin,
    const amrex::Real* cBoilT,
    amrex::Real* Y_skin,
    amrex::Real* X_vapor,
    amrex::Real* L_fuel,
    amrex::Real& B_M,
    amrex::Real& sumXVap,
    amrex::Real& cp_skin,
    amrex::Real& mw_skin)
  {
    amrex::Real mbar_vap = 0.; // mean molar mass of evaporated fuel in vapor
    auto eos = pele::physics::PhysicsType::eos(fdat.eosparm);

    // Calculate X_{v,n}- the vapour state

    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      const int fdspec = fdat.dep_indx[spf];
      const amrex::Real mw_li = fdat.liqprops.mw_i(spf);
      amrex::Real T_part = amrex::min<amrex::Real>(T_in, 0.99 * cBoilT[spf]);

      // Compute latent heat using correlation
      amrex::Real part_latent =
        fdat.liqprops.latent_correlation_i(fdat.liqprops.boilT_i(spf), spf);
      L_fuel[spf] = part_latent;

      // psat - any approach should work because we populated part_latent
      amrex::Real pres_sat = fdat.liqprops.psat_i(T_part, spf, part_latent);

      //
      amrex::Real X_fluid = gpv.Yspec_fluid[fdspec] * gpv.mw_mix / mw_li;

      if (pres_sat > X_fluid * gpv.p_fluid) {
        amrex::Real Xl = Y_l[spf] * mw_part / mw_li; // X_{d,n} in documentation
        amrex::Real Xv =
          Xl * pres_sat / gpv.p_fluid; // X_{v,n} in documentation
        mbar_vap += Xv * mw_li;
        sumXVap += Xv;
        X_vapor[spf] = Xv;
      } else {
        X_vapor[spf] = 0.;
      }
    }

    // Y_{g,i}^{pc}  is already stored in  gpv.Yspec_fluid
    // Calculate Y_{v,i}, Y_{r,i}^{pc}, and sum Y_{r,i}^{pc}
    amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> Yfv{0.0};
    amrex::Real sumYVap = 0.0;
    amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> Y_ri_pc{0.0};
    amrex::Real sumYripc = 0.0;
    amrex::Real sumYgipc = 0.0;
    amrex::Real totalmwx = mbar_vap + (1. - sumXVap) * gpv.mw_mix;
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      if (X_vapor[spf] > 0.) {
        const amrex::Real mw_li = fdat.liqprops.mw_i(spf);
        Yfv[spf] = X_vapor[spf] * mw_li / totalmwx;
        sumYVap += Yfv[spf]; //\sigma Y_{v,k}^{pc} in the document
        Y_ri_pc[spf] = Yfv[spf] + rule * (gpv.Yspec_fluid[spf] - Yfv[spf]);
        sumYripc += Y_ri_pc[spf];         //\sigma Y_{r,k}^{pc} in the document
        sumYgipc += gpv.Yspec_fluid[spf]; //\sigma Y_{g,k}^{pc} in the document
      }
    }

    // Now calculate theta
    amrex::Real theta = (1.0 - sumYripc) / (1.0 - sumYgipc);

    // Now calculate \zi_{r,j}
    amrex::GpuArray<amrex::Real, NUM_SPECIES> zi_rj{0.0};
    amrex::GpuArray<amrex::Real, NUM_SPECIES> WYgi{0.0}; // W_ij Y_{g,i}^{pc}
    amrex::GpuArray<amrex::Real, NUM_SPECIES> WYri{0.0}; // W_ij Y_{r,i}^{pc}

    // Alternative way of doing above without using W_ij assuming W_ij is some
    // row exchanged form of the identity matrix

    for (int nls = 0; nls < SPRAY_FUEL_NUM; nls++) {
      if (X_vapor[nls] > 0.) {
        WYgi[fdat.dep_indx[nls]] = gpv.Yspec_fluid[nls];
        WYri[fdat.dep_indx[nls]] = Y_ri_pc[nls];
      }
    }

    for (int n = 0; n < NUM_SPECIES; n++) {
      zi_rj[n] = WYri[n] + theta * (gpv.Y_fluid[n] - WYgi[n]);
      Y_skin[n] = zi_rj[n];
    }

    eos.Y2WBAR(Y_skin, mw_skin);
    mw_skin *= SprayUnits::mass_conv;

    // T_skin is just a dummy here
    eos.TY2Cp(T_skin, Y_skin, cp_skin);
    cp_skin *= SprayUnits::eng_conv;

    B_M = (sumYVap - sumYgipc) / amrex::max<amrex::Real>(C_eps, (1. - sumYVap));
    B_M = amrex::min<amrex::Real>(20., B_M);
  }
};

// Compute source terms and update particles
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
amrex::Real
calculateSpraySource(
  const amrex::Real flow_dt,
  GasPhaseVals& gpv,
  const SprayData& fdat,
  SprayParticleContainer::ParticleType& p,
  amrex::Real* cBoilT,
  pele::physics::transport::TransParm<
    pele::physics::EosType,
    pele::physics::TransportType> const* trans_parm)
{
  auto eos = pele::physics::PhysicsType::eos(fdat.eosparm);
  const amrex::Real rule = 1. / 3.;
  const amrex::Real C_eps = 1.E-15;
  const amrex::Real B_eps = 1.E-7;
  const amrex::Real min_mass = SprayUnits::min_mass;
  const int nSubMax = 100;
  bool get_xi = false;
  bool get_Ddiag = true;
  bool get_lambda = true;
  bool get_mu = true;
  bool get_chi = false;
  if (!fdat.mass_trans) {
    get_Ddiag = false;
    get_lambda = false;
  }
  amrex::GpuArray<amrex::Real, NUM_SPECIES> Y_skin;
  amrex::GpuArray<amrex::Real, NUM_SPECIES> Ddiag;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> L_fuel;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> mi_dot;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> Y_part; // Liquid mass fractions
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> X_vapor;

  Y_skin.fill(0.0);
  Ddiag.fill(0.0);
  mi_dot.fill(0.0);
  L_fuel.fill(0.0);
  Y_part.fill(0.0);
  X_vapor.fill(0.0);

  const amrex::Real num_ppp = p.rdata(SprayComps::pstateNumDens);
  amrex::RealVect vel_part(AMREX_D_DECL(
    p.rdata(SprayComps::pstateVel), p.rdata(SprayComps::pstateVel + 1),
    p.rdata(SprayComps::pstateVel + 2)));
  // TAB model distortion
  amrex::Real y_tab = 0.;
  if (fdat.do_breakup == 1) {
    y_tab = p.rdata(SprayComps::pstateBM1);
  }
  // If particle is fixed in place, make velocity zero
  if (fdat.fixed_parts) {
    vel_part = amrex::RealVect::TheZeroVector();
  }
  amrex::Real T_part = p.rdata(SprayComps::pstateT);
  amrex::Real dia_part = p.rdata(SprayComps::pstateDia);

  // Storing Y_{d,n}
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    Y_part[spf] = p.rdata(SprayComps::pstateY + spf);
  }
  // Calculating \rho_{d}
  amrex::Real rho_part = fdat.liqprops.rho_mix(Y_part, T_part, cBoilT);

  amrex::Real dt = flow_dt;
  int isub = 1;
  int nsub = 1;
  amrex::Real pmass = M_PI / 6. * rho_part * std::pow(dia_part, 3);
  amrex::Real startmass = pmass;
  amrex::Real Reyn;
  amrex::RealVect part_vel_src;
  while (isub <= nsub) {
    const amrex::Real cp_part = fdat.liqprops.cp_mix(Y_part, T_part, cBoilT);
    const amrex::Real Tboil = fdat.liqprops.boilT_mix(Y_part, cBoilT);
    // Average molar mass of liquid droplet
    const amrex::Real mw_part = fdat.liqprops.mw_mix(Y_part);
    // Model the fuel vapor using the one-third rule
    amrex::Real delT = amrex::max<amrex::Real>(gpv.T_fluid - T_part, 0.);
    amrex::Real T_skin = gpv.T_fluid;
    if (fdat.mass_trans) {
      T_skin = T_part + rule * delT;
    }
    Y_skin.fill(0.0);
    amrex::Real cp_skin = 0.; // Average C_p in modeled skin phase
    amrex::Real mw_skin = 0.; // Average molar mass of skin phase
    amrex::Real B_M = 0.;     // Mass Spalding number
    amrex::Real sumXVap = 0.; // Sum of Y_L Psat_f / mw_f

    if (fdat.mass_trans) {
      calcVaporState<pele::physics::PhysicsType::eos_type>()(
        fdat, gpv, rule, T_part, C_eps, mw_part, Y_part.data(), T_skin, cBoilT,
        Y_skin.data(), X_vapor.data(), L_fuel.data(), B_M, sumXVap, cp_skin,
        mw_skin);
    } else {
      for (int n = 0; n < NUM_SPECIES; ++n) {
        Y_skin[n] = gpv.Y_fluid[n];
      }
      eos.TY2Cp(T_skin, Y_skin.data(), cp_skin);
      cp_skin *= SprayUnits::eng_conv;
      mw_skin = gpv.mw_mix;
    }
    amrex::Real lambda_skin = 0.;
    amrex::Real mu_skin = 0.;
    amrex::Real xi_skin = 0.;
    amrex::Real rho_skin, rho_cgs;
    if (fdat.mass_trans) {
      eos.PYT2R(
        gpv.p_fluid / SprayUnits::pres_conv, Y_skin.data(), T_skin, rho_cgs);
      rho_skin = rho_cgs * SprayUnits::rho_conv;
    } else {
      rho_skin = gpv.rho_fluid;
      rho_cgs = rho_skin / SprayUnits::rho_conv;
    }
    auto trans = pele::physics::PhysicsType::transport();

    trans.transport(
      get_xi, get_mu, get_lambda, get_Ddiag, get_chi, T_skin, rho_cgs,
      Y_skin.data(), Ddiag.data(), nullptr, mu_skin, xi_skin, lambda_skin,
      trans_parm);

    mu_skin *= SprayUnits::mu_conv;
    lambda_skin *= SprayUnits::lambda_conv;
#ifdef USE_MANIFOLD_EOS
    // Hack: manifold model may inaccurately predict temperature in skin layer
    // to improve predictions, adjust computed transport coefficients to the
    // skin layer temperature using a simple temperature scaling based on
    // Sutherland's law. Rho also needs to be modified (uses ideal gas law)
    amrex::Real T_skin_predicted = 300.0;
    eos.HY2T(0.0, Y_skin.data(), T_skin_predicted);
    const amrex::Real temp_ratio = T_skin / T_skin_predicted;
    const amrex::Real sutherland_temp = 110.4;
    const amrex::Real sutherland_ratio = std::sqrt(temp_ratio) * temp_ratio *
                                         (T_skin_predicted + sutherland_temp) /
                                         (T_skin + sutherland_temp);
    mu_skin *= sutherland_ratio;
    lambda_skin *= sutherland_ratio;
    for (int spf = 0; spf < NUM_SPECIES; ++spf) {
      Ddiag[spf] *= sutherland_ratio;
    }
    rho_skin /= temp_ratio;
#endif
    amrex::RealVect diff_vel = gpv.vel_fluid - vel_part;
    amrex::Real diff_vel_mag = diff_vel.vectorLength();
    // Local Reynolds number
    Reyn = rho_skin * diff_vel_mag * dia_part / mu_skin;

    // Solve momentum source terms
    const amrex::Real inv_pmass = 1. / pmass;

    amrex::Real drag_coef = 0.;
    if (Reyn > 1000.) {
      drag_coef = 0.424;
    } else if (Reyn > 1.) {
      drag_coef = 24. / Reyn * (1. + std::cbrt(Reyn * Reyn) / 6.);
    } else if (Reyn > 0.) {
      drag_coef = 24. / Reyn;
    }
    drag_coef *= (1. + 2.6432 * y_tab);
    amrex::Real drag_force =
      0.125 * rho_skin * drag_coef * M_PI * dia_part * dia_part * diff_vel_mag;
    amrex::Real drag_force_p =
      0.75 * rho_skin * drag_coef * diff_vel_mag / (dia_part * rho_part);
    part_vel_src = drag_force_p * diff_vel + fdat.body_force;

    if (fdat.mom_trans) {
      gpv.fluid_mom_src += num_ppp * drag_force * diff_vel;
#ifndef PELELM_USE_SPRAY
      // s_d,mu dot u_d
      amrex::Real S_dmu_dot_u = diff_vel.dotProduct(vel_part);
      gpv.fluid_eng_src += num_ppp * drag_force * S_dmu_dot_u;
#endif
      if (isub == 1) {
        amrex::Real inv_tau_var = drag_force_p;
        nsub = amrex::min<int>(
          amrex::max<int>(nsub, static_cast<int>(flow_dt * inv_tau_var) + 1),
          nSubMax);
      }
    }

    // Solve mass and energy transfer source terms
    amrex::Real part_temp_src = 0.;
    amrex::Real m_dot = 0.;
    if (fdat.mass_trans) {
      amrex::Real Pr_skin = mu_skin * cp_skin / lambda_skin;

      amrex::Real powR = amrex::max<amrex::Real>(
        std::pow(amrex::min<amrex::Real>(400., Reyn), 0.077), 1.);
      amrex::Real Nu_num = 0.;
      amrex::Real sumL = 0.;
      amrex::Real Nu_0 = 1. + powR * std::cbrt(1. + Reyn * Pr_skin);
      // If gas phase is not saturated
      if (sumXVap > 0.) {
        amrex::Real rhoDtotal = 0.;
        amrex::Real cp_fuel = 0.;
        amrex::Real sumYSkin = 0.;
        for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
          const int fspec = fdat.indx[spf];
          // Convert mass diffusion coefficient from mixture average
          // to binary for fuel only, not concerned with other species

          Ddiag[fspec] *= mw_skin / gpv.mw[fspec] * SprayUnits::rhod_conv;
          // Get normalized vapor molar fraction
          amrex::Real X_norm = X_vapor[spf] / sumXVap;
          Ddiag[fspec] *= X_norm;
          rhoDtotal += Ddiag[fspec];
          sumYSkin += Y_skin[fspec];
        }

#ifndef USE_MANIFOLD_EOS
        amrex::GpuArray<amrex::Real, NUM_SPECIES> Y_fuel_skin;
        Y_fuel_skin.fill(0.0);
        for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
          Y_fuel_skin[fdat.indx[spf]] = Y_skin[fdat.indx[spf]] / sumYSkin;
        }
        eos.TY2Cp(T_skin, Y_fuel_skin.data(), cp_fuel);
#else
        eos.Y2GenericManifoldData(
          Y_skin.data(), &(fdat.eosparm->idx_Cp_fuel), &cp_fuel, 1);
        // TODO: account for variation between table temperature and skin
        // temperature
        //  cp_fuel = 1e7 + (cp_fuel - 1e7) * T_skin / T_skin_predicted;
#endif
        cp_fuel *= SprayUnits::eng_conv;

        amrex::Real Sc_skin = mu_skin / rhoDtotal;
        amrex::Real logB = std::log1p(B_M);
        amrex::Real invFM = B_M / (logB * std::pow(1. + B_M, 0.7));
        amrex::Real Sh_0 = 1. + powR * std::cbrt(1. + Reyn * Sc_skin);
        amrex::Real Sh_num = 2. + (Sh_0 - 2.) * invFM;
        amrex::Real mdotcoeff = M_PI * dia_part * Sh_num * logB;
        amrex::Real ratio = cp_fuel * Sh_num * rhoDtotal / lambda_skin;

        Nu_num = calcHeatCoeff(ratio, B_M, B_eps, C_eps, Nu_0);
        for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
          // Species index
          const int fspec = fdat.indx[spf];
          if (X_vapor[spf] > 0.) {
            mi_dot[spf] =
              -amrex::max<amrex::Real>(Ddiag[fspec] * mdotcoeff, 0.);
            m_dot += mi_dot[spf];
            sumL += mi_dot[spf] * L_fuel[spf];
          }
        }
      } else {
        // Still do convective heat transfer for saturated gas
        Nu_num = Nu_0;
      }

      amrex::Real conv_src = M_PI * lambda_skin * dia_part * delT * Nu_num;
      gpv.fluid_eng_src += num_ppp * conv_src;
      part_temp_src = (sumL + conv_src) * inv_pmass / cp_part;
      if (isub == 1 && delT > C_eps) {
        // Limit dt so change in mass does not exceed 10%
        amrex::Real inv_tau_d = -m_dot / (0.2 * pmass);
        amrex::Real inv_tau_T = conv_src * inv_pmass / (cp_part * delT);
        nsub = amrex::min<int>(
          amrex::max<int>(
            nsub, amrex::max<int>(
                    static_cast<int>(flow_dt * inv_tau_d),
                    static_cast<int>(flow_dt * inv_tau_T)) +
                    1),
          nSubMax);
      }
    }

    if (isub == 1) {
      dt = flow_dt / static_cast<amrex::Real>(nsub);
    }
    const amrex::Real part_dt = fdat.dtmod * dt;
    if (!fdat.fixed_parts) {
      // Update particle velocity
      AMREX_D_TERM(
        vel_part[0] += part_dt * part_vel_src[0];
        , vel_part[1] += part_dt * part_vel_src[1];
        , vel_part[2] += part_dt * part_vel_src[2];);
    }
    if (fdat.mass_trans) {
      amrex::Real newT = T_part + part_dt * part_temp_src;
      amrex::Real new_mass = pmass + m_dot * part_dt;
      // If droplet is still reasonable size and temperature
      if (new_mass > min_mass && newT < Tboil) {
        T_part = newT;
#if SPRAY_FUEL_NUM > 1
        amrex::Real sumY = 0.;
        for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
          Y_part[spf] = amrex::min<amrex::Real>(
            1.,
            amrex::max<amrex::Real>(
              0., (Y_part[spf] * pmass + mi_dot[spf] * part_dt) / new_mass));
          if (Y_part[spf] <= 1.E-12) {
            Y_part[spf] = 0.;
          }
          sumY += Y_part[spf];
        }
        for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
          Y_part[spf] /= sumY;
        }
        rho_part = fdat.liqprops.rho_mix(Y_part, T_part, cBoilT);
#else
        rho_part =
          fdat.liqprops.rho_i(amrex::min<amrex::Real>(T_part, cBoilT[0]), 0);
#endif
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
    gpv.fluid_eng_src /= static_cast<amrex::Real>(nsub);
    gpv.fluid_mom_src /= static_cast<amrex::Real>(nsub);
  }
  // Must add any mass related sources at the end in case
  // some species disappear completely
  amrex::Real mdot_total = (pmass - startmass) / (fdat.dtmod * flow_dt);
  gpv.fluid_mass_src = num_ppp * mdot_total;
#ifndef PELELM_USE_SPRAY
  amrex::Real part_ke = 0.5 * vel_part.radSquared();
  gpv.fluid_eng_src += num_ppp * part_ke * mdot_total;
#endif
  if (fdat.mom_trans) {
    gpv.fluid_mom_src += num_ppp * vel_part * mdot_total;
  }

#ifndef USE_MANIFOLD_EOS
  amrex::Real h_gas[NUM_SPECIES];
  eos.T2Hi(T_part, h_gas);
#endif
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    const int fdspec = fdat.dep_indx[spf];
    amrex::Real oldY = p.rdata(SprayComps::pstateY + spf);
    amrex::Real newY = Y_part[spf];
    amrex::Real midot =
      (newY * pmass - oldY * startmass) / (fdat.dtmod * flow_dt);
    gpv.fluid_Y_dot[spf] = num_ppp * midot;
#ifdef USE_MANIFOLD_EOS
    gpv.fluid_eng_src = 0.0;
#else
    gpv.fluid_eng_src += num_ppp * midot * h_gas[fdspec] * SprayUnits::eng_conv;
#endif
    p.rdata(SprayComps::pstateY + spf) = Y_part[spf];
  }
  AMREX_D_TERM(p.rdata(SprayComps::pstateVel) = vel_part[0];
               , p.rdata(SprayComps::pstateVel + 1) = vel_part[1];
               , p.rdata(SprayComps::pstateVel + 2) = vel_part[2];);
  p.rdata(SprayComps::pstateT) = T_part;
  p.rdata(SprayComps::pstateDia) = dia_part;
  return Reyn;
}

#endif
