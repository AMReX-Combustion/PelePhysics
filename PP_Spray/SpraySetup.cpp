
#include "SprayParticles.H"

using namespace amrex;

std::string SprayParticleContainer::m_sprayFuelNames[SPRAY_FUEL_NUM];
std::string SprayParticleContainer::m_sprayDepNames[SPRAY_FUEL_NUM];
Vector<std::string> SprayParticleContainer::m_sprayDeriveVars;

void
getInpCoef(
  Real* coef,
  const ParmParse& ppp,
  const std::string* fuel_names,
  const std::string& varname,
  bool is_required = false)
{
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    std::string var_read = fuel_names[spf] + "_" + varname;
    int numvals = ppp.countval(var_read.c_str());
    // If 4 values are specified, assume fit coefficients
    if (numvals == 4) {
      std::vector<Real> inp_coef(4, 0.);
      if (is_required) {
        ppp.getarr(var_read.c_str(), inp_coef);
      } else {
        ppp.queryarr(var_read.c_str(), inp_coef);
      }
      for (int i = 0; i < 4; ++i) {
        coef[4 * spf + i] = inp_coef[i];
      }
    } else if (numvals == 1) {
      // If 1 value is specified, assume constant value
      Real inp_coef = 0.;
      for (int i = 0; i < 4; ++i) {
        coef[4 * spf + i] = 0.;
      }
      if (is_required) {
        ppp.get(var_read.c_str(), inp_coef);
      } else {
        ppp.query(var_read.c_str(), inp_coef);
      }
      coef[4 * spf] = inp_coef;
    }
  }
}

void
getInpVal(
  Real* coef,
  const ParmParse& ppp,
  const std::string* fuel_names,
  const std::string& varname)
{
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    std::string var_read = fuel_names[spf] + "_" + varname;
    ppp.get(var_read.c_str(), coef[spf]);
  }
}

void
SprayParticleContainer::readSprayParams(
  int& particle_verbose,
  Real& particle_cfl,
  int& write_spray_ascii_files,
  int& plot_spray_src,
  int& init_function,
  std::string& init_file,
  SprayData& sprayData,
  const Real& max_cfl)
{
  ParmParse pp("particles");
  //
  // Control the verbosity of the Particle class
  pp.query("v", particle_verbose);

  pp.query("mass_transfer", sprayData.mass_trans);
  pp.query("mom_transfer", sprayData.mom_trans);
  pp.query("fixed_parts", sprayData.fixed_parts);
  pp.query("cfl", particle_cfl);
  if (particle_cfl > max_cfl) {
    Abort("particles.cfl must be <= " + std::to_string(max_cfl));
  }
  // Number of fuel species in spray droplets
  // Must match the number specified at compile time
  const int nfuel = pp.countval("fuel_species");
  if (nfuel != SPRAY_FUEL_NUM) {
    Abort("Number of fuel species in input file must match SPRAY_FUEL_NUM");
  }

  std::vector<std::string> fuel_names;
  std::vector<std::string> dep_fuel_names;
  bool has_dep_spec = false;
  {
    pp.getarr("fuel_species", fuel_names);
    if (pp.contains("dep_fuel_species")) {
      has_dep_spec = true;
      pp.getarr("dep_fuel_species", dep_fuel_names);
    }
    getInpVal(sprayData.critT.data(), pp, fuel_names.data(), "crit_temp");
    getInpVal(sprayData.boilT.data(), pp, fuel_names.data(), "boil_temp");
    getInpVal(sprayData.cp.data(), pp, fuel_names.data(), "cp");
    getInpVal(sprayData.ref_latent.data(), pp, fuel_names.data(), "latent");

    getInpCoef(sprayData.lambda_coef.data(), pp, fuel_names.data(), "lambda");
    getInpCoef(sprayData.psat_coef.data(), pp, fuel_names.data(), "psat");
    getInpCoef(sprayData.rho_coef.data(), pp, fuel_names.data(), "rho", true);
    getInpCoef(sprayData.mu_coef.data(), pp, fuel_names.data(), "mu");
    for (int i = 0; i < nfuel; ++i) {
      m_sprayFuelNames[i] = fuel_names[i];
      if (has_dep_spec) {
        m_sprayDepNames[i] = dep_fuel_names[i];
      } else {
        m_sprayDepNames[i] = m_sprayFuelNames[i];
      }
      sprayData.latent[i] = sprayData.ref_latent[i];
    }
  }

  Real spray_ref_T = 300.;
  bool splash_model = false;
  int breakup_model = 0;
  //
  // Set the number of particles per parcel
  //
  pp.query("max_parcel_size", max_num_ppp);
  pp.query("use_splash_model", splash_model);
  std::string breakup_model_str = "None";
  pp.query("use_breakup_model", breakup_model_str);
  if (breakup_model_str == "TAB") {
    breakup_model = 1;
    pp.query("max_parcel_size", max_num_ppp);
  } else if (breakup_model_str == "KHRT") {
    breakup_model = 2;
    pp.query("KHRT_B0", m_khrtB0);
    pp.query("KHRT_B1", m_khrtB1);
    pp.query("KHRT_C3", m_khrtC3);
  } else if (breakup_model_str == "None") {
    breakup_model = 0;
  } else {
    Abort("'use_breakup_model' input not recognized. Must be 'TAB', 'KHRT', or "
          "'None'");
  }
  if (splash_model || (breakup_model > 0)) {
    pp.query("breakup_parcel_factor", breakup_ppp_fact);
    if (breakup_ppp_fact > 1. || breakup_ppp_fact < 0.) {
      Abort("'breakup_parcel_factor' must be between 0 and 1");
    }
    bool wrong_data = false;
    for (int i = 0; i < nfuel; ++i) {
      std::string var_read = fuel_names[i] + "_mu";
      if (!pp.contains(var_read.c_str())) {
        wrong_data = true;
      }
    }
    if (wrong_data || !pp.contains("fuel_sigma")) {
      Abort(
        "fuel_sigma and mu coeffs must be set for splash or breakup model.");
    }
    if (splash_model) {
      // TODO: Have this retrieved from proper boundary data
      pp.get("wall_temp", sprayData.wall_T);
      Real theta_c_deg = -1.;
      pp.get("contact_angle", theta_c_deg);
      if (theta_c_deg < 0. || theta_c_deg > 180.) {
        Abort("'contact_angle' must be between 0 and 180");
      }
      sprayData.theta_c = theta_c_deg * M_PI / 180.;
    }
    // Set the fuel surface tension and contact angle
    pp.get("fuel_sigma", sprayData.sigma);
    sprayData.do_splash = splash_model;
    sprayData.do_breakup = breakup_model;
  }

  // Must use same reference temperature for all fuels
  pp.get("fuel_ref_temp", spray_ref_T);
  //
  // Set if spray ascii files should be written
  //
  pp.query("write_ascii_files", write_spray_ascii_files);
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

  sprayData.ref_T = spray_ref_T;

  // List of known derived spray quantities
  std::vector<std::string> derive_names = {
    "spray_mass",      // Total liquid mass in a cell
    "spray_density",   // Liquid mass divided by cell volume
    "spray_num",       // Number of spray droplets in a cell
    "spray_vol",       // Total liquid volume in a cell
    "spray_surf_area", // Total liquid surface area in a cell
    "spray_vol_frac",  // Volume fraction of liquid in cell
    "d10",             // Average diameter
    "d32",             // SMD
    "wall_film_hght",  // Wall film height
    "spray_temp",      // Mass-weighted average temperature
    "num_parcels",     // Number of parcels in a cell
    AMREX_D_DECL("spray_x_vel", "spray_y_vel", "spray_z_vel")};
  int derive_plot_vars = 1;
  pp.query("derive_plot_vars", derive_plot_vars);
  int derive_plot_species = 1;
  pp.query("derive_plot_species", derive_plot_species);
  // If derive_spray_vars if present, add above spray quantities in the same
  // order
  for (const auto& derive_name : derive_names) {
    m_sprayDeriveVars.push_back(derive_name);
  }
  if (derive_plot_species == 1 && SPRAY_FUEL_NUM > 1) {
    for (auto& fuel_name : m_sprayFuelNames) {
      m_sprayDeriveVars.push_back("spray_mass_" + fuel_name);
    }
  }

  if (particle_verbose >= 1 && ParallelDescriptor::IOProcessor()) {
    Print() << "Spray fuel species " << m_sprayFuelNames[0];
#if SPRAY_FUEL_NUM > 1
    for (int i = 1; i < SPRAY_FUEL_NUM; ++i) {
      Print() << ", " << m_sprayFuelNames[i];
    }
#endif
    Print() << std::endl;
    Print() << "Number of particles per parcel " << parcel_size << std::endl;
  }
  Gpu::streamSynchronize();
  ParallelDescriptor::Barrier();
}

void
SprayParticleContainer::spraySetup(SprayData& sprayData, const Real* body_force)
{
#if NUM_SPECIES > 1
  Vector<std::string> spec_names;
  pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(
    spec_names);
  for (int i = 0; i < SPRAY_FUEL_NUM; ++i) {
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
      std::string gas_spec = spec_names[ns];
      if (gas_spec == m_sprayFuelNames[i]) {
        sprayData.indx[i] = ns;
      }
      if (gas_spec == m_sprayDepNames[i]) {
        sprayData.dep_indx[i] = ns;
      }
    }
    if (sprayData.indx[i] < 0) {
      Abort("Fuel " + m_sprayFuelNames[i] + " not found in species list");
    }
    if (sprayData.dep_indx[i] < 0) {
      Abort("Fuel " + m_sprayDepNames[i] + " not found in species list");
    }
  }
#else
  sprayData.indx[0] = 0;
  sprayData.dep_indx[0] = 0;
#endif
  SprayUnits SPU;
  Vector<Real> fuelEnth(NUM_SPECIES);
  auto eos = pele::physics::PhysicsType::eos();
  eos.T2Hi(sprayData.ref_T, fuelEnth.data());
  for (int ns = 0; ns < SPRAY_FUEL_NUM; ++ns) {
    const int fspec = sprayData.indx[ns];
    sprayData.latent[ns] -= fuelEnth[fspec] * SPU.eng_conv;
  }
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    sprayData.body_force[dir] = body_force[dir];
  }
  Gpu::streamSynchronize();
  ParallelDescriptor::Barrier();
}
