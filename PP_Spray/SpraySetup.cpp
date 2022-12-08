
#include "SprayParticles.H"

using namespace amrex;

std::string SprayParticleContainer::spray_fuel_names[SPRAY_FUEL_NUM];
std::string SprayParticleContainer::spray_dep_names[SPRAY_FUEL_NUM];
Vector<std::string> SprayParticleContainer::spray_derive_vars;

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
  std::vector<std::string> dep_fuel_names;
  std::vector<Real> crit_T;
  std::vector<Real> boil_T;
  std::vector<Real> spraycp;
  std::vector<Real> latent;
  std::vector<Real> sprayrho;
  std::vector<Real> mu(nfuel, 0.);
  std::vector<Real> lambda(nfuel, 0.);
  bool has_dep_spec = false;
  {
    pp.getarr("fuel_species", fuel_names);
    pp.getarr("fuel_crit_temp", crit_T);
    pp.getarr("fuel_boil_temp", boil_T);
    pp.getarr("fuel_cp", spraycp);
    pp.getarr("fuel_latent", latent);
    pp.getarr("fuel_rho", sprayrho);
    pp.queryarr("fuel_mu", mu);
    pp.queryarr("fuel_lambda", lambda);
    if (pp.contains("dep_fuel_species")) {
      has_dep_spec = true;
      pp.getarr("dep_fuel_species", dep_fuel_names);
    }
    for (int i = 0; i < nfuel; ++i) {
      spray_fuel_names[i] = fuel_names[i];
      if (has_dep_spec) {
        spray_dep_names[i] = dep_fuel_names[i];
      } else {
        spray_dep_names[i] = spray_fuel_names[i];
      }
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
  bool splash_model = false;
  bool breakup_model = false;
  //
  // Set the number of particles per parcel
  //
  pp.query("parcel_size", parcel_size);
  pp.query("use_splash_model", splash_model);
  pp.query("use_breakup_model", breakup_model);
  if (splash_model || breakup_model) {
    if (!pp.contains("fuel_sigma") || !pp.contains("fuel_mu")) {
      Print()
        << "fuel_sigma and fuel_mu must be set for splash or breakup model. "
        << std::endl;
      Abort();
    }
    if (splash_model) {
      // TODO: Have this retrieved from proper boundary data
      pp.get("wall_temp", sprayData.wall_T);
    }
    // Set the fuel surface tension and contact angle
    pp.get("fuel_sigma", sprayData.sigma);
    sprayData.do_splash = splash_model;
    sprayData.do_breakup = breakup_model;
  }

  // Must use same reference temperature for all fuels
  // TODO: This means the reference temperature must be the same for all fuel
  // species
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

  sprayData.num_ppp = parcel_size;
  sprayData.ref_T = spray_ref_T;

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
    spray_derive_vars.push_back(derive_name);
  }
  if (derive_plot_species == 1 && SPRAY_FUEL_NUM > 1) {
    for (auto& fuel_name : spray_fuel_names) {
      spray_derive_vars.push_back("spray_mass_" + fuel_name);
    }
  }

  if (particle_verbose >= 1 && ParallelDescriptor::IOProcessor()) {
    amrex::Print() << "Spray fuel species " << spray_fuel_names[0];
    for (int i = 1; i < SPRAY_FUEL_NUM; ++i) {
      amrex::Print() << ", " << spray_fuel_names[i];
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
SprayParticleContainer::spraySetup(SprayData& sprayData)
{
#if NUM_SPECIES > 1
  Vector<std::string> spec_names;
  pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(
    spec_names);
  for (int i = 0; i < SPRAY_FUEL_NUM; ++i) {
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
      std::string gas_spec = spec_names[ns];
      if (gas_spec == spray_fuel_names[i]) {
        sprayData.indx[i] = ns;
      }
      if (gas_spec == spray_dep_names[i]) {
        sprayData.dep_indx[i] = ns;
      }
    }
    if (sprayData.indx[i] < 0) {
      amrex::Print() << "Fuel " << spray_fuel_names[i]
                     << " not found in species list" << std::endl;
      amrex::Abort();
    }
    if (sprayData.dep_indx[i] < 0) {
      amrex::Print() << "Fuel " << spray_dep_names[i]
                     << " not found in species list" << std::endl;
      amrex::Abort();
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
}
