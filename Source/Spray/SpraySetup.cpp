
#include "SprayParticles.H"

using namespace amrex;

std::string SprayParticleContainer::m_sprayFuelNames[SPRAY_FUEL_NUM];
std::string SprayParticleContainer::m_sprayDepNames[SPRAY_FUEL_NUM];
Vector<std::string> SprayParticleContainer::m_sprayDeriveVars;
SprayData* SprayParticleContainer::m_sprayData = nullptr;
SprayData* SprayParticleContainer::d_sprayData = nullptr;
SprayComps SprayParticleContainer::m_sprayIndx;
Real SprayParticleContainer::spray_cfl = 0.5;
bool SprayParticleContainer::write_ascii_files = false;
bool SprayParticleContainer::plot_spray_src = false;
Real SprayParticleContainer::m_maxNumPPP = 100.;
Real SprayParticleContainer::m_breakupPPPFact = 0.5;
Real SprayParticleContainer::m_khrtB0 = 0.61;
Real SprayParticleContainer::m_khrtB1 = 7.;
Real SprayParticleContainer::m_khrtC3 = 1.;
std::string SprayParticleContainer::spray_init_file;

void
SprayParticleContainer::readSprayParams(int& particle_verbose)
{
  amrex::Print() << "\n Reading spray model parameters ..." << std::endl;
#if AMREX_SPACEDIM == 1
  amrex::Abort("Spray model not valid in 1D");
#elif AMREX_SPACEDIM == 2
  amrex::Print()
    << " Warning: Spray model in 2D assumes narrow domain in z-direction (Lz = "
       "dz)!"
    << std::endl;
#endif
  m_sprayData = new SprayData{};
  d_sprayData =
    static_cast<SprayData*>(amrex::The_Arena()->alloc(sizeof(SprayData)));
  ParmParse pp("particles");
  // Control the verbosity of the Particle class
  pp.query("v", particle_verbose);

  pp.query("mass_transfer", m_sprayData->mass_trans);
  pp.query("mom_transfer", m_sprayData->mom_trans);
  pp.query("fixed_parts", m_sprayData->fixed_parts);
#ifdef PELELM_USE_SPRAY
  Real max_cfl = 2.;
#else
  Real max_cfl = 0.5;
#endif
  pp.query("cfl", spray_cfl);
  if (spray_cfl > max_cfl) {
    Abort("particles.cfl must be <= " + std::to_string(max_cfl));
  }
  // Number of fuel species in spray droplets
  // Must match the number specified at compile time
  const int nfuel = pp.countval("fuel_species");
  if (nfuel != SPRAY_FUEL_NUM) {
    amrex::Abort(
      "Error! Number of fuel species in input file must match "
      "SPRAY_FUEL_NUM");
  }

  std::vector<std::string> fuel_names;
  std::vector<std::string> dep_fuel_names;
  pele::physics::SprayProps::InitLiqProps<
    pele::physics::SprayProps::LiqPropType>
    init_liq_props;

  bool has_dep_spec = false;
  {
    pp.getarr("fuel_species", fuel_names);
    if (pp.contains("dep_fuel_species")) {
      has_dep_spec = true;
      pp.getarr("dep_fuel_species", dep_fuel_names);
    }

    // Read input parameters for liquid properties
    init_liq_props(&(m_sprayData->liqprops), fuel_names);

    // Set the fuel names
    for (int i = 0; i < nfuel; ++i) {
      m_sprayFuelNames[i] = fuel_names[i];
      if (has_dep_spec) {
        m_sprayDepNames[i] = dep_fuel_names[i];
      } else {
        m_sprayDepNames[i] = m_sprayFuelNames[i];
      }
    }
  }

  bool splash_model = false;
  int breakup_model = 0;
  //
  // Set the number of particles per parcel
  //
  pp.query("use_splash_model", splash_model);
  std::string breakup_model_str = "None";
  pp.query("use_breakup_model", breakup_model_str);
  if (breakup_model_str == "TAB") {
    breakup_model = 1;
    pp.query("max_parcel_size", m_maxNumPPP);
  } else if (breakup_model_str == "KHRT") {
    breakup_model = 2;
    pp.query("KHRT_B0", m_khrtB0);
    pp.query("KHRT_B1", m_khrtB1);
    pp.query("KHRT_C3", m_khrtC3);
  } else if (breakup_model_str == "None") {
    breakup_model = 0;
  } else {
    Abort(
      "'use_breakup_model' input not recognized. Must be 'TAB', 'KHRT', or "
      "'None'");
  }
  if (splash_model || (breakup_model > 0)) {
    pp.query("breakup_parcel_factor", m_breakupPPPFact);
    if (m_breakupPPPFact > 1. || m_breakupPPPFact < 0.) {
      Abort("'breakup_parcel_factor' must be between 0 and 1");
    }

    // Check proper input data for sigma and mu
    init_liq_props.init_breakupsplash(&(m_sprayData->liqprops), fuel_names);

    if (splash_model) {
      // TODO: Have this retrieved from proper boundary data
      pp.get("wall_temp", m_sprayData->wall_T);
      Real theta_c_deg = -1.;
      pp.get("contact_angle", theta_c_deg);
      if (theta_c_deg < 0. || theta_c_deg > 180.) {
        Abort("'contact_angle' must be between 0 and 180");
      }
      m_sprayData->theta_c = theta_c_deg * M_PI / 180.;
    }
    // Set the contact angle
    m_sprayData->do_splash = splash_model;
    m_sprayData->do_breakup = breakup_model;
  }

  // Set if spray ascii files should be written
  //
  pp.query("write_ascii_files", write_ascii_files);
  //
  // Set if gas phase spray source term should be written
  //
  pp.query("plot_src", plot_spray_src);
  //
  // Used in initData() on startup to read in a file of particles.
  //
  pp.query("init_file", spray_init_file);
#ifdef AMREX_USE_EB
  //
  // Spray source terms are only added to cells with a volume fraction higher
  // than this value
  //
  pp.query("min_eb_vfrac", m_sprayData->min_eb_vfrac);
#endif

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
    "wall_film_mass",  // Wall film mass
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
  }
  Gpu::streamSynchronize();
  ParallelDescriptor::Barrier();
}

void
SprayParticleContainer::spraySetup(
  const Real* body_force,
  pele::physics::eos::EosParm<pele::physics::PhysicsType::eos_type>* eosparms_h,
  const pele::physics::eos::EosParm<pele::physics::PhysicsType::eos_type>*
    eosparms_d)
{
#ifndef USE_MANIFOLD_EOS
#if NUM_SPECIES > 1
  Vector<std::string> spec_names;
  pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(
    spec_names, eosparms_h);

  for (int i = 0; i < SPRAY_FUEL_NUM; ++i) {
    for (int ns = 0; ns < NUM_SPECIES; ++ns) {
      std::string gas_spec = spec_names[ns];
      if (gas_spec == m_sprayFuelNames[i]) {
        m_sprayData->indx[i] = ns;
      }
      if (gas_spec == m_sprayDepNames[i]) {
        m_sprayData->dep_indx[i] = ns;
      }
    }
    if (m_sprayData->indx[i] < 0) {
      Abort("Fuel " + m_sprayFuelNames[i] + " not found in species list");
    }
    if (m_sprayData->dep_indx[i] < 0) {
      Abort("Fuel " + m_sprayDepNames[i] + " not found in species list");
    }
  }

  for (int i = 0; i < SPRAY_FUEL_NUM; ++i) {
    for (int j = i + 1; j < SPRAY_FUEL_NUM; ++j) {
      if (m_sprayData->dep_indx[i] == m_sprayData->dep_indx[j]) {
        m_sprayData->liquid_spec_share_gas_dep = true;
      }
    }
  }
#else
  m_sprayData->indx[0] = 0;
  m_sprayData->dep_indx[0] = 0;
#endif

  auto eos = pele::physics::PhysicsType::eos(eosparms_h);
  amrex::GpuArray<amrex::Real, NUM_SPECIES> mw;
  Vector<Real> fuelEnth(NUM_SPECIES);
  eos.molecular_weight(mw.data());
  m_sprayData->liqprops.init_mw(mw, m_sprayData->indx.data());
  eos.T2Hi(m_sprayData->liqprops.ref_T, fuelEnth.data());
  for (int ns = 0; ns < SPRAY_FUEL_NUM; ++ns) {
    const int fspec = m_sprayData->indx[ns];
    m_sprayData->liqprops.latentRef_minus_gasRefH_i[ns] =
      m_sprayData->liqprops.latent[ns] - fuelEnth[fspec] * SprayUnits::eng_conv;
  }

#else // USE_MANIFOLD_EOS is defined
  // Verify EOS can give molecular weights
  if (!eosparms_h->has_species_mw) {
    amrex::Error(
      "SpraySetup: Manifold EOS must contains spec molecular weights for "
      "Spray");
  }

  Vector<std::string> chemspec_names, manivar_names;
  // Manifold: For now, we require that each liquid/spray species
  // is cacuable from the Manifold model. We also require that
  // each species contributes to exactly one manifold variable with weight 1
  // which is specified through the "DepNames"
  std::set<std::string> unique_dep_names(
    m_sprayDepNames, m_sprayDepNames + SPRAY_FUEL_NUM);
  if (unique_dep_names.size() != SPRAY_FUEL_NUM) {
    amrex::Abort(
      "Each liquid spray species must uniquely contribute to one manifold "
      "parameter, as specified through dep_fuel_species");
  }

  pele::physics::eos::chemSpeciesNames<pele::physics::PhysicsType::eos_type>(
    chemspec_names, eosparms_h);
  pele::physics::eos::speciesNames<pele::physics::PhysicsType::eos_type>(
    manivar_names, eosparms_h);

  for (int i = 0; i < SPRAY_FUEL_NUM; ++i) {
    amrex::Print() << "\n SpraySpec = " << m_sprayFuelNames[i] << " "
                   << " dep species = " << m_sprayDepNames[i];
    for (int ns = 0; ns < chemspec_names.size(); ++ns) {
      amrex::Print() << "\n ChemSpec from manifold = " << chemspec_names[ns];
      std::string gas_spec = chemspec_names[ns];
      if (gas_spec == m_sprayFuelNames[i]) {
        m_sprayData->indx[i] = ns;
      }
    }
    if (m_sprayData->indx[i] < 0) {
      Abort(
        "Fuel " + m_sprayFuelNames[i] +
        " not found in species available in the manifold");
    }
    for (int ns = 0; ns < manivar_names.size(); ++ns) {
      std::string gas_spec = manivar_names[ns];
      amrex::Print() << "\n Manivar = " << manivar_names[ns];
      if (gas_spec == m_sprayDepNames[i]) {
        m_sprayData->dep_indx[i] = ns;
      }
    }
    if (m_sprayData->dep_indx[i] < 0) {
      Abort(
        "dep_fuel_species " + m_sprayDepNames[i] +
        " not found as a manifold parameter");
    }
  }

  // Initialize MW for spray model
  amrex::Vector<amrex::Real> mw(chemspec_names.size());
  auto eos = pele::physics::PhysicsType::eos(eosparms_h);
  eos.molecular_weight(mw.data());
  m_sprayData->liqprops.init_mw(mw.data(), m_sprayData->indx.data());

  // TODO: Handle latent heat for Manifold EOS
#endif

  // Stuff for both detailed chem and manifold
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    m_sprayData->body_force[dir] = body_force[dir];
  }
  m_sprayData->eosparm = eosparms_d;
  Gpu::copy(Gpu::hostToDevice, m_sprayData, m_sprayData + 1, d_sprayData);
  m_sprayData->eosparm = eosparms_h;
  Gpu::streamSynchronize();
  ParallelDescriptor::Barrier();
}

void
SprayParticleContainer::SprayInitialize(const std::string& restart_dir)
{
  bool init_sprays = false;
  if (restart_dir.empty() && spray_init_file.empty()) {
    init_sprays = true;
  }
  InitSprayParticles(init_sprays);
  if (!spray_init_file.empty()) {
    InitFromAsciiFile(spray_init_file, NSR_SPR);
  } else if (!restart_dir.empty()) {
    Restart(restart_dir, "particles");
  }
  PostInitRestart(restart_dir);
}
