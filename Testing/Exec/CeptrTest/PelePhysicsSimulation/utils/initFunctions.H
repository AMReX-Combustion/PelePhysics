void
parseInput(
  amrex::ParmParse& pp,
  amrex::ParmParse& ppode,
  std::string& fuel_name,
  std::string& chem_integrator,
  bool& do_plt,
  std::string& pltfile,
  int& initFromChk,
  std::string& chkfile,
  std::string& reactFormat,
  int& reactFunc,
  int& ode_ncells,
  amrex::Real& dt,
  int& ndt,
  int& ode_iE,
  amrex::Real& rtol,
  amrex::Real& atol,
  int& use_typ_vals,
  std::array<int, 3>& ncells,
  int& max_grid_size)
{
  pp.get("fuel_name", fuel_name);
  chem_integrator = "";
  pp.get("chem_integrator", chem_integrator);
  do_plt = false;
  if (pp.countval("plotfile") > 0) {
    pp.get("plotfile", pltfile);
    do_plt = true;
  }
  /* initial data */
  // Can either generate single level data
  // or read in data from a checkpoint like file
  initFromChk = 0;
  pp.query("initFromFile", initFromChk);
  chkfile = "";
  if (initFromChk) {
    pp.query("initFile", chkfile);
  }
  /* react() function version */
  // 1 -> Array4 version of react()  (Default)
  // 2 -> 1d raw pointer version of react()
  reactFormat = "Array4";
  pp.query("reactFormat", reactFormat);
  if (reactFormat == "Array4") {
    reactFunc = 1;
  } else if (reactFormat == "1dArray") {
    reactFunc = 2;
  } else {
    amrex::Abort(" --> reactFormat can only be 'Array4' or '1dArray' !");
  }

  /* ODE inputs */
  ode_ncells = 1;
  ppode.query("ode_ncells", ode_ncells); // number of cells to integrate per
                                         // call, used only if reactFunc = 2

  dt = 1.e-5;
  ppode.query("dt", dt);

  ndt = 1;
  ppode.query("ndt", ndt); // number of solver calls per dt

  ode_iE = 2;
  ppode.query(
    "reactor_type",
    ode_iE); // RHS type, 1: e (PeleC), !1: h (PeleLM)  <------ FIXME!

  rtol = 1e-10;
  ppode.query("rtol", rtol);

  atol = 1e-10;
  ppode.query("atol", atol);

  use_typ_vals = 0;
  ppode.query("use_typ_vals", use_typ_vals);

  amrex::Print() << "ODE solver: " << chem_integrator << std::endl;
  amrex::Print() << "Type of reactor: "
                 << (ode_iE == 1 ? "e (PeleC)" : "h (PeleLM)")
                 << std::endl; // <---- FIXME
  amrex::Print() << "Fuel: " << fuel_name << ", Oxy: O2" << std::endl;

  // ~~~~ Geometry
  ncells = {AMREX_D_DECL(1, 1, 1)};
  if (pp.countval("ncells") == 1) {
    pp.get("ncells", ncells[0]);
    ncells = {AMREX_D_DECL(ncells[0], 1, 1)};
  } else if (pp.countval("ncells") >= AMREX_SPACEDIM) {
    amrex::Vector<int> nc(AMREX_SPACEDIM);
    pp.getarr("ncells", nc, 0, AMREX_SPACEDIM);
    ncells = {AMREX_D_DECL(nc[0], nc[1], nc[2])};
  } else {
    amrex::Abort("ncells has to have length 1 or spacedim");
  }
  max_grid_size = 16;
  pp.query("max_grid_size", max_grid_size);
}

void
getFuelID(std::string& fuel_name, int& fuel_idx)
{

  /* Mixture info */
  fuel_idx = -1;
  if (fuel_name == "H2") {
    fuel_idx = H2_ID;
#ifdef CH4_ID
  } else if (fuel_name == "CH4") {
    fuel_idx = CH4_ID;
#endif
#ifdef NC12H26_ID
  } else if (fuel_name == "NC12H26") {
    fuel_idx = NC12H26_ID;
#endif
#ifdef IC8H18_ID
  } else if (fuel_name == "IC8H18") {
    fuel_idx = IC8H18_ID;
#endif
  }
}

void
initializeGeom(
  amrex::Vector<amrex::Geometry>& geoms,
  amrex::Vector<amrex::BoxArray>& grids,
  amrex::Vector<amrex::DistributionMapping>& dmaps,
  int& finest_level,
  std::array<int, 3>& ncells,
  int& ndt,
  amrex::Real& dt,
  int& max_grid_size)
{
  // -----------------------------------------------------------------------------
  // Initialize geom/data
  // When initFromChk = 0, default is single level data provided by
  // initialize_data with a -1:1 unit length realbox in each dir and cell
  // count provided by the user otherwise number of levels, geom and all are
  // read in from the chkfile
  // -----------------------------------------------------------------------------

  // -----------------------------------------------------------------------------
  // Resize vectors
  // -----------------------------------------------------------------------------
  geoms.resize(finest_level + 1);
  grids.resize(finest_level + 1);
  dmaps.resize(finest_level + 1);

  amrex::Box domain(
    amrex::IntVect(AMREX_D_DECL(0, 0, 0)),
    amrex::IntVect(AMREX_D_DECL(ncells[0] - 1, ncells[1] - 1, ncells[2] - 1)));

  amrex::RealBox real_box(
    {AMREX_D_DECL(-1.0, -1.0, -1.0)}, {AMREX_D_DECL(1.0, 1.0, 1.0)});

  int coord = 0;

  amrex::Array<int, AMREX_SPACEDIM> is_periodic{AMREX_D_DECL(1, 1, 1)};

  geoms[0] = amrex::Geometry(domain, real_box, coord, is_periodic);

  amrex::Print() << "Integrating " << domain.numPts() << " cells for: " << dt
                 << " seconds with " << ndt << " substeps \n";

  // -----------------------------------------------------------------------------
  // Define amrex::BoxArray / Dmap
  // -----------------------------------------------------------------------------
  grids[0] = amrex::BoxArray(domain);
  grids[0].maxSize(max_grid_size);
  dmaps[0] =
    amrex::DistributionMapping(grids[0], amrex::ParallelDescriptor::NProcs());
}

void
initializeData(
  int& num_grow,
  amrex::Vector<amrex::MultiFab>& mf,
  amrex::Vector<amrex::MultiFab>& rY_source_ext,
  amrex::Vector<amrex::MultiFab>& mfE,
  amrex::Vector<amrex::MultiFab>& rY_source_energy_ext,
  amrex::Vector<amrex::MultiFab>& fctCount,
  amrex::Vector<amrex::iMultiFab>& dummyMask,
  int& finest_level,
  amrex::Vector<amrex::Geometry>& geoms,
  amrex::Vector<amrex::BoxArray>& grids,
  amrex::Vector<amrex::DistributionMapping>& dmaps,
  int& fuel_idx,
  int& ode_iE)
{

  // -----------------------------------------------------------------------------
  // Create MFs and generate initial data
  // -----------------------------------------------------------------------------
  for (int lev = 0; lev <= finest_level; ++lev) {
    mf[lev].define(grids[lev], dmaps[lev], NUM_SPECIES + 1, num_grow);
    rY_source_ext[lev].define(grids[lev], dmaps[lev], NUM_SPECIES, num_grow);
    mfE[lev].define(grids[lev], dmaps[lev], 1, num_grow);
    rY_source_energy_ext[lev].define(grids[lev], dmaps[lev], 1, num_grow);
    fctCount[lev].define(grids[lev], dmaps[lev], 1, num_grow);
    dummyMask[lev].define(grids[lev], dmaps[lev], 1, num_grow);
    dummyMask[lev].setVal(1);
  }

  amrex::FabArrayBase::mfiter_tile_size =
    amrex::IntVect(AMREX_D_DECL(1024, 1024, 1024));

  // -----------------------------------------------------------------------------
  // Initialize data
  // -----------------------------------------------------------------------------
  for (int lev = 0; lev <= finest_level; ++lev) {
    const auto geomdata = geoms[lev].data();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(mf[lev], amrex::TilingIfNotGPU()); mfi.isValid();
         ++mfi) {

      const amrex::Box& box = mfi.tilebox();

      amrex::Array4<amrex::Real> const& rY_a = mf[lev].array(mfi);
      amrex::Array4<amrex::Real> const& rYs_a = rY_source_ext[lev].array(mfi);
      amrex::Array4<amrex::Real> const& E_a = mfE[lev].array(mfi);
      amrex::Array4<amrex::Real> const& rE_a =
        rY_source_energy_ext[lev].array(mfi);
      amrex::ParallelFor(
        box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
          initialize_data(
            i, j, k, fuel_idx, ode_iE, rY_a, rYs_a, E_a, rE_a, geomdata);
        });
    }
  }
}
