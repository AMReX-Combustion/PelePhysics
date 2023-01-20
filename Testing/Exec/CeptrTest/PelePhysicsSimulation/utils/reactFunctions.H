void
doReact_ode_iE1(
  int& lev,
  amrex::Real& dt,
  int& ndt,
  int& omp_thread,
  amrex::MFIter& mfi,
  amrex::Vector<amrex::MultiFab>& mf,
  amrex::Vector<amrex::MultiFab>& rY_source_ext,
  amrex::Vector<amrex::MultiFab>& mfE,
  amrex::Vector<amrex::MultiFab>& rY_source_energy_ext,
  amrex::Vector<amrex::MultiFab>& fctCount,
  amrex::Vector<amrex::iMultiFab>& dummyMask,
  std::unique_ptr<pele::physics::reactions::ReactorBase>& reactor,
  pele::physics::transport::TransportParams<
    pele::physics::PhysicsType::transport_type>& trans_parms)
{

  BL_PROFILE_VAR_NS("React", ReactInLoop);
  const amrex::Box& box = mfi.tilebox();
  const int nc = box.numPts();
  auto const& rhoY = mf[lev].array(mfi);
  auto const& T = mf[lev].array(mfi, NUM_SPECIES);
  auto const& rhoE = mfE[lev].array(mfi);
  auto const& frcExt = rY_source_ext[lev].array(mfi);
  auto const& frcEExt = rY_source_energy_ext[lev].array(mfi);
  auto const& fc = fctCount[lev].array(mfi);
  auto const& mask = dummyMask[lev].array(mfi);
  amrex::Real time = 0.0;
  amrex::Real dt_lev = dt / std::pow(2, lev);
  amrex::Real dt_incr = dt_lev / ndt;
  int tmp_fc;
  if (omp_thread == 0) {
    amrex::Print() << "  [" << lev << "]"
                   << " integrating " << nc << " cells" << std::endl;
  }

  // transport needs
  auto const* ltransparm = trans_parms.device_trans_parm();
  amrex::FArrayBox rho_fab(box, 1, amrex::The_Async_Arena());
  auto const& rho_a = rho_fab.array();
  amrex::FArrayBox Y_fab(box, NUM_SPECIES, amrex::The_Async_Arena());
  auto const& Y_a = Y_fab.array();
  amrex::FArrayBox D_fab(box, NUM_SPECIES, amrex::The_Async_Arena());
  auto const& D_a = D_fab.array();
  amrex::FArrayBox mu_fab(box, 1, amrex::The_Async_Arena());
  auto const& mu_a = mu_fab.array();
  amrex::FArrayBox xi_fab(box, 1, amrex::The_Async_Arena());
  auto const& xi_a = xi_fab.array();
  amrex::FArrayBox lam_fab(box, 1, amrex::The_Async_Arena());
  auto const& lam_a = lam_fab.array();
  amrex::FArrayBox chi_fab(box, NUM_SPECIES, amrex::The_Async_Arena());
  auto const& chi_a = chi_fab.array();
  // Transport
  amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    rho_a(i, j, k) = 0.0;
    for (int n = 0; n < NUM_SPECIES; n++) {
      rho_a(i, j, k) += rhoY(i, j, k, n);
    }
    for (int n = 0; n < NUM_SPECIES; n++) {
      Y_a(i, j, k, n) = rhoY(i, j, k, n) / rho_a(i, j, k);
    }
  });
  amrex::launch(box, [=] AMREX_GPU_DEVICE(amrex::Box const& tbx) {
    auto trans = pele::physics::PhysicsType::transport();
    trans.get_transport_coeffs(
      tbx, Y_a, T, rho_a, D_a, chi_a, mu_a, xi_a, lam_a, ltransparm);
  });

  /* Solve */
  BL_PROFILE_VAR_START(ReactInLoop);
  std::ofstream myfile;
  myfile.open("results.txt");
  myfile << "time,temperature,density,viscosity" << std::endl;
  myfile << 0 << "," << T.dataPtr()[0] << "," << rho_a.dataPtr()[0] << ","
         << mu_a.dataPtr()[0] << std::endl;
  double timeElapsed = 0;
  for (int ii = 0; ii < ndt; ++ii) {
    tmp_fc = reactor->react(
      box, rhoY, frcExt, T, rhoE, frcEExt, fc, mask, dt_incr, time
#ifdef AMREX_USE_GPU
      ,
      amrex::Gpu::gpuStream()
#endif
    );

    // Transport
    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      rho_a(i, j, k) = 0.0;
      for (int n = 0; n < NUM_SPECIES; n++) {
        rho_a(i, j, k) += rhoY(i, j, k, n);
      }
      for (int n = 0; n < NUM_SPECIES; n++) {
        Y_a(i, j, k, n) = rhoY(i, j, k, n) / rho_a(i, j, k);
      }
    });
    amrex::launch(box, [=] AMREX_GPU_DEVICE(amrex::Box const& tbx) {
      auto trans = pele::physics::PhysicsType::transport();
      trans.get_transport_coeffs(
        tbx, Y_a, T, rho_a, D_a, chi_a, mu_a, xi_a, lam_a, ltransparm);
    });

    dt_incr = dt_lev / ndt;
    timeElapsed += dt_incr;
    amrex::Gpu::Device::streamSynchronize();
    myfile << timeElapsed << "," << T.dataPtr()[0] << "," << rho_a.dataPtr()[0]
           << "," << mu_a.dataPtr()[0] << std::endl;
  }
  BL_PROFILE_VAR_STOP(ReactInLoop);
}

void
doReact_ode_iE2(
  int& lev,
  amrex::Real& dt,
  int& ndt,
  int& omp_thread,
  int& ode_ncells,
  amrex::MFIter& mfi,
  amrex::Vector<amrex::MultiFab>& mf,
  amrex::Vector<amrex::MultiFab>& rY_source_ext,
  amrex::Vector<amrex::MultiFab>& mfE,
  amrex::Vector<amrex::MultiFab>& rY_source_energy_ext,
  amrex::Vector<amrex::MultiFab>& fctCount,
  amrex::Vector<amrex::iMultiFab>& dummyMask,
  std::unique_ptr<pele::physics::reactions::ReactorBase>& reactor)
{

  BL_PROFILE_VAR_NS("Allocs", Allocs);
  BL_PROFILE_VAR_NS("Flatten", mainflatten);
  BL_PROFILE_VAR_NS("React", ReactInLoop);
  int extra_cells = 0;
  const amrex::Box& box = mfi.tilebox();
  const int nc = box.numPts();
  auto const& rhoY = mf[lev].array(mfi);
  auto const& T = mf[lev].array(mfi, NUM_SPECIES);
  auto const& rhoE = mfE[lev].array(mfi);
  auto const& frcExt = rY_source_ext[lev].array(mfi);
  auto const& frcEExt = rY_source_energy_ext[lev].array(mfi);
  auto const& fc = fctCount[lev].array(mfi);
  auto const& mask = dummyMask[lev].array(mfi);

#ifdef AMREX_USE_GPU
  ode_ncells = nc;
#endif
  extra_cells = nc - (nc / ode_ncells) * ode_ncells;
  if (omp_thread == 0) {
    amrex::Print() << " Integrating " << nc << " cells with a " << ode_ncells
                   << " ode cell buffer ";
    amrex::Print() << "(" << extra_cells << " extra cells) \n";
  }

  BL_PROFILE_VAR_START(Allocs);
  int nCells = nc + extra_cells;

#ifdef AMREX_USE_GPU
  auto tmp_vect_d = (amrex::Real*)amrex::The_Device_Arena()->alloc(
    nCells * (NUM_SPECIES + 1) * sizeof(amrex::Real));
  auto tmp_src_vect_d = (amrex::Real*)amrex::The_Device_Arena()->alloc(
    nCells * NUM_SPECIES * sizeof(amrex::Real));
  auto tmp_vect_energy_d = (amrex::Real*)amrex::The_Device_Arena()->alloc(
    nCells * sizeof(amrex::Real));
  auto tmp_src_vect_energy_d = (amrex::Real*)amrex::The_Device_Arena()->alloc(
    nCells * sizeof(amrex::Real));
  auto tmp_fc_d =
    (long int*)amrex::The_Device_Arena()->alloc(nCells * sizeof(amrex::Real));
  auto tmp_mask_d = (amrex::Real*)amrex::The_Device_Arena()->alloc(
    nCells * sizeof(amrex::Real));
#endif

  auto tmp_vect = new amrex::Real[nCells * (NUM_SPECIES + 1)];
  auto tmp_src_vect = new amrex::Real[nCells * NUM_SPECIES];
  auto tmp_vect_energy = new amrex::Real[nCells];
  auto tmp_src_vect_energy = new amrex::Real[nCells];
  auto tmp_fc = new long int[nCells];
  auto tmp_mask = new int[nCells];

  BL_PROFILE_VAR_STOP(Allocs);

  BL_PROFILE_VAR_START(mainflatten);
#ifndef AMREX_USE_GPU
  reactor->flatten(
    box, nCells, rhoY, frcExt, T, rhoE, frcEExt, tmp_vect, tmp_src_vect,
    tmp_vect_energy, tmp_src_vect_energy);

  for (int icell = nc; icell < nc + extra_cells; icell++) {
    for (int sp = 0; sp < NUM_SPECIES; sp++) {
      tmp_vect[icell * (NUM_SPECIES + 1) + sp] = rhoY(0, 0, 0, sp);
      tmp_src_vect[icell * NUM_SPECIES + sp] = frcExt(0, 0, 0, sp);
    }
    tmp_vect[icell * (NUM_SPECIES + 1) + NUM_SPECIES] = T(0, 0, 0);
    tmp_vect_energy[icell] = rhoE(0, 0, 0);
    tmp_src_vect_energy[icell] = frcEExt(0, 0, 0);
    tmp_mask[icell] = mask(0, 0, 0);
  }
#else
  reactor->flatten(
    box, nCells, rhoY, frcExt, T, rhoE, frcEExt, tmp_vect_d, tmp_src_vect_d,
    tmp_vect_energy_d, tmp_src_vect_energy_d);

  amrex::Gpu::copy(
    amrex::Gpu::deviceToHost, tmp_vect_d,
    tmp_vect_d + nCells * (NUM_SPECIES + 1), tmp_vect);
  amrex::Gpu::copy(
    amrex::Gpu::deviceToHost, tmp_src_vect_d,
    tmp_src_vect_d + nCells * NUM_SPECIES, tmp_src_vect);
  amrex::Gpu::copy(
    amrex::Gpu::deviceToHost, tmp_vect_energy_d, tmp_vect_energy_d + nCells,
    tmp_vect_energy);
  amrex::Gpu::copy(
    amrex::Gpu::deviceToHost, tmp_src_vect_energy_d,
    tmp_src_vect_energy_d + nCells, tmp_src_vect_energy);
#endif
  BL_PROFILE_VAR_STOP(mainflatten);

  /* Solve */
  BL_PROFILE_VAR_START(ReactInLoop);
  for (int i = 0; i < nCells; i += ode_ncells) {
    tmp_fc[i] = 0;
    amrex::Real time = 0.0;
    amrex::Real dt_lev = dt / std::pow(2, lev);
    amrex::Real dt_incr = dt_lev / ndt;
    for (int ii = 0; ii < ndt; ++ii) {
      tmp_fc[i] += reactor->react(
        &tmp_vect[i * (NUM_SPECIES + 1)], &tmp_src_vect[i * NUM_SPECIES],
        &tmp_vect_energy[i], &tmp_src_vect_energy[i], dt_incr, time, ode_ncells
#ifdef AMREX_USE_GPU
        ,
        amrex::Gpu::gpuStream()
#endif
      );

      dt_incr = dt_lev / ndt;
      for (int ic = i + 1; ic < i + ode_ncells; ++ic) {
        tmp_fc[ic] = tmp_fc[i];
      }
      amrex::Gpu::Device::streamSynchronize();
    }
  }
  BL_PROFILE_VAR_STOP(ReactInLoop);

  BL_PROFILE_VAR_START(mainflatten);
#ifndef AMREX_USE_GPU
  reactor->unflatten(
    box, nCells, rhoY, T, rhoE, frcEExt, fc, tmp_vect, tmp_vect_energy, tmp_fc,
    dt);
#else

  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, tmp_vect, tmp_vect + nCells * (NUM_SPECIES + 1),
    tmp_vect_d);
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, tmp_src_vect, tmp_src_vect + nCells * NUM_SPECIES,
    tmp_src_vect_d);
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, tmp_vect_energy, tmp_vect_energy + nCells,
    tmp_vect_energy_d);
  amrex::Gpu::copy(
    amrex::Gpu::hostToDevice, tmp_src_vect_energy, tmp_src_vect_energy + nCells,
    tmp_src_vect_energy_d);
  amrex::Gpu::copy(amrex::Gpu::hostToDevice, tmp_fc, tmp_fc + nCells, tmp_fc_d);

  reactor->unflatten(
    box, nCells, rhoY, T, rhoE, frcEExt, fc, tmp_vect_d, tmp_vect_energy_d,
    tmp_fc_d, dt);
#endif
  BL_PROFILE_VAR_STOP(mainflatten);

  delete[] tmp_vect;
  delete[] tmp_src_vect;
  delete[] tmp_vect_energy;
  delete[] tmp_src_vect_energy;
  delete[] tmp_fc;
  delete[] tmp_mask;
#ifdef AMREX_USE_GPU
  amrex::The_Device_Arena()->free(tmp_vect_d);
  amrex::The_Device_Arena()->free(tmp_src_vect_d);
  amrex::The_Device_Arena()->free(tmp_vect_energy_d);
  amrex::The_Device_Arena()->free(tmp_src_vect_energy_d);
  amrex::The_Device_Arena()->free(tmp_fc_d);
  amrex::The_Device_Arena()->free(tmp_mask_d);
#endif
}
