
#include "SprayParticles.H"

using namespace amrex;

void
SprayParticleContainer::SprayParticleIO(
  const int level,
  const bool is_checkpoint,
  const std::string& dir)
{
  Vector<std::string> real_comp_names(NSR_SPR);
  AMREX_D_TERM(real_comp_names[SprayComps::pstateVel] = "xvel";
               , real_comp_names[SprayComps::pstateVel + 1] = "yvel";
               , real_comp_names[SprayComps::pstateVel + 2] = "zvel";);
  real_comp_names[SprayComps::pstateT] = "temperature";
  real_comp_names[SprayComps::pstateDia] = "diam";
  for (int sp = 0; sp < SPRAY_FUEL_NUM; ++sp) {
    real_comp_names[SprayComps::pstateY + sp] =
      "spray_mf_" + m_sprayFuelNames[sp];
  }
  real_comp_names[SprayComps::pstateNumDens] = "number_density";
  if (m_sprayData->do_breakup == 1) {
    real_comp_names[SprayComps::pstateBM1] = "Y_TAB";
    real_comp_names[SprayComps::pstateBM2] = "Ydot_TAB";
  } else if (m_sprayData->do_breakup == 2) {
    real_comp_names[SprayComps::pstateBM1] = "d0";
    real_comp_names[SprayComps::pstateBM2] = "rt_time";
  } else {
    real_comp_names[SprayComps::pstateBM1] = "unused1";
    real_comp_names[SprayComps::pstateBM2] = "unused2";
  }
  real_comp_names[SprayComps::pstateFilmHght] = "wall_film_height";
  Vector<std::string> int_comp_names;
  Checkpoint(dir, "particles", is_checkpoint, real_comp_names, int_comp_names);
  // Here we write ascii information every time we write a plot file
  if (level == 0 && SprayParticleContainer::write_ascii_files) {
    size_t num_end_loc = dir.find_last_of("0123456789") + 1;
    // Remove anything following numbers, like .temp
    std::string dirout = dir.substr(0, num_end_loc);
    size_t num_start_loc = dirout.find_last_not_of("0123456789") + 1;
    std::string numstring = dirout.substr(num_start_loc, dirout.length());
    std::string dir_path = dir;
    size_t num_end_path = dir_path.find_last_of("/") + 1;
    std::string part_dir_path = dir_path.substr(0, num_end_path);
    std::string fname = part_dir_path + "spray" + numstring + ".p3d";
    WriteAsciiFile(fname);
  }
  // Since injection can occur over multiple time steps, we must write the
  // current status of each jet in a checkpoint to ensure injection isn't
  // interrupted during restart
  // File line 0: Number of jets.
  // Each line after lists the jet name then the oustanding mass and injection
  // time, then the minimum injection parcel
  if (is_checkpoint && !m_sprayJets.empty()) {
    int numjets = static_cast<int>(m_sprayJets.size());
    std::string filename = dir + "/particles/injection_data.log";
    if (ParallelDescriptor::IOProcessor()) {
      std::ofstream file;
      file.open(filename.c_str(), std::ios::out | std::ios::trunc);
      if (!file.good()) {
        FileOpenFailed(filename);
      }
      file << numjets << "\n";
      file.flush();
      file.close();
      if (!file.good()) {
        Abort("Problem writing injection file");
      }
    }
    ParallelDescriptor::Barrier();
    // Each jet contains name, injection number density, outstanding mass and
    // time, and minimum parcels for injection
    for (int jindx = 0; jindx < numjets; ++jindx) {
      SprayJet* js = m_sprayJets[jindx].get();
      int jet_proc = js->Proc();
      if (ParallelDescriptor::MyProc() == jet_proc) {
        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
        std::ofstream file;
        file.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
        file.open(filename.c_str(), std::ios::out | std::ios::app);
        file.precision(15);
        if (!file.good()) {
          FileOpenFailed(filename);
        }
        file << js->jet_name() << " " << js->num_ppp() << " "
             << js->m_sumInjMass << " " << js->m_sumInjTime << " "
             << js->m_minParcel << "\n";
        file.flush();
        file.close();
        if (!file.good()) {
          Abort("Problem writing injection file");
        }
      }
      ParallelDescriptor::Barrier();
    }
  }
}

void
SprayParticleContainer::PostInitRestart(const std::string& dir)
{
  std::string JetDataFileName = dir + "/particles/injection_data.log";
  if (!m_sprayJets.empty()) {
    int numjets = static_cast<int>(m_sprayJets.size());
    // Ensure all jet names are unique
    std::vector<std::string> jet_names(numjets);
    for (int i = 0; i < numjets; ++i) {
      jet_names[i] = m_sprayJets[i]->jet_name();
      if (jet_names[i].empty()) {
        Abort("Jet not provided name");
      }
    }
    auto count = std::distance(
      jet_names.begin(), std::unique(jet_names.begin(), jet_names.end()));
    if (static_cast<int>(count) != numjets) {
      Abort("Duplicate jet names detected");
    }
    if (FileSystem::Exists(JetDataFileName)) {
      Vector<char> fileCharPtr;
      ParallelDescriptor::ReadAndBcastFile(JetDataFileName, fileCharPtr);
      std::string fileCharPtrString(fileCharPtr.dataPtr());
      std::istringstream JetDataFile(fileCharPtrString, std::istringstream::in);
      int in_numjets;
      JetDataFile >> in_numjets;
      Vector<std::string> in_jet_names(in_numjets);
      Vector<Real> in_inj_mass(in_numjets);
      Vector<Real> in_inj_time(in_numjets);
      Vector<Real> in_min_parcel(in_numjets);
      for (int i = 0; i < in_numjets; ++i) {
        JetDataFile >> in_jet_names[i] >> in_inj_mass[i] >> in_inj_time[i] >>
          in_min_parcel[i];
      }
      for (int ijets = 0; ijets < in_numjets; ++ijets) {
        std::string in_name = in_jet_names[ijets];
        for (int mjets = 0; mjets < numjets; ++mjets) {
          SprayJet* js = m_sprayJets[mjets].get();
          if (js->jet_name() == in_name) {
            js->m_sumInjMass = in_inj_mass[ijets];
            js->m_sumInjTime = in_inj_time[ijets];
            js->m_minParcel = in_min_parcel[ijets];
          }
        }
      }
    }
  } else if (FileSystem::Exists(JetDataFileName)) {
    if (ParallelDescriptor::IOProcessor()) {
      Print() << "Warning: Restart file contains jet information but no "
                 "SprayJets have been initialized\n";
    }
  }
  Gpu::streamSynchronize();
}
