
#include "SprayParticles.H"

using namespace amrex;

void
SprayParticleContainer::SprayParticleIO(
  const int level,
  const bool is_checkpoint,
  const int write_ascii,
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
      "spray_mf_" + spray_fuel_names[sp];
  }
  Vector<std::string> int_comp_names;
  Checkpoint(dir, "particles", is_checkpoint, real_comp_names, int_comp_names);
  // Here we write ascii information every time we write a plot file
  if (level == 0 && write_ascii == 1) {
    // TODO: Would be nice to be able to use file_name_digits
    // instead of doing this
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
  if (
    is_checkpoint && m_sprayJets.size() > 0 &&
    ParallelDescriptor::IOProcessor()) {
    std::string filename = dir + "/particles/injection_data.log";
    std::ofstream file;
    file.open(filename.c_str(), std::ios::out | std::ios::trunc);
    if (!file.good()) {
      FileOpenFailed(filename);
    }
    int numjets = static_cast<int>(m_sprayJets.size());
    file << numjets << "\n";
    for (int jindx = 0; jindx < numjets; ++jindx) {
      SprayJet* js = m_sprayJets[jindx].get();
      file << js->jet_name() << " " << js->m_sumInjMass << " "
           << js->m_sumInjTime << " " << js->m_minParcel << "\n";
    }
    file.flush();
    file.close();
    if (!file.good()) {
      Abort("Problem writing injection file");
    }
  }
}

void
SprayParticleContainer::PostInitRestart(const std::string& dir)
{
  int numjets = static_cast<int>(m_sprayJets.size());
  // Ensure all jet names are unique
  if (numjets > 0) {
    std::vector<std::string> jet_names(numjets);
    for (int i = 0; i < numjets; ++i) {
      jet_names[i] = m_sprayJets[i]->jet_name();
      if (jet_names[i] == "") {
        Abort("Jet not provided proper name");
      }
    }
    int count = std::distance(
      jet_names.begin(), std::unique(jet_names.begin(), jet_names.end()));
    if (count > 0) {
      Abort("Duplicate jet names detected");
    }
  }
  std::string JetDataFileName = dir + "/particles/injection_data.log";
  if (FileSystem::Exists(JetDataFileName)) {
    if (numjets == 0 && ParallelDescriptor::IOProcessor()) {
      Print() << "Warning: Restart file contains jet information but no "
                 "SprayJets have been initialized"
              << std::endl;
      return;
    }
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
}
