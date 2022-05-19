void
plotResult(
  bool& do_plt,
  std::string& pltfile,
  int& finest_level,
  amrex::Vector<amrex::MultiFab>& mf,
  amrex::Vector<amrex::Geometry>& geoms)
{

  if (do_plt) {

    BL_PROFILE_VAR_START(PlotFile);
    std::string outfile = amrex::Concatenate(pltfile, 1);
    // TODO: add fct count to this output
    amrex::Vector<int> isteps(finest_level + 1, 0);
    amrex::Vector<amrex::IntVect> refRatios(
      finest_level, {AMREX_D_DECL(2, 2, 2)});
    amrex::Vector<std::string> plt_VarsName;
    for (int k = 0; k < NUM_SPECIES; ++k) {
      plt_VarsName.push_back("SPEC" + std::to_string(k));
    }
    plt_VarsName.push_back("TEMP");

    amrex::WriteMultiLevelPlotfile(
      outfile, finest_level + 1, GetVecOfConstPtrs(mf), plt_VarsName, geoms,
      0.0, isteps, refRatios);
    BL_PROFILE_VAR_STOP(PlotFile);
  }
}
