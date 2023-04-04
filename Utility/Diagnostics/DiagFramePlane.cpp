#include "DiagFramePlane.H"
#include "AMReX_VisMF.H"
#include <AMReX_FPC.H>
#include "AMReX_PlotFileUtil.H"
#include <regex>
#include <stdio.h>

void
printLowerDimIntVect(
  std::ostream& a_File, const amrex::IntVect& a_IntVect, int skipDim)
{
  int doneDim = 0;
  a_File << '(';
  for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
    if (idim != skipDim) {
      a_File << a_IntVect[idim];
      doneDim++;
      if (doneDim < AMREX_SPACEDIM - 1)
        a_File << ",";
    }
  }
  a_File << ')';
}

void
printLowerDimBox(std::ostream& a_File, const amrex::Box& a_box, int skipDim)
{
  a_File << '(';
  printLowerDimIntVect(a_File, a_box.smallEnd(), skipDim);
  a_File << ' ';
  printLowerDimIntVect(a_File, a_box.bigEnd(), skipDim);
  a_File << ' ';
  printLowerDimIntVect(a_File, a_box.type(), skipDim);
  a_File << ')';
}

void
DiagFramePlane::init(const std::string& a_prefix, std::string_view a_diagName)
{
  DiagBase::init(a_prefix, a_diagName);

  if (m_filters.size() != 0) {
    amrex::Print() << " Filters are not available on DiagFramePlane and will "
                      "be discarded \n";
  }

  amrex::ParmParse pp(a_prefix);

  // Outputed variables
  int nOutFields = pp.countval("field_names");
  AMREX_ASSERT(nOutFields > 0);
  m_fieldNames.resize(nOutFields);
  m_fieldIndices_d.resize(nOutFields);
  for (int f{0}; f < nOutFields; ++f) {
    pp.get("field_names", m_fieldNames[f], f);
  }

  // Plane normal
  pp.get("normal", m_normal);
  AMREX_ASSERT(m_normal >= 0 && m_normal < AMREX_SPACEDIM);

  // Plane center
  amrex::Vector<amrex::Real> center;
  pp.getarr("center", center, 0, pp.countval("center"));
  if (center.size() == AMREX_SPACEDIM) {
    for (int idim = 0; idim < AMREX_SPACEDIM; idim++) {
      m_center[idim] = center[idim];
    }
  } else if (center.size() == 1) {
    m_center[m_normal] = center[0];
  }

  // Interpolation
  std::string intType = "Quadratic";
  pp.query("interpolation", intType);
  if (intType == "Linear") {
    m_interpType = Linear;
  } else if (intType == "Quadratic") {
    m_interpType = Quadratic;
  } else {
    amrex::Abort("Unknown interpolation type for " + a_prefix);
  }
}

void
DiagFramePlane::addVars(amrex::Vector<std::string>& a_varList)
{
  DiagBase::addVars(a_varList);
  for (const auto& v : m_fieldNames) {
    a_varList.push_back(v);
  }
}

void
DiagFramePlane::prepare(
  int a_nlevels,
  const amrex::Vector<amrex::Geometry>& a_geoms,
  const amrex::Vector<amrex::BoxArray>& a_grids,
  const amrex::Vector<amrex::DistributionMapping>& a_dmap,
  const amrex::Vector<std::string>& a_varNames)
{
  if (first_time) {
    // Store the level0 geometry
    auto initDomain = a_geoms[0].Domain();
    auto initRealBox = a_geoms[0].ProbDomain();
    const amrex::Real* dxlcl = a_geoms[0].CellSize();
    int cdim = 0;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
      if (idim != m_normal) {
        initDomain.setRange(
          cdim, initDomain.smallEnd(idim), initDomain.bigEnd(idim) + 1);
        initRealBox.setLo(cdim, a_geoms[0].ProbLo(idim));
        initRealBox.setHi(cdim, a_geoms[0].ProbHi(idim));
        cdim += 1;
      }
    }
    initDomain.setRange(2, 0, 1);
    initRealBox.setLo(2, 0.0);
    initRealBox.setHi(2, dxlcl[2]);
    m_geomLev0.define(
      initDomain, initRealBox, a_geoms[0].Coord(),
      amrex::Array<int, AMREX_SPACEDIM>({AMREX_D_DECL(0, 0, 0)}));

    int nOutFields = m_fieldIndices_d.size();
    amrex::Vector<int> m_fieldIndices(nOutFields, 0);
    for (int f{0}; f < nOutFields; ++f) {
      m_fieldIndices[f] = getFieldIndex(m_fieldNames[f], a_varNames);
    }
    amrex::Gpu::copy(
      amrex::Gpu::hostToDevice, m_fieldIndices.begin(), m_fieldIndices.end(),
      m_fieldIndices_d.begin());

    first_time = false;
  }

  // Resize internal vectors
  m_intwgt.resize(a_nlevels);
  m_k0.resize(a_nlevels);

  // On each level, find the k0 where the plane lays
  // and the weight of the directionnal interpolation
  for (int lev = 0; lev < a_nlevels; lev++) {
    const amrex::Real* dx = a_geoms[lev].CellSize();
    const amrex::Real* problo = a_geoms[lev].ProbLo();
    // How many dx away from the lowest cell-center ?
    amrex::Real dist =
      (m_center[m_normal] - (problo[m_normal] + 0.5 * dx[m_normal])) /
      dx[m_normal];
    int k0 = static_cast<int>(std::round(dist));
    dist -= static_cast<amrex::Real>(k0);
    m_k0[lev] = k0;
    if (m_interpType == Quadratic) {
      // Quadratic interp. weights on k0-1, k0, k0+1
      m_intwgt[lev][0] = 0.5 * (dist - 1.0) * (dist - 2.0);
      m_intwgt[lev][1] = dist * (2.0 - dist);
      m_intwgt[lev][2] = 0.5 * dist * (dist - 1.0);
      ;
    } else if (m_interpType == Linear) {
      // linear interp. weights on k0-1, k0, k0+1
      if (dist > 0.0) {
        m_intwgt[lev][0] = 0.0;
        m_intwgt[lev][1] = 1.0 - dist;
        m_intwgt[lev][2] = dist;
      } else if (dist < 0.0) {
        m_intwgt[lev][0] = -dist;
        m_intwgt[lev][1] = 1.0 + dist;
        m_intwgt[lev][2] = 0.0;
      } else {
        m_intwgt[lev][0] = 0.0;
        m_intwgt[lev][1] = 1.0;
        m_intwgt[lev][2] = 0.0;
      }
    }
  }

  // Assemble the 2D slice boxArray
  m_sliceBA.resize(a_nlevels);
  m_sliceDM.resize(a_nlevels);
  m_dmConvert.resize(a_nlevels);
  for (int lev = 0; lev < a_nlevels; lev++) {
    amrex::Vector<int> pmap;
    amrex::BoxList bl(a_grids[lev].ixType());
    bl.reserve(a_grids[lev].size());
    amrex::Vector<int> dmConvertLev;
    for (int i = 0; i < a_grids[lev].size(); ++i) {
      auto cBox = a_grids[lev][i];
      amrex::IntVect ploc(
        AMREX_D_DECL(cBox.smallEnd(0), cBox.smallEnd(1), cBox.smallEnd(2)));
      ploc[m_normal] = m_k0[lev];
      if (cBox.contains(ploc)) {
        amrex::Box zNormalBax;
        int idx = 0;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
          if (idim != m_normal) {
            zNormalBax.setRange(idx, cBox.smallEnd(idim), cBox.length(idim));
            idx++;
          }
        }
        zNormalBax.setRange(AMREX_SPACEDIM - 1, 0, 1);
        bl.push_back(zNormalBax);
        pmap.push_back(a_dmap[lev][i]);
        dmConvertLev.push_back(i);
      }
      m_dmConvert[lev] = dmConvertLev;
    }
    m_sliceBA[lev].define(bl);
    m_sliceDM[lev].define(pmap);
  }
}

void
DiagFramePlane::processDiag(
  int a_nstep,
  const amrex::Real& a_time,
  const amrex::Vector<const amrex::MultiFab*>& a_state,
  const amrex::Vector<std::string>& a_stateVar)
{
  // Interpolate data to slice
  amrex::Vector<amrex::MultiFab> planeData(a_state.size());
  for (int lev = 0; lev < a_state.size(); ++lev) {
    planeData[lev].define(
      m_sliceBA[lev], m_sliceDM[lev], m_fieldNames.size(), 0);
    int p0 = m_k0[lev];
#ifdef AMREX_USE_OMP
#pragma omp parallel if (amrex::Gpu::notInLaunchRegion())
#endif
    for (amrex::MFIter mfi(planeData[lev], amrex::TilingIfNotGPU());
         mfi.isValid(); ++mfi) {
      const auto& bx = mfi.tilebox();
      const int state_idx = m_dmConvert[lev][mfi.index()];
      auto const& state = a_state[lev]->const_array(state_idx, 0);
      auto const& plane = planeData[lev].array(mfi);
      auto const& intwgt = m_intwgt[lev];
      auto* idx_d_p = m_fieldIndices_d.dataPtr();
      if (m_normal == 0) {
        amrex::ParallelFor(
          bx, m_fieldNames.size(),
          [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
            int stIdx = idx_d_p[n];
            plane(i, j, k, n) = intwgt[0] * state(p0 - 1, i, j, stIdx) +
                                intwgt[1] * state(p0, i, j, stIdx) +
                                intwgt[2] * state(p0 + 1, i, j, stIdx);
          });
      } else if (m_normal == 1) {
        amrex::ParallelFor(
          bx, m_fieldNames.size(),
          [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
            int stIdx = idx_d_p[n];
            plane(i, j, k, n) = intwgt[0] * state(i, p0 - 1, j, stIdx) +
                                intwgt[1] * state(i, p0, j, stIdx) +
                                intwgt[2] * state(i, p0 + 1, j, stIdx);
          });
      } else if (m_normal == 2) {
        amrex::ParallelFor(
          bx, m_fieldNames.size(),
          [=] AMREX_GPU_DEVICE(int i, int j, int k, int n) noexcept {
            int stIdx = idx_d_p[n];
            plane(i, j, k, n) = intwgt[0] * state(i, j, p0 - 1, stIdx) +
                                intwgt[1] * state(i, j, p0, stIdx) +
                                intwgt[2] * state(i, j, p0 + 1, stIdx);
          });
      }
    }
  }

  // Count the number of level where the cut exists
  int nlevs = 0;
  for (int lev = 0; lev < a_state.size(); lev++) {
    if (!m_sliceBA[lev].empty())
      nlevs += 1;
  }

  // Build up a z-normal 2D Geom
  amrex::Vector<amrex::Geometry> pltGeoms(nlevs);
  pltGeoms[0] = m_geomLev0;
  amrex::Vector<amrex::IntVect> ref_ratio;
  amrex::IntVect rref(AMREX_D_DECL(2, 2, 1));
  for (int lev = 1; lev < nlevs; ++lev) {
    pltGeoms[lev] = amrex::refine(pltGeoms[lev - 1], rref);
    ref_ratio.push_back(rref);
  }

  // File name based on tep or time
  std::string diagfile;
  if (m_interval > 0) {
    diagfile = amrex::Concatenate(m_diagfile, a_nstep, 6);
  }
  if (m_per > 0.0) {
    diagfile = m_diagfile + std::to_string(a_time);
  }
  amrex::Vector<int> step_array(nlevs, a_nstep);
  Write2DMultiLevelPlotfile(
    diagfile, nlevs, GetVecOfConstPtrs(planeData), m_fieldNames, pltGeoms,
    a_time, step_array, ref_ratio);
}

void
DiagFramePlane::Write2DMultiLevelPlotfile(
  const std::string& a_pltfile,
  int a_nlevels,
  const amrex::Vector<const amrex::MultiFab*>& a_slice,
  const amrex::Vector<std::string>& a_varnames,
  const amrex::Vector<amrex::Geometry>& a_geoms,
  const amrex::Real& a_time,
  const amrex::Vector<int>& a_steps,
  const amrex::Vector<amrex::IntVect>& a_rref)
{
  const std::string levelPrefix = "Level_";
  const std::string mfPrefix = "Cell";
  const std::string versionName = "HyperCLaw-V1.1";

  bool callBarrier(false);
  amrex::PreBuildDirectorHierarchy(
    a_pltfile, levelPrefix, a_nlevels, callBarrier);
  amrex::ParallelDescriptor::Barrier();

  if (
    amrex::ParallelDescriptor::MyProc() ==
    amrex::ParallelDescriptor::NProcs() - 1) {
    amrex::Vector<amrex::BoxArray> boxArrays(a_nlevels);
    for (int level(0); level < boxArrays.size(); ++level) {
      boxArrays[level] = a_slice[level]->boxArray();
    }

    amrex::VisMF::IO_Buffer io_buffer(amrex::VisMF::IO_Buffer_Size);
    std::string HeaderFileName(a_pltfile + "/Header");
    std::ofstream HeaderFile;
    HeaderFile.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());
    HeaderFile.open(
      HeaderFileName.c_str(),
      std::ofstream::out | std::ofstream::trunc | std::ofstream::binary);
    if (!HeaderFile.good())
      amrex::FileOpenFailed(HeaderFileName);
    Write2DPlotfileHeader(
      HeaderFile, a_nlevels, boxArrays, a_varnames, a_geoms, a_time, a_steps,
      a_rref, versionName, levelPrefix, mfPrefix);
  }
  for (int level = 0; level < a_nlevels; ++level) {
    VisMF2D(
      *a_slice[level], amrex::MultiFabFileFullPrefix(
                         level, a_pltfile, levelPrefix, mfPrefix + "2D"));
    amrex::VisMF::Write(
      *a_slice[level],
      amrex::MultiFabFileFullPrefix(level, a_pltfile, levelPrefix, mfPrefix));
    amrex::ParallelDescriptor::Barrier();
    ReWriteLevelVisMFHeader(
      amrex::MultiFabFileFullPrefix(level, a_pltfile, levelPrefix, ""));
  }
}

void
DiagFramePlane::Write2DPlotfileHeader(
  std::ostream& HeaderFile,
  int nlevels,
  const amrex::Vector<amrex::BoxArray>& bArray,
  const amrex::Vector<std::string>& varnames,
  const amrex::Vector<amrex::Geometry>& geom,
  const amrex::Real& time,
  const amrex::Vector<int>& level_steps,
  const amrex::Vector<amrex::IntVect>& ref_ratio,
  const std::string& versionName,
  const std::string& levelPrefix,
  const std::string& mfPrefix)
{
  int finest_level(nlevels - 1);
  HeaderFile.precision(17);

  int lowerSpaceDim = AMREX_SPACEDIM - 1;

  HeaderFile << versionName << '\n';
  HeaderFile << varnames.size() << '\n';
  for (int ivar = 0; ivar < varnames.size(); ++ivar) {
    HeaderFile << varnames[ivar] << "\n";
  }
  HeaderFile << lowerSpaceDim << '\n';
  HeaderFile << time << '\n';
  HeaderFile << finest_level << '\n';
  for (int idim = 0; idim < lowerSpaceDim; ++idim) {
    HeaderFile << geom[0].ProbLo(idim) << ' ';
  }
  HeaderFile << '\n';
  for (int idim = 0; idim < lowerSpaceDim; ++idim) {
    HeaderFile << geom[0].ProbHi(idim) << ' ';
  }
  HeaderFile << '\n';
  for (int i = 0; i < finest_level; ++i) {
    HeaderFile << ref_ratio[i][0] << ' ';
  }
  HeaderFile << '\n';
  for (int i = 0; i <= finest_level; ++i) {
    printLowerDimBox(HeaderFile, geom[i].Domain(), 2);
    HeaderFile << ' ';
  }
  HeaderFile << '\n';
  for (int i = 0; i <= finest_level; ++i) {
    HeaderFile << level_steps[i] << ' ';
  }
  HeaderFile << '\n';
  for (int i = 0; i <= finest_level; ++i) {
    for (int idim = 0; idim < lowerSpaceDim; ++idim) {
      HeaderFile << geom[i].CellSize()[idim] << ' ';
    }
    HeaderFile << '\n';
  }
  HeaderFile << (int)geom[0].Coord() << '\n';
  HeaderFile << "0\n";

  for (int level = 0; level <= finest_level; ++level) {
    HeaderFile << level << ' ' << bArray[level].size() << ' ' << time << '\n';
    HeaderFile << level_steps[level] << '\n';

    const amrex::IntVect& domain_lo = geom[level].Domain().smallEnd();
    for (int i = 0; i < bArray[level].size(); ++i) {
      // Need to shift because the RealBox ctor we call takes the
      // physical location of index (0,0,0).  This does not affect
      // the usual cases where the domain index starts with 0.
      const amrex::Box& b = amrex::shift(bArray[level][i], -domain_lo);
      amrex::RealBox loc =
        amrex::RealBox(b, geom[level].CellSize(), geom[level].ProbLo());
      for (int idim = 0; idim < lowerSpaceDim; ++idim) {
        HeaderFile << loc.lo(idim) << ' ' << loc.hi(idim) << '\n';
      }
    }
    HeaderFile << amrex::MultiFabHeaderPath(level, levelPrefix, mfPrefix)
               << '\n';
  }
}

void
DiagFramePlane::ReWriteLevelVisMFHeader(const std::string& a_HeaderPath)
{
  std::string OldHeaderFileName(a_HeaderPath + "Cell_H");
  amrex::Vector<char> oldfileCharPtr;
  amrex::ParallelDescriptor::ReadAndBcastFile(
    OldHeaderFileName, oldfileCharPtr);

  if (amrex::ParallelDescriptor::IOProcessor()) {
    amrex::VisMF::IO_Buffer io_buffer_new(amrex::VisMF::IO_Buffer_Size);
    std::string HeaderFileName(a_HeaderPath + "Cell_temp_H");
    std::ofstream HeaderFile;
    HeaderFile.rdbuf()->pubsetbuf(
      io_buffer_new.dataPtr(), io_buffer_new.size());
    HeaderFile.open(
      HeaderFileName.c_str(),
      std::ofstream::out | std::ofstream::trunc | std::ofstream::binary);
    if (!HeaderFile.good())
      amrex::FileOpenFailed(HeaderFileName);

    std::string fileCharPtrString(oldfileCharPtr.dataPtr());
    std::istringstream is(fileCharPtrString, std::istringstream::in);
    std::string line, word;

    // Version, How, nComp, nGrow
    int version, how, nComp, nGrow;
    is >> version;
    is.ignore(std::streamsize(10000), '\n');
    is >> how;
    is.ignore(std::streamsize(10000), '\n');
    is >> nComp;
    is.ignore(std::streamsize(10000), '\n');
    is >> nGrow;

    HeaderFile << version << "\n";
    HeaderFile << how << "\n";
    HeaderFile << nComp << "\n";
    HeaderFile << nGrow << "\n";

    // BoxArray
    int nbox = -1, dummy = 0;
    is.ignore(std::streamsize(10000), '(') >> nbox >> dummy;
    is.ignore(std::streamsize(10000), '\n');
    HeaderFile << '(' << nbox << ' ' << 0 << '\n';
    amrex::Vector<amrex::Box> fullDimBox(nbox);
    for (int i = 0; i < nbox; ++i) {
      is >> fullDimBox[i];
      printLowerDimBox(HeaderFile, fullDimBox[i], 2);
      HeaderFile << '\n';
    }
    is.ignore(std::streamsize(10000), '\n');
    std::getline(is, line);
    HeaderFile << line << "\n";

    //------------------------------------------------------------------
    // Get list of Datafiles and correct for FAB header length in offset
    //------------------------------------------------------------------
    std::getline(is, line);
    HeaderFile << line << "\n";
    int nFabs = std::stoi(line);

    // Stash away lines and count the number of files
    amrex::Vector<std::string> dataFiles;
    amrex::Vector<std::string> fablines(nFabs);
    for (int i = 0; i < nFabs; ++i) {
      std::getline(is, fablines[i]);
      std::istringstream lis(fablines[i]);
      std::string w1, w2, w3;
      lis >> w1;
      lis >> w2;
      lis >> w3;
      int offset = std::stoi(w3);
      if (offset == 0)
        dataFiles.push_back(w2);
    }

    // Map files and FABs
    int nDataFiles = dataFiles.size();
    amrex::Vector<int> nFabsInFiles(nDataFiles);
    amrex::Vector<amrex::Vector<int>> FabsInFiles(nDataFiles);
    for (int i = 0; i < nFabs; ++i) {
      std::istringstream lis(fablines[i]);
      std::string w1, w2, w3;
      lis >> w1;
      lis >> w2;
      lis >> w3;
      for (int f = 0; f < nDataFiles; ++f) {
        if (w2 == dataFiles[f]) {
          nFabsInFiles[f] += 1;
          FabsInFiles[f].push_back(i);
          break;
        }
      }
    }
    // Check
    for (int f = 0; f < nDataFiles; ++f) {
      if (nFabsInFiles[f] != FabsInFiles[f].size()) {
        amrex::Abort("Something went wrong while counting FABs in files");
      }
    }

    // Write out the FAB on disk list, corrected
    for (int i = 0; i < nFabs; ++i) {
      std::istringstream lis(fablines[i]);
      std::string w1, w2, w3;
      lis >> w1;
      lis >> w2;
      lis >> w3;
      int offset = std::stoi(w3);
      for (int f = 0; f < nDataFiles; ++f) {
        if (w2 == dataFiles[f]) {
          for (int fa = 0; fa < FabsInFiles[f].size(); ++fa) {
            if (FabsInFiles[f][fa] == i) {
              if (fa > 0) {
                for (int fback = fa - 1; fback >= 0; --fback) {
                  offset -= diff_2D3D_header(fullDimBox[fback], nComp);
                }
              }
            }
          }
        }
      }
      HeaderFile << w1 << " " << w2 << " " << offset << "\n";
    }

    // Just pass from istream to ostream `till the end
    while (is) {
      std::getline(is, line);
      HeaderFile << line << "\n";
    }
    HeaderFile.close();

    // Replace header file
    std::rename(HeaderFileName.c_str(), OldHeaderFileName.c_str());
  }

  // Replace 3D data file by 2D ones
  std::string oldname = amrex::Concatenate(
    a_HeaderPath + "Cell2D_D_", amrex::ParallelDescriptor::MyProc(), 5);
  if (amrex::FileSystem::Exists(oldname.c_str())) {
    std::string newname =
      std::regex_replace(oldname, std::regex("Cell2D_"), "Cell_");
    std::remove(newname.c_str());
    std::rename(oldname.c_str(), newname.c_str());
  }
}

void
DiagFramePlane::VisMF2D(
  const amrex::MultiFab& a_mf, const std::string& a_mf_name)
{
  auto whichRD = amrex::FArrayBox::getDataDescriptor();
  bool doConvert(*whichRD != amrex::FPC::NativeRealDescriptor());

  amrex::Long bytesWritten(0);

  std::string filePrefix(a_mf_name + "_D_");

  bool calcMinMax = false;
  amrex::VisMF::Header::Version currentVersion =
    amrex::VisMF::Header::Version_v1;
  amrex::VisMF::How how = amrex::VisMF::How::NFiles;
  amrex::VisMF::Header hdr(a_mf, how, currentVersion, calcMinMax);

  int nOutFiles =
    std::max(1, std::min(amrex::ParallelDescriptor::NProcs(), 256));
  bool groupSets = false;
  bool setBuf = true;

  amrex::NFilesIter nfi(nOutFiles, filePrefix, groupSets, setBuf);

  // Check if mf has sparse data
  bool useSparseFPP = false;
  const amrex::Vector<int>& pmap = a_mf.DistributionMap().ProcessorMap();
  std::set<int> procsWithData;
  amrex::Vector<int> procsWithDataVector;
  for (int i(0); i < pmap.size(); ++i) {
    procsWithData.insert(pmap[i]);
  }
  if (static_cast<int>(procsWithData.size()) < nOutFiles) {
    useSparseFPP = true;
    for (std::set<int>::iterator it = procsWithData.begin();
         it != procsWithData.end(); ++it) {
      procsWithDataVector.push_back(*it);
    }
  }

  if (useSparseFPP) {
    nfi.SetSparseFPP(procsWithDataVector);
  } else {
    nfi.SetDynamic();
  }
  for (; nfi.ReadyToWrite(); ++nfi) {
    int whichRDBytes(whichRD->numBytes()), nFABs(0);
    amrex::Long writeDataItems(0), writeDataSize(0);
    for (amrex::MFIter mfi(a_mf); mfi.isValid(); ++mfi) {
      const amrex::FArrayBox& fab = a_mf[mfi];
      std::stringstream hss;
      write_2D_header(hss, fab, fab.nComp());
      bytesWritten += static_cast<std::streamoff>(hss.tellp());
      bytesWritten += fab.box().numPts() * a_mf.nComp() * whichRDBytes;
      ++nFABs;
    }
    char* allFabData = nullptr;
    bool canCombineFABs = false;
    if ((nFABs > 1 || doConvert) && amrex::VisMF::GetUseSingleWrite()) {
      allFabData = new (std::nothrow) char[bytesWritten];
    } // ---- else { no need to make a copy for one fab }
    if (allFabData == nullptr) {
      canCombineFABs = false;
    } else {
      canCombineFABs = true;
    }

    if (canCombineFABs) {
      amrex::Long writePosition = 0;
      for (amrex::MFIter mfi(a_mf); mfi.isValid(); ++mfi) {
        int hLength(0);
        const amrex::FArrayBox& fab = a_mf[mfi];
        writeDataItems = fab.box().numPts() * a_mf.nComp();
        writeDataSize = writeDataItems * whichRDBytes;
        char* afPtr = allFabData + writePosition;
        std::stringstream hss;
        write_2D_header(hss, fab, fab.nComp());
        hLength = static_cast<std::streamoff>(hss.tellp());
        auto tstr = hss.str();
        std::memcpy(afPtr, tstr.c_str(), hLength); // ---- the fab header
        amrex::Real const* fabdata = fab.dataPtr();
#ifdef AMREX_USE_GPU
        std::unique_ptr<amrex::FArrayBox> hostfab;
        if (fab.arena()->isManaged() || fab.arena()->isDevice()) {
          hostfab = std::make_unique<amrex::FArrayBox>(
            fab.box(), fab.nComp(), amrex::The_Pinned_Arena());
          amrex::Gpu::dtoh_memcpy_async(
            hostfab->dataPtr(), fab.dataPtr(),
            fab.size() * sizeof(amrex::Real));
          amrex::Gpu::streamSynchronize();
          fabdata = hostfab->dataPtr();
        }
#endif
        memcpy(afPtr + hLength, fabdata, writeDataSize);
        writePosition += hLength + writeDataSize;
      }
      nfi.Stream().write(allFabData, bytesWritten);
      nfi.Stream().flush();
      delete[] allFabData;
    } else {
      for (amrex::MFIter mfi(a_mf); mfi.isValid(); ++mfi) {
        int hLength = 0;
        const amrex::FArrayBox& fab = a_mf[mfi];
        writeDataItems = fab.box().numPts() * a_mf.nComp();
        writeDataSize = writeDataItems * whichRDBytes;
        std::stringstream hss;
        write_2D_header(hss, fab, fab.nComp());
        hLength = static_cast<std::streamoff>(hss.tellp());
        auto tstr = hss.str();
        nfi.Stream().write(tstr.c_str(), hLength); // ---- the fab header
        nfi.Stream().flush();
        amrex::Real const* fabdata = fab.dataPtr();
#ifdef AMREX_USE_GPU
        std::unique_ptr<amrex::FArrayBox> hostfab;
        if (fab.arena()->isManaged() || fab.arena()->isDevice()) {
          hostfab = std::make_unique<amrex::FArrayBox>(
            fab.box(), fab.nComp(), amrex::The_Pinned_Arena());
          amrex::Gpu::dtoh_memcpy_async(
            hostfab->dataPtr(), fab.dataPtr(),
            fab.size() * sizeof(amrex::Real));
          amrex::Gpu::streamSynchronize();
          fabdata = hostfab->dataPtr();
        }
#endif
        nfi.Stream().write((char*)fabdata, writeDataSize);
        nfi.Stream().flush();
      }
    }
  }
}

void
DiagFramePlane::write_2D_header(
  std::ostream& os, const amrex::FArrayBox& f, int nvar)
{
  os << "FAB " << amrex::FPC::NativeRealDescriptor();
  amrex::StreamRetry sr(os, "FABio_write_header", 4);
  while (sr.TryOutput()) {
    printLowerDimBox(os, f.box(), 2);
    os << ' ' << nvar << '\n';
  }
}

int
DiagFramePlane::diff_2D3D_header(const amrex::Box a_box, int nvar)
{
  std::stringstream hss2D;
  std::stringstream hss3D;
  hss2D << "FAB " << amrex::FPC::NativeRealDescriptor();
  amrex::StreamRetry sr(hss2D, "FABio_write_header", 4);
  while (sr.TryOutput()) {
    printLowerDimBox(hss2D, a_box, 2); // Skipping z-dir
    hss2D << ' ' << nvar << '\n';
  }

  hss3D << "FAB " << amrex::FPC::NativeRealDescriptor();
  amrex::StreamRetry sr3(hss3D, "FABio_write_header", 4);
  while (sr3.TryOutput()) {
    hss3D << a_box;
    hss3D << ' ' << nvar << '\n';
  }

  int Length2D = static_cast<std::streamoff>(hss2D.tellp());
  int Length3D = static_cast<std::streamoff>(hss3D.tellp());

  return Length3D - Length2D;
}
