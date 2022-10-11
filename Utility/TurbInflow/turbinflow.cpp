#include <turbinflow.H>

namespace pele::physics::turbinflow {
void
TurbInflow::init(amrex::Geometry const& /*geom*/)
{
  amrex::ParmParse ppr;

  int n_tp = 0;
  n_tp = ppr.countval("turbinflows");
  amrex::Vector<std::string> tp_list;
  if (n_tp > 0) {
    tp.resize(n_tp);
    tp_list.resize(n_tp);
    for (int n = 0; n < n_tp; n++) {
      ppr.get("turbinflows", tp_list[n], n);
    }
  }

  for (int n = 0; n < n_tp; n++) {

    amrex::ParmParse pp("turbinflow." + tp_list[n]);
    if (pp.countval("turb_file") > 0) {

      // Query data
      pp.query("turb_file", tp[n].m_turb_file);
      tp[n].dir = -1;
      pp.query("dir", tp[n].dir);
      AMREX_ASSERT_WITH_MESSAGE(
        tp[n].dir >= 0 && tp[n].dir < AMREX_SPACEDIM,
        "Injection direction is needed: 0, 1 or 2");
      std::string side;
      pp.query("side", side);
      if (side == "low") {
        tp[n].side = amrex::Orientation::low;
      } else if (side == "high") {
        tp[n].side = amrex::Orientation::high;
      } else {
        amrex::Abort("turbinflow.side can only be low or high");
      }
      pp.query("time_offset", tp[n].time_shift);
      pp.query("turb_scale_loc", tp[n].turb_scale_loc);
      pp.query("turb_scale_vel", tp[n].turb_scale_vel);
      amrex::Print() << "Initializing turbInflow " << tp_list[n]
                     << " with file " << tp[n].m_turb_file
                     << " (location coordinates in will be scaled by "
                     << tp[n].turb_scale_loc
                     << " and velocity out to be scaled by "
                     << tp[n].turb_scale_vel << ") \n";

      // Get the turbcenter on the injection face
      amrex::Vector<amrex::Real> turb_center(AMREX_SPACEDIM - 1, 0);
      pp.getarr("turb_center", turb_center);
      AMREX_ASSERT_WITH_MESSAGE(
        turb_center.size() == AMREX_SPACEDIM - 1,
        "turb_center must have AMREX_SPACEDIM-1 elements");
      for (int idim = 0; idim < turb_center.size(); ++idim) {
        turb_center[idim] *= tp[n].turb_scale_loc;
      }

      pp.query("turb_nplane", tp[n].nplane);
      AMREX_ASSERT(tp[n].nplane > 0);
      pp.query("turb_conv_vel", tp[n].turb_conv_vel);
      AMREX_ASSERT(tp[n].turb_conv_vel > 0);

      // Set other stuff
      std::string turb_header = tp[n].m_turb_file + "/HDR";
      std::ifstream is(turb_header.c_str());
      if (!is.is_open()) {
        amrex::Abort("Unable to open input file " + turb_header);
      }
      amrex::Array<int, AMREX_SPACEDIM> npts = {{0}};
      amrex::Array<amrex::Real, AMREX_SPACEDIM> probsize = {{0}};
      amrex::Array<int, AMREX_SPACEDIM> iper = {{0}};
      is >> npts[0] >> npts[1] >> npts[2];
      is >> probsize[0] >> probsize[1] >> probsize[2];
      is >> iper[0] >> iper[1] >>
        iper[2]; // Unused - we assume it is always fully periodic

      for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        tp[n].dx[idim] = probsize[idim] / amrex::Real(npts[idim] - 1);
        tp[n].dxinv[idim] = 1.0 / tp[n].dx[idim];
      }

      // The following is relative to the injection face:
      // 0 and 1 are transverse directions, 2 is normal
      // one ghost point on each side, tangential to inflow face
      tp[n].pboxsize[0] = probsize[0] - 2.0 * tp[n].dx[0];
      tp[n].pboxsize[1] = probsize[1] - 2.0 * tp[n].dx[1];
      tp[n].pboxsize[2] = probsize[2];

      tp[n].npboxcells[0] = npts[0] - 3;
      tp[n].npboxcells[1] = npts[1] - 3;
      tp[n].npboxcells[2] = npts[2];

      // Center the turbulence
      tp[n].pboxlo[0] = turb_center[0] - 0.5 * tp[n].pboxsize[0];
      tp[n].pboxlo[1] = turb_center[1] - 0.5 * tp[n].pboxsize[1];
      tp[n].pboxlo[2] = 0.;

      amrex::Box sbx(
        amrex::IntVect(AMREX_D_DECL(1, 1, 1)),
        amrex::IntVect(AMREX_D_DECL(npts[0], npts[1], tp[n].nplane)));

      tp[n].sdata = new amrex::FArrayBox(sbx, 3, amrex::The_Async_Arena());

      tp[n].kmax = npts[2];

      if (tp[n].isswirltype) {
        for (int i = 0; i < tp[n].kmax; i++) {
          amrex::Real rdummy = 0.0;
          is >> rdummy; // Time for each plane - unused at the moment
        }
      }

      // Offset for each plane in Binary TurbFile
      tp[n].m_offset.resize(tp[n].kmax * AMREX_SPACEDIM);
      tp[n].offset = tp[n].m_offset.data();
      tp[n].offset_size = tp[n].m_offset.size();
      for (int i = 0; i < tp[n].offset_size; i++) {
        is >> tp[n].offset[i];
      }
      is.close();
    }
    turbinflow_initialized = true;
  }
}

void
TurbInflow::add_turb(
  amrex::Box const& bx,
  amrex::FArrayBox& data,
  const int dcomp,
  amrex::Geometry const& geom,
  const amrex::Real time,
  const int dir,
  const amrex::Orientation::Side& side)
{
  AMREX_ASSERT(turbinflow_initialized);

  // Box on which we will access data
  amrex::Box bvalsBox = bx;
  int planeLoc =
    (side == amrex::Orientation::low ? geom.Domain().smallEnd()[dir] - 1
                                     : geom.Domain().bigEnd()[dir] + 1);
  bvalsBox.setSmall(dir, planeLoc);
  bvalsBox.setBig(dir, planeLoc);

  // Define box that we will fill with turb: need to be z-normal
  // Get transverse directions
  int tdir1 = (dir != 0) ? 0 : 1;
  int tdir2 = (dir != 0) ? ((dir == 2) ? 1 : 2) : 2;
  int tr1Lo = bvalsBox.smallEnd()[tdir1];
  int tr1Hi = bvalsBox.bigEnd()[tdir1];
  int tr2Lo = bvalsBox.smallEnd()[tdir2];
  int tr2Hi = bvalsBox.bigEnd()[tdir2];
  const amrex::IntVect lo(AMREX_D_DECL(tr1Lo, tr2Lo, planeLoc));
  const amrex::IntVect hi(AMREX_D_DECL(tr1Hi, tr2Hi, planeLoc));
  amrex::Box turbBox(lo, hi);
  amrex::FArrayBox v(turbBox, 3, amrex::The_Async_Arena());
  v.setVal<amrex::RunOn::Device>(0);

  // Add turbulence from all the tp acting on this face
  for (int n = 0; n < tp.size(); n++) {

    if (tp[n].dir == dir && tp[n].side == side) {

      // 0 and 1 are the two transverse directions
      amrex::Vector<amrex::Real> x(turbBox.size()[0]), y(turbBox.size()[1]);
      for (int i = turbBox.smallEnd()[0]; i <= turbBox.bigEnd()[0]; ++i) {
        x[i - turbBox.smallEnd()[0]] =
          (geom.ProbLo()[tdir1] + (i + 0.5) * geom.CellSize(tdir1)) *
          tp[n].turb_scale_loc;
      }
      for (int j = turbBox.smallEnd()[1]; j <= turbBox.bigEnd()[1]; ++j) {
        y[j - turbBox.smallEnd()[1]] =
          (geom.ProbLo()[tdir2] + (j + 0.5) * geom.CellSize(tdir2)) *
          tp[n].turb_scale_loc;
      }

      // Get the turbulence
      amrex::Real z =
        (time + tp[n].time_shift) * tp[n].turb_conv_vel * tp[n].turb_scale_loc;
      fill_turb_plane(tp[n], x, y, z, v);
    }
  }

  // Moving it into data
  set_turb(dir, tdir1, tdir2, v, data, dcomp);

#if 0
  std::string junk = "TurbV_AftTP"+std::to_string(n)+"_D";
  std::ofstream os;
  os.precision(15);
  os.open(junk.c_str());
  data.writeOn(os);
  os.close();
  amrex::Abort();
#endif
}

void
TurbInflow::set_turb(
  int normDir,
  int transDir1,
  int transDir2,
  amrex::FArrayBox& v,
  amrex::FArrayBox& data,
  const int dcomp)
{
  // copy velocity fluctuations from plane into data
  const auto& box = v.box(); // z-normal plane
  const auto& v_in = v.array();
  const auto& v_out = data.array(dcomp);

  amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    // From z-normal box index to data box index
    int idx[3] = {0};
    idx[transDir1] = i;
    idx[transDir2] = j;
    idx[normDir] = k;
    v_out(idx[0], idx[1], idx[2], transDir1) =
      v_in(i, j, k, 0); // transverse velocity 1
    v_out(idx[0], idx[1], idx[2], transDir2) =
      v_in(i, j, k, 1); // transverse velocity 2
    v_out(idx[0], idx[1], idx[2], normDir) =
      v_in(i, j, k, 2); // normal velocity
  });
}

void
TurbInflow::read_one_turb_plane(TurbParm& a_tp, int iplane, int k)
{
  // There are AMREX_SPACEDIM * kmax planes of FABs.
  // The first component are in the first kmax planes,
  // the second component in the next kmax planes, ....
  // Note also that both (*plane) and (*ncomp) start from
  // 1 not 0 since they're passed from Fortran.

  std::string turb_data = a_tp.m_turb_file + "/DAT";
  std::ifstream ifs(turb_data.c_str());
  if (!ifs.is_open()) {
    amrex::Abort("Unable to open input file " + turb_data);
  }

  amrex::Box dstBox = a_tp.sdata->box();
  dstBox.setSmall(2, iplane);
  dstBox.setBig(2, iplane);

  for (int n = 0; n < AMREX_SPACEDIM; ++n) {

    const long offset_idx = (k + 1) + (n * a_tp.kmax);
    AMREX_ASSERT_WITH_MESSAGE(
      offset_idx < a_tp.offset_size, "Bad turb fab offset idx");

    const long start = a_tp.offset[offset_idx];

    ifs.seekg(start, std::ios::beg);

    if (!ifs.good()) {
      amrex::Abort("getplane(): seekg() failed");
    }

    amrex::FArrayBox tmp;
    tmp.readFrom(ifs);
    amrex::Box srcBox = tmp.box();

    a_tp.sdata->copy<amrex::RunOn::Device>(tmp, srcBox, 0, dstBox, n, 1);
  }
  ifs.close();
}

void
TurbInflow::read_turb_planes(TurbParm& a_tp, amrex::Real z)
{
  int izlo = (int)(round(z * a_tp.dxinv[2])) - 1;
  int izhi = izlo + a_tp.nplane - 1;
  a_tp.szlo = static_cast<amrex::Real>(izlo) * a_tp.dx[2];
  a_tp.szhi = static_cast<amrex::Real>(izhi) * a_tp.dx[2];

#if 0
  amrex::AllPrint() << "read_turb_planes filling " << izlo << " to " << izhi
                 << " covering " << a_tp.szlo + 0.5 * a_tp.dx[2]
                 << " to "       << a_tp.szhi - 0.5 * a_tp.dx[2] << " for z = " << z << std::endl;
#endif

  for (int iplane = 1; iplane <= a_tp.nplane; ++iplane) {
    int k = (izlo + iplane - 1) % (a_tp.npboxcells[2] - 2);
    read_one_turb_plane(a_tp, iplane, k);
  }
#if 0
  int myproc = amrex::ParallelDescriptor::MyProc();
  std::string junk = "TurbData_proc"+std::to_string(myproc)+"_D";
  std::ofstream os;
  os.precision(15);
  os.open(junk.c_str());
  a_tp.sdata->writeOn(os);
  os.close();
  //amrex::Abort();
#endif
}

void
TurbInflow::fill_turb_plane(
  TurbParm& a_tp,
  const amrex::Vector<amrex::Real>& x,
  const amrex::Vector<amrex::Real>& y,
  amrex::Real z,
  amrex::FArrayBox& v)
{
  if (
    (z < a_tp.szlo + 0.5 * a_tp.dx[2]) || (z > a_tp.szhi - 0.5 * a_tp.dx[2])) {
#if 0
    {
      amrex::AllPrint() << "Reading new data because z " << z << " is outside " << a_tp.szlo + 0.5 * a_tp.dx[2] << " and "
                     << a_tp.szhi - 0.5 * a_tp.dx[2] << std::endl;
    }
#endif
    read_turb_planes(a_tp, z);
  }

  const auto& bx = v.box();
  const auto& vd = v.array();

  amrex::Gpu::DeviceVector<amrex::Real> x_dev(x.size());
  amrex::Gpu::DeviceVector<amrex::Real> y_dev(y.size());
  amrex::Gpu::copyAsync(
    amrex::Gpu::hostToDevice, x.begin(), x.end(), x_dev.begin());
  amrex::Gpu::copyAsync(
    amrex::Gpu::hostToDevice, y.begin(), y.end(), y_dev.begin());
  amrex::Real* xd = x_dev.data();
  amrex::Real* yd = y_dev.data();

  amrex::Real velScale = (a_tp.side == amrex::Orientation::high)
                           ? -a_tp.turb_scale_vel
                           : a_tp.turb_scale_vel;
  const auto& npboxcells = a_tp.npboxcells;
  const auto& pboxlo = a_tp.pboxlo;
  const auto& szlo = a_tp.szlo;
  const auto& dxinv = a_tp.dxinv;
  const auto& sd = a_tp.sdata->array();

  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
    amrex::Real cx[3], cy[3], cz[3], ydata[3];
    amrex::Real zdata[3][3];

    amrex::Real zz =
      (z - szlo) * dxinv[2];        // How many dz away from the left side ?
    int k0 = (int)(std::round(zz)); // What's the closest point ?
    zz -= amrex::Real(k0);
    cz[0] = 0.5 * (zz - 1.0) * (zz - 2.0); // Weight of k0 - 1
    cz[1] = zz * (2.0 - zz);               // Weight of k0
    cz[2] = 0.5 * zz * (zz - 1.0);         // Weight of k0 + 1
    k0 += 1;                               // Index starting at 1

    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
      amrex::Real xx = (xd[i - bx.smallEnd(0)] - pboxlo[0]) * dxinv[0];
      amrex::Real yy = (yd[j - bx.smallEnd(1)] - pboxlo[1]) * dxinv[1];
      int i0 = (int)(std::round(xx));
      int j0 = (int)(std::round(yy));
      xx -= amrex::Real(i0);
      yy -= amrex::Real(j0);
      cx[0] = 0.5 * (xx - 1.0) * (xx - 2.0);
      cy[0] = 0.5 * (yy - 1.0) * (yy - 2.0);
      cx[1] = xx * (2.0 - xx);
      cy[1] = yy * (2.0 - yy);
      cx[2] = 0.5 * xx * (xx - 1.0);
      cy[2] = 0.5 * yy * (yy - 1.0);

      if (i0 >= 0 && i0 < npboxcells[0] && j0 >= 0 && j0 < npboxcells[1]) {
        i0 += 2;
        j0 += 2;
        for (int ii = 0; ii <= 2; ++ii) {
          for (int jj = 0; jj <= 2; ++jj) {
            zdata[ii][jj] = cz[0] * sd(i0 + ii, j0 + jj, k0 - 1, n) +
                            cz[1] * sd(i0 + ii, j0 + jj, k0, n) +
                            cz[2] * sd(i0 + ii, j0 + jj, k0 + 1, n);
          }
        }
        for (int ii = 0; ii <= 2; ++ii) {
          ydata[ii] =
            cy[0] * zdata[ii][0] + cy[1] * zdata[ii][1] + cy[2] * zdata[ii][2];
        }
        vd(i, j, k, n) = cx[0] * ydata[0] + cx[1] * ydata[1] + cx[2] * ydata[2];
        vd(i, j, k, n) *= velScale;
      }
    }
  });
  amrex::Gpu::synchronize(); // Ensure that DeviceVector's don't leave scope
                             // early
}
} // namespace pele::physics::turbinflow
