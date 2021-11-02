#include <turbinflow.H>

void
init_turbinflow(
  const std::string& turb_file,
  amrex::Real turb_scale_loc,
  amrex::Real turb_scale_vel,
  const amrex::Vector<amrex::Real>& turb_center,
  amrex::Real turb_conv_vel,
  int nplane,
  TurbParm& tp)
{
  amrex::Print() << "Initializing turbulence file: " << turb_file
                 << " (location coordinates in will be scaled by "
                 << turb_scale_loc << " and velocity out to be scaled by "
                 << turb_scale_vel << ")" << std::endl;

  tp.tph->turb_file = turb_file;
  tp.nplane = nplane;
  tp.turb_conv_vel = turb_conv_vel;
  tp.turb_scale_loc = turb_scale_loc;
  tp.turb_scale_vel = turb_scale_vel;

  std::string turb_header = tp.tph->turb_file + "/HDR";
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

  for (int n = 0; n < AMREX_SPACEDIM; ++n) {
    tp.dx[n] = probsize[n] / amrex::Real(npts[n] - 1);
    tp.dxinv[n] = 1.0 / tp.dx[n];
  }

  // one ghost point on each side, tangential to inflow face
  tp.pboxsize[0] = probsize[0] - 2.0 * tp.dx[0];
  tp.pboxsize[1] = probsize[1] - 2.0 * tp.dx[1];
  tp.pboxsize[2] = probsize[2];

  tp.npboxcells[0] = npts[0] - 3;
  tp.npboxcells[1] = npts[1] - 3;
  tp.npboxcells[2] = npts[2];

  // Center the turbulence)
  tp.pboxlo[0] = turb_center[0] - 0.5 * tp.pboxsize[0];
  tp.pboxlo[1] = turb_center[1] - 0.5 * tp.pboxsize[1];
  tp.pboxlo[2] = 0.;

  amrex::Box sbx(
    amrex::IntVect(AMREX_D_DECL(1, 1, 1)),
    amrex::IntVect(AMREX_D_DECL(npts[0], npts[1], tp.nplane)));

  tp.sdata = new amrex::FArrayBox(sbx, 3);

  tp.kmax = npts[2];

  amrex::Real rdummy;
  if (tp.isswirltype) {
    for (int i = 0; i < tp.kmax; i++) {
      is >> rdummy; // Time for each plane - unused at the moment
    }
  }
  tp.tph->offset_dv.resize(tp.kmax * AMREX_SPACEDIM);
  tp.offset = tp.tph->offset_dv.data();
  tp.offset_size = tp.tph->offset_dv.size();
  for (int i = 0; i < tp.offset_size; i++) {
    is >> tp.offset[i];
  }
  is.close();

  tp.turbinflow_initialized = true;
}

void
read_one_turb_plane(int iplane, int k, TurbParm& tp)
{
  //
  // There are AMREX_SPACEDIM * kmax planes of FABs.
  // The first component are in the first kmax planes,
  // the second component in the next kmax planes, ....
  // Note also that both (*plane) and (*ncomp) start from
  // 1 not 0 since they're passed from Fortran.
  //
  std::string turb_data = tp.tph->turb_file + "/DAT";
  std::ifstream ifs(turb_data.c_str());
  if (!ifs.is_open()) {
    amrex::Abort("Unable to open input file " + turb_data);
  }

  amrex::Box dstBox = tp.sdata->box();
  dstBox.setSmall(2, iplane);
  dstBox.setBig(2, iplane);

  for (int n = 0; n < AMREX_SPACEDIM; ++n) {

    const long offset_idx = (k + 1) + (n * tp.kmax);
    AMREX_ASSERT_WITH_MESSAGE(
      offset_idx < tp.offset_size, "Bad turb fab offset idx");

    const long start = tp.offset[offset_idx];

    ifs.seekg(start, std::ios::beg);

    if (!ifs.good())
      amrex::Abort("getplane(): seekg() failed");

    amrex::FArrayBox tmp;
    tmp.readFrom(ifs);
    amrex::Box srcBox = tmp.box();

    tp.sdata->copy<amrex::RunOn::Device>(tmp, srcBox, 0, dstBox, n, 1);
  }
  ifs.close();
}

void
read_turb_planes(amrex::Real z, TurbParm& tp)
{
  int izlo = (int)(round(z * tp.dxinv[2])) - 1;
  int izhi = izlo + tp.nplane - 1;
  tp.szlo = izlo * tp.dx[2];
  tp.szhi = izhi * tp.dx[2];

#if 0
  amrex::Print() << "read_turb_planes filling " << izlo << " to " << izhi
                 << " covering " << tp.szlo + 0.5 * tp.dx[2]
                 << " to "       << tp.szhi - 0.5 * tp.dx[2] << " for z = " << z << std::endl;
#endif

  for (int iplane = 1; iplane <= tp.nplane; ++iplane) {
    int k = (izlo + iplane - 1) % (tp.npboxcells[2] - 2);
    read_one_turb_plane(iplane, k, tp);
  }
  tp.turbinflow_planes_initialized = true;
}

void
fill_turb_plane(
  const amrex::Vector<amrex::Real>& x,
  const amrex::Vector<amrex::Real>& y,
  amrex::Real z,
  amrex::FArrayBox& v,
  TurbParm& tp)
{
  if (
    (!tp.turbinflow_planes_initialized) || (z < tp.szlo + 0.5 * tp.dx[2]) ||
    (z > tp.szhi - 0.5 * tp.dx[2])) {
#if 0
    if (!tp.turbinflow_planes_initialized) {
      amrex::Print() << "Reading new data because planes uninitialized at z: " << z << std::endl;
    }
    else {
      amrex::Print() << "Reading new data because z " << z << " is outside " << tp.szlo + 0.5 * tp.dx[2] << " and "
                     << tp.szhi - 0.5 * tp.dx[2] << std::endl;
    }
#endif
    read_turb_planes(z, tp);
  }

  const auto& bx = v.box();
  const auto& vd = v.array();
  const auto& sd = tp.sdata->array();

  amrex::Gpu::DeviceVector<amrex::Real> x_dev(x.size());
  amrex::Gpu::DeviceVector<amrex::Real> y_dev(y.size());
  amrex::Gpu::copyAsync(
    amrex::Gpu::hostToDevice, x.begin(), x.end(), x_dev.begin());
  amrex::Gpu::copyAsync(
    amrex::Gpu::hostToDevice, y.begin(), y.end(), y_dev.begin());
  amrex::Real* xd = x_dev.data();
  amrex::Real* yd = y_dev.data();
  amrex::ParallelFor(
    bx,
    [z, sd, vd, tp, bx, xd, yd] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      amrex::Real cx[3], cy[3], cz[3], ydata[3];
      amrex::Real zdata[3][3];

      amrex::Real zz = (z - tp.szlo) * tp.dxinv[2];
      int k0 = (int)(std::round(zz)) - 1;
      zz -= amrex::Real(k0);
      cz[0] = 0.5 * (zz - 1.0) * (zz - 2.0);
      cz[1] = zz * (2.0 - zz);
      cz[2] = 0.5 * zz * (zz - 1.0);
      k0 += 2;
      k0 = amrex::min(amrex::max(k0, 1), tp.nplane - 2);

      for (int n = 0; n < 3; ++n) {
        amrex::Real xx = (xd[i - bx.smallEnd(0)] - tp.pboxlo[0]) * tp.dxinv[0];
        amrex::Real yy = (yd[j - bx.smallEnd(1)] - tp.pboxlo[1]) * tp.dxinv[1];
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

        if (
          i0 >= 0 && i0 < tp.npboxcells[0] && j0 >= 0 &&
          j0 < tp.npboxcells[1]) {
          i0 += 2;
          j0 += 2;
          for (int ii = 0; ii <= 2; ++ii) {
            for (int jj = 0; jj <= 2; ++jj) {
              zdata[ii][jj] = cz[0] * sd(i0 + ii, j0 + jj, k0, n) +
                              cz[1] * sd(i0 + ii, j0 + jj, k0 + 1, n) +
                              cz[2] * sd(i0 + ii, j0 + jj, k0 + 2, n);
            }
          }
          for (int ii = 0; ii <= 2; ++ii) {
            ydata[ii] = cy[0] * zdata[ii][0] + cy[1] * zdata[ii][1] +
                        cy[2] * zdata[ii][2];
          }
          vd(i, j, k, n) =
            cx[0] * ydata[0] + cx[1] * ydata[1] + cx[2] * ydata[2];
        } else {
          vd(i, j, k, n) = 0.0;
        }
      }
    });
  amrex::Gpu::synchronize(); // Ensure that DeviceVector's don't leave scope
                             // early
}

void
add_turb(amrex::Box const& bx,
         amrex::FArrayBox& data,
         const int dcomp,
         amrex::Geometry const& geom,
         const amrex::Real time,
         const int dir,
         const amrex::Orientation::Side& side,
         TurbParm& tp)
{
  AMREX_ASSERT_WITH_MESSAGE(
    dir == 2,
    "Sadly, the fluctuation code currently only works in the third dimension");
  AMREX_ASSERT(tp.turbinflow_initialized);

  amrex::Box bvalsBox = bx;
  int planeLoc =
    (side == amrex::Orientation::low ? geom.Domain().smallEnd()[dir] - 1
                                     : geom.Domain().bigEnd()[dir] + 1);
  bvalsBox.setSmall(dir, planeLoc);
  bvalsBox.setBig(dir, planeLoc);

  amrex::FArrayBox v(bvalsBox, 3);

  amrex::Vector<amrex::Real> x(bvalsBox.size()[0]), y(bvalsBox.size()[1]);
  for (int i = bvalsBox.smallEnd()[0]; i <= bvalsBox.bigEnd()[0]; ++i) {
    x[i - bvalsBox.smallEnd()[0]] =
      (geom.ProbLo()[0] + (i + 0.5) * geom.CellSize(0)) * tp.turb_scale_loc;
  }
  for (int j = bvalsBox.smallEnd()[1]; j <= bvalsBox.bigEnd()[1]; ++j) {
    y[j - bvalsBox.smallEnd()[1]] =
      (geom.ProbLo()[1] + (j + 0.5) * geom.CellSize(1)) * tp.turb_scale_loc;
  }

  v.setVal<amrex::RunOn::Device>(0);
  amrex::Real z = time * tp.turb_conv_vel * tp.turb_scale_loc;
  fill_turb_plane(x, y, z, v, tp);
  if (side == amrex::Orientation::high) {
    v.mult<amrex::RunOn::Device>(-tp.turb_scale_vel);
  } else {
    v.mult<amrex::RunOn::Device>(tp.turb_scale_vel);
  }
  amrex::Box ovlp = bvalsBox & data.box();
  set_turb(v,data,dcomp);
}

void set_turb(amrex::FArrayBox& v,
              amrex::FArrayBox& data,
              const int dcomp)
{
  //copy velocity fluctuations from plane into data
  const auto& box   = v.box();
  const auto& v_in  = v.array();
  const auto& v_out = data.array(dcomp);

  amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
                     v_out(i,j,k,0) =  v_in(i,j,k,0); //UMX
                     v_out(i,j,k,1) =  v_in(i,j,k,1); //UMY
#if AMREX_SPACEDIM == 3
                     v_out(i,j,k,2) =  v_in(i,j,k,2); //UMZ
#endif
    });
}

void
fill_with_turb(
  amrex::Box const& bx,
  amrex::FArrayBox& data,
  const int dcomp,
  amrex::Geometry const& geom,
  TurbParm& tp)
{
  int dir = 2;
  AMREX_ASSERT(tp.turbinflow_initialized);

  for (int planeloc = bx.smallEnd()[2]; planeloc <= bx.bigEnd()[2];
       ++planeloc) {
    amrex::Box bvalsBox = bx;
    bvalsBox.setSmall(dir, planeloc);
    bvalsBox.setBig(dir, planeloc);

    amrex::FArrayBox v(bvalsBox, 3);
    v.setVal<amrex::RunOn::Device>(0);

    amrex::Vector<amrex::Real> x(bvalsBox.size()[0]), y(bvalsBox.size()[1]);
    for (int i = bvalsBox.smallEnd()[0]; i <= bvalsBox.bigEnd()[0]; ++i) {
      x[i - bvalsBox.smallEnd()[0]] =
        (geom.ProbLo()[0] + (i + 0.5) * geom.CellSize(0)) * tp.turb_scale_loc;
    }
    for (int j = bvalsBox.smallEnd()[1]; j <= bvalsBox.bigEnd()[1]; ++j) {
      y[j - bvalsBox.smallEnd()[1]] =
        (geom.ProbLo()[1] + (j + 0.5) * geom.CellSize(1)) * tp.turb_scale_loc;
    }

    amrex::Real z = (geom.ProbLo()[2] + (planeloc + 0.5) * geom.CellSize(2)) *
                    tp.turb_scale_loc;
    fill_turb_plane(x, y, z, v, tp);
    v.mult<amrex::RunOn::Device>(tp.turb_scale_vel);
    amrex::Box ovlp = bvalsBox & data.box();
    data.copy<amrex::RunOn::Device>(v, ovlp, 0, ovlp, dcomp, AMREX_SPACEDIM);
  }
}
