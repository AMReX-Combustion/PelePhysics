
#include "SprayJet.H"
#include <AMReX_ParmParse.H>
#include "SprayParticles.H"

// Constructor where parameters are set from input file
SprayJet::SprayJet(const std::string& jet_name, const amrex::Geometry& geom)
  : m_jetName(std::move(jet_name))
{
  std::string ppspray = "spray." + m_jetName;
  amrex::ParmParse ps(ppspray);
  std::string dist_type;
  ps.get("dist_type", dist_type);
  m_dropDist = DistBase::create(dist_type);
  m_dropDist->init(ppspray);
  m_avgDia = m_dropDist->get_avg_dia();
  std::vector<amrex::Real> jcent(AMREX_SPACEDIM);
  ps.getarr("jet_cent", jcent);
  std::vector<amrex::Real> jnorm(AMREX_SPACEDIM);
  ps.getarr("jet_norm", jnorm);
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    m_cent[dir] = jcent[dir];
    m_norm[dir] = jnorm[dir];
  }
  amrex::Real mag = m_norm.vectorLength();
  m_norm /= mag;
  check_jet_cent(geom);
  ps.get("jet_dia", m_jetDia);
  ps.get("spread_angle", m_spreadAngle);
  if (m_spreadAngle < 0. || m_spreadAngle > 180.) {
    amrex::Abort("'spread_angle' must be between 0 and 180");
  }
  m_spreadAngle *= M_PI / 180.; // Assumes spread angle is in degrees
  ps.query("swirl_angle", m_swirlAngle);
  if (m_swirlAngle < 0. || m_swirlAngle > 90.) {
    amrex::Abort("'swirl_angle' must be between 0 and 90");
  }
  m_swirlAngle *= M_PI / 180.;
  ps.get("T", m_jetT);
  amrex::Real rho_part = 0.;
  SprayData* fdat = SprayParticleContainer::getSprayData();
  if (SPRAY_FUEL_NUM == 1) {
    m_jetY[0] = 1.;
    rho_part = fdat->rhoL(m_jetT, 0);
  } else {
    std::vector<amrex::Real> in_Y_jet(SPRAY_FUEL_NUM, 0.);
    ps.getarr("Y", in_Y_jet);
    amrex::Real sumY = 0.;
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      m_jetY[spf] = in_Y_jet[spf];
      sumY += in_Y_jet[spf];
      rho_part += m_jetY[spf] / fdat->rhoL(m_jetT, spf);
    }
    rho_part = 1. / rho_part;
    if (std::abs(sumY - 1.) > 1.E-8) {
      amrex::Abort(ppspray + ".Y must sum to 1");
    }
  }
  // If a rate shape profile is generated at
  // https://www.cmt.upv.es/#/ecn/download/InjectionRateGenerator/InjectionRateGenerator,
  // it can be provided here for use
  if (ps.contains("roi_file")) {
    amrex::Real cd = 0.;
    ps.get("discharge_coeff", cd);
    std::string roi_file;
    ps.get("roi_file", roi_file);
    readROI(roi_file, rho_part, cd);
    m_useROI = true;
  } else {
    ps.query("start_time", m_startTime);
    ps.query("end_time", m_endTime);
    ps.get("jet_vel", m_jetVel);
    m_maxJetVel = m_jetVel;
    ps.get("mass_flow_rate", m_massFlow);
    if (m_jetVel < 0. || m_massFlow < 0.) {
      amrex::Abort("Jet u, T, and mdot must be greater than 0");
    }
  }
  ps.query("hollow_spray", m_hollowSpray);
  if (m_hollowSpray) {
    ps.query("hollow_spread", m_hollowSpread);
    if (m_hollowSpread < 0. || m_hollowSpread > 180.) {
      amrex::Abort("'hollow_spread' must be between 0 and 180");
    }
    m_hollowSpread *= M_PI / 180.;
  }
}

// Constructor for assigning parameters directly
SprayJet::SprayJet(
  const std::string& jet_name,
  const amrex::Geometry& geom,
  const amrex::RealVect jet_cent,
  const amrex::RealVect jet_norm,
  const amrex::Real spread_angle,
  const amrex::Real jet_dia,
  const amrex::Real jet_vel,
  const amrex::Real mass_flow,
  const amrex::Real jet_temp,
  const amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> jet_Y,
  const std::string& dist_type,
  const amrex::Real start_time,
  const amrex::Real end_time,
  const amrex::Real phi_swirl,
  bool hollow_spray,
  const amrex::Real hollow_spread)
  : m_jetName(std::move(jet_name)),
    m_norm(jet_norm),
    m_cent(jet_cent),
    m_spreadAngle(spread_angle * M_PI / 180.),
    m_swirlAngle(phi_swirl * M_PI / 180.),
    m_jetDia(jet_dia),
    m_jetVel(jet_vel),
    m_maxJetVel(jet_vel),
    m_massFlow(mass_flow),
    m_jetT(jet_temp),
    m_startTime(start_time),
    m_endTime(end_time),
    m_hollowSpray(hollow_spray),
    m_hollowSpread(hollow_spread * M_PI / 180.)
{
  if (m_spreadAngle < 0. || m_spreadAngle > M_PI) {
    amrex::Abort("Spread angle must be between 0 and 180");
  }
  if (m_swirlAngle < 0. || m_swirlAngle > 0.5 * M_PI) {
    amrex::Abort("Swirl angle must be between 0 and 90");
  }
  if (m_hollowSpread < 0. || m_hollowSpread > M_PI) {
    amrex::Abort("Hollow spread must be between 0 and 180");
  }
  if (m_jetVel < 0. || m_massFlow < 0. || m_jetT < 0.) {
    amrex::Abort("Jet u, T, and mdot must be greater than 0");
  }
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    m_jetY[spf] = jet_Y[spf];
  }
  m_dropDist = DistBase::create(dist_type);
  std::string ppspray = "spray";
  m_dropDist->init(ppspray);
  m_avgDia = m_dropDist->get_avg_dia();
  amrex::Real mag = m_norm.vectorLength();
  m_norm /= mag;
  check_jet_cent(geom);
}

// Default get_new_particle
bool
SprayJet::get_new_particle(
  const amrex::Real time,
  const amrex::Real& phi_radial,
  const amrex::Real& cur_radius,
  amrex::Real& umag,
  amrex::Real& theta_spread,
  amrex::Real& phi_swirl,
  amrex::Real& dia_part,
  amrex::Real& T_part,
  amrex::Real* Y_part)
{
#if AMREX_SPACEDIM == 3
  // Percentage of droplet along radius
  // In 3D, cur_radius is from [0, jetDia/2]
  amrex::Real radp = 2. * cur_radius / m_jetDia;
  theta_spread = radp * m_spreadAngle / 2.;
  if (m_hollowSpray) {
    amrex::Real rand = amrex::Random() - 0.5;
    theta_spread += m_hollowSpread * rand;
  }
#else
  // In 2D, cur_radius is from [-jetDia/2, jetDia/2]
  amrex::Real radp = cur_radius / m_jetDia + 0.5;
  theta_spread = -(radp - 0.5) * m_spreadAngle;
#endif
  // Provide no swirling component
  phi_swirl = m_swirlAngle;
  dia_part = m_dropDist->get_dia();
  T_part = m_jetT;
  umag = jet_vel(time);
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    Y_part[spf] = m_jetY[spf];
  }
  amrex::ignore_unused(phi_radial);
  return true;
}

std::string
read_inject_file(std::ifstream& in)
{
  return static_cast<std::stringstream const&>(
           std::stringstream() << in.rdbuf())
    .str();
}

void
SprayJet::readROI(
  const std::string roi_file, const amrex::Real rho_part, const amrex::Real cd)
{
  amrex::Real jet_area = 0.25 * M_PI * std::pow(m_jetDia, 2);
  std::string firstline, remaininglines;

  std::ifstream infile(roi_file);
  const std::string memfile = read_inject_file(infile);
  if (!amrex::FileSystem::Exists(roi_file)) {
    amrex::Abort("ROI file does not exists");
  }
  infile.close();
  std::istringstream iss(memfile);

  std::getline(iss, firstline);
  // Time is given in ms
  amrex::Real tconv = 1.E-3;
#ifdef PELELM_USE_SPRAY
  amrex::Real mconv = 1.E-3;
#else
  amrex::Real mconv = 1.;
#endif
  m_maxJetVel = 0.;
  bool jet_started = false;
  amrex::Real prev_time = 0.;
  while (std::getline(iss, remaininglines)) {
    // Replace any semi-colons with spaces
    std::size_t pos = remaininglines.find(";");
    while (pos != std::string::npos) {
      remaininglines.replace(pos, 1, " ");
      pos = remaininglines.find(";");
    }
    std::istringstream sinput(remaininglines);
    amrex::Real time, mdot;
    sinput >> time;
    sinput >> mdot;
    time *= tconv;
    mdot *= mconv;
    if (mdot > 0.) {
      if (!jet_started) {
        jet_started = true;
        m_startTime = time;
        inject_time.push_back(prev_time);
        inject_mass.push_back(0.);
        inject_vel.push_back(0.);
      }
      amrex::Real vel = mdot / (rho_part * jet_area * cd);
      m_maxJetVel = amrex::max(m_maxJetVel, vel);
      inject_time.push_back(time);
      inject_mass.push_back(mdot);
      inject_vel.push_back(vel);
      m_endTime = time;
    } else if (jet_started) {
      m_endTime = time;
      inject_time.push_back(time);
      inject_mass.push_back(0.);
      inject_vel.push_back(0.);
      break;
    }
    prev_time = time;
  }
}
