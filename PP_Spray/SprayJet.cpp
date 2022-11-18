
#include "SprayJet.H"
#include <AMReX_ParmParse.H>

// Constructor where parameters are set from inputfile
SprayJet::SprayJet(const std::string& jet_name, const amrex::Geometry& geom)
{
  std::string ppspray = "spray." + jet_name;
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
  m_spreadAngle *= M_PI / 180.; // Assumes spread angle is in degrees
  ps.query("swirl_angle", m_swirlAngle);
  m_swirlAngle *= M_PI / 180.;
  ps.get("T", m_jetT);
  if (SPRAY_FUEL_NUM == 1) {
    m_jetY[0] = 1.;
  } else {
    std::vector<amrex::Real> in_Y_jet(SPRAY_FUEL_NUM, 0.);
    ps.getarr("Y", in_Y_jet);
    amrex::Real sumY = 0.;
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      m_jetY[spf] = in_Y_jet[spf];
      sumY += in_Y_jet[spf];
    }
    if (amrex::Math::abs(sumY - 1.) > 1.E-8) {
      amrex::Abort(ppspray + ".Y must sum to 1");
    }
  }
  ps.query("start_time", m_startTime);
  ps.query("end_time", m_endTime);
  ps.get("jet_vel", m_jetVel);
  ps.get("mass_flow_rate", m_massFlow);
  ps.query("hollow_spray", m_hollowSpray);
}

// Constructor for assigning parameters directly
SprayJet::SprayJet(
  const amrex::Geometry& geom,
  const amrex::RealVect jet_cent,
  const amrex::RealVect jet_norm,
  const amrex::Real spread_angle,
  const amrex::Real jet_dia,
  const amrex::Real jet_vel,
  const amrex::Real mass_flow,
  const amrex::Real jet_temp,
  const amrex::Real jet_Y[SPRAY_FUEL_NUM],
  const std::string& dist_type,
  const amrex::Real start_time,
  const amrex::Real end_time,
  const amrex::Real phi_swirl,
  bool hollow_spray)
  : m_norm(jet_norm),
    m_cent(jet_cent),
    m_spreadAngle(spread_angle * M_PI / 180.),
    m_swirlAngle(phi_swirl * M_PI / 180.),
    m_jetDia(jet_dia),
    m_jetVel(jet_vel),
    m_massFlow(mass_flow),
    m_jetT(jet_temp),
    m_startTime(start_time),
    m_endTime(end_time),
    m_hollowSpray(hollow_spray)
{
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    m_jetY[spf] = jet_Y[spf];
  }
  m_dropDist = DistBase::create(dist_type);
  std::string ppspray = "spray";
  m_dropDist->init(ppspray);
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
#else
  // In 2D, cur_radius is from [-jetDia/2, jetDia/2]
  amrex::Real radp = cur_radius / m_jetDia + 0.5;
  theta_spread = -(radp - 0.5) * m_spreadAngle;
#endif
  // Provide no swirling component
  phi_swirl = m_swirlAngle;
  dia_part = m_dropDist->get_dia();
  T_part = m_jetT;
  umag = m_jetVel;
  for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
    Y_part[spf] = m_jetY[spf];
  }
  amrex::ignore_unused(phi_radial, time);
  return true;
}
