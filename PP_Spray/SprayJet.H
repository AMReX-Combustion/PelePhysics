#ifndef SPRAYJET_H
#define SPRAYJET_H

#include "DistBase.H"

class SprayJet
{
public:
  // Constructor for assigning parameters directly
  SprayJet(
    const amrex::Geometry& geom,
    const amrex::RealVect jet_cent,
    const amrex::RealVect jet_norm,
    const amrex::Real spread_angle,
    const amrex::Real jet_dia,
    const amrex::Real jet_vel,
    const amrex::Real mass_flow,
    const amrex::Real jet_temp,
    const amrex::Real jet_Y[SPRAY_FUEL_NUM],
    const std::string dist_type,
    const amrex::Real start_time = 0.,
    const amrex::Real end_time = 1.E8,
    bool hollow_spray = false)
    : m_norm(jet_norm),
      m_cent(jet_cent),
      m_spreadAngle(spread_angle),
      m_jetDia(jet_dia),
      m_jetVel(jet_vel),
      m_massFlow(mass_flow),
      m_jetT(jet_temp),
      m_startTime(start_time),
      m_endTime(end_time),
      m_hollowSpray(hollow_spray)
  {
    m_spreadAngle *= M_PI / 180.;
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      m_jetY[spf] = jet_Y[spf];
    }
    m_dropDist = DistBase::create(dist_type);
    std::string ppspray = "spray";
    m_dropDist->init(ppspray);
    check_jet_cent(geom);
  }

  // Constructor assumes parameters are set from input file
  SprayJet(const std::string jet_name, const amrex::Geometry& geom)
  {
    std::string ppspray = "spray." + jet_name;
    amrex::ParmParse ps(ppspray);
    /*
      Some parameters that likely remain constant over time are required as
      inputs: jet normal, jet center, jet diameter, spread angle, injection
      droplet temperature, and injection droplet mass fraction Other parameters
      can be set as constant but might vary with time: jet velocity, mass flow
      rate
    */
    std::string dist_type;
    ps.get("dist_type", dist_type);
    m_dropDist = DistBase::create(dist_type);
    m_dropDist->init(ppspray);
    std::vector<amrex::Real> jcent(AMREX_SPACEDIM);
    ps.getarr("jet_cent", jcent);
    std::vector<amrex::Real> jnorm(AMREX_SPACEDIM);
    ps.getarr("jet_norm", jnorm);
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      m_cent[dir] = jcent[dir];
      m_norm[dir] = jnorm[dir];
    }
    check_jet_cent(geom);
    ps.get("jet_dia", m_jetDia);
    ps.get("spread_angle", m_spreadAngle);
    ps.get("T", m_jetT);
    m_spreadAngle *= M_PI / 180.; // Assumes spread angle is in degrees
    std::vector<amrex::Real> in_Y_jet(SPRAY_FUEL_NUM, 0.);
    if (SPRAY_FUEL_NUM == 1) {
      m_jetY[0] = 1.;
    } else {
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
    ps.query("jet_vel", m_jetVel);
    ps.query("mass_flow_rate", m_massFlow);
    ps.query("hollow_spray", m_hollowSpray);
  }

  SprayJet(SprayJet const&) = delete;
  SprayJet& operator=(SprayJet const&) = delete;

  inline const amrex::RealVect& jet_norm() const { return m_norm; }
  inline const amrex::RealVect& jet_cent() const { return m_cent; }
  inline const amrex::Real& jet_vel() const { return m_jetVel; }
  inline const amrex::Real& jet_dia() const { return m_jetDia; }
  inline const amrex::Real& spread_angle() const { return m_spreadAngle; }
  inline const amrex::Real& mass_flow_rate() const { return m_massFlow; }
  inline const amrex::Real& Y(const int spf) const { return m_jetY[spf]; }
  inline const amrex::Real& T() const { return m_jetT; }
  inline const amrex::Real& start_time() const { return m_startTime; }
  inline const amrex::Real& end_time() const { return m_endTime; }

  // Call this before using spray jet
  bool jet_active(const amrex::Real time) const
  {
    if (m_jetVel < 0.) {
      amrex::Abort("Jet velocity is not set");
    } else if (m_massFlow < 0.) {
      amrex::Abort("Jet mass flow rate is not set");
    }
    bool check = true;
    if (time < m_startTime || time > m_endTime) {
      check = false;
    }
    return check;
  }

  inline bool hollow_spray() const { return m_hollowSpray; }
  inline amrex::Real get_dia() const { return m_dropDist->get_dia(); }
  inline amrex::Real get_avg_dia() const { return m_dropDist->get_avg_dia(); }

  void set_jet_vel(amrex::Real jet_vel) { m_jetVel = jet_vel; }
  void set_mass_flow(amrex::Real mass_flow) { m_massFlow = mass_flow; }
  void set_jet_norm(amrex::RealVect jet_norm) { m_norm = jet_norm; }
  void set_jet_cent(amrex::RealVect jet_cent, const amrex::Geometry& geom)
  {
    m_cent = jet_cent;
    check_jet_cent(geom);
  }

  void check_jet_cent(const amrex::Geometry& geom) const
  {
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      if (m_cent[dir] > geom.ProbHi(dir) || m_cent[dir] < geom.ProbLo(dir)) {
        amrex::Abort("Jet center outside of problem domain");
      }
    }
  }

  void set_mass_frac(const amrex::Real* Y_jet)
  {
    amrex::Real sumtest = 0.;
    for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
      m_jetY[spf] = Y_jet[spf];
      sumtest += m_jetY[spf];
    }
    if (amrex::Math::abs(1. - sumtest) > 1.E-8) {
      amrex::Abort("Jet mass fractions must sum to 1");
    }
  }

  void set_start_time(amrex::Real start_time) { m_startTime = start_time; }
  void set_end_time(amrex::Real end_time) { m_endTime = end_time; }

  void reset_sum()
  {
    m_sumInjMass = -1;
    m_sumInjTime = -1;
  }

  amrex::Real m_sumInjMass = -1.;
  amrex::Real m_sumInjTime = -1.;

protected:
  // Member data
  amrex::RealVect m_norm = amrex::RealVect::TheZeroVector();
  amrex::RealVect m_cent = amrex::RealVect::TheZeroVector();
  amrex::Real m_spreadAngle = 0.;
  amrex::Real m_jetDia = -1.;
  amrex::Real m_jetVel = -1.;
  amrex::Real m_massFlow = -1.;
  amrex::Real m_jetT = -1.;
  amrex::Real m_jetY[SPRAY_FUEL_NUM] = {0.};
  std::unique_ptr<DistBase> m_dropDist;
  amrex::Real m_startTime = 0.;
  amrex::Real m_endTime = 1.E8;
  bool m_hollowSpray = false;
};
#endif
