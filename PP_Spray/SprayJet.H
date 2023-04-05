#ifndef SPRAYJET_H
#define SPRAYJET_H

#include "DistBase.H"
#include <AMReX_RealVect.H>
#include <AMReX_Geometry.H>

class SprayJet
{
public:
  // Default constructor
  SprayJet() = default;

  // Constructor for assigning parameters directly
  SprayJet(
    const std::string jet_name,
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
    const amrex::Real start_time = 0.,
    const amrex::Real end_time = 1.E8,
    const amrex::Real phi_swirl = 0.,
    bool hollow_spray = false,
    const amrex::Real hollow_spread = 0.);

  // Constructor where parameters are set from input file
  SprayJet(const std::string jet_name, const amrex::Geometry& geom);

  virtual ~SprayJet() = default;

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
  inline const std::string& jet_name() const { return m_jetName; }
  inline int Proc() const { return m_proc; }

  // Call this before using spray jet
  bool jet_active(const amrex::Real time) const
  {
    if (m_massFlow < 0.) {
      amrex::Abort("Jet mass flow rate is not set");
    } else if (m_jetDia < 0.) {
      amrex::Abort("Jet diameter is not set");
    } else if (m_avgDia < 0.) {
      amrex::Abort("Average droplet diameter is not set");
    } else if (m_jetT < 0.) {
      amrex::Abort("Droplet temperature is not set");
    }
    bool check = true;
    if (time < m_startTime || time > m_endTime) {
      check = false;
    }
    return check;
  }

  inline bool hollow_spray() const { return m_hollowSpray; }

  /// Returns the average values given for the jet
  amrex::Real get_avg_dia() const { return m_avgDia; }

  amrex::Real get_avg_T() const { return m_jetT; }

  amrex::Real get_avg_Y(const int spf) const { return m_jetY[spf]; }

  void set_jet_vel(amrex::Real jet_vel) { m_jetVel = jet_vel; }

  void set_mass_flow(amrex::Real mass_flow) { m_massFlow = mass_flow; }

  void set_jet_norm(amrex::RealVect jet_norm)
  {
    m_norm = jet_norm;
    amrex::Real mag = m_norm.vectorLength();
    m_norm /= mag;
  }

  void set_inj_proc(int inj_proc)
  {
    m_proc = inj_proc;
    if (m_proc < 0 || m_proc > amrex::ParallelDescriptor::NProcs()) {
      amrex::Abort("SprayJet proc not valid");
    }
  }

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
  void set_num_ppp(amrex::Real num_ppp) { m_numPPP = num_ppp; }

  void reset_sum()
  {
    m_sumInjMass = -1;
    m_sumInjTime = -1;
  }

  /**
     Function for creating new particle parameters based on jet location. Note
     all angles must be in radians.
     @param[in] time Current solution time
     @param[in] phi_radial Azimuthal angle of particle location in jet CS
     @param[in] cur_radius Current radial location in jet CS
     @param[out] umag Magnitude of velocity
     @param[out] theta_spread Spread angle for velocity relative to jet norm
     @param[out] phi_swirl Azimuthal angle for tangential velocity component
     @param[out] dia_part Droplet diameter
     @param[out] T_part Droplet temperature
     @param[out] Y_part Droplet mass fractions
  */
  virtual bool get_new_particle(
    const amrex::Real time,
    const amrex::Real& phi_radial,
    const amrex::Real& cur_radius,
    amrex::Real& umag,
    amrex::Real& theta_spread,
    amrex::Real& phi_swirl,
    amrex::Real& dia_part,
    amrex::Real& T_part,
    amrex::Real* Y_part);

  /**
     Solve for transformed location and velocity based on provided angles and
     radius.
     @param[in] theta_spread Spread angle for velocity relative to jet norm
     @param[in] phi_radial Azimuthal angle of particle location in jet CS
     @param[in] cur_radius Current radial location of particle
     @param[in] umag Particle velocity magnitude
     @param[in] phi_swirl Azimuthal angle for tangential velocity component
     @param[out] part_vel Particle velocity in jet CS to be solved for
     @param[out] part_loc Particle location in jet CS to be solved for
   */
  void transform_loc_vel(
    const amrex::Real theta_spread,
    const amrex::Real phi_radial,
    const amrex::Real cur_radius,
    const amrex::Real umag,
    const amrex::Real phi_swirl,
    amrex::RealVect& part_vel,
    amrex::RealVect& part_loc)
  {
#if AMREX_SPACEDIM == 3
    amrex::Real norm_mag = m_norm.vectorLength();
    amrex::Real theta_jet = std::acos(m_norm[2] / norm_mag);
    amrex::Real phi_jet =
      std::atan2(m_norm[1] / norm_mag, m_norm[0] / norm_mag);
    amrex::Real sp1 = std::sin(phi_jet);
    amrex::Real cp1 = std::cos(phi_jet);
    amrex::Real sp2 = std::sin(phi_radial);
    amrex::Real cp2 = std::cos(phi_radial);
#else
    amrex::Real theta_jet = std::atan2(m_norm[1], m_norm[0]) + M_PI / 2.;
    amrex::ignore_unused(phi_radial, phi_swirl);
#endif
    amrex::Real st1 = std::sin(theta_jet);
    amrex::Real ct1 = std::cos(theta_jet);
    amrex::Real st2 = std::sin(theta_spread);
    amrex::Real ct2 = std::cos(theta_spread);
#if AMREX_SPACEDIM == 3
    amrex::RealVect dp(AMREX_D_DECL(
      cp1 * cp2 * ct1 - sp1 * sp2, sp1 * cp2 * ct1 + cp1 * sp2,
      -std::sin(theta_jet) * cp2));
    // Add phi_swirl for velocity
    amrex::Real phivel = phi_radial + phi_swirl;
    sp2 = std::sin(phivel);
    cp2 = std::cos(phivel);
    amrex::Real v1 = st1 * ct2 + st2 * cp2 * ct1;
    part_vel = {
      cp1 * v1 - sp1 * sp2 * st2, sp1 * v1 + sp2 * st2 * cp1,
      ct1 * ct2 - st1 * st2 * cp2};
#else
    amrex::RealVect dp(ct1, st1);
    part_vel = {st1 * ct2 - st2 * ct1, -ct1 * ct2 - st1 * st2};
#endif
    part_loc = m_cent + cur_radius * dp;
    part_vel *= umag;
  }

  amrex::Real m_sumInjMass = -1.;
  amrex::Real m_sumInjTime = -1.;
  // Minimum parcels to inject at a time
  amrex::Real m_minParcel = 1.;

protected:
  // Member data
  std::string m_jetName;
  amrex::RealVect m_norm = amrex::RealVect::TheZeroVector();
  amrex::RealVect m_cent = amrex::RealVect::TheZeroVector();
  amrex::Real m_spreadAngle = 0.;
  amrex::Real m_swirlAngle = 0.;
  amrex::Real m_jetDia = -1.;
  amrex::Real m_jetVel = -1.;
  amrex::Real m_massFlow = -1.;
  amrex::Real m_jetT = -1.;
  amrex::GpuArray<amrex::Real, SPRAY_FUEL_NUM> m_jetY = {{0.}};
  amrex::Real m_avgDia = -1.;
  std::unique_ptr<DistBase> m_dropDist;
  amrex::Real m_startTime = 0.;
  amrex::Real m_endTime = 1.E8;
  bool m_hollowSpray = false;
  amrex::Real m_hollowSpread = 0.;
  int m_proc = 0;
};

#endif
