
#include "SprayParticles.H"
#include "SBData.H"
#include "Distributions.H"

using namespace amrex;

// Find tangents along surface
void
find_tangents(
  const RealVect& pvel,
  RealVect& tanPsi,
  const RealVect& norm,
  RealVect& tanBeta)
{
#if AMREX_SPACEDIM == 3
  amrex::RealVect testvec = -pvel;
  // Check if directions of norm and velocity are the same
  if (testvec.crossProduct(norm) == RealVect::TheZeroVector()) {
    // If so, pick an arbitrary direction
    testvec = {norm[1], norm[2], norm[0]};
  }
  tanPsi = testvec.crossProduct(norm);
  tanPsi /= tanPsi.vectorLength();
  tanBeta = tanPsi.crossProduct(norm);
#else
  amrex::ignore_unused(pvel, tanPsi);
  tanBeta[0] = -norm[1];
  tanBeta[1] = norm[0];
#endif
}

void
SprayParticleContainer::CreateSBDroplets(
  const int Np,
  const splash_breakup* N_SB_h,
  const SBPtrs& rfh,
  const int level)
{
  ParticleLocData pld;
  SprayComps SPI = m_sprayIndx;
  const SprayData* fdat = m_sprayData;
  std::map<std::pair<int, int>, Gpu::HostVector<ParticleType>> host_particles;
  std::pair<int, int> ind(pld.m_grid, pld.m_tile);
  for (int n = 0; n < Np; n++) {
    if (N_SB_h[n] != splash_breakup::no_change) {
      RealVect normal;
      RealVect loc0;
      RealVect vel0;
      const int vn = AMREX_SPACEDIM * n;
      for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        normal[dir] = rfh.norm[vn + dir];
        loc0[dir] = rfh.loc[vn + dir];
        vel0[dir] = rfh.vel[vn + dir];
      }
      Real d0 = rfh.d0[n];
      Real dtpp = rfh.dtpp[n];
      Real phi1 = rfh.phi1[n];
      Real phi2 = rfh.phi2[n];
      Real T0 = rfh.T0[n];
      Array<Real, SPRAY_FUEL_NUM> Y0;
#if SPRAY_FUEL_NUM > 1
      Real rho_part = 0.;
      Real mu_part = 0.;
      const int vy = SPRAY_FUEL_NUM * n;
      for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
        Y0[spf] = rfh.Y0[vy + spf];
        rho_part += Y0[spf] / fdat->rho[spf];
        mu_part += Y0[spf] * fdat->mu[spf];
      }
      rho_part = 1. / rho_part;
#else
      Real rho_part = fdat->rho[0];
      Real mu_part = fdat->mu[0];
      Y0[0] = 1.;
#endif
      const Real sigma = fdat->sigma;
      Real pmass = M_PI / 6. * rho_part * std::pow(d0, 3);
      Real U0mag = vel0.vectorLength();

      if (
        N_SB_h[n] == splash_breakup::splash_splash ||
        N_SB_h[n] == splash_breakup::splash_thermal_breakup) {
        // tanPsi: vector tangent to wall normal and velocity direction
        // tanBeta: vector tangent to wall normal in plane with velocity
        RealVect tanPsi, tanBeta;
        find_tangents(vel0, tanPsi, normal, tanBeta);
        Real Kv = phi1;
        Real ms = phi2;
        Real U0norm = normal.dotProduct(vel0);
        Real alpha = std::asin(amrex::Math::abs(U0norm) / U0mag);
        Real U0tan = std::sqrt(U0mag * U0mag - U0norm * U0norm);
        Real Unmag = 0.0028 * Kv * std::exp(-0.0062 * Kv) * std::abs(U0norm);
        Real Utmag = 0.0065 * Kv * std::exp(-0.004 * Kv) * std::abs(U0norm);

        // Number of splashed droplets from the Marmanis-Thoroddsen model
        amrex::Real Nsdrops = std::ceil(
          std::abs(U0norm) / (20. * std::sqrt(mu_part / rho_part)) *
          std::pow(6. * M_PI * pmass / sigma, 0.25));
        int Nsint = static_cast<int>(Nsdrops);
        Real tanterm = 0.4 * std::tan(M_PI / 2. - alpha);
        // Average velocity of all splashed droplets
        RealVect avg_vel(RealVect::TheZeroVector());
        // Sum splashed mass
        Real sum_mass = 0.;
        for (int new_parts = 0; new_parts < Nsint; ++new_parts) {
          Real psi = amrex::Random() * M_PI;
          ParticleType p;
          p.id() = ParticleType::NextID();
          p.cpu() = ParallelDescriptor::MyProc();
          Real new_mass = ms * (1. + tanterm * std::cos(psi)) / Nsdrops;
          Real dia_part = std::cbrt(6. * new_mass / (M_PI * rho_part));
          p.rdata(SPI.pstateDia) = dia_part;
          Real psi_sign = std::copysign(1., 0.5 - amrex::Random());
          AMREX_D_TERM(
            Real un = Unmag;
            , Real utBeta = std::cos(psi / 3.) * U0tan + Utmag * std::cos(psi);
            , Real utPsi = psi_sign * std::sin(psi) * (0.2 * U0tan + Utmag););
          for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            Real pvel = AMREX_D_TERM(
              un * normal[dir], +utBeta * tanBeta[dir], +utPsi * tanPsi[dir]);
            avg_vel[dir] += pvel / Nsdrops;
            p.pos(dir) = loc0[dir] + dtpp * pvel;
            p.rdata(SPI.pstateVel + dir) = pvel;
          }
          for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
            p.rdata(SPI.pstateY + spf) = Y0[spf];
          }
          p.rdata(SPI.pstatePb) = 0.;
          p.rdata(SPI.pstatePbdot) = 0.;
          bool where = Where(p, pld);
          if (!where) {
            amrex::Abort("Bad reflected particle");
          }
          p.rdata(SPI.pstateT) = T0;
          host_particles[ind].push_back(p);
          sum_mass += new_mass;
        }
        // Remaining unsplashed mass
        Real rem_mass = pmass - sum_mass;
        // Original droplet is recreated with new values
        ParticleType p;
        p.id() = ParticleType::NextID();
        p.cpu() = ParallelDescriptor::MyProc();
        for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
          p.rdata(SPI.pstateY + spf) = Y0[spf];
        }
        p.rdata(SPI.pstateT) = T0;
        p.rdata(SPI.pstatePb) = 0.;
        p.rdata(SPI.pstatePbdot) = 0.;
        // If droplet splashing is not thermally breakup, center droplet also
        // reflects
        if (N_SB_h[n] == splash_breakup::splash_splash) {
          Real dia_rem = std::cbrt(6. * rem_mass / (M_PI * rho_part));
          p.rdata(SPI.pstateDia) = dia_rem;
          for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
            p.pos(dir) = loc0[dir] + dtpp * avg_vel[dir];
            p.rdata(SPI.pstateVel + dir) = avg_vel[dir];
          }
          bool where = Where(p, pld);
          if (!where) {
            amrex::Abort("Bad reflected particle");
          }
          host_particles[ind].push_back(p);
        } else {
          // Otherwise, droplet forms a film
          // TODO: Add wall film
          // for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
          //   p.pos(dir) = loc0[dir];
          //   p.rdata(SPI.pstateVel + dir) = 0.;
          // }
          // host_particles[ind].push_back(p);
        }
      } else if (N_SB_h[n] == splash_breakup::breakup) {
        // TODO: Add distribution for radii
        Real r32 = phi1;
        Real d32 = 2. * r32;
        Real Utan = phi2;
        Real bmass = M_PI / 6. * rho_part * std::pow(d32, 3);
        int Nsint = static_cast<int>(pmass / bmass);
        auto newbmass = pmass / static_cast<Real>(Nsint);
        Real newd32 = std::cbrt(6. * newbmass / (M_PI * rho_part));
#if AMREX_SPACEDIM == 3
        RealVect testvec(normal[1], normal[2], normal[0]);
        RealVect tanPsi = testvec.crossProduct(normal);
        RealVect tanBeta = tanPsi.crossProduct(normal);
#else
        RealVect tanBeta(normal[1], normal[0]);
#endif
        for (int new_parts = 0; new_parts < Nsint; ++new_parts) {
          Real rand = amrex::Random();
          ParticleType p;
          p.id() = ParticleType::NextID();
          p.cpu() = ParallelDescriptor::MyProc();
          p.rdata(SPI.pstateDia) = d32;
          p.rdata(SPI.pstateT) = T0;
          for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
            p.rdata(SPI.pstateY + spf) = Y0[spf];
          }
          p.rdata(SPI.pstatePb) = 0.;
          p.rdata(SPI.pstatePbdot) = 0.;
          for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
#if AMREX_SPACEDIM == 3
            Real psi = rand * 2. * M_PI;
            Real pvel = vel0[dir] + Utan * (std::sin(psi) * tanPsi[dir] +
                                            std::cos(psi) * tanBeta[dir]);
#else
            Real sgn = std::copysign(1., 0.5 - rand);
            Real pvel = vel0[dir] + sgn * Utan * tanBeta[dir];
#endif
            p.pos(dir) = loc0[dir] + dtpp * pvel;
            p.rdata(SPI.pstateVel + dir) = pvel;
          }
          bool where = Where(p, pld);
          if (!where) {
            amrex::Abort("Bad breakup particle");
          }
          host_particles[ind].push_back(p);
        }
      }
    }
  }
  for (auto& kv : host_particles) {
    auto grid = kv.first.first;
    auto tile = kv.first.second;
    const auto& src_tile = kv.second;
    auto& dst_tile = GetParticles(level)[std::make_pair(grid, tile)];
    auto old_size = dst_tile.GetArrayOfStructs().size();
    auto new_size = old_size + src_tile.size();
    dst_tile.resize(new_size);
    // Copy the AoS part of the host particles to the GPU
    Gpu::copy(
      Gpu::hostToDevice, src_tile.begin(), src_tile.end(),
      dst_tile.GetArrayOfStructs().begin() + old_size);
  }
}
