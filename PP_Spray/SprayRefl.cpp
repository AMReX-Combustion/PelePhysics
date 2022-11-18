
#include "SprayParticles.H"
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
  tanPsi = -pvel.crossProduct(norm);
  tanBeta = tanPsi.crossProduct(norm);
  tanPsi /= tanPsi.vectorLength();
  tanBeta /= tanBeta.vectorLength();
#else
  amrex::ignore_unused(pvel, tanPsi);
  tanBeta[0] = -norm[1];
  tanBeta[1] = norm[0];
#endif
}

void
SprayParticleContainer::CreateReflectedDroplets(
  const int Np, const int* N_refl_h, const ReflPtrs& rfh, const int level)
{
  ParticleLocData pld;
  SprayComps SPI = m_sprayIndx;
  const SprayData* fdat = m_sprayData;
  std::map<std::pair<int, int>, Gpu::HostVector<ParticleType>> host_particles;
  std::pair<int, int> ind(pld.m_grid, pld.m_tile);
  for (int n = 0; n < Np; n++) {
    if (N_refl_h[n] > 0) {
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
      Real Kv = rfh.Kv[n];
      Real ms = rfh.ms[n];
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
      Real alpha = rfh.alpha[n];

      Real U0mag = vel0.vectorLength();
      Real U0norm = normal.dotProduct(vel0);
      Real U0tan = std::sqrt(U0mag * U0mag - U0norm * U0norm);
      Real Unmag = 0.0028 * Kv * std::exp(-0.0062 * Kv) * std::abs(U0norm);
      Real Utmag = 0.0065 * Kv * std::exp(-0.004 * Kv) * std::abs(U0norm);

      // tanPsi: vector tangent to wall normal and velocity direction
      // tanBeta: vector tangent to wall normal in plane with velocity
      RealVect tanPsi, tanBeta;
      find_tangents(vel0, tanPsi, normal, tanBeta);

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
          p.rdata(SPI.pstateY) = Y0[spf];
        }
        p.rdata(SPI.pstateTime) = 0.;
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
      // If droplet breaksup, center droplet also reflects
      if (N_refl_h[n] == 1) {
        Real dia_rem = std::cbrt(6. * rem_mass / (M_PI * rho_part));
        p.rdata(SPI.pstateDia) = dia_rem;
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
          p.pos(dir) = loc0[dir] + dtpp * avg_vel[dir];
          p.rdata(SPI.pstateVel + dir) = avg_vel[dir];
        }
        p.rdata(SPI.pstateT) = T0;
        for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
          p.rdata(SPI.pstateY) = Y0[spf];
        }
        p.rdata(SPI.pstateTime) = 0.;
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
