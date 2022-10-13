
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
  // LogNormal un_dist;
  // Distribution for secondary droplet diameter
  LogNormal dia_dist;
  Normal v_dist;
  std::pair<int, int> ind(pld.m_grid, pld.m_tile);
  for (int n = 0; n < Np; n++) {
    if (N_refl_h[n] > 0) {
      // Represents the D32/D10 of the secondary droplets. Varies from 1.28-2 in
      // literature
      const Real nu32 = 2.;
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
      Real We = rfh.We[n];
      Real Tstar = rfh.Tstar[n];
      Real T0 = rfh.T0[n];
      Array<Real, SPRAY_FUEL_NUM> Y0;
#if SPRAY_FUEL_NUM > 1
      Real rho_part = 0.;
      Real mu_part = 0.;
      const int vy = SPRAY_FUEL_NUM * n;
      for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
        Y0[spf] = rfh.Y0[vy + spf];
        rho_part += Y_part[spf] / fdat->rho[spf];
        mu_part += Y_part[spf] * fdat->mu[spf];
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
      Real Unmax = 0.3 * U0norm;

      Real Re_L = rho_part * std::abs(U0norm) * d0 / mu_part;
      Real Kv = std::sqrt(We * std::sqrt(Re_L));
      const Real expon = 3.6 * (alpha / M_PI) * (alpha / M_PI);
      // Mean diameter of secondary droplets
      // Real d10 = d0 * 3.3 * std::exp(expon) * std::pow(We, -0.65);
      Real d10 = d0 * 2.2 * std::exp(expon) * std::pow(We, -0.36);
      dia_dist.init(1.16, 0.5);
      Real omega = M_PI * std::sqrt(std::cos(alpha) / (1. - std::cos(alpha)));
      Real expomega = 1. - std::exp(-omega);
      v_dist.init(0., 0.1 * std::abs(U0norm));
      // un_dist.init(1.16, 0.5);
      // tanPsi: vector tangent to wall normal and velocity direction
      // tanBeta: vector tangent to wall normal in plane with velocity
      RealVect tanPsi, tanBeta;
      find_tangents(vel0, tanPsi, normal, tanBeta);

      Real splash_mass = pmass;
      // If droplet splashes instead of breakup
      if (N_refl_h[n] == 2) {
        const Real B = 0.2 + 0.6 * amrex::Random();
        splash_mass *= min(1., (Tstar - 0.8) / 0.3 * (1. - B) + B);
        // TODO: Make actual wall film model
        Real mass_depo = pmass - splash_mass;
        Real dia_depo = std::cbrt(mass_depo * 6. / (M_PI * rho_part));
        // ParticleType p;
        // p.id() = ParticleType::NextID();
        // p.cpu() = ParallelDescriptor::MyProc();
        // for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        //   p.pos(dir) = loc0[dir];
        //   p.rdata(SPI.pstateVel + dir) = 0.;
        // }
        // p.rdata(SPI.pstateT) = T0;
        // p.rdata(SPI.pstateDia) = dia_depo;
        // for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
        //   p.rdata(SPI.pstateY + spf) = Y0[spf];
        // }
        // host_particles[ind].push_back(p);
      }
      Real refl_mass = 0.;
      while (refl_mass < splash_mass) {
        ParticleType p;
        p.id() = ParticleType::NextID();
        p.cpu() = ParallelDescriptor::MyProc();
        Real ln_val = dia_dist.get_dia();
        Real dia_part = d10 * ln_val;
        Real new_mass = M_PI / 6. * rho_part * std::pow(dia_part, 3);
        if (refl_mass + new_mass > splash_mass) {
          new_mass = splash_mass - refl_mass;
          dia_part = std::cbrt(6. * new_mass / (rho_part * M_PI));
        }
        p.rdata(SPI.pstateDia) = dia_part;
        Real psi = 0.;
        if (alpha * 180. / M_PI > 3. && AMREX_SPACEDIM == 3) {
          Real psi_sign = std::copysign(1., 0.5 - amrex::Random());
          Real rand1 = amrex::Random();
          psi = psi_sign * M_PI / omega * std::log(1. - rand1 * expomega);
          if (alpha * 180. / M_PI > 80.) {
            psi = psi_sign * rand1 * M_PI;
          }
        }
        Real Unorm = -Unmax * ln_val;
        Real Utan = std::abs(0.12 * U0norm) + std::abs(v_dist.get_dia());
        AMREX_D_TERM(Real un = Unorm;
                     , Real utBeta = Utan * (std::cos(psi) + 0.8 * U0tan);
                     , Real utPsi = Utan * std::sin(psi););
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
          Real pvel = AMREX_D_TERM(
            un * normal[dir], +utBeta * tanBeta[dir], +utPsi * tanPsi[dir]);
          p.pos(dir) = loc0[dir] + dtpp * pvel;
          p.rdata(SPI.pstateVel + dir) = pvel;
        }
        for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
          p.rdata(SPI.pstateY) = Y0[spf];
        }
        bool where = Where(p, pld);
        if (!where) {
          amrex::Abort("Bad reflected particle");
        }
        p.rdata(SPI.pstateT) = T0;
        host_particles[ind].push_back(p);
        refl_mass += M_PI / 6. * rho_part * std::pow(dia_part, 3);
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
