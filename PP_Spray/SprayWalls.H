#ifndef _SPRAYWALLS_H_
#define _SPRAYWALLS_H_

enum splash_type
{
  rebound = 0,
  deposit,
  splash,
  thermal_breakup
};

// Determine the type of particle impact at the wall
// Return 0 - rebound, 1 - deposit, 2 - splash, 3 - thermal breakup
splash_type splash_criteria (const Real Kv,
                             const Real Tstar,
                             const Real alpha)
{
  Real Kcrit = 20. + 2.*alpha/M_PI*20.;
  if (Tstar < 1.1) {
    Kcrit = 130.;
    if (Tstar < 1.)
      Kcrit = 54. + 76.*std::exp(13.*(Tstar - 1.));
    if (Kv < Kcrit) return splash_type::deposit;
    else return splash_type::splash;
  }
  if (Kv < Kcrit) return splash_type::rebound;
  return splash_type::thermal_breakup;
}

// Check if tile is adjacent to non-periodic boundaries
bool tile_at_bndry (const Box& in_box,
                    const IntVect& bndry_lo,
                    const IntVect& bndry_hi,
                    const Box& domain)
{
  Box testBox(in_box);
  testBox.grow(1);
  if (domain.contains(testBox)) return false;
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    Box boxlo(in_box);
    boxlo.growLo(dir, 1);
    Box boxhi(in_box);
    boxhi.growHi(dir, 1);
    if (!domain.contains(boxlo) && bndry_lo[dir] != 0) {
      return true;
    } else if (!domain.contains(boxhi) && bndry_hi[dir] != 0) {
      return true;
    }
  }
  return false;
}

// Find tangents along surface
void find_tangents (const RealVect& testvec,
                    const RealVect& norm,
                    RealVect& tanBeta,
                    RealVect& tanPsi)
{
#if AMREX_SPACEDIM == 3
  tanPsi = testvec.crossProduct(norm);
  tanBeta = tanPsi.crossProduct(norm);
  tanPsi /= tanPsi.vectorLength();
  tanBeta /= tanBeta.vectorLength();
#else
  tanBeta[0] = -norm[1];
  tanBeta[1] = norm[0];
#endif
}

// Compute the angles of reflecting droplets
void compute_angles (const Real& alpha,
                     const Real& Tstar,
                     const Real& We,
                     const bool dry,
                     SprayRefl& SPRF)
{
  // Inclination angle in degrees
  const Real alphad = alpha*180./M_PI;
  Real omega = 0.;
  if (alphad <= 80. && AMREX_SPACEDIM == 3)
    omega =
      std::sqrt((1. + 8.872*std::cos(1.152*alpha))/(1. - std::cos(alpha)));
  SPRF.omega = omega;
  SPRF.expomega = 1. - std::exp(-omega);
  Real beta_mean;
  if (dry)
    beta_mean = 9.3 + 0.22*alphad;
  else if (Tstar > 1.1)
    beta_mean = 0.225*alphad*
      std::exp(std::pow(0.017*alphad - 0.937, 2));
  else
    beta_mean = 0.96*alphad*std::exp(-4.5E-3*We);
  Real stdev = 4.; // 4 degrees
  // Now convert mean and stdev to log terms
  Real mean2 = beta_mean*beta_mean;
  Real st2 = stdev*stdev;
  Real term1 = std::log(beta_mean);
  Real term2 = std::log(mean2 + st2);
  SPRF.beta_mean = 2.*term1 - 0.5*term2;
  SPRF.beta_stdv = std::sqrt(amrex::max(-2.*term1 + term2, 0.));
}

void impose_wall (SprayParticleContainer::ParticleType& p,
                  SprayComps SPI,
                  SprayUnits SPU,
                  SprayData fdat,
                  const IntVect& ijkp,
                  const RealVect& dx,
                  const RealVect& dxi,
                  const RealVect& plo,
                  const RealVect& phi,
#ifdef AMREX_USE_EB
                  const bool use_EB,
                  const EBCellFlagFab& flags,
                  const CutFab& bcent,
                  const CutFab& bnorm,
#endif
                  const IntVect& bndry_lo,
                  const IntVect& bndry_hi,
                  const Real T_wall,
                  SprayRefl& SPRF,
                  RealVect& pos)
{
  const Real tolerance = std::numeric_limits<Real>::epsilon();
  Real sigma = fdat.sigma;
  // Determine the location of the particle after moving
  const RealVect lx = (p.pos() - plo)*dxi;
  const IntVect ijk = lx.floor();
  AMREX_D_TERM(const int ipn = ijk[0];, const int jpn = ijk[1];, const int kpn = ijk[2];);
  AMREX_D_TERM(const int ip = ijkp[0];, const int jp = ijkp[1];, const int kp = ijkp[2];);
  RealVect normal(RealVect::TheZeroVector()); // Normal vector to boundary face
  RealVect bcentv(RealVect::TheZeroVector()); // Point on boundary face
  IntVect bloc(ijk); // Index of nearest cell center
  bool check_part = false; // Check if particle is reflected
  // First check if particle has exited the domain through a Cartesian boundary
  for (int hilo = 0; hilo < 2; ++hilo) {
    Real fact = 1.;
    if (hilo == 1) fact = -1.;
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      Real dom_loc = plo[dir];
      int bndry_flag = bndry_lo[dir];
      if (hilo == 1) {
        dom_loc = phi[dir];
        bndry_flag = bndry_hi[dir];
      }
      if (fact*(p.pos(dir) - dom_loc) < 0.) {
        switch (bndry_flag) {
        case -1 : // Non-reflective, remove particle
          p.id() = -1;
          return;
        case 1: // Reflective, set vectors and locations
          normal[dir] = -1.;
          bcentv[dir] = fact*0.5;
          bloc[dir] = ijk[dir];
          check_part = true;
          break;
        default:
          break;
        }
      }
    }
  }
#ifdef AMREX_USE_EB
  if (use_EB) {
    // If particle has moved into a cut-cell
    if (flags(ijk).isSingleValued() && !check_part) {
      // Use the normal and boundary centroid from the new cell
      bloc = {ipn, jpn, kpn};
      normal = {-bnorm(ijk,0), -bnorm(ijk,1), -bnorm(ijk,2)};
      bcentv = {bcent(ijk,0), bcent(ijk,1), bcent(ijk,2)};
      check_part = true;
      // If the particle has moved into a covered cell
    } else if (flags(ijk).isCovered() && !check_part) {
      // and if the cell the particle came from is a cut-cell
      if (flags(ijkp).isSingleValued()) {
        // Use the normal and boundary centroid from the previous cell
        bloc = {ip, jp, kp};
        normal = {-bnorm(ijkp,0), -bnorm(ijkp,1), -bnorm(ijkp,2)};
        bcentv = {bcent(ijkp,0), bcent(ijkp,1), bcent(ijkp,2)};
        // otherwise an error has occurred
      } else {
        Abort("SprayInterpolation::reflect_wall: Particle is outside EB");
      }
      check_part = true;
    }
  }
#endif
  if (check_part) {
    // Projection of vector pointing from EB centroid to particle onto EB normal
    const Real par_dot_EB =
      AMREX_D_TERM((p.pos(0) - (bloc[0] + 0.5 + bcentv[0])*dx[0] - plo[0])*normal[0], +
                   (p.pos(1) - (bloc[1] + 0.5 + bcentv[1])*dx[1] - plo[1])*normal[1], +
                   (p.pos(2) - (bloc[2] + 0.5 + bcentv[2])*dx[2] - plo[2])*normal[2]);
    // Check if particle is on the wrong side of the EB
    if (par_dot_EB < tolerance) {
      // Plane ref point
      Real Nw_Vp = 0.;
      for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
        Nw_Vp += normal[dir]*p.rdata(SPI.pstateVel+dir);
      // Parcel normal velocity
      RealVect Vpn =
        {AMREX_D_DECL(Nw_Vp*normal[0],
                      Nw_Vp*normal[1],
                      Nw_Vp*normal[2])};
      // Parcel tangential velocity
      RealVect Vpt =
        {AMREX_D_DECL(p.rdata(SPI.pstateVel) - Vpn[0],
                      p.rdata(SPI.pstateVel+1) - Vpn[1],
                      p.rdata(SPI.pstateVel+2) - Vpn[2])};
      // dis_fact = 2 -> Particle is reflected off of the wall
      // dis_fact = 1 -> Particle is placed on the wall
      Real dis_fact = 2.;
      // If using a splash model
      if (sigma > 0.) {
        Real mu_part = 0.;
        Real rho_part = 0.;
        // TODO: Determine correct method for handling multi-component liquids
        Real Tstar = 0.; // Average sum Y_i T_wall/T_boil,i
        for (int spf = 0; spf != SPRAY_FUEL_NUM; ++spf) {
          const int ysp = SPI.pstateY + spf;
          mu_part += p.rdata(ysp)*fdat.mu[spf];
          rho_part += p.rdata(ysp)*fdat.rho[spf];
          Tstar += T_wall*p.rdata(ysp)/fdat.boilT[spf];
        }
        const Real dia_part = p.rdata(SPI.pstateDia);
        const Real d3 = pow(dia_part, 3);
        const Real pmass = M_PI/6.*rho_part*d3;
        // Weber number
        const Real We = rho_part*dia_part*Nw_Vp*Nw_Vp/sigma;
        const Real Re_L = std::abs(Nw_Vp)*dia_part*rho_part/mu_part;
        const Real sqReyn = std::sqrt(Re_L);
        const Real Kv = std::sqrt(We)*std::pow(Re_L, 0.25);
        Real umag = 0.;
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
          umag += std::pow(p.rdata(SPI.pstateVel+dir), 2);
        umag = std::sqrt(umag);
        // Inclination of on-coming particle to surface in radians
        const Real alpha = std::asin(std::abs(Nw_Vp)/umag);
        // TODO: Add functionality for wet walls
        bool dry = true;
        splash_type part_crit = splash_criteria(Kv, Tstar, alpha);
        // Boundary layer thickness of particle estimate provided by
        // Pasandideh-Fard et al. 1996
        const Real delta = 2.*dia_part/sqReyn;
        // part_crit == 0, then particle rebounds, no breakup, no deposit
        if (part_crit != splash_type::rebound) {
          dis_fact = 1.; // Place particle at wall
          // Only retain the tangential velocities
          Vpn = RealVect::TheZeroVector();
          // If part_crit == 1, particle is completely deposited onto the wall
          if (part_crit == splash_type::deposit) {
            // Make droplet stationary
            Vpt = RealVect::TheZeroVector();
            // Solve for the max diameter using the volume of a dome shape of height delta
            p.rdata(SPI.pstateDia) = 2.*std::sqrt((6.*pmass/(M_PI*rho_part*delta) - delta*delta)/3.);
            // TODO: Particle should become wall film, not disappear
            p.id() = -1;
            // Otherwise, secondary drops are formed from splashing or thermal breakup
          } else {
            // Determine diameter of secondary droplets
            Real expon = 3.6*(alpha/M_PI)*(alpha/M_PI);
            if (dry)
              SPRF.dia_refl = dia_part*3.3*std::exp(expon)*std::pow(We, -0.65);
            else
              SPRF.dia_refl = dia_part*2.2*std::exp(expon)*std::pow(We, -0.36);
            // Compute angles of reflecting droplets
            compute_angles(alpha, Tstar, We, dry, SPRF);
            // If thermal breakup, all mass is converted to secondary droplets
            Real splash_mass = pmass;
            // If splash, some mass rebounds and some deposits
            if (part_crit == splash_type::splash) {
              // Determine the fraction of the droplet mass
              // that forms secondary drops from Kuhnke 2004
              // TODO: This is incomplete by always assuming a dry wall
              const Real B = 0.2 + 0.6*amrex::Random();
              splash_mass *= amrex::min(1., (Tstar - 0.8)/0.3*(1. - B) + B);
              // Mass deposited into wall film
              Real depot_mass = pmass - splash_mass;
              // Solve for the max diameter using the volume of a dome shape of height delta
              p.rdata(SPI.pstateDia) = 2.*std::sqrt((6.*depot_mass/(M_PI*rho_part*delta) - delta*delta)/3.);
              // TODO: Particle should become wall film if splash occurs
              p.id() = -1;
            } else {
              p.id() = -1;
            }
            Real mass_refl = M_PI/6.*rho_part*std::pow(SPRF.dia_refl, 3);
            Real Ns = splash_mass/mass_refl;
            SPRF.Ns_refl = int(Ns);
            // TODO: Not sure what this variable is but it ranges from 1.45-2
            Real nu32 = 1.45;
            Real We_out = SPRF.dia_refl/dia_part*
              (We*(1. - 0.85*std::pow(std::sin(alpha*M_PI/180.), 2)) + 12.) - 12./nu32;
            SPRF.Unorm = std::sqrt(sigma*We_out/(rho_part*SPRF.dia_refl));
            // This is the dt to move the secondary droplets
            SPRF.dt_pp = par_dot_EB/Nw_Vp;
            // Test vector for finding tangents
            RealVect testvec = {AMREX_D_DECL(-p.rdata(SPI.pstateVel),
                                             -p.rdata(SPI.pstateVel+1),
                                             -p.rdata(SPI.pstateVel+2))};
            // Find two tangent vectors by taking the cross product with the velocity vector
            find_tangents(testvec, normal, SPRF.tanBeta, SPRF.tanPsi);
          }
        } // if (part_crit != splash_type::rebound)
      } // if (sigma > 0)
      // TODO: Make else so that the droplet sticks to the wall
      for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        SPRF.norm[dir] = normal[dir];
        pos[dir] -= dis_fact*par_dot_EB*normal[dir];
        p.pos(dir) -= dis_fact*par_dot_EB*normal[dir];
        p.rdata(SPI.pstateVel+dir) = -Vpn[dir] + Vpt[dir];
      }
    } // if (par_dot_EB < tolerance)
  } // if (check_part)
}

void create_splash_droplet (SprayParticleContainer::ParticleType& p,
                            SprayComps SPI,
                            SprayRefl SPRF)
{
  Real rand1 = amrex::Random();
  Real mean = SPRF.beta_mean;
  Real stdev = SPRF.beta_stdv;
  Real beta = amrex::RandomNormal(mean, stdev);
  beta = std::exp(beta)*M_PI/180.;
  Real omega = SPRF.omega;
  Real expb = SPRF.expomega;
  // If incoming angle is greater than 80 degrees, azimuthal
  // angle of outgoing drop can be between 0 to 2*pi
  Real psi = AMREX_D_PICK(0., 0., rand1*2.*M_PI);
  // Otherwise, follow we modify the recommendation by Naber and Reitz 1988,
  // so the azimuthal angle distribution favors the pre-splash drop path
  // as the inclination angle decreases
  if (omega > 0.) {
    Real rand2 = std::copysign(1., 0.5 - amrex::Random());
    psi = -rand2/omega*std::log(1. - rand1*expb)*M_PI;
  }
  Real costhetad = std::cos(beta);
  Real sinthetad = std::sin(beta);
  Real un = SPRF.Unorm*sinthetad;
  Real utBeta = SPRF.Unorm*costhetad*std::cos(psi);
  Real utPsi = SPRF.Unorm*costhetad*std::sin(psi);
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    Real pvel = AMREX_D_TERM(un*SPRF.norm[dir],
                             +utBeta*SPRF.tanBeta[dir],
                             +utPsi*SPRF.tanPsi[dir]);
    Gpu::Atomic::Add(&p.pos(dir), SPRF.dt_pp*pvel);
    Gpu::Atomic::Add(&p.rdata(SPI.pstateVel+dir), pvel);
  }
}

#endif // _SPRAYWALLS_H_
