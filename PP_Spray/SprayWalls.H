#ifndef _SPRAYWALLS_H_
#define _SPRAYWALLS_H_

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void impose_wall (SprayParticleContainer::ParticleType& p,
                  SprayComps SPI,
                  SprayUnits SPU,
                  GasPhaseVals gpv,
                  SprayData fdat,
                  const amrex::IntVect& ijkp,
                  const amrex::RealVect& dx,
                  const amrex::RealVect& dxi,
#ifdef AMREX_USE_EB
                  amrex::Array4<const EBCellFlag> const& flags,
                  amrex::Array4<const amrex::Real> const& ccent,
                  amrex::Array4<const amrex::Real> const& bcent,
                  amrex::Array4<const amrex::Real> const& bnorm,
#endif
                  const amrex::RealVect& plo,
                  const amrex::RealVect& phi,
                  const amrex::IntVect& bndry_lo,
                  const amrex::IntVect& bndry_hi,
                  int& Ns_pp,
                  amrex::Real& dt_pp,
                  amrex::Real& beta_max,
                  amrex::RealVect& norm_pp)
{
    const amrex::Real tolerance = std::numeric_limits<amrex::Real>::epsilon();
    // Splash threshold
    const amrex::Real K_splash = 57.7;
    amrex::Real sigma = fdat.sigma;
    amrex::Real cangle = fdat.cangle;
    // Determine the location of the particle after moving
    const amrex::RealVect lx = (p.pos() - plo)*dxi;
    const amrex::IntVect ijk = lx.floor();
    AMREX_D_TERM(const int ipn = ijk[0];, const int jpn = ijk[1];, const int kpn = ijk[2];);
    AMREX_D_TERM(const int ip = ijkp[0];, const int jp = ijkp[1];, const int kp = ijkp[2];);
    amrex::RealVect normal(RealVect::TheZeroVector()); // Normal vector to boundary face
    amrex::RealVect bcentv(RealVect::TheZeroVector()); // Point on boundary face
    amrex::IntVect bloc(ijk); // Index of nearest cell center
    bool check_part = false; // Check if particle is reflected
    // First check if particle has exited the domain through a Cartesian boundary
    for (int hilo = 0; hilo < 2; ++hilo) {
      amrex::Real fact = 1.;
      if (hilo == 1) fact = -1.;
      for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
        amrex::Real dom_loc = plo[dir];
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
    // If particle has moved into a cut-cell
    if (flags(ipn,jpn,kpn).isSingleValued() && !check_part) {
      // Use the normal and boundary centroid from the new cell
      bloc = {ipn, jpn, kpn};
      normal = {-bnorm(ipn,jpn,kpn,0), -bnorm(ipn,jpn,kpn,1), -bnorm(ipn,jpn,kpn,2)};
      bcentv = {bcent(ipn,jpn,kpn,0), bcent(ipn,jpn,kpn,1), bcent(ipn,jpn,kpn,2)};
      check_part = true;
      // If the particle has moved into a covered cell
    } else if (flags(ipn,jpn,kpn).isCovered() && !check_part) {
      // and if the cell the particle came from is a cut-cell
      if (flags(ip,jp,kp).isSingleValued()) {
        // Use the normal and boundary centroid from the previous cell
        bloc = {ip, jp, kp};
        normal = {-bnorm(ip,jp,kp,0), -bnorm(ip,jp,kp,1), -bnorm(ip,jp,kp,2)};
        bcentv = {bcent(ip,jp,kp,0), bcent(ip,jp,kp,1), bcent(ip,jp,kp,2)};
        // otherwise an error has occurred
      } else {
        amrex::Abort("SprayInterpolation::reflect_wall: Particle is outside EB");
      }
      check_part = true;
    }
#endif
    if (check_part) {
      // Projection of vector pointing from EB centroid to particle onto EB normal
      const amrex::Real par_dot_EB =
        AMREX_D_TERM((p.pos(0) - (bloc[0] + 0.5 + bcentv[0])*dx[0] - plo[0])*normal[0], +
                     (p.pos(1) - (bloc[1] + 0.5 + bcentv[1])*dx[1] - plo[1])*normal[1], +
                     (p.pos(2) - (bloc[2] + 0.5 + bcentv[2])*dx[2] - plo[2])*normal[2]);
      // Check if particle is on the wrong side of the EB
      if (par_dot_EB < tolerance) {
        // Plane ref point
        amrex::Real Nw_Vp = 0.;
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
          Nw_Vp += normal[dir]*p.rdata(SPI.pstateVel+dir);

        // Parcel normal velocity
        amrex::RealVect Vpn =
          {AMREX_D_DECL(Nw_Vp*normal[0],
                        Nw_Vp*normal[1],
                        Nw_Vp*normal[2])};
        // Parcel tangential velocity
        amrex::RealVect Vpt =
          {AMREX_D_DECL(p.rdata(SPI.pstateVel) - Vpn[0],
                        p.rdata(SPI.pstateVel+1) - Vpn[1],
                        p.rdata(SPI.pstateVel+2) - Vpn[2])};
        const amrex::Real rho_fluid = gpv.rho_fluid;
        const amrex::Real dia_part = p.rdata(SPI.pstateDia);
        amrex::Real rho_part = 0.;
        for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf)
          rho_part += p.rdata(SPI.pstateY+spf)*fdat.rho(spf);
        const amrex::Real d3 = pow(dia_part, 3);
        // Weber number
        amrex::Real mu_part = 0.;
        for (int spf = 0; spf != SPRAY_FUEL_NUM; ++spf)
          mu_part += p.rdata(SPI.pstateY+spf)*fdat.mu(spf);
        const amrex::Real We = rho_part*dia_part*Nw_Vp*Nw_Vp/sigma;
        const amrex::Real Re_L = std::abs(Nw_Vp)*dia_part*rho_part/mu_part;
        const amrex::Real sqReyn = sqrt(Re_L);
        const amrex::Real Kv = sqrt(We)*pow(Re_L, 0.25);
        // dis_fact = 2 -> Particle is reflected off of the wall
        // dis_fact = 1 -> Particle is placed on the wall
        amrex::Real dis_fact = 2.;
        if (Kv > K_splash && sigma > 0.) {
          dis_fact = 1.; // Place particle at wall
          beta_max =
            sqrt((12. + We)/(3.*(1. - std::cos(cangle)) + 4.*We/sqReyn));
          // Number of secondary drops formed
          amrex::Real Ns =
            std::abs(Nw_Vp)/(20.*std::sqrt(mu_part/rho_part))*std::pow(M_PI*M_PI*rho_fluid*d3/sigma, 0.25);
          Ns_pp = int(Ns);
          // This is the dt to move the secondary droplets
          dt_pp = par_dot_EB/Nw_Vp;
          // Only retain the tangential velocities
          Vpn = RealVect::TheZeroVector();
        }
        // TODO: Make else so that the droplet sticks to the wall
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
          norm_pp[dir] = normal[dir];
          p.pos(dir) -= dis_fact*par_dot_EB*normal[dir];
          p.rdata(SPI.pstateVel+dir) = -Vpn[dir] + Vpt[dir];
        }
      }
    }
}

void find_tangents (const amrex::RealVect& testvec,
                    const amrex::RealVect& norm,
                    amrex::RealVect& tan1,
                    amrex::RealVect& tan2)
{
#if AMREX_SPACEDIM == 3
  tan1 = testvec.crossProduct(norm);
  tan2 = tan1.crossProduct(norm);
  tan1 /= tan1.vectorLength();
  tan2 /= tan2.vectorLength();
#else
  tan1[0] = -norm[1];
  tan1[1] = norm[0];
#endif
}

void splash_vel_dia (const int Ns,
                     const amrex::Real sigma,
                     const amrex::Real cangle,
                     const amrex::Real dia_part,
                     const amrex::Real rho_part,
                     const amrex::Real beta_max,
                     amrex::Real& vel_pp,
                     amrex::Real& dia_pp)
{
  amrex::Real d_max = beta_max*dia_part;
  amrex::Real pmass = M_PI/6.*rho_part*dia_part*dia_part*dia_part;
  // Mass and diameter of secondary drops
  amrex::Real mass_pp = pmass/Ns;
  dia_pp = std::pow(mass_pp*6./(M_PI*rho_part), 1./3.);
  amrex::Real Es_p =
    M_PI/4.*d_max*d_max*(1. - std::cos(cangle))*sigma;
  amrex::Real Es_pp = sigma*M_PI*dia_pp*dia_pp;
  amrex::Real Eks_pp = Es_p/Ns - Es_pp;
  if (Eks_pp < 0.)
    amrex::Abort("splash_vel_dia contact angle resulted in negative droplet energy");
  vel_pp = std::sqrt(2.*Eks_pp/mass_pp);
}

void create_splash_droplet (SprayParticleContainer::ParticleType& p,
                            SprayComps SPI,
                            const amrex::Real Unorm,
                            const amrex::Real dt_pp,
                            const amrex::RealVect p_pos,
                            const amrex::RealVect normal,
                            const amrex::RealVect tangent1,
                            const amrex::RealVect tangent2)
{
  amrex::Real rand1 = amrex::Random();
  // This ensures the tangential velocity is {1/32,3/8}
  amrex::Real costhetad = 1./32. + 11./32.*rand1;
  amrex::Real sinthetad = std::sqrt(1. - costhetad*costhetad);
  amrex::Real psi = AMREX_D_PICK(0.,0.,amrex::Random()*2.*M_PI);
  amrex::Real un = Unorm*sinthetad;
  amrex::Real ut1 = Unorm*costhetad*std::cos(psi);
  amrex::Real ut2 = Unorm*costhetad*std::sin(psi);
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    amrex::Real pvel = AMREX_D_TERM(un*normal[dir],+ut1*tangent1[dir],+ut2*tangent2[dir]);
    p.pos(dir) = p_pos[dir] + dt_pp*pvel;
    p.rdata(SPI.pstateVel+dir) += pvel;
  }
}

#endif // _SPRAYWALLS_H_
