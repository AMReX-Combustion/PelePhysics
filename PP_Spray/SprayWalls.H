#ifndef _SPRAYWALLS_H_
#define _SPRAYWALLS_H_

AMREX_GPU_HOST_DEVICE AMREX_INLINE
void impose_wall(SprayParticleContainer::ParticleType& p,
                 SprayComps SPI,
                 SprayUnits SPU,
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
                 const amrex::IntVect& bndry_hi)
{
    const amrex::Real tolerance = std::numeric_limits<amrex::Real>::epsilon();
    // Determine the location of the particle after moving
    const amrex::RealVect lx = (p.pos() - plo)*dxi;
    const amrex::IntVect ijk = lx.floor();
    AMREX_D_TERM(const int ipn = ijk[0];, const int jpn = ijk[1];, const int kpn = ijk[2];);
    AMREX_D_TERM(const int ip = ijkp[0];, const int jp = ijkp[1];, const int kp = ijkp[2];);
    amrex::RealVect normal(RealVect::TheZeroVector()); // Normal vector to boundary face
    amrex::RealVect bcentv(RealVect::TheZeroVector()); // Point on boundary face
    amrex::IntVect bloc(ijk);    // Index of nearest cell center
    // If reflective
    bool check_part = false;
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
        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
          Nw_Vp += normal[dir]*p.rdata(SPI.pstateVel+dir);
          p.pos(dir) -= 2.*par_dot_EB*normal[dir];
        }

        // Parcel normal velocity
        const amrex::RealVect Vpn =
          {AMREX_D_DECL(Nw_Vp*normal[0],
                        Nw_Vp*normal[1],
                        Nw_Vp*normal[2])};
        // Parcel tangential velocity
        const amrex::RealVect Vpt =
          {AMREX_D_DECL(p.rdata(SPI.pstateVel) - Vpn[0],
                        p.rdata(SPI.pstateVel+1) - Vpn[1],
                        p.rdata(SPI.pstateVel+2) - Vpn[2])};
        // Rebound particle if it moving toward the wall
        if (Nw_Vp < 0.) {
          for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
            p.rdata(SPI.pstateVel+dir) = -Vpn[dir] + Vpt[dir];
        }
      }
    }
}

#endif // _SPRAYWALLS_H_
