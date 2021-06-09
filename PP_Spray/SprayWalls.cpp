
#include "SprayParticles.H"
#include "Transport.H"
#include "WallFunctions.H"

using namespace amrex;

// TODO: This file might not be necessary anymore

void
SprayParticleContainer::wallImpingement(
  const int& level,
  const Real& flow_dt,
  const Real& time,
#ifdef AMREX_USE_EB
  const FabArray<EBCellFlagFab>& flagmf,
  const MultiCutFab* bndrycent,
  const MultiCutFab* bndrynorm,
#endif
  const int state_ghosts,
  const int source_ghosts,
  const bool isActive)
{
  BL_PROFILE("ParticleContainer::wallImpingement()");
  const auto dxiarr = this->Geom(level).InvCellSizeArray();
  const auto dxarr = this->Geom(level).CellSizeArray();
  const auto ploarr = this->Geom(level).ProbLoArray();
  const auto phiarr = this->Geom(level).ProbHiArray();
  const RealVect dxi(AMREX_D_DECL(dxiarr[0], dxiarr[1], dxiarr[2]));
  const RealVect dx(AMREX_D_DECL(dxarr[0], dxarr[1], dxarr[2]));
  const RealVect plo(AMREX_D_DECL(ploarr[0], ploarr[1], ploarr[2]));
  const RealVect phi(AMREX_D_DECL(phiarr[0], phiarr[1], phiarr[2]));
  const auto domain = this->Geom(level).Domain();
  IntVect bndry_lo; // Designation for boundary types
  IntVect bndry_hi; // 0 - Periodic, 1 - Reflective, -1 - Non-reflective
  for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
    if (!this->Geom(level).isPeriodic(dir)) {
      if (reflect_lo[dir])
        bndry_lo[dir] = 1;
      else
        bndry_lo[dir] = -1;
      if (reflect_hi[dir])
        bndry_hi[dir] = 1;
      else
        bndry_hi[dir] = -1;
    } else {
      bndry_lo[dir] = 0;
      bndry_hi[dir] = 0;
    }
  }
  SprayComps SPI = m_sprayIndx;
  SprayUnits SPU;
  // Loop back over particles to see if any have interacted with walls
  // This should occur on the host
  for (MyParIter pti(*this, level); pti.isValid(); ++pti) {
    PairIndex index(pti.index(), pti.LocalTileIndex());
    const Box tile_box = pti.tilebox();
    const Box src_box = pti.growntilebox(source_ghosts);
    const Long Np = pti.numParticles();
    // Check if tile has walls
    bool at_bounds = tile_at_bndry(tile_box, bndry_lo, bndry_hi, domain);
#ifdef AMREX_USE_EB
    const EBCellFlagFab& flags = flagmf[pti];
    bool eb_in_box = false;
    if (flags.getType(src_box) != FabType::regular) {
      eb_in_box = true;
      at_bounds = true;
    }
#endif
    // Check if particles have gone outside of walls
    if (Np > 0 && at_bounds) {
      Vector<IntVect> film_locs;
      FArrayBox wall_film(src_box, SPI.wf_num);
      wall_film.setVal<RunOn::Host>(0., src_box, 0, SPI.wf_num);
      // This IAB notes the particle ID that will represent
      // the wall film at each location
      IArrayBox film_id(src_box, 1);
      film_id.setVal<RunOn::Host>(-1, src_box, 0, 1);
      const SprayData* fdat = m_sprayData;
      auto& ptile = GetParticles(level)[index];
      auto& pval = ptile.GetArrayOfStructs();
#ifdef AMREX_USE_EB
      Array4<const EBCellFlag> flags_fab = flags.array();
      Array4<const Real> bcent_fab;
      Array4<const Real> bnorm_fab;
      if (eb_in_box) {
        bcent_fab = bndrycent->array(pti);
        bnorm_fab = bndrynorm->array(pti);
      }
#endif
      for (int pid = 0; pid < Np; ++pid) {
        ParticleType& p = pval[pid];
        if (p.id() > 0) {
          const RealVect lx = (p.pos() - plo) * dxi;
          IntVect ijk = lx.floor(); // Closest cell center
          const Real T_part = p.rdata(SPI.pstateT);
          const Real dia_part = p.rdata(SPI.pstateDia);
          splash_type splash_flag = splash_type::no_impact;
          // Check if particle is already wall film
          if (T_part > 0.) {
            SprayRefl SPRF; // Structure holding data for reflected particles
            SPRF.pos_refl = p.pos();
            for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf)
              SPRF.Y_refl[spf] = p.rdata(SPI.pstateY + spf);
            bool dry_wall = true;
            if (film_id(ijk, 0) > 0)
              dry_wall = false;
//             splash_flag = impose_wall(
//               p, SPI, SPU, *fdat, ijk, dx, dxi, plo, phi,
// #ifdef AMREX_USE_EB
//               eb_in_box, flags_fab, bcent_fab, bnorm_fab,
// #endif
//               bndry_lo, bndry_hi, flow_dt, m_wallT, SPRF, isActive, dry_wall);
            // Only add active particles, not ghost or virtual
            if (SPRF.Ns_refl > 0 && isActive) {
              for (int nsp = 0; nsp < SPRF.Ns_refl; ++nsp) {
                ParticleType pnew;
                pnew.id() = ParticleType::NextID();
                pnew.cpu() = ParallelDescriptor::MyProc();
                pnew.rdata(SPI.pstateDia) = SPRF.dia_refl;
                pnew.rdata(SPI.pstateT) = T_part;
                for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf)
                  pnew.rdata(SPI.pstateY + spf) = SPRF.Y_refl[spf];
                for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
                  pnew.rdata(SPI.pstateVel + dir) = 0.;
                  pnew.pos(dir) = SPRF.pos_refl[dir];
                }
                create_splash_droplet(pnew, SPI, SPRF, SPU);
                ptile.push_back(pnew);
              }
            } // if (Ns_refl > 0)
          } else {
            splash_flag = splash_type::wall_film;
            p.rdata(SPI.pstateT) *= -1.;
          }
          // Check if droplet is deposited, splashes, or is already a wall film
          if (
            splash_flag == splash_type::deposit ||
            splash_flag == splash_type::splash ||
            splash_flag == splash_type::wall_film) {
            if (film_id(ijk, 0) < 0) {
              film_id(ijk, 0) = pid;
              film_locs.push_back(ijk);
            } else {
              p.id() = -1;
            }
            // Velocity component now holds the volume
            Real new_vol = p.rdata(SPI.pstateVol);
            wall_film(ijk, SPI.wf_vol) += new_vol;
            wall_film(ijk, SPI.wf_temp) += new_vol * p.rdata(SPI.pstateT);
            Real drop_height = p.rdata(SPI.pstateHt);
            wall_film(ijk, SPI.wf_ht) += drop_height;
            for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf)
              wall_film(ijk, SPI.wf_Y + spf) +=
                new_vol * p.rdata(SPI.pstateY + spf);
          }
        } // if (p.id() > 0)
      }   // for (int pid...
      for (int wfl = 0; wfl < film_locs.size(); ++wfl) {
        IntVect ijk = film_locs[wfl];
        int pid = film_id(ijk, 0);
        Real vol = wall_film(ijk, SPI.wf_vol);
        Real T = wall_film(ijk, SPI.wf_temp) / vol;
        ParticleType& p = pval[pid];
        for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf)
          p.rdata(SPI.pstateY + spf) = wall_film(ijk, SPI.wf_Y + spf) / vol;
        // Diameter index will hold height
        p.rdata(SPI.pstateHt) = wall_film(ijk, SPI.wf_ht);
        // Velocity index will hold volume
        p.rdata(SPI.pstateVol) = vol;
        // For wall films, the temperature is set to a negative value
        // TODO: Determine better way to model the wall film temperature
        p.rdata(SPI.pstateT) = -T;
      }
    } // if (do_move && Np > 0 && at_bounds)
  }   // for (MyParIter pti ...
}
