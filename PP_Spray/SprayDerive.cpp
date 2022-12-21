
#include "SprayParticles.H"
#ifdef AMREX_USE_EB
#include <AMReX_EBFArrayBox.H>
#endif

using namespace amrex;

void
SprayParticleContainer::computeDerivedVars(
  MultiFab& mf_var, const int level, const int start_indx)
{
  auto derivePlotVarCount = spray_derive_vars.size();
  AMREX_ALWAYS_ASSERT(mf_var.nComp() >= derivePlotVarCount);
  const auto dxiarr = this->Geom(level).InvCellSizeArray();
  const auto dxarr = this->Geom(level).CellSizeArray();
  const auto ploarr = this->Geom(level).ProbLoArray();
  const auto phiarr = this->Geom(level).ProbHiArray();
  const RealVect dxi(AMREX_D_DECL(dxiarr[0], dxiarr[1], dxiarr[2]));
  const RealVect dx(AMREX_D_DECL(dxarr[0], dxarr[1], dxarr[2]));
  const RealVect plo(AMREX_D_DECL(ploarr[0], ploarr[1], ploarr[2]));
  const RealVect phi(AMREX_D_DECL(phiarr[0], phiarr[1], phiarr[2]));
  const Real cell_vol = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);
#ifdef AMREX_USE_EB
  const auto& factory =
    dynamic_cast<EBFArrayBoxFactory const&>(mf_var.Factory());
  const auto* volfrac = &(factory.getVolFrac());
#endif
  int total_spec_indx = -1;
  for (auto ivar = 0; ivar < derivePlotVarCount; ++ivar) {
    if (spray_derive_vars[ivar] == "spray_mass_" + spray_fuel_names[0]) {
      total_spec_indx = ivar;
    }
  }
  const int mass_indx = 0;
  const int num_indx = mass_indx + 1;
  const int vol_indx = num_indx + 1;
  const int surf_indx = vol_indx + 1;
  const int volf_indx = surf_indx + 1;
  const int d10_indx = volf_indx + 1;
  const int d32_indx = d10_indx + 1;
  const int wfh_indx = d32_indx + 1;
  const int temp_indx = wfm_indx + 1;
  const int vel_indx = temp_indx + 1;
  for (MyParIter pti(*this, level); pti.isValid(); ++pti) {
    const Long Np = pti.numParticles();
    ParticleType* pstruct = pti.GetArrayOfStructs().data();
    const SprayData* fdat = d_sprayData;
    FArrayBox& varfab = mf_var[pti];
    Array4<Real> const& vararr = mf_var.array(pti, start_indx);
#ifdef AMREX_USE_EB
    Box box = pti.tilebox();
    const auto& interp_fab = static_cast<EBFArrayBox const&>(mf_var[pti]);
    const EBCellFlagFab& flags = interp_fab.getEBCellFlagFab();
    Array4<const Real> volfrac_fab;
    const auto& flags_array = flags.array();
    if (flags.getType(box) != FabType::regular) {
      volfrac_fab = volfrac->array(pti);
    }
#endif
    amrex::ParallelFor(
      Np, [pstruct, fdat, vararr, plo, dxi, cell_vol, total_spec_indx
#ifdef AMREX_USE_EB
           ,
           flags_array, volfrac_fab
#endif
    ] AMREX_GPU_DEVICE(int pid) noexcept {
        ParticleType& p = pstruct[pid];
        if (p.id() > 0) {
          RealVect lxc = (p.pos() - plo) * dxi;
          IntVect ijkc = lxc.floor(); // Cell with particle
          Real T_part = p.rdata(SprayComps::pstateT);
          Real dia_part = p.rdata(SprayComps::pstateDia);
          Real rho_part = 0.;
          for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
            rho_part +=
              p.rdata(SprayComps::pstateY + spf) / fdat->rhoL(T_part, spf);
          }
          rho_part = 1. / rho_part;
          Real surf = M_PI * dia_part * dia_part;
          Real vol = M_PI / 6. * std::pow(dia_part, 3);
          Real pmass = vol * rho_part;
          Real num_ppp = fdat->num_ppp; // Particles per parcel
          // TODO: Adjust face area for EB
          Real film_hght = p.rdata(SprayComps::pstateFilmVol) /
                           (AMREX_D_TERM(1., *dx[0], *dx[0]));
          if (film_hght == 0.) {
            Gpu::Atomic::Add(&vararr(ijkc, mass_indx), num_ppp * pmass);
            Gpu::Atomic::Add(&vararr(ijkc, num_indx), num_ppp);
            Gpu::Atomic::Add(&vararr(ijkc, vol_indx), num_ppp * vol);
            Gpu::Atomic::Add(&vararr(ijkc, surf_indx), num_ppp * surf);
            amrex::Real curvol = cell_vol;
#ifdef AMREX_USE_EB
            if (!flags_array(ijkc).isRegular()) {
              curvol /= volfrac_fab(ijkc);
            }
#endif
            Gpu::Atomic::Add(&vararr(ijkc, volf_indx), num_ppp * vol / curvol);
            Gpu::Atomic::Add(
              &vararr(ijkc, d10_indx),
              num_ppp * dia_part); // To be divided by num later
            Gpu::Atomic::Add(
              &vararr(ijkc, d32_indx),
              num_ppp * vol * 6.); // To be divided by surf later
            Gpu::Atomic::Add(
              &vararr(ijkc, temp_indx), num_ppp * pmass * T_part);
            for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
              Gpu::Atomic::Add(
                &vararr(ijkc, vel_indx + dir),
                num_ppp * pmass * p.rdata(SprayComps::pstateVel + dir));
            }
            if (total_spec_indx >= 0) {
              for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
                Gpu::Atomic::Add(
                  &vararr(ijkc, total_spec_indx + spf),
                  p.rdata(SprayComps::pstateY + spf) * pmass);
              }
            }
          } else {
            Gpu::Atomic::Add(&vararr(ijkc, wfh_indx), film_hght);
          }
        }
      });
#ifdef AMREX_USE_GPU
    varfab.protected_divide<RunOn::Device>(
      varfab, start_indx + num_indx, start_indx + d10_indx); // Get d10
    varfab.protected_divide<RunOn::Device>(
      varfab, start_indx + surf_indx, start_indx + d32_indx); // Get d32
    varfab.protected_divide<RunOn::Device>(
      varfab, start_indx + mass_indx, start_indx + temp_indx);
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      varfab.protected_divide<RunOn::Device>(
        varfab, start_indx + mass_indx,
        start_indx + vel_indx + dir); // Divide momentum by total mass
    }
#else
    varfab.protected_divide(
      varfab, start_indx + num_indx, start_indx + d10_indx); // Get d10
    varfab.protected_divide(
      varfab, start_indx + surf_indx, start_indx + d32_indx); // Get d32
    varfab.protected_divide(
      varfab, start_indx + mass_indx, start_indx + temp_indx);
    for (int dir = 0; dir < AMREX_SPACEDIM; ++dir) {
      varfab.protected_divide(
        varfab, start_indx + mass_indx,
        start_indx + vel_indx + dir); // Divide momentum by total mass
    }
#endif
  }
}
