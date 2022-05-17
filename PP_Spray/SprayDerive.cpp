
#include "SprayParticles.H"
#include "Drag.H"

using namespace amrex;

void
SprayParticleContainer::computeDerivedVars(
  MultiFab& mf_var,
  const int level,
  const int startComp,
  const Vector<std::string> derivePlotVars,
  std::string* sprayFuelNames)
{
  const int derivePlotVarCount = derivePlotVars.size();
  AMREX_ALWAYS_ASSERT(mf_var.nComp() >= derivePlotVarCount);
  SprayComps SPI = m_sprayIndx;
  const auto dxiarr = this->Geom(level).InvCellSizeArray();
  const auto dxarr = this->Geom(level).CellSizeArray();
  const auto ploarr = this->Geom(level).ProbLoArray();
  const auto phiarr = this->Geom(level).ProbHiArray();
  const RealVect dxi(AMREX_D_DECL(dxiarr[0], dxiarr[1], dxiarr[2]));
  const RealVect dx(AMREX_D_DECL(dxarr[0], dxarr[1], dxarr[2]));
  const RealVect plo(AMREX_D_DECL(ploarr[0], ploarr[1], ploarr[2]));
  const RealVect phi(AMREX_D_DECL(phiarr[0], phiarr[1], phiarr[2]));
  int totalSpecComp = -1;
  int avgSpecComp = -1;
  for (int ivar = 0; ivar < derivePlotVarCount; ++ivar) {
    if (derivePlotVars[ivar] == "spray_total_mass_" + sprayFuelNames[0]) {
      totalSpecComp = ivar;
    } else if (derivePlotVars[ivar] == "spray_avg_mass_" + sprayFuelNames[0]) {
      avgSpecComp = ivar;
    }
  }
  const auto domain = this->Geom(level).Domain();
  mf_var.setVal(0., startComp, derivePlotVarCount);
  for (MyParIter pti(*this, level); pti.isValid(); ++pti) {
    const Box tile_box = pti.tilebox();
    const Long Np = pti.numParticles();
    ParticleType* pstruct = &(pti.GetArrayOfStructs()[0]);
    const SprayData* fdat = d_sprayData;
    Array4<Real> const& vararr = mf_var.array(pti, startComp);
    amrex::ParallelFor(
      Np, [pstruct, SPI, fdat, vararr, plo,
           dxi] AMREX_GPU_DEVICE(int pid) noexcept {
        ParticleType& p = pstruct[pid];
        if (p.id() > 0) {
          RealVect lxc = (p.pos() - plo) * dxi;
          IntVect ijkc = lxc.floor(); // Cell with particle
          Real T_part = p.rdata(SPI.pstateT);
          Real dia_part = p.rdata(SPI.pstateDia);
          Real rho_part = 0.;
          for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
            rho_part += p.rdata(SPI.pstateY + spf) / rhoL(*fdat, T_part, spf);
          }
          rho_part = 1. / rho_part;
          Real pmass = M_PI / 6. * rho_part * std::pow(dia_part, 3);
          Gpu::Atomic::Add(&vararr(ijkc, 0), pmass);
          Gpu::Atomic::Add(&vararr(ijkc, 1), 1.);
        }
      });
    amrex::ParallelFor(
      Np, [pstruct, SPI, fdat, vararr, plo, dxi, totalSpecComp,
           avgSpecComp] AMREX_GPU_DEVICE(int pid) noexcept {
        ParticleType& p = pstruct[pid];
        if (p.id() > 0) {
          RealVect lxc = (p.pos() - plo) * dxi;
          IntVect ijkc = lxc.floor(); // Cell with particle
          Real cell_mass = vararr(ijkc, 0);
          Real cell_num = vararr(ijkc, 1);
          Real T_part = p.rdata(SPI.pstateT);
          Real dia_part = p.rdata(SPI.pstateDia);
          Real rho_part = 0.;
          for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
            rho_part += p.rdata(SPI.pstateY + spf) / rhoL(*fdat, T_part, spf);
          }
          rho_part = 1. / rho_part;
          Real pmass = M_PI / 6. * rho_part * std::pow(dia_part, 3);
          Real avg_temp = pmass * T_part / cell_mass;
          Real avg_mass = pmass / cell_num;
          Real avg_dia = pmass * dia_part / cell_mass;
          Gpu::Atomic::Add(&vararr(ijkc, 2), avg_mass);
          Gpu::Atomic::Add(&vararr(ijkc, 3), avg_dia);
          Gpu::Atomic::Add(&vararr(ijkc, 4), avg_temp);
          if (totalSpecComp >= 0) {
            for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
              Gpu::Atomic::Add(
                &vararr(ijkc, totalSpecComp + spf),
                p.rdata(SPI.pstateY + spf) * pmass);
            }
          }
          if (avgSpecComp >= 0) {
            for (int spf = 0; spf < SPRAY_FUEL_NUM; ++spf) {
              Gpu::Atomic::Add(
                &vararr(ijkc, avgSpecComp + spf),
                p.rdata(SPI.pstateY + spf) * pmass / cell_num);
            }
          }
        }
      });
  }
}
