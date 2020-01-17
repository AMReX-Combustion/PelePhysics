
#include <SprayParticles.H>
#include <AMReX_Particles.H>
#include <Transport_F.H>
#include <drag_F.H>
#include <Drag.H>

using namespace amrex;

void
SprayParticleContainer::init_bcs()
{
  for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
      if (phys_bc->lo(dir) == Symmetry   ||
	  phys_bc->lo(dir) == SlipWall   ||
	  phys_bc->lo(dir) == NoSlipWall)
	{
          reflect_lo[dir] = 1;
	}
      else
	{
	  reflect_lo[dir] = 0;
	}
      if (phys_bc->hi(dir) == Symmetry   ||
	  phys_bc->hi(dir) == SlipWall   ||
	  phys_bc->hi(dir) == NoSlipWall)
	{
          reflect_hi[dir] = 1;
	}
      else
	{
          reflect_hi[dir] = 0;
	}
    }
}

void
SprayParticleContainer::SetAll (Real val, int pstate_idx, int lev)
{
  BL_ASSERT(lev >= 0 && lev < GetParticles().size());
  ParticleLevel& plev = GetParticles(lev);
  for (auto& kv : plev)
    {
      AoS& particles = kv.second.GetArrayOfStructs();
      for (auto& p : particles)
	{
	  if (p.id() > 0)
	    {
	      p.rdata(pstate_idx) = val;
	    }
	}
    }
  return;
}

void
SprayParticleContainer::moveKickDrift (MultiFab& state,
				       MultiFab& source,
				       int       lev,
				       Real      dt,
				       int       tmp_src_width,
				       int       where_width)
{
  const Real* dx = Geom(lev).CellSize();
  const RealVect dx_vect(AMREX_D_DECL(dx[0], dx[1], dx[2]));
  ParticleLevel& pmap = this->GetParticles(lev);
  BL_PROFILE("ParticleContainer::moveKickDrift()");
  BL_ASSERT(lev >= 0);
  BL_ASSERT(state.nGrow() >= 2);

  //If there are no particles at this level
  if (lev >= this->GetParticles().size())
    return;

  const Real strttime = ParallelDescriptor::second();

  MultiFab* state_ptr;
  MultiFab* tmp_src_ptr;

  // ********************************************************************************
  // We only make a new state_ptr if the boxArray of the state differs from the
  // boxArray onto which the particles are decomposed
  // ********************************************************************************

  if (this->OnSameGrids(lev, state))
    {
      state_ptr = &state;
    }
  else
    {
      state_ptr = new MultiFab(this->m_gdb->ParticleBoxArray(lev),
			       this->m_gdb->ParticleDistributionMap(lev),
			       state.nComp(), state.nGrow());
      state_ptr->setVal(0.);
      state_ptr->copy(state,0,0,state.nComp());
      state_ptr->FillBoundary(Geom(lev).periodicity());
    }

  // ********************************************************************************
  // We make a temporary MultiFab for the source here and initialize it to zero
  // because if we use one that already has values in it, the SumBoundary call
  // will end up duplicating those values with each call
  // ********************************************************************************
  tmp_src_width = 3; // HACK?

  tmp_src_ptr = new MultiFab(this->m_gdb->ParticleBoxArray(lev),
			     this->m_gdb->ParticleDistributionMap(lev),
			     source.nComp(), tmp_src_width);
  tmp_src_ptr->setVal(0.);

  const Real* plo = Geom(lev).ProbLo();
  const Real* phi = Geom(lev).ProbHi();
  const Box& domain = Geom(lev).Domain();
  int do_move = 1;

  for (MyParIter pti(*this, lev); pti.isValid(); ++pti)
    {
      AoS& particles = pti.GetArrayOfStructs();
      int Np = particles.size();
      if (Np > 0)
	{
	  const Box& state_box = (*state_ptr)[pti].box();
	  const Box&   src_box = (*tmp_src_ptr)[pti].box();
	  update_particlesCPP(lev, particles, state[pti], (*tmp_src_ptr)[pti],
	  		      domain, Geom(lev).ProbDomain(), reflect_lo,
	  		      reflect_hi, dx_vect, dt, do_move);
// 	  update_particles(&Np,&lev, particles.data(),
// 	                   (*state_ptr)[pti].dataPtr(),
// 	                   state_box.loVect(), state_box.hiVect(),
// 	                   (*tmp_src_ptr)[pti].dataPtr(),
// 	                   src_box.loVect(), src_box.hiVect(),
// 	                   domain.loVect(), domain.hiVect(),
// 	                   plo,phi,reflect_lo.getVect(),reflect_hi.getVect(),dx,dt,&do_move);
	}
    }

  // ********************************************************************************
  // Make sure the momentum put into ghost cells of each grid is added to both 
  //      valid regions AND the ghost cells of other grids.  If at level = 0 
  //      we can accomplish this with SumBoundary; however if this level has ghost
  //      cells not covered at this level then we need to do more.
  // ********************************************************************************
  if (lev > 0)
    {
      int ncomp = tmp_src_ptr->nComp();
      MultiFab tmp(this->m_gdb->ParticleBoxArray(lev),
		   this->m_gdb->ParticleDistributionMap(lev),
		   ncomp, tmp_src_ptr->nGrow());
      tmp.setVal(0.);
      tmp.copy(*tmp_src_ptr, 0, 0, ncomp, tmp_src_ptr->nGrow(), tmp_src_ptr->nGrow(),
	       Geom(lev).periodicity(), FabArrayBase::ADD);
      tmp_src_ptr->copy(tmp, 0, 0, ncomp, tmp_src_width, tmp_src_width,
			Geom(lev).periodicity(), FabArrayBase::COPY);
    }
  else
    {
      tmp_src_ptr->SumBoundary(Geom(lev).periodicity());
    }

  // Add new sources into source *after* we have called SumBoundary
  MultiFab::Add(source,*tmp_src_ptr,0,0,source.nComp(),
		std::min(source.nGrow(),tmp_src_ptr->nGrow()));
  delete tmp_src_ptr;

  // Fill ghost cells after we've synced up ..
  source.FillBoundary(Geom(lev).periodicity());

  // Only delete this if in fact we created it.  Note we didn't change state_ptr so 
  //  we don't need to copy anything back into state
  if (state_ptr != &state ) delete state_ptr;

  // ********************************************************************************

  if (lev > 0 && sub_cycle)
    {
      ParticleLocData pld;
      for (auto& kv : pmap)
	{
	  AoS& pbox = kv.second.GetArrayOfStructs();
	  const int n = pbox.size();
#ifdef _OPENMP
#pragma omp parallel for private(pld)
#endif
	  for (int i = 0; i < n; i++)
	    {
	      ParticleType& p = pbox[i];
	      if (p.id() <= 0) continue;
	      // Move the particle to the proper ghost cell.
	      //      and remove any *ghost* particles that have gone too far
	      // Note that this should only negate ghost particles, not real particles.
	      if (!this->Where(p, pld, lev, lev, where_width))
		{
		  // Assert that the particle being removed is a ghost particle;
		  // the ghost particle is no longer in relevant ghost cells for this grid.
		  if (p.id() == GhostParticleID)
		    {
		      p.id() = -1;
		    }
		  else
		    {
		      std::cout << "Oops -- removing particle " << p.id() << " " << p.pos(0) << " " << p.pos(1) << std::endl;
		      p.id() = -1;
		    }
		}
	    }
	}
    }

  // ********************************************************************************

  if (this->m_verbose > 1)
    {
      Real stoptime = ParallelDescriptor::second() - strttime;
      ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());
      if (ParallelDescriptor::IOProcessor())
	{
	  Print() << "SprayParticleContainer::moveKickDrift() time: "
		  << stoptime << '\n';
	}
    }
}

void
SprayParticleContainer::moveKick(MultiFab& state,
				 MultiFab& source,
				 int       lev,
				 Real      dt,
				 int       tmp_src_width)
{
  const Real* dx = Geom(lev).CellSize();
  const RealVect dx_vect(AMREX_D_DECL(dx[0], dx[1], dx[2]));
  ParticleLevel& pmap = this->GetParticles(lev);

  BL_PROFILE("ParticleContainer::moveKick()");
  BL_ASSERT(lev >= 0);
  BL_ASSERT(state.nGrow() >= 2);

  //If there are no particles at this level
  if (lev >= this->GetParticles().size())
    return;

  const Real strttime = ParallelDescriptor::second();

  MultiFab* state_ptr;
  MultiFab* tmp_src_ptr;

  // ********************************************************************************
  // We only make a new state_ptr if the boxArray of the state differs from the
  // boxArray onto which the particles are decomposed
  // ********************************************************************************

  if (this->OnSameGrids(lev, state))
    {
      state_ptr = &state;
    }
  else
    {
      state_ptr = new MultiFab(this->m_gdb->ParticleBoxArray(lev),
			       this->m_gdb->ParticleDistributionMap(lev),
			       state.nComp(),state.nGrow());
      state_ptr->setVal(0.);
      state_ptr->copy(state,0,0,state.nComp());
      state_ptr->FillBoundary(Geom(lev).periodicity());
    }

  // ********************************************************************************
  // We make a temporary MultiFab for the source here and initialize it to zero
  // because if we use one that already has values in it, the SumBoundary call
  // will end up duplicating those values with each call
  // ********************************************************************************

  tmp_src_width = 3; // HACK?

  tmp_src_ptr = new MultiFab(this->m_gdb->ParticleBoxArray(lev),
			     this->m_gdb->ParticleDistributionMap(lev),
			     source.nComp(),tmp_src_width);
  tmp_src_ptr->setVal(0.);

  // ********************************************************************************

  const Real* plo = Geom(lev).ProbLo();
  const Real* phi = Geom(lev).ProbHi();

  const Box& domain = Geom(lev).Domain();

  int do_move = 0;

  for (MyParIter pti(*this, lev); pti.isValid(); ++pti)
    {
      AoS& particles = pti.GetArrayOfStructs();
      int Np = particles.size();
      if (Np > 0)
	{
	  const Box& state_box = (*state_ptr)[pti].box();
	  const Box&   src_box = (*tmp_src_ptr)[pti].box();
	  update_particlesCPP(lev, particles, state[pti], (*tmp_src_ptr)[pti],
	  		      domain, Geom(lev).ProbDomain(), reflect_lo,
	  		      reflect_hi, dx_vect, dt, do_move);
//  	  update_particles(&Np,&lev, particles.data(),
//  	                   (*state_ptr)[pti].dataPtr(),
//  	                   state_box.loVect(), state_box.hiVect(),
//  	                   (*tmp_src_ptr)[pti].dataPtr(),
//  	                   src_box.loVect(), src_box.hiVect(),
//  	                   domain.loVect(), domain.hiVect(),
//  	                   plo,phi,reflect_lo.getVect(),reflect_hi.getVect(),
//  	                   dx,dt,&do_move);
	}
    }

  // ********************************************************************************
  // Make sure the momentum put into ghost cells of each grid is added to both
  //      valid regions AND the ghost cells of other grids.  If at level = 0
  //      we can accomplish this with SumBoundary; however if this level has ghost
  //      cells not covered at this level then we need to do more.
  // ********************************************************************************
  if (lev > 0)
    {
      int ncomp = tmp_src_ptr->nComp();
      int ngrow = source.nGrow();
      MultiFab tmp(this->m_gdb->ParticleBoxArray(lev),
		   this->m_gdb->ParticleDistributionMap(lev),
		   ncomp, ngrow);
      tmp.setVal(0.);
      tmp.copy(*tmp_src_ptr, 0, 0, ncomp, ngrow, ngrow, Geom(lev).periodicity(),
	       FabArrayBase::ADD);
      tmp_src_ptr->copy(tmp, 0, 0, ncomp, ngrow, ngrow, Geom(lev).periodicity(),
			FabArrayBase::COPY);
    }
  else
    {
      tmp_src_ptr->SumBoundary(Geom(lev).periodicity());
    }

  // Add new sources into source *after* we have called SumBoundary
  MultiFab::Add(source,*tmp_src_ptr,0,0,source.nComp(),source.nGrow());
  delete tmp_src_ptr;

  // Fill ghost cells after we've synced up ..
  source.FillBoundary(Geom(lev).periodicity());

  // Only delete this if in fact we created it.  Note we didn't change state_ptr so
  //  we don't need to copy anything back into state
  if (state_ptr != &state ) delete state_ptr;

  // ********************************************************************************

  if (this->m_verbose > 1)
    {
      Real stoptime = ParallelDescriptor::second() - strttime;

      ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

      if (ParallelDescriptor::IOProcessor())
        {
	  std::cout << "SprayParticleContainer::moveKick() time: " << stoptime << '\n';
        }
    }
}

Real
SprayParticleContainer::estTimestep (int lev, Real cfl) const
{
  Real dt = 1.e50;

  if (this->GetParticles().size() == 0)
    return dt;

  const Real      strttime         = ParallelDescriptor::second();
  const Geometry& geom             = this->m_gdb->Geom(lev);
  const Real*     dx               = geom.CellSize();
  int   tnum                       = 1;

#ifdef _OPENMP
  tnum = omp_get_max_threads();
#endif

  Vector<Real> ldt(tnum,1e50);

  long num_particles_at_level = 0;
  Real max_mag_vel_over_dx;

  for (MyParConstIter pti(*this, lev); pti.isValid(); ++pti)
    {
      auto& particles = pti.GetArrayOfStructs();
      size_t Np = pti.numParticles();
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (unsigned i = 0; i < Np; ++i)
	{
          const ParticleType& p = particles[i];
          if (p.id() > 0)
	    {
	      max_mag_vel_over_dx = std::max({D_DECL(std::abs(p.rdata(PeleC::pstate_vel))/dx[0],
						     std::abs(p.rdata(PeleC::pstate_vel+1))/dx[1],
						     std::abs(p.rdata(PeleC::pstate_vel+2))/dx[2])});
	    }

	  Real dt_part = (max_mag_vel_over_dx > 0) ? (cfl / max_mag_vel_over_dx) : 1e50;
	  int tid = 0;
#ifdef _OPENMP
	  tid = omp_get_thread_num();
#endif
	  ldt[tid] = std::min(dt_part, ldt[tid]);
	}
      num_particles_at_level += Np;
    }

  for (int i = 0; i < ldt.size(); i++)
    dt = std::min(dt, ldt[i]);

  ParallelDescriptor::ReduceRealMin(dt);

  //
  // Set dt negative if there are no particles at this level.
  //
  ParallelDescriptor::ReduceLongSum(num_particles_at_level);

  if (num_particles_at_level == 0) dt = -1.e50;

  if (this->m_verbose > 1)
    {
      Real stoptime = ParallelDescriptor::second() - strttime;
      ParallelDescriptor::ReduceRealMax(stoptime,
					ParallelDescriptor::IOProcessorNumber());
      if (ParallelDescriptor::IOProcessor())
	std::cout << "SprayParticleContainer::estTimestep() time: "
		  << stoptime << '\n';
    }
  return dt;
}
