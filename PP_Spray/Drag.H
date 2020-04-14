#ifndef DRAG_H_
#define DRAG_H_

#include <AMReX_REAL.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_AmrParticles.H>
#include <AMReX_Particles.H>
#include <PeleC.H>
#include <SprayParticles.H>

using namespace amrex;

AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
void AdjIndexWeights(const RealVect& len,
                     IntVect*        indx_array,
                     Real*           coef,
                     const IntVect&  dom_lo,
                     const IntVect&  dom_hi)
{
  // Compute distances for trilinear interpolation
  IntVect indx(AMREX_D_DECL(amrex::Math::floor(len[0]),
                            amrex::Math::floor(len[1]),
                            amrex::Math::floor(len[2])));
  RealVect delL(AMREX_D_DECL(len[0] - indx[0], len[1] - indx[1], len[2] - indx[2]));
  Real ssx[2] = {1._rt - delL[0], delL[0]};
  int aindx = 0;
  IntVect cur_indx;
#if AMREX_SPACEDIM > 1
  Real ssy[2] = {1._rt - delL[1], delL[1]};
#if AMREX_SPACEDIM == 3
  Real ssz[2] = {1._rt - delL[2], delL[2]};
  for (int kk = 0; kk != 2; ++kk) {
    cur_indx[2] =
      amrex::min(dom_hi[2], amrex::max(dom_lo[2], indx[2] + kk - 1));
#endif
    for (int jj = 0; jj != 2; ++jj) {
      cur_indx[1] =
        amrex::min(dom_hi[1], amrex::max(dom_lo[1], indx[1] + jj - 1));
#endif
      for (int ii = 0; ii != 2; ++ii) {
        cur_indx[0] =
          amrex::min(dom_hi[0], amrex::max(dom_lo[0], indx[0] + ii - 1));
        AMREX_D_TERM(indx_array[aindx][0] = cur_indx[0];,
                     indx_array[aindx][1] = cur_indx[1];,
                     indx_array[aindx][2] = cur_indx[2];);
	coef[aindx] = AMREX_D_TERM(ssx[ii],*ssy[jj],*ssz[kk]);
	aindx++;
      }
#if AMREX_SPACEDIM > 1
    }
#if AMREX_SPACEDIM == 3
  }
#endif
#endif
}

// Compute the inverse of the correction factors
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
Real calcInvCorrectionFactor(const Real& B,
			     const Real& logB)
{
  return B/(logB*std::pow(1. + B, 0.7));
}

// Calculate the fuel vapor mass fraction
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
Real calcVaporMF(const Real& h_L,
		 const Real& T_part,
		 const Real& p_fluid,
		 const Real& mw_mix,
		 const Real& mw_fuel,
		 const Real& inv_boil_temp,
		 const Real& inv_Ru,
		 const Real& p_ref,
		 const Real& C_eps)
{
  Real pres_sat = p_ref*std::exp(h_L*inv_Ru*mw_fuel*
				 (inv_boil_temp - 1./T_part)) + C_eps;
  Real Yfv = mw_fuel*pres_sat/(mw_mix*p_fluid +
			       (mw_fuel - mw_mix)* pres_sat);
  return amrex::min(1. - C_eps, Yfv);
}

// Calculate evaporation rate and Sherwood number
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
Real calcSpecEvapRate(const Real& dia_part,
		      const Real& B_M_num,
		      const Real& logB,
		      const Real& Reyn,
		      const Real& powR,
		      const Real& Sc_skin,
		      const Real& D_skin)
{
  Real invFM = calcInvCorrectionFactor(B_M_num, logB);
  Real Sh_0 = 1. + powR*std::cbrt(1. + Reyn*Sc_skin);
  return 2. + (Sh_0 - 2.)*invFM;
}

// Compute the heat transfer coefficient using the
// corrected Nusselt number and B_T value
AMREX_GPU_DEVICE
AMREX_FORCE_INLINE
Real calcHeatCoeff(const Real& ratio,
		   const Real& B_M,
		   const Real& B_eps,
		   const Real& Nu_0)
{
  Real Nu_num = Nu_0;
  Real phi = ratio/Nu_0;
  Real BM1 = 1. + B_M;
  Real B_T = std::pow(BM1, phi) - 1.;
  const int maxIter = 100;
  int k = 0;
  Real logB;
  Real B_T_old;
  Real invFT;
  Real error = B_eps*1000.;
  while (k < maxIter && error > B_eps) {
    B_T_old = B_T;
    logB = std::log(1. + B_T);
    invFT = calcInvCorrectionFactor(B_T, logB);
    Nu_num = 2. + (Nu_0 - 2.)*invFT;
    phi = ratio/Nu_num;
    B_T = std::pow(BM1, phi) - 1.;
    error = std::abs(B_T - B_T_old)/B_T;
    k++;
  }
  logB = std::log(1. + B_T);
  invFT = calcInvCorrectionFactor(B_T, logB);
  Nu_num = 2. + (Nu_0 - 2.)*invFT;
  return Nu_num*logB/B_T;
}

// Previous wall functionality
//   if (do_move)
//     {
//       for (int dir = 0; dir != AMREX_SPACEDIM; ++dir)
// 	{
// 	  if (reflect_lo[dir])
// 	    {
// 	      // Loop over all particles
// 	      for (int n = 0; n != Np; ++n)
// 		{
// 		  if (particles[n].id() != -1 && 
// 		      particles[n].pos(dir) < plo[dir])
// 		    {
// 		      particles[n].pos(dir) *= -1.;
// 		      particles[n].pos(dir) += 2.*plo[dir];
// 		      particles[n].rdata(PeleC::pstate_vel + dir) *= -1.;
// 		    }
// 		}
// 	    }
// 	  if (reflect_hi[dir])
// 	    {
// 	      for (int n = 0; n != Np; ++n)
// 		{
// 		  if (particles[n].id() != -1 && 
// 		      particles[n].pos(dir) > phi[dir])
// 		    {
// 		      particles[n].pos(dir) *= -1.;
// 		      particles[n].pos(dir) += 2.*phi[dir];
// 		      particles[n].rdata(PeleC::pstate_vel + dir) *= -1.;
// 		    }
// 		}
// 	    }
// 	}
//     }
//   for (int dir = 0; dir != AMREX_SPACEDIM; ++dir)
//     {
//       if (reflect_lo[dir])
// 	{
// 	  Box extDomBox = adjCellLo(domain, dir);
// 	  // Extract box that is outside the domain but part of src_box
// 	  Box extCells = src_box & extDomBox;
// 	  if (!extCells.isEmpty())
// 	    {
// 	      // Extract cells adjacent to the domain edges
// 	      Box intCells = extCells;
// 	      intCells.shift(dir, 1);
// 	      source.atomicAdd(source, extCells, intCells,
// 			       0, 0, PeleC::NUM_STATE);
// 	    }
// 	}
//       if (reflect_hi[dir])
// 	{
// 	  Box extDomBox = adjCellHi(domain, dir);
// 	  // Extract box that is outside the domain but part of src_box
// 	  Box extCells = src_box & extDomBox;
// 	  if (!extCells.isEmpty())
// 	    {
// 	      // Extract cells adjacent to the domain edges
// 	      Box intCells = extCells;
// 	      intCells.shift(dir, -1);
// 	      source.atomicAdd(source, extCells, intCells,
// 			       0, 0, PeleC::NUM_STATE);
// 	    }
// 	}
//     }
// }

#endif
