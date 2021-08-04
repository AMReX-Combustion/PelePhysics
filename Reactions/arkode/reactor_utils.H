#ifndef _REACTOR_UTILS_H_
#define _REACTOR_UTILS_H_

AMREX_GPU_DEVICE AMREX_FORCE_INLINE int
sunvec_index(const int specie_idx, const int icell, const int ncells)
{
  return specie_idx * ncells + icell;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE int
sunspec_index(const int specie_idx, const int icell, const int ncells)
{
  return sunvec_index(specie_idx, icell, ncells);
}
#endif
