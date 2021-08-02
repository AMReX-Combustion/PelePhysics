#ifndef _REACTOR_UTILS_H_
#define _REACTOR_UTILS_H_

AMREX_GPU_DEVICE AMREX_FORCE_INLINE int
sunvec_index(const int specie_idx, const int icell, const int /*ncells*/)
{
  return icell * (NUM_SPECIES + 1) + specie_idx;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE int
sunspec_index(const int specie_idx, const int icell, const int /*ncells*/)
{
  return icell * NUM_SPECIES + specie_idx;
}
#endif
