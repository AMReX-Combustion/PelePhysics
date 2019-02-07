#include <AMReX_BLFort.H>

#ifdef __cplusplus
extern "C"
{
#endif
    void update_particles(const int* np, const int* lev, void* particles,
                          const amrex::Real* state, const int* state_lo, const int* state_hi,
                          const amrex::Real* source, const int* source_lo, const int* source_hi,
                          const int* domain_lo, const int* domain_hi,
                          const amrex::Real* prob_lo, const amrex::Real* prob_hi, 
                          const int* reflect_lo, const int* reflect_hi,
                          const amrex::Real* dx, const amrex::Real& dt, const int* do_move);
#ifdef __cplusplus
}
#endif
