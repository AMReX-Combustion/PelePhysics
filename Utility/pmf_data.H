#ifndef PMF_DATA_H_
#define PMF_DATA_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuContainers.H>
#include <mechanism.h>

// Define a data structure to hold the PMF data,
struct PmfData {
  unsigned int pmf_N;
  unsigned int pmf_M;
  int pmf_do_average;
  amrex::Real* pmf_X;
  amrex::Real* pmf_Y;
};

extern PmfData* pmf_data_g;
#endif
