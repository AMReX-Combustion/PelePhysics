#ifndef PMF_H_
#define PMF_H_

#include <AMReX_REAL.H>
#include <AMReX_GpuContainers.H>
#include <mechanism.h>

namespace PMF {
  
  extern AMREX_GPU_DEVICE_MANAGED unsigned int pmf_N;
  extern AMREX_GPU_DEVICE_MANAGED unsigned int pmf_M;
  extern AMREX_GPU_DEVICE_MANAGED bool pmf_do_average;

  extern amrex::Gpu::ManagedVector<amrex::Real> pmf_X;
  extern amrex::Gpu::ManagedVector<amrex::Real> pmf_Y;

  extern amrex::Vector<std::string> pmf_names;

  void read_pmf(const std::string& myfile);
  
  AMREX_GPU_HOST_DEVICE
  AMREX_FORCE_INLINE
  void
  pmf(
    amrex::Real xlo,
    amrex::Real xhi,
    amrex::GpuArray<amrex::Real, NUM_SPECIES + 4>& y_vector)
  {
    amrex::Real sum = 0.0, xmid = 0.0;
    int lo_loside = 0, lo_hiside = 0;
    int hi_loside = 0, hi_hiside = 0;
    int loside = 0, hiside = 0;
    amrex::Real ylo = 0.0, yhi = 0.0, x1 = 0.0, y1 = 0.0, x2 = 0.0, y2 = 0.0,
      dydx = 0.0;

    auto d_pmf_X = PMF::pmf_X.dataPtr();
    auto d_pmf_Y = PMF::pmf_Y.dataPtr();
    
    if (pmf_do_average) {
      if (xlo < d_pmf_X[0]) {
        lo_loside = 0;
        lo_hiside = 0;
      }
      if (xhi < d_pmf_X[0]) {
        hi_loside = 0;
        hi_hiside = 0;
      }
      if (xlo > d_pmf_X[pmf_N - 1]) {
        lo_loside = pmf_N - 1;
        lo_hiside = pmf_N - 1;
      }
      if (xhi > d_pmf_X[pmf_N - 1]) {
        hi_loside = pmf_N - 1;
        hi_hiside = pmf_N - 1;
      }
      if (lo_loside == 0) {
        for (int i = 0; i < pmf_N - 1; i++) {
          if ((xlo > d_pmf_X[i]) && (xlo < d_pmf_X[i + 1])) {
            lo_loside = i;
            lo_hiside = i + 1;
          }
        }
      }
      if (hi_loside == 0) {
        for (int i = 0; i < pmf_N - 1; i++) {
          if ((xhi > d_pmf_X[i]) && (xhi < d_pmf_X[i + 1])) {
            hi_loside = i;
            hi_hiside = i + 1;
          }
        }
      }
      for (int j = 0; j < pmf_M; j++) {
        x1 = d_pmf_X[lo_loside];
        y1 = d_pmf_Y[pmf_N * j + lo_loside];
        x2 = d_pmf_X[lo_hiside];
        y2 = d_pmf_Y[pmf_N * j + lo_hiside];
        dydx = lo_loside == lo_hiside ? 0 : (y2 - y1) / (x2 - x1);
        ylo = y1 + dydx * (xlo - x1);
        if (lo_loside == hi_loside) {
          yhi = y1 + dydx * (xhi - x1);
          y_vector[j] = 0.5 * (ylo + yhi);
        } else {
          sum = (x2 - xlo) * 0.5 * (ylo + y2);
          x1 = d_pmf_X[hi_loside];
          y1 = d_pmf_Y[pmf_N * j + hi_loside];
          x2 = d_pmf_X[hi_hiside];
          y2 = d_pmf_Y[pmf_N * j + hi_hiside];
          dydx = hi_loside == hi_hiside ? 0.0 : (y2 - y1) / (x2 - x1);
          yhi = y1 + dydx * (xhi - x1);
          sum = sum + (xhi - x1) * 0.5 * (yhi + y1);
          for (int k = lo_hiside; k < hi_loside - 1; k++) {
            sum = sum + (d_pmf_X[k + 1] - d_pmf_X[k]) * 0.5 *
              (d_pmf_Y[pmf_N * j + k] +
               d_pmf_Y[pmf_N * j + k + 1]);
          }
          y_vector[j] = sum / (xhi - xlo);
        }
      }
    } else {
      xmid = 0.5 * (xlo + xhi);
      loside = -1;
      hiside = -1;
      if (xmid < d_pmf_X[0]) {
        loside = 0;
        hiside = 0;
      }
      if (xmid > d_pmf_X[pmf_N - 1]) {
        loside = pmf_N - 1;
        hiside = pmf_N - 1;
      }
      if (loside == -1) {
        for (int i = 0; i < pmf_N - 1; i++) {
          if (
            (xmid >= d_pmf_X[i]) &&
            (xmid <= d_pmf_X[i + 1])) {
            loside = i;
            hiside = i + 1;
          }
        }
      }
      for (int j = 0; j < pmf_M; j++) {
        x1 = d_pmf_X[loside];
        y1 = d_pmf_Y[pmf_N * j + loside];
        x2 = d_pmf_X[hiside];
        y2 = d_pmf_Y[pmf_N * j + hiside];
        dydx = loside == hiside  ?  0.0  :  (y2 - y1) / (x2 - x1);
        y_vector[j] = y1 + dydx * (xlo - x1);
      }
    }
  }
}
#endif
