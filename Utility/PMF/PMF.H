#ifndef PMF_H
#define PMF_H

#include <AMReX_REAL.H>
#include <AMReX_GpuContainers.H>
#include <mechanism.H>
#include <PMFData.H>

namespace pele {
namespace physics {
namespace PMF {

AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
pmf(
  PmfData::DataContainer const* pmf_data,
  amrex::Real xlo,
  amrex::Real xhi,
  amrex::GpuArray<amrex::Real, NUM_SPECIES + 4>& y_vector)
{
  // Average the PMF data between xlo and xhi
  // to get a true Finite Volume, cell-centered value
  if (pmf_data->m_doAverage) {
    int lo_loside = -1;
    int lo_hiside = -1;
    if (xlo < pmf_data->pmf_X[0]) {
      lo_loside = 0;
      lo_hiside = 0;
    }
    if (xlo > pmf_data->pmf_X[pmf_data->m_nPoint - 1]) {
      lo_loside = pmf_data->m_nPoint - 1;
      lo_hiside = pmf_data->m_nPoint - 1;
    }
    int hi_loside = -1;
    int hi_hiside = -1;
    if (xhi < pmf_data->pmf_X[0]) {
      hi_loside = 0;
      hi_hiside = 0;
    }
    if (xhi > pmf_data->pmf_X[pmf_data->m_nPoint - 1]) {
      hi_loside = pmf_data->m_nPoint - 1;
      hi_hiside = pmf_data->m_nPoint - 1;
    }
    if (lo_loside == -1) {
      for (unsigned int i = 0; i < pmf_data->m_nPoint - 1; i++) {
        if ((xlo >= pmf_data->pmf_X[i]) && (xlo <= pmf_data->pmf_X[i + 1])) {
          lo_loside = i;
          lo_hiside = i + 1;
        }
      }
    }
    if (hi_loside == -1) {
      for (unsigned int i = 0; i < pmf_data->m_nPoint - 1; i++) {
        if ((xhi >= pmf_data->pmf_X[i]) && (xhi <= pmf_data->pmf_X[i + 1])) {
          hi_loside = i;
          hi_hiside = i + 1;
        }
      }
    }
    for (unsigned int j = 0; j < pmf_data->m_nVar; j++) {
      amrex::Real x1 = pmf_data->pmf_X[lo_loside];
      amrex::Real y1 = pmf_data->pmf_Y[pmf_data->m_nPoint * j + lo_loside];
      amrex::Real x2 = pmf_data->pmf_X[lo_hiside];
      amrex::Real y2 = pmf_data->pmf_Y[pmf_data->m_nPoint * j + lo_hiside];
      amrex::Real dydx = 0.0;
      if (lo_loside != lo_hiside) {
        dydx = (y2 - y1) / (x2 - x1);
      }
      amrex::Real ylo = y1 + dydx * (xlo - x1);
      if (lo_loside == hi_loside) {
        yhi = y1 + dydx * (xhi - x1);
        y_vector[j] = 0.5 * (ylo + yhi);
      } else {
        amrex::Real sum = (x2 - xlo) * 0.5 * (ylo + y2);
        x1 = pmf_data->pmf_X[hi_loside];
        y1 = pmf_data->pmf_Y[pmf_data->m_nPoint * j + hi_loside];
        x2 = pmf_data->pmf_X[hi_hiside];
        y2 = pmf_data->pmf_Y[pmf_data->m_nPoint * j + hi_hiside];
        if (lo_loside == lo_hiside) {
          dydx = 0.0;
        } else {
          dydx = (y2 - y1) / (x2 - x1);
        }
        yhi = y1 + dydx * (xhi - x1);
        sum += (xhi - x1) * 0.5 * (yhi + y1);
        for (int k = lo_hiside; k <= hi_loside - 1; k++) {
          sum += (pmf_data->pmf_X[k + 1] - pmf_data->pmf_X[k]) * 0.5 *
                 (pmf_data->pmf_Y[pmf_data->m_nPoint * j + k] +
                  pmf_data->pmf_Y[pmf_data->m_nPoint * j + k + 1]);
        }
        y_vector[j] = sum / (xhi - xlo);
      }
    }

    // Get the data from interpolating PMF data
    // to the center of [xlo:xli]
  } else {
    amrex::Real xmid = 0.5 * (xlo + xhi);
    int loside = -1;
    int hiside = -1;
    if (xmid < pmf_data->pmf_X[0]) {
      loside = 0;
      hiside = 0;
    }
    if (xmid > pmf_data->pmf_X[pmf_data->m_nPoint - 1]) {
      loside = pmf_data->m_nPoint - 1;
      hiside = pmf_data->m_nPoint - 1;
    }
    if (loside == -1) {
      for (unsigned int i = 0; i < pmf_data->m_nPoint - 1; i++) {
        if ((xmid >= pmf_data->pmf_X[i]) && (xmid <= pmf_data->pmf_X[i + 1])) {
          loside = i;
          hiside = i + 1;
        }
      }
    }
    for (unsigned int j = 0; j < pmf_data->m_nVar; j++) {
      amrex::Real x1 = pmf_data->pmf_X[loside];
      amrex::Real y1 = pmf_data->pmf_Y[pmf_data->m_nPoint * j + loside];
      amrex::Real x2 = pmf_data->pmf_X[hiside];
      amrex::Real y2 = pmf_data->pmf_Y[pmf_data->m_nPoint * j + hiside];
      amrex::Real dydx = 0.0;
      if (loside != hiside) {
        dydx = (y2 - y1) / (x2 - x1);
      }
      y_vector[j] = y1 + dydx * (xlo - x1);
    }
  }
}
} // namespace PMF
} // namespace physics
} // namespace pele
#endif
