#ifndef _PHYSICSCONSTANTS_H_
#define _PHYSICSCONSTANTS_H_

#include <AMReX_REAL.H>

namespace pele {
namespace physics {

struct Constants
{
  static constexpr amrex::Real gamma = 1.4;
  static constexpr amrex::Real RU = 8.31446261815324e7;
  static constexpr amrex::Real RUC = 1.98721558317399615845;
  static constexpr amrex::Real PATM = 1.01325e+06;
  static constexpr amrex::Real AIRMW = 28.97;
  static constexpr amrex::Real Avna = 6.022140857e23;
};
} // namespace physics
} // namespace pele
#endif
