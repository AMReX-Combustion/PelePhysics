#ifndef UTILITIES_H
#define UTILITIES_H

#include "UnitConversions.H"
namespace pele::physics::utilities {

// -----------------------------------------------------------
// Compute the overlap area between an arbitrary rectangle and circle
// code modified from:
// https://stackoverflow.com/questions/622287/area-of-intersection-between-circle-and-rectangle
// INPUTS:
// x0, x1        => low and high x values of rectangle
// y0, y1        => low and high y values of rectangle
// cx, cy        => x and y values of circle center
// r             => radius of circle
// OUTPUT:
// area of intersection in unit/coordinate systemm of inputs
// -----------------------------------------------------------

// helper function
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
section_area(
  const amrex::Real x0,
  const amrex::Real x1,
  const amrex::Real h,
  const amrex::Real r)
{
  const amrex::Real sec = (h < r) ? std::sqrt(r * r - h * h) : 0;
  const amrex::Real arg1 =
    amrex::max<amrex::Real>(-sec, amrex::min<amrex::Real>(sec, x1));
  const amrex::Real arg2 =
    amrex::max<amrex::Real>(-sec, amrex::min<amrex::Real>(sec, x0));
  const amrex::Real term1 =
    0.5 * (std::sqrt(1.0 - arg1 * arg1 / (r * r)) * arg1 * r +
           r * r * std::asin(arg1 / r) - 2.0 * h * arg1);
  const amrex::Real term2 =
    0.5 * (std::sqrt(1.0 - arg2 * arg2 / (r * r)) * arg2 * r +
           r * r * std::asin(arg2 / r) - 2.0 * h * arg2);
  return term1 - term2;
}

// Area of a square and circle intersection
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
amrex::Real
rectangle_circle_intersection_area(
  const amrex::Real x0,
  const amrex::Real x1,
  const amrex::Real y0,
  const amrex::Real y1,
  const amrex::Real cx,
  const amrex::Real cy,
  const amrex::Real r)
{
  // We assume coordinates of the box are passed in sorted order x0 = xlo, x1 =
  // xhi, etc
  AMREX_ASSERT(x1 > x0);
  AMREX_ASSERT(y1 > y0);
  AMREX_ASSERT(r >= 0.0);

  // Convert to coordinates centered on circle
  const amrex::Real x0new = x0 - cx;
  const amrex::Real x1new = x1 - cx;
  const amrex::Real y0new = y0 - cy;
  const amrex::Real y1new = y1 - cy;

  // Return 0 if fully outside the square that circumscribes the circle
  if (x0new >= r || x1new <= -r || y0new >= r || y1new <= -r) {
    return 0.0;
  }

  // Compute the area
  const amrex::Real sign = copysign(1.0, y0new) * copysign(1.0, y1new);
  const amrex::Real y0mod =
    amrex::min<amrex::Real>(std::abs(y0new), std::abs(y1new));
  const amrex::Real y1mod =
    amrex::max<amrex::Real>(std::abs(y0new), std::abs(y1new));
  const amrex::Real area0 = section_area(x0new, x1new, 0.0, r);
  const amrex::Real area1 = area0 - section_area(x0new, x1new, y0mod, r);
  const amrex::Real area2 = area0 - section_area(x0new, x1new, y1mod, r);

  return area2 - sign * area1;
}

// -----------------------------------------------------------
// Search for the closest index in an array to a given value
// using the bisection technique.
// INPUTS/OUTPUTS:
// xtable(0:n-1) => array to search in (ascending order)
// n             => array size
// x             => x location
// idxlo        <=> output st. xtable(idxlo) <= x < xtable(idxlo+1)
// -----------------------------------------------------------
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
void
locate(const amrex::Real* xtable, const int n, const amrex::Real& x, int& idxlo)
{
  // If x is out of bounds, return boundary index
  if (x >= xtable[n - 1]) {
    idxlo = n - 1;
    return;
  }
  if (x <= xtable[0]) {
    idxlo = 0;
    return;
  }

  // Do the bisection
  idxlo = 0;
  int idxhi = n - 1;
  bool notdone = true;
  while (notdone) {
    if (idxhi - idxlo <= 1) {
      notdone = false;
    } else {
      const int idxmid = (idxhi + idxlo) / 2;
      if (x >= xtable[idxmid]) {
        idxlo = idxmid;
      } else {
        idxhi = idxmid;
      }
    }
  }
}

} // namespace pele::physics::utilities

#endif
