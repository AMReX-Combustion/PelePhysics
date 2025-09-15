#include <AMReX.H>
#include <AMReX_Print.H>
#include "Utilities.H"
#include "Table.H"
#include "BlackBoxFunction.H"

namespace m2c = pele::physics::utilities::mks2cgs;
namespace c2m = pele::physics::utilities::cgs2mks;

int
main(int argc, char* argv[])
{
  amrex::Initialize(argc, argv);

  // Test table interpolation
  pele::physics::TabulatedFunctionParams tabfunc_par;
  tabfunc_par.host_only_parm().parm_parse_prefix = "test";
  tabfunc_par.initialize();
  pele::physics::TabulatedFunction tabfunc(
    static_cast<const pele::physics::TabulatedFunctionData*>(
      &tabfunc_par.host_parm()));
  amrex::Real outval;
  const amrex::GpuArray<amrex::Real, 2> invals{1.0, 1.0};
  const std::string varname = "T";
  const int idxT =
    pele::physics::get_var_index(varname.c_str(), &tabfunc_par.host_parm());
  tabfunc.get_value(idxT, invals.data(), outval);
  AMREX_ALWAYS_ASSERT(amrex::almostEqual(outval, 3.60));
  tabfunc_par.deallocate();
  amrex::Print() << "TabulatedFunction Tests passed" << std::endl;

  // Test Unit Conversions
  AMREX_ALWAYS_ASSERT(m2c::Length(1.0) == 100.0);
  AMREX_ALWAYS_ASSERT(c2m::Length(1.0, 2) == 1.0e-4);
  AMREX_ALWAYS_ASSERT(m2c::Mass(c2m::Mass(1.0)) == 1.0);
  AMREX_ALWAYS_ASSERT(m2c::Mu(1.0) == m2c::Mass(1.0) / m2c::Length(1.0));
  amrex::Print() << "Unit Conversion Tests Passed" << std::endl;

  // Test locate function
  amrex::Array<amrex::Real, 5> sorted_data{0.0, 0.2, 0.5, 0.9, 1.0};
  int length = sorted_data.size();
  int idxlo = -1;
  amrex::Real x = -50.0;
  pele::physics::utilities::locate(sorted_data.data(), length, x, idxlo);
  AMREX_ALWAYS_ASSERT(idxlo == 0);
  x = 1.0;
  pele::physics::utilities::locate(sorted_data.data(), length, x, idxlo);
  AMREX_ALWAYS_ASSERT(idxlo == length - 1);
  x = 0.25;
  pele::physics::utilities::locate(sorted_data.data(), length, x, idxlo);
  AMREX_ALWAYS_ASSERT(idxlo == 1);
  amrex::Print() << "Locate Tests Passed" << std::endl;

  // Test rectangle circle intersection area
  // Quarter circle intersection
  amrex::Real x0 = 0.0;
  amrex::Real x1 = 1.0;
  amrex::Real y0 = 0.0;
  amrex::Real y1 = 1.0;
  amrex::Real cx = 1.0;
  amrex::Real cy = 1.0;
  amrex::Real r = 1.0;
  amrex::Real A1 = pele::physics::utilities::rectangle_circle_intersection_area(
    x0, x1, y0, y1, cx, cy, r);
  AMREX_ALWAYS_ASSERT(std::abs(A1 - 0.25 * M_PI) < 1e-14);
  // No intersection
  cx = 3.0;
  cy = 3.0;
  amrex::Real A2 = pele::physics::utilities::rectangle_circle_intersection_area(
    x0, x1, y0, y1, cx, cy, r);
  AMREX_ALWAYS_ASSERT(std::abs(A2) < 1e-14);
  // Square covered by circle
  r = 20.0;
  amrex::Real A3 = pele::physics::utilities::rectangle_circle_intersection_area(
    x0, x1, y0, y1, cx, cy, r);
  AMREX_ALWAYS_ASSERT(std::abs(A3 - 1.0) < 1e-14);
  amrex::Print() << "Intersection Area Tests Passed" << std::endl;

  amrex::Finalize();
  return 0;
}
