#include <AMReX_PhysBCFunct.H>
#include <AMReX_Geometry.H>

namespace pele {
namespace physics {
namespace pltfilemanager {

struct FillExtDirDummy
{

  AMREX_GPU_HOST
  constexpr FillExtDirDummy() = default;

  AMREX_GPU_DEVICE
  void operator()(
    const amrex::IntVect& /*iv*/,
    amrex::Array4<amrex::Real> const& /*dummy*/,
    const int /*dcomp*/,
    const int /*numcomp*/,
    amrex::GeometryData const& /*geom*/,
    const amrex::Real /*time*/,
    const amrex::BCRec* /*bcr*/,
    const int /*bcomp*/,
    const int /*orig_comp*/) const
  {
  }
};

} // namespace pltfilemanager
} // namespace physics
} // namespace pele
