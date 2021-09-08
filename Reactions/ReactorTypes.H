#ifndef _REACTORTYPES_H_
#define _REACTORTYPES_H_

#include <AMReX.H>

namespace pele {
namespace physics {
namespace reactions {

struct ReactorTypes
{
  static constexpr int e_reactor_type = 1;
  static constexpr int h_reactor_type = 2;
  static void check_reactor_type(const int reactor_type)
  {
    if (!((reactor_type == ReactorTypes::e_reactor_type) ||
          (reactor_type == ReactorTypes::h_reactor_type))) {
      amrex::Abort("Wrong reactor type. Choose between 1 (e) or 2 (h).");
    }
  };
};
} // namespace reactions
} // namespace physics
} // namespace pele
#endif
