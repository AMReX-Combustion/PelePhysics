#include "ReactorCvodeUtils.H"

namespace pele {
namespace physics {
namespace reactions {
namespace cvode {

// Error function for CVODE
void
cvodeErrHandler(
  int error_code,
  const char* /*module*/,
  const char* /*function*/,
  char* msg,
  void* /*eh_data*/)
{
  if (error_code != CV_WARNING) {
    std::cout << "From CVODE: " << msg << std::endl;
    amrex::Abort("Aborting from CVODE");
  }
}

} // namespace cvode
} // namespace reactions
} // namespace physics
} // namespace pele
