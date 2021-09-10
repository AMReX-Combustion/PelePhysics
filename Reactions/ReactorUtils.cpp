#include "ReactorUtils.H"

namespace pele {
namespace physics {
namespace reactions {
namespace utils {

// Check function return value...
//     opt == 0 means SUNDIALS function allocates memory so check if
//              returned NULL pointer
//     opt == 1 means SUNDIALS function returns a flag so check if
//              flag >= 0
//     opt == 2 means function allocates memory so check if returned
//              NULL pointer
int
check_flag(void* flagvalue, const char* funcname, int opt)
{
  int* errflag;

  if (opt == 0 && flagvalue == nullptr) {
    if (amrex::ParallelDescriptor::IOProcessor()) {
      fprintf(
        stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
        funcname);
      amrex::Abort("abort");
    }
    return (1);
  }
  if (opt == 1) {
    errflag = (int*)flagvalue;
    if (*errflag < 0) {
      if (amrex::ParallelDescriptor::IOProcessor()) {
        fprintf(
          stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n", funcname,
          *errflag);
        amrex::Abort("abort");
      }
      return (1);
    }
  } else if (opt == 2 && flagvalue == NULL) {
    if (amrex::ParallelDescriptor::IOProcessor()) {
      fprintf(
        stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
        funcname);
      amrex::Abort("abort");
    }
    return (1);
  }

  return (0);
}
} // namespace utils
} // namespace reactions
} // namespace physics
} // namespace pele
