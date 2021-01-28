#include <AMReX.H>
#include <sundials/sundials_memory.h>

namespace amrex {
namespace sundials {

class MemoryHelper {
public:
  MemoryHelper();
  ~MemoryHelper();
  operator SUNMemoryHelper() const { return helper; }
  static void Initialize();
  static void Finalize();
private:
  SUNMemoryHelper helper;
};

MemoryHelper* The_SUNMemory_Helper ();

}
}
