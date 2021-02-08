#include "EOS.H"

namespace EOS {

void
init()
{
  CKINIT();
}

void 
close()
{
  CKFINALIZE();
}

void
speciesNames(amrex::Vector<std::string>& spn) {
  spn.resize(1);
  spn[0] = "AIR";
}

void
atomic_weightsCHON(amrex::Real atwCHON[])
{
  //CHON
  for (int i = 0; i < 4; i++) {
      atwCHON[i] = 0.0;
  }
}

void
element_compositionCHON(int ecompCHON[])
{
  for (int k = 0; k < 4; k++) {
      ecompCHON[k] = 0;
  }
}

} // namespace EOS
