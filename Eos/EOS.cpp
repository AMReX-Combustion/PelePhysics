#include "EOS.H"

namespace pele {
namespace physics {
namespace eos {

template <typename EOSType>
void
atomic_weightsCHON(amrex::Real* /*atwCHON*/)
{
}

template <typename EOSType>
void
element_compositionCHON(int* /*ecompCHON*/)
{
}

template <typename EOSType>
void
speciesNames(amrex::Vector<std::string>& /*spn*/)
{
}

template <>
void
atomic_weightsCHON<GammaLaw>(amrex::Real atwCHON[])
{
  // CHON
  for (int i = 0; i < 4; i++) {
    atwCHON[i] = 0.0;
  }
}

template <>
void
element_compositionCHON<GammaLaw>(int ecompCHON[])
{
  for (int k = 0; k < 4; k++) {
    ecompCHON[k] = 0;
  }
}

template <>
void
speciesNames<GammaLaw>(amrex::Vector<std::string>& spn)
{
  spn.resize(1);
  spn[0] = "AIR";
}

template <>
void
atomic_weightsCHON<Fuego>(amrex::Real atwCHON[])
{
  amrex::Vector<std::string> ename;
  CKSYME_STR(ename);
  amrex::Real atw[NUM_ELEMENTS];
  CKAWT(atw);
  // CHON
  for (int i = 0; i < 4; i++) {
    atwCHON[i] = 0.0;
  }
  for (int i = 0; i < NUM_ELEMENTS; i++) {
    if (ename[i] == "C") {
      atwCHON[0] = atw[i];
    }
    if (ename[i] == "H") {
      atwCHON[1] = atw[i];
    }
    if (ename[i] == "O") {
      atwCHON[2] = atw[i];
    }
    if (ename[i] == "N") {
      atwCHON[3] = atw[i];
    }
  }
}

template <>
void
element_compositionCHON<Fuego>(int ecompCHON[])
{
  amrex::Vector<std::string> ename;
  CKSYME_STR(ename);
  // CHON
  int CHON[4] = {-1};
  for (int i = 0; i < NUM_ELEMENTS; i++) {
    if (ename[i] == "C") {
      CHON[0] = i;
    }
    if (ename[i] == "H") {
      CHON[1] = i;
    }
    if (ename[i] == "O") {
      CHON[2] = i;
    }
    if (ename[i] == "N") {
      CHON[3] = i;
    }
  }
  int ecomp[NUM_SPECIES * NUM_ELEMENTS];
  CKNCF(ecomp);
  for (int i = 0; i < NUM_SPECIES; i++) {
    for (int k = 0; k < 4; k++) {
      if (CHON[k] > -1) {
        ecompCHON[i * 4 + k] = ecomp[i * NUM_ELEMENTS + CHON[k]];
      } else {
        ecompCHON[i * 4 + k] = 0;
      }
    }
  }
}

template <>
void
speciesNames<Fuego>(amrex::Vector<std::string>& spn)
{
  CKSYMS_STR(spn);
}

template <>
void
atomic_weightsCHON<SRK>(amrex::Real atwCHON[])
{
  atomic_weightsCHON<Fuego>(atwCHON);
}

template <>
void
element_compositionCHON<SRK>(int ecompCHON[])
{
  element_compositionCHON<Fuego>(ecompCHON);
}

template <>
void
speciesNames<SRK>(amrex::Vector<std::string>& spn)
{
  speciesNames<Fuego>(spn);
}

} // namespace eos
} // namespace physics
} // namespace pele
