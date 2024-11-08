#ifndef EOS_H
#define EOS_H

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <AMReX_Utility.H>
#include <AMReX_Arena.H>
#include <AMReX_ParmParse.H>

#include "mechanism.H"
#include "PhysicsConstants.H"
#include "EosParams.H"
#include "GammaLaw.H"
#include "Fuego.H"
#include "SRK.H"
#ifndef AMREX_USE_SYCL
#include "Manifold.H"
#endif

namespace pele::physics {
#ifdef USE_GAMMALAW_EOS
using EosType = eos::GammaLaw;
#elif USE_FUEGO_EOS
using EosType = eos::Fuego;
#elif USE_SRK_EOS
using EosType = eos::SRK;
#elif USE_MANIFOLD_EOS && !defined(AMREX_USE_SYCL)
using EosType = eos::Manifold;
#else
static_assert(false, "Invalid EOS specified");
#endif
namespace eos {

template <typename EOSType>
void
atomic_weightsCHON(amrex::Real atwCHON[])
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

template <typename EOSType>
void
element_compositionCHON(int ecompCHON[])
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
inline void
atomic_weightsCHON<GammaLaw>(amrex::Real* /*atwCHON[]*/)
{
  amrex::Error("atomic_weightsCHON not supported for GammaLaw EOS");
}

template <>
inline void
element_compositionCHON<GammaLaw>(int* /*ecompCHON[]*/)
{
  amrex::Error("element_compositionCHON not supported for GammaLaw EOS");
}

template <>
inline void
atomic_weightsCHON<Manifold>(amrex::Real* /*atwCHON[]*/)
{
  amrex::Error("atomic_weightsCHON not supported for Manifold EOS");
}

template <>
inline void
element_compositionCHON<Manifold>(int* /*ecompCHON[]*/)
{
  amrex::Error("element_compositionCHON not supported for Manifold EOS");
}

template <typename EOSType>
void
speciesNames(
  amrex::Vector<std::string>& spn, const EosParm<EOSType>* /*eparm*/ = nullptr)
{
  CKSYMS_STR(spn);
}

#ifndef AMREX_USE_SYCL
template <>
inline void
speciesNames<Manifold>(
  amrex::Vector<std::string>& spn, const EosParm<Manifold>* eparm)
{
  if (eparm == nullptr) {
    amrex::Abort("Manifold EOS requires an EosParm input in order to determine "
                 "species names");
  }
  spn.resize(NUM_SPECIES);
  const auto* mani_data = eparm->manf_data;
  for (int n = 0; n < MANIFOLD_DIM; n++) {
    std::string nametmp = std::string(
      &(mani_data->dimnames)[n * mani_data->len_str], mani_data->len_str);
    spn[n] = amrex::trim(nametmp);
  }
  spn[MANIFOLD_DIM] = "XRHO";
}
#endif

} // namespace eos
} // namespace pele::physics
#endif
