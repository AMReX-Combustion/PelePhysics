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
#include "Manifold.H"

namespace pele::physics {
#ifdef USE_GAMMALAW_EOS
using EosType = eos::GammaLaw;
#elif USE_FUEGO_EOS
using EosType = eos::Fuego;
#elif USE_SRK_EOS
using EosType = eos::SRK;
#elif USE_MANIFOLD_EOS
using EosType = eos::Manifold;
#else
static_assert(false, "Invalid EOS specified");
#endif
namespace eos {

template <typename EOSType>
void atomic_weightsCHON(amrex::Real atwCHON[]);

template <typename EOSType>
void element_compositionCHON(int ecompCHON[]);

template <typename EOSType>
void speciesNames(amrex::Vector<std::string>& spn);

} // namespace eos
} // namespace pele::physics
#endif
