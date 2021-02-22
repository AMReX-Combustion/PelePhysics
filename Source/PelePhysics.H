#ifndef _PELEPHYSICS_H_
#define _PELEPHYSICS_H_

#include <type_traits>
#include <iostream>
#include "PhysicsConstants.H"
#include "EOS.H"
#include "Transport.H"

namespace pele {
namespace physics {

struct LESModel
{
  using eos_type = LESModel;
  using transport_type = LESModel;

  static std::string identifier() { return "Bruce's AwesomeLESModel"; }

  // template<typename Args...>
  // AMREX_GPU_DEVICE
  // AMREX_FORCE_INLINE
  // void RYP2T(amrex::Real R, amrex::Real Y[], amrex::Real P, amrex::Real& T,
  // const Args&&... args){
  //   eos.RYP2T(R, Y, P, T, std::forward<Args>(args)...);
  // }
};

template <typename T>
struct is_eos_type
{
  static constexpr bool value = std::is_base_of<typename T::eos_type, T>::value;
};

template <typename T>
struct is_transport_type
{
  static constexpr bool value =
    std::is_base_of<typename T::transport_type, T>::value;
};

template <typename EosModel, typename TransportModel>
struct PelePhysics
{
  static_assert(
    is_eos_type<EosModel>::value,
    "PelePhysics must have EOS model as its first argument");
  static_assert(
    is_transport_type<TransportModel>::value,
    "PelePhysics must have Transport model as its second argument");

  using eos_type = typename EosModel::eos_type;
  using transport_type = typename TransportModel::transport_type;

  static std::string identifier()
  {
    if (std::is_same<EosModel, TransportModel>::value) {
      return EosModel::identifier();
    } else {
      return EosModel::identifier() + "-" + TransportModel::identifier();
    }
  }

  template <class... Args>
  AMREX_GPU_HOST_DEVICE
  static eos_type eos(Args... args)
  {
    return eos_type(std::forward<Args>(args)...);
  }

  template <class... Args>
  AMREX_GPU_HOST_DEVICE
  static transport_type transport(Args... args)
  {
    return transport_type(std::forward<Args>(args)...);
  }
};

using PhysicsType = PelePhysics<EosType, TransportType>;

} // namespace physics
} // namespace pele

#endif
