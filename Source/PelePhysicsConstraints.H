#include <type_traits>

namespace pele::physics {

// Forward declarations
namespace eos {
struct Fuego;
struct GammaLaw;
struct SRK;
} // namespace eos

// Forward declarations
namespace transport {
struct ConstTransport;
struct SimpleTransport;
struct SutherlandTransport;
} // namespace transport

template <typename EosModel, typename TransportModel>
struct is_valid_physics_combination : public std::true_type
{
};

template <>
struct is_valid_physics_combination<eos::GammaLaw, transport::SimpleTransport>
  : public std::false_type
{
};

template <>
struct is_valid_physics_combination<
  eos::GammaLaw,
  transport::SutherlandTransport> : public std::false_type
{
};

template <>
struct is_valid_physics_combination<eos::SRK, transport::ConstTransport>
  : public std::false_type
{
};

template <>
struct is_valid_physics_combination<eos::SRK, transport::SutherlandTransport>
  : public std::false_type
{
};

} // namespace pele::physics
